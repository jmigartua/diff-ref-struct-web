//! The least-squares engine (P42–P44, L053), run as a session that can be stepped, paused and resumed (W11, L109) and
//! steered (W12, `steering`).
//!
//! Minimises S = Σ_i w_i (y_i − y_c,i)² over the points outside the excluded regions, with the weights of Iwg
//! (`Weighting`; the maximum-likelihood weights 1/σ²_calc are recomputed from Ycal each cycle), by Gauss–Newton steps
//! (FullProf manual eqs. 3.5–3.6, p. 43) with the analytic derivatives of `jacobian`, damped Levenberg–Marquardt style
//! only when a step would raise S:
//! - the normal matrix A = JᵀWJ is scaled to unit diagonal (Â = D⁻¹AD⁻¹, D² = diag A) and decomposed into its
//!   eigenvectors (Jacobi rotations): Â = VΛVᵀ, the singular-value decomposition of the scaled weighted Jacobian
//!   squared. Directions with λ < 10⁻¹² λ_max (complete correlations) are a singular matrix: the session stops and
//!   names the combination as `check` does (R24), or, with the guard off, leaves them out of the step and reports them;
//! - the step is δ = D⁻¹V(Λ + λ)⁻¹VᵀD⁻¹g, λ = the steering's start (0: Gauss–Newton) first; if S rises, λ = 10⁻³ and
//!   ×10 until it falls (up to 10⁸), each try written in the cycle's notes; then ÷10 at each accepted cycle. With
//!   `Steering::marquardt` off every step is applied, as FullProf applies them;
//! - the step applied is δ times each parameter's damping, shortened to the values' maximum shifts (W12);
//! - limits (P39): a step that would take a bounded field past its bound is cut back to the bound for the parameter
//!   that moves only that field, and the parameter is held there (left out of the solve) while the gradient pushes it
//!   out; such parameters are reported as at their limit. A field whose bound reflects or is periodic is reflected
//!   back inside or wrapped round instead (`steering::Behaviour`);
//! - convergence: every |computed shift|/esd below ε (the `.pcr`'s Eps, FullProf's criterion, manual p. 77), on the
//!   shift computed, before damping;
//! - divergence (`steering::Guards`): χ² rising k cycles in a row, a computed shift of many esd's, an impossible value
//!   (`impossible`);
//! - esd's: σ_j = √((A⁻¹)_jj χ²_ν), χ²_ν = S/(N − P + C) over all points used (eq. 3.7–3.8); FullProf takes χ²_ν over
//!   the points with Bragg contributions only (p. 43), so its esd's can be larger; the esd of a field is
//!   √(m_fᵀ C m_f) with its multipliers m_f. Damping never changes an esd; FullProf's printed esd of a parameter with a
//!   multiplier a is a·σ (its eq. 3.7): `fullprof_esd`.
//!
//! Every cycle is an event (`Cycle`): for each free parameter its value, shift, esd, shift/esd and computed shift; the
//! cycle's λ, Ycal, the R-factors, the restraints with their weights and contributions, and notes. Every stop has its
//! reason as data (`Stop`).

use crate::data::Pattern;
use crate::jacobian::Field;
use crate::powder::{Calculated, Model, Weighting};
use crate::refinable::ParameterSet;
use crate::steering::Behaviour;

/// One parameter's state after a cycle.
#[derive(Clone, Debug, PartialEq)]
pub struct ParamState {
    pub name: String,
    pub value: f64,
    /// The shift applied (damped, shortened to the maximum shifts, reflected or wrapped at a bound).
    pub shift: f64,
    pub esd: f64,
    /// shift/esd.
    pub shift_esd: f64,
    /// The shift the least squares computed, before damping (convergence is judged on it).
    pub computed: f64,
    pub at_limit: bool,
}

/// The agreement of a calculated pattern with the data (P43).
#[derive(Clone, Debug, PartialEq)]
pub struct Statistics {
    /// Points used (outside the excluded regions, with a positive variance), parameters, constraints.
    pub n: usize,
    pub p: usize,
    pub c: usize,
    pub rp: f64,
    pub rwp: f64,
    pub rexp: f64,
    /// χ² = S/(N − P + C).
    pub chi2: f64,
    /// The Durbin–Watson statistic of the weighted residuals (Hill & Flack 1987, as FullProf prints it, p. 44).
    pub durbin_watson: f64,
    /// Per phase: R_Bragg and R_F (%), with the observed intensities apportioned as Rietveld's (FullProf, p. 44).
    pub phases: Vec<(f64, f64)>,
    /// Rp, Rwp and Rexp with the background taken out of the denominators (yᵢ − bᵢ for yᵢ: FullProf's "conventional"
    /// R factors, manual p. 44; McCusker et al. 1999, §11), % (P2).
    pub rp_background_free: f64,
    pub rwp_background_free: f64,
    pub rexp_background_free: f64,
}

/// A restraint's state.
#[derive(Clone, Debug, PartialEq)]
pub struct RestraintState {
    pub label: String,
    pub target: f64,
    pub value: f64,
    pub sigma: f64,
    /// Its weight in the minimised sum (`restraint_states_at`).
    pub weight: f64,
    /// What it adds to χ² (`restraint_states_at`).
    pub contribution: f64,
}

/// A step tried and refused because it raised χ²: the λ it was solved with and the χ² it gave (W12).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Try {
    pub lambda: f64,
    pub chi2: f64,
}

/// What the steering did to a parameter's step in a cycle (E15c), as data: the kind, the parameter's address (the
/// dictionary's names), and the value it applies to (the shift's bound, the value reached, the λ and χ² of a refused
/// step). `Cycle::notes` keep the same in words.
#[derive(Clone, Debug, PartialEq)]
pub struct SteeringNote {
    pub kind: NoteKind,
    /// The free parameter's name (for a maximum shift or a limit) or the field's address (reflected, wrapped); none
    /// for a refused step.
    pub address: Option<String>,
    /// The maximum shift's field, the bound held at, the value reflected or wrapped to; none for a refused step.
    pub value: Option<f64>,
    /// A refused step's λ and χ².
    pub lambda: Option<f64>,
    pub chi2: Option<f64>,
}

/// The kinds of `SteeringNote`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NoteKind {
    /// The step shortened to a field's maximum shift.
    MaximumShift,
    /// The parameter held at a limit its step would cross.
    AtLimit,
    /// A field's value reflected at its bound.
    Reflected,
    /// A field's value wrapped round its period.
    Wrapped,
    /// A step refused because χ² rose; λ raised.
    Refused,
}

impl NoteKind {
    /// A stable word: maximum-shift, at-limit, reflected, wrapped, refused.
    pub fn word(self) -> &'static str {
        match self {
            NoteKind::MaximumShift => "maximum-shift",
            NoteKind::AtLimit => "at-limit",
            NoteKind::Reflected => "reflected",
            NoteKind::Wrapped => "wrapped",
            NoteKind::Refused => "refused",
        }
    }

    /// The kind written as `word`.
    pub fn from_word(w: &str) -> Option<NoteKind> {
        [NoteKind::MaximumShift, NoteKind::AtLimit, NoteKind::Reflected, NoteKind::Wrapped, NoteKind::Refused].into_iter().find(|k| k.word() == w)
    }
}

/// Why a run paused (E15a).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PauseWhy {
    /// A step of n cycles ended (`Session::step`).
    Step,
    /// The caller paused it (`Session::pause`: the user's Pause).
    Caller,
}

impl PauseWhy {
    /// "step" or "caller".
    pub fn word(self) -> &'static str {
        match self {
            PauseWhy::Step => "step",
            PauseWhy::Caller => "caller",
        }
    }
}

/// What was changed during a pause (E15a).
#[derive(Clone, Debug, PartialEq)]
pub enum Change {
    /// A value set by its address.
    Value { address: String, value: f64 },
    /// The free parameters changed.
    Parameters,
}

impl Change {
    pub fn words(&self) -> String {
        match self {
            Change::Value { address, value } => format!("{address} set to {value}"),
            Change::Parameters => "the free parameters changed".into(),
        }
    }
}

/// A pause of a run (E15a): after which cycle, why, λ at that moment, what was changed before the run went on and,
/// when something was, the λ the damping restarted at.
#[derive(Clone, Debug, PartialEq)]
pub struct Pause {
    pub after: usize,
    pub why: PauseWhy,
    pub lambda: f64,
    pub changes: Vec<Change>,
    pub restarted: Option<f64>,
}

/// How a run stopped and after which cycle (E15b): the stop's facts are its variant's fields.
#[derive(Clone, Debug, PartialEq)]
pub struct Stopped {
    pub stop: Stop,
    pub after: usize,
}

/// One cycle (Q5).
#[derive(Clone, Debug, PartialEq)]
pub struct Cycle {
    pub cycle: usize,
    pub params: Vec<ParamState>,
    /// Ycal at the cycle's values (empty in a session's `history`, which keeps the rest).
    pub ycalc: Vec<f64>,
    pub stats: Statistics,
    /// χ² at the start of the cycle.
    pub chi2_start: f64,
    /// λ, the Marquardt factor of the accepted step.
    pub marquardt: f64,
    /// The steps tried and refused before it, each with its λ and χ².
    pub tries: Vec<Try>,
    /// Directions left out as singular (complete correlations), when the singular guard is off.
    pub singular: usize,
    /// The restraints at the cycle's values, with their weights and contributions.
    pub restraints: Vec<RestraintState>,
    /// What the steering did: steps shortened, maximum shifts reached, values reflected or wrapped.
    pub notes: Vec<String>,
    /// The same as data, with each note's kind and address, and the parameters held at a limit (E15c); empty for a
    /// cycle read back from a file.
    pub steering: Vec<SteeringNote>,
    /// The covariance of the free parameters at the cycle's values, C χ² as `Refinement::covariance` (E1: the
    /// correlations at every cycle, `correlations_of`); empty for a cycle read back from a file.
    pub covariance: Vec<Vec<f64>>,
    /// The eigenvalues of the normal matrix scaled to unit diagonal at the start of the cycle, over the largest,
    /// largest first (E1, E5; `check::directions` gives the vectors at any state).
    pub spectrum: Vec<f64>,
    /// The time the cycle took, in seconds (E1).
    pub seconds: f64,
}

/// The covariance ellipse of two refined values (Explore's pair, L138): from their esd's σa, σb and correlation r, the
/// covariance C = [[σa², r σa σb], [r σa σb, σb²]], its eigenvalues λ₁ ≥ λ₂ and the angle of λ₁'s axis from a's axis
/// (radians, the 2 × 2 closed form). The kσ ellipse about the end has semi-axes k√λ₁ and k√λ₂ along that angle; at
/// k = 1 and 2 it is the contour Δχ² = 1 and 4 of the pair's marginal. None unless both esd's are positive.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Ellipse {
    pub l1: f64,
    pub l2: f64,
    pub angle: f64,
}

pub fn covariance_ellipse(sa: f64, sb: f64, r: f64) -> Option<Ellipse> {
    if !(sa > 0.0 && sb > 0.0) {
        return None;
    }
    let (cxx, cyy, cxy) = (sa * sa, sb * sb, r * sa * sb);
    let (tr, det) = (cxx + cyy, cxx * cyy - cxy * cxy);
    let disc = (tr * tr / 4.0 - det).max(0.0).sqrt();
    let (l1, l2) = (tr / 2.0 + disc, (tr / 2.0 - disc).max(0.0));
    let angle = if cxy.abs() < 1e-300 { if cxx >= cyy { 0.0 } else { std::f64::consts::FRAC_PI_2 } } else { (l1 - cxx).atan2(cxy) };
    Some(Ellipse { l1, l2, angle })
}

/// Correlations from a covariance (0 where a variance is 0).
pub fn correlations_of(cov: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = cov.len();
    (0..n).map(|j| (0..n).map(|k| if cov[j][j] > 0.0 && cov[k][k] > 0.0 { cov[j][k] / (cov[j][j] * cov[k][k]).sqrt() } else { 0.0 }).collect()).collect()
}

/// The path of a refinement (E1): the model it started from (links set), the free parameters, and every cycle run,
/// from which the model at any cycle is rebuilt (`state_at`). A session's `trajectory()`, or a results file's cycles
/// read back with the control file's model and parameters.
#[derive(Clone, Debug)]
pub struct Trajectory {
    pub start: Model,
    pub set: ParameterSet,
    pub cycles: Vec<Cycle>,
    /// The run's pauses (E15a), in order.
    pub pauses: Vec<Pause>,
    /// How it stopped and after which cycle (E15b); None while it can go on.
    pub stop: Option<Stopped>,
}

impl Trajectory {
    /// A path from its start, parameters and cycles, without pauses or stop.
    pub fn new(start: Model, set: ParameterSet, cycles: Vec<Cycle>) -> Trajectory {
        Trajectory { start, set, cycles, pauses: Vec::new(), stop: None }
    }

    /// The path of a run read back (E1): the model at its end (a results file's control file) and its cycles; the
    /// start is the end with every applied shift taken back, last cycle first. Floating-point subtraction does not undo
    /// the additions exactly, so the start rebuilt may differ from the run's own start in its last bits (Y₂O₃'s zero
    /// shift: 0.05 started, 0.05000000000000001 rebuilt on macOS, L187); which bits depends on the run's numbers, and
    /// those on the platform's libm (sin, cos, exp are not correctly rounded and differ between glibc and macOS). The
    /// same end and cycles always rebuild the same start.
    pub fn from_end(end: &Model, set: &ParameterSet, cycles: Vec<Cycle>) -> Trajectory {
        let mut m = end.clone();
        for c in cycles.iter().rev() {
            for q in &c.params {
                let Some(p) = set.free.iter().find(|p| p.name == q.name) else { continue };
                let d = if p.value.is_some() { q.shift } else { q.shift / p.fields[0].1 };
                for &(f, k) in &p.fields {
                    m.set(f, m.get(f) - k * d);
                }
            }
        }
        Trajectory::new(m, set.clone(), cycles)
    }

    /// The model after k cycles (0: the start): `state_at` over the first k.
    pub fn state_at(&self, k: usize) -> Model {
        state_at(&self.start, &self.set, &self.cycles[..k.min(self.cycles.len())])
    }
}

/// A path in units of the esd at its end (Explore's traces, DESIGN.md §2.1; L139: every number shown is the core's):
/// for each state, the start then each cycle, and each free parameter in the set's order, (value − value at the end)/
/// esd at the end (`from_end`) and the cycle's own shift/esd (`shift_esd`, none at the start); each parameter's
/// `reach`, its largest |distance| along the path; the value and esd at the end.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PathInEsd {
    pub names: Vec<String>,
    /// The state each row is: 0 the start, then each cycle's number.
    pub states: Vec<usize>,
    pub from_end: Vec<Vec<Option<f64>>>,
    pub shift_esd: Vec<Vec<Option<f64>>>,
    pub reach: Vec<f64>,
    pub end_value: Vec<f64>,
    pub end_esd: Vec<f64>,
}

impl Trajectory {
    /// The path in units of the esd at its end (`PathInEsd`). The end is the last cycle's values; its esd's are
    /// `end_esd` when given (a finished run's, from its covariance), else the last cycle's. A parameter without a
    /// positive esd at the end has no distances.
    pub fn in_esd(&self, end_esd: Option<&[f64]>) -> PathInEsd {
        let names: Vec<String> = self.set.free.iter().map(|p| p.name.clone()).collect();
        let at = |c: &Cycle, n: &str| c.params.iter().find(|q| q.name == n).cloned();
        let start: Vec<f64> = self.set.free.iter().map(|p| p.value_in(&self.start)).collect();
        let last = self.cycles.last();
        let end_value: Vec<f64> = names.iter().enumerate().map(|(j, n)| last.and_then(|c| at(c, n)).map_or(start[j], |q| q.value)).collect();
        let end_esd: Vec<f64> = match end_esd {
            Some(e) => e.to_vec(),
            None => names.iter().map(|n| last.and_then(|c| at(c, n)).map_or(0.0, |q| q.esd)).collect(),
        };
        let dist = |j: usize, v: f64| (end_esd.get(j).copied().unwrap_or(0.0) > 0.0).then(|| (v - end_value[j]) / end_esd[j]);
        let mut states = vec![0];
        let mut from_end = vec![(0..names.len()).map(|j| dist(j, start[j])).collect::<Vec<_>>()];
        let mut shift_esd = vec![vec![None; names.len()]];
        for c in &self.cycles {
            states.push(c.cycle);
            from_end.push(names.iter().enumerate().map(|(j, n)| at(c, n).and_then(|q| dist(j, q.value))).collect());
            shift_esd.push(names.iter().map(|n| at(c, n).and_then(|q| (q.esd > 0.0).then_some(q.shift_esd))).collect());
        }
        let reach = (0..names.len()).map(|j| from_end.iter().filter_map(|row| row[j]).fold(0.0_f64, |a, d| a.max(d.abs()))).collect();
        PathInEsd { names, states, from_end, shift_esd, reach, end_value, end_esd }
    }
}

/// The half-width in esd's of the χ² map that holds a pair's whole path (Explore's span "the path"): the larger reach
/// of the two, at least 4σ, with 15 % to spare.
pub fn pair_span(reach_a: f64, reach_b: f64) -> f64 {
    reach_a.max(reach_b).max(4.0) * 1.15
}

/// The mean time of the cycles run, in seconds (Refine's "s per cycle"; L146); none before a cycle.
pub fn mean_seconds(cycles: &[Cycle]) -> Option<f64> {
    (!cycles.is_empty()).then(|| cycles.iter().map(|c| c.seconds).sum::<f64>() / cycles.len() as f64)
}

impl RestraintState {
    /// How far the restraint is from its target, in its σ: (value − target)/σ (the restraints' table; L146); none
    /// without a positive σ.
    pub fn deviation(&self) -> Option<f64> {
        (self.sigma > 0.0).then(|| (self.value - self.target) / self.sigma)
    }
}

/// The share of a run done after its `cycle`-th cycle of `max_cycles` asked for, in [0, 1] (the Dock's gauge, §9;
/// L139): a resumed run counts on from its earlier cycles; no cycles asked for is done.
pub fn fraction_done(cycle: usize, max_cycles: usize) -> f64 {
    if max_cycles == 0 { 1.0 } else { (cycle as f64 / max_cycles as f64).min(1.0) }
}

impl Cycle {
    /// The largest |shift|/esd of the cycle and its parameter (the cycles' table, W11); none without a positive esd.
    pub fn largest_shift(&self) -> Option<(f64, &str)> {
        self.params.iter().filter(|p| p.esd > 0.0).map(|p| (p.shift_esd.abs(), p.name.as_str())).fold(None, |a: Option<(f64, &str)>, b| match a {
            Some(x) if x.0 >= b.0 => Some(x),
            _ => Some(b),
        })
    }
}

/// The model after the given cycles from `start` (E1): each cycle's applied shift of each parameter (its
/// `ParamState::shift` over its first multiplier, a variable's own) moved through every field it moves, by its
/// multiplier, as the engine applies a step. Parameters are matched by name, so a set changed during a pause is
/// followed. A Le Bail phase's extracted intensities are not rebuilt (the engine extracts them again each cycle).
pub fn state_at(start: &Model, set: &ParameterSet, cycles: &[Cycle]) -> Model {
    let mut m = start.clone();
    for c in cycles {
        for q in &c.params {
            let Some(p) = set.free.iter().find(|p| p.name == q.name) else { continue };
            let d = if p.value.is_some() { q.shift } else { q.shift / p.fields[0].1 };
            for &(f, k) in &p.fields {
                m.set(f, m.get(f) + k * d);
            }
        }
    }
    m
}

/// Why a refinement stopped (W11), as data: `code` is stable, `text` is the sentence shown.
#[derive(Clone, Debug, PartialEq)]
pub enum Stop {
    /// Every computed shift below `eps` times its esd.
    Converged { eps: f64 },
    /// No step, however short, lowers χ²: a minimum within the limits.
    Minimum,
    /// The cycles asked for were all run (Le Bail fits run them all).
    Cycles { n: usize, le_bail: bool },
    Diverging(Divergence),
    /// The normal matrix is singular (R24): the combination of free parameters the data cannot fix (name, component),
    /// its eigenvalue over the largest, and why and what to do, from the dictionary (`check::null_direction`).
    Singular { combination: Vec<(String, f64)>, ratio: f64, why: String, tip: String },
    /// Stopped by the user (or the caller).
    User,
}

/// How a refinement diverged (`steering::Guards`).
#[derive(Clone, Debug, PartialEq)]
pub enum Divergence {
    /// χ² rose `cycles` cycles in a row, from `from` (at the start of the first) to `to`.
    Rising { cycles: usize, from: f64, to: f64 },
    /// A parameter's computed shift was `ratio` esd's in cycle `cycle`.
    Jump { parameter: String, ratio: f64, cycle: usize },
    /// A value no structure or profile can have (`impossible`).
    Impossible { what: String, value: f64 },
}

impl Stop {
    /// A stable code: converged, minimum, cycles, diverging-rising, diverging-jump, diverging-impossible, singular,
    /// user.
    pub fn code(&self) -> &'static str {
        match self {
            Stop::Converged { .. } => "converged",
            Stop::Minimum => "minimum",
            Stop::Cycles { .. } => "cycles",
            Stop::Diverging(Divergence::Rising { .. }) => "diverging-rising",
            Stop::Diverging(Divergence::Jump { .. }) => "diverging-jump",
            Stop::Diverging(Divergence::Impossible { .. }) => "diverging-impossible",
            Stop::Singular { .. } => "singular",
            Stop::User => "user",
        }
    }

    /// The stop's facts (E15b), as (name, value) pairs: converged eps; cycles n and whether Le Bail; diverging-rising
    /// cycles, from, to; diverging-jump parameter, ratio, cycle; diverging-impossible what, value; singular ratio and
    /// the combination ("name coefficient; …").
    pub fn facts(&self) -> Vec<(&'static str, String)> {
        match self {
            Stop::Converged { eps } => vec![("eps", eps.to_string())],
            Stop::Minimum | Stop::User => Vec::new(),
            Stop::Cycles { n, le_bail } => vec![("count", n.to_string()), ("le_bail", if *le_bail { "yes" } else { "no" }.into())],
            Stop::Diverging(Divergence::Rising { cycles, from, to }) => vec![("cycles", cycles.to_string()), ("from", from.to_string()), ("to", to.to_string())],
            Stop::Diverging(Divergence::Jump { parameter, ratio, cycle }) => vec![("parameter", parameter.clone()), ("ratio", ratio.to_string()), ("cycle", cycle.to_string())],
            Stop::Diverging(Divergence::Impossible { what, value }) => vec![("what", what.clone()), ("value", value.to_string())],
            Stop::Singular { combination, ratio, .. } => vec![("ratio", ratio.to_string()), ("combination", combination.iter().map(|(n, c)| format!("{n} {c}")).collect::<Vec<_>>().join("; "))],
        }
    }

    /// A stop from its code and facts (`facts`, read back from a results file); the singular stop's why and tip are
    /// not kept. None for an unknown code or a missing fact.
    pub fn from_facts(code: &str, fact: impl Fn(&str) -> Option<String>) -> Option<Stop> {
        let num = |k: &str| fact(k).and_then(|v| v.parse::<f64>().ok());
        let int = |k: &str| fact(k).and_then(|v| v.parse::<usize>().ok());
        Some(match code {
            "converged" => Stop::Converged { eps: num("eps")? },
            "minimum" => Stop::Minimum,
            "user" => Stop::User,
            "cycles" => Stop::Cycles { n: int("count")?, le_bail: fact("le_bail").as_deref() == Some("yes") },
            "diverging-rising" => Stop::Diverging(Divergence::Rising { cycles: int("cycles")?, from: num("from")?, to: num("to")? }),
            "diverging-jump" => Stop::Diverging(Divergence::Jump { parameter: fact("parameter")?, ratio: num("ratio")?, cycle: int("cycle")? }),
            "diverging-impossible" => Stop::Diverging(Divergence::Impossible { what: fact("what")?, value: num("value")? }),
            "singular" => {
                let combination = fact("combination")?
                    .split("; ")
                    .filter(|t| !t.is_empty())
                    .map(|t| t.rsplit_once(' ').and_then(|(n, c)| c.parse().ok().map(|c| (n.to_string(), c))))
                    .collect::<Option<Vec<_>>>()?;
                Stop::Singular { combination, ratio: num("ratio")?, why: String::new(), tip: String::new() }
            }
            _ => return None,
        })
    }

    /// The sentence shown.
    pub fn text(&self) -> String {
        match self {
            Stop::Converged { eps } => format!("every shift below {eps} esd"),
            Stop::Minimum => "no step lowers χ² (a minimum within the limits)".into(),
            Stop::Cycles { n, le_bail: false } => format!("{n} cycles"),
            Stop::Cycles { n, le_bail: true } => format!("{n} cycles, all run (Le Bail)"),
            Stop::Diverging(Divergence::Rising { cycles, from, to }) => format!("diverging: χ² rose {cycles} cycles in a row, from {from:.4} to {to:.4}"),
            Stop::Diverging(Divergence::Jump { parameter, ratio, cycle }) => format!("diverging: {parameter} would move {ratio:.0} esd's in cycle {cycle}"),
            Stop::Diverging(Divergence::Impossible { what, value }) => format!("diverging: an impossible value, {what} = {value}"),
            Stop::Singular { combination, ratio, .. } => {
                let names: Vec<String> = combination.iter().map(|(n, c)| format!("{c:+.2}·{n}")).collect();
                format!("the matrix is singular: the data cannot fix the combination {} (eigenvalue {ratio:.1e} of the largest)", names.join(" "))
            }
            Stop::User => STOPPED.into(),
        }
    }

    /// Whether the run ended at a minimum of χ² (converged, or no step lowers it).
    pub fn converged(&self) -> bool {
        matches!(self, Stop::Converged { .. } | Stop::Minimum)
    }
}

/// A finished refinement.
#[derive(Clone, Debug)]
pub struct Refinement {
    pub cycles: usize,
    pub converged: bool,
    pub params: Vec<ParamState>,
    pub stats: Statistics,
    /// Correlation matrix of the free parameters.
    pub correlations: Vec<Vec<f64>>,
    /// Covariance of the free parameters, C χ² (the esd's are the roots of its diagonal, before the multipliers).
    pub covariance: Vec<Vec<f64>>,
    pub fields: Vec<(String, f64, f64)>,
    /// The restraints at the values reached: label, target, value, σ, weight and contribution to χ².
    pub restraints: Vec<RestraintState>,
    /// Why it stopped, as text (`reason.text()`).
    pub stop: String,
    /// Why it stopped, as data.
    pub reason: Stop,
    /// The largest |shift|/esd of the last cycle run (the CIF's `_refine_ls_shift/su_max`; 0 before any cycle).
    pub shift_su_max: f64,
    /// Every cycle run, without its Ycal (E1: the results file's `_drs_cycle_*` loops).
    pub history: Vec<Cycle>,
    /// The run's pauses (E15a: the results file's `_drs_pause_*` loop).
    pub pauses: Vec<Pause>,
    /// A plan run's (`run_plan`, P37): per cycle of `history`, the plan step that ran it, its number (1, 2, …) and
    /// name (the results file's `_drs_cycle_step` and `_drs_cycle_step_name`); empty for a run without a plan.
    pub cycle_steps: Vec<(usize, String)>,
    /// Bérar & Lelann's factor for local correlations (`berar_factor`, P43): the esd's times it estimate the probable
    /// errors.
    pub berar_factor: f64,
}

/// `Refinement::stop` of a refinement stopped by its caller (`refine_until`).
pub const STOPPED: &str = "stopped by the user";

/// Options of `refine`.
#[derive(Clone, Copy, Debug)]
pub struct Options {
    pub max_cycles: usize,
    /// Convergence: every |shift|/esd below eps.
    pub eps: f64,
}

/// The cycles of a run when the control file states none (L145; the app's Settings ▸ Default cycles replaces it).
pub const DEFAULT_CYCLES: usize = 10;
/// The convergence when the control file states none: every |shift| below 0.1 esd (FullProf's tutorial value).
pub const DEFAULT_EPS: f64 = 0.1;

impl Options {
    /// The options of a run from what the control file (or the caller) states, the defaults otherwise; at least one
    /// cycle.
    pub fn stated(cycles: Option<usize>, eps: Option<f64>) -> Options {
        Options { max_cycles: cycles.unwrap_or(DEFAULT_CYCLES).max(1), eps: eps.unwrap_or(DEFAULT_EPS) }
    }
}

/// |r| above which two parameters are reported as strongly correlated (the Refine panel's list, the heat map's
/// outlined cells; L145).
pub const STRONG_CORRELATION: f64 = 0.8;

/// The pairs (j, k), j < k, of a correlation matrix with |r| above `STRONG_CORRELATION`, the strongest first.
pub fn strong_correlations(r: &[Vec<f64>]) -> Vec<(usize, usize, f64)> {
    let mut out: Vec<(usize, usize, f64)> = (0..r.len()).flat_map(|j| (j + 1..r.len()).map(move |k| (j, k))).filter(|&(j, k)| r[j][k].abs() > STRONG_CORRELATION).map(|(j, k)| (j, k, r[j][k])).collect();
    out.sort_by(|a, b| b.2.abs().total_cmp(&a.2.abs()));
    out
}

/// The threads a parallel computation of the core uses (the χ² map): half the logical processors, at least one, so
/// that the computer stays usable while it runs (L150: a MacBook Pro with 16 logical processors was brought to a halt
/// when all 16 computed patterns at once).
pub fn worker_threads() -> usize {
    std::thread::available_parallelism().map_or(1, |t| t.get()).div_ceil(2).max(1)
}

/// The χ² range of a map's computed points for its colour bar: the smallest, the geometric middle (the bar is
/// logarithmic) and the largest (Explore's χ² map; L148). None without a positive value.
pub fn chi2_range(values: impl IntoIterator<Item = f64>) -> Option<(f64, f64, f64)> {
    let (lo, hi) = values.into_iter().filter(|v| v.is_finite() && *v > 0.0).fold((f64::INFINITY, 0.0_f64), |(lo, hi), v| (lo.min(v), hi.max(v)));
    (hi > 0.0).then(|| (lo, (lo * hi).sqrt(), hi))
}

/// The most strongly correlated pair (j, k), j < k, among the parameters `counted` (Explore's default pair: the
/// background points left out; L148).
pub fn strongest_pair(r: &[Vec<f64>], counted: impl Fn(usize) -> bool) -> Option<(usize, usize)> {
    (0..r.len()).flat_map(|j| (j + 1..r.len()).map(move |k| (j, k))).filter(|&(j, k)| counted(j) && counted(k)).max_by(|&(a, b), &(c, d)| r[a][b].abs().total_cmp(&r[c][d].abs()))
}

/// The half-widths, in the parameters' units, of a χ² map `span_esd` esd's wide on each side (Explore's map; L145).
pub fn half_widths(span_esd: f64, esd_a: f64, esd_b: f64) -> (f64, f64) {
    (span_esd * esd_a, span_esd * esd_b)
}

/// The weight of each point in Σw(y − yc)² for a calculated pattern: 0 for a point left out (in an excluded region,
/// by `Model::excluded_points`, or without a positive variance), else 1/variance under the model's `Weighting`.
pub fn weights(model: &Model, data: &Pattern, ycalc: &[f64]) -> Vec<f64> {
    let excluded = model.excluded_points(&data.x);
    (0..data.x.len())
        .map(|i| {
            if excluded[i] {
                return 0.0;
            }
            let (y, s2) = (data.y[i], data.sigma[i] * data.sigma[i]);
            let var = match model.weighting {
                Weighting::Variance => s2,
                Weighting::MaximumLikelihood => ycalc[i] * if y > 0.0 && s2 > 0.0 { s2 / y } else { 1.0 },
                Weighting::Unit => 1.0,
            };
            if var > 0.0 {
                1.0 / var
            } else {
                0.0
            }
        })
        .collect()
}

pub(crate) fn objective(data: &Pattern, ycalc: &[f64], w: &[f64]) -> f64 {
    (0..w.len()).map(|i| w[i] * (data.y[i] - ycalc[i]).powi(2)).sum()
}

/// Bérar & Lelann's factor for local correlations (J. Appl. Cryst. 24, 1–5, 1991, in `manuals/papers/`; P43, L244):
/// the esd's of a refinement whose residuals are serially correlated, multiplied by it, estimate the probable errors.
/// With aᵢ = √wᵢ (yᵢ − ycᵢ), the residuals of correlated neighbours add linearly, the others quadratically:
/// S″ = Σ(1 − zᵢ²)aᵢ² + Σ_groups (Σ zᵢ|aᵢ|)² (eq. 10), factor √(S″/S) (eq. 11, S = Σaᵢ²), where a point and its
/// neighbour of the same sign are correlated with zᵢ = q/(2 + q), q = √(2π(aᵢ² + a²_next)) (eq. 12 with the appendix's
/// 2π; eq. 12 writes 2), zero for opposite signs, and a group is a run of correlated points. FullProf's reading, which
/// reproduces its printed factor to its four decimals in `berar_factor` (L244):
/// - each pair's z is given to its first point at the low angle side, as the appendix's loop does (eq. 10 says "with
///   its predecessor"), the loop running from the high-angle end;
/// - only the points with Bragg contributions count (inside a peak's window; the appendix skips points without a
///   reflection, `KR(I) = 0`), and a point left out breaks a run.
pub fn berar_factor(model: &Model, data: &Pattern, ycalc: &[f64], peaks: &[crate::powder::Peak]) -> f64 {
    let w = weights(model, data, ycalc);
    let bragg = bragg_points(data, peaks);
    let a: Vec<Option<f64>> = (0..w.len()).map(|i| (w[i] > 0.0 && bragg[i]).then(|| w[i].sqrt() * (data.y[i] - ycalc[i]))).collect();
    let (mut s2, mut s, mut group) = (0.0, 0.0, 0.0);
    for i in (0..a.len()).rev() {
        let Some(x) = a[i] else {
            s2 += group * group;
            group = 0.0;
            continue;
        };
        s += x * x;
        let z = match a.get(i + 1).copied().flatten() {
            Some(y) if x * y > 0.0 => {
                let q = (2.0 * std::f64::consts::PI * (x * x + y * y)).sqrt();
                q / (2.0 + q)
            }
            _ => 0.0,
        };
        if z == 0.0 {
            s2 += group * group;
            group = 0.0;
        }
        s2 += (1.0 - z * z) * x * x;
        group += z * x.abs();
    }
    s2 += group * group;
    if s > 0.0 {
        (s2 / s).sqrt()
    } else {
        1.0
    }
}

/// The statistics of a calculated pattern with p free parameters (P43).
pub fn statistics(model: &Model, data: &Pattern, ycalc: &[f64], p: usize) -> Statistics {
    statistics_over(model, data, ycalc, p, None)
}

/// The points with Bragg contributions: inside some peak's window (FullProf's second block of statistics, its global
/// χ² "(Bragg contrib.)" and the Bérar factor's points).
pub fn bragg_points(data: &Pattern, peaks: &[crate::powder::Peak]) -> Vec<bool> {
    let mut bragg = vec![false; data.x.len()];
    for pk in peaks {
        let lo = data.x.partition_point(|&t| t < pk.two_theta - pk.reach);
        let hi = data.x.partition_point(|&t| t <= pk.two_theta + pk.reach);
        bragg[lo..hi.max(lo)].iter_mut().for_each(|b| *b = true);
    }
    bragg
}

/// `statistics` over the points of `mask` only (all the points used when None); R_Bragg and R_F are the phases'
/// either way.
pub fn statistics_over(model: &Model, data: &Pattern, ycalc: &[f64], p: usize, mask: Option<&[bool]>) -> Statistics {
    let mut w = weights(model, data, ycalc);
    if let Some(m) = mask {
        w.iter_mut().zip(m).filter(|(_, &keep)| !keep).for_each(|(wi, _)| *wi = 0.0);
    }
    let c = model.calculate(&data.x);
    let (mut sa, mut sy, mut swd, mut swy, mut n) = (0.0, 0.0, 0.0, 0.0, 0usize);
    let (mut syb, mut swyb) = (0.0, 0.0);
    let (mut dw_num, mut prev): (f64, Option<f64>) = (0.0, None);
    for i in 0..w.len() {
        if w[i] == 0.0 {
            continue;
        }
        let (y, d) = (data.y[i], data.y[i] - ycalc[i]);
        sa += d.abs();
        sy += y.abs();
        swd += w[i] * d * d;
        swy += w[i] * y * y;
        // The net counts summed with their signs, as the manual's yᵢ − bᵢ reads and FullProf sums them (LaB6: 26.7
        // with the signs, 23.9 with |yᵢ − bᵢ|; P2, L114).
        let net = y - c.background[i];
        syb += net;
        swyb += w[i] * net * net;
        n += 1;
        let r = w[i].sqrt() * d;
        if let Some(q) = prev {
            dw_num += (r - q) * (r - q);
        }
        prev = Some(r);
    }
    let dof = (n.saturating_sub(p)).max(1) as f64;
    // R_Bragg and R_F: I_obs,h = I_calc,h Σ_i Ω_h,i (y_i − b_i)/(y_c,i − b_i).
    let excluded = model.excluded_points(&data.x);
    let mut phases = vec![(0.0, 0.0, 0.0, 0.0); model.phases.len()];
    for pk in &c.peaks {
        // Both sums are over the peak's points, so that they are in the same units (counts).
        let (mut iobs, mut icalc) = (0.0, 0.0);
        model.peak_profile(&data.x, &excluded, pk, |i, v| {
            let net = c.ycalc[i] - c.background[i];
            if net > 0.0 && w[i] > 0.0 {
                iobs += v * (data.y[i] - c.background[i]) / net;
                icalc += v;
            }
        });
        let e = &mut phases[pk.phase];
        e.0 += (iobs - icalc).abs();
        e.1 += iobs.abs();
        e.2 += (iobs.max(0.0).sqrt() - icalc.max(0.0).sqrt()).abs();
        e.3 += iobs.max(0.0).sqrt();
    }
    Statistics {
        n,
        p,
        c: 0,
        rp: 100.0 * sa / sy,
        rwp: 100.0 * (swd / swy).sqrt(),
        rexp: 100.0 * (dof / swy).sqrt(),
        chi2: swd / dof,
        durbin_watson: if swd > 0.0 { dw_num / swd } else { 0.0 },
        phases: phases.iter().map(|e| (100.0 * e.0 / e.1.max(1e-300), 100.0 * e.2 / e.3.max(1e-300))).collect(),
        // Σ(y − b) as it is (manual p. 44: "in the denominators the quantity yᵢ is changed to yᵢ − bᵢ"): negative where the
        // background lies above the counts, as FullProf prints it (y2o3_poly12_b7.sum, cRp −145; L353).
        rp_background_free: if syb == 0.0 { 0.0 } else { 100.0 * sa / syb },
        rwp_background_free: 100.0 * (swd / swyb.max(1e-300)).sqrt(),
        rexp_background_free: 100.0 * (dof / swyb.max(1e-300)).sqrt(),
    }
}

/// The restraints of a model and a parameter set: (target − value, σ, ∂value/∂parameter).
fn restraint_rows(model: &Model, set: &ParameterSet) -> Vec<(f64, f64, Vec<f64>)> {
    let np = set.free.len();
    let mut out = Vec::new();
    for r in &model.restraints {
        let (v, df) = r.derivatives(model);
        let mut gr = vec![0.0; np];
        for (j, p) in set.free.iter().enumerate() {
            for &(f, m) in &p.fields {
                if let Some(&(_, d)) = df.iter().find(|q| q.0 == f) {
                    gr[j] += m * d;
                }
            }
        }
        // FullProf's angle restraints pull with π/180 of their stated weight (`Convention::angle_restraints`, L058).
        let k = if model.convention.angle_restraints && matches!(r.kind, crate::restraints::Kind::Angle { .. }) { std::f64::consts::PI / 180.0 } else { 1.0 };
        out.push((r.target - v, r.sigma / k.sqrt(), gr));
    }
    for l in &set.linear {
        let (v, d) = l.value(model, set);
        out.push((l.target - v, l.sigma, d));
    }
    out
}

/// Σ ((target − value)/σ)².
fn penalty(rows: &[(f64, f64, Vec<f64>)]) -> f64 {
    rows.iter().map(|r| (r.0 / r.1).powi(2)).sum()
}

/// The restraints' rows added to the normal equations with weight c/σ².
fn add_restraint_rows(a: &mut [Vec<f64>], g: &mut [f64], rows: &[(f64, f64, Vec<f64>)], c: f64) {
    for (res, sigma, gr) in rows {
        let wt = c / (sigma * sigma);
        for j in 0..gr.len() {
            if gr[j] == 0.0 {
                continue;
            }
            g[j] += wt * res * gr[j];
            for k in 0..gr.len() {
                a[j][k] += wt * gr[j] * gr[k];
            }
        }
    }
}

/// The restraints' labels, targets, values and σ, each with its weight 1/σ² (as pulled: FullProf's angle restraints
/// with π/180 of it) and its term ((target − value)/σ)² of the penalty (`restraint_states_at` scales both).
pub fn restraint_states(model: &Model, set: &ParameterSet) -> Vec<RestraintState> {
    restraint_states_at(model, set, 1.0, 1)
}

/// The restraints' states with their weight in the minimised sum, c/σ² (c the data's χ²_ν, as the engine weights them,
/// P40), and their contribution to χ²: c((target − value)/σ)²/(N − P), what each adds to the minimised sum over the
/// degrees of freedom (`dof` = N − P).
pub fn restraint_states_at(model: &Model, set: &ParameterSet, c: f64, dof: usize) -> Vec<RestraintState> {
    let rows = restraint_rows(model, set);
    let mut out: Vec<RestraintState> = model.restraints.iter().map(|r| RestraintState { label: r.label(model), target: r.target, value: r.value(model), sigma: r.sigma, weight: 0.0, contribution: 0.0 }).collect();
    for l in &set.linear {
        out.push(RestraintState { label: l.name.clone(), target: l.target, value: l.value(model, set).0, sigma: l.sigma, weight: 0.0, contribution: 0.0 });
    }
    for (st, (res, sigma, _)) in out.iter_mut().zip(&rows) {
        st.weight = c / (sigma * sigma);
        st.contribution = st.weight * res * res / dof.max(1) as f64;
    }
    out
}

/// A parameter's path over the cycles: (cycle, value, esd) per cycle.
pub type Path = Vec<(usize, f64, f64)>;

/// The restraints' part of χ², the sum of their contributions (the app's request 27).
pub fn restraints_part(states: &[RestraintState]) -> f64 {
    states.iter().map(|r| r.contribution).sum()
}

/// Fill the lower triangle of a matrix from its upper triangle.
pub(crate) fn symmetrize(a: &mut [Vec<f64>]) {
    for j in 1..a.len() {
        let (top, rest) = a.split_at_mut(j);
        for (k, row) in top.iter().enumerate() {
            rest[0][k] = row[j];
        }
    }
}

/// The size above which `symmetric_eigen` reduces to tridiagonal form instead of rotating (P36): Jacobi's sweeps
/// cost of order n³ each, and a Pawley fit's 400 intensities took half a minute a cycle. Below it, where every
/// comparison with FullProf was made and frozen (L037), the rotations are kept, to the last bit.
pub const JACOBI_LIMIT: usize = 120;

/// Eigen-decomposition of a symmetric matrix: (eigenvalues, eigenvectors as columns). Cyclic Jacobi rotations up to
/// `JACOBI_LIMIT` (the eigenvalues in the diagonal's order), Householder's reduction and the implicit QL iteration
/// above it (`householder_ql`, the eigenvalues ascending).
pub fn symmetric_eigen(a: &[Vec<f64>]) -> (Vec<f64>, Vec<Vec<f64>>) {
    if a.len() > JACOBI_LIMIT {
        return householder_ql(a);
    }
    jacobi(a)
}

/// Householder's reduction to tridiagonal form and the implicit QL iteration with shifts, as EISPACK's tred2 and
/// tql2 (Bowdler, Martin, Reinsch & Wilkinson 1968; Smith et al. 1976) in JAMA's port (public domain): eigenvalues
/// ascending, eigenvectors as columns.
pub fn householder_ql(a: &[Vec<f64>]) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = a.len();
    if n == 0 {
        return (Vec::new(), Vec::new());
    }
    let mut v = a.to_vec();
    let mut d = vec![0.0; n];
    let mut e = vec![0.0; n];
    // tred2.
    d.copy_from_slice(&v[n - 1]);
    for i in (1..n).rev() {
        let mut scale = 0.0;
        let mut h = 0.0;
        for dk in d.iter().take(i) {
            scale += dk.abs();
        }
        if scale == 0.0 {
            e[i] = d[i - 1];
            for j in 0..i {
                d[j] = v[i - 1][j];
                v[i][j] = 0.0;
                v[j][i] = 0.0;
            }
        } else {
            for dk in d.iter_mut().take(i) {
                *dk /= scale;
                h += *dk * *dk;
            }
            let f = d[i - 1];
            let mut g = h.sqrt();
            if f > 0.0 {
                g = -g;
            }
            e[i] = scale * g;
            h -= f * g;
            d[i - 1] = f - g;
            for ej in e.iter_mut().take(i) {
                *ej = 0.0;
            }
            for j in 0..i {
                let f = d[j];
                v[j][i] = f;
                let mut g = e[j] + v[j][j] * f;
                for k in j + 1..i {
                    g += v[k][j] * d[k];
                    e[k] += v[k][j] * f;
                }
                e[j] = g;
            }
            let mut f = 0.0;
            for j in 0..i {
                e[j] /= h;
                f += e[j] * d[j];
            }
            let hh = f / (h + h);
            for j in 0..i {
                e[j] -= hh * d[j];
            }
            for j in 0..i {
                let (f, g) = (d[j], e[j]);
                for k in j..i {
                    v[k][j] -= f * e[k] + g * d[k];
                }
                d[j] = v[i - 1][j];
                v[i][j] = 0.0;
            }
        }
        d[i] = h;
    }
    for i in 0..n - 1 {
        v[n - 1][i] = v[i][i];
        v[i][i] = 1.0;
        let h = d[i + 1];
        if h != 0.0 {
            for k in 0..=i {
                d[k] = v[k][i + 1] / h;
            }
            for j in 0..=i {
                let g: f64 = v.iter().take(i + 1).map(|row| row[i + 1] * row[j]).sum();
                for (row, dk) in v.iter_mut().zip(&d).take(i + 1) {
                    row[j] -= g * dk;
                }
            }
        }
        for row in v.iter_mut().take(i + 1) {
            row[i + 1] = 0.0;
        }
    }
    for j in 0..n {
        d[j] = v[n - 1][j];
        v[n - 1][j] = 0.0;
    }
    v[n - 1][n - 1] = 1.0;
    e[0] = 0.0;
    // tql2.
    for i in 1..n {
        e[i - 1] = e[i];
    }
    e[n - 1] = 0.0;
    let (mut f, mut tst1) = (0.0f64, 0.0f64);
    let eps = f64::EPSILON;
    for l in 0..n {
        tst1 = tst1.max(d[l].abs() + e[l].abs());
        let mut m = l;
        while m < n {
            if e[m].abs() <= eps * tst1 {
                break;
            }
            m += 1;
        }
        if m > l {
            for _ in 0..200 {
                let mut g = d[l];
                let mut p = (d[l + 1] - g) / (2.0 * e[l]);
                let mut r = p.hypot(1.0);
                if p < 0.0 {
                    r = -r;
                }
                d[l] = e[l] / (p + r);
                d[l + 1] = e[l] * (p + r);
                let dl1 = d[l + 1];
                let mut h = g - d[l];
                for di in d.iter_mut().skip(l + 2) {
                    *di -= h;
                }
                f += h;
                p = d[m];
                let (mut c, mut c2, mut c3) = (1.0, 1.0, 1.0);
                let el1 = e[l + 1];
                let (mut s, mut s2) = (0.0, 0.0);
                for i in (l..m).rev() {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[i];
                    h = c * p;
                    r = p.hypot(e[i]);
                    e[i + 1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = c * d[i] - s * g;
                    d[i + 1] = h + s * (c * g + s * d[i]);
                    for row in v.iter_mut() {
                        let h = row[i + 1];
                        row[i + 1] = s * row[i] + c * h;
                        row[i] = c * row[i] - s * h;
                    }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;
                if e[l].abs() <= eps * tst1 {
                    break;
                }
            }
        }
        d[l] += f;
        e[l] = 0.0;
    }
    // Ascending, with the vectors.
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&x, &y| d[x].total_cmp(&d[y]));
    let lam = order.iter().map(|&k| d[k]).collect();
    let vec = v.iter().map(|row| order.iter().map(|&k| row[k]).collect()).collect();
    (lam, vec)
}

/// Eigen-decomposition of a symmetric matrix by cyclic Jacobi rotations: (eigenvalues, eigenvectors as columns).
fn jacobi(a: &[Vec<f64>]) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = a.len();
    let mut m = a.to_vec();
    let mut v = vec![vec![0.0; n]; n];
    for (i, r) in v.iter_mut().enumerate() {
        r[i] = 1.0;
    }
    for _ in 0..100 {
        let off: f64 = (0..n).flat_map(|i| (0..n).filter(move |&j| j != i).map(move |j| (i, j))).map(|(i, j)| m[i][j] * m[i][j]).sum();
        let diag: f64 = (0..n).map(|i| m[i][i] * m[i][i]).sum();
        if off <= 1e-30 * diag.max(1e-300) {
            break;
        }
        for p in 0..n {
            for q in p + 1..n {
                if m[p][q].abs() < 1e-300 {
                    continue;
                }
                let theta = (m[q][q] - m[p][p]) / (2.0 * m[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for row in m.iter_mut() {
                    let (mkp, mkq) = (row[p], row[q]);
                    row[p] = c * mkp - s * mkq;
                    row[q] = s * mkp + c * mkq;
                }
                let (rp, rq) = (m[p].clone(), m[q].clone());
                for (k, (mpk, mqk)) in rp.iter().zip(&rq).enumerate() {
                    m[p][k] = c * mpk - s * mqk;
                    m[q][k] = s * mpk + c * mqk;
                }
                for row in v.iter_mut() {
                    let (vkp, vkq) = (row[p], row[q]);
                    row[p] = c * vkp - s * vkq;
                    row[q] = s * vkp + c * vkq;
                }
            }
        }
    }
    ((0..n).map(|i| m[i][i]).collect(), v)
}

/// The scaled normal matrix's decomposition and the step for a Marquardt factor, over the parameters `active`.
pub(crate) struct Solver {
    d: Vec<f64>,
    lam: Vec<f64>,
    vec: Vec<Vec<f64>>,
    idx: Vec<usize>,
    singular: usize,
}

impl Solver {
    /// The null directions (eigenvalue below 10⁻¹² of the largest), each as (eigenvalue over the largest, components
    /// on the free parameters: index into the full set, value).
    fn null_directions(&self) -> Vec<(f64, Vec<(usize, f64)>)> {
        let top = self.lam.iter().copied().fold(0.0, f64::max);
        (0..self.idx.len()).filter(|&k| self.lam[k] < 1e-12 * top).map(|k| (self.lam[k] / top, (0..self.idx.len()).map(|i| (self.idx[i], self.vec[i][k])).collect())).collect()
    }

    pub(crate) fn new(a: &[Vec<f64>], idx: Vec<usize>) -> Solver {
        let d: Vec<f64> = idx.iter().map(|&j| a[j][j].sqrt()).collect();
        let n = idx.len();
        let s: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|k| a[idx[i]][idx[k]] / (d[i] * d[k])).collect()).collect();
        let (lam, vec) = symmetric_eigen(&s);
        let top = lam.iter().copied().fold(0.0, f64::max);
        let singular = lam.iter().filter(|&&l| l < 1e-12 * top).count();
        Solver { d, lam, vec, idx, singular }
    }

    /// (A + μ D²)⁻¹ applied to g, as a full-length vector.
    pub(crate) fn solve(&self, g: &[f64], mu: f64, len: usize) -> Vec<f64> {
        let n = self.idx.len();
        let top = self.lam.iter().copied().fold(0.0, f64::max);
        let gs: Vec<f64> = (0..n).map(|i| g[self.idx[i]] / self.d[i]).collect();
        let mut out = vec![0.0; len];
        for k in 0..n {
            if self.lam[k] < 1e-12 * top {
                continue;
            }
            let proj: f64 = (0..n).map(|i| self.vec[i][k] * gs[i]).sum::<f64>() / (self.lam[k] + mu);
            for i in 0..n {
                out[self.idx[i]] += self.vec[i][k] * proj / self.d[i];
            }
        }
        out
    }

    /// The pseudo-inverse of A (zero outside `idx`).
    pub(crate) fn inverse(&self, len: usize) -> Vec<Vec<f64>> {
        let n = self.idx.len();
        let top = self.lam.iter().copied().fold(0.0, f64::max);
        let mut inv = vec![vec![0.0; len]; len];
        for i in 0..n {
            for j in 0..n {
                let mut s = 0.0;
                for k in 0..n {
                    if self.lam[k] >= 1e-12 * top {
                        s += self.vec[i][k] * self.vec[j][k] / self.lam[k];
                    }
                }
                inv[self.idx[i]][self.idx[j]] = s / (self.d[i] * self.d[j]);
            }
        }
        inv
    }
}

/// The derivative matrix of the free parameters (one row per parameter) and the calculated pattern.
pub(crate) fn jacobian(model: &Model, data: &Pattern, set: &ParameterSet, fields: &[Field]) -> (Vec<f64>, Vec<Vec<f64>>) {
    let (c, jf) = model.derivatives(&data.x, fields);
    let jq = set
        .free
        .iter()
        .map(|p| {
            let mut col = vec![0.0; data.x.len()];
            for &(f, m) in &p.fields {
                let k = fields.iter().position(|g| *g == f).unwrap();
                for (a, b) in col.iter_mut().zip(&jf[k]) {
                    *a += m * b;
                }
            }
            col
        })
        .collect();
    (c.ycalc, jq)
}

/// The normal matrix the engine would invert at the model's values (data weighted as the engine weighs them, the
/// restraints' rows at χ²_ν), and each free parameter's derivative column: what `check` reads before a run (R24, L076).
pub fn normal_matrix(model: &Model, data: &Pattern, set: &ParameterSet) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let fields = set.fields();
    let np = set.free.len();
    let (yc, jq) = jacobian(model, data, set, &fields);
    let w = weights(model, data, &yc);
    let mut a = vec![vec![0.0; np]; np];
    for i in 0..w.len() {
        if w[i] == 0.0 {
            continue;
        }
        for j in 0..np {
            for k in j..np {
                a[j][k] += w[i] * jq[j][i] * jq[k][i];
            }
        }
    }
    symmetrize(&mut a);
    let stats = statistics(model, data, &yc, np);
    add_restraint_rows(&mut a, &mut vec![0.0; np], &restraint_rows(model, set), stats.chi2);
    // The columns as the weights see them (a column zero at every weighted point moves nothing that is fitted).
    let cols = jq.into_iter().map(|c| c.iter().zip(&w).map(|(d, wi)| d * wi.sqrt()).collect()).collect();
    (a, cols)
}

fn apply(model: &mut Model, set: &ParameterSet, delta: &[f64]) {
    for (p, &d) in set.free.iter().zip(delta) {
        for &(f, m) in &p.fields {
            model.set(f, model.get(f) + m * d);
        }
    }
}

/// How many free parameters move each field.
fn movers(set: &ParameterSet) -> std::collections::BTreeMap<Field, usize> {
    let mut count = std::collections::BTreeMap::new();
    for p in &set.free {
        for &(f, _) in &p.fields {
            *count.entry(f).or_insert(0) += 1;
        }
    }
    count
}

/// The range of the shift of each parameter allowed by the limits of the fields that only it moves, where the field
/// stops at its bounds (`Behaviour::Stop`; reflected and periodic bounds are applied after the step).
fn shift_bounds(model: &Model, set: &ParameterSet) -> Vec<(f64, f64)> {
    let count = movers(set);
    set.free
        .iter()
        .map(|p| {
            let (mut lo, mut hi) = (f64::NEG_INFINITY, f64::INFINITY);
            for &(f, m) in &p.fields {
                if count[&f] != 1 || m == 0.0 || set.behaviour_of(f) != Behaviour::Stop {
                    continue;
                }
                let (fl, fh) = set.limits_of(f, model);
                let v = model.get(f);
                for (bound, is_low) in [(fl, true), (fh, false)] {
                    let Some(b) = bound else { continue };
                    let s = (b - v) / m;
                    if (is_low && m > 0.0) || (!is_low && m < 0.0) {
                        lo = lo.max(s);
                    } else {
                        hi = hi.min(s);
                    }
                }
            }
            (lo.min(0.0), hi.max(0.0))
        })
        .collect()
}

/// A value of the model no structure or profile can have (the divergence guard, W11): a value that is not a finite
/// number, a cell length not above 0 or an angle outside (0°, 180°) or a cell of no volume, or a negative Gaussian H²
/// at a reflection (a width no peak can have; FullProf stops on it). Returns what and the value.
pub fn impossible(model: &Model, calc: &Calculated) -> Option<(String, f64)> {
    for f in model.all_fields() {
        let v = model.get(f);
        if !v.is_finite() {
            return Some((crate::refinable::address(f, model), v));
        }
    }
    for (p, ph) in model.phases.iter().enumerate() {
        let c = &ph.cell;
        for (k, v) in [c.a, c.b, c.c].into_iter().enumerate() {
            if v <= 0.0 {
                return Some((crate::refinable::address(Field::Cell(p, k), model), v));
            }
        }
        for (k, v) in [c.alpha, c.beta, c.gamma].into_iter().enumerate() {
            if v <= 0.0 || v >= 180.0 {
                return Some((crate::refinable::address(Field::Cell(p, 3 + k), model), v));
            }
        }
        let vol = c.volume();
        if !(vol.is_finite() && vol > 0.0) {
            return Some((format!("phase{}'s cell volume", p + 1), vol));
        }
    }
    for pk in &calc.peaks {
        let h2 = model.gaussian_h2(&model.phases[pk.phase], pk.two_theta, pk.wavelength);
        if h2 < 0.0 {
            return Some((format!("phase{}'s Gaussian H² at 2θ = {:.3}° ({} {} {})", pk.phase + 1, pk.two_theta, pk.hkl[0], pk.hkl[1], pk.hkl[2]), h2));
        }
    }
    if calc.ycalc.iter().any(|y| !y.is_finite()) {
        return Some(("Ycal".into(), f64::NAN));
    }
    None
}

/// FullProf's printed esd of free parameter j (its eq. 3.7, p. 43: σ·a, a the codeword's multiplier, which the core
/// imports as the parameter's own damping): the FullProf-compatibility reading of an esd, for comparisons. The core's
/// esd is σ.
pub fn fullprof_esd(set: &ParameterSet, j: usize, esd: f64) -> f64 {
    let p = &set.free[j];
    esd * set.steering.damping.iter().find(|d| d.0 == p.name).map_or(1.0, |d| d.1)
}

/// A refinement run cycle by cycle (W11): `step` runs n cycles and pauses, `run` goes on to the end; between calls
/// the session keeps its state (the model at the cycle's values, λ, the cycle count, the history of every cycle, the
/// count of cycles χ² has risen), so a paused run resumes where it was. Values or free parameters changed during a
/// pause (`set_value`, `set_parameters`) restart the damping: λ goes back to its start, and the notes say so.
/// `to_text` and `resume` save and restore the state for scripts (`drs refine --session`).
#[derive(Clone, Debug)]
pub struct Session {
    pub model: Model,
    pub data: Pattern,
    pub set: ParameterSet,
    pub opts: Options,
    /// The Marquardt factor λ the next cycle starts from.
    pub lambda: f64,
    /// The cycles run so far.
    pub cycle: usize,
    /// Every cycle run, without its Ycal.
    pub history: Vec<Cycle>,
    /// Why it stopped; None while it can go on (paused, or not started).
    pub stop: Option<Stop>,
    /// What the session did outside the cycles: links set at the start, λ restarted after a change.
    pub notes: Vec<String>,
    /// The model the session started from, links set (E1: `trajectory`).
    pub start: Model,
    /// The run's pauses (E15a): each step's end and each caller's pause, with what was changed before going on.
    pub pauses: Vec<Pause>,
    rising: usize,
    rising_from: f64,
    /// Values were set from outside since the last cycle: the next cycle checks them first.
    edited: bool,
}

impl Session {
    /// A session on a model, its data and the parameters to refine; the links' values are set from the values they
    /// follow (`ParameterSet::enforce_links`).
    pub fn new(mut model: Model, data: Pattern, set: ParameterSet, opts: Options) -> Session {
        let mut notes = Vec::new();
        for (f, old, new) in set.enforce_links(&mut model) {
            notes.push(format!("{} set from {old} to {new} by its link", crate::refinable::address(f, &model)));
        }
        let lambda = set.steering.lambda;
        Session { start: model.clone(), model, data, set, opts, lambda, cycle: 0, history: Vec::new(), stop: None, notes, pauses: Vec::new(), rising: 0, rising_from: 0.0, edited: true }
    }

    /// The path so far (E1): the start, the parameters, every cycle run.
    pub fn trajectory(&self) -> Trajectory {
        Trajectory { start: self.start.clone(), set: self.set.clone(), cycles: self.history.clone(), pauses: self.pauses.clone(), stop: self.stop.clone().map(|stop| Stopped { stop, after: self.cycle }) }
    }

    /// Run up to n cycles, calling `on_cycle` after each (false stops the run as the user's Stop does), then pause;
    /// returns the stop if the run ended. The cycles asked for (`opts.max_cycles`) do not limit a step.
    pub fn step(&mut self, n: usize, mut on_cycle: impl FnMut(&Cycle) -> bool) -> Option<&Stop> {
        for _ in 0..n {
            if self.stop.is_some() {
                break;
            }
            self.one_cycle(&mut on_cycle);
        }
        // The step's end is a pause (E15a), unless the run ended.
        if self.stop.is_none() {
            self.mark_pause(PauseWhy::Step);
        }
        self.stop.as_ref()
    }

    /// The caller's pause (E15a): recorded after the cycle run last, in place of the Stop a caller's `on_cycle` gives
    /// (which is taken back: the run can go on).
    pub fn pause(&mut self) {
        if self.stop == Some(Stop::User) {
            self.stop = None;
        }
        self.mark_pause(PauseWhy::Caller);
    }

    /// The pause after the current cycle, made the caller's when there is one already (a step's end paused by the user).
    fn mark_pause(&mut self, why: PauseWhy) {
        match self.pauses.last_mut() {
            Some(p) if p.after == self.cycle => {
                if why == PauseWhy::Caller {
                    p.why = why;
                }
            }
            _ => self.pauses.push(Pause { after: self.cycle, why, lambda: self.lambda, changes: Vec::new(), restarted: None }),
        }
    }

    /// Run until the run ends: converged, diverging, singular, stopped by `on_cycle`, or `opts.max_cycles` cycles run
    /// in all.
    pub fn run(&mut self, mut on_cycle: impl FnMut(&Cycle) -> bool) -> &Stop {
        while self.stop.is_none() {
            if self.cycle >= self.opts.max_cycles {
                self.stop = Some(Stop::Cycles { n: self.cycle, le_bail: self.model.has_extracted() });
                break;
            }
            self.one_cycle(&mut on_cycle);
        }
        self.stop.as_ref().unwrap()
    }

    fn restart(&mut self, change: Change) {
        let why = change.words();
        // The change belongs to the pause after this cycle (a caller's, if none was marked).
        if self.pauses.last().is_none_or(|p| p.after != self.cycle) {
            self.pauses.push(Pause { after: self.cycle, why: PauseWhy::Caller, lambda: self.lambda, changes: Vec::new(), restarted: None });
        }
        self.lambda = self.set.steering.lambda;
        self.rising = 0;
        self.stop = None;
        self.edited = true;
        let p = self.pauses.last_mut().unwrap();
        p.changes.push(change);
        p.restarted = Some(self.lambda);
        self.notes.push(format!("after cycle {}: {why}; λ restarted at {}", self.cycle, self.lambda));
    }

    /// Set a value during a pause (by its address); the damping restarts.
    pub fn set_value(&mut self, address: &str, value: f64) -> Result<(), String> {
        let f = crate::refinable::field_of(address, &self.model)?;
        self.model.set(f, value);
        self.set.enforce_links(&mut self.model);
        self.restart(Change::Value { address: address.to_string(), value });
        Ok(())
    }

    /// Change what is refined during a pause; the damping restarts.
    pub fn set_parameters(&mut self, set: ParameterSet) {
        self.set = set;
        self.set.enforce_links(&mut self.model);
        self.restart(Change::Parameters);
    }

    fn one_cycle(&mut self, on_cycle: &mut dyn FnMut(&Cycle) -> bool) {
        let clock = std::time::Instant::now();
        let np = self.set.free.len();
        let le_bail = self.model.has_extracted();
        let guards = self.set.steering.guards;
        if self.edited {
            self.edited = false;
            if let Some((what, value)) = impossible(&self.model, &self.model.calculate(&self.data.x)) {
                self.stop = Some(Stop::Diverging(Divergence::Impossible { what, value }));
                return;
            }
        }
        let cycle = self.cycle + 1;
        if le_bail {
            self.model.le_bail(&self.data);
            if np == 0 {
                self.cycle = cycle;
                let yc = self.model.calculate(&self.data.x).ycalc;
                let stats = statistics(&self.model, &self.data, &yc, 0);
                let c = Cycle { cycle, params: Vec::new(), ycalc: yc, chi2_start: stats.chi2, stats, marquardt: 0.0, tries: Vec::new(), singular: 0, restraints: Vec::new(), notes: Vec::new(), steering: Vec::new(), covariance: Vec::new(), spectrum: Vec::new(), seconds: clock.elapsed().as_secs_f64() };
                let go_on = on_cycle(&c);
                self.history.push(Cycle { ycalc: Vec::new(), ..c });
                if !go_on {
                    self.stop = Some(Stop::User);
                }
                return;
            }
        }
        let (model, data, set) = (&self.model, &self.data, &self.set);
        let fields = set.fields();
        let (yc, jq) = jacobian(model, data, set, &fields);
        let w = weights(model, data, &yc);
        let s0 = objective(data, &yc, &w);
        let mut a = vec![vec![0.0; np]; np];
        let mut g = vec![0.0; np];
        for i in 0..w.len() {
            if w[i] == 0.0 {
                continue;
            }
            let r = data.y[i] - yc[i];
            for j in 0..np {
                let wj = w[i] * jq[j][i];
                if wj == 0.0 {
                    continue;
                }
                g[j] += wj * r;
                for k in j..np {
                    a[j][k] += wj * jq[k][i];
                }
            }
        }
        symmetrize(&mut a);
        // Restraints (P40): rows weighted by χ²_ν of the data at the start of the cycle, as FullProf's (`restraints`).
        let used = w.iter().filter(|&&x| x > 0.0).count();
        let dof = used.saturating_sub(np).max(1);
        let c = s0 / dof as f64;
        let chi2_start = c;
        let rows = restraint_rows(model, set);
        add_restraint_rows(&mut a, &mut g, &rows, c);
        let s0 = s0 + c * penalty(&rows);
        // A singular matrix stops the run (R24), unless the guard is off: then the null directions are left out.
        let full = Solver::new(&a, (0..np).filter(|&j| a[j][j] > 0.0).collect());
        if guards.singular {
            if let Some((ratio, parts)) = full.null_directions().into_iter().next() {
                let dg = crate::check::null_direction(model, set, &parts, ratio);
                let mut combination: Vec<(String, f64)> = parts.iter().filter(|x| x.1.abs() > 0.2).map(|&(j, v)| (set.free[j].name.clone(), v)).collect();
                combination.sort_by(|x, y| y.1.abs().total_cmp(&x.1.abs()));
                self.stop = Some(Stop::Singular { combination, ratio, why: dg.why, tip: dg.suggestion.unwrap_or_default() });
                return;
            }
        }
        // The steering of each parameter: damping, and the largest step its values allow.
        let damping: Vec<f64> = set.free.iter().map(|p| p.damping).collect();
        let max_step: Vec<(f64, Option<String>)> = set
            .free
            .iter()
            .map(|p| {
                p.fields
                    .iter()
                    .filter_map(|&(f, m)| set.steering.max_shift_of(f, model).filter(|_| m != 0.0).map(|s| (s / m.abs(), crate::refinable::address(f, model))))
                    .fold((f64::INFINITY, None), |acc, (s, name)| if s < acc.0 { (s, Some(name)) } else { acc })
            })
            .collect();
        let bounds = shift_bounds(model, set);
        let count = movers(set);
        let mut mu = self.lambda;
        let mut held = vec![false; np];
        let mut tries = Vec::new();
        let mut notes = Vec::new();
        let mut steering = Vec::new();
        let mut accepted = None;
        let attempts = if set.steering.marquardt { 13 } else { 1 };
        for _ in 0..attempts {
            // The active set: parameters held at a limit while the step pushes them out.
            let mut active: Vec<bool> = (0..np).map(|j| a[j][j] > 0.0).collect();
            // The applied shift of each parameter held at a limit: onto the bound its step would cross.
            let mut fixed = vec![None::<f64>; np];
            let mut clipped: Vec<Option<String>>;
            let (mut computed, mut applied, mut solver);
            loop {
                // The others are solved with the held shifts taken out of the gradient.
                let mut ge = g.clone();
                for (h, fh) in fixed.iter().enumerate() {
                    if let Some(dh) = fh {
                        for (j, gj) in ge.iter_mut().enumerate() {
                            *gj -= a[j][h] * dh;
                        }
                    }
                }
                solver = Solver::new(&a, (0..np).filter(|&j| active[j]).collect());
                computed = solver.solve(&ge, mu, np);
                applied = (0..np).map(|j| computed[j] * damping[j]).collect::<Vec<f64>>();
                clipped = vec![None; np];
                for j in 0..np {
                    if applied[j].abs() > max_step[j].0 {
                        applied[j] = max_step[j].0.copysign(applied[j]);
                        clipped[j] = max_step[j].1.clone();
                    }
                }
                let mut changed = false;
                for j in 0..np {
                    if active[j] && (applied[j] < bounds[j].0 || applied[j] > bounds[j].1) {
                        active[j] = false;
                        held[j] = true;
                        fixed[j] = Some(if applied[j] < bounds[j].0 { bounds[j].0 } else { bounds[j].1 });
                        changed = true;
                    }
                }
                if !changed {
                    break;
                }
            }
            let mut step_notes = Vec::new();
            let mut step_data = Vec::new();
            for (j, fj) in fixed.iter().enumerate() {
                if let Some(dj) = fj {
                    applied[j] = *dj;
                    step_data.push(SteeringNote { kind: NoteKind::AtLimit, address: Some(set.free[j].name.clone()), value: Some(*dj), lambda: None, chi2: None });
                } else if let Some(f) = &clipped[j] {
                    step_notes.push(format!("{}: the step shortened to {f}'s maximum shift", set.free[j].name));
                    step_data.push(SteeringNote { kind: NoteKind::MaximumShift, address: Some(set.free[j].name.clone()), value: Some(applied[j].abs()), lambda: None, chi2: None });
                }
            }
            // Reflected and periodic bounds, on the fields only one parameter moves.
            for (j, p) in set.free.iter().enumerate() {
                for &(f, m) in &p.fields {
                    let b = set.behaviour_of(f);
                    if count[&f] != 1 || m == 0.0 || b == Behaviour::Stop {
                        continue;
                    }
                    let (lo, hi) = set.limits_of(f, model);
                    let v = model.get(f);
                    let mut t = v + m * applied[j];
                    let name = crate::refinable::address(f, model);
                    match b {
                        Behaviour::Reflect => {
                            let before = t;
                            for _ in 0..4 {
                                if let Some(h) = hi.filter(|&h| t > h) {
                                    t = 2.0 * h - t;
                                } else if let Some(l) = lo.filter(|&l| t < l) {
                                    t = 2.0 * l - t;
                                }
                            }
                            t = t.clamp(lo.unwrap_or(f64::NEG_INFINITY), hi.unwrap_or(f64::INFINITY));
                            if t != before {
                                step_notes.push(format!("{name}: reflected at its bound, to {t}"));
                                step_data.push(SteeringNote { kind: NoteKind::Reflected, address: Some(name.clone()), value: Some(t), lambda: None, chi2: None });
                            }
                        }
                        Behaviour::Periodic => {
                            if let (Some(l), Some(h)) = (lo, hi) {
                                let wrapped = l + (t - l).rem_euclid(h - l);
                                if wrapped != t {
                                    step_notes.push(format!("{name}: wrapped round its period, to {wrapped}"));
                                    step_data.push(SteeringNote { kind: NoteKind::Wrapped, address: Some(name.clone()), value: Some(wrapped), lambda: None, chi2: None });
                                    t = wrapped;
                                }
                            }
                        }
                        Behaviour::Stop => {}
                    }
                    applied[j] = (t - v) / m;
                }
            }
            let mut trial = model.clone();
            apply(&mut trial, set, &applied);
            let mut trial_set = set.clone();
            for (p, d) in trial_set.free.iter_mut().zip(&applied) {
                if let Some(v) = p.value.as_mut() {
                    *v += d;
                }
            }
            let yt = trial.calculate(&data.x).ycalc;
            let wt = weights(&trial, data, &yt);
            let s1 = objective(data, &yt, &wt) + c * penalty(&restraint_rows(&trial, &trial_set));
            if !set.steering.marquardt || (s1.is_finite() && s1 <= s0) {
                notes.extend(step_notes);
                steering.extend(step_data);
                accepted = Some((trial, trial_set, computed, applied, solver));
                break;
            }
            tries.push(Try { lambda: mu, chi2: s1 / dof as f64 });
            notes.push(format!("χ² rose to {:.4} with λ = {mu:e}: the step shortened", s1 / dof as f64));
            steering.push(SteeringNote { kind: NoteKind::Refused, address: None, value: None, lambda: Some(mu), chi2: Some(s1 / dof as f64) });
            mu = if mu == 0.0 { 1e-3 } else { mu * 10.0 };
            if mu > 1e8 {
                break;
            }
            held = vec![false; np];
        }
        let Some((trial, trial_set, computed, applied, solver)) = accepted else {
            if le_bail {
                // The intensities may still move: the next cycle extracts them again.
                self.cycle = cycle;
                return;
            }
            self.stop = Some(Stop::Minimum);
            return;
        };
        self.cycle = cycle;
        self.model = trial;
        self.set = trial_set;
        let (model, data, set) = (&self.model, &self.data, &self.set);
        let calc = model.calculate(&data.x);
        let stats = statistics(model, data, &calc.ycalc, np);
        let inv = solver.inverse(np);
        let esd: Vec<f64> = (0..np).map(|j| (inv[j][j] * stats.chi2).max(0.0).sqrt()).collect();
        let params: Vec<ParamState> = set
            .free
            .iter()
            .enumerate()
            .map(|(j, p)| {
                let (shift, e) = (p.shift_of(applied[j]), p.esd_of(esd[j]));
                ParamState { name: p.name.clone(), value: p.value_in(model), shift, esd: e, shift_esd: if e > 0.0 { shift / e } else { 0.0 }, computed: p.shift_of(computed[j]), at_limit: held[j] }
            })
            .collect();
        let restraints = restraint_states_at(model, set, stats.chi2, stats.n.saturating_sub(np));
        // Divergence (W11), after an impossible value: a computed shift of many esd's, χ² rising cycle after cycle.
        let mut diverging = impossible(model, &calc).map(|(what, value)| Divergence::Impossible { what, value });
        if diverging.is_none() && guards.jump > 0.0 && cycle >= guards.jump_from {
            let worst = (0..np).filter(|&j| !held[j] && esd[j] > 0.0).map(|j| (j, (computed[j] / esd[j]).abs())).max_by(|x, y| x.1.total_cmp(&y.1));
            if let Some((j, ratio)) = worst.filter(|w| w.1 > guards.jump) {
                diverging = Some(Divergence::Jump { parameter: set.free[j].name.clone(), ratio, cycle });
            }
        }
        if stats.chi2 > chi2_start * (1.0 + 1e-9) {
            if self.rising == 0 {
                self.rising_from = chi2_start;
            }
            self.rising += 1;
        } else {
            self.rising = 0;
        }
        if diverging.is_none() && guards.rising > 0 && self.rising >= guards.rising {
            diverging = Some(Divergence::Rising { cycles: self.rising, from: self.rising_from, to: stats.chi2 });
        }
        let done = (0..np).all(|j| held[j] || esd[j] == 0.0 || (computed[j] / esd[j]).abs() < self.opts.eps);
        let covariance: Vec<Vec<f64>> = inv.iter().map(|r| r.iter().map(|v| v * stats.chi2).collect()).collect();
        let top = full.lam.iter().copied().fold(0.0, f64::max);
        let mut spectrum: Vec<f64> = full.lam.iter().map(|l| if top > 0.0 { l / top } else { 0.0 }).collect();
        spectrum.sort_by(|x, y| y.total_cmp(x));
        let cyc = Cycle { cycle, params, ycalc: calc.ycalc, chi2_start, stats, marquardt: mu, tries, singular: solver.singular, restraints, notes, steering, covariance, spectrum, seconds: clock.elapsed().as_secs_f64() };
        let go_on = on_cycle(&cyc);
        self.history.push(Cycle { ycalc: Vec::new(), ..cyc });
        self.lambda = if mu < 1e-6 { 0.0 } else { mu / 10.0 };
        if done && !le_bail {
            self.stop = Some(Stop::Converged { eps: self.opts.eps });
        } else if let Some(d) = diverging {
            self.stop = Some(Stop::Diverging(d));
        } else if !go_on {
            self.stop = Some(Stop::User);
        }
    }

    /// The refinement at the values reached: the final covariance, esd's and correlations, every field moved with its
    /// esd, the restraints; `reason` the stop (a paused session: the cycles run).
    pub fn finish(&self) -> Refinement {
        // The pauses of the run; a step's end after the last cycle, with nothing changed, is the run's end, not a pause.
        let pauses: Vec<Pause> = self.pauses.iter().filter(|p| p.after < self.cycle || !p.changes.is_empty()).cloned().collect();
        let (model, data, set) = (&self.model, &self.data, &self.set);
        let fields = set.fields();
        let np = set.free.len();
        let (yc, jq) = jacobian(model, data, set, &fields);
        let w = weights(model, data, &yc);
        let mut a = vec![vec![0.0; np]; np];
        for i in 0..w.len() {
            if w[i] == 0.0 {
                continue;
            }
            for j in 0..np {
                for k in j..np {
                    a[j][k] += w[i] * jq[j][i] * jq[k][i];
                }
            }
        }
        symmetrize(&mut a);
        let stats = statistics(model, data, &yc, np);
        let berar = berar_factor(model, data, &yc, &model.calculate(&data.x).peaks);
        let rows = restraint_rows(model, set);
        add_restraint_rows(&mut a, &mut vec![0.0; np], &rows, stats.chi2);
        let restraints = restraint_states_at(model, set, stats.chi2, stats.n.saturating_sub(np));
        let solver = Solver::new(&a, (0..np).filter(|&j| a[j][j] > 0.0).collect());
        let cov = solver.inverse(np);
        let esd: Vec<f64> = (0..np).map(|j| (cov[j][j] * stats.chi2).max(0.0).sqrt()).collect();
        let correlations: Vec<Vec<f64>> = (0..np).map(|j| (0..np).map(|k| if cov[j][j] > 0.0 && cov[k][k] > 0.0 { cov[j][k] / (cov[j][j] * cov[k][k]).sqrt() } else { 0.0 }).collect()).collect();
        let bounds = shift_bounds(model, set);
        let params: Vec<ParamState> = set
            .free
            .iter()
            .enumerate()
            .map(|(j, p)| {
                let v = p.value_in(model);
                let at_limit = bounds[j].0.abs() < 1e-12 * (1.0 + v.abs()) || bounds[j].1.abs() < 1e-12 * (1.0 + v.abs());
                ParamState { name: p.name.clone(), value: v, shift: 0.0, esd: p.esd_of(esd[j]), shift_esd: 0.0, computed: 0.0, at_limit }
            })
            .collect();
        // Every field moved, with its esd √(mᵀ C m).
        let mut fields_out = Vec::new();
        for f in &fields {
            let m: Vec<f64> = set.free.iter().map(|p| p.fields.iter().filter(|q| q.0 == *f).map(|q| q.1).sum()).collect();
            let var: f64 = (0..np).map(|j| (0..np).map(|k| m[j] * cov[j][k] * m[k]).sum::<f64>()).sum::<f64>() * stats.chi2;
            fields_out.push((crate::refinable::address(*f, model), model.get(*f), var.max(0.0).sqrt()));
        }
        let covariance = cov.iter().map(|r| r.iter().map(|v| v * stats.chi2).collect()).collect();
        let reason = self.stop.clone().unwrap_or(Stop::Cycles { n: self.cycle, le_bail: model.has_extracted() });
        let shift_su_max = self.history.last().map_or(0.0, |c| c.params.iter().map(|q| q.shift_esd.abs()).fold(0.0, f64::max));
        Refinement { cycles: self.cycle, converged: reason.converged(), params, stats, correlations, covariance, fields: fields_out, restraints, stop: reason.text(), reason, shift_su_max, history: self.history.clone(), pauses, cycle_steps: Vec::new(), berar_factor: berar }
    }

    /// Each free parameter's path over the cycles run (the app's request 25): its name and, per cycle, the cycle, its
    /// value and its esd.
    pub fn paths(&self) -> Vec<(String, Path)> {
        self.set
            .free
            .iter()
            .map(|p| (p.name.clone(), self.history.iter().filter_map(|c| c.params.iter().find(|q| q.name == p.name).map(|q| (c.cycle, q.value, q.esd))).collect()))
            .collect()
    }

    /// The session's state as text, for a script to resume it (`resume`): the free parameters, every value they move
    /// and the variables' values (to the last digit), λ, the cycle count, the count of cycles χ² has risen, and the
    /// history (each cycle's statistics, λ, and each parameter's value, shift, esd and computed shift). A Le Bail
    /// phase's extracted intensities are not kept: they are extracted again at the next cycle.
    pub fn to_text(&self) -> String {
        use std::fmt::Write;
        let mut o = String::new();
        let _ = writeln!(o, "# drs refinement session (W11): resumed by `drs refine <control file> --session <this file>`");
        let _ = writeln!(o, "session 1");
        let _ = writeln!(o, "cycle {}", self.cycle);
        let _ = writeln!(o, "lambda {:?}", self.lambda);
        let _ = writeln!(o, "rising {} {:?}", self.rising, self.rising_from);
        if let Some(s) = &self.stop {
            let _ = writeln!(o, "stop {} {}", s.code(), s.text());
        }
        for p in &self.set.free {
            match p.value {
                Some(v) => {
                    let _ = writeln!(o, "variable {} {v:?}", p.name);
                }
                None => {
                    let _ = writeln!(o, "free {}", p.name);
                }
            }
        }
        for f in self.set.fields() {
            let _ = writeln!(o, "value {} {:?}", crate::refinable::address(f, &self.model), self.model.get(f));
        }
        for l in &self.set.links {
            let _ = writeln!(o, "value {} {:?}", crate::refinable::address(l.value, &self.model), self.model.get(l.value));
        }
        for c in &self.history {
            let s = &c.stats;
            let _ = writeln!(o, "history {} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {} {} {}", c.cycle, c.chi2_start, s.chi2, s.rp, s.rwp, s.rexp, s.durbin_watson, c.marquardt, s.n, s.p, c.singular);
            for p in &c.params {
                let _ = writeln!(o, "param {} {} {:?} {:?} {:?} {:?} {}", c.cycle, p.name, p.value, p.shift, p.esd, p.computed, p.at_limit);
            }
            for n in &c.notes {
                let _ = writeln!(o, "note {} {n}", c.cycle);
            }
        }
        // The pauses (E15a), each with what was changed before the run went on.
        for p in &self.pauses {
            let _ = writeln!(o, "pause {} {} {:?} {}", p.after, p.why.word(), p.lambda, p.restarted.map_or("-".to_string(), |v| format!("{v:?}")));
            for c in &p.changes {
                match c {
                    Change::Value { address, value } => {
                        let _ = writeln!(o, "change {} value {address} {value:?}", p.after);
                    }
                    Change::Parameters => {
                        let _ = writeln!(o, "change {} parameters", p.after);
                    }
                }
            }
        }
        o
    }

    /// A session resumed from `to_text`, on the same model, data and parameters it was saved with (the free
    /// parameters must be the same, in the same order). The run goes on from the saved values, λ and cycle count; a
    /// stopped session resumes paused, so that it can be run again.
    pub fn resume(text: &str, model: Model, data: Pattern, set: ParameterSet, opts: Options) -> Result<Session, String> {
        let mut s = Session::new(model, data, set, opts);
        s.notes.clear();
        let mut names = Vec::new();
        let num = |w: Option<&str>, what: &str| -> Result<f64, String> { w.and_then(|v| v.parse::<f64>().ok()).ok_or(format!("session: {what} is not a number")) };
        for (n, line) in text.lines().enumerate() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let mut w = line.split_whitespace();
            let key = w.next().unwrap_or("");
            let err = |e: String| format!("session line {}: {e}", n + 1);
            match key {
                "session" => {}
                "cycle" => s.cycle = num(w.next(), "cycle").map_err(err)? as usize,
                "lambda" => s.lambda = num(w.next(), "lambda").map_err(err)?,
                "rising" => {
                    s.rising = num(w.next(), "rising").map_err(err)? as usize;
                    s.rising_from = num(w.next(), "rising").map_err(err)?;
                }
                "stop" => {}
                "free" => names.push(w.next().unwrap_or("").to_string()),
                "variable" => {
                    let name = w.next().unwrap_or("").to_string();
                    let v = num(w.next(), "a variable's value").map_err(err)?;
                    if let Some(p) = s.set.free.iter_mut().find(|p| p.name == name && p.value.is_some()) {
                        p.value = Some(v);
                    }
                    names.push(name);
                }
                "value" => {
                    let addr = w.next().unwrap_or("");
                    let f = crate::refinable::field_of(addr, &s.model).map_err(err)?;
                    let v = num(w.next(), addr).map_err(err)?;
                    s.model.set(f, v);
                }
                "history" => {
                    let v: Vec<&str> = w.collect();
                    let x = |i: usize| v.get(i).and_then(|t| t.parse::<f64>().ok()).ok_or(err(format!("history field {}", i + 1)));
                    let stats = Statistics { n: x(8)? as usize, p: x(9)? as usize, c: 0, rp: x(3)?, rwp: x(4)?, rexp: x(5)?, chi2: x(2)?, durbin_watson: x(6)?, phases: Vec::new(), rp_background_free: 0.0, rwp_background_free: 0.0, rexp_background_free: 0.0 };
                    s.history.push(Cycle { cycle: x(0)? as usize, params: Vec::new(), ycalc: Vec::new(), stats, chi2_start: x(1)?, marquardt: x(7)?, tries: Vec::new(), singular: x(10)? as usize, restraints: Vec::new(), notes: Vec::new(), steering: Vec::new(), covariance: Vec::new(), spectrum: Vec::new(), seconds: 0.0 });
                }
                "param" => {
                    let v: Vec<&str> = w.collect();
                    let x = |i: usize| v.get(i).and_then(|t| t.parse::<f64>().ok()).ok_or(err(format!("param field {}", i + 1)));
                    let (value, shift, esd, computed) = (x(2)?, x(3)?, x(4)?, x(5)?);
                    let c = s.history.last_mut().ok_or(err("a param before its history line".into()))?;
                    c.params.push(ParamState { name: v.get(1).unwrap_or(&"").to_string(), value, shift, esd, shift_esd: if esd > 0.0 { shift / esd } else { 0.0 }, computed, at_limit: v.get(6) == Some(&"true") });
                }
                "note" => {
                    let rest = line.splitn(3, ' ').nth(2).unwrap_or("").to_string();
                    if let Some(c) = s.history.last_mut() {
                        c.notes.push(rest);
                    }
                }
                "pause" => {
                    let after = num(w.next(), "a pause's cycle").map_err(err)? as usize;
                    let why = match w.next() {
                        Some("step") => PauseWhy::Step,
                        Some("caller") => PauseWhy::Caller,
                        o => return Err(err(format!("a pause's why is step or caller, not {o:?}"))),
                    };
                    let lambda = num(w.next(), "a pause's λ").map_err(err)?;
                    let restarted = match w.next() {
                        Some("-") | None => None,
                        v => Some(num(v, "a pause's restarted λ").map_err(err)?),
                    };
                    s.pauses.push(Pause { after, why, lambda, changes: Vec::new(), restarted });
                }
                "change" => {
                    let _after = w.next();
                    let c = match w.next() {
                        Some("value") => {
                            let address = w.next().unwrap_or("").to_string();
                            Change::Value { address, value: num(w.next(), "a change's value").map_err(err)? }
                        }
                        Some("parameters") => Change::Parameters,
                        o => return Err(err(format!("a change is value or parameters, not {o:?}"))),
                    };
                    s.pauses.last_mut().ok_or(err("a change before its pause line".into()))?.changes.push(c);
                }
                other => return Err(err(format!("unknown line '{other}'"))),
            }
        }
        let ours: Vec<&String> = s.set.free.iter().map(|p| &p.name).collect();
        if names.iter().collect::<Vec<_>>() != ours {
            return Err(format!("the session was saved with the free parameters {names:?}, not these ({ours:?})"));
        }
        s.edited = true;
        Ok(s)
    }
}

/// Refine the free parameters of `set` in `model` against `data`; `on_cycle` receives every cycle.
pub fn refine(model: &mut Model, data: &Pattern, set: &ParameterSet, opts: Options, mut on_cycle: impl FnMut(&Cycle)) -> Refinement {
    refine_until(model, data, set, opts, |c| {
        on_cycle(c);
        true
    })
}

/// As `refine`, with `on_cycle` saying after every cycle whether to go on (the app's request 6, L070): when it returns
/// false the refinement stops there, with the model at that cycle's values, their esd's and correlations, and
/// `stop` = "stopped by the user". A caller with a Stop button keeps an `AtomicBool` and returns its negation. A
/// `Session` run to its end.
pub fn refine_until(model: &mut Model, data: &Pattern, set: &ParameterSet, opts: Options, on_cycle: impl FnMut(&Cycle) -> bool) -> Refinement {
    let mut s = Session::new(model.clone(), data.clone(), set.clone(), opts);
    s.run(on_cycle);
    *model = s.model.clone();
    s.finish()
}

/// A step of a plan run (`run_plan`) as it starts (`stop` None) and as it ends.
#[derive(Clone, Copy, Debug)]
pub struct PlanStepRun<'a> {
    /// The step's index in the plan (0, 1, …) and its name.
    pub index: usize,
    pub name: &'a str,
    /// The parameters it frees and the cycles it may run.
    pub parameters: usize,
    pub max_cycles: usize,
    pub stop: Option<&'a Stop>,
}

/// A plan run (L110's `[[refine.step]]`, W12; P37): the steps `steps` of `rf`'s plan in that order, each a `Session` of
/// the parameters it frees (`Refinement::for_step`), started from the model the step before reached, the variables
/// carrying their values to the next step. The run stops after a step stopped by the user (`on_cycle` returning
/// false), a divergence or a singular matrix. `on_step` receives each step as it starts and as it ends; `on_cycle` each
/// cycle as the step's session numbers it (1, 2, … within the step).
///
/// The refinement returned is the last step's (its free parameters, esd's, correlations, statistics and stop) with the
/// cycles of every step: `history` in the order run, numbered on from one step to the next (1, 2, … over the run),
/// `cycle_steps` the step of each, `cycles` their count, and the pauses after the cycles so numbered. The last step's
/// session is returned with it (its parameter set and trajectory). `model` is left at the values reached.
pub fn run_plan(model: &mut Model, data: &Pattern, rf: &crate::refinable::Refinement, steps: &[usize], mut on_step: impl FnMut(&PlanStepRun), mut on_cycle: impl FnMut(&Cycle) -> bool) -> Result<(Refinement, Session), String> {
    let mut rf = rf.clone();
    let (mut history, mut cycle_steps, mut pauses) = (Vec::new(), Vec::new(), Vec::new());
    let mut last: Option<(Refinement, Session)> = None;
    for &k in steps {
        let rk = rf.for_step(k, model)?;
        let set = rk.build(model)?;
        let opts = Options::stated(rk.cycles, rk.eps);
        let name = rf.plan[k].name.clone();
        on_step(&PlanStepRun { index: k, name: &name, parameters: set.free.len(), max_cycles: opts.max_cycles, stop: None });
        let mut session = Session::new(model.clone(), data.clone(), set, opts);
        let stop = session.run(&mut on_cycle).clone();
        on_step(&PlanStepRun { index: k, name: &name, parameters: session.set.free.len(), max_cycles: opts.max_cycles, stop: Some(&stop) });
        *model = session.model.clone();
        for v in rf.variables.iter_mut() {
            if let Some(x) = session.set.free.iter().find(|p| p.name == v.name).and_then(|p| p.value) {
                v.value = x;
            }
        }
        let r = session.finish();
        let offset = history.len();
        history.extend(r.history.iter().map(|c| Cycle { cycle: c.cycle + offset, ..c.clone() }));
        cycle_steps.extend(r.history.iter().map(|_| (k + 1, name.clone())));
        pauses.extend(r.pauses.iter().map(|p| Pause { after: p.after + offset, ..p.clone() }));
        last = Some((r, session));
        if matches!(stop, Stop::User | Stop::Diverging(_) | Stop::Singular { .. }) {
            break;
        }
    }
    let (mut r, session) = last.ok_or("the plan has no step to run")?;
    r.cycles = history.len();
    r.history = history;
    r.cycle_steps = cycle_steps;
    r.pauses = pauses;
    Ok((r, session))
}

/// A state at given values of free parameters (Explore's E20a): the model, its pattern, Yobs − Ycal and the statistics.
#[derive(Clone, Debug)]
pub struct PointState {
    pub model: Model,
    pub ycalc: Vec<f64>,
    pub difference: Vec<f64>,
    pub statistics: Statistics,
}

/// The state with the free parameters named at the values given (E20a, the dragged point of Explore's §2.7), each
/// set as a step of the engine sets it: its value (`Free::value_in`: a variable's own, else its first field's) moved
/// to the one given, every field it moves by its coefficient times the parameter's shift, then the links enforced
/// (`ParameterSet::enforce_links`); the others at the model's values. χ² over N − P, P the set's free parameters (the
/// restraints' penalty not added), as `chi2_grid`'s nodes, which are computed by it.
pub fn point(model: &Model, data: &Pattern, set: &ParameterSet, values: &[(&str, f64)]) -> Result<PointState, String> {
    let mut delta = vec![0.0; set.free.len()];
    for &(name, v) in values {
        let j = set.free.iter().position(|p| p.name == name).ok_or_else(|| format!("{name} is not a free parameter of the refinement"))?;
        let p = &set.free[j];
        let dv = v - p.value_in(model);
        // The parameter's shift for a change of its value (`Free::shift_of` inverted).
        delta[j] = if p.value.is_some() { dv } else { dv / p.fields[0].1 };
    }
    let mut m = model.clone();
    apply(&mut m, set, &delta);
    set.enforce_links(&mut m);
    let ycalc = m.calculate(&data.x).ycalc;
    let statistics = statistics(&m, data, &ycalc, set.free.len());
    let difference = crate::diagnostics::difference(&data.y, &ycalc);
    Ok(PointState { model: m, ycalc, difference, statistics })
}

/// A contour line of a χ² grid at one level (E21): the points of its pieces in the grid's units (a along x, b along
/// y), NaN between pieces, as the shared plot draws a map's contours.
#[derive(Clone, Debug, PartialEq)]
pub struct Contour {
    pub level: f64,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
}

/// χ² over a grid of two free parameters' values (E2).
#[derive(Clone, Debug, PartialEq)]
pub struct Grid {
    /// The two parameters' names.
    pub a: String,
    pub b: String,
    /// Their values along the grid (the first field's value, as `ParamState::value`).
    pub xa: Vec<f64>,
    pub xb: Vec<f64>,
    /// χ² at each point, `chi2[i][j]` at (xa[i], xb[j]); NaN where not computed (the run cancelled).
    pub chi2: Vec<Vec<f64>>,
    /// The points in the order computed (`coarse_to_fine`), and how many were.
    pub order: Vec<(usize, usize)>,
    pub computed: usize,
    /// N − P + C of χ² (`Statistics`), which scales the esd's (0 when not known): the ellipses' levels in χ² (E21).
    pub dof: usize,
}

impl Grid {
    /// The lowest χ² computed: (i, j, χ²).
    pub fn minimum(&self) -> Option<(usize, usize, f64)> {
        let mut best: Option<(usize, usize, f64)> = None;
        for (i, row) in self.chi2.iter().enumerate() {
            for (j, &v) in row.iter().enumerate() {
                if v.is_finite() && best.is_none_or(|b| v < b.2) {
                    best = Some((i, j, v));
                }
            }
        }
        best
    }

    /// The contours of the conditional cut at the levels of the pair's kσ ellipses (E21; `covariance_ellipse`, the
    /// marginal's): the esd's are scaled by χ² (S/(N − P + C)), so the kσ contour of S is a rise of k²χ² in S, of
    /// k²χ²/(N − P + C) in χ², above the grid's minimum: χ²_min(1 + k²/dof). k = 1 and 2 for Δχ² = 1 and 4 in units of
    /// the esd's. None of them without the grid's `dof`.
    pub fn contours_at_sigma(&self, ks: &[f64]) -> Vec<Contour> {
        match (self.minimum(), self.dof) {
            (Some((_, _, m)), d) if d > 0 => self.contours(&ks.iter().map(|k| m * (1.0 + k * k / d as f64)).collect::<Vec<_>>()),
            _ => Vec::new(),
        }
    }

    /// The contour lines at each level (E21), by marching squares over the cells whose four nodes are computed, the
    /// crossing on each edge by linear interpolation, a cell crossed four times (a saddle) split by its centre's mean;
    /// the segments joined into pieces (a closed piece ends at its start), a cell with a node not computed leaving a
    /// break.
    pub fn contours(&self, levels: &[f64]) -> Vec<Contour> {
        levels.iter().map(|&l| self.contour(l)).collect()
    }

    fn contour(&self, level: f64) -> Contour {
        let (na, nb) = (self.xa.len(), self.xb.len());
        // An edge by its nodes: (i, j, 0) from (i, j) to (i + 1, j), (i, j, 1) from (i, j) to (i, j + 1).
        type Edge = (usize, usize, u8);
        let z = |i: usize, j: usize| self.chi2[i][j];
        let at = |e: Edge| -> (f64, f64) {
            let (i, j, d) = e;
            let (i2, j2) = if d == 0 { (i + 1, j) } else { (i, j + 1) };
            let t = (level - z(i, j)) / (z(i2, j2) - z(i, j));
            (self.xa[i] + t * (self.xa[i2] - self.xa[i]), self.xb[j] + t * (self.xb[j2] - self.xb[j]))
        };
        let mut segments: Vec<(Edge, Edge)> = Vec::new();
        for i in 0..na.saturating_sub(1) {
            for j in 0..nb.saturating_sub(1) {
                let v = [z(i, j), z(i + 1, j), z(i + 1, j + 1), z(i, j + 1)];
                if v.iter().any(|x| !x.is_finite()) {
                    continue;
                }
                // The edges bottom, right, top, left, each between corners k and k + 1.
                let edges: [Edge; 4] = [(i, j, 0), (i + 1, j, 1), (i, j + 1, 0), (i, j, 1)];
                let above = v.map(|x| x > level);
                let crossed: Vec<usize> = (0..4).filter(|&k| above[k] != above[(k + 1) % 4]).collect();
                match crossed.len() {
                    2 => segments.push((edges[crossed[0]], edges[crossed[1]])),
                    4 => {
                        let centre = v.iter().sum::<f64>() / 4.0 > level;
                        if centre == above[0] {
                            // Corners 1 and 3 cut off.
                            segments.push((edges[0], edges[1]));
                            segments.push((edges[2], edges[3]));
                        } else {
                            segments.push((edges[3], edges[0]));
                            segments.push((edges[1], edges[2]));
                        }
                    }
                    _ => {}
                }
            }
        }
        // Join the segments at their shared edges.
        let mut by_edge: std::collections::HashMap<Edge, Vec<usize>> = std::collections::HashMap::new();
        for (k, &(a, b)) in segments.iter().enumerate() {
            by_edge.entry(a).or_default().push(k);
            by_edge.entry(b).or_default().push(k);
        }
        let mut used = vec![false; segments.len()];
        let (mut x, mut y) = (Vec::new(), Vec::new());
        let next = |e: Edge, used: &[bool]| by_edge.get(&e).and_then(|v| v.iter().copied().find(|&k| !used[k]));
        for start in 0..segments.len() {
            if used[start] {
                continue;
            }
            used[start] = true;
            // Walk backwards first, so that an open piece is followed from one end.
            let (mut head, mut tail) = segments[start];
            let mut piece = std::collections::VecDeque::from([head, tail]);
            while let Some(k) = next(tail, &used) {
                used[k] = true;
                let (a, b) = segments[k];
                tail = if a == tail { b } else { a };
                piece.push_back(tail);
            }
            while let Some(k) = next(head, &used) {
                used[k] = true;
                let (a, b) = segments[k];
                head = if a == head { b } else { a };
                piece.push_front(head);
            }
            if !x.is_empty() {
                x.push(f64::NAN);
                y.push(f64::NAN);
            }
            for e in piece {
                let (px, py) = at(e);
                x.push(px);
                y.push(py);
            }
        }
        Contour { level, x, y }
    }

    /// The grid as a map (E13a): columns at `xa`, rows at `xb` with their edges as `plots::map` makes them (midpoints,
    /// half a step out at the ends), χ² on the scale shown (None where not computed, or ≤ 0 on log), the colour range
    /// over the points computed. The same numbers as a plot set of one member per value of b, without the set.
    pub fn map(&self, scale: crate::plots::YScale) -> Result<crate::plots::Map, String> {
        if self.xa.is_empty() || self.xb.is_empty() {
            return Err("the grid is empty".into());
        }
        let rows = (0..self.xb.len()).map(|j| (j, self.xb[j], self.chi2.iter().map(|row| Some(row[j])).collect())).collect();
        Ok(crate::plots::map_of_rows(&self.xa, rows, Vec::new(), None, scale, None))
    }
}

impl Ellipse {
    /// n + 1 points of the ellipse at k σ about `centre` (E13b), the last the first again: at t = 2π i/n, the point
    /// k(√λ₁ cos t, √λ₂ sin t) turned by `angle`, as the pair's view drew it (`draw.js`, `ellipseOutline`).
    pub fn outline(&self, centre: (f64, f64), k: f64, n: usize) -> Vec<(f64, f64)> {
        let n = n.max(1);
        let (c, s) = (self.angle.cos(), self.angle.sin());
        let (ra, rb) = (k * self.l1.sqrt(), k * self.l2.sqrt());
        let mut out: Vec<(f64, f64)> = (0..n)
            .map(|i| {
                let t = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
                let (u, v) = (ra * t.cos(), rb * t.sin());
                (centre.0 + u * c - v * s, centre.1 + u * s + v * c)
            })
            .collect();
        out.push(out[0]);
        out
    }
}

/// The points of an n × n grid from coarse to fine (E2): the centre, the four corners, then the points of the grid of
/// stride s for s = the power of two at or above n − 1, halved down to 1 (the last row and column with each), each
/// point once. The first (n−1)/s + 1 squared make a picture the later ones refine.
pub fn coarse_to_fine(n: usize) -> Vec<(usize, usize)> {
    let mut seen = vec![vec![false; n]; n];
    let mut out = Vec::with_capacity(n * n);
    let mut push = |i: usize, j: usize, out: &mut Vec<(usize, usize)>| {
        if i < n && j < n && !seen[i][j] {
            seen[i][j] = true;
            out.push((i, j));
        }
    };
    if n == 0 {
        return out;
    }
    let c = (n - 1) / 2;
    push(c, c, &mut out);
    for (i, j) in [(0, 0), (n - 1, 0), (0, n - 1), (n - 1, n - 1)] {
        push(i, j, &mut out);
    }
    let mut stride = (n - 1).max(1).next_power_of_two();
    loop {
        let along = |s: usize| -> Vec<usize> {
            let mut v: Vec<usize> = (0..n).step_by(s).collect();
            if *v.last().unwrap() != n - 1 {
                v.push(n - 1);
            }
            v
        };
        for &i in &along(stride) {
            for &j in &along(stride) {
                push(i, j, &mut out);
            }
        }
        if stride == 1 {
            break;
        }
        stride /= 2;
    }
    out
}

/// χ² on an n × n grid of two free parameters (E2), indices `a` and `b` into `set.free`, each over its value ± its
/// span (in the value's units, `ParamState::value`), the others at the model's values: at each point the state of
/// `point` (the two parameters' fields moved by their multipliers as a step moves them, the links enforced; χ² over
/// N − P, P the set's free parameters; the restraints' penalty not added). The points in `coarse_to_fine`'s
/// order, on all the machine's threads in batches; after each batch `on_progress(done, total)`, false to stop (the
/// points not computed are NaN).
pub fn chi2_grid(model: &Model, data: &Pattern, set: &ParameterSet, (a, b): (usize, usize), span: (f64, f64), n: usize, mut on_progress: impl FnMut(usize, usize) -> bool) -> Grid {
    let n = n.max(2);
    let (pa, pb) = (&set.free[a], &set.free[b]);
    let (va, vb) = (pa.value_in(model), pb.value_in(model));
    let along = |v: f64, s: f64| -> Vec<f64> { (0..n).map(|k| v - s + 2.0 * s * k as f64 / (n - 1) as f64).collect() };
    let (xa, xb) = (along(va, span.0), along(vb, span.1));
    // Each node as `point` sets it (E20a): the two parameters as a step moves them, the links enforced.
    let point = |i: usize, j: usize| -> (f64, usize) {
        point(model, data, set, &[(pa.name.as_str(), xa[i]), (pb.name.as_str(), xb[j])]).map_or((f64::NAN, 0), |p| (p.statistics.chi2, (p.statistics.n + p.statistics.c).saturating_sub(p.statistics.p)))
    };
    let order = coarse_to_fine(n);
    let total = order.len();
    let mut chi2 = vec![vec![f64::NAN; n]; n];
    let threads = worker_threads();
    let (mut done, mut dof) = (0, 0);
    for batch in order.chunks(threads) {
        let values: Vec<(f64, usize)> = std::thread::scope(|sc| {
            let handles: Vec<_> = batch.iter().map(|&(i, j)| sc.spawn(move || point(i, j))).collect();
            handles.into_iter().map(|h| h.join().unwrap_or((f64::NAN, 0))).collect()
        });
        for (&(i, j), (v, d)) in batch.iter().zip(values) {
            chi2[i][j] = v;
            dof = dof.max(d);
        }
        done += batch.len();
        if !on_progress(done, total) {
            break;
        }
    }
    Grid { a: pa.name.clone(), b: pb.name.clone(), xa, xb, chi2, order, computed: done, dof }
}

#[cfg(test)]
mod tests {
    use super::{covariance_ellipse, householder_ql, jacobi, symmetric_eigen};

    /// The pair's ellipse (L138): uncorrelated, the axes are the variances along a and b; equal esd's with r = 0.9, the
    /// axes are σ²(1 ± r) at 45°; r = −0.9 at −45° (135°); λ₁λ₂ = det C always; no ellipse without both esd's.
    /// L145: the options stated or defaulted; the strong correlations, strongest first; the map's half-widths.
    /// Householder + QL against Jacobi's rotations on a symmetric matrix of 150 (pseudo-random, positive definite,
    /// with a near-null direction as a Pawley fit's accidental overlaps leave): the same eigenvalues to 10⁻¹¹ of the
    /// largest, each pair A v = λ v to 10⁻¹¹, the vectors orthonormal to 10⁻¹².
    #[test]
    fn householder_ql_as_jacobi() {
        let n = 150;
        let mut seed = 12345u64;
        let mut rnd = || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((seed >> 11) as f64 / (1u64 << 53) as f64) - 0.5
        };
        let b: Vec<Vec<f64>> = (0..n).map(|_| (0..n).map(|_| rnd()).collect()).collect();
        let mut a: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|j| (0..n).map(|k| b[i][k] * b[j][k]).sum()).collect()).collect();
        // Two columns nearly equal: a near-null direction.
        for row in a.iter_mut() {
            row[1] = row[0] * (1.0 + 1e-9);
        }
        let first = a[0].clone();
        for (x, y) in a[1].iter_mut().zip(&first) {
            *x = y * (1.0 + 1e-9);
        }
        a[1][1] = a[0][0] * (1.0 + 1e-9) * (1.0 + 1e-9);
        let (l1, v1) = householder_ql(&a);
        let (mut l2, _) = jacobi(&a);
        l2.sort_by(f64::total_cmp);
        let top = l1.last().unwrap().abs();
        for k in 0..n {
            assert!((l1[k] - l2[k]).abs() < 1e-11 * top, "eigenvalue {k}: {} against {}", l1[k], l2[k]);
            for i in 0..n {
                let av: f64 = (0..n).map(|j| a[i][j] * v1[j][k]).sum();
                assert!((av - l1[k] * v1[i][k]).abs() < 1e-11 * top, "A v = λ v, {k}, row {i}");
            }
            for q in 0..n {
                let dot: f64 = (0..n).map(|i| v1[i][k] * v1[i][q]).sum();
                assert!((dot - if k == q { 1.0 } else { 0.0 }).abs() < 1e-12, "orthonormal {k} {q}: {dot}");
            }
        }
        assert!(symmetric_eigen(&a).0 == l1);
    }

    #[test]
    fn options_correlations_and_half_widths() {
        let o = super::Options::stated(None, None);
        assert_eq!((o.max_cycles, o.eps), (super::DEFAULT_CYCLES, super::DEFAULT_EPS));
        let o = super::Options::stated(Some(0), Some(0.0));
        assert_eq!((o.max_cycles, o.eps), (1, 0.0));
        let r = vec![vec![1.0, 0.85, -0.95], vec![0.85, 1.0, 0.1], vec![-0.95, 0.1, 1.0]];
        assert_eq!(super::strong_correlations(&r), vec![(0, 2, -0.95), (0, 1, 0.85)]);
        assert_eq!(super::half_widths(4.0, 0.5, 2.0), (2.0, 8.0));
        assert_eq!(super::chi2_range([4.0, f64::NAN, 1.0, 9.0]), Some((1.0, 3.0, 9.0)));
        assert_eq!(super::strongest_pair(&r, |_| true), Some((0, 2)));
        assert_eq!(super::strongest_pair(&r, |j| j != 2), Some((0, 1)));
    }

    /// L146: a restraint's deviation in σ.
    #[test]
    fn a_restraints_deviation() {
        let r = super::RestraintState { label: "d".into(), target: 2.0, value: 2.3, sigma: 0.1, weight: 100.0, contribution: 9.0 };
        assert!((r.deviation().unwrap() - 3.0).abs() < 1e-12);
        assert_eq!(super::RestraintState { sigma: 0.0, ..r }.deviation(), None);
    }

    /// L150: half the logical processors, at least one.
    #[test]
    fn the_cores_worker_threads() {
        let all = std::thread::available_parallelism().map_or(1, |t| t.get());
        assert_eq!(super::worker_threads(), all.div_ceil(2).max(1));
        assert!(super::worker_threads() >= 1 && super::worker_threads() <= all);
    }

    #[test]
    fn the_share_of_a_run_done() {
        assert_eq!(super::fraction_done(3, 12), 0.25);
        assert_eq!(super::fraction_done(15, 10), 1.0, "a resumed run beyond its count");
        assert_eq!(super::fraction_done(0, 0), 1.0);
    }

    #[test]
    fn the_covariance_ellipse_of_a_pair() {
        let e = covariance_ellipse(2.0, 1.0, 0.0).unwrap();
        assert_eq!((e.l1, e.l2, e.angle), (4.0, 1.0, 0.0));
        assert_eq!(covariance_ellipse(1.0, 2.0, 0.0).unwrap().angle, std::f64::consts::FRAC_PI_2);
        let e = covariance_ellipse(0.5, 0.5, 0.9).unwrap();
        assert!((e.l1 - 0.25 * 1.9).abs() < 1e-15 && (e.l2 - 0.25 * 0.1).abs() < 1e-15);
        assert!((e.angle - std::f64::consts::FRAC_PI_4).abs() < 1e-12);
        // an axis is a direction modulo π: −45° is 135°
        assert!((covariance_ellipse(0.5, 0.5, -0.9).unwrap().angle.rem_euclid(std::f64::consts::PI) - 3.0 * std::f64::consts::FRAC_PI_4).abs() < 1e-12);
        let (sa, sb, r) = (3e-4, 7e-3, 0.857);
        let e = covariance_ellipse(sa, sb, r).unwrap();
        assert!((e.l1 * e.l2 / (sa * sa * sb * sb * (1.0 - r * r)) - 1.0).abs() < 1e-9);
        assert!(covariance_ellipse(0.0, 1.0, 0.5).is_none());
    }

    #[test]
    fn eigen_decomposition_rebuilds_the_matrix() {
        let a = vec![vec![4.0, 1.0, -2.0], vec![1.0, 2.0, 0.5], vec![-2.0, 0.5, 3.0]];
        let (l, v) = symmetric_eigen(&a);
        for i in 0..3 {
            for j in 0..3 {
                let r: f64 = (0..3).map(|k| v[i][k] * l[k] * v[j][k]).sum();
                assert!((r - a[i][j]).abs() < 1e-12, "({i}, {j}): {r} against {}", a[i][j]);
                let o: f64 = (0..3).map(|k| v[k][i] * v[k][j]).sum();
                assert!((o - if i == j { 1.0 } else { 0.0 }).abs() < 1e-12);
            }
        }
        // The trace is the sum of the eigenvalues.
        assert!((l.iter().sum::<f64>() - 9.0).abs() < 1e-12);
    }
}

// E11: two refinements' paths compared ------------------------------------------------------------------------------

/// Two ends agree when they differ by at most this many combined σ, √(σ_A² + σ_B²): the usual test of two independent
/// estimates of one value (E11).
pub const ENDS_AGREE: f64 = 3.0;

/// Which run of a comparison is the reference (its end is the common origin, its end esd's the common unit).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Side {
    #[default]
    A,
    B,
}

/// A run to compare: its path, its esd's at the end and the correlations there (both by its own `set.free` order: a
/// finished run's final covariance, `Refinement::params` and `correlations`; without them the last cycle's), and why
/// it stopped (`Stop::code`, `Stop::text`) when known.
#[derive(Clone, Copy, Debug)]
pub struct Run<'a> {
    pub path: &'a Trajectory,
    pub end_esd: Option<&'a [f64]>,
    pub end_correlations: Option<&'a [Vec<f64>]>,
    pub stop: Option<(&'a str, &'a str)>,
}

/// One common parameter's ends compared (E11): B's end less A's in its unit, in the reference's σ, and in the combined
/// σ; whether they agree (`ENDS_AGREE`).
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct EndDifference {
    pub address: String,
    pub delta: f64,
    pub in_scale: Option<f64>,
    pub combined: Option<f64>,
    pub agree: bool,
}

/// The comparison in numbers (E11).
#[derive(Clone, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct CompareSummary {
    pub common: usize,
    pub only_a: usize,
    pub only_b: usize,
    pub fixed_differences: usize,
    pub agree: usize,
    pub disagree: usize,
    /// The largest |combined| and its address.
    pub largest: Option<(String, f64)>,
    pub last_chi2_a: Option<f64>,
    pub last_chi2_b: Option<f64>,
    /// Each run's stop: (code, text), when known.
    pub stop_a: Option<(String, String)>,
    pub stop_b: Option<(String, String)>,
}

/// Two refinements' paths compared (Explore's E11, R31): every number the compare view shows, by parameter address.
/// `[state][j]` and `[step][j]` run over `common`.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PathsCompared {
    pub reference: Side,
    /// Free in both runs, in A's order.
    pub common: Vec<String>,
    /// Free in one run only: (address, start, end) in that run.
    pub only_a: Vec<(String, f64, f64)>,
    pub only_b: Vec<(String, f64, f64)>,
    /// Free in neither run and different at the two starts: (address, A's, B's).
    pub fixed_differences: Vec<(String, f64, f64)>,
    /// Each common parameter's unit, from the dictionary ("" for a variable).
    pub units: Vec<String>,
    pub states_a: Vec<usize>,
    pub states_b: Vec<usize>,
    pub chi2_a: Vec<f64>,
    pub chi2_b: Vec<f64>,
    pub values_a: Vec<Vec<f64>>,
    pub values_b: Vec<Vec<f64>>,
    pub esd_a: Vec<Vec<Option<f64>>>,
    pub esd_b: Vec<Vec<Option<f64>>>,
    pub end_esd_a: Vec<f64>,
    pub end_esd_b: Vec<f64>,
    /// The reference's end value and end esd of each common parameter: the common origin and unit of both paths.
    pub origin: Vec<f64>,
    pub scale: Vec<f64>,
    pub a_scaled: Vec<Vec<Option<f64>>>,
    pub b_scaled: Vec<Vec<Option<f64>>>,
    /// The comparison's time axis: per step, A's state and B's state shown; the shorter run held at its end.
    pub steps: Vec<(usize, usize)>,
    pub held_a: Vec<bool>,
    pub held_b: Vec<bool>,
    /// (B − A)/scale at each step's two states.
    pub difference: Vec<Vec<Option<f64>>>,
    pub end: Vec<EndDifference>,
    pub summary: CompareSummary,
    /// Each run's correlations at its end over `common` (for `pair`), when known.
    pub end_correlations_a: Option<Vec<Vec<f64>>>,
    pub end_correlations_b: Option<Vec<Vec<f64>>>,
}

/// A pair of common parameters in both runs (E11a): each run's path in the parameters' units, its end, its σ of the
/// two at its end (the view's extent along each axis) and its covariance ellipse there.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PairCompared {
    pub a: String,
    pub b: String,
    pub path_a: Vec<(f64, f64)>,
    pub path_b: Vec<(f64, f64)>,
    pub centre_a: (f64, f64),
    pub centre_b: (f64, f64),
    pub sigmas_a: (f64, f64),
    pub sigmas_b: (f64, f64),
    pub correlation_a: Option<f64>,
    pub correlation_b: Option<f64>,
    pub ellipse_a: Option<Ellipse>,
    pub ellipse_b: Option<Ellipse>,
}

/// One run's numbers by its own free parameters: states, χ², values and esd's per state, end esd's, end correlations.
struct RunTable {
    names: Vec<String>,
    states: Vec<usize>,
    chi2: Vec<f64>,
    values: Vec<Vec<f64>>,
    esd: Vec<Vec<Option<f64>>>,
    end_esd: Vec<f64>,
    correlations: Option<Vec<Vec<f64>>>,
}

fn run_table(r: &Run) -> RunTable {
    let t = r.path;
    let names: Vec<String> = t.set.free.iter().map(|p| p.name.clone()).collect();
    let mut values = vec![t.set.free.iter().map(|p| p.value_in(&t.start)).collect::<Vec<f64>>()];
    let mut esd = vec![vec![None; names.len()]];
    let mut states = vec![0];
    let mut chi2 = vec![t.cycles.first().map_or(f64::NAN, |c| c.chi2_start)];
    for c in &t.cycles {
        let prev = values.last().unwrap().clone();
        values.push(names.iter().enumerate().map(|(j, n)| c.params.iter().find(|q| &q.name == n).map_or(prev[j], |q| q.value)).collect());
        esd.push(names.iter().map(|n| c.params.iter().find(|q| &q.name == n).map(|q| q.esd)).collect());
        states.push(c.cycle);
        chi2.push(c.stats.chi2);
    }
    let end_esd = match r.end_esd {
        Some(e) => e.to_vec(),
        None => esd.last().unwrap().iter().map(|e| e.unwrap_or(0.0)).collect(),
    };
    let correlations = match r.end_correlations {
        Some(c) => Some(c.to_vec()),
        None => t.cycles.last().filter(|c| c.covariance.len() == names.len() && !names.is_empty()).map(|c| correlations_of(&c.covariance)),
    };
    RunTable { names, states, chi2, values, esd, end_esd, correlations }
}

/// Every field of a model by address, with its value.
fn fields_by_address(m: &Model) -> Vec<(String, f64)> {
    m.all_fields().into_iter().map(|f| (crate::refinable::address(f, m), m.get(f))).collect()
}

/// Two refinements' paths compared (E11; see `PathsCompared`): the parameters in common and apart, the fixed values
/// changed between the two starts, both paths on the reference's scale (its end the origin, its end esd's the unit),
/// the difference at each step of the longer run (the shorter held at its end), the ends compared, a summary.
pub fn compare(a: &Run, b: &Run, reference: Side) -> PathsCompared {
    let (ta, tb) = (run_table(a), run_table(b));
    let ja: Vec<usize> = (0..ta.names.len()).filter(|&j| tb.names.contains(&ta.names[j])).collect();
    let jb: Vec<usize> = ja.iter().map(|&j| tb.names.iter().position(|n| n == &ta.names[j]).unwrap()).collect();
    let common: Vec<String> = ja.iter().map(|&j| ta.names[j].clone()).collect();
    let apart = |t: &RunTable, other: &RunTable| -> Vec<(String, f64, f64)> {
        (0..t.names.len()).filter(|&j| !other.names.contains(&t.names[j])).map(|j| (t.names[j].clone(), t.values[0][j], t.values[t.values.len() - 1][j])).collect()
    };
    let (only_a, only_b) = (apart(&ta, &tb), apart(&tb, &ta));
    // The values free in neither run, compared between the two starts.
    let moved = |t: &Trajectory| -> Vec<String> { t.set.free.iter().flat_map(|p| p.fields.iter().map(|&(f, _)| crate::refinable::address(f, &t.start))).collect() };
    let (moved_a, moved_b) = (moved(a.path), moved(b.path));
    let fb = fields_by_address(&b.path.start);
    let fixed_differences: Vec<(String, f64, f64)> = fields_by_address(&a.path.start)
        .into_iter()
        .filter(|(n, _)| !moved_a.contains(n) && !moved_b.contains(n))
        .filter_map(|(n, va)| {
            let vb = fb.iter().find(|x| x.0 == n)?.1;
            ((va - vb).abs() > 1e-12 * va.abs().max(vb.abs())).then_some((n, va, vb))
        })
        .collect();
    let units: Vec<String> = ja
        .iter()
        .map(|&j| {
            let p = &a.path.set.free[j];
            if p.value.is_some() {
                String::new()
            } else {
                crate::parameters::lookup(crate::refinable::dictionary_id(p.fields[0].0, &a.path.start)).map_or(String::new(), |d| d.unit.to_string())
            }
        })
        .collect();
    type Picked = (Vec<Vec<f64>>, Vec<Vec<Option<f64>>>, Vec<f64>);
    let pick = |t: &RunTable, js: &[usize]| -> Picked {
        (t.values.iter().map(|row| js.iter().map(|&j| row[j]).collect()).collect(), t.esd.iter().map(|row| js.iter().map(|&j| row[j]).collect()).collect(), js.iter().map(|&j| t.end_esd.get(j).copied().unwrap_or(0.0)).collect())
    };
    let (values_a, esd_a, end_esd_a) = pick(&ta, &ja);
    let (values_b, esd_b, end_esd_b) = pick(&tb, &jb);
    let (end_a, end_b) = (values_a[values_a.len() - 1].clone(), values_b[values_b.len() - 1].clone());
    let (origin, scale) = match reference {
        Side::A => (end_a.clone(), end_esd_a.clone()),
        Side::B => (end_b.clone(), end_esd_b.clone()),
    };
    let scaled = |vs: &[Vec<f64>]| -> Vec<Vec<Option<f64>>> { vs.iter().map(|row| row.iter().enumerate().map(|(j, &v)| (scale[j] > 0.0).then(|| (v - origin[j]) / scale[j])).collect()).collect() };
    let (a_scaled, b_scaled) = (scaled(&values_a), scaled(&values_b));
    let (na, nb) = (values_a.len(), values_b.len());
    let steps: Vec<(usize, usize)> = (0..na.max(nb)).map(|k| (k.min(na - 1), k.min(nb - 1))).collect();
    let held_a = (0..steps.len()).map(|k| k > na - 1).collect();
    let held_b = (0..steps.len()).map(|k| k > nb - 1).collect();
    let difference = steps.iter().map(|&(ka, kb)| (0..common.len()).map(|j| (scale[j] > 0.0).then(|| (values_b[kb][j] - values_a[ka][j]) / scale[j])).collect()).collect();
    let end: Vec<EndDifference> = (0..common.len())
        .map(|j| {
            let delta = end_b[j] - end_a[j];
            // √(σA² + σB²) as products, a sum and a root, each correctly rounded (IEEE 754) and the sum commutative:
            // the same bits for A against B and B against A on every platform. `f64::hypot` is the C library's, which
            // need not be symmetric in its last bit: macOS's gave hypot(σA, σB) ≠ hypot(σB, σA) (L187, L192).
            let (ea, eb) = (end_esd_a[j], end_esd_b[j]);
            let s = (ea * ea + eb * eb).sqrt();
            let combined = (s > 0.0).then(|| delta / s);
            EndDifference { address: common[j].clone(), delta, in_scale: (scale[j] > 0.0).then(|| delta / scale[j]), combined, agree: combined.map_or(delta == 0.0, |c| c.abs() <= ENDS_AGREE) }
        })
        .collect();
    let agree = end.iter().filter(|e| e.agree).count();
    let largest = end.iter().filter_map(|e| e.combined.map(|c| (e.address.clone(), c))).max_by(|x, y| x.1.abs().total_cmp(&y.1.abs()));
    let sub = |c: &Option<Vec<Vec<f64>>>, js: &[usize]| c.as_ref().map(|m| js.iter().map(|&j| js.iter().map(|&k| m[j][k]).collect()).collect());
    let summary = CompareSummary {
        common: common.len(),
        only_a: only_a.len(),
        only_b: only_b.len(),
        fixed_differences: fixed_differences.len(),
        agree,
        disagree: end.len() - agree,
        largest,
        last_chi2_a: ta.chi2.last().copied().filter(|v| v.is_finite()),
        last_chi2_b: tb.chi2.last().copied().filter(|v| v.is_finite()),
        stop_a: a.stop.map(|(c, t)| (c.to_string(), t.to_string())),
        stop_b: b.stop.map(|(c, t)| (c.to_string(), t.to_string())),
    };
    PathsCompared {
        reference,
        common,
        only_a,
        only_b,
        fixed_differences,
        units,
        states_a: ta.states.clone(),
        states_b: tb.states.clone(),
        chi2_a: ta.chi2.clone(),
        chi2_b: tb.chi2.clone(),
        values_a,
        values_b,
        esd_a,
        esd_b,
        end_esd_a,
        end_esd_b,
        origin,
        scale,
        a_scaled,
        b_scaled,
        steps,
        held_a,
        held_b,
        difference,
        end,
        summary,
        end_correlations_a: sub(&ta.correlations, &ja),
        end_correlations_b: sub(&tb.correlations, &jb),
    }
}

impl PathsCompared {
    /// A pair of common parameters in both runs (E11a).
    pub fn pair(&self, a: &str, b: &str) -> Result<PairCompared, String> {
        let at = |n: &str| self.common.iter().position(|c| c == n).ok_or(format!("{n} is not free in both runs"));
        let (i, k) = (at(a)?, at(b)?);
        let path = |vs: &[Vec<f64>]| vs.iter().map(|row| (row[i], row[k])).collect::<Vec<_>>();
        let (path_a, path_b) = (path(&self.values_a), path(&self.values_b));
        let r = |c: &Option<Vec<Vec<f64>>>| c.as_ref().map(|m| m[i][k]);
        let (correlation_a, correlation_b) = (r(&self.end_correlations_a), r(&self.end_correlations_b));
        let (sigmas_a, sigmas_b) = ((self.end_esd_a[i], self.end_esd_a[k]), (self.end_esd_b[i], self.end_esd_b[k]));
        Ok(PairCompared {
            a: a.to_string(),
            b: b.to_string(),
            centre_a: *path_a.last().unwrap(),
            centre_b: *path_b.last().unwrap(),
            path_a,
            path_b,
            sigmas_a,
            sigmas_b,
            correlation_a,
            correlation_b,
            ellipse_a: correlation_a.and_then(|r| covariance_ellipse(sigmas_a.0, sigmas_a.1, r)),
            ellipse_b: correlation_b.and_then(|r| covariance_ellipse(sigmas_b.0, sigmas_b.1, r)),
        })
    }
}
