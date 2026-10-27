//! Joint refinement (P06): several patterns of the same sample, X-rays and neutrons, several instruments, refined
//! together, as FullProf's multi-pattern files hold them (FullProf_Manual.pdf pp. 64–66, LINE2–3 NPATT and W_PAT;
//! p. 88, LINE 19-2, the contributions of each phase to each pattern).
//!
//! Each pattern keeps its own model and data, as a single-pattern refinement has them: its instrument (wavelengths,
//! zero, background, excluded regions, Lp), and of each phase that contributes to it its own scale, profile, cell and
//! asymmetry, the structure (atoms) shared. A parameter of the joint refinement moves the fields of every pattern it
//! names: a `.pcr`'s codewords tie the patterns where they share a parameter number (pbb.pcr's cell, 351 361 371 in
//! both patterns; its atoms, written once), and the multipliers are each pattern's own.
//!
//! The quantity minimised is Q = Σ_k W_k S_k, S_k = Σ w (y − y_c)² over pattern k's points and W_k its W_PAT
//! (normalised to sum 1, manual p. 64): the patterns' normal equations added with the factors W_k, solved by the
//! single-pattern engine's Marquardt step (`lsq::Solver`). It is FullProf's: from the starting values of its run pbb,
//! where FullProf's first cycle shifts nothing, the engine's minimum of Q lies within 0.25 of FullProf's esd from
//! FullProf's values for 51 of the 63 parameters, every atom's within 0.11, where weighting the patterns' χ²,
//! Σ_k W_k S_k/(N_k − P_k), puts the atoms up to 2.8 esd away (Pb's β33; test `joint_refinement`, L322).
//!
//! The esd's (`Esd`): by default √(C_jj Q/(N − P)), C the inverse of the joint normal matrix, N the points of the
//! patterns refined and P the joint parameters: the points weighted by W_k w as known up to one scale, which for one
//! pattern is the single-pattern engine's √(C_jj χ²). FullProf's are √(C_jj Σ_k W_k χ²_k), the patterns' χ² over
//! their points with Bragg contributions, its single-pattern rule (`fullprof_refinements`) added over the patterns
//! with their weights: for patterns of equal weights about √NPATT times the default (pbb: 1/0.694, reproduced to 0.1%
//! for every parameter of the neutron pattern).
//!
//! FullProf's global χ² is the patterns' χ² over their points with Bragg contributions weighted by W_PAT ("Global
//! user-weigthed Chi2 (Bragg contrib.)": pbb_calc.out, 0.5·3.34 + 0.5·3.00 = 3.17); `Statistics::chi2` gives it with
//! the patterns' χ² over all their points beside it.
//!
//! Not yet in a joint refinement: Le Bail phases, restraints, the steering of the single-pattern engine (damping,
//! maximum shifts, bounds), plans of steps.

use crate::data::Pattern;
use crate::jacobian::Field;
use crate::lsq::{self, Options, Statistics};
use crate::pcr::MultiPcr;
use crate::powder::{Calculated, Model};
use crate::refinable::ParameterSet;
use std::path::Path;

/// One pattern of a joint refinement.
#[derive(Clone, Debug)]
pub struct Member {
    /// Its name: the data file's name without its folder and extension.
    pub label: String,
    /// The data file as the control file names it.
    pub file: String,
    pub model: Model,
    pub data: Pattern,
    /// W_PAT, normalised to sum 1 over the patterns.
    pub weight: f64,
    /// Refined (FullProf's flag 1), or calculated only (0, excluded).
    pub refined: bool,
    /// The joint phase (0, 1, …) of each of the model's phases.
    pub phases: Vec<usize>,
    /// The pattern's free parameters, and FullProf's parameter number of each (`refinable::pcr_numbered_set`; None
    /// for a pattern read from a control file, whose patterns the joint control file ties).
    pub set: ParameterSet,
    pub numbers: Vec<Option<u32>>,
    /// What the pattern's own control file holds besides the model: its title, its refinement as the control file
    /// writes it (`set` is built from it), its data file with its reader and operations (relative to `Joint::dir`),
    /// what an imported `.pcr` wrote for the export, and its background points' file.
    pub title: String,
    pub refinement: crate::refinable::Refinement,
    pub source: Option<crate::control::DataSource>,
    pub fullprof: crate::control::FullProfRecord,
    pub background_file: Option<String>,
}

/// A parameter of the joint refinement: the free parameters of the patterns it moves, (pattern, index into its set).
#[derive(Clone, Debug, PartialEq)]
pub struct Parameter {
    /// The first pattern's name of it, with "pattern{k}." before it when it moves one pattern of several.
    pub name: String,
    /// FullProf's parameter number, when the patterns were tied by it.
    pub number: Option<u32>,
    pub members: Vec<(usize, usize)>,
}

/// A joint refinement.
#[derive(Clone, Debug)]
pub struct Joint {
    pub title: String,
    /// The phases' names (each phase once, whatever patterns it contributes to).
    pub phases: Vec<String>,
    pub patterns: Vec<Member>,
    pub parameters: Vec<Parameter>,
    /// The cycles and convergence the control file asks for (its first pattern's NCY and Eps).
    pub cycles: Option<usize>,
    pub eps: Option<f64>,
    /// The folder the patterns' files are named from.
    pub dir: std::path::PathBuf,
}

/// The esd's of a joint refinement.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Esd {
    /// √(C_jj Q/(N − P)) (the default).
    Default,
    /// FullProf's: √(C_jj Σ_k W_k χ²_k,Bragg).
    FullProf,
}

/// The agreement of one pattern: over its points, over its points with Bragg contributions, and the joint
/// parameters that move it.
#[derive(Clone, Debug, PartialEq)]
pub struct PatternStatistics {
    pub all: Statistics,
    pub bragg: Statistics,
    pub parameters: usize,
}

/// The agreement of a joint refinement.
#[derive(Clone, Debug, PartialEq)]
pub struct JointStatistics {
    pub patterns: Vec<PatternStatistics>,
    /// Σ_k W_k χ²_k over the points with Bragg contributions (FullProf's global χ²), and over all the points.
    pub chi2: f64,
    pub chi2_all: f64,
    /// R_p, R_wp, R_exp (%) over every pattern's points, each pattern's sums weighted by W_k:
    /// R_wp = √(Σ W_k S_k / Σ W_k Σ w y²), R_exp = √(Σ W_k (N_k − P_k) / Σ W_k Σ w y²), R_p = Σ W_k Σ|y − y_c| / Σ W_k Σ y.
    pub rp: f64,
    pub rwp: f64,
    pub rexp: f64,
}

/// One cycle of a joint refinement.
#[derive(Clone, Debug, PartialEq)]
pub struct JointCycle {
    pub cycle: usize,
    /// The global χ² (Bragg points) and the objective at the start of the cycle and after its step.
    pub chi2_start: f64,
    pub chi2: f64,
    pub objective: f64,
    pub marquardt: f64,
    /// Per parameter: value after the step, shift, esd.
    pub params: Vec<(String, f64, f64, f64)>,
}

/// The end of a joint refinement.
#[derive(Clone, Debug, PartialEq)]
pub struct JointRefinement {
    pub cycles: usize,
    pub converged: bool,
    /// Why it stopped, as text and as data.
    pub stop: String,
    pub reason: lsq::Stop,
    /// Per parameter: name, value, esd.
    pub params: Vec<(String, f64, f64)>,
    /// The covariance of the parameters, C Q (`esd = √diag`).
    pub covariance: Vec<Vec<f64>>,
    pub statistics: JointStatistics,
    pub history: Vec<JointCycle>,
}

fn label_of(file: &str) -> String {
    Path::new(file.rsplit(['\\', '/']).next().unwrap_or(file)).file_stem().map_or(file.to_string(), |s| s.to_string_lossy().to_string())
}

impl Joint {
    /// A joint refinement from a multi-pattern `.pcr` (`pcr::split_patterns`; a single-pattern file gives a joint
    /// refinement of its one pattern, of weight 1): each pattern's model and data as the
    /// single-pattern import builds them (`Model::from_pcr`, or with `fullprof` true `Model::from_pcr_fullprof`: the
    /// values under FullProf's rules), its free parameters from its codewords, the patterns tied by their numbers.
    /// `dir` is the file's folder, `stem` its name without extension.
    pub fn from_pcr(text: &str, dir: &Path, stem: &str, fullprof: bool) -> Result<Joint, String> {
        let m: MultiPcr = if crate::pcr::is_multipattern(text) {
            crate::pcr::split_patterns(text).map_err(|e| e.to_string())?
        } else {
            // A single-pattern file: a joint refinement of one pattern, weight 1.
            let p = crate::pcr::read(text).map_err(|e| e.to_string())?;
            let file = p.data_file.clone().unwrap_or_else(|| format!("{stem}.dat"));
            MultiPcr { title: p.title.clone(), refined: vec![true], weights: vec![1.0], files: vec![file], texts: vec![text.to_string()], phases: p.phases.iter().map(|ph| ph.name.clone()).collect(), contributions: vec![vec![true]; p.phases.len()] }
        };
        let mut patterns = Vec::new();
        let (mut cycles, mut eps) = (None, None);
        for (k, t) in m.texts.iter().enumerate() {
            let what = |e: String| format!("pattern {} ({}): {e}", k + 1, m.files[k]);
            let p = crate::pcr::read(t).map_err(|e| what(e.to_string()))?;
            if k == 0 {
                (cycles, eps) = (Some(p.ncy.max(0) as usize), Some(p.eps));
            }
            let l = if fullprof { Model::from_pcr_fullprof(&p, dir, stem) } else { Model::from_pcr(&p, dir, stem) }.map_err(what)?;
            let data = l.data.ok_or_else(|| what(format!("the data file {} is not found", m.files[k])))?;
            if l.model.has_extracted() {
                return Err(what("a Le Bail phase in a joint refinement is not refined yet".into()));
            }
            let (set, numbers) = crate::refinable::pcr_numbered_set(&p, &l.model).map_err(what)?;
            let phases = (0..m.phases.len()).filter(|&q| m.contributions[q][k]).collect();
            // The control file's parts of the pattern, as `drs convert` imports a .pcr (its data relative to `dir`).
            let c = crate::control::Control::from_pcr(&p, dir, stem, dir).map_err(what)?;
            let refinement = crate::refinable::Refinement::from_pcr(&p, &l.model).map_err(what)?;
            patterns.push(Member {
                label: label_of(&m.files[k]),
                file: m.files[k].clone(),
                model: l.model,
                data,
                weight: m.weights[k],
                refined: m.refined[k],
                phases,
                set,
                numbers,
                title: c.title,
                refinement,
                source: c.data,
                fullprof: c.fullprof,
                background_file: c.background_file,
            });
        }
        let parameters = tie(&patterns);
        Ok(Joint { title: m.title, phases: m.phases, patterns, parameters, cycles, eps, dir: dir.to_path_buf() })
    }

    /// Each pattern's calculated pattern at its data's points.
    pub fn calculate(&self) -> Vec<Calculated> {
        self.patterns.iter().map(|p| p.model.calculate(&p.data.x)).collect()
    }

    /// The joint parameters that move pattern k.
    pub fn parameters_of(&self, k: usize) -> usize {
        self.parameters.iter().filter(|q| q.members.iter().any(|m| m.0 == k)).count()
    }

    /// The agreement at the models' values.
    pub fn statistics(&self) -> JointStatistics {
        let mut pats = Vec::new();
        let (mut chi2, mut chi2_all) = (0.0, 0.0);
        let (mut sa, mut sy, mut swd, mut swy, mut dof) = (0.0, 0.0, 0.0, 0.0, 0.0);
        for (k, p) in self.patterns.iter().enumerate() {
            let c = p.model.calculate(&p.data.x);
            let np = self.parameters_of(k);
            let all = lsq::statistics(&p.model, &p.data, &c.ycalc, np);
            let mask = lsq::bragg_points(&p.data, &c.peaks);
            let bragg = lsq::statistics_over(&p.model, &p.data, &c.ycalc, np, Some(&mask));
            chi2 += p.weight * bragg.chi2;
            chi2_all += p.weight * all.chi2;
            let w = lsq::weights(&p.model, &p.data, &c.ycalc);
            for ((&wi, &y), &yc) in w.iter().zip(&p.data.y).zip(&c.ycalc) {
                if wi > 0.0 {
                    let d = y - yc;
                    sa += p.weight * d.abs();
                    sy += p.weight * y.abs();
                    swd += p.weight * wi * d * d;
                    swy += p.weight * wi * y * y;
                }
            }
            dof += p.weight * (all.n.saturating_sub(np)).max(1) as f64;
            pats.push(PatternStatistics { all, bragg, parameters: np });
        }
        JointStatistics { patterns: pats, chi2, chi2_all, rp: 100.0 * sa / sy.max(1e-300), rwp: 100.0 * (swd / swy.max(1e-300)).sqrt(), rexp: 100.0 * (dof / swy.max(1e-300)).sqrt() }
    }

    /// Each joint parameter's value (its first member's).
    pub fn values(&self) -> Vec<f64> {
        self.parameters.iter().map(|q| { let (k, i) = q.members[0]; self.patterns[k].set.free[i].value_in(&self.patterns[k].model) }).collect()
    }

    fn apply(&mut self, delta: &[f64]) {
        for (q, &d) in self.parameters.iter().zip(delta) {
            for &(k, i) in &q.members {
                let p = &mut self.patterns[k];
                for &(f, m) in &p.set.free[i].fields {
                    p.model.set(f, p.model.get(f) + m * d);
                }
            }
        }
    }

    /// A pattern's S = Σ w (y − y_c)² at its model's values (refined or not).
    pub fn objective_of_pattern(&self, p: &Member) -> f64 {
        let yc = p.model.calculate(&p.data.x).ycalc;
        lsq::objective(&p.data, &yc, &lsq::weights(&p.model, &p.data, &yc))
    }

    /// The objective Q = Σ_k W_k S_k over the patterns refined, at the models' values.
    pub fn objective(&self) -> f64 {
        self.patterns
            .iter()
            .filter(|p| p.refined)
            .map(|p| {
                let yc = p.model.calculate(&p.data.x).ycalc;
                p.weight * lsq::objective(&p.data, &yc, &lsq::weights(&p.model, &p.data, &yc))
            })
            .sum()
    }

    /// Refine: Marquardt cycles over the patterns refined (FullProf's flag 1), up to `opts.max_cycles` or until every
    /// |shift|/esd is below `opts.eps`; `on_cycle` receives each cycle and stops the run when it returns false.
    pub fn refine(&mut self, opts: Options, esd: Esd, mut on_cycle: impl FnMut(&JointCycle) -> bool) -> JointRefinement {
        let np = self.parameters.len();
        let mut mu = self.patterns.first().map_or(1e-3, |p| p.set.steering.lambda);
        let mut history = Vec::new();
        let (mut converged, mut stop) = (false, String::new());
        let mut reason = lsq::Stop::Cycles { n: 0, le_bail: false };
        let mut cov = vec![vec![0.0; np]; np];
        let mut cycle = 0;
        while cycle < opts.max_cycles {
            cycle += 1;
            let (a, g, q0) = self.normal_equations();
            let solver = lsq::Solver::new(&a, (0..np).filter(|&j| a[j][j] > 0.0).collect());
            let chi2_start = self.statistics().chi2;
            let start: Vec<Model> = self.patterns.iter().map(|p| p.model.clone()).collect();
            let mut accepted = None;
            for _ in 0..13 {
                let delta = solver.solve(&g, mu, np);
                self.apply(&delta);
                let q = self.objective();
                if q.is_finite() && q <= q0 {
                    accepted = Some((delta, q));
                    mu = (mu / 10.0).max(1e-12);
                    break;
                }
                for (p, m) in self.patterns.iter_mut().zip(&start) {
                    p.model = m.clone();
                }
                mu *= 10.0;
            }
            let Some((delta, q)) = accepted else {
                // At the minimum to the objective's rounding: the step of the smallest λ, against the esd's.
                cov = self.covariance(esd);
                let delta = solver.solve(&g, 0.0, np);
                cycle -= 1;
                converged = (0..np).all(|j| cov[j][j] <= 0.0 || (delta[j] / cov[j][j].sqrt()).abs() < opts.eps);
                reason = if converged { lsq::Stop::Converged { eps: opts.eps } } else { lsq::Stop::Minimum };
                stop = if converged {
                    format!("converged after {cycle} cycles: no step lowers the objective, and every |shift|/esd is below {}", opts.eps)
                } else {
                    format!("no step lowers the objective after cycle {cycle}")
                };
                break;
            };
            // The esd's at the values reached.
            cov = self.covariance(esd);
            let values = self.values();
            let esd: Vec<f64> = (0..np).map(|j| cov[j][j].max(0.0).sqrt()).collect();
            let params = (0..np).map(|j| (self.parameters[j].name.clone(), values[j], delta[j], esd[j])).collect();
            let c = JointCycle { cycle, chi2_start, chi2: self.statistics().chi2, objective: q, marquardt: mu, params };
            let go_on = on_cycle(&c);
            history.push(c);
            if (0..np).all(|j| esd[j] == 0.0 || (delta[j] / esd[j]).abs() < opts.eps) {
                converged = true;
                reason = lsq::Stop::Converged { eps: opts.eps };
                stop = format!("converged after {cycle} cycles: every |shift|/esd below {}", opts.eps);
                break;
            }
            if !go_on {
                reason = lsq::Stop::User;
                stop = lsq::STOPPED.to_string();
                break;
            }
        }
        if stop.is_empty() {
            reason = lsq::Stop::Cycles { n: cycle, le_bail: false };
            stop = format!("{cycle} cycles run (the cycles asked for)");
        }
        if history.is_empty() {
            cov = self.covariance(esd);
        }
        let values = self.values();
        let params = (0..np).map(|j| (self.parameters[j].name.clone(), values[j], cov[j][j].max(0.0).sqrt())).collect();
        JointRefinement { cycles: cycle, converged, stop, reason, params, covariance: cov, statistics: self.statistics(), history }
    }

    /// The joint normal matrix, gradient and objective at the models' values.
    fn normal_equations(&self) -> (Vec<Vec<f64>>, Vec<f64>, f64) {
        let np = self.parameters.len();
        let mut a = vec![vec![0.0; np]; np];
        let mut g = vec![0.0; np];
        let mut q = 0.0;
        for (k, p) in self.patterns.iter().enumerate().filter(|(_, p)| p.refined) {
            let fields = p.set.fields();
            let (yc, jq) = lsq::jacobian(&p.model, &p.data, &p.set, &fields);
            let w = lsq::weights(&p.model, &p.data, &yc);
            q += p.weight * lsq::objective(&p.data, &yc, &w);
            // Each joint parameter's column in this pattern: the sum of its members' columns here.
            let cols: Vec<(usize, Vec<f64>)> = self
                .parameters
                .iter()
                .enumerate()
                .filter_map(|(j, par)| {
                    let mine: Vec<usize> = par.members.iter().filter(|m| m.0 == k).map(|m| m.1).collect();
                    if mine.is_empty() {
                        return None;
                    }
                    let mut col = vec![0.0; w.len()];
                    for i in mine {
                        for (c, d) in col.iter_mut().zip(&jq[i]) {
                            *c += d;
                        }
                    }
                    Some((j, col))
                })
                .collect();
            for i in 0..w.len() {
                if w[i] == 0.0 {
                    continue;
                }
                let r = p.data.y[i] - yc[i];
                let wi = p.weight * w[i];
                for (x, (j, cj)) in cols.iter().enumerate() {
                    let wj = wi * cj[i];
                    if wj == 0.0 {
                        continue;
                    }
                    g[*j] += wj * r;
                    for (l, ck) in &cols[x..] {
                        a[*j][*l] += wj * ck[i];
                    }
                }
            }
        }
        lsq::symmetrize(&mut a);
        (a, g, q)
    }

    /// The pattern (0, 1, …) and field FullProf's symbolic name names in a joint refinement (`output_names`): `_pat<k>`
    /// names the pattern, `_ph<n>` the joint phase; a name without `_pat` (an atom's, X_Pb_ph1) the first pattern its
    /// phase contributes to, whose field the joint parameter moves with the others'.
    pub fn fullprof_field(&self, name: &str) -> Option<(usize, Field)> {
        let n = name.trim();
        let (base, k) = match n.rfind("_pat") {
            Some(i) => (&n[..i], n[i + 4..].parse::<usize>().ok()?.checked_sub(1)?),
            None => (n, usize::MAX),
        };
        let (stem, phase) = match base.rfind("_ph") {
            Some(i) => (&base[..i], Some(base[i + 3..].parse::<usize>().ok()?.checked_sub(1)?)),
            None => (base, None),
        };
        let k = match (k, phase) {
            (usize::MAX, Some(g)) => self.patterns.iter().position(|p| p.phases.contains(&g))?,
            (usize::MAX, None) => 0,
            (k, _) => k,
        };
        let p = self.patterns.get(k)?;
        let local = match phase {
            Some(g) => format!("{stem}_ph{}", p.phases.iter().position(|&x| x == g)? + 1),
            None => stem.to_string(),
        };
        crate::output_names::fullprof_field(&local, &p.model).map(|f| (k, f))
    }

    /// The fields of a joint parameter, by pattern: (pattern, field, multiplier).
    pub fn fields_of(&self, j: usize) -> Vec<(usize, Field, f64)> {
        self.parameters[j].members.iter().flat_map(|&(k, i)| self.patterns[k].set.free[i].fields.iter().map(move |&(f, m)| (k, f, m))).collect()
    }
}

fn scaled(inv: &[Vec<f64>], q: f64) -> Vec<Vec<f64>> {
    inv.iter().map(|r| r.iter().map(|v| v * q).collect()).collect()
}

impl Joint {
    /// N − P of the refinement: the points of the patterns refined, less the joint parameters.
    pub fn degrees_of_freedom(&self) -> usize {
        let n: usize = self
            .patterns
            .iter()
            .filter(|p| p.refined)
            .map(|p| {
                let yc = p.model.calculate(&p.data.x).ycalc;
                lsq::weights(&p.model, &p.data, &yc).iter().filter(|&&w| w > 0.0).count()
            })
            .sum();
        n.saturating_sub(self.parameters.len()).max(1)
    }

    /// The covariance at the models' values (`Esd`): C Q/(N − P), or FullProf's C Σ_k W_k χ²_k,Bragg.
    pub fn covariance(&self, esd: Esd) -> Vec<Vec<f64>> {
        let np = self.parameters.len();
        let (a, _, q) = self.normal_equations();
        let inv = lsq::Solver::new(&a, (0..np).filter(|&j| a[j][j] > 0.0).collect()).inverse(np);
        let factor = match esd {
            Esd::Default => q / self.degrees_of_freedom() as f64,
            Esd::FullProf => self.statistics().chi2,
        };
        scaled(&inv, factor)
    }
}

/// A pattern's name of a parameter in the joint refinement's: its phases numbered as the joint refinement's
/// (`Member::phases`), and, when it moves this pattern only of several, "pattern{k}." for the single pattern's
/// "pattern1." or before a phase's name.
fn joint_name(name: &str, k: usize, phases: &[usize], own: bool) -> String {
    let renumbered = match name.strip_prefix("phase").and_then(|r| r.split_once('.')) {
        Some((n, rest)) => match n.parse::<usize>().ok().and_then(|n| phases.get(n.wrapping_sub(1))) {
            Some(g) => format!("phase{}.{rest}", g + 1),
            None => name.to_string(),
        },
        None => name.to_string(),
    };
    if !own {
        return renumbered;
    }
    match renumbered.strip_prefix("pattern1.") {
        Some(rest) => format!("pattern{}.{rest}", k + 1),
        None => format!("pattern{}.{renumbered}", k + 1),
    }
}

/// The joint parameters: the patterns' free parameters with the same FullProf number are one; the others their own.
fn tie(patterns: &[Member]) -> Vec<Parameter> {
    let mut out: Vec<Parameter> = Vec::new();
    for (k, p) in patterns.iter().enumerate() {
        for (i, f) in p.set.free.iter().enumerate() {
            match p.numbers[i].and_then(|n| out.iter_mut().find(|q| q.number == Some(n))) {
                Some(q) => q.members.push((k, i)),
                None => out.push(Parameter { name: f.name.clone(), number: p.numbers[i], members: vec![(k, i)] }),
            }
        }
    }
    for q in out.iter_mut() {
        let ks: std::collections::BTreeSet<usize> = q.members.iter().map(|m| m.0).collect();
        let (k, _) = q.members[0];
        q.name = joint_name(&q.name, k, &patterns[k].phases, patterns.len() > 1 && ks.len() == 1);
    }
    out
}

/// What `read_results` reads of a joint results file's pattern block.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PatternResults {
    /// Its number (1, 2, …), name (the data file's), `_pd_block_id`, W_PAT and whether it was refined.
    pub id: usize,
    pub label: String,
    pub block: String,
    pub weight: f64,
    pub refined: bool,
    /// Its agreement: R_p, R_wp, R_exp (%), χ² over its points and over those with Bragg contributions.
    pub rp: f64,
    pub rwp: f64,
    pub rexp: f64,
    pub chi2: f64,
    pub chi2_bragg: f64,
    /// The joint phases (1, 2, …) that contribute to it, the profile's points and the reflections (phase, pattern).
    pub phases: Vec<usize>,
    pub points: usize,
    pub reflections: Vec<(usize, usize)>,
    /// Its abscissa: 2θ, or a TOF bank's time of flight (P46: a multi-bank refinement is a joint one, a block per bank).
    pub axis: crate::formats::Axis,
    /// A TOF bank's angle 2θ (deg).
    pub two_theta_fixed: Option<f64>,
    /// Its control file's text (`Joint::pattern_control`), empty in a file written before core 0.61.0 (P50).
    pub control_file: String,
}

/// What `read_results` reads of a joint results file.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct JointResults {
    pub patterns: Vec<PatternResults>,
    /// The joint agreement: χ² (Σ W_k χ²_k over all points, and over the Bragg points: FullProf's global χ²), R_p,
    /// R_wp, R_exp (%), the objective Q, the points used and the parameters.
    pub chi2: f64,
    pub chi2_bragg: f64,
    pub rp: f64,
    pub rwp: f64,
    pub rexp: f64,
    pub objective: f64,
    pub points: usize,
    pub parameters: usize,
    /// The refined parameters: name, value, su.
    pub refined: Vec<(String, f64, f64)>,
    /// The phases' blocks: (phase id, `_pd_block_id`).
    pub phase_blocks: Vec<(usize, String)>,
    /// The overall block's `_drs_control_file`: the source the refinement ran, a multi-pattern `.pcr` or a joint
    /// control file; and the name the file was written under (`_pd_block_id` before `|overall`).
    pub control_file: String,
    pub name: String,
}

impl JointResults {
    /// Each pattern's control file, in the patterns' order: its model (with its phases, `PatternResults::phases`
    /// naming them in the joint refinement), refinement and data source, read without the data (Publish's P50). From
    /// each pattern block's own control file; in a file written before core 0.61.0, from the overall block's source:
    /// a multi-pattern `.pcr` split by pattern and imported (`dir`, the folder of the `.pcr`'s files, is where its
    /// resolution and hkl files are looked for), or a joint control file's patterns' control files read from `dir`.
    pub fn pattern_controls(&self, dir: &Path) -> Result<Vec<crate::control::Control>, String> {
        use crate::control::Control;
        let diag = |d: Vec<crate::control::syntax::Diagnostic>| d.iter().map(|x| format!("line {}: {}", x.span.line, x.message)).collect::<Vec<_>>().join("; ");
        if self.patterns.iter().all(|p| !p.control_file.is_empty()) {
            return self.patterns.iter().map(|p| Control::from_text(&p.control_file).map_err(|d| format!("pattern {}: {}", p.id, diag(d)))).collect();
        }
        let src = &self.control_file;
        if src.trim().is_empty() {
            return Err("the file holds no control file".into());
        }
        if is_joint_control(src) {
            let doc = crate::control::syntax::parse(src);
            return crate::control::relocate::paths(&doc)
                .into_iter()
                .filter(|p| p.table == "joint.pattern")
                .map(|p| {
                    let f = dir.join(&p.path);
                    let text = std::fs::read_to_string(&f).map_err(|e| format!("{}: {e}", f.display()))?;
                    Control::from_text_in(&text, dir).map_err(|d| format!("{}: {}", f.display(), diag(d)))
                })
                .collect();
        }
        let texts = if crate::pcr::is_multipattern(src) { crate::pcr::split_patterns(src).map_err(|e| e.to_string())?.texts } else { vec![src.clone()] };
        let stem = self.name.rsplit_once(' ').map_or(self.name.as_str(), |(s, _)| s);
        texts
            .iter()
            .enumerate()
            .map(|(k, t)| {
                let p = crate::pcr::read(t).map_err(|e| format!("pattern {}: {e}", k + 1))?;
                Control::from_pcr(&p, dir, stem, dir).map_err(|e| format!("pattern {}: {e}", k + 1))
            })
            .collect()
    }
}

impl Joint {
    /// The joint parameter of pattern k's free parameter i.
    fn joint_index(&self, k: usize, i: usize) -> usize {
        self.parameters.iter().position(|q| q.members.contains(&(k, i))).expect("every free parameter is in a joint one")
    }

    /// Pattern k's part of a joint refinement as the single-pattern engine reports a refinement (`lsq::Refinement`),
    /// for the writers that take one (the phases' CIF blocks, the pattern block): its free parameters with the joint
    /// covariance (`JointRefinement::covariance`) restricted to them, each field's su from it, its own agreement and
    /// Bérar factor, the run's cycles and stop.
    pub fn pattern_refinement(&self, k: usize, r: &JointRefinement) -> lsq::Refinement {
        let p = &self.patterns[k];
        let np = p.set.free.len();
        let idx: Vec<usize> = (0..np).map(|i| self.joint_index(k, i)).collect();
        let cov: Vec<Vec<f64>> = (0..np).map(|a| (0..np).map(|b| r.covariance.get(idx[a]).and_then(|row| row.get(idx[b])).copied().unwrap_or(0.0)).collect()).collect();
        let correlations = (0..np).map(|a| (0..np).map(|b| if cov[a][a] > 0.0 && cov[b][b] > 0.0 { cov[a][b] / (cov[a][a] * cov[b][b]).sqrt() } else { 0.0 }).collect()).collect();
        let params = p
            .set
            .free
            .iter()
            .enumerate()
            .map(|(i, f)| lsq::ParamState { name: f.name.clone(), value: f.value_in(&p.model), shift: 0.0, esd: f.esd_of(cov[i][i].max(0.0).sqrt()), shift_esd: 0.0, computed: 0.0, at_limit: false })
            .collect();
        let fields = p
            .set
            .fields()
            .into_iter()
            .map(|f| {
                let m: Vec<f64> = p.set.free.iter().map(|q| q.fields.iter().filter(|x| x.0 == f).map(|x| x.1).sum()).collect();
                let var: f64 = (0..np).map(|a| (0..np).map(|b| m[a] * cov[a][b] * m[b]).sum::<f64>()).sum();
                (crate::refinable::address(f, &p.model), p.model.get(f), var.max(0.0).sqrt())
            })
            .collect();
        let c = p.model.calculate(&p.data.x);
        let stats = r.statistics.patterns.get(k).map_or_else(|| lsq::statistics(&p.model, &p.data, &c.ycalc, np), |s| s.all.clone());
        let last = r.history.last();
        let shift_su_max = last.map_or(0.0, |c| c.params.iter().filter(|q| q.3 > 0.0).map(|q| (q.2 / q.3).abs()).fold(0.0, f64::max));
        lsq::Refinement {
            cycles: r.cycles,
            converged: r.converged,
            params,
            stats,
            correlations,
            covariance: cov,
            fields,
            restraints: Vec::new(),
            stop: r.stop.clone(),
            reason: r.reason.clone(),
            shift_su_max,
            history: Vec::new(),
            pauses: Vec::new(),
            cycle_steps: Vec::new(),
            berar_factor: lsq::berar_factor(&p.model, &p.data, &c.ycalc, &c.peaks),
        }
    }

    /// The results file of a joint refinement (P41): a powder CIF as the single-pattern one (`results::write`) with
    /// - `data_overall` (`{name}|overall`): the joint agreement (R_p, R_wp, R_exp over every pattern's points, each
    ///   pattern's sums weighted by W_k, `JointStatistics`), `_drs_chi2` (Σ W_k χ²_k over all points),
    ///   `_drs_chi2_bragg` (over the Bragg points: FullProf's global χ²) and `_drs_joint_objective` (Q = Σ W_k S_k),
    ///   the points and parameters, the cycles and the stop; the patterns (`_drs_pattern_*`: id, block, data file,
    ///   weight, refined), the phases with their blocks, which phase contributes to which pattern
    ///   (`_drs_contribution_*`), the diffractograms' blocks, the joint parameters with their values and su's (the
    ///   `Esd` the refinement used), and the source `.pcr` (`_drs_control_file`);
    /// - one block per phase (`{name}|phase n`, n the joint phase), the deposition block of the first pattern it
    ///   contributes to, its `_pd_block_diffractogram_id` a loop of the patterns it contributes to;
    /// - one block per pattern, `data_pattern_k` (`{name}|pattern k`): its weight, refined flag, agreement (R_p,
    ///   R_wp, R_exp, χ² over all its points and over the Bragg points, the parameters that move it) and each of its
    ///   phases' scale, R_Bragg, R_F and cell (`_drs_pattern_phase_*`: a phase has its own per pattern), its own control
    ///   file (`_drs_control_file`, `pattern_control`: its model read without the data, Publish's P50), then as a
    ///   single pattern's: radiation, wavelengths, excluded regions, the profile and the reflections, these with the
    ///   pattern's id first (`_drs_refln_pattern_id`) and the joint phase's id.
    ///
    /// Without `r` the file holds the values at the models' state and no uncertainties.
    pub fn results(&self, r: Option<&JointRefinement>, name: &str, source: &str) -> Result<String, String> {
        use crate::results::{num, quoted, PREFIX};
        use std::fmt::Write;
        let st = r.map_or_else(|| self.statistics(), |r| r.statistics.clone());
        let refinements: Vec<Option<lsq::Refinement>> = (0..self.patterns.len()).map(|k| r.map(|r| self.pattern_refinement(k, r))).collect();
        let version = env!("CARGO_PKG_VERSION");
        let mut out = String::new();
        let w = &mut out;
        let _ = writeln!(w, "# The results of a joint refinement of {} patterns, written by diff-ref-struct {version} (drs) as a powder CIF (R28, P41).", self.patterns.len());
        let _ = writeln!(w, "# Items the IUCr's core and powder dictionaries do not name carry the local prefix {PREFIX}.");
        let _ = writeln!(w, "\ndata_overall");
        let _ = writeln!(w, "_pd_block_id                      {}", quoted(&format!("{name}|overall")));
        let _ = writeln!(w, "_audit_creation_method            'diff-ref-struct {version} (drs)'");
        let _ = writeln!(w, "_computing_structure_refinement   'diff-ref-struct {version}'");
        let _ = writeln!(w, "_pd_proc_ls_prof_R_factor         {}", num(st.rp / 100.0));
        let _ = writeln!(w, "_pd_proc_ls_prof_wR_factor        {}", num(st.rwp / 100.0));
        let _ = writeln!(w, "_pd_proc_ls_prof_wR_expected      {}", num(st.rexp / 100.0));
        let _ = writeln!(w, "_refine_ls_goodness_of_fit_all    {}", num(st.chi2_all.sqrt()));
        let _ = writeln!(w, "{PREFIX}chi2                         {}", num(st.chi2_all));
        let _ = writeln!(w, "{PREFIX}chi2_bragg                   {}", num(st.chi2));
        let _ = writeln!(w, "{PREFIX}joint_objective              {}", num(self.objective()));
        let points: usize = st.patterns.iter().zip(&self.patterns).filter(|(_, p)| p.refined).map(|(s, _)| s.all.n).sum();
        let _ = writeln!(w, "_pd_proc_number_of_points         {points}");
        let _ = writeln!(w, "_refine_ls_number_parameters      {}", self.parameters.len());
        let _ = writeln!(w, "_refine_ls_number_restraints      0");
        if let Some(r) = r {
            let _ = writeln!(w, "{PREFIX}cycles                       {}", r.cycles);
            let _ = writeln!(w, "{PREFIX}stop_code                    {}", r.reason.code());
            let _ = writeln!(w, "{PREFIX}stop_text                    {}", quoted(&r.stop));
        }
        let _ = writeln!(w, "loop_\n{PREFIX}pattern_id\n{PREFIX}pattern_block_id\n{PREFIX}pattern_data_file\n{PREFIX}pattern_weight\n{PREFIX}pattern_refined");
        for (k, p) in self.patterns.iter().enumerate() {
            let _ = writeln!(w, "{} {} {} {} {}", k + 1, quoted(&format!("{name}|pattern {}", k + 1)), quoted(&p.file), num(p.weight), if p.refined { "yes" } else { "no" });
        }
        let _ = writeln!(w, "loop_\n_pd_phase_id\n_pd_phase_block_id\n{PREFIX}phase_name");
        for (g, n) in self.phases.iter().enumerate() {
            let _ = writeln!(w, "{} {} {}", g + 1, quoted(&format!("{name}|phase {}", g + 1)), quoted(n));
        }
        let _ = writeln!(w, "loop_\n{PREFIX}contribution_phase_id\n{PREFIX}contribution_pattern_id");
        for (k, p) in self.patterns.iter().enumerate() {
            for &g in &p.phases {
                let _ = writeln!(w, "{} {}", g + 1, k + 1);
            }
        }
        let _ = writeln!(w, "loop_\n_pd_block_diffractogram_id");
        for k in 0..self.patterns.len() {
            let _ = writeln!(w, "{}", quoted(&format!("{name}|pattern {}", k + 1)));
        }
        let _ = writeln!(w, "loop_\n{PREFIX}parameter_name\n{PREFIX}parameter_value\n{PREFIX}parameter_su");
        let values = self.values();
        for (j, q) in self.parameters.iter().enumerate() {
            let su = r.and_then(|r| r.covariance.get(j).and_then(|row| row.get(j))).map_or(f64::NAN, |v| v.max(0.0).sqrt());
            let _ = writeln!(w, "{} {} {}", quoted(&q.name), num(values[j]), if su.is_nan() { ".".to_string() } else { num(su) });
        }
        let _ = writeln!(w, "{PREFIX}control_file\n;\n{}\n;", source.trim_end().replace("\n;", "\n ;"));
        // The phases, each from the first pattern it contributes to.
        for (g, _) in self.phases.iter().enumerate() {
            let Some(k) = self.patterns.iter().position(|p| p.phases.contains(&g)) else { continue };
            let p = &self.patterns[k];
            let q = p.phases.iter().position(|&x| x == g).unwrap();
            if !matches!(p.model.phases[q].intensities, crate::powder::Intensities::Structure(..)) {
                continue;
            }
            let rk = refinements[k].as_ref();
            let dep = crate::cif_write::deposition_in(&p.model, &p.data, Some(&p.set), rk, &self.title, Some(name))?;
            let blocks = crate::results::phase_blocks(&dep);
            let before = (0..q).filter(|&x| matches!(p.model.phases[x].intensities, crate::powder::Intensities::Structure(..))).count();
            let (_, text) = &blocks[before];
            let mut lines = text.lines();
            let head = lines.next().unwrap_or("");
            let _ = writeln!(w, "\n{head}");
            for l in lines {
                if l.starts_with("_pd_block_id") {
                    let _ = writeln!(w, "_pd_block_id                      {}", quoted(&format!("{name}|phase {}", g + 1)));
                } else if l.starts_with("_pd_block_diffractogram_id") {
                    let _ = writeln!(w, "loop_\n_pd_block_diffractogram_id");
                    for (kk, pp) in self.patterns.iter().enumerate() {
                        if pp.phases.contains(&g) {
                            let _ = writeln!(w, "{}", quoted(&format!("{name}|pattern {}", kk + 1)));
                        }
                    }
                } else {
                    let _ = writeln!(w, "{l}");
                }
            }
            crate::results::modes_block(w, &p.model, q, Some(&p.set), rk);
            crate::results::structure_block(w, &p.model, q, Some(&p.set), rk);
        }
        // The patterns.
        for (k, p) in self.patterns.iter().enumerate() {
            let ps = &st.patterns[k];
            let extra = |w: &mut String| {
                let _ = writeln!(w, "{PREFIX}pattern_id                   {}", k + 1);
                let _ = writeln!(w, "{PREFIX}pattern_data_file            {}", quoted(&p.file));
                let _ = writeln!(w, "{PREFIX}pattern_weight               {}", num(p.weight));
                let _ = writeln!(w, "{PREFIX}pattern_refined              {}", if p.refined { "yes" } else { "no" });
                let _ = writeln!(w, "_pd_proc_ls_prof_R_factor         {}", num(ps.all.rp / 100.0));
                let _ = writeln!(w, "_pd_proc_ls_prof_wR_factor        {}", num(ps.all.rwp / 100.0));
                let _ = writeln!(w, "_pd_proc_ls_prof_wR_expected      {}", num(ps.all.rexp / 100.0));
                let _ = writeln!(w, "{PREFIX}chi2                         {}", num(ps.all.chi2));
                let _ = writeln!(w, "{PREFIX}chi2_bragg                   {}", num(ps.bragg.chi2));
                let _ = writeln!(w, "_pd_proc_number_of_points         {}", ps.all.n);
                let _ = writeln!(w, "{PREFIX}pattern_parameters           {}", ps.parameters);
                let _ = writeln!(w, "loop_\n{PREFIX}pattern_phase_id\n{PREFIX}pattern_phase_scale\n{PREFIX}pattern_phase_R_Bragg\n{PREFIX}pattern_phase_R_F");
                let _ = writeln!(w, "{PREFIX}pattern_phase_cell_a\n{PREFIX}pattern_phase_cell_b\n{PREFIX}pattern_phase_cell_c\n{PREFIX}pattern_phase_cell_alpha\n{PREFIX}pattern_phase_cell_beta\n{PREFIX}pattern_phase_cell_gamma");
                for (q, ph) in p.model.phases.iter().enumerate() {
                    let (rb, rf) = ps.all.phases.get(q).copied().unwrap_or((f64::NAN, f64::NAN));
                    let c = ph.cell;
                    let _ = writeln!(w, "{} {} {} {} {} {} {} {} {} {}", p.phases[q] + 1, num(ph.scale), num(rb / 100.0), num(rf / 100.0), num(c.a), num(c.b), num(c.c), num(c.alpha), num(c.beta), num(c.gamma));
                }
                // The pattern's own control file, its files named from the joint refinement's folder: its model read
                // without the data (Publish's P50).
                let _ = writeln!(w, "{PREFIX}control_file\n;\n{}\n;", self.pattern_control(k, &self.dir).to_text().trim_end().replace("\n;", "\n ;"));
            };
            crate::results::pattern_block(w, &p.model, &p.data, Some(&p.set), refinements[k].as_ref(), &format!("pattern_{}", k + 1), &format!("{name}|pattern {}", k + 1), Some(k), &p.phases, &extra);
        }
        Ok(out)
    }
}

/// A joint results file read back (`Joint::results`): its overall agreement and parameters, its phases' blocks, and
/// per pattern its block, weight, agreement, phases, points and reflections.
pub fn read_results(text: &str) -> Result<JointResults, String> {
    let blocks = crate::cif::parse(text);
    let o = blocks.iter().find(|b| b.get("_pd_block_id").is_some_and(|v| v.ends_with("|overall"))).ok_or("no overall block")?;
    let f = |b: &crate::cif::Block, tag: &str| b.get(tag).and_then(crate::cif::number).ok_or(format!("{tag} missing"));
    let col = |b: &crate::cif::Block, tag: &str| -> Result<Vec<String>, String> { Ok(b.column(tag).ok_or(format!("{tag} missing"))?.iter().map(|s| s.to_string()).collect()) };
    let n = |s: &str| crate::cif::number(s).unwrap_or(f64::NAN);
    let mut r = JointResults {
        chi2: f(o, "_drs_chi2")?,
        chi2_bragg: f(o, "_drs_chi2_bragg")?,
        rp: 100.0 * f(o, "_pd_proc_ls_prof_r_factor")?,
        rwp: 100.0 * f(o, "_pd_proc_ls_prof_wr_factor")?,
        rexp: 100.0 * f(o, "_pd_proc_ls_prof_wr_expected")?,
        objective: f(o, "_drs_joint_objective")?,
        points: f(o, "_pd_proc_number_of_points")? as usize,
        parameters: f(o, "_refine_ls_number_parameters")? as usize,
        control_file: o.get("_drs_control_file").unwrap_or("").to_string(),
        name: o.get("_pd_block_id").and_then(|v| v.strip_suffix("|overall")).unwrap_or("").to_string(),
        ..Default::default()
    };
    let (names, values, sus) = (col(o, "_drs_parameter_name")?, col(o, "_drs_parameter_value")?, col(o, "_drs_parameter_su")?);
    r.refined = names.iter().zip(&values).zip(&sus).map(|((a, v), s)| (a.clone(), n(v), n(s))).collect();
    let (ids, pblocks) = (col(o, "_pd_phase_id")?, col(o, "_pd_phase_block_id")?);
    r.phase_blocks = ids.iter().zip(&pblocks).map(|(i, b)| (i.parse().unwrap_or(0), b.clone())).collect();
    let (cp, ck) = (col(o, "_drs_contribution_phase_id")?, col(o, "_drs_contribution_pattern_id")?);
    let (kid, kblock, kfile, kw, kref) = (col(o, "_drs_pattern_id")?, col(o, "_drs_pattern_block_id")?, col(o, "_drs_pattern_data_file")?, col(o, "_drs_pattern_weight")?, col(o, "_drs_pattern_refined")?);
    for i in 0..kid.len() {
        let id: usize = kid[i].parse().map_err(|_| format!("_drs_pattern_id {}", kid[i]))?;
        let b = blocks.iter().find(|b| b.get("_pd_block_id") == Some(kblock[i].as_str())).ok_or(format!("no block {}", kblock[i]))?;
        let label = std::path::Path::new(&kfile[i]).file_stem().map_or(kfile[i].clone(), |s| s.to_string_lossy().to_string());
        let reflections = match (b.column("_drs_refln_pattern_id"), b.column("_pd_refln_phase_id")) {
            (Some(k), Some(ph)) => ph.iter().zip(k).filter(|(p, _)| **p != ".").map(|(p, k)| (p.parse().unwrap_or(0), k.parse().unwrap_or(0))).collect(),
            _ => Vec::new(),
        };
        r.patterns.push(PatternResults {
            id,
            label,
            block: kblock[i].clone(),
            weight: n(&kw[i]),
            refined: kref[i] == "yes",
            rp: 100.0 * f(b, "_pd_proc_ls_prof_r_factor")?,
            rwp: 100.0 * f(b, "_pd_proc_ls_prof_wr_factor")?,
            rexp: 100.0 * f(b, "_pd_proc_ls_prof_wr_expected")?,
            chi2: f(b, "_drs_chi2")?,
            chi2_bragg: f(b, "_drs_chi2_bragg")?,
            phases: cp.iter().zip(&ck).filter(|(_, k)| k.parse() == Ok(id)).map(|(p, _)| p.parse().unwrap_or(0)).collect(),
            points: b.column("_pd_proc_point_id").map_or(0, |c| c.len()),
            reflections,
            axis: if b.column("_pd_meas_time_of_flight").is_some() { crate::formats::Axis::TimeOfFlight } else { crate::formats::Axis::TwoTheta },
            two_theta_fixed: b.get("_pd_meas_2theta_fixed").and_then(crate::cif::number),
            control_file: b.get("_drs_control_file").unwrap_or("").to_string(),
        });
    }
    Ok(r)
}

/// A joint parameter by its FullProf number and its members' keys, (pattern, the first field its free parameter
/// moves).
type Keyed = (Option<u32>, Vec<(usize, Field)>);

/// A joint parameter by its FullProf number and its members, (pattern, index into its free parameters).
type Resolved = (Option<u32>, Vec<(usize, usize)>);

/// The format line of a joint control file.
pub const JOINT_FORMAT: &str = "diff-ref-struct joint 1";

/// Whether a control file's text is a joint refinement's (its format line, or a `[[joint.pattern]]` table).
pub fn is_joint_control(text: &str) -> bool {
    text.lines().any(|l| {
        let t = l.trim();
        t.starts_with("[[joint.pattern]]") || (t.starts_with("format") && t.contains(JOINT_FORMAT))
    })
}

/// A refinement's description with every field mapped (`map`), the terms whose field maps to None left out, the
/// parameters left empty dropped.
fn remap_refinement(rf: &mut crate::refinable::Refinement, map: &dyn Fn(Field) -> Option<Field>) {
    for p in rf.parameters.iter_mut() {
        *p = p.iter().filter_map(|&(f, c)| map(f).map(|g| (g, c))).collect();
    }
    rf.parameters.retain(|p| !p.is_empty());
    rf.limits = rf.limits.iter().filter_map(|&(f, a, b)| map(f).map(|g| (g, a, b))).collect();
    rf.behaviour = rf.behaviour.iter().filter_map(|&(f, b)| map(f).map(|g| (g, b))).collect();
    rf.links = rf.links.iter().filter_map(|l| Some(crate::steering::Link { value: map(l.value)?, follows: map(l.follows)?, ..*l })).collect();
    for v in rf.variables.iter_mut() {
        v.terms = v.terms.iter().filter_map(|&(f, c)| map(f).map(|g| (g, c))).collect();
    }
    rf.variables.retain(|v| !v.terms.is_empty());
}

impl Joint {
    /// Each joint parameter's members as (pattern, the first field its free parameter moves): what identifies them
    /// when a pattern's parameters are built again.
    fn member_keys(&self) -> Vec<Vec<(usize, Field)>> {
        self.parameters.iter().map(|q| q.members.iter().map(|&(k, i)| (k, self.patterns[k].set.free[i].fields[0].0)).collect()).collect()
    }

    /// The joint parameters from their members' keys (`member_keys`), each pattern's free parameters found by their
    /// first field; a free parameter no key names is a parameter of its own, in its pattern's order.
    fn tie_keys(&mut self, keys: Vec<Keyed>) {
        let find = |p: &Member, f: Field| p.set.free.iter().position(|x| x.fields[0].0 == f);
        let mut out: Vec<Parameter> = Vec::new();
        let mut taken: std::collections::BTreeSet<(usize, usize)> = Default::default();
        let resolved: Vec<Resolved> = keys
            .into_iter()
            .map(|(n, ms)| (n, ms.into_iter().filter_map(|(k, f)| self.patterns.get(k).and_then(|p| find(p, f)).map(|i| (k, i))).collect::<Vec<_>>()))
            .filter(|(_, ms)| !ms.is_empty())
            .collect();
        let owner: std::collections::BTreeMap<(usize, usize), usize> = resolved.iter().enumerate().flat_map(|(j, (_, ms))| ms.iter().map(move |&m| (m, j))).collect();
        let mut placed = vec![false; resolved.len()];
        for (k, p) in self.patterns.iter().enumerate() {
            for (i, f) in p.set.free.iter().enumerate() {
                if taken.contains(&(k, i)) {
                    continue;
                }
                match owner.get(&(k, i)) {
                    Some(&j) if !placed[j] => {
                        placed[j] = true;
                        taken.extend(resolved[j].1.iter().copied());
                        out.push(Parameter { name: f.name.clone(), number: resolved[j].0, members: resolved[j].1.clone() });
                    }
                    Some(_) => {}
                    None => {
                        taken.insert((k, i));
                        out.push(Parameter { name: f.name.clone(), number: None, members: vec![(k, i)] });
                    }
                }
            }
        }
        for q in out.iter_mut() {
            let ks: std::collections::BTreeSet<usize> = q.members.iter().map(|m| m.0).collect();
            let (k, _) = q.members[0];
            q.name = joint_name(&q.name, k, &self.patterns[k].phases, self.patterns.len() > 1 && ks.len() == 1);
        }
        self.parameters = out;
    }

    /// Set pattern k's weight against the others' as they are, then every weight normalised to sum 1 (as W_PAT is
    /// read; app request 43). Two patterns of 0.5 and pattern 1 set to 1 give 2/3 and 1/3.
    pub fn set_weight(&mut self, k: usize, w: f64) -> Result<(), String> {
        let mut ws: Vec<f64> = self.patterns.iter().map(|p| p.weight).collect();
        *ws.get_mut(k).ok_or(format!("the joint refinement has {} patterns", self.patterns.len()))? = w;
        self.set_weights(&ws)
    }

    /// Every pattern's weight (any positive numbers), normalised to sum 1.
    pub fn set_weights(&mut self, w: &[f64]) -> Result<(), String> {
        if w.len() != self.patterns.len() {
            return Err(format!("{} weights for {} patterns", w.len(), self.patterns.len()));
        }
        if let Some(x) = w.iter().find(|x| !(x.is_finite() && **x > 0.0)) {
            return Err(format!("{x}: a pattern's weight is a number above 0"));
        }
        let total: f64 = w.iter().sum();
        for (p, x) in self.patterns.iter_mut().zip(w) {
            p.weight = x / total;
        }
        Ok(())
    }

    /// Refine pattern k, or calculate it only (FullProf's flag 0); at least one pattern stays refined.
    pub fn set_refined(&mut self, k: usize, on: bool) -> Result<(), String> {
        let n = self.patterns.len();
        if k >= n {
            return Err(format!("the joint refinement has {n} patterns"));
        }
        if !on && self.patterns.iter().enumerate().all(|(j, p)| j == k || !p.refined) {
            return Err("at least one pattern is refined".into());
        }
        self.patterns[k].refined = on;
        Ok(())
    }

    /// Phase g (0, 1, …) contributing to pattern k or not (app request 43). Switched off, the pattern's model is
    /// relieved of it, with its free parameters there. Switched on, the pattern's model is given the phase as the
    /// first pattern it contributes to has it (its structure, scale, profile and cell, the atoms' scattering taken
    /// again for this pattern's radiation), with that pattern's free parameters of the phase: those of its structure
    /// (atoms, modes) tied to that pattern's, the others (scale, profile, cell) this pattern's own. The values kept,
    /// the joint parameters tied again.
    pub fn set_contribution(&mut self, k: usize, g: usize, on: bool) -> Result<(), String> {
        if k >= self.patterns.len() || g >= self.phases.len() {
            return Err(format!("the joint refinement has {} patterns and {} phases", self.patterns.len(), self.phases.len()));
        }
        let has = self.patterns[k].phases.contains(&g);
        if has == on {
            return Ok(());
        }
        let mut keys: Vec<Keyed> = self.parameters.iter().map(|q| q.number).zip(self.member_keys()).collect();
        if !on {
            let p = &mut self.patterns[k];
            if p.phases.len() == 1 {
                return Err(format!("phase {} is pattern {}'s only phase", g + 1, k + 1));
            }
            let q = p.phases.iter().position(|&x| x == g).unwrap();
            p.model.remove_phase(q)?;
            p.phases.remove(q);
            let map = move |f: Field| match f.phase() {
                Some(r) if r == q => None,
                Some(r) if r > q => Some(f.with_phase(r - 1)),
                _ => Some(f),
            };
            remap_refinement(&mut p.refinement, &map);
            for (_, ms) in keys.iter_mut() {
                *ms = ms.iter().filter_map(|&(kk, f)| if kk == k { map(f).map(|g| (kk, g)) } else { Some((kk, f)) }).collect();
            }
        } else {
            let s = self.patterns.iter().position(|p| p.phases.contains(&g)).ok_or(format!("phase {} contributes to no pattern to take it from", g + 1))?;
            let qs = self.patterns[s].phases.iter().position(|&x| x == g).unwrap();
            let mut phase = self.patterns[s].model.phases[qs].clone();
            let source_rf: Vec<Vec<(Field, f64)>> = self.patterns[s].refinement.parameters.iter().filter(|t| t.iter().all(|(f, _)| f.phase() == Some(qs))).cloned().collect();
            let p = &mut self.patterns[k];
            if let crate::powder::Intensities::Structure(st, sc) = &mut phase.intensities {
                let (radiation, user) = (p.model.radiation(), p.model.user_scatterers());
                *sc = st.atoms.iter().map(|a| crate::scattering::resolve(&a.species, radiation, &user)).collect::<Result<Vec<_>, _>>()?;
            }
            p.model.phases.push(phase);
            let q = p.model.phases.len() - 1;
            p.phases.push(g);
            for t in &source_rf {
                p.refinement.parameters.push(t.iter().map(|&(f, c)| (f.with_phase(q), c)).collect());
                let (first_s, first_k) = (t[0].0, t[0].0.with_phase(q));
                match keys.iter_mut().find(|(_, ms)| first_s.is_structural() && ms.contains(&(s, first_s))) {
                    Some((_, ms)) => ms.push((k, first_k)),
                    None => keys.push((None, vec![(k, first_k)])),
                }
            }
        }
        let p = &mut self.patterns[k];
        p.set = p.refinement.build(&p.model)?;
        p.numbers = vec![None; p.set.free.len()];
        self.tie_keys(keys);
        Ok(())
    }

    /// Pattern k's control file: its data, model and free parameters, the files it names (the data, the files the
    /// operations read, the background points) relative to `folder` (app 45).
    pub fn pattern_control(&self, k: usize, folder: &Path) -> crate::control::Control {
        let p = &self.patterns[k];
        let rel = |f: &str| crate::control::relative(&self.dir.join(f), folder);
        let source = p.source.clone().map(|mut d| {
            d.file = rel(&d.file);
            for op in &mut d.operations {
                if let Some(f) = op.file_mut() {
                    *f = rel(f);
                }
            }
            d
        });
        crate::control::Control { title: p.title.clone(), model: p.model.clone(), data: source, refinement: p.refinement.clone(), fullprof: p.fullprof.clone(), background_file: p.background_file.as_deref().map(rel) }
    }

    /// The control files of the joint refinement (P06, app request 43), each (name, text): first the joint one,
    /// `{stem}.drc`, then each pattern's, `{stem} pattern {k}.drc` (the file dictionary's `pattern_control`), their
    /// files named relative to `folder`, where they are to be written. The joint file:
    ///
    /// ```text
    /// format = "diff-ref-struct joint 1"
    /// title = "…"
    /// [[joint.pattern]]
    /// control = "pbb pattern 1.drc"   # the pattern's control file: its data, model and free parameters
    /// weight = 0.5                    # W_k, normalised to sum 1 when read
    /// refined = true
    /// phases = [1]                    # the joint phases its model holds, in its order
    /// [[joint.phase]]
    /// name = "PbSO4"
    /// [[joint.tie]]                   # a parameter several patterns share: each member "k: the address of the first
    /// name = "phase1.cell_a"          # value its free parameter moves" in pattern k's control file
    /// members = ["1: phase1.cell_a", "2: phase1.cell_a"]
    /// ```
    pub fn control_files(&self, stem: &str, folder: &Path) -> Result<Vec<(String, String)>, String> {
        use crate::control::syntax::{float, string};
        let mut out = Vec::new();
        let mut j = format!("# diff-ref-struct control file of a joint refinement (P06; https://jmigartua.github.io/diff-ref-struct-site/)\nformat = {}\ntitle = {}\n", string(JOINT_FORMAT), string(&self.title));
        let mut patterns = Vec::new();
        for (k, p) in self.patterns.iter().enumerate() {
            let name = crate::formats::files::path("pattern_control", &[("control", stem), ("n", &(k + 1).to_string())])?;
            patterns.push((name.clone(), self.pattern_control(k, folder).to_text()));
            j.push_str(&format!(
                "\n[[joint.pattern]]\ncontrol = {}\nweight = {}\nrefined = {}\nphases = [{}]\n",
                string(&name),
                float(p.weight),
                p.refined,
                p.phases.iter().map(|g| (g + 1).to_string()).collect::<Vec<_>>().join(", ")
            ));
        }
        for n in &self.phases {
            j.push_str(&format!("\n[[joint.phase]]\nname = {}\n", string(n.trim())));
        }
        for q in self.parameters.iter().filter(|q| q.members.iter().map(|m| m.0).collect::<std::collections::BTreeSet<_>>().len() > 1) {
            let members: Vec<String> = q.members.iter().map(|&(k, i)| string(&format!("{}: {}", k + 1, crate::refinable::address(self.patterns[k].set.free[i].fields[0].0, &self.patterns[k].model)))).collect();
            j.push_str(&format!("\n[[joint.tie]]\nname = {}\nmembers = [{}]\n", string(&q.name), members.join(", ")));
        }
        out.push((format!("{stem}.drc"), j));
        out.extend(patterns);
        Ok(out)
    }

    /// Write the control files (`control_files`) to `path`'s folder, the joint file at `path`; the paths written.
    pub fn save(&self, path: &Path) -> Result<Vec<std::path::PathBuf>, String> {
        let folder = path.parent().map(|d| d.to_path_buf()).unwrap_or_default();
        let stem = path.file_stem().map(|s| s.to_string_lossy().to_string()).ok_or(format!("{}: no file name", path.display()))?;
        let mut written = Vec::new();
        for (name, text) in self.control_files(&stem, &folder)? {
            let full = folder.join(&name);
            std::fs::write(&full, text).map_err(|e| format!("{}: {e}", full.display()))?;
            written.push(full);
        }
        Ok(written)
    }

    /// Open a joint control file (`control_files`): each pattern's control file and data, its free parameters from
    /// its refinement, the patterns tied as its `[[joint.tie]]` tables say.
    pub fn open(path: &Path) -> Result<Joint, String> {
        use crate::control::syntax::Value;
        let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
        let doc = crate::control::syntax::parse(&text);
        let at = |sp: crate::control::syntax::Span| format!("{}:{}:{}", path.display(), sp.line, sp.column);
        if let Some(d) = doc.diagnostics.iter().find(|d| !d.warning) {
            return Err(format!("{}: {}", at(d.span), d.message));
        }
        let dir = path.parent().map(|d| d.to_path_buf()).unwrap_or_default();
        let root = &doc.tables[0];
        let string_of = |t: &crate::control::syntax::Table, key: &str| -> Result<String, String> {
            match t.get(key).map(|e| &e.value.value) {
                Some(Value::String(s)) => Ok(s.clone()),
                _ => Err(format!("{}: {} needs {key} = \"…\"", t.header.map_or(path.display().to_string(), at), t.name())),
            }
        };
        let format = string_of(root, "format")?;
        if format != JOINT_FORMAT {
            return Err(format!("{}: format \"{format}\" is not a joint control file's (\"{JOINT_FORMAT}\")", path.display()));
        }
        let title = string_of(root, "title").unwrap_or_default();
        let tables = |name: &'static str| doc.tables.iter().filter(move |t| t.array && t.path.join(".") == name);
        let phases: Vec<String> = tables("joint.phase").map(|t| string_of(t, "name")).collect::<Result<_, _>>()?;
        let mut patterns = Vec::new();
        for t in tables("joint.pattern") {
            let where_ = t.header.map_or(path.display().to_string(), at);
            let control = string_of(t, "control")?;
            let o = crate::control::open(&dir.join(&control))?;
            let num = |key: &str| match t.get(key).map(|e| &e.value.value) {
                Some(Value::Float(v)) => Some(*v),
                Some(Value::Integer(v)) => Some(*v as f64),
                _ => None,
            };
            let weight = num("weight").ok_or(format!("{where_}: weight = <a number above 0>"))?;
            let refined = match t.get("refined").map(|e| &e.value.value) {
                Some(Value::Bool(b)) => *b,
                None => true,
                _ => return Err(format!("{where_}: refined = true or false")),
            };
            let ph: Vec<usize> = match t.get("phases").map(|e| &e.value.value) {
                Some(Value::Array(a)) => a.iter().map(|i| if let Value::Integer(v) = i.value { (v >= 1 && (v as usize) <= phases.len()).then(|| v as usize - 1).ok_or(format!("{}: phase {v} is not one of the {} [[joint.phase]]", at(i.span), phases.len())) } else { Err(format!("{}: a phase's number", at(i.span))) }).collect::<Result<_, _>>()?,
                _ => (0..o.control.model.phases.len()).collect(),
            };
            if ph.len() != o.control.model.phases.len() {
                return Err(format!("{where_}: {} phases named for the {} of {control}", ph.len(), o.control.model.phases.len()));
            }
            let data = o.data.ok_or(format!("{control}: no data"))?;
            let c = o.control;
            let set = c.refinement.build(&c.model).map_err(|e| format!("{control}: {e}"))?;
            let file = c.data.as_ref().map_or(String::new(), |d| d.file.clone());
            let n = set.free.len();
            patterns.push(Member { label: label_of(&file), file, model: c.model, data, weight, refined, phases: ph, set, numbers: vec![None; n], title: c.title, refinement: c.refinement, source: c.data, fullprof: c.fullprof, background_file: c.background_file });
        }
        if patterns.is_empty() {
            return Err(format!("{}: no [[joint.pattern]]", path.display()));
        }
        let (cycles, eps) = (patterns[0].refinement.cycles, patterns[0].refinement.eps);
        let mut j = Joint { title, phases, patterns, parameters: Vec::new(), cycles, eps, dir };
        let ws: Vec<f64> = j.patterns.iter().map(|p| p.weight).collect();
        j.set_weights(&ws).map_err(|e| format!("{}: {e}", path.display()))?;
        let mut keys = Vec::new();
        for t in tables("joint.tie") {
            let Some(Value::Array(a)) = t.get("members").map(|e| &e.value.value) else { return Err(format!("{}: members = [\"k: address\", …]", t.header.map_or(path.display().to_string(), at))) };
            let mut ms = Vec::new();
            for i in a {
                let Value::String(m) = &i.value else { return Err(format!("{}: a member is \"k: address\"", at(i.span))) };
                let (k, addr) = m.split_once(':').ok_or(format!("{}: {m}: \"k: address\"", at(i.span)))?;
                let k: usize = k.trim().parse::<usize>().ok().filter(|&k| k >= 1 && k <= j.patterns.len()).ok_or(format!("{}: {m}: no pattern {}", at(i.span), k.trim()))? - 1;
                let f = crate::refinable::field_of(addr.trim(), &j.patterns[k].model).map_err(|e| format!("{}: {e}", at(i.span)))?;
                if !j.patterns[k].set.free.iter().any(|x| x.fields[0].0 == f) {
                    return Err(format!("{}: {m}: pattern {}'s control file does not refine it", at(i.span), k + 1));
                }
                ms.push((k, f));
            }
            keys.push((None, ms));
        }
        j.tie_keys(keys);
        Ok(j)
    }
}
