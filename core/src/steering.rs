//! Steering a refinement (W12, L107–L110): how each free parameter is allowed to move, written in the control file
//! (`[refine]`, `[[refine.link]]`, `[[refine.variable]]`, `[[refine.step]]`) and read by the engine (`lsq::Session`).
//!
//! - **Damping** (`Steering::damping`, `Steering::class_damping`): the shift the least squares compute for a parameter
//!   is multiplied by its damping before it is applied: d = d_parameter × d_class(first field). The esd is never
//!   changed: it is that of the parameter's value, whatever step reached it; convergence is judged on the computed
//!   (undamped) shift. FullProf's codeword multiplier and relaxation factors both damp its shifts (manual p. 13, p. 43,
//!   LINE 9 on p. 77); its classes are `Class`. FullProf also multiplies the printed esd by the multiplier (its eq. 3.7,
//!   p. 43; measured in `manuals/fullprof/runs/shift_control/`); the core does not (`lsq::fullprof_esd`).
//! - **Maximum shift** (`Steering::max_shift`, the dictionary's `max_shift`): per cycle, each field moves at most this
//!   much, in its unit; a parameter whose step would move a field further has its whole step shortened, and the cycle
//!   says so.
//! - **λ** (`Steering::lambda`): the Levenberg–Marquardt factor the first cycle starts from (0: a Gauss–Newton step
//!   first); each cycle reports the λ it used. When a step raises χ², λ is raised and the step shortened, and the
//!   cycle's notes record each try. `marquardt = false` applies every computed step as FullProf does, damped only by
//!   the damping: then χ² may rise, and the divergence guards stop the run.
//! - **Bounds** (`Behaviour`): at a bound a parameter stops there (held, as P39's limits), is reflected back inside,
//!   or, between two bounds, wraps round (periodic).
//! - **Links** (`Link`): a value follows another as x = c·p + d (two species sharing a site as p and s − p, equal B's,
//!   one value across phases).
//! - **Variables** (`Variable`): a new parameter v refined in place of the fields it moves, f_i = f_i⁰ + c_i v; the
//!   amplitude of a symmetry mode is one, the c_i being the mode's displacements per unit amplitude.
//! - **Divergence** (`Guards`): χ² rising for k cycles, a shift of many esd's, an impossible value, a singular matrix.
//! - **A plan** (`PlanStep`): steps, each freeing part of the parameters with its damping, cycles and stop rule.

use crate::jacobian::Field;
use crate::powder::Model;

/// FullProf's four classes of parameters, each with its relaxation factor (manual, LINE 9, p. 77): the class of a
/// parameter is the dictionary's `class` (`class_of`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Class {
    /// Coordinates, magnetic moments, occupancies, isotropic B's (R_at).
    Atomic,
    /// Anisotropic displacement parameters (R_an).
    Anisotropic,
    /// Profile, asymmetry, overall B, cell, preferred orientation, strain, size (R_pr).
    Profile,
    /// Zero, background, displacement, transparency (R_gl).
    Global,
}

impl Class {
    pub const ALL: [Class; 4] = [Class::Atomic, Class::Anisotropic, Class::Profile, Class::Global];

    /// The word the control file and the dictionary use.
    pub fn word(self) -> &'static str {
        match self {
            Class::Atomic => "atomic",
            Class::Anisotropic => "anisotropic",
            Class::Profile => "profile",
            Class::Global => "global",
        }
    }

    pub fn from_word(w: &str) -> Option<Class> {
        Class::ALL.into_iter().find(|c| c.word() == w)
    }

    pub fn index(self) -> usize {
        self as usize
    }
}

/// The class of a field (the dictionary's `class`, checked in `parameter_dictionary`). Where FullProf's manual names
/// no class (scale, ratio, mode amplitudes) the dictionary says what is believed.
pub fn class_of(f: Field) -> Class {
    match f {
        Field::AtomX(..) | Field::AtomB(..) | Field::AtomOcc(..) | Field::Mode(..) | Field::Rigid(..) => Class::Atomic,
        Field::AtomBeta(..) => Class::Anisotropic,
        Field::Zero | Field::Sycos | Field::Sysin | Field::Ratio | Field::Wavelength | Field::Micro(_) | Field::Background(_) | Field::BackgroundTerm(_) => Class::Global,
        // A TOF bank's TOF(d) is global, as Zero (manual LINE 9, p. 77: R_gl); a phase's TOF profile is its profile.
        Field::Tof(_) => Class::Global,
        Field::TofProfile(..) => Class::Profile,
        Field::Scale(_) | Field::Bov(_) | Field::Cell(..) | Field::U(_) | Field::V(_) | Field::W(_) | Field::X(_) | Field::Y(_) | Field::GaussSize(_) | Field::Eta0(_) | Field::Asym(..) | Field::SL(_) | Field::DL(_) | Field::Pref(..) | Field::Texture(..) => Class::Profile,
        // Not named by the manual's LINE 9 (p. 77); believed profile, as LINES 26–38 are the phase's profile (L211).
        Field::Quartic(..) | Field::QuarticShare(_) | Field::ShapeSize(_) => Class::Profile,
        // A Pawley fit's intensities (P36): FullProf, which has no Pawley fit, names no class; profile, as the
        // dictionary has it for a list's intensities.
        Field::Intensity(..) => Class::Profile,
    }
}

/// What a parameter does at a bound.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Behaviour {
    /// It stops at the bound and is held there while the step pushes it out (P39's limits).
    #[default]
    Stop,
    /// A step that would cross the bound is reflected back inside by as much as it would have crossed.
    Reflect,
    /// Between its two bounds, the value wraps round: past the upper bound it comes back from the lower one (a
    /// coordinate taken modulo 1, an angle modulo 360°).
    Periodic,
}

impl Behaviour {
    pub fn word(self) -> &'static str {
        match self {
            Behaviour::Stop => "stop",
            Behaviour::Reflect => "reflect",
            Behaviour::Periodic => "periodic",
        }
    }

    pub fn from_word(w: &str) -> Option<Behaviour> {
        [Behaviour::Stop, Behaviour::Reflect, Behaviour::Periodic].into_iter().find(|b| b.word() == w)
    }
}

/// A value that follows another: value = coefficient × follows + offset. If `follows` is refined, `value` moves with
/// it (coefficient × its shift); if not, the value is still set from it at the start of a run.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Link {
    pub value: Field,
    pub follows: Field,
    pub coefficient: f64,
    pub offset: f64,
}

/// A new variable, refined in place of the fields it moves: each field f_i = f_i⁰ + c_i·v, with `terms` the (f_i, c_i)
/// and `value` the variable's value at the model's present values.
#[derive(Clone, Debug, PartialEq)]
pub struct Variable {
    pub name: String,
    pub value: f64,
    pub terms: Vec<(Field, f64)>,
}

/// The divergence guards (W11), each with its threshold; a threshold of 0 turns its guard off.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Guards {
    /// Stop when χ² has risen this many cycles in a row (default 3).
    pub rising: usize,
    /// Stop when a parameter's computed shift is more than this many esd's (default 100) …
    pub jump: f64,
    /// … from this cycle on (default 3: the first cycles from rough values legitimately move many esd's).
    pub jump_from: usize,
    /// Stop when the normal matrix is singular (a combination of parameters the data cannot fix, named as `check`
    /// names it, R24) rather than leave the combination out and go on (default true).
    pub singular: bool,
}

impl Default for Guards {
    fn default() -> Guards {
        Guards { rising: 3, jump: 100.0, jump_from: 3, singular: true }
    }
}

/// How a refinement is steered: what the control file's `[refine]` table says beyond what is free.
#[derive(Clone, Debug, PartialEq)]
pub struct Steering {
    /// λ's start (the Marquardt factor of the first cycle): 0 for a Gauss–Newton step first.
    pub lambda: f64,
    /// Raise λ and shorten the step when χ² rises (true), or apply every step (FullProf's way, false).
    pub marquardt: bool,
    /// Damping of each class (`Class::index`): FullProf's R_at, R_an, R_pr, R_gl.
    pub class_damping: [f64; 4],
    /// Damping of a free parameter, by its name (the address of its first value, or a variable's name).
    pub damping: Vec<(String, f64)>,
    /// The user's maximum shift of a value per cycle, in its unit (over the dictionary's default).
    pub max_shift: Vec<(Field, f64)>,
    pub guards: Guards,
}

impl Default for Steering {
    fn default() -> Steering {
        Steering { lambda: 0.0, marquardt: true, class_damping: [1.0; 4], damping: Vec::new(), max_shift: Vec::new(), guards: Guards::default() }
    }
}

impl Steering {
    /// The damping of a free parameter named `name` whose first value is `first`: its own times its class's.
    pub fn damping_of(&self, name: &str, first: Field) -> f64 {
        let own = self.damping.iter().find(|d| d.0 == name).map_or(1.0, |d| d.1);
        own * self.class_damping[class_of(first).index()]
    }

    /// The maximum shift of a value per cycle: the user's, else the dictionary's (`max_shift_default`); None for no
    /// limit.
    pub fn max_shift_of(&self, f: Field, m: &Model) -> Option<f64> {
        self.max_shift.iter().find(|s| s.0 == f).map(|s| s.1).or_else(|| max_shift_default(f, m))
    }
}

/// The dictionary's maximum shift per cycle of a field (`max_shift` in `parameters.yml`), None where it gives none.
pub fn max_shift_default(f: Field, m: &Model) -> Option<f64> {
    let id = crate::refinable::dictionary_id(f, m);
    crate::parameters::lookup(id).and_then(|p| p.max_shift.parse::<f64>().ok())
}

/// One step of a refinement plan: which parameters are free (addresses, prefixes of addresses ending in `*`, variable
/// names, class words, or "all"), with the step's damping (over the steering's), its cycles and its convergence.
#[derive(Clone, Debug, PartialEq)]
pub struct PlanStep {
    pub name: String,
    pub free: Vec<String>,
    pub damping: Option<f64>,
    pub cycles: Option<usize>,
    pub eps: Option<f64>,
}

impl PlanStep {
    /// Whether the step frees a parameter named `name` whose first value is `first`.
    pub fn frees(&self, name: &str, first: Field) -> bool {
        self.free.iter().any(|s| selects(s, name, first))
    }
}

/// Whether a selector of a plan step names a parameter: "all", a class word, the parameter's name, a prefix of it
/// ending in `*`, or a scope of it (`phase1.O1` selects `phase1.O1.atom_position.x`).
pub fn selects(sel: &str, name: &str, first: Field) -> bool {
    if sel == "all" {
        return true;
    }
    if let Some(c) = Class::from_word(sel) {
        return class_of(first) == c;
    }
    if let Some(prefix) = sel.strip_suffix('*') {
        return name.starts_with(prefix);
    }
    name == sel || name.strip_prefix(sel).is_some_and(|rest| rest.starts_with('.'))
}
