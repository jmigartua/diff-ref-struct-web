//! Named parameters, constraints and limits (P37–P39, L053).
//!
//! Every refinable value of the model is a `jacobian::Field`, addressed by its scope and its parameter-dictionary id
//! (`ecosystem/data/parameters.yml`): `pattern1.zero_shift`, `phase1.cell_a`, `phase2.O1.atom_position.x`,
//! `phase1.atom_betas.11` … (`address`).
//!
//! A refinement moves a set of free parameters q_j; each field moves with them linearly, f = f₀ + Σ_j M_fj δq_j
//! (`Free::fields` holds the column of M). The set is built from a `.pcr` (`ParameterSet::from_pcr`):
//! - codewords (FullProf manual p. 13): C = sign(a)(10p + |a|) refines a field with parameter p and multiplier a;
//!   fields sharing p are one free parameter, each moving by its multiplier (equality constraints, P38); a bare
//!   "1.00" (|C| < 10, Aut = 1) is a parameter of its own; `VARY ratio` refines the Kα2/Kα1 ratio;
//! - the multipliers are carried into the default model the import translated the values into (L050): the X-ray scale
//!   ×2, P1 + 3P2 and P3 + 3P4, the axial pair, −G1 for the exponential preferred orientation;
//! - the automatic constraints of symmetry (P38) then replace the codewords' where they apply: the cell parameters
//!   tied or fixed by the crystal system, and each atom's coordinates and betas by its site symmetry (`site`), as
//!   FullProf does with Aut = 1;
//! - limits (P39, `limits`): the bounds of the dictionary, which the engine respects.

use crate::jacobian::Field;
use crate::pcr::Pcr;
use crate::powder::{BackgroundModel, InstrumentalWidths, Intensities, Model, Shape};
use crate::site::Site;
use crate::steering::{Behaviour, Link, PlanStep, Steering, Variable};
use std::collections::BTreeMap;

/// A free parameter: its name (the address of its first field, or a variable's name) and the fields it moves, with
/// their multipliers; the damping of its shifts (its own times its class's, `Steering::damping_of`); and, for a
/// variable (`steering::Variable`), its value at the model's present values.
#[derive(Clone, Debug, PartialEq)]
pub struct Free {
    pub name: String,
    pub fields: Vec<(Field, f64)>,
    pub damping: f64,
    pub value: Option<f64>,
}

impl Free {
    /// Its value in a model: a variable's own, else its first field's.
    pub fn value_in(&self, m: &Model) -> f64 {
        self.value.unwrap_or_else(|| m.get(self.fields[0].0))
    }

    /// The value's shift for a shift δ of the parameter.
    pub fn shift_of(&self, delta: f64) -> f64 {
        if self.value.is_some() {
            delta
        } else {
            delta * self.fields[0].1
        }
    }

    /// The value's esd for the parameter's esd σ.
    pub fn esd_of(&self, sigma: f64) -> f64 {
        if self.value.is_some() {
            sigma
        } else {
            sigma * self.fields[0].1.abs()
        }
    }
}

/// A linear restraint on free parameters (P40): Σ c·p ≈ target ± sigma, p the value of the named parameter (that of
/// its first field over its multiplier).
#[derive(Clone, Debug, PartialEq)]
pub struct LinearRestraint {
    pub name: String,
    pub terms: Vec<(String, f64)>,
    pub target: f64,
    pub sigma: f64,
}

impl LinearRestraint {
    /// Its value in a model, and its derivative with respect to each free parameter of `set`.
    pub fn value(&self, m: &Model, set: &ParameterSet) -> (f64, Vec<f64>) {
        let mut d = vec![0.0; set.free.len()];
        let mut v = 0.0;
        for (name, c) in &self.terms {
            if let Some(j) = set.free.iter().position(|f| &f.name == name) {
                let p = &set.free[j];
                let (f, mult) = p.fields[0];
                v += c * p.value.unwrap_or_else(|| m.get(f) / mult);
                d[j] += c;
            }
        }
        (v, d)
    }
}

/// The free parameters of a refinement, the linear restraints on them, the user's limits and what each does at its
/// bound, the links whose values follow others, and the steering (W12).
#[derive(Clone, Debug, PartialEq)]
pub struct ParameterSet {
    pub free: Vec<Free>,
    pub linear: Vec<LinearRestraint>,
    /// Limits given by the user, which replace `limits` for their field.
    pub user_limits: Vec<(Field, Option<f64>, Option<f64>)>,
    /// What a field does at its bounds, where it is not `Behaviour::Stop`.
    pub behaviour: Vec<(Field, Behaviour)>,
    /// Values that follow others (x = c·p + d), set at the start of a run (`enforce_links`) and moved with p.
    pub links: Vec<Link>,
    pub steering: Steering,
}

/// What is refined, as the control file writes it (L066): each free parameter as the values it moves before the
/// symmetry's ties, with their multipliers (value = multiplier × parameter; the first value names the parameter, a
/// value may move with several parameters and then the terms add), the linear restraints, the user's limits, and the
/// number of cycles and the convergence. `build` applies the crystal system's and the sites' constraints, as
/// `ParameterSet::from_pcr` does, so that the file need not repeat what the symmetry implies.
///
/// W12 adds the steering (λ, damping, maximum shifts, guards), what each bounded value does at its bound, the links,
/// the new variables and a plan of steps (`steering`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Refinement {
    pub parameters: Vec<Vec<(Field, f64)>>,
    pub linear: Vec<LinearRestraint>,
    pub limits: Vec<(Field, Option<f64>, Option<f64>)>,
    pub cycles: Option<usize>,
    pub eps: Option<f64>,
    pub steering: Steering,
    pub behaviour: Vec<(Field, Behaviour)>,
    pub links: Vec<Link>,
    pub variables: Vec<Variable>,
    pub plan: Vec<PlanStep>,
}

impl Refinement {
    /// The parameter set: one column per parameter, the fields in the order they are written (each parameter's,
    /// then the next's), the symmetry applied.
    pub fn build(&self, m: &Model) -> Result<ParameterSet, String> {
        let mut rows: Rows = BTreeMap::new();
        let mut order = Vec::new();
        let mut namer = Vec::new();
        for (col, terms) in self.parameters.iter().enumerate() {
            let Some(&(first, _)) = terms.first() else { return Err(format!("parameter {} moves no value", col + 1)) };
            namer.push(first);
            for &(f, c) in terms {
                if rows.entry(f).or_default().insert(col, c).is_some() {
                    return Err(format!("'{}' is written twice in the parameter '{}'", address(f, m), address(first, m)));
                }
                order.push(f);
            }
        }
        // The new variables: columns of their own, after the parameters', through the symmetry like them.
        let first_var = self.parameters.len();
        for (k, v) in self.variables.iter().enumerate() {
            let Some(&(first, _)) = v.terms.first() else { return Err(format!("the variable '{}' moves no value", v.name)) };
            if v.name.is_empty() || v.name.contains(char::is_whitespace) || m.all_fields().iter().any(|&f| address(f, m) == v.name) {
                return Err(format!("'{}' cannot name a variable: a variable's name is one word and not a value's address", v.name));
            }
            if self.variables[..k].iter().any(|w| w.name == v.name) {
                return Err(format!("the variable '{}' is written twice", v.name));
            }
            namer.push(first);
            for &(f, c) in &v.terms {
                if let Field::AtomX(p, ..) = f {
                    if m.phases[p].modes.is_some() {
                        return Err(format!("the variable '{}' moves '{}', which follows from the symmetry modes", v.name, address(f, m)));
                    }
                }
                // Variables may move the same value (two modes the same coordinate: the terms add); a parameter may not.
                if let Some(other) = rows.get(&f).and_then(|r| r.keys().find(|&&c| c < first_var)) {
                    return Err(format!("the variable '{}' moves '{}', which the parameter '{}' moves too: a variable is refined in its values' place", v.name, address(f, m), address(namer[*other], m)));
                }
                if rows.get(&f).is_some_and(|r| r.contains_key(&(first_var + k))) {
                    return Err(format!("the variable '{}' moves '{}' twice", v.name, address(f, m)));
                }
                if self.links.iter().any(|l| l.value == f) {
                    return Err(format!("the variable '{}' moves '{}', which follows a link", v.name, address(f, m)));
                }
                rows.entry(f).or_default().insert(first_var + k, c);
                order.push(f);
            }
        }
        let (mut free, col_name) = columns(m, rows, order, &namer, &self.links)?;
        for (p, &col) in free.iter_mut().zip(col_name.keys()) {
            if col >= first_var {
                let v = &self.variables[col - first_var];
                p.name = v.name.clone();
                p.value = Some(v.value);
            }
        }
        for v in &self.variables {
            if !free.iter().any(|p| p.name == v.name) {
                return Err(format!("the variable '{}' moves nothing the symmetry leaves free", v.name));
            }
        }
        for p in free.iter_mut() {
            p.damping = self.steering.damping_of(&p.name, p.fields[0].0);
        }
        for l in &self.linear {
            for (name, _) in &l.terms {
                if !free.iter().any(|p| &p.name == name) {
                    return Err(format!("the linear restraint '{}' names '{name}', which is not a free parameter", l.name));
                }
            }
        }
        for (name, _) in &self.steering.damping {
            if !free.iter().any(|p| &p.name == name) {
                return Err(format!("damping is given for '{name}', which is not a free parameter"));
            }
        }
        for &(f, b) in &self.behaviour {
            let (lo, hi) = self.limits.iter().find(|l| l.0 == f).map_or_else(|| limits(f, m), |l| (l.1, l.2));
            if b == Behaviour::Periodic && (lo.is_none() || hi.is_none()) {
                return Err(format!("'{}' is periodic but has not two bounds", address(f, m)));
            }
        }
        Ok(ParameterSet { free, linear: self.linear.clone(), user_limits: self.limits.clone(), behaviour: self.behaviour.clone(), links: self.links.clone(), steering: self.steering.clone() })
    }

    /// The refinement of one step of the plan (`steering::PlanStep`): the parameters and variables it frees, the
    /// step's damping over every parameter's, its cycles and convergence; the restraints on parameters it leaves
    /// fixed are dropped. An error if the step frees nothing.
    pub fn for_step(&self, k: usize, m: &Model) -> Result<Refinement, String> {
        let step = self.plan.get(k).ok_or(format!("the plan has {} steps", self.plan.len()))?;
        let mut r = self.clone();
        r.parameters.retain(|terms| step.frees(&address(terms[0].0, m), terms[0].0));
        r.variables.retain(|v| step.frees(&v.name, v.terms[0].0));
        if r.parameters.is_empty() && r.variables.is_empty() {
            return Err(format!("step {} ('{}') frees no parameter", k + 1, step.name));
        }
        let names: Vec<String> = r.parameters.iter().map(|t| address(t[0].0, m)).chain(r.variables.iter().map(|v| v.name.clone())).collect();
        r.linear.retain(|l| l.terms.iter().all(|(n, _)| names.contains(n)));
        r.steering.damping.retain(|(n, _)| names.contains(n));
        if let Some(d) = step.damping {
            for c in r.steering.class_damping.iter_mut() {
                *c *= d;
            }
        }
        r.cycles = step.cycles.or(self.cycles);
        r.eps = step.eps.or(self.eps);
        r.plan.clear();
        Ok(r)
    }

    /// The refinement a `.pcr` codes, with its cycles and convergence.
    pub fn from_pcr(pcr: &Pcr, m: &Model) -> Result<Refinement, String> {
        let (rows, order, namer, numbered, damping) = pcr_rows(pcr, m)?;
        let mut parameters: BTreeMap<usize, Vec<(Field, f64)>> = BTreeMap::new();
        let mut seen = std::collections::BTreeSet::new();
        for f in &order {
            if !seen.insert(*f) {
                continue;
            }
            for (&col, &c) in rows.get(f).into_iter().flatten() {
                parameters.entry(col).or_default().push((*f, c));
            }
        }
        for (col, terms) in parameters.iter_mut() {
            if let Some(k) = terms.iter().position(|t| t.0 == namer[*col]) {
                let t = terms.remove(k);
                terms.insert(0, t);
            }
        }
        let names: BTreeMap<usize, String> = parameters.iter().map(|(&c, t)| (c, address(t[0].0, m))).collect();
        let mut linear = Vec::new();
        for l in &pcr.linear_restraints {
            let terms = l
                .terms
                .iter()
                .map(|&(c, n)| numbered.get(&n).and_then(|col| names.get(col)).map(|name| (name.clone(), c)).ok_or(format!("linear restraint {}: parameter {n} is not refined", l.name)))
                .collect::<Result<Vec<_>, _>>()?;
            linear.push(LinearRestraint { name: l.name.clone(), terms, target: l.value, sigma: l.sigma });
        }
        let steering = pcr_steering(pcr, &damping, &names);
        Ok(Refinement { parameters: parameters.into_values().collect(), linear, limits: Vec::new(), cycles: Some(pcr.ncy.max(0) as usize), eps: Some(pcr.eps), steering, ..Default::default() })
    }
}

/// The field an address names in a model (the inverse of `address`); an address that names no field of the model is
/// an error with the address in it.
pub fn field_of(addr: &str, m: &Model) -> Result<Field, String> {
    m.all_fields().into_iter().find(|&f| address(f, m) == addr).ok_or(format!("'{addr}' names no value of this model"))
}

/// How a field is shown (the app's request 39, L276, L292): what the app needs to label a value without a list of the
/// core's fields of its own.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct FieldLabel {
    /// The parameter dictionary's id and name (`parameter_table`): "atom_position", "Atom position".
    pub id: &'static str,
    pub name: &'static str,
    /// The phase it belongs to (`Field::phase`), from 0; none for the pattern's values.
    pub phase: Option<usize>,
    /// The atom it belongs to: its index in the phase and its label.
    pub atom: Option<(usize, String)>,
    /// The component, where fields share the dictionary's id: x, y, z; β₁₁ … β₂₃; P1 … P4; S400 …; a mode's name; a
    /// background point's 2θ ("at 23.45°") or a series' order ("order 3").
    pub component: Option<String>,
}

impl Field {
    /// How this field of `m` is shown (`FieldLabel`): its dictionary id and name, phase, atom and component.
    pub fn describe(self, m: &Model) -> FieldLabel {
        let id = dictionary_id(self, m);
        let name = crate::parameter_table::PARAMETERS.iter().find(|p| p.id == id).map_or(id, |p| p.name);
        let atom = match self {
            Field::AtomX(p, a, _) | Field::AtomB(p, a) | Field::AtomOcc(p, a) | Field::AtomBeta(p, a, _) => match m.phases.get(p).map(|ph| &ph.intensities) {
                Some(Intensities::Structure(st, _)) => st.atoms.get(a).map(|at| (a, at.label.clone())),
                _ => None,
            },
            _ => None,
        };
        let component = match self {
            Field::Background(k) => Some(match &m.background {
                BackgroundModel::Points(pts) | BackgroundModel::Spline(pts) => pts.get(k).map_or(format!("point {}", k + 1), |p| if m.tof.is_some() { format!("at {:.1} µs", p.0) } else { format!("at {:.2}°", p.0) }),
                _ => format!("order {k}"),
            }),
            Field::BackgroundTerm(k) => Some(format!("term {}", m.background_terms.id(k).1)),
            Field::Asym(_, k) => Some(format!("P{}", k + 1)),
            Field::AtomX(_, _, k) => Some(["x", "y", "z"][k].to_string()),
            Field::AtomBeta(_, _, k) => Some(["β₁₁", "β₂₂", "β₃₃", "β₁₂", "β₁₃", "β₂₃"][k].to_string()),
            Field::Mode(p, k) => m.phases.get(p).and_then(|ph| ph.modes.as_ref()).and_then(|ms| ms.modes.get(k)).map(|md| md.name.clone()),
            Field::Rigid(p, g, k) => m.phases.get(p).and_then(|ph| ph.rigid.get(g)).map(|gr| rigid_scope(m, p, g, k).1 + " " + crate::rigid_body::value_component(gr, k)),
            Field::Quartic(p, k) => m.phases.get(p).and_then(|ph| ph.quartic_strain.as_ref()).and_then(|q| q.terms.get(k)).map(|(t, _)| format!("S{}{}{}", t[0], t[1], t[2])),
            Field::Intensity(p, k) => list_hkl(m, p, k).map(|h| format!("{} {} {}", h[0], h[1], h[2])),
            Field::Texture(p, k) => m.phases.get(p).and_then(|ph| ph.texture.as_ref()).map(|t| format!("axis {}", t.id(k).1)),
            _ => None,
        };
        FieldLabel { id, name, phase: self.phase(), atom, component }
    }
}

/// The parameter-dictionary id of a field (`parameter_table`; checked in `parameter_dictionary`).
pub fn dictionary_id(f: Field, m: &Model) -> &'static str {
    match f {
        Field::Zero => "zero_shift",
        Field::Sycos => "sample_displacement",
        Field::Sysin => "sample_transparency",
        Field::Ratio => "wavelength_ratio",
        Field::Wavelength => "wavelength_1",
        Field::Micro(k) => ["microabsorption_p0", "microabsorption_cp", "microabsorption_tau"][k],
        Field::Background(_) => match m.background {
            BackgroundModel::Points(_) | BackgroundModel::Spline(_) => "background_points",
            BackgroundModel::Polynomial { .. } => "background_polynomial",
            BackgroundModel::Chebyshev { .. } => "background_chebyshev",
            BackgroundModel::Fourier { .. } => "background_fourier",
        },
        Field::BackgroundTerm(k) => m.background_terms.id(k).0,
        Field::Scale(_) => "scale_factor",
        Field::Bov(_) => "overall_b",
        Field::Cell(_, k) => ["cell_a", "cell_b", "cell_c", "cell_alpha", "cell_beta", "cell_gamma"][k],
        Field::U(_) => "caglioti_u",
        Field::V(_) => "caglioti_v",
        Field::W(_) => "caglioti_w",
        Field::X(p) => {
            if m.phases[p].shape == Shape::PseudoVoigt {
                "pv_eta_slope"
            } else {
                "lorentzian_strain"
            }
        }
        Field::Y(_) => "lorentzian_size",
        Field::GaussSize(_) => "gaussian_size",
        Field::Eta0(_) => "pv_eta0",
        Field::Asym(..) => "berar_baldinozzi",
        Field::SL(_) => "axial_divergence_s_l",
        Field::DL(_) => "axial_divergence_d_l",
        Field::Texture(p, k) => m.phases.get(p).and_then(|ph| ph.texture.as_ref()).map_or("po_axis_ratio", |t| t.id(k).0),
        Field::Pref(_, 0) => "po_g1",
        Field::Pref(..) => "po_g2",
        Field::AtomX(..) => "atom_position",
        Field::AtomB(..) => "atom_biso",
        Field::AtomOcc(..) => "atom_occupancy",
        Field::AtomBeta(..) => "atom_betas",
        Field::Mode(..) => "mode_amplitude",
        Field::Rigid(p, g, k) => m.phases.get(p).and_then(|ph| ph.rigid.get(g)).map_or("rigid_origin", |gr| crate::rigid_body::value_id(gr, k)),
        Field::Quartic(..) => "quartic_strain",
        Field::QuarticShare(_) => "quartic_strain_lorentzian",
        Field::ShapeSize(_) => "shape_size",
        Field::Intensity(..) => "reflection_intensity",
        Field::Tof(k) => crate::tof::BANK_FIELDS[k],
        Field::TofProfile(_, k) => crate::tof::PROFILE_FIELDS[k],
    }
}

/// The indices of reflection k of a phase's list.
fn list_hkl(m: &Model, p: usize, k: usize) -> Option<[i32; 3]> {
    match &m.phases.get(p)?.intensities {
        Intensities::List(l) => l.get(k).map(|r| r.0),
        Intensities::Structure(..) => None,
    }
}

fn atom_label(m: &Model, p: usize, a: usize) -> String {
    match &m.phases[p].intensities {
        Intensities::Structure(st, _) => {
            let l = st.atoms[a].label.clone();
            // A label that an earlier atom of the phase already has takes its position: "Y:2".
            if st.atoms[..a].iter().any(|b| b.label == l) {
                format!("{l}:{}", a + 1)
            } else {
                l
            }
        }
        Intensities::List(_) => format!("atom{}", a + 1),
    }
}

/// The scope of a rigid group's k-th value: (address part, label): `group1` for its origin, orientation and TLS
/// values, its atom's label for an internal coordinate.
fn rigid_scope(m: &Model, p: usize, g: usize, k: usize) -> (String, String) {
    if crate::rigid_body::is_group_value(&m.phases[p].rigid[g], k) {
        (format!("group{}", g + 1), format!("group {}", g + 1))
    } else {
        let a = atom_label(m, p, m.phases[p].rigid[g].members[(k - 6) / 3].index);
        (a.clone(), a)
    }
}

/// The address of a field: scope (`pattern1`, `phaseN`, `phaseN.<atom label>`), dictionary id, and the component
/// where the id names several values (x, y, z; 11 … 23; P1 … P4; the background's k).
pub fn address(f: Field, m: &Model) -> String {
    let id = dictionary_id(f, m);
    match f {
        Field::Zero | Field::Sycos | Field::Sysin | Field::Ratio | Field::Wavelength | Field::Micro(_) | Field::Tof(_) => format!("pattern1.{id}"),
        Field::TofProfile(p, _) => format!("phase{}.{id}", p + 1),
        Field::Background(k) => format!("pattern1.{id}.{k}"),
        // A further term of the background by its number from 1: pattern1.background_debye_distance.2 (P23).
        Field::BackgroundTerm(k) => format!("pattern1.{id}.{}", m.background_terms.id(k).1),
        Field::Asym(p, k) => format!("phase{}.{id}.{}", p + 1, k + 1),
        Field::AtomX(p, a, k) => format!("phase{}.{}.{id}.{}", p + 1, atom_label(m, p, a), ["x", "y", "z"][k]),
        Field::AtomB(p, a) | Field::AtomOcc(p, a) => format!("phase{}.{}.{id}", p + 1, atom_label(m, p, a)),
        Field::AtomBeta(p, a, k) => format!("phase{}.{}.{id}.{}", p + 1, atom_label(m, p, a), ["11", "22", "33", "12", "13", "23"][k]),
        Field::Mode(p, k) => format!("phase{}.{}.{id}", p + 1, m.phases[p].modes.as_ref().map_or(String::new(), |ms| ms.modes[k].name.clone())),
        // A group's origin and orientation by its number, an internal coordinate by its atom: phase1.group1.
        // rigid_orientation.theta, phase1.A3.rigid_angle.angle (P41).
        Field::Rigid(p, g, k) => format!("phase{}.{}.{id}.{}", p + 1, rigid_scope(m, p, g, k).0, crate::rigid_body::value_component(&m.phases[p].rigid[g], k)),
        // The term by its exponents: phase1.quartic_strain.400.
        Field::Quartic(p, k) => format!("phase{}.{id}.{}", p + 1, m.phases[p].quartic_strain.as_ref().map_or(String::new(), |q| q.terms[k].0.iter().map(|e| e.to_string()).collect())),
        Field::QuarticShare(p) | Field::ShapeSize(p) => format!("phase{}.{id}", p + 1),
        // An axis of a texture by its number from 1: phase1.po_axis_fraction.1 (P30).
        Field::Texture(p, k) => format!("phase{}.{id}.{}", p + 1, m.phases[p].texture.as_ref().map_or(String::new(), |t| t.id(k).1)),
        // The reflection by its indices: phase1.reflection_intensity.1_-1_0 (P36).
        Field::Intensity(p, k) => format!("phase{}.{id}.{}", p + 1, list_hkl(m, p, k).map_or(String::new(), |h| format!("{}_{}_{}", h[0], h[1], h[2]))),
        Field::Scale(p) | Field::Bov(p) | Field::Cell(p, _) | Field::U(p) | Field::V(p) | Field::W(p) | Field::X(p) | Field::Y(p) | Field::GaussSize(p) | Field::Eta0(p) | Field::SL(p) | Field::DL(p) | Field::Pref(p, _) => {
            format!("phase{}.{id}", p + 1)
        }
    }
}

/// The bounds a field keeps during a refinement (P39), where the value has a physical sign: scale, B, occupancy, the
/// Gaussian and Lorentzian size terms, S/L and H/L, the ratio, G2 within [0, 1], a March ratio G1 above zero, the cell
/// lengths; and, over a resolution file (whose U, V, W, X, Y are the instrument's), the phase's U and X, which are
/// then the sample's strain. Listed in the dictionary's `limits`.
/// Where a value stands against its limits (the Parameters list's "at its limit", "outside its limits"; L145).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum LimitState {
    Inside,
    /// On a limit, within 10⁻¹² relative.
    AtLimit,
    /// Outside them: a value imported so, which may only move towards them (L053).
    Outside,
}

impl LimitState {
    pub fn of(v: f64, (lo, hi): (Option<f64>, Option<f64>)) -> LimitState {
        let on = |b: f64| (v - b).abs() < 1e-12 * (1.0 + v.abs());
        if lo.is_some_and(|b| v < b && !on(b)) || hi.is_some_and(|b| v > b && !on(b)) {
            LimitState::Outside
        } else if lo.is_some_and(on) || hi.is_some_and(on) {
            LimitState::AtLimit
        } else {
            LimitState::Inside
        }
    }
    /// The words the window shows; none inside.
    pub fn words(self) -> Option<&'static str> {
        match self {
            LimitState::Inside => None,
            LimitState::AtLimit => Some("at its limit"),
            LimitState::Outside => Some("outside its limits"),
        }
    }
}

pub fn limits(f: Field, m: &Model) -> (Option<f64>, Option<f64>) {
    let irf = !matches!(m.instrumental, InstrumentalWidths::None);
    match f {
        Field::AtomOcc(p, a) => (Some(0.0), full_site(m, p, a)),
        Field::Scale(_) | Field::AtomB(..) | Field::GaussSize(_) | Field::Y(_) | Field::SL(_) | Field::DL(_) | Field::Ratio | Field::ShapeSize(_) | Field::Intensity(..) => (Some(0.0), None),
        // A TOF phase's extinction and absorption (P05): none below zero.
        Field::TofProfile(_, 17 | 18) => (Some(0.0), None),
        Field::U(p) | Field::X(p) if irf && m.phases[p].shape == Shape::ThompsonCoxHastings => (Some(0.0), None),
        Field::Pref(_, 1) | Field::QuarticShare(_) => (Some(0.0), Some(1.0)),
        Field::Pref(_, 0) if m.march => (Some(1e-3), None),
        // A texture's fractions within [0, 1], its March ratios above zero (P30): an odd value is a fraction.
        Field::Texture(_, k) if k % 2 == 1 => (Some(0.0), Some(1.0)),
        Field::Texture(..) => (Some(1e-3), None),
        Field::Cell(_, k) if k < 3 => (Some(1e-3), None),
        Field::Wavelength => (Some(1e-3), None),
        // FullProf's Max_Amplitude: the box |Q| ≤ value on every amplitude (note of 16 January 2009).
        Field::Mode(p, _) => match m.phases[p].modes.as_ref().and_then(|ms| ms.max_amplitude) {
            Some(b) => (Some(-b), Some(b)),
            None => (None, None),
        },
        _ => (None, None),
    }
}

/// The steering a `.pcr` codes: each lone multiplier's |a| as its parameter's damping, the relaxation factors as the
/// classes' (L107).
fn pcr_steering(pcr: &Pcr, damping: &BTreeMap<usize, f64>, names: &BTreeMap<usize, String>) -> Steering {
    let mut s = Steering::default();
    if let Some(r) = pcr.relax {
        s.class_damping = r;
    }
    s.damping = damping.iter().filter_map(|(c, &d)| names.get(c).map(|n| (n.clone(), d))).collect();
    s
}

/// The Occ of a full site for atom a of phase p (W12): FullProf's Occ = o·m/M (o the chemical occupancy, m the site's
/// multiplicity, M the general position's; manual p. 100), so o ≤ 1 is Occ ≤ m/M; for a file whose occupancies follow
/// the primitive cell's count (`bonds::chemical_occupancies`), m/c. None when the structure's occupancies follow neither
/// (files that scale them): no upper bound is then imposed.
fn full_site(m: &Model, p: usize, a: usize) -> Option<f64> {
    let Intensities::Structure(st, _) = &m.phases[p].intensities else { return None };
    let chem = crate::bonds::chemical_occupancies(st).ok()?;
    if chem.iter().any(|&o| o > 1.0 + 1e-3) || st.atoms[a].occ <= 0.0 || chem[a] <= 0.0 {
        return None;
    }
    Some(st.atoms[a].occ / chem[a])
}

/// Decode a codeword (manual p. 13): (p, a), p = 0 for a parameter numbered automatically (Aut = 1, |C| < 10).
fn codeword(c: f64) -> Option<(u32, f64)> {
    if c == 0.0 {
        return None;
    }
    let ac = c.abs();
    let p = (ac / 10.0 + 1e-9).floor();
    let a = c.signum() * (ac - 10.0 * p);
    if a.abs() < 1e-9 {
        return None;
    }
    Some((p as u32, a))
}

/// The rows of M: field → (column → coefficient).
type Rows = BTreeMap<Field, BTreeMap<usize, f64>>;

fn row(rows: &Rows, f: Field) -> BTreeMap<usize, f64> {
    rows.get(&f).cloned().unwrap_or_default()
}

fn combine(terms: &[(f64, &BTreeMap<usize, f64>)]) -> BTreeMap<usize, f64> {
    let mut out = BTreeMap::new();
    for (c, r) in terms {
        for (&k, &v) in r.iter() {
            *out.entry(k).or_insert(0.0) += c * v;
        }
    }
    out.retain(|_, v| v.abs() > 1e-14);
    out
}

/// The rows of a `.pcr`'s codewords in the default model's values (field → column → coefficient), the fields in the
/// order they were coded, the field that names each column, FullProf's parameter numbers → columns, and the damping of
/// the columns whose multiplier is a lone one.
type PcrRows = (Rows, Vec<Field>, Vec<Field>, BTreeMap<u32, usize>, BTreeMap<usize, f64>);

fn pcr_rows(pcr: &Pcr, m: &Model) -> Result<PcrRows, String> {
        let mut coded: Vec<(Field, f64)> = Vec::new();
        let mut refused = Vec::new();
        let z = &pcr.zero_codes;
        if m.tof.is_some() {
            // A TOF pattern (P05): line 14's Zero, Dtt1, Dtt2, Dtt_1overd and line 16's Zerot, Dtt1t, Dtt2t, x-cross,
            // Width (Npr = 10).
            for (k, f) in [Field::Zero, Field::Tof(0), Field::Tof(1), Field::Tof(2)].into_iter().enumerate() {
                coded.push((f, z.get(k).copied().unwrap_or(0.0)));
            }
            for (k, &c) in pcr.thermal_codes.iter().enumerate().take(5) {
                coded.push((Field::Tof(3 + k), c));
            }
        } else {
            for (k, f) in [Field::Zero, Field::Sycos, Field::Sysin].into_iter().enumerate() {
                coded.push((f, z.get(k).copied().unwrap_or(0.0)));
            }
            // Line 14's Lambda (P18): λ₁, with λ₂ when the pattern has one wavelength (`Field::Wavelength`).
            if let Some(&c) = z.get(3).filter(|&&c| c != 0.0) {
                coded.push((Field::Wavelength, c));
            }
        }
        for (k, &c) in pcr.microabsorption_codes.iter().enumerate() {
            coded.push((Field::Micro(k), c));
        }
        if pcr.flag("Nba") == 1 {
            // Nba = 1: the codes are those of b1–b4 of the .bac background (eq. 3.6), not of its points: the file's
            // scale, offset, stretch and shift (the app's request 47), the last four values of the further terms.
            let first = m.background_terms.len().saturating_sub(4);
            if m.background_terms.file.is_some() {
                for (k, &c) in pcr.background_codes.iter().take(4).enumerate() {
                    coded.push((Field::BackgroundTerm(first + k), c));
                }
            }
        } else {
            // Nba = −1 and −3 (P23): the codes past the polynomial's are those of the further terms, in their order
            // (the Debye-like terms' amplitudes, then their distances; the inverse powers).
            let n = match pcr.flag("Nba") {
                -1 => 6,
                -3 => 9,
                _ => usize::MAX,
            };
            let nd = m.background_terms.debye.len();
            for (k, &c) in pcr.background_codes.iter().enumerate() {
                let f = match k.checked_sub(n) {
                    None => Field::Background(k),
                    Some(j) if j < nd => Field::BackgroundTerm(2 * j),
                    Some(j) if j < 2 * nd => Field::BackgroundTerm(2 * (j - nd) + 1),
                    Some(j) => Field::BackgroundTerm(j),
                };
                coded.push((f, c));
            }
            // MORE = 2, 3: the profiles' coefficients.
            for (j, &c) in pcr.profile_codes.iter().enumerate() {
                coded.push((Field::BackgroundTerm(2 * nd + m.background_terms.inverse.len() + j), c));
            }
            // FullProf's own end point of a spline (Nba < 0, L081) takes the last point's value: the same codeword.
            if let (crate::powder::BackgroundModel::Spline(q), Some(&c)) = (&m.background, pcr.background_codes.last()) {
                if pcr.flag("Nba") < 0 && q.len() == pcr.background_codes.len() + 1 {
                    coded.push((Field::Background(q.len() - 1), c));
                }
            }
        }
        for cmd in &pcr.commands {
            let t = cmd.trim();
            let mut w = t.split_whitespace();
            match w.next().map(str::to_ascii_uppercase).as_deref() {
                Some("VARY") => {
                    for k in w {
                        // `VARY back_n` of a background term its codeword already refines changes nothing: FullProf's
                        // CeO2_PEARL.out lists the .pcr's 20 coded parameters, Bck_0 … Bck_6, after VARY back_7
                        // (back_n read as the n-th term, Bck_{n−1}, believed). Any other is refused.
                        let back = k.to_ascii_lowercase().strip_prefix("back_").and_then(|n| n.parse::<usize>().ok());
                        if k.eq_ignore_ascii_case("ratio") {
                            coded.push((Field::Ratio, 1.0));
                        } else if back.is_some_and(|n| n >= 1 && pcr.background_codes.get(n - 1).is_some_and(|&c| c != 0.0)) {
                        } else {
                            refused.push(format!("VARY {k}"));
                        }
                    }
                }
                Some("FIX") if w.next().is_some() => refused.push(t.to_string()),
                _ => {}
            }
        }
        if m.phases.len() != pcr.phases.len() {
            return Err(format!("{} phases in the model for {} in the .pcr", m.phases.len(), pcr.phases.len()));
        }
        for (p, (pp, mp)) in pcr.phases.iter().zip(&m.phases).enumerate() {
            let c = |v: &Vec<f64>, k: usize| v.get(k).copied().unwrap_or(0.0);
            coded.push((Field::Scale(p), c(&pp.scale_codes, 0)));
            if let (Some(t), Some(_)) = (&pp.tof, &mp.tof) {
                // A TOF phase (P05, pp. 110–113): Extinc in Shape1's place; Sig-2 Sig-1 Sig-0 Sigma-Q (line 27-1),
                // Gam-2 Gam-1 Gam-0 (27-3), alph0 beta0 alph1 beta1 alphQ betaQ after Pref1 Pref2 (30-1), alph0t
                // beta0t alph1t beta1t (31-1), ABSCOR1 (33).
                coded.push((Field::TofProfile(p, 17), c(&pp.scale_codes, 1)));
                for (k, f) in [2, 1, 0, 3].into_iter().enumerate() {
                    coded.push((Field::TofProfile(p, f), c(&t.sigma_codes, k)));
                }
                for (k, f) in [6, 5, 4].into_iter().enumerate() {
                    coded.push((Field::TofProfile(p, f), c(&t.gamma_codes, k)));
                }
                for (k, f) in [7, 10, 8, 11, 9, 12].into_iter().enumerate() {
                    coded.push((Field::TofProfile(p, f), c(&t.shape_codes, 2 + k)));
                }
                if mp.tof.as_ref().is_some_and(|x| x.thermal.is_some()) {
                    for k in 0..4 {
                        coded.push((Field::TofProfile(p, 13 + k), c(&t.thermal_codes, k)));
                    }
                }
                coded.push((Field::TofProfile(p, 18), t.absorption_codes[0]));
                let unread = [(&t.sigma_codes, 4..7, "Iso-GStrain, Iso-GSize, Ani-LSize"), (&t.gamma_codes, 3..5, "Iso-LorStrain, Iso-LorSize")];
                for (v, r, names) in unread {
                    if r.into_iter().any(|k| c(v, k) != 0.0) {
                        refused.push(format!("phase {}: the TOF size and strain terms ({names})", p + 1));
                    }
                }
                if t.absorption_codes[1] != 0.0 {
                    refused.push(format!("phase {}: ABSCOR2", p + 1));
                }
            } else if c(&pp.scale_codes, 1) != 0.0 {
                if mp.shape == Shape::PseudoVoigt {
                    coded.push((Field::Eta0(p), c(&pp.scale_codes, 1)));
                } else {
                    refused.push(format!("phase {}: Shape1 of a non-pseudo-Voigt profile", p + 1));
                }
            }
            coded.push((Field::Bov(p), c(&pp.scale_codes, 2)));
            let hexagonal = pp.flag("Str") == 0 && pp.strain_model == 8;
            if hexagonal {
                // Strain model 8 as the quartic form (`powder::hexagonal_strain`): Saa, Scc, Cac refined as S_400,
                // S_004 and S_112, one to one where Saa, Scc > 0 (the same minimum, other parameters).
                if let Some(q) = &mp.quartic_strain {
                    for (k, t) in [[4, 0, 0], [0, 0, 4], [1, 1, 2]].iter().enumerate() {
                        if let Some(i) = q.terms.iter().position(|x| x.0 == *t) {
                            coded.push((Field::Quartic(p, i), c(&pp.scale_codes, 3 + k)));
                        }
                    }
                }
            } else if (3..6).any(|k| c(&pp.scale_codes, k) != 0.0) {
                refused.push(format!("phase {}: strain parameters Str1–Str3 (FullProf's strain models, not calculated: L211)", p + 1));
            }
            // The quartic form of strain (P25): each S_HKL in the model's order of terms, and ξ.
            if let (Some(q), Some(order)) = (&mp.quartic_strain, crate::microstructure::fullprof_strain_model(pp.strain_model)) {
                for (k, (t, _)) in q.terms.iter().enumerate() {
                    // 4/m's S_310, which FullProf has not, is fixed (L267).
                    let code = order.iter().position(|o| o == t).map_or(0.0, |i| c(&pp.quartic_codes, i));
                    coded.push((Field::Quartic(p, k), code));
                }
                coded.push((Field::QuarticShare(p), pp.xi.map_or(0.0, |x| x.1)));
            }
            for (k, f) in [Field::U(p), Field::V(p), Field::W(p), Field::X(p), Field::Y(p), Field::GaussSize(p)].into_iter().enumerate() {
                coded.push((f, c(&pp.width_codes, k)));
            }
            // LorSiz of the anisotropic size (L267); refused without one (Model::from_pcr refuses the phase then).
            if mp.shape_size.is_some() {
                coded.push((Field::ShapeSize(p), c(&pp.width_codes, 6)));
            } else if c(&pp.width_codes, 6) != 0.0 {
                refused.push(format!("phase {}: LorSiz without a size model", p + 1));
            }
            for k in 0..6 {
                coded.push((Field::Cell(p, k), c(&pp.cell_codes, k)));
            }
            coded.push((Field::Pref(p, 0), c(&pp.pref_codes, 0)));
            coded.push((Field::Pref(p, 1), c(&pp.pref_codes, 1)));
            // Nor = 2, 3 (P30): each axis's r and fraction, the last axis's r alone.
            for (k, &code) in pp.texture_codes.iter().enumerate() {
                coded.push((Field::Texture(p, k), code));
            }
            for k in 0..4 {
                coded.push((Field::Asym(p, k), c(&pp.pref_codes, 2 + k)));
            }
            coded.push((Field::SL(p), c(&pp.pref_codes, 6)));
            coded.push((Field::DL(p), c(&pp.pref_codes, 7)));
            // Rigid groups (P41): θ, φ (codes 6, 7 of line 25-2), the internal coordinates, X0 Y0 Z0 and χ (line 25-4);
            // a group atom's x, y, z are FullProf's information only (manual p. 35).
            for (g, gr) in mp.rigid.iter().enumerate() {
                let head = &pp.atoms[gr.members[0].index];
                for (k, code) in [(3, c(&head.codes, 6)), (4, c(&head.codes, 7)), (0, c(&head.rigid_codes, 3)), (1, c(&head.rigid_codes, 4)), (2, c(&head.rigid_codes, 5)), (5, c(&head.rigid_codes, 6))] {
                    coded.push((Field::Rigid(p, g, k), code));
                }
                for (mi, m) in gr.members.iter().enumerate() {
                    for cc in 0..3 {
                        coded.push((Field::Rigid(p, g, 6 + 3 * mi + cc), c(&pp.atoms[m.index].rigid_codes, cc)));
                    }
                }
                // TLS (L367): the codewords of the phase's lines "t11 value code" … (ties such as s22 = −s11 by sign).
                if gr.displacements.is_some() {
                    if let Some((_, codes)) = crate::rigid_body::tls_of(pp) {
                        let n = 6 + 3 * gr.members.len();
                        coded.extend(codes.iter().enumerate().map(|(i, &code)| (Field::Rigid(p, g, n + i), code)));
                    }
                }
            }
            for (a, at) in pp.atoms.iter().enumerate() {
                let grouped = mp.rigid.iter().any(|g| g.members.iter().any(|m| m.index == a));
                for k in 0..3 {
                    coded.push((Field::AtomX(p, a, k), if grouped { 0.0 } else { c(&at.codes, k) }));
                }
                coded.push((Field::AtomB(p, a), c(&at.codes, 3)));
                coded.push((Field::AtomOcc(p, a), c(&at.codes, 4)));
                if let Some(bc) = &at.beta_codes {
                    for k in 0..6 {
                        coded.push((Field::AtomBeta(p, a, k), c(bc, k)));
                    }
                }
            }
            if let (Some(pm), Some(_)) = (&pp.modes, &mp.modes) {
                for (k, &code) in pm.codes.iter().enumerate() {
                    coded.push((Field::Mode(p, k), code));
                }
            }
        }
        if !refused.is_empty() {
            return Err(format!("not refined by the core yet: {}", refused.join("; ")));
        }
        // Columns: one per parameter number, one per automatic codeword.
        let mut rows: Rows = BTreeMap::new();
        let mut order: Vec<Field> = Vec::new();
        let mut numbered: BTreeMap<u32, usize> = BTreeMap::new();
        // The field each column's codeword was first written on: it names the free parameter.
        let mut namer: Vec<Field> = Vec::new();
        let mut ncol = 0;
        for &(f, cw) in &coded {
            let Some((pn, a)) = codeword(cw) else { continue };
            let col = if pn == 0 {
                ncol += 1;
                namer.push(f);
                ncol - 1
            } else {
                *numbered.entry(pn).or_insert_with(|| {
                    ncol += 1;
                    namer.push(f);
                    ncol - 1
                })
            };
            rows.entry(f).or_default().insert(col, a);
            order.push(f);
        }
        // Multipliers by their meaning (W12, L107): a column carried by one value alone is that value refined with its
        // shift multiplied by |a| (FullProf, manual p. 13: "the factor by which the computed shift … will be
        // multiplied"): the core takes it as damping, the coefficient keeping a's sign. Values the symmetry ties to it
        // (a cubic cell's b and c, coded with a's codeword) do not make it shared. Shared multipliers are the
        // coefficients of a link between the values, and stay so.
        let mut damping = BTreeMap::new();
        for col in 0..ncol {
            let carriers: Vec<Field> = rows.iter().filter(|(_, r)| r.contains_key(&col)).map(|(f, _)| *f).collect();
            let free: Vec<Field> = carriers.iter().copied().filter(|&f| m.fixed_by(f).is_none()).collect();
            if let [f] = free[..] {
                let a = rows[&f][&col];
                if (a.abs() - 1.0).abs() > 1e-9 {
                    for g in carriers {
                        let c = rows.get_mut(&g).unwrap().get_mut(&col).unwrap();
                        *c /= a.abs();
                    }
                    damping.insert(col, a.abs());
                }
            }
        }
        // Into the default model's values (L050).
        let conv = m.convention;
        for (p, ph) in m.phases.iter().enumerate() {
            if !conv.lorentz_polarisation && !m.neutrons && !matches!(m.ilo, -1 | 2 | 3) && matches!(ph.intensities, Intensities::Structure(..)) {
                let r = row(&rows, Field::Scale(p));
                rows.insert(Field::Scale(p), combine(&[(2.0, &r)]));
            }
            if !conv.berar_baldinozzi {
                for (lo, hi) in [(0, 1), (2, 3)] {
                    let (a, b) = (row(&rows, Field::Asym(p, lo)), row(&rows, Field::Asym(p, hi)));
                    rows.insert(Field::Asym(p, lo), combine(&[(1.0, &a), (3.0, &b)]));
                }
            }
            let pd_l = pcr.phases[p].d_l;
            if !conv.axial_divergence && pd_l < 0.0002 && (pcr.phases[p].s_l != 0.0 || pd_l != 0.0) {
                let (s, d) = (row(&rows, Field::SL(p)), row(&rows, Field::DL(p)));
                rows.insert(Field::SL(p), combine(&[(0.5, &s), (-0.5, &d)]));
                rows.insert(Field::DL(p), combine(&[(0.5, &s), (0.5, &d)]));
            }
            if !conv.preferred_orientation && !m.march {
                let g = row(&rows, Field::Pref(p, 0));
                rows.insert(Field::Pref(p, 0), combine(&[(-1.0, &g)]));
            }
        }
        // A value the translations give a row keeps it, though it had no codeword of its own (D_L moving with S_L when
        // D_L < 0.0002; P1 with P2): before 0.24.0 such a value was left out of its parameter (L110).
        for f in rows.keys() {
            if !order.contains(f) {
                order.push(*f);
            }
        }
        // A column that moves one value alone after the translations is that value refined: its coefficient is only a
        // change of variable (f = c·q), written ±1 so that the set is the same however the file wrote it.
        for col in 0..ncol {
            let carriers: Vec<Field> = rows.iter().filter(|(_, r)| r.contains_key(&col)).map(|(f, _)| *f).collect();
            if let [f] = carriers[..] {
                let c = rows.get_mut(&f).unwrap().get_mut(&col).unwrap();
                *c = c.signum();
            }
        }
        Ok((rows, order, namer, numbered, damping))
}

/// FullProf's parameter numbers that a `.pcr`'s codewords use (the p of C = sign(a)(10p + |a|)), each with the fields
/// that carry it, in the default model built from the file; the ratio of `VARY ratio` and automatic codewords (Aut = 1)
/// have no number and are not listed. FullProf needs them to be 1…N with no hole, N the number line 13 declares
/// ("Hole in the matrix!, missing element", L071).
pub fn pcr_parameter_numbers(pcr: &Pcr, m: &Model) -> Result<BTreeMap<u32, Vec<Field>>, String> {
    let (rows, _, _, numbered, _) = pcr_rows(pcr, m)?;
    let mut out: BTreeMap<u32, Vec<Field>> = BTreeMap::new();
    for (&n, &col) in &numbered {
        let fields: Vec<Field> = rows.iter().filter(|(_, r)| r.contains_key(&col)).map(|(f, _)| *f).collect();
        out.insert(n, fields);
    }
    Ok(out)
}

/// The free parameters of a `.pcr` (`ParameterSet::from_pcr`) with FullProf's parameter number of each (the p of its
/// codeword C = sign(a)(10p + |a|); None for one without: the ratio of `VARY ratio`, automatic codewords). The patterns
/// of a joint refinement share a parameter where they share its number (P06, `joint`).
pub fn pcr_numbered_set(pcr: &Pcr, m: &Model) -> Result<(ParameterSet, Vec<Option<u32>>), String> {
    let set = ParameterSet::from_pcr(pcr, m)?;
    let (rows, order, namer, numbered, _) = pcr_rows(pcr, m)?;
    let (_, col_name) = columns(m, rows, order, &namer, &[])?;
    let number_of: BTreeMap<usize, u32> = numbered.iter().map(|(&n, &c)| (c, n)).collect();
    let numbers = set.free.iter().map(|p| col_name.iter().find(|(_, name)| **name == p.name).and_then(|(c, _)| number_of.get(c).copied())).collect();
    Ok((set, numbers))
}

impl ParameterSet {
    /// A phase's cell volume (Å³) and its su from `covariance`, the covariance of the free parameters (`Refinement::
    /// covariance`, C χ²): the cell's correlations and its symmetry's ties included (a cubic cell's volume moves as
    /// 3a²·δa). The one su of a cell volume the core gives (P43): the results file's `_cell_volume`, the series
    /// table's, the summary's. 0 without a covariance of the set's size.
    pub fn cell_volume(&self, model: &Model, phase: usize, covariance: &[Vec<f64>]) -> (f64, f64) {
        let c = model.phases[phase].cell;
        let v0 = c.volume();
        if covariance.len() != self.free.len() || self.free.is_empty() {
            return (v0, 0.0);
        }
        // ∂V/∂(cell parameter k), by central differences.
        let grad: Vec<f64> = (0..6)
            .map(|k| {
                let (mut cp, mut cm) = (c, c);
                let x = [c.a, c.b, c.c, c.alpha, c.beta, c.gamma][k];
                let h = 1e-6 * x.abs().max(1.0);
                *[&mut cp.a, &mut cp.b, &mut cp.c, &mut cp.alpha, &mut cp.beta, &mut cp.gamma][k] = x + h;
                *[&mut cm.a, &mut cm.b, &mut cm.c, &mut cm.alpha, &mut cm.beta, &mut cm.gamma][k] = x - h;
                (cp.volume() - cm.volume()) / (2.0 * h)
            })
            .collect();
        let mut var = 0.0;
        for j in 0..6 {
            for k in 0..6 {
                var += grad[j] * grad[k] * self.field_covariance(covariance, Field::Cell(phase, j), Field::Cell(phase, k));
            }
        }
        (v0, var.max(0.0).sqrt())
    }

    /// The free parameters of a `.pcr` file, for the model `Model::from_pcr` (or `from_pcr_fullprof`) built from it
    /// with the same phases (magnetic phases removed from both). Refused with the list of what the core does not
    /// refine yet: wavelength, strain and size models, LorSiz, Shape1 of a TCH phase, VARY keywords other than ratio.
    pub fn from_pcr(pcr: &Pcr, m: &Model) -> Result<ParameterSet, String> {
        let (rows, order, namer, numbered, damping) = pcr_rows(pcr, m)?;
        let (mut free, col_name) = columns(m, rows, order, &namer, &[])?;
        let steering = pcr_steering(pcr, &damping, &col_name);
        for p in free.iter_mut() {
            p.damping = steering.damping_of(&p.name, p.fields[0].0);
        }
        // Linear restraints (NLI): FullProf's parameter numbers are its codewords' p.
        let mut linear = Vec::new();
        for l in &pcr.linear_restraints {
            let terms = l
                .terms
                .iter()
                .map(|&(c, n)| numbered.get(&n).and_then(|col| col_name.get(col)).map(|name| (name.clone(), c)).ok_or(format!("linear restraint {}: parameter {n} is not refined", l.name)))
                .collect::<Result<Vec<_>, _>>()?;
            linear.push(LinearRestraint { name: l.name.clone(), terms, target: l.value, sigma: l.sigma });
        }
        Ok(ParameterSet { free, linear, user_limits: Vec::new(), behaviour: Vec::new(), links: Vec::new(), steering })
    }

    /// A parameter set from the user's choices (the app's request 2, L058): the addresses set free, equalities
    /// `(a, m, b)` that make field a move as m times free field b (several add up), and limits that replace the
    /// engine's for a field.
    /// The crystal system's and the sites' constraints are applied as `from_pcr` applies them; an address that names
    /// no field, a value the symmetry fixes or ties, or an equality to an address that is not free is an error.
    pub fn from_addresses(free: &[String], equal: &[(String, f64, String)], limits: &[(String, Option<f64>, Option<f64>)], m: &Model) -> Result<ParameterSet, String> {
        let mut rows: Rows = BTreeMap::new();
        let mut order = Vec::new();
        let mut namer = Vec::new();
        for (col, a) in free.iter().enumerate() {
            let f = field_of(a, m)?;
            if let Some(by) = m.fixed_by(f) {
                return Err(format!("'{a}' follows from the {by}: free the value it follows from"));
            }
            if rows.contains_key(&f) {
                return Err(format!("'{a}' is given twice"));
            }
            rows.entry(f).or_default().insert(col, 1.0);
            order.push(f);
            namer.push(f);
        }
        for (a, mult, b) in equal {
            let fa = field_of(a, m)?;
            field_of(b, m)?;
            let col = free.iter().position(|x| x == b).ok_or(format!("'{a}' is made equal to '{b}', which is not free"))?;
            let row = rows.entry(fa).or_default();
            if row.contains_key(&col) {
                return Err(format!("'{a}' already moves with '{b}'"));
            }
            // A field may move with several free parameters (the import's P1 + 3·P2): the terms add.
            row.insert(col, *mult);
            order.push(fa);
        }
        let mut user_limits = Vec::new();
        for (a, lo, hi) in limits {
            if let (Some(l), Some(h)) = (lo, hi) {
                if l > h {
                    return Err(format!("'{a}': its lower limit {l} is above its upper limit {h}"));
                }
            }
            user_limits.push((field_of(a, m)?, *lo, *hi));
        }
        let (free, _) = columns(m, rows, order, &namer, &[])?;
        Ok(ParameterSet { free, linear: Vec::new(), user_limits, behaviour: Vec::new(), links: Vec::new(), steering: Steering::default() })
    }

    /// The covariance of two fields' values from that of the free parameters (`lsq::Refinement::covariance`): each
    /// field moves as Σ_j a_j q_j, so cov(f, g) = Σ_jk a_fj a_gk C_jk. Zero for a field no parameter moves.
    pub fn field_covariance(&self, cov: &[Vec<f64>], f: Field, g: Field) -> f64 {
        let terms = |x: Field| -> Vec<(usize, f64)> { self.free.iter().enumerate().flat_map(|(j, p)| p.fields.iter().filter(move |(h, _)| *h == x).map(move |(_, a)| (j, *a))).collect() };
        let (tf, tg) = (terms(f), terms(g));
        tf.iter().flat_map(|&(j, a)| tg.iter().map(move |&(k, b)| a * b * cov[j][k])).sum()
    }

    /// Set each linked value from the value it follows (x = c·p + d), in the order the links are written; returns the
    /// values changed, with their old and new values.
    pub fn enforce_links(&self, m: &mut Model) -> Vec<(Field, f64, f64)> {
        let mut out = Vec::new();
        for l in &self.links {
            let (old, new) = (m.get(l.value), l.coefficient * m.get(l.follows) + l.offset);
            if (old - new).abs() > 1e-12 * (1.0 + new.abs()) {
                m.set(l.value, new);
                out.push((l.value, old, new));
            }
        }
        out
    }

    /// What a bounded field does at its bounds.
    pub fn behaviour_of(&self, f: Field) -> Behaviour {
        self.behaviour.iter().find(|b| b.0 == f).map_or(Behaviour::Stop, |b| b.1)
    }

    /// Each free parameter as the app shows it (W12): its fields with their coefficients, the links that follow its
    /// fields, each field's bounds and behaviour, its damping and each field's maximum shift per cycle.
    pub fn describe(&self, m: &Model) -> Vec<FreeDescription> {
        self.free
            .iter()
            .map(|p| FreeDescription {
                name: p.name.clone(),
                class: crate::steering::class_of(p.fields[0].0).word(),
                variable: p.value.is_some(),
                damping: p.damping,
                fields: p
                    .fields
                    .iter()
                    .map(|&(f, c)| {
                        let (min, max) = self.limits_of(f, m);
                        FieldDescription {
                            address: address(f, m),
                            coefficient: c,
                            min,
                            max,
                            behaviour: self.behaviour_of(f).word(),
                            max_shift: self.steering.max_shift_of(f, m),
                            link: self.links.iter().find(|l| l.value == f).map(|l| (address(l.follows, m), l.coefficient, l.offset)),
                        }
                    })
                    .collect(),
            })
            .collect()
    }

    /// The limits of a field in this set: the user's if given, else the engine's (`limits`).
    pub fn limits_of(&self, f: Field, m: &Model) -> (Option<f64>, Option<f64>) {
        self.user_limits.iter().find(|u| u.0 == f).map_or_else(|| limits(f, m), |u| (u.1, u.2))
    }

    /// Every field the set moves.
    pub fn fields(&self) -> Vec<Field> {
        let mut v: Vec<Field> = self.free.iter().flat_map(|p| p.fields.iter().map(|f| f.0)).collect();
        v.sort();
        v.dedup();
        v
    }
}

/// A free parameter as `ParameterSet::describe` gives it.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct FreeDescription {
    pub name: String,
    /// FullProf's relaxation class of its first value.
    pub class: &'static str,
    /// A new variable (`steering::Variable`) rather than a value of the model.
    pub variable: bool,
    /// Its own damping times its class's.
    pub damping: f64,
    pub fields: Vec<FieldDescription>,
}

/// A value a free parameter moves.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct FieldDescription {
    pub address: String,
    /// How much it moves per unit of the parameter.
    pub coefficient: f64,
    pub min: Option<f64>,
    pub max: Option<f64>,
    /// "stop", "reflect" or "periodic".
    pub behaviour: &'static str,
    /// The largest shift per cycle, in its unit (None: no limit).
    pub max_shift: Option<f64>,
    /// The link it follows: (the value followed, coefficient, offset).
    pub link: Option<(String, f64, f64)>,
}

/// The free parameters from the rows of the fields (field → column → coefficient): the crystal system's cell rules and
/// each atom's site constraints applied, then one free parameter per column that moves some field, named after the
/// field `namer[column]` (its first). Also returns each column's parameter name.
fn columns(m: &Model, mut rows: Rows, mut order: Vec<Field>, namer: &[Field], links: &[Link]) -> Result<(Vec<Free>, BTreeMap<usize, String>), String> {
    // FullProf refines a single wavelength (manual p. 83); with a doublet, which line moves and whether λ₂/λ₁ is kept
    // is not said, so it is not refined (P18).
    if rows.get(&Field::Wavelength).is_some_and(|r| !r.is_empty()) && m.lambda[1] != m.lambda[0] {
        return Err(format!("the wavelength is refined for a pattern of one wavelength only: this one has λ₁ = {} and λ₂ = {} (FullProf refines a single wavelength, manual p. 83)", m.lambda[0], m.lambda[1]));
    }
    // Symmetry: the cell, and each atom's site.
    for (p, ph) in m.phases.iter().enumerate() {
        for (k, rule) in cell_rules(ph).into_iter().enumerate() {
            match rule {
                CellRule::Free => {}
                CellRule::Fixed => {
                    rows.remove(&Field::Cell(p, k));
                }
                CellRule::Equal(src) => {
                    let r = row(&rows, Field::Cell(p, src));
                    rows.insert(Field::Cell(p, k), r);
                    order.push(Field::Cell(p, k));
                }
            }
        }
        if let Intensities::Structure(st, _) = &ph.intensities {
            for (a, at) in st.atoms.iter().enumerate() {
                let site = Site::new(&ph.group, at.x).map_err(|e| format!("phase {} atom {}: {e}", p + 1, at.label))?;
                let apply = |rows: &mut Rows, order: &mut Vec<Field>, cons: &crate::site::Constraints, field: &dyn Fn(usize) -> Field| {
                    let old: Vec<BTreeMap<usize, f64>> = (0..cons.names.len()).map(|i| row(rows, field(i))).collect();
                    for i in 0..cons.names.len() {
                        let terms: Vec<(f64, &BTreeMap<usize, f64>)> = cons.coef[i].iter().zip(&cons.free).map(|(c, &fv)| (c.value(), &old[fv])).collect();
                        let r = combine(&terms);
                        if r.is_empty() {
                            rows.remove(&field(i));
                        } else {
                            rows.insert(field(i), r);
                            order.push(field(i));
                        }
                    }
                };
                apply(&mut rows, &mut order, &site.coordinates, &|i| Field::AtomX(p, a, i));
                if at.beta.is_some() {
                    apply(&mut rows, &mut order, &site.adp, &|i| Field::AtomBeta(p, a, i));
                }
            }
        }
    }
    // Links (W12): the value moves as c times the value it follows, after the symmetry's ties.
    for l in links {
        let (a, b) = (address(l.value, m), address(l.follows, m));
        if l.value == l.follows {
            return Err(format!("'{a}' is linked to itself"));
        }
        if rows.contains_key(&l.value) {
            return Err(format!("'{a}' follows '{b}' by a link and is also refined (or tied by the symmetry)"));
        }
        if let Some(by) = m.fixed_by(l.value) {
            return Err(format!("'{a}' follows from the {by}: it cannot follow '{b}' by a link"));
        }
        let r = combine(&[(l.coefficient, &row(&rows, l.follows))]);
        if !r.is_empty() {
            rows.insert(l.value, r);
            order.push(l.value);
        }
    }
    // The free parameters: one per column that moves some field, named after the field its codeword was written
    // on, which comes first in its list (its value and esd are reported).
    let mut cols: BTreeMap<usize, Vec<(Field, f64)>> = BTreeMap::new();
    let mut seen = std::collections::BTreeSet::new();
    for f in &order {
        if !seen.insert(*f) {
            continue;
        }
        let Some(r) = rows.get(f) else { continue };
        for (&col, &coef) in r {
            cols.entry(col).or_default().push((*f, coef));
        }
    }
    let mut free = Vec::new();
    let mut col_name: BTreeMap<usize, String> = BTreeMap::new();
    for (col, mut fields) in cols {
        if let Some(k) = fields.iter().position(|q| q.0 == namer[col]) {
            fields.swap(0, k);
        }
        col_name.insert(col, address(fields[0].0, m));
        free.push(Free { name: address(fields[0].0, m), fields, damping: 1.0, value: None });
    }
    Ok((free, col_name))
}

pub(crate) enum CellRule {
    Free,
    Fixed,
    Equal(usize),
}

/// The cell parameters the crystal system ties or fixes (a, b, c, α, β, γ).
pub(crate) fn cell_rules(ph: &crate::powder::Phase) -> [CellRule; 6] {
    use CellRule::*;
    let c = &ph.cell;
    match ph.group.crystal_system {
        "cubic" => [Free, Equal(0), Equal(0), Fixed, Fixed, Fixed],
        "tetragonal" => [Free, Equal(0), Free, Fixed, Fixed, Fixed],
        "hexagonal" => [Free, Equal(0), Free, Fixed, Fixed, Fixed],
        "trigonal" if (c.gamma - 120.0).abs() < 1e-6 => [Free, Equal(0), Free, Fixed, Fixed, Fixed],
        "trigonal" => [Free, Equal(0), Equal(0), Free, Equal(3), Equal(3)],
        "orthorhombic" => [Free, Free, Free, Fixed, Fixed, Fixed],
        "monoclinic" => {
            // The unique axis: the axis of the 2-fold rotation (or the normal of the mirror).
            let unique = ph
                .group
                .ops
                .iter()
                .find_map(|op| {
                    let tr = op.r[0][0] + op.r[1][1] + op.r[2][2];
                    (tr == -1 || tr == 1).then(|| (0..3).find(|&i| op.r[i][i] == if tr == -1 { 1 } else { -1 })).flatten()
                })
                .unwrap_or(1);
            let mut r = [Free, Free, Free, Fixed, Fixed, Fixed];
            r[3 + unique] = Free;
            r
        }
        _ => [Free, Free, Free, Free, Free, Free],
    }
}

#[cfg(test)]
mod limit_state_tests {
    use super::LimitState;
    /// L145: inside, on a limit, outside.
    #[test]
    fn where_a_value_stands() {
        assert_eq!(LimitState::of(0.5, (Some(0.0), Some(1.0))), LimitState::Inside);
        assert_eq!(LimitState::of(0.0, (Some(0.0), None)), LimitState::AtLimit);
        assert_eq!(LimitState::of(-0.1, (Some(0.0), None)), LimitState::Outside);
        assert_eq!(LimitState::of(2.0, (None, None)).words(), None);
        assert_eq!(LimitState::of(1.0, (None, Some(1.0))).words(), Some("at its limit"));
    }
}

#[cfg(test)]
mod tests {
    use super::codeword;

    #[test]
    fn codewords_decode_as_the_manual_says() {
        // C = sign(a)(10p + |a|), manual p. 13.
        assert_eq!(codeword(41.0), Some((4, 1.0)));
        assert_eq!(codeword(-20.5).map(|(p, a)| (p, (a * 1e6).round() / 1e6)), Some((2, -0.5)));
        assert_eq!(codeword(161.0), Some((16, 1.0)));
        // Aut = 1: a bare multiplier is a parameter numbered automatically.
        assert_eq!(codeword(1.0), Some((0, 1.0)));
        // Fixed: no code, or a parameter number with a zero multiplier.
        assert_eq!(codeword(0.0), None);
        assert_eq!(codeword(10.0), None);
    }
}
