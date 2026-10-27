//! Guided refinement (P55, L302): the documented order of refinement offered as a plan of steps (L110's
//! `[[refine.step]]`), each step naming the values it frees.
//!
//! The order is FullProf's manual's "sensible sequence of refinement of a crystal structure" (FullProf_Manual.pdf,
//! "The Rietveld Method in Practice", p. 14), which refers the reader to McCusker et al. (1999, J. Appl. Cryst. 32,
//! 36–50, in `manuals/papers/`) to complete it; from McCusker et al., §7 ("Refining the structure", pp. 45–46): the
//! positions of the heavier atoms before those of the lighter, refinement "in sets of two to five cycles", and at the
//! end all the parameters together, to a largest shift/e.s.d. of 0.10, "to obtain correct estimated standard
//! deviations". The steps are cumulative, as the manual's "add" says:
//!
//! 1. the scale factor;
//! 2. the zero point, the first background term and the cell (the manual: "or better fix the background using linear
//!    interpolation between a set of fixed points": background points are not refined);
//! 3. the atomic positions, the heavier atoms' first (McCusker), then all of them with the overall B; a phase described
//!    by symmetry modes frees its mode amplitudes;
//! 4. the peak shape and the asymmetry (U, V, W, X, Y, η; Bérar–Baldinozzi P1–P4, or S/L and H/L for the
//!    Thompson–Cox–Hastings profile); over an instrumental resolution file the phase's U, X, Y are its microstructure
//!    and wait for step 10;
//! 5. the occupancies, "if required": the atoms of a partly occupied site (one shared by several atoms, or whose
//!    Occ/m is below that of the phase's fullest site);
//! 6. the overall B turned into each atom's isotropic B;
//! 7. the rest of the background series;
//! 8. the anisotropic betas, of the atoms that have them;
//! 9. Sycos and Sysin (constant wavelength);
//! 10. the microstructure: size and strain (the Gaussian size; the quartic strain and its share; the anisotropic size;
//!     over a resolution file, U, X and Y);
//! 11. all of them together, to a largest shift/e.s.d. of 0.10 (McCusker).
//!
//! A step that would add nothing to the one before is left out. "Heavier" is the scattering power: the X-ray form
//! factor at s = 0 (Σaᵢ + c), or |b| for neutrons; an atom is heavier when its power is at least half the phase's
//! largest, and the split is made only when the phase has atoms on both sides. What the manual's list leaves out
//! (preferred orientation, the wavelength, the micro-absorption) is refined only in the last step, and only when the
//! refinement it is added to already refines it.

use crate::jacobian::Field;
use crate::powder::{BackgroundModel, InstrumentalWidths, Intensities, Model, Shape};
use crate::refinable::{address, Refinement};
use crate::scattering::Scatterer;
use crate::steering::PlanStep;

/// One step of the documented order: the plan step (its name, the addresses it frees, its cycles and convergence), the
/// source of the step, and what it adds to and takes from the step before.
#[derive(Clone, Debug, PartialEq)]
pub struct GuidedStep {
    pub step: PlanStep,
    pub source: &'static str,
    pub adds: Vec<String>,
    pub drops: Vec<String>,
}

const MANUAL: &str = "FullProf manual, p. 14";
const MCCUSKER: &str = "McCusker et al. 1999, §7";

fn power(sc: &Scatterer) -> f64 {
    match sc {
        Scatterer::Xray { a, c, .. } => a.iter().sum::<f64>() + c,
        Scatterer::Neutron { b, .. } => b.abs(),
        Scatterer::XrayTable { points, .. } => points.first().map_or(0.0, |p| p.1),
    }
}

/// The documented order for a model (`Model::fields`, the values the symmetry leaves free).
pub fn documented_plan(m: &Model) -> Vec<GuidedStep> {
    let free: Vec<Field> = m.fields().into_iter().filter(|f| f.fixed_by.is_none()).map(|f| f.field).collect();
    let pick = |keep: &dyn Fn(Field) -> bool| -> Vec<Field> { free.iter().copied().filter(|&f| keep(f)).collect() };
    let structure = |p: usize| matches!(m.phases[p].intensities, Intensities::Structure(..)) && !m.phases[p].extract;
    let series = !matches!(m.background, BackgroundModel::Points(_) | BackgroundModel::Spline(_));
    let resolution = !matches!(m.instrumental, InstrumentalWidths::None);
    let tch = |p: usize| m.phases[p].shape == Shape::ThompsonCoxHastings;

    // Heavier atoms: per phase, those with at least half its largest scattering power, when lighter ones exist too.
    let heavier = |p: usize, a: usize| -> bool {
        let Intensities::Structure(_, sc) = &m.phases[p].intensities else { return true };
        let top = sc.iter().map(power).fold(0.0, f64::max);
        let light = sc.iter().any(|s| power(s) < 0.5 * top);
        !light || sc.get(a).is_some_and(|s| power(s) >= 0.5 * top)
    };
    // A partly occupied site, whatever the file's normalisation of Occ (FullProf's Occ is the chemical occupancy times
    // m/M, but a file may scale it, as Pd3MnD's does to m; the scale factor takes the rest): an atom that shares its
    // position with another (a mixed site), or whose Occ/m is below that of the phase's fullest site.
    let partial = |p: usize, a: usize| -> bool {
        let ph = &m.phases[p];
        let Intensities::Structure(st, _) = &ph.intensities else { return false };
        let per_site = |x: [f64; 3], occ: f64| occ / crate::site::Site::new(&ph.group, x).map_or(ph.group.multiplicity(), |s| s.multiplicity) as f64;
        let fullest = st.atoms.iter().map(|at| per_site(at.x, at.occ)).fold(0.0, f64::max);
        let at = &st.atoms[a];
        let shared = st.atoms.iter().enumerate().any(|(b, o)| b != a && (0..3).all(|k| (o.x[k] - at.x[k]).abs() < 1e-4));
        shared || per_site(at.x, at.occ) < 0.995 * fullest
    };
    let modes = |p: usize| m.phases[p].modes.is_some();

    let groups: Vec<(&str, &'static str, Vec<Field>, Vec<Field>)> = vec![
        ("scale factor", MANUAL, pick(&|f| matches!(f, Field::Scale(p) if structure(p))), vec![]),
        ("zero point, first background term and cell", MANUAL, pick(&|f| matches!(f, Field::Zero | Field::Cell(..)) || (series && f == Field::Background(0))), vec![]),
        ("positions of the heavier atoms", MCCUSKER, pick(&|f| matches!(f, Field::AtomX(p, a, _) if structure(p) && !modes(p) && heavier(p, a)) || matches!(f, Field::Mode(p, _) if structure(p))), vec![]),
        ("all positions and the overall B", MANUAL, pick(&|f| matches!(f, Field::AtomX(p, _, _) if structure(p) && !modes(p)) || matches!(f, Field::Bov(p) if structure(p))), vec![]),
        (
            "peak shape and asymmetry",
            MANUAL,
            pick(&|f| match f {
                Field::U(_) | Field::V(_) | Field::W(_) | Field::X(_) | Field::Y(_) => !resolution,
                Field::Eta0(_) => true,
                Field::SL(p) | Field::DL(p) => tch(p),
                Field::Asym(p, _) => !tch(p) || m.phases[p].asym.iter().any(|&v| v != 0.0),
                _ => false,
            }),
            vec![],
        ),
        ("occupancies of the partly occupied sites", MANUAL, pick(&|f| matches!(f, Field::AtomOcc(p, a) if structure(p) && partial(p, a))), vec![]),
        ("each atom's isotropic B for the overall B", MANUAL, pick(&|f| matches!(f, Field::AtomB(p, _) if structure(p))), pick(&|f| matches!(f, Field::Bov(_)))),
        ("the rest of the background", MANUAL, pick(&|f| series && matches!(f, Field::Background(_))), vec![]),
        ("anisotropic betas", MANUAL, pick(&|f| matches!(f, Field::AtomBeta(p, _, _) if structure(p))), vec![]),
        ("Sycos and Sysin", MANUAL, pick(&|f| matches!(f, Field::Sycos | Field::Sysin)), vec![]),
        (
            "microstructure: size and strain",
            MANUAL,
            pick(&|f| match f {
                Field::GaussSize(_) | Field::Quartic(..) | Field::QuarticShare(_) | Field::ShapeSize(_) => true,
                Field::U(_) | Field::X(_) | Field::Y(_) => resolution,
                _ => false,
            }),
            vec![],
        ),
    ];
    let mut out = Vec::new();
    let mut current: Vec<Field> = Vec::new();
    for (name, source, add, drop) in groups {
        let adds: Vec<Field> = add.into_iter().filter(|f| !current.contains(f)).collect();
        let drops: Vec<Field> = drop.into_iter().filter(|f| current.contains(f)).collect();
        if adds.is_empty() && drops.is_empty() {
            continue;
        }
        current.retain(|f| !drops.contains(f));
        current.extend(adds.iter().copied());
        out.push(GuidedStep {
            step: PlanStep { name: format!("{} {name}", out.len() + 1), free: current.iter().map(|&f| address(f, m)).collect(), damping: None, cycles: Some(5), eps: None },
            source,
            adds: adds.iter().map(|&f| address(f, m)).collect(),
            drops: drops.iter().map(|&f| address(f, m)).collect(),
        });
    }
    if !current.is_empty() {
        out.push(GuidedStep {
            step: PlanStep { name: format!("{} all together, to convergence", out.len() + 1), free: current.iter().map(|&f| address(f, m)).collect(), damping: None, cycles: None, eps: Some(0.1) },
            source: MCCUSKER,
            adds: Vec::new(),
            drops: Vec::new(),
        });
    }
    out
}

/// A refinement guided by the documented order: the parameters of `base` (its ties, constraints and steering kept) and
/// every value a step frees that `base` does not refine, each a parameter of its own; the plan its steps. The last
/// step frees every parameter of the refinement, those `base` refines outside the documented order too (preferred
/// orientation, the wavelength …), but what the plan drops and `base` does not refine: the overall B, kept at its
/// value once the atoms' B's are free (the model's B is the overall B plus the atom's, so this is the manual's turning
/// of the one into the others; refining both would be degenerate).
pub fn guided(m: &Model, base: &Refinement) -> Refinement {
    let steps = documented_plan(m);
    let mut r = base.clone();
    let have: Vec<Field> = r.parameters.iter().flatten().map(|t| t.0).collect();
    let mut dropped: Vec<String> = Vec::new();
    for s in &steps {
        for a in &s.step.free {
            if let Ok(f) = crate::refinable::field_of(a, m) {
                if !r.parameters.iter().flatten().any(|t| t.0 == f) {
                    r.parameters.push(vec![(f, 1.0)]);
                }
            }
        }
        dropped.extend(s.drops.iter().cloned());
    }
    let all: Vec<String> = r
        .parameters
        .iter()
        .filter(|t| have.contains(&t[0].0) || !dropped.contains(&address(t[0].0, m)))
        .map(|t| address(t[0].0, m))
        .collect();
    r.plan = steps.into_iter().map(|s| s.step).collect();
    if let Some(last) = r.plan.last_mut() {
        last.free = all;
    }
    r
}
