//! A state's derived quantities in one bundle (Explore E3, L115): what a module showing a refinement's state, at any
//! cycle (`lsq::Trajectory::state_at`), needs in one call, from the calls that compute each (`Model::gaussian_h2`,
//! `diagnostics::reflections`, `geometry::contacts_with_su`, `bonds::valence_sums`, `Modes::irreps`,
//! `Model::breadths`, `lsq::restraint_states`), with the plausible ranges of the bonds (`plausible`) from the
//! bond-valence and ionic-radius tables.

use crate::bonds::{self, Contact, Instability, ValenceSum};
use crate::cell::Cell;
use crate::data::Pattern;
use crate::diagnostics::ReflectionRow;
use crate::geometry::Propagator;
use crate::lsq::{RestraintState, Statistics};
use crate::microstructure::Breadths;
use crate::powder::{Intensities, Model};
use crate::refinable::ParameterSet;

/// The peak width of a phase at one 2θ: FWHM and η as the profile uses them, and the Gaussian H_G² before it is
/// clamped at 0 (negative where U, V, W give a width no peak can have) with the Lorentzian H_L.
#[derive(Clone, Debug, PartialEq)]
pub struct Width {
    pub two_theta: f64,
    pub fwhm: f64,
    pub eta: f64,
    pub h_g2: f64,
    pub h_l: f64,
}

/// A mode's or an irrep's amplitude with its su (Å).
pub type Amplitude = (String, f64, f64);

/// One phase's part of the state.
#[derive(Clone, Debug, PartialEq)]
pub struct PhaseState {
    pub name: String,
    pub cell: Cell,
    pub volume: f64,
    /// a, b, c (Å), α, β, γ (degrees) and V (Å³), each with its su propagated from the covariance given (0 without
    /// one), the cell's ties counted (a cubic cell's V has 3a²σ(a)): Explore's E18a.
    pub cell_su: [(f64, f64); 7],
    /// Against the reference cell given (`Options::reference_cells`): a/a₀ − 1, b/b₀ − 1, c/c₀ − 1 and V/V₀ − 1.
    pub strain: Option<[f64; 4]>,
    pub widths: Vec<Width>,
    /// The atoms' labels (unique), and each one's contacts within `Options::dmax` (atom, operator, translation, with
    /// propagated su's when a covariance is given), as `geometry::contacts_with_su` lists them.
    pub labels: Vec<String>,
    pub contacts: Vec<Vec<Contact>>,
    pub angles: Vec<bonds::Angle>,
    /// The bond-valence sums, when every atom has a charge (`Atom::charge`, or a species such as "Y+3").
    pub valence: Option<(Vec<ValenceSum>, Instability)>,
    pub modes: Vec<Amplitude>,
    pub irreps: Vec<Amplitude>,
    /// The size and strain breadths of each reflection (TCH profiles only).
    pub breadths: Vec<Breadths>,
}

/// The state: the agreement, the reflections, each phase's part, the restraints.
#[derive(Clone, Debug, PartialEq)]
pub struct State {
    pub stats: Statistics,
    pub reflections: Vec<ReflectionRow>,
    pub phases: Vec<PhaseState>,
    pub restraints: Vec<RestraintState>,
}

/// What `summary` computes beyond the defaults.
/// The 2θ range and step at which the widths are sampled over the data's range [lo, hi]: 120 steps, at least 0.1°
/// (Explore's peak-shape view; L145).
pub fn width_sampling(lo: f64, hi: f64) -> (f64, f64, f64) {
    (lo, hi, ((hi - lo) / 120.0).max(0.1))
}

/// The 2θ range over which the Gaussian width's square H_G² = U tan²θ + V tanθ + W is negative at some state: from
/// sampled points (2θ, H_G²) of every state (Explore's widths view, hatched there; L148). None where it never is.
pub fn negative_gaussian_range(points: impl IntoIterator<Item = (f64, f64)>) -> Option<(f64, f64)> {
    points.into_iter().filter(|&(_, hg2)| hg2 < 0.0).fold(None, |acc: Option<(f64, f64)>, (t, _)| Some(acc.map_or((t, t), |(a, b)| (a.min(t), b.max(t)))))
}

#[cfg(test)]
mod sampling_tests {
    /// L145: 120 steps over the range, at least 0.1°.
    #[test]
    fn the_widths_sampling() {
        assert_eq!(super::width_sampling(10.0, 130.0), (10.0, 130.0, 1.0));
        assert_eq!(super::width_sampling(10.0, 16.0), (10.0, 16.0, 0.1));
        assert_eq!(super::negative_gaussian_range([(10.0, 1.0), (12.0, -0.1), (11.0, -0.2), (15.0, 2.0)]), Some((11.0, 12.0)));
        assert_eq!(super::negative_gaussian_range([(10.0, 1.0)]), None);
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Options {
    /// The 2θ range and step of the widths (default: the data's range, 1°).
    pub widths: Option<(f64, f64, f64)>,
    /// Contacts up to this distance (Å; 3.0), angles between neighbours up to `dangl` (0: none).
    pub dmax: f64,
    pub dangl: f64,
    /// A reference cell per phase, for the strain.
    pub reference_cells: Vec<Option<Cell>>,
}

impl Default for Options {
    fn default() -> Options {
        Options { widths: None, dmax: 3.0, dangl: 0.0, reference_cells: Vec::new() }
    }
}

/// The state of `model` against `data` (E3). With the parameters and the covariance of a refinement
/// (`Refinement::covariance`), the su's of the distances, angles, sums and amplitudes are propagated (`geometry`);
/// without, they are 0.
pub fn summary(model: &Model, data: &Pattern, set: Option<&ParameterSet>, refinement: Option<&crate::lsq::Refinement>, opts: &Options) -> State {
    summary_at(model, data, set, refinement.map(|r| r.covariance.as_slice()), opts)
}

/// `summary` with the covariance of the free parameters given directly (E14a): a cycle's (`Trajectory::cycles[k −
/// 1].covariance`, at that cycle's values), so that the su's of the modes, irreps, distances and angles follow the
/// path. An empty covariance (a cycle read back from a results file, which keeps none) gives no su's, as None.
pub fn summary_at(model: &Model, data: &Pattern, set: Option<&ParameterSet>, covariance: Option<&[Vec<f64>]>, opts: &Options) -> State {
    let covariance = covariance.filter(|c| !c.is_empty());
    let calc = model.calculate(&data.x);
    let np = set.map_or(0, |s| s.free.len());
    let stats = crate::lsq::statistics(model, data, &calc.ycalc, np);
    let reflections = crate::diagnostics::reflections(model, data);
    let prop = match (set, covariance) {
        (Some(s), Some(c)) => Some(Propagator::from_covariance(model, s, c)),
        _ => None,
    };
    let (lo, hi) = (data.x.first().copied().unwrap_or(0.0), data.x.last().copied().unwrap_or(0.0));
    let (w0, w1, step) = opts.widths.unwrap_or((lo, hi, 1.0));
    let mut phases = Vec::new();
    for (p, ph) in model.phases.iter().enumerate() {
        let volume = ph.cell.volume();
        let cell_of = |m: &Model| -> Vec<f64> {
            let c = m.phases[p].cell;
            vec![c.a, c.b, c.c, c.alpha, c.beta, c.gamma, c.volume()]
        };
        let cs = match &prop {
            Some(pr) => pr.of_phase(p, &cell_of),
            None => cell_of(model).into_iter().map(|v| (v, 0.0)).collect(),
        };
        let cell_su: [(f64, f64); 7] = std::array::from_fn(|k| cs[k]);
        let strain = opts.reference_cells.get(p).copied().flatten().map(|r| [ph.cell.a / r.a - 1.0, ph.cell.b / r.b - 1.0, ph.cell.c / r.c - 1.0, volume / r.volume() - 1.0]);
        let mut widths = Vec::new();
        let mut t = w0;
        while t <= w1 + 1e-9 && step > 0.0 {
            let (fwhm, eta) = model.width(ph, t);
            let (_, h_l) = model.width_components(ph, t);
            widths.push(Width { two_theta: t, fwhm, eta, h_g2: model.gaussian_h2(ph, t, 1), h_l });
            t += step;
        }
        let (mut labels, mut contacts, mut angles, mut valence) = (Vec::new(), Vec::new(), Vec::new(), None);
        if matches!(ph.intensities, Intensities::Structure(..)) {
            if let Ok((st, cs, an)) = crate::geometry::contacts_with_su(model, p, prop.as_ref(), opts.dmax, opts.dangl) {
                labels = crate::geometry::unique_labels(&st);
                let mut cs = cs;
                if st.atoms.iter().all(|a| bonds::atom_ion(a).is_some()) {
                    let range = ph.valence_sum_range.unwrap_or(bonds::VALENCE_RANGE);
                    valence = bonds::chemical_occupancies(&st).ok().and_then(|occ| bonds::valence_sums(&st, &cs, &occ, 0.2).ok()).map(|(v, i)| (v.into_iter().map(|s| s.with_range(range)).collect::<Vec<_>>(), i));
                    set_plausible(&st, &mut cs, valence.as_ref().map(|v| v.0.as_slice()));
                }
                contacts = cs;
                angles = an;
            }
        }
        let (mut modes, mut irreps) = (Vec::new(), Vec::new());
        if let Some(ms) = &ph.modes {
            let cov = |a: usize, b: usize| match (set, covariance) {
                (Some(s), Some(c)) => s.field_covariance(c, crate::jacobian::Field::Mode(p, a), crate::jacobian::Field::Mode(p, b)),
                _ => 0.0,
            };
            modes = ms.modes.iter().enumerate().map(|(k, m)| (m.name.clone(), m.amplitude, cov(k, k).max(0.0).sqrt())).collect();
            irreps = ms.irreps().iter().map(|ir| (ir.name.clone(), ir.amplitude, ms.irrep_esd(ir, cov))).collect();
        }
        let breadths = reflections.iter().filter(|r| r.phase == p).filter_map(|r| model.breadths(p, r.hkl, r.wavelength)).collect();
        phases.push(PhaseState { name: crate::cif_write::phase_name(model, p), cell: ph.cell, volume, cell_su, strain, widths, labels, contacts, angles, valence, modes, irreps, breadths });
    }
    let restraints = set.map(|s| crate::lsq::restraint_states(model, s)).unwrap_or_default();
    State { stats, reflections, phases, restraints }
}

/// Each cation–anion contact's plausible distances (E18b): `plausible` for the pair, the cation's ion and the anion's,
/// at the cation's coordination at this state (its bonds in the valence sum when there is one, else its contacts to
/// atoms of the other sign).
fn set_plausible(st: &crate::structure_factor::Structure, contacts: &mut [Vec<Contact>], sums: Option<&[ValenceSum]>) {
    let ions: Vec<Option<(String, i32)>> = st.atoms.iter().map(bonds::atom_ion).collect();
    let coordination = |cation: usize, contacts: &[Vec<Contact>]| -> usize {
        sums.and_then(|s| s.iter().find(|v| v.atom == cation)).map(|v| v.coordination).unwrap_or_else(|| {
            let q = ions[cation].as_ref().map_or(0, |i| i.1);
            contacts[cation].iter().filter(|c| ions[c.atom].as_ref().is_some_and(|i| i.1.signum() == -q.signum())).count()
        })
    };
    let n: Vec<usize> = (0..contacts.len()).map(|i| coordination(i, contacts)).collect();
    for (i, list) in contacts.iter_mut().enumerate() {
        for c in list.iter_mut() {
            let (Some(a), Some(b)) = (&ions[i], &ions[c.atom]) else { continue };
            if a.1.signum() * b.1.signum() >= 0 {
                continue;
            }
            let (cat, an, k) = if a.1 > 0 { (a, b, i) } else { (b, a, c.atom) };
            c.plausible = Some(plausible((&cat.0, cat.1), (&an.0, an.1), &[n[k]]));
        }
    }
}

/// One reflection's profile over its window (E3): the points of the data inside the window of the phase's reflection
/// hkl (wavelength 1), the peak's own contribution at each (I·Ω, scale and corrections included) and the background
/// there, apart. None when the model has no such reflection in the data's range.
pub fn reflection_profile(model: &Model, data: &Pattern, phase: usize, hkl: [i32; 3]) -> Option<(Vec<f64>, Vec<f64>, Vec<f64>)> {
    let calc = model.calculate(&data.x);
    let pk = calc.peaks.iter().find(|k| k.phase == phase && k.hkl == hkl && k.wavelength == 1)?;
    let excluded = model.excluded_points(&data.x);
    let mut peak = vec![0.0; data.x.len()];
    model.peak_profile(&data.x, &excluded, pk, |i, v| peak[i] += v);
    let idx: Vec<usize> = (0..data.x.len()).filter(|&i| (data.x[i] - pk.two_theta).abs() <= pk.reach).collect();
    Some((idx.iter().map(|&i| data.x[i]).collect(), idx.iter().map(|&i| peak[i]).collect(), idx.iter().map(|&i| calc.background[i]).collect()))
}

/// The plausible distances of a cation–anion pair (E3): for each coordination number of `coordinations`, the
/// bond-valence single-bond distance d₀ − B₀ ln(V/N) (V the cation's valence; Brown 2002, *The chemical bond in
/// inorganic chemistry*, eq. 3.1 with s = V/N), and the sum of the ionic radii (`bonds::ionic_radius`, CrysFML's
/// table). None for the bond valence when the table has no d₀ for the pair.
pub fn plausible(cation: (&str, i32), anion: (&str, i32), coordinations: &[usize]) -> bonds::Plausible {
    let name = |el: &str, q: i32| format!("{}{}{}", el.to_uppercase(), if q < 0 { "-" } else { "+" }, q.abs());
    let pars = bonds::parameters(&name(cation.0, cation.1), &name(anion.0, anion.1));
    let v = cation.1.abs() as f64;
    let bv = coordinations.iter().map(|&n| (n, pars.map(|(d0, b0, _)| d0 - b0 * (v / n as f64).ln()))).collect();
    let ionic = bonds::ionic_radius(&cation.0.to_uppercase(), cation.1) + bonds::ionic_radius(&anion.0.to_uppercase(), anion.1);
    (bv, ionic)
}
