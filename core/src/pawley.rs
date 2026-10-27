//! The Pawley fit (P36; Pawley 1981, J. Appl. Cryst. 14, 357–361): a phase without a structure whose reflection
//! intensities are least-squares parameters, refined with the cell, the profile and the background. FullProf has no
//! Pawley fit (its manual names the method beside Le Bail's, p. 14); GSAS-II has one (`PWLref`, its reflections' F²),
//! and the core follows it in what it ties:
//!
//! - A phase is a Pawley phase when it has a list of reflections (`Intensities::List`) and is not a Le Bail phase
//!   (`extract` false): its intensities are then fields of the model (`Field::Intensity`, addressed
//!   `phase1.reflection_intensity.1_-1_0`) that a refinement frees. The list's intensity carries no scale factor, so
//!   the phase's scale is not refined with them.
//! - Reflections that coincide for every cell the crystal system allows (330 and 411 of a cubic phase, h k 0 and k h 0
//!   of 4/m) cannot be told apart by the pattern: they are tied, with equal |F|², so that their intensities stay in
//!   proportion to their multiplicities (GSAS-II's `GetPawleyConstr` equivalences their F², by the same sums of
//!   indices). Whether two reflections coincide is found without the crystal system's name: the metric averaged over
//!   the group's rotations, R S Rᵀ for a generic S, is a metric of the lattice the group allows, and two reflections
//!   coincide for every such metric when they have the same 1/d² for two of them. GSAS-II also ties, in −1, 2/m and
//!   mmm, reflections whose d differ by less than 10⁻⁴ (accidental overlaps); the core does not: they part as the
//!   cell is refined, and the engine leaves out the directions the data do not determine.
//! - An intensity is held at zero where the data push it below (`refinable::limits`); GSAS-II penalises negative ones
//!   instead ('Pawley neg wt').
//! - The intensities start from Le Bail's partition (`Model::le_bail`, a few iterations at the starting cell and
//!   profile), the tied ones then shared out by multiplicity.

use crate::data::Pattern;
use crate::jacobian::Field;
use crate::powder::{Intensities, Model};
use crate::refinable::Refinement;

/// The Le Bail iterations that give a Pawley fit its starting intensities, when the caller states none.
pub const START_ITERATIONS: usize = 10;

/// The metric of the lattice the group allows, averaged from `s` over the group's rotations (reflections transform as
/// h·R, so the average is Σ R s Rᵀ).
fn averaged_metric(m: &Model, p: usize, s: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let ops = &m.phases[p].group.ops;
    let mut g = [[0.0; 3]; 3];
    for op in ops {
        let r = op.r.map(|row| row.map(|v| v as f64));
        for i in 0..3 {
            for j in 0..3 {
                g[i][j] += (0..3).map(|k| (0..3).map(|l| r[i][k] * s[k][l] * r[j][l]).sum::<f64>()).sum::<f64>() / ops.len() as f64;
            }
        }
    }
    g
}

fn quad(g: &[[f64; 3]; 3], h: [i32; 3]) -> f64 {
    let h = h.map(|v| v as f64);
    (0..3).map(|i| (0..3).map(|j| h[i] * g[i][j] * h[j]).sum::<f64>()).sum()
}

/// The reflections of phase `p`'s list that coincide for every cell its group allows, as groups of indices into the
/// list (each of two or more, in the list's order); empty for a phase with a structure.
pub fn coincident(m: &Model, p: usize) -> Vec<Vec<usize>> {
    let Intensities::List(list) = &m.phases[p].intensities else { return Vec::new() };
    // Two generic symmetric matrices (fixed, so that the groups are the same at every call).
    let s1 = [[1.0, 0.173, 0.291], [0.173, 1.37, 0.113], [0.291, 0.113, 0.83]];
    let s2 = [[0.71, -0.237, 0.059], [-0.237, 1.13, 0.311], [0.059, 0.311, 1.51]];
    let (g1, g2) = (averaged_metric(m, p, s1), averaged_metric(m, p, s2));
    let key: Vec<(f64, f64)> = list.iter().map(|r| (quad(&g1, r.0), quad(&g2, r.0))).collect();
    let same = |a: (f64, f64), b: (f64, f64)| (a.0 - b.0).abs() <= 1e-9 * a.0.abs().max(b.0.abs()) && (a.1 - b.1).abs() <= 1e-9 * a.1.abs().max(b.1.abs());
    let mut seen = vec![false; list.len()];
    let mut out = Vec::new();
    for i in 0..list.len() {
        if seen[i] {
            continue;
        }
        let group: Vec<usize> = (i..list.len()).filter(|&j| !seen[j] && same(key[i], key[j])).collect();
        for &j in &group {
            seen[j] = true;
        }
        if group.len() > 1 {
            out.push(group);
        }
    }
    out
}

/// The parameters of phase `p`'s Pawley fit: one per reflection of the list whose peak's window reaches the data,
/// the coincident ones (`coincident`) one parameter whose intensities move in proportion to their multiplicities
/// (the multiplier of each, m_k/m_first).
pub fn parameters(m: &Model, p: usize, data: &Pattern) -> Vec<Vec<(Field, f64)>> {
    let Intensities::List(list) = &m.phases[p].intensities else { return Vec::new() };
    let (lo, hi) = (data.x.first().copied().unwrap_or(0.0), data.x.last().copied().unwrap_or(0.0));
    // A reflection whose window reaches the data, its centre within or not.
    let seen: std::collections::BTreeSet<[i32; 3]> = m.peaks(lo - 10.0, hi + 10.0).into_iter().filter(|pk| pk.phase == p && pk.two_theta + pk.reach >= lo && pk.two_theta - pk.reach <= hi).map(|pk| pk.hkl).collect();
    let groups = coincident(m, p);
    let mut taken = vec![false; list.len()];
    let mut out = Vec::new();
    for k in 0..list.len() {
        if taken[k] || !seen.contains(&list[k].0) {
            continue;
        }
        let group = groups.iter().find(|g| g.contains(&k)).cloned().unwrap_or_else(|| vec![k]);
        for &j in &group {
            taken[j] = true;
        }
        let first = list[group[0]].1 as f64;
        out.push(group.iter().map(|&j| (Field::Intensity(p, j), list[j].1 as f64 / first)).collect());
    }
    out
}

/// What `setup` did.
#[derive(Clone, Debug, PartialEq)]
pub struct Setup {
    /// The phase (0, 1, …), the reflections of its list, and the intensities refined (parameters).
    pub phase: usize,
    pub reflections: usize,
    pub parameters: usize,
    /// The groups of coincident reflections tied, by their indices h k l.
    pub ties: Vec<Vec<[i32; 3]>>,
}

/// Make phase `p` a Pawley phase and `refinement` its fit: a phase with a structure takes the list of its reflections
/// to 5° past the data (`start::le_bail_list`), a Le Bail phase keeps its list; the intensities start from
/// `iterations` Le Bail partitions at the model's cell and profile (from the structure's intensities, or the list's),
/// the coincident ones are shared out by multiplicity, and the refinement frees every intensity of a reflection
/// within the data (`parameters`) and no longer the phase's scale factor. The refinement's other parameters are kept.
pub fn setup(m: &mut Model, data: &Pattern, refinement: &mut Refinement, p: usize, iterations: usize) -> Result<Setup, String> {
    let ph = m.phases.get(p).ok_or(format!("the model has {} phases", m.phases.len()))?;
    if matches!(ph.intensities, Intensities::Structure(..)) {
        // The structure's intensities at the start (scale included): the peaks of the first wavelength, the list's
        // being the total of both wavelengths.
        let (lo, hi) = (data.x.first().copied().unwrap_or(0.0), data.x.last().copied().unwrap_or(180.0));
        let lambda_min = m.lambda[0].min(m.lambda[1]);
        let mut list = crate::start::le_bail_list(&ph.group, &ph.cell, lambda_min, hi);
        let two = m.ratio != 0.0 && (m.lambda[1] - m.lambda[0]).abs() > 1e-9;
        let peaks = m.peaks(lo - 5.0, hi + 5.0);
        for r in list.iter_mut() {
            r.2 = peaks.iter().find(|pk| pk.phase == p && pk.wavelength == 1 && pk.hkl == r.0).map_or(1.0, |pk| pk.intensity * if two { 1.0 + m.ratio } else { 1.0 });
        }
        m.phases[p].intensities = Intensities::List(list);
        m.phases[p].modes = None;
    }
    m.phases[p].extract = true;
    for _ in 0..iterations {
        m.le_bail(data);
    }
    m.phases[p].extract = false;
    let groups = coincident(m, p);
    if let Intensities::List(list) = &mut m.phases[p].intensities {
        for g in &groups {
            let (sum, mult): (f64, f64) = g.iter().fold((0.0, 0.0), |(s, n), &k| (s + list[k].2, n + list[k].1 as f64));
            for &k in g {
                list[k].2 = sum * list[k].1 as f64 / mult;
            }
        }
    }
    let params = parameters(m, p, data);
    refinement.parameters.retain(|t| !t.iter().any(|(f, _)| matches!(f, Field::Scale(q) | Field::Intensity(q, _) if *q == p)));
    refinement.parameters.extend(params.iter().cloned());
    // Reflections that overlap by accident at this cell (three of LaMnO₃'s near 2θ = 150°, 1/d² within 3·10⁻⁴) leave
    // combinations the data cannot fix: they are left out of each cycle's step rather than stop the run (the guard
    // `singular`, R24), as Pawley's slack constraints and GSAS-II's ties do in their ways.
    refinement.steering.guards.singular = false;
    let Intensities::List(list) = &m.phases[p].intensities else { unreachable!() };
    let ties = params.iter().filter(|t| t.len() > 1).map(|t| t.iter().map(|(f, _)| if let Field::Intensity(_, k) = f { list[*k].0 } else { [0; 3] }).collect()).collect();
    Ok(Setup { phase: p, reflections: list.len(), parameters: params.len(), ties })
}

/// A tie of a Pawley fit: the phase (0, 1, …) and its reflections (h k l) with their intensities' multipliers (the
/// first 1).
pub type Tie = (usize, Vec<([i32; 3], f64)>);

/// The ties of a Pawley fit as a refinement holds them, from its parameters that move several intensities.
pub fn ties(m: &Model, r: &Refinement) -> Vec<Tie> {
    r.parameters
        .iter()
        .filter(|t| t.len() > 1 && t.iter().all(|(f, _)| matches!(f, Field::Intensity(..))))
        .filter_map(|t| {
            let Field::Intensity(p, _) = t[0].0 else { return None };
            let Intensities::List(list) = &m.phases.get(p)?.intensities else { return None };
            let rows = t.iter().filter_map(|(f, c)| if let Field::Intensity(q, k) = f { (*q == p).then(|| (list.get(*k).map(|e| e.0), *c)) } else { None }).filter_map(|(h, c)| h.map(|h| (h, c))).collect();
            Some((p, rows))
        })
        .collect()
}

/// Whether phase `p` of a refinement is fitted by Pawley's method: some of its intensities are refined.
pub fn is_pawley(r: &Refinement, p: usize) -> bool {
    r.parameters.iter().flatten().any(|(f, _)| matches!(f, Field::Intensity(q, _) if *q == p))
}
