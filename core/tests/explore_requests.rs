//! The Explore stream's requests to the core (explore/REQUESTS.md on branch explore, E1–E7; L115), on case A of the
//! mode course (`manuals/fullprof/runs/modes/A_modes`, Sr2ScSbO6 in P2₁/n by seven irreps of Fm-3m):
//! - `e1_trajectory`: three cycles; the model at every cycle rebuilt from the start and the shifts (`state_at`) is the
//!   engine's to the bit; each cycle's covariance gives its esd's, its spectrum runs from 1 down, its time is
//!   measured; the spectrum at the start is `check::directions`'s (E5); the cycles written in the results file read
//!   back to the bit, and the states rebuilt from the control file written there (`Trajectory::from_end`) are the
//!   engine's; the record as JSON (`drs refine --trace`);
//! - `e2_chi2_grid`: χ² on a 5 × 5 grid of Γ₄⁺ and X₃⁺ about the converged values, ± 1 esd: every point once, centre
//!   and corners first, the centre χ² at the values, the minimum there, and the rest the normal matrix's quadratic
//!   form within 5 %; stopped after a batch, the rest NaN;
//! - `e3_state`: the bundle: unclamped H_G² (negative where W makes it so), strain against a reference cell, the irreps
//!   with su's, Sb1's six contacts with su's, the bond-valence distances of Sb⁵⁺–O and Sc³⁺–O against case A's, a
//!   reflection's profile (the profiles at a point add up to the phase's curve);
//! - `e4_scene`: the cell's 20 atoms, Cartesian through the cell's matrix, Sb1's octahedron (8 faces, its volume the
//!   regular one's within 1 %), the irreps' displacements adding up to each image's distance from its reference; the
//!   hull of a cube;
//! - `e6_e7_modes`: the irreps' displacements add up to the whole; the 4 amplitude sets of case A's structure under
//!   Fm-3m's operations (24 before P21's fix, L218), the Explore stream's mirror solution among them (by −x,−y,z), GM1+'s sign never reversed,
//!   and the nearest set to a reference;
//! - `e14a_esds_at_every_cycle`, `e14b_the_sign_reference`, `e14c_fullprof_reported_amplitudes` (L208): the modes'
//!   and irreps' esd's from a cycle's covariance; the control file's sign reference and, per state, the equivalent set
//!   nearest it with the operation named; FullProf's reported amplitudes read from its `.sum`.

use drs::powder::Model;
use drs::refinable::ParameterSet;
use std::path::{Path, PathBuf};

fn root(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)
}

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display())).iter().map(|&b| b as char).collect()
}

/// Case A (Sr2ScSbO6, P2₁/n, seven irreps from Fm-3m) at FullProf's values, with its parent from AMPLIMODES's
/// analysis 1 of the course.
fn case_a() -> (Model, drs::data::Pattern, ParameterSet, Vec<drs::modes::Analysis>) {
    let d = root("../manuals/fullprof/runs/modes/A_modes");
    let p = drs::pcr::read(&latin1(&d.join("a_modes.pcr"))).unwrap();
    let l = Model::from_pcr(&p, &d, "a_modes").unwrap();
    let mut m = l.model;
    let a = drs::modes::read_amplimodes(&latin1(&root("../manuals/bcs/mode-crystallography-manual/course-um6p-2024/04_tools_materials/0_1_P21nToCubicTransData.txt")));
    m.phases[0].modes.as_mut().unwrap().parent = a[0].mode_parent();
    let set = ParameterSet::from_pcr(&p, &m).unwrap();
    (m, l.data.unwrap(), set, a)
}

#[test]
fn e6_e7_modes() {
    let (m, _, _, analyses) = case_a();
    let a = &analyses[0];
    let ph = &m.phases[0];
    let ms = ph.modes.as_ref().unwrap();
    let drs::powder::Intensities::Structure(st, _) = &ph.intensities else { panic!() };
    let species: Vec<String> = st.atoms.iter().map(|x| x.species.clone()).collect();
    // E6: the irreps' displacements sum to the whole, and displacement_of the amplitudes is displacement().
    let mut sum = vec![[0.0; 3]; ms.reference.len()];
    for ir in ms.irreps() {
        for (s, u) in sum.iter_mut().zip(ms.irrep_displacement(&ir.name)) {
            for k in 0..3 {
                s[k] += u[k];
            }
        }
    }
    let q: Vec<f64> = ms.modes.iter().map(|x| x.amplitude).collect();
    let all = ms.displacement_of(&q);
    for (j, (s, a)) in sum.iter().zip(&all).enumerate() {
        let whole = ms.displacement(j);
        assert!((0..3).all(|k| (s[k] - whole[k]).abs() < 1e-15 && (a[k] - whole[k]).abs() < 1e-15));
    }
    // E7: the sets equivalent to FullProf's amplitudes.
    let metric = a.reference_cell().unwrap_or(ph.cell);
    let eq = ms.equivalent_amplitudes(&species, &ph.group, &metric, &q).unwrap();
    for e in &eq {
        println!("{:24} {:?}: {}", e.operation, e.translation, e.amplitudes.iter().map(|v| format!("{v:+.4}")).collect::<Vec<_>>().join(" "));
    }
    // 4 sets (P21, L218; 24 before the whole cell was checked, 20 of them other structures): the identity, FullProf's
    // own amplitudes, first.
    assert_eq!(eq.len(), 4);
    assert!(eq[0].operation == "x,y,z" && eq[0].translation == [0, 0, 0] && eq[0].amplitudes.iter().zip(&q).all(|(x, y)| (x - y).abs() < 1e-6), "{:?}", eq[0]);
    // The Explore stream's mirror solution (from +0.1 Å starts): Γ₄⁺ +0.8288 for FullProf's −0.8286, the X₅⁺ pair
    // reversed, the rest as FullProf's; the same structure, by the parent's 2-fold rotation about c.
    let name = |n: &str| ms.modes.iter().position(|x| x.name == n).unwrap();
    let mut mirror = q.clone();
    for n in ["A3_GM4+", "A10_X5+", "A11_X5+"] {
        mirror[name(n)] = -mirror[name(n)];
    }
    let same = ms.same_structure(&species, &ph.group, &metric, &q, &mirror, 1e-4).unwrap().expect("the mirror solution is the same structure");
    println!("the mirror solution: {} {:?}", same.operation, same.translation);
    // GM1+, the totally symmetric irrep, keeps its sign under every operation: its sign reversed is another structure.
    let mut other = q.clone();
    other[name("A1_GM1+")] = -other[name("A1_GM1+")];
    assert!(ms.same_structure(&species, &ph.group, &metric, &q, &other, 1e-4).unwrap().is_none());
    // A reference's signs as the convention (`nearest_to`): from the mirror solution, FullProf's amplitudes as the
    // reference give them back. (AMPLIMODES's own block for case A is in its P2₁/c setting, whose X2+ sign and GM5+
    // and X5+ bases differ from this phase's, `symmetry_modes`: a reference must be in the phase's basis, as the
    // analysis the control file's modes come from is.)
    let near = ms.nearest_to(&species, &ph.group, &metric, &mirror, &q).unwrap();
    println!("from the mirror, nearest FullProf's: {} {:?}", near.operation, near.translation);
    assert!(near.amplitudes.iter().zip(&q).all(|(a, b)| (a - b).abs() < 1e-6), "{near:?}");
}

/// The largest difference of two models' values, relative to the value (1 below 1).
fn largest_difference(a: &Model, b: &Model) -> (f64, String) {
    a.fields().iter().zip(b.fields()).map(|(x, y)| ((x.value - y.value).abs() / x.value.abs().max(1.0), x.address.clone())).fold((0.0, String::new()), |acc, v| if v.0 > acc.0 { v } else { acc })
}

#[test]
fn e1_trajectory() {
    let (m, data, set, _) = case_a();
    // Three cycles from the file's values with every cycle's model kept as the engine reached it.
    let mut s = drs::lsq::Session::new(m, data.clone(), set.clone(), drs::lsq::Options { max_cycles: 3, eps: 1e-6 });
    let mut models = vec![s.model.clone()];
    for _ in 0..3 {
        s.step(1, |_| true);
        models.push(s.model.clone());
    }
    let t = s.trajectory();
    assert_eq!(t.cycles.len(), 3);
    // (a) The state at every cycle, rebuilt from the start and the shifts, is the engine's.
    let mut worst: f64 = 0.0;
    assert!(models[3].fields().len() > 20 && largest_difference(&models[0], &models[3]).0 > 1e-6, "the cycles move the model");
    for (k, mk) in models.iter().enumerate() {
        let (d, at) = largest_difference(mk, &t.state_at(k));
        println!("cycle {k}: state rebuilt within {d:.1e} ({at})");
        worst = worst.max(d);
    }
    assert!(worst < 1e-12, "{worst}");
    // (b), (c) The covariance, spectrum and time of every cycle: each esd the root of its covariance's diagonal (with
    // the multiplier), the spectrum from 1 down, the time measured.
    for c in &t.cycles {
        assert_eq!(c.covariance.len(), set.free.len());
        for (j, (q, p)) in c.params.iter().zip(&set.free).enumerate() {
            assert!((p.esd_of(c.covariance[j][j].sqrt()) - q.esd).abs() <= 1e-12 * q.esd.max(1e-300), "{}: {} against {}", q.name, p.esd_of(c.covariance[j][j].sqrt()), q.esd);
        }
        assert!((c.spectrum[0] - 1.0).abs() < 1e-12 && c.spectrum.windows(2).all(|w| w[0] >= w[1]) && c.seconds > 0.0);
        let r = drs::lsq::correlations_of(&c.covariance);
        assert!((0..r.len()).all(|j| (r[j][j] - 1.0).abs() < 1e-12 || r[j][j] == 0.0));
    }
    // E5: the spectrum at the start, from `check::directions`, is the first cycle's.
    let dirs = drs::check::directions(&t.start, &data, &set);
    assert_eq!(dirs.len(), set.free.len());
    for (d, l) in dirs.iter().zip(&t.cycles[0].spectrum) {
        assert!((d.ratio - l).abs() <= 1e-9 * l.max(1e-12), "{} against {l}", d.ratio);
        assert!((d.components.iter().map(|c| c.1 * c.1).sum::<f64>() - 1.0).abs() < 1e-9);
    }
    println!("spectrum at the start: {:.3e} … {:.3e}", dirs[0].ratio, dirs.last().unwrap().ratio);
    // (d) Written in the results file and read back: the cycles as run, to the last bit; the states rebuilt from the
    // control file written (the end) and the cycles read are the engine's.
    let r = s.finish();
    let d = root("../manuals/fullprof/runs/modes/A_modes");
    let p = drs::pcr::read(&latin1(&d.join("a_modes.pcr"))).unwrap();
    let control = drs::control::Control::from_pcr(&p, &d, "a_modes", &d).unwrap();
    let control = drs::control::Control { model: s.model.clone(), ..control };
    let (text, _) = drs::results::write(&control, &data, Some(&set), Some(&r), "a_modes trajectory").unwrap();
    let back = drs::results::read(&text).unwrap();
    assert_eq!(back.cycles.len(), 3);
    for (a, b) in back.cycles.iter().zip(&t.cycles) {
        assert!(a.cycle == b.cycle && a.stats.chi2 == b.stats.chi2 && a.stats.rwp == b.stats.rwp && a.marquardt == b.marquardt && a.singular == b.singular && a.chi2_start == b.chi2_start);
        assert!(a.params.len() == b.params.len() && a.params.iter().zip(&b.params).all(|(x, y)| x.name == y.name && x.value == y.value && x.shift == y.shift && x.esd == y.esd && x.computed == y.computed && x.at_limit == y.at_limit));
    }
    let end = drs::control::Control::from_text(&back.control_file).unwrap().model;
    let from_file = drs::lsq::Trajectory::from_end(&end, &set, back.cycles.clone());
    let mut worst_file: f64 = 0.0;
    for (k, mk) in models.iter().enumerate() {
        worst_file = worst_file.max(largest_difference(mk, &from_file.state_at(k)).0);
    }
    println!("the states from the results file within {worst_file:.1e}");
    assert!(worst_file < 1e-9, "{worst_file}");
    // (e) The record as JSON.
    let json = drs::results::trace_json(&t, s.stop.as_ref().unwrap_or(&drs::lsq::Stop::User));
    assert!(json.contains("\"cycles\": [") && json.matches("\"cycle\": ").count() == 3 && json.contains("\"covariance\": [["));
}

#[test]
fn e2_chi2_grid() {
    let (mut m, data, set, _) = case_a();
    drs::lsq::refine(&mut m, &data, &set, drs::lsq::Options { max_cycles: 20, eps: 0.01 }, |_| {});
    let r = drs::lsq::Session::new(m.clone(), data.clone(), set.clone(), drs::lsq::Options { max_cycles: 0, eps: 0.1 }).finish();
    let (a, b) = (set.free.iter().position(|p| p.name == "phase1.A3_GM4+.mode_amplitude").unwrap(), set.free.iter().position(|p| p.name == "phase1.A9_X3+.mode_amplitude").unwrap());
    let (sa, sb) = (r.params[a].esd, r.params[b].esd);
    let n = 5;
    let g = drs::lsq::chi2_grid(&m, &data, &set, (a, b), (sa, sb), n, |_, _| true);
    // Every point once, coarse to fine: the centre, then the corners.
    let mut seen = std::collections::BTreeSet::new();
    assert!(g.order.iter().all(|p| seen.insert(*p)) && seen.len() == n * n && g.computed == n * n);
    assert_eq!(&g.order[..5], &[(2, 2), (0, 0), (4, 0), (0, 4), (4, 4)]);
    // The centre is χ² at the model's values; the converged minimum is the grid's lowest.
    let c = drs::lsq::statistics(&m, &data, &m.calculate(&data.x).ycalc, set.free.len()).chi2;
    assert_eq!(g.chi2[2][2], c);
    assert!(g.chi2.iter().flatten().all(|v| *v >= c - 1e-9));
    // Near the minimum χ² is the quadratic form of the normal matrix: Δχ² = δᵀ A δ /(N − P) over the two parameters
    // (the others fixed), at ± 1 esd within 5 %.
    let (na, _) = drs::lsq::normal_matrix(&m, &data, &set);
    let dof = (r.stats.n - r.stats.p) as f64;
    let mut worst: f64 = 0.0;
    for (i, j) in [(1, 2), (3, 2), (2, 1), (2, 3), (1, 1), (3, 3), (1, 3), (3, 1)] {
        let (da, db) = ((g.xa[i] - g.xa[2]) / set.free[a].fields[0].1, (g.xb[j] - g.xb[2]) / set.free[b].fields[0].1);
        let pred = (da * da * na[a][a] + 2.0 * da * db * na[a][b] + db * db * na[b][b]) / dof;
        let got = g.chi2[i][j] - c;
        worst = worst.max((got / pred - 1.0).abs());
    }
    println!("χ² on the grid as the normal matrix predicts within {:.2} %", 100.0 * worst);
    assert!(worst < 0.05, "{worst}");
    // Stopped after the first batch: the rest NaN.
    let h = drs::lsq::chi2_grid(&m, &data, &set, (a, b), (sa, sb), n, |_, _| false);
    assert!(h.computed < n * n && h.chi2.iter().flatten().filter(|v| v.is_nan()).count() == n * n - h.computed);
    // The coarse-to-fine order on larger grids: the 3 × 3 of stride 4 first on a 9 × 9, every point once on 21 × 21.
    let o = drs::lsq::coarse_to_fine(9);
    assert!(o[..9].iter().all(|(i, j)| i % 4 == 0 && j % 4 == 0));
    let o21 = drs::lsq::coarse_to_fine(21);
    assert_eq!(o21.iter().collect::<std::collections::BTreeSet<_>>().len(), 441);
}

#[test]
fn e3_state() {
    let (m, data, set, _) = case_a();
    let r = drs::lsq::Session::new(m.clone(), data.clone(), set.clone(), drs::lsq::Options { max_cycles: 0, eps: 0.1 }).finish();
    let o = drs::state::Options { widths: Some((20.0, 150.0, 10.0)), dmax: 2.2, dangl: 0.0, reference_cells: vec![Some(drs::cell::Cell::new(5.7, 5.7, 8.0, 90.0, 90.0, 90.0).unwrap())] };
    let st = drs::state::summary(&m, &data, Some(&set), Some(&r), &o);
    let ph = &st.phases[0];
    assert_eq!(ph.name, "Sr2ScSbO6");
    // The widths: H_G² unclamped is the model's, and the strain against the reference.
    for w in &ph.widths {
        assert_eq!(w.h_g2, m.gaussian_h2(&m.phases[0], w.two_theta, 1));
    }
    let s = ph.strain.unwrap();
    assert!((s[2] - (m.phases[0].cell.c / 8.0 - 1.0)).abs() < 1e-15);
    // A negative H_G² is kept negative where the profile clamps it (W = −1).
    let mut neg = m.clone();
    neg.phases[0].w = -1.0;
    let stn = drs::state::summary(&neg, &data, None, None, &o);
    assert!(stn.phases[0].widths.iter().any(|w| w.h_g2 < 0.0));
    // The irreps with their su's, as `drs refine` prints them.
    let ms = m.phases[0].modes.as_ref().unwrap();
    for ((n, a, e), ir) in ph.irreps.iter().zip(ms.irreps()) {
        assert!(n == &ir.name && *a == ir.amplitude && *e >= 0.0);
    }
    assert!(ph.irreps.iter().any(|x| x.2 > 0.0));
    // The contacts: Sb1's six O at 1.98 Å, with su's.
    let sb = ph.labels.iter().position(|l| l == "Sb1").unwrap();
    assert_eq!(ph.contacts[sb].len(), 6);
    assert!(ph.contacts[sb].iter().all(|c| (c.distance - 1.98).abs() < 0.01 && c.sigma > 0.0));
    // The plausible ranges: Sb⁵⁺–O²⁻ and Sc³⁺–O²⁻ in six-fold coordination by bond valence, against case A's.
    let (bv, ionic) = drs::state::plausible(("Sb", 5), ("O", -2), &[4, 6, 8]);
    let d6 = bv[1].1.unwrap();
    let sc = ph.labels.iter().position(|l| l == "Sc1").unwrap();
    let (bv_sc, _) = drs::state::plausible(("Sc", 3), ("O", -2), &[6]);
    let mean = |cs: &[drs::bonds::Contact]| cs.iter().map(|c| c.distance).sum::<f64>() / cs.len() as f64;
    println!("Sb–O: bond valence {d6:.4} (N = 6), ionic radii {ionic:.3}, case A {:.4}; Sc–O: {:.4}, case A {:.4}", mean(&ph.contacts[sb]), bv_sc[0].1.unwrap(), mean(&ph.contacts[sc]));
    assert!((mean(&ph.contacts[sb]) - d6).abs() < 0.05 && (mean(&ph.contacts[sc]) - bv_sc[0].1.unwrap()).abs() < 0.05);
    // A reflection's profile: the profiles of the reflections whose windows hold a point add up to the phase's curve.
    let c = m.calculate(&data.x);
    let pk = c.peaks.iter().filter(|k| k.wavelength == 1).max_by(|x, y| x.intensity.total_cmp(&y.intensity)).unwrap();
    let i0 = data.x.iter().position(|&x| x >= pk.two_theta).unwrap();
    let x0 = data.x[i0];
    let mut sum = 0.0;
    for k in c.peaks.iter().filter(|k| (k.two_theta - x0).abs() <= k.reach) {
        let (x, y, _) = drs::state::reflection_profile(&m, &data, k.phase, k.hkl).unwrap();
        if k.wavelength == 1 {
            if let Some(i) = x.iter().position(|&v| v == x0) {
                sum += y[i];
            }
        }
    }
    let lambda2 = c.peaks.iter().any(|k| k.wavelength == 2);
    println!("at {x0}°: the reflections' profiles {sum:.6}, the phase's curve {:.6} (a second wavelength: {lambda2})", c.per_phase[0][i0]);
    assert!(lambda2 || (sum - c.per_phase[0][i0]).abs() <= 1e-9 * c.per_phase[0][i0]);
}

#[test]
fn e4_scene() {
    let (m, _, _, _) = case_a();
    // The cell's content: 2 Sb + 2 Sc + 4 Sr + 12 O in [0, 1).
    let sc = drs::geometry::scene(&m, 0, [0.0; 3], [1.0; 3], &drs::geometry::SceneRule::Within(2.2)).unwrap();
    let inside = sc.images.iter().filter(|im| im.fract.iter().all(|v| *v < 1.0 - 1e-6)).count();
    assert_eq!(inside, 20);
    // The Cartesian positions are the fractions through the cell's matrix.
    let o = m.phases[0].cell.crystal_to_orthonormal();
    for im in &sc.images {
        for (c, row) in im.cart.iter().zip(&o) {
            assert!((c - (0..3).map(|j| row[j] * im.fract[j]).sum::<f64>()).abs() < 1e-12);
        }
    }
    // Sb1's octahedron: 6 vertices, 8 faces, its volume within 1 % of the regular one of its mean bond.
    let oct = sc.polyhedra.iter().find(|p| sc.images[p.centre].label == "Sb1").unwrap();
    let mean = oct.vertices.iter().map(|&v| sc.bonds.iter().find(|b| (b.0 == oct.centre && b.1 == v) || (b.1 == oct.centre && b.0 == v)).unwrap().2).sum::<f64>() / 6.0;
    println!("Sb1: {} faces, V {:.4} Å³ (regular {:.4}), distortion {:.5}", oct.faces.len(), oct.volume, 4.0 / 3.0 * mean.powi(3), oct.distortion);
    assert!(oct.faces.len() == 8 && (oct.volume / (4.0 / 3.0 * mean.powi(3)) - 1.0).abs() < 0.01);
    // The irreps' displacements add up to each image's distance from its reference.
    for (i, im) in sc.images.iter().enumerate() {
        let total: [f64; 3] = std::array::from_fn(|k| sc.displacements.iter().map(|(_, u)| u[i][k]).sum());
        assert!((0..3).all(|k| (im.cart[k] - sc.reference[i][k] - total[k]).abs() < 1e-12));
    }
    // The hull of a cube: 12 triangles, outward, volume 1.
    let cube: Vec<[f64; 3]> = (0..8).map(|k| [(k & 1) as f64, ((k >> 1) & 1) as f64, ((k >> 2) & 1) as f64]).collect();
    let faces = drs::geometry::hull(&cube);
    let centre = [0.5; 3];
    let vol: f64 = faces.iter().map(|(v, n)| {
        let (a, b, c) = (cube[v[0]], cube[v[1]], cube[v[2]]);
        let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let w = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let cr = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]];
        assert!(cr[0] * n[0] + cr[1] * n[1] + cr[2] * n[2] > 0.0, "a face turned inward");
        let r = [a[0] - centre[0], a[1] - centre[1], a[2] - centre[2]];
        (r[0] * cr[0] + r[1] * cr[1] + r[2] * cr[2]) / 6.0
    }).sum();
    assert!(faces.len() == 12 && (vol - 1.0).abs() < 1e-12, "{} faces, V {vol}", faces.len());
}

/// The path in units of the esd at its end (L139: Explore's traces, its table and the χ² map's span are the core's):
/// on three cycles of case A, every distance is (value − end)/esd_end from the cycles' own numbers, the start's from
/// the start model; the last state is 0; shift/esd is the cycle's; reach is the largest |distance|; pair_span holds
/// the path with 15 % to spare and at least 4σ; a cycle's largest shift is its largest |shift/esd|.
#[test]
fn the_path_in_units_of_the_esd_at_its_end() {
    let (m, data, set, _) = case_a();
    let mut s = drs::lsq::Session::new(m, data, set.clone(), drs::lsq::Options { max_cycles: 3, eps: 1e-6 });
    s.step(3, |_| true);
    let t = s.trajectory();
    let p = t.in_esd(None);
    assert_eq!(p.names.len(), set.free.len());
    assert_eq!(p.states, vec![0, 1, 2, 3]);
    let last = t.cycles.last().unwrap();
    for (j, free) in set.free.iter().enumerate() {
        let q = last.params.iter().find(|q| q.name == free.name).unwrap();
        assert_eq!((p.end_value[j], p.end_esd[j]), (q.value, q.esd));
        if q.esd > 0.0 {
            assert_eq!(p.from_end[3][j], Some(0.0), "{} ends at 0", free.name);
            let start = (free.value_in(&t.start) - q.value) / q.esd;
            assert!((p.from_end[0][j].unwrap() - start).abs() <= 1e-12 * start.abs().max(1.0));
            let c1 = t.cycles[0].params.iter().find(|x| x.name == free.name).unwrap();
            assert_eq!(p.from_end[1][j], Some((c1.value - q.value) / q.esd));
            assert_eq!(p.shift_esd[1][j], Some(c1.shift_esd));
            let reach = p.from_end.iter().filter_map(|r| r[j]).fold(0.0_f64, |a, d| a.max(d.abs()));
            assert_eq!(p.reach[j], reach);
        }
        assert_eq!(p.shift_esd[0][j], None, "no shift at the start");
    }
    // the esd's at the end may be given (a finished run's)
    let given: Vec<f64> = p.end_esd.iter().map(|e| 2.0 * e).collect();
    let p2 = t.in_esd(Some(&given));
    let j = (0..p.names.len()).find(|&j| p.end_esd[j] > 0.0).unwrap();
    assert!((p2.from_end[0][j].unwrap() - p.from_end[0][j].unwrap() / 2.0).abs() < 1e-12 * p.from_end[0][j].unwrap().abs().max(1.0));
    assert_eq!(drs::lsq::pair_span(1.0, 2.0), 4.0 * 1.15);
    assert_eq!(drs::lsq::pair_span(10.0, 2.0), 10.0 * 1.15);
    let (w, name) = last.largest_shift().unwrap();
    assert!(last.params.iter().filter(|q| q.esd > 0.0).all(|q| q.shift_esd.abs() <= w));
    assert!(last.params.iter().any(|q| q.name == name && q.shift_esd.abs() == w));
}

#[test]
fn e14a_esds_at_every_cycle() {
    let (m, data, set, _) = case_a();
    let mut s = drs::lsq::Session::new(m, data.clone(), set.clone(), drs::lsq::Options { max_cycles: 3, eps: 1e-6 });
    s.step(3, |_| true);
    let t = s.trajectory();
    let r = s.finish();
    let opts = drs::state::Options::default();
    let irreps = |st: &drs::state::State| st.phases[0].irreps.clone();
    // At the end, from the refinement and from its covariance given directly: the same, to the bit.
    let end = drs::state::summary(&s.model, &data, Some(&set), Some(&r), &opts);
    let direct = drs::state::summary_at(&s.model, &data, Some(&set), Some(&r.covariance), &opts);
    assert_eq!(irreps(&end), irreps(&direct));
    assert_eq!(end.phases[0].modes, direct.phases[0].modes);
    // At every cycle, from that cycle's covariance (the normal matrix at the cycle's start, times its χ²).
    for k in 1..=3 {
        let st = drs::state::summary_at(&t.state_at(k), &data, Some(&set), Some(&t.cycles[k - 1].covariance), &opts);
        assert!(st.phases[0].irreps.iter().filter(|i| i.1 != 0.0).all(|i| i.2 > 0.0), "cycle {k}: every irrep refined has an esd");
        if k == 3 {
            // Against the end's: the cycle's normal matrix is the one at its start values, the end's at its last; the
            // difference is that of the matrix between the two (measured: below 10⁻³ of each esd here).
            for (a, b) in irreps(&st).iter().zip(irreps(&end)) {
                let rel = if b.2 > 0.0 { (a.2 - b.2).abs() / b.2 } else { a.2 };
                println!("{:6} esd at cycle 3 {:.6e}, at the end {:.6e} ({rel:.1e})", a.0, a.2, b.2);
                assert!(rel < 1e-3, "{}: {rel}", a.0);
            }
        }
    }
    // A cycle read back from a results file keeps no covariance: no esd's.
    let none = drs::state::summary_at(&s.model, &data, Some(&set), Some(&[]), &opts);
    assert!(none.phases[0].irreps.iter().all(|i| i.2 == 0.0) && none.phases[0].modes.iter().all(|m| m.2 == 0.0));
}

#[test]
fn e14b_the_sign_reference() {
    let (m, _, _, _) = case_a();
    let d = root("../manuals/fullprof/runs/modes/A_modes");
    let ph = &m.phases[0];
    let drs::powder::Intensities::Structure(st, _) = &ph.intensities else { panic!() };
    let species: Vec<String> = st.atoms.iter().map(|x| x.species.clone()).collect();
    // The control file's key, its file FullProf's .sum beside the run: the reference is FullProf's amplitudes.
    let pcr = drs::pcr::read(&latin1(&d.join("a_modes.pcr"))).unwrap();
    let c = drs::control::Control::from_pcr(&pcr, &d, "a_modes", &d).unwrap();
    let with = |key: &str| c.to_text().replacen("[[phase]]\n", &format!("[[phase]]\n{key}\n"), 1);
    let c2 = drs::control::Control::from_text_in(&with("mode_reference = { file = \"a_modes.sum\" }"), &d).unwrap();
    let ms2 = c2.model.phases[0].modes.as_ref().unwrap();
    let sr = ms2.sign_reference.as_ref().unwrap();
    assert_eq!(sr.amplitudes.len(), 12);
    assert_eq!(sr.amplitudes[2], ("A3_GM4+".to_string(), -0.828632));
    // Written back and read again, the same.
    let c3 = drs::control::Control::from_text_in(&c2.to_text(), &d).unwrap();
    assert_eq!(c3.model.phases[0].modes.as_ref().unwrap().sign_reference, ms2.sign_reference);
    // An AMPLIMODES analysis as the reference: its amplitudes by mode name and the parent it states (a .pcr states
    // none). Case A's analysis 3 is in its P2₁/c setting, analysis 1 (the phase's P2₁/n) has no block for FullProf.
    let ampli = "../../../../bcs/mode-crystallography-manual/course-um6p-2024/04_tools_materials/0_1_P21nToCubicTransData.txt";
    let c4 = drs::control::Control::from_text_in(&with(&format!("mode_reference = {{ file = \"{ampli}\", analysis = 3 }}")), &d).unwrap();
    let ms4 = c4.model.phases[0].modes.as_ref().unwrap();
    assert_eq!(ms4.parent.as_ref().map(|p| p.group.as_str()), Some("225"), "the parent from the analysis");
    assert_eq!(ms4.sign_reference.as_ref().unwrap().amplitudes[0], ("A1_GM1+".to_string(), -0.242338));
    let e = drs::control::Control::from_text_in(&with(&format!("mode_reference = {{ file = \"{ampli}\", analysis = 1 }}")), &d).unwrap_err();
    assert!(e.iter().any(|x| x.message.contains("analysis 1 has no block for FullProf")), "{e:?}");
    // Per state: the set nearest FullProf's amplitudes and the operation. FullProf's own: its setting. The mirror
    // solution of Explore's +0.1 Å starts (e6_e7_modes): brought back by the parent's 2-fold about c.
    let mut ms = ph.modes.clone().unwrap();
    ms.sign_reference = ms2.sign_reference.clone();
    let metric = drs::modes::read_amplimodes(&latin1(&root("../manuals/bcs/mode-crystallography-manual/course-um6p-2024/04_tools_materials/0_1_P21nToCubicTransData.txt")))[0].reference_cell().unwrap();
    let q: Vec<f64> = ms.modes.iter().map(|x| x.amplitude).collect();
    let own = ms.reference_signs(&species, &ph.group, &metric, &q).unwrap().unwrap();
    assert!(own.is_identity() && own.words() == "the reference's setting", "{own:?}");
    let name = |n: &str| ms.modes.iter().position(|x| x.name == n).unwrap();
    let mut mirror = q.clone();
    for n in ["A3_GM4+", "A10_X5+", "A11_X5+"] {
        mirror[name(n)] = -mirror[name(n)];
    }
    let back = ms.reference_signs(&species, &ph.group, &metric, &mirror).unwrap().unwrap();
    println!("the mirror: {}", back.words());
    assert!(back.words().contains("2 about [0 0 1]"), "{}", back.words());
    assert!(back.amplitudes.iter().zip(&q).all(|(a, b)| (a - b).abs() < 1e-6));
    // The command line: drs modes on the control file names the state against its reference.
    let f = root("target/test-e14b.drc");
    let text = with(&format!("mode_reference = {{ file = \"{}\" }}", d.join("a_modes.sum").display())).replace("\"SrScLRRT.dat\"", &format!("\"{}\"", d.join("SrScLRRT.dat").display()));
    std::fs::write(&f, text).unwrap();
    let o = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(["modes", f.to_str().unwrap()]).output().unwrap();
    assert!(String::from_utf8_lossy(&o.stdout).contains("sign reference"), "{}{}", String::from_utf8_lossy(&o.stdout), String::from_utf8_lossy(&o.stderr));
}

#[test]
fn e14c_fullprof_reported_amplitudes() {
    let d = root("../manuals/fullprof/runs/modes/A_modes");
    let b = drs::modes::read_fullprof_amplitudes(&latin1(&d.join("a_modes.sum")));
    assert_eq!(b.len(), 1);
    let a = &b[0];
    // Twelve modes (A1–A12) and seven irreps, as a_modes.sum prints them (lines 58–71, 147–153).
    assert_eq!((a.modes.len(), a.irreps.len()), (12, 7));
    assert_eq!(a.modes[0], ("A1_GM1+".to_string(), -0.283144, 0.005271));
    assert_eq!(a.modes[2], ("A3_GM4+".to_string(), -0.828632, 0.003009));
    assert_eq!(a.modes[10], ("A11_X5+".to_string(), -0.051054, 0.008631));
    let ir = |n: &str| a.irreps.iter().find(|x| x.0 == n).unwrap().clone();
    assert_eq!(ir("GM4+").1, 0.829);
    assert!((ir("GM4+").2.unwrap() - 0.003).abs() < 1e-15);
    assert_eq!(ir("X3+").1, 0.4);
    assert!((ir("X5+").2.unwrap() - 0.006).abs() < 1e-15);
    assert_eq!(ir("GM3+"), ("GM3+".to_string(), 0.0, None));
    // The calc pass's .sum prints the phase's results twice: one entry, the same values.
    let c = drs::modes::read_fullprof_amplitudes(&latin1(&d.join("a_modes_calc.sum")));
    assert_eq!(c.len(), 1);
    assert_eq!((c[0].modes.clone(), c[0].irreps.clone()), (a.modes.clone(), a.irreps.clone()));
    // Beside the core's: the .pcr's modes by name. The .pcr holds the values FullProf read (it is a_modes.pcr.in),
    // the .sum those after its last cycle, whose shifts FullProf stops below Eps·σ (Eps = 0.2 in the .pcr; measured:
    // A10_X5+ 0.110418 read, 0.110325 after, 0.019σ, the largest); the core's irreps from the .pcr's amplitudes
    // within that and the .sum's print.
    let (m, _, _, _) = case_a();
    let ms = m.phases[0].modes.as_ref().unwrap();
    assert_eq!(ms.modes.iter().map(|x| x.name.clone()).collect::<Vec<_>>(), a.modes.iter().map(|x| x.0.clone()).collect::<Vec<_>>());
    for (x, y) in ms.modes.iter().zip(&a.modes) {
        assert!((x.amplitude - y.1).abs() <= 0.2 * y.2 + 5e-7, "{}: {} against {} ({})", x.name, x.amplitude, y.1, y.2);
    }
    let worst = ms.modes.iter().zip(&a.modes).filter(|(_, y)| y.2 > 0.0).map(|(x, y)| (x.amplitude - y.1).abs() / y.2).fold(0.0, f64::max);
    println!("the .pcr's amplitudes against the .sum's: at most {worst:.3}σ");
    for i in ms.irreps() {
        let (_, v, s) = ir(&i.name);
        let half = s.map_or(5e-6, |e| 5e-4 + 0.2 * e);
        assert!((i.amplitude - v).abs() <= half, "{}: {} against {v}", i.name, i.amplitude);
    }
    // The command line: drs modes on the .pcr prints them.
    let o = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(["modes", d.join("a_modes.pcr").to_str().unwrap()]).output().unwrap();
    let out = String::from_utf8_lossy(&o.stdout);
    assert!(out.contains("FullProf's reported amplitudes") && out.contains("A3_GM4+       -0.828632 (0.003009)"), "{out}");
}

/// Publish's P21 (L218): every set `equivalent_amplitudes` offers keeps each irrep's amplitude (a parent operation maps
/// each irrep's subspace onto itself), and case A has 4. Before the whole cell was checked it offered 24, 20 of them
/// moving amplitude between irreps (-y,x,z: GM4+ 0.8288 Å to 0.0363): the asymmetric unit of another domain,
/// symmetrised again by the phase's operations, which is another structure.
#[test]
fn p21_equivalent_sets_keep_each_irreps_amplitude() {
    let (m, _, _, analyses) = case_a();
    let ph = &m.phases[0];
    let ms = ph.modes.as_ref().unwrap();
    let drs::powder::Intensities::Structure(st, _) = &ph.intensities else { panic!() };
    let species: Vec<String> = st.atoms.iter().map(|x| x.species.clone()).collect();
    let metric = analyses[0].reference_cell().unwrap();
    let q: Vec<f64> = ms.modes.iter().map(|x| x.amplitude).collect();
    let eq = ms.equivalent_amplitudes(&species, &ph.group, &metric, &q).unwrap();
    let ops: Vec<&str> = eq.iter().map(|e| e.operation.as_str()).collect();
    assert_eq!(ops, ["x,y,z", "x,y+1/2,z+1/2", "-x,-y,z", "-y,-x,-z"]);
    let own: Vec<(String, f64)> = ms.irreps().iter().map(|i| (i.name.clone(), i.amplitude)).collect();
    for e in &eq {
        let mut other = ms.clone();
        for (mode, a) in other.modes.iter_mut().zip(&e.amplitudes) {
            mode.amplitude = *a;
        }
        for (i, (name, a)) in other.irreps().iter().zip(&own) {
            assert_eq!(&i.name, name);
            assert!((i.amplitude - a).abs() < 1e-6, "{}: {} keeps {name}'s amplitude {a}, not {}", e.operation, name, i.amplitude);
        }
        // and each is the same structure: the moved cell, rebuilt from its amplitudes, is case A's own cell moved.
        assert!(ms.same_structure(&species, &ph.group, &metric, &q, &e.amplitudes, 1e-6).unwrap().is_some());
    }
}
