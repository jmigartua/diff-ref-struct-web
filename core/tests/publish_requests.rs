//! The Publish stream's requests to the core (publish/REQUESTS.md P2–P11; L114), each against FullProf where FullProf
//! prints the number:
//! - P2: the background-free R factors (FullProf's "conventional" Rp, Rwp: yᵢ → yᵢ − bᵢ in the denominators, manual
//!   p. 44, the net counts summed with their signs as FullProf does) at FullProf's converged values of Y2O3, LaB6 and
//!   rutile–anatase, over all non-excluded points, to FullProf's print and the profile's known difference; the results file of case
//!   A holds R_F, the background-free R's, the mode amplitudes and irreps with su's, the reference structure, the parent
//!   cell, shift/su_max and how the su's were obtained;
//! - P3: case A's positions reduced to the cell with su's, each within FullProf's CFL print (its esd's those of the
//!   amplitudes taken as independent), and its Sb1–O11 distance against FullProf's Bond_Str (1.9778(17),
//!   `a_modes_1.dis`): Bond_Str's number is the CFL's, positions and esd's rounded, which the core reproduces to the
//!   digit; the core's own, 1.9777(16) at unrounded positions with the covariance, differs from it by less than the
//!   CFL's rounding moves the distance;
//! - P5: charges from the control file or the caller give valence sums; the species keeps its scattering;
//! - P7: Y2O3's O–Y1–O angle of 180°, fixed by the inversion centre, without su, in the CIF and (LaMnO3's) from a
//!   CFL; each bond once;
//! - P8: the CIF's items (profile, background and preferred-orientation functions, constraints, R_F², weighting,
//!   shift/su_max, block ids, computing) and unique labels (Y1, Y2);
//! - P9: AMPLIMODES's summary per irrep (case A's analysis: k, irrep, direction, isotropy subgroup, dimension);
//! - P10: the phases' proposed names (Y2O3, Sr2ScSbO6), with no ':' or FullProf's keywords;
//! - P11: case A's cell reduced to the parent's.

use drs::control::Control;
use drs::powder::Model;
use drs::refinable::ParameterSet;
use std::path::{Path, PathBuf};

fn runs() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs")
}

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display())).iter().map(|&b| b as char).collect()
}

fn load(dir: &Path, file: &str, code: &str) -> (drs::pcr::Pcr, drs::powder::Loaded) {
    let mut p = drs::pcr::read(&latin1(&dir.join(file))).unwrap();
    p.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
    let l = Model::from_pcr(&p, dir, code).unwrap_or_else(|e| panic!("{file}: {e}"));
    (p, l)
}

/// A state with its refinement (n cycles from the file's values) and propagator.
fn refined(dir: &Path, file: &str, code: &str, cycles: usize) -> (Model, drs::data::Pattern, ParameterSet, drs::lsq::Refinement) {
    let (p, l) = load(dir, file, code);
    let data = l.data.clone().unwrap();
    let mut m = l.model.clone();
    let set = ParameterSet::from_pcr(&p, &m).unwrap();
    let r = drs::lsq::refine(&mut m, &data, &set, drs::lsq::Options { max_cycles: cycles, eps: 0.1 }, |_| {});
    (m, data, set, r)
}

/// Case A's parent (Fm-3m, P from AMPLIMODES's analysis 1, as the control file's mode_parent_* give it).
fn case_a_parent() -> drs::modes::ModeParent {
    let f = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/bcs/mode-crystallography-manual/course-um6p-2024/04_tools_materials/0_1_P21nToCubicTransData.txt");
    drs::modes::read_amplimodes(&latin1(&f))[0].mode_parent().unwrap()
}

#[test]
fn p2_background_free_r_factors() {
    // FullProf's "Conventional Rietveld Rp,Rwp,Re and Chi2" of the last cycle of each converged run, at the values it
    // wrote (the core at the same values, 0 cycles).
    let mut worst: f64 = 0.0;
    for (run, code) in [("Y2O3", "y2o3_structure"), ("LaB6", "LaB6rx"), ("rutile-anatase", "rutana_tch")] {
        let d = runs().join(run);
        let out = latin1(&d.join(format!("{code}_conv.out")));
        let all = out.rfind("ALL NON-EXCLUDED POINTS").unwrap();
        let line = out[all..].lines().find(|l| l.contains("Conventional Rietveld Rp,Rwp,Re and Chi2:")).unwrap();
        let v: Vec<f64> = line.split(':').nth(1).unwrap().split_whitespace().map(|x| x.parse().unwrap()).collect();
        let (_, l) = load(&d, &format!("{code}_conv.pcr"), code);
        let data = l.data.unwrap();
        let s = drs::lsq::statistics(&l.model, &data, &l.model.calculate(&data.x).ycalc, 0);
        let (rp, rwp) = (s.rp_background_free, s.rwp_background_free);
        // The core's Rp and Rwp themselves against FullProf's at these values (the profile's known difference).
        let plain = out.lines().rev().find(|l| l.contains("=> R-Factors:") && l.contains("Patt#:  1")).unwrap();
        let w: Vec<&str> = plain.split_whitespace().collect();
        let (fp_rp, fp_rwp): (f64, f64) = (w[2].parse().unwrap(), w[3].parse().unwrap());
        let d_plain = (s.rp - fp_rp).abs().max((s.rwp - fp_rwp).abs());
        let d_conv = (rp - v[0]).abs().max((rwp - v[1]).abs());
        println!("{run}: background-free Rp {rp:.2} Rwp {rwp:.2}, FullProf's conventional {} {}; Rp {:.2} Rwp {:.2}, FullProf's {fp_rp} {fp_rwp}", v[0], v[1], s.rp, s.rwp);
        // Within FullProf's print (0.05) and as far as the plain R's differ, scaled by the denominators' ratio.
        worst = worst.max(d_conv - 0.05 - d_plain * (rp / s.rp.max(1e-9)));
    }
    assert!(worst <= 0.0, "background-free R's beyond FullProf's by {worst}");
    // The results file of case A (one cycle).
    let d = runs().join("modes/A_modes");
    let (p, _) = load(&d, "a_modes.pcr", "a_modes");
    let control = Control::from_pcr(&p, &d, "a_modes", &d).unwrap();
    let (mut m, data, set, r) = refined(&d, "a_modes.pcr", "a_modes", 1);
    m.phases[0].modes.as_mut().unwrap().parent = Some(case_a_parent());
    let control = Control { model: m, ..control };
    let (text, structures) = drs::results::write(&control, &data, Some(&set), Some(&r), "a_modes 2026-10-11 12.00.00").unwrap();
    for item in [
        "_drs_phase_R_F",
        "_drs_prof_R_factor_background_free",
        "_drs_prof_wR_factor_background_free",
        "_refine_ls_shift/su_max",
        "_refine_special_details",
        "_drs_mode_amplitude_su",
        "_drs_irrep_amplitude_su",
        "_drs_mode_reference_fract_x",
        "_drs_parent_cell_a",
        "_refine_ls_R_factor_all",
    ] {
        assert!(text.contains(item), "the results file lacks {item}");
    }
    let back = drs::results::read(&text).unwrap();
    assert_eq!(back.refined.len(), r.fields.len());
    assert_eq!(structures[0].0, "Sr2ScSbO6");
}

#[test]
fn p3_mode_positions_and_distances() {
    let d = runs().join("modes/A_modes");
    let (m, _, set, r) = refined(&d, "a_modes.pcr", "a_modes", 0);
    let prop = drs::geometry::Propagator::new(&m, &set, &r);
    let g = drs::geometry::geometry(&m, 0, Some(&prop), drs::geometry::BondRule::Within(2.1)).unwrap();
    for (l, x) in g.labels.iter().zip(&g.positions) {
        assert!(x.iter().all(|v| (0.0..1.0).contains(&v.0)), "{l}: {x:?}");
    }
    // FullProf's CFL (`a_modes_1.cfl`), which Bond_Str reads: positions rounded to their esd's last digit, the esd's
    // from the amplitudes taken as independent (`symmetry_modes`), rounded.
    let cfl = drs::cfl::read(&latin1(&d.join("a_modes_1.cfl"))).unwrap();
    let at = |l: &str| cfl.structure.atoms.iter().position(|a| a.label == l).unwrap();
    let o11 = g.labels.iter().position(|l| l == "O11").unwrap();
    let gi = drs::geometry::geometry(&m, 0, Some(&prop.independent()), drs::geometry::BondRule::Within(2.1)).unwrap();
    for k in 0..3 {
        let (x, e) = (cfl.structure.atoms[at("O11")].x[k], cfl.uncertainties.x[at("O11")][k]);
        let places = (-e.log10()).round() as i32;
        let h = 0.5 * 10f64.powi(-places);
        assert!(((g.positions[o11][k].0 - x + 0.5).rem_euclid(1.0) - 0.5).abs() <= h + 1e-9, "O11 {k}: {:?} against {x}", g.positions[o11][k]);
        assert!(g.positions[o11][k].1 > 0.0 && (gi.positions[o11][k].1 - e).abs() <= h + 1e-9, "O11 {k}: su {:?} against {e}", gi.positions[o11][k]);
    }
    // Sb1–O11: FullProf's Bond_Str 1.9778(17) is the CFL's (the core's CFL path gives it to the digit, `bond_str`);
    // the core at its unrounded positions differs from it by what the CFL's rounding does to the distance, and its su,
    // through the covariance, from the CFL's by the rounding of O11's esd's (0.00026 printed 0.0003).
    let (sb, o) = (at("Sb1"), at("O11"));
    let c = drs::bonds::contacts(&cfl.structure, Some(&cfl.uncertainties), 2.1)[sb].iter().find(|c| c.atom == o).unwrap().clone();
    let from_cfl = (c.distance, drs::bonds::bond_str_distance_sigma(&cfl.structure, &cfl.uncertainties, sb, &c));
    let sbg = g.labels.iter().position(|l| l == "Sb1").unwrap();
    let b = g.bonds.iter().find(|b| (b.i == sbg && b.j == o11) || (b.i == o11 && b.j == sbg)).unwrap();
    println!("case A: O11 at {:?}; Sb1–O11 {:.6}({:.6}) = {}; from FullProf's CFL {:.6}({:.6}); FullProf's Bond_Str 1.9778(17)", g.positions[o11], b.distance.0, b.distance.1, drs::cif_write::with_su(b.distance.0, b.distance.1, 4), from_cfl.0, from_cfl.1);
    assert!((from_cfl.0 - 1.9778).abs() <= 0.5e-4 && (from_cfl.1 - 0.0017).abs() <= 0.5e-4, "{from_cfl:?}");
    // The CFL's rounding: O11 within half a unit of its printed digits moves the distance by at most Σ|∂d/∂x|·h.
    let grad: f64 = (0..3).map(|k| {
        let mut mv = m.clone();
        let h = 1e-6;
        let mut st2 = drs::geometry::structure(&m, 0).unwrap();
        st2.atoms[o11].x[k] += h;
        if let drs::powder::Intensities::Structure(s, _) = &mut mv.phases[0].intensities {
            **s = st2;
        }
        let g2 = drs::geometry::geometry(&mv, 0, None, drs::geometry::BondRule::Within(2.1)).unwrap();
        let b2 = g2.bonds.iter().find(|x| x.i == b.i && x.j == b.j && x.op == b.op && x.translation == b.translation).unwrap();
        let places = (-cfl.uncertainties.x[at("O11")][k].log10()).round() as i32;
        ((b2.distance.0 - b.distance.0) / h).abs() * 0.5 * 10f64.powi(-places)
    }).sum();
    println!("the CFL's rounding moves Sb1–O11 by up to {grad:.6} Å; here {:.6}", (b.distance.0 - from_cfl.0).abs());
    assert!((b.distance.0 - from_cfl.0).abs() <= grad, "{} against the CFL's {}, rounding bound {grad}", b.distance.0, from_cfl.0);
    assert!(b.distance.1 < from_cfl.1 && b.distance.1 > 0.0);
}

#[test]
fn p5_charges() {
    let d = runs().join("Y2O3");
    let (p, l) = load(&d, "y2o3_structure_conv.pcr", "y2o3_structure");
    let st = drs::geometry::structure(&l.model, 0).unwrap();
    assert!(st.atoms.iter().all(|a| drs::bonds::atom_ion(a).is_none()));
    let charged = drs::bonds::with_charges(&st, &[("Y".into(), 3), ("O".into(), -2)]);
    let contacts = drs::bonds::contacts(&charged, None, 3.0);
    let occ = drs::bonds::bond_str_occupancies(&charged).unwrap();
    let (sums, _) = drs::bonds::valence_sums(&charged, &contacts, &occ, 0.2).unwrap();
    println!("Y2O3 valence sums: {:?}", sums.iter().map(|v| (charged.atoms[v.atom].label.clone(), v.sum.0)).collect::<Vec<_>>());
    assert!(sums.iter().all(|v| (v.sum.0.abs() - v.valence.abs()).abs() < 0.3));
    // In the control file: the charge is kept and read back, the species unchanged (it scatters as Y).
    let mut c = Control::from_pcr(&p, &d, "y2o3_structure", &d).unwrap();
    if let drs::powder::Intensities::Structure(st, _) = &mut c.model.phases[0].intensities {
        for (a, q) in st.atoms.iter_mut().zip([3, 3, -2]) {
            a.charge = Some(q);
        }
    }
    let text = c.to_text();
    assert!(text.contains("charge = 3") && text.contains("charge = -2"));
    let back = Control::from_text(&text).unwrap();
    let st2 = drs::geometry::structure(&back.model, 0).unwrap();
    assert_eq!(st2.atoms.iter().map(|a| (a.species.clone(), a.charge)).collect::<Vec<_>>(), vec![("Y".to_string(), Some(3)), ("Y".to_string(), Some(3)), ("O".to_string(), Some(-2))]);
    assert_eq!(c.model.calculate(&[30.0, 31.0]).ycalc, back.model.calculate(&[30.0, 31.0]).ycalc);
}

#[test]
fn p7_p8_p10_the_cif() {
    // P7 on a CFL (`drs bonds`): LaMnO3's O1–Mn1–O1, fixed at 180° by the inversion at Mn1, without su (it was
    // 180.0000(9) from the arc cosine's rounding).
    let cfl = drs::cfl::read(&latin1(&runs().join("../7-FullProf-Tutorial-SymmetryModes/LaMnO3/LaMnO3_cub_1.cfl"))).unwrap();
    let contacts = drs::bonds::contacts(&cfl.structure, Some(&cfl.uncertainties), cfl.dmax.0);
    let straight: Vec<_> = drs::bonds::angles(&cfl.structure, Some(&cfl.uncertainties), &contacts, 2.5).into_iter().filter(|a| (a.angle.0 - 180.0).abs() < 1e-6).collect();
    assert!(!straight.is_empty() && straight.iter().all(|a| a.angle.0 == 180.0 && a.angle.1 == 0.0), "{:?}", straight.iter().map(|a| a.angle).collect::<Vec<_>>());
    let d = runs().join("Y2O3");
    let (mut m, data, set, r) = refined(&d, "y2o3_structure_conv.pcr", "y2o3_structure", 0);
    // Charges, so that the bonds are Y–O.
    if let drs::powder::Intensities::Structure(st, _) = &mut m.phases[0].intensities {
        for (a, q) in st.atoms.iter_mut().zip([3, 3, -2]) {
            a.charge = Some(q);
        }
    }
    let text = drs::cif_write::deposition(&m, &data, Some(&set), Some(&r), "Y2O3").unwrap();
    let b = drs::cif::parse(&text);
    let blk = b.iter().find(|x| x.column("_atom_site_label").is_some()).unwrap();
    // P10: the data block's name; P8: unique labels.
    assert!(text.contains("\ndata_Y2O3\n"), "{}", &text[..300]);
    assert_eq!(blk.column("_atom_site_label").unwrap(), vec!["Y1", "Y2", "O"]);
    // P7: the angle O–Y1–O of 180° without su, every other with its own (Y2O3's positions are not refined: the
    // angles carry none, the distances the cell's); each pair once, from the cation.
    let angles = blk.column("_geom_angle").unwrap();
    assert!(angles.contains(&"180"), "{angles:?}");
    assert!(angles.iter().all(|a| !a.starts_with("180.") ), "{angles:?}");
    let l1 = blk.column("_geom_bond_atom_site_label_1").unwrap();
    assert!(l1.iter().all(|l| l.starts_with('Y')), "{l1:?}");
    assert!(blk.column("_geom_bond_distance").unwrap().iter().all(|d| d.contains('(')));
    // P8: the items.
    for item in ["_pd_proc_ls_profile_function", "_pd_proc_ls_background_function", "_pd_proc_ls_pref_orient_corr"] {
        let v = blk.get(item).unwrap();
        assert!(v != "?" && !v.is_empty(), "{item}: {v}");
    }
    for item in ["_refine_ls_number_constraints", "_refine_ls_r_fsqd_factor", "_refine_ls_weighting_scheme", "_refine_ls_shift/su_max", "_pd_block_id", "_pd_block_diffractogram_id", "_computing_publication_material", "_refine_ls_r_factor_all", "_refine_special_details"] {
        assert!(blk.get(item).is_some(), "{item} missing");
    }
    println!("Y2O3: profile '{}'; background '{}'; preferred orientation '{}'; constraints {}; R_F² {}", blk.get("_pd_proc_ls_profile_function").unwrap(), blk.get("_pd_proc_ls_background_function").unwrap(), blk.get("_pd_proc_ls_pref_orient_corr").unwrap(), blk.get("_refine_ls_number_constraints").unwrap(), blk.get("_refine_ls_r_fsqd_factor").unwrap());
    // P10: case A's name, without FullProf's keywords.
    let (m, ..) = refined(&runs().join("modes/A_modes"), "a_modes.pcr", "a_modes", 0);
    assert_eq!(drs::cif_write::phase_name(&m, 0), "Sr2ScSbO6");
}

#[test]
fn p9_amplimodes_summary() {
    let f = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/bcs/mode-crystallography-manual/course-um6p-2024/04_tools_materials/0_1_P21nToCubicTransData.txt");
    let a = drs::modes::read_amplimodes(&latin1(&f));
    let s = &a[0].summary;
    let row = |irrep: &str| s.iter().find(|x| x.irrep == irrep).unwrap();
    assert_eq!(s.len(), 7);
    let gm4 = row("GM4+");
    assert_eq!((gm4.k.as_str(), gm4.direction.as_str(), gm4.isotropy.as_str(), gm4.dimension), ("(0,0,0)", "(a,a,0)", "C2/m (12)", 1));
    let x5 = row("X5+");
    assert_eq!((x5.k.as_str(), x5.direction.as_str(), x5.isotropy.as_str(), x5.dimension), ("(0,1,0)", "(a,a,0,0,a,-a)", "Pnnm (58)", 3));
}

#[test]
fn p11_parent_cell() {
    let d = runs().join("modes/A_modes");
    let (mut m, _, set, r) = refined(&d, "a_modes.pcr", "a_modes", 0);
    m.phases[0].modes.as_mut().unwrap().parent = Some(case_a_parent());
    let prop = drs::geometry::Propagator::new(&m, &set, &r);
    let c = drs::geometry::parent_cell(&m, 0, Some(&prop)).unwrap();
    let cell = m.phases[0].cell;
    println!("case A reduced to the parent: {:?} (the phase {:.5} {:.5} {:.5} {:.3})", c, cell.a, cell.b, cell.c, cell.beta);
    // The parent's F cell holds 2 × 2 × 2 primitive perovskites: 2a_p = √2·a_m ≈ √2·b_m ≈ c_m.
    let ap = [cell.a * 2f64.sqrt(), cell.b * 2f64.sqrt(), cell.c];
    assert!(c[..3].iter().all(|x| ap.iter().any(|y| (x.0 - y).abs() < 0.05)) && c.iter().all(|x| x.1 > 0.0 || x.0 == 90.0));
}
