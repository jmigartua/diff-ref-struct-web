//! The open requests to the core of L153's step 2 (L168):
//! - Publish P15: U_eq of an anisotropic atom, a third of the trace of U in Cartesian axes, against that trace
//!   computed directly on a triclinic cell; PbSO₄'s results file writes it, with its su, where it wrote 0;
//! - Publish P17: the cell reduced to the parent's for a phase refined by its coordinates, its parent given in the
//!   control file (case A's structure as coordinates): the same cell as the mode phase's, kept by the control file,
//!   written with su's in the results file;
//! - Publish P18: the coordination, mean distance and bond-valence sum of each atom with a charge, with su's through
//!   the covariance, in the results file; the values those of `bonds::valence_sums`;
//! - app 31: a file named by its place in the project before its content (`describe_in`): a results file is a result;
//! - app 32: FullProf's exports in `Exports/FullProf/`, FullProf/ left to the original files;
//! - Explore E7 on the command line: `drs modes --equivalent --reference`.

use drs::control::Control;
use drs::powder::{Intensities, Model};
use drs::refinable::ParameterSet;
use std::path::{Path, PathBuf};

mod common;
use common::latin1;

fn root(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)
}

const AMPLIMODES_A: &str = "../manuals/bcs/mode-crystallography-manual/course-um6p-2024/04_tools_materials/0_1_P21nToCubicTransData.txt";

/// Case A (Sr2ScSbO6 by modes) with its parent from AMPLIMODES's analysis 1, its control file and data.
fn case_a() -> (Control, drs::data::Pattern) {
    let d = root("../manuals/fullprof/runs/modes/A_modes");
    let p = drs::pcr::read(&latin1(&d.join("a_modes.pcr"))).unwrap();
    let mut c = Control::from_pcr(&p, &d, "a_modes", &d).unwrap();
    let a = drs::modes::read_amplimodes(&latin1(&root(AMPLIMODES_A)));
    c.model.phases[0].modes.as_mut().unwrap().parent = a[0].mode_parent();
    let data = Model::from_pcr(&p, &d, "a_modes").unwrap().data.unwrap();
    (c, data)
}

#[test]
fn p15_u_equivalent() {
    // A triclinic cell and a β: U_eq against ⅓ tr(A N U N Aᵀ), A the orthonormalisation, N = diag(a*), U^ij from β.
    let cell = drs::cell::Cell::new(5.1, 6.3, 7.7, 81.0, 97.5, 104.0).unwrap();
    let beta = [0.0041, 0.0029, 0.0018, 0.0006, -0.0004, 0.0003];
    let r = cell.reciprocal();
    let n = [r.a, r.b, r.c];
    let pi2 = std::f64::consts::PI * std::f64::consts::PI;
    let idx = [[0, 3, 4], [3, 1, 5], [4, 5, 2]];
    let u = |i: usize, j: usize| beta[idx[i][j]] / (2.0 * pi2 * n[i] * n[j]);
    let a = cell.crystal_to_orthonormal();
    let mut trace = 0.0;
    for row in &a {
        for i in 0..3 {
            for j in 0..3 {
                trace += row[i] * n[i] * u(i, j) * n[j] * row[j];
            }
        }
    }
    let ueq = drs::geometry::u_equivalent(&cell, &beta);
    println!("U_eq {ueq:.8}, ⅓ tr U_cart {:.8}", trace / 3.0);
    assert!((ueq - trace / 3.0).abs() < 1e-14);

    // PbSO₄ (FP_tests/pbso4_cheb.pcr, every atom anisotropic): U_eq with its su in the results file, not 0.
    let d = root("../manuals/fullprof/FP_tests");
    let p = drs::pcr::read(&latin1(&d.join("pbso4_cheb.pcr"))).unwrap();
    let c = Control::from_pcr(&p, &d, "pbso4_cheb", &d).unwrap();
    let data = Model::from_pcr(&p, &d, "pbso4_cheb").unwrap().data.unwrap();
    let mut m = c.model.clone();
    let set = ParameterSet::from_pcr(&p, &m).unwrap();
    let r = drs::lsq::refine(&mut m, &data, &set, drs::lsq::Options { max_cycles: 2, eps: 0.1 }, |_| {});
    let control = Control { model: m.clone(), ..c };
    let (text, _) = drs::results::write(&control, &data, Some(&set), Some(&r), "pbso4_cheb 2026-10-12 12.00.00").unwrap();
    let prop = drs::geometry::Propagator::new(&m, &set, &r);
    let ueq = drs::geometry::u_eq(&m, 0, Some(&prop));
    let Intensities::Structure(st, _) = &m.phases[0].intensities else { panic!() };
    for (a, e) in st.atoms.iter().zip(&ueq) {
        let (v, s) = e.expect("anisotropic");
        // Orthorhombic: U_eq is the mean of U11, U22, U33.
        let b = a.beta.unwrap();
        let mean = (0..3).map(|k| b[k] / (2.0 * pi2 * [m.phases[0].cell.reciprocal().a, m.phases[0].cell.reciprocal().b, m.phases[0].cell.reciprocal().c][k].powi(2))).sum::<f64>() / 3.0;
        assert!((v - mean).abs() < 1e-15 && s > 0.0, "{}: {v} {s} against {mean}", a.label);
        let line = text.lines().find(|l| l.starts_with(&format!("{} ", a.label)) && l.contains(" Uani ")).unwrap();
        assert!(line.contains(&drs::cif_write::with_su(v, s, 5)), "{line}");
        println!("{line}");
    }
}

#[test]
fn p17_parent_cell_of_a_coordinate_phase() {
    let (c, data) = case_a();
    let modes_cell = drs::geometry::parent_cell(&c.model, 0, None).unwrap();
    // The same structure as coordinates (the positions the modes give), the parent given to the phase.
    let mut m = c.model.clone();
    let parent = m.phases[0].modes.take().unwrap().parent;
    m.phases[0].parent = parent;
    let coords = drs::geometry::parent_cell(&m, 0, None).unwrap();
    for k in 0..6 {
        assert!((coords[k].0 - modes_cell[k].0).abs() < 1e-12, "{k}: {:?} {:?}", coords[k], modes_cell[k]);
    }
    println!("parent cell {:?}", coords.map(|x| x.0));
    // Kept by the control file.
    let control = Control { model: m.clone(), refinement: Default::default(), ..c.clone() };
    let text = control.to_text();
    assert!(text.contains("mode_parent_group") && !text.contains("[[phase.mode]]"));
    let back = Control::from_text(&text).unwrap();
    assert_eq!(back.model.phases[0].parent, m.phases[0].parent);
    // Refined by its cell: the reduced cell with su's in the results file.
    let set = ParameterSet::from_addresses(&["phase1.scale_factor", "phase1.cell_a", "phase1.cell_b", "phase1.cell_c", "phase1.cell_beta"].map(String::from), &[], &[], &m).unwrap();
    let r = drs::lsq::refine(&mut m, &data, &set, drs::lsq::Options { max_cycles: 2, eps: 0.1 }, |_| {});
    let control = Control { model: m.clone(), refinement: Default::default(), ..c };
    let (text, _) = drs::results::write(&control, &data, Some(&set), Some(&r), "a_coordinates 2026-10-12 12.00.00").unwrap();
    let line = text.lines().find(|l| l.starts_with("_drs_parent_cell_a")).expect("the parent cell");
    println!("{line}");
    assert!(line.contains('('), "an su: {line}");
}

#[test]
fn p18_valence_sums_in_the_results_file() {
    let (mut c, data) = case_a();
    // Charges: Sr²⁺, Sc³⁺, Sb⁵⁺, O²⁻ (P5).
    if let Intensities::Structure(st, _) = &mut c.model.phases[0].intensities {
        for a in st.atoms.iter_mut() {
            a.charge = Some(match a.species.to_uppercase().chars().take(2).collect::<String>().as_str() {
                "SR" => 2,
                "SC" => 3,
                "SB" => 5,
                _ => -2,
            });
        }
    }
    let mut m = c.model.clone();
    let p = drs::pcr::read(&latin1(&root("../manuals/fullprof/runs/modes/A_modes/a_modes.pcr"))).unwrap();
    let set = ParameterSet::from_pcr(&p, &m).unwrap();
    let r = drs::lsq::refine(&mut m, &data, &set, drs::lsq::Options { max_cycles: 1, eps: 0.1 }, |_| {});
    let prop = drs::geometry::Propagator::new(&m, &set, &r);
    let vs = drs::geometry::valence(&m, 0, Some(&prop), 3.0, 0.2).unwrap().expect("every atom has a charge");
    // The values are bonds::valence_sums's on the same contacts.
    let (st, contacts, _) = drs::geometry::contacts_with_su(&m, 0, None, 3.0, 0.0).unwrap();
    let occ = drs::bonds::bond_str_occupancies(&st).unwrap();
    let (sums, _) = drs::bonds::valence_sums(&st, &contacts, &occ, 0.2).unwrap();
    assert_eq!(vs.len(), sums.len());
    for (v, s) in vs.iter().zip(&sums) {
        println!("{:6} CN {:2} ⟨d⟩ {:.4}({:.4}) BVS {:.4}({:.4}); independent distances: ({:.4})", st.atoms[v.atom].label, v.coordination, v.mean_distance.0, v.mean_distance.1, v.sum.0, v.sum.1, s.sum.1);
        assert_eq!((v.atom, v.coordination), (s.atom, s.coordination));
        assert!((v.mean_distance.0 - s.mean_distance.0).abs() < 1e-12 && (v.sum.0 - s.sum.0).abs() < 1e-12);
        assert!(v.sum.1 > 0.0 && v.sum.1.is_finite() && v.mean_distance.1 > 0.0);
    }
    let control = Control { model: m.clone(), ..c };
    let (text, _) = drs::results::write(&control, &data, Some(&set), Some(&r), "a_modes 2026-10-12 12.00.00").unwrap();
    let rows: Vec<&str> = text.lines().skip_while(|l| !l.starts_with("_drs_valence_sum")).skip(1).take_while(|l| !l.trim().is_empty() && !l.starts_with('_') && !l.starts_with("loop_") && !l.starts_with("data_")).collect();
    assert_eq!(rows.len(), vs.len(), "{rows:?}");
    let labels = drs::geometry::unique_labels(&st);
    for (row, v) in rows.iter().zip(&vs) {
        assert!(row.trim_start_matches('\'').starts_with(&labels[v.atom]) && row.contains(&drs::cif_write::with_su(v.sum.0, v.sum.1, 3)), "{row}");
    }
}

#[test]
fn app31_a_file_named_by_its_place() {
    use drs::formats::project::{describe, describe_in, Description, Role};
    use drs::formats::Probe;
    let text = std::fs::read(root("../publish/examples/Results/y2o3_atoms 2026-10-11 12.00.00 results.cif")).unwrap();
    let rel = "Results/y2o3_atoms 2026-10-11 12.00.00 results.cif";
    let p = Probe::new(rel, &text);
    // By its content it was a powder CIF the registry reads as a pattern (another results file, a control file: it
    // holds one); its first line now names it, and by its place it is the results of a refinement (L199).
    assert!(matches!(describe(&p), Description::Project(ref k) if k.id == "results"));
    match describe_in(rel, &p) {
        Description::Project(k) => assert_eq!((k.id, k.role), ("results", Role::Output)),
        other => panic!("{other:?}"),
    }
    // Measured data keep the content's description, whatever their name.
    let dat = std::fs::read(root("../manuals/fullprof/runs/Y2O3/y2o3.dat")).unwrap();
    assert!(matches!(describe_in("Data/y2o3.dat", &Probe::new("Data/y2o3.dat", &dat)), Description::Data { .. }));
    // A file outside the dictionary's places, by its content.
    assert_eq!(describe_in("notes/y2o3.dat", &Probe::new("notes/y2o3.dat", &dat)), describe(&Probe::new("notes/y2o3.dat", &dat)));
    // The command line.
    let d = root("target/test-open-requests/app31");
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(d.join("Results")).unwrap();
    std::fs::write(d.join(rel), &text).unwrap();
    let o = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(["data", "info", d.join(rel).to_str().unwrap(), "--in", d.to_str().unwrap()]).output().unwrap();
    let out = String::from_utf8_lossy(&o.stdout);
    assert!(o.status.success() && out.contains("(results; Output)"), "{out}{}", String::from_utf8_lossy(&o.stderr));
}

#[test]
fn app32_fullprof_exports_in_exports_fullprof() {
    use drs::formats::files::{kind_of, path};
    for id in ["fullprof_pcr_export", "fullprof_background_export", "fullprof_reflections_export", "fullprof_resolution_export", "fullprof_data_export"] {
        let p = path(id, &[("control", "y2o3"), ("n", "1")]).unwrap();
        assert!(p.starts_with("Exports/FullProf/"), "{p}");
        assert_eq!(kind_of(&p).unwrap().id, id);
    }
    // FullProf/ holds FullProf's original files: the dictionary names none of them an export.
    for f in ["FullProf/y2o3.pcr", "FullProf/y2o3.bac", "FullProf/y2o31.hkl", "FullProf/y2o3.irf"] {
        assert!(kind_of(f).is_none_or(|k| !k.id.ends_with("_export")), "{f}");
    }
    // drs convert writes there.
    let d = root("target/test-open-requests/app32");
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    let run = root("../manuals/fullprof/runs/Y2O3");
    for f in ["y2o3.dat", "y2o3_structure.pcr"] {
        std::fs::copy(run.join(f), d.join(f)).unwrap();
    }
    let drs_bin = env!("CARGO_BIN_EXE_drs");
    let o = std::process::Command::new(drs_bin).args(["convert", d.join("y2o3_structure.pcr").to_str().unwrap()]).output().unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let o = std::process::Command::new(drs_bin).args(["convert", d.join("y2o3_structure.drc").to_str().unwrap()]).output().unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(d.join("Exports/FullProf/y2o3_structure.pcr").exists());
}

#[test]
fn e7_on_the_command_line() {
    let d = root("../manuals/fullprof/runs/modes/A_modes");
    let pcr = d.join("a_modes.pcr");
    let o = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(["modes", pcr.to_str().unwrap(), "--amplimodes", root(AMPLIMODES_A).to_str().unwrap(), "--equivalent", "--reference", pcr.to_str().unwrap()]).output().unwrap();
    let out = String::from_utf8_lossy(&o.stdout);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    // The 4 sets of `explore_requests` (24 before P21's fix, L218), and the reference's own amplitudes as the nearest
    // (the identity).
    assert!(out.contains("4 amplitude sets describe this structure"), "{out}");
    assert!(out.contains("nearest") && out.lines().any(|l| l.contains("nearest") && l.contains("x,y,z [0, 0, 0]")), "{out}");
    assert!(out.lines().any(|l| l.trim_start().starts_with("A3_GM4+") && l.contains("-0.828624")), "{out}");
}
