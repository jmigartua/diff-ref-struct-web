//! The small requests of L170's first step (L171):
//! - Publish P20: B with its su, `b_from_u` and `u_from_b` exact inverses; the results file's B_iso for each atom
//!   (B_eq for an anisotropic one) with its su, the refined B's su that of the engine's;
//! - app 33: `control::open_in` opens a version in `Versions/` with its files named from the project's folder, and
//!   `--base` does it on the command line;
//! - Publish P19: case A's structure as coordinates, described by its own mode file, gives the mode phase's amplitudes
//!   (residual below 10⁻⁹ Å) and the same scene (reference, displacements by irrep); every image's reference position
//!   in fractions, its Cartesian position the scene's reference; `drs structure --modes`.

use drs::control::Control;
use drs::powder::{Intensities, Model};
use drs::refinable::ParameterSet;
use std::path::{Path, PathBuf};

mod common;
use common::latin1;

fn root(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)
}

#[test]
fn p20_b_with_its_su() {
    let (b, sb) = drs::structure_factor::b_from_u(0.0125, 0.0004);
    let (u, su) = drs::structure_factor::u_from_b(b, sb);
    assert!((u - 0.0125).abs() < 1e-17 && (su - 0.0004).abs() < 1e-18);
    assert!((b - 8.0 * std::f64::consts::PI.powi(2) * 0.0125).abs() < 1e-15);
    // Y2O3's B's refined (the atoms' B's freed): the results file's B and su against the engine's.
    let d = root("../manuals/fullprof/runs/Y2O3");
    let p = drs::pcr::read(&latin1(&d.join("y2o3_structure.pcr"))).unwrap();
    let c = Control::from_pcr(&p, &d, "y2o3_structure", &d).unwrap();
    let data = Model::from_pcr(&p, &d, "y2o3_structure").unwrap().data.unwrap();
    let mut m = c.model.clone();
    let free: Vec<String> = ["phase1.scale_factor", "phase1.Y.atom_biso", "phase1.Y:2.atom_biso", "phase1.O.atom_biso"].map(String::from).to_vec();
    let set = ParameterSet::from_addresses(&free, &[], &[], &m).unwrap_or_else(|e| panic!("{e}: {:?}", m.fields().iter().filter(|f| f.address.contains("biso")).map(|f| &f.address).collect::<Vec<_>>()));
    let r = drs::lsq::refine(&mut m, &data, &set, drs::lsq::Options { max_cycles: 3, eps: 0.1 }, |_| {});
    let control = Control { model: m.clone(), refinement: Default::default(), ..c };
    let (text, _) = drs::results::write(&control, &data, Some(&set), Some(&r), "y2o3 2026-10-13 12.00.00").unwrap();
    let rows: Vec<&str> = text.lines().skip_while(|l| !l.starts_with("_drs_atom_site_B_iso_or_equiv")).skip(1).take(3).collect();
    let Intensities::Structure(st, _) = &m.phases[0].intensities else { panic!() };
    let labels = drs::geometry::unique_labels(st);
    for (a, (row, label)) in rows.iter().zip(&labels).enumerate() {
        let addr = drs::refinable::address(drs::jacobian::Field::AtomB(0, a), &m);
        let (_, v, e) = r.fields.iter().find(|x| x.0 == addr).cloned().unwrap();
        println!("{row}   (engine {addr} {v} ± {e})");
        assert_eq!(*row, format!("'{label}' {}", drs::cif_write::with_su(v, e, 4)));
    }
}

#[test]
fn app33_a_version_opened_from_the_projects_folder() {
    let d = root("target/test-small-requests/app33");
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(d.join("Versions")).unwrap();
    let run = root("../manuals/fullprof/runs/Y2O3");
    std::fs::copy(run.join("y2o3.dat"), d.join("y2o3.dat")).unwrap();
    let p = drs::pcr::read(&latin1(&run.join("y2o3_structure.pcr"))).unwrap();
    let c = Control::from_pcr(&p, &d.join("x"), "y2o3_structure", &d).unwrap();
    let mut c = c;
    c.data.as_mut().unwrap().file = "y2o3.dat".into();
    let version = d.join("Versions/y2o3 2026-10-13 12.00.00.drc");
    std::fs::write(&version, c.to_text()).unwrap();
    assert!(drs::control::open(&version).is_err(), "beside the version there is no data");
    let o = drs::control::open_in(&version, &d).unwrap();
    assert_eq!(o.data.unwrap().len(), 13091);
    assert_eq!(o.data_path.unwrap(), d.join("y2o3.dat"));
    // --base on the command line.
    let run = |args: &[&str]| std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(args).output().unwrap();
    let o = run(&["calc", version.to_str().unwrap()]);
    assert!(!o.status.success());
    let o = run(&["calc", version.to_str().unwrap(), "--base", d.to_str().unwrap()]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(String::from_utf8_lossy(&o.stdout).contains("Rp "));
}

const AMPLIMODES_A: &str = "../manuals/bcs/mode-crystallography-manual/course-um6p-2024/04_tools_materials/0_1_P21nToCubicTransData.txt";

#[test]
fn p19_a_coordinate_phase_described_by_modes() {
    let d = root("../manuals/fullprof/runs/modes/A_modes");
    let p = drs::pcr::read(&latin1(&d.join("a_modes.pcr"))).unwrap();
    let mut modes = Model::from_pcr(&p, &d, "a_modes").unwrap().model;
    let a = drs::modes::read_amplimodes(&latin1(&root(AMPLIMODES_A)));
    modes.phases[0].modes.as_mut().unwrap().parent = a[0].mode_parent();
    let metric = a[0].reference_cell();
    // The same structure as coordinates: the positions the modes give, no modes; its parent given (P17).
    let mut coords = modes.clone();
    let ms = coords.phases[0].modes.take().unwrap();
    coords.phases[0].parent = ms.parent.clone();
    let (dm, residual) = drs::geometry::decompose(&coords, 0, &modes, 0, metric, [0.0; 3]).unwrap();
    println!("residual {residual:.2e} Å");
    assert!(residual < 1e-9);
    let (q0, q1): (Vec<f64>, Vec<f64>) = (ms.modes.iter().map(|m| m.amplitude).collect(), dm.phases[0].modes.as_ref().unwrap().modes.iter().map(|m| m.amplitude).collect());
    for (x, y) in q0.iter().zip(&q1) {
        assert!((x - y).abs() < 1e-6, "{q0:?} {q1:?}");
    }
    let rule = drs::geometry::SceneRule::Ionic(1.2);
    let (a0, b0) = ([0.0; 3], [1.0; 3]);
    let s0 = drs::geometry::scene(&modes, 0, a0, b0, &rule).unwrap();
    let s1 = drs::geometry::scene(&dm, 0, a0, b0, &rule).unwrap();
    assert_eq!(s0.images.len(), s1.images.len());
    let close = |u: &[[f64; 3]], v: &[[f64; 3]]| u.iter().zip(v).all(|(x, y)| (0..3).all(|k| (x[k] - y[k]).abs() < 1e-6));
    assert!(close(&s0.reference, &s1.reference));
    assert_eq!(s0.displacements.len(), s1.displacements.len());
    for ((n0, u0), (n1, u1)) in s0.displacements.iter().zip(&s1.displacements) {
        assert!(n0 == n1 && close(u0, u1), "{n0}");
    }
    // The reference in fractions: the Cartesian reference of each image.
    let m3 = modes.phases[0].cell.crystal_to_orthonormal();
    for (f, c) in s0.reference_fract.iter().zip(&s0.reference) {
        let x: [f64; 3] = std::array::from_fn(|i| (0..3).map(|j| m3[i][j] * f[j]).sum());
        assert!((0..3).all(|k| (x[k] - c[k]).abs() < 1e-9));
    }
    // An image whose reference sits inside a slab its refined position leaves (P19's Sr₂ScSbO₆ case): told apart.
    let moved = s0.images.iter().zip(&s0.reference_fract).filter(|(im, r)| (0..3).any(|k| (im.fract[k] - r[k]).abs() > 0.01)).count();
    println!("{} images, {moved} more than 0.01 from their reference in some fraction", s0.images.len());
    assert!(moved > 0);
    // The command line, on a control file of the coordinates.
    let dir = root("target/test-small-requests/p19");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::copy(d.join("SrScLRRT.dat"), dir.join("SrScLRRT.dat")).unwrap();
    let mut c = Control::from_pcr(&p, &d, "a_modes", &dir).unwrap();
    c.model = coords;
    c.refinement = Default::default();
    std::fs::write(dir.join("coords.drc"), c.to_text()).unwrap();
    let o = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(["structure", dir.join("coords.drc").to_str().unwrap(), "--cartesian", "--modes", d.join("a_modes.pcr").to_str().unwrap(), "--amplimodes", root(AMPLIMODES_A).to_str().unwrap()]).output().unwrap();
    let out = String::from_utf8_lossy(&o.stdout);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(out.contains("the part they do not describe") && out.lines().any(|l| l.starts_with("reference ")) && out.lines().any(|l| l.starts_with("displacement ")), "{out}");
}
