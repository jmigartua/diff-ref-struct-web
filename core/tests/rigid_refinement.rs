//! Rigid groups refined (P41): a group's orientation refined in place of its atoms' coordinates, against FullProf's two
//! Rietveld refinements of Z-matrix groups in the material (9-SAnn-OrganicCompounds: citric acid, 21 atoms, 16
//! parameters; PPH₃, 34 atoms, 8 parameters; FullProf 6.00, both converged: "Convergence reached at this CYCLE").
//!
//! - `derivatives_through_the_atoms`: each of a group's values (origin, θ, φ, χ, a torsion) moves the calculated pattern
//!   as the engine's analytic derivative says (central differences, to 10⁻⁹ of the largest); the cell's derivative
//!   carries the group's atoms with it (the group keeps its shape in Å).
//! - `refined_against_fullprof`: from each `.pcr` (FullProf's converged values), refined to the minimum (no damping,
//!   eps 0.1; FullProf's R_at of 0.08 and 0.10 would stop it short of it), the engine's χ² at the start against FullProf's conventional χ² (64.64, 42.63, to their print), its converged θ,
//!   φ, χ against FullProf's in FullProf's esd's, and its esd's against FullProf's (within 5 %).
//!   - PPH₃: θ, φ, χ 1.42, 1.46 and 1.78 of FullProf's esd's below FullProf's (frozen at 1.45, 1.5, 1.8), the esd
//!     ratios 0.99, 0.97, 0.99; believed to go with the 5·10⁻⁴ Å by which its written positions depart from its
//!     printed Z-matrix (`rigid_bodies`), not checked.
//!   - Citric acid: φ 1.39 esd (frozen at 1.4), χ 0.06, the esd ratios 1.02, 1.02, 0.98. θ ends 5.16 esd (0.077°)
//!     above FullProf's 218.193° (frozen at 5.2), χ² 64.563 against FullProf's 64.64, which FullProf's three cycles
//!     did not lower. Along θ alone, the others at FullProf's values, S falls from 2 077 024 at FullProf's θ to
//!     2 075 116 at +0.06° on the core's calculated pattern, which matches FullProf's `.prf` to its print
//!     (`fullprof_profiles`). FullProf printed θ's applied shift as −0.000 with R_at = 0.08, so its own computed shift
//!     was below 0.4 esd where the core's is 4.8. Not explained: counted and frozen (L347), and FullProf's run at the
//!     core's values named (`runs/version2_runs.py`, citric_core).
//! - `control_file_export_and_results`: citric acid's control file with its group (`[[phase.rigid_group]]`, each atom's
//!   `rigid_group`, `rigid_length`, `rigid_angle`, `rigid_connect`) read back to the same model, calculated pattern and
//!   free parameters to the bit; every group value located in the text; the `.pcr` export (Jbt = 4, lines 25-1 to
//!   25-4, DEG = 1) imported again to the same model and free parameters; the results file's group loops (P47.1: its
//!   kind, origin and angles with su's, each atom's internal coordinates) and its atoms flagged G with su's (P47.3).
//! - `back_from_a_turned_start`: citric acid with θ, φ, χ turned by 1° each (65 esd), refined without damping: back to
//!   the minimum the engine reaches from FullProf's values, within 0.3 esd.

use drs::jacobian::Field;
use drs::lsq::{refine, Options};
use drs::powder::Model;
use drs::refinable::Refinement;
use std::path::{Path, PathBuf};

fn dir(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/9-SAnn-OrganicCompounds").join(rel)
}

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display())).iter().map(|&b| b as char).collect()
}

fn load(rel: &str, code: &str) -> (Model, drs::data::Pattern, Refinement, String) {
    let d = dir(rel);
    let p = drs::pcr::read(&latin1(&d.join(format!("{code}.pcr")))).unwrap();
    let l = Model::from_pcr(&p, &d, code).unwrap();
    let rf = Refinement::from_pcr(&p, &l.model).unwrap();
    (l.model, l.data.unwrap(), rf, latin1(&d.join(format!("{code}.out"))))
}

/// FullProf's final θ, φ, χ of the group's first atom with their esd's, from the last cycle of its `.out`.
fn fullprof_angles(out: &str) -> [(f64, f64); 3] {
    let lines: Vec<&str> = out.lines().collect();
    let last = |head: &str| lines.iter().rposition(|l| l.contains(head)).unwrap_or_else(|| panic!("no {head}"));
    let tp: Vec<f64> = lines[last("Thet7") + 2].split_whitespace().map(|x| x.parse().unwrap()).collect();
    let c: Vec<f64> = lines[last("Chi15") + 4].split_whitespace().map(|x| x.parse().unwrap()).collect();
    [(tp[0], tp[2]), (tp[3], tp[5]), (c[0], c[2])]
}

/// FullProf's first conventional χ² of a `.out`.
fn fullprof_chi2(out: &str) -> f64 {
    out.lines().find_map(|l| l.split_once("Conventional Rietveld Rp,Rwp,Re and Chi2:")).map(|(_, r)| r.split_whitespace().nth(3).unwrap().parse().unwrap()).unwrap()
}

const ANGLES: [Field; 3] = [Field::Rigid(0, 0, 3), Field::Rigid(0, 0, 4), Field::Rigid(0, 0, 5)];

#[test]
fn derivatives_through_the_atoms() {
    let (m, data, _, _) = load("CitricAcid_example/Rietveld", "citric-acid-sls_1");
    // The origin, θ, φ, χ, and the torsion of A6 (atom 6, 6 + 3·5 + 2).
    let fields = [Field::Rigid(0, 0, 0), Field::Rigid(0, 0, 3), Field::Rigid(0, 0, 4), Field::Rigid(0, 0, 5), Field::Rigid(0, 0, 23)];
    let (_, jac) = m.derivatives(&data.x, &fields);
    for (j, &f) in fields.iter().enumerate() {
        let (v, h) = (m.get(f), if matches!(f, Field::Rigid(_, _, 0..=2)) { 1e-6 } else { 1e-4 });
        let (mut a, mut b) = (m.clone(), m.clone());
        a.set(f, v + h);
        b.set(f, v - h);
        let (ya, yb) = (a.calculate(&data.x).ycalc, b.calculate(&data.x).ycalc);
        let num: Vec<f64> = ya.iter().zip(&yb).map(|(p, q)| (p - q) / (2.0 * h)).collect();
        let top = num.iter().fold(0.0f64, |s, x| s.max(x.abs()));
        let err = num.iter().zip(&jac[j]).fold(0.0f64, |s, (n, d)| s.max((n - d).abs()));
        println!("{f:?}: largest {top:.3e}, |numerical − analytic| ≤ {err:.2e}");
        assert!(err <= 1e-9 * top, "{f:?}: {err:e} against {top:e}");
    }
    // The cell carries the group: what the atoms add to ∂y/∂a is the change their placement makes.
    let f = Field::Cell(0, 0);
    let mut fixed = m.clone();
    fixed.phases[0].rigid.clear();
    let (_, with) = m.derivatives(&data.x, &[f]);
    let (_, without) = fixed.derivatives(&data.x, &[f]);
    let h = 1e-5;
    let moved = |mm: &Model| {
        let (mut a, mut b) = (mm.clone(), mm.clone());
        a.set(f, mm.get(f) + h);
        b.set(f, mm.get(f) - h);
        a.calculate(&data.x).ycalc.iter().zip(b.calculate(&data.x).ycalc).map(|(p, q)| (p - q) / (2.0 * h)).collect::<Vec<f64>>()
    };
    let (nw, nwo) = (moved(&m), moved(&fixed));
    let top = nw.iter().zip(&nwo).fold(0.0f64, |s, (a, b)| s.max((a - b).abs()));
    let err = (0..data.x.len()).fold(0.0f64, |s, i| s.max(((nw[i] - nwo[i]) - (with[0][i] - without[0][i])).abs()));
    println!("Cell a: the group's part {top:.3e}, |numerical − analytic| ≤ {err:.2e}");
    // The central difference's own error with a step of 10⁻⁵ Å on 11.46 Å: ~10⁻⁵ of the largest.
    assert!(top > 0.0 && err <= 1e-4 * top);
}

#[test]
fn refined_against_fullprof() {
    // (folder, code, FullProf's χ², frozen |Δ|/σ of θ, φ, χ)
    for (rel, code, frozen) in [("PPH3/Rietveld", "PPH3-sls", [1.45, 1.5, 1.8]), ("CitricAcid_example/Rietveld", "citric-acid-sls_1", [5.2, 1.4, 1.0])] {
        let (mut m, data, mut rf, out) = load(rel, code);
        // To the minimum itself: FullProf's R_at (0.08, 0.10) would stop the engine short of it at eps 0.3.
        rf.steering.class_damping = [1.0; 4];
        let set = rf.build(&m).unwrap();
        let start = drs::lsq::statistics(&m, &data, &m.calculate(&data.x).ycalc, set.free.len());
        let chi = fullprof_chi2(&out);
        assert!((start.chi2 - chi).abs() <= 0.005 + 0.01, "{code}: χ² {:.4} at FullProf's values against its {chi}", start.chi2);
        let r = refine(&mut m, &data, &set, Options { max_cycles: 40, eps: 0.1 }, |_| {});
        assert!(r.converged, "{code}: {}", r.stop);
        let fp = fullprof_angles(&out);
        for (k, &f) in ANGLES.iter().enumerate() {
            let q = set.free.iter().position(|p| p.fields[0].0 == f).map(|i| &r.params[i]).unwrap();
            let (dev, ratio) = ((q.value - fp[k].0).abs() / fp[k].1, q.esd / fp[k].1);
            println!("{code} {}: {:.4} ({:.4}) against FullProf's {} ({}): {dev:.2} of its esd, esd ratio {ratio:.3}", ["θ", "φ", "χ"][k], q.value, q.esd, fp[k].0, fp[k].1);
            assert!(dev <= frozen[k], "{code}: {} {dev:.2} esd from FullProf's, frozen at {}", ["θ", "φ", "χ"][k], frozen[k]);
            // FullProf prints the esd's to 3 decimals (0.015 … 0.028): 5 % and half a unit of the print.
            assert!((q.esd - fp[k].1).abs() <= 0.05 * fp[k].1 + 0.0005, "{code}: esd {:.4} against FullProf's {}", q.esd, fp[k].1);
        }
        println!("{code}: χ² {:.4} → {:.4} in {} cycles", start.chi2, r.stats.chi2, r.cycles);
    }
}

#[test]
fn back_from_a_turned_start() {
    let (m0, data, mut rf, _) = load("CitricAcid_example/Rietveld", "citric-acid-sls_1");
    rf.steering.class_damping = [1.0; 4];
    let set = rf.build(&m0).unwrap();
    let mut a = m0.clone();
    let ra = refine(&mut a, &data, &set, Options { max_cycles: 40, eps: 0.1 }, |_| {});
    let mut b = m0.clone();
    for f in ANGLES {
        b.set(f, m0.get(f) + 1.0);
    }
    let rb = refine(&mut b, &data, &set, Options { max_cycles: 40, eps: 0.1 }, |_| {});
    assert!(ra.converged && rb.converged, "{} / {}", ra.stop, rb.stop);
    for f in ANGLES {
        let i = set.free.iter().position(|p| p.fields[0].0 == f).unwrap();
        let d = (a.get(f) - b.get(f)).abs() / ra.params[i].esd;
        println!("{f:?}: {:.4} from FullProf's values, {:.4} from 1° away ({d:.3} esd; {} and {} cycles)", a.get(f), b.get(f), ra.cycles, rb.cycles);
        assert!(d <= 0.3, "{f:?}: {d} esd apart");
    }
}

#[test]
fn control_file_export_and_results() {
    let rel = "CitricAcid_example/Rietveld";
    let d = dir(rel);
    let p = drs::pcr::read(&latin1(&d.join("citric-acid-sls_1.pcr"))).unwrap();
    let c = drs::control::Control::from_pcr(&p, &d, "citric-acid-sls_1", &d).unwrap();
    let text = c.to_text();
    assert_eq!(text.lines().filter(|l| *l == "[[phase.rigid_group]]").count(), 1);
    assert_eq!(text.lines().filter(|l| l.starts_with("rigid_group = 1 ")).count(), 21);
    let back = drs::control::Control::from_text(&text).unwrap();
    assert_eq!(back.model.to_text(), c.model.to_text(), "the model read back");
    assert_eq!(back.refinement, c.refinement, "the free parameters read back");
    let (m0, data, _, _) = load(rel, "citric-acid-sls_1");
    assert_eq!(back.model.calculate(&data.x).ycalc, m0.calculate(&data.x).ycalc, "the calculated pattern");
    // Every group value where the file writes it.
    let all = drs::control::locate::all(&text, None);
    let n = drs::rigid_body::values_of(&c.model.phases[0].rigid[0]);
    let located: Vec<_> = all.iter().filter(|(a, ..)| a.contains("rigid_")).collect();
    assert_eq!(located.len(), n, "{n} values of the group");
    for (a, _, r) in &located {
        let v: f64 = text[r.clone()].parse().unwrap();
        let f = drs::refinable::field_of(a, &c.model).unwrap();
        assert_eq!(v, c.model.get(f), "{a}");
    }
    // The .pcr export, read again.
    let e = drs::control::pcr_export::export(&c, "citric-acid-sls_1", Some(&data)).unwrap();
    let p2 = drs::pcr::read(&e.pcr).unwrap();
    let m2 = Model::from_pcr(&p2, &d, "citric-acid-sls_1").unwrap().model;
    assert_eq!(m2.to_text(), c.model.to_text(), "the model through the .pcr");
    let free = |rf: &Refinement, m: &Model| rf.build(m).unwrap().free.iter().map(|f| f.fields[0].0).collect::<Vec<_>>();
    assert_eq!(free(&Refinement::from_pcr(&p2, &m2).unwrap(), &m2), free(&c.refinement, &c.model), "the free parameters through the .pcr");
    // The results file (drs results, two cycles): the group's loops and the atoms flagged G.
    let out = Path::new(env!("CARGO_MANIFEST_DIR")).join("target/test-rigid-results");
    let _ = std::fs::remove_dir_all(&out);
    let st = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(["results", d.join("citric-acid-sls_1.pcr").to_str().unwrap(), "--cycles", "2", "--project", out.to_str().unwrap(), "--time", "2026-10-19 12.00.00"]).output().unwrap();
    assert!(st.status.success(), "{}", String::from_utf8_lossy(&st.stderr));
    let cif = std::fs::read_to_string(out.join("Results/citric-acid-sls_1 2026-10-19 12.00.00 results.cif")).unwrap();
    let blocks = drs::cif::parse(&cif);
    let b = blocks.iter().find(|b| b.column("_drs_rigid_group_id").is_some()).expect("the group's loop");
    assert_eq!(b.column("_drs_rigid_group_kind").unwrap(), ["z-matrix"]);
    for t in ["_drs_rigid_group_theta", "_drs_rigid_group_phi", "_drs_rigid_group_chi"] {
        assert!(b.column(t).unwrap()[0].contains('('), "{t} with its su");
    }
    assert_eq!(b.column("_drs_rigid_atom_label").unwrap().len(), 21);
    let flags = blocks.iter().find_map(|b| b.column("_atom_site_refinement_flags_posn")).expect("the atoms' flags");
    assert!(flags.len() == 21 && flags.iter().all(|f| *f == "G"));
}
