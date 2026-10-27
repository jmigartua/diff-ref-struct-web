//! A joint refinement's settings and its control files (P06; the app's request 43, L329).
//!
//! - `control_files_round_trip`: scor_rxn's joint refinement (two phases in two patterns) written as its control files
//!   (the joint `.drc` and one per pattern, `Joint::save`) and read back (`Joint::open`): the same patterns, weights,
//!   flags, phases and joint parameters, the same models to the bit, the same objective; written again, the same
//!   text; refined two cycles, the same values as the joint refinement read from the `.pcr`. Each pattern's file opens
//!   on its own as a control file of one pattern.
//! - `weights_and_flags`: `set_weight` and `set_weights` normalise to sum 1 and the objective follows, Q = Σ W_k S_k;
//!   a pattern calculated only leaves Q; the last refined pattern cannot be switched off.
//! - `contributions`: scor_rxn's Sr₂Co₂O₅ taken out of its X-ray pattern (that pattern's Ycal the old one less the
//!   phase's curve; the phase's parameters there gone, its structure's now pattern 1's alone) and given back (the
//!   phase as pattern 1 has it, its atoms' scattering the X-rays' again, the same as the X-ray pattern had; its atoms'
//!   parameters tied to pattern 1's again, its scale, profile and cell the pattern's own).

use drs::joint::{Esd, Joint};
use drs::lsq::Options;
use std::path::{Path, PathBuf};

fn root(rel: &str) -> String {
    format!("{}/../{rel}", env!("CARGO_MANIFEST_DIR"))
}

fn scor() -> Joint {
    let dir = root("manuals/fullprof/runs/scor_rxn");
    let text: String = std::fs::read(format!("{dir}/scor_rxn.pcr")).unwrap().iter().map(|&b| b as char).collect();
    Joint::from_pcr(&text, Path::new(&dir), "scor_rxn", false).unwrap()
}

fn temp(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("drs_joint_control_{name}_{}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn shape(j: &Joint) -> Vec<(String, Vec<(usize, String)>)> {
    j.parameters.iter().map(|q| (q.name.clone(), q.members.iter().map(|&(k, i)| (k, j.patterns[k].set.free[i].name.clone())).collect())).collect()
}

#[test]
fn control_files_round_trip() {
    let j = scor();
    let d = temp("rt");
    let written = j.save(&d.join("scor.drc")).unwrap();
    assert_eq!(written.len(), 3);
    let back = Joint::open(&d.join("scor.drc")).unwrap();
    assert_eq!(back.patterns.len(), 2);
    for (a, b) in j.patterns.iter().zip(&back.patterns) {
        assert_eq!((a.weight, a.refined, &a.phases, a.data.x.len()), (b.weight, b.refined, &b.phases, b.data.x.len()));
        assert_eq!(a.model.to_text(), b.model.to_text(), "the model read back");
        assert_eq!(a.set.free.iter().map(|f| &f.name).collect::<Vec<_>>(), b.set.free.iter().map(|f| &f.name).collect::<Vec<_>>());
    }
    assert_eq!(shape(&j), shape(&back), "the joint parameters");
    assert_eq!(j.objective(), back.objective());
    let again = back.control_files("scor", &d).unwrap();
    let first = j.control_files("scor", &d).unwrap();
    assert_eq!(first, again, "written again, the same text");
    // Each pattern's file is a control file of one pattern.
    for k in 1..=2 {
        let o = drs::control::open(&d.join(format!("scor pattern {k}.drc"))).unwrap();
        assert_eq!(o.control.model.to_text(), back.patterns[k - 1].model.to_text());
    }
    // Refined, the same as from the .pcr.
    let (mut a, mut b) = (scor(), back);
    let opts = Options { max_cycles: 2, eps: 1e-9 };
    let (ra, rb) = (a.refine(opts, Esd::Default, |_| true), b.refine(opts, Esd::Default, |_| true));
    for (x, y) in ra.params.iter().zip(&rb.params) {
        assert_eq!(x.0, y.0);
        assert!((x.1 - y.1).abs() <= 1e-9 * x.1.abs().max(1e-6) && (x.2 - y.2).abs() <= 1e-6 * x.2.max(1e-12), "{}: {} ({}) against {} ({})", x.0, x.1, x.2, y.1, y.2);
    }
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn weights_and_flags() {
    let mut j = scor();
    let s: Vec<f64> = j.patterns.iter().map(|p| j.objective_of_pattern(p)).collect();
    j.set_weight(0, 1.0).unwrap();
    let w: Vec<f64> = j.patterns.iter().map(|p| p.weight).collect();
    assert!((w[0] - 2.0 / 3.0).abs() < 1e-15 && (w[1] - 1.0 / 3.0).abs() < 1e-15, "{w:?}");
    assert!((j.objective() - (w[0] * s[0] + w[1] * s[1])).abs() < 1e-9 * j.objective());
    j.set_weights(&[3.0, 1.0]).unwrap();
    assert_eq!(j.patterns.iter().map(|p| p.weight).collect::<Vec<_>>(), vec![0.75, 0.25]);
    assert!(j.set_weights(&[1.0, 0.0]).is_err() && j.set_weight(5, 1.0).is_err());
    j.set_refined(1, false).unwrap();
    assert!((j.objective() - 0.75 * s[0]).abs() < 1e-9 * j.objective());
    assert!(j.set_refined(0, false).is_err(), "the last refined pattern");
}

#[test]
fn contributions() {
    let mut j = scor();
    let before = j.patterns[1].model.calculate(&j.patterns[1].data.x);
    let shared_atoms = |j: &Joint| j.parameters.iter().filter(|q| q.name.contains("atom_") && q.members.len() == 2).count();
    let atoms = shared_atoms(&j);
    assert!(atoms > 0);
    let sc0 = match &j.patterns[1].model.phases[0].intensities {
        drs::powder::Intensities::Structure(_, sc) => sc.clone(),
        _ => unreachable!(),
    };
    // Out of the X-ray pattern.
    j.set_contribution(1, 0, false).unwrap();
    assert_eq!(j.patterns[1].phases, vec![1]);
    let after = j.patterns[1].model.calculate(&j.patterns[1].data.x);
    for i in 0..after.ycalc.len() {
        let want = before.ycalc[i] - before.per_phase[0][i];
        assert!((after.ycalc[i] - want).abs() <= 1e-9 * want.abs().max(1.0), "point {i}: {} against {want}", after.ycalc[i]);
    }
    assert_eq!(shared_atoms(&j), 0, "the structure's parameters pattern 1's alone");
    assert!(j.patterns[1].set.free.iter().all(|f| f.fields.iter().all(|(g, _)| g.phase() != Some(1))), "nothing of a phase 2 left");
    let n_out = j.parameters.len();
    // Given back.
    j.set_contribution(1, 0, true).unwrap();
    assert_eq!(j.patterns[1].phases, vec![1, 0]);
    let sc = match &j.patterns[1].model.phases[1].intensities {
        drs::powder::Intensities::Structure(_, sc) => sc.clone(),
        _ => unreachable!(),
    };
    assert_eq!(sc, sc0, "the atoms' X-ray scattering");
    assert_eq!(shared_atoms(&j), atoms, "its atoms tied to pattern 1's again");
    let own: Vec<&str> = j.parameters.iter().filter(|q| q.members.iter().all(|m| m.0 == 1) && q.name.contains("phase1.")).map(|q| q.name.as_str()).collect();
    assert!(!own.is_empty() && own.iter().all(|n| !n.contains("atom_")), "{own:?}");
    assert!(j.parameters.len() > n_out);
    assert!(j.set_contribution(1, 1, false).is_ok() && j.set_contribution(1, 0, false).is_err(), "a pattern's only phase stays");
    // A refinement runs on the joint refinement so changed.
    let r = j.refine(Options { max_cycles: 1, eps: 0.01 }, Esd::Default, |_| true);
    assert!(r.statistics.chi2.is_finite());
}
