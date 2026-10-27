//! A parameter across the versions of a control file (W7, L141): `trend::trend` and `drs trend` on two versions of the
//! Y2O3 control file, the second with a refined a and the results file of its refinement.

use drs::control::Control;
use std::path::{Path, PathBuf};

fn root(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)
}
fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap().iter().map(|&b| b as char).collect()
}

#[test]
fn a_parameter_across_versions() {
    let tmp = std::env::temp_dir().join(format!("drs_trend_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&tmp);
    std::fs::create_dir_all(tmp.join("Versions")).unwrap();
    let run = root("../manuals/fullprof/runs/Y2O3");
    std::fs::copy(run.join("y2o3.dat"), tmp.join("y2o3.dat")).unwrap();
    let p = drs::pcr::read(&latin1(&run.join("y2o3_structure.pcr"))).unwrap();
    let c = Control::from_pcr(&p, &run, "y2o3_structure", &tmp).unwrap();
    let a = drs::refinable::field_of("phase1.cell_a", &c.model).unwrap();
    // version 1: a moved off by 0.01 Å; version 2: the file as FullProf left it, refined by `drs results` (3 cycles)
    let a0 = c.model.get(a);
    let mut c1 = c.clone();
    c1.model.set(a, a0 + 0.01);
    let v1 = tmp.join("Versions/y2o3 2026-10-12 10.00.00.drc");
    std::fs::write(&v1, c1.to_text()).unwrap();
    let v2 = tmp.join("y2o3.drc");
    std::fs::write(&v2, c.to_text()).unwrap();
    let out = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(["results", v2.to_str().unwrap(), "--cycles", "3", "--time", "2026-10-12 11.00.00"]).output().unwrap();
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let r2 = tmp.join("Results/y2o3 2026-10-12 11.00.00 results.cif");
    assert!(r2.exists());
    let refined = drs::results::read(&std::fs::read_to_string(&r2).unwrap()).unwrap();
    let (_, _, su_a) = refined.refined.iter().find(|x| x.0 == "phase1.cell_a").cloned().unwrap();
    let t = drs::trend::trend(&[(v1.clone(), None), (v2.clone(), Some(r2.clone()))], Some(&tmp), "phase1.cell_a");
    let (p1, p2) = (t[0].as_ref().unwrap(), t[1].as_ref().unwrap());
    assert_eq!((p1.value, p1.su, p1.refined), (a0 + 0.01, None, false));
    assert_eq!((p2.value, p2.su, p2.refined), (a0, Some(su_a), true));
    assert!(su_a > 0.0 && su_a < 1e-3);
    // a value the model does not have is an error at its place, not a failure of the trend
    let bad = drs::trend::trend(&[(v1.clone(), None)], Some(&tmp), "phase2.cell_a");
    assert!(bad[0].as_ref().unwrap_err().contains("phase2.cell_a"));
    // the command line prints the same
    let out = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(["trend", "phase1.cell_a", v1.to_str().unwrap(), &format!("{}:{}", v2.display(), r2.display()), "--base", tmp.to_str().unwrap()]).output().unwrap();
    let text = String::from_utf8_lossy(&out.stdout);
    assert!(out.status.success() && text.contains(&format!("\t{}\t-", a0 + 0.01)) && text.contains(&format!("\t{a0}\t{su_a}")), "{text}");
    std::fs::remove_dir_all(&tmp).unwrap();
}
