//! A model whose pattern is not a number is refused by `drs calc` and `drs simulate` (the tutorial's TU6): Y2O3's
//! control file with the March–Dollase G₁ = 0, below its limit 10⁻³, where the function is undefined. Before core
//! 0.27.0 `drs simulate` never ended on it (the Poisson sampler's rejection loop never accepts a NaN mean).

use std::path::{Path, PathBuf};

fn root(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)
}

#[test]
fn a_pattern_that_is_not_a_number_is_refused() {
    let d = root("target/test-simulate-nan");
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    let run = root("../manuals/fullprof/runs/Y2O3");
    for f in ["y2o3.dat", "y2o3_structure.pcr"] {
        std::fs::copy(run.join(f), d.join(f)).unwrap();
    }
    let pcr = drs::pcr::read(&std::fs::read(d.join("y2o3_structure.pcr")).unwrap().iter().map(|&b| b as char).collect::<String>()).unwrap();
    let mut c = drs::control::Control::from_pcr(&pcr, &d, "y2o3_structure", &d).unwrap();
    c.model.phases[0].pref[0] = 0.0;
    // No peak reaches 20°; the first NaN is at a reflection's point.
    let calc = c.model.calculate(&[20.0, 29.15, 33.8]);
    let e = c.model.finite(&calc).unwrap_err();
    assert!(e.contains("phase1.po_g1 = 0") && e.contains("point 3 (2θ = 33.8)"), "{e}");
    let drc = d.join("g0.drc");
    std::fs::write(&drc, c.to_text()).unwrap();
    for cmd in ["calc", "simulate"] {
        let t = std::time::Instant::now();
        let o = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args([cmd, drc.to_str().unwrap()]).output().unwrap();
        let err = String::from_utf8_lossy(&o.stderr);
        assert!(!o.status.success() && err.contains("not a number") && err.contains("phase1.po_g1"), "{cmd}: {err}");
        assert!(t.elapsed().as_secs() < 60, "{cmd} took {:?}", t.elapsed());
    }
    // The same file with G₁ = 1 is calculated and simulated.
    c.model.phases[0].pref[0] = 1.0;
    assert!(c.model.finite(&c.model.calculate(&[20.0, 29.15, 33.8])).is_ok());
    std::fs::write(&drc, c.to_text()).unwrap();
    let o = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(["simulate", drc.to_str().unwrap(), "--out", d.join("sim.xy").to_str().unwrap()]).output().unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
}
