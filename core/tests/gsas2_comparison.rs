//! A refinement compared with GSAS-II's (the app's request 50, `compare::gsas2_prepare`, `gsas2_read`; test
//! `gsas2_comparison`, L399), on the reference run `manuals/gsas2/runs/y2o3_tch` (its README): a pattern the core
//! simulated from a Thompson–Cox–Hastings model of Y₂O₃, fitted from the same start by the core and by GSAS-II.
//!
//! - **Prepared**: the run's folder for `fit.drc`: the instrument file in GSAS-II's terms (U, V, W as σ² in centideg²,
//!   GSAS-II's X the core's Y and its Y the core's X, in centideg; P = 1/(1 + Cthm)), the CIF (U_iso = B/8π², the
//!   chemical occupancy), the flags of the input, the script, `run.json`.
//! - **Read and compared, on GSAS-II's stored output**: χ² within 0.005, every free parameter but the scale and the
//!   background compared, all within GSAS-II's esd but the cell, frozen at 1.01 of it.
//! - **Run** (only with `GSAS2_PYTHON` set to GSAS-II's Python and `GSAS2_PATH` to the folder holding `GSASII`, on
//!   the Mac): GSAS-II run on the prepared folder and read, with the same checks.

#![cfg(feature = "serde")]

use drs::compare::{self, Run};
use std::path::{Path, PathBuf};

fn run_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/gsas2/runs/y2o3_tch")
}

/// The checks of a comparison of `fit.drc` with GSAS-II.
fn check(c: &compare::Comparison, what: &str, fails: &mut Vec<String>) {
    let chi2 = &c.agreement[0];
    println!("{what}: χ² {:.4} against GSAS-II's {:?}; {} of {} parameters agree", chi2.core, chi2.other, c.parameters.iter().filter(|p| p.agrees).count(), c.parameters.len());
    for p in &c.parameters {
        println!("    {:32} {:12.6} ({:.1e}) GSAS-II {:12.6} ({:.1e}): {:+.2} σ", p.address, p.core_value, p.core_esd, p.other_value, p.other_esd, p.delta_sigma);
    }
    if chi2.difference.is_none_or(|d| d.abs() > 0.005) {
        fails.push(format!("{what}: χ² {:.4} against {:?}", chi2.core, chi2.other));
    }
    if c.parameters.len() != 12 {
        fails.push(format!("{what}: {} parameters compared, 12 expected; not compared {:?}", c.parameters.len(), c.not_compared));
    }
    // The cell: 4·10⁻⁵ Å apart, 1.00 of GSAS-II's esd, whose esd is 0.80 of the core's (the zero shift's is 1.17):
    // frozen at 1.01.
    for p in &c.parameters {
        let allowed = if p.address == "phase1.cell_a" { 1.01 } else { 1.0 };
        if p.delta_sigma.abs() > allowed {
            fails.push(format!("{what}: {} {:.6} against GSAS-II's {:.6}: {:+.2} σ", p.address, p.core_value, p.other_value, p.delta_sigma));
        }
    }
    let names: Vec<&str> = c.not_compared.iter().map(|n| n.name.as_str()).collect();
    if names != ["phase1.scale_factor", "the background"] {
        fails.push(format!("{what}: not compared {names:?}"));
    }
    if c.patterns.ycalc_other.iter().filter(|v| v.is_finite()).count() != c.patterns.x.len() {
        fails.push(format!("{what}: GSAS-II's Ycal not at every point"));
    }
}

#[test]
fn gsas2_comparison() {
    let mut fails: Vec<String> = Vec::new();
    let tmp = std::env::temp_dir().join(format!("drs_gsas2_comparison_{}", std::process::id()));
    let ours = drs::control::open(&run_dir().join("fit.drc")).unwrap();

    // Prepared.
    let r = compare::gsas2_prepare(&ours, &tmp, None).unwrap();
    let read = |f: &str| std::fs::read_to_string(tmp.join(f)).unwrap_or_else(|e| panic!("{f}: {e}"));
    let ins: Vec<(String, f64)> = read("compare.instprm").lines().filter_map(|l| l.split_once(':')).filter_map(|(k, v)| Some((k.to_string(), v.trim().parse().ok()?))).collect();
    let get = |k: &str| ins.iter().find(|x| x.0 == k).map(|x| x.1).unwrap_or(f64::NAN);
    let k = 1e4 / (8.0 * std::f64::consts::LN_2);
    for (key, want) in [("U", 0.02 * k), ("V", -0.002 * k), ("W", 0.008 * k), ("X", 100.0 * 0.03), ("Y", 100.0 * 0.02), ("Zero", 0.0), ("Polariz.", 1.0 / 1.78), ("Lam1", 1.5406), ("I(L2)/I(L1)", 0.5)] {
        if (get(key) - want).abs() > 1e-9 * want.abs().max(1.0) {
            fails.push(format!("instprm {key}: {} against {want}", get(key)));
        }
    }
    let cif = read("phase1.cif");
    let o = cif.lines().find(|l| l.starts_with("O ")).unwrap_or("").split_whitespace().collect::<Vec<_>>();
    let u: f64 = o.get(5).and_then(|v| v.parse().ok()).unwrap_or(f64::NAN);
    if (u - 0.5 / (8.0 * std::f64::consts::PI.powi(2))).abs() > 1e-12 || o.get(6).and_then(|v| v.parse::<f64>().ok()).is_none_or(|v| (v - 1.0).abs() > 1e-9) {
        fails.push(format!("the CIF's O: {o:?}"));
    }
    let y1 = cif.lines().find(|l| l.starts_with("Y ")).unwrap_or("").split_whitespace().nth(6).and_then(|v| v.parse::<f64>().ok()).unwrap_or(f64::NAN);
    if (y1 - 0.16666 * 6.0).abs() > 1e-9 {
        fails.push(format!("the CIF's first Y: chemical occupancy {y1}, FullProf's 0.16666 × 48/8"));
    }
    let input: serde_json::Value = serde_json::from_str(&read("compare_input.json")).unwrap();
    let inst: Vec<&str> = input["instrument_refine"].as_array().unwrap().iter().map(|v| v.as_str().unwrap()).collect();
    if inst != ["U", "V", "W", "X", "Y", "Zero"] || input["phases"][0]["refine_cell"] != true || input["phases"][0]["atoms"][2]["refine"] != "UX" || input["background"]["terms"] != 1 || input["cycles"] != 15 {
        fails.push(format!("the input's flags: {input}"));
    }
    if read("compare.py") != std::fs::read_to_string(run_dir().join("gsas2/compare.py")).unwrap() {
        fails.push("the script differs from the reference run's: run it again (manuals/gsas2/runs/README.md)".into());
    }
    if Run::load(&tmp).unwrap() != r || r.args != ["compare.py"] || r.env[0].0 != "PYTHONPATH" {
        fails.push(format!("the run: {r:?}"));
    }

    // Read and compared, on GSAS-II's stored output.
    let stored = run_dir().join("gsas2");
    let run = Run::load(&stored).unwrap();
    let c = compare::gsas2_read(&stored, &run, &ours).unwrap();
    if c.other.version != "GSAS-II 5843 (Python 3.13.16)" || !c.other.converged {
        fails.push(format!("GSAS-II's version and convergence: {} {}", c.other.version, c.other.converged));
    }
    check(&c, "stored", &mut fails);

    // Run, where GSAS-II is installed.
    match (std::env::var("GSAS2_PYTHON"), std::env::var("GSAS2_PATH")) {
        (Ok(py), Ok(path)) => {
            let out = std::process::Command::new(&py).args(&r.args).current_dir(&r.dir).env("PYTHONPATH", &path).stdin(std::process::Stdio::null()).output().unwrap();
            if !out.status.success() {
                fails.push(format!("GSAS-II: {}\n{}", out.status, String::from_utf8_lossy(&out.stderr).lines().rev().take(5).collect::<Vec<_>>().join("\n")));
            } else {
                let c = compare::gsas2_read(&r.dir, &r, &ours).unwrap();
                check(&c, "run", &mut fails);
            }
        }
        _ => println!("run: GSAS2_PYTHON and GSAS2_PATH not set, GSAS-II not run"),
    }
    let _ = std::fs::remove_dir_all(&tmp);
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}
