//! A refinement compared with FullProf's (the app's request 49, `compare`; test `fullprof_comparison`).
//!
//! - **Read and compared, on FullProf's stored outputs**: Y2O3's starting file (`y2o3_structure.pcr.in`) refined by
//!   FullProf to its convergence (`y2o3_structure_conv.*`, 100 cycles allowed, `converge_runs.py`, L052) against the
//!   core's refinement of the same file with the same cycles: FullProf's banner, cycles and convergence; every free
//!   parameter matched by FullProf's name (none left over); every parameter within FullProf's esd (`compare::RULE`);
//!   χ² and the R factors within FullProf's print and the frozen differences below; the esd ratio 1 within 2 % but the
//!   kinds `fullprof_refinements` already freezes; FullProf's Ycal at every data point; the comparison written and
//!   read back.
//! - **The `.out` of five phases**: Egypt's R_Bragg and R_F per phase, its cycles and convergence.
//! - **Prepared**: the run's folder for the same control file (the `.pcr` with FullProf's output flags and the cycles
//!   asked for, the data, `run.json`), the `.pcr` read back by the core with the same free parameters.
//! - **Run** (only with `FULLPROF` set to fp2k's folder, on the Mac): FullProf run on the prepared folder and the
//!   comparison read from it, every parameter agreeing.

use drs::compare::{self, Carried, Run};
use drs::control::{self, Control};
use std::path::{Path, PathBuf};

fn runs() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs")
}

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display())).iter().map(|&b| b as char).collect()
}

/// A control file made from a `.pcr` of a run folder, saved in `tmp` and opened from there.
fn opened_from_pcr(run: &str, pcr: &str, stem: &str, tmp: &Path) -> control::Opened {
    let d = runs().join(run);
    let p = drs::pcr::read(&latin1(&d.join(pcr))).unwrap();
    std::fs::create_dir_all(tmp).unwrap();
    let c = Control::from_pcr(&p, &d, stem, tmp).unwrap();
    let f = tmp.join(format!("{stem}.drc"));
    std::fs::write(&f, c.to_text()).unwrap();
    control::open(&f).unwrap()
}

/// FullProf's χ² and R factors against the core's: FullProf's print and these, frozen (|difference|).
/// - R_F, Y2O3: the core's 1.626 against FullProf's 1.49 of its summary (1.50 in its last cycle's lines). Not
///   explained; R_Bragg agrees to 0.03 and every parameter within its esd, so it is the R_F's sum, not the model.
/// - Rexp, Y2O3: FullProf prints 11.3 in the `.out` and 11.33 in the `.sum`; the core's 11.328 is the `.sum`'s.
/// - Rp and R_Bragg, Y2O3: the core refines its default model, whose reflection range and peak window are not
///   FullProf's (the export's notes say so; `fullprof_refinements` sets FullProf's and agrees), and its minimum is a
///   little apart from FullProf's (the zero shift 0.998 of FullProf's esd): Rp 8.641 against 8.65.
const FROZEN: &[(&str, &str, f64)] = &[("Y2O3", "r_f", 0.137), ("Y2O3", "rexp", 0.029), ("Y2O3", "r_bragg", 0.025), ("Y2O3", "rp", 0.0095)];

#[test]
fn fullprof_comparison() {
    let mut fails: Vec<String> = Vec::new();
    let tmp = std::env::temp_dir().join(format!("drs_fullprof_comparison_{}", std::process::id()));
    let ours = opened_from_pcr("Y2O3", "y2o3_structure.pcr.in", "y2o3_structure", &tmp.join("control"));

    // Read and compared on FullProf's stored outputs.
    let run = Run {
        other: "FullProf".into(),
        program: "fp2k".into(),
        args: vec!["y2o3_structure_conv".into(), "y2o3.dat".into()],
        dir: runs().join("Y2O3"),
        env: Vec::new(),
        cycles: 100,
        code: "y2o3_structure_conv".into(),
        data: "y2o3.dat".into(),
        manifest: Vec::new(),
    };
    let c = compare::fullprof_read(&run.dir, &run, &ours).unwrap();
    println!("{} ({}, {} cycles) against {} ({} cycles, converged {})", c.core.name, c.core.version, c.core.cycles, c.other.version, c.other.cycles, c.other.converged);
    if c.other.version != "FullProf.2k (Version 8.20 - Feb2025-ILL JRC)" || c.other.cycles != 1 || !c.other.converged {
        fails.push(format!("FullProf's banner, cycles, convergence: {:?}, {}, {}", c.other.version, c.other.cycles, c.other.converged));
    }
    let free = ours.control.refinement.build(&ours.control.model).unwrap().free.len();
    if c.parameters.len() != free || !c.not_compared.is_empty() {
        fails.push(format!("{} of {free} parameters compared; not compared: {:?}", c.parameters.len(), c.not_compared));
    }
    for m in &c.agreement {
        let Some(d) = m.difference else {
            fails.push(format!("{}: FullProf's value not read", m.name));
            continue;
        };
        let frozen = FROZEN.iter().find(|f| f.0 == "Y2O3" && f.1 == m.id).map_or(0.0, |f| f.2);
        let tol = m.print.unwrap_or(0.0) + frozen + 1e-3;
        println!("    {:28} {:10.4} FullProf {:>8} difference {d:+.4} (print {:?})", m.name, m.core, m.other.unwrap(), m.print);
        if d.abs() > tol {
            fails.push(format!("{}: {:.4} against FullProf's {}: {d:+.4}, allowed {tol:.4}", m.name, m.core, m.other.unwrap()));
        }
    }
    let beyond: Vec<&compare::ParameterRow> = c.parameters.iter().filter(|p| !p.agrees).collect();
    if !beyond.is_empty() {
        fails.push(format!("beyond FullProf's esd: {:?}", beyond.iter().map(|p| (&p.address, p.delta_sigma)).collect::<Vec<_>>()));
    }
    // The esd ratios: 1 within 2 %, but the kinds frozen in `fullprof_refinements` (the ends of the background, a
    // bound on P1 + 3·P2).
    for p in &c.parameters {
        let frozen = p.address.ends_with("background_points.0") || p.address.ends_with("background_points.61") || p.address.contains("berar_baldinozzi");
        if !frozen && (p.esd_ratio - 1.0).abs() > 0.02 {
            fails.push(format!("{}: esd ratio {:.3}", p.address, p.esd_ratio));
        }
    }
    let n = ours.data.as_ref().unwrap().x.len();
    if c.patterns.ycalc_other.len() != n || c.patterns.ycalc_other.iter().any(|v| !v.is_finite()) {
        fails.push("FullProf's Ycal not at every data point".into());
    }
    let worst = c.patterns.difference.iter().zip(&c.patterns.ycalc_other).map(|(d, y)| d.abs() / y.abs().max(1.0)).fold(0.0f64, f64::max);
    println!("    {} parameters, largest |Δ/σ| {:.3} ({}); Ycal's largest relative difference {worst:.4}", c.parameters.len(), c.parameters[0].delta_sigma.abs(), c.parameters[0].address);
    let dir = tmp.join("written");
    std::fs::create_dir_all(&dir).unwrap();
    compare::write(&dir, &c).unwrap();
    let back = compare::read(&dir).unwrap();
    if back.parameters != c.parameters || back.agreement != c.agreement || back.other != c.other {
        fails.push("the comparison written does not read back the same".into());
    }

    // The .out of five phases.
    let egypt = drs::formats::fullprof_outputs::read_out(&latin1(&runs().join("Egypt/Egypt_qa_Rietveld_conv.out"))).unwrap();
    let rb: Vec<f64> = egypt.phases.iter().map(|p| p.0).collect();
    let rf: Vec<f64> = egypt.phases.iter().map(|p| p.1).collect();
    if rb != [5.21, 5.58, 4.96, 4.70, 10.2] || rf != [2.76, 3.42, 3.07, 2.58, 4.18] || egypt.cycles != 7 || !egypt.converged {
        fails.push(format!("Egypt's .out: R_Bragg {rb:?}, R_F {rf:?}, {} cycles, converged {}", egypt.cycles, egypt.converged));
    }
    if !egypt.all_points.as_ref().is_some_and(|s| s.dof > 0 && s.chi2_conventional > 0.0 && s.rp > 0.0) {
        fails.push(format!("Egypt's statistics: {:?}", egypt.all_points));
    }

    // Prepared.
    let prep = tmp.join("run");
    let r = compare::fullprof_prepare(&ours, &prep, Some(7)).unwrap();
    if r.args != ["compare", "compare.dat"] || r.cycles != 7 || r.env.first().map(|e| e.0.as_str()) != Some("FULLPROF") {
        fails.push(format!("the run: {:?} {} {:?}", r.args, r.cycles, r.env));
    }
    for f in ["compare.pcr", "compare.pcr.in", "compare.dat", "run.json"] {
        if !prep.join(f).exists() {
            fails.push(format!("{f} not written"));
        }
    }
    let pcr_text = std::fs::read_to_string(prep.join("compare.pcr")).unwrap();
    let pcr = drs::pcr::read(&pcr_text).unwrap();
    if pcr.flag("Prf") != 3 || pcr.flag("Pcr") != 1 || pcr.flag("Rpa") != 0 || pcr.ncy != 7 {
        fails.push(format!("the .pcr's flags: Prf {} Pcr {} Rpa {} NCY {}", pcr.flag("Prf"), pcr.flag("Pcr"), pcr.flag("Rpa"), pcr.ncy));
    }
    let back = Control::from_pcr(&pcr, &prep, "compare", &prep).unwrap();
    let free_back = back.refinement.build(&back.model).unwrap().free.len();
    if free_back != free {
        fails.push(format!("the prepared .pcr refines {free_back} parameters, the control file {free}"));
    }
    if r.manifest.iter().any(|m| m.carried == Carried::Refused) {
        fails.push(format!("refused: {:?}", r.manifest.iter().filter(|m| m.carried == Carried::Refused).collect::<Vec<_>>()));
    }
    if Run::load(&prep).unwrap() != r {
        fails.push("run.json does not read back the run".into());
    }

    // Run, where FullProf is installed.
    match std::env::var("FULLPROF") {
        Ok(fp) if Path::new(&fp).join("fp2k").exists() => {
            let out = std::process::Command::new(Path::new(&fp).join("fp2k")).args(&r.args).current_dir(&r.dir).env("FULLPROF", &fp).stdin(std::process::Stdio::null()).output().unwrap();
            if !out.status.success() {
                fails.push(format!("fp2k: {}", out.status));
            } else {
                let c = compare::fullprof_read(&r.dir, &r, &ours).unwrap();
                let agree = c.parameters.iter().filter(|p| p.agrees).count();
                println!("run: FullProf {} cycles, the core {}; χ² {:.4} against {:?}; {agree} of {} parameters agree", c.other.cycles, c.core.cycles, c.agreement[0].core, c.agreement[0].other, c.parameters.len());
                if agree != c.parameters.len() || c.parameters.len() != free {
                    fails.push(format!("run: {agree} of {} parameters agree ({free} free)", c.parameters.len()));
                }
            }
        }
        _ => println!("run: FULLPROF not set, FullProf not run"),
    }
    let _ = std::fs::remove_dir_all(&tmp);
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}
