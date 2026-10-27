//! The results file (R28, pdCIF chosen in L092; `results`, L101), on case A of the mode course refined one cycle from
//! FullProf's values (`manuals/fullprof/runs/modes/A_modes/a_modes.pcr`, as `drs results … --cycles 1` runs it):
//! - written and read back: the profile (2θ, Yobs and its σ, the weights, the background, Ycal, the phase's curve), the
//!   reflections (`diagnostics::reflections`), the agreement factors, the refined parameters with their uncertainties
//!   and the control file (read again, it writes the same text) are the core's, to the seven figures written;
//! - the structure CIF of the phase, read by the core's CIF reader, gives the refined cell and atoms;
//! - the committed file `tests/data/results/a_modes_results.cif` (and its phase's CIF) is what the core writes today,
//!   to the figures written; GSAS-II's reading of it (`tests/data/gsas2_results.txt`, written by
//!   `make_gsas2_results.py` with GSAS-II's own powder and phase importers, commit 305a159) is the core's: the points,
//!   Σx, Σy, Σw with GSAS-II's weights (1/max(σ, 1)², σ from the su of Yobs, else 1/y; it does not read
//!   `_pd_proc_ls_weight`), the first and last points, the wavelength, and each phase's cell, space group and atoms.

use drs::control::Control;
use drs::results;
use std::path::Path;

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display())).iter().map(|&b| b as char).collect()
}

fn close(a: f64, b: f64, rel: f64) -> bool {
    (a - b).abs() <= rel * a.abs().max(b.abs()).max(1e-12) || (a.is_nan() && b.is_nan())
}

#[test]
fn results_file() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let dir = root.join("../manuals/fullprof/runs/modes/A_modes");
    let p = drs::pcr::read(&latin1(&dir.join("a_modes.pcr"))).unwrap();
    let mut control = Control::from_pcr(&p, &dir, "a_modes", &dir).unwrap();
    let data = drs::powder::Model::from_pcr(&p, &dir, "a_modes").unwrap().data.unwrap();
    let set = control.refinement.build(&control.model).unwrap();
    let r = drs::lsq::refine(&mut control.model, &data, &set, drs::lsq::Options { max_cycles: 1, eps: 0.1 }, |_| {});
    let (text, structures) = results::write(&control, &data, Some(&set), Some(&r), "a_modes 2026-10-11 12.00.00").unwrap();
    let mut fails = Vec::new();
    let back = results::read(&text).unwrap();
    // The profile.
    let c = control.model.calculate(&data.x);
    let w = drs::lsq::weights(&control.model, &data, &c.ycalc);
    let n = data.x.len();
    let col_ok = |name: &str, got: &[f64], want: &[f64], rel: f64| -> Option<String> {
        if got.len() != want.len() {
            return Some(format!("{name}: {} values, {} expected", got.len(), want.len()));
        }
        (0..got.len()).find(|&i| !close(got[i], want[i], rel)).map(|i| format!("{name}[{i}]: {} against {}", got[i], want[i]))
    };
    fails.extend(col_ok("2θ", &back.x, &data.x, 1e-12));
    fails.extend(col_ok("Yobs", &back.yobs, &data.y, 1e-12));
    fails.extend(col_ok("σ", &back.sigma, &data.sigma, 1e-3));
    fails.extend(col_ok("weight", &back.weight, &w, 1e-6));
    fails.extend(col_ok("background", &back.background, &c.background, 1e-6));
    fails.extend(col_ok("Ycal", &back.ycalc, &c.ycalc, 1e-6));
    fails.extend(col_ok("phase 1", &back.per_phase[0], &c.per_phase[0], 1e-6));
    // The reflections.
    let rows = drs::diagnostics::reflections(&control.model, &data);
    if back.reflections.len() != rows.len() {
        fails.push(format!("{} reflections read, {} written", back.reflections.len(), rows.len()));
    }
    for (b, r) in back.reflections.iter().zip(&rows) {
        let want = [r.d, r.two_theta, r.fwhm, r.f2.unwrap_or(f64::NAN), r.icalc, r.iobs, r.sigma];
        if b.0 != r.phase + 1 || b.1 != r.wavelength || b.2 != r.hkl || b.3 != r.multiplicity || (0..7).any(|j| !close(b.4[j], want[j], 1e-6)) {
            fails.push(format!("reflection {:?}: {b:?}", r.hkl));
            break;
        }
    }
    // The agreement factors, the parameters, the control file.
    let s = &r.stats;
    if !(close(back.rp, s.rp, 1e-6) && close(back.rwp, s.rwp, 1e-6) && close(back.rexp, s.rexp, 1e-6) && close(back.chi2, s.chi2, 1e-6) && back.points_used == s.n && back.parameters == s.p) {
        fails.push(format!("agreement: {} {} {} {} {} {} against {s:?}", back.rp, back.rwp, back.rexp, back.chi2, back.points_used, back.parameters));
    }
    if back.refined.len() != r.fields.len() || back.refined.iter().zip(&r.fields).any(|(a, b)| a.0 != b.0 || !close(a.1, b.1, 1e-6) || !close(a.2, b.2, 1e-6)) {
        fails.push(format!("{} parameters read, {} written", back.refined.len(), r.fields.len()));
    }
    match Control::from_text(&back.control_file) {
        Ok(c2) if c2.to_text() == control.to_text() => {}
        Ok(_) => fails.push("the control file read back writes another text".into()),
        Err(e) => fails.push(format!("the control file read back: {e:?}")),
    }
    // The structure CIF of the phase.
    assert_eq!(structures.len(), 1);
    let st = drs::cif::structure(&structures[0].1).unwrap();
    // Written to their standard uncertainty (the IUCr's rule): within half of it.
    let cell = control.model.phases[0].cell;
    let su = |name: &str| r.fields.iter().find(|f| f.0 == name).map_or(0.0, |f| f.2);
    let within = |a: f64, b: f64, name: &str| (a - b).abs() <= 0.5 * su(name) + 1e-9;
    if !(within(st.cell.a, cell.a, "phase1.cell_a") && within(st.cell.b, cell.b, "phase1.cell_b") && within(st.cell.c, cell.c, "phase1.cell_c") && within(st.cell.beta, cell.beta, "phase1.cell_beta")) || st.atoms.len() != 6 {
        fails.push(format!("structure CIF: {:?}, {} atoms", st.cell, st.atoms.len()));
    }
    // The committed file is what the core writes, to the figures written.
    let committed = results::read(&latin1(&root.join("tests/data/results/a_modes_results.cif"))).unwrap();
    fails.extend(col_ok("committed Ycal", &committed.ycalc, &back.ycalc, 1e-6));
    fails.extend(col_ok("committed weight", &committed.weight, &back.weight, 1e-6));
    if committed.reflections.len() != back.reflections.len() || !close(committed.rwp, back.rwp, 1e-6) {
        fails.push("the committed results file is not what the core writes: regenerate it and gsas2_results.txt".into());
    }
    // GSAS-II's reading of the committed files.
    let g2 = latin1(&root.join("tests/data/gsas2_results.txt"));
    let fsum = |v: &[f64]| v.iter().sum::<f64>();
    let g2w: Vec<f64> = committed.sigma.iter().zip(&committed.yobs).map(|(&s, &y)| if s > 0.0 { 1.0 / s.max(1.0).powi(2) } else { 1.0 / y }).collect();
    for l in g2.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = l.split_whitespace().collect();
        let v = |i: usize| f[i].parse::<f64>().unwrap();
        match f[0] {
            "profile" => {
                let last = committed.x.len() - 1;
                let ours = [committed.x.len() as f64, fsum(&committed.x), fsum(&committed.yobs), fsum(&g2w), committed.x[0], committed.yobs[0], g2w[0], committed.x[last], committed.yobs[last], g2w[last], control.model.lambda[0]];
                if (1..12).any(|i| !close(v(i), ours[i - 1], 1e-9)) {
                    fails.push(format!("GSAS-II reads the profile as {l}; the core {ours:?}"));
                }
            }
            "phase" => {
                let ours = drs::cif::structure(&latin1(&root.join("tests/data/results").join(f[1]))).unwrap();
                let oc = [ours.cell.a, ours.cell.b, ours.cell.c, ours.cell.alpha, ours.cell.beta, ours.cell.gamma];
                let sum: f64 = ours.atoms.iter().map(|a| a.x[0] + a.x[1] + a.x[2]).sum();
                let hm = ours.group.setting.map(|s| s.hm.replace([' ', '_'], "")).unwrap_or_default();
                if (0..6).any(|k| !close(v(2 + k), oc[k], 1e-9)) || !hm.starts_with(f[8].trim_end_matches('1')) && hm != f[8] || f[9].parse::<usize>().unwrap() != ours.atoms.len() || !close(v(10), sum, 1e-9) {
                    fails.push(format!("GSAS-II reads {} as {l}; the core: {oc:?} {hm} {} atoms, Σ {sum}", f[1], ours.atoms.len()));
                }
            }
            _ => fails.push(format!("unknown line {l}")),
        }
    }
    println!("{n} points, {} reflections, {} parameters; Rwp {:.3} %; read back, and as GSAS-II reads it", back.reflections.len(), back.refined.len(), back.rwp);
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}
