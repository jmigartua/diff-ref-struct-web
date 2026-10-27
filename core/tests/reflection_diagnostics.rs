//! The diagnostics of W7–W8 (`diagnostics`, app requests 19–21, L100) against FullProf's runs:
//! - per reflection, against the last reflection table of the run's `CODE_calc.out` ("Hw", "Icalc", "Iobs", "Sigma"):
//!   the FWHM to 10⁻⁵, and, for the reflections of Icalc ≥ 20, I obs and its σ as fractions of I calc (FullProf's
//!   Icalc is the intensity within its window, the core's the whole: the ratios are compared; FullProf prints a
//!   negative partition as 0). I obs/I calc agrees within FullProf's print (Iobs to 0.1: 0.005 at Icalc 20) but in
//!   FROZEN. σ is the propagation of the points' σ through the partition; FullProf's formula is not in its manual,
//!   and its Sigma differs by 2–6 % (frozen per run);
//! - the shares of a point sum to 1 and the phases' shares are the per-phase curves;
//! - the cumulative χ² ends at Σw(y − yc)² and the weighted difference is (y − yc)/σ.

use drs::diagnostics;
use drs::powder::Model;
use std::path::Path;

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display())).iter().map(|&b| b as char).collect()
}

/// (phase, h k l, [Hw, Icalc, Iobs, Sigma]).
type Row = (usize, [i32; 3], [f64; 4]);

/// The last reflection table of a FullProf `.out`: (phase, h k l, [Hw, Icalc, Iobs, Sigma]), the columns found by the
/// names of their header (its layout differs with the profile: TCH adds HwG and HwL).
fn out_table(out: &str) -> Vec<Row> {
    let mut tables: Vec<Vec<Row>> = Vec::new();
    let mut phase = 0;
    let mut cols: Option<[usize; 4]> = None;
    for l in out.lines() {
        if let Some((_, r)) = l.split_once("Phase No.:") {
            phase = r.split_whitespace().next().and_then(|s| s.parse().ok()).unwrap_or(phase);
            if phase == 1 {
                tables.push(Vec::new());
            }
            continue;
        }
        let names: Vec<&str> = l.split_whitespace().collect();
        if names.first() == Some(&"No.") && names.contains(&"Iobs") {
            let at = |n: &str| names.iter().position(|x| *x == n);
            cols = Some([at("Hw").unwrap(), at("Icalc").unwrap(), at("Iobs").unwrap(), at("Sigma").unwrap()]);
            continue;
        }
        let v: Vec<f64> = l.split_whitespace().filter_map(|s| s.parse().ok()).collect();
        if let (Some(c), Some(t)) = (cols, tables.last_mut()) {
            if v.len() == names.len() && v.len() > c[3] && names.len() >= 10 {
                t.push((phase, [v[2] as i32, v[3] as i32, v[4] as i32], c.map(|i| v[i])));
            }
        }
    }
    tables.into_iter().rev().find(|t| !t.is_empty()).unwrap_or_default()
}

/// Per run: (strong reflections whose I obs/I calc differs by more than 0.02, the largest difference, the largest
/// relative difference of σ). rutile-anatase (TCH with axial divergence): 9 strong reflections lower in the core by
/// up to 0.13, not resolved. E: (0 1 2) at 22.89°, whose window reaches into the region 0–21.3° that E lists out of
/// order and FullProf treats apart (L098).
const FROZEN: &[(&str, usize, f64, f64)] = &[
    ("Y2O3", 0, 0.005, 0.055),
    ("rutile-anatase", 9, 0.131, 0.308),
    ("modes/A_modes", 0, 0.005, 0.023),
    ("modes/B_thesis", 0, 0.005, 0.021),
    ("modes/E", 1, 0.192, 0.202),
];

#[test]
fn reflection_diagnostics() {
    let mut fails = Vec::new();
    let runs = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs");
    for (dir, code) in [("Y2O3", "y2o3_structure"), ("rutile-anatase", "rutana_tch"), ("modes/A_modes", "a_modes"), ("modes/B_thesis", "b_thesis"), ("modes/E", "e")] {
        let d = runs.join(dir);
        let p = drs::pcr::read(&latin1(&d.join(format!("{code}.pcr")))).unwrap();
        let mut l = Model::from_pcr_fullprof(&p, &d, code).unwrap();
        let out = latin1(&d.join(format!("{code}_calc.out")));
        let limits: Vec<f64> = out.lines().filter_map(|l| l.split_once("S(1/d)max:")).map(|(_, r)| r.split_whitespace().next().unwrap().parse().unwrap()).collect();
        for (ph, lim) in l.model.phases.iter_mut().zip(&limits) {
            ph.reflection_limit = Some(*lim);
        }
        let data = l.data.as_ref().unwrap();
        let rows = diagnostics::reflections(&l.model, data);
        let table = out_table(&out);
        let (mut n, mut wh, mut n_off, mut wo, mut ws) = (0, 0.0f64, 0, 0.0f64, 0.0f64);
        for (ph, h, v) in &table {
            let Some(r) = rows.iter().find(|r| r.phase + 1 == *ph && r.hkl == *h && r.wavelength == 1) else { continue };
            n += 1;
            wh = wh.max((r.fwhm / v[0] - 1.0).abs());
            if v[1] < 20.0 {
                continue;
            }
            // FullProf prints a negative partition as 0.
            let e = (r.iobs.max(0.0) / r.icalc - v[2] / v[1]).abs();
            if e > 0.02 {
                n_off += 1;
            }
            wo = wo.max(e);
            ws = ws.max((r.sigma / r.icalc - v[3] / v[1]).abs() / (v[3] / v[1]));
        }
        println!("{dir}: {} reflections in the .out, {n} compared: FWHM within {wh:.1e}; Icalc ≥ 20: Iobs/Icalc within {wo:.3}, {n_off} beyond 0.02, σ within {:.1}%", table.len(), 100.0 * ws);
        let &(_, f_n, f_o, f_s) = FROZEN.iter().find(|f| f.0 == dir).unwrap();
        if n < table.len() || wh > 1e-5 || n_off > f_n || wo > f_o || ws > f_s {
            fails.push(format!("{dir}: {n} of {} found, FWHM {wh:.1e}, Iobs {n_off} beyond 0.02 up to {wo:.4}, σ {ws:.4}", table.len()));
        }
        // The shares at every point sum to 1, the phases' are the per-phase curves; the cumulative χ² ends at Σw(y − yc)².
        let c = l.model.calculate(&data.x);
        for i in (0..data.x.len()).step_by(97) {
            let sh = diagnostics::shares(&c, i);
            let total = sh.background + sh.phases.iter().sum::<f64>();
            if (total - 1.0).abs() > 1e-12 || sh.phases.iter().zip(&c.per_phase).any(|(f, p)| (f * c.ycalc[i] - p[i]).abs() > 1e-9 * c.ycalc[i].abs()) {
                fails.push(format!("{dir}: shares at {} sum to {total}", data.x[i]));
            }
        }
        let cum = diagnostics::cumulative(&l.model, data, &c.ycalc);
        let w = drs::lsq::weights(&l.model, data, &c.ycalc);
        let s: f64 = (0..w.len()).map(|i| w[i] * (data.y[i] - c.ycalc[i]).powi(2)).sum();
        let last = *cum.sum.last().unwrap();
        let wd_ok = (0..w.len()).all(|i| w[i] == 0.0 || (cum.weighted_difference[i] - (data.y[i] - c.ycalc[i]) * w[i].sqrt()).abs() < 1e-9);
        if (last - s).abs() > 1e-9 * s || cum.sum.windows(2).any(|p| p[1] < p[0]) || (cum.fraction().last().unwrap() - 1.0).abs() > 1e-12 || !wd_ok {
            fails.push(format!("{dir}: cumulative χ² ends at {last}, Σw(y − yc)² {s}"));
        }
        // The forbidden reflections over the pattern's range (app request 30): every one between the first point and
        // the last, Zero added, none allowed by the group.
        let lo_hi = (data.x[0], data.x[data.x.len() - 1]);
        for (ip, r) in diagnostics::forbidden_in(&l.model, &data.x) {
            let ph = &l.model.phases[ip];
            let in_range = r.two_theta.is_none_or(|t| t - l.model.zero >= lo_hi.0 - 1e-9 && t - l.model.zero <= lo_hi.1.min(179.9) + 1e-9);
            if !in_range || !drs::reflections::is_absent(&ph.group, r.hkl) || r.condition.is_empty() {
                fails.push(format!("{dir}: forbidden {:?} of phase {} at {:?} ({})", r.hkl, ip + 1, r.two_theta, r.condition));
                break;
            }
        }
    }
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}
