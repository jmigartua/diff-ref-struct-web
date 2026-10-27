//! Size and strain against FullProf (test `fullprof_microstructure`).
//!
//! Two kinds of reference, both from the FullProf runs of `manuals/fullprof/runs/` (L034, L038):
//!
//! 1. The widths of every reflection in the `.out` of each run with isotropic broadening (Hw, HwG, HwL, ETA), against
//!    the core's widths of the phase read from `CODE.pcr` (the parameters FullProf read, L035). Tolerance: half a unit
//!    of the printed decimal and 2·10⁻⁵ of the value (FullProf works in single precision and prints 2θ to three
//!    decimals; the worst found is 1.9·10⁻⁵, YVO3).
//!
//! 2. The microstructure files (`CODE_calc1.mic`) of the three runs with a resolution file and isotropic size or strain
//!    (YVO3: strain X, size Y; Yb2ZrWO8: size Y and I_G; NdNi-HTO, nuclear phase: strain U and X), row by row:
//!    - HGo, HGi, HLo, HLi, betaG, betaL, Beta_TCH, HpVo and eta, against the phase of `CODE.pcr`, to half a unit of
//!      the fourth decimal printed and 2·10⁻⁵ of the value;
//!    - App-size and Max-strain, against the phase of `CODE_calc.pcr`: FullProf computes them after its refinement
//!      cycle, from the shifted values (YVO3, whose cycle moves U from 0 to −0.0914 and X from 0.645 to 0.821, prints
//!      a strain of 60.8810 that the values read would give as 44.2415, and the shifted ones as 60.8810). Tolerance:
//!      the effect of half a unit of the sixth decimal of U, X, Y, I_G as `CODE_calc.pcr` prints them, plus half a
//!      unit of the printed value. An apparent size above 99999 Å is printed as 99999.
//!    - where the resolution file gives a negative instrumental H_G² (NdNi-HTO above 150°), FullProf prints HGi and
//!      betaG as NaN and Beta_TCH as 0; the core must do the same;
//!    - the column "beta" (the de Keijser formula) and the observed intensities are not compared.
//!
//! The magnetic phase of NdNi-HTO (Jbt = ±1) is not calculated by the core yet and is left out of its model.
//!
//! 3. `williamson_hall_axes` (Explore's E19, L376): the Williamson–Hall plot's axes the core gives with each
//!    reflection's breadths, d* = 1/d against the `.mic`'s s and β* in Å⁻¹ against its betaG and betaL (×1000), the
//!    reflection's hkl its own, on NdNi-HTO's first phase (FullProf's NdNi_HTO_calc1.mic).

use drs::powder::Model;
use std::path::Path;

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display())).iter().map(|&b| b as char).collect()
}

fn dir(run: &str) -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs").join(run)
}

/// The model of a run's `.pcr`, its magnetic phases left out.
fn model(run: &str, pcr: &str) -> Model {
    let d = dir(run);
    let mut p = drs::pcr::read(&latin1(&d.join(format!("{pcr}.pcr")))).unwrap();
    p.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
    // The .pcr's values under FullProf's own rules (λ1 in β*, √|U|; L046, L050).
    Model::from_pcr_fullprof(&p, &d, pcr).unwrap_or_else(|e| panic!("{run}: {e}")).model
}

/// Half a unit of the last decimal of a printed number.
fn half(s: &str) -> f64 {
    0.5 * 10f64.powi(-(s.split_once('.').map_or(0, |(_, d)| d.len()) as i32))
}

/// Differences not resolved, frozen at what was measured (L039; CLAUDE.md, the gate): the run and the largest ratio
/// to the tolerance allowed. Y2O3 (pseudo-Voigt, Npr = 5): FullProf's η goes as η0 + X·2θ with a slope 7·10⁻⁵ below
/// the X = 0.003508 it read (from its own rows: 0.605175 at 20.497°, 1.040431 at 144.581°, slope 0.0035077), up to
/// 4·10⁻⁵ in η at 145°; the cause is not known. The profile is not affected beyond its print (fullprof_profiles).
const FROZEN: [(&str, f64); 1] = [("Y2O3", 1.89)];

/// Isotropic runs and their `.pcr` code.
const WIDTH_RUNS: [(&str, &str); 9] = [
    ("Y2O3", "y2o3_structure"),
    ("Co-YBCO", "Cocus-asym"),
    ("MgAl2O4-s", "MgAl2O4s_Rietveld"),
    ("rutile-anatase", "rutana_tch"),
    ("Egypt", "Egypt_qa_Rietveld"),
    ("YVO3", "YVO3a"),
    ("Yb2ZrWO8", "ybfluo_Riet_IRF"),
    ("NdNi-HTO", "NdNi_HTO"),
    ("LaB6", "LaB6rx"),
];

#[test]
fn fullprof_microstructure() {
    let mut fails = Vec::new();
    // 1. Widths of the .out reflection lists.
    let mut rows_total = 0;
    for (run, code) in WIDTH_RUNS {
        let m = model(run, code);
        let out = latin1(&dir(run).join(format!("{code}_calc.out")));
        let (mut ph, mut n, mut worst, mut at) = (usize::MAX, 0, 0.0f64, String::new());
        let mut hdr: Vec<String> = Vec::new();
        for (no, l) in out.lines().enumerate() {
            if l.contains("Phase No.:") && l.contains("Pattern#") {
                ph = l.split("Phase No.:").nth(1).unwrap().split_whitespace().next().unwrap().parse::<usize>().unwrap() - 1;
            }
            if l.trim_start().starts_with("No.  Code") {
                hdr = l.split_whitespace().map(String::from).collect();
                continue;
            }
            let f: Vec<&str> = l.split_whitespace().collect();
            if hdr.is_empty() || f.len() + 1 < hdr.len() || f[0].parse::<usize>().is_err() || ph >= m.phases.len() {
                continue;
            }
            let col = |name: &str| hdr.iter().position(|h| h == name).and_then(|i| f.get(i).copied());
            let (Some(tt), Some(hw)) = (col("2theta").or(col("2theta/TOF")), col("Hw")) else { continue };
            // At the Bragg angle of FullProf's d (six decimals; its 2θ has three, and its cell may differ from the
            // .pcr's in the sixth digit, L028).
            // Only the final reflection lists, which print d; Code 2 marks the λ2 reflections.
            let Some(d) = col("d-hkl").and_then(|s| s.parse::<f64>().ok()) else { continue };
            let lambda = if f[1] == "2" { m.lambda[1] } else { m.lambda[0] };
            let _ = tt;
            let tt = 2.0 * (lambda / (2.0 * d)).asin().to_degrees();
            let (h, eta) = m.width(&m.phases[ph], tt);
            let hg = m.gauss_width(&m.phases[ph], tt);
            let mut check = |name: &str, s: Option<&str>, v: f64| {
                if let Some(s) = s {
                    let r: f64 = s.parse().unwrap();
                    let e = (v - r).abs() / (half(s) + 2e-5 * r.abs());
                    if e > worst {
                        (worst, at) = (e, format!("{name} {v:.6} against {s} (line {})", no + 1));
                    }
                }
            };
            check("Hw", Some(hw), h);
            if m.phases[ph].shape == drs::powder::Shape::ThompsonCoxHastings {
                check("HwG", col("HwG").or(col("HWG")), hg);
                check("HwL", col("HwL").or(col("HWL")), m.width_components(&m.phases[ph], tt).1);
            }
            check("ETA", col("ETA").or(col("ETA/M")), eta);
            n += 1;
        }
        println!("{run}: {n} reflections of the .out, widths at {worst:.2} of the tolerance ({at})");
        let allowed = FROZEN.iter().find(|(r, _)| *r == run).map_or(1.0, |&(_, a)| a);
        if n == 0 || worst > allowed {
            fails.push(format!("{run}: widths at {worst:.2} of the tolerance: {at}"));
        }
        rows_total += n;
    }

    // 2. The microstructure files.
    for (run, code) in [("YVO3", "YVO3a"), ("Yb2ZrWO8", "ybfluo_Riet_IRF"), ("NdNi-HTO", "NdNi_HTO")] {
        let read = model(run, code);
        let refined = model(run, &format!("{code}_calc"));
        let mic = latin1(&dir(run).join(format!("{code}_calc1.mic")));
        let (mut n, mut worst, mut at, mut worst_ss, mut at_ss) = (0, 0.0f64, String::new(), 0.0f64, String::new());
        for (no, l) in mic.lines().enumerate() {
            let f: Vec<&str> = l.split_whitespace().collect();
            if f.len() < 20 || l.starts_with('!') || f[0].parse::<f64>().is_err() {
                continue;
            }
            let hkl = [f[11].parse::<i32>().unwrap(), f[12].parse().unwrap(), f[13].parse().unwrap()];
            let tt: f64 = f[14].parse().unwrap();
            // The wavelength whose Bragg angle is the printed 2θ.
            let (wl, b) = (1..=2)
                .filter_map(|w| read.breadths(0, hkl, w).map(|b| (w, b)))
                .min_by(|a, b| (a.1.two_theta - tt).abs().partial_cmp(&(b.1.two_theta - tt).abs()).unwrap())
                .unwrap_or_else(|| panic!("{run}: {hkl:?} at {tt}"));
            if (b.two_theta - tt).abs() > 2e-3 {
                fails.push(format!("{run} line {}: {hkl:?} at {tt}, core {:.4}", no + 1, b.two_theta));
                continue;
            }
            let cols = [(1, b.hg_obs, "HGo"), (2, b.hg_ins, "HGi"), (3, b.hl_obs, "HLo"), (4, b.hl_ins, "HLi"), (5, b.beta_g, "betaG"), (6, b.beta_l, "betaL"), (15, b.h_obs, "HpVo"), (16, b.eta_obs, "eta"), (19, b.beta_tch, "Beta_TCH")];
            for (i, v, name) in cols {
                // NaN where FullProf's instrumental H_G² is negative: the core must say the same.
                let e = if f[i] == "NaN" {
                    if v.is_nan() { 0.0 } else { f64::INFINITY }
                } else {
                    let r: f64 = f[i].parse().unwrap();
                    (v - r).abs() / (half(f[i]) + 2e-5 * r.abs())
                };
                if e > worst {
                    (worst, at) = (e, format!("{name} of {hkl:?} {v:.5} against {} (line {})", f[i], no + 1));
                }
            }
            // Apparent size and strain from the refined parameters, with the print of those parameters.
            let after = refined.breadths(0, hkl, wl).unwrap();
            let spread = |get: &dyn Fn(&drs::microstructure::Breadths) -> Option<f64>| -> f64 {
                let mut s = 0.0;
                for k in 0..4 {
                    let mut m2 = refined.clone();
                    let p = &mut m2.phases[0];
                    match k {
                        0 => p.u += 5e-7,
                        1 => p.x += 5e-7,
                        2 => p.y += 5e-7,
                        _ => p.gauss_size += 5e-7,
                    }
                    let up = get(&m2.breadths(0, hkl, wl).unwrap()).unwrap_or(0.0);
                    s += (up - get(&after).unwrap_or(0.0)).abs();
                }
                s
            };
            for (i, v, name, g) in [
                (9, after.size, "App-size", (|b: &drs::microstructure::Breadths| b.size) as fn(&drs::microstructure::Breadths) -> Option<f64>),
                (10, after.strain, "Max-strain", |b: &drs::microstructure::Breadths| b.strain),
            ] {
                let r: f64 = f[i].parse().unwrap();
                let v = match v {
                    Some(v) if v < 99999.0 => v,
                    _ => 99999.0,
                };
                let r_none = if name == "App-size" { r >= 99999.0 } else { r == 0.0 };
                let core_none = v >= 99999.0 || (name == "Max-strain" && after.strain.is_none());
                let e = if r_none || core_none {
                    if r_none == core_none { 0.0 } else { f64::INFINITY }
                } else {
                    (v - r).abs() / (half(f[i]) + spread(&g))
                };
                if e > worst_ss {
                    (worst_ss, at_ss) = (e, format!("{name} of {hkl:?} {v:.4} against {} (line {})", f[i], no + 1));
                }
            }
            n += 1;
        }
        println!("{run}: {n} rows of {code}_calc1.mic; widths and breadths at {worst:.2} of the tolerance ({at}); size and strain at {worst_ss:.2} ({at_ss})");
        if n == 0 || worst > 1.0 || worst_ss > 1.0 {
            fails.push(format!("{run}: .mic at {worst:.2} ({at}) and {worst_ss:.2} ({at_ss})"));
        }
    }
    println!("{rows_total} reflections of the .out files compared");
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}

#[test]
fn williamson_hall_axes() {
    let m = model("NdNi-HTO", "NdNi_HTO");
    let mic = latin1(&dir("NdNi-HTO").join("NdNi_HTO_calc1.mic"));
    let (mut n, mut worst, mut at) = (0, 0.0f64, String::new());
    for l in mic.lines() {
        let f: Vec<&str> = l.split_whitespace().collect();
        if f.len() < 20 || l.starts_with('!') || f[0].parse::<f64>().is_err() {
            continue;
        }
        let hkl = [f[11].parse::<i32>().unwrap(), f[12].parse().unwrap(), f[13].parse().unwrap()];
        let tt: f64 = f[14].parse().unwrap();
        let b = (1..=2).filter_map(|w| m.breadths(0, hkl, w)).min_by(|a, b| (a.two_theta - tt).abs().partial_cmp(&(b.two_theta - tt).abs()).unwrap()).unwrap();
        assert_eq!(b.hkl, hkl);
        for (i, v, name) in [(0, b.d_star, "s = 1/d"), (5, 1000.0 * b.beta_g_star, "betaG"), (6, 1000.0 * b.beta_l_star, "betaL")] {
            if f[i] == "NaN" {
                assert!(v.is_nan(), "{name} of {hkl:?}: NaN in the .mic");
                continue;
            }
            let r: f64 = f[i].parse().unwrap();
            let e = (v - r).abs() / (half(f[i]) + 2e-5 * r.abs());
            if e > worst {
                (worst, at) = (e, format!("{name} of {hkl:?}: {v:.6} against {}", f[i]));
            }
        }
        n += 1;
    }
    println!("NdNi-HTO: {n} reflections, d* and β* at {worst:.2} of the tolerance ({at})");
    assert!(n > 100 && worst <= 1.0, "{n} rows; {at}");
}
