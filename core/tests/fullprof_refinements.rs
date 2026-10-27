//! The least-squares engine against FullProf's converged refinements (test `fullprof_refinements`, L052, L053).
//!
//! FullProf refined each tutorial's starting file (`CODE.pcr.in`) to its own convergence (`CODE_conv.*`, L052). The
//! engine refines the same starting file, imported into the default model (L050) with the parameters of its codewords
//! (`refinable`), until every shift is below the `.pcr`'s Eps times its esd (FullProf's criterion), with FullProf's
//! reflection limit (S(1/d)max of `CODE_conv.out`). Compared:
//! - χ² over N − P + C against the last "Conventional Rietveld … Chi2" of `CODE_conv.out` (printed to 3 decimals);
//! - every refined parameter against FullProf's converged value (`CODE_conv.pcr`, imported into the default model the
//!   same way), within FullProf's esd (`CODE_conv.out`, translated: the X-ray scale ×2; P1 + 3P2 with σ₁ + 3σ₂, a
//!   bound; (S − D)/2 and (S + D)/2 with (σ_S + σ_D)/2);
//! - the esd's: the ratio of the engine's to FullProf's, printed.

use drs::jacobian::Field;
use drs::lsq::{refine, statistics, Options};
use drs::powder::{Intensities, Model, Shape};
use drs::refinable::ParameterSet;
use std::collections::BTreeMap;
use std::path::Path;

mod common;
use common::{fullprof_field, fullprof_params, latin1};

/// FullProf's last statistics over all non-excluded points (the block "RELIABILITY FACTORS WITH ALL NON-EXCLUDED
/// POINTS"; the next block is over the points with Bragg contributions only).
struct FpStats {
    dof: usize,
    residual: f64,
    chi2: f64,
    rp: f64,
    rwp: f64,
    dw: f64,
    rbragg: Vec<f64>,
}

fn fullprof_stats(out: &str) -> FpStats {
    let lines: Vec<&str> = out.lines().collect();
    let start = lines.iter().rposition(|l| l.contains("RELIABILITY FACTORS WITH ALL NON-EXCLUDED POINTS")).expect("statistics block");
    let num = |l: &str, key: &str, k: usize| -> f64 { l.split_once(key).unwrap().1.split_whitespace().nth(k).unwrap().parse().unwrap() };
    let mut s = FpStats { dof: 0, residual: 0.0, chi2: 0.0, rp: 0.0, rwp: 0.0, dw: 0.0, rbragg: Vec::new() };
    let mut i = start;
    while i < lines.len() && !lines[i].contains("FOR POINTS WITH BRAGG") {
        let l = lines[i];
        if l.contains("=> R-Factors:") {
            s.rp = num(l, "R-Factors:", 0);
            s.rwp = num(l, "R-Factors:", 1);
            s.dw = num(l, "DW-Stat.:", 0);
        } else if l.contains("=> N-P+C:") {
            s.dof = num(l, "N-P+C:", 0) as usize;
        } else if l.contains("SumYdif") {
            s.residual = lines[i + 1].split_whitespace().nth(4).unwrap().parse().unwrap();
        } else if l.contains("Conventional Rietveld") {
            s.chi2 = num(l, "Chi2:", 3);
        }
        i += 1;
    }
    s.rbragg = lines[i..].iter().filter_map(|l| l.split_once("Bragg R-factor:")).map(|(_, r)| r.split_whitespace().next().unwrap().parse().unwrap()).collect();
    s
}

fn load(d: &Path, file: &str, code: &str, data: Option<&Path>) -> (drs::pcr::Pcr, drs::powder::Loaded) {
    let mut p = drs::pcr::read(&latin1(&d.join(file))).unwrap();
    p.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
    let l = Model::from_pcr_with_data(&p, d, code, data).unwrap_or_else(|e| panic!("{file}: {e}"));
    (p, l)
}

const RUNS: [(&str, &str); 9] = [
    ("Y2O3", "y2o3_structure"),
    ("Yb2ZrWO8", "ybfluo_Riet_IRF"),
    ("MgAl2O4-s", "MgAl2O4s_Rietveld"),
    ("LaB6", "LaB6rx"),
    ("Egypt", "Egypt_qa_Rietveld"),
    ("rutile-anatase", "rutana_tch"),
    ("Co-YBCO", "Cocus-asym"),
    ("NdNi-HTO", "NdNi_HTO"),
    ("YVO3", "YVO3a"),
];

/// χ² tolerance, as in `fullprof_profiles`: FullProf's print (5·10⁻⁴) and 0.15 % of χ².
fn chi2_tol(chi2: f64) -> f64 {
    5e-4 + 1.5e-3 * chi2
}

/// Runs whose function differs from FullProf's at FullProf's converged values: frozen |χ² − χ²_FullProf|/χ²_FullProf.
/// - LaB6: FullProf's peak windows in its 8-cycle run are not those of its 1-cycle run at the same values (its own
///   `LaB6rx_calc.prf` and `LaB6rx_conv.prf`, from identical parameters, give Σw(y − yc)² 107888 and 108157); neither
///   Wdt·H_noX nor Wdt·H reproduces the second; the rule was not found.
/// - YVO3: the window of the parametric resolution file (L039, L041), frozen in `fullprof_profiles`.
/// - NdNi-HTO: its magnetic phase is not calculated.
const FUNCTION_FROZEN: &[(&str, f64)] = &[("LaB6", 0.0064), ("YVO3", 0.0256), ("NdNi-HTO", 10.13)];

/// Runs whose refined χ² differs from FullProf's: frozen (χ² − χ²_FullProf)/χ²_FullProf, signed.
/// - Co-YBCO: FullProf's last point is not a minimum of its own χ² (the function agrees there to 10⁻⁴ and the engine,
///   started from it, lowers χ² from 1.5239 to 1.4834); its 18 cycles show χ² rising from 1.56 to 1.74 in cycles 1–3
///   and falling by 0.01 per cycle at the end.
/// - LaB6, YVO3, NdNi-HTO: the function differs (above). YVO3 also stops FullProf without convergence (L052).
const ENGINE_FROZEN: &[(&str, f64)] = &[("Co-YBCO", -0.0268), ("LaB6", 0.0025), ("YVO3", -0.3764), ("NdNi-HTO", 0.0877)];

/// Refined parameters more than one FullProf esd away from FullProf's value: frozen (count, largest in esd's).
/// Co-YBCO and NdNi-HTO for the reasons above; LaB6's cell (1.29 esd of FullProf's esd, 0.71 of the engine's, which
/// is twice FullProf's for the cell, see below).
const PARAMS_FROZEN: &[(&str, usize, f64)] = &[("Co-YBCO", 17, 25.7), ("LaB6", 1, 1.30), ("NdNi-HTO", 11, 67.4)];

/// The esd's, compared as q = (σ_engine/σ_FullProf)/√(χ²_all/χ²_Bragg): FullProf scales by the χ² over the points with
/// Bragg contributions (manual eqs. 3.7–3.8), the engine by that over all the points used, so q = 1 for most
/// parameters, within 2 %. Kinds that differ, frozen:
/// - the cell parameters and I_G of a Thompson–Cox–Hastings phase (Npr = 7): q = 2 (1.98–2.02, in the five runs
///   compared): FullProf prints half of √(C_kk χ²), where the engine's derivative is checked against finite
///   differences of the calculated pattern (`engine_derivatives`; at Yb2ZrWO8's values with FullProf's rules too,
///   Σw(∂y/∂a)² to 10⁻⁶). With the pseudo-Voigt (Y2O3, Npr = 5) q = 1.001. Why FullProf's is half is not known;
/// - Bérar–Baldinozzi P1 + 3P2: FullProf's σ₁ + 3σ₂ is an upper bound (q ≤ 1.02);
/// - the parameters named in `ESD_FROZEN`, by their q.
const ESD_FROZEN: &[(&str, &str, f64)] = &[
    // The ends of a linearly interpolated background.
    ("Y2O3", "pattern1.background_points.0", 1.050),
    ("Y2O3", "pattern1.background_points.61", 0.682),
    ("Yb2ZrWO8", "pattern1.background_points.0", 1.250),
    ("Yb2ZrWO8", "pattern1.background_points.1", 1.086),
    ("MgAl2O4-s", "pattern1.background_points.0", 1.028),
    ("MgAl2O4-s", "pattern1.background_points.31", 0.671),
    // LaB6's wavelength ratio (VARY ratio) and the scale it correlates with.
    ("LaB6", "pattern1.wavelength_ratio", 0.562),
    ("LaB6", "phase1.scale_factor", 0.475),
];

/// FullProf's χ² over the points with Bragg contributions (the second statistics block of the last cycle).
fn fullprof_chi2_bragg(out: &str) -> f64 {
    let lines: Vec<&str> = out.lines().collect();
    let start = lines.iter().rposition(|l| l.contains("FOR POINTS WITH BRAGG CONTRIBUTIONS")).expect("Bragg block");
    lines[start..].iter().find_map(|l| l.split_once("Conventional Rietveld Rp,Rwp,Re and Chi2:")).map(|(_, r)| r.split_whitespace().nth(3).unwrap().parse().unwrap()).unwrap()
}

#[test]
fn fullprof_refinements() {
    let mut fails: Vec<String> = Vec::new();
    for (run, code) in RUNS {
        let d = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs").join(run);
        // The converged file names the data FullProf was run on (the starting files may name another, L052).
        let conv_file = if run == "YVO3" { format!("{code}_conv.new") } else { format!("{code}_conv.pcr") };
        let (_, conv) = load(&d, &conv_file, code, None);
        let (p, start) = load(&d, &format!("{code}.pcr.in"), code, conv.data_path.as_deref());
        let out = latin1(&d.join(format!("{code}_conv.out")));
        let mut m = start.model.clone();
        let lim: Vec<f64> = out.lines().filter_map(|l| l.split_once("S(1/d)max:")).map(|(_, r)| r.split_whitespace().next().unwrap().parse().unwrap()).collect();
        if lim.len() >= m.phases.len() {
            for (ph, l) in m.phases.iter_mut().zip(&lim) {
                ph.reflection_limit = Some(*l);
            }
        }
        // FullProf's windows in its refinement cycles: rutile-anatase's 11-cycle run cuts each peak at Wdt·H, H with
        // X tan θ (its `.prf` at the converged values: Σw(y − yc)² 10552.2, the core with Wdt·H 10551.5, with its
        // 1-cycle rule Wdt·H_noX 10810.3), the default model's window. Its 1-cycle runs cut at Wdt·H_noX (L041).
        if run == "rutile-anatase" {
            m.convention.peak_window = false;
        }
        let set = ParameterSet::from_pcr(&p, &m).unwrap_or_else(|e| panic!("{run}: {e}"));
        let data = start.data.as_ref().unwrap();
        let t = std::time::Instant::now();
        let r = refine(&mut m, data, &set, Options { max_cycles: 100, eps: p.eps }, |_| {});
        let fp = fullprof_stats(&out);
        println!(
            "{run:15} {} parameters, {} cycles ({}; {:.0} s): χ² {:.4} over N − P = {}, FullProf {:.4} over N − P + C = {}; Rp {:.2} Rwp {:.2} DW {:.4} (FullProf {:.2} {:.2} {:.4}); R_Bragg {:?} (FullProf {:?})",
            set.free.len(),
            r.cycles,
            r.stop,
            t.elapsed().as_secs_f64(),
            r.stats.chi2,
            r.stats.n - r.stats.p,
            fp.chi2,
            fp.dof,
            r.stats.rp,
            r.stats.rwp,
            r.stats.durbin_watson,
            fp.rp,
            fp.rwp,
            fp.dw,
            r.stats.phases.iter().map(|q| format!("{:.2}", q.0)).collect::<Vec<_>>(),
            fp.rbragg
        );
        // Egypt: FullProf counts one point more than it uses (7001 points, 171 in 5.00–6.70°, and 75.00°, the last,
        // which its `.prf` marks excluded too: 6829 − 44 = 6785, FullProf 6786); believed to be that last point.
        if r.stats.n - r.stats.p + usize::from(run == "Egypt") != fp.dof {
            fails.push(format!("{run}: N − P {} against FullProf's N − P + C {}", r.stats.n - r.stats.p, fp.dof));
        }
        // The same function at FullProf's converged values, and the engine started there.
        let mut at_fp = conv.model.clone();
        at_fp.convention = m.convention;
        for (ph, q) in at_fp.phases.iter_mut().zip(&m.phases) {
            ph.reflection_limit = q.reflection_limit;
        }
        let s_fp = statistics(&at_fp, data, &at_fp.calculate(&data.x).ycalc, set.free.len());
        let r_fp = refine(&mut at_fp, data, &set, Options { max_cycles: 100, eps: p.eps }, |_| {});
        println!("    χ² at FullProf's values {:.4}; the engine from there {:.4} in {} cycles ({})", s_fp.chi2, r_fp.stats.chi2, r_fp.cycles, r_fp.stop);
        let df = (s_fp.chi2 - fp.chi2).abs();
        match FUNCTION_FROZEN.iter().find(|f| f.0 == run) {
            Some(f) if df / fp.chi2 > f.1 => fails.push(format!("{run}: χ² at FullProf's values {:.4} against {}: {:.4} of it, frozen {}", s_fp.chi2, fp.chi2, df / fp.chi2, f.1)),
            None if df > chi2_tol(fp.chi2) => fails.push(format!("{run}: χ² at FullProf's values {:.4} against {}", s_fp.chi2, fp.chi2)),
            _ => {}
        }
        let de = r.stats.chi2 - fp.chi2;
        match ENGINE_FROZEN.iter().find(|f| f.0 == run) {
            Some(f) if (de / fp.chi2).abs() > f.1.abs() || de.signum() != f.1.signum() => {
                fails.push(format!("{run}: refined χ² {:.4} against FullProf's {}: {:+.4} of it, frozen {:+}", r.stats.chi2, fp.chi2, de / fp.chi2, f.1))
            }
            None if de.abs() > chi2_tol(fp.chi2) => fails.push(format!("{run}: refined χ² {:.4} against FullProf's {}", r.stats.chi2, fp.chi2)),
            _ => {}
        }
        // Where FullProf's point is not the minimum, the engine reaches the same χ² from both starts.
        if run == "Co-YBCO" && (r_fp.stats.chi2 - r.stats.chi2).abs() > chi2_tol(r.stats.chi2) {
            fails.push(format!("{run}: the engine from FullProf's values reaches {:.4}, from the start {:.4}", r_fp.stats.chi2, r.stats.chi2));
        }
        // FullProf's esd's in the default model's fields.
        let labels: Vec<Vec<String>> = p.phases.iter().map(|ph| ph.atoms.iter().map(|a| a.label.clone()).collect()).collect();
        let fpp = fullprof_params(&out);
        let mut sig: BTreeMap<Field, f64> = BTreeMap::new();
        for (name, &(_, e)) in &fpp {
            if let Some(f) = fullprof_field(name, &labels) {
                sig.insert(f, e);
            }
        }
        let translated = |f: Field| -> Option<f64> {
            match f {
                Field::Scale(q) if !m.neutrons && matches!(m.phases[q].intensities, Intensities::Structure(..)) => sig.get(&f).map(|e| 2.0 * e),
                Field::Asym(q, k) if k % 2 == 0 => Some(sig.get(&f).copied().unwrap_or(0.0) + 3.0 * sig.get(&Field::Asym(q, k + 1)).copied().unwrap_or(0.0)).filter(|&e| e > 0.0),
                _ => sig.get(&f).copied(),
            }
        };
        let scaling = (fp.chi2 / fullprof_chi2_bragg(&out)).sqrt();
        let (mut beyond, mut worst, mut qs) = (0, (0.0f64, String::new()), Vec::new());
        for (fpar, ours) in set.free.iter().zip(&r.params) {
            let f = fpar.fields[0].0;
            let (v, vf) = (m.get(f), conv.model.get(f));
            let Some(e) = translated(f) else {
                fails.push(format!("{run}: {} has no FullProf esd", fpar.name));
                continue;
            };
            let dev = (v - vf).abs() / e;
            if dev > worst.0 {
                worst = (dev, fpar.name.clone());
            }
            if dev > 1.0 {
                beyond += 1;
                println!("    {:40} {v:14.6} ({:.1e}) FullProf {vf:14.6} ({e:.1e}): {dev:.2} esd{}", fpar.name, ours.esd, if ours.at_limit { ", held at its limit" } else { "" });
            }
            // The esd's, where the function and the minimum are FullProf's.
            if FUNCTION_FROZEN.iter().any(|x| x.0 == run && x.0 != "LaB6") || ENGINE_FROZEN.iter().any(|x| x.0 == run && x.0 != "LaB6") {
                continue;
            }
            let q = ours.esd / e / scaling;
            qs.push(q);
            let ok = match (f, ESD_FROZEN.iter().find(|x| x.0 == run && x.1 == fpar.name)) {
                (_, Some(x)) => (q - 1.0).abs() <= (x.2 - 1.0f64).abs() + 1e-3,
                (Field::Cell(k, _) | Field::GaussSize(k), _) if m.phases[k].shape == Shape::ThompsonCoxHastings => (1.98..=2.02).contains(&q),
                (Field::Asym(..), _) => q <= 1.02,
                _ => (0.98..=1.02).contains(&q),
            };
            if !ok {
                fails.push(format!("{run}: {} esd {:.3e} against FullProf's {e:.3e}: q = {q:.3}", fpar.name, ours.esd));
            }
        }
        qs.sort_by(f64::total_cmp);
        println!(
            "    {beyond} parameters beyond 1 esd, the largest {:.3} ({}); esd's: q = {} (FullProf's scaled by {scaling:.4})",
            worst.0,
            worst.1,
            if qs.is_empty() { "not compared".to_string() } else { format!("{:.3} … {:.3}, median {:.3}", qs[0], qs[qs.len() - 1], qs[qs.len() / 2]) }
        );
        // Quantitative phase analysis (P33): the weight fractions against FullProf's "Fract(%)" of its last cycle,
        // each within FullProf's esd. The esd's: FullProf's are reproduced by treating the numerator and the sum as
        // independent, σ²(W_φ) = W_φ²[(σ_φ/S_φ)² + Σ_ψ W_ψ²(σ_ψ/S_ψ)²] (rutile: 0.521 % against the printed 0.52), so
        // a single phase gets an esd (Y2O3: 100 ± 0.23 %). The core propagates the covariance, which gives 0 for one
        // phase and smaller esd's where the scales are correlated with the sum; FullProf's formula is checked here
        // from the engine's scale esd's (scaled as above), within 2 % and the print.
        let fract: Vec<(f64, f64)> = out
            .lines()
            .filter_map(|l| l.split_once("Fract(%):"))
            .map(|(_, r)| {
                let v: Vec<f64> = r.replace(['(', ')'], " ").split_whitespace().filter_map(|t| t.parse().ok()).collect();
                (v[0], v[1])
            })
            .collect();
        if !fract.is_empty() {
            let last = &fract[fract.len() - m.phases.len()..];
            let ours = drs::qpa::weight_fractions_esd(&m, &set, &r.covariance);
            for (k, (o, &(wf, ef))) in ours.iter().zip(last).enumerate() {
                let Some((w, e)) = o else {
                    fails.push(format!("{run}: phase {} has no weight fraction", k + 1));
                    continue;
                };
                let (w, e) = (100.0 * w, 100.0 * e);
                // FullProf's formula from the engine's relative scale esd's.
                let rel = |q: usize| -> f64 {
                    let a = drs::refinable::address(Field::Scale(q), &m);
                    r.fields.iter().find(|x| x.0 == a).map_or(0.0, |x| x.2 / x.1)
                };
                let wall: Vec<f64> = ours.iter().map(|o| o.map_or(0.0, |x| x.0)).collect();
                let sum: f64 = (0..m.phases.len()).map(|q| (wall[q] * rel(q)).powi(2)).sum();
                let e_fp = w * (rel(k).powi(2) + sum).sqrt();
                println!("    phase {} weight fraction {w:.3} ({e:.3}) %, FullProf {wf} ({ef}): {:.2} esd; FullProf's formula {e_fp:.3}, scaled {:.3}", k + 1, (w - wf).abs() / ef, e_fp / scaling);
                // FullProf prints the fraction to 0.01 % and its esd to 0.01 %: half a digit is added to each.
                if (w - wf).abs() > ef + 0.005 {
                    fails.push(format!("{run}: phase {} weight fraction {w:.3} % against FullProf's {wf} ± {ef}", k + 1));
                }
                // Where the esd's are compared (the same minimum and function, LaB6's scale esd frozen above).
                let compared = !FUNCTION_FROZEN.iter().any(|x| x.0 == run) && !ENGINE_FROZEN.iter().any(|x| x.0 == run);
                if compared && (e_fp / scaling - ef).abs() > 0.02 * ef + 0.005 {
                    fails.push(format!("{run}: phase {} weight fraction esd by FullProf's formula {:.4} % against FullProf's {ef}", k + 1, e_fp / scaling));
                }
                if e > e_fp + 1e-9 {
                    fails.push(format!("{run}: phase {} the propagated esd {e:.4} % exceeds FullProf's formula's {e_fp:.4}", k + 1));
                }
            }
        }
        if run == "YVO3" {
            // Not compared: FullProf does not converge (L052). Reported: the engine holds U at its limit 0 where
            // FullProf takes it to −0.091.
            for ps in r.params.iter().filter(|q| q.at_limit) {
                println!("    held at its limit: {} = {}", ps.name, ps.value);
            }
            if !r.params.iter().any(|q| q.at_limit && q.name == "phase1.caglioti_u" && q.value == 0.0) {
                fails.push(format!("{run}: U is not held at its limit 0"));
            }
            continue;
        }
        match PARAMS_FROZEN.iter().find(|x| x.0 == run) {
            Some(x) if beyond > x.1 || worst.0 > x.2 => fails.push(format!("{run}: {beyond} parameters beyond 1 esd, the largest {:.2}, frozen {} and {}", worst.0, x.1, x.2)),
            None if beyond > 0 => fails.push(format!("{run}: {beyond} parameters beyond 1 esd, the largest {:.3} ({})", worst.0, worst.1)),
            _ => {}
        }
    }
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}
