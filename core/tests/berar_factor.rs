//! Bérar & Lelann's factor (P43, L244) against the one FullProf prints in its `.out` files ("Standard deviations have
//! to be multiplied by:", which cites J.-F. Bérar & P. Lelann, J. Appl. Cryst. 24, 1–5, 1991, in `manuals/papers/`).
//!
//! For every `.out` in `manuals/` that prints the factor (the patterns of a joint run each, 176 in all), the
//! core's `lsq::berar_factor` is computed from FullProf's own Ycal (its `.prf`) and the data as FullProf reads them
//! (`Model::from_pcr_fullprof`), with the weights and the points with Bragg contributions of the core's FullProf
//! rules. A run is compared when those are FullProf's: Σw(y − yc)² over all the points and over the points with Bragg contributions equal the "Residual" of the
//! last cycle's two statistics blocks of its `.out` within 0.05 %. The others are left out by reason, counted:
//! - no `.prf` beside the `.out`;
//! - the model is not calculated by the core (a magnetic phase, a size or strain model; L211);
//! - the weights differ (S over all the points is not FullProf's: data read differently);
//! - the points with Bragg contributions differ (the core's peak windows are not FullProf's: the window rules not
//!   found in L042, frozen in `fullprof_profiles`).
//!
//! Compared runs agree within 0.05 %; or FullProf's value lies within the range the core's factor takes when the
//! `.prf`'s Ycal moves by up to half a unit of its last printed decimal (40 draws: the Y₂O₃ runs print Ycal to one
//! decimal), counted; or they are frozen in `FROZEN` with their deviation and its reason.

use std::path::{Path, PathBuf};

mod common;
use common::latin1;

/// Compared runs beyond 0.05 % and outside the .prf's rounding: (run, the largest |deviation| allowed, in %), as
/// measured in L244:
/// - the narrow-window experiments of L042: the core's peak windows differ from FullProf's at their edges there (the
///   rules not found, frozen in `fullprof_profiles`), so the points with Bragg contributions differ in a few points
///   though S over them agrees within 0.05 %; the factor depends on which runs of residuals they cut;
/// - LBF_lmopr (a Le Bail fit, Ycal printed to three decimals) and sr-ox-rest_conv (a restrained refinement): not
///   explained; listed for Josu;
/// - prba10k (L267): the window edges of its needle-and-strain phase, as the narrow windows above.
const FROZEN: &[(&str, f64)] = &[
    ("fullprof/7-FullProf-Tutorial-SymmetryModes/LaMnO3/LBF_lmopr.out", 0.087),
    ("fullprof/runs/Sr-ox-rest/sr-ox-rest_conv.out", 0.880),
    ("fullprof/runs/window_experiments/wdt10_x0_narrow/ybfluo_Riet_IRF.out", 14.965),
    ("fullprof/runs/window_experiments/wdt20_x0_narrow/ybfluo_Riet_IRF.out", 12.617),
    ("fullprof/runs/window_experiments/wdt40_x0.645_narrow/ybfluo_Riet_IRF.out", 0.091),
    ("fullprof/runs/window_experiments/wdt40_x0_narrow/ybfluo_Riet_IRF.out", 0.086),
    ("fullprof/runs/window_experiments/yb_irf1_x0.645/ybfluo_Riet_IRF.out", 6.024),
    // L267: prba10k's windows end within half a step of FullProf's on either side (592 edge points frozen in
    // fullprof_profiles); both passes print the same factor.
    ("fullprof/runs/prba10k/prba10k.out", 0.156),
    ("fullprof/runs/prba10k/prba10k_calc.out", 0.156),
];

/// (printing the factor, compared, of them within the .prf's rounding, no .prf, not calculated, other weights, other
/// Bragg points). L267 added FullProf's runs of prba10k, cuf1k and lac (6 outputs): prba10k's two are compared, and
/// Pd3MnD's two, calculated since; cuf1k's and lac's four are not (their magnetic phases). L322 added the joint runs
/// pbb, scor_rxn and prbaco (6 outputs, 12 patterns): 10 patterns compared, all within 0.016 %; prbaco's neutron
/// pattern in both outputs has other weights (its Σ 0.13% from FullProf's, frozen in `joint_refinement`). L334–L341
/// added P05's TOF runs (27 outputs, gabriela's three patterns each): arg_si_invd compared (−0.007 %); arg_si_abs2 and
/// abs2b without a .prf (FullProf stopped, "Singular matrix"); arg_si_plain not calculated (FullProf rewrote its .pcr
/// with NaN); 13 with other weights, each a χ² difference frozen in `fullprof_profiles` (CeO2_PEARL's RALF σ, 2;
/// gabriela's σ = 0 points, 6; lamn_pol's and nac-osiris' Iwg = 1, 5); 10 with other Bragg points (arg_si and six of
/// its variants, alum and alum_noq's: FullProf's TOF window rule, not found, frozen there).
/// L353 added L346's runs (42 outputs): the 16 of background_terms but y2o3_filter compared, every one within 0.05 %
/// or within the .prf's rounding (y2o3_poly12_b2's −1.5 % among them); ureaz's two; PPH3-sls's two with other weights
/// (S 1 192 400 against FullProf's 1 268 000, as the tutorial's own PPH3-sls.out in this count); the 15 with TLS not
/// calculated (P41's next step, L367) and the four of texture/ not read (P30). L354 added L351's citric runs (3 outputs:
/// citric-acid-sls_1's two passes, citric_core), each with other weights, as the tutorial's own citric-acid-sls_1.out.
/// L360 reads texture/'s four (P30): compared, within 0.005 % (yb_nor3's +0.0045 % among them).
/// L367 calculates the 18 outputs of rigid/'s TLS runs (P41's TLS): 13 compared (pdd's five and urea_conv, urea_test,
/// ureac, ureacyl with their _calc), within 0.05 % or, ureac.out (−0.19 %), the .prf's rounding; 5 with other Bragg
/// points (hmt's two, urea's, urea_calc's and urea_notls'), S over all the points agreeing.
const COUNTS: (usize, usize, usize, usize, usize, usize, usize) = (221, 113, 29, 6, 15, 31, 56);

fn outs(d: &Path, v: &mut Vec<PathBuf>) {
    for e in std::fs::read_dir(d).unwrap() {
        let p = e.unwrap().path();
        if p.is_dir() {
            outs(&p, v);
        } else if p.extension().is_some_and(|x| x == "out") {
            v.push(p);
        }
    }
}

/// The "Residual" of the last two statistics blocks of a `.out` (of pattern k of the n a joint run's last cycle
/// reports in turn): all the points, then those with Bragg contributions.
fn residuals(out: &str, n: usize, k: usize) -> Option<(f64, f64)> {
    let ls: Vec<&str> = out.lines().collect();
    let idx: Vec<usize> = ls.iter().enumerate().filter(|(_, l)| l.contains("SumYdif")).map(|(i, _)| i).collect();
    let r = |k: usize| -> Option<f64> { ls.get(k + 1)?.split_whitespace().nth(4)?.parse().ok() };
    let first = idx.len().checked_sub(2 * n)? + 2 * k;
    Some((r(idx[first])?, r(idx[first + 1])?))
}

#[test]
fn berar_factor() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals");
    let mut files = Vec::new();
    outs(&root, &mut files);
    files.sort();
    let (mut printed, mut no_prf, mut not_calculated, mut weights_differ, mut bragg_differ) = (0, 0, 0, 0, 0);
    let mut compared = Vec::new();
    for out in files {
        let text = latin1(&out);
        let factors: Vec<f64> = text.lines().filter_map(|l| l.split_once("Standard deviations have to be multiplied by:")).filter_map(|(_, r)| r.trim().parse::<f64>().ok()).collect();
        if factors.is_empty() {
            continue;
        }
        let file = out.strip_prefix(&root).unwrap().display().to_string();
        let stem = out.file_stem().unwrap().to_str().unwrap().to_string();
        if !out.with_extension("pcr").exists() {
            // A .out without its .pcr: counted as before, without a .prf (asserted).
            printed += 1;
            assert!(!out.with_extension("prf").exists(), "{file}: a .prf without its .pcr");
            no_prf += 1;
            continue;
        }
        let pcr_text = latin1(&out.with_extension("pcr"));
        // A joint run's .out (P06, L322) prints a factor and two statistics blocks per pattern, in order, and its
        // .prf is one per pattern (CODE_k.prf): each pattern is compared as a run of its own ("CODE.out#k").
        // A file the reader refuses: not calculated (none since P30's step read runs/texture, L360).
        let pcrs = match drs::pcr::read_all(&pcr_text) {
            Ok(p) => p,
            Err(e) if e.message.contains("not read yet") => {
                println!("{file:<62} left out: not read ({e})");
                printed += 1;
                not_calculated += 1;
                continue;
            }
            Err(e) => panic!("{file}: {e}"),
        };
        let n = pcrs.len();
        for (k, pcr) in pcrs.into_iter().enumerate() {
        let fp = factors[factors.len() - n + k];
        printed += 1;
        let (name, prf) = if n == 1 { (file.clone(), out.with_extension("prf")) } else { (format!("{file}#{}", k + 1), out.with_file_name(format!("{stem}_{}.prf", k + 1))) };
        if !prf.exists() {
            println!("{name:<62} left out: no .prf");
            no_prf += 1;
            continue;
        }
        let p = drs::formats::fullprof_outputs::read_prf(&latin1(&prf)).unwrap_or_else(|e| panic!("{name}: {e:?}"));
        let l = match drs::powder::Model::from_pcr_fullprof(&pcr, out.parent().unwrap(), &stem) {
            Ok(l) => l,
            Err(e) => {
                println!("{name:<62} left out: not calculated ({e})");
                not_calculated += 1;
                continue;
            }
        };
        let d = l.data.as_ref().unwrap_or_else(|| panic!("{name}: no data"));
        let m = &l.model;
        // FullProf's Ycal at the data's points (the .prf lists the same points, the excluded ones too; a TOF .prf
        // prints x in single precision, 2·10⁻³ µs at 24 000 µs: matched within 10⁻⁴ + 10⁻⁷|x|, as `fullprof_profiles`,
        // L341; 2θ within 10⁻⁴).
        let tof = l.model.tof.is_some();
        let ycalc: Vec<f64> = d
            .x
            .iter()
            .map(|&x| {
                // (The search for a 2θ pattern starts 10⁻⁵ below x, as before L341: from 10⁻⁴ below, five more
                // outputs would be compared, PPH3-sls at −1.77 %; left for a later step.)
                let (tol, from) = if tof { (1e-4 + 1e-7 * x.abs(), 1e-4 + 1e-7 * x.abs()) } else { (1e-4, 1e-5) };
                let j = p.x.partition_point(|&t| t < x - from);
                if j < p.x.len() && (p.x[j] - x).abs() <= tol {
                    p.ycalc[j]
                } else {
                    f64::NAN
                }
            })
            .collect();
        let w = drs::lsq::weights(m, d, &ycalc);
        let peaks = m.calculate(&d.x).peaks;
        let mut bragg = vec![false; d.x.len()];
        for pk in &peaks {
            let lo = d.x.partition_point(|&t| t < pk.two_theta - pk.reach);
            let hi = d.x.partition_point(|&t| t <= pk.two_theta + pk.reach);
            bragg[lo..hi.max(lo)].iter_mut().for_each(|b| *b = true);
        }
        let term = |i: usize| if w[i] > 0.0 && ycalc[i].is_finite() { w[i] * (d.y[i] - ycalc[i]).powi(2) } else { 0.0 };
        let s_all: f64 = (0..w.len()).map(term).sum();
        let s_bragg: f64 = (0..w.len()).filter(|&i| bragg[i]).map(term).sum();
        let (r_all, r_bragg) = residuals(&text, n, k).unwrap_or_else(|| panic!("{name}: no statistics blocks"));
        if (s_all / r_all - 1.0).abs() > 5e-4 {
            println!("{name:<62} left out: other weights (S {s_all:.6e}, FullProf {r_all:.6e})");
            weights_differ += 1;
            continue;
        }
        if (s_bragg / r_bragg - 1.0).abs() > 5e-4 {
            println!("{name:<62} left out: other Bragg points (S {s_bragg:.6e}, FullProf {r_bragg:.6e})");
            bragg_differ += 1;
            continue;
        }
        let f = drs::lsq::berar_factor(m, d, &ycalc, &peaks);
        let dev = 100.0 * (f / fp - 1.0);
        // The factor for Ycal within half a unit of its print: does FullProf's value lie in that range?
        let rounding = dev.abs() > 0.05 && {
            let h = 0.5 * 10f64.powi(-(p.ycalc_decimals as i32));
            let mut seed: u64 = 7;
            let (mut lo, mut hi) = (f, f);
            for _ in 0..40 {
                let y: Vec<f64> = ycalc
                    .iter()
                    .map(|&v| {
                        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                        v + (2.0 * ((seed >> 11) as f64 / (1u64 << 53) as f64) - 1.0) * h
                    })
                    .collect();
                let g = drs::lsq::berar_factor(m, d, &y, &peaks);
                (lo, hi) = (lo.min(g), hi.max(g));
            }
            fp >= lo - 5e-5 && fp <= hi + 5e-5
        };
        println!("{name:<62} FullProf {fp:.4}  core {f:.4}  {dev:+.4} %{}", if rounding { "  (within the .prf's rounding)" } else { "" });
        compared.push((name, dev, rounding));
        }
    }
    let within_rounding = compared.iter().filter(|c| c.2).count();
    println!("{printed} outputs print the factor: {} compared ({within_rounding} beyond 0.05 % within the .prf's rounding), {no_prf} without a .prf, {not_calculated} not calculated, {weights_differ} with other weights, {bragg_differ} with other Bragg points", compared.len());
    let mut beyond = Vec::new();
    for (name, dev, rounding) in &compared {
        if dev.abs() > 0.05 && !rounding {
            match FROZEN.iter().find(|f| f.0 == name) {
                Some(&(_, allowed)) => assert!(dev.abs() <= allowed, "{name}: {dev:+.4} %, frozen at {allowed} %"),
                None => beyond.push(format!("{name} {dev:+.4} %")),
            }
        }
    }
    assert!(beyond.is_empty(), "beyond 0.05 % and not frozen:\n{}", beyond.join("\n"));
    assert_eq!(compared.iter().filter(|c| c.1.abs() > 0.05 && !c.2).count(), FROZEN.len(), "a frozen run now agrees: unfreeze it");
    assert_eq!((printed, compared.len(), within_rounding, no_prf, not_calculated, weights_differ, bragg_differ), COUNTS, "the counts");
}

/// The factor of a run: `Refinement::berar_factor`, in the results file (`_drs_berar_factor`, read back) and printed by
/// `drs refine` (case A, from FullProf's converged state).
#[test]
fn berar_factor_in_the_record() {
    let src = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs/modes/A_modes");
    let d = Path::new(env!("CARGO_MANIFEST_DIR")).join("target/test-berar");
    std::fs::create_dir_all(&d).unwrap();
    for f in ["a_modes.pcr", "SrScLRRT.dat"] {
        std::fs::copy(src.join(f), d.join(f)).unwrap();
    }
    let pcr = drs::pcr::read(&latin1(&d.join("a_modes.pcr"))).unwrap();
    let c = drs::control::Control::from_pcr(&pcr, &d, "a_modes", &d).unwrap();
    let l = drs::powder::Model::from_pcr(&pcr, &d, "a_modes").unwrap();
    let data = l.data.unwrap();
    let mut m = c.model.clone();
    let set = c.refinement.build(&m).unwrap();
    let r = drs::lsq::refine(&mut m, &data, &set, drs::lsq::Options { max_cycles: 1, eps: 0.1 }, |_| {});
    let yc = m.calculate(&data.x);
    assert_eq!(r.berar_factor, drs::lsq::berar_factor(&m, &data, &yc.ycalc, &yc.peaks));
    assert!(r.berar_factor > 3.0 && r.berar_factor < 4.0, "{}", r.berar_factor);
    let mut cc = c.clone();
    cc.model = m.clone();
    let (text, _) = drs::results::write(&cc, &data, Some(&set), Some(&r), "run 2026-10-15 12.00.00").unwrap();
    assert_eq!(drs::results::read(&text).unwrap().berar_factor.map(|f| (f * 1e4).round()), Some((r.berar_factor * 1e4).round()));
    let o = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(["refine", d.join("a_modes.pcr").to_str().unwrap(), "--cycles", "1"]).output().unwrap();
    let out = String::from_utf8_lossy(&o.stdout);
    assert!(o.status.success() && out.contains("Bérar factor 3."), "{out}");
}
