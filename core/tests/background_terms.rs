//! The further terms of the background (P23, L353) against FullProf 8.20's runs of `runs/version2_runs.py`
//! (`manuals/fullprof/runs/background_terms`, run on the Mac: L346).
//!
//! - `against_fullprof_backgrounds`: the core's background at the points of each run's `.prf` against FullProf's
//!   Backg, under FullProf's rules, within the precision FullProf prints it to, by the bound `background_points`
//!   derives (one offset c with |ours − c − printed| ≤ ½·10^(−decimals) plus single precision at every point, L081),
//!   widened by FullProf's single-precision evaluation of the function (`single_precision`, derived), and the offset a
//!   quarter of the least background (L082, L087):
//!   - Nba = −1 (manual eq. 3.4, Q = 4π sin θ/λ₁ at the data's 2θ): nac_3t2_calc and nac_3t2_debye (its Debye-like
//!     terms alone), 2 508 of 2 510 points; the first two, in the excluded region, frozen (`RUNS`);
//!   - Nba = −3: nine coefficients of the polynomial, then c₁/2θ + c₂/2θ² + c₃/2θ³, 2θ in degrees (y2o3_poly12_b1, b2,
//!     b7, b9, b10, b11, b12, each with one coefficient besides B1 = 20);
//!   - Nba = 1 with b1–b4 (manual eq. 3.6): y2o3_bac_b1b2 (b1 = 0.8, b2 = 5) and y2o3_bac_b3b4 (b3 = 0.01, b4 = 0.5);
//!   - MORE = 2, 3 (the News of 28 July 2003): y2o3_more2, y2o3_more3 and y2o3_more2_two (y2o3_1.bac × 0.5,
//!     y2o3_2.bac × 2);
//!   - FullProf's background is 0 where the function falls below 0 or above the data's largest count
//!     (`Convention::background_bounds`): nac_3t2_debye (1 504 points below 0) and y2o3_poly12_b9 (96 points above
//!     Y2O3's 11 248 counts, from 144.95°).
//! - `fourier_filter_refused`: Nba = −2 (y2o3_filter, FWINDOW 3000) is refused by the import: the filtering of manual
//!   eq. 3.7 (a smoothing of Yobs − Ybragg by FFT) did not reproduce FullProf's background from its `.prf` (L353:
//!   Numerical Recipes' SMOOFT at FWINDOW and at other windows, rms 0.5–1.9 counts against corrections up to 6.3).
//! - `debye_refinement`: nac_3t2's refinement (50 parameters, the polynomial's six and five Debye-like amplitudes among
//!   them) from FullProf's start, to the engine's convergence, against FullProf's converged values in its `.out`.
//! - `control_file_and_export`: five runs' control files, written in another folder, read back to the same model,
//!   pattern and free parameters (the profiles' files named relative to the control file); every further term located
//!   in the text; the `.pcr` export (Nba = −1, −3, 1, MORE = 2 with the profiles written at the data's points) read
//!   again to the same background and free parameters; the results file's background parts (P49.2) adding up to its
//!   background.
//! - `bac_with_b1_b4_written_as_held` (the app's request 48, L374): a `.bac` with b1–b4 written as its points are held,
//!   not moved by b1–b4, and left alone when it already holds them.
//! - `paths_relative_across_the_root` (Publish's P61, L374): a project sharing only "/" with the run gets relative
//!   paths, `_drs_bkg_file` among them.

use drs::jacobian::Field;
use drs::lsq::{refine, Options};
use drs::powder::Model;
use drs::refinable::Refinement;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

fn dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs/background_terms")
}

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display())).iter().map(|&b| b as char).collect()
}

/// The rows of a `.prf`: (2θ, Backg, the Backg field as printed).
fn prf_background(text: &str) -> Vec<(f64, f64, String)> {
    let n: usize = text.lines().nth(1).unwrap().split_whitespace().nth(1).unwrap().parse().unwrap();
    let start = text.lines().position(|l| l.trim_start().starts_with("2Theta")).unwrap() + 1;
    text.lines()
        .skip(start)
        .take(n)
        .map(|l| {
            let f: Vec<&str> = l.split('\t').map(str::trim).collect();
            (f[0].parse().unwrap(), f[4].parse().unwrap(), f[4].to_string())
        })
        .collect()
}

/// The interval of offsets c that put ours within the print of FullProf's at every point (`background_points`, L081),
/// each point's bound widened by `extra`, FullProf's single-precision rounding of the function itself
/// (`single_precision`).
fn offset_interval(ours: &[f64], prf: &[&(f64, f64, String)], extra: &[f64]) -> (f64, f64, f64) {
    let decimals = prf[0].2.split('.').nth(1).map_or(0, str::len) as i32;
    let h = 0.5 * 10f64.powi(-decimals);
    let (mut lo, mut hi) = (f64::NEG_INFINITY, f64::INFINITY);
    let mut spread = (f64::INFINITY, f64::NEG_INFINITY);
    for ((o, r), e) in ours.iter().zip(prf).zip(extra) {
        let d = o - r.1;
        let u = (o.abs() + r.1.abs()) * 2f64.powi(-24) + e;
        lo = lo.max(d - h - u);
        hi = hi.min(d + h + u);
        spread = (spread.0.min(d), spread.1.max(d));
    }
    (lo, hi, spread.1 - spread.0)
}

/// A bound on the error of FullProf's single-precision evaluation (ε = 2⁻²⁴) of the background at 2θ = t, beyond the
/// rounding of its value that `offset_interval` counts: the polynomial Σ B_m x^m, x = t/Bkpos − 1, from t rounded
/// (ε|t|), divided (ε|t|/Bkpos) and less 1 (ε|x|), so δx ≤ ε(2|t|/Bkpos + |x|), carried into each term as
/// m|B_m||x|^(m−1)δx; each power by m − 1 products and its coefficient's (m + 1 roundings, ε(m + 1)|B_m x^m|); the sum's
/// roundings, ε Σ|B_m x^m| per term added. The inverse powers likewise (j + 2 roundings each); the Debye-like terms
/// from Q r_j rounded (ε·Q r_j·|A_j|, sin x/x changing by at most |x| times x's relative error) and their quotient (2ε|A_j|).
/// For y2o3_poly12_b9 at 143.35° (B9·x⁸ = 9 951) it is 0.024 counts, against the print's 0.05; FullProf's own sequence
/// is not known, the bound being the worst case of any order of the same operations (L353).
fn single_precision(m: &Model, t: f64) -> f64 {
    let eps = 2f64.powi(-24);
    let mut e = 0.0;
    if let drs::powder::BackgroundModel::Polynomial { origin, coefficients } = &m.background {
        if *origin != 0.0 {
            let x = t / origin - 1.0;
            let dx = eps * (2.0 * t.abs() / origin.abs() + x.abs());
            let n = coefficients.len() as f64;
            for (k, b) in coefficients.iter().enumerate() {
                let term = (b * x.powi(k as i32)).abs();
                e += eps * (k as f64 + 1.0 + n) * term;
                if k > 0 {
                    e += k as f64 * b.abs() * x.abs().powi(k as i32 - 1) * dx;
                }
            }
        }
    }
    let bt = &m.background_terms;
    if t > 0.0 {
        e += bt.inverse.iter().enumerate().map(|(j, c)| eps * (j as f64 + 3.0) * (c * t.powi(-(j as i32 + 1))).abs()).sum::<f64>();
    }
    let q = m.background_q(t);
    e += bt.debye.iter().map(|d| eps * d.amplitude.abs() * (q * d.distance + 2.0)).sum::<f64>();
    e
}

/// The run's model as FullProf read it: its `.pcr.in` for the NCY = 0 variants (FullProf rewrote their `.pcr`), its
/// `.pcr` for nac_3t2_calc (pass 2 of the refinement, `run_fullprof.py`, L035).
fn load(run: &str) -> drs::powder::Loaded {
    let d = dir();
    let (file, stem) = match run.strip_suffix("_calc") {
        Some(base) => (d.join(format!("{base}.pcr")), base),
        None => (d.join(format!("{run}.pcr.in")), run),
    };
    let p = drs::pcr::read(&latin1(&file)).unwrap_or_else(|e| panic!("{run}: {e}"));
    Model::from_pcr_fullprof(&p, &d, stem).unwrap_or_else(|e| panic!("{run}: {e}"))
}

/// (run, points compared, points left out and frozen, the largest |ours − FullProf's| there): nac_3t2's first two
/// points, in the excluded region 0–5° (no χ² counts them). At 2θ = 0 (Q = 0) FullProf's background lacks the
/// Debye-like terms (31.033 = −Σ A_j of the refined terms: it leaves sin(Qr)/(Qr) out at Q = 0, where its limit is 1;
/// the core keeps the limit). At 0.05° (Q = 0.0045 Å⁻¹, sin(Qr)/(Qr) within 10⁻⁴ of 1) FullProf's is 0.557 below the
/// core's: not explained, frozen (L353). In nac_3t2_debye the background at both is 0 under either rule (Σ A_j < 0).
const RUNS: [(&str, usize, usize, f64); 14] = [
    ("nac_3t2_calc", 2508, 2, 31.033),
    ("nac_3t2_debye", 2508, 2, 0.0),
    ("y2o3_poly12_b1", 13091, 0, 0.0),
    ("y2o3_poly12_b2", 13091, 0, 0.0),
    ("y2o3_poly12_b7", 13091, 0, 0.0),
    ("y2o3_poly12_b9", 13091, 0, 0.0),
    ("y2o3_poly12_b10", 13091, 0, 0.0),
    ("y2o3_poly12_b11", 13091, 0, 0.0),
    ("y2o3_poly12_b12", 13091, 0, 0.0),
    ("y2o3_bac_b1b2", 13091, 0, 0.0),
    ("y2o3_bac_b3b4", 13091, 0, 0.0),
    ("y2o3_more2", 13091, 0, 0.0),
    ("y2o3_more3", 13091, 0, 0.0),
    ("y2o3_more2_two", 13091, 0, 0.0),
];

#[test]
fn against_fullprof_backgrounds() {
    let mut fails = Vec::new();
    for &(run, compared, frozen, worst) in &RUNS {
        let l = load(run);
        let rows = prf_background(&latin1(&dir().join(format!("{run}.prf"))));
        let x: Vec<f64> = rows.iter().map(|r| r.0).collect();
        let ours = l.model.calculate(&x).background;
        // The points left out: nac_3t2's first two, 2θ = 0 and 0.05° (`RUNS`).
        let keep: Vec<bool> = rows.iter().map(|r| !run.starts_with("nac_3t2") || r.0 > 0.075).collect();
        let (o, p): (Vec<f64>, Vec<&(f64, f64, String)>) = ours.iter().zip(&rows).zip(&keep).filter(|x| *x.1).map(|((o, r), _)| (*o, r)).unzip();
        let extra: Vec<f64> = p.iter().map(|r| single_precision(&l.model, r.0)).collect();
        let (lo, hi, spread) = offset_interval(&o, &p, &extra);
        let least = o.iter().copied().fold(f64::INFINITY, f64::min);
        let quarter = least / 4.0;
        let delta = least.abs() * 2f64.powi(-24) / 4.0 + quarter.abs() * 2f64.powi(-24);
        let outside = (lo - quarter).max(quarter - hi).max(0.0);
        // The points left out: ours against FullProf's with the quarter added back.
        let off: Vec<f64> = ours.iter().zip(&rows).zip(&keep).filter(|x| !*x.1).map(|((o, r), _)| (o - r.1 - quarter).abs()).collect();
        let off_max = off.iter().copied().fold(0.0, f64::max);
        println!(
            "{run}: {} points compared, ours − FullProf's spans {spread:.5}, one offset in [{lo:.6}, {hi:.6}], a quarter of the least {quarter:.6} (outside by {outside:.1e}); {} left out, up to {off_max:.3} from FullProf's",
            o.len(),
            off.len()
        );
        if lo > hi {
            fails.push(format!("{run}: no offset puts the background within FullProf's print (spread {spread:.4})"));
        }
        if outside > delta {
            fails.push(format!("{run}: the offset is not a quarter of the least background ({quarter} outside [{lo}, {hi}] by {outside:.2e})"));
        }
        if o.len() != compared || off.len() != frozen || off_max > worst + 5e-4 {
            fails.push(format!("{run}: {} compared and {} left out up to {off_max}, frozen at {compared}, {frozen} and {worst}", o.len(), off.len()));
        }
    }
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}

#[test]
fn fourier_filter_refused() {
    let p = drs::pcr::read(&latin1(&dir().join("y2o3_filter.pcr.in"))).unwrap();
    let e = Model::from_pcr(&p, &dir(), "y2o3_filter").map(|_| ()).expect_err("Nba = −2 is refused");
    assert!(e.contains("Nba = -2") && e.contains("Fourier filtering"), "{e}");
}


/// FullProf's refined parameters in its last cycle: name → (value, esd) (as `common::fullprof_params`).
fn fullprof_params(out: &str) -> BTreeMap<String, (f64, f64)> {
    let mut m = BTreeMap::new();
    for l in out.lines().filter(|l| l.contains("Parameter number") && l.contains("+/-")) {
        let Some((_, r)) = l.split_once(':') else { continue };
        let mut t = r.split_whitespace();
        let name = t.next().unwrap().to_string();
        let rest = t.collect::<Vec<_>>().join(" ").replace(['(', ')'], " ").replace("+/-", " ");
        let v: Vec<f64> = rest.split_whitespace().filter_map(|s| s.parse().ok()).collect();
        if v.len() >= 2 {
            m.insert(name, (v[0], v[1]));
        }
    }
    m
}

/// FullProf's residual S and N − P + C over all non-excluded points, in its first (`first`) or last statistics block.
fn fullprof_residual(out: &str, first: bool) -> (f64, usize) {
    let lines: Vec<&str> = out.lines().collect();
    let heads: Vec<usize> = lines.iter().enumerate().filter(|(_, l)| l.contains("RELIABILITY FACTORS WITH ALL NON-EXCLUDED POINTS")).map(|(i, _)| i).collect();
    let h = if first { heads[0] } else { *heads.last().unwrap() };
    let dof = lines[h..].iter().find_map(|l| l.split_once("N-P+C:")).unwrap().1.trim().parse().unwrap();
    let k = h + lines[h..].iter().position(|l| l.contains("SumYdif")).unwrap();
    (lines[k + 1].split_whitespace().nth(4).unwrap().parse().unwrap(), dof)
}

#[test]
fn debye_refinement() {
    // nac_3t2 (3T2 neutrons, Na₂Ca₃Al₂F₁₄ and a second phase): FullProf refined 51 parameters from nac_3t2.pcr (it left
    // the file as read and wrote its values to nac_3t2.new, rounded to its print) and converged in 8 cycles. Its 51st,
    // β23 of Ca1, it freed itself (the .pcr's codeword is 0; −1.6(1.2)·10⁻⁴ at the end); the core refines the
    // codewords' 50, β23 kept at 0 (with FullProf's value the residual below is 7 564.6, further from FullProf's). The
    // engine refines the same start without damping to its own convergence.
    let d = dir();
    let p = drs::pcr::read(&latin1(&d.join("nac_3t2.pcr"))).unwrap();
    let l = Model::from_pcr(&p, &d, "nac_3t2").unwrap();
    let (mut m, data) = (l.model, l.data.unwrap());
    let pm = m.clone();
    let mut rf = Refinement::from_pcr(&p, &m).unwrap();
    rf.steering.class_damping = [1.0; 4];
    let set = rf.build(&m).unwrap();
    assert_eq!(set.free.len(), 50);
    // At the start: the residual against FullProf's at the same values (nac_3t2_calc.out, pass 2), within the bound
    // of fullprof_profiles for χ² (5·10⁻⁴ + 1.5·10⁻³ χ²).
    // The residuals are FullProf's, so under FullProf's rules (the start's model imported so; the default model's
    // tapered peak windows add 0.07 % to S).
    let pf = Model::from_pcr_fullprof(&p, &d, "nac_3t2").unwrap().model;
    let s0 = drs::lsq::statistics(&pf, &data, &pf.calculate(&data.x).ycalc, set.free.len());
    let (fp_s0, fp_dof) = fullprof_residual(&latin1(&d.join("nac_3t2_calc.out")), true);
    let our_s0 = s0.chi2 * (s0.n + s0.c - s0.p) as f64;
    println!("start: S {our_s0:.1} against FullProf's {fp_s0} ({:+.3} %), N − P + C {} against its {fp_dof}", 100.0 * (our_s0 / fp_s0 - 1.0), s0.n + s0.c - s0.p);
    assert_eq!(s0.n + s0.c - s0.p, fp_dof + 1, "one parameter less than FullProf's 51");
    // +0.199 %: the 157 points just outside a peak window's edge of a phase with X ≠ 0, where FullProf's window is
    // not the core's (L035); `fullprof_profiles` takes them from FullProf's .prf and then has χ² 3.6238 against
    // FullProf's 3.619, within its bound. Frozen at 0.2 %.
    assert!((our_s0 / fp_s0 - 1.0).abs() <= 0.002, "S at the start, frozen at 0.2 %");
    let r = refine(&mut m, &data, &set, Options { max_cycles: 40, eps: 0.1 }, |_| {});
    assert!(r.converged, "{}", r.stop);
    println!("χ² {:.4} → {:.4} in {} cycles", s0.chi2, r.stats.chi2, r.cycles);
    // Each refined value against FullProf's last, in FullProf's esd's: the background's (the polynomial's six and the
    // five Debye-like amplitudes) within 0.11, the others within 0.37; the esd's 0.927–0.990 of FullProf's (the others 0.955–0.992). The cell's a
    // of phase 1 ends 2.97 of FullProf's esd's from it, 1.50 of the core's (frozen at 1.51): FullProf prints the cell's esd's at half of
    // what its fit supports (L053, L113), its ratio here 1.98.
    let fp = fullprof_params(&latin1(&d.join("nac_3t2.out")));
    let (mut background, mut others) = (0, 0);
    for (name, (v, e)) in &fp {
        let Some(f) = drs::output_names::fullprof_field(name, &pm) else { continue };
        let Some(i) = set.free.iter().position(|q| q.fields[0].0 == f) else { continue };
        let q = &r.params[i];
        let (dev, ratio) = ((q.value - v) / e, q.esd / e);
        println!("{name:<20} {:.6} ({:.6}) against FullProf's {v} ({e}): {dev:+.2} of its esd, esd ratio {ratio:.3}", q.value, q.esd);
        match f {
            Field::Background(_) | Field::BackgroundTerm(_) => {
                background += 1;
                assert!(dev.abs() <= 0.11 && (0.92..=0.995).contains(&ratio), "{name}: {dev:.2} esd, ratio {ratio:.3}");
            }
            Field::Cell(..) => assert!((q.value - v).abs() <= 1.51 * q.esd && (1.97..=1.99).contains(&ratio), "{name}: {dev:.2} esd, ratio {ratio:.3}"),
            _ => {
                others += 1;
                assert!(dev.abs() <= 0.37 && (0.95..=0.995).contains(&ratio), "{name}: {dev:.2} esd, ratio {ratio:.3}");
            }
        }
    }
    assert_eq!((background, others), (11, 37));
    // At FullProf's last values (in full, from its .out; each moving its parameter's fields): the residual 0.64 % below
    // FullProf's last printed one, 7 624, which its pass-1 .prf gives too (7 623.4 from its Ycal with the core's
    // weights). Believed, not checked: the .prf and that residual belong to the parameters of the cycle before the last
    // (as every pass-1 .prf, L035), not to the values the .out prints after it; the core's background at the printed
    // values is the .prf's to 0.005 at 10° but 1.01 below it at 124.85°, where the polynomial's high powers dominate,
    // as a change of its coefficients would make it. Frozen at 0.65 % (L353).
    let mut mf = pf.clone();
    for (name, (v, _)) in &fp {
        let Some(f) = drs::output_names::fullprof_field(name, &pm) else { continue };
        let Some(q) = set.free.iter().find(|q| q.fields[0].0 == f) else { continue };
        let delta = (v - pf.get(f)) / q.fields[0].1;
        for &(g, c) in &q.fields {
            mf.set(g, pf.get(g) + delta * c);
        }
    }
    let s1 = drs::lsq::statistics(&mf, &data, &mf.calculate(&data.x).ycalc, set.free.len());
    let (fp_s1, _) = fullprof_residual(&latin1(&d.join("nac_3t2.out")), false);
    let our_s1 = s1.chi2 * (s1.n + s1.c - s1.p) as f64;
    let gap = our_s1 / fp_s1 - 1.0;
    println!("at FullProf's last values: S {our_s1:.1} against FullProf's {fp_s1} ({:+.3} %); the engine's minimum {:.1}", 100.0 * gap, r.stats.chi2 * (s1.n + s1.c - s1.p) as f64);
    assert!(gap.abs() <= 0.0065, "the residual at FullProf's last values, frozen at 0.65 %: {:+.3} %", 100.0 * gap);
}


#[test]
fn control_file_and_export() {
    let d = dir();
    let tmp = Path::new(env!("CARGO_MANIFEST_DIR")).join("target/test-background-terms");
    let _ = std::fs::remove_dir_all(&tmp);
    std::fs::create_dir_all(&tmp).unwrap();
    // (run, its file, Nba and MORE of the export)
    for (run, file, nba, more) in [
        ("nac_3t2", "nac_3t2.pcr", -1, 0),
        ("y2o3_poly12_b10", "y2o3_poly12_b10.pcr.in", -3, 0),
        ("y2o3_more2_two", "y2o3_more2_two.pcr.in", 0, 2),
        ("y2o3_more3", "y2o3_more3.pcr.in", 0, 2),
        ("y2o3_bac_b3b4", "y2o3_bac_b3b4.pcr.in", 1, 0),
    ] {
        let p = drs::pcr::read(&latin1(&d.join(file))).unwrap();
        let l = Model::from_pcr(&p, &d, run).unwrap();
        let data = l.data.unwrap();
        // The control file in another folder: the data and the profiles named relative to it.
        let c = drs::control::Control::from_pcr(&p, &d, run, &tmp).unwrap();
        let text = c.to_text();
        let back = drs::control::Control::from_text_in(&text, &tmp).unwrap_or_else(|e| panic!("{run}: {e:?}"));
        assert_eq!(back.model.to_text(), c.model.to_text(), "{run}: the model read back");
        assert_eq!(back.model.calculate(&data.x).ycalc, l.model.calculate(&data.x).ycalc, "{run}: the pattern read back");
        assert_eq!(back.refinement, c.refinement, "{run}: the free parameters read back");
        // Every further term where the file writes it.
        let n = c.model.background_terms.len();
        let located: Vec<_> = drs::control::locate::all(&text, Some(&tmp)).into_iter().filter(|(a, ..)| a.contains("background_debye") || a.contains("background_inverse") || a.contains("background_profile") || a.contains("background_file_")).collect();
        assert_eq!(located.len(), n, "{run}: {n} further terms located");
        for (a, _, r) in &located {
            let f = drs::refinable::field_of(a, &c.model).unwrap();
            assert!(matches!(f, Field::BackgroundTerm(_)), "{a}");
            assert_eq!(text[r.clone()].parse::<f64>().unwrap(), c.model.get(f), "{a}");
        }
        // The .pcr export, read again from its folder: the same background at the data's points (a profile's .bac
        // written at them), and the same free parameters.
        let e = drs::control::pcr_export::export(&c, run, Some(&data)).unwrap();
        std::fs::write(tmp.join(format!("{run}.pcr")), &e.pcr).unwrap();
        for (name, t) in &e.files {
            std::fs::write(tmp.join(name), t).unwrap();
        }
        let p2 = drs::pcr::read(&e.pcr).unwrap_or_else(|er| panic!("{run}: {er}\n{}", e.pcr));
        assert_eq!((p2.flag("Nba"), !p2.profiles.is_empty()), (nba, more == 2), "{run}: Nba and MORE of the export");
        let m2 = Model::from_pcr(&p2, &tmp, run).unwrap_or_else(|er| panic!("{run}: {er}")).model;
        let (b1, b2) = (c.model.background_at(&data.x), m2.background_at(&data.x));
        let worst = b1.iter().zip(&b2).map(|(a, b)| (a - b).abs() / a.abs().max(1.0)).fold(0.0, f64::max);
        assert!(worst < 1e-9, "{run}: the background through the .pcr, {worst:.2e}");
        let free = |rf: &Refinement, m: &Model| rf.build(m).unwrap().free.iter().map(|f| f.fields[0].0).collect::<Vec<_>>();
        assert_eq!(free(&Refinement::from_pcr(&p2, &m2).unwrap(), &m2), free(&c.refinement, &c.model), "{run}: the free parameters through the .pcr");
        // The results file: the background's parts (P49.2) add up to it.
        let (cif, _) = drs::results::write(&c, &data, None, None, run).unwrap();
        let rs = drs::results::read(&cif).unwrap();
        let names: Vec<&str> = rs.background_parts.iter().map(|p| p.0.as_str()).collect();
        let bt = &c.model.background_terms;
        let mut expected = vec!["smooth".to_string()];
        expected.extend((1..=bt.debye.len()).map(|j| format!("debye_{j}")));
        expected.extend((1..=bt.profiles.len()).map(|j| format!("file_{j}")));
        if bt.debye.is_empty() && bt.profiles.is_empty() {
            expected.clear();
        }
        assert_eq!(names, expected, "{run}: the background's parts");
        for i in (0..data.x.len()).step_by(97) {
            if !rs.background_parts.is_empty() {
                let sum: f64 = rs.background_parts.iter().map(|p| p.1[i]).sum();
                assert!((sum - rs.background[i]).abs() <= 1e-6 * rs.background[i].abs().max(1.0), "{run}: the parts at {}", data.x[i]);
            }
        }
        println!("{run}: control file {} kB, {n} further terms located, .pcr Nba {nba} MORE {more}, {} parts in the results", text.len() / 1024, names.len());
    }
}


/// Nba = 1's b1–b4 as values (the app's request 47, Publish's P57): the `.bac`'s points kept as read, b1–b4 four values
/// of the background (`background_terms::FileScale`, `Field::BackgroundTerm`), their derivatives against differences,
/// refined back from moved values, written to the `.pcr` with their codewords and read again, and in the results file
/// (`_drs_bkg_file*`, the function's text).
#[test]
fn file_scale_values() {
    let d = dir();
    let tmp = Path::new(env!("CARGO_MANIFEST_DIR")).join("target/test-background-file");
    let _ = std::fs::remove_dir_all(&tmp);
    std::fs::create_dir_all(&tmp).unwrap();
    for (run, values) in [("y2o3_bac_b1b2", [0.8, 5.0, 0.0, 0.0]), ("y2o3_bac_b3b4", [1.0, 0.0, 0.01, 0.5])] {
        let p = drs::pcr::read(&latin1(&d.join(format!("{run}.pcr.in")))).unwrap();
        let l = Model::from_pcr(&p, &d, run).unwrap();
        let (m, data) = (l.model, l.data.unwrap());
        let f = m.background_terms.file.as_ref().unwrap();
        assert_eq!((f.file.clone(), f.values()), (format!("{run}.bac"), values), "{run}");
        // The points as the file holds them.
        let file = drs::background::read_path(&d.join(format!("{run}.bac"))).unwrap();
        assert!(matches!(&m.background, drs::powder::BackgroundModel::Points(q) if *q == file.points), "{run}: the file's points as read");
        let first = m.background_terms.len() - 4;
        let fields: Vec<Field> = (0..4).map(|k| Field::BackgroundTerm(first + k)).collect();
        let names: Vec<String> = fields.iter().map(|&f| drs::refinable::address(f, &m)).collect();
        assert_eq!(names, ["pattern1.background_file_scale.1", "pattern1.background_file_offset.1", "pattern1.background_file_stretch.1", "pattern1.background_file_shift.1"]);
        // The derivatives of the pattern against central differences.
        let (_, jac) = m.derivatives(&data.x, &fields);
        for (j, &f) in fields.iter().enumerate() {
            let h = 1e-6;
            let (mut a, mut b) = (m.clone(), m.clone());
            a.set(f, m.get(f) + h);
            b.set(f, m.get(f) - h);
            let (ya, yb) = (a.calculate(&data.x).ycalc, b.calculate(&data.x).ycalc);
            let top = jac[j].iter().fold(0.0f64, |s, v| s.max(v.abs()));
            // At a .bac point the slope changes: a difference across it averages the two sides, as the derivative does.
            let worst = (0..data.x.len()).map(|i| ((ya[i] - yb[i]) / (2.0 * h) - jac[j][i]).abs()).fold(0.0, f64::max);
            assert!(top > 0.0 && worst < 1e-4 * top, "{run} {f:?}: off by {worst:.3e} of {top:.3e}");
        }
        // Refined back from moved values on the core's own pattern.
        let y = m.calculate(&data.x).ycalc;
        let synthetic = drs::data::Pattern { y: y.clone(), sigma: y.iter().map(|v| v.max(1.0).sqrt()).collect(), ..data.clone() };
        let mut moved = m.clone();
        for (k, dv) in [0.05, 3.0, 0.002, 0.05].into_iter().enumerate() {
            moved.set(fields[k], m.get(fields[k]) + dv);
        }
        let mut rf = Refinement { parameters: fields.iter().map(|&f| vec![(f, 1.0)]).collect(), ..Default::default() };
        rf.steering.guards.jump = f64::INFINITY;
        let set = rf.build(&moved).unwrap();
        let r = refine(&mut moved, &synthetic, &set, Options { max_cycles: 30, eps: 1e-6 }, |_| {});
        let got: Vec<f64> = fields.iter().map(|&f| moved.get(f)).collect();
        println!("{run}: b1–b4 back to {got:?} in {} cycles (χ² {:.2e})", r.cycles, r.stats.chi2);
        for (g, w) in got.iter().zip(values) {
            assert!((g - w).abs() < 1e-5 * w.abs().max(1.0), "{run}: {got:?}");
        }
        // The control file and the .pcr export with b1–b4 free: their codewords, the same pattern read again.
        let mut c = drs::control::Control::from_pcr(&p, &d, run, &tmp).unwrap();
        c.refinement.parameters.extend(fields.iter().map(|&f| vec![(f, 1.0)]));
        let text = c.to_text();
        assert!(text.contains("background_file_scale") && text.contains("file = "), "{run}");
        let back = drs::control::Control::from_text_in(&text, &tmp).unwrap_or_else(|e| panic!("{run}: {e:?}"));
        assert_eq!(back.model.background_terms.file.as_ref().map(|f| f.values()), Some(values), "{run}: the control file");
        let e = drs::control::pcr_export::export(&c, run, Some(&data)).unwrap();
        for (name, t) in &e.files {
            std::fs::write(tmp.join(name), t).unwrap();
        }
        let p2 = drs::pcr::read(&e.pcr).unwrap_or_else(|er| panic!("{run}: {er}\n{}", e.pcr));
        assert_eq!(p2.flag("Nba"), 1, "{run}");
        let m2 = Model::from_pcr(&p2, &tmp, run).unwrap().model;
        assert_eq!(m2.background_terms.file.as_ref().map(|f| f.values()), Some(values), "{run}: through the .pcr");
        let (b1, b2) = (c.model.background_at(&data.x), m2.background_at(&data.x));
        let worst = b1.iter().zip(&b2).map(|(a, b)| (a - b).abs() / a.abs().max(1.0)).fold(0.0, f64::max);
        assert!(worst < 1e-6, "{run}: the background through the .pcr, {worst:.2e}");
        let free = |rf: &Refinement, m: &Model| rf.build(m).unwrap().free.iter().map(|f| f.fields[0].0).collect::<Vec<_>>();
        let f2 = free(&Refinement::from_pcr(&p2, &m2).unwrap(), &m2);
        assert!(fields.iter().all(|f| f2.contains(f)), "{run}: b1–b4 free through the .pcr: {f2:?}");
        // The results file (P57): the file and b1–b4 with su's, and the function's text.
        let mut cr = c.clone();
        let set = c.refinement.build(&cr.model).unwrap();
        let r = refine(&mut cr.model, &data, &set, Options { max_cycles: 1, eps: 1e-3 }, |_| {});
        let (cif, _) = drs::results::write(&cr, &data, Some(&set), Some(&r), run).unwrap();
        assert!(cif.lines().any(|l| l.starts_with("_drs_bkg_file ") && l.contains(&format!("{run}.bac"))), "{run}");
        for item in ["scale", "offset", "stretch", "shift"] {
            let l = cif.lines().find(|l| l.starts_with(&format!("_drs_bkg_file_{item}"))).unwrap_or_else(|| panic!("{run}: {item}"));
            assert!(l.contains('('), "{run}: {l}");
        }
        assert!(cif.contains(&format!("read from {run}.bac, as b~1~ B[(1 + b~3~)2\\q + b~4~] + b~2~")), "{run}");
    }
}

/// The app's request 48: a save of a control file whose background is a `.bac` with b1–b4 writes the file's points as
/// held, not the background with b1–b4 at the data's points (which the next reading would move again), and leaves a
/// file that holds those points as it is: FullProf's own y2o3_bac_b3b4.bac is not rewritten by a save that changed
/// nothing in it.
#[test]
fn bac_with_b1_b4_written_as_held() {
    let d = dir();
    let tmp = Path::new(env!("CARGO_MANIFEST_DIR")).join("target/test-background-held");
    let _ = std::fs::remove_dir_all(&tmp);
    std::fs::create_dir_all(&tmp).unwrap();
    for run in ["y2o3_bac_b1b2", "y2o3_bac_b3b4"] {
        let p = drs::pcr::read(&latin1(&d.join(format!("{run}.pcr.in")))).unwrap();
        let c = drs::control::Control::from_pcr(&p, &d, run, &d).unwrap();
        let data = Model::from_pcr(&p, &d, run).unwrap().data.unwrap();
        let file = drs::background::read_path(&d.join(format!("{run}.bac"))).unwrap();
        // Given the data's points or not, the text is the points as read (the first 19.1834, not b1–b4's 19.0873).
        for x in [Some(&data.x[..]), None] {
            let (name, text) = c.background_file_text(x).unwrap().unwrap();
            assert_eq!(name, format!("{run}.bac"));
            assert_eq!(drs::background::read(&text).unwrap().points, file.points, "{run}: the points as held");
        }
        // FullProf's file, unchanged, is left alone.
        assert_eq!(c.background_file_to_write(&d, Some(&data.x)).unwrap().unwrap(), None, "{run}: nothing to write");
        // A copy saved and read again keeps its points and its pattern, save after save.
        let mut copy = c.clone();
        let (name, text) = copy.background_file_text(Some(&data.x)).unwrap().unwrap();
        std::fs::write(tmp.join(&name), &text).unwrap();
        let back = drs::control::Control::from_text_in(&c.to_text().replace(&format!("{run}.dat"), &d.join(format!("{run}.dat")).to_string_lossy()), &tmp).unwrap();
        assert_eq!(back.model.background, c.model.background, "{run}: read again");
        assert_eq!(back.model.calculate(&data.x).ycalc, c.model.calculate(&data.x).ycalc, "{run}: the pattern read again");
        assert_eq!(back.background_file_to_write(&tmp, Some(&data.x)).unwrap().unwrap(), None, "{run}: a second save writes nothing");
        // A point edited is written.
        if let drs::powder::BackgroundModel::Points(q) = &mut copy.model.background {
            q[3].1 += 1.0;
        }
        let w = copy.background_file_to_write(&tmp, Some(&data.x)).unwrap().unwrap().expect("an edited point is written");
        assert_eq!(drs::background::read(&w.1).unwrap().points[3].1, file.points[3].1 + 1.0, "{run}");
    }
}

/// Publish's P61: a results file written for a project that shares only the root with the run (a project in the
/// system's temporary folder, the run in the repository) names the data and the background file relative to the
/// control file, as `relocate` writes them, not by full paths: `_drs_bkg_file` is the control file's own `file`.
#[test]
fn paths_relative_across_the_root() {
    let d = dir();
    let base = std::env::temp_dir().join("drs-p61-project");
    let _ = std::fs::remove_dir_all(&base);
    std::fs::create_dir_all(&base).unwrap();
    let run = "y2o3_bac_b3b4";
    let p = drs::pcr::read(&latin1(&d.join(format!("{run}.pcr.in")))).unwrap();
    let c = drs::control::Control::from_pcr(&p, &d, run, &base).unwrap();
    let data = Model::from_pcr(&p, &d, run).unwrap().data.unwrap();
    let bac = c.background_file.clone().unwrap();
    let dat = c.data.as_ref().unwrap().file.clone();
    for f in [&bac, &dat] {
        assert!(f.starts_with("../"), "{f}: relative");
        let there = std::fs::canonicalize(base.join(f)).unwrap();
        assert_eq!(there.parent().unwrap(), std::fs::canonicalize(&d).unwrap(), "{f}: the run's file");
    }
    assert_eq!(c.model.background_terms.file.as_ref().unwrap().file, bac);
    let (cif, _) = drs::results::write(&c, &data, None, None, run).unwrap();
    let item = cif.lines().find(|l| l.starts_with("_drs_bkg_file ")).unwrap();
    assert_eq!(item.split_whitespace().nth(1).unwrap().trim_matches('\''), bac, "_drs_bkg_file");
    // The control file written into the results file reads its files from the control file's folder.
    let back = drs::control::Control::from_text_in(&c.to_text(), &base).unwrap();
    assert_eq!(back.model.background, c.model.background);
    // Paths without a common root stay absolute; two relative names relative to each other.
    assert_eq!(drs::control::relative(Path::new("a/b/c/x.dat"), Path::new("a/b")), "c/x.dat");
}
