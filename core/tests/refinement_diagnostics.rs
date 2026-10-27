//! P45, refinement diagnostics (`impact`, L272): each parameter's impact on χ² and on the calculated pattern, and the
//! normal probability plot of the weighted residuals.
//!
//! - **The normal distribution.** ln Φ against its values to 10⁻¹³ (both branches, and across x = −3 where they meet);
//!   Φ⁻¹ round-trips Φ; the expected normal order statistics are exact for n = 2 and 3 (±1/√π, ±3/(2√π), to 10⁻⁹),
//!   and for n = 10 as mpmath's quadrature at 30 digits gives them (`tests/data/expected_normal_order_statistics.py`),
//!   to 10⁻⁹.
//! - **The plot** on Y₂O₃ (FullProf's tutorial, its converged values): its own data, against GSAS-II's expected
//!   values (`PlotDeltSig`), and on data made from the calculated pattern with Gaussian noise of the stated σ: slope
//!   1 and intercept 0 within their statistical spread; with σ stated as half the noise's, slope 2. It bites.
//! - **The impact.** The analytic slope ∂χ²/∂p against GSAS-II's finite-difference slope of the rms weighted residual
//!   (`AllPrmDerivs`: steps of max(10⁻⁴|p|, 10⁻⁶) down and up), for every parameter; at the engine's convergence every
//!   refined parameter's first-order shift is a small fraction of its esd; a zero shift of 0.01° left unrefined ranks
//!   first, with its shift and the fall of χ² as refining it alone gives them; for a linear parameter (the scale) the
//!   predicted fall of χ² is the engine's to 10⁻⁶.

use drs::impact::{expected_normal_order_statistics, impact_of, impacts, ln_normal_cdf, normal_probability, normal_quantile};
use drs::jacobian::Field;
use drs::powder::Model;
use std::path::Path;

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap().iter().map(|&b| b as char).collect()
}

#[test]
fn normal_distribution() {
    // Φ(x) to 16 digits (from the error function's tables; Φ(−10) = ½ erfc(10/√2)).
    for (x, v) in [(-10.0, 7.619853024160527e-24), (-5.0, 2.866515718791939e-7), (-3.0, 1.3498980316300946e-3), (-1.0, 0.15865525393145707), (0.0, 0.5), (1.96, 0.9750021048517795), (3.5, 0.9997673709209645)] {
        let got = ln_normal_cdf(x).exp();
        assert!((got - v).abs() <= 1e-13 * v, "Φ({x}) = {got:e}, not {v:e}");
    }
    // The two branches meet at −3.
    let (a, b) = (ln_normal_cdf(-3.0 - 1e-12), ln_normal_cdf(-3.0 + 1e-12));
    assert!((a - b).abs() < 1e-11, "{a} {b}");
    for p in [1e-300, 1e-100, 1e-20, 1e-6, 0.001, 0.025, 0.1, 0.3, 0.5, 0.7, 0.975, 0.999999] {
        let x = normal_quantile(p);
        assert!((ln_normal_cdf(x) - p.ln()).abs() <= 1e-12 * p.ln().abs().max(1.0), "Φ(Φ⁻¹({p}))");
    }
    assert!((normal_quantile(0.975) - 1.959963984540054).abs() < 1e-13);

    // Expected normal order statistics.
    let pi = std::f64::consts::PI;
    let e2 = expected_normal_order_statistics(2);
    assert!((e2[1] - 1.0 / pi.sqrt()).abs() < 1e-9 && (e2[0] + e2[1]).abs() < 1e-15, "{e2:?}");
    let e3 = expected_normal_order_statistics(3);
    assert!((e3[2] - 1.5 / pi.sqrt()).abs() < 1e-9 && e3[1] == 0.0, "{e3:?}");
    let e10 = expected_normal_order_statistics(10);
    for (k, v) in [(9, 1.53875273083517), (8, 1.00135704457581), (7, 0.656059105364761), (6, 0.375764696997878), (5, 0.122667752284338)] {
        assert!((e10[k] - v).abs() <= 1e-9, "n = 10, order {}: {} against {v}", k + 1, e10[k]);
    }
    // A large n: sorted, symmetric; the extremes as mpmath gives them (orders 1 and 101 of 13000), and Blom's
    // approximation Φ⁻¹((i − 3/8)/(n + 1/4)) within 10⁻³ in the middle (it is 0.4 % off at the extremes).
    let n = 13000;
    let e = expected_normal_order_statistics(n);
    assert!(e.windows(2).all(|w| w[0] < w[1]) && (e[0] + e[n - 1]).abs() < 1e-15);
    for (i, v) in [(0, -3.91565922316484), (100, -2.42118140450983)] {
        assert!((e[i] - v).abs() < 1e-9, "order {} of {n}: {} against {v}", i + 1, e[i]);
    }
    for i in [3000, 6499] {
        let blom = normal_quantile((i as f64 + 1.0 - 0.375) / (n as f64 + 0.25));
        assert!((e[i] - blom).abs() < 1e-3 * blom.abs().max(1.0), "order {}: {} against Blom's {blom}", i + 1, e[i]);
    }
}

fn y2o3() -> (drs::pcr::Pcr, drs::powder::Loaded) {
    let d = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs/Y2O3");
    let p = drs::pcr::read(&latin1(&d.join("y2o3_structure_conv.pcr"))).unwrap();
    let l = Model::from_pcr_with_data(&p, &d, "y2o3_structure", None).unwrap();
    (p, l)
}

/// GSAS-II's expected values (PlotDeltSig): Abramowitz & Stegun 26.2.23 at equally spaced probabilities 0…1.
fn gsas2_expected(n: usize) -> Vec<f64> {
    (0..n)
        .map(|k| {
            let dx = k as f64 / (n - 1) as f64;
            let t = (1.0 / (dx * dx)).ln().sqrt();
            let top = 2.515517 + 0.802853 * t + 0.010328 * t * t;
            let bot = 1.0 + 1.432788 * t + 0.189269 * t * t + 0.001308 * t * t * t;
            if dx > 0.0 { -(t - top / bot) } else { t - top / bot }
        })
        .collect()
}

#[test]
fn normal_probability_plot() {
    let (_, l) = y2o3();
    let (m, data) = (l.model, l.data.unwrap());
    let yc = m.calculate(&data.x).ycalc;
    let np = normal_probability(&m, &data, &yc);
    // GSAS-II's line over its own expected values, for the same residuals.
    let g = gsas2_expected(np.residuals.len());
    let (gs, gi, gn) = drs::impact::central_line(&g, &np.residuals);
    println!("Y2O3: {} points; slope {:.4}, intercept {:.4} over {} central points; GSAS-II's expected values: {gs:.4}, {gi:.4} over {gn}", np.residuals.len(), np.slope, np.intercept, np.central);
    assert!((np.slope - gs).abs() <= 0.02 * gs && (np.intercept - gi).abs() <= 0.02 * gs.max(1.0), "the core's and GSAS-II's lines differ by more than 2 %");
    assert!((np.central as f64 / np.residuals.len() as f64 - 0.6827).abs() < 0.001, "the central part is the 68 % within ±1");

    // Gaussian noise of the stated σ: slope 1, intercept 0, within 4 times their spread (≈ 1/√(2N_c) and 1/√N_c).
    let mut state = 0x2545F4914F6CDD1Du64;
    let mut uniform = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    let mut noisy = data.clone();
    for ((y, s), c) in noisy.y.iter_mut().zip(noisy.sigma.iter_mut()).zip(&yc) {
        *s = c.max(1.0).sqrt();
        let (u1, u2) = (uniform().max(1e-300), uniform());
        let z = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
        *y = c + *s * z;
    }
    let np = normal_probability(&m, &noisy, &yc);
    let spread = 1.0 / (np.central as f64).sqrt();
    println!("noise of the stated σ: slope {:.4}, intercept {:.4} (spread {spread:.4})", np.slope, np.intercept);
    assert!((np.slope - 1.0).abs() < 4.0 * spread && np.intercept.abs() < 4.0 * spread);
    // σ stated as half the noise's: slope 2.
    noisy.sigma.iter_mut().for_each(|s| *s /= 2.0);
    let np2 = normal_probability(&m, &noisy, &yc);
    assert!((np2.slope - 2.0 * np.slope).abs() < 1e-9 && (np2.intercept - 2.0 * np.intercept).abs() < 1e-9);
}

#[test]
fn parameter_impact() {
    let (p, l) = y2o3();
    let (mut m, data) = (l.model, l.data.unwrap());
    let set = drs::refinable::ParameterSet::from_pcr(&p, &m).unwrap();
    // To the engine's own convergence.
    let r = drs::lsq::refine(&mut m, &data, &set, drs::lsq::Options { max_cycles: 30, eps: 1e-4 }, |_| {});
    let im = impacts(&m, &data, &set).unwrap();
    let (refined, free): (Vec<_>, Vec<_>) = im.iter().partition(|i| i.refined);
    println!("Y2O3 at the engine's convergence ({}): {} refined, {} not; χ² {:.4}", r.stop, refined.len(), free.len(), r.stats.chi2);
    for i in im.iter().take(8) {
        println!("  {:40} refined {:5}  shift {:+.4e} = {:+.3} esd  χ² −{:.3e}; largest change {:.1} at {:.2}°; {} regions", i.name, i.refined, i.shift, i.esds, i.gain, i.largest, i.at, i.regions.len());
    }
    assert_eq!(refined.len(), set.free.len());
    // A column of rounding noise moves nothing: Y₂O₃'s March–Dollase G1 at 1 (no preferred orientation, a cubic phase).
    let g1 = im.iter().find(|i| i.name == "phase1.po_g1").unwrap();
    assert!(g1.shift == 0.0 && g1.gain == 0.0 && g1.regions.is_empty(), "{g1:?}");
    let worst = refined.iter().map(|i| i.esds.abs()).fold(0.0, f64::max);
    assert!(worst < 0.05, "a refined parameter would still move {worst:.3} esd");

    // The slope against GSAS-II's: −d rms/dp by central differences at its step, rms = √(S/N) over the points used.
    let rms = |mm: &Model| {
        let yc = mm.calculate(&data.x).ycalc;
        let w = drs::lsq::weights(mm, &data, &yc);
        let (s, n) = (0..w.len()).filter(|&i| w[i] > 0.0).fold((0.0, 0usize), |(s, n), i| (s + w[i] * (data.y[i] - yc[i]).powi(2), n + 1));
        ((s / n as f64).sqrt(), n)
    };
    let (rms0, n) = rms(&m);
    let stats = drs::lsq::statistics(&m, &data, &m.calculate(&data.x).ycalc, set.free.len());
    let dof = (stats.n + stats.c - stats.p) as f64;
    let mut compared = 0;
    for i in im.iter().filter(|i| !i.refined) {
        let one = drs::refinable::ParameterSet::from_addresses(std::slice::from_ref(&i.name), &[], &[], &m).unwrap();
        let fields = one.free[0].fields.clone();
        let delta = (1e-4 * i.value.abs()).max(1e-6);
        let at = |d: f64| {
            let mut mm = m.clone();
            for (f, mult) in &fields {
                mm.set(*f, mm.get(*f) + mult * d);
            }
            rms(&mm).0
        };
        let g2 = 0.5 * (at(-delta) - at(delta)) / delta;
        let ours = -dof / (2.0 * n as f64 * rms0) * i.slope;
        let tol = 2e-3 * g2.abs() + 1e-9 * rms0 / delta.max(1e-6);
        assert!((ours - g2).abs() <= tol.max(1e-7 * rms0), "{}: slope {ours:e} against GSAS-II's {g2:e}", i.name);
        compared += 1;
    }
    println!("{compared} unrefined slopes as GSAS-II's finite differences give them");
    assert!(compared > 10);

    // A zero shift of 0.01° left unrefined: first, with its shift and the fall of χ² refining it alone gives.
    let zero = m.fields().into_iter().find(|f| f.field == Field::Zero).unwrap().address;
    let mut off = m.clone();
    off.zero += 0.01;
    let mut without = set.clone();
    without.free.retain(|f| f.fields.iter().all(|x| x.0 != Field::Zero));
    let im = impacts(&off, &data, &without).unwrap();
    let first = im.iter().find(|i| !i.refined).unwrap();
    assert_eq!(first.name, zero, "the zero first among the unrefined: {:?}", im.iter().take(3).map(|i| &i.name).collect::<Vec<_>>());
    assert!((first.shift + 0.01).abs() < 0.03 * 0.01, "shift {} for the 0.01° planted", first.shift);
    let alone = drs::refinable::ParameterSet::from_addresses(std::slice::from_ref(&zero), &[], &[], &off).unwrap();
    let chi_before = drs::lsq::statistics(&off, &data, &off.calculate(&data.x).ycalc, without.free.len()).chi2;
    let mut fit = off.clone();
    drs::lsq::refine(&mut fit, &data, &alone, drs::lsq::Options { max_cycles: 20, eps: 1e-6 }, |_| {});
    let chi_after = drs::lsq::statistics(&fit, &data, &fit.calculate(&data.x).ycalc, without.free.len()).chi2;
    println!("zero: predicted shift {:.5}°, χ² −{:.4}; refined alone: {:.5}°, χ² −{:.4}", first.shift, first.gain, fit.zero - off.zero, chi_before - chi_after);
    assert!((first.gain - (chi_before - chi_after)).abs() < 0.05 * (chi_before - chi_after));
    // Where it moves the pattern: on the peaks' flanks, the largest at a reflection.
    let one = impact_of(&off, &data, &without, &zero).unwrap();
    assert_eq!(one.change.len(), data.x.len());
    assert!(!one.regions.is_empty() && one.regions.iter().all(|r| r.1 - r.0 < 2.0), "{:?}", &one.regions[..one.regions.len().min(5)]);
    let peaks = off.calculate(&data.x).peaks;
    assert!(peaks.iter().any(|pk| (pk.two_theta - one.at).abs() < 2.0 * pk.fwhm), "largest change at {}°", one.at);

    // A linear parameter (the scale), 1 % off and left unrefined: the fall of χ² as the engine gives it.
    let scale = m.fields().into_iter().find(|f| f.field == Field::Scale(0)).unwrap().address;
    let mut off = m.clone();
    off.phases[0].scale *= 1.01;
    let mut without = set.clone();
    without.free.retain(|f| f.fields.iter().all(|x| x.0 != Field::Scale(0)));
    let i = impact_of(&off, &data, &without, &scale).unwrap();
    let alone = drs::refinable::ParameterSet::from_addresses(std::slice::from_ref(&scale), &[], &[], &off).unwrap();
    let chi_before = drs::lsq::statistics(&off, &data, &off.calculate(&data.x).ycalc, without.free.len()).chi2;
    let mut fit = off.clone();
    drs::lsq::refine(&mut fit, &data, &alone, drs::lsq::Options { max_cycles: 20, eps: 1e-8 }, |_| {});
    let chi_after = drs::lsq::statistics(&fit, &data, &fit.calculate(&data.x).ycalc, without.free.len()).chi2;
    println!("scale: predicted χ² −{:.6}, refined alone −{:.6}", i.gain, chi_before - chi_after);
    assert!((i.gain - (chi_before - chi_after)).abs() < 1e-6 * (chi_before - chi_after));
    assert!((fit.phases[0].scale - off.phases[0].scale - i.shift).abs() < 1e-6 * i.shift.abs());
}
