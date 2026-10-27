//! Single-peak fitting (P51, `peakfit`).
//!
//! - **The esd's are calibrated:** 400 synthetic pseudo-Voigt peaks over a sloping background, each with Gaussian
//!   noise of variance y (counting statistics); each parameter lies within one esd of the truth in 62–74 % of the
//!   trials (68.3 % expected; ±6 % is about 2.5 standard errors of a fraction over 400) and the mean deviation is
//!   below 0.2 esd. The same for a Kα₁/Kα₂ doublet (Cu, ratio 0.5).
//! - A Gaussian and a Lorentzian fitted with η held: η's esd is None, the others recovered.
//! - Real peaks against FullProf: see `fullprof_peaks`.

use drs::peakfit::{fit, PeakRequest, PeakShape};

/// A deterministic normal deviate (Box–Muller over a 64-bit LCG).
struct Noise(u64);
impl Noise {
    fn uniform(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    }
    fn normal(&mut self) -> f64 {
        (-2.0 * self.uniform().ln()).sqrt() * (2.0 * std::f64::consts::PI * self.uniform()).cos()
    }
}

/// (x₀, A, H, η, b₀, b₁): the truth of the synthetic peaks.
const TRUTH: [f64; 6] = [40.0, 3000.0, 0.15, 0.4, 200.0, -20.0];

fn synthetic(noise: &mut Noise, doublet: Option<(f64, f64)>) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let x: Vec<f64> = (0..161).map(|i| 39.3 + 0.01 * i as f64).collect();
    let xc = 0.5 * (x[0] + x[160]);
    let [x0, a, h, eta, b0, b1] = TRUTH;
    let ideal: Vec<f64> = x
        .iter()
        .map(|&t| {
            let mut y = a * drs::profile::pseudo_voigt(t - x0, eta, h) + b0 + b1 * (t - xc);
            if let Some((rl, ratio)) = doublet {
                let x2 = 2.0 * (rl * (x0.to_radians() / 2.0).sin()).asin().to_degrees();
                y += ratio * a * drs::profile::pseudo_voigt(t - x2, eta, h);
            }
            y
        })
        .collect();
    let y: Vec<f64> = ideal.iter().map(|&v| v + v.sqrt() * noise.normal()).collect();
    let s: Vec<f64> = y.iter().map(|v| v.max(1.0).sqrt()).collect();
    (x, y, s)
}

fn coverage(doublet: Option<(f64, f64)>, label: &str) {
    let mut noise = Noise(2026);
    let trials = 400;
    let (mut within, mut bias) = ([0usize; 4], [0.0f64; 4]);
    for _ in 0..trials {
        let (x, y, s) = synthetic(&mut noise, doublet);
        let r = fit(&x, &y, &s, &PeakRequest { from: 39.3, to: 40.9, shape: PeakShape::PseudoVoigt, doublet }).unwrap();
        assert!(r.converged, "{label}: not converged");
        for (k, e) in [r.position, r.intensity, r.fwhm, r.eta].iter().enumerate() {
            let z = (e.value - TRUTH[k]) / e.esd.unwrap();
            if z.abs() <= 1.0 {
                within[k] += 1;
            }
            bias[k] += z / trials as f64;
        }
    }
    let frac: Vec<f64> = within.iter().map(|&w| w as f64 / trials as f64).collect();
    println!("{label}: within one esd (x₀, A, H, η) {frac:.3?}; mean deviation in esd's {bias:.3?}");
    for k in 0..4 {
        assert!((0.62..=0.74).contains(&frac[k]), "{label}: parameter {k} within one esd in {:.3}", frac[k]);
        assert!(bias[k].abs() < 0.2, "{label}: parameter {k} biased by {:.3} esd", bias[k]);
    }
}

#[test]
fn esds_calibrated() {
    coverage(None, "single peak");
    coverage(Some((1.544390 / 1.540560, 0.5)), "Kα doublet");
}

#[test]
fn shapes_held() {
    let x: Vec<f64> = (0..121).map(|i| 19.4 + 0.01 * i as f64).collect();
    for (shape, eta) in [(PeakShape::Gaussian, 0.0), (PeakShape::Lorentzian, 1.0)] {
        let y: Vec<f64> = x.iter().map(|&t| 500.0 * drs::profile::pseudo_voigt(t - 20.0, eta, 0.2) + 50.0).collect();
        let s: Vec<f64> = y.iter().map(|v| v.sqrt()).collect();
        let r = fit(&x, &y, &s, &PeakRequest { from: 19.4, to: 20.6, shape, doublet: None }).unwrap();
        println!("{shape:?}: x₀ {:?} A {:?} H {:?} η {:?}", r.position, r.intensity, r.fwhm, r.eta);
        assert_eq!(r.eta.esd, None);
        assert_eq!(r.eta.value, eta);
        assert!((r.position.value - 20.0).abs() < 1e-6 && (r.intensity.value / 500.0 - 1.0).abs() < 1e-6 && (r.fwhm.value / 0.2 - 1.0).abs() < 1e-6, "{r:?}");
    }
    // Too few points, and no peak.
    assert!(fit(&x, &x, &x, &PeakRequest { from: 19.4, to: 19.45, shape: PeakShape::PseudoVoigt, doublet: None }).is_err());
    let flat = vec![100.0; x.len()];
    assert!(fit(&x, &flat, &vec![10.0; x.len()], &PeakRequest { from: 19.4, to: 20.6, shape: PeakShape::PseudoVoigt, doublet: None }).is_err());
}

/// FullProf's last reflection table of a `.out`: (2θ, Hw, η, ΣIcalc) per position (reflections at one 2θ summed).
fn fullprof_positions(out: &str) -> Vec<(f64, f64, f64, f64)> {
    let lines: Vec<&str> = out.lines().collect();
    // The table comes in pages, each under its header: from the last page that starts with reflection 1.
    let first_row = |i: usize| lines[i + 1..].iter().find_map(|l| l.split_whitespace().next()?.parse::<usize>().ok());
    let start = (0..lines.len()).rev().find(|&i| lines[i].contains("2theta/TOF") && lines[i].contains("Hw") && first_row(i) == Some(1)).unwrap();
    let mut v: Vec<(f64, f64, f64, f64)> = Vec::new();
    for l in &lines[start + 1..] {
        if l.contains("-----") || l.contains("=>") {
            break;
        }
        let f: Vec<f64> = l.split_whitespace().filter_map(|s| s.parse().ok()).collect();
        if f.len() < 10 {
            continue;
        }
        let (hw, eta, tt, ic) = (f[6], f[7], f[8], f[9]);
        match v.last_mut() {
            Some(last) if (last.0 - tt).abs() < 1e-3 => last.3 += ic,
            _ => v.push((tt, hw, eta, ic)),
        }
    }
    v
}

/// Y₂O₃'s isolated reflections (`y2o3_structure_calc.out`, Kα₁/Kα₂ with ratio 0.5) fitted one by one, against
/// FullProf's whole-pattern profile at its refined values, as measured in L260: (peaks with A at 20 esd's or more,
/// the largest |x₀ − (2θ + Zero)| in esd's, the largest relative difference of H from FullProf's Hw, of A from ΣIcalc).
/// FullProf's peaks are its pseudo-Voigt with η up to 1.0076 (above 1, which the fit does not allow), the
/// Bérar–Baldinozzi asymmetry below 30° and widths from U, V, W across the pattern; the fit is one free pseudo-Voigt.
const FULLPROF_PEAKS: (usize, f64, f64, f64) = (23, 4.2, 0.195, 0.088);

#[test]
fn fullprof_peaks() {
    let d = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs/Y2O3");
    let text: String = std::fs::read(d.join("y2o3.dat")).unwrap().iter().map(|&b| b as char).collect();
    let data = drs::data::read(&text, 0).unwrap();
    let out: String = std::fs::read(d.join("y2o3_structure_calc.out")).unwrap().iter().map(|&b| b as char).collect();
    let pos = fullprof_positions(&out);
    println!("{} positions, {} points", pos.len(), data.x.len());
    let (l1, l2, ratio) = (1.540600, 1.544400, 0.5);
    // FullProf's table lists the Bragg positions; the peaks are at 2θ + Zero (its .pcr's 0.00250).
    let zero = 0.0025;
    let shift = |t: f64| 2.0 * ((l2 / l1) * (t.to_radians() / 2.0).sin()).asin().to_degrees();
    let mut strong: Vec<(f64, f64, f64, f64)> = Vec::new();
    for (k, &(tt, hw, eta, ic)) in pos.iter().enumerate() {
        let (lo, hi) = (tt - 4.0 * hw, shift(tt) + 4.0 * hw);
        let crowded = pos.iter().enumerate().any(|(j, p)| j != k && p.0 > lo - 4.0 * p.1 && p.0 < hi + 4.0 * p.1);
        if crowded || tt < 30.0 {
            continue;
        }
        let r = match fit(&data.x, &data.y, &data.sigma, &PeakRequest { from: lo, to: hi, shape: PeakShape::PseudoVoigt, doublet: Some((l2 / l1, ratio)) }) {
            Ok(r) => r,
            Err(e) => {
                println!("{tt}: {e}");
                continue;
            }
        };
        let tt = tt + zero;
        println!(
            "2θ {tt:8.3}: x₀ − 2θ {:+.4} ({:.1} esd); H {:.4} against Hw {hw:.4} ({:+.1} %); η {:.3}({:.3}) against {eta:.3}; A {:.1} ΣIcalc {ic:.1}; χ² {:.2}",
            r.position.value - tt,
            (r.position.value - tt) / r.position.esd.unwrap(),
            r.fwhm.value,
            100.0 * (r.fwhm.value / hw - 1.0),
            r.eta.value,
            r.eta.esd.unwrap(),
            r.intensity.value,
            r.chi2
        );
        // The peaks the data determine well (A at 20 esd's or more): position in esd's, width and intensity.
        if r.intensity.value > 20.0 * r.intensity.esd.unwrap() {
            strong.push((tt, (r.position.value - tt) / r.position.esd.unwrap(), r.fwhm.value / hw - 1.0, r.intensity.value / ic - 1.0));
        }
    }
    let worst = |k: usize| strong.iter().map(|s| [s.1, s.2, s.3][k].abs()).fold(0.0f64, f64::max);
    println!("{} peaks at 20 esd's or more: x₀ within {:.2} esd of FullProf's 2θ + Zero, H within {:.3}, A within {:.3} of FullProf's Hw and ΣIcalc", strong.len(), worst(0), worst(1), worst(2));
    assert_eq!(strong.len(), FULLPROF_PEAKS.0);
    assert!(worst(0) <= FULLPROF_PEAKS.1 && worst(1) <= FULLPROF_PEAKS.2 && worst(2) <= FULLPROF_PEAKS.3, "{strong:?}");
}

/// The app's request 38: the doublet from the scan's wavelengths. A PANalytical scan (GSAS-II's synthetic .xrdml)
/// states Cu Kα₁, Kα₂ and their ratio; a scan with one wavelength, or a λ₂ without its ratio, gives one peak.
#[test]
fn doublet_from_the_scan() {
    use drs::formats::Wavelength;
    use drs::peakfit::doublet_of;
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/gsas2/synthetic/scan.xrdml");
    let d = drs::formats::read_path(&p, None).unwrap();
    let w = d.scans[0].wavelength.unwrap();
    println!("scan.xrdml: {w:?}");
    let (rl, ratio) = doublet_of(&w).expect("the scan states λ2 and the ratio");
    assert_eq!((rl, ratio), (w.lambda2.unwrap() / w.lambda1, w.ratio.unwrap()));
    assert!((rl - 1.0025).abs() < 1e-3 && ratio > 0.0);
    assert_eq!(doublet_of(&Wavelength::single(1.5406)), None);
    assert_eq!(doublet_of(&Wavelength { lambda1: 1.5406, lambda2: Some(1.5444), ratio: None }), None);
    let r = PeakRequest::with_wavelength(39.3, 40.9, PeakShape::PseudoVoigt, &Wavelength { lambda1: 1.540560, lambda2: Some(1.544390), ratio: Some(0.5) });
    assert_eq!(r.doublet, Some((1.544390 / 1.540560, 0.5)));
    // The same fit as the doublet given by hand.
    let mut noise = Noise(7);
    let (x, y, s) = synthetic(&mut noise, r.doublet);
    let a = fit(&x, &y, &s, &r).unwrap();
    let b = fit(&x, &y, &s, &PeakRequest { from: 39.3, to: 40.9, shape: PeakShape::PseudoVoigt, doublet: Some((1.544390 / 1.540560, 0.5)) }).unwrap();
    assert_eq!(a, b);
}
