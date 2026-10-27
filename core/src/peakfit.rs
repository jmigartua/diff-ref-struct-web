//! Single-peak fitting (P51; WinPLOTR's interactive peak fitting, FullProf Suite): one peak over a linear background,
//! fitted to the points of a 2θ range by weighted least squares, each parameter with its standard deviation.
//!
//! The peak is the pseudo-Voigt the core's profiles use (`profile::pseudo_voigt`, FullProf manual eq. 3.13, unit area
//! in 2θ): y(x) = A·pV(x − x₀; η, H) + b₀ + b₁(x − x_c), x_c the range's centre; for a Kα₁/Kα₂ doublet, plus R·A·pV(x −
//! x₂; η, H) with x₂ = 2 arcsin((λ₂/λ₁) sin(x₀/2)), the same H and η, R the intensity ratio. The free parameters are
//! the position x₀, the integrated intensity A, the FWHM H, η (or η fixed at 0, a Gaussian, or 1, a Lorentzian) and
//! the background b₀, b₁. Levenberg–Marquardt with exact derivatives (dual numbers); weights 1/σ²; the esd's are
//! √(C_kk χ²_ν), C the inverse normal matrix and χ²_ν = S/(N − P), as the engine's (`lsq`).
//!
//! Checked in `peak_fit`: on synthetic peaks with Poisson noise, the parameters within one esd of the truth in 68 % of
//! the trials (the esd's calibrated); on LaB6's and Y₂O₃'s isolated reflections, the positions and widths against
//! FullProf's.

use crate::dual::{Dual, Real};
use crate::profile::pseudo_voigt;

/// The peak's shape.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum PeakShape {
    /// η refined in [0, 1].
    PseudoVoigt,
    Gaussian,
    Lorentzian,
}

/// What a fit is asked: the range, the shape, and the second wavelength of a doublet (λ₂/λ₁ and the ratio).
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PeakRequest {
    pub from: f64,
    pub to: f64,
    pub shape: PeakShape,
    pub doublet: Option<(f64, f64)>,
}

/// The doublet a scan's wavelengths make (the app's request 38): (λ₂/λ₁, the intensity ratio I(λ₂)/I(λ₁)) when the
/// scan states λ₂ and the ratio, both positive and λ₂ ≠ λ₁; None for one wavelength, or a λ₂ without its ratio (no
/// ratio is assumed).
pub fn doublet_of(w: &crate::formats::Wavelength) -> Option<(f64, f64)> {
    let (l2, ratio) = (w.lambda2?, w.ratio?);
    (w.lambda1 > 0.0 && l2 > 0.0 && ratio > 0.0 && (l2 - w.lambda1).abs() > 1e-9).then(|| (l2 / w.lambda1, ratio))
}

impl PeakRequest {
    /// A request with the doublet of the scan's wavelengths (`doublet_of`), or one peak when the scan states one.
    pub fn with_wavelength(from: f64, to: f64, shape: PeakShape, w: &crate::formats::Wavelength) -> PeakRequest {
        PeakRequest { from, to, shape, doublet: doublet_of(w) }
    }
}

/// A value and its standard deviation (None when the parameter was held).
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Estimate {
    pub value: f64,
    pub esd: Option<f64>,
}

/// The fitted peak.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PeakFit {
    /// x₀ (degrees 2θ).
    pub position: Estimate,
    /// A, the integrated intensity (counts·degrees) of the first wavelength's peak.
    pub intensity: Estimate,
    /// H (degrees).
    pub fwhm: Estimate,
    pub eta: Estimate,
    /// b₀ and b₁ about the range's centre.
    pub background: [Estimate; 2],
    pub centre: f64,
    /// The points used, χ² = S/(N − P), the iterations and whether the shifts fell below 10⁻⁴ of their esd's.
    pub points: usize,
    pub chi2: f64,
    pub iterations: usize,
    pub converged: bool,
    /// The fit at the range's points.
    pub x: Vec<f64>,
    pub ycalc: Vec<f64>,
}

const N: usize = 6;

fn model(p: &[Dual<N>; N], x: f64, xc: f64, doublet: Option<(f64, f64)>) -> Dual<N> {
    let [x0, a, h, eta, b0, b1] = *p;
    let mut y = a * pseudo_voigt(Dual::<N>::cst(x) - x0, eta, h) + b0 + b1 * (x - xc);
    if let Some((rl, ratio)) = doublet {
        let half = (x0 * (std::f64::consts::PI / 360.0)).sin() * rl;
        if half.val().abs() < 1.0 {
            let x2 = half.asin() * (360.0 / std::f64::consts::PI);
            y = y + a * ratio * pseudo_voigt(Dual::<N>::cst(x) - x2, eta, h);
        }
    }
    y
}

fn invert(mut a: Vec<Vec<f64>>) -> Option<Vec<Vec<f64>>> {
    let n = a.len();
    let mut inv: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect()).collect();
    for c in 0..n {
        let p = (c..n).max_by(|&i, &j| a[i][c].abs().total_cmp(&a[j][c].abs()))?;
        if a[p][c].abs() < 1e-300 {
            return None;
        }
        a.swap(c, p);
        inv.swap(c, p);
        let d = a[c][c];
        for j in 0..n {
            a[c][j] /= d;
            inv[c][j] /= d;
        }
        for r in 0..n {
            if r != c {
                let f = a[r][c];
                for j in 0..n {
                    a[r][j] -= f * a[c][j];
                    inv[r][j] -= f * inv[c][j];
                }
            }
        }
    }
    Some(inv)
}

/// Fits one peak to the points of x, y, σ (σ > 0) within the request's range.
pub fn fit(x: &[f64], y: &[f64], sigma: &[f64], req: &PeakRequest) -> Result<PeakFit, String> {
    let idx: Vec<usize> = (0..x.len()).filter(|&i| x[i] >= req.from && x[i] <= req.to && sigma[i] > 0.0).collect();
    let n = idx.len();
    let free: Vec<usize> = match req.shape {
        PeakShape::PseudoVoigt => vec![0, 1, 2, 3, 4, 5],
        _ => vec![0, 1, 2, 4, 5],
    };
    if n < free.len() + 3 {
        return Err(format!("{n} points between {} and {}: too few for a peak", req.from, req.to));
    }
    let xs: Vec<f64> = idx.iter().map(|&i| x[i]).collect();
    let ys: Vec<f64> = idx.iter().map(|&i| y[i]).collect();
    let w: Vec<f64> = idx.iter().map(|&i| 1.0 / (sigma[i] * sigma[i])).collect();
    let xc = 0.5 * (xs[0] + xs[n - 1]);
    // Starting values: the background through the ends, the highest point, its half width, η = 0.5.
    let k = (n / 10).max(1);
    let (ya, yb) = (ys[..k].iter().sum::<f64>() / k as f64, ys[n - k..].iter().sum::<f64>() / k as f64);
    let (xa, xb) = (xs[..k].iter().sum::<f64>() / k as f64, xs[n - k..].iter().sum::<f64>() / k as f64);
    let b1 = (yb - ya) / (xb - xa);
    let b0 = ya + b1 * (xc - xa);
    let net: Vec<f64> = (0..n).map(|i| ys[i] - b0 - b1 * (xs[i] - xc)).collect();
    let top = (0..n).max_by(|&i, &j| net[i].total_cmp(&net[j])).unwrap();
    if net[top] <= 0.0 {
        return Err(format!("no peak above the background between {} and {}", req.from, req.to));
    }
    let above: Vec<f64> = (0..n).filter(|&i| net[i] >= 0.5 * net[top]).map(|i| xs[i]).collect();
    let step = (xs[n - 1] - xs[0]) / (n - 1) as f64;
    let h0 = (above[above.len() - 1] - above[0]).max(2.0 * step);
    let eta0 = match req.shape {
        PeakShape::PseudoVoigt => 0.5,
        PeakShape::Gaussian => 0.0,
        PeakShape::Lorentzian => 1.0,
    };
    let mut p = [xs[top], net[top] * h0 * 1.2, h0, eta0, b0, b1];
    let chi = |p: &[f64; N]| -> (f64, Vec<f64>, Vec<Vec<f64>>, Vec<f64>) {
        let d: [Dual<N>; N] = std::array::from_fn(|i| Dual::var(p[i], i));
        let (mut s, mut g, mut a, mut yc) = (0.0, vec![0.0; free.len()], vec![vec![0.0; free.len()]; free.len()], Vec::with_capacity(n));
        for i in 0..n {
            let m = model(&d, xs[i], xc, req.doublet);
            let r = ys[i] - m.v;
            s += w[i] * r * r;
            for (u, &fu) in free.iter().enumerate() {
                g[u] += w[i] * r * m.d[fu];
                for (v, &fv) in free.iter().enumerate() {
                    a[u][v] += w[i] * m.d[fu] * m.d[fv];
                }
            }
            yc.push(m.v);
        }
        (s, g, a, yc)
    };
    let (mut s, mut g, mut a, mut yc) = chi(&p);
    let mut lambda = 1e-3;
    let (mut iterations, mut converged) = (0, false);
    let dof = (n - free.len()) as f64;
    while iterations < 200 {
        iterations += 1;
        let mut aa = a.clone();
        for (u, row) in aa.iter_mut().enumerate() {
            row[u] *= 1.0 + lambda;
        }
        let Some(inv) = invert(aa) else { return Err("the normal matrix is singular".into()) };
        let delta: Vec<f64> = (0..free.len()).map(|u| (0..free.len()).map(|v| inv[u][v] * g[v]).sum()).collect();
        let mut q = p;
        for (u, &fu) in free.iter().enumerate() {
            q[fu] += delta[u];
        }
        q[2] = q[2].max(1e-6);
        q[3] = q[3].clamp(0.0, 1.0);
        let (s2, g2, a2, yc2) = chi(&q);
        if s2 <= s {
            let small = invert(a.clone()).is_some_and(|c| (0..free.len()).all(|u| delta[u].abs() <= 1e-4 * (c[u][u] * s / dof).sqrt().max(1e-300)));
            (p, s, g, a, yc) = (q, s2, g2, a2, yc2);
            lambda = (lambda * 0.3).max(1e-12);
            if small {
                converged = true;
                break;
            }
        } else {
            lambda *= 10.0;
            if lambda > 1e12 {
                converged = true;
                break;
            }
        }
    }
    let c = invert(a).ok_or("the normal matrix is singular at the end")?;
    let chi2 = s / dof;
    let est = |k: usize| Estimate { value: p[k], esd: free.iter().position(|&f| f == k).map(|u| (c[u][u] * chi2).sqrt()) };
    Ok(PeakFit {
        position: est(0),
        intensity: est(1),
        fwhm: est(2),
        eta: est(3),
        background: [est(4), est(5)],
        centre: xc,
        points: n,
        chi2,
        iterations,
        converged,
        x: xs,
        ycalc: yc,
    })
}
