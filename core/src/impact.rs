//! Refinement diagnostics (P45, L272): what each parameter does to the fit, and how the weighted residuals are
//! distributed.
//!
//! **A parameter's impact.** With the weighted residuals rᵢ = yᵢ − ycᵢ, weights wᵢ, S = Σwᵢrᵢ² and χ² = S/(N − P + C)
//! (`lsq::statistics`), and gᵢ = ∂ycᵢ/∂p the parameter's column of the Jacobian (its symmetry ties moving with it):
//! - the slope ∂χ²/∂p = −2 Σwᵢrᵢgᵢ/(N − P + C), analytic. GSAS-II's "Parameter Impact" (Toby 2024, J. Appl. Cryst.
//!   57, 175–180; its `AllPrmDerivs`, GSASIIstrMain.py) ranks the parameters by the finite-difference slope of the
//!   rms weighted residual, −∂√(S/N)/∂p, at the step max(10⁻⁴|p|, 10⁻⁶); that is −(N − P + C)/(2N·rms) times this
//!   slope (test `refinement_diagnostics`);
//! - a slope has the parameter's unit, so the core ranks by what refining the parameter alone would do, to first
//!   order (one Gauss–Newton step): it would move by δ = Σwrg/Σwg², which is δ/σ₁ of its esd's alone,
//!   σ₁ = √(χ²/Σwg²), and χ² would fall by (Σwrg)²/Σwg²/(N − P + C) = (δ/σ₁)²·χ²/(N − P + C). A refined parameter
//!   at convergence has δ = 0 (the slope is zero at a minimum);
//! - where it moves the calculated pattern: gᵢδ point by point, its largest |gᵢδ| and the 2θ ranges where |gᵢ| is
//!   at least a tenth of its largest; and the largest change for one esd, max|gᵢ|σ₁.
//!
//! The restraints' terms (P40) are not in S here: a restrained parameter's impact is that of the pattern alone.
//!
//! **The normal probability plot** (Abrahams & Keve 1971, Acta Cryst. A27, 157–165; GSAS-II's "Error Analysis",
//! `PlotDeltSig` in GSASIIplot.py): the weighted residuals √wᵢ(yᵢ − ycᵢ) of the points used, sorted, against the
//! expected values of the order statistics of as many standard normal values (`expected_normal_order_statistics`,
//! exact to 10⁻⁹). For residuals that are normal with the stated σ, the points lie on a line of slope 1 through 0;
//! the slope over the central part (expected value between −1 and 1, about 68 % of the points, as GSAS-II fits it)
//! is the goodness of fit of the well-fitted points, and the intercept a bias. GSAS-II takes its expected values from
//! a rational approximation at equally spaced probabilities (Abramowitz & Stegun 26.2.23); the difference this makes
//! is measured in `refinement_diagnostics`. Neither paper is in `manuals/papers/` (L272): the definitions are those
//! GSAS-II's code and help implement and cite.

use crate::data::Pattern;
use crate::jacobian::Field;
use crate::powder::Model;
use crate::refinable::ParameterSet;

/// What refining one parameter would do, to first order.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Impact {
    /// The parameter's address (`phase1.cell_a`), or the refined parameter's name.
    pub name: String,
    pub refined: bool,
    pub value: f64,
    /// ∂χ²/∂p.
    pub slope: f64,
    /// σ₁ = √(χ²/Σwg²): its esd if it were refined alone.
    pub esd_alone: f64,
    /// δ = Σwrg/Σwg²: the shift refining it alone would make.
    pub shift: f64,
    /// δ/σ₁: that shift in its esd's (signed).
    pub esds: f64,
    /// The fall of χ² that shift gives, to first order.
    pub gain: f64,
    /// max|gᵢδ| (counts) and its 2θ.
    pub largest: f64,
    pub at: f64,
    /// max|gᵢ|σ₁: the largest change of the calculated pattern for one esd.
    pub per_esd: f64,
    /// The 2θ ranges where |gᵢ| ≥ max|g|/10.
    pub regions: Vec<(f64, f64)>,
    /// gᵢδ point by point (`impact_of` fills it; `impacts` leaves it empty).
    pub change: Vec<f64>,
}

/// A parameter weighed: its name, whether it is refined, and the values it moves with their multipliers.
type Candidate = (String, bool, Vec<(Field, f64)>);

/// The parameters to weigh: each free parameter of the set (refined, with its fields and multipliers) and each value
/// of the model that is not refined and not fixed by the symmetry (its ties moving with it).
fn candidates(model: &Model, set: &ParameterSet) -> Vec<Candidate> {
    let mut v: Vec<Candidate> = set.free.iter().map(|f| (f.name.clone(), true, f.fields.clone())).collect();
    let refined = set.fields();
    for fi in model.fields() {
        if fi.fixed_by.is_some() || refined.contains(&fi.field) {
            continue;
        }
        if let Ok(one) = ParameterSet::from_addresses(std::slice::from_ref(&fi.address), &[], &[], model) {
            if let Some(f) = one.free.first() {
                v.push((fi.address.clone(), false, f.fields.clone()));
            }
        }
    }
    v
}

fn regions(x: &[f64], g: &[f64], w: &[f64]) -> Vec<(f64, f64)> {
    let top = g.iter().zip(w).filter(|(_, &w)| w > 0.0).map(|(g, _)| g.abs()).fold(0.0, f64::max);
    let mut out: Vec<(f64, f64)> = Vec::new();
    if top <= 0.0 {
        return out;
    }
    let mut start: Option<usize> = None;
    for i in 0..=x.len() {
        let on = i < x.len() && w[i] > 0.0 && g[i].abs() >= 0.1 * top;
        match (on, start) {
            (true, None) => start = Some(i),
            (false, Some(s)) => {
                out.push((x[s], x[i - 1]));
                start = None;
            }
            _ => {}
        }
    }
    out
}

fn impacts_with(model: &Model, data: &Pattern, set: &ParameterSet, only: Option<&str>) -> Result<Vec<Impact>, String> {
    let mut cands = candidates(model, set);
    if let Some(name) = only {
        cands.retain(|c| c.0 == name);
        if cands.is_empty() {
            return Err(format!("{name}: no such parameter, or fixed by the symmetry"));
        }
    }
    let mut fields: Vec<Field> = cands.iter().flat_map(|c| c.2.iter().map(|f| f.0)).collect();
    fields.sort();
    fields.dedup();
    let (c, jac) = model.derivatives(&data.x, &fields);
    let w = crate::lsq::weights(model, data, &c.ycalc);
    let stats = crate::lsq::statistics(model, data, &c.ycalc, set.free.len());
    let dof = (stats.n + stats.c).saturating_sub(stats.p).max(1) as f64;
    let chi2 = stats.chi2;
    let r: Vec<f64> = data.y.iter().zip(&c.ycalc).map(|(y, yc)| y - yc).collect();
    let ymax = c.ycalc.iter().map(|v| v.abs()).fold(0.0, f64::max);
    let mut out = Vec::new();
    for (name, refined, fs) in cands {
        let mut g = vec![0.0; data.x.len()];
        for (f, m) in &fs {
            let col = &jac[fields.iter().position(|x| x == f).unwrap()];
            g.iter_mut().zip(col).for_each(|(gi, ci)| *gi += m * ci);
        }
        let value = model.get(fs[0].0) / fs[0].1;
        // A column that is rounding noise (below 10⁻⁹ of the largest Ycal for a change of max(|p|, 1)) moves nothing:
        // Y₂O₃'s March–Dollase G1 at 1, whose derivative the average over a cubic reflection's equivalents cancels.
        let gmax_all = g.iter().zip(&w).filter(|(_, &w)| w > 0.0).map(|(g, _)| g.abs()).fold(0.0, f64::max);
        if gmax_all * value.abs().max(1.0) <= 1e-9 * ymax {
            g.iter_mut().for_each(|x| *x = 0.0);
        }
        let (b, a): (f64, f64) = (0..g.len()).fold((0.0, 0.0), |(b, a), i| (b + w[i] * r[i] * g[i], a + w[i] * g[i] * g[i]));
        let (shift, esd_alone) = if a > 0.0 { (b / a, (chi2 / a).sqrt()) } else { (0.0, f64::INFINITY) };
        let esds = if esd_alone.is_finite() && esd_alone > 0.0 { shift / esd_alone } else { 0.0 };
        let (mut largest, mut at, mut gmax) = (0.0, f64::NAN, 0.0f64);
        for i in 0..g.len() {
            if w[i] > 0.0 && (g[i] * shift).abs() > largest {
                (largest, at) = ((g[i] * shift).abs(), data.x[i]);
            }
            if w[i] > 0.0 {
                gmax = gmax.max(g[i].abs());
            }
        }
        let gain = if a > 0.0 { b * b / a / dof } else { 0.0 };
        out.push(Impact {
            regions: regions(&data.x, &g, &w),
            change: if only.is_some() { g.iter().map(|gi| gi * shift).collect() } else { Vec::new() },
            name,
            refined,
            value,
            slope: -2.0 * b / dof,
            esd_alone,
            shift,
            esds,
            gain,
            largest,
            at,
            per_esd: if esd_alone.is_finite() { gmax * esd_alone } else { 0.0 },
        });
    }
    out.sort_by(|p, q| q.gain.partial_cmp(&p.gain).unwrap_or(std::cmp::Ordering::Equal));
    Ok(out)
}

/// Every parameter's impact (P45), refined or not, the largest fall of χ² first; without the point-by-point change.
pub fn impacts(model: &Model, data: &Pattern, set: &ParameterSet) -> Result<Vec<Impact>, String> {
    impacts_with(model, data, set, None)
}

/// One parameter's impact, by its name or address, with the change of the calculated pattern point by point.
pub fn impact_of(model: &Model, data: &Pattern, set: &ParameterSet, name: &str) -> Result<Impact, String> {
    Ok(impacts_with(model, data, set, Some(name))?.remove(0))
}

/// The normal probability plot of the weighted residuals.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct NormalProbability {
    /// √wᵢ(yᵢ − ycᵢ) of the points used, sorted.
    pub residuals: Vec<f64>,
    /// The expected values of the normal order statistics, as many.
    pub expected: Vec<f64>,
    /// The straight line through the central points (expected value in [−1, 1]): residual = slope·expected + intercept.
    pub slope: f64,
    pub intercept: f64,
    /// How many points the line is fitted to.
    pub central: usize,
}

/// The normal probability plot of a calculated pattern against its data (P45).
pub fn normal_probability(model: &Model, data: &Pattern, ycalc: &[f64]) -> NormalProbability {
    let w = crate::lsq::weights(model, data, ycalc);
    let mut d: Vec<f64> = (0..w.len()).filter(|&i| w[i] > 0.0).map(|i| w[i].sqrt() * (data.y[i] - ycalc[i])).collect();
    d.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let e = expected_normal_order_statistics(d.len());
    let (slope, intercept, central) = central_line(&e, &d);
    NormalProbability { residuals: d, expected: e, slope, intercept, central }
}

/// The least-squares line y = slope·x + intercept through the points with |x| ≤ 1.
pub fn central_line(x: &[f64], y: &[f64]) -> (f64, f64, usize) {
    let pts: Vec<(f64, f64)> = x.iter().zip(y).filter(|(x, _)| x.abs() <= 1.0).map(|(a, b)| (*a, *b)).collect();
    let n = pts.len() as f64;
    if pts.len() < 2 {
        return (f64::NAN, f64::NAN, pts.len());
    }
    let (mx, my) = (pts.iter().map(|p| p.0).sum::<f64>() / n, pts.iter().map(|p| p.1).sum::<f64>() / n);
    let sxx: f64 = pts.iter().map(|p| (p.0 - mx).powi(2)).sum();
    let sxy: f64 = pts.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum();
    let slope = sxy / sxx;
    (slope, my - slope * mx, pts.len())
}

/// ln Φ(x), the standard normal distribution function, to about 10⁻¹⁴ relative: Marsaglia's series
/// Φ(x) = ½ + φ(x)(x + x³/3 + x⁵/(3·5) + …) (J. Stat. Softw. 11(4), 2004) for x ≥ −3, and for x < −3 the continued
/// fraction of the Mills ratio, Φ(x) = φ(x)/(z + 1/(z + 2/(z + 3/(z + …)))), z = −x.
pub fn ln_normal_cdf(x: f64) -> f64 {
    let ln_phi = -0.5 * x * x - 0.5 * (2.0 * std::f64::consts::PI).ln();
    if x >= -3.0 {
        if x > 8.0 {
            // Φ = 1 − Φ(−x), ln(1 − ε) ≈ −ε.
            return -ln_normal_cdf(-x).exp();
        }
        let (mut term, mut sum) = (x, x);
        let mut k = 1.0;
        while term.abs() > 1e-17 * sum.abs().max(1e-300) {
            k += 2.0;
            term *= x * x / k;
            sum += term;
        }
        (0.5 + ln_phi.exp() * sum).ln()
    } else {
        let z = -x;
        // Lentz's evaluation of z + 1/(z + 2/(z + 3/(z + …))), from the bottom: 200 terms are many at z ≥ 3.
        let mut f = z;
        for k in (1..=200).rev() {
            f = z + k as f64 / f;
        }
        ln_phi - f.ln()
    }
}

/// Φ⁻¹(p): Newton's method on `ln_normal_cdf` from a rough start, to its accuracy.
pub fn normal_quantile(p: f64) -> f64 {
    if p <= 0.0 {
        return f64::NEG_INFINITY;
    }
    if p >= 1.0 {
        return f64::INFINITY;
    }
    if p > 0.5 {
        return -normal_quantile(1.0 - p);
    }
    // A start from the tail's asymptote, then Newton on ln Φ(x) = ln p: x ← x − (lnΦ(x) − ln p)·Φ(x)/φ(x).
    let lp = p.ln();
    let mut x = if p > 0.3 { (p - 0.5) * (2.0 * std::f64::consts::PI).sqrt() } else { -(-2.0 * lp - (-2.0 * lp).ln() - (2.0 * std::f64::consts::PI).ln()).max(0.0).sqrt() };
    for _ in 0..100 {
        let lc = ln_normal_cdf(x);
        let ln_phi = -0.5 * x * x - 0.5 * (2.0 * std::f64::consts::PI).ln();
        let step = (lc - lp) * (lc - ln_phi).exp();
        x -= step;
        if step.abs() <= 1e-15 * x.abs().max(1.0) {
            break;
        }
    }
    x
}

/// The expected values of the order statistics of n standard normal values, smallest first:
/// E[X₍ᵢ₎] = ∫ x · n!/((i−1)!(n−i)!) Φ(x)^{i−1} (1 − Φ(x))^{n−i} φ(x) dx, by the trapezoidal rule on a grid fine for
/// n (spacing min(0.002, 0.1·1.25/√n), where the central order statistic's spread is 1.25/√n), each density
/// normalised by its own sum and taken where its logarithm is within 60 of its largest. Exact for n = 2 and 3
/// (±1/√π, ±3/(2√π)) and for n = 10 against mpmath's quadrature at 30 digits, to 10⁻⁹ (test `refinement_diagnostics`).
pub fn expected_normal_order_statistics(n: usize) -> Vec<f64> {
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![0.0];
    }
    let h = (0.125 / (n as f64).sqrt()).min(0.002);
    let lim = 10.0;
    let m = (2.0 * lim / h).ceil() as usize + 1;
    let xs: Vec<f64> = (0..m).map(|k| -lim + k as f64 * h).collect();
    let lc: Vec<f64> = xs.iter().map(|&x| ln_normal_cdf(x)).collect();
    let lu: Vec<f64> = xs.iter().map(|&x| ln_normal_cdf(-x)).collect();
    let mut e = vec![0.0; n];
    let half = n.div_ceil(2);
    for i in 1..=half {
        let (a, b) = ((i - 1) as f64, (n - i) as f64);
        let lw = |k: usize| a * lc[k] + b * lu[k] - 0.5 * xs[k] * xs[k];
        // Start at the median of the order statistic, about Φ⁻¹(i/(n + 1)), and walk out both ways.
        let centre = normal_quantile(i as f64 / (n as f64 + 1.0)).clamp(-lim, lim);
        let kc = (((centre + lim) / h).round() as usize).min(m - 1);
        let mut top = lw(kc);
        let mut terms: Vec<(f64, f64)> = vec![(xs[kc], top)];
        for dir in [-1i64, 1] {
            let mut k = kc as i64 + dir;
            while k >= 0 && (k as usize) < m {
                let v = lw(k as usize);
                top = top.max(v);
                if v < top - 60.0 {
                    break;
                }
                terms.push((xs[k as usize], v));
                k += dir;
            }
        }
        let (s0, s1) = terms.iter().fold((0.0, 0.0), |(s0, s1), &(x, v)| {
            let wv = (v - top).exp();
            (s0 + wv, s1 + wv * x)
        });
        e[i - 1] = s1 / s0;
        e[n - i] = -e[i - 1];
    }
    if n % 2 == 1 {
        e[n / 2] = 0.0;
    }
    e
}
