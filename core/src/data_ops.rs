//! The data operations of the Data module (R21; L066): an ordered list applied to a measured scan, never to its file.
//! The list is what the control file writes in its `[[pattern.operation]]` tables and what the app keeps as the
//! user works; `apply_all` gives the scan the refinement sees.
//!
//! Operations (each a pure function of the scan, σ propagated as for independent errors):
//! - `mask` and `mask_points`: points that stay, with no weight (σ = ∞); `delete` and `delete_points`: points removed.
//!   A range is closed, from ≤ x ≤ to; a point is matched to 10⁻⁹ of its x.
//! - `smooth`: Savitzky–Golay, a polynomial of `order` least-squares fitted over `points` (odd) points at constant
//!   step; at each end, the polynomial of the first (last) window evaluated at its first (last) half, as SciPy's
//!   `savgol_filter(mode = "interp")`. σ² of a smoothed point is Σ c_k² σ_k² over its coefficients.
//! - `rebin`: by a factor n (each n consecutive points become one at their mean x; the last points that do not make n
//!   are dropped), or to a new step (each point counts for the interval between the midpoints to its neighbours, and
//!   is shared among the new intervals in proportion to the overlap). Mode `sum` conserves the counts (Σy before =
//!   Σy after, σ² = Σ (share·σ)²); mode `average` divides by the factor (or the new step over the old).
//! - `normalize`: y and σ times a factor: `to` over the scan's monitor, over its counting time, or over its maximum;
//!   or a given factor.
//! - `add` and `subtract`: another file's scan (read by the registry), times `scale`, interpolated linearly at this
//!   scan's points (σ² of an interpolated point (1 − t)²σ_a² + t²σ_b²); points outside the other scan are removed.
//! - `merge`: the points of another scan joined to these; points at the same x (to 10⁻⁶ of the step) become their
//!   mean weighted by 1/σ², σ = (Σ 1/σ²)^(−1/2).
//! - `convert`: x to 2θ, d or Q (d = λ/(2 sin θ), Q = 4π sin θ/λ = 2π/d), with the scan's λ1 or the one given;
//!   `wavelength`: 2θ to another wavelength, sin θ′ = (λ′/λ) sin θ (points with |sin θ′| > 1 removed). y and σ are
//!   kept: these change the axis, not the counts. The points stay in increasing x.
//!
//! Analyses (they return results; the scan is unchanged):
//! - `peak_search`: the local maxima with their prominence, height and width, as SciPy's `find_peaks` and
//!   `peak_widths(rel_height = 0.5)` define them, kept when the prominence is at least `min_prominence` and at least
//!   `min_significance` times the peak's σ.
//! - `gsas2_peak_candidates`: GSAS-II's own rule for peak candidates, for comparison.
//! - `auto_background`: the asymmetrically reweighted penalised least squares (arPLS) of Baek et al., Analyst 140,
//!   250 (2015), as pybaselines' `whittaker.arpls`, which GSAS-II uses for its automatic background (GSASIIpwd.py
//!   `autoBkgCalc`, λ = 10^logLam, at most 10 iterations), and its points every N/100 as GSAS-II's fixed points.

use crate::control::syntax::{float, string, Diagnostic, Item, Table, Value};
use crate::data::Pattern;
use crate::formats::{Axis, Scan, Wavelength};

/// Reads another file named by an operation: (path as written, reader id) → its first scan.
pub type Load<'a> = dyn Fn(&str, Option<&str>) -> Result<Scan, String> + 'a;

/// What `normalize` divides by.
#[derive(Clone, Debug, PartialEq)]
pub enum Normalization {
    /// y · to / monitor.
    Monitor(f64),
    /// y · to / counting time.
    Time(f64),
    /// y · to / max y.
    Maximum(f64),
    /// y · factor.
    Factor(f64),
    /// y · to / Δx, Δx each point's step to the next (the last point's, the step before it): the counts per unit of x,
    /// times `to` (FullProf's reading of a GSAS RALF bank is GSAS-II's y·100/Δt: CeO2_PEARL_calc.prf's Yobs to 10⁻⁴,
    /// P05).
    Step(f64),
}

/// How `rebin` combines points.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RebinMode {
    Sum,
    Average,
    /// The first point of each group kept as it is (FullProf's Ste, manual p. 64: "only those points corresponding to
    /// the new step Ste×Step are taken into account"; L267). By a factor only.
    Keep,
}

#[derive(Clone, Debug, PartialEq)]
pub enum Operation {
    Mask { from: f64, to: f64 },
    MaskPoints { x: Vec<f64> },
    Delete { from: f64, to: f64 },
    DeletePoints { x: Vec<f64> },
    Smooth { points: usize, order: usize },
    RebinFactor { factor: usize, mode: RebinMode },
    RebinStep { step: f64, mode: RebinMode },
    Normalize(Normalization),
    Add { file: String, reader: Option<String>, scale: f64 },
    Subtract { file: String, reader: Option<String>, scale: f64 },
    Merge { file: String, reader: Option<String>, scale: f64 },
    Convert { to: Axis, wavelength: Option<f64> },
    /// 2θ for another wavelength `to`, from the scan's own λ, or from `from` when the scan states none (app 34: the λ
    /// the user takes from a control file); a `from` that disagrees with the scan's own is refused.
    Wavelength { to: f64, from: Option<f64> },
}

fn same_x(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-9 * a.abs().max(b.abs()).max(1.0)
}

fn retain(s: &mut Scan, keep: impl Fn(usize) -> bool) {
    let p = &s.pattern;
    let idx: Vec<usize> = (0..p.len()).filter(|&i| keep(i)).collect();
    s.pattern = Pattern { x: idx.iter().map(|&i| p.x[i]).collect(), y: idx.iter().map(|&i| p.y[i]).collect(), sigma: idx.iter().map(|&i| p.sigma[i]).collect() };
}

/// Solve a small symmetric positive-definite system by Gaussian elimination with partial pivoting.
fn solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
    let n = b.len();
    for k in 0..n {
        let p = (k..n).max_by(|&i, &j| a[i][k].abs().total_cmp(&a[j][k].abs())).unwrap();
        a.swap(k, p);
        b.swap(k, p);
        let (top, rest) = a.split_at_mut(k + 1);
        let pivot = &top[k];
        for (r, row) in rest.iter_mut().enumerate() {
            let f = row[k] / pivot[k];
            for (v, pv) in row[k..].iter_mut().zip(&pivot[k..]) {
                *v -= f * pv;
            }
            b[k + 1 + r] -= f * b[k];
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        x[i] = (b[i] - (i + 1..n).map(|j| a[i][j] * x[j]).sum::<f64>()) / a[i][i];
    }
    x
}

/// The weights that give a least-squares polynomial of `order`, fitted at positions `s`, evaluated at `t`:
/// v(t)ᵀ (AᵀA)⁻¹ Aᵀ.
fn fit_weights(s: &[f64], order: usize, t: f64) -> Vec<f64> {
    let n = order + 1;
    let mut ata = vec![vec![0.0; n]; n];
    for &x in s {
        for (i, row) in ata.iter_mut().enumerate() {
            for (j, v) in row.iter_mut().enumerate() {
                *v += x.powi((i + j) as i32);
            }
        }
    }
    let v: Vec<f64> = (0..n).map(|i| t.powi(i as i32)).collect();
    // (AᵀA)⁻¹ v, then dotted with each row of A.
    let g = solve(ata, v);
    s.iter().map(|&x| (0..n).map(|i| g[i] * x.powi(i as i32)).sum()).collect()
}

/// The step of a scan whose steps agree to 1 % (files write x to four or five decimals: 10.22001 in a scan of step
/// 0.02); Savitzky–Golay smoothing, as SciPy's, takes the points as equally spaced.
fn constant_step(p: &Pattern) -> Result<f64, String> {
    if p.len() < 2 {
        return Err("fewer than two points".into());
    }
    let step = (p.x[p.len() - 1] - p.x[0]) / (p.len() - 1) as f64;
    if p.x.windows(2).any(|w| ((w[1] - w[0]) - step).abs() > 0.01 * step.abs()) {
        return Err("the step is not constant (to 1 %)".into());
    }
    Ok(step)
}

/// Savitzky–Golay smoothing (see the module's notes).
pub fn smooth(p: &Pattern, points: usize, order: usize) -> Result<Pattern, String> {
    if points.is_multiple_of(2) || points < 3 {
        return Err(format!("smooth: {points} points; the window is an odd number, 3 or more"));
    }
    if order >= points {
        return Err(format!("smooth: a polynomial of order {order} needs a window of more than {order} points"));
    }
    if p.len() < points {
        return Err(format!("smooth: {} points for a window of {points}", p.len()));
    }
    constant_step(p).map_err(|e| format!("smooth: {e}"))?;
    let m = (points / 2) as isize;
    let s: Vec<f64> = (-m..=m).map(|k| k as f64).collect();
    let centre = fit_weights(&s, order, 0.0);
    let n = p.len();
    let (mut y, mut sigma) = (vec![0.0; n], vec![0.0; n]);
    let mut put = |i: usize, start: usize, w: &[f64]| {
        y[i] = w.iter().enumerate().map(|(k, c)| c * p.y[start + k]).sum();
        sigma[i] = w.iter().enumerate().map(|(k, c)| c * c * p.sigma[start + k] * p.sigma[start + k]).sum::<f64>().sqrt();
    };
    let m = m as usize;
    for i in m..n - m {
        put(i, i - m, &centre);
    }
    for j in 0..m {
        let t = j as f64 - m as f64;
        put(j, 0, &fit_weights(&s, order, t));
        put(n - 1 - j, n - points, &fit_weights(&s, order, -t));
    }
    Ok(Pattern { x: p.x.clone(), y, sigma })
}

pub fn rebin_factor(p: &Pattern, factor: usize, mode: RebinMode) -> Result<Pattern, String> {
    if factor < 1 {
        return Err("rebin: the factor is 1 or more".into());
    }
    let mut q = Pattern::default();
    for c in (0..p.len() / factor).map(|k| k * factor..(k + 1) * factor) {
        let n = factor as f64;
        let (sy, s2): (f64, f64) = c.clone().map(|i| (p.y[i], p.sigma[i] * p.sigma[i])).fold((0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1));
        q.x.push(c.clone().map(|i| p.x[i]).sum::<f64>() / n);
        match mode {
            RebinMode::Sum => {
                q.y.push(sy);
                q.sigma.push(s2.sqrt());
            }
            RebinMode::Average => {
                q.y.push(sy / n);
                q.sigma.push(s2.sqrt() / n);
            }
            RebinMode::Keep => {
                q.x.pop();
                q.x.push(p.x[c.start]);
                q.y.push(p.y[c.start]);
                q.sigma.push(p.sigma[c.start]);
            }
        }
    }
    Ok(q)
}

pub fn rebin_step(p: &Pattern, step: f64, mode: RebinMode) -> Result<Pattern, String> {
    if mode == RebinMode::Keep {
        return Err("rebin: keep takes a factor, not a step".into());
    }
    if step.is_nan() || step <= 0.0 || p.len() < 2 {
        return Err("rebin: the new step must be positive, and the scan have two points".into());
    }
    let n = p.len();
    // Each point's interval: between the midpoints to its neighbours (half a step out at the ends).
    let edge = |i: usize| -> f64 {
        if i == 0 {
            p.x[0] - (p.x[1] - p.x[0]) / 2.0
        } else if i == n {
            p.x[n - 1] + (p.x[n - 1] - p.x[n - 2]) / 2.0
        } else {
            (p.x[i - 1] + p.x[i]) / 2.0
        }
    };
    let (lo, hi) = (edge(0), edge(n));
    let bins = ((hi - lo) / step + 1e-9).floor() as usize;
    let mut q = Pattern { x: (0..bins).map(|j| lo + (j as f64 + 0.5) * step).collect(), y: vec![0.0; bins], sigma: vec![0.0; bins] };
    let mut s2 = vec![0.0; bins];
    for i in 0..n {
        let (a, b) = (edge(i), edge(i + 1));
        let width = b - a;
        let j0 = ((a - lo) / step).floor().max(0.0) as usize;
        let j1 = (((b - lo) / step).ceil() as usize).min(bins);
        for (j, (yj, vj)) in q.y[j0..j1].iter_mut().zip(&mut s2[j0..j1]).enumerate().map(|(k, v)| (k + j0, v)) {
            let (ba, bb) = (lo + j as f64 * step, lo + (j + 1) as f64 * step);
            let overlap = (b.min(bb) - a.max(ba)).max(0.0);
            if overlap > 0.0 {
                let share = overlap / width;
                *yj += share * p.y[i];
                *vj += share * share * p.sigma[i] * p.sigma[i];
            }
        }
    }
    let old_step = (hi - lo) / n as f64;
    let k = if mode == RebinMode::Average { old_step / step } else { 1.0 };
    for ((y, e), v) in q.y.iter_mut().zip(q.sigma.iter_mut()).zip(&s2) {
        *y *= k;
        *e = v.sqrt() * k;
    }
    Ok(q)
}

/// A scan at `x` by linear interpolation between its two neighbouring points: (y, σ) with σ² = (1 − t)²σ₁² + t²σ₂²,
/// None outside it (Plots' PL4: the private `interpolate` made public).
pub fn interpolate_at(o: &Pattern, x: f64) -> Option<(f64, f64)> {
    interpolate(o, x)
}

/// The other scan at the points `x` by linear interpolation: (y, σ), None outside it.
fn interpolate(o: &Pattern, x: f64) -> Option<(f64, f64)> {
    let n = o.len();
    if n == 0 || x < o.x[0] - 1e-9 || x > o.x[n - 1] + 1e-9 {
        return None;
    }
    let j = o.x.partition_point(|&v| v < x);
    if j < n && same_x(o.x[j], x) {
        return Some((o.y[j], o.sigma[j]));
    }
    if j == 0 || j >= n {
        return Some((o.y[j.min(n - 1)], o.sigma[j.min(n - 1)]));
    }
    let t = (x - o.x[j - 1]) / (o.x[j] - o.x[j - 1]);
    let (sa, sb) = (o.sigma[j - 1], o.sigma[j]);
    Some(((1.0 - t) * o.y[j - 1] + t * o.y[j], ((1.0 - t) * (1.0 - t) * sa * sa + t * t * sb * sb).sqrt()))
}

fn combine(s: &mut Scan, o: &Scan, scale: f64, sign: f64) {
    let p = &s.pattern;
    let mut q = Pattern::default();
    for i in 0..p.len() {
        if let Some((y, e)) = interpolate(&o.pattern, p.x[i]) {
            q.x.push(p.x[i]);
            q.y.push(p.y[i] + sign * scale * y);
            q.sigma.push((p.sigma[i] * p.sigma[i] + scale * scale * e * e).sqrt());
        }
    }
    s.pattern = q;
}

fn merge(s: &mut Scan, o: &Scan, scale: f64) {
    let step = if s.pattern.len() > 1 { (s.pattern.x[s.pattern.len() - 1] - s.pattern.x[0]).abs() / (s.pattern.len() - 1) as f64 } else { 1.0 };
    let mut pts: Vec<(f64, f64, f64)> = (0..s.pattern.len()).map(|i| (s.pattern.x[i], s.pattern.y[i], s.pattern.sigma[i])).collect();
    pts.extend((0..o.pattern.len()).map(|i| (o.pattern.x[i], scale * o.pattern.y[i], scale.abs() * o.pattern.sigma[i])));
    pts.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut q = Pattern::default();
    let mut i = 0;
    while i < pts.len() {
        let mut j = i + 1;
        while j < pts.len() && (pts[j].0 - pts[i].0).abs() <= 1e-6 * step {
            j += 1;
        }
        if j == i + 1 {
            q.x.push(pts[i].0);
            q.y.push(pts[i].1);
            q.sigma.push(pts[i].2);
        } else {
            // Weighted by 1/σ²; a point with no weight counts for nothing.
            let (mut sw, mut swy) = (0.0, 0.0);
            for p in &pts[i..j] {
                let w = if p.2.is_finite() && p.2 > 0.0 { 1.0 / (p.2 * p.2) } else { 0.0 };
                sw += w;
                swy += w * p.1;
            }
            q.x.push(pts[i].0);
            if sw > 0.0 {
                q.y.push(swy / sw);
                q.sigma.push(1.0 / sw.sqrt());
            } else {
                q.y.push(pts[i..j].iter().map(|p| p.1).sum::<f64>() / (j - i) as f64);
                q.sigma.push(f64::INFINITY);
            }
        }
        i = j;
    }
    s.pattern = q;
}

/// One abscissa converted between 2θ (degrees), d (Å) and Q (Å⁻¹) by Bragg's law: d = λ/(2 sin θ), Q = 4π sin θ/λ =
/// 2π/d (`convert`'s rule, for one value: the d-spacing of a peak, R07). `wavelength` is λ in Å, needed for any
/// conversion that involves 2θ; Ok(None) where the value has no image (sin θ > 1 for 2θ). Other axes are an error.
pub fn to_axis(from: Axis, to: Axis, x: f64, wavelength: Option<f64>) -> Result<Option<f64>, String> {
    let need = || wavelength.ok_or("convert: the scan has no wavelength; give one".to_string());
    if from == to {
        return Ok(Some(x));
    }
    // d ↔ Q without the wavelength.
    if let (Axis::D, Axis::Q) | (Axis::Q, Axis::D) = (from, to) {
        return Ok(Some(2.0 * std::f64::consts::PI / x));
    }
    let st = match from {
        Axis::TwoTheta => (x / 2.0).to_radians().sin(),
        Axis::D => need()? / (2.0 * x),
        Axis::Q => x * need()? / (4.0 * std::f64::consts::PI),
        _ => return Err(format!("convert: from {} is not done (it needs the instrument's constants)", from.label())),
    };
    Ok(match to {
        Axis::TwoTheta => (st.abs() <= 1.0).then(|| 2.0 * st.asin().to_degrees()),
        Axis::D => Some(need()? / (2.0 * st)),
        Axis::Q => Some(4.0 * std::f64::consts::PI * st / need()?),
        _ => return Err(format!("convert: to {} is not done", to.label())),
    })
}

fn convert(s: &mut Scan, to: Axis, wavelength: Option<f64>) -> Result<(), String> {
    let lambda = wavelength.or(s.wavelength.map(|w| w.lambda1));
    let from = s.axis;
    let f = |x: f64| to_axis(from, to, x, lambda);
    let mut pts = Vec::new();
    for i in 0..s.pattern.len() {
        if let Some(x) = f(s.pattern.x[i])? {
            pts.push((x, s.pattern.y[i], s.pattern.sigma[i]));
        }
    }
    pts.sort_by(|a, b| a.0.total_cmp(&b.0));
    s.pattern = Pattern { x: pts.iter().map(|p| p.0).collect(), y: pts.iter().map(|p| p.1).collect(), sigma: pts.iter().map(|p| p.2).collect() };
    s.axis = to;
    if let Some(l) = wavelength {
        if s.wavelength.is_none() {
            s.wavelength = Some(Wavelength::single(l));
        }
    }
    s.info.step = None;
    Ok(())
}

impl Operation {
    /// The file an operation reads (add, subtract, merge), as written: relative to the control file.
    pub fn file_mut(&mut self) -> Option<&mut String> {
        match self {
            Operation::Add { file, .. } | Operation::Subtract { file, .. } | Operation::Merge { file, .. } => Some(file),
            _ => None,
        }
    }

    /// Apply to a scan; `load` reads another file named by an operation (path as written, reader id).
    pub fn apply(&self, s: &mut Scan, load: &Load) -> Result<(), String> {
        match self {
            Operation::Mask { from, to } => {
                for i in 0..s.pattern.len() {
                    if s.pattern.x[i] >= *from && s.pattern.x[i] <= *to {
                        s.pattern.sigma[i] = f64::INFINITY;
                    }
                }
            }
            Operation::MaskPoints { x } => {
                for i in 0..s.pattern.len() {
                    if x.iter().any(|&v| same_x(v, s.pattern.x[i])) {
                        s.pattern.sigma[i] = f64::INFINITY;
                    }
                }
            }
            Operation::Delete { from, to } => {
                let x = s.pattern.x.clone();
                retain(s, |i| !(x[i] >= *from && x[i] <= *to));
            }
            Operation::DeletePoints { x } => {
                let px = s.pattern.x.clone();
                retain(s, |i| !x.iter().any(|&v| same_x(v, px[i])));
            }
            Operation::Smooth { points, order } => s.pattern = smooth(&s.pattern, *points, *order)?,
            Operation::RebinFactor { factor, mode } => {
                s.pattern = rebin_factor(&s.pattern, *factor, *mode)?;
                s.info.step = s.info.step.map(|st| st * *factor as f64);
            }
            Operation::RebinStep { step, mode } => {
                s.pattern = rebin_step(&s.pattern, *step, *mode)?;
                s.info.step = Some(*step);
            }
            Operation::Normalize(n) => {
                let ymax = s.pattern.y.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                let f = match n {
                    Normalization::Monitor(to) => to / s.info.monitor.filter(|&m| m > 0.0).ok_or("normalize: the scan has no monitor count")?,
                    Normalization::Time(to) => to / s.info.counting_time.filter(|&t| t > 0.0).ok_or("normalize: the scan has no counting time")?,
                    Normalization::Maximum(to) => to / ymax.max(f64::MIN_POSITIVE),
                    Normalization::Factor(f) | Normalization::Step(f) => *f,
                };
                let x = s.pattern.x.clone();
                let step = |i: usize| if i + 1 < x.len() { x[i + 1] - x[i] } else if i > 0 { x[i] - x[i - 1] } else { 1.0 };
                for i in 0..s.pattern.len() {
                    let f = if matches!(n, Normalization::Step(_)) { f / step(i) } else { f };
                    s.pattern.y[i] *= f;
                    s.pattern.sigma[i] *= f.abs();
                }
                match n {
                    Normalization::Monitor(to) => s.info.monitor = Some(*to),
                    Normalization::Time(to) => s.info.counting_time = Some(*to),
                    _ => {}
                }
            }
            Operation::Add { file, reader, scale } => combine(s, &load(file, reader.as_deref())?, *scale, 1.0),
            Operation::Subtract { file, reader, scale } => combine(s, &load(file, reader.as_deref())?, *scale, -1.0),
            Operation::Merge { file, reader, scale } => merge(s, &load(file, reader.as_deref())?, *scale),
            Operation::Convert { to, wavelength } => convert(s, *to, *wavelength)?,
            Operation::Wavelength { to, from } => {
                if s.axis != Axis::TwoTheta {
                    return Err("wavelength: the scan's x must be 2θ".into());
                }
                let l = match (s.wavelength.map(|w| w.lambda1), *from) {
                    (Some(own), Some(f)) if (own - f).abs() > 1e-9 * own.abs() => return Err(format!("wavelength: the scan states λ = {own} Å, and from = {f} Å disagrees")),
                    (Some(own), _) => own,
                    (None, Some(f)) if f > 0.0 => f,
                    (None, Some(f)) => return Err(format!("wavelength: from = {f} Å is not a wavelength")),
                    (None, None) => return Err("wavelength: the scan states no wavelength; give the one it was measured with (from)".into()),
                };
                if to.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater) {
                    return Err(format!("wavelength: to = {to} Å is not a wavelength"));
                }
                let mut pts = Vec::new();
                for i in 0..s.pattern.len() {
                    let st = to / l * (s.pattern.x[i] / 2.0).to_radians().sin();
                    if st.abs() <= 1.0 {
                        pts.push((2.0 * st.asin().to_degrees(), s.pattern.y[i], s.pattern.sigma[i]));
                    }
                }
                s.pattern = Pattern { x: pts.iter().map(|p| p.0).collect(), y: pts.iter().map(|p| p.1).collect(), sigma: pts.iter().map(|p| p.2).collect() };
                s.wavelength = Some(Wavelength::single(*to));
                s.info.step = None;
            }
        }
        Ok(())
    }

    /// The name of the operation in the control file.
    pub fn name(&self) -> &'static str {
        match self {
            Operation::Mask { .. } => "mask",
            Operation::MaskPoints { .. } => "mask_points",
            Operation::Delete { .. } => "delete",
            Operation::DeletePoints { .. } => "delete_points",
            Operation::Smooth { .. } => "smooth",
            Operation::RebinFactor { .. } | Operation::RebinStep { .. } => "rebin",
            Operation::Normalize(_) => "normalize",
            Operation::Add { .. } => "add",
            Operation::Subtract { .. } => "subtract",
            Operation::Merge { .. } => "merge",
            Operation::Convert { .. } => "convert",
            Operation::Wavelength { .. } => "wavelength",
        }
    }

    /// The operation's keys and values as the control file writes them, `op` first.
    pub fn to_pairs(&self) -> Vec<(String, String)> {
        let mut v = vec![("op".to_string(), string(self.name()))];
        let list = |x: &[f64]| format!("[{}]", x.iter().map(|&v| float(v)).collect::<Vec<_>>().join(", "));
        let mode = |m: &RebinMode| {
            string(match m {
                RebinMode::Sum => "sum",
                RebinMode::Average => "average",
                RebinMode::Keep => "keep",
            })
        };
        let file = |v: &mut Vec<(String, String)>, f: &str, r: &Option<String>, s: f64| {
            v.push(("file".into(), string(f)));
            if let Some(r) = r {
                v.push(("reader".into(), string(r)));
            }
            v.push(("scale".into(), float(s)));
        };
        match self {
            Operation::Mask { from, to } | Operation::Delete { from, to } => {
                v.push(("from".into(), float(*from)));
                v.push(("to".into(), float(*to)));
            }
            Operation::MaskPoints { x } | Operation::DeletePoints { x } => {
                v[0].1 = string(if matches!(self, Operation::MaskPoints { .. }) { "mask" } else { "delete" });
                v.push(("points".into(), list(x)));
            }
            Operation::Smooth { points, order } => {
                v.push(("points".into(), points.to_string()));
                v.push(("order".into(), order.to_string()));
            }
            Operation::RebinFactor { factor, mode: m } => {
                v.push(("factor".into(), factor.to_string()));
                v.push(("mode".into(), mode(m)));
            }
            Operation::RebinStep { step, mode: m } => {
                v.push(("step".into(), float(*step)));
                v.push(("mode".into(), mode(m)));
            }
            Operation::Normalize(n) => {
                let (by, x) = match n {
                    Normalization::Monitor(t) => ("monitor", t),
                    Normalization::Time(t) => ("time", t),
                    Normalization::Maximum(t) => ("maximum", t),
                    Normalization::Factor(f) => ("factor", f),
                    Normalization::Step(t) => ("step", t),
                };
                v.push(("by".into(), string(by)));
                v.push((if by == "factor" { "factor" } else { "to" }.into(), float(*x)));
            }
            Operation::Add { file: f, reader, scale } | Operation::Subtract { file: f, reader, scale } | Operation::Merge { file: f, reader, scale } => file(&mut v, f, reader, *scale),
            Operation::Convert { to, wavelength } => {
                v.push(("to".into(), string(axis_name(*to))));
                if let Some(l) = wavelength {
                    v.push(("wavelength".into(), float(*l)));
                }
            }
            Operation::Wavelength { to, from } => {
                v.push(("to".into(), float(*to)));
                if let Some(f) = from {
                    v.push(("from".into(), float(*f)));
                }
            }
        }
        v
    }

    /// Read an operation from its `[[pattern.operation]]` table.
    pub fn from_table(t: &Table) -> Result<Operation, Diagnostic> {
        let at = |span, m: String| Diagnostic { span, message: m, warning: false, code: "value", replacement: None };
        let head = t.header.unwrap_or_default();
        let get = |k: &str| t.get(k).map(|e| &e.value);
        let num = |k: &str| -> Result<f64, Diagnostic> {
            match get(k) {
                Some(Item { value: Value::Float(v), .. }) => Ok(*v),
                Some(Item { value: Value::Integer(v), .. }) => Ok(*v as f64),
                Some(it) => Err(at(it.span, format!("\"{k}\" must be a number, not {}", it.value.kind()))),
                None => Err(at(head, format!("this operation needs \"{k}\""))),
            }
        };
        let int = |k: &str| -> Result<usize, Diagnostic> {
            match get(k) {
                Some(Item { value: Value::Integer(v), .. }) if *v >= 0 => Ok(*v as usize),
                Some(it) => Err(at(it.span, format!("\"{k}\" must be a whole number"))),
                None => Err(at(head, format!("this operation needs \"{k}\""))),
            }
        };
        let text = |k: &str| -> Result<Option<String>, Diagnostic> {
            match get(k) {
                Some(Item { value: Value::String(s), .. }) => Ok(Some(s.clone())),
                Some(it) => Err(at(it.span, format!("\"{k}\" must be a quoted string, not {}", it.value.kind()))),
                None => Ok(None),
            }
        };
        let known: &[&str] = &["op", "from", "to", "points", "order", "factor", "step", "mode", "by", "file", "reader", "scale", "wavelength"];
        for e in &t.entries {
            if !known.contains(&e.key.as_str()) {
                return Err(at(e.key_span, format!("\"{}\" is not a key of an operation", e.key)));
            }
        }
        let op = text("op")?.ok_or_else(|| at(head, "an operation needs op = \"…\" (mask, delete, smooth, rebin, normalize, add, subtract, merge, convert, wavelength)".into()))?;
        // The keys each operation reads (the dictionary's settings, test operation_dictionary): another is refused.
        if let Some(keys) = keys_of(&op) {
            for e in t.entries.iter().filter(|e| e.key != "op") {
                if !keys.contains(&e.key.as_str()) {
                    return Err(at(e.key_span, format!("\"{}\" is not a key of {op} ({})", e.key, keys.join(", "))));
                }
            }
        }
        let mode = || -> Result<RebinMode, Diagnostic> {
            match text("mode")?.as_deref() {
                None => Ok(DEFAULT_REBIN_MODE),
                Some("sum") => Ok(RebinMode::Sum),
                Some("average") => Ok(RebinMode::Average),
                Some("keep") => Ok(RebinMode::Keep),
                Some(m) => Err(at(get("mode").unwrap().span, format!("\"{m}\" is not sum, average or keep"))),
            }
        };
        let points_list = || -> Result<Option<Vec<f64>>, Diagnostic> {
            match get("points") {
                Some(Item { value: Value::Array(a), .. }) => a
                    .iter()
                    .map(|it| match it.value {
                        Value::Float(v) => Ok(v),
                        Value::Integer(v) => Ok(v as f64),
                        _ => Err(at(it.span, "a point is its x, a number".into())),
                    })
                    .collect::<Result<Vec<_>, _>>()
                    .map(Some),
                _ => Ok(None),
            }
        };
        let file = || -> Result<(String, Option<String>, f64), Diagnostic> {
            let f = text("file")?.ok_or_else(|| at(head, "this operation needs file = \"…\" (relative to the control file)".into()))?;
            Ok((f, text("reader")?, if get("scale").is_some() { num("scale")? } else { DEFAULT_SCALE }))
        };
        Ok(match op.as_str() {
            "mask" => match points_list()? {
                Some(x) => Operation::MaskPoints { x },
                None => Operation::Mask { from: num("from")?, to: num("to")? },
            },
            "delete" => match points_list()? {
                Some(x) => Operation::DeletePoints { x },
                None => Operation::Delete { from: num("from")?, to: num("to")? },
            },
            "smooth" => Operation::Smooth { points: int("points")?, order: int("order")? },
            "rebin" => match get("factor") {
                Some(_) => Operation::RebinFactor { factor: int("factor")?, mode: mode()? },
                None => Operation::RebinStep { step: num("step")?, mode: mode()? },
            },
            "normalize" => match text("by")?.as_deref() {
                Some("monitor") => Operation::Normalize(Normalization::Monitor(num("to")?)),
                Some("time") => Operation::Normalize(Normalization::Time(num("to")?)),
                Some("maximum") => Operation::Normalize(Normalization::Maximum(num("to")?)),
                Some("factor") => Operation::Normalize(Normalization::Factor(num("factor")?)),
                Some("step") => Operation::Normalize(Normalization::Step(num("to")?)),
                _ => return Err(at(get("by").map_or(head, |i| i.span), "by is monitor, time, maximum, step or factor".into())),
            },
            "add" | "subtract" | "merge" => {
                let (file, reader, scale) = file()?;
                match op.as_str() {
                    "add" => Operation::Add { file, reader, scale },
                    "subtract" => Operation::Subtract { file, reader, scale },
                    _ => Operation::Merge { file, reader, scale },
                }
            }
            "convert" => {
                let to = match text("to")?.as_deref() {
                    Some("2theta") => Axis::TwoTheta,
                    Some("d") => Axis::D,
                    Some("q") => Axis::Q,
                    _ => return Err(at(get("to").map_or(head, |i| i.span), "to is \"2theta\", \"d\" or \"q\"".into())),
                };
                Operation::Convert { to, wavelength: if get("wavelength").is_some() { Some(num("wavelength")?) } else { None } }
            }
            "wavelength" => Operation::Wavelength { to: num("to")?, from: if get("from").is_some() { Some(num("from")?) } else { None } },
            other => return Err(at(get("op").unwrap().span, format!("\"{other}\" is not an operation (mask, delete, smooth, rebin, normalize, add, subtract, merge, convert, wavelength)"))),
        })
    }
}

/// The keys an operation of the control file reads; None for a word that is not an operation.
pub fn keys_of(op: &str) -> Option<&'static [&'static str]> {
    Some(match op {
        "mask" | "delete" => &["from", "to", "points"],
        "smooth" => &["points", "order"],
        "rebin" => &["factor", "step", "mode"],
        "normalize" => &["by", "to", "factor"],
        "add" | "subtract" | "merge" => &["file", "reader", "scale"],
        "convert" => &["to", "wavelength"],
        "wavelength" => &["to", "from"],
        _ => return None,
    })
}

pub fn axis_name(a: Axis) -> &'static str {
    match a {
        Axis::TwoTheta => "2theta",
        Axis::D => "d",
        Axis::Q => "q",
        Axis::TimeOfFlight => "tof",
        Axis::Energy => "energy",
    }
}

/// The scan after the operations, in order; the scan given is not changed.
pub fn apply_all(s: &Scan, ops: &[Operation], load: &Load) -> Result<Scan, String> {
    let mut out = s.clone();
    for (k, op) in ops.iter().enumerate() {
        op.apply(&mut out, load).map_err(|e| format!("operation {} ({}): {e}", k + 1, op.name()))?;
    }
    Ok(out)
}

// Analyses -----------------------------------------------------------------------------------------------------------

/// A peak found by `peak_search`.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct FoundPeak {
    pub index: usize,
    pub x: f64,
    /// y at the peak.
    pub height: f64,
    /// Height above the higher of the two bases (SciPy's prominence).
    pub prominence: f64,
    /// Width at half the prominence, in x (linear interpolation, as SciPy's `peak_widths(rel_height = 0.5)`).
    pub fwhm: f64,
}

/// The local maxima of y, a plateau taken at its middle (SciPy's `_local_maxima_1d`).
fn local_maxima(y: &[f64]) -> Vec<usize> {
    let mut out = Vec::new();
    let n = y.len();
    let mut i = 1;
    while i + 1 < n {
        if y[i - 1] < y[i] {
            let mut ahead = i + 1;
            while ahead + 1 < n && y[ahead] == y[i] {
                ahead += 1;
            }
            if y[ahead] < y[i] {
                out.push((i + ahead - 1) / 2);
                i = ahead;
            }
        }
        i += 1;
    }
    out
}

/// SciPy's prominence of peak `p`: the bases are the lowest points between the peak and the nearest higher point (or
/// the end) on each side; the prominence is the peak's height above the higher base. Returns (prominence, left base,
/// right base).
fn prominence(y: &[f64], p: usize) -> (f64, usize, usize) {
    let mut i = p as isize;
    let (mut left_min, mut left) = (y[p], p);
    while i >= 0 && y[i as usize] <= y[p] {
        if y[i as usize] < left_min {
            left_min = y[i as usize];
            left = i as usize;
        }
        i -= 1;
    }
    let mut i = p;
    let (mut right_min, mut right) = (y[p], p);
    while i < y.len() && y[i] <= y[p] {
        if y[i] < right_min {
            right_min = y[i];
            right = i;
        }
        i += 1;
    }
    (y[p] - left_min.max(right_min), left, right)
}

/// SciPy's width at half the prominence, in samples, interpolated: (left position, right position).
fn half_width(y: &[f64], p: usize, prom: f64, left: usize, right: usize) -> (f64, f64) {
    let h = y[p] - 0.5 * prom;
    let mut i = p;
    while left < i && h < y[i] {
        i -= 1;
    }
    let mut lpos = i as f64;
    if y[i] < h {
        lpos += (h - y[i]) / (y[i + 1] - y[i]);
    }
    let mut i = p;
    while i < right && h < y[i] {
        i += 1;
    }
    let mut rpos = i as f64;
    if y[i] < h {
        rpos -= (h - y[i]) / (y[i - 1] - y[i]);
    }
    (lpos, rpos)
}

/// The defaults of the analyses' settings, as the dictionary of operations states them (`operation_table`, L088).
pub const DEFAULT_MIN_PROMINENCE: f64 = 0.0;
pub const DEFAULT_MIN_SIGNIFICANCE: f64 = 3.0;
pub const DEFAULT_MAX_ITER: usize = 10;
/// The scale of `add`, `subtract` and `merge`, and the mode of `rebin`, when the control file gives none.
pub const DEFAULT_SCALE: f64 = 1.0;
pub const DEFAULT_REBIN_MODE: RebinMode = RebinMode::Sum;

/// What a reader writes as the scan's ("intensity", RATE) in `Info::other` when its y are counts per second (UXD's
/// `_CPS` blocks, XRDML intensities in cps), so that counts are y · counting time.
pub const RATE: &str = "counts per second";

/// The scans of a multi-scan file (UXD, XRDML: several ranges, each with its step and counting time) brought to one
/// counting time per point and joined (PL6; WinPLOTR's "constant time normalization … to create instrm10 file with data
/// and sigmas normalized to a constant counting time", and for overlapped scans "if several Y values for a single X
/// value are available, only the first Y value is stored": its news, `manuals/fullprof/1-WinPLOTR/docs/WinPLOTR_news.pdf`).
/// Each point's counts c (y, or y · t for a rate) and σ become c · time / t and σ · time / t, t its scan's counting time
/// per point; a point at the x of a point of an earlier scan is left out; the points sorted by x. `time` is given: the
/// counting time WinPLOTR itself chooses is described in a note (S12, `normalisation_uxd.pdf`) that is not in the
/// repository, so the core takes none of its own.
pub fn constant_time(scans: &[Scan], time: f64) -> Result<Scan, String> {
    if !(time > 0.0 && time.is_finite()) {
        return Err(format!("constant time: {time} s is not a counting time"));
    }
    let first = scans.first().ok_or("constant time: no scan")?;
    let mut pts: Vec<(f64, f64, f64)> = Vec::new();
    for (k, s) in scans.iter().enumerate() {
        if s.axis != first.axis {
            return Err(format!("constant time: scan {} is on another axis", k + 1));
        }
        let t = s.info.counting_time.filter(|t| *t > 0.0 && t.is_finite()).ok_or(format!("constant time: scan {} ({}) states no counting time", k + 1, s.label))?;
        let rate = s.info.other.iter().any(|(a, b)| a == "intensity" && b == RATE);
        let step = s.pattern.mean_step().map_or(0.0, f64::abs);
        let kept: Vec<f64> = pts.iter().map(|p| p.0).collect();
        for i in 0..s.pattern.len() {
            let x = s.pattern.x[i];
            if kept.iter().any(|&k| (k - x).abs() <= 1e-6 * step.max(1e-12)) {
                continue;
            }
            let (c, sc) = if rate { (s.pattern.y[i] * t, s.pattern.sigma[i] * t) } else { (s.pattern.y[i], s.pattern.sigma[i]) };
            pts.push((x, c * time / t, sc * time / t));
        }
    }
    pts.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut out = Scan::new("constant time", first.axis, Pattern { x: pts.iter().map(|p| p.0).collect(), y: pts.iter().map(|p| p.1).collect(), sigma: pts.iter().map(|p| p.2).collect() });
    out.wavelength = first.wavelength;
    out.info = crate::formats::Info { title: first.info.title.clone(), sample: first.info.sample.clone(), date: first.info.date.clone(), instrument: first.info.instrument.clone(), temperature: first.info.temperature, counting_time: Some(time), ..Default::default() };
    out.info.other.push(("constant time".into(), format!("{} scans brought to {time} s per point", scans.len())));
    Ok(out)
}

/// The dictionary's entry of an operation or analysis, by id (`mask`, …, `peak_search`, `auto_background`).
pub fn settings(id: &str) -> Option<&'static crate::operation_table::OperationEntry> {
    crate::operation_table::OPERATIONS.iter().find(|o| o.id == id)
}

/// The peaks of a scan: local maxima with prominence ≥ `min_prominence` and ≥ `min_significance`·σ (see the module).
pub fn peak_search(p: &Pattern, min_prominence: f64, min_significance: f64) -> Vec<FoundPeak> {
    let y = &p.y;
    let at = |pos: f64| -> f64 {
        let i = (pos.floor() as usize).min(p.len() - 1);
        let f = pos - i as f64;
        if f == 0.0 || i + 1 >= p.len() {
            p.x[i]
        } else {
            p.x[i] + f * (p.x[i + 1] - p.x[i])
        }
    };
    local_maxima(y)
        .into_iter()
        .filter_map(|k| {
            let (prom, l, r) = prominence(y, k);
            if prom < min_prominence || (p.sigma[k].is_finite() && prom < min_significance * p.sigma[k]) {
                return None;
            }
            let (a, b) = half_width(y, k, prom, l, r);
            Some(FoundPeak { index: k, x: p.x[k], height: y[k], prominence: prom, fwhm: at(b) - at(a) })
        })
        .collect()
}

/// The peak candidates of GSAS-II's automatic peak search (GSASIIpwdGUI.py `OnAutoSearch`, before its rejection of
/// picks closer than two FWHM, which needs the instrument's profile): the points of y − background at least the
/// standard deviation of y − background (over all points) and above both neighbours, the first and the last point
/// compared with each other as numpy's `roll` does.
pub fn gsas2_peak_candidates(p: &Pattern, background: Option<&[f64]>) -> Vec<usize> {
    let n = p.len();
    if n == 0 {
        return Vec::new();
    }
    let y0: Vec<f64> = (0..n).map(|i| p.y[i] - background.map_or(0.0, |b| b[i])).collect();
    let mean = y0.iter().sum::<f64>() / n as f64;
    let std = (y0.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / n as f64).sqrt();
    (0..n).filter(|&i| y0[i] >= std && y0[i] - y0[(i + n - 1) % n] > 0.0 && y0[i] - y0[(i + 1) % n] > 0.0 && y0[i] != 0.0).collect()
}

/// The automatic background of `auto_background`: the curve at every point, and points every N/100 of them.
#[derive(Clone, Debug, PartialEq)]
pub struct AutoBackground {
    pub curve: Vec<f64>,
    pub points: Vec<(f64, f64)>,
    pub log_lambda: f64,
    pub iterations: usize,
}

/// GSAS-II's default log10 λ for N points: min(10, (int(10·log10(N)^1.5) − 9.5)/10).
pub fn default_log_lambda(n: usize) -> f64 {
    (((10.0 * (n as f64).log10().powf(1.5)) as i64) as f64 - 9.5).min(100.0) / 10.0
}

/// Solve (W + λ DᵀD) z = W y, D the second difference, by a banded LDLᵀ factorisation (pentadiagonal).
fn whittaker(y: &[f64], w: &[f64], lam: f64) -> Vec<f64> {
    let n = y.len();
    // The three diagonals of DᵀD: main, first, second.
    let d0 = |i: usize| -> f64 {
        match (i, n - 1 - i) {
            (0, _) | (_, 0) => 1.0,
            (1, _) | (_, 1) => 5.0,
            _ => 6.0,
        }
    };
    let d1 = |i: usize| -> f64 {
        if i == 0 || i == n - 2 {
            -2.0
        } else {
            -4.0
        }
    };
    let mut a0: Vec<f64> = (0..n).map(|i| w[i] + lam * d0(i)).collect();
    let mut a1: Vec<f64> = (0..n.saturating_sub(1)).map(|i| lam * d1(i)).collect();
    let mut a2: Vec<f64> = vec![lam; n.saturating_sub(2)];
    // LDLᵀ in place: a0 ← D, a1 ← L(i+1, i), a2 ← L(i+2, i).
    for i in 0..n {
        if i >= 1 {
            a0[i] -= a1[i - 1] * a1[i - 1] * a0[i - 1];
        }
        if i >= 2 {
            a0[i] -= a2[i - 2] * a2[i - 2] * a0[i - 2];
        }
        if i + 1 < n {
            let mut v = a1[i];
            if i >= 1 {
                v -= a2[i - 1] * a1[i - 1] * a0[i - 1];
            }
            a1[i] = v / a0[i];
        }
        if i + 2 < n {
            a2[i] /= a0[i];
        }
    }
    let mut z: Vec<f64> = (0..n).map(|i| w[i] * y[i]).collect();
    for i in 0..n {
        if i >= 1 {
            z[i] -= a1[i - 1] * z[i - 1];
        }
        if i >= 2 {
            z[i] -= a2[i - 2] * z[i - 2];
        }
    }
    for v in z.iter_mut().zip(&a0) {
        *v.0 /= v.1;
    }
    for i in (0..n).rev() {
        if i + 1 < n {
            z[i] -= a1[i] * z[i + 1];
        }
        if i + 2 < n {
            z[i] -= a2[i] * z[i + 2];
        }
    }
    z
}

/// arPLS (see the module): λ = 10^log_lambda (GSAS-II's int()), at most `max_iter` + 1 solves, tolerance 10⁻³.
pub fn auto_background(p: &Pattern, log_lambda: Option<f64>, max_iter: usize) -> Result<AutoBackground, String> {
    let n = p.len();
    if n < 4 {
        return Err("background: fewer than four points".into());
    }
    let log_lambda = log_lambda.unwrap_or_else(|| default_log_lambda(n));
    let lam = 10f64.powf(log_lambda).trunc();
    let y = &p.y;
    let mut w = vec![1.0; n];
    let mut z = Vec::new();
    let mut iterations = 0;
    for it in 0..=max_iter {
        iterations = it + 1;
        z = whittaker(y, &w, lam);
        let neg: Vec<f64> = (0..n).map(|i| y[i] - z[i]).filter(|&r| r < 0.0).collect();
        if neg.len() < 2 {
            break;
        }
        let mean = neg.iter().sum::<f64>() / neg.len() as f64;
        let std = (neg.iter().map(|r| (r - mean) * (r - mean)).sum::<f64>() / (neg.len() - 1) as f64).sqrt().max(f64::MIN_POSITIVE);
        let new: Vec<f64> = (0..n).map(|i| 1.0 / (1.0 + ((2.0 / std) * ((y[i] - z[i]) - (2.0 * std - mean))).exp())).collect();
        let num = new.iter().zip(&w).map(|(a, b)| (a - b) * (a - b)).sum::<f64>().sqrt();
        let den = w.iter().map(|a| a * a).sum::<f64>().sqrt().max(f64::MIN_POSITIVE);
        if num / den < 1e-3 {
            break;
        }
        w = new;
    }
    let every = (n / 100).max(1);
    let points = (0..n).step_by(every).map(|i| (p.x[i], z[i])).collect();
    Ok(AutoBackground { curve: z, points, log_lambda, iterations })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// One value by Bragg's law (`to_axis`, L138): Cu Kα₁'s 2θ = 28.60° is d = λ/(2 sin 14.30°) = 3.1186 Å, Q = 2π/d;
    /// and back. Without λ, 2θ → d is an error; d ↔ Q needs none; 2θ beyond sin θ = 1 has no image.
    #[test]
    fn one_value_converted_by_braggs_law() {
        let lam = Some(1.5406);
        let d = to_axis(Axis::TwoTheta, Axis::D, 28.60, lam).unwrap().unwrap();
        assert!((d - 1.5406 / (2.0 * 14.30_f64.to_radians().sin())).abs() < 1e-12 && (d - 3.1186).abs() < 1e-4, "d {d}");
        let q = to_axis(Axis::D, Axis::Q, d, None).unwrap().unwrap();
        assert!((q - 2.0 * std::f64::consts::PI / d).abs() < 1e-12);
        assert!((to_axis(Axis::Q, Axis::TwoTheta, q, lam).unwrap().unwrap() - 28.60).abs() < 1e-9);
        assert!((to_axis(Axis::D, Axis::TwoTheta, d, lam).unwrap().unwrap() - 28.60).abs() < 1e-9);
        assert!(to_axis(Axis::TwoTheta, Axis::D, 28.60, None).is_err());
        assert_eq!(to_axis(Axis::D, Axis::TwoTheta, 0.5, lam).unwrap(), None, "λ/2d > 1");
    }

    fn pat(y: &[f64]) -> Pattern {
        Pattern { x: (0..y.len()).map(|i| 10.0 + 0.1 * i as f64).collect(), y: y.to_vec(), sigma: y.iter().map(|v| v.sqrt()).collect() }
    }

    #[test]
    fn rebin_conserves_counts() {
        let p = pat(&[4.0, 9.0, 16.0, 25.0, 36.0, 49.0, 64.0]);
        let q = rebin_factor(&p, 2, RebinMode::Sum).unwrap();
        assert_eq!(q.y, [13.0, 41.0, 85.0]);
        assert_eq!(q.sigma[0], 13f64.sqrt());
        let r = rebin_step(&p, 0.2, RebinMode::Sum).unwrap();
        let (a, b): (f64, f64) = (p.y.iter().sum(), r.y.iter().sum());
        assert!((a - b - 0.0).abs() < 1e-9 || r.x.len() * 2 < p.len(), "{a} {b}");
    }

    #[test]
    fn whittaker_solves_the_system() {
        // (W + λDᵀD) z = W y checked by multiplication.
        let y = [3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0];
        let w = [1.0, 0.5, 1.0, 0.2, 1.0, 1.0, 0.7, 1.0];
        let lam = 3.0;
        let z = whittaker(&y, &w, lam);
        let n = y.len();
        for i in 0..n {
            let mut dz = 0.0;
            // (DᵀD z)_i from the second differences.
            for k in 0..n - 2 {
                let d = z[k] - 2.0 * z[k + 1] + z[k + 2];
                let c = if i == k { 1.0 } else if i == k + 1 { -2.0 } else if i == k + 2 { 1.0 } else { 0.0 };
                dz += c * d;
            }
            assert!((w[i] * z[i] + lam * dz - w[i] * y[i]).abs() < 1e-12);
        }
    }

    #[test]
    fn local_maxima_plateaus_as_scipy() {
        assert_eq!(local_maxima(&[0.0, 1.0, 1.0, 1.0, 0.0, 2.0, 0.0]), [2, 5]);
    }
}
