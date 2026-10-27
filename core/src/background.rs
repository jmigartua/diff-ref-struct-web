//! Background points and functions (R26, APP_WORKFLOWS W4; L078): the files that hold background points, their
//! interpolation, and a function fitted to them.
//!
//! - Files: our own (`.drb`, a strict subset of TOML as the control file: `format = "diff-ref-struct background 1"`,
//!   `interpolation`, `points = [[2θ, counts], …]`), and FullProf's `.bac` in its two formats (manual p. 8): (1) a first
//!   line "start step end comment" and the intensities after it; (2) a comment line and pairs "x intensity". Both read
//!   and written; FullProf 8.20 reads only the first, with a value at every data point (L079, L081), so a `.bac` is
//!   written in the first format at the data's points.
//! - Interpolation: linear between the points (FullProf's Nba ≥ 2, manual p. 69), or a natural cubic spline through
//!   them (FullProf's cubic splines, Nba < 0 with |Nba| > 4 on p. 69 and Nba ≤ −4 on p. 78; −4 and −5 are the Fourier
//!   and Chebyshev series in the reader, News 17 June 2011). FullProf's spline (Nba = −6, −17, measured in L081 with
//!   `runs/background_runs.py`) is the natural one through the points and one point of its own at the data's last
//!   point with the last value, which the `.pcr` import writes in (`powder`); none when the last point lies past
//!   the data, and the data's end rather than Thmax (L087, `runs/background`: spl6beyond, spl6thmax). Outside the points the background is
//!   the first or last point's value, as FullProf's linear interpolation.
//! - Fitting: a Chebyshev series over [from, to] or FullProf's polynomial Σ B_m (2θ/origin − 1)^m (manual eq. 3.4),
//!   least squares through the points, as starting values for a refined function.

use crate::control::syntax::{self, float, Diagnostic, Value};

/// How the points are joined.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Interpolation {
    Linear,
    Spline,
}

impl Interpolation {
    pub fn name(self) -> &'static str {
        match self {
            Interpolation::Linear => "linear",
            Interpolation::Spline => "spline",
        }
    }

    pub fn from_name(s: &str) -> Option<Interpolation> {
        match s {
            "linear" => Some(Interpolation::Linear),
            "spline" => Some(Interpolation::Spline),
            _ => None,
        }
    }
}

/// The format line of our background file.
pub const FORMAT: &str = "diff-ref-struct background 1";

/// A file of background points.
#[derive(Clone, Debug, PartialEq)]
pub struct PointsFile {
    pub points: Vec<(f64, f64)>,
    /// Our format says it; a `.bac` does not (FullProf takes it from Nba).
    pub interpolation: Option<Interpolation>,
}

/// The second derivatives of the natural cubic spline through (x, y) (zero at both ends).
fn second_derivatives(x: &[f64], y: &[f64]) -> Vec<f64> {
    let n = x.len();
    let mut m = vec![0.0; n];
    if n < 3 {
        return m;
    }
    // Tridiagonal system for the interior points (Thomas algorithm).
    let mut c = vec![0.0; n];
    let mut d = vec![0.0; n];
    for i in 1..n - 1 {
        let (h0, h1) = (x[i] - x[i - 1], x[i + 1] - x[i]);
        let a = h0 / 6.0;
        let b = (h0 + h1) / 3.0;
        let cc = h1 / 6.0;
        let r = (y[i + 1] - y[i]) / h1 - (y[i] - y[i - 1]) / h0;
        let den = b - a * c[i - 1];
        c[i] = cc / den;
        d[i] = (r - a * d[i - 1]) / den;
    }
    for i in (1..n - 1).rev() {
        m[i] = d[i] - c[i] * m[i + 1];
    }
    m
}

/// The natural cubic spline through (x, y) at t; the first or last value outside.
pub fn spline_at(x: &[f64], y: &[f64], m: &[f64], t: f64) -> f64 {
    let n = x.len();
    if n == 0 {
        return 0.0;
    }
    if t <= x[0] {
        return y[0];
    }
    if t >= x[n - 1] {
        return y[n - 1];
    }
    let k = x.partition_point(|&v| v <= t).clamp(1, n - 1);
    let h = x[k] - x[k - 1];
    let a = (x[k] - t) / h;
    let b = (t - x[k - 1]) / h;
    a * y[k - 1] + b * y[k] + ((a * a * a - a) * m[k - 1] + (b * b * b - b) * m[k]) * h * h / 6.0
}

/// The background through the points at each t.
pub fn interpolate(points: &[(f64, f64)], how: Interpolation, t: &[f64]) -> Vec<f64> {
    let x: Vec<f64> = points.iter().map(|p| p.0).collect();
    let y: Vec<f64> = points.iter().map(|p| p.1).collect();
    match how {
        Interpolation::Linear => t.iter().map(|&v| linear_at(points, v)).collect(),
        Interpolation::Spline => {
            let m = second_derivatives(&x, &y);
            t.iter().map(|&v| spline_at(&x, &y, &m, v)).collect()
        }
    }
}

/// Linear interpolation between the points; the first or last value outside.
pub fn linear_at(p: &[(f64, f64)], t: f64) -> f64 {
    if p.is_empty() {
        return 0.0;
    }
    if t <= p[0].0 {
        return p[0].1;
    }
    for w in p.windows(2) {
        if t <= w[1].0 {
            return w[0].1 + (w[1].1 - w[0].1) * (t - w[0].0) / (w[1].0 - w[0].0);
        }
    }
    p[p.len() - 1].1
}

/// d(background at t)/d(y_k) for every point k: the interpolation is linear in the y's.
pub fn weights(points: &[(f64, f64)], how: Interpolation, t: f64) -> Vec<f64> {
    let n = points.len();
    let x: Vec<f64> = points.iter().map(|p| p.0).collect();
    (0..n)
        .map(|k| {
            let e: Vec<f64> = (0..n).map(|i| f64::from(u8::from(i == k))).collect();
            match how {
                Interpolation::Linear => {
                    let p: Vec<(f64, f64)> = x.iter().zip(&e).map(|(a, b)| (*a, *b)).collect();
                    linear_at(&p, t)
                }
                Interpolation::Spline => spline_at(&x, &e, &second_derivatives(&x, &e), t),
            }
        })
        .collect()
}

// Files ---------------------------------------------------------------------------------------------------------------

/// Read a file of background points: our format, or FullProf's `.bac` (either of its two formats).
pub fn read(text: &str) -> Result<PointsFile, String> {
    if text.lines().any(|l| l.trim_start().starts_with("format") && l.contains("diff-ref-struct background")) {
        return read_ours(text).map_err(|d| d.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("; "));
    }
    read_bac(text)
}

fn read_ours(text: &str) -> Result<PointsFile, Vec<Diagnostic>> {
    let doc = syntax::parse(text);
    if doc.has_errors() {
        return Err(doc.diagnostics);
    }
    let root = doc.tables.first().ok_or_else(|| vec![Diagnostic { message: "empty file".into(), ..Default::default() }])?;
    let err = |span, m: &str| vec![Diagnostic { span, message: m.into(), code: "value", ..Default::default() }];
    let mut interpolation = None;
    if let Some(e) = root.get("interpolation") {
        match &e.value.value {
            Value::String(s) => interpolation = Some(Interpolation::from_name(s).ok_or_else(|| err(e.value.span, "interpolation is \"linear\" or \"spline\""))?),
            _ => return Err(err(e.value.span, "interpolation is \"linear\" or \"spline\"")),
        }
    }
    let e = root.get("points").ok_or_else(|| err(root.header.unwrap_or_default(), "the file needs points = [[2θ, counts], …]"))?;
    let Value::Array(rows) = &e.value.value else { return Err(err(e.value.span, "points is a list of [2θ, counts]")) };
    let mut points = Vec::new();
    for r in rows {
        let Value::Array(v) = &r.value else { return Err(err(r.span, "a point is [2θ, counts]")) };
        let nums: Vec<f64> = v.iter().filter_map(|i| match i.value {
            Value::Float(f) => Some(f),
            Value::Integer(n) => Some(n as f64),
            _ => None,
        }).collect();
        if nums.len() != 2 {
            return Err(err(r.span, "a point is [2θ, counts]"));
        }
        points.push((nums[0], nums[1]));
    }
    if points.windows(2).any(|w| w[1].0 <= w[0].0) {
        return Err(err(e.value.span, "the points must be in increasing 2θ"));
    }
    Ok(PointsFile { points, interpolation })
}

fn numbers(l: &str) -> Vec<f64> {
    l.split_whitespace().map_while(|t| t.parse::<f64>().ok()).collect()
}

fn read_bac(text: &str) -> Result<PointsFile, String> {
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let first = lines.next().ok_or("empty .bac file")?;
    let head = numbers(first);
    let rest: Vec<&str> = lines.collect();
    // Format 1 starts with "start step end"; format 2 with a comment (manual p. 8).
    let points = if head.len() >= 3 && !first.trim_start().starts_with('!') {
        // Format 1: start step end, then the intensities.
        let (start, step) = (head[0], head[1]);
        let y: Vec<f64> = rest.iter().flat_map(|l| numbers(l)).collect();
        y.iter().enumerate().map(|(i, &v)| (start + i as f64 * step, v)).collect::<Vec<_>>()
    } else {
        // Format 2: a comment line, then pairs.
        rest.iter().map(|l| numbers(l)).filter(|v| v.len() >= 2).map(|v| (v[0], v[1])).collect()
    };
    if points.is_empty() {
        return Err("no background points in the .bac file".into());
    }
    if points.windows(2).any(|w| w[1].0 <= w[0].0) {
        return Err("the .bac points are not in increasing 2θ".into());
    }
    Ok(PointsFile { points, interpolation: None })
}

/// Our format's text.
pub fn write_ours(points: &[(f64, f64)], how: Interpolation) -> String {
    let mut s = String::from("# diff-ref-struct background points (R26)\n");
    s.push_str(&format!("format = {}\n", syntax::string(FORMAT)));
    s.push_str(&format!("interpolation = {}          # linear or spline\n", syntax::string(how.name())));
    s.push_str("points = [                         # [2θ (deg), counts]\n");
    for (x, y) in points {
        s.push_str(&format!("  [{}, {}],\n", float(*x), float(*y)));
    }
    s.push_str("]\n");
    s
}

/// FullProf's `.bac`, format 1 (manual p. 8): a first line "start step end comment", then the background at every
/// point of the data, ten to a line. FullProf reads a `.bac` this way only (its format 2 is "not available at
/// present", p. 8 footnote) and needs a value at each data point: with fewer it stops at "end of file *.BAC" (L079).
pub fn write_bac(values: &[f64], start: f64, step: f64, title: &str) -> String {
    let end = start + step * (values.len().max(1) - 1) as f64;
    let mut s = format!("{} {} {} ! {title} (written by drs)\n", float(start), float(step), float(end));
    for row in values.chunks(10) {
        s.push_str(&row.iter().map(|v| float(*v)).collect::<Vec<_>>().join(" "));
        s.push('\n');
    }
    s
}

/// The background of `points` joined `how`, written as a `.bac` (format 1) at the points `x` of the data, which must
/// be equally spaced as FullProf's `.bac` is.
pub fn write_bac_at(points: &[(f64, f64)], how: Interpolation, x: &[f64], title: &str) -> Result<String, String> {
    let (start, step) = regular(x).ok_or("a .bac holds the background at equally spaced points: the data's are not")?;
    Ok(write_bac(&interpolate(points, how, x), start, step, title))
}

/// (start, step) when `x` is equally spaced (to 10⁻⁶ of a step), with at least two points.
pub fn regular(x: &[f64]) -> Option<(f64, f64)> {
    if x.len() < 2 {
        return None;
    }
    let step = (x[x.len() - 1] - x[0]) / (x.len() - 1) as f64;
    (step > 0.0 && x.iter().enumerate().all(|(i, &v)| (v - x[0] - i as f64 * step).abs() <= 1e-6 * step)).then_some((x[0], step))
}

/// FullProf's `.bac` format 2 (a comment line, then "2θ intensity" pairs). FullProf 8.20 does not read it (manual
/// p. 8 footnote): it takes the numbers in order as the intensities of format 1 (L081). Written for other programs.
pub fn write_bac_pairs(points: &[(f64, f64)], title: &str) -> String {
    let mut s = format!("! Background points: {title} (written by drs)\n");
    for (x, y) in points {
        s.push_str(&format!("{:>14} {:>14}\n", float(*x), float(*y)));
    }
    s
}

/// Read a points file by its path (our format or a `.bac`).
pub fn read_path(path: &std::path::Path) -> Result<PointsFile, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let text: String = String::from_utf8(bytes.clone()).unwrap_or_else(|_| bytes.iter().map(|&b| b as char).collect());
    read(&text).map_err(|e| format!("{}: {e}", path.display()))
}

/// Write a points file in the format its extension names: `.bac` FullProf's format 1 at the data's points `x` (or at
/// the points themselves when they are equally spaced), otherwise ours.
pub fn write_path(path: &std::path::Path, points: &[(f64, f64)], how: Interpolation, x: Option<&[f64]>) -> Result<(), String> {
    let is_bac = path.extension().is_some_and(|e| e.eq_ignore_ascii_case("bac"));
    let title = path.file_stem().map(|s| s.to_string_lossy().to_string()).unwrap_or_default();
    let text = if is_bac { bac_text(points, how, x, &title)? } else { write_ours(points, how) };
    std::fs::write(path, text).map_err(|e| format!("{}: {e}", path.display()))
}

/// A `.bac` of these points: at the data's points `x` when given, else at the points themselves when they are
/// equally spaced (a background already given at every point); otherwise an error, as FullProf needs a value at every
/// data point.
pub fn bac_text(points: &[(f64, f64)], how: Interpolation, x: Option<&[f64]>, title: &str) -> Result<String, String> {
    match x {
        Some(x) => write_bac_at(points, how, x, title),
        None => {
            let px: Vec<f64> = points.iter().map(|p| p.0).collect();
            let (start, step) = regular(&px).ok_or("a .bac needs the background at every data point: give the data")?;
            Ok(write_bac(&points.iter().map(|p| p.1).collect::<Vec<_>>(), start, step, title))
        }
    }
}

// Fitting -------------------------------------------------------------------------------------------------------------

/// A function to fit to background points.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Function {
    /// T_0 … T_{terms−1} over [from, to].
    Chebyshev { from: f64, to: f64, terms: usize },
    /// Σ B_m (2θ/origin − 1)^m, m < terms (FullProf's eq. 3.4: six terms).
    Polynomial { origin: f64, terms: usize },
}

/// A Chebyshev series over the data's range unless the range is given (L145).
pub fn chebyshev_over(x: &[f64], from: Option<f64>, to: Option<f64>, terms: usize) -> Function {
    Function::Chebyshev { from: from.unwrap_or(x.first().copied().unwrap_or(0.0)), to: to.unwrap_or(x.last().copied().unwrap_or(1.0)), terms }
}

/// A polynomial about the middle of the data's range unless the origin is given (L145).
pub fn polynomial_about(x: &[f64], origin: Option<f64>, terms: usize) -> Function {
    let (lo, hi) = (x.first().copied().unwrap_or(0.0), x.last().copied().unwrap_or(1.0));
    Function::Polynomial { origin: origin.unwrap_or((lo + hi) / 2.0), terms }
}

impl Function {
    fn basis(&self, t: f64) -> Vec<f64> {
        match *self {
            Function::Chebyshev { from, to, terms } => {
                let x = (2.0 * (t - from) / (to - from) - 1.0).clamp(-1.0, 1.0);
                let mut v = Vec::with_capacity(terms);
                for j in 0..terms {
                    v.push(match j {
                        0 => 1.0,
                        1 => x,
                        _ => 2.0 * x * v[j - 1] - v[j - 2],
                    });
                }
                v
            }
            Function::Polynomial { origin, terms } => {
                let u = t / origin - 1.0;
                (0..terms).map(|m| u.powi(m as i32)).collect()
            }
        }
    }

    /// The function with these coefficients at t.
    pub fn at(&self, c: &[f64], t: f64) -> f64 {
        self.basis(t).iter().zip(c).map(|(b, c)| b * c).sum()
    }
}

/// The coefficients that fit the function to the points by least squares, and the rms of the residuals.
pub fn fit(points: &[(f64, f64)], f: Function) -> Result<(Vec<f64>, f64), String> {
    let terms = match f {
        Function::Chebyshev { terms, .. } | Function::Polynomial { terms, .. } => terms,
    };
    if terms == 0 || points.len() < terms {
        return Err(format!("{} points for {terms} terms: at least as many points as terms", points.len()));
    }
    if let Function::Chebyshev { from, to, .. } = f {
        if to <= from {
            return Err("the Chebyshev range needs to > from".into());
        }
    }
    // Least squares by Householder QR of the design matrix (well conditioned for Chebyshev terms).
    let mut a: Vec<Vec<f64>> = points.iter().map(|p| f.basis(p.0)).collect();
    let mut b: Vec<f64> = points.iter().map(|p| p.1).collect();
    let (m, n) = (a.len(), terms);
    for k in 0..n {
        let norm = (k..m).map(|i| a[i][k] * a[i][k]).sum::<f64>().sqrt();
        if norm == 0.0 {
            return Err(format!("term {k} is zero at every point"));
        }
        let alpha = if a[k][k] > 0.0 { -norm } else { norm };
        let mut v: Vec<f64> = (k..m).map(|i| a[i][k]).collect();
        v[0] -= alpha;
        let vv: f64 = v.iter().map(|x| x * x).sum();
        if vv == 0.0 {
            continue;
        }
        for j in k..n {
            let s: f64 = a[k..].iter().zip(&v).map(|(row, vi)| vi * row[j]).sum::<f64>() * 2.0 / vv;
            for (row, vi) in a[k..].iter_mut().zip(&v) {
                row[j] -= s * vi;
            }
        }
        let s: f64 = (k..m).map(|i| v[i - k] * b[i]).sum::<f64>() * 2.0 / vv;
        for i in k..m {
            b[i] -= s * v[i - k];
        }
    }
    let mut c = vec![0.0; n];
    for k in (0..n).rev() {
        c[k] = (b[k] - (k + 1..n).map(|j| a[k][j] * c[j]).sum::<f64>()) / a[k][k];
    }
    let rms = (points.iter().map(|p| (p.1 - f.at(&c, p.0)).powi(2)).sum::<f64>() / points.len() as f64).sqrt();
    Ok((c, rms))
}

#[cfg(test)]
mod default_tests {
    /// L145: the function's range or origin from the data when not given.
    #[test]
    fn a_functions_defaults_from_the_data() {
        let x = [10.0, 20.0, 50.0];
        assert_eq!(super::chebyshev_over(&x, None, Some(40.0), 4), super::Function::Chebyshev { from: 10.0, to: 40.0, terms: 4 });
        assert_eq!(super::polynomial_about(&x, None, 6), super::Function::Polynomial { origin: 30.0, terms: 6 });
        assert_eq!(super::polynomial_about(&x, Some(40.0), 6), super::Function::Polynomial { origin: 40.0, terms: 6 });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spline_through_a_cubic_and_a_line() {
        // A natural spline reproduces a straight line exactly.
        let p: Vec<(f64, f64)> = (0..8).map(|i| (10.0 + 3.0 * i as f64, 5.0 - 0.5 * i as f64)).collect();
        for t in [10.0, 11.3, 20.7, 31.0] {
            let v = interpolate(&p, Interpolation::Spline, &[t])[0];
            assert!((v - (5.0 - 0.5 * (t - 10.0) / 3.0)).abs() < 1e-12);
        }
        // It passes through every point; the weights add up to 1.
        let q = [(10.0, 3.0), (15.0, 9.0), (21.0, 4.0), (30.0, 7.0)];
        for (x, y) in q {
            assert!((interpolate(&q, Interpolation::Spline, &[x])[0] - y).abs() < 1e-12);
        }
        let w = weights(&q, Interpolation::Spline, 17.5);
        assert!((w.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        let direct = interpolate(&q, Interpolation::Spline, &[17.5])[0];
        assert!((w.iter().zip(&q).map(|(w, p)| w * p.1).sum::<f64>() - direct).abs() < 1e-12);
    }

    #[test]
    fn bac_formats_read() {
        let f1 = "   10.0  0.5  11.0  Background of: x\n 1 2 3\n\n";
        assert_eq!(read(f1).unwrap().points, [(10.0, 1.0), (10.5, 2.0), (11.0, 3.0)]);
        let f2 = "! Difference pattern for : x\n 10.0  4.0\n 10.5  5.0\n";
        assert_eq!(read(f2).unwrap().points, [(10.0, 4.0), (10.5, 5.0)]);
        let ours = write_ours(&[(10.0, 4.0), (20.0, 5.5)], Interpolation::Spline);
        let r = read(&ours).unwrap();
        assert_eq!((r.points, r.interpolation), (vec![(10.0, 4.0), (20.0, 5.5)], Some(Interpolation::Spline)));
        assert_eq!(read(&write_bac(&[4.0, 5.5, 6.0], 10.0, 0.5, "t")).unwrap().points, [(10.0, 4.0), (10.5, 5.5), (11.0, 6.0)]);
        assert_eq!(read(&write_bac_pairs(&[(10.0, 4.0), (20.0, 5.5)], "t")).unwrap().points, [(10.0, 4.0), (20.0, 5.5)]);
        assert!(bac_text(&[(10.0, 4.0), (20.0, 5.5), (25.0, 1.0)], Interpolation::Linear, None, "t").is_err());
    }

    #[test]
    fn fit_recovers_a_series() {
        let f = Function::Chebyshev { from: 10.0, to: 100.0, terms: 5 };
        let c = [50.0, -12.0, 3.0, 0.5, -0.25];
        let p: Vec<(f64, f64)> = (0..30).map(|i| (10.0 + 3.0 * i as f64, 0.0)).map(|(t, _)| (t, f.at(&c, t))).collect();
        let (got, rms) = fit(&p, f).unwrap();
        assert!(rms < 1e-10 && got.iter().zip(&c).all(|(a, b)| (a - b).abs() < 1e-9), "{got:?}");
        let g = Function::Polynomial { origin: 60.0, terms: 4 };
        let d = [20.0, 4.0, -1.0, 0.5];
        let p: Vec<(f64, f64)> = (0..20).map(|i| 15.0 + 7.0 * i as f64).map(|t| (t, g.at(&d, t))).collect();
        let (got, _) = fit(&p, g).unwrap();
        assert!(got.iter().zip(&d).all(|(a, b)| (a - b).abs() < 1e-9), "{got:?}");
    }
}
