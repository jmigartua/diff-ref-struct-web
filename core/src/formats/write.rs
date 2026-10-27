//! Writing a scan to a file the registry reads back (app request 12; L088), and the line of each point of a file
//! (app request 15; L088).
//!
//! Formats written:
//! - `xys` (FullProf's Ins = 10, manual p. 136): `XYDATA`, a `TITLE` line and four comment lines, then x, y, σ; the
//!   wavelength in a comment (FullProf takes it from the `.pcr`).
//! - `xye` (Topas and GSAS-II's `G2pwd_xye.py`): '#' comment lines with the title and the wavelength, then x, y, σ;
//!   x in 2θ, d or Q, the quantity named in the header line ("x = d (Angstrom)"), which the core's reader takes back.
//! - `qye` (GSAS-II's Q, y, σ): as xye, x in Q.
//!
//! Numbers are written to the shortest text that reads back to the same value. A point with no weight (σ = ∞, a masked
//! point) has no form in either: it is written with σ = 10³⁰, which both readers take as σ and the core reads back as
//! no weight (`NO_WEIGHT`). xys holds x as 2θ only. Every value reads back to the bit,
//! but for one rule of the xye reader (GSAS-II's): a point with y ≤ 0 is read as y = 0 with no weight; xys keeps it.

use super::{Axis, Scan};
use crate::control::syntax::float;

/// The σ written for a point with no weight; a σ this large or larger is read back as no weight.
pub const NO_WEIGHT: f64 = 1e30;

/// The formats `write` writes, by the registry's ids.
pub const WRITTEN: [&str; 3] = ["xys", "xye", "qye"];

/// The words for x in the header line, which the core's xye reader takes back as the scan's axis (PL7).
pub fn x_words(a: Axis) -> Option<&'static str> {
    match a {
        Axis::TwoTheta => Some("2theta (deg)"),
        Axis::D => Some("d (Angstrom)"),
        Axis::Q => Some("Q (1/Angstrom)"),
        _ => None,
    }
}

/// A scan as the text of a file in format `id` (one of `WRITTEN`): xys in 2θ only (FullProf takes an Ins = 10 file's
/// x as the pattern's, and its header has no word for d or Q: INTER and TEMP only, manual p. 136); xye in 2θ, d or Q
/// (PL7), the quantity in its header line; qye (GSAS-II's Q, x, y, σ) in Q only.
pub fn write(s: &Scan, id: &str) -> Result<String, String> {
    let ok = match id {
        "xys" => s.axis == Axis::TwoTheta,
        "xye" => matches!(s.axis, Axis::TwoTheta | Axis::D | Axis::Q),
        "qye" => s.axis == Axis::Q,
        _ => true,
    };
    if !ok {
        return Err(format!("{id} does not hold x as {} (xys: 2θ; xye: 2θ, d or Q; qye: Q)", s.axis.label()));
    }
    let xw = x_words(s.axis).unwrap_or("2theta (deg)");
    let title = s.info.title.clone().unwrap_or_else(|| s.label.clone());
    let lambda = s.wavelength.map(|w| format!("wavelength {}{}{}", float(w.lambda1), w.lambda2.map_or(String::new(), |l| format!(" {}", float(l))), w.ratio.map_or(String::new(), |r| format!(" ratio {}", float(r)))));
    let mut t = match id {
        "xys" => format!(
            "XYDATA\nTITLE: {title} (written by drs {})\n! {}\n! x = 2theta (deg), y, sigma; sigma {} marks a point with no weight\n!\n!\n",
            env!("CARGO_PKG_VERSION"),
            lambda.as_deref().unwrap_or("wavelength not stated"),
            float(NO_WEIGHT)
        ),
        "xye" | "qye" => format!(
            "# {title} (written by drs {})\n# {}\n# x = {xw}, y, sigma; sigma {} marks a point with no weight\n",
            env!("CARGO_PKG_VERSION"),
            lambda.as_deref().unwrap_or("wavelength not stated"),
            float(NO_WEIGHT)
        ),
        other => return Err(format!("{other}: not a format the core writes ({})", WRITTEN.join(", "))),
    };
    let p = &s.pattern;
    for i in 0..p.len() {
        let sigma = if p.sigma[i].is_finite() { float(p.sigma[i]) } else { float(NO_WEIGHT) };
        t.push_str(&format!("{} {} {}\n", float(p.x[i]), float(p.y[i]), sigma));
    }
    Ok(t)
}

/// The format `write` takes for a file name: its extension when it is one of `WRITTEN`.
pub fn format_of(name: &str) -> Option<&'static str> {
    let ext = name.rsplit_once('.').map(|(_, e)| e.to_lowercase())?;
    WRITTEN.iter().copied().find(|w| *w == ext)
}

// The line of each point ---------------------------------------------------------------------------------------------

fn numbers(l: &str) -> Vec<f64> {
    l.split(|c: char| c.is_whitespace() || c == ',' || c == ';').filter(|t| !t.is_empty()).map_while(|t| t.parse::<f64>().ok()).collect()
}

fn same(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-9 * a.abs().max(b.abs()).max(1.0)
}

/// The line (from 1) each point of the scan was read from, when it can be established from the text; None otherwise
/// (a binary file, a format whose values are transformed on reading, a file not laid out either way below):
/// - one point per row: each point on the next line that starts with two or more numbers, whose first is the point's x
///   and which holds its y (column files: xy, xye, xys, csv, …);
/// - steps (FullProf's free format and the like): the line holding Ti, step and Tf, then every number after it, in
///   order, one point's y each, no number skipped.
pub fn point_lines(text: &str, s: &Scan) -> Option<Vec<usize>> {
    let p = &s.pattern;
    if p.is_empty() {
        return Some(Vec::new());
    }
    let ls: Vec<&str> = text.lines().collect();
    rows(&ls, s).or_else(|| steps(&ls, s))
}

fn rows(ls: &[&str], s: &Scan) -> Option<Vec<usize>> {
    let p = &s.pattern;
    let mut out = Vec::with_capacity(p.len());
    let mut k = 0;
    // The first row: the first line whose numbers start with the first point's x and hold its y.
    while k < ls.len() {
        let v = numbers(ls[k]);
        if v.len() >= 2 && same(v[0], p.x[0]) && v[1..].iter().any(|&y| same(y, p.y[0])) {
            break;
        }
        k += 1;
    }
    for i in 0..p.len() {
        // The next line with two numbers or more (others, comments and blanks, skipped).
        while k < ls.len() && numbers(ls[k]).len() < 2 {
            k += 1;
        }
        let v = numbers(ls.get(k)?);
        if !(same(v[0], p.x[i]) && v[1..].iter().any(|&y| same(y, p.y[i]))) {
            return None;
        }
        out.push(k + 1);
        k += 1;
    }
    Some(out)
}

fn steps(ls: &[&str], s: &Scan) -> Option<Vec<usize>> {
    let p = &s.pattern;
    let n = p.len();
    let step = if n > 1 { (p.x[n - 1] - p.x[0]) / (n - 1) as f64 } else { 0.0 };
    // The header: the first line whose first three numbers are Ti, the step and a Tf at or past the last point.
    let k = ls.iter().position(|l| {
        let v = numbers(l);
        v.len() >= 3 && same(v[0], p.x[0]) && (n < 2 || (v[1] - step).abs() <= 1e-6 * step.abs().max(1e-12) + 1e-9) && v[2] >= p.x[n - 1] - step.abs() * 0.5
    })?;
    let mut out = Vec::with_capacity(n);
    for (j, l) in ls.iter().enumerate().skip(k + 1) {
        for v in numbers(l) {
            if out.len() == n {
                break;
            }
            if !same(v, p.y[out.len()]) {
                return None;
            }
            out.push(j + 1);
        }
        if out.len() == n {
            break;
        }
    }
    (out.len() == n).then_some(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::Pattern;

    fn scan() -> Scan {
        Scan::new("t", Axis::TwoTheta, Pattern { x: vec![10.0, 10.02, 10.04], y: vec![5.0, 7.5, 4.0], sigma: vec![2.0, f64::INFINITY, 2.0] })
    }

    #[test]
    fn written_and_read_back() {
        for id in WRITTEN {
            // qye holds Q only (PL7): its scan on Q, the same numbers.
            let s = if id == "qye" { Scan::new("t", Axis::Q, scan().pattern) } else { scan() };
            let t = write(&s, id).unwrap();
            let d = crate::formats::read(&crate::formats::Probe::new(&format!("x.{id}"), t.as_bytes()), Some(id)).unwrap();
            let p = &d.scans[0].pattern;
            assert_eq!((p.x.clone(), p.y.clone()), (scan().pattern.x, scan().pattern.y), "{id}");
            assert_eq!((p.sigma[0], p.sigma[2]), (2.0, 2.0), "{id}");
            assert!(p.sigma[1].is_infinite(), "{id}: no weight");
            assert_eq!(point_lines(&t, &d.scans[0]), Some(vec![t.lines().count() - 2, t.lines().count() - 1, t.lines().count()]), "{id}");
        }
        // xys holds 2θ only and qye Q only (PL7, L195); xye holds 2θ, d or Q.
        assert!(write(&Scan::new("d", Axis::D, scan().pattern), "xys").is_err());
        assert!(write(&scan(), "qye").is_err());
        assert!(write(&Scan::new("d", Axis::D, scan().pattern), "xye").is_ok());
    }

    #[test]
    fn lines_of_steps() {
        let t = "title\n10.0 0.02 10.04\n5 7.5\n4\n";
        let s = Scan::new("t", Axis::TwoTheta, Pattern { x: vec![10.0, 10.02, 10.04], y: vec![5.0, 7.5, 4.0], sigma: vec![1.0; 3] });
        assert_eq!(point_lines(t, &s), Some(vec![3, 3, 4]));
        // A number between that is not the next point's: not established.
        assert_eq!(point_lines("10.0 0.02 10.04\n5 9 7.5 4\n", &s), None);
    }
}
