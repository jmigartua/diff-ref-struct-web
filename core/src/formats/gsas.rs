//! GSAS standard powder data files (`.gsas`, `.gss`, `.gsa`, `.gda`, `.fxye`, `.raw`, `.xra`, `.cwn`; FullProf's
//! Ins = 12): a title line, an optional "Instrument parameter file:" line, '#' comments, and banks, each opened by a
//! BANK record `BANK n NCHAN NREC BINTYP BCOEF1 BCOEF2 ... TYPE`, read as GSAS-II's `G2pwd_fxye.py` reads them:
//! - STD (the default): ten (I2, F6) fields per line, counters and counts; σ² = y/counters; y ≤ 0 has no weight.
//! - ESD: five (F8, F8) fields per line, y and σ; y ≤ 0 has no weight.
//! - FXY, FXYE: x, y [, σ] per line, x in centidegrees (CW) or µs (TOF); y ≤ 0 (or σ ≤ 0) has no weight; FXY σ² = y.
//! - ALT: four (F8, F7, F6) fields per line, x·32, y·1000, σ·1000.
//! - TIME_MAP records give the TOF of each channel (clock width in ns) for the banks of BINTYP TIME_MAP.
//! - EDS in the BANK record: energy-dispersive, E = A + B·c + C·c² + D·c³ keV for channel c.
//!
//! x is 2θ in degrees for constant wavelength (BINTYP CONST or CONS), the TOF in µs for RALF, SLOG and TIME_MAP banks
//! (GSAS-II's reader divides these by 100 too and multiplies back once instrument parameters say TOF), keV for EDS.
//! A point with no weight has σ = ∞.
//!
//! Two departures from GSAS-II's reading, both of its loops: a STD bank keeps NCHAN points (GSAS-II's `j < Nch` keeps
//! NCHAN − 1), and the empty fields of the last line of an ESD bank are not points (GSAS-II's test for eight blanks
//! does not see the end of a short line, and adds zero-weight points there).

use super::fullprof::fnum;
use super::{lines, Axis, Dataset, Info, Probe, Scan, SIGNATURE};
use crate::data::{DataError, Pattern};

fn err<T>(m: impl Into<String>) -> Result<T, DataError> {
    Err(DataError(m.into()))
}

/// Python's `S[a:b]` of a line.
fn sub(s: &str, a: usize, b: usize) -> &str {
    let n = s.len();
    s.get(a.min(n)..b.min(n)).unwrap_or("")
}

/// GSAS-II's sfloat: blank is 0.
fn sfloat(s: &str) -> Result<f64, DataError> {
    let t = s.trim();
    if t.is_empty() {
        Ok(0.0)
    } else {
        fnum(t).ok_or(DataError(format!("GSAS: \"{t}\" is not a number")))
    }
}

/// A TIME_MAP: rows (first channel, TOF in clock ticks, tick step), and the clock width in µs.
struct TimeMap {
    rows: Vec<[i64; 3]>,
    clock: f64,
}

impl TimeMap {
    /// GSAS-II's GetTimeMap and Tmap2TOF: the TOF of every channel.
    fn tof(&self) -> Vec<f64> {
        let mut m = self.rows.clone();
        let n = m.len();
        if n < 2 {
            return Vec::new();
        }
        let tmax = m[n - 1][0];
        let nch = m[n - 2][0] + (tmax - m[n - 2][1] + m[n - 2][2] - 1) / m[n - 2][2].max(1);
        m[n - 1] = [nch + 1, tmax, 0];
        for r in &mut m {
            r[0] -= 1;
        }
        let mut out = Vec::new();
        let (mut tch, mut t, mut step) = (m[0][0], m[0][1], m[0][2]);
        for r in &m[1..] {
            for i in tch..r[0] {
                out.push((t + step * (i - tch)) as f64 * self.clock);
            }
            (tch, t, step) = (r[0], r[1], r[2]);
        }
        out
    }
}

struct Bank<'a> {
    record: &'a str,
    body: Vec<&'a str>,
    comments: Vec<String>,
}

pub(crate) fn parse(text: &str) -> Result<Dataset, DataError> {
    let ls = lines(text);
    let title = ls.first().map(|l| l.trim().to_string()).unwrap_or_default();
    let mut banks: Vec<Bank> = Vec::new();
    let mut maps: Vec<(String, TimeMap)> = Vec::new();
    let mut comments: Vec<String> = Vec::new();
    let mut instprm = None;
    let mut i = 1;
    while i < ls.len() {
        let s = ls[i];
        if i == 1 && s.to_lowercase().starts_with("inst") && s.contains(':') {
            instprm = s.split_once(':').map(|(_, v)| v.trim().trim_matches(['[', ']']).trim().to_string());
            i += 1;
            continue;
        }
        if let Some(c) = s.strip_prefix('#') {
            comments.push(c.trim().to_string());
            i += 1;
            continue;
        }
        if let Some(rest) = s.strip_prefix("TIME_MAP") {
            let cons: Vec<&str> = rest.split_whitespace().collect();
            let nrec: usize = cons.get(2).and_then(|v| v.parse().ok()).ok_or(DataError("GSAS: TIME_MAP record".into()))?;
            let clock = cons.get(4).and_then(|v| fnum(v)).ok_or(DataError("GSAS: TIME_MAP clock width".into()))? / 1000.0;
            let no = if s.contains("TIME_MAP10") { "1".to_string() } else { cons[0].to_string() };
            let vals: Vec<i64> = ls[i + 1..(i + 1 + nrec).min(ls.len())].iter().flat_map(|l| l.split_whitespace().filter_map(|v| v.parse().ok())).collect();
            let rows = vals.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
            maps.push((no, TimeMap { rows, clock }));
            i += 1 + nrec;
            continue;
        }
        if s.starts_with("BANK") {
            let mut body = Vec::new();
            i += 1;
            while i < ls.len() && !ls[i].starts_with("BANK") && !ls[i].starts_with('#') && !ls[i].starts_with("TIME_MAP") {
                body.push(ls[i]);
                i += 1;
            }
            banks.push(Bank { record: s, body, comments: std::mem::take(&mut comments) });
            continue;
        }
        i += 1;
    }
    if banks.is_empty() {
        return err("GSAS: no BANK record");
    }
    let mut scans = Vec::new();
    for b in &banks {
        let cons: Vec<&str> = b.record.split_whitespace().collect();
        let bintyp = cons.get(4).copied().unwrap_or("");
        let kind = ["FXYE", "FXY", "ESD", "STD", "ALT"].into_iter().find(|t| b.record.contains(t)).unwrap_or("STD");
        let eds = b.record.contains("EDS");
        let tof = matches!(bintyp, "RALF" | "SLOG" | "TIME_MAP");
        let nch: usize = cons.get(2).and_then(|v| v.parse().ok()).unwrap_or(usize::MAX);
        let map = if bintyp == "TIME_MAP" || !maps.is_empty() {
            let key = cons.get(5).copied().unwrap_or("1");
            maps.iter().find(|(k, _)| k == key).or(maps.first()).map(|(_, m)| m)
        } else {
            None
        };
        // The scale from the file's x to the scan's: centidegrees for CW, µs as written for TOF.
        let xs = if tof { 1.0 } else { 0.01 };
        let (start, step) = (cons.get(5).map_or(Ok(0.0), |v| sfloat(v))? * xs, cons.get(6).map_or(Ok(0.0), |v| sfloat(v))? * xs);
        let ecoef: Vec<f64> = if eds { cons.iter().skip(5).take(4).filter_map(|v| fnum(v)).collect() } else { Vec::new() };
        let mut p = Pattern::default();
        let mut push = |x: f64, y: f64, var_or_sigma: Result<f64, f64>| {
            p.x.push(x);
            match var_or_sigma {
                Ok(sigma) if y > 0.0 && sigma > 0.0 => {
                    p.y.push(y);
                    p.sigma.push(sigma);
                }
                Err(var) if y > 0.0 && var > 0.0 => {
                    p.y.push(y);
                    p.sigma.push(var.sqrt());
                }
                _ => {
                    p.y.push(0.0);
                    p.sigma.push(f64::INFINITY);
                }
            }
        };
        match kind {
            "FXYE" | "FXY" => {
                for l in &b.body {
                    let v: Vec<&str> = l.split_whitespace().collect();
                    if v.len() < 2 {
                        continue;
                    }
                    let (x, y) = (sfloat(v[0])? * xs, sfloat(v[1])?);
                    if kind == "FXYE" {
                        push(x, y, Ok(sfloat(v.get(2).copied().unwrap_or(""))?));
                    } else {
                        push(x, y, Err(y));
                    }
                }
            }
            "ESD" => {
                let mut j = 0;
                'lines: for l in &b.body {
                    for i in (0..80).step_by(16) {
                        if sub(l, i, i + 8).trim().is_empty() || j >= nch {
                            continue 'lines;
                        }
                        push(start + step * j as f64, sfloat(sub(l, i, i + 8))?, Ok(sfloat(sub(l, i + 8, i + 16))?));
                        j += 1;
                    }
                }
            }
            "ALT" => {
                for l in &b.body {
                    for i in (0..80).step_by(20) {
                        let xi = sfloat(sub(l, i, i + 8))?;
                        if xi == 0.0 {
                            break;
                        }
                        let x = if tof { xi / 32.0 } else { xi / 3200.0 };
                        push(x, sfloat(sub(l, i + 8, i + 15))? / 1000.0, Ok(sfloat(sub(l, i + 15, i + 21))? / 1000.0));
                    }
                }
            }
            _ => {
                let mut j = 0usize;
                'lines: for l in &b.body {
                    if l.contains('\u{1a}') {
                        break;
                    }
                    for i in (0..80).step_by(8) {
                        if sub(l, i, i + 8).trim().is_empty() || j >= nch {
                            continue 'lines;
                        }
                        let ni = sfloat(sub(l, i, i + 2))?.max(1.0);
                        let yi = sfloat(sub(l, i + 2, i + 8))?.max(0.0);
                        let x = if eds {
                            let c = j as f64;
                            ecoef.first().unwrap_or(&0.0) + (ecoef.get(1).unwrap_or(&0.0) + ecoef.get(2).unwrap_or(&0.0) * c) * c + ecoef.get(3).unwrap_or(&0.0) * c * c * c
                        } else {
                            start + step * j as f64
                        };
                        push(x, yi, Err(yi / ni));
                        j += 1;
                    }
                }
            }
        }
        if let Some(m) = map {
            let t = m.tof();
            for (x, tx) in p.x.iter_mut().zip(t) {
                *x = tx;
            }
        }
        let axis = if eds {
            Axis::Energy
        } else if tof {
            Axis::TimeOfFlight
        } else {
            Axis::TwoTheta
        };
        let mut s = Scan::new(cons.iter().take(2).copied().collect::<Vec<_>>().join(" "), axis, p);
        s.info = Info { title: Some(title.clone()), ..Default::default() };
        if !tof && !eds && matches!(kind, "STD" | "ESD") {
            s.info.step = Some(step);
        }
        if let Some(f) = &instprm {
            s.info.other.push(("instrument parameter file".into(), f.clone()));
        }
        for c in &b.comments {
            let (k, v) = c.split_once('=').unwrap_or((c, ""));
            let num = v.trim().parse::<f64>().ok();
            if k.contains("Temp") {
                s.info.temperature = num.or(s.info.temperature);
            } else if k.to_lowercase().contains("instrument") {
                s.info.instrument = Some(v.trim().to_string()).filter(|v| !v.is_empty()).or_else(|| c.split_once(':').map(|(_, v)| v.trim().to_string()));
            } else if k.to_lowercase().contains("monitor") {
                s.info.monitor = num;
            }
            s.info.other.push(("comment".into(), c.clone()));
        }
        scans.push(s);
    }
    Ok(Dataset { format: "gsas", scans })
}

pub fn recognise(p: &Probe) -> f64 {
    // No test for text: GSAS files may hold stray bytes in their title (LaMnO3_50k.gsas has a NUL); a BANK record
    // as GSAS writes it: BANK as GSAS writes it: BANK, the bank number, NCHAN and NREC, and a BINTYP.
    let record = p.text().lines().any(|l| {
        let c: Vec<&str> = l.split_whitespace().collect();
        l.starts_with("BANK") && c.len() >= 5 && c[1].parse::<u32>().is_ok() && c[2].parse::<u64>().is_ok() && matches!(c[4], "CONST" | "CONS" | "RALF" | "SLOG" | "TIME_MAP" | "LPSD" | "LOG6")
    });
    match parse(p.text()) {
        Ok(d) if record && d.scans.iter().any(|s| s.pattern.len() >= 4) => SIGNATURE,
        _ => 0.0,
    }
}

pub fn read(p: &Probe) -> Result<Dataset, DataError> {
    parse(p.text())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn std_esd_fxye() {
        let t = "title\nBANK 1 3 1 CONST 1000 5 0 0 STD\n 1   100 2   200 1     0\nBANK 2 3 1 CONST 1000 5 0 0 ESD\n   100.0    10.0   200.0    20.0\nBANK 3 2 2 CONS 1000 5 0 0 FXYE\n 1000 4.0 2.0\n 1005 9.0 0\n";
        let d = parse(t).unwrap();
        assert_eq!(d.scans.len(), 3);
        let s = &d.scans[0].pattern;
        assert_eq!(s.x, [10.0, 10.05, 10.1]);
        assert_eq!(s.y, [100.0, 200.0, 0.0]);
        assert!((s.sigma[1] - 10.0).abs() < 1e-12 && s.sigma[2].is_infinite());
        assert_eq!(d.scans[1].pattern.sigma, [10.0, 20.0]);
        assert_eq!(d.scans[2].pattern.sigma, [2.0, f64::INFINITY]);
    }
}
