//! Rigaku's formats, read as GSAS-II's `G2pwd_rigaku.py` reads them.
//! - RAS: a first line `*RAS_DATA_START`; each `*RAS_HEADER_START … *RAS_INT_START … *RAS_INT_END` is a scan, its
//!   lines x, y (and an attenuator factor, not used); σ² = max(y, 1). The header's `*HW_XG_WAVE_LENGTH_ALPHA1`,
//!   `…_ALPHA2`, `*MEAS_SCAN_STEP`, `*FILE_SAMPLE`, `*MEAS_SCAN_START_TIME` are kept (GSAS-II does not read them).
//! - RASX: a zip holding `Data0/Profile0.txt` (x, y per line after a three-byte mark) and
//!   `Data0/MesurementConditions0.xml` (sic); σ² = max(y, 1). The wavelength is read from the XML's
//!   `WavelengthKalpha1` when present.
//! - TXT: comment lines (any line not wholly numbers), then lines of an angle and a fixed number k ≥ 1 of intensities
//!   at a constant step: the line's j-th intensity is at angle + j·step, the step the difference of two lines' angles
//!   over k; σ² = max(y, 1). GSAS-II requires more than 30 such lines.

use super::xml;
use super::zip::Archive;
use super::{all_numbers, lines, Axis, Dataset, Probe, Scan, Wavelength, LOOSE, SIGNATURE, STRUCTURE};
use crate::data::{DataError, Pattern};

fn err<T>(m: impl Into<String>) -> Result<T, DataError> {
    Err(DataError(m.into()))
}

fn pattern(pairs: &[(f64, f64)]) -> Pattern {
    Pattern { x: pairs.iter().map(|p| p.0).collect(), y: pairs.iter().map(|p| p.1).collect(), sigma: pairs.iter().map(|p| p.1.max(1.0).sqrt()).collect() }
}

fn parse_ras(p: &Probe) -> Result<Dataset, DataError> {
    let ls = lines(p.text());
    if ls.first().map(|l| l.trim_end()) != Some("*RAS_DATA_START") {
        return err("RAS: the first line is not *RAS_DATA_START");
    }
    let mut scans = Vec::new();
    let mut header: Vec<(String, String)> = Vec::new();
    let mut i = 1;
    while i < ls.len() {
        let l = ls[i].trim_end();
        i += 1;
        if l == "*RAS_HEADER_START" {
            header.clear();
            continue;
        }
        if let Some(rest) = l.strip_prefix('*') {
            if let Some((k, v)) = rest.split_once(' ') {
                header.push((k.to_string(), v.trim().trim_matches('"').to_string()));
            }
        }
        if l != "*RAS_INT_START" {
            continue;
        }
        let mut pairs = Vec::new();
        while i < ls.len() && ls[i].trim_end() != "*RAS_INT_END" {
            let v: Vec<f64> = ls[i].split_whitespace().take(2).map(|t| t.parse().map_err(|_| DataError(format!("RAS: line {} is not x, y", i + 1)))).collect::<Result<_, _>>()?;
            if v.len() == 2 {
                pairs.push((v[0], v[1]));
            }
            i += 1;
        }
        let h = |k: &str| header.iter().find(|(a, _)| a == k).map(|(_, v)| v.clone()).filter(|v| !v.is_empty());
        let hf = |k: &str| h(k).and_then(|v| v.parse::<f64>().ok());
        let mut s = Scan::new(format!("sample {}", scans.len() + 1), Axis::TwoTheta, pattern(&pairs));
        s.info.sample = h("FILE_SAMPLE");
        s.info.date = h("MEAS_SCAN_START_TIME");
        s.info.step = hf("MEAS_SCAN_STEP");
        s.info.instrument = h("HW_UNIT_NAME").or_else(|| h("FILE_MD_NAME"));
        s.wavelength = hf("HW_XG_WAVE_LENGTH_ALPHA1").map(|l1| Wavelength { lambda1: l1, lambda2: hf("HW_XG_WAVE_LENGTH_ALPHA2"), ratio: None });
        scans.push(s);
    }
    if scans.is_empty() {
        return err("RAS: no *RAS_INT_START block");
    }
    Ok(Dataset { format: "rigaku_ras", scans })
}

pub fn recognise_ras(p: &Probe) -> f64 {
    if p.is_text() && p.text().starts_with("*RAS_DATA_START") && parse_ras(p).is_ok() {
        SIGNATURE
    } else {
        0.0
    }
}

pub fn read_ras(p: &Probe) -> Result<Dataset, DataError> {
    parse_ras(p)
}

fn parse_rasx(p: &Probe) -> Result<Dataset, DataError> {
    let a = Archive::open(p.bytes).map_err(|e| DataError(format!("RASX: {e}")))?;
    let prof = a.read("Data0/Profile0.txt").map_err(|e| DataError(format!("RASX: {e}")))?;
    // GSAS-II seeks past the first three bytes (a UTF-8 byte-order mark).
    let text = String::from_utf8_lossy(prof.get(3..).unwrap_or(&[])).into_owned();
    let mut pairs = Vec::new();
    for l in text.lines() {
        let v: Vec<&str> = l.split_whitespace().collect();
        if v.len() < 2 {
            continue;
        }
        pairs.push((v[0].parse().map_err(|_| DataError("RASX: x is not a number".into()))?, v[1].parse().map_err(|_| DataError("RASX: y is not a number".into()))?));
    }
    let mut s = Scan::new("bank 1", Axis::TwoTheta, pattern(&pairs));
    if let Ok(c) = a.read("Data0/MesurementConditions0.xml") {
        if let Ok(root) = xml::parse(&String::from_utf8_lossy(&c)) {
            let num = |n: &str| root.text_of(n).and_then(|t| t.parse::<f64>().ok());
            s.wavelength = num("WavelengthKalpha1").map(|l1| Wavelength { lambda1: l1, lambda2: num("WavelengthKalpha2"), ratio: None });
            s.info.sample = root.text_of("SampleName").map(String::from);
            s.info.step = num("Step");
        }
    }
    if s.pattern.is_empty() {
        return err("RASX: no data in Data0/Profile0.txt");
    }
    Ok(Dataset { format: "rigaku_rasx", scans: vec![s] })
}

pub fn recognise_rasx(p: &Probe) -> f64 {
    match p.zip_names() {
        Some(n) if n.iter().any(|n| n == "Data0/Profile0.txt") && parse_rasx(p).is_ok() => SIGNATURE,
        _ => 0.0,
    }
}

pub fn read_rasx(p: &Probe) -> Result<Dataset, DataError> {
    parse_rasx(p)
}

fn parse_txt(p: &Probe) -> Result<(Dataset, usize), DataError> {
    let ls = lines(p.text());
    let skip = ls.iter().position(|l| all_numbers(l).is_some()).ok_or(DataError("Rigaku TXT: no line of numbers".into()))?;
    let rows: Vec<Vec<f64>> = ls[skip..].iter().filter(|l| !l.trim().is_empty()).map(|l| all_numbers(l).ok_or(DataError("Rigaku TXT: a line after the header is not all numbers".into()))).collect::<Result<_, _>>()?;
    let k = rows[0].len();
    if k < 2 || rows.iter().any(|r| r.len() != k) || rows.len() < 2 {
        return err("Rigaku TXT: lines with different numbers of values");
    }
    let step = (rows[1][0] - rows[0][0]) / (k - 1) as f64;
    // The step GSAS-II warns about when it changes; here a changing step means another format.
    if step.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater) || rows.windows(2).any(|w| ((w[1][0] - w[0][0]) / (k - 1) as f64 - step).abs() > 1e-4 * step.abs().max(1e-12) + 1e-9) {
        return err("Rigaku TXT: the angle does not advance by a constant step");
    }
    let mut pairs = Vec::new();
    for r in &rows {
        let mut a = r[0];
        for &y in &r[1..] {
            pairs.push((a, y));
            a += step;
        }
    }
    let mut s = Scan::new("pattern", Axis::TwoTheta, pattern(&pairs));
    s.info.step = Some(step);
    Ok((Dataset { format: "rigaku_txt", scans: vec![s] }, k))
}

pub fn recognise_txt(p: &Probe) -> f64 {
    if !p.is_text() {
        return 0.0;
    }
    match parse_txt(p) {
        // GSAS-II's validator wants more than 30 lines; three or more intensities a line are this format's mark (two
        // are x, y, σ as well).
        Ok((d, k)) if d.scans[0].pattern.len() > 30 => {
            if k > 3 {
                STRUCTURE
            } else {
                LOOSE
            }
        }
        _ => 0.0,
    }
}

pub fn read_txt(p: &Probe) -> Result<Dataset, DataError> {
    parse_txt(p).map(|r| r.0)
}
