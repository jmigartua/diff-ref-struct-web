//! Powder CIF (the pdCIF dictionary), read as GSAS-II's `G2pwd_CIF.py` reads it: in each data block, an x item and
//! an intensity item with the same number of values make a scan.
//! - x, in this order of preference: `_pd_meas_2theta_range_min/max/inc`, `_pd_proc_2theta_range_min/max/inc` (a
//!   range of 1 + int(0.5 + (max − min)/inc) points), `_pd_meas_2theta_scan`, `_pd_meas_time_of_flight`,
//!   `_pd_proc_2theta_corrected`, `_pd_proc_recip_len_q`;
//! - y: `_pd_meas_counts_total`, `_pd_meas_intensity_total`, `_pd_proc_intensity_net`, `_pd_proc_intensity_total`,
//!   `_pd_calc_intensity_net`, `_pd_calc_intensity_total`;
//! - σ: the standard uncertainty written with each intensity, at least 1 ("58(4)" is 58 ± 4), else √y.
//!
//! Where GSAS-II asks the user (several x or y items of the same length), the first in the order above is taken.
//! The wavelength (`_diffrn_radiation_wavelength`, one or two values) and the temperature (`_diffrn_ambient_
//! temperature`) are read always; GSAS-II reads them only in the branch where it has asked. Q is kept as Q (GSAS-II
//! converts it to 2θ with the wavelength).

use super::{Axis, Dataset, Probe, Scan, Wavelength, SIGNATURE};
use crate::cif::{self, Block};
use crate::data::{DataError, Pattern};

const X_RANGES: [[&str; 3]; 2] = [
    ["_pd_meas_2theta_range_min", "_pd_meas_2theta_range_max", "_pd_meas_2theta_range_inc"],
    ["_pd_proc_2theta_range_min", "_pd_proc_2theta_range_max", "_pd_proc_2theta_range_inc"],
];
const X_LOOPS: [(&str, Axis); 4] = [
    ("_pd_meas_2theta_scan", Axis::TwoTheta),
    ("_pd_meas_time_of_flight", Axis::TimeOfFlight),
    ("_pd_proc_2theta_corrected", Axis::TwoTheta),
    ("_pd_proc_recip_len_q", Axis::Q),
];
const Y_ITEMS: [&str; 6] = ["_pd_meas_counts_total", "_pd_meas_intensity_total", "_pd_proc_intensity_net", "_pd_proc_intensity_total", "_pd_calc_intensity_net", "_pd_calc_intensity_total"];

/// A CIF number and its standard uncertainty: "58(4)" → (58, Some(4)), "1.234(12)" → (1.234, Some(0.012)).
pub fn number_su(s: &str) -> Option<(f64, Option<f64>)> {
    let s = s.trim();
    match s.find('(') {
        Some(i) => {
            let v: f64 = s[..i].parse().ok()?;
            let digits = s[i + 1..].trim_end_matches(')');
            let su: f64 = digits.parse().ok()?;
            let decimals = s[..i].split_once('.').map_or(0, |(_, d)| d.trim_end_matches(|c: char| !c.is_ascii_digit()).len());
            let exp = s[..i].find(['e', 'E']).and_then(|k| s[k + 1..i].parse::<i32>().ok()).unwrap_or(0);
            Some((v, Some(su * 10f64.powi(exp - decimals as i32))))
        }
        None => s.parse().ok().map(|v| (v, None)),
    }
}

fn scans_of(b: &Block) -> Vec<Scan> {
    // x candidates by length, in order of preference.
    let mut xs: Vec<(usize, Vec<f64>, Axis)> = Vec::new();
    for r in X_RANGES {
        let v: Option<Vec<f64>> = r.iter().map(|t| b.items.get(*t).and_then(|v| number_su(v)).map(|n| n.0)).collect();
        if let Some(v) = v.filter(|v| v[2] != 0.0) {
            let n = 1 + (0.5 + (v[1] - v[0]) / v[2]) as usize;
            xs.push((n, (0..n).map(|i| v[0] + i as f64 * v[2]).collect(), Axis::TwoTheta));
        }
    }
    for (t, axis) in X_LOOPS {
        if let Some(c) = b.column(t) {
            xs.push((c.len(), c.iter().map(|v| number_su(v).map_or(f64::NAN, |n| n.0)).collect(), axis));
        }
    }
    let wl: Vec<f64> = match b.column("_diffrn_radiation_wavelength") {
        Some(c) => c.iter().filter_map(|v| number_su(v).map(|n| n.0)).collect(),
        None => b.get("_diffrn_radiation_wavelength").and_then(number_su).map(|n| vec![n.0]).unwrap_or_default(),
    };
    let temperature = b.get("_diffrn_ambient_temperature").and_then(number_su).map(|n| n.0);
    let mut out = Vec::new();
    let mut used = Vec::new();
    for (n, x, axis) in &xs {
        if used.contains(n) {
            continue;
        }
        let Some((yt, yc)) = Y_ITEMS.iter().find_map(|t| b.column(t).filter(|c| c.len() == *n).map(|c| (*t, c))) else { continue };
        used.push(*n);
        let (mut y, mut sigma) = (Vec::with_capacity(*n), Vec::with_capacity(*n));
        for v in yc {
            match number_su(v) {
                Some((v, Some(e))) => {
                    y.push(v);
                    sigma.push(e.max(1.0));
                }
                Some((v, None)) => {
                    y.push(v);
                    sigma.push(v.sqrt());
                }
                None => {
                    y.push(f64::NAN);
                    sigma.push(f64::NAN);
                }
            }
        }
        let mut s = Scan::new(format!("{} ({yt})", b.name), *axis, Pattern { x: x.clone(), y, sigma });
        s.wavelength = wl.first().map(|&l1| Wavelength { lambda1: l1, lambda2: wl.get(1).copied(), ratio: None });
        s.info.temperature = temperature;
        s.info.title = Some(b.name.clone());
        s.info.instrument = b.get("_pd_instr_location").filter(|v| *v != "?" && *v != ".").map(String::from);
        out.push(s);
    }
    out
}

fn parse(p: &Probe) -> Result<Dataset, DataError> {
    let scans: Vec<Scan> = cif::parse(p.text()).iter().flat_map(scans_of).collect();
    if scans.is_empty() {
        return Err(DataError("powder CIF: no block with an x item and an intensity item of the same length".into()));
    }
    Ok(Dataset { format: "powder_cif", scans })
}

pub fn recognise(p: &Probe) -> f64 {
    if !p.is_text() || !p.text().contains("_pd_") {
        return 0.0;
    }
    match parse(p) {
        Ok(d) if d.scans.iter().any(|s| s.pattern.len() >= 4) => SIGNATURE,
        _ => 0.0,
    }
}

pub fn read(p: &Probe) -> Result<Dataset, DataError> {
    parse(p)
}

#[cfg(test)]
mod tests {
    #[test]
    fn standard_uncertainties() {
        assert_eq!(super::number_su("58(4)"), Some((58.0, Some(4.0))));
        let (v, e) = super::number_su("1.234(12)").unwrap();
        assert_eq!(v, 1.234);
        assert!((e.unwrap() - 0.012).abs() < 1e-15);
    }
}
