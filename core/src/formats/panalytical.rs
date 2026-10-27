//! PANalytical XRDML, read as GSAS-II's `G2pwd_Panalytical.py` reads it: the root `xrdMeasurements`, each `scan` of
//! each `xrdMeasurement` a scan, x spread evenly from `startPosition` to `endPosition` of its `positions`, y the
//! `intensities` (or `counts`), σ² = y (σ = 1 for y ≤ 0, as GSAS-II's w = 1). The wavelengths are those of
//! `usedWavelength` (kAlpha1, kAlpha2, ratioKAlpha2KAlpha1); the goniometer radius, sample id, times, tube, tension
//! and current, and a non-ambient temperature, are kept as GSAS-II keeps them.
//!
//! GSAS-II accepts any XML file with a root element; here the root must be `xrdMeasurements`.

use super::xml;
use super::{Axis, Dataset, Info, Probe, Scan, Wavelength, SIGNATURE};
use crate::data::{DataError, Pattern};

fn parse(p: &Probe) -> Result<Dataset, DataError> {
    let root = xml::parse(p.text()).map_err(|e| DataError(format!("XRDML: {e}")))?;
    if root.name != "xrdMeasurements" {
        return Err(DataError("XRDML: the root element is not xrdMeasurements".into()));
    }
    let num = |e: Option<&xml::Element>| e.and_then(|e| e.text.trim().parse::<f64>().ok());
    let mut scans = Vec::new();
    for m in root.children_named("xrdMeasurement") {
        let wave = m.child("usedWavelength");
        let wl = num(wave.and_then(|w| w.child("kAlpha1"))).map(|l1| Wavelength {
            lambda1: l1,
            lambda2: num(wave.and_then(|w| w.child("kAlpha2"))),
            ratio: num(wave.and_then(|w| w.child("ratioKAlpha2KAlpha1"))),
        });
        let incident = m.child("incidentBeamPath");
        let radius = num(incident.and_then(|i| i.child("radius"))).unwrap_or(300.0);
        let tube = incident.and_then(|i| i.child("xRayTube"));
        for (k, sc) in m.children_named("scan").enumerate() {
            let Some(dp) = sc.child("dataPoints") else { continue };
            let Some(pos) = dp.child("positions") else { continue };
            let (Some(a), Some(b)) = (num(pos.child("startPosition")), num(pos.child("endPosition"))) else { continue };
            let Some(sec) = dp.child("intensities").or_else(|| dp.child("counts")) else { continue };
            let y: Vec<f64> = sec.text.split_whitespace().map(|t| t.parse().map_err(|_| DataError("XRDML: an intensity is not a number".into()))).collect::<Result<_, _>>()?;
            let n = y.len();
            let x: Vec<f64> = (0..n).map(|i| if n > 1 { a + (b - a) * i as f64 / (n - 1) as f64 } else { a }).collect();
            let sigma = y.iter().map(|&v| if v > 0.0 { v.sqrt() } else { 1.0 }).collect();
            let mut s = Scan::new(format!("scan {}", k + 1), Axis::TwoTheta, Pattern { x, y, sigma });
            let header = sc.child("header");
            s.info = Info {
                sample: root.child("sample").and_then(|e| e.text_of("id")).map(String::from),
                date: header.and_then(|h| h.text_of("startTimeStamp")).map(String::from),
                step: (n > 1).then(|| (b - a) / (n - 1) as f64),
                counting_time: num(dp.child("commonCountingTime")),
                ..Default::default()
            };
            s.info.other.push(("goniometer radius (mm)".into(), format!("{radius}")));
            if sec.attribute("unit").is_some_and(|u| u.eq_ignore_ascii_case("cps")) {
                s.info.other.push(("intensity".into(), crate::data_ops::RATE.into()));
            }
            if let Some(t) = tube {
                s.info.instrument = t.attribute("name").map(String::from);
                for key in ["tension", "current"] {
                    if let Some(v) = t.text_of(key) {
                        s.info.other.push((key.into(), v.to_string()));
                    }
                }
            }
            if let Some(c) = sc.child("nonAmbientPoints").filter(|c| c.attribute("type") == Some("Temperature")) {
                s.info.temperature = c.text_of("nonAmbientValues").and_then(|v| v.split_whitespace().last()).and_then(|v| v.parse().ok());
            }
            s.wavelength = wl;
            scans.push(s);
        }
    }
    if scans.is_empty() {
        return Err(DataError("XRDML: no scan with positions and intensities".into()));
    }
    Ok(Dataset { format: "panalytical_xrdml", scans })
}

pub fn recognise(p: &Probe) -> f64 {
    if !p.is_text() || !p.text().contains("xrdMeasurements") {
        return 0.0;
    }
    match parse(p) {
        Ok(_) => SIGNATURE,
        Err(_) => 0.0,
    }
}

pub fn read(p: &Probe) -> Result<Dataset, DataError> {
    parse(p)
}
