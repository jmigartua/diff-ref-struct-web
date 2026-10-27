//! Bruker's binary RAW files and zipped BRML files, read as GSAS-II's `G2pwd_BrukerRAW.py` and
//! `G2pwd_BrukerBRML.py` read them (offsets and record types from there; Bruker does not publish them).
//!
//! RAW versions, by their first bytes: "RAW " version 1 (DIFFRAC-AT; recognised, not read: GSAS-II does not read
//! it either and its layout is not known to us), "RAW2" version 2, "RAW1.01" version 3 (DIFFRAC plus), "RAW4.00"
//! version 4 (DIFFRAC.SUITE). Every range of the file is a scan. Counts below 1 are taken as 1 and σ² = y, as
//! GSAS-II does. The wavelengths (Kα1, Kα2, their ratio) are read from the header of versions 2 and 3 and from the
//! X-ray source record (type 30) of version 4.
//!
//! BRML is a zip of XML files; each `RawData<n>.xml` is a scan: the `Datum` lines of its data route hold, as
//! comma-separated values, the counting time, ..., 2θ (third value) and the counts (fifth, or the fourth when the
//! fifth does not vary), y = counts·time/effective time; σ² = y (no weight for y ≤ 0). The wavelength is read from
//! the `WaveLengthAlpha1/2` elements when present (GSAS-II does not read it).

use super::xml::{self, Element};
use super::zip::Archive;
use super::{Axis, Dataset, Info, Probe, Scan, Wavelength, SIGNATURE};
use crate::data::{DataError, Pattern};

fn err<T>(m: impl Into<String>) -> Result<T, DataError> {
    Err(DataError(m.into()))
}

struct Bin<'a> {
    b: &'a [u8],
}

impl Bin<'_> {
    fn bytes(&self, at: usize, n: usize) -> Result<&[u8], DataError> {
        self.b.get(at..at + n).ok_or(DataError(format!("Bruker RAW: the file ends before byte {}", at + n)))
    }
    fn u16(&self, at: usize) -> Result<u16, DataError> {
        Ok(u16::from_le_bytes(self.bytes(at, 2)?.try_into().unwrap()))
    }
    fn i32(&self, at: usize) -> Result<i32, DataError> {
        Ok(i32::from_le_bytes(self.bytes(at, 4)?.try_into().unwrap()))
    }
    fn u32(&self, at: usize) -> Result<u32, DataError> {
        Ok(u32::from_le_bytes(self.bytes(at, 4)?.try_into().unwrap()))
    }
    fn f32(&self, at: usize) -> Result<f64, DataError> {
        Ok(f32::from_le_bytes(self.bytes(at, 4)?.try_into().unwrap()) as f64)
    }
    fn f64(&self, at: usize) -> Result<f64, DataError> {
        Ok(f64::from_le_bytes(self.bytes(at, 8)?.try_into().unwrap()))
    }
    fn text(&self, at: usize, n: usize) -> Result<String, DataError> {
        Ok(self.bytes(at, n)?.iter().take_while(|&&c| c != 0).map(|&c| c as char).collect::<String>().trim().to_string())
    }
}

fn counts(bin: &Bin, at: usize, n: usize) -> Result<Vec<f64>, DataError> {
    (0..n).map(|i| bin.f32(at + 4 * i).map(|v| v.max(1.0))).collect()
}

fn scan(label: String, start: f64, step: f64, y: Vec<f64>) -> Scan {
    let x = (0..y.len()).map(|i| start + i as f64 * step).collect();
    let sigma = y.iter().map(|v| v.sqrt()).collect();
    let mut s = Scan::new(label, Axis::TwoTheta, Pattern { x, y, sigma });
    s.info.step = Some(step);
    s
}

fn wavelength(ka1: f64, ka2: f64, ratio: f64) -> Option<Wavelength> {
    (ka1 > 0.0).then_some(Wavelength { lambda1: ka1, lambda2: (ka2 > 0.0).then_some(ka2), ratio: (ka2 > 0.0 && ratio > 0.0).then_some(ratio) })
}

fn version(b: &[u8]) -> Option<u8> {
    match b.get(..7)? {
        [b'R', b'A', b'W', b' ', ..] => Some(1),
        [b'R', b'A', b'W', b'2', ..] => Some(2),
        b"RAW1.01" => Some(3),
        b"RAW4.00" => Some(4),
        _ => None,
    }
}

pub fn recognise_raw_v1(p: &Probe) -> f64 {
    if version(p.bytes) == Some(1) {
        SIGNATURE
    } else {
        0.0
    }
}

fn read_v2(b: &[u8]) -> Result<Dataset, DataError> {
    let bin = Bin { b };
    let nblock = bin.i32(4)?.max(0) as usize;
    let mut info = Info { date: Some(bin.text(168, 20)?), ..Default::default() };
    info.other.push(("anode".into(), bin.text(188, 2)?));
    let wl = wavelength(bin.f32(190)?, bin.f32(194)?, bin.f32(198)?);
    let mut pos = 256;
    let mut scans = Vec::new();
    for k in 0..nblock {
        let head = bin.u16(pos)? as usize;
        let n = bin.u16(pos + 2)? as usize;
        let step = bin.f32(pos + 12)?;
        let start = bin.f32(pos + 16)?;
        let y = counts(&bin, pos + head, n)?;
        let mut s = scan(format!("range {}", k + 1), start, step, y);
        s.info = Info { step: Some(step), ..info.clone() };
        s.wavelength = wl;
        scans.push(s);
        pos += head + 4 * n;
    }
    Ok(Dataset { format: "bruker_raw_v2", scans })
}

fn read_v3(b: &[u8]) -> Result<Dataset, DataError> {
    let bin = Bin { b };
    let nblock = bin.i32(12)?.max(0) as usize;
    let mut info = Info { date: Some(format!("{} {}", bin.text(16, 10)?, bin.text(26, 10)?)), sample: Some(bin.text(326, 60)?), ..Default::default() };
    info.other.push(("goniometer radius (mm)".into(), format!("{}", bin.f32(564)?)));
    info.other.push(("anode".into(), bin.text(608, 4)?));
    let wl = wavelength(bin.f64(624)?, bin.f64(632)?, bin.f64(648)?);
    let mut pos = 712;
    let mut scans = Vec::new();
    for k in 0..nblock {
        let mut head = bin.i32(pos)?.max(0) as usize;
        let n = bin.i32(pos + 4)?.max(0) as usize;
        if n == 0 {
            break;
        }
        head += if nblock > 1 { bin.i32(pos + 256)?.max(0) as usize } else { 40 };
        let start = bin.f64(pos + 16)?;
        let temp = bin.f32(pos + 212)?;
        let step = bin.f64(pos + 176)?;
        // GSAS-II falls back to 40 bytes earlier when the block runs past the end ("this is absurd").
        let y = counts(&bin, pos + head, n).or_else(|_| counts(&bin, pos + head - 40, n))?;
        let mut s = scan(format!("range {}", k + 1), start, step, y);
        s.info = Info { step: Some(step), temperature: (temp > 0.0).then_some(temp), ..info.clone() };
        s.wavelength = wl;
        scans.push(s);
        pos += head + 4 * n;
    }
    Ok(Dataset { format: "bruker_raw_v3", scans })
}

fn read_v4(b: &[u8]) -> Result<Dataset, DataError> {
    let bin = Bin { b };
    let nbanks = b.windows(6).filter(|w| w == b"2Theta").count();
    let mut info = Info { date: Some(format!("{} {}", bin.text(12, 12)?, bin.text(24, 10)?)), ..Default::default() };
    let mut wl = None;
    let mut pos = 61;
    let mut scans = Vec::new();
    while scans.len() < nbanks.max(1) && pos + 8 <= b.len() {
        // Header segments up to the start of a range (type 0 or 160); the data follow the range's own header (GSAS-II
        // skips the datum size at +104 and reads four bytes a point).
        let mut seg;
        loop {
            seg = bin.u32(pos)?;
            pos += 4;
            if seg == 0 || seg == 160 {
                break;
            }
            let len = (bin.u32(pos)? as usize).max(8);
            pos += 4;
            match seg {
                10 => {
                    let key = bin.text(pos + 4, 24)?;
                    let val = bin.text(pos + 28, len.saturating_sub(36))?;
                    info.other.push((key, val));
                }
                30 => {
                    let base = pos + 64;
                    wl = wavelength(bin.f64(base + 8)?, bin.f64(base + 16)?, bin.f64(base + 32)?);
                    info.other.push(("anode".into(), bin.text(base + 44, 4)?));
                }
                _ => {}
            }
            pos += len - 8;
        }
        let at = pos + 28;
        let scan_type = bin.text(at, 24)?;
        let start = bin.f64(at + 40)?;
        let step = bin.f64(at + 48)?;
        let n = bin.u32(at + 56)? as usize;
        let step_time = bin.f32(at + 60)?;
        let used_wave = bin.f64(at + 80)?;
        let mut hdr = bin.u32(at + 108)? as i64;
        let mut q = at + 128;
        if !matches!(scan_type.as_str(), "Locked Coupled" | "Unlocked Coupled" | "Detector Scan") {
            return err(format!("Bruker RAW 4: scan type \"{scan_type}\" is not read (GSAS-II reads Locked Coupled, Unlocked Coupled and Detector Scan)"));
        }
        while hdr > 0 {
            let len = (bin.u32(q + 4)? as usize).max(8);
            q += len;
            hdr -= len as i64;
        }
        let temp = bin.f32(q - 16)?;
        let y = counts(&bin, q, n)?;
        let mut s = scan(format!("range {}", scans.len() + 1), start, step, y);
        s.info = Info { step: Some(step), counting_time: Some(step_time / 1000.0), temperature: (temp > 7.0).then_some(temp), ..info.clone() };
        s.info.other.push(("scan type".into(), scan_type));
        s.wavelength = wl.or_else(|| (used_wave > 0.0).then_some(Wavelength::single(used_wave)));
        scans.push(s);
        pos = q + 4 * n;
    }
    if scans.is_empty() {
        return err("Bruker RAW 4: no range");
    }
    Ok(Dataset { format: "bruker_raw_v4", scans })
}

macro_rules! raw_version {
    ($rec:ident, $read:ident, $v:expr, $f:ident) => {
        pub fn $rec(p: &Probe) -> f64 {
            match version(p.bytes) {
                Some(v) if v == $v && $f(p.bytes).is_ok_and(|d| d.scans.iter().any(|s| s.pattern.len() >= 4)) => SIGNATURE,
                _ => 0.0,
            }
        }
        pub fn $read(p: &Probe) -> Result<Dataset, DataError> {
            if version(p.bytes) != Some($v) {
                return err(concat!("not a Bruker RAW version ", $v, " file"));
            }
            $f(p.bytes)
        }
    };
}

raw_version!(recognise_raw_v2, read_raw_v2, 2, read_v2);
raw_version!(recognise_raw_v3, read_raw_v3, 3, read_v3);
raw_version!(recognise_raw_v4, read_raw_v4, 4, read_v4);

// BRML ---------------------------------------------------------------------------------------------------------------

fn raw_data_names(names: &[String]) -> Vec<String> {
    let mut v: Vec<(i64, String)> = names
        .iter()
        .filter(|n| n.contains("RawData") && n.ends_with(".xml"))
        .map(|n| (n.split("RawData").nth(1).and_then(|r| r.trim_end_matches(".xml").parse().ok()).unwrap_or(i64::MAX), n.clone()))
        .collect();
    v.sort();
    v.into_iter().map(|(_, n)| n).collect()
}

fn f(e: &Element, name: &str) -> Option<f64> {
    e.text_of(name).and_then(|t| t.parse().ok())
}

fn brml_scan(root: &Element, label: String) -> Result<Option<Scan>, DataError> {
    let mut routes = Vec::new();
    root.find_all("DataRoute", &mut routes);
    // GSAS-II takes the second data route when there are several, else the only one.
    let route = match routes.len() {
        0 => return Ok(None),
        1 => routes[0],
        _ => routes[1],
    };
    let info = route.find("ScanInformation");
    let mut data = Vec::new();
    route.find_all("Datum", &mut data);
    let wl = root.find("WaveLengthAlpha1").and_then(|e| e.attribute("Value").and_then(|v| v.parse().ok()).or_else(|| e.text.trim().parse().ok())).map(|l1: f64| {
        let l2 = root.find("WaveLengthAlpha2").and_then(|e| e.attribute("Value").and_then(|v| v.parse().ok()));
        let r = root.find("WaveLengthRatio").and_then(|e| e.attribute("Value").and_then(|v| v.parse().ok()));
        Wavelength { lambda1: l1, lambda2: l2, ratio: r }
    });
    // A stationary detector: one Datum holding all counts from the fourth value on (GSAS-II's XtractXMLNoscan).
    if let Some(axis) = info.and_then(|i| i.find("ScaleAxisInfo")) {
        let (Some(start), Some(stop), Some(inc)) = (f(axis, "Start"), f(axis, "Stop"), f(axis, "Increment")) else { return Ok(None) };
        let n = (0.5 + (stop - start) / inc) as usize;
        if n <= 10 {
            return Ok(None);
        }
        let Some(d) = data.first() else { return Ok(None) };
        let v: Vec<&str> = d.text.split(',').collect();
        let temp: Option<f64> = v.get(2).and_then(|t| t.trim().parse().ok());
        let y: Vec<f64> = v.iter().skip(3).take(n).map(|t| t.trim().parse().map_err(|_| DataError("BRML: a count is not a number".into()))).collect::<Result<_, _>>()?;
        let x = (0..y.len()).map(|i| start + i as f64 * inc).collect();
        let sigma = y.iter().map(|&v| if v > 0.0 { v.sqrt() } else { f64::INFINITY }).collect();
        let mut s = Scan::new(label, Axis::TwoTheta, Pattern { x, y, sigma });
        s.info.temperature = temp.map(|t| 273.0 + t);
        s.info.step = Some(inc);
        s.wavelength = wl;
        return Ok(Some(s));
    }
    let n: usize = info.and_then(|i| i.text_of("MeasurementPoints")).and_then(|t| t.parse().ok()).unwrap_or(0);
    if n <= 10 {
        return Ok(None);
    }
    let eff: f64 = info.and_then(|i| f(i, "TimePerStepEffective")).ok_or(DataError("BRML: no TimePerStepEffective".into()))?;
    let (mut x, mut y, mut y3) = (Vec::with_capacity(n), Vec::with_capacity(n), Vec::with_capacity(n));
    let mut temps: Vec<f64> = Vec::new();
    for d in data.iter().take(n) {
        let v: Vec<f64> = d.text.split(',').map(|t| t.trim().parse().unwrap_or(f64::NAN)).collect();
        if v.len() < 5 {
            return err("BRML: a Datum with fewer than five values");
        }
        x.push(v[2]);
        y.push(v[4] * v[0] / eff);
        y3.push(v[3]);
        // GSAS-II looks for a temperature (°C) in the ninth, seventh or sixth value, the first that is > 0.
        for &c in &[8usize, 6, 5] {
            if let Some(&t) = v.get(c).filter(|t| **t > 0.0) {
                temps.push(t);
                break;
            }
        }
    }
    let (lo, hi) = y.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &v| (a.min(v), b.max(v)));
    if hi.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater) || (1.0 - lo / hi).abs() < 1e-4 {
        y = y3;
    }
    let sigma = y.iter().map(|&v| if v > 0.0 { v.sqrt() } else { f64::INFINITY }).collect();
    let mut s = Scan::new(label, Axis::TwoTheta, Pattern { x, y, sigma });
    s.info.counting_time = Some(eff);
    if !temps.is_empty() {
        s.info.temperature = Some(273.15 + temps.iter().sum::<f64>() / temps.len() as f64);
    }
    s.wavelength = wl;
    Ok(Some(s))
}

fn parse_brml(p: &Probe) -> Result<Dataset, DataError> {
    let a = Archive::open(p.bytes).map_err(|e| DataError(format!("BRML: {e}")))?;
    let names: Vec<String> = a.names().map(String::from).collect();
    let mut scans = Vec::new();
    for n in raw_data_names(&names) {
        let bytes = a.read(&n).map_err(|e| DataError(format!("BRML: {e}")))?;
        let root = xml::parse(&String::from_utf8_lossy(&bytes)).map_err(|e| DataError(format!("BRML {n}: {e}")))?;
        let label = n.rsplit('/').next().unwrap_or(&n).to_string();
        if let Some(mut s) = brml_scan(&root, label)? {
            s.info.date = root.text_of("TimeStampStarted").map(String::from);
            s.info.sample = root.find("Identifier").and_then(|i| i.text_of("SampleName").or_else(|| i.text_of("Name"))).map(String::from);
            scans.push(s);
        }
    }
    if scans.is_empty() {
        return err("BRML: no RawData*.xml with a scan of more than ten points");
    }
    Ok(Dataset { format: "bruker_brml", scans })
}

pub fn recognise_brml(p: &Probe) -> f64 {
    match p.zip_names() {
        Some(names) if !raw_data_names(names).is_empty() && parse_brml(p).is_ok() => SIGNATURE,
        _ => 0.0,
    }
}

pub fn read_brml(p: &Probe) -> Result<Dataset, DataError> {
    parse_brml(p)
}
