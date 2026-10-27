//! The formats of FullProf's list (`Ins` of the `.pcr` file), read as the Fortran that the manual prints for each of
//! them reads them (manual pp. 133–137, section "FILE.dat"), FullProf's quirks kept, and Ins = 14 as the FullProf
//! News of 24 September 2003 describes it. The readers of Ins = 0, 1, 3, 5, 6, 8 and 10 were in `data.rs` since
//! 0.6.0 and are checked against FullProf's `.prf` files (`fullprof_profiles`, `data_formats`).
//!
//! - `free_steps` (Ins = 0): up to seven lines are searched for Ti, step, Tf (a line starting "BANK" gives npts and
//!   Ti, step instead); then NPTS = int((Tf − Ti)/step + 1.5) intensities in free format. σ² = max(y, 1).
//! - `ill_d1a_d2b_old` (Ins = 1): Ti, step, Tf; NPTS pairs (counters, counts) in (10(I2,I6)); y ≤ 10⁻⁵ → 1,
//!   counters 0 → 1, σ² = y/counters.
//! - `ill_d1b_old` (Ins = 2): date, title, "io numor", two text lines, then 40 lines (F7.3, 10F6.0), the first value of
//!   the first line a1: 400 points from 2θ = 2·a1 in steps of 0.2°, y(1) = y(2), y(400) = y(399), σ² = y.
//! - `ill_d1b_d20` (Ins = 3): nset; date and title; nset1, numor, nsumad; the angles (x0 the third value, the step the
//!   ninth); NPTS; NPTS pairs (counters, counts); σ² = y/counters as for Ins = 1.
//! - `brookhaven` (Ins = 4): Ti, step, Tf; then NPTS/10 + 1 pairs of lines (10F8): intensities, then their σ; for
//!   i ≥ 2, y < 10⁻⁵ takes the previous y and σ < 10⁻⁵ is 1; σ² is the square.
//! - `two_axis` (Ins = 5): three title lines; NPTS, TSamp, TReg, Ivari, Rmon1, Rmon2; Ti, step, Tf; NPTS intensities,
//!   then their σ if Ivari ≠ 0 (y < 10⁻⁴ → 10⁻⁴, σ² < 10⁻⁶ → 1), else σ² = y·Rmon1/Rmon2 (or y).
//! - `ill_d2b` (Ins = 6): title; step in columns 17–24; Ti; Rmon, Rmon_old; then pairs (10(I2,I6)) up to the line whose
//!   first count is −1000. As in FullProf's loop, the line read before that one is not counted: NPTS is ten times the
//!   number of lines read before the line that holds −1000, less ten. σ² = cnorm·y/counters, cnorm = Rmon/Rmon_old (1
//!   if Rmon_old < 1).
//! - `ill_d4_d20l` (Ins = 7): NPTS, TSamp, title (I5,1X,F5.0,1X,A72); then NPTS triplets x, y, σ in
//!   (5(F7.2,1X,F8.1,1X,F5.1)); σ² = 0.5·σ², as the manual's code has it.
//! - `psi_dmc` (Ins = 8): two title lines; Ti, step, Tf; NPTS = int((Tf − Ti)/step + 1.005) intensities (10F8.0), then
//!   their σ (10F8.0).
//! - `socabim_uxd` (Ins = 9): the manual names the format and gives no layout. Read as the UXD keyword file of Bruker's
//!   DIFFRAC-AT: `_KEY=value` lines (`_2THETA` or `_START`, `_STEPSIZE`, `_STEPTIME`, `_WL1`, `_WL2`, `_WLRATIO`), then
//!   a data block after `_COUNTS` or `_CPS` (intensities, several per line) or `_2THETACOUNTS` or `_2THETACPS`
//!   (pairs); one scan per block; σ² = counts, or cps/steptime for count rates. Believed, not checked: no UXD file is
//!   in manuals/.
//! - `xys` (Ins = 10): the first six lines are comments (or, if the first starts with XYDATA, a header whose INTER line
//!   gives fac_x, fac_y and whose TEMP line the temperature); then x, y and optionally σ (σ = √y when absent), lines
//!   starting with '!' skipped; variable step allowed. The manual does not print this reader.
//! - `variable_time` (Ins = 11): four comment lines; Ti, step, Tf; NPTS pairs (time, counts) in (5(F6.0,F10.0));
//!   y = counts·⟨time⟩/time, σ² = y.
//! - `isis_tof_banks` (Ins = 14, News of 24 September 2003): a title; comment lines starting with '#' (one may hold
//!   "fac_y" and a factor, 1000 by default); each line starting with BANK opens a pattern of lines x, y, σ; y and σ
//!   are multiplied by fac_y and, when the BANK line holds RALF, divided by the step. The News page gives that
//!   division as an equation lost from our copy; the step is taken as x(i+1) − x(i) (the last point's as the one before
//!   it). Believed, not checked: FP_tests has Ins = 14 runs but no FullProf output for them.
//!
//! Header lines are also scanned for what they say about the measurement (`header_info`).

use super::{all_numbers, lines, Axis, Dataset, Info, Probe, Scan, Wavelength, LOOSE, SIGNATURE, STRUCTURE};
use crate::data::{DataError, Pattern};

fn err<T>(m: impl Into<String>) -> Result<T, DataError> {
    Err(DataError(m.into()))
}

/// A Fortran real: "1.5", ".02", "5.", "1.0D+03".
pub(crate) fn fnum(s: &str) -> Option<f64> {
    let s = s.trim();
    if s.is_empty() {
        return None;
    }
    s.parse::<f64>().ok().or_else(|| s.replace(['D', 'd'], "E").parse().ok())
}

/// The leading `n` reals of a line read list-directed (blanks and commas separate).
fn reals(line: &str, n: usize) -> Option<Vec<f64>> {
    let v: Vec<f64> = line.split(|c: char| c.is_whitespace() || c == ',').filter(|s| !s.is_empty()).take(n).map(fnum).collect::<Option<_>>()?;
    (v.len() == n).then_some(v)
}

/// A fixed-width field (columns from 0); blanks inside are ignored as a Fortran formatted read does (BLANK='NULL'):
/// in Cocuc.dat the pair " 10   60" read as (I2,I6) is 1 and 60. Blank → None; Some(Err) if not a number.
fn field(line: &str, from: usize, width: usize) -> Option<Result<f64, ()>> {
    let b = line.as_bytes();
    if from >= b.len() {
        return None;
    }
    let s: String = String::from_utf8_lossy(&b[from..(from + width).min(b.len())]).chars().filter(|c| !c.is_whitespace()).collect();
    if s.is_empty() {
        None
    } else {
        Some(fnum(&s).ok_or(()))
    }
}

/// Fixed-width fields of a line, `count` of `width` columns from `from`: None for a blank field, Err if one is not a
/// number.
fn fields(line: &str, from: usize, width: usize, count: usize) -> Result<Vec<Option<f64>>, ()> {
    (0..count).map(|k| field(line, from + k * width, width).transpose()).collect()
}

/// An integer field (Fortran I format): blank → None; Some(Err) if it is not an integer ("5." is not).
fn ifield(line: &str, from: usize, width: usize) -> Option<Result<f64, ()>> {
    let b = line.as_bytes();
    if from >= b.len() {
        return None;
    }
    let s: String = String::from_utf8_lossy(&b[from..(from + width).min(b.len())]).chars().filter(|c| !c.is_whitespace()).collect();
    if s.is_empty() {
        None
    } else {
        Some(s.parse::<i64>().map(|v| v as f64).map_err(|_| ()))
    }
}

/// The (counters, counts) pairs of a line in (10(I2,I6)), and how many counter fields were not blank; Err if a field
/// is not an integer.
fn pairs_i2i6(line: &str) -> Result<(Vec<(f64, f64)>, usize), ()> {
    let mut written = 0;
    let v = (0..10)
        .map(|k| {
            let n = ifield(line, 8 * k, 2).transpose()?;
            written += n.is_some() as usize;
            Ok((n.unwrap_or(0.0), ifield(line, 8 * k + 2, 6).transpose()?.unwrap_or(0.0)))
        })
        .collect::<Result<Vec<_>, ()>>()?;
    Ok((v, written))
}

/// Every token of the lines as a number, in order; None at the first token that is not one (a list-directed read
/// fails there).
fn free_numbers<'a>(ls: impl Iterator<Item = &'a str>, want: usize) -> (Vec<f64>, bool) {
    let mut v = Vec::new();
    for l in ls {
        for t in l.split(|c: char| c.is_whitespace() || c == ',').filter(|s| !s.is_empty()) {
            match fnum(t) {
                Some(x) => v.push(x),
                None => return (v, false),
            }
            if v.len() > want + 1_000_000 {
                return (v, true);
            }
        }
    }
    (v, true)
}

fn scan_from_steps(label: &str, ti: f64, step: f64, y: Vec<f64>, var: Vec<f64>) -> Scan {
    let x = (0..y.len()).map(|i| ti + i as f64 * step).collect();
    let mut s = Scan::new(label, Axis::TwoTheta, Pattern { x, y, sigma: var.iter().map(|v| v.max(0.0).sqrt()).collect() });
    s.info.step = Some(step);
    s
}

/// What header lines say: "lambda=", "WL:" (wavelength), "T=", "TEMP:" (temperature), "Date=", "M:" or
/// "Monitor/counts" (monitor), "Counting time" (time), and the title (the first line).
pub(crate) fn header_info(header: &[&str], info: &mut Info, wl: &mut Option<Wavelength>) {
    let after = |l: &str, key: &str| -> Option<f64> {
        let i = l.find(key)? + key.len();
        l[i..].split(|c: char| c.is_whitespace() || c == ',' || c == '+').find(|s| !s.is_empty()).and_then(fnum)
    };
    for (k, l) in header.iter().enumerate() {
        if k == 0 && info.title.is_none() && !l.trim().is_empty() {
            info.title = Some(l.trim().trim_start_matches(['!', '#']).trim().to_string());
        }
        let low = l.to_lowercase();
        if wl.is_none() {
            if let Some(v) = after(&low, "lambda=").or_else(|| after(&low, "wl:")).or_else(|| after(&low, "wavelength=")) {
                if v > 0.0 && v < 20.0 {
                    *wl = Some(Wavelength::single(v));
                }
            }
        }
        if info.temperature.is_none() {
            info.temperature = after(&low, "temp:").or_else(|| after(&low, " t=")).or_else(|| after(&low, "temp="));
        }
        if info.date.is_none() {
            if let Some(i) = low.find("date=") {
                info.date = Some(l[i + 5..].trim().to_string());
            }
        }
        if info.monitor.is_none() {
            if let Some(i) = low.find("monitor/counts, counting time/sec:") {
                let v: Vec<f64> = l[i + 34..].split_whitespace().filter_map(fnum).collect();
                info.monitor = v.first().copied();
                info.counting_time = v.get(1).copied();
            } else if let Some(v) = after(l, " M:") {
                info.monitor = Some(v);
            }
        }
    }
}

fn dataset(id: &'static str, s: Scan) -> Dataset {
    Dataset { format: id, scans: vec![s] }
}

fn confidence(r: Result<(Dataset, bool), DataError>) -> f64 {
    match r {
        Ok((d, c)) if d.scans.iter().any(|s| s.pattern.len() >= 4) => {
            if c {
                STRUCTURE
            } else {
                LOOSE
            }
        }
        _ => 0.0,
    }
}

/// Lines after `k` that are not blank, ILL's end marks (a line holding only −1000 or −10000) not counted.
fn rest_nonblank(ls: &[&str], k: usize) -> usize {
    ls.iter().skip(k).filter(|l| !l.trim().is_empty() && !is_end_mark(l)).count()
}

fn is_end_mark(l: &str) -> bool {
    matches!(l.trim(), "-1000" | "-10000")
}

// Ins = 0 ------------------------------------------------------------------------------------------------------------

fn parse_free(p: &Probe) -> Result<(Dataset, bool), DataError> {
    let ls = lines(p.text());
    for (k, line) in ls.iter().enumerate().take(7) {
        let (ti, step, npts, tof) = if line.starts_with("BANK") {
            // READ(aline(5:41),*) inum,npts ; READ(aline(47:90),*) thmin,step
            let a: Vec<f64> = line.get(4..41.min(line.len())).unwrap_or("").split_whitespace().filter_map(fnum).collect();
            let b: Vec<f64> = line.get(46..90.min(line.len())).unwrap_or("").split_whitespace().filter_map(fnum).collect();
            if a.len() < 2 || b.len() < 2 {
                continue;
            }
            (b[0], b[1], a[1] as usize, true)
        } else {
            let Some(v) = reals(line, 3) else { continue };
            if v[1] <= 0.0 || v[2] <= v[0] {
                continue;
            }
            (v[0], v[1], ((v[2] - v[0]) / v[1] + 1.5) as usize, false)
        };
        if !(4..=50_000_000).contains(&npts) {
            continue;
        }
        let (y, ok) = free_numbers(ls[k + 1..].iter().copied(), npts);
        if y.len() < npts {
            return err(format!("free format: {} intensities for {npts} points", y.len()));
        }
        // FullProf reads NPTS values and ignores what follows (y2o3.dat has 9 more, b8p.dat runs on past Tf); about
        // twice NPTS numbers are pairs (Ins = 1, 3) more likely than intensities.
        // Lines of (I2,I6) pairs with their counters written are pairs (Ins = 1) even when the free reading finds
        // numbers enough: there the counter and the count run together ("19    9819   10719").
        let consistent = ok && (y.len() as f64) < 1.9 * npts as f64 && !paired(&ls[k + 1..], npts);
        let y: Vec<f64> = y.into_iter().take(npts).collect();
        let var = y.iter().map(|&v| v.max(1.0)).collect();
        let mut s = scan_from_steps("pattern", ti, step, y, var);
        if tof {
            s.axis = Axis::TimeOfFlight;
        }
        header_info(&ls[..=k], &mut s.info, &mut s.wavelength);
        return Ok((dataset("free_steps", s), consistent));
    }
    err("free format: no line with Ti, step, Tf in the first seven")
}

pub fn recognise_free(p: &Probe) -> f64 {
    if !p.is_text() {
        return 0.0;
    }
    confidence(parse_free(p))
}

pub fn read_free(p: &Probe) -> Result<Dataset, DataError> {
    parse_free(p).map(|r| r.0)
}

/// Do the lines hold (I2,I6) pairs, their counters written in nine fields of ten or more?
fn paired(ls: &[&str], npts: usize) -> bool {
    let nlines = npts.div_ceil(10);
    let mut counters = 0;
    for l in ls.iter().take(nlines) {
        match pairs_i2i6(l) {
            Ok((_, c)) => counters += c,
            Err(()) => return false,
        }
    }
    counters as f64 >= 0.9 * npts as f64
}

// Ins = 1 ------------------------------------------------------------------------------------------------------------

/// σ² = y/counters with FullProf's guards (y ≤ 10⁻⁵ → 1, counters 0 → 1).
fn counted(pairs: &[(f64, f64)], cnorm: f64) -> (Vec<f64>, Vec<f64>) {
    let mut y = Vec::with_capacity(pairs.len());
    let mut var = Vec::with_capacity(pairs.len());
    for &(n, c) in pairs {
        let c = if c <= 0.00001 { 1.0 } else { c };
        let n = if n == 0.0 { 1.0 } else { n };
        y.push(c);
        var.push(cnorm * c / n);
    }
    (y, var)
}

fn parse_d1a(p: &Probe) -> Result<(Dataset, bool), DataError> {
    let ls = lines(p.text());
    let h = ls.first().and_then(|l| reals(l, 3)).ok_or(DataError("D1A/D2B format: first line must be Ti, step, Tf".into()))?;
    if h[1] <= 0.0 || h[2] <= h[0] {
        return err("D1A/D2B format: Ti, step, Tf");
    }
    let npts = ((h[2] - h[0]) / h[1] + 1.5) as usize;
    let nlines = npts.div_ceil(10);
    let mut pairs = Vec::new();
    let mut counters = 0;
    for l in ls.iter().skip(1).take(nlines) {
        let (p, c) = pairs_i2i6(l).map_err(|_| DataError("D1A/D2B format: a field is not an integer pair (I2,I6)".into()))?;
        pairs.extend(p);
        counters += c;
    }
    pairs.truncate(npts);
    if pairs.len() < npts {
        return err("D1A/D2B format: fewer points than Ti, step, Tf give");
    }
    // The counters written (a free-format file of integers in eight-column fields fits (I2,I6) with blank counters).
    // Files often hold pairs past Tf (ndni_ltm.dat: 12 lines), which FullProf does not read.
    let past = ls.iter().skip(1 + nlines).filter(|l| !l.trim().is_empty() && !is_end_mark(l));
    let consistent = past.clone().all(|l| pairs_i2i6(l).is_ok_and(|(_, c)| c >= 1)) && counters as f64 >= 0.9 * npts as f64 && ls.iter().skip(1).take(nlines).all(|l| l.len() <= 80);
    let (y, var) = counted(&pairs, 1.0);
    let mut s = scan_from_steps("pattern", h[0], h[1], y, var);
    header_info(&ls[..1], &mut s.info, &mut s.wavelength);
    Ok((dataset("ill_d1a_d2b_old", s), consistent))
}

pub fn recognise_d1a(p: &Probe) -> f64 {
    if !p.is_text() {
        return 0.0;
    }
    confidence(parse_d1a(p))
}

pub fn read_d1a(p: &Probe) -> Result<Dataset, DataError> {
    parse_d1a(p).map(|r| r.0)
}

// Ins = 2 ------------------------------------------------------------------------------------------------------------

fn parse_d1b_old(p: &Probe) -> Result<(Dataset, bool), DataError> {
    let ls = lines(p.text());
    if ls.len() < 45 {
        return err("D1B old format: 5 header lines and 40 data lines expected");
    }
    let io: Vec<i64> = ls[2].split_whitespace().take(2).filter_map(|t| t.parse().ok()).collect();
    if io.len() < 2 || rest_nonblank(&ls, 45) != 0 {
        return err("D1B old format: line 3 must hold io, numor, and 40 data lines end the file");
    }
    let mut y = Vec::with_capacity(400);
    let mut a1 = 0.0;
    for (k, l) in ls[5..45].iter().enumerate() {
        // (F7.3, 10F6.0) on the first line, (10F6.0) on the others: the manual's READ of lines 2–40 has the same
        // format, so the F7.3 field is skipped there too.
        let f0 = field(l, 0, 7).transpose().map_err(|_| DataError("D1B old format: bad angle field".into()))?;
        if k == 0 {
            a1 = f0.ok_or(DataError("D1B old format: a1 missing".into()))?;
        }
        let v = fields(l, 7, 6, 10).map_err(|_| DataError("D1B old format: bad intensity field".into()))?;
        y.extend(v.into_iter().map(|x| x.unwrap_or(0.0)));
    }
    y[0] = y[1];
    y[399] = y[398];
    let var = y.clone();
    let mut s = scan_from_steps("pattern", 2.0 * a1, 0.2, y, var);
    s.info.date = Some(ls[0].trim().to_string());
    s.info.title = Some(ls[1].trim().to_string());
    s.info.other.push(("numor".into(), format!("{}", io[1])));
    let consistent = ls[5..45].iter().all(|l| l.len() <= 67) && rest_nonblank(&ls, 45) == 0;
    Ok((dataset("ill_d1b_old", s), consistent))
}

pub fn recognise_d1b_old(p: &Probe) -> f64 {
    if !p.is_text() {
        return 0.0;
    }
    confidence(parse_d1b_old(p))
}

pub fn read_d1b_old(p: &Probe) -> Result<Dataset, DataError> {
    parse_d1b_old(p).map(|r| r.0)
}

// Ins = 3 ------------------------------------------------------------------------------------------------------------

fn parse_d20(p: &Probe) -> Result<(Dataset, bool), DataError> {
    let ls = lines(p.text());
    if ls.len() < 6 {
        return err("D1B/D20 format: header too short");
    }
    // nset alone on the first line; the angles' line is not a line of (I2,I6) pairs.
    if ls[0].split_whitespace().next().is_none_or(|t| t.parse::<i64>().is_err()) || pairs_i2i6(ls[3]).is_ok_and(|(_, c)| c >= 9) {
        return err("D1B/D20 format: nset on line 1 and the angles on line 4");
    }
    let angles: Vec<f64> = ls[3].split_whitespace().map_while(fnum).collect();
    if angles.len() < 9 {
        return err("D1B/D20 format: line 4 must hold xm, ti, x0, ome, chi, phi, tr1, tr2, step");
    }
    let (x0, step) = (angles[2], angles[8]);
    if step <= 0.0 {
        return err("D1B/D20 format: step");
    }
    let npts = ls[4].split_whitespace().next().and_then(|s| s.parse::<usize>().ok()).filter(|&n| n > 0).ok_or(DataError("D1B/D20 format: NPTS".into()))?;
    let (v, ok) = free_numbers(ls[5..].iter().copied(), 2 * npts);
    if v.len() < 2 * npts {
        return err("D1B/D20 format: fewer points than NPTS");
    }
    let pairs: Vec<(f64, f64)> = (0..npts).map(|i| (v[2 * i], v[2 * i + 1])).collect();
    let (y, var) = counted(&pairs, 1.0);
    let mut s = scan_from_steps("pattern", x0, step, y, var);
    let date_title = ls[1].trim();
    s.info.date = date_title.split_whitespace().take(2).collect::<Vec<_>>().join(" ").into();
    s.info.title = Some(date_title.to_string());
    if let Some(n) = ls[2].split_whitespace().nth(1) {
        s.info.other.push(("numor".into(), n.to_string()));
    }
    s.info.monitor = angles.first().copied();
    if let Some(&t) = angles.get(11) {
        s.info.temperature = Some(t);
    }
    // Then FullProf's inum (READ format 20), or ILL's end marks −1000 and −10000.
    let consistent = ok && v.len() >= 2 * npts && v.len() <= 2 * npts + 2;
    Ok((dataset("ill_d1b_d20", s), consistent))
}

pub fn recognise_d20(p: &Probe) -> f64 {
    if !p.is_text() {
        return 0.0;
    }
    confidence(parse_d20(p))
}

pub fn read_d20(p: &Probe) -> Result<Dataset, DataError> {
    parse_d20(p).map(|r| r.0)
}

// Ins = 4 ------------------------------------------------------------------------------------------------------------

fn parse_brookhaven(p: &Probe) -> Result<(Dataset, bool), DataError> {
    let ls = lines(p.text());
    let h = ls.first().and_then(|l| reals(l, 3)).ok_or(DataError("Brookhaven format: first line must be Ti, step, Tf".into()))?;
    if h[1] <= 0.0 || h[2] <= h[0] {
        return err("Brookhaven format: Ti, step, Tf");
    }
    let npts = ((h[2] - h[0]) / h[1] + 1.5) as usize;
    let nlines = npts / 10 + 1;
    if ls.len() < 1 + 2 * nlines {
        return err("Brookhaven format: fewer lines than NPTS needs");
    }
    let (mut y, mut var) = (Vec::new(), Vec::new());
    for k in 0..nlines {
        let a = fields(ls[1 + 2 * k], 0, 8, 10).map_err(|_| DataError("Brookhaven format: bad field".into()))?;
        let b = fields(ls[2 + 2 * k], 0, 8, 10).map_err(|_| DataError("Brookhaven format: bad field".into()))?;
        y.extend(a.into_iter().map(|v| v.unwrap_or(0.0)));
        var.extend(b.into_iter().map(|v| v.unwrap_or(0.0)));
    }
    y.truncate(npts);
    var.truncate(npts);
    var[0] *= var[0];
    for i in 1..npts {
        if y[i] < 0.00001 {
            y[i] = y[i - 1];
        }
        if var[i] < 0.00001 {
            var[i] = 1.0;
        }
        var[i] *= var[i];
    }
    let consistent = rest_nonblank(&ls, 1) == 2 * nlines && ls[1..1 + 2 * nlines].iter().all(|l| l.len() <= 80);
    let s = scan_from_steps("pattern", h[0], h[1], y, var);
    Ok((dataset("brookhaven", s), consistent))
}

pub fn recognise_brookhaven(p: &Probe) -> f64 {
    if !p.is_text() {
        return 0.0;
    }
    confidence(parse_brookhaven(p))
}

pub fn read_brookhaven(p: &Probe) -> Result<Dataset, DataError> {
    parse_brookhaven(p).map(|r| r.0)
}

// Ins = 5 ------------------------------------------------------------------------------------------------------------

fn parse_two_axis(p: &Probe) -> Result<(Dataset, bool), DataError> {
    let ls = lines(p.text());
    if ls.len() < 6 {
        return err("two-axis format: header too short");
    }
    let h: Vec<f64> = ls[3].split_whitespace().map_while(fnum).collect();
    let r = reals(ls[4], 3).ok_or(DataError("two-axis format: Ti, step, Tf".into()))?;
    if h.len() < 6 || r[1] <= 0.0 {
        return err("two-axis format: NPTS, TSamp, TReg, Ivari, Rmon1, Rmon2");
    }
    let (npts, ivari, rmon1, rmon2) = (h[0] as usize, h[3], h[4], h[5]);
    let announced = ((r[2] - r[0]) / r[1] + 1.0).round() as usize;
    if npts < 1 || h[0].fract() != 0.0 || h[3].fract() != 0.0 || announced.abs_diff(npts) > 1 {
        return err("two-axis format: NPTS must be an integer and agree with Ti, step, Tf");
    }
    let want = if ivari != 0.0 { 2 * npts } else { npts };
    let (v, ok) = free_numbers(ls[5..].iter().copied(), want);
    if v.len() < want {
        return err("two-axis format: fewer values than NPTS (and σ) need");
    }
    let y: Vec<f64> = v[..npts].to_vec();
    let var: Vec<f64> = if ivari != 0.0 {
        v[npts..2 * npts].iter().map(|s| if s * s < 0.000001 { 1.0 } else { s * s }).collect()
    } else {
        let c = if rmon1 > 1.0 && rmon2 > 1.0 { rmon1 / rmon2 } else { 1.0 };
        y.iter().map(|v| v * c).collect()
    };
    let y = if ivari != 0.0 { y.iter().map(|&v| if v < 0.0001 { 0.0001 } else { v }).collect() } else { y };
    let consistent = ok && v.len() == want;
    let mut s = scan_from_steps("pattern", r[0], r[1], y, var);
    s.info.title = Some(ls[0].trim().to_string());
    s.info.temperature = Some(h[1]);
    s.info.monitor = Some(rmon1);
    Ok((dataset("two_axis", s), consistent))
}

pub fn recognise_two_axis(p: &Probe) -> f64 {
    if !p.is_text() {
        return 0.0;
    }
    confidence(parse_two_axis(p))
}

pub fn read_two_axis(p: &Probe) -> Result<Dataset, DataError> {
    parse_two_axis(p).map(|r| r.0)
}

// Ins = 6 ------------------------------------------------------------------------------------------------------------

fn parse_d1a_standard(p: &Probe) -> Result<(Dataset, bool), DataError> {
    let ls = lines(p.text());
    if ls.len() < 5 {
        return err("D1A/D2B standard format: header too short");
    }
    // Lines 2–4 are a header, not (I2,I6) pairs: an Ins = 1 file ending in ILL's −1000 mark is not this format.
    if ls[1..4].iter().any(|l| pairs_i2i6(l).is_ok_and(|(_, c)| c >= 9)) {
        return err("D1A/D2B standard format: lines 2–4 hold data pairs, not the step, Ti and monitors");
    }
    let num = |l: &str, a, w| field(l, a, w).and_then(|r| r.ok());
    let step = num(ls[1], 16, 8).filter(|&s| s > 0.0).ok_or(DataError("D1A/D2B standard format: step (columns 17-24)".into()))?;
    let ti = num(ls[2], 0, 8).ok_or(DataError("D1A/D2B standard format: Ti".into()))?;
    let rmoni = num(ls[3], 0, 8).unwrap_or(0.0);
    let rmoniold = num(ls[3], 8, 8).unwrap_or(0.0);
    let cnorm = if rmoniold < 1.0 { 1.0 } else { rmoni / rmoniold };
    let mut pairs: Vec<(f64, f64)> = Vec::new();
    let mut found = false;
    for l in &ls[4..] {
        let (p, _) = pairs_i2i6(l).map_err(|_| DataError("D1A/D2B standard format: a field is not an integer pair (I2,I6)".into()))?;
        if (p[0].1 + 1000.0).abs() < 1e-3 {
            found = true;
            break;
        }
        pairs.extend(p);
    }
    if !found {
        return err("D1A/D2B standard format: no −1000 end mark");
    }
    // j = j − 10 after the loop, as FullProf's code has it.
    let npts = pairs.len().saturating_sub(10);
    pairs.truncate(npts);
    let (y, var) = counted(&pairs, cnorm);
    let mut s = scan_from_steps("pattern", ti, step, y, var);
    s.info.title = Some(ls[0].trim().to_string());
    s.info.monitor = Some(rmoni);
    header_info(&ls[..1], &mut s.info, &mut s.wavelength);
    Ok((dataset("ill_d2b", s), true))
}

pub fn recognise_d1a_standard(p: &Probe) -> f64 {
    if !p.is_text() {
        return 0.0;
    }
    // The −1000 end mark in the first count of a line is this format's alone.
    match parse_d1a_standard(p) {
        Ok((d, _)) if d.scans[0].pattern.len() >= 4 => SIGNATURE,
        _ => 0.0,
    }
}

pub fn read_d1a_standard(p: &Probe) -> Result<Dataset, DataError> {
    parse_d1a_standard(p).map(|r| r.0)
}

// Ins = 7 ------------------------------------------------------------------------------------------------------------

fn parse_d4(p: &Probe) -> Result<(Dataset, bool), DataError> {
    let ls = lines(p.text());
    let first = ls.first().ok_or(DataError("D4 format: empty".into()))?;
    let npts = ifield(first, 0, 5).and_then(|r| r.ok()).filter(|&n| n >= 1.0).ok_or(DataError("D4 format: NPTS (I5) in columns 1–5".into()))? as usize;
    let tsamp = field(first, 6, 5).and_then(|r| r.ok());
    let mut t = Vec::with_capacity(3 * npts);
    let mut used = 0;
    for l in &ls[1..] {
        if t.len() >= 3 * npts {
            break;
        }
        used += 1;
        for k in 0..5 {
            let o = 25 * k;
            for (a, w) in [(o, 7), (o + 8, 8), (o + 17, 5)] {
                match field(l, a, w) {
                    Some(Ok(v)) => t.push(v),
                    Some(Err(())) => return err("D4 format: a field is not a number"),
                    None => {}
                }
            }
        }
    }
    if t.len() < 3 * npts {
        return err("D4 format: fewer triplets than NPTS");
    }
    t.truncate(3 * npts);
    let x: Vec<f64> = t.iter().step_by(3).copied().collect();
    let y: Vec<f64> = t.iter().skip(1).step_by(3).copied().collect();
    let sigma: Vec<f64> = t.iter().skip(2).step_by(3).map(|s| (0.5 * s * s).sqrt()).collect();
    if rest_nonblank(&ls, 1 + used) != 0 || !x.windows(2).all(|w| w[1] > w[0]) {
        return err("D4 format: lines after the NPTS triplets, or x not increasing");
    }
    let consistent = true;
    let mut s = Scan::new("pattern", Axis::TwoTheta, Pattern { x, y, sigma });
    s.info.temperature = tsamp;
    s.info.title = first.get(12..).map(|t| t.trim().to_string());
    Ok((dataset("ill_d4_d20l", s), consistent))
}

pub fn recognise_d4(p: &Probe) -> f64 {
    if !p.is_text() {
        return 0.0;
    }
    confidence(parse_d4(p))
}

pub fn read_d4(p: &Probe) -> Result<Dataset, DataError> {
    parse_d4(p).map(|r| r.0)
}

// Ins = 8 ------------------------------------------------------------------------------------------------------------

fn parse_dmc(p: &Probe) -> Result<(Dataset, bool), DataError> {
    let ls = lines(p.text());
    let r = ls.get(2).and_then(|l| reals(l, 3)).ok_or(DataError("DMC format: Ti, step, Tf on line 3".into()))?;
    if r[1] <= 0.0 || r[2] <= r[0] {
        return err("DMC format: Ti, step, Tf");
    }
    let npts = ((r[2] - r[0]) / r[1] + 1.005) as usize;
    let mut v = Vec::with_capacity(2 * npts);
    let mut used = 0;
    for l in &ls[3..] {
        if v.len() >= 2 * npts {
            break;
        }
        used += 1;
        let f = fields(l, 0, 8, 10).map_err(|_| DataError("DMC format: a field is not a number".into()))?;
        v.extend(f.into_iter().flatten());
    }
    if v.len() < 2 * npts {
        return err("DMC format: fewer values than 2·NPTS");
    }
    // PSI files end with lines of keys and values (Numor=…, Mon=…): no more numbers after the 2·NPTS.
    let trailing = &ls[(3 + used).min(ls.len())..];
    let consistent = v.len() == 2 * npts && trailing.iter().all(|l| l.trim().is_empty() || all_numbers(l).is_none());
    let var = v[npts..2 * npts].iter().map(|s| s * s).collect();
    let mut s = scan_from_steps("pattern", r[0], r[1], v[..npts].to_vec(), var);
    header_info(&ls[..3], &mut s.info, &mut s.wavelength);
    for l in trailing.iter().filter(|l| l.contains('=')) {
        for kv in l.split(';') {
            if let Some((k, val)) = kv.split_once('=') {
                s.info.other.push((k.trim().to_string(), val.trim().to_string()));
            }
        }
    }
    Ok((dataset("psi_dmc", s), consistent))
}

pub fn recognise_dmc(p: &Probe) -> f64 {
    if !p.is_text() {
        return 0.0;
    }
    confidence(parse_dmc(p))
}

pub fn read_dmc(p: &Probe) -> Result<Dataset, DataError> {
    parse_dmc(p).map(|r| r.0)
}

// Ins = 9 ------------------------------------------------------------------------------------------------------------

fn parse_uxd(p: &Probe) -> Result<Dataset, DataError> {
    let ls = lines(p.text());
    let mut keys: Vec<(String, String)> = Vec::new();
    let get = |keys: &[(String, String)], k: &str| keys.iter().rev().find(|(a, _)| a == k).and_then(|(_, v)| fnum(v.trim_matches('\'')));
    let mut scans = Vec::new();
    let mut i = 0;
    while i < ls.len() {
        let l = ls[i].trim();
        i += 1;
        if l.is_empty() || l.starts_with(';') {
            continue;
        }
        let Some(rest) = l.strip_prefix('_') else { continue };
        let (k, v) = rest.split_once('=').map_or((rest.trim(), ""), |(a, b)| (a.trim(), b.trim()));
        let k = k.to_uppercase();
        let block = matches!(k.as_str(), "COUNTS" | "CPS" | "2THETACOUNTS" | "2THETACPS");
        if !block {
            keys.push((k, v.to_string()));
            continue;
        }
        let mut vals = Vec::new();
        while i < ls.len() {
            let t = ls[i].trim();
            if t.starts_with('_') || t.starts_with(';') {
                break;
            }
            i += 1;
            if t.is_empty() {
                continue;
            }
            vals.extend(all_numbers(t).ok_or(DataError("UXD: a data line holds something that is not a number".into()))?);
        }
        let steptime = get(&keys, "STEPTIME").filter(|&t| t > 0.0);
        let rate = k.ends_with("CPS");
        let (x, y): (Vec<f64>, Vec<f64>) = if k.starts_with("2THETA") {
            vals.chunks_exact(2).map(|c| (c[0], c[1])).unzip()
        } else {
            let start = get(&keys, "2THETA").or_else(|| get(&keys, "START")).ok_or(DataError("UXD: a _COUNTS block with no _2THETA or _START".into()))?;
            let step = get(&keys, "STEPSIZE").filter(|&s| s > 0.0).ok_or(DataError("UXD: a _COUNTS block with no _STEPSIZE".into()))?;
            ((0..vals.len()).map(|j| start + j as f64 * step).collect(), vals)
        };
        let sigma = y.iter().map(|&v| if rate { (v.max(0.0) / steptime.unwrap_or(1.0)).sqrt() } else { v.max(0.0).sqrt() }).collect();
        let mut s = Scan::new(format!("range {}", scans.len() + 1), Axis::TwoTheta, Pattern { x, y, sigma });
        s.info.step = get(&keys, "STEPSIZE");
        s.info.counting_time = steptime;
        if rate {
            s.info.other.push(("intensity".into(), crate::data_ops::RATE.into()));
        }
        if let Some(l1) = get(&keys, "WL1") {
            s.wavelength = Some(Wavelength { lambda1: l1, lambda2: get(&keys, "WL2").filter(|&v| v > 0.0), ratio: get(&keys, "WLRATIO").filter(|&v| v > 0.0) });
        }
        let txt = |k: &str| keys.iter().rev().find(|(a, _)| a == k).map(|(_, v)| v.trim_matches('\'').to_string());
        s.info.sample = txt("SAMPLE");
        s.info.date = txt("DATEMEASURED").or_else(|| txt("DATE"));
        s.info.instrument = txt("SITE");
        s.info.temperature = get(&keys, "TEMPERATURE");
        scans.push(s);
    }
    if scans.is_empty() {
        return err("UXD: no _COUNTS, _CPS, _2THETACOUNTS or _2THETACPS block");
    }
    Ok(Dataset { format: "socabim_uxd", scans })
}

pub fn recognise_uxd(p: &Probe) -> f64 {
    if !p.is_text() {
        return 0.0;
    }
    let t = p.text();
    let keyed = t.lines().take(200).any(|l| {
        let l = l.trim_start().to_uppercase();
        l.starts_with("_FILEVERSION") || l.starts_with("_2THETA=") || l.starts_with("_STEPSIZE")
    });
    match parse_uxd(p) {
        Ok(d) if keyed && d.scans.iter().any(|s| s.pattern.len() >= 4) => SIGNATURE,
        _ => 0.0,
    }
}

pub fn read_uxd(p: &Probe) -> Result<Dataset, DataError> {
    parse_uxd(p)
}

// Ins = 10 -----------------------------------------------------------------------------------------------------------

fn parse_xys(p: &Probe) -> Result<(Dataset, bool), DataError> {
    let ls = lines(p.text());
    let (mut fx, mut fy) = (1.0, 1.0);
    let xydata = ls.first().is_some_and(|l| l.starts_with("XYDATA"));
    let mut s = Scan::new("pattern", Axis::TwoTheta, Pattern::default());
    if xydata {
        for l in ls.iter().take(6).skip(1) {
            let t = l.trim_start();
            if let Some(rest) = t.strip_prefix("INTER") {
                let v: Vec<f64> = rest.split_whitespace().filter_map(fnum).collect();
                if v.len() >= 2 {
                    fx = v[0];
                    fy = v[1];
                }
            } else if let Some(rest) = t.strip_prefix("TEMP") {
                s.info.temperature = rest.split_whitespace().next().and_then(fnum);
            }
        }
        s.info.title = ls.first().map(|l| l[6..].trim().to_string()).filter(|t| !t.is_empty());
    }
    let header: Vec<&str> = ls.iter().take(6).copied().collect();
    header_info(&header[1.min(header.len())..], &mut s.info, &mut s.wavelength);
    // The first '!' line's first word names the instrument ("! D1A …"), except in WinPLOTR's header ("! File created
    // by WinPLOTR", "! File saved from WinPLOTR-2006"), whose first word is "File".
    let instrument_line = header.iter().find(|l| l.starts_with('!')).filter(|l| !l.to_lowercase().contains("winplotr"));
    if let Some(inst) = instrument_line.and_then(|l| l[1..].split_whitespace().next()) {
        if inst.len() <= 6 && inst.chars().all(|c| c.is_ascii_alphanumeric()) {
            s.info.instrument = Some(inst.to_uppercase());
        }
    }
    let mut columns_ok = true;
    for l in ls.iter().skip(6) {
        if l.starts_with('!') || l.trim().is_empty() {
            continue;
        }
        let v: Vec<f64> = l.split_whitespace().map_while(fnum).collect();
        if v.len() < 2 {
            columns_ok = false;
            continue;
        }
        if v.len() > 3 {
            columns_ok = false;
        }
        s.pattern.x.push(v[0] * fx);
        s.pattern.y.push(v[1] * fy);
        s.pattern.sigma.push(if v.len() > 2 { v[2] * fy } else { (v[1] * fy).max(0.0).sqrt() });
    }
    if s.pattern.is_empty() {
        return err("x, y, σ format: no data after the six header lines");
    }
    Ok((dataset("xys", s), xydata && columns_ok))
}

pub fn recognise_xys(p: &Probe) -> f64 {
    if !p.is_text() {
        return 0.0;
    }
    match parse_xys(p) {
        Ok((d, true)) if d.scans[0].pattern.len() >= 4 => SIGNATURE,
        Ok((d, false)) if d.scans[0].pattern.len() >= 4 => {
            // Without XYDATA, FullProf takes any six lines as comments: columns after them read loosely.
            let ls = lines(p.text());
            let data: Vec<&&str> = ls.iter().skip(6).filter(|l| !l.starts_with('!') && !l.trim().is_empty()).collect();
            if data.iter().all(|l| all_numbers(l).is_some_and(|v| v.len() >= 2)) {
                LOOSE
            } else {
                0.0
            }
        }
        _ => 0.0,
    }
}

pub fn read_xys(p: &Probe) -> Result<Dataset, DataError> {
    parse_xys(p).map(|r| r.0)
}

// Ins = 11 -----------------------------------------------------------------------------------------------------------

fn parse_variable_time(p: &Probe) -> Result<(Dataset, bool), DataError> {
    let ls = lines(p.text());
    let r = ls.get(4).and_then(|l| reals(l, 3)).ok_or(DataError("variable-time format: Ti, step, Tf on line 5".into()))?;
    if r[1] <= 0.0 || r[2] <= r[0] {
        return err("variable-time format: Ti, step, Tf");
    }
    let npts = ((r[2] - r[0]) / r[1] + 1.5) as usize;
    let mut v = Vec::with_capacity(2 * npts);
    let mut used = 0;
    for l in &ls[5..] {
        if v.len() >= 2 * npts {
            break;
        }
        used += 1;
        for k in 0..5 {
            for (a, w) in [(16 * k, 6), (16 * k + 6, 10)] {
                match field(l, a, w) {
                    Some(Ok(x)) => v.push(x),
                    Some(Err(())) => return err("variable-time format: a field is not a number"),
                    None => v.push(0.0),
                }
            }
        }
    }
    if v.len() < 2 * npts {
        return err("variable-time format: fewer pairs than NPTS");
    }
    let time: Vec<f64> = v.iter().step_by(2).take(npts).copied().collect();
    let counts: Vec<f64> = v.iter().skip(1).step_by(2).take(npts).copied().collect();
    if time.iter().any(|&t| t < 1.0e-6) {
        return err("variable-time format: zero time");
    }
    let mean = time.iter().sum::<f64>() / npts as f64;
    let y: Vec<f64> = counts.iter().zip(&time).map(|(c, t)| c * mean / t).collect();
    let var = y.clone();
    let consistent = rest_nonblank(&ls, 5 + used) == 0;
    let mut s = scan_from_steps("pattern", r[0], r[1], y, var);
    s.info.title = Some(ls[0].trim().to_string());
    s.info.counting_time = Some(mean);
    Ok((dataset("variable_time", s), consistent))
}

pub fn recognise_variable_time(p: &Probe) -> f64 {
    if !p.is_text() {
        return 0.0;
    }
    confidence(parse_variable_time(p))
}

pub fn read_variable_time(p: &Probe) -> Result<Dataset, DataError> {
    parse_variable_time(p).map(|r| r.0)
}

// Ins = 14 -----------------------------------------------------------------------------------------------------------

fn parse_isis(p: &Probe) -> Result<(Dataset, bool), DataError> {
    let ls = lines(p.text());
    let mut fac = 1000.0;
    let mut scans: Vec<Scan> = Vec::new();
    let mut ralf = false;
    let mut columns_ok = true;
    let mut info = Info { title: ls.first().map(|l| l.trim().to_string()), ..Default::default() };
    for l in ls.iter().skip(1) {
        if let Some(c) = l.strip_prefix('#') {
            if let Some(i) = c.find("fac_y") {
                if let Some(v) = c[i + 5..].split_whitespace().next().and_then(fnum) {
                    fac = v;
                }
            }
            info.other.push(("comment".into(), c.trim().to_string()));
            continue;
        }
        if l.starts_with("BANK") {
            ralf = l.contains("RALF");
            let axis = if ralf || l.contains("SLOG") || l.contains("TIME_MAP") { Axis::TimeOfFlight } else { Axis::TwoTheta };
            let mut s = Scan::new(l.split_whitespace().take(2).collect::<Vec<_>>().join(" "), axis, Pattern::default());
            s.info = info.clone();
            scans.push(s);
            continue;
        }
        if l.trim().is_empty() {
            continue;
        }
        let Some(s) = scans.last_mut() else { return err("ISIS multibank: data before the first BANK line") };
        let v = all_numbers(l).ok_or(DataError("ISIS multibank: a data line holds something that is not a number".into()))?;
        if v.len() != 3 {
            columns_ok = false;
            if v.len() < 2 {
                continue;
            }
        }
        s.pattern.x.push(v[0]);
        s.pattern.y.push(v[1] * fac);
        s.pattern.sigma.push(v.get(2).copied().unwrap_or(v[1].abs().sqrt()) * fac);
        if ralf {
            s.info.other.retain(|(k, _)| k != "ralf");
            s.info.other.push(("ralf".into(), "1".into()));
        }
    }
    for s in &mut scans {
        if s.info.other.iter().any(|(k, _)| k == "ralf") {
            s.info.other.retain(|(k, _)| k != "ralf");
            let x = &s.pattern.x;
            let n = x.len();
            if n >= 2 {
                for i in 0..n {
                    let dx = if i + 1 < n { x[i + 1] - x[i] } else { x[i] - x[i - 1] };
                    if dx > 0.0 {
                        s.pattern.y[i] /= dx;
                        s.pattern.sigma[i] /= dx;
                    }
                }
            }
        }
    }
    if scans.is_empty() || scans.iter().all(|s| s.pattern.is_empty()) {
        return err("ISIS multibank: no BANK with data");
    }
    Ok((Dataset { format: "isis_tof_banks", scans }, columns_ok))
}

pub fn recognise_isis(p: &Probe) -> f64 {
    if !p.is_text() {
        return 0.0;
    }
    // GSAS FXYE banks have this layout too; the GSAS reader recognises them by their signature, this one by structure.
    confidence(parse_isis(p))
}

pub fn read_isis(p: &Probe) -> Result<Dataset, DataError> {
    parse_isis(p).map(|r| r.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn probe(t: &str) -> Probe<'_> {
        Probe::new("t.dat", t.as_bytes())
    }

    #[test]
    fn free_and_d1a() {
        let t = "!COMM test\n 10.0 0.5 12.0\n 1 2 3\n 4 5\n";
        let d = read_free(&probe(t)).unwrap();
        assert_eq!(d.scans[0].pattern.x, [10.0, 10.5, 11.0, 11.5, 12.0]);
        assert_eq!(d.scans[0].pattern.y, [1.0, 2.0, 3.0, 4.0, 5.0]);
        assert_eq!(recognise_free(&probe(t)), STRUCTURE);
        let t = " 1.0 0.1 1.4 title\n 1    21 2    40 1     0 1    10 1    12\n";
        let d = read_d1a(&probe(t)).unwrap();
        assert_eq!(d.scans[0].pattern.y, [21.0, 40.0, 1.0, 10.0, 12.0]);
        assert!((d.scans[0].pattern.sigma[1] - 20f64.sqrt()).abs() < 1e-12);
        assert_eq!(recognise_d1a(&probe(t)), STRUCTURE);
        // The D1A file read as free format: twice as many numbers as points.
        assert_eq!(recognise_free(&probe(t)), LOOSE);
    }

    #[test]
    fn xys_with_and_without_sigma() {
        let t = "a\nb\nc\nd\ne\nf\n1.0 4.0\n! c\n2.0 9.0 2.0\n";
        let d = read_xys(&probe(t)).unwrap();
        assert_eq!(d.scans[0].pattern.x, [1.0, 2.0]);
        assert_eq!(d.scans[0].pattern.sigma, [2.0, 2.0]);
    }

    #[test]
    fn uxd_counts_and_pairs() {
        let t = "; comment\n_FILEVERSION=1\n_WL1=1.540600\n_WL2=1.544390\n_WLRATIO=0.5\n_2THETA=10.000\n_STEPSIZE=0.020\n_STEPTIME=2.0\n_COUNTS\n 100 121 144\n 169 196\n_2THETACPS\n20.0 50.0\n20.1 72.0\n20.2 8.0\n20.3 2.0\n";
        let d = read_uxd(&probe(t)).unwrap();
        assert_eq!(d.scans.len(), 2);
        assert_eq!(d.scans[0].pattern.y, [100.0, 121.0, 144.0, 169.0, 196.0]);
        assert!((d.scans[0].pattern.x[4] - 10.08).abs() < 1e-12);
        assert_eq!(d.scans[0].wavelength.unwrap().lambda2, Some(1.54439));
        assert_eq!(d.scans[1].pattern.sigma[1], 6.0);
        assert_eq!(recognise_uxd(&probe(t)), SIGNATURE);
    }

    #[test]
    fn isis_banks_scaled() {
        let t = "title\n# fac_y 10\nBANK 1 3 3 RALF 1 1 1 1 FXYE\n 100 1.0 0.1\n 102 2.0 0.2\n 106 3.0 0.3\n 110 4.0 0.4\n";
        let d = read_isis(&probe(t)).unwrap();
        let s = &d.scans[0];
        assert_eq!(s.axis, Axis::TimeOfFlight);
        assert_eq!(s.pattern.y, [5.0, 5.0, 7.5, 10.0]);
    }
}
