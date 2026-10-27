//! Columns of numbers: the general table reader (L061: delimiter, decimal mark, header lines and columns recognised,
//! or chosen by the user), and GSAS-II's column formats read by GSAS-II's rules (`G2pwd_xye.py`, `G2pwd_csv.py`).
//!
//! The general reader, `table`, finds:
//! - comment lines: those starting with '#', '!', '\'' or '%', and C comments `/* … */`;
//! - the delimiter: the one of ',', ';', tab and blanks that splits the most rows into the same number (≥ 2) of
//!   numbers;
//! - the decimal mark: ',' when the delimiter is not ',' and numbers are written "12,5";
//! - the header: the lines before the first row of numbers; the last of them, split by the delimiter, names the
//!   columns when it has as many fields as the rows;
//! - the columns: x, y and σ from their names (2θ, 2theta, tth, angle; d; Q; TOF, time; intensity, counts, I, y;
//!   sigma, esd, error, su, e), else the first three (σ only when there are three or more columns);
//! - the quantity on x from the x column's name (2θ, d, Q, TOF), else 2θ, noted as assumed.
//!
//! σ, when no column gives it, is √y (1 for y ≤ 0). Every choice can be given instead (`TableOptions`).
//!
//! GSAS-II's rules for its column formats:
//! - `xye` (Topas .xye and the like): comments '#', '!', '\'' and `/* … */` at the top, a first line "XYDATA" or a
//!   line "TITLE" skipped; rows of two or three numbers (separated by blanks, ',' or ';'). y ≤ 0 has no weight
//!   (σ = ∞); σ from the third column when it is not zero, else √y. GSAS-II drops the y of a row whose σ is negative
//!   but keeps its x, putting the arrays out of step; here the whole row is dropped.
//! - `qye`: the same, x being Q and the wavelength in a comment ("wavelength" and the value, third word).
//!   GSAS-II converts Q to 2θ and multiplies y by 100 on reading; here Q and y are kept as in the file.
//! - `fit2d_chi`: Fit2D's four header lines (title, "2-Theta Angle (Degrees)" or "Q (...)", "Intensity", the number
//!   of points), or pyFAI's '#' header naming "chi_2theta"; then x, y [, σ] as for xye. The Q header of a .qchi
//!   gives the wavelength after "Q".

use super::{lines, Axis, Dataset, Probe, Scan, Wavelength, LOOSE, SIGNATURE, STRUCTURE};
use crate::data::{DataError, Pattern};

fn err<T>(m: impl Into<String>) -> Result<T, DataError> {
    Err(DataError(m.into()))
}

/// How a table is read; None means "recognise it".
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TableOptions {
    /// ',', ';', '\t', or ' ' for blanks.
    pub delimiter: Option<char>,
    /// '.' or ','.
    pub decimal: Option<char>,
    /// The number of lines before the first row of data.
    pub header_lines: Option<usize>,
    /// Column numbers from 0.
    pub x: Option<usize>,
    pub y: Option<usize>,
    /// Some(None): no σ column.
    pub sigma: Option<Option<usize>>,
    pub axis: Option<Axis>,
}

/// What the reader found (or was given): the options in effect.
#[derive(Clone, Debug, PartialEq)]
pub struct TableLayout {
    pub delimiter: char,
    pub decimal: char,
    pub header_lines: usize,
    pub columns: usize,
    pub names: Vec<String>,
    pub x: usize,
    pub y: usize,
    pub sigma: Option<usize>,
    pub axis: Axis,
    pub axis_assumed: bool,
    /// Several y columns named (PL8): each y column with the σ column its name pairs it with; empty for one scan.
    pub series: Vec<(usize, Option<usize>)>,
    /// A long table's member column (PL7's long CSV: "member, x, y, sigma"): one scan per member.
    pub member: Option<usize>,
    /// The columns not read, so that the reader says so (PL8).
    pub unread: Vec<usize>,
}

fn is_comment(l: &str) -> bool {
    let t = l.trim_start();
    t.starts_with('#') || t.starts_with('!') || t.starts_with('\'') || t.starts_with('%')
}

fn split(l: &str, d: char) -> Vec<&str> {
    if d == ' ' {
        l.split_whitespace().collect()
    } else {
        l.split(d).map(str::trim).collect()
    }
}

fn number(s: &str, decimal: char) -> Option<f64> {
    let s = s.trim().trim_matches('"');
    if s.is_empty() {
        return None;
    }
    if decimal == ',' {
        s.replace(',', ".").parse().ok()
    } else {
        s.parse().ok()
    }
}

/// A row whose cells are numbers or empty (a member without a point there, PL7's wide CSV on a union grid): the first
/// cell and one more a number. Blanks never give an empty cell.
fn cells(l: &str, d: char, decimal: char) -> Option<Vec<Option<f64>>> {
    let f = split(l, d);
    if f.len() < 2 {
        return None;
    }
    let v: Vec<Option<f64>> = f.iter().map(|s| if s.trim().is_empty() && d != ' ' { Ok(None) } else { number(s, decimal).map(Some).ok_or(()) }).collect::<Result<_, _>>().ok()?;
    (v[0].is_some() && v.iter().flatten().count() >= 2).then_some(v)
}

const Y_WORDS: [&str; 9] = ["intensity", "counts", "count", "cps", "yobs", "iobs", "int", "y", "i"];
const SIGMA_WORDS: [&str; 12] = ["uncertainty", "uncert", "sigma", "error", "err", "esd", "std", "dy", "su", "σ", "e", "s"];

/// What pairs a y column with its σ column (PL8): the name without its leading y or σ word and separators ("y3" and
/// "sigma3" → "3"; "I_300K" and "esd_300K" → "300k").
fn pairing_key(name: &str, words: &[&str]) -> String {
    let n = name.trim().to_lowercase();
    let rest = words.iter().find(|w| n.starts_with(**w)).map_or(n.as_str(), |w| &n[w.len()..]);
    rest.trim_matches(|c: char| c == '_' || c == '-' || c == ' ' || c == '.' || c == '(' || c == ')').to_string()
}

/// The lines that are not comments (C comments at the top included), with their line numbers.
fn content(ls: &[&str]) -> Vec<(usize, String)> {
    let mut out = Vec::new();
    let mut in_c = false;
    for (k, l) in ls.iter().enumerate() {
        let t = l.trim();
        if in_c {
            if t.contains("*/") {
                in_c = false;
            }
            continue;
        }
        if t.starts_with("/*") {
            in_c = !t.contains("*/");
            continue;
        }
        if t.is_empty() || is_comment(l) {
            continue;
        }
        out.push((k, l.to_string()));
    }
    out
}

fn axis_of(name: &str) -> Option<Axis> {
    let n = name.to_lowercase().replace([' ', '_', '-', '(', ')', '°'], "");
    if n.starts_with("2th") || n.starts_with("twoth") || n.starts_with("2θ") || n == "tth" || n.starts_with("angle") {
        Some(Axis::TwoTheta)
    } else if n.starts_with("tof") || n.starts_with("time") {
        Some(Axis::TimeOfFlight)
    } else if n == "q" || n.starts_with("q[") || n.starts_with("qå") || n.starts_with("q1/") || n.starts_with("qa") || n.starts_with("qnm") || n == "qinverseangstrom" {
        Some(Axis::Q)
    } else if n == "d" || n.starts_with("dspac") || n.starts_with("då") || n.starts_with("d[") || n.starts_with("da") && n.len() <= 3 || n == "dangstrom" {
        Some(Axis::D)
    } else {
        None
    }
}

fn is_y(name: &str) -> bool {
    let n = name.to_lowercase();
    ["int", "count", "cps", "y", "i", "yobs", "iobs"].iter().any(|k| n == *k || n.starts_with(k) && *k != "i" && *k != "y")
}

fn is_sigma(name: &str) -> bool {
    let n = name.to_lowercase();
    ["sigma", "σ", "esd", "err", "su", "e", "s", "std", "dy", "uncert"].iter().any(|k| n == *k || n.starts_with(k) && k.len() > 2)
}

/// The layout of a table, recognised where the options leave it open.
pub fn layout(text: &str, o: &TableOptions) -> Result<TableLayout, DataError> {
    let ls = lines(text);
    let body = content(&ls);
    if body.is_empty() {
        return err("table: no lines of data");
    }
    let candidates: Vec<char> = match o.delimiter {
        Some(d) => vec![d],
        None => vec![',', ';', '\t', ' '],
    };
    // For each delimiter and decimal mark, the longest run of rows with the same number of numeric fields. A delimiter
    // absent from the data lines splits none of them (`cells` takes two fields or more), and a decimal comma without a
    // comma reads as the point does (the point's pass, first, keeps a tie): those passes are not run (W9, a large file).
    let has = |c: char| body.iter().any(|(_, l)| l.contains(c));
    let mut best: Option<(usize, char, char, usize, usize)> = None; // (rows, delimiter, decimal, first line, columns)
    for &d in &candidates {
        if d != ' ' && o.delimiter.is_none() && !has(d) {
            continue;
        }
        let decimals: Vec<char> = match o.decimal {
            Some(c) => vec![c],
            None if d == ',' => vec!['.'],
            None if !has(',') => vec!['.'],
            None => vec!['.', ','],
        };
        for &dec in &decimals {
            let mut run: Option<(usize, usize, usize)> = None; // (first line, columns, rows)
            let mut longest: Option<(usize, usize, usize)> = None;
            for (k, l) in &body {
                match cells(l, d, dec) {
                    Some(v) => match &mut run {
                        Some((_, c, n)) if *c == v.len() => *n += 1,
                        _ => run = Some((*k, v.len(), 1)),
                    },
                    None => run = None,
                }
                if let Some(r) = run {
                    if longest.is_none_or(|b| r.2 > b.2) {
                        longest = Some(r);
                    }
                }
            }
            if let Some((first, cols, rows)) = longest {
                if best.is_none_or(|b| rows > b.0) {
                    best = Some((rows, d, dec, first, cols));
                }
            }
        }
    }
    let (_, delimiter, decimal, first, columns) = best.ok_or(DataError("table: no rows of two or more numbers".into()))?;
    let header_lines = o.header_lines.unwrap_or(first);
    // Column names: the last non-blank line before the data, comment marks removed.
    let names: Vec<String> = ls[..header_lines.min(ls.len())]
        .iter()
        .rev()
        .find(|l| !l.trim().is_empty())
        .map(|l| split(l.trim().trim_start_matches(['#', '!', '\'', '%']).trim(), delimiter).iter().map(|s| s.trim_matches('"').to_string()).collect())
        .filter(|v: &Vec<String>| v.len() == columns && v.iter().any(|s| number(s, decimal).is_none()))
        .unwrap_or_default();
    let named = |f: &dyn Fn(&str) -> bool| names.iter().position(|n| f(n));
    let x = o.x.or_else(|| named(&|n| axis_of(n).is_some())).unwrap_or(0);
    let y = o.y.or_else(|| names.iter().enumerate().position(|(i, n)| i != x && is_y(n))).unwrap_or(if x == 0 { 1 } else { 0 });
    let sigma = match o.sigma {
        Some(s) => s,
        None => names.iter().enumerate().position(|(i, n)| i != x && i != y && is_sigma(n)).or_else(|| (columns >= 3 && names.is_empty()).then_some(2)),
    };
    if x >= columns || y >= columns || sigma.is_some_and(|s| s >= columns) {
        return err(format!("table: {columns} columns, the columns chosen are out of range"));
    }
    let named_axis = names.get(x).and_then(|n| axis_of(n));
    let axis = o.axis.or(named_axis).unwrap_or(Axis::TwoTheta);
    // PL8: a long table (a "member" column), or several y columns named, each σ column paired by name.
    let free = o.y.is_none() && o.sigma.is_none();
    let member = if free { names.iter().position(|n| n.trim().eq_ignore_ascii_case("member")).filter(|&m| m != x) } else { None };
    let (y, sigma) = match member {
        Some(m) => {
            let y = names.iter().enumerate().position(|(i, n)| i != x && i != m && is_y(n)).unwrap_or(y);
            (y, names.iter().enumerate().position(|(i, n)| i != x && i != m && i != y && is_sigma(n)))
        }
        None => (y, sigma),
    };
    let mut series = Vec::new();
    if free && member.is_none() && !names.is_empty() {
        let ys: Vec<usize> = (0..columns).filter(|&i| i != x && !is_sigma(&names[i])).collect();
        if ys.len() >= 2 {
            let sigmas: Vec<usize> = (0..columns).filter(|&i| i != x && is_sigma(&names[i])).collect();
            for &yi in &ys {
                let key = pairing_key(&names[yi], &Y_WORDS);
                series.push((yi, sigmas.iter().copied().find(|&si| pairing_key(&names[si], &SIGMA_WORDS) == key)));
            }
        }
    }
    let used: Vec<usize> = if series.is_empty() {
        [Some(x), Some(y), sigma, member].into_iter().flatten().collect()
    } else {
        std::iter::once(x).chain(series.iter().flat_map(|(a, b)| [Some(*a), *b]).flatten()).collect()
    };
    let unread = (0..columns).filter(|i| !used.contains(i)).collect();
    Ok(TableLayout { delimiter, decimal, header_lines, columns, names, x, y, sigma, axis, axis_assumed: o.axis.is_none() && named_axis.is_none(), series, member, unread })
}

/// What PL7's `#` lines say of a member: "file (scan 1 of 1, format), T = 4 K (source), …" (or the prototype's
/// "file (format), T = 4 K (source)", "T = not stated"): the file, the format, and the variable's name, value, unit and
/// source.
/// A member line's variable: name, value, unit, source.
type LineVariable = (String, f64, String, String);

fn member_line(t: &str) -> (String, Option<String>, Option<LineVariable>) {
    let t = t.trim();
    let Some(open) = t.find(" (") else { return (t.to_string(), None, None) };
    let file = t[..open].to_string();
    let rest = &t[open + 2..];
    let close = rest.find(')').unwrap_or(rest.len());
    let format = rest[..close].split(", ").map(str::trim).filter(|i| !i.starts_with("scan ")).last().map(String::from);
    let after = rest.get(close + 1..).unwrap_or("").trim_start_matches(',').trim();
    let var = (|| {
        let (name, v) = after.split_once(" = ")?;
        let mut w = v.splitn(2, ' ');
        let value: f64 = w.next()?.parse().ok()?;
        let tail = w.next()?;
        let (unit, src) = tail.split_once(" (")?;
        let mut depth = 1;
        let end = src.char_indices().find(|&(_, c)| {
            depth += match c {
                '(' => 1,
                ')' => -1,
                _ => 0,
            };
            depth == 0
        })?;
        Some((name.trim().to_string(), value, unit.trim().to_string(), src[..end.0].to_string()))
    })();
    (file, format, var)
}

fn describe_member(s: &mut Scan, line: Option<&String>) {
    let Some(line) = line else { return };
    let (file, format, var) = member_line(line);
    s.label = std::path::Path::new(&file).file_name().map_or(file.clone(), |f| f.to_string_lossy().into_owned());
    s.info.other.push(("member file".into(), file));
    if let Some(f) = format {
        s.info.other.push(("member format".into(), f));
    }
    if let Some((name, value, unit, source)) = var {
        s.info.other.push(("member variable".into(), name));
        s.info.other.push(("member value".into(), value.to_string()));
        s.info.other.push(("member unit".into(), unit));
        s.info.other.push(("member source".into(), source));
    }
}

/// Read a table with the given options (None: recognised). Several y columns named give one scan each (PL8), labelled
/// by PL7's `#` lines when the file has them (each member's file, format, and value of the variable with its source,
/// in `Info::other` as "member …"), else by the column's name; a long table one scan per member. A column not read is
/// said, in each scan's `Info::other` ("columns not read").
pub fn read_with(text: &str, o: &TableOptions) -> Result<(Dataset, TableLayout), DataError> {
    let lay = layout(text, o)?;
    let ls = lines(text);
    let rows: Vec<Vec<Option<f64>>> = content(&ls[lay.header_lines.min(ls.len())..]).iter().filter_map(|(_, l)| cells(l, lay.delimiter, lay.decimal)).filter(|v| v.len() == lay.columns).collect();
    // PL7's # lines: "# y3: …" (wide), "# member 3: …" (long).
    let mut said: std::collections::HashMap<String, String> = Default::default();
    for l in &ls[..lay.header_lines.min(ls.len())] {
        let t = l.trim_start().trim_start_matches('#').trim();
        if let Some((k, v)) = t.split_once(": ") {
            let k = k.trim();
            if (k.starts_with('y') && k[1..].parse::<usize>().is_ok()) || (k.starts_with("member ") && k[7..].parse::<usize>().is_ok()) {
                said.insert(k.to_string(), v.to_string());
            }
        }
    }
    let title = (lay.header_lines > 0).then(|| ls.iter().take(lay.header_lines).find(|l| !l.trim().is_empty()).map(|l| l.trim().trim_start_matches(['#', '!', '\'', '%']).trim().to_string())).flatten();
    let one = |label: String, x: usize, y: usize, sig: Option<usize>, keep: &dyn Fn(&Vec<Option<f64>>) -> bool| -> Scan {
        let mut p = Pattern::default();
        for v in rows.iter().filter(|v| keep(v)) {
            let (Some(xv), Some(yv)) = (v[x], v[y]) else { continue };
            p.x.push(xv);
            p.y.push(yv);
            p.sigma.push(match sig.and_then(|s| v[s]) {
                Some(s) => s,
                None if yv > 0.0 => yv.sqrt(),
                None => 1.0,
            });
        }
        let mut s = Scan::new(label, lay.axis, p);
        s.info.title = title.clone();
        s
    };
    let mut scans = Vec::new();
    if let Some(m) = lay.member {
        let mut ids: Vec<f64> = rows.iter().filter_map(|v| v[m]).collect();
        ids.dedup();
        let mut seen: Vec<f64> = Vec::new();
        for id in ids {
            if seen.contains(&id) {
                continue;
            }
            seen.push(id);
            let mut s = one(format!("member {id}"), lay.x, lay.y, lay.sigma, &|v| v[m] == Some(id));
            describe_member(&mut s, said.get(&format!("member {id}")));
            scans.push(s);
        }
    } else if !lay.series.is_empty() {
        for &(yi, si) in &lay.series {
            let name = lay.names[yi].clone();
            let mut s = one(name.clone(), lay.x, yi, si, &|_| true);
            describe_member(&mut s, said.get(name.trim()));
            if si.is_none() {
                s.info.other.push(("sigma".into(), format!("√y: no σ column pairs with \"{name}\"")));
            }
            scans.push(s);
        }
    } else {
        scans.push(one("table".into(), lay.x, lay.y, lay.sigma, &|_| true));
    }
    scans.retain(|s| !s.pattern.is_empty());
    if scans.is_empty() {
        return err("table: no rows of data");
    }
    for s in scans.iter_mut() {
        if lay.axis_assumed {
            s.info.other.push(("axis".into(), "assumed 2θ: no column name says what x is".into()));
        }
        if !lay.unread.is_empty() {
            let named = |i: &usize| lay.names.get(*i).map_or(format!("column {}", i + 1), |n| format!("column {} \"{n}\"", i + 1));
            s.info.other.push(("columns not read".into(), format!("{} of the file's {} columns: {}", lay.unread.len(), lay.columns, lay.unread.iter().map(named).collect::<Vec<_>>().join(", "))));
        }
    }
    Ok((Dataset { format: "table", scans }, lay))
}

pub fn recognise_table(p: &Probe) -> f64 {
    if !p.is_text() {
        return 0.0;
    }
    match read_with(p.text(), &TableOptions::default()) {
        Ok((d, lay)) if d.scans[0].pattern.len() >= 4 => {
            // A delimited table that names its columns is a table by structure; blank-separated columns are loose.
            if lay.delimiter != ' ' && !lay.names.is_empty() {
                STRUCTURE
            } else {
                LOOSE - 0.05
            }
        }
        _ => 0.0,
    }
}

pub fn read_table(p: &Probe) -> Result<Dataset, DataError> {
    read_with(p.text(), &TableOptions::default()).map(|r| r.0)
}

// GSAS-II's column formats -------------------------------------------------------------------------------------------

/// GSAS-II's xye rows: (x, y, σ) with y ≤ 0 unweighted (σ = ∞), σ = 0 → √y, negative σ dropped.
fn xye_rows<'a>(rows: impl Iterator<Item = &'a str>) -> Result<Pattern, DataError> {
    let mut p = Pattern::default();
    for l in rows {
        let v: Vec<&str> = l.split(|c: char| c.is_whitespace() || c == ',' || c == ';').filter(|s| !s.is_empty()).collect();
        if v.len() < 2 {
            continue;
        }
        let num = |s: &str| s.parse::<f64>().map_err(|_| DataError(format!("xye: not a number in \"{}\"", l.trim())));
        let (x, y) = (num(v[0])?, num(v[1])?);
        let s = if v.len() == 3 { Some(num(v[2])?) } else { None };
        let sigma = if y <= 0.0 {
            f64::INFINITY
        } else {
            match s {
                Some(s) if s < 0.0 => continue,
                Some(s) if s != 0.0 => s,
                _ => y.sqrt(),
            }
        };
        p.x.push(x);
        p.y.push(if y > 0.0 { y } else { 0.0 });
        p.sigma.push(sigma);
    }
    Ok(p)
}

/// GSAS-II's header rules for .xye: the lines of data after the comments at the top; None if a data line does not
/// hold two or three values (GSAS-II's validator).
fn xye_body(ls: &[&str]) -> Option<(Vec<String>, usize)> {
    let mut begin = true;
    let mut in_c = false;
    let mut comments = Vec::new();
    let mut start = ls.len();
    for (i, s) in ls.iter().enumerate() {
        if begin {
            let t = s.trim();
            if in_c {
                comments.push(s.to_string());
                if t.contains("*/") {
                    in_c = false;
                }
                continue;
            }
            if t.starts_with("/*") {
                in_c = !t.contains("*/");
                comments.push(s.to_string());
                continue;
            }
            if s.starts_with(['\'', '#', '!']) || (i == 0 && s.to_lowercase().contains("xydata")) || s.starts_with("TITLE") {
                comments.push(s.to_string());
                continue;
            }
            if t.is_empty() {
                continue;
            }
            begin = false;
            start = i;
        }
        let n = s.split(|c: char| c.is_whitespace() || c == ',' || c == ';').filter(|s| !s.is_empty()).count();
        if !(n == 2 || n == 3 || n == 0) {
            return None;
        }
    }
    Some((comments, start))
}

fn parse_xye(p: &Probe, q: bool) -> Result<Dataset, DataError> {
    let ls = lines(p.text());
    let (comments, start) = xye_body(&ls).ok_or(DataError("xye: a line of data does not hold two or three values".into()))?;
    let pattern = xye_rows(ls[start..].iter().copied())?;
    if pattern.is_empty() {
        return err("xye: no data");
    }
    let mut s = Scan::new("pattern", if q { Axis::Q } else { Axis::TwoTheta }, pattern);
    for c in &comments {
        let low = c.to_lowercase();
        if let Some((k, v)) = low.split_once('=') {
            if k.contains("temp") {
                s.info.temperature = v.trim().parse().ok();
            }
        }
        if q && s.wavelength.is_none() && low.contains("wavelength") {
            if let Some(w) = c.split_whitespace().nth(2).and_then(|w| w.parse().ok()) {
                s.wavelength = Some(Wavelength::single(w));
            }
        }
        if low.contains("x-axis unit is: time-of-flight") {
            s.axis = Axis::TimeOfFlight;
        }
        // The core's own header line (formats::write, PL7): "x = d (Angstrom), y, sigma", "x = Q (1/Angstrom), …".
        if low.contains("x = d (angstrom), y, sigma") {
            s.axis = Axis::D;
        } else if low.contains("x = q (1/angstrom), y, sigma") {
            s.axis = Axis::Q;
        }
        if low.contains("instrument:") {
            s.info.instrument = c.split_once(':').map(|(_, v)| v.trim().to_string());
        }
    }
    s.info.title = comments.first().map(|c| c.trim().trim_start_matches(['#', '!', '\'']).trim().to_string()).filter(|t| !t.is_empty());
    Ok(Dataset { format: if q { "qye" } else { "xye" }, scans: vec![s] })
}

pub fn recognise_xye(p: &Probe) -> f64 {
    if !p.is_text() {
        return 0.0;
    }
    match parse_xye(p, false) {
        Ok(d) if d.scans[0].pattern.len() >= 4 => LOOSE,
        _ => 0.0,
    }
}

pub fn read_xye(p: &Probe) -> Result<Dataset, DataError> {
    parse_xye(p, false)
}

pub fn recognise_qye(p: &Probe) -> f64 {
    if !p.is_text() {
        return 0.0;
    }
    match parse_xye(p, true) {
        Ok(d) if d.scans[0].pattern.len() >= 4 && d.scans[0].wavelength.is_some() => LOOSE,
        _ => 0.0,
    }
}

pub fn read_qye(p: &Probe) -> Result<Dataset, DataError> {
    parse_xye(p, true)
}

fn parse_chi(p: &Probe) -> Result<(Dataset, f64), DataError> {
    let ls = lines(p.text());
    let low: Vec<String> = ls.iter().take(4).map(|l| l.to_lowercase()).collect();
    let pyfai = low.iter().any(|l| l.contains("chi_2theta"));
    let two_theta = low.iter().skip(1).any(|l| l.contains("2-theta"));
    let q_line = low.iter().skip(1).position(|l| l.trim_start().starts_with("q ") || l.trim() == "q").map(|k| k + 1);
    let (start, conf) = if pyfai {
        (ls.iter().position(|l| !l.trim_start().starts_with('#')).unwrap_or(ls.len()), STRUCTURE)
    } else if two_theta || q_line.is_some() {
        // Fit2D: four header lines, the fourth the number of points.
        let n = ls.get(3).and_then(|l| l.trim().parse::<usize>().ok());
        let rows = ls.iter().skip(4).filter(|l| !l.trim().is_empty()).count();
        (4, if n == Some(rows) { SIGNATURE } else { STRUCTURE })
    } else {
        return err("chi: no Fit2D or pyFAI header");
    };
    let pattern = xye_rows(ls[start.min(ls.len())..].iter().copied())?;
    if pattern.is_empty() {
        return err("chi: no data");
    }
    let mut s = Scan::new("pattern", if q_line.is_some() { Axis::Q } else { Axis::TwoTheta }, pattern);
    s.info.title = ls.first().map(|l| l.trim().trim_start_matches('#').trim().to_string());
    if let Some(k) = q_line {
        if let Some(w) = ls[k].split_whitespace().nth(1).and_then(|w| w.parse().ok()) {
            s.wavelength = Some(Wavelength::single(w));
        }
    }
    Ok((Dataset { format: "fit2d_chi", scans: vec![s] }, conf))
}

pub fn recognise_chi(p: &Probe) -> f64 {
    if !p.is_text() {
        return 0.0;
    }
    match parse_chi(p) {
        Ok((d, c)) if d.scans[0].pattern.len() >= 4 => c,
        _ => 0.0,
    }
}

pub fn read_chi(p: &Probe) -> Result<Dataset, DataError> {
    parse_chi(p).map(|r| r.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn csv_with_names_and_decimal_commas() {
        let t = "# exported\n\"2Theta\";\"Intensity\";\"Error\"\n10,0;100;10\n10,1;121;11\n10,2;144;12\n10,3;169;13\n";
        let (d, lay) = read_with(t, &TableOptions::default()).unwrap();
        assert_eq!((lay.delimiter, lay.decimal, lay.header_lines), (';', ',', 2));
        assert_eq!((lay.x, lay.y, lay.sigma, lay.axis), (0, 1, Some(2), Axis::TwoTheta));
        assert_eq!(d.scans[0].pattern.x, [10.0, 10.1, 10.2, 10.3]);
        assert_eq!(d.scans[0].pattern.sigma, [10.0, 11.0, 12.0, 13.0]);
        assert_eq!(recognise_table(&Probe::new("a.csv", t.as_bytes())), STRUCTURE);
    }

    #[test]
    fn columns_by_name_and_chosen() {
        let t = "d,sigma,counts\n2.0,3,9\n1.9,4,16\n1.8,5,25\n1.7,6,36\n";
        let (d, lay) = read_with(t, &TableOptions::default()).unwrap();
        assert_eq!((lay.x, lay.y, lay.sigma, lay.axis), (0, 2, Some(1), Axis::D));
        assert_eq!(d.scans[0].pattern.y, [9.0, 16.0, 25.0, 36.0]);
        let o = TableOptions { y: Some(1), sigma: Some(None), axis: Some(Axis::Q), ..Default::default() };
        let (d, lay) = read_with(t, &o).unwrap();
        assert_eq!((lay.y, lay.sigma, lay.axis), (1, None, Axis::Q));
        assert_eq!(d.scans[0].pattern.sigma[0], 3f64.sqrt());
    }

    #[test]
    fn xye_rules() {
        let t = "' title\n1.0 4.0 0.5\n2.0 0.0 1.0\n3.0 9.0 0\n4.0 16.0 -1\n5.0 25.0\n";
        let d = read_xye(&Probe::new("a.xye", t.as_bytes())).unwrap();
        let p = &d.scans[0].pattern;
        assert_eq!(p.x, [1.0, 2.0, 3.0, 5.0]);
        assert_eq!(p.sigma, [0.5, f64::INFINITY, 3.0, 5.0]);
    }
}
