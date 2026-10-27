//! Plots (R32, L154, L169): several data files as one set, the counterpart of WinPLOTR's buffers. The app's Plots module
//! draws what this module returns and computes nothing (L139).
//!
//! - `PlotSet`, `Member`, `Facts` (PL1): the members of a set, each one scan of one file read by the registry, with
//!   the facts the reader gives; opened from a list of files (`open`) or from WinPLOTR's `.buf` (`read_buf`,
//!   `open_buf`).
//! - `Variable`, `Found` (PL2): each member's value of the set's variable (the temperature, K) and where it was found:
//!   the reader's header, a `.buf` comment, the file's title (FullProf's free-format titles write "T.=.4.K"), the file
//!   name ("_300K", "_T300"), or typed by the user. None when nothing says it: never a guess (Josu, L169). The members
//!   sorted by it, those without one after the others.
//! - `Displacement`, `YScale`, `displace` (PL3): each member's curve displaced by the set's rule (none, a step per
//!   member, or in proportion to the variable) or its own override, on the y scale shown (linear, √, log₁₀), with its
//!   baseline and its background displaced alike.
//! - `same_grid`, `common_grid`, `resample`, `on_grid`, `map` (PL4): the members on one grid (their own points when
//!   they share them; otherwise interpolated where the grid is finer than a member, averaged over each bin where it is
//!   coarser, σ propagated), and the map: rows at the variable's values, the colour range, the difference map.
//! - `derive`, `derived_text` (PL9): the sum, average (plain or σ-weighted) and difference of members with σ, written as
//!   an XYDATA (.xys) file that names its members (the file dictionary's `derived_data`).
//! - `Window`, `WindowBackground`, `windows`, `window_value` (PL5): for named x windows, per member, the integral above
//!   a background (the automatic one, given points, or the line between the window's ends) with σ, the centroid with
//!   its σ, the strongest peak of `peak_search`; the background at chosen x: the series view (F5).
//! - `project` (PL17): the projection with hidden parts (F2): what each curve leaves visible of the curves behind it,
//!   given the displacements, the order front to back and the floor of each curve's curtain.

use crate::background::Interpolation;
use crate::data::Pattern;
use crate::data_ops::FoundPeak;
use crate::formats::{Axis, Scan};
use std::path::{Path, PathBuf};

/// What the reader says of a member's scan (PL1).
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Facts {
    /// The registry's format id.
    pub format: String,
    pub axis: Axis,
    pub points: usize,
    /// The smallest and the largest x.
    pub range: (f64, f64),
    /// The mean step, (x_last − x₀)/(n − 1) (`Pattern::mean_step`).
    pub step: Option<f64>,
    /// λ1 (Å) when the file states it.
    pub wavelength: Option<f64>,
    pub title: Option<String>,
    /// The temperature the reader found in the header (K), the monitor and the counting time.
    pub temperature: Option<f64>,
    pub monitor: Option<f64>,
    pub counting_time: Option<f64>,
}

/// Where a member's value of the variable comes from (PL2).
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Found {
    /// The reader's header (the D20 header's temperature, a TEMP line, a CIF's `_diffrn_ambient_temperature`).
    Header,
    /// The comment written after the file's name in a `.buf`, the text quoted.
    BufComment(String),
    /// The file's title, the text quoted (FullProf's free format: "T.=.4.K").
    Title(String),
    /// The file's name, the part read quoted ("300K", "T300").
    FileName(String),
    /// Typed by the user, with their note of its source.
    Typed(String),
    /// FullProf's temperatures file of a sequential run (PL21), by its name: the line "file TSAMP value" for this
    /// member's file.
    TemperatureFile(String),
    /// A CSV's `#` line for this member (PL7 written, PL8 read back), the source it states quoted.
    CsvLine(String),
}

impl Found {
    /// In words, for a table or a tooltip.
    pub fn words(&self) -> String {
        match self {
            Found::Header => "the file's header".to_string(),
            Found::BufComment(t) => format!("the buffer's comment \"{t}\""),
            Found::Title(t) => format!("the title's \"{t}\""),
            Found::FileName(t) => format!("the file name's \"{t}\""),
            Found::Typed(n) => format!("typed: {n}"),
            Found::TemperatureFile(f) => format!("the temperatures file {f}"),
            Found::CsvLine(src) => format!("the CSV's # line, which gives as its source: {src}"),
        }
    }
}

impl Found {
    /// The unit of a temperature from this source: °C for a name's "T0900C", "900C" or "900°C" (W8); kelvin for the
    /// header (the readers give kelvin, a header's °C converted), a buffer's comment, the title, a name's "300K" or
    /// "T300", FullProf's temperatures file; none for a typed value or a CSV's line, whose numbers are as given.
    pub fn unit(&self) -> Option<&'static str> {
        match self {
            Found::FileName(t) if t.ends_with('C') => Some("°C"),
            Found::Header | Found::BufComment(_) | Found::Title(_) | Found::FileName(_) | Found::TemperatureFile(_) => Some("K"),
            Found::Typed(_) | Found::CsvLine(_) => None,
        }
    }
}

/// A member's value of the set's variable, and its source.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Variable {
    pub value: f64,
    pub found: Found,
}

impl Variable {
    /// The value in `unit` (K or °C), converted (± 273.15) when its source gives the other; as it is when the source
    /// gives no unit or `unit` is neither.
    pub fn in_unit(&self, unit: &str) -> Variable {
        let value = match (self.found.unit(), unit) {
            (Some("°C"), "K") => self.value + 273.15,
            (Some("K"), "°C") => self.value - 273.15,
            _ => self.value,
        };
        Variable { value, found: self.found.clone() }
    }
}

impl PlotSet {
    /// The temperatures in `unit` (K or °C): each member's value and candidates whose source gives the other unit,
    /// converted; the set's unit set. A typed value is left as typed.
    pub fn temperatures_in(&mut self, unit: &str) {
        for m in self.members.iter_mut() {
            m.variable = m.variable.as_ref().map(|v| v.in_unit(unit));
            m.candidates = m.candidates.iter().map(|v| v.in_unit(unit)).collect();
        }
        self.variable_unit = unit.to_string();
    }
}

impl PlotSet {
    /// Every member rebinned to `step` in its x (W9): the counts averaged over each new bin, σ propagated
    /// (`data_ops::rebin_step`, as `drs data op rebin step=… mode=average`); the facts' points, range and step those of
    /// the rebinned scan.
    pub fn rebin(&mut self, step: f64) -> Result<(), String> {
        for m in self.members.iter_mut() {
            m.pattern = crate::data_ops::rebin_step(&m.pattern, step, crate::data_ops::RebinMode::Average).map_err(|e| format!("{}: {e}", m.label))?;
            m.facts.points = m.pattern.len();
            m.facts.range = crate::data::extent(&m.pattern.x).unwrap_or((0.0, 0.0));
            m.facts.step = m.pattern.mean_step();
        }
        Ok(())
    }
}

/// The unit of a set's temperatures: °C when every member's value whose source gives a unit gives °C (the ESRF names,
/// W8), else kelvin.
fn temperature_unit(members: &[Member]) -> &'static str {
    let units: Vec<&str> = members.iter().filter_map(|m| m.variable.as_ref().and_then(|v| v.found.unit())).collect();
    if !units.is_empty() && units.iter().all(|u| *u == "°C") {
        "°C"
    } else {
        "K"
    }
}

/// One member of a set: one scan of one file.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Member {
    /// The file as given (relative to the `.buf`, or as named).
    pub file: String,
    /// The scan of the file (from 0), and how many the file has.
    pub scan: usize,
    pub scans: usize,
    /// What the legend says: the file's name, with the scan's number when the file has several.
    pub label: String,
    /// The text written after the name in a `.buf`.
    pub comment: Option<String>,
    pub facts: Facts,
    /// The value of the set's variable, when something says it.
    pub variable: Option<Variable>,
    /// Every source that gave a value, in the order they are taken (the first is `variable`), so that a disagreement
    /// is seen.
    pub candidates: Vec<Variable>,
    pub pattern: Pattern,
    /// Other variables than the set's, each with its source (W3): a temperature beside the time a set is drawn
    /// against, so that a map labels its rows with both. Empty unless a source gives one.
    #[cfg_attr(feature = "serde", serde(default))]
    pub others: Vec<Other>,
    /// The file's modification time, seconds since the epoch, when the host gives it (the app reads it from the file
    /// system; the browser from the file it was handed; W2, W3).
    #[cfg_attr(feature = "serde", serde(default))]
    pub modified: Option<f64>,
}

/// A value of a variable other than the set's (W3): its name ("T", "t"), the value in the variable's unit, and where
/// it was found (`Found`: the file's header, typed in the set file …).
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Other {
    pub name: String,
    pub value: f64,
    pub found: Found,
    /// The value's unit ("K", "s", "h"; empty when the source does not say it, as a value typed in a set file).
    #[cfg_attr(feature = "serde", serde(default))]
    pub unit: String,
}

/// A file given by its name and its content (W1): what a browser reads, with no path. `comment` is a buffer's text
/// after the name (PL2's source); `variable` a value typed for the set's variable, with the note of where it came from.
#[derive(Clone, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct OpenedFile {
    pub name: String,
    pub text: String,
    #[cfg_attr(feature = "serde", serde(default))]
    pub comment: Option<String>,
    #[cfg_attr(feature = "serde", serde(default))]
    pub variable: Option<(f64, String)>,
    /// The file's modification time, seconds since the epoch, as the host read it (a browser's `File.lastModified`).
    #[cfg_attr(feature = "serde", serde(default))]
    pub modified: Option<f64>,
}

/// A set of members (PL1): its title, the variable's name and unit, the members in their order.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PlotSet {
    pub title: Option<String>,
    pub variable_name: String,
    pub variable_unit: String,
    pub members: Vec<Member>,
}

/// A WinPLOTR buffer file (`.buf`), as WinPLOTR's documents give it (`manuals/fullprof/1-WinPLOTR/docs/`):
/// - one data file per line, "without the directory name", optionally followed by "informations about the files
///   (description of the experiment, temperatures)": `file_1.dat      T=5K`, `file_2.dat      T = 10 K` (the user's
///   guide, March 2003, p. 19, §5.1.4, examples 6 and 7);
/// - "the title of the plot in the first line if begun by COMMENT or TITLE" (WinPLOTR's news, July 2013);
/// - empty lines allowed (the news, May 1998).
///
/// A name is relative to the buffer. The two buffers of Tb₂Pd₂Sn hold names only, which the documents allow. One rule
/// is ours: a line is a file's name whole when that file exists (a name with a blank), else its first word.
#[derive(Clone, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Buffer {
    pub title: Option<String>,
    pub entries: Vec<(String, Option<String>)>,
}

/// Read a buffer's text; `exists` says whether a name (relative to the buffer) is a file.
pub fn read_buf(text: &str, exists: &dyn Fn(&str) -> bool) -> Buffer {
    let mut b = Buffer::default();
    for line in text.lines() {
        let l = line.trim();
        if l.is_empty() {
            continue;
        }
        // The first line that is not empty (the news's "first line"; empty lines are allowed).
        if b.title.is_none() && b.entries.is_empty() {
            let up = l.to_uppercase();
            if let Some(rest) = ["TITLE", "COMMENT"].iter().find_map(|w| up.starts_with(w).then(|| l[w.len()..].trim_start_matches([':', '=', ' ']).trim())) {
                b.title = Some(rest.to_string());
                continue;
            }
        }
        if exists(l) {
            b.entries.push((l.to_string(), None));
            continue;
        }
        match l.split_once(char::is_whitespace) {
            Some((name, rest)) if !rest.trim().is_empty() => b.entries.push((name.to_string(), Some(rest.trim().to_string()))),
            _ => b.entries.push((l.to_string(), None)),
        }
    }
    b
}

/// A temperature written "T = 4 K", "T.=.4.K" (FullProf's free-format titles, dots for blanks), "Temp: 300K" or
/// "temperature = 1.5 K": the letter T (or "Temp", "Temperature") not inside a word, then '=' or ':', then a number,
/// then K. None when the text says none.
pub fn temperature_in(text: &str) -> Option<(f64, String)> {
    let c: Vec<char> = text.chars().collect();
    let sep = |ch: char| ch == '.' || ch == ' ' || ch == '_';
    for i in 0..c.len() {
        if !(c[i] == 'T' || c[i] == 't') || (i > 0 && c[i - 1].is_alphanumeric()) {
            continue;
        }
        let mut j = i + 1;
        // "Temp", "Temperature".
        for word in ["emperature", "emp"] {
            let w: Vec<char> = word.chars().collect();
            if c.len() >= j + w.len() && c[j..j + w.len()].iter().zip(&w).all(|(a, b)| a.to_ascii_lowercase() == *b) {
                j += w.len();
                break;
            }
        }
        while j < c.len() && sep(c[j]) {
            j += 1;
        }
        if j >= c.len() || !(c[j] == '=' || c[j] == ':') {
            continue;
        }
        j += 1;
        while j < c.len() && sep(c[j]) {
            j += 1;
        }
        let Some((v, end)) = number_at(&c, j) else { continue };
        let mut k = end;
        while k < c.len() && sep(c[k]) {
            k += 1;
        }
        if k < c.len() && (c[k] == 'K' || c[k] == 'k') && (k + 1 == c.len() || !c[k + 1].is_alphanumeric()) {
            return Some((v, c[i..=k].iter().collect()));
        }
    }
    None
}

/// A number at `j`: optional sign, digits, and a decimal part only when a digit follows the point (so that "4.K"
/// reads 4).
fn number_at(c: &[char], j: usize) -> Option<(f64, usize)> {
    let mut k = j;
    if k < c.len() && (c[k] == '-' || c[k] == '+') {
        k += 1;
    }
    let d0 = k;
    while k < c.len() && c[k].is_ascii_digit() {
        k += 1;
    }
    if k == d0 {
        return None;
    }
    if k + 1 < c.len() && (c[k] == '.' || c[k] == 'p') && c[k + 1].is_ascii_digit() {
        k += 1;
        while k < c.len() && c[k].is_ascii_digit() {
            k += 1;
        }
    }
    let s: String = c[j..k].iter().map(|&ch| if ch == 'p' { '.' } else { ch }).collect();
    s.parse().ok().map(|v| (v, k))
}

/// The elapsed seconds a file's name gives (W3's t): a number followed by s, "sec" or "secs" ("_3600s", "_12.5sec"),
/// or a lower-case t and such a number ("_t00000s", the ESRF run's names, W8), not inside a longer word; the part read
/// quoted.
pub fn elapsed_in_name(name: &str) -> Option<(f64, String)> {
    let stem = Path::new(name).file_stem().map_or(name.to_string(), |s| s.to_string_lossy().into_owned());
    let c: Vec<char> = stem.chars().collect();
    for j in 0..c.len() {
        // a lower-case t before the number, itself not inside a word, says a time (as T says a temperature)
        let before_t = j > 0 && c[j - 1] == 't' && (j == 1 || !c[j - 2].is_alphanumeric());
        if !c[j].is_ascii_digit() || (j > 0 && (c[j - 1].is_ascii_digit() || c[j - 1] == '.' || (c[j - 1].is_alphabetic() && !before_t))) {
            continue;
        }
        let Some((v, end)) = number_at(&c, j) else { continue };
        for unit in ["secs", "sec", "s"] {
            let u: Vec<char> = unit.chars().collect();
            let stop = end + u.len();
            if stop <= c.len() && c[end..stop].iter().map(|x| x.to_ascii_lowercase()).eq(u.iter().copied()) && (stop == c.len() || !c[stop].is_alphabetic()) {
                return Some((v, c[if before_t { j - 1 } else { j }..stop].iter().collect()));
            }
        }
    }
    None
}

impl PlotSet {
    /// W3's t for each member without one: the seconds since the set's first member, from the elapsed seconds its name
    /// gives (`elapsed_in_name`) when every member's name gives them, else, with `modified`, from the modification
    /// times when every member has one; nothing otherwise (never a mix of the two sources). A set opened takes the
    /// names' only: a file's modification time says when it was copied as often as when it was measured, so it is
    /// taken only where the host says the times are the run's (W2's `propose`).
    pub fn fill_elapsed(&mut self, modified: bool) {
        let from_names: Option<Vec<(f64, String)>> = self.members.iter().map(|m| elapsed_in_name(&m.file)).collect();
        let times: Option<Vec<(f64, Found)>> = match from_names {
            Some(v) if !v.is_empty() => Some(v.into_iter().map(|(s, t)| (s, Found::FileName(t))).collect()),
            _ if modified => self.members.iter().map(|m| m.modified.map(|t| (t, Found::Typed("the file's modification time".into())))).collect(),
            _ => None,
        };
        let Some(times) = times.filter(|t| !t.is_empty()) else { return };
        let first = times.iter().map(|t| t.0).fold(f64::INFINITY, f64::min);
        for (m, (t, found)) in self.members.iter_mut().zip(times) {
            if !m.others.iter().any(|o| o.name == "t") {
                m.others.push(Other { name: "t".into(), value: t - first, found, unit: "s".into() });
            }
        }
    }
}

/// The set's y extent on a scale (PL41, W4's shared y interval): the smallest and the largest of the members' counts
/// (those named, by their index in the set, else all) put through the scale, finite values only (a log of y ≤ 0 is
/// none); none when no value is left.
pub fn extent(set: &PlotSet, scale: YScale, members: Option<&[usize]>) -> Option<(f64, f64)> {
    let ys: Vec<f64> = set.members.iter().enumerate().filter(|(k, _)| members.is_none_or(|ms| ms.contains(k))).flat_map(|(_, m)| m.pattern.y.iter().map(|&y| scale.apply(y))).collect();
    crate::data::extent(&ys)
}

/// A temperature in a file's name: a number followed by K ("_300K", "1p5K" for 1.5 K) not inside a longer word, or an
/// upper-case T and a number ("_T300"), in kelvin; a number followed by an upper-case C or °C, with or without the T
/// ("T0900C", "900C", "900°C"), in °C (W8: the ESRF run's names; `Found::unit` says which). A lower-case "t01" is a
/// file's index, not a temperature, and is not read.
pub fn temperature_in_name(name: &str) -> Option<(f64, String)> {
    let stem = Path::new(name).file_stem().map_or(name.to_string(), |s| s.to_string_lossy().into_owned());
    let c: Vec<char> = stem.chars().collect();
    for j in 0..c.len() {
        if !c[j].is_ascii_digit() || (j > 0 && (c[j - 1].is_ascii_digit() || c[j - 1] == '.')) {
            continue;
        }
        let Some((v, end)) = number_at(&c, j) else { continue };
        let before_t = j > 0 && c[j - 1] == 'T' && (j == 1 || !c[j - 2].is_alphanumeric());
        let k_after = end < c.len() && (c[end] == 'K' || c[end] == 'k') && (end + 1 == c.len() || !c[end + 1].is_alphabetic());
        // degrees Celsius (W8): "T0900C", "900C", "900°C", an upper-case C only, not inside a longer word
        let c_end = if end < c.len() && c[end] == 'C' { Some(end) } else if end + 1 < c.len() && c[end] == '°' && c[end + 1] == 'C' { Some(end + 1) } else { None };
        let c_after = c_end.is_some_and(|e| e + 1 == c.len() || !c[e + 1].is_alphabetic());
        let free_before = j == 0 || !c[j - 1].is_alphabetic() || before_t;
        if k_after && free_before {
            let s = if before_t { j - 1 } else { j };
            return Some((v, c[s..=end].iter().collect()));
        }
        if c_after && free_before {
            let s = if before_t { j - 1 } else { j };
            return Some((v, c[s..=c_end.unwrap()].iter().collect()));
        }
        if before_t && (end == c.len() || !c[end].is_alphanumeric()) {
            return Some((v, c[j - 1..end].iter().collect()));
        }
    }
    None
}

/// Every value of the variable that the member's sources give, in the order they are taken: the reader's header, the
/// buffer's comment, the title, the file name.
pub fn temperature_candidates(file: &str, comment: Option<&str>, scan: &Scan) -> Vec<Variable> {
    let mut out = Vec::new();
    // A member of a CSV that describes it (PL8): what its # line says, and its own file's name; not the CSV's title
    // nor the CSV's name, which are the set's.
    let other = |k: &str| scan.info.other.iter().find(|o| o.0 == k).map(|o| o.1.clone());
    if let Some(own) = other("member file") {
        if let (Some(v), Some(src)) = (other("member value").and_then(|v| v.parse::<f64>().ok()), other("member source")) {
            out.push(Variable { value: v, found: Found::CsvLine(src) });
        }
        if let Some((v, t)) = temperature_in_name(&own) {
            out.push(Variable { value: v, found: Found::FileName(t) });
        }
        return out;
    }
    if let Some(t) = scan.info.temperature {
        out.push(Variable { value: t, found: Found::Header });
    }
    if let Some((v, _)) = comment.and_then(temperature_in) {
        out.push(Variable { value: v, found: Found::BufComment(comment.unwrap().to_string()) });
    }
    if let Some((v, t)) = scan.info.title.as_deref().and_then(temperature_in) {
        out.push(Variable { value: v, found: Found::Title(t) });
    }
    if let Some((v, t)) = temperature_in_name(file) {
        out.push(Variable { value: v, found: Found::FileName(t) });
    }
    out
}

fn facts(format: &str, s: &Scan) -> Facts {
    let p = &s.pattern;
    Facts {
        format: format.to_string(),
        axis: s.axis,
        points: p.len(),
        range: crate::data::extent(&p.x).unwrap_or((0.0, 0.0)),
        step: p.mean_step(),
        wavelength: s.wavelength.map(|w| w.lambda1),
        title: s.info.title.clone(),
        temperature: s.info.temperature,
        monitor: s.info.monitor,
        counting_time: s.info.counting_time,
    }
}

/// Open files as a set (PL1, PL2): each file read by the registry (`reader` when given), every scan a member, each
/// member's temperature from its sources; the members sorted by it, those without one after the others in the order
/// given. `base` is the folder the names are relative to; `comments` the buffer's texts by file.
pub fn open(base: &Path, files: &[(String, Option<String>)], reader: Option<&str>, title: Option<String>) -> Result<PlotSet, String> {
    open_with(base, files, reader, title, &|_| Vec::new())
}

/// The operations Data keeps for a file (PL20), by the file's name as the set gives it: the caller's (the app keeps
/// them in its own state; a script names them), applied by the core as `data_ops::apply_all` applies them.
pub type Operations<'a> = dyn Fn(&str) -> Vec<crate::data_ops::Operation> + 'a;

/// `open`, each member drawn after the operations Data keeps for its file (PL20): masked, rebinned, normalised …
/// as Data shows it. A file another operation names (add, subtract, merge) is read relative to `base`. The member's
/// facts are the scan's after the operations (points, range, step); its header's (format, λ, title, temperature …)
/// as read.
pub fn open_with(base: &Path, files: &[(String, Option<String>)], reader: Option<&str>, title: Option<String>, operations: &Operations) -> Result<PlotSet, String> {
    let mut given = Vec::new();
    for (file, comment) in files {
        let path: PathBuf = base.join(file);
        let bytes = std::fs::read(&path).map_err(|e| format!("{file}: {}: {e}", path.display()))?;
        // the file system's modification time, as the app's members' list shows it (W2's time, W3's t)
        let modified = std::fs::metadata(&path).ok().and_then(|m| m.modified().ok()).and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok()).map(|d| d.as_secs_f64());
        given.push((file.clone(), comment.clone(), None, bytes, modified));
    }
    let load = |f: &str| std::fs::read(base.join(f)).map_err(|e| format!("{f}: {}: {e}", base.join(f).display()));
    members_of(given, reader, title, operations, &load)
}

/// The set from files given by name and content (W1), as `open` reads them from paths: the registry identifies each
/// on its text and name (the name's extension only), nothing is read from disk. A file a member's operations name
/// (add, subtract, merge) is taken from `files` by its name. A typed `variable` is the member's value, its note the
/// source.
pub fn open_from(files: &[OpenedFile], reader: Option<&str>, title: Option<String>) -> Result<PlotSet, String> {
    open_from_with(files, reader, title, &|_| Vec::new())
}

/// `open_from` with Data's operations per file (PL20).
pub fn open_from_with(files: &[OpenedFile], reader: Option<&str>, title: Option<String>, operations: &Operations) -> Result<PlotSet, String> {
    let given = files.iter().map(|f| (f.name.clone(), f.comment.clone(), f.variable.clone(), f.text.as_bytes().to_vec(), f.modified)).collect();
    let load = |name: &str| files.iter().find(|f| f.name == name).map(|f| f.text.as_bytes().to_vec()).ok_or_else(|| format!("{name}: an operation names it, and no file of that name was given"));
    members_of(given, reader, title, operations, &load)
}

/// The members of files given as (name, comment, typed value, bytes), each after its operations, the files those name
/// read by `load`; sorted by the variable.
#[allow(clippy::type_complexity)]
fn members_of(given: Vec<(String, Option<String>, Option<(f64, String)>, Vec<u8>, Option<f64>)>, reader: Option<&str>, title: Option<String>, operations: &Operations, load: &dyn Fn(&str) -> Result<Vec<u8>, String>) -> Result<PlotSet, String> {
    let scan_of = |f: &str, r: Option<&str>| -> Result<Scan, String> {
        let bytes = load(f)?;
        crate::formats::read(&crate::formats::Probe::new(f, &bytes), r).map_err(|e| format!("{f}: {e}"))?.scans.into_iter().next().ok_or(format!("{f}: no scan"))
    };
    let mut members = Vec::new();
    for (file, comment, typed, bytes, modified) in given {
        let ds = crate::formats::read(&crate::formats::Probe::new(&file, &bytes), reader).map_err(|e| format!("{file}: {e}"))?;
        let ops = operations(&file);
        let n = ds.scans.len();
        for (k, s) in ds.scans.into_iter().enumerate() {
            let s = if ops.is_empty() { s } else { crate::data_ops::apply_all(&s, &ops, &scan_of).map_err(|e| format!("{file}: {e}"))? };
            let mut candidates = temperature_candidates(&file, comment.as_deref(), &s);
            let name = Path::new(&file).file_name().map_or(file.clone(), |f| f.to_string_lossy().into_owned());
            let mut variable = candidates.first().cloned();
            if let Some((v, note)) = &typed {
                let t = Variable { value: *v, found: Found::Typed(note.clone()) };
                candidates.insert(0, t.clone());
                variable = Some(t);
            }
            members.push(Member {
                file: file.clone(),
                scan: k,
                scans: n,
                // A table of several members (PL8) labels each by its # line's file or its column's name.
                label: if n > 1 && ds.format == "table" { s.label.clone() } else if n > 1 { format!("{name} ({})", k + 1) } else { name },
                comment: comment.clone(),
                facts: facts(ds.format, &s),
                variable,
                candidates,
                pattern: s.pattern,
                others: Vec::new(),
                modified,
            });
        }
    }
    let unit = temperature_unit(&members);
    let mut set = PlotSet { title, variable_name: "T".into(), variable_unit: "K".into(), members };
    set.temperatures_in(unit);
    set.sort_by_variable();
    // W3's t from the elapsed seconds the names give, when they all give them
    set.fill_elapsed(false);
    Ok(set)
}

/// Open a WinPLOTR buffer as a set (PL1): its names relative to it, its comments, its title.
pub fn open_buf(path: &Path, reader: Option<&str>) -> Result<PlotSet, String> {
    open_buf_with(path, reader, &|_| Vec::new())
}

/// `open_buf` with Data's operations per file (PL20).
pub fn open_buf_with(path: &Path, reader: Option<&str>, operations: &Operations) -> Result<PlotSet, String> {
    let text = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let text: String = String::from_utf8(text.clone()).unwrap_or_else(|_| text.iter().map(|&b| b as char).collect());
    let base = path.parent().map(Path::to_path_buf).unwrap_or_default();
    let b = read_buf(&text, &|name| base.join(name).is_file());
    if b.entries.is_empty() {
        return Err(format!("{}: no file named", path.display()));
    }
    open_with(&base, &b.entries, reader, b.title, operations)
}

// PL21: FullProf's temperatures file -----------------------------------------------------------------------------------

/// FullProf's temperatures file of a sequential run (PL21; FullProf's news of 25 September 2013,
/// `manuals/fullprof/news/FullProf_News_2013.htm`): one data file's name per line, the keyword TSAMP and the temperature
/// ("data_23.dat   TSAMP 1.5"); "the temperature is attributed only if the name of the file coincides (including the
/// extension)". The news's own example names data_23.dat twice, with 1.5 and 5.0, and does not say which FullProf takes:
/// a name given twice with different values is listed in `conflicts` and attributes nothing (never a guess, L169).
#[derive(Clone, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TemperatureFile {
    /// (file name, temperature, line from 1), in the file's order.
    pub entries: Vec<(String, f64, usize)>,
    /// The names given twice or more with different values, with every value.
    pub conflicts: Vec<(String, Vec<f64>)>,
    /// The lines that are not "name TSAMP value" (line from 1, text), kept, not dropped.
    pub unread: Vec<(usize, String)>,
}

/// Read a temperatures file (see `TemperatureFile`). The keyword is matched without regard to case.
pub fn read_temperatures(text: &str) -> TemperatureFile {
    let mut t = TemperatureFile::default();
    for (n, line) in text.lines().enumerate() {
        let w: Vec<&str> = line.split_whitespace().collect();
        if w.is_empty() {
            continue;
        }
        match (w.len(), w.get(1).map(|k| k.eq_ignore_ascii_case("TSAMP")), w.get(2).and_then(|v| v.replace(['d', 'D'], "e").parse::<f64>().ok())) {
            (3, Some(true), Some(v)) if v.is_finite() => t.entries.push((w[0].to_string(), v, n + 1)),
            _ => t.unread.push((n + 1, line.to_string())),
        }
    }
    for (name, _, _) in &t.entries {
        let mut values: Vec<f64> = t.entries.iter().filter(|e| &e.0 == name).map(|e| e.1).collect();
        values.dedup();
        if values.len() > 1 && !t.conflicts.iter().any(|c| &c.0 == name) {
            t.conflicts.push((name.clone(), values));
        }
    }
    t
}

impl TemperatureFile {
    /// The temperature of a data file, by its name with its extension (not its folder); None when the file is not
    /// named or is named with different values.
    pub fn value_of(&self, file_name: &str) -> Option<f64> {
        if self.conflicts.iter().any(|c| c.0 == file_name) {
            return None;
        }
        self.entries.iter().find(|e| e.0 == file_name).map(|e| e.1)
    }
}

impl PlotSet {
    /// Each member whose file the temperatures file names (by its name with its extension, as FullProf matches it)
    /// takes that value, its source the temperatures file (`Found::TemperatureFile(name)`, also added to its
    /// candidates); a typed value is kept unless `replace_typed`. The members are sorted again. Returns how many took
    /// a value.
    pub fn variable_from_file(&mut self, name: &str, t: &TemperatureFile, replace_typed: bool) -> usize {
        let mut n = 0;
        for m in self.members.iter_mut() {
            let file_name = Path::new(&m.file).file_name().map_or(m.file.clone(), |f| f.to_string_lossy().into_owned());
            let Some(v) = t.value_of(&file_name) else { continue };
            let found = Variable { value: v, found: Found::TemperatureFile(name.to_string()) };
            if !m.candidates.contains(&found) {
                m.candidates.push(found.clone());
            }
            if !replace_typed && matches!(m.variable.as_ref().map(|v| &v.found), Some(Found::Typed(_))) {
                continue;
            }
            m.variable = Some(found);
            n += 1;
        }
        self.sort_by_variable();
        n
    }
}

impl PlotSet {
    /// The members with a value first, by it (a stable sort: equal values keep their order), then the others.
    pub fn sort_by_variable(&mut self) {
        self.members.sort_by(|a, b| match (&a.variable, &b.variable) {
            (Some(x), Some(y)) => x.value.total_cmp(&y.value),
            (Some(_), None) => std::cmp::Ordering::Less,
            (None, Some(_)) => std::cmp::Ordering::Greater,
            (None, None) => std::cmp::Ordering::Equal,
        });
    }

    /// The user's value for a member, with their note of its source (L169: what a file does not say is typed); the
    /// members sorted again.
    pub fn set_variable(&mut self, member: usize, value: f64, note: &str) -> Result<(), String> {
        let m = self.members.get_mut(member).ok_or(format!("no member {}", member + 1))?;
        if !value.is_finite() {
            return Err(format!("{value} is not a number"));
        }
        m.variable = Some(Variable { value, found: Found::Typed(note.to_string()) });
        self.sort_by_variable();
        Ok(())
    }

    /// The members without a value of the variable (drawn apart until one is typed, L169).
    pub fn without_variable(&self) -> Vec<usize> {
        (0..self.members.len()).filter(|&k| self.members[k].variable.is_none()).collect()
    }
}

/// The y scale a set is shown on (PL3).
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum YScale {
    Linear,
    /// √y (y < 0 shown as −√|y|, so that a difference keeps its sign).
    Sqrt,
    /// log₁₀ y; not a number where y ≤ 0.
    Log,
}

impl YScale {
    pub fn apply(self, y: f64) -> f64 {
        match self {
            YScale::Linear => y,
            YScale::Sqrt => y.signum() * y.abs().sqrt(),
            YScale::Log if y > 0.0 => y.log10(),
            YScale::Log => f64::NAN,
        }
    }
}

/// The rule of a set's displacements (PL3).
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Displacement {
    None,
    /// The k-th member (from 0, in the set's order) moved by k·dx and k·dy: a stack (dx = 0) or a projection.
    Step { dx: f64, dy: f64 },
    /// Moved by (v − v₀)·dx and (v − v₀)·dy, v the member's value of the variable and v₀ the first member's: unequal
    /// steps in T show as unequal gaps. A member without a value is not placed.
    Variable { dx: f64, dy: f64 },
}

/// One member displaced: its x and y on the scale shown, the shift applied, its baseline (the scale's zero moved by
/// the shift, where the stack's label sits) and its background displaced alike.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Displaced {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub shift: (f64, f64),
    pub baseline: f64,
    pub background: Option<Vec<f64>>,
}

/// Each member displaced (PL3), in the set's order: the rule's shift, or the member's own override `(dx, dy)`, added to
/// x and to the scaled y (in √ the shift is added to √y, so that a stack keeps equal gaps on the screen). None for a
/// member the rule cannot place (the `Variable` rule without a value). `backgrounds`, when given, one per member.
pub fn displace(set: &PlotSet, rule: Displacement, scale: YScale, overrides: &[(usize, (f64, f64))], backgrounds: Option<&[Vec<f64>]>) -> Vec<Option<Displaced>> {
    let v0 = set.members.iter().find_map(|m| m.variable.as_ref().map(|v| v.value));
    set.members
        .iter()
        .enumerate()
        .map(|(k, m)| {
            let shift = match overrides.iter().find(|o| o.0 == k) {
                Some(o) => Some(o.1),
                None => match rule {
                    Displacement::None => Some((0.0, 0.0)),
                    Displacement::Step { dx, dy } => Some((k as f64 * dx, k as f64 * dy)),
                    Displacement::Variable { dx, dy } => match (&m.variable, v0) {
                        (Some(v), Some(v0)) => Some(((v.value - v0) * dx, (v.value - v0) * dy)),
                        _ => None,
                    },
                },
            }?;
            let p = &m.pattern;
            Some(Displaced {
                x: p.x.iter().map(|x| x + shift.0).collect(),
                y: p.y.iter().map(|&y| scale.apply(y) + shift.1).collect(),
                shift,
                baseline: shift.1,
                background: backgrounds.and_then(|b| b.get(k)).map(|b| b.iter().map(|&y| scale.apply(y) + shift.1).collect()),
            })
        })
        .collect()
}

// PL4: a common grid and the map ----------------------------------------------------------------------------------------

/// The range of a common grid: where every member has points (the overlap), or where any has (the union).
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum GridRange {
    Overlap,
    Union,
}

/// The step of a common grid: the finest or the coarsest of the members' mean steps, or a given one.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum GridStep {
    Finest,
    Coarsest,
    Given(f64),
}

fn same_x(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-9 * a.abs().max(b.abs()).max(1.0)
}

/// Whether every member has the same points (WinPLOTR's sums need "same X grid and same npts", its guide p. 30).
pub fn same_grid(set: &PlotSet) -> bool {
    let Some(first) = set.members.first() else { return true };
    set.members.iter().all(|m| m.pattern.len() == first.pattern.len() && m.pattern.x.iter().zip(&first.pattern.x).all(|(a, b)| same_x(*a, *b)))
}

/// A common grid for the set (PL4): the members' own points when they all share them and the rule asks for nothing
/// else (the overlap at the finest step); otherwise lo + k·step over the range the rule gives.
pub fn common_grid(set: &PlotSet, range: GridRange, step: GridStep) -> Result<Vec<f64>, String> {
    if set.members.is_empty() {
        return Err("the set has no member".into());
    }
    if same_grid(set) && range == GridRange::Overlap && step == GridStep::Finest {
        return Ok(set.members[0].pattern.x.clone());
    }
    let ranges: Vec<(f64, f64)> = set.members.iter().filter_map(|m| crate::data::extent(&m.pattern.x)).collect();
    let (lo, hi) = match range {
        GridRange::Overlap => (ranges.iter().map(|r| r.0).fold(f64::NEG_INFINITY, f64::max), ranges.iter().map(|r| r.1).fold(f64::INFINITY, f64::min)),
        GridRange::Union => (ranges.iter().map(|r| r.0).fold(f64::INFINITY, f64::min), ranges.iter().map(|r| r.1).fold(f64::NEG_INFINITY, f64::max)),
    };
    if hi.partial_cmp(&lo) != Some(std::cmp::Ordering::Greater) {
        return Err(format!("the members' ranges do not overlap ({lo} to {hi})"));
    }
    let steps: Vec<f64> = set.members.iter().filter_map(|m| m.pattern.mean_step()).map(f64::abs).filter(|s| *s > 0.0).collect();
    let h = match step {
        GridStep::Finest => steps.iter().copied().fold(f64::INFINITY, f64::min),
        GridStep::Coarsest => steps.iter().copied().fold(0.0, f64::max),
        GridStep::Given(h) => h,
    };
    if !(h > 0.0 && h.is_finite()) {
        return Err(format!("the step {h} is not positive"));
    }
    let n = ((hi - lo) / h + 1e-9).floor() as usize + 1;
    Ok((0..n).map(|k| lo + k as f64 * h).collect())
}

/// How a member was put on the grid.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Resampling {
    /// The grid is the member's own points.
    AsMeasured,
    /// Linear interpolation between the two neighbouring points, σ² = (1 − t)²σ₁² + t²σ₂² (`data_ops::interpolate`).
    Interpolated,
    /// The grid is coarser than the member: each grid point the mean of the member's points over its bin
    /// [x − h/2, x + h/2], each point weighted by the part of its own interval (between the midpoints to its
    /// neighbours) inside the bin, σ² = Σwᵢ²σᵢ²/(Σwᵢ)²; counts per point are kept, not summed.
    BinAverage,
}

/// A member on a common grid: y and σ at each grid point, None outside the member's range.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Resampled {
    pub y: Vec<Option<f64>>,
    pub sigma: Vec<Option<f64>>,
    pub how: Resampling,
}

/// One pattern on a grid (PL4).
pub fn resample(p: &Pattern, grid: &[f64]) -> Resampled {
    let n = p.len();
    if n == grid.len() && p.x.iter().zip(grid).all(|(a, b)| same_x(*a, *b)) {
        return Resampled { y: p.y.iter().map(|&v| Some(v)).collect(), sigma: p.sigma.iter().map(|&v| Some(v)).collect(), how: Resampling::AsMeasured };
    }
    let h = if grid.len() > 1 { (grid[grid.len() - 1] - grid[0]) / (grid.len() - 1) as f64 } else { 0.0 };
    let own = p.mean_step().map_or(0.0, f64::abs);
    let (lo, hi) = crate::data::extent(&p.x).unwrap_or((f64::NAN, f64::NAN));
    let inside = |x: f64| x >= lo - 1e-9 * lo.abs().max(1.0) && x <= hi + 1e-9 * hi.abs().max(1.0);
    if h > own * (1.0 + 1e-6) && n > 1 {
        // Each point's interval, between the midpoints to its neighbours.
        let edge = |i: usize| -> f64 {
            if i == 0 {
                p.x[0] - (p.x[1] - p.x[0]) / 2.0
            } else if i == n {
                p.x[n - 1] + (p.x[n - 1] - p.x[n - 2]) / 2.0
            } else {
                (p.x[i - 1] + p.x[i]) / 2.0
            }
        };
        let mut y = Vec::with_capacity(grid.len());
        let mut s = Vec::with_capacity(grid.len());
        for &g in grid {
            if !inside(g) {
                y.push(None);
                s.push(None);
                continue;
            }
            let (a, b) = (g - h / 2.0, g + h / 2.0);
            let i0 = p.x.partition_point(|&v| v < a - 2.0 * own).saturating_sub(1);
            let (mut sw, mut swy, mut sw2s2) = (0.0, 0.0, 0.0);
            for i in i0..n {
                let (ea, eb) = (edge(i), edge(i + 1));
                if ea > b {
                    break;
                }
                let w = (eb.min(b) - ea.max(a)).max(0.0);
                if w > 0.0 && p.sigma[i].is_finite() {
                    sw += w;
                    swy += w * p.y[i];
                    sw2s2 += w * w * p.sigma[i] * p.sigma[i];
                }
            }
            if sw > 0.0 {
                y.push(Some(swy / sw));
                s.push(Some(sw2s2.sqrt() / sw));
            } else {
                y.push(None);
                s.push(None);
            }
        }
        return Resampled { y, sigma: s, how: Resampling::BinAverage };
    }
    let (y, s) = grid.iter().map(|&g| crate::data_ops::interpolate_at(p, g).map_or((None, None), |(v, e)| (Some(v), Some(e)))).unzip();
    Resampled { y, sigma: s, how: Resampling::Interpolated }
}

/// Every member on the grid.
pub fn on_grid(set: &PlotSet, grid: &[f64]) -> Vec<Resampled> {
    set.members.iter().map(|m| resample(&m.pattern, grid)).collect()
}

/// A row of the map: its member, its value of the variable, and its lower and upper edges.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MapRow {
    pub member: usize,
    pub value: f64,
    pub edges: (f64, f64),
    /// The row's member's other variables (W3), each with its source in words, for the map's right-hand axis.
    #[cfg_attr(feature = "serde", serde(default))]
    pub others: Vec<RowOther>,
}

/// A row's other variable (W3): its name, value and source in words (`Found::words`).
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct RowOther {
    pub name: String,
    pub value: f64,
    pub source: String,
    #[cfg_attr(feature = "serde", serde(default))]
    pub unit: String,
}

/// The map of a set (PL4): rows at the variable's values, columns at the grid's points, the values on the scale shown,
/// and the colour range.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Map {
    pub x: Vec<f64>,
    /// The columns' edges: midpoints between the grid's points, half a step out at the ends.
    pub x_edges: Vec<f64>,
    /// The rows in increasing value, each spanning from the midpoint to its lower neighbour to the midpoint to its upper
    /// one (half the gap out at the ends; ±0.5 for a single row).
    pub rows: Vec<MapRow>,
    /// The values, one row per `rows` entry, on the scale shown; None where a member has no point or the scale has no
    /// value (log of y ≤ 0).
    pub z: Vec<Vec<Option<f64>>>,
    /// The members without a value of the variable, drawn apart (L169), and their rows on the grid and the scale, less
    /// the reference for a difference map (PL18), in `apart`'s order.
    pub apart: Vec<usize>,
    #[cfg_attr(feature = "serde", serde(default))]
    pub apart_z: Vec<Vec<Option<f64>>>,
    pub scale: YScale,
    /// The colour range: the smallest and the largest value shown, the rows apart included; for a difference map, −m to
    /// +m with m the largest |difference| (a diverging scale centred on zero).
    pub range: (f64, f64),
    /// The member subtracted, for a difference map.
    pub reference: Option<usize>,
}

/// The map of a set on a grid (PL4); with a reference member, the difference map (each member less the reference, on
/// the linear or √ scale; log₁₀ of a difference has no meaning and is refused). The members without a value of the
/// variable have their rows too, apart (PL18), less the reference alike.
pub fn map(set: &PlotSet, grid: &[f64], scale: YScale, reference: Option<usize>) -> Result<Map, String> {
    if grid.is_empty() {
        return Err("the grid is empty".into());
    }
    if reference.is_some() && scale == YScale::Log {
        return Err("a difference map is drawn on the linear or the √ scale, not on log".into());
    }
    let res = on_grid(set, grid);
    let r = match reference {
        Some(k) => Some(res.get(k).ok_or(format!("no member {}", k + 1))?.y.clone()),
        None => None,
    };
    let rows: Vec<(usize, f64, Vec<Option<f64>>)> = set.members.iter().enumerate().filter_map(|(k, m)| m.variable.as_ref().map(|v| (k, v.value, res[k].y.clone()))).collect();
    let apart: Vec<(usize, Vec<Option<f64>>)> = set.without_variable().into_iter().map(|k| (k, res[k].y.clone())).collect();
    let mut m = map_of_rows(grid, rows, apart, r.as_deref(), scale, reference);
    // Each row's other variables with their sources (W3), for the map's right-hand axis.
    for row in m.rows.iter_mut() {
        row.others = set.members[row.member].others.iter().map(|o| RowOther { name: o.name.clone(), value: o.value, source: o.found.words(), unit: o.unit.clone() }).collect();
    }
    Ok(m)
}

/// A map from rows already on the grid (PL4, PL18, E13a): `rows` (member, value of the variable, y at each of the
/// grid's points), sorted here by value (a stable sort); `apart` the rows without a value; `less` the y subtracted from
/// every row for a difference map. The values on the scale, None where a row has no point or the scale has no value
/// (log of y ≤ 0); the colour range over every value drawn, the rows apart included.
pub(crate) fn map_of_rows(grid: &[f64], mut rows_in: Vec<(usize, f64, Vec<Option<f64>>)>, apart_in: Vec<(usize, Vec<Option<f64>>)>, less: Option<&[Option<f64>]>, scale: YScale, reference: Option<usize>) -> Map {
    rows_in.sort_by(|a, b| a.1.total_cmp(&b.1));
    let n = rows_in.len();
    let rows: Vec<MapRow> = (0..n)
        .map(|i| {
            let v = rows_in[i].1;
            let lo = if i > 0 { (rows_in[i - 1].1 + v) / 2.0 } else if n > 1 { v - (rows_in[1].1 - v) / 2.0 } else { v - 0.5 };
            let hi = if i + 1 < n { (v + rows_in[i + 1].1) / 2.0 } else if n > 1 { v + (v - rows_in[n - 2].1) / 2.0 } else { v + 0.5 };
            MapRow { member: rows_in[i].0, value: v, edges: (lo, hi), others: Vec::new() }
        })
        .collect();
    let on_scale = |y: &[Option<f64>]| -> Vec<Option<f64>> {
        (0..grid.len())
            .map(|j| {
                let y = y.get(j).copied().flatten()?;
                let d = match less {
                    Some(r) => y - r.get(j).copied().flatten()?,
                    None => y,
                };
                let v = scale.apply(d);
                v.is_finite().then_some(v)
            })
            .collect()
    };
    let z: Vec<Vec<Option<f64>>> = rows_in.iter().map(|r| on_scale(&r.2)).collect();
    let apart_z: Vec<Vec<Option<f64>>> = apart_in.iter().map(|r| on_scale(&r.1)).collect();
    let finite = z.iter().chain(&apart_z).flatten().flatten().copied();
    let range = if reference.is_some() {
        let m = finite.fold(0.0_f64, |a, v| a.max(v.abs()));
        (-m, m)
    } else {
        finite.fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), v| (a.min(v), b.max(v)))
    };
    let g = grid;
    let x_edges: Vec<f64> = (0..=g.len())
        .map(|j| {
            if g.len() == 1 {
                g[0] + if j == 0 { -0.5 } else { 0.5 }
            } else if j == 0 {
                g[0] - (g[1] - g[0]) / 2.0
            } else if j == g.len() {
                g[j - 1] + (g[j - 1] - g[j - 2]) / 2.0
            } else {
                (g[j - 1] + g[j]) / 2.0
            }
        })
        .collect();
    Map { x: grid.to_vec(), x_edges, rows, z, apart: apart_in.iter().map(|r| r.0).collect(), apart_z, scale, range, reference }
}

// PL9: derived patterns ---------------------------------------------------------------------------------------------------

/// What a derived pattern is.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Derivation {
    /// Σyᵢ, σ² = Σσᵢ².
    Sum,
    /// Σyᵢ/n, σ = √(Σσᵢ²)/n.
    Average,
    /// Σwᵢyᵢ/Σwᵢ with wᵢ = 1/σᵢ², σ = 1/√Σwᵢ.
    WeightedAverage,
    /// y₁ − y₂ of two members, σ² = σ₁² + σ₂².
    Difference,
}

impl Derivation {
    pub fn word(self) -> &'static str {
        match self {
            Derivation::Sum => "sum",
            Derivation::Average => "average",
            Derivation::WeightedAverage => "weighted average",
            Derivation::Difference => "difference",
        }
    }
}

/// A derived pattern (PL9): x, y, σ on the grid at the points where every member chosen has a point, and what it was
/// made from.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Derived {
    pub derivation: Derivation,
    pub members: Vec<usize>,
    pub pattern: Pattern,
    /// The grid's points left out because a member has none there.
    pub dropped: usize,
}

/// The sum, average or difference of members on a grid (PL9). A difference takes exactly two members, the first less
/// the second. A point where a member has no point, or has one with no weight, is left out (and counted).
pub fn derive(set: &PlotSet, members: &[usize], derivation: Derivation, grid: &[f64]) -> Result<Derived, String> {
    if members.is_empty() {
        return Err("give the members to combine".into());
    }
    if derivation == Derivation::Difference && members.len() != 2 {
        return Err(format!("a difference takes two members, not {}", members.len()));
    }
    if let Some(k) = members.iter().find(|&&k| k >= set.members.len()) {
        return Err(format!("no member {}", k + 1));
    }
    let res: Vec<Resampled> = members.iter().map(|&k| resample(&set.members[k].pattern, grid)).collect();
    let mut p = Pattern::default();
    let mut dropped = 0;
    for (j, &x) in grid.iter().enumerate() {
        let pts: Option<Vec<(f64, f64)>> = res.iter().map(|r| Some((r.y[j]?, r.sigma[j]?)).filter(|v| v.1.is_finite())).collect();
        let Some(pts) = pts else {
            dropped += 1;
            continue;
        };
        let n = pts.len() as f64;
        let (y, s) = match derivation {
            Derivation::Sum => (pts.iter().map(|v| v.0).sum(), pts.iter().map(|v| v.1 * v.1).sum::<f64>().sqrt()),
            Derivation::Average => (pts.iter().map(|v| v.0).sum::<f64>() / n, pts.iter().map(|v| v.1 * v.1).sum::<f64>().sqrt() / n),
            Derivation::WeightedAverage => {
                if pts.iter().any(|v| v.1 <= 0.0) {
                    dropped += 1;
                    continue;
                }
                let w: f64 = pts.iter().map(|v| 1.0 / (v.1 * v.1)).sum();
                (pts.iter().map(|v| v.0 / (v.1 * v.1)).sum::<f64>() / w, 1.0 / w.sqrt())
            }
            Derivation::Difference => (pts[0].0 - pts[1].0, (pts[0].1 * pts[0].1 + pts[1].1 * pts[1].1).sqrt()),
        };
        p.x.push(x);
        p.y.push(y);
        p.sigma.push(s);
    }
    Ok(Derived { derivation, members: members.to_vec(), pattern: p, dropped })
}

/// A derived pattern as FullProf's XYDATA file (`.xys`, Ins = 10: x, y, σ; read back by the registry, negative y kept,
/// which an xye file would lose: GSAS-II reads y ≤ 0 as no weight). The six header lines are those `formats::write`
/// writes (XYDATA, the title, the derivation, the columns); the members follow as '!' lines, each with its file and its
/// value of the variable with the source, then the data. The core's reader skips '!' lines after the header; whether
/// FullProf does is not stated in its manual (not checked).
pub fn derived_text(set: &PlotSet, d: &Derived, title: &str) -> String {
    derived_text_with(set, d, title, &[])
}

/// `derived_text` with lines saying how the members were shown (PL25, `PlotSet::as_shown`), written as '!' lines after
/// the header, before the members.
pub fn derived_text_with(set: &PlotSet, d: &Derived, title: &str, notes: &[String]) -> String {
    let n = d.members.len();
    let mut t = format!(
        "XYDATA\nTITLE: {title} (written by drs {})\n! the {} of {n} member{}{}\n! x = 2theta (deg), y, sigma\n!\n!\n",
        env!("CARGO_PKG_VERSION"),
        d.derivation.word(),
        if n == 1 { "" } else { "s" },
        if d.derivation == Derivation::Difference { ", the first less the second" } else { "" }
    );
    for l in notes {
        t.push_str(&format!("! {l}\n"));
    }
    for (i, &k) in d.members.iter().enumerate() {
        let m = &set.members[k];
        let v = m.variable.as_ref().map_or(format!("{}: nothing says it", set.variable_name), |v| format!("{} = {} {} ({})", set.variable_name, v.value, set.variable_unit, v.found.words()));
        t.push_str(&format!("! member {}: {} (scan {}), {v}\n", i + 1, m.file, m.scan + 1));
    }
    t.push_str(&format!("! {} points{}\n", d.pattern.len(), if d.dropped > 0 { format!(", {} grid points left out where a member has none", d.dropped) } else { String::new() }));
    for i in 0..d.pattern.len() {
        t.push_str(&format!("{} {} {}\n", crate::control::syntax::float(d.pattern.x[i]), crate::control::syntax::float(d.pattern.y[i]), crate::control::syntax::float(d.pattern.sigma[i])));
    }
    t
}

// PL5: windows ------------------------------------------------------------------------------------------------------------

/// A named x window (PL5).
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Window {
    pub name: String,
    pub from: f64,
    pub to: f64,
}

/// The background a window's integral is taken above (PL5).
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum WindowBackground {
    /// Each member's automatic background (`data_ops::auto_background`, GSAS-II's arPLS), at this log₁₀ λ or GSAS-II's
    /// default for its number of points.
    Auto { log_lambda: Option<f64> },
    /// Given points, joined as `background::interpolate` joins them, the same for every member.
    Points { points: Vec<(f64, f64)>, interpolation: Interpolation },
    /// The straight line between the member's own values at the window's two ends (each the linear interpolation of
    /// its two nearest points).
    Line,
}

/// A window's numbers for one member (PL5).
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct WindowValue {
    pub name: String,
    pub from: f64,
    pub to: f64,
    /// The points whose interval (between the midpoints to their neighbours) overlaps the window, and those of them
    /// left out because they have no weight (σ infinite: no count, or y ≤ 0 in an xye file).
    pub points: usize,
    pub unweighted: usize,
    /// ∫ (y − b) dx over [from, to], each point's y − b taken over the part of its interval inside the window, with σ
    /// from the points' σ through every y the sum depends on (the line's two ends included); the automatic background's
    /// and given points' own uncertainty left out.
    pub integral: f64,
    pub sigma: f64,
    /// Σ x (y − b) Δx / Σ (y − b) Δx with its σ (propagated as the integral's), when the integral is positive.
    pub centroid: Option<f64>,
    pub centroid_sigma: Option<f64>,
    /// The strongest peak in the window (by prominence) of the member's `peak_search` at the proposed prominence (a
    /// fiftieth of the largest count) and the dictionary's significance (3σ): its position, height and width.
    pub peak: Option<FoundPeak>,
    /// Whether the integral stands above `significance` times its σ (PL22): a value below it is drawn open, as the
    /// design's F5 draws it.
    pub significant: bool,
    pub significance: f64,
}

/// The integral of a window counts as measured when it stands above this many σ (PL22; F5's "open: below 2σ"): the
/// usual two-sided 95 % level for a single value. The threshold is the core's, as `lsq::strong_correlations`'s 0.8.
pub const WINDOW_SIGNIFICANCE: f64 = 2.0;

/// The windows of one member, and its background at the chosen x (PL5).
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MemberWindows {
    pub member: usize,
    pub windows: Vec<WindowValue>,
    /// The background at each chosen x: None outside the member's range (automatic) and for the line, which a window
    /// defines only inside itself.
    pub background_at: Vec<Option<f64>>,
}

/// The windows of a set (PL5): the series view's integrals, centroids and peaks against the variable.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Windows {
    pub background: WindowBackground,
    pub at: Vec<f64>,
    pub members: Vec<MemberWindows>,
}

/// The two points x lies between and its fraction of the way from the first: y(x) = (1 − u) y_k + u y_{k+1}.
fn bracket(x: &[f64], t: f64) -> Option<(usize, f64)> {
    let n = x.len();
    if n < 2 || !(t >= x[0] && t <= x[n - 1]) {
        return None;
    }
    let k = x.partition_point(|&v| v <= t).clamp(1, n - 1) - 1;
    Some((k, (t - x[k]) / (x[k + 1] - x[k])))
}

/// One window of one pattern (PL5): `b` the background at each point (unused for the line).
pub fn window_value(p: &Pattern, b: &[f64], line: bool, w: &Window, peaks: &[FoundPeak]) -> Result<WindowValue, String> {
    let n = p.len();
    if n < 2 || p.x.windows(2).any(|v| v[1] <= v[0]) {
        return Err("a window needs two points or more, x increasing".into());
    }
    if w.to.partial_cmp(&w.from) != Some(std::cmp::Ordering::Greater) {
        return Err(format!("window {}: from {} to {} is empty", w.name, w.from, w.to));
    }
    let x = &p.x;
    let edge = |i: usize| -> f64 {
        if i == 0 {
            x[0] - (x[1] - x[0]) / 2.0
        } else if i == n {
            x[n - 1] + (x[n - 1] - x[n - 2]) / 2.0
        } else {
            (x[i - 1] + x[i]) / 2.0
        }
    };
    // The line through the member's values at the two ends, each from its two nearest points.
    let ends = if line {
        let a = bracket(x, w.from).ok_or(format!("window {}: {} lies outside the data", w.name, w.from))?;
        let c = bracket(x, w.to).ok_or(format!("window {}: {} lies outside the data", w.name, w.to))?;
        let v = |(k, u): (usize, f64)| (1.0 - u) * p.y[k] + u * p.y[k + 1];
        Some((a, c, v(a), v(c)))
    } else {
        None
    };
    // d(integral)/dy_j and d(moment)/dy_j, for σ.
    let (mut a, mut m) = (vec![0.0; n], vec![0.0; n]);
    let (mut integral, mut moment, mut points, mut unweighted) = (0.0, 0.0, 0, 0);
    let mut s = [0.0; 4]; // Σ Δx (1 − t), Σ Δx t, Σ Δx x (1 − t), Σ Δx x t over the window's points, t along the line.
    for i in 0..n {
        let dx = (edge(i + 1).min(w.to) - edge(i).max(w.from)).max(0.0);
        if dx <= 0.0 {
            continue;
        }
        points += 1;
        if !p.sigma[i].is_finite() {
            unweighted += 1;
            continue;
        }
        let t = (x[i] - w.from) / (w.to - w.from);
        let bi = match ends {
            Some((_, _, ya, yb)) => ya + (yb - ya) * t,
            None => b[i],
        };
        integral += dx * (p.y[i] - bi);
        moment += dx * x[i] * (p.y[i] - bi);
        a[i] += dx;
        m[i] += dx * x[i];
        s[0] += dx * (1.0 - t);
        s[1] += dx * t;
        s[2] += dx * x[i] * (1.0 - t);
        s[3] += dx * x[i] * t;
    }
    if let Some(((ka, ua), (kb, ub), _, _)) = ends {
        for (k, u, si, sm) in [(ka, ua, s[0], s[2]), (kb, ub, s[1], s[3])] {
            a[k] -= si * (1.0 - u);
            a[k + 1] -= si * u;
            m[k] -= sm * (1.0 - u);
            m[k + 1] -= sm * u;
        }
    }
    let spread = |c: &dyn Fn(usize) -> f64| -> f64 { (0..n).filter(|&j| c(j) != 0.0).map(|j| (c(j) * p.sigma[j]).powi(2)).sum::<f64>().sqrt() };
    let sigma = spread(&|j| a[j]);
    let (centroid, centroid_sigma) = if integral > 0.0 {
        let c = moment / integral;
        (Some(c), Some(spread(&|j| (m[j] - c * a[j]) / integral)))
    } else {
        (None, None)
    };
    let peak = peaks.iter().filter(|q| q.x >= w.from && q.x <= w.to).max_by(|q, r| q.prominence.total_cmp(&r.prominence)).cloned();
    let significant = integral > WINDOW_SIGNIFICANCE * sigma;
    Ok(WindowValue { name: w.name.clone(), from: w.from, to: w.to, points, unweighted, integral, sigma, centroid, centroid_sigma, peak, significant, significance: WINDOW_SIGNIFICANCE })
}

/// The windows of every member of a set (PL5), each on its own points, and the background at the chosen x.
pub fn windows(set: &PlotSet, list: &[Window], background: &WindowBackground, at: &[f64]) -> Result<Windows, String> {
    let mut members = Vec::new();
    for (k, mb) in set.members.iter().enumerate() {
        let p = &mb.pattern;
        let peaks = crate::data_ops::peak_search(p, crate::data::proposed_prominence(p).unwrap_or(0.0), crate::data_ops::DEFAULT_MIN_SIGNIFICANCE);
        let (b, background_at): (Vec<f64>, Vec<Option<f64>>) = match background {
            WindowBackground::Auto { log_lambda } => {
                let curve = crate::data_ops::auto_background(p, *log_lambda, crate::data_ops::DEFAULT_MAX_ITER).map_err(|e| format!("{}: {e}", mb.label))?.curve;
                let at_x = at.iter().map(|&t| bracket(&p.x, t).map(|(i, u)| (1.0 - u) * curve[i] + u * curve[i + 1])).collect();
                (curve, at_x)
            }
            WindowBackground::Points { points, interpolation } => {
                (crate::background::interpolate(points, *interpolation, &p.x), crate::background::interpolate(points, *interpolation, at).into_iter().map(Some).collect())
            }
            WindowBackground::Line => (Vec::new(), vec![None; at.len()]),
        };
        let line = matches!(background, WindowBackground::Line);
        let windows = list.iter().map(|w| window_value(p, &b, line, w, &peaks)).collect::<Result<Vec<_>, _>>().map_err(|e| format!("{}: {e}", mb.label))?;
        members.push(MemberWindows { member: k, windows, background_at });
    }
    Ok(Windows { background: background.clone(), at: at.to_vec(), members })
}

// PL17: the projection's hidden parts ----------------------------------------------------------------------------------

/// One member's curve in the projection (PL17): the pieces of it that the curves in front leave visible, each a
/// polyline in the displaced coordinates, and the share of its x length they cover.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Projected {
    pub member: usize,
    pub pieces: Vec<Vec<(f64, f64)>>,
    pub visible: f64,
}

/// The projection with hidden parts (PL17): the curves front to back, each with its visible pieces, and the oblique
/// axis: each member's last x at the foot of its curtain (the floor raised by its shift, or its baseline).
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Projection {
    pub order: Vec<usize>,
    pub floor: Option<f64>,
    pub curves: Vec<Projected>,
    pub axis: Vec<(usize, (f64, f64))>,
}

struct Curtain {
    x: Vec<f64>,
    y: Vec<f64>,
    bottom: f64,
}

impl Curtain {
    fn at(&self, t: f64) -> f64 {
        let k = self.x.partition_point(|&v| v <= t).clamp(1, self.x.len() - 1) - 1;
        self.y[k] + (self.y[k + 1] - self.y[k]) * (t - self.x[k]) / (self.x[k + 1] - self.x[k])
    }
}

/// The part of [lo, hi] where g ≤ 0, g linear on [u, v] with g(u) = gu, g(v) = gv.
fn where_not_above(lo: f64, hi: f64, u: f64, v: f64, gu: f64, gv: f64) -> (f64, f64) {
    match (gu <= 0.0, gv <= 0.0) {
        (true, true) => (lo, hi),
        (false, false) => (hi, lo),
        (true, false) => (lo, hi.min(u + (v - u) * gu / (gu - gv))),
        (false, true) => (lo.max(u + (v - u) * gu / (gu - gv)), hi),
    }
}

/// The x intervals of the segment (xa, ya)–(xb, yb) inside a curtain: within its x range, at or below its curve and at
/// or above its bottom.
fn hidden_by(xa: f64, ya: f64, xb: f64, yb: f64, f: &Curtain, out: &mut Vec<(f64, f64)>) {
    let (lo, hi) = (xa.max(f.x[0]), xb.min(f.x[f.x.len() - 1]));
    if hi <= lo {
        return;
    }
    let line = |t: f64| ya + (yb - ya) * (t - xa) / (xb - xa);
    let mut cuts = vec![lo];
    cuts.extend(f.x[f.x.partition_point(|&v| v <= lo)..].iter().take_while(|&&v| v < hi));
    cuts.push(hi);
    for c in cuts.windows(2) {
        let (u, v) = (c[0], c[1]);
        if v <= u {
            continue;
        }
        let (p, q) = where_not_above(u, v, u, v, line(u) - f.at(u), line(v) - f.at(v));
        let (p, q) = if f.bottom.is_finite() && q > p { where_not_above(p, q, u, v, f.bottom - line(u), f.bottom - line(v)) } else { (p, q) };
        if q > p {
            out.push((p, q));
        }
    }
}

/// The projection with hidden parts (PL17): the members in `order`, front to back, each displaced (`displace`); each
/// curve hides what lies under it, down to its floor (`floor` raised by the member's y shift) or without end when no
/// floor is given, across its own x range. Returned: the visible pieces of each curve, which the app draws as they are.
pub fn project(displaced: &[Option<Displaced>], order: &[usize], floor: Option<f64>) -> Result<Projection, String> {
    let mut seen = vec![false; displaced.len()];
    let mut front: Vec<Curtain> = Vec::new();
    let (mut curves, mut axis) = (Vec::new(), Vec::new());
    for &k in order {
        let d = displaced.get(k).ok_or(format!("no member {}", k + 1))?.as_ref().ok_or(format!("member {} is not placed (no value of the variable)", k + 1))?;
        if std::mem::replace(&mut seen[k], true) {
            return Err(format!("member {} is twice in the order", k + 1));
        }
        let (mut x, mut y) = (d.x.clone(), d.y.clone());
        if x.len() > 1 && x[0] > x[x.len() - 1] {
            x.reverse();
            y.reverse();
        }
        if x.len() < 2 || x.windows(2).any(|v| v[1] <= v[0]) || y.iter().any(|v| !v.is_finite()) {
            return Err(format!("member {}: the projection needs x strictly monotonic and every y a number on the scale shown", k + 1));
        }
        let bottom = floor.map_or(f64::NEG_INFINITY, |f| f + d.shift.1);
        let (mut pieces, mut seen_len): (Vec<Vec<(f64, f64)>>, f64) = (Vec::new(), 0.0);
        for i in 0..x.len() - 1 {
            let (xa, ya, xb, yb) = (x[i], y[i], x[i + 1], y[i + 1]);
            let mut hid = Vec::new();
            for f in &front {
                hidden_by(xa, ya, xb, yb, f, &mut hid);
            }
            hid.sort_by(|p, q| p.0.total_cmp(&q.0));
            let eps = 1e-9 * (xb - xa);
            let mut start = xa;
            let mut visible = Vec::new();
            for (p, q) in hid {
                if p > start + eps {
                    visible.push((start, p));
                }
                start = start.max(q);
            }
            if xb > start + eps {
                visible.push((start, xb));
            }
            let line = |t: f64| ya + (yb - ya) * (t - xa) / (xb - xa);
            for (p, q) in visible {
                seen_len += q - p;
                let (a, b) = ((p, if p == xa { ya } else { line(p) }), (q, if q == xb { yb } else { line(q) }));
                match pieces.last_mut() {
                    Some(last) if last.last().is_some_and(|e| e.0 == p) => last.push(b),
                    _ => pieces.push(vec![a, b]),
                }
            }
        }
        axis.push((k, (x[x.len() - 1], if floor.is_some() { bottom } else { d.baseline })));
        curves.push(Projected { member: k, pieces, visible: seen_len / (x[x.len() - 1] - x[0]) });
        front.push(Curtain { x, y, bottom });
    }
    Ok(Projection { order: order.to_vec(), floor, curves, axis })
}

// PL19: a proposed step for a stack -----------------------------------------------------------------------------------

/// A step for a stack proposed from the data (PL19), on the scale shown, with its rule and each member's spread.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ProposedStep {
    pub dy: f64,
    pub scale: YScale,
    pub rule: String,
    /// Each member's spread on the scale: its 90th less its 10th percentile (None without two values on the scale).
    pub spreads: Vec<Option<f64>>,
}

/// The value at quantile q of sorted values, by nearest rank: v[round(q (n − 1))].
fn quantile(sorted: &[f64], q: f64) -> f64 {
    sorted[((q * (sorted.len() - 1) as f64).round() as usize).min(sorted.len() - 1)]
}

/// The smallest of 1, 2 and 5 × 10ᵏ at or above v (a step the axis labels well).
fn round_up_125(v: f64) -> f64 {
    let p = 10f64.powf(v.log10().floor());
    [1.0, 2.0, 5.0, 10.0].iter().map(|m| m * p).find(|&s| s >= v * (1.0 - 1e-12)).unwrap_or(10.0 * p)
}

/// A step for a stack (PL19): the median over the members of their spread on the scale shown (the 90th less the 10th
/// percentile of the scaled y, nearest rank), rounded up to 1, 2 or 5 × 10ᵏ. A member's bulk (its background and the
/// feet of its peaks) then stays clear of the next one's while its strong peaks may reach above it, as the mock-up F1
/// and WinPLOTR's stacks show (its guide gives no rule). None without a member with two values on the scale.
pub fn proposed_step(set: &PlotSet, scale: YScale) -> Option<ProposedStep> {
    let spreads: Vec<Option<f64>> = set
        .members
        .iter()
        .map(|m| {
            let mut v: Vec<f64> = m.pattern.y.iter().map(|&y| scale.apply(y)).filter(|v| v.is_finite()).collect();
            if v.len() < 2 {
                return None;
            }
            v.sort_by(f64::total_cmp);
            Some(quantile(&v, 0.9) - quantile(&v, 0.1))
        })
        .collect();
    let mut known: Vec<f64> = spreads.iter().flatten().copied().filter(|s| *s > 0.0).collect();
    if known.is_empty() {
        return None;
    }
    known.sort_by(f64::total_cmp);
    let n = known.len();
    let median = if n % 2 == 1 { known[n / 2] } else { (known[n / 2 - 1] + known[n / 2]) / 2.0 };
    let dy = round_up_125(median);
    let rule = format!("the median over {n} members of their 90th less 10th percentile on the {} scale ({}), rounded up to 1, 2 or 5 × 10^k", scale.word(), crate::control::syntax::float(median));
    Some(ProposedStep { dy, scale, rule, spreads })
}

impl YScale {
    /// The scale's word in the plot set file.
    pub fn word(self) -> &'static str {
        match self {
            YScale::Linear => "linear",
            YScale::Sqrt => "sqrt",
            YScale::Log => "log",
        }
    }

    pub fn from_word(w: &str) -> Option<YScale> {
        match w {
            "linear" => Some(YScale::Linear),
            "sqrt" => Some(YScale::Sqrt),
            "log" => Some(YScale::Log),
            _ => None,
        }
    }
}

// PL12: the plot set file ---------------------------------------------------------------------------------------------

/// The format line of a plot set file.
pub const SET_FORMAT: &str = "drs plot set 1";

/// The members' order in a set file (PL12): by the variable (those without one last), by file name, or as given.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Order {
    #[default]
    Variable,
    Name,
    Given,
}

/// One `[[member]]` of a set file: the file (relative to the project's folder), the scan (from 1; every scan of the
/// file when not given), a value of the variable the user chose with its source (typed, or read from a temperatures
/// file, PL23), and how it is drawn. What is not written is
/// the core's default: the variable as read from the file, visible, opacity 1, the set's palette.
#[derive(Clone, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SetMember {
    pub file: String,
    pub scan: Option<usize>,
    /// The value and its source: `Found::Typed` (`source = "…"`) or `Found::TemperatureFile` (`temperatures = "…"`),
    /// the two a file read again does not give.
    pub value: Option<(f64, Found)>,
    pub visible: Option<bool>,
    pub opacity: Option<f64>,
    pub colour: Option<String>,
    /// Values of other variables than the set's, each typed with its source (W3): `variables = { T = { value = 1270,
    /// source = "…" } }`, in the file's order.
    #[cfg_attr(feature = "serde", serde(default))]
    pub variables: Vec<(String, f64, String)>,
}

/// A plot set as a file of the project (PL12, `plots/DESIGN.md` §2.2): `Plots/{set}.drv`, the control file's TOML
/// subset. The x quantity is 2θ until the conversion of a set (PL10).
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SetFile {
    pub title: Option<String>,
    pub variable_name: String,
    pub variable_unit: String,
    pub x: String,
    pub y_scale: YScale,
    pub order: Order,
    /// The stack's displacement, and the y scale its step was set on (PL27; none: `y_scale`).
    pub displacement: Displacement,
    #[cfg_attr(feature = "serde", serde(default))]
    pub displacement_scale: Option<YScale>,
    /// The projection's displacement apart from the stack's, and the y scale its steps were set on (PL27).
    #[cfg_attr(feature = "serde", serde(default = "no_displacement"))]
    pub projection: Displacement,
    #[cfg_attr(feature = "serde", serde(default))]
    pub projection_scale: Option<YScale>,
    /// The map's colour scale apart from the curves' y scale (PL27; none: `y_scale`).
    #[cfg_attr(feature = "serde", serde(default))]
    pub map_scale: Option<YScale>,
    /// The curves' and the row strip's y interval (PL41, W4): fitted to the whole set (the default) or to the member
    /// chosen or shown.
    #[cfg_attr(feature = "serde", serde(default))]
    pub y_interval: YInterval,
    /// The members rebinned on opening to this step in their x (W9, a large set: `rebin = 0.01` in 2θ), the counts
    /// averaged over each new bin and σ propagated, as `drs data op rebin step=… mode=average`; none by default.
    #[cfg_attr(feature = "serde", serde(default))]
    pub rebin: Option<f64>,
    /// The view's normalisation (PL6, PL25, PL27); none by default.
    #[cfg_attr(feature = "serde", serde(default))]
    pub normalisation: Normalisation,
    /// The reference member's file (differences, the difference map).
    pub reference: Option<String>,
    pub windows: Vec<Window>,
    /// The map's colour range as the user clipped it (PL37, I08), on the map's scale (`map_scale`, else `y_scale`):
    /// values beyond are drawn saturated at the end colours; none: the core's own range (`Map::range`). A view, never an
    /// operation on the data.
    #[cfg_attr(feature = "serde", serde(default))]
    pub map_range: Option<Range>,
    /// The difference map's, apart (its range is −m to m on its own scale).
    #[cfg_attr(feature = "serde", serde(default))]
    pub difference_range: Option<Range>,
    /// The regions of the variable–x plane saved by name (PL37), `[[regions]]`.
    #[cfg_attr(feature = "serde", serde(default))]
    pub regions: Vec<Region>,
    pub members: Vec<SetMember>,
}

/// The y interval of the curves and the row strip (PL41): the whole set's (`plots::extent`), or each member's own.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(rename_all = "lowercase"))]
pub enum YInterval {
    #[default]
    Shared,
    Own,
}

impl YInterval {
    pub fn word(self) -> &'static str {
        match self {
            YInterval::Shared => "shared",
            YInterval::Own => "own",
        }
    }
}

/// A colour range chosen by the user (PL37): `from` below `to`, on the scale the map is drawn on.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Range {
    pub from: f64,
    pub to: f64,
}

/// A region of the map saved by name (PL37, I08): `x` from its smaller end in the quantity it was drawn in (none: the
/// set's `x`), `variable` from its smaller end in the variable's unit, a colour range of its own when "range from the
/// view" was on, and the window it was drawn from, if any (a name, nothing more: a region is a view of the map, a window
/// a 2θ interval whose quantities are computed per member).
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Region {
    pub name: String,
    pub x: [f64; 2],
    #[cfg_attr(feature = "serde", serde(default))]
    pub x_quantity: Option<XQuantity>,
    pub variable: [f64; 2],
    #[cfg_attr(feature = "serde", serde(default))]
    pub map_range: Option<Range>,
    /// Its own colour range on the difference map, on that map's scale (P65): a region saved on the difference map with
    /// "range from the view" on keeps it here, one saved on the map in `map_range`, so that each range is applied to its
    /// own view only.
    #[cfg_attr(feature = "serde", serde(default))]
    pub difference_range: Option<Range>,
    #[cfg_attr(feature = "serde", serde(default))]
    pub window: Option<String>,
}

/// The x a map is drawn in (PL37): its quantity (none: TOF or energy, not converted) and, for each member as the map's
/// rows name it, its λ (the one the view used, else its file's) and its label, so that a region's bounds are converted
/// row by row as the view converted the member.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MapX {
    pub quantity: Option<XQuantity>,
    pub lambda: Vec<Option<f64>>,
    pub labels: Vec<String>,
}

impl MapX {
    /// A set drawn as its files give it: the quantity of its members' axis (none for TOF, energy, or members in
    /// different quantities), each member's λ from its file.
    pub fn of(set: &PlotSet) -> MapX {
        let q = |a: Axis| match a {
            Axis::TwoTheta => Some(XQuantity::TwoTheta),
            Axis::D => Some(XQuantity::D),
            Axis::Q => Some(XQuantity::Q),
            _ => None,
        };
        let first = set.members.first().and_then(|m| q(m.facts.axis));
        let quantity = first.filter(|f| set.members.iter().all(|m| q(m.facts.axis) == Some(*f)));
        MapX { quantity, lambda: set.members.iter().map(|m| m.facts.wavelength).collect(), labels: set.members.iter().map(|m| m.label.clone()).collect() }
    }

    /// A set shown in another quantity (PL10, `x_view`): the view's quantity, each member's λ the view's, else its file's.
    pub fn of_view(set: &PlotSet, view: &XView) -> MapX {
        let lambda = set.members.iter().enumerate().map(|(k, m)| view.members.iter().find(|x| x.member == k).and_then(|x| x.lambda).or(m.facts.wavelength)).collect();
        MapX { quantity: Some(view.quantity), lambda, labels: set.members.iter().map(|m| m.label.clone()).collect() }
    }
}

/// The block of a map inside a region (PL37): the rows inside (indices into `Map::rows`, the row's value of the
/// variable within the region's), the columns inside (indices into `Map::x`, the column's x within the region's bounds
/// converted to the map's quantity with each row's λ; their union when the rows' λ differ), those bounds (their union),
/// and the smallest and largest value of the block's cells, each row's own columns only.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MapBlock {
    pub rows: Vec<usize>,
    pub columns: Vec<usize>,
    pub x: (f64, f64),
    pub range: (f64, f64),
    /// Where `rows` and `columns` count from: always 0, the map's arrays' own indices (PL38). A person's count, as the
    /// page and `drs plots region` print it, is one more.
    #[cfg_attr(feature = "serde", serde(default))]
    pub counted_from: usize,
}

/// A value of x in sin θ/λ (PL10's route): from 2θ sin(θ)/λ, d 1/(2d), Q Q/4π, s s/2.
fn stl_of(q: XQuantity, v: f64, lambda: Option<f64>) -> Option<f64> {
    let t = match q {
        XQuantity::TwoTheta => (v / 2.0).to_radians().sin() / lambda?,
        XQuantity::D => 1.0 / (2.0 * v),
        XQuantity::Q => v / (4.0 * std::f64::consts::PI),
        XQuantity::Stl => v,
        XQuantity::S => v / 2.0,
    };
    t.is_finite().then_some(t)
}

/// A value of sin θ/λ in a quantity (`x_view`'s): d 1/(2t), Q 4πt, s 2t, 2θ 2 asin(λt); none where it has no value.
fn of_stl(q: XQuantity, t: f64, lambda: Option<f64>) -> Option<f64> {
    let v = match q {
        XQuantity::Stl => t,
        XQuantity::S => 2.0 * t,
        XQuantity::Q => 4.0 * std::f64::consts::PI * t,
        XQuantity::D => (t > 0.0).then(|| 1.0 / (2.0 * t))?,
        XQuantity::TwoTheta => {
            let a = lambda? * t;
            (a.abs() <= 1.0).then(|| 2.0 * a.asin().to_degrees())?
        }
    };
    v.is_finite().then_some(v)
}

/// Whether a region drawn in `from` needs a λ to be shown in `to`: when the two differ and one is 2θ.
fn needs_lambda(from: XQuantity, to: XQuantity) -> bool {
    from != to && (from == XQuantity::TwoTheta || to == XQuantity::TwoTheta)
}

/// A region's x bounds in another quantity, with a λ (the smaller end first).
fn bounds_in(x: [f64; 2], from: XQuantity, to: XQuantity, lambda: Option<f64>) -> Result<(f64, f64), String> {
    if from == to {
        return Ok((x[0], x[1]));
    }
    let one = |v: f64| stl_of(from, v, lambda).and_then(|t| of_stl(to, t, lambda)).ok_or(format!("{} {v} has no value in {}", from.word(), to.word()));
    let (a, b) = (one(x[0])?, one(x[1])?);
    Ok((a.min(b), a.max(b)))
}

/// The block of a map inside a region (PL37, I08), the map drawn in `x`: a reading of the matrix the core computed, so
/// that "range from the view" and the zoom take the block and its range from the core (L139, PL14). Refused, with the
/// members named, when the region's quantity is not the map's and the conversion needs a λ some members lack; when the
/// region holds no row, no column, or no value. The block is the same whatever colour ranges the region keeps
/// (`map_range`, `difference_range`, P65): they are views of the block, applied by the host to their own map.
pub fn map_region(map: &Map, region: &Region, x: &MapX) -> Result<MapBlock, String> {
    let name = &region.name;
    let to = x.quantity;
    let from = region.x_quantity.or(to);
    let (from, to) = match (from, to) {
        (Some(f), Some(t)) => (f, t),
        (None, None) => (XQuantity::TwoTheta, XQuantity::TwoTheta),
        (Some(f), None) => return Err(format!("region \"{name}\" is in {}, and the map's x is not converted (TOF or energy)", f.word())),
        (None, Some(t)) => (t, t),
    };
    let rows: Vec<usize> = (0..map.rows.len()).filter(|&i| (region.variable[0]..=region.variable[1]).contains(&map.rows[i].value)).collect();
    if rows.is_empty() {
        return Err(format!("region \"{name}\": no row of the map has its value of the variable from {} to {}", region.variable[0], region.variable[1]));
    }
    if needs_lambda(from, to) {
        let lacking: Vec<String> = rows.iter().map(|&i| map.rows[i].member).filter(|&k| x.lambda.get(k).copied().flatten().is_none()).map(|k| x.labels.get(k).cloned().unwrap_or(format!("member {}", k + 1))).collect();
        if !lacking.is_empty() {
            return Err(format!("region \"{name}\" is in {} and the map in {}: converting it needs each member's wavelength, and these state none: {}", from.word(), to.word(), lacking.join(", ")));
        }
    }
    let mut columns: Vec<usize> = Vec::new();
    let mut lo_hi = (f64::INFINITY, f64::NEG_INFINITY);
    let mut range = (f64::INFINITY, f64::NEG_INFINITY);
    for &i in &rows {
        let lambda = x.lambda.get(map.rows[i].member).copied().flatten();
        let (lo, hi) = bounds_in(region.x, from, to, lambda).map_err(|e| format!("region \"{name}\": {e}"))?;
        lo_hi = (lo_hi.0.min(lo), lo_hi.1.max(hi));
        for (j, &c) in map.x.iter().enumerate() {
            if c >= lo && c <= hi {
                if !columns.contains(&j) {
                    columns.push(j);
                }
                if let Some(v) = map.z[i].get(j).copied().flatten() {
                    range = (range.0.min(v), range.1.max(v));
                }
            }
        }
    }
    columns.sort_unstable();
    if columns.is_empty() {
        return Err(format!("region \"{name}\": no column of the map lies from {} to {} in {}", lo_hi.0, lo_hi.1, to.word()));
    }
    if range.0 > range.1 {
        return Err(format!("region \"{name}\": its cells hold no value on the map's scale"));
    }
    Ok(MapBlock { rows, columns, x: lo_hi, range, counted_from: 0 })
}

/// A colour range checked against a map's values (PL37): a warning when it lies wholly outside them, every cell then
/// drawn saturated at one end; none otherwise.
pub fn range_warning(what: &str, r: &Range, map: &Map) -> Option<String> {
    (r.to < map.range.0 || r.from > map.range.1).then(|| format!("{what} {} to {} lies wholly outside the map's values, {} to {}: every cell is drawn saturated", r.from, r.to, map.range.0, map.range.1))
}

// PL38: the page's data -------------------------------------------------------------------------------------------------

/// The data of the page without the app (PL38, I08 step 5): everything the Map viewer draws, in the shape its header
/// names (`app/ui/viewer/map.js`, the first comment, the contract). The page computes nothing: every number is the
/// core's, and a region's block is the one written here.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PageData {
    pub title: Option<String>,
    pub path: String,
    /// The core's version and the day the page was written (UTC, YYYY-MM-DD).
    pub core: String,
    pub written: String,
    pub variable: PageVariable,
    pub x: PageX,
    /// "counts", or "normalised counts" when the set file names a normalisation.
    pub y_word: String,
    pub members: Vec<PageMember>,
    pub maps: PageMaps,
    pub clip: PageClip,
    pub regions: Vec<PageRegion>,
    pub regions_refused: Option<String>,
    pub same_grid: bool,
    /// The view the page opens on, "map" or "difference" (`drs plots show --view`); the page's shell of Plots' step 18
    /// opens on the map and does not yet read it.
    pub view: String,
}

#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PageVariable {
    pub name: String,
    pub unit: String,
}

/// The x the map is drawn in: the label the app shows ("2θ (°)", "d (Å)") and its word (none for TOF, energy or mixed).
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PageX {
    pub label: String,
    pub quantity: Option<String>,
}

/// A member by its index in the set, as the map's rows name it: its label, its value of the variable with the source
/// in the core's words, and its pattern as the set holds it (left out with `--no-patterns`).
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PageMember {
    pub label: String,
    pub variable: Option<PageValue>,
    pub pattern: Option<PagePattern>,
}

#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PageValue {
    pub value: f64,
    pub source: String,
}

#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PagePattern {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
}

#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PageMaps {
    pub map: Option<PageMap>,
    pub difference: Option<PageMap>,
}

/// A map as the app's `plots_map` gives it to the viewer: the cells, the rows apart as `{ member, z }`, the scale as its
/// word, the reference member.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PageMap {
    pub map: PageCells,
    pub apart: Vec<PageApart>,
    pub scale: String,
    pub reference: Option<usize>,
}

#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PageCells {
    pub x: Vec<f64>,
    pub x_edges: Vec<f64>,
    pub rows: Vec<MapRow>,
    pub z: Vec<Vec<Option<f64>>>,
    pub range: (f64, f64),
}

#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PageApart {
    pub member: usize,
    pub z: Option<Vec<Option<f64>>>,
}

#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PageClip {
    pub map: Option<Range>,
    pub difference: Option<Range>,
}

/// A saved region with its blocks of both maps precomputed (the page has no core to ask), or why it has none.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PageRegion {
    pub name: String,
    pub x: [f64; 2],
    pub x_quantity: Option<String>,
    pub variable: [f64; 2],
    pub map_range: Option<Range>,
    /// The region's own range on the difference map (P65), beside `map_range` on the map.
    pub difference_range: Option<Range>,
    pub window: Option<String>,
    pub blocks: PageBlocks,
    pub refused: Option<String>,
}

#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PageBlocks {
    pub map: Option<MapBlock>,
    pub difference: Option<MapBlock>,
}

/// The scale the difference map is drawn on for a map scale (the module's rule, `diffScaleOf`): the square root for a
/// log map scale, since log₁₀ of a difference has no meaning (`map` refuses it); the map's own otherwise.
pub fn difference_scale(map_scale: YScale) -> YScale {
    if map_scale == YScale::Log {
        YScale::Sqrt
    } else {
        map_scale
    }
}

fn page_map(m: &Map) -> PageMap {
    PageMap {
        map: PageCells { x: m.x.clone(), x_edges: m.x_edges.clone(), rows: m.rows.clone(), z: m.z.clone(), range: m.range },
        apart: m.apart.iter().zip(&m.apart_z).map(|(&member, z)| PageApart { member, z: Some(z.clone()) }).collect(),
        scale: m.scale.word().to_string(),
        reference: m.reference,
    }
}

/// The page's data (PL38): `set` as the view shows it (in the set file's x, normalised as the file says), drawn in `x`
/// (`MapX::of`, or `MapX::of_view` for a set shown in another quantity), on the common grid (the overlap at the finest
/// step, as `drs plots map`); the map on the file's map scale (else its y scale) and, with `reference`, the difference
/// map on `difference_scale`; the file's clip and regions, each region's blocks by `map_region` or its refusal in the
/// core's words; `path` the file's name as the page says it; the members' patterns only with `patterns`. Refused when
/// fewer than two members have a value of the variable: a map has one row per value.
pub fn page_data(set: &PlotSet, file: &SetFile, x: &MapX, reference: Option<usize>, path: &str, patterns: bool) -> Result<PageData, String> {
    let with_value: Vec<&str> = set.members.iter().filter(|m| m.variable.is_some()).map(|m| m.label.as_str()).collect();
    if with_value.len() < 2 {
        return Err(format!(
            "a map has one row per value of {}, and this set has {} member{} with one{}: no map, so no page",
            set.variable_name,
            with_value.len(),
            if with_value.len() == 1 { "" } else { "s" },
            if with_value.is_empty() { String::new() } else { format!(" ({})", with_value.join(", ")) }
        ));
    }
    if let Some(k) = reference {
        if k >= set.members.len() {
            return Err(format!("the reference is member {}, and the set has {}", k + 1, set.members.len()));
        }
    }
    let scale = file.map_scale.unwrap_or(file.y_scale);
    let grid = common_grid(set, GridRange::Overlap, GridStep::Finest)?;
    let map = map(set, &grid, scale, None)?;
    let difference = match reference {
        Some(k) => Some(crate::plots::map(set, &grid, difference_scale(scale), Some(k))?),
        None => None,
    };
    let regions = file
        .regions
        .iter()
        .map(|r| {
            let on_map = map_region(&map, r, x);
            let on_difference = difference.as_ref().map(|d| map_region(d, r, x));
            let refused = match (&on_map, &on_difference) {
                (Err(e), _) => Some(e.clone()),
                (Ok(_), Some(Err(e))) => Some(format!("on the difference map: {e}")),
                _ => None,
            };
            PageRegion {
                name: r.name.clone(),
                x: r.x,
                x_quantity: r.x_quantity.map(|q| q.word().to_string()),
                variable: r.variable,
                map_range: r.map_range,
                difference_range: r.difference_range,
                window: r.window.clone(),
                blocks: PageBlocks { map: on_map.ok(), difference: on_difference.and_then(|b| b.ok()) },
                refused,
            }
        })
        .collect();
    let members = set
        .members
        .iter()
        .map(|m| PageMember {
            label: m.label.clone(),
            variable: m.variable.as_ref().map(|v| PageValue { value: v.value, source: v.found.words() }),
            pattern: patterns.then(|| PagePattern { x: m.pattern.x.clone(), y: m.pattern.y.clone() }),
        })
        .collect();
    Ok(PageData {
        title: set.title.clone().or_else(|| file.title.clone()),
        path: path.to_string(),
        core: env!("CARGO_PKG_VERSION").to_string(),
        written: crate::formats::files::now_text()[..10].to_string(),
        variable: PageVariable { name: set.variable_name.clone(), unit: set.variable_unit.clone() },
        x: PageX { label: x.quantity.map_or_else(|| set.members.first().map_or("x", |m| m.facts.axis.label()).to_string(), |q| q.label().to_string()), quantity: x.quantity.map(|q| q.word().to_string()) },
        y_word: if file.normalisation == Normalisation::None { "counts" } else { "normalised counts" }.to_string(),
        members,
        maps: PageMaps { map: Some(page_map(&map)), difference: difference.as_ref().map(page_map) },
        clip: PageClip { map: file.map_range, difference: file.difference_range },
        regions,
        regions_refused: file.check_regions(set).err(),
        same_grid: same_grid(set),
        view: "map".into(),
    })
}

/// The page's shell and the viewer's bundle, generated by the app's bundler (`app/ui/tools/bundle-viewer.mjs`, L409,
/// L412) into `core/assets/viewer/`, never edited here.
pub const PAGE_SHELL: &str = include_str!("../assets/viewer/page.html");
pub const VIEWER_BUNDLE: &str = include_str!("../assets/viewer/map.js");
pub const VIEWER_SOURCES: &str = include_str!("../assets/viewer/SOURCES.txt");
/// The shell's two markers, where the JSON and the bundle go.
pub const DATA_MARKER: &str = "<!--DRS_PAGE_DATA-->";
pub const BUNDLE_MARKER: &str = "<!--DRS_VIEWER_BUNDLE-->";

/// The self-contained page (PL38): the shell with the page's JSON and the viewer's bundle in its two script elements,
/// nothing fetched. Inside a script element the text must not hold `</` (it would close the element): the JSON writes
/// every `<` as `\u003c`, the bundle `</` as `<\/` (its `</` are in template literals, where both read back the same).
/// Refused in words when the shell lacks a marker.
#[cfg(feature = "serde")]
pub fn page_html(data: &PageData) -> Result<String, String> {
    let json = serde_json::to_string(data).map_err(|e| e.to_string())?.replace('<', "\\u003c");
    let bundle = VIEWER_BUNDLE.replace("</", "<\\/");
    let (Some(d), Some(b)) = (PAGE_SHELL.find(DATA_MARKER), PAGE_SHELL.find(BUNDLE_MARKER)) else {
        return Err(format!("the page's shell (core/assets/viewer/page.html) lacks the marker {}: regenerate it with the app's bundler", if PAGE_SHELL.contains(DATA_MARKER) { BUNDLE_MARKER } else { DATA_MARKER }));
    };
    // Each marker's text replaced in place, whatever the other's content holds.
    let mut parts = vec![(d, DATA_MARKER.len(), json), (b, BUNDLE_MARKER.len(), bundle)];
    parts.sort_by_key(|p| p.0);
    let mut out = String::with_capacity(PAGE_SHELL.len() + parts.iter().map(|p| p.2.len()).sum::<usize>());
    let mut at = 0;
    for (pos, len, text) in parts {
        out.push_str(&PAGE_SHELL[at..pos]);
        out.push_str(&text);
        at = pos + len;
    }
    out.push_str(&PAGE_SHELL[at..]);
    Ok(out)
}

#[cfg(feature = "serde")]
fn no_displacement() -> Displacement {
    Displacement::None
}

impl Default for SetFile {
    fn default() -> SetFile {
        SetFile {
            title: None,
            variable_name: "T".into(),
            variable_unit: "K".into(),
            x: "2theta".into(),
            y_scale: YScale::Linear,
            order: Order::Variable,
            displacement: Displacement::None,
            displacement_scale: None,
            projection: Displacement::None,
            projection_scale: None,
            map_scale: None,
            y_interval: YInterval::Shared,
            rebin: None,
            normalisation: Normalisation::None,
            reference: None,
            windows: Vec::new(),
            map_range: None,
            difference_range: None,
            regions: Vec::new(),
            members: Vec::new(),
        }
    }
}

use crate::control::syntax::{self as toml, Diagnostic, Entry, Value};

fn bad(span: toml::Span, m: impl Into<String>) -> Vec<Diagnostic> {
    vec![Diagnostic { span, message: m.into(), code: "value", ..Default::default() }]
}

fn number_of(e: &Entry) -> Result<f64, Vec<Diagnostic>> {
    match e.value.value {
        Value::Float(v) => Ok(v),
        Value::Integer(v) => Ok(v as f64),
        _ => Err(bad(e.value.span, format!("{} is a number", e.key))),
    }
}

fn text_of(e: &Entry) -> Result<String, Vec<Diagnostic>> {
    match &e.value.value {
        Value::String(v) => Ok(v.clone()),
        _ => Err(bad(e.value.span, format!("{} is a string", e.key))),
    }
}

fn inline(e: &Entry) -> Result<&[Entry], Vec<Diagnostic>> {
    match &e.value.value {
        Value::Table(t) => Ok(t),
        _ => Err(bad(e.value.span, format!("{} is a table {{ … }}", e.key))),
    }
}

fn known_keys(entries: &[Entry], keys: &[&str], what: &str) -> Result<(), Vec<Diagnostic>> {
    match entries.iter().find(|e| !keys.contains(&e.key.as_str())) {
        Some(e) => Err(bad(e.key_span, format!("\"{}\" is not a key of {what} ({})", e.key, keys.join(", ")))),
        None => Ok(()),
    }
}

/// A normalisation as the set file writes it (PL27): `{ by = "area", from = 20, until = 30, to = 1 }`; a peak's
/// background `"auto"` (with `log_lambda` when given), `"line"`, or `"points"` with `points` and `interpolation`; a
/// reference member counted from 1.
fn normalisation_text(n: &Normalisation) -> String {
    use crate::control::syntax::{float, string};
    match n {
        Normalisation::None => "{ by = \"none\" }".into(),
        Normalisation::Monitor { to } => format!("{{ by = \"monitor\", to = {} }}", float(*to)),
        Normalisation::Time { to } => format!("{{ by = \"time\", to = {} }}", float(*to)),
        Normalisation::Maximum { to } => format!("{{ by = \"maximum\", to = {} }}", float(*to)),
        Normalisation::MaximumIn { from, until, to } => format!("{{ by = \"maximum\", from = {}, until = {}, to = {} }}", float(*from), float(*until), float(*to)),
        Normalisation::Area { from, until, to } => format!("{{ by = \"area\", from = {}, until = {}, to = {} }}", float(*from), float(*until), float(*to)),
        Normalisation::Peak { window, background, to } => {
            let b = match background {
                WindowBackground::Auto { log_lambda: None } => "background = \"auto\"".to_string(),
                WindowBackground::Auto { log_lambda: Some(l) } => format!("background = \"auto\", log_lambda = {}", float(*l)),
                WindowBackground::Line => "background = \"line\"".to_string(),
                WindowBackground::Points { points, interpolation } => format!(
                    "background = \"points\", interpolation = {}, points = [{}]",
                    string(interpolation.name()),
                    points.iter().map(|(x, y)| format!("[{}, {}]", float(*x), float(*y))).collect::<Vec<_>>().join(", ")
                ),
            };
            format!("{{ by = \"peak\", name = {}, from = {}, until = {}, to = {}, {b} }}", string(&window.name), float(window.from), float(window.to), float(*to))
        }
        Normalisation::Reference { member, from, until } => format!("{{ by = \"reference\", member = {}, from = {}, until = {} }}", member + 1, float(*from), float(*until)),
    }
}

/// A colour range as the set file writes it (PL37): `{ from = 120.0, to = 900.0 }`.
fn range_text(r: &Range) -> String {
    use crate::control::syntax::float;
    format!("{{ from = {}, to = {} }}", float(r.from), float(r.to))
}

/// A colour range read, refused when `from` is not below `to`.
fn range_of(e: &Entry, what: &str) -> Result<Range, Vec<Diagnostic>> {
    let t = inline(e)?;
    known_keys(t, &["from", "to"], what)?;
    let need = |k: &str| t.iter().find(|x| x.key == k).map(number_of).transpose()?.ok_or_else(|| bad(e.value.span, format!("{what} needs {k}")));
    let (from, to) = (need("from")?, need("to")?);
    if from >= to {
        return Err(bad(e.value.span, format!("{what}: from = {from} is not below to = {to}; a range is written from its lower end")));
    }
    Ok(Range { from, to })
}

/// A pair `[a, b]` read, refused when its ends are reversed (PL37's region).
fn pair_of(e: &Entry, what: &str) -> Result<[f64; 2], Vec<Diagnostic>> {
    let Value::Array(items) = &e.value.value else { return Err(bad(e.value.span, format!("{what} is a pair [from, to]"))) };
    let v: Vec<f64> = items.iter().filter_map(|i| match i.value {
        Value::Float(v) => Some(v),
        Value::Integer(v) => Some(v as f64),
        _ => None,
    }).collect();
    if v.len() != 2 || items.len() != 2 {
        return Err(bad(e.value.span, format!("{what} is a pair of numbers [from, to]")));
    }
    if v[0] >= v[1] {
        return Err(bad(e.value.span, format!("{what} = [{}, {}] runs backwards or is empty: the smaller end first", v[0], v[1])));
    }
    Ok([v[0], v[1]])
}

/// The set file's `normalisation = { … }` read (`normalisation_text`).
fn normalisation_of(e: &Entry) -> Result<Normalisation, Vec<Diagnostic>> {
    let t = inline(e)?;
    known_keys(t, &["by", "to", "from", "until", "member", "name", "background", "log_lambda", "interpolation", "points"], "normalisation")?;
    let num = |k: &str| t.iter().find(|x| x.key == k).map(number_of).transpose();
    let need = |k: &str| num(k)?.ok_or_else(|| bad(e.value.span, format!("this normalisation needs {k}")));
    let to = num("to")?.unwrap_or(1.0);
    let by = t.iter().find(|x| x.key == "by").map(text_of).transpose()?.ok_or_else(|| bad(e.value.span, "a normalisation says by = \"…\""))?;
    Ok(match by.as_str() {
        "none" => Normalisation::None,
        "monitor" => Normalisation::Monitor { to },
        "time" => Normalisation::Time { to },
        "maximum" => match (num("from")?, num("until")?) {
            (Some(from), Some(until)) => Normalisation::MaximumIn { from, until, to },
            _ => Normalisation::Maximum { to },
        },
        "area" => Normalisation::Area { from: need("from")?, until: need("until")?, to },
        "peak" => {
            let name = t.iter().find(|x| x.key == "name").map(text_of).transpose()?.unwrap_or_else(|| "peak".into());
            let background = match t.iter().find(|x| x.key == "background").map(text_of).transpose()?.as_deref() {
                None | Some("auto") => WindowBackground::Auto { log_lambda: num("log_lambda")? },
                Some("line") => WindowBackground::Line,
                Some("points") => {
                    let interpolation = match t.iter().find(|x| x.key == "interpolation").map(text_of).transpose()? {
                        Some(w) => crate::background::Interpolation::from_name(&w).ok_or_else(|| bad(e.value.span, "interpolation is \"linear\" or \"spline\""))?,
                        None => crate::background::Interpolation::Linear,
                    };
                    let p = t.iter().find(|x| x.key == "points").ok_or_else(|| bad(e.value.span, "a background of points needs points = [[x, y], …]"))?;
                    let Value::Array(rows) = &p.value.value else { return Err(bad(p.value.span, "points = [[x, y], …]")) };
                    let mut points = Vec::new();
                    for r in rows {
                        let Value::Array(xy) = &r.value else { return Err(bad(r.span, "a point is [x, y]")) };
                        let v: Vec<f64> = xy.iter().filter_map(|i| match i.value {
                            Value::Float(v) => Some(v),
                            Value::Integer(v) => Some(v as f64),
                            _ => None,
                        }).collect();
                        if v.len() != 2 {
                            return Err(bad(r.span, "a point is [x, y]"));
                        }
                        points.push((v[0], v[1]));
                    }
                    WindowBackground::Points { points, interpolation }
                }
                Some(o) => return Err(bad(e.value.span, format!("background = \"{o}\": auto, line or points"))),
            };
            Normalisation::Peak { window: Window { name, from: need("from")?, to: need("until")? }, background, to }
        }
        "reference" => {
            let k = need("member")?;
            if k < 1.0 || k.fract() != 0.0 {
                return Err(bad(e.value.span, "member counts from 1"));
            }
            Normalisation::Reference { member: k as usize - 1, from: need("from")?, until: need("until")? }
        }
        o => return Err(bad(e.value.span, format!("by = \"{o}\": none, monitor, time, maximum, area, peak or reference"))),
    })
}

/// How `SetFile::open_by` gets the set of the files a set file names (by path, or by name and content).
type Opener<'a> = dyn Fn(&[(String, Option<String>)]) -> Result<PlotSet, String> + 'a;

impl SetFile {
    /// Read a set file; what is wrong, with its place.
    pub fn read(text: &str) -> Result<SetFile, Vec<Diagnostic>> {
        let doc = toml::parse(text);
        if doc.has_errors() {
            return Err(doc.diagnostics);
        }
        let root = doc.tables.first().ok_or_else(|| bad(Default::default(), "an empty file"))?;
        known_keys(&root.entries, &["format", "title", "variable", "x", "y_scale", "order", "displacement", "projection", "map_scale", "y_interval", "rebin", "normalisation", "reference", "windows", "map_range", "difference_range"], "the top of a plot set")?;
        match root.get("format") {
            Some(e) if text_of(e)? == SET_FORMAT => {}
            Some(e) => return Err(bad(e.value.span, format!("format is \"{SET_FORMAT}\""))),
            None => return Err(bad(Default::default(), format!("a plot set starts with format = \"{SET_FORMAT}\""))),
        }
        let mut f = SetFile::default();
        if let Some(e) = root.get("title") {
            f.title = Some(text_of(e)?);
        }
        if let Some(e) = root.get("variable") {
            let t = inline(e)?;
            known_keys(t, &["name", "unit"], "variable")?;
            for x in t {
                match x.key.as_str() {
                    "name" => f.variable_name = text_of(x)?,
                    _ => f.variable_unit = text_of(x)?,
                }
            }
        }
        if let Some(e) = root.get("x") {
            f.x = text_of(e)?;
            if !["2theta", "d", "q", "stl", "s", "tof", "energy"].contains(&f.x.as_str()) {
                return Err(bad(e.value.span, "x is \"2theta\", \"d\", \"q\", \"stl\", \"s\", \"tof\" or \"energy\""));
            }
        }
        if let Some(e) = root.get("y_scale") {
            f.y_scale = YScale::from_word(&text_of(e)?).ok_or_else(|| bad(e.value.span, "y_scale is \"linear\", \"sqrt\" or \"log\""))?;
        }
        if let Some(e) = root.get("order") {
            f.order = match text_of(e)?.as_str() {
                "variable" => Order::Variable,
                "name" => Order::Name,
                "given" => Order::Given,
                _ => return Err(bad(e.value.span, "order is \"variable\", \"name\" or \"given\"")),
            };
        }
        // The stack's displacement and the projection's (PL27), each with the y scale its steps were set on.
        let steps = |e: &Entry, what: &str| -> Result<(Displacement, Option<YScale>), Vec<Diagnostic>> {
            let t = inline(e)?;
            known_keys(t, &["rule", "dx", "dy", "scale"], what)?;
            let get = |k: &str| t.iter().find(|x| x.key == k).map(number_of).transpose().map(|v| v.unwrap_or(0.0));
            let (dx, dy) = (get("dx")?, get("dy")?);
            let rule = t.iter().find(|x| x.key == "rule").map(text_of).transpose()?.unwrap_or_else(|| "none".into());
            let d = match rule.as_str() {
                "none" => Displacement::None,
                "step" => Displacement::Step { dx, dy },
                "variable" => Displacement::Variable { dx, dy },
                _ => return Err(bad(e.value.span, format!("the {what}'s rule is \"none\", \"step\" or \"variable\""))),
            };
            let scale = match t.iter().find(|x| x.key == "scale") {
                Some(x) => Some(YScale::from_word(&text_of(x)?).ok_or_else(|| bad(x.value.span, "scale is \"linear\", \"sqrt\" or \"log\""))?),
                None => None,
            };
            Ok((d, scale))
        };
        if let Some(e) = root.get("displacement") {
            (f.displacement, f.displacement_scale) = steps(e, "displacement")?;
        }
        if let Some(e) = root.get("projection") {
            (f.projection, f.projection_scale) = steps(e, "projection")?;
        }
        if let Some(e) = root.get("map_scale") {
            f.map_scale = Some(YScale::from_word(&text_of(e)?).ok_or_else(|| bad(e.value.span, "map_scale is \"linear\", \"sqrt\" or \"log\""))?);
        }
        if let Some(e) = root.get("y_interval") {
            f.y_interval = match text_of(e)?.as_str() {
                "shared" => YInterval::Shared,
                "own" => YInterval::Own,
                _ => return Err(bad(e.value.span, "y_interval is \"shared\" (the whole set's) or \"own\" (each member's)")),
            };
        }
        if let Some(e) = root.get("rebin") {
            let v = number_of(e)?;
            if v.is_nan() || v <= 0.0 {
                return Err(bad(e.value.span, "rebin is the new step, a positive number"));
            }
            f.rebin = Some(v);
        }
        if let Some(e) = root.get("normalisation") {
            f.normalisation = normalisation_of(e)?;
        }
        if let Some(e) = root.get("reference") {
            f.reference = Some(text_of(e)?);
        }
        if let Some(e) = root.get("windows") {
            let Value::Array(items) = &e.value.value else { return Err(bad(e.value.span, "windows is a list of { name, from, to }")) };
            for it in items {
                let Value::Table(t) = &it.value else { return Err(bad(it.span, "a window is { name = …, from = …, to = … }")) };
                known_keys(t, &["name", "from", "to"], "a window")?;
                let find = |k: &str| t.iter().find(|x| x.key == k).ok_or_else(|| bad(it.span, format!("a window needs {k}")));
                f.windows.push(Window { name: text_of(find("name")?)?, from: number_of(find("from")?)?, to: number_of(find("to")?)? });
            }
        }
        if let Some(e) = root.get("map_range") {
            f.map_range = Some(range_of(e, "map_range")?);
        }
        if let Some(e) = root.get("difference_range") {
            f.difference_range = Some(range_of(e, "difference_range")?);
        }
        for t in doc.tables.iter().skip(1) {
            if t.array && t.path == ["regions"] {
                // A region of the map saved by name (PL37).
                let at = t.header.unwrap_or_default();
                known_keys(&t.entries, &["name", "x", "x_quantity", "variable", "map_range", "difference_range", "window"], "[[regions]]")?;
                let name = text_of(t.get("name").ok_or_else(|| bad(at, "a region needs name"))?)?;
                if f.regions.iter().any(|r| r.name == name) {
                    return Err(bad(at, format!("two regions are named \"{name}\": a region is known by its name")));
                }
                let x = pair_of(t.get("x").ok_or_else(|| bad(at, format!("region \"{name}\" needs x = [from, to]")))?, &format!("region \"{name}\": x"))?;
                let variable = pair_of(t.get("variable").ok_or_else(|| bad(at, format!("region \"{name}\" needs variable = [from, to]")))?, &format!("region \"{name}\": variable"))?;
                let x_quantity = match t.get("x_quantity") {
                    Some(e) => Some(XQuantity::from_word(&text_of(e)?).ok_or_else(|| bad(e.value.span, "x_quantity is \"2theta\", \"d\", \"q\", \"stl\" or \"s\""))?),
                    None => None,
                };
                let map_range = t.get("map_range").map(|e| range_of(e, &format!("region \"{name}\": map_range"))).transpose()?;
                let difference_range = t.get("difference_range").map(|e| range_of(e, &format!("region \"{name}\": difference_range"))).transpose()?;
                let window = t.get("window").map(text_of).transpose()?;
                f.regions.push(Region { name, x, x_quantity, variable, map_range, difference_range, window });
                continue;
            }
            if !(t.array && t.path == ["member"]) {
                return Err(bad(t.header.unwrap_or_default(), format!("{} is not a table of a plot set ([[member]], [[regions]])", t.name())));
            }
            known_keys(&t.entries, &["file", "scan", "variable", "variables", "visible", "opacity", "colour"], "[[member]]")?;
            let mut m = SetMember { file: text_of(t.get("file").ok_or_else(|| bad(t.header.unwrap_or_default(), "a member needs file"))?)?, ..Default::default() };
            if let Some(e) = t.get("scan") {
                match e.value.value {
                    Value::Integer(k) if k >= 1 => m.scan = Some(k as usize),
                    _ => return Err(bad(e.value.span, "scan counts from 1")),
                }
            }
            if let Some(e) = t.get("variable") {
                let v = inline(e)?;
                known_keys(v, &["value", "source", "temperatures"], "a member's variable")?;
                let value = number_of(v.iter().find(|x| x.key == "value").ok_or_else(|| bad(e.value.span, "a member's variable needs value"))?)?;
                let source = v.iter().find(|x| x.key == "source").map(text_of).transpose()?;
                let file = v.iter().find(|x| x.key == "temperatures").map(text_of).transpose()?;
                let found = match (source, file) {
                    (Some(s), None) => Found::Typed(s),
                    (None, Some(f)) => Found::TemperatureFile(f),
                    (Some(_), Some(_)) => return Err(bad(e.value.span, "a value has one source: source = \"…\" (typed) or temperatures = \"…\" (a temperatures file), not both")),
                    (None, None) => return Err(bad(e.value.span, "a value says its source (source = \"…\", or temperatures = \"…\"): never a value without one (L169)")),
                };
                m.value = Some((value, found));
            }
            if let Some(e) = t.get("visible") {
                match e.value.value {
                    Value::Bool(b) => m.visible = Some(b),
                    _ => return Err(bad(e.value.span, "visible is true or false")),
                }
            }
            if let Some(e) = t.get("opacity") {
                let o = number_of(e)?;
                if !(0.0..=1.0).contains(&o) {
                    return Err(bad(e.value.span, "opacity runs from 0 to 1"));
                }
                m.opacity = Some(o);
            }
            if let Some(e) = t.get("colour") {
                m.colour = Some(text_of(e)?);
            }
            // Other variables than the set's, each { value, source } (W3); never a value without its source (L169).
            if let Some(e) = t.get("variables") {
                for v in inline(e)? {
                    let w = inline(v)?;
                    known_keys(w, &["value", "source"], &format!("the member's variable {}", v.key))?;
                    let value = number_of(w.iter().find(|x| x.key == "value").ok_or_else(|| bad(v.value.span, format!("the variable {} needs value", v.key)))?)?;
                    let source = text_of(w.iter().find(|x| x.key == "source").ok_or_else(|| bad(v.value.span, format!("the variable {} says its source (source = \"…\"): never a value without one (L169)", v.key)))?)?;
                    if v.key == f.variable_name {
                        return Err(bad(v.key_span, format!("{} is the set's variable: its value is the member's variable = {{ … }}, not one of its other variables", v.key)));
                    }
                    if m.variables.iter().any(|x| x.0 == v.key) {
                        return Err(bad(v.key_span, format!("the variable {} is given twice", v.key)));
                    }
                    m.variables.push((v.key.clone(), value, source));
                }
            }
            f.members.push(m);
        }
        if f.members.is_empty() {
            return Err(bad(Default::default(), "a plot set needs a [[member]]"));
        }
        Ok(f)
    }

    /// The file's text: what differs from the defaults, in the order `read` takes it.
    pub fn to_text(&self) -> String {
        use crate::control::syntax::{float, string};
        let mut t = format!("format = {}\n", string(SET_FORMAT));
        if let Some(x) = &self.title {
            t.push_str(&format!("title = {}\n", string(x)));
        }
        t.push_str(&format!("variable = {{ name = {}, unit = {} }}\n", string(&self.variable_name), string(&self.variable_unit)));
        t.push_str(&format!("x = {}\ny_scale = {}\n", string(&self.x), string(self.y_scale.word())));
        t.push_str(&format!("order = {}\n", string(match self.order {
            Order::Variable => "variable",
            Order::Name => "name",
            Order::Given => "given",
        })));
        let steps = |key: &str, d: Displacement, scale: Option<YScale>| -> String {
            let sc = scale.map_or(String::new(), |s| format!(", scale = {}", string(s.word())));
            match d {
                Displacement::None if scale.is_none() => String::new(),
                Displacement::None => format!("{key} = {{ rule = \"none\"{sc} }}\n"),
                Displacement::Step { dx, dy } => format!("{key} = {{ rule = \"step\", dx = {}, dy = {}{sc} }}\n", float(dx), float(dy)),
                Displacement::Variable { dx, dy } => format!("{key} = {{ rule = \"variable\", dx = {}, dy = {}{sc} }}\n", float(dx), float(dy)),
            }
        };
        t.push_str(&steps("displacement", self.displacement, self.displacement_scale));
        t.push_str(&steps("projection", self.projection, self.projection_scale));
        if let Some(m) = self.map_scale {
            t.push_str(&format!("map_scale = {}\n", string(m.word())));
        }
        if self.y_interval != YInterval::Shared {
            t.push_str(&format!("y_interval = {}\n", string(self.y_interval.word())));
        }
        if let Some(r) = self.rebin {
            t.push_str(&format!("rebin = {}\n", float(r)));
        }
        if self.normalisation != Normalisation::None {
            t.push_str(&format!("normalisation = {}\n", normalisation_text(&self.normalisation)));
        }
        if let Some(r) = &self.reference {
            t.push_str(&format!("reference = {}\n", string(r)));
        }
        if !self.windows.is_empty() {
            let w: Vec<String> = self.windows.iter().map(|w| format!("{{ name = {}, from = {}, to = {} }}", string(&w.name), float(w.from), float(w.to))).collect();
            t.push_str(&format!("windows = [{}]\n", w.join(", ")));
        }
        if let Some(r) = &self.map_range {
            t.push_str(&format!("map_range = {}\n", range_text(r)));
        }
        if let Some(r) = &self.difference_range {
            t.push_str(&format!("difference_range = {}\n", range_text(r)));
        }
        for r in &self.regions {
            t.push_str(&format!("\n[[regions]]\nname = {}\nx = [{}, {}]\n", string(&r.name), float(r.x[0]), float(r.x[1])));
            if let Some(q) = r.x_quantity {
                t.push_str(&format!("x_quantity = {}\n", string(q.word())));
            }
            t.push_str(&format!("variable = [{}, {}]\n", float(r.variable[0]), float(r.variable[1])));
            if let Some(m) = &r.map_range {
                t.push_str(&format!("map_range = {}\n", range_text(m)));
            }
            if let Some(m) = &r.difference_range {
                t.push_str(&format!("difference_range = {}\n", range_text(m)));
            }
            if let Some(w) = &r.window {
                t.push_str(&format!("window = {}\n", string(w)));
            }
        }
        for m in &self.members {
            t.push_str(&format!("\n[[member]]\nfile = {}\n", string(&m.file)));
            if let Some(k) = m.scan {
                t.push_str(&format!("scan = {k}\n"));
            }
            match &m.value {
                Some((v, Found::TemperatureFile(f))) => t.push_str(&format!("variable = {{ value = {}, temperatures = {} }}\n", float(*v), string(f))),
                Some((v, Found::Typed(src))) => t.push_str(&format!("variable = {{ value = {}, source = {} }}\n", float(*v), string(src))),
                _ => {}
            }
            if let Some(b) = m.visible {
                t.push_str(&format!("visible = {b}\n"));
            }
            if let Some(o) = m.opacity {
                t.push_str(&format!("opacity = {}\n", float(o)));
            }
            if let Some(c) = &m.colour {
                t.push_str(&format!("colour = {}\n", string(c)));
            }
            if !m.variables.is_empty() {
                let v: Vec<String> = m.variables.iter().map(|(k, v, src)| format!("{k} = {{ value = {}, source = {} }}", float(*v), string(src))).collect();
                t.push_str(&format!("variables = {{ {} }}\n", v.join(", ")));
            }
        }
        t
    }

    /// The set it describes, its files read relative to `base` (the project's folder), each after Data's operations
    /// (PL20), each value typed or read from a temperatures file set with its source, in the file's order: by the variable, by name, or as given.
    pub fn open(&self, base: &Path, reader: Option<&str>, operations: &Operations) -> Result<PlotSet, String> {
        self.open_by(&|files| open_with(base, files, reader, self.title.clone(), operations))
    }

    /// The set it describes from files given by name and content (W1): each member's file taken from `files` by the
    /// name the set file gives it (relative to the project's folder, as written); refused, naming it, when a member's
    /// file is not given. As `open` otherwise.
    pub fn open_from(&self, files: &[OpenedFile], reader: Option<&str>, operations: &Operations) -> Result<PlotSet, String> {
        self.open_by(&|wanted| {
            let given: Vec<OpenedFile> = wanted
                .iter()
                .map(|(name, _)| files.iter().find(|f| &f.name == name).cloned().ok_or_else(|| format!("{name}: the set file names it, and no file of that name was given")))
                .collect::<Result<_, _>>()?;
            // the files the members' operations name are looked up among all those given
            let mut all = given;
            for f in files {
                if !all.iter().any(|g| g.name == f.name) {
                    all.push(f.clone());
                }
            }
            let mut set = open_from_with(&all, reader, self.title.clone(), operations)?;
            set.members.retain(|m| wanted.iter().any(|w| w.0 == m.file));
            Ok(set)
        })
    }

    fn open_by(&self, opener: &Opener) -> Result<PlotSet, String> {
        // `x` is the view's quantity (PL10, PL27): the members open as their files give them; `x_view` shows them so.
        let mut files: Vec<(String, Option<String>)> = Vec::new();
        for m in &self.members {
            if !files.iter().any(|f| f.0 == m.file) {
                files.push((m.file.clone(), None));
            }
        }
        let all = opener(&files)?;
        let mut members = Vec::new();
        for m in &self.members {
            let mine: Vec<&Member> = all.members.iter().filter(|x| x.file == m.file && m.scan.is_none_or(|k| x.scan + 1 == k)).collect();
            if mine.is_empty() {
                return Err(format!("{}: no scan {}", m.file, m.scan.unwrap_or(1)));
            }
            for x in mine {
                let mut x = x.clone();
                if let Some((v, found)) = &m.value {
                    let found = Variable { value: *v, found: found.clone() };
                    // A temperatures file's value is one of the member's candidates too, as when it was read (PL21).
                    if matches!(found.found, Found::TemperatureFile(_)) && !x.candidates.contains(&found) {
                        x.candidates.push(found.clone());
                    }
                    x.variable = Some(found);
                }
                // The other variables (W3): those the file types, each with its source; a temperature the member's own
                // sources give (PL2: the header, the buffer's comment, the title, the name) when the set's variable is
                // another one and the file types none.
                x.others = m.variables.iter().map(|(name, value, src)| Other { name: name.clone(), value: *value, found: Found::Typed(src.clone()), unit: String::new() }).collect();
                if self.variable_name != "T" && !x.others.iter().any(|o| o.name == "T") {
                    if let Some(t) = x.candidates.iter().find(|c| !matches!(c.found, Found::Typed(_) | Found::TemperatureFile(_))) {
                        x.others.push(Other { name: "T".into(), value: t.value, found: t.found.clone(), unit: all.variable_unit.clone() });
                    }
                }
                members.push(x);
            }
        }
        let mut set = PlotSet { title: self.title.clone(), variable_name: self.variable_name.clone(), variable_unit: self.variable_unit.clone(), members };
        if let Some(step) = self.rebin {
            set.rebin(step)?;
        }
        // a temperature in the file's unit (K or °C), whatever unit its source gives (W8)
        if self.variable_name == "T" && (self.variable_unit == "K" || self.variable_unit == "°C") {
            let unit = self.variable_unit.clone();
            set.temperatures_in(&unit);
        }
        // W3's t: the time since the set's first member, from the seconds a name gives or the modification time, when
        // the set's variable is another one and the file types none
        if self.variable_name != "t" {
            set.fill_elapsed(false);
        }
        match self.order {
            Order::Variable => set.sort_by_variable(),
            Order::Name => set.members.sort_by(|a, b| a.label.cmp(&b.label)),
            Order::Given => {}
        }
        Ok(set)
    }

    /// The reference member's index in an opened set.
    pub fn reference_in(&self, set: &PlotSet) -> Option<usize> {
        let r = self.reference.as_ref()?;
        set.members.iter().position(|m| &m.file == r)
    }

    /// The set's x quantity as the file says it (`x`); none for TOF and energy, which are not converted.
    pub fn x_quantity(&self) -> Option<XQuantity> {
        XQuantity::from_word(&self.x)
    }

    /// The regions checked against the set opened from the file (PL37), the set shown in the file's `x` with each
    /// member's λ from its file: a region in another quantity whose conversion needs a λ (one of the two is 2θ) is
    /// refused, with the members that state none named; a region in a quantity on a set in TOF or energy likewise.
    pub fn check_regions(&self, set: &PlotSet) -> Result<(), String> {
        let to = self.x_quantity();
        for r in &self.regions {
            let Some(from) = r.x_quantity else { continue };
            let Some(to) = to else { return Err(format!("region \"{}\" is in {}, and the set's x, {}, is not converted", r.name, from.word(), self.x)) };
            if needs_lambda(from, to) {
                let lacking: Vec<&str> = set.members.iter().filter(|m| m.facts.wavelength.is_none()).map(|m| m.label.as_str()).collect();
                if !lacking.is_empty() {
                    return Err(format!("region \"{}\" is in {} and the set in {}: converting it needs each member's wavelength, and these state none: {}", r.name, from.word(), to.word(), lacking.join(", ")));
                }
            }
        }
        Ok(())
    }

    /// What the file says beyond its members, in words (`drs plots describe`; PL12, PL27, PL37).
    pub fn words(&self) -> Vec<String> {
        let mut w = vec![format!("{} members; variable {} ({}); x {}; y on the {} scale", self.members.len(), self.variable_name, self.variable_unit, self.x, self.y_scale.word())];
        if let Some(t) = &self.title {
            w.insert(0, format!("title \"{t}\""));
        }
        let map_scale = self.map_scale.unwrap_or(self.y_scale).word();
        w.push(format!("the map on the {map_scale} scale{}", if self.map_scale.is_some() { "" } else { " (the curves')" }));
        w.push(match self.y_interval {
            YInterval::Shared => "the curves' y interval: the whole set's".into(),
            YInterval::Own => "the curves' y interval: each member's own".into(),
        });
        if let Some(r) = self.rebin {
            w.push(format!("the members rebinned on opening to a step of {r} in their x, the counts averaged, σ propagated"));
        }
        w.push(match &self.map_range {
            Some(r) => format!("the map's colour range clipped to {} to {} on the {map_scale} scale, values beyond saturated", r.from, r.to),
            None => "the map's colour range: the core's (no clip)".into(),
        });
        w.push(match &self.difference_range {
            Some(r) => format!("the difference map's colour range clipped to {} to {}, values beyond saturated", r.from, r.to),
            None => "the difference map's colour range: the core's (no clip)".into(),
        });
        if let Some(r) = &self.reference {
            w.push(format!("reference {r}"));
        }
        if !self.windows.is_empty() {
            w.push(format!("windows {}", self.windows.iter().map(|x| format!("{} {}–{}", x.name, x.from, x.to)).collect::<Vec<_>>().join(", ")));
        }
        if self.regions.is_empty() {
            w.push("no region saved".into());
        }
        let (map_scale, difference_scale) = (self.map_scale.unwrap_or(self.y_scale), difference_scale(self.map_scale.unwrap_or(self.y_scale)));
        for r in &self.regions {
            let q = r.x_quantity.map_or(format!("{} (the set's x)", self.x), |q| q.word().to_string());
            // Each range with the map it is on and that map's scale (P65).
            let ranges = [(r.map_range, "the map", map_scale), (r.difference_range, "the difference map", difference_scale)]
                .iter()
                .filter_map(|(m, on, sc)| m.map(|m| format!(", its colour range on {on} {} to {} ({} scale)", m.from, m.to, sc.word())))
                .collect::<String>();
            w.push(format!(
                "region \"{}\": {} {} to {}, {} {} to {} {}{}{}",
                r.name,
                q,
                r.x[0],
                r.x[1],
                self.variable_name,
                r.variable[0],
                r.variable[1],
                self.variable_unit,
                ranges,
                r.window.as_ref().map_or(String::new(), |x| format!(", drawn from the window {x}"))
            ));
        }
        w
    }

    /// A set file for a set (members as they are, each value typed or read from a temperatures file kept with its
    /// source, PL23), the rest default.
    pub fn of(set: &PlotSet) -> SetFile {
        let members = set
            .members
            .iter()
            .map(|m| SetMember {
                file: m.file.clone(),
                scan: (m.scans > 1).then_some(m.scan + 1),
                value: match &m.variable {
                    Some(Variable { value, found: f @ (Found::Typed(_) | Found::TemperatureFile(_)) }) => Some((*value, f.clone())),
                    _ => None,
                },
                // the other variables typed (W3), those a file read again does not give
                variables: m.others.iter().filter_map(|o| match &o.found {
                    Found::Typed(src) => Some((o.name.clone(), o.value, src.clone())),
                    _ => None,
                }).collect(),
                ..Default::default()
            })
            .collect();
        SetFile { title: set.title.clone(), variable_name: set.variable_name.clone(), variable_unit: set.variable_unit.clone(), members, ..Default::default() }
    }
}

/// A file's name relative to a folder, '/' between its parts, ".." where it lies beside it; the absolute name when the
/// two share no root (another disk).
pub fn relative_to(file: &Path, folder: &Path) -> String {
    let abs = |p: &Path| -> PathBuf {
        let p = if p.is_absolute() { p.to_path_buf() } else { std::env::current_dir().unwrap_or_default().join(p) };
        let mut out = PathBuf::new();
        for c in p.components() {
            match c {
                std::path::Component::ParentDir => {
                    out.pop();
                }
                std::path::Component::CurDir => {}
                c => out.push(c),
            }
        }
        out
    };
    let (f, d) = (abs(file), abs(folder));
    let (fc, dc): (Vec<_>, Vec<_>) = (f.components().collect(), d.components().collect());
    let common = fc.iter().zip(&dc).take_while(|(a, b)| a == b).count();
    if common == 0 {
        return f.to_string_lossy().into_owned();
    }
    let mut parts: Vec<String> = vec!["..".to_string(); dc.len() - common];
    parts.extend(fc[common..].iter().map(|c| c.as_os_str().to_string_lossy().into_owned()));
    parts.join("/")
}

/// The project's folder for a set file: the folder above `Plots/` when the file is in one so named (the file
/// dictionary's place), else the file's own folder.
pub fn set_base(path: &Path) -> PathBuf {
    let dir = path.parent().map(Path::to_path_buf).unwrap_or_default();
    if dir.file_name().is_some_and(|n| n == "Plots") {
        dir.parent().map(Path::to_path_buf).unwrap_or_default()
    } else {
        dir
    }
}

/// The name of a set's x column in a CSV (PL7).
fn x_column(a: Axis) -> &'static str {
    match a {
        Axis::TwoTheta => "2theta_deg",
        Axis::D => "d_angstrom",
        Axis::Q => "q_inverse_angstrom",
        Axis::TimeOfFlight => "tof_microseconds",
        Axis::Energy => "energy_kev",
    }
}

fn member_line(set: &PlotSet, k: usize) -> String {
    let m = &set.members[k];
    let v = m.variable.as_ref().map_or(format!("{}: nothing says it", set.variable_name), |v| format!("{} = {} {} ({})", set.variable_name, v.value, set.variable_unit, v.found.words()));
    format!("{} (scan {} of {}, {}), {v}", m.file, m.scan + 1, m.scans, m.facts.format)
}

/// The set with every member's x in another quantity (PL7's CSV in d or Q): each member converted as Data's convert
/// converts it, with its own λ or `wavelength` when the file states none; refused, naming the member, when neither
/// gives one. The facts' range and step follow; the rest as read.
pub fn converted(set: &PlotSet, to: Axis, wavelength: Option<f64>) -> Result<PlotSet, String> {
    let mut out = set.clone();
    for m in out.members.iter_mut() {
        if m.facts.axis == to {
            continue;
        }
        let lambda = m.facts.wavelength.or(wavelength).ok_or(format!("{}: the file states no wavelength; give one", m.label))?;
        let mut s = Scan::new(m.label.clone(), m.facts.axis, m.pattern.clone());
        crate::data_ops::Operation::Convert { to, wavelength: Some(lambda) }.apply(&mut s, &|_, _| Err("no file".into())).map_err(|e| format!("{}: {e}", m.label))?;
        m.pattern = s.pattern;
        m.facts.axis = to;
        m.facts.range = crate::data::extent(&m.pattern.x).unwrap_or((f64::NAN, f64::NAN));
        m.facts.step = m.pattern.mean_step();
        m.facts.points = m.pattern.len();
        m.facts.wavelength = Some(lambda);
    }
    Ok(out)
}

/// A set as CSV, wide (PL7, PL12's export): x on the grid, then y and σ of each member, what each column is in `#`
/// lines above the header (its file, scan, format and value of the variable with its source). An empty cell where a
/// member has no point; a σ of no weight is written inf. The x quantity is the members' (`converted` for another).
pub fn csv(set: &PlotSet, grid: &[f64]) -> String {
    csv_with(set, grid, &[])
}

/// `csv` with lines saying how the members were shown (PL25, `PlotSet::as_shown`), as `#` lines after the first.
pub fn csv_with(set: &PlotSet, grid: &[f64], notes: &[String]) -> String {
    let res = on_grid(set, grid);
    let axis = set.members.first().map_or(Axis::TwoTheta, |m| m.facts.axis);
    let mut t = format!("# {}: {} members on a grid of {} points, x = {} (written by drs {})\n", set.title.clone().unwrap_or_else(|| "a plot set".into()), set.members.len(), grid.len(), axis.label(), env!("CARGO_PKG_VERSION"));
    for l in notes {
        t.push_str(&format!("# {l}\n"));
    }
    for (k, r) in res.iter().enumerate() {
        t.push_str(&format!("# y{}: {}, {:?} on the grid\n", k + 1, member_line(set, k), r.how));
    }
    t.push_str(x_column(axis));
    for k in 1..=set.members.len() {
        t.push_str(&format!(",y{k},sigma{k}"));
    }
    t.push('\n');
    let cell = |v: Option<f64>| v.map_or(String::new(), |v| if v.is_infinite() { "inf".into() } else { format!("{v}") });
    for (j, x) in grid.iter().enumerate() {
        t.push_str(&format!("{x}"));
        for r in &res {
            t.push_str(&format!(",{},{}", cell(r.y[j]), cell(r.sigma[j])));
        }
        t.push('\n');
    }
    t
}

/// A set as CSV, long (PL7): one row per member and point, on each member's own points: member (from 1), x, y, σ; what
/// each member is in `#` lines.
pub fn csv_long(set: &PlotSet) -> String {
    csv_long_with(set, &[])
}

/// `csv_long` with lines saying how the members were shown (PL25), as `#` lines after the first.
pub fn csv_long_with(set: &PlotSet, notes: &[String]) -> String {
    let axis = set.members.first().map_or(Axis::TwoTheta, |m| m.facts.axis);
    let mut t = format!("# {}: {} members, each on its own points, x = {} (written by drs {})\n", set.title.clone().unwrap_or_else(|| "a plot set".into()), set.members.len(), axis.label(), env!("CARGO_PKG_VERSION"));
    for l in notes {
        t.push_str(&format!("# {l}\n"));
    }
    for k in 0..set.members.len() {
        t.push_str(&format!("# member {}: {}\n", k + 1, member_line(set, k)));
    }
    t.push_str(&format!("member,{},y,sigma\n", x_column(axis)));
    for (k, m) in set.members.iter().enumerate() {
        let p = &m.pattern;
        for i in 0..p.len() {
            let s = if p.sigma[i].is_infinite() { "inf".to_string() } else { format!("{}", p.sigma[i]) };
            t.push_str(&format!("{},{},{},{s}\n", k + 1, p.x[i], p.y[i]));
        }
    }
    t
}

// PL6: normalisation as a view -----------------------------------------------------------------------------------------

/// How the members of a set are normalised for display (PL6): a view, never an operation written into a file's list
/// (Data's `Operation::Normalize` is that). Each member's y and σ are multiplied by its factor; the factor is taken as
/// exact (its own uncertainty, from the area or the integral it divides by, left out).
#[derive(Clone, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Normalisation {
    #[default]
    None,
    /// y · to / monitor (the reader's monitor count).
    Monitor { to: f64 },
    /// y · to / counting time (the reader's).
    Time { to: f64 },
    /// y · to / max y.
    Maximum { to: f64 },
    /// y · to / max y over [from, until].
    MaximumIn { from: f64, until: f64, to: f64 },
    /// The same area over [from, until] for every member (Josu, L169): y · to / ∫ y dx there, the integral by
    /// `window_value`'s rule (each point over the part of its interval inside the range) above no background.
    Area { from: f64, until: f64, to: f64 },
    /// A chosen peak's integral: y · to / the window's integral above its background (`window_value`, PL5).
    Peak { window: Window, background: WindowBackground, to: f64 },
    /// Every member's area over [from, until] made the reference member's (its factor 1).
    Reference { member: usize, from: f64, until: f64 },
}

/// One member normalised (PL6): its factor and its pattern, or why it has none.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Normalised {
    pub member: usize,
    pub factor: Option<f64>,
    pub pattern: Option<Pattern>,
    pub why_not: Option<String>,
}

fn area_over(p: &Pattern, from: f64, until: f64) -> Result<f64, String> {
    let w = Window { name: "area".into(), from, to: until };
    Ok(window_value(p, &vec![0.0; p.len()], false, &w, &[])?.integral)
}

/// Each member normalised (PL6): its factor and its pattern scaled (σ alike); a member the rule cannot normalise (no
/// monitor, no counting time, a zero or negative area or integral, a range outside it) has none, with the reason.
pub fn normalise(set: &PlotSet, rule: &Normalisation) -> Result<Vec<Normalised>, String> {
    let reference_area = match rule {
        Normalisation::Reference { member, from, until } => {
            let m = set.members.get(*member).ok_or(format!("no member {}", member + 1))?;
            let a = area_over(&m.pattern, *from, *until).map_err(|e| format!("the reference {}: {e}", m.label))?;
            if a.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater) {
                return Err(format!("the reference {} has no positive area over {from} to {until}", m.label));
            }
            Some(a)
        }
        _ => None,
    };
    Ok(set
        .members
        .iter()
        .enumerate()
        .map(|(k, m)| {
            let p = &m.pattern;
            let max_of = |from: f64, until: f64| p.x.iter().zip(&p.y).filter(|(x, _)| **x >= from && **x <= until).map(|(_, y)| *y).fold(f64::NEG_INFINITY, f64::max);
            let positive = |v: f64, what: &str| -> Result<f64, String> { if v > 0.0 && v.is_finite() { Ok(v) } else { Err(format!("{what} is {v}, not positive")) } };
            let factor: Result<f64, String> = match rule {
                Normalisation::None => Ok(1.0),
                Normalisation::Monitor { to } => m.facts.monitor.ok_or("the file states no monitor".to_string()).and_then(|v| positive(v, "the monitor")).map(|v| to / v),
                Normalisation::Time { to } => m.facts.counting_time.ok_or("the file states no counting time".to_string()).and_then(|v| positive(v, "the counting time")).map(|v| to / v),
                Normalisation::Maximum { to } => positive(max_of(f64::NEG_INFINITY, f64::INFINITY), "the largest y").map(|v| to / v),
                Normalisation::MaximumIn { from, until, to } => positive(max_of(*from, *until), "the largest y in the range").map(|v| to / v),
                Normalisation::Area { from, until, to } => area_over(p, *from, *until).and_then(|a| positive(a, "the area")).map(|a| to / a),
                Normalisation::Peak { window, background, to } => {
                    let b = match background {
                        WindowBackground::Auto { log_lambda } => crate::data_ops::auto_background(p, *log_lambda, crate::data_ops::DEFAULT_MAX_ITER).map(|b| b.curve),
                        WindowBackground::Points { points, interpolation } => Ok(crate::background::interpolate(points, *interpolation, &p.x)),
                        WindowBackground::Line => Ok(Vec::new()),
                    };
                    b.and_then(|b| window_value(p, &b, matches!(background, WindowBackground::Line), window, &[])).and_then(|v| positive(v.integral, "the peak's integral")).map(|i| to / i)
                }
                Normalisation::Reference { from, until, .. } => area_over(p, *from, *until).and_then(|a| positive(a, "the area")).map(|a| reference_area.unwrap() / a),
            };
            match factor {
                Ok(f) => Normalised { member: k, factor: Some(f), pattern: Some(Pattern { x: p.x.clone(), y: p.y.iter().map(|y| y * f).collect(), sigma: p.sigma.iter().map(|s| s * f.abs()).collect() }), why_not: None },
                Err(e) => Normalised { member: k, factor: None, pattern: None, why_not: Some(e) },
            }
        })
        .collect())
}

impl PlotSet {
    /// The set as the view shows it after a normalisation (PL6): every member scaled, the files and their operations
    /// untouched; refused, naming them, when some member cannot be normalised by the rule.
    pub fn normalised(&self, rule: &Normalisation) -> Result<PlotSet, String> {
        let n = normalise(self, rule)?;
        let missing: Vec<String> = n.iter().filter(|x| x.factor.is_none()).map(|x| format!("{}: {}", self.members[x.member].label, x.why_not.clone().unwrap_or_default())).collect();
        if !missing.is_empty() {
            return Err(format!("not every member can be normalised so: {}", missing.join("; ")));
        }
        let mut out = self.clone();
        for (m, x) in out.members.iter_mut().zip(n) {
            m.pattern = x.pattern.unwrap();
        }
        Ok(out)
    }

    /// The set as the view shows it after a normalisation, and the lines a file written from it carries to say so
    /// (PL25): the rule in words, then each member's factor. A file of normalised numbers that did not say how they
    /// were normalised would mislead whoever reads it. With `Normalisation::None`, the set as read and no lines.
    pub fn as_shown(&self, rule: &Normalisation) -> Result<(PlotSet, Vec<String>), String> {
        if *rule == Normalisation::None {
            return Ok((self.clone(), Vec::new()));
        }
        let shown = self.normalised(rule)?;
        let n = normalise(self, rule)?;
        let mut lines = vec![format!("y normalised as the view shows it (the files are as read): {}", rule.words(self))];
        for x in &n {
            lines.push(format!("factor of {}: {}", self.members[x.member].label, x.factor.unwrap()));
        }
        Ok((shown, lines))
    }
}

impl Normalisation {
    /// The rule in words, in the set's x (PL25).
    pub fn words(&self, set: &PlotSet) -> String {
        let x = set.members.first().map_or(Axis::TwoTheta, |m| m.facts.axis).label();
        match self {
            Normalisation::None => "none".into(),
            Normalisation::Monitor { to } => format!("y × {to} / the monitor count"),
            Normalisation::Time { to } => format!("y × {to} / the counting time"),
            Normalisation::Maximum { to } => format!("the largest y made {to}"),
            Normalisation::MaximumIn { from, until, to } => format!("the largest y over {x} {from} to {until} made {to}"),
            Normalisation::Area { from, until, to } => format!("the same area, {to}, over {x} {from} to {until}"),
            Normalisation::Peak { window, background, to } => {
                let b = match background {
                    WindowBackground::Auto { .. } => "the automatic background",
                    WindowBackground::Line => "the line between its ends",
                    WindowBackground::Points { .. } => "the given background points",
                };
                format!("the integral of the peak over {x} {} to {} above {b} made {to}", window.from, window.to)
            }
            Normalisation::Reference { member, from, until } => {
                format!("the area over {x} {from} to {until} made that of {}", set.members.get(*member).map_or(format!("member {}", member + 1), |m| m.label.clone()))
            }
        }
    }
}

// PL10: a set in another quantity, as a view -----------------------------------------------------------------------------

/// The x quantities a set is shown in (PL10; WinPLOTR's x spaces, deck p. 43: 2THETA, Q, STL, D, S).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum XQuantity {
    TwoTheta,
    D,
    Q,
    /// sin θ/λ (Å⁻¹).
    Stl,
    /// s = 1/d = 2 sin θ/λ (Å⁻¹).
    S,
}

impl XQuantity {
    pub fn word(self) -> &'static str {
        match self {
            XQuantity::TwoTheta => "2theta",
            XQuantity::D => "d",
            XQuantity::Q => "q",
            XQuantity::Stl => "stl",
            XQuantity::S => "s",
        }
    }

    pub fn from_word(w: &str) -> Option<XQuantity> {
        [XQuantity::TwoTheta, XQuantity::D, XQuantity::Q, XQuantity::Stl, XQuantity::S].into_iter().find(|q| q.word() == w)
    }

    pub fn label(self) -> &'static str {
        match self {
            XQuantity::TwoTheta => "2θ (°)",
            XQuantity::D => "d (Å)",
            XQuantity::Q => "Q (Å⁻¹)",
            XQuantity::Stl => "sin θ/λ (Å⁻¹)",
            XQuantity::S => "s = 1/d (Å⁻¹)",
        }
    }
}

/// Where a member's λ comes from (PL10).
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum LambdaFrom {
    /// The data file states it.
    File,
    /// A control file's λ1, by its name.
    ControlFile(String),
    /// Typed, with the user's note.
    Typed(String),
}

/// One member in another quantity (PL10): the λ used and its source (none where the conversion needs none), its
/// pattern, or why it has none.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct XMember {
    pub member: usize,
    pub lambda: Option<f64>,
    pub from: Option<LambdaFrom>,
    pub pattern: Option<Pattern>,
    pub why_not: Option<String>,
}

/// A set shown in another quantity (PL10): every member, and those without a λ listed.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct XView {
    pub quantity: XQuantity,
    pub members: Vec<XMember>,
    pub without: Vec<usize>,
}

/// The set shown in another x quantity (PL10), as a view: no file and no operation changed. Each member's λ is the one
/// given for it, else its file's, else `fallback` (a control file's or a typed one); a member that needs one (in 2θ,
/// or shown in 2θ) and has none is listed in `without`. Through sin θ/λ: from 2θ sin(θ)/λ, from d 1/(2d), from Q Q/4π;
/// to d 1/(2 sin θ/λ), Q 4π sin θ/λ, s 2 sin θ/λ, 2θ 2 asin(λ sin θ/λ). Points with no value there (d of sin θ/λ ≤ 0,
/// 2θ of λ sin θ/λ > 1) are left out; the points sorted by x; y and σ unchanged. TOF and energy are not converted.
pub fn x_view(set: &PlotSet, quantity: XQuantity, given: &[(usize, f64, LambdaFrom)], fallback: Option<(f64, LambdaFrom)>) -> XView {
    let mut members = Vec::new();
    let mut without = Vec::new();
    for (k, m) in set.members.iter().enumerate() {
        let lambda: Option<(f64, LambdaFrom)> = given.iter().find(|g| g.0 == k).map(|g| (g.1, g.2.clone())).or(m.facts.wavelength.map(|l| (l, LambdaFrom::File))).or(fallback.clone());
        let same = matches!((m.facts.axis, quantity), (Axis::TwoTheta, XQuantity::TwoTheta) | (Axis::D, XQuantity::D) | (Axis::Q, XQuantity::Q));
        let needs = !same && (m.facts.axis == Axis::TwoTheta || quantity == XQuantity::TwoTheta);
        let mut x = XMember { member: k, lambda: None, from: None, pattern: None, why_not: None };
        if same {
            x.pattern = Some(m.pattern.clone());
        } else if !matches!(m.facts.axis, Axis::TwoTheta | Axis::D | Axis::Q) {
            x.why_not = Some(format!("its x is {}, which is not converted", m.facts.axis.label()));
        } else if needs && lambda.is_none() {
            x.why_not = Some("no wavelength: the file states none, and none is given".into());
            without.push(k);
        } else {
            let l = lambda.as_ref().map(|l| l.0);
            if needs {
                (x.lambda, x.from) = (l, lambda.map(|l| l.1));
            }
            let stl = |v: f64| -> f64 {
                match m.facts.axis {
                    Axis::TwoTheta => (v / 2.0).to_radians().sin() / l.unwrap(),
                    Axis::D => 1.0 / (2.0 * v),
                    _ => v / (4.0 * std::f64::consts::PI),
                }
            };
            let to = |t: f64| -> Option<f64> {
                let v = match quantity {
                    XQuantity::Stl => t,
                    XQuantity::S => 2.0 * t,
                    XQuantity::Q => 4.0 * std::f64::consts::PI * t,
                    XQuantity::D => (t > 0.0).then(|| 1.0 / (2.0 * t))?,
                    XQuantity::TwoTheta => {
                        let a = l.unwrap() * t;
                        (a.abs() <= 1.0).then(|| 2.0 * a.asin().to_degrees())?
                    }
                };
                v.is_finite().then_some(v)
            };
            let p = &m.pattern;
            let mut pts: Vec<(f64, f64, f64)> = (0..p.len()).filter_map(|i| to(stl(p.x[i])).map(|x| (x, p.y[i], p.sigma[i]))).collect();
            pts.sort_by(|a, b| a.0.total_cmp(&b.0));
            x.pattern = Some(Pattern { x: pts.iter().map(|q| q.0).collect(), y: pts.iter().map(|q| q.1).collect(), sigma: pts.iter().map(|q| q.2).collect() });
        }
        members.push(x);
    }
    XView { quantity, members, without }
}

// W2: the series a run's files offer ----------------------------------------------------------------------------------

/// A set the files offer (W2): what it is (`kind` "by_temperature", "hold" or "against_time"), a title that says it and
/// how many scans, the exposure class it is of (a hold's, when the run has two or more; none for a set across classes),
/// the set file the host opens as it is, the scans it leaves out, with why, and the reader every member was read with
/// (W9: given back when the host opens the set, the files are not identified again).
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Proposal {
    pub kind: String,
    pub title: String,
    pub class: Option<usize>,
    pub set: SetFile,
    pub left_out: Vec<(String, String)>,
    #[cfg_attr(feature = "serde", serde(default))]
    pub reader: Option<String>,
}

/// An exposure class (W2, W9): scans whose counts' scale, the median of σ²/y, lies within `CLASS_STEP` of its
/// neighbours' once sorted; the classes in increasing scale, the longest exposure (the most counts per unit of y)
/// first.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ExposureClass {
    pub scans: usize,
    /// The smallest and the largest median of σ²/y of its scans.
    pub median_scale: (f64, f64),
}

/// What `propose` found: the sets offered, the exposure classes, and the scans in no set, with why.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Proposals {
    pub proposals: Vec<Proposal>,
    pub classes: Vec<ExposureClass>,
    pub left_out: Vec<(String, String)>,
}

/// One scan's facts for `propose`.
struct Fact {
    file: String,
    scan: usize,
    scans: usize,
    label: String,
    /// In the run's unit (`unit`).
    temperature: Option<Variable>,
    /// Seconds, and where they come from.
    time: Option<(f64, Found)>,
    scale: f64,
    class: usize,
    points: usize,
    two_theta: bool,
    format: String,
}

impl Fact {
    /// The scan as `left_out` names it: the file as given, and its scan when the file holds several.
    fn named(&self) -> String {
        if self.scans > 1 {
            format!("{} ({})", self.file, self.scan + 1)
        } else {
            self.file.clone()
        }
    }
}

/// Two files' modification times closer than this (seconds) are one write, a copy's or a checkout's, not two scans'.
pub const SAME_WRITE: f64 = 1.0;

/// Two temperatures are one when within this (in the run's unit, K or °C): a hold's scans read 1270.0 and 1270.3.
pub const SAME_TEMPERATURE: f64 = 0.5;

/// A new exposure class where a scan's counts' scale (the median of σ²/y) is more than this times the one before, once
/// sorted (W9, L428). For counts N scaled by c (y = cN, σ = c√N), σ²/y = c whatever the intensity, so a class is one
/// scale. Measured: in a team's ESRF run (101 scans) the step between the long and the short
/// exposure is 1.265 and between the short and scan 085 2.195, and no step inside a class exceeds 1.0075; Tb₂Pd₂Sn's
/// 20 raw-count scans have σ²/y = 1 each, one class. The relative error σ/y that F02 compares (0.040 against 0.047) is
/// not used: it moves with the intensity, and Tb₂Pd₂Sn's 20 scans of one exposure step by up to 1.107 in it, against
/// the ESRF classes' 1.118.
pub const CLASS_STEP: f64 = 1.1;

/// A proposal whose members have more points than this is rebinned on opening (W9), to `LARGE_STEP` in 2θ.
pub const LARGE_SCAN: usize = 5000;

/// The step a large set is rebinned to, in degrees 2θ (W9; the step the ESRF run's maps were made on).
pub const LARGE_STEP: f64 = 0.01;

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(|a, b| a.total_cmp(b));
    let n = v.len();
    if n == 0 {
        f64::NAN
    } else if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    }
}

fn rounded(v: f64) -> f64 {
    (v * 1000.0).round() / 1000.0
}

/// The unit a span of time is shown in: seconds under 2 min, minutes under 5 h, hours above (the ESRF run's sets
/// built by hand: the 2.8 h hold in minutes, the 10 h run in hours).
fn time_unit(span: f64) -> (&'static str, f64) {
    if span < 120.0 {
        ("s", 1.0)
    } else if span < 18000.0 {
        ("min", 60.0)
    } else {
        ("h", 3600.0)
    }
}

/// A scan's counts' scale: the median of σ²/y over its points with y > 0 and σ > 0, above 2° in 2θ (below, the direct
/// beam's tail) or all of them on another axis.
fn counts_scale(m: &Member) -> f64 {
    let p = &m.pattern;
    let tt = m.facts.axis == crate::formats::Axis::TwoTheta;
    let mut v: Vec<f64> = (0..p.len()).filter(|&i| p.y[i] > 0.0 && p.sigma[i] > 0.0 && p.sigma[i].is_finite() && (!tt || p.x[i] > 2.0)).map(|i| p.sigma[i] * p.sigma[i] / p.y[i]).collect();
    median(&mut v)
}

/// The series a run's files offer (W2, I10 step 3; W8, W9): from each file's facts as the registry reads them (one scan
/// per member; its temperature from PL2's sources, in °C when every one comes from names that say °C, else in K; its
/// time from the elapsed seconds its name gives, else the modification time the host gives, `OpenedFile::modified`,
/// unless two files' fall within `SAME_WRITE` (a copy's); its exposure class by its counts' scale, `CLASS_STEP`):
/// - by temperature: one scan per temperature (within `SAME_TEMPERATURE`), of the longest exposure that has one there
///   and the last of those in time order (the nearest to equilibrium; in the files' order without times), when three
///   or more temperatures are found;
/// - against time: every scan with a time, all classes, against the time since the first; when most scans give a
///   temperature, one that gives none (a background) is left out and named;
/// - the hold: in each class, the longest run of three or more scans in time order at one temperature, against the
///   time since its first; the class's other holds named in its title.
///
/// Each proposal is a set file the host opens as it is: its title says what it is and how many scans; a temperature
/// that a re-read gives is not typed, a time is typed with its source, and the temperature is a time set's members'
/// other variable (W3) with its source. A proposal whose members have more than `LARGE_SCAN` points in 2θ is rebinned
/// to `LARGE_STEP` on opening (`SetFile::rebin`), its title saying so. A file the registry cannot read, or a scan in
/// no set, is named with why: by its file as given, and its scan in brackets when the file holds several.
pub fn propose(files: &[OpenedFile]) -> Proposals {
    let mut left_out: Vec<(String, String)> = Vec::new();
    let mut facts: Vec<Fact> = Vec::new();
    for f in files {
        match open_from(std::slice::from_ref(f), None, None) {
            Err(e) => left_out.push((f.name.clone(), format!("not read: {e}"))),
            Ok(set) => {
                for m in &set.members {
                    let time = elapsed_in_name(&m.file).map(|(v, t)| (v, Found::FileName(t))).or(m.modified.map(|t| (t, Found::Typed("the file's modification time".into()))));
                    facts.push(Fact {
                        file: m.file.clone(),
                        scan: m.scan,
                        scans: m.scans,
                        label: m.label.clone(),
                        temperature: m.variable.clone(),
                        time,
                        scale: counts_scale(m),
                        class: 0,
                        points: m.pattern.len(),
                        two_theta: m.facts.axis == crate::formats::Axis::TwoTheta,
                        format: m.facts.format.clone(),
                    });
                }
            }
        }
    }
    // The run's unit: °C when every temperature whose source gives a unit gives °C (W8), else K, each converted.
    let units: Vec<&str> = facts.iter().filter_map(|f| f.temperature.as_ref().and_then(|v| v.found.unit())).collect();
    let unit = if !units.is_empty() && units.iter().all(|u| *u == "°C") { "°C" } else { "K" };
    for f in facts.iter_mut() {
        f.temperature = f.temperature.as_ref().map(|v| v.in_unit(unit));
    }
    // Modification times say when a scan was written only when no two files' fall within `SAME_WRITE` of each other:
    // files checked out or copied together carry the time of the copy.
    let mut written: Vec<(&str, f64)> = Vec::new();
    for f in &facts {
        if let Some((t, Found::Typed(_))) = &f.time {
            if !written.iter().any(|w| w.0 == f.file) {
                written.push((&f.file, *t));
            }
        }
    }
    written.sort_by(|a, b| a.1.total_cmp(&b.1));
    if written.windows(2).any(|w| w[1].1 - w[0].1 < SAME_WRITE) {
        for f in facts.iter_mut() {
            if matches!(f.time, Some((_, Found::Typed(_)))) {
                f.time = None;
            }
        }
    }
    // The exposure classes: the scales sorted, a new class where one is more than `CLASS_STEP` times the one before.
    let mut order: Vec<usize> = (0..facts.len()).filter(|&i| facts[i].scale.is_finite()).collect();
    order.sort_by(|&a, &b| facts[a].scale.total_cmp(&facts[b].scale));
    let mut classes: Vec<ExposureClass> = Vec::new();
    let mut last = f64::NAN;
    for &i in &order {
        let m = facts[i].scale;
        if classes.is_empty() || m > CLASS_STEP * last {
            classes.push(ExposureClass { scans: 0, median_scale: (m, m) });
        }
        let c = classes.len() - 1;
        classes[c].scans += 1;
        classes[c].median_scale.1 = m;
        facts[i].class = c;
        last = m;
    }
    let several = classes.len() > 1;
    let class_words = |c: usize| if several { format!(" (exposure class {} of {}, σ²/y {:.3} to {:.3})", c + 1, classes.len(), classes[c].median_scale.0, classes[c].median_scale.1) } else { String::new() };
    let across = if several { format!(" ({} exposure classes)", classes.len()) } else { String::new() };
    // The scans with an esd, in time order when every one has a time, else in the files' order.
    let mut seq: Vec<usize> = order.clone();
    seq.sort();
    if seq.iter().all(|&i| facts[i].time.is_some()) {
        seq.sort_by(|&a, &b| facts[a].time.as_ref().unwrap().0.total_cmp(&facts[b].time.as_ref().unwrap().0));
    }
    let time_of = |i: usize| facts[i].time.as_ref().unwrap().0;
    // Each member against the time since `t0`, its temperature the other variable with its source.
    let timed_member = |i: usize, t0: f64, first: &str, per: f64| -> SetMember {
        let f = &facts[i];
        let (t, found) = f.time.clone().unwrap();
        let mut m = SetMember { file: f.file.clone(), scan: (f.scans > 1).then_some(f.scan + 1), ..Default::default() };
        m.value = Some(((t - t0) / per, Found::Typed(format!("the time since {first}, from {}", found.words()))));
        if let Some(v) = &f.temperature {
            m.variables.push(("T".into(), v.value, format!("{} ({unit})", v.found.words())));
        }
        m
    };
    let time_set = |members: Vec<SetMember>, per: &str, title: String| SetFile { title: Some(title), variable_name: "t".into(), variable_unit: per.into(), order: Order::Variable, members, ..Default::default() };
    let mut proposals: Vec<Proposal> = Vec::new();
    // By temperature: the temperatures in the order first met, each scan in the group of the first within
    // `SAME_TEMPERATURE` of it; the longest exposure there, its last scan.
    let mut groups: Vec<(f64, Vec<usize>)> = Vec::new();
    let mut by_t_out: Vec<(String, String)> = Vec::new();
    for &i in &seq {
        match facts[i].temperature.as_ref().map(|v| v.value) {
            None => by_t_out.push((facts[i].named(), "no temperature: neither its header, a buffer's comment, its title nor its name gives one".into())),
            Some(v) => match groups.iter_mut().find(|g| (g.0 - v).abs() <= SAME_TEMPERATURE) {
                Some(g) => g.1.push(i),
                None => groups.push((v, vec![i])),
            },
        }
    }
    if groups.len() >= 3 {
        let mut picked: Vec<usize> = Vec::new();
        for (_, g) in &groups {
            let best = g.iter().map(|&i| facts[i].class).min().unwrap();
            let p = *g.iter().rev().find(|&&i| facts[i].class == best).unwrap();
            for &i in g.iter().filter(|&&i| i != p) {
                let v = facts[i].temperature.as_ref().unwrap().value;
                let why = if several { "the last scan of the longest exposure there is taken" } else { "the last scan there is taken" };
                by_t_out.push((facts[i].named(), format!("{v} {unit} is {}'s: {why}", facts[p].label)));
            }
            picked.push(p);
        }
        picked.sort_by(|&a, &b| facts[a].temperature.as_ref().unwrap().value.total_cmp(&facts[b].temperature.as_ref().unwrap().value));
        let (lo_t, hi_t) = (facts[picked[0]].temperature.as_ref().unwrap().value, facts[picked[picked.len() - 1]].temperature.as_ref().unwrap().value);
        let members: Vec<SetMember> = picked
            .iter()
            .map(|&i| {
                let f = &facts[i];
                let mut m = SetMember { file: f.file.clone(), scan: (f.scans > 1).then_some(f.scan + 1), ..Default::default() };
                // a temperature a re-read gives is not typed; one typed for the file keeps its source
                if let Some(Variable { value, found: found @ Found::Typed(_) }) = &f.temperature {
                    m.value = Some((*value, found.clone()));
                }
                m
            })
            .collect();
        let classes_in: Vec<usize> = picked.iter().map(|&i| facts[i].class).collect();
        let one = classes_in.iter().all(|&c| c == classes_in[0]);
        let words = if !several {
            String::new()
        } else if one {
            class_words(classes_in[0])
        } else {
            format!(" (the longest exposure at each temperature, of {} exposure classes)", classes.len())
        };
        let title = format!("{} scans by temperature, {lo_t} to {hi_t} {unit}{words}", members.len());
        proposals.push(Proposal { kind: "by_temperature".into(), title: title.clone(), class: (several && one).then_some(classes_in[0]), set: SetFile { title: Some(title), variable_unit: unit.into(), members, ..Default::default() }, left_out: by_t_out, reader: None });
    }
    // Against time: every scan with a time; when most give a temperature, those that give none left out.
    let with_temperature = seq.iter().filter(|&&i| facts[i].temperature.is_some()).count();
    let most = 2 * with_temperature > seq.len();
    let run: Vec<usize> = seq.iter().copied().filter(|&i| facts[i].time.is_some() && (facts[i].temperature.is_some() || !most)).collect();
    if run.len() >= 3 {
        let t0 = time_of(run[0]);
        let span = time_of(run[run.len() - 1]) - t0;
        let per = time_unit(span);
        let first = facts[run[0]].label.clone();
        let members: Vec<SetMember> = run.iter().map(|&i| timed_member(i, t0, &first, per.1)).collect();
        let mut lo: Vec<(String, String)> = Vec::new();
        for &i in &seq {
            if facts[i].time.is_none() {
                lo.push((facts[i].named(), "no time: its name gives no elapsed seconds and the host gave no modification time".into()));
            } else if !run.contains(&i) {
                lo.push((facts[i].named(), "no temperature, where the run's other scans give one: left out of the run".into()));
            }
        }
        let title = format!("{} scans against time, {} {} from {first}{across}", members.len(), rounded(span / per.1), per.0);
        let one = run.iter().all(|&i| facts[i].class == facts[run[0]].class);
        proposals.push(Proposal { kind: "against_time".into(), title: title.clone(), class: (several && one).then_some(facts[run[0]].class), set: time_set(members, per.0, title), left_out: lo, reader: None });
    }
    // The hold of each class: the longest run at one temperature, in time order.
    for class in 0..classes.len() {
        let mine: Vec<usize> = seq.iter().copied().filter(|&i| facts[i].class == class && facts[i].time.is_some()).collect();
        let mut holds: Vec<(usize, usize)> = Vec::new();
        let mut k = 0;
        while k < mine.len() {
            let Some(t) = facts[mine[k]].temperature.as_ref().map(|v| v.value) else {
                k += 1;
                continue;
            };
            let mut e = k + 1;
            while e < mine.len() && facts[mine[e]].temperature.as_ref().is_some_and(|v| (v.value - t).abs() <= SAME_TEMPERATURE) {
                e += 1;
            }
            if e - k >= 3 {
                holds.push((k, e));
            }
            k = e;
        }
        let Some(&(k, e)) = holds.iter().max_by(|a, b| (a.1 - a.0).cmp(&(b.1 - b.0)).then(b.0.cmp(&a.0))) else { continue };
        let held = &mine[k..e];
        let t = facts[held[0]].temperature.as_ref().unwrap().value;
        let (hs, he) = (time_of(held[0]), time_of(held[held.len() - 1]));
        let per = time_unit(he - hs);
        let first = facts[held[0]].label.clone();
        let members: Vec<SetMember> = held.iter().map(|&i| timed_member(i, hs, &first, per.1)).collect();
        let others: Vec<String> = holds.iter().filter(|h| h.0 != k).map(|&(a, b)| format!("{} {unit} ({} scans)", facts[mine[a]].temperature.as_ref().unwrap().value, b - a)).collect();
        let others = if others.is_empty() { String::new() } else { format!("; the class's other holds: {}", others.join(", ")) };
        let title = format!("the hold at {t} {unit}: {} scans against time, {} {} from {first}{}{others}", members.len(), rounded((he - hs) / per.1), per.0, class_words(class));
        proposals.push(Proposal { kind: "hold".into(), title: title.clone(), class: several.then_some(class), set: time_set(members, per.0, title), left_out: Vec::new(), reader: None });
    }
    // A large set rebinned on opening (W9), its title saying so; the reader its members share.
    for p in proposals.iter_mut() {
        let formats: Vec<&str> = p.set.members.iter().filter_map(|m| facts.iter().find(|f| f.file == m.file).map(|f| f.format.as_str())).collect();
        p.reader = formats.first().filter(|f| formats.iter().all(|g| g == *f)).map(|f| f.to_string());
        let big = p.set.members.iter().any(|m| facts.iter().any(|f| f.file == m.file && f.points > LARGE_SCAN));
        let tt = p.set.members.iter().all(|m| facts.iter().any(|f| f.file == m.file && f.two_theta));
        if big && tt {
            p.set.rebin = Some(LARGE_STEP);
            p.title = format!("{}, rebinned to {LARGE_STEP}°", p.title);
            p.set.title = Some(p.title.clone());
        }
    }
    for f in &facts {
        if !proposals.iter().any(|p| p.set.members.iter().any(|m| m.file == f.file && m.scan.is_none_or(|k| k == f.scan + 1))) {
            let why = if !f.scale.is_finite() {
                "no esd: its points carry no weight".to_string()
            } else if f.temperature.is_none() && f.time.is_none() {
                "no temperature and no time: in no set".to_string()
            } else if f.temperature.is_none() && most {
                "no temperature, where the run's other scans give one: in no set".to_string()
            } else {
                "in no set: too few scans with a time or temperatures".to_string()
            };
            left_out.push((f.named(), why));
        }
    }
    Proposals { proposals, classes, left_out }
}

