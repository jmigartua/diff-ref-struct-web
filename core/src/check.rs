//! The core's checks (R24, L076): what is wrong with a control file or a `.pcr`, its model and its refinement,
//! found before a refinement runs, each said where it is, what is wrong, why, and what to do.
//!
//! Four levels, each run when the one before it leaves something to check:
//! - **file**: the syntax (the control file's parser, L066), unknown names (with the nearest known one as a
//!   replacement), the units written in comments against the dictionary's, values outside the engine's limits (P39);
//! - **model**: the data and resolution files, the wavelength, the angular range against the data, occupancies beyond
//!   a site, values a site's symmetry fixes but the refinement sets free, phases without reflections in the data;
//! - **refinement**: holes and wrong counts in a `.pcr`'s parameter numbers ("Hole in the matrix", L071), ties or
//!   constraints that do not build, free parameters that do not change the pattern (a zero column of the Jacobian);
//! - **run**: the normal matrix at the starting values (`lsq::normal_matrix`), its null and near-null directions
//!   (eigen-decomposition of the matrix scaled to unit diagonal, as the engine solves it), each named as the
//!   parameters it combines, and the correlations it predicts above a threshold. A pair or group is explained by the
//!   parameter dictionary's `correlations` (why, tip, source; `parameters.yml`), never by a text written here.
//!
//! Each diagnostic has a stable `code`; the app calls `check_path` (or `check_control`, `check_pcr`) as `drs check`
//! does, so that both report the same problems.

use crate::control::syntax::{self, Span};
use crate::control::schema;
use crate::data::Pattern;
use crate::jacobian::Field;
use crate::parameter_table::Correlation;
use crate::powder::{Intensities, Model};
use crate::refinable::{address, dictionary_id, limits, ParameterSet};
use std::path::Path;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Severity {
    Error,
    Warning,
    Hint,
}

impl Severity {
    pub fn name(self) -> &'static str {
        match self {
            Severity::Error => "error",
            Severity::Warning => "warning",
            Severity::Hint => "hint",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Level {
    File,
    Model,
    Refinement,
    Run,
}

impl Level {
    pub fn name(self) -> &'static str {
        match self {
            Level::File => "file",
            Level::Model => "model",
            Level::Refinement => "refinement",
            Level::Run => "run",
        }
    }
}

/// One problem: where it is (the file's line and column, or the model's addresses), what is wrong, why, what to do
/// (and the text to put in place of the span, when the fix is mechanical), and the dictionary entry it concerns.
#[derive(Clone, Debug, PartialEq)]
pub struct Diagnostic {
    pub severity: Severity,
    pub level: Level,
    /// Stable: `syntax`, `unknown-name`, `value`, `missing`, `unit`, `range`, `data-missing`, `irf-missing`,
    /// `wavelength`, `range-data`, `occupancy`, `symmetry-fixed`, `no-reflections`, `model`, `parameter-hole`,
    /// `parameter-count`, `refinement`, `no-effect`, `singular`, `near-singular`, `correlation`.
    pub code: &'static str,
    pub span: Option<Span>,
    pub addresses: Vec<String>,
    pub what: String,
    pub why: String,
    pub suggestion: Option<String>,
    pub replacement: Option<String>,
    /// The dictionary ids it concerns.
    pub entries: Vec<&'static str>,
}

impl Diagnostic {
    pub(crate) fn new(severity: Severity, level: Level, code: &'static str, what: impl Into<String>, why: impl Into<String>) -> Diagnostic {
        Diagnostic { severity, level, code, span: None, addresses: Vec::new(), what: what.into(), why: why.into(), suggestion: None, replacement: None, entries: Vec::new() }
    }

    pub(crate) fn at(mut self, span: Option<Span>) -> Self {
        self.span = span;
        self
    }

    pub(crate) fn addr(mut self, a: impl Into<String>) -> Self {
        self.addresses.push(a.into());
        self
    }

    pub(crate) fn fix(mut self, s: impl Into<String>) -> Self {
        self.suggestion = Some(s.into());
        self
    }

    pub(crate) fn entry(mut self, id: &'static str) -> Self {
        if !id.is_empty() && !self.entries.contains(&id) {
            self.entries.push(id);
        }
        self
    }

    /// `file:line:column` (or `file` and the addresses).
    pub fn place(&self, file: &str) -> String {
        match self.span {
            Some(s) if s.line > 0 => format!("{file}:{}:{}", s.line, s.column),
            _ if !self.addresses.is_empty() => format!("{file}: {}", self.addresses.join(", ")),
            _ => file.to_string(),
        }
    }

    /// The diagnostic as `drs check` prints it.
    pub fn to_text(&self, file: &str) -> String {
        let mut s = format!("{}: {}[{}]: {}", self.place(file), self.severity.name(), self.code, self.what);
        if !self.why.is_empty() {
            s.push_str(&format!("\n    why: {}", self.why));
        }
        if let Some(f) = &self.suggestion {
            s.push_str(&format!("\n    fix: {f}"));
        }
        if let Some(r) = &self.replacement {
            s.push_str(&format!("\n    replace with: {r}"));
        }
        if !self.entries.is_empty() {
            s.push_str(&format!("\n    see: {}", self.entries.join(", ")));
        }
        s
    }

    /// The diagnostic as a JSON object (`drs check --json`; the app reads the same fields).
    pub fn to_json(&self, file: &str) -> String {
        let q = |s: &str| {
            let mut o = String::from("\"");
            for c in s.chars() {
                match c {
                    '"' => o.push_str("\\\""),
                    '\\' => o.push_str("\\\\"),
                    '\n' => o.push_str("\\n"),
                    c if (c as u32) < 0x20 => o.push_str(&format!("\\u{:04x}", c as u32)),
                    c => o.push(c),
                }
            }
            o.push('"');
            o
        };
        let opt = |v: &Option<String>| v.as_deref().map_or("null".to_string(), q);
        let list = |v: &[String]| format!("[{}]", v.iter().map(|s| q(s)).collect::<Vec<_>>().join(", "));
        let span = self.span.map_or("null".to_string(), |s| format!("{{\"line\": {}, \"column\": {}, \"start\": {}, \"end\": {}}}", s.line, s.column, s.start, s.end));
        format!(
            "{{\"file\": {}, \"severity\": {}, \"level\": {}, \"code\": {}, \"span\": {span}, \"addresses\": {}, \"what\": {}, \"why\": {}, \"suggestion\": {}, \"replacement\": {}, \"entries\": {}}}",
            q(file),
            q(self.severity.name()),
            q(self.level.name()),
            q(self.code),
            list(&self.addresses),
            q(&self.what),
            q(&self.why),
            opt(&self.suggestion),
            opt(&self.replacement),
            list(&self.entries.iter().map(|s| s.to_string()).collect::<Vec<_>>())
        )
    }
}

/// What the checks take.
#[derive(Clone, Copy, Debug)]
pub struct Options {
    /// Predicted correlations at or above this (absolute value) are reported.
    pub correlation: f64,
    /// Run the run level (the normal matrix): it calculates the derivatives once.
    pub run: bool,
}

impl Default for Options {
    fn default() -> Self {
        Options { correlation: 0.9, run: true }
    }
}

/// An eigenvalue of the scaled normal matrix below this fraction of the largest is a null direction (the engine
/// leaves such directions out, `lsq::Solver`); below `NEAR` a near-null one.
pub const NULL: f64 = 1e-12;
pub const NEAR: f64 = 1e-8;

/// Check a control file (`.drc`) or a `.pcr`, by its extension.
pub fn check_path(path: &Path, opts: &Options) -> Vec<Diagnostic> {
    check_path_with_data(path, None, opts)
}

/// As `check_path`, a `.pcr` with the data file FullProf would be given on its command line (`fp2k CODFIL DATFIL`,
/// L079); a control file names its own data, and `data` is not used for it.
pub fn check_path_with_data(path: &Path, data: Option<&Path>, opts: &Options) -> Vec<Diagnostic> {
    let ext = path.extension().map(|e| e.to_string_lossy().to_lowercase()).unwrap_or_default();
    let bytes = match std::fs::read(path) {
        Ok(b) => b,
        Err(e) => return vec![Diagnostic::new(Severity::Error, Level::File, "missing", format!("{}: {e}", path.display()), "the file cannot be read")],
    };
    let text: String = String::from_utf8(bytes.clone()).unwrap_or_else(|_| bytes.iter().map(|&b| b as char).collect());
    let dir = path.parent().map(|d| d.to_path_buf()).unwrap_or_default();
    let dir = std::fs::canonicalize(&dir).unwrap_or(dir);
    match ext.as_str() {
        "pcr" => check_pcr_with_data(&text, &dir, path.file_stem().and_then(|s| s.to_str()).unwrap_or(""), data, opts),
        _ => check_control(&text, &dir, opts),
    }
}

// Level (a): the file ------------------------------------------------------------------------------------------------

/// `tof`: the file's pattern is a time-of-flight bank, whose 2θ quantities (its zero, its excluded regions) are in µs
/// (the app's request 46).
fn file_level(doc: &syntax::Document, read: &[syntax::Diagnostic], tof: bool) -> Vec<Diagnostic> {
    let mut out = Vec::new();
    for d in read {
        let (code, why) = match d.code {
            "syntax" => ("syntax", "the control file is a strict subset of TOML 1.0 (L066): one key = value a line, tables in [ ], strings in quotes"),
            "unknown-name" => ("unknown-name", "the names are the parameter dictionary's ids and the format's tables (drs check lists the nearest one)"),
            "missing" => ("missing", "the format needs this key to build the model"),
            _ => ("value", "the value does not have the type or the form the format gives this key"),
        };
        let mut x = Diagnostic::new(if d.warning { Severity::Warning } else { Severity::Error }, Level::File, code, d.message.clone(), why).at(Some(d.span));
        x.suggestion = Some(
            match code {
                "syntax" => "write the line as key = value, the value a number (1.5, 2e-3), true or false, a string in quotes, or a list in [ ]",
                "unknown-name" => "use a name of the parameter dictionary (the site's Parameters page) or of the format's tables",
                "missing" => "add the key to the table",
                _ => "give the value the type and form the message names",
            }
            .to_string(),
        );
        if let Some(r) = &d.replacement {
            x.replacement = Some(r.clone());
            x.suggestion = Some(format!("write \"{r}\""));
        }
        out.push(x);
    }
    // Units: a comment after a dictionary value that starts with a unit of the dictionary other than the entry's.
    let units: Vec<&str> = crate::parameters::dictionary().iter().map(|p| p.unit).filter(|u| !u.is_empty() && *u != "1").collect();
    for t in &doc.tables {
        for e in &t.entries {
            let Some(p) = crate::parameters::lookup(&e.key) else { continue };
            if p.unit.is_empty() || p.unit == "1" {
                continue;
            }
            // A TOF bank's zero and excluded regions are times of flight (`schema` writes them "# µs"); its other 2θ
            // values (displacement, transparency, the asymmetry's limit) are not used and keep their unit.
            let expected = if tof && matches!(p.id, "zero_shift" | "excluded_regions") { "µs" } else { p.unit };
            let rest = &doc.source[e.value.span.end..e.lines.1.min(doc.source.len())];
            let Some(h) = rest.find('#') else { continue };
            let comment = rest[h + 1..].trim();
            let written = units.iter().filter(|u| comment.starts_with(**u)).max_by_key(|u| u.len());
            if let Some(w) = written {
                if !comment.starts_with(expected) {
                    let start = e.value.span.end + h;
                    let end = e.value.span.end + rest[..].find('\n').unwrap_or(rest.len());
                    let col = e.value.span.column + doc.source[e.value.span.start..start].chars().count();
                    let span = Span { start, end, line: e.value.span.line, column: col };
                    let mut x = Diagnostic::new(Severity::Warning, Level::File, "unit", format!("\"{}\" is written in {w}; the file's unit for it is {expected}", e.key), "values are read in the dictionary's unit, whatever the comment says: a value meant in another unit is read wrongly").at(Some(span)).entry(p.id);
                    x.suggestion = Some(format!("give the value in {expected} and write the comment \"# {expected}\""));
                    x.replacement = Some(format!("# {expected}"));
                    out.push(x);
                }
            }
        }
    }
    out
}

/// The span of a field's value in a control file, where it can be found.
fn locate(doc: &syntax::Document, f: Field, m: &Model) -> Option<Span> {
    let id = dictionary_id(f, m);
    // The p-th [[phase]] table, and the a-th [[phase.atom]] after it.
    let (mut phase, mut atom) = (None::<usize>, 0usize);
    let mut table = None;
    for t in &doc.tables {
        if t.array && t.path == ["phase"] {
            phase = Some(phase.map_or(0, |k| k + 1));
            atom = 0;
        }
        let hit = match f {
            Field::Zero | Field::Sycos | Field::Sysin | Field::Ratio | Field::Wavelength | Field::Micro(_) => t.path == ["pattern"],
            Field::Background(_) | Field::BackgroundTerm(_) => t.path == ["pattern", "background"],
            Field::AtomX(p, a, _) | Field::AtomB(p, a) | Field::AtomOcc(p, a) | Field::AtomBeta(p, a, _) => t.array && t.path == ["phase", "atom"] && phase == Some(p) && atom == a,
            _ => f.phase().is_some_and(|p| t.array && t.path == ["phase"] && phase == Some(p)),
        };
        if t.array && t.path == ["phase", "atom"] {
            atom += 1;
        }
        if hit {
            table = Some(t);
            break;
        }
    }
    table?.get(id).map(|e| e.value.span)
}

// Level (b): the model -----------------------------------------------------------------------------------------------

fn model_level(m: &Model, data: Option<&Pattern>, out: &mut Vec<Diagnostic>) {
    // The wavelength; a TOF bank has none (its d from the time of flight, `tof::Bank`; the app's request 46).
    let tof = m.tof.is_some();
    let l1 = m.lambda[0];
    if tof {
    } else if l1.is_nan() || l1 <= 0.0 {
        out.push(Diagnostic::new(Severity::Error, Level::Model, "wavelength", format!("the wavelength λ1 is {l1} Å"), "every peak position follows from λ = 2d sin θ: λ must be positive").addr("pattern1.wavelength_1").entry("wavelength_1").fix("give the wavelength of the radiation in Å (Cu Kα1 1.5406, Mo Kα1 0.7093)"));
    } else if !(0.05..=10.0).contains(&l1) {
        out.push(Diagnostic::new(Severity::Warning, Level::Model, "wavelength", format!("the wavelength λ1 is {l1} Å"), "laboratory, synchrotron and neutron powder diffraction use 0.1–5 Å; a value outside suggests another unit (nm, pm) or a typing error").addr("pattern1.wavelength_1").entry("wavelength_1").fix("give the wavelength in Å"));
    }
    // Phases without reflections in the data's range, and the range of the data.
    if let Some(d) = data {
        if !d.is_empty() {
            let (lo, hi) = (d.x.iter().copied().fold(f64::MAX, f64::min), d.x.iter().copied().fold(f64::MIN, f64::max));
            if tof && lo <= 0.0 {
                out.push(Diagnostic::new(Severity::Error, Level::Model, "range-data", format!("the data run from {lo} to {hi} µs, not all positive"), "a time-of-flight pattern's x is the time of flight in µs, positive").addr("pattern1.data").fix("read the data with the reader of its format (drs data identify); a GSAS bank's CONST times are in units of 10 µs"));
            } else if !tof && (lo < 0.0 || hi > 180.0) {
                out.push(Diagnostic::new(Severity::Error, Level::Model, "range-data", format!("the data run from {lo} to {hi}, outside 0–180° 2θ"), "a constant-wavelength pattern's x is 2θ in degrees").addr("pattern1.data").fix("read the data with the reader of its format (drs data identify), or convert its x to 2θ (drs data op … convert)"));
            } else if tof || l1 > 0.0 {
                let peaks = m.peaks(lo, hi);
                for (p, ph) in m.phases.iter().enumerate() {
                    if !peaks.iter().any(|k| k.phase == p && k.two_theta >= lo && k.two_theta <= hi) {
                        let d = Diagnostic::new(Severity::Error, Level::Model, "no-reflections", format!("phase {} ({}) has no reflection between {lo:.3}{u} and {hi:.3}{u}", p + 1, ph.name.trim(), u = if tof { " µs" } else { "°" }), if tof { "with the cell and the bank's constants given, every d of the phase falls outside the data: its scale and structure cannot be refined" } else { "with the cell and the wavelength given, every d of the phase falls outside the data: its scale and structure cannot be refined" }).addr(format!("phase{}", p + 1)).entry("cell_a");
                        out.push(if tof { d.entry("tof_dtt1").fix("check the cell (Å) and Dtt1 (µs/Å) against the data") } else { d.entry("wavelength_1").fix("check the cell (Å) and the wavelength (Å) against the data") });
                    }
                }
            }
        }
    }
    // Occupancies beyond a site: the chemical occupancies of the atoms at one position add up to more than 1.
    for (p, ph) in m.phases.iter().enumerate() {
        let Intensities::Structure(st, _) = &ph.intensities else { continue };
        // FullProf's Occ is the chemical occupancy times m/M, or atoms per primitive cell when that is the only
        // reading that keeps every atom at 1 or below (`bonds::chemical_occupancies`; LaB6rx.pcr, ybfluo_Riet_IRF.pcr).
        let Ok(c1) = crate::bonds::chemical_occupancies(st) else { continue };
        let mut groups: Vec<(usize, [f64; 3], f64, Vec<usize>)> = Vec::new();
        for (a, at) in st.atoms.iter().enumerate() {
            let c = c1[a];
            let near = |u: [f64; 3], v: [f64; 3]| (0..3).all(|k| {
                let d = (u[k] - v[k]).rem_euclid(1.0);
                d.min(1.0 - d) < 1e-4
            });
            match groups.iter_mut().find(|g| near(g.1, at.x)) {
                Some(g) => {
                    g.2 += c;
                    g.3.push(a);
                }
                None => groups.push((a, at.x, c, vec![a])),
            }
        }
        for (_, x, sum, atoms) in groups {
            if sum > 1.0 + 1e-3 {
                let labels: Vec<String> = atoms.iter().map(|&a| st.atoms[a].label.clone()).collect();
                let mut d = Diagnostic::new(
                    Severity::Error,
                    Level::Model,
                    "occupancy",
                    format!("the site at ({:.4}, {:.4}, {:.4}) of phase {} holds {} with a total occupancy of {sum:.4}", x[0], x[1], x[2], p + 1, labels.join(", ")),
                    "a site holds at most one atom: the chemical occupancies of the atoms that share it add up to 1 at most (FullProf's Occ is the chemical occupancy times m/M, manual p. 102)",
                )
                .entry("atom_occupancy")
                .fix("lower the occupancies, or tie them so that they add up to 1 (a linear constraint)");
                for &a in &atoms {
                    d = d.addr(address(Field::AtomOcc(p, a), m));
                }
                out.push(d);
            }
        }
    }
}

/// Values a site's symmetry (or the crystal system) fixes, set free by the refinement: a parameter that moves only
/// such values (a cubic cell's a, b, c under one parameter is the symmetry's own tie and is not reported).
fn symmetry_level(m: &Model, coded: &[(Field, String)], out: &mut Vec<Diagnostic>) {
    let mut params: Vec<&String> = coded.iter().map(|c| &c.1).collect();
    params.dedup();
    for parameter in params {
        let fields: Vec<Field> = coded.iter().filter(|c| &c.1 == parameter).map(|c| c.0).collect();
        if fields.is_empty() || !fields.iter().all(|f| m.fixed_by(*f).is_some()) {
            continue;
        }
        for f in fields {
            let by = m.fixed_by(f).unwrap_or("symmetry");
            out.push(
                Diagnostic::new(Severity::Warning, Level::Model, "symmetry-fixed", format!("{} is fixed by the {by}, but the refinement sets it free ({parameter})", address(f, m)), format!("the {by} fixes this value or ties it to another: it cannot move on its own, and the engine keeps it where the symmetry puts it"))
                    .addr(address(f, m))
                    .entry(dictionary_id(f, m))
                    .fix("take it out of the refined parameters, or refine the value it follows"),
            );
        }
    }
}

// Levels (c) and (d): the refinement and the run --------------------------------------------------------------------

/// The dictionary's explanation of a group of parameters: the entry whose correlation (the entry and its `with`)
/// covers the most members of the group, at least two.
fn explain(ids: &[Vec<&'static str>]) -> Option<(&'static str, &'static Correlation)> {
    let mut best: Option<(usize, &'static str, &'static Correlation)> = None;
    for p in crate::parameters::dictionary() {
        for c in p.correlations {
            let covered = ids.iter().filter(|v| v.iter().any(|id| *id == p.id || c.with.contains(id))).count();
            if covered >= 2 && best.is_none_or(|b| covered > b.0) {
                best = Some((covered, p.id, c));
            }
        }
    }
    best.map(|b| (b.1, b.2))
}

/// The diagnostic of a null or near-null direction of the scaled normal matrix (R24): `parts` the free parameters'
/// components of the eigenvector (index into `set.free`, component), `ratio` its eigenvalue over the largest. Names the
/// parameters whose components exceed 0.2, largest first, and takes why and what to do from the dictionary's
/// correlations. The engine stops with it on a singular matrix (`lsq::Stop::Singular`, W11).
pub fn null_direction(m: &Model, set: &ParameterSet, parts: &[(usize, f64)], ratio: f64) -> Diagnostic {
    let ids: Vec<Vec<&'static str>> = set.free.iter().map(|p| {
        let mut v: Vec<&'static str> = p.fields.iter().map(|q| dictionary_id(q.0, m)).collect();
        v.dedup();
        v
    }).collect();
    let mut parts: Vec<(usize, f64)> = parts.iter().copied().filter(|x| x.1.abs() > 0.2).collect();
    parts.sort_by(|x, y| y.1.abs().total_cmp(&x.1.abs()));
    let names: Vec<String> = parts.iter().map(|&(j, c)| format!("{:+.2}·{}", c, set.free[j].name)).collect();
    let group: Vec<Vec<&'static str>> = parts.iter().map(|&(j, _)| ids[j].clone()).collect();
    let (code, sev) = if ratio < NULL { ("singular", Severity::Error) } else { ("near-singular", Severity::Warning) };
    let mut dg = Diagnostic::new(sev, Level::Run, code, format!("the data cannot fix the combination {} (eigenvalue {:.1e} of the largest)", names.join(" "), ratio), "");
    for &(j, _) in &parts {
        dg = dg.addr(set.free[j].name.clone()).entry(ids[j][0]);
    }
    match explain(&group) {
        Some((id, c)) => {
            dg.why = c.why.to_string();
            dg.suggestion = Some(c.tip.to_string());
            dg = dg.entry(id);
        }
        None => {
            dg.why = "these parameters change the pattern in nearly the same way; the dictionary has no entry for this group yet".into();
            dg.suggestion = Some("fix one of them, or tie them".into());
        }
    }
    dg
}

/// One direction of the scaled normal matrix (E5): its eigenvalue over the largest, and its components on the free
/// parameters (name, component; a unit vector in the scaled parameters).
#[derive(Clone, Debug, PartialEq)]
pub struct Direction {
    pub ratio: f64,
    pub components: Vec<(String, f64)>,
}

/// The whole eigen-spectrum of the normal matrix scaled to unit diagonal, as the engine solves it (E5): every
/// direction, largest eigenvalue first, at the model's values. `refinement_level` turns the null and near-null ones
/// (ratio below 10⁻⁶) into diagnostics; this returns them all. A parameter that moves nothing is left out (its
/// diagonal is zero).
pub fn directions(m: &Model, data: &Pattern, set: &ParameterSet) -> Vec<Direction> {
    let (a, _) = crate::lsq::normal_matrix(m, data, set);
    let live: Vec<usize> = (0..set.free.len()).filter(|&j| a[j][j] > 0.0).collect();
    let n = live.len();
    if n == 0 {
        return Vec::new();
    }
    let d: Vec<f64> = live.iter().map(|&j| a[j][j].sqrt()).collect();
    let s: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|k| a[live[i]][live[k]] / (d[i] * d[k])).collect()).collect();
    let (lam, vec) = crate::lsq::symmetric_eigen(&s);
    let top = lam.iter().copied().fold(0.0, f64::max);
    let mut out: Vec<Direction> = (0..n).map(|k| Direction { ratio: lam[k] / top, components: (0..n).map(|i| (set.free[live[i]].name.clone(), vec[i][k])).collect() }).collect();
    out.sort_by(|x, y| y.ratio.total_cmp(&x.ratio));
    out
}

fn refinement_level(m: &Model, data: &Pattern, set: &ParameterSet, opts: &Options, out: &mut Vec<Diagnostic>) {
    let np = set.free.len();
    if np == 0 {
        return;
    }
    let (a, cols) = crate::lsq::normal_matrix(m, data, set);
    let ids: Vec<Vec<&'static str>> = set.free.iter().map(|p| {
        let mut v: Vec<&'static str> = p.fields.iter().map(|q| dictionary_id(q.0, m)).collect();
        v.dedup();
        v
    }).collect();
    // A free parameter that does not change the fitted pattern.
    let mut live = Vec::new();
    for (j, p) in set.free.iter().enumerate() {
        let moves = cols[j].iter().any(|v| *v != 0.0) || a[j][j] > 0.0;
        if moves {
            live.push(j);
        } else {
            out.push(
                Diagnostic::new(Severity::Warning, Level::Refinement, "no-effect", format!("{} does not change the pattern at any point used in the fit", p.name), "its derivative is zero at every weighted point: the data say nothing about it, and the engine leaves it where it is")
                    .addr(p.name.clone())
                    .entry(ids[j][0])
                    .fix("take it out of the refined parameters (or check that the points it acts on are not excluded)"),
            );
        }
    }
    if !opts.run || live.len() < 2 {
        return;
    }
    // The run: the matrix scaled to unit diagonal, as the engine solves it.
    let n = live.len();
    let d: Vec<f64> = live.iter().map(|&j| a[j][j].sqrt()).collect();
    let s: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|k| a[live[i]][live[k]] / (d[i] * d[k])).collect()).collect();
    let (lam, vec) = crate::lsq::symmetric_eigen(&s);
    let top = lam.iter().copied().fold(0.0, f64::max);
    for k in 0..n {
        let ratio = lam[k] / top;
        if ratio >= NEAR {
            continue;
        }
        let parts: Vec<(usize, f64)> = (0..n).map(|i| (live[i], vec[i][k])).collect();
        let dg = null_direction(m, set, &parts, ratio);
        out.push(dg);
    }
    // The correlations the engine would report at these values (pseudo-inverse over the non-null directions).
    let mut inv = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            inv[i][j] = (0..n).filter(|&k| lam[k] >= NULL * top).map(|k| vec[i][k] * vec[j][k] / lam[k]).sum();
        }
    }
    let mut pairs: Vec<(usize, usize, f64)> = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            if inv[i][i] <= 0.0 || inv[j][j] <= 0.0 {
                continue;
            }
            let c = inv[i][j] / (inv[i][i] * inv[j][j]).sqrt();
            if c.abs() >= opts.correlation {
                pairs.push((live[i], live[j], c));
            }
        }
    }
    // Pairs that share parameters form groups (zero, displacement and cell; a background's terms): one diagnostic for
    // each group of more than two, one for each lone pair. Pairs within a group already reported as a null or
    // near-null direction are its consequence and are not repeated.
    let mut group_of: Vec<usize> = (0..np).collect();
    fn root(g: &mut [usize], i: usize) -> usize {
        let mut r = i;
        while g[r] != r {
            r = g[r];
        }
        g[i] = r;
        r
    }
    let in_null = |i: usize, j: usize| out.iter().any(|d| matches!(d.code, "singular" | "near-singular") && d.addresses.contains(&set.free[i].name) && d.addresses.contains(&set.free[j].name));
    let pairs: Vec<(usize, usize, f64)> = pairs.into_iter().filter(|&(i, j, _)| !in_null(i, j)).collect();
    for &(i, j, _) in &pairs {
        let (ri, rj) = (root(&mut group_of, i), root(&mut group_of, j));
        group_of[ri] = rj;
    }
    type Pair = (usize, usize, f64);
    let mut groups: Vec<(usize, Vec<Pair>)> = Vec::new();
    for &(i, j, c) in &pairs {
        let r = root(&mut group_of, i);
        match groups.iter_mut().find(|g| g.0 == r) {
            Some(g) => g.1.push((i, j, c)),
            None => groups.push((r, vec![(i, j, c)])),
        }
    }
    for (_, g) in groups {
        let mut members: Vec<usize> = g.iter().flat_map(|p| [p.0, p.1]).collect();
        members.sort_unstable();
        members.dedup();
        let what = if g.len() == 1 {
            let (i, j, c) = g[0];
            format!("{} and {} will be {:.0} % correlated", set.free[i].name, set.free[j].name, 100.0 * c)
        } else {
            let worst = g.iter().map(|p| p.2.abs()).fold(0.0, f64::max);
            let list: Vec<String> = g.iter().take(6).map(|&(i, j, c)| format!("{}–{} {:.0} %", set.free[i].name, set.free[j].name, 100.0 * c)).collect();
            format!("{} parameters will be correlated with one another: {} pairs at {:.0} % or more, up to {:.1} % ({}{})", members.len(), g.len(), 100.0 * opts.correlation, 100.0 * worst, list.join(", "), if g.len() > 6 { ", …" } else { "" })
        };
        let mut dg = Diagnostic::new(Severity::Hint, Level::Run, "correlation", what, "");
        for &j in &members {
            dg = dg.addr(set.free[j].name.clone()).entry(ids[j][0]);
        }
        let group: Vec<Vec<&'static str>> = members.iter().map(|&j| ids[j].clone()).collect();
        match explain(&group) {
            Some((id, cr)) => {
                dg.why = cr.why.to_string();
                dg.suggestion = Some(cr.tip.to_string());
                dg.severity = Severity::Warning;
                dg = dg.entry(id);
            }
            None => dg.why = "they change the pattern in similar ways at these values; the dictionary has no entry for them yet".into(),
        }
        out.push(dg);
    }
}

/// The predicted correlation matrix of the free parameters at the model's values (what the run level reports), for
/// the test against the engine's after a refinement.
pub fn predicted_correlations(m: &Model, data: &Pattern, set: &ParameterSet) -> Vec<Vec<f64>> {
    let np = set.free.len();
    let (a, _) = crate::lsq::normal_matrix(m, data, set);
    let live: Vec<usize> = (0..np).filter(|&j| a[j][j] > 0.0).collect();
    let n = live.len();
    let d: Vec<f64> = live.iter().map(|&j| a[j][j].sqrt()).collect();
    let s: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|k| a[live[i]][live[k]] / (d[i] * d[k])).collect()).collect();
    let (lam, vec) = crate::lsq::symmetric_eigen(&s);
    let top = lam.iter().copied().fold(0.0, f64::max);
    let mut out = vec![vec![0.0; np]; np];
    let mut inv = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            inv[i][j] = (0..n).filter(|&k| lam[k] >= NULL * top).map(|k| vec[i][k] * vec[j][k] / lam[k]).sum();
        }
    }
    for i in 0..n {
        for j in 0..n {
            if inv[i][i] > 0.0 && inv[j][j] > 0.0 {
                out[live[i]][live[j]] = inv[i][j] / (inv[i][i] * inv[j][j]).sqrt();
            }
        }
    }
    out
}

/// Values outside the engine's limits (P39).
fn ranges(m: &Model, doc: Option<&syntax::Document>, out: &mut Vec<Diagnostic>) {
    // The wavelength has its own diagnostic (λ must be positive), not this one as well.
    for f in m.all_fields().into_iter().filter(|&f| f != Field::Wavelength) {
        let v = m.get(f);
        let (lo, hi) = limits(f, m);
        let below = lo.is_some_and(|l| v < l);
        let above = hi.is_some_and(|h| v > h);
        if below || above {
            let id = dictionary_id(f, m);
            let bound = if below { format!("below its lower limit {}", lo.unwrap()) } else { format!("above its upper limit {}", hi.unwrap()) };
            out.push(
                Diagnostic::new(Severity::Warning, Level::File, "range", format!("{} = {v} is {bound}; the engine will only move it towards the limit", address(f, m)), crate::parameters::lookup(id).map_or(String::new(), |p| format!("{} (the dictionary's limits: {})", p.meaning, p.limits)))
                    .at(doc.and_then(|d| locate(d, f, m)))
                    .addr(address(f, m))
                    .entry(id)
                    .fix(match (lo, hi) {
                        (Some(l), Some(h)) => format!("give a value from {l} to {h}"),
                        (Some(l), None) => format!("give a value of at least {l}"),
                        (None, Some(h)) => format!("give a value of at most {h}"),
                        _ => String::new(),
                    }),
            );
        }
    }
}

/// Check a control file's text; its data file and the operations' files are relative to `dir`.
pub fn check_control(text: &str, dir: &Path, opts: &Options) -> Vec<Diagnostic> {
    let doc = syntax::parse(text);
    let (control, read) = schema::read_with(&doc, Some(dir));
    let tof = control.as_ref().is_some_and(|c| c.model.tof.is_some());
    let mut out = file_level(&doc, &read, tof);
    let Some(c) = control else { return out };
    ranges(&c.model, Some(&doc), &mut out);
    // The data.
    let span_of = |key: &str| doc.tables.iter().find(|t| t.path == ["pattern"]).and_then(|t| t.get(key)).map(|e| e.value.span);
    let mut data = None;
    if let Some(src) = &c.data {
        let path = dir.join(&src.file);
        if !path.is_file() {
            out.push(Diagnostic::new(Severity::Error, Level::Model, "data-missing", format!("the data file {} is not found", src.file), format!("it is read relative to the control file's folder, {}", dir.display())).at(span_of("data")).addr("pattern1.data").fix("give its path relative to the control file"));
        } else {
            let load = |f: &str, r: Option<&str>| -> Result<crate::formats::Scan, String> {
                crate::formats::read_path(&dir.join(f), r).map_err(|e| e.to_string())?.scans.into_iter().next().ok_or(format!("{f}: no scan"))
            };
            // A TOF bank's data on its time axis, as the calculation reads them (`control::data_scan`).
            match crate::control::data_scan_with(dir, &src.file, src.reader.as_deref(), &src.reader_options, tof).and_then(|s| crate::data_ops::apply_all(&s, &src.operations, &load)) {
                Ok(s) => data = Some(s.pattern),
                Err(e) => out.push(Diagnostic::new(Severity::Error, Level::Model, "data-missing", format!("the data file {} does not read: {e}", src.file), "the reader named (or the one the registry recognises) cannot read it, or an operation fails").at(span_of("data")).addr("pattern1.data").fix("check the reader with drs data identify")),
            }
        }
    }
    model_level(&c.model, data.as_ref(), &mut out);
    if let (Some(d), Some([a, _, b])) = (data.as_ref(), c.fullprof.range) {
        range_against(d, a, b, if tof { " µs" } else { "°" }, &mut out);
    }
    // The refinement.
    let coded: Vec<(Field, String)> = c.refinement.parameters.iter().flat_map(|terms| {
        let name = terms.first().map(|t| address(t.0, &c.model)).unwrap_or_default();
        terms.iter().map(move |t| (t.0, name.clone()))
    }).collect();
    symmetry_level(&c.model, &coded, &mut out);
    // A value that names two parameters (the first value of each): the two are one parameter written twice.
    let span_params = doc.tables.iter().find(|t| t.path == ["refine"]).and_then(|t| t.get("parameters")).map(|e| e.key_span);
    let firsts: Vec<Field> = c.refinement.parameters.iter().filter_map(|t| t.first().map(|x| x.0)).collect();
    for (k, f) in firsts.iter().enumerate() {
        if firsts[..k].contains(f) {
            out.push(Diagnostic::new(Severity::Error, Level::Refinement, "refinement", format!("{} names two parameters", address(*f, &c.model)), "a parameter is named by its first value; two parameters on the same value move it twice over and the data cannot tell them apart (the normal matrix is singular)").at(span_params).addr(address(*f, &c.model)).entry(dictionary_id(*f, &c.model)).fix("keep one of the two, or tie the values in one parameter"));
        }
    }
    if out.iter().any(|x| x.code == "refinement") {
        return out;
    }
    let set = match c.refinement.build(&c.model) {
        Ok(s) => s,
        Err(e) => {
            let span = doc.tables.iter().find(|t| t.path == ["refine"]).and_then(|t| t.get("parameters")).map(|e| e.key_span);
            out.push(Diagnostic::new(Severity::Error, Level::Refinement, "refinement", e, "a value can be moved by one parameter only, and a linear restraint names free parameters").at(span).fix("give each value to one parameter; name in restraints only refined parameters"));
            return out;
        }
    };
    if let Some(d) = &data {
        if !out.iter().any(|x| x.severity == Severity::Error) {
            refinement_level(&c.model, d, &set, opts, &mut out);
        }
    }
    out
}

/// `unit`: "°" for 2θ, " µs" for a TOF bank (TOF-min, TOF-max).
fn range_against(d: &Pattern, thmin: f64, thmax: f64, unit: &str, out: &mut Vec<Diagnostic>) {
    if d.len() < 2 {
        return;
    }
    let (lo, hi) = (d.x[0], d.x[d.len() - 1]);
    let step = (hi - lo) / (d.len() - 1) as f64;
    // A TOF bank's TOF-max is the period of its Fourier or Chebyshev background (`tof::Bank::tmax`, P05), not the
    // data's end (CeO2_PEARL.pcr: 19380.78 µs, the data to 19400.19): its TOF-min alone is compared.
    let tof = unit.contains("µs");
    if (thmin - lo).abs() > step.abs() || (!tof && (thmax - hi).abs() > step.abs()) {
        out.push(Diagnostic::new(Severity::Warning, Level::Model, "range-data", format!("the range {thmin}–{thmax}{unit} is not the data's, {lo}–{hi}{unit}"), "FullProf reads the points of the data file; a range that does not match them is a file of another pattern or a mistyped Thmin, Thmax").addr("pattern1.range").fix(format!("write Thmin {lo} and Thmax {hi}, or check the data file")));
    }
}

/// Check a `.pcr` in `dir` (its data found as FullProf finds it, its resolution and reflection files as `drs calc`
/// finds them).
pub fn check_pcr(text: &str, dir: &Path, stem: &str, opts: &Options) -> Vec<Diagnostic> {
    check_pcr_with_data(text, dir, stem, None, opts)
}

/// Check a `.pcr` with its data as FullProf takes them (L079): the file given (FullProf's command line, `fp2k CODFIL
/// DATFIL`), else CODFIL.dat next to the `.pcr`. Without either, a file the core finds by the `.pcr`'s comment or in
/// the folders around it is used, and a hint says that FullProf would not find it so.
pub fn check_pcr_with_data(text: &str, dir: &Path, stem: &str, data: Option<&Path>, opts: &Options) -> Vec<Diagnostic> {
    let mut out = Vec::new();
    let mut pcr = match crate::pcr::read(text) {
        Ok(p) => p,
        Err(e) => {
            let span = Span { start: 0, end: 0, line: e.line, column: 1 };
            out.push(Diagnostic::new(Severity::Error, Level::File, "syntax", e.message.clone(), "FullProf's .pcr is read line by line in a fixed order (manual pp. 80–126): a line out of place stops the reading").at(Some(span)));
            return out;
        }
    };
    pcr.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
    if let Some(f) = &pcr.resolution_file {
        if crate::powder::find_file(dir, f).is_none() {
            out.push(Diagnostic::new(Severity::Error, Level::Model, "irf-missing", format!("the resolution file {f} is not found"), "Res ≠ 0 names a file of the instrument's widths, looked for next to the .pcr and in the folders around it").addr("pattern1.resolution").entry("instrument_resolution").fix("put the file next to the .pcr, or set Res = 0"));
            return out;
        }
    }
    let own = dir.join(format!("{stem}.dat"));
    let given = match data {
        Some(d) => Some(d.to_path_buf()),
        None => own.exists().then_some(own),
    };
    let l = match crate::powder::Model::from_pcr_with_data(&pcr, dir, stem, given.as_deref()) {
        Ok(l) => l,
        Err(e) => {
            out.push(Diagnostic::new(Severity::Error, Level::Model, "model", e, "the core builds no model from this file"));
            return out;
        }
    };
    let m = &l.model;
    ranges(m, None, &mut out);
    // A .bac in pairs: FullProf 8.20 reads every .bac as its first format (manual p. 8, footnote; L082).
    if pcr.flag("Nba") == 1 {
        if let Some(bac) = crate::powder::find_file(dir, &format!("{stem}.bac")) {
            let first = std::fs::read(&bac).ok().map(|b| b.iter().map(|&c| c as char).collect::<String>()).and_then(|t| t.lines().find(|l| !l.trim().is_empty()).map(str::to_string)).unwrap_or_default();
            if first.trim_start().starts_with('!') {
                out.push(Diagnostic::new(Severity::Warning, Level::Model, "bac-format", format!("{stem}.bac is written as a comment and pairs (the manual's second format)"), "FullProf reads every .bac as the first format, start step end and then the intensities, and takes the pairs' numbers in order as intensities (the manual's footnote, p. 8: the second format \"is not available at present\"; L082)").addr("pattern1.background").fix("write it in the first format at the data's points: drs background convert <points> --to <stem>.bac --data <data>"));
            }
        }
    }
    match (&given, &l.data_path) {
        (None, Some(found)) => {
            let name = found.file_name().map(|n| n.to_string_lossy().to_string()).unwrap_or_default();
            out.push(Diagnostic::new(Severity::Hint, Level::Model, "data-elsewhere", format!("the data are read from {}, found by the .pcr's comment or in the folders around it", found.display()), format!("FullProf takes the data file from its command line, else {stem}.dat next to the .pcr: run without a file it would not find this one")).addr("pattern1.data").fix(format!("run FullProf as fp2k {stem} {name}, and check as drs check {stem}.pcr {name}")));
        }
        (_, None) => {
            out.push(Diagnostic::new(Severity::Error, Level::Model, "data-missing", format!("no data file for {stem}.pcr is found"), format!("FullProf reads the data file given on its command line (fp2k {stem} <data>), else {stem}.dat next to the .pcr")).addr("pattern1.data").fix(format!("give the data: drs check {stem}.pcr <data>, or put them next to the .pcr as {stem}.dat")));
        }
        _ => {}
    }
    model_level(m, l.data.as_ref(), &mut out);
    if let Some(d) = &l.data {
        if pcr.step > 0.0 {
            range_against(d, pcr.thmin, pcr.thmax, if pcr.flag("Uni") == 1 { " µs" } else { "°" }, &mut out);
        }
    }
    // The parameter numbers: 1…N with no hole, N as declared.
    match crate::refinable::pcr_parameter_numbers(&pcr, m) {
        Ok(nums) => {
            let coded: Vec<(Field, String)> = nums.iter().flat_map(|(n, fs)| fs.iter().map(move |f| (*f, format!("parameter {n}")))).collect();
            symmetry_level(m, &coded, &mut out);
            let max = nums.keys().copied().max().unwrap_or(0) as usize;
            let holes: Vec<u32> = (1..=max as u32).filter(|k| !nums.contains_key(k)).collect();
            // With Aut = 1 FullProf numbers the codewords itself: no hole is left (manual pp. 68–69).
            if !holes.is_empty() && pcr.flag("Aut") == 0 {
                out.push(Diagnostic::new(Severity::Error, Level::Refinement, "parameter-hole", format!("no codeword uses parameter number{} {}", if holes.len() > 1 { "s" } else { "" }, holes.iter().map(|h| h.to_string()).collect::<Vec<_>>().join(", ")), "FullProf numbers the refined parameters 1…N through the codewords (C = sign(a)(10p + |a|), manual p. 13); a number no codeword uses leaves an empty row and column in its matrix, and FullProf stops with \"Hole in the matrix!, missing element\" (L071)").addr("refinement").fix(format!("renumber the codewords 1…{} without gaps, or set Aut = 1 to let FullProf number them", nums.len())));
            }
            let declared = pcr.declared_parameters;
            if declared != max && pcr.flag("Aut") == 0 {
                out.push(Diagnostic::new(Severity::Error, Level::Refinement, "parameter-count", format!("the file declares {declared} refined parameters and its codewords number up to {max}"), "line 13's count (MAXS) is the size of FullProf's matrix: fewer than the highest number is out of range, more leaves empty rows (\"Hole in the matrix\", manual pp. 68–69)").addr("refinement").fix(format!("write {max} on the line \"Number of refined parameters\"")).at(None));
            }
        }
        Err(e) => out.push(Diagnostic::new(Severity::Warning, Level::Refinement, "refinement", e, "the core does not refine some of what the file codes; the run levels are not checked")),
    }
    if out.iter().any(|x| x.severity == Severity::Error) {
        return out;
    }
    match ParameterSet::from_pcr(&pcr, m) {
        Ok(set) => {
            if let Some(d) = &l.data {
                refinement_level(m, d, &set, opts, &mut out);
            }
        }
        Err(e) => out.push(Diagnostic::new(Severity::Warning, Level::Refinement, "refinement", e, "the refinement is not built: the run level is not checked")),
    }
    out
}

/// Is any diagnostic an error?
pub fn has_errors(d: &[Diagnostic]) -> bool {
    d.iter().any(|x| x.severity == Severity::Error)
}
