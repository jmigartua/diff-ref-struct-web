//! Per-format reading options (the app's request 51, app/DESIGN.md §11.2): what a reader can be told that its file
//! does not say, or says ambiguously. Each format declares its options (`Format::options`), with a kind, a default
//! and a line of help, so that the app's Settings ▸ Data Formats and Data ▸ Read As… are generated from the
//! declarations; `read_with` reads with them, and a control file keeps them with its data (`[pattern]
//! reader_options`), so that it reads the same whatever the Settings say.
//!
//! The options:
//! - `scan`: which scan of a file of several is read (a GSAS file's bank, a Bruker file's range, an XRDML or RAS
//!   file's scan, a CIF's block, a table's series), 1 the first. Without it a control file reads the first.
//! - `wavelength` (Å): for a file that states none; a wavelength the file states is kept.
//! - `axis`: the quantity on x of a column file (2θ, TOF, Q, d, energy); empty: from the column's name, else 2θ.
//! - `tof_unit`: µs or ms for a time-of-flight x (ms read as 1000 µs).
//! - the general table's layout (`table::TableOptions`): `delimiter`, `decimal`, `header_lines`, and the columns
//!   `x_column`, `y_column`, `sigma_column`, counted from 1 (`sigma_column` 0: no σ column, σ = √y). Empty: as the
//!   reader recognises them.

use super::table::TableOptions;
use super::{Axis, Dataset, Format};
use crate::data::DataError;
use std::collections::BTreeMap;

/// The kind of an option's value.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
pub enum OptionKind {
    /// An integer, at least `min`.
    Integer { min: i64 },
    /// One of these words.
    Choice(&'static [&'static str]),
    /// A number in this unit.
    Number { unit: &'static str },
    /// A column of the file, counted from 1 (0 where `none` is allowed: no such column).
    Column { none: bool },
}

/// An option a reader takes.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize))]
pub struct FormatOption {
    pub id: &'static str,
    pub name: &'static str,
    pub kind: OptionKind,
    /// The default as text; empty: the reader decides from the file.
    pub default: &'static str,
    pub help: &'static str,
}

/// An option's value.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(untagged))]
pub enum OptionValue {
    Integer(i64),
    Number(f64),
    Text(String),
}

impl OptionValue {
    pub fn as_f64(&self) -> Option<f64> {
        match self {
            OptionValue::Integer(v) => Some(*v as f64),
            OptionValue::Number(v) => Some(*v),
            OptionValue::Text(_) => None,
        }
    }

    pub fn as_text(&self) -> String {
        match self {
            OptionValue::Integer(v) => v.to_string(),
            OptionValue::Number(v) => v.to_string(),
            OptionValue::Text(s) => s.clone(),
        }
    }
}

/// Options by id.
pub type Options = BTreeMap<String, OptionValue>;

const fn scan(name: &'static str, help: &'static str) -> FormatOption {
    FormatOption { id: "scan", name, kind: OptionKind::Integer { min: 1 }, default: "1", help }
}

pub const BANK: FormatOption = scan("Bank", "The bank of the file read, 1 the first.");
pub const RANGE: FormatOption = scan("Range", "The range (scan) of the file read, 1 the first.");
pub const SCAN: FormatOption = scan("Scan", "The scan of the file read, 1 the first.");
pub const BLOCK: FormatOption = scan("Block", "The data block of the file read, 1 the first.");
pub const SERIES: FormatOption = scan("Series", "The series (a y column, or a member of a long table) read, 1 the first.");
pub const WAVELENGTH: FormatOption = FormatOption {
    id: "wavelength",
    name: "Wavelength",
    kind: OptionKind::Number { unit: "Å" },
    default: "",
    help: "The wavelength, for a file that states none; a wavelength the file states is kept.",
};
pub const AXIS: FormatOption = FormatOption {
    id: "axis",
    name: "Quantity on x",
    kind: OptionKind::Choice(&["2theta", "tof", "q", "d", "energy"]),
    default: "",
    help: "2θ (°), time of flight (µs), Q (Å⁻¹), d (Å) or energy (keV); empty: from the x column's name, else 2θ.",
};
pub const TOF_UNIT: FormatOption = FormatOption { id: "tof_unit", name: "Time-of-flight unit", kind: OptionKind::Choice(&["µs", "ms"]), default: "µs", help: "The unit of a time-of-flight x: ms is read as 1000 µs." };
pub const DELIMITER: FormatOption = FormatOption { id: "delimiter", name: "Delimiter", kind: OptionKind::Choice(&[",", ";", "tab", "blanks"]), default: "", help: "What separates the columns; empty: the one that splits the most rows alike." };
pub const DECIMAL: FormatOption = FormatOption { id: "decimal", name: "Decimal mark", kind: OptionKind::Choice(&[".", ","]), default: "", help: "The decimal mark; empty: ',' when numbers are written so and ',' does not separate the columns." };
pub const HEADER_LINES: FormatOption = FormatOption { id: "header_lines", name: "Header lines", kind: OptionKind::Integer { min: 0 }, default: "", help: "The lines before the first row of numbers; empty: up to the first row of numbers." };
pub const X_COLUMN: FormatOption = FormatOption { id: "x_column", name: "x column", kind: OptionKind::Column { none: false }, default: "", help: "The column of x, from 1; empty: by its name, else the first." };
pub const Y_COLUMN: FormatOption = FormatOption { id: "y_column", name: "y column", kind: OptionKind::Column { none: false }, default: "", help: "The column of the counts, from 1; empty: by its name, else the second." };
pub const SIGMA_COLUMN: FormatOption = FormatOption { id: "sigma_column", name: "σ column", kind: OptionKind::Column { none: true }, default: "", help: "The column of σ, from 1; 0: none, σ = √y; empty: by its name, else the third if there is one." };

/// The general table's options.
pub const TABLE: &[FormatOption] = &[SERIES, AXIS, DELIMITER, DECIMAL, HEADER_LINES, X_COLUMN, Y_COLUMN, SIGMA_COLUMN, WAVELENGTH, TOF_UNIT];

/// Check options against a format's declarations: each known, of its kind and range.
pub fn validate(f: &Format, o: &Options) -> Result<(), String> {
    for (id, v) in o {
        let Some(d) = f.options.iter().find(|d| d.id == id) else {
            let known: Vec<&str> = f.options.iter().map(|d| d.id).collect();
            return Err(format!("{}: no option \"{id}\" ({})", f.id, if known.is_empty() { "it takes none".to_string() } else { format!("it takes {}", known.join(", ")) }));
        };
        let ok = match (d.kind, v) {
            (OptionKind::Integer { min }, OptionValue::Integer(n)) => *n >= min,
            (OptionKind::Column { none }, OptionValue::Integer(n)) => *n >= 1 || (none && *n == 0),
            (OptionKind::Number { .. }, x) => x.as_f64().is_some_and(|x| x.is_finite() && x > 0.0),
            (OptionKind::Choice(words), OptionValue::Text(s)) => words.contains(&s.as_str()),
            _ => false,
        };
        if !ok {
            let want = match d.kind {
                OptionKind::Integer { min } => format!("an integer from {min}"),
                OptionKind::Column { none } => format!("a column from 1{}", if none { ", or 0 for none" } else { "" }),
                OptionKind::Number { unit } => format!("a positive number in {unit}"),
                OptionKind::Choice(words) => format!("one of {}", words.iter().map(|w| format!("\"{w}\"")).collect::<Vec<_>>().join(", ")),
            };
            return Err(format!("{}: option \"{id}\" = {} must be {want}", f.id, v.as_text()));
        }
    }
    Ok(())
}

fn axis_of(word: &str) -> Option<Axis> {
    Some(match word {
        "2theta" => Axis::TwoTheta,
        "tof" => Axis::TimeOfFlight,
        "q" => Axis::Q,
        "d" => Axis::D,
        "energy" => Axis::Energy,
        _ => return None,
    })
}

/// The general table's layout from the options (validated).
pub fn table_options(o: &Options) -> TableOptions {
    let int = |k: &str| match o.get(k) {
        Some(OptionValue::Integer(n)) => Some(*n),
        _ => None,
    };
    let text = |k: &str| match o.get(k) {
        Some(OptionValue::Text(s)) => Some(s.as_str()),
        _ => None,
    };
    TableOptions {
        delimiter: text("delimiter").map(|d| match d {
            "tab" => '\t',
            "blanks" => ' ',
            d => d.chars().next().unwrap_or(' '),
        }),
        decimal: text("decimal").and_then(|d| d.chars().next()),
        header_lines: int("header_lines").map(|n| n as usize),
        x: int("x_column").map(|n| n as usize - 1),
        y: int("y_column").map(|n| n as usize - 1),
        sigma: int("sigma_column").map(|n| (n > 0).then(|| n as usize - 1)),
        axis: text("axis").and_then(axis_of),
    }
}

/// Is any of the table's layout given?
pub fn has_table_layout(o: &Options) -> bool {
    ["delimiter", "decimal", "header_lines", "x_column", "y_column", "sigma_column", "axis"].iter().any(|k| o.contains_key(*k))
}

/// The options that apply after the reading: the scan chosen, the quantity on x, the time-of-flight unit and the
/// wavelength where the file states none.
pub fn apply(d: &mut Dataset, o: &Options) -> Result<(), DataError> {
    if let Some(OptionValue::Integer(n)) = o.get("scan") {
        let n = *n as usize;
        if n > d.scans.len() {
            return Err(DataError(format!("scan {n} asked for: the file holds {}", d.scans.len())));
        }
        let s = d.scans.swap_remove(n - 1);
        d.scans = vec![s];
    }
    if let Some(a) = o.get("axis").and_then(|v| axis_of(&v.as_text())) {
        for s in &mut d.scans {
            s.axis = a;
        }
    }
    if o.get("tof_unit").map(OptionValue::as_text).as_deref() == Some("ms") {
        for s in d.scans.iter_mut().filter(|s| s.axis == Axis::TimeOfFlight) {
            for x in &mut s.pattern.x {
                *x *= 1000.0;
            }
        }
    }
    if let Some(l) = o.get("wavelength").and_then(OptionValue::as_f64) {
        for s in d.scans.iter_mut().filter(|s| s.wavelength.is_none()) {
            s.wavelength = Some(super::Wavelength::single(l));
        }
    }
    Ok(())
}
