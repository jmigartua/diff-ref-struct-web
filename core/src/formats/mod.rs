//! The reader registry (R21; L061): one list of the powder data formats that FullProf and GSAS-II read, used by `drs`
//! and the app.
//!
//! Each format has a readable id and name (`bruker_raw_v4`, "Bruker RAW (version 4)"), the extensions it is usually
//! found with, where its layout comes from (FullProf's manual, pp. 133–137, or the GSAS-II importer that reads it), and
//! two functions: one gives a confidence that a file is in the format, from its content; the other reads it.
//! FullProf's numeric `Ins` codes appear only as a translation (`by_ins`, `Format::ins`), for `.pcr` import and export.
//!
//! Recognition is by content first, then by extension. A reader's confidence is one of:
//! - `SIGNATURE` (1.0): something only this format has (magic bytes, a mandatory keyword or structure) and the file
//!   reads;
//! - `STRUCTURE` (0.8): the file reads with the format's whole structure, its counts consistent (the points that
//!   Ti, step and Tf announce are all there), but other formats share that structure;
//! - `LOOSE` (0.5): the file reads, but only as columns of numbers that many formats are;
//! - 0: it does not read.
//!
//! An extension the format is usually found with adds `EXTENSION` (0.1). A file is identified when the best score
//! leads the next by at least `MARGIN` (0.1); otherwise it is ambiguous and the candidates are returned, best first.
//! Columns of numbers in a `.dat` file, which FullProf's Ins = 10 (after six comment lines), GSAS-II's xye rules and
//! the general table all read loosely, are such a case: the `.pcr` file's `Ins` or the user decides.
//!
//! What a reader returns is a `Dataset`: one `Scan` per bank, range or scan of the file, each with its pattern
//! (x, y, σ), the quantity on x, the wavelengths the file gives, and what else the file says about the measurement
//! (`Info`: title, instrument, date, sample, temperature, monitor, counting time, step). A point that the reading
//! program gives no weight (GSAS-II's w = 0, for y ≤ 0) has σ = ∞.

pub mod bruker;
pub mod cif_powder;
pub mod file_table;
pub mod files;
pub mod fullprof;
pub mod fullprof_outputs;
pub mod gsas;
pub mod hdf;
pub mod options;
pub mod panalytical;
pub mod project;
pub mod rigaku;
pub mod table;
pub mod write;
pub mod xml;
pub mod zip;

use crate::data::{DataError, Pattern};
use std::cell::OnceCell;

pub const SIGNATURE: f64 = 1.0;
pub const STRUCTURE: f64 = 0.8;
pub const LOOSE: f64 = 0.5;
pub const EXTENSION: f64 = 0.1;
pub const MARGIN: f64 = 0.1;

/// The quantity on the x axis of a scan, in the unit given.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Axis {
    /// 2θ in degrees (constant wavelength).
    #[default]
    TwoTheta,
    /// Time of flight in µs.
    TimeOfFlight,
    /// Q = 4π sin θ / λ in Å⁻¹.
    Q,
    /// d-spacing in Å.
    D,
    /// Energy in keV (energy-dispersive).
    Energy,
}

impl Axis {
    pub fn label(self) -> &'static str {
        match self {
            Axis::TwoTheta => "2θ (°)",
            Axis::TimeOfFlight => "TOF (µs)",
            Axis::Q => "Q (Å⁻¹)",
            Axis::D => "d (Å)",
            Axis::Energy => "E (keV)",
        }
    }
}

/// The wavelengths a file gives (Å): λ1, and λ2 with I(λ2)/I(λ1) for a doublet.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Wavelength {
    pub lambda1: f64,
    pub lambda2: Option<f64>,
    pub ratio: Option<f64>,
}

impl Wavelength {
    pub fn single(l: f64) -> Wavelength {
        Wavelength { lambda1: l, lambda2: None, ratio: None }
    }
}

/// What a file says about its measurement, as far as the reader finds it. Numbers are in the file's units:
/// temperature usually in K (BRML's °C are converted, as GSAS-II does), counting time in s, step in the unit of x.
#[derive(Clone, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Info {
    pub title: Option<String>,
    pub instrument: Option<String>,
    pub date: Option<String>,
    pub sample: Option<String>,
    pub temperature: Option<f64>,
    pub monitor: Option<f64>,
    pub counting_time: Option<f64>,
    /// The step the file states (Ti, step, Tf; a step-size keyword), not one computed from x.
    pub step: Option<f64>,
    /// Anything else, as key and text.
    pub other: Vec<(String, String)>,
}

/// One measured pattern of a file: a bank, a range or a scan.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Scan {
    pub label: String,
    pub axis: Axis,
    pub pattern: Pattern,
    pub wavelength: Option<Wavelength>,
    pub info: Info,
}

impl Scan {
    pub fn new(label: impl Into<String>, axis: Axis, pattern: Pattern) -> Scan {
        Scan { label: label.into(), axis, pattern, wavelength: None, info: Info::default() }
    }
}

/// What a reader returns: the scans of a file, and the format it was read as.
#[derive(Clone, Debug, PartialEq)]
pub struct Dataset {
    pub format: &'static str,
    pub scans: Vec<Scan>,
}

/// A file offered to the readers: its name (for the extension) and its bytes, with the text and zip directory
/// decoded once, when first asked for.
pub struct Probe<'a> {
    pub name: &'a str,
    pub bytes: &'a [u8],
    text: OnceCell<String>,
    archive: OnceCell<Option<Vec<String>>>,
}

impl<'a> Probe<'a> {
    pub fn new(name: &'a str, bytes: &'a [u8]) -> Probe<'a> {
        Probe { name, bytes, text: OnceCell::new(), archive: OnceCell::new() }
    }

    /// The bytes as text, each byte one character (ISO 8859-1), as FullProf and GSAS-II open their text files;
    /// UTF-8 when the bytes are valid UTF-8.
    pub fn text(&self) -> &str {
        self.text.get_or_init(|| match std::str::from_utf8(self.bytes) {
            Ok(s) => s.to_string(),
            Err(_) => self.bytes.iter().map(|&b| b as char).collect(),
        })
    }

    /// The names of the entries when the file is a zip archive.
    pub fn zip_names(&self) -> Option<&[String]> {
        self.archive.get_or_init(|| zip::Archive::open(self.bytes).ok().map(|a| a.names().map(String::from).collect())).as_deref()
    }

    /// The file's extension, lower case, without the dot ("" if none). For "x.zarr.zip" it is "zarr.zip".
    pub fn extension(&self) -> String {
        let base = self.name.rsplit(['/', '\\']).next().unwrap_or(self.name).to_lowercase();
        if base.ends_with(".zarr.zip") {
            return "zarr.zip".into();
        }
        base.rsplit_once('.').map_or(String::new(), |(_, e)| e.to_string())
    }

    /// Does the content look like text (no NUL byte in the first 4 KiB)?
    pub fn is_text(&self) -> bool {
        !self.bytes.iter().take(4096).any(|&b| b == 0)
    }
}

/// A data format of the registry.
pub struct Format {
    pub id: &'static str,
    pub name: &'static str,
    /// The extensions it is usually found with, lower case, without the dot.
    pub extensions: &'static [&'static str],
    /// FullProf's `Ins` for this format, if FullProf reads it.
    pub ins: Option<i32>,
    /// Where the layout and the reading rules come from.
    pub source: &'static str,
    /// None if the format is read; otherwise why it is only recognised.
    pub not_read: Option<&'static str>,
    pub recognise: fn(&Probe) -> f64,
    pub read: fn(&Probe) -> Result<Dataset, DataError>,
    /// The options its reader takes (the app's request 51, `options`).
    pub options: &'static [options::FormatOption],
}

impl Format {
    pub fn reads(&self) -> bool {
        self.not_read.is_none()
    }
}

/// A format a file may be in, with its confidence.
#[derive(Clone, Debug, PartialEq)]
pub struct Candidate {
    pub id: &'static str,
    pub name: &'static str,
    /// From the content alone (0, `LOOSE`, `STRUCTURE` or `SIGNATURE`).
    pub content: f64,
    /// Whether the file's extension is one the format is usually found with.
    pub extension: bool,
    pub score: f64,
}

/// What `identify` finds.
#[derive(Clone, Debug, PartialEq)]
pub struct Identification {
    /// Every format whose reader reads the file, best first.
    pub candidates: Vec<Candidate>,
    /// The reader an order of preference chose among those tied with the best (`identify_with`); None when one led
    /// by `MARGIN`, or no preference named one of the tied.
    pub preferred: Option<&'static str>,
}

impl Identification {
    /// The candidates tied with the best (within `MARGIN`) when none leads; empty when one does.
    fn tied(&self) -> Vec<&Candidate> {
        let Some(first) = self.candidates.first() else { return Vec::new() };
        match self.candidates.get(1) {
            Some(second) if first.score - second.score < MARGIN - 1e-12 => self.candidates.iter().filter(|c| first.score - c.score < MARGIN - 1e-12).collect(),
            _ => Vec::new(),
        }
    }

    /// The format: the one that leads the others by `MARGIN` or more, or else the one the preference chose.
    pub fn format(&self) -> Option<&Candidate> {
        let first = self.candidates.first()?;
        if self.tied().is_empty() {
            return Some(first);
        }
        self.preferred.and_then(|id| self.candidates.iter().find(|c| c.id == id))
    }

    /// Was it the order of preference, not the content, that decided?
    pub fn by_preference(&self) -> bool {
        self.preferred.is_some()
    }

    /// The candidates tied with the best (within `MARGIN`), when the file is ambiguous and no preference decided.
    pub fn ambiguous(&self) -> Vec<&Candidate> {
        if self.format().is_some() {
            return Vec::new();
        }
        self.tied()
    }
}

macro_rules! entry {
    ($id:literal, $name:literal, [$($e:literal),*], $ins:expr, $src:literal, $rec:path, $read:path) => {
        entry!($id, $name, [$($e),*], $ins, $src, $rec, $read, &[])
    };
    ($id:literal, $name:literal, [$($e:literal),*], $ins:expr, $src:literal, $rec:path, $read:path, $opts:expr) => {
        Format { id: $id, name: $name, extensions: &[$($e),*], ins: $ins, source: $src, not_read: None, recognise: $rec, read: $read, options: $opts }
    };
}

use options::{AXIS, BANK, BLOCK, RANGE, SCAN, TOF_UNIT, WAVELENGTH};

/// The registry: FullProf's list first (by `Ins`), then the general table, then the formats GSAS-II reads that
/// FullProf does not (L061).
pub static FORMATS: &[Format] = &[
    entry!("free_steps", "Free format: Ti, step, Tf, then the intensities", ["dat"], Some(0), "FullProf manual p. 133 (Ins = 0); GSAS-II reads the same layout (G2pwd_FP.py)", fullprof::recognise_free, fullprof::read_free, &[WAVELENGTH]),
    entry!("ill_d1a_d2b_old", "ILL D1A, D2B (original Rietveld–Hewat format)", ["dat"], Some(1), "FullProf manual p. 133 (Ins = 1)", fullprof::recognise_d1a, fullprof::read_d1a),
    entry!("ill_d1b_old", "ILL D1B (old DEC-10 format)", ["dat"], Some(2), "FullProf manual p. 133 (Ins = 2)", fullprof::recognise_d1b_old, fullprof::read_d1b_old),
    entry!("ill_d1b_d20", "ILL D1B, D20", ["dat"], Some(3), "FullProf manual p. 134 (Ins = 3)", fullprof::recognise_d20, fullprof::read_d20),
    entry!("brookhaven", "Brookhaven synchrotron (intensities and σ in alternate lines)", ["dat"], Some(4), "FullProf manual p. 134 (Ins = 4)", fullprof::recognise_brookhaven, fullprof::read_brookhaven, &[WAVELENGTH]),
    Format {
        id: "dbws",
        name: "DBWS synchrotron (DBW3.2S-8711)",
        extensions: &["dat"],
        ins: Some(-4),
        source: "FullProf manual p. 134 (Ins = −4)",
        not_read: Some("the manual names the format but gives neither its layout nor its code, and no example is in manuals/"),
        recognise: never,
        read: unread,
        options: &[],
    },
    entry!("two_axis", "General two-axis instrument (title, NPTS, Ti/step/Tf, y and σ)", ["dat"], Some(5), "FullProf manual pp. 134–135 (Ins = 5)", fullprof::recognise_two_axis, fullprof::read_two_axis, &[WAVELENGTH]),
    entry!("ill_d2b", "ILL D1A, D2B (D1ASUM, D2BSUM, ADDET, MPDSUM)", ["dat"], Some(6), "FullProf manual p. 135 (Ins = 6)", fullprof::recognise_d1a_standard, fullprof::read_d1a_standard),
    entry!("ill_d4_d20l", "ILL D4, D20L", ["dat"], Some(7), "FullProf manual p. 135 (Ins = 7)", fullprof::recognise_d4, fullprof::read_d4),
    entry!("psi_dmc", "PSI DMC, HRPT", ["dat"], Some(8), "FullProf manual p. 136 (Ins = 8)", fullprof::recognise_dmc, fullprof::read_dmc),
    entry!("socabim_uxd", "Socabim UXD (Bruker DIFFRAC-AT text)", ["uxd"], Some(9), "FullProf manual p. 136 (Ins = 9) names it without a layout; read as the UXD keyword file, believed, not checked", fullprof::recognise_uxd, fullprof::read_uxd, &[RANGE, WAVELENGTH]),
    entry!("xys", "x, y, σ columns (FullProf, six header lines)", ["dat", "xys", "xy"], Some(10), "FullProf manual p. 136 (Ins = 10)", fullprof::recognise_xys, fullprof::read_xys, &[AXIS, WAVELENGTH, TOF_UNIT]),
    entry!("variable_time", "X-rays with variable counting time (time, counts)", ["dat"], Some(11), "FullProf manual p. 136 (Ins = 11)", fullprof::recognise_variable_time, fullprof::read_variable_time, &[WAVELENGTH]),
    entry!("gsas", "GSAS powder data (STD, ESD, FXY, FXYE, ALT banks)", ["gsas", "gss", "gsa", "gda", "fxye", "fxy", "raw", "dat", "xra", "cwn"], Some(12), "FullProf manual p. 136 (Ins = 12); GSAS-II G2pwd_fxye.py", gsas::recognise, gsas::read, &[BANK, WAVELENGTH]),
    entry!("isis_tof_banks", "ISIS TOF multibank (x, y, σ under BANK lines)", ["gss", "dat"], Some(14), "FullProf News, 24 September 2003 (Ins = 14)", fullprof::recognise_isis, fullprof::read_isis, &[BANK, TOF_UNIT]),
    entry!("table", "Table of columns (CSV, TSV, blanks; delimiter, decimal mark, header and columns recognised)", ["csv", "tsv", "txt", "dat", "xy"], None, "the general reader of L061; GSAS-II's G2pwd_csv.py reads a subset", table::recognise_table, table::read_table, options::TABLE),
    entry!("bruker_raw_v1", "Bruker RAW (version 1, DIFFRAC-AT)", ["raw"], None, "GSAS-II G2pwd_BrukerRAW.py recognises it and does not read it", bruker::recognise_raw_v1, unread),
    entry!("bruker_raw_v2", "Bruker RAW (version 2)", ["raw"], None, "GSAS-II G2pwd_BrukerRAW.py", bruker::recognise_raw_v2, bruker::read_raw_v2, &[RANGE]),
    entry!("bruker_raw_v3", "Bruker RAW (version 3, DIFFRAC plus)", ["raw"], None, "GSAS-II G2pwd_BrukerRAW.py", bruker::recognise_raw_v3, bruker::read_raw_v3, &[RANGE]),
    entry!("bruker_raw_v4", "Bruker RAW (version 4, DIFFRAC.SUITE)", ["raw"], None, "GSAS-II G2pwd_BrukerRAW.py", bruker::recognise_raw_v4, bruker::read_raw_v4, &[RANGE]),
    entry!("bruker_brml", "Bruker BRML (DIFFRAC.SUITE, zip of XML)", ["brml"], None, "GSAS-II G2pwd_BrukerBRML.py", bruker::recognise_brml, bruker::read_brml, &[RANGE]),
    entry!("panalytical_xrdml", "PANalytical XRDML", ["xrdml", "xml"], None, "GSAS-II G2pwd_Panalytical.py", panalytical::recognise, panalytical::read, &[SCAN]),
    entry!("rigaku_ras", "Rigaku RAS", ["ras"], None, "GSAS-II G2pwd_rigaku.py", rigaku::recognise_ras, rigaku::read_ras, &[SCAN]),
    entry!("rigaku_rasx", "Rigaku RASX (zip)", ["rasx"], None, "GSAS-II G2pwd_rigaku.py", rigaku::recognise_rasx, rigaku::read_rasx, &[SCAN]),
    entry!("rigaku_txt", "Rigaku text export (angle and intensities on each line)", ["txt"], None, "GSAS-II G2pwd_rigaku.py", rigaku::recognise_txt, rigaku::read_txt, &[WAVELENGTH]),
    entry!("xye", "x, y, σ columns (Topas .xye, GSAS-II's rules)", ["xye", "xy", "dat", "txt"], None, "GSAS-II G2pwd_xye.py", table::recognise_xye, table::read_xye, &[AXIS, WAVELENGTH, TOF_UNIT]),
    entry!("qye", "Q, y, σ columns, wavelength in a comment (.qye)", ["qye"], None, "GSAS-II G2pwd_xye.py", table::recognise_qye, table::read_qye),
    entry!("fit2d_chi", "Fit2D or pyFAI .chi (2θ or Q)", ["chi", "qchi"], None, "GSAS-II G2pwd_xye.py", table::recognise_chi, table::read_chi, &[WAVELENGTH]),
    entry!("powder_cif", "Powder CIF (pdCIF)", ["cif"], None, "GSAS-II G2pwd_CIF.py; the pdCIF dictionary", cif_powder::recognise, cif_powder::read, &[BLOCK, WAVELENGTH]),
    Format {
        id: "nexus_nxazint",
        name: "NeXus/HDF5 azimuthal integration (MAX IV NXazint1d, NXazint2d)",
        extensions: &["h5", "hdf", "hdf5", "nxs"],
        ins: None,
        source: "GSAS-II G2pwd_HDF5.py",
        not_read: Some("an HDF5 decoder is not part of the core yet: the file is recognised by its HDF5 signature only"),
        recognise: hdf::recognise_hdf5,
        read: unread,
        options: &[],
    },
    Format {
        id: "midas_zarr",
        name: "MIDAS zarr (zipped lineouts)",
        extensions: &["zarr.zip"],
        ins: None,
        source: "GSAS-II G2pwd_MIDAS.py",
        not_read: Some("zarr's Blosc codecs are not part of the core yet: the file is recognised by its zip entries only"),
        recognise: hdf::recognise_zarr,
        read: unread,
        options: &[],
    },
];

fn never(_: &Probe) -> f64 {
    0.0
}

fn unread(p: &Probe) -> Result<Dataset, DataError> {
    let f = FORMATS.iter().find(|f| f.not_read.is_some() && (f.recognise)(p) > 0.0);
    Err(DataError(match f {
        Some(f) => format!("{} ({}) is recognised but not read: {}", f.name, f.id, f.not_read.unwrap_or("")),
        None => "this format is recognised but not read".into(),
    }))
}

/// The format with this id.
pub use write::{format_of, point_lines, write, NO_WEIGHT};

pub fn by_id(id: &str) -> Option<&'static Format> {
    FORMATS.iter().find(|f| f.id == id)
}

/// The format FullProf reads with this `Ins` (the `.pcr` translation).
pub fn by_ins(ins: i32) -> Option<&'static Format> {
    FORMATS.iter().find(|f| f.ins == Some(ins))
}

/// The candidates for a file, by content and then extension.
pub fn identify(p: &Probe) -> Identification {
    let ext = p.extension();
    let mut candidates: Vec<Candidate> = FORMATS
        .iter()
        .filter_map(|f| {
            let content = (f.recognise)(p);
            (content > 0.0).then(|| {
                let extension = f.extensions.contains(&ext.as_str());
                Candidate { id: f.id, name: f.name, content, extension, score: content + if extension { EXTENSION } else { 0.0 } }
            })
        })
        .collect();
    candidates.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap_or(std::cmp::Ordering::Equal));
    Identification { candidates, preferred: None }
}

/// The candidates for a file with an order of preference among the readers (the app's request 51): when no reader
/// leads by `MARGIN`, the first of the tied readers in `preference` reads the file, and `Identification::preferred`
/// says so. A preference never overrides a reader that leads.
pub fn identify_with(p: &Probe, preference: &[&str]) -> Identification {
    let mut id = identify(p);
    let tied: Vec<&'static str> = id.tied().iter().map(|c| c.id).collect();
    id.preferred = preference.iter().find_map(|w| tied.iter().find(|t| *t == w).copied());
    id
}

/// Read a file as the given format, or as the format it is identified as.
pub fn read(p: &Probe, format: Option<&str>) -> Result<Dataset, DataError> {
    read_with(p, format, &options::Options::new(), &[])
}

/// Read a file as the given format, or as the format it is identified as with an order of preference among the
/// readers (`identify_with`), with the reader's options (`options`; checked against the format's declarations).
pub fn read_with(p: &Probe, format: Option<&str>, opts: &options::Options, preference: &[&str]) -> Result<Dataset, DataError> {
    let f = match format {
        Some(id) => by_id(id).ok_or_else(|| DataError(format!("no data format with id {id} (drs data formats lists them)")))?,
        None => {
            let id = identify_with(p, preference);
            match id.format() {
                Some(c) => by_id(c.id).expect("candidate in the registry"),
                None if id.candidates.is_empty() => return Err(DataError("no reader recognises this file".into())),
                None => {
                    let names: Vec<String> = id.ambiguous().iter().map(|c| format!("{} ({:.2})", c.id, c.score)).collect();
                    return Err(DataError(format!("ambiguous: the file reads as {}; choose one", names.join(", "))));
                }
            }
        }
    };
    if let Some(why) = f.not_read {
        return Err(DataError(format!("{} ({}) is recognised but not read: {why}", f.name, f.id)));
    }
    options::validate(f, opts).map_err(DataError)?;
    // The general table's layout, when given, is the reader's own choice; the rest applies after the reading.
    let mut d = if f.id == "table" && options::has_table_layout(opts) {
        table::read_with(p.text(), &options::table_options(opts))?.0
    } else {
        (f.read)(p)?
    };
    options::apply(&mut d, opts)?;
    // A WinPLOTR-made header's temperature (PL32), whichever reader took the columns after it.
    if p.is_text() {
        let w = winplotr_temperature(p.text());
        for s in d.scans.iter_mut().filter(|s| s.info.temperature.is_none()) {
            match &w {
                Some(WinplotrTemperature { kelvin: Some(k), said }) => {
                    s.info.temperature = Some(*k);
                    s.info.other.push(("temperature as written".into(), said.clone()));
                }
                Some(WinplotrTemperature { kelvin: None, said }) => s.info.other.push(("temperature not read".into(), said.clone())),
                None => {}
            }
        }
    }
    // A σ of 10³⁰ or more is a point with no weight, as `write` writes it (L088).
    for s in &mut d.scans {
        for v in &mut s.pattern.sigma {
            if *v >= NO_WEIGHT {
                *v = f64::INFINITY;
            }
        }
    }
    Ok(d)
}

/// The temperature line of a header that WinPLOTR wrote (PL32), and what is known of its unit.
#[derive(Clone, Debug, PartialEq)]
pub struct WinplotrTemperature {
    /// In K; None when the unit is not known, for the line is never read as kelvin by default.
    pub kelvin: Option<f64>,
    /// The line's value as written, with its unit and where the unit comes from, or why it was not read.
    pub said: String,
}

/// A header written by WinPLOTR ("File created by WinPLOTR", "File saved from WinPLOTR-2006"): its leading comment
/// lines (those before the first line of numbers, up to 40) searched for `! TEMPERATURE <value> [unit]`.
///
/// WinPLOTR writes the value without a unit (`manuals/fullprof/runs/modes/E/t7.dat`: `! TEMPERATURE 179.99`). Its
/// unit is taken from the line when one follows the number (K, C, °C, Celsius), else from the header's source file:
/// one named `! INITIAL .UXD FILE:` is Bruker's, whose temperatures are degrees Celsius (Josu, L243: case E is
/// 179.99 °C = 453.14 K; Bruker's BRML states °C, `bruker.rs`). Any other header's value is not read, and `said` says
/// why: no kelvin by default (PL32).
pub fn winplotr_temperature(text: &str) -> Option<WinplotrTemperature> {
    let header: Vec<&str> = lines(text).into_iter().take(40).take_while(|l| all_numbers(l).is_none()).collect();
    let low: Vec<String> = header.iter().map(|l| l.to_lowercase()).collect();
    if !low.iter().any(|l| l.contains("winplotr") && (l.contains("created by") || l.contains("saved from"))) {
        return None;
    }
    let line = header.iter().map(|l| l.trim_start_matches('!').trim()).find(|l| l.to_uppercase().starts_with("TEMPERATURE"))?;
    let rest = line["TEMPERATURE".len()..].trim_start_matches([' ', '\t', ':', '=']);
    let mut words = rest.split_whitespace();
    let value = words.next().and_then(fullprof::fnum)?;
    let unit = words.next().map(|u| u.trim_start_matches('°').to_uppercase());
    let uxd = header.iter().find(|l| l.to_uppercase().contains("INITIAL .UXD FILE"));
    let said = |unit: &str, from: &str| format!("{value} {unit} (WinPLOTR's header{from})");
    Some(match unit.as_deref() {
        Some("K") => WinplotrTemperature { kelvin: Some(value), said: said("K", "") },
        Some("C") | Some("CELSIUS") => WinplotrTemperature { kelvin: Some(value + 273.15), said: said("°C", "") },
        Some(u) => WinplotrTemperature { kelvin: None, said: format!("{value} {u}: WinPLOTR's header gives a unit the core does not know") },
        None if uxd.is_some() => WinplotrTemperature { kelvin: Some(value + 273.15), said: said("°C", ", from a Bruker .UXD file, whose temperatures are °C: L243") },
        None => WinplotrTemperature {
            kelvin: None,
            said: format!("{value}: WinPLOTR's header gives no unit and names no source that gives one, and the value is not taken as kelvin (PL32)"),
        },
    })
}

/// Read a file from disk (see `read`).
pub fn read_path(path: &std::path::Path, format: Option<&str>) -> Result<Dataset, DataError> {
    let bytes = std::fs::read(path).map_err(|e| DataError(format!("{}: {e}", path.display())))?;
    let name = path.to_string_lossy();
    read(&Probe::new(&name, &bytes), format)
}

/// `read_with` on a file.
pub fn read_path_with(path: &std::path::Path, format: Option<&str>, opts: &options::Options, preference: &[&str]) -> Result<Dataset, DataError> {
    let bytes = std::fs::read(path).map_err(|e| DataError(format!("{}: {e}", path.display())))?;
    let name = path.to_string_lossy();
    read_with(&Probe::new(&name, &bytes), format, opts, preference)
}

/// Lines of text, the line ends of any system removed.
pub(crate) fn lines(t: &str) -> Vec<&str> {
    t.lines().map(|l| l.trim_end_matches('\r')).collect()
}

/// The numbers of a whitespace-separated line, or None if any field is not a number.
pub(crate) fn all_numbers(line: &str) -> Option<Vec<f64>> {
    let v: Option<Vec<f64>> = line.split_whitespace().map(|s| s.parse::<f64>().ok()).collect();
    v.filter(|v| !v.is_empty())
}
