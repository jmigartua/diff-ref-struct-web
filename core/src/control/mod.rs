//! The control file of diff-ref-struct (P47, Q3; L066): the model, its data and its refinement written as text a
//! person reads and diffs, read back exactly. Files end `.drc`; the syntax and the extension accepted by Josu (L070).
//!
//! - `syntax`: the file's syntax, a strict subset of TOML 1.0, parsed with the place of every key and value, its
//!   comments and its order kept; every error with its line, column and span (the editor's diagnostics) and an
//!   outline of its tables (folding).
//! - `schema`: what the file holds (the parameter dictionary's names, explicit units, readable choices) and the
//!   canonical writer.
//! - `update`: a file written again after the model changed, its comments, order and formatting kept: only the values
//!   that changed are rewritten, what is new is added to its table, what is gone is removed.
//! - `pcr_export`: a `.pcr` file for FullProf, where FullProf can express the model.
//! - `relocate`: a file moved to another folder, every path it holds rewritten to reach the same files (app 45).

pub mod complete;
pub mod locate;
pub mod pcr_export;
pub mod relocate;
pub mod schema;
pub mod syntax;
pub mod update;

pub use relocate::relocate;
use crate::data_ops::Operation;
use crate::powder::Model;
use crate::refinable::Refinement;

/// The pattern's data: the file (as written, relative to the control file), the registry's reader with the options
/// it read the file with (the app's request 51: `[pattern] reader_options`), and the operations applied to it, in
/// order.
#[derive(Clone, Debug, PartialEq)]
pub struct DataSource {
    pub file: String,
    pub reader: Option<String>,
    pub reader_options: crate::formats::options::Options,
    pub operations: Vec<Operation>,
}

/// What an imported `.pcr` wrote that the default model does not use, kept for the export back to FullProf (L070):
/// the angular range (Thmin, Step, Thmax) and the excluded regions as written. The import turns FullProf's regions
/// into the data points FullProf leaves out (Co-YBCO's 0–5 becomes 0.069–5.019), which cannot be turned back into
/// the numbers that were written; the export writes these while the model's regions are still the ones they give.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FullProfRecord {
    pub range: Option<[f64; 3]>,
    pub excluded: Vec<(f64, f64)>,
}

impl FullProfRecord {
    pub fn is_empty(&self) -> bool {
        self.range.is_none() && self.excluded.is_empty()
    }
}

/// Everything a control file holds.
#[derive(Clone, Debug)]
pub struct Control {
    pub title: String,
    pub model: Model,
    pub data: Option<DataSource>,
    pub refinement: Refinement,
    pub fullprof: FullProfRecord,
    /// The file that holds the background points, relative to the control file, when they are not written in it
    /// (R26, L078).
    pub background_file: Option<String>,
}

impl Control {
    /// The canonical text (`schema::write`).
    pub fn to_text(&self) -> String {
        schema::write(self)
    }

    /// Read a control file's text; the diagnostics when it does not read.
    pub fn from_text(text: &str) -> Result<Control, Vec<syntax::Diagnostic>> {
        match schema::read_text(text) {
            (Some(c), _) => Ok(c),
            (None, d) => Err(d),
        }
    }

    /// Read a control file's text in its folder (a background-points file it names is read from there).
    pub fn from_text_in(text: &str, dir: &std::path::Path) -> Result<Control, Vec<syntax::Diagnostic>> {
        match schema::read_with(&syntax::parse(text), Some(dir)) {
            (Some(c), _) => Ok(c),
            (None, d) => Err(d),
        }
    }

    /// The background points as their own file: (path relative to the control file, text), when the control file
    /// names one; the format is the extension's (`.bac` FullProf's first format at the data's points `x`, as FullProf
    /// needs a value at every data point; otherwise ours). A file's b1–b4 (Nba = 1, `BackgroundTerms::file`) apply
    /// where it is read: its points are written as held, at their own spacing, whatever `x` (the app's request 48;
    /// written at `x` they would carry b1–b4 and be moved by them again when read).
    pub fn background_file_text(&self, x: Option<&[f64]>) -> Option<Result<(String, String), String>> {
        let x = if self.model.background_terms.file.is_some() { None } else { x };
        let f = self.background_file.as_ref()?;
        let (pts, how) = match &self.model.background {
            crate::powder::BackgroundModel::Points(p) => (p, crate::background::Interpolation::Linear),
            crate::powder::BackgroundModel::Spline(p) => (p, crate::background::Interpolation::Spline),
            _ => return None,
        };
        let text = match (f.to_lowercase().ends_with(".bac"), x) {
            // At the data's points, the background as the model gives it there (with `Convention::background_edge`).
            (true, Some(x)) => crate::background::regular(x)
                .ok_or_else(|| "a .bac holds the background at equally spaced points: the data's are not".to_string())
                .map(|(start, step)| crate::background::write_bac(&self.model.background_at(x), start, step, &self.title)),
            (true, None) => crate::background::bac_text(pts, how, None, &self.title),
            (false, _) => Ok(crate::background::write_ours(pts, how)),
        };
        Some(text.map(|t| (f.clone(), t)).map_err(|e| format!("{f}: {e}")))
    }

    /// The background file a save writes, as `background_file_text`, or `Ok(None)` when the file already in `dir` (the
    /// control file's folder) holds the same points: a save that changed nothing in it leaves it, FullProf's own `.bac`
    /// among them, as it is (the app's request 48).
    pub fn background_file_to_write(&self, dir: &std::path::Path, x: Option<&[f64]>) -> Option<Result<Option<(String, String)>, String>> {
        let r = self.background_file_text(x)?;
        Some(r.map(|(f, t)| {
            let on_disk = crate::background::read_path(&dir.join(&f)).ok();
            let new = crate::background::read(&t).ok();
            let same = matches!((&on_disk, &new), (Some(a), Some(b)) if a.points == b.points && (a.interpolation.is_none() || a.interpolation == b.interpolation));
            (!same).then_some((f, t))
        }))
    }

    /// A `.pcr` imported into the default model (as `drs calc` imports it), with its data file named relative to
    /// `base` (the folder the control file will be in) and its refinement.
    pub fn from_pcr(pcr: &crate::pcr::Pcr, dir: &std::path::Path, stem: &str, base: &std::path::Path) -> Result<Control, String> {
        let l = Model::from_pcr(pcr, dir, stem)?;
        let refinement = Refinement::from_pcr(pcr, &l.model)?;
        let data = match l.data_path.as_ref() {
            Some(p) => {
                // FullProf's Ste > 1: every Ste-th point kept (L267).
                let mut operations = if pcr.flag("Ste") > 1 { vec![crate::data_ops::Operation::RebinFactor { factor: pcr.flag("Ste") as usize, mode: crate::data_ops::RebinMode::Keep }] } else { Vec::new() };
                // A TOF pattern's counts as FullProf reads them (P05, `powder::read_tof_data`).
                if l.model.tof.is_some() {
                    let text: String = std::fs::read(p).map_err(|e| format!("{}: {e}", p.display()))?.iter().map(|&b| b as char).collect();
                    operations.extend(crate::powder::read_tof_data(&text, pcr.flag("Ins"))?.1);
                }
                Some(DataSource { file: relative(p, base), reader: crate::formats::by_ins(pcr.flag("Ins")).map(|f| f.id.to_string()), reader_options: Default::default(), operations })
            }
            None => None,
        };
        let fullprof = FullProfRecord { range: (pcr.step > 0.0).then_some([pcr.thmin, pcr.step, pcr.thmax]), excluded: pcr.excluded.clone() };
        // Nba = 1: the points stay in the .bac file the .pcr reads, as read, which the control file names; b1–b4 are
        // the background's values `background_file_scale` … `_shift` (eq. 3.6; the app's request 47).
        let background_file = (pcr.flag("Nba") == 1).then(|| crate::powder::find_file(dir, &format!("{stem}.bac")).map(|p| relative(&p, base))).flatten();
        // MORE = 2, 3: the profiles' files relative to the control file (P23).
        let mut model = l.model;
        if let (Some(f), Some(name)) = (&mut model.background_terms.file, &background_file) {
            f.file = name.clone();
        }
        for p in &mut model.background_terms.profiles {
            p.file = relative(std::path::Path::new(&p.file), base);
        }
        Ok(Control { title: pcr.title.trim().to_string(), model, data, refinement, fullprof, background_file })
    }
}

/// `path` relative to the folder `base` ("../data/x.dat"), with '/' as separator; absolute when they share no root
/// (another drive).
pub fn relative(path: &std::path::Path, base: &std::path::Path) -> String {
    let abs = |p: &std::path::Path| std::fs::canonicalize(p).unwrap_or_else(|_| p.to_path_buf());
    let (p, b) = (abs(path), abs(base));
    let pc: Vec<_> = p.components().collect();
    let bc: Vec<_> = b.components().collect();
    let common = pc.iter().zip(&bc).take_while(|(a, b)| a == b).count();
    // Absolute only without a common root (Windows' drives): two paths that share "/" alone, a run in /tmp and a
    // project under /Users, are relative too (Publish's P61, L374).
    let shared_root = pc[..common].iter().any(|c| matches!(c, std::path::Component::RootDir));
    if if p.has_root() || b.has_root() { !shared_root } else { common == 0 } {
        return p.to_string_lossy().replace('\\', "/");
    }
    let mut parts: Vec<String> = vec!["..".to_string(); bc.len() - common];
    parts.extend(pc[common..].iter().map(|c| c.as_os_str().to_string_lossy().into_owned()));
    parts.join("/")
}

/// A control file opened: the control, and its data (the reader's first scan, its operations applied) with the data
/// file's path.
pub struct Opened {
    pub control: Control,
    pub data: Option<crate::data::Pattern>,
    /// The scan the data come from, after the operations (its axis, wavelength and what the file says).
    pub scan: Option<crate::formats::Scan>,
    pub data_path: Option<std::path::PathBuf>,
}

/// Open a control file and its data (`data` relative to the control file's folder; the operations' files too).
pub fn open(path: &std::path::Path) -> Result<Opened, String> {
    open_in(path, &path.parent().map(|d| d.to_path_buf()).unwrap_or_default())
}

/// Open a control file whose files are named relative to another folder, `base` (the app's request 33): a version in
/// `Versions/` is the working file as it was, its data, operations' and background files named from the project's
/// folder, the parent of `Versions/`.
pub fn open_in(path: &std::path::Path, base: &std::path::Path) -> Result<Opened, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let control = Control::from_text_in(&text, base).map_err(|d| d.iter().map(|x| format!("{}:{}:{}: {}", path.display(), x.span.line, x.span.column, x.message)).collect::<Vec<_>>().join("\n"))?;
    let dir = base.to_path_buf();
    let read = |file: &str, reader: Option<&str>| -> Result<crate::formats::Scan, String> {
        let d = crate::formats::read_path(&dir.join(file), reader).map_err(|e| e.to_string())?;
        d.scans.into_iter().next().ok_or(format!("{file}: no scan"))
    };
    let (scan, data_path) = match &control.data {
        Some(src) => {
            let scan = data_scan_with(&dir, &src.file, src.reader.as_deref(), &src.reader_options, control.model.tof.is_some())?;
            (Some(crate::data_ops::apply_all(&scan, &src.operations, &read)?), Some(dir.join(&src.file)))
        }
        None => (None, None),
    };
    Ok(Opened { control, data: scan.as_ref().map(|s| s.pattern.clone()), scan, data_path })
}

/// A control file's data file read (its first scan), before its operations: for a TOF pattern (`tof`) a GSAS bank put on
/// the time axis in µs as GSAS-II reads it (`tof::time_axis`, P05); a file of another format keeps its x as written
/// (alum's and lamn_pol's `xys` are in µs already: the app's request 46, L354).
pub fn data_scan(dir: &std::path::Path, file: &str, reader: Option<&str>, tof: bool) -> Result<crate::formats::Scan, String> {
    data_scan_with(dir, file, reader, &Default::default(), tof)
}

/// `data_scan` with the reader's options (the app's request 51): the scan they choose, else the first.
pub fn data_scan_with(dir: &std::path::Path, file: &str, reader: Option<&str>, options: &crate::formats::options::Options, tof: bool) -> Result<crate::formats::Scan, String> {
    let d = crate::formats::read_path_with(&dir.join(file), reader, options, &[]).map_err(|e| format!("{file}: {e}"))?;
    let gsas = d.format == "gsas";
    let mut scan = d.scans.into_iter().next().ok_or(format!("{file}: no scan"))?;
    if tof && gsas {
        crate::tof::time_axis(&mut scan);
    }
    Ok(scan)
}

/// Is this file a control file (by its name: `.drc`, L070)?
pub fn is_control_file(path: &std::path::Path) -> bool {
    path.extension().is_some_and(|e| e.eq_ignore_ascii_case("drc"))
}
