//! Sequential refinement (P54, L317): a series of patterns refined in turn with one control file, as FullProf's
//! sequential mode runs them (its log of the D20 run in `manuals/fullprof/1-WinPLOTR/2_FP_sequential_refinements/
//! LaMnO3_D20/la0.log`: each data file refined from the values the one before ended at, to convergence or the
//! cycles allowed, its refined state written as `la0_NNN.pcr` and its summary appended to the `.seq`). A series may
//! also start every pattern from the control file's values (`Start::Each`), so that one pattern's failure does not
//! carry on into the next.
//!
//! Each pattern is read as the control file reads its own (`control::open_in`): the same reader and the same data
//! operations; its external parameter is the temperature its header gives (`formats::Info::temperature`, FullProf's
//! PARVAL), or one the caller states. A control file with a plan of steps (L110) runs its plan on every pattern
//! (`lsq::run_plan`). The outputs of a series (P59): `series_csv`, a table of every point's agreement and refined
//! values with their su's and each phase's cell volume, and `series_cif`, every point's results file in one CIF.

use crate::control::Control;
use crate::data::Pattern;
use crate::lsq::{Cycle, Options, Refinement};
use crate::powder::Model;
use crate::refinable::ParameterSet;
use std::fmt::Write;
use std::path::Path;

/// Where each pattern of a series starts from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Start {
    /// From the values the pattern before ended at (FullProf's sequential mode); the first from the control file's.
    Previous,
    /// Every pattern from the control file's values.
    Each,
}

/// One pattern of a series: its label (the data file's name without its folder and extension), its external parameter
/// (a temperature, a pressure), and what it ended at.
#[derive(Clone, Debug)]
pub struct Point {
    pub label: String,
    pub parameter: Option<f64>,
    pub data: Pattern,
    /// The model at the values reached, its free parameters, and the refinement (statistics, esd's, covariance).
    pub model: Model,
    pub set: ParameterSet,
    pub refinement: Refinement,
}

/// A pattern of the series, read: its label, the pattern, and the temperature its header gives.
#[derive(Clone, Debug)]
pub struct Input {
    pub label: String,
    pub data: Pattern,
    pub parameter: Option<f64>,
}

/// Read a data file as the control file reads its own: its reader and data operations (`control::DataSource`), the
/// other files the operations name relative to `base`. The label is the file's stem; the parameter, the temperature
/// its header gives.
pub fn read(control: &Control, base: &Path, file: &Path) -> Result<Input, String> {
    let reader = control.data.as_ref().and_then(|d| d.reader.clone());
    let ops = control.data.as_ref().map_or(Vec::new(), |d| d.operations.clone());
    let read = |f: &str, r: Option<&str>| -> Result<crate::formats::Scan, String> {
        let d = crate::formats::read_path(&base.join(f), r).map_err(|e| e.to_string())?;
        d.scans.into_iter().next().ok_or(format!("{f}: no scan"))
    };
    let path = if file.is_absolute() { file.to_path_buf() } else { base.join(file) };
    let d = crate::formats::read_path(&path, reader.as_deref()).map_err(|e| format!("{}: {e}", path.display()))?;
    let scan = d.scans.into_iter().next().ok_or(format!("{}: no scan", path.display()))?;
    let scan = crate::data_ops::apply_all(&scan, &ops, &read)?;
    let label = path.file_stem().map_or(String::new(), |s| s.to_string_lossy().to_string());
    Ok(Input { label, parameter: scan.info.temperature, data: scan.pattern })
}

/// Refine the patterns in turn with `control`'s model and refinement (its plan, when it has one), each from the
/// values the pattern before ended at or from the control file's (`start`). `on_cycle` receives the point's index
/// and each cycle, and stops the series when it returns false (the point then ends there, and no other is run).
pub fn run(control: &Control, inputs: &[Input], start: Start, mut on_cycle: impl FnMut(usize, &Cycle) -> bool) -> Result<Vec<Point>, String> {
    let mut out: Vec<Point> = Vec::new();
    let rf = &control.refinement;
    for (k, input) in inputs.iter().enumerate() {
        let mut model = match (start, out.last()) {
            (Start::Previous, Some(p)) => p.model.clone(),
            _ => control.model.clone(),
        };
        let mut go_on = true;
        let (refinement, set) = if rf.plan.is_empty() {
            let set = rf.build(&model)?;
            let opts = Options::stated(rf.cycles, rf.eps);
            let r = crate::lsq::refine_until(&mut model, &input.data, &set, opts, |c| {
                go_on = on_cycle(k, c);
                go_on
            });
            (r, set)
        } else {
            let steps: Vec<usize> = (0..rf.plan.len()).collect();
            let (r, session) = crate::lsq::run_plan(&mut model, &input.data, rf, &steps, |_| {}, |c| {
                go_on = on_cycle(k, c);
                go_on
            })?;
            (r, session.set)
        };
        out.push(Point { label: input.label.clone(), parameter: input.parameter, data: input.data.clone(), model, set, refinement });
        if !go_on {
            break;
        }
    }
    Ok(out)
}

/// A phase's cell volume (Å³) and its su from the covariance of the free parameters (`ParameterSet::cell_volume`, the
/// core's one su of a volume, P43).
pub fn cell_volume(p: &Point, phase: usize) -> (f64, f64) {
    p.set.cell_volume(&p.model, phase, &p.refinement.covariance)
}

fn csv_field(s: &str) -> String {
    if s.contains([',', '"', '\n']) {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}

/// The series as CSV (P59): one row per point, its label, external parameter, cycles, stop, χ², Rp, Rwp and Rexp (%),
/// each phase's cell volume and su (`phase1.cell_volume`, Å³), and every free parameter's value and su, by name, in
/// the order they first appear (empty where a point does not refine it). Numbers to the last bit.
pub fn series_csv(points: &[Point]) -> String {
    let mut names: Vec<String> = Vec::new();
    for p in points {
        for q in &p.refinement.params {
            if !names.contains(&q.name) {
                names.push(q.name.clone());
            }
        }
    }
    let phases = points.iter().map(|p| p.model.phases.len()).max().unwrap_or(0);
    let mut o = String::from("point,label,parameter,cycles,stop,chi2,rp,rwp,rexp");
    for k in 1..=phases {
        let _ = write!(o, ",phase{k}.cell_volume,phase{k}.cell_volume su");
    }
    for n in &names {
        let _ = write!(o, ",{},{} su", csv_field(n), csv_field(n));
    }
    o.push('\n');
    let num = |v: f64| format!("{v:?}");
    for (i, p) in points.iter().enumerate() {
        let s = &p.refinement.stats;
        let _ = write!(o, "{},{},{},{},{},{},{},{},{}", i + 1, csv_field(&p.label), p.parameter.map_or(String::new(), num), p.refinement.cycles, p.refinement.reason.code(), num(s.chi2), num(s.rp), num(s.rwp), num(s.rexp));
        for k in 0..phases {
            if k < p.model.phases.len() {
                let (v, e) = cell_volume(p, k);
                let _ = write!(o, ",{},{}", num(v), num(e));
            } else {
                o.push_str(",,");
            }
        }
        for n in &names {
            match p.refinement.params.iter().find(|q| &q.name == n) {
                Some(q) => {
                    let _ = write!(o, ",{},{}", num(q.value), num(q.esd));
                }
                None => o.push_str(",,"),
            }
        }
        o.push('\n');
    }
    o
}

/// The marker before each point's results file in `series_cif`.
const POINT_MARK: &str = "#==== point ";

/// The series as one CIF (P59): an index block (`data_series`: per point its label, external parameter, χ², Rwp and
/// the `_pd_block_id` of its overall block), then each point's results file (`results::write`, named "{name} {label}"),
/// its data blocks renamed `data_p{n}_…` so that the names stay unique. `split_series` gives each point's text back,
/// which `results::read` reads.
pub fn series_cif(control: &Control, points: &[Point], name: &str) -> Result<String, String> {
    let mut o = format!("# A series refined by diff-ref-struct {} (drs): {} points\n\ndata_series\n_drs_series_name                  '{}'\n", env!("CARGO_PKG_VERSION"), points.len(), name.replace('\'', ""));
    let _ = writeln!(o, "loop_\n_drs_series_point_id\n_drs_series_label\n_drs_series_parameter\n_drs_series_chi2\n_drs_series_wR_factor\n_drs_series_block_id");
    for (i, p) in points.iter().enumerate() {
        let _ = writeln!(o, "{} '{}' {} {:?} {:?} '{} {}|overall'", i + 1, p.label.replace('\'', ""), p.parameter.map_or("?".to_string(), |v| format!("{v:?}")), p.refinement.stats.chi2, p.refinement.stats.rwp / 100.0, name, p.label);
    }
    for (i, p) in points.iter().enumerate() {
        let c = Control { model: p.model.clone(), ..control.clone() };
        let (text, _) = crate::results::write(&c, &p.data, Some(&p.set), Some(&p.refinement), &format!("{name} {}", p.label))?;
        let _ = writeln!(o, "\n{POINT_MARK}{} {}", i + 1, p.label);
        for line in text.lines() {
            match line.strip_prefix("data_") {
                Some(rest) => {
                    let _ = writeln!(o, "data_p{}_{rest}", i + 1);
                }
                None => {
                    let _ = writeln!(o, "{line}");
                }
            }
        }
    }
    Ok(o)
}

/// Each point's results file from a `series_cif` text, in order.
pub fn split_series(text: &str) -> Vec<String> {
    let mut out: Vec<String> = Vec::new();
    for line in text.lines() {
        if line.starts_with(POINT_MARK) {
            out.push(String::new());
        } else if let Some(s) = out.last_mut() {
            s.push_str(line);
            s.push('\n');
        }
    }
    out
}
