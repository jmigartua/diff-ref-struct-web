//! A control file pointed at new data (R25, APP_WORKFLOWS W3; L078): the refinement of other data taken as the start
//! of this one. The control is copied, never edited where it was; its data become the new file, and what was written
//! for the old data and does not fit the new is reported as the checks of R24 report (`check::Diagnostic`), each with
//! a suggestion: another wavelength or axis, the range, another instrument or step, excluded regions and background
//! points outside the new data or of other counts, a Chebyshev range of other angles, the old data's operations.
//! Then the checks themselves run on the result (`check::check_path`).

use crate::check::{Diagnostic, Level, Severity};
use crate::control::{relative, Control, DataSource, FullProfRecord};
use crate::formats::{Axis, Scan};
use crate::powder::BackgroundModel;
use std::path::Path;

/// The rebased control and what does not fit.
pub struct Rebased {
    pub control: Control,
    pub diagnostics: Vec<Diagnostic>,
}

fn step(s: &Scan) -> Option<f64> {
    let p = &s.pattern;
    (p.len() >= 2).then(|| (p.x[p.len() - 1] - p.x[0]) / (p.len() - 1) as f64)
}

/// Point a control at a new scan: `file` and `reader` as the control file will name them; `old` the scan the control
/// was written for, when it is at hand (for the range, instrument and step).
pub fn rebase(c: &Control, old: Option<&Scan>, new: &Scan, file: String, reader: Option<String>) -> Rebased {
    let mut out = Vec::new();
    let mut r = c.clone();
    let m = &c.model;
    let p = &new.pattern;
    let (lo, hi) = (p.x.first().copied().unwrap_or(0.0), p.x.last().copied().unwrap_or(0.0));
    let st = step(new).unwrap_or(0.0).abs();
    // The data and its operations.
    if let Some(d) = &c.data {
        if !d.operations.is_empty() {
            let names: Vec<&str> = d.operations.iter().map(|o| o.name()).collect();
            out.push(Diagnostic::new(Severity::Warning, Level::Model, "rebase-operations", format!("the old data's {} operation(s) ({}) are not carried over", names.len(), names.join(", ")), "masks, deleted points, smoothing and the rest were chosen for the old file's points").addr("pattern1.data").fix("apply again those the new data need (drs data op …)"));
        }
    }
    r.data = Some(DataSource { file, reader, reader_options: Default::default(), operations: Vec::new() });
    r.fullprof = FullProfRecord::default();
    // The axis and the wavelength.
    if new.axis != Axis::TwoTheta {
        out.push(Diagnostic::new(Severity::Error, Level::Model, "rebase-axis", format!("the new data's x is {}, not 2θ", new.axis.label()), "the model calculates a constant-wavelength pattern in 2θ").addr("pattern1.data").fix("convert the data to 2θ (drs data op … convert to=2theta wavelength=…)"));
    }
    match new.wavelength {
        Some(w) if (w.lambda1 - m.lambda[0]).abs() > 1e-4 => {
            let mut d = Diagnostic::new(Severity::Warning, Level::Model, "rebase-wavelength", format!("the new data were measured at λ1 = {} Å; the control file has {} Å", w.lambda1, m.lambda[0]), "every peak position follows from λ = 2d sin θ: with the old wavelength the peaks fall in the wrong places").addr("pattern1.wavelength_1").entry("wavelength_1").fix(format!("write wavelength_1 = {}{}", w.lambda1, w.lambda2.map_or(String::new(), |l2| format!(", wavelength_2 = {l2}"))));
            d.replacement = Some(format!("wavelength_1 = {}", w.lambda1));
            out.push(d);
        }
        None => out.push(Diagnostic::new(Severity::Hint, Level::Model, "rebase-wavelength", format!("the new data file does not state its wavelength; the control file's λ1 = {} Å is kept", m.lambda[0]), "few data formats carry the wavelength (FullProf's do not)").addr("pattern1.wavelength_1").entry("wavelength_1").fix("check it against the instrument the data come from")),
        _ => {}
    }
    // The range, the instrument and the step.
    if let Some(o) = old {
        let (olo, ohi) = (o.pattern.x.first().copied().unwrap_or(0.0), o.pattern.x.last().copied().unwrap_or(0.0));
        if (olo - lo).abs() > st || (ohi - hi).abs() > st {
            out.push(Diagnostic::new(Severity::Hint, Level::Model, "rebase-range", format!("the new data run from {lo}° to {hi}°, the old from {olo}° to {ohi}°"), "the excluded regions, the background and the refined profile were set for the old range").addr("pattern1.data").fix("look at the ends of the new pattern before refining"));
        }
        if let (Some(a), Some(b)) = (&o.info.instrument, &new.info.instrument) {
            if a.trim() != b.trim() {
                out.push(Diagnostic::new(Severity::Warning, Level::Model, "rebase-instrument", format!("the new data come from {b}, the old from {a}"), "the zero, the resolution (U V W, a resolution file) and the asymmetry belong to an instrument").addr("pattern1.data").entry("zero_shift").entry("caglioti_u").fix("take the instrument's zero and resolution from a standard measured on it"));
            }
        }
        if let (Some(a), Some(b)) = (step(o), step(new)) {
            if ((a - b) / a).abs() > 0.01 {
                out.push(Diagnostic::new(Severity::Hint, Level::Model, "rebase-step", format!("the new data's step is {b:.5}°, the old {a:.5}°"), "another step often means another instrument setting: slits, detector, counting").addr("pattern1.data").fix("check the instrument's parameters"));
            }
        }
    }
    // Excluded regions outside the new data.
    for (k, &(a, b)) in m.excluded.iter().enumerate() {
        if b < lo || a > hi {
            out.push(Diagnostic::new(Severity::Warning, Level::Model, "rebase-excluded", format!("excluded region {} ({a}–{b}°) is outside the new data ({lo}–{hi}°)", k + 1), "it was written for the old data and leaves out nothing here").addr("pattern1.excluded_regions").entry("excluded_regions").fix("remove it, or exclude what the new data need"));
        } else if a <= lo && b >= hi {
            out.push(Diagnostic::new(Severity::Error, Level::Model, "rebase-excluded", format!("excluded region {} ({a}–{b}°) covers all the new data", k + 1), "no point would be left to fit").addr("pattern1.excluded_regions").entry("excluded_regions").fix("remove it"));
        }
    }
    // The background.
    match &m.background {
        BackgroundModel::Points(pts) | BackgroundModel::Spline(pts) if !pts.is_empty() => {
            let (p0, pn) = (pts[0].0, pts[pts.len() - 1].0);
            // What the points leave uncovered at each end, beyond what they left of the old data (or, without the old
            // data, beyond 1 % of the new range).
            let (olo, ohi) = old.map_or((p0, pn), |o| (o.pattern.x.first().copied().unwrap_or(p0), o.pattern.x.last().copied().unwrap_or(pn)));
            let allow = 0.01 * (hi - lo).abs() + st;
            if (p0 - lo) - (p0 - olo).max(0.0) > allow || (hi - pn) - (ohi - pn).max(0.0) > allow {
                out.push(Diagnostic::new(Severity::Warning, Level::Model, "rebase-background", format!("the background points run from {p0}° to {pn}°, the new data from {lo}° to {hi}°"), "outside its points the background is held at the first or last point's value").addr("pattern1.background_points").entry("background_points").fix("add points at the ends, or pick them anew on the new data (drs data background)"));
            }
            // The points' counts against the new data's low points near each of them.
            let ratios: Vec<f64> = pts
                .iter()
                .filter(|q| q.0 >= lo && q.0 <= hi && q.1 > 0.0)
                .filter_map(|q| {
                    let near: Vec<f64> = p.x.iter().zip(&p.y).filter(|(x, _)| (*x - q.0).abs() <= 10.0 * st.max(1e-6)).map(|(_, y)| *y).collect();
                    let low = near.iter().copied().fold(f64::INFINITY, f64::min);
                    low.is_finite().then(|| low / q.1)
                })
                .collect();
            if !ratios.is_empty() {
                let mut s = ratios.clone();
                s.sort_by(f64::total_cmp);
                let med = s[s.len() / 2];
                if !(0.5..=2.0).contains(&med) {
                    out.push(Diagnostic::new(Severity::Warning, Level::Model, "rebase-background", format!("the new data are about {med:.2} times the background points' counts"), "the points were picked on other data: other counting, sample or instrument").addr("pattern1.background_points").entry("background_points").fix("pick the points anew on the new data, or scale them (the automatic background gives a start: drs data background)"));
                }
            }
        }
        BackgroundModel::Chebyshev { from, to, .. } if (from - lo).abs() > st || (to - hi).abs() > st => {
            let mut d = Diagnostic::new(Severity::Warning, Level::Model, "rebase-background", format!("the Chebyshev background spans {from}–{to}°, the new data {lo}–{hi}°"), "its terms are independent over the range they span; over other angles they become nearly dependent (L076)").addr("pattern1.background_chebyshev").entry("background_chebyshev").fix(format!("write from = {lo} and to = {hi}, and refine the terms again"));
            d.replacement = Some(format!("from = {lo}\nto = {hi}"));
            out.push(d);
        }
        _ => {}
    }
    Rebased { control: r, diagnostics: out }
}

/// `drs rebase`: the control file (`.drc` or `.pcr`) pointed at `data` and written to `to`, its background file
/// copied next to it, and every diagnostic: those of the rebase, then the checks of the written file.
pub fn rebase_path(file: &Path, data: &Path, reader: Option<&str>, to: &Path) -> Result<Rebased, String> {
    let to_dir = to.parent().filter(|d| !d.as_os_str().is_empty()).map(|d| d.to_path_buf()).unwrap_or_else(|| std::path::PathBuf::from("."));
    std::fs::create_dir_all(&to_dir).map_err(|e| format!("{}: {e}", to_dir.display()))?;
    let ext = file.extension().map(|e| e.to_string_lossy().to_lowercase()).unwrap_or_default();
    let (control, old, from_dir) = if ext == "pcr" {
        let text: String = std::fs::read(file).map_err(|e| format!("{}: {e}", file.display()))?.iter().map(|&b| b as char).collect();
        let mut pcr = crate::pcr::read(&text).map_err(|e| format!("{}: {e}", file.display()))?;
        pcr.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
        let dir = std::fs::canonicalize(file).ok().and_then(|p| p.parent().map(|d| d.to_path_buf())).unwrap_or_default();
        let stem = file.file_stem().and_then(|s| s.to_str()).unwrap_or("");
        let c = Control::from_pcr(&pcr, &dir, stem, &dir)?;
        let old = c.data.as_ref().and_then(|d| crate::formats::read_path(&dir.join(&d.file), d.reader.as_deref()).ok()).and_then(|d| d.scans.into_iter().next());
        (c, old, dir)
    } else {
        let o = crate::control::open(file)?;
        let dir = file.parent().map(|d| d.to_path_buf()).unwrap_or_default();
        (o.control, o.scan, dir)
    };
    let bytes = std::fs::read(data).map_err(|e| format!("{}: {e}", data.display()))?;
    let name = data.to_string_lossy().to_string();
    let probe = crate::formats::Probe::new(&name, &bytes);
    let id = match reader {
        Some(r) => r.to_string(),
        None => crate::formats::identify(&probe).format().map(|c| c.id.to_string()).ok_or(format!("{name}: the format is not recognised or is ambiguous (drs data identify); give it with --as"))?,
    };
    let new = crate::formats::read(&probe, Some(&id)).map_err(|e| format!("{name}: {e}"))?.scans.into_iter().next().ok_or(format!("{name}: no scan"))?;
    let mut r = rebase(&control, old.as_ref(), &new, relative(data, &to_dir), Some(id));
    // The background points' file, copied next to the new control file under its own name.
    if let Some(f) = r.control.background_file.clone() {
        let src = from_dir.join(&f);
        let name = Path::new(&f).file_name().map(|n| n.to_string_lossy().to_string()).unwrap_or(f.clone());
        std::fs::copy(&src, to_dir.join(&name)).map_err(|e| format!("{}: {e}", src.display()))?;
        r.control.background_file = Some(name);
    }
    std::fs::write(to, r.control.to_text()).map_err(|e| format!("{}: {e}", to.display()))?;
    let mut checks = crate::check::check_path(to, &crate::check::Options::default());
    // A check that says what a rebase diagnostic already says is not repeated.
    checks.retain(|d| !(d.code == "wavelength" && r.diagnostics.iter().any(|x| x.code == "rebase-wavelength")));
    r.diagnostics.extend(checks);
    Ok(r)
}
