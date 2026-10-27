//! A refinement checked against another program (the app's requests 49 and 50, app/DESIGN.md §11.1): the run
//! prepared in a folder of its own, the other program's outputs read, and the two refinements compared in one
//! `Comparison`. The app finds the program, runs the command `Run` gives and shows what this module returns; it writes
//! no input file, reads no output and judges no difference (L139).
//!
//! **FullProf** (`fullprof_prepare`, `fullprof_read`). The control file is exported as a `.pcr` (`pcr_export`, the
//! export FullProf reproduces to its print in eleven runs, `manuals/fullprof/runs/export_check.txt`), with the cycles
//! asked for and FullProf's output flags of `manuals/fullprof/runs/run_fullprof.py` (Prf = 3; Pcr = 1, the refined
//! values written back into the `.pcr`). FullProf is run as `fp2k CODE DATA` in the folder, `FULLPROF` set to its
//! folder and nothing on its input. Read back: the `.out` (the banner, the parameters with their esd's, the
//! statistics, R_Bragg and R_F of the last cycle), the `.sum` (kept for the listing), the `.prf` (the calculated
//! pattern) and the refined `.pcr`, imported into the default model as any `.pcr` is (L050), so that each value is
//! compared in the core's terms. FullProf's names are matched to the core's parameters by `output_names`.
//!
//! **The esd's in the default model.** FullProf's printed esd's are put into the default model's terms by the inverse
//! of the export's translations: the X-ray scale's ×2 (L047), the Bérar–Baldinozzi P1 + 3·P2's σ₁ + 3σ₂ (an upper
//! bound), the axial pair's (σ_S + σ_D)/2, and a codeword's multiplier a taken out (FullProf prints σ·a, manual eq.
//! 3.7, `lsq::fullprof_esd`). Two conventions are left as FullProf has them and said in the notes, both measured in
//! `fullprof_refinements` (L052, L053): FullProf scales its esd's by the χ² over the points with Bragg contributions,
//! the core by that over all the points used; and it prints half the esd of a Thompson–Cox–Hastings phase's cell
//! parameters and I_G, for a reason not known.
//!
//! **GSAS-II** (`gsas2_prepare`, `gsas2_read`, in `compare/gsas2.rs` with the model's translation into GSAS-II's
//! terms): a script for GSASIIscriptable (`compare/gsas2.py`, fixed) that reads the core's choices from a JSON, builds
//! the project from the data, an instrument file and each phase's CIF, refines the cycles and writes what it found as
//! JSON. Both need the `serde` feature.
//!
//! **The verdict.** A parameter agrees when the two values differ by at most the other program's esd (`RULE`): the
//! test `fullprof_refinements` holds the engine to it on FullProf's converged tutorials.

use crate::control::{pcr_export, Opened};
use crate::formats::fullprof_outputs::{read_out, read_prf, Out};
use crate::jacobian::Field;
use crate::lsq::{refine, Options};
use crate::powder::{Model, Shape};
use std::path::{Path, PathBuf};

#[cfg(feature = "serde")]
mod gsas2;
#[cfg(feature = "serde")]
pub use gsas2::{gsas2_prepare, gsas2_read, GSAS2_CODE};

/// The rule of the verdict, written once.
pub const RULE: &str = "|Δ| ≤ 1 σ, σ the other program's esd in the default model's terms";

/// The name of the comparison file in a run's folder (the file dictionary's `comparison`).
pub const COMPARISON_FILE: &str = "comparison.json";

/// The name of the run's record in its folder (`Run::save`, `Run::load`).
pub const RUN_FILE: &str = "run.json";

/// The code (FullProf's CODFIL) of a FullProf run: the `.pcr` and every file FullProf writes start with it.
pub const FULLPROF_CODE: &str = "compare";

/// What became of an item of the control file in the other program's input.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Carried {
    /// Written as the control file has it.
    Carried,
    /// Written, in the other program's nearest terms: the results may differ by it.
    Approximated,
    /// Not written (a free parameter written fixed, a feature the other program has not).
    Refused,
}

/// An item of the manifest: what (a parameter's address, or a feature in words), and why.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ManifestItem {
    pub name: String,
    pub carried: Carried,
    pub why: String,
}

/// A comparison run prepared: the command the app runs and the manifest it shows before it does.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Run {
    /// "FullProf" or "GSAS-II".
    pub other: String,
    /// The program to run: `fp2k`, found by the app; GSAS-II's Python, given by the app.
    pub program: String,
    pub args: Vec<String>,
    /// The folder to run it in.
    pub dir: PathBuf,
    /// The environment the program needs: (variable, what the app sets it to).
    pub env: Vec<(String, String)>,
    /// The refinement cycles both programs run.
    pub cycles: usize,
    /// FullProf's CODFIL, or the GSAS-II script's stem.
    pub code: String,
    /// The data file in the folder, as the program reads it.
    pub data: String,
    pub manifest: Vec<ManifestItem>,
}

#[cfg(feature = "serde")]
impl Run {
    /// Write the run's record into its folder (`RUN_FILE`).
    pub fn save(&self) -> Result<(), String> {
        let p = self.dir.join(RUN_FILE);
        std::fs::write(&p, serde_json::to_string_pretty(self).map_err(|e| e.to_string())?).map_err(|e| format!("{}: {e}", p.display()))
    }

    /// The run prepared in `dir`, from its record; its folder is `dir` wherever the folder was moved.
    pub fn load(dir: &Path) -> Result<Run, String> {
        let p = dir.join(RUN_FILE);
        let t = std::fs::read_to_string(&p).map_err(|e| format!("{}: {e}", p.display()))?;
        let mut r: Run = serde_json::from_str(&t).map_err(|e| format!("{}: {e}", p.display()))?;
        r.dir = dir.to_path_buf();
        Ok(r)
    }
}

/// A program and its version.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Program {
    pub name: String,
    pub version: String,
    /// The cycles it ran, and whether it stopped by its own convergence test.
    pub cycles: usize,
    pub converged: bool,
}

/// A measure of agreement of both refinements: χ², an R factor, or a phase's R_Bragg or R_F.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Measure {
    /// "chi2", "rp", "rwp", "rexp", "rp_background_free", "rwp_background_free", "rexp_background_free",
    /// "durbin_watson", "r_bragg", "r_f".
    pub id: String,
    /// "χ²", "Rp (%)", …
    pub name: String,
    /// The phase (0, 1, …) of R_Bragg and R_F.
    pub phase: Option<usize>,
    pub core: f64,
    /// As the other program prints it; None where it prints none.
    pub other: Option<f64>,
    /// core − other.
    pub difference: Option<f64>,
    /// Half a unit of the other program's last printed digit: a difference within it is the print's.
    pub print: Option<f64>,
}

/// A parameter refined by both, compared.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ParameterRow {
    /// The core's address ("phase1.cell_a").
    pub address: String,
    /// The other program's name for it ("Cell_A_ph1_pat1").
    pub other_name: String,
    pub core_value: f64,
    pub core_esd: f64,
    /// The other's value and esd in the default model's terms.
    pub other_value: f64,
    pub other_esd: f64,
    /// (core − other)/σ_other.
    pub delta_sigma: f64,
    /// σ_core/σ_other.
    pub esd_ratio: f64,
    /// By `RULE`.
    pub agrees: bool,
    /// A convention of the other program's that bears on this row, in words.
    pub note: Option<String>,
}

/// A parameter refined by one program and not compared, and why.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct NotCompared {
    /// The core's address, or the other program's name when the core has no field for it.
    pub name: String,
    /// The program that refined it.
    pub refined_by: String,
    pub why: String,
}

/// Both calculated patterns at the data's points.
#[derive(Clone, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Patterns {
    pub x: Vec<f64>,
    pub yobs: Vec<f64>,
    pub ycalc_core: Vec<f64>,
    /// The other's Ycal at the data's points (interpolated where its points are not the data's; NaN outside them).
    pub ycalc_other: Vec<f64>,
    /// ycalc_core − ycalc_other.
    pub difference: Vec<f64>,
}

/// A refinement compared with another program's (the app's request 49).
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Comparison {
    pub core: Program,
    pub other: Program,
    pub agreement: Vec<Measure>,
    /// Sorted by |Δ/σ|, the largest first.
    pub parameters: Vec<ParameterRow>,
    pub rule: String,
    pub not_compared: Vec<NotCompared>,
    pub manifest: Vec<ManifestItem>,
    /// What bears on the comparison as a whole (conventions, scalings), in words.
    pub notes: Vec<String>,
    pub patterns: Patterns,
}

/// Write a comparison into its run's folder (`COMPARISON_FILE`); the path written.
#[cfg(feature = "serde")]
pub fn write(dir: &Path, c: &Comparison) -> Result<PathBuf, String> {
    let p = dir.join(COMPARISON_FILE);
    std::fs::write(&p, serde_json::to_string(c).map_err(|e| e.to_string())?).map_err(|e| format!("{}: {e}", p.display()))?;
    Ok(p)
}

/// A comparison written by `write`: the comparison file, or the run's folder holding it.
#[cfg(feature = "serde")]
pub fn read(path: &Path) -> Result<Comparison, String> {
    let p = if path.is_dir() { path.join(COMPARISON_FILE) } else { path.to_path_buf() };
    let t = std::fs::read_to_string(&p).map_err(|e| format!("{}: {e}", p.display()))?;
    serde_json::from_str(&t).map_err(|e| format!("{}: {e}", p.display()))
}

fn latin1(p: &Path) -> Result<String, String> {
    Ok(std::fs::read(p).map_err(|e| format!("{}: {e}", p.display()))?.iter().map(|&b| b as char).collect())
}

/// The cycles of a run: those asked for, or the control file's (10 when it gives none).
fn cycles_of(opened: &Opened, cycles: Option<usize>) -> usize {
    cycles.or(opened.control.refinement.cycles).unwrap_or(10).max(1)
}

/// Set FullProf's output flags of `run_fullprof.py` on the "!Ipr … Hkl" line: Prf = 3 (one layout of `.prf`), Pcr = 1
/// (the refined values written back into the `.pcr`), Rpa = 0 (no CIF output, during which FullProf 8.20 stops on
/// some files, `run_fullprof.py`).
fn fullprof_flags(pcr: &str) -> String {
    let mut lines: Vec<String> = pcr.lines().map(str::to_string).collect();
    if let Some(i) = lines.iter().position(|l| l.trim_start().starts_with("!Ipr") && l.contains("Hkl")) {
        let cols: Vec<String> = lines[i].trim_start_matches(['!', ' ']).split_whitespace().map(str::to_string).collect();
        if let Some(next) = lines.get(i + 1) {
            let mut vals: Vec<String> = next.split_whitespace().map(str::to_string).collect();
            if vals.len() == cols.len() {
                for (k, v) in [("Prf", "3"), ("Pcr", "1"), ("Rpa", "0")] {
                    if let Some(j) = cols.iter().position(|c| c == k) {
                        vals[j] = v.to_string();
                    }
                }
                lines[i + 1] = vals.iter().map(|v| format!("{v:>4}")).collect();
            }
        }
    }
    let mut s = lines.join("\n");
    s.push('\n');
    s
}

/// The manifest of an export: the free parameters carried as codewords, those written fixed, and the export's notes
/// (each an approximation, or a refusal where it names a parameter written fixed).
fn manifest(opened: &Opened, e: &pcr_export::Export) -> Result<Vec<ManifestItem>, String> {
    let set = opened.control.refinement.build(&opened.control.model)?;
    let mut m: Vec<ManifestItem> = Vec::new();
    for p in &set.free {
        for &(f, _) in &p.fields {
            let a = crate::refinable::address(f, &opened.control.model);
            if !e.fixed.contains(&a) {
                m.push(ManifestItem { name: a, carried: Carried::Carried, why: "refined by a codeword".into() });
            }
        }
    }
    for n in &e.notes {
        let fixed = e.fixed.iter().find(|a| n.contains(&format!("'{a}'")));
        match fixed {
            Some(a) => m.push(ManifestItem { name: a.clone(), carried: Carried::Refused, why: n.clone() }),
            None => m.push(ManifestItem { name: String::new(), carried: Carried::Approximated, why: n.clone() }),
        }
    }
    for a in &e.fixed {
        if !m.iter().any(|x| &x.name == a) {
            m.push(ManifestItem { name: a.clone(), carried: Carried::Refused, why: "written fixed: FullProf cannot refine it as the control file does".into() });
        }
    }
    if !opened.control.refinement.plan.is_empty() {
        m.push(ManifestItem { name: String::new(), carried: Carried::Approximated, why: "the control file's plan of steps is not carried: both programs refine every free parameter from the first cycle".into() });
    }
    Ok(m)
}

/// Prepare a FullProf run of a control file in `dir` (created if needed): the `.pcr` and the files it names, the
/// data, and the `.pcr` as written kept as `CODE.pcr.in` (FullProf writes its refined values over `CODE.pcr`). The
/// app runs `Run::program` with `Run::args` in `Run::dir`, `FULLPROF` set to fp2k's folder, nothing on its input.
/// `cycles`: those of the control file when None. Refused: what the export refuses.
pub fn fullprof_prepare(opened: &Opened, dir: &Path, cycles: Option<usize>) -> Result<Run, String> {
    let cycles = cycles_of(opened, cycles);
    std::fs::create_dir_all(dir).map_err(|e| format!("{}: {e}", dir.display()))?;
    let mut c = opened.control.clone();
    c.refinement.cycles = Some(cycles);
    // The data file copied into the folder under the code's name, the extension kept (FullProf reads it by Ins).
    let data_name = match (&c.data, &opened.data_path) {
        (Some(_), Some(p)) => {
            let ext = p.extension().map(|e| e.to_string_lossy().into_owned()).unwrap_or_else(|| "dat".into());
            format!("{FULLPROF_CODE}.{ext}")
        }
        _ => return Err("a comparison needs the control file's data".into()),
    };
    if let Some(d) = c.data.as_mut() {
        d.file = data_name.clone();
    }
    let e = pcr_export::export(&c, FULLPROF_CODE, opened.data.as_ref())?;
    let pcr = fullprof_flags(&e.pcr);
    // The data as FullProf reads them: the processed file the export writes, or the data file copied.
    let data = e.files.iter().find(|(n, _)| n.to_lowercase().ends_with(".xys") && n.starts_with(FULLPROF_CODE)).map_or(data_name.clone(), |(n, _)| n.clone());
    for (n, t) in &e.files {
        std::fs::write(dir.join(n), t).map_err(|er| format!("{n}: {er}"))?;
    }
    if data == data_name {
        let src = opened.data_path.as_ref().unwrap();
        std::fs::copy(src, dir.join(&data_name)).map_err(|er| format!("{}: {er}", src.display()))?;
    }
    for name in [format!("{FULLPROF_CODE}.pcr"), format!("{FULLPROF_CODE}.pcr.in")] {
        std::fs::write(dir.join(&name), &pcr).map_err(|er| format!("{name}: {er}"))?;
    }
    let run = Run {
        other: "FullProf".into(),
        program: "fp2k".into(),
        args: vec![FULLPROF_CODE.into(), data.clone()],
        dir: dir.to_path_buf(),
        env: vec![("FULLPROF".into(), "the folder holding fp2k".into())],
        cycles,
        code: FULLPROF_CODE.into(),
        data,
        manifest: manifest(opened, &e)?,
    };
    #[cfg(feature = "serde")]
    run.save()?;
    Ok(run)
}

/// FullProf's printed esd of `f` in the default model's terms (see the module's notes); None if FullProf printed none.
fn fullprof_esd_in_default(f: Field, m: &Model, sig: &dyn Fn(Field) -> Option<f64>) -> Option<f64> {
    let c = m.convention;
    match f {
        Field::Scale(q) if pcr_export::scale_doubled(m, q) => sig(f).map(|e| 2.0 * e),
        Field::Asym(q, k) if !c.berar_baldinozzi && k % 2 == 0 => {
            let (a, b) = (sig(f), sig(Field::Asym(q, k + 1)));
            (a.is_some() || b.is_some()).then(|| a.unwrap_or(0.0) + 3.0 * b.unwrap_or(0.0))
        }
        Field::SL(q) | Field::DL(q) if !c.axial_divergence && m.phases[q].d_l < 0.0002 => {
            let (s, d) = (sig(Field::SL(q)), sig(Field::DL(q)));
            (s.is_some() || d.is_some()).then(|| (s.unwrap_or(0.0) + d.unwrap_or(0.0)) / 2.0)
        }
        _ => sig(f),
    }
}

/// A convention of FullProf's that bears on the esd of `f` (measured in `fullprof_refinements`).
fn fullprof_esd_note(f: Field, m: &Model) -> Option<String> {
    match f {
        Field::Cell(q, _) | Field::GaussSize(q) if m.phases[q].shape == Shape::ThompsonCoxHastings => {
            Some("FullProf prints half the esd of a Thompson–Cox–Hastings phase's cell parameters and I_G (σ_core/σ_FullProf = 2 in five runs, test fullprof_refinements; why is not known)".into())
        }
        Field::Asym(_, k) if !m.convention.berar_baldinozzi && k % 2 == 0 => Some("FullProf's σ₁ + 3σ₂ of P1 + 3·P2 is an upper bound of the esd".into()),
        _ => None,
    }
}

/// The other program's Ycal at the data's points: as printed where its points are the data's, otherwise
/// interpolated linearly (NaN outside its range).
fn at_points(x: &[f64], ox: &[f64], oy: &[f64]) -> Vec<f64> {
    if ox.len() == x.len() && ox.iter().zip(x).all(|(a, b)| (a - b).abs() <= 1e-6 * (1.0 + b.abs())) {
        return oy.to_vec();
    }
    x.iter()
        .map(|&v| {
            let i = ox.partition_point(|&a| a < v);
            match i {
                _ if ox.is_empty() => f64::NAN,
                0 if (ox[0] - v).abs() <= 1e-9 => oy[0],
                0 => f64::NAN,
                i if i == ox.len() => f64::NAN,
                i => {
                    let t = (v - ox[i - 1]) / (ox[i] - ox[i - 1]);
                    oy[i - 1] + t * (oy[i] - oy[i - 1])
                }
            }
        })
        .collect()
}

/// The core's own refinement of a copy of the control file, as the comparison runs it: the same free parameters
/// and cycles, the control file's ε.
struct Ours {
    model: Model,
    refinement: crate::lsq::Refinement,
    set: crate::refinable::ParameterSet,
}

fn refine_ours(opened: &Opened, cycles: usize) -> Result<Ours, String> {
    let mut model = opened.control.model.clone();
    let set = opened.control.refinement.build(&model)?;
    let data = opened.data.as_ref().ok_or("the control file's data are not read")?;
    let eps = opened.control.refinement.eps.unwrap_or(0.1);
    let refinement = refine(&mut model, data, &set, Options { max_cycles: cycles, eps }, |_| {});
    Ok(Ours { model, refinement, set })
}

fn measure(id: &str, name: &str, phase: Option<usize>, core: f64, other: Option<f64>, print: Option<f64>) -> Measure {
    let other = other.filter(|v| v.is_finite());
    Measure { id: id.into(), name: name.into(), phase, core, other, difference: other.map(|o| core - o), print }
}

/// Read a FullProf run prepared by `fullprof_prepare` and run, refine `ours` (a copy) with the same free parameters
/// and cycles, and compare the two.
pub fn fullprof_read(dir: &Path, run: &Run, ours: &Opened) -> Result<Comparison, String> {
    let code = &run.code;
    let out = read_out(&latin1(&dir.join(format!("{code}.out")))?).map_err(|e| e.to_string())?;
    // FullProf writes the refined values into CODE.pcr (Pcr = 1), or into CODE.new when it stops unconverged on some
    // files (YVO3, L052): the newer of the two.
    let pcr_path = [format!("{code}.new"), format!("{code}.pcr")]
        .into_iter()
        .map(|n| dir.join(n))
        .filter(|p| p.exists())
        .max_by_key(|p| p.metadata().and_then(|m| m.modified()).ok())
        .ok_or(format!("{code}.pcr: not found in {}", dir.display()))?;
    let pcr = crate::pcr::read(&latin1(&pcr_path)?).map_err(|e| format!("{}: {e}", pcr_path.display()))?;
    let fp = Model::from_pcr_with_data(&pcr, dir, code, Some(&dir.join(&run.data)))?.model;
    let prf = std::fs::metadata(dir.join(format!("{code}.prf"))).ok().and_then(|_| latin1(&dir.join(format!("{code}.prf"))).ok()).and_then(|t| read_prf(&t).ok());
    let o = refine_ours(ours, run.cycles)?;
    compare_fullprof(&out, &fp, prf.as_ref(), &o, ours, run)
}

fn compare_fullprof(out: &Out, fp: &Model, prf: Option<&crate::formats::fullprof_outputs::Prf>, o: &Ours, opened: &Opened, run: &Run) -> Result<Comparison, String> {
    let m = &o.model;
    let r = &o.refinement;
    let data = opened.data.as_ref().ok_or("the control file's data are not read")?;
    let mut notes: Vec<String> = Vec::new();
    // Agreement.
    let s = &r.stats;
    let a = out.all_points.clone().unwrap_or_default();
    // FullProf prints its R factors to three figures, χ² (conventional) to three decimals, d to four.
    let fig3 = |v: f64| Some(if v >= 100.0 { 0.5 } else if v >= 10.0 { 0.05 } else if v >= 1.0 { 0.005 } else { 0.0005 });
    let mut agreement = vec![
        measure("chi2", "χ²", None, s.chi2, Some(a.chi2_conventional), Some(5e-4)),
        measure("rp", "Rp (%)", None, s.rp, Some(a.rp), fig3(a.rp)),
        measure("rwp", "Rwp (%)", None, s.rwp, Some(a.rwp), fig3(a.rwp)),
        measure("rexp", "Rexp (%)", None, s.rexp, Some(a.rexp), fig3(a.rexp)),
        measure("rp_background_free", "Rp, background-free (%)", None, s.rp_background_free, Some(a.rp_conventional), fig3(a.rp_conventional)),
        measure("rwp_background_free", "Rwp, background-free (%)", None, s.rwp_background_free, Some(a.rwp_conventional), fig3(a.rwp_conventional)),
        measure("rexp_background_free", "Rexp, background-free (%)", None, s.rexp_background_free, Some(a.rexp_conventional), fig3(a.rexp_conventional)),
        measure("durbin_watson", "Durbin–Watson d", None, s.durbin_watson, Some(a.durbin_watson), Some(5e-5)),
    ];
    for (k, &(rb, rf)) in s.phases.iter().enumerate() {
        let (orb, orf) = out.phases.get(k).copied().map_or((None, None), |(x, y)| (Some(x), Some(y)));
        agreement.push(measure("r_bragg", &format!("R_Bragg, phase {} (%)", k + 1), Some(k), rb, orb, orb.and_then(fig3)));
        agreement.push(measure("r_f", &format!("R_F, phase {} (%)", k + 1), Some(k), rf, orf, orf.and_then(fig3)));
    }
    if s.n - s.p + s.c != a.dof && a.dof > 0 {
        notes.push(format!("the points less the parameters differ: N − P + C = {} here, {} in FullProf", s.n - s.p + s.c, a.dof));
    }
    // FullProf scales its esd's by the χ² over the points with Bragg contributions.
    if let Some(b) = &out.bragg_points {
        if b.chi2_conventional > 0.0 && (a.chi2_conventional - b.chi2_conventional).abs() > 5e-4 {
            notes.push(format!(
                "FullProf scales its esd's by χ² over the points with Bragg contributions ({}), the core by χ² over all the points used ({}): FullProf's esd's are smaller by √ of their ratio, {:.4} (manual eqs. 3.7–3.8)",
                b.chi2_conventional,
                a.chi2_conventional,
                (b.chi2_conventional / a.chi2_conventional).sqrt()
            ));
        }
    }
    // FullProf's parameters as fields of its model (the refined .pcr imported), names shared by atoms of one label
    // taken in order of their numbers.
    let mut seen: std::collections::BTreeMap<String, usize> = std::collections::BTreeMap::new();
    let mut fp_fields: Vec<(Field, &crate::formats::fullprof_outputs::OutParameter)> = Vec::new();
    let mut not_compared: Vec<NotCompared> = Vec::new();
    for p in &out.parameters {
        let fs = crate::output_names::fullprof_fields(&p.name, fp);
        let k = seen.entry(p.name.clone()).or_insert(0);
        match fs.get(*k).or(fs.first()) {
            Some(&f) => fp_fields.push((f, p)),
            None => not_compared.push(NotCompared { name: p.name.clone(), refined_by: "FullProf".into(), why: "the core has no field for this FullProf name".into() }),
        }
        *k += 1;
    }
    let sig_of = |f: Field| -> Option<f64> { fp_fields.iter().find(|x| x.0 == f).map(|x| x.1.esd) };
    let mut parameters: Vec<ParameterRow> = Vec::new();
    let mut matched: Vec<Field> = Vec::new();
    for (j, free) in o.set.free.iter().enumerate() {
        // The value the parameter is named after (a parameter may move others with it).
        let f = free.fields.iter().map(|x| x.0).find(|&f| crate::refinable::address(f, m) == free.name).unwrap_or(free.fields[0].0);
        let address = crate::refinable::address(f, m);
        let Some(&(_, fpp)) = fp_fields.iter().find(|x| x.0 == f) else {
            let why = run.manifest.iter().find(|x| x.name == address && x.carried == Carried::Refused).map_or("FullProf lists no refined parameter for it".to_string(), |x| x.why.clone());
            not_compared.push(NotCompared { name: address, refined_by: "this app".into(), why });
            continue;
        };
        matched.push(f);
        let Some(printed) = fullprof_esd_in_default(f, fp, &sig_of) else { continue };
        // FullProf prints σ·a, a the codeword's multiplier (imported as the parameter's damping).
        let other_esd = printed / free.damping.abs().max(f64::MIN_POSITIVE);
        let core_value = m.get(f);
        let core_esd = r.params.get(j).map_or(f64::NAN, |p| p.esd);
        let other_value = fp.get(f);
        let delta_sigma = if other_esd > 0.0 { (core_value - other_value) / other_esd } else if core_value == other_value { 0.0 } else { f64::INFINITY };
        parameters.push(ParameterRow {
            address,
            other_name: fpp.name.clone(),
            core_value,
            core_esd,
            other_value,
            other_esd,
            delta_sigma,
            esd_ratio: if other_esd > 0.0 { core_esd / other_esd } else { f64::NAN },
            agrees: delta_sigma.abs() <= 1.0,
            note: fullprof_esd_note(f, m),
        });
    }
    for (f, p) in &fp_fields {
        if !matched.contains(f) {
            not_compared.push(NotCompared { name: p.name.clone(), refined_by: "FullProf".into(), why: format!("the control file does not refine {}", crate::refinable::address(*f, m)) });
        }
    }
    parameters.sort_by(|a, b| b.delta_sigma.abs().total_cmp(&a.delta_sigma.abs()));
    // The patterns.
    let ycalc_core = m.calculate(&data.x).ycalc;
    let ycalc_other = match prf {
        Some(p) => {
            // A .prf's Yobs and Ycal may be the data's divided by a power of ten.
            let (sd, sp): (f64, f64) = (data.y.iter().sum(), p.yobs.iter().sum());
            let factor = if sp > 0.0 && sd > 0.0 { 10f64.powi((sd / sp).log10().round() as i32) } else { 1.0 };
            let y: Vec<f64> = p.ycalc.iter().map(|v| v * factor).collect();
            at_points(&data.x, &p.x, &y)
        }
        None => {
            notes.push("FullProf wrote no .prf: its calculated pattern is not compared".into());
            vec![f64::NAN; data.x.len()]
        }
    };
    let difference = ycalc_core.iter().zip(&ycalc_other).map(|(a, b)| a - b).collect();
    let patterns = Patterns { x: data.x.clone(), yobs: data.y.clone(), ycalc_core, ycalc_other, difference };
    Ok(Comparison {
        core: Program { name: "diff-ref-struct".into(), version: format!("core {}", env!("CARGO_PKG_VERSION")), cycles: r.cycles, converged: r.converged },
        other: Program { name: "FullProf".into(), version: out.banner.clone(), cycles: out.cycles, converged: out.converged },
        agreement,
        parameters,
        rule: RULE.into(),
        not_compared,
        manifest: run.manifest.clone(),
        notes,
        patterns,
    })
}
