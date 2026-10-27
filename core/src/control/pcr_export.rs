//! A `.pcr` file for FullProf from a control (L066): the single-pattern ("COMM") format that the core's `.pcr` reader
//! and FullProf read, with the files it names (resolution, reflection lists, data).
//!
//! The model is put back into FullProf's rules where the default model departs from them, by the exact inverses of
//! the import's translations (`Model::translate_from_fullprof`): the X-ray scale halved, P1 − 3·P2 and P3 − 3·P4,
//! the axial pair written so that FullProf reads the paper's S/L and H/L, the exponential preferred orientation's sign,
//! and the excluded regions as the imported `.pcr` wrote them while they still leave out the same points (otherwise
//! each run of points half a step inside its ends, a step beyond the data at its ends). The angular range is the
//! imported `.pcr`'s (`FullProfRecord`) unless the data are written processed. The refinement is
//! written as FullProf's codewords, one parameter number per free parameter: C = sign(a)·(10p + |a|).
//!
//! What FullProf cannot express is reported in `Export::notes`, and the rest is written:
//! - a value that moves with two parameters (codewords have one), a multiplier of 10 or more;
//! - the user's limits (the classic format has none);
//! - the taper of the default peak window (FullProf cuts each peak at ±Wdt·H) and the default's reflections;
//! - angle restraints at their stated weight (FullProf weighs them at π/180 of it, L058);
//! - data operations: the data are then written as processed (x, y, σ; Ins = 10), as is data read by a reader that has
//!   no Ins;
//! - a space group with no symbol FullProf reads (its own table's, or one of the 530 settings').
//!
//! FullProf holds its numbers in single precision; the file writes every value in full.

use super::syntax::float;
use super::Control;
use crate::data::Pattern;
use crate::jacobian::Field;
use crate::powder::{BackgroundModel, InstrumentalWidths, Intensities, Model, Shape};
use crate::restraints::{Kind, Operator};
use std::collections::BTreeMap;

/// The files of an export: the `.pcr` and those it names, and what FullProf cannot express.
#[derive(Clone, Debug, Default)]
pub struct Export {
    pub pcr: String,
    /// (file name, contents) to write next to the `.pcr`.
    pub files: Vec<(String, String)>,
    pub notes: Vec<String>,
    /// The addresses of the free parameters written fixed, which FullProf cannot refine as the control does (each
    /// also has its note): the manifest of a comparison run (the app's request 49).
    pub fixed: Vec<String>,
}

/// The model in FullProf's rules, and the rows (field → column → coefficient) of the refinement in FullProf's values.
type Rows = BTreeMap<Field, BTreeMap<usize, f64>>;

pub(crate) fn scale_doubled(m: &Model, p: usize) -> bool {
    !m.convention.lorentz_polarisation && !m.neutrons && !matches!(m.ilo, -1 | 2 | 3) && matches!(m.phases[p].intensities, Intensities::Structure(..))
}

/// The inverse of `translate_from_fullprof`, on the values and on the refinement's rows.
fn to_fullprof(m: &Model, rows: &mut Rows, x: Option<&[f64]>, record: &super::FullProfRecord, notes: &mut Vec<String>) -> Model {
    let mut f = m.clone();
    let c = m.convention;
    let scale = |rows: &mut Rows, field: Field, k: f64| {
        if let Some(r) = rows.get_mut(&field) {
            for v in r.values_mut() {
                *v *= k;
            }
        }
    };
    // r(a) − k·r(b).
    let minus = |rows: &mut Rows, a: Field, b: Field, k: f64| {
        let rb = rows.get(&b).cloned().unwrap_or_default();
        let ra = rows.entry(a).or_default();
        for (col, v) in rb {
            *ra.entry(col).or_insert(0.0) -= k * v;
        }
        ra.retain(|_, v| *v != 0.0);
        if ra.is_empty() {
            rows.remove(&a);
        }
    };
    for p in 0..m.phases.len() {
        if scale_doubled(m, p) {
            f.phases[p].scale /= 2.0;
            scale(rows, Field::Scale(p), 0.5);
        }
        if !c.berar_baldinozzi {
            let a = &mut f.phases[p].asym;
            a[0] -= 3.0 * a[1];
            a[2] -= 3.0 * a[3];
            minus(rows, Field::Asym(p, 0), Field::Asym(p, 1), 3.0);
            minus(rows, Field::Asym(p, 2), Field::Asym(p, 3), 3.0);
        }
        if !c.axial_divergence {
            let (s, h) = (m.phases[p].s_l, m.phases[p].d_l);
            if (s != 0.0 || h != 0.0) && h < 0.0002 {
                // FullProf reads a pair with D_L < 0.0002 as (H/L + S/L, H/L − S/L).
                f.phases[p].s_l = s + h;
                f.phases[p].d_l = h - s;
                let (rs, rh) = (rows.get(&Field::SL(p)).cloned().unwrap_or_default(), rows.get(&Field::DL(p)).cloned().unwrap_or_default());
                let mut ns = rs.clone();
                let mut nd = rh.clone();
                for (&col, &v) in &rh {
                    *ns.entry(col).or_insert(0.0) += v;
                }
                for (&col, &v) in &rs {
                    *nd.entry(col).or_insert(0.0) -= v;
                }
                ns.retain(|_, v| *v != 0.0);
                nd.retain(|_, v| *v != 0.0);
                for (field, r) in [(Field::SL(p), ns), (Field::DL(p), nd)] {
                    if r.is_empty() {
                        rows.remove(&field);
                    } else {
                        rows.insert(field, r);
                    }
                }
            }
        }
        if !c.preferred_orientation && !m.march {
            f.phases[p].pref[0] = -m.phases[p].pref[0];
            scale(rows, Field::Pref(p, 0), -1.0);
        }
    }
    if !c.excluded_points {
        match x {
            Some(x) if x.len() >= 2 => {
                // The points the default leaves out (closed regions).
                let out: Vec<bool> = x.iter().map(|&v| m.excluded.iter().any(|&(a, b)| v >= a && v <= b)).collect();
                // The regions as the imported .pcr wrote them, while they still leave out these points.
                let mut written = m.clone();
                written.excluded = record.excluded.clone();
                if !record.excluded.is_empty() && written.fullprof_excluded(x) == out {
                    f.excluded = record.excluded.clone();
                } else {
                    // Regions changed or new: FullProf leaves out a point whose step reaches into a region (L033), so
                    // each run of points becomes a region half a step inside its first and last points; at an end of
                    // the data the region goes a step beyond it, never clamped to the data.
                    let mut regions = Vec::new();
                    let mut i = 0;
                    while i < x.len() {
                        if out[i] {
                            let first = i;
                            while i + 1 < x.len() && out[i + 1] {
                                i += 1;
                            }
                            let step = if first + 1 < x.len() { x[first + 1] - x[first] } else { x[first] - x[first - 1] };
                            let back = if i > 0 { x[i] - x[i - 1] } else { step };
                            let a = if first == 0 { x[0] - step } else { x[first] + step / 2.0 };
                            let b = if i + 1 == x.len() { x[i] + back } else { x[i] - back / 2.0 };
                            regions.push((a, b.max(a)));
                        }
                        i += 1;
                    }
                    if !regions.is_empty() {
                        notes.push("excluded regions changed since the import are written half a step inside their first and last points, as FullProf's rule needs (L033)".into());
                    }
                    f.excluded = regions;
                }
            }
            _ if !m.excluded.is_empty() => notes.push("the excluded regions are written as they are: without the data's points they cannot be narrowed, and FullProf leaves out also the points a step from them".into()),
            _ => {}
        }
    }
    if !c.window_edge {
        notes.push("the default peak window tapers each peak to zero beyond ±Wdt·H; FullProf cuts it there (the pattern differs by the tails, L058)".into());
    }
    if !c.background_edge && matches!(m.background, BackgroundModel::Points(_)) {
        notes.push("FullProf reads the data point at or just below the first background point on the line through the first two; the default holds the first point's value there (L098)".into());
    }
    // FullProf's background is 0 below 0 and above the data's largest count (`Convention::background_bounds`, L353).
    if let (false, Some(x), BackgroundModel::Polynomial { .. }) = (c.background_bounds, x, &m.background) {
        let ceiling = m.background_ceiling;
        if m.background_unbounded(x).iter().any(|&b| b < 0.0 || ceiling.is_some_and(|top| b > top)) {
            notes.push("the polynomial background falls below 0 (or above the data's largest count) at some points, where FullProf sets it to 0 (L353)".into());
        }
    }
    if !c.reflections {
        notes.push("the default's reflection range (to the data's end) is FullProf's own rule there: FullProf adds reflections up to 5° past the data".into());
    }
    if !c.angle_restraints && m.restraints.iter().any(|r| matches!(r.kind, Kind::Angle { .. })) {
        notes.push("angle restraints are written with their σ; FullProf weighs them at π/180 of the default's weight (L058)".into());
    }
    f.convention = crate::powder::Convention::FULLPROF;
    f
}

fn npr(s: Shape) -> i32 {
    match s {
        Shape::Gaussian => 0,
        Shape::Lorentzian => 1,
        Shape::PseudoVoigt => 5,
        Shape::ThompsonCoxHastings => 7,
    }
}

/// FullProf's ITnum for an operator's rotation (manual p. 123), the frame chosen by the group: None if the rotation is
/// not one of FullProf's table.
fn itnum(op: &Operator, hexagonal: bool) -> Option<i32> {
    let range: Vec<i32> = if hexagonal { (25..=36).collect() } else { (1..=24).collect() };
    range.into_iter().flat_map(|n| [n, -n]).find(|&n| Operator::itnum(n, op.t).is_some_and(|o| o.rot == op.rot))
}

fn line(out: &mut String, s: &str) {
    out.push_str(s);
    out.push('\n');
}

fn nums(v: &[f64]) -> String {
    v.iter().map(|&x| format!("{:>12}", float(x))).collect::<Vec<_>>().join(" ")
}

/// The end point FullProf adds to a spline (L081): the last point, at the data's last x (or beyond it) with the value
/// of the point before it.
fn spline_end(p: &[(f64, f64)], data: Option<&Pattern>) -> Option<usize> {
    let n = p.len();
    let end = data.and_then(|d| d.x.last().copied())?;
    (n >= 2 && p[n - 1].0 >= end && p[n - 1].1 == p[n - 2].1).then_some(n - 1)
}

/// Write a control as a `.pcr` named `stem` (the files it names start with the stem), with the data's points when
/// known (the excluded regions and the angular range need them).
pub fn export(c: &Control, stem: &str, data: Option<&Pattern>) -> Result<Export, String> {
    let m = &c.model;
    let mut e = Export::default();
    let set = c.refinement.build(m)?;
    let mut rows: Rows = BTreeMap::new();
    for (col, p) in set.free.iter().enumerate() {
        for &(f, k) in &p.fields {
            rows.entry(f).or_default().insert(col, k);
        }
    }
    // A Pawley fit's intensities (P36): FullProf has none (its manual names Pawley's method only, p. 14); the list is
    // written with the values reached, refined by no codeword.
    let pawley = rows.keys().filter(|f| matches!(f, Field::Intensity(..))).count();
    if pawley > 0 {
        e.fixed.extend(rows.keys().filter(|f| matches!(f, Field::Intensity(..))).map(|f| crate::refinable::address(*f, m)));
        rows.retain(|f, _| !matches!(f, Field::Intensity(..)));
        e.notes.push(format!("{pawley} reflection intensities refined by a Pawley fit: FullProf has no Pawley fit, so they are written at the values reached, fixed"));
    }
    let x = data.map(|d| d.x.as_slice());
    let fm = to_fullprof(m, &mut rows, x, &c.fullprof, &mut e.notes);
    // Codewords. The fields that a codeword can carry, by parameter (column).
    let mut coded: Vec<(Field, usize, f64)> = Vec::new();
    for (f, r) in &rows {
        let name = crate::refinable::address(*f, m);
        if r.len() > 1 {
            e.notes.push(format!("'{name}' moves with {} parameters in FullProf's values; a codeword has one: left fixed", r.len()));
            e.fixed.push(name);
            continue;
        }
        let (&col, &a) = r.iter().next().unwrap();
        if a.abs() >= 10.0 {
            e.notes.push(format!("'{name}' moves as {a} × its parameter; a codeword's multiplier is below 10: left fixed"));
            e.fixed.push(name);
            continue;
        }
        coded.push((*f, col, a));
    }
    // FullProf numbers the parameters of the codewords 1…N with no hole ("Hole in the matrix", L071), and the ratio of
    // `VARY ratio` itself: the columns that carry a codeword are numbered in order, the ratio's left out.
    let mut used: Vec<usize> = coded.iter().filter(|c| c.0 != Field::Ratio).map(|c| c.1).collect();
    used.sort_unstable();
    used.dedup();
    let number: BTreeMap<usize, usize> = used.iter().enumerate().map(|(k, &col)| (col, k + 1)).collect();
    let mut codes: BTreeMap<Field, f64> = BTreeMap::new();
    // A parameter's own damping (W12) is FullProf's multiplier of a lone codeword (manual p. 13); on a shared one the
    // multipliers are the link's coefficients, and the damping cannot be written.
    let own: Vec<f64> = set.free.iter().map(|p| c.refinement.steering.damping.iter().find(|d| d.0 == p.name).map_or(1.0, |d| d.1)).collect();
    for (col, p) in set.free.iter().enumerate() {
        if own[col] != 1.0 && coded.iter().filter(|x| x.1 == col).count() != 1 {
            e.notes.push(format!("the damping {} of '{}' is not written: its codeword is shared, and FullProf's multipliers are then the link's coefficients", own[col], p.name));
        }
    }
    for &(f, col, a) in &coded {
        // A lone value moves as its parameter whatever the coefficient (f = c·q refines f itself): FullProf's
        // multiplier is then the damping alone, with the coefficient's sign.
        let a = if coded.iter().filter(|x| x.1 == col).count() == 1 { a.signum() * own[col] } else { a };
        if f == Field::Ratio {
            if coded.iter().any(|c| c.1 == col && c.0 != Field::Ratio) {
                e.notes.push("the ratio λ2/λ1 moves with other values; FullProf's VARY ratio refines it alone: left fixed".into());
                e.fixed.push(crate::refinable::address(Field::Ratio, m));
            } else {
                codes.insert(f, 1.0);
            }
            continue;
        }
        codes.insert(f, a.signum() * (10.0 * number[&col] as f64 + a.abs()));
    }
    let declared = number.len();
    if !c.refinement.limits.is_empty() {
        e.notes.push(format!("{} limits of the user are not written (the .pcr has none)", c.refinement.limits.len()));
    }
    // The steering FullProf has no place for (W12).
    let rf = &c.refinement;
    let d = crate::steering::Steering::default();
    for (missing, what) in [
        (rf.steering.lambda != d.lambda || rf.steering.marquardt != d.marquardt, "λ's start and the Marquardt switch"),
        (rf.steering.guards != d.guards, "the divergence guards"),
        (!rf.steering.max_shift.is_empty(), "the maximum shifts"),
        (!rf.behaviour.is_empty(), "the behaviour at the bounds"),
        (!rf.plan.is_empty(), "the plan's steps"),
    ] {
        if missing {
            e.notes.push(format!("{what} not written (the .pcr has no place for them)"));
        }
    }
    let code = |f: Field| -> String { format!("{:>10}", float(codes.get(&f).copied().unwrap_or(0.0))) };
    // FullProf's Ste: a lone rebin that keeps every n-th point is written as Ste = n, the data as read.
    let ste = match c.data.as_ref().map(|d| d.operations.as_slice()) {
        Some([crate::data_ops::Operation::RebinFactor { factor, mode: crate::data_ops::RebinMode::Keep }]) if *factor > 1 => *factor,
        _ => 1,
    };
    // A TOF pattern's counts as FullProf reads a GSAS file (P05, `powder::read_tof_data`): the data as read.
    let tof_read = m.tof.is_some()
        && c.data.as_ref().is_some_and(|d| {
            use crate::data_ops::{Normalization, Operation};
            d.reader.as_deref() == Some("gsas") && matches!(d.operations.as_slice(), [Operation::Normalize(Normalization::Step(n) | Normalization::Factor(n))] if *n == 100.0 || *n == 1000.0)
        });
    // The data: as read, or processed when there are operations (other than Ste's) or the reader has no Ins.
    let (data_name, ins) = match &c.data {
        Some(d) => match d.reader.as_deref().and_then(crate::formats::by_id).and_then(|f| f.ins) {
            Some(ins) if d.operations.is_empty() || ste > 1 || tof_read => (d.file.clone(), ins),
            _ => {
                let name = format!("{stem}.xys");
                if let Some(p) = data {
                    let axis = if m.tof.is_some() { crate::formats::Axis::TimeOfFlight } else { crate::formats::Axis::TwoTheta };
                    let mut s = crate::formats::Scan::new("processed data", axis, p.clone());
                    s.info.title = Some(format!("{} (processed data)", c.title));
                    let t = crate::formats::write(&s, "xys")?;
                    e.files.push((name.clone(), t));
                    e.notes.push(format!("the data are written processed in {name} (Ins = 10): FullProf does not read their format or does not apply the operations"));
                } else {
                    e.notes.push("the data's format has no FullProf Ins and the data are not at hand: no data file is written".into());
                }
                (name, 10)
            }
        },
        None => (format!("{stem}.dat"), 0),
    };
    // The range: as the imported .pcr gave it, unless the data are written processed; otherwise the data's.
    let processed = ins == 10 && c.data.as_ref().is_some_and(|d| (!d.operations.is_empty() && ste == 1 && !tof_read) || d.reader.as_deref().and_then(crate::formats::by_id).and_then(|f| f.ins).is_none());
    let (thmin, step, mut thmax) = match (c.fullprof.range, data) {
        (Some([a, s, b]), _) if !processed => (a, s, b),
        (_, Some(p)) if p.len() >= 2 => (p.x[0], (p.x[p.len() - 1] - p.x[0]) / (p.len() - 1) as f64, p.x[p.len() - 1]),
        _ => (0.0, 0.0, 0.0),
    };
    // A TOF pattern's TOF-max is the Fourier background's period (`tof::Bank::tmax`).
    if let Some(b) = &m.tof {
        if thmax != b.tmax {
            e.notes.push(format!("TOF-max is written as the bank's {} µs, the Fourier background's period (the data end at {thmax})", b.tmax));
            thmax = b.tmax;
        }
    }
    let mut out = String::new();
    line(&mut out, &format!("COMM {}", c.title));
    line(&mut out, "! Written by drs from a diff-ref-struct control file");
    line(&mut out, &format!("! Files => DAT-file: {data_name},  PCR-file: {stem}"));
    // Background points in a .bac (Nba = 1, b1 = 1): when the control file keeps them in one, or when there are too many
    // for the .pcr's field (FullProf reads a .bac in its first format, a value at every data point: L082). Not when
    // the points are refined (FullProf refines a, b, c, d of a .bac, not its values).
    let points_refined = codes.keys().any(|f| matches!(f, Field::Background(_)));
    let bac = match &fm.background {
        BackgroundModel::Points(p) | BackgroundModel::Spline(p) if !points_refined && (c.background_file.as_ref().is_some_and(|f| f.to_lowercase().ends_with(".bac")) || p.len() > 9999 || fm.background_terms.file.is_some()) => {
            let how = if matches!(fm.background, BackgroundModel::Spline(_)) { crate::background::Interpolation::Spline } else { crate::background::Interpolation::Linear };
            match crate::background::bac_text(p, how, x, &c.title) {
                Ok(t) => {
                    e.files.push((format!("{stem}.bac"), t));
                    true
                }
                Err(er) => {
                    e.notes.push(format!("the background points stay in the .pcr: {er}"));
                    false
                }
            }
        }
        _ => false,
    };
    let (nba, bkpos) = match &fm.background {
        _ if bac => (1, 0.0),
        BackgroundModel::Points(p) => (p.len() as i32, 0.0),
        BackgroundModel::Spline(p) if spline_end(p, data).is_some() => {
            // The point FullProf adds itself (at the last data point, with the last value: L081) is not written.
            let n = p.len() - 1;
            if n <= 5 {
                e.notes.push(format!("a spline through {n} background points: FullProf reads Nba = −{n} as another background; written as linear interpolation"));
                (n as i32, 0.0)
            } else {
                (-(n as i32), 0.0)
            }
        }
        BackgroundModel::Spline(p) => {
            e.notes.push("a spline whose last point is not FullProf's own end point: FullProf adds a point of its own at the last data point with the last value (L081), which changes its curve, most beyond the last point but one".into());
            // FullProf's cubic splines: Nba < 0 with |Nba| > 4 (manual p. 69); −4 and −5 are its Fourier and Chebyshev
            // series, so a spline of five points or fewer cannot be written (L078; Nba = −6 measured a spline in L081).
            if p.len() <= 5 {
                e.notes.push(format!("a spline through {} background points: FullProf reads Nba = −{} as another background; written as linear interpolation", p.len(), p.len()));
                (p.len() as i32, 0.0)
            } else {
                (-(p.len() as i32), 0.0)
            }
        }
        BackgroundModel::Polynomial { origin, .. } => (0, *origin),
        BackgroundModel::Chebyshev { from, to, .. } => {
            // FullProf's spans the data's range, or a TOF pattern's TOF-min to TOF-max (P05, CeO2_PEARL_calc.prf).
            if m.tof.is_some() {
                if (*from, *to) != (thmin, thmax) {
                    e.notes.push(format!("FullProf's Chebyshev background of a TOF pattern spans TOF-min to TOF-max, {thmin}–{thmax}; this one spans {from}–{to}"));
                }
            } else if data.is_none_or(|p| p.x.first() != Some(from) || p.x.last() != Some(to)) {
                e.notes.push(format!("FullProf's Chebyshev background spans the data's range; this one spans {from}–{to}"));
            }
            (-5, 0.0)
        }
        BackgroundModel::Fourier { .. } => (-4, 0.0),
    };
    // Further terms of the background (P23): Debye-like terms with the polynomial of six coefficients (Nba = −1, eq. 3.4;
    // fewer terms written as six, the others' amplitudes 0), inverse powers of 2θ with that of nine (Nba = −3), and
    // profiles from files (MORE = 2, 3, from {data}_n.bac), for a constant-wavelength pattern.
    let bt = &fm.background_terms;
    let poly = match &fm.background {
        BackgroundModel::Polynomial { coefficients, .. } if !bac && fm.tof.is_none() => Some(coefficients.len()),
        _ => None,
    };
    let term_nba = match (bt.debye.is_empty(), bt.inverse.is_empty(), poly) {
        (true, true, _) => None,
        (false, true, Some(n)) if n <= 6 && bt.debye.len() <= 6 => Some(-1),
        (true, false, Some(n)) if n <= 9 => Some(-3),
        _ => {
            e.notes.push("Debye-like background terms and inverse powers of 2θ are FullProf's with its polynomial alone (Nba = −1: six coefficients and six terms; −3: nine coefficients and three powers), one kind at a time, and not in a TOF pattern: not written".into());
            None
        }
    };
    let nba = term_nba.unwrap_or(nba);
    let data_stem = std::path::Path::new(&data_name).file_stem().map_or(stem.to_string(), |s| s.to_string_lossy().to_string());
    let mut profile_files = Vec::new();
    if !bt.profiles.is_empty() {
        if fm.tof.is_some() {
            e.notes.push("background profiles from files are FullProf's for constant-wavelength patterns (MORE on line 14): not written".into());
        } else {
            for (k, p) in bt.profiles.iter().enumerate() {
                match crate::background::bac_text(&p.points, crate::background::Interpolation::Linear, x, &format!("{} (background profile {})", c.title, k + 1)) {
                    Ok(t) => profile_files.push((format!("{data_stem}_{}.bac", k + 1), t)),
                    Err(er) => {
                        e.notes.push(format!("the background profiles are not written: {er}"));
                        profile_files.clear();
                        break;
                    }
                }
            }
        }
    }
    let profiles_written = !profile_files.is_empty();
    e.files.extend(profile_files);
    let instrument = fm.tof.as_ref().and_then(|b| b.instrument.as_ref());
    let res = match &fm.instrumental {
        InstrumentalWidths::None if instrument.is_some() => 5,
        InstrumentalWidths::None => 0,
        InstrumentalWidths::Function { res, .. } => *res,
        InstrumentalWidths::Table(_) => 4,
    };
    // A TOF pattern (P05): Job = −1, Npr 9 (10 with Hodges' thermal terms), Uni = 1.
    let tof = fm.tof.as_ref();
    let job = if tof.is_some() { -1 } else if fm.neutrons { 1 } else { 0 };
    let tof_npr = tof.map(|b| if b.thermal.is_some() { 10 } else { 9 });
    let npr_of = |sh: Shape| tof_npr.unwrap_or_else(|| npr(sh));
    let iwg = match fm.weighting {
        crate::powder::Weighting::Variance => 0,
        crate::powder::Weighting::MaximumLikelihood => 1,
        crate::powder::Weighting::Unit => 2,
    };
    let nli = c.refinement.linear.len();
    // The user's scattering factors (P11) as line 12: Nsc < 0 when one is a table (manual p. 70), and then a form factor
    // given by its coefficients cannot be written (ITY = 0 reads a table there).
    let user = fm.user_scatterers();
    let tables = user.iter().any(|u| !u.form_factor_table.is_empty());
    let mut line12 = Vec::new();
    for u in &user {
        let name = if fm.neutrons { u.species.clone() } else { u.species.to_lowercase() };
        if fm.neutrons {
            line12.push(format!("{name:<6} {} {} 0", float(u.scattering_length.unwrap_or(0.0)), float(u.scattering_length_imaginary)));
        } else if let (Some(k), false) = (u.form_factor, tables) {
            line12.push(format!("{name:<6} {} {} 0", float(u.f_prime), float(u.f_double_prime)));
            line12.push(format!("  {}", nums(&k)));
        } else if !u.form_factor_table.is_empty() {
            line12.push(format!("{name:<6} {} {} 0", float(u.f_prime), float(u.f_double_prime)));
            for (x, f) in &u.form_factor_table {
                line12.push(format!("  {} {}", float(*x), float(*f)));
            }
            line12.push("  -100".to_string());
        } else if u.form_factor.is_some() {
            e.notes.push(format!("the form factor of '{}' given by its coefficients is not written: FullProf reads tables only when another species gives one (Nsc < 0)", u.species));
        } else {
            line12.push(format!("{name:<6} {} {} 2", float(u.f_prime), float(u.f_double_prime)));
        }
    }
    let nsc = user.len() as i32 * if tables { -1 } else { 1 };
    // Nor = 2, 3 (P30): several axes, which every phase then gives; a phase of one March–Dollase axis without G2 gives
    // it as one axis.
    let textured: Vec<&crate::texture::MultiAxial> = fm.phases.iter().filter_map(|p| p.texture.as_ref()).collect();
    let nor = if textured.is_empty() {
        i32::from(fm.march)
    } else if textured.iter().all(|t| t.steps > 0) {
        3
    } else if textured.iter().all(|t| t.steps == 0) {
        2
    } else {
        return Err("phases of a flat sample and of an anvil cell (po_axis_steps) together: FullProf's Nor is one for the pattern".into());
    };
    if nor >= 2 {
        if let Some(p) = fm.phases.iter().position(|p| p.texture.is_none() && !(fm.march && p.pref[1] == 0.0)) {
            return Err(format!("phase {}: FullProf's Nor = {nor} gives every phase March–Dollase axes, and this phase's preferred orientation is not one (exponential, or G2 ≠ 0)", p + 1));
        }
    }
    let axes_lines = |out: &mut String, p: usize| {
        if nor < 2 {
            return;
        }
        let ph = &fm.phases[p];
        let one = crate::texture::MultiAxial { axes: vec![crate::texture::Axis { direction: ph.pref_direction, r: ph.pref[0], fraction: 1.0 }], steps: textured.first().map_or(0, |t| t.steps) };
        let t = ph.texture.as_ref().unwrap_or(&one);
        if nor == 3 {
            line(out, "!  Number of preferred orientation directions and integration steps");
            line(out, &format!("    {}   {}", t.axes.len(), t.steps));
        } else {
            line(out, "!  Number of preferred orientation directions");
            line(out, &format!("    {}", t.axes.len()));
        }
        line(out, "!  Pref-hkl       MD-r     Code   Fraction  Code");
        let n = t.axes.len();
        for (j, a) in t.axes.iter().enumerate() {
            let rc = if ph.texture.is_some() { code(Field::Texture(p, 2 * j)) } else { code(Field::Pref(p, 0)) };
            let fc = if j + 1 < n { code(Field::Texture(p, 2 * j + 1)) } else { format!("{:>10}", "0.000") };
            line(out, &format!("  {} {} {} {} {rc} {} {fc}", float(a.direction[0]), float(a.direction[1]), float(a.direction[2]), float(a.r), float(a.fraction)));
        }
    };
    line(&mut out, "!Job Npr Nph Nba Nex Nsc Nor Dum Iwg Ilo Ias Res Ste Nre Cry Uni Cor Opt Aut");
    let npr0 = fm.phases.first().map_or(tof_npr.unwrap_or(5), |p| npr_of(p.shape));
    line(&mut out, &format!("{job:>4}{npr0:>4}{:>4}{nba:>4}{:>4}{nsc:>4}{:>4}   1{iwg:>4}{:>4}   1{res:>4}{ste:>4}   0   0{:>4}   0   0   0", fm.phases.len(), fm.excluded.len(), nor, fm.ilo, i32::from(tof.is_some())));
    if res != 0 {
        let name = format!("{stem}.irf");
        let mut t = String::from("Instrumental resolution written by drs\n");
        match &fm.instrumental {
            InstrumentalWidths::Function { uvwxyz, .. } => {
                t.push_str(&format!("{}\n{}\n", nums(uvwxyz), nums(uvwxyz)));
            }
            InstrumentalWidths::Table(rows) => {
                t.push_str(&format!("{}\n", rows.len()));
                for r in rows {
                    t.push_str(&format!("{}\n", nums(r)));
                }
            }
            InstrumentalWidths::None => {
                // A TOF resolution file's keywords (Res = 5; FullProf's News of 18 April 2011), as `tof::Instrument::read`
                // reads them, the bank's angle with them (TWOTH replaces the .pcr's).
                if let (Some(i), Some(b)) = (instrument, tof) {
                    t.push_str(&format!("ZD2TOF {}\nZD2TOT {}\nTWOTH {}\n", nums(&[i.zero, i.dtt[0], i.dtt[1]]), nums(&i.thermal), float(b.two_theta)));
                    t.push_str(&format!("SIGMA {}\nGAMMA {}\n", nums(&[i.sigma[2], i.sigma[1], i.sigma[0]]), nums(&[i.gamma[2], i.gamma[1], i.gamma[0]])));
                    t.push_str(&format!("ALFBE {}\nALFBT {}\n", nums(&[i.alpha[0], i.beta[0], i.alpha[1], i.beta[1]]), nums(&i.thermal_decay)));
                }
            }
        }
        // The name alone on its line: FullProf takes the whole line as the file's name (L069).
        line(&mut out, "!  Resolution file for Pattern#   1");
        line(&mut out, &name);
        e.files.push((name, t));
    }
    line(&mut out, "!");
    line(&mut out, "!Ipr Ppl Ioc Mat Pcr Ls1 Ls2 Ls3 NLI Prf Ins Rpa Sym Hkl Fou Sho Ana");
    line(&mut out, &format!("   0   0   1   0   1   0   4   0{nli:>4}   3{ins:>4}   0   0   0   0   0   0"));
    line(&mut out, "!");
    // Ratio < 0 tells FullProf that each phase gives the second wavelength's U, V, W (line 32).
    let own_uvw2 = fm.phases.iter().any(|p| p.uvw2.is_some());
    if let Some(b) = tof {
        line(&mut out, "!  Bkpos        Wdt    Iabscor for Pattern#  1");
        line(&mut out, &format!(" {} {} {}", float(bkpos), float(fm.wdt), b.absorption.iabscor()));
    } else {
        line(&mut out, "! Lambda1  Lambda2    Ratio    Bkpos    Wdt    Cthm     muR   AsyLim   Rpolarz  2nd-muR -> Patt# 1");
        let ratio = if own_uvw2 { -fm.ratio.abs() } else { fm.ratio };
        line(&mut out, &format!(" {} {} {} {} {} {} {} {} {} 0.0", float(fm.lambda[0]), float(fm.lambda[1]), float(ratio), float(bkpos), float(fm.wdt), float(fm.cthm), float(fm.mur), float(fm.asym_limit), float(fm.rpolarz)));
    }
    // Commands straight after the wavelength line, as FullProf's own files have them (L069).
    if codes.contains_key(&Field::Ratio) {
        line(&mut out, "VARY ratio");
    }
    line(&mut out, "!");
    let r = c.refinement.steering.class_damping;
    let ncy = format!("{:>3} {} {} {} {} {} {} {} {}", c.refinement.cycles.unwrap_or(10), float(c.refinement.eps.unwrap_or(0.1)), float(r[0]), float(r[1]), float(r[2]), float(r[3]), float(thmin), float(step), float(thmax));
    if tof.is_some() {
        // A TOF pattern's line 9 ends at TOF-max (p. 76).
        line(&mut out, "!NCY  Eps  R_at  R_an  R_pr  R_gl    TOF-min      <Step>       TOF-max");
        line(&mut out, &ncy);
    } else {
        line(&mut out, "!NCY  Eps  R_at  R_an  R_pr  R_gl     Thmin       Step       Thmax    PSD    Sent0");
        line(&mut out, &format!("{ncy} 0.000 0.000"));
    }
    line(&mut out, "!");
    if let (false, BackgroundModel::Points(p) | BackgroundModel::Spline(p)) = (bac, &fm.background) {
        let p = match &fm.background {
            BackgroundModel::Spline(q) if spline_end(q, data).is_some() => &q[..q.len() - 1],
            _ => &p[..],
        };
        line(&mut out, "!2Theta/TOF/E(Kev)   Background  for Pattern#  1");
        for (k, (a, b)) in p.iter().enumerate() {
            line(&mut out, &format!("{:>14} {:>14} {}", float(*a), float(*b), code(Field::Background(k))));
        }
        line(&mut out, "!");
    }
    if !fm.excluded.is_empty() {
        line(&mut out, "!       Excluded regions (LowT  HighT) for Pattern#  1");
        for (a, b) in &fm.excluded {
            line(&mut out, &format!("{:>14} {:>14}", float(*a), float(*b)));
        }
        line(&mut out, "!");
    }
    if !line12.is_empty() {
        line(&mut out, "! Additional scattering factors for Pattern#  1");
        for l in &line12 {
            line(&mut out, l);
        }
        line(&mut out, "!");
    }
    line(&mut out, &format!("{declared:>6}    !Number of refined parameters"));
    line(&mut out, "!");
    if let Some(b) = tof {
        // Line 14 (p. 82, with Dtt_1overd as FullProf's newer files write it) and line 16 for Npr = 10 (p. 84).
        line(&mut out, "!    Zero       Code      Dtt1      Code       Dtt2     Code  Dtt_1overd    Code  2ThetaBank  -> Patt#  1");
        line(&mut out, &format!(" {} {} {} {} {} {} {} {} {}", float(fm.zero), code(Field::Zero), float(b.dtt1), code(Field::Tof(0)), float(b.dtt2), code(Field::Tof(1)), float(b.dtt_1overd), code(Field::Tof(2)), float(b.two_theta)));
        if let Some(t) = &b.thermal {
            line(&mut out, "!   Zerot   Dtt1t   Dtt2t   x-cross   Width");
            line(&mut out, &format!(" {}", nums(&[t.zero, t.dtt1, t.dtt2, t.cross, t.width])));
            line(&mut out, &format!(" {}", (3..8).map(|k| code(Field::Tof(k))).collect::<Vec<_>>().join(" ")));
        }
    } else {
        line(&mut out, "!  Zero    Code    SyCos    Code   SySin    Code  Lambda     Code MORE ->Patt# 1");
        // Lambda and its codeword when the wavelength is refined (P18): FullProf writes λ₁ there (manual p. 83). More = 1 and
        // line 15 for the surface micro-absorption of a Bragg–Brentano X-ray pattern (P29, p. 84).
        let lambda = if codes.contains_key(&Field::Wavelength) { format!("{} {}", float(fm.lambda[0]), code(Field::Wavelength)) } else { "0.0 0.0".to_string() };
        let micro = !fm.neutrons && fm.ilo == 0 && ((0..3).any(|k| codes.contains_key(&Field::Micro(k))) || fm.microabsorption.iter().any(|&v| v != 0.0));
        // MORE = 2 with background profiles, 3 with both (P23).
        let more = match (micro, profiles_written) {
            (false, false) => 0,
            (true, false) => 1,
            (false, true) => 2,
            (true, true) => 3,
        };
        line(&mut out, &format!(" {} {} {} {} {} {} {lambda} {more}", float(fm.zero), code(Field::Zero), float(fm.sycos), code(Field::Sycos), float(fm.sysin), code(Field::Sysin)));
        if micro {
            line(&mut out, "!  P0    Cod_P0    Cp    Cod_Cp    Tau    Cod_Tau");
            let m = fm.microabsorption;
            line(&mut out, &format!(" {} {} {} {} {} {}", float(m[0]), code(Field::Micro(0)), float(m[1]), code(Field::Micro(1)), float(m[2]), code(Field::Micro(2))));
        } else if fm.microabsorption.iter().any(|&v| v != 0.0) {
            e.notes.push("the surface micro-absorption is FullProf's for Bragg–Brentano X-ray patterns only (line 15): not written".into());
        }
        if profiles_written {
            let n = bt.profiles.len();
            let first = 2 * bt.debye.len() + bt.inverse.len();
            line(&mut out, "! Number and coefficients of profiles contributing to background of Pattern#  1");
            line(&mut out, &format!("  {n} {}", nums(&bt.profiles.iter().map(|p| p.coefficient).collect::<Vec<_>>())));
            line(&mut out, &format!("  {}", (0..n).map(|j| code(Field::BackgroundTerm(first + j))).collect::<Vec<_>>().join(" ")));
        }
    }
    if codes.contains_key(&Field::Ratio) {
        e.notes.push("the ratio λ2/λ1 refined is written as FullProf's command VARY ratio".into());
    }
    let coef: Option<&Vec<f64>> = match &fm.background {
        BackgroundModel::Polynomial { coefficients, .. } | BackgroundModel::Chebyshev { coefficients, .. } | BackgroundModel::Fourier { coefficients } => Some(coefficients),
        BackgroundModel::Points(_) | BackgroundModel::Spline(_) => None,
    };
    if bac {
        // Back(2θ) = a·BackRead[(1 + c)2θ + d] + b, as FullProf's .out writes it: b1–b4 of the file (eq. 3.6; the app's
        // request 47) with their codewords, a = 1 and the others 0 without them.
        line(&mut out, "!");
        line(&mut out, "!   Background coefficients/codes  for Pattern#  1");
        match &fm.background_terms.file {
            Some(f) => {
                let first = fm.background_terms.len() - 4;
                line(&mut out, &format!("  {} 0.0 0.0", nums(&f.values())));
                line(&mut out, &format!("  {} 0.0 0.0", (0..4).map(|k| code(Field::BackgroundTerm(first + k))).collect::<Vec<_>>().join(" ")));
            }
            None => {
                line(&mut out, "  1.0 0.0 0.0 0.0 0.0 0.0");
                line(&mut out, "  0.0 0.0 0.0 0.0 0.0 0.0");
            }
        }
    } else if let Some(f) = fm.background_terms.file.as_ref().filter(|f| !f.is_identity()) {
        e.notes.push(format!("the background's points are refined: FullProf refines a .bac's b1–b4, not its points, so b1–b4 ({}) are not written", nums(&f.values())));
    }
    if let Some(coef) = coef.filter(|_| !bac) {
        line(&mut out, "!");
        line(&mut out, "!   Background coefficients/codes  for Pattern#  1");
        let padded = |v: &[f64], n: usize| (0..n).map(|i| v.get(i).copied().unwrap_or(0.0)).collect::<Vec<f64>>();
        match term_nba {
            // Nba = −1: the six coefficients, then the six terms' amplitudes and distances (an absent term: A = 0, r = 1).
            Some(-1) => {
                line(&mut out, &nums(&padded(coef, 6)));
                line(&mut out, &(0..6).map(|i| code(Field::Background(i))).collect::<Vec<_>>().join(" "));
                line(&mut out, "!  Additional background coefficients/codes for Pattern#  1");
                let d = &bt.debye;
                line(&mut out, &nums(&padded(&d.iter().map(|t| t.amplitude).collect::<Vec<_>>(), 6)));
                line(&mut out, &(0..6).map(|j| if j < d.len() { code(Field::BackgroundTerm(2 * j)) } else { format!("{:>10}", float(0.0)) }).collect::<Vec<_>>().join(" "));
                line(&mut out, &nums(&(0..6).map(|j| d.get(j).map_or(1.0, |t| t.distance)).collect::<Vec<_>>()));
                line(&mut out, &(0..6).map(|j| if j < d.len() { code(Field::BackgroundTerm(2 * j + 1)) } else { format!("{:>10}", float(0.0)) }).collect::<Vec<_>>().join(" "));
            }
            // Nba = −3: nine coefficients and three inverse powers, as two pairs of lines.
            Some(_) => {
                let v: Vec<f64> = padded(coef, 9).into_iter().chain(padded(&bt.inverse, 3)).collect();
                let f = |i: usize| if i < 9 { Field::Background(i) } else { Field::BackgroundTerm(i - 9) };
                for k in 0..2 {
                    line(&mut out, &nums(&v[6 * k..6 * k + 6]));
                    line(&mut out, &(6 * k..6 * k + 6).map(|i| code(f(i))).collect::<Vec<_>>().join(" "));
                }
            }
            None => {
                for (k, chunk) in coef.chunks(6).enumerate() {
                    line(&mut out, &nums(chunk));
                    line(&mut out, &(0..chunk.len()).map(|i| code(Field::Background(6 * k + i))).collect::<Vec<_>>().join(" "));
                }
            }
        }
    }
    for (p, ph) in fm.phases.iter().enumerate() {
        line(&mut out, "!-------------------------------------------------------------------------------");
        line(&mut out, &format!("!  Data for PHASE number: {:>3}", p + 1));
        line(&mut out, "!-------------------------------------------------------------------------------");
        line(&mut out, &ph.name);
        line(&mut out, "!");
        let (nat, atoms) = match &ph.intensities {
            Intensities::Structure(st, _) => (st.atoms.len(), Some(st)),
            Intensities::List(_) => (0, None),
        };
        let restraints: Vec<&crate::restraints::Restraint> = fm.restraints.iter().filter(|r| r.phase == p).collect();
        let dis = restraints.iter().filter(|r| matches!(r.kind, Kind::Distance { .. })).count();
        let ang = restraints.len() - dis;
        // Symmetry modes (R29): Jbt = 6 with Furth the number of amplitudes (FullProf's note of 29 August 2008).
        let (jbt, irf) = if ph.extract || matches!(ph.intensities, Intensities::List(_)) {
            (2, 2)
        } else if ph.modes.is_some() {
            (6, 0)
        } else if !ph.rigid.is_empty() {
            // Rigid groups (P41): Jbt = 4.
            (4, 0)
        } else {
            (0, 0)
        };
        // Furth: the modes' amplitudes, or a TLS group's 21 values (L367; LINE 43, manual p. 38).
        let tls = ph.rigid.iter().enumerate().find(|(_, g)| g.displacements.is_some());
        if ph.rigid.iter().filter(|g| g.displacements.is_some()).count() > 1 {
            return Err(format!("phase {}: TLS on several rigid groups: FullProf's Furth lines hold one group's", p + 1));
        }
        let furth = ph.modes.as_ref().map_or(if tls.is_some() { 21 } else { 0 }, |m| m.modes.len());
        let more = i32::from(ph.brindley != 1.0);
        line(&mut out, "!Nat Dis Ang Pr1 Pr2 Pr3 Jbt Irf Isy Str Furth       ATZ    Nvk Npr More");
        let d = ph.pref_direction;
        // The quartic form of strain (P25): Str = 1 and the Strain-Model of its Laue class (manual pp. 110, 119–120).
        // FullProf's strain model 8 (Str = 0; `powder::hexagonal_strain`): a trigonal or hexagonal form with S_400,
        // S_004 and S_112 alone, Gaussian, written back as Saa, Scc, Cac with their codewords (L322).
        let hexagonal = ph.quartic_strain.as_ref().and_then(|q| {
            let get = |t: [u8; 3]| q.terms.iter().position(|x| x.0 == t);
            let (a, c, ac) = (get([4, 0, 0])?, get([0, 0, 4])?, get([1, 1, 2])?);
            let others_zero = q.terms.iter().enumerate().all(|(k, t)| k == a || k == c || k == ac || t.1 == 0.0);
            let (s400, s004, s112) = (q.terms[a].1, q.terms[c].1, q.terms[ac].1);
            let (saa, scc) = (s400.max(0.0).sqrt(), s004.max(0.0).sqrt());
            let cac = if saa * scc > 0.0 { s112 / (2.0 * saa * scc) } else { 0.0 };
            let ok = matches!(ph.group.laue, "-3" | "-3m" | "6/m" | "6/mmm") && others_zero && q.lorentzian == 0.0 && s400 >= 0.0 && s004 >= 0.0 && (saa * scc > 0.0 || s112 == 0.0) && crate::microstructure::fullprof_model_of(&q.terms.iter().map(|t| t.0).collect::<Vec<_>>()).is_none();
            ok.then_some(([saa, scc, cac], [a, c, ac]))
        });
        let strain_model = match &ph.quartic_strain {
            Some(_) if hexagonal.is_some() => None,
            Some(q) => {
                let terms: Vec<[u8; 3]> = q.terms.iter().map(|t| t.0).collect();
                let m = crate::microstructure::fullprof_model_of(&terms).ok_or(format!("phase {}: a quartic form of strain with the terms {} has no Strain-Model FullProf reads in a known order (−1: L211)", p + 1, terms.iter().map(|t| crate::microstructure::term_name(*t)).collect::<Vec<_>>().join(", ")))?;
                // 4/m is FullProf's 4, 4/mmm its 5 (the same terms); 4/m's S_310 FullProf has not (L267).
                if q.terms.iter().any(|t| t.0 == [3, 1, 0] && t.1 != 0.0) {
                    return Err(format!("phase {}: S_310 of 4/m is not in FullProf's Strain-Model 4 (manual p. 120): the .pcr cannot hold it (L267)", p + 1));
                }
                Some(if m == 5 && ph.group.laue == "4/m" { 4 } else { m })
            }
            None => None,
        };
        let str_flag = i32::from(strain_model.is_some());
        line(&mut out, &format!("{nat:>4}{dis:>4}{ang:>4} {} {} {}{jbt:>4}{irf:>4}   0{str_flag:>4} {furth:>3} {}   0{:>4}{more:>4}", float(d[0]), float(d[1]), float(d[2]), float(ph.atz), npr_of(ph.shape)));
        if more == 1 {
            line(&mut out, "!Jvi Jdi Hel Sol Mom Ter  Brind   RMua    RMub    RMuc");
            line(&mut out, &format!("   0   0   0   0   0   0 {}  0.0 0.0 0.0", float(ph.brindley)));
        }
        line(&mut out, "!");
        let symbol = ph.group.fullprof_symbol().ok_or(format!("phase {}: its space group has no symbol FullProf reads (FullProf's own symbols, or the 530 settings')", p + 1))?;
        line(&mut out, &format!("{symbol:<25}<--Space group symbol"));
        if let (Some(st), false) = (atoms, ph.rigid.is_empty()) {
            // Lines 25-1 to 25-4 of a rigid-body phase (manual pp. 100, 104; RBdoc pp. 9–11): x, y, z as placed (FullProf's
            // information), P6 θ φ of a group's first atom, a Z-matrix atom's N1 N2 N3; its internal coordinates, and the
            // first atom's X0 Y0 Z0, χ, P16 (a satellite's N1.N2), DEG = 1 and KIND.
            line(&mut out, "!Atom   Typ       x        y        z         B      Occ       P6     THETA     PHI  Spc");
            line(&mut out, "!    r/xc/rho the/yc/phi phi/zc/z   X0       Y0       Z0      CHI    P16:SAT   DEG  KIND");
            for (a, at) in st.atoms.iter().enumerate() {
                let Some((g, gr, k)) = ph.rigid.iter().enumerate().find_map(|(g, gr)| gr.members.iter().position(|m| m.index == a).map(|k| (g, gr, k))) else {
                    return Err(format!("phase {}: atom {} is in no rigid group: a free atom of a rigid-body phase is not written yet", p + 1, at.label));
                };
                let m = &gr.members[k];
                let r = |i: usize| code(Field::Rigid(p, g, i));
                let kind = [crate::rigid_body::Coordinates::Spherical, crate::rigid_body::Coordinates::Cartesian, crate::rigid_body::Coordinates::Cylindrical, crate::rigid_body::Coordinates::ZMatrix].iter().position(|c| *c == gr.coordinates).unwrap_or(0);
                let base = format!("{:<6} {:<6} {} {} {} {} {}", at.label, at.species, float(at.x[0]), float(at.x[1]), float(at.x[2]), float(at.biso), float(at.occ));
                let i0 = 6 + 3 * k;
                if k == 0 {
                    line(&mut out, &format!("{base} {} {} {}   0", float(gr.option), float(gr.theta), float(gr.phi)));
                    line(&mut out, &format!("              0.0 0.0 0.0 {} {} 0.0 {} {}", code(Field::AtomB(p, a)), code(Field::AtomOcc(p, a)), r(3), r(4)));
                    let p16 = gr.satellite.map_or(0.0, |s| s.atom as f64 + s.from as f64 / 100.0);
                    line(&mut out, &format!("      {} {} {} {} {} {} {} {} 1 {kind}", float(m.internal[0]), float(m.internal[1]), float(m.internal[2]), float(gr.origin[0]), float(gr.origin[1]), float(gr.origin[2]), float(gr.chi), float(p16)));
                    line(&mut out, &format!("      {} {} {} {} {} {} {}", r(i0), r(i0 + 1), r(i0 + 2), r(0), r(1), r(2), r(5)));
                } else {
                    line(&mut out, &format!("{base} {} {} {}   0", m.connect[0], m.connect[1], m.connect[2]));
                    line(&mut out, &format!("              0.0 0.0 0.0 {} {}", code(Field::AtomB(p, a)), code(Field::AtomOcc(p, a))));
                    line(&mut out, &format!("      {} {} {}", float(m.internal[0]), float(m.internal[1]), float(m.internal[2])));
                    line(&mut out, &format!("      {} {} {}", r(i0), r(i0 + 1), r(i0 + 2)));
                }
            }
        } else if let Some(st) = atoms {
            line(&mut out, "!Atom   Typ       X        Y        Z     Biso       Occ     In Fin N_t Spc /Codes");
            for (a, at) in st.atoms.iter().enumerate() {
                let nt = if at.beta.is_some() { 2 } else { 0 };
                // A mode phase's atoms are written at the reference positions; FullProf adds the modes.
                let x = ph.modes.as_ref().map_or(at.x, |m| m.reference[a]);
                line(&mut out, &format!("{:<6} {:<6} {} {} {} {} {}   0   0 {nt:>3}    0", at.label, at.species, float(x[0]), float(x[1]), float(x[2]), float(at.biso), float(at.occ)));
                line(&mut out, &format!("              {} {} {} {} {}", code(Field::AtomX(p, a, 0)), code(Field::AtomX(p, a, 1)), code(Field::AtomX(p, a, 2)), code(Field::AtomB(p, a)), code(Field::AtomOcc(p, a))));
                if let Some(b) = at.beta {
                    line(&mut out, &format!("              {}", nums(&b)));
                    line(&mut out, &format!("              {}", (0..6).map(|k| code(Field::AtomBeta(p, a, k))).collect::<Vec<_>>().join(" ")));
                }
            }
        }
        if let (Some(ms), Some(st)) = (&ph.modes, atoms) {
            let nv: usize = ms.modes.iter().map(|m| m.vectors.len()).sum();
            line(&mut out, "! Polarisation Vectors of Symmetry Modes for each atom");
            line(&mut out, &format!("V_MODES {nv:>4}"));
            line(&mut out, "! Nm  Atm     Irrep               Vx        Vy        Vz         Coeff");
            for (k, m) in ms.modes.iter().enumerate() {
                for v in &m.vectors {
                    line(&mut out, &format!("{:>4}  {:<7} {:<12} {} {} {} {}", k + 1, st.atoms[v.atom].label, m.irrep, float(v.v[0]), float(v.v[1]), float(v.v[2]), float(v.coeff)));
                }
            }
            line(&mut out, "! Amplitudes of Symmetry Modes");
            let fst = if ms.fst.is_empty() { String::new() } else { format!(" {}", ms.fst.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(" ")) };
            line(&mut out, &format!("A_MODES {:>4}{fst}", ms.modes.len()));
            if let Some(b) = ms.max_amplitude {
                line(&mut out, &format!("Max_Amplitude {}", float(b)));
            }
            for (k, m) in ms.modes.iter().enumerate() {
                line(&mut out, &format!("      {:<16} {} {}", m.name, float(m.amplitude), code(Field::Mode(p, k))));
            }
        }
        if let Intensities::List(list) = &ph.intensities {
            let name = format!("{stem}{}.hkl", p + 1);
            let mut t = format!("Reflections of phase {} written by drs\n  h   k   l   mult   intensity\n", p + 1);
            for (h, mult, i) in list {
                t.push_str(&format!("{:>4}{:>4}{:>4}{:>6} {}\n", h[0], h[1], h[2], mult, float(*i)));
            }
            e.files.push((name, t));
        }
        line(&mut out, "!-------> Profile Parameters for Pattern #  1");
        line(&mut out, "!  Scale        Shape1      Bov      Str1      Str2      Str3   Strain-Model");
        match hexagonal {
            Some((v, k)) => {
                line(&mut out, &format!(" {} {} {} {} {} {} 8", float(ph.scale), float(ph.eta0), float(ph.bov), float(v[0]), float(v[1]), float(v[2])));
                line(&mut out, &format!(" {} {} {} {}", code(Field::Scale(p)), code(Field::Eta0(p)), code(Field::Bov(p)), k.iter().map(|&k| code(Field::Quartic(p, k))).collect::<Vec<_>>().join(" ")));
            }
            None => {
                // A TOF phase's Shape1 is its extinction (p. 110).
                let (shape1, shape1_code) = match &ph.tof {
                    Some(t) => (t.extinction, code(Field::TofProfile(p, 17))),
                    None => (ph.eta0, code(Field::Eta0(p))),
                };
                line(&mut out, &format!(" {} {} {} 0.0 0.0 0.0 {}", float(ph.scale), float(shape1), float(ph.bov), strain_model.unwrap_or(0)));
                line(&mut out, &format!(" {} {shape1_code} {} 0.0 0.0 0.0", code(Field::Scale(p)), code(Field::Bov(p))));
            }
        }
        if let Some(t) = &ph.tof {
            // Lines 27–33 of a TOF phase (pp. 110–113): the widths, the cell, the rise and decay, Hodges' terms, the
            // absorption.
            let tp = |k: usize| code(Field::TofProfile(p, k));
            line(&mut out, "!  Sig-2     Sig-1     Sig-0    Sigma-Q  Iso-GStrain Iso-GSize Ani-LSize  Size-Model");
            line(&mut out, &format!(" {} {} {} {} 0.0 0.0 0.0 0", float(t.sigma[2]), float(t.sigma[1]), float(t.sigma[0]), float(t.sigma_q)));
            line(&mut out, &format!(" {} {} {} {} 0.0 0.0 0.0", tp(2), tp(1), tp(0), tp(3)));
            line(&mut out, "!  Gam-2     Gam-1     Gam-0  Iso-LorStrain Iso-LorSize");
            line(&mut out, &format!(" {} {} {} 0.0 0.0", float(t.gamma[2]), float(t.gamma[1]), float(t.gamma[0])));
            line(&mut out, &format!(" {} {} {} 0.0 0.0", tp(6), tp(5), tp(4)));
            line(&mut out, "!     a          b         c        alpha      beta       gamma      #Cell Info");
            let cl = &ph.cell;
            line(&mut out, &format!(" {} {} {} {} {} {}", float(cl.a), float(cl.b), float(cl.c), float(cl.alpha), float(cl.beta), float(cl.gamma)));
            line(&mut out, &(0..6).map(|k| code(Field::Cell(p, k))).collect::<Vec<_>>().join(" "));
            line(&mut out, "!  Pref1    Pref2    alph0    beta0    alph1    beta1    alphQ    betaQ");
            line(&mut out, &format!(" {} {} {} {} {} {} {} {}", float(ph.pref[0]), float(ph.pref[1]), float(t.alpha[0]), float(t.beta[0]), float(t.alpha[1]), float(t.beta[1]), float(t.alpha[2]), float(t.beta[2])));
            line(&mut out, &format!(" {} {} {} {} {} {} {} {}", code(Field::Pref(p, 0)), code(Field::Pref(p, 1)), tp(7), tp(10), tp(8), tp(11), tp(9), tp(12)));
            if let Some(th) = t.thermal {
                line(&mut out, "!  alph0t   beta0t   alph1t   beta1t");
                line(&mut out, &format!(" {}", nums(&th)));
                line(&mut out, &format!(" {}", (13..17).map(tp).collect::<Vec<_>>().join(" ")));
            }
            axes_lines(&mut out, p);
            line(&mut out, "!  Abs1    Code    Abs2   Code");
            line(&mut out, &format!(" {} {} 0.0 0.0", float(t.absorption), tp(18)));
        } else {
            line(&mut out, "!       U         V          W           X          Y        GauSiz   LorSiz Size-Model");
            // LorSiz and its Size-Model (L267): ±1 with the vector on LINE 36, −14 the antiphase rule.
            let (lorsiz, size_model, size_code) = ph.shape_size.map_or((0.0, 0, "0.0".to_string()), |sz| (sz.value, sz.fullprof_model(), code(Field::ShapeSize(p))));
            line(&mut out, &format!(" {} {} {} {} {} {} {} {size_model}", float(ph.u), float(ph.v), float(ph.w), float(ph.x), float(ph.y), float(ph.gauss_size), float(lorsiz)));
            line(&mut out, &format!(" {} {} {} {} {} {} {size_code}", code(Field::U(p)), code(Field::V(p)), code(Field::W(p)), code(Field::X(p)), code(Field::Y(p)), code(Field::GaussSize(p))));
            line(&mut out, "!     a          b         c        alpha      beta       gamma      #Cell Info");
            let cl = &ph.cell;
            line(&mut out, &format!(" {} {} {} {} {} {}", float(cl.a), float(cl.b), float(cl.c), float(cl.alpha), float(cl.beta), float(cl.gamma)));
            line(&mut out, &(0..6).map(|k| code(Field::Cell(p, k))).collect::<Vec<_>>().join(" "));
            line(&mut out, "!  Pref1    Pref2      Asy1     Asy2     Asy3     Asy4      S_L      D_L");
            line(&mut out, &format!(" {} {} {} {} {} {} {} {}", float(ph.pref[0]), float(ph.pref[1]), float(ph.asym[0]), float(ph.asym[1]), float(ph.asym[2]), float(ph.asym[3]), float(ph.s_l), float(ph.d_l)));
            line(&mut out, &format!("{} {} {} {} {} {} {} {}", code(Field::Pref(p, 0)), code(Field::Pref(p, 1)), code(Field::Asym(p, 0)), code(Field::Asym(p, 1)), code(Field::Asym(p, 2)), code(Field::Asym(p, 3)), code(Field::SL(p)), code(Field::DL(p))));
            axes_lines(&mut out, p);
        }
        if own_uvw2 {
            line(&mut out, "!Additional U,V,W parameters for Lambda2");
            line(&mut out, &format!(" {}   <--  U2,V2,W2 for lambda(2)", nums(&ph.uvw(2))));
            line(&mut out, "   0.0 0.0 0.0");
        }
        if let Some(crate::microstructure::ShapeSize { shape: crate::microstructure::SizeShape::Platelet(t) | crate::microstructure::SizeShape::Needle(t), .. }) = ph.shape_size {
            line(&mut out, &format!(" {} {} {}             Platelet-Needle vector (Size)", float(t[0]), float(t[1]), float(t[2])));
        }
        if let (Some(q), Some(model)) = (&ph.quartic_strain, strain_model) {
            // LINE 37 in FullProf's order for the Strain-Model, LINE 38 ξ (manual pp. 119–121).
            let order = crate::microstructure::fullprof_strain_model(model).unwrap();
            let k_of = |t: &[u8; 3]| q.terms.iter().position(|x| &x.0 == t).unwrap();
            for chunk in order.chunks(if order.len() > 6 { 5 } else { 6 }) {
                line(&mut out, &format!("! {}", chunk.iter().map(|t| crate::microstructure::term_name(*t)).collect::<Vec<_>>().join("  ")));
                line(&mut out, &format!(" {}", chunk.iter().map(|t| float(q.terms[k_of(t)].1)).collect::<Vec<_>>().join(" ")));
                line(&mut out, &format!(" {}", chunk.iter().map(|t| code(Field::Quartic(p, k_of(t)))).collect::<Vec<_>>().join(" ")));
            }
            line(&mut out, "!  Lorentzian strain coeff.+ code");
            line(&mut out, &format!(" {} {}", float(q.lorentzian), code(Field::QuarticShare(p))));
        }
        if hexagonal.is_some() {
            // As FullProf writes it after strain model 8 (scor_rxn.pcr): ξ, 0 here.
            line(&mut out, "!  Lorentzian strain coeff.+ code");
            line(&mut out, " 0.0 0.0");
        }
        if !restraints.is_empty() {
            let Some(st) = atoms else { return Err(format!("phase {}: restraints on a phase without atoms", p + 1)) };
            let hex = ph.group.crystal_system == "hexagonal" || ph.group.crystal_system == "trigonal" && (ph.cell.gamma - 120.0).abs() < 1e-6;
            let lab = |a: usize| st.atoms[a].label.clone();
            line(&mut out, "! Soft distance constraints:");
            for r in &restraints {
                if let Kind::Distance { atoms, op } = &r.kind {
                    let n = itnum(op, hex).ok_or(format!("phase {}: a restraint's operator is not one of FullProf's ITnum", p + 1))?;
                    line(&mut out, &format!("{:<6} {:<6} {n:>3} {} {} {} {} {}", lab(atoms[0]), lab(atoms[1]), float(op.t[0]), float(op.t[1]), float(op.t[2]), float(r.target), float(r.sigma)));
                }
            }
            line(&mut out, "! Soft angle constraints:");
            for r in &restraints {
                if let Kind::Angle { atoms, ops } = &r.kind {
                    let n1 = itnum(&ops[0], hex).ok_or(format!("phase {}: a restraint's operator is not one of FullProf's ITnum", p + 1))?;
                    let n2 = itnum(&ops[1], hex).ok_or(format!("phase {}: a restraint's operator is not one of FullProf's ITnum", p + 1))?;
                    line(
                        &mut out,
                        &format!(
                            "{:<6} {:<6} {:<6} {n1:>3} {n2:>3} {} {} {} {} {} {} {} {}",
                            lab(atoms[0]),
                            lab(atoms[1]),
                            lab(atoms[2]),
                            float(ops[0].t[0]),
                            float(ops[0].t[1]),
                            float(ops[0].t[2]),
                            float(ops[1].t[0]),
                            float(ops[1].t[1]),
                            float(ops[1].t[2]),
                            float(r.target),
                            float(r.sigma)
                        ),
                    );
                }
            }
        }
        // LINE 43 (Furth lines): NAMEPAR VALUEPAR CODEPAR, T in Å², L in rad², S in Å·rad (L367). FullProf reads the
        // first TLS group's; its centre is FullProf's own (P6 = 5.01, the centre of mass; 5.0, the origin).
        if let Some((g, gr)) = tls {
            let n = 6 + 3 * gr.members.len();
            line(&mut out, "! TLS parameters:");
            for (i, name) in crate::rigid_body::TLS_NAMES.iter().enumerate() {
                line(&mut out, &format!("{name:<6} {} {}", float(gr.value(n + i)), code(Field::Rigid(p, g, n + i))));
            }
        }
    }
    line(&mut out, "!  2Th1/TOF1    2Th2/TOF2  Pattern to plot");
    line(&mut out, &format!(" {} {}       1", float(thmin), float(thmax)));
    if nli > 0 {
        let number = |name: &str| set.free.iter().position(|f| f.name == name).and_then(|k| number.get(&k).copied());
        for l in &c.refinement.linear {
            line(&mut out, &format!("{} {} {} {}", l.name.replace(' ', "_"), l.terms.len(), float(l.target), float(l.sigma)));
            let mut t = Vec::new();
            for (n, k) in &l.terms {
                t.push(format!("{} {}", float(*k), number(n).ok_or(format!("linear restraint {}: '{n}' is not a free parameter", l.name))?));
            }
            line(&mut out, &t.join(" "));
        }
    }
    e.pcr = out;
    Ok(e)
}
