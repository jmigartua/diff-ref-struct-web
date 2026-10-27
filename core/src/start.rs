//! A refinement started from a data file and what is known of the material, with no previous control file (L153; the
//! workflows W2 and W3 of `APP_WORKFLOWS.md`, the tutorial's start, L152): `drs new`.
//!
//! What is known of each phase is a CIF (its structure: Rietveld), or a formula, a space group and a cell (no atoms:
//! the intensities are extracted, Le Bail). From these and the data file the core writes a control file:
//! - the data file, its reader (the registry's) and its range, the wavelength the file states or the one given;
//! - one phase per structure, its profile started from the data: W = H², H the median FWHM of the peaks the data show
//!   (`data_ops::peak_search` at the proposed prominence, `data::proposed_prominence`), U = V = 0;
//! - the automatic background (arPLS, `data_ops::auto_background`) as points, kept in their own file
//!   (`{stem} background.drb`, the file dictionary's `background_points`);
//! - each Rietveld phase's scale by linear least squares against the data less that background;
//! - free: the scale of each phase with a structure and every background point; the rest fixed;
//! - every value not given and not taken from the data a default, each one listed (`Start::stated`) and written as a
//!   comment at the head of the file, with where it comes from.

use crate::cell::Cell;
use crate::control::{relative, Control, DataSource, FullProfRecord};
use crate::data::Pattern;
use crate::formats::Axis;
use crate::jacobian::Field;
use crate::powder::{BackgroundModel, Convention, InstrumentalWidths, Intensities, Model, Shape, Weighting};
use crate::refinable::Refinement;
use crate::scattering::{Radiation, Scatterer};
use crate::spacegroup::SpaceGroup;
use std::path::Path;

/// What is known of one phase.
#[derive(Clone, Debug, PartialEq)]
pub enum Known {
    /// A CIF's text (`cif::structure`), the phase named `name`.
    Cif { name: String, text: String },
    /// A formula ("Y2O3"), a space group (a symbol `SpaceGroup::from_symbol` reads) and the cell's free values in
    /// the crystal system's order (`cell_from`): a phase without atoms, its intensities extracted (Le Bail).
    Cell { formula: String, group: String, cell: Vec<f64> },
}

/// What is said of the instrument; None is "not given": the data file's, else a default.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Instrument {
    pub neutrons: bool,
    /// λ1, and λ2 with the ratio I(λ2)/I(λ1) when there are two lines.
    pub wavelength: Option<(f64, Option<(f64, f64)>)>,
    /// An X-ray anode ("Cu", "Mo", …): Kα1 and Kα2 of the core's table at the ratio `KALPHA_RATIO`.
    pub anode: Option<String>,
    /// Cthm (FullProf), cos² 2θ_m of the monochromator.
    pub monochromator_polarisation: Option<f64>,
    /// Wdt: each peak calculated within ±Wdt FWHM.
    pub peak_window: Option<f64>,
    /// The asymmetry correction applies below this 2θ (degrees).
    pub asymmetry_limit: Option<f64>,
    pub profile: Option<Shape>,
}

/// I(Kα2)/I(Kα1) taken with a named anode: 0.5, the ratio of the lines' statistical weights (2p3/2 : 2p1/2 = 2 : 1);
/// the parameter dictionary's `wavelength_ratio` ("Cu Kα2/Kα1 about 0.5", 0.43–0.50 in FullProf's files).
pub const KALPHA_RATIO: f64 = 0.5;

/// Wdt when none is given: 8 FWHM, the median of FullProf's files (the parameter dictionary's `peak_window`).
pub const DEFAULT_PEAK_WINDOW: f64 = 8.0;

/// The asymmetry limit when none is given: 170°, the median of FullProf's files (the dictionary's
/// `asymmetry_limit`).
pub const DEFAULT_ASYMMETRY_LIMIT: f64 = 170.0;

/// One value of the new control file that was neither given nor read from the data: its key in the control file,
/// the value as written, and where it comes from.
#[derive(Clone, Debug, PartialEq)]
pub struct Stated {
    pub key: String,
    pub value: String,
    pub source: String,
}

/// A control file started (`start`): the control, its text (the canonical text with the defaults stated at its head),
/// the background points' file (name relative to the control file, text) and what was taken by default.
#[derive(Clone, Debug)]
pub struct Start {
    pub control: Control,
    pub text: String,
    pub background_file: (String, String),
    pub stated: Vec<Stated>,
    /// From the data: the range, the number of peaks found and their median FWHM.
    pub range: (f64, f64),
    pub peaks: usize,
    pub fwhm: f64,
}

/// The six cell parameters of a group's crystal system from its free values, in this order: cubic a; tetragonal and
/// hexagonal (trigonal in hexagonal axes) a, c; rhombohedral axes a, α; orthorhombic a, b, c; monoclinic a, b, c and
/// the angle of the unique axis; triclinic a, b, c, α, β, γ. Six values are taken as they are.
pub fn cell_from(group: &SpaceGroup, v: &[f64]) -> Result<Cell, String> {
    // Rhombohedral axes: the Hall symbols of the R groups in them carry the 3-fold along the body diagonal, "P 3*".
    let rhombohedral = group.crystal_system == "trigonal" && group.hall.contains('*');
    let c = match (group.crystal_system, v.len()) {
        (_, 6) => [v[0], v[1], v[2], v[3], v[4], v[5]],
        ("cubic", 1) => [v[0], v[0], v[0], 90.0, 90.0, 90.0],
        ("trigonal", 2) if rhombohedral => [v[0], v[0], v[0], v[1], v[1], v[1]],
        ("tetragonal", 2) => [v[0], v[0], v[1], 90.0, 90.0, 90.0],
        ("hexagonal" | "trigonal", 2) => [v[0], v[0], v[1], 90.0, 90.0, 120.0],
        ("orthorhombic", 3) => [v[0], v[1], v[2], 90.0, 90.0, 90.0],
        ("monoclinic", 4) => {
            // The unique axis: that of the 2-fold rotation (or the normal of the mirror), as `refinable::cell_rules`.
            let unique = group
                .ops
                .iter()
                .find_map(|op| {
                    let tr = op.r[0][0] + op.r[1][1] + op.r[2][2];
                    (tr == -1 || tr == 1).then(|| (0..3).find(|&i| op.r[i][i] == if tr == -1 { 1 } else { -1 })).flatten()
                })
                .unwrap_or(1);
            let mut c = [v[0], v[1], v[2], 90.0, 90.0, 90.0];
            c[3 + unique] = v[3];
            c
        }
        (system, n) => {
            let want = match system {
                "cubic" => "1 (a)",
                "tetragonal" | "hexagonal" => "2 (a, c)",
                "trigonal" if rhombohedral => "2 (a, α)",
                "trigonal" => "2 (a, c)",
                "orthorhombic" => "3 (a, b, c)",
                "monoclinic" => "4 (a, b, c and the unique axis's angle)",
                _ => "6 (a, b, c, α, β, γ)",
            };
            return Err(format!("a {system} cell is given by {want} values or all six, not {n}"));
        }
    };
    Cell::new(c[0], c[1], c[2], c[3], c[4], c[5]).map_err(|e| format!("the cell does not make sense: {e}"))
}

/// A formula's elements and counts ("Y2O3" → Y 2, O 3); each element must be in the scattering tables.
pub fn formula_elements(formula: &str) -> Result<Vec<(String, f64)>, String> {
    let mut out = Vec::new();
    let c: Vec<char> = formula.chars().filter(|c| !c.is_whitespace()).collect();
    let mut i = 0;
    while i < c.len() {
        if !c[i].is_ascii_uppercase() {
            return Err(format!("formula '{formula}': '{}' does not start an element (a capital letter); brackets are not read", c[i]));
        }
        let mut el = c[i].to_string();
        i += 1;
        while i < c.len() && c[i].is_ascii_lowercase() {
            el.push(c[i]);
            i += 1;
        }
        let start = i;
        while i < c.len() && (c[i].is_ascii_digit() || c[i] == '.') {
            i += 1;
        }
        let n: f64 = if i > start { c[start..i].iter().collect::<String>().parse().map_err(|_| format!("formula '{formula}': the count of {el}"))? } else { 1.0 };
        Scatterer::lookup(&el, Radiation::Neutron).map_err(|_| format!("formula '{formula}': {el} is not an element"))?;
        out.push((el, n));
    }
    if out.is_empty() {
        return Err("the formula is empty".into());
    }
    Ok(out)
}

/// The reflections a phase without atoms starts its extraction from (Le Bail): every one the group allows up to 5° of
/// 2θ past the data, each with intensity 1, as an imported profile-matching phase without its `.hkl` (`Model::from_pcr`).
pub fn le_bail_list(group: &SpaceGroup, cell: &Cell, lambda_min: f64, two_theta_max: f64) -> Vec<([i32; 3], usize, f64)> {
    let top = (two_theta_max + 5.0).min(179.0);
    let dmin = lambda_min / (2.0 * (top / 2.0).to_radians().sin());
    crate::reflections::generate(group, cell, dmin, 1e9, None).iter().map(|r| (r.hkl, r.multiplicity, 1.0)).collect()
}

/// The median FWHM of the peaks of the data (degrees), and how many there are: the peaks of `data_ops::peak_search`
/// at the proposed prominence and its default significance (3σ).
pub fn median_fwhm(p: &Pattern) -> Option<(f64, usize)> {
    let prom = crate::data::proposed_prominence(p)?;
    let mut w: Vec<f64> = crate::data_ops::peak_search(p, prom, 3.0).iter().map(|k| k.fwhm).filter(|w| *w > 0.0).collect();
    if w.is_empty() {
        return None;
    }
    w.sort_by(f64::total_cmp);
    Some((w[w.len() / 2], w.len()))
}

/// The scale of each Rietveld phase by linear least squares, the others' and the background held: s = Σ w r p / Σ w p²
/// with p the phase's profile at scale 1, r the data less the background and the other phases, w = 1/σ². A phase
/// whose profile is zero keeps its scale.
pub fn least_squares_scales(m: &mut Model, data: &Pattern) {
    for p in 0..m.phases.len() {
        if !matches!(m.phases[p].intensities, Intensities::Structure(..)) {
            continue;
        }
        m.phases[p].scale = 1.0;
        let calc = m.calculate(&data.x);
        let (mut num, mut den) = (0.0, 0.0);
        for i in 0..data.len() {
            let w = if data.sigma[i] > 0.0 && data.sigma[i].is_finite() { 1.0 / (data.sigma[i] * data.sigma[i]) } else { 0.0 };
            let own = calc.per_phase.get(p).map_or(0.0, |c| c[i]);
            let r = data.y[i] - (calc.ycalc[i] - own);
            num += w * r * own;
            den += w * own * own;
        }
        if den > 0.0 && num > 0.0 {
            m.phases[p].scale = num / den;
        }
    }
}

fn stated(out: &mut Vec<Stated>, key: &str, value: impl std::fmt::Display, source: &str) {
    out.push(Stated { key: key.to_string(), value: value.to_string(), source: source.to_string() });
}

/// Start a control file at `out` (it is not written: `Start::write`) from the data file `data` (read by `reader`, or
/// the reader the registry identifies) and what is known of each phase.
pub fn start(data: &Path, reader: Option<&str>, known: &[Known], inst: &Instrument, out: &Path) -> Result<Start, String> {
    if known.is_empty() {
        return Err("give at least one phase: a CIF, or a formula, a space group and a cell".into());
    }
    let ds = crate::formats::read_path(data, reader).map_err(|e| format!("{}: {e}", data.display()))?;
    let scan = ds.scans.into_iter().next().ok_or(format!("{}: no scan", data.display()))?;
    if scan.axis != Axis::TwoTheta {
        return Err(format!("{}: the scan's x is {}, and a refinement needs 2θ (time of flight is version 2, P05)", data.display(), crate::data_ops::axis_name(scan.axis)));
    }
    let pattern = scan.pattern.clone();
    let (lo, hi) = crate::data::extent(&pattern.x).ok_or(format!("{}: no points", data.display()))?;
    let mut defaults = Vec::new();

    // The radiation.
    let anode = |a: &str| crate::scattering_table::ANODES.iter().find(|x| x.0.eq_ignore_ascii_case(a)).copied();
    let (lambda, ratio, said) = match (&inst.wavelength, &inst.anode, &scan.wavelength) {
        (Some((l1, two)), _, _) => (*l1, *two, "given"),
        (None, Some(a), _) => {
            let (_, k1, k2, _) = anode(a).ok_or(format!("no anode '{a}' in the table (Cr, Fe, Co, Ni, Cu, Mo, Ag)"))?;
            stated(&mut defaults, "wavelength_ratio", KALPHA_RATIO, "I(Kα2)/I(Kα1) of the anode's lines: 2 : 1 by the levels' weights (the parameter dictionary's wavelength_ratio)");
            (k1, Some((k2, KALPHA_RATIO)), "the anode's Kα1 and Kα2 (core/src/scattering_table.rs, ANODES)")
        }
        (None, None, Some(w)) => (w.lambda1, w.lambda2.map(|l2| (l2, w.ratio.unwrap_or(KALPHA_RATIO))), "the data file's"),
        (None, None, None) => return Err(format!("{} does not say its wavelength: give it (--wavelength) or the anode (--anode)", data.display())),
    };
    if matches!(&scan.wavelength, Some(w) if w.lambda2.is_some() && w.ratio.is_none()) && inst.wavelength.is_none() && inst.anode.is_none() {
        stated(&mut defaults, "wavelength_ratio", KALPHA_RATIO, "the data file gives two lines without their ratio: I(Kα2)/I(Kα1) = 0.5");
    }
    let lambda2 = ratio.map_or(lambda, |r| r.0);
    let ratio = ratio.map_or(0.0, |r| r.1);
    if inst.neutrons && inst.anode.is_some() {
        return Err("an anode is an X-ray source: neutrons take --wavelength".into());
    }

    // Cthm = cos² 2θ_M (FullProf's manual, eq. 3.26, p. 53): without a monochromator 2θ_M = 0 and Cthm = 1, the
    // polarisation of an unpolarised beam, ½(1 + cos² 2θ). Neutrons are not polarised by it (K = 0).
    let cthm = inst.monochromator_polarisation.unwrap_or_else(|| {
        if inst.neutrons {
            0.0
        } else {
            stated(&mut defaults, "monochromator_polarisation", 1, "no monochromator: cos² 2θ_M with 2θ_M = 0 (FullProf's manual, eq. 3.26, p. 53); graphite with Cu Kα is 0.7998");
            1.0
        }
    });
    let wdt = inst.peak_window.unwrap_or_else(|| {
        stated(&mut defaults, "peak_window", DEFAULT_PEAK_WINDOW, "FWHM: the median of FullProf's files (the parameter dictionary's peak_window)");
        DEFAULT_PEAK_WINDOW
    });
    let asym_limit = inst.asymmetry_limit.unwrap_or_else(|| {
        stated(&mut defaults, "asymmetry_limit", DEFAULT_ASYMMETRY_LIMIT, "deg 2θ: the median of FullProf's files (the parameter dictionary's asymmetry_limit)");
        DEFAULT_ASYMMETRY_LIMIT
    });
    stated(&mut defaults, "geometry", "debye-scherrer", "the Lorentz–polarisation of Bragg–Brentano and Debye–Scherrer geometry, FullProf's Ilo = 0 (its manual, eq. 3.26, p. 53)");
    stated(&mut defaults, "zero_shift, sample_displacement, sample_transparency, absorption_mur, polarisation_fraction", 0, "the parameter dictionary's defaults");
    stated(&mut defaults, "weighting", "variance", "w = 1/σ², σ as the data file gives it or √counts");
    stated(&mut defaults, "[convention]", "all false", "the papers' rules where FullProf departs from them (L050)");

    let mut model = Model {
        convention: Convention::DEFAULT,
        neutrons: inst.neutrons,
        lambda: [lambda, lambda2],
        ratio,
        zero: 0.0,
        sycos: 0.0,
        sysin: 0.0,
        cthm,
        rpolarz: 0.0,
        ilo: 0,
        mur: 0.0,
        wdt,
        asym_limit,
        march: true,
        background: BackgroundModel::Points(Vec::new()),
        instrumental: InstrumentalWidths::None,
        phases: Vec::new(),
        excluded: Vec::new(),
        weighting: Weighting::from_iwg(0),
        restraints: Vec::new(),
        fullprof_step: None,
        microabsorption: [0.0; 3],
        tof: None,
        background_terms: Default::default(),
        background_ceiling: None,
    };

    // The background: arPLS's points.
    let auto = crate::data_ops::auto_background(&pattern, None, crate::data_ops::DEFAULT_MAX_ITER)?;
    stated(&mut defaults, "background", format!("{} points, linear", auto.points.len()), &format!("arPLS (Baek et al. 2015, as GSAS-II) with log₁₀ λ = {} for {} points (GSAS-II's default) and {} iterations at most; a point every N/100 of the scan", auto.log_lambda, pattern.len(), crate::data_ops::DEFAULT_MAX_ITER));
    model.set_background(BackgroundModel::Points(auto.points.clone()))?;

    // The widths from the data.
    let (fwhm, peaks) = median_fwhm(&pattern).ok_or(format!("{}: no peak stands out of the data", data.display()))?;
    let shape = inst.profile.unwrap_or_else(|| {
        stated(&mut defaults, "profile", "thompson-cox-hastings", "the pseudo-Voigt of Thompson, Cox & Hastings (1987), a new phase's (Model::add_phase)");
        Shape::ThompsonCoxHastings
    });
    stated(&mut defaults, "caglioti_w", format!("{:.6}", fwhm * fwhm), &format!("deg²: H², H = {fwhm:.4}° the median FWHM of the {peaks} peaks found (peak_search at the proposed prominence, 3σ); caglioti_u = caglioti_v = 0"));
    stated(&mut defaults, "pv_eta0, pv_eta_slope, lorentzian_size, lorentzian_strain, gaussian_size", 0, "a Gaussian start");
    stated(&mut defaults, "berar_baldinozzi", "[0, 0, 0, 0]", "no asymmetry at the start");
    stated(&mut defaults, "po_g1", 1, "no preferred orientation (March–Dollase)");

    let lambda_min = model.lambda[0].min(model.lambda[1]);
    let mut le_bail = Vec::new();
    for (k, kn) in known.iter().enumerate() {
        let p = match kn {
            Known::Cif { name, text } => model.add_phase_from_cif(name, text).map_err(|e| format!("{name}: {e}"))?,
            Known::Cell { formula, group, cell } => {
                formula_elements(formula)?;
                let (g, others) = SpaceGroup::from_symbol(group)?;
                if !others.is_empty() {
                    stated(&mut defaults, &format!("phase {} space_group", k + 1), g.setting.map_or(g.hall.clone(), |s| s.hm.to_string()), &format!("the first of {} settings '{group}' names; another is chosen by its symbol with ':' (\"P n n n:2\")", others.len() + 1));
                }
                let c = cell_from(&g, cell)?;
                let p = model.add_phase(formula, g.clone(), c)?;
                model.phases[p].intensities = Intensities::List(le_bail_list(&g, &c, lambda_min, hi));
                model.phases[p].extract = true;
                le_bail.push(p);
                p
            }
        };
        let ph = &mut model.phases[p];
        ph.shape = shape;
        ph.w = fwhm * fwhm;
    }
    if !le_bail.is_empty() {
        stated(&mut defaults, "le_bail", "true", "a phase without atoms: its intensities are extracted every cycle (Le Bail), from 1 for every reflection to 5° past the data; its scale fixed");
    }
    least_squares_scales(&mut model, &pattern);
    for (p, ph) in model.phases.iter().enumerate() {
        if matches!(ph.intensities, Intensities::Structure(..)) {
            stated(&mut defaults, &format!("phase {} scale_factor", p + 1), format!("{:.6e}", ph.scale), "linear least squares against the data less the background");
        }
    }

    // Free: each Rietveld phase's scale and every background point.
    let mut parameters: Vec<Vec<(Field, f64)>> = (0..auto.points.len()).map(|k| vec![(Field::Background(k), 1.0)]).collect();
    parameters.extend((0..model.phases.len()).filter(|p| !le_bail.contains(p)).map(|p| vec![(Field::Scale(p), 1.0)]));
    let refinement = Refinement { parameters, cycles: Some(crate::lsq::DEFAULT_CYCLES), eps: Some(crate::lsq::DEFAULT_EPS), ..Default::default() };
    stated(&mut defaults, "cycles, convergence", format!("{}, {}", crate::lsq::DEFAULT_CYCLES, crate::lsq::DEFAULT_EPS), "the core's (lsq::DEFAULT_CYCLES, DEFAULT_EPS)");

    let dir = out.parent().filter(|d| !d.as_os_str().is_empty()).map_or_else(|| std::path::PathBuf::from("."), Path::to_path_buf);
    let stem = out.file_stem().map(|s| s.to_string_lossy().into_owned()).ok_or("the control file needs a name")?;
    let bg_name = crate::formats::file_table::FILES.iter().find(|f| f.id == "background_points").map_or_else(|| format!("{stem} background.drb"), |f| f.pattern.replace("{stem}", &stem));
    let title = known.iter().map(|k| match k {
        Known::Cif { name, .. } => name.clone(),
        Known::Cell { formula, .. } => formula.clone(),
    }).collect::<Vec<_>>().join(" + ");
    let control = Control {
        title: format!("{title}, {}", data.file_name().map_or(String::new(), |n| n.to_string_lossy().into_owned())),
        model,
        data: Some(DataSource { file: relative(data, &dir), reader: Some(ds.format.to_string()), reader_options: Default::default(), operations: Vec::new() }),
        refinement,
        fullprof: FullProfRecord::default(),
        background_file: Some(bg_name.clone()),
    };
    let (_, bg_text) = control.background_file_text(None).ok_or("the background is not points")??;
    let text = with_stated(&control.to_text(), &defaults, &known_sources(data, known, said));
    Ok(Start { control, text, background_file: (bg_name, bg_text), stated: defaults, range: (lo, hi), peaks, fwhm })
}

fn known_sources(data: &Path, known: &[Known], wavelength: &str) -> String {
    let phases: Vec<String> = known.iter().map(|k| match k {
        Known::Cif { name, .. } => format!("{name} (a CIF)"),
        Known::Cell { formula, group, .. } => format!("{formula} ({group}, a cell: Le Bail)"),
    }).collect();
    format!("{} and {}; the wavelength: {wavelength}", data.file_name().map_or(String::new(), |n| n.to_string_lossy().into_owned()), phases.join(", "))
}

/// The canonical text with the stated defaults as comments after its first line.
fn with_stated(text: &str, defaults: &[Stated], from: &str) -> String {
    let mut head = vec![format!("# Started by drs new (L153) from {from}."), "# Free: the scale of each phase with atoms and every background point; the rest is fixed.".to_string(), "# Not given and not read from the data (each a default, with where it comes from):".to_string()];
    head.extend(defaults.iter().map(|d| format!("#   {} = {}: {}", d.key, d.value, d.source)));
    let mut lines: Vec<&str> = text.lines().collect();
    let at = 1.min(lines.len());
    let head_refs: Vec<&str> = head.iter().map(String::as_str).collect();
    lines.splice(at..at, head_refs);
    let mut out = lines.join("\n");
    out.push('\n');
    out
}

impl Start {
    /// Write the control file at `out` and the background points beside it; refused if either exists.
    pub fn write(&self, out: &Path) -> Result<(), String> {
        let dir = out.parent().filter(|d| !d.as_os_str().is_empty()).map_or_else(|| std::path::PathBuf::from("."), Path::to_path_buf);
        let bg = dir.join(&self.background_file.0);
        for p in [out, bg.as_path()] {
            if p.exists() {
                return Err(format!("{} exists: nothing written", p.display()));
            }
        }
        std::fs::write(&bg, &self.background_file.1).map_err(|e| format!("{}: {e}", bg.display()))?;
        std::fs::write(out, &self.text).map_err(|e| format!("{}: {e}", out.display()))
    }
}
