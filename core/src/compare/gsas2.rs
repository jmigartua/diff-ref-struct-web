//! The comparison with GSAS-II (the app's request 50; `compare`'s notes).
//!
//! The run's folder holds the data as `compare.xye` (x, y, σ: GSAS-II's "Topas xye" reader), the instrument file
//! `compare.instprm`, one CIF per phase (`phase1.cif`, …), `compare_input.json` (every choice the core made: values the
//! CIFs do not hold, the refinement flags, the cycles) and the script `compare.py` (`compare/gsas2.py`, fixed: it
//! decides nothing of its own). GSAS-II's Python runs it with the GSAS-II folder on PYTHONPATH; it writes
//! `compare.gpx`, GSAS-II's listing `compare.lst` and `compare_output.json`, which `gsas2_read` reads.
//!
//! The model in GSAS-II's terms (the parameter dictionary's `gsas2` fields, checked against GSAS-II 5843's source):
//! - widths: GSAS-II's U, V, W are a Gaussian variance σ² in centideg² (GSASIIstrMath.py `GetReflSigGamCW`), the
//!   core's FWHM² in deg², so U_GSAS = U·10⁴/(8 ln 2); GSAS-II's Lorentzian γ = X/cos θ + Y tan θ in centideg, the
//!   core's H_L = X tan θ + Y/cos θ in deg, so X_GSAS = 100·Y and Y_GSAS = 100·X; GSAS-II's profile is always
//!   Thompson–Cox–Hastings's pseudo-Voigt, shared by the phases, with Finger–Cox–Jephcoat's SH/L = S/L + H/L;
//! - Zero is added to 2θ in both (GSASIIstrMath.py `GetReflPos`);
//! - Bragg–Brentano's sample shift and transparency, Δ2θ = −(0.36/πR)·Shift·cos θ − (9/πR)·T·sin 2θ (R the
//!   goniometer radius in mm, Shift in µm), are FullProf's SyCos and SySin;
//! - the polarisation (1 − P)cos²2θ + P (GSASIIpwd.py `Polarization`, azimuth 0) is the core's 1 + Cthm·cos²2θ for
//!   P = 1/(1 + Cthm), and its 1 − K + K·m·cos²2θ (Ilo = 3) for P = (1 − K)/(1 − K + K·m);
//! - atoms: U_iso = (B + B_overall)/8π², the occupancy chemical (FullProf's Occ × M/m), coordinates as they are;
//! - the background: GSAS-II's own Chebyshev polynomial, refined; the core's background is not carried.

use super::*;
use std::path::Path;

/// The code of a GSAS-II run: its script, data, instrument file and outputs start with it.
pub const GSAS2_CODE: &str = "compare";

const GSAS2_SCRIPT: &str = include_str!("gsas2.py");

/// σ² in centideg² per FWHM² in deg² (GSAS-II's U, V, W from the core's).
fn gsas2_uvw() -> f64 {
    1e4 / (8.0 * std::f64::consts::LN_2)
}

/// GSAS-II's goniometer radius (mm), its default, in which the shift and transparency are written.
const GONIO_RADIUS: f64 = 200.0;

/// An element's symbol with its charge as GSAS-II writes it: "TI+4" → "Ti+4".
fn gsas2_species(s: &str) -> String {
    let (el, rest): (String, String) = {
        let k = s.find(|c: char| !c.is_ascii_alphabetic()).unwrap_or(s.len());
        (s[..k].to_string(), s[k..].to_string())
    };
    let mut c = el.chars();
    let el = match c.next() {
        Some(f) => f.to_ascii_uppercase().to_string() + &c.as_str().to_ascii_lowercase(),
        None => String::new(),
    };
    el + &rest
}

/// A phase as a CIF that GSAS-II reads: the cell, the Hermann–Mauguin symbol and the operators, the atoms with
/// U_iso = (B + B_overall)/8π² and the chemical occupancy.
fn gsas2_cif(m: &Model, p: usize) -> Result<String, String> {
    use std::fmt::Write;
    let ph = &m.phases[p];
    let st = match &ph.intensities {
        crate::powder::Intensities::Structure(st, _) => st,
        crate::powder::Intensities::List(_) => return Err(format!("phase {}: a phase without atoms (Le Bail, Pawley) is not carried to GSAS-II yet", p + 1)),
    };
    let mut w = String::new();
    let _ = writeln!(w, "data_phase{}\n_pd_phase_name '{}'", p + 1, ph.name.replace('\'', ""));
    let c = &st.cell;
    for (n, v) in ["length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma"].iter().zip([c.a, c.b, c.c, c.alpha, c.beta, c.gamma]) {
        let _ = writeln!(w, "_cell_{n} {v}");
    }
    if let Some(s) = &st.group.setting {
        let _ = writeln!(w, "_space_group_name_H-M_alt '{}'", s.hm.replace('_', "").split(':').next().unwrap_or("").trim());
    }
    let _ = writeln!(w, "loop_\n_space_group_symop_operation_xyz");
    for op in &st.group.ops {
        let _ = writeln!(w, "'{}'", op.to_xyz());
    }
    let _ = writeln!(w, "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_U_iso_or_equiv\n_atom_site_occupancy");
    let b8 = 8.0 * std::f64::consts::PI * std::f64::consts::PI;
    for a in &st.atoms {
        let chem = crate::structure_factor::chemical_from_occupancy(&st.group, a.x, a.occ)?;
        let _ = writeln!(w, "{} {} {} {} {} {} {}", a.label, gsas2_species(&a.species), a.x[0], a.x[1], a.x[2], (a.biso + ph.bov) / b8, chem);
    }
    Ok(w)
}

/// GSAS-II's polarisation P for the core's Lorentz–polarisation (see the notes above); None where it has none.
fn gsas2_polarisation(m: &Model) -> Option<f64> {
    match m.ilo {
        0 | 1 => Some(1.0 / (1.0 + m.cthm)),
        3 => {
            let mm = if m.cthm == 0.0 { 1.0 } else { m.cthm };
            Some((1.0 - m.rpolarz) / (1.0 - m.rpolarz + m.rpolarz * mm))
        }
        _ => None,
    }
}

/// Prepare a GSAS-II run of a control file in `dir` (see the notes above). The app runs `Run::program` (GSAS-II's
/// Python, which it gives) with `Run::args` in `Run::dir`, `PYTHONPATH` set to the GSAS-II folder that holds
/// `GSASII`. Refused: a time-of-flight pattern, a phase without atoms, data without points.
pub fn gsas2_prepare(opened: &Opened, dir: &Path, cycles: Option<usize>) -> Result<Run, String> {
    use crate::powder::{BackgroundModel, Intensities};
    let cycles = cycles_of(opened, cycles);
    let m = &opened.control.model;
    if m.tof.is_some() {
        return Err("a time-of-flight pattern is not carried to GSAS-II yet".into());
    }
    let data = opened.data.as_ref().ok_or("a comparison needs the control file's data")?;
    std::fs::create_dir_all(dir).map_err(|e| format!("{}: {e}", dir.display()))?;
    let set = opened.control.refinement.build(m)?;
    let mut manifest: Vec<ManifestItem> = Vec::new();
    let mut note = |name: &str, carried: Carried, why: &str| manifest.push(ManifestItem { name: name.into(), carried, why: why.into() });
    // The data: the points with a weight.
    let keep: Vec<usize> = (0..data.x.len()).filter(|&i| data.sigma[i].is_finite() && data.sigma[i] > 0.0).collect();
    if keep.len() < data.x.len() {
        note("", Carried::Approximated, &format!("{} points without a weight are left out of GSAS-II's data", data.x.len() - keep.len()));
    }
    let pat = crate::data::Pattern { x: keep.iter().map(|&i| data.x[i]).collect(), y: keep.iter().map(|&i| data.y[i]).collect(), sigma: keep.iter().map(|&i| data.sigma[i]).collect() };
    let xye = crate::formats::write(&crate::formats::Scan::new("data", crate::formats::Axis::TwoTheta, pat), "xye")?;
    std::fs::write(dir.join(format!("{GSAS2_CODE}.xye")), xye).map_err(|e| e.to_string())?;
    // The instrument file.
    let ph0 = m.phases.first().ok_or("no phase")?;
    let k = gsas2_uvw();
    let two = m.ratio > 0.0 && m.lambda[1] > 0.0;
    let mut ins = String::from("#GSAS-II instrument parameter file; do not add/delete items!\n");
    ins += &format!("Type:{}\nBank:1.0\n", if m.neutrons { "PNC" } else { "PXC" });
    if two {
        ins += &format!("Lam1:{}\nLam2:{}\nI(L2)/I(L1):{}\n", m.lambda[0], m.lambda[1], m.ratio);
    } else {
        ins += &format!("Lam:{}\n", m.lambda[0]);
    }
    let pol = if m.neutrons { 0.5 } else { gsas2_polarisation(m).unwrap_or(0.5) };
    if !m.neutrons && gsas2_polarisation(m).is_none() {
        note("", Carried::Approximated, &format!("the Lorentz–polarisation of Ilo = {} is written as GSAS-II's unpolarised beam (P = 0.5)", m.ilo));
    }
    ins += &format!("Zero:{}\nPolariz.:{pol}\nAzimuth:0.0\n", m.zero);
    // X and Y are Lorentzian widths only in the Thompson–Cox–Hastings shape.
    let (lx, ly) = if ph0.shape == Shape::ThompsonCoxHastings { (100.0 * ph0.y, 100.0 * ph0.x) } else { (0.0, 0.0) };
    ins += &format!("U:{}\nV:{}\nW:{}\nX:{lx}\nY:{ly}\nZ:0.0\nSH/L:{}\n", ph0.u * k, ph0.v * k, ph0.w * k, ph0.s_l + ph0.d_l);
    std::fs::write(dir.join(format!("{GSAS2_CODE}.instprm")), ins).map_err(|e| e.to_string())?;
    // What GSAS-II cannot hold as the core does.
    if ph0.shape != Shape::ThompsonCoxHastings {
        note("", Carried::Approximated, &format!("the peak shape ({:?}) is written as GSAS-II's Thompson–Cox–Hastings pseudo-Voigt from U, V, W, X, Y; η is not carried", ph0.shape));
    }
    if m.phases.len() > 1 && m.phases.iter().any(|q| (q.u, q.v, q.w, q.x, q.y) != (ph0.u, ph0.v, ph0.w, ph0.x, ph0.y)) {
        note("", Carried::Approximated, "GSAS-II's U, V, W, X, Y are the instrument's, shared by the phases: phase 1's are written");
    }
    if (ph0.s_l != 0.0 || ph0.d_l != 0.0) && ph0.s_l != ph0.d_l {
        note("", Carried::Approximated, "the axial divergence's S/L and H/L are written as GSAS-II's one SH/L = S/L + H/L");
    }
    for (p, q) in m.phases.iter().enumerate() {
        if q.asym.iter().any(|&a| a != 0.0) {
            note(&format!("phase{}", p + 1), Carried::Refused, "Bérar–Baldinozzi's asymmetry: GSAS-II has none");
        }
        if q.pref.iter().any(|&a| a != 0.0) && !(m.march && q.pref[0] == 1.0) {
            note(&format!("phase{}", p + 1), Carried::Refused, "preferred orientation is not carried to GSAS-II yet");
        }
        if let Intensities::Structure(st, _) = &q.intensities {
            if st.atoms.iter().any(|a| a.beta.is_some()) {
                note(&format!("phase{}", p + 1), Carried::Approximated, "anisotropic displacements are written as their atoms' B only");
            }
        }
    }
    if !matches!(m.instrumental, crate::powder::InstrumentalWidths::None) {
        note("", Carried::Approximated, "the resolution file's widths are not carried: GSAS-II's U, V, W, X, Y are phase 1's alone");
    }
    if m.mur != 0.0 {
        note("", Carried::Refused, "the absorption correction μR is not carried to GSAS-II yet");
    }
    if !matches!(m.background, BackgroundModel::Polynomial { .. } | BackgroundModel::Points(_) | BackgroundModel::Spline(_)) || !m.background_terms.profiles.is_empty() {
        note("", Carried::Approximated, "the background is GSAS-II's Chebyshev polynomial, refined");
    }
    // The flags.
    let geometry = if m.neutrons { "Debye-Scherrer" } else { "Bragg-Brentano" };
    let mut instrument: Vec<&str> = Vec::new();
    let mut sample: Vec<&str> = Vec::new();
    let mut refine_cell = vec![false; m.phases.len()];
    let mut refine_scale = vec![false; m.phases.len()];
    let mut atom_codes: Vec<Vec<String>> = m.phases.iter().map(|q| match &q.intensities { Intensities::Structure(st, _) => vec![String::new(); st.atoms.len()], _ => Vec::new() }).collect();
    let mut background_terms = 0usize;
    for fr in &set.free {
        for &(f, _) in &fr.fields {
            let a = crate::refinable::address(f, m);
            let (carried, why): (Carried, String) = match f {
                Field::Zero => (Carried::Carried, "Zero".into()),
                Field::U(0) | Field::V(0) | Field::W(0) => {
                    instrument.push(match f { Field::U(_) => "U", Field::V(_) => "V", _ => "W" });
                    (Carried::Carried, "U, V, W, as σ² in centideg²".into())
                }
                Field::X(0) | Field::Y(0) if ph0.shape != Shape::ThompsonCoxHastings => (Carried::Refused, "not a Lorentzian width outside the Thompson–Cox–Hastings shape (the pseudo-Voigt's η slope): GSAS-II has no counterpart; left fixed".into()),
                Field::X(0) => (Carried::Carried, "GSAS-II's Y (its tan θ term), in centideg".into()),
                Field::Y(0) => (Carried::Carried, "GSAS-II's X (its 1/cos θ term), in centideg".into()),
                Field::SL(0) | Field::DL(0) => (Carried::Approximated, "GSAS-II's SH/L = S/L + H/L".into()),
                Field::Sycos if !m.neutrons => (Carried::Carried, "GSAS-II's Shift (µm)".into()),
                Field::Sysin if !m.neutrons => (Carried::Carried, "GSAS-II's Transparency".into()),
                Field::Cell(p, _) => {
                    refine_cell[p] = true;
                    (Carried::Carried, "the cell, refined by GSAS-II under the group's symmetry".into())
                }
                Field::Scale(p) => {
                    refine_scale[p] = true;
                    (Carried::Carried, "GSAS-II's phase fraction (values not compared: another normalisation)".into())
                }
                Field::AtomX(p, i, _) => {
                    atom_codes[p][i].push('X');
                    (Carried::Carried, "GSAS-II's X: all of the atom's free coordinates".into())
                }
                Field::AtomB(p, i) => {
                    atom_codes[p][i].push('U');
                    (Carried::Carried, "GSAS-II's U_iso".into())
                }
                Field::AtomOcc(p, i) => {
                    atom_codes[p][i].push('F');
                    (Carried::Carried, "GSAS-II's chemical occupancy".into())
                }
                Field::Background(_) | Field::BackgroundTerm(_) => {
                    background_terms += 1;
                    (Carried::Approximated, "the background is GSAS-II's Chebyshev polynomial, refined; not compared".into())
                }
                _ => (Carried::Refused, "GSAS-II refines no counterpart here: left fixed".into()),
            };
            match f {
                Field::Zero => instrument.push("Zero"),
                Field::X(0) if ph0.shape == Shape::ThompsonCoxHastings => instrument.push("Y"),
                Field::Y(0) if ph0.shape == Shape::ThompsonCoxHastings => instrument.push("X"),
                Field::SL(0) | Field::DL(0) => instrument.push("SH/L"),
                Field::Sycos if !m.neutrons => sample.push("Shift"),
                Field::Sysin if !m.neutrons => sample.push("Transparency"),
                _ => {}
            }
            if !matches!(f, Field::Background(_) | Field::BackgroundTerm(_)) && !manifest.iter().any(|x| x.name == a) {
                manifest.push(ManifestItem { name: a, carried, why });
            }
        }
    }
    if background_terms > 0 {
        manifest.push(ManifestItem { name: String::new(), carried: Carried::Approximated, why: format!("the background's {background_terms} free values: GSAS-II refines its own Chebyshev polynomial instead; not compared") });
    }
    instrument.sort_unstable();
    instrument.dedup();
    // The CIFs and the input.
    let mut phases = Vec::new();
    for (p, q) in m.phases.iter().enumerate() {
        let cif = format!("phase{}.cif", p + 1);
        std::fs::write(dir.join(&cif), gsas2_cif(m, p)?).map_err(|e| e.to_string())?;
        let labels: Vec<String> = match &q.intensities {
            Intensities::Structure(st, _) => st.atoms.iter().map(|a| a.label.clone()).collect(),
            _ => Vec::new(),
        };
        let atoms: Vec<serde_json::Value> = labels.iter().zip(&atom_codes[p]).map(|(l, c)| {
            let mut c: Vec<char> = c.chars().collect();
            c.sort_unstable();
            c.dedup();
            serde_json::json!({ "label": l, "refine": c.into_iter().collect::<String>() })
        }).collect();
        phases.push(serde_json::json!({ "cif": cif, "name": format!("phase{}", p + 1), "refine_cell": refine_cell[p], "refine_scale": refine_scale[p], "atoms": atoms }));
    }
    // Sample shift (µm) and transparency from SyCos and SySin.
    let r = GONIO_RADIUS;
    let mut sample_values = serde_json::Map::new();
    if !m.neutrons {
        sample_values.insert("Shift".into(), serde_json::json!(-m.sycos * std::f64::consts::PI * r / 0.36));
        sample_values.insert("Transparency".into(), serde_json::json!(-m.sysin * std::f64::consts::PI * r / 9.0));
    }
    let input = serde_json::json!({
        "data": format!("{GSAS2_CODE}.xye"),
        "instprm": format!("{GSAS2_CODE}.instprm"),
        "geometry": geometry,
        "goniometer_radius": r,
        "sample_values": sample_values,
        "excluded": m.excluded.iter().map(|&(a, b)| vec![a, b]).collect::<Vec<_>>(),
        "background": { "type": "chebyschev-1", "terms": background_terms.clamp(1, 12) },
        "instrument_refine": instrument,
        "sample_refine": sample,
        "phases": phases,
        "cycles": cycles,
    });
    std::fs::write(dir.join(format!("{GSAS2_CODE}_input.json")), serde_json::to_string_pretty(&input).map_err(|e| e.to_string())?).map_err(|e| e.to_string())?;
    std::fs::write(dir.join(format!("{GSAS2_CODE}.py")), GSAS2_SCRIPT).map_err(|e| e.to_string())?;
    let run = Run {
        other: "GSAS-II".into(),
        program: "python".into(),
        args: vec![format!("{GSAS2_CODE}.py")],
        dir: dir.to_path_buf(),
        env: vec![("PYTHONPATH".into(), "the GSAS-II folder that holds GSASII".into())],
        cycles,
        code: GSAS2_CODE.into(),
        data: format!("{GSAS2_CODE}.xye"),
        manifest,
    };
    run.save()?;
    Ok(run)
}

/// Read a GSAS-II run prepared by `gsas2_prepare` and run, refine `ours` (a copy) with the same free parameters and
/// cycles, and compare the two.
pub fn gsas2_read(dir: &Path, run: &Run, ours: &Opened) -> Result<Comparison, String> {
    use serde_json::Value;
    let path = dir.join(format!("{}_output.json", run.code));
    let t = std::fs::read_to_string(&path).map_err(|e| format!("{}: {e} (did GSAS-II run?)", path.display()))?;
    let g: Value = serde_json::from_str(&t).map_err(|e| format!("{}: {e}", path.display()))?;
    let o = refine_ours(ours, run.cycles)?;
    let m = &o.model;
    let r = &o.refinement;
    let data = ours.data.as_ref().ok_or("the control file's data are not read")?;
    let num = |v: &Value| v.as_f64();
    let varied = |name: &str| -> Option<(f64, f64)> {
        let a = g["varied"].get(name)?.as_array()?;
        Some((a.first()?.as_f64()?, a.get(1)?.as_f64()?))
    };
    let mut notes: Vec<String> = Vec::new();
    // Agreement: GSAS-II's reduced χ² is GOF²; R, wR, wRmin are Rp, Rwp, Rexp; Rb, wRb the background-free.
    let s = &r.stats;
    let res = &g["residuals"];
    let gof = num(&g["gof"]);
    let mut agreement = vec![
        measure("chi2", "χ²", None, s.chi2, gof.map(|v| v * v), None),
        measure("rp", "Rp (%)", None, s.rp, num(&res["R"]), None),
        measure("rwp", "Rwp (%)", None, s.rwp, num(&res["wR"]), None),
        measure("rexp", "Rexp (%)", None, s.rexp, num(&res["wRmin"]), None),
        measure("rp_background_free", "Rp, background-free (%)", None, s.rp_background_free, num(&res["Rb"]), None),
        measure("rwp_background_free", "Rwp, background-free (%)", None, s.rwp_background_free, num(&res["wRb"]), None),
    ];
    let gphases = g["phases"].as_array().cloned().unwrap_or_default();
    for (k, &(rb, rf)) in s.phases.iter().enumerate() {
        let pid = gphases.get(k).and_then(|p| p["id"].as_i64()).unwrap_or(k as i64);
        agreement.push(measure("r_bragg", &format!("R_Bragg, phase {} (%)", k + 1), Some(k), rb, num(&res[format!("{pid}:0:Rf^2")]), None));
        agreement.push(measure("r_f", &format!("R_F, phase {} (%)", k + 1), Some(k), rf, num(&res[format!("{pid}:0:Rf")]), None));
    }
    notes.push("GSAS-II's R_Bragg here is its R on F² (Rf^2) and its R_F its Rf; both over its own partition of the observed intensities".into());
    let n_core = s.n;
    if let Some(nobs) = g["nobs"].as_u64() {
        if nobs as usize != n_core {
            notes.push(format!("the points used differ: {n_core} here, {nobs} in GSAS-II"));
        }
    }
    // Parameters.
    let k = gsas2_uvw();
    let b8 = 8.0 * std::f64::consts::PI * std::f64::consts::PI;
    let r_mm = GONIO_RADIUS;
    let mut parameters = Vec::new();
    let mut not_compared = Vec::new();
    let mut background = 0usize;
    let cell_keys = ["length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma"];
    for (j, free) in o.set.free.iter().enumerate() {
        let f = free.fields.iter().map(|x| x.0).find(|&f| crate::refinable::address(f, m) == free.name).unwrap_or(free.fields[0].0);
        let address = crate::refinable::address(f, m);
        // (GSAS-II's name, value, esd) in the default model's terms.
        let other: Result<(String, f64, f64), String> = match f {
            Field::Zero => varied(":0:Zero").map(|(v, e)| (":0:Zero".into(), v, e)).ok_or_default(),
            Field::U(0) => varied(":0:U").map(|(v, e)| (":0:U".into(), v / k, e / k)).ok_or_default(),
            Field::V(0) => varied(":0:V").map(|(v, e)| (":0:V".into(), v / k, e / k)).ok_or_default(),
            Field::W(0) => varied(":0:W").map(|(v, e)| (":0:W".into(), v / k, e / k)).ok_or_default(),
            Field::X(0) | Field::Y(0) if m.phases[0].shape != Shape::ThompsonCoxHastings => Err("not a Lorentzian width in this shape: GSAS-II has no counterpart".into()),
            Field::X(0) => varied(":0:Y").map(|(v, e)| (":0:Y".into(), v / 100.0, e / 100.0)).ok_or_default(),
            Field::Y(0) => varied(":0:X").map(|(v, e)| (":0:X".into(), v / 100.0, e / 100.0)).ok_or_default(),
            Field::Sycos if !m.neutrons => varied(":0:Shift").map(|(v, e)| (":0:Shift".into(), -0.36 * v / (std::f64::consts::PI * r_mm), 0.36 * e / (std::f64::consts::PI * r_mm))).ok_or_default(),
            Field::Sysin if !m.neutrons => varied(":0:Transparency").map(|(v, e)| (":0:Transparency".into(), -9.0 * v / (std::f64::consts::PI * r_mm), 9.0 * e / (std::f64::consts::PI * r_mm))).ok_or_default(),
            Field::Cell(p, c) => {
                let ph = gphases.get(p).ok_or("no such phase in GSAS-II's output".to_string())?;
                match (num(&ph["cell"][cell_keys[c]]), num(&ph["cell_esd"][cell_keys[c]])) {
                    (Some(v), Some(e)) => Ok((format!("{}::{}", ph["id"], cell_keys[c]), v, e)),
                    _ => Err("GSAS-II printed no esd for it".into()),
                }
            }
            Field::AtomX(p, i, _) | Field::AtomB(p, i) | Field::AtomOcc(p, i) => {
                let ph = gphases.get(p).ok_or("no such phase in GSAS-II's output".to_string())?;
                let pid = ph["id"].as_i64().unwrap_or(p as i64);
                let atom = &ph["atoms"][i];
                match f {
                    Field::AtomX(..) => {
                        let c = c_of(f);
                        let name = format!("{pid}::dA{}:{i}", ["x", "y", "z"][c]);
                        match (num(&atom["xyz"][c]), varied(&name)) {
                            (Some(v), Some((_, e))) => Ok((name, v, e)),
                            _ => Err("GSAS-II did not refine this coordinate (fixed by the site's symmetry, or tied)".into()),
                        }
                    }
                    Field::AtomB(..) => {
                        let name = format!("{pid}::AUiso:{i}");
                        match (num(&atom["uiso"]), varied(&name)) {
                            (Some(u), Some((_, e))) => Ok((name, b8 * u - m.phases[p].bov, b8 * e)),
                            _ => Err("GSAS-II did not refine U_iso".into()),
                        }
                    }
                    _ => {
                        let name = format!("{pid}::Afrac:{i}");
                        let st = match &m.phases[p].intensities {
                            crate::powder::Intensities::Structure(st, _) => st,
                            _ => return Err("no atoms".into()),
                        };
                        let q = crate::structure_factor::occupancy_from_chemical(&st.group, st.atoms[i].x, 1.0)?;
                        match (num(&atom["frac"]), varied(&name)) {
                            (Some(v), Some((_, e))) => Ok((name, v * q, e * q)),
                            _ => Err("GSAS-II did not refine the occupancy".into()),
                        }
                    }
                }
            }
            Field::Scale(_) => Err("GSAS-II's phase fraction is on its own normalisation: not compared".into()),
            Field::Background(_) | Field::BackgroundTerm(_) => Err("GSAS-II refines its own Chebyshev background: not compared".into()),
            _ => Err(run.manifest.iter().find(|x| x.name == address && x.carried == Carried::Refused).map_or("GSAS-II has no counterpart here".into(), |x| x.why.clone())),
        };
        match other {
            Ok((other_name, other_value, other_esd)) => {
                let core_value = m.get(f);
                let core_esd = r.params.get(j).map_or(f64::NAN, |p| p.esd);
                let delta_sigma = if other_esd > 0.0 { (core_value - other_value) / other_esd } else if core_value == other_value { 0.0 } else { f64::INFINITY };
                parameters.push(ParameterRow { address, other_name, core_value, core_esd, other_value, other_esd, delta_sigma, esd_ratio: if other_esd > 0.0 { core_esd / other_esd } else { f64::NAN }, agrees: delta_sigma.abs() <= 1.0, note: None });
            }
            Err(_) if matches!(f, Field::Background(_) | Field::BackgroundTerm(_)) => background += 1,
            Err(why) => not_compared.push(NotCompared { name: address, refined_by: "this app".into(), why }),
        }
    }
    if background > 0 {
        not_compared.push(NotCompared { name: "the background".into(), refined_by: "both".into(), why: format!("{background} values here; GSAS-II refines its own Chebyshev polynomial: not compared") });
    }
    parameters.sort_by(|a, b| b.delta_sigma.abs().total_cmp(&a.delta_sigma.abs()));
    // The patterns.
    let to_vec = |v: &Value| -> Vec<f64> { v.as_array().map_or(Vec::new(), |a| a.iter().map(|x| x.as_f64().unwrap_or(f64::NAN)).collect()) };
    let (gx, gy) = (to_vec(&g["x"]), to_vec(&g["ycalc"]));
    let ycalc_core = m.calculate(&data.x).ycalc;
    let ycalc_other = if gx.is_empty() { vec![f64::NAN; data.x.len()] } else { at_points(&data.x, &gx, &gy) };
    let difference = ycalc_core.iter().zip(&ycalc_other).map(|(a, b)| a - b).collect();
    Ok(Comparison {
        core: Program { name: "diff-ref-struct".into(), version: format!("core {}", env!("CARGO_PKG_VERSION")), cycles: r.cycles, converged: r.converged },
        other: Program {
            name: "GSAS-II".into(),
            version: format!("GSAS-II {} (Python {})", g["version"].as_str().unwrap_or("?"), g["python"].as_str().unwrap_or("?")),
            cycles: run.cycles,
            converged: g["converged"].as_bool().unwrap_or(false),
        },
        agreement,
        parameters,
        rule: RULE.into(),
        not_compared,
        manifest: run.manifest.clone(),
        notes,
        patterns: Patterns { x: data.x.clone(), yobs: data.y.clone(), ycalc_core, ycalc_other, difference },
    })
}

/// A coordinate's index (0 x, 1 y, 2 z).
fn c_of(f: Field) -> usize {
    match f {
        Field::AtomX(_, _, c) => c,
        _ => 0,
    }
}

trait OkOrDefault {
    fn ok_or_default(self) -> Result<(String, f64, f64), String>;
}

impl OkOrDefault for Option<(String, f64, f64)> {
    /// None: GSAS-II did not vary it.
    fn ok_or_default(self) -> Result<(String, f64, f64), String> {
        self.ok_or_else(|| "GSAS-II did not vary it".into())
    }
}
