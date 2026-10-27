//! The results file of a refinement (R28; the format chosen by Josu, L092; written L101): a powder CIF (pdCIF) named
//! "{control} {time} results.cif" in `Results/` (the file dictionary's `results`), and beside it one structure CIF per
//! phase ("{control} {phase}.cif", `deposition_cif`).
//!
//! The results file holds every phase, the profile, the reflections, the agreement factors and the model with its
//! standard uncertainties, in data blocks linked by `_pd_block_id` as the powder dictionary links them:
//! - `data_overall`: the agreement factors (`_pd_proc_ls_prof_R_factor`, `_pd_proc_ls_prof_wR_factor`,
//!   `_pd_proc_ls_prof_wR_expected`, `_refine_ls_goodness_of_fit_all`), the numbers of points, parameters and
//!   restraints, the blocks of the phases and of the pattern, each refined parameter with its value and standard
//!   uncertainty, and the control file itself;
//! - one block per structure phase: the deposition CIF's block (`cif_write::deposition`: cell, symmetry, atoms with
//!   their uncertainties, R_Bragg, geometry), with its `_pd_block_id` and `_pd_phase_name`;
//! - `data_pattern`: the radiation and wavelengths, the excluded regions, and the profile point by point (2θ, Yobs, its
//!   weight, the background, Ycal and each phase's curve; Yobs with its σ as su, to four significant digits, since
//!   GSAS-II takes σ from there and not from the weight) and the reflections (phase, wavelength, h k l, multiplicity,
//!   d, 2θ, FWHM, |F|², I calc, I obs partitioned by the model and its σ: `diagnostics::reflections`).
//!
//! Items are named by the IUCr's DDL1 core and powder dictionaries where a name exists; the rest under the local
//! prefix `_drs_` (each phase's curve, the reflections' 2θ and FWHM, the parameters, χ², the control file). GSAS-II's
//! powder CIF importer reads the profile (`G2pwd_CIF.py`: x from `_pd_meas_2theta_scan`, y from
//! `_pd_meas_intensity_total`, w = 1/max(σ, 1)² from its su, else 1/y: at commit 305a159 its line reading
//! `_pd_proc_ls_weight` into w is commented out) and its phase importer the structures; the test
//! `results_file` checks both readings against the core's (`tests/data/make_gsas2_results.py`).

use crate::control::Control;
use crate::data::Pattern;
use crate::lsq::Refinement;
use crate::powder::Intensities;
use crate::refinable::ParameterSet;
use std::fmt::Write;

/// The local prefix of the items the IUCr's dictionaries do not name.
pub const PREFIX: &str = "_drs_";

pub(crate) fn quoted(s: &str) -> String {
    if s.contains('\'') {
        format!("\"{s}\"")
    } else {
        format!("'{s}'")
    }
}

/// A number as the loops write it: seven significant figures, no trailing zeros.
pub(crate) fn num(v: f64) -> String {
    if !v.is_finite() {
        return "?".into();
    }
    if v == 0.0 {
        return "0".into();
    }
    let d = (6 - v.abs().log10().floor() as i32).clamp(0, 15) as usize;
    let s = format!("{v:.d$}");
    let s = if s.contains('.') { s.trim_end_matches('0').trim_end_matches('.').to_string() } else { s };
    if s == "-0" { "0".into() } else { s }
}

/// A value with its standard uncertainty to four significant digits ("212.00(1456)"), so that the weights it gives
/// survive the print; without an uncertainty (0), the value alone.
fn with_sigma(v: f64, su: f64) -> String {
    if !(su > 0.0 && su.is_finite()) {
        return num(v);
    }
    let d = (3 - su.log10().floor() as i32).max(0) as usize;
    let digits = (su * 10f64.powi(d as i32)).round() as i64;
    format!("{v:.d$}({digits})")
}

/// A mode phase's amplitudes, its irreps' and their su's, its reference structure and its cell reduced to the
/// parent's (P2, P11), under the local prefix; nothing for a phase without modes.
pub(crate) fn modes_block(w: &mut String, m: &crate::powder::Model, p: usize, set: Option<&ParameterSet>, refinement: Option<&Refinement>) {
    let Some(ms) = &m.phases[p].modes else { return };
    let field_su = |f: crate::jacobian::Field| refinement.map_or(0.0, |r| {
        let a = crate::refinable::address(f, m);
        r.fields.iter().find(|x| x.0 == a).map_or(0.0, |x| x.2)
    });
    let _ = writeln!(w, "loop_\n{PREFIX}mode_name\n{PREFIX}mode_irrep\n{PREFIX}mode_amplitude\n{PREFIX}mode_amplitude_su");
    for (k, md) in ms.modes.iter().enumerate() {
        let _ = writeln!(w, "{} {} {} {}", quoted(&md.name), quoted(&md.irrep), num(md.amplitude), num(field_su(crate::jacobian::Field::Mode(p, k))));
    }
    let _ = writeln!(w, "loop_\n{PREFIX}irrep_name\n{PREFIX}irrep_amplitude\n{PREFIX}irrep_amplitude_su");
    for ir in ms.irreps() {
        let su = match (set, refinement) {
            (Some(s), Some(r)) => ms.irrep_esd(&ir, |a, b| s.field_covariance(&r.covariance, crate::jacobian::Field::Mode(p, a), crate::jacobian::Field::Mode(p, b))),
            _ => 0.0,
        };
        let _ = writeln!(w, "{} {} {}", quoted(&ir.name), num(ir.amplitude), num(su));
    }
    if let crate::powder::Intensities::Structure(st, _) = &m.phases[p].intensities {
        let _ = writeln!(w, "loop_\n{PREFIX}mode_reference_label\n{PREFIX}mode_reference_fract_x\n{PREFIX}mode_reference_fract_y\n{PREFIX}mode_reference_fract_z");
        for (a, x) in st.atoms.iter().zip(&ms.reference) {
            let _ = writeln!(w, "{} {} {} {}", quoted(&a.label), num(x[0]), num(x[1]), num(x[2]));
        }
    }
}

/// A rigid-body phase's groups (P41, Publish's P47.1), under the local prefix: per group its kind, origin and θ, φ, χ
/// (degrees) with their su's and a satellite's anchor (parent group, atoms N1 and N2); per group atom its internal
/// coordinates with su's (in the order of the kind: r θ φ, x y z, ρ φ z, or distance, angle, torsion) and a
/// Z-matrix's N1 N2 N3. Nothing for a phase without groups.
pub(crate) fn rigid_block(w: &mut String, m: &crate::powder::Model, p: usize, refinement: Option<&Refinement>) {
    let ph = &m.phases[p];
    if ph.rigid.is_empty() {
        return;
    }
    let su = |f: crate::jacobian::Field| refinement.map_or(0.0, |r| {
        let a = crate::refinable::address(f, m);
        r.fields.iter().find(|x| x.0 == a).map_or(0.0, |x| x.2)
    });
    let v = |g: usize, k: usize| crate::cif_write::with_su(ph.rigid[g].value(k), su(crate::jacobian::Field::Rigid(p, g, k)), 5);
    let _ = writeln!(w, "loop_\n{PREFIX}rigid_group_id\n{PREFIX}rigid_group_kind\n{PREFIX}rigid_group_origin_x\n{PREFIX}rigid_group_origin_y\n{PREFIX}rigid_group_origin_z\n{PREFIX}rigid_group_theta\n{PREFIX}rigid_group_phi\n{PREFIX}rigid_group_chi\n{PREFIX}rigid_group_satellite_of\n{PREFIX}rigid_group_satellite_atom\n{PREFIX}rigid_group_satellite_from");
    for (g, gr) in ph.rigid.iter().enumerate() {
        let sat = gr.satellite.map_or(". . .".to_string(), |s| format!("{} {} {}", s.parent + 1, s.atom, s.from));
        let _ = writeln!(w, "{} {} {} {} {} {} {} {} {sat}", g + 1, quoted(gr.coordinates.name()), v(g, 0), v(g, 1), v(g, 2), v(g, 3), v(g, 4), v(g, 5));
    }
    let labels = match &ph.intensities {
        crate::powder::Intensities::Structure(st, _) => crate::geometry::unique_labels(st),
        crate::powder::Intensities::List(_) => Vec::new(),
    };
    let _ = writeln!(w, "loop_\n{PREFIX}rigid_atom_group_id\n{PREFIX}rigid_atom_label\n{PREFIX}rigid_atom_internal_1\n{PREFIX}rigid_atom_internal_2\n{PREFIX}rigid_atom_internal_3\n{PREFIX}rigid_atom_connect_n1\n{PREFIX}rigid_atom_connect_n2\n{PREFIX}rigid_atom_connect_n3");
    for (g, gr) in ph.rigid.iter().enumerate() {
        for (k, mm) in gr.members.iter().enumerate() {
            let c = |n: usize| if gr.coordinates == crate::rigid_body::Coordinates::ZMatrix && n > 0 { n.to_string() } else { ".".to_string() };
            let label = labels.get(mm.index).cloned().unwrap_or_else(|| mm.label.clone());
            let _ = writeln!(w, "{} {} {} {} {} {} {} {}", g + 1, quoted(&label), v(g, 6 + 3 * k), v(g, 7 + 3 * k), v(g, 8 + 3 * k), c(mm.connect[0]), c(mm.connect[1]), c(mm.connect[2]));
        }
    }
}

/// A phase's TLS groups (L367, Publish's P47.2 and T16), under the local prefix: per group its centre (fractions, the
/// group's origin plus `Tls::centre`, with the origin's su's); per value its tensor (T, L, S), component (11 … 23, and
/// S's 21 31 32), value with su in FullProf's printed units (T Å², L deg², S Å·deg) and the free parameter that moves it
/// with its multiplier (`.` for a fixed value), so that the values one parameter moves are the ones the refinement ties
/// (FullProf's codewords, e.g. s22 = −s11); per atom its U11 … U23 (Å²) from T, L and S with su's through the
/// covariance. Nothing for a phase without TLS.
fn tls_block(w: &mut String, m: &crate::powder::Model, p: usize, set: Option<&ParameterSet>, refinement: Option<&Refinement>, prop: Option<&crate::geometry::Propagator>) {
    use crate::jacobian::Field;
    let ph = &m.phases[p];
    if !ph.rigid.iter().any(|g| g.displacements.is_some()) {
        return;
    }
    let su = |f: Field| refinement.map_or(0.0, |r| {
        let a = crate::refinable::address(f, m);
        r.fields.iter().find(|x| x.0 == a).map_or(0.0, |x| x.2)
    });
    let _ = writeln!(w, "loop_\n{PREFIX}rigid_tls_group_id\n{PREFIX}rigid_tls_centre_x\n{PREFIX}rigid_tls_centre_y\n{PREFIX}rigid_tls_centre_z");
    for (g, gr) in ph.rigid.iter().enumerate() {
        let Some(c) = crate::rigid_body::tls_centre(gr) else { continue };
        let cs: Vec<String> = (0..3).map(|k| crate::cif_write::with_su(c[k], su(Field::Rigid(p, g, k)), 5)).collect();
        let _ = writeln!(w, "{} {}", g + 1, cs.join(" "));
    }
    let _ = writeln!(w, "loop_\n{PREFIX}rigid_tls_value_group_id\n{PREFIX}rigid_tls_value_tensor\n{PREFIX}rigid_tls_value_component\n{PREFIX}rigid_tls_value\n{PREFIX}rigid_tls_value_parameter\n{PREFIX}rigid_tls_value_multiplier");
    let deg = 180.0 / std::f64::consts::PI;
    for (g, gr) in ph.rigid.iter().enumerate().filter(|(_, g)| g.displacements.is_some()) {
        let n = 6 + 3 * gr.members.len();
        for i in 0..21 {
            let f = Field::Rigid(p, g, n + i);
            let (tensor, unit) = match i {
                0..=5 => ("T", 1.0),
                6..=11 => ("L", deg * deg),
                _ => ("S", deg),
            };
            let free = set.and_then(|s| s.free.iter().find_map(|fr| fr.fields.iter().find(|x| x.0 == f).map(|x| (fr.name.clone(), x.1))));
            let (name, mult) = free.map_or((".".to_string(), ".".to_string()), |(n, c)| (quoted(&n), num(c)));
            let _ = writeln!(w, "{} {tensor} {} {} {name} {mult}", g + 1, crate::rigid_body::value_component(gr, n + i), crate::cif_write::with_su(gr.value(n + i) * unit, su(f) * unit, 6));
        }
    }
    let labels = match &ph.intensities {
        crate::powder::Intensities::Structure(st, _) => crate::geometry::unique_labels(st),
        crate::powder::Intensities::List(_) => Vec::new(),
    };
    let _ = writeln!(w, "loop_\n{PREFIX}rigid_tls_atom_label\n{PREFIX}rigid_tls_atom_U_11\n{PREFIX}rigid_tls_atom_U_22\n{PREFIX}rigid_tls_atom_U_33\n{PREFIX}rigid_tls_atom_U_12\n{PREFIX}rigid_tls_atom_U_13\n{PREFIX}rigid_tls_atom_U_23");
    for (a, u) in crate::geometry::tls_u(m, p, prop) {
        let us: Vec<String> = u.iter().map(|&(v, s)| crate::cif_write::with_su(v, s, 5)).collect();
        let _ = writeln!(w, "{} {}", quoted(labels.get(a).map_or("?", |l| l.as_str())), us.join(" "));
    }
}

/// A phase's preferred orientation as numbers (P30, Publish's P48.1–2), under the local prefix: the function
/// (`none`, `march-dollase`, `exponential`, `march-dollase axes`), and for one axis its direction with G1 and G2, for
/// several its geometry (`reflection`, or `anvil cell` with the integration steps) and per axis h k l, March's r and
/// the fraction with su's (the last fraction's through the covariance of the others).
pub(crate) fn texture_block(w: &mut String, m: &crate::powder::Model, p: usize, prop: Option<&crate::geometry::Propagator>) {
    use crate::jacobian::Field;
    let ph = &m.phases[p];
    let su = |vals: &dyn Fn(&crate::powder::Model) -> Vec<f64>, moves: &dyn Fn(Field) -> bool| -> Vec<(f64, f64)> {
        match prop {
            Some(pr) => pr.propagate(vals, moves),
            None => vals(m).into_iter().map(|v| (v, 0.0)).collect(),
        }
    };
    let s = |v: (f64, f64)| crate::cif_write::with_su(v.0, v.1, 5);
    match &ph.texture {
        Some(t) => {
            let _ = writeln!(w, "{PREFIX}pref_orient_function          'march-dollase axes'");
            let _ = writeln!(w, "{PREFIX}pref_orient_geometry          {}", if t.steps > 0 { "'anvil cell'" } else { "reflection" });
            if t.steps > 0 {
                let _ = writeln!(w, "{PREFIX}pref_orient_steps             {}", t.steps);
            }
            let v = su(&|mm: &crate::powder::Model| mm.phases[p].texture.as_ref().map_or(Vec::new(), |t| t.axes.iter().flat_map(|a| [a.r, a.fraction]).collect()), &|f| matches!(f, Field::Texture(q, _) if q == p));
            let _ = writeln!(w, "loop_\n{PREFIX}pref_orient_axis_id\n{PREFIX}pref_orient_axis_index_h\n{PREFIX}pref_orient_axis_index_k\n{PREFIX}pref_orient_axis_index_l\n{PREFIX}pref_orient_axis_ratio\n{PREFIX}pref_orient_axis_fraction");
            for (j, a) in t.axes.iter().enumerate() {
                let _ = writeln!(w, "{} {} {} {} {} {}", j + 1, num(a.direction[0]), num(a.direction[1]), num(a.direction[2]), s(v[2 * j]), s(v[2 * j + 1]));
            }
        }
        None => {
            let none = if m.march { ph.pref[0] == 1.0 && ph.pref[1] == 0.0 } else { ph.pref[0] == 0.0 };
            let function = if none { "none" } else if m.march { "march-dollase" } else { "exponential" };
            let _ = writeln!(w, "{PREFIX}pref_orient_function          {function}");
            if !none {
                let v = su(&|mm: &crate::powder::Model| mm.phases[p].pref.to_vec(), &|f| matches!(f, Field::Pref(q, _) if q == p));
                for (k, item) in ["h", "k", "l"].iter().enumerate() {
                    let _ = writeln!(w, "{PREFIX}pref_orient_index_{item}             {}", num(ph.pref_direction[k]));
                }
                let _ = writeln!(w, "{PREFIX}pref_orient_G1                {}", s(v[0]));
                let _ = writeln!(w, "{PREFIX}pref_orient_G2                {}", s(v[1]));
            }
        }
    }
}

/// A structure phase's cell reduced to its parent's (P11; P17 for a phase refined by its coordinates, its parent
/// given), each atom's B_iso (B_eq for an anisotropic atom) with its su (P20), and its atoms' coordination, mean distance and bond-valence sum with su's when every atom has a charge
/// (P18: contacts within 3 Å, Bond_Str's tolerance 0.20), under the local prefix.
pub(crate) fn structure_block(w: &mut String, m: &crate::powder::Model, p: usize, set: Option<&ParameterSet>, refinement: Option<&Refinement>) {
    rigid_block(w, m, p, refinement);
    let prop = match (set, refinement) {
        (Some(s), Some(r)) => Some(crate::geometry::Propagator::new(m, s, r)),
        _ => None,
    };
    tls_block(w, m, p, set, refinement, prop.as_ref());
    texture_block(w, m, p, prop.as_ref());
    if let Some(c) = crate::geometry::parent_cell(m, p, prop.as_ref()) {
        for (k, item) in ["a", "b", "c", "alpha", "beta", "gamma"].iter().enumerate() {
            let _ = writeln!(w, "{PREFIX}parent_cell_{item:<6}            {}", crate::cif_write::with_su(c[k].0, c[k].1, 5));
        }
    }
    // B_iso, or B_eq = 8π²U_eq for an anisotropic atom, with its su through the covariance (P20), by the CIF's labels.
    if let Some(st) = crate::geometry::structure(m, p) {
        let of = |mm: &crate::powder::Model| -> Vec<f64> {
            let s = crate::geometry::structure(mm, p).unwrap();
            s.atoms.iter().map(|a| a.beta.map_or(a.biso, |b| crate::structure_factor::b_from_u(crate::geometry::u_equivalent(&s.cell, &b), 0.0).0)).collect()
        };
        let values: Vec<(f64, f64)> = match &prop {
            Some(pr) => pr.of_phase(p, &of),
            None => of(m).into_iter().map(|v| (v, 0.0)).collect(),
        };
        let labels = crate::geometry::unique_labels(&st);
        let _ = writeln!(w, "loop_\n{PREFIX}atom_site_label\n{PREFIX}atom_site_B_iso_or_equiv");
        for (l, (v, s)) in labels.iter().zip(&values) {
            let _ = writeln!(w, "{} {}", quoted(l), crate::cif_write::with_su(*v, *s, 4));
        }
    }
    if let (Ok(Some(vs)), Some(st)) = (crate::geometry::valence(m, p, prop.as_ref(), 3.0, 0.2), crate::geometry::structure(m, p)) {
        let labels = crate::geometry::unique_labels(&st);
        let _ = writeln!(w, "loop_\n{PREFIX}valence_atom_site_label\n{PREFIX}valence_coordination\n{PREFIX}valence_formal\n{PREFIX}valence_mean_distance\n{PREFIX}valence_sum");
        for v in &vs {
            let _ = writeln!(w, "{} {} {} {} {}", quoted(&labels[v.atom]), v.coordination, num(v.valence), crate::cif_write::with_su(v.mean_distance.0, v.mean_distance.1, 4), crate::cif_write::with_su(v.sum.0, v.sum.1, 3));
        }
    }
}

/// A number to the last bit (Rust's shortest representation that reads back the same), for the cycles' loops.
fn exact(v: f64) -> String {
    if v.is_finite() { format!("{v:?}") } else { "?".into() }
}

/// Every cycle of the refinement (E1), under the local prefix: per cycle χ² at its start and end, Rp, Rwp, Rexp, λ,
/// the singular directions left out and the time; per cycle and parameter the value, the shift applied, the esd, the
/// shift computed and whether it was held at a limit. The numbers to the last bit, so that the states are rebuilt
/// exactly (`lsq::Trajectory::from_end`). Nothing without cycles.
fn cycles_block(w: &mut String, refinement: Option<&Refinement>) {
    let Some(r) = refinement.filter(|r| !r.history.is_empty()) else { return };
    // A plan run's cycles carry the step that ran them (P37): its number and name.
    let steps = r.cycle_steps.len() == r.history.len();
    let _ = writeln!(w, "loop_\n{PREFIX}cycle_number\n{PREFIX}cycle_chi2_start\n{PREFIX}cycle_chi2\n{PREFIX}cycle_R_factor\n{PREFIX}cycle_wR_factor\n{PREFIX}cycle_wR_expected\n{PREFIX}cycle_marquardt\n{PREFIX}cycle_singular\n{PREFIX}cycle_seconds{}", if steps { format!("\n{PREFIX}cycle_step\n{PREFIX}cycle_step_name") } else { String::new() });
    for (i, c) in r.history.iter().enumerate() {
        let s = &c.stats;
        let step = if steps { format!(" {} {}", r.cycle_steps[i].0, quoted(&r.cycle_steps[i].1)) } else { String::new() };
        let _ = writeln!(w, "{} {} {} {} {} {} {} {} {}{step}", c.cycle, exact(c.chi2_start), exact(s.chi2), exact(s.rp), exact(s.rwp), exact(s.rexp), exact(c.marquardt), c.singular, exact(c.seconds));
    }
    let _ = writeln!(w, "loop_\n{PREFIX}cycle_parameter_cycle\n{PREFIX}cycle_parameter_name\n{PREFIX}cycle_parameter_value\n{PREFIX}cycle_parameter_shift\n{PREFIX}cycle_parameter_esd\n{PREFIX}cycle_parameter_computed\n{PREFIX}cycle_parameter_at_limit");
    for c in &r.history {
        for q in &c.params {
            let _ = writeln!(w, "{} {} {} {} {} {} {}", c.cycle, quoted(&q.name), exact(q.value), exact(q.shift), exact(q.esd), exact(q.computed), if q.at_limit { "yes" } else { "no" });
        }
    }
    // The steering's notes as data (E15c): per cycle each note's kind, address, value, and a refused step's λ and χ².
    if r.history.iter().any(|c| !c.steering.is_empty()) {
        let _ = writeln!(w, "loop_\n{PREFIX}steering_cycle\n{PREFIX}steering_kind\n{PREFIX}steering_address\n{PREFIX}steering_value\n{PREFIX}steering_lambda\n{PREFIX}steering_chi2");
        let opt = |x: Option<f64>| x.map_or(".".to_string(), exact);
        for c in &r.history {
            for n in &c.steering {
                let _ = writeln!(w, "{} {} {} {} {} {}", c.cycle, n.kind.word(), n.address.as_deref().map_or(".".to_string(), quoted), opt(n.value), opt(n.lambda), opt(n.chi2));
            }
        }
    }
}

/// The free parameters at the end, each with its value and esd (`Refinement::params`; a field's su above differs from
/// the esd of the parameter named after it when another parameter also moves the field, as Y₂O₃'s second asymmetry
/// parameter moves P1), and their correlations (E11b): the upper triangle of `Refinement::correlations`, by the
/// parameters' names, so that a run read back has each pair's ellipse.
fn correlations_block(w: &mut String, set: Option<&ParameterSet>, refinement: Option<&Refinement>) {
    let (Some(set), Some(r)) = (set, refinement) else { return };
    if r.params.len() != set.free.len() || r.params.is_empty() {
        return;
    }
    let _ = writeln!(w, "loop_\n{PREFIX}free_name\n{PREFIX}free_value\n{PREFIX}free_esd");
    for p in &r.params {
        let _ = writeln!(w, "{} {} {}", quoted(&p.name), exact(p.value), exact(p.esd));
    }
    let n = set.free.len().min(r.correlations.len());
    if n < 2 {
        return;
    }
    let _ = writeln!(w, "loop_\n{PREFIX}correlation_a\n{PREFIX}correlation_b\n{PREFIX}correlation_r");
    for j in 0..n {
        for k in j + 1..n {
            let _ = writeln!(w, "{} {} {}", quoted(&set.free[j].name), quoted(&set.free[k].name), exact(r.correlations[j][k]));
        }
    }
}

/// A string as JSON writes it.
fn json_str(s: &str) -> String {
    let mut o = String::from("\"");
    for ch in s.chars() {
        match ch {
            '"' => o.push_str("\\\""),
            '\\' => o.push_str("\\\\"),
            '\n' => o.push_str("\\n"),
            c if (c as u32) < 0x20 => o.push_str(&format!("\\u{:04x}", c as u32)),
            c => o.push(c),
        }
    }
    o.push('"');
    o
}

/// A number as JSON writes it, to the last bit (null when not finite).
fn json_num(v: f64) -> String {
    if v.is_finite() { format!("{v:?}") } else { "null".into() }
}

/// The record of a run as JSON (E1, `drs refine --trace`): the free parameters with their values at the start, why
/// it stopped, and per cycle χ² at its start and end, Rp, Rwp, Rexp, the Durbin–Watson d, λ, the steps refused (λ,
/// χ²), the singular directions left out, the time, each parameter's value, shift, esd, computed shift and limit,
/// the spectrum of the scaled normal matrix (`Cycle::spectrum`), the covariance and the steering's notes, in words and
/// as data (E15c); the stop's facts and the cycle after which it happened (E15b); the pauses (E15a).
pub fn trace_json(t: &crate::lsq::Trajectory, stop: &crate::lsq::Stop) -> String {
    let list = |v: &[f64]| format!("[{}]", v.iter().map(|x| json_num(*x)).collect::<Vec<_>>().join(", "));
    let mut o = String::from("{\n  \"drs\": ");
    o.push_str(&json_str(env!("CARGO_PKG_VERSION")));
    o.push_str(",\n  \"parameters\": [");
    o.push_str(&t.set.free.iter().map(|p| format!("{{\"name\": {}, \"start\": {}}}", json_str(&p.name), json_num(p.value_in(&t.start)))).collect::<Vec<_>>().join(", "));
    // The stop's facts (E15b), and the cycle after which it happened when the trajectory knows it.
    let facts: String = stop.facts().iter().map(|(k, v)| format!(", {}: {}", json_str(k), json_str(v))).collect();
    let after = t.stop.as_ref().map_or("null".to_string(), |s| s.after.to_string());
    let _ = write!(o, "],\n  \"stop\": {{\"code\": {}, \"text\": {}, \"after\": {after}{facts}}},\n  \"cycles\": [", json_str(stop.code()), json_str(&stop.text()));
    for (k, c) in t.cycles.iter().enumerate() {
        let s = &c.stats;
        let _ = write!(o, "{}\n    {{\"cycle\": {}, \"chi2_start\": {}, \"chi2\": {}, \"rp\": {}, \"rwp\": {}, \"rexp\": {}, \"durbin_watson\": {}, \"lambda\": {}, \"singular\": {}, \"seconds\": {},", if k > 0 { "," } else { "" }, c.cycle, json_num(c.chi2_start), json_num(s.chi2), json_num(s.rp), json_num(s.rwp), json_num(s.rexp), json_num(s.durbin_watson), json_num(c.marquardt), c.singular, json_num(c.seconds));
        let _ = write!(o, "\n     \"tries\": [{}],", c.tries.iter().map(|x| format!("{{\"lambda\": {}, \"chi2\": {}}}", json_num(x.lambda), json_num(x.chi2))).collect::<Vec<_>>().join(", "));
        let _ = write!(o, "\n     \"params\": [{}],", c.params.iter().map(|q| format!("{{\"name\": {}, \"value\": {}, \"shift\": {}, \"esd\": {}, \"computed\": {}, \"at_limit\": {}}}", json_str(&q.name), json_num(q.value), json_num(q.shift), json_num(q.esd), json_num(q.computed), q.at_limit)).collect::<Vec<_>>().join(", "));
        let _ = write!(o, "\n     \"spectrum\": {},", list(&c.spectrum));
        let _ = write!(o, "\n     \"covariance\": [{}],", c.covariance.iter().map(|r| list(r)).collect::<Vec<_>>().join(", "));
        let _ = write!(o, "\n     \"notes\": [{}],", c.notes.iter().map(|n| json_str(n)).collect::<Vec<_>>().join(", "));
        let opt = |x: Option<f64>| x.map_or("null".to_string(), json_num);
        let _ = write!(o, "\n     \"steering\": [{}]}}", c.steering.iter().map(|n| format!("{{\"kind\": {}, \"address\": {}, \"value\": {}, \"lambda\": {}, \"chi2\": {}}}", json_str(n.kind.word()), n.address.as_deref().map_or("null".to_string(), json_str), opt(n.value), opt(n.lambda), opt(n.chi2))).collect::<Vec<_>>().join(", "));
    }
    o.push_str("\n  ],\n  \"pauses\": [");
    o.push_str(&t.pauses.iter().map(|p| format!("{{\"after\": {}, \"why\": {}, \"lambda\": {}, \"restarted\": {}, \"changes\": [{}]}}", p.after, json_str(p.why.word()), json_num(p.lambda), p.restarted.map_or("null".to_string(), json_num), p.changes.iter().map(|c| json_str(&c.words())).collect::<Vec<_>>().join(", "))).collect::<Vec<_>>().join(", "));
    o.push_str("]\n}\n");
    o
}

/// The run's pauses (E15a): after which cycle, why, λ then, the λ the damping restarted at, what was changed (in
/// words); and the changes as data in a second loop.
fn pauses_block(w: &mut String, refinement: Option<&Refinement>) {
    let Some(r) = refinement.filter(|r| !r.pauses.is_empty()) else { return };
    let _ = writeln!(w, "loop_\n{PREFIX}pause_after\n{PREFIX}pause_why\n{PREFIX}pause_lambda\n{PREFIX}pause_restarted\n{PREFIX}pause_changes");
    for p in &r.pauses {
        let changes = if p.changes.is_empty() { ".".to_string() } else { quoted(&p.changes.iter().map(|c| c.words()).collect::<Vec<_>>().join("; ")) };
        let _ = writeln!(w, "{} {} {} {} {changes}", p.after, p.why.word(), num(p.lambda), p.restarted.map_or(".".to_string(), num));
    }
    if r.pauses.iter().any(|p| !p.changes.is_empty()) {
        let _ = writeln!(w, "loop_\n{PREFIX}pause_change_after\n{PREFIX}pause_change_kind\n{PREFIX}pause_change_address\n{PREFIX}pause_change_value");
        for p in &r.pauses {
            for c in &p.changes {
                match c {
                    crate::lsq::Change::Value { address, value } => {
                        let _ = writeln!(w, "{} value {} {}", p.after, quoted(address), num(*value));
                    }
                    crate::lsq::Change::Parameters => {
                        let _ = writeln!(w, "{} parameters . .", p.after);
                    }
                }
            }
        }
    }
}

/// The block names of the phases (as the deposition CIF writes them).
pub(crate) fn phase_blocks(text: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    for b in text.split("\ndata_").skip(1) {
        let name = b.lines().next().unwrap_or("").trim().to_string();
        out.push((name, format!("data_{b}")));
    }
    out
}

/// A Pawley fit's ties (P36, `pawley::ties`): per tie, its number, the phase, and each reflection's indices and the
/// multiplier of its intensity (the first 1; the others the ratio of their multiplicities, equal |F|²). Nothing
/// without ties.
fn pawley_ties_block(w: &mut String, control: &Control) {
    let ties = crate::pawley::ties(&control.model, &control.refinement);
    if ties.is_empty() {
        return;
    }
    let _ = writeln!(w, "loop_\n{PREFIX}pawley_tie_id\n{PREFIX}pawley_tie_phase_id\n{PREFIX}pawley_tie_index_h\n{PREFIX}pawley_tie_index_k\n{PREFIX}pawley_tie_index_l\n{PREFIX}pawley_tie_multiplier");
    for (t, (p, rows)) in ties.iter().enumerate() {
        for (h, c) in rows {
            let _ = writeln!(w, "{} {} {} {} {} {}", t + 1, p + 1, h[0], h[1], h[2], num(*c));
        }
    }
}

/// The results file of a refinement of `control` on `data` (`set` and `refinement` give the standard uncertainties;
/// without them the file holds the values at the control file's state and no uncertainties), and the structure CIFs
/// of its phases ((phase name, text)). `name` is "{control} {time}", the results' id.
pub fn write(control: &Control, data: &Pattern, set: Option<&ParameterSet>, refinement: Option<&Refinement>, name: &str) -> Result<(String, Vec<(String, String)>), String> {
    let m = &control.model;
    let c = m.calculate(&data.x);
    let stats = match refinement {
        Some(r) => r.stats.clone(),
        None => crate::lsq::statistics(m, data, &c.ycalc, set.map_or(0, |s| s.free.len())),
    };
    let has_structure = m.phases.iter().any(|p| matches!(p.intensities, Intensities::Structure(..)));
    let deposition = if has_structure { crate::cif_write::deposition_in(m, data, set, refinement, &control.title, Some(name))? } else { String::new() };
    let blocks = phase_blocks(&deposition);
    let version = env!("CARGO_PKG_VERSION");
    let mut out = String::new();
    let w = &mut out;
    let _ = writeln!(w, "# The results of a refinement, written by diff-ref-struct {version} (drs) as a powder CIF (R28, L092).");
    let _ = writeln!(w, "# Items the IUCr's core and powder dictionaries do not name carry the local prefix {PREFIX}.");
    // The overall block.
    let _ = writeln!(w, "\ndata_overall");
    let _ = writeln!(w, "_pd_block_id                      {}", quoted(&format!("{name}|overall")));
    let _ = writeln!(w, "_audit_creation_method            'diff-ref-struct {version} (drs)'");
    let _ = writeln!(w, "_computing_structure_refinement   'diff-ref-struct {version}'");
    let _ = writeln!(w, "_pd_proc_ls_prof_R_factor         {}", num(stats.rp / 100.0));
    let _ = writeln!(w, "_pd_proc_ls_prof_wR_factor        {}", num(stats.rwp / 100.0));
    let _ = writeln!(w, "_pd_proc_ls_prof_wR_expected      {}", num(stats.rexp / 100.0));
    let _ = writeln!(w, "_refine_ls_goodness_of_fit_all    {}", num(stats.chi2.sqrt()));
    // McCusker et al. (1999), §11: the R factors also with the background taken out (FullProf's "conventional").
    let _ = writeln!(w, "{PREFIX}prof_R_factor_background_free   {}", num(stats.rp_background_free / 100.0));
    let _ = writeln!(w, "{PREFIX}prof_wR_factor_background_free  {}", num(stats.rwp_background_free / 100.0));
    let _ = writeln!(w, "{PREFIX}prof_wR_expected_background_free {}", num(stats.rexp_background_free / 100.0));
    let _ = writeln!(w, "{PREFIX}chi2                         {}", num(stats.chi2));
    let _ = writeln!(w, "{PREFIX}durbin_watson                {}", num(stats.durbin_watson));
    let _ = writeln!(w, "_pd_proc_number_of_points         {}", stats.n);
    let _ = writeln!(w, "_refine_ls_number_parameters      {}", stats.p);
    let _ = writeln!(w, "_refine_ls_number_restraints      {}", m.restraints.len() + set.map_or(0, |s| s.linear.len()));
    if let Some(r) = refinement {
        let _ = writeln!(w, "{PREFIX}cycles                       {}", r.cycles);
        let _ = writeln!(w, "{PREFIX}stop_code                    {}", r.reason.code());
        let _ = writeln!(w, "{PREFIX}stop_text                    {}", quoted(&r.stop));
        // The stop's facts with its cycle (E15b).
        let _ = writeln!(w, "{PREFIX}stop_after                   {}", r.cycles);
        for (k, v) in r.reason.facts() {
            let _ = writeln!(w, "{PREFIX}stop_{k:<24} {}", quoted(&v));
        }
        let _ = writeln!(w, "_refine_ls_shift/su_max           {}", num(r.shift_su_max));
        // Bérar & Lelann's factor for local correlations (P43): the esd's times it estimate the probable errors.
        let _ = writeln!(w, "{PREFIX}berar_factor                 {}", num(r.berar_factor));
    }
    // Each phase with its block, the name the core proposes (P10), its R_Bragg and R_F (P2).
    let _ = writeln!(w, "loop_\n_pd_phase_id\n_pd_phase_block_id\n{PREFIX}phase_name\n{PREFIX}phase_R_Bragg\n{PREFIX}phase_R_F");
    for (p, ph) in m.phases.iter().enumerate() {
        let id = if matches!(ph.intensities, Intensities::Structure(..)) { format!("{name}|phase {}", p + 1) } else { ".".to_string() };
        let (rb, rf) = stats.phases.get(p).copied().unwrap_or((f64::NAN, f64::NAN));
        let _ = writeln!(w, "{} {} {} {} {}", p + 1, quoted(&id), quoted(&crate::cif_write::phase_name(m, p)), num(rb / 100.0), num(rf / 100.0));
    }
    let _ = writeln!(w, "_pd_block_diffractogram_id        {}", quoted(&format!("{name}|pattern")));
    let _ = writeln!(w, "loop_\n{PREFIX}parameter_name\n{PREFIX}parameter_value\n{PREFIX}parameter_su");
    match refinement {
        Some(r) => {
            for f in &r.fields {
                let _ = writeln!(w, "{} {} {}", quoted(&f.0), num(f.1), num(f.2));
            }
        }
        None => {
            let _ = writeln!(w, ". . .");
        }
    }
    cycles_block(w, refinement);
    pauses_block(w, refinement);
    correlations_block(w, set, refinement);
    pawley_ties_block(w, control);
    let _ = writeln!(w, "{PREFIX}control_file\n;\n{}\n;", control.to_text().trim_end());
    // The phases.
    let mut structures = Vec::new();
    let mut k = 0;
    for (p, ph) in m.phases.iter().enumerate() {
        if !matches!(ph.intensities, Intensities::Structure(..)) {
            continue;
        }
        let (bname, text) = &blocks[k];
        k += 1;
        let mut lines = text.lines();
        let head = lines.next().unwrap_or("");
        let body: Vec<&str> = lines.collect();
        let _ = writeln!(w, "\n{head}");
        let _ = writeln!(w, "{}", body.join("\n").trim_end());
        modes_block(w, m, p, set, refinement);
        structure_block(w, m, p, set, refinement);
        structures.push((crate::cif_write::phase_name(m, p), format!("# CIF written by diff-ref-struct {version} (drs): phase {} of {name}\n\n{head}\n{}\n", bname, body.join("\n").trim_end())));
    }
    // The pattern.
    let ids: Vec<usize> = (0..m.phases.len()).collect();
    pattern_block(w, m, data, set, refinement, "pattern", &format!("{name}|pattern"), None, &ids, &|_| {});
    Ok((out, structures))
}

/// A TOF pattern's instrument items (P46): the bank's angle (`_pd_meas_2theta_fixed`), the TOF–d equation with its
/// values (`_pd_calibration_conversion_eqn`), and each TOF value of the bank and of the phases with its standard
/// uncertainty when refined, the resolution file's part apart and the total the calculation takes
/// (`_drs_tof_parameter_*`). The flight path (`_pd_instr_dist_spec/detc`) is not written: neither the `.pcr` nor the
/// FullProf runs' data and resolution files give it (Dtt1 holds it with the angle).
fn tof_items(w: &mut String, m: &crate::powder::Model, b: &crate::tof::Bank, refinement: Option<&Refinement>) {
    use crate::jacobian::Field;
    let _ = writeln!(w, "_pd_meas_2theta_fixed             {}", num(b.two_theta));
    let _ = writeln!(w, "_pd_calibration_conversion_eqn\n;\n{}\n;", b.equation(m.zero));
    let su = |f: Field| refinement.and_then(|r| r.fields.iter().find(|x| x.0 == crate::refinable::address(f, m))).map_or(0.0, |x| x.2);
    let _ = writeln!(w, "loop_\n{PREFIX}tof_parameter_id\n{PREFIX}tof_parameter_phase_id\n{PREFIX}tof_parameter_value\n{PREFIX}tof_parameter_su\n{PREFIX}tof_parameter_instrument\n{PREFIX}tof_parameter_total");
    let mut row = |id: &str, phase: Option<usize>, f: Field, instrument: f64| {
        let v = m.get(f);
        // Values as they are (shortest round trip): Dtt1's eight or nine digits matter at long d.
        let _ = writeln!(w, "{id} {} {v:?} {} {instrument:?} {:?}", phase.map_or(".".to_string(), |p| (p + 1).to_string()), num(su(f)), v + instrument);
    };
    row("zero_shift", None, Field::Zero, b.instrument_part(None));
    let n = if b.thermal.is_some() { 8 } else { 3 };
    for (k, id) in crate::tof::BANK_FIELDS.iter().enumerate().take(n) {
        row(id, None, Field::Tof(k), b.instrument_part(Some(k)));
    }
    for (p, ph) in m.phases.iter().enumerate() {
        let Some(t) = &ph.tof else { continue };
        for (k, id) in crate::tof::PROFILE_FIELDS.iter().enumerate() {
            if t.thermal.is_some() || !(13..=16).contains(&k) {
                row(id, Some(p), Field::TofProfile(p, k), b.profile_instrument_part(k));
            }
        }
    }
}

/// The block of a pattern (`data_{block}`, `_pd_block_id` `id`): `extra` writes items after the id; the radiation,
/// the excluded regions, the profile and the reflections follow. A joint refinement's (P41): `pattern` its pattern
/// (0, 1, …), written before each reflection (`_drs_refln_pattern_id`), and `phase_ids` the joint phase of each of
/// the model's phases (its phase curves' and reflections' phase ids); a single pattern's: None and 0, 1, ….
#[allow(clippy::too_many_arguments)]
pub(crate) fn pattern_block(w: &mut String, m: &crate::powder::Model, data: &Pattern, set: Option<&ParameterSet>, refinement: Option<&Refinement>, block: &str, id: &str, pattern: Option<usize>, phase_ids: &[usize], extra: &dyn Fn(&mut String)) {
    let c = m.calculate(&data.x);
    let w_pts = crate::lsq::weights(m, data, &c.ycalc);
    let _ = writeln!(w, "\ndata_{block}");
    let _ = writeln!(w, "_pd_block_id                      {}", quoted(id));
    extra(w);
    let _ = writeln!(w, "_diffrn_radiation_probe           {}", if m.neutrons { "neutron" } else { "x-ray" });
    let two = m.ratio != 0.0 && (m.lambda[1] - m.lambda[0]).abs() > 1e-9;
    if let Some(b) = &m.tof {
        tof_items(w, m, b, refinement);
    } else {
        let _ = writeln!(w, "_diffrn_radiation_wavelength      {}", num(m.lambda[0]));
        if two {
            let _ = writeln!(w, "loop_\n_diffrn_radiation_wavelength_id\n_diffrn_radiation_wavelength_wt\n{PREFIX}wavelength\n1 1 {}\n2 {} {}", num(m.lambda[0]), num(m.ratio), num(m.lambda[1]));
        }
    }
    let _ = writeln!(w, "_pd_meas_number_of_points         {}", data.x.len());
    let unit = if m.tof.is_some() { "microseconds TOF" } else { "degrees 2theta" };
    let regions: Vec<String> = m.excluded.iter().map(|(a, b)| format!("{} to {} {unit}", num(*a), num(*b))).collect();
    let _ = writeln!(w, "_pd_proc_info_excluded_regions\n;\n{}\n;", if regions.is_empty() { "none".to_string() } else { regions.join("\n") });
    // A file background's source and b1–b4 with su's (FullProf's Nba = 1, eq. 3.6; Publish's P57.2).
    if let Some(f) = &m.background_terms.file {
        let first = m.background_terms.len() - 4;
        let su = |k: usize| refinement.map_or(0.0, |r| {
            let a = crate::refinable::address(crate::jacobian::Field::BackgroundTerm(first + k), m);
            r.fields.iter().find(|x| x.0 == a).map_or(0.0, |x| x.2)
        });
        let _ = writeln!(w, "{PREFIX}bkg_file                      {}", quoted(&f.file));
        for (k, (item, v)) in ["scale", "offset", "stretch", "shift"].iter().zip(f.values()).enumerate() {
            let _ = writeln!(w, "{PREFIX}bkg_file_{item:<8}             {}", crate::cif_write::with_su(v, su(k), 5));
        }
    }
    // The abscissa: 2θ, or a TOF pattern's time of flight in µs (P46).
    let abscissa = if m.tof.is_some() { "_pd_meas_time_of_flight" } else { "_pd_meas_2theta_scan" };
    let _ = writeln!(w, "loop_\n_pd_proc_point_id\n{abscissa}\n_pd_meas_intensity_total\n_pd_proc_ls_weight\n_pd_proc_intensity_bkg_calc\n_pd_calc_intensity_total");
    for id in &phase_ids[..m.phases.len()] {
        let _ = writeln!(w, "{PREFIX}calc_intensity_phase_{}", id + 1);
    }
    // The background by its parts (P49.2), when it has Debye-like terms or profiles from files (P23): the function
    // (with the inverse powers of 2θ), each Debye-like term, each profile; they add up to _pd_proc_intensity_bkg_calc.
    let bt = &m.background_terms;
    let parts = if bt.debye.is_empty() && bt.profiles.is_empty() { Vec::new() } else { m.background_parts(&data.x) };
    for (name, _) in &parts {
        let _ = writeln!(w, "{PREFIX}bkg_calc_{name}");
    }
    for i in 0..data.x.len() {
        let _ = write!(w, "{} {} {} {} {} {}", i + 1, num(data.x[i]), with_sigma(data.y[i], data.sigma[i]), num(w_pts[i]), num(c.background[i]), num(c.ycalc[i]));
        for pc in &c.per_phase {
            let _ = write!(w, " {}", num(pc[i]));
        }
        for (_, v) in &parts {
            let _ = write!(w, " {}", num(v[i]));
        }
        let _ = writeln!(w);
    }
    let rows = crate::diagnostics::reflections(m, data);
    // A Pawley phase's reflections (P36): the intensity measured is the one refined, its share of each wavelength as the
    // calculated one's, and its su from the covariance (the reflection's own, or that of the parameter tying it).
    let two = m.ratio != 0.0 && (m.lambda[1] - m.lambda[0]).abs() > 1e-9;
    let pawley: std::collections::BTreeSet<crate::jacobian::Field> = match (set, refinement) {
        (Some(s), Some(r)) if r.covariance.len() == s.free.len() => s.fields().into_iter().filter(|f| matches!(f, crate::jacobian::Field::Intensity(..))).collect(),
        _ => Default::default(),
    };
    let refined = |r: &crate::diagnostics::ReflectionRow| -> Option<(f64, f64)> {
        let Intensities::List(l) = &m.phases[r.phase].intensities else { return None };
        let f = crate::jacobian::Field::Intensity(r.phase, l.iter().position(|e| e.0 == r.hkl)?);
        if !pawley.contains(&f) {
            return None;
        }
        let share = if !two { 1.0 } else if r.wavelength == 1 { 1.0 / (1.0 + m.ratio) } else { m.ratio / (1.0 + m.ratio) };
        let var = set?.field_covariance(&refinement?.covariance, f, f);
        Some((m.get(f) * share, var.max(0.0).sqrt() * share))
    };
    if pattern.is_some() {
        let _ = writeln!(w, "loop_\n{PREFIX}refln_pattern_id");
    } else {
        let _ = writeln!(w, "loop_");
    }
    // A reflection's position: its 2θ, or its time of flight (µs; the FWHM in µs too) in a TOF pattern (P46).
    let position = if m.tof.is_some() { "refln_tof" } else { "refln_2theta" };
    let _ = writeln!(w, "_pd_refln_phase_id\n_pd_refln_wavelength_id\n_refln_index_h\n_refln_index_k\n_refln_index_l\n_refln_symmetry_multiplicity\n_refln_d_spacing\n{PREFIX}{position}\n{PREFIX}refln_fwhm\n_refln_F_squared_calc\n_refln_intensity_calc\n_refln_intensity_meas\n_refln_intensity_sigma");
    for r in &rows {
        if let Some(k) = pattern {
            let _ = write!(w, "{} ", k + 1);
        }
        let _ = writeln!(
            w,
            "{} {} {} {} {} {} {} {} {} {} {} {} {}",
            phase_ids[r.phase] + 1,
            r.wavelength,
            r.hkl[0],
            r.hkl[1],
            r.hkl[2],
            r.multiplicity,
            num(r.d),
            num(r.two_theta),
            num(r.fwhm),
            r.f2.map_or(".".into(), num),
            num(r.icalc),
            num(refined(r).map_or(r.iobs, |x| x.0)),
            num(refined(r).map_or(r.sigma, |x| x.1))
        );
    }
    if rows.is_empty() {
        let _ = writeln!(w, "{}", vec!["."; if pattern.is_some() { 14 } else { 13 }].join(" "));
    }
}

/// The results file `text` with the user's notes on its reflections (W8, L142; `notes`) in the pattern block, the last
/// one, beside the reflections: a loop of the phase and h k l, or the 2θ, and the text. No notes, no loop.
pub fn with_notes(text: &str, notes: &[crate::notes::Note]) -> String {
    if notes.is_empty() {
        return text.to_string();
    }
    let mut w = text.trim_end().to_string();
    let _ = writeln!(w, "\n\nloop_\n{PREFIX}note_phase_id\n{PREFIX}note_index_h\n{PREFIX}note_index_k\n{PREFIX}note_index_l\n{PREFIX}note_2theta\n{PREFIX}note_text");
    for n in notes {
        let (h, k, l) = n.hkl.map_or((".".to_string(), ".".to_string(), ".".to_string()), |[h, k, l]| (h.to_string(), k.to_string(), l.to_string()));
        let _ = writeln!(w, "{} {h} {k} {l} {} {}", n.phase.map_or(".".into(), |p| p.to_string()), n.two_theta.map_or(".".into(), num), quoted(&n.text.replace(['\n', '\r'], " ")));
    }
    w
}

/// A reflection of a results file: phase (1, 2, …), wavelength, h k l, multiplicity, and d, 2θ, FWHM, |F|² (NaN if
/// none), I calc, I obs, σ.
pub type ReflectionLine = (usize, usize, [i32; 3], usize, [f64; 7]);

/// A results file read back.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Results {
    pub rp: f64,
    pub rwp: f64,
    pub rexp: f64,
    pub chi2: f64,
    pub points_used: usize,
    pub parameters: usize,
    /// Each refined parameter: name, value, standard uncertainty.
    pub refined: Vec<(String, f64, f64)>,
    pub control_file: String,
    /// The phases: name and block id ("." for a phase without a structure).
    pub phases: Vec<(String, String)>,
    pub x: Vec<f64>,
    pub yobs: Vec<f64>,
    /// Yobs's σ, as its su (0 where none is written).
    pub sigma: Vec<f64>,
    pub weight: Vec<f64>,
    pub background: Vec<f64>,
    /// The background by its parts (P49.2, P23): ("smooth", the function with the inverse powers of 2θ), ("debye_1",
    /// a Debye-like term) …, ("file_1", a profile from a file) …, in the file's order; empty when the background has
    /// none of the terms, and in a file written before core 0.63.0.
    pub background_parts: Vec<(String, Vec<f64>)>,
    pub ycalc: Vec<f64>,
    pub per_phase: Vec<Vec<f64>>,
    /// The reflections: phase (1, 2, …), wavelength, h k l, multiplicity, d, 2θ, FWHM, |F|² (NaN if none), I calc,
    /// I obs, σ.
    pub reflections: Vec<ReflectionLine>,
    /// The cycles run (E1): number, χ² at the start, the statistics (Rp, Rwp, Rexp, χ²), λ, singular directions,
    /// time, and each parameter's state; covariance, spectrum, restraints and notes are not kept.
    pub cycles: Vec<crate::lsq::Cycle>,
    /// A plan run's (P37): per cycle of `cycles`, the plan step that ran it, its number (1, 2, …) and name; empty for a
    /// run without a plan and in a file written before core 0.54.0.
    pub cycle_steps: Vec<(usize, String)>,
    /// A Pawley fit's ties (P36): per tie, the phase (0, 1, …) and its reflections with their intensities'
    /// multipliers; empty without ties and in a file written before core 0.55.0.
    pub pawley_ties: Vec<crate::pawley::Tie>,
    /// The user's notes on the reflections (W8, `with_notes`).
    pub notes: Vec<crate::notes::Note>,
    /// Why the run stopped (E11c): `Stop::code` and `Stop::text`; a file written before core 0.32.0 has the text only.
    pub stop_code: Option<String>,
    pub stop_text: Option<String>,
    /// The stop with its facts and its cycle (E15b); none in a file written before core 0.40.0 (or by a run without a
    /// stop).
    pub stop: Option<crate::lsq::Stopped>,
    /// The run's pauses (E15a); empty in a file written before core 0.40.0.
    pub pauses: Vec<crate::lsq::Pause>,
    /// Bérar & Lelann's factor (P43); none in a file written before core 0.42.0.
    pub berar_factor: Option<f64>,
    /// The free parameters at the end: name, value and esd (E11b; empty in a file written before core 0.32.0).
    pub free: Vec<(String, f64, f64)>,
    /// The correlations of the free parameters at the end (E11b), the upper triangle: (a, b, r).
    pub correlations: Vec<(String, String, f64)>,
    /// The profile's abscissa (P46): 2θ (deg) or the time of flight (µs); a reflection's position (`reflections`, the
    /// 2θ of its seven numbers) and FWHM are then in µs too.
    pub axis: crate::formats::Axis,
    /// A TOF pattern's bank angle 2θ (deg), its TOF–d equation as written, and its TOF values (P46); none otherwise
    /// and in a file written before core 0.59.0.
    pub two_theta_fixed: Option<f64>,
    pub conversion: Option<String>,
    pub tof_parameters: Vec<TofParameter>,
}

/// A TOF value of a results file (P46): its dictionary id, the phase (0, 1, …; none for the bank's), the control
/// file's value with its standard uncertainty (0 when not refined), the resolution file's part and the total the
/// calculation takes.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TofParameter {
    pub id: String,
    pub phase: Option<usize>,
    pub value: f64,
    pub su: f64,
    pub instrument: f64,
    pub total: f64,
}

impl Results {
    /// The profile's abscissa in d (Å) or Q (Å⁻¹; `to` is `Axis::D` or `Axis::Q`), with the refined constants of the
    /// model the file holds (its control file): a TOF pattern's T(d) inverted (`tof::Bank::d_spacing`), a constant
    /// wavelength's λ₁ with the zero (`powder::Model::d_spacing_at`). NaN where a point has none (P46).
    pub fn abscissa_in(&self, to: crate::formats::Axis) -> Result<Vec<f64>, String> {
        use crate::formats::Axis;
        let c = crate::control::Control::from_text(&self.control_file).map_err(|d| d.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("; "))?;
        let d: Vec<f64> = self.x.iter().map(|&x| c.model.d_spacing_at(x).unwrap_or(f64::NAN)).collect();
        match to {
            Axis::D => Ok(d),
            Axis::Q => Ok(d.iter().map(|d| 2.0 * std::f64::consts::PI / d).collect()),
            a => Err(format!("the abscissa in {} is not done: d or Q", a.label())),
        }
    }

    /// The correlation matrix of the parameters named, in their order (1 on the diagonal, 0 for a pair the file does
    /// not hold); None if the file holds no correlations.
    pub fn correlation_matrix(&self, names: &[String]) -> Option<Vec<Vec<f64>>> {
        if self.correlations.is_empty() {
            return None;
        }
        let index: std::collections::HashMap<&str, usize> = names.iter().enumerate().map(|(i, n)| (n.as_str(), i)).collect();
        let mut m: Vec<Vec<f64>> = (0..names.len()).map(|j| (0..names.len()).map(|k| if j == k { 1.0 } else { 0.0 }).collect()).collect();
        for (a, b, r) in &self.correlations {
            if let (Some(&j), Some(&k)) = (index.get(a.as_str()), index.get(b.as_str())) {
                m[j][k] = *r;
                m[k][j] = *r;
            }
        }
        Some(m)
    }
}

/// Read a results file (`write`'s).
pub fn read(text: &str) -> Result<Results, String> {
    let blocks = crate::cif::parse(text);
    let find = |id: &str| blocks.iter().find(|b| b.get("_pd_block_id").is_some_and(|v| v.ends_with(id))).ok_or(format!("no block {id}"));
    let o = find("|overall")?;
    let f = |b: &crate::cif::Block, tag: &str| b.get(tag).and_then(crate::cif::number).ok_or(format!("{tag} missing"));
    let col = |b: &crate::cif::Block, tag: &str| -> Result<Vec<f64>, String> { Ok(b.column(tag).ok_or(format!("{tag} missing"))?.iter().map(|v| crate::cif::number(v).unwrap_or(f64::NAN)).collect()) };
    let mut r = Results {
        rp: 100.0 * f(o, "_pd_proc_ls_prof_r_factor")?,
        rwp: 100.0 * f(o, "_pd_proc_ls_prof_wr_factor")?,
        rexp: 100.0 * f(o, "_pd_proc_ls_prof_wr_expected")?,
        chi2: f(o, "_drs_chi2")?,
        points_used: f(o, "_pd_proc_number_of_points")? as usize,
        parameters: f(o, "_refine_ls_number_parameters")? as usize,
        control_file: o.get("_drs_control_file").unwrap_or("").to_string(),
        ..Default::default()
    };
    if let Some(names) = o.column("_drs_parameter_name") {
        let (v, s) = (col(o, "_drs_parameter_value")?, col(o, "_drs_parameter_su")?);
        r.refined = names.iter().zip(v).zip(s).filter(|((n, _), _)| **n != ".").map(|((n, v), s)| (n.to_string(), v, s)).collect();
    }
    if let Some(nums) = o.column("_drs_cycle_number") {
        let c = |t: &str| col(o, &format!("_drs_cycle_{t}"));
        let (c0, c1, rp, rwp, rexp, mu, sing, sec) = (c("chi2_start")?, c("chi2")?, c("r_factor")?, c("wr_factor")?, c("wr_expected")?, c("marquardt")?, c("singular")?, c("seconds")?);
        for i in 0..nums.len() {
            let stats = crate::lsq::Statistics { n: 0, p: 0, c: 0, rp: rp[i], rwp: rwp[i], rexp: rexp[i], chi2: c1[i], durbin_watson: 0.0, phases: Vec::new(), rp_background_free: 0.0, rwp_background_free: 0.0, rexp_background_free: 0.0 };
            r.cycles.push(crate::lsq::Cycle { cycle: nums[i].parse().map_err(|_| format!("_drs_cycle_number {}", nums[i]))?, params: Vec::new(), ycalc: Vec::new(), stats, chi2_start: c0[i], marquardt: mu[i], tries: Vec::new(), singular: sing[i] as usize, restraints: Vec::new(), notes: Vec::new(), steering: Vec::new(), covariance: Vec::new(), spectrum: Vec::new(), seconds: sec[i] });
        }
        if let (Some(k), Some(names)) = (o.column("_drs_cycle_step"), o.column("_drs_cycle_step_name")) {
            for (k, n) in k.iter().zip(names) {
                r.cycle_steps.push((k.parse().map_err(|_| format!("_drs_cycle_step {k}"))?, n.to_string()));
            }
        }
        if let Some(cyc) = o.column("_drs_cycle_parameter_cycle") {
            let names = o.column("_drs_cycle_parameter_name").ok_or("_drs_cycle_parameter_name missing")?;
            let limit = o.column("_drs_cycle_parameter_at_limit").ok_or("_drs_cycle_parameter_at_limit missing")?;
            let c = |t: &str| col(o, &format!("_drs_cycle_parameter_{t}"));
            let (v, sh, e, cp) = (c("value")?, c("shift")?, c("esd")?, c("computed")?);
            for i in 0..cyc.len() {
                let k: usize = cyc[i].parse().map_err(|_| format!("_drs_cycle_parameter_cycle {}", cyc[i]))?;
                let target = r.cycles.iter_mut().find(|x| x.cycle == k).ok_or(format!("a parameter of cycle {k}, which is not in the cycles' loop"))?;
                target.params.push(crate::lsq::ParamState { name: names[i].to_string(), value: v[i], shift: sh[i], esd: e[i], shift_esd: if e[i] > 0.0 { sh[i] / e[i] } else { 0.0 }, computed: cp[i], at_limit: limit[i] == "yes" });
            }
        }
        if let Some(cyc) = o.column("_drs_steering_cycle") {
            let col_of = |t: &str| o.column(&format!("_drs_steering_{t}")).ok_or(format!("_drs_steering_{t} missing"));
            let (kind, address, value, lambda, chi2) = (col_of("kind")?, col_of("address")?, col_of("value")?, col_of("lambda")?, col_of("chi2")?);
            for i in 0..cyc.len() {
                let k: usize = cyc[i].parse().map_err(|_| format!("_drs_steering_cycle {}", cyc[i]))?;
                let kind = crate::lsq::NoteKind::from_word(kind[i]).ok_or(format!("_drs_steering_kind {}", kind[i]))?;
                let address = (address[i] != "." && address[i] != "?").then(|| address[i].to_string());
                let note = crate::lsq::SteeringNote { kind, address, value: crate::cif::number(value[i]), lambda: crate::cif::number(lambda[i]), chi2: crate::cif::number(chi2[i]) };
                r.cycles.iter_mut().find(|x| x.cycle == k).ok_or(format!("a steering note of cycle {k}, which is not in the cycles' loop"))?.steering.push(note);
            }
        }
    }
    if let Some(ids) = o.column("_drs_pawley_tie_id") {
        let int = |t: &str| -> Result<Vec<i64>, String> { o.column(t).ok_or(format!("{t} missing"))?.iter().map(|v| v.parse::<i64>().map_err(|_| format!("{t} {v}"))).collect() };
        let (ph, h, k, l, c) = (int("_drs_pawley_tie_phase_id")?, int("_drs_pawley_tie_index_h")?, int("_drs_pawley_tie_index_k")?, int("_drs_pawley_tie_index_l")?, col(o, "_drs_pawley_tie_multiplier")?);
        let mut last = None;
        for i in 0..ids.len() {
            if last != Some(ids[i]) {
                r.pawley_ties.push(((ph[i] - 1).max(0) as usize, Vec::new()));
                last = Some(ids[i]);
            }
            r.pawley_ties.last_mut().unwrap().1.push(([h[i] as i32, k[i] as i32, l[i] as i32], c[i]));
        }
    }
    r.stop_code = o.get("_drs_stop_code").map(str::to_string);
    r.berar_factor = o.get("_drs_berar_factor").and_then(crate::cif::number);
    r.stop_text = o.get("_drs_stop_text").or(o.get("_drs_stop")).map(str::to_string);
    if let (Some(code), Some(after)) = (&r.stop_code, o.get("_drs_stop_after").and_then(|v| v.parse::<usize>().ok())) {
        r.stop = crate::lsq::Stop::from_facts(code, |k| o.get(&format!("_drs_stop_{k}")).map(str::to_string)).map(|stop| crate::lsq::Stopped { stop, after });
    }
    if let Some(after) = o.column("_drs_pause_after") {
        let why = o.column("_drs_pause_why").ok_or("_drs_pause_why missing")?;
        let restarted = o.column("_drs_pause_restarted").ok_or("_drs_pause_restarted missing")?;
        let lambda = col(o, "_drs_pause_lambda")?;
        for i in 0..after.len() {
            let why = match why[i] {
                "step" => crate::lsq::PauseWhy::Step,
                "caller" => crate::lsq::PauseWhy::Caller,
                w => return Err(format!("_drs_pause_why {w}")),
            };
            r.pauses.push(crate::lsq::Pause { after: after[i].parse().map_err(|_| format!("_drs_pause_after {}", after[i]))?, why, lambda: lambda[i], changes: Vec::new(), restarted: crate::cif::number(restarted[i]) });
        }
        if let Some(at) = o.column("_drs_pause_change_after") {
            let kind = o.column("_drs_pause_change_kind").ok_or("_drs_pause_change_kind missing")?;
            let address = o.column("_drs_pause_change_address").ok_or("_drs_pause_change_address missing")?;
            let value = col(o, "_drs_pause_change_value")?;
            for i in 0..at.len() {
                let k: usize = at[i].parse().map_err(|_| format!("_drs_pause_change_after {}", at[i]))?;
                let c = match kind[i] {
                    "value" => crate::lsq::Change::Value { address: address[i].to_string(), value: value[i] },
                    "parameters" => crate::lsq::Change::Parameters,
                    w => return Err(format!("_drs_pause_change_kind {w}")),
                };
                r.pauses.iter_mut().find(|p| p.after == k).ok_or(format!("a change after cycle {k}, which has no pause"))?.changes.push(c);
            }
        }
    }
    if let Some(names) = o.column("_drs_free_name") {
        let (v, e) = (col(o, "_drs_free_value")?, col(o, "_drs_free_esd")?);
        r.free = names.iter().zip(v).zip(e).map(|((n, v), e)| (n.to_string(), v, e)).collect();
    }
    if let Some(a) = o.column("_drs_correlation_a") {
        let b = o.column("_drs_correlation_b").ok_or("_drs_correlation_b missing")?;
        let v = col(o, "_drs_correlation_r")?;
        r.correlations = a.iter().zip(b).zip(v).map(|((a, b), v)| (a.to_string(), b.to_string(), v)).collect();
    }
    let names = o.column("_drs_phase_name").unwrap_or_default();
    let ids = o.column("_pd_phase_block_id").unwrap_or_default();
    r.phases = names.iter().zip(ids).map(|(n, i)| (n.to_string(), i.to_string())).collect();
    let p = find("|pattern")?;
    // The abscissa: 2θ or a TOF pattern's time of flight (P46).
    (r.x, r.axis) = match col(p, "_pd_meas_time_of_flight") {
        Ok(x) => (x, crate::formats::Axis::TimeOfFlight),
        Err(_) => (col(p, "_pd_meas_2theta_scan")?, crate::formats::Axis::TwoTheta),
    };
    if r.axis == crate::formats::Axis::TimeOfFlight {
        r.two_theta_fixed = p.get("_pd_meas_2theta_fixed").and_then(crate::cif::number);
        r.conversion = p.get("_pd_calibration_conversion_eqn").map(|s| s.trim().to_string());
        if let Some(ids) = p.column(&format!("{PREFIX}tof_parameter_id")) {
            let phases = p.column(&format!("{PREFIX}tof_parameter_phase_id")).ok_or("_drs_tof_parameter_phase_id missing")?;
            let c = |t: &str| col(p, &format!("{PREFIX}tof_parameter_{t}"));
            let (v, su, ins, tot) = (c("value")?, c("su")?, c("instrument")?, c("total")?);
            for i in 0..ids.len() {
                r.tof_parameters.push(TofParameter { id: ids[i].to_string(), phase: phases[i].parse::<usize>().ok().and_then(|p| p.checked_sub(1)), value: v[i], su: su[i], instrument: ins[i], total: tot[i] });
            }
        }
    }
    let y = p.column("_pd_meas_intensity_total").ok_or("_pd_meas_intensity_total missing")?;
    (r.yobs, r.sigma) = y.iter().map(|v| crate::cif::value_su(v).unwrap_or((f64::NAN, 0.0))).unzip();
    r.weight = col(p, "_pd_proc_ls_weight")?;
    r.background = col(p, "_pd_proc_intensity_bkg_calc")?;
    r.ycalc = col(p, "_pd_calc_intensity_total")?;
    for k in 1..=r.phases.len() {
        r.per_phase.push(col(p, &format!("_drs_calc_intensity_phase_{k}"))?);
    }
    // The background's parts (P49.2), in the order written.
    let tags: Vec<String> = p.loops.iter().find(|(t, _)| t.iter().any(|x| x == "_pd_proc_point_id")).map_or(Vec::new(), |(t, _)| t.clone());
    for name in tags.iter().filter_map(|n| n.strip_prefix(&format!("{PREFIX}bkg_calc_")).map(str::to_string)) {
        r.background_parts.push((name.clone(), col(p, &format!("{PREFIX}bkg_calc_{name}"))?));
    }
    if let Some(texts) = p.column(&format!("{PREFIX}note_text")) {
        let c = |t: &str| p.column(&format!("{PREFIX}note_{t}")).unwrap_or_default();
        let (ph, h, k, l, tt) = (c("phase_id"), c("index_h"), c("index_k"), c("index_l"), c("2theta"));
        for i in 0..texts.len() {
            let int = |v: &[&str]| v.get(i).and_then(|x| x.parse::<i32>().ok());
            let hkl = match (int(&h), int(&k), int(&l)) {
                (Some(a), Some(b), Some(c)) => Some([a, b, c]),
                _ => None,
            };
            r.notes.push(crate::notes::Note { phase: int(&ph).map(|p| p as usize), hkl, two_theta: tt.get(i).and_then(|v| crate::cif::number(v)), text: texts[i].to_string() });
        }
    }
    if let Some(ph) = p.column("_pd_refln_phase_id") {
        let ints = |tag: &str| -> Result<Vec<i32>, String> { Ok(p.column(tag).ok_or(format!("{tag} missing"))?.iter().map(|v| v.parse().unwrap_or(0)).collect()) };
        let (wl, h, k, l, mult) = (ints("_pd_refln_wavelength_id")?, ints("_refln_index_h")?, ints("_refln_index_k")?, ints("_refln_index_l")?, ints("_refln_symmetry_multiplicity")?);
        let position = if r.axis == crate::formats::Axis::TimeOfFlight { "_drs_refln_tof" } else { "_drs_refln_2theta" };
        let cols: Vec<Vec<f64>> = ["_refln_d_spacing", position, "_drs_refln_fwhm", "_refln_f_squared_calc", "_refln_intensity_calc", "_refln_intensity_meas", "_refln_intensity_sigma"].iter().map(|t| col(p, t)).collect::<Result<_, _>>()?;
        for i in 0..ph.len() {
            if ph[i] == "." {
                continue;
            }
            r.reflections.push((ph[i].parse().unwrap_or(0), wl[i] as usize, [h[i], k[i], l[i]], mult[i] as usize, std::array::from_fn(|j| cols[j][i])));
        }
    }
    Ok(r)
}

/// The reflections as CSV (the file dictionary's `export_reflections`, written on demand): one line per reflection.
pub fn reflections_csv(rows: &[crate::diagnostics::ReflectionRow]) -> String {
    let mut s = String::from("phase,wavelength,h,k,l,multiplicity,two_theta,d,fwhm,eta,F2,I_calc,I_obs,sigma_I_obs,residual\n");
    for r in rows {
        let _ = writeln!(
            s,
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            r.phase + 1, r.wavelength, r.hkl[0], r.hkl[1], r.hkl[2], r.multiplicity, num(r.two_theta), num(r.d), num(r.fwhm), num(r.eta),
            r.f2.map_or(String::new(), num), num(r.icalc), num(r.iobs), num(r.sigma), num(r.residual)
        );
    }
    s
}

/// The profile as CSV (`export_profile`, written on demand): 2θ, Yobs, σ, Ycal, the background, each phase's curve,
/// the weighted difference and the cumulative χ² as a fraction of the whole.
pub fn profile_csv(model: &crate::powder::Model, data: &Pattern) -> String {
    let c = model.calculate(&data.x);
    let cum = crate::diagnostics::cumulative(model, data, &c.ycalc);
    let frac = cum.fraction();
    let mut s = String::from("two_theta,y_obs,sigma,y_calc,background");
    for p in 1..=model.phases.len() {
        let _ = write!(s, ",phase_{p}");
    }
    s.push_str(",weighted_difference,cumulative_chi2_fraction\n");
    for i in 0..data.x.len() {
        let _ = write!(s, "{},{},{},{},{}", num(data.x[i]), num(data.y[i]), num(data.sigma[i]), num(c.ycalc[i]), num(c.background[i]));
        for pc in &c.per_phase {
            let _ = write!(s, ",{}", num(pc[i]));
        }
        let _ = writeln!(s, ",{},{}", num(cum.weighted_difference[i]), num(frac[i]));
    }
    s
}
