//! A CIF of a refinement, for deposition (P12, R-list "outputs of phase 2"; L081): one data block per structure phase
//! with what the IUCr's core and powder dictionaries ask of a Rietveld refinement, the numbers with their standard
//! uncertainties from the engine's covariance, and "?" where only the author can say (as FullProf's own CIF, "a
//! template with filled structural items", leaves them).
//!
//! - Chemistry: the sum formula and Z from the atoms' chemical occupancies (`bonds::chemical_occupancies`) times their
//!   multiplicities: Z is the largest integer that divides every element's count in the cell.
//! - Symmetry: the space group's number, Hermann–Mauguin and Hall symbols, and every operator; the site-symmetry codes
//!   of the geometry ("n_klm") number the operators in that list.
//! - Cell: a, b, c, α, β, γ with their uncertainties, the volume with its own (first order, from those of the cell).
//! - Experiment and refinement: radiation and wavelength, the 2θ range and step of the data, Rp, Rwp, Rexp, the
//!   goodness of fit √χ², R_Bragg, the numbers of points, parameters, restraints and reflections.
//! - Atoms: label, type, x, y, z, U_iso (B/8π²), chemical occupancy, multiplicity and Wyckoff letter; U_ij where the
//!   displacements are anisotropic (U_ij = β_ij/(2π² a*_i a*_j)).
//! - Geometry: the bonds (cation–anion pairs Bond_Str counts, or every distance up to 3 Å when the species carry no
//!   charge) and the angles between the bonds of each atom, with their uncertainties (`bonds`).
//!
//! Uncertainties are written by the IUCr's rule: two digits when the first two are 19 or less, one otherwise.

use crate::bonds;
use crate::data::Pattern;
use crate::jacobian::Field;
use crate::lsq::{Refinement, Statistics};
use crate::powder::{Intensities, Model};
use crate::refinable::ParameterSet;
use crate::site::Site;
use std::f64::consts::PI;
use std::fmt::Write;

/// A value with its standard uncertainty as a CIF writes it: "5.74819(5)", "0.0118(2)", "90"; without one, the value
/// to `places` decimals with the trailing zeros dropped. The one rule of the crate for an su (P6, Josu's decision,
/// L116): the IUCr's Statistical descriptors (Rec. 12: "in general, uncertainties should be rounded up"; Schwarzenbach
/// et al. 1995, Acta Cryst. A51, 565) with its range: the su is written in units of the value's last digit as 2 to 19
/// (two digits when the leading ones are 10 to 19, one otherwise), **rounded up**, and the value rounded half away
/// from zero to that digit: 2.2865 ± 0.0032 is 2.287(4), 5.680030 ± 0.000103 is 5.68003(11).
pub fn with_su(v: f64, su: f64, places: usize) -> String {
    if su.is_nan() || su <= 0.0 || !su.is_finite() {
        let s = format!("{v:.places$}");
        let s = if s.contains('.') { s.trim_end_matches('0').trim_end_matches('.').to_string() } else { s };
        return if s == "-0" { "0".into() } else { s };
    }
    // Rounded up, the binary representation's last bits not counting (0.004 is 4, not 5).
    let up = |x: f64| (x * (1.0 - 1e-12)).ceil().max(1.0) as i64;
    // The decimal k at which su is 1.0 to 9.99… units of 10⁻ᵏ: below 2 it is written with two digits at the next
    // decimal (10 to 20; 20 is 2 at this one), from 2 with one (2 to 10; 10 is written as two digits, "10").
    let mut k = (-su.log10().floor()) as i32;
    let x = su * 10f64.powi(k);
    let units = if x < 2.0 {
        k += 1;
        match up(su * 10f64.powi(k)) {
            u if u >= 20 => {
                k -= 1;
                2
            }
            u => u,
        }
    } else {
        up(x)
    };
    // The value rounded half away from zero to 10⁻ᵏ (a half written in decimal stays a half in binary's last bits).
    let scaled = v * 10f64.powi(k);
    let r = scaled.signum() * (scaled.abs() + 0.5 + 1e-9 * scaled.abs().max(1.0)).floor();
    if k <= 0 {
        let unit = 10f64.powi(-k);
        return format!("{}({})", (r * unit) as i64, units * unit as i64);
    }
    let s = format!("{:.*}", k as usize, r / 10f64.powi(k));
    let s = if s.starts_with("-0") && r == 0.0 { s[1..].to_string() } else { s };
    format!("{s}({units})")
}

/// The uncertainty of each field: from a refinement's covariance (the fields of `set`, in the order of
/// `Refinement::fields`); 0 for a field not refined.
pub struct Uncertainty {
    fields: Vec<Field>,
    values: Vec<f64>,
}

impl Uncertainty {
    pub fn new(set: Option<&ParameterSet>, refinement: Option<&Refinement>) -> Uncertainty {
        match (set, refinement) {
            (Some(s), Some(r)) => Uncertainty { fields: s.fields(), values: r.fields.iter().map(|f| f.2).collect() },
            _ => Uncertainty { fields: Vec::new(), values: Vec::new() },
        }
    }

    pub fn of(&self, f: Field) -> f64 {
        self.fields.iter().position(|g| *g == f).and_then(|k| self.values.get(k).copied()).unwrap_or(0.0)
    }
}

fn element(species: &str) -> String {
    let s: String = species.chars().take_while(|c| c.is_ascii_alphabetic()).collect();
    let mut c = s.chars();
    match c.next() {
        Some(f) => f.to_ascii_uppercase().to_string() + &c.as_str().to_lowercase(),
        None => species.to_string(),
    }
}

/// Text quoted for CIF when it holds blanks or starts with a reserved character.
fn quoted(s: &str) -> String {
    if s.is_empty() {
        "?".into()
    } else if s.contains(char::is_whitespace) || s.starts_with(['_', '#', '$', '\'', '"', '[', ']', ';']) {
        format!("'{}'", s.replace('\'', "’"))
    } else {
        s.to_string()
    }
}

/// The sum formula ("La Mn O3") and Z of a structure, from its chemical occupancies.
pub fn formula(st: &crate::structure_factor::Structure) -> Result<(String, usize), String> {
    let occ = bonds::chemical_occupancies(st)?;
    let mut counts: Vec<(String, f64)> = Vec::new();
    for (a, o) in st.atoms.iter().zip(&occ) {
        let m = Site::new(&st.group, a.x)?.multiplicity as f64;
        let el = element(&a.species);
        match counts.iter_mut().find(|c| c.0 == el) {
            Some(c) => c.1 += o * m,
            None => counts.push((el, o * m)),
        }
    }
    // Hill order: C, H, then alphabetical; without C, all alphabetical.
    let has_c = counts.iter().any(|c| c.0 == "C");
    let rank = |e: &str| match e {
        "C" if has_c => (0, String::new()),
        "H" if has_c => (1, String::new()),
        _ => (2, e.to_string()),
    };
    counts.sort_by_key(|c| rank(&c.0));
    // Z: the largest integer that leaves every count in the cell an integer.
    let least = counts.iter().map(|c| c.1).fold(f64::INFINITY, f64::min).round().max(1.0) as usize;
    let z = (1..=least).rev().find(|z| counts.iter().all(|c| {
        let q = c.1 / *z as f64;
        (q - q.round()).abs() < 1e-3
    })).unwrap_or(1);
    let parts: Vec<String> = counts.iter().map(|(e, n)| {
        let q = n / z as f64;
        let s = with_su(q, 0.0, 4);
        if s == "1" { e.clone() } else { format!("{e}{s}") }
    }).collect();
    Ok((parts.join(" "), z))
}

/// The volume's uncertainty from the cell's (first order, independent).
fn volume_su(c: &crate::cell::Cell, su: &[f64; 6]) -> f64 {
    let v = |p: &[f64; 6]| crate::cell::Cell { a: p[0], b: p[1], c: p[2], alpha: p[3], beta: p[4], gamma: p[5] }.volume();
    let p0 = [c.a, c.b, c.c, c.alpha, c.beta, c.gamma];
    let mut var = 0.0;
    for k in 0..6 {
        if su[k] > 0.0 {
            let h = 1e-6 * p0[k].abs().max(1.0);
            let (mut p, mut m) = (p0, p0);
            p[k] += h;
            m[k] -= h;
            var += ((v(&p) - v(&m)) / (2.0 * h) * su[k]).powi(2);
        }
    }
    var.sqrt()
}

/// The CIF of a refinement: `model` at its refined values, its `data`, and the refinement (its statistics and
/// covariance; without one the statistics are calculated and no uncertainties are written).
pub fn deposition(model: &Model, data: &Pattern, set: Option<&ParameterSet>, refinement: Option<&Refinement>, title: &str) -> Result<String, String> {
    deposition_in(model, data, set, refinement, title, None)
}

/// As `deposition`, its blocks named for the results file `results` ("{control} {time}", R28): each phase block's
/// `_pd_block_id` is "{results}|phase n" and its `_pd_block_diffractogram_id` "{results}|pattern" (P8). Without one,
/// the ids are the phase's name and "pattern".
pub fn deposition_in(model: &Model, data: &Pattern, set: Option<&ParameterSet>, refinement: Option<&Refinement>, title: &str, results: Option<&str>) -> Result<String, String> {
    let u = Uncertainty::new(set, refinement);
    let prop = match (set, refinement) {
        (Some(s), Some(r)) => Some(crate::geometry::Propagator::new(model, s, r)),
        _ => None,
    };
    let stats: Statistics = match refinement {
        Some(r) => r.stats.clone(),
        None => crate::lsq::statistics(model, data, &model.calculate(&data.x).ycalc, 0),
    };
    let mut out = String::new();
    let w = &mut out;
    let _ = writeln!(w, "# CIF written by diff-ref-struct {} (drs) from a Rietveld refinement", env!("CARGO_PKG_VERSION"));
    let _ = writeln!(w, "# Items marked ? are for the authors to give before deposition.");
    let n_phases = model.phases.iter().filter(|p| matches!(p.intensities, Intensities::Structure(..))).count();
    if n_phases == 0 {
        return Err("no structure phase to describe".into());
    }
    for (p, ph) in model.phases.iter().enumerate() {
        let Intensities::Structure(st0, _) = &ph.intensities else { continue };
        let mut st = (**st0).clone();
        st.cell = ph.cell;
        let name = phase_name(model, p);
        let _ = writeln!(w, "\ndata_{name}");
        let (block, pattern) = match results {
            Some(r) => (format!("{r}|phase {}", p + 1), format!("{r}|pattern")),
            None => (name.clone(), "pattern".to_string()),
        };
        let _ = writeln!(w, "_pd_block_id                      {}", quoted(&block));
        let _ = writeln!(w, "_pd_block_diffractogram_id        {}", quoted(&pattern));
        let _ = writeln!(w, "_pd_phase_name                    {}", quoted(&name));
        let _ = writeln!(w, "_audit_creation_method            'diff-ref-struct {} (drs)'", env!("CARGO_PKG_VERSION"));
        let _ = writeln!(w, "_computing_structure_refinement   'diff-ref-struct {}'", env!("CARGO_PKG_VERSION"));
        let _ = writeln!(w, "_computing_publication_material   'diff-ref-struct {}'", env!("CARGO_PKG_VERSION"));
        let _ = writeln!(w, "_computing_data_collection        ?\n_computing_structure_solution     ?\n_computing_molecular_graphics     ?");
        let _ = writeln!(w, "_publ_section_title\n;\n{}\n;", if title.trim().is_empty() { "?" } else { title.trim() });
        let _ = writeln!(w, "_publ_contact_author_name         ?\n_publ_contact_author_email        ?\n_journal_name_full                ?");
        // Chemistry.
        let (sum, z) = formula(&st)?;
        let _ = writeln!(w, "_chemical_name_systematic         ?");
        let _ = writeln!(w, "_chemical_formula_sum             '{sum}'");
        let _ = writeln!(w, "_cell_formula_units_Z             {z}");
        // Symmetry.
        let g = &st.group;
        if let Some(s) = &g.setting {
            let _ = writeln!(w, "_space_group_IT_number            {}", s.number);
            let _ = writeln!(w, "_space_group_name_H-M_alt         '{}'", s.hm.replace('_', ""));
        }
        let _ = writeln!(w, "_space_group_name_Hall            '{}'", g.hall);
        let _ = writeln!(w, "_space_group_crystal_system       {}", crystal_system(g.point_group));
        let _ = writeln!(w, "loop_\n_space_group_symop_id\n_space_group_symop_operation_xyz");
        for (k, op) in g.ops.iter().enumerate() {
            let _ = writeln!(w, "{} '{}'", k + 1, op.to_xyz());
        }
        // Cell.
        let c = &st.cell;
        let su: [f64; 6] = std::array::from_fn(|k| u.of(Field::Cell(p, k)));
        let names = ["length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma"];
        let vals = [c.a, c.b, c.c, c.alpha, c.beta, c.gamma];
        for k in 0..6 {
            let _ = writeln!(w, "_cell_{:<28} {}", names[k], with_su(vals[k], su[k], if k < 3 { 5 } else { 4 }));
        }
        // The volume's su through the covariance, the cell's correlations and ties included (P43; as the series
        // table's, `ParameterSet::cell_volume`); from the cell's su's alone only without one.
        let vsu = match (set, refinement) {
            (Some(s), Some(r)) if r.covariance.len() == s.free.len() && !s.free.is_empty() => s.cell_volume(model, p, &r.covariance).1,
            _ => volume_su(c, &su),
        };
        let _ = writeln!(w, "_cell_volume                      {}", with_su(c.volume(), vsu, 3));
        let _ = writeln!(w, "_cell_measurement_temperature     ?");
        // Experiment.
        let _ = writeln!(w, "_diffrn_radiation_probe           {}", if model.neutrons { "neutron" } else { "x-ray" });
        // A TOF pattern (P05) has no wavelength: its bank's angle instead.
        match &model.tof {
            Some(b) => {
                let _ = writeln!(w, "_pd_meas_2theta_fixed             {}", with_su(b.two_theta, 0.0, 4));
            }
            None => {
                let _ = writeln!(w, "_diffrn_radiation_wavelength      {}", with_su(model.lambda[0], 0.0, 6));
            }
        }
        if model.tof.is_none() && !model.neutrons && model.lambda[1] > 0.0 && (model.lambda[1] - model.lambda[0]).abs() > 1e-6 {
            let _ = writeln!(w, "_diffrn_radiation_type            'K\\a~1~ and K\\a~2~ ({}, {}; ratio {})'", model.lambda[0], model.lambda[1], with_su(model.ratio, 0.0, 4));
        }
        let _ = writeln!(w, "_diffrn_ambient_temperature       ?\n_diffrn_measurement_device_type   ?");
        if let (Some(x0), Some(x1), None) = (data.x.first(), data.x.last(), &model.tof) {
            let step = if data.x.len() > 1 { (x1 - x0) / (data.x.len() - 1) as f64 } else { 0.0 };
            let _ = writeln!(w, "_pd_meas_number_of_points         {}", data.x.len());
            let _ = writeln!(w, "_pd_meas_2theta_range_min         {}", with_su(*x0, 0.0, 4));
            let _ = writeln!(w, "_pd_meas_2theta_range_max         {}", with_su(*x1, 0.0, 4));
            let _ = writeln!(w, "_pd_meas_2theta_range_inc         {}", with_su(step, 0.0, 5));
        }
        // Refinement.
        // The reflections down to the shortest d the data reach: a TOF pattern's first time of flight, a constant
        // wavelength's last angle.
        let dmin = match &model.tof {
            Some(_) => data.x.first().and_then(|&t| model.d_spacing_at(t)).unwrap_or(f64::INFINITY),
            None => model.lambda[0] / (2.0 * (data.x.last().copied().unwrap_or(180.0).min(179.9) / 2.0).to_radians().sin()),
        };
        let reflections = crate::reflections::generate(g, &st.cell, dmin, 1e9, None).len();
        let _ = writeln!(w, "_refine_ls_structure_factor_coef  Inet");
        let _ = writeln!(w, "_pd_proc_ls_prof_R_factor         {}", with_su(stats.rp / 100.0, 0.0, 4));
        let _ = writeln!(w, "_pd_proc_ls_prof_wR_factor        {}", with_su(stats.rwp / 100.0, 0.0, 4));
        let _ = writeln!(w, "_pd_proc_ls_prof_wR_expected      {}", with_su(stats.rexp / 100.0, 0.0, 4));
        let _ = writeln!(w, "_refine_ls_goodness_of_fit_all    {}", with_su(stats.chi2.sqrt(), 0.0, 3));
        if let Some(&(rb, rf)) = stats.phases.get(p) {
            let _ = writeln!(w, "_refine_ls_R_I_factor             {}", with_su(rb / 100.0, 0.0, 4));
            let _ = writeln!(w, "_refine_ls_R_factor_all           {}", with_su(rf / 100.0, 0.0, 4));
        }
        if let Some(rf) = r_fsqd(model, data, p) {
            let _ = writeln!(w, "_refine_ls_R_Fsqd_factor          {}", with_su(rf, 0.0, 4));
        }
        let _ = writeln!(w, "_refine_ls_number_parameters      {}", set.map_or(0, |s| s.free.len()));
        let _ = writeln!(w, "_refine_ls_number_restraints      {}", model.restraints.iter().filter(|r| r.phase == p).count());
        let _ = writeln!(w, "_refine_ls_number_constraints     {}", set.map_or(0, constraints));
        let _ = writeln!(w, "_refine_ls_number_reflns          {reflections}");
        let _ = writeln!(w, "_refine_ls_weighting_scheme       {}", match model.weighting {
            crate::powder::Weighting::Variance => "sigma",
            crate::powder::Weighting::MaximumLikelihood => "calc",
            crate::powder::Weighting::Unit => "unit",
        });
        let _ = writeln!(w, "_refine_ls_weighting_details      {}", quoted(match model.weighting {
            crate::powder::Weighting::Variance => "w = 1/sigma^2^(y~obs~)",
            crate::powder::Weighting::MaximumLikelihood => "w = 1/sigma^2^(y~calc~), the variance recalculated from y~calc~ at each cycle",
            crate::powder::Weighting::Unit => "w = 1",
        }));
        if let Some(r) = refinement {
            let _ = writeln!(w, "_refine_ls_shift/su_max           {}", with_su(r.shift_su_max, 0.0, 3));
        }
        let _ = writeln!(w, "_refine_special_details\n;\n{}\n;", special_details(stats.chi2));
        let _ = writeln!(w, "_pd_proc_ls_profile_function      {}", quoted(&profile_function(model, ph)));
        let _ = writeln!(w, "_pd_proc_ls_background_function   {}", quoted(&background_function(model)));
        let _ = writeln!(w, "_pd_proc_ls_pref_orient_corr      {}", quoted(&preferred_orientation(model, ph)));
        // Atoms (P3, P8): labels unique, positions reduced to the cell with su's through the covariance (a mode
        // phase's through its amplitudes'), none on a value the symmetry fixes (P7).
        let geo = crate::geometry::geometry(model, p, prop.as_ref(), crate::geometry::BondRule::Default)?;
        let labels = &geo.labels;
        let st = &geo.structure;
        let occ = bonds::chemical_occupancies(st)?;
        let ueq = crate::geometry::u_eq(model, p, prop.as_ref());
        // A rigid-body phase's atoms flagged G, their positions those of their group (P41, Publish's P47.3).
        let grouped = |a: usize| ph.rigid.iter().any(|g| g.members.iter().any(|m| m.index == a));
        let flags = if ph.rigid.is_empty() { "" } else { "\n_atom_site_refinement_flags_posn" };
        let _ = writeln!(w, "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_U_iso_or_equiv\n_atom_site_adp_type\n_atom_site_occupancy\n_atom_site_site_symmetry_multiplicity\n_atom_site_Wyckoff_symbol{flags}");
        for (a, at) in st.atoms.iter().enumerate() {
            let site = Site::new(g, at.x)?;
            let xs: Vec<String> = (0..3).map(|k| with_su(geo.positions[a][k].0, geo.positions[a][k].1, 5)).collect();
            // U_iso from B, or U_eq for an anisotropic atom (P15).
            let (uiso, usu) = match ueq.get(a).copied().flatten() {
                Some(e) => e,
                None => (at.biso / (8.0 * PI * PI), u.of(Field::AtomB(p, a)) / (8.0 * PI * PI)),
            };
            let ratio = occ[a] / at.occ.max(1e-300);
            let osu = u.of(Field::AtomOcc(p, a)) * ratio;
            let _ = writeln!(
                w,
                "{} {} {} {} {} {} {} {} {} {}",
                quoted(&labels[a]),
                element(&at.species),
                xs[0],
                xs[1],
                xs[2],
                with_su(uiso, usu, 5),
                if at.beta.is_some() { "Uani" } else { "Uiso" },
                with_su(occ[a], osu, 4),
                site.multiplicity,
                site.wyckoff.map_or("?".to_string(), |e| e.letter.to_string()) + if ph.rigid.is_empty() { "" } else if grouped(a) { " G" } else { " ." }
            );
        }
        if st.atoms.iter().any(|a| a.beta.is_some()) {
            let r = st.cell.reciprocal();
            let rs = [r.a, r.b, r.c];
            let idx = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];
            let _ = writeln!(w, "loop_\n_atom_site_aniso_label\n_atom_site_aniso_U_11\n_atom_site_aniso_U_22\n_atom_site_aniso_U_33\n_atom_site_aniso_U_12\n_atom_site_aniso_U_13\n_atom_site_aniso_U_23");
            // A TLS group's atoms' U's with su's through T, L, S and the group's values (L367).
            let tls = crate::geometry::tls_u(model, p, prop.as_ref());
            for (a, at) in st.atoms.iter().enumerate() {
                let Some(b) = at.beta else { continue };
                let us: Vec<String> = (0..6).map(|k| {
                    let f = 2.0 * PI * PI * rs[idx[k].0] * rs[idx[k].1];
                    match tls.iter().find(|t| t.0 == a) {
                        Some((_, t)) => with_su(t[k].0, t[k].1, 5),
                        None => with_su(b[k] / f, u.of(Field::AtomBeta(p, a, k)) / f, 5),
                    }
                }).collect();
                let _ = writeln!(w, "{} {}", quoted(&labels[a]), us.join(" "));
            }
        }
        // Geometry (P7): each bond once, from the cation (or the atom listed first), its angles; all flagged for
        // publication, being the bonds of the rule.
        let _ = writeln!(w, "loop_\n_geom_bond_atom_site_label_1\n_geom_bond_atom_site_label_2\n_geom_bond_site_symmetry_2\n_geom_bond_distance\n_geom_bond_publ_flag");
        for b in &geo.bonds {
            let _ = writeln!(w, "{} {} {} {} yes", quoted(&labels[b.i]), quoted(&labels[b.j]), crate::geometry::symmetry_code(b), with_su(b.distance.0, b.distance.1, 4));
        }
        if geo.bonds.is_empty() {
            let _ = writeln!(w, "? ? ? ? ?");
        }
        let _ = writeln!(w, "loop_\n_geom_angle_atom_site_label_1\n_geom_angle_atom_site_label_2\n_geom_angle_atom_site_label_3\n_geom_angle_site_symmetry_1\n_geom_angle_site_symmetry_3\n_geom_angle\n_geom_angle_publ_flag");
        for a in &geo.angles {
            let _ = writeln!(w, "{} {} {} {} {} {} yes", quoted(&labels[a.first.j]), quoted(&labels[a.centre]), quoted(&labels[a.second.j]), crate::geometry::symmetry_code(&a.first), crate::geometry::symmetry_code(&a.second), with_su(a.angle.0, a.angle.1, 2));
        }
        if geo.angles.is_empty() {
            let _ = writeln!(w, "? ? ? ? ? ? ?");
        }
    }
    Ok(out)
}

/// Pauling electronegativities (CRC Handbook of Chemistry and Physics, "Electronegativity", the Pauling scale; noble
/// gases and elements without a value left out), for the order of a formula's elements (`phase_name`).
const ELECTRONEGATIVITY: &[(&str, f64)] = &[
    ("H", 2.20), ("Li", 0.98), ("Be", 1.57), ("B", 2.04), ("C", 2.55), ("N", 3.04), ("O", 3.44), ("F", 3.98),
    ("Na", 0.93), ("Mg", 1.31), ("Al", 1.61), ("Si", 1.90), ("P", 2.19), ("S", 2.58), ("Cl", 3.16),
    ("K", 0.82), ("Ca", 1.00), ("Sc", 1.36), ("Ti", 1.54), ("V", 1.63), ("Cr", 1.66), ("Mn", 1.55), ("Fe", 1.83),
    ("Co", 1.88), ("Ni", 1.91), ("Cu", 1.90), ("Zn", 1.65), ("Ga", 1.81), ("Ge", 2.01), ("As", 2.18), ("Se", 2.55),
    ("Br", 2.96), ("Rb", 0.82), ("Sr", 0.95), ("Y", 1.22), ("Zr", 1.33), ("Nb", 1.6), ("Mo", 2.16), ("Tc", 1.9),
    ("Ru", 2.2), ("Rh", 2.28), ("Pd", 2.20), ("Ag", 1.93), ("Cd", 1.69), ("In", 1.78), ("Sn", 1.96), ("Sb", 2.05),
    ("Te", 2.1), ("I", 2.66), ("Cs", 0.79), ("Ba", 0.89), ("La", 1.10), ("Ce", 1.12), ("Pr", 1.13), ("Nd", 1.14),
    ("Sm", 1.17), ("Gd", 1.20), ("Dy", 1.22), ("Ho", 1.23), ("Er", 1.24), ("Tm", 1.25), ("Lu", 1.27), ("Hf", 1.3),
    ("Ta", 1.5), ("W", 2.36), ("Re", 1.9), ("Os", 2.2), ("Ir", 2.20), ("Pt", 2.28), ("Au", 2.54), ("Hg", 2.00),
    ("Tl", 1.62), ("Pb", 1.87), ("Bi", 2.02), ("Po", 2.0), ("At", 2.2), ("Fr", 0.7), ("Ra", 0.9), ("Ac", 1.1),
    ("Th", 1.3), ("Pa", 1.5), ("U", 1.38), ("Np", 1.36), ("Pu", 1.28),
];

/// The name the core proposes for phase p (P10, L114), so that the CIF's data block, the tables, the file names and the
/// app use one: the phase's own name when it is one word of letters, digits, '-' or '_' (not "phase1"), else the
/// formula of its content, Z taken out, its elements by increasing electronegativity, the electropositive first (IUPAC,
/// Nomenclature of Inorganic Chemistry 2005, IR-4.4.2: "Y2O3", "Sr2ScSbO6"; an element without a value goes last, in
/// the order its atoms are written); "phase{n}" for a phase without a structure and no usable name.
pub fn phase_name(model: &Model, p: usize) -> String {
    let ph = &model.phases[p];
    let own = ph.name.trim().trim_start_matches("Name:").trim();
    let word = !own.is_empty() && own.chars().all(|c| c.is_ascii_alphanumeric() || c == '-' || c == '_') && !own.to_ascii_lowercase().starts_with("phase");
    if word {
        return own.to_string();
    }
    let Intensities::Structure(st, _) = &ph.intensities else { return format!("phase{}", p + 1) };
    let Ok(occ) = bonds::chemical_occupancies(st) else { return format!("phase{}", p + 1) };
    let mut counts: Vec<(String, f64)> = Vec::new();
    for (a, o) in st.atoms.iter().zip(&occ) {
        let Ok(site) = Site::new(&st.group, a.x) else { return format!("phase{}", p + 1) };
        let el = element(&a.species);
        match counts.iter_mut().find(|c| c.0 == el) {
            Some(c) => c.1 += o * site.multiplicity as f64,
            None => counts.push((el, o * site.multiplicity as f64)),
        }
    }
    let en = |e: &str| ELECTRONEGATIVITY.iter().find(|x| x.0 == e).map_or(f64::INFINITY, |x| x.1);
    counts.sort_by(|a, b| en(&a.0).total_cmp(&en(&b.0)));
    let least = counts.iter().map(|c| c.1).fold(f64::INFINITY, f64::min).round().max(1.0) as usize;
    let z = (1..=least).rev().find(|z| counts.iter().all(|c| ((c.1 / *z as f64) - (c.1 / *z as f64).round()).abs() < 1e-3)).unwrap_or(1);
    let name: String = counts
        .iter()
        .map(|(e, n)| {
            let q = with_su(n / z as f64, 0.0, 3);
            if q == "1" { e.clone() } else { format!("{e}{q}") }
        })
        .collect();
    if name.is_empty() { format!("phase{}", p + 1) } else { name }
}

/// R_F² of phase p (P8): Σ|F²_obs − F²_calc| / ΣF²_obs over its reflections, F²_obs = I_obs F²_calc / I_calc with I_obs
/// apportioned by the model (`diagnostics::reflections`); None without a structure.
fn r_fsqd(model: &Model, data: &Pattern, p: usize) -> Option<f64> {
    let (mut num, mut den) = (0.0, 0.0);
    for r in crate::diagnostics::reflections(model, data).iter().filter(|r| r.phase == p && r.icalc > 0.0) {
        let f2 = r.f2?;
        let obs = r.iobs * f2 / r.icalc;
        num += (obs - f2).abs();
        den += obs.abs();
    }
    (den > 0.0).then(|| num / den)
}

/// The constraints of a parameter set (P8): every value that moves with a free parameter but is not its own (tied by
/// the symmetry, a link, a shared codeword or a variable), counted once.
pub fn constraints(set: &ParameterSet) -> usize {
    set.fields().len().saturating_sub(set.free.iter().filter(|q| q.value.is_none()).count())
}

/// How the su's were obtained (P2; the IUCr's Statistical descriptors, Rec. 8).
fn special_details(chi2: f64) -> String {
    format!(
        "Standard uncertainties from the inverse of the least-squares normal matrix at the final values, multiplied by\n\
         chi^2^ = S/(N - P) = {:.4} over all the points used (Prince, 2004, eq. 6.4.3); those of derived quantities\n\
         (positions of a mode phase, distances, angles, the volume) propagated through the full covariance. Values\n\
         fixed by symmetry or by a constraint carry none. su's rounded up, the IUCr's 2-19 rule (Schwarzenbach et al.,\n\
         1995).",
        chi2
    )
}

/// The profile function of a phase, in words (P8).
fn profile_function(model: &Model, ph: &crate::powder::Phase) -> String {
    use crate::powder::Shape;
    let shape = match ph.shape {
        Shape::Gaussian => "Gaussian, H^2^ = U tan^2^\\q + V tan\\q + W + I~G~/cos^2^\\q".to_string(),
        Shape::Lorentzian => "Lorentzian".to_string(),
        Shape::PseudoVoigt => "pseudo-Voigt, H^2^ = U tan^2^\\q + V tan\\q + W + I~G~/cos^2^\\q, \\h = \\h~0~ + X 2\\q".to_string(),
        Shape::ThompsonCoxHastings => "Thompson-Cox-Hastings pseudo-Voigt (Thompson, Cox & Hastings, 1987), H~G~^2^ = U tan^2^\\q + V tan\\q + W + I~G~/cos^2^\\q, H~L~ = X tan\\q + Y/cos\\q".to_string(),
    };
    let irf = if matches!(model.instrumental, crate::powder::InstrumentalWidths::None) { "" } else { "; the instrument's widths from a resolution function added" };
    let asym = if ph.s_l != 0.0 || ph.d_l != 0.0 {
        "; axial-divergence asymmetry (Finger, Cox & Jephcoat, 1994)"
    } else if ph.asym.iter().any(|a| *a != 0.0) {
        "; asymmetry (Berar & Baldinozzi, 1993)"
    } else {
        ""
    };
    format!("{shape}{irf}{asym}")
}

/// The background function, in words (P8).
/// The background function in words, with its further terms (Publish's P57.1): a file's b1–b4, Debye-like terms,
/// inverse powers of 2θ, profiles from files.
fn background_function(model: &Model) -> String {
    use crate::powder::BackgroundModel as B;
    let bt = &model.background_terms;
    let mut s = match &model.background {
        B::Points(p) => match &bt.file {
            Some(f) => format!(
                "linear interpolation between {} points read from {}, as b~1~ B[(1 + b~3~)2\\q + b~4~] + b~2~ with b~1~ = {}, b~2~ = {}, b~3~ = {}, b~4~ = {}",
                p.len(),
                f.file.rsplit(['/', '\\']).next().unwrap_or(&f.file),
                num_text(f.scale),
                num_text(f.offset),
                num_text(f.stretch),
                num_text(f.shift)
            ),
            None => format!("linear interpolation between {} points", p.len()),
        },
        B::Spline(p) => format!("natural cubic spline through {} points", p.len()),
        B::Polynomial { origin, coefficients } => format!("polynomial of degree {} in (2\\q/{} - 1)", coefficients.len().saturating_sub(1), num_text(*origin)),
        B::Chebyshev { from, to, coefficients } => format!("Chebyshev polynomial of the first kind, {} terms, on {} to {} degrees 2\\q", coefficients.len(), num_text(*from), num_text(*to)),
        B::Fourier { coefficients } => format!("Fourier cosine series, {} terms", coefficients.len()),
    };
    if !bt.inverse.is_empty() {
        s += &format!(", plus {} inverse powers of 2\\q", bt.inverse.len());
    }
    if !bt.debye.is_empty() {
        let r: Vec<String> = bt.debye.iter().map(|d| num_text(d.distance)).collect();
        s += &format!(", plus {} Debye-like terms A sin(Qr)/(Qr) at r = {} \\%A", bt.debye.len(), r.join(", "));
    }
    if !bt.profiles.is_empty() {
        let f: Vec<String> = bt.profiles.iter().map(|p| format!("{} x {}", p.file.rsplit(['/', '\\']).next().unwrap_or(&p.file), num_text(p.coefficient))).collect();
        s += &format!(", plus {} profiles from files ({})", bt.profiles.len(), f.join(", "));
    }
    s
}

fn num_text(v: f64) -> String {
    with_su(v, 0.0, 4)
}

/// The preferred-orientation correction of a phase, in words (P8).
fn preferred_orientation(model: &Model, ph: &crate::powder::Phase) -> String {
    // Several axes (P30): each with its direction, March's r and fraction.
    if let Some(t) = &ph.texture {
        let axes: Vec<String> = t.axes.iter().map(|a| format!("[{}] r = {}, fraction {}", a.direction.map(|v| with_su(v, 0.0, 3)).join(" "), with_su(a.r, 0.0, 4), with_su(a.fraction, 0.0, 4))).collect();
        let geometry = if t.steps > 0 { format!(", in an anvil cell with the beam along the force (Wright et al., 1996; {} steps)", t.steps) } else { String::new() };
        return format!("March-Dollase (Dollase, 1986) along {} axes{geometry}: {}", t.axes.len(), axes.join("; "));
    }
    let d = ph.pref_direction.map(|v| with_su(v, 0.0, 3)).join(" ");
    let none = if model.march { ph.pref[0] == 1.0 && ph.pref[1] == 0.0 } else { ph.pref[0] == 0.0 };
    if none {
        "none".into()
    } else if model.march {
        format!("March-Dollase (Dollase, 1986) along [{d}], G~1~ = {}, untextured fraction G~2~ = {}", with_su(ph.pref[0], 0.0, 4), with_su(ph.pref[1], 0.0, 4))
    } else {
        format!("exponential, exp(-G~1~ \\a^2^), along [{d}], G~1~ = {}, G~2~ = {}", with_su(ph.pref[0], 0.0, 4), with_su(ph.pref[1], 0.0, 4))
    }
}

fn crystal_system(point_group: &str) -> &'static str {
    match point_group {
        "1" | "-1" => "triclinic",
        "2" | "m" | "2/m" => "monoclinic",
        "222" | "mm2" | "mmm" => "orthorhombic",
        "4" | "-4" | "4/m" | "422" | "4mm" | "-42m" | "-4m2" | "4/mmm" => "tetragonal",
        "3" | "-3" | "32" | "3m" | "-3m" | "321" | "312" | "3m1" | "31m" | "-3m1" | "-31m" => "trigonal",
        "6" | "-6" | "6/m" | "622" | "6mm" | "-6m2" | "-62m" | "6/mmm" => "hexagonal",
        _ => "cubic",
    }
}

#[cfg(test)]
mod tests {
    use super::with_su;

    #[test]
    fn uncertainties_as_the_iucr_writes_them() {
        // Rounded up (P6): the su is never written smaller than it is.
        assert_eq!(with_su(5.748193, 0.0000512, 5), "5.74819(6)");
        assert_eq!(with_su(0.01183, 0.000197, 5), "0.0118(2)");
        assert_eq!(with_su(0.30685, 0.00019, 5), "0.30685(19)");
        assert_eq!(with_su(244.9112, 0.0031, 3), "244.911(4)");
        assert_eq!(with_su(90.0, 0.0, 4), "90");
        assert_eq!(with_su(0.25, 0.0, 5), "0.25");
        assert_eq!(with_su(0.0045, 0.00019, 5), "0.00450(19)");
        assert_eq!(with_su(1234.5, 12.0, 2), "1235(12)");
        assert_eq!(with_su(1234.5, 25.0, 2), "1230(30)");
        // Publish's examples (P6): 2.2865 ± 0.0032 and case A's b.
        assert_eq!(with_su(2.2865, 0.0032, 4), "2.287(4)");
        assert_eq!(with_su(5.680030, 0.000103, 5), "5.68003(11)");
        // Exact su's stay as they are; 1.96 units of a place become 20 of the next, written 2; 9.2 becomes 10.
        assert_eq!(with_su(1.0, 0.004, 4), "1.000(4)");
        assert_eq!(with_su(1.0, 0.0019, 4), "1.0000(19)");
        assert_eq!(with_su(1.0, 0.00196, 4), "1.000(2)");
        assert_eq!(with_su(1.0, 0.0092, 4), "1.000(10)");
        assert_eq!(with_su(1.0, 0.0015, 4), "1.0000(15)");
        assert_eq!(with_su(-0.00004, 0.0021, 4), "0.000(3)");
        assert_eq!(with_su(180.0, 0.0017, 4), "180.0000(17)");
    }
}
