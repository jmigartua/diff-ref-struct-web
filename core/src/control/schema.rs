//! What the control file holds and how it is written (L066). The keys are the parameter dictionary's ids
//! (`sample_displacement`, `caglioti_u`, `atom_biso`), with the dictionary's unit in a comment after each value; the
//! choices are readable words (`geometry = "debye-scherrer"`, `profile = "thompson-cox-hastings"`); FullProf's numeric
//! codes appear nowhere.
//!
//! ```toml
//! # diff-ref-struct control file
//! format = "diff-ref-struct 1"
//! title = "Y2O3"
//! radiation = "x-ray"                  # x-ray or neutron
//!
//! [pattern]                            # the instrument, the data, the excluded regions
//! data = "y2o3.dat"                    # relative to this file
//! reader = "free_steps"                # the reader registry's id (drs data formats)
//! wavelength_1 = 1.5406                # Å
//! ...
//! [pattern.background]                 # function = points, polynomial, chebyshev or fourier
//! [pattern.resolution]                 # the instrument's widths, when given
//! [pattern.fullprof]                   # an imported .pcr's range and excluded regions as written (L070)
//! [[pattern.operation]]                # the data operations, in order (the file is never changed)
//! [[pattern.scatterer]]                # a scattering factor of the user's for a species (FullProf's line 12, P11)
//! [[phase]]                            # one per phase: cell, profile, scale, preferred orientation...
//! [[phase.atom]]                       # one per atom of the phase above
//! [[phase.mode]]                       # symmetry modes (R29): name, irrep, amplitude, vectors on the atoms above
//! [[phase.rigid_group]]                # a rigid group (P41): its kind, origin, orientation; its atoms name it
//! [[restraint]]                        # distances and angles
//! [refine]                             # what is refined: parameters, limits, linear restraints, cycles, steering
//! [[refine.link]]                      # values that follow others: value = coefficient × follows + offset (W12)
//! [[refine.variable]]                  # new variables refined in place of the values they move (W12)
//! [[refine.step]]                      # a refinement plan: what each step frees, its damping, cycles, stop rule
//! [convention]                         # the rules where FullProf departs from its papers (true: FullProf's)
//! ```
//!
//! Reading checks every key and value and reports what is wrong where it is: an unknown key (with the nearest known
//! one), a value of the wrong kind, a missing value, an atom or address that does not exist, a space group that does
//! not read.

use super::syntax::{self, float, string, Diagnostic, Document, Entry, Item, Span, Table, Value};
use super::Control;
use crate::cell::Cell;
use crate::data_ops::Operation;
use crate::jacobian::Field;
use crate::powder::{BackgroundModel, Convention, InstrumentalWidths, Intensities, Model, Phase, Shape, Weighting};
use crate::refinable::{address, dictionary_id, field_of, LinearRestraint, Refinement};
use crate::steering::{Behaviour, Class, Guards, Link, PlanStep, Steering, Variable};
use crate::restraints::{Kind, Operator, Restraint};
use crate::scattering::Radiation;
use crate::spacegroup::SpaceGroup;
use crate::structure_factor::{Atom, Structure};
use crate::symop::SymOp;

/// The format line every file starts with.
pub const FORMAT: &str = "diff-ref-struct 1";

const ROOT_KEYS: &[&str] = &["format", "title", "radiation"];
const PATTERN_KEYS: &[&str] = &[
    "data",
    "reader",
    "reader_options",
    "wavelength_1",
    "wavelength_2",
    "wavelength_ratio",
    "zero_shift",
    "sample_displacement",
    "sample_transparency",
    "geometry",
    "monochromator_polarisation",
    "polarisation_fraction",
    "absorption_mur",
    "microabsorption_p0",
    "microabsorption_cp",
    "microabsorption_tau",
    "peak_window",
    "asymmetry_limit",
    "preferred_orientation_function",
    "weighting",
    "background_ceiling",
    "fullprof_step",
    "excluded_regions",
    "tof_dtt1",
    "tof_dtt2",
    "tof_dtt_1overd",
    "tof_thermal_zero",
    "tof_thermal_dtt1",
    "tof_thermal_dtt2",
    "tof_crossover",
    "tof_crossover_width",
    "tof_bank_two_theta",
    "tof_absorption_correction",
    "tof_max",
];
const SCATTERER_KEYS: &[&str] = &["species", "scattering_length", "scattering_length_imaginary", "f_prime", "f_double_prime", "form_factor", "form_factor_table"];
const BACKGROUND_KEYS: &[&str] = &[
    "function",
    "interpolation",
    "file",
    "background_points",
    "background_origin",
    "background_polynomial",
    "from",
    "to",
    "background_chebyshev",
    "background_fourier",
    "background_debye_amplitude",
    "background_debye_distance",
    "background_inverse_power",
    "profile_files",
    "background_profile",
    "background_file_scale",
    "background_file_offset",
    "background_file_stretch",
    "background_file_shift",
];
const RESOLUTION_KEYS: &[&str] = &["function", "instrument_resolution"];
const FULLPROF_KEYS: &[&str] = &["range", "excluded_regions"];
const PHASE_KEYS: &[&str] = &[
    "name",
    "space_group",
    "hall",
    "operators",
    "cell_a",
    "cell_b",
    "cell_c",
    "cell_alpha",
    "cell_beta",
    "cell_gamma",
    "profile",
    "scale_factor",
    "overall_b",
    "pv_eta0",
    "caglioti_u",
    "caglioti_v",
    "caglioti_w",
    "lorentzian_strain",
    "pv_eta_slope",
    "lorentzian_size",
    "gaussian_size",
    "quartic_strain",
    "quartic_strain_lorentzian",
    "shape_size",
    "shape_size_kind",
    "shape_size_direction",
    "berar_baldinozzi",
    "axial_divergence_s_l",
    "axial_divergence_d_l",
    "po_g1",
    "po_g2",
    "po_direction",
    "po_axis_direction",
    "po_axis_ratio",
    "po_axis_fraction",
    "po_axis_steps",
    "qpa_atz",
    "brindley_coefficient",
    "internal_standard_weight_fraction",
    "valence_sum_range",
    "le_bail",
    "fullprof_reflection_limit",
    "reflection_intensity",
    "mode_reference_cell",
    "mode_reference",
    "mode_max_amplitude",
    "mode_fst",
    "mode_parent_group",
    "mode_parent_basis",
    "mode_parent_origin",
    "lambda2_uvw",
    "tof_sigma0",
    "tof_sigma1",
    "tof_sigma2",
    "tof_sigma_q",
    "tof_gamma0",
    "tof_gamma1",
    "tof_gamma2",
    "tof_alpha0",
    "tof_alpha1",
    "tof_alpha_q",
    "tof_beta0",
    "tof_beta1",
    "tof_beta_q",
    "tof_alpha0_thermal",
    "tof_beta0_thermal",
    "tof_alpha1_thermal",
    "tof_beta1_thermal",
    "tof_extinction",
    "tof_absorption",
];
const MODE_KEYS: &[&str] = &["name", "irrep", "mode_amplitude", "mode_vector"];
const ATOM_KEYS: &[&str] = &["label", "species", "charge", "atom_position", "atom_biso", "atom_occupancy", "atom_betas", "rigid_group", "rigid_length", "rigid_angle", "rigid_connect"];
const RIGID_KEYS: &[&str] = &["kind", "rigid_origin", "rigid_orientation", "satellite", "fullprof_option", "rigid_tls_t", "rigid_tls_l", "rigid_tls_s", "tls_centre"];
const RESTRAINT_KEYS: &[&str] = &["phase", "distance", "angle", "operator", "operators", "target", "sigma"];
const REFINE_KEYS: &[&str] = &["cycles", "convergence", "lambda", "marquardt", "damping", "divergence", "parameters", "limits", "steering"];
const LINK_KEYS: &[&str] = &["value", "follows", "coefficient", "offset"];
const VARIABLE_KEYS: &[&str] = &["name", "value", "terms"];
const STEP_KEYS: &[&str] = &["name", "free", "damping", "cycles", "convergence"];
const LINEAR_KEYS: &[&str] = &["name", "terms", "target", "sigma"];
const CONVENTION_KEYS: &[&str] = &[
    "lorentz_polarisation",
    "berar_baldinozzi",
    "axial_divergence",
    "preferred_orientation",
    "excluded_points",
    "peak_window",
    "reflections",
    "microstructure",
    "angle_restraints",
    "window_edge",
    "background_edge",
    "displacement",
    "background_bounds",
];
const GEOMETRY: &[(i32, &str)] = &[(-1, "corrected"), (0, "debye-scherrer"), (1, "flat-plate-psd"), (2, "transmission"), (3, "synchrotron")];
const SHAPES: &[(Shape, &str)] = &[(Shape::Gaussian, "gaussian"), (Shape::Lorentzian, "lorentzian"), (Shape::PseudoVoigt, "pseudo-voigt"), (Shape::ThompsonCoxHastings, "thompson-cox-hastings")];
const TOF_ABSORPTION: &[(crate::tof::Absorption, &str)] = &[(crate::tof::Absorption::FlatPlate, "flat-plate"), (crate::tof::Absorption::Cylinder, "cylinder"), (crate::tof::Absorption::Exponential, "exponential")];
const WEIGHTS: &[(Weighting, &str)] = &[(Weighting::Variance, "variance"), (Weighting::MaximumLikelihood, "maximum-likelihood"), (Weighting::Unit, "unit")];

/// The tables a control file may hold and the keys each takes, as the reader checks them (P49's completion):
/// (path, an array of tables, its keys). `[[pattern.operation]]` takes `op` and the settings of the operation it
/// names (`operation_table`).
pub const TABLES: &[(&[&str], bool, &[&str])] = &[
    (&[], false, ROOT_KEYS),
    (&["pattern"], false, PATTERN_KEYS),
    (&["pattern", "background"], false, BACKGROUND_KEYS),
    (&["pattern", "resolution"], false, RESOLUTION_KEYS),
    (&["pattern", "fullprof"], false, FULLPROF_KEYS),
    (&["pattern", "operation"], true, &["op"]),
    (&["pattern", "scatterer"], true, SCATTERER_KEYS),
    (&["phase"], true, PHASE_KEYS),
    (&["phase", "atom"], true, ATOM_KEYS),
    (&["phase", "mode"], true, MODE_KEYS),
    (&["phase", "rigid_group"], true, RIGID_KEYS),
    (&["restraint"], true, RESTRAINT_KEYS),
    (&["refine"], false, REFINE_KEYS),
    (&["refine", "linear"], true, LINEAR_KEYS),
    (&["refine", "link"], true, LINK_KEYS),
    (&["refine", "variable"], true, VARIABLE_KEYS),
    (&["refine", "step"], true, STEP_KEYS),
    (&["convention"], false, CONVENTION_KEYS),
];

/// A `[[phase]]` table with its atoms', modes' and rigid groups' tables.
type PhaseTables<'a> = (&'a Table, Vec<&'a Table>, Vec<&'a Table>, Vec<&'a Table>);

/// A phase's rigid groups (P41): each `[[phase.rigid_group]]` table's kind, origin, orientation and satellite anchor,
/// its atoms the phase's `[[phase.atom]]` tables that name it (`rigid_group`), in their order, with their internal
/// coordinates (`rigid_length`, `rigid_angle`) and a Z-matrix atom's `rigid_connect`.
fn read_rigid(r: &mut R, tables: &[&Table], atoms: &[&Table]) -> Vec<crate::rigid_body::Group> {
    use crate::rigid_body::{Coordinates, Group, Member, Satellite};
    let mut out: Vec<Group> = Vec::new();
    for t in tables {
        r.check_keys(t, RIGID_KEYS);
        let coordinates = r.choice(t, "kind", &Coordinates::NAMES, Coordinates::Spherical);
        let o = r.nums_key(t, "rigid_origin", Some(3), None);
        let a = r.nums_key(t, "rigid_orientation", Some(3), None);
        let satellite = t.get("satellite").and_then(|e| r.nums(&e.value, "\"satellite\"", Some(3))).map(|v| Satellite { parent: (v[0] as usize).saturating_sub(1), atom: v[1] as usize, from: v[2] as usize });
        if let (Some(sat), Some(e)) = (satellite, t.get("satellite")) {
            if sat.parent >= out.len() {
                r.err(e.value.span, format!("satellite: group {} is not a group before this one", sat.parent + 1));
            }
        }
        let option = t.get("fullprof_option").map_or(if coordinates == Coordinates::ZMatrix { 4.0 } else if satellite.is_some() { 2.0 + (satellite.map_or(0, |s| s.parent + 1) as f64) / 100.0 } else { 1.0 }, |_| r.f(t, "fullprof_option", None));
        // TLS displacements (P41, L367): T (Å²), L (rad²), S (Å·rad), the centre less the origin (fractions).
        let displacements = t.get("rigid_tls_t").map(|_| {
            let tt = r.nums_key(t, "rigid_tls_t", Some(6), None);
            let l = r.nums_key(t, "rigid_tls_l", Some(6), Some(vec![0.0; 6]));
            let s = r.nums_key(t, "rigid_tls_s", Some(9), Some(vec![0.0; 9]));
            let c = r.nums_key(t, "tls_centre", Some(3), Some(vec![0.0; 3]));
            crate::rigid_body::Tls { t: std::array::from_fn(|i| tt[i]), l: std::array::from_fn(|i| l[i]), s: std::array::from_fn(|i| s[i]), centre: [c[0], c[1], c[2]] }
        });
        if displacements.is_none() {
            for k in ["rigid_tls_l", "rigid_tls_s", "tls_centre"] {
                if let Some(e) = t.get(k) {
                    r.err(e.key_span, format!("\"{k}\" belongs to a TLS group, one with rigid_tls_t"));
                }
            }
        }
        out.push(Group { option, coordinates, degrees: true, origin: [o[0], o[1], o[2]], theta: a[0], phi: a[1], chi: a[2], satellite, members: Vec::new(), displacements });
    }
    for (i, at) in atoms.iter().enumerate() {
        let Some(e) = at.get("rigid_group") else {
            for k in ["rigid_length", "rigid_angle", "rigid_connect"] {
                if let Some(e) = at.get(k) {
                    r.err(e.key_span, format!("\"{k}\" belongs to an atom of a rigid group (rigid_group = n)"));
                }
            }
            continue;
        };
        let g = match &e.value.value {
            Value::Integer(n) if *n >= 1 && (*n as usize) <= out.len() => *n as usize - 1,
            _ => {
                r.err(e.value.span, format!("\"rigid_group\" is the number of one of the phase's {} [[phase.rigid_group]] tables", out.len()));
                continue;
            }
        };
        let (len, ang) = out[g].coordinates.split();
        let l = if len.is_empty() { Vec::new() } else { r.nums_key(at, "rigid_length", Some(len.len()), None) };
        let an = if ang.is_empty() { Vec::new() } else { r.nums_key(at, "rigid_angle", Some(ang.len()), None) };
        let mut internal = [0.0; 3];
        for (k, &c) in len.iter().enumerate() {
            internal[c] = l.get(k).copied().unwrap_or(0.0);
        }
        for (k, &c) in ang.iter().enumerate() {
            internal[c] = an.get(k).copied().unwrap_or(0.0);
        }
        let n = out[g].members.len();
        let connect = if out[g].coordinates == Coordinates::ZMatrix && n > 0 {
            let c = r.nums_key(at, "rigid_connect", Some(3), None);
            let c: [usize; 3] = std::array::from_fn(|k| c.get(k).copied().unwrap_or(0.0) as usize);
            if c.iter().any(|&k| k > n) || c[0] == 0 || (n >= 2 && c[1] == 0) {
                r.err(at.get("rigid_connect").map_or(at.header.unwrap_or_default(), |e| e.value.span), format!("rigid_connect: N1 N2 N3 are atoms of the group before this one (1 … {n}), N1 and, from the third atom, N2 given"));
            }
            c
        } else {
            [0; 3]
        };
        out[g].members.push(Member { label: r.s(at, "label").map_or(String::new(), |s| s.0), index: i, internal, connect });
    }
    for (k, g) in out.iter().enumerate() {
        if g.members.is_empty() {
            r.err(tables[k].header.unwrap_or_default(), format!("[[phase.rigid_group]] {}: no atom names it (rigid_group = {})", k + 1, k + 1));
        }
    }
    out.retain(|g| !g.members.is_empty());
    out
}

/// The values a key takes from a closed list, as the reader accepts them: (table path, key, values).
pub fn choices() -> Vec<(&'static [&'static str], &'static str, Vec<&'static str>)> {
    let mut v: Vec<(&'static [&'static str], &'static str, Vec<&'static str>)> = vec![
        (&[], "radiation", vec!["x-ray", "neutron"]),
        (&["pattern"], "geometry", GEOMETRY.iter().map(|g| g.1).collect()),
        (&["pattern"], "weighting", WEIGHTS.iter().map(|w| w.1).collect()),
        (&["phase", "rigid_group"], "kind", crate::rigid_body::Coordinates::NAMES.iter().map(|n| n.1).collect()),
        (&["phase"], "profile", SHAPES.iter().map(|s| s.1).collect()),
        (&["pattern"], "preferred_orientation_function", vec!["march-dollase", "exponential"]),
        (&["pattern"], "tof_absorption_correction", TOF_ABSORPTION.iter().map(|a| a.1).collect()),
        (&["phase"], "le_bail", vec!["true", "false"]),
        (&["phase"], "shape_size_kind", vec!["platelet", "needle", "antiphase"]),
        (&["pattern", "background"], "function", vec!["points", "polynomial", "chebyshev", "fourier"]),
        (&["pattern", "background"], "interpolation", vec!["linear", "spline"]),
        (&["pattern", "operation"], "op", crate::operation_table::OPERATIONS.iter().map(|o| o.id).collect()),
    ];
    for k in CONVENTION_KEYS {
        v.push((&["convention"], k, vec!["true", "false"]));
    }
    v
}

fn geometry_name(ilo: i32) -> String {
    GEOMETRY.iter().find(|g| g.0 == ilo).map_or(format!("fullprof-ilo-{ilo}"), |g| g.1.to_string())
}

/// The dictionary's unit of an id, as a comment.
fn unit(id: &str) -> String {
    match crate::parameters::lookup(id).map(|p| p.unit) {
        Some("1") | None => String::new(),
        Some(u) => format!("# {u}"),
    }
}

// The writer ---------------------------------------------------------------------------------------------------------

struct W {
    out: String,
}

impl W {
    fn line(&mut self, s: &str) {
        self.out.push_str(s);
        self.out.push('\n');
    }

    /// `key = value`, the comment aligned at column 41 when there is one.
    fn kv(&mut self, key: &str, value: &str, comment: &str) {
        let mut l = format!("{key} = {value}");
        if !comment.is_empty() {
            let pad = 40usize.saturating_sub(l.chars().count()).max(1);
            l.push_str(&" ".repeat(pad));
            l.push_str(comment);
        }
        self.line(&l);
    }

    fn num(&mut self, id: &str, v: f64) {
        self.kv(id, &float(v), &unit(id));
    }
}

fn floats(v: &[f64]) -> String {
    format!("[{}]", v.iter().map(|&x| float(x)).collect::<Vec<_>>().join(", "))
}

/// A multi-line array of rows, one row a line.
fn rows(rows: &[Vec<String>]) -> String {
    if rows.is_empty() {
        return "[]".into();
    }
    let mut s = String::from("[\n");
    for r in rows {
        s.push_str(&format!("  [{}],\n", r.join(", ")));
    }
    s.push(']');
    s
}

/// An atom's name in a restraint or an address: its label, or "label:n" when an earlier atom has the same label.
fn atom_name(st: &Structure, a: usize) -> String {
    let l = &st.atoms[a].label;
    if st.atoms[..a].iter().any(|b| &b.label == l) {
        format!("{l}:{}", a + 1)
    } else {
        l.clone()
    }
}

fn operator_text(op: &Operator) -> String {
    let rotation = SymOp { r: op.rot, t: [0; 3] }.to_xyz();
    format!("{{ rotation = {}, translation = {} }}", string(&rotation), floats(&op.t))
}

/// The canonical text of a control file: every value written, in the order of `schema`'s tables.
pub fn write(c: &Control) -> String {
    let m = &c.model;
    let mut w = W { out: String::new() };
    w.line("# diff-ref-struct control file (https://jmigartua.github.io/diff-ref-struct-site/)");
    w.kv("format", &string(FORMAT), "");
    w.kv("title", &string(&c.title), "");
    w.kv("radiation", &string(if m.neutrons { "neutron" } else { "x-ray" }), "# x-ray or neutron");
    w.line("");
    w.line("[pattern]");
    if let Some(d) = &c.data {
        w.kv("data", &string(&d.file), "# relative to this file");
        if let Some(r) = &d.reader {
            w.kv("reader", &string(r), "# drs data formats");
        }
        if !d.reader_options.is_empty() {
            use crate::formats::options::OptionValue;
            let o: Vec<String> = d
                .reader_options
                .iter()
                .map(|(k, v)| match v {
                    OptionValue::Integer(n) => format!("{k} = {n}"),
                    OptionValue::Number(x) => format!("{k} = {}", float(*x)),
                    OptionValue::Text(t) => format!("{k} = {}", string(t)),
                })
                .collect();
            w.kv("reader_options", &format!("{{ {} }}", o.join(", ")), "# the reader's options (drs data formats)");
        }
    }
    w.num("wavelength_1", m.lambda[0]);
    w.num("wavelength_2", m.lambda[1]);
    w.num("wavelength_ratio", m.ratio);
    if m.tof.is_some() {
        w.kv("zero_shift", &float(m.zero), "# µs: the time-of-flight zero");
    } else {
        w.num("zero_shift", m.zero);
    }
    w.num("sample_displacement", m.sycos);
    w.num("sample_transparency", m.sysin);
    w.kv("geometry", &string(&geometry_name(m.ilo)), "# corrected, debye-scherrer, flat-plate-psd, transmission, synchrotron");
    w.num("monochromator_polarisation", m.cthm);
    w.num("polarisation_fraction", m.rpolarz);
    w.num("absorption_mur", m.mur);
    // Bragg–Brentano surface micro-absorption (P29): written only when it is there.
    if m.microabsorption.iter().any(|&v| v != 0.0) {
        for (k, v) in ["microabsorption_p0", "microabsorption_cp", "microabsorption_tau"].iter().zip(m.microabsorption) {
            w.num(k, v);
        }
    }
    w.num("peak_window", m.wdt);
    w.num("asymmetry_limit", m.asym_limit);
    w.kv("preferred_orientation_function", &string(if m.march { "march-dollase" } else { "exponential" }), "# march-dollase or exponential");
    w.kv("weighting", &string(WEIGHTS.iter().find(|x| x.0 == m.weighting).map_or("variance", |x| x.1)), "# variance, maximum-likelihood or unit");
    if let Some(c) = m.background_ceiling {
        w.kv("background_ceiling", &float(c), "# counts: FullProf's background is 0 above it (convention background_bounds)");
    }
    if let Some(s) = m.fullprof_step {
        w.kv("fullprof_step", &float(s), "# deg: FullProf's step, as its .out prints it, for its tests");
    }
    w.kv("excluded_regions", &rows(&m.excluded.iter().map(|&(a, b)| vec![float(a), float(b)]).collect::<Vec<_>>()), &if m.tof.is_some() { "# µs".to_string() } else { unit("excluded_regions") });
    // A time-of-flight bank (P05): T(d), its angle, the sample's absorption, the last time of flight.
    if let Some(b) = &m.tof {
        let n = if b.thermal.is_some() { 8 } else { 3 };
        for (k, id) in crate::tof::BANK_FIELDS.iter().enumerate().take(n) {
            w.num(id, b.get(k));
        }
        w.num("tof_bank_two_theta", b.two_theta);
        w.kv("tof_absorption_correction", &string(TOF_ABSORPTION.iter().find(|a| a.0 == b.absorption).unwrap().1), "# flat-plate, cylinder or exponential");
        w.num("tof_max", b.tmax);
    }
    w.line("");
    w.line("[pattern.background]");
    match &m.background {
        BackgroundModel::Points(p) | BackgroundModel::Spline(p) => {
            w.kv("function", &string("points"), "# points, polynomial, chebyshev or fourier");
            let how = if matches!(m.background, BackgroundModel::Spline(_)) { "spline" } else { "linear" };
            w.kv("interpolation", &string(how), "# linear or spline");
            match &c.background_file {
                Some(f) => w.kv("file", &string(f), "# the points, relative to this file (ours or FullProf's .bac)"),
                None => w.kv("background_points", &rows(&p.iter().map(|&(x, y)| vec![float(x), float(y)]).collect::<Vec<_>>()), "# [2θ (deg), counts]"),
            }
        }
        BackgroundModel::Polynomial { origin, coefficients } => {
            w.kv("function", &string("polynomial"), "# points, polynomial, chebyshev or fourier");
            w.num("background_origin", *origin);
            w.kv("background_polynomial", &floats(coefficients), &unit("background_polynomial"));
        }
        BackgroundModel::Chebyshev { from, to, coefficients } => {
            w.kv("function", &string("chebyshev"), "# points, polynomial, chebyshev or fourier");
            w.kv("from", &float(*from), "# deg 2θ");
            w.kv("to", &float(*to), "# deg 2θ");
            w.kv("background_chebyshev", &floats(coefficients), &unit("background_chebyshev"));
        }
        BackgroundModel::Fourier { coefficients } => {
            w.kv("function", &string("fourier"), "# points, polynomial, chebyshev or fourier");
            w.kv("background_fourier", &floats(coefficients), &unit("background_fourier"));
        }
    }
    // Further terms, added to the function (P23).
    let bt = &m.background_terms;
    if !bt.debye.is_empty() {
        w.kv("background_debye_amplitude", &floats(&bt.debye.iter().map(|d| d.amplitude).collect::<Vec<_>>()), &unit("background_debye_amplitude"));
        w.kv("background_debye_distance", &floats(&bt.debye.iter().map(|d| d.distance).collect::<Vec<_>>()), &unit("background_debye_distance"));
    }
    if !bt.inverse.is_empty() {
        w.kv("background_inverse_power", &floats(&bt.inverse), "# c1, c2, c3 of c1/2θ + c2/2θ² + c3/2θ³, 2θ in deg");
    }
    if !bt.profiles.is_empty() {
        w.kv("profile_files", &format!("[{}]", bt.profiles.iter().map(|p| string(&p.file)).collect::<Vec<_>>().join(", ")), "# profiles added to the background, relative to this file");
        w.kv("background_profile", &floats(&bt.profiles.iter().map(|p| p.coefficient).collect::<Vec<_>>()), "# the coefficient of each profile");
    }
    // A file's b1–b4 (FullProf's Nba = 1, eq. 3.6; the app's request 47): b1·B[(1 + b3)2θ + b4] + b2.
    if let Some(f) = &bt.file {
        w.num("background_file_scale", f.scale);
        w.num("background_file_offset", f.offset);
        w.num("background_file_stretch", f.stretch);
        w.num("background_file_shift", f.shift);
    }
    match &m.instrumental {
        InstrumentalWidths::None => {
            if let Some(i) = m.tof.as_ref().and_then(|b| b.instrument.as_ref()) {
                w.line("");
                w.line("[pattern.resolution]");
                w.kv("function", &string("tof"), "# a TOF resolution file's values, added to the bank's and the phases'");
                w.kv("instrument_resolution", &floats(&i.to_values()), "# Zero Dtt1 Dtt2, Zerot Dtt1t Dtt2t x-cross Width, σ₀ σ₁ σ₂, γ₀ γ₁ γ₂, α₀ α₁ β₀ β₁, α₀ₜ β₀ₜ α₁ₜ β₁ₜ");
            }
        }
        InstrumentalWidths::Function { res, uvwxyz } => {
            w.line("");
            w.line("[pattern.resolution]");
            w.kv("function", &res.to_string(), "# FullProf's resolution function 1–3 (U V W X Y Z)");
            w.kv("instrument_resolution", &floats(uvwxyz), &unit("instrument_resolution"));
        }
        InstrumentalWidths::Table(t) => {
            w.line("");
            w.line("[pattern.resolution]");
            w.kv("function", &string("table"), "# the widths at given angles");
            w.kv("instrument_resolution", &rows(&t.iter().map(|r| r.iter().map(|&v| float(v)).collect()).collect::<Vec<_>>()), "# [2θ, H_G, H_L] (deg)");
        }
    }
    if !c.fullprof.is_empty() {
        w.line("");
        w.line("[pattern.fullprof]                     # as the imported .pcr wrote them, for the export back to FullProf");
        if let Some(r) = c.fullprof.range {
            w.kv("range", &floats(&r), if m.tof.is_some() { "# TOF-min, Step, TOF-max (µs)" } else { "# Thmin, Step, Thmax (deg 2θ)" });
        }
        if !c.fullprof.excluded.is_empty() {
            w.kv("excluded_regions", &rows(&c.fullprof.excluded.iter().map(|&(a, b)| vec![float(a), float(b)]).collect::<Vec<_>>()), if m.tof.is_some() { "# µs, under FullProf's rule" } else { "# deg 2θ, under FullProf's rule" });
        }
    }
    // The user's scattering factors: each species whose scatterer is not the tables' (P11).
    let radiation = if m.neutrons { Radiation::Neutron } else { Radiation::Xray { lambda: m.lambda[0], anomalous: true } };
    let mut written: Vec<String> = Vec::new();
    for ph in &m.phases {
        if let Intensities::Structure(st, sc) = &ph.intensities {
            for (a, s) in st.atoms.iter().zip(sc) {
                if written.iter().any(|x| x.eq_ignore_ascii_case(&a.species)) {
                    continue;
                }
                if let Some(u) = crate::scattering::UserScatterer::of(&a.species, s, radiation) {
                    written.push(a.species.clone());
                    w.line("");
                    w.line("[[pattern.scatterer]]");
                    w.kv("species", &string(&u.species), "");
                    if let Some(b) = u.scattering_length {
                        w.kv("scattering_length", &float(b), "# 10⁻¹² cm");
                        if u.scattering_length_imaginary != 0.0 {
                            w.kv("scattering_length_imaginary", &float(u.scattering_length_imaginary), "# 10⁻¹² cm, an absorbing nucleus");
                        }
                    } else {
                        w.kv("f_prime", &float(u.f_prime), "# electrons");
                        w.kv("f_double_prime", &float(u.f_double_prime), "# electrons");
                    }
                    if let Some(k) = u.form_factor {
                        w.kv("form_factor", &floats(&k), "# a₁ b₁ a₂ b₂ a₃ b₃ a₄ b₄ c of f(s), s = sin θ/λ");
                    }
                    if !u.form_factor_table.is_empty() {
                        w.kv("form_factor_table", &rows(&u.form_factor_table.iter().map(|&(x, f)| vec![float(x), float(f)]).collect::<Vec<_>>()), "# [sin θ/λ (Å⁻¹), f], linearly between");
                    }
                }
            }
        }
    }
    if let Some(d) = &c.data {
        for op in &d.operations {
            w.line("");
            w.line("[[pattern.operation]]");
            for (k, v) in op.to_pairs() {
                w.kv(&k, &v, "");
            }
        }
    }
    for (p, ph) in m.phases.iter().enumerate() {
        w.line("");
        w.line("[[phase]]");
        w.kv("name", &string(&ph.name), "");
        if let Some(s) = &ph.group.setting {
            w.kv("space_group", &string(s.hm), "# Hermann–Mauguin");
        }
        w.kv("hall", &string(&ph.group.hall), "# Hall symbol: the operators");
        if SpaceGroup::from_hall(&ph.group.hall).map(|g| g.ops) != Ok(ph.group.ops.clone()) {
            w.kv("operators", &format!("[{}]", ph.group.ops.iter().map(|o| string(&o.to_xyz())).collect::<Vec<_>>().join(", ")), "");
        }
        let c = &ph.cell;
        for (id, v) in [("cell_a", c.a), ("cell_b", c.b), ("cell_c", c.c), ("cell_alpha", c.alpha), ("cell_beta", c.beta), ("cell_gamma", c.gamma)] {
            w.num(id, v);
        }
        w.kv("profile", &string(SHAPES.iter().find(|s| s.0 == ph.shape).map_or("pseudo-voigt", |s| s.1)), "# gaussian, lorentzian, pseudo-voigt, thompson-cox-hastings");
        w.num("scale_factor", ph.scale);
        w.num("overall_b", ph.bov);
        w.num("pv_eta0", ph.eta0);
        w.num("caglioti_u", ph.u);
        w.num("caglioti_v", ph.v);
        w.num("caglioti_w", ph.w);
        w.num(dictionary_id(Field::X(p), m), ph.x);
        w.num("lorentzian_size", ph.y);
        w.num("gaussian_size", ph.gauss_size);
        if let Some(t) = &ph.tof {
            for (k, id) in crate::tof::PROFILE_FIELDS.iter().enumerate() {
                if t.thermal.is_some() || !(13..=16).contains(&k) {
                    w.num(id, t.get(k));
                }
            }
        }
        if let Some(q) = &ph.quartic_strain {
            w.kv("quartic_strain", &rows(&q.terms.iter().map(|(t, v)| vec![t[0].to_string(), t[1].to_string(), t[2].to_string(), float(*v)]).collect::<Vec<_>>()), "# [H, K, L, S_HKL]: σ²(1/d²) = 10⁻⁸ Σ S_HKL h^H k^K l^L Å⁻⁴ (Stephens)");
            w.num("quartic_strain_lorentzian", q.lorentzian);
        }
        if let Some(sz) = &ph.shape_size {
            use crate::microstructure::SizeShape;
            w.num("shape_size", sz.value);
            let (kind, t) = match sz.shape {
                SizeShape::Platelet(t) => ("platelet", Some(t)),
                SizeShape::Needle(t) => ("needle", Some(t)),
                SizeShape::Antiphase => ("antiphase", None),
            };
            w.kv("shape_size_kind", &string(kind), "# platelet (|cos φ|), needle (sin φ) or antiphase (h, k, l of mixed parity)");
            if let Some(t) = t {
                w.kv("shape_size_direction", &floats(&t), "# [u, v, w]: the platelets' normal or the needles' axis");
            }
        }
        if let Some(u) = ph.uvw2 {
            w.kv("lambda2_uvw", &floats(&u), "# deg²: U, V, W of the second wavelength's peaks (FullProf's Ratio < 0)");
        }
        w.kv("berar_baldinozzi", &floats(&ph.asym), "# P1 P2 P3 P4");
        w.num("axial_divergence_s_l", ph.s_l);
        w.num("axial_divergence_d_l", ph.d_l);
        w.num("po_g1", ph.pref[0]);
        w.num("po_g2", ph.pref[1]);
        w.kv("po_direction", &floats(&ph.pref_direction), &unit("po_direction"));
        // Preferred orientation along several axes (P30): it takes the place of po_g1, po_g2 and po_direction.
        if let Some(t) = &ph.texture {
            w.kv("po_axis_direction", &rows(&t.axes.iter().map(|a| a.direction.iter().map(|v| float(*v)).collect()).collect::<Vec<_>>()), "# [h, k, l] of each axis: March–Dollase along several axes, in place of po_g1, po_g2, po_direction");
            w.kv("po_axis_ratio", &floats(&t.axes.iter().map(|a| a.r).collect::<Vec<_>>()), "# March's r of each axis");
            w.kv("po_axis_fraction", &floats(&t.axes.iter().map(|a| a.fraction).collect::<Vec<_>>()), "# the fraction of each axis: they sum to 1, the last is 1 minus the others");
            if t.steps > 0 {
                w.kv("po_axis_steps", &t.steps.to_string(), "# an anvil cell, the beam along the force: steps of the integration about each cone (FullProf's Nor = 3)");
            }
        }
        w.num("qpa_atz", ph.atz);
        w.num("brindley_coefficient", ph.brindley);
        if let Some(r) = ph.valence_sum_range {
            w.kv(crate::bonds::VALENCE_RANGE_KEY, &float(r), "# a valence sum is plausible within |q|(1 ± this) (default 0.10)");
        }
        if let Some(s) = ph.internal_standard {
            w.kv("internal_standard_weight_fraction", &float(s), "# weighed into the sample: the amorphous content follows (P33)");
        }
        w.kv("le_bail", if ph.extract { "true" } else { "false" }, "# intensities extracted (profile matching)");
        if let Some(l) = ph.reflection_limit {
            w.kv("fullprof_reflection_limit", &float(l), "# FullProf's S(1/d)max (1/Å), for its tests");
        }
        match &ph.intensities {
            Intensities::List(list) => {
                w.kv("reflection_intensity", &rows(&list.iter().map(|(h, mult, i)| vec![h[0].to_string(), h[1].to_string(), h[2].to_string(), mult.to_string(), float(*i)]).collect::<Vec<_>>()), "# [h, k, l, multiplicity, intensity]");
            }
            Intensities::Structure(st, _) => {
                // A phase described by symmetry modes (R29): its atoms at the reference structure's positions, the modes
                // after them; the positions the modes give are computed on reading.
                if let (None, Some(pa)) = (&ph.modes, &ph.parent) {
                    w.kv("mode_parent_group", &string(&pa.group), "# the parent's space group (P17: for the cell reduced to the parent's)");
                    w.kv("mode_parent_basis", &rows(&pa.basis.iter().map(|r| r.iter().map(|v| float(*v)).collect()).collect::<Vec<_>>()), "# P by rows: its columns are this phase's basis vectors in the parent's");
                    w.kv("mode_parent_origin", &floats(&pa.origin), "# p: this phase's origin in the parent's basis");
                }
                if let Some(ms) = &ph.modes {
                    if let Some(pa) = &ms.parent {
                        w.kv("mode_parent_group", &string(&pa.group), "# the parent's space group");
                        w.kv("mode_parent_basis", &rows(&pa.basis.iter().map(|r| r.iter().map(|v| float(*v)).collect()).collect::<Vec<_>>()), "# P by rows: its columns are this phase's basis vectors in the parent's");
                        w.kv("mode_parent_origin", &floats(&pa.origin), "# p: this phase's origin in the parent's basis");
                    }
                    if let Some(c) = &ms.reference_cell {
                        w.kv("mode_reference_cell", &floats(&[c.a, c.b, c.c, c.alpha, c.beta, c.gamma]), "# the reference (parent) cell in this basis: the modes' normalisation");
                    }
                    if let Some(b) = ms.max_amplitude {
                        w.kv("mode_max_amplitude", &float(b), "# Å: |amplitude| at most this (FullProf's Max_Amplitude)");
                    }
                    if let Some(sr) = &ms.sign_reference {
                        let a = sr.analysis.map_or(String::new(), |k| format!(", analysis = {k}"));
                        w.kv("mode_reference", &format!("{{ file = {}{a} }}", string(&sr.file)), "# the sign reference: an AMPLIMODES analysis or FullProf's .sum (E14b)");
                    }
                    if !ms.fst.is_empty() {
                        w.kv("mode_fst", &format!("[{}]", ms.fst.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(", ")), "# FullProf's FullProf Studio output (A_MODES line)");
                    }
                }
                for (j, a) in st.atoms.iter().enumerate() {
                    w.line("");
                    w.line("[[phase.atom]]");
                    w.kv("label", &string(&a.label), "");
                    w.kv("species", &string(&a.species), "");
                    if let Some(q) = a.charge {
                        w.kv("charge", &q.to_string(), "# formal oxidation state, for bond valences (the species keeps its scattering)");
                    }
                    let grouped = ph.rigid.iter().enumerate().find_map(|(g, gr)| gr.members.iter().position(|m| m.index == j).map(|k| (g, gr, k)));
                    match (&ph.modes, grouped) {
                        (Some(ms), _) => w.kv("atom_position", &floats(&ms.reference[j]), "# fractions of the cell: the reference position, the modes add their displacement"),
                        (None, Some(_)) => w.kv("atom_position", &floats(&a.x), "# fractions of the cell: where its rigid group places it (computed on reading)"),
                        (None, None) => w.kv("atom_position", &floats(&a.x), &unit("atom_position")),
                    }
                    // A rigid group's atom (P41): its group and internal coordinates, lengths and angles apart.
                    if let Some((g, gr, k)) = grouped {
                        let (len, ang) = gr.coordinates.split();
                        let m = &gr.members[k];
                        w.kv("rigid_group", &(g + 1).to_string(), &format!("# its atoms placed by [[phase.rigid_group]] {}", g + 1));
                        if !len.is_empty() {
                            w.kv("rigid_length", &floats(&len.iter().map(|&c| m.internal[c]).collect::<Vec<_>>()), "# Å");
                        }
                        if !ang.is_empty() {
                            w.kv("rigid_angle", &floats(&ang.iter().map(|&c| m.internal[c]).collect::<Vec<_>>()), "# deg");
                        }
                        if gr.coordinates == crate::rigid_body::Coordinates::ZMatrix && k > 0 {
                            w.kv("rigid_connect", &format!("[{}]", m.connect.iter().map(|c| c.to_string()).collect::<Vec<_>>().join(", ")), "# the Z-matrix's N1 N2 N3, by their order in the group");
                        }
                    }
                    w.num("atom_biso", a.biso);
                    w.num("atom_occupancy", a.occ);
                    if let Some(b) = a.beta {
                        w.kv("atom_betas", &floats(&b), "# β11 β22 β33 β12 β13 β23");
                    }
                }
                for gr in &ph.rigid {
                    w.line("");
                    w.line("[[phase.rigid_group]]");
                    w.kv("kind", &string(gr.coordinates.name()), "# spherical, cartesian, cylindrical or z-matrix: its atoms' internal coordinates");
                    w.kv("rigid_origin", &floats(&gr.origin), "# fractions of the cell: X0 Y0 Z0");
                    w.kv("rigid_orientation", &floats(&[gr.theta, gr.phi, gr.chi]), "# deg: θ φ χ, R = Rz(φ) Ry(θ) Rz(χ)");
                    if let Some(sat) = gr.satellite {
                        w.kv("satellite", &format!("[{}, {}, {}]", sat.parent + 1, sat.atom, sat.from), "# [parent group, its atom N1, its atom N2 (0: its origin)]: origin, θ and φ follow");
                    }
                    w.kv("fullprof_option", &float(gr.option), "# FullProf's P6, as read");
                    if let Some(t) = &gr.displacements {
                        w.kv("rigid_tls_t", &floats(&t.t), "# Å²: T11 T22 T33 T12 T13 T23");
                        w.kv("rigid_tls_l", &floats(&t.l), "# rad²: L11 L22 L33 L12 L13 L23");
                        w.kv("rigid_tls_s", &floats(&t.s), "# Å·rad: S11 S22 S33 S12 S13 S23 S21 S31 S32");
                        w.kv("tls_centre", &floats(&t.centre), "# fractions: the centre T, L, S refer to, less rigid_origin");
                    }
                }
                for mode in ph.modes.iter().flat_map(|ms| &ms.modes) {
                    w.line("");
                    w.line("[[phase.mode]]");
                    w.kv("name", &string(&mode.name), "");
                    w.kv("irrep", &string(&mode.irrep), "");
                    w.num("mode_amplitude", mode.amplitude);
                    let rows_v: Vec<Vec<String>> = mode.vectors.iter().map(|v| {
                        let mut r = vec![string(&atom_name(st, v.atom))];
                        r.extend(v.v.iter().map(|x| float(*x)));
                        r.push(float(v.coeff));
                        r
                    }).collect();
                    w.kv("mode_vector", &rows(&rows_v), "# [atom, Vx, Vy, Vz (fractions of the cell), coefficient]");
                }
            }
        }
    }
    for r in &m.restraints {
        let Intensities::Structure(st, _) = &m.phases[r.phase].intensities else { continue };
        w.line("");
        w.line("[[restraint]]");
        w.kv("phase", &(r.phase + 1).to_string(), "");
        match &r.kind {
            Kind::Distance { atoms, op } => {
                w.kv("distance", &format!("[{}]", atoms.iter().map(|&a| string(&atom_name(st, a))).collect::<Vec<_>>().join(", ")), "");
                w.kv("operator", &operator_text(op), "# applied to the second atom");
                w.kv("target", &float(r.target), "# Å");
                w.kv("sigma", &float(r.sigma), "# Å");
            }
            Kind::Angle { atoms, ops } => {
                w.kv("angle", &format!("[{}]", atoms.iter().map(|&a| string(&atom_name(st, a))).collect::<Vec<_>>().join(", ")), "");
                w.kv("operators", &format!("[{}, {}]", operator_text(&ops[0]), operator_text(&ops[1])), "# applied to the second and third atoms");
                w.kv("target", &float(r.target), "# deg");
                w.kv("sigma", &float(r.sigma), "# deg");
            }
        }
    }
    let rf = &c.refinement;
    w.line("");
    w.line("[refine]");
    if let Some(n) = rf.cycles {
        w.kv("cycles", &n.to_string(), "");
    }
    if let Some(e) = rf.eps {
        w.kv("convergence", &float(e), "# every |shift| below this times its esd");
    }
    let st = &rf.steering;
    let dst = Steering::default();
    if st.lambda != dst.lambda {
        w.kv("lambda", &float(st.lambda), "# λ's start, the Marquardt factor of the first cycle (0: Gauss–Newton first)");
    }
    if st.marquardt != dst.marquardt {
        w.kv("marquardt", &st.marquardt.to_string(), "# false: every step applied, as FullProf applies them");
    }
    if st.class_damping != dst.class_damping {
        let d: Vec<String> = Class::ALL.iter().map(|c| format!("{} = {}", c.word(), float(st.class_damping[c.index()]))).collect();
        w.kv("damping", &format!("{{ {} }}", d.join(", ")), "# each class's shifts multiplied by this (FullProf's R_at, R_an, R_pr, R_gl)");
    }
    if st.guards != dst.guards {
        let g = &st.guards;
        w.kv("divergence", &format!("{{ rising = {}, jump = {}, jump_from = {}, singular = {} }}", g.rising, float(g.jump), g.jump_from, g.singular), "# 0 turns a guard off");
    }
    w.line("# Each parameter: the values it moves and their multipliers (value = multiplier × parameter); it is named by");
    w.line("# its first value. The values the symmetry ties to these follow by themselves.");
    w.kv(
        "parameters",
        &rows(&rf.parameters.iter().map(|t| t.iter().flat_map(|&(f, c)| [string(&address(f, m)), float(c)]).collect()).collect::<Vec<_>>()),
        "",
    );
    // Limits, and the behaviour at a bound (with the dictionary's bounds where the user gives none).
    let mut bounded: Vec<Field> = rf.limits.iter().map(|l| l.0).collect();
    bounded.extend(rf.behaviour.iter().map(|b| b.0).filter(|f| !rf.limits.iter().any(|l| l.0 == *f)));
    if !bounded.is_empty() {
        let l: Vec<Vec<String>> = bounded
            .iter()
            .map(|&f| {
                let mut v = vec![format!("value = {}", string(&address(f, m)))];
                if let Some(&(_, lo, hi)) = rf.limits.iter().find(|l| l.0 == f) {
                    if let Some(lo) = lo {
                        v.push(format!("min = {}", float(lo)));
                    }
                    if let Some(hi) = hi {
                        v.push(format!("max = {}", float(hi)));
                    }
                }
                if let Some(b) = rf.behaviour.iter().find(|b| b.0 == f) {
                    v.push(format!("behaviour = {}", string(b.1.word())));
                }
                v
            })
            .collect();
        let mut s = String::from("[\n");
        for r in l {
            s.push_str(&format!("  {{ {} }},\n", r.join(", ")));
        }
        s.push(']');
        w.kv("limits", &s, "# behaviour at a bound: stop (the default), reflect or periodic");
    }
    // Per parameter (named by its first value) its damping, per value its largest shift per cycle.
    let mut steered: Vec<String> = st.damping.iter().map(|d| d.0.clone()).collect();
    for (f, _) in &st.max_shift {
        let a = address(*f, m);
        if !steered.contains(&a) {
            steered.push(a);
        }
    }
    if !steered.is_empty() {
        let mut s = String::from("[\n");
        for a in &steered {
            let mut v = vec![format!("value = {}", string(a))];
            if let Some(d) = st.damping.iter().find(|d| &d.0 == a) {
                v.push(format!("damping = {}", float(d.1)));
            }
            if let Some(x) = st.max_shift.iter().find(|x| &address(x.0, m) == a) {
                v.push(format!("max_shift = {}", float(x.1)));
            }
            s.push_str(&format!("  {{ {} }},\n", v.join(", ")));
        }
        s.push(']');
        w.kv("steering", &s, "# damping of a parameter (its shifts × damping, never its esd); largest shift of a value per cycle");
    }
    for l in &rf.links {
        w.line("");
        w.line("[[refine.link]]");
        w.kv("value", &string(&address(l.value, m)), "# = coefficient × follows + offset");
        w.kv("follows", &string(&address(l.follows, m)), "");
        w.kv("coefficient", &float(l.coefficient), "");
        w.kv("offset", &float(l.offset), "");
    }
    for v in &rf.variables {
        w.line("");
        w.line("[[refine.variable]]");
        w.kv("name", &string(&v.name), "");
        w.kv("value", &float(v.value), "# at the values written above");
        w.kv("terms", &rows(&v.terms.iter().map(|&(f, c)| vec![string(&address(f, m)), float(c)]).collect::<Vec<_>>()), "# [value, how much it moves per unit of the variable]");
    }
    for p in &rf.plan {
        w.line("");
        w.line("[[refine.step]]");
        w.kv("name", &string(&p.name), "");
        w.kv("free", &format!("[{}]", p.free.iter().map(|x| string(x)).collect::<Vec<_>>().join(", ")), "# addresses, prefixes ending in *, variables, classes, or \"all\"");
        if let Some(d) = p.damping {
            w.kv("damping", &float(d), "");
        }
        if let Some(n) = p.cycles {
            w.kv("cycles", &n.to_string(), "");
        }
        if let Some(e) = p.eps {
            w.kv("convergence", &float(e), "");
        }
    }
    for l in &rf.linear {
        w.line("");
        w.line("[[refine.linear]]");
        w.kv("name", &string(&l.name), "");
        w.kv("terms", &rows(&l.terms.iter().map(|(n, c)| vec![string(n), float(*c)]).collect::<Vec<_>>()), "# [parameter, coefficient]");
        w.kv("target", &float(l.target), "");
        w.kv("sigma", &float(l.sigma), "");
    }
    w.line("");
    w.line("[convention]");
    w.line("# true: FullProf's rule where it departs from the papers (the compatibility layer); false: the papers' (L050).");
    let cv = m.convention;
    for (k, v) in [
        ("lorentz_polarisation", cv.lorentz_polarisation),
        ("berar_baldinozzi", cv.berar_baldinozzi),
        ("axial_divergence", cv.axial_divergence),
        ("preferred_orientation", cv.preferred_orientation),
        ("excluded_points", cv.excluded_points),
        ("peak_window", cv.peak_window),
        ("reflections", cv.reflections),
        ("microstructure", cv.microstructure),
        ("angle_restraints", cv.angle_restraints),
        ("window_edge", cv.window_edge),
        ("background_edge", cv.background_edge),
        ("displacement", cv.displacement),
        ("background_bounds", cv.background_bounds),
    ] {
        w.kv(k, if v { "true" } else { "false" }, "");
    }
    w.out
}

// The reader ---------------------------------------------------------------------------------------------------------

struct R {
    diags: Vec<Diagnostic>,
}

/// The nearest of the known keys (by edit distance), for "did you mean".
fn nearest<'k>(key: &str, known: &[&'k str]) -> Option<&'k str> {
    let d = |a: &str, b: &str| {
        let (a, b): (Vec<char>, Vec<char>) = (a.chars().collect(), b.chars().collect());
        let mut prev: Vec<usize> = (0..=b.len()).collect();
        for i in 1..=a.len() {
            let mut cur = vec![i; b.len() + 1];
            for j in 1..=b.len() {
                cur[j] = (prev[j] + 1).min(cur[j - 1] + 1).min(prev[j - 1] + usize::from(a[i - 1] != b[j - 1]));
            }
            prev = cur;
        }
        prev[b.len()]
    };
    known.iter().map(|k| (d(key, k), *k)).filter(|(n, k)| *n <= 3.max(k.len() / 3)).min().map(|x| x.1)
}

impl R {
    fn err(&mut self, span: Span, msg: impl Into<String>) {
        self.diags.push(Diagnostic { span, message: msg.into(), warning: false, code: "value", replacement: None });
    }

    /// A name that is not known, with the nearest known one as its replacement.
    fn unknown(&mut self, span: Span, msg: String, nearest: Option<&str>) {
        let hint = nearest.map_or(String::new(), |n| format!("; did you mean \"{n}\"?"));
        self.diags.push(Diagnostic { span, message: format!("{msg}{hint}"), warning: false, code: "unknown-name", replacement: nearest.map(str::to_string) });
    }

    fn check_keys(&mut self, t: &Table, known: &[&str]) {
        for e in &t.entries {
            if !known.contains(&e.key.as_str()) {
                self.unknown(e.key_span, format!("\"{}\" is not a key of {}", e.key, t.name()), nearest(&e.key, known));
            }
        }
    }

    fn missing(&mut self, t: &Table, key: &str) {
        let span = t.header.unwrap_or_default();
        self.diags.push(Diagnostic { span, message: format!("{} needs \"{key}\"", t.name()), warning: false, code: "missing", replacement: None });
    }

    fn num_item(&mut self, it: &Item, what: &str) -> Option<f64> {
        match it.value {
            Value::Float(v) => Some(v),
            Value::Integer(v) => Some(v as f64),
            ref v => {
                self.err(it.span, format!("{what} must be a number, not {}", v.kind()));
                None
            }
        }
    }

    fn f(&mut self, t: &Table, key: &str, default: Option<f64>) -> f64 {
        match t.get(key) {
            Some(e) => self.num_item(&e.value, &format!("\"{key}\"")).unwrap_or(0.0),
            None => default.unwrap_or_else(|| {
                self.missing(t, key);
                0.0
            }),
        }
    }

    fn s(&mut self, t: &Table, key: &str) -> Option<(String, Span)> {
        let e = t.get(key)?;
        match &e.value.value {
            Value::String(s) => Some((s.clone(), e.value.span)),
            v => {
                self.err(e.value.span, format!("\"{key}\" must be a quoted string, not {}", v.kind()));
                None
            }
        }
    }

    fn b(&mut self, t: &Table, key: &str, default: bool) -> bool {
        match t.get(key) {
            None => default,
            Some(e) => match e.value.value {
                Value::Bool(v) => v,
                ref v => {
                    self.err(e.value.span, format!("\"{key}\" must be true or false, not {}", v.kind()));
                    default
                }
            },
        }
    }

    fn array<'a>(&mut self, it: &'a Item, what: &str) -> Option<&'a [Item]> {
        match &it.value {
            Value::Array(v) => Some(v),
            v => {
                self.err(it.span, format!("{what} must be an array [ … ], not {}", v.kind()));
                None
            }
        }
    }

    /// An array of numbers, of length `n` when given.
    fn nums(&mut self, it: &Item, what: &str, n: Option<usize>) -> Option<Vec<f64>> {
        let a = self.array(it, what)?;
        if let Some(n) = n {
            if a.len() != n {
                self.err(it.span, format!("{what} must hold {n} numbers, not {}", a.len()));
                return None;
            }
        }
        a.iter().map(|x| self.num_item(x, what)).collect()
    }

    fn nums_key(&mut self, t: &Table, key: &str, n: Option<usize>, default: Option<Vec<f64>>) -> Vec<f64> {
        match t.get(key) {
            Some(e) => self.nums(&e.value, &format!("\"{key}\""), n).unwrap_or_else(|| vec![0.0; n.unwrap_or(0)]),
            None => default.unwrap_or_else(|| {
                self.missing(t, key);
                vec![0.0; n.unwrap_or(0)]
            }),
        }
    }

    /// An array of rows of `n` numbers.
    fn table_rows(&mut self, t: &Table, key: &str, n: usize) -> Vec<Vec<f64>> {
        let Some(e) = t.get(key) else { return Vec::new() };
        let Some(a) = self.array(&e.value, &format!("\"{key}\"")) else { return Vec::new() };
        a.iter().filter_map(|r| self.nums(r, &format!("a row of \"{key}\""), Some(n))).collect()
    }

    fn choice<T: Copy>(&mut self, t: &Table, key: &str, options: &[(T, &str)], default: T) -> T {
        let Some((v, span)) = self.s(t, key) else { return default };
        match options.iter().find(|o| o.1 == v) {
            Some(o) => o.0,
            None => {
                self.err(span, format!("\"{v}\" is not one of: {}", options.iter().map(|o| o.1).collect::<Vec<_>>().join(", ")));
                default
            }
        }
    }
}

fn index_of_atom(st: &Structure, name: &str) -> Option<usize> {
    if let Some((l, n)) = name.rsplit_once(':') {
        if let Ok(n) = n.parse::<usize>() {
            if st.atoms.get(n.wrapping_sub(1)).is_some_and(|a| a.label == l) {
                return Some(n - 1);
            }
        }
    }
    st.atoms.iter().position(|a| a.label == name)
}

fn operator(r: &mut R, it: &Item) -> Option<Operator> {
    let Value::Table(entries) = &it.value else {
        r.err(it.span, "an operator is an inline table { rotation = \"x,y,z\", translation = [0.0, 0.0, 0.0] }");
        return None;
    };
    let get = |k: &str| entries.iter().find(|e| e.key == k);
    let rot = match get("rotation").map(|e| &e.value) {
        Some(Item { value: Value::String(s), span }) => match SymOp::parse_xyz(s) {
            Ok(op) if op.t == [0; 3] => op.r,
            Ok(_) => {
                r.err(*span, "the rotation has no translation: give it in \"translation\"");
                return None;
            }
            Err(e) => {
                r.err(*span, format!("the rotation \"{s}\" does not read: {e}"));
                return None;
            }
        },
        _ => {
            r.err(it.span, "the operator needs rotation = \"x,y,z\"");
            return None;
        }
    };
    let t = match get("translation") {
        Some(e) => r.nums(&e.value, "the translation", Some(3))?,
        None => vec![0.0; 3],
    };
    Some(Operator { rot, t: [t[0], t[1], t[2]] })
}

/// A value named by its address in the control file.
fn addr_item(r: &mut R, it: &Item, m: &Model) -> Option<Field> {
    match &it.value {
        Value::String(s) => field_of(s, m).map_err(|e| r.err(it.span, e)).ok(),
        v => {
            r.err(it.span, format!("a value is named by its address (\"phase1.cell_a\"), not {}", v.kind()));
            None
        }
    }
}

/// A damping is above 0; above 1 the step overshoots (over-relaxation), which is allowed with a warning.
fn check_damping(r: &mut R, span: Span, d: f64) {
    if d <= 0.0 {
        r.err(span, "a damping is above 0 (1: the whole shift)");
    } else if d > 1.0 {
        r.diags.push(Diagnostic { span, message: format!("damping {d} is above 1: each step overshoots by {:.0} %, and χ² may rise (over-relaxation)", (d - 1.0) * 100.0), warning: true, code: "steering", replacement: None });
    }
}

/// Read a parsed control file into a control (None when an error leaves nothing to build); the diagnostics are the
/// document's and the schema's.
pub fn read(doc: &Document) -> (Option<Control>, Vec<Diagnostic>) {
    read_with(doc, None)
}

/// As `read`, with the folder of the control file: a background-points file it names (`file = …`) is read from there.
/// Without it the points of such a file are left out (and the refinement is not checked against them).
pub fn read_with(doc: &Document, base: Option<&std::path::Path>) -> (Option<Control>, Vec<Diagnostic>) {
    let mut r = R { diags: doc.diagnostics.clone() };
    let mut background_file: Option<String> = None;
    let mut points_unread = false;
    let root = &doc.tables[0];
    r.check_keys(root, ROOT_KEYS);
    match r.s(root, "format") {
        Some((f, _)) if f == FORMAT => {}
        Some((f, span)) => r.err(span, format!("this file is \"{f}\"; this version of the core reads \"{FORMAT}\"")),
        None => r.err(Span { start: 0, end: 0, line: 1, column: 1 }, format!("the file must start with format = \"{FORMAT}\"")),
    }
    let title = r.s(root, "title").map_or(String::new(), |s| s.0);
    let neutrons = r.choice(root, "radiation", &[(false, "x-ray"), (true, "neutron")], false);
    let find = |path: &[&str]| doc.tables.iter().find(|t| !t.array && t.path == path);
    let empty = Table { path: vec!["pattern".into()], array: false, header: None, header_lines: (0, 0), entries: Vec::new() };
    let pt = match find(&["pattern"]) {
        Some(t) => t,
        None => {
            r.err(Span { start: 0, end: 0, line: 1, column: 1 }, "the file needs a [pattern] table (the wavelength at least)");
            &empty
        }
    };
    r.check_keys(pt, PATTERN_KEYS);
    let data_file = r.s(pt, "data").map(|s| s.0);
    let reader = r.s(pt, "reader");
    if let Some((id, span)) = &reader {
        if crate::formats::by_id(id).is_none() {
            let ids: Vec<&str> = crate::formats::FORMATS.iter().map(|f| f.id).collect();
            r.unknown(*span, format!("\"{id}\" is not a reader of the registry (drs data formats)"), nearest(id, &ids));
        }
    }
    // The reader's options (the app's request 51), checked against the reader's declarations.
    let mut reader_options = crate::formats::options::Options::new();
    if let Some(e) = pt.get("reader_options") {
        use crate::formats::options::OptionValue;
        match &e.value.value {
            Value::Table(entries) => {
                for x in entries {
                    let v = match &x.value.value {
                        Value::Integer(n) => OptionValue::Integer(*n),
                        Value::Float(f) => OptionValue::Number(*f),
                        Value::String(t) => OptionValue::Text(t.clone()),
                        v => {
                            r.err(x.value.span, format!("a reader's option is a number or a quoted word, not {}", v.kind()));
                            continue;
                        }
                    };
                    reader_options.insert(x.key.clone(), v);
                }
                match reader.as_ref().and_then(|(id, _)| crate::formats::by_id(id)) {
                    Some(f) => {
                        if let Err(m) = crate::formats::options::validate(f, &reader_options) {
                            r.err(e.value.span, m);
                        }
                    }
                    None if reader.is_none() => r.err(e.key_span, "reader_options needs the reader they are for: reader = \"…\""),
                    None => {}
                }
            }
            _ => r.err(e.value.span, "reader_options is { option = value, … } (drs data formats lists each reader's)"),
        }
    }
    let l1 = r.f(pt, "wavelength_1", None);
    let ilo = match r.s(pt, "geometry") {
        None => 0,
        Some((g, span)) => match GEOMETRY.iter().find(|x| x.1 == g) {
            Some(x) => x.0,
            None => match g.strip_prefix("fullprof-ilo-").and_then(|n| n.parse().ok()) {
                Some(n) => n,
                None => {
                    r.err(span, format!("\"{g}\" is not one of: {}", GEOMETRY.iter().map(|x| x.1).collect::<Vec<_>>().join(", ")));
                    0
                }
            },
        },
    };
    let tof = pt.get("tof_dtt1").map(|_| {
        let mut b = crate::tof::Bank {
            dtt1: 0.0,
            dtt2: 0.0,
            dtt_1overd: 0.0,
            two_theta: r.f(pt, "tof_bank_two_theta", None),
            thermal: None,
            absorption: r.choice(pt, "tof_absorption_correction", TOF_ABSORPTION, crate::tof::Absorption::FlatPlate),
            instrument: None,
            tmax: r.f(pt, "tof_max", None),
        };
        let n = if pt.get("tof_thermal_dtt1").is_some() { 8 } else { 3 };
        for (k, id) in crate::tof::BANK_FIELDS.iter().enumerate().take(n) {
            b.set(k, r.f(pt, id, if k == 0 { None } else { Some(0.0) }));
        }
        b
    });
    if tof.is_none() {
        for k in &PATTERN_KEYS[PATTERN_KEYS.len() - 11..] {
            if let Some(e) = pt.get(k) {
                r.err(e.key_span, format!("{k} belongs to a time-of-flight pattern, one with tof_dtt1"));
            }
        }
    }
    let mut model = Model {
        convention: Convention::DEFAULT,
        neutrons,
        lambda: [l1, r.f(pt, "wavelength_2", Some(l1))],
        ratio: r.f(pt, "wavelength_ratio", Some(0.0)),
        zero: r.f(pt, "zero_shift", Some(0.0)),
        sycos: r.f(pt, "sample_displacement", Some(0.0)),
        sysin: r.f(pt, "sample_transparency", Some(0.0)),
        cthm: r.f(pt, "monochromator_polarisation", Some(0.0)),
        rpolarz: r.f(pt, "polarisation_fraction", Some(0.0)),
        ilo,
        mur: r.f(pt, "absorption_mur", Some(0.0)),
        wdt: r.f(pt, "peak_window", Some(8.0)),
        asym_limit: r.f(pt, "asymmetry_limit", Some(0.0)),
        march: r.choice(pt, "preferred_orientation_function", &[(true, "march-dollase"), (false, "exponential")], true),
        background: BackgroundModel::Polynomial { origin: 0.0, coefficients: Vec::new() },
        instrumental: InstrumentalWidths::None,
        phases: Vec::new(),
        excluded: r.table_rows(pt, "excluded_regions", 2).into_iter().map(|v| (v[0], v[1])).collect(),
        weighting: r.choice(pt, "weighting", WEIGHTS, Weighting::Variance),
        restraints: Vec::new(),
        fullprof_step: pt.get("fullprof_step").map(|_| r.f(pt, "fullprof_step", None)),
        microabsorption: [r.f(pt, "microabsorption_p0", Some(0.0)), r.f(pt, "microabsorption_cp", Some(0.0)), r.f(pt, "microabsorption_tau", Some(0.0))],
        tof,
        background_terms: Default::default(),
        background_ceiling: pt.get("background_ceiling").map(|_| r.f(pt, "background_ceiling", None)),
    };
    if let Some(bt) = find(&["pattern", "background"]) {
        r.check_keys(bt, BACKGROUND_KEYS);
        model.background = match r.s(bt, "function").map(|s| (s.0, s.1)) {
            Some((f, _)) if f == "points" => {
                // The points: in a file the control file names (ours or FullProf's .bac), or written here (R26).
                let stated = r.s(bt, "interpolation").map(|(s, span)| (crate::background::Interpolation::from_name(&s), span));
                let mut how = match stated {
                    Some((Some(h), _)) => Some(h),
                    Some((None, span)) => {
                        r.err(span, "interpolation is \"linear\" or \"spline\"");
                        None
                    }
                    None => None,
                };
                let points = match r.s(bt, "file") {
                    Some((file, span)) => {
                        background_file = Some(file.clone());
                        match base {
                            Some(dir) => match crate::background::read_path(&dir.join(&file)) {
                                Ok(pf) => {
                                    how = how.or(pf.interpolation);
                                    pf.points
                                }
                                Err(e) => {
                                    r.diags.push(Diagnostic { span, message: format!("the background file does not read: {e}"), warning: false, code: "missing", replacement: None });
                                    Vec::new()
                                }
                            },
                            None => {
                                points_unread = true;
                                Vec::new()
                            }
                        }
                    }
                    None => r.table_rows(bt, "background_points", 2).into_iter().map(|v| (v[0], v[1])).collect(),
                };
                match how.unwrap_or(crate::background::Interpolation::Linear) {
                    crate::background::Interpolation::Linear => BackgroundModel::Points(points),
                    crate::background::Interpolation::Spline => BackgroundModel::Spline(points),
                }
            }
            Some((f, _)) if f == "polynomial" => BackgroundModel::Polynomial { origin: r.f(bt, "background_origin", Some(0.0)), coefficients: r.nums_key(bt, "background_polynomial", None, Some(Vec::new())) },
            Some((f, _)) if f == "chebyshev" => BackgroundModel::Chebyshev { from: r.f(bt, "from", None), to: r.f(bt, "to", None), coefficients: r.nums_key(bt, "background_chebyshev", None, Some(Vec::new())) },
            Some((f, _)) if f == "fourier" => BackgroundModel::Fourier { coefficients: r.nums_key(bt, "background_fourier", None, Some(Vec::new())) },
            Some((f, span)) => {
                r.err(span, format!("\"{f}\" is not one of: points, polynomial, chebyshev, fourier"));
                model.background.clone()
            }
            None => {
                r.missing(bt, "function");
                model.background.clone()
            }
        };
        // Further terms (P23): Debye-like terms, inverse powers of 2θ, profiles read from files.
        let span = |key: &str| bt.get(key).map(|e| e.value.span).unwrap_or_default();
        let amplitude = r.nums_key(bt, "background_debye_amplitude", None, Some(Vec::new()));
        let distance = r.nums_key(bt, "background_debye_distance", None, Some(Vec::new()));
        if amplitude.len() != distance.len() {
            r.err(span("background_debye_distance"), format!("background_debye_distance must hold a distance for each of the {} amplitudes", amplitude.len()));
        }
        let terms = &mut model.background_terms;
        terms.debye = amplitude.iter().zip(&distance).map(|(&a, &d)| crate::background_terms::DebyeTerm { amplitude: a, distance: d }).collect();
        terms.inverse = r.nums_key(bt, "background_inverse_power", None, Some(Vec::new()));
        if terms.inverse.len() > 3 {
            r.err(span("background_inverse_power"), format!("background_inverse_power holds at most 3 coefficients, not {}", terms.inverse.len()));
            terms.inverse.truncate(3);
        }
        let mut files = Vec::new();
        if let Some(e) = bt.get("profile_files") {
            if let Some(a) = r.array(&e.value, "\"profile_files\"") {
                for it in a {
                    match &it.value {
                        Value::String(s) => files.push((s.clone(), it.span)),
                        v => r.err(it.span, format!("a file of \"profile_files\" must be a quoted string, not {}", v.kind())),
                    }
                }
            }
        }
        let coefficients = r.nums_key(bt, "background_profile", None, Some(vec![1.0; files.len()]));
        if coefficients.len() != files.len() {
            r.err(span("background_profile"), format!("background_profile must hold a coefficient for each of the {} profile files", files.len()));
        }
        for ((file, fspan), &c) in files.into_iter().zip(&coefficients) {
            let points = match base {
                Some(dir) => match crate::background::read_path(&dir.join(&file)) {
                    Ok(pf) => pf.points,
                    Err(e) => {
                        r.diags.push(Diagnostic { span: fspan, message: format!("the background profile does not read: {e}"), warning: false, code: "missing", replacement: None });
                        Vec::new()
                    }
                },
                None => {
                    points_unread = true;
                    Vec::new()
                }
            };
            model.background_terms.profiles.push(crate::background_terms::Profile { file, coefficient: c, points });
        }
        // A file's b1–b4 (eq. 3.6; the app's request 47): with points, the four or none.
        let keys = ["background_file_scale", "background_file_offset", "background_file_stretch", "background_file_shift"];
        if keys.iter().any(|k| bt.get(k).is_some()) {
            if !matches!(model.background, BackgroundModel::Points(_)) {
                r.err(span(keys.iter().find(|k| bt.get(k).is_some()).unwrap()), "background_file_scale … _shift move a background of points joined linearly (FullProf's Nba = 1)");
            }
            let d = crate::background_terms::FileScale::default();
            let f = crate::background_terms::FileScale {
                file: background_file.clone().unwrap_or_default(),
                scale: r.f(bt, keys[0], Some(d.scale)),
                offset: r.f(bt, keys[1], Some(d.offset)),
                stretch: r.f(bt, keys[2], Some(d.stretch)),
                shift: r.f(bt, keys[3], Some(d.shift)),
            };
            if 1.0 + f.stretch <= 0.0 {
                r.err(span(keys[2]), format!("background_file_stretch {} reads the points nowhere (1 + b3 ≤ 0)", f.stretch));
            }
            model.background_terms.file = Some(f);
        }
    }
    if let Some(rt) = find(&["pattern", "resolution"]) {
        r.check_keys(rt, RESOLUTION_KEYS);
        match rt.get("function").map(|e| &e.value) {
            Some(Item { value: Value::String(s), span }) if s == "tof" => {
                let v = r.nums_key(rt, "instrument_resolution", Some(22), None);
                match model.tof.as_mut() {
                    Some(b) => b.instrument = Some(crate::tof::Instrument::from_values(&v)),
                    None => r.err(*span, "a \"tof\" resolution belongs to a time-of-flight pattern, one with tof_dtt1"),
                }
            }
            Some(Item { value: Value::String(s), .. }) if s == "table" => model.instrumental = InstrumentalWidths::Table(r.table_rows(rt, "instrument_resolution", 3).into_iter().map(|v| [v[0], v[1], v[2]]).collect()),
            Some(Item { value: Value::Integer(res), span }) => {
                if !(1..=3).contains(res) {
                    r.err(*span, "the resolution function is 1, 2 or 3 (or \"table\")");
                }
                let v = r.nums_key(rt, "instrument_resolution", Some(6), None);
                model.instrumental = InstrumentalWidths::Function { res: *res as i32, uvwxyz: [v[0], v[1], v[2], v[3], v[4], v[5]] };
            }
            Some(it) => r.err(it.span, "\"function\" is 1, 2, 3, \"table\" or \"tof\""),
            None => r.missing(rt, "function"),
        }
    }
    let mut fullprof = super::FullProfRecord::default();
    if let Some(ft) = find(&["pattern", "fullprof"]) {
        r.check_keys(ft, FULLPROF_KEYS);
        if ft.get("range").is_some() {
            let v = r.nums_key(ft, "range", Some(3), None);
            fullprof.range = (v.len() == 3).then(|| [v[0], v[1], v[2]]);
        }
        fullprof.excluded = r.table_rows(ft, "excluded_regions", 2).into_iter().map(|v| (v[0], v[1])).collect();
    }
    let mut operations = Vec::new();
    for t in doc.tables.iter().filter(|t| t.array && t.path == ["pattern", "operation"]) {
        match Operation::from_table(t) {
            Ok(op) => operations.push(op),
            Err(d) => r.diags.push(d),
        }
    }
    // The user's scattering factors (P11): matched to the atoms' species case aside.
    let mut user = Vec::new();
    for t in doc.tables.iter().filter(|t| t.array && t.path == ["pattern", "scatterer"]) {
        r.check_keys(t, SCATTERER_KEYS);
        let Some((species, _)) = r.s(t, "species") else {
            r.missing(t, "species");
            continue;
        };
        let form_factor = t.get("form_factor").and_then(|e| r.nums(&e.value, "\"form_factor\"", Some(9))).map(|k| [k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8]]);
        let form_factor_table = r.table_rows(t, "form_factor_table", 2).into_iter().map(|v| (v[0], v[1])).collect();
        let scattering_length = t.get("scattering_length").map(|_| r.f(t, "scattering_length", None));
        let scattering_length_imaginary = r.f(t, "scattering_length_imaginary", Some(0.0));
        user.push(crate::scattering::UserScatterer { species, scattering_length, scattering_length_imaginary, f_prime: r.f(t, "f_prime", Some(0.0)), f_double_prime: r.f(t, "f_double_prime", Some(0.0)), form_factor, form_factor_table });
    }
    // Phases and their atoms.
    let radiation = if neutrons { Radiation::Neutron } else { Radiation::Xray { lambda: l1, anomalous: true } };
    let mut phase_tables: Vec<PhaseTables> = Vec::new();
    for t in &doc.tables {
        if t.array && t.path == ["phase"] {
            phase_tables.push((t, Vec::new(), Vec::new(), Vec::new()));
        } else if t.array && t.path == ["phase", "atom"] {
            match phase_tables.last_mut() {
                Some(p) => p.1.push(t),
                None => r.err(t.header.unwrap_or_default(), "[[phase.atom]] comes after the [[phase]] it belongs to"),
            }
        } else if t.array && t.path == ["phase", "mode"] {
            match phase_tables.last_mut() {
                Some(p) => p.2.push(t),
                None => r.err(t.header.unwrap_or_default(), "[[phase.mode]] comes after the [[phase]] it belongs to"),
            }
        } else if t.array && t.path == ["phase", "rigid_group"] {
            match phase_tables.last_mut() {
                Some(p) => p.3.push(t),
                None => r.err(t.header.unwrap_or_default(), "[[phase.rigid_group]] comes after the [[phase]] it belongs to"),
            }
        } else if t.header.is_some() && !matches!(t.path.iter().map(String::as_str).collect::<Vec<_>>()[..], ["pattern"] | ["pattern", "background"] | ["pattern", "resolution"] | ["pattern", "fullprof"] | ["pattern", "operation"] | ["pattern", "scatterer"] | ["phase"] | ["phase", "atom"] | ["phase", "mode"] | ["phase", "rigid_group"] | ["restraint"] | ["refine"] | ["refine", "linear"] | ["refine", "link"] | ["refine", "variable"] | ["refine", "step"] | ["convention"]) {
            let known = ["pattern", "pattern.background", "pattern.resolution", "pattern.fullprof", "pattern.operation", "pattern.scatterer", "phase", "phase.atom", "phase.mode", "phase.rigid_group", "restraint", "refine", "refine.linear", "refine.link", "refine.variable", "refine.step", "convention"];
            let name = t.path.join(".");
            let near = nearest(&name, &known);
            let hint = near.map_or(String::new(), |n| format!("; did you mean \"{n}\"?"));
            let replacement = near.map(|n| if t.array { format!("[[{n}]]") } else { format!("[{n}]") });
            r.diags.push(Diagnostic { span: t.header.unwrap_or_default(), message: format!("{} is not a table of this format{hint}", t.name()), warning: false, code: "unknown-name", replacement });
        }
    }
    for (pt, atoms, mode_tables, rigid_tables) in &phase_tables {
        r.check_keys(pt, PHASE_KEYS);
        let name = r.s(pt, "name").map_or(String::new(), |s| s.0);
        let group = match (r.s(pt, "hall"), r.s(pt, "space_group")) {
            (Some((h, span)), _) => SpaceGroup::from_hall(&h).map_err(|e| r.err(span, format!("the Hall symbol \"{h}\" does not read: {e}"))).ok(),
            (None, Some((s, span))) => SpaceGroup::from_fullprof_symbol(&s).map_err(|e| r.err(span, format!("the space group \"{s}\" does not read: {e}"))).ok(),
            (None, None) => {
                r.missing(pt, "hall");
                None
            }
        };
        let group = match (group, pt.get("operators")) {
            (g, Some(e)) => {
                let ops: Option<Vec<SymOp>> = r.array(&e.value, "\"operators\"").map(|a| {
                    a.iter()
                        .filter_map(|it| match &it.value {
                            Value::String(s) => SymOp::parse_xyz(s).map_err(|er| r.err(it.span, format!("\"{s}\" does not read: {er}"))).ok(),
                            v => {
                                r.err(it.span, format!("an operator is a string \"x,y,z\", not {}", v.kind()));
                                None
                            }
                        })
                        .collect()
                });
                match ops {
                    Some(ops) => match SpaceGroup::from_operators(&ops) {
                        Ok(mut g2) => {
                            if let Some(g) = g {
                                g2.hall = g.hall;
                            }
                            Some(g2)
                        }
                        Err(er) => {
                            r.err(e.value.span, format!("these operators do not make a group: {er}"));
                            None
                        }
                    },
                    None => g,
                }
            }
            (g, None) => g,
        };
        let cell_vals: Vec<f64> = ["cell_a", "cell_b", "cell_c", "cell_alpha", "cell_beta", "cell_gamma"].iter().map(|k| r.f(pt, k, None)).collect();
        let cell = Cell::new(cell_vals[0], cell_vals[1], cell_vals[2], cell_vals[3], cell_vals[4], cell_vals[5]).map_err(|e| r.err(pt.get("cell_a").map_or(pt.header.unwrap_or_default(), |e| e.key_span), format!("the cell does not make sense: {e}"))).ok();
        let shape = r.choice(pt, "profile", SHAPES, Shape::PseudoVoigt);
        let x = match (pt.get("lorentzian_strain"), pt.get("pv_eta_slope")) {
            (Some(_), Some(e)) => {
                r.err(e.key_span, "give one of \"lorentzian_strain\" (Thompson–Cox–Hastings) and \"pv_eta_slope\" (pseudo-Voigt)");
                0.0
            }
            (Some(_), None) => r.f(pt, "lorentzian_strain", None),
            (None, Some(_)) => r.f(pt, "pv_eta_slope", None),
            (None, None) => 0.0,
        };
        let asym = r.nums_key(pt, "berar_baldinozzi", Some(4), Some(vec![0.0; 4]));
        let po_dir = r.nums_key(pt, "po_direction", Some(3), Some(vec![0.0, 0.0, 1.0]));
        // Preferred orientation along several axes (P30): directions, ratios and fractions of as many axes, 1 to 5.
        let texture = pt.get("po_axis_direction").map(|e| {
            let span = e.value.span;
            let dirs = r.table_rows(pt, "po_axis_direction", 3);
            let n = dirs.len();
            let ratio = r.nums_key(pt, "po_axis_ratio", Some(n), None);
            let fraction = r.nums_key(pt, "po_axis_fraction", Some(n), None);
            if !(1..=5).contains(&n) {
                r.err(span, format!("{n} axes: 1 to 5 are taken"));
            }
            if let Some(x) = ratio.iter().find(|&&x| x <= 0.0) {
                r.err(pt.get("po_axis_ratio").map_or(span, |e| e.value.span), format!("a March ratio must be above 0, not {x}"));
            }
            let sum: f64 = fraction.iter().sum();
            if n > 0 && (sum - 1.0).abs() > 1e-4 {
                r.err(pt.get("po_axis_fraction").map_or(span, |e| e.value.span), format!("the fractions sum to {sum}, not 1"));
            }
            let steps = pt.get("po_axis_steps").map_or(0.0, |_| r.f(pt, "po_axis_steps", None));
            if steps < 0.0 || steps.fract() != 0.0 {
                r.err(pt.get("po_axis_steps").map_or(span, |e| e.value.span), "po_axis_steps is a whole number of steps");
            }
            let mut t = crate::texture::MultiAxial {
                axes: dirs.iter().zip(&ratio).zip(&fraction).map(|((d, &r), &f)| crate::texture::Axis { direction: [d[0], d[1], d[2]], r, fraction: f }).collect(),
                steps: steps.max(0.0) as usize,
            };
            t.close();
            t
        });
        for k in ["po_axis_ratio", "po_axis_fraction", "po_axis_steps"] {
            if texture.is_none() {
                if let Some(e) = pt.get(k) {
                    r.err(e.key_span, format!("{k} belongs to a phase with po_axis_direction"));
                }
            }
        }
        let list = pt.get("reflection_intensity").map(|e| {
            r.array(&e.value, "\"reflection_intensity\"")
                .unwrap_or(&[])
                .iter()
                .filter_map(|row| {
                    let v = r.nums(row, "a reflection [h, k, l, multiplicity, intensity]", Some(5))?;
                    if v[..4].iter().any(|x| x.fract() != 0.0) {
                        r.err(row.span, "h, k, l and the multiplicity are whole numbers");
                        return None;
                    }
                    Some(([v[0] as i32, v[1] as i32, v[2] as i32], v[3] as usize, v[4]))
                })
                .collect::<Vec<_>>()
        });
        let (Some(group), Some(cell)) = (group, cell) else { continue };
        // The quartic form of strain (P25): each row a term the group allows, H + K + L = 4.
        let quartic_strain = pt.get("quartic_strain").map(|e| {
            let span = e.value.span;
            let allowed = crate::microstructure::quartic_terms(&group).unwrap_or_else(|m| {
                r.err(span, m);
                Vec::new()
            });
            let mut terms: Vec<([u8; 3], f64)> = Vec::new();
            for v in r.table_rows(pt, "quartic_strain", 4) {
                let ok = v[..3].iter().all(|x| x.fract() == 0.0 && *x >= 0.0) && v[0] + v[1] + v[2] == 4.0;
                let t = [v[0] as u8, v[1] as u8, v[2] as u8];
                if !ok {
                    r.err(span, format!("[{}, {}, {}] is not a term of the quartic form: H, K, L are whole numbers ≥ 0 with H + K + L = 4", v[0], v[1], v[2]));
                } else if !allowed.is_empty() && !allowed.contains(&t) {
                    let names: Vec<String> = allowed.iter().map(|t| crate::microstructure::term_name(*t)).collect();
                    r.err(span, format!("{} is not allowed by the phase's Laue class {} ({})", crate::microstructure::term_name(t), group.laue, names.join(", ")));
                } else if terms.iter().any(|x| x.0 == t) {
                    r.err(span, format!("{} is given twice", crate::microstructure::term_name(t)));
                } else {
                    terms.push((t, v[3]));
                }
            }
            crate::microstructure::QuarticStrain::for_group(&group, terms, r.f(pt, "quartic_strain_lorentzian", Some(0.0)))
        });
        if quartic_strain.is_none() && pt.get("quartic_strain_lorentzian").is_some() {
            r.err(pt.get("quartic_strain_lorentzian").unwrap().key_span, "quartic_strain_lorentzian belongs to a phase with quartic_strain");
        }
        let mut modes = None;
        let mut rigid: Vec<crate::rigid_body::Group> = Vec::new();
        let mut parent = None;
        let intensities = match list {
            Some(list) => {
                if let Some(a) = atoms.first().or(mode_tables.first()) {
                    r.err(a.header.unwrap_or_default(), "a phase with reflection_intensity has no atoms or modes");
                }
                Intensities::List(list)
            }
            None => {
                let mut st = Structure { group: group.clone(), cell, atoms: Vec::new() };
                let mut sc = Vec::new();
                for at in atoms {
                    r.check_keys(at, ATOM_KEYS);
                    let label = r.s(at, "label").map_or_else(
                        || {
                            r.missing(at, "label");
                            String::new()
                        },
                        |s| s.0,
                    );
                    let (species, sspan) = r.s(at, "species").unwrap_or_else(|| {
                        r.missing(at, "species");
                        (String::new(), at.header.unwrap_or_default())
                    });
                    let x = r.nums_key(at, "atom_position", Some(3), None);
                    let beta = at.get("atom_betas").and_then(|e| r.nums(&e.value, "\"atom_betas\"", Some(6))).map(|b| [b[0], b[1], b[2], b[3], b[4], b[5]]);
                    let charge = match at.get("charge").map(|e| &e.value) {
                        Some(Item { value: Value::Integer(q), .. }) if *q != 0 && q.abs() <= 8 => Some(*q as i32),
                        Some(it) => {
                            r.err(it.span, "\"charge\" is the formal oxidation state, a whole number from −8 to 8 other than 0 (3 for Y³⁺, −2 for O²⁻)");
                            None
                        }
                        None => None,
                    };
                    let atom = Atom { label, species: species.clone(), x: [x[0], x[1], x[2]], biso: r.f(at, "atom_biso", Some(0.0)), occ: r.f(at, "atom_occupancy", None), beta, charge };
                    match crate::scattering::resolve(&species, radiation, &user) {
                        Ok(s) => sc.push(s),
                        Err(e) => r.err(sspan, format!("\"{species}\": {e}")),
                    }
                    st.atoms.push(atom);
                }
                // Rigid groups (P41): each table's kind and values, its atoms those that name it, in their order.
                rigid = read_rigid(&mut r, rigid_tables, atoms);
                if !rigid.is_empty() {
                    let mut x: Vec<[f64; 3]> = st.atoms.iter().map(|a| a.x).collect();
                    crate::rigid_body::place(&mut rigid, &st.cell, &mut x);
                    for (a, x) in st.atoms.iter_mut().zip(x) {
                        a.x = x;
                    }
                }
                // The parent: a mode phase's (in its modes), or that of a phase refined by its coordinates (P17).
                if let Some((group, _)) = r.s(pt, "mode_parent_group") {
                    let b = r.table_rows(pt, "mode_parent_basis", 3);
                    let o = if pt.get("mode_parent_origin").is_some() { r.nums_key(pt, "mode_parent_origin", Some(3), None) } else { vec![0.0; 3] };
                    if b.len() == 3 {
                        parent = Some(crate::modes::ModeParent { group, basis: [[b[0][0], b[0][1], b[0][2]], [b[1][0], b[1][1], b[1][2]], [b[2][0], b[2][1], b[2][2]]], origin: [o[0], o[1], o[2]] });
                    } else if let Some(e) = pt.get("mode_parent_group") {
                        r.err(e.value.span, "mode_parent_group needs mode_parent_basis, three rows of three numbers (the columns of P are this phase's basis vectors in the parent's)");
                    }
                }
                // Symmetry modes: the atoms' positions above are the reference; the structure takes those the modes give.
                if !mode_tables.is_empty() {
                    let mut ms = crate::modes::Modes { reference: st.atoms.iter().map(|a| a.x).collect(), modes: Vec::new(), reference_cell: None, max_amplitude: None, fst: Vec::new(), parent: parent.take(), sign_reference: None };
                    // The sign reference (E14b): its file read relative to the control file's folder; the parent it
                    // states taken when the control file states none.
                    if let Some(e) = pt.get("mode_reference") {
                        let span = e.value.span;
                        match &e.value.value {
                            Value::Table(t) => {
                                let file = t.iter().find(|x| x.key == "file").and_then(|x| if let Value::String(v) = &x.value.value { Some(v.clone()) } else { None });
                                let analysis = t.iter().find(|x| x.key == "analysis").and_then(|x| if let Value::Integer(v) = x.value.value { Some(v as usize) } else { None });
                                if let Some(x) = t.iter().find(|x| x.key != "file" && x.key != "analysis") {
                                    r.err(x.key_span, format!("\"{}\" is not a key of mode_reference (file, analysis)", x.key));
                                }
                                match file {
                                    None => r.err(span, "mode_reference needs file = \"…\" (an AMPLIMODES output or FullProf's .sum)"),
                                    Some(f) => {
                                        let path = base.map_or_else(|| std::path::PathBuf::from(&f), |b| b.join(&f));
                                        match std::fs::read(&path) {
                                            Err(err) => r.err(span, format!("{f}: {err}")),
                                            Ok(bytes) => {
                                                let text: String = String::from_utf8(bytes.clone()).unwrap_or_else(|_| bytes.iter().map(|&b| b as char).collect());
                                                match crate::modes::SignReference::read(&f, &text, analysis, model.phases.len()) {
                                                    Ok((sr, pa)) => {
                                                        if ms.parent.is_none() {
                                                            ms.parent = pa;
                                                        }
                                                        ms.sign_reference = Some(sr);
                                                    }
                                                    Err(m) => r.err(span, m),
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            _ => r.err(span, "mode_reference is { file = \"…\", analysis = k }"),
                        }
                    }
                    if pt.get("mode_reference_cell").is_some() {
                        let c = r.nums_key(pt, "mode_reference_cell", Some(6), None);
                        ms.reference_cell = Cell::new(c[0], c[1], c[2], c[3], c[4], c[5]).map_err(|e| r.err(pt.get("mode_reference_cell").unwrap().value.span, format!("the reference cell does not make sense: {e}"))).ok();
                    }
                    if pt.get("mode_max_amplitude").is_some() {
                        ms.max_amplitude = Some(r.f(pt, "mode_max_amplitude", None));
                    }
                    if pt.get("mode_fst").is_some() {
                        ms.fst = r.nums_key(pt, "mode_fst", None, None).iter().map(|&v| v as i32).collect();
                    }
                    for mt in mode_tables {
                        r.check_keys(mt, MODE_KEYS);
                        let name = r.s(mt, "name").map_or_else(|| { r.missing(mt, "name"); String::new() }, |s| s.0);
                        let irrep = r.s(mt, "irrep").map(|s| s.0).or_else(|| crate::modes::irrep_of_name(&name).map(str::to_string)).unwrap_or_default();
                        let amplitude = r.f(mt, "mode_amplitude", Some(0.0));
                        let mut vectors = Vec::new();
                        match mt.get("mode_vector") {
                            Some(e) => {
                                for row in r.array(&e.value, "\"mode_vector\"").unwrap_or(&[]) {
                                    let Value::Array(items) = &row.value else {
                                        r.err(row.span, "a vector is [atom, Vx, Vy, Vz, coefficient]");
                                        continue;
                                    };
                                    let atom = match items.first().map(|i| &i.value) {
                                        Some(Value::String(l)) => index_of_atom(&st, l).or_else(|| {
                                            r.err(items[0].span, format!("the phase has no atom \"{l}\""));
                                            None
                                        }),
                                        _ => {
                                            r.err(row.span, "a vector is [atom, Vx, Vy, Vz, coefficient], the atom by its label");
                                            None
                                        }
                                    };
                                    let v: Vec<f64> = items.iter().skip(1).filter_map(|i| r.num_item(i, "a component")).collect();
                                    if let (Some(atom), 4) = (atom, v.len()) {
                                        vectors.push(crate::modes::ModeVector { atom, v: [v[0], v[1], v[2]], coeff: v[3] });
                                    } else if atom.is_some() {
                                        r.err(row.span, "a vector is [atom, Vx, Vy, Vz, coefficient]");
                                    }
                                }
                            }
                            None => r.missing(mt, "mode_vector"),
                        }
                        ms.modes.push(crate::modes::Mode { name, irrep, amplitude, vectors });
                    }
                    for (a, x) in st.atoms.iter_mut().zip(ms.positions()) {
                        a.x = x;
                    }
                    modes = Some(ms);
                } else if pt.get("mode_reference_cell").is_some() || pt.get("mode_max_amplitude").is_some() {
                    r.err(pt.header.unwrap_or_default(), "mode_reference_cell and mode_max_amplitude belong to a phase with [[phase.mode]] tables");
                }
                Intensities::Structure(Box::new(st), sc)
            }
        };
        let tof = model.tof.as_ref().map(|_| {
            let mut t = crate::tof::Profile::default();
            for (k, id) in crate::tof::PROFILE_FIELDS.iter().enumerate() {
                if (13..=16).contains(&k) && pt.get("tof_alpha0_thermal").is_none() {
                    continue;
                }
                t.set(k, r.f(pt, id, Some(0.0)));
            }
            t
        });
        if tof.is_none() {
            for k in crate::tof::PROFILE_FIELDS {
                if let Some(e) = pt.get(k) {
                    r.err(e.key_span, format!("{k} belongs to a phase of a time-of-flight pattern, one with tof_dtt1"));
                }
            }
        }
        model.phases.push(Phase {
            name,
            group,
            cell,
            intensities,
            scale: r.f(pt, "scale_factor", None),
            bov: r.f(pt, "overall_b", Some(0.0)),
            shape,
            eta0: r.f(pt, "pv_eta0", Some(0.0)),
            u: r.f(pt, "caglioti_u", Some(0.0)),
            v: r.f(pt, "caglioti_v", Some(0.0)),
            w: r.f(pt, "caglioti_w", None),
            x,
            y: r.f(pt, "lorentzian_size", Some(0.0)),
            gauss_size: r.f(pt, "gaussian_size", Some(0.0)),
            asym: [asym[0], asym[1], asym[2], asym[3]],
            s_l: r.f(pt, "axial_divergence_s_l", Some(0.0)),
            d_l: r.f(pt, "axial_divergence_d_l", Some(0.0)),
            pref: [r.f(pt, "po_g1", Some(0.0)), r.f(pt, "po_g2", Some(0.0))],
            pref_direction: [po_dir[0], po_dir[1], po_dir[2]],
            reflection_limit: pt.get("fullprof_reflection_limit").map(|_| r.f(pt, "fullprof_reflection_limit", None)),
            extract: r.b(pt, "le_bail", false),
            atz: r.f(pt, "qpa_atz", Some(0.0)),
            brindley: r.f(pt, "brindley_coefficient", Some(1.0)),
            internal_standard: pt.get("internal_standard_weight_fraction").map(|_| r.f(pt, "internal_standard_weight_fraction", None)),
            valence_sum_range: pt.get(crate::bonds::VALENCE_RANGE_KEY).map(|_| r.f(pt, crate::bonds::VALENCE_RANGE_KEY, None)),
            modes,
            rigid,
            parent,
            uvw2: pt.get("lambda2_uvw").map(|_| {
                let v = r.nums_key(pt, "lambda2_uvw", Some(3), None);
                [v[0], v[1], v[2]]
            }),
            quartic_strain,
            tof,
            shape_size: pt.get("shape_size").map(|_| {
                use crate::microstructure::{ShapeSize, SizeShape};
                let value = r.f(pt, "shape_size", None);
                let kind = r.choice(pt, "shape_size_kind", &[(0, "platelet"), (1, "needle"), (2, "antiphase")], 1);
                let shape = if kind == 2 {
                    SizeShape::Antiphase
                } else {
                    let t = r.nums_key(pt, "shape_size_direction", Some(3), None);
                    let t = [t[0], t[1], t[2]];
                    if kind == 0 { SizeShape::Platelet(t) } else { SizeShape::Needle(t) }
                };
                ShapeSize { shape, value }
            }),
            texture,
        });
    }
    // Restraints.
    for t in doc.tables.iter().filter(|t| t.array && t.path == ["restraint"]) {
        r.check_keys(t, RESTRAINT_KEYS);
        let pe = t.get("phase");
        let p = match pe.map(|e| &e.value) {
            Some(Item { value: Value::Integer(n), span }) => {
                if *n < 1 || *n as usize > model.phases.len() {
                    r.err(*span, format!("there is no phase {n} ({} phases)", model.phases.len()));
                    continue;
                }
                *n as usize - 1
            }
            Some(it) => {
                r.err(it.span, "\"phase\" is the phase's number, from 1");
                continue;
            }
            None => {
                r.missing(t, "phase");
                continue;
            }
        };
        let Intensities::Structure(st, _) = &model.phases[p].intensities else {
            r.err(pe.unwrap().value.span, "a restraint needs a phase with atoms");
            continue;
        };
        let atoms_of = |r: &mut R, e: &Entry, n: usize| -> Option<Vec<usize>> {
            let a = r.array(&e.value, &format!("\"{}\"", e.key))?;
            if a.len() != n {
                r.err(e.value.span, format!("\"{}\" names {n} atoms", e.key));
                return None;
            }
            a.iter()
                .map(|it| match &it.value {
                    Value::String(s) => index_of_atom(st, s).or_else(|| {
                        r.err(it.span, format!("phase {} has no atom \"{s}\"", p + 1));
                        None
                    }),
                    v => {
                        r.err(it.span, format!("an atom is named by its label, not {}", v.kind()));
                        None
                    }
                })
                .collect()
        };
        let kind = match (t.get("distance"), t.get("angle")) {
            (Some(e), None) => {
                let atoms = atoms_of(&mut r, e, 2);
                let op = t.get("operator").map_or(Some(Operator::identity()), |o| operator(&mut r, &o.value));
                match (atoms, op) {
                    (Some(a), Some(op)) => Kind::Distance { atoms: [a[0], a[1]], op },
                    _ => continue,
                }
            }
            (None, Some(e)) => {
                let atoms = atoms_of(&mut r, e, 3);
                let ops = match t.get("operators") {
                    None => Some([Operator::identity(), Operator::identity()]),
                    Some(o) => r.array(&o.value, "\"operators\"").and_then(|a| {
                        if a.len() != 2 {
                            r.err(o.value.span, "an angle has two operators (second and third atoms)");
                            return None;
                        }
                        Some([operator(&mut r, &a[0])?, operator(&mut r, &a[1])?])
                    }),
                };
                match (atoms, ops) {
                    (Some(a), Some(ops)) => Kind::Angle { atoms: [a[0], a[1], a[2]], ops },
                    _ => continue,
                }
            }
            _ => {
                r.err(t.header.unwrap_or_default(), "a restraint holds one of \"distance\" and \"angle\"");
                continue;
            }
        };
        model.restraints.push(Restraint { phase: p, kind, target: r.f(t, "target", None), sigma: r.f(t, "sigma", None) });
    }
    // Convention.
    if let Some(ct) = find(&["convention"]) {
        r.check_keys(ct, CONVENTION_KEYS);
        let d = Convention::DEFAULT;
        model.convention = Convention {
            lorentz_polarisation: r.b(ct, "lorentz_polarisation", d.lorentz_polarisation),
            berar_baldinozzi: r.b(ct, "berar_baldinozzi", d.berar_baldinozzi),
            axial_divergence: r.b(ct, "axial_divergence", d.axial_divergence),
            preferred_orientation: r.b(ct, "preferred_orientation", d.preferred_orientation),
            excluded_points: r.b(ct, "excluded_points", d.excluded_points),
            peak_window: r.b(ct, "peak_window", d.peak_window),
            reflections: r.b(ct, "reflections", d.reflections),
            microstructure: r.b(ct, "microstructure", d.microstructure),
            angle_restraints: r.b(ct, "angle_restraints", d.angle_restraints),
            window_edge: r.b(ct, "window_edge", d.window_edge),
            background_edge: r.b(ct, "background_edge", d.background_edge),
            displacement: r.b(ct, "displacement", d.displacement),
            background_bounds: r.b(ct, "background_bounds", d.background_bounds),
        };
    }
    // The refinement: addresses need the model.
    let mut refinement = Refinement::default();
    if let Some(rt) = find(&["refine"]) {
        r.check_keys(rt, REFINE_KEYS);
        match rt.get("cycles").map(|e| &e.value) {
            Some(Item { value: Value::Integer(n), .. }) if *n >= 0 => refinement.cycles = Some(*n as usize),
            Some(it) => r.err(it.span, "\"cycles\" is a whole number"),
            None => {}
        }
        if rt.get("convergence").is_some() {
            refinement.eps = Some(r.f(rt, "convergence", None));
        }
        let addr = addr_item;
        // The steering (W12), each value checked.
        if let Some(e) = rt.get("lambda") {
            let l = r.f(rt, "lambda", None);
            if l < 0.0 {
                r.err(e.value.span, "λ's start is 0 (Gauss–Newton first) or above");
            }
            refinement.steering.lambda = l.max(0.0);
        }
        refinement.steering.marquardt = r.b(rt, "marquardt", true);
        if let Some(e) = rt.get("damping") {
            match &e.value.value {
                Value::Table(entries) => {
                    for x in entries {
                        match Class::from_word(&x.key) {
                            Some(c) => {
                                if let Some(d) = r.num_item(&x.value, "a damping") {
                                    check_damping(&mut r, x.value.span, d);
                                    refinement.steering.class_damping[c.index()] = d;
                                }
                            }
                            None => {
                                let words: Vec<&str> = Class::ALL.iter().map(|c| c.word()).collect();
                                r.unknown(x.key_span, format!("\"{}\" is not a class (atomic, anisotropic, profile, global)", x.key), nearest(&x.key, &words));
                            }
                        }
                    }
                }
                _ => r.err(e.value.span, "damping is { atomic = …, anisotropic = …, profile = …, global = … }"),
            }
        }
        if let Some(e) = rt.get("divergence") {
            match &e.value.value {
                Value::Table(entries) => {
                    let mut g = Guards::default();
                    for x in entries {
                        match (x.key.as_str(), &x.value.value) {
                            ("rising", Value::Integer(n)) if *n >= 0 => g.rising = *n as usize,
                            ("jump_from", Value::Integer(n)) if *n >= 1 => g.jump_from = *n as usize,
                            ("jump", _) => {
                                if let Some(j) = r.num_item(&x.value, "jump") {
                                    if j < 0.0 {
                                        r.err(x.value.span, "jump is a number of esd's, 0 (off) or above");
                                    }
                                    g.jump = j.max(0.0);
                                }
                            }
                            ("singular", Value::Bool(b)) => g.singular = *b,
                            ("rising" | "jump_from", _) => r.err(x.value.span, format!("{} is a whole number{}", x.key, if x.key == "jump_from" { " from 1" } else { "" })),
                            ("singular", _) => r.err(x.value.span, "singular is true or false"),
                            (k, _) => r.unknown(x.key_span, format!("\"{k}\" is not a divergence guard (rising, jump, jump_from, singular)"), nearest(k, &["rising", "jump", "jump_from", "singular"])),
                        }
                    }
                    refinement.steering.guards = g;
                }
                _ => r.err(e.value.span, "divergence is { rising = 3, jump = 100, jump_from = 3, singular = true }"),
            }
        }
        if let Some(e) = rt.get("steering") {
            for it in r.array(&e.value, "\"steering\"").unwrap_or(&[]).to_vec() {
                let Value::Table(entries) = &it.value else {
                    r.err(it.span, "a steering entry is { value = \"address\", damping = …, max_shift = … }");
                    continue;
                };
                let get = |k: &str| entries.iter().find(|e| e.key == k);
                for x in entries {
                    if !["value", "damping", "max_shift"].contains(&x.key.as_str()) {
                        r.unknown(x.key_span, format!("\"{}\" is not a key of a steering entry", x.key), nearest(&x.key, &["value", "damping", "max_shift"]));
                    }
                }
                let Some(v) = get("value") else {
                    r.err(it.span, "a steering entry needs value = \"address\"");
                    continue;
                };
                let Value::String(name) = &v.value.value else {
                    r.err(v.value.span, "value is an address (or a variable's name)");
                    continue;
                };
                if let Some(d) = get("damping").and_then(|x| r.num_item(&x.value, "damping").map(|d| (d, x.value.span))) {
                    check_damping(&mut r, d.1, d.0);
                    refinement.steering.damping.push((name.clone(), d.0));
                }
                if let Some(x) = get("max_shift") {
                    if let (Some(f), Some(ms)) = (addr(&mut r, &v.value, &model), r.num_item(&x.value, "max_shift")) {
                        if ms <= 0.0 {
                            r.err(x.value.span, "a maximum shift is above 0, in the value's unit");
                        }
                        refinement.steering.max_shift.push((f, ms));
                    }
                }
            }
        }
        if let Some(e) = rt.get("parameters") {
            for row in r.array(&e.value, "\"parameters\"").unwrap_or(&[]).to_vec() {
                let Some(items) = r.array(&row, "a parameter [address, multiplier, …]").map(<[Item]>::to_vec) else { continue };
                if items.is_empty() || items.len() % 2 != 0 {
                    r.err(row.span, "a parameter is written [address, multiplier, address, multiplier, …]");
                    continue;
                }
                let terms: Option<Vec<(Field, f64)>> = items.chunks(2).map(|c| Some((addr(&mut r, &c[0], &model)?, r.num_item(&c[1], "a multiplier")?))).collect();
                if let Some(t) = terms {
                    refinement.parameters.push(t);
                }
            }
        }
        if let Some(e) = rt.get("limits") {
            for it in r.array(&e.value, "\"limits\"").unwrap_or(&[]).to_vec() {
                let Value::Table(entries) = &it.value else {
                    r.err(it.span, "a limit is { value = \"address\", min = …, max = … }");
                    continue;
                };
                let get = |k: &str| entries.iter().find(|e| e.key == k);
                let Some(v) = get("value") else {
                    r.err(it.span, "a limit needs value = \"address\"");
                    continue;
                };
                let Some(f) = addr(&mut r, &v.value, &model) else { continue };
                let lo = get("min").and_then(|e| r.num_item(&e.value, "min"));
                let hi = get("max").and_then(|e| r.num_item(&e.value, "max"));
                if let (Some(l), Some(h)) = (lo, hi) {
                    if l > h {
                        r.err(it.span, format!("the lower limit {l} is above the upper limit {h}"));
                    }
                }
                if lo.is_some() || hi.is_some() {
                    refinement.limits.push((f, lo, hi));
                }
                if let Some(b) = get("behaviour") {
                    match &b.value.value {
                        Value::String(w) => match Behaviour::from_word(w) {
                            Some(Behaviour::Stop) => {}
                            Some(bh) => refinement.behaviour.push((f, bh)),
                            None => r.unknown(b.value.span, format!("\"{w}\" is not a behaviour at a bound (stop, reflect, periodic)"), nearest(w, &["stop", "reflect", "periodic"])),
                        },
                        _ => r.err(b.value.span, "behaviour is \"stop\", \"reflect\" or \"periodic\""),
                    }
                }
            }
        }
    }
    for t in doc.tables.iter().filter(|t| t.array && t.path == ["refine", "linear"]) {
        r.check_keys(t, LINEAR_KEYS);
        let name = r.s(t, "name").map_or(String::new(), |s| s.0);
        let mut terms = Vec::new();
        if let Some(e) = t.get("terms") {
            for row in r.array(&e.value, "\"terms\"").unwrap_or(&[]).to_vec() {
                let Some(items) = r.array(&row, "a term [parameter, coefficient]").map(<[Item]>::to_vec) else { continue };
                match (items.first().map(|i| &i.value), items.get(1)) {
                    (Some(Value::String(n)), Some(c)) if items.len() == 2 => {
                        if let Some(c) = r.num_item(c, "a coefficient") {
                            terms.push((n.clone(), c));
                        }
                    }
                    _ => r.err(row.span, "a term is [\"parameter\", coefficient]"),
                }
            }
        }
        refinement.linear.push(LinearRestraint { name, terms, target: r.f(t, "target", None), sigma: r.f(t, "sigma", None) });
    }
    for t in doc.tables.iter().filter(|t| t.array && t.path == ["refine", "link"]) {
        r.check_keys(t, LINK_KEYS);
        let mut get = |k: &str| match t.get(k) {
            Some(e) => addr_item(&mut r, &e.value, &model),
            None => {
                r.missing(t, k);
                None
            }
        };
        let (value, follows) = (get("value"), get("follows"));
        let coefficient = r.f(t, "coefficient", Some(1.0));
        let offset = r.f(t, "offset", Some(0.0));
        if let (Some(value), Some(follows)) = (value, follows) {
            refinement.links.push(Link { value, follows, coefficient, offset });
        }
    }
    for t in doc.tables.iter().filter(|t| t.array && t.path == ["refine", "variable"]) {
        r.check_keys(t, VARIABLE_KEYS);
        let Some((name, span)) = r.s(t, "name") else {
            r.missing(t, "name");
            continue;
        };
        if name.is_empty() || name.contains(char::is_whitespace) {
            r.err(span, "a variable's name is one word");
        }
        let value = r.f(t, "value", Some(0.0));
        let mut terms = Vec::new();
        if let Some(e) = t.get("terms") {
            for row in r.array(&e.value, "\"terms\"").unwrap_or(&[]).to_vec() {
                let Some(items) = r.array(&row, "a term [address, coefficient]").map(<[Item]>::to_vec) else { continue };
                if items.len() != 2 {
                    r.err(row.span, "a term is [\"address\", coefficient]");
                    continue;
                }
                if let (Some(f), Some(c)) = (addr_item(&mut r, &items[0], &model), r.num_item(&items[1], "a coefficient")) {
                    terms.push((f, c));
                }
            }
        } else {
            r.missing(t, "terms");
        }
        refinement.variables.push(Variable { name, value, terms });
    }
    for t in doc.tables.iter().filter(|t| t.array && t.path == ["refine", "step"]) {
        r.check_keys(t, STEP_KEYS);
        let name = r.s(t, "name").map_or(format!("step {}", refinement.plan.len() + 1), |s| s.0);
        let mut free = Vec::new();
        match t.get("free") {
            Some(e) => {
                for it in r.array(&e.value, "\"free\"").unwrap_or(&[]).to_vec() {
                    match &it.value {
                        Value::String(x) => free.push(x.clone()),
                        _ => r.err(it.span, "free lists addresses, prefixes ending in *, variables, classes or \"all\""),
                    }
                }
            }
            None => r.missing(t, "free"),
        }
        let damping = t.get("damping").map(|e| {
            let d = r.f(t, "damping", None);
            check_damping(&mut r, e.value.span, d);
            d
        });
        let cycles = match t.get("cycles").map(|e| &e.value) {
            Some(Item { value: Value::Integer(n), .. }) if *n >= 1 => Some(*n as usize),
            Some(it) => {
                r.err(it.span, "\"cycles\" is a whole number from 1");
                None
            }
            None => None,
        };
        let eps = t.get("convergence").map(|_| r.f(t, "convergence", None));
        refinement.plan.push(PlanStep { name, free, damping, cycles, eps });
    }
    // The parameter set must build (symmetry, names of linear restraints, links, variables, damping, bounds), and
    // each step of the plan must free something.
    if (!refinement.parameters.is_empty() || !refinement.linear.is_empty() || !refinement.variables.is_empty() || !refinement.links.is_empty()) && !points_unread {
        let span = find(&["refine"]).and_then(|t| t.get("parameters")).map_or(Span::default(), |e| e.key_span);
        if let Err(e) = refinement.build(&model) {
            r.err(span, e);
        }
        let steps: Vec<&Table> = doc.tables.iter().filter(|t| t.array && t.path == ["refine", "step"]).collect();
        for k in 0..refinement.plan.len() {
            if let Err(e) = refinement.for_step(k, &model).and_then(|s| s.build(&model).map(|_| ())) {
                r.err(steps.get(k).and_then(|t| t.header).unwrap_or(span), e);
            }
        }
    }
    let data = data_file.map(|file| super::DataSource { file, reader: reader.map(|r| r.0), reader_options, operations });
    let errors = r.diags.iter().any(|d| !d.warning);
    let control = (!errors).then_some(Control { title, model, data, refinement, fullprof, background_file });
    let mut diags = r.diags;
    // In the order of the text, as an editor lists them.
    diags.sort_by_key(|d| d.span.start);
    (control, diags)
}

/// Parse and read a control file's text.
pub fn read_text(text: &str) -> (Option<Control>, Vec<Diagnostic>) {
    read(&syntax::parse(text))
}
