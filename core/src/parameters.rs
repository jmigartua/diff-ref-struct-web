//! The core's model fields and the parameter dictionary (L044, L046).
//!
//! `ecosystem/data/parameters.yml` is the single source of the parameters' names and theory; `parameter_table`
//! is generated from it. `FIELDS` names, for every field of the model's structures (`Model`, `Phase`, `Atom`), the
//! dictionary ids it holds, or why it holds none. The test `parameter_dictionary` checks that every field of those
//! structures is listed here, that every id exists in the dictionary, and that the dictionary's own `core` lists say the
//! same; so a field cannot be added to the model without its entry in the dictionary.

pub use crate::parameter_table::{Parameter, GROUPS, PARAMETERS};

/// The whole dictionary: every entry with every field of `parameters.yml` (the app's request 4, L058).
pub fn dictionary() -> &'static [Parameter] {
    PARAMETERS
}

/// The dictionary's groups, (key, title), in order.
pub fn groups() -> &'static [(&'static str, &'static str)] {
    GROUPS
}

/// A field of the model, the dictionary ids it holds, and, when it holds none, why.
pub struct Field {
    pub field: &'static str,
    pub ids: &'static [&'static str],
    pub not_a_parameter: &'static str,
}

const fn p(field: &'static str, ids: &'static [&'static str]) -> Field {
    Field { field, ids, not_a_parameter: "" }
}

const fn not(field: &'static str, why: &'static str) -> Field {
    Field { field, ids: &[], not_a_parameter: why }
}

pub const FIELDS: &[Field] = &[
    not("Model.convention", "the rules followed where FullProf departs from its manual (powder::Convention)"),
    not("Model.neutrons", "the radiation, a choice"),
    p("Model.lambda", &["wavelength_1", "wavelength_2"]),
    p("Model.ratio", &["wavelength_ratio"]),
    p("Model.zero", &["zero_shift"]),
    p("Model.sycos", &["sample_displacement"]),
    p("Model.sysin", &["sample_transparency"]),
    p("Model.cthm", &["monochromator_polarisation"]),
    p("Model.rpolarz", &["polarisation_fraction"]),
    not("Model.ilo", "the geometry, a choice (FullProf's Ilo)"),
    p("Model.mur", &["absorption_mur"]),
    p("Model.microabsorption", &["microabsorption_p0", "microabsorption_cp", "microabsorption_tau"]),
    p("Model.tof", &["tof_dtt1", "tof_dtt2", "tof_dtt_1overd", "tof_thermal_zero", "tof_thermal_dtt1", "tof_thermal_dtt2", "tof_crossover", "tof_crossover_width", "tof_bank_two_theta", "tof_max", "instrument_resolution"]),
    p("Model.wdt", &["peak_window"]),
    p("Model.asym_limit", &["asymmetry_limit"]),
    not("Model.march", "the preferred-orientation function, a choice (FullProf's Nor)"),
    p("Model.background", &["background_points", "background_polynomial", "background_origin", "background_chebyshev", "background_fourier"]),
    p("Model.background_terms", &["background_debye_amplitude", "background_debye_distance", "background_inverse_power", "background_profile", "background_file_scale", "background_file_offset", "background_file_stretch", "background_file_shift"]),
    p("Model.instrumental", &["instrument_resolution"]),
    not("Model.phases", "the list of phases"),
    p("Model.excluded", &["excluded_regions"]),
    not("Model.weighting", "the weights of the agreement factors, a choice (FullProf's Iwg)"),
    not("Model.restraints", "distance and angle restraints: observations added to the fit, not parameters (P40)"),
    not("Model.fullprof_step", "FullProf's printed step, used by its rule for the background's low edge in the FullProf tests (L098)"),
    not("Model.background_ceiling", "the data's largest count, above which FullProf's background is 0 (Convention::background_bounds, L353): read from the data, not refined"),
    not("Phase.name", "a label"),
    not("Phase.group", "the space group, a choice"),
    p("Phase.cell", &["cell_a", "cell_b", "cell_c", "cell_alpha", "cell_beta", "cell_gamma"]),
    p("Phase.intensities", &["reflection_intensity"]),
    p("Phase.scale", &["scale_factor"]),
    p("Phase.bov", &["overall_b"]),
    not("Phase.shape", "the profile function, a choice (FullProf's Npr)"),
    p("Phase.eta0", &["pv_eta0"]),
    p("Phase.u", &["caglioti_u"]),
    p("Phase.v", &["caglioti_v"]),
    p("Phase.w", &["caglioti_w"]),
    p("Phase.x", &["lorentzian_strain", "pv_eta_slope"]),
    p("Phase.y", &["lorentzian_size"]),
    p("Phase.gauss_size", &["gaussian_size"]),
    p("Phase.quartic_strain", &["quartic_strain", "quartic_strain_lorentzian"]),
    p("Phase.shape_size", &["shape_size"]),
    p(
        "Phase.tof",
        &[
            "tof_sigma0", "tof_sigma1", "tof_sigma2", "tof_sigma_q", "tof_gamma0", "tof_gamma1", "tof_gamma2", "tof_alpha0", "tof_alpha1", "tof_alpha_q", "tof_beta0", "tof_beta1", "tof_beta_q",
            "tof_alpha0_thermal", "tof_beta0_thermal", "tof_alpha1_thermal", "tof_beta1_thermal", "tof_extinction", "tof_absorption",
        ],
    ),
    p("Phase.asym", &["berar_baldinozzi"]),
    p("Phase.s_l", &["axial_divergence_s_l"]),
    p("Phase.d_l", &["axial_divergence_d_l"]),
    p("Phase.pref", &["po_g1", "po_g2"]),
    p("Phase.pref_direction", &["po_direction"]),
    p("Phase.texture", &["po_axis_ratio", "po_axis_fraction", "po_axis_steps"]),
    not("Phase.reflection_limit", "FullProf's printed S(1/d)max, used by the FullProf tests only (L035)"),
    p("Phase.modes", &["mode_amplitude", "mode_vector"]),
    p("Phase.rigid", &["rigid_origin", "rigid_orientation", "rigid_length", "rigid_angle", "rigid_tls_t", "rigid_tls_l", "rigid_tls_s"]),
    not("Phase.parent", "the parent structure's group and transformation of a phase refined by its coordinates (the control file's mode_parent_*; P17): what the cell is reduced to, not a refined value"),
    not("Phase.uvw2", "U, V, W of the second wavelength's peaks, read from and written to FullProf's line 32 (Ratio < 0, manual p. 75) and held fixed: the Caglioti U, V, W of λ1 are the refined ones (L098)"),
    not("Phase.extract", "whether the phase's intensities are extracted by Le Bail's iteration, a choice (FullProf's Jbt = 2)"),
    p("Phase.atz", &["qpa_atz"]),
    p("Phase.brindley", &["brindley_coefficient"]),
    not("Phase.valence_sum_range", "the relative range about |q| within which a valence sum is drawn as plausible (Explore's E18c), a reading of the result, not a value of the model"),
    not("Phase.internal_standard", "the weight fraction of the phase weighed into the sample as an internal standard, a measured input that is not refined (P33; Orlhac et al. 2001)"),
    not("Atom.label", "a label"),
    not("Atom.species", "the scattering species, a choice"),
    not("Atom.charge", "the formal charge for bond-valence sums (P5), a choice"),
    p("Atom.x", &["atom_position"]),
    p("Atom.biso", &["atom_biso"]),
    p("Atom.occ", &["atom_occupancy"]),
    p("Atom.beta", &["atom_betas"]),
];

/// The dictionary entry of an id.
pub fn lookup(id: &str) -> Option<&'static Parameter> {
    PARAMETERS.iter().find(|p| p.id == id)
}
