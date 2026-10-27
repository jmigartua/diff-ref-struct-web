//! The parameter a line of an output names (the app's request 36, P49's jump from an output to its parameter):
//! FullProf's symbolic names, as its `.out` and `.sum` print them ("Parameter number n -> Symbolic Name:", manual
//! p. 9), and the standard CIF tags of a results file or a structure CIF (`_cell_length_a`, an `_atom_site_` row's
//! label and column), each as the model's field and its address (`phase1.cell_a`).
//!
//! FullProf's grammar: a name, then `_ph<n>` for a phase's parameter and `_pat<n>` for a pattern's (one pattern
//! here); an atom's parameters are `<kind>_<label>` (X_Y1, Biso_O1, Occ_Y2, bet11_O …), a symmetry mode's amplitude
//! `A<k>_<irrep>`, the micro-absorption's `P0_mabs`, `Cp_mabs`, `Tau_mabs`. Checked on every `.out` of
//! `manuals/` whose `.pcr` the core calculates: each name maps to one of the fields its parameter number carries in
//! the model built from that `.pcr` (test `output_names`).

use crate::jacobian::Field;
use crate::powder::{Intensities, Model};

fn atom_index(model: &Model, phase: usize, label: &str) -> Option<usize> {
    match &model.phases.get(phase)?.intensities {
        Intensities::Structure(st, _) => st.atoms.iter().position(|a| a.label.eq_ignore_ascii_case(label)),
        Intensities::List(_) => None,
    }
}

/// The field FullProf's symbolic name names in `model` (the model built from the same `.pcr`); None for a name the
/// core does not refine (magnetic, size and strain models, a second pattern) or that names no atom of the model.
pub fn fullprof_field(name: &str, model: &Model) -> Option<Field> {
    let n = name.trim();
    let n = match n.rfind("_pat") {
        Some(i) if n[i + 4..].parse::<usize>().ok()? == 1 => &n[..i],
        Some(_) => return None,
        None => n,
    };
    let (base, phase) = match n.rfind("_ph") {
        Some(i) => (&n[..i], n[i + 3..].parse::<usize>().ok()?.checked_sub(1)?),
        None => (n, 0),
    };
    if phase >= model.phases.len().max(1) {
        return None;
    }
    Some(match base {
        "Zero" => Field::Zero,
        "SyCos" | "Sycos" => Field::Sycos,
        "SySin" | "Sysin" => Field::Sysin,
        "Lambda" => Field::Wavelength,
        "Lambda_rat" => Field::Ratio,
        "Scale" => Field::Scale(phase),
        "Bover" => Field::Bov(phase),
        "U-Cagl" => Field::U(phase),
        "V-Cagl" => Field::V(phase),
        "W-Cagl" => Field::W(phase),
        "X-tan" => Field::X(phase),
        "Y-cos" => Field::Y(phase),
        "G-Size" => Field::GaussSize(phase),
        "EtaPV" => Field::Eta0(phase),
        "SL_asym" => Field::SL(phase),
        "DL_asym" => Field::DL(phase),
        // FullProf 8.20 names Pref1 Or1 in rigid/'s pdd runs (its value is the .pcr's Pref1, 1.58325).
        "Pref1" | "Or1" => Field::Pref(phase, 0),
        "Pref2" => Field::Pref(phase, 1),
        b if b.starts_with("Asym") => Field::Asym(phase, b[4..].parse::<usize>().ok()?.checked_sub(1)?),
        b if b.starts_with("Cell_") => Field::Cell(phase, "ABCDEF".find(&b[5..])?),
        // Nba = −1 (P23): FullProf's Bck_12 … Bck_17 are the Debye-like terms' amplitudes (runs/background_terms/
        // nac_3t2.out: Bck_12_pat1 −37.464, BC1 of its .pcr) and Bck_18 … Bck_23 their distances (believed: the runs
        // refine no distance).
        b if b.starts_with("Bck_") => {
            let k: usize = b[4..].parse().ok()?;
            let nd = model.background_terms.debye.len();
            match k.checked_sub(12) {
                Some(j) if j < nd => Field::BackgroundTerm(2 * j),
                Some(j) if (6..6 + nd).contains(&j) => Field::BackgroundTerm(2 * (j - 6) + 1),
                _ => Field::Background(k),
            }
        }
        // The quartic form of strain (P25): Strain<n> is the n-th S_HKL in FullProf's order for the Strain-Model, LorStr
        // its Lorentzian share ξ; Ani-LSiz the anisotropic size (L267).
        b if b.starts_with("Strain") && b[6..].parse::<usize>().is_ok() => {
            let q = model.phases[phase].quartic_strain.as_ref()?;
            let terms: Vec<[u8; 3]> = q.terms.iter().map(|t| t.0).collect();
            // FullProf's strain model 8 of a trigonal or hexagonal phase: Str1–Str3 (Saa, Scc, Cac) are refined as
            // S_400, S_004 and S_112 (`powder::hexagonal_strain`, L322).
            let hexagonal: &[[u8; 3]] = &[[4, 0, 0], [0, 0, 4], [1, 1, 2]];
            let order = match crate::microstructure::fullprof_model_of(&terms) {
                Some(m) => crate::microstructure::fullprof_strain_model(m)?,
                None if matches!(model.phases[phase].group.laue, "-3" | "-3m" | "6/m" | "6/mmm") => hexagonal,
                None => return None,
            };
            let t = order.get(b[6..].parse::<usize>().ok()?.checked_sub(1)?)?;
            Field::Quartic(phase, terms.iter().position(|x| x == t)?)
        }
        "LorStr" if model.phases[phase].quartic_strain.is_some() => Field::QuarticShare(phase),
        "Ani-LSiz" if model.phases[phase].shape_size.is_some() => Field::ShapeSize(phase),
        // A TOF pattern's names (P05; as the .out of L334's runs print them): Dtt1, Dtt2, Dtt_1overD of the bank, and
        // per phase Sig-n, Gam-n, Alpha0/1/Q, Beta0/1/Q, ExtTOF (Extinc) and TOF_Abs1.
        "Dtt1" if model.tof.is_some() => Field::Tof(0),
        "Dtt2" if model.tof.is_some() => Field::Tof(1),
        // With Hodges' terms (Npr = 10) FullProf names Zerot so: gabriela.out's parameter 9, Dtt_1overD_pat1 = 23.96,
        // is line 16's Zerot (codeword 91), the line 14 of gabriela.pcr having no Dtt_1overd.
        "Dtt_1overD" if model.tof.as_ref().is_some_and(|b| b.thermal.is_some()) => Field::Tof(3),
        "Dtt_1overD" if model.tof.is_some() => Field::Tof(2),
        b if model.phases[phase].tof.is_some() && tof_profile_name(b).is_some() => Field::TofProfile(phase, tof_profile_name(b)?),
        // Bragg–Brentano surface micro-absorption (P29).
        "P0_mabs" => Field::Micro(0),
        "Cp_mabs" => Field::Micro(1),
        "Tau_mabs" => Field::Micro(2),
        // A symmetry mode's amplitude: A<k>_<irrep>, the k-th mode of the phase (Jbt = 6).
        b if b.starts_with('A') && b[1..].split('_').next().is_some_and(|k| k.parse::<usize>().is_ok()) && model.phases[phase].modes.is_some() => {
            let k: usize = b[1..].split('_').next()?.parse().ok()?;
            let count = model.phases[phase].modes.as_ref()?.modes.len();
            (1..=count).contains(&k).then_some(Field::Mode(phase, k - 1))?
        }
        // A TLS value (L367): t11 … s32, or T11 … (pdd's runs), of the phase's TLS group.
        b if crate::rigid_body::TLS_NAMES.iter().any(|n| n.eq_ignore_ascii_case(b)) => {
            let i = crate::rigid_body::TLS_NAMES.iter().position(|n| n.eq_ignore_ascii_case(b))?;
            let g = model.phases[phase].rigid.iter().position(|g| g.displacements.is_some())?;
            Field::Rigid(phase, g, 6 + 3 * model.phases[phase].rigid[g].members.len() + i)
        }
        b => {
            // A group's origin by its first atom (X_RGB_u1) and an atom's internal coordinates (Dm_In_u2, TheIn_u3,
            // PhiIn_…), as FullProf 8.20 names them in rigid/'s runs (L367).
            let b = b.replacen("_RGB_", "RGB_", 1).replacen("Dm_In_", "DmIn_", 1);
            let (kind, label) = b.split_once('_')?;
            let a = atom_index(model, phase, label)?;
            let internal = |c: usize| model.phases[phase].rigid.iter().enumerate().find_map(|(g, gr)| gr.members.iter().position(|m| m.index == a).map(|m| Field::Rigid(phase, g, 6 + 3 * m + c)));
            // A rigid group's θ, φ, χ by its first atom (P41: ThRGB_A1, citric-acid-sls_1.out, PPH3-sls.out).
            let rigid = |k: usize| model.phases[phase].rigid.iter().position(|g| g.members[0].index == a).map(|g| Field::Rigid(phase, g, k));
            match kind {
                "ThRGB" => rigid(3)?,
                "PhRGB" => rigid(4)?,
                "ChRGB" => rigid(5)?,
                "XRGB" => rigid(0)?,
                "YRGB" => rigid(1)?,
                "ZRGB" => rigid(2)?,
                "DmIn" => internal(0)?,
                "TheIn" => internal(1)?,
                "PhiIn" => internal(2)?,
                "X" => Field::AtomX(phase, a, 0),
                "Y" => Field::AtomX(phase, a, 1),
                "Z" => Field::AtomX(phase, a, 2),
                "Biso" => Field::AtomB(phase, a),
                "Occ" => Field::AtomOcc(phase, a),
                k if k.starts_with("bet") => Field::AtomBeta(phase, a, ["11", "22", "33", "12", "13", "23"].iter().position(|s| *s == &k[3..])?),
                _ => return None,
            }
        }
    })
}

/// A TOF phase's value FullProf's name names (`tof::PROFILE_FIELDS`' order).
fn tof_profile_name(b: &str) -> Option<usize> {
    Some(match b {
        "Sig-0" => 0,
        "Sig-1" => 1,
        "Sig-2" => 2,
        "Gam-0" => 4,
        "Gam-1" => 5,
        "Gam-2" => 6,
        "Alpha0" => 7,
        "Alpha1" => 8,
        "AlphaQ" => 9,
        "Beta0" => 10,
        "Beta1" => 11,
        "BetaQ" => 12,
        "ExtTOF" => 17,
        "TOF_Abs1" => 18,
        _ => return None,
    })
}

/// Every field FullProf's symbolic name can name in `model`: one, except where atoms share a label (FullProf names an
/// atom's parameters by its label, so `Biso_MN_ph1` names both B's of Pd₃MnD's two MN atoms, pdfinal.out's parameters
/// 4 and 6, L267); `fullprof_field` gives the first, and the parameter number tells them apart.
pub fn fullprof_fields(name: &str, model: &Model) -> Vec<Field> {
    let Some(f) = fullprof_field(name, model) else { return Vec::new() };
    let (phase, a) = match f {
        Field::AtomX(p, a, _) | Field::AtomB(p, a) | Field::AtomOcc(p, a) | Field::AtomBeta(p, a, _) => (p, a),
        _ => return vec![f],
    };
    let Some(Intensities::Structure(st, _)) = model.phases.get(phase).map(|p| &p.intensities) else { return vec![f] };
    let label = &st.atoms[a].label;
    (0..st.atoms.len())
        .filter(|&b| st.atoms[b].label.eq_ignore_ascii_case(label))
        .map(|b| match f {
            Field::AtomX(p, _, k) => Field::AtomX(p, b, k),
            Field::AtomB(p, _) => Field::AtomB(p, b),
            Field::AtomOcc(p, _) => Field::AtomOcc(p, b),
            Field::AtomBeta(p, _, k) => Field::AtomBeta(p, b, k),
            other => other,
        })
        .collect()
}

/// FullProf's symbolic name as the model's address (`refinable::address`).
pub fn address_of_fullprof(name: &str, model: &Model) -> Option<String> {
    fullprof_field(name, model).map(|f| crate::refinable::address(f, model))
}

/// The field a CIF tag names in `model`: `phase` is the index of the phase whose block holds the tag, `label` the
/// `_atom_site_label` (or `_atom_site_aniso_label`) of its row for the atom loops. The tags the core's CIFs write with
/// a standard uncertainty, and their B and β variants: `_cell_length_a` … `_cell_angle_gamma`, `_atom_site_fract_x`
/// … `_z`, `_atom_site_U_iso_or_equiv` and `_B_iso_or_equiv`, `_atom_site_occupancy`, `_atom_site_aniso_U_11` …
/// (and `_B_`, `_beta_`), `_diffrn_radiation_wavelength`. Case is ignored, as in CIF.
pub fn cif_field(tag: &str, phase: usize, label: Option<&str>, model: &Model) -> Option<Field> {
    let t = tag.trim().to_ascii_lowercase();
    let t = t.trim_start_matches('_');
    if t == "diffrn_radiation_wavelength" {
        return Some(Field::Wavelength);
    }
    if phase >= model.phases.len() {
        return None;
    }
    let cell = ["cell_length_a", "cell_length_b", "cell_length_c", "cell_angle_alpha", "cell_angle_beta", "cell_angle_gamma"];
    if let Some(k) = cell.iter().position(|c| *c == t) {
        return Some(Field::Cell(phase, k));
    }
    let a = atom_index(model, phase, label?)?;
    Some(match t {
        "atom_site_fract_x" => Field::AtomX(phase, a, 0),
        "atom_site_fract_y" => Field::AtomX(phase, a, 1),
        "atom_site_fract_z" => Field::AtomX(phase, a, 2),
        "atom_site_u_iso_or_equiv" | "atom_site_b_iso_or_equiv" => Field::AtomB(phase, a),
        "atom_site_occupancy" => Field::AtomOcc(phase, a),
        _ => {
            let ij = t.strip_prefix("atom_site_aniso_u_").or_else(|| t.strip_prefix("atom_site_aniso_b_")).or_else(|| t.strip_prefix("atom_site_aniso_beta_"))?;
            Field::AtomBeta(phase, a, ["11", "22", "33", "12", "13", "23"].iter().position(|s| *s == ij)?)
        }
    })
}

/// A CIF tag as the model's address.
pub fn address_of_cif(tag: &str, phase: usize, label: Option<&str>, model: &Model) -> Option<String> {
    cif_field(tag, phase, label, model).map(|f| crate::refinable::address(f, model))
}
