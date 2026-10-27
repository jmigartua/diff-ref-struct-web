//! Every editable value of a model, enumerated, and the model's structural edits (the app's request 3, L058).
//!
//! `Model::fields` lists each field with its address (`refinable::address`), dictionary id, value, limits and whether
//! the crystal system or the site symmetry fixes or ties it; `refinable::field_of` turns an address back into a field;
//! `Model::set_field` sets a value with those checks (the engine's `Model::set` stays unchecked). The structural edits
//! (phases, atoms, background, excluded regions, profile, asymmetry and preferred orientation) each return a `Result`
//! and leave the model unchanged when they fail.

use crate::cell::Cell;
use crate::jacobian::Field;
use crate::powder::{BackgroundModel, Intensities, Model, Phase, Shape};
use crate::refinable::{address, cell_rules, dictionary_id, limits, CellRule};
use crate::restraints::Kind;
use crate::scattering::Radiation;
use crate::site::Site;
use crate::spacegroup::SpaceGroup;
use crate::structure_factor::{Atom, Structure};

/// One editable value of a model.
#[derive(Clone, Debug, PartialEq)]
pub struct FieldInfo {
    pub field: Field,
    /// `pattern1.zero_shift`, `phase1.O1.atom_position.x`, …
    pub address: String,
    /// The parameter dictionary's id.
    pub id: &'static str,
    pub value: f64,
    /// The bounds the engine holds it within (P39).
    pub limits: (Option<f64>, Option<f64>),
    /// "crystal system" or "site symmetry" when the value follows from others (fixed, or tied to another); None when
    /// it can be set by itself.
    pub fixed_by: Option<&'static str>,
}

/// The asymmetry function of a phase.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Asymmetry {
    None,
    /// Bérar & Baldinozzi's P1–P4 (the paper's F_b, L050).
    BerarBaldinozzi([f64; 4]),
    /// Finger, Cox & Jephcoat's axial divergence S/L and H/L (TCH profile).
    AxialDivergence { s_l: f64, h_l: f64 },
}

fn background_len(b: &BackgroundModel) -> usize {
    match b {
        BackgroundModel::Points(p) | BackgroundModel::Spline(p) => p.len(),
        BackgroundModel::Polynomial { coefficients, .. } | BackgroundModel::Chebyshev { coefficients, .. } | BackgroundModel::Fourier { coefficients } => coefficients.len(),
    }
}

impl Model {
    /// Every field of the model that a value can be given to, in a fixed order: the pattern's, then each phase's, then
    /// each atom's.
    pub fn all_fields(&self) -> Vec<Field> {
        let mut f = vec![Field::Zero];
        // A TOF bank (P05): its TOF(d) in place of the displacement, transparency and wavelengths.
        if let Some(b) = &self.tof {
            f.extend((0..if b.thermal.is_some() { 8 } else { 3 }).map(Field::Tof));
        } else {
            f.extend([Field::Sycos, Field::Sysin]);
            if (self.lambda[1] - self.lambda[0]).abs() > 1e-9 {
                f.push(Field::Ratio);
            }
            f.push(Field::Wavelength);
        }
        // Surface micro-absorption, for X-rays in Bragg–Brentano geometry (manual p. 84).
        if !self.neutrons && self.ilo == 0 {
            f.extend((0..3).map(Field::Micro));
        }
        f.extend((0..background_len(&self.background)).map(Field::Background));
        f.extend((0..self.background_terms.len()).map(Field::BackgroundTerm));
        for (p, ph) in self.phases.iter().enumerate() {
            f.extend([Field::Scale(p), Field::Bov(p)]);
            f.extend((0..6).map(|k| Field::Cell(p, k)));
            // A TOF phase's profile (P05): its σ, γ, α, β (and Hodges' thermal ones), extinction and absorption, in
            // place of U … Y and the asymmetry.
            if let Some(t) = &ph.tof {
                let n = if t.thermal.is_some() { 17 } else { 13 };
                f.extend((0..n).chain([17, 18]).map(|k| Field::TofProfile(p, k)));
                f.extend(pref_fields(ph, p));
                if let Intensities::Structure(st, _) = &ph.intensities {
                    for (a, at) in st.atoms.iter().enumerate() {
                        f.extend((0..3).map(|k| Field::AtomX(p, a, k)));
                        f.extend([Field::AtomB(p, a), Field::AtomOcc(p, a)]);
                        if at.beta.is_some() {
                            f.extend((0..6).map(|k| Field::AtomBeta(p, a, k)));
                        }
                    }
                }
                continue;
            }
            f.extend([Field::U(p), Field::V(p), Field::W(p), Field::X(p)]);
            if ph.shape == Shape::PseudoVoigt {
                f.push(Field::Eta0(p));
            } else {
                f.extend([Field::Y(p), Field::GaussSize(p)]);
            }
            f.extend((0..4).map(|k| Field::Asym(p, k)));
            if ph.shape == Shape::ThompsonCoxHastings {
                f.extend([Field::SL(p), Field::DL(p)]);
            }
            f.extend(pref_fields(ph, p));
            if let Some(q) = &ph.quartic_strain {
                f.extend((0..q.terms.len()).map(|k| Field::Quartic(p, k)));
                f.push(Field::QuarticShare(p));
            }
            if ph.shape_size.is_some() {
                f.push(Field::ShapeSize(p));
            }
            if let Intensities::Structure(st, _) = &ph.intensities {
                for (a, at) in st.atoms.iter().enumerate() {
                    f.extend((0..3).map(|k| Field::AtomX(p, a, k)));
                    f.extend([Field::AtomB(p, a), Field::AtomOcc(p, a)]);
                    if at.beta.is_some() {
                        f.extend((0..6).map(|k| Field::AtomBeta(p, a, k)));
                    }
                }
            }
            if let Some(ms) = &ph.modes {
                f.extend((0..ms.modes.len()).map(|m| Field::Mode(p, m)));
            }
            for (g, gr) in ph.rigid.iter().enumerate() {
                f.extend((0..crate::rigid_body::values_of(gr)).map(|k| Field::Rigid(p, g, k)));
            }
            // A list's intensities, but a Le Bail phase's, which are extracted, not refined (P36).
            if let (Intensities::List(l), false) = (&ph.intensities, ph.extract) {
                f.extend((0..l.len()).map(|k| Field::Intensity(p, k)));
            }
        }
        f
    }

    /// What fixes or ties a field: the crystal system for the cell, the site symmetry for an atom's coordinates and
    /// betas.
    pub fn fixed_by(&self, f: Field) -> Option<&'static str> {
        match f {
            Field::Cell(p, k) => match cell_rules(&self.phases[p])[k] {
                CellRule::Free => None,
                _ => Some("crystal system"),
            },
            Field::AtomX(p, ..) if self.phases[p].modes.is_some() => Some("symmetry modes"),
            Field::AtomX(p, a, _) if self.phases[p].rigid.iter().any(|g| g.members.iter().any(|m| m.index == a)) => Some("rigid group"),
            // A satellite's origin, θ and φ follow its parent (its χ is its own).
            Field::Rigid(p, g, 0..=4) if self.phases[p].rigid[g].satellite.is_some() => Some("parent group"),
            // A TLS group's atoms take their β from its T, L and S (L367).
            Field::AtomBeta(p, a, _) if (0..self.phases[p].rigid.len()).any(|g| crate::rigid_body::tls_members(&self.phases[p].rigid, g).iter().any(|&(h, m)| self.phases[p].rigid[h].members[m].index == a)) => Some("TLS"),
            Field::AtomX(p, a, k) | Field::AtomBeta(p, a, k) => {
                let ph = &self.phases[p];
                let Intensities::Structure(st, _) = &ph.intensities else { return None };
                let site = Site::new(&ph.group, st.atoms[a].x).ok()?;
                let cons = if matches!(f, Field::AtomX(..)) { &site.coordinates } else { &site.adp };
                (!cons.free.contains(&k)).then_some("site symmetry")
            }
            _ => None,
        }
    }

    /// Every editable value with its address, dictionary id, value, limits and what fixes it.
    pub fn fields(&self) -> Vec<FieldInfo> {
        self.all_fields().into_iter().map(|f| FieldInfo { field: f, address: address(f, self), id: dictionary_id(f, self), value: self.get(f), limits: limits(f, self), fixed_by: self.fixed_by(f) }).collect()
    }

    /// Set one value, checked: the field must exist in this model, the value be finite, and the field not follow
    /// from others. The values tied to it follow: the cell parameters its crystal system makes equal, and the
    /// coordinates or betas its site expresses through it.
    pub fn set_field(&mut self, f: Field, v: f64) -> Result<(), String> {
        if !self.all_fields().contains(&f) {
            return Err(format!("{f:?} is not a field of this model"));
        }
        if !v.is_finite() {
            return Err(format!("{}: {v} is not a number", address(f, self)));
        }
        if let Some(by) = self.fixed_by(f) {
            return Err(format!("{} follows from the {by}: set the value it follows from", address(f, self)));
        }
        self.set(f, v);
        match f {
            Field::Cell(p, k) => {
                for (j, rule) in cell_rules(&self.phases[p]).into_iter().enumerate() {
                    if let CellRule::Equal(src) = rule {
                        if src == k {
                            self.set(Field::Cell(p, j), v);
                        }
                    }
                }
            }
            Field::AtomX(p, a, _) | Field::AtomBeta(p, a, _) => {
                let ph = &self.phases[p];
                let Intensities::Structure(st, _) = &ph.intensities else { return Ok(()) };
                let site = Site::new(&ph.group, st.atoms[a].x).map_err(|e| e.to_string())?;
                let is_x = matches!(f, Field::AtomX(..));
                let cons = if is_x { site.coordinates } else { site.adp };
                let field = |i: usize| if is_x { Field::AtomX(p, a, i) } else { Field::AtomBeta(p, a, i) };
                let free: Vec<f64> = cons.free.iter().map(|&i| self.get(field(i))).collect();
                for i in 0..cons.names.len() {
                    if cons.free.contains(&i) {
                        continue;
                    }
                    let value = cons.constant[i].value() + cons.coef[i].iter().zip(&free).map(|(c, x)| c.value() * x).sum::<f64>();
                    self.set(field(i), value);
                }
            }
            _ => {}
        }
        Ok(())
    }

    pub(crate) fn radiation(&self) -> Radiation {
        if self.neutrons {
            Radiation::Neutron
        } else {
            Radiation::Xray { lambda: self.lambda[0], anomalous: true }
        }
    }

    /// Add a phase from a space-group setting and a cell, with no atoms; returns its index. Its profile is the
    /// Thompson–Cox–Hastings with W = 0.01 deg² and the rest zero, its scale 1.
    pub fn add_phase(&mut self, name: &str, group: SpaceGroup, cell: Cell) -> Result<usize, String> {
        self.add_structure(name, Structure { group, cell, atoms: Vec::new() })
    }

    /// Add a phase from a CIF (`cif::structure`); returns its index.
    pub fn add_phase_from_cif(&mut self, name: &str, text: &str) -> Result<usize, String> {
        let st = crate::cif::structure(text)?;
        self.add_structure(name, st)
    }

    /// The user's scattering factors of the pattern (P11): every species whose scatterer is not the tables'. A phase or
    /// atom added takes them, as FullProf's line 12 holds for the whole pattern.
    pub fn user_scatterers(&self) -> Vec<crate::scattering::UserScatterer> {
        let radiation = self.radiation();
        let mut out: Vec<crate::scattering::UserScatterer> = Vec::new();
        for ph in &self.phases {
            if let Intensities::Structure(st, sc) = &ph.intensities {
                for (a, s) in st.atoms.iter().zip(sc) {
                    if !out.iter().any(|u| u.species.eq_ignore_ascii_case(&a.species)) {
                        out.extend(crate::scattering::UserScatterer::of(&a.species, s, radiation));
                    }
                }
            }
        }
        out
    }

    fn add_structure(&mut self, name: &str, st: Structure) -> Result<usize, String> {
        let radiation = self.radiation();
        let user = self.user_scatterers();
        let sc = st.atoms.iter().map(|a| crate::scattering::resolve(&a.species, radiation, &user)).collect::<Result<Vec<_>, _>>()?;
        let cell = st.cell;
        self.phases.push(Phase {
            name: name.to_string(),
            group: st.group.clone(),
            cell,
            intensities: Intensities::Structure(Box::new(st), sc),
            scale: 1.0,
            bov: 0.0,
            shape: Shape::ThompsonCoxHastings,
            eta0: 0.0,
            u: 0.0,
            v: 0.0,
            w: 0.01,
            x: 0.0,
            y: 0.0,
            gauss_size: 0.0,
            asym: [0.0; 4],
            s_l: 0.0,
            d_l: 0.0,
            pref: [if self.march { 1.0 } else { 0.0 }, 0.0],
            pref_direction: [0.0, 0.0, 1.0],
            reflection_limit: None,
            extract: false,
            atz: 0.0,
            brindley: 1.0,
            internal_standard: None,
            valence_sum_range: None,
            modes: None,
            rigid: Vec::new(),
            uvw2: None,
            parent: None,
            quartic_strain: None,
            shape_size: None,
            tof: None,
            texture: None,
        });
        Ok(self.phases.len() - 1)
    }

    /// Remove a phase, with its restraints.
    pub fn remove_phase(&mut self, p: usize) -> Result<(), String> {
        if p >= self.phases.len() {
            return Err(format!("no phase {}", p + 1));
        }
        self.phases.remove(p);
        self.restraints.retain(|r| r.phase != p);
        for r in &mut self.restraints {
            if r.phase > p {
                r.phase -= 1;
            }
        }
        Ok(())
    }

    /// Add an atom to a phase (species as the tables write it: "O", "Ti+4"; FullProf's Occ); returns its index.
    pub fn add_atom(&mut self, p: usize, label: &str, species: &str, x: [f64; 3], biso: f64, occ: f64) -> Result<usize, String> {
        let radiation = self.radiation();
        let user = self.user_scatterers();
        let ph = self.phases.get_mut(p).ok_or(format!("no phase {}", p + 1))?;
        let Intensities::Structure(st, sc) = &mut ph.intensities else { return Err(format!("phase {} has no structure (profile matching)", p + 1)) };
        if st.atoms.iter().any(|a| a.label == label) {
            return Err(format!("phase {}: there is already an atom {label}", p + 1));
        }
        Site::new(&st.group, x).map_err(|e| format!("{label}: {e}"))?;
        let s = crate::scattering::resolve(species, radiation, &user)?;
        st.atoms.push(Atom { label: label.to_string(), species: species.to_string(), x, biso, occ, beta: None, charge: None });
        sc.push(s);
        Ok(st.atoms.len() - 1)
    }

    /// Remove an atom; refused while a restraint names it.
    pub fn remove_atom(&mut self, p: usize, a: usize) -> Result<(), String> {
        let ph = self.phases.get_mut(p).ok_or(format!("no phase {}", p + 1))?;
        let Intensities::Structure(st, sc) = &mut ph.intensities else { return Err(format!("phase {} has no structure", p + 1)) };
        if a >= st.atoms.len() {
            return Err(format!("phase {}: no atom {}", p + 1, a + 1));
        }
        if self.restraints.iter().any(|r| r.phase == p && r.atoms().contains(&a)) {
            return Err(format!("phase {}: atom {} is in a restraint; remove the restraint first", p + 1, st.atoms[a].label));
        }
        st.atoms.remove(a);
        sc.remove(a);
        for r in self.restraints.iter_mut().filter(|r| r.phase == p) {
            let shift = |i: &mut usize| {
                if *i > a {
                    *i -= 1;
                }
            };
            match &mut r.kind {
                Kind::Distance { atoms, .. } => atoms.iter_mut().for_each(shift),
                Kind::Angle { atoms, .. } => atoms.iter_mut().for_each(shift),
            }
        }
        Ok(())
    }

    /// Replace the background model (points in increasing 2θ; at least one value).
    pub fn set_background(&mut self, b: BackgroundModel) -> Result<(), String> {
        match &b {
            BackgroundModel::Points(p) | BackgroundModel::Spline(p) if p.windows(2).any(|w| w[1].0 <= w[0].0) => return Err("background points must be in increasing 2θ".into()),
            BackgroundModel::Chebyshev { from, to, .. } if from >= to => return Err("Chebyshev background: its range is empty".into()),
            _ => {}
        }
        if background_len(&b) == 0 {
            return Err("a background needs at least one value".into());
        }
        self.background = b;
        Ok(())
    }

    /// The number of coefficients of a polynomial, Chebyshev or Fourier background: new ones are 0, extra ones
    /// dropped. Points are added and removed one by one (`set_background`).
    pub fn set_background_terms(&mut self, n: usize) -> Result<(), String> {
        if n == 0 {
            return Err("a background needs at least one value".into());
        }
        match &mut self.background {
            BackgroundModel::Points(_) | BackgroundModel::Spline(_) => Err("a background by points changes by its points, not by a number of terms".into()),
            BackgroundModel::Polynomial { coefficients, .. } | BackgroundModel::Chebyshev { coefficients, .. } | BackgroundModel::Fourier { coefficients } => {
                coefficients.resize(n, 0.0);
                Ok(())
            }
        }
    }

    /// Add an excluded region [a, b] (degrees 2θ); regions are kept in order.
    pub fn add_excluded(&mut self, a: f64, b: f64) -> Result<(), String> {
        if a.partial_cmp(&b) != Some(std::cmp::Ordering::Less) {
            return Err(format!("excluded region: {a} is not below {b}"));
        }
        self.excluded.push((a, b));
        self.excluded.sort_by(|x, y| x.0.total_cmp(&y.0));
        Ok(())
    }

    pub fn remove_excluded(&mut self, i: usize) -> Result<(), String> {
        if i >= self.excluded.len() {
            return Err(format!("no excluded region {}", i + 1));
        }
        self.excluded.remove(i);
        Ok(())
    }

    /// The profile function of a phase. The axial divergence needs the TCH profile: leaving it drops S/L and H/L.
    pub fn set_shape(&mut self, p: usize, shape: Shape) -> Result<(), String> {
        let ph = self.phases.get_mut(p).ok_or(format!("no phase {}", p + 1))?;
        if shape != Shape::ThompsonCoxHastings {
            ph.s_l = 0.0;
            ph.d_l = 0.0;
        }
        ph.shape = shape;
        Ok(())
    }

    /// The asymmetry function of a phase (the axial divergence with the TCH profile only).
    pub fn set_asymmetry(&mut self, p: usize, a: Asymmetry) -> Result<(), String> {
        let ph = self.phases.get_mut(p).ok_or(format!("no phase {}", p + 1))?;
        match a {
            Asymmetry::None => (ph.asym, ph.s_l, ph.d_l) = ([0.0; 4], 0.0, 0.0),
            Asymmetry::BerarBaldinozzi(q) => (ph.asym, ph.s_l, ph.d_l) = (q, 0.0, 0.0),
            Asymmetry::AxialDivergence { s_l, h_l } => {
                if ph.shape != Shape::ThompsonCoxHastings {
                    return Err(format!("phase {}: the axial divergence needs the Thompson–Cox–Hastings profile", p + 1));
                }
                if s_l < 0.0 || h_l < 0.0 {
                    return Err("the axial divergence's S/L and H/L are not negative".into());
                }
                (ph.asym, ph.s_l, ph.d_l) = ([0.0; 4], s_l, h_l);
            }
        }
        Ok(())
    }

    /// The preferred-orientation function of the pattern: March–Dollase (`march`) or the exponential; each phase's G1
    /// is set to its function's "none" (1 for March, 0 for the exponential).
    pub fn set_preferred_orientation_function(&mut self, march: bool) {
        if self.march != march {
            for ph in &mut self.phases {
                ph.pref = [if march { 1.0 } else { 0.0 }, 0.0];
            }
        }
        self.march = march;
    }

    /// A phase's preferred orientation: the axis (reciprocal-lattice components), G1 and G2.
    pub fn set_preferred_orientation(&mut self, p: usize, direction: [f64; 3], g1: f64, g2: f64) -> Result<(), String> {
        let march = self.march;
        let ph = self.phases.get_mut(p).ok_or(format!("no phase {}", p + 1))?;
        if direction.iter().all(|&v| v == 0.0) {
            return Err("the preferred-orientation axis is zero".into());
        }
        if march && g1 <= 0.0 {
            return Err("March's G1 is positive (1: no preferred orientation)".into());
        }
        if !(0.0..=1.0).contains(&g2) {
            return Err("G2, the untextured fraction, is between 0 and 1".into());
        }
        (ph.pref_direction, ph.pref) = (direction, [g1, g2]);
        Ok(())
    }
}

/// A phase's preferred-orientation fields: G1, G2, or its texture's values when it has several axes (P30).
fn pref_fields(ph: &crate::powder::Phase, p: usize) -> Vec<Field> {
    match &ph.texture {
        Some(t) => (0..t.len()).map(|k| Field::Texture(p, k)).collect(),
        None => vec![Field::Pref(p, 0), Field::Pref(p, 1)],
    }
}
