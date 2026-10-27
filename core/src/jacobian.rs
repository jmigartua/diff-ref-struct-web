//! The calculated pattern's analytic derivatives with respect to the model's refinable values (L053).
//!
//! `Field` names one value of the model (the zero, a cell edge, an atom's x, a background height …); `Model::get` and
//! `Model::set` read and write it. `Model::derivatives` returns the pattern of `Model::calculate` and, for each field
//! asked, the derivative of Ycal at every point:
//! - per peak, the position, FWHM, η and integrated intensity are recomputed with dual numbers (`dual::Dual`) over the
//!   global and phase variables (zero, Sycos, Sysin, Ratio, the six cell parameters, U V W X Y I_G η0, G1 G2, B_ov,
//!   scale), following `Model::peaks` formula for formula: the cell enters through the reciprocal metric, so d, θ,
//!   the widths at the Bragg angle, Lp, absorption, the form factors and Debye–Waller terms at s = 1/2d and the
//!   preferred orientation all carry their cell derivatives;
//! - the structure factor's derivatives with respect to s and to each atom's x, y, z, B, Occ and betas are analytic
//!   (`Structure::structure_factor_derivatives`);
//! - at each point, the profile (pseudo-Voigt, TCH, Finger–Cox–Jephcoat, Bérar–Baldinozzi) is evaluated with dual
//!   numbers over its position, FWHM, η, P1–P4, S/L and H/L;
//! - the background's derivatives are its basis functions.
//!
//! The reflections, windows and excluded points are those `calculate` chooses at the current values (their
//! boundaries move discontinuously, so they carry no derivative). Test `engine_derivatives` checks every kind of field
//! against a central finite difference, and the dual values against `calculate`.

use crate::dual::{Dual, Real};
use crate::powder::{BackgroundModel, Calculated, InstrumentalWidths, Intensities, Model, Phase, Shape};
use crate::profile::{berar_baldinozzi, gaussian, lorentzian, pseudo_voigt, taper, tch, Fcj};

/// One refinable value of the model. Phases, atoms and components count from 0.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Field {
    Zero,
    Sycos,
    Sysin,
    Ratio,
    /// The k-th background value: a point's height, or a coefficient of the polynomial, Chebyshev or Fourier series.
    Background(usize),
    /// Value k of the further terms of the background (P23, `background_terms::BackgroundTerms::part`): a Debye-like
    /// term's amplitude and distance, an inverse power's coefficient, a profile's coefficient.
    BackgroundTerm(usize),
    Scale(usize),
    Bov(usize),
    /// Cell parameter k (a, b, c, α, β, γ) of a phase.
    Cell(usize, usize),
    U(usize),
    V(usize),
    W(usize),
    X(usize),
    Y(usize),
    GaussSize(usize),
    Eta0(usize),
    /// Bérar–Baldinozzi P_k (k = 0 … 3).
    Asym(usize, usize),
    SL(usize),
    DL(usize),
    /// Preferred orientation G_k (k = 0, 1).
    Pref(usize, usize),
    /// Value k of a phase's preferred orientation along several axes (P30, `texture::MultiAxial::part`): axis j's
    /// March ratio r_j (2j) and fraction f_j (2j + 1; the last axis's is 1 minus the others).
    Texture(usize, usize),
    /// Atom a's fractional coordinate k of phase p.
    AtomX(usize, usize, usize),
    AtomB(usize, usize),
    AtomOcc(usize, usize),
    /// Atom a's β_k (11, 22, 33, 12, 13, 23).
    AtomBeta(usize, usize, usize),
    /// The amplitude of symmetry mode m of phase p (R29, L090): it moves the atoms' coordinates linearly.
    Mode(usize, usize),
    /// Value k of rigid group g of phase p (P41, `rigid_body::Group::value`): its origin (0–2), θ, φ, χ (3–5), the
    /// internal coordinates of its atom m (6 + 3m + c); it moves the group's atoms.
    Rigid(usize, usize, usize),
    /// The k-th S_HKL of a phase's quartic form of strain, in the order of its terms (P25, L211).
    Quartic(usize, usize),
    /// The Lorentzian share ξ of a phase's quartic strain.
    QuarticShare(usize),
    /// LorSiz of a phase's anisotropic size by shape (platelets, needles) or the antiphase rule (P25, L267).
    ShapeSize(usize),
    /// The wavelength λ₁, and λ₂ with it when the pattern has one wavelength (P18, L233).
    Wavelength,
    /// Bragg–Brentano surface micro-absorption: P0, Cp, τ (k = 0, 1, 2; FullProf's eq. 3.36, P29, L233).
    Micro(usize),
    /// The intensity of reflection k of a phase's list (no structure): a Pawley fit refines it (P36).
    Intensity(usize, usize),
    /// A TOF bank's value k (P05, `tof::BANK_FIELDS`): Dtt1, Dtt2, Dtt_1overd (0–2), and Hodges' Zerot, Dtt1t, Dtt2t,
    /// x-cross, Width (3–7). Its zero is `Zero`.
    Tof(usize),
    /// A phase's TOF profile value k (P05, `tof::PROFILE_FIELDS`): σ₀, σ₁, σ₂, σ_Q (0–3), γ₀, γ₁, γ₂ (4–6), α₀, α₁, α_Q
    /// (7–9), β₀, β₁, β_Q (10–12), α₀ₜ, β₀ₜ, α₁ₜ, β₁ₜ (13–16), the extinction (17), the absorption A₁ (18).
    TofProfile(usize, usize),
}

impl Field {
    /// The same field of phase q (a pattern's fields unchanged).
    pub fn with_phase(self, q: usize) -> Field {
        use Field::*;
        match self {
            Scale(_) => Scale(q),
            Bov(_) => Bov(q),
            Cell(_, k) => Cell(q, k),
            U(_) => U(q),
            V(_) => V(q),
            W(_) => W(q),
            X(_) => X(q),
            Y(_) => Y(q),
            GaussSize(_) => GaussSize(q),
            Eta0(_) => Eta0(q),
            Asym(_, k) => Asym(q, k),
            SL(_) => SL(q),
            DL(_) => DL(q),
            Pref(_, k) => Pref(q, k),
            Texture(_, k) => Texture(q, k),
            AtomX(_, a, k) => AtomX(q, a, k),
            AtomB(_, a) => AtomB(q, a),
            AtomOcc(_, a) => AtomOcc(q, a),
            AtomBeta(_, a, k) => AtomBeta(q, a, k),
            Mode(_, m) => Mode(q, m),
            Rigid(_, g, k) => Rigid(q, g, k),
            Quartic(_, k) => Quartic(q, k),
            QuarticShare(_) => QuarticShare(q),
            ShapeSize(_) => ShapeSize(q),
            Intensity(_, k) => Intensity(q, k),
            TofProfile(_, k) => TofProfile(q, k),
            f => f,
        }
    }

    /// Whether the field is of a phase's structure (its atoms or symmetry modes), which a phase keeps whatever pattern
    /// it is seen in, against its scale, profile and cell, which are each pattern's own (P06).
    pub fn is_structural(self) -> bool {
        matches!(self, Field::AtomX(..) | Field::AtomB(..) | Field::AtomOcc(..) | Field::AtomBeta(..) | Field::Mode(..) | Field::Rigid(..))
    }

    /// The phase a field belongs to (none for the pattern's fields).
    pub fn phase(self) -> Option<usize> {
        use Field::*;
        match self {
            Zero | Sycos | Sysin | Ratio | Background(_) | BackgroundTerm(_) | Wavelength | Micro(_) | Tof(_) => None,
            Scale(p) | Bov(p) | Cell(p, _) | U(p) | V(p) | W(p) | X(p) | Y(p) | GaussSize(p) | Eta0(p) | Asym(p, _) | SL(p) | DL(p) | Pref(p, _) | Texture(p, _) | AtomX(p, ..) | AtomB(p, _) | AtomOcc(p, _) | AtomBeta(p, ..) | Mode(p, _) | Rigid(p, ..) | Quartic(p, _) | QuarticShare(p) | ShapeSize(p) | Intensity(p, _) | TofProfile(p, _) => Some(p),
        }
    }
}

/// The global and phase variables carried as duals: their index in the phase dual.
const ZERO: usize = 0;
const SYCOS: usize = 1;
const SYSIN: usize = 2;
const RATIO: usize = 3;
const CELL: usize = 4;
const U: usize = 10;
const V: usize = 11;
const W: usize = 12;
const X: usize = 13;
const Y: usize = 14;
const IG: usize = 15;
const ETA0: usize = 16;
const G1: usize = 17;
const G2: usize = 18;
const BOV: usize = 19;
const SCALE: usize = 20;
/// The quartic form's sum Σ S_HKL h^H k^K l^L of the reflection in hand (a variable per peak: ∂/∂S_HKL is its
/// monomial times ∂/∂SUM), and the Lorentzian share ξ.
const SUM: usize = 21;
const XI: usize = 22;
/// The wavelength λ₁ and the micro-absorption's P0, Cp, τ.
const WAVE: usize = 23;
const MICRO: usize = 24;
/// LorSiz of the anisotropic size (P25, L267).
const SIZE: usize = 27;
const NV: usize = 28;
type PD = Dual<NV>;

/// The profile's variables: position, FWHM, η, P1–P4, S/L, H/L, and the Bragg angle 2θ_B at which the asymmetry is
/// taken (L367).
const NP: usize = 10;
type QD = Dual<NP>;

/// A TOF pattern's variables carried as duals (P05): the zero, the bank's eight values (`tof::BANK_FIELDS`), the cell,
/// scale, Bov, G1 and G2, the phase's nineteen TOF values (`tof::PROFILE_FIELDS`), and the reflection's F² (through
/// which the atoms reach the extinction).
const T_ZERO: usize = 0;
const T_BANK: usize = 1;
const T_CELL: usize = 9;
const T_SCALE: usize = 15;
const T_BOV: usize = 16;
const T_G: usize = 17;
const T_PROFILE: usize = 19;
const T_F2: usize = 38;
const NT: usize = 39;
type TD = Dual<NT>;

impl Field {
    /// The TOF variable this field is, for a phase p of a TOF pattern.
    fn tof_var(&self, p: usize) -> Option<usize> {
        match *self {
            Field::Zero => Some(T_ZERO),
            Field::Tof(k) => Some(T_BANK + k),
            Field::Cell(q, k) if q == p => Some(T_CELL + k),
            Field::Scale(q) if q == p => Some(T_SCALE),
            Field::Bov(q) if q == p => Some(T_BOV),
            Field::Pref(q, k) if q == p => Some(T_G + k),
            Field::TofProfile(q, k) if q == p => Some(T_PROFILE + k),
            _ => None,
        }
    }

    /// The phase variable this field is, for a phase p (None for the others).
    fn phase_var(&self, p: usize) -> Option<usize> {
        match *self {
            Field::Zero => Some(ZERO),
            Field::Sycos => Some(SYCOS),
            Field::Sysin => Some(SYSIN),
            Field::Ratio => Some(RATIO),
            Field::Wavelength => Some(WAVE),
            Field::Micro(k) => Some(MICRO + k),
            Field::Cell(q, k) if q == p => Some(CELL + k),
            Field::U(q) if q == p => Some(U),
            Field::V(q) if q == p => Some(V),
            Field::W(q) if q == p => Some(W),
            Field::X(q) if q == p => Some(X),
            Field::Y(q) if q == p => Some(Y),
            Field::GaussSize(q) if q == p => Some(IG),
            Field::Eta0(q) if q == p => Some(ETA0),
            Field::Pref(q, k) if q == p => Some(G1 + k),
            Field::Bov(q) if q == p => Some(BOV),
            Field::Scale(q) if q == p => Some(SCALE),
            Field::QuarticShare(q) if q == p => Some(XI),
            Field::ShapeSize(q) if q == p => Some(SIZE),
            _ => None,
        }
    }

    /// The profile variable this field is, for a phase p.
    fn profile_var(&self, p: usize) -> Option<usize> {
        match *self {
            Field::Asym(q, k) if q == p => Some(3 + k),
            Field::SL(q) if q == p => Some(7),
            Field::DL(q) if q == p => Some(8),
            _ => None,
        }
    }
}

fn cell_params(ph: &Phase) -> [f64; 6] {
    let c = &ph.cell;
    [c.a, c.b, c.c, c.alpha, c.beta, c.gamma]
}

fn background_values(b: &BackgroundModel) -> Vec<f64> {
    match b {
        BackgroundModel::Points(p) | BackgroundModel::Spline(p) => p.iter().map(|q| q.1).collect(),
        BackgroundModel::Polynomial { coefficients, .. } | BackgroundModel::Chebyshev { coefficients, .. } | BackgroundModel::Fourier { coefficients } => coefficients.clone(),
    }
}

/// The derivative of the background at t with respect to its k-th value (the background is linear in its values).
fn background_basis(b: &BackgroundModel, k: usize, t: f64) -> f64 {
    match b {
        BackgroundModel::Spline(p) => crate::background::weights(p, crate::background::Interpolation::Spline, t)[k],
        BackgroundModel::Points(p) => {
            if p.is_empty() {
                return 0.0;
            }
            if t <= p[0].0 {
                return if k == 0 { 1.0 } else { 0.0 };
            }
            for (i, w) in p.windows(2).enumerate() {
                if t <= w[1].0 {
                    let f = (t - w[0].0) / (w[1].0 - w[0].0);
                    return if k == i { 1.0 - f } else if k == i + 1 { f } else { 0.0 };
                }
            }
            if k == p.len() - 1 {
                1.0
            } else {
                0.0
            }
        }
        BackgroundModel::Polynomial { origin, .. } => (t / origin - 1.0).powi(k as i32),
        BackgroundModel::Fourier { .. } => (k as f64 * t.to_radians()).cos(),
        BackgroundModel::Chebyshev { from, to, .. } => {
            let x = (2.0 * (t - from) / (to - from) - 1.0).clamp(-1.0, 1.0);
            let (mut t0, mut t1) = (1.0, x);
            match k {
                0 => 1.0,
                1 => x,
                _ => {
                    let mut tk = 0.0;
                    for _ in 2..=k {
                        tk = 2.0 * x * t1 - t0;
                        t0 = t1;
                        t1 = tk;
                    }
                    tk
                }
            }
        }
    }
}

/// The direct metric from the six cell parameters (Å, degrees).
fn direct_metric<T: Real>(c: &[T; 6]) -> [[T; 3]; 3] {
    let (a, b, cc) = (c[0], c[1], c[2]);
    let (ca, cb, cg) = (c[3].to_radians().cos(), c[4].to_radians().cos(), c[5].to_radians().cos());
    [[a * a, a * b * cg, a * cc * cb], [a * b * cg, b * b, b * cc * ca], [a * cc * cb, b * cc * ca, cc * cc]]
}

/// The reciprocal metric from the six cell parameters (Å, degrees).
fn reciprocal_metric<T: Real>(c: &[T; 6]) -> [[T; 3]; 3] {
    let g = direct_metric(c);
    let co = |i: usize, j: usize| {
        let (r0, r1) = ([1, 0, 0][i], [2, 2, 1][i]);
        let (c0, c1) = ([1, 0, 0][j], [2, 2, 1][j]);
        g[r0][c0] * g[r1][c1] - g[r0][c1] * g[r1][c0]
    };
    let det = g[0][0] * co(0, 0) - g[0][1] * co(0, 1) + g[0][2] * co(0, 2);
    let mut inv = [[T::cst(0.0); 3]; 3];
    for (i, row) in inv.iter_mut().enumerate() {
        for (j, v) in row.iter_mut().enumerate() {
            let sign = if (i + j) % 2 == 0 { 1.0 } else { -1.0 };
            *v = co(j, i) * sign / det;
        }
    }
    inv
}

fn quad<T: Real>(g: &[[T; 3]; 3], a: [f64; 3], b: [f64; 3]) -> T {
    let mut s = T::cst(0.0);
    for i in 0..3 {
        for j in 0..3 {
            if a[i] != 0.0 && b[j] != 0.0 {
                s = s + g[i][j] * (a[i] * b[j]);
            }
        }
    }
    s
}

/// Bragg–Brentano surface micro-absorption S_r = 1 − P(h), P(h) = P0 + Cp (τ/sin θ)(1 − τ/sin θ) (FullProf's eqs. 3.35
/// and 3.36, manual p. 55), as `Model::micro_absorption`.
fn micro_g(th: PD, [p0, cp, tau]: [PD; 3]) -> PD {
    let r = tau / th.sin();
    -(p0 + cp * r * (-r + 1.0)) + 1.0
}

/// The phase variables as duals, seeded for the fields asked.
struct Vars {
    v: [PD; NV],
}

impl Model {
    /// The value of a field.
    pub fn get(&self, f: Field) -> f64 {
        let ph = |p: usize| &self.phases[p];
        let atom = |p: usize, a: usize| match &self.phases[p].intensities {
            Intensities::Structure(st, _) => &st.atoms[a],
            Intensities::List(_) => panic!("phase {} has no atoms", p + 1),
        };
        match f {
            Field::Zero => self.zero,
            Field::Sycos => self.sycos,
            Field::Sysin => self.sysin,
            Field::Ratio => self.ratio,
            Field::Background(k) => background_values(&self.background)[k],
            Field::BackgroundTerm(k) => self.background_terms.value(k),
            Field::Scale(p) => ph(p).scale,
            Field::Bov(p) => ph(p).bov,
            Field::Cell(p, k) => cell_params(ph(p))[k],
            Field::U(p) => ph(p).u,
            Field::V(p) => ph(p).v,
            Field::W(p) => ph(p).w,
            Field::X(p) => ph(p).x,
            Field::Y(p) => ph(p).y,
            Field::GaussSize(p) => ph(p).gauss_size,
            Field::Eta0(p) => ph(p).eta0,
            Field::Asym(p, k) => ph(p).asym[k],
            Field::SL(p) => ph(p).s_l,
            Field::DL(p) => ph(p).d_l,
            Field::Pref(p, k) => ph(p).pref[k],
            Field::Texture(p, k) => ph(p).texture.as_ref().map_or(0.0, |t| t.value(k)),
            Field::AtomX(p, a, k) => atom(p, a).x[k],
            Field::AtomB(p, a) => atom(p, a).biso,
            Field::AtomOcc(p, a) => atom(p, a).occ,
            Field::AtomBeta(p, a, k) => atom(p, a).beta.map_or(0.0, |b| b[k]),
            Field::Mode(p, m) => ph(p).modes.as_ref().map_or(0.0, |ms| ms.modes[m].amplitude),
            Field::Rigid(p, g, k) => ph(p).rigid.get(g).map_or(0.0, |gr| gr.value(k)),
            Field::Quartic(p, k) => ph(p).quartic_strain.as_ref().map_or(0.0, |q| q.terms[k].1),
            Field::QuarticShare(p) => ph(p).quartic_strain.as_ref().map_or(0.0, |q| q.lorentzian),
            Field::ShapeSize(p) => ph(p).shape_size.map_or(0.0, |s| s.value),
            Field::Wavelength => self.lambda[0],
            Field::Micro(k) => self.microabsorption[k],
            Field::Tof(k) => self.tof.as_ref().map_or(0.0, |b| b.get(k)),
            Field::TofProfile(p, k) => ph(p).tof.as_ref().map_or(0.0, |t| t.get(k)),
            Field::Intensity(p, k) => match &ph(p).intensities {
                Intensities::List(l) => l[k].2,
                Intensities::Structure(..) => panic!("phase {} has a structure: no list of intensities", p + 1),
            },
        }
    }

    /// Set a field. A cell parameter is set in the phase's cell and in its structure's.
    pub fn set(&mut self, f: Field, v: f64) {
        match f {
            Field::Zero => self.zero = v,
            Field::Sycos => self.sycos = v,
            Field::Sysin => self.sysin = v,
            Field::Ratio => self.ratio = v,
            // One wavelength: λ₂ = λ₁ moves with it, as FullProf writes both (FP_tests/pb_ho.pcr and eight more patterns
            // whose line 14 refines Lambda: λ₁ = λ₂ = Lambda).
            Field::Wavelength => {
                if self.lambda[1] == self.lambda[0] {
                    self.lambda[1] = v;
                }
                self.lambda[0] = v;
            }
            Field::Micro(k) => self.microabsorption[k] = v,
            Field::Tof(k) => self.tof.as_mut().unwrap_or_else(|| panic!("the pattern is not TOF")).set(k, v),
            Field::TofProfile(p, k) => self.phases[p].tof.as_mut().unwrap_or_else(|| panic!("phase {} has no TOF profile", p + 1)).set(k, v),
            Field::Intensity(p, k) => match &mut self.phases[p].intensities {
                Intensities::List(l) => l[k].2 = v,
                Intensities::Structure(..) => panic!("phase {} has a structure: no list of intensities", p + 1),
            },
            Field::Background(k) => match &mut self.background {
                BackgroundModel::Points(p) | BackgroundModel::Spline(p) => p[k].1 = v,
                BackgroundModel::Polynomial { coefficients, .. } | BackgroundModel::Chebyshev { coefficients, .. } | BackgroundModel::Fourier { coefficients } => coefficients[k] = v,
            },
            Field::BackgroundTerm(k) => self.background_terms.set(k, v),
            Field::Scale(p) => self.phases[p].scale = v,
            Field::Bov(p) => self.phases[p].bov = v,
            Field::Cell(p, k) => {
                let ph = &mut self.phases[p];
                let c = &mut ph.cell;
                *[&mut c.a, &mut c.b, &mut c.c, &mut c.alpha, &mut c.beta, &mut c.gamma][k] = v;
                if let Intensities::Structure(st, _) = &mut ph.intensities {
                    st.cell = ph.cell;
                }
                // A rigid group keeps its shape in Å: its atoms' fractions follow the cell (P41).
                place_rigid(ph);
            }
            Field::U(p) => self.phases[p].u = v,
            Field::V(p) => self.phases[p].v = v,
            Field::W(p) => self.phases[p].w = v,
            Field::X(p) => self.phases[p].x = v,
            Field::Y(p) => self.phases[p].y = v,
            Field::GaussSize(p) => self.phases[p].gauss_size = v,
            Field::Eta0(p) => self.phases[p].eta0 = v,
            Field::Asym(p, k) => self.phases[p].asym[k] = v,
            Field::SL(p) => self.phases[p].s_l = v,
            Field::DL(p) => self.phases[p].d_l = v,
            Field::Pref(p, k) => self.phases[p].pref[k] = v,
            Field::Texture(p, k) => self.phases[p].texture.as_mut().unwrap_or_else(|| panic!("phase {} has no preferred orientation along several axes", p + 1)).set(k, v),
            Field::Quartic(p, k) => self.phases[p].quartic_strain.as_mut().unwrap_or_else(|| panic!("phase {} has no quartic strain", p + 1)).terms[k].1 = v,
            Field::QuarticShare(p) => self.phases[p].quartic_strain.as_mut().unwrap_or_else(|| panic!("phase {} has no quartic strain", p + 1)).lorentzian = v,
            Field::ShapeSize(p) => self.phases[p].shape_size.as_mut().unwrap_or_else(|| panic!("phase {} has no anisotropic size", p + 1)).value = v,
            Field::Mode(p, m) => {
                let ph = &mut self.phases[p];
                let Some(ms) = &mut ph.modes else { panic!("phase {} has no symmetry modes", p + 1) };
                ms.modes[m].amplitude = v;
                if let Intensities::Structure(st, _) = &mut ph.intensities {
                    for (at, x) in st.atoms.iter_mut().zip(ms.positions()) {
                        at.x = x;
                    }
                }
            }
            Field::Rigid(p, g, k) => {
                let ph = &mut self.phases[p];
                ph.rigid[g].set_value(k, v);
                place_rigid(ph);
            }
            Field::AtomX(p, a, _) | Field::AtomB(p, a) | Field::AtomOcc(p, a) | Field::AtomBeta(p, a, _) => {
                let Intensities::Structure(st, _) = &mut self.phases[p].intensities else { panic!("phase {} has no atoms", p + 1) };
                let at = &mut st.atoms[a];
                match f {
                    Field::AtomX(_, _, k) => at.x[k] = v,
                    Field::AtomB(..) => at.biso = v,
                    Field::AtomOcc(..) => at.occ = v,
                    Field::AtomBeta(_, _, k) => at.beta.get_or_insert([0.0; 6])[k] = v,
                    _ => unreachable!(),
                }
            }
        }
    }

    fn vars(&self, p: usize, fields: &[Field]) -> Vars {
        let ph = &self.phases[p];
        let c = cell_params(ph);
        let mut val = [0.0; NV];
        val[ZERO] = self.zero;
        val[SYCOS] = self.sycos;
        val[SYSIN] = self.sysin;
        val[RATIO] = self.ratio;
        val[CELL..CELL + 6].copy_from_slice(&c);
        (val[U], val[V], val[W], val[X], val[Y], val[IG], val[ETA0]) = (ph.u, ph.v, ph.w, ph.x, ph.y, ph.gauss_size, ph.eta0);
        (val[G1], val[G2], val[BOV], val[SCALE]) = (ph.pref[0], ph.pref[1], ph.bov, ph.scale);
        val[XI] = ph.quartic_strain.as_ref().map_or(0.0, |q| q.lorentzian);
        val[SIZE] = ph.shape_size.map_or(0.0, |s| s.value);
        val[WAVE] = self.lambda[0];
        val[MICRO..MICRO + 3].copy_from_slice(&self.microabsorption);
        let mut v = val.map(PD::cst);
        for f in fields {
            if let Some(i) = f.phase_var(p) {
                v[i] = PD::var(val[i], i);
            }
        }
        Vars { v }
    }

    /// X-ray or neutron Lorentz–polarisation at 2θ (degrees), as `lorentz_polarisation`.
    fn lp_g<T: Real>(&self, two_theta: T) -> T {
        let th = (two_theta / 2.0).to_radians();
        let c2 = (th * 2.0).cos();
        match self.ilo {
            -1 => T::cst(1.0),
            2 => {
                let k = self.rpolarz;
                let (c1, c2m) = (self.cthm.sqrt(), self.cthm);
                let p = if self.neutrons { T::cst(1.0) } else { (c2 * c2 * c1 + 1.0) * (k / (1.0 + c1)) + (c2 * c2 * c2m + 1.0) * ((1.0 - k) / (1.0 + c2m)) };
                p / (th * 2.0).sin()
            }
            _ => {
                let pol = if self.neutrons {
                    T::cst(1.0)
                } else if self.ilo == 3 {
                    // Cthm = 0 taken as cos²2θ_M = 1, as `lorentz_polarisation` (L322).
                    let m = if self.cthm == 0.0 { 1.0 } else { self.cthm };
                    c2 * c2 * (self.rpolarz * m) + (1.0 - self.rpolarz)
                } else if self.convention.lorentz_polarisation {
                    c2 * c2 * self.cthm + 1.0
                } else {
                    (c2 * c2 * self.cthm + 1.0) * 0.5
                };
                pol / (th.sin().powi(2) * th.cos() * 2.0)
            }
        }
    }

    fn absorption_g<T: Real>(&self, th: T) -> T {
        // (The surface micro-absorption, which multiplies it, is `micro_g`, over its own variables.)
        if self.mur == 0.0 {
            return T::cst(1.0);
        }
        let s2 = th.sin().powi(2);
        if self.ilo == 2 {
            (T::cst(-self.mur) / th.cos()).exp() / th.cos()
        } else {
            ((s2 * 0.0368 - 1.7133) * self.mur + (s2 * 0.375 + 0.0927) * (self.mur * self.mur)).exp()
        }
    }

    /// The instrument's (H_G², H_L) at 2θ, as `InstrumentalWidths::at`.
    fn instrumental_g<T: Real>(&self, two_theta: T) -> (T, T) {
        let th = (two_theta / 2.0).to_radians();
        match &self.instrumental {
            InstrumentalWidths::None => (T::cst(0.0), T::cst(0.0)),
            InstrumentalWidths::Function { res, uvwxyz: [u, v, w, x, y, z] } => match res {
                1 => ((th.tan() * *u + *v) * th.tan() + *w, th.tan() * *x + T::cst(*y) / th.cos() + *z),
                2 => ((th.tan() * *u + *v) * th.tan() + *w, (two_theta * *x + *y) * two_theta + *z),
                _ => ((two_theta * *u + *v) * two_theta + *w, (two_theta * *x + *y) * two_theta + *z),
            },
            InstrumentalWidths::Table(t) => {
                let i = t.iter().position(|r| r[0] >= two_theta.val()).unwrap_or(t.len() - 1).max(1);
                let (a, b) = (t[i - 1], t[i]);
                let f = if b[0] > a[0] { (two_theta - a[0]) / (b[0] - a[0]) } else { T::cst(0.0) };
                let hg = f * (b[1] - a[1]) + a[1];
                (hg * hg, f * (b[2] - a[2]) + a[2])
            }
        }
    }

    /// FWHM and η at the Bragg angle 2θ, as `peak_width`, with the anisotropic strain's (H_G², H_L) added.
    fn width_g(&self, shape: Shape, v: &[PD; NV], two_theta: PD, uvw2: Option<[f64; 3]>, extra: (PD, PD)) -> (PD, PD) {
        let th = (two_theta / 2.0).to_radians();
        let (t, c) = (th.tan(), th.cos());
        let (ig2, il) = self.instrumental_g(two_theta);
        // The second wavelength's own U, V, W (FullProf's Ratio < 0) are constants: not refined.
        let (u, vv, w) = match uvw2 {
            Some([a, b, c2]) => (PD::cst(a), PD::cst(b), PD::cst(c2)),
            None => (v[U], v[V], v[W]),
        };
        let hg2 = (u * t + vv) * t + w + v[IG] / (c * c) + ig2 + extra.0;
        match shape {
            Shape::ThompsonCoxHastings => tch(hg2.max0().sqrt(), v[X] * t + v[Y] / c + il + extra.1),
            Shape::PseudoVoigt => (hg2.max0().sqrt(), v[ETA0] + v[X] * two_theta),
            Shape::Gaussian => (hg2.max0().sqrt(), PD::cst(0.0)),
            Shape::Lorentzian => (hg2.max0().sqrt(), PD::cst(1.0)),
        }
    }

    /// The preferred-orientation factor, as `preferred_orientation_at`, with the metric and G1, G2 as duals; θ (radians)
    /// for an anvil cell's texture (P30).
    fn pref_g<T: Real>(&self, ph: &Phase, gs: &[[T; 3]; 3], g1: T, g2: T, h: [f64; 3], theta: f64) -> T {
        let p = ph.pref_direction;
        let eq = crate::reflections::equivalents(&ph.group, h.map(|x| x.round() as i32));
        if let Some(t) = &ph.texture {
            let psi = t.psi(theta);
            let mut sum = T::cst(0.0);
            for k in &eq {
                let k = k.map(|x| x as f64);
                let kk = quad(gs, k, k);
                let cos: Vec<T> = t.axes.iter().map(|a| quad(gs, k, a.direction) / (kk * quad(gs, a.direction, a.direction)).sqrt()).collect();
                sum = sum + t.factor_g(&cos, psi);
            }
            return sum / eq.len() as f64;
        }
        let pp = quad(gs, p, p);
        let mut sum = T::cst(0.0);
        for k in &eq {
            let k = k.map(|x| x as f64);
            let cosa = quad(gs, k, p) / (quad(gs, k, k) * pp).sqrt();
            let cosa = if cosa.val() > 1.0 { T::cst(1.0) } else if cosa.val() < -1.0 { T::cst(-1.0) } else { cosa };
            let term = if self.march {
                let sin2 = -(cosa * cosa) + 1.0;
                ((g1 * cosa).powi(2) + sin2 / g1).powf(-1.5)
            } else {
                let c = cosa.abs();
                let alpha = if c.val() > 1.0 - 1e-12 { T::cst(0.0) } else { c.acos() };
                let g = if self.convention.preferred_orientation { -g1 } else { g1 };
                (-g * alpha * alpha).exp()
            };
            sum = sum + (-g2 + 1.0) * term + g2;
        }
        sum / eq.len() as f64
    }

    /// (column, dP/dt ÷ P) for the texture's fields of phase ip among the fields, at reflection h and θ (P30).
    fn texture_derivatives(&self, ph: &Phase, ip: usize, fields: &[Field], h: [f64; 3], theta: f64) -> Vec<(usize, f64)> {
        let Some(t) = &ph.texture else { return Vec::new() };
        let cols: Vec<(usize, usize)> = fields.iter().enumerate().filter_map(|(j, f)| match *f {
            Field::Texture(p, k) if p == ip => Some((j, k)),
            _ => None,
        }).collect();
        if cols.is_empty() {
            return Vec::new();
        }
        let psi = t.psi(theta);
        let cos = crate::powder::texture_cosines(ph, h);
        let po = cos.iter().map(|c| t.factor(c, psi)).sum::<f64>();
        if po == 0.0 {
            return Vec::new();
        }
        cols.into_iter().map(|(j, k)| (j, cos.iter().map(|c| t.derivative(k, c, psi)).sum::<f64>() / po)).collect()
    }

    /// The pattern of `calculate` at x, and the derivative of Ycal at every point with respect to each field.
    pub fn derivatives(&self, x: &[f64], fields: &[Field]) -> (Calculated, Vec<Vec<f64>>) {
        let c = self.calculate(x);
        let mut jac = vec![vec![0.0; x.len()]; fields.len()];
        // The background; none where FullProf's bounds set it to 0 (`Convention::background_bounds`).
        let bounded: Vec<bool> = if self.background_bounded() { self.background_unbounded(x).iter().map(|&b| self.background_out_of_bounds(b)).collect() } else { vec![false; x.len()] };
        // A file's b1 scales the function (eq. 3.6; b3, b4 in `background_x`).
        let scale = self.background_terms.file.as_ref().map_or(1.0, |f| f.scale);
        for (j, f) in fields.iter().enumerate() {
            if let Field::Background(k) = f {
                for (i, &t) in x.iter().enumerate() {
                    let u = self.background_x(t);
                    jac[j][i] = scale
                        * match self.background_edge_segment(x, t) {
                            Some([(x0, _), (x1, _)]) => {
                                let f = (u - x0) / (x1 - x0);
                                [1.0 - f, f].get(*k).copied().unwrap_or(0.0)
                            }
                            None => background_basis(&self.background, *k, u),
                        };
                }
            }
            if let Field::BackgroundTerm(k) = f {
                let file = match self.background_terms.part(*k) {
                    Some((crate::background_terms::Part::File, _, b)) => Some(b),
                    _ => None,
                };
                for (i, &t) in x.iter().enumerate() {
                    jac[j][i] = match file {
                        Some(b) => self.background_file_derivative(b, x, t),
                        None => self.background_terms.derivative(*k, t, self.background_q(t)),
                    };
                }
            }
            if matches!(f, Field::Background(_) | Field::BackgroundTerm(_)) {
                for (i, _) in bounded.iter().enumerate().filter(|b| *b.1) {
                    jac[j][i] = 0.0;
                }
            }
        }
        if let Some(bank) = &self.tof {
            self.tof_derivatives(bank, x, fields, &c, &mut jac);
            return (c, jac);
        }
        let conv = self.convention;
        let excluded = self.excluded_points(x);
        for (ip, ph) in self.phases.iter().enumerate() {
            let vars = self.vars(ip, fields);
            let v = &vars.v;
            let gs = reciprocal_metric(&[v[CELL], v[CELL + 1], v[CELL + 2], v[CELL + 3], v[CELL + 4], v[CELL + 5]]);
            // The fields that touch this phase: (column, phase variable), (column, profile variable), atoms.
            let pv: Vec<(usize, usize)> = fields.iter().enumerate().filter_map(|(j, f)| f.phase_var(ip).map(|k| (j, k))).collect();
            let qv: Vec<(usize, usize)> = fields.iter().enumerate().filter_map(|(j, f)| f.profile_var(ip).map(|k| (j, k))).collect();
            // (column, term) of the quartic form's S_HKL: through SUM, times the reflection's monomial.
            let sv: Vec<(usize, usize)> = fields.iter().enumerate().filter_map(|(j, f)| match *f {
                Field::Quartic(p, k) if p == ip => Some((j, k)),
                _ => None,
            }).collect();
            // (column, atom, atom parameter, coefficient): a mode amplitude moves the coordinates by c·V (chain rule).
            let mut av: Vec<(usize, usize, usize, f64)> = fields
                .iter()
                .enumerate()
                .flat_map(|(j, f)| match *f {
                    Field::AtomX(p, a, k) if p == ip => vec![(j, a, k, 1.0)],
                    Field::AtomB(p, a) if p == ip => vec![(j, a, 3, 1.0)],
                    Field::AtomOcc(p, a) if p == ip => vec![(j, a, 4, 1.0)],
                    Field::AtomBeta(p, a, k) if p == ip => vec![(j, a, 5 + k, 1.0)],
                    Field::Mode(p, m) if p == ip => ph.modes.as_ref().map_or(Vec::new(), |ms| ms.derivative(m).into_iter().map(|(a, k, c)| (j, a, k, c)).collect()),
                    // A rigid group's value moves its atoms and its satellites'; the cell moves them too, the groups keeping their
                    // shape (P41, `rigid_body::placement_derivative`).
                    Field::Rigid(p, g, k) if p == ip => crate::rigid_body::placement_derivative(&ph.rigid, &ph.cell, Some((g, k)), None).into_iter().map(|(a, c, d)| (j, a, c, d)).collect(),
                    Field::Cell(p, k) if p == ip && !ph.rigid.is_empty() => crate::rigid_body::placement_derivative(&ph.rigid, &ph.cell, None, Some(k)).into_iter().map(|(a, c, d)| (j, a, c, d)).collect(),
                    _ => Vec::new(),
                })
                .collect();
            // A TLS group's atoms' β (P41, L367): through their positions, the cell and the TLS values.
            av.extend(tls_terms(ph, ip, fields, &av));
            for pk in c.peaks.iter().filter(|p| p.phase == ip) {
                let lambda = if pk.wavelength == 1 { v[WAVE] } else { PD::cst(self.lambda[1]) };
                let h = pk.hkl.map(|q| q as f64);
                let inv_d2 = quad(&gs, h, h);
                let s = inv_d2.sqrt() * 0.5;
                let th = (s * lambda).asin();
                let pos = th.to_degrees() * 2.0 + v[ZERO] + v[SYCOS] * (if self.convention.displacement { (th * 2.0).cos() } else { th.cos() }) + v[SYSIN] * (th * 2.0).sin();
                let bragg = th.to_degrees() * 2.0;
                // The quartic strain (P25): D_ST tan θ = c·√SUM/M·tan θ, M = 1/d² with the cell's derivatives.
                let (extra, monomials) = match &ph.quartic_strain {
                    Some(q) => {
                        let sum = q.sum(h);
                        let sum = if sv.is_empty() || sum <= 0.0 { PD::cst(sum) } else { PD::var(sum, SUM) };
                        // The Gaussian part, [(1 − ξ)c tan θ/M]²·SUM, is linear in SUM and keeps its derivative where SUM
                        // is 0 (an S_HKL at 0 for the (0 0 l) of mmm); the Lorentzian part ξ c √SUM tan θ/M has none
                        // there (an infinite slope on one side), and is taken as constant.
                        let k = th.tan() / inv_d2 * crate::microstructure::STRAIN_WIDTH;
                        let gk = (-v[XI] + 1.0) * k;
                        if sum.val() < 0.0 {
                            ((PD::cst(0.0), PD::cst(0.0)), Vec::new())
                        } else {
                            let l = if sum.val() > 0.0 { v[XI] * k * sum.sqrt() } else { PD::cst(0.0) };
                            ((gk * gk * sum, l), (0..q.terms.len()).map(|i| q.term_value(i, h)).collect::<Vec<f64>>())
                        }
                    }
                    None => ((PD::cst(0.0), PD::cst(0.0)), Vec::new()),
                };
                // The anisotropic size (L267): k_s λ LorSiz g / cos θ, λ the reflection's (2 d sin θ).
                let extra = match &ph.shape_size {
                    Some(sz) => {
                        // g with the cell's derivatives: cos φ = h·t / (|h*| |t|), |t|² from the direct metric.
                        let g = match sz.shape {
                            crate::microstructure::SizeShape::Platelet(t) | crate::microstructure::SizeShape::Needle(t) => {
                                let ht: f64 = (0..3).map(|i| h[i] * t[i]).sum();
                                let tt = quad(&direct_metric(&[v[CELL], v[CELL + 1], v[CELL + 2], v[CELL + 3], v[CELL + 4], v[CELL + 5]]), t, t);
                                let cos = PD::cst(ht) / (inv_d2 * tt).sqrt();
                                let sin2 = PD::cst(1.0) - cos * cos;
                                match sz.shape {
                                    crate::microstructure::SizeShape::Platelet(_) => cos.abs(),
                                    // sin φ = 0 (h along the needles) has no slope from one side: none taken there.
                                    _ if sin2.val() <= 1e-12 => PD::cst(0.0),
                                    _ => sin2.sqrt(),
                                }
                            }
                            crate::microstructure::SizeShape::Antiphase => PD::cst(sz.factor(h, &ph.cell)),
                        };
                        (extra.0, extra.1 + v[SIZE] * lambda * g * crate::microstructure::SIZE_WIDTH / th.cos())
                    }
                    None => extra,
                };
                let (fwhm, eta) = self.width_g(ph.shape, v, bragg, if pk.wavelength == 2 { ph.uvw2 } else { None }, extra);
                let two = self.ratio != 0.0 && (self.lambda[1] - self.lambda[0]).abs() > 1e-9;
                let weight = if pk.wavelength == 1 { PD::cst(1.0) } else { v[RATIO] };
                // Intensity, and its derivatives with respect to the atoms' parameters.
                let (intensity, datoms): (PD, Vec<(usize, f64)>) = match &ph.intensities {
                    Intensities::Structure(st, sc) => {
                        let (f2, df2ds, dfa) = st.intensity_derivatives(sc, pk.hkl, s.val());
                        let pre = v[SCALE] * weight * pk.multiplicity as f64 * self.lp_g(bragg) * (-v[BOV] * s * s * 2.0).exp() * self.pref_g(ph, &gs, v[G1], v[G2], h, th.val()) * self.absorption_g(th) * micro_g(th, [v[MICRO], v[MICRO + 1], v[MICRO + 2]]);
                        let mut da: Vec<(usize, f64)> = av.iter().map(|&(j, a, k, c)| (j, c * pre.val() * dfa[a][k])).collect();
                        // A texture's values (P30): the intensity times dP/dt over P.
                        let iv = pre.val() * f2;
                        da.extend(self.texture_derivatives(ph, ip, fields, h, th.val()).into_iter().map(|(j, r)| (j, iv * r)));
                        (pre * s.map(f2, df2ds), da)
                    }
                    Intensities::List(l) => {
                        let at = l.iter().position(|e| e.0 == pk.hkl);
                        let given = at.map_or(0.0, |k| l[k].2);
                        let w = if two { weight / (v[RATIO] + 1.0) } else { weight };
                        // A Pawley fit's intensity (P36): the peak is linear in it.
                        let di = fields.iter().enumerate().filter(|(_, f)| at.is_some_and(|k| **f == Field::Intensity(ip, k))).map(|(j, _)| (j, w.val())).collect();
                        (w * given, di)
                    }
                };
                // The window and the profile, as `calculate`.
                // The asymmetry at the Bragg angle (as `peak_profile`): its test and Bérar & Baldinozzi's θ.
                let asym = bragg.val() < self.asym_limit;
                let use_fcj = asym && matches!(ph.shape, Shape::ThompsonCoxHastings | Shape::PseudoVoigt) && (ph.s_l != 0.0 || ph.d_l != 0.0);
                // The Bérar–Baldinozzi multiplier is 1 where every P is zero; it is kept when a P is asked for, so that its
                // derivative is there.
                let berar_too = asym && (ph.asym.iter().any(|&a| a != 0.0) || !qv.iter().all(|&(_, k)| k > 6));
                let (low, high) = (pk.two_theta - pk.reach, pk.two_theta + pk.reach);
                // As `calculate`: the extra low point with a Lorentzian strain, X or the quartic form's ξ (L267).
                let lorentzian_strain = ph.x != 0.0 || ph.quartic_strain.as_ref().is_some_and(|q| q.lorentzian != 0.0);
                let extra = conv.peak_window && ph.shape == Shape::ThompsonCoxHastings && lorentzian_strain;
                let start = if extra { x.partition_point(|&t| t <= low).saturating_sub(1) } else { x.partition_point(|&t| t < low) };
                let q0 = QD::var(pos.v, 0);
                // The Bragg angle, a variable of its own: it moves with the cell and λ, not with the zero or the shifts.
                let qb = QD::var(bragg.v, 9);
                let q1 = QD::var(fwhm.v, 1);
                let q2 = QD::var(eta.v, 2);
                let qa = [QD::var(ph.asym[0], 3), QD::var(ph.asym[1], 4), QD::var(ph.asym[2], 5), QD::var(ph.asym[3], 6)];
                let (qs, qd) = (QD::var(ph.s_l, 7), QD::var(ph.d_l, 8));
                let (qs, qd) = if conv.axial_divergence && ph.d_l < 0.0002 { ((qs - qd) * 0.5, (qs + qd) * 0.5) } else { (qs, qd) };
                let fcj = use_fcj.then(|| Fcj::new(q2, q1, qs, qd, q0));
                for i in start..x.len() {
                    if x[i] > high {
                        break;
                    }
                    if excluded[i] {
                        continue;
                    }
                    let dx = -q0 + x[i];
                    let omega: QD = match &fcj {
                        Some(f) if berar_too => f.at(x[i]) * berar_baldinozzi(dx / q1, qb, qa, conv.berar_baldinozzi),
                        Some(f) => f.at(x[i]),
                        None => {
                            let base = match ph.shape {
                                Shape::Gaussian => gaussian(dx, q1),
                                Shape::Lorentzian => lorentzian(dx, q1),
                                _ => pseudo_voigt(dx, q2, q1),
                            };
                            if berar_too {
                                base * berar_baldinozzi(dx / q1, qb, qa, conv.berar_baldinozzi)
                            } else {
                                base
                            }
                        }
                    };
                    // The default window's taper, a function of the position and the FWHM (reach = Wdt·H/0.8).
                    let omega = if conv.window_edge {
                        omega
                    } else if conv.peak_window {
                        omega * taper(dx, QD::cst(pk.reach))
                    } else {
                        omega * taper(dx, q1 * (self.wdt / crate::profile::TAPER_START))
                    };
                    let iv = intensity.v;
                    for &(j, k) in &pv {
                        jac[j][i] += omega.v * intensity.d[k] + iv * (omega.d[0] * pos.d[k] + omega.d[1] * fwhm.d[k] + omega.d[2] * eta.d[k] + omega.d[9] * bragg.d[k]);
                    }
                    for &(j, k) in &qv {
                        jac[j][i] += iv * omega.d[k];
                    }
                    for &(j, k) in &sv {
                        jac[j][i] += iv * (omega.d[1] * fwhm.d[SUM] + omega.d[2] * eta.d[SUM]) * monomials[k];
                    }
                    for &(j, di) in &datoms {
                        jac[j][i] += omega.v * di;
                    }
                }
            }
        }
        (c, jac)
    }

    /// The derivatives of a TOF pattern's peaks (P05), as `tof_peaks` and `peak_profile` calculate them: the position,
    /// rise, decay and widths of `tof::Parameters::at`, the intensity scale·m·F²·d⁴ sin θ_B·A·E·e^{−2B s²}·PO, the
    /// back-to-back profile (`tof::back_to_back`, with its own derivatives) and the taper.
    fn tof_derivatives(&self, bank: &crate::tof::Bank, x: &[f64], fields: &[Field], c: &Calculated, jac: &mut [Vec<f64>]) {
        use crate::tof::{back_to_back_g, extinction, gauss_fwhm, Parameters};
        let excluded = self.excluded_points(x);
        let th = (bank.two_theta / 2.0).to_radians();
        let (sin_th, s2) = (th.sin(), th.sin().powi(2));
        for (ip, ph) in self.phases.iter().enumerate() {
            let Some(prof) = &ph.tof else { continue };
            let tv: Vec<(usize, usize)> = fields.iter().enumerate().filter_map(|(j, f)| f.tof_var(ip).map(|k| (j, k))).collect();
            // (column, atom, atom parameter, coefficient), as for a constant-wavelength pattern.
            let mut av: Vec<(usize, usize, usize, f64)> = fields
                .iter()
                .enumerate()
                .flat_map(|(j, f)| match *f {
                    Field::AtomX(p, a, k) if p == ip => vec![(j, a, k, 1.0)],
                    Field::AtomB(p, a) if p == ip => vec![(j, a, 3, 1.0)],
                    Field::AtomOcc(p, a) if p == ip => vec![(j, a, 4, 1.0)],
                    Field::AtomBeta(p, a, k) if p == ip => vec![(j, a, 5 + k, 1.0)],
                    Field::Mode(p, m) if p == ip => ph.modes.as_ref().map_or(Vec::new(), |ms| ms.derivative(m).into_iter().map(|(a, k, c)| (j, a, k, c)).collect()),
                    // A rigid group's value moves its atoms and its satellites'; the cell moves them too, the groups keeping their
                    // shape (P41, `rigid_body::placement_derivative`).
                    Field::Rigid(p, g, k) if p == ip => crate::rigid_body::placement_derivative(&ph.rigid, &ph.cell, Some((g, k)), None).into_iter().map(|(a, c, d)| (j, a, c, d)).collect(),
                    Field::Cell(p, k) if p == ip && !ph.rigid.is_empty() => crate::rigid_body::placement_derivative(&ph.rigid, &ph.cell, None, Some(k)).into_iter().map(|(a, c, d)| (j, a, c, d)).collect(),
                    _ => Vec::new(),
                })
                .collect();
            // A TLS group's atoms' β (P41, L367): through their positions, the cell and the TLS values.
            av.extend(tls_terms(ph, ip, fields, &av));
            // The values as `tof_peaks` takes them, the instrument's part added (a sum: the same derivatives).
            let pf = Parameters::of(self.zero, bank, prof);
            let var = |v: f64, i: usize| TD::var(v, i);
            let par = Parameters::<TD> {
                zero: var(pf.zero, T_ZERO),
                dtt: std::array::from_fn(|k| var(pf.dtt[k], T_BANK + k)),
                thermal: pf.thermal.map(|t| std::array::from_fn(|k| var(t[k], T_BANK + 3 + k))),
                sigma: std::array::from_fn(|k| var(pf.sigma[k], T_PROFILE + k)),
                sigma_q: var(pf.sigma_q, T_PROFILE + 3),
                gamma: std::array::from_fn(|k| var(pf.gamma[k], T_PROFILE + 4 + k)),
                alpha: std::array::from_fn(|k| var(pf.alpha[k], T_PROFILE + 7 + k)),
                beta: std::array::from_fn(|k| var(pf.beta[k], T_PROFILE + 10 + k)),
                thermal_decay: pf.thermal_decay.map(|t| std::array::from_fn(|k| var(t[k], T_PROFILE + 13 + k))),
            };
            let cl: [TD; 6] = std::array::from_fn(|k| var(cell_params(ph)[k], T_CELL + k));
            let gs = reciprocal_metric(&cl);
            let g = direct_metric(&cl);
            let volume = (g[0][0] * (g[1][1] * g[2][2] - g[1][2] * g[2][1]) - g[0][1] * (g[1][0] * g[2][2] - g[1][2] * g[2][0]) + g[0][2] * (g[1][0] * g[2][1] - g[1][1] * g[2][0])).sqrt();
            let (scale, bov) = (var(ph.scale, T_SCALE), var(ph.bov, T_BOV));
            let (g1, g2) = (var(ph.pref[0], T_G), var(ph.pref[1], T_G + 1));
            let (ext, abs) = (var(prof.extinction, T_PROFILE + 17), var(prof.absorption, T_PROFILE + 18));
            for pk in c.peaks.iter().filter(|p| p.phase == ip) {
                let h = pk.hkl.map(|q| q as f64);
                let d = TD::cst(1.0) / quad(&gs, h, h).sqrt();
                let shape = par.at(d);
                let (fwhm, eta) = tch(gauss_fwhm(shape.sigma2), shape.lorentz);
                let (intensity, datoms): (TD, Vec<(usize, f64)>) = match &ph.intensities {
                    Intensities::Structure(st, sc) => {
                        let s = TD::cst(0.5) / d;
                        let (f2, df2ds, dfa) = st.intensity_derivatives(sc, pk.hkl, s.val());
                        // F² with its cell derivatives through s, and as a variable of its own for the atoms.
                        let f2d = s.map(f2, df2ds) + (TD::var(f2, T_F2) - f2);
                        let lambda = d * (2.0 * sin_th);
                        let a = bank.absorption.at(abs * lambda, s2);
                        let e = extinction(ext, f2d, lambda, volume, s2);
                        let po = self.pref_g(ph, &gs, g1, g2, h, std::f64::consts::FRAC_PI_2);
                        let i = scale * pk.multiplicity as f64 * f2d * d.powi(4) * sin_th * a * e * (-bov * s * s * 2.0).exp() * po;
                        let mut da: Vec<(usize, f64)> = av.iter().map(|&(j, a, k, c)| (j, c * i.d[T_F2] * dfa[a][k])).collect();
                        da.extend(self.texture_derivatives(ph, ip, fields, h, std::f64::consts::FRAC_PI_2).into_iter().map(|(j, r)| (j, i.val() * r)));
                        (i, da)
                    }
                    Intensities::List(l) => {
                        let at = l.iter().position(|e| e.0 == pk.hkl);
                        let di = fields.iter().enumerate().filter(|(_, f)| at.is_some_and(|k| **f == Field::Intensity(ip, k))).map(|(j, _)| (j, 1.0)).collect();
                        (TD::cst(at.map_or(0.0, |k| l[k].2)), di)
                    }
                };
                let reach = (!self.convention.window_edge).then(|| fwhm * (self.wdt / crate::profile::TAPER_START));
                let (low, high) = (pk.two_theta - pk.reach, pk.two_theta + pk.reach);
                for i in x.partition_point(|&t| t < low)..x.len() {
                    if x[i] > high {
                        break;
                    }
                    if excluded[i] {
                        continue;
                    }
                    let dt = -shape.position + x[i];
                    let omega = back_to_back_g(dt, shape.alpha, shape.beta, fwhm, eta);
                    let omega = match reach {
                        Some(r) => omega * taper(dt, r),
                        None => omega,
                    };
                    let y = intensity * omega;
                    for &(j, k) in &tv {
                        jac[j][i] += y.d[k];
                    }
                    for &(j, di) in &datoms {
                        jac[j][i] += omega.val() * di;
                    }
                }
            }
        }
    }

    /// The pattern recomputed through the dual-number path (values only): equal to `calculate`'s Ycal, which the test
    /// `engine_derivatives` checks.
    pub fn calculate_dual_values(&self, x: &[f64]) -> Vec<f64> {
        let fields: Vec<Field> = (0..self.phases.len()).map(Field::Scale).collect();
        let (c, jac) = self.derivatives(x, &fields);
        // Ycal = background + Σ_phases scale·∂Ycal/∂scale.
        let mut y = c.background.clone();
        for (p, col) in jac.iter().enumerate() {
            for (yi, d) in y.iter_mut().zip(col) {
                *yi += self.phases[p].scale * d;
            }
        }
        y
    }
}

/// Place a phase's rigid groups and write their atoms' positions into its structure (P41).
pub(crate) fn place_rigid(ph: &mut crate::powder::Phase) {
    if ph.rigid.is_empty() {
        return;
    }
    if let Intensities::Structure(st, _) = &mut ph.intensities {
        let mut x: Vec<[f64; 3]> = st.atoms.iter().map(|a| a.x).collect();
        crate::rigid_body::place(&mut ph.rigid, &ph.cell, &mut x);
        // A TLS group's atoms take its displacements (P41, L367).
        for (a, b) in crate::rigid_body::betas(&ph.rigid, &ph.cell, &x) {
            st.atoms[a].beta = Some(b);
        }
        for (a, x) in st.atoms.iter_mut().zip(x) {
            a.x = x;
        }
    }
}

/// The TLS terms of the atoms' parameters (column, atom, β component 5 + b, coefficient) for phase ip (P41, L367):
/// each position term of a TLS atom carried on to its β (∂β/∂x), each cell parameter's own effect on the β (the atoms'
/// fractions fixed; their move is already a position term), each TLS value (β linear in it) and a TLS group's origin
/// (its centre moves with it, the atoms' terms coming through their positions).
fn tls_terms(ph: &crate::powder::Phase, ip: usize, fields: &[Field], av: &[(usize, usize, usize, f64)]) -> Vec<(usize, usize, usize, f64)> {
    use crate::rigid_body::{beta_derivative, Wrt};
    if !ph.rigid.iter().any(|g| g.displacements.is_some()) {
        return Vec::new();
    }
    let Intensities::Structure(st, _) = &ph.intensities else { return Vec::new() };
    let x: Vec<[f64; 3]> = st.atoms.iter().map(|a| a.x).collect();
    let (groups, cell) = (&ph.rigid, &ph.cell);
    let tls_atoms: std::collections::BTreeSet<usize> = crate::rigid_body::betas(groups, cell, &x).into_iter().map(|(a, _)| a).collect();
    // (atom, coordinate) → the β terms its move gives (atom, component, ∂β/∂x).
    type Terms = Vec<(usize, usize, f64)>;
    let mut dx: std::collections::BTreeMap<(usize, usize), Terms> = std::collections::BTreeMap::new();
    let mut out = Vec::new();
    for &(j, a, c, d) in av.iter().filter(|e| e.2 < 3 && tls_atoms.contains(&e.1)) {
        let terms = dx.entry((a, c)).or_insert_with(|| beta_derivative(groups, cell, &x, Wrt::Position(a, c)));
        out.extend(terms.iter().map(|&(b, k, v)| (j, b, 5 + k, d * v)));
    }
    for (j, f) in fields.iter().enumerate() {
        match *f {
            Field::Cell(p, k) if p == ip => out.extend(beta_derivative(groups, cell, &x, Wrt::Cell(k)).into_iter().map(|(b, k, v)| (j, b, 5 + k, v))),
            Field::Rigid(p, g, k) if p == ip => {
                if let Some(i) = crate::rigid_body::tls_index(&groups[g], k) {
                    out.extend(beta_derivative(groups, cell, &x, Wrt::Tls(g, i)).into_iter().map(|(b, k, v)| (j, b, 5 + k, v)));
                } else if k < 3 && groups[g].displacements.is_some() {
                    out.extend(beta_derivative(groups, cell, &x, Wrt::Centre(g, k)).into_iter().map(|(b, k, v)| (j, b, 5 + k, v)));
                }
            }
            _ => {}
        }
    }
    out
}

