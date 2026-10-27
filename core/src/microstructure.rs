//! Isotropic size and strain: the integral breadths, apparent sizes and maximum strains FullProf writes, for each
//! reflection, in its microstructure file (`.mic`) when a resolution file is given.
//!
//! The widths are those of the pattern (`Model::width`, at the Bragg angle without Zero): the instrument's H_G², H_L
//! from the resolution file, and the phase's (Rodríguez-Carvajal, "Study of micro-structural effects by powder
//! diffraction using the program FullProf", `manuals/fullprof/docs/Microstructural_effects.pdf`, eq. 19, p. 4)
//!   H_G² = U tan²θ + V tan θ + W + I_G/cos²θ,   H_L = X tan θ + Y/cos θ.
//! The `.mic` file states how FullProf splits them (its header, e.g. `runs/YVO3/YVO3a_calc1.mic`, lines 31–41):
//! - size: H_G = √I_G / cos θ and H_L = Y / cos θ (LorSiz and the size models add to Y; not handled here);
//! - strain: H_G = √U tan θ and H_L = X tan θ;
//! - each pair is turned into (H, η) by the TCH formulas (eqs. 9–10, p. 2) and into the integral breadth of the
//!   pseudo-Voigt, β = ½H / (η/π + (1 − η)/√(π/ln 2)) (eq. 7), in degrees 2θ, then into reciprocal units
//!   β* = β (radians) cos θ / λ1 (Å⁻¹);
//! - apparent size D = 1/β*_size (Scherrer, eq. 14); maximum strain e = ½ β*_strain d (eq. 17), printed × 10⁴.
//!
//! Found on FullProf's runs, not stated in the papers (L039), and followed only with `Convention::microstructure` (L046, L050):
//! - β* is taken with λ1 for the λ2 reflections too, so their strains are larger by λ2/λ1 (YVO3a_calc1.mic: 61.0327
//!   against 60.8810, λ2/λ1 = 1.00249); the default takes the reflection's own λ, as eq. 14 has it;
//! - a negative U enters the strain as √|U| (YVO3a_calc.pcr after its cycle, U = −0.091446: e = 60.8810, as printed);
//!   the default leaves the strain undefined (`None`), as √U tan θ is;
//! - where the resolution file's H_G² is negative, HGi and betaG are undefined and Beta_TCH is 0, as FullProf prints;
//! - FullProf prints the apparent size and strain of the parameters after its refinement cycle, and the widths and
//!   breadths of those it read (`fullprof_microstructure`).

use crate::powder::{Model, Shape};
use crate::profile::tch;
use crate::spacegroup::SpaceGroup;
use std::f64::consts::{LN_2, PI};

/// Anisotropic microstrain as a quartic form in hkl (P25, L211): Stephens' notation, which FullProf's manual quotes
/// (p. 22, Figure 1; p. 121, LINE 37):
///   σ²(M_hkl) = 10⁻⁸ Σ_{H+K+L=4} S_HKL h^H k^K l^L   (Å⁻⁴),   M_hkl = 1/d²,
/// the variance of M over the crystallites (the fluctuations and correlations of the metric). The strain's FWHM is
/// (manual p. 22)
///   D_ST = √(8 ln 2) (180/π) σ(M)/M   (degrees, times tan θ),
/// shared between the TCH profile's Gaussian and Lorentzian parts by ξ, the Lorentzian share (the `.mic` header,
/// `runs/NdNi-LTOa/NdNi_LTOa_calc1.mic`, lines 36–38):
///   H_G² += [(1 − ξ) D_ST]² tan²θ,   H_L += ξ D_ST tan θ.
/// Each S_HKL multiplies the invariant of its monomial h^H k^K l^L under the Laue class: the sum of the monomial's
/// images under the class's rotations, scaled so that the monomial itself has coefficient 1 (`quartic_forms`): one
/// monomial for −1, 2/m and mmm; h⁴ + k⁴ for S_400 of 4/mmm, (h² + k²)l² for its S_202 (L267). Stephens' paper (J.
/// Appl. Cryst. 32, 281, 1999) is not in the repository: the form, its 10⁻⁸ and the split by ξ are as the manual and
/// the `.mic` header state them, checked on FullProf's runs of Nd₂NiO₄ (mmm), CuF₂ (2/m), Pr₀.₅Ba₀.₅CoO₃ and La₂NiO₄
/// (4/mmm) (tests `fullprof_anisotropic_strain`, `fullprof_anisotropic_size`).
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct QuarticStrain {
    /// The terms (H, K, L) with H + K + L = 4 and their S_HKL (10⁻⁸ Å⁻⁴), in the order of the phase's Laue class
    /// (`quartic_terms`).
    pub terms: Vec<([u8; 3], f64)>,
    /// Each term's invariant polynomial (monomials and coefficients, `quartic_forms`); empty for a term that is its
    /// monomial alone (the classes −1, 2/m, mmm, and a form built before L267).
    #[cfg_attr(feature = "serde", serde(default))]
    pub forms: Vec<Vec<([u8; 3], f64)>>,
    /// ξ, the Lorentzian share of the strain's width (FullProf's "Lorentzian strain coeff.", XI, LINE 38).
    pub lorentzian: f64,
}

/// Anisotropic Lorentzian size by the crystallites' shape, or by a selective rule (P25, L267): FullProf's LorSiz with
/// Size-Model ±1 or −14 (manual pp. 23–24, LINE 27 and LINE 36).
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum SizeShape {
    /// Platelets of normal t (a direction [u v w] of the crystal): g = |cos φ| (Size-Model 1).
    Platelet([f64; 3]),
    /// Needles along t: g = sin φ (Size-Model −1).
    Needle([f64; 3]),
    /// Antiphase domains of the Cu₃Au kind: the reflections whose h, k, l are not all of one parity, g = 1, the others
    /// 0 (Size-Model −14, as FullProf's own run states it: `runs/Pd3MnD/pdfinal_calc.out`, line 344).
    Antiphase,
}

/// A phase's anisotropic size: its shape and LorSiz, the Lorentzian width coefficient.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ShapeSize {
    pub shape: SizeShape,
    pub value: f64,
}

/// k_s/λ of the manual's size models: 360/π² · 10⁻³ (degrees 2θ per Å of wavelength; p. 23).
pub const SIZE_WIDTH: f64 = 360.0 / (PI * PI) * 1e-3;

impl ShapeSize {
    /// FullProf's Size-Model of the shape.
    pub fn fullprof_model(&self) -> i32 {
        match self.shape {
            SizeShape::Platelet(_) => 1,
            SizeShape::Needle(_) => -1,
            SizeShape::Antiphase => -14,
        }
    }

    /// g(h): |cos φ| for platelets, sin φ for needles, φ the angle between the reciprocal vector h* and the direction
    /// t (cos φ = h·t / (|h*| |t|), since a*_i · a_j = δ_ij); 1 or 0 for the antiphase rule.
    pub fn factor(&self, h: [f64; 3], cell: &crate::cell::Cell) -> f64 {
        let cos = |t: &[f64; 3]| {
            let ht: f64 = (0..3).map(|i| h[i] * t[i]).sum();
            let tt = crate::mat3::quad(&cell.metric(), t);
            let hh = cell.inv_d_squared(h);
            if tt <= 0.0 || hh <= 0.0 { 0.0 } else { (ht / (hh * tt).sqrt()).clamp(-1.0, 1.0) }
        };
        match &self.shape {
            SizeShape::Platelet(t) => cos(t).abs(),
            // h along the needles: 1 − cos²φ is rounding noise there, taken as 0 (as the Jacobian does).
            SizeShape::Needle(t) => {
                let s2 = 1.0 - cos(t).powi(2);
                if s2 <= 1e-12 { 0.0 } else { s2.sqrt() }
            }
            SizeShape::Antiphase => {
                let p = h.map(|v| (v.round() as i64).rem_euclid(2));
                f64::from(!(p[0] == p[1] && p[1] == p[2]))
            }
        }
    }

    /// The Lorentzian FWHM it adds at the Bragg angle (degrees): k_s LorSiz g / cos θ, k_s = 360/π² λ 10⁻³ with the
    /// reflection's λ = 2 d sin θ (manual p. 23; checked on prba10k and Pd3MnD, test `fullprof_anisotropic_size`).
    pub fn width(&self, h: [f64; 3], cell: &crate::cell::Cell, two_theta: f64) -> f64 {
        let th = (two_theta / 2.0).to_radians();
        let lambda = 2.0 * th.sin() / cell.inv_d_squared(h).sqrt();
        SIZE_WIDTH * lambda * self.value * self.factor(h, cell) / th.cos()
    }
}

/// The constant of D_ST: √(8 ln 2)·180/π, with the manual's 10⁻⁸ of σ² as 10⁻⁴ of σ (p. 22).
pub const STRAIN_WIDTH: f64 = 2.354_820_045_030_949_3 * (180.0 / PI) * 1e-4;

/// The monomial h^H k^K l^L of a term.
pub fn monomial(term: [u8; 3], h: [f64; 3]) -> f64 {
    h[0].powi(term[0] as i32) * h[1].powi(term[1] as i32) * h[2].powi(term[2] as i32)
}

impl QuarticStrain {
    /// Σ S_HKL P_HKL(h) (in 10⁻⁸ Å⁻⁴): σ²(M) / 10⁻⁸, P_HKL the term's invariant.
    pub fn sum(&self, h: [f64; 3]) -> f64 {
        (0..self.terms.len()).map(|i| self.terms[i].1 * self.term_value(i, h)).sum()
    }

    /// P_HKL(h) of term i: its invariant polynomial, or its monomial when it has none.
    pub fn term_value(&self, i: usize, h: [f64; 3]) -> f64 {
        match self.forms.get(i).filter(|f| !f.is_empty()) {
            Some(f) => f.iter().map(|(t, c)| c * monomial(*t, h)).sum(),
            None => monomial(self.terms[i].0, h),
        }
    }

    /// A form for a group: its terms with their invariants (`quartic_forms`).
    pub fn for_group(group: &SpaceGroup, terms: Vec<([u8; 3], f64)>, lorentzian: f64) -> QuarticStrain {
        let forms = quartic_forms(group, &terms.iter().map(|t| t.0).collect::<Vec<_>>());
        QuarticStrain { terms, forms, lorentzian }
    }

    /// D_ST of reflection h whose 1/d² is `inv_d2` (degrees; the width is D_ST tan θ). A negative quadratic sum, which
    /// no variance can be, gives no broadening.
    pub fn width(&self, h: [f64; 3], inv_d2: f64) -> f64 {
        STRAIN_WIDTH * self.sum(h).max(0.0).sqrt() / inv_d2
    }

    /// (Gaussian H_G², Lorentzian H_L) that the strain adds at θ for reflection h.
    pub fn widths(&self, h: [f64; 3], inv_d2: f64, tan_theta: f64) -> (f64, f64) {
        let d = self.width(h, inv_d2) * tan_theta;
        (((1.0 - self.lorentzian) * d).powi(2), self.lorentzian * d)
    }
}

/// The 15 monomials of degree 4 in h, k, l.
const MONOMIALS: [[u8; 3]; 15] = [[4, 0, 0], [0, 4, 0], [0, 0, 4], [2, 2, 0], [2, 0, 2], [0, 2, 2], [2, 1, 1], [1, 2, 1], [1, 1, 2], [3, 1, 0], [3, 0, 1], [1, 3, 0], [1, 0, 3], [0, 1, 3], [0, 3, 1]];

/// The quartic polynomial of a monomial's images under a group's rotations: Σ_R (hR)^H (kR)^K (lR)^L (reciprocal indices
/// transform as h' = h·R), as (monomial, coefficient) over the 15 monomials.
fn orbit_sum(group: &SpaceGroup, term: [u8; 3]) -> [f64; 15] {
    let mut rotations: Vec<[[i32; 3]; 3]> = group.ops.iter().map(|o| o.r).collect();
    rotations.sort();
    rotations.dedup();
    let mut sum = [0.0; 15];
    for r in &rotations {
        // The linear forms (hR)_j = Σ_i h_i R_ij, multiplied out.
        let mut poly: Vec<([u8; 3], f64)> = vec![([0, 0, 0], 1.0)];
        for (j, &e) in term.iter().enumerate() {
            for _ in 0..e {
                let mut next: Vec<([u8; 3], f64)> = Vec::new();
                for (m, c) in &poly {
                    for i in 0..3 {
                        if r[i][j] != 0 {
                            let mut m2 = *m;
                            m2[i] += 1;
                            match next.iter_mut().find(|x| x.0 == m2) {
                                Some(x) => x.1 += c * r[i][j] as f64,
                                None => next.push((m2, c * r[i][j] as f64)),
                            }
                        }
                    }
                }
                poly = next;
            }
        }
        for (m, c) in poly {
            let k = MONOMIALS.iter().position(|x| *x == m).unwrap();
            sum[k] += c;
        }
    }
    sum
}

/// The invariant of each term under the group: its orbit sum scaled so that the term's monomial has coefficient 1.
pub fn quartic_forms(group: &SpaceGroup, terms: &[[u8; 3]]) -> Vec<Vec<([u8; 3], f64)>> {
    terms
        .iter()
        .map(|t| {
            let s = orbit_sum(group, *t);
            let k = MONOMIALS.iter().position(|x| x == t).unwrap();
            if s[k].abs() < 1e-9 {
                return vec![(*t, 1.0)];
            }
            MONOMIALS.iter().zip(s).filter(|(_, c)| c.abs() > 1e-9).map(|(m, c)| (*m, c / s[k])).collect()
        })
        .collect()
}

/// FullProf's order of the terms for a Laue class, where it has a Strain-Model (manual pp. 119–120), followed by the
/// monomials in order.
fn preferred_order(laue: &str) -> Vec<[u8; 3]> {
    let first: &[[u8; 3]] = match laue {
        "4/m" | "4/mmm" => &[[4, 0, 0], [0, 0, 4], [2, 2, 0], [2, 0, 2]],
        "-3" | "-3m" | "6/m" | "6/mmm" => &[[4, 0, 0], [0, 0, 4], [1, 1, 2], [2, 1, 1]],
        "m-3" | "m-3m" => &[[4, 0, 0], [2, 2, 0]],
        _ => &[],
    };
    first.iter().chain(MONOMIALS.iter().filter(|m| !first.contains(m))).copied().collect()
}

/// The terms of the quartic form a space group allows: the monomials whose invariants under its rotations are not
/// zero and are independent of those before them, in FullProf's order for its Laue class (manual pp. 119–120), the
/// others after. For −1, 2/m and mmm each invariant is one monomial; for 4/mmm they are S_400 (h⁴ + k⁴), S_004 (l⁴),
/// S_220 (h²k²) and S_202 ((h² + k²)l²), and 4/m adds S_310 (h³k − hk³), which FullProf's Strain-Model 4 does not have
/// (L267). The trigonal, hexagonal and cubic invariants follow from the same rule, unchecked against FullProf.
pub fn quartic_terms(group: &SpaceGroup) -> Result<Vec<[u8; 3]>, String> {
    let mut basis: Vec<[f64; 15]> = Vec::new();
    let mut out = Vec::new();
    for t in preferred_order(group.laue) {
        let mut v = orbit_sum(group, t);
        let k = MONOMIALS.iter().position(|x| *x == t).unwrap();
        if v[k].abs() < 1e-9 {
            continue;
        }
        for b in &basis {
            let d: f64 = (0..15).map(|i| v[i] * b[i]).sum();
            for i in 0..15 {
                v[i] -= d * b[i];
            }
        }
        let n = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        if n > 1e-9 {
            basis.push(v.map(|x| x / n));
            out.push(t);
        }
    }
    // The diagonal classes in FullProf's order of the monomials (as before L267).
    if group.ops.iter().all(|o| (0..3).all(|i| (0..3).all(|j| i == j || o.r[i][j] == 0))) {
        out.sort_by_key(|t| MONOMIALS.iter().position(|m| m == t));
    }
    Ok(out)
}

/// FullProf's Strain-Model for a Laue class of diagonal rotations, and the order of its S_HKL in the `.pcr` (manual
/// pp. 119–120, LINE 37): 2 for 1 2/m 1, −2 for 1 1 2/m, 3 for mmm. Model 1 (−1) is not read: the manual lists S_301
/// twice and no S_310 (p. 119), so the order of its 15 terms is not known from the repository.
pub fn fullprof_strain_model(model: i32) -> Option<&'static [[u8; 3]]> {
    match model {
        2 => Some(&[[4, 0, 0], [0, 4, 0], [0, 0, 4], [2, 2, 0], [2, 0, 2], [0, 2, 2], [1, 2, 1], [3, 0, 1], [1, 0, 3]]),
        -2 => Some(&[[4, 0, 0], [0, 4, 0], [0, 0, 4], [2, 2, 0], [2, 0, 2], [0, 2, 2], [1, 1, 2], [3, 1, 0], [1, 3, 0]]),
        3 => Some(&[[4, 0, 0], [0, 4, 0], [0, 0, 4], [2, 2, 0], [2, 0, 2], [0, 2, 2]]),
        // 4/m and 4/mmm (L267): checked on prba10k and lac.
        4 | 5 => Some(&[[4, 0, 0], [0, 0, 4], [2, 2, 0], [2, 0, 2]]),
        _ => None,
    }
}

/// The Strain-Model FullProf reads a quartic form's terms with (2, −2, 3, or 5 for the tetragonal classes), from the
/// terms a group allows; None for −1, whose order in the `.pcr` is not known (`fullprof_strain_model`), and for the
/// classes not checked against FullProf. 4/m's S_310, which FullProf has not, is left out: the export refuses it
/// when it is not zero.
pub fn fullprof_model_of(terms: &[[u8; 3]]) -> Option<i32> {
    let mut have: Vec<[u8; 3]> = terms.iter().copied().filter(|t| *t != [3, 1, 0]).collect();
    have.sort();
    [2, -2, 3, 5].into_iter().find(|&m| {
        let mut want = fullprof_strain_model(m).unwrap().to_vec();
        want.sort();
        want == have
    })
}

/// How many S_HKL FullProf reads on LINE 37 for a Strain-Model (manual pp. 119–120), to read past them.
pub fn fullprof_strain_terms(model: i32) -> Option<usize> {
    match model {
        1 => Some(15),
        2 | -2 => Some(9),
        3 => Some(6),
        4..=7 => Some(4),
        8..=12 => Some(3),
        13 | 14 => Some(2),
        _ => None,
    }
}

/// The term's name, "S_400".
pub fn term_name(t: [u8; 3]) -> String {
    format!("S_{}{}{}", t[0], t[1], t[2])
}

/// The integral breadth of the pseudo-Voigt of FWHM h and mixing η (eq. 7), in the units of h.
pub fn pseudo_voigt_breadth(h: f64, eta: f64) -> f64 {
    if h <= 0.0 {
        return 0.0;
    }
    0.5 * h / (eta / PI + (1.0 - eta) / (PI / LN_2).sqrt())
}

/// The integral breadth of the Voigt of Gaussian and Lorentzian FWHM hg, hl through its TCH pseudo-Voigt.
pub fn voigt_breadth(hg: f64, hl: f64) -> f64 {
    let (h, eta) = tch(hg, hl);
    pseudo_voigt_breadth(h, eta)
}

/// What the `.mic` file lists for one reflection. Widths in degrees 2θ; breadths in Å⁻¹ × 1000; size in Å; strain
/// × 10⁴.
#[derive(Clone, Debug, PartialEq)]
pub struct Breadths {
    /// The reflection, its Bragg angle (degrees, without Zero), d (Å) and d* = 1/d (Å⁻¹, the `.mic`'s s), the
    /// Williamson–Hall plot's abscissa (Explore's E19).
    pub hkl: [i32; 3],
    pub two_theta: f64,
    pub d: f64,
    pub d_star: f64,
    /// Gaussian and Lorentzian FWHM, observed (instrument and sample) and instrumental: HGo, HGi, HLo, HLi.
    pub hg_obs: f64,
    pub hg_ins: f64,
    pub hl_obs: f64,
    pub hl_ins: f64,
    /// Integral breadths of the sample's Gaussian and Lorentzian parts (betaG, betaL), and of their Voigt (Beta_TCH).
    pub beta_g: f64,
    pub beta_l: f64,
    pub beta_tch: f64,
    /// The same three as β* = β cos θ/λ in Å⁻¹, the units of `d_star`: the Williamson–Hall plot's ordinates (E19).
    pub beta_g_star: f64,
    pub beta_l_star: f64,
    pub beta_tch_star: f64,
    /// FWHM, η and integral breadth of the whole profile (HpVo, eta, betaVo).
    pub h_obs: f64,
    pub eta_obs: f64,
    pub beta_obs: f64,
    /// Apparent size (Å) and maximum strain (× 10⁴); `None` without size or strain broadening.
    pub size: Option<f64>,
    pub strain: Option<f64>,
}

impl Model {
    /// The `.mic` quantities of reflection hkl of a phase, for wavelength 1 or 2. `None` if the phase does not use the
    /// TCH profile or the reflection does not exist at that wavelength.
    pub fn breadths(&self, phase: usize, hkl: [i32; 3], wavelength: usize) -> Option<Breadths> {
        let ph = self.phases.get(phase)?;
        if ph.shape != Shape::ThompsonCoxHastings {
            return None;
        }
        let d = ph.cell.d_spacing(hkl.map(|v| v as f64))?;
        let sin_th = self.lambda[wavelength - 1] / (2.0 * d);
        if sin_th >= 1.0 {
            return None;
        }
        let th = sin_th.asin();
        let two_theta = 2.0 * th.to_degrees();
        let (t, c) = (th.tan(), th.cos());
        let (ig2, il) = self.instrumental.at(two_theta);
        // The quartic strain's Gaussian and Lorentzian parts, [(1 − ξ)D_ST]² tan²θ and ξ D_ST tan θ (P25).
        let (ag2, al) = self.strain_widths(ph, hkl, two_theta);
        let gs2 = (ph.u * t + ph.v) * t + ph.w + ph.gauss_size / (c * c) + ag2;
        let ls = ph.x * t + ph.y / c + al;
        // Degrees 2θ to Å⁻¹ × 1000, with the reflection's λ (eq. 14); FullProf uses λ1 for both wavelengths.
        let fullprof = self.convention.microstructure;
        let lambda = if fullprof { self.lambda[0] } else { self.lambda[wavelength - 1] };
        let star = |beta_deg: f64| beta_deg.to_radians() * c / lambda * 1000.0;
        let (hg_obs, hl_obs) = ((gs2 + ig2).max(0.0).sqrt(), ls + il);
        let (h_obs, eta_obs) = tch(hg_obs, hl_obs);
        // FullProf takes the sample's H_G as √(HGo² − HGi²): where the resolution file gives a negative H_G², HGi,
        // betaG are undefined (it prints NaN) and Beta_TCH is printed as 0 (NdNi_HTO_calc1.mic above 150°).
        let gs = if ig2 < 0.0 { f64::NAN } else { gs2.max(0.0).sqrt() };
        let size_b = star(voigt_breadth(ph.gauss_size.max(0.0).sqrt() / c, ph.y / c));
        // The strain's H_G = √(U tan²θ + [(1 − ξ)D_ST tan θ]²), H_L = X tan θ + ξ D_ST tan θ (the `.mic` header).
        let u = if fullprof { ph.u.abs() } else { ph.u };
        let strain_g2 = u * t * t + ag2;
        let strain_b = star(voigt_breadth(strain_g2.max(0.0).sqrt(), ph.x * t + al));
        let (beta_g, beta_l) = (star(0.5 * gs * (PI / LN_2).sqrt()), star(0.5 * PI * ls));
        let beta_tch = if gs.is_nan() { 0.0 } else { star(voigt_breadth(gs, ls)) };
        Some(Breadths {
            hkl,
            two_theta,
            d,
            d_star: 1.0 / d,
            hg_obs,
            hg_ins: if ig2 < 0.0 { f64::NAN } else { ig2.sqrt() },
            hl_obs,
            hl_ins: il,
            beta_g,
            beta_l,
            beta_tch,
            beta_g_star: beta_g / 1000.0,
            beta_l_star: beta_l / 1000.0,
            beta_tch_star: beta_tch / 1000.0,
            h_obs,
            eta_obs,
            beta_obs: star(pseudo_voigt_breadth(h_obs, eta_obs)),
            size: (size_b > 0.0).then(|| 1000.0 / size_b),
            strain: (strain_b > 0.0 && strain_g2 >= 0.0).then(|| 0.5 * strain_b / 1000.0 * d * 1e4),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn breadth_limits() {
        // A Gaussian's integral breadth is ½H√(π/ln2), a Lorentzian's πH/2 (eqs. 1–2, p. 1).
        assert!((pseudo_voigt_breadth(0.2, 0.0) - 0.1 * (PI / LN_2).sqrt()).abs() < 1e-12);
        assert!((pseudo_voigt_breadth(0.2, 1.0) - 0.1 * PI).abs() < 1e-12);
        assert!((voigt_breadth(0.0, 0.2) - 0.1 * PI).abs() < 1e-4, "TCH η is 1.00000 for a pure Lorentzian");
    }
}
