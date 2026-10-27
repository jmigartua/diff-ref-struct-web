//! Peak-shape functions for constant-wavelength powder patterns, normalised to unit area in the variable x = 2θ − 2θ_h
//! (degrees), with H the full width at half maximum.
//!
//! - Gaussian, Lorentzian, pseudo-Voigt: FullProf manual eqs. 3.8, 3.9, 3.13 (pp. 47–48). The Gaussian is set to zero
//!   beyond |x| > 6H, as CrysFML's `Gaussian` does (`CFML_PowderProfiles_Finger.f90`).
//! - Thompson–Cox–Hastings: H and η from the Gaussian and Lorentzian widths, eqs. 3.16–3.17 (p. 48).
//! - Finger–Cox–Jephcoat axial divergence: the pseudo-Voigt convoluted with the axial-divergence function of Finger,
//!   Cox & Jephcoat (J. Appl. Cryst. 27, 892, 1994), ported from CrysFML's `Prof_Val` (commit 8ed8670, the version of
//!   May 2013 following J. Hester, J. Appl. Cryst. 46, 1219, 2013), value only, with its Gauss–Legendre tables
//!   (`fcj_table`).
//! - Bérar–Baldinozzi empirical asymmetry: the multiplier of eq. 3.29 (p. 54), as FullProf computes it, which is not
//!   as the manual prints it (L035): z of the opposite sign, and F_b = (4z² − 3)·F_a where the manual has
//!   2(2z² − 3)·F_a.

use crate::dual::Real;
use crate::fcj_table::{WP, XP};
use std::f64::consts::{LN_2, PI};

/// Normalised Gaussian of FWHM h at x (zero beyond |x| > 6h). Generic over `Real` (f64, or dual numbers for the
/// engine's derivatives), as are the other shapes.
pub fn gaussian<T: Real>(x: T, h: T) -> T {
    if x.val().abs() / h.val() > 6.0 {
        return T::cst(0.0);
    }
    let c = 2.0 * LN_2.sqrt();
    let cg = 2.0 * (LN_2 / PI).sqrt();
    (-(x * c / h).powi(2)).exp() * cg / h
}

/// Normalised Lorentzian of FWHM h at x.
pub fn lorentzian<T: Real>(x: T, h: T) -> T {
    h * (2.0 / PI) / (x * x * 4.0 + h * h)
}

/// Normalised pseudo-Voigt η L + (1 − η) G of FWHM h.
pub fn pseudo_voigt<T: Real>(x: T, eta: T, h: T) -> T {
    eta * lorentzian(x, h) + (-eta + 1.0) * gaussian(x, h)
}

/// FWHM and η of the pseudo-Voigt that approximates the Voigt of Gaussian and Lorentzian FWHM hg, hl (TCH, eqs.
/// 3.16–3.17; Thompson, Cox & Hastings 1987, eqs. 2–3, p. 80).
pub fn tch<T: Real>(hg: T, hl: T) -> (T, T) {
    let h5 = hg.powi(5) + hg.powi(4) * hl * 2.69269 + hg.powi(3) * hl.powi(2) * 2.42843 + hg.powi(2) * hl.powi(3) * 4.47163 + hg * hl.powi(4) * 0.07842 + hl.powi(5);
    let h = if h5.val() > 0.0 { h5.powf(0.2) } else { T::cst(0.0) };
    let r = if h.val() > 0.0 { hl / h } else { T::cst(0.0) };
    (h, r * 1.36603 - r * r * 0.47719 + r * r * r * 0.11116)
}

/// Bérar–Baldinozzi asymmetry multiplier at z = (2θ − 2θ_h)/H, θ_h the Bragg angle (degrees for 2θ_h).
///
/// Bérar & Baldinozzi (1993), J. Appl. Cryst. 26, 128–129 (`manuals/papers/Berar_Baldinozzi_1993.pdf`):
/// g = g₀ [1 + A(θ) F_a(u) + B(θ) F_b(u)] (eq. 14, p. 129) with u = dθ/w the deviation-to-width ratio,
/// F_a(u) = 2u exp(−u²) (eq. 12), F_b(u) = (8u³ − 12u) exp(−u²) (eq. 13) and A(θ) = A₀/tan θ + A₁/tan 2θ (eq. 15).
/// The paper's deviation dθ is positive on the low-angle side: its transfer functions are nonzero for dθ > 0 (eqs.
/// 1–3, p. 128), and axial divergence sends intensity only below the Bragg angle (Finger, Cox & Jephcoat 1994,
/// p. 895). So u = −z here, as FullProf has it; FullProf's manual prints z (eq. 3.29, p. 54), the opposite sign. The
/// width w is not named in the paper; the FWHM is taken, as FullProf does.
///
/// `fullprof = false`: the paper's F_b. `fullprof = true`: FullProf's F_b = (4u² − 3)·F_a = (8u³ − 6u) exp(−u²)
/// (found from its output, L035; the EDI project found the same, enhantica/edi 855a426). FullProf's differs from the
/// paper's by exactly 3·F_a, so FullProf's (P1, P2) are the paper's (P1 + 3·P2, P2), and the same for (P3, P4)
/// (`powder::Model::translate_from_fullprof`).
pub fn berar_baldinozzi<T: Real>(z: T, two_theta_h: T, p: [T; 4], fullprof: bool) -> T {
    let th = (two_theta_h / 2.0).to_radians();
    let u = -z;
    let e = (-u * u).exp();
    let fa = u * e * 2.0;
    let fb = if fullprof { (u * u * 4.0 - 3.0) * fa } else { (u * u * u * 8.0 - u * 12.0) * e };
    (p[0] * fa + p[1] * fb) / th.tan() + (p[2] * fa + p[3] * fb) / (th * 2.0).tan() + 1.0
}

/// Where the default window's taper begins, as a fraction of the half-window.
pub const TAPER_START: f64 = 0.8;

/// The default model's peak window (L058): 1 within 0.8 of the half-window `reach`, then a cosine down to 0 at the
/// edge, with zero slope there. The default's half-window is Wdt·H/0.8, so the profile is as before within ±Wdt·H
/// (FullProf's window, manual p. 75) and its tail is brought to zero beyond: S does not jump when a window edge crosses
/// a point as the widths and positions move (FullProf cuts the profile at the edge, `Convention::window_edge`).
pub fn taper<T: Real>(dx: T, reach: T) -> T {
    let u = dx.abs() / reach;
    if u.val() <= TAPER_START {
        T::cst(1.0)
    } else if u.val() >= 1.0 {
        T::cst(0.0)
    } else {
        (((u - TAPER_START) * (std::f64::consts::PI / (1.0 - TAPER_START))).cos() + 1.0) * 0.5
    }
}

/// FullProf's (CrysFML's) reading of its axial-divergence pair (S_L, D_L) as the paper's (S/L, H/L): when D_L is below
/// 0.0002 the pair is taken as (H/L + S/L, H/L − S/L) (`CFML_PowderProfiles_Finger.f90`, `Prof_Val`, lines 1289–1294,
/// commit 8ed8670), so S_L alone means S/L = H/L = S_L/2; otherwise S/L = S_L and H/L = D_L.
pub fn fullprof_axial(s_l: f64, d_l: f64) -> (f64, f64) {
    if d_l < 0.0002 {
        (0.5 * (s_l - d_l), 0.5 * (s_l + d_l))
    } else {
        (s_l, d_l)
    }
}

const NTERMS: [usize; 14] = [6, 10, 20, 40, 60, 80, 100, 150, 200, 300, 400, 600, 800, 1000];
const FSTTERM: [usize; 14] = [0, 3, 8, 18, 38, 68, 108, 158, 233, 333, 483, 683, 983, 1383];

/// The pseudo-Voigt (η, FWHM h) of a peak at two_theta0 convoluted with the Finger–Cox–Jephcoat axial divergence,
/// with s_l = S/L and h_l = H/L the sample and detector half-heights over the goniometer radius (Finger, Cox &
/// Jephcoat 1994, J. Appl. Cryst. 27, 892–900, `manuals/papers/Finger_Cox_Jephcoat_1994.pdf`: 2φ_min, eq. 4, and
/// 2φ_infl, eq. 5, p. 895; the Gauss–Legendre sum of eq. 9, p. 895), evaluated at two_theta. Port of CrysFML's
/// `Prof_Val` (value only), with S/L and H/L taken as given; FullProf's reading of its pair is `fullprof_axial`.
pub struct Fcj<T: Real = f64> {
    two_theta0: T,
    s_l: T,
    d_l: T,
    emin: T,
    einfl: T,
    cstwoth: T,
    ngt: usize,
    it: isize,
    eta: T,
    h: T,
    plain: bool,
}

impl<T: Real> Fcj<T> {
    pub fn new(eta: T, h: T, s_l: T, h_l: T, two_theta0: T) -> Fcj<T> {
        let z = T::cst(0.0);
        let use_asym = s_l.val().abs() + h_l.val().abs() >= 0.00001;
        let mut f = Fcj { two_theta0, s_l: z, d_l: z, emin: z, einfl: z, cstwoth: z, ngt: 0, it: 0, eta, h, plain: !use_asym || (two_theta0.val() - 90.0).abs() < 0.4 };
        if f.plain {
            return f;
        }
        let eps_asym2 = 0.0002;
        let eps = 1.0e-6;
        let twoth0r = two_theta0.to_radians();
        let cstwoth = twoth0r.cos();
        let d_l = h_l;
        let apb = s_l + d_l;
        let amb = s_l - d_l;
        let s_eq_d = amb.val().abs() <= eps_asym2;
        let apb2 = apb * apb;
        let tmp = (amb * amb + 1.0).sqrt() * cstwoth;
        let mut einfl = if tmp.val().abs() > 1.0 || tmp.val().abs() <= cstwoth.val().abs() { two_theta0 } else { tmp.acos().to_degrees() };
        let tmp = (apb2 + 1.0).sqrt() * cstwoth;
        if s_l.val().abs() <= eps || d_l.val().abs() <= eps || s_eq_d {
            einfl = two_theta0;
        }
        let emin = if tmp.val().abs() <= 1.0 {
            tmp.acos().to_degrees()
        } else if tmp.val() > 0.0 {
            T::cst(0.0)
        } else {
            T::cst(180.0)
        };
        let k = (1000.0 * (two_theta0.val() - emin.val())) as i64;
        let mut arraynum = 0;
        while arraynum < 13 && k > NTERMS[arraynum] as i64 {
            arraynum += 1;
        }
        let ngt = NTERMS[arraynum];
        f.it = FSTTERM[arraynum] as isize - (ngt / 2) as isize;
        (f.s_l, f.d_l, f.emin, f.einfl, f.cstwoth, f.ngt) = (s_l, d_l, emin, einfl, cstwoth, ngt);
        f
    }

    /// The profile at two_theta.
    pub fn at(&self, two_theta: f64) -> T {
        let x = -self.two_theta0 + two_theta;
        if self.plain {
            return pseudo_voigt(x, self.eta, self.h);
        }
        let (t0, emin) = (self.two_theta0, self.emin);
        let apb = self.s_l + self.d_l;
        let mut sumwg = T::cst(0.0);
        let mut sumwrg = T::cst(0.0);
        for k in self.ngt / 2 + 1..=self.ngt {
            let idx = (k as isize + self.it) as usize;
            for side in [1.0, -1.0] {
                let delta = (t0 + emin) / 2.0 + (t0 - emin) * (side * XP[idx] / 2.0);
                let mut cosdelta = delta.to_radians().cos();
                if cosdelta.val().abs() < 1.0e-15 {
                    cosdelta = T::cst(1.0e-15);
                }
                let rcosdelta = T::cst(1.0) / cosdelta;
                let tmp = cosdelta * cosdelta - self.cstwoth * self.cstwoth;
                let f = if tmp.val() > 0.00001 { self.cstwoth.abs() / tmp.sqrt() } else { T::cst(0.0) };
                let g = if (delta - emin).val().abs() > (self.einfl - emin).val().abs() {
                    (if self.s_l.val() > self.d_l.val() { self.d_l } else { self.s_l }) * f * rcosdelta * 2.0
                } else {
                    ((apb * f - 1.0) * rcosdelta).abs()
                };
                let pv = pseudo_voigt(-delta + two_theta, self.eta, self.h);
                sumwg = sumwg + g * WP[idx];
                sumwrg = sumwrg + g * pv * WP[idx];
            }
        }
        if sumwg.val().abs() <= 1.0e-6 {
            T::cst(0.0)
        } else {
            sumwrg / sumwg
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn area(f: impl Fn(f64) -> f64, from: f64, to: f64) -> f64 {
        let n = 100000;
        let dx = (to - from) / n as f64;
        (0..n).map(|i| f(from + (i as f64 + 0.5) * dx)).sum::<f64>() * dx
    }

    #[test]
    fn normalised_shapes() {
        assert!((area(|x| gaussian(x, 0.3), -3.0, 3.0) - 1.0).abs() < 1e-9);
        assert!((pseudo_voigt(0.15, 0.4, 0.3) / pseudo_voigt(0.0, 0.4, 0.3) - 0.5).abs() < 1e-12, "half maximum at H/2");
        // The Lorentzian's tails are long: its area within ±1000 H is 1 − 2/(π·2000).
        assert!((area(|x| lorentzian(x, 0.3), -300.0, 300.0) - (1.0 - 1.0 / (PI * 1000.0))).abs() < 1e-6);
    }

    #[test]
    fn tch_limits() {
        let (h, eta) = tch(0.2, 0.0);
        assert!((h - 0.2).abs() < 1e-12 && eta.abs() < 1e-12);
        let (h, eta) = tch(0.0, 0.2);
        assert!((h - 0.2).abs() < 1e-12 && (eta - (1.36603 - 0.47719 + 0.11116)).abs() < 1e-12);
    }

    #[test]
    fn fcj_keeps_area_and_shifts_low_angle_peaks_down() {
        let f = Fcj::new(0.3, 0.2, 0.03, 0.02, 15.0);
        let a = area(|t| f.at(t), 10.0, 20.0);
        let plain = area(|t| pseudo_voigt(t - 15.0, 0.3, 0.2), 10.0, 20.0); // the same window (Lorentzian tails cut)
        assert!((a - plain).abs() < 2e-4, "area {a} against {plain}");
        let mean = area(|t| t * f.at(t), 10.0, 20.0) / a;
        assert!(mean < 15.0, "axial divergence moves low-angle peaks to lower angles: {mean}");
    }
}
