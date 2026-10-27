//! Time-of-flight neutron patterns (P05): one detector bank, the time of flight of a reflection from its d-spacing,
//! the peak as back-to-back exponentials convolved with a pseudo-Voigt, and the bank's Lorentz factor, absorption and
//! extinction.
//!
//! The equations are the papers' and FullProf's documents', each checked against FullProf's runs of
//! `manuals/fullprof/runs/{arg_si,CeO2_PEARL,alum,lamn_pol,nac_osiris,gabriela}` (L334, L341), reflection by reflection
//! on their `.out` lists and point by point on their `.prf`:
//! - the position T = Zero + Dtt1 d + Dtt2 d² + Dtt_1overd/d (FullProf_Manual.pdf p. 82; its `.out`: "TOF=TOF+
//!   dtt_1overD/dsp"), or Hodges' epithermal and thermal positions mixed by n = ½ erfc(w (x_c − 1/d)) (p. 84–85): Tₑ =
//!   Zero + Dtt1 d, Tₜ = Zerot + Dtt1t d − Dtt2t/d, T = n Tₑ + (1 − n) Tₜ, to 0.003 µs on gabriela's 1288 reflections;
//! - the rise and decay α = α₀ + α₁/d + α_Q/√d, β = β₀ + β₁/d⁴ + β_Q/d² (p. 112 and "Using FullProf to analyze Time of
//!   Flight Neutron Powder Diffraction data", `docs/TOF_FullProf.pdf`, Panel 8, for α₀ … β₁; α_Q and β_Q found from
//!   FullProf's `.out` of alum, to its six digits, and GSAS-II's beta-q is the same term,
//!   `GSASII/GSASIIpwd.py#L2484`); with Hodges' (Npr = 10), 1/α = n(α₀ + α₁ d) + (1 − n)(α₀ₜ − α₁ₜ/d) and the same for
//!   β, as the 2003 document's Panel 10 writes it in µs, where the manual (p. 112) mixes α itself: FullProf's `.out`
//!   follows the document, to six digits;
//! - the widths σ² = σ₀ + σ₁ d² + σ₂ d⁴ + σ_Q d (p. 110; σ_Q as GSAS-II's sig-q, `GSASIIpwd.py#L2510`, believed: alum's
//!   0.001 is below the print of FullProf's `.out`), H_G = √(8 ln 2 σ²), H_L = γ₀ + γ₁ d + γ₂ d², and H, η by the
//!   Thompson–Cox–Hastings relations (`profile::tch`), as the `.out`'s HwG, HwL, Hw and ETA;
//! - the profile, Jorgensen and Von Dreele's convolution as CrysFML codes it (`Tof_Jorgensen_Vondreele`,
//!   deps-pycrysfml 2258e77 `CFML/Src/CFML_Profiles/Profile_TOF_Jorg_Vondreele.f90`; the 2003 document's Panels 1–5):
//!   Ω = (1 − η) N [eᵘ erfc y + eᵛ erfc z] − η N (2/π) [Im eᵖ E₁(p) + Im e^q E₁(q)], N = αβ/(2(α + β));
//! - the intensity scale·m·F²·d⁴ sin θ_B·A·E, θ_B half the bank's 2θ: d⁴ sin θ_B is the Lorentz factor (FullProf's
//!   overall factor against the core's, with d⁴ sin θ_B: 0.9531/0.9533 for arg_si's 144.8°, 0.9997 for CeO2's 90°,
//!   0.9943 for lamn_pol's 29.0°, 1.0001 for alum without its Q terms);
//! - the absorption A(λ), λ = 2 d sin θ_B, x = A₁λ (`Absorption`) and Sabine's extinction (`extinction`).

use crate::dual::Real;

/// The parameter dictionary's ids of a bank's values, in the order of `Field::Tof`.
pub const BANK_FIELDS: [&str; 8] = ["tof_dtt1", "tof_dtt2", "tof_dtt_1overd", "tof_thermal_zero", "tof_thermal_dtt1", "tof_thermal_dtt2", "tof_crossover", "tof_crossover_width"];

/// The parameter dictionary's ids of a phase's TOF profile values, in the order of `Field::TofProfile`.
pub const PROFILE_FIELDS: [&str; 19] = [
    "tof_sigma0", "tof_sigma1", "tof_sigma2", "tof_sigma_q", "tof_gamma0", "tof_gamma1", "tof_gamma2", "tof_alpha0", "tof_alpha1", "tof_alpha_q", "tof_beta0", "tof_beta1", "tof_beta_q",
    "tof_alpha0_thermal", "tof_beta0_thermal", "tof_alpha1_thermal", "tof_beta1_thermal", "tof_extinction", "tof_absorption",
];

impl Bank {
    pub fn get(&self, k: usize) -> f64 {
        let t = self.thermal.unwrap_or(Thermal { zero: 0.0, dtt1: 0.0, dtt2: 0.0, cross: 0.0, width: 0.0 });
        [self.dtt1, self.dtt2, self.dtt_1overd, t.zero, t.dtt1, t.dtt2, t.cross, t.width][k]
    }

    pub fn set(&mut self, k: usize, v: f64) {
        match k {
            0 => self.dtt1 = v,
            1 => self.dtt2 = v,
            2 => self.dtt_1overd = v,
            _ => {
                let t = self.thermal.get_or_insert(Thermal { zero: 0.0, dtt1: 0.0, dtt2: 0.0, cross: 0.0, width: 0.0 });
                *[&mut t.zero, &mut t.dtt1, &mut t.dtt2, &mut t.cross, &mut t.width][k - 3] = v;
            }
        }
    }
}

impl Profile {
    pub fn get(&self, k: usize) -> f64 {
        let th = self.thermal.unwrap_or_default();
        match k {
            0..=2 => self.sigma[k],
            3 => self.sigma_q,
            4..=6 => self.gamma[k - 4],
            7..=9 => self.alpha[k - 7],
            10..=12 => self.beta[k - 10],
            13..=16 => th[k - 13],
            17 => self.extinction,
            _ => self.absorption,
        }
    }

    pub fn set(&mut self, k: usize, v: f64) {
        match k {
            0..=2 => self.sigma[k] = v,
            3 => self.sigma_q = v,
            4..=6 => self.gamma[k - 4] = v,
            7..=9 => self.alpha[k - 7] = v,
            10..=12 => self.beta[k - 10] = v,
            13..=16 => self.thermal.get_or_insert([0.0; 4])[k - 13] = v,
            17 => self.extinction = v,
            _ => self.absorption = v,
        }
    }
}

/// The absorption of the sample, FullProf's Iabscor (line 8, p. 76), with x = A₁ λ, A₁ the phase's ABSCOR1.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Absorption {
    /// Iabscor = 1, a flat plate perpendicular to the beam: (1 − e^{−2x})/(2x) (FullProf's CORR in `arg_si_abs1.out`
    /// to its six digits; the manual gives no formula).
    FlatPlate,
    /// Iabscor = 2, a cylinder: Hewat's exp(−(1.7133 − 0.0368 sin²θ) x + (0.0927 + 0.375 sin²θ) x²), the expression of
    /// the constant-wavelength cylinder (FullProf manual eq. 3.33) with μR = x (FullProf's CORR in `arg_si_calc.out`
    /// divided by the extinction, to 8·10⁻⁶).
    Cylinder,
    /// Iabscor = 3: exp(−x) (FullProf's CORR in `arg_si_abs3.out` to its six digits; the manual, p. 76, prints
    /// exp(−c λ²)).
    Exponential,
}

impl Absorption {
    pub fn from_iabscor(i: i32) -> Result<Absorption, String> {
        match i {
            1 => Ok(Absorption::FlatPlate),
            2 => Ok(Absorption::Cylinder),
            3 => Ok(Absorption::Exponential),
            n => Err(format!("Iabscor = {n}: FullProf's TOF absorption is 1 (flat plate), 2 (cylinder) or 3 (exponential)")),
        }
    }

    pub fn iabscor(self) -> i32 {
        match self {
            Absorption::FlatPlate => 1,
            Absorption::Cylinder => 2,
            Absorption::Exponential => 3,
        }
    }

    /// The transmission for x = A₁ λ, sin²θ_B = s2, with its derivative also at x = 0 (where an ABS1 refined from 0
    /// starts).
    pub fn at<T: Real>(self, x: T, s2: f64) -> T {
        match self {
            // Near 0 by its series, 1 − x + 2x²/3 − x³/3 (the quotient is 0/0 there).
            Absorption::FlatPlate if x.val().abs() < 1e-4 => -x + x * x * (2.0 / 3.0) - x * x * x * (1.0 / 3.0) + 1.0,
            Absorption::FlatPlate => (-(x * -2.0).exp() + 1.0) / (x * 2.0),
            Absorption::Cylinder => (x * (0.0368 * s2 - 1.7133) + x * x * (0.0927 + 0.375 * s2)).exp(),
            Absorption::Exponential => (-x).exp(),
        }
    }
}

/// Hodges' formulation of TOF(d) (FullProf's Npr = 10, line 16, pp. 84–85): the thermal position Tₜ = zero + dtt1·d −
/// dtt2/d and the crossover n = ½ erfc(width (cross − 1/d)) between it and the epithermal Zero + Dtt1·d.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Thermal {
    pub zero: f64,
    pub dtt1: f64,
    pub dtt2: f64,
    pub cross: f64,
    pub width: f64,
}

/// A TOF pattern's detector bank (P05): how d becomes a time of flight, the bank's angle, the sample's absorption and
/// the range that sets the Fourier background's period. The zero point is the model's `zero`.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Bank {
    /// T = Zero + dtt1·d + dtt2·d² + dtt_1overd/d (µs, d in Å).
    pub dtt1: f64,
    pub dtt2: f64,
    pub dtt_1overd: f64,
    /// The bank's scattering angle 2θ (degrees), in the Lorentz factor, the absorption and the extinction.
    pub two_theta: f64,
    /// Hodges' epithermal–thermal TOF(d), when the bank has it (FullProf's Npr = 10).
    pub thermal: Option<Thermal>,
    pub absorption: Absorption,
    /// The instrument's part of the profile, from a resolution file, added to the bank's and the phases' values
    /// (FullProf's Res = 5, 6): what the refinement varies is the rest.
    #[cfg_attr(feature = "serde", serde(default))]
    pub instrument: Option<Instrument>,
    /// The pattern's last time of flight (FullProf's TOF-max): the Fourier background is Σ B_j cos(jπT/T_max)
    /// (FullProf manual eq. 3.4′ with C = T_max/π, where the manual, p. 47, says C = T_max: arg_si_calc.prf's
    /// background is met to its print, 0.03, so and missed by 49 with C = T_max).
    pub tmax: f64,
}

/// The instrument's part of a TOF profile: a resolution file's values, added to the `.pcr`'s ("Using FullProf to analyze
/// Time of Flight Neutron Powder Diffraction data", `docs/TOF_FullProf.pdf`, "External instrumental resolution function
/// files": "the values of the major part of parameters are added to those given in the PCR file"; FullProf's News of
/// 18 April 2011 for ZD2TOF and ZD2TOT). The bank's angle (TWOTH) replaces the `.pcr`'s instead.
#[derive(Clone, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Instrument {
    /// Zero, Dtt1, Dtt2 (D2TOF Dtt1 Dtt2, or ZD2TOF Zero Dtt1 Dtt2).
    pub zero: f64,
    pub dtt: [f64; 2],
    /// Zerot, Dtt1t, Dtt2t, x-cross, Width (D2TOT Dtt1t Dtt2t x-cross Width, or ZD2TOT with Zerot first).
    pub thermal: [f64; 5],
    /// σ₀, σ₁, σ₂ (SIGMA Sig-2 Sig-1 Sig-0), γ₀, γ₁, γ₂ (GAMMA Gam-2 Gam-1 Gam-0), α₀, α₁ and β₀, β₁ (ALFBE alph0 beta0
    /// alph1 beta1), and α₀ₜ, β₀ₜ, α₁ₜ, β₁ₜ (ALFBT).
    pub sigma: [f64; 3],
    pub gamma: [f64; 3],
    pub alpha: [f64; 2],
    pub beta: [f64; 2],
    pub thermal_decay: [f64; 4],
}

impl Instrument {
    /// Its 22 values in one row, as the control file writes them: Zero Dtt1 Dtt2, Zerot Dtt1t Dtt2t x-cross Width,
    /// σ₀ σ₁ σ₂, γ₀ γ₁ γ₂, α₀ α₁ β₀ β₁, α₀ₜ β₀ₜ α₁ₜ β₁ₜ.
    pub fn to_values(&self) -> Vec<f64> {
        let mut v = vec![self.zero, self.dtt[0], self.dtt[1]];
        for part in [&self.thermal[..], &self.sigma, &self.gamma, &self.alpha, &self.beta, &self.thermal_decay] {
            v.extend_from_slice(part);
        }
        v
    }

    /// The inverse of `to_values` (22 values).
    pub fn from_values(v: &[f64]) -> Instrument {
        let a = |k: usize| v.get(k).copied().unwrap_or(0.0);
        Instrument {
            zero: a(0),
            dtt: [a(1), a(2)],
            thermal: [a(3), a(4), a(5), a(6), a(7)],
            sigma: [a(8), a(9), a(10)],
            gamma: [a(11), a(12), a(13)],
            alpha: [a(14), a(15)],
            beta: [a(16), a(17)],
            thermal_decay: [a(18), a(19), a(20), a(21)],
        }
    }

    /// A TOF resolution file (Res = 5, 6): its keyword lines (the first line is the title; '!' and '#' lines are
    /// comments). The bank's angle of TWOTH is returned apart.
    pub fn read(text: &str) -> Result<(Instrument, Option<f64>), String> {
        let mut i = Instrument::default();
        let mut twoth = None;
        let mut any = false;
        for l in text.lines().skip(1) {
            let t = l.trim_start();
            if t.starts_with('!') || t.starts_with('#') || t.is_empty() {
                continue;
            }
            let mut w = t.split_whitespace();
            let key = w.next().unwrap_or("").to_uppercase();
            let v: Vec<f64> = w.map_while(|x| x.parse().ok()).collect();
            let g = |k: usize| v.get(k).copied().unwrap_or(0.0);
            any |= matches!(key.as_str(), "D2TOF" | "ZD2TOF" | "SIGMA" | "GAMMA" | "ALFBE");
            match key.as_str() {
                "D2TOF" => (i.dtt, i.zero) = ([g(0), g(1)], g(2)),
                "ZD2TOF" => (i.zero, i.dtt) = (g(0), [g(1), g(2)]),
                "D2TOT" => i.thermal = [0.0, g(0), g(1), g(2), g(3)],
                "ZD2TOT" => i.thermal = [g(0), g(1), g(2), g(3), g(4)],
                "TWOTH" => twoth = v.first().copied(),
                "SIGMA" => i.sigma = [g(2), g(1), g(0)],
                "GAMMA" => i.gamma = [g(2), g(1), g(0)],
                "ALFBE" => (i.alpha, i.beta) = ([g(0), g(2)], [g(1), g(3)]),
                "ALFBT" => i.thermal_decay = [g(0), g(1), g(2), g(3)],
                _ => {}
            }
        }
        if !any {
            return Err("a TOF resolution file needs D2TOF, SIGMA, GAMMA or ALFBE lines".into());
        }
        Ok((i, twoth))
    }
}

/// A phase's peak shape in a TOF pattern (P05). Without Hodges' terms: α = α₀ + α₁/d + α_Q/√d, β = β₀ + β₁/d⁴ + β_Q/d²
/// (µs⁻¹). With them (`thermal`, the bank's `Thermal`): the epithermal α₀ + α₁ d and the thermal α₀ₜ − α₁ₜ/d are times
/// (µs), mixed as 1/α (α_Q, β_Q unused).
#[derive(Clone, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Profile {
    /// σ₀, σ₁, σ₂: σ² = σ₀ + σ₁ d² + σ₂ d⁴ + σ_Q d (µs²).
    pub sigma: [f64; 3],
    pub sigma_q: f64,
    /// γ₀, γ₁, γ₂: H_L = γ₀ + γ₁ d + γ₂ d² (µs).
    pub gamma: [f64; 3],
    /// α₀, α₁, α_Q and β₀, β₁, β_Q.
    pub alpha: [f64; 3],
    pub beta: [f64; 3],
    /// α₀ₜ, β₀ₜ, α₁ₜ, β₁ₜ (FullProf's line 31, Npr = 10).
    pub thermal: Option<[f64; 4]>,
    /// The extinction parameter (Sabine's, `extinction`; FullProf's Extinc).
    pub extinction: f64,
    /// A₁ of the absorption (FullProf's ABSCOR1), per Å of wavelength.
    pub absorption: f64,
}

/// What a reflection's TOF peak is: its position (µs), rise and decay (µs⁻¹), Gaussian variance (µs²) and Lorentzian
/// FWHM (µs), from its d-spacing, generic for the derivatives.
pub struct PeakShape<T> {
    pub position: T,
    pub alpha: T,
    pub beta: T,
    pub sigma2: T,
    pub lorentz: T,
}

impl Bank {
    /// The d-spacing (Å) whose time of flight is t (µs), with the model's zero and the instrument's part added: the
    /// inverse of T(d) (`Parameters::at`), by bisection over 10⁻³–10³ Å, to 10⁻¹² Å relative; None outside, or where
    /// T(d) is not increasing (P46: the profile and the reflections in d and Q).
    pub fn d_spacing(&self, zero: f64, t: f64) -> Option<f64> {
        let par = Parameters::of(zero, self, &Profile::default());
        let at = |d: f64| par.at(d).position;
        let (mut lo, mut hi) = (1e-3_f64, 1e3_f64);
        if !(at(lo) <= t && t <= at(hi)) {
            return None;
        }
        for _ in 0..200 {
            let m = 0.5 * (lo + hi);
            if at(m) < t {
                lo = m;
            } else {
                hi = m;
            }
            if hi - lo < 1e-12 * hi {
                break;
            }
        }
        Some(0.5 * (lo + hi))
    }

    /// The TOF–d equation as text, with its values (pdCIF's `_pd_calibration_conversion_eqn`; P46): the instrument's
    /// part added, as the calculation takes it.
    pub fn equation(&self, zero: f64) -> String {
        let p = Parameters::of(zero, self, &Profile::default());
        match p.thermal {
            Some([zt, d1t, d2t, cross, width]) => format!(
                "TOF = n (Zero + Dtt1 d) + (1 - n) (Zerot + Dtt1t d - Dtt2t/d), n = erfc(Width (x_cross - 1/d))/2; Zero = {} us, Dtt1 = {} us/A, Zerot = {} us, Dtt1t = {} us/A, Dtt2t = {} us A, x_cross = {} 1/A, Width = {} A (Hodges; FullProf's Npr = 10)",
                p.zero, p.dtt[0], zt, d1t, d2t, cross, width
            ),
            None => format!("TOF = Zero + Dtt1 d + Dtt2 d^2 + Dtt_1overd/d; Zero = {} us, Dtt1 = {} us/A, Dtt2 = {} us/A^2, Dtt_1overd = {} us A", p.zero, p.dtt[0], p.dtt[1], p.dtt[2]),
        }
    }

    /// The instrument's part (a resolution file's, added to the control file's values) of the bank's field k
    /// (`BANK_FIELDS`) and of the zero (k = None).
    pub fn instrument_part(&self, k: Option<usize>) -> f64 {
        let Some(i) = &self.instrument else { return 0.0 };
        match k {
            None => i.zero,
            Some(0) => i.dtt[0],
            Some(1) => i.dtt[1],
            Some(2) => 0.0,
            Some(k) => i.thermal[k - 3],
        }
    }

    /// The instrument's part of a phase's TOF value k (`PROFILE_FIELDS`).
    pub fn profile_instrument_part(&self, k: usize) -> f64 {
        let Some(i) = &self.instrument else { return 0.0 };
        match k {
            0..=2 => i.sigma[k],
            4..=6 => i.gamma[k - 4],
            7 | 8 => i.alpha[k - 7],
            10 | 11 => i.beta[k - 10],
            13..=16 => i.thermal_decay[k - 13],
            _ => 0.0,
        }
    }
}

/// The values the shape depends on, as the bank and the phase give them (generic for the derivatives).
pub struct Parameters<T> {
    pub zero: T,
    pub dtt: [T; 3],
    pub thermal: Option<[T; 5]>,
    pub sigma: [T; 3],
    pub sigma_q: T,
    pub gamma: [T; 3],
    pub alpha: [T; 3],
    pub beta: [T; 3],
    pub thermal_decay: Option<[T; 4]>,
}

impl Parameters<f64> {
    pub fn of(zero: f64, bank: &Bank, p: &Profile) -> Parameters<f64> {
        let i = bank.instrument.clone().unwrap_or_default();
        let add = |a: [f64; 3], b: [f64; 3]| [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
        Parameters {
            zero: zero + i.zero,
            dtt: [bank.dtt1 + i.dtt[0], bank.dtt2 + i.dtt[1], bank.dtt_1overd],
            thermal: bank.thermal.map(|t| [t.zero + i.thermal[0], t.dtt1 + i.thermal[1], t.dtt2 + i.thermal[2], t.cross + i.thermal[3], t.width + i.thermal[4]]),
            sigma: add(p.sigma, i.sigma),
            sigma_q: p.sigma_q,
            gamma: add(p.gamma, i.gamma),
            alpha: [p.alpha[0] + i.alpha[0], p.alpha[1] + i.alpha[1], p.alpha[2]],
            beta: [p.beta[0] + i.beta[0], p.beta[1] + i.beta[1], p.beta[2]],
            thermal_decay: p.thermal.map(|t| [t[0] + i.thermal_decay[0], t[1] + i.thermal_decay[1], t[2] + i.thermal_decay[2], t[3] + i.thermal_decay[3]]),
        }
    }
}

/// A scan read as 2θ from a GSAS CONST bank put on the time-of-flight axis: its x × 100 in µs, as GSAS-II multiplies it
/// back once the instrument says TOF (`formats::gsas`). A scan already in TOF is left as it is.
pub fn time_axis(scan: &mut crate::formats::Scan) {
    if scan.axis == crate::formats::Axis::TwoTheta {
        for x in &mut scan.pattern.x {
            *x *= 100.0;
        }
        scan.axis = crate::formats::Axis::TimeOfFlight;
    }
}

/// erfc for the dual numbers: d/dx erfc x = −(2/√π) e^{−x²}.
pub fn erfc<T: Real>(x: T) -> T {
    let v = x.val();
    x.map(libm_erfc(v), -std::f64::consts::FRAC_2_SQRT_PI * (-v * v).exp())
}

/// The complementary error function (W. J. Cody's rational approximations are not in std; this is the continued
/// fraction and series of Numerical Recipes' erfc with the precision of f64 to 10⁻¹⁵ relative).
pub fn libm_erfc(x: f64) -> f64 {
    if x < 0.0 {
        return 2.0 - libm_erfc(-x);
    }
    if x < 2.0 {
        // erf by its Taylor series.
        let (mut term, mut sum, x2) = (x, x, x * x);
        for n in 1..200 {
            term *= -x2 / n as f64;
            let t = term / (2 * n + 1) as f64;
            sum += t;
            if t.abs() < 1e-17 * sum.abs() {
                break;
            }
        }
        return 1.0 - std::f64::consts::FRAC_2_SQRT_PI * sum;
    }
    // The continued fraction erfc x = e^{−x²}/√π · 1/(x + ½/(x + 1/(x + 3/2/(x + …)))), evaluated backwards.
    let mut f = 0.0;
    for k in (1..120).rev() {
        f = (k as f64 / 2.0) / (x + f);
    }
    (-x * x).exp() / std::f64::consts::PI.sqrt() / (x + f)
}

impl<T: Real> Parameters<T> {
    /// The shape at d (Å).
    pub fn at(&self, d: T) -> PeakShape<T> {
        let d2 = d * d;
        let sigma2 = self.sigma[0] + self.sigma[1] * d2 + self.sigma[2] * d2 * d2 + self.sigma_q * d;
        let lorentz = self.gamma[0] + self.gamma[1] * d + self.gamma[2] * d2;
        match (self.thermal, self.thermal_decay) {
            (Some([zt, d1t, d2t, cross, width]), Some([a0t, b0t, a1t, b1t])) => {
                let n = erfc(width * (cross - T::cst(1.0) / d)) * 0.5;
                let te = self.zero + self.dtt[0] * d;
                let tt = zt + d1t * d - d2t / d;
                let position = n * te + (-n + 1.0) * tt;
                let ia = n * (self.alpha[0] + self.alpha[1] * d) + (-n + 1.0) * (a0t - a1t / d);
                let ib = n * (self.beta[0] + self.beta[1] * d) + (-n + 1.0) * (b0t - b1t / d);
                PeakShape { position, alpha: T::cst(1.0) / ia, beta: T::cst(1.0) / ib, sigma2, lorentz }
            }
            _ => PeakShape {
                position: self.zero + self.dtt[0] * d + self.dtt[1] * d2 + self.dtt[2] / d,
                alpha: self.alpha[0] + self.alpha[1] / d + self.alpha[2] / d.sqrt(),
                beta: self.beta[0] + self.beta[1] / (d2 * d2) + self.beta[2] / d2,
                sigma2,
                lorentz,
            },
        }
    }
}

/// The Gaussian FWHM from the variance: H_G = √(8 ln 2 σ²).
pub fn gauss_fwhm<T: Real>(sigma2: T) -> T {
    (sigma2.max0() * (8.0 * std::f64::consts::LN_2)).sqrt()
}

/// Sabine's extinction of a reflection (Sabine, Von Dreele & Jørgensen 1988, Acta Cryst. A44, 374–379, as GSAS has it):
/// x = Ext·(λF/V)², E = E_L cos²θ_B + E_B sin²θ_B, E_B = 1/√(1 + x), E_L = 1 − x/2 + x²/4 − 5x³/48 + 7x⁴/192 for x ≤ 1
/// and √(2/(πx))(1 − 1/(8x) − 3/(128x²) − 15/(1024x³)) above. FullProf's TOF Extinc is this Ext: its CORR in
/// `arg_si_ext10.out` and `arg_si_ext100.out` (absorption 0) gives x/(Ext·F²λ²) = 3.8953·10⁻⁵ = 1/V² (V = 160.224 Å³)
/// on every reflection (the manual, p. 63, names only SHELX's empirical correction). The paper is not in
/// `manuals/papers/`: the form is cited from it as GSAS-II codes it, and checked on FullProf's runs.
pub fn extinction<T: Real>(ext: T, f2: T, lambda: T, volume: T, sin2_theta: f64) -> T {
    // E = 1 at x = 0, and its slope there, −½ ∂x, kept (an Extinc refined from 0 starts there).
    let x = ext * f2 * lambda * lambda / (volume * volume);
    let eb = T::cst(1.0) / (x + 1.0).sqrt();
    let el = if x.val() <= 1.0 {
        let (x2, x3) = (x * x, x * x * x);
        -x * 0.5 + x2 * 0.25 - x3 * (5.0 / 48.0) + x3 * x * (7.0 / 192.0) + 1.0
    } else {
        let ix = T::cst(1.0) / x;
        (ix * (2.0 / std::f64::consts::PI)).sqrt() * (-ix * 0.125 - ix * ix * (3.0 / 128.0) - ix * ix * ix * (15.0 / 1024.0) + 1.0)
    };
    el * (1.0 - sin2_theta) + eb * sin2_theta
}

/// e^z E₁(z) for complex z (CrysFML's `Expi_E1`, `CFML_Profiles.f90`, from Zhang & Jin, Computation of Special
/// Functions, 1996): the series and the continued fraction. CrysFML takes the series for |z| ≤ 10 (or Re z ≤ 0 and
/// |z| < 20); for Re z > 0 the series' terms reach 10³ where E₁ ≈ e^{−z}/z is 10⁻⁶, and its relative error grows to
/// 1.2·10⁻⁸ at |z| = 10 (against mpmath's e^z E₁(z) at 40 digits, Im z > 0), noise that the profile's derivatives by
/// finite differences magnify to 10⁻⁴ (L341). Here the continued fraction is taken for Re z ≥ 0 from |z| = 2 on, where
/// it holds 3·10⁻¹⁶ (the series 1.2·10⁻¹⁴ at |z| = 2); for Re z < 0, CrysFML's choice, the better of the two there
/// (the series 2.9·10⁻⁹ at |z| = 20, the continued fraction 1.9·10⁻⁷; 6.4·10⁻⁸ and 1.6·10⁻⁸ at 22, for Im z ≥ 0.1).
fn exp_e1(re: f64, im: f64) -> (f64, f64) {
    use std::f64::consts::PI;
    const EULER: f64 = 0.577_215_664_901_532_8;
    let mul = |a: (f64, f64), b: (f64, f64)| (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0);
    let div = |a: (f64, f64), b: (f64, f64)| {
        let n = b.0 * b.0 + b.1 * b.1;
        ((a.0 * b.0 + a.1 * b.1) / n, (a.1 * b.0 - a.0 * b.1) / n)
    };
    let z = (re, im);
    let a0 = re.hypot(im);
    let ez = {
        let m = re.exp();
        (m * im.cos(), m * im.sin())
    };
    if a0 == 0.0 {
        return (1e300, 0.0);
    }
    if if re >= 0.0 { a0 < 2.0 } else { a0 < 20.0 } {
        let (mut e, mut cr) = ((1.0, 0.0), (1.0, 0.0));
        for k in 1..=150 {
            let km = k as f64;
            let t = mul(cr, z);
            cr = (-t.0 * km / ((km + 1.0) * (km + 1.0)), -t.1 * km / ((km + 1.0) * (km + 1.0)));
            e = (e.0 + cr.0, e.1 + cr.1);
            if cr.0.hypot(cr.1) <= e.0.hypot(e.1) * 1e-15 {
                break;
            }
        }
        let ln = (a0.ln(), im.atan2(re));
        let ze = mul(z, e);
        let s = (-EULER - ln.0 + ze.0, -ln.1 + ze.1);
        mul(s, ez)
    } else {
        let mut ct = (0.0, 0.0);
        for k in (1..=120).rev() {
            let km = k as f64;
            let q = div((km, 0.0), (z.0 + ct.0, z.1 + ct.1));
            ct = div((km, 0.0), (1.0 + q.0, q.1));
        }
        let mut r = div((1.0, 0.0), (z.0 + ct.0, z.1 + ct.1));
        if re <= 0.0 && im.abs() <= 1e-9 {
            r = (r.0 + PI * ez.1, r.1 - PI * ez.0);
        }
        r
    }
}

/// The value of the back-to-back exponentials ⊗ pseudo-Voigt at dt = T − T_h (µs) and its derivatives with respect to
/// dt, α, β, H and η (CrysFML's `Tof_Jorgensen_Vondreele` and its `deriv`), H the pseudo-Voigt's FWHM.
pub fn back_to_back(dt: f64, alpha: f64, beta: f64, h: f64, eta: f64) -> (f64, [f64; 5]) {
    use std::f64::consts::{FRAC_2_SQRT_PI, PI};
    let inv_8ln2 = 1.0 / (8.0 * std::f64::consts::LN_2);
    let sigma = h * h * inv_8ln2;
    let u = 0.5 * alpha * (alpha * sigma + 2.0 * dt);
    let v = 0.5 * beta * (beta * sigma - 2.0 * dt);
    let norm = 0.5 * alpha * beta / (alpha + beta);
    let deno = (2.0 * sigma).sqrt();
    let y = (alpha * sigma + dt) / deno;
    let z = (beta * sigma - dt) / deno;
    // e^u erfc(y), computed as exp(u + ln erfc y) where erfc underflows, so that a large u does not overflow.
    let e_erfc = |u: f64, y: f64| {
        let e = libm_erfc(y);
        if e > 0.0 {
            (u + e.ln()).exp()
        } else {
            0.0
        }
    };
    let omega = e_erfc(u, y);
    let omegb = e_erfc(v, z);
    let omeg = norm * (omega + omegb);
    let (fz1, fz2) = (exp_e1(alpha * dt, 0.5 * alpha * h), exp_e1(-beta * dt, 0.5 * beta * h));
    let two_over_pi = 2.0 / PI;
    let (oml_a, oml_b) = (-fz1.1 * two_over_pi, -fz2.1 * two_over_pi);
    let omel = norm * (oml_a + oml_b);
    let value = (1.0 - eta) * omeg + eta * omel;
    // Derivatives (CrysFML's, term for term).
    let (a2, b2) = (alpha * alpha, beta * beta);
    let denoinv = 1.0 / deno;
    let d3 = 2.0 * denoinv * denoinv * denoinv;
    let ca = d3 * (alpha * sigma - dt);
    let cb = d3 * (beta * sigma + dt);
    let (domg_t, domg_a, domg_b, domg_g) = if omeg <= 1e-30 {
        (0.0, 0.0, 0.0, 0.0)
    } else {
        // erfc′(x) = −(2/√π) e^{−x²}, times e^u.
        let a = -FRAC_2_SQRT_PI * (u - y * y).exp();
        let b = -FRAC_2_SQRT_PI * (v - z * z).exp();
        (
            norm * (alpha * omega - beta * omegb + (a - b) * denoinv),
            norm * (2.0 * omeg / a2 + deno * (y * omega + 0.5 * a)),
            norm * (2.0 * omeg / b2 + deno * (z * omegb + 0.5 * b)),
            inv_8ln2 * norm * h * (a2 * omega + b2 * omegb + a * ca + b * cb),
        )
    };
    let (exper1, exper2) = (-fz1.0 * two_over_pi, -fz2.0 * two_over_pi);
    let doml_t = norm * (alpha * oml_a - beta * oml_b);
    let doml_a = norm * (2.0 * omel / a2 + dt * oml_a + 0.5 * h * exper1);
    let doml_b = norm * (2.0 * omel / b2 - dt * oml_b + 0.5 * h * exper2);
    let doml_g = 0.5 * norm * (alpha * exper1 + beta * exper2);
    let one_e = 1.0 - eta;
    (value, [one_e * domg_t + eta * doml_t, one_e * domg_a + eta * doml_a, one_e * domg_b + eta * doml_b, one_e * domg_g + eta * doml_g, omel - omeg])
}

/// `back_to_back` for dual numbers: its value with the chain rule through its five arguments.
pub fn back_to_back_g<T: Real>(dt: T, alpha: T, beta: T, h: T, eta: T) -> T {
    let (v, d) = back_to_back(dt.val(), alpha.val(), beta.val(), h.val(), eta.val());
    let lin = |x: T, k: usize| (x - x.val()) * d[k];
    lin(dt, 0) + lin(alpha, 1) + lin(beta, 2) + lin(h, 3) + lin(eta, 4) + v
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn erfc_and_the_profile() {
        // erfc against its known values (Abramowitz & Stegun, Table 7.1).
        for (x, e) in [(0.0, 1.0), (0.5, 0.479_500_122_186_953_5), (1.0, 0.157_299_207_050_285_1), (2.0, 0.004_677_734_981_047_266), (3.0, 2.209_049_699_858_544e-5), (-1.0, 1.842_700_792_949_715)] {
            assert!((libm_erfc(x) - e).abs() < 1e-14 * e.max(1e-3), "erfc({x}) = {} for {e}", libm_erfc(x));
        }
        // The profile integrates to 1 and its derivatives are its finite differences.
        let (a, b, h, eta) = (0.19, 0.0423, 47.15, 0.218);
        // Over ±600 µs, plus the Lorentzian's tails beyond, η·H/(π·600) (the exponentials' share there is below 1e-9).
        let s: f64 = (-6000..6000).map(|i| back_to_back(i as f64 * 0.1, a, b, h, eta).0 * 0.1).sum();
        let tails = eta * h / (std::f64::consts::PI * 600.0);
        assert!((s + tails - 1.0).abs() < 2e-4, "integral {s} + tails {tails}");
        let x = [12.0, a, b, h, eta];
        let f = |x: [f64; 5]| back_to_back(x[0], x[1], x[2], x[3], x[4]).0;
        let (_, d) = back_to_back(x[0], x[1], x[2], x[3], x[4]);
        for k in 0..5 {
            let e = 1e-6 * x[k].abs().max(1e-3);
            let (mut p, mut m) = (x, x);
            p[k] += e;
            m[k] -= e;
            let fd = (f(p) - f(m)) / (2.0 * e);
            assert!((fd - d[k]).abs() < 1e-6 * fd.abs().max(1e-6), "derivative {k}: {} analytic, {fd} numeric", d[k]);
        }
    }
}
