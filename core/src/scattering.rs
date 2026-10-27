//! Scattering of one atom: X-ray form factors with the anomalous terms f′ and f″, neutron scattering lengths.
//!
//! The numbers come from CrysFML's tables (`scattering_table`, generated from `CFML_Scattering_Chemical_Tables.f90`),
//! the tables FullProf prints in its `.out` files ("Scattering coefficients from internal table"):
//! - X-rays: f(s) = Σ_{i=1..4} a_i exp(−b_i s²) + c, s = sinθ/λ (Å⁻¹), in electrons (International Tables Vol. C,
//!   Table 6.1.1.4). The species is looked up as written, case aside: "Ti+4" is the Ti⁴⁺ ion, "Ti" the neutral atom;
//!   a species that is not in the table is an error, as in CrysFML (`Create_Table_fabc_Xray`).
//! - Anomalous terms: f′ and f″ of the element, tabulated at the Kα1 lines of Cr, Fe, Cu, Mo and Ag; for a wavelength λ
//!   the anode whose Kα1 is nearest to λ is used (CrysFML, `Create_Table_fabc_Xray`: min |λ − Kα1|). FullProf does
//!   the same: at λ = 1.001 Å it prints the Mo values (test `fullprof_structure_factors`).
//! - Neutrons: the coherent scattering length b of the element (charge removed), in 10⁻¹² cm, as FullProf uses it.
//! - The user's own (P11, L233; FullProf's line 12, manual pp. 70 and 79): a neutron b, or for X-rays f′ and f″ with
//!   the tabulated f(s), the nine coefficients of f(s), or a table of (sin θ/λ, f) interpolated linearly (`UserScatterer`).

use crate::scattering_table::{ANODES, ANOMALOUS, NEUTRON_B, XRAY_FORM};
use std::fmt;
use std::ops::{Add, AddAssign, Mul, Sub};

/// A complex number, enough for structure factors (the core has no external dependencies).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Complex {
    pub re: f64,
    pub im: f64,
}

impl Complex {
    pub const ZERO: Complex = Complex { re: 0.0, im: 0.0 };

    pub fn new(re: f64, im: f64) -> Complex {
        Complex { re, im }
    }

    /// exp(iφ).
    pub fn cis(phi: f64) -> Complex {
        let (s, c) = phi.sin_cos();
        Complex { re: c, im: s }
    }

    pub fn conj(self) -> Complex {
        Complex { re: self.re, im: -self.im }
    }

    /// |z|².
    pub fn norm_sqr(self) -> f64 {
        self.re * self.re + self.im * self.im
    }

    pub fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }

    /// The phase in degrees, −180 < φ ≤ 180.
    pub fn phase_deg(self) -> f64 {
        self.im.atan2(self.re).to_degrees()
    }
}

impl Add for Complex {
    type Output = Complex;
    fn add(self, o: Complex) -> Complex {
        Complex { re: self.re + o.re, im: self.im + o.im }
    }
}

impl Sub for Complex {
    type Output = Complex;
    fn sub(self, o: Complex) -> Complex {
        Complex { re: self.re - o.re, im: self.im - o.im }
    }
}

impl AddAssign for Complex {
    fn add_assign(&mut self, o: Complex) {
        self.re += o.re;
        self.im += o.im;
    }
}

impl Mul for Complex {
    type Output = Complex;
    fn mul(self, o: Complex) -> Complex {
        Complex { re: self.re * o.re - self.im * o.im, im: self.re * o.im + self.im * o.re }
    }
}

impl Mul<f64> for Complex {
    type Output = Complex;
    fn mul(self, k: f64) -> Complex {
        Complex { re: self.re * k, im: self.im * k }
    }
}

impl fmt::Display for Complex {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}{:+}i", self.re, self.im)
    }
}

/// The radiation a structure factor is computed for.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Radiation {
    /// X-rays of wavelength λ (Å). With `anomalous`, f′ and f″ of the nearest anode are added.
    Xray { lambda: f64, anomalous: bool },
    /// Neutrons: scattering lengths, independent of λ and s.
    Neutron,
}

/// The anode names of the table, with Kα1, Kα2 and Kβ (Å): Cr, Fe, Cu, Mo, Ag, Co, Ni.
pub fn anode(name: &str) -> Option<(&'static str, f64, f64, f64)> {
    ANODES.iter().find(|a| a.0.eq_ignore_ascii_case(name)).copied()
}

/// The column of the anomalous table for wavelength λ: the anode among Cr, Fe, Cu, Mo, Ag whose Kα1 is nearest to λ
/// (the first one on a tie, as CrysFML's loop keeps it).
pub fn nearest_anode(lambda: f64) -> usize {
    let mut best = 0;
    for (i, a) in ANODES.iter().enumerate().take(5) {
        if (lambda - a.1).abs() < (lambda - ANODES[best].1).abs() {
            best = i;
        }
    }
    best
}

/// The element of a species: its leading letters ("Ti+4" → "Ti", "O-2" → "O", "CU" → "CU").
pub fn element_of(species: &str) -> &str {
    let s = species.trim();
    let n = s.find(|c: char| !c.is_ascii_alphabetic()).unwrap_or(s.len());
    &s[..n]
}

/// The scattering of one species for one radiation, ready to be evaluated at any s = sinθ/λ.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Scatterer {
    /// X-ray form factor coefficients and the anomalous terms (zero when not wanted or not tabulated).
    Xray { symbol: String, a: [f64; 4], b: [f64; 4], c: f64, fp: f64, fpp: f64 },
    /// A neutron scattering length (10⁻¹² cm), with the imaginary part of an absorbing nucleus (0 in the tables; the
    /// user's DFPP, as FullProf uses it: nabd4_calc.out prints "Scattering lengths (real & imaginary parts) … For B
    /// 0.530000 −0.213000", though its manual, p. 79, says DFPP is ignored for neutrons; L251).
    Neutron {
        symbol: String,
        b: f64,
        #[cfg_attr(feature = "serde", serde(default))]
        b_imaginary: f64,
    },
    /// An X-ray form factor given as a table of (sin θ/λ, f), interpolated linearly (FullProf's Nsc < 0, manual p. 70),
    /// held at its end values outside it, with the anomalous terms.
    XrayTable { symbol: String, points: Vec<(f64, f64)>, fp: f64, fpp: f64 },
}

/// A scattering factor the user gives for a species (FullProf's line 12; the control file's `[[pattern.scatterer]]`):
/// for neutrons the scattering length b; for X-rays f′ and f″ with the tabulated f(s) of the species, or with the nine
/// coefficients a₁ b₁ … a₄ b₄ c of f(s) (FullProf's order), or with a table of (sin θ/λ, f).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct UserScatterer {
    /// The species it stands for, matched case aside (FullProf turns X-ray names to lower case, p. 80).
    pub species: String,
    pub scattering_length: Option<f64>,
    /// The imaginary part of b (an absorbing nucleus; FullProf's DFPP for neutrons).
    pub scattering_length_imaginary: f64,
    pub f_prime: f64,
    pub f_double_prime: f64,
    pub form_factor: Option<[f64; 9]>,
    pub form_factor_table: Vec<(f64, f64)>,
}

impl UserScatterer {
    /// The scatterer it gives for a radiation.
    pub fn scatterer(&self, radiation: Radiation) -> Result<Scatterer, String> {
        let sp = &self.species;
        match radiation {
            Radiation::Neutron => match self.scattering_length {
                Some(b) => Ok(Scatterer::Neutron { symbol: sp.clone(), b, b_imaginary: self.scattering_length_imaginary }),
                None => Err(format!("the scatterer '{sp}' gives no scattering length for neutrons")),
            },
            Radiation::Xray { .. } => {
                let (fp, fpp) = (self.f_prime, self.f_double_prime);
                if let Some(k) = self.form_factor {
                    // "If A2 is set to zero, only A1 and B1 are used" (manual p. 79).
                    let (a, b, c) = if k[2] == 0.0 { ([k[0], 0.0, 0.0, 0.0], [k[1], 0.0, 0.0, 0.0], 0.0) } else { ([k[0], k[2], k[4], k[6]], [k[1], k[3], k[5], k[7]], k[8]) };
                    Ok(Scatterer::Xray { symbol: sp.clone(), a, b, c, fp, fpp })
                } else if !self.form_factor_table.is_empty() {
                    let mut points = self.form_factor_table.clone();
                    points.sort_by(|x, y| x.0.total_cmp(&y.0));
                    Ok(Scatterer::XrayTable { symbol: sp.clone(), points, fp, fpp })
                } else {
                    match Scatterer::lookup_or_element(sp, Radiation::Xray { lambda: 1.0, anomalous: false })?.0 {
                        Scatterer::Xray { symbol, a, b, c, .. } => Ok(Scatterer::Xray { symbol, a, b, c, fp, fpp }),
                        other => Ok(other),
                    }
                }
            }
        }
    }

    /// What a scatterer differs in from the tables, as the user's (None when it is the tables' own, `lookup_or_element`).
    pub fn of(species: &str, s: &Scatterer, radiation: Radiation) -> Option<UserScatterer> {
        let table = Scatterer::lookup_or_element(species, radiation).ok().map(|t| t.0);
        if table.as_ref() == Some(s) {
            return None;
        }
        let mut u = UserScatterer { species: species.to_string(), ..Default::default() };
        match s {
            Scatterer::Neutron { b, b_imaginary, .. } => (u.scattering_length, u.scattering_length_imaginary) = (Some(*b), *b_imaginary),
            Scatterer::Xray { a, b, c, fp, fpp, .. } => {
                (u.f_prime, u.f_double_prime) = (*fp, *fpp);
                let same_f = matches!(&table, Some(Scatterer::Xray { a: ta, b: tb, c: tc, .. }) if ta == a && tb == b && tc == c);
                if !same_f {
                    u.form_factor = Some([a[0], b[0], a[1], b[1], a[2], b[2], a[3], b[3], *c]);
                }
            }
            Scatterer::XrayTable { points, fp, fpp, .. } => {
                (u.f_prime, u.f_double_prime) = (*fp, *fpp);
                u.form_factor_table = points.clone();
            }
        }
        Some(u)
    }
}

/// The scatterer of a species: the user's when one names it (case aside), else the tables' (`lookup_or_element`).
pub fn resolve(species: &str, radiation: Radiation, user: &[UserScatterer]) -> Result<Scatterer, String> {
    match user.iter().find(|u| u.species.trim().eq_ignore_ascii_case(species.trim())) {
        Some(u) => u.scatterer(radiation),
        None => Scatterer::lookup_or_element(species, radiation).map(|s| s.0),
    }
}

impl Scatterer {
    /// The table entry for a species ("Ti+4", "O", "pb", "CU") and a radiation.
    pub fn lookup(species: &str, radiation: Radiation) -> Result<Scatterer, String> {
        match radiation {
            Radiation::Xray { lambda, anomalous } => {
                let key = species.trim().to_ascii_lowercase();
                let &(symbol, _, a, b, c) = XRAY_FORM
                    .iter()
                    .find(|e| e.0 == key)
                    .ok_or_else(|| format!("species '{species}' is not in the X-ray form-factor table (CrysFML Xray_Form)"))?;
                let (mut fp, mut fpp) = (0.0, 0.0);
                if anomalous {
                    let el = element_of(&key);
                    if let Some(e) = ANOMALOUS.iter().find(|e| e.0 == el) {
                        let k = nearest_anode(lambda);
                        fp = e.1[k];
                        fpp = e.2[k];
                    }
                }
                Ok(Scatterer::Xray { symbol: symbol.to_string(), a, b, c, fp, fpp })
            }
            Radiation::Neutron => {
                let el = element_of(species).to_ascii_uppercase();
                let &(symbol, _, b) = NEUTRON_B
                    .iter()
                    .find(|e| e.0 == el)
                    .ok_or_else(|| format!("element '{el}' of species '{species}' is not in the neutron table (CrysFML Chem_Info)"))?;
                Ok(Scatterer::Neutron { symbol: symbol.to_string(), b, b_imaginary: 0.0 })
            }
        }
    }

    /// As `lookup`, but an ion that is not in the X-ray table falls back to its neutral element (the second value is then
    /// true). CrysFML and FullProf stop instead; the command line uses this for CIF type symbols such as "S2-".
    pub fn lookup_or_element(species: &str, radiation: Radiation) -> Result<(Scatterer, bool), String> {
        match Scatterer::lookup(species, radiation) {
            Ok(s) => Ok((s, false)),
            Err(e) => match Scatterer::lookup(element_of(species), radiation) {
                Ok(s) if element_of(species) != species.trim() => Ok((s, true)),
                _ => Err(e),
            },
        }
    }

    /// The scattering amplitude at s = sinθ/λ (Å⁻¹): f(s) + f′ + i f″ for X-rays, b for neutrons.
    pub fn at(&self, s: f64) -> Complex {
        match self {
            Scatterer::Xray { a, b, c, fp, fpp, .. } => {
                let s2 = s * s;
                let f0: f64 = (0..4).map(|i| a[i] * (-b[i] * s2).exp()).sum::<f64>() + c;
                Complex::new(f0 + fp, *fpp)
            }
            Scatterer::Neutron { b, b_imaginary, .. } => Complex::new(*b, *b_imaginary),
            Scatterer::XrayTable { points, fp, fpp, .. } => Complex::new(table_at(points, s).0 + fp, *fpp),
        }
    }

    /// d f/d s at s (the anomalous terms and neutron lengths do not depend on s).
    pub fn ds(&self, s: f64) -> f64 {
        match self {
            Scatterer::Xray { a, b, .. } => -2.0 * s * (0..4).map(|i| a[i] * b[i] * (-b[i] * s * s).exp()).sum::<f64>(),
            Scatterer::Neutron { .. } => 0.0,
            Scatterer::XrayTable { points, .. } => table_at(points, s).1,
        }
    }
}

/// A table of (s, f) at s, linearly: the value and the slope of the segment (0 outside the table, held at its ends).
fn table_at(points: &[(f64, f64)], s: f64) -> (f64, f64) {
    match points {
        [] => (0.0, 0.0),
        [p] => (p.1, 0.0),
        _ if s <= points[0].0 => (points[0].1, 0.0),
        _ if s >= points[points.len() - 1].0 => (points[points.len() - 1].1, 0.0),
        _ => {
            let i = points.partition_point(|p| p.0 <= s).max(1);
            let (a, b) = (points[i - 1], points[i]);
            let slope = (b.1 - a.1) / (b.0 - a.0);
            (a.1 + slope * (s - a.0), slope)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn form_factor_at_zero_is_close_to_the_electron_count() {
        // f(0) = Σa + c; FullProf prints f(0) next to Z: 14.9993 for P, 5.9992 for C, 1.0000 for H (PPH3-sls.out:423-425).
        for (sp, f0) in [("P", 14.9993), ("C", 5.9992), ("H", 1.0000)] {
            let s = Scatterer::lookup(sp, Radiation::Xray { lambda: 1.0, anomalous: false }).unwrap();
            assert!((s.at(0.0).re - f0).abs() < 1e-4, "{sp}");
        }
    }

    #[test]
    fn ions_and_case() {
        let a = Scatterer::lookup("TI+4", Radiation::Xray { lambda: 1.54, anomalous: false }).unwrap();
        let b = Scatterer::lookup("ti+4", Radiation::Xray { lambda: 1.54, anomalous: false }).unwrap();
        assert_eq!(a, b);
        assert!(Scatterer::lookup("Ti+9", Radiation::Xray { lambda: 1.54, anomalous: false }).is_err());
        // Neutrons take the element: Ti+4 scatters as Ti.
        assert_eq!(
            Scatterer::lookup("Ti+4", Radiation::Neutron).unwrap(),
            Scatterer::lookup("Ti", Radiation::Neutron).unwrap()
        );
    }

    #[test]
    fn nearest_anode_columns() {
        // Columns: Cr, Fe, Cu, Mo, Ag.
        assert_eq!(nearest_anode(1.5406), 2);
        assert_eq!(nearest_anode(1.001), 3);
        assert_eq!(nearest_anode(0.5), 4);
        assert_eq!(nearest_anode(1.9095), 1);
        assert_eq!(anode("cu").unwrap().1, 1.54059);
    }

    #[test]
    fn complex_arithmetic() {
        let z = Complex::new(1.0, 2.0) * Complex::new(3.0, -1.0);
        assert_eq!(z, Complex::new(5.0, 5.0));
        assert!((Complex::cis(std::f64::consts::FRAC_PI_2).im - 1.0).abs() < 1e-15);
    }
}
