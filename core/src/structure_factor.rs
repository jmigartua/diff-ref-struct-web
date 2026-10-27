//! Structure factors, with FullProf's conventions.
//!
//! F(h) = Σ_j Occ_j · f_j(s) · exp(−B_j s²) · Σ_{g ∈ G} T_jg(h) · exp(2πi h·(R_g x_j + t_g)),   s = sinθ/λ = 1/(2d),
//!
//! the inner sum over all the operators of the group, centring and inversion included (M of them, the multiplicity of
//! the general position), and T_jg(h) = exp(−(h R_g) β_j (h R_g)ᵀ) the anisotropic factor of the image at R_g x_j + t_g,
//! whose betas are R_g β_j R_gᵀ.
//!
//! - Occ is FullProf's occupancy (manual p. 102: "Occ = chemical occupancy × site multiplicity (can be normalised to
//!   the multiplicity of the general position)"): an atom fully occupying a site of multiplicity m has Occ = m/M.
//!   Summing over all M operators counts each of its m distinct positions M/m times, so Occ = m/M weighs every
//!   position once. `occupancy_from_chemical` converts.
//! - B_j is FullProf's Biso (Å²); the betas are FullProf's, exp(−(β11h² + β22k² + β33l² + 2β12hk + 2β13hl + 2β23kl));
//!   `beta_from_u` converts the U_ij of a CIF.
//! - f_j(s) is the X-ray form factor with f′ and f″ (electrons), or the neutron scattering length (10⁻¹² cm); see
//!   `scattering`.
//!
//! FullProf's printed F² (column `F2` of the `.hkl` file written with Jbt/Hkl output, `StrFactor^2` of the `.out`) is
//! |F|² of this sum, anomalous terms included, without the overall temperature factor B_ov: on PPH3 (P 21/c, 7538
//! reflections) and citric acid (P 21/c, 4074 reflections) it agrees with every reflection (test
//! `fullprof_structure_factors`). The manual's eq. 3.44 (p. 57) sums over the reduced set of operators instead, and puts
//! the lattice and centrosymmetry factor NLAT = (LAT·ICEN)² on F² (p. 58); the full sum here contains both.

use crate::cell::Cell;
use crate::scattering::{Complex, Radiation, Scatterer};
use crate::site::Site;
use crate::spacegroup::SpaceGroup;
use crate::symop::DEN;
use std::f64::consts::PI;

/// An atom of the asymmetric unit, as in a FullProf `.pcr` file.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Atom {
    pub label: String,
    /// Scattering species: element or ion as the tables write it ("Ti+4", "O", "Cu").
    pub species: String,
    /// Fractional coordinates.
    pub x: [f64; 3],
    /// Isotropic displacement B (Å²).
    pub biso: f64,
    /// FullProf's Occ: chemical occupancy × m/M (m the multiplicity of the site, M of the general position).
    pub occ: f64,
    /// Anisotropic betas β11 β22 β33 β12 β13 β23 (FullProf's), applied on top of Biso.
    pub beta: Option<[f64; 6]>,
    /// The formal oxidation state, for bond-valence sums and the bonds' rule (P5, L114); never the scattering species,
    /// which keeps its own name ("Y" with charge 3 scatters as Y). None: taken from the species ("SR+2"), if it has one.
    #[cfg_attr(feature = "serde", serde(default))]
    pub charge: Option<i32>,
}

/// A crystal structure: space group, cell, atoms.
#[derive(Clone, Debug)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Structure {
    pub group: SpaceGroup,
    pub cell: Cell,
    pub atoms: Vec<Atom>,
}

/// FullProf's Occ of an atom at x with a chemical occupancy: chemical × m/M.
pub fn occupancy_from_chemical(g: &SpaceGroup, x: [f64; 3], chemical: f64) -> Result<f64, String> {
    let m = Site::new(g, x)?.multiplicity;
    Ok(chemical * m as f64 / g.multiplicity() as f64)
}

/// The chemical occupancy of an atom at x with FullProf's Occ: Occ × M/m.
pub fn chemical_from_occupancy(g: &SpaceGroup, x: [f64; 3], occ: f64) -> Result<f64, String> {
    let m = Site::new(g, x)?.multiplicity;
    Ok(occ * g.multiplicity() as f64 / m as f64)
}

/// FullProf's betas from the U_ij (Å²) of a CIF: β_ij = 2π² a*_i a*_j U_ij, ordered 11 22 33 12 13 23.
pub fn beta_from_u(cell: &Cell, u: [f64; 6]) -> [f64; 6] {
    let r = cell.reciprocal();
    let s = [r.a, r.b, r.c];
    let pairs = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];
    let mut b = [0.0; 6];
    for (k, &(i, j)) in pairs.iter().enumerate() {
        b[k] = 2.0 * PI * PI * s[i] * s[j] * u[k];
    }
    b
}

/// B (Å²) from U (Å²): B = 8π²U.
pub fn biso_from_uiso(u: f64) -> f64 {
    8.0 * PI * PI * u
}

/// B and its su from U and its su (Publish P20): B = 8π²U, σ_B = 8π²σ_U.
pub fn b_from_u(u: f64, su: f64) -> (f64, f64) {
    (8.0 * PI * PI * u, 8.0 * PI * PI * su)
}

/// U and its su from B and its su: U = B/(8π²), σ_U = σ_B/(8π²).
pub fn u_from_b(b: f64, su: f64) -> (f64, f64) {
    (b / (8.0 * PI * PI), su / (8.0 * PI * PI))
}

impl Structure {
    /// The scatterer of every atom for a radiation, in the order of `atoms`.
    pub fn scatterers(&self, radiation: Radiation) -> Result<Vec<Scatterer>, String> {
        self.atoms.iter().map(|a| Scatterer::lookup(&a.species, radiation)).collect()
    }

    /// sinθ/λ = 1/(2d) of a reflection.
    pub fn sin_theta_over_lambda(&self, hkl: [i32; 3]) -> f64 {
        0.5 * self.cell.inv_d_squared([hkl[0] as f64, hkl[1] as f64, hkl[2] as f64]).sqrt()
    }

    /// The structure factor F(h), with one scatterer per atom (from `scatterers`, or given).
    pub fn structure_factor(&self, scatterers: &[Scatterer], hkl: [i32; 3]) -> Complex {
        assert_eq!(scatterers.len(), self.atoms.len(), "one scatterer per atom");
        let s = self.sin_theta_over_lambda(hkl);
        let h = hkl.map(|v| v as f64);
        let mut total = Complex::ZERO;
        for (atom, sc) in self.atoms.iter().zip(scatterers) {
            let mut sum = Complex::ZERO;
            for op in &self.group.ops {
                // h·(R x + t), and h' = h·R for the betas of the image.
                let mut phase = 0.0;
                let mut hr = [0.0; 3];
                for ((hi, row), ti) in h.iter().zip(&op.r).zip(&op.t) {
                    let rx: f64 = row.iter().zip(&atom.x).map(|(r, x)| *r as f64 * x).sum();
                    phase += hi * (rx + *ti as f64 / DEN as f64);
                    for (v, r) in hr.iter_mut().zip(row) {
                        *v += hi * *r as f64;
                    }
                }
                let t = match atom.beta {
                    Some(b) => (-(b[0] * hr[0] * hr[0] + b[1] * hr[1] * hr[1] + b[2] * hr[2] * hr[2]
                        + 2.0 * (b[3] * hr[0] * hr[1] + b[4] * hr[0] * hr[2] + b[5] * hr[1] * hr[2])))
                        .exp(),
                    None => 1.0,
                };
                sum += Complex::cis(2.0 * PI * phase) * t;
            }
            total += sc.at(s) * sum * (atom.occ * (-atom.biso * s * s).exp());
        }
        total
    }
}

impl Structure {
    /// The powder intensity |F|² of reflection h: the average of |F(h)|² and |F(−h)|², which fall in the same peak
    /// (L251). They differ only when a scatterer has an imaginary part (f″, or the b″ of an absorbing nucleus) in a
    /// group without a centre of symmetry; otherwise this is |F(h)|² exactly. FullProf's own intensities of nabd4
    /// (P −4 2₁ c, B with b = 0.530 − 0.213i) are reproduced only so.
    pub fn intensity(&self, scatterers: &[Scatterer], hkl: [i32; 3]) -> f64 {
        let f = self.structure_factor(scatterers, hkl).norm_sqr();
        if scatterers.iter().all(|s| s.at(0.0).im == 0.0) {
            return f;
        }
        0.5 * (f + self.structure_factor(scatterers, hkl.map(|v| -v)).norm_sqr())
    }

    /// `intensity` with its derivatives: (|F|², d|F|²/ds, d|F|²/d(atom parameter) per atom), each the average over h and −h.
    pub fn intensity_derivatives(&self, scatterers: &[Scatterer], hkl: [i32; 3], s: f64) -> (f64, f64, Vec<[f64; 11]>) {
        let one = |fd: &FactorDerivatives| {
            let f2 = fd.f.norm_sqr();
            let ds = 2.0 * (fd.f.re * fd.ds.re + fd.f.im * fd.ds.im);
            let atoms: Vec<[f64; 11]> = fd.atoms.iter().map(|d| d.map(|x| 2.0 * (fd.f.re * x.re + fd.f.im * x.im))).collect();
            (f2, ds, atoms)
        };
        let a = one(&self.structure_factor_derivatives(scatterers, hkl, s));
        if scatterers.iter().all(|sc| sc.at(0.0).im == 0.0) {
            return a;
        }
        let b = one(&self.structure_factor_derivatives(scatterers, hkl.map(|v| -v), s));
        let atoms = a.2.iter().zip(&b.2).map(|(x, y)| std::array::from_fn(|k| 0.5 * (x[k] + y[k]))).collect();
        (0.5 * (a.0 + b.0), 0.5 * (a.1 + b.1), atoms)
    }
}

/// The structure factor and its analytic derivatives (L053): with respect to s = sinθ/λ, and for each atom with respect
/// to x, y, z, Biso, Occ and β11 … β23 (in that order; the betas' derivatives are zero for an isotropic atom).
#[derive(Clone, Debug)]
pub struct FactorDerivatives {
    pub f: Complex,
    pub ds: Complex,
    pub atoms: Vec<[Complex; 11]>,
}

impl Structure {
    /// F(h) at s with its derivatives (`FactorDerivatives`). Agrees with `structure_factor` when s is the reflection's.
    pub fn structure_factor_derivatives(&self, scatterers: &[Scatterer], hkl: [i32; 3], s: f64) -> FactorDerivatives {
        let h = hkl.map(|v| v as f64);
        let (mut f, mut fds) = (Complex::ZERO, Complex::ZERO);
        let mut atoms = Vec::with_capacity(self.atoms.len());
        for (atom, sc) in self.atoms.iter().zip(scatterers) {
            let mut sum = Complex::ZERO;
            let mut dx = [Complex::ZERO; 3];
            let mut db = [Complex::ZERO; 6];
            for op in &self.group.ops {
                let mut phase = 0.0;
                let mut hr = [0.0; 3];
                for ((hi, row), ti) in h.iter().zip(&op.r).zip(&op.t) {
                    let rx: f64 = row.iter().zip(&atom.x).map(|(r, x)| *r as f64 * x).sum();
                    phase += hi * (rx + *ti as f64 / DEN as f64);
                    for (v, r) in hr.iter_mut().zip(row) {
                        *v += hi * *r as f64;
                    }
                }
                let q = [hr[0] * hr[0], hr[1] * hr[1], hr[2] * hr[2], 2.0 * hr[0] * hr[1], 2.0 * hr[0] * hr[2], 2.0 * hr[1] * hr[2]];
                let t = match atom.beta {
                    Some(b) => (-(0..6).map(|k| b[k] * q[k]).sum::<f64>()).exp(),
                    None => 1.0,
                };
                let e = Complex::cis(2.0 * PI * phase) * t;
                sum += e;
                // d/dx_k of exp(2πi h·(Rx + t)) = 2πi (hR)_k exp(…).
                for k in 0..3 {
                    dx[k] += Complex::new(-e.im, e.re) * (2.0 * PI * hr[k]);
                }
                if atom.beta.is_some() {
                    for k in 0..6 {
                        db[k] += e * (-q[k]);
                    }
                }
            }
            let fs = sc.at(s);
            let dw = (-atom.biso * s * s).exp();
            let common = fs * (atom.occ * dw);
            let term = common * sum;
            f += term;
            fds += (Complex::new(sc.ds(s), 0.0) * dw + fs * (-2.0 * atom.biso * s * dw)) * sum * atom.occ;
            let mut d = [Complex::ZERO; 11];
            for k in 0..3 {
                d[k] = common * dx[k];
            }
            d[3] = term * (-s * s);
            d[4] = fs * dw * sum;
            for k in 0..6 {
                d[5 + k] = common * db[k];
            }
            atoms.push(d);
        }
        FactorDerivatives { f, ds: fds, atoms }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn atom(species: &str, x: [f64; 3], occ: f64) -> Atom {
        Atom { label: species.into(), species: species.into(), x, biso: 0.0, occ, beta: None, charge: None }
    }

    #[test]
    fn rock_salt_by_hand() {
        // NaCl-type with neutron lengths b1, b2: F(111) = 4(b1 − b2), F(200) = 4(b1 + b2) for full occupancy, Occ = 4/192.
        let (g, _) = SpaceGroup::from_symbol("F m -3 m").unwrap();
        let cell = Cell::new(5.64, 5.64, 5.64, 90.0, 90.0, 90.0).unwrap();
        let o = occupancy_from_chemical(&g, [0.0; 3], 1.0).unwrap();
        assert!((o - 4.0 / 192.0).abs() < 1e-15);
        let st = Structure { group: g, cell, atoms: vec![atom("Na", [0.0; 3], o), atom("Cl", [0.5; 3], o)] };
        let sc = st.scatterers(Radiation::Neutron).unwrap();
        let (b1, b2) = (sc[0].at(0.0).re, sc[1].at(0.0).re);
        let f111 = st.structure_factor(&sc, [1, 1, 1]);
        let f200 = st.structure_factor(&sc, [2, 0, 0]);
        assert!((f111.re - 4.0 * (b1 - b2)).abs() < 1e-12 && f111.im.abs() < 1e-12);
        assert!((f200.re - 4.0 * (b1 + b2)).abs() < 1e-12);
        assert!(st.structure_factor(&sc, [1, 0, 0]).abs() < 1e-12);
    }

    #[test]
    fn beta_of_isotropic_u() {
        // An isotropic U on a cubic cell: β11 = 2π²U/a², and exp(−hβh) = exp(−B s²) with B = 8π²U.
        let cell = Cell::new(4.0, 4.0, 4.0, 90.0, 90.0, 90.0).unwrap();
        let u = 0.01;
        let b = beta_from_u(&cell, [u, u, u, 0.0, 0.0, 0.0]);
        let h = [1.0, 2.0, 3.0];
        let hbh = b[0] * h[0] * h[0] + b[1] * h[1] * h[1] + b[2] * h[2] * h[2];
        let s2 = cell.inv_d_squared(h) / 4.0;
        assert!((hbh - biso_from_uiso(u) * s2).abs() < 1e-14);
    }
}
