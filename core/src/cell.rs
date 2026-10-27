//! The unit cell: direct and reciprocal parameters, metric tensors, volume, d-spacings, and the
//! Cartesian frame FullProf uses (x along a, y in the ab-plane, z = x × y).
//!
//! Lengths in Å, angles in degrees, reciprocal quantities without the 2π factor (a* = 1/a for an
//! orthogonal cell), as in FullProf's `.out` files and the International Tables, Vol. B.

use crate::mat3::{self, Mat3, Vec3};
use std::fmt;

/// Direct-cell parameters a, b, c (Å) and α, β, γ (degrees).
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Cell {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub alpha: f64,
    pub beta: f64,
    pub gamma: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum CellError {
    /// A length that is not a positive finite number.
    Length(&'static str, f64),
    /// An angle outside 0° < angle < 180°.
    Angle(&'static str, f64),
    /// Angles that cannot close a cell (for instance α + β < γ): the volume would be imaginary.
    Angles { alpha: f64, beta: f64, gamma: f64 },
}

impl fmt::Display for CellError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            CellError::Length(n, v) => {
                write!(f, "cell length {n} = {v} must be a positive number (Å)")
            }
            CellError::Angle(n, v) => write!(
                f,
                "cell angle {n} = {v} must lie strictly between 0° and 180°"
            ),
            CellError::Angles { alpha, beta, gamma } => {
                write!(
                    f,
                    "the angles α = {alpha}°, β = {beta}°, γ = {gamma}° do not close a cell"
                )
            }
        }
    }
}

impl std::error::Error for CellError {}

impl Cell {
    pub fn new(
        a: f64,
        b: f64,
        c: f64,
        alpha: f64,
        beta: f64,
        gamma: f64,
    ) -> Result<Cell, CellError> {
        for (n, v) in [("a", a), ("b", b), ("c", c)] {
            if !(v.is_finite() && v > 0.0) {
                return Err(CellError::Length(n, v));
            }
        }
        for (n, v) in [("α", alpha), ("β", beta), ("γ", gamma)] {
            if !(v.is_finite() && v > 0.0 && v < 180.0) {
                return Err(CellError::Angle(n, v));
            }
        }
        let cell = Cell {
            a,
            b,
            c,
            alpha,
            beta,
            gamma,
        };
        if cell.volume_factor_squared() <= 1e-12 {
            return Err(CellError::Angles { alpha, beta, gamma });
        }
        Ok(cell)
    }

    fn cosines(&self) -> Vec3 {
        [self.alpha, self.beta, self.gamma].map(|x| x.to_radians().cos())
    }

    /// 1 − cos²α − cos²β − cos²γ + 2 cosα cosβ cosγ, so that V = abc √(this).
    fn volume_factor_squared(&self) -> f64 {
        let [ca, cb, cg] = self.cosines();
        1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg
    }

    /// Direct metric tensor G, with G_ij = a_i · a_j.
    pub fn metric(&self) -> Mat3 {
        let [ca, cb, cg] = self.cosines();
        let (a, b, c) = (self.a, self.b, self.c);
        [
            [a * a, a * b * cg, a * c * cb],
            [a * b * cg, b * b, b * c * ca],
            [a * c * cb, b * c * ca, c * c],
        ]
    }

    /// Cell volume in Å³.
    pub fn volume(&self) -> f64 {
        self.a * self.b * self.c * self.volume_factor_squared().sqrt()
    }

    /// Reciprocal metric tensor G* = G⁻¹ (Å⁻²).
    pub fn reciprocal_metric(&self) -> Mat3 {
        mat3::inverse(&self.metric()).expect("a valid cell has a non-singular metric")
    }

    /// The reciprocal cell a*, b*, c* (Å⁻¹) and α*, β*, γ* (degrees).
    pub fn reciprocal(&self) -> Cell {
        let g = self.reciprocal_metric();
        let (a, b, c) = (g[0][0].sqrt(), g[1][1].sqrt(), g[2][2].sqrt());
        let ang = |gij: f64, x: f64, y: f64| (gij / (x * y)).clamp(-1.0, 1.0).acos().to_degrees();
        Cell {
            a,
            b,
            c,
            alpha: ang(g[1][2], b, c),
            beta: ang(g[0][2], a, c),
            gamma: ang(g[0][1], a, b),
        }
    }

    /// 1/d² = hᵀ G* h for the reflection (h, k, l).
    pub fn inv_d_squared(&self, hkl: Vec3) -> f64 {
        mat3::quad(&self.reciprocal_metric(), &hkl)
    }

    /// The interplanar spacing d(hkl) in Å; `None` for (0, 0, 0).
    pub fn d_spacing(&self, hkl: Vec3) -> Option<f64> {
        let s2 = self.inv_d_squared(hkl);
        (s2 > 0.0).then(|| 1.0 / s2.sqrt())
    }

    /// The Bragg angle 2θ in degrees for wavelength λ (Å); `None` if the reflection is beyond
    /// the limiting sphere (λ / 2d > 1) or is (0, 0, 0).
    pub fn two_theta(&self, hkl: Vec3, lambda: f64) -> Option<f64> {
        let s = lambda / (2.0 * self.d_spacing(hkl)?);
        // a rounding margin, so that back-scattering at exactly 2θ = 180° is not refused
        (s <= 1.0 + 1e-12).then(|| 2.0 * s.min(1.0).asin().to_degrees())
    }

    /// Crystal-to-orthonormal matrix M (FullProf's `Cr_Orth_cel`): its columns are a, b, c in the
    /// Cartesian frame x ∥ a, y in the ab-plane, z = x × y, so that r_cart = M · r_frac.
    pub fn crystal_to_orthonormal(&self) -> Mat3 {
        let [ca, cb, cg] = self.cosines();
        let sg = self.gamma.to_radians().sin();
        let (a, b, c) = (self.a, self.b, self.c);
        [
            [a, b * cg, c * cb],
            [0.0, b * sg, c * (ca - cb * cg) / sg],
            [0.0, 0.0, self.volume() / (a * b * sg)],
        ]
    }

    /// The cell's twelve edges, each from one corner to another, in Cartesian Å (`crystal_to_orthonormal`): the box
    /// drawn to scale for any phase, one with atoms or a Le Bail or profile phase (Explore's E22). The edges of
    /// `geometry::scene` are these.
    pub fn edges(&self) -> Vec<(Vec3, Vec3)> {
        const CORNERS: [([f64; 3], [f64; 3]); 12] = [
            ([0., 0., 0.], [1., 0., 0.]),
            ([0., 0., 0.], [0., 1., 0.]),
            ([0., 0., 0.], [0., 0., 1.]),
            ([1., 0., 0.], [1., 1., 0.]),
            ([1., 0., 0.], [1., 0., 1.]),
            ([0., 1., 0.], [1., 1., 0.]),
            ([0., 1., 0.], [0., 1., 1.]),
            ([0., 0., 1.], [1., 0., 1.]),
            ([0., 0., 1.], [0., 1., 1.]),
            ([1., 1., 0.], [1., 1., 1.]),
            ([1., 0., 1.], [1., 1., 1.]),
            ([0., 1., 1.], [1., 1., 1.]),
        ];
        let m = self.crystal_to_orthonormal();
        CORNERS.iter().map(|(a, b)| (mat3::mul_vec(&m, a), mat3::mul_vec(&m, b))).collect()
    }

    /// Orthonormal-to-crystal matrix M⁻¹ (FullProf's `Orth_Cr_cel`).
    pub fn orthonormal_to_crystal(&self) -> Mat3 {
        mat3::inverse(&self.crystal_to_orthonormal())
            .expect("a valid cell has a non-singular frame")
    }

    /// Busing–Levy B matrix (Acta Cryst. 22 (1967) 457; FullProf's `BL_M`): h_cart = B · h, with
    /// B upper triangular, x ∥ a*, y in the a*b*-plane.
    pub fn busing_levy_b(&self) -> Mat3 {
        let r = self.reciprocal();
        let (cbs, cgs) = (r.beta.to_radians().cos(), r.gamma.to_radians().cos());
        let (sbs, sgs) = (r.beta.to_radians().sin(), r.gamma.to_radians().sin());
        let ca = self.alpha.to_radians().cos();
        [
            [r.a, r.b * cgs, r.c * cbs],
            [0.0, r.b * sgs, -r.c * sbs * ca],
            [0.0, 0.0, 1.0 / self.c],
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(x: f64, y: f64, tol: f64) -> bool {
        (x - y).abs() <= tol
    }

    #[test]
    fn cubic_cell() {
        let c = Cell::new(4.0, 4.0, 4.0, 90.0, 90.0, 90.0).unwrap();
        assert!(close(c.volume(), 64.0, 1e-12));
        assert!(close(
            c.d_spacing([1.0, 1.0, 1.0]).unwrap(),
            4.0 / 3f64.sqrt(),
            1e-12
        ));
        assert!(close(c.reciprocal().a, 0.25, 1e-12));
    }

    #[test]
    fn hexagonal_d100_is_a_sqrt3_over_2() {
        let c = Cell::new(3.0, 3.0, 5.0, 90.0, 90.0, 120.0).unwrap();
        assert!(close(
            c.d_spacing([1.0, 0.0, 0.0]).unwrap(),
            3.0 * 3f64.sqrt() / 2.0,
            1e-12
        ));
        assert!(close(c.reciprocal().gamma, 60.0, 1e-9));
    }

    #[test]
    fn triclinic_volume_matches_metric_determinant_and_frame() {
        let c = Cell::new(5.1, 6.3, 7.7, 81.0, 97.5, 104.2).unwrap();
        let v = c.volume();
        assert!(close(mat3::det(&c.metric()).sqrt(), v, 1e-9 * v));
        assert!(close(mat3::det(&c.crystal_to_orthonormal()), v, 1e-9 * v));
        // the reciprocal of the reciprocal is the direct cell
        let rr = c.reciprocal().reciprocal();
        for (x, y) in [
            (rr.a, c.a),
            (rr.b, c.b),
            (rr.c, c.c),
            (rr.alpha, c.alpha),
            (rr.beta, c.beta),
            (rr.gamma, c.gamma),
        ] {
            assert!(close(x, y, 1e-9), "{rr:?} vs {c:?}");
        }
        // |B h| = 1/d for any h
        let h = [2.0, -1.0, 3.0];
        let bh = mat3::mul_vec(&c.busing_levy_b(), &h);
        let len = (bh[0] * bh[0] + bh[1] * bh[1] + bh[2] * bh[2]).sqrt();
        assert!(close(len, 1.0 / c.d_spacing(h).unwrap(), 1e-12));
    }

    #[test]
    fn two_theta_exact_cases_and_beyond_the_limiting_sphere() {
        // λ = a√2 puts (100) at 2θ = 90° exactly (sin θ = λ / 2a = 1/√2); λ = 2a puts it at 180°.
        let si = Cell::new(5.43102, 5.43102, 5.43102, 90.0, 90.0, 90.0).unwrap();
        assert!(close(
            si.two_theta([1.0, 0.0, 0.0], 5.43102 * 2f64.sqrt())
                .unwrap(),
            90.0,
            1e-9
        ));
        // asin is ill-conditioned at 1: a rounding of 1e-16 in sin θ moves 2θ by ~1e-6°
        assert!(close(
            si.two_theta([1.0, 0.0, 0.0], 2.0 * 5.43102).unwrap(),
            180.0,
            1e-5
        ));
        assert!(si.two_theta([20.0, 0.0, 0.0], 1.540593).is_none());
        assert!(si.d_spacing([0.0, 0.0, 0.0]).is_none());
    }

    #[test]
    fn impossible_cells_are_refused() {
        assert!(matches!(
            Cell::new(-1.0, 1.0, 1.0, 90.0, 90.0, 90.0),
            Err(CellError::Length("a", _))
        ));
        assert!(matches!(
            Cell::new(1.0, 1.0, 1.0, 180.0, 90.0, 90.0),
            Err(CellError::Angle("α", _))
        ));
        assert!(matches!(
            Cell::new(1.0, 1.0, 1.0, 30.0, 30.0, 100.0),
            Err(CellError::Angles { .. })
        ));
    }
}
