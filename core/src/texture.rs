//! Preferred orientation along several axes (P30, L360): the multi-axial March–Dollase model, FullProf's Nor = 2 and 3
//! (FullProf's News of 12 October 2003, `manuals/fullprof/news/Fullprof_News_2003.htm`).
//!
//! A phase's crystallites are a mixture: a fraction f_j of them has its axis j (a reciprocal-lattice direction)
//! distributed about the sample's axis by March's density with ratio r_j (March 1932; Dollase 1986, Table 1, no. 6):
//!
//!   D(r, ρ) = (r² cos²ρ + sin²ρ/r)^(−3/2),  P_h = Σ_j f_j D_j,  Σ_j f_j = 1,
//!
//! ρ the angle between the axis and the sample's axis, averaged over the reflections equivalent to h as the one-axis
//! function is (`Model::preferred_orientation`). The last fraction is not a value of its own: it is 1 minus the others
//! (the News: "The fraction corresponding to the last POD cannot be refined").
//!
//! - A flat sample in reflection (Nor = 2): the scattering vector lies along the sample's axis, ρ is the angle α between
//!   h and axis j. FullProf's G2 (Pref2) does not enter: yb_nor2_g2 (G2 = 0.3) has yb_nor2's CORR (L360).
//! - An anvil cell with the beam along the force and the image plate perpendicular to it (Nor = 3; Wright, Nelmes,
//!   Belmonte & McMahon (1996), J. Synchrotron Rad. 3, 112–119, eqs. 4–6, which the News names; the paper is not in
//!   `manuals/papers`): the scattering vector makes ψ = 90° − θ with the sample's axis, so axis j, at α from h, lies on a
//!   cone about the scattering vector, cos ρ = cos α cos ψ + sin α sin ψ cos φ, and D_j is averaged over φ with the
//!   `steps` of the integration. That form is believed, not checked: FullProf's yb_nor3 CORR is up to 3.6 % below it
//!   (test `texture`, frozen and listed for Josu, L360).
//!
//! The values are numbered for `Field::Texture(p, k)`: r_j at 2j, f_j at 2j + 1, the last axis's r at 2n − 2 (no
//! fraction of its own).

use crate::dual::Real;

/// One axis of the model: its direction (h k l, a reciprocal-lattice vector), March's r and the fraction.
#[derive(Clone, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Axis {
    pub direction: [f64; 3],
    pub r: f64,
    pub fraction: f64,
}

/// The multi-axial model of a phase.
#[derive(Clone, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MultiAxial {
    pub axes: Vec<Axis>,
    /// 0: a flat sample in reflection (Nor = 2). Above 0: an anvil cell (Nor = 3), the steps of the integration about
    /// each reflection's cone (FullProf's NOrstp).
    pub steps: usize,
}

/// March–Dollase's density at the angle ρ to the axis, given cos ρ: (r² cos²ρ + sin²ρ/r)^(−3/2), in multiples of a
/// random distribution (MRD; its mean over the sphere is 1). For Publish's F10 (P48.4); against ρ itself,
/// `march_dollase_at` and `march_dollase_curve` (P59).
pub fn march_dollase(r: f64, cos_rho: f64) -> f64 {
    march_g(r, cos_rho)
}

/// March–Dollase's density at ρ degrees from the axis (Publish's P59, L374).
pub fn march_dollase_at(r: f64, rho_degrees: f64) -> f64 {
    march_g(r, rho_degrees.to_radians().cos())
}

/// March–Dollase's density against the angle, for Publish's F10 (P59, L374): `n` + 1 points (ρ in degrees, D) at ρ =
/// 0, 90/n, …, 90°, from D(0) = r⁻³ along the axis to D(90°) = r^{3/2} across it.
pub fn march_dollase_curve(r: f64, n: usize) -> Vec<(f64, f64)> {
    let n = n.max(1);
    (0..=n).map(|k| 90.0 * k as f64 / n as f64).map(|rho| (rho, march_dollase_at(r, rho))).collect()
}

fn march_g<T: Real>(r: T, c: T) -> T {
    let c2 = c * c;
    let c2 = if c2.val() > 1.0 { T::cst(1.0) } else { c2 };
    (r * r * c2 + (-c2 + 1.0) / r).powf(-1.5)
}

/// dD/dr at cos ρ.
fn march_dr(r: f64, c: f64) -> f64 {
    let c2 = (c * c).min(1.0);
    let s2 = 1.0 - c2;
    -1.5 * (r * r * c2 + s2 / r).powf(-2.5) * (2.0 * r * c2 - s2 / (r * r))
}

/// cos ρ at the cone's steps for cos α and ψ (one value, cos α, when ψ is 0 or there are no steps).
fn cone<T: Real>(cos_alpha: T, psi: f64, steps: usize) -> Vec<T> {
    if steps == 0 || psi == 0.0 {
        return vec![cos_alpha];
    }
    let c2 = cos_alpha * cos_alpha;
    let sin_alpha = if c2.val() >= 1.0 { T::cst(0.0) } else { (-c2 + 1.0).sqrt() };
    (0..steps).map(|k| cos_alpha * psi.cos() + sin_alpha * (psi.sin() * (std::f64::consts::PI * (k as f64 + 0.5) / steps as f64).cos())).collect()
}

impl MultiAxial {
    /// The number of values (`Field::Texture`).
    pub fn len(&self) -> usize {
        (2 * self.axes.len()).saturating_sub(1)
    }

    pub fn is_empty(&self) -> bool {
        self.axes.is_empty()
    }

    /// What value k is: (axis, 0 for r or 1 for the fraction).
    pub fn part(&self, k: usize) -> Option<(usize, usize)> {
        (k < self.len()).then_some((k / 2, k % 2))
    }

    pub fn value(&self, k: usize) -> f64 {
        match self.part(k) {
            Some((j, 0)) => self.axes[j].r,
            Some((j, _)) => self.axes[j].fraction,
            None => 0.0,
        }
    }

    /// Sets value k; a fraction moves the last axis's, which keeps the sum 1.
    pub fn set(&mut self, k: usize, v: f64) {
        match self.part(k) {
            Some((j, 0)) => self.axes[j].r = v,
            Some((j, _)) => {
                self.axes[j].fraction = v;
                self.close();
            }
            None => {}
        }
    }

    /// The last axis's fraction: 1 minus the others.
    pub fn close(&mut self) {
        let n = self.axes.len();
        if n > 0 {
            let others: f64 = self.axes[..n - 1].iter().map(|a| a.fraction).sum();
            self.axes[n - 1].fraction = 1.0 - others;
        }
    }

    /// The dictionary id and the component (the axis, from 1) of value k.
    pub fn id(&self, k: usize) -> (&'static str, String) {
        match self.part(k) {
            Some((j, 1)) => ("po_axis_fraction", (j + 1).to_string()),
            Some((j, _)) => ("po_axis_ratio", (j + 1).to_string()),
            None => ("po_axis_ratio", (k / 2 + 1).to_string()),
        }
    }

    /// ψ, the angle between the scattering vector and the sample's axis, at Bragg angle θ (radians): 0 in reflection,
    /// 90° − θ in the anvil cell.
    pub fn psi(&self, theta: f64) -> f64 {
        if self.steps == 0 {
            0.0
        } else {
            std::f64::consts::FRAC_PI_2 - theta
        }
    }

    /// Axis j's density for a reflection at cos α from it, at ψ.
    pub fn density(&self, j: usize, cos_alpha: f64, psi: f64) -> f64 {
        self.density_g(j, cos_alpha, psi)
    }

    fn density_g<T: Real>(&self, j: usize, cos_alpha: T, psi: f64) -> T {
        let c = cone(cos_alpha, psi, self.steps);
        let n = c.len() as f64;
        c.into_iter().fold(T::cst(0.0), |s, x| s + march_g(T::cst(self.axes[j].r), x)) / n
    }

    /// The factor for one reflection, given cos α to each axis: Σ_j f_j D_j.
    pub fn factor_g<T: Real>(&self, cos_alpha: &[T], psi: f64) -> T {
        (0..self.axes.len()).fold(T::cst(0.0), |s, j| s + self.density_g(j, cos_alpha[j], psi) * self.axes[j].fraction)
    }

    pub fn factor(&self, cos_alpha: &[f64], psi: f64) -> f64 {
        self.factor_g(cos_alpha, psi)
    }

    /// The derivative of `factor` with respect to value k: f_j dD_j/dr_j, or D_j − D_last for a fraction.
    pub fn derivative(&self, k: usize, cos_alpha: &[f64], psi: f64) -> f64 {
        match self.part(k) {
            Some((j, 0)) => {
                let c = cone(cos_alpha[j], psi, self.steps);
                let n = c.len() as f64;
                self.axes[j].fraction * c.into_iter().map(|x| march_dr(self.axes[j].r, x)).sum::<f64>() / n
            }
            Some((j, _)) => {
                let last = self.axes.len() - 1;
                self.density(j, cos_alpha[j], psi) - self.density(last, cos_alpha[last], psi)
            }
            None => 0.0,
        }
    }
}

/// A phase's preferred orientation as axes, for Publish's F10 (P48.4): its texture, or, under March–Dollase
/// (`Model::march`), its one axis with G1 as r and the untextured fraction G2 as a second axis of r = 1 (D = 1
/// everywhere): G2 + (1 − G2) D(G1) = Σ f_j D_j. None under the exponential function.
pub fn of_phase(ph: &crate::powder::Phase, march: bool) -> Option<MultiAxial> {
    if let Some(t) = &ph.texture {
        return Some(t.clone());
    }
    if !march {
        return None;
    }
    let mut axes = vec![Axis { direction: ph.pref_direction, r: ph.pref[0], fraction: 1.0 - ph.pref[1] }];
    if ph.pref[1] != 0.0 {
        axes.push(Axis { direction: ph.pref_direction, r: 1.0, fraction: ph.pref[1] });
    }
    Some(MultiAxial { axes, steps: 0 })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn model(steps: usize) -> MultiAxial {
        MultiAxial { axes: vec![Axis { direction: [1.0, 1.0, 0.0], r: 0.85, fraction: 0.6 }, Axis { direction: [0.0, 0.0, 1.0], r: 1.2, fraction: 0.4 }], steps }
    }

    #[test]
    fn derivatives_by_differences() {
        for steps in [0, 20] {
            let t = model(steps);
            for (c, psi) in [([0.3, 0.9], 0.0), ([0.7, 0.1], 0.6), ([1.0, 0.0], 1.2)] {
                for k in 0..t.len() {
                    let h = 1e-6;
                    let (mut a, mut b) = (t.clone(), t.clone());
                    a.set(k, t.value(k) + h);
                    b.set(k, t.value(k) - h);
                    let num = (a.factor(&c, psi) - b.factor(&c, psi)) / (2.0 * h);
                    assert!((num - t.derivative(k, &c, psi)).abs() < 1e-6, "k {k}: {num} {}", t.derivative(k, &c, psi));
                }
            }
        }
    }

    #[test]
    fn march_dollase_is_normalised_on_the_sphere() {
        // ∫ D sin ρ dρ over [0, π/2] is 1 for every r (Dollase 1986, p. 268).
        for r in [0.5, 0.85, 1.0, 1.2, 2.0] {
            let n = 20000;
            let s: f64 = (0..n).map(|k| (k as f64 + 0.5) / n as f64).map(|u| march_dollase(r, u)).sum::<f64>() / n as f64;
            assert!((s - 1.0).abs() < 1e-6, "r {r}: {s}");
        }
    }

    #[test]
    fn the_curve_against_the_angle() {
        // P59: the ends r⁻³ and r^{3/2}, and ∫ D sin ρ dρ = 1 over the curve (Simpson's rule on 90 intervals).
        for r in [0.6, 1.0, 1.4] {
            let c = march_dollase_curve(r, 90);
            assert_eq!(c.len(), 91);
            assert!((c[0].1 - r.powi(-3)).abs() < 1e-12 && (c[90].1 - r.powf(1.5)).abs() < 1e-12 && c[90].0 == 90.0);
            let h = (1.0f64).to_radians();
            let f = |k: usize| c[k].1 * c[k].0.to_radians().sin();
            let s: f64 = h / 3.0 * (f(0) + f(90) + (1..90).map(|k| if k % 2 == 1 { 4.0 * f(k) } else { 2.0 * f(k) }).sum::<f64>());
            assert!((s - 1.0).abs() < 1e-4, "r {r}: {s}");
            assert!((march_dollase_at(r, 37.0) - march_dollase(r, 37f64.to_radians().cos())).abs() < 1e-15);
        }
    }

    #[test]
    fn the_cone_keeps_the_pole_figure_normalised() {
        // At any ψ, the density averaged over the axis's directions on the sphere is still 1 (the cone is a rotation).
        let t = MultiAxial { axes: vec![Axis { direction: [0.0, 0.0, 1.0], r: 0.7, fraction: 1.0 }], steps: 200 };
        for psi in [0.3, 0.9, 1.4] {
            let n = 4000;
            let s: f64 = (0..n).map(|k| (k as f64 + 0.5) / n as f64).map(|u| t.density(0, u, psi)).sum::<f64>() / n as f64;
            assert!((s - 1.0).abs() < 1e-4, "ψ {psi}: {s}");
        }
    }
}
