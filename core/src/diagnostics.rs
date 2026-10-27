//! What the app's Refine and Reflections views draw beside the pattern (W7–W8; app requests 19–21, L100): each
//! reflection's width, calculated and observed intensity; each phase's share of the calculated pattern at a point; the
//! cumulative χ² and the weighted difference along the pattern.
//!
//! **I obs is partitioned by the model** (Rietveld 1969; FullProf's "Iobs", manual p. 11): the measured counts above
//! the background are shared among the reflections in proportion to their calculated profiles,
//!
//!   I_obs,h = I_calc,h · Σ_i Ω_h,i (y_i − b_i)/(y_c,i − b_i) / Σ_i Ω_h,i,
//!
//! the sums over the points of the reflection's window that count in χ² and where y_c − b > 0; Ω_h,i its calculated
//! profile at point i. It is not a measurement of the reflection: where reflections overlap, it repeats the model's
//! ratios. Its σ is the propagation of the points' σ_i through the same sum,
//!
//!   σ(I_obs,h) = I_calc,h / Σ_i Ω_h,i · √(Σ_i (Ω_h,i σ_i / (y_c,i − b_i))²),
//!
//! and the residual (I_obs − I_calc)/σ. FullProf's `.out` prints Iobs and Sigma per reflection; the test
//! `reflection_diagnostics` compares them.

use crate::data::Pattern;
use crate::powder::{Calculated, Intensities, Model};

/// A reflection of a calculated pattern, as the Reflections view lists it.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ReflectionRow {
    pub phase: usize,
    pub hkl: [i32; 3],
    pub multiplicity: usize,
    /// 1 or 2.
    pub wavelength: usize,
    /// 2θ of the peak (degrees), Zero and the displacement terms included.
    pub two_theta: f64,
    pub d: f64,
    /// The profile's FWHM H (degrees) and η.
    pub fwhm: f64,
    pub eta: f64,
    /// |F|² of the phase's structure (None for intensities given as a list).
    pub f2: Option<f64>,
    /// The integrated intensity, scale and corrections included.
    pub icalc: f64,
    /// I obs partitioned by the model, its σ, and (I obs − I calc)/σ (0 where σ is 0: no point counts).
    pub iobs: f64,
    pub sigma: f64,
    pub residual: f64,
}

/// The reflections of the model's pattern at the data's points, with I obs partitioned by the model.
pub fn reflections(model: &Model, data: &Pattern) -> Vec<ReflectionRow> {
    let c = model.calculate(&data.x);
    let w = crate::lsq::weights(model, data, &c.ycalc);
    let excluded = model.excluded_points(&data.x);
    c.peaks
        .iter()
        .map(|pk| {
            let (mut omega, mut part, mut var) = (0.0, 0.0, 0.0);
            model.peak_profile(&data.x, &excluded, pk, |i, v| {
                let net = c.ycalc[i] - c.background[i];
                if net > 0.0 && w[i] > 0.0 {
                    omega += v;
                    part += v * (data.y[i] - c.background[i]) / net;
                    var += (v * data.sigma[i] / net).powi(2);
                }
            });
            let ph = &model.phases[pk.phase];
            let (iobs, sigma) = if omega > 0.0 { (pk.intensity * part / omega, pk.intensity / omega * var.sqrt()) } else { (0.0, 0.0) };
            let f2 = match &ph.intensities {
                Intensities::Structure(st, sc) => Some(st.intensity(sc, pk.hkl)),
                Intensities::List(_) => None,
            };
            ReflectionRow {
                phase: pk.phase,
                hkl: pk.hkl,
                multiplicity: pk.multiplicity,
                wavelength: pk.wavelength,
                two_theta: pk.two_theta,
                d: ph.cell.d_spacing(pk.hkl.map(|v| v as f64)).unwrap_or(f64::NAN),
                fwhm: pk.fwhm,
                eta: pk.eta,
                f2,
                icalc: pk.intensity,
                iobs,
                sigma,
                residual: if sigma > 0.0 { (iobs - pk.intensity) / sigma } else { 0.0 },
            }
        })
        .collect()
}

/// What the calculated pattern is made of at a point: the background's share and each phase's, summing to 1.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Shares {
    pub ycalc: f64,
    pub background: f64,
    pub phases: Vec<f64>,
}

/// The shares at point i of a calculated pattern (all 0 where Ycal is 0).
pub fn shares(c: &Calculated, i: usize) -> Shares {
    let y = c.ycalc[i];
    let f = |v: f64| if y != 0.0 { v / y } else { 0.0 };
    Shares { ycalc: y, background: f(c.background[i]), phases: c.per_phase.iter().map(|p| f(p[i])).collect() }
}

/// A reflection's profile on the background, point by point (Explore's peak-shape view; L148).
pub fn on_background(peak: &[f64], background: &[f64]) -> Vec<f64> {
    peak.iter().zip(background).map(|(p, b)| p + b).collect()
}

/// The strongest of some reflections: the index of the largest I calc (Explore's default reflection; L148).
pub fn strongest(icalc: &[f64]) -> Option<usize> {
    (0..icalc.len()).filter(|&i| icalc[i].is_finite()).max_by(|&a, &b| icalc[a].total_cmp(&icalc[b]))
}

/// Yobs − Ycal at each point (the Rietveld plot's difference strip; L146).
pub fn difference(yobs: &[f64], ycalc: &[f64]) -> Vec<f64> {
    yobs.iter().zip(ycalc).map(|(o, c)| o - c).collect()
}

impl Shares {
    /// The background's and each phase's share in per cent ("who is this peak"; L146).
    pub fn in_percent(&self) -> (f64, Vec<f64>) {
        (100.0 * self.background, self.phases.iter().map(|p| 100.0 * p).collect())
    }
}

/// The data point nearest 2θ = t.
pub fn nearest_point(x: &[f64], t: f64) -> Option<usize> {
    (0..x.len()).min_by(|&a, &b| (x[a] - t).abs().total_cmp(&(x[b] - t).abs()))
}

/// The shares at the data point nearest 2θ = t ("who is this peak", W7).
pub fn shares_at(model: &Model, data: &Pattern, t: f64) -> Option<Shares> {
    let i = nearest_point(&data.x, t)?;
    Some(shares(&model.calculate(&data.x), i))
}

/// How far a reflection's peak reaches, in its FWHM ("who is this peak": the reflections whose peaks reach a point).
pub const REACH_FWHM: f64 = 1.5;
/// The reach of a reflection whose width is not known (degrees).
pub const REACH_UNKNOWN: f64 = 0.05;

/// Whether the peak of a reflection at `two_theta` with `fwhm` reaches 2θ = t: within 1.5 FWHM, or 0.05° when the
/// width is not known (L145).
pub fn reaches(two_theta: f64, fwhm: f64, t: f64) -> bool {
    (two_theta - t).abs() <= REACH_FWHM * if fwhm.is_finite() && fwhm > 0.0 { fwhm } else { REACH_UNKNOWN }
}

#[cfg(test)]
mod reach_tests {
    /// L146: the difference; the cumulative χ² on the counts' scale; the shares in per cent.
    #[test]
    fn the_difference_the_cumulative_on_counts_and_shares_in_percent() {
        assert_eq!(super::difference(&[10.0, 5.0], &[8.0, 6.0]), vec![2.0, -1.0]);
        let c = super::Cumulative { sum: vec![0.0, 1.0, 4.0], weighted_difference: vec![0.0; 3] };
        assert_eq!(c.on_counts(&[100.0, 50.0, 20.0], &[90.0, 200.0, 10.0]), vec![0.0, 50.0, 200.0]);
        let s = super::Shares { ycalc: 10.0, background: 0.25, phases: vec![0.5, 0.25] };
        assert_eq!(s.in_percent(), (25.0, vec![50.0, 25.0]));
    }

    /// L148: a profile on its background; the strongest reflection.
    #[test]
    fn a_profile_on_its_background_and_the_strongest_reflection() {
        assert_eq!(super::on_background(&[1.0, 5.0], &[10.0, 10.0]), vec![11.0, 15.0]);
        assert_eq!(super::strongest(&[3.0, 9.0, f64::NAN, 4.0]), Some(1));
        assert_eq!(super::strongest(&[]), None);
    }

    /// L145: the nearest point; a peak reaches 1.5 FWHM, 0.05° without a width.
    #[test]
    fn the_nearest_point_and_a_peaks_reach() {
        assert_eq!(super::nearest_point(&[10.0, 10.1, 10.2], 10.14), Some(1));
        assert!(super::reaches(30.0, 0.1, 30.15) && !super::reaches(30.0, 0.1, 30.16));
        assert!(super::reaches(30.0, f64::NAN, 30.07) && !super::reaches(30.0, 0.0, 30.08));
    }
}

/// The cumulative χ² and the weighted difference along the pattern (W7).
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Cumulative {
    /// Σ w(y − yc)² from the first point to each (the points left out add nothing).
    pub sum: Vec<f64>,
    /// √w (y − yc) at each point: (y − yc)/σ under variance weights; 0 for a point left out.
    pub weighted_difference: Vec<f64>,
}

impl Cumulative {
    /// The sum as a fraction of the whole (0 to 1).
    /// The cumulative χ² on the counts' scale of a Rietveld plot: the fraction times the largest of Yobs and Ycal, so
    /// that it rises from 0 to the top of the frame (Refine's "cumulative χ² (to the observed maximum)"; L146).
    pub fn on_counts(&self, yobs: &[f64], ycalc: &[f64]) -> Vec<f64> {
        let top = yobs.iter().chain(ycalc).copied().filter(|v| v.is_finite()).fold(0.0_f64, f64::max);
        self.fraction().into_iter().map(|f| f * top).collect()
    }

    pub fn fraction(&self) -> Vec<f64> {
        let total = self.sum.last().copied().unwrap_or(0.0);
        self.sum.iter().map(|s| if total > 0.0 { s / total } else { 0.0 }).collect()
    }

    /// The sum scaled as GSAS-II draws it under the pattern: its end at the observed maximum (GSASIIpwdplot.py,
    /// `CY = np.cumsum(W**2*(Y-Z)**2); scale = np.max(CY)/np.max(Y)`, GSAS-II 305a159).
    pub fn gsas2(&self, ymax: f64) -> Vec<f64> {
        self.fraction().iter().map(|f| f * ymax).collect()
    }
}

/// The cumulative χ² of a calculated pattern against the data, with the model's weights (`lsq::weights`).
pub fn cumulative(model: &Model, data: &Pattern, ycalc: &[f64]) -> Cumulative {
    let w = crate::lsq::weights(model, data, ycalc);
    let mut s = 0.0;
    let mut sum = Vec::with_capacity(w.len());
    let mut weighted_difference = Vec::with_capacity(w.len());
    for i in 0..w.len() {
        let r = w[i].sqrt() * (data.y[i] - ycalc[i]);
        s += r * r;
        sum.push(s);
        weighted_difference.push(r);
    }
    Cumulative { sum, weighted_difference }
}

/// The forbidden reflections of every phase within the pattern's 2θ range (the app's request 30; W8): (phase,
/// reflection), each reflection's 2θ with Zero added, as `drs reflections --forbidden` lists them. The d range is
/// that of the first point and the last (below 179.9°) through Bragg's law with λ1.
pub fn forbidden_in(model: &Model, x: &[f64]) -> Vec<(usize, crate::reflections::Forbidden)> {
    let (Some(&lo), Some(&hi)) = (x.first(), x.last()) else { return Vec::new() };
    let s = |t: f64| (t / 2.0).to_radians().sin() / model.lambda[0];
    let (dmin, dmax) = (0.5 / s(hi.min(179.9)), if lo > 0.0 { 0.5 / s(lo) } else { 1e9 });
    let mut out = Vec::new();
    for (ip, ph) in model.phases.iter().enumerate() {
        for mut r in crate::reflections::forbidden(&ph.group, &ph.cell, dmin, dmax, Some(model.lambda[0])) {
            r.two_theta = r.two_theta.map(|t| t + model.zero);
            out.push((ip, r));
        }
    }
    out
}
