//! Quantitative phase analysis (P33): weight fractions from the scale factors (Hill & Howard 1987, as FullProf's
//! manual writes it, p. 31): W_φ = S_φ·ATZ_φ·V_φ/t_φ / Σ_ψ S_ψ·ATZ_ψ·V_ψ/t_ψ, ATZ = Z·M_w·f²/t (p. 91), t the
//! Brindley coefficient. Only phases with ATZ > 0 take part.
//!
//! The esd's propagate the covariance of a refinement through the fractions' derivatives with respect to the
//! scales and, when the cell is refined, to the cell parameters (through V): σ²(W_φ) = gᵀ C g.
//!
//! Amorphous content with an internal standard (P33; X. Orlhac et al., J. Appl. Cryst. 34, 114–118, 2001, in
//! `manuals/fullprof/5-Quantitative_Phase_Analysis/Quant-amorphous.pdf`, the paper FullProf's QPA tutorial carries;
//! FullProf has no option for it): a known weight fraction W_s of a crystalline standard is mixed into the sample. The
//! fractions above, (W_φ)_x, assume a fully crystalline mixture; with the standard they become (W_φ)_real =
//! (W_φ)_x·W_s/(W_s)_x in the mixture (eq. 3), (W_φ)_raw = (W_φ)_real/(1 − W_s) in the sample (eq. 4), and the
//! amorphous content is what the crystalline phases leave: 1 − Σ_φ (W_φ)_real in the mixture, 1 − Σ_{φ≠s} (W_φ)_raw
//! in the sample (§4.3; their Table 3). Brindley's particle absorption contrast τ (the coefficient t above) from the
//! table the FullProf manual reproduces (p. 31; Brindley 1945): `brindley_tau`.

use crate::jacobian::Field;
use crate::powder::Model;
use crate::refinable::ParameterSet;

/// The weight fraction of each phase (None where ATZ = 0).
pub fn weight_fractions(model: &Model) -> Vec<Option<f64>> {
    let a: Vec<f64> = model.phases.iter().map(|p| if p.atz > 0.0 { p.scale * p.atz * p.cell.volume() / p.brindley } else { 0.0 }).collect();
    let total: f64 = a.iter().sum();
    model.phases.iter().zip(&a).map(|(p, &ai)| (p.atz > 0.0 && total > 0.0).then(|| ai / total)).collect()
}

/// The weight fractions with their esd's from the covariance `cov` of the free parameters of `set` (the
/// `covariance` of `lsq::Refinement`).
pub fn weight_fractions_esd(model: &Model, set: &ParameterSet, cov: &[Vec<f64>]) -> Vec<Option<(f64, f64)>> {
    with_esd(model, set, cov, weight_fractions)
}

/// The values of `f` with their esd's from the covariance: σ² = gᵀ C g, g the derivatives of each value with respect
/// to the free parameters through the fields they move that the fractions depend on (scales and cell parameters).
fn with_esd(model: &Model, set: &ParameterSet, cov: &[Vec<f64>], values: impl Fn(&Model) -> Vec<Option<f64>>) -> Vec<Option<(f64, f64)>> {
    let w0 = values(model);
    let np = set.free.len();
    // ∂W_φ/∂(free parameter j), through the fields each parameter moves: scales and cell parameters.
    let mut g = vec![vec![0.0; np]; w0.len()];
    for (j, p) in set.free.iter().enumerate() {
        for &(f, m) in &p.fields {
            if !matches!(f, Field::Scale(_) | Field::Cell(..)) {
                continue;
            }
            // A central difference: W is a smooth ratio of products, and V of the cell is a closed form.
            let v = model.get(f);
            let h = 1e-6 * v.abs().max(1e-6);
            let mut up = model.clone();
            up.set(f, v + h);
            let mut dn = model.clone();
            dn.set(f, v - h);
            let (wu, wd) = (values(&up), values(&dn));
            for (k, gk) in g.iter_mut().enumerate() {
                if let (Some(a), Some(b)) = (wu[k], wd[k]) {
                    gk[j] += m * (a - b) / (2.0 * h);
                }
            }
        }
    }
    w0.iter()
        .zip(&g)
        .map(|(w, gk)| {
            w.map(|w| {
                let var: f64 = (0..np).map(|j| (0..np).map(|k| gk[j] * cov[j][k] * gk[k]).sum::<f64>()).sum();
                (w, var.max(0.0).sqrt())
            })
        })
        .collect()
}

/// The amorphous content found with an internal standard (Orlhac et al. 2001, eqs. 3–4), each value with its esd
/// propagated from the covariance (the paper adds the phases' esd's instead, an upper bound: Table 3's amorphous esd's
/// are the sums of the phases').
#[derive(Clone, Debug, PartialEq)]
pub struct Amorphous {
    /// The phase that is the standard, and the weight fraction of it weighed into the sample.
    pub standard: usize,
    pub standard_fraction: f64,
    /// (W_φ)_real: each phase's weight fraction in the mixture with the standard (eq. 3).
    pub in_mixture: Vec<Option<(f64, f64)>>,
    /// (W_φ)_raw: each phase's weight fraction in the sample without the standard (eq. 4); None for the standard.
    pub in_sample: Vec<Option<(f64, f64)>>,
    /// The amorphous weight fraction in the mixture and in the sample.
    pub amorphous_in_mixture: (f64, f64),
    pub amorphous_in_sample: (f64, f64),
}

/// The fractions of eqs. 3–4 from the fully crystalline fractions (W_φ)_x: (W_φ)_real for every phase, (W_φ)_raw for
/// every phase but the standard, then the amorphous content in the mixture and in the sample (the last two entries).
pub fn with_standard(w_x: &[Option<f64>], standard: usize, w_s: f64) -> Vec<Option<f64>> {
    let k = w_x[standard].map(|ws_x| w_s / ws_x);
    let real: Vec<Option<f64>> = w_x.iter().map(|w| Some((*w)? * k?)).collect();
    let raw: Vec<Option<f64>> = real.iter().enumerate().map(|(i, w)| if i == standard { None } else { Some((*w)? / (1.0 - w_s)) }).collect();
    let am_mix = k.map(|_| 1.0 - real.iter().flatten().sum::<f64>());
    let am_raw = k.map(|_| 1.0 - raw.iter().flatten().sum::<f64>());
    real.into_iter().chain(raw).chain([am_mix, am_raw]).collect()
}

/// The amorphous content of a refined model whose phases include an internal standard (`Phase::internal_standard`):
/// None when no phase is one; an error when more than one is, when the standard takes no part in the analysis
/// (ATZ = 0) or when its weight fraction is not within (0, 1).
pub fn amorphous_content(model: &Model, set: &ParameterSet, cov: &[Vec<f64>]) -> Result<Option<Amorphous>, String> {
    let standards: Vec<(usize, f64)> = model.phases.iter().enumerate().filter_map(|(i, p)| p.internal_standard.map(|w| (i, w))).collect();
    let (standard, w_s) = match standards[..] {
        [] => return Ok(None),
        [one] => one,
        _ => return Err(format!("{} phases are given as the internal standard; one is", standards.len())),
    };
    let name = &model.phases[standard].name;
    if !(w_s > 0.0 && w_s < 1.0) {
        return Err(format!("phase {} ({name}): the internal standard's weight fraction {w_s} is not between 0 and 1", standard + 1));
    }
    if model.phases[standard].atz <= 0.0 {
        return Err(format!("phase {} ({name}) is the internal standard but takes no part in the quantitative analysis (ATZ = 0)", standard + 1));
    }
    let v = with_esd(model, set, cov, |m| with_standard(&weight_fractions(m), standard, w_s));
    let n = model.phases.len();
    Ok(Some(Amorphous {
        standard,
        standard_fraction: w_s,
        in_mixture: v[..n].to_vec(),
        in_sample: v[n..2 * n].to_vec(),
        amorphous_in_mixture: v[2 * n].expect("the standard has a fraction"),
        amorphous_in_sample: v[2 * n + 1].expect("the standard has a fraction"),
    }))
}

/// Brindley's table of τ against (μ_j − μ̄)R at 2θ = 0°, 90° and 180°, as the FullProf manual reproduces it (p. 31).
const BRINDLEY: [(f64, [f64; 3]); 29] = [
    (-0.50, [2.068, 2.036, 2.029]),
    (-0.40, [1.813, 1.807, 1.827]),
    (-0.30, [1.508, 1.573, 1.585]),
    (-0.20, [1.352, 1.353, 1.362]),
    (-0.10, [1.159, 1.162, 1.163]),
    (-0.09, [1.142, 1.143, 1.144]),
    (-0.08, [1.124, 1.125, 1.125]),
    (-0.07, [1.107, 1.108, 1.108]),
    (-0.06, [1.090, 1.091, 1.091]),
    (-0.05, [1.074, 1.073, 1.074]),
    (-0.04, [1.059, 1.058, 1.059]),
    (-0.03, [1.043, 1.042, 1.042]),
    (-0.02, [1.028, 1.027, 1.027]),
    (-0.01, [1.014, 1.014, 1.014]),
    (0.00, [1.000, 1.000, 1.000]),
    (0.01, [0.986, 0.986, 0.986]),
    (0.02, [0.972, 0.973, 0.973]),
    (0.03, [0.959, 0.960, 0.960]),
    (0.04, [0.945, 0.946, 0.947]),
    (0.05, [0.932, 0.933, 0.934]),
    (0.06, [0.918, 0.919, 0.921]),
    (0.07, [0.905, 0.906, 0.908]),
    (0.08, [0.892, 0.893, 0.895]),
    (0.09, [0.878, 0.879, 0.882]),
    (0.10, [0.865, 0.866, 0.870]),
    (0.20, [0.742, 0.753, 0.760]),
    (0.30, [0.640, 0.653, 0.671]),
    (0.40, [0.545, 0.569, 0.587]),
    (0.50, [0.468, 0.496, 0.529]),
];

/// Brindley's particle absorption contrast τ for (μ_j − μ̄)R (μ_j the linear absorption of the phase's particles, μ̄
/// the powder's mean, R the particle radius) at 2θ (degrees), interpolated linearly in Brindley's table (FullProf's
/// manual, p. 31) in both. None outside the table, |(μ_j − μ̄)R| > 0.5.
pub fn brindley_tau(mu_r: f64, two_theta: f64) -> Option<f64> {
    let rows = &BRINDLEY;
    if !(mu_r >= rows[0].0 - 1e-12 && mu_r <= rows[rows.len() - 1].0 + 1e-12) || !(0.0..=180.0).contains(&two_theta) {
        return None;
    }
    let j = rows.partition_point(|r| r.0 < mu_r - 1e-12).clamp(1, rows.len() - 1);
    let (a, b) = (rows[j - 1], rows[j]);
    let t = if b.0 > a.0 { ((mu_r - a.0) / (b.0 - a.0)).clamp(0.0, 1.0) } else { 0.0 };
    let at = |c: usize| a.1[c] + t * (b.1[c] - a.1[c]);
    let (c, u) = if two_theta <= 90.0 { (0, two_theta / 90.0) } else { (1, (two_theta - 90.0) / 90.0) };
    Some(at(c) + u * (at(c + 1) - at(c)))
}
