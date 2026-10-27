//! Le Bail's extraction of integrated intensities (P35; Le Bail, Duroy & Fourquet 1988, as FullProf's profile
//! matching with constant scale, Jbt = 2, manual p. 14).
//!
//! Each reflection's intensity is replaced by its "observed" one, the Rietveld formula (manual eq. 3.15, p. 44):
//! I_obs,h = I_calc,h Σ_i Ω_hi (y_i − b_i)/(y_c,i − b_i), with Ω_hi the reflection's normalised profile at point i,
//! y_c the whole calculated pattern and b the background, over the points used. The profile is normalised over the
//! points: I_obs,h = I_calc,h Σ_i v_hi r_i / Σ_i v_hi, v_hi the reflection's contribution at point i (both
//! wavelengths summed), r_i = (y_i − b_i)/(y_c,i − b_i). A reflection whose calculated intensity is 0 keeps 0, as
//! FullProf's (p. 44); a negative result is set to 0. The engine calls it at
//! the start of every cycle (`lsq::refine`), so the intensities follow the profile and cell as they are refined.

use crate::data::Pattern;
use crate::powder::{Intensities, Model};
use std::collections::BTreeMap;

impl Model {
    /// Whether any phase has its intensities extracted.
    pub fn has_extracted(&self) -> bool {
        self.phases.iter().any(|p| p.extract && matches!(p.intensities, Intensities::List(_)))
    }

    /// One Le Bail iteration over the extracted phases. Returns the largest change of an intensity over the largest
    /// intensity of its phase.
    pub fn le_bail(&mut self, data: &Pattern) -> f64 {
        if !self.has_extracted() {
            return 0.0;
        }
        let c = self.calculate(&data.x);
        let excluded = self.excluded_points(&data.x);
        // Per reflection (both wavelengths): Σ v_i r_i and Σ v_i over the points used.
        let mut obs: BTreeMap<(usize, [i32; 3]), (f64, f64)> = BTreeMap::new();
        for pk in c.peaks.iter().filter(|pk| self.phases[pk.phase].extract) {
            let e = obs.entry((pk.phase, pk.hkl)).or_insert((0.0, 0.0));
            self.peak_profile(&data.x, &excluded, pk, |i, v| {
                let net = c.ycalc[i] - c.background[i];
                if net > 0.0 && data.sigma[i] > 0.0 {
                    e.0 += v * (data.y[i] - c.background[i]) / net;
                    e.1 += v;
                }
            });
        }
        let mut change = 0.0f64;
        for (k, ph) in self.phases.iter_mut().enumerate() {
            if !ph.extract {
                continue;
            }
            if let Intensities::List(list) = &mut ph.intensities {
                let top = list.iter().map(|r| r.2.abs()).fold(0.0, f64::max).max(1e-300);
                for r in list.iter_mut() {
                    // Reflections outside the data keep their value.
                    if let Some(&(num, den)) = obs.get(&(k, r.0)).filter(|e| e.1 > 0.0) {
                        let new = (r.2 * num / den).max(0.0);
                        change = change.max((new - r.2).abs() / top);
                        r.2 = new;
                    }
                }
            }
        }
        change
    }
}

/// One reflection of a file of extracted intensities (P35): h k l, its 2θ (Zero included) and FWHM (degrees), the
/// intensity extracted, G = j·Lp·|F|² as a Le Bail phase holds it (FullProf's "Obs = j LP F^2", the intensity its
/// `.hkl` lists), and its σ.
#[derive(Clone, Debug, PartialEq)]
pub struct Extracted {
    pub hkl: [i32; 3],
    pub two_theta: f64,
    pub fwhm: f64,
    pub g: f64,
    pub sigma: f64,
}

/// FullProf's Rmub and Rmuc when the `.pcr` gives zeros (manual, line 19 of a phase, note on Jvi = 11).
pub const CLUSTER_DEFAULTS: (f64, f64) = (1.0, 0.2);

/// The clusters of overlapped reflections, FullProf's rule for its `CODFILn_cltr.int` (manual, note on Jvi = 11 and
/// Jbt = 2): reflections i and i + 1, in the order of 2θ, belong to the same cluster if T(i+1) − T(i) < ½(H_i +
/// H_{i+1})·Rmub, or if T(i+1) − T(i) < ½(H_i + H_{i+1}) and G(i+1) < G_sum·Rmuc, T the Bragg position, H the FWHM, G the
/// intensity and G_sum the intensity of the cluster so far. Each range is a cluster.
pub fn clusters(rows: &[Extracted], rmub: f64, rmuc: f64) -> Vec<std::ops::Range<usize>> {
    let mut out = Vec::new();
    let (mut start, mut gsum) = (0, 0.0);
    for i in 0..rows.len() {
        gsum += rows[i].g;
        let joined = rows.get(i + 1).is_some_and(|n| {
            let (dt, half) = (n.two_theta - rows[i].two_theta, 0.5 * (rows[i].fwhm + n.fwhm));
            dt < half * rmub || (dt < half && n.g < gsum * rmuc)
        });
        if !joined {
            out.push(start..i + 1);
            (start, gsum) = (i + 1, 0.0);
        }
    }
    out
}

/// The extracted intensities of a Le Bail phase at the model's values: for each reflection of the first wavelength,
/// its intensity partitioned from the data (`diagnostics::reflections`, I obs and its σ) in the phase's units, times
/// the list's intensity over the reflection's calculated one; in the order of 2θ.
pub fn extracted(model: &Model, data: &Pattern, phase: usize) -> Result<Vec<Extracted>, String> {
    let ph = model.phases.get(phase).ok_or(format!("no phase {}", phase + 1))?;
    let Intensities::List(list) = &ph.intensities else {
        return Err(format!("phase {} ({}) has a structure: extracted intensities belong to a Le Bail phase", phase + 1, ph.name));
    };
    let of: BTreeMap<[i32; 3], f64> = list.iter().map(|r| (r.0, r.2)).collect();
    let mut rows: Vec<Extracted> = crate::diagnostics::reflections(model, data)
        .into_iter()
        .filter(|r| r.phase == phase && r.wavelength == 1)
        .map(|r| {
            let ratio = if r.icalc > 0.0 { of.get(&r.hkl).copied().unwrap_or(0.0) / r.icalc } else { 0.0 };
            Extracted { hkl: r.hkl, two_theta: r.two_theta, fwhm: r.fwhm, g: r.iobs * ratio, sigma: r.sigma * ratio }
        })
        .collect();
    rows.sort_by(|a, b| a.two_theta.total_cmp(&b.two_theta));
    Ok(rows)
}

/// FullProf's `CODFILn_cltr.int` (manual: "Files suitable as input for integrated intensity refinements", Jbt = 2,
/// More = 1, Jvi = 11), read by FullProf with Irf = 4 and Ipow = 2: a title, the format `(3i4,2f16.5,i4,3f14.4)`, λ
/// with Itypdata = 0 and Ipow = 2, then per reflection h k l, G, σ, the code 1, c1 = c2 = 0 and 2θ. Within a cluster
/// every reflection but the last has G = −1 and its own σ; the last carries the cluster's sum and the root of the sum of
/// its σ²'s. A σ of 0 (no calculated intensity to partition) is written as FullProf writes its near-zero observations,
/// G = 0.0001 with σ = 1.0 (its files in `manuals/`; the threshold at which FullProf does so is not in its manual, so
/// only a σ of 0 is replaced), so that a weighted refinement can read the line.
pub fn int_file(title: &str, lambda: f64, rows: &[Extracted], (rmub, rmuc): (f64, f64)) -> String {
    use std::fmt::Write;
    let mut s = format!("! {title:<46} Overlapped reflections re-grouped  -> Obs = j LP F^2\n(3i4,2f16.5,i4,3f14.4)\n{lambda:9.5}{:4}{:4}{:10.4}\n", 0, 2, 0.0);
    for c in clusters(rows, rmub, rmuc) {
        let (g, var): (f64, f64) = rows[c.clone()].iter().fold((0.0, 0.0), |(g, v), r| (g + r.g, v + r.sigma * r.sigma));
        for i in c.clone() {
            let r = &rows[i];
            let (g, sigma) = if i + 1 == c.end { (g, var.sqrt()) } else { (-1.0, r.sigma) };
            let (g, sigma) = if sigma > 0.0 { (g, sigma) } else { (if g == -1.0 { g } else { g.max(1e-4) }, 1.0) };
            let _ = writeln!(s, "{:4}{:4}{:4}{g:16.5}{sigma:16.5}{:4}{:14.4}{:14.4}{:14.4}", r.hkl[0], r.hkl[1], r.hkl[2], 1, 0.0, 0.0, r.two_theta);
        }
    }
    s
}
