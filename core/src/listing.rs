//! The summary and the listing of a run (P57), FullProf's `.sum` and `.out` in the core's own words: what was refined
//! and how well, phase by phase, and, for the listing, cycle by cycle.
//!
//! The agreement of a block of points (FullProf's manual, eqs. 3.9–3.20, pp. 43–45), each value checked against
//! FullProf's `.sum` on its runs in `fullprof_summary`:
//! - Rp, Rwp, Rexp and χ² over N − P + C, and the conventional Rp, Rwp, Rexp with the background taken out of the
//!   denominators (`lsq::statistics`);
//! - the deviance D = 2 Σ{yᵢ ln(yᵢ/ycᵢ) − (yᵢ − ycᵢ)} (eq. 3.17; Antoniadis, Berruyer & Filhol 1990), and D per
//!   degree of freedom, D/(N − P + C) (FullProf's "Dev*" is not defined in its manual; this reading gives its print
//!   in 50 of 68 runs, D/N in Y₂O₃'s, `fullprof_summary`);
//! - the Durbin–Watson d (eq. 3.19) and its 99.9 % bound Q_D = 2[(N − 1)/(N − P) − 3.0901/√(N + 2)] (eq. 3.20);
//! - the goodness of fit in standard deviations of its expectation, (χ² − 1)/√(2/(N − P + C)) (FullProf's "N-sigma
//!   of the GoF", not defined in its manual: this reading gives its printed values).
//!
//! Two blocks, as FullProf's: all the points used, and the points with Bragg contributions (`lsq::bragg_points`);
//! FullProf's global χ² "(Bragg contrib.)" is S over all the points over N − P + C of the second.

use crate::data::Pattern;
use crate::lsq::{bragg_points, statistics_over, weights, Refinement, Statistics};
use crate::powder::Model;
use crate::refinable::ParameterSet;
use std::fmt::Write;

/// The agreement over a block of points.
#[derive(Clone, Debug, PartialEq)]
pub struct Block {
    /// The core's statistics over the block (R's, χ², the conventional R's, d).
    pub stats: Statistics,
    pub deviance: f64,
    /// D/(N − P + C).
    pub deviance_per_dof: f64,
    /// Q_D, the Durbin–Watson bound at 99.9 %.
    pub durbin_watson_bound: f64,
    /// (χ² − 1)/√(2/(N − P + C)).
    pub gof_sigmas: f64,
}

/// A phase's lines of the summary.
#[derive(Clone, Debug, PartialEq)]
pub struct PhaseLine {
    pub name: String,
    pub group: String,
    /// The reflections of the first wavelength calculated within the data.
    pub reflections: usize,
    pub r_bragg: f64,
    pub r_f: f64,
    pub volume: f64,
    /// Its su from the covariance (`ParameterSet::cell_volume`, P43); 0 without a refinement.
    pub volume_su: f64,
    /// The weight fraction and its esd (QPA, two phases or more with ATZ > 0).
    pub weight_fraction: Option<(f64, f64)>,
}

/// The summary of a run: its agreement, its phases and its parameters.
#[derive(Clone, Debug, PartialEq)]
pub struct Summary {
    pub all: Block,
    pub bragg: Block,
    /// S over all the points over N − P + C of the points with Bragg contributions.
    pub global_bragg_chi2: f64,
    pub berar_factor: f64,
    pub phases: Vec<PhaseLine>,
    pub amorphous: Option<crate::qpa::Amorphous>,
    /// The refined parameters: name, value, esd (none before a refinement).
    pub parameters: Vec<(String, f64, Option<f64>)>,
    /// Cycles run, and the stop, when there was a refinement.
    pub cycles: Option<(usize, String)>,
}

fn block(model: &Model, data: &Pattern, ycalc: &[f64], p: usize, mask: Option<&[bool]>) -> Block {
    let stats = statistics_over(model, data, ycalc, p, mask);
    let w = weights(model, data, ycalc);
    let used = |i: usize| w[i] > 0.0 && mask.is_none_or(|m| m[i]);
    let deviance: f64 = (0..ycalc.len())
        .filter(|&i| used(i))
        .map(|i| {
            let (y, yc) = (data.y[i], ycalc[i]);
            2.0 * (if y > 0.0 && yc > 0.0 { y * (y / yc).ln() } else { 0.0 } - (y - yc))
        })
        .sum();
    let (n, dof) = (stats.n as f64, (stats.n + stats.c).saturating_sub(stats.p).max(1) as f64);
    Block {
        deviance,
        deviance_per_dof: deviance / dof,
        durbin_watson_bound: 2.0 * ((n - 1.0) / (n - stats.p as f64).max(1.0) - 3.0901 / (n + 2.0).sqrt()),
        gof_sigmas: (stats.chi2 - 1.0) / (2.0 / dof).sqrt(),
        stats,
    }
}

/// The summary at the model's values: after a refinement (its free parameters, esd's, cycles and stop) or of a
/// calculation with `set`'s parameters counted as refined (P), as FullProf counts them in a run of no cycles.
pub fn summary(model: &Model, data: &Pattern, set: &ParameterSet, r: Option<&Refinement>) -> Summary {
    let c = model.calculate(&data.x);
    let p = set.free.len();
    let all = block(model, data, &c.ycalc, p, None);
    let mask = bragg_points(data, &c.peaks);
    let bragg = block(model, data, &c.ycalc, p, Some(&mask));
    let s_all = all.stats.chi2 * ((all.stats.n + all.stats.c).saturating_sub(p).max(1)) as f64;
    let global_bragg_chi2 = s_all / ((bragg.stats.n + bragg.stats.c).saturating_sub(p).max(1)) as f64;
    let cov = r.map(|r| r.covariance.clone()).unwrap_or_default();
    let fractions = if cov.is_empty() { crate::qpa::weight_fractions(model).into_iter().map(|w| w.map(|w| (w, 0.0))).collect() } else { crate::qpa::weight_fractions_esd(model, set, &cov) };
    let qpa = fractions.iter().flatten().count() > 1;
    let (lo, hi) = (data.x.first().copied().unwrap_or(0.0), data.x.last().copied().unwrap_or(0.0));
    let phases = model
        .phases
        .iter()
        .enumerate()
        .map(|(k, ph)| PhaseLine {
            name: ph.name.clone(),
            group: ph.group.setting.as_ref().map_or_else(|| ph.group.hall.clone(), |st| st.hm.to_string()),
            reflections: c.peaks.iter().filter(|q| q.phase == k && q.wavelength == 1 && q.two_theta >= lo && q.two_theta <= hi).count(),
            r_bragg: all.stats.phases.get(k).map_or(0.0, |x| x.0),
            r_f: all.stats.phases.get(k).map_or(0.0, |x| x.1),
            volume: ph.cell.volume(),
            volume_su: set.cell_volume(model, k, &cov).1,
            weight_fraction: if qpa { fractions[k] } else { None },
        })
        .collect();
    let amorphous = if cov.is_empty() { None } else { crate::qpa::amorphous_content(model, set, &cov).ok().flatten() };
    let parameters = match r {
        Some(r) => r.params.iter().map(|q| (q.name.clone(), q.value, Some(q.esd))).collect(),
        None => set.free.iter().map(|f| (f.name.clone(), model.get(f.fields[0].0) / f.fields[0].1, None)).collect(),
    };
    Summary {
        global_bragg_chi2,
        berar_factor: crate::lsq::berar_factor(model, data, &c.ycalc, &c.peaks),
        all,
        bragg,
        phases,
        amorphous,
        parameters,
        cycles: r.map(|r| (r.cycles, r.stop.clone())),
    }
}

fn block_text(s: &mut String, title: &str, b: &Block) {
    let t = &b.stats;
    let _ = writeln!(s, "{title}");
    let _ = writeln!(s, "  N − P + C        {} ({} points, {} parameters, {} constraints)", t.n + t.c - t.p.min(t.n + t.c), t.n, t.p, t.c);
    let _ = writeln!(s, "  Rp {:.2}  Rwp {:.2}  Rexp {:.2}  χ² {:.4}", t.rp, t.rwp, t.rexp, t.chi2);
    let _ = writeln!(s, "  conventional (background taken out): Rp {:.2}  Rwp {:.2}  Rexp {:.2}", t.rp_background_free, t.rwp_background_free, t.rexp_background_free);
    let _ = writeln!(s, "  deviance D {:.4e}  D/(N − P + C) {:.4}", b.deviance, b.deviance_per_dof);
    let _ = writeln!(s, "  Durbin–Watson d {:.4}  (99.9 % bound Q_D {:.4}: {})", t.durbin_watson, b.durbin_watson_bound, if t.durbin_watson < b.durbin_watson_bound {
        "positive serial correlation"
    } else if t.durbin_watson > 4.0 - b.durbin_watson_bound {
        "negative serial correlation"
    } else {
        "no serial correlation"
    });
    let _ = writeln!(s, "  goodness of fit: χ² {:.3} standard deviations of its expectation above 1", b.gof_sigmas);
}

/// The summary as text: `title` (the control file and the run), the agreement in its two blocks, each phase, the
/// amorphous content, and the parameters with their esd's.
pub fn summary_text(title: &str, sm: &Summary) -> String {
    let mut s = String::new();
    let _ = writeln!(s, "drs — summary of a run (FullProf's .sum in the core's words)\n{title}\n");
    if let Some((n, stop)) = &sm.cycles {
        let _ = writeln!(s, "Refinement: {n} cycles; {stop}\n");
    } else {
        let _ = writeln!(s, "Calculation at the file's values: no cycle run\n");
    }
    block_text(&mut s, "Agreement, all the points used", &sm.all);
    let _ = writeln!(s);
    block_text(&mut s, "Agreement, the points with Bragg contributions", &sm.bragg);
    let _ = writeln!(s, "\nGlobal χ² (S over all the points, N − P + C of the points with Bragg contributions) {:.4}", sm.global_bragg_chi2);
    let _ = writeln!(s, "Bérar factor {:.4}: the esd's times it estimate the probable errors (Bérar & Lelann 1991)\n", sm.berar_factor);
    for (k, ph) in sm.phases.iter().enumerate() {
        let _ = writeln!(s, "Phase {} {} ({})", k + 1, ph.name, ph.group);
        let _ = writeln!(s, "  {} reflections in the data; R_Bragg {:.2}  R_F {:.2}; volume {:.4}{} Å³", ph.reflections, ph.r_bragg, ph.r_f, ph.volume, if ph.volume_su > 0.0 { format!(" ± {:.4}", ph.volume_su) } else { String::new() });
        if let Some((w, e)) = ph.weight_fraction {
            let _ = writeln!(s, "  weight fraction {} %", crate::cif_write::with_su(100.0 * w, 100.0 * e, 2));
        }
    }
    if let Some(a) = &sm.amorphous {
        let pc = |v: (f64, f64)| crate::cif_write::with_su(100.0 * v.0, 100.0 * v.1, 2);
        let _ = writeln!(s, "Amorphous content {} % of the sample, {} % of the mixture (phase {} the internal standard at {} %)", pc(a.amorphous_in_sample), pc(a.amorphous_in_mixture), a.standard + 1, 100.0 * a.standard_fraction);
    }
    let _ = writeln!(s, "\n{:<45} {:>16} {:>14}", "parameter", "value", "esd");
    for (name, v, e) in &sm.parameters {
        let _ = writeln!(s, "{name:<45} {v:16.6} {:>14}", e.map_or("-".into(), |e| format!("{e:.6}")));
    }
    s
}

/// The listing of a refinement, cycle by cycle (FullProf's `.out` in the core's words): each cycle's agreement, the
/// steps refused, the steering's notes, then every parameter's value, shift, esd and shift/esd; the summary at the end.
pub fn listing_text(title: &str, cycles: &[crate::lsq::Cycle], sm: &Summary) -> String {
    let mut s = String::new();
    let _ = writeln!(s, "drs — listing of a run (FullProf's .out in the core's words)\n{title}\n");
    for c in cycles {
        let t = &c.stats;
        let _ = writeln!(s, "Cycle {}: χ² {:.4} (from {:.4})  Rp {:.2}  Rwp {:.2}  Rexp {:.2}  d {:.4}  Marquardt λ {:.3e}", c.cycle, t.chi2, c.chi2_start, t.rp, t.rwp, t.rexp, t.durbin_watson, c.marquardt);
        for tr in &c.tries {
            let _ = writeln!(s, "  step refused: λ {:.3e} gave χ² {:.4}", tr.lambda, tr.chi2);
        }
        for n in &c.notes {
            let _ = writeln!(s, "  {n}");
        }
        let _ = writeln!(s, "  {:<45} {:>16} {:>14} {:>12} {:>9}", "parameter", "value", "shift", "esd", "shift/esd");
        for q in &c.params {
            let _ = writeln!(s, "  {:<45} {:16.6} {:14.6} {:12.6} {:9.3}{}", q.name, q.value, q.shift, q.esd, q.shift_esd, if q.at_limit { "  at its limit" } else { "" });
        }
        for (k, (rb, rf)) in t.phases.iter().enumerate() {
            let _ = writeln!(s, "  phase {}: R_Bragg {rb:.2}  R_F {rf:.2}", k + 1);
        }
        let _ = writeln!(s);
    }
    s.push_str(&summary_text(title, sm));
    s
}
