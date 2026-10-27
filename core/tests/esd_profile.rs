//! The engine's esd's against the profile of the fit (test `esd_profile`, L054, L055).
//!
//! An esd is a statement about the χ² surface: with σ_k = √(C_kk χ²_ν), fixing parameter k at its best value ± σ_k
//! and refining all the others again raises S = Σw(y − y_c)² by χ²_ν, exactly for a model linear in its
//! parameters. This test makes that measurement on the five runs where the engine reaches FullProf's minimum
//! (L053), one parameter of a different kind in each: S_min is the engine's minimum from FullProf's converged values
//! (Eps = 10⁻³), and ΔS± that of the refinement with the parameter fixed at v ± σ; the ratio ΔS±/χ²_ν is 1 when the
//! esd is the one the surface supports. FullProf's cell esd's under TCH, half the engine's (L053), would give 4.
//!
//! The peaks are calculated over ±400 FWHM (Wdt = 400, default window). At the runs' own Wdt (10–40) S jumps
//! whenever a window edge crosses a point, and the two sides disagree: Y2O3's W gave 1.38 and 0.32, rutile's
//! anatase O z 1.00 and 0.53, Egypt's scale 1.09 and 0.52; at Wdt = 100 all but Y2O3 (a pseudo-Voigt with η ≈ 0.5,
//! whose Lorentzian tails are long) are symmetric, at 400 all are, to 0.011.
//!
//! The ratios left (L055): Y2O3's W 0.88–0.89, Yb2ZrWO8's cell a 0.95, MgAl2O4's O x 1.00, Egypt's scale 1.01,
//! rutile + anatase's O z 0.97. At ±σ/2 they are a quarter of these to 0.01 (Y2O3 0.223 and 0.217), so the surface
//! is quadratic over ±σ, and what is short of 1 is the curvature of the model itself, Σ r_i w_i ∂²y_i, which the
//! Gauss–Newton matrix JᵀWJ (the engine's, as FullProf's, manual eq. 3.7) leaves out: the profile esd is larger by
//! √(1/ratio), 6 % for Y2O3's W, 3 % for Yb2ZrWO8's cell. Frozen per run: the test fails if |ratio − 1| grows.
//!
//! The default window (L058) tapers each peak to zero beyond ±Wdt·H instead of cutting it, so that S is continuous.
//! A second pass makes the same measurement at each run's own Wdt (10–40) with it: the two sides agree to 0.010
//! (Y2O3 0.883 and 0.893, where FullProf's cut gave 1.38 and 0.32) and the ratios are those of Wdt = 400 to 0.004
//! (Y2O3 0.888, Yb2ZrWO8 0.948, MgAl2O4 0.999, Egypt 1.011, rutile + anatase 0.969). Frozen too, with the difference
//! between the sides.

use drs::lsq::{refine, statistics, Options};
use drs::powder::Model;
use drs::refinable::ParameterSet;
use std::path::Path;

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display())).iter().map(|&b| b as char).collect()
}

/// Run, code, the parameter profiled, and the frozen |ratio − 1| (largest of the two sides).
/// Then the same at the run's own Wdt with the default's tapered window edge (L058): frozen |ratio − 1| and the
/// largest difference between the two sides.
const CASES: &[(&str, &str, &str, f64, f64, f64)] = &[
    ("Y2O3", "y2o3_structure", "phase1.caglioti_w", 0.122, 0.118, 0.011),
    ("Yb2ZrWO8", "ybfluo_Riet_IRF", "phase1.cell_a", 0.054, 0.053, 0.001),
    ("MgAl2O4-s", "MgAl2O4s_Rietveld", "phase1.O.atom_position.x", 0.004, 0.004, 0.006),
    ("Egypt", "Egypt_qa_Rietveld", "phase1.scale_factor", 0.014, 0.014, 0.007),
    ("rutile-anatase", "rutana_tch", "phase2.O.atom_position.z", 0.036, 0.035, 0.008),
];

#[test]
fn esd_profile() {
    let mut fails = Vec::new();
    let opts = |eps| Options { max_cycles: 100, eps };
    for &(run, code, name, frozen_400, frozen_taper, frozen_sides) in CASES {
      for tapered in [false, true] {
        let d = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs").join(run);
        let mut p = drs::pcr::read(&latin1(&d.join(format!("{code}_conv.pcr")))).unwrap();
        p.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
        let loaded = Model::from_pcr(&p, &d, code).unwrap_or_else(|e| panic!("{run}: {e}"));
        let data = loaded.data.as_ref().unwrap();
        let mut m = loaded.model.clone();
        let out = latin1(&d.join(format!("{code}_conv.out")));
        let lim: Vec<f64> = out.lines().filter_map(|l| l.split_once("S(1/d)max:")).map(|(_, r)| r.split_whitespace().next().unwrap().parse().unwrap()).collect();
        for (ph, l) in m.phases.iter_mut().zip(&lim) {
            ph.reflection_limit = Some(*l);
        }
        // A window wide enough that S is smooth (see above); or the run's own Wdt with the default's tapered edge.
        m.convention.peak_window = false;
        if tapered {
            m.convention.window_edge = false;
        } else {
            m.wdt = 400.0;
        }
        let mode = if tapered { format!("Wdt {} tapered", m.wdt) } else { "Wdt 400".to_string() };
        let set = ParameterSet::from_pcr(&p, &m).unwrap();
        let r = refine(&mut m, data, &set, opts(1e-3), |_| {});
        let s_min = r.stats.chi2 * (r.stats.n - r.stats.p) as f64;
        let chi2 = r.stats.chi2;
        let j = set.free.iter().position(|f| f.name == name).unwrap_or_else(|| panic!("{run}: no parameter {name}"));
        let (sigma, mult) = (r.covariance[j][j].sqrt(), set.free[j].fields.clone());
        // The others, with parameter j left out (fixed).
        let mut rest = set.clone();
        rest.free.remove(j);
        let mut ratios = Vec::new();
        for sign in [1.0, -1.0] {
            let mut mk = m.clone();
            for &(f, mf) in &mult {
                mk.set(f, mk.get(f) + sign * sigma * mf);
            }
            let rk = refine(&mut mk, data, &rest, opts(1e-3), |_| {});
            // S over the same points, with N − P of the full set (the fixed parameter is still a parameter).
            let s = statistics(&mk, data, &mk.calculate(&data.x).ycalc, set.free.len());
            let ds = s.chi2 * (s.n - s.p) as f64 - s_min;
            ratios.push(ds / chi2);
            println!("{run:15} {mode:16} {name:28} {}σ = {:+.3e}: ΔS {ds:.4} ({} cycles), χ²_ν {chi2:.4}: ratio {:.3}", if sign > 0.0 { "+" } else { "−" }, sign * sigma * mult[0].1, rk.cycles, ds / chi2);
        }
        let worst = ratios.iter().map(|q| (q - 1.0).abs()).fold(0.0, f64::max);
        let sides = (ratios[0] - ratios[1]).abs();
        println!("{run:15} {mode:16} S_min {s_min:.3}; mean ratio {:.3}, largest |ratio − 1| {worst:.3}, sides differ by {sides:.3}", (ratios[0] + ratios[1]) / 2.0);
        let frozen = if tapered { frozen_taper } else { frozen_400 };
        if worst > frozen || (tapered && sides > frozen_sides) {
            fails.push(format!("{run} ({mode}): {name}: ΔS/χ²_ν {:.3} and {:.3}, |ratio − 1| {worst:.3} (frozen {frozen}), sides {sides:.3}", ratios[0], ratios[1]));
        }
      }
    }
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}
