//! The default model against FullProf's rules, and the import of a `.pcr` into it (tests `default_rules`,
//! `pcr_import`; L046, L047, L050).
//!
//! `default_rules` checks the default's formulas against the papers in `manuals/papers/` and the exact identities the
//! import relies on. `pcr_import` imports every FullProf run of `manuals/fullprof/runs/` (the ten runs, the
//! preferred-orientation run and the 35 window experiments) into the default model and compares its χ² with that of
//! the `.pcr`'s values under FullProf's rules, over the same points and weights:
//! - translated as imported (FullProf's peak window and reflections past the data kept, having no exact
//!   translation), the χ² must agree to 10⁻⁹ of its value: every other rule is translated exactly;
//! - with those two rules also switched to the default, the χ² moves by what they explain; that change is printed per
//!   run and frozen (the test fails if it grows).

use drs::powder::{agreement_weighted, Convention, Model};
use drs::profile::{berar_baldinozzi, fullprof_axial};
use std::path::Path;

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display())).iter().map(|&b| b as char).collect()
}

#[test]
fn default_rules() {
    let p = [0.1, 0.02, 0.03, 0.004];
    let th = 10f64.to_radians();
    for &z in &[-2.0, -0.7, 0.0, 0.3, 1.5] {
        // Bérar & Baldinozzi (1993): u = dθ/w positive on the low-angle side, F_a = 2u e^{−u²} (eq. 12),
        // F_b = (8u³ − 12u) e^{−u²} (eq. 13), A(θ) = A0/tan θ + A1/tan 2θ (eq. 15).
        let u = -z;
        let fa = 2.0 * u * f64::exp(-u * u);
        let fb = (8.0 * u * u * u - 12.0 * u) * f64::exp(-u * u);
        let want = 1.0 + (p[0] * fa + p[1] * fb) / th.tan() + (p[2] * fa + p[3] * fb) / (2.0 * th).tan();
        assert!((berar_baldinozzi(z, 20.0, p, false) - want).abs() < 1e-14);
        // FullProf's F_b = (4u² − 3) F_a is the paper's + 3 F_a: FullProf's (P1, P2) are the paper's (P1 + 3 P2, P2).
        let t = [p[0] + 3.0 * p[1], p[1], p[2] + 3.0 * p[3], p[3]];
        assert!((berar_baldinozzi(z, 20.0, p, true) - berar_baldinozzi(z, 20.0, t, false)).abs() < 1e-14);
    }
    // CrysFML's reading of a small D_L, and the pair as given otherwise.
    assert_eq!(fullprof_axial(0.06, 0.0), (0.03, 0.03));
    assert_eq!(fullprof_axial(0.06, 0.02), (0.06, 0.02));
}

/// The FullProf runs: name, folder, the .pcr's name (CODE, or `input` for the window experiments).
fn runs() -> Vec<(String, String, String)> {
    let base = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs");
    let mut v: Vec<(String, String, String)> = [
        ("Y2O3", "y2o3_structure"),
        ("Co-YBCO", "Cocus-asym"),
        ("MgAl2O4-s", "MgAl2O4s_Rietveld"),
        ("rutile-anatase", "rutana_tch"),
        ("Egypt", "Egypt_qa_Rietveld"),
        ("YVO3", "YVO3a"),
        ("Yb2ZrWO8", "ybfluo_Riet_IRF"),
        ("LaB6", "LaB6rx"),
        ("NdNi-HTO", "NdNi_HTO"),
        ("NdNi-LTOa", "NdNi_LTOa"),
        ("NdNi-LTOb", "NdNi_LTOb"),
        ("Pd3MnD", "pdfinal"),
        ("Yb2ZrWO8-PO-March", "ybfluo_Riet_IRF"),
        // P25's rest (L267).
        ("prba10k", "prba10k"),
        ("cuf1k", "cuf1k"),
        ("lac", "lac"),
    ]
    .iter()
    .map(|(r, c)| (r.to_string(), base.join(r).to_string_lossy().into_owned(), c.to_string()))
    .collect();
    let mut exps: Vec<_> = std::fs::read_dir(base.join("window_experiments")).unwrap().flatten().map(|e| e.path()).filter(|p| p.join("input.pcr").is_file()).collect();
    exps.sort();
    for e in exps {
        v.push((format!("window_experiments/{}", e.file_name().unwrap().to_string_lossy()), e.to_string_lossy().into_owned(), "input".into()));
    }
    v
}

/// Frozen |Δχ²|/χ² of the fully default model (peak window, its tapered edge and reflections too) against FullProf's
/// rules, per run (L050: the value measured, + 1%; re-measured in 0.14.0, L058, when the default window gained its
/// taper beyond ±Wdt·H: the tails it adds are what moves most runs). The window experiments are modified copies, not fits: there the default
/// window (±Wdt·H with the Lorentzian strain term) often fits better than FullProf's.
const UNTRANSLATED: &[(&str, f64)] = &[
    ("Y2O3", 3.030e-04),
    ("Co-YBCO", 1.434e-03),
    ("MgAl2O4-s", 1.667e-03),
    ("rutile-anatase", 2.586e-02),
    ("Egypt", 1.212e-03),
    ("YVO3", 6.242e-02),
    ("Yb2ZrWO8", 3.404e-03),
    ("LaB6", 4.060e-02),
    // The default's Sycos·cos θ against FullProf's cos 2θ with Sycos = 0.145 (L267): measured 9.105, + 1%. Before
    // L267 both took cos θ (3.949e-3).
    ("NdNi-HTO", 9.196e+00),
    // Imported since 0.37.0 (P25, its nuclear phase with Stephens' quartic form; L211): measured 6.990e-3, + 1%.
    ("NdNi-LTOa", 7.060e-03),
    ("Yb2ZrWO8-PO-March", 3.464e-03),
    // Imported since 0.47.0 (P25's rest, L267): measured, + 1%. cuf1k's includes its Sycos (−0.003) under cos θ.
    ("Pd3MnD", 1.907e-03),
    ("prba10k", 3.889e-02),
    ("cuf1k", 3.009e-03),
    ("lac", 2.292e-03),
    ("window_experiments/irf1_narrow_x0.645_wdt20", 1.677e-01),
    ("window_experiments/irf1_narrow_x0_noXins_ref", 3.656e-04),
    ("window_experiments/irf1_narrow_x0_ref", 5.050e-04),
    ("window_experiments/w4_x0", 1.091e-01),
    ("window_experiments/w4_x0.1", 1.162e-01),
    ("window_experiments/w4_x0.3", 5.757e-02),
    ("window_experiments/w4_x0.645", 1.162e-02),
    ("window_experiments/w4_x0.645_y0", 2.363e-01),
    ("window_experiments/w4_x0.645_y0_g0", 6.474e-01),
    ("window_experiments/w4_x0_y0", 5.727e-04),
    ("window_experiments/w4_x1.5", 5.262e-02),
    ("window_experiments/w8_x0", 3.899e-02),
    ("window_experiments/w8_x0.645", 1.909e-02),
    ("window_experiments/wdt10_x0.645_narrow", 3.252e-01),
    ("window_experiments/wdt10_x0_narrow", 1.818e-04),
    ("window_experiments/wdt20_x0.645_narrow", 1.040e-01),
    ("window_experiments/wdt20_x0_narrow", 5.141e-05),
    ("window_experiments/wdt40_x0.645_narrow", 1.454e-02),
    ("window_experiments/wdt40_x0_narrow", 7.868e-06),
    ("window_experiments/wdt4_x0.645_narrow", 6.474e-01),
    ("window_experiments/wdt4_x0_narrow", 7.565e-04),
    ("window_experiments/yb_FCJ_x0", 7.232e-02),
    ("window_experiments/yb_FCJ_x0.645", 1.081e-02),
    ("window_experiments/yb_K2_x0", 4.909e-03),
    ("window_experiments/yb_K2_x0.645", 4.959e-02),
    ("window_experiments/yb_all_x0", 1.656e-03),
    ("window_experiments/yb_all_x0.645", 3.485e-02),
    ("window_experiments/yb_irf1_x0", 3.929e-03),
    ("window_experiments/yb_irf1_x0.645", 4.909e-02),
    ("window_experiments/yb_plain_x0", 1.091e-01),
    ("window_experiments/yb_plain_x0.645", 1.162e-02),
    ("window_experiments/yvo3_X0", 1.545e-02),
    ("window_experiments/yvo3_base_w4", 9.272e-01),
    ("window_experiments/yvo3_noK2", 6.393e-01),
    ("window_experiments/yvo3_noK2_noasym", 6.474e-01),
    ("window_experiments/yvo3_noXins", 9.423e-01),
    ("window_experiments/yvo3_noasym", 8.868e-01),
];

#[test]
fn pcr_import() {
    let mut fails = Vec::new();
    let (mut n, mut refused) = (0, Vec::new());
    for (name, dir, code) in runs() {
        let d = Path::new(&dir);
        let mut p = drs::pcr::read(&latin1(&d.join(format!("{code}.pcr")))).unwrap();
        p.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
        let stem = if code == "input" { "ybfluo_Riet_IRF" } else { code.as_str() };
        let fp = match Model::from_pcr_fullprof(&p, d, stem) {
            Ok(l) => l,
            Err(e) => {
                refused.push(format!("{name}: {e}"));
                continue;
            }
        };
        let imp = Model::from_pcr(&p, d, stem).unwrap();
        let data = fp.data.as_ref().expect("data");
        let x = &data.x;
        let (cf, ci) = (fp.model.calculate(x), imp.model.calculate(x));
        // The same points and weights for both: FullProf's excluded regions.
        let chi2 = |y: &[f64]| agreement_weighted(data, y, &fp.model.excluded, 0, fp.model.weighting).3;
        let (xf, xi) = (chi2(&cf.ycalc), chi2(&ci.ycalc));
        let kept = imp.model.convention.fullprof_rules();
        let mut all = imp.model.clone();
        all.convention = Convention::DEFAULT;
        let xa = chi2(&all.calculate(x).ycalc);
        let (di, da) = ((xi - xf).abs() / xf, (xa - xf).abs() / xf);
        println!("{name:45} χ² {xf:12.6} imported {xi:12.6} (Δ {di:.1e}; kept: {}); all default {xa:12.6} (Δ {da:.2e})", kept.join(", "));
        if di > 1e-9 {
            fails.push(format!("{name}: imported χ² {xi} against FullProf's rules {xf}"));
        }
        let frozen = UNTRANSLATED.iter().find(|u| u.0 == name).map_or(0.0, |u| u.1);
        if da > frozen {
            fails.push(format!("{name}: the untranslated rules move χ² by {da:.3e}, more than the {frozen:.3e} frozen"));
        }
        n += 1;
    }
    for r in &refused {
        println!("not imported: {r}");
    }
    println!("{n} runs imported, {} refused", refused.len());
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}
