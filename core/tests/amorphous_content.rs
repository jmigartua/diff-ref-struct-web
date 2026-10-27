//! Amorphous content with an internal standard (P33; X. Orlhac, C. Fillet, P. Deniard, A. M. Dulac & R. Brec, J. Appl.
//! Cryst. 34, 114–118, 2001: `manuals/fullprof/5-Quantitative_Phase_Analysis/Quant-amorphous.pdf`, the paper FullProf's
//! QPA tutorial carries; FullProf's manual has no option for it) and Brindley's τ (FullProf's manual, p. 31).
//!
//! - The paper's Table 3 (four glasses, 8 wt% TiO₂ as the standard): from the crystalline phases' (W_φ)_raw, the
//!   fully crystalline fractions (W_φ)_x a refinement would give are built back with eqs. 3–4, and `qpa::with_standard`
//!   returns each (W_φ)_raw and the amorphous content of the table to its 0.01 wt%. The table's amorphous esd's are
//!   the sums of the phases' (0.10 + 0.10 = 0.20, …): the core propagates the covariance instead.
//! - FullProf's rutile + anatase run (`runs/rutile-anatase/rutana_tch`), refined by the engine, with anatase as the
//!   standard: given at its refined fraction the amorphous content is 0; given at 0.9 of it, the closed form
//!   1 − k(1 − W)/W with k = W_s/(1 − W_s); the esd's against the closed form's derivative times W's esd.
//! - The control file's `internal_standard_weight_fraction` written and read back; `drs refine` prints the content;
//!   the errors (two standards, ATZ = 0, a fraction outside (0, 1)).
//! - Brindley's τ: the table's nodes, its interpolation, and the paper's eq. 8, τ = 1 − 1.450x + 1.426x² for
//!   x = (μ_j − μ̄)R, within 0.0043 of the table for |x| ≤ 0.1, where the paper uses it (largest at x = 0.1, 0°).

use drs::qpa::{amorphous_content, brindley_tau, with_standard};
use std::path::Path;

mod common;
use common::latin1;

/// A sample of the paper's Table 3: its name, [(phase, (W_φ)_raw AMC in wt%, its esd)], the amorphous wt% and its esd.
type Sample = (&'static str, &'static [(&'static str, f64, f64)], f64, f64);

const TABLE_3: [Sample; 4] = [
    ("quenched glass", &[("RuO2", 1.69, 0.10), ("Pd-Te", 0.27, 0.10)], 98.04, 0.20),
    ("20 h", &[("RuO2", 1.97, 0.11), ("Pd-Te", 0.39, 0.10), ("CaMoO4", 1.50, 0.25), ("CeO2", 0.06, 0.07), ("ZnCr2O4", 0.10, 0.13)], 95.98, 0.66),
    ("50 h", &[("RuO2", 1.09, 0.10), ("Pd-Te", 0.22, 0.10), ("CaMoO4", 1.63, 0.10), ("CeO2", 0.11, 0.03), ("ZnCr2O4", 0.19, 0.05)], 96.76, 0.38),
    ("1000 h", &[("RuO2", 1.48, 0.10), ("Pd-Te", 0.23, 0.06), ("CaMoO4", 1.92, 0.26), ("CeO2", 0.60, 0.08), ("ZnCr2O4", 0.78, 0.14)], 94.99, 0.64),
];

/// The weight of TiO₂ added (§3.1: "8 wt% TiO2 was added").
const W_STANDARD: f64 = 0.08;

#[test]
fn orlhac_table_3() {
    for (sample, phases, amorphous, esd) in TABLE_3 {
        // Eqs. 3–4 backwards: (W_φ)_real = (W_φ)_raw (1 − W_s), the standard W_s; the refinement's (W_φ)_x are these
        // over their sum (a fully crystalline mixture).
        let mut real: Vec<f64> = phases.iter().map(|p| p.1 / 100.0 * (1.0 - W_STANDARD)).collect();
        real.push(W_STANDARD);
        let sum: f64 = real.iter().sum();
        let w_x: Vec<Option<f64>> = real.iter().map(|w| Some(w / sum)).collect();
        let s = phases.len();
        let v = with_standard(&w_x, s, W_STANDARD);
        let n = w_x.len();
        for (k, p) in phases.iter().enumerate() {
            let raw = 100.0 * v[n + k].unwrap();
            assert!((raw - p.1).abs() < 1e-9, "{sample} {}: {raw} against {}", p.0, p.1);
        }
        assert_eq!(v[n + s], None, "the standard has no fraction in the sample");
        let am = 100.0 * v[2 * n + 1].unwrap();
        println!("{sample}: amorphous {am:.3} wt% of the sample (Table 3: {amorphous}), {:.3} of the mixture", 100.0 * v[2 * n].unwrap());
        assert!((am - amorphous).abs() <= 0.005 + 1e-9, "{sample}: amorphous {am} against {amorphous}");
        // The table's esd of the amorphous content is the sum of the phases' (an upper bound, not the propagation).
        let sum_esd: f64 = phases.iter().map(|p| p.2).sum();
        assert!((sum_esd - esd).abs() < 1e-9, "{sample}: Σ esd {sum_esd} against {esd}");
    }
}

fn rutana() -> (drs::powder::Model, drs::refinable::ParameterSet, drs::lsq::Refinement) {
    let d = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs/rutile-anatase");
    let p = drs::pcr::read(&latin1(&d.join("rutana_tch.pcr"))).unwrap();
    let l = drs::powder::Model::from_pcr(&p, &d, "rutana_tch").unwrap();
    let data = l.data.unwrap();
    let mut m = l.model;
    let set = drs::refinable::ParameterSet::from_pcr(&p, &m).unwrap();
    let r = drs::lsq::refine(&mut m, &data, &set, drs::lsq::Options { max_cycles: 3, eps: 0.1 }, |_| {});
    (m, set, r)
}

#[test]
fn rutile_anatase_with_a_standard() {
    let (mut m, set, r) = rutana();
    assert!(amorphous_content(&m, &set, &r.covariance).unwrap().is_none(), "no standard, no amorphous content");
    let w = drs::qpa::weight_fractions_esd(&m, &set, &r.covariance);
    let (wa, ea) = w[1].unwrap();
    println!("anatase (W)_x {wa:.5} ± {ea:.5}");
    for f in [1.0, 0.9] {
        let w_s = f * wa;
        m.phases[1].internal_standard = Some(w_s);
        let a = amorphous_content(&m, &set, &r.covariance).unwrap().unwrap();
        let k = w_s / (1.0 - w_s);
        let (am, e) = a.amorphous_in_sample;
        let closed = 1.0 - k * (1.0 - wa) / wa;
        // ∂/∂W of 1 − k(1 − W)/W = k/W²: W's esd carried through.
        let e_closed = k / (wa * wa) * ea;
        println!("standard at {f} × (W)_x: amorphous {am:.6} ± {e:.6} of the sample (closed form {closed:.6} ± {e_closed:.6}), {:.6} of the mixture", a.amorphous_in_mixture.0);
        assert!((am - closed).abs() < 1e-12, "amorphous {am} against {closed}");
        assert!((e / e_closed - 1.0).abs() < 1e-4, "esd {e} against {e_closed}");
        // The mixture: 1 − Σ (W_φ)_real, and (W_φ)_raw = (W_φ)_real/(1 − W_s) for rutile.
        let (rr, _) = a.in_mixture[0].unwrap();
        assert!((a.amorphous_in_mixture.0 - (1.0 - rr - w_s)).abs() < 1e-12);
        assert!((a.in_sample[0].unwrap().0 - rr / (1.0 - w_s)).abs() < 1e-12);
        assert_eq!(a.in_sample[1], None);
        if f == 1.0 {
            assert!(am.abs() < 1e-12, "the standard at its refined fraction leaves no amorphous content: {am}");
        }
    }
    // The errors.
    m.phases[0].internal_standard = Some(0.5);
    assert!(amorphous_content(&m, &set, &r.covariance).unwrap_err().contains("2 phases"));
    m.phases[0].internal_standard = None;
    m.phases[1].internal_standard = Some(1.2);
    assert!(amorphous_content(&m, &set, &r.covariance).unwrap_err().contains("not between 0 and 1"));
    m.phases[1].internal_standard = Some(0.1);
    m.phases[1].atz = 0.0;
    assert!(amorphous_content(&m, &set, &r.covariance).unwrap_err().contains("ATZ = 0"));
}

/// The control file's key, written and read back, and `drs refine` printing the content.
#[test]
fn internal_standard_in_the_control_file() {
    let src = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs/rutile-anatase");
    let d = Path::new(env!("CARGO_MANIFEST_DIR")).join("target/test-amorphous");
    std::fs::create_dir_all(&d).unwrap();
    for f in ["rutana_tch.pcr", "rutana.dat", "rx_cukalpha.irf"] {
        std::fs::copy(src.join(f), d.join(f)).unwrap();
    }
    let p = drs::pcr::read(&latin1(&d.join("rutana_tch.pcr"))).unwrap();
    let mut c = drs::control::Control::from_pcr(&p, &d, "rutana_tch", &d).unwrap();
    assert!(!c.to_text().contains("internal_standard_weight_fraction"), "written only when given");
    c.model.phases[1].internal_standard = Some(0.2);
    let text = c.to_text();
    assert!(text.contains("internal_standard_weight_fraction = 0.2"), "{text}");
    let back = drs::control::Control::from_text(&text).unwrap();
    assert_eq!(back.model.phases[1].internal_standard, Some(0.2));
    assert_eq!(back.model.phases[0].internal_standard, None);
    let drc = d.join("rutana_std.drc");
    std::fs::write(&drc, &text).unwrap();
    let o = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(["refine", drc.to_str().unwrap(), "--cycles", "2"]).output().unwrap();
    let out = String::from_utf8_lossy(&o.stdout);
    assert!(o.status.success(), "{out}{}", String::from_utf8_lossy(&o.stderr));
    let line = out.lines().find(|l| l.contains("amorphous content")).unwrap_or_else(|| panic!("{out}"));
    println!("{line}");
    assert!(line.contains("% of the sample") && line.contains("20 % of phase 2") && !line.contains("negative"), "{line}");
    // Weighed above its refined fraction (0.22), the standard gives a negative content, and the report says why.
    std::fs::write(&drc, text.replace("internal_standard_weight_fraction = 0.2", "internal_standard_weight_fraction = 0.25")).unwrap();
    let o = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(["refine", drc.to_str().unwrap(), "--cycles", "2"]).output().unwrap();
    let out = String::from_utf8_lossy(&o.stdout);
    assert!(out.contains("negative: the standard's weighed fraction"), "{out}");
}

#[test]
fn brindley_table() {
    // The nodes, at 0°, 90° and 180°.
    assert_eq!(brindley_tau(0.0, 45.0), Some(1.0));
    assert!((brindley_tau(-0.5, 0.0).unwrap() - 2.068).abs() < 1e-12);
    assert!((brindley_tau(0.3, 90.0).unwrap() - 0.653).abs() < 1e-12);
    assert!((brindley_tau(0.5, 180.0).unwrap() - 0.529).abs() < 1e-12);
    // Linear between the nodes, in (μ − μ̄)R and in 2θ.
    assert!((brindley_tau(0.25, 90.0).unwrap() - 0.5 * (0.753 + 0.653)).abs() < 1e-12);
    assert!((brindley_tau(0.1, 135.0).unwrap() - 0.5 * (0.866 + 0.870)).abs() < 1e-12);
    // Outside the table.
    assert_eq!(brindley_tau(0.6, 90.0), None);
    assert_eq!(brindley_tau(-0.51, 90.0), None);
    // The paper's eq. 8 against the table where it applies (|x| ≤ 0.1, §2): within 0.0043 at every angle.
    let mut worst = 0.0f64;
    for i in -100..=100 {
        let x = i as f64 * 1e-3;
        let eq8 = 1.0 - 1.450 * x + 1.426 * x * x;
        for tt in [0.0, 45.0, 90.0, 135.0, 180.0] {
            worst = worst.max((eq8 - brindley_tau(x, tt).unwrap()).abs());
        }
    }
    println!("eq. 8 against Brindley's table, |x| ≤ 0.1: largest difference {worst:.4}");
    assert!(worst <= 0.0043, "{worst}");
}
