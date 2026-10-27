//! FullProf's sample displacement (L267): FullProf 8.20 shifts a reflection by Sycos·cos 2θ, not by its manual's
//! Sycos·cos θ (eq. 3.37, p. 56); `Convention::displacement` is that rule. The two runs with Sycos ≠ 0 show it:
//! - the ticks of `NdNi_HTO_calc.prf` and `cuf1k_calc.prf` are their `.hkl`'s 2T shifted by Zero + Sycos·cos 2θ +
//!   Sysin·sin 2θ to the print (`fullprof_outputs::hkl_files`); a fit of the three terms to the ticks gives back line 2's
//!   values (−0.12737, 0.14539, 0.06421 and −0.12701, −0.00316, 0.13128) within 10⁻⁴°, where cos θ leaves 0.019° and
//!   2·10⁻⁴°;
//! - the profile of NdNi_HTO's nuclear phase, FullProf's `NdNi_HTO_calc1.sub` (the phase and the background, in whole
//!   counts), against the core's under FullProf's rules: here, and 5 % away with cos θ.
//!
//! What is left, frozen (L037): 203 of the 2990 points differ by more than the print (half a count), by up to 1.95
//! counts in about 120; 193 of them lie within 0.15° of a peak window's edge, where this run's windows (Lorentzian strain
//! X and the D2B resolution file) are not FullProf's to the point. Not resolved.

use drs::formats::fullprof_outputs::read_prf;
use drs::powder::Model;
use std::path::Path;

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap().iter().map(|&b| b as char).collect()
}

/// (Σ(y − y_FullProf)² / Σ(y_FullProf − background)², the largest |y − y_FullProf|, points beyond half a count).
fn against_sub(m: &Model, x: &[f64], sub: &[f64]) -> (f64, f64, usize) {
    let c = m.calculate(x);
    let (mut s2, mut n2, mut worst, mut beyond) = (0.0, 0.0, 0.0f64, 0);
    for ((p, b), s) in c.per_phase[0].iter().zip(&c.background).zip(sub) {
        let r = p + b - s;
        s2 += r * r;
        n2 += (s - b).powi(2);
        worst = worst.max(r.abs());
        beyond += usize::from(r.abs() > 0.5);
    }
    (s2 / n2, worst, beyond)
}

#[test]
fn fullprof_displacement() {
    let d = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs/NdNi-HTO");
    let mut p = drs::pcr::read(&latin1(&d.join("NdNi_HTO_calc.pcr"))).unwrap();
    // The magnetic phase (Jbt = 1) is not calculated; FullProf's .sub of phase 1 holds the nuclear phase alone.
    p.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
    let mut m = Model::from_pcr_fullprof(&p, &d, "NdNi_HTO_calc").unwrap().model;
    assert!(m.convention.displacement && m.sycos == 0.14539);
    let prf = read_prf(&latin1(&d.join("NdNi_HTO_calc.prf"))).unwrap();
    let text = latin1(&d.join("NdNi_HTO_calc1.sub"));
    assert!(text.lines().next().unwrap().contains("Phase No:   1"));
    let sub: Vec<f64> = text.lines().skip(1).flat_map(|l| l.split_whitespace().map(|v| v.parse::<f64>().unwrap()).collect::<Vec<_>>()).collect();
    assert_eq!(sub.len(), prf.x.len());
    let (rel, worst, beyond) = against_sub(&m, &prf.x, &sub);
    println!("NdNi_HTO phase 1 with Sycos·cos 2θ: relative {rel:.3e}, worst {worst:.2} counts, {beyond} of {} points beyond the print", prf.x.len());
    assert!(rel <= 4.22e-6 && worst <= 1.95 && beyond <= 203, "frozen at 4.22e-6, 1.95 and 203 (L267)");
    // It bites: the manual's cos θ.
    m.convention.displacement = false;
    let (rel, worst, _) = against_sub(&m, &prf.x, &sub);
    println!("  with Sycos·cos θ: relative {rel:.3e}, worst {worst:.0} counts");
    assert!(rel > 0.05 && worst > 400.0);
    // The translation into the default rules keeps FullProf's where Sycos ≠ 0, and drops it where Sycos = 0.
    let mut t = m.clone();
    t.convention.displacement = true;
    t.translate_from_fullprof(None);
    assert!(t.convention.displacement);
    t.sycos = 0.0;
    t.translate_from_fullprof(None);
    assert!(!t.convention.displacement);
}
