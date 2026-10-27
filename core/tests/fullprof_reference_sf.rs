//! Structure factors against FullProf's own, for the reference refinements (test `fullprof_reference_sf`).
//!
//! `manuals/fullprof/runs/` holds FullProf's runs of the five reference refinements (L034). With Hkl = −3 FullProf
//! writes, in `CODE_calc.hkl`, the real and imaginary parts of F of every reflection of every phase (conventional
//! cell, f′ and f″ included, temperature factors included, the overall B not). The core reads each phase from the
//! `.pcr` and computes F for the same h k l.
//!
//! The `.pcr` the structures are read from is `CODE.pcr`, not `CODE_calc.pcr` (L035): asked for NCY = 0, FullProf
//! still runs one cycle ("Number of cycles: 1" in `CODE_calc.out`) and writes the shifted values back into
//! `CODE_calc.pcr`, while its `.hkl` and `.prf` belong to the values it read, those of `CODE.pcr`. With
//! `CODE_calc.pcr` anatase's F differ by up to 1.8% at high angle (its Biso of Ti moved by 0.05 Å² in that cycle);
//! with `CODE.pcr` every F agrees to FullProf's print.
//!
//! Tolerance on |F_core − F_FullProf| (complex): half a unit of the fourth decimal FullProf prints in both parts, plus
//! what FullProf's single precision (unit roundoff u = 2⁻²⁴ ≈ 6·10⁻⁸) can do to each term f_j Occ_j T_j exp(2πi h·x)
//! of the sum: the phase 2π h·x carries the rounding of x (u per coordinate, times |h|+|k|+|l|) and of the product
//! (u of the phase), at most 4π(|h|+|k|+|l|)u in all, and f, Occ and T a few u (4u). So the bound is
//! 5·10⁻⁵·√2 + Σ_j |f_j| |Occ_j| N_ops · u·(4π(|h|+|k|+|l|) + 4), a worst case. The worst ratio to it is printed.

use drs::powder::{Intensities, Model};
use std::path::Path;

const RUNS: [(&str, &str); 5] = [
    ("Y2O3", "y2o3_structure"),
    ("Co-YBCO", "Cocus-asym"),
    ("MgAl2O4-s", "MgAl2O4s_Rietveld"),
    ("rutile-anatase", "rutana_tch"),
    ("Egypt", "Egypt_qa_Rietveld"),
];

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display())).iter().map(|&b| b as char).collect()
}

#[test]
fn fullprof_reference_sf() {
    let mut fails = Vec::new();
    let (mut phases, mut total) = (0, 0);
    for (dir, code) in RUNS {
        let d = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs").join(dir);
        let pcr = drs::pcr::read(&latin1(&d.join(format!("{code}.pcr")))).unwrap();
        let model = Model::from_pcr(&pcr, &d, code).unwrap().model;
        let hkl = latin1(&d.join(format!("{code}_calc.hkl")));
        // Per phase: (reflections, worst ratio to the bound, worst |ΔF|, where).
        let mut stats: Vec<(usize, f64, f64, String)> = vec![(0, 0.0, 0.0, String::new()); model.phases.len()];
        let mut ph = usize::MAX;
        for (line_no, l) in hkl.lines().enumerate() {
            if let Some((_, r)) = l.split_once("Phase No.:") {
                ph = r.split_whitespace().next().unwrap().parse::<usize>().unwrap() - 1;
                continue;
            }
            let f: Vec<&str> = l.split_whitespace().collect();
            if f.len() < 8 {
                continue;
            }
            let Ok(h) = f[0].parse::<i32>() else { continue };
            let hkl3 = [h, f[1].parse().unwrap(), f[2].parse().unwrap()];
            let (fr, fi): (f64, f64) = (f[4].parse().unwrap(), f[5].parse().unwrap());
            let Intensities::Structure(st, sc) = &model.phases[ph].intensities else { panic!("{dir}: phase {} has no structure", ph + 1) };
            let c = st.structure_factor(sc, hkl3);
            let s = st.sin_theta_over_lambda(hkl3);
            let nops = (st.group.ops.len() * st.group.centring.len()) as f64;
            let scale: f64 = st.atoms.iter().zip(sc).map(|(a, x)| x.at(s).abs() * a.occ.abs() * nops).sum();
            let u = 2f64.powi(-24);
            let order = hkl3.iter().map(|v| v.abs() as f64).sum::<f64>();
            let bound = 5e-5 * 2f64.sqrt() + scale * u * (4.0 * std::f64::consts::PI * order + 4.0);
            let df = ((c.re - fr).powi(2) + (c.im - fi).powi(2)).sqrt();
            let e = &mut stats[ph];
            e.0 += 1;
            if df / bound > e.1 {
                (e.1, e.2) = (df / bound, df);
                e.3 = format!("{hkl3:?} at line {}: core {:.4} {:+.4}i, FullProf {fr} {fi:+}i", line_no + 1, c.re, c.im);
            }
        }
        for (k, (n, r, df, at)) in stats.iter().enumerate() {
            println!("{dir} phase {} ({}): {n} reflections, worst |ΔF| {df:.1e} at {r:.2} of the bound ({at})", k + 1, model.phases[k].name.trim());
            if *n == 0 || *r > 1.0 {
                fails.push(format!("{dir} phase {}: {n} reflections, worst at {r:.2} of the bound: {at}", k + 1));
            }
            phases += 1;
            total += n;
        }
    }
    println!("{phases} phases, {total} structure factors compared with FullProf's");
    assert_eq!(phases, 11, "the five runs hold 11 phases");
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}
