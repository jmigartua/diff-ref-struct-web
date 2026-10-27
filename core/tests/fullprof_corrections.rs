//! Lorentz–polarisation, absorption and preferred orientation against FullProf (test `fullprof_corrections`).
//!
//! FullProf's final reflection lists in the `.out` of its runs (`manuals/fullprof/runs/`, L034, L038) end with a
//! column CORR. Found on the runs (L039): CORR = m · Lp(θ) · A(θ) · P(h), with m the multiplicity, and × Ratio for the
//! λ2 reflections (Code 2). The core's product is compared reflection by reflection, θ from the d-hkl FullProf prints;
//! tolerance: the effect of half a unit of the sixth decimal of d and of the last decimal of CORR, plus 10⁻⁶ of the
//! value (single precision).
//!
//! What the runs exercise: X-ray Lp with a monochromator (Cthm ≠ 0: Y2O3, rutile-anatase, YVO3, LaB6), without
//! (Yb2ZrWO8), neutron Lp, λ2 with Ratio, and cylinder absorption with μR = 0.386 (the NdNi runs; without it the
//! product would be 5–6% off across the pattern), and March–Dollase preferred orientation (Yb2ZrWO8-PO-March, G1 =
//! 0.85 about [1 1 0], L042: P is averaged over the reflections equivalent to h). The exponential form is in no run.
//!
//! Left out: Y2O3, whose pseudo-Voigt list (Npr = 5) prints CORR = 1.000000 for every reflection (another quantity;
//! its Lp and intensities are checked by its profile in `fullprof_profiles`), and the Egyptian make-up (μR = 3.0). Its CORR is not the factor FullProf applies to its profile: its
//! profile, whose peaks carry m·Lp·A, agrees with the core's to 1.6·10⁻⁵ (`fullprof_profiles`), while its CORR
//! differs from m·Lp·A by 6% at 9.6° and by a factor 17 at 75°. What it holds there is not known (listed in L039).
//! Magnetic phases (NdNi-HTO, phase 2) are not calculated by the core yet.

use drs::powder::Model;
use std::path::Path;

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display())).iter().map(|&b| b as char).collect()
}

const RUNS: [(&str, &str); 11] = [
    ("Co-YBCO", "Cocus-asym"),
    ("MgAl2O4-s", "MgAl2O4s_Rietveld"),
    ("rutile-anatase", "rutana_tch"),
    ("YVO3", "YVO3a"),
    ("Yb2ZrWO8", "ybfluo_Riet_IRF"),
    ("NdNi-HTO", "NdNi_HTO"),
    ("NdNi-LTOa", "NdNi_LTOa"),
    ("NdNi-LTOb", "NdNi_LTOb"),
    ("LaB6", "LaB6rx"),
    ("Pd3MnD", "pdfinal"),
    ("Yb2ZrWO8-PO-March", "ybfluo_Riet_IRF"),
];

#[test]
fn fullprof_corrections() {
    let mut fails = Vec::new();
    let mut total = 0;
    for (run, code) in RUNS {
        let d = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs").join(run);
        let mut p = drs::pcr::read(&latin1(&d.join(format!("{code}.pcr")))).unwrap();
        p.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
        // Anisotropic size and strain do not enter this product; the model is built without them.
        for ph in p.phases.iter_mut() {
            ph.flags.insert("Str".into(), 0.0);
            ph.size_model = 0;
            ph.strain_model = 0;
            ph.strain = [0.0; 3];
            ph.lorentz_size = 0.0;
        }
        // The .pcr's values under FullProf's own rules (its X-ray Lp is twice eq. 3.26; L046, L050).
        let m = Model::from_pcr_fullprof(&p, &d, code).unwrap_or_else(|e| panic!("{run}: {e}")).model;
        let out = latin1(&d.join(format!("{code}_calc.out")));
        let (mut ph, mut n, mut worst, mut at) = (usize::MAX, 0, 0.0f64, String::new());
        let mut hdr: Vec<String> = Vec::new();
        for (no, l) in out.lines().enumerate() {
            if l.contains("Phase No.:") && l.contains("Pattern#") {
                ph = l.split("Phase No.:").nth(1).unwrap().split_whitespace().next().unwrap().parse::<usize>().unwrap() - 1;
            }
            if l.trim_start().starts_with("No.  Code") {
                hdr = l.split_whitespace().map(String::from).collect();
                continue;
            }
            let f: Vec<&str> = l.split_whitespace().collect();
            if hdr.is_empty() || ph >= m.phases.len() || f.len() + 1 < hdr.len() || f[0].parse::<usize>().is_err() {
                continue;
            }
            let col = |name: &str| hdr.iter().position(|h| h == name).and_then(|i| f.get(i).copied());
            let (Some(ds), Some(cs)) = (col("d-hkl"), col("CORR")) else { continue };
            let num = |s: &str| s.parse::<f64>().ok();
            let (Some(dh), Some(corr), Some(mult), Some(h), Some(k), Some(l)) = (num(ds), num(cs), num(f[5]), num(f[2]), num(f[3]), num(f[4])) else {
                continue;
            };
            let hkl = [h, k, l];
            let two = f[1] == "2";
            let product = |dd: f64| {
                let th = ((if two { m.lambda[1] } else { m.lambda[0] }) / (2.0 * dd)).asin();
                mult * m.lorentz_polarisation(2.0 * th.to_degrees()) * m.absorption(th) * m.preferred_orientation(&m.phases[ph], hkl) * if two { m.ratio } else { 1.0 }
            };
            let c = product(dh);
            let bound = (product(dh + 5e-7) - c).abs().max((product(dh - 5e-7) - c).abs()) + 0.5 * 10f64.powi(-(cs.split_once('.').map_or(0, |(_, x)| x.len()) as i32)) + 1e-6 * c;
            let e = (c - corr).abs() / bound;
            if e > worst {
                (worst, at) = (e, format!("{} {} {} (line {}): core {c:.6}, FullProf {cs}", f[2], f[3], f[4], no + 1));
            }
            n += 1;
        }
        println!("{run}: {n} reflections, CORR at {worst:.2} of the tolerance ({at})");
        if n == 0 || worst > 1.0 {
            fails.push(format!("{run}: {n} reflections, CORR at {worst:.2} of the tolerance: {at}"));
        }
        total += n;
    }
    println!("{total} reflections compared");
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}
