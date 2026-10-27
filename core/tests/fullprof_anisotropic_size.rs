//! The rest of P25 (L267) against FullProf 8.20's runs of its own test files (L216) and of tutorial 6 (L038): the
//! widths FullProf drew each reflection with (Hw, HwG, HwL, ETA of the reflection list of the calculating pass's
//! `.out`), at the Bragg angle of FullProf's d, against `Model::peak_width` of the phase read from the run's `.pcr`, to
//! half a unit of the printed decimal and 2·10⁻⁵ of the value (as `fullprof_anisotropic_strain`):
//!
//! - `prba10k` (Pr₀.₅Ba₀.₅CoO₃₋ₓ, P 4/m m m, synchrotron λ ≈ 0.44 Å): needles along [0 0 1] (Size-Model −1, LorSiz)
//!   and the tetragonal quartic strain (Strain-Model 4: S_400, S_004, S_220, S_202, with ξ);
//! - `lac` (La₂NiO₄, P 4₂/n c m, D1B): the tetragonal quartic strain of phase 1 (phase 2 is magnetic);
//! - `cuf1k` (CuF₂, P 2₁/n, D1A over d1a_llb.irf): the monoclinic quartic strain (Strain-Model 2) of phase 1 (phase 2
//!   is magnetic);
//! - `Pd3MnD` (Pd₃MnD₀.₈, P m −3 m, 3T2 over its resolution file): the antiphase domains of Cu₃Au (Size-Model −14):
//!   only the reflections of mixed parity broadened.
//!
//! Each check bites: the same rows without the anisotropic term must miss.

use drs::powder::Model;
use std::path::{Path, PathBuf};

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display())).iter().map(|&b| b as char).collect()
}

fn dir(run: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs").join(run)
}

/// The nuclear phases of a run's `.pcr`, under FullProf's rules (L046).
fn model(run: &str, pcr: &str) -> Model {
    let mut p = drs::pcr::read(&latin1(&dir(run).join(format!("{pcr}.pcr")))).unwrap();
    p.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
    Model::from_pcr_fullprof(&p, &dir(run), pcr).unwrap_or_else(|e| panic!("{pcr}: {e}")).model
}

fn half(s: &str) -> f64 {
    0.5 * 10f64.powi(-(s.split_once('.').map_or(0, |(_, d)| d.len()) as i32))
}

/// Phase 1's rows of the `.out`'s last reflection list: (hkl, d, the printed Hw, HwG, HwL, ETA).
fn rows(out: &str) -> Vec<([i32; 3], f64, [String; 4])> {
    let (mut ph, mut hdr) = (usize::MAX, Vec::<String>::new());
    let mut v: Vec<([i32; 3], f64, [String; 4])> = Vec::new();
    for l in out.lines() {
        if l.contains("Phase No.:") && l.contains("Pattern#") {
            ph = l.split("Phase No.:").nth(1).unwrap().split_whitespace().next().unwrap().parse::<usize>().unwrap();
        }
        if l.trim_start().starts_with("No.  Code") {
            hdr = l.split_whitespace().map(String::from).collect();
            continue;
        }
        let f: Vec<&str> = l.split_whitespace().collect();
        if hdr.is_empty() || f.len() + 1 < hdr.len() || f[0].parse::<usize>().is_err() || ph != 1 {
            continue;
        }
        let col = |name: &str| hdr.iter().position(|h| h == name).and_then(|i| f.get(i).copied());
        let (Some(d), Some(h), Some(k), Some(ll)) = (col("d-hkl"), col("H"), col("K"), col("L")) else { continue };
        // A new list starts again at row 1: keep the last.
        if f[0] == "1" {
            v.clear();
        }
        let hkl = [h.parse().unwrap(), k.parse().unwrap(), ll.parse().unwrap()];
        v.push((hkl, d.parse().unwrap(), ["Hw", "HwG", "HwL", "ETA"].map(|n| col(n).unwrap().to_string())));
    }
    v
}

/// The worst row at the tolerance, and how many rows the model `bare` (the anisotropic term left out) misses.
fn compare(m: &Model, bare: &Model, out: &str) -> (usize, f64, String, usize) {
    let (mut n, mut worst, mut at, mut missed) = (0, 0.0f64, String::new(), 0);
    for (hkl, d, printed) in rows(out) {
        let tt = 2.0 * (m.lambda[0] / (2.0 * d)).asin().to_degrees();
        let p = &m.phases[0];
        let (hw, eta) = m.peak_width(p, hkl, tt, 1);
        let (hg, hl) = m.peak_width_components(p, hkl, tt);
        for (s, v, name) in [(&printed[0], hw, "Hw"), (&printed[1], hg, "HwG"), (&printed[2], hl, "HwL"), (&printed[3], eta, "ETA")] {
            let r: f64 = s.parse().unwrap();
            let e = (v - r).abs() / (half(s) + 2e-5 * r.abs());
            if e > worst {
                (worst, at) = (e, format!("{name} of {hkl:?} {v:.6} against {s}"));
            }
        }
        let (hw0, _) = bare.peak_width(&bare.phases[0], hkl, tt, 1);
        if (hw0 - printed[0].parse::<f64>().unwrap()).abs() > half(&printed[0]) + 2e-5 * hw0 {
            missed += 1;
        }
        n += 1;
    }
    (n, worst, at, missed)
}

/// (run, .pcr, .out, the reflections of phase 1 and those the term misses, as measured in L267).
const RUNS: [(&str, &str, &str, usize, usize); 4] = [
    ("prba10k", "prba10k", "prba10k_calc.out", 207, 195),
    ("lac", "lac", "lac_calc.out", 28, 14),
    ("cuf1k", "cuf1k", "cuf1k_calc.out", 57, 57),
    ("Pd3MnD", "pdfinal", "pdfinal_calc.out", 33, 22),
];

#[test]
fn fullprof_anisotropic_size() {
    for (run, pcr, out, rows_expected, missed_expected) in RUNS {
        let m = model(run, pcr);
        let p = &m.phases[0];
        println!("{run}: {} {}; quartic {:?}; size {:?}", p.group.laue, p.name.trim(), p.quartic_strain.as_ref().map(|q| q.terms.clone()), p.shape_size);
        // Without the anisotropic terms: the size alone where there is one, else the strain.
        let mut bare = m.clone();
        if bare.phases[0].shape_size.is_some() {
            bare.phases[0].shape_size = None;
        } else {
            bare.phases[0].quartic_strain = None;
        }
        let (n, worst, at, missed) = compare(&m, &bare, &latin1(&dir(run).join(out)));
        println!("  {n} reflections of phase 1 in {out}: widths at {worst:.2} of the tolerance ({at}); without the anisotropic term {missed} miss");
        assert_eq!((n, missed), (rows_expected, missed_expected), "{run}: rows, and rows missed without the term");
        assert!(worst <= 1.0, "{run}: widths at {worst:.2} of the tolerance: {at}");
    }
}

/// The anisotropic size in the files and the engine (L267): prba10k's needles and Pd3MnD's antiphase domains are
/// free from their codewords (51, 101), written to the control file and read back, exported to a `.pcr` and read back
/// with the same codewords; and the engine refines LorSiz back from 20 % off on a pattern calculated with the run's.
#[test]
fn shape_size_files_and_refinement() {
    for (run, code, kind) in [("prba10k", "prba10k", "\"needle\""), ("Pd3MnD", "pdfinal", "\"antiphase\"")] {
        let d = dir(run);
        let pcr = drs::pcr::read(&latin1(&d.join(format!("{code}.pcr")))).unwrap();
        let loaded = Model::from_pcr(&pcr, &d, code).unwrap();
        let set = drs::refinable::ParameterSet::from_pcr(&pcr, &loaded.model).unwrap();
        assert!(set.free.iter().any(|p| p.name == "phase1.shape_size"), "{run}: LorSiz free");

        // The control file: written, read back, the same size.
        let c = drs::control::Control::from_pcr(&pcr, &d, code, &d).unwrap();
        let text = c.to_text();
        assert!(text.lines().any(|l| l.starts_with("shape_size_kind = ") && l.contains(kind)), "{run}: {}", text.lines().filter(|l| l.contains("shape_size")).collect::<Vec<_>>().join("\n"));
        let back = drs::control::Control::from_text_in(&text, &d).unwrap();
        assert_eq!(back.model.phases[0].shape_size, loaded.model.phases[0].shape_size, "{run}");
        assert_eq!(back.model.phases[0].quartic_strain, loaded.model.phases[0].quartic_strain, "{run}");

        // The .pcr export, read back: the same size, form and free parameters.
        let ex = drs::control::pcr_export::export(&back, "size", loaded.data.as_ref()).unwrap();
        let out = Path::new(env!("CARGO_MANIFEST_DIR")).join("target/test-anisotropic-size").join(run);
        let _ = std::fs::remove_dir_all(&out);
        std::fs::create_dir_all(&out).unwrap();
        std::fs::write(out.join("size.pcr"), &ex.pcr).unwrap();
        for (name, t) in &ex.files {
            std::fs::write(out.join(name), t).unwrap();
        }
        let again = drs::pcr::read(&ex.pcr).unwrap();
        let m2 = Model::from_pcr(&again, &out, "size").unwrap().model;
        assert_eq!(m2.phases[0].shape_size, loaded.model.phases[0].shape_size, "{run}");
        assert_eq!(m2.phases[0].quartic_strain, loaded.model.phases[0].quartic_strain, "{run}");
        let set2 = drs::refinable::ParameterSet::from_pcr(&again, &m2).unwrap();
        assert_eq!(set2.free.iter().map(|p| p.name.clone()).collect::<Vec<_>>(), set.free.iter().map(|p| p.name.clone()).collect::<Vec<_>>(), "{run}");

        // The engine: LorSiz from 20 % off, on noiseless data calculated with the run's value.
        let truth = loaded.model.clone();
        let mut data = loaded.data.clone().unwrap();
        data.y = truth.calculate(&data.x).ycalc;
        data.sigma = data.y.iter().map(|y| y.abs().max(1.0).sqrt()).collect();
        let mut m = truth.clone();
        m.phases[0].shape_size.as_mut().unwrap().value *= 0.8;
        let set = drs::refinable::ParameterSet::from_addresses(&["phase1.shape_size".to_string()], &[], &[], &m).unwrap();
        let r = drs::lsq::refine(&mut m, &data, &set, drs::lsq::Options { max_cycles: 30, eps: 1e-6 }, |_| {});
        let (a, b) = (truth.phases[0].shape_size.unwrap().value, m.phases[0].shape_size.unwrap().value);
        println!("{run}: LorSiz refined back to {b:.6} against {a:.6} ({:?})", r.stop);
        assert!((a - b).abs() <= 1e-5 * a, "{run}: LorSiz {b} against {a} ({:?})", r.stop);
    }
}
