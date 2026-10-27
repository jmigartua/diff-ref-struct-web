//! Anisotropic strain, Stephens' quartic form, against FullProf (P25, L211).
//!
//! The reference is FullProf 8.20's run of Nd₂NiO₄ at 1.5 K (D2B, `manuals/fullprof/runs/NdNi-LTOa`, from tutorial 6,
//! L038): Pccn, Laue class mmm, Str = 1, Strain-Model 3, S_400 S_040 S_004 S_220 S_202 S_022 = 20.490307 16.256466
//! 0 −36.733093 0 0, ξ = 0.13621, over the resolution file D2B_Old-5min.irf, with an isotropic U and X besides. Its
//! second phase is magnetic (Jbt = −1), which the core does not calculate: the whole pattern (`.prf`) cannot be
//! compared, so the reference is FullProf's own widths, reflection by reflection, of the nuclear phase:
//!
//! 1. the reflection list of `NdNi_LTOa_calc.out` (Hw, HwG, HwL, ETA, the widths FullProf drew each peak with), at the
//!    Bragg angle of FullProf's d, against `Model::peak_width` of the phase read from `NdNi_LTOa.pcr`, to half a unit
//!    of the printed decimal and 2·10⁻⁵ of the value (as `fullprof_microstructure`);
//! 2. the microstructure file `NdNi_LTOa_calc1.mic` (HGo, HGi, HLo, HLi, betaG, betaL, HpVo, eta, Beta_TCH, and the
//!    maximum strain of each reflection), as `fullprof_microstructure` reads it.
//!
//! The isotropic widths alone (the quartic form left out) must miss most rows: the test bites on the form.

use drs::microstructure::{quartic_terms, QuarticStrain};
use drs::powder::Model;
use drs::spacegroup::SpaceGroup;
use std::path::{Path, PathBuf};

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display())).iter().map(|&b| b as char).collect()
}

fn dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs/NdNi-LTOa")
}

/// The nuclear phase of a `.pcr` of the run, under FullProf's rules (λ1 in β*, √|U|; L046).
fn model(pcr: &str) -> Model {
    let mut p = drs::pcr::read(&latin1(&dir().join(format!("{pcr}.pcr")))).unwrap();
    p.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
    Model::from_pcr_fullprof(&p, &dir(), pcr).unwrap_or_else(|e| panic!("{pcr}: {e}")).model
}

fn half(s: &str) -> f64 {
    0.5 * 10f64.powi(-(s.split_once('.').map_or(0, |(_, d)| d.len()) as i32))
}

fn drs(args: &[&str]) -> std::process::Output {
    std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(args).output().unwrap()
}

#[test]
fn fullprof_anisotropic_strain() {
    let m = model("NdNi_LTOa");
    let q = m.phases[0].quartic_strain.clone().expect("the quartic form read");
    let s = |t: [u8; 3]| q.terms.iter().find(|x| x.0 == t).unwrap().1;
    assert_eq!((s([4, 0, 0]), s([0, 4, 0]), s([0, 0, 4]), s([2, 2, 0]), s([2, 0, 2]), s([0, 2, 2])), (20.490307, 16.256466, 0.0, -36.733093, 0.0, 0.0));
    assert_eq!(q.lorentzian, 0.13621);
    let mut iso = m.clone();
    iso.phases[0].quartic_strain = None;

    // 1. The widths of the .out's reflection list, phase 1.
    let out = latin1(&dir().join("NdNi_LTOa_calc.out"));
    let (mut ph, mut hdr) = (usize::MAX, Vec::<String>::new());
    let (mut n, mut worst, mut at, mut missed) = (0, 0.0f64, String::new(), 0);
    for (no, l) in out.lines().enumerate() {
        if l.contains("Phase No.:") && l.contains("Pattern#") {
            ph = l.split("Phase No.:").nth(1).unwrap().split_whitespace().next().unwrap().parse::<usize>().unwrap() - 1;
        }
        if l.trim_start().starts_with("No.  Code") {
            hdr = l.split_whitespace().map(String::from).collect();
            continue;
        }
        let f: Vec<&str> = l.split_whitespace().collect();
        if hdr.is_empty() || f.len() + 1 < hdr.len() || f[0].parse::<usize>().is_err() || ph != 0 {
            continue;
        }
        let col = |name: &str| hdr.iter().position(|h| h == name).and_then(|i| f.get(i).copied());
        let (Some(d), Some(h), Some(k), Some(ll)) = (col("d-hkl"), col("H"), col("K"), col("L")) else { continue };
        let hkl = [h.parse().unwrap(), k.parse().unwrap(), ll.parse().unwrap()];
        let tt = 2.0 * (m.lambda[0] / (2.0 * d.parse::<f64>().unwrap())).asin().to_degrees();
        let p = &m.phases[0];
        let (hw, eta) = m.peak_width(p, hkl, tt, 1);
        let (hg, hl) = m.peak_width_components(p, hkl, tt);
        let mut row = 0.0f64;
        for (name, v) in [("Hw", hw), ("HwG", hg), ("HwL", hl), ("ETA", eta)] {
            let s = col(name).unwrap();
            let r: f64 = s.parse().unwrap();
            let e = (v - r).abs() / (half(s) + 2e-5 * r.abs());
            row = row.max(e);
            if e > worst {
                (worst, at) = (e, format!("{name} of {hkl:?} {v:.6} against {s} (line {})", no + 1));
            }
        }
        // The same row without the quartic form.
        let (hw0, _) = iso.peak_width(&iso.phases[0], hkl, tt, 1);
        let s = col("Hw").unwrap();
        if (hw0 - s.parse::<f64>().unwrap()).abs() > half(s) + 2e-5 * hw0 {
            missed += 1;
        }
        let _ = row;
        n += 1;
    }
    println!("{n} reflections of phase 1 in NdNi_LTOa_calc.out: widths at {worst:.2} of the tolerance ({at}); without the quartic form {missed} of them miss");
    assert!(n > 300, "{n} reflections read");
    assert!(worst <= 1.0, "widths at {worst:.2} of the tolerance: {at}");
    assert!(missed * 10 > n * 9, "the isotropic widths alone should miss nearly every row: {missed} of {n}");

    // 2. The microstructure file. The widths from the parameters read; the maximum strain from those after the cycle
    //    (`NdNi_LTOa_calc.pcr`, as fullprof_microstructure explains: the calc pass runs one cycle, NCY = 1), at the
    //    reflections' d of the cell read. Found here: with the cell after the cycle (a, b, c moved by 4·10⁻⁶ to 10⁻⁵ Å)
    //    the strains of the strongest rows are 1.5·10⁻⁴ above FullProf's print, 2.4 times the effect of the other
    //    parameters' print; with the cell read, every row is within it (L211).
    let mut refined = model("NdNi_LTOa_calc");
    refined.set(drs::jacobian::Field::Cell(0, 0), m.phases[0].cell.a);
    refined.set(drs::jacobian::Field::Cell(0, 1), m.phases[0].cell.b);
    refined.set(drs::jacobian::Field::Cell(0, 2), m.phases[0].cell.c);
    let mic = latin1(&dir().join("NdNi_LTOa_calc1.mic"));
    let (mut n, mut worst, mut at, mut worst_s, mut at_s) = (0, 0.0f64, String::new(), 0.0f64, String::new());
    for (no, l) in mic.lines().enumerate() {
        let f: Vec<&str> = l.split_whitespace().collect();
        if f.len() < 20 || l.starts_with('!') || f[0].parse::<f64>().is_err() {
            continue;
        }
        let hkl = [f[11].parse::<i32>().unwrap(), f[12].parse().unwrap(), f[13].parse().unwrap()];
        let tt: f64 = f[14].parse().unwrap();
        let b = m.breadths(0, hkl, 1).unwrap();
        assert!((b.two_theta - tt).abs() < 2e-3, "{hkl:?} at {tt}: core {:.4}", b.two_theta);
        for (i, v, name) in [(1, b.hg_obs, "HGo"), (2, b.hg_ins, "HGi"), (3, b.hl_obs, "HLo"), (4, b.hl_ins, "HLi"), (5, b.beta_g, "betaG"), (6, b.beta_l, "betaL"), (15, b.h_obs, "HpVo"), (16, b.eta_obs, "eta"), (19, b.beta_tch, "Beta_TCH")] {
            let r: f64 = f[i].parse().unwrap();
            let e = (v - r).abs() / (half(f[i]) + 2e-5 * r.abs());
            if e > worst {
                (worst, at) = (e, format!("{name} of {hkl:?} {v:.5} against {} (line {})", f[i], no + 1));
            }
        }
        // The maximum strain, with the effect of half a unit of the last printed decimal of the refined U, X, S_HKL
        // (six decimals) and ξ (five).
        let after = refined.breadths(0, hkl, 1).unwrap();
        let strain = after.strain.unwrap();
        let mut spread = 0.0;
        let nq = refined.phases[0].quartic_strain.as_ref().unwrap().terms.len();
        for k in 0..nq + 3 {
            let mut m2 = refined.clone();
            let p = &mut m2.phases[0];
            let q = p.quartic_strain.as_mut().unwrap();
            match k {
                0 => p.u += 5e-7,
                1 => p.x += 5e-7,
                2 => q.lorentzian += 5e-6,
                _ => q.terms[k - 3].1 += 5e-7,
            }
            spread += (m2.breadths(0, hkl, 1).unwrap().strain.unwrap() - strain).abs();
        }
        let r: f64 = f[10].parse().unwrap();
        let e = (strain - r).abs() / (half(f[10]) + spread);
        if e > worst_s {
            (worst_s, at_s) = (e, format!("Max-strain of {hkl:?} {strain:.4} against {} (line {})", f[10], no + 1));
        }
        n += 1;
    }
    println!("{n} rows of NdNi_LTOa_calc1.mic: widths and breadths at {worst:.2} of the tolerance ({at}); strain at {worst_s:.2} ({at_s})");
    assert!(n > 300, "{n} rows read");
    assert!(worst <= 1.0 && worst_s <= 1.0, "widths at {worst:.2} ({at}), strain at {worst_s:.2} ({at_s})");

    // The command line: drs mic lists the same rows, the magnetic phase left out and said.
    let o = drs(&["mic", dir().join("NdNi_LTOa.pcr").to_str().unwrap(), "--fullprof"]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(String::from_utf8_lossy(&o.stderr).contains("phase 2 is magnetic (Jbt = -1): left out"));
    let text = String::from_utf8_lossy(&o.stdout);
    let row = text.lines().find(|l| l.split_whitespace().skip(9).take(3).collect::<Vec<_>>() == ["2", "0", "0"]).expect("the (2 0 0) row");
    assert_eq!(row.split_whitespace().nth(1), Some("0.5899"), "HGo of (2 0 0) as FullProf prints it: {row}");
}

/// The terms each Laue class allows, from the group's rotations (the invariants of its orbit sums, L267), against
/// FullProf's LINE 37 of each Strain-Model (manual pp. 119–120); each form takes one value over a reflection's
/// equivalents.
#[test]
fn quartic_terms_by_symmetry() {
    let names = |g: &str| quartic_terms(&SpaceGroup::from_fullprof_symbol(g).unwrap()).map(|t| t.iter().map(|t| drs::microstructure::term_name(*t)).collect::<Vec<_>>().join(" "));
    assert_eq!(names("P c c n").unwrap(), "S_400 S_040 S_004 S_220 S_202 S_022");
    assert_eq!(names("P 1 21/c 1").unwrap(), "S_400 S_040 S_004 S_220 S_202 S_022 S_121 S_301 S_103");
    assert_eq!(names("P 1 1 21/b").unwrap(), "S_400 S_040 S_004 S_220 S_202 S_022 S_112 S_310 S_130");
    assert_eq!(names("P -1").unwrap().split(' ').count(), 15);
    // As FullProf's Strain-Models 5, 7, 8–12 and 13–14.
    for (g, terms) in [
        ("P 4/m m m", "S_400 S_004 S_220 S_202"),
        ("R -3 m", "S_400 S_004 S_112 S_211"),
        ("P -3 m 1", "S_400 S_004 S_112 S_211"),
        ("P 6/m", "S_400 S_004 S_112"),
        ("P 63/m m c", "S_400 S_004 S_112"),
        ("F m -3 m", "S_400 S_220"),
        ("P a -3", "S_400 S_220"),
    ] {
        assert_eq!(names(g).unwrap(), terms, "{g}");
    }
    // Where FullProf leaves an invariant out: 4/m's S_310 (h³k − hk³; Strain-Model 4 reads 4/mmm's four) and −3's S_121
    // (Strain-Model 6 reads −3m's four, and 8 the hexagonal three).
    assert_eq!(names("P 4/m").unwrap(), "S_400 S_004 S_220 S_202 S_310");
    assert_eq!(names("I 41/a").unwrap(), "S_400 S_004 S_220 S_202 S_310");
    assert_eq!(names("R -3").unwrap(), "S_400 S_004 S_112 S_211 S_121");
    for g in ["P 4/m", "R -3 m", "R -3", "P 63/m m c", "F m -3 m", "P a -3"] {
        let sg = SpaceGroup::from_fullprof_symbol(g).unwrap();
        let terms = quartic_terms(&sg).unwrap();
        let q = QuarticStrain::for_group(&sg, terms.iter().enumerate().map(|(i, t)| (*t, 1.0 + i as f64)).collect(), 0.0);
        for hkl in [[3, 1, 2], [2, -1, 4], [1, 2, 3], [4, 1, 0]] {
            let v = q.sum(hkl.map(|x| x as f64));
            for e in drs::reflections::equivalents(&sg, hkl) {
                assert!((q.sum(e.map(|x| x as f64)) - v).abs() <= 1e-9 * v.abs(), "{g}: {hkl:?} and {e:?}");
            }
        }
    }
    assert_eq!(drs::microstructure::fullprof_model_of(&quartic_terms(&SpaceGroup::from_fullprof_symbol("P 4/m m m").unwrap()).unwrap()), Some(5));
    // FullProf's Strain-Model of each class read, and its order.
    assert_eq!(drs::microstructure::fullprof_model_of(&quartic_terms(&SpaceGroup::from_fullprof_symbol("P c c n").unwrap()).unwrap()), Some(3));
    assert_eq!(drs::microstructure::fullprof_model_of(&quartic_terms(&SpaceGroup::from_fullprof_symbol("P 1 21/c 1").unwrap()).unwrap()), Some(2));
    assert_eq!(drs::microstructure::fullprof_model_of(&quartic_terms(&SpaceGroup::from_fullprof_symbol("P -1").unwrap()).unwrap()), None);
    // Equivalent reflections of Pccn have the same width; inequivalent ones of the same d need not.
    let m = model("NdNi_LTOa");
    let p = &m.phases[0];
    for hkl in [[2, 1, 3], [1, 2, 3], [3, 0, 4]] {
        let tt = p.cell.two_theta(hkl.map(|v| v as f64), m.lambda[0]).unwrap();
        let w = m.peak_width(p, hkl, tt, 1);
        for e in drs::reflections::equivalents(&p.group, hkl) {
            assert_eq!(m.peak_width(p, e, tt, 1), w, "{hkl:?} and {e:?}");
        }
    }
    // The form by hand for (2 0 0): σ²(M) = 10⁻⁸ S_400 2⁴, M = 4/a², D_ST = √(8 ln 2)(180/π)σ/M.
    let q = p.quartic_strain.as_ref().unwrap();
    let a = p.cell.a;
    let dst = (8.0 * 2f64.ln()).sqrt() * (180.0 / std::f64::consts::PI) * (1e-8 * 20.490307 * 16.0f64).sqrt() / (4.0 / (a * a));
    assert!((q.width([2.0, 0.0, 0.0], 4.0 / (a * a)) - dst).abs() < 1e-12 * dst);
    let _ = QuarticStrain { terms: Vec::new(), forms: Vec::new(), lorentzian: 0.0 };
}

/// The quartic form through the control file and the `.pcr` export, its codewords kept; the engine recovers it from a
/// pattern the core calculates, and the classes not read are refused with their reason.
#[test]
fn quartic_strain_files_and_refinement() {
    let d = dir();
    let mut pcr = drs::pcr::read(&latin1(&d.join("NdNi_LTOa.pcr"))).unwrap();
    pcr.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
    let loaded = Model::from_pcr(&pcr, &d, "NdNi_LTOa").unwrap();
    let set = drs::refinable::ParameterSet::from_pcr(&pcr, &loaded.model).unwrap();
    let names: Vec<&str> = set.free.iter().map(|p| p.name.as_str()).collect();
    for n in ["phase1.quartic_strain.400", "phase1.quartic_strain.040", "phase1.quartic_strain.220", "phase1.quartic_strain_lorentzian"] {
        assert!(names.contains(&n), "{n} free (codes 201, 151, 221, 211): {names:?}");
    }
    assert!(!names.contains(&"phase1.quartic_strain.004"), "S_004 has code 0");

    // The control file: written, read back, the same model.
    let c = drs::control::Control::from_pcr(&pcr, &d, "NdNi_LTOa", &d).unwrap();
    let text = c.to_text();
    assert!(text.contains("quartic_strain = [") && text.contains("[4, 0, 0, 20.490307]"), "{}", text.lines().filter(|l| l.contains("quartic")).collect::<Vec<_>>().join("\n"));
    let back = drs::control::Control::from_text_in(&text, &d).unwrap();
    assert_eq!(back.model.phases[0].quartic_strain, loaded.model.phases[0].quartic_strain);
    // A term the group does not allow, and the share without a form, are refused with their place.
    let wrong = text.replace("[2, 0, 2, 0.0]", "[2, 1, 1, 0.0]");
    let e = drs::control::Control::from_text_in(&wrong, &d).unwrap_err();
    assert!(e.iter().any(|x| x.message.contains("S_211 is not allowed by the phase's Laue class mmm")), "{e:?}");

    // The .pcr export, read back: the same S_HKL, ξ and codewords.
    let ex = drs::control::pcr_export::export(&back, "ndni", loaded.data.as_ref()).unwrap();
    let out = Path::new(env!("CARGO_MANIFEST_DIR")).join("target/test-anisotropic-strain");
    let _ = std::fs::remove_dir_all(&out);
    std::fs::create_dir_all(&out).unwrap();
    std::fs::write(out.join("ndni.pcr"), &ex.pcr).unwrap();
    for (name, t) in &ex.files {
        std::fs::write(out.join(name), t).unwrap();
    }
    let again = drs::pcr::read(&ex.pcr).unwrap();
    assert_eq!((again.phases[0].flag("Str"), again.phases[0].strain_model), (1, 3));
    let m2 = Model::from_pcr(&again, &out, "ndni").unwrap().model;
    assert_eq!(m2.phases[0].quartic_strain, loaded.model.phases[0].quartic_strain);
    let set2 = drs::refinable::ParameterSet::from_pcr(&again, &m2).unwrap();
    assert_eq!(set2.free.iter().map(|p| p.name.clone()).collect::<Vec<_>>(), set.free.iter().map(|p| p.name.clone()).collect::<Vec<_>>());

    // The engine: a pattern calculated with the run's form, the S_HKL and ξ moved by 20 %, refined back. The data
    // have no noise, so χ² goes to 0 with the esd's, and the divergence guard stops the run when the last shift is
    // many of those vanishing esd's: the values are then those of the cycle before, to 10⁻⁵.
    let data0 = loaded.data.clone().unwrap();
    let truth = loaded.model.clone();
    let mut data = data0.clone();
    data.y = truth.calculate(&data.x).ycalc;
    data.sigma = data.y.iter().map(|y| y.abs().max(1.0).sqrt()).collect();
    let mut m = truth.clone();
    let q = m.phases[0].quartic_strain.as_mut().unwrap();
    for t in q.terms.iter_mut() {
        t.1 *= 0.8;
    }
    q.lorentzian *= 1.2;
    let free: Vec<String> = ["400", "040", "220"].iter().map(|t| format!("phase1.quartic_strain.{t}")).chain(["phase1.quartic_strain_lorentzian".to_string()]).collect();
    let set = drs::refinable::ParameterSet::from_addresses(&free, &[], &[], &m).unwrap();
    let r = drs::lsq::refine(&mut m, &data, &set, drs::lsq::Options { max_cycles: 30, eps: 1e-6 }, |_| {});
    let (qt, qm) = (truth.phases[0].quartic_strain.as_ref().unwrap(), m.phases[0].quartic_strain.as_ref().unwrap());
    for ((t, a), (_, b)) in qt.terms.iter().zip(&qm.terms) {
        assert!((a - b).abs() <= 1e-5 * a.abs().max(1.0), "S_{t:?}: {b} against {a} ({:?})", r.stop);
    }
    assert!((qt.lorentzian - qm.lorentzian).abs() < 1e-5 * qt.lorentzian, "ξ {} against {}", qm.lorentzian, qt.lorentzian);
    println!("refined back: S_400 {:.6}, ξ {:.6} ({:?})", qm.terms[0].1, qm.lorentzian, r.stop);

    // Refused, with the reason: NdNi-LTOb's strain model 6 (Str = 0), highpress's hkl-dependent shifts (P26; Pd3MnD's
    // size model −14 is calculated since L267).
    let refuse = |run: &str, code: &str| {
        let d = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs").join(run);
        let mut p = drs::pcr::read(&latin1(&d.join(format!("{code}.pcr")))).unwrap();
        p.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
        Model::from_pcr(&p, &d, code).err().unwrap()
    };
    assert!(refuse("NdNi-LTOb", "NdNi_LTOb").contains("strain model 6 (Str = 0) is described in words"));
    assert!(refuse("highpress", "highpress").contains("hkl-dependent peak shifts (Sol = 1"));
}
