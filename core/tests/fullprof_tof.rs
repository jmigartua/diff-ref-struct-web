//! Time-of-flight patterns (P05) in the engine (L341). Their calculated patterns are checked against FullProf's
//! `_calc.prf` by `fullprof_profiles`; here:
//!
//! - `tof_derivatives`: `Model::derivatives` of every refinable field of seven TOF runs (Npr 9 with each absorption
//!   function, with the Q terms, with a resolution file and a RALF bank's Chebyshev background, two phases, and
//!   Hodges' Npr 10 of gabriela's first pattern) against central finite differences of `Model::calculate`, as
//!   `engine_derivatives` checks the constant-wavelength fields: steps h and h/2, h = 10⁻⁵·max(|value|, 10⁻²) (the
//!   scale 10⁻⁵·value, the cell 10⁻⁶·value), Richardson's extrapolation, 300 points; the error, the largest
//!   |analytic − difference| over the largest of either, below 2·10⁻⁵. The models are taken in the default convention
//!   at the run's own Wdt: the window's tapered edge (L058) is smooth, so that it may lie between the evaluations and
//!   inside the data, and its derivative is checked with the rest (a window of 10⁴ H, as `engine_derivatives` widens
//!   them, reaches every point from each of arg_si's 1500 reflections: too slow). Also: the dual-number path
//!   reproduces `calculate` to 10⁻¹⁰.
//! - `tof_refinement`: arg_si's refinement as its `.pcr` codes it (16 parameters), from the values FullProf started
//!   from (`arg_si.pcr.in`), run by the core in FullProf's convention to the `.pcr`'s Eps: its χ² at or below the
//!   core's χ² at FullProf's final values (`arg_si.pcr` as FullProf rewrote it), and every parameter within FullProf's
//!   esd of FullProf's value (`arg_si.out`); the esd's ratio to FullProf's printed. The core's χ², 2.675, is 1.5 % above
//!   FullProf's 2.637: the 586 window-edge points of arg_si, where FullProf's TOF window rule was not found and the
//!   difference is frozen in `fullprof_profiles` (which compares χ² with FullProf's Ycal there, to 2.6370). The esd's
//!   are FullProf's times one ratio, 0.726 = √(2.675/5.08): FullProf's esd's take the χ² over the points with Bragg
//!   contributions, 5.08 (as `esd_profile` found for constant wavelength), the core's the χ² over all points. ExtTOF's
//!   departs: its ratio is 1.908 times the others' (FullProf's esd of Extinc smaller), not explained; frozen (L341),
//!   accepted by Josu as it is (L343).
//! - `tof_results`: P46 (below).
//! - `tof_checks`: the app's request 46 (L354). Each single-pattern TOF run's control file, written in another folder,
//!   is checked (`check::check_path`, which Refine runs before every Run) with no error and none of the diagnostics
//!   that read a TOF bank as constant wavelength (a wavelength, the data's range in 2θ, a unit of 2θ for its zero or
//!   excluded regions); and opened (`control::open`) to the pattern its `.pcr` gives, on the same axis: a GSAS bank's
//!   CONST times ×100 in µs, an `xys` file's times as written (alum and lamn_pol were read ×100 before, with no
//!   reflection in the data).

use drs::jacobian::Field;
use drs::powder::{Intensities, Model};
use std::path::{Path, PathBuf};

mod common;

fn runs() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs")
}

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display())).iter().map(|&b| b as char).collect()
}

/// The default model of pattern k of `folder/file`.
fn load(folder: &str, file: &str, k: usize) -> drs::powder::Loaded {
    let d = runs().join(folder);
    let p = drs::pcr::read_all(&latin1(&d.join(file))).unwrap().remove(k);
    let stem = file.split('.').next().unwrap();
    Model::from_pcr(&p, &d, stem).unwrap_or_else(|e| panic!("{folder}/{file}: {e}"))
}

const RUNS: [(&str, &str, usize); 8] = [
    ("arg_si", "arg_si.pcr", 0),
    ("arg_si", "arg_si_abs1.pcr", 0),
    ("arg_si", "arg_si_abs3.pcr", 0),
    ("alum", "alum.pcr", 0),
    ("CeO2_PEARL", "CeO2_PEARL.pcr", 0),
    ("lamn_pol", "lamn_pol_irf.pcr", 0),
    ("nac_osiris", "nac-osiris.pcr", 0),
    ("gabriela", "gabriela.pcr", 0),
];

/// Every refinable field of the model (`Model::all_fields`), the atoms' of the first two atoms of each phase.
fn fields(m: &Model) -> Vec<Field> {
    m.all_fields()
        .into_iter()
        .filter(|f| match *f {
            Field::AtomX(_, a, _) | Field::AtomB(_, a) | Field::AtomOcc(_, a) | Field::AtomBeta(_, a, _) => a < 2,
            _ => true,
        })
        .collect()
}

#[test]
fn tof_derivatives() {
    let mut fails = Vec::new();
    let mut checked = 0;
    let mut kinds = std::collections::BTreeSet::new();
    for (folder, file, k) in RUNS {
        let loaded = load(folder, file, k);
        let data = loaded.data.as_ref().expect("data");
        let mut m = loaded.model.clone();
        assert!(m.tof.is_some(), "{file}: not a TOF pattern");
        if m.phases.iter().all(|ph| matches!(ph.intensities, Intensities::Structure(..))) {
            let (a, b) = (m.calculate(&data.x).ycalc, m.calculate_dual_values(&data.x));
            let worst = a.iter().zip(&b).map(|(u, v)| (u - v).abs() / u.abs().max(1.0)).fold(0.0, f64::max);
            if worst > 1e-10 {
                fails.push(format!("{file}: the dual path differs from calculate by {worst:.2e}"));
            }
        }
        m.convention = drs::powder::Convention::DEFAULT;
        let n = data.x.len();
        let x: Vec<f64> = (0..300).map(|i| data.x[i * (n - 1) / 299]).collect();
        let fl = fields(&m);
        let (c0, jac) = m.derivatives(&x, &fl);
        let ymax = c0.ycalc.iter().map(|v| v.abs()).fold(0.0, f64::max);
        let mut worst = (0.0f64, String::new());
        for (j, &f) in fl.iter().enumerate() {
            let v = m.get(f);
            let h = match f {
                Field::Scale(_) => 1e-5 * v.abs(),
                Field::Cell(..) => 1e-6 * v.abs(),
                _ => 1e-5 * v.abs().max(1e-2),
            };
            let central = |h: f64| -> Vec<f64> {
                let mut up = m.clone();
                up.set(f, v + h);
                let mut dn = m.clone();
                dn.set(f, v - h);
                let (yu, yd) = (up.calculate(&x).ycalc, dn.calculate(&x).ycalc);
                yu.iter().zip(&yd).map(|(a, b)| (a - b) / (2.0 * h)).collect()
            };
            let (d1, d2) = (central(h), central(h / 2.0));
            let fd: Vec<f64> = d1.iter().zip(&d2).map(|(a, b)| (4.0 * b - a) / 3.0).collect();
            let top = fd.iter().chain(&jac[j]).map(|v| v.abs()).fold(0.0, f64::max);
            checked += 1;
            if top * h < 1e-7 * ymax * h.max(1e-2) / 1e-2 || top < 1e-7 * ymax {
                continue;
            }
            kinds.insert(drs::refinable::dictionary_id(f, &m));
            let err = jac[j].iter().zip(&fd).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max) / top;
            if err > worst.0 {
                worst = (err, format!("{f:?}"));
            }
            if err > 2e-5 {
                fails.push(format!("{file}: {f:?}: analytic against finite difference {err:.2e} (largest derivative {top:.3e})"));
            }
        }
        println!("{file:22} {} fields, worst {:.1e} ({})", fl.len(), worst.0, worst.1);
    }
    println!("{checked} derivatives checked; non-zero ones of {}", kinds.into_iter().collect::<Vec<_>>().join(", "));
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}

/// P46: a TOF pattern's results file as Publish reads it. arg_si (Npr 9) and lamn_pol_irf (its resolution file's part
/// apart) written and read back: the abscissa in µs, the bank's angle, the TOF–d equation, each TOF value with the
/// resolution file's part and the total the calculation takes, the reflections at their times of flight, and the d of
/// each reflection's time of flight by the inverse of T(d) equal to its d to 10⁻⁹ Å. gabriela (a TOF bank and two
/// constant-wavelength patterns) as a joint results file, a block per pattern with its own abscissa; FullProf's
/// four-bank FP_tests/Ho2BaCuO5/TOF/test-tof.pcr refused for its magnetic phase, with the reason.
#[test]
fn tof_results() {
    use drs::formats::Axis;
    for (folder, file, instrument) in [("arg_si", "arg_si.pcr", false), ("lamn_pol", "lamn_pol_irf.pcr", true)] {
        let d = runs().join(folder);
        let p = drs::pcr::read(&latin1(&d.join(file))).unwrap();
        let stem = file.trim_end_matches(".pcr");
        let c = drs::control::Control::from_pcr(&p, &d, stem, &d).unwrap();
        let loaded = Model::from_pcr(&p, &d, stem).unwrap();
        let data = loaded.data.as_ref().unwrap();
        let (text, _) = drs::results::write(&c, data, None, None, stem).unwrap();
        let r = drs::results::read(&text).unwrap();
        let bank = c.model.tof.as_ref().unwrap();
        assert_eq!(r.axis, Axis::TimeOfFlight, "{file}");
        assert_eq!(r.x.len(), data.x.len());
        assert!(r.x.iter().zip(&data.x).all(|(a, b)| (a - b).abs() <= 1e-6 * b.abs()), "{file}: the abscissa");
        assert_eq!(r.two_theta_fixed, Some(bank.two_theta), "{file}");
        assert!(r.conversion.as_deref().is_some_and(|e| e.starts_with("TOF = Zero + Dtt1 d")), "{file}: {:?}", r.conversion);
        let dtt1 = r.tof_parameters.iter().find(|t| t.id == "tof_dtt1").unwrap();
        assert_eq!(dtt1.value, bank.dtt1, "{file}");
        assert_eq!(dtt1.instrument, bank.instrument.as_ref().map_or(0.0, |i| i.dtt[0]), "{file}");
        assert!((dtt1.total - (dtt1.value + dtt1.instrument)).abs() <= 1e-6 * dtt1.total.abs(), "{file}");
        assert_eq!(dtt1.instrument != 0.0, instrument, "{file}: the resolution file's part");
        assert_eq!(r.tof_parameters.iter().filter(|t| t.phase == Some(0)).count(), 15, "{file}: the phase's 15 TOF values (19 less Hodges' four)");
        // The reflections at their times of flight, and their d back from them.
        let c0 = c.model.calculate(&data.x);
        assert!(!r.reflections.is_empty());
        let mut worst = 0f64;
        for (_, _, hkl, _, v) in &r.reflections {
            let pk = c0.peaks.iter().find(|q| q.hkl == *hkl).unwrap();
            assert!((v[1] - pk.two_theta).abs() <= 1e-6 * pk.two_theta, "{file} {hkl:?}: TOF {} against {}", v[1], pk.two_theta);
            let back = c.model.d_spacing_at(pk.two_theta).unwrap();
            worst = worst.max((back - v[0]).abs());
        }
        assert!(worst < 1e-6, "{file}: d from T(d) inverted, worst {worst:.2e} Å against the reflections' d (printed to 7 digits)");
        let exact = c0.peaks.iter().map(|pk| (c.model.d_spacing_at(pk.two_theta).unwrap() - c.model.phases[0].cell.d_spacing(pk.hkl.map(|h| h as f64)).unwrap()).abs()).fold(0.0, f64::max);
        assert!(exact < 1e-9, "{file}: the inverse of T(d), worst {exact:.2e} Å");
        let q = r.abscissa_in(Axis::Q).unwrap();
        let dd = r.abscissa_in(Axis::D).unwrap();
        assert!(q.iter().zip(&dd).all(|(q, d)| (q * d - 2.0 * std::f64::consts::PI).abs() < 1e-9), "{file}: Q = 2π/d");
        println!("{file}: {} points in µs, {} reflections, {} TOF values; d from T(d) to {exact:.1e} Å", r.x.len(), r.reflections.len(), r.tof_parameters.len());
    }
    // A joint refinement with a TOF bank: each pattern's block with its abscissa.
    let d = runs().join("gabriela");
    let j = drs::joint::Joint::from_pcr(&latin1(&d.join("gabriela.pcr")), &d, "gabriela", false).unwrap();
    let text = j.results(None, "gabriela", "gabriela.pcr").unwrap();
    let r = drs::joint::read_results(&text).unwrap();
    let axes: Vec<Axis> = r.patterns.iter().map(|p| p.axis).collect();
    assert_eq!(axes, [Axis::TimeOfFlight, Axis::TwoTheta, Axis::TwoTheta], "gabriela's patterns");
    assert!(r.patterns[0].two_theta_fixed.is_some_and(|t| (t - 144.845).abs() < 1e-9));
    // Ho2BaCuO5's four banks: its first phase is magnetic (a magnetic space group, Jbt = −6).
    let d = runs().join("../FP_tests/Ho2BaCuO5/TOF");
    let Err(e) = drs::joint::Joint::from_pcr(&latin1(&d.join("test-tof.pcr")), &d, "test-tof", false) else { panic!("Ho2BaCuO5's test-tof.pcr not refused") };
    assert!(e.contains("magnetic space group P2_1'/c") && e.contains("magnetic scattering is not calculated"), "{e}");
    println!("Ho2BaCuO5 test-tof.pcr: {e}");
}

#[test]
fn tof_refinement() {
    let d = runs().join("arg_si");
    let p = drs::pcr::read(&latin1(&d.join("arg_si.pcr.in"))).unwrap();
    let mut l = Model::from_pcr_fullprof(&p, &d, "arg_si").unwrap();
    let data = l.data.clone().unwrap();
    let set = drs::refinable::ParameterSet::from_pcr(&p, &l.model).unwrap();
    assert_eq!(set.free.len(), 16);
    let r = drs::lsq::refine(&mut l.model, &data, &set, drs::lsq::Options { max_cycles: 30, eps: p.eps }, |_| {});
    // The core's χ² at FullProf's final values.
    let pf = drs::pcr::read(&latin1(&d.join("arg_si.pcr"))).unwrap();
    let at_fp = Model::from_pcr_fullprof(&pf, &d, "arg_si").unwrap().model;
    let chi2_fp = drs::lsq::statistics(&at_fp, &data, &at_fp.calculate(&data.x).ycalc, 16).chi2;
    assert!(r.stats.chi2 <= chi2_fp * (1.0 + 1e-9), "χ² {} above the core's {chi2_fp} at FullProf's values", r.stats.chi2);
    let fp = common::fullprof_params(&latin1(&d.join("arg_si.out")));
    assert_eq!(fp.len(), 16);
    let esd = |f: Field| r.fields.iter().find(|x| x.0 == drs::refinable::address(f, &l.model)).map(|x| x.2).unwrap();
    let mut ratios: Vec<(String, f64)> = Vec::new();
    for (name, (v, s)) in &fp {
        let f = drs::output_names::fullprof_field(name, &l.model).unwrap_or_else(|| panic!("{name}"));
        let mine = l.model.get(f);
        println!("{name:24} FullProf {v:>14.6} ({s:.2e})  core {mine:>14.6} ({:.2e})  Δ/σ {:+.2}  esd ratio {:.3}", esd(f), (mine - v) / s, esd(f) / s);
        assert!((mine - v).abs() <= *s, "{name}: core {mine}, FullProf {v} ± {s}");
        ratios.push((name.clone(), esd(f) / s));
    }
    let mut sorted: Vec<f64> = ratios.iter().map(|r| r.1).collect();
    sorted.sort_by(f64::total_cmp);
    let median = sorted[sorted.len() / 2];
    for (name, ratio) in &ratios {
        let rel = ratio / median;
        if name.starts_with("ExtTOF") {
            assert!((1.90..=1.915).contains(&rel), "{name}: esd ratio {rel:.4} × the others' (frozen at 1.908)");
        } else {
            assert!((rel - 1.0).abs() < 5e-3, "{name}: esd ratio {ratio:.4}, the others' {median:.4}");
        }
    }
    let chi2_bragg = 5.08;
    assert!((median - (r.stats.chi2 / chi2_bragg).sqrt()).abs() < 3e-3, "the esd ratio {median:.4} against √(χ²/χ²_Bragg)");
    println!("arg_si refined in {} cycles: χ² {:.4}, the core's at FullProf's values {chi2_fp:.4} (FullProf's own 2.637); {}", r.cycles, r.stats.chi2, r.stop);
}

#[test]
fn tof_checks() {
    let tmp = Path::new(env!("CARGO_MANIFEST_DIR")).join("target/test-tof-checks");
    let _ = std::fs::remove_dir_all(&tmp);
    std::fs::create_dir_all(&tmp).unwrap();
    let mut checked = 0;
    for (folder, file) in [("arg_si", "arg_si.pcr"), ("arg_si", "arg_si_abs1.pcr"), ("arg_si", "arg_si_abs3.pcr"), ("arg_si", "arg_si_invd.pcr"), ("alum", "alum.pcr"), ("CeO2_PEARL", "CeO2_PEARL.pcr"), ("lamn_pol", "lamn_pol.pcr"), ("lamn_pol", "lamn_pol_irf.pcr"), ("nac_osiris", "nac-osiris.pcr")] {
        let d = runs().join(folder);
        let stem = file.trim_end_matches(".pcr");
        let p = drs::pcr::read(&latin1(&d.join(file))).unwrap();
        let l = Model::from_pcr(&p, &d, stem).unwrap_or_else(|e| panic!("{file}: {e}"));
        let c = drs::control::Control::from_pcr(&p, &d, stem, &tmp).unwrap();
        let path = tmp.join(format!("{stem}.drc"));
        std::fs::write(&path, c.to_text()).unwrap();
        let diags = drs::check::check_path(&path, &drs::check::Options::default());
        let bad: Vec<String> = diags.iter().filter(|x| x.severity == drs::check::Severity::Error || matches!(x.code, "wavelength" | "range-data" | "unit")).map(|x| format!("{}: {}", x.code, x.what)).collect();
        assert!(bad.is_empty(), "{file}:\n{}", bad.join("\n"));
        // Opened as Refine opens it: the same points and pattern as the .pcr's.
        let o = drs::control::open(&path).unwrap();
        let (data, d0) = (o.data.as_ref().unwrap(), l.data.as_ref().unwrap());
        assert_eq!(data.x.len(), d0.x.len(), "{file}");
        assert!(data.x.iter().zip(&d0.x).all(|(a, b)| (a - b).abs() <= 1e-9 * b.abs()), "{file}: the axis");
        let (a, b) = (o.control.model.calculate(&data.x).ycalc, l.model.calculate(&d0.x).ycalc);
        assert!(a.iter().zip(&b).all(|(u, v)| (u - v).abs() <= 1e-9 * v.abs().max(1.0)), "{file}: the pattern");
        checked += 1;
    }
    assert_eq!(checked, 9);
}
