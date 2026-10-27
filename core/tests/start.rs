//! A refinement started from a data file and a structure, with no previous control file (`drs new`, L153).
//!
//! Y2O3's data with a CIF holding the atoms of FullProf's run (`tests/data/start/y2o3.cif`: its positions and B, which
//! that run does not refine) and the cell rounded to 10.60 Å, and the instrument as that run states it (Cu Kα1/Kα2
//! 1.5406/1.5444 at 0.5, Cthm 0.78, Wdt 40, the asymmetry below 30°, a pseudo-Voigt). The control file `start::start`
//! writes is checked as written (its data, reader, radiation, the phase, what is free, the background file, the
//! stated defaults), then refined: its own refinement first (the scale and the background), then the other
//! parameters FullProf's run refines, in three steps (zero and cell; U, V, W, η0, X; the asymmetry), each to
//! FullProf's convergence (Eps 0.2). The result is compared with FullProf's converged run (`y2o3_structure_conv.*`):
//! - with FullProf's own 62 background points in place of the automatic ones, the same file and steps reach FullProf's
//!   minimum: χ² within the tolerance of `fullprof_refinements` and every parameter within 1.2 FullProf esd's (one
//!   esd, and the 0.2 esd at which both refinements stop);
//! - with its automatic background (101 points, arPLS, all refined) the minimum moves: frozen below (L037), the
//!   difference being the background's, as the run with FullProf's points shows.

use drs::jacobian::Field;
use drs::lsq::{refine, Options};
use drs::powder::{BackgroundModel, Intensities, Model, Shape};
use drs::refinable::{address, Refinement};
use drs::start::{start, Instrument, Known};
use std::path::{Path, PathBuf};

mod common;
use common::{fullprof_field, fullprof_params, latin1};

fn root(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)
}

fn scratch(name: &str) -> PathBuf {
    let d = root("target/test-start").join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    std::fs::copy(root("../manuals/fullprof/runs/Y2O3/y2o3.dat"), d.join("y2o3.dat")).unwrap();
    d
}

/// The instrument of FullProf's Y2O3 run (y2o3_structure.pcr, line 8 and the profile's Npr = 5).
fn y2o3_instrument() -> Instrument {
    Instrument {
        wavelength: Some((1.5406, Some((1.5444, 0.5)))),
        monochromator_polarisation: Some(0.78),
        peak_window: Some(40.0),
        asymmetry_limit: Some(30.0),
        profile: Some(Shape::PseudoVoigt),
        ..Default::default()
    }
}

fn y2o3_cif() -> Known {
    Known::Cif { name: "Y2O3".into(), text: std::fs::read_to_string(root("tests/data/start/y2o3.cif")).unwrap() }
}

fn points(m: &Model) -> &[(f64, f64)] {
    match &m.background {
        BackgroundModel::Points(p) => p,
        _ => &[],
    }
}

/// The file's refinement, then the parameters of FullProf's run in three steps; the last result.
fn refine_in_steps(m: &mut Model, data: &drs::data::Pattern, first: &Refinement) -> drs::lsq::Refinement {
    let opts = Options { max_cycles: 100, eps: 0.2 };
    let mut free = first.parameters.clone();
    let set = first.build(m).unwrap();
    let mut r = refine(m, data, &set, opts, |_| {});
    for more in [vec![Field::Zero, Field::Cell(0, 0)], vec![Field::W(0), Field::U(0), Field::V(0), Field::Eta0(0), Field::X(0)], vec![Field::Asym(0, 0), Field::Asym(0, 1)]] {
        free.extend(more.into_iter().map(|f| vec![(f, 1.0)]));
        let set = Refinement { parameters: free.clone(), ..Default::default() }.build(m).unwrap();
        r = refine(m, data, &set, opts, |_| {});
        assert!(r.reason.converged(), "{}", r.stop);
    }
    r
}

/// FullProf's χ² over all non-excluded points in its last cycle (as `fullprof_refinements` reads it).
fn fullprof_chi2(out: &str) -> f64 {
    let lines: Vec<&str> = out.lines().collect();
    let start = lines.iter().rposition(|l| l.contains("RELIABILITY FACTORS WITH ALL NON-EXCLUDED POINTS")).unwrap();
    lines[start..].iter().find_map(|l| l.split_once("Conventional Rietveld").and_then(|(_, r)| r.split_once("Chi2:"))).map(|(_, r)| r.split_whitespace().nth(3).unwrap().parse().unwrap()).unwrap()
}

/// Each parameter FullProf's run refines (its background aside): (address, ours − FullProf's in FullProf's esd's), the
/// esd's translated as `fullprof_refinements` does (the X-ray scale ×2; P1 + 3P2 with σ₁ + 3σ₂).
fn against_fullprof(m: &Model) -> Vec<(String, f64)> {
    let fp = root("../manuals/fullprof/runs/Y2O3");
    let pcr = drs::pcr::read(&latin1(&fp.join("y2o3_structure_conv.pcr"))).unwrap();
    let conv = Model::from_pcr(&pcr, &fp, "y2o3_structure").unwrap().model;
    let labels: Vec<Vec<String>> = pcr.phases.iter().map(|ph| ph.atoms.iter().map(|a| a.label.clone()).collect()).collect();
    let fpp = fullprof_params(&latin1(&fp.join("y2o3_structure_conv.out")));
    let sig = |f: Field| fpp.iter().find(|(n, _)| fullprof_field(n, &labels) == Some(f)).map(|(_, v)| v.1);
    let mut out = Vec::new();
    for name in fpp.keys() {
        let Some(f) = fullprof_field(name, &labels) else { continue };
        let e = match f {
            Field::Background(_) | Field::Asym(_, 1) => continue,
            Field::Scale(_) => 2.0 * sig(f).unwrap(),
            Field::Asym(q, 0) => sig(f).unwrap() + 3.0 * sig(Field::Asym(q, 1)).unwrap(),
            _ => sig(f).unwrap(),
        };
        out.push((address(f, m), (m.get(f) - conv.get(f)) / e));
    }
    out
}

/// With the automatic background, the parameters more than one FullProf esd from FullProf's values (address, ours −
/// FullProf's in esd's, signed), frozen: the background's 101 refined points against FullProf's 62 take up part of
/// the peaks' Lorentzian tails (η0 and X) and of the intensity (the scale); with FullProf's points every one is within
/// 1.04 esd. The test fails if one is added, grows, or changes sign.
const FROZEN: &[(&str, f64)] = &[
    ("pattern1.zero_shift", 1.51),
    ("phase1.pv_eta0", 7.49),
    ("phase1.caglioti_w", -2.21),
    ("phase1.pv_eta_slope", -3.20),
    ("phase1.scale_factor", 8.40),
];

/// χ² with the automatic background against FullProf's, (ours − FullProf's)/FullProf's, frozen: lower, the 101
/// refined points following the data more closely than FullProf's 62.
const CHI2_FROZEN: f64 = -0.02643;

#[test]
fn start_from_data_and_a_structure() {
    let d = scratch("y2o3");
    let out = d.join("y2o3.drc");
    let s = start(&d.join("y2o3.dat"), None, &[y2o3_cif()], &y2o3_instrument(), &out).unwrap();
    s.write(&out).unwrap();
    assert!(s.write(&out).unwrap_err().contains("exists"), "nothing is overwritten");

    // The file as written, opened as the app and the command line open it.
    let o = drs::control::open(&out).unwrap();
    let c = &o.control;
    let m0 = &c.model;
    let data = o.data.as_ref().unwrap();
    let src = c.data.as_ref().unwrap();
    assert_eq!((src.file.as_str(), src.reader.as_deref()), ("y2o3.dat", Some("free_steps")));
    assert_eq!((m0.lambda, m0.ratio, m0.cthm, m0.wdt, m0.asym_limit, m0.neutrons), ([1.5406, 1.5444], 0.5, 0.78, 40.0, 30.0, false));
    assert_eq!(m0.convention, drs::powder::Convention::DEFAULT);
    assert_eq!(m0.phases.len(), 1);
    let ph = &m0.phases[0];
    assert_eq!((ph.group.hall.as_str(), ph.cell.a, ph.shape, ph.extract), ("-I 2b 2c 3", 10.60, Shape::PseudoVoigt, false));
    let Intensities::Structure(st, _) = &ph.intensities else { panic!("a structure") };
    // FullProf's Occ, the chemical occupancy times m/M: 8/48, 24/48, 48/48.
    assert_eq!(st.atoms.iter().map(|a| a.label.as_str()).collect::<Vec<_>>(), ["Y1", "Y2", "O1"]);
    for (a, want) in st.atoms.iter().zip([8.0 / 48.0, 0.5, 1.0]) {
        assert!((a.occ - want).abs() < 1e-12, "{} against {want}", a.occ);
    }
    assert_eq!((ph.u, ph.v, ph.eta0, ph.x, ph.asym, ph.pref[0]), (0.0, 0.0, 0.0, 0.0, [0.0; 4], 1.0));
    // W from the data's peaks.
    let (h, n) = drs::start::median_fwhm(data).unwrap();
    assert_eq!((ph.w, s.fwhm, s.peaks), (h * h, h, n));
    println!("{n} peaks, median FWHM {h:.4}°: W = {:.5} deg²", ph.w);
    // The background: arPLS's points, in their own file, read back exactly.
    let auto = drs::data_ops::auto_background(data, None, drs::data_ops::DEFAULT_MAX_ITER).unwrap();
    assert_eq!(c.background_file.as_deref(), Some("y2o3 background.drb"));
    assert!(d.join("y2o3 background.drb").exists());
    assert_eq!(points(m0), auto.points.as_slice());
    assert_eq!(points(m0).len(), 101);
    // Free: the scale and every background point, nothing else.
    let set = c.refinement.build(m0).unwrap();
    let mut names: Vec<String> = set.free.iter().map(|p| p.name.clone()).collect();
    names.sort();
    let mut want: Vec<String> = (0..101).map(|k| format!("pattern1.background_points.{k}")).chain(["phase1.scale_factor".to_string()]).collect();
    want.sort();
    assert_eq!(names, want);
    // The scale is the least-squares one: χ² rises either side of it, the background held.
    let chi2 = |s: f64| {
        let mut m = m0.clone();
        m.phases[0].scale *= s;
        drs::lsq::statistics(&m, data, &m.calculate(&data.x).ycalc, 0).chi2
    };
    assert!(chi2(1.01) > chi2(1.0) && chi2(0.99) > chi2(1.0));
    // Every default stated, at the head of the file; what was given is not among them.
    let text = std::fs::read_to_string(&out).unwrap();
    for st in &s.stated {
        assert!(text.contains(&format!("#   {} = {}: {}", st.key, st.value, st.source)), "{} not stated", st.key);
    }
    let keys: Vec<&str> = s.stated.iter().map(|x| x.key.as_str()).collect();
    for k in ["geometry", "weighting", "[convention]", "background", "caglioti_w", "phase 1 scale_factor", "cycles, convergence"] {
        assert!(keys.contains(&k), "{k}");
    }
    for k in ["monochromator_polarisation", "peak_window", "asymmetry_limit", "profile", "wavelength_ratio"] {
        assert!(!keys.contains(&k), "{k} was given");
    }
    // The checks find no error in it.
    let diags = drs::check::check_path(&out, &drs::check::Options::default());
    assert!(!drs::check::has_errors(&diags), "{:?}", diags.iter().map(|x| &x.what).collect::<Vec<_>>());

    // Refined in steps from the file, with FullProf's background points: FullProf's minimum.
    let fp = root("../manuals/fullprof/runs/Y2O3");
    let fp_chi2 = fullprof_chi2(&latin1(&fp.join("y2o3_structure_conv.out")));
    let pin = drs::pcr::read(&latin1(&fp.join("y2o3_structure.pcr.in"))).unwrap();
    let fp_start = Model::from_pcr(&pin, &fp, "y2o3_structure").unwrap().model;
    let mut m = m0.clone();
    m.background = fp_start.background.clone();
    let mut first = c.refinement.clone();
    first.parameters = (0..points(&m).len()).map(|k| vec![(Field::Background(k), 1.0)]).chain([vec![(Field::Scale(0), 1.0)]]).collect();
    let r = refine_in_steps(&mut m, data, &first);
    let tol = 5e-4 + 1.5e-3 * fp_chi2;
    println!("with FullProf's background points: χ² {:.4}, FullProf {fp_chi2} (tolerance {tol:.4})", r.stats.chi2);
    assert!((r.stats.chi2 - fp_chi2).abs() <= tol);
    for (name, dev) in against_fullprof(&m) {
        println!("    {name:28} {dev:+.2} esd");
        assert!(dev.abs() <= 1.2, "{name}: {dev:+.2} esd from FullProf's value");
    }

    // Refined in steps from the file as written: its automatic background moves the minimum (frozen).
    let mut m = m0.clone();
    let r = refine_in_steps(&mut m, data, &c.refinement);
    let rel = (r.stats.chi2 - fp_chi2) / fp_chi2;
    println!("with the automatic background: χ² {:.4} ({rel:+.6} of FullProf's)", r.stats.chi2);
    assert!(rel.signum() == CHI2_FROZEN.signum() && rel.abs() <= CHI2_FROZEN.abs(), "χ² {:.4}: {rel:+.4} of FullProf's, frozen {CHI2_FROZEN}", r.stats.chi2);
    let mut beyond = Vec::new();
    for (name, dev) in against_fullprof(&m) {
        println!("    {name:28} {dev:+.4} esd");
        if dev.abs() > 1.0 {
            beyond.push(name.clone());
            match FROZEN.iter().find(|f| f.0 == name) {
                Some(f) => assert!(dev.signum() == f.1.signum() && dev.abs() <= f.1.abs(), "{name}: {dev:+.2} esd, frozen {}", f.1),
                None => panic!("{name}: {dev:+.2} esd from FullProf's value, not frozen"),
            }
        }
    }
    assert_eq!(beyond.len(), FROZEN.len(), "fewer parameters beyond one esd than frozen: freeze the smaller count");
}

/// What is not given is a default, stated: with only an anode and a CIF.
#[test]
fn the_defaults_are_stated() {
    let d = scratch("defaults");
    let out = d.join("y2o3.drc");
    let s = start(&d.join("y2o3.dat"), None, &[y2o3_cif()], &Instrument { anode: Some("cu".into()), ..Default::default() }, &out).unwrap();
    let m = &s.control.model;
    let (_, k1, k2, _) = drs::scattering_table::ANODES[2];
    assert_eq!((m.lambda, m.ratio), ([k1, k2], drs::start::KALPHA_RATIO));
    // No monochromator: Cthm = cos² 0 = 1 (FullProf's manual, eq. 3.26).
    assert_eq!((m.cthm, m.wdt, m.asym_limit), (1.0, drs::start::DEFAULT_PEAK_WINDOW, drs::start::DEFAULT_ASYMMETRY_LIMIT));
    assert_eq!(m.phases[0].shape, Shape::ThompsonCoxHastings);
    let stated: Vec<(&str, &str)> = s.stated.iter().map(|x| (x.key.as_str(), x.value.as_str())).collect();
    for kv in [("wavelength_ratio", "0.5"), ("monochromator_polarisation", "1"), ("peak_window", "8"), ("asymmetry_limit", "170"), ("profile", "thompson-cox-hastings")] {
        assert!(stated.contains(&kv), "{kv:?} not stated");
    }
    // A data file that does not say its wavelength, and none given: refused.
    let e = start(&d.join("y2o3.dat"), None, &[y2o3_cif()], &Instrument::default(), &out).unwrap_err();
    assert!(e.contains("wavelength"), "{e}");
    assert!(start(&d.join("y2o3.dat"), None, &[], &Instrument { anode: Some("Cu".into()), ..Default::default() }, &out).is_err());
}

/// A formula, a space group and a cell: a Le Bail phase, its reflections those the group allows to 5° past the
/// data, its scale fixed; the file refines. A Le Bail fit runs every cycle asked for (`lsq::Stop::Cycles`), as
/// FullProf's profile matching does.
#[test]
fn a_cell_without_atoms_is_a_le_bail_phase() {
    let d = scratch("lebail");
    let out = d.join("y2o3.drc");
    let known = [Known::Cell { formula: "Y2O3".into(), group: "I a -3".into(), cell: vec![10.60] }];
    let s = start(&d.join("y2o3.dat"), None, &known, &y2o3_instrument(), &out).unwrap();
    s.write(&out).unwrap();
    let o = drs::control::open(&out).unwrap();
    let ph = &o.control.model.phases[0];
    assert!(ph.extract);
    let Intensities::List(list) = &ph.intensities else { panic!("a list") };
    let g = drs::spacegroup::SpaceGroup::from_symbol("I a -3").unwrap().0;
    let dmin = 1.5406 / (2.0 * (150.9f64 / 2.0).to_radians().sin());
    assert_eq!(list.len(), drs::reflections::generate(&g, &ph.cell, dmin, 1e9, None).len());
    let set = o.control.refinement.build(&o.control.model).unwrap();
    assert!(set.free.iter().all(|p| p.name.starts_with("pattern1.background_points.")), "only the background is free");
    // Refined as the Rietveld file is: the background, then zero and cell, then the widths.
    let data = o.data.as_ref().unwrap();
    let mut m = o.control.model.clone();
    let before = drs::lsq::statistics(&m, data, &m.calculate(&data.x).ycalc, 0).chi2;
    let mut free = o.control.refinement.parameters.clone();
    let mut last = None;
    for more in [vec![], vec![Field::Zero, Field::Cell(0, 0)], vec![Field::W(0), Field::U(0), Field::V(0), Field::Eta0(0), Field::X(0)]] {
        free.extend(more.into_iter().map(|f| vec![(f, 1.0)]));
        let set = Refinement { parameters: free.clone(), ..Default::default() }.build(&m).unwrap();
        last = Some(refine(&mut m, data, &set, Options { max_cycles: 20, eps: 0.2 }, |_| {}));
    }
    let r = last.unwrap();
    println!("Le Bail: χ² {before:.3} → {:.4} ({}), a = {:.6} Å", r.stats.chi2, r.stop, m.phases[0].cell.a);
    assert!(matches!(r.reason, drs::lsq::Stop::Cycles { n: 20, le_bail: true }) && r.stats.chi2 < before);
    // A cubic cell is one value; a wrong count is refused, as an unknown element.
    assert!(drs::start::cell_from(&g, &[10.6, 10.6]).unwrap_err().contains("cubic"));
    assert!(drs::start::formula_elements("Y2Q3").is_err());
    assert_eq!(drs::start::formula_elements("Y2O3").unwrap(), vec![("Y".to_string(), 2.0), ("O".to_string(), 3.0)]);
}

/// `drs new` writes what `start::start` gives.
#[test]
fn the_command_line() {
    let d = scratch("cli");
    let cif = root("tests/data/start/y2o3.cif");
    let run = |args: &[&str]| std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(args).output().unwrap();
    let dat = d.join("y2o3.dat");
    let o = run(&["new", dat.to_str().unwrap(), "--cif", cif.to_str().unwrap(), "--wavelength", "1.5406,1.5444,0.5", "--polarisation", "0.78", "--window", "40", "--asymmetry-limit", "30", "--profile", "pseudo-voigt"]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let written = std::fs::read_to_string(d.join("y2o3.drc")).unwrap();
    let s = start(&dat, None, &[Known::Cif { name: "y2o3".into(), text: std::fs::read_to_string(&cif).unwrap() }], &y2o3_instrument(), &d.join("y2o3.drc")).unwrap();
    assert_eq!(written, s.text);
    assert_eq!(std::fs::read_to_string(d.join("y2o3 background.drb")).unwrap(), s.background_file.1);
    let again = run(&["new", dat.to_str().unwrap(), "--cif", cif.to_str().unwrap(), "--anode", "Cu"]);
    assert!(!again.status.success() && String::from_utf8_lossy(&again.stderr).contains("exists"));
    let lb = run(&["new", dat.to_str().unwrap(), "--phase", "Y2O3", "I a -3", "10.6", "--anode", "Cu", "--out", d.join("lb.drc").to_str().unwrap()]);
    assert!(lb.status.success(), "{}", String::from_utf8_lossy(&lb.stderr));
    assert!(d.join("lb background.drb").exists());
}
