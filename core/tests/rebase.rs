//! A control file pointed at new data (R25, APP_WORKFLOWS W3; L078).
//! - `onto_its_own_data`: Y2O3 rebased onto its own data keeps its model, pattern and χ², and reports only that the
//!   data file does not state its wavelength.
//! - `onto_other_data`: Y2O3 onto LaB6's data (another range and step; background points of other angles) and
//!   Co-YBCO onto Y2O3's (an excluded region outside the data, a Chebyshev range of other angles): each diagnostic
//!   reported, with what, why and a suggestion; the file written checks with no error.
//! - `what_does_not_fit`: on a scan made for the purpose, another wavelength (with the line to write), another axis,
//!   an excluded region over all the data, the old data's operations, and background points of other counts.

use drs::check::Severity;
use drs::control::Control;
use drs::formats::{Axis, Scan, Wavelength};
use drs::rebase::{rebase, rebase_path};
use std::path::{Path, PathBuf};

fn root(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)
}

fn tmp(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("drs_rebase_{name}_{}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn codes(r: &drs::rebase::Rebased) -> Vec<&'static str> {
    r.diagnostics.iter().map(|d| d.code).collect()
}

#[test]
fn onto_its_own_data() {
    let run = root("../manuals/fullprof/runs/Y2O3");
    let d = tmp("self");
    let r = rebase_path(&run.join("y2o3_structure.pcr"), &run.join("y2o3.dat"), None, &d.join("y.drc")).unwrap();
    assert_eq!(codes(&r), ["rebase-wavelength"], "{:?}", r.diagnostics);
    assert_eq!(r.diagnostics[0].severity, Severity::Hint);
    let o = drs::control::open(&d.join("y.drc")).unwrap();
    let p = drs::pcr::read(&std::fs::read(run.join("y2o3_structure.pcr")).unwrap().iter().map(|&b| b as char).collect::<String>()).unwrap();
    let l = drs::powder::Model::from_pcr(&p, &run, "y2o3_structure").unwrap();
    let data = l.data.unwrap();
    assert_eq!(o.control.model.calculate(&data.x).ycalc, l.model.calculate(&data.x).ycalc);
    assert_eq!(o.data.unwrap(), data);
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn onto_other_data() {
    let runs = root("../manuals/fullprof/runs");
    let d = tmp("other");
    let r = rebase_path(&runs.join("Y2O3/y2o3_structure.pcr"), &runs.join("LaB6/LaB6.dat"), Some("xys"), &d.join("y.drc")).unwrap();
    let c = codes(&r);
    for want in ["rebase-range", "rebase-step", "rebase-background"] {
        assert!(c.contains(&want), "Y2O3 onto LaB6: {c:?}");
    }
    assert!(!drs::check::has_errors(&r.diagnostics), "{c:?}");
    let r = rebase_path(&runs.join("Co-YBCO/Cocus-asym.pcr"), &runs.join("Y2O3/y2o3.dat"), None, &d.join("c.drc")).unwrap();
    let c = codes(&r);
    assert!(c.contains(&"rebase-excluded") && c.contains(&"rebase-background"), "Co-YBCO onto Y2O3: {c:?}");
    let cheb = r.diagnostics.iter().find(|x| x.code == "rebase-background").unwrap();
    assert_eq!(cheb.replacement.as_deref(), Some("from = 15\nto = 145.9"));
    for x in r.diagnostics.iter().filter(|x| x.code.starts_with("rebase-")) {
        assert!(!x.what.is_empty() && !x.why.is_empty() && x.suggestion.is_some(), "{x:?}");
    }
    println!("Co-YBCO onto Y2O3: {c:?}");
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn what_does_not_fit() {
    let c = Control::from_text(&std::fs::read_to_string(root("tests/data/broken/base.drc")).unwrap()).unwrap();
    let x: Vec<f64> = (0..1000).map(|i| 20.0 + 0.1 * i as f64).collect();
    let pattern = drs::data::Pattern { y: vec![170.0; x.len()], sigma: vec![13.0; x.len()], x };
    // Another wavelength and another axis.
    let mut s = Scan::new("new", Axis::D, pattern.clone());
    s.wavelength = Some(Wavelength::single(0.7093));
    let r = rebase(&c, None, &s, "new.xye".into(), Some("xye".into()));
    let w = r.diagnostics.iter().find(|d| d.code == "rebase-wavelength").unwrap();
    assert_eq!((w.severity, w.replacement.as_deref()), (Severity::Warning, Some("wavelength_1 = 0.7093")));
    assert!(r.diagnostics.iter().any(|d| d.code == "rebase-axis" && d.severity == Severity::Error));
    assert_eq!(r.control.data.as_ref().unwrap().file, "new.xye");
    assert!(r.control.fullprof.is_empty(), "the old .pcr's range is dropped");
    // An excluded region over all the data; the old data's operations; points of ten times fewer counts.
    let mut c2 = c.clone();
    c2.model.excluded = vec![(0.0, 180.0)];
    c2.data.as_mut().unwrap().operations = vec![drs::data_ops::Operation::Mask { from: 30.0, to: 31.0 }];
    let s2 = Scan::new("new", Axis::TwoTheta, pattern);
    let r = rebase(&c2, None, &s2, "new.xye".into(), None);
    let get = |code: &str| r.diagnostics.iter().filter(|d| d.code == code).collect::<Vec<_>>();
    assert!(get("rebase-excluded").iter().any(|d| d.severity == Severity::Error));
    assert_eq!(get("rebase-operations").len(), 1);
    assert!(get("rebase-background").iter().any(|d| d.what.contains("times the background points' counts")), "{:?}", r.diagnostics);
    assert!(r.control.data.as_ref().unwrap().operations.is_empty());
}
