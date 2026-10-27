//! The core's checks (R24, L076).
//! - `broken_files`: every file of `tests/data/broken/` (written by `make_broken.py`, each a run with one mistake)
//!   raises exactly one diagnostic, the one its name gives; the bases they come from raise none.
//! - `good_runs`: the eight runs FullProf reproduces from the core's export (L071), as `.pcr` and as control files,
//!   raise no error.
//! - `pcr_data`: a `.pcr`'s data taken as FullProf takes them (L079): the file given, else CODFIL.dat next to it;
//!   a file found only by the comment is used with a hint, none an error; a `.bac` in pairs warned (L082).
//! - `predicted_correlations`: on each of those runs, the correlations the checks predict at the starting values
//!   against those the engine reports after refining from them.

use drs::check::{self, Options, Severity};
use drs::powder::Model;
use drs::refinable::ParameterSet;
use std::path::{Path, PathBuf};

fn root(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)
}

/// The eight runs FullProf gives its own numbers on from the core's export (L071): folder, .pcr stem.
const GOOD: &[(&str, &str)] = &[
    ("Co-YBCO", "Cocus-asym"),
    ("Egypt", "Egypt_qa_Rietveld"),
    ("MgAl2O4-s", "MgAl2O4s_Rietveld"),
    ("Y2O3", "y2o3_structure"),
    ("YVO3", "YVO3a"),
    ("Yb2ZrWO8", "ybfluo_Riet_IRF"),
    ("Yb2ZrWO8-PO-March", "ybfluo_Riet_IRF"),
    ("rutile-anatase", "rutana_tch"),
];

fn summary(d: &[check::Diagnostic], file: &str) -> String {
    d.iter().map(|x| x.to_text(file)).collect::<Vec<_>>().join("\n")
}

#[test]
fn broken_files() {
    let dir = root("tests/data/broken");
    let opts = Options::default();
    // The .pcr cases are checked with the Y2O3 data given, as FullProf is run on them (fp2k <case> y2o3.dat,
    // broken_runs.py; L079).
    let y2o3 = root("../manuals/fullprof/runs/Y2O3/y2o3.dat");
    let mut names: Vec<PathBuf> = std::fs::read_dir(&dir).unwrap().map(|e| e.unwrap().path()).collect();
    names.sort();
    let mut checked = Vec::new();
    for path in names {
        let name = path.file_name().unwrap().to_string_lossy().to_string();
        let stem = name.split('.').next().unwrap().to_string();
        let ext = path.extension().map(|e| e.to_string_lossy().to_string()).unwrap_or_default();
        let diags = match ext.as_str() {
            "drc" => check::check_path(&path, &opts),
            "pcr" => {
                // The .pcr exports of the .drc cases are for FullProf (L076); the .pcr cases themselves are checked, and
                // the base export too.
                if dir.join(format!("{stem}.drc")).exists() && stem != "base" {
                    continue;
                }
                check::check_path_with_data(&path, Some(&y2o3), &opts)
            }
            _ => continue,
        };
        if stem.starts_with("base") {
            assert!(diags.is_empty(), "{name}: the base raises\n{}", summary(&diags, &name));
            continue;
        }
        let code = stem.split("--").next().unwrap();
        assert_eq!(diags.len(), 1, "{name}: expected one diagnostic [{code}], got\n{}", summary(&diags, &name));
        assert_eq!(diags[0].code, code, "{name}:\n{}", summary(&diags, &name));
        // Every diagnostic says what, why and what to do.
        assert!(!diags[0].what.is_empty() && !diags[0].why.is_empty() && diags[0].suggestion.is_some(), "{name}: {:?}", diags[0]);
        checked.push(format!("{name}: {}[{}]", diags[0].severity.name(), diags[0].code));
    }
    println!("{} broken files, each its own diagnostic:\n  {}", checked.len(), checked.join("\n  "));
    assert!(checked.len() >= 18, "{} broken files", checked.len());
}

#[test]
fn good_runs() {
    let opts = Options::default();
    let tmp = std::env::temp_dir().join(format!("drs_good_{}", std::process::id()));
    std::fs::create_dir_all(&tmp).unwrap();
    for (run, code) in GOOD {
        let d = root(&format!("../manuals/fullprof/runs/{run}"));
        let pcr = d.join(format!("{code}.pcr"));
        let diags = check::check_path(&pcr, &opts);
        let errors: Vec<_> = diags.iter().filter(|x| x.severity == Severity::Error).collect();
        assert!(errors.is_empty(), "{run}: {}", summary(&diags, code));
        // The same run as a control file.
        let text: String = std::fs::read(&pcr).unwrap().iter().map(|&b| b as char).collect();
        let mut p = drs::pcr::read(&text).unwrap();
        p.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
        let c = drs::control::Control::from_pcr(&p, &d, code, &tmp).unwrap();
        let drc = tmp.join(format!("{run}.drc"));
        std::fs::write(&drc, c.to_text()).unwrap();
        let d2 = check::check_path(&drc, &opts);
        assert!(!check::has_errors(&d2), "{run}.drc: {}", summary(&d2, run));
        let codes = |v: &[check::Diagnostic]| v.iter().map(|x| x.code).collect::<Vec<_>>();
        println!("{run:18} .pcr: {:?}  .drc: {:?}", codes(&diags), codes(&d2));
    }
    let _ = std::fs::remove_dir_all(&tmp);
}

/// Run, the largest |predicted − after refining| measured (rounded up in the fifth decimal), the pairs above 0.9 in one and not the other (L076).
const FROZEN: &[(&str, f64, usize)] = &[
    // Re-measured in core 0.24.0 after the fix of L108 (FullProf's S_L moves the default model's H/L too, which the
    // import had left out): 0.03323, from 0.03169 with the old set (checked with the fix disabled: 0.031683).
    ("Co-YBCO", 0.03324, 1),
    ("Egypt", 0.00042, 0),
    ("MgAl2O4-s", 0.00052, 0),
    ("Y2O3", 0.00339, 0),
    ("YVO3", 0.00979, 0),
    ("Yb2ZrWO8", 0.00031, 0),
    ("Yb2ZrWO8-PO-March", 0.00514, 0),
    ("rutile-anatase", 0.00611, 0),
];

#[test]
fn predicted_correlations() {
    for (run, code) in GOOD {
        let d = root(&format!("../manuals/fullprof/runs/{run}"));
        let text: String = std::fs::read(d.join(format!("{code}.pcr"))).unwrap().iter().map(|&b| b as char).collect();
        let mut p = drs::pcr::read(&text).unwrap();
        p.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
        let l = Model::from_pcr(&p, &d, code).unwrap();
        let data = l.data.as_ref().unwrap();
        let set = ParameterSet::from_pcr(&p, &l.model).unwrap();
        let predicted = check::predicted_correlations(&l.model, data, &set);
        let mut m = l.model.clone();
        let r = drs::lsq::refine(&mut m, data, &set, drs::lsq::Options { max_cycles: 30, eps: p.eps }, |_| {});
        let n = set.free.len();
        let (mut worst, mut at) = (0f64, (0, 0));
        let (mut both, mut only_p, mut only_e) = (0, 0, 0);
        for (i, (pr, er)) in predicted.iter().zip(&r.correlations).enumerate() {
            for j in i + 1..n {
                let (a, b) = (pr[j], er[j]);
                if (a - b).abs() > worst {
                    worst = (a - b).abs();
                    at = (i, j);
                }
                match (a.abs() >= 0.9, b.abs() >= 0.9) {
                    (true, true) => both += 1,
                    (true, false) => only_p += 1,
                    (false, true) => only_e += 1,
                    _ => {}
                }
            }
        }
        println!(
            "{run:18} {n:3} parameters, {} cycles: largest |Δr| {worst:.6} ({} – {}); pairs ≥ 0.9: {both} both, {only_p} predicted only, {only_e} after refining only",
            r.cycles,
            set.free.get(at.0).map_or("", |q| q.name.as_str()),
            set.free.get(at.1).map_or("", |q| q.name.as_str())
        );
        // Frozen (L076): measured on these runs from their refined .pcr, they must not grow.
        // Frozen (L076), measured from each run's refined .pcr: the largest difference and the pairs that cross 0.9
        // must not grow. Co-YBCO's cell moves during its 4 cycles (a – c 0.032 apart, one pair crossing 0.9).
        let (_, bound, crossing) = FROZEN.iter().find(|f| f.0 == *run).unwrap();
        assert!(worst <= *bound, "{run}: predicted and refined correlations differ by {worst} (frozen at {bound})");
        assert!(only_p + only_e <= *crossing, "{run}: {} pairs cross 0.9 (frozen at {crossing})", only_p + only_e);
    }
}

#[test]
fn pcr_data() {
    let opts = Options { run: false, ..Options::default() };
    let base = root("tests/data/broken/base.pcr");
    let y2o3 = root("../manuals/fullprof/runs/Y2O3/y2o3.dat");
    let codes = |d: &[check::Diagnostic]| d.iter().map(|x| x.code).collect::<Vec<_>>();
    // Given, as on FullProf's command line: no diagnostic.
    assert_eq!(codes(&check::check_path_with_data(&base, Some(&y2o3), &opts)), Vec::<&str>::new());
    // Not given and no base.dat next to it: the comment's file is found in the folders around, with a hint.
    let d = check::check_path(&base, &opts);
    let hint = d.iter().find(|x| x.code == "data-elsewhere");
    println!("{}", summary(&d, "base.pcr"));
    // CODFIL.dat next to the .pcr: read without a word.
    // Three folders down, so that the folders around it, where the core looks, hold nothing else.
    let top = std::env::temp_dir().join(format!("drs_pcr_data_{}", std::process::id()));
    let tmp = top.join("a/b/c");
    std::fs::create_dir_all(&tmp).unwrap();
    std::fs::copy(&base, tmp.join("base.pcr")).unwrap();
    let alone = check::check_path(&tmp.join("base.pcr"), &opts);
    std::fs::copy(&y2o3, tmp.join("base.dat")).unwrap();
    assert_eq!(codes(&check::check_path(&tmp.join("base.pcr"), &opts)), Vec::<&str>::new());
    let _ = std::fs::remove_dir_all(&top);
    match hint {
        Some(h) => assert_eq!(h.severity, Severity::Hint),
        None => assert!(d.iter().any(|x| x.code == "data-missing")),
    }
    // Alone in a folder: no data, an error that says how FullProf finds them.
    let missing = alone.iter().find(|x| x.code == "data-missing").expect("data-missing");
    assert_eq!(missing.severity, Severity::Error);
    assert!(missing.why.contains("fp2k base <data>"), "{}", missing.why);
    // A .bac in pairs, which FullProf reads as its first format (L082): warned; the first format is not.
    let b = root("../manuals/fullprof/runs/background");
    let dat = b.join("y2o3.dat");
    assert_eq!(codes(&check::check_path_with_data(&b.join("y2o3_bacpairs.pcr"), Some(&dat), &opts)), ["bac-format"]);
    assert_eq!(codes(&check::check_path_with_data(&b.join("y2o3_bacfirst.pcr"), Some(&dat), &opts)), Vec::<&str>::new());
}
