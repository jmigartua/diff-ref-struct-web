//! The user's notes on reflections (W8, L142): kept in `<control> notes.toml` beside the control file, carried by
//! `drs results` into the results file's pattern block and read back with it; `drs notes` lists them.

use drs::control::Control;
use drs::notes::Note;
use std::path::{Path, PathBuf};

fn root(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)
}

#[test]
fn notes_go_into_the_results_file() {
    let tmp = std::env::temp_dir().join(format!("drs_notes_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&tmp);
    std::fs::create_dir_all(&tmp).unwrap();
    let run = root("../manuals/fullprof/runs/Y2O3");
    std::fs::copy(run.join("y2o3.dat"), tmp.join("y2o3.dat")).unwrap();
    let text: String = std::fs::read(run.join("y2o3_structure.pcr")).unwrap().iter().map(|&b| b as char).collect();
    let c = Control::from_pcr(&drs::pcr::read(&text).unwrap(), &run, "y2o3_structure", &tmp).unwrap();
    let drc = tmp.join("y2o3.drc");
    std::fs::write(&drc, c.to_text()).unwrap();
    let notes = vec![
        Note { phase: Some(1), hkl: Some([2, 2, 2]), two_theta: None, text: "the strongest; a shoulder at low angle".into() },
        Note { phase: None, hkl: None, two_theta: Some(31.25), text: "no phase's: an impurity?\nto check".into() },
    ];
    let file = tmp.join(drs::formats::files::path("reflection_notes", &[("control", "y2o3")]).unwrap());
    assert_eq!(file, tmp.join("y2o3 notes.toml"));
    std::fs::write(&file, drs::notes::write(&notes, "y2o3.drc")).unwrap();
    let out = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(["results", drc.to_str().unwrap(), "--time", "2026-10-12 12.00.00"]).output().unwrap();
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let r = drs::results::read(&std::fs::read_to_string(tmp.join("Results/y2o3 2026-10-12 12.00.00 results.cif")).unwrap()).unwrap();
    assert_eq!(r.notes.len(), 2);
    assert_eq!(r.notes[0], notes[0]);
    assert_eq!((r.notes[1].two_theta, r.notes[1].text.as_str()), (Some(31.25), "no phase's: an impurity? to check"), "one line in the CIF");
    assert!(!r.reflections.is_empty(), "the reflections are still read");
    let out = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(["notes", drc.to_str().unwrap()]).output().unwrap();
    let t = String::from_utf8_lossy(&out.stdout);
    assert!(t.contains("# 2 notes") && t.contains("phase 1 (2 2 2)\tthe strongest") && t.contains("2θ 31.25°"), "{t}");
    // no notes file: no loop, nothing read
    std::fs::remove_file(&file).unwrap();
    let out = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(["results", drc.to_str().unwrap(), "--time", "2026-10-12 12.00.01"]).output().unwrap();
    assert!(out.status.success());
    let t = std::fs::read_to_string(tmp.join("Results/y2o3 2026-10-12 12.00.01 results.cif")).unwrap();
    assert!(!t.contains("_drs_note_text") && drs::results::read(&t).unwrap().notes.is_empty());
    std::fs::remove_dir_all(&tmp).unwrap();
}
