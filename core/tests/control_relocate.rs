//! A control file moved to another folder, its paths rewritten (`control::relocate`, the app's request 45, L337).
//!
//! - Y₂O₃'s control file in `project/a/`, its data in `project/data/`, its background points in a file beside it and a
//!   data operation adding another file, a comment and a value written in the person's way: moved to
//!   `project/b/c/`, every path reaches the same file (`../../data/…`, `../../a/…`), the rest of the text is the same
//!   line for line, and opened there it gives the same data and the same model as before.
//! - The mode case A with its sign reference (`mode_reference = { file = … }`), moved one folder down.
//! - scor_rxn's joint control file, its patterns' control files named in it, moved to another folder and opened there:
//!   the same patterns.
//! - Refused, the file named: a path to a file that is not there, and a text that does not read.

use drs::control::relocate::paths;
use drs::control::{relocate, Control};
use std::path::{Path, PathBuf};

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap().iter().map(|&b| b as char).collect()
}

fn runs() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs")
}

fn fresh(name: &str) -> PathBuf {
    let d = Path::new(env!("CARGO_MANIFEST_DIR")).join("target").join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

#[test]
fn moved_with_its_files() {
    let root = fresh("test-relocate");
    for sub in ["a", "data", "b/c"] {
        std::fs::create_dir_all(root.join(sub)).unwrap();
    }
    let d = runs().join("Y2O3");
    std::fs::copy(d.join("y2o3.dat"), root.join("data/y2o3.dat")).unwrap();
    std::fs::copy(d.join("y2o3.dat"), root.join("data/empty holder.dat")).unwrap();
    let a = root.join("a");
    let p = drs::pcr::read(&latin1(&d.join("y2o3_structure.pcr"))).unwrap();
    let mut c = Control::from_pcr(&p, &d, "y2o3_structure", &a).unwrap();
    let src = c.data.as_mut().unwrap();
    src.file = "../data/y2o3.dat".into();
    src.operations.push(drs::data_ops::Operation::Subtract { file: "../data/empty holder.dat".into(), reader: src.reader.clone(), scale: 0.0 });
    c.background_file = Some("y2o3 background.drb".into());
    let (name, points) = c.background_file_text(None).unwrap().unwrap();
    std::fs::write(a.join(&name), points).unwrap();
    // The person's own text: a comment, and a value written in their way.
    let canonical = c.to_text();
    let first = canonical.lines().find(|l| l.starts_with("peak_window = ")).unwrap().to_string();
    let text = format!("# Y2O3, my notes\n{}", canonical.replacen(&first, "peak_window = 8.00   # as I like it", 1));
    std::fs::write(a.join("y2o3.drc"), &text).unwrap();

    let found: Vec<(String, String, String)> = paths(&drs::control::syntax::parse(&text)).into_iter().map(|p| (p.table, p.key, p.path)).collect();
    assert_eq!(found, [("pattern".to_string(), "data".to_string(), "../data/y2o3.dat".to_string()), ("pattern.background".into(), "file".into(), "y2o3 background.drb".into()), ("pattern.operation".into(), "file".into(), "../data/empty holder.dat".into())]);

    let to = root.join("b/c");
    let moved = relocate(&text, &a, &to).unwrap();
    std::fs::write(to.join("y2o3.drc"), &moved).unwrap();
    let new: Vec<String> = paths(&drs::control::syntax::parse(&moved)).into_iter().map(|p| p.path).collect();
    assert_eq!(new, ["../../data/y2o3.dat", "../../a/y2o3 background.drb", "../../data/empty holder.dat"]);
    // The rest of the text as it was, line for line.
    let (l0, l1): (Vec<&str>, Vec<&str>) = (text.lines().collect(), moved.lines().collect());
    assert_eq!(l0.len(), l1.len());
    let changed: Vec<usize> = (0..l0.len()).filter(|&i| l0[i] != l1[i]).collect();
    assert_eq!(changed.len(), 3, "only the three paths' lines change");
    assert!(moved.starts_with("# Y2O3, my notes\n") && moved.contains("peak_window = 8.00   # as I like it"));
    // Opened in its new folder: the same data and model.
    let (o0, o1) = (drs::control::open(&a.join("y2o3.drc")).unwrap(), drs::control::open(&to.join("y2o3.drc")).unwrap());
    let (d0, d1) = (o0.data.unwrap(), o1.data.unwrap());
    assert!(d0.x == d1.x && d0.y == d1.y && d0.sigma == d1.sigma, "the same data");
    assert_eq!(o0.control.model.to_text(), o1.control.model.to_text(), "the same model, its background points read");
    // Moving it back gives the text it had.
    assert_eq!(relocate(&moved, &to, &a).unwrap(), text);

    // A path to nothing, and a text that does not read: refused.
    let gone = text.replace("../data/y2o3.dat", "../data/none.dat");
    let e = relocate(&gone, &a, &to).unwrap_err();
    assert!(e.contains("none.dat") && e.contains("[pattern] data"), "{e}");
    assert!(relocate("[pattern\n", &a, &to).is_err());
    assert!(relocate(&text, &a, &root.join("no such folder")).unwrap_err().contains("is not a folder"));
    // `drs move`: the file written in the folder with its paths rewritten, the old one gone.
    let cli = root.join("b/cli");
    std::fs::create_dir_all(&cli).unwrap();
    std::fs::write(a.join("again.drc"), &text).unwrap();
    let st = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(["move", a.join("again.drc").to_str().unwrap(), cli.to_str().unwrap()]).status().unwrap();
    assert!(st.success() && !a.join("again.drc").exists());
    assert_eq!(std::fs::read_to_string(cli.join("again.drc")).unwrap(), relocate(&text, &a, &cli).unwrap());
    println!("moved: {new:?}; refused: {e}");
}

#[test]
fn the_sign_reference_moved() {
    let d = runs().join("modes/A_modes");
    let root = fresh("test-relocate-modes");
    let to = root.join("deeper");
    std::fs::create_dir_all(&to).unwrap();
    let pcr = drs::pcr::read(&latin1(&d.join("a_modes.pcr"))).unwrap();
    let c = Control::from_pcr(&pcr, &d, "a_modes", &d).unwrap();
    let text = c.to_text().replacen("[[phase]]\n", "[[phase]]\nmode_reference = { file = \"a_modes.sum\" }\n", 1);
    let moved = relocate(&text, &d, &to).unwrap();
    let p: Vec<(String, String)> = paths(&drs::control::syntax::parse(&moved)).into_iter().map(|p| (p.key, p.path)).collect();
    assert_eq!(p.len(), 2);
    assert_eq!(p[1].0, "mode_reference.file");
    let c2 = Control::from_text_in(&moved, &to).unwrap();
    let sr = c2.model.phases[0].modes.as_ref().unwrap().sign_reference.as_ref().unwrap();
    assert_eq!(sr.amplitudes.len(), 12, "the reference read from its new path {}", p[1].1);
}

#[test]
fn a_joint_control_file_moved() {
    let d = runs().join("scor_rxn");
    let root = fresh("test-relocate-joint");
    let (from, to) = (root.join("here"), root.join("there/deeper"));
    std::fs::create_dir_all(&from).unwrap();
    std::fs::create_dir_all(&to).unwrap();
    let j = drs::joint::Joint::from_pcr(&latin1(&d.join("scor_rxn.pcr")), &d, "scor_rxn", false).unwrap();
    j.save(&from.join("scor_rxn.drc")).unwrap();
    let text = std::fs::read_to_string(from.join("scor_rxn.drc")).unwrap();
    let moved = relocate(&text, &from, &to).unwrap();
    let p: Vec<String> = paths(&drs::control::syntax::parse(&moved)).into_iter().map(|p| p.path).collect();
    assert_eq!(p.len(), j.patterns.len());
    assert!(p.iter().all(|x| x.starts_with("../../here/")), "{p:?}");
    std::fs::write(to.join("scor_rxn.drc"), &moved).unwrap();
    let back = drs::joint::Joint::open(&to.join("scor_rxn.drc")).unwrap();
    assert_eq!(back.patterns.len(), j.patterns.len());
    assert_eq!(back.objective(), j.objective(), "the same patterns and data from the new folder");
}
