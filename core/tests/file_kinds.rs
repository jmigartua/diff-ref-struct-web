//! Every file of a FullProf project named by what it is (W5, L080, L081): the classifier (`formats::project`) against
//! the expected list `tests/data/file_kinds.txt` (written by `make_file_kinds.py` from a survey of the files of
//! `manuals/fullprof/runs`, the tutorials and FP_tests, apart from the classifier). Each kind carries the source that
//! says what it is, the manual's section (Input or Output files) for FullProf's own; a file nothing documents is
//! named "unknown"; no file with one of FullProf's extensions is read as a table of columns.

use drs::formats::project::{self, Description};
use drs::formats::Probe;
use std::path::Path;

/// Extensions of FullProf and its suite (the manual's Input and Output files, the News).
const FULLPROF_EXTENSIONS: &[&str] = &[
    "pcr", "new", "bac", "hkl", "int", "irf", "out", "sum", "log", "prf", "sym", "fou", "inp", "ins", "mic", "dis", "sub", "atm", "rpa", "sav", "sim", "sch",
    "shp", "cor", "hlp", "fst", "mcif", "vesta", "ggrid", "pgrid", "bin", "str", "siz", "cfl", "bvs", "tpcr", "spr", "seq", "so1", "so2", "anl", "cfg", "bgr",
    "buf", "aps", "dic", "ind", "ord", "ucp", "cry", "cnd", "spg", "fch", "fhz", "red", "flts", "msf", "ort", "mlf",
];

#[test]
fn file_kinds() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof");
    let text = std::fs::read_to_string(Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/file_kinds.txt")).unwrap();
    let mut sources = std::collections::BTreeMap::new();
    let (mut n, mut wrong, mut missing) = (0, Vec::new(), 0);
    let mut tally = std::collections::BTreeMap::<&str, usize>::new();
    for l in text.lines().filter(|l| !l.starts_with('#') && !l.is_empty()) {
        if let Some(k) = l.strip_prefix("kind ") {
            let (id, src) = k.split_once('\t').unwrap();
            sources.insert(id.to_string(), src.to_string());
            continue;
        }
        let (path, want) = l.split_once('\t').unwrap();
        let Ok(bytes) = std::fs::read(root.join(path)) else {
            missing += 1;
            continue;
        };
        let p = Probe::new(path, &bytes);
        let d = project::describe(&p);
        let ext = path.rsplit('.').next().unwrap_or("").to_lowercase();
        if let Description::Data { id, .. } = &d {
            assert!(!(FULLPROF_EXTENSIONS.contains(&ext.as_str()) && *id == "table"), "{path}: a FullProf file read as a table");
        }
        if d.id() != want {
            wrong.push(format!("{path}: {} (expected {want})", d.id()));
        }
        if let Description::Project(k) = d {
            let src = sources.get(k.id).unwrap_or_else(|| panic!("{path}: no source listed for {}", k.id));
            assert_eq!(k.source, src, "{path}: the source of {}", k.id);
            assert!(!k.name.is_empty());
        }
        *tally.entry(d.id()).or_default() += 1;
        n += 1;
    }
    // FullProf's own kinds cite the manual's section.
    for (id, src) in &sources {
        if ["fullprof_listing", "fullprof_summary", "fullprof_profile", "fullprof_symmetry", "fullprof_fourier", "fullprof_microstructure", "fullprof_distances", "fullprof_reflections"].contains(&id.as_str()) {
            assert!(src.contains("Output files"), "{id}: {src}");
        }
        if ["fullprof_pcr", "fullprof_background", "fullprof_resolution", "fullprof_integrated_intensities"].contains(&id.as_str()) {
            assert!(src.contains("Input files"), "{id}: {src}");
        }
    }
    println!("{n} files, {} kinds: {tally:?}", tally.len());
    if missing > 0 {
        println!("{missing} files of the list are not on disk");
    }
    assert!(wrong.is_empty(), "{} files named otherwise:\n{}", wrong.len(), wrong.join("\n"));
    assert!(n > 2600, "{n}");
}
