//! The web stream's W1 and W3 (I10, R34; L422, L423, L424):
//! - W1: a set from files given by name and content (`plots::open_from`), what a browser reads, with no path: the
//!   Tb₂Pd₂Sn files by name and text give the same set and the same map as by path, a buffer's comments the same values,
//!   a typed value its source; a file the members' operations name taken by its name; a set file's members by name
//!   (`SetFile::open_from`), a missing one refused by name;
//! - W3: a member's other variables with their sources (`variables = { T = { value, source } }` in the set file),
//!   read and written back, carried by the map's rows; a temperature the member's own sources give filling `T` when
//!   the set's variable is another one.

use drs::data_ops::Operation;
use drs::plots::{self, Found, OpenedFile, SetFile, YScale};
use std::path::{Path, PathBuf};

fn root(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)
}

fn tb() -> PathBuf {
    root("../manuals/fullprof/1-WinPLOTR/1_Tb2pd2sn")
}

fn names() -> Vec<String> {
    (1..=20).map(|k| format!("tb_t{k:02}.dat")).collect()
}

fn given(prefix: &str) -> Vec<OpenedFile> {
    names().into_iter().map(|n| OpenedFile { text: std::fs::read_to_string(tb().join(&n)).unwrap(), name: format!("{prefix}{n}"), ..Default::default() }).collect()
}

/// The set without the files' modification times: a set opened by path takes each file's from the disk, one from
/// names and text only those the host gives (W4, L426).
fn untimed(mut set: plots::PlotSet) -> plots::PlotSet {
    for m in set.members.iter_mut() {
        m.modified = None;
    }
    set
}

fn map_of(set: &plots::PlotSet) -> plots::Map {
    let grid = plots::common_grid(set, plots::GridRange::Overlap, plots::GridStep::Finest).unwrap();
    plots::map(set, &grid, YScale::Linear, None).unwrap()
}

#[test]
fn w1_a_set_from_names_and_text_is_the_set_from_paths() {
    let files: Vec<(String, Option<String>)> = names().into_iter().map(|n| (n, None)).collect();
    let by_path = plots::open(&tb(), &files, None, None).unwrap();
    let by_text = plots::open_from(&given(""), None, None).unwrap();
    assert!(by_path.members.iter().all(|m| m.modified.is_some()) && by_text.members.iter().all(|m| m.modified.is_none()), "the disk's times by path, none given by text");
    let by_path = untimed(by_path);
    assert_eq!(by_text, by_path, "every member: file, label, facts, variable and its candidates, pattern");
    assert_eq!(by_text.members.len(), 20);
    assert_eq!(map_of(&by_text), map_of(&by_path));
    // A buffer's comments are the members' comments, the same values as open_buf gives.
    let buf = std::fs::read_to_string(tb().join("Tb_vs_T.buf")).unwrap();
    let b = plots::read_buf(&buf, &|n| tb().join(n).is_file());
    let with_comments: Vec<OpenedFile> = b.entries.iter().map(|(n, c)| OpenedFile { name: n.clone(), text: std::fs::read_to_string(tb().join(n)).unwrap(), comment: c.clone(), variable: None, modified: None }).collect();
    assert_eq!(plots::open_from(&with_comments, None, b.title.clone()).unwrap(), untimed(plots::open_buf(&tb().join("Tb_vs_T.buf"), None).unwrap()));
    // A value typed for a member is its variable, its note the source, first among its candidates.
    let mut typed = given("");
    typed[0].variable = Some((1.5, "WinPLOTR's demo 4: BT is 1.5 K".into()));
    let s = plots::open_from(&typed, None, None).unwrap();
    let bt = s.members.iter().find(|m| m.file == "tb_t01.dat").unwrap();
    assert_eq!(bt.variable.as_ref().map(|v| (v.value, v.found.clone())), Some((1.5, Found::Typed("WinPLOTR's demo 4: BT is 1.5 K".into()))));
    assert_eq!(s.members[0].file, "tb_t01.dat", "sorted by the variable, BT at 1.5 K first");
    // A file the operations name, taken from the files given by its name; refused by name when absent.
    let ops = |f: &str| if f == "tb_t10.dat" { vec![Operation::Subtract { file: "tb_t20.dat".into(), reader: None, scale: 1.0 }] } else { Vec::new() };
    let d = plots::open_from_with(&given(""), None, None, &ops).unwrap();
    let p = untimed(plots::open_with(&tb(), &files, None, None, &ops).unwrap());
    assert_eq!(d, p);
    let t10 = d.members.iter().find(|m| m.file == "tb_t10.dat").unwrap();
    let raw10 = by_text.members.iter().find(|m| m.file == "tb_t10.dat").unwrap();
    assert!(t10.pattern.y.iter().zip(&raw10.pattern.y).any(|(a, b)| a != b), "the subtraction applied");
    let only_ten: Vec<OpenedFile> = given("").into_iter().filter(|f| f.name == "tb_t10.dat").collect();
    assert_eq!(plots::open_from_with(&only_ten, None, None, &ops).unwrap_err(), "tb_t10.dat: operation 1 (subtract): tb_t20.dat: an operation names it, and no file of that name was given");
}

const SET: &str = r#"format = "drs plot set 1"
title = "Tb2Pd2Sn by name"
variable = { name = "T", unit = "K" }

[[member]]
file = "Data/tb_t20.dat"

[[member]]
file = "Data/tb_t01.dat"
variable = { value = 1.5, source = "WinPLOTR demo 4" }

[[member]]
file = "Data/tb_t02.dat"
"#;

#[test]
fn w1_a_set_file_from_names_and_text() {
    let f = SetFile::read(SET).unwrap();
    let by_text = f.open_from(&given("Data/"), None, &|_| Vec::new()).unwrap();
    // The same set from the project's folder on disk: Data/ beside Plots/.
    let d = root("target/test-plots-open-from/project");
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(d.join("Data")).unwrap();
    for n in names() {
        std::fs::copy(tb().join(&n), d.join("Data").join(&n)).unwrap();
    }
    let by_path = untimed(f.open(&d, None, &|_| Vec::new()).unwrap());
    assert_eq!(by_text, by_path);
    assert_eq!(by_text.members.iter().map(|m| m.label.as_str()).collect::<Vec<_>>(), ["tb_t01.dat", "tb_t02.dat", "tb_t20.dat"]);
    assert_eq!(by_text.members[0].variable.as_ref().map(|v| v.value), Some(1.5));
    // Only the files the set names are members, whatever else is given; a member's file not given is refused by name.
    let few: Vec<OpenedFile> = given("Data/").into_iter().filter(|g| g.name != "Data/tb_t02.dat").collect();
    assert_eq!(f.open_from(&few, None, &|_| Vec::new()).unwrap_err(), "Data/tb_t02.dat: the set file names it, and no file of that name was given");
}

/// W3: a set against time, the hold's temperatures as the members' other variable.
const AGAINST_TIME: &str = r#"format = "drs plot set 1"
title = "a hold against time"
variable = { name = "t", unit = "h" }

[[member]]
file = "Data/tb_t02.dat"
variable = { value = 0.5, source = "the run's log" }
variables = { T = { value = 1270.0, source = "the pyrometer, °C" } }

[[member]]
file = "Data/tb_t03.dat"
variable = { value = 1.0, source = "the run's log" }
variables = { T = { value = 1272.5, source = "the pyrometer, °C" }, p = { value = 1.0, source = "the cell, bar" } }

[[member]]
file = "Data/tb_t04.dat"
variable = { value = 1.5, source = "the run's log" }
"#;

#[test]
fn w3_other_variables_with_their_sources() {
    let f = SetFile::read(AGAINST_TIME).unwrap();
    assert_eq!(f.members[0].variables, vec![("T".to_string(), 1270.0, "the pyrometer, °C".to_string())]);
    assert_eq!(f.members[1].variables.len(), 2);
    assert!(f.members[2].variables.is_empty());
    // Written and read back: the same file, the same text again; a file without them writes none.
    let text = f.to_text();
    assert_eq!(SetFile::read(&text).unwrap(), f);
    assert!(text.contains("variables = { T = { value = 1270.0, source = \"the pyrometer, °C\" } }\n"), "{text}");
    assert!(!SetFile::read(SET).unwrap().to_text().contains("variables"));
    // Opened: the typed others with their sources; a member typing none takes the temperature its title gives (PL2).
    let set = f.open_from(&given("Data/"), None, &|_| Vec::new()).unwrap();
    assert_eq!((set.variable_name.as_str(), set.variable_unit.as_str()), ("t", "h"));
    let others = |k: usize| set.members[k].others.iter().map(|o| (o.name.clone(), o.value, o.found.words())).collect::<Vec<_>>();
    assert_eq!(others(0), [("T".to_string(), 1270.0, "typed: the pyrometer, °C".to_string())]);
    assert_eq!(others(1), [("T".to_string(), 1272.5, "typed: the pyrometer, °C".to_string()), ("p".to_string(), 1.0, "typed: the cell, bar".to_string())]);
    assert_eq!(others(2), [("T".to_string(), 8.0, "the title's \"T.=.8.K\"".to_string())]);
    // The map's rows carry them, each with its source in words, for the right-hand axis.
    let m = map_of(&set);
    assert_eq!(m.rows.iter().map(|r| r.value).collect::<Vec<_>>(), [0.5, 1.0, 1.5]);
    assert_eq!((m.rows[0].others[0].name.as_str(), m.rows[0].others[0].value, m.rows[0].others[0].source.as_str()), ("T", 1270.0, "typed: the pyrometer, °C"));
    assert_eq!(m.rows[2].others[0].source, "the title's \"T.=.8.K\"");
    // The set written from the set keeps the typed others, not the ones a file read again gives.
    let back = SetFile::of(&set);
    assert_eq!(back.members.iter().map(|m| m.variables.len()).collect::<Vec<_>>(), [1, 2, 0]);
    // A set whose variable is T takes no other T from the members' sources: the map's rows carry nothing.
    let by_t = plots::open_from(&given(""), None, None).unwrap();
    assert!(by_t.members.iter().all(|m| m.others.is_empty()));
    assert!(map_of(&by_t).rows.iter().all(|r| r.others.is_empty()));
    // Refused: an other variable without its source, the set's own variable as an other, one given twice.
    let wrong = |a: &str, b: &str| SetFile::read(&AGAINST_TIME.replace(a, b)).unwrap_err()[0].message.clone();
    assert!(wrong("T = { value = 1270.0, source = \"the pyrometer, °C\" }", "T = { value = 1270.0 }").contains("says its source"));
    assert!(wrong("T = { value = 1270.0, source = \"the pyrometer, °C\" }", "t = { value = 3.0, source = \"x\" }").contains("is the set's variable"));
    assert!(wrong("p = { value = 1.0, source = \"the cell, bar\" }", "T = { value = 1.0, source = \"again\" }").contains("given twice"));
    // The JSON of a member and of a map row: the others beside the rest, absent ones as empty lists.
    let row = serde_json::to_value(&m.rows[1]).unwrap();
    assert_eq!(row["others"][1], serde_json::json!({ "name": "p", "value": 1.0, "source": "typed: the cell, bar", "unit": "" }));
    let mut old = serde_json::to_value(&m.rows[1]).unwrap();
    old.as_object_mut().unwrap().remove("others");
    assert!(serde_json::from_value::<plots::MapRow>(old).unwrap().others.is_empty());
}
