//! Plots' PL8 and PL10 and the app's requests 31 and 32 (L199):
//! - PL8: a CSV of several y columns read as several scans, σ paired by name; PL7's `#` lines read back as each member's
//!   file and variable with its source (the mock-up's 41-column series CSV, read before as one scan with 19 dropped
//!   silently, and the core's own wide and long CSV); a column not read said;
//! - PL10: a set shown in d, Q, sin θ/λ and s = 1/d as a view, each member with its own λ (its file's, a control file's,
//!   typed), the members without one listed;
//! - app 31: a results file named by its place on the command line without `--in`; app 32: FullProf's exports under
//!   `Exports/FullProf/`;
//! - PL25 (L208): a CSV or a derived file written as the view shows the set, its normalisation said (the rule and each
//!   member's factor);
//! - PL26 (L208): the core's own CSV read back in d and Q, wide and long, as the set.

use drs::formats::{self, Axis};
use drs::plots::{self, Found, LambdaFrom, XQuantity};
use std::path::{Path, PathBuf};

fn root(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)
}

fn tb() -> PathBuf {
    root("../manuals/fullprof/1-WinPLOTR/1_Tb2pd2sn")
}

fn drs(args: &[&str]) -> std::process::Output {
    std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(args).output().unwrap()
}

fn dir(name: &str) -> PathBuf {
    let d = root("target/test-plots-tables").join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn other<'a>(s: &'a formats::Scan, k: &str) -> Option<&'a str> {
    s.info.other.iter().find(|o| o.0 == k).map(|o| o.1.as_str())
}

#[test]
fn pl8_the_series_csv_as_a_set() {
    // The mock-up's CSV (plots/mockups/out/tb2pd2sn_series.csv, the prototype of PL7): 41 columns, 20 files.
    let csv = root("../plots/mockups/out/tb2pd2sn_series.csv");
    let d = formats::read_path(&csv, None).unwrap();
    assert_eq!((d.format, d.scans.len()), ("table", 20));
    let set = plots::open_buf(&tb().join("Tb_vs_T.buf"), None).unwrap();
    for s in &d.scans {
        let m = set.members.iter().find(|m| m.label == s.label).unwrap_or_else(|| panic!("{}", s.label));
        assert_eq!(s.pattern, m.pattern, "{}: the CSV's numbers are the file's", s.label);
        assert_eq!(other(s, "member format"), Some("free_steps"));
        assert!(other(s, "columns not read").is_none());
    }
    let t02 = d.scans.iter().find(|s| s.label == "tb_t02.dat").unwrap();
    assert_eq!((other(t02, "member value"), other(t02, "member source")), (Some("4"), Some("from the title")));
    assert_eq!(other(&d.scans[0], "member value"), None, "tb_t01.dat: T = not stated");
    // Opened as a set: each member's temperature from its # line, with that line's source; BT has none.
    let s = plots::open(csv.parent().unwrap(), &[("tb2pd2sn_series.csv".into(), None)], None, None).unwrap();
    assert_eq!(s.members.len(), 20);
    let t4 = s.members.iter().find(|m| m.label == "tb_t02.dat").unwrap();
    assert_eq!(t4.variable.as_ref().map(|v| (v.value, v.found.clone())), Some((4.0, Found::CsvLine("from the title".into()))));
    let bt = s.members.iter().find(|m| m.label == "tb_t01.dat").unwrap();
    assert!(bt.variable.is_none() && bt.candidates.is_empty(), "the CSV's title says 1.5–30 K for the set, not BT's temperature");
    assert_eq!(s.members.last().unwrap().label, "tb_t01.dat");
}

#[test]
fn pl8_the_cores_own_csv_read_back() {
    let set = plots::open_buf(&tb().join("Tb_vs_T.buf"), None).unwrap();
    let d = dir("roundtrip");
    // Wide on the common grid.
    let grid = plots::common_grid(&set, plots::GridRange::Overlap, plots::GridStep::Finest).unwrap();
    std::fs::write(d.join("wide.csv"), plots::csv(&set, &grid)).unwrap();
    let back = plots::open(&d, &[("wide.csv".into(), None)], None, None).unwrap();
    assert_eq!(back.members.len(), 20);
    for (a, b) in back.members.iter().zip(&set.members) {
        assert_eq!((&a.label, &a.pattern), (&b.label, &b.pattern));
        match (&a.variable, &b.variable) {
            (Some(x), Some(y)) => assert_eq!((x.value, x.found.clone()), (y.value, Found::CsvLine(y.found.words()))),
            (None, None) => {}
            o => panic!("{o:?}"),
        }
    }
    // Long.
    std::fs::write(d.join("long.csv"), plots::csv_long(&set)).unwrap();
    let long = formats::read_path(&d.join("long.csv"), None).unwrap();
    assert_eq!(long.scans.len(), 20);
    for (s, m) in long.scans.iter().zip(&set.members) {
        assert_eq!((&s.label, &s.pattern), (&m.label, &m.pattern));
    }
    // On the union with a D20 file: empty cells are points a member has not, the rest kept.
    let mixed = plots::open(&root("../manuals/fullprof/1-WinPLOTR"), &[("1_Tb2pd2sn/tb_t02.dat".into(), None), ("2_FP_sequential_refinements/LaMnO3_D20/la0128.dat".into(), None)], None, None).unwrap();
    let u = plots::common_grid(&mixed, plots::GridRange::Union, plots::GridStep::Finest).unwrap();
    std::fs::write(d.join("union.csv"), plots::csv(&mixed, &u)).unwrap();
    let r = formats::read_path(&d.join("union.csv"), None).unwrap();
    let on = plots::on_grid(&mixed, &u);
    for (s, o) in r.scans.iter().zip(&on) {
        assert_eq!(s.pattern.len(), o.y.iter().filter(|v| v.is_some()).count(), "{}", s.label);
    }
    assert_eq!(r.scans[0].pattern.len(), 300);
}

#[test]
fn pl8_columns_paired_and_said() {
    let d = dir("columns");
    // σ paired by name; a σ column with no y said; with no names, a fourth and fifth column said.
    let named = "2theta,I_300K,esd_300K,I_10K,esd_10K,sigma_5K\n10,100,10,400,20,1\n10.1,121,11,441,21,1\n10.2,144,12,484,22,1\n10.3,169,13,529,23,1\n";
    std::fs::write(d.join("named.csv"), named).unwrap();
    let r = formats::read_path(&d.join("named.csv"), None).unwrap();
    assert_eq!(r.scans.iter().map(|s| s.label.as_str()).collect::<Vec<_>>(), ["I_300K", "I_10K"]);
    assert_eq!((r.scans[0].pattern.sigma[1], r.scans[1].pattern.sigma[1]), (11.0, 21.0));
    assert_eq!(other(&r.scans[0], "columns not read"), Some("1 of the file's 6 columns: column 6 \"sigma_5K\""));
    let plain = "10 100 10 7 8\n10.1 121 11 7 8\n10.2 144 12 7 8\n10.3 169 13 7 8\n";
    std::fs::write(d.join("plain.dat"), plain).unwrap();
    let r = formats::read_path(&d.join("plain.dat"), Some("table")).unwrap();
    assert_eq!(r.scans.len(), 1);
    assert_eq!(other(&r.scans[0], "columns not read"), Some("2 of the file's 5 columns: column 4, column 5"));
    // y given by the user in a file of several: one scan, and the others said.
    let o = formats::table::TableOptions { y: Some(3), sigma: Some(Some(4)), ..Default::default() };
    let (ds, lay) = formats::table::read_with(named, &o).unwrap();
    assert_eq!((ds.scans.len(), lay.unread.clone()), (1, vec![1, 2, 5]));
    // The command line names the scans.
    let o = drs(&["data", "info", root("../plots/mockups/out/tb2pd2sn_series.csv").to_str().unwrap()]);
    assert!(String::from_utf8_lossy(&o.stdout).contains("scan 20: tb_t20.dat, 301 points"), "{}", String::from_utf8_lossy(&o.stdout));
}

#[test]
fn pl10_a_set_in_another_quantity() {
    let d = dir("quantities");
    // A UXD stating its λ, a d-spacing xye (no λ needed for d → Q), and Tb's tb_t02.dat (no λ).
    std::fs::write(d.join("cu.uxd"), "_FILEVERSION=1\n_WL1=1.5406\n_STEPTIME=1\n_STEPSIZE=1\n_2THETA=20\n_COUNTS\n10 20 30 40 50\n").unwrap();
    let mut ds = formats::Scan::new("d", Axis::D, drs::data::Pattern { x: vec![1.0, 2.0, 4.0, 5.0, 8.0], y: vec![5.0, 6.0, 7.0, 8.0, 9.0], sigma: vec![1.0; 5] });
    ds.info.title = Some("d".into());
    std::fs::write(d.join("d.xye"), formats::write(&ds, "xye").unwrap()).unwrap();
    std::fs::copy(tb().join("tb_t02.dat"), d.join("tb_t02.dat")).unwrap();
    let set = plots::open(&d, &[("cu.uxd".into(), None), ("d.xye".into(), None), ("tb_t02.dat".into(), None)], None, None).unwrap();
    let k = |f: &str| set.members.iter().position(|m| m.file == f).unwrap();
    // No λ given: tb_t02 listed without; the UXD from its file; d → Q needs none.
    let v = plots::x_view(&set, XQuantity::Q, &[], None);
    assert_eq!(v.without, vec![k("tb_t02.dat")]);
    let cu = &v.members[k("cu.uxd")];
    assert_eq!((cu.lambda, cu.from.clone()), (Some(1.5406), Some(LambdaFrom::File)));
    let q20 = 4.0 * std::f64::consts::PI * 10f64.to_radians().sin() / 1.5406;
    assert!((cu.pattern.as_ref().unwrap().x[0] - q20).abs() < 1e-12);
    let dm = &v.members[k("d.xye")];
    let pi = std::f64::consts::PI;
    assert_eq!((dm.lambda, dm.pattern.as_ref().unwrap().x.clone()), (None, vec![4.0 * pi / 16.0, 4.0 * pi / 10.0, 4.0 * pi / 8.0, 4.0 * pi / 4.0, 4.0 * pi / 2.0]));
    // A control file's λ for those that state none; typed for one member, over its file's.
    let c = plots::x_view(&set, XQuantity::D, &[(k("cu.uxd"), 1.54, LambdaFrom::Typed("a note".into()))], Some((2.426, LambdaFrom::ControlFile("g41.drc".into()))));
    assert!(c.without.is_empty());
    let tb = &c.members[k("tb_t02.dat")];
    assert_eq!(tb.from, Some(LambdaFrom::ControlFile("g41.drc".into())));
    let x = &tb.pattern.as_ref().unwrap().x;
    assert!(x.windows(2).all(|w| w[0] < w[1]), "sorted by d");
    assert!((x[0] - 2.426 / (2.0 * 23f64.to_radians().sin())).abs() < 1e-12 && (x[300] - 2.426 / (2.0 * 8f64.to_radians().sin())).abs() < 1e-12);
    assert_eq!(c.members[k("cu.uxd")].from, Some(LambdaFrom::Typed("a note".into())));
    // sin θ/λ and s; back to 2θ the points as they were.
    let stl = plots::x_view(&set, XQuantity::Stl, &[], Some((2.426, LambdaFrom::Typed("G41".into()))));
    let s = plots::x_view(&set, XQuantity::S, &[], Some((2.426, LambdaFrom::Typed("G41".into()))));
    let (a, b) = (stl.members[k("tb_t02.dat")].pattern.as_ref().unwrap(), s.members[k("tb_t02.dat")].pattern.as_ref().unwrap());
    assert!(a.x.iter().zip(&b.x).all(|(p, q)| *q == 2.0 * p));
    let mut back = set.clone();
    back.members[k("tb_t02.dat")].pattern = c.members[k("tb_t02.dat")].pattern.clone().unwrap();
    back.members[k("tb_t02.dat")].facts.axis = Axis::D;
    let tt = plots::x_view(&back, XQuantity::TwoTheta, &[], Some((2.426, LambdaFrom::Typed("G41".into()))));
    let p = tt.members[k("tb_t02.dat")].pattern.as_ref().unwrap();
    for (x, x0) in p.x.iter().zip(&set.members[k("tb_t02.dat")].pattern.x) {
        assert!((x - x0).abs() < 1e-12, "{x} {x0}");
    }
    // The command line, the λ from a control file.
    let o = drs(&["plots", "view", d.join("tb_t02.dat").to_str().unwrap(), "--x", "stl", "--lambda-from", root("../explore/compare/cases/y2o3_a.drc").to_str().unwrap(), "--json"]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let v: plots::XView = serde_json::from_slice(&o.stdout).unwrap();
    assert_eq!((v.members[0].lambda, v.members[0].from.clone()), (Some(1.5406), Some(LambdaFrom::ControlFile("y2o3_a.drc".into()))));
}

#[test]
fn app31_and_32_on_the_command_line() {
    use drs::formats::project::{describe, describe_path, Description, Role};
    let d = dir("project");
    std::fs::create_dir_all(d.join("Results")).unwrap();
    let f = d.join("Results/a_modes 2026-10-11 12.00.00 results.cif");
    std::fs::copy(root("tests/data/results/a_modes_results.cif"), &f).unwrap();
    let bytes = std::fs::read(&f).unwrap();
    let name = f.to_string_lossy().into_owned();
    let p = formats::Probe::new(&name, &bytes);
    // By its content alone it was a control file (L199: it holds its control file inside); its first line now says.
    match describe(&p) {
        Description::Project(k) => assert_eq!((k.id, k.role), ("results", Role::Output)),
        o => panic!("{o:?}"),
    }
    match describe_path(&f, &p) {
        Description::Project(k) => assert_eq!((k.id, k.role), ("results", Role::Output)),
        o => panic!("{o:?}"),
    }
    let o = drs(&["data", "info", f.to_str().unwrap()]);
    assert!(String::from_utf8_lossy(&o.stdout).contains("(results; Output)"), "{}", String::from_utf8_lossy(&o.stdout));
    // A data file stays data, wherever it lies.
    let dat = root("../manuals/fullprof/runs/Y2O3/y2o3.dat");
    let b = std::fs::read(&dat).unwrap();
    let dname = dat.to_string_lossy().into_owned();
    assert!(matches!(describe_path(&dat, &formats::Probe::new(&dname, &b)), Description::Data { .. }));
    // App 32: the exports for FullProf in Exports/FullProf/, FullProf/ the originals'.
    for (kind, rel) in [("fullprof_pcr_export", "Exports/FullProf/y2o3.pcr"), ("fullprof_reflections_export", "Exports/FullProf/y2o3_structure1.hkl")] {
        assert_eq!(formats::files::kind_of(rel).map(|k| k.id), Some(kind));
    }
    assert!(formats::files::kind_of("FullProf/y2o3.pcr").is_none_or(|k| !k.id.contains("export")));
}

#[test]
fn pl25_the_view_written_as_shown() {
    let set = plots::open_buf(&tb().join("Tb_vs_T.buf"), None).unwrap();
    let rule = plots::Normalisation::Area { from: 20.0, until: 30.0, to: 1.0 };
    let (shown, notes) = set.as_shown(&rule).unwrap();
    assert_eq!(notes[0], "y normalised as the view shows it (the files are as read): the same area, 1, over 2θ (°) 20 to 30");
    assert_eq!(notes.len(), set.members.len() + 1);
    let n = plots::normalise(&set, &rule).unwrap();
    assert_eq!(notes[3], format!("factor of {}: {}", set.members[2].label, n[2].factor.unwrap()));
    // Without a rule: the set as read and no lines.
    let (same, none) = set.as_shown(&plots::Normalisation::None).unwrap();
    assert!(none.is_empty() && same == set);
    // The long CSV as shown: its # lines say how, and read back it is the normalised set, each member's area 1.
    let d = dir("pl25");
    let text = plots::csv_long_with(&shown, &notes);
    assert!(text.lines().nth(1).unwrap().starts_with("# y normalised as the view shows it"));
    std::fs::write(d.join("shown.csv"), &text).unwrap();
    let back = plots::open(&d, &[("shown.csv".into(), None)], None, None).unwrap();
    assert_eq!(back.members.len(), set.members.len());
    for (b, m) in back.members.iter().zip(&shown.members) {
        let o = shown.members.iter().find(|x| x.label == b.label).unwrap();
        assert_eq!((&b.pattern.x, &b.pattern.y), (&o.pattern.x, &o.pattern.y), "{}", b.label);
        let _ = m;
    }
    for m in &shown.members {
        let a = plots::window_value(&m.pattern, &vec![0.0; m.pattern.len()], false, &plots::Window { name: "a".into(), from: 20.0, to: 30.0 }, &[]).unwrap().integral;
        assert!((a - 1.0).abs() < 1e-12, "{}: area {a}", m.label);
    }
    // The wide CSV and a derived file say the same.
    let grid = plots::common_grid(&shown, plots::GridRange::Overlap, plots::GridStep::Finest).unwrap();
    assert!(plots::csv_with(&shown, &grid, &notes).contains(&format!("# {}", notes[1])));
    let dv = plots::derive(&shown, &[0, 1], plots::Derivation::Sum, &grid).unwrap();
    let t = plots::derived_text_with(&shown, &dv, "sum", &notes);
    assert!(t.lines().nth(6).unwrap() == format!("! {}", notes[0]), "after the six header lines");
    // The command line: export and derive with --by, as `drs plots normalise` takes it.
    let buf = tb().join("Tb_vs_T.buf");
    let o = drs(&["plots", "export", buf.to_str().unwrap(), "--csv", "--long", "--by", "area", "--within", "20:30"]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    assert_eq!(String::from_utf8_lossy(&o.stdout), text, "the CLI's CSV as shown is the core's");
    let out = d.join("sum.xys");
    let o = drs(&["plots", "derive", buf.to_str().unwrap(), "--sum", "--members", "1,2", "--by", "area", "--within", "20:30", "--out", out.to_str().unwrap()]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(std::fs::read_to_string(&out).unwrap().contains("! y normalised as the view shows it"));
    // A rule some member cannot follow is refused, naming it: no file of half-normalised numbers.
    assert!(set.as_shown(&plots::Normalisation::Monitor { to: 1.0 }).unwrap_err().contains("the file states no monitor"));
}

#[test]
fn pl26_the_cores_csv_read_back_in_d_and_q() {
    let set = plots::open_buf(&tb().join("Tb_vs_T.buf"), None).unwrap();
    let n = set.members.len();
    let d = dir("pl26");
    for (axis, word) in [(Axis::D, "d"), (Axis::Q, "q")] {
        let c = plots::converted(&set, axis, Some(1.5406)).unwrap();
        let grid = plots::common_grid(&c, plots::GridRange::Overlap, plots::GridStep::Finest).unwrap();
        for (long, text) in [(false, plots::csv(&c, &grid)), (true, plots::csv_long(&c))] {
            let f = format!("{}-{word}.csv", if long { "long" } else { "wide" });
            std::fs::write(d.join(&f), &text).unwrap();
            let back = plots::open(&d, &[(f.clone(), None)], None, None).unwrap();
            assert_eq!(back.members.len(), n, "{f}: the members");
            assert!(back.members.iter().all(|m| m.facts.axis == axis), "{f}: read in {word}");
            if long {
                for b in &back.members {
                    let o = c.members.iter().find(|x| x.label == b.label).unwrap();
                    assert_eq!((&b.pattern.x, &b.pattern.y), (&o.pattern.x, &o.pattern.y), "{f}: {}", b.label);
                }
            }
            let t02 = back.members.iter().find(|m| m.label == "tb_t02.dat").unwrap();
            assert_eq!(t02.variable.as_ref().map(|v| v.value), Some(4.0), "{f}");
        }
    }
}
