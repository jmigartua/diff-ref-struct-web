//! Plots' requests PL1–PL3 (R32, L170 step 2, L172), on the two series of `manuals/fullprof/1-WinPLOTR/`: Tb₂Pd₂Sn
//! (20 FullProf free-format files, two WinPLOTR buffers) and LaMnO₃ on D20 (71 ILL files).
//! - PL1: a buffer opened as a set, every member with the reader's facts; a buffer's title and comments;
//! - PL2: each member's temperature with its source (Tb₂Pd₂Sn's titles "T.=.4.K", D20's header), none for the file
//!   that says "BT", never a guess; the members sorted by it, a typed value with its note;
//! - PL3: the displacements of a stack and of a projection, on the linear, √ and log scales, by member or in
//!   proportion to T, with overrides and the background displaced alike;
//! - the command line, `drs plots open` and `drs plots displace`.

use drs::plots::{self, Displacement, Found, YScale};
use std::path::{Path, PathBuf};

fn root(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)
}

fn tb() -> PathBuf {
    root("../manuals/fullprof/1-WinPLOTR/1_Tb2pd2sn")
}

fn d20() -> PathBuf {
    root("../manuals/fullprof/1-WinPLOTR/2_FP_sequential_refinements/LaMnO3_D20")
}

#[test]
fn pl1_a_buffer_is_a_set() {
    for buf in ["Tb_vs_T.buf", "data_buffer_file.buf"] {
        let set = plots::open_buf(&tb().join(buf), None).unwrap();
        assert_eq!(set.members.len(), 20, "{buf}");
        for m in &set.members {
            let f = &m.facts;
            assert_eq!((f.format.as_str(), f.points, f.range, m.scans), ("free_steps", 301, (16.0, 46.0), 1), "{}", m.file);
            assert!((f.step.unwrap() - 0.1).abs() < 1e-12);
            assert_eq!(m.pattern.len(), 301);
        }
    }
    // A buffer's title and comments: the user's guide's example 7 (p. 19, §5.1.4: "file_1.dat      T=5K",
    // "file_2.dat      T = 10 K"), the news's title line (July 2013), empty lines (May 1998); a name with a blank
    // kept whole when that file exists (ours).
    let exists = |n: &str| ["file_1.dat", "file_2.dat"].contains(&n);
    let b = plots::read_buf("TITLE my plot\n\nfile_1.dat      T=5K\nfile_2.dat      T = 10 K\n", &exists);
    assert_eq!(b.title.as_deref(), Some("my plot"));
    assert_eq!(b.entries.iter().map(|e| (e.0.as_str(), e.1.as_deref().and_then(plots::temperature_in).map(|t| t.0))).collect::<Vec<_>>(), vec![("file_1.dat", Some(5.0)), ("file_2.dat", Some(10.0))]);
    let exists = |n: &str| ["tb_t01.dat", "my file.dat"].contains(&n);
    let b = plots::read_buf("TITLE Tb2Pd2Sn on G4.1\ntb_t01.dat   T = 1.5 K\nmy file.dat\n\ntb_t02.dat\n", &exists);
    assert_eq!(b.title.as_deref(), Some("Tb2Pd2Sn on G4.1"));
    assert_eq!(b.entries, vec![("tb_t01.dat".to_string(), Some("T = 1.5 K".to_string())), ("my file.dat".to_string(), None), ("tb_t02.dat".to_string(), None)]);
    let b = plots::read_buf("COMMENT: a series\nx.dat\n", &|_| true);
    assert_eq!(b.title.as_deref(), Some("a series"));
    // The set's shape, through JSON and back (the app and drs share it).
    let set = plots::open_buf(&tb().join("Tb_vs_T.buf"), None).unwrap();
    let text = serde_json::to_string(&set).unwrap();
    let back: plots::PlotSet = serde_json::from_str(&text).unwrap();
    assert_eq!(back, set);
}

#[test]
fn pl2_the_temperature_and_its_source() {
    // The parsers.
    for (t, v) in [("Tb2Pd2Sn................T.=.4.K...........29/06/96", Some(4.0)), ("T = 300 K", Some(300.0)), ("Temp: 1.5K", Some(1.5)), ("temperature=12.25 K", Some(12.25)), ("T=4.5K run 3", Some(4.5)), ("Tb2Pd2Sn..................BT..............29/06/96  INSTRM=0", None), ("Ti = 300 K", None), ("T = 300", None)] {
        assert_eq!(plots::temperature_in(t).map(|x| x.0), v, "{t}");
    }
    for (n, v) in [("sample_300K.dat", Some(300.0)), ("run_T300.dat", Some(300.0)), ("x_1p5K.xye", Some(1.5)), ("tb_t01.dat", None), ("LaB6.dat", None), ("p_10Kbar.dat", None), ("la0128.dat", None), ("T20K.dat", Some(20.0))] {
        assert_eq!(plots::temperature_in_name(n).map(|x| x.0), v, "{n}");
    }
    // Tb₂Pd₂Sn: the titles say 4, 6, … 30 K; the first file says "BT" and gets none, after the others.
    let mut set = plots::open_buf(&tb().join("Tb_vs_T.buf"), None).unwrap();
    let ts: Vec<Option<f64>> = set.members.iter().map(|m| m.variable.as_ref().map(|v| v.value)).collect();
    let want: Vec<Option<f64>> = [4., 6., 8., 10., 12., 14., 16., 18., 20., 21., 22., 23., 24., 25., 26., 27., 28., 29., 30.].iter().map(|&t| Some(t)).chain([None]).collect();
    assert_eq!(ts, want);
    assert!(set.members[..19].iter().all(|m| matches!(&m.variable.as_ref().unwrap().found, Found::Title(t) if t.starts_with("T.=.") && t.ends_with(".K"))));
    assert_eq!(set.members[19].file, "tb_t01.dat");
    assert_eq!(set.without_variable(), vec![19]);
    // Typed, with its note: sorted first (1.5 K).
    set.set_variable(19, 1.5, "WinPLOTR's demo 4 gives 1.5–30 K; the file says BT").unwrap();
    assert_eq!(set.members[0].file, "tb_t01.dat");
    assert!(matches!(&set.members[0].variable.as_ref().unwrap().found, Found::Typed(n) if n.contains("BT")));
    assert!(set.without_variable().is_empty());
    // A buffer's comment is taken before the title, and the disagreement kept.
    let dir = root("target/test-plots/comments");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    for f in ["tb_t01.dat", "tb_t02.dat"] {
        std::fs::copy(tb().join(f), dir.join(f)).unwrap();
    }
    std::fs::write(dir.join("two.buf"), "TITLE two files\ntb_t01.dat  T = 1.5 K\ntb_t02.dat  T = 4.2 K\n").unwrap();
    let set = plots::open_buf(&dir.join("two.buf"), None).unwrap();
    assert_eq!(set.title.as_deref(), Some("two files"));
    let m = &set.members[1];
    assert_eq!((m.file.as_str(), m.variable.as_ref().unwrap().value), ("tb_t02.dat", 4.2));
    assert!(matches!(m.variable.as_ref().unwrap().found, Found::BufComment(_)));
    assert_eq!(m.candidates.iter().map(|c| c.value).collect::<Vec<_>>(), vec![4.2, 4.0]);
    // LaMnO₃ on D20: the header's temperature, la0128 79.618 K (the request's example); every file has one, sorted.
    let mut files: Vec<String> = std::fs::read_dir(d20()).unwrap().filter_map(|e| e.ok()).map(|e| e.file_name().to_string_lossy().into_owned()).filter(|n| n.ends_with(".dat")).collect();
    files.sort();
    let set = plots::open(&d20(), &files.iter().map(|f| (f.clone(), None)).collect::<Vec<_>>(), None, None).unwrap();
    assert_eq!(set.members.len(), 71);
    assert!(set.members.iter().all(|m| m.variable.as_ref().is_some_and(|v| v.found == Found::Header) && m.facts.format == "ill_d1b_d20"));
    let la0128 = set.members.iter().find(|m| m.file == "la0128.dat").unwrap();
    assert_eq!(la0128.variable.as_ref().unwrap().value, 79.618);
    assert!(set.members.windows(2).all(|w| w[0].variable.as_ref().unwrap().value <= w[1].variable.as_ref().unwrap().value));
    let (lo, hi) = (set.members[0].variable.as_ref().unwrap().value, set.members[70].variable.as_ref().unwrap().value);
    println!("D20: 71 files from {lo} K to {hi} K");
}

#[test]
fn pl3_displacements() {
    let mut set = plots::open_buf(&tb().join("Tb_vs_T.buf"), None).unwrap();
    // A stack: k·dy on the linear scale; the baseline the shift.
    let d = plots::displace(&set, Displacement::Step { dx: 0.0, dy: 3000.0 }, YScale::Linear, &[], None);
    for (k, (m, o)) in set.members.iter().zip(&d).enumerate() {
        let o = o.as_ref().unwrap();
        assert_eq!((o.shift, o.baseline), ((0.0, 3000.0 * k as f64), 3000.0 * k as f64));
        assert!(o.y.iter().zip(&m.pattern.y).all(|(a, b)| *a == b + 3000.0 * k as f64) && o.x == m.pattern.x);
    }
    // On √: the shift added to √y, so equal gaps on the screen.
    let d = plots::displace(&set, Displacement::Step { dx: 0.0, dy: 10.0 }, YScale::Sqrt, &[], None);
    let o = d[3].as_ref().unwrap();
    assert!(o.y.iter().zip(&set.members[3].pattern.y).all(|(a, b)| (a - (b.sqrt() + 30.0)).abs() < 1e-12));
    // A projection: k·dx too; log₁₀ y.
    let d = plots::displace(&set, Displacement::Step { dx: 0.5, dy: 0.1 }, YScale::Log, &[], None);
    let o = d[2].as_ref().unwrap();
    assert!(o.x.iter().zip(&set.members[2].pattern.x).all(|(a, b)| (a - (b + 1.0)).abs() < 1e-12));
    assert!(o.y.iter().zip(&set.members[2].pattern.y).all(|(a, b)| (a - (b.log10() + 0.2)).abs() < 1e-12));
    assert!(YScale::Log.apply(0.0).is_nan() && YScale::Sqrt.apply(-4.0) == -2.0);
    // In proportion to T, from the first member's (4 K): 30 K is 26·dy up; the file without a value is not placed
    // until one is typed.
    let d = plots::displace(&set, Displacement::Variable { dx: 0.0, dy: 100.0 }, YScale::Linear, &[], None);
    assert_eq!(d[18].as_ref().unwrap().shift, (0.0, 2600.0));
    assert!(d[19].is_none());
    set.set_variable(19, 1.5, "typed").unwrap();
    let d = plots::displace(&set, Displacement::Variable { dx: 0.0, dy: 100.0 }, YScale::Linear, &[], None);
    assert_eq!(d[1].as_ref().unwrap().shift, (0.0, 250.0));
    // A member's own override; a background displaced alike.
    let bg: Vec<Vec<f64>> = set.members.iter().map(|m| vec![1000.0; m.pattern.len()]).collect();
    let d = plots::displace(&set, Displacement::Step { dx: 0.0, dy: 3000.0 }, YScale::Linear, &[(5, (0.0, -1.0))], Some(&bg));
    let o = d[5].as_ref().unwrap();
    assert_eq!(o.shift, (0.0, -1.0));
    assert!(o.background.as_ref().unwrap().iter().all(|&b| b == 999.0));
    assert!(d[6].as_ref().unwrap().background.as_ref().unwrap().iter().all(|&b| b == 1000.0 + 18000.0));
}

#[test]
fn the_command_line() {
    let run = |args: &[&str]| std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(args).output().unwrap();
    let buf = tb().join("Tb_vs_T.buf");
    let o = run(&["plots", "open", buf.to_str().unwrap(), "--set", "20=1.5:the demo"]);
    let out = String::from_utf8_lossy(&o.stdout);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(out.lines().nth(1).unwrap().contains("tb_t01.dat") && out.contains("T = 1.5 K (typed: the demo)") && out.contains("T = 4 K (the title's \"T.=.4.K\")"), "{out}");
    let o = run(&["plots", "open", buf.to_str().unwrap(), "--json"]);
    let set: plots::PlotSet = serde_json::from_slice(&o.stdout).unwrap();
    assert_eq!(set.members.len(), 20);
    let o = run(&["plots", "displace", buf.to_str().unwrap(), "--step", "3000", "--scale", "sqrt"]);
    let out = String::from_utf8_lossy(&o.stdout);
    assert!(o.status.success() && out.contains("shift x +0.0000 y +57000.0000"), "{out}");
    let o = run(&["plots", "displace", buf.to_str().unwrap(), "--per-unit", "100"]);
    assert!(String::from_utf8_lossy(&o.stdout).contains("not placed: no value of T"));
}

/// PL4: the grid, the resampling with σ, the map and the difference map.
#[test]
fn pl4_grid_and_map() {
    let set = plots::open_buf(&tb().join("Tb_vs_T.buf"), None).unwrap();
    assert!(plots::same_grid(&set));
    // The members' own points when they share them.
    let grid = plots::common_grid(&set, plots::GridRange::Overlap, plots::GridStep::Finest).unwrap();
    assert_eq!(grid, set.members[0].pattern.x);
    assert!(plots::on_grid(&set, &grid).iter().all(|r| r.how == plots::Resampling::AsMeasured));
    // A coarser grid (0.2°, every other point): each new point the mean over its bin, the middle point weighted 0.1
    // and its two neighbours 0.05 each (the parts of their intervals inside the bin); σ² = Σw²σ²/(Σw)².
    let g2 = plots::common_grid(&set, plots::GridRange::Overlap, plots::GridStep::Given(0.2)).unwrap();
    assert_eq!(g2.len(), 151);
    let p = &set.members[0].pattern;
    let r = plots::resample(p, &g2);
    assert_eq!(r.how, plots::Resampling::BinAverage);
    for j in 1..150 {
        let i = 2 * j;
        let y = (0.1 * p.y[i] + 0.05 * p.y[i - 1] + 0.05 * p.y[i + 1]) / 0.2;
        let s = (0.01 * p.sigma[i].powi(2) + 0.0025 * (p.sigma[i - 1].powi(2) + p.sigma[i + 1].powi(2))).sqrt() / 0.2;
        assert!((r.y[j].unwrap() - y).abs() < 1e-9 * y && (r.sigma[j].unwrap() - s).abs() < 1e-9 * s, "{j}");
    }
    // A finer grid: linear interpolation, exact at the grid's points that are the member's.
    let g3 = plots::common_grid(&set, plots::GridRange::Overlap, plots::GridStep::Given(0.05)).unwrap();
    let r = plots::resample(p, &g3);
    assert_eq!(r.how, plots::Resampling::Interpolated);
    assert!((r.y[2].unwrap() - p.y[1]).abs() < 1e-9 && (r.y[1].unwrap() - (p.y[0] + p.y[1]) / 2.0).abs() < 1e-9);
    assert!((r.sigma[1].unwrap() - ((p.sigma[0].powi(2) + p.sigma[1].powi(2)) / 4.0).sqrt()).abs() < 1e-9);
    // Two series of different ranges: the union with None outside each member, the overlap only where both are.
    let mixed = plots::open(&root("../manuals/fullprof/1-WinPLOTR"), &[("1_Tb2pd2sn/tb_t02.dat".into(), None), ("2_FP_sequential_refinements/LaMnO3_D20/la0128.dat".into(), None)], None, None).unwrap();
    assert!(!plots::same_grid(&mixed));
    let u = plots::common_grid(&mixed, plots::GridRange::Union, plots::GridStep::Finest).unwrap();
    let o = plots::common_grid(&mixed, plots::GridRange::Overlap, plots::GridStep::Finest).unwrap();
    assert!((u[0] + 3.101).abs() < 1e-9 && (o[0] - 16.0).abs() < 1e-9 && (o[o.len() - 1] - 46.0).abs() < 1e-6);
    let tb_on_u = plots::on_grid(&mixed, &u).into_iter().find(|r| r.y.iter().any(|v| v.is_none())).unwrap();
    // The union grid starts at D20's −3.101° with D20's step 0.1°, so its points sit 0.001° off Tb's (16.000 … 46.000):
    // 16.099 … 45.999 fall inside Tb's range, 300 of them, and 15.999 just outside.
    assert_eq!(tb_on_u.y.iter().filter(|v| v.is_some()).count(), 300);
    // The map: 19 rows at 4…30 K, edges at the midpoints (20 K from 19 to 20.5), tb_t01 apart; the colour range the
    // smallest and largest count.
    let mp = plots::map(&set, &grid, YScale::Linear, None).unwrap();
    assert_eq!((mp.rows.len(), mp.apart.clone()), (19, vec![19]));
    assert_eq!(mp.rows[0].edges, (3.0, 5.0));
    let r20 = mp.rows.iter().find(|r| r.value == 20.0).unwrap();
    assert_eq!(r20.edges, (19.0, 20.5));
    assert_eq!(mp.rows[18].edges, (29.5, 30.5));
    // The colour range over every value drawn, the row apart (BT, PL18, L192) included: its largest count, 68658, is
    // above the 19 rows' 68455.
    let all: Vec<f64> = set.members.iter().flat_map(|m| m.pattern.y.clone()).collect();
    assert_eq!(mp.range, (all.iter().copied().fold(f64::INFINITY, f64::min), all.iter().copied().fold(f64::NEG_INFINITY, f64::max)));
    assert_eq!((mp.x_edges.len(), mp.x_edges[0], mp.x_edges[301]), (302, 15.95, 46.05));
    // On √ and log the range is the scaled one.
    let ms = plots::map(&set, &grid, YScale::Sqrt, None).unwrap();
    assert!((ms.range.1 - mp.range.1.sqrt()).abs() < 1e-9);
    // The difference map less the 30 K pattern (the paramagnetic one, F4): its own row 0, a range symmetric about 0.
    let k30 = set.members.iter().position(|m| m.variable.as_ref().is_some_and(|v| v.value == 30.0)).unwrap();
    let dm = plots::map(&set, &grid, YScale::Linear, Some(k30)).unwrap();
    assert_eq!(dm.range.0, -dm.range.1);
    assert!(dm.z[18].iter().all(|v| *v == Some(0.0)));
    let m4 = dm.z[0].iter().zip(&set.members[0].pattern.y).zip(&set.members[k30].pattern.y).all(|((d, a), b)| *d == Some(a - b));
    assert!(m4);
    assert!(plots::map(&set, &grid, YScale::Log, Some(k30)).is_err());
    println!("map {} to {}; difference map ±{}", mp.range.0, mp.range.1, dm.range.1);
}

/// PL9: sum, averages and difference with σ, written as a file that names its members and read back.
#[test]
fn pl9_derived_patterns() {
    let set = plots::open_buf(&tb().join("Tb_vs_T.buf"), None).unwrap();
    let grid = set.members[0].pattern.x.clone();
    let (a, b) = (&set.members[0].pattern, &set.members[1].pattern);
    // The sum of two counted patterns: σ² = σ₁² + σ₂², here √(y₁ + y₂), the counts' own.
    let s = plots::derive(&set, &[0, 1], plots::Derivation::Sum, &grid).unwrap();
    assert_eq!((s.pattern.len(), s.dropped), (301, 0));
    for i in 0..301 {
        assert_eq!(s.pattern.y[i], a.y[i] + b.y[i]);
        assert!((s.pattern.sigma[i] - (a.y[i] + b.y[i]).sqrt()).abs() < 1e-9);
    }
    let av = plots::derive(&set, &[0, 1], plots::Derivation::Average, &grid).unwrap();
    assert!((0..301).all(|i| av.pattern.y[i] == (a.y[i] + b.y[i]) / 2.0 && (av.pattern.sigma[i] - s.pattern.sigma[i] / 2.0).abs() < 1e-12));
    let w = plots::derive(&set, &[0, 1], plots::Derivation::WeightedAverage, &grid).unwrap();
    for i in 0..301 {
        let (wa, wb) = (1.0 / a.sigma[i].powi(2), 1.0 / b.sigma[i].powi(2));
        assert!((w.pattern.y[i] - (wa * a.y[i] + wb * b.y[i]) / (wa + wb)).abs() < 1e-9 && (w.pattern.sigma[i] - 1.0 / (wa + wb).sqrt()).abs() < 1e-12);
    }
    let d = plots::derive(&set, &[0, 18], plots::Derivation::Difference, &grid).unwrap();
    let c = &set.members[18].pattern;
    assert!((0..301).all(|i| d.pattern.y[i] == a.y[i] - c.y[i]));
    assert!(plots::derive(&set, &[0, 1, 2], plots::Derivation::Difference, &grid).is_err());
    // On a grid wider than a member: the points outside it left out and counted.
    let wide: Vec<f64> = (0..=340).map(|k| 14.0 + 0.1 * k as f64).collect();
    let s2 = plots::derive(&set, &[0, 1], plots::Derivation::Sum, &wide).unwrap();
    assert_eq!((s2.pattern.len(), s2.dropped), (301, 40));
    // Written and read back by the registry: the same numbers, the members named with their temperatures.
    let dir = root("target/test-plots/derived");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let text = plots::derived_text(&set, &d, "Tb2Pd2Sn");
    assert!(text.starts_with("XYDATA\nTITLE: Tb2Pd2Sn") && text.contains("! the difference of 2 members, the first less the second"), "{text}");
    assert!(text.contains("! member 1: tb_t02.dat (scan 1), T = 4 K (the title's \"T.=.4.K\")") && text.contains("! member 2: tb_t20.dat (scan 1), T = 30 K"), "{text}");
    assert!(d.pattern.y.iter().any(|&y| y < 0.0), "a difference with negative points, which the file must keep");
    let path = dir.join(drs::formats::files::path("derived_data", &[("set", "Tb2Pd2Sn"), ("operation", "difference")]).unwrap());
    std::fs::create_dir_all(path.parent().unwrap()).unwrap();
    std::fs::write(&path, &text).unwrap();
    assert_eq!(drs::formats::files::kind_of("Plots/Tb2Pd2Sn difference.xys").unwrap().id, "derived_data");
    let back = drs::formats::read_path(&path, None).unwrap();
    assert_eq!(back.scans[0].pattern, d.pattern);
    // The command line.
    let o = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).current_dir(&dir).args(["plots", "derive", tb().join("Tb_vs_T.buf").to_str().unwrap(), "--sum", "--members", "1,2", "--name", "Tb"]).output().unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let back = drs::formats::read_path(&dir.join("Plots/Tb sum.xys"), None).unwrap();
    assert_eq!(back.scans[0].pattern, s.pattern);
    let o = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(["plots", "map", tb().join("Tb_vs_T.buf").to_str().unwrap(), "--reference", "19", "--json"]).output().unwrap();
    let m: plots::Map = serde_json::from_slice(&o.stdout).unwrap();
    assert_eq!(m.reference, Some(18));
}
