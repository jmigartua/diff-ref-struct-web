//! Plots' PL6, PL7, PL21 and PL22 (L170 step 5, L195):
//! - PL6: normalisation as a view (monitor, time, maximum, maximum in a range, the same area over a range, a peak's
//!   integral, a reference member's area), the files and their operations untouched; the constant counting time of a
//!   multi-scan UXD (counts and counts per second) and XRDML, overlaps keeping the first scan (WinPLOTR's news);
//! - PL7: a set as CSV, wide on the common grid or long, in 2θ, d or Q; the single-scan writers in d and Q, read back;
//! - PL21: FullProf's temperatures file (its news of 25 September 2013, its own example), each member its file names
//!   taking the value, a name given twice with two values taking none;
//! - PL22: a window's value flagged below 2σ, the threshold the core's;
//! - the command line: `drs plots normalise`, `export --long --x d`, `--temperatures`, `drs data constant-time`.

use drs::data::Pattern;
use drs::formats::{self, Axis};
use drs::plots::{self, Found, Normalisation, Window, WindowBackground};
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

fn drs(args: &[&str]) -> std::process::Output {
    std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(args).output().unwrap()
}

fn dir(name: &str) -> PathBuf {
    let d = root("target/test-plots-views").join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn area(p: &Pattern, from: f64, to: f64) -> f64 {
    plots::window_value(p, &vec![0.0; p.len()], false, &Window { name: "a".into(), from, to }, &[]).unwrap().integral
}

#[test]
fn pl6_normalisation_as_a_view() {
    let set = plots::open_buf(&tb().join("Tb_vs_T.buf"), None).unwrap();
    let before = set.clone();
    // The largest count 1.
    let n = plots::normalise(&set, &Normalisation::Maximum { to: 1.0 }).unwrap();
    for (x, m) in n.iter().zip(&set.members) {
        let p = x.pattern.as_ref().unwrap();
        assert_eq!(x.factor, Some(1.0 / m.pattern.max_count().unwrap()));
        assert!((p.y.iter().copied().fold(f64::MIN, f64::max) - 1.0).abs() < 1e-15);
        assert!(p.sigma.iter().zip(&m.pattern.sigma).all(|(s, s0)| *s == s0 * x.factor.unwrap()));
    }
    // The same area over 16–18° (the background before the first peak; Josu, L169): every member's area there 1000.
    let v = set.normalised(&Normalisation::Area { from: 16.0, until: 18.0, to: 1000.0 }).unwrap();
    for m in &v.members {
        assert!((area(&m.pattern, 16.0, 18.0) - 1000.0).abs() < 1e-9, "{}", m.label);
    }
    // To the 30 K member's area over 16–18°: its own factor 1, every other member's area its.
    let k30 = set.members.iter().position(|m| m.variable.as_ref().is_some_and(|v| v.value == 30.0)).unwrap();
    let r = plots::normalise(&set, &Normalisation::Reference { member: k30, from: 16.0, until: 18.0 }).unwrap();
    assert_eq!(r[k30].factor, Some(1.0));
    let a30 = area(&set.members[k30].pattern, 16.0, 18.0);
    for x in &r {
        assert!((area(x.pattern.as_ref().unwrap(), 16.0, 18.0) - a30).abs() < 1e-9 * a30);
    }
    // The magnetic peak at 19.1° above its background: each member's factor 1 over that integral.
    let w = Window { name: "m1".into(), from: 18.5, to: 19.8 };
    let peak = plots::normalise(&set, &Normalisation::Peak { window: w.clone(), background: WindowBackground::Line, to: 1.0 }).unwrap();
    let i4 = plots::window_value(&set.members[0].pattern, &[], true, &w, &[]).unwrap().integral;
    assert_eq!(peak[0].factor, Some(1.0 / i4));
    // At 30 K the peak is gone: its integral is not positive, and the member is said, not normalised.
    let gone: Vec<&plots::Normalised> = peak.iter().filter(|x| x.factor.is_none()).collect();
    println!("peak normalisation refused for: {:?}", gone.iter().map(|x| (&set.members[x.member].label, &x.why_not)).collect::<Vec<_>>());
    assert_eq!(gone.iter().map(|x| set.members[x.member].label.as_str()).collect::<Vec<_>>(), ["tb_t20.dat"]);
    assert!(gone[0].why_not.as_deref().unwrap().starts_with("the peak's integral is -"));
    assert!(set.normalised(&Normalisation::Peak { window: w, background: WindowBackground::Line, to: 1.0 }).unwrap_err().contains("tb_t20.dat"));
    // The largest count in a range.
    let mi = plots::normalise(&set, &Normalisation::MaximumIn { from: 26.0, until: 27.4, to: 100.0 }).unwrap();
    let m0 = set.members[0].pattern.x.iter().zip(&set.members[0].pattern.y).filter(|(x, _)| **x >= 26.0 && **x <= 27.4).map(|(_, y)| *y).fold(f64::MIN, f64::max);
    assert_eq!(mi[0].factor, Some(100.0 / m0));
    // No monitor in these files: each member says so, and the set is not normalised by it.
    let mon = plots::normalise(&set, &Normalisation::Monitor { to: 1e5 }).unwrap();
    assert!(mon.iter().all(|x| x.factor.is_none() && x.why_not.as_deref() == Some("the file states no monitor")));
    assert!(set.normalised(&Normalisation::Monitor { to: 1e5 }).unwrap_err().contains("tb_t02.dat: the file states no monitor"));
    // D20's files state their monitor: y · to / monitor.
    let d = plots::open(&d20(), &[("la0120.dat".into(), None), ("la0150.dat".into(), None)], None, None).unwrap();
    let dm = d.normalised(&Normalisation::Monitor { to: 1e5 }).unwrap();
    for (a, b) in dm.members.iter().zip(&d.members) {
        let f = 1e5 / b.facts.monitor.unwrap();
        assert!(a.pattern.y.iter().zip(&b.pattern.y).all(|(y, y0)| *y == y0 * f));
    }
    println!("D20 monitors {:?}, counting times {:?}", d.members.iter().map(|m| m.facts.monitor).collect::<Vec<_>>(), d.members.iter().map(|m| m.facts.counting_time).collect::<Vec<_>>());
    // A view: the set itself is as it was.
    assert_eq!(set, before);
}

/// A UXD of two ranges: 10.0–10.9° at 1 s per point, then 10.5–11.4° at 4 s per point, overlapping on 10.5–10.9°;
/// the second as counts or as counts per second.
fn uxd(cps: bool) -> String {
    let a: Vec<f64> = (0..10).map(|i| 100.0 + i as f64).collect();
    let b: Vec<f64> = (0..10).map(|i| 400.0 + 4.0 * i as f64).collect();
    let b_text: Vec<String> = b.iter().map(|v| if cps { format!("{}", v / 4.0) } else { format!("{v}") }).collect();
    format!(
        "_FILEVERSION=1\n_SAMPLE='test'\n_WL1=1.5406\n\n; range 1\n_STEPTIME=1\n_STEPSIZE=0.1\n_2THETA=10.0\n_COUNTS\n{}\n\n; range 2\n_STEPTIME=4\n_STEPSIZE=0.1\n_2THETA=10.5\n{}\n{}\n",
        a.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(" "),
        if cps { "_CPS" } else { "_COUNTS" },
        b_text.join(" ")
    )
}

#[test]
fn pl6_constant_counting_time() {
    let d = dir("constant");
    for cps in [false, true] {
        let f = d.join(if cps { "two_cps.uxd" } else { "two.uxd" });
        std::fs::write(&f, uxd(cps)).unwrap();
        let ds = formats::read_path(&f, None).unwrap();
        assert_eq!((ds.format, ds.scans.len(), ds.scans[0].info.counting_time, ds.scans[1].info.counting_time), ("socabim_uxd", 2, Some(1.0), Some(4.0)));
        let s = drs::data_ops::constant_time(&ds.scans, 2.0).unwrap();
        // 10.0–10.9 from the first range (100 … 109 counts in 1 s, ×2); 11.0–11.4 from the second (counts in 4 s, ×½).
        assert_eq!(s.pattern.len(), 15);
        for i in 0..10 {
            assert!((s.pattern.x[i] - (10.0 + 0.1 * i as f64)).abs() < 1e-12 && s.pattern.y[i] == 2.0 * (100.0 + i as f64), "{cps} {i}");
            assert!((s.pattern.sigma[i] - 2.0 * (100.0 + i as f64).sqrt()).abs() < 1e-12);
        }
        for i in 10..15 {
            let counts = 400.0 + 4.0 * (i - 5) as f64;
            assert!((s.pattern.y[i] - counts / 2.0).abs() < 1e-12 && (s.pattern.sigma[i] - counts.sqrt() / 2.0).abs() < 1e-12, "{cps} {i}: {} {}", s.pattern.y[i], s.pattern.sigma[i]);
        }
        assert_eq!(s.info.counting_time, Some(2.0));
    }
    // XRDML: two scans of 1 s and 2 s per point.
    let x = |a: f64, b: f64, t: f64, counts: &str| format!("<scan appendNumber=\"0\" mode=\"Continuous\" scanAxis=\"Gonio\" status=\"Completed\"><dataPoints><positions axis=\"2Theta\" unit=\"deg\" listType=\"StartEnd\"><startPosition>{a}</startPosition><endPosition>{b}</endPosition></positions><commonCountingTime unit=\"seconds\">{t}</commonCountingTime><counts unit=\"counts\">{counts}</counts></dataPoints></scan>");
    let text = format!("<?xml version=\"1.0\"?><xrdMeasurements><xrdMeasurement measurementType=\"Scan\" status=\"Completed\">{}{}</xrdMeasurement></xrdMeasurements>", x(20.0, 20.4, 1.0, "10 20 30 40 50"), x(20.5, 20.9, 2.0, "60 80 100 120 140"));
    let f = d.join("two.xrdml");
    std::fs::write(&f, text).unwrap();
    let ds = formats::read_path(&f, None).unwrap();
    let s = drs::data_ops::constant_time(&ds.scans, 1.0).unwrap();
    assert_eq!(s.pattern.y, vec![10.0, 20.0, 30.0, 40.0, 50.0, 30.0, 40.0, 50.0, 60.0, 70.0]);
    // A scan without a counting time is refused, not guessed.
    let mut no = ds.scans.clone();
    no[1].info.counting_time = None;
    assert!(drs::data_ops::constant_time(&no, 1.0).unwrap_err().contains("states no counting time"));
    // The command line: written as xys and read back.
    let out = d.join("two.xys");
    let o = drs(&["data", "constant-time", d.join("two.uxd").to_str().unwrap(), "--time", "2", "--out", out.to_str().unwrap()]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let back = formats::read_path(&out, None).unwrap().scans.remove(0);
    let direct = drs::data_ops::constant_time(&formats::read_path(&d.join("two.uxd"), None).unwrap().scans, 2.0).unwrap();
    assert_eq!(back.pattern, direct.pattern);
}

#[test]
fn pl7_csv_and_the_writers_in_d_and_q() {
    let set = plots::open_buf(&tb().join("Tb_vs_T.buf"), None).unwrap();
    // In d: the files state no λ, so it is given (G41's 2.426 Å, WinPLOTR's demo 4); without one, refused by name.
    assert!(plots::converted(&set, Axis::D, None).unwrap_err().contains("tb_t02.dat: the file states no wavelength"));
    let d = plots::converted(&set, Axis::D, Some(2.426)).unwrap();
    let m = &d.members[0];
    assert_eq!((m.facts.axis, m.facts.wavelength, m.pattern.len()), (Axis::D, Some(2.426), 301));
    let two_theta = 46.0_f64;
    let d_expect = 2.426 / (2.0 * (two_theta / 2.0).to_radians().sin());
    assert!((m.pattern.x[0] - d_expect).abs() < 1e-12, "the smallest d is at 46°: {} {d_expect}", m.pattern.x[0]);
    // Wide on the common grid in d: the header names d, every member's file and format.
    let grid = plots::common_grid(&d, plots::GridRange::Overlap, plots::GridStep::Finest).unwrap();
    let wide = plots::csv(&d, &grid);
    assert!(wide.lines().any(|l| l.starts_with("d_angstrom,y1,sigma1")));
    assert!(wide.lines().any(|l| l.starts_with("# y1: tb_t02.dat (scan 1 of 1, free_steps), T = 4 K (the title's \"T.=.4.K\")")), "{}", &wide[..400]);
    // Long: one row per member and point, each the member's own.
    let long = plots::csv_long(&set);
    let rows: Vec<Vec<&str>> = long.lines().filter(|l| !l.starts_with('#')).map(|l| l.split(',').collect()).collect();
    assert_eq!(rows[0], ["member", "2theta_deg", "y", "sigma"]);
    assert_eq!(rows.len(), 1 + 20 * 301);
    for r in &rows[1..] {
        let k: usize = r[0].parse::<usize>().unwrap() - 1;
        let i = set.members[k].pattern.x.iter().position(|x| x.to_string() == r[1]).unwrap();
        assert_eq!((r[2].parse::<f64>().unwrap(), r[3].parse::<f64>().unwrap()), (set.members[k].pattern.y[i], set.members[k].pattern.sigma[i]));
    }
    // A single scan written in d (xye) and in Q (xye, qye), read back on its axis with its numbers; xys refuses d.
    let dd = dir("writers");
    let scan = |set: &plots::PlotSet, axis: Axis| {
        let mut s = formats::Scan::new("tb_t02", axis, set.members[0].pattern.clone());
        s.wavelength = Some(formats::Wavelength::single(2.426));
        s
    };
    let q = plots::converted(&set, Axis::Q, Some(2.426)).unwrap();
    for (s, id, name) in [(scan(&d, Axis::D), "xye", "d.xye"), (scan(&q, Axis::Q), "xye", "q.xye"), (scan(&q, Axis::Q), "qye", "q.qye")] {
        let f = dd.join(name);
        std::fs::write(&f, formats::write(&s, id).unwrap()).unwrap();
        let back = formats::read_path(&f, None).unwrap().scans.remove(0);
        assert_eq!((back.axis, &back.pattern), (s.axis, &s.pattern), "{name}");
    }
    assert!(formats::write(&scan(&d, Axis::D), "xys").unwrap_err().contains("xys: 2θ"));
    // The command line: long, in d.
    let o = drs(&["plots", "export", tb().join("Tb_vs_T.buf").to_str().unwrap(), "--csv", "--long", "--x", "d", "--wavelength", "2.426"]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    assert_eq!(String::from_utf8_lossy(&o.stdout), plots::csv_long(&d));
}

#[test]
fn pl21_the_temperatures_file() {
    // FullProf's own example (news of 25 September 2013): data_23.dat twice, with 1.5 and 5.0; "....." in between.
    let news = "data_23.dat   TSAMP 1.5\ndata_24.dat   TSAMP 2.0\ndata_23.dat   TSAMP 5.0\n.....\ndata_120.dat  TSAMP 260.0\n";
    let t = plots::read_temperatures(news);
    assert_eq!(t.entries.len(), 4);
    assert_eq!(t.conflicts, vec![("data_23.dat".to_string(), vec![1.5, 5.0])]);
    assert_eq!(t.unread, vec![(4, ".....".to_string())]);
    assert_eq!((t.value_of("data_24.dat"), t.value_of("data_23.dat"), t.value_of("data_120.dat"), t.value_of("data_24")), (Some(2.0), None, Some(260.0), None));
    // On three D20 files: the name must coincide with its extension and its case.
    let text = "la0120.dat TSAMP 70.5\nla0121 TSAMP 71.5\nLA0122.DAT tsamp 72.5\n";
    let t = plots::read_temperatures(text);
    let mut set = plots::open(&d20(), &[("la0120.dat".into(), None), ("la0121.dat".into(), None), ("la0122.dat".into(), None)], None, None).unwrap();
    set.set_variable(set.members.iter().position(|m| m.file == "la0121.dat").unwrap(), 99.0, "typed").unwrap();
    assert_eq!(set.variable_from_file("temps.txt", &t, false), 1);
    let m = set.members.iter().find(|m| m.file == "la0120.dat").unwrap();
    assert_eq!(m.variable.as_ref().map(|v| (v.value, v.found.clone())), Some((70.5, Found::TemperatureFile("temps.txt".into()))));
    assert_eq!(m.variable.as_ref().unwrap().found.words(), "the temperatures file temps.txt");
    assert!(m.candidates.iter().any(|c| c.found == Found::Header), "the header's value kept among the candidates");
    assert_eq!(set.members.iter().find(|m| m.file == "la0121.dat").unwrap().variable.as_ref().map(|v| v.value), Some(99.0), "no line for la0121.dat");
    assert_eq!(set.members.iter().find(|m| m.file == "la0122.dat").unwrap().variable.as_ref().map(|v| v.found.clone()), Some(Found::Header), "LA0122.DAT is another name");
    // The command line.
    let d = dir("temps");
    std::fs::write(d.join("temps.txt"), text).unwrap();
    let o = drs(&["plots", "open", d20().join("la0120.dat").to_str().unwrap(), d20().join("la0121.dat").to_str().unwrap(), "--temperatures", d.join("temps.txt").to_str().unwrap(), "--json"]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let s: plots::PlotSet = serde_json::from_slice(&o.stdout).unwrap();
    assert_eq!(s.members.iter().filter(|m| matches!(m.variable.as_ref().map(|v| &v.found), Some(Found::TemperatureFile(_)))).count(), 1);
    assert!(String::from_utf8_lossy(&o.stderr).contains("1 of 2 members took their temperature from it"));
}

#[test]
fn pl22_a_value_below_two_sigma() {
    let set = plots::open_buf(&tb().join("Tb_vs_T.buf"), None).unwrap();
    let ws = [("m1", 18.5, 19.8), ("m2", 26.0, 27.4), ("m3", 31.9, 33.3), ("m4", 41.5, 42.9)].map(|(n, a, b)| Window { name: n.into(), from: a, to: b });
    let w = plots::windows(&set, &ws, &WindowBackground::Auto { log_lambda: None }, &[]).unwrap();
    for mw in &w.members {
        for v in &mw.windows {
            assert_eq!((v.significant, v.significance), (v.integral > plots::WINDOW_SIGNIFICANCE * v.sigma, 2.0));
        }
    }
    // At 30 K (member 18): 19.1° and 32.6° gone, below 2σ; 26.7° and 42.2° keep the nuclear peaks (L177).
    let at30: Vec<bool> = w.members[18].windows.iter().map(|v| v.significant).collect();
    assert_eq!(at30, [false, true, false, true]);
    assert!(w.members[0].windows.iter().all(|v| v.significant));
}
