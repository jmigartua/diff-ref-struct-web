//! The web stream's W2 (I10 step 3; L423, L426, L428): the series a run's files offer (`plots::propose`), and W3's `t`.
//! - Tb₂Pd₂Sn's buffer (20 files, the temperatures from the titles, no time, σ = √y: one exposure class): one set, by
//!   temperature, its 19 scans 4 to 30 K; BT (tb_t01, no temperature) left out with why; no set against time.
//! - A run written here (xys files whose names give the elapsed seconds, a TEMP line in each header): a ramp, a hold at
//!   1270 K and a cooling, and two scans of a shorter exposure (σ ten times larger): by temperature (the long exposure's
//!   last scan at each temperature), the whole run against time (both classes) and the long class's hold; a file the
//!   registry cannot read named.
//! - W3's `t`: from the elapsed seconds in the names when every name gives them, else, in `propose` only, from the
//!   modification times the host gives, unless two fall within a second (a copy's or a checkout's); never a mix.
//!
//! The ESRF run's 101 files are `plots_esrf`'s.

use drs::plots::{self, OpenedFile};
use std::path::{Path, PathBuf};

fn tb() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/1-WinPLOTR/1_Tb2pd2sn")
}

/// Tb₂Pd₂Sn's buffer as the host gives it: its files with the buffer's comments.
fn tb_files() -> Vec<OpenedFile> {
    let buf = std::fs::read_to_string(tb().join("Tb_vs_T.buf")).unwrap();
    plots::read_buf(&buf, &|n| tb().join(n).is_file()).entries.into_iter().map(|(n, c)| OpenedFile { text: std::fs::read_to_string(tb().join(&n)).unwrap(), name: n, comment: c, ..Default::default() }).collect()
}

#[test]
fn w2_the_tb_buffer_offers_one_set_by_temperature() {
    let p = plots::propose(&tb_files());
    assert_eq!(p.proposals.len(), 1, "{:?}", p.proposals.iter().map(|x| &x.title).collect::<Vec<_>>());
    let by_t = &p.proposals[0];
    assert_eq!(by_t.kind, "by_temperature");
    assert_eq!(by_t.title, "19 scans by temperature, 4 to 30 K");
    assert_eq!(by_t.class, None, "one exposure class");
    assert_eq!(by_t.set.members.len(), 19);
    assert_eq!((by_t.set.variable_name.as_str(), by_t.set.variable_unit.as_str()), ("T", "K"));
    assert!(by_t.set.members.iter().all(|m| m.value.is_none()), "a temperature the titles give is not typed: a re-read gives it");
    assert_eq!(by_t.left_out, vec![("tb_t01.dat".to_string(), "no temperature: neither its header, a buffer's comment, its title nor its name gives one".to_string())]);
    assert_eq!(p.classes.len(), 1);
    assert_eq!(p.left_out, vec![("tb_t01.dat".to_string(), "no temperature and no time: in no set".to_string())]);
    // The set file opens as it is: the 19 members, by temperature, as the buffer gives them.
    let set = by_t.set.open_from(&tb_files(), None, &|_| Vec::new()).unwrap();
    assert_eq!(set.members.len(), 19);
    let opened = plots::open_buf(&tb().join("Tb_vs_T.buf"), None).unwrap();
    for m in &set.members {
        let o = opened.members.iter().find(|x| x.file == m.file).unwrap();
        assert_eq!(m.variable.as_ref().map(|v| v.value), o.variable.as_ref().map(|v| v.value), "{}", m.file);
    }
    // Written and read back as text: the host may save it.
    assert_eq!(plots::SetFile::read(&by_t.set.to_text()).unwrap(), by_t.set);
}

/// One xys file of the run written here: its header's TEMP line, its points at σ = `esd`.
fn scan(name: &str, t: f64, esd: f64) -> OpenedFile {
    let mut text = format!("XYDATA a run written for the test\nTEMP {t}\n!\n!\n!\n!\n");
    for i in 0..20 {
        text.push_str(&format!("{} {} {esd}\n", 10.0 + 0.1 * i as f64, 100.0 + i as f64));
    }
    OpenedFile { name: name.into(), text, ..Default::default() }
}

/// The run: five scans of a ramp 300–700 K, four held at 1270 K, two of a cooling (1000, 800 K), each 600 s apart and
/// named by its elapsed seconds, σ ≈ 1; two scans of a short exposure (σ ≈ 10) at 1270 K; a file nothing reads.
fn run() -> Vec<OpenedFile> {
    let mut f = Vec::new();
    let temps = [300.0, 400.0, 500.0, 600.0, 700.0, 1270.0, 1270.2, 1269.9, 1270.1, 1000.0, 800.0];
    for (k, t) in temps.iter().enumerate() {
        f.push(scan(&format!("run_{:05}s.xys", 600 * k), *t, 1.0 + 0.01 * k as f64));
    }
    f.push(scan("run_06900s_short.xys", 1270.0, 10.0));
    f.push(scan("run_07500s_short.xys", 1270.0, 10.2));
    f.push(OpenedFile { name: "notes.txt".into(), text: "the run's notes, no data\n".into(), ..Default::default() });
    f
}

#[test]
fn w2_a_run_offers_its_series() {
    let p = plots::propose(&run());
    let titles: Vec<&str> = p.proposals.iter().map(|x| x.title.as_str()).collect();
    // Two exposure classes by the counts' scale σ²/y, the long first (the smaller scale).
    assert_eq!(p.classes.iter().map(|c| c.scans).collect::<Vec<_>>(), [11, 2]);
    let kinds: Vec<(&str, Option<usize>, usize)> = p.proposals.iter().map(|x| (x.kind.as_str(), x.class, x.set.members.len())).collect();
    assert_eq!(kinds, [("by_temperature", Some(0), 8), ("against_time", None, 13), ("hold", Some(0), 4)], "{titles:?}");
    assert!(titles[0].starts_with("8 scans by temperature, 300 to 1270.1 K (exposure class 1 of 2, σ²/y 0.009 to 0.011)"), "{}", titles[0]);
    assert_eq!(titles[1], "13 scans against time, 125 min from run_00000s.xys (2 exposure classes)");
    assert!(titles[2].starts_with("the hold at 1270 K: 4 scans against time, 30 min from run_03000s.xys (exposure class 1 of 2"), "{}", titles[2]);
    // The time typed with its source, the temperature the member's other variable with its source and unit.
    let t = &p.proposals[1].set;
    assert_eq!((t.variable_name.as_str(), t.variable_unit.as_str()), ("t", "min"));
    let m5 = &t.members[5];
    assert_eq!(m5.file, "run_03000s.xys");
    assert_eq!(m5.value.as_ref().map(|v| v.0), Some(50.0));
    assert_eq!(m5.value.as_ref().map(|v| v.1.words()), Some("typed: the time since run_00000s.xys, from the file name's \"03000s\"".to_string()));
    assert_eq!(m5.variables, vec![("T".to_string(), 1270.0, "the file's header (K)".to_string())]);
    // The hold against the time since its own first scan, as the ESRF hold was built by hand.
    assert_eq!(p.proposals[2].set.members.iter().map(|m| m.value.as_ref().unwrap().0).collect::<Vec<_>>(), [0.0, 10.0, 20.0, 30.0]);
    // By temperature: at each temperature the last scan of the longest exposure, the others left out naming it.
    let by_t = &p.proposals[0];
    assert!(by_t.set.members.iter().any(|m| m.file == "run_04800s.xys"));
    for f in ["run_03000s.xys", "run_06900s_short.xys"] {
        assert!(by_t.left_out.iter().any(|(x, why)| x == f && why.ends_with(" K is run_04800s.xys's: the last scan of the longest exposure there is taken")), "{f}: {:?}", by_t.left_out);
    }
    // The short scans are in the run against time; the notes are named with why.
    assert!(p.proposals[1].set.members.iter().any(|m| m.file == "run_07500s_short.xys"));
    assert_eq!(p.left_out.iter().map(|x| x.0.as_str()).collect::<Vec<_>>(), ["notes.txt"]);
    assert!(p.left_out[0].1.starts_with("not read:"), "{:?}", p.left_out);
    // Each proposal opens as it is, and its rows carry T (W3) on the map against time; nothing to rebin (20 points).
    let readable: Vec<OpenedFile> = run().into_iter().filter(|f| f.name != "notes.txt").collect();
    let set = p.proposals[1].set.open_from(&readable, None, &|_| Vec::new()).unwrap();
    assert_eq!(set.members.len(), 13);
    let grid = plots::common_grid(&set, plots::GridRange::Overlap, plots::GridStep::Finest).unwrap();
    let map = plots::map(&set, &grid, plots::YScale::Linear, None).unwrap();
    assert_eq!(map.rows.iter().map(|r| r.value).collect::<Vec<_>>(), [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 115.0, 125.0]);
    assert_eq!((map.rows[5].others[0].name.as_str(), map.rows[5].others[0].value, map.rows[5].others[0].source.as_str()), ("T", 1270.0, "typed: the file's header (K)"));
    for x in &p.proposals {
        assert_eq!(x.set.rebin, None, "{}", x.title);
        assert_eq!(plots::SetFile::read(&x.set.to_text()).unwrap(), x.set, "{}", x.title);
    }
}

#[test]
fn w2_times_from_the_host_when_the_names_give_none() {
    // The same run without seconds in the names, the host giving each file's modification time.
    let files: Vec<OpenedFile> = run()
        .into_iter()
        .enumerate()
        .map(|(k, f)| OpenedFile { name: f.name.replace("run_", "scan").replace('s', "x"), modified: Some(1_700_000_000.0 + 600.0 * k as f64), ..f })
        .collect();
    let p = plots::propose(&files);
    let t = p.proposals.iter().find(|x| x.kind == "against_time").unwrap();
    assert_eq!(t.set.members.len(), 13);
    assert_eq!(t.set.members[1].value.as_ref().map(|v| (v.0, v.1.words())), Some((10.0, format!("typed: the time since {}, from typed: the file's modification time", files[0].name))));
    // Files checked out or copied together carry one write's time (within `SAME_WRITE`): no time taken from any, no
    // set against time, the rest as before.
    let copied: Vec<OpenedFile> = files.iter().enumerate().map(|(k, f)| OpenedFile { modified: Some(1_761_400_000.0 + 0.01 * k as f64), ..f.clone() }).collect();
    let c = plots::propose(&copied);
    assert!(c.proposals.iter().all(|x| x.kind == "by_temperature"), "{:?}", c.proposals.iter().map(|x| &x.title).collect::<Vec<_>>());
    assert_eq!(c.proposals[0].set.members.len(), 8);
    let tb_on_disk: Vec<OpenedFile> = tb_files().into_iter().map(|f| OpenedFile { modified: Some(1_761_400_000.0), ..f }).collect();
    assert_eq!(plots::propose(&tb_on_disk).proposals.iter().map(|x| x.kind.as_str()).collect::<Vec<_>>(), ["by_temperature"], "a checkout's times");
}

#[test]
fn w3_t_from_the_names() {
    assert_eq!(plots::elapsed_in_name("run_03600s.xys"), Some((3600.0, "03600s".into())));
    // A lower-case t before the seconds says a time, as the ESRF names write it (W8); L426's guess that it was inside
    // a word is undone by those names.
    assert_eq!(plots::elapsed_in_name("hold_t12.5sec.dat"), Some((12.5, "t12.5sec".into())));
    assert_eq!(plots::elapsed_in_name("scan_032_T1270C_t16530s_short.xye"), Some((16530.0, "t16530s".into())));
    assert_eq!(plots::elapsed_in_name("hold_at12.5sec.dat"), None, "at12.5 is inside a word");
    assert_eq!(plots::elapsed_in_name("hold_12.5sec.dat"), Some((12.5, "12.5sec".into())));
    assert_eq!(plots::elapsed_in_name("scan_3600seconds.xye"), None, "seconds is a longer word");
    assert_eq!(plots::elapsed_in_name("tb_t01.dat"), None);
    // A set opened from files whose names give their seconds carries t on each member, since the first.
    let set = plots::open_from(&run()[..11], None, None).unwrap();
    let t_of = |f: &str| set.members.iter().find(|m| m.file == f).unwrap().others.iter().find(|o| o.name == "t").map(|o| (o.value, o.unit.clone(), o.found.words()));
    assert_eq!(t_of("run_00000s.xys"), Some((0.0, "s".into(), "the file name's \"00000s\"".into())));
    assert_eq!(t_of("run_06000s.xys"), Some((6000.0, "s".into(), "the file name's \"06000s\"".into())));
    // A set whose names do not all give seconds carries none; a modification time alone is never taken on opening.
    let mixed: Vec<OpenedFile> = run().into_iter().take(11).enumerate().map(|(k, f)| if k == 3 { OpenedFile { name: "other.xys".into(), modified: Some(1.0), ..f } } else { f }).collect();
    assert!(plots::open_from(&mixed, None, None).unwrap().members.iter().all(|m| m.others.iter().all(|o| o.name != "t")));
    let tb = plots::open_from(&tb_files(), None, None).unwrap();
    assert!(tb.members.iter().all(|m| m.others.is_empty()));
}

#[test]
fn w2_drs_plots_propose() {
    // The buffer's files on disk with their modification times (a checkout's here: no set against time; the text
    // checked does not depend on them).
    let o = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(["plots", "propose"]).arg(tb().join("Tb_vs_T.buf")).output().unwrap();
    let out = String::from_utf8_lossy(&o.stdout);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(out.contains(". by_temperature: 19 scans by temperature, 4 to 30 K\n   left out tb_t01.dat: no temperature: neither its header, a buffer's comment, its title nor its name gives one\n"), "{out}");
    let j = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(["plots", "propose"]).arg(tb().join("Tb_vs_T.buf")).arg("--json").output().unwrap();
    let p: plots::Proposals = serde_json::from_slice(&j.stdout).unwrap();
    assert_eq!(p, plots::propose(&tb_files().into_iter().map(|f| OpenedFile { modified: std::fs::metadata(tb().join(&f.name)).unwrap().modified().ok().map(|t| t.duration_since(std::time::UNIX_EPOCH).unwrap().as_secs_f64()), ..f }).collect::<Vec<_>>()));
}
