//! A team's ESRF ID22 in-situ run (10–11 June 2026, 101 scans; unpublished: W8, W9, L428; the names neutral, W13).
//! - W8, always: the 101 names, in the structure the run's renamed files have (`scan_<seq>_T<temp>C_t<seconds>s_<class>`),
//!   read: `T0900C` a temperature of 900 °C, `t00000s` 0 s elapsed, as the run's index gives them; the background
//!   scan's name gives no temperature; a set opened from files of these names is in °C.
//! - W2/W9, when the run's folder is on this machine (`DRS_ESRF_RENAMED` names it): `drs plots propose` on its 101
//!   files offers the three sets as they were built by hand: 25 by temperature, 100 against time, the 1270 °C hold of
//!   39, each rebinned to 0.01° on opening; the background scan left out and named. The files are never copied into
//!   this repository: without the folder the test says so and checks the names only.

use drs::plots::{self, OpenedFile};
use std::path::PathBuf;

/// The 101 names, neutral (W13), with the nominal temperature (°C) and the elapsed seconds the run's index gives, and
/// the exposure class it states. Names only: the data are not in this repository.
const ESRF: [(&str, Option<f64>, f64, &str); 101] = [
    ("scan_001_T0900C_t00000s_long.xye", Some(900.0), 0.0, "long"),
    ("scan_002_T0925C_t00810s_long.xye", Some(925.0), 810.0, "long"),
    ("scan_003_T0950C_t01622s_long.xye", Some(950.0), 1622.0, "long"),
    ("scan_004_T0975C_t02432s_long.xye", Some(975.0), 2432.0, "long"),
    ("scan_005_T1000C_t03242s_long.xye", Some(1000.0), 3242.0, "long"),
    ("scan_006_T1025C_t04050s_long.xye", Some(1025.0), 4050.0, "long"),
    ("scan_007_T1050C_t04862s_long.xye", Some(1050.0), 4862.0, "long"),
    ("scan_008_T1075C_t05672s_long.xye", Some(1075.0), 5672.0, "long"),
    ("scan_009_T1100C_t06482s_long.xye", Some(1100.0), 6482.0, "long"),
    ("scan_010_T1125C_t07296s_long.xye", Some(1125.0), 7296.0, "long"),
    ("scan_011_T1150C_t08102s_long.xye", Some(1150.0), 8102.0, "long"),
    ("scan_012_T1175C_t08910s_long.xye", Some(1175.0), 8910.0, "long"),
    ("scan_013_T1200C_t09720s_long.xye", Some(1200.0), 9720.0, "long"),
    ("scan_014_T1225C_t10528s_long.xye", Some(1225.0), 10528.0, "long"),
    ("scan_015_T1225C_t11156s_long.xye", Some(1225.0), 11156.0, "long"),
    ("scan_016_T1240C_t11778s_long.xye", Some(1240.0), 11778.0, "long"),
    ("scan_017_T1240C_t12422s_long.xye", Some(1240.0), 12422.0, "long"),
    ("scan_018_T1240C_t12672s_short.xye", Some(1240.0), 12672.0, "short"),
    ("scan_019_T1240C_t12934s_short.xye", Some(1240.0), 12934.0, "short"),
    ("scan_020_T1240C_t13196s_short.xye", Some(1240.0), 13196.0, "short"),
    ("scan_021_T1240C_t13458s_short.xye", Some(1240.0), 13458.0, "short"),
    ("scan_022_T1240C_t13732s_short.xye", Some(1240.0), 13732.0, "short"),
    ("scan_023_T1245C_t13992s_short.xye", Some(1245.0), 13992.0, "short"),
    ("scan_024_T1245C_t14252s_short.xye", Some(1245.0), 14252.0, "short"),
    ("scan_025_T1245C_t14524s_short.xye", Some(1245.0), 14524.0, "short"),
    ("scan_026_T1250C_t14786s_short.xye", Some(1250.0), 14786.0, "short"),
    ("scan_027_T1250C_t15046s_short.xye", Some(1250.0), 15046.0, "short"),
    ("scan_028_T1250C_t15310s_short.xye", Some(1250.0), 15310.0, "short"),
    ("scan_029_T1250C_t15570s_short.xye", Some(1250.0), 15570.0, "short"),
    ("scan_030_T1250C_t15832s_short.xye", Some(1250.0), 15832.0, "short"),
    ("scan_031_T1266C_t16268s_short.xye", Some(1266.0), 16268.0, "short"),
    ("scan_032_T1270C_t16530s_short.xye", Some(1270.0), 16530.0, "short"),
    ("scan_033_T1270C_t16792s_short.xye", Some(1270.0), 16792.0, "short"),
    ("scan_034_T1270C_t17054s_short.xye", Some(1270.0), 17054.0, "short"),
    ("scan_035_T1270C_t17314s_short.xye", Some(1270.0), 17314.0, "short"),
    ("scan_036_T1270C_t17578s_short.xye", Some(1270.0), 17578.0, "short"),
    ("scan_037_T1270C_t17838s_short.xye", Some(1270.0), 17838.0, "short"),
    ("scan_038_T1270C_t18100s_short.xye", Some(1270.0), 18100.0, "short"),
    ("scan_039_T1270C_t18362s_short.xye", Some(1270.0), 18362.0, "short"),
    ("scan_040_T1270C_t18624s_short.xye", Some(1270.0), 18624.0, "short"),
    ("scan_041_T1270C_t18884s_short.xye", Some(1270.0), 18884.0, "short"),
    ("scan_042_T1270C_t19146s_short.xye", Some(1270.0), 19146.0, "short"),
    ("scan_043_T1270C_t19408s_short.xye", Some(1270.0), 19408.0, "short"),
    ("scan_044_T1270C_t19686s_short.xye", Some(1270.0), 19686.0, "short"),
    ("scan_045_T1270C_t19948s_short.xye", Some(1270.0), 19948.0, "short"),
    ("scan_046_T1270C_t20210s_short.xye", Some(1270.0), 20210.0, "short"),
    ("scan_047_T1270C_t20472s_short.xye", Some(1270.0), 20472.0, "short"),
    ("scan_048_T1270C_t20734s_short.xye", Some(1270.0), 20734.0, "short"),
    ("scan_049_T1270C_t20996s_short.xye", Some(1270.0), 20996.0, "short"),
    ("scan_050_T1270C_t21256s_short.xye", Some(1270.0), 21256.0, "short"),
    ("scan_051_T1270C_t21518s_short.xye", Some(1270.0), 21518.0, "short"),
    ("scan_052_T1270C_t21778s_short.xye", Some(1270.0), 21778.0, "short"),
    ("scan_053_T1270C_t22040s_short.xye", Some(1270.0), 22040.0, "short"),
    ("scan_054_T1270C_t22300s_short.xye", Some(1270.0), 22300.0, "short"),
    ("scan_055_T1270C_t22562s_short.xye", Some(1270.0), 22562.0, "short"),
    ("scan_056_T1270C_t22824s_short.xye", Some(1270.0), 22824.0, "short"),
    ("scan_057_T1270C_t23086s_short.xye", Some(1270.0), 23086.0, "short"),
    ("scan_058_T1270C_t23346s_short.xye", Some(1270.0), 23346.0, "short"),
    ("scan_059_T1270C_t23610s_short.xye", Some(1270.0), 23610.0, "short"),
    ("scan_060_T1270C_t23870s_short.xye", Some(1270.0), 23870.0, "short"),
    ("scan_061_T1270C_t24132s_short.xye", Some(1270.0), 24132.0, "short"),
    ("scan_062_T1270C_t24394s_short.xye", Some(1270.0), 24394.0, "short"),
    ("scan_063_T1270C_t24656s_short.xye", Some(1270.0), 24656.0, "short"),
    ("scan_064_T1270C_t24918s_short.xye", Some(1270.0), 24918.0, "short"),
    ("scan_065_T1270C_t25180s_short.xye", Some(1270.0), 25180.0, "short"),
    ("scan_066_T1270C_t25442s_short.xye", Some(1270.0), 25442.0, "short"),
    ("scan_067_T1270C_t25704s_short.xye", Some(1270.0), 25704.0, "short"),
    ("scan_068_T1270C_t25964s_short.xye", Some(1270.0), 25964.0, "short"),
    ("scan_069_T1270C_t26228s_short.xye", Some(1270.0), 26228.0, "short"),
    ("scan_070_T1270C_t26516s_short.xye", Some(1270.0), 26516.0, "short"),
    ("scan_071_T1272C_t26806s_short.xye", Some(1272.0), 26806.0, "short"),
    ("scan_072_T1275C_t27066s_short.xye", Some(1275.0), 27066.0, "short"),
    ("scan_073_T1275C_t27328s_short.xye", Some(1275.0), 27328.0, "short"),
    ("scan_074_T1275C_t27590s_short.xye", Some(1275.0), 27590.0, "short"),
    ("scan_075_T1275C_t27850s_short.xye", Some(1275.0), 27850.0, "short"),
    ("scan_076_T1275C_t28110s_short.xye", Some(1275.0), 28110.0, "short"),
    ("scan_077_T1275C_t28372s_short.xye", Some(1275.0), 28372.0, "short"),
    ("scan_078_T1275C_t28634s_short.xye", Some(1275.0), 28634.0, "short"),
    ("scan_079_T1275C_t28904s_short.xye", Some(1275.0), 28904.0, "short"),
    ("scan_080_T1280C_t29164s_short.xye", Some(1280.0), 29164.0, "short"),
    ("scan_081_T1280C_t29426s_short.xye", Some(1280.0), 29426.0, "short"),
    ("scan_082_T1280C_t29692s_short.xye", Some(1280.0), 29692.0, "short"),
    ("scan_083_T1280C_t29950s_short.xye", Some(1280.0), 29950.0, "short"),
    ("scan_084_T1290C_t30726s_long.xye", Some(1290.0), 30726.0, "long"),
    ("scan_085_T1290C_t31630s_short.xye", Some(1290.0), 31630.0, "short"),
    ("scan_086_T1299C_t32088s_short.xye", Some(1299.0), 32088.0, "short"),
    ("scan_087_T1300C_t32348s_short.xye", Some(1300.0), 32348.0, "short"),
    ("scan_088_T1300C_t32610s_short.xye", Some(1300.0), 32610.0, "short"),
    ("scan_089_T1300C_t32872s_short.xye", Some(1300.0), 32872.0, "short"),
    ("scan_090_T1300C_t33132s_short.xye", Some(1300.0), 33132.0, "short"),
    ("scan_091_T1300C_t33394s_short.xye", Some(1300.0), 33394.0, "short"),
    ("scan_092_T1300C_t33658s_short.xye", Some(1300.0), 33658.0, "short"),
    ("scan_093_T1300C_t33920s_short.xye", Some(1300.0), 33920.0, "short"),
    ("scan_094_T1300C_t34182s_short.xye", Some(1300.0), 34182.0, "short"),
    ("scan_095_T1300C_t34444s_short.xye", Some(1300.0), 34444.0, "short"),
    ("scan_096_T1300C_t34704s_short.xye", Some(1300.0), 34704.0, "short"),
    ("scan_097_T1300C_t34966s_short.xye", Some(1300.0), 34966.0, "short"),
    ("scan_098_T1300C_t35228s_short.xye", Some(1300.0), 35228.0, "short"),
    ("scan_099_T1300C_t35490s_short.xye", Some(1300.0), 35490.0, "short"),
    ("scan_100_T1300C_t36128s_long.xye", Some(1300.0), 36128.0, "long"),
    ("scan_101_background_nosample_t38192s_long.xye", None, 38192.0, "long"),
];

#[test]
fn w8_the_esrf_names() {
    for (name, t_c, secs, _) in ESRF {
        let t = plots::temperature_in_name(name);
        match t_c {
            Some(v) => {
                let (value, text) = t.unwrap_or_else(|| panic!("{name}: no temperature"));
                assert_eq!(value, v, "{name}");
                assert_eq!(text, format!("T{:04}C", v as i64), "{name}");
                assert_eq!(plots::Found::FileName(text).unit(), Some("°C"), "{name}");
            }
            None => assert_eq!(t, None, "{name}"),
        }
        let (s, text) = plots::elapsed_in_name(name).unwrap_or_else(|| panic!("{name}: no time"));
        assert_eq!(s, secs, "{name}");
        assert_eq!(text, format!("t{:05}s", secs as i64), "{name}");
    }
    // The originals' names (file_index.csv's original_name) say °C without the T: read alike.
    assert_eq!(plots::temperature_in_name("scan_900C_w0003_b0001-0s.xye"), Some((900.0, "900C".into())));
    assert_eq!(plots::temperature_in_name("sample_1270°C.xye"), Some((1270.0, "1270°C".into())));
    // A C inside a word is not a temperature; a K name is in kelvin as before.
    assert_eq!(plots::temperature_in_name("run_12Cu.xye"), None);
    assert_eq!(plots::temperature_in_name("LaMnO3_300K.dat"), Some((300.0, "300K".into())));
    assert_eq!(plots::Found::FileName("300K".into()).unit(), Some("K"));
    // A set from files of these names is in °C; with one in K among them, in K, the °C converted.
    let xye = |name: &str| OpenedFile { name: name.into(), text: (0..40).map(|i| format!("{} {} {}\n", 1.0 + 0.1 * i as f64, 100.0 + i as f64, 10.0)).collect(), ..Default::default() };
    let three: Vec<OpenedFile> = ESRF[..3].iter().map(|e| xye(e.0)).collect();
    let set = plots::open_from(&three, Some("xye"), None).unwrap();
    assert_eq!(set.variable_unit, "°C");
    assert_eq!(set.members.iter().map(|m| m.variable.as_ref().unwrap().value).collect::<Vec<_>>(), [900.0, 925.0, 950.0]);
    assert_eq!(set.members[0].variable.as_ref().unwrap().found.words(), "the file name's \"T0900C\"");
    let mut mixed = three.clone();
    mixed.push(xye("other_1200K.xye"));
    let k = plots::open_from(&mixed, Some("xye"), None).unwrap();
    assert_eq!(k.variable_unit, "K");
    assert_eq!(k.members.iter().map(|m| m.variable.as_ref().unwrap().value).collect::<Vec<_>>(), [900.0 + 273.15, 925.0 + 273.15, 1200.0, 950.0 + 273.15], "sorted by T");
}

/// The folder of the 101 files, when this machine has it (`DRS_ESRF_RENAMED`): its `.xye` files in name order, which is
/// the order of acquisition.
fn esrf_folder() -> Option<Vec<String>> {
    let dir = PathBuf::from(std::env::var("DRS_ESRF_RENAMED").ok()?);
    let mut files: Vec<String> = std::fs::read_dir(&dir).ok()?.filter_map(|e| e.ok()).map(|e| e.path()).filter(|p| p.extension().is_some_and(|x| x == "xye")).map(|p| p.to_string_lossy().into_owned()).collect();
    files.sort();
    (files.len() == ESRF.len()).then_some(files)
}

/// The scans (their numbers in the names) of the sets built by hand (by-temperature-map.drv, all-scans.drv and
/// hold-1270.drv of the run's project).
const BY_TEMPERATURE: [usize; 25] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 17, 25, 30, 31, 70, 71, 79, 83, 84, 86, 100];

/// A scan's number in its name: the token before its temperature (`T0900C`) or before `background`, in the neutral
/// names and in the run's own alike.
fn number(file: &str) -> usize {
    let base = file.rsplit('/').next().unwrap_or(file);
    let parts: Vec<&str> = base.split('_').collect();
    let at = parts.iter().position(|p| *p == "background" || (p.starts_with('T') && p.ends_with('C') && p[1..p.len() - 1].chars().all(|c| c.is_ascii_digit()))).unwrap();
    parts[at - 1].parse().unwrap()
}

#[test]
fn w2_the_esrf_run_offers_l021s_three_sets() {
    let Some(files) = esrf_folder() else {
        eprintln!("the ESRF folder is not on this machine (set DRS_ESRF_RENAMED): its names are checked by w8_the_esrf_names only");
        return;
    };
    let o = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(["plots", "propose"]).args(&files).arg("--json").output().unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let p: plots::Proposals = serde_json::from_slice(&o.stdout).unwrap();
    let kinds: Vec<(&str, usize)> = p.proposals.iter().map(|x| (x.kind.as_str(), x.set.members.len())).collect();
    assert_eq!(kinds, [("by_temperature", 25), ("against_time", 100), ("hold", 39)], "{:?}", p.proposals.iter().map(|x| &x.title).collect::<Vec<_>>());
    let scans = |k: usize| p.proposals[k].set.members.iter().map(|m| number(&m.file)).collect::<Vec<_>>();
    assert_eq!(scans(0), BY_TEMPERATURE);
    assert_eq!(scans(1), (1..=100).collect::<Vec<_>>());
    assert_eq!(scans(2), (32..=70).collect::<Vec<_>>());
    // In °C; the hold against the minutes since its first scan, the run against the hours since scan 001, as built.
    assert_eq!(p.proposals[0].set.variable_unit, "°C");
    assert_eq!((p.proposals[1].set.variable_unit.as_str(), p.proposals[2].set.variable_unit.as_str()), ("h", "min"));
    let last = |k: usize| p.proposals[k].set.members.last().unwrap().value.as_ref().unwrap().0;
    assert!((last(1) - 10.0356).abs() < 1e-4 && (last(2) - 166.43).abs() < 1e-2, "{} {}", last(1), last(2));
    // Three exposure classes: the long scans with the background, the short ones, and scan 085 alone.
    assert_eq!(p.classes.iter().map(|c| c.scans).collect::<Vec<_>>(), [20, 80, 1]);
    assert_eq!(p.proposals[2].class, Some(1));
    for (name, _, _, class) in ESRF.iter().filter(|e| e.0.contains("scan_085_") || e.0.contains("scan_001_") || e.0.contains("scan_032_")) {
        println!("{name}: {class} in file_index.csv");
    }
    // The background named, out of every set; each set rebinned to 0.01° on opening, its title saying so.
    let bg = files[100].clone();
    assert_eq!(p.left_out, vec![(bg.clone(), "no temperature, where the run's other scans give one: in no set".to_string())]);
    assert!(p.proposals[1].left_out.contains(&(bg, "no temperature, where the run's other scans give one: left out of the run".to_string())));
    for x in &p.proposals {
        assert_eq!(x.set.rebin, Some(0.01), "{}", x.title);
        assert!(x.title.ends_with(", rebinned to 0.01°"), "{}", x.title);
        println!("{}", x.title);
    }
    // Opened as it is: the hold's 39 members rebinned, 26 992 points each to about 2 700.
    let hold = p.proposals[2].set.open(std::path::Path::new(""), None, &|_| Vec::new()).unwrap();
    assert_eq!(hold.members.len(), 39);
    assert!(hold.members.iter().all(|m| (2690..=2705).contains(&m.pattern.len())), "{}", hold.members[0].pattern.len());
    assert_eq!(hold.members[0].others.iter().map(|o| (o.name.as_str(), o.value)).collect::<Vec<_>>(), [("T", 1270.0)]);
}
