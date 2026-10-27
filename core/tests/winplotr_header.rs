//! Plots' PL32 (L243): a WinPLOTR-made header's temperature read with its unit, never as kelvin by default.
//!
//! La₂CoMnO₆'s case E (`manuals/fullprof/runs/modes/E/t7.dat`, "File created by WinPLOTR" from a Bruker `.UXD`) says
//! `! TEMPERATURE 179.99` without a unit; Josu gives it as 179.99 °C = 453.14 K (L243). The test checks it through the
//! three readers that take the file's columns (xys, xye, table), as a plot set's member (its source the file's header)
//! and as a sequential refinement's external parameter; and that a WinPLOTR header without a unit or a Bruker source
//! gives no temperature, but says why.

use drs::formats::{self, winplotr_temperature, Probe};
use drs::plots::{self, Found};
use std::path::{Path, PathBuf};

fn root(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)
}

fn case_e() -> PathBuf {
    root("../manuals/fullprof/runs/modes/E/t7.dat")
}

fn other<'a>(s: &'a formats::Scan, key: &str) -> Option<&'a str> {
    s.info.other.iter().find(|o| o.0 == key).map(|o| o.1.as_str())
}

#[test]
fn case_e_is_179_99_celsius() {
    for reader in ["xys", "xye", "table"] {
        let d = formats::read_path(&case_e(), Some(reader)).unwrap();
        let s = &d.scans[0];
        assert_eq!(s.info.temperature, Some(179.99 + 273.15), "{reader}");
        assert!((s.info.temperature.unwrap() - 453.14).abs() < 1e-9, "{reader}");
        let said = other(s, "temperature as written").unwrap();
        assert!(said.starts_with("179.99 °C") && said.contains("Bruker .UXD") && said.contains("L243"), "{reader}: {said}");
        assert_eq!(s.pattern.len(), 1576, "{reader}");
    }
    // The first word of WinPLOTR's first line ("! File created by WinPLOTR") is not an instrument.
    assert_eq!(formats::read_path(&case_e(), Some("xys")).unwrap().scans[0].info.instrument, None);
}

#[test]
fn units_known_and_unknown() {
    let file = |lines: &[&str]| format!("{}\n 10.0 100.0 10.0\n 10.1 110.0 10.5\n 10.2 120.0 11.0\n 10.3 100.0 10.0\n", lines.join("\n"));
    let t = |text: &str| winplotr_temperature(text);
    let uxd = file(&["! File created by WinPLOTR", "! INITIAL .UXD FILE: C:\\x.UXD", "! TEMPERATURE         179.99"]);
    assert_eq!(t(&uxd).unwrap().kelvin, Some(179.99 + 273.15));
    // A unit on the line decides, whatever the source.
    let k = file(&["! File created by WinPLOTR", "! INITIAL .UXD FILE: C:\\x.UXD", "! TEMPERATURE  300.0 K"]);
    assert_eq!(t(&k).unwrap().kelvin, Some(300.0));
    let c = file(&["! File saved from WinPLOTR-2006", "! TEMPERATURE: 25 °C"]);
    assert_eq!(t(&c).unwrap().kelvin, Some(298.15));
    // No unit and no source that gives one: not read, never kelvin by default (PL32), and the reason kept.
    let bare = file(&["! File created by WinPLOTR", "! INITIAL file: C:\\x.xy", "! TEMPERATURE 179.99"]);
    let w = t(&bare).unwrap();
    assert_eq!(w.kelvin, None);
    assert!(w.said.contains("not taken as kelvin"), "{}", w.said);
    let d = formats::read(&Probe::new("bare.xys", bare.as_bytes()), Some("xys")).unwrap();
    assert_eq!(d.scans[0].info.temperature, None);
    assert!(other(&d.scans[0], "temperature not read").is_some());
    let odd = file(&["! File created by WinPLOTR", "! TEMPERATURE 179.99 F"]);
    assert_eq!(t(&odd).unwrap().kelvin, None);
    // Not WinPLOTR's header, or WinPLOTR's without a temperature line: nothing.
    assert_eq!(t(&file(&["! some program", "! TEMPERATURE 179.99"])), None);
    let srsc = std::fs::read_to_string(root("../manuals/fullprof/runs/modes/A_modes/SrScLRRT.dat")).unwrap();
    assert_eq!(t(&srsc), None);
}

#[test]
fn a_plot_set_and_a_series_take_it() {
    let dir = case_e();
    let base = dir.parent().unwrap();
    let set = plots::open(base, &[("t7.dat".to_string(), None)], Some("xys"), None).unwrap();
    let v = set.members[0].variable.as_ref().unwrap();
    assert!((v.value - 453.14).abs() < 1e-9);
    assert_eq!(v.found, Found::Header);
    assert_eq!(v.found.words(), "the file's header");
    assert_eq!(set.members[0].facts.temperature, Some(179.99 + 273.15));
    // A sequential refinement's external parameter is the temperature its header gives.
    let text: String = std::fs::read(base.join("e.pcr")).unwrap().iter().map(|&b| b as char).collect();
    let control = drs::control::Control::from_pcr(&drs::pcr::read(&text).unwrap(), base, "e", base).unwrap();
    let input = drs::sequential::read(&control, base, Path::new("t7.dat")).unwrap();
    assert_eq!(input.parameter, Some(179.99 + 273.15));
}

#[test]
fn drs_data_read_says_it() {
    let out = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(["data", "read"]).arg(case_e()).args(["--as", "xys"]).output().unwrap();
    let text = String::from_utf8_lossy(&out.stdout);
    assert!(out.status.success());
    assert!(text.contains("# temperature: 453.14\n"), "{}", &text[..400.min(text.len())]);
    assert!(text.contains("# temperature as written: 179.99 °C (WinPLOTR's header, from a Bruker .UXD file"));
}
