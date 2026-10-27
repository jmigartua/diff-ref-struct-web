//! Per-format reading options and an order of preference among the readers (the app's request 51; test
//! `data_format_options`, L400).
//!
//! - **The preference**: the 16 example files that `data_formats::recognition` leaves ambiguous (its 17 file–format
//!   pairs: xys with xye, or with the general table, read their columns) are read by the first of the tied readers in
//!   the preference, and say so; a preference
//!   naming a reader that is not tied changes nothing; a file one reader leads is never moved by a preference.
//! - **The scan**: a GSAS file of three banks read at bank 2 is the second bank of the whole reading (which
//!   `data_formats` holds to GSAS-II's); a Bruker RAW file's range 2 is its 3501 points (GSAS-II's 3501); a scan
//!   beyond the file's is refused.
//! - **The table's layout**: a semicolon table with the decimal comma, x in its second column, y in its fourth and no
//!   σ, read with those options.
//! - **The declarations**: every option of every format has a kind, a help and a default of its kind; an unknown
//!   option, a value of the wrong kind or out of range is refused, naming the format's options.
//! - **The control file**: `reader_options` written, read back the same, applied when the control file is opened
//!   (the wavelength of a file that states none), and refused by the checks when they do not fit the reader.

use drs::formats::options::{OptionKind, OptionValue, Options};
use drs::formats::{self, Probe};
use std::path::{Path, PathBuf};

fn root(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..").join(rel)
}

fn opts(pairs: &[(&str, OptionValue)]) -> Options {
    pairs.iter().map(|(k, v)| (k.to_string(), v.clone())).collect()
}

#[test]
fn data_format_options() {
    let mut fails: Vec<String> = Vec::new();

    // The preference.
    let index = std::fs::read_to_string(root("core/tests/data/gsas2_reads.txt")).unwrap();
    let mut files: Vec<String> = index.lines().filter(|l| !l.starts_with('#')).map(|l| l.split('\t').next().unwrap().to_string()).collect();
    let fp = std::fs::read_to_string(root("core/tests/data/fullprof_reads.txt")).unwrap();
    files.extend(fp.lines().filter(|l| l.starts_with("READ") || l.starts_with("FILE")).map(|l| l.split_whitespace().nth(3).unwrap().to_string()));
    files.sort();
    files.dedup();
    let (mut tied, mut decided, mut leaders) = (0, 0, 0);
    for f in &files {
        let Ok(bytes) = std::fs::read(root(f)) else { continue };
        let p = Probe::new(f, &bytes);
        let plain = formats::identify(&p);
        if plain.candidates.is_empty() {
            continue;
        }
        match plain.format() {
            Some(lead) => {
                leaders += 1;
                // A preference for any other reader moves nothing.
                let others: Vec<&str> = plain.candidates.iter().map(|c| c.id).filter(|id| *id != lead.id).collect();
                let w = formats::identify_with(&p, &others);
                if w.format().map(|c| c.id) != Some(lead.id) || w.by_preference() {
                    fails.push(format!("{f}: {} leads, but the preference {others:?} gave {:?}", lead.id, w.format().map(|c| c.id)));
                }
            }
            None => {
                tied += 1;
                let ids: Vec<&str> = plain.ambiguous().iter().map(|c| c.id).collect();
                for pref in [["xye", "xys"], ["xys", "xye"]] {
                    let w = formats::identify_with(&p, &pref);
                    let want = pref.iter().find(|x| ids.contains(x)).copied();
                    if w.format().map(|c| c.id) != want || !w.by_preference() {
                        fails.push(format!("{f}: tied {ids:?}, preference {pref:?}: {:?}", w.format().map(|c| c.id)));
                    }
                }
                // Read so: the preferred reader's dataset.
                match formats::read_with(&p, None, &Options::new(), &["xye", "xys"]) {
                    Ok(d) if d.format != if ids.contains(&"xye") { "xye" } else { "xys" } => fails.push(format!("{f}: read as {} with xye, xys preferred", d.format)),
                    Err(e) => fails.push(format!("{f}: {e}")),
                    _ => decided += 1,
                }
                // A preference naming no tied reader leaves it ambiguous.
                if formats::identify_with(&p, &["panalytical_xrdml"]).format().is_some() {
                    fails.push(format!("{f}: a preference for a reader not tied decided"));
                }
            }
        }
    }
    println!("{} files: {leaders} led by a reader, {tied} tied, {decided} read by the preference", files.len());
    // `data_formats::recognition` freezes 17 ambiguous (file, format) pairs: 16 files.
    if tied != 16 || decided != 16 {
        fails.push(format!("{tied} tied and {decided} decided, 16 expected (data_formats::recognition)"));
    }

    // The scan.
    let read = |rel: &str, o: &Options| formats::read_path_with(&root(rel), None, o, &[]);
    let gsas = "manuals/gsas2/tutorials/RMCProfile-I/data/sf6_190gsas.dat";
    let all = read(gsas, &Options::new()).unwrap();
    let two = read(gsas, &opts(&[("scan", OptionValue::Integer(2))])).unwrap();
    if all.scans.len() != 3 || two.scans.len() != 1 || two.scans[0] != all.scans[1] {
        fails.push(format!("GSAS bank 2: {} scans of {}", two.scans.len(), all.scans.len()));
    }
    if read(gsas, &opts(&[("scan", OptionValue::Integer(4))])).is_ok() {
        fails.push("GSAS bank 4 of 3 read".into());
    }
    let raw = read("manuals/gsas2/tutorials/StackingFaults-II/data/Keokuk kaolinite.RAW", &opts(&[("scan", OptionValue::Integer(2))])).unwrap();
    if raw.scans.len() != 1 || raw.scans[0].pattern.len() != 3501 {
        fails.push(format!("Bruker range 2: {} scans, {} points (GSAS-II 3501)", raw.scans.len(), raw.scans.first().map_or(0, |s| s.pattern.len())));
    }

    // The table's layout.
    let text = "angle;x2θ;n;counts\n1;10,0;7;100,5\n2;10,5;7;200\n3;11,0;7;150,25\n";
    let o = opts(&[
        ("delimiter", OptionValue::Text(";".into())),
        ("decimal", OptionValue::Text(",".into())),
        ("header_lines", OptionValue::Integer(1)),
        ("x_column", OptionValue::Integer(2)),
        ("y_column", OptionValue::Integer(4)),
        ("sigma_column", OptionValue::Integer(0)),
    ]);
    match formats::read_with(&Probe::new("t.csv", text.as_bytes()), Some("table"), &o, &[]) {
        Ok(d) => {
            let s = &d.scans[0];
            if s.pattern.x != [10.0, 10.5, 11.0] || s.pattern.y != [100.5, 200.0, 150.25] || (s.pattern.sigma[1] - 200f64.sqrt()).abs() > 1e-12 {
                fails.push(format!("the table: {:?}", s.pattern));
            }
        }
        Err(e) => fails.push(format!("the table: {e}")),
    }

    // The declarations.
    for f in formats::FORMATS {
        for o in f.options {
            if o.help.is_empty() || o.name.is_empty() {
                fails.push(format!("{}: option {} without a name or help", f.id, o.id));
            }
            let fits = o.default.is_empty()
                || match o.kind {
                    OptionKind::Integer { min } => o.default.parse::<i64>().is_ok_and(|n| n >= min),
                    OptionKind::Choice(w) => w.contains(&o.default),
                    OptionKind::Number { .. } => o.default.parse::<f64>().is_ok(),
                    OptionKind::Column { .. } => o.default.parse::<i64>().is_ok(),
                };
            if !fits {
                fails.push(format!("{}: option {}'s default {:?} is not of its kind", f.id, o.id, o.default));
            }
        }
    }
    let gsas_f = formats::by_id("gsas").unwrap();
    for (o, why) in [
        (opts(&[("bank", OptionValue::Integer(2))]), "unknown"),
        (opts(&[("scan", OptionValue::Integer(0))]), "out of range"),
        (opts(&[("scan", OptionValue::Text("two".into()))]), "of the wrong kind"),
        (opts(&[("wavelength", OptionValue::Number(-1.0))]), "negative"),
    ] {
        match formats::options::validate(gsas_f, &o) {
            Ok(()) => fails.push(format!("an option {why} accepted: {o:?}")),
            Err(e) => println!("refused ({why}): {e}"),
        }
    }

    // The control file.
    let tmp = std::env::temp_dir().join(format!("drs_format_options_{}", std::process::id()));
    std::fs::create_dir_all(&tmp).unwrap();
    std::fs::copy(root("manuals/gsas2/runs/y2o3_tch/synth.xys"), tmp.join("synth.xys")).unwrap();
    let text = std::fs::read_to_string(root("manuals/gsas2/runs/y2o3_tch/fit.drc")).unwrap();
    let with = text.replacen("reader = \"xys\"", "reader = \"xys\"\nreader_options = { axis = \"2theta\", wavelength = 1.5406 }", 1);
    let f = tmp.join("fit.drc");
    std::fs::write(&f, &with).unwrap();
    match drs::control::open(&f) {
        Ok(o) => {
            let d = o.control.data.as_ref().unwrap();
            let want = opts(&[("axis", OptionValue::Text("2theta".into())), ("wavelength", OptionValue::Number(1.5406))]);
            if d.reader_options != want {
                fails.push(format!("reader_options read: {:?}", d.reader_options));
            }
            if o.scan.as_ref().and_then(|s| s.wavelength).map(|w| w.lambda1) != Some(1.5406) {
                fails.push(format!("the wavelength option not applied: {:?}", o.scan.as_ref().map(|s| s.wavelength)));
            }
            let back = drs::control::Control::from_text_in(&o.control.to_text(), &tmp).unwrap();
            if back.data != o.control.data || !o.control.to_text().contains("reader_options = { axis = \"2theta\", wavelength = 1.5406 }") {
                fails.push("reader_options not written back the same".into());
            }
        }
        Err(e) => fails.push(format!("the control file with reader_options: {e}")),
    }
    for (bad, why) in [("reader_options = { scan = 2 }", "an option xys has not"), ("reader_options = { axis = \"theta\" }", "a word not among the choices"), ("reader_options = 3", "not a table")] {
        let t = text.replacen("reader = \"xys\"", &format!("reader = \"xys\"\n{bad}"), 1);
        match drs::control::Control::from_text_in(&t, &tmp) {
            Ok(_) => fails.push(format!("{bad}: accepted ({why})")),
            Err(d) => println!("refused ({why}): {}", d[0].message),
        }
    }
    let _ = std::fs::remove_dir_all(&tmp);
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}
