//! The model written as text and read back (test `model_text`, the app's request 1, L058; feature `serde`, on by
//! default). Every FullProf run the core imports (the runs of `manuals/fullprof/runs`, the preferred-orientation run,
//! the window experiments, the restraint and Le Bail runs) is imported, written with `Model::to_text`, read with
//! `Model::from_text` and written again: the two texts are identical, the calculated patterns identical to the bit at
//! the data's points, and so are the lists of editable values (`Model::fields`).
#![cfg(feature = "serde")]

use drs::powder::Model;
use std::path::Path;

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display())).iter().map(|&b| b as char).collect()
}

#[test]
fn model_text_round_trip() {
    let base = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs");
    let mut pcrs = Vec::new();
    for e in std::fs::read_dir(&base).unwrap().flatten() {
        let d = e.path();
        if !d.is_dir() {
            continue;
        }
        if d.file_name().unwrap() == "window_experiments" {
            for x in std::fs::read_dir(&d).unwrap().flatten() {
                if x.path().join("input.pcr").is_file() {
                    pcrs.push((x.path(), "input".to_string(), "ybfluo_Riet_IRF".to_string()));
                }
            }
            continue;
        }
        for f in std::fs::read_dir(&d).unwrap().flatten() {
            let n = f.file_name().to_string_lossy().into_owned();
            if let Some(code) = n.strip_suffix(".pcr") {
                if !code.ends_with("_calc") {
                    pcrs.push((d.clone(), code.to_string(), code.to_string()));
                }
            }
        }
    }
    pcrs.sort();
    let (mut n, mut refused, mut bytes) = (0, Vec::new(), 0);
    for (d, code, stem) in &pcrs {
        // A joint run's patterns each (P06, L322).
        // A file the reader refuses is counted as refused (none since P30's step read runs/texture, L360).
        let all = match drs::pcr::read_all(&latin1(&d.join(format!("{code}.pcr")))) {
            Ok(a) => a,
            Err(e) => {
                assert!(e.message.contains("not read yet"), "{code}: {e}");
                refused.push(format!("{code}: {e}"));
                continue;
            }
        };
        let several = all.len() > 1;
        for (k, mut p) in all.into_iter().enumerate() {
            p.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
            let l = match Model::from_pcr(&p, d, stem) {
                Ok(l) => l,
                Err(e) => {
                    refused.push(format!("{}/{code}{}: {e}", d.file_name().unwrap().to_string_lossy(), if several { format!("#{}", k + 1) } else { String::new() }));
                    continue;
                }
            };
            let name = format!("{}/{code}{}", d.file_name().unwrap().to_string_lossy(), if several { format!("#{}", k + 1) } else { String::new() });
            let t1 = l.model.to_text();
            let m2 = Model::from_text(&t1).unwrap_or_else(|e| panic!("{name}: {e}"));
            let t2 = m2.to_text();
            assert_eq!(t1, t2, "{name}: written again, the text differs");
            if let Some(data) = &l.data {
                let (a, b) = (l.model.calculate(&data.x).ycalc, m2.calculate(&data.x).ycalc);
                assert!(a == b, "{name}: the pattern of the model read back differs");
            }
            assert_eq!(l.model.fields(), m2.fields(), "{name}: the editable values differ");
            bytes += t1.len();
            n += 1;
        }
    }
    for r in &refused {
        println!("not imported: {r}");
    }
    println!("{n} models written and read back ({} kB of text), {} not imported", bytes / 1024, refused.len());
    assert!(n >= 47, "only {n} models");
    assert!(Model::from_text("{").is_err());
}
