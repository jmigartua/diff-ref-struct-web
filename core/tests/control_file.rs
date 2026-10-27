//! The control file (L066).
//! - `control_round_trip`: every FullProf run the core imports (the runs of `manuals/fullprof/runs`, the window
//!   experiments) is imported into the default model with its refinement, written as a control file, and read back:
//!   the model is the same to the bit (its JSON text, `Model::to_text`, identical), so are the calculated pattern and
//!   χ² at the data's points, the refinement's parameter set is the one `ParameterSet::from_pcr` builds, the data file
//!   is found again from its relative path, and the text written again is identical.
//! - `comments_and_order_kept`: a file edited by hand (comments, keys moved, values written their own way) and then
//!   updated after the model changed keeps everything but the changed value, and reads back to the changed model.
//! - `errors_have_places`: every kind of mistake is reported with its line and column.
//! - `pcr_export_round_trip`: every run's control exported to a `.pcr` (with its resolution and reflection files) and
//!   imported again by the core's `.pcr` reader: the same model to the bit, the same parameter set; what FullProf
//!   cannot express is listed. FullProf itself is not run here (cloud).
#![cfg(feature = "serde")]

use drs::control::{schema, syntax, update, Control};
use drs::powder::Model;
use drs::refinable::ParameterSet;
use std::path::{Path, PathBuf};

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display())).iter().map(|&b| b as char).collect()
}

/// The runs imported: (folder, .pcr stem, data stem).
fn runs() -> Vec<(PathBuf, String, String)> {
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
    pcrs
}

#[test]
fn control_round_trip() {
    let (mut n, mut with_set, mut refused, mut bytes) = (0, 0, Vec::new(), 0);
    for (d, code, stem) in runs() {
        // A joint run's patterns each (P06, L322).
        // A file the reader refuses is counted as refused (none since P30's step read runs/texture, L360).
        let all = match drs::pcr::read_all(&latin1(&d.join(format!("{code}.pcr")))) {
            Ok(a) => a,
            Err(e) => {
                assert!(e.message.contains("not read yet"), "{code}: {e}");
                refused.push(format!("{}/{code}: {e}", d.file_name().unwrap().to_string_lossy()));
                continue;
            }
        };
        let several = all.len() > 1;
        for (k, mut p) in all.into_iter().enumerate() {
            let name = format!("{}/{code}{}", d.file_name().unwrap().to_string_lossy(), if several { format!("#{}", k + 1) } else { String::new() });
            p.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
            let l = match Model::from_pcr(&p, &d, &stem) {
                Ok(l) => l,
                Err(e) => {
                    refused.push(format!("{name}: {e}"));
                    continue;
                }
            };
            let c = match Control::from_pcr(&p, &d, &stem, &d) {
                Ok(c) => c,
                Err(e) => {
                    // A refinement the core does not refine yet (the model alone is checked by model_text).
                    refused.push(format!("{name}: {e}"));
                    continue;
                }
            };
            let t1 = c.to_text();
            let c2 = Control::from_text_in(&t1, &d).unwrap_or_else(|e| panic!("{name}: {}\n{t1}", e.iter().map(|d| d.to_string()).collect::<Vec<_>>().join("\n")));
            assert_eq!(c.model.to_text(), c2.model.to_text(), "{name}: the model read back differs");
            assert_eq!(t1, c2.to_text(), "{name}: written again, the text differs");
            assert_eq!(c.title, c2.title);
            // The data file from its relative path, and the pattern and χ² at its points.
            if let (Some(data), Some(src)) = (&l.data, &c2.data) {
                let path = d.join(&src.file);
                assert_eq!(std::fs::canonicalize(&path).unwrap(), std::fs::canonicalize(l.data_path.as_ref().unwrap()).unwrap(), "{name}");
                let (a, b) = (l.model.calculate(&data.x).ycalc, c2.model.calculate(&data.x).ycalc);
                assert!(a == b, "{name}: the pattern of the model read back differs");
                let s1 = drs::lsq::statistics(&l.model, data, &a, 0);
                let s2 = drs::lsq::statistics(&c2.model, data, &b, 0);
                assert!(s1.chi2.to_bits() == s2.chi2.to_bits(), "{name}: χ² {} and {}", s1.chi2, s2.chi2);
            }
            // The refinement: the parameter set of the .pcr, exactly.
            if let Ok(set) = ParameterSet::from_pcr(&p, &l.model) {
                assert_eq!(c2.refinement.build(&c2.model).unwrap(), set, "{name}: the parameter set differs");
                with_set += 1;
            }
            bytes += t1.len();
            n += 1;
        }
    }
    for r in &refused {
        println!("not converted: {r}");
    }
    println!("{n} runs written as control files and read back ({} kB), {with_set} with their parameter sets; {} not converted", bytes / 1024, refused.len());
    // The runs not converted, each for its reason: a strain model FullProf describes in words (L211), hkl-dependent
    // shifts (P26), and two of P05's TOF variants (L334, L341): arg_si_abs2b's ABSCOR2 (FullProf stopped on it,
    // "Singular matrix") and arg_si_plain, whose .pcr FullProf rewrote with NaN after its one cycle; of L346's runs, the
    // Fourier filtering of the background (Nba = −2, not reproduced: L353). The texture with several directions (Nor = 2,
    // 3) is converted since P30 (L360), the rigid groups with TLS (rigid/'s eleven) since P41's TLS (L367).
    let names: Vec<&str> = refused.iter().map(|r| r.split(':').next().unwrap()).collect();
    let expected = ["NdNi-LTOb/NdNi_LTOb", "arg_si/arg_si_abs2b", "arg_si/arg_si_plain", "background_terms/y2o3_filter", "highpress/highpress"];
    assert_eq!(names, expected, "the runs not converted");
    assert!(n >= 110, "{n} runs");
}

#[test]
fn comments_and_order_kept() {
    let d = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs/Y2O3");
    let p = drs::pcr::read(&latin1(&d.join("y2o3_structure.pcr"))).unwrap();
    let c = Control::from_pcr(&p, &d, "y2o3_structure", &d).unwrap();
    let canonical = c.to_text();
    // The person's edits: a comment, a value written their own way, two keys swapped, a key the canonical writer
    // writes (the transparency) removed.
    let mut text = canonical.replace("[pattern]\n", "[pattern]\n# measured on the D8, 2019\n");
    let zero_line = text.lines().find(|l| l.starts_with("zero_shift")).unwrap().to_string();
    let disp_line = text.lines().find(|l| l.starts_with("sample_displacement")).unwrap().to_string();
    text = text.replace(&format!("{zero_line}\n{disp_line}\n"), &format!("{disp_line}\n{zero_line}   # checked\n"));
    let transparency = text.lines().find(|l| l.starts_with("sample_transparency")).unwrap().to_string();
    text = text.replace(&format!("{transparency}\n"), "");
    let w = c.model.lambda[0];
    text = text.replacen(&format!("wavelength_1 = {}", syntax::float(w)), "wavelength_1 = 1.54060000", 1);
    let edited = Control::from_text(&text).unwrap();
    assert_eq!(edited.model.to_text(), c.model.to_text(), "the edits change no value");
    // The model changes: the zero shift and one atom's B.
    let mut changed = edited.clone();
    changed.model.zero = 0.0123;
    if let drs::powder::Intensities::Structure(st, _) = &mut changed.model.phases[0].intensities {
        st.atoms[0].biso = 0.4567;
    }
    let out = update::update(&text, &changed);
    assert!(out.contains("# measured on the D8, 2019\n"));
    assert!(out.contains("zero_shift = 0.0123"), "{out}");
    assert!(out.contains("   # checked\n"));
    assert!(out.contains("wavelength_1 = 1.54060000"), "a value written its own way stays when it did not change");
    assert!(out.find("sample_displacement").unwrap() < out.find("zero_shift").unwrap(), "the order is kept");
    assert!(out.contains("sample_transparency = "), "a key the person removed comes back at the end of its table");
    assert!(out.contains("atom_biso = 0.4567"));
    let back = Control::from_text(&out).unwrap();
    assert_eq!(back.model.to_text(), changed.model.to_text());
    // The lines that differ: the two changed values and the key added back.
    let diff = out.lines().filter(|l| !text.lines().any(|m| m == *l)).count();
    assert_eq!(diff, 3, "{out}");
    // Removing an atom removes its table only.
    let mut fewer = changed.clone();
    if let drs::powder::Intensities::Structure(st, sc) = &mut fewer.model.phases[0].intensities {
        st.atoms.pop();
        sc.pop();
    }
    fewer.refinement.parameters.retain(|t| t.iter().all(|&(f, _)| !matches!(f, drs::jacobian::Field::AtomX(0, 1, _) | drs::jacobian::Field::AtomB(0, 1) | drs::jacobian::Field::AtomOcc(0, 1))));
    let out2 = update::update(&out, &fewer);
    assert_eq!(out2.matches("[[phase.atom]]").count(), out.matches("[[phase.atom]]").count() - 1);
    assert_eq!(Control::from_text(&out2).unwrap().model.to_text(), fewer.model.to_text());
}

#[test]
fn errors_have_places() {
    let d = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs/Y2O3");
    let p = drs::pcr::read(&latin1(&d.join("y2o3_structure.pcr"))).unwrap();
    let good = Control::from_pcr(&p, &d, "y2o3_structure", &d).unwrap().to_text();
    let line_of = |t: &str, key: &str| t.lines().position(|l| l.starts_with(key)).unwrap() + 1;
    let cases: Vec<(String, usize, &str)> = vec![
        (good.replacen("cell_a =", "cel_a =", 1), line_of(&good, "cell_a"), "did you mean \"cell_a\""),
        (good.replacen("geometry = \"debye-scherrer\"", "geometry = \"bragg\"", 1), line_of(&good, "geometry"), "is not one of"),
        (good.replacen("species = \"Y\"", "species = \"Qq\"", 1), line_of(&good, "species = \"Y\""), "\"Qq\""),
        (good.replacen("reader = \"free_steps\"", "reader = \"free_step\"", 1), line_of(&good, "reader"), "did you mean \"free_steps\""),
        (good.replacen("hall = ", "hall = 3 #", 1), line_of(&good, "hall"), "must be a quoted string"),
        (good.replacen("\"phase1.scale_factor\"", "\"phase1.scale\"", 1), line_of(&good, "  [\"phase1.scale"), "names no value"),
        (good.replacen("format = \"diff-ref-struct 1\"", "format = \"diff-ref-struct 9\"", 1), 2, "this version of the core reads"),
        (good.replacen("[[phase.atom]]", "[[phase.atoms]]", 1), line_of(&good, "[[phase.atom]]"), "did you mean \"phase.atom\""),
    ];
    for (text, line, msg) in &cases {
        let (c, diags) = schema::read_text(text);
        assert!(c.is_none(), "{msg}: read without an error");
        assert!(diags.iter().any(|d| d.span.line == *line && d.message.contains(msg)), "{msg} at line {line}: {:?}", diags);
    }
    // The outline of a file: its phase and atoms with their lines.
    let doc = syntax::parse(&good);
    let o = doc.outline();
    assert!(o.iter().any(|i| i.title.starts_with("phase ") && i.first_line == line_of(&good, "[[phase]]")));
    println!("{} mistakes placed; outline of {} items", cases.len(), o.len());
}


/// The parameter set as a comparable form: each parameter's fields and multipliers, sorted; parameters sorted.
fn set_form(s: &ParameterSet) -> Vec<Vec<(drs::jacobian::Field, f64)>> {
    let mut v: Vec<Vec<(drs::jacobian::Field, f64)>> = s
        .free
        .iter()
        .map(|f| {
            let mut x = f.fields.clone();
            x.sort_by_key(|a| a.0);
            x
        })
        .collect();
    v.sort_by(|a, b| a[0].0.cmp(&b[0].0));
    v
}

/// The lines of the exported `.pcr` that FullProf refused or read differently (L069), against the original's fields,
/// numbers compared as numbers: the resolution file's name alone on its line; Thmin, Step and Thmax as the `.pcr` gave
/// them; the excluded regions as written; VARY and FIX commands straight after the wavelength line. (The codewords'
/// numbers, 1…N with no hole, are checked in `pcr_export_round_trip`.)
fn lines_as_written(name: &str, p: &drs::pcr::Pcr, p2: &drs::pcr::Pcr, e: &drs::control::pcr_export::Export) {
    let text: Vec<&str> = e.pcr.lines().collect();
    let data: Vec<(usize, &str)> = text.iter().enumerate().filter(|(_, l)| !l.trim_start().starts_with('!') && !l.trim().is_empty()).map(|(i, l)| (i, *l)).collect();
    // The resolution file: the third line that is not a comment (after COMM and the flags), one field, a file written.
    if p.resolution_file.is_some() {
        let l = data[2].1;
        let fields: Vec<&str> = l.split_whitespace().collect();
        assert_eq!(fields.len(), 1, "{name}: the resolution file's line is '{l}'");
        assert!(e.files.iter().any(|f| f.0 == fields[0]), "{name}: {} is not written", fields[0]);
        assert_eq!(p2.resolution_file.as_deref(), Some(fields[0]), "{name}");
    }
    // Thmin, Step, Thmax.
    assert!(p2.thmin == p.thmin && p2.step == p.step && p2.thmax == p.thmax, "{name}: range {} {} {} written {} {} {}", p.thmin, p.step, p.thmax, p2.thmin, p2.step, p2.thmax);
    // The excluded regions.
    assert_eq!(p2.excluded, p.excluded, "{name}: excluded regions");
    // Commands: the original's VARY ratio kept, every command line straight after the wavelength line.
    let ratio = |q: &drs::pcr::Pcr| q.commands.iter().any(|c| c.split_whitespace().map(str::to_uppercase).collect::<Vec<_>>() == ["VARY", "RATIO"]);
    assert_eq!(ratio(p2), ratio(p), "{name}: VARY ratio");
    // Line 8: the wavelengths, or a TOF pattern's Bkpos, Wdt and Iabscor (P05).
    let lambda = text.iter().position(|l| l.trim_start().starts_with("! Lambda1") || l.trim_start().starts_with("!  Bkpos")).expect("the wavelength header") + 1;
    for (i, l) in &data {
        let w = l.split_whitespace().next().unwrap_or("").to_uppercase();
        if matches!(w.as_str(), "VARY" | "FIX") {
            let k = text[lambda + 1..*i].iter().filter(|t| !t.trim_start().starts_with('!')).count();
            let commands_before = text[lambda + 1..*i].iter().all(|t| matches!(t.split_whitespace().next().map(str::to_uppercase).as_deref(), Some("VARY" | "FIX")));
            assert!(commands_before, "{name}: '{l}' at line {} is not straight after the wavelength line ({k} lines between)", i + 1);
        }
    }
}

#[test]
fn pcr_export_round_trip() {
    let tmp = std::env::temp_dir().join(format!("drs_pcr_export_{}", std::process::id()));
    let (mut n, mut exact, mut worst_y, mut worst_chi, mut notes) = (0, 0, 0f64, 0f64, std::collections::BTreeMap::<String, usize>::new());
    let mut rewritten: Vec<String> = Vec::new();
    let (mut fields_checked, mut numbers_checked) = (0, 0);
    for (d, code, stem) in runs() {
        // A joint run's patterns each (P06, L322).
        // A file the reader refuses is left out (none since P30's step read runs/texture, L360).
        let Ok(all) = drs::pcr::read_all(&latin1(&d.join(format!("{code}.pcr")))) else { continue };
        let several = all.len() > 1;
        for (k, mut p) in all.into_iter().enumerate() {
            let name = format!("{}/{code}{}", d.file_name().unwrap().to_string_lossy(), if several { format!("#{}", k + 1) } else { String::new() });
            p.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
            let Ok(l) = Model::from_pcr(&p, &d, &stem) else { continue };
            let Ok(c) = Control::from_pcr(&p, &d, &stem, &d) else { continue };
            // Through the control file's text, as `drs convert` goes.
            let c = Control::from_text_in(&c.to_text(), &d).unwrap_or_else(|_| panic!("{name}: the control file does not read"));
            let e = drs::control::pcr_export::export(&c, "rt", l.data.as_ref()).unwrap_or_else(|er| panic!("{name}: {er}"));
            for note in &e.notes {
                *notes.entry(note.split(':').next().unwrap().to_string()).or_default() += 1;
            }
            let dir = tmp.join(name.replace('/', "_"));
            std::fs::create_dir_all(&dir).unwrap();
            std::fs::write(dir.join("rt.pcr"), &e.pcr).unwrap();
            for (f, t) in &e.files {
                std::fs::write(dir.join(f), t).unwrap();
            }
            let p2 = drs::pcr::read(&e.pcr).unwrap_or_else(|er| panic!("{name}: the exported .pcr does not read: {er}\n{}", e.pcr));
            lines_as_written(&name, &p, &p2, &e);
            fields_checked += 1;
            let data_path = e.files.iter().find(|f| f.0.ends_with(".xys")).map(|f| dir.join(&f.0)).or(l.data_path.clone());
            let l2 = Model::from_pcr_with_data(&p2, &dir, "rt", data_path.as_deref()).unwrap_or_else(|er| panic!("{name}: the exported .pcr does not import: {er}"));
            // The codewords' parameter numbers: 1…N with no hole, N as line 13 declares (L071).
            let nums: Vec<u32> = drs::refinable::pcr_parameter_numbers(&p2, &l2.model).unwrap().keys().copied().collect();
            assert_eq!(nums, (1..=nums.len() as u32).collect::<Vec<_>>(), "{name}: the parameter numbers have a hole");
            assert_eq!(nums.len(), p2.declared_parameters, "{name}: {} numbered, {} declared", nums.len(), p2.declared_parameters);
            numbers_checked += 1;
            // The background's profiles (MORE = 2, 3; P23) and a .bac's b1–b4 (the app's request 47) are named by where
            // each model read them: their points compared.
            let unnamed = |m: &Model| {
                let mut m = m.clone();
                m.background_terms.profiles.iter_mut().for_each(|p| p.file.clear());
                if let Some(f) = &mut m.background_terms.file {
                    f.file.clear();
                }
                m.to_text()
            };
            if unnamed(&l2.model) == unnamed(&c.model) {
                exact += 1;
            } else if l.data.as_ref().is_some_and(|d| matches!(&l2.model.background, drs::powder::BackgroundModel::Points(q) if q.len() == d.x.len() && q.len() > 9999)) {
                // A .bac of points at other angles (pairs, FullProf's second format) is
                // written at the data's points, the only form FullProf reads (L082): the same pattern (checked below),
                // other points.
                exact += 1;
                rewritten.push(name.clone());
            } else {
                println!("{name}: not the same model");
            }
            if let Some(data) = &l.data {
                let (a, b) = (c.model.calculate(&data.x).ycalc, l2.model.calculate(&data.x).ycalc);
                for (u, v) in a.iter().zip(&b) {
                    worst_y = worst_y.max((u - v).abs() / u.abs().max(1e-300));
                }
                let (s1, s2) = (drs::lsq::statistics(&c.model, data, &a, 0), drs::lsq::statistics(&l2.model, data, &b, 0));
                worst_chi = worst_chi.max((s1.chi2 - s2.chi2).abs() / s1.chi2);
            }
            let s1 = set_form(&c.refinement.build(&c.model).unwrap());
            let s2 = set_form(&ParameterSet::from_pcr(&p2, &l2.model).unwrap_or_else(|er| panic!("{name}: {er}")));
            assert_eq!(s1.len(), s2.len(), "{name}: free parameters");
            for (x, y) in s1.iter().zip(&s2) {
                assert_eq!(x.len(), y.len(), "{name}: {x:?} / {y:?}");
                for (u, v) in x.iter().zip(y) {
                    assert!(u.0 == v.0 && (u.1 - v.1).abs() < 1e-12 * u.1.abs().max(1.0), "{name}: {x:?} / {y:?}");
                }
            }
            n += 1;
        }
    }
    let _ = std::fs::remove_dir_all(&tmp);
    println!("background points rewritten at the data's points: {rewritten:?}");
    println!("{n} runs exported to .pcr and imported again: {exact} models identical to the bit; pattern within {worst_y:.1e}, χ² within {worst_chi:.1e}; parameter sets the same");
    for (k, v) in &notes {
        println!("  note ({v} runs): {k}");
    }
    assert!(n >= 110, "{n} runs exported");
    assert_eq!(fields_checked, n, "the lines FullProf refused (L069) checked on every run");
    assert_eq!(numbers_checked, n, "the parameter numbers checked on every run");
    assert_eq!(exact, n, "every model exported and imported again is the same to the bit");
    // y2o3_bac_b3b4's points are kept as read since b1–b4 are values (the app's request 47, L364): no longer rewritten.
    assert_eq!(rewritten, ["background/y2o3_bac62", "background/y2o3_bacpairs"], "only the pairs .bac are rewritten");
    assert!(worst_y == 0.0 && worst_chi == 0.0, "pattern {worst_y:e}, χ² {worst_chi:e}");
}
