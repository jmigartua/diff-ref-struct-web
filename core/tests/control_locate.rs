//! Where a control file writes each value (`control::locate::locate`, the app's request 40, L299), on the control files
//! the core writes for eight FullProf runs (`Control::from_pcr`, `to_text`): for every value of each model that the
//! symmetry does not fix or tie, the bytes `locate` gives are replaced by the value plus 0.125, the file is read again,
//! and that value, and it alone of the values compared, has moved by 0.125. The line is the one of those bytes. A value
//! the file does not write (the micro-absorption at zero) and an address of no value give None.

use drs::control::locate::locate;
use drs::control::Control;
use drs::jacobian::Field;
use std::path::Path;

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap().iter().map(|&b| b as char).collect()
}

fn control_text(run: &str, code: &str) -> String {
    let d = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs").join(run);
    let mut p = drs::pcr::read(&latin1(&d.join(format!("{code}.pcr")))).unwrap();
    p.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
    Control::from_pcr(&p, &d, code, &d).unwrap_or_else(|e| panic!("{run}: {e}")).to_text()
}

#[test]
fn every_value_located() {
    let runs = [("Y2O3", "y2o3_structure"), ("NdNi-LTOa", "NdNi_LTOa"), ("NdNi-HTO", "NdNi_HTO"), ("prba10k", "prba10k"), ("Pd3MnD", "pdfinal"), ("modes/A_modes", "a_modes"), ("Co-YBCO", "Cocus-asym"), ("LaB6", "LaB6rx")];
    let mut total = 0;
    for (run, code) in runs {
        let text = control_text(run, code);
        let m = Control::from_text(&text).unwrap().model;
        let fields = m.fields();
        let mut located = 0;
        for fi in fields.iter().filter(|f| f.fixed_by.is_none()) {
            let Some((line, r)) = locate(&text, &fi.address) else {
                // Not written: only the micro-absorption at zero.
                assert!(matches!(fi.field, Field::Micro(_)) && fi.value == 0.0, "{run}: {} not located", fi.address);
                continue;
            };
            assert_eq!(line, text[..r.start].matches('\n').count() + 1, "{run}: {}", fi.address);
            let written: f64 = text[r.clone()].parse().unwrap_or_else(|_| panic!("{run}: {} at '{}'", fi.address, &text[r.clone()]));
            let moved = format!("{}{}{}", &text[..r.start], written + 0.125, &text[r.end..]);
            let m2 = Control::from_text(&moved).unwrap_or_else(|e| panic!("{run}: {}: {e:?}", fi.address)).model;
            let d = m2.get(fi.field) - m.get(fi.field);
            assert!((d - 0.125).abs() < 1e-9, "{run}: {} moved by {d} for {}", fi.address, &text[r.clone()]);
            // Nothing else of the same table moved (its neighbours in the list of fields).
            for other in fields.iter().filter(|o| o.field != fi.field && o.field.phase() == fi.field.phase() && o.fixed_by.is_none()).take(40) {
                assert!((m2.get(other.field) - m.get(other.field)).abs() < 1e-12, "{run}: {} moved {}", fi.address, other.address);
            }
            located += 1;
        }
        println!("{run:14} {located} values located in its control file");
        total += located;
    }
    assert!(total > 600, "{total}");
    // A value not written, an address of no value, and text that is not a control file.
    let text = control_text("Y2O3", "y2o3_structure");
    assert_eq!(locate(&text, "pattern1.microabsorption_p0"), None);
    assert_eq!(locate(&text, "phase2.cell_a"), None);
    assert_eq!(locate("not a control file", "pattern1.zero_shift"), None);
    // The second atom labelled Y, by its address with ':2'.
    let (line, r) = locate(&text, "phase1.Y:2.atom_position.x").unwrap();
    assert_eq!(&text[r], "-0.03252");
    println!("phase1.Y:2.atom_position.x on line {line}");
}

/// With the control file's folder (the app's request 41, L305): Y₂O₃'s control file with its background points in a
/// file of their own, as `drs new` writes them, and the points freed. Without the folder the file has no model and
/// nothing is located; with it every value is, and `all` (one reading) gives each value's place as `locate_in` does,
/// and `address_at` the address back from any byte of the value.
#[test]
fn with_the_folder_and_read_once() {
    let d = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs/Y2O3");
    let out = Path::new(env!("CARGO_MANIFEST_DIR")).join("target/test-control-locate");
    let _ = std::fs::remove_dir_all(&out);
    std::fs::create_dir_all(&out).unwrap();
    let p = drs::pcr::read(&latin1(&d.join("y2o3_structure.pcr"))).unwrap();
    let mut c = Control::from_pcr(&p, &d, "y2o3_structure", &out).unwrap();
    c.background_file = Some("y2o3 background.drb".into());
    let text = c.to_text();
    let (name, points) = c.background_file_text(None).unwrap().unwrap();
    std::fs::write(out.join(&name), points).unwrap();
    assert!(text.contains("file = \"y2o3 background.drb\"") && !text.contains("background_points = ["), "the points are not in the control file");
    assert!(text.contains("\"pattern1.background_points.0\""), "the points are refined");

    // Without the folder: no model, nothing located; with it, the cell.
    assert_eq!(locate(&text, "phase1.cell_a"), None);
    let (line, r) = drs::control::locate::locate_in(&text, "phase1.cell_a", Some(&out)).unwrap();
    assert!(text.lines().nth(line - 1).unwrap().starts_with("cell_a = "), "line {line}");
    assert_eq!(text[r.clone()].parse::<f64>().unwrap(), c.model.phases[0].cell.a);

    // Every value from one reading, as locate_in gives it, and the address back from its first, a middle and its
    // past-the-end byte.
    let t = std::time::Instant::now();
    let places = drs::control::locate::all(&text, Some(&out));
    let once = t.elapsed();
    let written = c.model.fields().into_iter().filter(|f| !matches!(f.field, Field::Background(_) | Field::Micro(_))).count();
    assert_eq!(places.len(), written, "every value but the points (in their own file) and the micro-absorption at zero");
    let t = std::time::Instant::now();
    for (a, line, r) in &places {
        assert_eq!(drs::control::locate::locate_in(&text, a, Some(&out)), Some((*line, r.clone())), "{a}");
    }
    let each = t.elapsed();
    for (a, _, r) in &places {
        for b in [r.start, (r.start + r.end) / 2, r.end] {
            assert_eq!(drs::control::locate::address_at(&text, b, Some(&out)).as_deref(), Some(a.as_str()), "{a} at byte {b}");
        }
    }
    let key = text.find("cell_a").unwrap();
    assert_eq!(drs::control::locate::address_at(&text, key, Some(&out)), None, "a key is not a value");
    println!("{} values located: all at once in {once:?}, one by one in {each:?}", places.len());
}
