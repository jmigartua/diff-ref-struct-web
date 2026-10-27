//! Writing a scan and the line of each point (app requests 12 and 15; L088), on every data file of the FullProf
//! material that the registry reads (the files `tests/data/file_kinds.txt` lists as data):
//! - `written_and_read_back`: each 2θ scan written as xye and as xys reads back with the same x, y and σ to the bit
//!   (a point with no weight stays with none);
//! - `lines_of_points`: where `point_lines` gives the lines, each point's line holds its y (and its x, for rows), the
//!   lines never go back, and the count of files per format with lines is reported; y2o3.dat (FullProf's free
//!   format) and a column file are checked line by line.

use drs::formats::{self, Axis, Probe};
use std::path::{Path, PathBuf};

fn data_files() -> Vec<PathBuf> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let list = std::fs::read_to_string(root.join("tests/data/file_kinds.txt")).unwrap();
    list.lines().filter_map(|l| l.split_once('\t')).filter(|(_, k)| *k == "data").map(|(p, _)| root.join("../manuals/fullprof").join(p)).filter(|p| p.exists()).collect()
}

#[test]
fn written_and_read_back() {
    let (mut files, mut scans, mut xye_rule) = (0, 0, 0);
    for f in data_files() {
        let bytes = std::fs::read(&f).unwrap();
        let Ok(d) = formats::read(&Probe::new(&f.to_string_lossy(), &bytes), None) else { continue };
        files += 1;
        for s in d.scans.iter().filter(|s| s.axis == Axis::TwoTheta && !s.pattern.is_empty()) {
            // The formats that hold 2θ (qye holds Q only, L195).
            for id in formats::write::WRITTEN.into_iter().filter(|id| *id != "qye") {
                let t = formats::write(s, id).unwrap();
                let back = formats::read(&Probe::new(&format!("x.{id}"), t.as_bytes()), Some(id)).unwrap_or_else(|e| panic!("{}: {id}: {e}", f.display()));
                let b = &back.scans[0].pattern;
                assert_eq!(&b.x, &s.pattern.x, "{}: {id}", f.display());
                // y and σ as written, except where the xye reader puts GSAS-II's rule on y ≤ 0: y = 0, no weight.
                for i in 0..b.len() {
                    let (y, sigma) = (s.pattern.y[i], s.pattern.sigma[i]);
                    let same = (b.y[i] == y && b.sigma[i] == sigma) || (id == "xye" && y <= 0.0 && b.y[i] == 0.0 && b.sigma[i].is_infinite());
                    assert!(same, "{}: {id}: (y, σ) at {} is ({y}, {sigma}) read back as ({}, {})", f.display(), s.pattern.x[i], b.y[i], b.sigma[i]);
                    xye_rule += usize::from(id == "xye" && !(b.y[i] == y && b.sigma[i] == sigma));
                }
            }
            scans += 1;
        }
    }
    println!("{files} data files read, {scans} 2θ scans written as xye and xys and read back to the bit; {xye_rule} points with y ≤ 0 read from xye as y = 0, no weight");
    assert!(scans > 100, "{scans}");
}

#[test]
fn lines_of_points() {
    let mut by_format: std::collections::BTreeMap<&str, (usize, usize)> = Default::default();
    for f in data_files() {
        let bytes = std::fs::read(&f).unwrap();
        let name = f.to_string_lossy();
        let p = Probe::new(&name, &bytes);
        let Ok(d) = formats::read(&p, None) else { continue };
        let text = p.text();
        let ls: Vec<&str> = text.lines().collect();
        for s in &d.scans {
            let e = by_format.entry(d.format).or_default();
            e.1 += 1;
            let Some(lines) = formats::point_lines(text, s) else { continue };
            e.0 += 1;
            assert_eq!(lines.len(), s.pattern.len());
            assert!(lines.windows(2).all(|w| w[0] <= w[1]), "{}: lines go back", f.display());
            for (i, &l) in lines.iter().enumerate() {
                let nums: Vec<f64> = ls[l - 1].split(|c: char| c.is_whitespace() || c == ',' || c == ';').filter_map(|t| t.parse().ok()).collect();
                assert!(nums.iter().any(|&v| (v - s.pattern.y[i]).abs() <= 1e-9 * v.abs().max(1.0)), "{}: point {i} not on line {l}", f.display());
            }
        }
    }
    println!("scans with the line of every point, by format: {by_format:?}");
    // FullProf's free format: y2o3.dat, a comment, Ti step Tf, then eight intensities a line.
    let y2o3 = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs/Y2O3/y2o3.dat");
    let bytes = std::fs::read(&y2o3).unwrap();
    let p = Probe::new("y2o3.dat", &bytes);
    let s = formats::read(&p, None).unwrap().scans.remove(0);
    let l = formats::point_lines(p.text(), &s).unwrap();
    assert_eq!((l[0], l[8], l[9], l[13090]), (3, 4, 4, 3 + 13090 / 8));
    // A column file written by the core: point i on line 4 + i (three header lines).
    let t = formats::write(&s, "xye").unwrap();
    let back = formats::read(&Probe::new("y.xye", t.as_bytes()), Some("xye")).unwrap().scans.remove(0);
    assert_eq!(formats::point_lines(&t, &back).unwrap()[..3], [4, 5, 6]);
    let free = by_format.get("free_steps").copied().unwrap_or_default();
    assert!(free.0 > 0 && free.0 == free.1, "every free-format scan has its lines: {free:?}");
}
