//! The unit cell against FullProf (test `fullprof_metric`).
//!
//! Every `.out` file under `manuals/fullprof/` has a "Metric information" block: the direct cell,
//! its volume, the reciprocal cell and volume, the metric tensors GD and GR, the Cartesian frame
//! (Cr_Orth_cel, Orth_Cr_cel) and the Busing–Levy matrix (BL_M, BL_Minv). From the printed cell
//! alone the core recomputes all of them, and each must agree with FullProf's number.
//!
//! Tolerance: half a unit of the last digit FullProf printed, plus 5·10⁻⁵ of the largest element of
//! the same quantity. The second term covers the rounding of the printed cell itself (lengths to
//! 10⁻⁴ Å, angles to 10⁻³°), from which the core has to start.

use drs::cell::Cell;
use drs::mat3::{self, Mat3};
use std::path::{Path, PathBuf};

const REL: f64 = 5e-5;

/// A number as printed, with the number of decimals it was printed with.
#[derive(Clone, Copy, Debug)]
struct Printed {
    v: f64,
    dec: i32,
}

fn printed(tok: &str) -> Option<Printed> {
    let v = tok.parse::<f64>().ok()?;
    let dec = tok.split_once('.').map_or(0, |(_, f)| f.len() as i32);
    Some(Printed { v, dec })
}

fn numbers(line: &str) -> Vec<Printed> {
    line.split(|c: char| c.is_whitespace() || c == '=')
        .filter_map(printed)
        .collect()
}

struct Block {
    cell: [Printed; 6],
    volume: Printed,
    recip: [Printed; 6],
    recip_volume: Printed,
    gd: [[Printed; 3]; 3],
    gr: [[Printed; 3]; 3],
    cr_orth: [[Printed; 3]; 3],
    orth_cr: [[Printed; 3]; 3],
    bl: [[Printed; 3]; 3],
    bl_inv: [[Printed; 3]; 3],
}

/// The six numbers after the first line that contains `key`, read from the next two lines with
/// numbers (three lengths, then three angles).
fn six_after(lines: &[&str], key: &str) -> [Printed; 6] {
    let i = lines
        .iter()
        .position(|l| l.contains(key))
        .unwrap_or_else(|| panic!("no '{key}'"));
    let v: Vec<Printed> = lines[i + 1..]
        .iter()
        .map(|l| numbers(l))
        .filter(|n| !n.is_empty())
        .take(2)
        .flatten()
        .collect();
    v.try_into()
        .unwrap_or_else(|v: Vec<_>| panic!("'{key}': expected 6 numbers, got {}", v.len()))
}

fn scalar_after(lines: &[&str], key: &str) -> Printed {
    let l = lines
        .iter()
        .find(|l| l.contains(key))
        .unwrap_or_else(|| panic!("no '{key}'"));
    *numbers(l).last().unwrap()
}

/// Two 3×3 matrices printed side by side in the three lines with six numbers after `key`.
fn pair_after(lines: &[&str], key: &str) -> ([[Printed; 3]; 3], [[Printed; 3]; 3]) {
    let i = lines
        .iter()
        .position(|l| l.contains(key))
        .unwrap_or_else(|| panic!("no '{key}'"));
    let rows: Vec<Vec<Printed>> = lines[i + 1..]
        .iter()
        .map(|l| numbers(l))
        .filter(|n| n.len() == 6)
        .take(3)
        .collect();
    assert_eq!(rows.len(), 3, "'{key}': expected three rows of six numbers");
    let m = |off: usize| [0, 1, 2].map(|r| [0, 1, 2].map(|c| rows[r][off + c]));
    (m(0), m(3))
}

fn read_block(text: &str) -> Block {
    let start = text
        .find("Metric information")
        .expect("no 'Metric information' block");
    let lines: Vec<&str> = text[start..].lines().take(60).collect();
    let (gd, gr) = pair_after(&lines, "GD");
    let (cr_orth, orth_cr) = pair_after(&lines, "Cr_Orth_cel");
    let (bl, bl_inv) = pair_after(&lines, "BL_M");
    Block {
        cell: six_after(&lines, "Direct cell parameters"),
        volume: scalar_after(&lines, "Direct Cell Volume"),
        recip: six_after(&lines, "Reciprocal cell parameters"),
        recip_volume: scalar_after(&lines, "Reciprocal Cell Volume"),
        gd,
        gr,
        cr_orth,
        orth_cr,
        bl,
        bl_inv,
    }
}

fn out_files(dir: &Path, found: &mut Vec<PathBuf>) {
    for e in std::fs::read_dir(dir).unwrap().flatten() {
        let p = e.path();
        if p.is_dir() {
            out_files(&p, found);
        } else if p.extension().is_some_and(|x| x.eq_ignore_ascii_case("out")) {
            found.push(p);
        }
    }
}

/// Compares one computed quantity with its printed counterpart; returns the worst excess, if any.
fn check(name: &str, computed: &[f64], printed: &[Printed], scale: f64, fails: &mut Vec<String>) {
    for (k, (c, p)) in computed.iter().zip(printed).enumerate() {
        let tol = 0.5 * 10f64.powi(-p.dec) + REL * scale;
        if (c - p.v).abs() > tol {
            fails.push(format!(
                "{name}[{k}]: core {c:.8}, FullProf {} (tolerance {tol:.2e})",
                p.v
            ));
        }
    }
}

fn flat(m: &Mat3) -> Vec<f64> {
    m.iter().flatten().copied().collect()
}

fn flat_p(m: &[[Printed; 3]; 3]) -> Vec<Printed> {
    m.iter().flatten().copied().collect()
}

fn scale(m: &[[Printed; 3]; 3]) -> f64 {
    m.iter().flatten().fold(0.0_f64, |s, p| s.max(p.v.abs()))
}

#[test]
fn fullprof_metric() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof");
    let mut files = Vec::new();
    out_files(&root, &mut files);
    files.sort();
    assert!(
        files.len() >= 13,
        "expected the 13 FullProf .out files in manuals/fullprof, found {}",
        files.len()
    );

    let mut report = Vec::new();
    for f in &files {
        let text = String::from_utf8_lossy(&std::fs::read(f).unwrap()).into_owned();
        // A run FullProf stopped before it read the phases has no metric block; left out only when its screen output
        // (CODE.run.log of the runs scripts) says why (runs/background/y2o3_bac62: "end of file *.BAC", L079).
        if !text.contains("Metric information") {
            let log = std::fs::read_to_string(f.with_extension("run.log")).unwrap_or_default();
            let why = log.lines().rev().find(|l| l.contains("end of file") || l.contains("ERROR") || l.contains("Error")).map(str::trim);
            let why = why.unwrap_or_else(|| panic!("{}: no 'Metric information' block and no stop in its run.log", f.display()));
            println!("left out {}: FullProf stopped before the phases ({why})", f.file_name().unwrap().to_string_lossy());
            continue;
        }
        let b = read_block(&text);
        let p = b.cell.map(|x| x.v);
        let cell = Cell::new(p[0], p[1], p[2], p[3], p[4], p[5]).unwrap();
        let r = cell.reciprocal();
        let mut fails = Vec::new();

        check(
            "volume",
            &[cell.volume()],
            &[b.volume],
            b.volume.v,
            &mut fails,
        );
        check(
            "reciprocal cell",
            &[r.a, r.b, r.c],
            &b.recip[..3],
            b.recip[0].v.max(b.recip[1].v).max(b.recip[2].v),
            &mut fails,
        );
        check(
            "reciprocal angles",
            &[r.alpha, r.beta, r.gamma],
            &b.recip[3..],
            180.0,
            &mut fails,
        );
        check(
            "reciprocal volume",
            &[r.volume()],
            &[b.recip_volume],
            b.recip_volume.v,
            &mut fails,
        );
        check(
            "GD",
            &flat(&cell.metric()),
            &flat_p(&b.gd),
            scale(&b.gd),
            &mut fails,
        );
        check(
            "GR",
            &flat(&cell.reciprocal_metric()),
            &flat_p(&b.gr),
            scale(&b.gr),
            &mut fails,
        );
        check(
            "Cr_Orth_cel",
            &flat(&cell.crystal_to_orthonormal()),
            &flat_p(&b.cr_orth),
            scale(&b.cr_orth),
            &mut fails,
        );
        check(
            "Orth_Cr_cel",
            &flat(&cell.orthonormal_to_crystal()),
            &flat_p(&b.orth_cr),
            scale(&b.orth_cr),
            &mut fails,
        );
        check(
            "BL_M",
            &flat(&cell.busing_levy_b()),
            &flat_p(&b.bl),
            scale(&b.bl),
            &mut fails,
        );
        let bl_inv = mat3::inverse(&cell.busing_levy_b()).unwrap();
        check(
            "BL_Minv",
            &flat(&bl_inv),
            &flat_p(&b.bl_inv),
            scale(&b.bl_inv),
            &mut fails,
        );

        let name = f.strip_prefix(&root).unwrap().display().to_string();
        println!(
            "{name}: a b c = {} {} {}, angles {} {} {}: {} disagreements",
            p[0],
            p[1],
            p[2],
            p[3],
            p[4],
            p[5],
            fails.len()
        );
        report.extend(fails.into_iter().map(|e| format!("{name}: {e}")));
    }
    assert!(
        report.is_empty(),
        "{} disagreements with FullProf:\n{}",
        report.len(),
        report.join("\n")
    );
}
