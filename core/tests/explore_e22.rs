//! Explore's E22: the cell's box for any phase, `Cell::edges` (test `explore_e22`, L401).
//!
//! - a cubic cell: twelve edges, each a long, those from a corner at right angles;
//! - a monoclinic cell: a, b, c long, β between a and c, b at right angles to both;
//! - Y₂O₃: the same edges as `geometry::scene`'s, to the bit;
//! - a Le Bail phase (MgAl₂O₄'s profile matching, no atoms): its edges, from its cell alone.

use drs::cell::Cell;
use std::path::Path;

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn len(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn angle(u: [f64; 3], v: [f64; 3]) -> f64 {
    ((u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (len(u) * len(v))).acos().to_degrees()
}

/// The three edges from the origin: a, b, c.
fn axes(e: &[([f64; 3], [f64; 3])]) -> [[f64; 3]; 3] {
    let from: Vec<[f64; 3]> = e.iter().filter(|(p, _)| len(*p) < 1e-12).map(|(p, q)| sub(*q, *p)).collect();
    assert_eq!(from.len(), 3, "three edges from the origin");
    [from[0], from[1], from[2]]
}

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap().iter().map(|&b| b as char).collect()
}

#[test]
fn explore_e22() {
    let mut fails: Vec<String> = Vec::new();
    // Cubic.
    let e = Cell::new(10.6049, 10.6049, 10.6049, 90.0, 90.0, 90.0).unwrap().edges();
    if e.len() != 12 || e.iter().any(|(p, q)| (len(sub(*q, *p)) - 10.6049).abs() > 1e-12) {
        fails.push(format!("cubic: {} edges, lengths {:?}", e.len(), e.iter().map(|(p, q)| len(sub(*q, *p))).collect::<Vec<_>>()));
    }
    let [a, b, c] = axes(&e);
    for (x, y, n) in [(a, b, "γ"), (a, c, "β"), (b, c, "α")] {
        if (angle(x, y) - 90.0).abs() > 1e-9 {
            fails.push(format!("cubic: {n} = {}", angle(x, y)));
        }
    }
    // Monoclinic.
    let e = Cell::new(5.1, 6.2, 7.3, 90.0, 101.5, 90.0).unwrap().edges();
    let [a, b, c] = axes(&e);
    let got = [len(a), len(b), len(c), angle(b, c), angle(a, c), angle(a, b)];
    let want = [5.1, 6.2, 7.3, 90.0, 101.5, 90.0];
    if got.iter().zip(&want).any(|(g, w)| (g - w).abs() > 1e-9) {
        fails.push(format!("monoclinic: a b c α β γ {got:?}"));
    }
    // Every edge joins two of the eight corners and the box closes: the twelve edges hold each corner three times.
    let corners: Vec<[f64; 3]> = e.iter().flat_map(|(p, q)| [*p, *q]).collect();
    for k in &corners {
        let n = corners.iter().filter(|x| len(sub(**x, *k)) < 1e-12).count();
        if n != 3 {
            fails.push(format!("monoclinic: corner {k:?} in {n} edges"));
        }
    }
    // Y₂O₃: the scene's edges.
    let d = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs/Y2O3");
    let p = drs::pcr::read(&latin1(&d.join("y2o3_structure.pcr"))).unwrap();
    let m = drs::powder::Model::from_pcr(&p, &d, "y2o3_structure").unwrap().model;
    let s = drs::geometry::scene(&m, 0, [0.0; 3], [1.0; 3], &drs::geometry::SceneRule::Within(2.5)).unwrap();
    if s.edges != m.phases[0].cell.edges() {
        fails.push("Y2O3: the scene's edges are not Cell::edges".into());
    }
    // A Le Bail phase.
    let d = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs/MgAl2O4-LeBail");
    let p = drs::pcr::read(&latin1(&d.join("MgAl2O4s_LeBail.pcr.in"))).unwrap();
    let m = drs::powder::Model::from_pcr(&p, &d, "MgAl2O4s_LeBail").unwrap().model;
    let ph = &m.phases[0];
    if !matches!(ph.intensities, drs::powder::Intensities::List(_)) {
        fails.push("MgAl2O4-LeBail: not a phase without atoms".into());
    }
    let e = ph.cell.edges();
    if e.len() != 12 || e.iter().any(|(p, q)| (len(sub(*q, *p)) - ph.cell.a).abs() > 1e-12) {
        fails.push(format!("Le Bail: {} edges, a = {}", e.len(), ph.cell.a));
    }
    println!("cubic, monoclinic, Y2O3's scene and MgAl2O4's Le Bail phase checked");
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}
