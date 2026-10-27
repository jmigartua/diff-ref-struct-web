//! Sites against spglib (test `spglib_wyckoff`).
//!
//! `tests/data/wyckoff_230.txt` was written by `tests/data/make_wyckoff_reference.py` with spglib 2.7.0: for every
//! Wyckoff position of the standard setting of each of the 230 types, at two sets of free parameters, a structure
//! with an atom on the position and a second species on a general position, and spglib's `get_symmetry_dataset`
//! for it: the Wyckoff letter, the site-symmetry symbol, and the multiplicity (atoms in the conventional cell).
//! 3462 positions. For each, the core must find:
//! - the same multiplicity, Wyckoff letter and site-symmetry symbol, also from four other images of the position;
//! - multiplicity × order of the site group = order of the group (the property test);
//! - a site point group equal to the class of spglib's site-symmetry symbol;
//! - as many free coordinates as the Wyckoff position has free parameters, and constraints that reproduce the
//!   position;
//! - as many free β as the dimension of the symmetric tensors invariant under the site group, counted
//!   independently by the character formula (1/|S|) Σ ((tr R)² + tr(R²))/2, and β built from the constraints
//!   invariant under every site rotation.
//!
//! A second test checks, over all 530 settings, that no rotation mixes two axes of different length, so that the
//! U_ij of the CIF convention obey the same relations as the β_ij.

use drs::sgtable::SETTINGS;
use drs::site::Site;
use drs::spacegroup::SpaceGroup;
use drs::wyckoff_table::positions;

/// The crystal class of an oriented site-symmetry symbol: dots removed, orientation variants named as the class.
fn class_of(symbol: &str) -> String {
    let s: String = symbol.chars().filter(|&c| c != '.').collect();
    match s.as_str() {
        "m2m" | "2mm" => "mm2".into(),
        "-4m2" => "-42m".into(),
        "-62m" => "-6m2".into(),
        "321" | "312" => "32".into(),
        "3m1" | "31m" => "3m".into(),
        "-3m1" | "-31m" => "-3m".into(),
        "m3" => "m-3".into(),
        "m3m" => "m-3m".into(),
        "42m" => "-42m".into(),
        _ => s,
    }
}

fn trace(m: &[[i32; 3]; 3]) -> i32 {
    m[0][0] + m[1][1] + m[2][2]
}

fn square(m: &[[i32; 3]; 3]) -> [[i32; 3]; 3] {
    let mut r = [[0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = (0..3).map(|k| m[i][k] * m[k][j]).sum();
        }
    }
    r
}

#[test]
fn spglib_wyckoff() {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/wyckoff_230.txt");
    let text = std::fs::read_to_string(path).unwrap();
    let mut fails = Vec::new();
    let (mut n, mut types) = (0, std::collections::BTreeSet::new());
    let mut groups: std::collections::HashMap<u16, SpaceGroup> = std::collections::HashMap::new();
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = line.split_whitespace().collect();
        let hall: u16 = f[0].parse().unwrap();
        let number: u16 = f[1].parse().unwrap();
        let (letter, mult, sym) = (f[2].chars().next().unwrap(), f[3].parse::<usize>().unwrap(), f[4]);
        let x = [f[5].parse::<f64>().unwrap(), f[6].parse().unwrap(), f[7].parse().unwrap()];
        n += 1;
        types.insert(number);
        let g = groups.entry(hall).or_insert_with(|| SpaceGroup::from_setting(&SETTINGS[hall as usize - 1]));
        let tag = format!("#{hall} {} {letter} {x:?}", SETTINGS[hall as usize - 1].hm);
        let s = match Site::new(g, x) {
            Ok(s) => s,
            Err(e) => {
                fails.push(format!("{tag}: {e}"));
                continue;
            }
        };
        let w = s.wyckoff;
        if s.multiplicity != mult || w.map(|w| (w.letter, w.site)) != Some((letter, sym)) {
            fails.push(format!("{tag}: core {} {:?}, spglib {mult} {letter} {sym}", s.multiplicity, w.map(|w| (w.letter, w.site))));
            continue;
        }
        if s.multiplicity * s.ops.len() != g.ops.len() {
            fails.push(format!("{tag}: {} × {} ≠ {}", s.multiplicity, s.ops.len(), g.ops.len()));
        }
        if s.point_group != class_of(sym) {
            fails.push(format!("{tag}: site point group {} against spglib's {sym}", s.point_group));
        }
        // The same answer from other images of the position.
        for k in [1usize, 3, 7, 11] {
            let op = g.ops[(k * 5) % g.ops.len()];
            let y = [0, 1, 2].map(|i| {
                let v: f64 = (0..3).map(|j| op.r[i][j] as f64 * x[j]).sum::<f64>() + op.t[i] as f64 / 12.0;
                v - v.floor()
            });
            match Site::new(g, y) {
                Ok(t) if t.multiplicity == mult && t.wyckoff.map(|w| w.letter) == Some(letter) => {}
                Ok(t) => fails.push(format!("{tag}: image {y:?} gives {} {:?}", t.multiplicity, t.wyckoff.map(|w| w.letter))),
                Err(e) => fails.push(format!("{tag}: image {y:?}: {e}")),
            }
        }
        // Free coordinates: as many as the parameters of the Wyckoff position's coordinate triplet.
        let wy = positions(hall).iter().find(|e| e.letter == letter).unwrap();
        let params = ['x', 'y', 'z'].iter().filter(|c| wy.coord.contains(**c)).count();
        let c = &s.coordinates;
        if c.free.len() != params {
            fails.push(format!("{tag}: {} free coordinates, Wyckoff '{}' has {params}", c.free.len(), wy.coord));
        }
        let back = c.apply(&c.free.iter().map(|&i| s.position[i]).collect::<Vec<_>>());
        if (0..3).any(|i| (back[i] - s.position[i]).abs() > 1e-9) {
            fails.push(format!("{tag}: constraints give {back:?} for {:?}", s.position));
        }
        // Free betas: dimension of the invariant symmetric tensors by the character formula.
        let sum: i32 = s.ops.iter().map(|o| (trace(&o.r).pow(2) + trace(&square(&o.r))) / 2).sum();
        let dim = sum as usize / s.ops.len();
        if s.adp.free.len() != dim {
            fails.push(format!("{tag}: {} free betas, character formula {dim}", s.adp.free.len()));
        }
        let free: Vec<f64> = (0..s.adp.free.len()).map(|k| 0.0013 + 0.0007 * k as f64).collect();
        let b = s.adp.apply(&free);
        let m = [[b[0], b[3], b[4]], [b[3], b[1], b[5]], [b[4], b[5], b[2]]];
        for o in &s.ops {
            for i in 0..3 {
                for j in 0..3 {
                    let v: f64 = (0..3).flat_map(|k| (0..3).map(move |l| (k, l))).map(|(k, l)| o.r[i][k] as f64 * m[k][l] * o.r[j][l] as f64).sum();
                    if (v - m[i][j]).abs() > 1e-12 {
                        fails.push(format!("{tag}: β not invariant"));
                    }
                }
            }
        }
    }
    let shown: Vec<_> = fails.iter().take(40).cloned().collect();
    assert!(fails.is_empty(), "{} failures:\n{}", fails.len(), shown.join("\n"));
    assert_eq!(n, 3462, "positions read");
    assert_eq!(types.len(), 230, "space-group types covered");
}

/// The axes of equal length in the conventional cell of a setting: (a = b, b = c, a = c).
fn equal_axes(number: u16, choice: &str) -> [bool; 3] {
    match number {
        75..=194 if choice != "R" => [true, false, false],
        143..=167 => [true, true, true],
        195..=230 => [true, true, true],
        _ => [false, false, false],
    }
}

#[test]
fn u_ij_follow_beta_relations() {
    let mut fails = Vec::new();
    for s in SETTINGS.iter() {
        let eq = equal_axes(s.number, s.choice);
        let same = |i: usize, j: usize| i == j || match (i.min(j), i.max(j)) {
            (0, 1) => eq[0],
            (1, 2) => eq[1],
            _ => eq[2],
        };
        let g = SpaceGroup::from_setting(s);
        for op in &g.ops {
            for i in 0..3 {
                for j in 0..3 {
                    if op.r[i][j] != 0 && !same(i, j) {
                        fails.push(format!("#{} {}: {} mixes axes {i} and {j}", s.hall_number, s.hm, op.to_xyz()));
                    }
                }
            }
        }
    }
    assert!(fails.is_empty(), "{}", fails[..fails.len().min(20)].join("\n"));
}
