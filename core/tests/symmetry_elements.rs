//! Symmetry elements (P09, the app's request 35) against cctbx's reading of every operator of the 530 settings,
//! centring translations included (`tests/data/symmetry_elements_530.txt`, written by `make_symmetry_elements.py`
//! from cctbx-base's `rot_mx_info` and `translation_part_info` over spglib's operators):
//! - the operators of each setting the core builds from its Hall symbol are cctbx's;
//! - for each operator: the rotation's order and sign; the axis (cctbx's eigenvector, up to its sign) and the sense
//!   of 3-, 4- and 6-fold axes taken with it; the intrinsic translation exactly; cctbx's point lies on the core's
//!   element (on its axis or plane, or at its point);
//! - the names: in each setting, every screw and glide letter of the full Hermann–Mauguin symbol is among the
//!   elements' symbols (2₁, 4₂ …, a, b, c, n, d; e and the double glides counted apart);
//! - the location strings of P 2₁/c, as International Tables Vol. A prints them.

use drs::spacegroup::SpaceGroup;
use drs::symmetry_elements::ElementKind;
use std::collections::HashMap;

fn data(f: &str) -> String {
    std::fs::read_to_string(format!("{}/tests/data/{f}", env!("CARGO_MANIFEST_DIR"))).unwrap()
}

/// (R, t) → (type, sense, ev, w in 1/24, p in 1/24).
type Ref = HashMap<([[i32; 3]; 3], [i32; 3]), (i32, i32, [i32; 3], [i32; 3], [i32; 3])>;

#[test]
fn symmetry_elements() {
    let refs = data("symmetry_elements_530.txt");
    let spg = data("spglib_530.txt");
    // The settings: Hall number, Hall symbol, full HM symbol (from spglib_530.txt's headers).
    let settings: Vec<(usize, String, String)> = spg
        .lines()
        .filter(|l| l.starts_with("# ") && l.contains('|') && l.as_bytes()[2].is_ascii_digit())
        .map(|l| {
            let f: Vec<&str> = l[2..].split('|').map(str::trim).collect();
            (f[0].split_whitespace().next().unwrap().parse().unwrap(), f[1].to_string(), f[3].to_string())
        })
        .collect();
    assert_eq!(settings.len(), 530);
    let mut table: HashMap<usize, Ref> = HashMap::new();
    let mut current = 0;
    for l in refs.lines() {
        if let Some(h) = l.strip_prefix("# ").and_then(|r| r.split('|').next()).and_then(|h| h.trim().parse().ok()) {
            current = h;
            continue;
        }
        if l.starts_with('#') {
            continue;
        }
        let parts: Vec<Vec<i32>> = l.split('|').map(|p| p.split_whitespace().map(|x| x.parse().unwrap()).collect()).collect();
        let r = [[parts[0][0], parts[0][1], parts[0][2]], [parts[0][3], parts[0][4], parts[0][5]], [parts[0][6], parts[0][7], parts[0][8]]];
        let t = [parts[0][9], parts[0][10], parts[0][11]];
        let v3 = |p: &Vec<i32>| [p[0], p[1], p[2]];
        table.entry(current).or_default().insert((r, t), (parts[1][0], parts[1][1], v3(&parts[2]), v3(&parts[3]), v3(&parts[4])));
    }
    assert_eq!(table.len(), 530);
    let (mut ops, mut fails) = (0usize, Vec::new());
    let (mut letters, mut letters_found, mut letters_apart) = (0, 0, Vec::new());
    for (h, hall, full) in &settings {
        let g = SpaceGroup::from_hall(hall).unwrap_or_else(|e| panic!("{h} {hall}: {e}"));
        let reference = &table[h];
        assert_eq!(g.ops.len(), reference.len(), "{h} {hall}: the operators");
        let elements = g.symmetry_elements();
        for (op, e) in &elements {
            ops += 1;
            let Some(&(ty, sense, ev, w24, p24)) = reference.get(&(op.r, op.t)) else {
                fails.push(format!("{h} {hall}: {} is not cctbx's", op.to_xyz()));
                continue;
            };
            let mut bad = Vec::new();
            if e.order != ty {
                bad.push(format!("order {} against {ty}", e.order));
            }
            // The axis, up to its sign, and the sense taken with it.
            let same = e.direction == ev;
            let opposite = e.direction == ev.map(|x| -x);
            if !(same || opposite) {
                bad.push(format!("axis {:?} against {ev:?}", e.direction));
            } else if sense != 0 && e.sense != if same { sense } else { -sense } {
                bad.push(format!("sense {} against {sense} (axis {:?}, cctbx's {ev:?})", e.sense, e.direction));
            }
            let w: Vec<i32> = e.intrinsic.iter().map(|x| (x * 24.0).round() as i32).collect();
            if w != w24 {
                bad.push(format!("intrinsic part {w:?}/24 against {w24:?}/24"));
            }
            // cctbx's point on the element: p − point in the span of its directions.
            if !matches!(e.kind, ElementKind::Identity | ElementKind::Translation) {
                let d: Vec<f64> = (0..3).map(|i| p24[i] as f64 / 24.0 - e.point[i]).collect();
                let off = match e.span.len() {
                    0 => d.iter().map(|x| x.abs()).fold(0.0, f64::max),
                    1 => {
                        let s = e.span[0];
                        let c = [d[1] * s[2] - d[2] * s[1], d[2] * s[0] - d[0] * s[2], d[0] * s[1] - d[1] * s[0]];
                        c.iter().map(|x| x.abs()).fold(0.0, f64::max)
                    }
                    2 => {
                        let (a, b) = (e.span[0], e.span[1]);
                        let n = [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
                        (0..3).map(|i| n[i] * d[i]).sum::<f64>().abs()
                    }
                    _ => 0.0,
                };
                if off > 1e-9 {
                    bad.push(format!("cctbx's point {p24:?}/24 is not on the element at {} ({:?})", e.location, e.point));
                }
            }
            if !bad.is_empty() {
                fails.push(format!("{h} {hall} {}: {} {}: {}", op.to_xyz(), e.symbol, e.location, bad.join("; ")));
            }
        }
        // The screw and glide letters of the full symbol (e.g. "P 1 21/c 1") among the elements' symbols.
        let symbols: Vec<&str> = elements.iter().map(|(_, e)| e.symbol.as_str()).collect();
        for part in full.split_whitespace().skip(1).flat_map(|p| p.split('/')) {
            let name = if part.len() == 2 && part.as_bytes()[0].is_ascii_digit() && part.as_bytes()[1].is_ascii_digit() && part.as_bytes()[1] < part.as_bytes()[0] {
                format!("{}_{}", &part[..1], &part[1..])
            } else if matches!(part, "a" | "b" | "c" | "n" | "d") {
                part.to_string()
            } else {
                continue;
            };
            letters += 1;
            if symbols.contains(&name.as_str()) {
                letters_found += 1;
            } else {
                letters_apart.push(format!("{h} {full}: {name}"));
            }
        }
    }
    println!("{ops} operators of 530 settings compared with cctbx; {} differ", fails.len());
    for f in fails.iter().take(20) {
        println!("  {f}");
    }
    println!("{letters} screw and glide letters of the full symbols, {letters_found} found among the elements; apart: {}", letters_apart.len());
    for f in letters_apart.iter().take(20) {
        println!("  {f}");
    }
    assert!(fails.is_empty(), "{} operators differ from cctbx", fails.len());
    // Apart, as measured in L260: R 3 c and R −3 c in rhombohedral axes, whose glide (1/2,1/2,1/2) the core names n
    // (ITA's diagonal glide, ½(a + b + c)) where the symbol keeps the c of the hexagonal axes.
    assert_eq!((letters, letters_found), (529, 527), "the letters of the symbols: {letters_apart:?}");
}

/// P 2₁/c (No. 14, unique axis b) as International Tables Vol. A lists its symmetry operations: (1) 1; (2) 2(0,1/2,0)
/// 0,y,1/4; (3) -1 0,0,0; (4) c x,1/4,z.
#[test]
fn p21c_as_ita() {
    let g = SpaceGroup::from_hall("-P 2ybc").unwrap();
    let e: Vec<(String, String, String)> = g.symmetry_elements().into_iter().map(|(o, e)| (o.to_xyz(), format!("{} {}", e.symbol, e.translation), e.location)).collect();
    println!("{e:#?}");
    let find = |xyz: &str| e.iter().find(|x| x.0 == xyz).unwrap_or_else(|| panic!("{xyz} in {e:?}")).clone();
    assert_eq!(find("-x,y+1/2,-z+1/2").1, "2_1 (0,1/2,0)");
    assert_eq!(find("-x,y+1/2,-z+1/2").2, "0,y,1/4");
    assert_eq!(find("-x,-y,-z").1.trim(), "-1");
    assert_eq!(find("-x,-y,-z").2, "0,0,0");
    assert_eq!(find("x,-y+1/2,z+1/2").1, "c (0,0,1/2)");
    assert_eq!(find("x,-y+1/2,z+1/2").2, "x,1/4,z");
}

/// Elements with their own locations, each against International Tables' line for the operator: the 3-fold along
/// [111] of P 2₁3, the −4 of P −4 2₁ c with its inversion point, and P n m a's mirror and n glide.
#[test]
fn ita_lines() {
    let check = |hm: &str, xyz: &str, symbol: &str, translation: &str, location: &str| {
        let g = SpaceGroup::from_symbol(hm).unwrap().0;
        let (_, e) = g.symmetry_elements().into_iter().find(|(o, _)| o.to_xyz() == xyz).unwrap_or_else(|| panic!("{hm}: {xyz}"));
        assert_eq!((e.symbol.as_str(), e.translation.as_str(), e.location.as_str()), (symbol, translation, location), "{hm} {xyz}");
    };
    // P 2₁3 (No. 198), (5) 3+ x,x,x.
    check("P 21 3", "z,x,y", "3", "", "x,x,x");
    // P −4 2₁ c (No. 114), (3) −4+ 0,0,z; 0,0,0.
    check("P -4 21 c", "y,-x,-z", "-4", "", "0,0,z; 0,0,0");
    // P n m a (No. 62): (6) m x,1/4,z and (8) n(0,1/2,1/2) 1/4,y,z.
    check("P n m a", "x,-y+1/2,z", "m", "", "x,1/4,z");
    check("P n m a", "-x+1/2,y+1/2,z+1/2", "n", "(0,1/2,1/2)", "1/4,y,z");
}
