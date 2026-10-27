//! Sites against FullProf (tests `fullprof_sites` and `crysfml_symbols`).
//!
//! `tests/data/fullprof_sites.txt` was written by `tests/data/make_fullprof_sites.py` from the files in
//! `manuals/fullprof/`, every value copied as printed, with its file and line:
//! - MULT records: the atoms of the 13 FullProf `.out` files with FullProf's site multiplicity (96 atoms);
//! - BETA records: the atoms of the `.pcr` files refined with anisotropic betas (N_t = 2), with the six betas and
//!   their refinement codes (69 atoms, in P n m a, F m m m, P 21/a, I 21 3, B m a b and F m 3 m).
//!
//! For each MULT atom the core's multiplicity must be FullProf's. For each BETA atom:
//! - the refined betas must satisfy the core's relations (to the five decimals printed);
//! - a beta the core fixes at zero must be zero in the file and carry the code 0;
//! - a beta the core ties to another (β_i = c·β_k), when both carry codes, must share its parameter number with
//!   the code multipliers in the ratio c: FullProf's codeword is ±(10·p + a), p the parameter, a the multiplier.
//!
//! `crysfml_symbols` reads FullProf's own symbol table (CrysFML, `sg_fullprof`, 612 entries). Each Hall symbol must
//! give a valid group of the right point group. When the group is one of the 530 settings of the main table, it must
//! have the same number, the same setting as the main table's entry of that name if there is one, and the lookup
//! must find it from CrysFML's name. CrysFML also lists cells outside the 530 (A 1, I -1, F 2/m, C 4 2 21, F 4/m m m,
//! the origin choices ":1" of groups that have none in the main table); these are counted and listed, not failed.

use drs::sg_fullprof::FULLPROF_SYMBOLS;
use drs::site::Site;
use drs::spacegroup::{find, setting_of, SpaceGroup};

/// The fields of a line, with the quoted symbol as one field.
fn fields(line: &str) -> (String, Vec<&str>) {
    let (_, rest) = line.split_once('"').unwrap();
    let (sym, rest) = rest.split_once('"').unwrap();
    (sym.to_string(), rest.split_whitespace().collect())
}

fn f(s: &str) -> f64 {
    s.parse().unwrap()
}

#[test]
fn fullprof_sites() {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/fullprof_sites.txt");
    let text = std::fs::read_to_string(path).unwrap();
    let mut fails = Vec::new();
    let (mut n_mult, mut n_beta, mut n_codes) = (0, 0, 0);
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        let kind = line.split_whitespace().next().unwrap();
        let (sym, v) = fields(line);
        let at = v.last().unwrap();
        let x = [f(v[1]), f(v[2]), f(v[3])];
        let g = match kind {
            "MULT" => SpaceGroup::from_hall(&sym),
            _ => SpaceGroup::from_symbol(&sym).map(|(g, _)| g),
        };
        let g = match g {
            Ok(g) => g,
            Err(e) => {
                fails.push(format!("{at}: {sym}: {e}"));
                continue;
            }
        };
        let s = match Site::new(&g, x) {
            Ok(s) => s,
            Err(e) => {
                fails.push(format!("{at}: {e}"));
                continue;
            }
        };
        if kind == "MULT" {
            n_mult += 1;
            let m: usize = v[4].parse().unwrap();
            if m != s.multiplicity {
                fails.push(format!("{at}: {} multiplicity {m}, core {}", v[0], s.multiplicity));
            }
            continue;
        }
        n_beta += 1;
        let b: Vec<f64> = v[4..10].iter().map(|t| f(t)).collect();
        let code: Vec<f64> = v[10..16].iter().map(|t| f(t)).collect();
        let c = &s.adp;
        let free: Vec<f64> = c.free.iter().map(|&i| b[i]).collect();
        let want = c.apply(&free);
        for i in 0..6 {
            if (want[i] - b[i]).abs() > 1.5e-5 {
                fails.push(format!("{at}: {} {} = {} but {} from the core's {} = {}", v[0], c.names[i], b[i], want[i], c.names[i], c.expression(i)));
            }
            if c.is_fixed(i) && code[i] != 0.0 {
                fails.push(format!("{at}: {} {} fixed by symmetry but refined with code {}", v[0], c.names[i], code[i]));
            }
            // A tie to one free beta: compare the codes.
            let terms: Vec<(usize, f64)> =
                c.coef[i].iter().enumerate().filter(|(_, r)| !r.is_zero()).map(|(k, r)| (c.free[k], r.value())).collect();
            if let [(k, ratio)] = terms[..] {
                if k != i && code[i] != 0.0 && code[k] != 0.0 {
                    n_codes += 1;
                    let (pi, ai) = ((code[i].abs() / 10.0).floor(), code[i].signum() * (code[i].abs() % 10.0));
                    let (pk, ak) = ((code[k].abs() / 10.0).floor(), code[k].signum() * (code[k].abs() % 10.0));
                    if pi != pk || (ai - ratio * ak).abs() > 1e-6 {
                        fails.push(format!("{at}: {} codes {} and {} do not encode {} = {}", v[0], code[i], code[k], c.names[i], c.expression(i)));
                    }
                }
            }
        }
    }
    assert!(fails.is_empty(), "{} disagreements:\n{}", fails.len(), fails.join("\n"));
    assert_eq!((n_mult, n_beta), (96, 69), "records read");
    println!("{n_mult} multiplicities, {n_beta} anisotropic atoms, {n_codes} tied betas with codes on both sides");
}

#[test]
fn crysfml_symbols() {
    let (mut fails, mut outside) = (Vec::new(), Vec::new());
    let key = |s: &str| -> String { s.chars().filter(|c| !c.is_whitespace() && *c != '_').collect::<String>().to_uppercase() };
    for (n, hm, hall) in FULLPROF_SYMBOLS.iter() {
        let g = match SpaceGroup::from_hall(hall) {
            Ok(g) => g,
            Err(e) => {
                fails.push(format!("{hm} '{hall}': {e}"));
                continue;
            }
        };
        match setting_of(&g.ops) {
            Some(s) if s.number == *n => {
                let named: Vec<_> = drs::sgtable::SETTINGS
                    .iter()
                    .filter(|t| [t.hm, t.hm_full].iter().chain(t.international.split(" = ").collect::<Vec<_>>().iter()).any(|a| key(a) == key(hm)))
                    .collect();
                if !named.is_empty() && !named.iter().any(|t| t.hall_number == s.hall_number) {
                    fails.push(format!("{hm} '{hall}' is table entry {} ({}), but the table's '{hm}' is entry {}", s.hall_number, s.hm, named[0].hall_number));
                }
                if find(hm).first().map(|t| t.hall_number).is_none() {
                    fails.push(format!("{hm}: not found by the lookup"));
                }
            }
            Some(s) => fails.push(format!("{hm} '{hall}': number {n}, table entry {} is no. {}", s.hall_number, s.number)),
            None => outside.push(format!("{hm} ({hall}, no. {n})")),
        }
    }
    assert!(fails.is_empty(), "{} of {} CrysFML symbols disagree:\n{}", fails.len(), FULLPROF_SYMBOLS.len(), fails.join("\n"));
    println!("{} CrysFML symbols, {} outside the 530 settings: {}", FULLPROF_SYMBOLS.len(), outside.len(), outside.join("; "));
}
