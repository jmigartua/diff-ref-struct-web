//! Space groups against spglib (test `spglib_530`).
//!
//! `tests/data/spglib_530.txt` holds, for each of the 530 settings, spglib's name of the setting and its stored
//! operators (spglib 2.7.0, written by `tests/data/make_spglib_reference.py`; spglib is not a dependency of the
//! core). For every setting the core generates the group from the Hall symbol of its own table (`sgtable`) and
//! checks, with no tolerance (integers and twelfths):
//! - that its table entry (number, Hall symbol, spglib's name string, full Hermann–Mauguin symbol) matches spglib's line;
//! - that the set of operators is spglib's set, operator for operator;
//! - its point group against spglib's, its crystal system against the number ranges of the International
//!   Tables (1–2 triclinic … 195–230 cubic), its centring vectors against the pure translations of spglib's set,
//!   and the presence of the inversion against spglib's set.

use drs::sgtable::SETTINGS;
use drs::spacegroup::SpaceGroup;
use drs::symop::SymOp;
use std::collections::BTreeSet;

struct Ref {
    hall_number: u16,
    number: u16,
    hall: String,
    international: String,
    hm_full: String,
    point_group: String,
    ops: BTreeSet<SymOp>,
}

fn read_reference() -> Vec<Ref> {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/spglib_530.txt");
    let text = std::fs::read_to_string(path).expect("tests/data/spglib_530.txt");
    let mut refs: Vec<Ref> = Vec::new();
    for line in text.lines() {
        if let Some(head) = line.strip_prefix("# ") {
            let mut f = head.split(" | ");
            let first: Vec<&str> = f.next().unwrap().split_whitespace().collect();
            let Ok(hall_number) = first[0].parse::<u16>() else {
                continue;
            };
            let hall = f.next().unwrap().to_string();
            let international = f.next().unwrap();
            refs.push(Ref {
                hall_number,
                number: first[1].parse().unwrap(),
                hall,
                international: international.to_string(),
                hm_full: f.next().unwrap().to_string(),
                point_group: f.next().unwrap().to_string(),
                ops: BTreeSet::new(),
            });
        } else if !line.trim().is_empty() {
            let v: Vec<i32> = line
                .split_whitespace()
                .map(|x| x.parse().unwrap())
                .collect();
            let r = [[v[0], v[1], v[2]], [v[3], v[4], v[5]], [v[6], v[7], v[8]]];
            refs.last_mut()
                .unwrap()
                .ops
                .insert(SymOp::new(r, [v[9], v[10], v[11]]));
        }
    }
    refs
}

fn system_of_number(n: u16) -> &'static str {
    match n {
        1..=2 => "triclinic",
        3..=15 => "monoclinic",
        16..=74 => "orthorhombic",
        75..=142 => "tetragonal",
        143..=167 => "trigonal",
        168..=194 => "hexagonal",
        _ => "cubic",
    }
}

#[test]
fn spglib_530() {
    let refs = read_reference();
    assert_eq!(refs.len(), 530);
    let mut failures = Vec::new();
    let mut n_ops = 0;
    for (s, r) in SETTINGS.iter().zip(&refs) {
        let tag = format!("#{} {} '{}'", s.hall_number, s.hm, s.hall);
        if (s.hall_number, s.number, s.hall, s.international, s.hm_full)
            != (
                r.hall_number,
                r.number,
                r.hall.as_str(),
                r.international.as_str(),
                r.hm_full.as_str(),
            )
        {
            failures.push(format!("{tag}: table entry differs from spglib's line"));
            continue;
        }
        let g = match SpaceGroup::from_hall(s.hall) {
            Ok(g) => g,
            Err(e) => {
                failures.push(format!("{tag}: {e}"));
                continue;
            }
        };
        let mine: BTreeSet<SymOp> = g.ops.iter().copied().collect();
        n_ops += mine.len();
        if mine.len() != g.ops.len() {
            failures.push(format!("{tag}: repeated operators"));
        }
        if mine != r.ops {
            let missing = r.ops.difference(&mine).count();
            let extra = mine.difference(&r.ops).count();
            failures.push(format!("{tag}: {} operators against spglib's {} ({missing} missing, {extra} not in spglib)", mine.len(), r.ops.len()));
            continue;
        }
        if g.point_group != r.point_group {
            failures.push(format!(
                "{tag}: point group {} against spglib's {}",
                g.point_group, r.point_group
            ));
        }
        if g.crystal_system != system_of_number(s.number) {
            failures.push(format!(
                "{tag}: crystal system {} for number {}",
                g.crystal_system, s.number
            ));
        }
        let pure: BTreeSet<[i32; 3]> = r
            .ops
            .iter()
            .filter(|o| o.r == drs::symop::IDENTITY)
            .map(|o| o.t)
            .collect();
        if pure != g.centring.iter().copied().collect() {
            failures.push(format!("{tag}: centring vectors differ"));
        }
        let inv = r
            .ops
            .iter()
            .any(|o| o.r == [[-1, 0, 0], [0, -1, 0], [0, 0, -1]]);
        if inv != g.centrosymmetric {
            failures.push(format!(
                "{tag}: centrosymmetric {} against spglib's {inv}",
                g.centrosymmetric
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "{} of 530 settings disagree:\n{}",
        failures.len(),
        failures.join("\n")
    );
    assert_eq!(n_ops, 7388, "operators compared");
}
