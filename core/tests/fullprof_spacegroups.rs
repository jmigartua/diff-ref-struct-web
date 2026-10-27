//! Space groups against FullProf (test `fullprof_spacegroups`).
//!
//! Every `.out` file under `manuals/fullprof/` (13 files) has an "Information on Space Group" block: the number,
//! the Hermann–Mauguin and Hall symbols, the setting choice (when not the standard one), crystal system, Laue
//! class, point group, general multiplicity, centrosymmetry and the list of operators without the inversion and
//! the centring translations (SYMM lines; one file lists every operator instead). The groups met are P n m a, P 21/c, P 21/a and P 2/m.
//!
//! For each block the core looks up FullProf's Hermann–Mauguin symbol in its table and must find FullProf's Hall
//! symbol, number and setting choice; it generates the group from that Hall symbol and must agree on crystal
//! system, Laue class, point group, multiplicity and centrosymmetry; and the full set of operators, rebuilt from
//! FullProf's SYMM list with the inversion (when centric) and the lattice centring, must be the core's set.

use drs::spacegroup::SpaceGroup;
use drs::symop::SymOp;
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

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

/// The value after "<key>:" on the first line at or after `from` that contains it.
fn field<'a>(lines: &[&'a str], from: usize, key: &str) -> Option<&'a str> {
    lines[from..]
        .iter()
        .find(|l| l.contains(key))
        .and_then(|l| l.split_once(':'))
        .map(|(_, v)| v.trim())
}

/// Two Hall symbols that generate the same set of operators.
fn same_group(a: &str, b: &str) -> bool {
    match (SpaceGroup::from_hall(a), SpaceGroup::from_hall(b)) {
        (Ok(x), Ok(y)) => x.ops.iter().collect::<BTreeSet<_>>() == y.ops.iter().collect::<BTreeSet<_>>() && x.centring.iter().collect::<BTreeSet<_>>() == y.centring.iter().collect::<BTreeSet<_>>(),
        _ => false,
    }
}

#[test]
fn fullprof_spacegroups() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof");
    let mut files = Vec::new();
    out_files(&root, &mut files);
    files.sort();
    // 13 tutorial files, plus the FullProf runs of the reference refinements in manuals/fullprof/runs
    assert!(
        files.len() >= 13,
        "expected at least the 13 FullProf .out files of the tutorials in manuals/fullprof, found {}",
        files.len()
    );
    let mut fails = Vec::new();
    let (mut blocks, mut groups) = (0, BTreeSet::new());
    for f in &files {
        let raw = std::fs::read(f).unwrap();
        let text = String::from_utf8_lossy(&raw);
        let lines: Vec<&str> = text.lines().collect();
        let name = f.strip_prefix(&root).unwrap().display().to_string();
        for start in lines
            .iter()
            .enumerate()
            .filter(|(_, l)| l.contains("Information on Space Group"))
            .map(|(i, _)| i)
        {
            blocks += 1;
            let get =
                |k: &str| field(&lines, start, k).unwrap_or_else(|| panic!("{name}: no '{k}'"));
            let (hm, hall) = (get("Hermann-Mauguin Symbol"), get("Hall Symbol"));
            // A magnetic phase given by its operators alone has no symbol (runs/NdNi-HTO, phase 2): left out.
            if hm == "unknown" {
                blocks -= 1;
                continue;
            }
            groups.insert(hm.to_string());
            let tag = format!("{name}: {hm}");
            let number: u16 = get("Number of Space group").parse().unwrap();
            // The setting choice is printed only for non-standard settings; read it within the block.
            // Two forms: the reduced list ("List of S.O. without inversion and lattice centring translations")
            // or, in LaMnO3_cub.out, every operator with its symbol ("List of all Symmetry Operators").
            let end = start
                + lines[start..]
                    .iter()
                    .position(|l| {
                        l.contains("List of S.O.") || l.contains("List of all Symmetry Operators")
                    })
                    .unwrap_or_else(|| panic!("{tag}: no list of operators"));
            let all_listed = lines[end].contains("List of all");
            let choice = field(&lines[..end], start, "Table Setting Choice").unwrap_or("");

            // FullProf's symbol read as FullProf reads it (L030): "F d -3 m" and "I 41/a m d" are origin
            // choice 2 there, where the International Tables' order puts choice 1 first.
            match SpaceGroup::from_fullprof_symbol(hm).map(|g| g.setting) {
                // The same group may be written by two Hall symbols that differ by a centring translation ("B m a b":
                // the table's '-B 2ab 2', FullProf's '-B 2bc 2', (½,½,0) + (½,0,½) ≡ (0,½,½)): the groups are compared.
                Ok(Some(s)) if same_group(s.hall, hall) && s.number == number && (choice.is_empty() || s.choice == choice) => {}
                Ok(Some(s)) => fails.push(format!("{tag}: the core reads '{}' no. {} choice '{}', FullProf '{hall}' no. {number} choice '{choice}'", s.hall, s.number, s.choice)),
                Ok(None) => fails.push(format!("{tag}: no table setting for the symbol")),
                Err(e) => fails.push(format!("{tag}: {e}")),
            }
            let g = match SpaceGroup::from_hall(hall) {
                Ok(g) => g,
                Err(e) => {
                    fails.push(format!("{tag}: {e}"));
                    continue;
                }
            };
            // FullProf writes a trigonal Laue class with its orientation on the hexagonal axes, −3m1 or −31m (alum.out,
            // R −3 c: "Laue Class: -3m1", and "Laue symmetry -3m1 will be used to generate HKL"; L341): the class is −3m,
            // as the core names it.
            let laue = match get("Laue Class") {
                "-3m1" | "-31m" => "-3m",
                l => l,
            };
            let fp = (
                get("Crystal System").to_lowercase(),
                laue.to_string(),
                get("Point Group").to_string(),
            );
            if fp
                != (
                    g.crystal_system.to_string(),
                    g.laue.to_string(),
                    g.point_group.to_string(),
                )
            {
                fails.push(format!(
                    "{tag}: FullProf {fp:?}, core ({}, {}, {})",
                    g.crystal_system, g.laue, g.point_group
                ));
            }
            let mult: usize = get("General multiplicity").parse().unwrap();
            if mult != g.multiplicity() {
                fails.push(format!(
                    "{tag}: multiplicity {mult} against {}",
                    g.multiplicity()
                ));
            }
            let centric = get("Centrosymmetry");
            let at_origin = centric.contains("-1 at origin");
            if centric.starts_with("Centric") != g.centrosymmetric
                || at_origin != g.inversion_at_origin
            {
                fails.push(format!(
                    "{tag}: '{centric}' against centrosymmetric {}, at origin {}",
                    g.centrosymmetric, g.inversion_at_origin
                ));
            }

            // FullProf's operators, then the inversion and the centring, as a full set.
            let mut reduced = Vec::new();
            for l in lines[end + 1..].iter().skip_while(|l| l.trim().is_empty()) {
                if !l.contains("SYMM(") {
                    break;
                }
                for part in l.split("SYMM(").skip(1) {
                    let xyz = part.split_once(':').unwrap().1.split("=>").next().unwrap();
                    let xyz = xyz.split("Symbol:").next().unwrap().trim();
                    reduced.push(SymOp::parse_xyz(xyz).unwrap_or_else(|e| panic!("{tag}: {e}")));
                }
            }
            let reduced_n: usize = if all_listed {
                mult
            } else {
                get("Reduced Number of S.O.").parse().unwrap()
            };
            if reduced.len() != reduced_n {
                fails.push(format!(
                    "{tag}: read {} SYMM lines, FullProf says {reduced_n}",
                    reduced.len()
                ));
            }
            let inv = SymOp::new([[-1, 0, 0], [0, -1, 0], [0, 0, -1]], [0, 0, 0]);
            let mut full = BTreeSet::new();
            for op in reduced.iter().filter(|_| !all_listed) {
                for c in &g.centring {
                    let t = SymOp::translation(*c);
                    full.insert(t.compose(op));
                    if at_origin {
                        full.insert(t.compose(&inv.compose(op)));
                    }
                }
            }
            if all_listed {
                full = reduced.iter().copied().collect();
            }
            let mine: BTreeSet<SymOp> = g.ops.iter().copied().collect();
            if full != mine {
                fails.push(format!(
                    "{tag}: FullProf's {} operators are not the core's {}",
                    full.len(),
                    mine.len()
                ));
            }
        }
    }
    assert!(fails.is_empty(), "{}", fails.join("\n"));
    // At least one block per file; the groups of the 13 tutorial files must all be among those read.
    assert!(
        blocks >= files.len(),
        "space-group blocks read: {blocks} for {} files",
        files.len()
    );
    let groups: Vec<_> = groups.into_iter().collect();
    println!(
        "{blocks} blocks in {} files, groups: {}",
        files.len(),
        groups.join(", ")
    );
    for g in ["P 2/m", "P 21/a", "P 21/c", "P n m a"] {
        assert!(
            groups.iter().any(|x| x == g),
            "{g} missing from the groups read"
        );
    }
}
