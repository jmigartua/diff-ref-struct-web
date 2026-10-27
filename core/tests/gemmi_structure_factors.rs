//! Structure factors against gemmi (test `gemmi_structure_factors`).
//!
//! `tests/data/gemmi_structure_factors.txt` holds the reference structures of L020 as their tutorial `.pcr` files give
//! them (Y2O3, MgAl2O4, Co-doped YBCO, and the phases of the two quantitative-analysis tutorials: rutile, anatase,
//! laurionite, phosgenite, cerussite, galena, gypsum), with the structure factors gemmi computed for them (X-rays
//! without anomalous terms, and neutrons), gemmi's form factors and scattering lengths, written by
//! `tests/data/make_gemmi_reference.py`. gemmi builds the group from the Hall symbol and places the atoms with their
//! chemical occupancy, on its own; the core reads the `.pcr` symbol (`from_fullprof_symbol`) and FullProf's Occ.
//!
//! The test separates the summation from the tables:
//! 1. With gemmi's coefficients, the core's F must equal gemmi's, X-rays and neutrons, to 10⁻⁶ of Σ_j Occ_j·M·|f_j|
//!    (gemmi stores its form-factor coefficients in single precision; the agreement found is 6·10⁻⁸ and 4·10⁻¹⁰); the unique reflections with d ≥ 0.9 Å must be the same set
//!    (up to equivalence), and the group the same operators.
//! 2. The core's own tables (CrysFML's) are compared with gemmi's, entry by entry; the entries that differ are listed,
//!    with the largest effect they have on |F| for each phase. For the 15 species of these structures the two X-ray
//!    tables agree to gemmi's single precision except O²⁻ ("O-2" of rutile and anatase): CrysFML's entry is a
//!    two-Gaussian fit (a = 4.758, 3.637; b = 7.831, 30.05; c = 1.594), gemmi's a four-Gaussian one; where each comes
//!    from was not checked. Its effect on |F| is at most 7·10⁻⁴ of Σ Occ·M·|f|. The neutron lengths of all 14 elements
//!    agree (CrysFML's `Chem_Info` in 10⁻¹² cm, gemmi's Neutron92 in fm).

use drs::cell::Cell;
use drs::reflections::{equivalents, generate};
use drs::scattering::{Radiation, Scatterer};
use drs::spacegroup::SpaceGroup;
use drs::structure_factor::{Atom, Structure};
use std::collections::BTreeSet;

struct Phase {
    name: String,
    structure: Structure,
    hall: String,
    xray: Vec<(String, Option<Vec<f64>>)>,
    neutron: Vec<(String, f64)>,
    refs: Vec<([i32; 3], [f64; 4])>,
}

fn read() -> Vec<Phase> {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/gemmi_structure_factors.txt");
    let text = std::fs::read_to_string(path).unwrap();
    let mut phases: Vec<Phase> = Vec::new();
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = line.split_whitespace().collect();
        match f[0] {
            "PHASE" => {
                let q: Vec<&str> = line.split('"').collect();
                let c: Vec<f64> = q[4].split_whitespace().take(6).map(|x| x.parse().unwrap()).collect();
                phases.push(Phase {
                    name: f[1].to_string(),
                    structure: Structure {
                        group: SpaceGroup::from_fullprof_symbol(q[1]).unwrap(),
                        cell: Cell::new(c[0], c[1], c[2], c[3], c[4], c[5]).unwrap(),
                        atoms: Vec::new(),
                    },
                    hall: q[3].to_string(),
                    xray: Vec::new(),
                    neutron: Vec::new(),
                    refs: Vec::new(),
                });
            }
            "ATOM" => {
                let v: Vec<f64> = f[3..8].iter().map(|x| x.parse().unwrap()).collect();
                phases.last_mut().unwrap().structure.atoms.push(Atom {
                    label: f[1].to_string(),
                    species: f[2].to_string(),
                    x: [v[0], v[1], v[2]],
                    biso: v[3],
                    occ: v[4],
                    beta: None,
                    charge: None,
                });
            }
            "XRAY" => {
                let c = if f[2] == "none" { None } else { Some(f[2..11].iter().map(|x| x.parse().unwrap()).collect()) };
                phases.last_mut().unwrap().xray.push((f[1].to_string(), c));
            }
            "NEUTRON" => phases.last_mut().unwrap().neutron.push((f[1].to_string(), f[2].parse().unwrap())),
            "REF" => {
                let h = [f[1].parse().unwrap(), f[2].parse().unwrap(), f[3].parse().unwrap()];
                let v: Vec<f64> = f[4..8].iter().map(|x| x.parse().unwrap()).collect();
                phases.last_mut().unwrap().refs.push((h, [v[0], v[1], v[2], v[3]]));
            }
            _ => panic!("unknown line {line}"),
        }
    }
    phases
}

#[test]
fn gemmi_structure_factors() {
    let phases = read();
    assert_eq!(phases.len(), 10);
    let mut fails = Vec::new();
    let mut table_diffs: BTreeSet<String> = BTreeSet::new();
    let mut n_refs = 0;
    for p in &phases {
        let st = &p.structure;
        let n = st.group.ops.len() as f64;
        // The same group as gemmi's Hall symbol.
        let hall = SpaceGroup::from_hall(&p.hall).unwrap();
        let a: BTreeSet<_> = st.group.ops.iter().map(|o| (o.r, o.t)).collect();
        let b: BTreeSet<_> = hall.ops.iter().map(|o| (o.r, o.t)).collect();
        if a != b {
            fails.push(format!("{}: the .pcr symbol does not give the operators of {}", p.name, p.hall));
        }
        // gemmi's scatterers, atom by atom; neutron lengths from fm to 10⁻¹² cm.
        let gx: Vec<Scatterer> = st
            .atoms
            .iter()
            .map(|at| {
                let c = p.xray.iter().find(|x| x.0 == at.species).unwrap().1.clone().expect("gemmi has the species");
                Scatterer::Xray { symbol: at.species.clone(), a: [c[0], c[1], c[2], c[3]], b: [c[4], c[5], c[6], c[7]], c: c[8], fp: 0.0, fpp: 0.0 }
            })
            .collect();
        let gn: Vec<Scatterer> = st
            .atoms
            .iter()
            .map(|at| Scatterer::Neutron { symbol: at.species.clone(), b_imaginary: 0.0, b: p.neutron.iter().find(|x| x.0 == at.species).unwrap().1 / 10.0 })
            .collect();
        let cx = st.scatterers(Radiation::Xray { lambda: 1.5406, anomalous: false }).unwrap();
        let cn = st.scatterers(Radiation::Neutron).unwrap();
        // The tables, entry by entry.
        for (sp, c) in &p.xray {
            let Scatterer::Xray { a, b, c: c0, .. } = Scatterer::lookup(sp, Radiation::Xray { lambda: 1.5406, anomalous: false }).unwrap() else { panic!() };
            let mine = [a[0], a[1], a[2], a[3], b[0], b[1], b[2], b[3], c0];
            let g = c.as_ref().unwrap();
            if mine.iter().zip(g).any(|(m, v)| (m - v).abs() > 1e-6 * v.abs().max(1.0)) {
                table_diffs.insert(format!("X-ray {sp}: CrysFML {mine:?}, gemmi {g:?}"));
            }
        }
        for (sp, bg) in &p.neutron {
            let Scatterer::Neutron { b, .. } = Scatterer::lookup(sp, Radiation::Neutron).unwrap() else { panic!() };
            if (b * 10.0 - bg).abs() > 1e-6 {
                table_diffs.insert(format!("neutron {}: CrysFML {b} × 10⁻¹² cm, gemmi {bg} fm", drs::scattering::element_of(sp)));
            }
        }
        // The reflections: the same unique set for d ≥ 0.9 Å.
        let mine = generate(&st.group, &st.cell, 0.9 - 1e-9, 1e9, None);
        if mine.len() != p.refs.len() {
            fails.push(format!("{}: {} unique reflections with d ≥ 0.9 Å, gemmi {}", p.name, mine.len(), p.refs.len()));
        }
        let reps: BTreeSet<[i32; 3]> = mine.iter().map(|r| r.hkl).collect();
        // Scale of F for the tolerance: Σ Occ·M·|f| at s = 0.
        let scale = |sc: &[Scatterer]| -> f64 { st.atoms.iter().zip(sc).map(|(a, s)| a.occ * n * s.at(0.0).abs()).sum() };
        let (sx, sn) = (scale(&gx), scale(&gn));
        let (mut worst_x, mut worst_n, mut table_x, mut table_n) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        for (h, v) in &p.refs {
            n_refs += 1;
            if equivalents(&st.group, *h).intersection(&reps).count() != 1 {
                fails.push(format!("{}: gemmi's {h:?} is not equivalent to exactly one of the core's reflections", p.name));
            }
            let fx = st.structure_factor(&gx, *h);
            let fnn = st.structure_factor(&gn, *h);
            let dx = (fx.re - v[0]).hypot(fx.im - v[1]) / sx;
            let dn = (fnn.re - v[2] / 10.0).hypot(fnn.im - v[3] / 10.0) / sn;
            worst_x = worst_x.max(dx);
            worst_n = worst_n.max(dn);
            if dx > 1e-6 || dn > 1e-6 {
                fails.push(format!("{}: {h:?} core X-ray {fx}, gemmi {}{:+}i; neutron {fnn}, gemmi/10 {}{:+}i", p.name, v[0], v[1], v[2] / 10.0, v[3] / 10.0));
            }
            // The effect of the tables: the core's own against gemmi's, relative to the scale.
            table_x = table_x.max((st.structure_factor(&cx, *h).abs() - fx.abs()).abs() / sx);
            table_n = table_n.max((st.structure_factor(&cn, *h).abs() - fnn.abs()).abs() / sn);
        }
        println!(
            "{:<11} {:>4} reflections: same tables {:.1e} (X) {:.1e} (n); own tables, largest |Δ|F|| {:.1e} (X) {:.1e} (n) of Σ Occ·M·|f|",
            p.name, p.refs.len(), worst_x, worst_n, table_x, table_n
        );
    }
    println!("table entries that differ ({}):\n  {}", table_diffs.len(), table_diffs.iter().cloned().collect::<Vec<_>>().join("\n  "));
    let expected = ["X-ray O-2"];
    let found: Vec<&str> = table_diffs.iter().map(|d| d.split(':').next().unwrap()).collect();
    assert_eq!(found, expected, "table entries that differ between CrysFML and gemmi");
    let shown: Vec<_> = fails.iter().take(40).cloned().collect();
    assert!(fails.is_empty(), "{} disagreements:\n{}", fails.len(), shown.join("\n"));
    assert_eq!(n_refs, 1435, "reflections compared");
}
