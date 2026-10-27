//! Reading structures from CIF files (test `cif_structures`).
//!
//! - The CIFs FullProf wrote for PPH3 (`PPH3-sls.cif`) and LaMnO3 (`Lmo_san.cif`) must give the structures of its `.out`
//!   files, as copied in `tests/data/fullprof_structure_factors.txt`: the same group, coordinates to the CIF's printed
//!   precision, Biso = 8π²·Uiso to the precision of the printed U, and FullProf's Occ from the CIF's chemical occupancy
//!   (LaMnO3: 0.5, 0.5, 1, 0.5 for the sites 4c, 4b, 8d, 4c of P n m a).
//! - An ICSD CIF (galena, `manuals/fullprof/5-Quantitative_Phase_Analysis/EgyptianMakeup/ICSD/galena.cif`, Latin-1, with
//!   symmetry operators and ionic type symbols) must give F m -3 m from its 192 operators, the species "Pb+2" and "S-2",
//!   and Occ = 4/192.

use drs::cif::structure;
use drs::spacegroup::SpaceGroup;
use std::collections::BTreeSet;
use std::f64::consts::PI;

fn root(rel: &str) -> String {
    format!("{}/../{rel}", env!("CARGO_MANIFEST_DIR"))
}

fn latin1(path: &str) -> String {
    std::fs::read(path).unwrap().iter().map(|&b| b as char).collect()
}

/// The ATOM lines of a phase of fullprof_structure_factors.txt: label, x, Biso, Occ.
fn fullprof_atoms(phase: &str) -> Vec<(String, [f64; 3], f64, f64)> {
    let text = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/fullprof_structure_factors.txt")).unwrap();
    let mut on = false;
    let mut out = Vec::new();
    for l in text.lines() {
        let f: Vec<&str> = l.split_whitespace().collect();
        match f.first() {
            Some(&"PHASE") => on = f[1] == phase,
            Some(&"ATOM") if on => {
                let v: Vec<f64> = f[3..8].iter().map(|x| x.parse().unwrap()).collect();
                out.push((f[1].to_string(), [v[0], v[1], v[2]], v[3], v[4]));
            }
            _ => {}
        }
    }
    out
}

fn same_ops(a: &SpaceGroup, hall: &str) -> bool {
    let b = SpaceGroup::from_hall(hall).unwrap();
    let x: BTreeSet<_> = a.ops.iter().map(|o| (o.r, o.t)).collect();
    let y: BTreeSet<_> = b.ops.iter().map(|o| (o.r, o.t)).collect();
    x == y
}

#[test]
fn cif_structures() {
    let cases = [
        ("PPH3", "manuals/fullprof/9-SAnn-OrganicCompounds/PPH3/Rietveld/PPH3-sls.cif", "-P 2ybc"),
        ("LaMnO3", "manuals/fullprof/7-FullProf-Tutorial-SymmetryModes/LaMnO3/SimulatedAnnealing/Lmo_san.cif", "-P 2ac 2n"),
    ];
    for (phase, file, hall) in cases {
        let st = structure(&latin1(&root(file))).unwrap();
        assert!(same_ops(&st.group, hall), "{file}: group");
        let fp = fullprof_atoms(phase);
        assert_eq!(st.atoms.len(), fp.len(), "{file}: atoms");
        for (a, (label, x, biso, occ)) in st.atoms.iter().zip(&fp) {
            assert_eq!(&a.label, label);
            for (k, (u, v)) in a.x.iter().zip(x).enumerate() {
                // The CIF prints 4 or 5 decimals.
                assert!((u - v).abs() <= 0.5e-4 + 1e-12, "{file}: {label} x{k} {u} against {v}");
            }
            // Uiso printed to five decimals: B to 8π²·0.5·10⁻⁵.
            assert!((a.biso - biso).abs() <= 8.0 * PI * PI * 0.5e-5 + 0.5e-5, "{file}: {label} Biso {} against {biso}", a.biso);
            assert!((a.occ - occ).abs() < 1e-12, "{file}: {label} Occ {} against {occ}", a.occ);
        }
    }
    let st = structure(&latin1(&root("manuals/fullprof/5-Quantitative_Phase_Analysis/EgyptianMakeup/ICSD/galena.cif"))).unwrap();
    assert!(same_ops(&st.group, "-F 4 2 3"), "galena: F m -3 m");
    let species: Vec<&str> = st.atoms.iter().map(|a| a.species.as_str()).collect();
    assert_eq!(species, ["Pb+2", "S-2"]);
    for a in &st.atoms {
        assert!((a.occ - 4.0 / 192.0).abs() < 1e-12);
    }
}
