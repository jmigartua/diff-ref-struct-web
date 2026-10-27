//! ISODISTORT's mode CIF (P13, `isodistort`).
//! - **Case A as ISODISTORT 6.12.2 wrote it** (`manuals/bcs/isodistort/`, made by the merging session, L275; L278):
//!   12 displacive modes in 7 irreps; the coordinate deltas from the modes, δ = M·(A × N), to the file's print; the
//!   amplitudes back from the deltas; and the core's own decomposition of ISODISTORT's input structures
//!   (`mode_decomposition::decompose`, basis (1/2,1/2,0), (−1/2,1/2,0), (0,0,1)): each irrep's amplitude, each orbit's
//!   and each mode's |amplitude| as ISODISTORT gives them, to its five decimals (GM4+ 0.82937 Å, the core 0.829371),
//!   and each irrep's dimension as the components of ISODISTORT's order-parameter direction.
//! - **The fixture** (`tests/data/isodistort_fixture.cif`, L263): the reader's mechanics (the loops by ID, the label's
//!   parts, the norms, the deltas and formulas), and a CIF without modes refused.

use drs::isodistort::{isodistort_cif, split_label};
use drs::mode_decomposition as md;

#[test]
fn fixture() {
    let text = std::fs::read_to_string(format!("{}/tests/data/isodistort_fixture.cif", env!("CARGO_MANIFEST_DIR"))).unwrap();
    let iso = isodistort_cif(&text).unwrap();
    assert_eq!(iso.modes.len(), 2);
    let m = &iso.modes[0];
    assert_eq!((m.id, m.parent.as_str(), m.k.as_str(), m.irrep.as_str(), m.direction.as_str()), (1, "Pm-3m", "[1/2,1/2,1/2]", "R4+", "(a,0,0)"));
    assert_eq!((m.atom.as_str(), m.wyckoff.as_str(), m.kind.as_str(), m.point_irrep.as_str()), ("O1", "d", "dsp", "Eu"));
    assert_eq!((m.value, m.norm), (0.1404, Some(0.0641)));
    assert_eq!(iso.modes[1].norm, Some(0.04533));
    assert_eq!(iso.irreps(), vec![("[1/2,1/2,1/2]R4+".to_string(), vec![0]), ("[0,0,0]GM1+".to_string(), vec![1])]);
    assert_eq!(iso.deltas, vec![(1, "O2_dx".to_string(), -0.009), (2, "O2_dy".to_string(), -0.009)]);
    assert_eq!(iso.formulas[0], ("O2_x".to_string(), "1/4 + O2_dx".to_string()));
    assert_eq!(iso.matrix.len(), 3);
    // A label without the parts, and a CIF without modes.
    assert!(split_label("GM1+").is_none());
    assert!(isodistort_cif("data_x\n_cell_length_a 5\n").unwrap_err().contains("not an ISODISTORT mode CIF"));
}

fn iso_dir() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/bcs/isodistort")
}

#[test]
fn case_a_as_isodistort_wrote_it() {
    let text = std::fs::read_to_string(iso_dir().join("caseA_Sr2ScSbO6_P21n_modes.cif")).unwrap();
    let iso = isodistort_cif(&text).unwrap();
    assert_eq!((iso.modes.len(), iso.irreps().len(), iso.deltas.len(), iso.formulas.len(), iso.matrix.len()), (12, 7, 12, 18, 34));
    assert!(iso.modes.iter().all(|m| m.norm.is_some() && m.parent == "Fm-3m" && m.kind == "dsp"));
    let gm4 = iso.modes.iter().find(|m| m.irrep == "GM4+").unwrap();
    assert_eq!((gm4.k.as_str(), gm4.direction.as_str(), gm4.atom.as_str(), gm4.wyckoff.as_str(), gm4.point_irrep.as_str(), gm4.value), ("[0,0,0]", "(a,-a,0)", "O1", "e", "E", 0.82937));

    // The deltas from the modes, to the print (amplitudes and norms to 5 decimals: ±2·10⁻⁵ at most).
    let d = iso.deltas_from_modes();
    let mut worst = 0.0f64;
    for (id, label, v) in &iso.deltas {
        let got = d.iter().find(|x| x.0 == *id).map_or(0.0, |x| x.1);
        worst = worst.max((got - v).abs());
        assert!((got - v).abs() <= 2e-5, "{label}: {got:.6} from the modes, {v} in the file");
    }
    // And back: the amplitudes from the deltas (the deltas' print, 5·10⁻⁶, over a norm of 0.036 or more).
    let a = iso.modes_from_deltas().unwrap();
    let mut worst_a = 0.0f64;
    for (id, v) in &a {
        let m = iso.modes.iter().find(|m| m.id == *id).unwrap();
        worst_a = worst_a.max((v - m.value).abs());
        assert!((v - m.value).abs() <= 3e-4, "{}: {v:.5} from the deltas, {} in the file", m.label, m.value);
    }
    println!("ISODISTORT's case A: deltas from the modes within {worst:.1e}, amplitudes from the deltas within {worst_a:.1e}");

    // The core's decomposition of ISODISTORT's own inputs.
    let read = |f: &str| md::from_cif(&std::fs::read_to_string(iso_dir().join(f)).unwrap()).unwrap();
    let (parent, distorted) = (read("input_caseA_parent_Fm-3m.cif"), read("input_caseA_distorted_P21n.cif"));
    let p = md::parse_basis("1/2a+1/2b,-1/2a+1/2b,c").unwrap();
    let dec = md::decompose(&parent, &distorted, &p, [0.0; 3]).unwrap();
    assert_eq!((dec.irreps.len(), dec.irreps.iter().map(|i| i.modes).sum::<usize>()), (7, 12));
    let mut used = vec![false; dec.irreps.len()];
    for (key, idx) in iso.irreps() {
        let modes: Vec<&drs::isodistort::IsoMode> = idx.iter().map(|&i| &iso.modes[i]).collect();
        let amp = modes.iter().map(|m| m.value * m.value).sum::<f64>().sqrt();
        let gamma = modes[0].k == "[0,0,0]";
        // The core's irrep of the same kind of k (Γ or not), the same number of modes and the nearest amplitude.
        let (j, ours) = dec
            .irreps
            .iter()
            .enumerate()
            .filter(|(j, ir)| !used[*j] && (ir.star == vec![[0.0; 3]]) == gamma && ir.modes == modes.len())
            .min_by(|a, b| (a.1.amplitude - amp).abs().total_cmp(&(b.1.amplitude - amp).abs()))
            .unwrap_or_else(|| panic!("{key}: no irrep of the core with {} modes", modes.len()));
        used[j] = true;
        let tol = 5e-6 * (modes.len() as f64).sqrt() + 2e-6;
        println!("  {key:12} {} modes  ISODISTORT {amp:.5} Å  the core's {} {:.6} Å", modes.len(), ours.name, ours.amplitude);
        assert!((ours.amplitude - amp).abs() <= tol, "{key}: {} against {amp}", ours.amplitude);
        // The irrep's dimension: the components of ISODISTORT's order-parameter direction, (a,-a,0) 3, (0,0;0,0;a,a) 6.
        let dim = modes[0].direction.trim_matches(|c| c == '(' || c == ')').split([',', ';']).count();
        assert_eq!(ours.dimension, dim, "{key}{}", modes[0].direction);
        // Orbit by orbit, and the modes' |amplitudes|.
        for o in &ours.orbits {
            let label = &parent.atoms[o.orbit].label;
            let mut theirs: Vec<f64> = modes.iter().filter(|m| &m.atom == label).map(|m| m.value.abs()).collect();
            let mut mine: Vec<f64> = o.modes.iter().map(|m| m.amplitude.abs()).collect();
            theirs.sort_by(f64::total_cmp);
            mine.sort_by(f64::total_cmp);
            assert_eq!(mine.len(), theirs.len(), "{key} on {label}");
            for (x, y) in mine.iter().zip(&theirs) {
                assert!((x - y).abs() <= 5.1e-6, "{key} on {label}: |amplitude| {x:.6} against {y}");
            }
        }
    }
}
