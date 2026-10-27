//! Symmetry modes computed from a parent and a distorted structure (P13, `mode_decomposition`), against AMPLIMODES
//! on the course's files (`manuals/bcs/mode-crystallography-manual/course-um6p-2024/04_tools_materials/0_1…0_4`):
//!
//! - **The summaries** (`amplimodes_summaries`): the 10 analyses that give both structures and the transformation.
//!   From the structures alone, the irreps the distorted group keeps are AMPLIMODES's, one to one: the same k (its
//!   k-vector in the computed star), parity, dimension (the length of its direction) and number of modes, on each
//!   orbit as its "Symmetry Modes Summary" lists them. Every amplitude and the global distortion agree to the 4
//!   decimals printed, in 54 of 56 irreps; the 2 others are analyses 1 and 2 of `0_4`, GM5+, 0.0229 against 0.0226
//!   and 0.0224 against 0.0221 (`FROM_THE_INPUT`): AMPLIMODES moved the input's Sr from x = −0.00003 to 0 and z
//!   from 0.2492 to 0.24923 (its "Atom Mappings"). From the positions of its atom mappings, all 56 agree.
//! - **The modes** (`amplimodes_modes`): analysis 2 of `0_1` decomposed and its irreps named by the 12 modes of
//!   analysis 3's block for FullProf (same parent and transformation): each mode has norm 1 and lies in the computed
//!   modes of its irrep and orbit; each amplitude of A_MODES is the computed one (up to the sign, a convention) where
//!   the irrep has one mode on the orbit, and their norm is the computed one where it has two.
//! - **The block for FullProf** (`fullprof_block`): the block the core writes, read back by the reader of AMPLIMODES's
//!   block, names the same irreps of a new decomposition and gives the amplitudes to its 6 decimals.
//! - **For refinement** (`modes_for_refinement`): the computed modes as `Modes` put every atom back where the
//!   distorted structure has it, and `Modes::project` (the refinement's own code) finds the same amplitudes.
//! - **A supercell** (`supercell_r_point`): SrTiO₃'s I4/mcm tilt from Pm-3m, a cell P⁻¹ with halves (n = 2): one
//!   mode at k = (1/2,1/2,1/2), of dimension 3, of amplitude 4δa (the closed form); odd with Sr at the parent's
//!   origin and even with Ti there, the inversion moving by (1/2,1/2,1/2) and e^{2πik·(1,1,1)} = −1 (the labels
//!   R5− and R4+ that the literature gives the two settings are believed, not checked here: no table of irreps is
//!   in `manuals/`).
//! - **Refused** (`refused`): an origin that makes the distorted group no subgroup; a distorted structure that does
//!   not fill the reference's sites.

use drs::mode_decomposition::{decompose, from_amplimodes, Atom, Decomposition, Structure};
use drs::modes::{read_amplimodes, Analysis};

fn files() -> Vec<(String, String)> {
    let d = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/bcs/mode-crystallography-manual/course-um6p-2024/04_tools_materials");
    let mut v: Vec<(String, String)> = std::fs::read_dir(&d)
        .unwrap()
        .filter_map(|e| {
            let p = e.ok()?.path();
            let name = p.file_name()?.to_str()?.to_string();
            (name.starts_with("0_") && name.ends_with(".txt")).then(|| (name, std::fs::read(&p).unwrap().iter().map(|&b| b as char).collect()))
        })
        .collect();
    v.sort();
    v
}

/// The "Symmetry Modes Summary" of each analysis of a text: (orbit label, irrep, count).
fn orbit_counts(text: &str) -> Vec<Vec<(String, String, usize)>> {
    text.split("Symmetry mode analysis")
        .skip(1)
        .map(|seg| {
            let ls: Vec<&str> = seg.lines().collect();
            let Some(i) = ls.iter().position(|l| l.trim() == "Symmetry Modes Summary") else { return Vec::new() };
            let mut v = Vec::new();
            for l in &ls[i + 2..] {
                let f: Vec<&str> = l.split_whitespace().collect();
                if f.len() < 3 || !f[2].contains('(') {
                    break;
                }
                for m in &f[2..] {
                    if let Some((ir, c)) = m.split_once('(') {
                        v.push((f[0].to_string(), ir.to_string(), c.trim_end_matches(')').parse().unwrap()));
                    }
                }
            }
            v
        })
        .collect()
}

fn parse_k(s: &str) -> [f64; 3] {
    let v: Vec<f64> = s.trim_matches(|c| c == '(' || c == ')').split(',').map(|t| t.trim().parse().unwrap()).collect();
    [v[0], v[1], v[2]]
}

/// Each AMPLIMODES summary line matched to one computed irrep: the same k, parity, dimension and number of modes; among
/// several such, the nearest in amplitude.
fn assign(d: &Decomposition, a: &Analysis) -> Vec<usize> {
    let mut taken = vec![false; d.irreps.len()];
    let mut out = Vec::new();
    for s in &a.summary {
        let k = parse_k(&s.k);
        let parity = if s.irrep.ends_with('+') { Some(1) } else if s.irrep.ends_with('-') { Some(-1) } else { None };
        let dim = s.direction.split(',').count();
        let cands: Vec<usize> = (0..d.irreps.len())
            .filter(|&r| {
                let ir = &d.irreps[r];
                !taken[r] && ir.star.iter().any(|q| (0..3).all(|c| ((q[c] - k[c]) / 2.0 - ((q[c] - k[c]) / 2.0).round()).abs() < 1e-9 || (q[c] - k[c]).abs() < 1e-9 || ((q[c] - k[c]) - (q[c] - k[c]).round()).abs() < 1e-9)) && ir.parity == parity && ir.dimension == dim && ir.modes == s.dimension
            })
            .collect();
        let r = *cands.iter().min_by(|&&x, &&y| (d.irreps[x].amplitude - s.amplitude).abs().total_cmp(&(d.irreps[y].amplitude - s.amplitude).abs())).unwrap_or_else(|| panic!("{} {}: no computed irrep with its k, parity, dimension and modes", s.k, s.irrep));
        taken[r] = true;
        out.push(r);
    }
    assert!(taken.iter().all(|&t| t), "computed irreps AMPLIMODES does not list: {:?}", d.irreps.iter().zip(&taken).filter(|(_, t)| !**t).map(|(r, _)| &r.name).collect::<Vec<_>>());
    out
}

/// Irreps whose amplitude from the input differs from AMPLIMODES's printed one, as measured in L263: (file, analysis,
/// irrep, computed, AMPLIMODES's). AMPLIMODES moved the input's Sr (its atom mappings).
const FROM_THE_INPUT: [(&str, usize, &str, f64, f64); 2] = [("0_4_co_te_data_012.txt", 1, "GM5+", 0.0229, 0.0226), ("0_4_co_te_data_012.txt", 2, "GM5+", 0.0224, 0.0221)];

/// Half a unit of the 4th decimal, and the rounding of exact halves (0.503650 printed 0.5036).
const PRINTED: f64 = 0.00005 + 1e-6;

#[test]
fn amplimodes_summaries() {
    let (mut analyses, mut irreps, mut agree, mut agree_mapped) = (0, 0, 0, 0);
    let mut apart = Vec::new();
    for (name, text) in files() {
        let counts = orbit_counts(&text);
        for (k, a) in read_amplimodes(&text).iter().enumerate() {
            let Ok((parent, distorted, p, o)) = from_amplimodes(a) else { continue };
            analyses += 1;
            let d = decompose(&parent, &distorted, &p, o).unwrap_or_else(|e| panic!("{name} {}: {e}", k + 1));
            assert!(d.residual < 1e-9 && d.translation.iter().all(|t| t.abs() < 1e-9), "{name} {}", k + 1);
            assert!((d.global - a.global.unwrap()).abs() <= PRINTED, "{name} {}: global {} against {:?}", k + 1, d.global, a.global);
            let map = assign(&d, a);
            // From the positions of AMPLIMODES's atom mappings.
            let mut dm = distorted.clone();
            dm.atoms = a.mappings.iter().map(|m| Atom { label: m.2.clone(), species: m.2.chars().take_while(|c| c.is_ascii_alphabetic()).collect(), position: m.3 }).collect();
            let dmap = decompose(&parent, &dm, &p, o).unwrap();
            let map2 = assign(&dmap, a);
            for (j, s) in a.summary.iter().enumerate() {
                irreps += 1;
                let r = &d.irreps[map[j]];
                if (r.amplitude - s.amplitude).abs() <= PRINTED {
                    agree += 1;
                } else {
                    apart.push((name.clone(), k + 1, s.irrep.clone(), (r.amplitude * 1e4).round() / 1e4, s.amplitude));
                }
                if (dmap.irreps[map2[j]].amplitude - s.amplitude).abs() <= PRINTED {
                    agree_mapped += 1;
                }
                // The modes on each orbit, as its "Symmetry Modes Summary".
                for (orbit, ir, c) in counts[k].iter().filter(|x| x.1 == s.irrep) {
                    let o = parent.atoms.iter().position(|at| at.label == *orbit).unwrap();
                    let got = r.orbits.iter().find(|x| x.orbit == o).map_or(0, |x| x.modes.len());
                    assert_eq!(got, *c, "{name} {}: {ir} on {orbit}", k + 1);
                }
            }
            println!("{name} {}: {} irreps, global {:.4}", k + 1, a.summary.len(), d.global);
        }
    }
    println!("{analyses} analyses, {irreps} irreps: {agree} agree from the input, {agree_mapped} from the mappings; apart {apart:?}");
    assert_eq!(analyses, 10);
    assert_eq!(irreps, 56);
    assert_eq!(agree_mapped, irreps);
    assert_eq!(apart.len(), FROM_THE_INPUT.len());
    for (x, y) in apart.iter().zip(FROM_THE_INPUT) {
        assert_eq!((x.0.as_str(), x.1, x.2.as_str()), (y.0, y.1, y.2));
        assert!((x.3 - y.3).abs() < 1e-9 && (x.4 - y.4).abs() < 1e-9, "{x:?}");
    }
}

fn analysis(file: &str, k: usize) -> (Analysis, Vec<Analysis>) {
    let (_, text) = files().into_iter().find(|f| f.0 == file).unwrap();
    let all = read_amplimodes(&text);
    (all[k - 1].clone(), all)
}

#[test]
fn amplimodes_modes() {
    let (a, all) = analysis("0_1_P21nToCubicTransData.txt", 2);
    let block = &all[2];
    assert_eq!(block.transformation, a.transformation);
    let (parent, distorted, p, o) = from_amplimodes(&a).unwrap();
    let mut d = decompose(&parent, &distorted, &p, o).unwrap();
    let modes = block.modes().unwrap();
    assert_eq!(modes.modes.len(), 12);
    assert_eq!(d.name_by_modes(&modes), 7);
    let mut names: Vec<&str> = d.irreps.iter().map(|r| r.name.as_str()).collect();
    names.sort();
    assert_eq!(names, ["GM1+", "GM3+", "GM4+", "GM5+", "X2+", "X3+", "X5+"]);
    // Each mode: norm 1, in its irrep's modes on its orbit.
    let mut blocks: std::collections::BTreeMap<(String, usize), Vec<f64>> = Default::default();
    for m in &modes.modes {
        let atoms: Vec<([f64; 3], [f64; 3])> = m.vectors.iter().map(|v| (modes.reference[v.atom], v.v)).collect();
        let f = d.field(&atoms).unwrap();
        let norm = d.inner(&f, &f).sqrt();
        let label = &block.fullprof.as_ref().unwrap().0[m.vectors[0].atom].0;
        let orbit = parent.atoms.iter().position(|at| at.label == label.split('_').next().unwrap()).unwrap();
        let r = d.irreps.iter().position(|x| x.name == m.irrep).unwrap();
        let out = d.outside(r, orbit, &f);
        println!("{} on {label}: norm {norm:.6}, outside {out:.2e}", m.name);
        assert!((norm - 1.0).abs() < 2e-5, "{}: norm {norm}", m.name);
        assert!(out < 2e-5, "{}: {out} outside", m.name);
        blocks.entry((m.irrep.clone(), orbit)).or_default().push(m.amplitude);
    }
    for ((irrep, orbit), q) in &blocks {
        let r = d.irreps.iter().find(|x| x.name == *irrep).unwrap();
        let part = r.orbits.iter().find(|x| x.orbit == *orbit).unwrap();
        assert_eq!(part.modes.len(), q.len(), "{irrep} on {orbit}");
        let theirs = q.iter().map(|x| x * x).sum::<f64>().sqrt();
        println!("{irrep} on {}: {:?} against {:?}", parent.atoms[*orbit].label, part.modes.iter().map(|m| m.amplitude).collect::<Vec<_>>(), q);
        assert!((part.amplitude - theirs).abs() < 2e-6, "{irrep}: {} against {theirs}", part.amplitude);
        if q.len() == 1 {
            assert!((part.modes[0].amplitude.abs() - q[0].abs()).abs() < 2e-6);
        }
    }
}

#[test]
fn fullprof_block() {
    let (a, all) = analysis("0_1_P21nToCubicTransData.txt", 2);
    let (parent, distorted, p, o) = from_amplimodes(&a).unwrap();
    let mut d = decompose(&parent, &distorted, &p, o).unwrap();
    d.name_by_modes(&all[2].modes().unwrap());
    let block = d.fullprof_block();
    let atoms: Vec<String> = parent.atoms.iter().map(|x| format!("{}\t1\t-\t{}\t{}\t{}", x.label, x.position[0], x.position[1], x.position[2])).collect();
    let text = format!("Symmetry mode analysis\nHigh symmetry structure\n225\n8.087481 8.087481 8.087481 90 90 90\n{}\n{}\nYou can copy and paste the following text on your .pcr file\n014                    <--Space group symbol\n{block}", atoms.len(), atoms.join("\n"));
    let back = read_amplimodes(&text);
    let m = back[0].modes().unwrap();
    let mut d2 = decompose(&parent, &distorted, &p, o).unwrap();
    assert_eq!(d2.name_by_modes(&m), 7);
    assert_eq!(d2.irreps.iter().map(|r| &r.name).collect::<Vec<_>>(), d.irreps.iter().map(|r| &r.name).collect::<Vec<_>>());
    let ours: Vec<f64> = d.irreps.iter().flat_map(|r| r.orbits.iter().flat_map(|o| o.modes.iter().map(|m| m.amplitude))).collect();
    assert_eq!(m.modes.len(), ours.len());
    for (x, y) in m.modes.iter().zip(&ours) {
        assert!((x.amplitude - y).abs() < 6e-7, "{}: {} against {y}", x.name, x.amplitude);
    }
}

#[test]
fn modes_for_refinement() {
    for (file, k) in [("0_1_P21nToCubicTransData.txt", 2), ("0_2_MonoclinicToCubicTransData.txt", 1), ("0_3_TetragonalToCubicTrans2.txt", 1)] {
        let (a, _) = analysis(file, k);
        let (parent, distorted, p, o) = from_amplimodes(&a).unwrap();
        let d = decompose(&parent, &distorted, &p, o).unwrap();
        let m = d.to_modes("225");
        // r⁰ + ΣQV is the distorted structure.
        let pos = m.positions();
        for (j, (r, u)) in d.reference.iter().zip(&d.displacements).enumerate() {
            for c in 0..3 {
                assert!((pos[j][c] - r.2[c] - u[c]).abs() < 1e-9, "{file}: atom {j}");
            }
        }
        // The refinement's own projection finds the same amplitudes.
        let (q, residual) = m.project(&distorted.group, &d.reference_cell, &d.displacements).unwrap();
        let worst = q.iter().zip(&m.modes).map(|(x, y)| (x - y.amplitude).abs()).fold(0.0, f64::max);
        println!("{file} {k}: {} modes, projection within {worst:.1e}, residual {residual:.1e}", m.modes.len());
        assert!(worst < 1e-9 && residual < 1e-9);
        // Each mode normalised as AMPLIMODES's: Σ over the primitive cell of |V|² = 1 Å².
        let gram = m.gram(&distorted.group, &d.reference_cell).unwrap();
        assert!((0..gram.len()).all(|i| (gram[i][i] - 1.0).abs() < 1e-9 && (0..gram.len()).all(|j| i == j || gram[i][j].abs() < 1e-9)));
    }
}

fn group(s: &str) -> drs::spacegroup::SpaceGroup {
    drs::spacegroup::SpaceGroup::from_symbol(s).unwrap().0
}

fn atom(l: &str, x: [f64; 3]) -> Atom {
    Atom { label: l.into(), species: l.chars().take_while(|c| c.is_ascii_alphabetic()).collect(), position: x }
}

#[test]
fn supercell_r_point() {
    let ap = 3.9;
    let delta = 0.009;
    let a2 = ap * 2f64.sqrt();
    let distorted = Structure {
        group: group("I 4/m c m"),
        cell: drs::cell::Cell::new(a2, a2, 2.0 * ap, 90.0, 90.0, 90.0).unwrap(),
        atoms: vec![atom("Sr1", [0.0, 0.5, 0.25]), atom("Ti1", [0.0, 0.0, 0.0]), atom("O1", [0.0, 0.0, 0.25]), atom("O2", [0.25 - delta, 0.75 - delta, 0.0])],
    };
    let p = drs::mode_decomposition::parse_basis("a-b,a+b,2c").unwrap();
    let cell = drs::cell::Cell::new(ap, ap, ap, 90.0, 90.0, 90.0).unwrap();
    // Sr at the parent's origin (Ti at 1/2,1/2,1/2), then Ti there.
    for (atoms, origin, parity) in [
        (vec![atom("Sr1", [0.0, 0.0, 0.0]), atom("Ti1", [0.5, 0.5, 0.5]), atom("O1", [0.0, 0.5, 0.5])], [0.5, 0.5, 0.5], -1),
        (vec![atom("Ti1", [0.0, 0.0, 0.0]), atom("Sr1", [0.5, 0.5, 0.5]), atom("O1", [0.5, 0.0, 0.0])], [0.0, 0.0, 0.0], 1),
    ] {
        let parent = Structure { group: group("P m -3 m"), cell, atoms };
        let d = decompose(&parent, &distorted, &p, origin).unwrap();
        println!("n {} |G| {} |H| {} N {}: {:?}", d.supercell, d.order, d.order_h, d.primitive_cells, d.irreps.iter().map(|r| (&r.name, &r.star, r.modes, r.amplitude)).collect::<Vec<_>>());
        assert_eq!((d.supercell, d.order, d.primitive_cells), (2, 48 * 8, 4.0));
        assert_eq!(d.irreps.len(), 1);
        let r = &d.irreps[0];
        assert_eq!((r.star.clone(), r.parity, r.dimension, r.modes), (vec![[0.5, 0.5, 0.5]], Some(parity), 3, 1));
        assert_eq!(r.name, format!("k(1/2,1/2,1/2){}[3]", if parity > 0 { "+" } else { "-" }));
        assert!((r.amplitude - 4.0 * delta * ap).abs() < 1e-12 && (d.global - r.amplitude).abs() < 1e-12, "{} against {}", r.amplitude, 4.0 * delta * ap);
    }
}

#[test]
fn refused() {
    let (a, _) = analysis("0_1_P21nToCubicTransData.txt", 2);
    let (parent, distorted, p, o) = from_amplimodes(&a).unwrap();
    let e = decompose(&parent, &distorted, &p, [o[0] + 0.25, o[1], o[2]]).unwrap_err();
    assert!(e.contains("not an operator of the parent"), "{e}");
    let mut fewer = distorted.clone();
    fewer.atoms.pop();
    let e = decompose(&parent, &fewer, &p, o).unwrap_err();
    assert!(e.contains("atoms in the supercell"), "{e}");
}
