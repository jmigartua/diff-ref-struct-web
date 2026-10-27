//! Refinement by symmetry modes (R29, P13; L090), against FullProf's mode refinements and AMPLIMODES's analyses of the
//! course material (`manuals/bcs/mode-crystallography-manual/course-um6p-2024/04_tools_materials/`, README of
//! `manuals/bcs`: the data files and FullProf's runs are the reference, not the book's chapters).
//! - `case_a_positions`: Sr₂ScSbO₆ P2₁/n (case A, `srsclrrt_modes_sb_to_origin`): the positions from FullProf's final
//!   amplitudes (`.new`) against its `.out`'s "Final Complete Structure" (lines 581–588) and each irrep's own
//!   displacements (lines 492–578), the distances in Å, the irrep amplitudes and their esd's, and the positions' esd's,
//!   these two by FullProf's rule (correlations between amplitudes neglected, its note of 29 August 2008).
//! - `case_a_refinement`: the engine refines the 17 parameters of `srsclrrt_modes_sb_to_origin.pcr` (5 amplitudes):
//!   χ² and N − P against the `.out` (2.069 over 2632; FullProf 8.20 reproduces the files, `manuals/bcs/README.md`), each
//!   amplitude within FullProf's esd and its esd within 2 % of FullProf's (scaled as `fullprof_refinements` does).
//!   The pattern itself is compared point by point with the `.prf` in `fullprof_profiles` (index of
//!   `make_fullprof_profiles.py`).
//! - `normalisation`: Σ_j mult_j |c V_mj|² over the distorted primitive cell is 1 Å² and the modes are orthogonal, in
//!   the reference metric (AMPLIMODES's normalisation): case A with AMPLIMODES's parent (0_1, Fm-3m a = 8.087481 Å and
//!   its transformation), case B with each of its two parents, case E (R-3, hexagonal axes) with the parent edge the
//!   La mode implies, the O modes then checked.
//! - `amplimodes_projection`: the distorted structures AMPLIMODES analysed (its atom mappings in 0_1 and 0_4, in the
//!   `.pcr` files' P2₁/n setting) projected onto the `.pcr` files' modes give AMPLIMODES's irrep amplitudes (its "Summary
//!   of Amplitudes", printed to 4 decimals) and, for case A, the per-mode amplitudes of its block for FullProf (0_1).
//! - `case_b_two_parents`: Sr₂CoTeO₆ described from the parent of a = 8.087481 Å (O at x = 0.2541) and of a = 7.9209 Å
//!   (x = 0.2471): the same positions, every amplitude but GM1+'s in the ratio of the two metrics.
//! - `coordinates_projected`: FullProf's coordinate refinement of the same data (`srsclrrt.new`, case A) projected onto
//!   the modes: its amplitudes against the mode refinement's, within their esd's.
//! - `amplimodes_reader`: the analyses of 0_1–0_4 read (structures, transformation, mappings, summary, block for
//!   FullProf).
//! - `control_file_round_trip`: case A through a control file and the `.pcr` exported from it.
//! - `cases_b_d_e_patterns`: B's two parents (their `.new` files) give one pattern.
//! - `parent_reflections`: each reflection of cases A and E in the parent's cell, with its k-vector (L100).
//! - `mode_runs_against_fullprof`: FullProf 8.20's runs of cases A, B, D and E (`manuals/fullprof/runs/modes_runs.py`):
//!   Σw(y − yc)² against their `.out`, the points at the edges of excluded regions where FullProf's choice differs from
//!   the core's named and frozen (L098).

use drs::cell::Cell;
use drs::lsq::{refine, Options};
use drs::modes::{self, Modes};
use drs::pcr;
use drs::powder::Model;
use drs::refinable::ParameterSet;
use drs::spacegroup::SpaceGroup;
use std::path::{Path, PathBuf};

fn tools() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/bcs/mode-crystallography-manual/course-um6p-2024/04_tools_materials")
}

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display())).iter().map(|&b| b as char).collect()
}

const CASE_A: &str = "1_sr2scsbo6_300K_structure_014_data_LR";
const CASE_B: &str = "2_sr2coteo6_300K_structure_014_data_NeutronRT";
const CASE_E: &str = "6_la2comno6_500K_modes_166_data_xrpd";

/// The modes and group of the first phase of a FullProf file.
fn modes_of(dir: &str, file: &str) -> (pcr::PcrPhase, Modes, SpaceGroup) {
    let p = pcr::read(&latin1(&tools().join(dir).join(file))).unwrap_or_else(|e| panic!("{file}: {e}"));
    let ph = p.phases[0].clone();
    let m = Modes::from_pcr(ph.modes.as_ref().expect("V_MODES"), &ph.atoms).unwrap();
    let g = SpaceGroup::from_fullprof_symbol(&ph.symbol).unwrap();
    (ph, m, g)
}

/// A number as FullProf prints it, "−0.55123(19)": value and esd (0 when none is printed).
fn printed(s: &str) -> (f64, f64) {
    match s.split_once('(') {
        Some((v, e)) => {
            let dec = v.split_once('.').map_or(0, |(_, d)| d.len());
            (v.parse().unwrap(), e.trim_end_matches(')').parse::<f64>().unwrap() * 10f64.powi(-(dec as i32)))
        }
        None => (s.parse().unwrap(), 0.0),
    }
}

/// Half a unit of the last decimal printed.
fn half_unit(s: &str) -> f64 {
    let v = s.split('(').next().unwrap();
    0.5 * 10f64.powi(-(v.split_once('.').map_or(0, |(_, d)| d.len()) as i32))
}

/// The amplitudes' esd's of FullProf's last cycle (`.out`, "Amplitudes of symmetry modes ==>": name, value, shift,
/// sigma).
fn fullprof_amplitudes(out: &str) -> Vec<(String, f64, f64)> {
    let ls: Vec<&str> = out.lines().collect();
    let i = ls.iter().rposition(|l| l.contains("=> Amplitudes of symmetry modes ==>")).expect("amplitudes block");
    ls[i + 2..].iter().map_while(|l| {
        let f: Vec<&str> = l.split_whitespace().collect();
        (f.len() == 4).then(|| (f[0].to_string(), f[1].parse().unwrap(), f[3].parse().unwrap()))
    }).collect()
}

/// The reference metric of an AMPLIMODES analysis (its parent cell transformed to the distorted phase's basis).
fn reference_metric(file: &str, k: usize) -> Cell {
    modes::read_amplimodes(&latin1(&tools().join(file)))[k].reference_cell().expect("parent and transformation")
}

#[test]
fn amplimodes_reader() {
    let mut fails = Vec::new();
    // (file, analyses, of which with mappings, with the block for FullProf)
    for (file, n, maps, fp) in [("0_1_P21nToCubicTransData.txt", 3, 2, 1), ("0_2_MonoclinicToCubicTransData.txt", 1, 1, 0), ("0_3_TetragonalToCubicTrans2.txt", 1, 1, 0), ("0_4_co_te_data_012.txt", 6, 6, 0)] {
        let a = modes::read_amplimodes(&latin1(&tools().join(file)));
        let got = (a.len(), a.iter().filter(|x| !x.mappings.is_empty()).count(), a.iter().filter(|x| x.fullprof.is_some()).count());
        println!("{file}: {} analyses; {} with atom mappings; {} with a block for FullProf; irreps {:?}", got.0, got.1, got.2, a.iter().map(|x| x.summary.len()).collect::<Vec<_>>());
        if got != (n, maps, fp) {
            fails.push(format!("{file}: {got:?} for {:?}", (n, maps, fp)));
        }
        for (k, x) in a.iter().enumerate() {
            // Each summary's irrep amplitudes square-sum to the global distortion (printed to 4 decimals).
            if let Some(g) = x.global {
                let s = x.summary.iter().map(|r| r.amplitude.powi(2)).sum::<f64>().sqrt();
                if (s - g).abs() > 2e-4 {
                    fails.push(format!("{file} analysis {}: √Σ A² = {s:.4}, global {g}", k + 1));
                }
            }
            if x.high.is_none() || x.reference.is_none() || x.transformation.is_none() || x.summary.is_empty() {
                fails.push(format!("{file} analysis {}: parent, reference, transformation or summary missing", k + 1));
            }
            // The reference cell AMPLIMODES prints is its parent transformed (P2₁/c or I2/m settings: 6 decimals).
            if let (Some(r), Some(c)) = (&x.reference, x.reference_cell()) {
                let ours = [c.a, c.b, c.c, c.alpha, c.beta, c.gamma];
                // AMPLIMODES's P2₁/n analyses print the standard setting's reference cell; compare when the settings agree.
                if x.low_setting.is_none() && ours.iter().zip(&r.cell).any(|(a, b)| (a - b).abs() > 5e-5 * b.abs().max(1.0)) {
                    fails.push(format!("{file} analysis {}: reference cell {ours:?}, printed {:?}", k + 1, r.cell));
                }
            }
        }
    }
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}

#[test]
fn normalisation() {
    let mut fails = Vec::new();
    let check = |name: &str, m: &Modes, g: &SpaceGroup, metric: &Cell| -> Option<String> {
        let gram = m.gram(g, metric).unwrap();
        let n = m.modes.len();
        let (mut worst_d, mut worst_o) = (0.0f64, 0.0f64);
        for (i, row) in gram.iter().enumerate() {
            for (j, g) in row.iter().enumerate() {
                if i == j {
                    worst_d = worst_d.max((g - 1.0).abs());
                } else {
                    worst_o = worst_o.max(g.abs());
                }
            }
        }
        // The vectors are printed to 6 decimals (|V| ≈ 0.03–0.09): a relative 10⁻⁵ on each, 2·10⁻⁵ on the sums.
        println!("{name}: {n} modes, |Σ|cV|² − 1| ≤ {worst_d:.1e} Å², |off-diagonal| ≤ {worst_o:.1e} Å², metric {:.6} {:.6} {:.6} {:.4} {:.4} {:.4}", metric.a, metric.b, metric.c, metric.alpha, metric.beta, metric.gamma);
        (worst_d > 5e-5 || worst_o > 5e-5).then(|| format!("{name}: diagonal off by {worst_d:.1e}, off-diagonal {worst_o:.1e}"))
    };
    // Case A: AMPLIMODES's analysis in the P2₁/n setting (0_1, the first), parent a = 8.087481 Å.
    let (_, m, g) = modes_of(CASE_A, "srsclrrt_modes_sb_to_origin.pcr");
    fails.extend(check("case A (Sr2ScSbO6, P21/n)", &m, &g, &reference_metric("0_1_P21nToCubicTransData.txt", 0)));
    // Case B: the two parents of 0_4 (analyses 1 and 2), the modes of each .new.
    for (file, k) in [("CoP21n_modes_parent_1.new", 0), ("CoP21n_modes_parent_2.new", 1)] {
        let (_, m, g) = modes_of(CASE_B, file);
        fails.extend(check(&format!("case B ({file})"), &m, &g, &reference_metric("0_4_co_te_data_012.txt", k)));
    }
    // Case E: R-3 on hexagonal axes from Fm-3m: a_h = (a − b)/2, b_h = (b − c)/2, c_h = a + b + c. No AMPLIMODES output
    // is in the material: the parent edge is the one that normalises the La mode (one atom, 2 in the primitive cell),
    // and the three O modes are then checked with it.
    let (_, m, g) = modes_of(CASE_E, "R-3.new");
    let p = [[0.5, 0.0, 1.0], [-0.5, 0.5, 1.0], [0.0, -0.5, 1.0]];
    let unit = modes::transform_cell(&Cell::new(1.0, 1.0, 1.0, 90.0, 90.0, 90.0).unwrap(), &p).unwrap();
    let la = m.modes.iter().position(|x| x.vectors.len() == 1 && m.reference[x.vectors[0].atom][2] == 0.25).expect("the La mode");
    let a_c = (1.0 / m.gram(&g, &unit).unwrap()[la][la]).sqrt();
    println!("case E: the La mode is normalised by a parent of a = {a_c:.5} Å (Sr2CoTeO6's parent in 0_4 is 7.9209 Å)");
    if (a_c - 7.9209).abs() > 2e-4 {
        fails.push(format!("case E: parent edge {a_c:.5}"));
    }
    let metric = modes::transform_cell(&Cell::new(a_c, a_c, a_c, 90.0, 90.0, 90.0).unwrap(), &p).unwrap();
    fails.extend(check("case E (La2CoMnO6, R-3 hexagonal)", &m, &g, &metric));
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}

#[test]
fn amplimodes_projection() {
    let mut fails = Vec::new();
    // (AMPLIMODES file, analysis, FullProf file whose modes and reference it shares)
    for (file, k, dir, pcr_file) in [("0_1_P21nToCubicTransData.txt", 0, CASE_A, "srsclrrt_modes_sb_to_origin.pcr"), ("0_4_co_te_data_012.txt", 0, CASE_B, "CoP21n_modes_parent_1.new"), ("0_4_co_te_data_012.txt", 1, CASE_B, "CoP21n_modes_parent_2.new")] {
        let an = &modes::read_amplimodes(&latin1(&tools().join(file)))[k];
        let (ph, mut m, g) = modes_of(dir, pcr_file);
        let metric = an.reference_cell().unwrap();
        // The mappings pair each reference atom (by label) with its distorted partner, both in P2₁/n.
        let mut u = Vec::new();
        for at in &ph.atoms {
            let Some(mp) = an.mappings.iter().find(|x| x.0 == at.label) else { panic!("{file}: no mapping for {}", at.label) };
            if mp.1.iter().zip(&at.x).any(|(a, b)| (a - b).abs() > 1e-5) {
                fails.push(format!("{file} {k}: {}'s reference {:?} is not the .pcr's {:?}", at.label, mp.1, at.x));
            }
            u.push([0, 1, 2].map(|i| {
                let d = mp.3[i] - mp.1[i];
                d - d.round()
            }));
        }
        let (q, residual) = m.project(&g, &metric, &u).unwrap();
        for (mode, a) in m.modes.iter_mut().zip(&q) {
            mode.amplitude = *a;
        }
        let ours = m.irreps();
        let mut line = Vec::new();
        for s in &an.summary {
            let a = ours.iter().find(|r| r.name == s.irrep).map_or(0.0, |r| r.amplitude);
            line.push(format!("{} {a:.4} ({})", s.irrep, s.amplitude));
            if (a - s.amplitude).abs() > 5e-5 + 1e-5 * a {
                fails.push(format!("{file} {k}: {} {a:.5} against AMPLIMODES's {}", s.irrep, s.amplitude));
            }
        }
        println!("{file} analysis {}: {}; the part no mode describes {residual:.1e} Å", k + 1, line.join(", "));
        if residual > 1e-4 {
            fails.push(format!("{file} {k}: residual {residual}"));
        }
        // Case A: the per-mode amplitudes of AMPLIMODES's block for FullProf (0_1, the analysis in P2₁/c): the same
        // modes for the non-degenerate irreps (GM1+, GM3+, GM4+, X2+, X3+: one mode each); within GM5+ and X5+ the basis
        // differs with the setting, and only the irrep's amplitude is compared.
        if k == 0 && dir == CASE_A {
            let fp = modes::read_amplimodes(&latin1(&tools().join(file)))[2].fullprof.clone().expect("block for FullProf");
            for (name, a) in &fp.2 {
                let irrep = modes::irrep_of_name(name).unwrap();
                // The sign of a one-mode irrep's basis vector is a convention of each setting (X2+ comes out with the
                // opposite sign in P2₁/c): the magnitudes are compared.
                if an.summary.iter().find(|s| s.irrep == irrep).is_some_and(|s| s.dimension == 1) {
                    let ours = m.modes.iter().find(|x| x.name == *name).unwrap().amplitude;
                    if (ours.abs() - a.abs()).abs() > 5e-6 + 1e-5 * a.abs() {
                        fails.push(format!("{name}: {ours:.6} against AMPLIMODES's {a}"));
                    }
                }
            }
        }
    }
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}

#[test]
fn case_b_two_parents() {
    let (ph1, m1, _) = modes_of(CASE_B, "CoP21n_modes_parent_1.new");
    let (_, m2, _) = modes_of(CASE_B, "CoP21n_modes_parent_2.new");
    let mut fails = Vec::new();
    // The same structure: positions to the files' 6 decimals of the amplitudes (a 5·10⁻⁷ Å error on Q moves x by
    // ≤ 5·10⁻⁸).
    let mut dpos = 0.0f64;
    for (j, (a, b)) in m1.positions().iter().zip(m2.positions()).enumerate() {
        let d = (0..3).map(|k| (a[k] - b[k]).abs()).fold(0.0, f64::max);
        dpos = dpos.max(d);
        if d > 2e-6 {
            fails.push(format!("{}: {a:?} and {b:?}", ph1.atoms[j].label));
        }
    }
    println!("case B: positions from the two parents differ by at most {dpos:.1e}");
    // Every amplitude but GM1+'s (which also absorbs the shift of the reference O, 0.2541 → 0.2471) in the ratio of
    // the parents' edges: the vectors scale by a₁/a₂, the amplitudes by a₂/a₁.
    let r = 7.9209 / 8.087481;
    let mut worst = 0.0f64;
    for (x, y) in m1.modes.iter().zip(&m2.modes) {
        if x.irrep == "GM1+" {
            println!("case B: GM1+ {} (parent 1) and {} (parent 2): the reference O moved by 0.007 in x", x.amplitude, y.amplitude);
            continue;
        }
        // The bound: the two amplitudes' rounding to 6 decimals, and the two files' positions, which agree to dpos:
        // a difference dpos along the mode's largest vector component is dpos/|V|max in its amplitude.
        let vmax = x.vectors.iter().flat_map(|v| v.v).fold(0.0f64, |a, b| a.max(b.abs()));
        let tol = 5e-7 * (1.0 + r) + dpos / vmax;
        let d = (y.amplitude - r * x.amplitude).abs();
        worst = worst.max(d / tol);
        if d > tol {
            fails.push(format!("{}: {} and {} for the ratio {r:.6}", x.name, x.amplitude, y.amplitude));
        }
    }
    println!("case B: positions alike from both parents; the 11 other amplitudes in the ratio 7.9209/8.087481 = {r:.6}, within {:.0} % of the rounding's bound", 100.0 * worst);
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}

#[test]
fn case_a_positions() {
    let dir = tools().join(CASE_A);
    let p = pcr::read(&latin1(&dir.join("srsclrrt_modes_sb_to_origin.new"))).unwrap();
    let ph = &p.phases[0];
    let m = Modes::from_pcr(ph.modes.as_ref().unwrap(), &ph.atoms).unwrap();
    let out = latin1(&dir.join("srsclrrt_modes_sb_to_origin.out"));
    let ls: Vec<&str> = out.lines().collect();
    let c = ph.cell;
    let cell = Cell::new(c[0], c[1], c[2], c[3], c[4], c[5]).unwrap();
    let g = cell.metric();
    let dist = |d: [f64; 3]| drs::mat3::quad(&g, &d).sqrt();
    let sig: Vec<f64> = fullprof_amplitudes(&out).iter().map(|x| x.2).collect();
    let diag = |a: usize, b: usize| if a == b { sig[a] * sig[a] } else { 0.0 };
    let mut fails = Vec::new();
    let (mut checked, mut esds) = (0, 0);
    // Each irrep's block: its displacements alone and their length.
    for r in m.irreps() {
        let head = format!("single Irreducible representation:{}", r.name);
        let i = ls.iter().position(|l| l.contains(&head)).unwrap_or_else(|| panic!("{head}"));
        let mut one = m.clone();
        for (k, mode) in one.modes.iter_mut().enumerate() {
            if !r.modes.contains(&k) {
                mode.amplitude = 0.0;
            }
        }
        for (j, at) in ph.atoms.iter().enumerate() {
            let l = ls[i + 2 + j];
            let f: Vec<&str> = l.split_whitespace().collect();
            assert_eq!(f[1], at.label, "{l}");
            let u = one.displacement(j);
            for k in 0..3 {
                if (u[k] - f[6 + k].parse::<f64>().unwrap()).abs() > half_unit(f[6 + k]) + 1e-9 {
                    fails.push(format!("{}: {} d{} {:.6} against {}", r.name, at.label, ["x", "y", "z"][k], u[k], f[6 + k]));
                }
                checked += 1;
            }
            if (dist(u) - f[9].parse::<f64>().unwrap()).abs() > half_unit(f[9]) + 1e-9 {
                fails.push(format!("{}: {} moves {:.6} Å, FullProf {}", r.name, at.label, dist(u), f[9]));
            }
        }
        // "Global Amplitude of Representation GM1+ : 0.283(5)".
        let l = ls[i + 2 + ph.atoms.len()..].iter().find(|l| l.contains("Global Amplitude of Representation")).unwrap();
        let s = l.rsplit(':').next().unwrap().trim();
        let (v, e) = printed(s);
        let e_ours = m.irrep_esd(&r, diag);
        if (r.amplitude - v).abs() > half_unit(s) + 1e-9 || (e > 0.0 && (e_ours - e).abs() > half_unit(s) + 1e-9) {
            fails.push(format!("{}: amplitude {:.5}({:.5}) against FullProf's {s}", r.name, r.amplitude, e_ours));
        }
        println!("{:5} A = {:.6} ({:.6}; FullProf {s}); polarisation {:?}", r.name, r.amplitude, e_ours, r.polarisation.iter().map(|c| format!("{c:.4}")).collect::<Vec<_>>());
    }
    // The final structure.
    let i = ls.iter().position(|l| l.contains("Final Complete Structure with contribution of all modes")).unwrap();
    for (j, (at, x)) in ph.atoms.iter().zip(m.positions()).enumerate() {
        let f: Vec<&str> = ls[i + 2 + j].split_whitespace().collect();
        assert_eq!(f[1], at.label);
        let e = m.position_esd(j, diag);
        for k in 0..3 {
            let (v, ef) = printed(f[3 + k]);
            if (x[k] - v).abs() > half_unit(f[3 + k]) + 1e-9 {
                fails.push(format!("{} {}: {:.6} against {}", at.label, ["x", "y", "z"][k], x[k], f[3 + k]));
            }
            if ef > 0.0 {
                esds += 1;
                if (e[k] - ef).abs() > half_unit(f[3 + k]) + 1e-12 {
                    fails.push(format!("{} σ{}: {:.6} against {}", at.label, ["x", "y", "z"][k], e[k], f[3 + k]));
                }
            }
        }
        let u = m.displacement(j);
        // Biso, then the occupancy left blank, dx dy dz and the distance (the last column).
        let dl = f[f.len() - 1];
        if (dist(u) - dl.parse::<f64>().unwrap()).abs() > half_unit(dl) + 1e-9 {
            fails.push(format!("{} moves {:.6} Å, FullProf {dl}", at.label, dist(u)));
        }
        println!("{:6} {:9.5} {:9.5} {:9.5}  σ {:.5} {:.5} {:.5}  |u| {:.5} Å   FullProf {} {dl}", at.label, x[0], x[1], x[2], e[0], e[1], e[2], dist(u), f[3..6].join(" "));
    }
    println!("{checked} displacements of the seven irreps' blocks, {} positions and {esds} of their esd's against FullProf's .out", 3 * ph.atoms.len());
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}

#[test]
fn case_a_refinement() {
    let dir = tools().join(CASE_A);
    let stem = "srsclrrt_modes_sb_to_origin";
    let p = pcr::read(&latin1(&dir.join(format!("{stem}.pcr")))).unwrap();
    let l = Model::from_pcr(&p, &dir, stem).unwrap();
    let data = l.data.as_ref().unwrap();
    let mut m = l.model.clone();
    let out = latin1(&dir.join(format!("{stem}.out")));
    // FullProf's reflections past the data, as far as S(1/d)max of its .out (L035; `fullprof_refinements` does the
    // same): without them the tails of the peaks beyond 160° are missing and χ² falls by 0.13 %.
    let lim: f64 = out.lines().rev().find_map(|l| l.split_once("S(1/d)max:")).map(|(_, r)| r.split_whitespace().next().unwrap().parse().unwrap()).unwrap();
    m.phases[0].reflection_limit = Some(lim);
    let set = ParameterSet::from_pcr(&p, &m).unwrap();
    let r = refine(&mut m, data, &set, Options { max_cycles: 100, eps: p.eps }, |_| {});
    let chi2_fp: f64 = out.lines().rfind(|l| l.contains("Conventional Rietveld Rp,Rwp,Re and Chi2:")).unwrap().split(':').nth(1).unwrap().split_whitespace().nth(3).unwrap().parse().unwrap();
    let dof_fp: usize = out.lines().rfind(|l| l.contains("=> N-P+C:")).unwrap().split(':').nth(1).unwrap().trim().parse().unwrap();
    // χ² at FullProf's values (its .pcr: the run was converged when written), against the .out's.
    let mut at_fp = l.model.clone();
    at_fp.phases[0].reflection_limit = Some(lim);
    let s_fp = drs::lsq::statistics(&at_fp, data, &at_fp.calculate(&data.x).ycalc, 17);
    let mut fails = Vec::new();
    println!("{} parameters, {} cycles ({}): χ² {:.4} over N − P = {}; FullProf {chi2_fp} over {dof_fp}", set.free.len(), r.cycles, r.stop, r.stats.chi2, r.stats.n - r.stats.p);
    if set.free.len() != 17 || r.stats.n - r.stats.p != dof_fp {
        fails.push(format!("{} parameters, N − P {}", set.free.len(), r.stats.n - r.stats.p));
    }
    println!("χ² at FullProf's values {:.4}", s_fp.chi2);
    if (s_fp.chi2 - chi2_fp).abs() > 5e-4 + 1.5e-3 * chi2_fp {
        fails.push(format!("χ² at FullProf's values {:.4} against {chi2_fp}", s_fp.chi2));
    }
    if (r.stats.chi2 - chi2_fp).abs() > 5e-4 + 1.5e-3 * chi2_fp {
        fails.push(format!("χ² {:.4} against FullProf's {chi2_fp}", r.stats.chi2));
    }
    // The amplitudes and their esd's (FullProf scales by χ² over the points with Bragg contributions, the engine over
    // all: here the two are equal, 2.07, so the scaling is 1).
    let fp = fullprof_amplitudes(&out);
    for (free, ours) in set.free.iter().zip(&r.params) {
        let drs::jacobian::Field::Mode(_, k) = free.fields[0].0 else { continue };
        let (name, v, e) = &fp[k];
        let dev = (ours.value - v).abs() / e;
        let q = ours.esd / e;
        println!("{name:8} {:10.6} ({:.6})  FullProf {v:10.6} ({e:.6}): {dev:.3} esd; esd ratio {q:.4}", ours.value, ours.esd);
        if dev > 1.0 || !(0.98..=1.02).contains(&q) {
            fails.push(format!("{name}: {:.6} ({:.6}) against FullProf's {v} ({e})", ours.value, ours.esd));
        }
    }
    // The irrep amplitudes reached, with the engine's covariance (correlations included).
    let ms = m.phases[0].modes.as_ref().unwrap();
    let cols: Vec<Option<usize>> = (0..ms.modes.len()).map(|k| set.free.iter().position(|f| f.fields[0].0 == drs::jacobian::Field::Mode(0, k))).collect();
    let cov = |a: usize, b: usize| match (cols[a], cols[b]) {
        (Some(i), Some(j)) => r.covariance[i][j],
        _ => 0.0,
    };
    for ir in ms.irreps().iter().filter(|x| x.amplitude > 0.0) {
        println!("  {:5} A = {:.4}({:.0})", ir.name, ir.amplitude, 1e4 * ms.irrep_esd(ir, cov));
    }
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}

#[test]
fn coordinates_projected() {
    // FullProf's coordinate refinement of the same data (srsclrrt.new: P 1 21/N 1, Sb at (0, 1/2, 0)) moved to the
    // modes' origin (Sb at 0: x, y + 1/2, z, AMPLIMODES's "Sb to the origin" in 0_1) and projected onto the modes.
    let dir = tools().join(CASE_A);
    let c = pcr::read(&latin1(&dir.join("srsclrrt.new"))).unwrap();
    let (ph, mut m, g) = modes_of(CASE_A, "srsclrrt_modes_sb_to_origin.new");
    let metric = reference_metric("0_1_P21nToCubicTransData.txt", 0);
    let reference: Vec<(String, [f64; 3])> = ph.atoms.iter().map(|a| (a.species.clone(), a.x)).collect();
    let distorted: Vec<(String, [f64; 3])> = c.phases[0].atoms.iter().map(|a| (a.species.clone(), a.x)).collect();
    let u = modes::displacements(&reference, &distorted, &g, &metric, [0.0, 0.5, 0.0], 0.5).unwrap();
    let (q, residual) = m.project(&g, &metric, &u).unwrap();
    let fp = fullprof_amplitudes(&latin1(&dir.join("srsclrrt_modes_sb_to_origin.out")));
    let mut fails = Vec::new();
    for (k, mode) in m.modes.iter_mut().enumerate() {
        let (name, v, e) = &fp[k];
        // Modes fixed at 0 in the mode refinement (e = 0) have only the coordinate refinement's value.
        let note = if *e > 0.0 { format!("{:.2} esd", (q[k] - v).abs() / e) } else { "fixed at 0 there".to_string() };
        println!("{name:8} coordinates {:9.5}   modes {v:9.5} ({e:.5})  {note}", q[k]);
        if *e > 0.0 && (q[k] - v).abs() > 3.0 * e {
            fails.push(format!("{name}: {:.5} from the coordinates, {v} ({e}) from the modes", q[k]));
        }
        mode.amplitude = q[k];
    }
    println!("the part of the coordinate refinement no mode describes: {residual:.1e} Å; irreps {:?}", m.irreps().iter().map(|r| format!("{} {:.4}", r.name, r.amplitude)).collect::<Vec<_>>());
    if residual > 1e-3 {
        fails.push(format!("residual {residual}"));
    }
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}

#[test]
fn control_file_round_trip() {
    // Case A's .pcr as a control file and back: the modes (reference, vectors, amplitudes, Max_Amplitude, the FST
    // integers), the positions they give and the parameter set (its addresses name the amplitudes) survive the control
    // file, its text written again is identical, and the .pcr exported from it reads back to the same modes and the
    // same 17 parameters.
    let dir = tools().join(CASE_A);
    let stem = "srsclrrt_modes_sb_to_origin";
    let p = pcr::read(&latin1(&dir.join(format!("{stem}.pcr")))).unwrap();
    let c = drs::control::Control::from_pcr(&p, &dir, stem, &dir).unwrap();
    let text = c.to_text();
    assert!(text.contains("[[phase.mode]]") && text.contains("phase1.A3_GM4+.mode_amplitude"), "the modes and their parameters are written");
    let back = drs::control::Control::from_text_in(&text, &dir).unwrap_or_else(|d| panic!("{d:?}"));
    let (m0, m1) = (c.model.phases[0].modes.as_ref().unwrap(), back.model.phases[0].modes.as_ref().unwrap());
    assert_eq!(m0, m1);
    let atoms = |m: &Model| match &m.phases[0].intensities {
        drs::powder::Intensities::Structure(st, _) => st.atoms.iter().map(|a| a.x).collect::<Vec<_>>(),
        _ => unreachable!(),
    };
    assert_eq!(atoms(&c.model), atoms(&back.model));
    assert_eq!(back.to_text(), text, "written again identically");
    let names = |c: &drs::control::Control| c.refinement.build(&c.model).unwrap().free.iter().map(|f| f.name.clone()).collect::<Vec<_>>();
    assert_eq!(names(&c), names(&back));
    let e = drs::control::pcr_export::export(&back, stem, None).unwrap();
    let q = pcr::read(&e.pcr).unwrap();
    let l = Model::from_pcr_with_data(&q, &dir, stem, Some(&dir.join("SrScLRRT.dat"))).unwrap();
    assert_eq!(l.model.phases[0].modes.as_ref().unwrap(), m0, "the exported .pcr's modes");
    let set = ParameterSet::from_pcr(&q, &l.model).unwrap();
    assert_eq!(set.free.len(), 17);
    println!("case A: {} modes through the control file and the exported .pcr; {} parameters; notes {:?}", m0.modes.len(), set.free.len(), e.notes);
    // The parent (AMPLIMODES's analysis 0_1) through the control file (mode_parent_group, mode_parent_basis,
    // mode_parent_origin; L100).
    let mut c2 = c.clone();
    c2.model.phases[0].modes.as_mut().unwrap().parent = modes::read_amplimodes(&latin1(&tools().join("0_1_P21nToCubicTransData.txt")))[0].mode_parent();
    let t2 = c2.to_text();
    assert!(t2.contains("mode_parent_basis"), "the parent is written");
    let b2 = drs::control::Control::from_text_in(&t2, &dir).unwrap_or_else(|d| panic!("{d:?}"));
    assert_eq!(b2.model.phases[0].modes.as_ref().unwrap().parent, c2.model.phases[0].modes.as_ref().unwrap().parent);
    assert_eq!(b2.to_text(), t2);
}

#[test]
fn cases_b_d_e_patterns() {
    // B's two parents describe one structure: their .new files (the same values from either parent) give one pattern.
    let d = tools().join(CASE_B);
    let chi2: Vec<f64> = ["CoP21n_modes_parent_1", "CoP21n_modes_parent_2"]
        .iter()
        .map(|stem| {
            let p = pcr::read(&latin1(&d.join(format!("{stem}.new")))).unwrap();
            let l = Model::from_pcr_fullprof_with_data(&p, &d, stem, Some(&d.join("Sr2CoTeO6_RT.dat"))).unwrap_or_else(|e| panic!("{stem}: {e}"));
            let obs = l.data.as_ref().unwrap();
            drs::lsq::statistics(&l.model, obs, &l.model.calculate(&obs.x).ycalc, 0).chi2
        })
        .collect();
    assert!((chi2[0] - chi2[1]).abs() < 1e-6 * chi2[0], "case B: the two parents give χ² {} and {}", chi2[0], chi2[1]);
}

/// Points at the edge of an excluded region that FullProf calculates and counts and the core leaves out (`true`), or
/// the reverse (`false`), by run (L098). FullProf turns a region into points with its own step, the "Step" of its
/// `.out`: for these files of x, y, σ (Ins = 10) it is not the data's (B 0.050018 for 0.05, D 0.100078 for 0.1, E
/// 0.066710 for 0.06667), and for the same data it changes between runs (Yb₂ZrWO₈: 0.040013 and 0.040016), so the
/// core cannot derive it and keeps its rule (L033). In B, 17.000 (region 0–17) has its peaks and counts in FullProf,
/// 149.95 (region 150–180) neither; in E, 21.3335 (region 0–21.3) has them. E lists its regions out of order (120–160,
/// then 0–21.3): FullProf calculates peaks within the second too (its `.prf`), though it leaves its points out of Σ.
/// Frozen: a change that adds a point here fails.
const EDGE_POINTS: &[(&str, &[(f64, bool)])] = &[
    ("B_article", &[(17.0, true), (149.95, false)]),
    ("B_parent_1", &[(17.0, true), (149.95, false)]),
    ("B_parent_2", &[(17.0, true), (149.95, false)]),
    ("B_thesis", &[(17.0, true), (149.95, false)]),
    ("E", &[(21.33346, true)]),
];

/// FullProf 8.20's runs of the course's cases (manuals/fullprof/runs/modes_runs.py, L093), each CODE.pcr (FullProf's
/// refined values) on its data with the reflection limits of CODE_calc.out (L035), against Σw(y − yc)² of
/// CODE_calc.out ("Residual", printed to 4 digits). The 1–6 % of L093 came from the reflection limits (the `.new`
/// files were taken without them), from D's U, V, W of λ/2 (Ratio < 0, line 32, ignored before 0.23.0) and from D's
/// linear background below its first point (`Convention::background_edge`); L098.
/// - The point set and the weights: Σ with FullProf's own Ycal (its `_calc.prf`) over the points the core counts,
///   EDGE_POINTS changed to FullProf's choice, within the print of the `.out` (half a unit of its 4th digit) and of
///   the `.prf`.
/// - The core's Σ within 0.15 % of FullProf's, as `fullprof_profiles::check_chi2` holds χ² (the window edges of A, L035,
///   and EDGE_POINTS); the profile itself is compared point by point in `fullprof_profiles`.
#[test]
fn mode_runs_against_fullprof() {
    let runs = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs/modes");
    let cases = [
        ("A_modes", "SrScLRRT.dat"),
        ("A_coords", "SrScLRRT.dat"),
        ("B_article", "Sr2CoTeO6_RT.dat"),
        ("B_parent_1", "Sr2CoTeO6_RT.dat"),
        ("B_parent_2", "Sr2CoTeO6_RT.dat"),
        ("B_thesis", "Sr2CoTeO6_RT.dat"),
        ("D", "L2CoMnO6_300K_130589.dat"),
        ("E", "t7.dat"),
    ];
    let mut fails = Vec::new();
    for (run, data) in cases {
        let dir = runs.join(run);
        let code = run.to_lowercase();
        let p = pcr::read(&latin1(&dir.join(format!("{code}.pcr")))).unwrap();
        let mut l = Model::from_pcr_fullprof_with_data(&p, &dir, &code, Some(&dir.join(data))).unwrap_or_else(|e| panic!("{run}: {e}"));
        let out = latin1(&dir.join(format!("{code}_calc.out")));
        let limits: Vec<f64> = out.lines().filter_map(|l| l.split_once("S(1/d)max:")).map(|(_, r)| r.split_whitespace().next().unwrap().parse().unwrap()).collect();
        for (ph, lim) in l.model.phases.iter_mut().zip(&limits) {
            ph.reflection_limit = Some(*lim);
        }
        l.model.fullprof_step = out.lines().find_map(|l| l.split_once("Step:")).and_then(|(_, r)| r.split_whitespace().next()?.parse().ok());
        let (m, obs) = (&l.model, l.data.as_ref().unwrap());
        let resid_text = out.lines().skip_while(|l| !l.contains("SumYdif")).nth(1).and_then(|l| l.split_whitespace().nth(4)).unwrap();
        let resid: f64 = resid_text.parse().unwrap();
        // Four significant digits ("0.8641E+05").
        let half_unit = 0.5 * 10f64.powi(resid.log10().floor() as i32 - 3);
        let _ = resid_text;
        // FullProf's Ycal at the data's points, in the data's units.
        let prf = latin1(&dir.join(format!("{code}_calc.prf")));
        let lines: Vec<&str> = prf.lines().collect();
        let start = lines.iter().position(|l| l.starts_with(" 2Theta") && l.contains("Yobs")).unwrap() + 1;
        let rows: Vec<[f64; 3]> = lines[start..]
            .iter()
            .map_while(|l| {
                let f: Vec<f64> = l.split('\t').take(3).filter_map(|t| t.trim().parse().ok()).collect();
                (f.len() == 3).then(|| [f[0], f[1], f[2]])
            })
            .collect();
        assert_eq!(rows.len(), obs.x.len(), "{run}: the .prf has the data's points");
        let k = rows.iter().zip(&obs.y).find(|(r, &y)| r[1] > 1.0 && y > 1.0).map(|(r, &y)| (y / r[1]).log10().round()).unwrap();
        let yc_fp: Vec<f64> = rows.iter().map(|r| r[2] * 10f64.powi(k as i32)).collect();
        // The weights of every point (none left out), and the points the core counts.
        let mut open = m.clone();
        open.excluded.clear();
        let w_fp = drs::lsq::weights(&open, obs, &yc_fp);
        let mut counted: Vec<bool> = m.excluded_points(&obs.x).iter().map(|e| !e).collect();
        let edges = EDGE_POINTS.iter().find(|e| e.0 == run).map_or(&[][..], |e| e.1);
        for &(t, fp_counts) in edges {
            let i = obs.x.iter().position(|&x| (x - t).abs() < 1e-3).unwrap_or_else(|| panic!("{run}: no point at {t}"));
            if counted[i] == fp_counts {
                fails.push(format!("{run}: the core now treats {t} as FullProf does: take it out of EDGE_POINTS"));
            }
            counted[i] = fp_counts;
        }
        let s_fp: f64 = (0..obs.x.len()).filter(|&i| counted[i]).map(|i| w_fp[i] * (obs.y[i] - yc_fp[i]).powi(2)).sum();
        let yc = m.calculate(&obs.x).ycalc;
        let w = drs::lsq::weights(m, obs, &yc);
        let s: f64 = (0..obs.x.len()).map(|i| w[i] * (obs.y[i] - yc[i]).powi(2)).sum();
        println!("{run:11} Σw(y − yc)²: FullProf {resid} (± {half_unit}); with its Ycal over the core's points {s_fp:.1}; the core's {s:.1} ({:+.2e})", s / resid - 1.0);
        if (s_fp - resid).abs() > half_unit + 1e-5 * resid {
            fails.push(format!("{run}: Σ with FullProf's Ycal {s_fp:.1} against its {resid}"));
        }
        if (s - resid).abs() > 1.5e-3 * resid {
            fails.push(format!("{run}: the core's Σ {s:.1} against FullProf's {resid}"));
        }
    }
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}

/// Each reflection of a mode phase in its parent's terms (W8, app request 20; L100): the parent indices h·P⁻¹ give h
/// back through P, a reflection is of k = 0 (Γ) exactly when the parent's lattice has it (integer indices its centring
/// allows), and the k-vectors are those of the irreps the phase's modes carry. Case A: P2₁/n from Fm-3m by AMPLIMODES's
/// analysis 0_1 (modes of GM and X: k = 0 and the star of X, (1 0 0), (0 1 0), (0 0 1)). Case E: R-3 on hexagonal axes
/// from Fm-3m, obverse (a_h = (b − a)/2, b_h = (c − b)/2, c_h = a + b + c, the metric of `normalisation`; its signs
/// there give the reverse setting, whose R centring is not the parent's F): modes of Γ only, so every reflection is the
/// parent's.
#[test]
fn parent_reflections() {
    let mut fails = Vec::new();
    let a = modes::read_amplimodes(&latin1(&tools().join("0_1_P21nToCubicTransData.txt")))[0].mode_parent().expect("analysis 0_1: parent and transformation");
    let p_e = modes::ModeParent { group: "F m -3 m".into(), basis: [[-0.5, 0.0, 1.0], [0.5, -0.5, 1.0], [0.0, 0.5, 1.0]], origin: [0.0; 3] };
    let (fm3m, _) = SpaceGroup::from_symbol("F m -3 m").unwrap();
    for (name, dir, file, parent, labels) in [("case A", CASE_A, "srsclrrt_modes_sb_to_origin.pcr", a, &["GM", "X"][..]), ("case E", CASE_E, "R-3.new", p_e, &["GM"][..])] {
        let (ph, m, g) = modes_of(dir, file);
        let irreps: std::collections::BTreeSet<String> = m.modes.iter().map(|x| x.irrep.trim_end_matches(['+', '-']).trim_end_matches(char::is_numeric).to_string()).collect();
        let c = ph.cell;
        let cell = Cell::new(c[0], c[1], c[2], c[3], c[4], c[5]).unwrap();
        let (mut gamma, mut x_star) = (0, 0);
        for r in drs::reflections::generate(&g, &cell, 1.2, 1e9, None) {
            let pr = parent.reflection(r.hkl).unwrap();
            let back: [f64; 3] = [0, 1, 2].map(|j| (0..3).map(|i| pr.hkl[i] * parent.basis[i][j]).sum());
            if (0..3).any(|i| (back[i] - r.hkl[i] as f64).abs() > 1e-9) {
                fails.push(format!("{name}: {:?} → {:?} → {back:?}", r.hkl, pr.hkl));
            }
            let integral = pr.hkl.iter().all(|v| (v - v.round()).abs() < 1e-9);
            let allowed = integral && !drs::reflections::is_absent(&fm3m, pr.hkl.map(|v| v.round() as i32));
            let is_gamma = pr.k == [0.0; 3];
            let is_x = pr.k.iter().filter(|v| v.abs() == 1.0).count() == 1 && pr.k.iter().filter(|&&v| v == 0.0).count() == 2;
            if is_gamma != allowed || !(is_gamma || is_x && labels.contains(&"X")) {
                fails.push(format!("{name}: {:?} is {:?} in the parent, k = {:?}", r.hkl, pr.hkl, pr.k));
            }
            gamma += usize::from(is_gamma);
            x_star += usize::from(is_x);
        }
        println!("{name}: irreps {irreps:?}; reflections to d = 1.2 Å: {gamma} of the parent's lattice (Γ), {x_star} of the star of X");
        if irreps.iter().any(|i| !labels.contains(&i.as_str())) || gamma == 0 || (labels.contains(&"X") && x_star == 0) {
            fails.push(format!("{name}: irreps {irreps:?}, {gamma} at Γ, {x_star} at X"));
        }
    }
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}
