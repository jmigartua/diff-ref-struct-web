//! Structure factors against FullProf (test `fullprof_structure_factors`).
//!
//! `tests/data/fullprof_structure_factors.txt` holds three structures and the F² FullProf printed for them, copied by
//! `tests/data/make_fullprof_structure_factors.py` with file and line from `manuals/fullprof/`:
//! - PPH3, P 21/c, 34 atoms, X-rays λ = 1.001 Å: 7538 reflections of `PPH3-sls.hkl`;
//! - citric acid, P 21/c, 21 atoms, X-rays λ = 1.001 Å: 4074 reflections of `citric-acid-sls_1.hkl`;
//! - LaMnO3, P n m a, atoms on 4c, 4b, 8d and 4c (Occ = m/M: 0.5, 0.5, 1, 0.5), neutrons: 414 reflections of
//!   `Lmo_san.out`.
//!
//! The core computes F with its own tables (CrysFML's, `scattering_table`) and must give, for every reflection:
//! - the multiplicity FullProf printed, and sinθ/λ to the printed five decimals (half a unit, plus 10⁻⁶ of s for the
//!   cell, printed to five or six decimals and used by FullProf in single precision);
//! - F² within the uncertainty that FullProf's printing leaves: the atoms are printed to 4–5 decimals (rigid bodies
//!   and symmetry modes compute them to more), so the bound carries half a unit of the last printed decimal of every
//!   x, y, z, Biso and Occ through F, plus the difference between the core's form factors and the coefficients FullProf
//!   printed (four decimals), plus half a unit of the last decimal of F², plus 10⁻⁶ of F² for FullProf's single
//!   precision. The worst ratio to the bound is printed.
//!
//! It also checks the core's tables against what FullProf printed: the X-ray coefficients to half a unit of the fourth
//! decimal, and f′, f″ (columns Dfp, Dfpp) exactly as the table of the anode nearest to λ (Mo at λ = 1.001 Å).

use drs::cell::Cell;
use drs::reflections::equivalents;
use drs::scattering::{Radiation, Scatterer};
use drs::spacegroup::SpaceGroup;
use drs::structure_factor::{Atom, Structure};
use std::f64::consts::PI;

/// A printed number and half a unit of its last decimal.
fn printed(s: &str) -> (f64, f64) {
    let dec = s.split_once('.').map_or(0, |(_, f)| f.len() as i32);
    (s.parse().unwrap(), 0.5 * 10f64.powi(-dec))
}

/// A reflection: h k l, multiplicity, sinθ/λ and F² (each value with half a unit of its last printed decimal), where.
type Ref = ([i32; 3], usize, (f64, f64), (f64, f64), String);

struct Phase {
    name: String,
    structure: Structure,
    radiation: Radiation,
    /// Half units of the printed x, y, z, Biso, Occ of each atom.
    rounding: Vec<[f64; 5]>,
    /// FullProf's printed X-ray coefficients by species: a1 b1 … a4 b4 c f0 Z dfp dfpp.
    coef: Vec<(String, Vec<f64>, String)>,
    /// h k l, mult, sinθ/λ, F² (value, half unit), where.
    refs: Vec<Ref>,
}

fn read() -> Vec<Phase> {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/fullprof_structure_factors.txt");
    let text = std::fs::read_to_string(path).unwrap();
    let mut phases: Vec<Phase> = Vec::new();
    let mut files = (String::new(), String::new());
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = line.split_whitespace().collect();
        match f[0] {
            "PHASE" => {
                let (_, r) = line.split_once('"').unwrap();
                let (hall, r) = r.split_once('"').unwrap();
                let g: Vec<&str> = r.split_whitespace().collect();
                let c: Vec<f64> = g[..6].iter().map(|x| x.parse().unwrap()).collect();
                let lambda: f64 = g[7].parse().unwrap();
                let radiation = if g[6] == "xray" { Radiation::Xray { lambda, anomalous: true } } else { Radiation::Neutron };
                files = (g[8].to_string(), g[9].to_string());
                phases.push(Phase {
                    name: f[1].to_string(),
                    structure: Structure {
                        group: SpaceGroup::from_hall(hall).unwrap(),
                        cell: Cell::new(c[0], c[1], c[2], c[3], c[4], c[5]).unwrap(),
                        atoms: Vec::new(),
                    },
                    radiation,
                    rounding: Vec::new(),
                    coef: Vec::new(),
                    refs: Vec::new(),
                });
            }
            "COEF" => {
                let p = phases.last_mut().unwrap();
                p.coef.push((f[1].to_string(), f[2..15].iter().map(|x| x.parse().unwrap()).collect(), format!("{}:{}", files.0, f[15])));
            }
            "ATOM" => {
                let p = phases.last_mut().unwrap();
                let v: Vec<(f64, f64)> = f[3..8].iter().map(|x| printed(x)).collect();
                p.structure.atoms.push(Atom {
                    label: f[1].to_string(),
                    species: f[2].to_string(),
                    x: [v[0].0, v[1].0, v[2].0],
                    biso: v[3].0,
                    occ: v[4].0,
                    beta: None,
                    charge: None,
                });
                p.rounding.push([v[0].1, v[1].1, v[2].1, v[3].1, v[4].1]);
            }
            "REF" => {
                let p = phases.last_mut().unwrap();
                p.refs.push((
                    [f[1].parse().unwrap(), f[2].parse().unwrap(), f[3].parse().unwrap()],
                    f[4].parse().unwrap(),
                    printed(f[5]),
                    printed(f[6]),
                    format!("{}:{}", files.1, f[7]),
                ));
            }
            _ => panic!("unknown line {line}"),
        }
    }
    phases
}

#[test]
fn fullprof_structure_factors() {
    let phases = read();
    assert_eq!(phases.len(), 3);
    let mut fails = Vec::new();
    let mut total = 0;
    for p in &phases {
        let st = &p.structure;
        let sc = st.scatterers(p.radiation).unwrap();
        // The tables against FullProf's printed coefficients.
        for (species, c, at) in &p.coef {
            let Scatterer::Xray { a, b, c: c0, fp, fpp, .. } = Scatterer::lookup(species, p.radiation).unwrap() else { panic!() };
            let mine = [a[0], b[0], a[1], b[1], a[2], b[2], a[3], b[3], c0];
            for (k, (m, v)) in mine.iter().zip(c).enumerate() {
                if (m - v).abs() > 0.5e-4 + 1e-12 {
                    fails.push(format!("{at}: {species} coefficient {k}: core {m}, FullProf {v}"));
                }
            }
            if (fp - c[11]).abs() > 1e-9 || (fpp - c[12]).abs() > 1e-9 {
                fails.push(format!("{at}: {species} f′, f″: core {fp} {fpp}, FullProf {} {}", c[11], c[12]));
            }
        }
        // FullProf's printed form factors, to bound the difference with the core's.
        let printed_f = |species: &str, s: f64| -> Option<f64> {
            p.coef.iter().find(|c| c.0.eq_ignore_ascii_case(species)).map(|(_, c, _)| {
                (0..4).map(|i| c[2 * i] * (-c[2 * i + 1] * s * s).exp()).sum::<f64>() + c[8] + c[11]
            })
        };
        // Controls: the same comparison without f′ and f″ (X-rays), or with chemical occupancies in place of FullProf's
        // Occ = m/M (LaMnO3), must fail for more than a tenth of the reflections, so that the test tells the conventions
        // apart (C, H, O of citric acid have small f′ and f″ at λ = 1.001 Å: 0.002 to 0.008 electrons).
        let no_anomalous = match p.radiation {
            Radiation::Xray { lambda, .. } => Some(st.scatterers(Radiation::Xray { lambda, anomalous: false }).unwrap()),
            Radiation::Neutron => None,
        };
        let mut chemical = st.clone();
        for a in chemical.atoms.iter_mut() {
            a.occ = drs::structure_factor::chemical_from_occupancy(&st.group, a.x, a.occ).unwrap();
        }
        let mut control_fails = 0;
        let mut worst = (0.0f64, String::new());
        for (hkl, mult, (s_fp, s_half), (f2_fp, f2_half), at) in &p.refs {
            total += 1;
            let m = equivalents(&st.group, *hkl).len();
            if m != *mult {
                fails.push(format!("{at}: {hkl:?} multiplicity {m}, FullProf {mult}"));
            }
            let s = st.sin_theta_over_lambda(*hkl);
            if (s - s_fp).abs() > s_half + 1e-6 * s {
                fails.push(format!("{at}: {hkl:?} sinθ/λ {s:.6}, FullProf {s_fp}"));
            }
            let f = st.structure_factor(&sc, *hkl);
            // Bound on |δF| from the printed precision of the atoms and the form factors.
            let mut df = 0.0;
            for ((atom, r), scj) in st.atoms.iter().zip(&p.rounding).zip(&sc) {
                let fj = scj.at(s).abs();
                let dfj = printed_f(&atom.species, s).map_or(0.0, |v| (v - scj.at(s).re).abs());
                let w = (-atom.biso * s * s).exp();
                let mut pos = 0.0;
                for op in &st.group.ops {
                    let hr: Vec<f64> = (0..3).map(|j| (0..3).map(|i| hkl[i] as f64 * op.r[i][j] as f64).sum()).collect();
                    pos += 2.0 * PI * (0..3).map(|j| hr[j].abs() * r[j]).sum::<f64>();
                }
                let n = st.group.ops.len() as f64;
                df += atom.occ * w * (fj * (pos + n * s * s * r[3]) + n * dfj) + r[4] * w * fj * n;
            }
            let fa = f.abs();
            let tol = 2.0 * fa * df + df * df + f2_half + 1e-6 * f2_fp;
            let d = (f.norm_sqr() - f2_fp).abs();
            if d / tol > worst.0 {
                worst = (d / tol, format!("{hkl:?} at {at}"));
            }
            if d > tol {
                fails.push(format!("{at}: {hkl:?} F² = {:.4}, FullProf {f2_fp} (bound {tol:.4})", f.norm_sqr()));
            }
            let control = match &no_anomalous {
                Some(sc0) => st.structure_factor(sc0, *hkl),
                None => chemical.structure_factor(&sc, *hkl),
            };
            if (control.norm_sqr() - f2_fp).abs() > tol {
                control_fails += 1;
            }
        }
        println!("{}: {} reflections, worst |ΔF²| at {:.3} of the bound ({})", p.name, p.refs.len(), worst.0, worst.1);
        let what = if no_anomalous.is_some() { "without f′ and f″" } else { "with chemical occupancies" };
        println!("  control {what}: {control_fails} reflections outside the bound");
        assert!(control_fails * 10 > p.refs.len(), "{}: the control {what} should fail for more than a tenth of the reflections", p.name);
    }
    let shown: Vec<_> = fails.iter().take(40).cloned().collect();
    assert!(fails.is_empty(), "{} disagreements:\n{}", fails.len(), shown.join("\n"));
    assert_eq!(total, 7538 + 4074 + 414, "reflections compared");
}
