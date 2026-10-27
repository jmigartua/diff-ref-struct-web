//! Properties of structure factors in every one of the 530 settings (test `structure_factor_properties`).
//!
//! Each setting gets a structure of two atoms: Fe on the general position (coordinates from a fixed pseudo-random
//! sequence, Biso 0.7 Å², arbitrary betas) and O on the most special Wyckoff position of the table (its free
//! coordinates from the same sequence, betas obeying the site's constraints, Occ = m/M). For every h with |h_i| ≤ 3,
//! and (1 3 5), (2 4 6):
//! - absences: F = 0 (|F| < 10⁻⁹ of Σ Occ·M·|f|) for every reflection the core calls absent, X-rays with f′ and f″ and
//!   anisotropic betas included;
//! - Friedel's law: without anomalous terms F(−h) = F(h)* (so |F(−h)| = |F(h)|), to 10⁻⁹; with them, in every
//!   non-centrosymmetric setting some pair differs by more than 10⁻⁶ (Fe and O have different f″/f), and in every
//!   centrosymmetric one none does;
//! - a centrosymmetric group with −1 at the origin gives real F without anomalous terms (|Im F| < 10⁻⁹);
//! - equivalent reflections h·R have the same |F| without anomalous terms (betas carried through R·β·Rᵀ);
//! - the occupancy convention: F equals that of the same atoms written out in P 1 (every distinct position of each orbit
//!   with Occ 1 and its rotated betas), which checks that FullProf's Occ = m/M counts each position once.

use drs::cell::Cell;
use drs::reflections::is_absent;
use drs::scattering::{Complex, Radiation};
use drs::sgtable::SETTINGS;
use drs::site::{eval_triplet, Site};
use drs::spacegroup::SpaceGroup;
use drs::structure_factor::{occupancy_from_chemical, Atom, Structure};
use drs::symop::DEN;
use drs::wyckoff_table::positions;

/// A cell of the right metric for a setting (unique axis and H/R axes from the setting choice).
fn cell_for(number: u16, choice: &str) -> Cell {
    let (a, b, c) = (6.1, 7.3, 8.9);
    let p = match number {
        1..=2 => (a, b, c, 95.0, 100.0, 105.0),
        3..=15 => match choice.trim_start_matches('-').chars().next() {
            Some('c') => (a, b, c, 90.0, 90.0, 103.0),
            Some('a') => (a, b, c, 103.0, 90.0, 90.0),
            _ => (a, b, c, 90.0, 103.0, 90.0),
        },
        16..=74 => (a, b, c, 90.0, 90.0, 90.0),
        75..=142 => (a, a, c, 90.0, 90.0, 90.0),
        143..=194 if choice == "R" => (a, a, a, 72.0, 72.0, 72.0),
        143..=194 => (a, a, c, 90.0, 90.0, 120.0),
        _ => (a, a, a, 90.0, 90.0, 90.0),
    };
    Cell::new(p.0, p.1, p.2, p.3, p.4, p.5).unwrap()
}

/// A fixed pseudo-random sequence in [0, 1) (a linear congruential generator), so the test is reproducible.
struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// R·β·Rᵀ for betas ordered 11 22 33 12 13 23.
fn rotate_beta(r: &[[i32; 3]; 3], b: [f64; 6]) -> [f64; 6] {
    let m = [[b[0], b[3], b[4]], [b[3], b[1], b[5]], [b[4], b[5], b[2]]];
    let mut o = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                for l in 0..3 {
                    o[i][j] += r[i][k] as f64 * m[k][l] * r[j][l] as f64;
                }
            }
        }
    }
    [o[0][0], o[1][1], o[2][2], o[0][1], o[0][2], o[1][2]]
}

/// The atoms written out in P 1: every distinct position of each orbit, Occ = chemical occupancy, rotated betas.
fn expand(st: &Structure, chemical: &[f64]) -> Vec<Atom> {
    let mut out = Vec::new();
    for (a, &occ) in st.atoms.iter().zip(chemical) {
        let mut seen: Vec<[f64; 3]> = Vec::new();
        for op in &st.group.ops {
            let y: Vec<f64> = (0..3)
                .map(|i| ((0..3).map(|j| op.r[i][j] as f64 * a.x[j]).sum::<f64>() + op.t[i] as f64 / DEN as f64).rem_euclid(1.0))
                .collect();
            let y = [y[0], y[1], y[2]];
            let same = |p: &[f64; 3]| (0..3).all(|i| {
                let d = (p[i] - y[i]).abs();
                d.min(1.0 - d) < 1e-6
            });
            if seen.iter().any(same) {
                continue;
            }
            seen.push(y);
            out.push(Atom { x: y, occ, beta: a.beta.map(|b| rotate_beta(&op.r, b)), ..a.clone() });
        }
    }
    out
}

#[test]
fn structure_factor_properties() {
    let mut rng = Lcg(2027);
    let mut fails = Vec::new();
    let p1 = SpaceGroup::from_hall("P 1").unwrap();
    let xray = Radiation::Xray { lambda: 1.5406, anomalous: true };
    let plain = Radiation::Xray { lambda: 1.5406, anomalous: false };
    let (mut n_absent, mut n_checked, mut bijvoet_settings) = (0usize, 0usize, 0usize);
    let lim = 3;
    for s in SETTINGS.iter() {
        let g = SpaceGroup::from_setting(s);
        let tag = format!("#{} {}", s.hall_number, s.hm);
        let cell = cell_for(s.number, s.choice);
        // Fe on the general position, arbitrary betas.
        let x = [rng.next(), rng.next(), rng.next()];
        let beta = [0.004, 0.003, 0.002, 0.0007, -0.0005, 0.0004].map(|v| v * (0.5 + rng.next()));
        // O on the most special position of the table, betas from the site's free parameters.
        let w = positions(s.hall_number).last().unwrap();
        let xo = eval_triplet(w.coord, [rng.next(), rng.next(), rng.next()]).unwrap();
        let site = Site::new(&g, xo).unwrap();
        let free: Vec<f64> = site.adp.free.iter().map(|_| 0.002 * (0.5 + rng.next())).collect();
        let bo = site.adp.apply(&free);
        let atoms = vec![
            Atom { label: "Fe1".into(), species: "Fe".into(), x, biso: 0.7, occ: occupancy_from_chemical(&g, x, 1.0).unwrap(), beta: Some(beta), charge: None },
            Atom {
                label: "O1".into(),
                species: "O".into(),
                x: site.position,
                biso: 0.4,
                occ: occupancy_from_chemical(&g, site.position, 1.0).unwrap(),
                beta: Some([bo[0], bo[1], bo[2], bo[3], bo[4], bo[5]]),
                charge: None,
            },
        ];
        let st = Structure { group: g.clone(), cell, atoms };
        let p1st = Structure { group: p1.clone(), cell, atoms: expand(&st, &[1.0, 1.0]) };
        let sc_x = st.scatterers(xray).unwrap();
        let sc_p = st.scatterers(plain).unwrap();
        let sc_p1 = p1st.scatterers(plain).unwrap();
        let scale: f64 = st.atoms.iter().zip(&sc_x).map(|(a, f)| a.occ * g.ops.len() as f64 * f.at(0.0).abs()).sum();
        let mut bijvoet = false;
        // The box, plus two reflections with three different indices, all odd and all even (in F 4 3 2 and F 41 3 2
        // every allowed reflection of the box has a zero or two equal indices, and is centric).
        let mut list: Vec<[i32; 3]> = Vec::new();
        for h in -lim..=lim {
            for k in -lim..=lim {
                for l in -lim..=lim {
                    if [h, k, l] != [0, 0, 0] {
                        list.push([h, k, l]);
                    }
                }
            }
        }
        list.extend([[1, 3, 5], [2, 4, 6]]);
        for &hkl in &list {
            let [h, k, l] = hkl;
            n_checked += 1;
            let fx = st.structure_factor(&sc_x, hkl);
            if is_absent(&g, hkl) {
                n_absent += 1;
                if fx.abs() > 1e-9 * scale {
                    fails.push(format!("{tag}: {hkl:?} absent, |F| = {:.3e}", fx.abs()));
                }
                continue;
            }
            let minus = [-h, -k, -l];
            let fp = st.structure_factor(&sc_p, hkl);
            let fm = st.structure_factor(&sc_p, minus);
            if (fm - fp.conj()).abs() > 1e-9 * scale {
                fails.push(format!("{tag}: {hkl:?} F(−h) = {fm} is not F(h)* = {}", fp.conj()));
            }
            if (st.structure_factor(&sc_x, minus).abs() - fx.abs()).abs() > 1e-6 * scale {
                bijvoet = true;
            }
            if g.inversion_at_origin && fp.im.abs() > 1e-9 * scale {
                fails.push(format!("{tag}: {hkl:?} −1 at the origin but Im F = {:.3e}", fp.im));
            }
            for op in &g.ops {
                let hr: [i32; 3] = std::array::from_fn(|j| (0..3).map(|i| hkl[i] * op.r[i][j]).sum());
                let f2 = st.structure_factor(&sc_p, hr);
                if (f2.abs() - fp.abs()).abs() > 1e-9 * scale {
                    fails.push(format!("{tag}: |F{hr:?}| = {:.6} against |F{hkl:?}| = {:.6}", f2.abs(), fp.abs()));
                    break;
                }
            }
            let f1: Complex = p1st.structure_factor(&sc_p1, hkl);
            if (f1 - fp).abs() > 1e-9 * scale {
                fails.push(format!("{tag}: {hkl:?} F = {fp} in the group, {f1} written out in P 1"));
            }
        }
        if bijvoet {
            bijvoet_settings += 1;
        }
        if bijvoet == g.centrosymmetric {
            fails.push(format!("{tag}: centrosymmetric {}, but Bijvoet differences found: {bijvoet}", g.centrosymmetric));
        }
    }
    let shown: Vec<_> = fails.iter().take(40).cloned().collect();
    assert!(fails.is_empty(), "{} failures:\n{}", fails.len(), shown.join("\n"));
    println!("{n_checked} reflections in 530 settings, {n_absent} absent; Bijvoet differences in {bijvoet_settings} non-centrosymmetric settings");
}
