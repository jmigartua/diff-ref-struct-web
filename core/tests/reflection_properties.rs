//! Properties of the reflections in every one of the 530 settings (test `reflection_properties`).
//!
//! - Absences: for an atom on a general position (coordinates from a fixed pseudo-random sequence), the structure
//!   factor F(h) = Σ exp(2πi h·(R·x + t)) over all the operators, computed directly, is zero (|F| < 10⁻⁹) for every
//!   reflection the core calls absent, and not zero (|F| > 10⁻⁶) for every other one, for all h with |h_i| ≤ 5.
//! - Multiplicities: in a sphere 1/d ≤ 1/dmin, the multiplicities of the core's unique reflections add up to the
//!   number of allowed reflections counted one by one, and every allowed reflection is equivalent to one of them.
//! - FullProf's region rule (CrysFML `Asu_Hkl`) picks exactly one representative per set of equivalents in the 230
//!   standard settings; the settings where it does not (non-standard orientations, where the fallback is used) are
//!   counted and printed.

use drs::cell::Cell;
use drs::reflections::{equivalents, generate, is_absent};
use drs::sgtable::SETTINGS;
use drs::spacegroup::SpaceGroup;
use std::collections::BTreeSet;

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

#[test]
fn reflection_properties() {
    let mut fails = Vec::new();
    let mut fallback_settings = Vec::new();
    let mut rng = Lcg(2026);
    let (mut n_absent, mut n_present) = (0usize, 0usize);
    let first_of_type: BTreeSet<u16> = (1..=230).map(|n| SETTINGS.iter().find(|s| s.number == n).unwrap().hall_number).collect();
    for s in SETTINGS.iter() {
        let g = SpaceGroup::from_setting(s);
        let tag = format!("#{} {}", s.hall_number, s.hm);
        // Absences against a direct structure factor.
        let x = [rng.next(), rng.next(), rng.next()];
        let lim = 5;
        for h in -lim..=lim {
            for k in -lim..=lim {
                for l in -lim..=lim {
                    let hkl = [h, k, l];
                    if hkl == [0, 0, 0] {
                        continue;
                    }
                    let (mut re, mut im) = (0.0f64, 0.0f64);
                    for op in &g.ops {
                        let phase: f64 = (0..3)
                            .map(|i| hkl[i] as f64 * ((0..3).map(|j| op.r[i][j] as f64 * x[j]).sum::<f64>() + op.t[i] as f64 / 12.0))
                            .sum();
                        re += (2.0 * std::f64::consts::PI * phase).cos();
                        im += (2.0 * std::f64::consts::PI * phase).sin();
                    }
                    let f = re.hypot(im);
                    if is_absent(&g, hkl) {
                        n_absent += 1;
                        if f > 1e-9 {
                            fails.push(format!("{tag}: {hkl:?} called absent, |F| = {f:.3e}"));
                        }
                    } else {
                        n_present += 1;
                        if f < 1e-6 {
                            fails.push(format!("{tag}: {hkl:?} allowed, but |F| = {f:.3e}"));
                        }
                    }
                }
            }
        }
        // Multiplicities against a count of every allowed reflection in a sphere.
        let cell = cell_for(s.number, s.choice);
        let dmin = 1.4;
        let list = generate(&g, &cell, dmin, 1e6, None);
        let sum: usize = list.iter().map(|r| r.multiplicity).sum();
        let reps: BTreeSet<[i32; 3]> = list.iter().map(|r| r.hkl).collect();
        let hmax = [cell.a, cell.b, cell.c].map(|v| (v / dmin) as i32 + 1);
        let mut count = 0;
        for h in -hmax[0]..=hmax[0] {
            for k in -hmax[1]..=hmax[1] {
                for l in -hmax[2]..=hmax[2] {
                    let hkl = [h, k, l];
                    if hkl == [0, 0, 0] || is_absent(&g, hkl) {
                        continue;
                    }
                    if cell.d_spacing([h as f64, k as f64, l as f64]).unwrap() >= dmin {
                        count += 1;
                        if equivalents(&g, hkl).intersection(&reps).count() != 1 {
                            fails.push(format!("{tag}: {hkl:?} not equivalent to exactly one listed reflection"));
                        }
                    }
                }
            }
        }
        if sum != count {
            fails.push(format!("{tag}: multiplicities add up to {sum}, {count} allowed reflections in the sphere"));
        }
        let fallback = list.iter().filter(|r| !r.fullprof_choice).count();
        if fallback > 0 {
            if first_of_type.contains(&s.hall_number) {
                fails.push(format!("{tag} (standard setting): {fallback} reflections not chosen by FullProf's rule"));
            }
            fallback_settings.push(format!("{} {}", s.hall_number, s.hm));
        }
    }
    let shown: Vec<_> = fails.iter().take(40).cloned().collect();
    assert!(fails.is_empty(), "{} failures:\n{}", fails.len(), shown.join("\n"));
    println!("{n_absent} absent and {n_present} allowed reflections checked against |F|");
    println!("settings where FullProf's region rule needs the fallback ({}): {}", fallback_settings.len(), fallback_settings.join(", "));
}
