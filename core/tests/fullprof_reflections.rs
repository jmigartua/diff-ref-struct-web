//! Reflections against FullProf (test `fullprof_reflections`).
//!
//! `tests/data/fullprof_reflections.txt` holds the reflection lists FullProf printed in the `.out` files of
//! `manuals/fullprof/` (9 phases, 14337 reflections, in P n m a, P 21/c, P 21/a and P 2/m), with each phase's Hall
//! symbol, refined cell and wavelength, written by `tests/data/make_fullprof_reflections.py` with file and line.
//!
//! For each phase the core generates the reflections between FullProf's smallest and largest d, from FullProf's
//! cell and wavelength, and must give, row for row in FullProf's order:
//! - the same h k l (FullProf's choice of representative, not only an equivalent one) and multiplicity;
//! - d to the printed six decimals and 2θ to the printed three (half a unit of the last digit, plus 2·10⁻⁷ of d for
//!   the cell, printed to six decimals and computed by FullProf in single precision);
//! - no reflection more or fewer: the same systematic absences.
//!
//! Two rows may come in the other order only where FullProf's own d values for them agree to 2·10⁻⁶ Å (FullProf
//! sorts single-precision values; such pairs are equal to its precision), or where FullProf's printed d values are
//! themselves out of order for the pair (FullProf's list is then not sorted by its own d; such pairs are listed).
//! The 2θ tolerance carries the d tolerance through 2θ: Δ2θ = 2·tanθ·Δd/d.

use drs::cell::Cell;
use drs::reflections::generate;
use drs::spacegroup::SpaceGroup;

/// A phase: Hall symbol, cell, wavelength, file, and its reflections.
type Phase = (String, Vec<f64>, f64, String, Vec<Row>);

struct Row {
    hkl: [i32; 3],
    mult: usize,
    tt: (f64, i32),
    d: (f64, i32),
    at: String,
}

fn printed(s: &str) -> (f64, i32) {
    (s.parse().unwrap(), s.split_once('.').map_or(0, |(_, f)| f.len() as i32))
}

#[test]
fn fullprof_reflections() {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/fullprof_reflections.txt");
    let text = std::fs::read_to_string(path).unwrap();
    let mut phases: Vec<Phase> = Vec::new();
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        if let Some(rest) = line.strip_prefix("PHASE ") {
            let (_, r) = rest.split_once('"').unwrap();
            let (hall, r) = r.split_once('"').unwrap();
            let f: Vec<&str> = r.split_whitespace().collect();
            let cell: Vec<f64> = f[..6].iter().map(|x| x.parse().unwrap()).collect();
            phases.push((hall.to_string(), cell, f[6].parse().unwrap(), f[7].to_string(), Vec::new()));
        } else {
            let f: Vec<&str> = line.split_whitespace().collect();
            phases.last_mut().unwrap().4.push(Row {
                hkl: [f[1].parse().unwrap(), f[2].parse().unwrap(), f[3].parse().unwrap()],
                mult: f[4].parse().unwrap(),
                tt: printed(f[5]),
                d: printed(f[6]),
                at: f[7].to_string(),
            });
        }
    }
    assert_eq!(phases.len(), 9, "phases");
    let mut fails = Vec::new();
    let (mut total, mut swaps, mut worst_d, mut worst_tt) = (0, 0, 0.0f64, 0.0f64);
    let mut unsorted: Vec<String> = Vec::new();
    for (hall, c, lambda, file, rows) in &phases {
        let g = SpaceGroup::from_hall(hall).unwrap();
        let cell = Cell::new(c[0], c[1], c[2], c[3], c[4], c[5]).unwrap();
        let dmin = rows.iter().map(|r| r.d.0).fold(f64::MAX, f64::min) - 2e-6;
        let dmax = rows.iter().map(|r| r.d.0).fold(0.0, f64::max) + 2e-6;
        let mine = generate(&g, &cell, dmin, dmax, Some(*lambda));
        total += rows.len();
        if mine.len() != rows.len() {
            fails.push(format!("{file} ({hall}): core {} reflections, FullProf {}", mine.len(), rows.len()));
        }
        let mut used = vec![false; mine.len()];
        for (i, r) in rows.iter().enumerate() {
            // The core's row at the same place, or a neighbour with the same d to FullProf's precision.
            let j = if i < mine.len() && mine[i].hkl == r.hkl {
                Some(i)
            } else {
                let near = (i.saturating_sub(3)..(i + 4).min(mine.len())).find(|&j| !used[j] && mine[j].hkl == r.hkl);
                if let Some(j) = near {
                    let other = rows.get(j).map_or(f64::MAX, |o| (o.d.0 - r.d.0).abs());
                    // FullProf's own order contradicts its printed d: the earlier of the two rows has the smaller d.
                    let (lo, hi) = (i.min(j), i.max(j));
                    let fp_unsorted = rows[lo].d.0 < rows[hi].d.0;
                    if other <= 2e-6 {
                        swaps += 1;
                    } else if fp_unsorted {
                        unsorted.push(format!("{} {:?} d {} before {} {:?} d {}", rows[lo].at, rows[lo].hkl, rows[lo].d.0, rows[hi].at, rows[hi].hkl, rows[hi].d.0));
                    } else {
                        fails.push(format!("{}: ({:?}) out of FullProf's order", r.at, r.hkl));
                    }
                }
                near
            };
            let Some(j) = j else {
                fails.push(format!("{}: FullProf's {:?} not in the core's list (core row: {:?})", r.at, r.hkl, mine.get(i).map(|m| m.hkl)));
                continue;
            };
            used[j] = true;
            let m = &mine[j];
            if m.multiplicity != r.mult {
                fails.push(format!("{}: {:?} multiplicity {} against FullProf's {}", r.at, r.hkl, m.multiplicity, r.mult));
            }
            if !m.fullprof_choice {
                fails.push(format!("{}: {:?} chosen by the fallback, not FullProf's rule", r.at, r.hkl));
            }
            let dd = (m.d - r.d.0).abs();
            let tol_d = 0.5 * 10f64.powi(-r.d.1) + 2e-7 * r.d.0;
            worst_d = worst_d.max(dd / tol_d);
            if dd > tol_d {
                fails.push(format!("{}: {:?} d = {:.7} against FullProf's {}", r.at, r.hkl, m.d, r.d.0));
            }
            let tt = m.two_theta.unwrap_or(f64::NAN);
            let tan = (r.tt.0.to_radians() / 2.0).tan();
            let tol_t = 0.5 * 10f64.powi(-r.tt.1) + 2.0 * tan * 2e-7 * 180.0 / std::f64::consts::PI;
            worst_tt = worst_tt.max((tt - r.tt.0).abs() / tol_t);
            if tt.is_nan() || (tt - r.tt.0).abs() > tol_t {
                fails.push(format!("{}: {:?} 2θ = {tt:.4} against FullProf's {}", r.at, r.hkl, r.tt.0));
            }
        }
    }
    let shown: Vec<_> = fails.iter().take(40).cloned().collect();
    assert!(fails.is_empty(), "{} disagreements:\n{}", fails.len(), shown.join("\n"));
    assert_eq!(total, 14337, "reflections compared");
    println!("{total} reflections, {swaps} pairs in the other order at equal d; worst d and 2θ at {worst_d:.2} and {worst_tt:.2} of the tolerance");
    unsorted.dedup(); // the pair is met from both of its rows
    println!("pairs that FullProf lists against its own printed d ({}): {}", unsorted.len(), unsorted.join("; "));
    assert!(unsorted.len() <= 1, "more pairs out of FullProf's own order than the one known (PPH3-sls_1.out:1763)");
}
