//! Reflections of a phase: the unique reflections of the Laue class between two d-spacings, their multiplicity,
//! d and 2θ, with the systematic absences of the space group (phase 1, item 4).
//!
//! As FullProf does it, in CrysFML's `CFML_Reflections_Utilities.f90` (J. Rodríguez-Carvajal, ILL, LGPL-3; read at
//! https://github.com/easyscience/crysfml, commit 8ed8670):
//! - **absences** come from the operators themselves (`Hkl_Absent`): h is absent when an operator {R|t} of the group
//!   has h·R = h and h·t not an integer, which covers centring, screw axes and glide planes. Reflections forbidden
//!   only by atoms on special positions are not absent here;
//! - **one reflection per set of equivalents** (`Hkl_Uni`): the member h of the set {h·R}, Friedel mates included,
//!   that lies in the index box `Hkl_Uni` loops over (h, k or l ≥ 0, by space-group number and unique axis) and for
//!   which `Asu_Hkl(h) = h`, the region of reciprocal space chosen for each Laue class. Both are ported here
//!   (`fullprof_box`, `fullprof_asu`); the box matters, since in −3m1 the region alone admits two members of some
//!   sets, e.g. (0 4 1) and (4 0 −1), and only the first is in the box (l ≥ 0). In settings
//!   where they do not single out one member (non-standard orientations), the representative is the member that
//!   comes first in a fixed order (largest h, then k, then l), and `Reflection::fullprof_choice` is false;
//! - **multiplicity** (`Hkl_Mult`): the number of distinct reflections h·R and −h·R;
//! - **order**: by increasing sinθ/λ = 1/(2d), and for equal d by h, k, l.

use crate::cell::Cell;
use crate::spacegroup::SpaceGroup;
use crate::symop::DEN;
use std::collections::BTreeSet;

/// One unique reflection.
#[derive(Clone, Debug, PartialEq)]
pub struct Reflection {
    pub hkl: [i32; 3],
    /// Number of equivalent reflections, Friedel mates included.
    pub multiplicity: usize,
    /// d-spacing, Å.
    pub d: f64,
    /// 2θ in degrees for the wavelength asked, if the reflection is within reach.
    pub two_theta: Option<f64>,
    /// Whether the representative is the one FullProf's region rule picks (false: the fallback order).
    pub fullprof_choice: bool,
}

/// h·R (h as a row vector): the reflection equivalent to h under the rotation R of an operator.
fn h_times_r(h: [i32; 3], r: &[[i32; 3]; 3]) -> [i32; 3] {
    [0, 1, 2].map(|j| (0..3).map(|i| h[i] * r[i][j]).sum())
}

/// Systematically absent: some operator {R|t} with h·R = h has h·t ≠ 0 (mod 1).
pub fn is_absent(g: &SpaceGroup, h: [i32; 3]) -> bool {
    g.ops.iter().any(|op| {
        h_times_r(h, &op.r) == h && (0..3).map(|i| h[i] * op.t[i]).sum::<i32>().rem_euclid(DEN) != 0
    })
}

/// The distinct reflections equivalent to h in the Laue class: h·R and −h·R for every rotation R.
pub fn equivalents(g: &SpaceGroup, h: [i32; 3]) -> BTreeSet<[i32; 3]> {
    let mut out = BTreeSet::new();
    for op in &g.ops {
        let k = h_times_r(h, &op.r);
        out.insert(k);
        out.insert(k.map(|v| -v));
    }
    out
}

/// The Laue-class rule FullProf uses for a space group (CrysFML `Asu_Hkl`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AsuRule {
    Triclinic,
    Monoclinic(char),
    Orthorhombic,
    Tetragonal4m,
    Tetragonal4mmm,
    Trigonal3,
    Trigonal31m,
    Trigonal3m1,
    Hexagonal6m,
    Hexagonal6mmm,
    Cubicm3,
    Cubicm3m,
}

/// The rule for the group: by its number, and for the monoclinic types by the unique axis of the setting
/// (CrysFML reads it from the setting information; here from the table's setting choice, b by default).
pub fn asu_rule(g: &SpaceGroup) -> AsuRule {
    let number = g.setting.map_or(0, |s| s.number);
    match number {
        1..=2 => AsuRule::Triclinic,
        // The unique axis is the first letter of the setting choice, after a '-' if there is one ("-c2"), as CrysFML
        // reads its setting information.
        3..=15 => AsuRule::Monoclinic(
            g.setting
                .and_then(|s| s.choice.trim_start_matches('-').chars().next())
                .filter(|c| matches!(c, 'a' | 'b' | 'c'))
                .unwrap_or('b'),
        ),
        16..=74 => AsuRule::Orthorhombic,
        75..=88 => AsuRule::Tetragonal4m,
        89..=142 => AsuRule::Tetragonal4mmm,
        143..=148 => AsuRule::Trigonal3,
        149 | 151 | 153 | 157 | 159 | 162 | 163 => AsuRule::Trigonal31m,
        150..=167 => AsuRule::Trigonal3m1,
        168..=176 => AsuRule::Hexagonal6m,
        177..=194 => AsuRule::Hexagonal6mmm,
        195..=206 => AsuRule::Cubicm3,
        207..=230 => AsuRule::Cubicm3m,
        // A group outside the table: by its Laue class.
        _ => match g.laue {
            "-1" => AsuRule::Triclinic,
            "2/m" => AsuRule::Monoclinic('b'),
            "mmm" => AsuRule::Orthorhombic,
            "4/m" => AsuRule::Tetragonal4m,
            "4/mmm" => AsuRule::Tetragonal4mmm,
            "-3" => AsuRule::Trigonal3,
            "-3m" => AsuRule::Trigonal3m1,
            "6/m" => AsuRule::Hexagonal6m,
            "6/mmm" => AsuRule::Hexagonal6mmm,
            "m-3" => AsuRule::Cubicm3,
            _ => AsuRule::Cubicm3m,
        },
    }
}

/// h or −h by CrysFML's sign normalisation: the first of the components (in the order given) that is not zero is
/// made positive; (0, 0, 0) stays.
fn normalise(h: [i32; 3], order: [usize; 3]) -> [i32; 3] {
    for &i in &order {
        if h[i] != 0 {
            return if h[i] > 0 { h } else { h.map(|v| -v) };
        }
    }
    h
}

/// CrysFML's `Asu_Hkl`: Some(hh) when h, after the sign normalisation, lies in the Laue class's region.
pub fn fullprof_asu(h: [i32; 3], rule: AsuRule) -> Option<[i32; 3]> {
    use AsuRule::*;
    let [a, b, c];
    match rule {
        Triclinic => Some(normalise(h, [2, 0, 1])),
        Orthorhombic => Some(h.map(i32::abs)),
        Monoclinic('c') => {
            let hh = normalise(h, [2, 1, 0]);
            [a, b, c] = hh;
            (if b == 0 { a >= 0 } else { b >= 0 && c >= 0 }).then_some(hh)
        }
        Monoclinic('a') => {
            let hh = normalise(h, [0, 2, 1]);
            [a, b, c] = hh;
            (if a == 0 { b >= 0 } else { a >= 0 && c >= 0 }).then_some(hh)
        }
        Monoclinic(_) => {
            let hh = normalise(h, [2, 1, 0]);
            [a, b, c] = hh;
            (if c == 0 { a >= 0 } else { b >= 0 && c >= 0 }).then_some(hh)
        }
        _ => {
            let hh = normalise(h, [0, 1, 2]);
            [a, b, c] = hh;
            let ok = match rule {
                Tetragonal4m => (a == 0 && b >= 0 && c >= 0) || (a > 0 && b > 0 && c >= 0),
                Tetragonal4mmm => a >= 0 && b >= 0 && c >= 0 && a >= b,
                Trigonal3 => (a == 0 && b == 0 && c > 0) || (a + b > 0 && c > 0) || (a > 0 && b >= 0 && c == 0),
                Trigonal3m1 => (a >= b && b >= 0) || (a >= 0 && b > 0 && c > 0) || (a >= 0 && b == a && c >= 0),
                Trigonal31m => (a >= b && a >= 0 && b > 0) || (a >= 0 && b == 0 && c >= 0),
                Hexagonal6m => (a > 0 && b > 0 && c >= 0) || (a == 0 && b >= 0 && c >= 0),
                Hexagonal6mmm => b >= 0 && a >= b && c >= 0,
                Cubicm3 => (c >= 0 && a >= c && b == c) || (c >= 0 && a > c && b > c),
                Cubicm3m => c >= 0 && b >= c && a >= b,
                _ => unreachable!(),
            };
            ok.then_some(hh)
        }
    }
}

/// The index box of CrysFML's `Hkl_Uni` loop for a rule: the signs h, k, l may take.
pub fn fullprof_box(h: [i32; 3], rule: AsuRule) -> bool {
    use AsuRule::*;
    let [a, b, c] = h;
    match rule {
        Triclinic => c >= 0,
        Monoclinic('a') => a >= 0 && c >= 0,
        Monoclinic(_) => b >= 0 && c >= 0,
        Trigonal3 => a >= 0 && c >= 0,
        Trigonal31m => a >= 0 && b >= 0,
        _ => a >= 0 && b >= 0 && c >= 0,
    }
}

/// The representative of the set of equivalents of h, and whether FullProf's rule picked it.
pub fn representative(g: &SpaceGroup, h: [i32; 3], rule: AsuRule) -> ([i32; 3], bool) {
    let eq = equivalents(g, h);
    let chosen: Vec<[i32; 3]> =
        eq.iter().copied().filter(|&k| fullprof_box(k, rule) && fullprof_asu(k, rule) == Some(k)).collect();
    match chosen[..] {
        [one] => (one, true),
        _ => (*eq.iter().next_back().unwrap(), false),
    }
}

/// The unique, allowed reflections with dmin ≤ d ≤ dmax, sorted by decreasing d (increasing 2θ), then h, k, l.
/// With a wavelength, 2θ is given for each reflection.
pub fn generate(g: &SpaceGroup, cell: &Cell, dmin: f64, dmax: f64, lambda: Option<f64>) -> Vec<Reflection> {
    let rule = asu_rule(g);
    let hmax = [cell.a, cell.b, cell.c].map(|x| (x / dmin).floor() as i32);
    let mut seen = BTreeSet::new();
    let mut out = Vec::new();
    for h in -hmax[0]..=hmax[0] {
        for k in -hmax[1]..=hmax[1] {
            for l in -hmax[2]..=hmax[2] {
                let hkl = [h, k, l];
                if hkl == [0, 0, 0] {
                    continue;
                }
                let Some(d) = cell.d_spacing([h as f64, k as f64, l as f64]) else { continue };
                if d < dmin * (1.0 - 1e-12) || d > dmax * (1.0 + 1e-12) || is_absent(g, hkl) {
                    continue;
                }
                let (rep, fullprof_choice) = representative(g, hkl, rule);
                if seen.insert(rep) {
                    let two_theta = lambda.and_then(|w| cell.two_theta([rep[0] as f64, rep[1] as f64, rep[2] as f64], w));
                    out.push(Reflection { hkl: rep, multiplicity: equivalents(g, rep).len(), d, two_theta, fullprof_choice });
                }
            }
        }
    }
    out.sort_by(|x, y| y.d.total_cmp(&x.d).then(x.hkl.cmp(&y.hkl)));
    out
}

/// A systematically absent reflection (W8: drawn faintly on demand, "a peak on one says the space group or the cell is
/// wrong"; app request 20, L100).
#[derive(Clone, Debug, PartialEq)]
pub struct Forbidden {
    pub hkl: [i32; 3],
    pub multiplicity: usize,
    pub d: f64,
    pub two_theta: Option<f64>,
    /// The condition it breaks, as the tables write them: "h0l: l = 2n", "hkl: h + k = 2n", "00l: l = 4n".
    pub condition: String,
}

/// The unique, systematically absent reflections with dmin ≤ d ≤ dmax, sorted as `generate` sorts the allowed ones,
/// each with the condition that forbids it (`absence_condition`).
pub fn forbidden(g: &SpaceGroup, cell: &Cell, dmin: f64, dmax: f64, lambda: Option<f64>) -> Vec<Forbidden> {
    let rule = asu_rule(g);
    let hmax = [cell.a, cell.b, cell.c].map(|x| (x / dmin).floor() as i32);
    let mut seen = BTreeSet::new();
    let mut out = Vec::new();
    for h in -hmax[0]..=hmax[0] {
        for k in -hmax[1]..=hmax[1] {
            for l in -hmax[2]..=hmax[2] {
                let hkl = [h, k, l];
                if hkl == [0, 0, 0] || !is_absent(g, hkl) {
                    continue;
                }
                let Some(d) = cell.d_spacing([h as f64, k as f64, l as f64]) else { continue };
                if d < dmin * (1.0 - 1e-12) || d > dmax * (1.0 + 1e-12) {
                    continue;
                }
                let (rep, _) = representative(g, hkl, rule);
                if seen.insert(rep) {
                    let two_theta = lambda.and_then(|w| cell.two_theta([rep[0] as f64, rep[1] as f64, rep[2] as f64], w));
                    let condition = absence_condition(g, rep).unwrap_or_default();
                    out.push(Forbidden { hkl: rep, multiplicity: equivalents(g, rep).len(), d, two_theta, condition });
                }
            }
        }
    }
    out.sort_by(|x, y| y.d.total_cmp(&x.d).then(x.hkl.cmp(&y.hkl)));
    out
}

/// The reflection condition h breaks, or None if it is allowed: from an operator {R|t} with h·R = h and h·t not
/// integral, the class of reflections R leaves fixed (its zone: "hkl" for a centring translation, "h0l" for a glide
/// plane normal to b, "00l" for a screw axis along c, "hhl", …) and the condition h·t integral on that class, written
/// as "a·h + b·k + c·l = m·n". Operators that leave every reflection fixed (centrings) are tried first, then those
/// with a plane of fixed reflections, then a line.
pub fn absence_condition(g: &SpaceGroup, h: [i32; 3]) -> Option<String> {
    let identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    let mut best: Option<(usize, String)> = None;
    for op in &g.ops {
        if h_times_r(h, &op.r) != h || (0..3).map(|i| h[i] * op.t[i]).sum::<i32>().rem_euclid(DEN) == 0 {
            continue;
        }
        let basis = if op.r == identity { vec![[1, 0, 0], [0, 1, 0], [0, 0, 1]] } else { fixed_basis(&op.r) };
        if best.as_ref().is_some_and(|b| b.0 >= basis.len()) {
            continue;
        }
        const L: [&str; 3] = ["h", "k", "l"];
        // Each basis vector stands for the letter of its first non-zero index.
        let letter: Vec<&str> = basis.iter().map(|b| L[b.iter().position(|&v| v != 0).unwrap()]).collect();
        let term = |c: i32, x: &str| match c {
            1 => x.to_string(),
            -1 => format!("-{x}"),
            c => format!("{c}{x}"),
        };
        let comps: Vec<String> = (0..3)
            .map(|i| {
                let t: Vec<String> = basis.iter().zip(&letter).filter(|(b, _)| b[i] != 0).map(|(b, x)| term(b[i], x)).collect();
                if t.is_empty() { "0".to_string() } else { t.join("+").replace("+-", "-") }
            })
            .collect();
        let zone = if comps.iter().all(|c| c.len() <= 2) { comps.concat() } else { comps.join(",") };
        let a: Vec<i32> = basis.iter().map(|b| (0..3).map(|i| b[i] * op.t[i]).sum::<i32>().rem_euclid(DEN)).collect();
        let gcd = a.iter().fold(DEN, |x, &y| gcd(x, y));
        let m = DEN / gcd;
        let form: Vec<String> = a.iter().zip(&letter).filter(|(&c, _)| c != 0).map(|(&c, x)| term(c / gcd, x)).collect();
        let text = format!("{zone}: {} = {m}n", form.join(" + "));
        best = Some((basis.len(), text));
    }
    best.map(|b| b.1)
}

fn gcd(a: i32, b: i32) -> i32 {
    if b == 0 { a.abs() } else { gcd(b, a % b) }
}

/// A basis of the reflections a rotation leaves fixed (h·R = h), in echelon form with small integers.
fn fixed_basis(r: &[[i32; 3]; 3]) -> Vec<[i32; 3]> {
    let mut fixed: Vec<[i32; 3]> = Vec::new();
    for a in -2..=2 {
        for b in -2..=2 {
            for c in -2..=2 {
                let v = [a, b, c];
                if v != [0, 0, 0] && h_times_r(v, r) == v {
                    fixed.push(v);
                }
            }
        }
    }
    // Echelon: pivots first by position, then the smallest vectors, the first index positive.
    fixed.retain(|v| v.iter().find(|&&x| x != 0).is_some_and(|&x| x > 0));
    fixed.sort_by_key(|v| (v.iter().position(|&x| x != 0).unwrap(), v.iter().map(|x| x.abs()).sum::<i32>(), std::cmp::Reverse(*v)));
    let mut basis: Vec<[i32; 3]> = Vec::new();
    for v in fixed {
        let p = v.iter().position(|&x| x != 0).unwrap();
        if basis.iter().any(|b| b.iter().position(|&x| x != 0).unwrap() == p) {
            continue;
        }
        // Independent of the basis so far, and zero at the earlier pivots.
        if basis.iter().all(|b| v[b.iter().position(|&x| x != 0).unwrap()] == 0) {
            basis.push(v);
        }
    }
    basis
}

#[cfg(test)]
mod tests {
    #[test]
    fn absence_conditions_as_the_tables_write_them() {
        let cond = |sym: &str, h: [i32; 3]| absence_condition(&SpaceGroup::from_symbol(sym).unwrap().0, h).unwrap_or_default();
        assert_eq!(cond("P 21/c", [0, 1, 0]), "0k0: k = 2n");
        assert_eq!(cond("P 21/c", [1, 0, 1]), "h0l: l = 2n");
        assert_eq!(cond("I m -3 m", [1, 0, 0]), "hkl: h + k + l = 2n");
        assert_eq!(cond("P 41 21 2", [0, 0, 2]), "00l: l = 4n");
        assert_eq!(cond("P n m a", [0, 1, 0]), "0kl: k + l = 2n");
        assert_eq!(absence_condition(&SpaceGroup::from_symbol("P 21/c").unwrap().0, [1, 0, 2]), None);
        let (g, _) = SpaceGroup::from_symbol("F d -3 m").unwrap();
        let c = Cell::new(8.0, 8.0, 8.0, 90.0, 90.0, 90.0).unwrap();
        let f = forbidden(&g, &c, 2.0, 10.0, None);
        // The d glide of origin choice 2 has translations off the tables' origin: (2 0 0) breaks "hk0: h + 3k = 4n",
        // the same lattice of reflections as the tables' "hk0: h + k = 4n" on (h k 0) with h, k even.
        assert!(f.iter().any(|r| r.hkl == [2, 0, 0] && r.condition.ends_with("= 4n")), "{f:?}");
        assert!(f.iter().all(|r| is_absent(&g, r.hkl)) && !f.iter().any(|r| r.hkl == [2, 2, 0]));
    }

    use super::*;

    #[test]
    fn centring_and_glides_from_operators() {
        let (g, _) = SpaceGroup::from_symbol("F m -3 m").unwrap();
        assert!(is_absent(&g, [1, 0, 0]) && is_absent(&g, [1, 1, 0]) && !is_absent(&g, [1, 1, 1]) && !is_absent(&g, [2, 0, 0]));
        let (g, _) = SpaceGroup::from_symbol("P 21/c").unwrap();
        assert!(is_absent(&g, [0, 1, 0]) && !is_absent(&g, [0, 2, 0]));
        assert!(is_absent(&g, [1, 0, 1]) && !is_absent(&g, [1, 0, 2]), "c glide: h0l with l odd");
    }

    #[test]
    fn cubic_list() {
        let (g, _) = SpaceGroup::from_symbol("F m -3 m").unwrap();
        let cell = Cell::new(4.0, 4.0, 4.0, 90.0, 90.0, 90.0).unwrap();
        let r = generate(&g, &cell, 1.0, 10.0, None);
        let hkl: Vec<_> = r.iter().map(|x| (x.hkl, x.multiplicity)).collect();
        assert_eq!(hkl[..4], [([1, 1, 1], 8), ([2, 0, 0], 6), ([2, 2, 0], 12), ([3, 1, 1], 24)]);
        assert!(r.iter().all(|x| x.fullprof_choice));
    }
}
