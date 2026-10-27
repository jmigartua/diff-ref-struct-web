//! Sites in a space group: site-symmetry group, multiplicity, orbit, Wyckoff position, and the constraints that the
//! site symmetry puts on the coordinates and on the anisotropic displacement parameters (phase 1, item 3).
//!
//! The site-symmetry group of a position x is the set of operators {R|t} of the group, with the lattice translation
//! that makes it exact, such that R·x + t = x. Its order times the multiplicity of x is the order of the group
//! (centring included). The constraints follow from the same operators:
//! - coordinates: (R − I)·x = −t for every site operator, solved exactly in rationals;
//! - anisotropic displacements β (FullProf's betas, exp(−hᵀβh)): R·β·Rᵀ = β for every site rotation.
//!
//! Both are returned as linear relations, value_i = c_i + Σ_k a_ik·p_k over free parameters p_k, each free
//! parameter being one of the variables itself (x, y or z; β11 … β23), the first ones in the usual order. In every
//! conventional setting the site rotations mix only axes of equal length, so the U_ij of the CIF convention
//! (β_ij = 2π² a*_i a*_j U_ij) obey the same relations as the β_ij (checked by the test `spglib_wyckoff`).
//! Wyckoff letters and site-symmetry symbols are names from spglib's table (`wyckoff_table`).

use crate::spacegroup::{point_group, SpaceGroup};
use crate::symop::{SymOp, DEN};
use crate::wyckoff_table::{positions, WyckoffEntry};
use std::fmt;

/// Default tolerance, in fractions of a lattice vector, for a position to be fixed by an operator.
pub const TOL: f64 = 1e-3;

// ----------------------------------------------------------------------------------------------- rationals

/// An exact fraction, always reduced, denominator positive.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Rational {
    pub num: i64,
    pub den: i64,
}

fn gcd(a: i64, b: i64) -> i64 {
    if b == 0 {
        a.abs()
    } else {
        gcd(b, a % b)
    }
}

impl Rational {
    pub fn new(num: i64, den: i64) -> Rational {
        assert!(den != 0, "zero denominator");
        let g = gcd(num, den).max(1) * den.signum();
        Rational { num: num / g, den: den / g }
    }
    pub fn int(n: i64) -> Rational {
        Rational { num: n, den: 1 }
    }
    pub fn zero() -> Rational {
        Rational::int(0)
    }
    pub fn is_zero(&self) -> bool {
        self.num == 0
    }
    pub fn value(&self) -> f64 {
        self.num as f64 / self.den as f64
    }
}

impl std::ops::Add for Rational {
    type Output = Rational;
    fn add(self, o: Rational) -> Rational {
        Rational::new(self.num * o.den + o.num * self.den, self.den * o.den)
    }
}

impl std::ops::Sub for Rational {
    type Output = Rational;
    fn sub(self, o: Rational) -> Rational {
        Rational::new(self.num * o.den - o.num * self.den, self.den * o.den)
    }
}

impl std::ops::Mul for Rational {
    type Output = Rational;
    fn mul(self, o: Rational) -> Rational {
        Rational::new(self.num * o.num, self.den * o.den)
    }
}

impl std::ops::Div for Rational {
    type Output = Rational;
    fn div(self, o: Rational) -> Rational {
        Rational::new(self.num * o.den, self.den * o.num)
    }
}

impl fmt::Display for Rational {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.den == 1 {
            write!(f, "{}", self.num)
        } else {
            write!(f, "{}/{}", self.num, self.den)
        }
    }
}

// ----------------------------------------------------------------------------------------------- constraints

/// Linear constraints on a set of named variables: variable i = `constant[i]` + Σ_k `coef[i][k]`·(free variable k),
/// where the free variables are the variables `free[k]` themselves.
#[derive(Clone, Debug, PartialEq)]
pub struct Constraints {
    /// Names of the variables, e.g. ["x", "y", "z"] or ["β11", …, "β23"].
    pub names: Vec<&'static str>,
    /// Indices of the free variables.
    pub free: Vec<usize>,
    /// Constant part of each variable.
    pub constant: Vec<Rational>,
    /// coef[i][k]: coefficient of free variable k in variable i.
    pub coef: Vec<Vec<Rational>>,
}

impl Constraints {
    /// Solves the homogeneous-or-not linear system rows·v = rhs for the variables `names`, keeping the earliest
    /// variables free: the system is reduced to row-echelon form with the columns in reverse order, so that pivots
    /// (dependent variables) fall on the later variables.
    fn solve(names: Vec<&'static str>, rows: &[(Vec<Rational>, Rational)]) -> Result<Constraints, String> {
        let n = names.len();
        // Columns reversed: column j of the matrix is variable n-1-j.
        let mut m: Vec<Vec<Rational>> = rows
            .iter()
            .map(|(r, b)| {
                let mut v: Vec<Rational> = (0..n).map(|j| r[n - 1 - j]).collect();
                v.push(*b);
                v
            })
            .collect();
        let mut pivots: Vec<usize> = Vec::new();
        let mut row = 0;
        for col in 0..n {
            let Some(p) = (row..m.len()).find(|&i| !m[i][col].is_zero()) else { continue };
            m.swap(row, p);
            let lead = m[row][col];
            for v in m[row].iter_mut() {
                *v = *v / lead;
            }
            let prow = m[row].clone();
            for (i, r) in m.iter_mut().enumerate() {
                if i != row && !r[col].is_zero() {
                    let f = r[col];
                    for (v, p) in r.iter_mut().zip(&prow) {
                        *v = *v - *p * f;
                    }
                }
            }
            pivots.push(col);
            row += 1;
            if row == m.len() {
                break;
            }
        }
        if m[row..].iter().any(|r| !r[n].is_zero() && r[..n].iter().all(Rational::is_zero)) {
            return Err("inconsistent constraints".into());
        }
        let var = |col: usize| n - 1 - col;
        let free: Vec<usize> = (0..n).filter(|&v| !pivots.iter().any(|&c| var(c) == v)).collect();
        let mut constant = vec![Rational::zero(); n];
        let mut coef = vec![vec![Rational::zero(); free.len()]; n];
        for (k, &v) in free.iter().enumerate() {
            coef[v][k] = Rational::int(1);
        }
        for (r, &c) in pivots.iter().enumerate() {
            let v = var(c);
            constant[v] = m[r][n];
            for (k, &fv) in free.iter().enumerate() {
                let col = n - 1 - fv;
                coef[v][k] = Rational::zero() - m[r][col];
            }
        }
        Ok(Constraints { names, free, constant, coef })
    }

    /// The variable i as an expression, e.g. "1/4", "x", "2x", "-x+1/2", "β11/2", "0".
    pub fn expression(&self, i: usize) -> String {
        let mut s = String::new();
        for (k, c) in self.coef[i].iter().enumerate() {
            if c.is_zero() {
                continue;
            }
            let name = self.names[self.free[k]];
            let (neg, mag) = if c.num < 0 { (true, Rational::new(-c.num, c.den)) } else { (false, *c) };
            s.push_str(if neg { "-" } else if s.is_empty() { "" } else { "+" });
            if mag.num != 1 {
                s.push_str(&mag.num.to_string());
            }
            s.push_str(name);
            if mag.den != 1 {
                s.push_str(&format!("/{}", mag.den));
            }
        }
        let c = self.constant[i];
        if !c.is_zero() || s.is_empty() {
            if !s.is_empty() && c.num > 0 {
                s.push('+');
            }
            s.push_str(&c.to_string());
        }
        s
    }

    /// Whether variable i is fixed (no free part).
    pub fn is_fixed(&self, i: usize) -> bool {
        self.coef[i].iter().all(Rational::is_zero)
    }

    /// Values of the variables from values of the free variables.
    pub fn apply(&self, free_values: &[f64]) -> Vec<f64> {
        (0..self.names.len())
            .map(|i| self.constant[i].value() + self.coef[i].iter().zip(free_values).map(|(c, v)| c.value() * v).sum::<f64>())
            .collect()
    }
}

// ----------------------------------------------------------------------------------------------- the site

/// A site-symmetry operator: the rotation and the exact translation (twelfths, not reduced) that fix the site.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SiteOp {
    pub r: [[i32; 3]; 3],
    pub t: [i32; 3],
}

impl SiteOp {
    pub fn apply(&self, x: [f64; 3]) -> [f64; 3] {
        let mut y = [0.0; 3];
        for (i, yi) in y.iter_mut().enumerate() {
            *yi = (0..3).map(|j| self.r[i][j] as f64 * x[j]).sum::<f64>() + self.t[i] as f64 / DEN as f64;
        }
        y
    }
}

/// A position in a space group and what its symmetry implies.
#[derive(Clone, Debug)]
pub struct Site {
    /// The position on its site: the given one averaged over the site group, so that it is exactly fixed.
    pub position: [f64; 3],
    /// The site-symmetry operators, identity first.
    pub ops: Vec<SiteOp>,
    /// Number of equivalent positions in the conventional cell.
    pub multiplicity: usize,
    /// The equivalent positions in the cell, in [0, 1), the site's own first.
    pub orbit: Vec<[f64; 3]>,
    /// Point group of the site (crystal-class symbol), e.g. "m", "-1", "mm2", "-3m".
    pub point_group: &'static str,
    /// The Wyckoff position of the table, when the group is one of the 530 settings.
    pub wyckoff: Option<WyckoffEntry>,
    /// Constraints on x, y, z.
    pub coordinates: Constraints,
    /// Constraints on β11, β22, β33, β12, β13, β23 (and on the U_ij of the CIF convention).
    pub adp: Constraints,
}

fn frac(v: f64) -> f64 {
    v - v.floor()
}

fn same_mod1(a: [f64; 3], b: [f64; 3], tol: f64) -> bool {
    (0..3).all(|i| {
        let d = a[i] - b[i];
        (d - d.round()).abs() < tol
    })
}

/// The site-symmetry operators of x in g: the operators that map x onto itself up to a lattice translation, with
/// that translation absorbed so that they fix x exactly (within `tol`).
pub fn site_ops(g: &SpaceGroup, x: [f64; 3], tol: f64) -> Vec<SiteOp> {
    let mut out = Vec::new();
    for op in &g.ops {
        let y = SiteOp { r: op.r, t: op.t }.apply(x);
        let l: Vec<f64> = (0..3).map(|i| (y[i] - x[i]).round()).collect();
        if (0..3).all(|i| (y[i] - x[i] - l[i]).abs() < tol) {
            out.push(SiteOp { r: op.r, t: [0, 1, 2].map(|i| op.t[i] - DEN * l[i] as i32) });
        }
    }
    out
}

impl Site {
    /// Everything about position x in group g, with the default tolerance.
    pub fn new(g: &SpaceGroup, x: [f64; 3]) -> Result<Site, String> {
        Site::with_tolerance(g, x, TOL)
    }

    pub fn with_tolerance(g: &SpaceGroup, x: [f64; 3], tol: f64) -> Result<Site, String> {
        let ops = site_ops(g, x, tol);
        if !g.ops.len().is_multiple_of(ops.len()) {
            return Err(format!("{} site operators do not divide the {} of the group: tolerance too large?", ops.len(), g.ops.len()));
        }
        let mut position = [0.0; 3];
        for op in &ops {
            let y = op.apply(x);
            for i in 0..3 {
                position[i] += y[i] / ops.len() as f64;
            }
        }
        let multiplicity = g.ops.len() / ops.len();
        let mut orbit: Vec<[f64; 3]> = Vec::new();
        for op in &g.ops {
            let y = SiteOp { r: op.r, t: op.t }.apply(position).map(frac);
            if !orbit.iter().any(|o| same_mod1(*o, y, tol)) {
                orbit.push(y);
            }
        }
        if orbit.len() != multiplicity {
            return Err(format!("orbit of {} positions for multiplicity {multiplicity}: tolerance too large?", orbit.len()));
        }
        let symops: Vec<SymOp> = ops.iter().map(|o| SymOp::new(o.r, o.t)).collect();
        let (point_group, _, _) = point_group(&symops)?;
        let coordinates = coordinate_constraints(&ops)?;
        let adp = adp_constraints(&ops)?;
        let wyckoff = g.setting.and_then(|s| wyckoff_of(g, s.hall_number, position, multiplicity, tol));
        Ok(Site { position, ops, multiplicity, orbit, point_group, wyckoff, coordinates, adp })
    }
}

/// Constraints on x, y, z: (R − I)·x = −t for every site operator.
pub fn coordinate_constraints(ops: &[SiteOp]) -> Result<Constraints, String> {
    let mut rows = Vec::new();
    for op in ops {
        for i in 0..3 {
            let r: Vec<Rational> = (0..3).map(|j| Rational::int((op.r[i][j] - i32::from(i == j)) as i64)).collect();
            rows.push((r, Rational::new(-op.t[i] as i64, DEN as i64)));
        }
    }
    Constraints::solve(vec!["x", "y", "z"], &rows)
}

/// The six independent components of a symmetric tensor, in the order β11, β22, β33, β12, β13, β23.
const PAIRS: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];

/// Constraints on the anisotropic displacement parameters: R·β·Rᵀ = β for every site rotation.
pub fn adp_constraints(ops: &[SiteOp]) -> Result<Constraints, String> {
    let index = |k: usize, l: usize| PAIRS.iter().position(|&p| p == (k.min(l), k.max(l))).unwrap();
    let mut rows = Vec::new();
    for op in ops {
        for &(i, j) in &PAIRS {
            let mut r = vec![0i64; 6];
            for k in 0..3 {
                for l in 0..3 {
                    r[index(k, l)] += (op.r[i][k] * op.r[j][l]) as i64;
                }
            }
            r[index(i, j)] -= 1;
            rows.push((r.into_iter().map(Rational::int).collect(), Rational::zero()));
        }
    }
    Constraints::solve(vec!["β11", "β22", "β33", "β12", "β13", "β23"], &rows)
}

// ----------------------------------------------------------------------------------------------- Wyckoff letters

/// An affine coordinate triplet of the Wyckoff table, e.g. "x,2x,z" or "-x+1/2,y+1/4,1/8", evaluated at given x, y, z.
pub fn eval_triplet(s: &str, v: [f64; 3]) -> Result<[f64; 3], String> {
    let parts: Vec<&str> = s.split(',').collect();
    if parts.len() != 3 {
        return Err(format!("'{s}': three parts expected"));
    }
    let mut out = [0.0; 3];
    for (i, p) in parts.iter().enumerate() {
        let mut rest = p.trim();
        while !rest.is_empty() {
            let (sign, body) = match rest.as_bytes()[0] {
                b'+' => (1.0, &rest[1..]),
                b'-' => (-1.0, &rest[1..]),
                _ => (1.0, rest),
            };
            let end = body.find(['+', '-']).unwrap_or(body.len());
            let term = &body[..end];
            rest = &body[end..];
            let (k, var) = match term.chars().last() {
                Some(c @ ('x' | 'y' | 'z')) => (&term[..term.len() - 1], Some(c)),
                _ => (term, None),
            };
            let num = |t: &str| -> Result<f64, String> {
                match t.split_once('/') {
                    Some((a, b)) => Ok(a.parse::<f64>().map_err(|e| e.to_string())? / b.parse::<f64>().map_err(|e| e.to_string())?),
                    None if t.is_empty() => Ok(1.0),
                    None => t.parse::<f64>().map_err(|e| format!("'{s}': {e}")),
                }
            };
            let k = num(k)?;
            out[i] += sign
                * match var {
                    Some('x') => k * v[0],
                    Some('y') => k * v[1],
                    Some('z') => k * v[2],
                    _ => k,
                };
        }
    }
    Ok(out)
}

/// Generic values for the free parameters of a Wyckoff representative.
const GENERIC: [f64; 3] = [0.1123, 0.2357, 0.3571];

/// The Wyckoff position of x: the entry of the table with the same multiplicity whose representative's site group
/// fixes some image of x (an operator of the group and a lattice translation).
fn wyckoff_of(g: &SpaceGroup, hall_number: u16, x: [f64; 3], mult: usize, tol: f64) -> Option<WyckoffEntry> {
    for w in positions(hall_number) {
        if w.multiplicity as usize != mult {
            continue;
        }
        let Ok(rep) = eval_triplet(w.coord, GENERIC) else { continue };
        let sw = site_ops(g, rep, 1e-6);
        for op in &g.ops {
            let y = SiteOp { r: op.r, t: op.t }.apply(x);
            for l in 0..125 {
                let lv = [(l / 25) as f64 - 2.0, ((l / 5) % 5) as f64 - 2.0, (l % 5) as f64 - 2.0];
                let p = [y[0] + lv[0], y[1] + lv[1], y[2] + lv[2]];
                if sw.iter().all(|s| {
                    let q = s.apply(p);
                    (0..3).all(|i| (q[i] - p[i]).abs() < tol)
                }) {
                    return Some(*w);
                }
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn site(sym: &str, x: [f64; 3]) -> Site {
        let (g, _) = SpaceGroup::from_symbol(sym).unwrap();
        Site::new(&g, x).unwrap()
    }

    #[test]
    fn general_position_is_free() {
        let s = site("P n m a", [0.11, 0.17, 0.31]);
        assert_eq!(s.multiplicity, 8);
        assert_eq!(s.coordinates.free, vec![0, 1, 2]);
        assert_eq!(s.adp.free.len(), 6);
    }

    #[test]
    fn mirror_in_pnma_fixes_y() {
        let s = site("P n m a", [0.31, 0.25, 0.62]);
        assert_eq!((s.multiplicity, s.point_group), (4, "m"));
        assert_eq!((0..3).map(|i| s.coordinates.expression(i)).collect::<Vec<_>>(), ["x", "1/4", "z"]);
        // A mirror perpendicular to b: β12 = β23 = 0.
        assert_eq!((0..6).map(|i| s.adp.expression(i)).collect::<Vec<_>>(), ["β11", "β22", "β33", "0", "β13", "0"]);
    }

    #[test]
    fn hexagonal_and_cubic_relations() {
        // P 63/m m c, (x, 2x, z): β22 = 2β12-type relations.
        let s = site("P 63/m m c", [0.17, 0.34, 0.06]);
        assert_eq!((0..3).map(|i| s.coordinates.expression(i)).collect::<Vec<_>>(), ["x", "2x", "z"]);
        let s = site("F m -3 m", [0.0, 0.0, 0.0]);
        assert_eq!((0..3).map(|i| s.coordinates.expression(i)).collect::<Vec<_>>(), ["0", "0", "0"]);
        assert_eq!((0..6).map(|i| s.adp.expression(i)).collect::<Vec<_>>(), ["β11", "β11", "β11", "0", "0", "0"]);
        let s = site("F m -3 m", [0.21, 0.21, 0.21]);
        assert_eq!((0..3).map(|i| s.coordinates.expression(i)).collect::<Vec<_>>(), ["x", "x", "x"]);
        assert_eq!((0..6).map(|i| s.adp.expression(i)).collect::<Vec<_>>(), ["β11", "β11", "β11", "β12", "β12", "β12"]);
    }

    #[test]
    fn rational_arithmetic() {
        assert_eq!(Rational::new(6, -12), Rational::new(-1, 2));
        assert_eq!(Rational::new(1, 3) + Rational::new(1, 6), Rational::new(1, 2));
        assert_eq!(Rational::new(2, 3).to_string(), "2/3");
    }
}
