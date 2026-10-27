//! Restraints (P40): distances and angles between atoms, and linear combinations of parameters, held near a value
//! within a standard deviation as extra observations (FullProf's soft constraints; manual pp. 123–124, and the
//! tutorial "How to use restraints and the rigid body description", `manuals/fullprof/docs`, pp. 2–4).
//!
//! The cost minimised is S = Σ_i w_i (y_i − y_c,i)² + χ²_ν Σ_j ((R_j − R_c,j)/σ_j)², the restraint term multiplied by
//! χ²_ν = S_data/(N − P) of the current cycle, as FullProf's (tutorial, p. 2): near the solution the restraints weigh
//! as observations of the same quality as the profile. The statistics (Rp, Rwp, χ²) are those of the data alone, and
//! P is not reduced by restraints (they are observations, not constraints).
//!
//! An atom of a restraint is taken through an operator {R|T}: R from FullProf's table of rotations, numbered ITnum
//! (1–24 in non-hexagonal frames, 25–36 in hexagonal ones; negative: combined with the inversion; manual p. 123), and
//! T the operator's whole translation. Distances and angles are computed in the cell's metric, generically over
//! `dual::Real`, so that their derivatives with respect to the atoms' coordinates and the cell are exact.

use crate::dual::{Dual, Real};
use crate::jacobian::Field;
use crate::powder::{Intensities, Model};
use crate::symop::SymOp;

/// FullProf's ITnum rotations (manual p. 123): non-hexagonal frames, 1–24.
const ITNUM_CUBIC: [&str; 24] = [
    "x,y,z", "-x,-y,z", "-x,y,-z", "x,-y,-z", "z,x,y", "z,-x,-y", "-z,-x,y", "-z,x,-y", "y,z,x", "-y,z,-x", "y,-z,-x", "-y,-z,x", "y,x,-z", "-y,-x,-z", "y,-x,z",
    "-y,x,z", "x,z,-y", "-x,z,y", "-x,-z,-y", "x,-z,y", "z,y,-x", "z,-y,x", "-z,y,x", "-z,-y,-x",
];
/// Hexagonal frames, 25–36.
const ITNUM_HEX: [&str; 12] = ["x,y,z", "-y,x-y,z", "-x+y,-x,z", "-x,-y,z", "y,-x+y,z", "x-y,x,z", "y,x,-z", "x-y,-y,-z", "-x,-x+y,-z", "-y,-x,-z", "-x+y,y,-z", "x,x-y,-z"];

/// An operator {R|T} applied to an atom of a restraint.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Operator {
    pub rot: [[i32; 3]; 3],
    pub t: [f64; 3],
}

impl Operator {
    /// FullProf's operator ITnum (negative: times the inversion) with translation T.
    pub fn itnum(n: i32, t: [f64; 3]) -> Option<Operator> {
        let k = n.unsigned_abs() as usize;
        let s = match k {
            1..=24 => ITNUM_CUBIC[k - 1],
            25..=36 => ITNUM_HEX[k - 25],
            _ => return None,
        };
        let op = SymOp::parse_xyz(s).ok()?;
        let sign = if n < 0 { -1 } else { 1 };
        Some(Operator { rot: op.r.map(|r| r.map(|v| sign * v)), t })
    }

    pub fn identity() -> Operator {
        Operator { rot: [[1, 0, 0], [0, 1, 0], [0, 0, 1]], t: [0.0; 3] }
    }

    fn apply<T: Real>(&self, x: &[T; 3]) -> [T; 3] {
        let mut y = [T::cst(0.0); 3];
        for (i, yi) in y.iter_mut().enumerate() {
            *yi = x[0] * self.rot[i][0] as f64 + x[1] * self.rot[i][1] as f64 + x[2] * self.rot[i][2] as f64 + self.t[i];
        }
        y
    }
}

/// What a restraint holds.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Kind {
    /// The distance (Å) from atom a to atom b taken through the operator.
    Distance { atoms: [usize; 2], op: Operator },
    /// The angle (degrees) at atom b taken through op1 between atom a and atom c taken through op2.
    Angle { atoms: [usize; 3], ops: [Operator; 2] },
}

/// A distance or angle restraint of a phase.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Restraint {
    pub phase: usize,
    pub kind: Kind,
    pub target: f64,
    pub sigma: f64,
}

/// The scalar product in the cell's metric of two fractional vectors.
fn dot<T: Real>(cell: &[T; 6], u: &[T; 3], v: &[T; 3]) -> T {
    let (a, b, c) = (cell[0], cell[1], cell[2]);
    let (ca, cb, cg) = (cell[3].to_radians().cos(), cell[4].to_radians().cos(), cell[5].to_radians().cos());
    u[0] * v[0] * a * a + u[1] * v[1] * b * b + u[2] * v[2] * c * c + (u[0] * v[1] + u[1] * v[0]) * a * b * cg + (u[0] * v[2] + u[2] * v[0]) * a * c * cb + (u[1] * v[2] + u[2] * v[1]) * b * c * ca
}

fn sub<T: Real>(u: &[T; 3], v: &[T; 3]) -> [T; 3] {
    [u[0] - v[0], u[1] - v[1], u[2] - v[2]]
}

impl Restraint {
    /// The atoms it involves.
    pub fn atoms(&self) -> Vec<usize> {
        match &self.kind {
            Kind::Distance { atoms, .. } => atoms.to_vec(),
            Kind::Angle { atoms, .. } => atoms.to_vec(),
        }
    }

    /// The value from a cell and the atoms' positions (in the order of `atoms()`).
    fn value_g<T: Real>(&self, cell: &[T; 6], x: &[[T; 3]]) -> T {
        match &self.kind {
            Kind::Distance { op, .. } => {
                let d = sub(&x[0], &op.apply(&x[1]));
                dot(cell, &d, &d).sqrt()
            }
            Kind::Angle { ops, .. } => {
                let centre = ops[0].apply(&x[1]);
                let (u, v) = (sub(&x[0], &centre), sub(&ops[1].apply(&x[2]), &centre));
                let c = dot(cell, &u, &v) / (dot(cell, &u, &u) * dot(cell, &v, &v)).sqrt();
                c.acos().to_degrees()
            }
        }
    }

    fn inputs(&self, m: &Model) -> ([f64; 6], Vec<[f64; 3]>) {
        let ph = &m.phases[self.phase];
        let c = &ph.cell;
        let Intensities::Structure(st, _) = &ph.intensities else { panic!("restraint on phase {} without atoms", self.phase + 1) };
        ([c.a, c.b, c.c, c.alpha, c.beta, c.gamma], self.atoms().iter().map(|&a| st.atoms[a].x).collect())
    }

    /// The calculated value (Å or degrees).
    pub fn value(&self, m: &Model) -> f64 {
        let (cell, x) = self.inputs(m);
        self.value_g(&cell, &x)
    }

    /// The calculated value and its derivatives with respect to the fields it depends on: the six cell parameters and
    /// the coordinates of its atoms.
    pub fn derivatives(&self, m: &Model) -> (f64, Vec<(Field, f64)>) {
        let (cell, x) = self.inputs(m);
        let p = self.phase;
        let atoms = self.atoms();
        // Up to 6 + 3·3 = 15 variables.
        let cd: [Dual<15>; 6] = std::array::from_fn(|k| Dual::var(cell[k], k));
        let xd: Vec<[Dual<15>; 3]> = x.iter().enumerate().map(|(j, xj)| std::array::from_fn(|k| Dual::var(xj[k], 6 + 3 * j + k))).collect();
        let v = self.value_g(&cd, &xd);
        let mut out: Vec<(Field, f64)> = (0..6).map(|k| (Field::Cell(p, k), v.d[k])).collect();
        for (j, &a) in atoms.iter().enumerate() {
            for k in 0..3 {
                out.push((Field::AtomX(p, a, k), v.d[6 + 3 * j + k]));
            }
        }
        // An atom used twice (C2–C2 through the inversion) contributes through both positions.
        out.sort_by_key(|a| a.0);
        out.dedup_by(|b, a| {
            if a.0 == b.0 {
                a.1 += b.1;
                true
            } else {
                false
            }
        });
        (v.v, out)
    }

    /// A label: "phase1 Ow3–Hw1", "phase1 Hw1–Ow3–Hw2".
    pub fn label(&self, m: &Model) -> String {
        let names: Vec<String> = match &m.phases[self.phase].intensities {
            Intensities::Structure(st, _) => self.atoms().iter().map(|&a| st.atoms[a].label.clone()).collect(),
            Intensities::List(_) => self.atoms().iter().map(|a| format!("atom{}", a + 1)).collect(),
        };
        format!("phase{} {}", self.phase + 1, names.join("–"))
    }
}

/// The distance and angle restraints of a `.pcr`'s phases, with the atoms found by their labels.
pub fn from_pcr(p: &crate::pcr::Pcr) -> Result<Vec<Restraint>, String> {
    let mut out = Vec::new();
    for (k, ph) in p.phases.iter().enumerate() {
        let find = |l: &str| ph.atoms.iter().position(|a| a.label == l).ok_or(format!("phase {}: restraint on {l}, not an atom of the phase", k + 1));
        let op = |n: i32, t: [f64; 3]| Operator::itnum(n, t).ok_or(format!("phase {}: ITnum {n} is not in FullProf's table (1–36)", k + 1));
        for d in &ph.distance_restraints {
            out.push(Restraint { phase: k, kind: Kind::Distance { atoms: [find(&d.atoms[0])?, find(&d.atoms[1])?], op: op(d.itnum, d.t)? }, target: d.value, sigma: d.sigma });
        }
        for a in &ph.angle_restraints {
            let atoms = [find(&a.atoms[0])?, find(&a.atoms[1])?, find(&a.atoms[2])?];
            out.push(Restraint { phase: k, kind: Kind::Angle { atoms, ops: [op(a.itnum[0], a.t[0])?, op(a.itnum[1], a.t[1])?] }, target: a.value, sigma: a.sigma });
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::Operator;

    #[test]
    fn itnum_table() {
        // ITnum 3: 2 (0, y, 0) → (−x, y, −z); negative: times the inversion.
        assert_eq!(Operator::itnum(3, [0.0; 3]).unwrap().rot, [[-1, 0, 0], [0, 1, 0], [0, 0, -1]]);
        assert_eq!(Operator::itnum(-3, [0.0; 3]).unwrap().rot, [[1, 0, 0], [0, -1, 0], [0, 0, 1]]);
        // ITnum 26: 3+ (0, 0, z) → (−y, x − y, z).
        assert_eq!(Operator::itnum(26, [0.0; 3]).unwrap().rot, [[0, -1, 0], [1, -1, 0], [0, 0, 1]]);
        assert!(Operator::itnum(37, [0.0; 3]).is_none());
    }
}
