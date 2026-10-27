//! Rigid bodies (P41, Publish's P47): a group of atoms placed by its internal coordinates, an orientation and an
//! origin, as FullProf's STRMOD places them (V. Rodríguez's "General Rigid-Body-Constraints/TLS Subroutine", Jbt = 4;
//! FullProf's manual pp. 34–36, 100, 104; `docs/FullProf-Restraints-RigidBodies.pdf` pp. 5–15).
//!
//! - A group starts at an atom whose P6 (line 25-1, after Occ) is not 0: 1.0x a rigid group, 2.xx a satellite of the
//!   xx-th group (in the order of the file), 4 a Z-matrix group, 5.0x a group whose displacements are TLS (Schomaker &
//!   Trueblood 1968); P6 < 0 asks FullProf for its matrices in the `.out`. Its line 25-3 holds the first atom's internal
//!   coordinates, then X0 Y0 Z0 (the origin, fractional), CHI, P16 (a satellite's N1.N2), DEG (1: angles in degrees, 0:
//!   radians) and KIND; THETA and PHI are P7 and P8 of line 25-1. The group's other atoms follow, each with its three
//!   internal coordinates on line 25-3 (and, in a Z-matrix group, its connectivity N1 N2 N3 in place of P6–P8).
//! - KIND 0 spherical (r, θ, φ: x = r sin θ cos φ, y = r sin θ sin φ, z = r cos θ), 1 Cartesian (x, y, z, Å),
//!   2 cylindrical (ρ, φ, z), 3 Z-matrix (distance, bond angle, torsion to N1, N2, N3; RBdoc p. 6: the first atom at
//!   the origin, the second on +x, the third in the xy plane).
//! - The molecule's frame to the crystal's Cartesian frame (x ∥ a, y in the ab plane, z ∥ c*; manual p. 34):
//!   r_C = R(θ, φ, χ) r_m + M X0, R = R_z(φ) R_y(θ) R_z(χ) (RBdoc p. 8), M the cell's crystal-to-orthonormal matrix;
//!   then fractional, M⁻¹ r_C. FullProf's printed Euler matrix of citric acid (`citric-acid-sls_1.out`, lines 413–416)
//!   is this R to its five decimals.
//! - Found on FullProf's own outputs, where the papers say nothing: in a Z-matrix the third atom (and any atom whose N3
//!   is 0) lies on the +y side of the xy plane, and the torsion is the opposite of IUPAC's sense (with IUPAC's, or with
//!   the third atom on −y, citric acid's and PPH₃'s atoms are 0.5–1.2 Å from FullProf's).
//! - A satellite: its origin is the parent's atom N1 and its z axis the direction from the parent's atom N2 (N2 = 0:
//!   the parent's origin) to N1, FullProf writing the X0 and θ, φ this gives back into the file (`satellite_axis`);
//!   only its χ is its own. The satellite is placed from the values written, as the parent is.
//!
//! In the powder model (`Phase::rigid`, L347) the groups place their atoms, and their values are refined in place of
//! the atoms' coordinates (`Field::Rigid`): the derivatives through the atoms exact (dual numbers, `Group::derivative`),
//! the cell's too (`Group::cell_derivative`: a group keeps its shape in Å). A satellite's origin, θ and φ follow its
//! parent (`place`); that rule is checked only on the urea files, whose parents are not turned. The Jacobian takes the
//! satellites along (`placement_derivative`, L367: urea's columns against differences, test `tls`).
//!
//! TLS displacements (P6 = 5.0x, `Tls`, P41 and Publish's P47.2, L367): a group's atoms (and its satellites') take
//! anisotropic displacements U = T + A L Aᵀ + A S + Sᵀ Aᵀ (Schomaker & Trueblood 1968, as FullProf's manual p. 38
//! names them), A = [[0, z, −y], [−z, 0, x], [y, −x, 0]] from r = M (x − centre) in the crystal's Cartesian frame,
//! T in Å², L in rad², S in Å·rad as FullProf's `.pcr` holds them (S's rows 1, 2, 3, S21 S31 S32 its last three), and
//! FullProf's β = 2π² M⁻¹ U M⁻ᵀ. Against the β of FullProf's `.ort` files (`runs/rigid`): hmt's and pdd's five runs
//! within their print (10⁻⁵); urea's within it only with the centre at z = 0.33034, where FullProf prints its centre of
//! mass at 0.326 and the group's atoms' masses (CrysFML's `AtWe`) put it at 0.32957: frozen, for Josu (test
//! `tls`, L367). The centre is the model's (`Tls::centre`, from the group's origin); the `.pcr` import sets it at the
//! main group's centre of mass when P6 = 5.01 asks for it (`tls_of`), at the origin otherwise.
//!
//! Not calculated: P6 = 3 (internal coordinates generated from x, y, z) and the chains of P6 = 4.xx.

use crate::cell::Cell;
use crate::dual::{Dual, Real};
use crate::mat3::{self, Mat3};
use crate::pcr::{PcrAtom, PcrPhase};

/// How a group gives its atoms' internal coordinates (KIND).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Coordinates {
    Spherical,
    Cartesian,
    Cylindrical,
    ZMatrix,
}

/// An atom of a group: its index among the phase's atoms, its internal coordinates (in the group's unit) and, in a
/// Z-matrix, the atoms N1 N2 N3 (1-based within the group, 0 for none) they refer to.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Member {
    pub label: String,
    pub index: usize,
    pub internal: [f64; 3],
    pub connect: [usize; 3],
}

/// A satellite group's anchor (P16 = N1.N2): its parent group (0-based) and the parent's atoms N1 and N2 (1-based,
/// N2 = 0 the parent's origin).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Satellite {
    pub parent: usize,
    pub atom: usize,
    pub from: usize,
}

/// A rigid group as the `.pcr` gives it.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Group {
    /// P6 as written (its sign asks for FullProf's verbose output).
    pub option: f64,
    pub coordinates: Coordinates,
    /// DEG: the angles (θ, φ, χ and the internal ones) in degrees, else radians.
    pub degrees: bool,
    /// X0 Y0 Z0, fractional.
    pub origin: [f64; 3],
    /// THETA, PHI, CHI in the group's unit.
    pub theta: f64,
    pub phi: f64,
    pub chi: f64,
    pub satellite: Option<Satellite>,
    pub members: Vec<Member>,
    /// TLS displacements of the group and its satellites (P6 = 5.0x). None: the atoms' own B.
    #[cfg_attr(feature = "serde", serde(default))]
    pub displacements: Option<Tls>,
}

/// A group's TLS tensors (Schomaker & Trueblood 1968) and the centre they refer to.
#[derive(Clone, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Tls {
    /// T11 T22 T33 T12 T13 T23 (Å²).
    pub t: [f64; 6],
    /// L11 L22 L33 L12 L13 L23 (rad²).
    pub l: [f64; 6],
    /// S11 S22 S33 S12 S13 S23 S21 S31 S32 (Å·rad).
    pub s: [f64; 9],
    /// The centre the tensors refer to, less the group's origin (fractions): it moves with the origin.
    pub centre: [f64; 3],
}

/// The TLS values' names in FullProf's order (the `.pcr`'s and `Tls`'s).
pub const TLS_NAMES: [&str; 21] = ["t11", "t22", "t33", "t12", "t13", "t23", "l11", "l22", "l33", "l12", "l13", "l23", "s11", "s22", "s33", "s12", "s13", "s23", "s21", "s31", "s32"];

impl Tls {
    pub fn value(&self, i: usize) -> f64 {
        match i {
            0..=5 => self.t[i],
            6..=11 => self.l[i - 6],
            _ => self.s[i - 12],
        }
    }

    pub fn set_value(&mut self, i: usize, v: f64) {
        match i {
            0..=5 => self.t[i] = v,
            6..=11 => self.l[i - 6] = v,
            _ => self.s[i - 12] = v,
        }
    }
}

/// FullProf's β (11 22 33 12 13 23) of atoms at fractions x for the TLS values `v` (`TLS_NAMES`' order), the centre
/// (fractions) and the cell's a, b, c, α, β, γ, over any scalar: U = T + A L Aᵀ + A S + Sᵀ Aᵀ, β = 2π² M⁻¹ U M⁻ᵀ.
pub fn tls_betas<T: Real>(v: &[T; 21], centre: [T; 3], c: &[T; 6], x: &[[T; 3]]) -> Vec<[T; 6]> {
    let sym = |a: &[T]| [[a[0], a[3], a[4]], [a[3], a[1], a[5]], [a[4], a[5], a[2]]];
    let (tt, ll) = (sym(&v[0..6]), sym(&v[6..12]));
    let s = &v[12..21];
    let ss = [[s[0], s[3], s[4]], [s[6], s[1], s[5]], [s[7], s[8], s[2]]];
    let m = crystal_to_orthonormal(c);
    let mi = orthonormal_to_crystal(c);
    let mul = |a: &[[T; 3]; 3], b: &[[T; 3]; 3]| -> [[T; 3]; 3] { std::array::from_fn(|i| std::array::from_fn(|j| a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j])) };
    let tr = |a: &[[T; 3]; 3]| -> [[T; 3]; 3] { std::array::from_fn(|i| std::array::from_fn(|j| a[j][i])) };
    let two_pi2 = 2.0 * std::f64::consts::PI * std::f64::consts::PI;
    x.iter()
        .map(|p| {
            let r = mat_vec(&m, &sub(*p, centre));
            let z = T::cst(0.0);
            let a = [[z, r[2], -r[1]], [-r[2], z, r[0]], [r[1], -r[0], z]];
            let (at, st) = (tr(&a), tr(&ss));
            let (alat, as_, sat) = (mul(&mul(&a, &ll), &at), mul(&a, &ss), mul(&st, &at));
            let u: [[T; 3]; 3] = std::array::from_fn(|i| std::array::from_fn(|j| tt[i][j] + alat[i][j] + as_[i][j] + sat[i][j]));
            let b = mul(&mul(&mi, &u), &tr(&mi));
            [b[0][0] * two_pi2, b[1][1] * two_pi2, b[2][2] * two_pi2, b[0][1] * two_pi2, b[0][2] * two_pi2, b[1][2] * two_pi2]
        })
        .collect()
}

/// R = R_z(φ) R_y(θ) R_z(χ), the angles in radians.
pub fn rotation(theta: f64, phi: f64, chi: f64) -> Mat3 {
    rotation_of(theta, phi, chi)
}

fn sub<T: Real>(a: [T; 3], b: [T; 3]) -> [T; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn add<T: Real>(a: [T; 3], b: [T; 3]) -> [T; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale<T: Real>(a: [T; 3], s: T) -> [T; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn cross<T: Real>(a: [T; 3], b: [T; 3]) -> [T; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

fn unit<T: Real>(a: [T; 3]) -> [T; 3] {
    let n = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
    [a[0] / n, a[1] / n, a[2] / n]
}

fn mat_vec<T: Real>(m: &[[T; 3]; 3], v: &[T; 3]) -> [T; 3] {
    std::array::from_fn(|i| m[i][0] * v[0] + m[i][1] * v[1] + m[i][2] * v[2])
}

/// R_z(φ) R_y(θ) R_z(χ) over any scalar.
fn rotation_of<T: Real>(t: T, p: T, c: T) -> [[T; 3]; 3] {
    let (st, ct, sp, cp, sc, cc) = (t.sin(), t.cos(), p.sin(), p.cos(), c.sin(), c.cos());
    let z = T::cst(0.0);
    let rz = |s: T, c: T| [[c, -s, z], [s, c, z], [z, z, T::cst(1.0)]];
    let ry = [[ct, z, st], [z, T::cst(1.0), z], [-st, z, ct]];
    let mul = |a: &[[T; 3]; 3], b: &[[T; 3]; 3]| -> [[T; 3]; 3] { std::array::from_fn(|i| std::array::from_fn(|j| a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j])) };
    mul(&mul(&rz(sp, cp), &ry), &rz(sc, cc))
}

/// The number of a group's values (`Group::value`): the origin's three, θ, φ, χ, three per atom, then the 21 TLS
/// values of a TLS group (`TLS_NAMES`).
pub fn values_of(g: &Group) -> usize {
    6 + 3 * g.members.len() + if g.displacements.is_some() { 21 } else { 0 }
}

/// Which TLS value the group's k-th value is, if one.
pub fn tls_index(g: &Group, k: usize) -> Option<usize> {
    let first = 6 + 3 * g.members.len();
    (g.displacements.is_some() && k >= first && k < first + 21).then(|| k - first)
}

/// Whether the group's k-th value belongs to the group (origin, orientation, TLS) rather than to one of its atoms.
pub fn is_group_value(g: &Group, k: usize) -> bool {
    k < 6 || tls_index(g, k).is_some()
}

impl Group {
    fn radians<T: Real>(&self, v: T) -> T {
        if self.degrees {
            v * (std::f64::consts::PI / 180.0)
        } else {
            v
        }
    }

    /// Whether its displacements are TLS (P6 = ±5.0x).
    pub fn tls(&self) -> bool {
        self.option.abs().trunc() == 5.0
    }

    /// Its k-th value (`values_of`): 0–2 the origin X0 Y0 Z0 (fractional), 3 θ, 4 φ, 5 χ (the group's unit), then the
    /// internal coordinates of atom m at 6 + 3m + c.
    pub fn value(&self, k: usize) -> f64 {
        if let (Some(i), Some(t)) = (tls_index(self, k), &self.displacements) {
            return t.value(i);
        }
        match k {
            0..=2 => self.origin[k],
            3 => self.theta,
            4 => self.phi,
            5 => self.chi,
            _ => self.members[(k - 6) / 3].internal[(k - 6) % 3],
        }
    }

    /// Set its k-th value.
    pub fn set_value(&mut self, k: usize, v: f64) {
        if let Some(i) = tls_index(self, k) {
            if let Some(t) = &mut self.displacements {
                t.set_value(i, v);
            }
            return;
        }
        match k {
            0..=2 => self.origin[k] = v,
            3 => self.theta = v,
            4 => self.phi = v,
            5 => self.chi = v,
            _ => self.members[(k - 6) / 3].internal[(k - 6) % 3] = v,
        }
    }

    /// Whether the k-th value is an angle (in the group's unit; an internal coordinate by its kind).
    pub fn is_angle(&self, k: usize) -> bool {
        if tls_index(self, k).is_some() {
            return false;
        }
        match k {
            0..=2 => false,
            3..=5 => true,
            _ => match (self.coordinates, (k - 6) % 3) {
                (Coordinates::Cartesian, _) => false,
                (Coordinates::Cylindrical, c) => c == 1,
                (_, c) => c > 0,
            },
        }
    }

    /// The group's orientation R(θ, φ, χ).
    pub fn rotation(&self) -> Mat3 {
        rotation(self.radians(self.theta), self.radians(self.phi), self.radians(self.chi))
    }

    /// Its atoms in its own frame (Å), its values given as any scalar (`vals`, in `values_of`'s order).
    fn molecular_of<T: Real>(&self, vals: &[T]) -> Vec<[T; 3]> {
        let n = self.members.len();
        let internal = |m: usize| [vals[6 + 3 * m], vals[7 + 3 * m], vals[8 + 3 * m]];
        match self.coordinates {
            Coordinates::Spherical => (0..n)
                .map(|m| {
                    let [r, t, p] = internal(m);
                    let (t, p) = (self.radians(t), self.radians(p));
                    [r * t.sin() * p.cos(), r * t.sin() * p.sin(), r * t.cos()]
                })
                .collect(),
            Coordinates::Cartesian => (0..n).map(internal).collect(),
            Coordinates::Cylindrical => (0..n)
                .map(|m| {
                    let [rho, p, z] = internal(m);
                    let p = self.radians(p);
                    [rho * p.cos(), rho * p.sin(), z]
                })
                .collect(),
            Coordinates::ZMatrix => {
                let z = T::cst(0.0);
                let mut pos: Vec<[T; 3]> = Vec::with_capacity(n);
                for (k, m) in self.members.iter().enumerate() {
                    let [d, ang, tor] = internal(k);
                    let (t, w) = (self.radians(ang), self.radians(tor));
                    let at = |i: usize| pos.get(i.wrapping_sub(1)).copied().unwrap_or([z, z, z]);
                    let p = match k {
                        0 => [z, z, z],
                        1 => [d, z, z],
                        _ if k == 2 || m.connect[2] == 0 => {
                            // In the xy plane, on its +y side.
                            let (a, b) = (at(m.connect[0]), at(m.connect[1]));
                            let u = unit(sub(b, a));
                            let mut v = [-u[1], u[0], z];
                            if v[1].val() < 0.0 {
                                v = [-v[0], -v[1], z];
                            }
                            add(a, scale(add(scale(u, t.cos()), scale(v, t.sin())), d))
                        }
                        _ => {
                            // From N1 (bonded), N2 (the angle) and N3 (the torsion, FullProf's sense).
                            let (a, b, c) = (at(m.connect[0]), at(m.connect[1]), at(m.connect[2]));
                            let bc = unit(sub(a, b));
                            let nn = unit(cross(sub(b, c), bc));
                            let mm = cross(nn, bc);
                            let (l0, l1, l2) = (-(d * t.cos()), d * t.sin() * w.cos(), -(d * t.sin() * w.sin()));
                            add(a, add(add(scale(bc, l0), scale(mm, l1)), scale(nn, l2)))
                        }
                    };
                    pos.push(p);
                }
                pos
            }
        }
    }

    /// Its atoms' fractional coordinates from its values and the cell's a, b, c, α, β, γ (Å, degrees), given as any
    /// scalar: x = X0 + M⁻¹ R r_m, M the crystal-to-orthonormal matrix (`Cell::crystal_to_orthonormal`).
    fn fractional_of<T: Real>(&self, vals: &[T], c: &[T; 6]) -> Vec<[T; 3]> {
        let r = rotation_of(self.radians(vals[3]), self.radians(vals[4]), self.radians(vals[5]));
        let mi = orthonormal_to_crystal(c);
        self.molecular_of(vals).iter().map(|x| add([vals[0], vals[1], vals[2]], mat_vec(&mi, &mat_vec(&r, x)))).collect()
    }

    fn values(&self) -> Vec<f64> {
        (0..values_of(self)).map(|k| self.value(k)).collect()
    }

    /// Its atoms in its own frame (Å).
    pub fn molecular(&self) -> Vec<[f64; 3]> {
        self.molecular_of(&self.values())
    }

    /// Its atoms in the crystal's Cartesian frame (Å).
    pub fn cartesian(&self, cell: &Cell) -> Vec<[f64; 3]> {
        let m = cell.crystal_to_orthonormal();
        self.fractional(cell).iter().map(|x| mat3::mul_vec(&m, x)).collect()
    }

    /// Its atoms' fractional coordinates.
    pub fn fractional(&self, cell: &Cell) -> Vec<[f64; 3]> {
        self.fractional_of(&self.values(), &cell_values(cell))
    }

    /// The derivatives of its atoms' fractional coordinates with respect to its k-th value (per unit of the value as
    /// written): (member, coordinate, ∂x/∂value), the zeros left out. Exact (dual numbers).
    pub fn derivative(&self, k: usize, cell: &Cell) -> Vec<(usize, usize, f64)> {
        let vals: Vec<Dual<1>> = self.values().iter().enumerate().map(|(i, &v)| if i == k { Dual::var(v, 0) } else { Dual::cst(v) }).collect();
        let c = cell_values(cell).map(Dual::cst);
        nonzero(&self.fractional_of(&vals, &c))
    }

    /// The derivatives of its atoms' fractional coordinates with respect to the cell's k-th parameter (a, b, c in Å,
    /// α, β, γ in degrees): the group keeps its shape in Å, so its atoms' fractions move with the cell.
    pub fn cell_derivative(&self, k: usize, cell: &Cell) -> Vec<(usize, usize, f64)> {
        let vals: Vec<Dual<1>> = self.values().iter().map(|&v| Dual::cst(v)).collect();
        let c: [Dual<1>; 6] = std::array::from_fn(|i| if i == k { Dual::var(cell_values(cell)[i], 0) } else { Dual::cst(cell_values(cell)[i]) });
        nonzero(&self.fractional_of(&vals, &c))
    }
}

/// The derivatives of every group's atoms' fractional coordinates as `place` sets them, each satellite following its
/// parent (its origin, θ and φ from the parent's atoms, `satellite_axis`), as (atom index, coordinate, ∂x), the zeros
/// left out: with respect to group g's value k (`value` Some((g, k)), per unit as written) or the cell's k-th parameter
/// (`cell_k`). Exact (dual numbers).
pub fn placement_derivative(groups: &[Group], cell: &Cell, value: Option<(usize, usize)>, cell_k: Option<usize>) -> Vec<(usize, usize, f64)> {
    let cv = cell_values(cell);
    let c: [Dual<1>; 6] = std::array::from_fn(|i| if cell_k == Some(i) { Dual::var(cv[i], 0) } else { Dual::cst(cv[i]) });
    let (m, mi) = (crystal_to_orthonormal(&c), orthonormal_to_crystal(&c));
    let mut placed: Vec<(Vec<[Dual<1>; 3]>, [Dual<1>; 3])> = Vec::with_capacity(groups.len());
    let mut out = Vec::new();
    for (h, gr) in groups.iter().enumerate() {
        let mut vals: Vec<Dual<1>> = gr.values().iter().enumerate().map(|(i, &v)| if value == Some((h, i)) { Dual::var(v, 0) } else { Dual::cst(v) }).collect();
        if let Some(s) = gr.satellite.filter(|s| s.parent < h) {
            let (at, origin) = &placed[s.parent];
            let n1 = s.atom.checked_sub(1).and_then(|i| at.get(i)).map(|x| mat_vec(&m, x));
            let n2 = if s.from == 0 { Some(mat_vec(&m, origin)) } else { at.get(s.from - 1).map(|x| mat_vec(&m, x)) };
            if let (Some(n1), Some(n2)) = (n1, n2) {
                let v = unit(sub(n1, n2));
                let z = if v[2].val().abs() > 1.0 { Dual::cst(v[2].val().signum()) } else { v[2] };
                let (t, p) = (z.acos(), atan2(v[1], v[0]));
                let unit_of = |x: Dual<1>| if gr.degrees { x.to_degrees() } else { x };
                let o = mat_vec(&mi, &n1);
                vals[..5].copy_from_slice(&[o[0], o[1], o[2], unit_of(t), unit_of(p)]);
            }
        }
        let x = gr.fractional_of(&vals, &c);
        out.extend(nonzero(&x).into_iter().map(|(k, i, d)| (gr.members[k].index, i, d)));
        placed.push((x, [vals[0], vals[1], vals[2]]));
    }
    out
}

/// atan2(y, x) over any scalar: atan(y/x) or π/2 − atan(x/y) locally, whichever is better conditioned.
fn atan2<T: Real>(y: T, x: T) -> T {
    let f = y.val().atan2(x.val());
    if x.val() == 0.0 && y.val() == 0.0 {
        T::cst(f)
    } else if x.val().abs() >= y.val().abs() {
        let t = y / x;
        t.map(f, 1.0 / (1.0 + t.val() * t.val()))
    } else {
        let t = x / y;
        t.map(f, -1.0 / (1.0 + t.val() * t.val()))
    }
}

fn nonzero(x: &[[Dual<1>; 3]]) -> Vec<(usize, usize, f64)> {
    let mut out = Vec::new();
    for (m, p) in x.iter().enumerate() {
        for (c, d) in p.iter().enumerate() {
            if d.d[0] != 0.0 {
                out.push((m, c, d.d[0]));
            }
        }
    }
    out
}

fn cell_values(cell: &Cell) -> [f64; 6] {
    [cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma]
}

/// M⁻¹, the orthonormal-to-crystal matrix of a, b, c, α, β, γ (degrees), over any scalar: the inverse of the upper
/// triangular M = [[a, b cos γ, c cos β], [0, b sin γ, c (cos α − cos β cos γ)/sin γ], [0, 0, V/(a b sin γ)]].
fn orthonormal_to_crystal<T: Real>(c: &[T; 6]) -> [[T; 3]; 3] {
    let r = std::f64::consts::PI / 180.0;
    let (ca, cb, cg, sg) = ((c[3] * r).cos(), (c[4] * r).cos(), (c[5] * r).cos(), (c[5] * r).sin());
    let v = c[0] * c[1] * c[2] * (-(ca * ca) - cb * cb - cg * cg + ca * cb * cg * 2.0 + 1.0).sqrt();
    let (m00, m01, m02) = (c[0], c[1] * cg, c[2] * cb);
    let (m11, m12) = (c[1] * sg, c[2] * (ca - cb * cg) / sg);
    let m22 = v / (c[0] * c[1] * sg);
    let z = T::cst(0.0);
    // The inverse of an upper triangular matrix.
    let (i00, i11, i22) = (T::cst(1.0) / m00, T::cst(1.0) / m11, T::cst(1.0) / m22);
    let i01 = -(m01 * i00 * i11);
    let i12 = -(m12 * i11 * i22);
    let i02 = (m01 * m12 - m02 * m11) * i00 * i11 * i22;
    [[i00, i01, i02], [z, i11, i12], [z, z, i22]]
}

/// M, the crystal-to-orthonormal matrix of a, b, c, α, β, γ (degrees), over any scalar (`orthonormal_to_crystal`'s
/// inverse): x ∥ a, y in the ab plane, z ∥ c*.
fn crystal_to_orthonormal<T: Real>(c: &[T; 6]) -> [[T; 3]; 3] {
    let r = std::f64::consts::PI / 180.0;
    let (ca, cb, cg, sg) = ((c[3] * r).cos(), (c[4] * r).cos(), (c[5] * r).cos(), (c[5] * r).sin());
    let v = c[0] * c[1] * c[2] * (-(ca * ca) - cb * cb - cg * cg + ca * cb * cg * 2.0 + 1.0).sqrt();
    let z = T::cst(0.0);
    [[c[0], c[1] * cg, c[2] * cb], [z, c[1] * sg, c[2] * (ca - cb * cg) / sg], [z, z, v / (c[0] * c[1] * sg)]]
}

/// The atoms whose displacements group g's TLS gives: its members and its satellites' (by their group's index), as
/// (group, member).
pub fn tls_members(groups: &[Group], g: usize) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    if groups.get(g).and_then(|gr| gr.displacements.as_ref()).is_none() {
        return out;
    }
    for (h, gr) in groups.iter().enumerate() {
        let mut top = h;
        while let Some(s) = groups[top].satellite {
            if s.parent >= top {
                break;
            }
            top = s.parent;
        }
        if top == g {
            out.extend((0..gr.members.len()).map(|m| (h, m)));
        }
    }
    out
}

/// The TLS centre of group g (fractions): its origin plus `Tls::centre`.
pub fn tls_centre(g: &Group) -> Option<[f64; 3]> {
    g.displacements.as_ref().map(|t| add(g.origin, t.centre))
}

/// Every atom's TLS β, as (atom index, β), from the groups placed at x (the phase's atoms).
pub fn betas(groups: &[Group], cell: &Cell, x: &[[f64; 3]]) -> Vec<(usize, [f64; 6])> {
    let c = cell_values(cell);
    let mut out = Vec::new();
    for (g, gr) in groups.iter().enumerate() {
        let (Some(t), Some(centre)) = (&gr.displacements, tls_centre(gr)) else { continue };
        let v: [f64; 21] = std::array::from_fn(|i| t.value(i));
        let atoms: Vec<usize> = tls_members(groups, g).into_iter().map(|(h, m)| groups[h].members[m].index).collect();
        let pos: Vec<[f64; 3]> = atoms.iter().map(|&a| x[a]).collect();
        out.extend(atoms.into_iter().zip(tls_betas(&v, centre, &c, &pos)));
    }
    out
}

/// The derivatives of the TLS atoms' β, as (atom, β component 0–5, ∂β): `Wrt::Tls(g, i)` (linear), `Wrt::Position(a, c)`
/// (atom a's fraction c, its centre fixed), `Wrt::Centre(g, c)` (the centre's fraction c, the atoms fixed) and
/// `Wrt::Cell(k)` (the atoms' fractions fixed). Exact (dual numbers).
pub fn beta_derivative(groups: &[Group], cell: &Cell, x: &[[f64; 3]], wrt: Wrt) -> Vec<(usize, usize, f64)> {
    let cv = cell_values(cell);
    let mut out = Vec::new();
    for (g, gr) in groups.iter().enumerate() {
        let (Some(t), Some(centre)) = (&gr.displacements, tls_centre(gr)) else { continue };
        let atoms: Vec<usize> = tls_members(groups, g).into_iter().map(|(h, m)| groups[h].members[m].index).collect();
        let d = |v: f64, on: bool| if on { Dual::<1>::var(v, 0) } else { Dual::cst(v) };
        let v: [Dual<1>; 21] = std::array::from_fn(|i| d(t.value(i), wrt == Wrt::Tls(g, i)));
        let ce: [Dual<1>; 3] = std::array::from_fn(|i| d(centre[i], wrt == Wrt::Centre(g, i)));
        let c: [Dual<1>; 6] = std::array::from_fn(|i| d(cv[i], wrt == Wrt::Cell(i)));
        let pos: Vec<[Dual<1>; 3]> = atoms.iter().map(|&a| std::array::from_fn(|i| d(x[a][i], wrt == Wrt::Position(a, i)))).collect();
        for (&a, b) in atoms.iter().zip(tls_betas(&v, ce, &c, &pos)) {
            for (k, bk) in b.iter().enumerate() {
                if bk.d[0] != 0.0 {
                    out.push((a, k, bk.d[0]));
                }
            }
        }
    }
    out
}

/// What a TLS β derivative is taken with respect to (`beta_derivative`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Wrt {
    Tls(usize, usize),
    Position(usize, usize),
    Centre(usize, usize),
    Cell(usize),
}

/// A phase's TLS values from its `.pcr` (P6 = 5.0x with Furth ≠ 0): the further lines "name value code" after the
/// phase's line 30 ("! TLS parameters:", `PcrPhase::extra`), t11 … s32 by name (FullProf writes 12, T and L, or 21),
/// and their codewords in `TLS_NAMES`' order. None without them.
pub fn tls_of(ph: &PcrPhase) -> Option<([f64; 21], [f64; 21])> {
    let (mut v, mut c, mut any) = ([0.0; 21], [0.0; 21], false);
    for l in &ph.extra {
        let t: Vec<&str> = l.split_whitespace().collect();
        if t.len() >= 2 {
            if let Some(i) = TLS_NAMES.iter().position(|n| n.eq_ignore_ascii_case(t[0])) {
                v[i] = t[1].parse().ok()?;
                c[i] = t.get(2).and_then(|x| x.parse().ok()).unwrap_or(0.0);
                any = true;
            }
        }
    }
    any.then_some((v, c))
}

/// A `.pcr` phase's TLS (P6 = 5.0x with Furth ≠ 0; P41, L367) onto its groups placed at x: the further lines' T, L, S
/// to the first TLS group, its centre the main group's centre of mass when P6 = ±5.01 asks for it (FullProf's manual
/// p. 38: "x = 1"), its origin otherwise. Nothing when Furth is 0.
pub fn read_tls(ph: &PcrPhase, groups: &mut [Group], species: &[String], x: &[[f64; 3]], cell: &Cell) -> Result<(), String> {
    if ph.flag("Furth") == 0 {
        return Ok(());
    }
    let (v, _) = tls_of(ph).ok_or(format!("Furth = {} without the TLS lines t11 … s32", ph.flag("Furth")))?;
    let g = groups.iter().position(|g| g.tls()).ok_or(format!("TLS values (Furth = {}) without a TLS group (P6 = 5.0x)", ph.flag("Furth")))?;
    let centre = if (groups[g].option.abs().fract() * 100.0).round() == 1.0 {
        let c = centre_of_mass(&groups[g], species, x, cell).ok_or(format!("the centre of mass of group {}: an atom without an atomic weight", g + 1))?;
        sub(c, groups[g].origin)
    } else {
        [0.0; 3]
    };
    groups[g].displacements = Some(Tls { t: std::array::from_fn(|i| v[i]), l: std::array::from_fn(|i| v[6 + i]), s: std::array::from_fn(|i| v[12 + i]), centre });
    Ok(())
}

/// The centre of mass (fractions) of a group's own atoms, by CrysFML's atomic weights (`scattering_table::ATOMIC_WEIGHT`,
/// the species' letters), placed at x; None when a species has no weight.
pub fn centre_of_mass(g: &Group, species: &[String], x: &[[f64; 3]], cell: &Cell) -> Option<[f64; 3]> {
    let m = cell.crystal_to_orthonormal();
    let (mut sum, mut w) = ([0.0; 3], 0.0);
    for mm in &g.members {
        let sym: String = species.get(mm.index)?.chars().take_while(|c| c.is_ascii_alphabetic()).collect::<String>().to_uppercase();
        let wt = crate::scattering_table::ATOMIC_WEIGHT.iter().find(|e| e.0 == sym)?.1;
        let r = mat3::mul_vec(&m, &x[mm.index]);
        for i in 0..3 {
            sum[i] += wt * r[i];
        }
        w += wt;
    }
    (w > 0.0).then(|| mat3::mul_vec(&cell.orthonormal_to_crystal(), &[sum[0] / w, sum[1] / w, sum[2] / w]))
}

impl Coordinates {
    /// The control file's word for it.
    pub const NAMES: [(Coordinates, &'static str); 4] = [(Coordinates::Spherical, "spherical"), (Coordinates::Cartesian, "cartesian"), (Coordinates::Cylindrical, "cylindrical"), (Coordinates::ZMatrix, "z-matrix")];

    pub fn name(self) -> &'static str {
        Coordinates::NAMES.iter().find(|n| n.0 == self).map_or("spherical", |n| n.1)
    }

    /// Which of an atom's three internal coordinates are lengths (Å) and which angles (degrees), in order: what the
    /// control file writes as `rigid_length` and `rigid_angle`.
    pub fn split(self) -> (&'static [usize], &'static [usize]) {
        match self {
            Coordinates::Spherical | Coordinates::ZMatrix => (&[0], &[1, 2]),
            Coordinates::Cartesian => (&[0, 1, 2], &[]),
            Coordinates::Cylindrical => (&[0, 2], &[1]),
        }
    }
}

/// The parameter dictionary's id of a group's k-th value: `rigid_origin` (fractions), `rigid_orientation` (θ, φ, χ),
/// `rigid_length` (Å) or `rigid_angle` (degrees) for an internal coordinate.
pub fn value_id(g: &Group, k: usize) -> &'static str {
    if let Some(i) = tls_index(g, k) {
        return ["rigid_tls_t", "rigid_tls_l", "rigid_tls_s"][(i / 6).min(2)];
    }
    match k {
        0..=2 => "rigid_origin",
        3..=5 => "rigid_orientation",
        _ if g.is_angle(k) => "rigid_angle",
        _ => "rigid_length",
    }
}

/// The component a group's k-th value is within its id: x0, y0, z0; theta, phi, chi; an internal coordinate by the
/// group's kind (r, theta, phi; x, y, z; rho, phi, z; distance, angle, torsion).
pub fn value_component(g: &Group, k: usize) -> &'static str {
    if let Some(i) = tls_index(g, k) {
        return ["11", "22", "33", "12", "13", "23", "11", "22", "33", "12", "13", "23", "11", "22", "33", "12", "13", "23", "21", "31", "32"][i];
    }
    match k {
        0..=2 => ["x0", "y0", "z0"][k],
        3..=5 => ["theta", "phi", "chi"][k - 3],
        _ => {
            let c = (k - 6) % 3;
            match g.coordinates {
                Coordinates::Spherical => ["r", "theta", "phi"][c],
                Coordinates::Cartesian => ["x", "y", "z"][c],
                Coordinates::Cylindrical => ["rho", "phi", "z"][c],
                Coordinates::ZMatrix => ["distance", "angle", "torsion"][c],
            }
        }
    }
}

impl Group {
    /// The group with its angles in degrees (DEG = 1), its atoms where they were.
    pub fn in_degrees(mut self) -> Group {
        if !self.degrees {
            for k in 3..values_of(&self) {
                if self.is_angle(k) {
                    let v = self.value(k).to_degrees();
                    self.set_value(k, v);
                }
            }
            self.degrees = true;
        }
        self
    }
}

/// A satellite's origin (fractional) and the polar angles θ, φ of its z axis (in its unit), from its parent: the
/// parent's atom N1, and the direction from N2 (or the parent's origin) to N1. None for a group that is not a
/// satellite.
pub fn satellite_axis(groups: &[Group], g: usize, cell: &Cell) -> Option<([f64; 3], f64, f64)> {
    let s = groups[g].satellite?;
    let parent = groups.get(s.parent)?;
    let at = parent.cartesian(cell);
    let n1 = *at.get(s.atom.checked_sub(1)?)?;
    let n2 = if s.from == 0 { mat3::mul_vec(&cell.crystal_to_orthonormal(), &parent.origin) } else { *at.get(s.from - 1)? };
    let v = unit(sub(n1, n2));
    let (t, p) = (v[2].clamp(-1.0, 1.0).acos(), v[1].atan2(v[0]));
    let unit_of = |x: f64| if groups[g].degrees { x.to_degrees() } else { x };
    Some((mat3::mul_vec(&cell.orthonormal_to_crystal(), &n1), unit_of(t), unit_of(p)))
}

/// Place every group: each satellite first anchored to its parent (`satellite_axis`: its origin and its z axis; its χ
/// its own), then its atoms' fractional coordinates written into `x` (indexed by the phase's atoms).
pub fn place(groups: &mut [Group], cell: &Cell, x: &mut [[f64; 3]]) {
    for g in 0..groups.len() {
        if let Some((o, t, p)) = satellite_axis(groups, g, cell) {
            (groups[g].origin, groups[g].theta, groups[g].phi) = (o, t, p);
        }
        for (m, p) in groups[g].members.iter().zip(groups[g].fractional(cell)) {
            x[m.index] = p;
        }
    }
}

fn is_head(a: &PcrAtom) -> bool {
    a.rigid.1.len() >= 10 && a.rigid.0.first().is_some_and(|&p6| p6 != 0.0)
}

/// The rigid groups of a phase (|Jbt| = 4): each from its first atom to the next group's or to an atom whose label
/// ends without a number (a free atom, manual p. 35). Refused, the atom named, where a group's line cannot be read.
pub fn groups(ph: &PcrPhase) -> Result<Vec<Group>, String> {
    let mut out: Vec<Group> = Vec::new();
    let mut current: Option<Group> = None;
    for (i, a) in ph.atoms.iter().enumerate() {
        let free = !a.label.chars().last().is_some_and(|c| c.is_ascii_digit());
        if is_head(a) || free {
            out.extend(current.take());
        }
        if free && !is_head(a) {
            continue;
        }
        if is_head(a) {
            let (l1, l3) = (&a.rigid.0, &a.rigid.1);
            let option = l1[0];
            let kind = l3[9] as i32;
            let coordinates = match kind {
                0 => Coordinates::Spherical,
                1 => Coordinates::Cartesian,
                2 => Coordinates::Cylindrical,
                3 => Coordinates::ZMatrix,
                k => return Err(format!("atom {}: KIND = {k} (0 spherical, 1 Cartesian, 2 cylindrical, 3 Z-matrix)", a.label)),
            };
            let k = option.abs().trunc() as i32;
            if !matches!(k, 1 | 2 | 4 | 5) {
                return Err(format!("atom {}: P6 = {option} is not calculated (1 rigid, 2 satellite, 4 Z-matrix, 5 TLS)", a.label));
            }
            if k == 4 && coordinates != Coordinates::ZMatrix {
                return Err(format!("atom {}: P6 = {option}, the chains and rings of FullProf's manual (p. 36), is not calculated", a.label));
            }
            let satellite = if k == 2 {
                let parent = ((option.abs().fract() * 100.0).round() as usize).checked_sub(1).ok_or(format!("atom {}: P6 = {option} names no parent group", a.label))?;
                if parent >= out.len() {
                    return Err(format!("atom {}: P6 = {option}: the satellite of group {} comes before it", a.label, parent + 1));
                }
                let p16 = l3[7];
                Some(Satellite { parent, atom: p16.trunc() as usize, from: (p16.fract() * 100.0).round() as usize })
            } else {
                None
            };
            current = Some(Group {
                option,
                coordinates,
                degrees: l3[8] == 1.0,
                origin: [l3[3], l3[4], l3[5]],
                theta: l1.get(1).copied().unwrap_or(0.0),
                phi: l1.get(2).copied().unwrap_or(0.0),
                chi: l3[6],
                satellite,
                members: Vec::new(),
                displacements: None,
            });
        }
        let g = current.as_mut().ok_or(format!("atom {}: a rigid-body phase's atom outside a group", a.label))?;
        if a.rigid.1.len() < 3 {
            return Err(format!("atom {}: its internal coordinates (line 25-3) are missing", a.label));
        }
        let first = g.members.is_empty();
        let connect = if g.coordinates == Coordinates::ZMatrix && !first { std::array::from_fn(|k| a.rigid.0.get(k).copied().unwrap_or(0.0) as usize) } else { [0; 3] };
        if connect.iter().any(|&c| c > g.members.len()) {
            return Err(format!("atom {}: its connectivity {connect:?} names an atom after it", a.label));
        }
        g.members.push(Member { label: a.label.clone(), index: i, internal: [a.rigid.1[0], a.rigid.1[1], a.rigid.1[2]], connect });
    }
    out.extend(current);
    Ok(out)
}

/// Every atom's fractional coordinates: a group's placed from its values, a free atom's as written.
pub fn positions(ph: &PcrPhase) -> Result<Vec<[f64; 3]>, String> {
    let c = ph.cell;
    let cell = Cell::new(c[0], c[1], c[2], c[3], c[4], c[5]).map_err(|e| format!("the cell: {e:?}"))?;
    let mut x: Vec<[f64; 3]> = ph.atoms.iter().map(|a| a.x).collect();
    for g in groups(ph)? {
        for (m, p) in g.members.iter().zip(g.fractional(&cell)) {
            x[m.index] = p;
        }
    }
    Ok(x)
}
