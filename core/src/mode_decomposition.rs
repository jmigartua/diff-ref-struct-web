//! Symmetry modes computed from a parent and a distorted structure (P13, R29; L263): what AMPLIMODES (Orayech et al.
//! 2015, J. Appl. Cryst. 48, 318–333) and ISODISTORT (Campbell et al. 2006, J. Appl. Cryst. 39, 607–614) do, the
//! decomposition of a distortion into the irreducible representations (irreps) of the parent's space group G that the
//! distorted phase's group H keeps.
//!
//! The method, which needs no table of irreps:
//! 1. **A supercell the whole of G keeps.** The distorted phase's lattice T need not be kept by G's rotations (a
//!    tetragonal cell in a cubic parent), so the displacements are taken periodic in nZ³, n times the parent's
//!    conventional lattice, n the least integer with n·P⁻¹ integral (every nZ³ vector a lattice vector of the
//!    distorted phase). G acts on the displacements of the supercell's atoms through the finite group G/nZ³ (the
//!    parent's operators times the supercell's translations), as orthogonal matrices D(g) in Cartesian coordinates of
//!    the parent's metric.
//! 2. **The reference structure:** the parent's atoms in the supercell. The distorted structure is brought to the
//!    parent's basis by the transformation (P, p) (ITA's convention, x_parent = P x + p), expanded by H over the
//!    supercell, and each reference atom is paired with the nearest atom of the same element: u its displacement.
//! 3. **The irreps.** A sum Σ_C c_C Σ_{g∈C} D(g) over the conjugacy classes C of G/nZ³, with fixed pseudo-random c_C,
//!    commutes with every D(g); made symmetric (adding its transpose, the classes of the inverses), its eigenspaces
//!    are the isotypic components of the displacement representation, one per physically irreducible representation
//!    (an irrep and its complex conjugate together, as ISOTROPY and AMPLIMODES count them). The same c_C on every
//!    parent orbit (Wyckoff orbit) gives the same eigenvalue to the same irrep, so the components of the orbits are
//!    matched by their eigenvalue.
//! 4. **The modes:** in each orbit's component of each irrep, the vectors H leaves fixed (the projector
//!    (1/|H|)Σ_{h∈H} D(h), which commutes with the component's), an orthonormal basis of them chosen canonically
//!    (the projections of the unit displacements of the supercell's atoms, in order, made orthonormal). Their number is
//!    the irrep's number of modes on that orbit (AMPLIMODES's "Symmetry Modes Summary").
//! 5. **Amplitudes,** with AMPLIMODES's normalisation (`modes`): the inner product (1/N)Σ u·v over the supercell's
//!    atoms, N the number of the distorted phase's primitive cells in it, in Å in the reference metric. A mode's
//!    amplitude is ⟨u, mode⟩; an irrep's is the norm of u's projection on its modes; the global distortion is |u|.
//! 6. **What names the irrep:** its k-vector star (the k of the parent's reciprocal lattice whose projector
//!    (1/|T|)Σ_τ cos(2πk·τ)D(τ), over the supercell's translations, is non-zero on the component), its parity under the
//!    inversion at the parent's origin when G has it, and its dimension. The dimension comes from the component's
//!    character χ: m the multiplicity, from Frobenius–Schur's (1/|G|)Σχ(g²) (m for a real irrep, 0 for a complex
//!    pair, −2m for a quaternionic one) and (1/|G|)Σχ(g)²; the dimension is the component's dimension over m.
//!    Without tables of irreps the number in a label such as GM4+ is not computed: an irrep is named from its k, parity
//!    and dimension ("GM+[3]a"), and takes AMPLIMODES's or ISODISTORT's name when their modes are given
//!    (`name_by_modes`).
//!
//! Checked in `mode_decomposition` against every analysis of AMPLIMODES in the course's files (`manuals/bcs/…/
//! 04_tools_materials/0_1…0_4`): the irreps, their k, parity, dimension, number of modes and amplitude, the global
//! distortion; and AMPLIMODES's mode vectors in the computed subspaces.

use crate::cell::Cell;
use crate::mat3::{self, Mat3};
use crate::modes::{Mode, ModeParent, ModeVector, Modes};
use crate::spacegroup::SpaceGroup;
use std::collections::HashMap;

/// Translations in 1/24 of the parent's conventional edges (1/2, 1/3, 1/4, 1/6, 1/8, 1/12).
const U: i64 = 24;

/// An atom of a structure's asymmetric unit.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Atom {
    pub label: String,
    /// The element ("O", "Sr"); compared without case or charge.
    pub species: String,
    pub position: [f64; 3],
}

/// A structure: its group (in the setting its coordinates use), cell and asymmetric unit.
#[derive(Clone, Debug)]
pub struct Structure {
    pub group: SpaceGroup,
    pub cell: Cell,
    pub atoms: Vec<Atom>,
}

/// One computed mode: its vectors on the reference structure's atoms (fractions of the distorted phase's cell, so
/// that Σ over its primitive cell of |V|² is 1 Å² in the reference metric, AMPLIMODES's normalisation) and the
/// distortion's amplitude along it (Å).
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ComputedMode {
    pub vectors: Vec<[f64; 3]>,
    pub amplitude: f64,
}

/// An irrep's modes on one parent orbit.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct OrbitPart {
    /// The parent atom (asymmetric unit) whose orbit this is.
    pub orbit: usize,
    pub modes: Vec<ComputedMode>,
    /// The norm of the distortion's projection on these modes (Å).
    pub amplitude: f64,
    /// The modes as fields over the supercell (Cartesian, Å; orthonormal in ⟨u, v⟩ = (1/N)Σu·v), for comparisons.
    #[cfg_attr(feature = "serde", serde(skip))]
    pub fields: Vec<Vec<[f64; 3]>>,
}

/// Whether an irrep is real, a complex pair or quaternionic (Frobenius–Schur).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Reality {
    Real,
    Complex,
    Quaternionic,
}

/// An irrep of the parent's group that the distorted phase's group keeps (at least one mode).
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct IrrepPart {
    /// Its name: computed ("GM+[3]a": k, parity, dimension, a letter among equals) or taken from a reference.
    pub name: String,
    /// The arms of its k-vector star, in the parent's conventional reciprocal basis, each in [0, 1) less the
    /// reciprocal lattice.
    pub star: Vec<[f64; 3]>,
    /// +1 or −1 under the inversion at the parent's origin, when G holds it and the irrep is even or odd.
    pub parity: Option<i32>,
    /// The dimension of the (physically irreducible) representation.
    pub dimension: usize,
    pub reality: Reality,
    pub modes: usize,
    pub amplitude: f64,
    pub orbits: Vec<OrbitPart>,
}

/// The decomposition of a distorted structure.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Decomposition {
    /// The reference structure in the distorted phase's basis: label, element, position, parent orbit.
    pub reference: Vec<(String, String, [f64; 3], usize)>,
    /// FullProf's Occ of each reference atom in the distorted group, m/M (a full site).
    pub occupancy: Vec<f64>,
    /// The parent's cell in the distorted phase's basis (the metric of the amplitudes).
    pub reference_cell: Cell,
    /// Each reference atom's displacement (fractions of the distorted phase's cell).
    pub displacements: Vec<[f64; 3]>,
    /// |u| (Å), and the largest atomic displacement (Å).
    pub global: f64,
    pub largest: f64,
    /// The uniform translation taken off the displacements (Å, Cartesian of the parent), non-zero only for a polar H.
    pub translation: [f64; 3],
    /// The part of u that H does not keep (Å): the rounding of the given coordinates.
    pub residual: f64,
    pub irreps: Vec<IrrepPart>,
    /// The supercell nZ³ and the orders of G/nZ³ and of H in it.
    pub supercell: usize,
    pub order: usize,
    pub order_h: usize,
    /// N, the distorted phase's primitive cells in the supercell.
    pub primitive_cells: f64,
    #[cfg_attr(feature = "serde", serde(skip))]
    work: Work,
}

/// What the comparisons need: the supercell's atoms and H.
#[derive(Clone, Debug, Default, PartialEq)]
struct Work {
    /// Each supercell atom: parent-frame position in [0, n)³ and its orbit.
    atoms: Vec<([f64; 3], usize)>,
    /// H's elements: rotation and translation (fractions) in the parent's frame.
    h: Vec<([[i32; 3]; 3], [f64; 3])>,
    a: Mat3,
    p: Mat3,
    origin: [f64; 3],
    n: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct El {
    r: [[i32; 3]; 3],
    t: [i64; 3],
}

fn element(s: &str) -> String {
    s.chars().take_while(|c| c.is_ascii_alphabetic()).collect::<String>().to_ascii_lowercase()
}

fn wrap(x: f64, n: f64) -> f64 {
    let y = x.rem_euclid(n);
    if (n - y).abs() < 1e-7 {
        0.0
    } else {
        y
    }
}

fn near(a: [f64; 3], b: [f64; 3], n: f64, tol: f64) -> bool {
    (0..3).all(|k| {
        let d = (a[k] - b[k]).rem_euclid(n);
        d.min(n - d) < tol
    })
}

fn imul(a: &[[i32; 3]; 3], b: &[[i32; 3]; 3]) -> [[i32; 3]; 3] {
    std::array::from_fn(|i| std::array::from_fn(|j| (0..3).map(|k| a[i][k] * b[k][j]).sum()))
}

fn ivec(a: &[[i32; 3]; 3], t: [i64; 3]) -> [i64; 3] {
    std::array::from_fn(|i| (0..3).map(|k| a[i][k] as i64 * t[k]).sum())
}

fn to_f(m: &[[i32; 3]; 3]) -> Mat3 {
    m.map(|r| r.map(|v| v as f64))
}

/// The Cartesian frame of a metric: A with AᵀA = G (upper triangular, Cholesky).
fn frame(g: &Mat3) -> Mat3 {
    let a00 = g[0][0].sqrt();
    let a01 = g[0][1] / a00;
    let a02 = g[0][2] / a00;
    let a11 = (g[1][1] - a01 * a01).sqrt();
    let a12 = (g[1][2] - a01 * a02) / a11;
    let a22 = (g[2][2] - a02 * a02 - a12 * a12).sqrt();
    [[a00, a01, a02], [0.0, a11, a12], [0.0, 0.0, a22]]
}

/// An eigenvector of an orbit's class-sum operator: (orbit, eigenvalue, the vector as a field over the supercell).
type Space = (usize, f64, Vec<[f64; 3]>);

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        0.5 + ((self.0 >> 11) as f64) / (1u64 << 53) as f64
    }
}

/// The decomposition of `distorted` (its coordinates and group in the distorted phase's setting) against `parent`, the
/// distorted phase's basis and origin being (a', b', c') = (a, b, c)P and p in the parent's (ITA; AMPLIMODES's
/// "Transformation matrix").
pub fn decompose(parent: &Structure, distorted: &Structure, p: &Mat3, origin: [f64; 3]) -> Result<Decomposition, String> {
    let pinv = mat3::inverse(p).ok_or("the transformation is singular")?;
    let det = mat3::det(p).abs();
    let n = (1..=24).find(|&n| pinv.iter().flatten().all(|v| (v * n as f64 - (v * n as f64).round()).abs() < 1e-6)).ok_or("the transformation's inverse has no denominator up to 24")?;
    let nf = n as f64;
    let nu = U * n as i64;
    let red = |t: [i64; 3]| t.map(|v| v.rem_euclid(nu));
    // G/nZ³.
    let mut els: Vec<El> = Vec::new();
    for op in &parent.group.ops {
        for v in 0..n * n * n {
            let v = [v / (n * n), (v / n) % n, v % n].map(|x| x as i64);
            els.push(El { r: op.r, t: red([0, 1, 2].map(|k| op.t[k] as i64 * (U / crate::symop::DEN as i64) + U * v[k])) });
        }
    }
    let index: HashMap<El, usize> = els.iter().enumerate().map(|(i, e)| (*e, i)).collect();
    if index.len() != els.len() {
        return Err("the parent's operators repeat".into());
    }
    let order = els.len();
    let mul = |a: &El, b: &El| El { r: imul(&a.r, &b.r), t: red(std::array::from_fn(|k| ivec(&a.r, b.t)[k] + a.t[k])) };
    let inv = |a: &El| {
        let ri = mat3::inverse(&to_f(&a.r)).unwrap().map(|r| r.map(|v| v.round() as i32));
        El { r: ri, t: red(ivec(&ri, a.t).map(|v| -v)) }
    };
    // H in the parent's frame: the distorted group's operators and its lattice's translations modulo nZ³.
    let to_parent_t = |t: [f64; 3], r: &Mat3| -> [f64; 3] {
        let pt = mat3::mul_vec(p, &t);
        let rp = mat3::mul_vec(r, &origin);
        std::array::from_fn(|k| pt[k] + origin[k] - rp[k])
    };
    let lattice: Vec<[f64; 3]> = {
        let mut v: Vec<[f64; 3]> = Vec::new();
        let m = (pinv.iter().flatten().fold(0.0f64, |a, b| a.max(b.abs())).ceil() as i64 + 1) * n as i64 + 1;
        for i in -m..=m {
            for j in -m..=m {
                for k in -m..=m {
                    let x = mat3::mul_vec(p, &[i as f64, j as f64, k as f64]).map(|c| wrap(c, nf));
                    if !v.iter().any(|y| near(*y, x, nf, 1e-6)) {
                        v.push(x);
                    }
                }
            }
        }
        v
    };
    let expected = nf.powi(3) / det;
    if (lattice.len() as f64 - expected).abs() > 1e-6 {
        return Err(format!("{} translations of the distorted lattice in the supercell, {expected} expected", lattice.len()));
    }
    let mut h: Vec<usize> = Vec::new();
    let mut h_ops: Vec<([[i32; 3]; 3], [f64; 3])> = Vec::new();
    for op in &distorted.group.ops {
        let rp = mat3::mul(p, &mat3::mul(&to_f(&op.r), &pinv));
        let ri = rp.map(|r| r.map(|v| v.round() as i32));
        if rp.iter().flatten().zip(ri.iter().flatten()).any(|(a, b)| (a - *b as f64).abs() > 1e-6) {
            return Err(format!("the distorted group's operator {} is no rotation of the parent's lattice in this transformation", op.to_xyz()));
        }
        let t0 = to_parent_t(op.t.map(|v| v as f64 / crate::symop::DEN as f64), &to_f(&ri));
        for l in &lattice {
            let t: [f64; 3] = std::array::from_fn(|k| t0[k] + l[k]);
            let tu = t.map(|v| v * U as f64);
            if tu.iter().any(|v| (v - v.round()).abs() > 1e-5) {
                return Err(format!("the distorted group's operator {} has a translation off the 1/24 grid in the parent's frame", op.to_xyz()));
            }
            let e = El { r: ri, t: red(tu.map(|v| v.round() as i64)) };
            let i = *index.get(&e).ok_or(format!("the distorted group's operator {} is not an operator of the parent: check the transformation and origin", op.to_xyz()))?;
            if !h.contains(&i) {
                h.push(i);
                h_ops.push((ri, std::array::from_fn(|k| wrap(t[k], nf))));
            }
        }
    }
    // The reference: the parent's atoms over the supercell.
    let mut atoms: Vec<([f64; 3], usize)> = Vec::new();
    for (o, at) in parent.atoms.iter().enumerate() {
        for e in &els {
            let y: [f64; 3] = std::array::from_fn(|i| wrap((0..3).map(|j| e.r[i][j] as f64 * at.position[j]).sum::<f64>() + e.t[i] as f64 / U as f64, nf));
            if !atoms.iter().any(|(x, _)| near(*x, y, nf, 1e-5)) {
                atoms.push((y, o));
            }
        }
    }
    let na = atoms.len();
    let find = |y: [f64; 3], o: usize| atoms.iter().position(|(x, oo)| *oo == o && near(*x, y, nf, 1e-4));
    // The permutation of each element.
    let mut perm: Vec<Vec<usize>> = Vec::with_capacity(order);
    for e in &els {
        let mut pm = Vec::with_capacity(na);
        for (x, o) in &atoms {
            let y: [f64; 3] = std::array::from_fn(|i| (0..3).map(|j| e.r[i][j] as f64 * x[j]).sum::<f64>() + e.t[i] as f64 / U as f64);
            pm.push(find(y, *o).ok_or("an operator maps a reference atom to no atom: the parent's atoms do not close")?);
        }
        perm.push(pm);
    }
    // The distorted structure over the supercell.
    let mut dist: Vec<([f64; 3], String)> = Vec::new();
    for at in &distorted.atoms {
        for op in &distorted.group.ops {
            let y = crate::site::SiteOp { r: op.r, t: op.t }.apply(at.position);
            let py = mat3::mul_vec(p, &y);
            let x: [f64; 3] = std::array::from_fn(|k| py[k] + origin[k]);
            for l in &lattice {
                let z: [f64; 3] = std::array::from_fn(|k| wrap(x[k] + l[k], nf));
                if !dist.iter().any(|(w, _)| near(*w, z, nf, 1e-5)) {
                    dist.push((z, element(&at.species)));
                }
            }
        }
    }
    if dist.len() != na {
        return Err(format!("the distorted structure has {} atoms in the supercell, the reference {na}", dist.len()));
    }
    let g = parent.cell.metric();
    let a = frame(&g);
    let ainv = mat3::inverse(&a).ok_or("singular metric")?;
    let mut u: Vec<[f64; 3]> = Vec::with_capacity(na);
    let mut used = vec![false; na];
    for (x, o) in &atoms {
        let sp = element(&parent.atoms[*o].species);
        let mut best: Option<(f64, usize, [f64; 3])> = None;
        for (j, (z, s)) in dist.iter().enumerate() {
            if *s != sp {
                continue;
            }
            let d: [f64; 3] = std::array::from_fn(|k| {
                let d = (z[k] - x[k]).rem_euclid(nf);
                if d > nf / 2.0 {
                    d - nf
                } else {
                    d
                }
            });
            let dd = mat3::quad(&g, &d).sqrt();
            if best.is_none_or(|b| dd < b.0) {
                best = Some((dd, j, d));
            }
        }
        let (_, j, d) = best.ok_or(format!("no {sp} in the distorted structure"))?;
        if used[j] {
            return Err(format!("two {sp} reference atoms pair with the same distorted atom: the structures do not match in this transformation"));
        }
        used[j] = true;
        u.push(mat3::mul_vec(&a, &d));
    }
    let primitive_cells = nf.powi(3) / (det / distorted.group.centring.len() as f64);
    let dot = |x: &[[f64; 3]], y: &[[f64; 3]]| x.iter().zip(y).map(|(p, q)| p[0] * q[0] + p[1] * q[1] + p[2] * q[2]).sum::<f64>() / primitive_cells;
    // The uniform translation (the origin of a polar group).
    let translation: [f64; 3] = std::array::from_fn(|k| u.iter().map(|v| v[k]).sum::<f64>() / na as f64);
    for v in u.iter_mut() {
        for k in 0..3 {
            v[k] -= translation[k];
        }
    }
    // The Cartesian rotations.
    let rc: Vec<Mat3> = els.iter().map(|e| mat3::mul(&a, &mat3::mul(&to_f(&e.r), &ainv))).collect();
    let act = |gi: usize, f: &[[f64; 3]]| -> Vec<[f64; 3]> {
        let mut out = vec![[0.0; 3]; na];
        for (i, v) in f.iter().enumerate() {
            out[perm[gi][i]] = mat3::mul_vec(&rc[gi], v);
        }
        out
    };
    let project_h = |f: &[[f64; 3]]| -> Vec<[f64; 3]> {
        let mut s = vec![[0.0; 3]; na];
        for &hi in &h {
            for (i, v) in act(hi, f).iter().enumerate() {
                for k in 0..3 {
                    s[i][k] += v[k] / h.len() as f64;
                }
            }
        }
        s
    };
    let uh = project_h(&u);
    let residual = dot(&u.iter().zip(&uh).map(|(x, y)| [x[0] - y[0], x[1] - y[1], x[2] - y[2]]).collect::<Vec<_>>(), &u.iter().zip(&uh).map(|(x, y)| [x[0] - y[0], x[1] - y[1], x[2] - y[2]]).collect::<Vec<_>>()).sqrt();
    // Conjugacy classes and the class weights.
    let mut class = vec![usize::MAX; order];
    let mut nclass = 0;
    for i in 0..order {
        if class[i] != usize::MAX {
            continue;
        }
        for gg in &els {
            let c = mul(&mul(gg, &els[i]), &inv(gg));
            class[index[&c]] = nclass;
        }
        nclass += 1;
    }
    let mut rng = Rng(0x5eed_2026);
    let weight: Vec<f64> = (0..nclass).map(|_| rng.next()).collect();
    // Per orbit: the symmetric class-sum operator and its eigenspaces.
    let orbits = parent.atoms.len();
    let mut spaces: Vec<Space> = Vec::new();
    let mut all_values: Vec<f64> = Vec::new();
    for o in 0..orbits {
        let idx: Vec<usize> = (0..na).filter(|&i| atoms[i].1 == o).collect();
        let pos: HashMap<usize, usize> = idx.iter().enumerate().map(|(k, &i)| (i, k)).collect();
        let d = 3 * idx.len();
        let mut m = vec![vec![0.0; d]; d];
        for gi in 0..order {
            let c = weight[class[gi]];
            for (ka, &ia) in idx.iter().enumerate() {
                let kb = pos[&perm[gi][ia]];
                for r in 0..3 {
                    for s in 0..3 {
                        let v = c * rc[gi][r][s];
                        m[3 * kb + r][3 * ka + s] += v;
                        m[3 * ka + s][3 * kb + r] += v;
                    }
                }
            }
        }
        let (vals, vecs) = crate::lsq::symmetric_eigen(&m);
        for (j, &lam) in vals.iter().enumerate() {
            let mut f = vec![[0.0; 3]; na];
            for (k, &i) in idx.iter().enumerate() {
                f[i] = [vecs[3 * k][j], vecs[3 * k + 1][j], vecs[3 * k + 2][j]];
            }
            spaces.push((o, lam, f));
            all_values.push(lam);
        }
    }
    // Cluster the eigenvalues: one cluster per irrep.
    let scale = all_values.iter().fold(1.0f64, |a, b| a.max(b.abs()));
    let mut sorted = all_values.clone();
    sorted.sort_by(f64::total_cmp);
    let mut centres: Vec<f64> = Vec::new();
    for v in sorted {
        if centres.last().is_none_or(|c| (v - c).abs() > 1e-6 * scale) {
            centres.push(v);
        }
    }
    let cluster = |v: f64| centres.iter().position(|c| (v - c).abs() <= 1e-6 * scale + 1e-12).unwrap();
    let mut irreps: Vec<IrrepPart> = Vec::new();
    let translations: Vec<usize> = (0..order).filter(|&i| els[i].r == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]).collect();
    let inversion = index.get(&El { r: [[-1, 0, 0], [0, -1, 0], [0, 0, -1]], t: [0, 0, 0] }).copied();
    // Elements squared, for Frobenius–Schur.
    let square: Vec<usize> = els.iter().map(|e| index[&mul(e, e)]).collect();
    let cand_k: Vec<[f64; 3]> = {
        let cen: Vec<[f64; 3]> = parent.group.centring.iter().map(|t| t.map(|v| v as f64 / crate::symop::DEN as f64)).collect();
        let in_recip = |d: [f64; 3]| d.iter().all(|v| (v - v.round()).abs() < 1e-6) && cen.iter().all(|t| {
            let s: f64 = (0..3).map(|k| d[k] * t[k]).sum();
            (s - s.round()).abs() < 1e-6
        });
        let mut ks: Vec<[f64; 3]> = Vec::new();
        let m = 6 * n;
        for i in 0..m {
            for j in 0..m {
                for k in 0..m {
                    let kv = [i, j, k].map(|q| q as f64 / nf);
                    if !ks.iter().any(|q| in_recip(std::array::from_fn(|c| kv[c] - q[c]))) {
                        ks.push(kv);
                    }
                }
            }
        }
        ks
    };
    for (ci, _) in centres.iter().enumerate() {
        let comp: Vec<&Space> = spaces.iter().filter(|s| cluster(s.1) == ci).collect();
        let dim_e = comp.len();
        // The character over G.
        let chi: Vec<f64> = (0..order).map(|gi| comp.iter().map(|s| dot(&act(gi, &s.2), &s.2) * primitive_cells).sum()).collect();
        let fs: f64 = (0..order).map(|gi| chi[square[gi]]).sum::<f64>() / order as f64;
        let nn: f64 = chi.iter().map(|c| c * c).sum::<f64>() / order as f64;
        let (reality, m) = if fs > 0.5 {
            (Reality::Real, fs.round())
        } else if fs < -0.5 {
            (Reality::Quaternionic, (-fs / 2.0).round())
        } else {
            (Reality::Complex, (nn / 2.0).sqrt().round())
        };
        let expect_nn = match reality {
            Reality::Real => m * m,
            Reality::Complex => 2.0 * m * m,
            Reality::Quaternionic => 4.0 * m * m,
        };
        if (nn - expect_nn).abs() > 1e-4 * expect_nn.max(1.0) {
            return Err(format!("an eigenspace of the class sums is not one irrep (⟨χ,χ⟩ {nn:.4}, {expect_nn} expected): two irreps share an eigenvalue"));
        }
        let dimension = (dim_e as f64 / m).round() as usize;
        // The fixed vectors of H, orbit by orbit.
        let mut parts: Vec<OrbitPart> = Vec::new();
        for o in 0..orbits {
            let e: Vec<&Vec<[f64; 3]>> = comp.iter().filter(|s| s.0 == o).map(|s| &s.2).collect();
            if e.is_empty() {
                continue;
            }
            // Canonical basis: the projections of the unit displacements onto the fixed subspace, orthonormalised.
            let ph: Vec<Vec<[f64; 3]>> = e.iter().map(|f| project_h(f)).collect();
            // The fixed subspace is spanned by ph (P_H maps the component onto its fixed part).
            let mut basis: Vec<Vec<[f64; 3]>> = Vec::new();
            let span_proj = |f: &[[f64; 3]], basis: &[Vec<[f64; 3]>]| -> Vec<[f64; 3]> {
                let mut r = f.to_vec();
                for b in basis {
                    let c = dot(&r, b);
                    for (x, y) in r.iter_mut().zip(b) {
                        for k in 0..3 {
                            x[k] -= c * y[k];
                        }
                    }
                }
                r
            };
            // An orthonormal basis of span(ph) first.
            let mut fixed: Vec<Vec<[f64; 3]>> = Vec::new();
            for f in &ph {
                let r = span_proj(f, &fixed);
                let nr = dot(&r, &r).sqrt();
                if nr > 1e-6 {
                    fixed.push(r.iter().map(|v| v.map(|x| x / nr)).collect());
                }
            }
            if fixed.is_empty() {
                continue;
            }
            // The canonical basis: unit displacements in order, projected on the fixed subspace.
            for i in 0..na {
                if atoms[i].1 != o || basis.len() == fixed.len() {
                    continue;
                }
                for k in 0..3 {
                    if basis.len() == fixed.len() {
                        break;
                    }
                    let mut unit = vec![[0.0; 3]; na];
                    unit[i][k] = 1.0;
                    let mut q = vec![[0.0; 3]; na];
                    for b in &fixed {
                        let c = dot(&unit, b);
                        for (x, y) in q.iter_mut().zip(b) {
                            for kk in 0..3 {
                                x[kk] += c * y[kk];
                            }
                        }
                    }
                    let r = span_proj(&q, &basis);
                    let nr = dot(&r, &r).sqrt();
                    if nr > 1e-6 {
                        basis.push(r.iter().map(|v| v.map(|x| x / nr)).collect());
                    }
                }
            }
            let modes: Vec<ComputedMode> = basis.iter().map(|b| ComputedMode { vectors: Vec::new(), amplitude: dot(&u, b) }).collect();
            let amplitude = modes.iter().map(|m| m.amplitude * m.amplitude).sum::<f64>().sqrt();
            parts.push(OrbitPart { orbit: o, modes, amplitude, fields: basis });
        }
        if parts.is_empty() {
            continue;
        }
        let modes = parts.iter().map(|p| p.modes.len()).sum();
        let amplitude = parts.iter().map(|p| p.amplitude * p.amplitude).sum::<f64>().sqrt();
        // The star of k.
        let star: Vec<[f64; 3]> = cand_k
            .iter()
            .filter(|k| {
                let w: f64 = translations
                    .iter()
                    .map(|&ti| {
                        let t = els[ti].t.map(|v| v as f64 / U as f64);
                        (2.0 * std::f64::consts::PI * (0..3).map(|c| k[c] * t[c]).sum::<f64>()).cos() * chi[ti]
                    })
                    .sum::<f64>()
                    / translations.len() as f64;
                w > 0.5
            })
            .map(|k| k.map(|v| if v.abs() < 1e-12 { 0.0 } else { v }))
            .collect();
        let parity = inversion.and_then(|ii| {
            let s = chi[ii] / dim_e as f64;
            if (s - 1.0).abs() < 1e-6 {
                Some(1)
            } else if (s + 1.0).abs() < 1e-6 {
                Some(-1)
            } else {
                None
            }
        });
        irreps.push(IrrepPart { name: String::new(), star, parity, dimension, reality, modes, amplitude, orbits: parts });
    }
    // Names: k, parity, dimension, and a letter among irreps that share them.
    irreps.sort_by(|x, y| {
        let kx = x.star.iter().map(|k| k[0] + k[1] + k[2]).fold(f64::MAX, f64::min);
        let ky = y.star.iter().map(|k| k[0] + k[1] + k[2]).fold(f64::MAX, f64::min);
        kx.total_cmp(&ky).then(y.parity.cmp(&x.parity)).then(x.dimension.cmp(&y.dimension)).then(y.amplitude.total_cmp(&x.amplitude))
    });
    let base: Vec<String> = irreps.iter().map(|r| format!("{}{}[{}]", k_name(&r.star), r.parity.map_or("", |p| if p > 0 { "+" } else { "-" }), r.dimension)).collect();
    for i in 0..irreps.len() {
        let same: Vec<usize> = (0..irreps.len()).filter(|&j| base[j] == base[i]).collect();
        irreps[i].name = if same.len() > 1 { format!("{}{}", base[i], (b'a' + same.iter().position(|&j| j == i).unwrap() as u8) as char) } else { base[i].clone() };
    }
    // The reference structure in the distorted phase's basis: one atom per H-orbit, and the modes' vectors on them.
    let mut rep_of: Vec<Option<usize>> = vec![None; na];
    let mut reference: Vec<(String, String, [f64; 3], usize)> = Vec::new();
    let mut ref_atom: Vec<usize> = Vec::new();
    let mut count: HashMap<usize, usize> = HashMap::new();
    for i in 0..na {
        if rep_of[i].is_some() {
            continue;
        }
        let r = reference.len();
        for &hi in &h {
            rep_of[perm[hi][i]] = Some(r);
        }
        let o = atoms[i].1;
        let c = count.entry(o).or_insert(0);
        *c += 1;
        let label = if *c == 1 { parent.atoms[o].label.clone() } else { format!("{}_{}", parent.atoms[o].label, c) };
        let x = atoms[i].0;
        let y = mat3::mul_vec(&pinv, &std::array::from_fn(|k| x[k] - origin[k])).map(|v| wrap(v, 1.0)).map(|v| if v.abs() < 1e-12 { 0.0 } else { v });
        reference.push((label, parent.atoms[o].species.clone(), y, o));
        ref_atom.push(i);
    }
    let to_child = |v: [f64; 3]| mat3::mul_vec(&pinv, &mat3::mul_vec(&ainv, &v));
    for r in irreps.iter_mut() {
        for part in r.orbits.iter_mut() {
            for (m, f) in part.modes.iter_mut().zip(&part.fields) {
                // Normalised in (1/N)Σ over the supercell, as AMPLIMODES's V over the primitive cell.
                m.vectors = ref_atom.iter().map(|&i| to_child(f[i])).collect();
            }
        }
    }
    let displacements: Vec<[f64; 3]> = ref_atom.iter().map(|&i| to_child(u[i])).collect();
    let occupancy = reference.iter().map(|r| crate::structure_factor::occupancy_from_chemical(&distorted.group, r.2, 1.0)).collect::<Result<Vec<f64>, String>>()?;
    let global = dot(&u, &u).sqrt();
    let largest = u.iter().map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()).fold(0.0, f64::max);
    let reference_cell = crate::modes::transform_cell(&parent.cell, p)?;
    Ok(Decomposition {
        reference,
        occupancy,
        reference_cell,
        displacements,
        global,
        largest,
        translation,
        residual,
        irreps,
        supercell: n,
        order,
        order_h: h.len(),
        primitive_cells,
        work: Work { atoms, h: h_ops, a, p: *p, origin, n },
    })
}

/// The name of a k-vector star: GM for Γ, else its first arm.
fn k_name(star: &[[f64; 3]]) -> String {
    if star.len() == 1 && star[0] == [0.0, 0.0, 0.0] {
        return "GM".into();
    }
    let k = star.iter().min_by(|a, b| (a[0] + a[1] + a[2]).total_cmp(&(b[0] + b[1] + b[2])).then(b[1].total_cmp(&a[1]))).copied().unwrap_or([0.0; 3]);
    format!("k({})", k.map(fmt_frac).join(","))
}

fn fmt_frac(v: f64) -> String {
    for d in [1, 2, 3, 4, 6, 8, 12] {
        let n = v * d as f64;
        if (n - n.round()).abs() < 1e-6 {
            return if d == 1 { format!("{}", n.round() as i64) } else { format!("{}/{d}", n.round() as i64) };
        }
    }
    format!("{v:.4}")
}

impl Decomposition {
    /// A field over the supercell (Cartesian, Å) from vectors on atoms of the distorted phase (positions and vectors
    /// in its basis, fractions), expanded by H: AMPLIMODES's V_MODES or an ISODISTORT mode, for comparison. None when
    /// a position is no atom of the reference.
    pub fn field(&self, atoms: &[([f64; 3], [f64; 3])]) -> Option<Vec<[f64; 3]>> {
        let w = &self.work;
        let nf = w.n as f64;
        let mut f = vec![[0.0; 3]; w.atoms.len()];
        let mut set = vec![false; w.atoms.len()];
        for (x, v) in atoms {
            let xp: [f64; 3] = {
                let y = mat3::mul_vec(&w.p, x);
                std::array::from_fn(|k| y[k] + w.origin[k])
            };
            let vp = mat3::mul_vec(&w.p, v);
            for (r, t) in &w.h {
                let y: [f64; 3] = std::array::from_fn(|i| (0..3).map(|j| r[i][j] as f64 * xp[j]).sum::<f64>() + t[i]);
                let i = w.atoms.iter().position(|(z, _)| near(*z, y, nf, 1e-4))?;
                let rv = mat3::mul_vec(&to_f(r), &vp);
                if !set[i] {
                    f[i] = mat3::mul_vec(&w.a, &rv);
                    set[i] = true;
                }
            }
        }
        Some(f)
    }

    /// ⟨x, y⟩ = (1/N)Σ x·y over the supercell (Å²).
    pub fn inner(&self, x: &[[f64; 3]], y: &[[f64; 3]]) -> f64 {
        x.iter().zip(y).map(|(p, q)| p[0] * q[0] + p[1] * q[1] + p[2] * q[2]).sum::<f64>() / self.primitive_cells
    }

    /// The part of a field outside the modes of irrep `r` on orbit `o` (Å): 0 when the field is one of their
    /// combinations.
    pub fn outside(&self, r: usize, o: usize, f: &[[f64; 3]]) -> f64 {
        let part = self.irreps[r].orbits.iter().find(|p| p.orbit == o);
        let mut rest = f.to_vec();
        if let Some(part) = part {
            for b in &part.fields {
                let c = self.inner(f, b);
                for (x, y) in rest.iter_mut().zip(b) {
                    for k in 0..3 {
                        x[k] -= c * y[k];
                    }
                }
            }
        }
        self.inner(&rest, &rest).sqrt()
    }

    /// Each irrep named after the modes given (AMPLIMODES's block for FullProf, a `.pcr`'s V_MODES or an ISODISTORT
    /// file), when each of its modes lies in exactly one computed irrep's modes; the names kept otherwise. Returns the
    /// number of irreps named.
    pub fn name_by_modes(&mut self, modes: &Modes) -> usize {
        let mut named = 0;
        let mut names: Vec<Option<String>> = vec![None; self.irreps.len()];
        for m in &modes.modes {
            let atoms: Vec<([f64; 3], [f64; 3])> = m.vectors.iter().map(|v| (modes.reference[v.atom], v.v)).collect();
            let Some(f) = self.field(&atoms) else { continue };
            let norm = self.inner(&f, &f).sqrt();
            if norm < 1e-9 {
                continue;
            }
            let hits: Vec<usize> = (0..self.irreps.len())
                .filter(|&r| {
                    let total: f64 = self.irreps[r].orbits.iter().flat_map(|p| p.fields.iter()).map(|b| self.inner(&f, b).powi(2)).sum();
                    (total.sqrt() - norm).abs() < 1e-3 * norm.max(1e-3)
                })
                .collect();
            if let [r] = hits[..] {
                if names[r].is_none() {
                    names[r] = Some(m.irrep.clone());
                    named += 1;
                }
            }
        }
        for (r, nm) in names.into_iter().enumerate() {
            if let Some(nm) = nm {
                self.irreps[r].name = nm;
            }
        }
        named
    }

    /// The modes as the model refines them (`Modes`): the reference structure, each computed mode named
    /// "A<m>_<irrep>" with its vectors and amplitude, the reference cell and the parent.
    pub fn to_modes(&self, parent_group: &str) -> Modes {
        let mut modes = Vec::new();
        for r in &self.irreps {
            for part in &r.orbits {
                for m in &part.modes {
                    let k = modes.len() + 1;
                    let vectors = m.vectors.iter().enumerate().filter(|(_, v)| v.iter().any(|x| x.abs() > 1e-12)).map(|(atom, v)| ModeVector { atom, v: *v, coeff: 1.0 }).collect();
                    modes.push(Mode { name: format!("A{k}_{}", r.name), irrep: r.name.clone(), amplitude: m.amplitude, vectors });
                }
            }
        }
        Modes {
            reference: self.reference.iter().map(|r| r.2).collect(),
            modes,
            reference_cell: Some(self.reference_cell),
            max_amplitude: None,
            fst: Vec::new(),
            parent: Some(ModeParent { group: parent_group.to_string(), basis: self.work.p, origin: self.work.origin }),
            sign_reference: None,
        }
    }
}

/// The structures and transformation of an AMPLIMODES analysis (`modes::read_amplimodes`): the parent (its number's
/// default setting), the distorted structure (in the setting the analysis names, else its number's default) and
/// (P, p).
pub fn from_amplimodes(a: &crate::modes::Analysis) -> Result<(Structure, Structure, Mat3, [f64; 3]), String> {
    let hs = a.high.as_ref().ok_or("the analysis has no high-symmetry structure")?;
    let ls = a.low.as_ref().ok_or("the analysis has no low-symmetry structure")?;
    let (p, o) = a.transformation.ok_or("the analysis has no transformation matrix")?;
    let st = |s: &crate::modes::AmpliStructure, symbol: &str| -> Result<Structure, String> {
        let (group, _) = SpaceGroup::from_symbol(symbol)?;
        let c = s.cell;
        let cell = Cell::new(c[0], c[1], c[2], c[3], c[4], c[5]).map_err(|e| e.to_string())?;
        let atoms = s.atoms.iter().map(|(l, _, x)| Atom { label: l.clone(), species: l.chars().take_while(|c| c.is_ascii_alphabetic()).collect(), position: *x }).collect();
        Ok(Structure { group, cell, atoms })
    };
    let parent = st(hs, &hs.group.to_string())?;
    let distorted = st(ls, a.low_setting.as_deref().unwrap_or(&ls.group.to_string()))?;
    Ok((parent, distorted, p, o))
}

/// A basis written as ITA and AMPLIMODES write it, the new vectors in the old: "c,b,-a-c", "1/2a+1/2b,-1/2a+1/2b,c",
/// "a-b,a+b,2c". The columns of the matrix returned are the new vectors' components.
pub fn parse_basis(s: &str) -> Result<Mat3, String> {
    let parts: Vec<&str> = s.split(',').map(str::trim).collect();
    if parts.len() != 3 {
        return Err(format!("\"{s}\": three vectors separated by commas, as \"a-b,a+b,c\""));
    }
    let mut p = [[0.0; 3]; 3];
    for (j, e) in parts.iter().enumerate() {
        let e: String = e.chars().filter(|c| !c.is_whitespace()).collect();
        let mut terms: Vec<String> = Vec::new();
        for (i, ch) in e.char_indices() {
            if (ch == '+' || ch == '-') && i > 0 {
                terms.push(String::new());
            }
            match terms.last_mut() {
                Some(t) => t.push(ch),
                None => terms.push(ch.to_string()),
            }
        }
        for t in terms {
            let axis = t.chars().last().ok_or(format!("\"{e}\": an empty term"))?;
            let i = "abc".find(axis).ok_or(format!("\"{e}\": a term must end in a, b or c"))?;
            let coef = &t[..t.len() - 1];
            let coef = coef.trim_end_matches('*');
            let v = match coef {
                "" | "+" => 1.0,
                "-" => -1.0,
                c => {
                    let (sign, c) = match c.strip_prefix('-') {
                        Some(r) => (-1.0, r),
                        None => (1.0, c.trim_start_matches('+')),
                    };
                    sign * match c.split_once('/') {
                        Some((n, d)) => n.parse::<f64>().map_err(|_| format!("\"{t}\""))? / d.parse::<f64>().map_err(|_| format!("\"{t}\""))?,
                        None => c.parse::<f64>().map_err(|_| format!("\"{t}\""))?,
                    }
                }
            };
            p[i][j] += v;
        }
    }
    if mat3::det(&p).abs() < 1e-9 {
        return Err(format!("\"{s}\": the vectors are not independent"));
    }
    Ok(p)
}

/// A structure for the decomposition from a CIF (`cif::structure`).
pub fn from_cif(text: &str) -> Result<Structure, String> {
    let s = crate::cif::structure(text)?;
    Ok(Structure { group: s.group, cell: s.cell, atoms: s.atoms.into_iter().map(|a| Atom { label: a.label, species: a.species, position: a.x }).collect() })
}


impl Decomposition {
    /// The block for FullProf's `.pcr` as AMPLIMODES writes it: the reference atoms, `V_MODES` (mode, atom, irrep,
    /// vector, coefficient) and `A_MODES` (each mode's amplitude).
    pub fn fullprof_block(&self) -> String {
        let mut s = String::from("!Atom   Typ       X         Y         Z     Biso      Occ     In Fin N_t Spc /Codes\n");
        for ((label, species, x, _), occ) in self.reference.iter().zip(&self.occupancy) {
            s += &format!("{label:<6} {:<6} {:9.6} {:9.6} {:9.6} 0.500000 {occ:8.6}  0   0   0    1\n                   0.00      0.00      0.00     0.00     0.00\n", species.to_uppercase(), x[0], x[1], x[2]);
        }
        // Names FullProf reads (letters, digits, +, −, _): "k(1/2,1/2,1/2)+[3]" as "k121212+_3".
        let mut m = self.to_modes("");
        let safe = |x: &str| x.replace('[', "_").chars().filter(|c| c.is_ascii_alphanumeric() || matches!(c, '+' | '-' | '_')).collect::<String>();
        for mode in m.modes.iter_mut() {
            (mode.name, mode.irrep) = (safe(&mode.name), safe(&mode.irrep));
        }
        let lines: usize = m.modes.iter().map(|x| x.vectors.len()).sum();
        s += &format!("! Polarisation Vectors of Symmetry Modes for each atom\nV_MODES {lines:>4}\n! Nm Atm     Irrep      Vx        Vy        Vz         Coeff\n");
        for (k, mode) in m.modes.iter().enumerate() {
            for v in &mode.vectors {
                s += &format!("{:>3} {:<6} {:<8} {:9.6} {:9.6} {:9.6} 1.00\n", k + 1, self.reference[v.atom].0, mode.irrep, v.v[0], v.v[1], v.v[2]);
            }
        }
        s += &format!("!Amplitudes of Symmetry Modes\nA_MODES {:>4}   2\n", m.modes.len());
        for mode in &m.modes {
            s += &format!("    {:<10} {:10.6}     1.00\n                0.00\n", mode.name, mode.amplitude);
        }
        s
    }
}

#[cfg(test)]
mod basis_tests {
    #[test]
    fn bases() {
        assert_eq!(super::parse_basis("c,b,-a-c").unwrap(), [[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, -1.0]]);
        assert_eq!(super::parse_basis("1/2a+1/2b, -1/2a+1/2b, c").unwrap(), [[0.5, -0.5, 0.0], [0.5, 0.5, 0.0], [0.0, 0.0, 1.0]]);
        assert_eq!(super::parse_basis("a-b,a+b,2c").unwrap(), [[1.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 2.0]]);
        assert!(super::parse_basis("a,a,c").is_err());
        assert!(super::parse_basis("a,b").is_err());
    }
}
