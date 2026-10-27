//! Distances, angles and bond-valence sums (P12; L081): what FullProf lists after a refinement through its program
//! Bond_Str (`CODFILn.dis`, `CODFILn_sum.bvs`; manual p. 11), from a structure and the standard deviations of its cell
//! and coordinates.
//!
//! - Distances: from each atom of the asymmetric unit to every atom the group's operators and the lattice
//!   translations place within `dmax`, each image once (CrysFML's `Calc_Dist_Angle_Sigma`).
//! - Angles: at each atom, between each pair of its neighbours within `dmax_angle`, with the triangle's other two.
//! - Standard deviations: first order, from the standard deviations of the cell and of the coordinates taken as
//!   independent (no covariances are given with them), each quantity differentiated as a whole. Bond_Str's differ
//!   (L081), kept for comparison: it leaves out the cell's standard deviations, adds the two ends of a distance as
//!   independent even when they are the same atom (`bond_str_sigma`), and takes an angle's from the triangle's three
//!   distances as if independent (`bond_str_angle_sigma`, CrysFML's `Angle_and_Sigma`).
//! - Bond valence: s = exp((d0 − d)/B0) (Brown & Altermatt 1985), d0 and B0 from CrysFML's table (`bv_table`). As
//!   Bond_Str (`Calc_BVS`): a cation–anion pair is a bond when d ≤ (r_cation + r_anion)(1 + tol), the ionic radii of
//!   the same table, tol 0.20 by default; each bond's valence is weighted by the neighbour's occupancy; the sum's
//!   standard deviation is from the distances'. The occupancies are the chemical ones (`chemical_occupancies`); Bond_Str
//!   normalises FullProf's Occ to the first atom of the list (`bond_str_occupancies`), which is the same when the first
//!   atom's site is full.

use crate::bv_table::{ANIONS, ANION_RADII, CATIONS, IONIC_RADII, REFERENCES};
use crate::cell::Cell;
use crate::site::Site;
use crate::structure_factor::Structure;
use crate::symop::DEN;

/// Standard deviations of a structure: the cell (a, b, c in Å; α, β, γ in degrees) and each atom's coordinates.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct Uncertainties {
    pub cell: [f64; 6],
    pub x: Vec<[f64; 3]>,
}

/// A neighbour of an atom: the atom of the asymmetric unit it is an image of, the operator (an index into the
/// group's `ops`) and lattice translation that place it, its position, and the distance with its standard deviation.
#[derive(Clone, Debug, PartialEq)]
pub struct Contact {
    pub atom: usize,
    pub op: usize,
    pub translation: [i32; 3],
    pub position: [f64; 3],
    pub distance: f64,
    pub sigma: f64,
    /// A cation–anion contact's plausible distances, as `state::plausible` gives them for the pair at the cation's
    /// coordination at that state (bond-valence single-bond distances by coordination, the sum of the ionic radii);
    /// None for other contacts, or where `contacts` alone is asked (Explore's E18b, `state::summary_at` sets it).
    pub plausible: Option<Plausible>,
}

/// A pair's plausible distances (`state::plausible`): per coordination number the bond-valence single-bond distance
/// (None without parameters for the pair), and the sum of the ionic radii.
pub type Plausible = (Vec<(usize, Option<f64>)>, f64);

/// An angle at an atom between two of its neighbours (indices into that atom's contacts), with the triangle's other
/// two angles (at the first and at the second neighbour) and the distance between the neighbours, each with its
/// standard deviation.
#[derive(Clone, Debug, PartialEq)]
pub struct Angle {
    pub centre: usize,
    pub first: usize,
    pub second: usize,
    pub angle: (f64, f64),
    pub at_first: (f64, f64),
    pub at_second: (f64, f64),
    pub between: (f64, f64),
}

/// The cell's parameters in the order of `Uncertainties::cell`.
fn cell_values(c: &Cell) -> [f64; 6] {
    [c.a, c.b, c.c, c.alpha, c.beta, c.gamma]
}

fn cell_from(v: &[f64; 6]) -> Cell {
    Cell { a: v[0], b: v[1], c: v[2], alpha: v[3], beta: v[4], gamma: v[5] }
}

fn length(c: &Cell, d: [f64; 3]) -> f64 {
    let g = c.metric();
    let mut s = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            s += d[i] * g[i][j] * d[j];
        }
    }
    s.max(0.0).sqrt()
}

/// One point of a geometric quantity: an atom of the asymmetric unit taken through a rotation R and a translation.
#[derive(Clone, Copy)]
struct Point {
    atom: usize,
    rot: [[i32; 3]; 3],
    shift: [f64; 3],
}

impl Point {
    fn at(&self, x: &[[f64; 3]]) -> [f64; 3] {
        let a = x[self.atom];
        let mut y = self.shift;
        for (i, yi) in y.iter_mut().enumerate() {
            *yi += self.rot[i][0] as f64 * a[0] + self.rot[i][1] as f64 * a[1] + self.rot[i][2] as f64 * a[2];
        }
        y
    }
}

/// The first-order standard deviation of f(cell, coordinates) with independent standard deviations: central
/// differences in each parameter the points depend on.
fn propagate(st: &Structure, esd: Option<&Uncertainties>, atoms: &[usize], f: &dyn Fn(&Cell, &[[f64; 3]]) -> f64) -> f64 {
    let Some(u) = esd else { return 0.0 };
    let cell = cell_values(&st.cell);
    let x: Vec<[f64; 3]> = st.atoms.iter().map(|a| a.x).collect();
    let mut var = 0.0;
    for k in 0..6 {
        if u.cell[k] > 0.0 {
            let h = 1e-6 * cell[k].abs().max(1.0);
            let (mut p, mut m) = (cell, cell);
            p[k] += h;
            m[k] -= h;
            let d = (f(&cell_from(&p), &x) - f(&cell_from(&m), &x)) / (2.0 * h);
            var += (d * u.cell[k]).powi(2);
        }
    }
    let mut seen: Vec<usize> = atoms.to_vec();
    seen.sort_unstable();
    seen.dedup();
    for &a in &seen {
        for k in 0..3 {
            let s = u.x.get(a).map_or(0.0, |v| v[k]);
            if s > 0.0 {
                let h = 1e-7;
                let (mut p, mut m) = (x.clone(), x.clone());
                p[a][k] += h;
                m[a][k] -= h;
                let d = (f(&st.cell, &p) - f(&st.cell, &m)) / (2.0 * h);
                var += (d * s).powi(2);
            }
        }
    }
    var.sqrt()
}

fn point(st: &Structure, c: &Contact) -> Point {
    let op = &st.group.ops[c.op];
    Point { atom: c.atom, rot: op.r, shift: std::array::from_fn(|i| op.t[i] as f64 / DEN as f64 + c.translation[i] as f64) }
}

fn origin(atom: usize) -> Point {
    Point { atom, rot: [[1, 0, 0], [0, 1, 0], [0, 0, 1]], shift: [0.0; 3] }
}

fn distance_between(p: Point, q: Point) -> impl Fn(&Cell, &[[f64; 3]]) -> f64 {
    move |c: &Cell, x: &[[f64; 3]]| {
        let (a, b) = (p.at(x), q.at(x));
        length(c, [b[0] - a[0], b[1] - a[1], b[2] - a[2]])
    }
}

fn angle_at(centre: Point, p: Point, q: Point) -> impl Fn(&Cell, &[[f64; 3]]) -> f64 {
    move |c: &Cell, x: &[[f64; 3]]| {
        let (o, a, b) = (centre.at(x), p.at(x), q.at(x));
        // In Cartesian axes, by atan2 (|u × v|, u·v): exact at 0° and 180°, where the arc cosine's slope is infinite and
        // an angle fixed by symmetry would get an su from rounding (P7, L114).
        let m = c.crystal_to_orthonormal();
        let cart = |w: [f64; 3]| -> [f64; 3] { std::array::from_fn(|i| (0..3).map(|j| m[i][j] * w[j]).sum()) };
        let u = cart([a[0] - o[0], a[1] - o[1], a[2] - o[2]]);
        let v = cart([b[0] - o[0], b[1] - o[1], b[2] - o[2]]);
        let cross = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
        let n = (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt();
        n.atan2(u[0] * v[0] + u[1] * v[1] + u[2] * v[2]).to_degrees()
    }
}

/// The neighbours of each atom of the asymmetric unit within `dmax` (Å), in the order CrysFML lists them: by atom, then
/// operator, then translation.
pub fn contacts(st: &Structure, esd: Option<&Uncertainties>, dmax: f64) -> Vec<Vec<Contact>> {
    let cell = &st.cell;
    let lengths = [cell.a, cell.b, cell.c];
    let reach: [i32; 3] = std::array::from_fn(|k| (dmax / lengths[k] + 2.5).round() as i32);
    // Each component of a vector is at most the distance times the reciprocal length (CrysFML's first filter).
    let r = cell.reciprocal();
    let recip = [r.a, r.b, r.c];
    let mut out = Vec::new();
    for i in 0..st.atoms.len() {
        let xo = st.atoms[i].x;
        let mut list = Vec::new();
        for k in 0..st.atoms.len() {
            let mut seen: Vec<[f64; 3]> = vec![xo];
            for (j, op) in st.group.ops.iter().enumerate() {
                let p = Point { atom: k, rot: op.r, shift: std::array::from_fn(|m| op.t[m] as f64 / DEN as f64) };
                let xx = p.at(&st.atoms.iter().map(|a| a.x).collect::<Vec<_>>());
                for n1 in -reach[0]..=reach[0] {
                    for n2 in -reach[1]..=reach[1] {
                        'n3: for n3 in -reach[2]..=reach[2] {
                            let t = [n1, n2, n3];
                            let x1: [f64; 3] = std::array::from_fn(|m| xx[m] + t[m] as f64);
                            for m in 0..3 {
                                if (x1[m] - xo[m]).abs() / recip[m] > dmax {
                                    continue 'n3;
                                }
                            }
                            if seen.iter().any(|s| (0..3).map(|m| (s[m] - x1[m]).abs()).sum::<f64>() <= 1e-5) {
                                continue;
                            }
                            let d = length(cell, [x1[0] - xo[0], x1[1] - xo[1], x1[2] - xo[2]]);
                            if !(0.001..=dmax).contains(&d) {
                                continue;
                            }
                            seen.push(x1);
                            let c = Contact { atom: k, op: j, translation: t, position: x1, distance: d, sigma: 0.0, plausible: None };
                            let sigma = propagate(st, esd, &[i, k], &distance_between(origin(i), point(st, &c)));
                            list.push(Contact { sigma, ..c });
                        }
                    }
                }
            }
        }
        out.push(list);
    }
    out
}

/// The angles at each atom between each pair of its neighbours within `dmax_angle` (Å), in CrysFML's order.
pub fn angles(st: &Structure, esd: Option<&Uncertainties>, contacts: &[Vec<Contact>], dmax_angle: f64) -> Vec<Angle> {
    let mut out = Vec::new();
    for (i, list) in contacts.iter().enumerate() {
        for j in 0..list.len() {
            if list[j].distance > dmax_angle {
                continue;
            }
            for k in j + 1..list.len() {
                if list[k].distance > dmax_angle {
                    continue;
                }
                let (o, p, q) = (origin(i), point(st, &list[j]), point(st, &list[k]));
                let atoms = [i, list[j].atom, list[k].atom];
                let between = distance_between(p, q);
                let d12 = between(&st.cell, &st.atoms.iter().map(|a| a.x).collect::<Vec<_>>());
                if d12 < 1e-4 {
                    continue;
                }
                let value = |f: &dyn Fn(&Cell, &[[f64; 3]]) -> f64| (f(&st.cell, &st.atoms.iter().map(|a| a.x).collect::<Vec<_>>()), propagate(st, esd, &atoms, f));
                out.push(Angle {
                    centre: i,
                    first: j,
                    second: k,
                    angle: value(&angle_at(o, p, q)),
                    at_first: value(&angle_at(p, o, q)),
                    at_second: value(&angle_at(q, o, p)),
                    between: value(&between),
                });
            }
        }
    }
    out
}

/// Bond_Str's standard deviation of a distance (CrysFML's `Distance_and_Sigma` as FullProf's Bond_Str prints it):
/// the two ends' coordinates as independent, each with its own atom's standard deviations, unrotated, and no term
/// for the cell (an Mn–Mn distance of c/2 is printed without one, L081). Apart from the cell's small terms it is the
/// default's except for an atom bonded to its own image (the same coordinates move both ends) or through an
/// operator that exchanges axes.
pub fn bond_str_distance_sigma(st: &Structure, esd: &Uncertainties, centre: usize, c: &Contact) -> f64 {
    bond_str_sigma(st, esd, (centre, st.atoms[centre].x), (c.atom, c.position))
}

/// Bond_Str's standard deviation of the distance between two points, each an atom's image at a position.
pub fn bond_str_sigma(st: &Structure, esd: &Uncertainties, a: (usize, [f64; 3]), b: (usize, [f64; 3])) -> f64 {
    let d: [f64; 3] = std::array::from_fn(|m| b.1[m] - a.1[m]);
    let dist = length(&st.cell, d);
    let mut var = 0.0;
    let g = st.cell.metric();
    for m in 0..3 {
        let df = (0..3).map(|n| g[m][n] * d[n]).sum::<f64>() / dist;
        let (s0, s1) = (esd.x.get(a.0).map_or(0.0, |v| v[m]), esd.x.get(b.0).map_or(0.0, |v| v[m]));
        var += df * df * (s0 * s0 + s1 * s1);
    }
    var.sqrt()
}

/// Bond_Str's standard deviation (degrees) of the angle opposite the side d3 in a triangle of sides d1, d2, d3 with
/// standard deviations s1, s2, s3 (CrysFML's `Angle_and_Sigma`): √((s1/d1)² + (s2/d2)² + (s3 d3/(d1 d2))²) in radians.
pub fn bond_str_angle_sigma(d1: (f64, f64), d2: (f64, f64), d3: (f64, f64)) -> f64 {
    ((d1.1 / d1.0).powi(2) + (d2.1 / d2.0).powi(2) + (d3.1 * d3.0 / (d1.0 * d2.0)).powi(2)).sqrt().to_degrees()
}

// Bond valence ---------------------------------------------------------------------------------------------------------

/// A species as CrysFML writes it ("La+3" → ("LA", 3); "O-2" → ("O", −2)); None without a charge.
pub fn ion(species: &str) -> Option<(String, i32)> {
    let s = species.trim().to_uppercase();
    let at = s.find(['+', '-'])?;
    let (el, rest) = s.split_at(at);
    let n: i32 = rest[1..].trim().parse().ok().filter(|n| *n > 0).unwrap_or(1);
    let el: String = el.chars().filter(|c| c.is_ascii_alphabetic()).collect();
    (!el.is_empty()).then(|| (el, if rest.starts_with('-') { -n } else { n }))
}

/// An atom's ion: its element and formal charge, from `Atom::charge` (P5) or, without one, from its species ("SR+2").
pub fn atom_ion(a: &crate::structure_factor::Atom) -> Option<(String, i32)> {
    match a.charge {
        Some(q) => {
            let el: String = a.species.chars().take_while(|c| c.is_ascii_alphabetic()).collect::<String>().to_uppercase();
            (!el.is_empty()).then_some((el, q))
        }
        None => ion(&a.species),
    }
}

/// The charges of a structure's atoms set from a list by label or by element ("Y" → 3, "O" → −2; the caller's, as
/// a publication set's `[sample] charges`); a label wins over an element. Atoms named in neither keep theirs.
pub fn with_charges(st: &Structure, charges: &[(String, i32)]) -> Structure {
    let mut out = st.clone();
    for a in out.atoms.iter_mut() {
        let el: String = a.species.chars().take_while(|c| c.is_ascii_alphabetic()).collect();
        if let Some((_, q)) = charges.iter().find(|(k, _)| *k == a.label).or_else(|| charges.iter().find(|(k, _)| k.eq_ignore_ascii_case(&el))) {
            a.charge = Some(*q);
        }
    }
    out
}

pub(crate) fn ion_name(el: &str, q: i32) -> String {
    format!("{el}{}{}", if q < 0 { "-" } else { "+" }, q.abs())
}

/// The ionic radius (Å) of an ion in CrysFML's table: an anion's from `ANION_RADII`, a cation's by its oxidation state;
/// 0.8 Å when the table has none (as `Set_Atoms_Conf` puts it).
pub fn ionic_radius(el: &str, q: i32) -> f64 {
    if q < 0 {
        if let Some(i) = ANIONS.iter().position(|a| *a == ion_name(el, q)) {
            return ANION_RADII[i];
        }
    }
    let r = IONIC_RADII.iter().find(|e| e.0 == el).and_then(|e| e.1.iter().position(|&o| o as i32 == q).map(|k| e.2[k])).unwrap_or(0.0);
    if r < 0.01 {
        0.8
    } else {
        r
    }
}

/// The bond-valence parameters (d0, B0, reference) of a cation–anion pair in CrysFML's table.
pub fn parameters(cation: &str, anion: &str) -> Option<(f64, f64, &'static str)> {
    let (ce, cq) = ion(cation)?;
    let (ae, aq) = ion(anion)?;
    let a = ANIONS.iter().position(|x| *x == ion_name(&ae, aq))?;
    let c = CATIONS.iter().find(|e| e.0 == ion_name(&ce, cq))?;
    (c.1[a] > 0.0).then(|| (c.1[a], c.2[a], REFERENCES[c.3[a] as usize]))
}

/// The chemical occupancy of each atom from FullProf's Occ: Occ × M/m (m the site's multiplicity, M the general
/// position's; manual p. 102), unless that puts an atom above 1 and Occ × c/m does not (c the number of centring
/// vectors): the Occ written as atoms per primitive cell, the manual's "normalised to the multiplicity of the general
/// position" (LaB6rx.pcr: La 1.0, B 6.0 in P m -3 m; ybfluo_Riet_IRF.pcr: Yb 0.5, O 2.0 in F m -3 m, L081). The scale
/// factor absorbs either reading in the pattern; distances do not depend on it; bond-valence sums and a CIF do.
pub fn chemical_occupancies(st: &Structure) -> Result<Vec<f64>, String> {
    let m_g = st.group.multiplicity() as f64;
    let c = st.group.centring.len().max(1) as f64;
    let m: Vec<f64> = st.atoms.iter().map(|a| Site::new(&st.group, a.x).map(|s| s.multiplicity as f64)).collect::<Result<_, _>>()?;
    let documented: Vec<f64> = st.atoms.iter().zip(&m).map(|(a, m)| a.occ * m_g / m).collect();
    let primitive: Vec<f64> = st.atoms.iter().zip(&m).map(|(a, m)| a.occ * c / m).collect();
    let full = |v: &[f64]| v.iter().all(|x| *x <= 1.0 + 1e-3);
    Ok(if !full(&documented) && full(&primitive) { primitive } else { documented })
}

/// Bond_Str's occupancies: the chemical ones divided by the first atom's (CrysFML's `Set_Atoms_Conf`, `fac1`): the
/// first atom is taken as fully occupied, whatever its Occ says.
pub fn bond_str_occupancies(st: &Structure) -> Result<Vec<f64>, String> {
    let c = chemical_occupancies(st)?;
    let f = c.first().copied().filter(|v| *v > 0.0).unwrap_or(1.0);
    Ok(c.iter().map(|v| v / f).collect())
}

/// One bond of an atom's valence sum: the contact (an index into the atom's contacts), its valence and standard
/// deviation, both weighted by the neighbour's occupancy.
#[derive(Clone, Debug, PartialEq)]
pub struct Bond {
    pub contact: usize,
    pub valence: f64,
    pub sigma: f64,
}

/// An atom's bond-valence sum, with Bond_Str's figures: the bonds, the coordination (their number and Σ of the
/// neighbours' occupancies), the mean distance and its standard deviation (√Σσ²/n), the distortion 10⁴(⟨d²⟩/⟨d⟩² − 1),
/// the formal valence, and the sum with its standard deviation.
#[derive(Clone, Debug, PartialEq)]
pub struct ValenceSum {
    pub atom: usize,
    pub bonds: Vec<Bond>,
    pub coordination: usize,
    pub effective_coordination: f64,
    pub mean_distance: (f64, f64),
    pub distortion: f64,
    pub valence: f64,
    pub sum: (f64, f64),
    /// What the sum is compared with, |q| (Bond_Str compares the sum with the formal valence's magnitude), and the range
    /// taken as plausible about it: |q|(1 ∓ r), r = `VALENCE_RANGE` unless the phase's control-file key
    /// `VALENCE_RANGE_KEY` gives another (Josu's answer 5 to Explore's §8; E18c).
    pub expected: f64,
    pub plausible: (f64, f64),
}

/// The relative range about |q| within which a valence sum is plausible, by default: V ± 10 % (Josu's answer 5 to
/// Explore's DESIGN.md §8, L230).
pub const VALENCE_RANGE: f64 = 0.10;

/// The phase's control-file key that overrides `VALENCE_RANGE` (`Phase::valence_sum_range`).
pub const VALENCE_RANGE_KEY: &str = "valence_sum_range";

impl ValenceSum {
    /// The plausible range for a relative range r about |q|.
    pub fn with_range(mut self, r: f64) -> ValenceSum {
        self.plausible = (self.expected * (1.0 - r), self.expected * (1.0 + r));
        self
    }
}

/// The global instability indices of Bond_Str (×100): the old one √(Σ(V_i − |q_i|)²/N) over the atoms of the
/// asymmetric unit, and the three normalised by the atoms in the cell (a: Σ|Δ|·occ·m; b: Σ 100|Δ/q|·occ·m, in %;
/// c: √(Σ Δ²·occ·m)).
#[derive(Clone, Debug, PartialEq)]
pub struct Instability {
    pub old: f64,
    pub a: f64,
    pub b: f64,
    pub c: f64,
}

/// The bond-valence sums of every atom with a charge, from its contacts, with these occupancies and the tolerance on
/// the sum of the ionic radii (0.20 is Bond_Str's default); and the instability indices.
pub fn valence_sums(st: &Structure, contacts: &[Vec<Contact>], occupancy: &[f64], tol: f64) -> Result<(Vec<ValenceSum>, Instability), String> {
    let ions: Vec<(String, i32)> = st.atoms.iter().map(|a| atom_ion(a).ok_or_else(|| format!("{}: the species {} has no charge (\"O-2\", \"La+3\", or the control file's charge = …)", a.label, a.species))).collect::<Result<_, _>>()?;
    let mut sums = Vec::new();
    let (mut del2, mut ga, mut gb, mut gc, mut total) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for (i, a) in st.atoms.iter().enumerate() {
        total += occupancy[i] * Site::new(&st.group, a.x)?.multiplicity as f64;
    }
    for (i, list) in contacts.iter().enumerate() {
        let (el, q) = &ions[i];
        let mut bonds = Vec::new();
        let (mut eff, mut s, mut var, mut dsum, mut d2, mut sd2) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        for (j, c) in list.iter().enumerate() {
            let (en, qn) = &ions[c.atom];
            if q.signum() == qn.signum() {
                continue;
            }
            if c.distance > (ionic_radius(el, *q) + ionic_radius(en, *qn)) * (1.0 + tol) {
                continue;
            }
            let (cat, an) = if *q > 0 { (ion_name(el, *q), ion_name(en, *qn)) } else { (ion_name(en, *qn), ion_name(el, *q)) };
            let (d0, b0, _) = parameters(&cat, &an).ok_or_else(|| format!("no bond-valence parameters for {cat}–{an} in CrysFML's table"))?;
            let occ = occupancy[c.atom];
            let v = (d0 - c.distance) / b0;
            let (val, sig) = (v.exp() * occ, v.exp() * c.sigma / b0 * occ);
            eff += occ;
            s += val;
            var += sig * sig;
            dsum += c.distance;
            d2 += c.distance * c.distance;
            sd2 += c.sigma * c.sigma;
            bonds.push(Bond { contact: j, valence: val, sigma: sig });
        }
        let n = bonds.len();
        if n == 0 {
            continue;
        }
        let nf = n as f64;
        let mean = dsum / nf;
        let distortion = 1e4 * (d2 / nf / (mean * mean) - 1.0);
        let del = s - q.abs() as f64;
        del2 += del * del;
        let fact = occupancy[i] * Site::new(&st.group, st.atoms[i].x)?.multiplicity as f64;
        ga += del.abs() * fact;
        gb += 100.0 * (del / *q as f64).abs() * fact;
        gc += del * del * fact;
        sums.push(ValenceSum { atom: i, bonds, coordination: n, effective_coordination: eff, mean_distance: (mean, sd2.sqrt() / nf), distortion, valence: *q as f64, sum: (s, var.sqrt()), expected: q.abs() as f64, plausible: (0.0, 0.0) }.with_range(VALENCE_RANGE));
    }
    let gii = Instability { old: (del2 / st.atoms.len() as f64).sqrt() * 100.0, a: ga * 100.0 / total, b: gb / total, c: (gc / total).sqrt() * 100.0 };
    Ok((sums, gii))
}
