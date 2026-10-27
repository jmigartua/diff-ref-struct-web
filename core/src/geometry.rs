//! A refined state's derived quantities with their standard uncertainties (Publish P3, P4, P7, P11, L114; Explore E3, E4,
//! L115): the positions reduced to the cell, the bonds and angles, the parent cell of a mode phase.
//!
//! **One way to every su** (`Propagator`): a quantity q computed from the model has σ²(q) = gᵀ C g, C the covariance
//! of the free parameters (`lsq::Refinement::covariance`, χ² included) and g its gradient with respect to them, taken
//! by central differences: each parameter moved by a thousandth of its esd, its values with their coefficients, as the
//! engine moves them. So a mode phase's positions and distances carry the amplitudes' covariance (P3), a value tied to
//! another carries the tie, and a value the symmetry fixes (an angle of 180° at an inversion centre, a coordinate on a
//! special position) has no gradient and is written without su (P7).
//!
//! **The bonds** (`geometry`): every contact within the rule's reach, each pair listed once (P7): between ions, from the
//! cation to the anion within 1.2 times the sum of their ionic radii (charges from `Atom::charge` or the species,
//! P5); without charges, every contact up to 3 Å, from the atom listed first. The angles at each atom between its
//! bonds. Labels unique within a phase (P8: two "Y" are Y1 and Y2).

use crate::bonds;
use crate::jacobian::Field;
use crate::lsq::Refinement;
use crate::powder::{Intensities, Model};
use crate::refinable::ParameterSet;
use crate::structure_factor::Structure;
use crate::symop::DEN;

/// The su's of a refined state's quantities, through the covariance of its free parameters.
#[derive(Clone, Debug)]
pub struct Propagator {
    model: Model,
    set: ParameterSet,
    cov: Vec<Vec<f64>>,
}

impl Propagator {
    pub fn new(model: &Model, set: &ParameterSet, refinement: &Refinement) -> Propagator {
        Propagator::from_covariance(model, set, &refinement.covariance)
    }

    /// With a covariance of the free parameters given directly: a cycle's (`Cycle::covariance`, E14a).
    pub fn from_covariance(model: &Model, set: &ParameterSet, covariance: &[Vec<f64>]) -> Propagator {
        Propagator { model: model.clone(), set: set.clone(), cov: covariance.to_vec() }
    }

    /// With the covariance's diagonal only: the parameters taken as independent, as FullProf propagates the
    /// amplitudes of a mode phase to its positions (L090); for comparisons.
    pub fn independent(&self) -> Propagator {
        let n = self.cov.len();
        let cov = (0..n).map(|j| (0..n).map(|k| if j == k { self.cov[j][j] } else { 0.0 }).collect()).collect();
        Propagator { model: self.model.clone(), set: self.set.clone(), cov }
    }

    /// The values of `f` at the state and their su's, differentiating only by the parameters that move a value for
    /// which `moves` is true (all others have no effect on `f`). An su below 10⁻⁸ of its value (10⁻¹² absolute) is
    /// 0: the value is fixed by the symmetry or a constraint.
    pub fn propagate(&self, f: &dyn Fn(&Model) -> Vec<f64>, moves: &dyn Fn(Field) -> bool) -> Vec<(f64, f64)> {
        let v0 = f(&self.model);
        let n = self.set.free.len();
        let mut grads: Vec<(usize, Vec<f64>)> = Vec::new();
        for j in 0..n {
            let p = &self.set.free[j];
            let var = self.cov.get(j).and_then(|r| r.get(j)).copied().unwrap_or(0.0);
            if var <= 0.0 || !p.fields.iter().any(|(fd, _)| moves(*fd)) {
                continue;
            }
            let h = 1e-3 * var.sqrt();
            let moved = |s: f64| {
                let mut m = self.model.clone();
                for &(fd, c) in &p.fields {
                    m.set(fd, m.get(fd) + s * c * h);
                }
                f(&m)
            };
            let (up, down) = (moved(1.0), moved(-1.0));
            grads.push((j, up.iter().zip(&down).map(|(a, b)| (a - b) / (2.0 * h)).collect()));
        }
        v0.iter()
            .enumerate()
            .map(|(i, &v)| {
                let mut var = 0.0;
                for (j, gj) in &grads {
                    for (k, gk) in &grads {
                        var += gj[i] * self.cov[*j][*k] * gk[i];
                    }
                }
                let su = var.max(0.0).sqrt();
                (v, if su < 1e-8 * v.abs().max(1e-4) { 0.0 } else { su })
            })
            .collect()
    }

    /// As `propagate`, for a quantity of phase p: only the parameters that move its cell, atoms or modes count.
    pub fn of_phase(&self, p: usize, f: &dyn Fn(&Model) -> Vec<f64>) -> Vec<(f64, f64)> {
        self.propagate(f, &|fd| match fd {
            Field::Cell(q, _) | Field::AtomX(q, ..) | Field::AtomB(q, _) | Field::AtomOcc(q, _) | Field::AtomBeta(q, ..) | Field::Mode(q, _) | Field::Rigid(q, ..) => q == p,
            _ => false,
        })
    }
}

/// U_eq of an anisotropic atom (P15): a third of the trace of U in Cartesian axes, U_eq = ⅓ Σᵢⱼ U^ij a*ᵢ a*ⱼ aᵢ·aⱼ
/// (Fischer & Tillmanns 1988, Acta Cryst. C44, 775–776; the CIF dictionary's `_atom_site_U_iso_or_equiv`). With
/// FullProf's βᵢⱼ = 2π² a*ᵢ a*ⱼ Uᵢⱼ it is Σᵢⱼ βᵢⱼ Gᵢⱼ / (6π²), G the direct metric; β holds β11 β22 β33 β12 β13 β23.
pub fn u_equivalent(cell: &crate::cell::Cell, beta: &[f64; 6]) -> f64 {
    let g = cell.metric();
    let pi2 = std::f64::consts::PI * std::f64::consts::PI;
    (beta[0] * g[0][0] + beta[1] * g[1][1] + beta[2] * g[2][2] + 2.0 * (beta[3] * g[0][1] + beta[4] * g[0][2] + beta[5] * g[1][2])) / (6.0 * pi2)
}

/// Each atom's U_eq of phase p with its su through the covariance (P15): None for an isotropic atom.
pub fn u_eq(model: &Model, p: usize, prop: Option<&Propagator>) -> Vec<Option<(f64, f64)>> {
    let Some(st) = structure(model, p) else { return Vec::new() };
    let of = |m: &Model| -> Vec<f64> {
        let s = structure(m, p).unwrap();
        s.atoms.iter().map(|a| a.beta.map_or(0.0, |b| u_equivalent(&s.cell, &b))).collect()
    };
    let values: Vec<(f64, f64)> = match prop {
        Some(pr) => pr.of_phase(p, &of),
        None => of(model).into_iter().map(|v| (v, 0.0)).collect(),
    };
    st.atoms.iter().zip(values).map(|(a, v)| a.beta.map(|_| v)).collect()
}

/// The anisotropic U's (U11 U22 U33 U12 U13 U23 in Å², the CIF's) of the atoms a TLS group gives them (L367), with su's
/// through the covariance (T, L, S, the group's values and the cell): (atom, values), in the atoms' order.
pub fn tls_u(model: &Model, p: usize, prop: Option<&Propagator>) -> Vec<(usize, [(f64, f64); 6])> {
    let atoms = tls_atoms(model, p);
    if atoms.is_empty() || structure(model, p).is_none() {
        return Vec::new();
    }
    let idx = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];
    let of = |m: &Model| -> Vec<f64> {
        let s = structure(m, p).unwrap();
        let r = s.cell.reciprocal();
        let rs = [r.a, r.b, r.c];
        atoms.iter().flat_map(|&a| { let b = s.atoms[a].beta.unwrap_or([0.0; 6]); (0..6).map(move |k| b[k] / (2.0 * std::f64::consts::PI * std::f64::consts::PI * rs[idx[k].0] * rs[idx[k].1])) }).collect()
    };
    let values: Vec<(f64, f64)> = match prop {
        Some(pr) => pr.of_phase(p, &of),
        None => of(model).into_iter().map(|v| (v, 0.0)).collect(),
    };
    atoms.iter().enumerate().map(|(i, &a)| (a, std::array::from_fn(|k| values[6 * i + k]))).collect()
}

/// The atoms of phase p's TLS groups (L367).
fn tls_atoms(model: &Model, p: usize) -> Vec<usize> {
    let Some(gs) = model.phases.get(p).map(|ph| &ph.rigid) else { return Vec::new() };
    let mut out: Vec<usize> = (0..gs.len()).flat_map(|g| crate::rigid_body::tls_members(gs, g).into_iter().map(|(h, m)| gs[h].members[m].index)).collect();
    out.sort_unstable();
    out
}

/// The bond-valence sum of one atom (Publish P18): its coordination, mean distance and sum with su's.
#[derive(Clone, Debug, PartialEq)]
pub struct Valence {
    pub atom: usize,
    pub coordination: usize,
    /// The formal valence |q|.
    pub valence: f64,
    pub mean_distance: (f64, f64),
    pub sum: (f64, f64),
}

/// Phase p's bond-valence sums, when every atom has a charge (`bonds::atom_ion`; the control file's `charge =`): the
/// contacts within `dmax` (Å) as `bonds::valence_sums` takes them (tolerance `tol` on the sum of the ionic radii,
/// Bond_Str's occupancies). The su's are propagated through the covariance with the contacts held (`Propagator`),
/// where `bonds::valence_sums` takes the distances as independent. None when an atom has no charge.
pub fn valence(model: &Model, p: usize, prop: Option<&Propagator>, dmax: f64, tol: f64) -> Result<Option<Vec<Valence>>, String> {
    let (st, contacts, _) = contacts_with_su(model, p, None, dmax, 0.0)?;
    if !st.atoms.iter().all(|a| bonds::atom_ion(a).is_some()) {
        return Ok(None);
    }
    let occ = bonds::bond_str_occupancies(&st)?;
    let (sums, _) = bonds::valence_sums(&st, &contacts, &occ, tol)?;
    // The same bonds at another state: the contacts' distances recomputed, the list held.
    let at = |m: &Model| -> Vec<f64> {
        let Some(s) = structure(m, p) else { return vec![f64::NAN; 2 * sums.len()] };
        let shift: Vec<[f64; 3]> = s.atoms.iter().map(|a| a.x.map(|v| (v + 1e-12).floor())).collect();
        let mut red = s.clone();
        for (a, sh) in red.atoms.iter_mut().zip(&shift) {
            a.x = std::array::from_fn(|k| a.x[k] - sh[k]);
        }
        let occ = bonds::bond_str_occupancies(&red).unwrap_or_else(|_| vec![f64::NAN; red.atoms.len()]);
        let mut out = Vec::new();
        for v in &sums {
            let ion = bonds::atom_ion(&red.atoms[v.atom]);
            let (mut d, mut sum) = (0.0, 0.0);
            for b in &v.bonds {
                let c = &contacts[v.atom][b.contact];
                let bond = Bond { i: v.atom, j: c.atom, op: c.op, translation: c.translation, distance: (c.distance, 0.0) };
                let dist = norm(vector(&red, &bond));
                d += dist;
                let other = bonds::atom_ion(&red.atoms[c.atom]);
                let par = match (&ion, &other) {
                    // Cation first, as `bonds::valence_sums` looks the pair up.
                    (Some((e1, q1)), Some((e2, q2))) if *q1 > 0 => bonds::parameters(&bonds::ion_name(e1, *q1), &bonds::ion_name(e2, *q2)),
                    (Some((e1, q1)), Some((e2, q2))) => bonds::parameters(&bonds::ion_name(e2, *q2), &bonds::ion_name(e1, *q1)),
                    _ => None,
                };
                sum += par.map_or(f64::NAN, |(d0, b0, _)| ((d0 - dist) / b0).exp() * occ[c.atom]);
            }
            out.push(d / v.bonds.len() as f64);
            out.push(sum);
        }
        out
    };
    let values: Vec<(f64, f64)> = match prop {
        Some(pr) => pr.of_phase(p, &at),
        None => at(model).into_iter().map(|v| (v, 0.0)).collect(),
    };
    Ok(Some(sums.iter().enumerate().map(|(k, v)| Valence { atom: v.atom, coordination: v.coordination, valence: v.valence.abs(), mean_distance: values[2 * k], sum: values[2 * k + 1] }).collect()))
}

/// Phase p's structure as the model holds it (the cell the phase's, positions from its modes).
pub fn structure(model: &Model, p: usize) -> Option<Structure> {
    match &model.phases.get(p)?.intensities {
        Intensities::Structure(st, _) => {
            let mut s = (**st).clone();
            s.cell = model.phases[p].cell;
            Some(s)
        }
        Intensities::List(_) => None,
    }
}

/// Labels unique within a structure (P8): a label two atoms share is numbered on each ("Y", "Y" → "Y1", "Y2"; "O1",
/// "O1" → "O1_1", "O1_2"); a label used once stays.
pub fn unique_labels(st: &Structure) -> Vec<String> {
    let labels: Vec<&str> = st.atoms.iter().map(|a| a.label.trim()).collect();
    let mut out = Vec::new();
    for (i, l) in labels.iter().enumerate() {
        let same: Vec<usize> = (0..labels.len()).filter(|&k| labels[k] == *l).collect();
        if same.len() == 1 {
            out.push(l.to_string());
            continue;
        }
        let n = same.iter().position(|&k| k == i).unwrap() + 1;
        let sep = if l.ends_with(|c: char| c.is_ascii_digit()) { "_" } else { "" };
        out.push(format!("{l}{sep}{n}"));
    }
    // A numbered label that another atom already has ("Y1" given and "Y", "Y"): primes until unique.
    for i in 0..out.len() {
        while (0..out.len()).any(|k| k != i && out[k] == out[i]) {
            out[i].push('\'');
        }
    }
    out
}

/// A bond: from atom i to the image of atom j by operator `op` and lattice `translation`, its length (Å) with su.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Bond {
    pub i: usize,
    pub j: usize,
    pub op: usize,
    pub translation: [i32; 3],
    pub distance: (f64, f64),
}

/// An angle at atom `centre` between two of its bonds (the images of `first` and `second`), in degrees with su.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct BondAngle {
    pub centre: usize,
    pub first: Bond,
    pub second: Bond,
    pub angle: (f64, f64),
}

/// A phase's geometry at a refined state.
#[derive(Clone, Debug)]
pub struct Geometry {
    /// The structure with every position reduced to [0, 1) (the codes of the bonds refer to it).
    pub structure: Structure,
    pub labels: Vec<String>,
    /// Each atom's reduced position with its su's.
    pub positions: Vec<[(f64, f64); 3]>,
    /// The bonds, each pair once.
    pub bonds: Vec<Bond>,
    pub angles: Vec<BondAngle>,
    /// Whether every atom has a charge (the bonds then follow the ionic rule).
    pub charged: bool,
}

/// The image of atom j by operator `op` and `translation` in a structure (fractions).
pub fn image(st: &Structure, j: usize, op: usize, t: [i32; 3]) -> [f64; 3] {
    let o = &st.group.ops[op];
    let x = st.atoms[j].x;
    std::array::from_fn(|m| (0..3).map(|k| o.r[m][k] as f64 * x[k]).sum::<f64>() + o.t[m] as f64 / DEN as f64 + t[m] as f64)
}

fn cartesian(st: &Structure, f: [f64; 3]) -> [f64; 3] {
    let m = st.cell.crystal_to_orthonormal();
    std::array::from_fn(|r| (0..3).map(|c| m[r][c] * f[c]).sum())
}

fn vector(st: &Structure, b: &Bond) -> [f64; 3] {
    let (o, x) = (st.atoms[b.i].x, image(st, b.j, b.op, b.translation));
    cartesian(st, std::array::from_fn(|m| x[m] - o[m]))
}

fn norm(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// The angle between two vectors (degrees), by atan2 so that 0° and 180° are exact.
fn angle(u: [f64; 3], v: [f64; 3]) -> f64 {
    let cross = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
    norm(cross).atan2(u[0] * v[0] + u[1] * v[1] + u[2] * v[2]).to_degrees()
}

/// The rule that makes a contact a bond.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BondRule {
    /// Between ions, cation–anion within 1.2 times the sum of the ionic radii; without charges, up to 3 Å.
    Default,
    /// Every contact up to this distance (Å).
    Within(f64),
}

/// Phase p's geometry at the model's state: positions reduced to the cell, the bonds by `rule` each once, the angles
/// between each atom's bonds, all with su's when a propagator is given.
pub fn geometry(model: &Model, p: usize, prop: Option<&Propagator>, rule: BondRule) -> Result<Geometry, String> {
    let st0 = structure(model, p).ok_or(format!("phase {} has no structure", p + 1))?;
    // The reduction to [0, 1): an integer translation per atom, kept fixed while the su's are differentiated.
    let shift: Vec<[f64; 3]> = st0.atoms.iter().map(|a| a.x.map(|v| (v + 1e-12).floor())).collect();
    let reduce = |st: &Structure| -> Structure {
        let mut s = st.clone();
        for (a, sh) in s.atoms.iter_mut().zip(&shift) {
            a.x = std::array::from_fn(|k| a.x[k] - sh[k]);
        }
        s
    };
    let st = reduce(&st0);
    let ions: Vec<Option<(String, i32)>> = st.atoms.iter().map(bonds::atom_ion).collect();
    let charged = ions.iter().all(Option::is_some);
    let reach = match rule {
        BondRule::Within(d) => d,
        BondRule::Default if charged => 4.0,
        BondRule::Default => 3.0,
    };
    let is_bond = |i: usize, j: usize, d: f64| -> bool {
        match rule {
            BondRule::Within(r) => d <= r,
            BondRule::Default if charged => {
                let ((e1, q1), (e2, q2)) = (ions[i].clone().unwrap(), ions[j].clone().unwrap());
                q1.signum() != q2.signum() && d <= 1.2 * (bonds::ionic_radius(&e1, q1) + bonds::ionic_radius(&e2, q2))
            }
            BondRule::Default => d <= 3.0,
        }
    };
    let contacts = bonds::contacts(&st, None, reach);
    let as_bond = |i: usize, c: &bonds::Contact| Bond { i, j: c.atom, op: c.op, translation: c.translation, distance: (c.distance, 0.0) };
    // Each pair once (P7): from the cation (ions) or from the atom listed first; an atom to its own images by length.
    let mut list: Vec<Bond> = Vec::new();
    for (i, cs) in contacts.iter().enumerate() {
        for c in cs.iter().filter(|c| is_bond(i, c.atom, c.distance)) {
            let keep = match (&ions[i], &ions[c.atom], charged) {
                (Some((_, q1)), _, true) => *q1 > 0,
                _ => i < c.atom || (i == c.atom && !list.iter().any(|b| b.i == i && b.j == i && (b.distance.0 - c.distance).abs() < 1e-5)),
            };
            if keep {
                list.push(as_bond(i, c));
            }
        }
    }
    // The angles at every atom between its bonds (from either end), each pair of bonds once.
    let mut angs: Vec<BondAngle> = Vec::new();
    for (i, cs) in contacts.iter().enumerate() {
        let mine: Vec<Bond> = cs.iter().filter(|c| is_bond(i, c.atom, c.distance)).map(|c| as_bond(i, c)).collect();
        for a in 0..mine.len() {
            for b in a + 1..mine.len() {
                let v = angle(vector(&st, &mine[a]), vector(&st, &mine[b]));
                angs.push(BondAngle { centre: i, first: mine[a].clone(), second: mine[b].clone(), angle: (v, 0.0) });
            }
        }
    }
    // Every value at once, so that one model per moved parameter serves them all.
    let n_at = st.atoms.len();
    let values = |m: &Model| -> Vec<f64> {
        let s = reduce(&structure(m, p).unwrap());
        let mut v: Vec<f64> = s.atoms.iter().flat_map(|a| a.x).collect();
        v.extend(list.iter().map(|b| norm(vector(&s, b))));
        v.extend(angs.iter().map(|a| angle(vector(&s, &a.first), vector(&s, &a.second))));
        v
    };
    let vs: Vec<(f64, f64)> = match prop {
        Some(pr) => pr.of_phase(p, &values),
        None => values(model).into_iter().map(|v| (v, 0.0)).collect(),
    };
    let positions = (0..n_at).map(|a| std::array::from_fn(|k| vs[3 * a + k])).collect();
    for (b, v) in list.iter_mut().zip(&vs[3 * n_at..]) {
        b.distance = *v;
    }
    for (a, v) in angs.iter_mut().zip(&vs[3 * n_at + list.len()..]) {
        a.angle = *v;
    }
    Ok(Geometry { labels: unique_labels(&st), structure: st, positions, bonds: list, angles: angs, charged })
}

/// The CIF's site-symmetry code of a bond's far end ("." for the atom itself, else "op_klm" with 5 for no
/// translation).
pub fn symmetry_code(b: &Bond) -> String {
    let t = b.translation;
    if b.op == 0 && t == [0, 0, 0] {
        ".".into()
    } else {
        format!("{}_{}{}{}", b.op + 1, 5 + t[0], 5 + t[1], 5 + t[2])
    }
}

/// The cell of a phase reduced to its parent's (P11; P17 for a phase refined by its coordinates, its parent given): the phase's metric in the parent's basis (P⁻¹ of the
/// parent's transformation, `modes::ModeParent`), a, b, c, α, β, γ with su's when a propagator is given (case A's
/// P2₁/n cell, √2a_p × √2a_p × 2a_p, gives the pseudo-cubic a_p three times). None without a parent.
pub fn parent_cell(model: &Model, p: usize, prop: Option<&Propagator>) -> Option<[(f64, f64); 6]> {
    let ph = model.phases.get(p)?;
    let parent = ph.modes.as_ref().and_then(|m| m.parent.as_ref()).or(ph.parent.as_ref())?;
    let inv = crate::mat3::inverse(&parent.basis)?;
    let f = |m: &Model| -> Vec<f64> {
        crate::modes::transform_cell(&m.phases[p].cell, &inv).map_or(vec![f64::NAN; 6], |c| vec![c.a, c.b, c.c, c.alpha, c.beta, c.gamma])
    };
    let v = match prop {
        Some(pr) => pr.of_phase(p, &f),
        None => f(model).into_iter().map(|x| (x, 0.0)).collect(),
    };
    v.iter().all(|x| x.0.is_finite()).then(|| std::array::from_fn(|k| v[k]))
}

/// A structure reduced to the cell, its atoms' contacts and the angles between them (`contacts_with_su`).
pub type ContactsWithSu = (Structure, Vec<Vec<bonds::Contact>>, Vec<bonds::Angle>);

/// Every contact of phase p within `dmax` (Å) as `bonds::contacts` lists them (each atom's neighbours, Bond_Str's
/// order), and the angles between neighbours within `dangl`, at the model's state with positions reduced to the cell,
/// each distance's and angle's su propagated through the covariance when a propagator is given (Publish P4: `drs bonds`
/// on the refined state). The structure returned is the reduced one the contacts refer to.
pub fn contacts_with_su(model: &Model, p: usize, prop: Option<&Propagator>, dmax: f64, dangl: f64) -> Result<ContactsWithSu, String> {
    let st0 = structure(model, p).ok_or(format!("phase {} has no structure", p + 1))?;
    let shift: Vec<[f64; 3]> = st0.atoms.iter().map(|a| a.x.map(|v| (v + 1e-12).floor())).collect();
    let reduce = |st: &Structure| -> Structure {
        let mut s = st.clone();
        for (a, sh) in s.atoms.iter_mut().zip(&shift) {
            a.x = std::array::from_fn(|k| a.x[k] - sh[k]);
        }
        s
    };
    let st = reduce(&st0);
    let mut contacts = bonds::contacts(&st, None, dmax);
    let mut angles = if dangl > 0.0 { bonds::angles(&st, None, &contacts, dangl) } else { Vec::new() };
    let as_bond = |i: usize, c: &bonds::Contact| Bond { i, j: c.atom, op: c.op, translation: c.translation, distance: (c.distance, 0.0) };
    let list: Vec<Bond> = contacts.iter().enumerate().flat_map(|(i, cs)| cs.iter().map(move |c| as_bond(i, c))).collect();
    let pairs: Vec<(Bond, Bond)> = angles.iter().map(|a| (as_bond(a.centre, &contacts[a.centre][a.first]), as_bond(a.centre, &contacts[a.centre][a.second]))).collect();
    if let Some(pr) = prop {
        let values = |m: &Model| -> Vec<f64> {
            let s = reduce(&structure(m, p).unwrap());
            let mut v: Vec<f64> = list.iter().map(|b| norm(vector(&s, b))).collect();
            v.extend(pairs.iter().map(|(a, b)| angle(vector(&s, a), vector(&s, b))));
            v
        };
        let vs = pr.of_phase(p, &values);
        let mut k = 0;
        for cs in contacts.iter_mut() {
            for c in cs.iter_mut() {
                c.sigma = vs[k].1;
                k += 1;
            }
        }
        for a in angles.iter_mut() {
            a.angle.1 = vs[k].1;
            k += 1;
        }
    }
    Ok((st, contacts, angles))
}

// The scene of the 3D view (Explore E4) ---------------------------------------------------------------------------

/// An atom's image in the scene: its atom of the asymmetric unit, the operator (index into the group's `ops`) and
/// lattice translation that place it, its position in fractions and in Cartesian Å (`Cell::crystal_to_orthonormal`).
#[derive(Clone, Debug, PartialEq)]
pub struct Image {
    pub atom: usize,
    pub label: String,
    pub op: usize,
    pub translation: [i32; 3],
    pub fract: [f64; 3],
    pub cart: [f64; 3],
}

/// A face of a polyhedron: three images (indices into `Scene::images`), counter-clockwise seen from outside, and its
/// outward unit normal.
#[derive(Clone, Debug, PartialEq)]
pub struct Face {
    pub vertices: [usize; 3],
    pub normal: [f64; 3],
}

/// A coordination polyhedron: its centre and vertices (images), the faces of their convex hull (coplanar vertices'
/// faces fanned), its volume (Å³) and the distortion of its bonds, Baur's Δ = (1/n) Σ |dᵢ − ⟨d⟩| / ⟨d⟩ (Baur 1974,
/// Acta Cryst. B30, 1195).
#[derive(Clone, Debug, PartialEq)]
pub struct Polyhedron {
    pub centre: usize,
    pub vertices: Vec<usize>,
    pub faces: Vec<Face>,
    pub volume: f64,
    pub distortion: f64,
}

/// The rule that draws a bond in the scene.
#[derive(Clone, Debug, PartialEq)]
pub enum SceneRule {
    /// Between ions, cation–anion within this factor times the sum of the ionic radii (`bonds::ionic_radius`; 1.2 as
    /// `geometry`'s); without charges nothing.
    Ionic(f64),
    /// Every pair of images up to this distance (Å).
    Within(f64),
    /// By element pair (either order), a distance range (Å).
    Pairs(Vec<(String, String, f64, f64)>),
}

/// The scene of phase p (E4): every image whose fractions lie in the box [lo, hi] (the cell: [0, 0, 0]–[1, 1, 1]), each
/// once; the bonds between them by the rule; the polyhedra of the images with four bonds or more; the cell's twelve
/// edges; for a mode phase, each image's reference position and each irrep's displacement on it (Cartesian, Å: the
/// displacement of its atom turned by its operator), so that r = r⁰ + t·Σ u and the arrows are the core's.
#[derive(Clone, Debug, PartialEq)]
pub struct Scene {
    pub images: Vec<Image>,
    pub bonds: Vec<(usize, usize, f64)>,
    pub polyhedra: Vec<Polyhedron>,
    pub edges: Vec<([f64; 3], [f64; 3])>,
    pub reference: Vec<[f64; 3]>,
    /// Each image's reference position in fractions (P19), so that a slab or a region is chosen on the structure drawn.
    pub reference_fract: Vec<[f64; 3]>,
    pub displacements: Vec<(String, Vec<[f64; 3]>)>,
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(u: [f64; 3], v: [f64; 3]) -> [f64; 3] {
    [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
}

fn dot(u: [f64; 3], v: [f64; 3]) -> f64 {
    u[0] * v[0] + u[1] * v[1] + u[2] * v[2]
}

/// The faces of the convex hull of a few points (a coordination polyhedron): every plane through three points with
/// all the others on one side; the points on that plane ordered round their centre and fanned into triangles, each
/// plane once. Outward normals. Nothing for fewer than four points or all in a plane.
pub fn hull(points: &[[f64; 3]]) -> Vec<([usize; 3], [f64; 3])> {
    let n = points.len();
    let centre = points.iter().fold([0.0; 3], |a, p| [a[0] + p[0] / n as f64, a[1] + p[1] / n as f64, a[2] + p[2] / n as f64]);
    let scale = points.iter().map(|p| norm(sub(*p, centre))).fold(0.0, f64::max).max(1e-12);
    let eps = 1e-6 * scale;
    let mut planes: Vec<([f64; 3], f64)> = Vec::new();
    let mut faces = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            for k in j + 1..n {
                let nv = cross(sub(points[j], points[i]), sub(points[k], points[i]));
                let len = norm(nv);
                if len < 1e-9 * scale * scale {
                    continue;
                }
                let mut u = nv.map(|v| v / len);
                let mut off = dot(u, points[i]);
                if dot(u, centre) > off {
                    u = u.map(|v| -v);
                    off = -off;
                }
                if points.iter().any(|p| dot(u, *p) - off > eps) || planes.iter().any(|(q, o)| dot(*q, u) > 1.0 - 1e-9 && (o - off).abs() < eps) {
                    continue;
                }
                planes.push((u, off));
                // The plane's points round their centre, fanned from the first.
                let on: Vec<usize> = (0..n).filter(|&m| (dot(u, points[m]) - off).abs() <= eps).collect();
                let c = on.iter().fold([0.0; 3], |a, &m| [a[0] + points[m][0], a[1] + points[m][1], a[2] + points[m][2]]).map(|v| v / on.len() as f64);
                let e1 = sub(points[on[0]], c);
                let e1 = e1.map(|v| v / norm(e1).max(1e-300));
                let e2 = cross(u, e1);
                let mut ring: Vec<(f64, usize)> = on.iter().map(|&m| (dot(sub(points[m], c), e2).atan2(dot(sub(points[m], c), e1)), m)).collect();
                ring.sort_by(|a, b| a.0.total_cmp(&b.0));
                for t in 1..ring.len() - 1 {
                    faces.push(([ring[0].1, ring[t].1, ring[t + 1].1], u));
                }
            }
        }
    }
    faces
}

/// A phase refined by its coordinates described by symmetry modes (Publish P19): the modes of another file's mode
/// phase `source.phases[q]`, in the same setting (its group, its reference structure, its modes and parent). The
/// refined positions are compared with the reference (`modes::displacements`: each reference atom's nearest image of
/// the same element, after `shift`, within 1 Å) and projected on the modes (`Modes::project`, in `metric`: the given
/// one, else the source's reference cell, else the phase's cell). Returns the phase as a mode phase with those
/// amplitudes (its cell the refined one, its parent the source's or its own), so that `scene` gives its reference
/// and the irreps' displacements, and the residual (Å): the part of the displacements the modes do not describe.
pub fn decompose(model: &Model, p: usize, source: &Model, q: usize, metric: Option<crate::cell::Cell>, shift: [f64; 3]) -> Result<(Model, f64), String> {
    let st = structure(model, p).ok_or(format!("phase {} has no structure", p + 1))?;
    let sph = source.phases.get(q).ok_or(format!("the source has no phase {}", q + 1))?;
    let ms = sph.modes.as_ref().ok_or(format!("the source's phase {} has no modes", q + 1))?;
    let Intensities::Structure(sst, ssc) = &sph.intensities else { return Err("the source's phase has no structure".into()) };
    let metric = metric.or(ms.reference_cell).unwrap_or(model.phases[p].cell);
    let reference: Vec<(String, [f64; 3])> = sst.atoms.iter().zip(&ms.reference).map(|(a, x)| (a.species.clone(), *x)).collect();
    let distorted: Vec<(String, [f64; 3])> = st.atoms.iter().map(|a| (a.species.clone(), a.x)).collect();
    let u = crate::modes::displacements(&reference, &distorted, &sph.group, &metric, shift, 1.0)?;
    let (amps, residual) = ms.project(&sph.group, &metric, &u)?;
    let mut nms = ms.clone();
    for (m, a) in nms.modes.iter_mut().zip(&amps) {
        m.amplitude = *a;
    }
    if nms.parent.is_none() {
        nms.parent = model.phases[p].parent.clone();
    }
    let mut nst = (**sst).clone();
    nst.cell = model.phases[p].cell;
    for (a, x) in nst.atoms.iter_mut().zip(nms.positions()) {
        a.x = x;
    }
    let mut out = model.clone();
    let ph = &mut out.phases[p];
    ph.group = sph.group.clone();
    ph.intensities = Intensities::Structure(Box::new(nst), ssc.clone());
    ph.modes = Some(nms);
    Ok((out, residual))
}

/// Phase p's scene in the box [lo, hi] of fractions (E4).
pub fn scene(model: &Model, p: usize, lo: [f64; 3], hi: [f64; 3], rule: &SceneRule) -> Result<Scene, String> {
    let st = structure(model, p).ok_or(format!("phase {} has no structure", p + 1))?;
    let m = st.cell.crystal_to_orthonormal();
    let cart = |x: [f64; 3]| -> [f64; 3] { std::array::from_fn(|i| (0..3).map(|j| m[i][j] * x[j]).sum()) };
    let labels = unique_labels(&st);
    let tol = 1e-6;
    let mut images: Vec<Image> = Vec::new();
    for (a, atom) in st.atoms.iter().enumerate() {
        for (o, op) in st.group.ops.iter().enumerate() {
            let y: [f64; 3] = std::array::from_fn(|i| (0..3).map(|j| op.r[i][j] as f64 * atom.x[j]).sum::<f64>() + op.t[i] as f64 / DEN as f64);
            let range = |k: usize| ((lo[k] - y[k] - tol).ceil() as i32)..=((hi[k] - y[k] + tol).floor() as i32);
            for t0 in range(0) {
                for t1 in range(1) {
                    for t2 in range(2) {
                        let t = [t0, t1, t2];
                        let f: [f64; 3] = std::array::from_fn(|k| y[k] + t[k] as f64);
                        if images.iter().any(|im| im.atom == a && (0..3).all(|k| (im.fract[k] - f[k]).abs() < tol)) {
                            continue;
                        }
                        images.push(Image { atom: a, label: labels[a].clone(), op: o, translation: t, fract: f, cart: cart(f) });
                    }
                }
            }
        }
    }
    // The bonds.
    let ions: Vec<Option<(String, i32)>> = st.atoms.iter().map(bonds::atom_ion).collect();
    let element = |a: usize| st.atoms[a].species.chars().take_while(|c| c.is_ascii_alphabetic()).collect::<String>().to_ascii_lowercase();
    let mut bonds_out = Vec::new();
    for i in 0..images.len() {
        for j in i + 1..images.len() {
            let d = norm(sub(images[i].cart, images[j].cart));
            if d < 1e-6 {
                continue;
            }
            let (a, b) = (images[i].atom, images[j].atom);
            let ok = match rule {
                SceneRule::Within(r) => d <= *r,
                SceneRule::Ionic(f) => match (&ions[a], &ions[b]) {
                    (Some((ea, qa)), Some((eb, qb))) if qa.signum() != qb.signum() => d <= f * (bonds::ionic_radius(ea, *qa) + bonds::ionic_radius(eb, *qb)),
                    _ => false,
                },
                SceneRule::Pairs(list) => list.iter().any(|(x, y, d0, d1)| {
                    let (x, y) = (x.to_ascii_lowercase(), y.to_ascii_lowercase());
                    ((element(a) == x && element(b) == y) || (element(a) == y && element(b) == x)) && d >= *d0 && d <= *d1
                }),
            };
            if ok {
                bonds_out.push((i, j, d));
            }
        }
    }
    // The polyhedra: an image with four bonds or more, round its neighbours.
    let mut polyhedra = Vec::new();
    for c in 0..images.len() {
        let nb: Vec<(usize, f64)> = bonds_out.iter().filter_map(|&(i, j, d)| if i == c { Some((j, d)) } else if j == c { Some((i, d)) } else { None }).collect();
        if nb.len() < 4 {
            continue;
        }
        let pts: Vec<[f64; 3]> = nb.iter().map(|&(k, _)| images[k].cart).collect();
        let hf = hull(&pts);
        if hf.is_empty() {
            continue;
        }
        let centre = images[c].cart;
        let faces: Vec<Face> = hf.iter().map(|(v, n)| Face { vertices: v.map(|k| nb[k].0), normal: *n }).collect();
        let volume: f64 = hf.iter().map(|(v, _)| dot(sub(pts[v[0]], centre), cross(sub(pts[v[1]], centre), sub(pts[v[2]], centre))).abs() / 6.0).sum();
        let mean = nb.iter().map(|x| x.1).sum::<f64>() / nb.len() as f64;
        let distortion = nb.iter().map(|x| (x.1 - mean).abs()).sum::<f64>() / (nb.len() as f64 * mean);
        polyhedra.push(Polyhedron { centre: c, vertices: nb.iter().map(|x| x.0).collect(), faces, volume, distortion });
    }
    // The cell's edges (E22: `Cell::edges`, for any phase).
    let edges = st.cell.edges();
    // A mode phase: the reference images and each irrep's displacement on them.
    let (mut reference, mut reference_fract, mut displacements) = (Vec::new(), Vec::new(), Vec::new());
    if let Some(ms) = &model.phases[p].modes {
        let rotate = |im: &Image, u: [f64; 3]| -> [f64; 3] {
            let r = st.group.ops[im.op].r;
            std::array::from_fn(|i| (0..3).map(|j| r[i][j] as f64 * u[j]).sum())
        };
        let turn = |im: &Image, u: [f64; 3]| -> [f64; 3] { cart(rotate(im, u)) };
        // The reference image: the image's position less its atom's displacement turned by its operator.
        for im in &images {
            let u = ms.displacement(im.atom);
            reference.push(sub(im.cart, turn(im, u)));
            reference_fract.push(sub(im.fract, rotate(im, u)));
        }
        for ir in ms.irreps() {
            let u = ms.irrep_displacement(&ir.name);
            displacements.push((ir.name.clone(), images.iter().map(|im| turn(im, u[im.atom])).collect()));
        }
    }
    Ok(Scene { images, bonds: bonds_out, polyhedra, edges, reference, reference_fract, displacements })
}
