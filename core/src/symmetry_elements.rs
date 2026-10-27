//! Symmetry elements (P09, the app's request 35): the geometric reading of a Seitz operator {R|t} as International
//! Tables, Vol. A, lists it in each group's "Symmetry operations" block (ITA 2016, §1.2.2.4 and §1.4.2): its kind and
//! symbol (1, t, −1, 2, 2₁, 3₁, m, a, b, c, n, d, g, −3, −4, −6 …), its direction, its intrinsic translation (the
//! screw or glide component) and its location, the points it leaves fixed.
//!
//! For an operator of order k (Rᵏ = I) the intrinsic part is w = (1/k) Σ_{j<k} Rʲ t, the screw or glide component
//! (ITA §1.2.2.4); the element is the set of fixed points of the reduced operator {R|t − w}: the solutions p of
//! (R − I) p = −(t − w), a point plus the kernel of R − I (an axis, a plane, or a point for −1 and the
//! rotoinversions, whose axis is the proper part −R's). The sense of a 3-, 4- or 6-fold axis is the sign of
//! det[u, v, Rv] (R's proper part) for the axis u and any v off it, as in a right-handed cell. Everything is exact:
//! the core keeps t in twelfths.
//!
//! Checked operator by operator for the 530 settings, centring translations included, against cctbx's
//! `rot_mx_info` and `translation_part_info` (test `symmetry_elements`, its reference written by
//! `tests/data/make_symmetry_elements.py`).

use crate::spacegroup::SpaceGroup;
use crate::symop::{SymOp, DEN};

/// An exact rational number.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Q(i64, i64);

fn gcd(a: i64, b: i64) -> i64 {
    if b == 0 {
        a.abs().max(1)
    } else {
        gcd(b, a % b)
    }
}

impl Q {
    fn new(n: i64, d: i64) -> Q {
        let g = gcd(n, d) * d.signum();
        Q(n / g, d / g)
    }
    fn int(n: i64) -> Q {
        Q(n, 1)
    }
    fn add(self, o: Q) -> Q {
        Q::new(self.0 * o.1 + o.0 * self.1, self.1 * o.1)
    }
    fn sub(self, o: Q) -> Q {
        self.add(Q(-o.0, o.1))
    }
    fn mul(self, o: Q) -> Q {
        Q::new(self.0 * o.0, self.1 * o.1)
    }
    fn div(self, o: Q) -> Q {
        Q::new(self.0 * o.1, self.1 * o.0)
    }
    fn is_zero(self) -> bool {
        self.0 == 0
    }
    fn f(self) -> f64 {
        self.0 as f64 / self.1 as f64
    }
}

/// The kind of a symmetry element.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum ElementKind {
    Identity,
    /// A lattice (centring) translation.
    Translation,
    InversionCentre,
    RotationAxis,
    ScrewAxis,
    MirrorPlane,
    GlidePlane,
    /// A rotoinversion axis (−3, −4, −6) with its inversion point.
    RotoinversionAxis,
}

/// An operator's symmetry element.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SymmetryElement {
    pub kind: ElementKind,
    /// The symbol ITA prints: 1, t, -1, 2, 2_1, 3, 3_1, 3_2, 4, 4_1, 4_2, 4_3, 6, 6_1 … 6_5, m, a, b, c, n, d, g, -3,
    /// -4, -6.
    pub symbol: String,
    /// The rotation's order with its sign (cctbx's type): 1, 2, 3, 4, 6, negative for an improper rotation (−2 a
    /// mirror or glide).
    pub order: i32,
    /// The sense of a 3-, 4- or 6-fold rotation or rotoinversion: +1 or −1; 0 otherwise.
    pub sense: i32,
    /// The axis (for a mirror or glide, the direction R reverses; none for 1, t, −1), as integers with no common
    /// factor.
    pub direction: [i32; 3],
    /// The intrinsic translation: the screw or glide component, or the translation of t (fractions of the cell).
    pub intrinsic: [f64; 3],
    /// A point of the element: the inversion point of −1 and of a rotoinversion.
    pub point: [f64; 3],
    /// The directions the element spans: one for an axis, two for a plane, none for a point.
    pub span: Vec<[f64; 3]>,
    /// The location as ITA writes it, e.g. "x,1/4,z", "0,y,0", "x,x,x", "0,0,0"; empty for 1 and t.
    pub location: String,
    /// The intrinsic translation as ITA writes it, e.g. "(0,1/2,0)"; empty when none.
    pub translation: String,
}

fn frac_text(q: Q) -> String {
    if q.1 == 1 {
        q.0.to_string()
    } else {
        format!("{}/{}", q.0, q.1)
    }
}

/// Solves A p = b (3×3, rational) for one solution with the free coordinates at 0, and the kernel of A in reduced
/// echelon form; None when there is no solution.
fn solve(a: [[Q; 3]; 3], b: [Q; 3]) -> Option<([Q; 3], Vec<[Q; 3]>)> {
    let mut m: Vec<[Q; 4]> = (0..3).map(|i| [a[i][0], a[i][1], a[i][2], b[i]]).collect();
    let mut pivots = Vec::new();
    let mut row = 0;
    for col in 0..3 {
        let Some(p) = (row..3).find(|&r| !m[r][col].is_zero()) else { continue };
        m.swap(row, p);
        let lead = m[row][col];
        for x in m[row].iter_mut() {
            *x = x.div(lead);
        }
        for r in 0..3 {
            if r != row && !m[r][col].is_zero() {
                let (f, pivot) = (m[r][col], m[row]);
                for (x, p) in m[r].iter_mut().zip(pivot) {
                    *x = x.sub(f.mul(p));
                }
            }
        }
        pivots.push(col);
        row += 1;
    }
    if (row..3).any(|r| !m[r][3].is_zero()) {
        return None;
    }
    let mut p = [Q::int(0); 3];
    for (r, &c) in pivots.iter().enumerate() {
        p[c] = m[r][3];
    }
    let free: Vec<usize> = (0..3).filter(|c| !pivots.contains(c)).collect();
    let kernel = free
        .iter()
        .map(|&f| {
            let mut v = [Q::int(0); 3];
            v[f] = Q::int(1);
            for (r, &c) in pivots.iter().enumerate() {
                v[c] = Q::int(0).sub(m[r][f]);
            }
            v
        })
        .collect();
    Some((p, kernel))
}

/// A rational vector scaled to integers with no common factor (sign kept).
fn primitive(v: [Q; 3]) -> [i64; 3] {
    let l = v.iter().fold(1, |a, q| a / gcd(a, q.1) * q.1);
    let n: Vec<i64> = v.iter().map(|q| q.0 * (l / q.1)).collect();
    let g = n.iter().fold(0, |a, &x| gcd(a, x));
    [n[0] / g, n[1] / g, n[2] / g]
}

/// The location as ITA writes it: each coordinate as its constant plus the parameters x, y, z of the directions
/// spanned (each named by the first coordinate it moves, its coefficients made whole: "2x,x,0").
fn location_text(p: [Q; 3], span: &[[i64; 3]]) -> String {
    let names = ["x", "y", "z"];
    let mut used = Vec::new();
    let params: Vec<(&str, [i64; 3])> = span
        .iter()
        .map(|d| {
            let first = (0..3).find(|&i| d[i] != 0 && !used.contains(&i)).or_else(|| (0..3).find(|&i| d[i] != 0)).unwrap();
            used.push(first);
            // The parameter's own coordinate positive.
            let s = d[first].signum();
            (names[first], [d[0] * s, d[1] * s, d[2] * s])
        })
        .collect();
    (0..3)
        .map(|i| {
            let mut s = String::new();
            for (name, d) in &params {
                let c = d[i];
                if c == 0 {
                    continue;
                }
                if c < 0 {
                    s.push('-');
                } else if !s.is_empty() {
                    s.push('+');
                }
                if c.abs() != 1 {
                    s.push_str(&c.abs().to_string());
                }
                s.push_str(name);
            }
            let k = p[i];
            if !k.is_zero() || s.is_empty() {
                if s.is_empty() {
                    s = frac_text(k);
                } else if k.0 < 0 {
                    s.push_str(&frac_text(k));
                } else {
                    s.push('+');
                    s.push_str(&frac_text(k));
                }
            }
            s
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn mat_vec(r: &[[i64; 3]; 3], v: [Q; 3]) -> [Q; 3] {
    let mut o = [Q::int(0); 3];
    for (i, oi) in o.iter_mut().enumerate() {
        for (j, vj) in v.iter().enumerate() {
            *oi = oi.add(Q::int(r[i][j]).mul(*vj));
        }
    }
    o
}

fn det3(c: [[i64; 3]; 3]) -> i64 {
    c[0][0] * (c[1][1] * c[2][2] - c[1][2] * c[2][1]) - c[0][1] * (c[1][0] * c[2][2] - c[1][2] * c[2][0]) + c[0][2] * (c[1][0] * c[2][1] - c[1][1] * c[2][0])
}

/// The symmetry element of an operator; `lattice` holds the group's centring translations (in twelfths, (0,0,0)
/// first), which decide a screw's subscript when the shortest lattice vector along the axis is a centring one.
pub fn element_in(op: &SymOp, lattice: &[[i32; 3]]) -> SymmetryElement {
    let r: [[i64; 3]; 3] = op.r.map(|row| row.map(|x| x as i64));
    let t: [Q; 3] = op.t.map(|x| Q::new(x as i64, DEN as i64));
    let det = op.det();
    let proper: [[i64; 3]; 3] = r.map(|row| row.map(|x| x * det as i64));
    // The rotation's order as a matrix, and the type: 1, 2, 3, 4, 6, negative when improper.
    let order = {
        let mut m = r;
        let mut k = 1;
        while m != [[1, 0, 0], [0, 1, 0], [0, 0, 1]] {
            let mut n = [[0i64; 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    n[i][j] = (0..3).map(|l| m[i][l] * r[l][j]).sum();
                }
            }
            m = n;
            k += 1;
        }
        k
    };
    let trace: i64 = proper[0][0] + proper[1][1] + proper[2][2];
    let n_fold: i32 = match trace {
        3 => 1,
        -1 => 2,
        0 => 3,
        1 => 4,
        2 => 6,
        _ => 0,
    };
    let signed = if det > 0 { n_fold } else { -n_fold };
    // The intrinsic part w = (1/k) Σ Rʲ t.
    let mut sum = [Q::int(0); 3];
    let mut v = t;
    for _ in 0..order {
        for i in 0..3 {
            sum[i] = sum[i].add(v[i]);
        }
        v = mat_vec(&r, v);
    }
    let w = sum.map(|x| x.div(Q::int(order)));
    // The axis: the kernel of (proper − I); for a mirror, the direction R reverses.
    let id = |i: usize, j: usize| if i == j { 1 } else { 0 };
    let axis: Option<[i64; 3]> = if n_fold > 1 {
        let a = std::array::from_fn(|i| std::array::from_fn(|j| Q::int(proper[i][j] - id(i, j))));
        solve(a, [Q::int(0); 3]).and_then(|(_, k)| k.first().map(|v| primitive(*v)))
    } else {
        None
    };
    // The element: fixed points of {R|t − w}.
    let a: [[Q; 3]; 3] = std::array::from_fn(|i| std::array::from_fn(|j| Q::int(r[i][j] - id(i, j))));
    let rhs: [Q; 3] = std::array::from_fn(|i| Q::int(0).sub(t[i].sub(w[i])));
    let (p, kernel) = solve(a, rhs).expect("the reduced operator has fixed points");
    let span: Vec<[i64; 3]> = kernel.iter().map(|v| primitive(*v)).collect();
    // The sense, for 3-, 4- and 6-fold proper parts.
    let sense = match (n_fold >= 3, axis) {
        (true, Some(u)) => {
            let vv = if u[0] == 0 && u[1] == 0 { [1, 0, 0] } else { [0, 0, 1] };
            let rv: [i64; 3] = std::array::from_fn(|i| (0..3).map(|j| proper[i][j] * vv[j]).sum());
            det3([u, vv, rv]).signum() as i32
        }
        _ => 0,
    };
    let w_zero = w.iter().all(|x| x.is_zero());
    let wq = w;
    let translation = if w_zero { String::new() } else { format!("({})", wq.iter().map(|q| frac_text(*q)).collect::<Vec<_>>().join(",")) };
    let is_lattice = |v: [Q; 3]| -> bool {
        lattice.iter().any(|c| (0..3).all(|i| v[i].sub(Q::new(c[i] as i64, DEN as i64)).1 == 1))
    };
    let (kind, symbol) = match signed {
        1 if w_zero => (ElementKind::Identity, "1".to_string()),
        1 => (ElementKind::Translation, "t".to_string()),
        -1 => (ElementKind::InversionCentre, "-1".to_string()),
        -3 | -4 | -6 => (ElementKind::RotoinversionAxis, signed.to_string()),
        2 | 3 | 4 | 6 => {
            if w_zero {
                (ElementKind::RotationAxis, n_fold.to_string())
            } else {
                // w = α·L for the shortest lattice vector L along the axis; the subscript is α·n.
                let u = axis.expect("an axis");
                let along = |m: i64| -> [Q; 3] { u.map(|x| Q::new(x, m)) };
                let m = [6, 4, 3, 2, 1].into_iter().find(|&m| is_lattice(along(m))).unwrap_or(1);
                let l = along(m);
                let i = (0..3).find(|&i| !l[i].is_zero()).unwrap();
                let alpha = wq[i].div(l[i]);
                let k = alpha.mul(Q::int(n_fold as i64));
                let k = if k.1 == 1 { k.0.rem_euclid(n_fold as i64) } else { -1 };
                if k == 0 {
                    (ElementKind::RotationAxis, n_fold.to_string())
                } else {
                    (ElementKind::ScrewAxis, if k > 0 { format!("{n_fold}_{k}") } else { format!("{n_fold}_?") })
                }
            }
        }
        -2 => {
            if w_zero {
                (ElementKind::MirrorPlane, "m".to_string())
            } else {
                let halves: Vec<usize> = (0..3).filter(|&i| !wq[i].is_zero()).collect();
                let all = |d: i64| halves.iter().all(|&i| wq[i].1 == d && wq[i].0.abs() == 1);
                let letter = match halves.len() {
                    1 if all(2) => ["a", "b", "c"][halves[0]],
                    2 | 3 if all(2) => "n",
                    _ if halves.iter().all(|&i| wq[i].1 == 4) => "d",
                    _ => "g",
                };
                (ElementKind::GlidePlane, letter.to_string())
            }
        }
        _ => (ElementKind::Identity, "?".to_string()),
    };
    let location = match kind {
        ElementKind::Identity | ElementKind::Translation => String::new(),
        ElementKind::RotoinversionAxis => format!("{}; {}", location_text(p, &[axis.unwrap()]), location_text(p, &[])),
        _ => location_text(p, &span),
    };
    SymmetryElement {
        kind,
        symbol,
        order: signed,
        sense,
        direction: axis.map_or([0, 0, 0], |u| u.map(|x| x as i32)),
        intrinsic: w.map(|x| x.f()),
        point: p.map(|x| x.f()),
        span: if matches!(kind, ElementKind::Identity | ElementKind::Translation) { Vec::new() } else { span.iter().map(|d| d.map(|x| x as f64)).collect() },
        location,
        translation,
    }
}

/// The symmetry element of an operator in a primitive lattice.
pub fn element(op: &SymOp) -> SymmetryElement {
    element_in(op, &[[0, 0, 0]])
}

impl SpaceGroup {
    /// Each operator with its symmetry element, centring translations included, in the order of `ops`.
    pub fn symmetry_elements(&self) -> Vec<(SymOp, SymmetryElement)> {
        self.ops.iter().map(|o| (*o, element_in(o, &self.centring))).collect()
    }
}
