//! The Explore stream's requests E18–E21 (explore/REQUESTS.md; L376), on the Y₂O₃ tutorial refined from its `.pcr`:
//! - `e18a_the_cell_with_its_sus`: the state's a, b, c, α, β, γ and V with su's through the covariance; a's su is the
//!   refined cell parameter's esd, V's 3a²σ(a) for the cubic cell (b and c tied to a), the angles 90° with none;
//!   without a covariance, none;
//! - `e18b_e18c_plausible_ranges`: each Y–O contact's plausible distances at Y's coordination 6 (`state::plausible`'s),
//!   none for Y–Y or O–O; each valence sum's |q| and its plausible range |q|(1 ∓ 0.10), +1.97 for O against 2 and not
//!   −2; the phase's `valence_sum_range` (0.05) through the control file;
//! - `e20a_a_point`: a state at given values of the free parameters, at a grid node equal to `chi2_grid`'s χ² there
//!   (and the grid's N − P + C its own; its 1σ contour one closed line);
//!   a parameter that moves two fields with coefficients 1 and ½, a variable, and a value linked to another by ½, set
//!   as a step sets them; a name that is no free parameter refused;
//! - `e21_contours`: a quadratic χ² on a 41 × 41 grid whose contours are ellipses: every point of the 1σ and 2σ
//!   contours (Δχ² = k²χ²_min/dof) on its ellipse within linear interpolation's bound h²Q''/8, each a closed line; a node not computed breaks the line there.

use drs::jacobian::Field;
use drs::lsq::{refine, Options};
use drs::powder::{Intensities, Model};
use drs::refinable::ParameterSet;
use std::path::{Path, PathBuf};

fn dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs/Y2O3")
}

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display())).iter().map(|&b| b as char).collect()
}

/// Y₂O₃ from the tutorial's `.pcr`, its charges set (Y³⁺, O²⁻), refined 4 cycles.
fn y2o3() -> (Model, drs::data::Pattern, ParameterSet, drs::lsq::Refinement) {
    let p = drs::pcr::read(&latin1(&dir().join("y2o3_structure.pcr"))).unwrap();
    let l = Model::from_pcr(&p, &dir(), "y2o3_structure").unwrap();
    let (mut m, data) = (l.model, l.data.unwrap());
    if let Intensities::Structure(st, _) = &mut m.phases[0].intensities {
        for a in &mut st.atoms {
            a.charge = Some(if a.species.starts_with('Y') { 3 } else { -2 });
        }
    }
    let set = ParameterSet::from_pcr(&p, &m).unwrap();
    let r = refine(&mut m, &data, &set, Options { max_cycles: 4, eps: 0.01 }, |_| {});
    (m, data, set, r)
}

#[test]
fn e18a_the_cell_with_its_sus() {
    let (m, data, set, r) = y2o3();
    let st = drs::state::summary(&m, &data, Some(&set), Some(&r), &drs::state::Options::default());
    let c = st.phases[0].cell_su;
    let esd = r.params.iter().find(|p| p.name == "phase1.cell_a").unwrap().esd;
    println!("a {:.6}({:.1e}), V {:.4}({:.2e}); the cell parameter's esd {esd:.1e}", c[0].0, c[0].1, c[6].0, c[6].1);
    // The propagator differentiates by central differences: to 10⁻⁶ of the su.
    assert!(esd > 0.0 && (c[0].1 - esd).abs() < 1e-6 * esd, "a's su {} against the esd {esd}", c[0].1);
    assert!((c[1].1 - esd).abs() < 1e-6 * esd && (c[2].1 - esd).abs() < 1e-6 * esd, "b and c tied to a");
    let a = c[0].0;
    assert!((c[6].0 - a * a * a).abs() < 1e-9 * a * a * a);
    assert!((c[6].1 - 3.0 * a * a * esd).abs() < 1e-6 * 3.0 * a * a * esd, "V's su {} against 3a²σ(a) {}", c[6].1, 3.0 * a * a * esd);
    for (v, s) in &c[3..6] {
        assert_eq!((*v, *s), (90.0, 0.0));
    }
    let none = drs::state::summary(&m, &data, None, None, &drs::state::Options::default());
    assert!(none.phases[0].cell_su.iter().all(|x| x.1 == 0.0));
}

#[test]
fn e18b_e18c_plausible_ranges() {
    let (m, data, set, r) = y2o3();
    let st = drs::state::summary(&m, &data, Some(&set), Some(&r), &drs::state::Options::default());
    let ph = &st.phases[0];
    let want = drs::state::plausible(("Y", 3), ("O", -2), &[6]);
    let species = |i: usize| match &m.phases[0].intensities {
        Intensities::Structure(s, _) => s.atoms[i].species.clone(),
        _ => unreachable!(),
    };
    let (mut yo, mut others) = (0, 0);
    for (i, list) in ph.contacts.iter().enumerate() {
        for c in list {
            if species(i).starts_with('Y') != species(c.atom).starts_with('Y') {
                assert_eq!(c.plausible.as_ref(), Some(&want), "{} – {}", ph.labels[i], ph.labels[c.atom]);
                yo += 1;
            } else {
                assert!(c.plausible.is_none());
                others += 1;
            }
        }
    }
    println!("{yo} Y–O contacts with the plausible distances {want:?}, {others} others without");
    assert!(yo > 0 && others > 0);
    let (sums, _) = ph.valence.as_ref().unwrap();
    for s in sums {
        let q = if species(s.atom).starts_with('Y') { 3.0 } else { 2.0 };
        assert_eq!(s.expected, q);
        assert!((s.plausible.0 - 0.9 * q).abs() < 1e-12 && (s.plausible.1 - 1.1 * q).abs() < 1e-12);
    }
    let o = sums.iter().find(|s| species(s.atom).starts_with('O')).unwrap();
    println!("O: sum {:.3} against {} within {:?}", o.sum.0, o.expected, o.plausible);
    assert!(o.valence == -2.0 && o.sum.0 > 0.0 && o.plausible.0 < o.sum.0 && o.sum.0 < o.plausible.1);
    // The phase's own range, through the control file.
    let p = drs::pcr::read(&latin1(&dir().join("y2o3_structure.pcr"))).unwrap();
    let mut c = drs::control::Control::from_pcr(&p, &dir(), "y2o3_structure", &dir()).unwrap();
    c.model.phases[0].valence_sum_range = Some(0.05);
    let text = c.to_text();
    assert!(text.contains(&format!("{} = 0.05", drs::bonds::VALENCE_RANGE_KEY)), "the key written");
    let back = drs::control::Control::from_text_in(&text, &dir()).unwrap();
    assert_eq!(back.model.phases[0].valence_sum_range, Some(0.05));
    let mut m2 = m.clone();
    m2.phases[0].valence_sum_range = Some(0.05);
    let st2 = drs::state::summary(&m2, &data, Some(&set), Some(&r), &drs::state::Options::default());
    let o2 = st2.phases[0].valence.as_ref().unwrap().0.iter().find(|s| s.atom == o.atom).unwrap();
    assert!((o2.plausible.0 - 1.9).abs() < 1e-12 && (o2.plausible.1 - 2.1).abs() < 1e-12, "{:?}", o2.plausible);
}

#[test]
fn e20a_a_point() {
    let (m, data, set, r) = y2o3();
    let a = set.free.iter().position(|p| p.name == "phase1.cell_a").unwrap();
    let b = set.free.iter().position(|p| p.name == "phase1.scale_factor").unwrap();
    let (sa, sb) = (r.params[a].esd, r.params[b].esd);
    let g = drs::lsq::chi2_grid(&m, &data, &set, (a, b), (sa, sb), 5, |_, _| true);
    let (na, nb) = (set.free[a].name.as_str(), set.free[b].name.as_str());
    for (i, j) in [(0, 0), (1, 3), (4, 2)] {
        let p = drs::lsq::point(&m, &data, &set, &[(na, g.xa[i]), (nb, g.xb[j])]).unwrap();
        assert_eq!(p.statistics.chi2, g.chi2[i][j], "node ({i}, {j})");
        assert_eq!(g.dof, p.statistics.n + p.statistics.c - p.statistics.p, "the grid's N − P + C");
        assert_eq!(p.model.get(set.free[a].fields[0].0), g.xa[i]);
        assert_eq!(p.difference.len(), data.x.len());
        assert!((p.difference[10] - (data.y[10] - p.ycalc[10])).abs() < 1e-9);
    }
    // E21 on the grid itself: ± 3σ on 11 × 11, the 1σ contour one closed line inside it.
    let g11 = drs::lsq::chi2_grid(&m, &data, &set, (a, b), (3.0 * sa, 3.0 * sb), 11, |_, _| true);
    let c1 = &g11.contours_at_sigma(&[1.0])[0];
    println!("Y2O3's 1σ contour: {} points, a over {:.2}σ", c1.x.len(), (c1.x.iter().cloned().fold(f64::MIN, f64::max) - c1.x.iter().cloned().fold(f64::MAX, f64::min)) / sa);
    assert!(c1.x.len() > 8 && c1.x.iter().all(|v| v.is_finite()) && c1.x[0] == *c1.x.last().unwrap(), "one closed line");
    // A parameter moving two fields, the second by ½; a variable moving one by 2.
    let (fy, fo) = (Field::AtomB(0, 0), Field::AtomB(0, 2));
    let mut s2 = set.clone();
    s2.free.retain(|p| !p.fields.iter().any(|f| f.0 == fy || f.0 == fo));
    s2.free.push(drs::refinable::Free { name: "pair".into(), fields: vec![(fy, 1.0), (fo, 0.5)], damping: 1.0, value: None });
    s2.free.push(drs::refinable::Free { name: "v".into(), fields: vec![(Field::AtomB(0, 1), 2.0)], damping: 1.0, value: Some(7.0) });
    let p = drs::lsq::point(&m, &data, &s2, &[("pair", m.get(fy) + 0.2), ("v", 7.1)]).unwrap();
    assert!((p.model.get(fy) - m.get(fy) - 0.2).abs() < 1e-12 && (p.model.get(fo) - m.get(fo) - 0.1).abs() < 1e-12, "the pair by 1 and ½");
    assert!((p.model.get(Field::AtomB(0, 1)) - m.get(Field::AtomB(0, 1)) - 0.2).abs() < 1e-12, "the variable's field by 2 × 0.1");
    assert!(drs::lsq::point(&m, &data, &s2, &[("nothing", 1.0)]).is_err());
    // A value linked to another with coefficient ½ (x = ½·p): it follows the parameter set.
    let mut s3 = set.clone();
    s3.free.push(drs::refinable::Free { name: "y".into(), fields: vec![(fy, 1.0)], damping: 1.0, value: None });
    s3.links.push(drs::steering::Link { value: fo, follows: fy, coefficient: 0.5, offset: 0.0 });
    let p = drs::lsq::point(&m, &data, &s3, &[("y", 0.5)]).unwrap();
    assert_eq!((p.model.get(fy), p.model.get(fo)), (0.5, 0.25), "the link enforced");
}

#[test]
fn e21_contours() {
    // χ² = 10 + Q(x, y), Q = u C⁻¹ u with σa = 0.5, σb = 2, r = 0.6 about (1, 3): the contour Q = k² is the ellipse
    // of k σ.
    let (ma, mb, sa, sb, r) = (1.0, 3.0, 0.5, 2.0, 0.6);
    let q = |x: f64, y: f64| {
        let (u, v) = ((x - ma) / sa, (y - mb) / sb);
        (u * u - 2.0 * r * u * v + v * v) / (1.0 - r * r)
    };
    let n = 41;
    let xa: Vec<f64> = (0..n).map(|k| ma - 3.0 * sa + 6.0 * sa * k as f64 / (n - 1) as f64).collect();
    let xb: Vec<f64> = (0..n).map(|k| mb - 3.0 * sb + 6.0 * sb * k as f64 / (n - 1) as f64).collect();
    let chi2: Vec<Vec<f64>> = xa.iter().map(|&x| xb.iter().map(|&y| 10.0 + q(x, y)).collect()).collect();
    // χ²_min/dof = 1: the kσ contour at χ²_min + k².
    let mut g = drs::lsq::Grid { a: "a".into(), b: "b".into(), xa, xb, chi2, order: Vec::new(), computed: n * n, dof: 10 };
    let cs = g.contours_at_sigma(&[1.0, 2.0]);
    let min = g.minimum().unwrap().2;
    // Linear interpolation along an edge of h (in σ, 6/40) misses the quadratic by at most h²Q''/8, Q'' = 2/(1 − r²)
    // along either axis.
    let bound = (6.0f64 / 40.0).powi(2) * 2.0 / (1.0 - r * r) / 8.0;
    for (c, k2) in cs.iter().zip([1.0, 4.0]) {
        assert!((c.level - min - k2).abs() < 1e-12);
        assert!(c.x.iter().all(|v| v.is_finite()), "one piece");
        let worst = c.x.iter().zip(&c.y).map(|(&x, &y)| (q(x, y) + 10.0 - min - k2).abs()).fold(0.0, f64::max);
        println!("Δχ² = {k2}: {} points, worst |ΔQ| {worst:.2e} (bound {bound:.2e})", c.x.len());
        assert!(c.x.len() > 40 && worst <= bound * (1.0 + 1e-9), "on the ellipse: {worst}");
        assert_eq!((c.x[0], c.y[0]), (*c.x.last().unwrap(), *c.y.last().unwrap()), "closed");
    }
    // A node not computed, on the Δχ² = 4 line: broken there.
    let i = 20;
    let j = (0..n).min_by(|&p, &s| ((g.chi2[i][p] - min - 4.0).abs()).partial_cmp(&(g.chi2[i][s] - min - 4.0).abs()).unwrap()).unwrap();
    g.chi2[i][j] = f64::NAN;
    let c = &g.contours(&[min + 4.0])[0];
    assert!(c.x.iter().any(|v| v.is_nan()) || c.x.first() != c.x.last(), "broken where the node is missing");
}
