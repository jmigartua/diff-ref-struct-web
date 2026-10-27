//! Every editable value of a model, and its edits (test `model_editing`, the app's requests 2 and 3, L058).
//!
//! - `Model::fields` on every FullProf run the core imports: each address unique, `refinable::field_of` its inverse,
//!   the value `Model::get`'s, and what fixes a value: the cubic cell's b, c and angles, the coordinates of atoms on
//!   special positions.
//! - `Model::set_field`: the values tied to the one set follow (a cubic cell's b and c, MgAl2O4's O y and z = x), and a
//!   value that follows from others is refused.
//! - `ParameterSet::from_addresses` reproduces `from_pcr`: on five runs the set built from the addresses of the
//!   `.pcr`'s parameters and their equalities refines to the same χ² and values (10⁻⁶), and its errors name the
//!   address; a user's limit is held.
//! - The structural edits: a phase added from a setting and a cell, and from a CIF, with atoms, calculated; removed,
//!   the pattern is the original's again; background, excluded regions, profile, asymmetry and preferred orientation,
//!   with the edits that must be refused.

use drs::edit::Asymmetry;
use drs::jacobian::Field;
use drs::lsq::{refine, Options};
use drs::powder::{BackgroundModel, Model, Shape};
use drs::refinable::{address, field_of, ParameterSet};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

fn root(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof").join(rel)
}

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display())).iter().map(|&b| b as char).collect()
}

fn load(run: &str, code: &str) -> (drs::pcr::Pcr, drs::powder::Loaded) {
    let d = root("runs").join(run);
    let mut p = drs::pcr::read(&latin1(&d.join(format!("{code}.pcr")))).unwrap();
    p.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
    let l = Model::from_pcr(&p, &d, code).unwrap_or_else(|e| panic!("{run}: {e}"));
    (p, l)
}

const RUNS: &[(&str, &str)] = &[
    ("Y2O3", "y2o3_structure"),
    ("Co-YBCO", "Cocus-asym"),
    ("MgAl2O4-s", "MgAl2O4s_Rietveld"),
    ("rutile-anatase", "rutana_tch"),
    ("Egypt", "Egypt_qa_Rietveld"),
    ("YVO3", "YVO3a"),
    ("Yb2ZrWO8", "ybfluo_Riet_IRF"),
    ("LaB6", "LaB6rx"),
    ("NdNi-HTO", "NdNi_HTO"),
    ("Sr-ox-test", "sr-ox-test_conv"),
    ("K2NaAgCN", "k2na1p5_conv"),
];

#[test]
fn fields_and_addresses() {
    let mut total = 0;
    for &(run, code) in RUNS {
        let (_, l) = load(run, code);
        let m = &l.model;
        let fields = m.fields();
        let addrs: BTreeSet<&str> = fields.iter().map(|f| f.address.as_str()).collect();
        assert_eq!(addrs.len(), fields.len(), "{run}: an address given twice");
        for fi in &fields {
            assert_eq!(field_of(&fi.address, m).unwrap(), fi.field, "{run}: {}", fi.address);
            assert_eq!(fi.value, m.get(fi.field), "{run}: {}", fi.address);
            assert!(drs::parameters::lookup(fi.id).is_some(), "{run}: {} has id {}", fi.address, fi.id);
        }
        let fixed = fields.iter().filter(|f| f.fixed_by.is_some()).count();
        println!("{run:15} {} values, {fixed} fixed or tied by the symmetry", fields.len());
        total += fields.len();
    }
    // Y2O3 (I a −3): the cell's b, c and angles follow from a; Y1 at 8b (1/4, 1/4, 1/4) has no free coordinate.
    let (_, l) = load("Y2O3", "y2o3_structure");
    let m = &l.model;
    let by = |a: &str| m.fields().into_iter().find(|f| f.address == a).unwrap_or_else(|| panic!("{a}")).fixed_by;
    assert_eq!(by("phase1.cell_a"), None);
    for a in ["phase1.cell_b", "phase1.cell_c", "phase1.cell_alpha", "phase1.cell_gamma"] {
        assert_eq!(by(a), Some("crystal system"), "{a}");
    }
    let err = field_of("phase1.Q9.atom_biso", m).unwrap_err();
    assert!(err.contains("phase1.Q9.atom_biso"), "{err}");
    println!("{total} values in {} runs", RUNS.len());
}

#[test]
fn set_field_follows_the_symmetry() {
    // Cubic: a sets b and c; b is refused.
    let (_, l) = load("Y2O3", "y2o3_structure");
    let mut m = l.model.clone();
    m.set_field(Field::Cell(0, 0), 10.61).unwrap();
    assert_eq!((m.get(Field::Cell(0, 1)), m.get(Field::Cell(0, 2))), (10.61, 10.61));
    assert!(m.set_field(Field::Cell(0, 1), 10.7).unwrap_err().contains("crystal system"));
    assert!(m.set_field(Field::Zero, f64::NAN).is_err());
    assert!(m.set_field(Field::AtomB(0, 99), 1.0).is_err());
    // MgAl2O4 (F d −3 m): O at 32e (x, x, x) moves as one.
    let (_, l) = load("MgAl2O4-s", "MgAl2O4s_Rietveld");
    let mut m = l.model.clone();
    let o = m.fields().into_iter().find(|f| f.address == "phase1.O.atom_position.x").unwrap().field;
    let Field::AtomX(p, a, _) = o else { panic!() };
    m.set_field(o, 0.2620).unwrap();
    assert_eq!((m.get(Field::AtomX(p, a, 1)), m.get(Field::AtomX(p, a, 2))), (0.2620, 0.2620));
    assert!(m.set_field(Field::AtomX(p, a, 1), 0.27).unwrap_err().contains("site symmetry"));
}

#[test]
fn from_addresses_reproduces_from_pcr() {
    for &(run, code) in &[("Y2O3", "y2o3_structure"), ("MgAl2O4-s", "MgAl2O4s_Rietveld"), ("rutile-anatase", "rutana_tch"), ("Egypt", "Egypt_qa_Rietveld"), ("Sr-ox-test", "sr-ox-test_conv")] {
        let (p, l) = load(run, code);
        let m0 = l.model.clone();
        let data = l.data.as_ref().unwrap();
        let set = ParameterSet::from_pcr(&p, &m0).unwrap();
        // The same choices as addresses: each parameter's first field free, its other fields equal to it (those the
        // symmetry ties follow by themselves).
        let free: Vec<String> = set.free.iter().map(|q| q.name.clone()).collect();
        let mut equal = Vec::new();
        for q in &set.free {
            let (f0, m_0) = q.fields[0];
            for &(f, mk) in &q.fields[1..] {
                if m0.fixed_by(f).is_none() {
                    equal.push((address(f, &m0), mk / m_0, address(f0, &m0)));
                }
            }
        }
        let set2 = ParameterSet::from_addresses(&free, &equal, &[], &m0).unwrap_or_else(|e| panic!("{run}: {e}"));
        assert_eq!(set2.free.len(), set.free.len(), "{run}");
        let opts = Options { max_cycles: 30, eps: p.eps };
        let (mut a, mut b) = (m0.clone(), m0.clone());
        let (ra, rb) = (refine(&mut a, data, &set, opts, |_| {}), refine(&mut b, data, &set2, opts, |_| {}));
        let worst = a.all_fields().into_iter().map(|f| (a.get(f) - b.get(f)).abs() / a.get(f).abs().max(1e-3)).fold(0.0, f64::max);
        println!("{run:15} {} parameters: χ² {:.6} (from_pcr) {:.6} (from_addresses); values differ by {worst:.1e}", set.free.len(), ra.stats.chi2, rb.stats.chi2);
        assert!((ra.stats.chi2 - rb.stats.chi2).abs() < 1e-6 * ra.stats.chi2, "{run}");
        assert!(worst < 1e-6, "{run}: values differ by {worst}");
    }
    // Errors, and a user's limit.
    let (_, l) = load("Y2O3", "y2o3_structure");
    let m = &l.model;
    let s = |v: &[&str]| v.iter().map(|x| x.to_string()).collect::<Vec<_>>();
    assert!(ParameterSet::from_addresses(&s(&["phase1.cell_x"]), &[], &[], m).unwrap_err().contains("phase1.cell_x"));
    assert!(ParameterSet::from_addresses(&s(&["phase1.cell_b"]), &[], &[], m).unwrap_err().contains("crystal system"));
    assert!(ParameterSet::from_addresses(&s(&["phase1.caglioti_u"]), &[("phase1.caglioti_v".into(), 1.0, "phase1.caglioti_w".into())], &[], m).unwrap_err().contains("not free"));
    // W starts above a lower limit set above its best value (w0): the refinement takes it down to the limit.
    let w0 = m.phases[0].w;
    let lo = w0 + 0.0003;
    let set = ParameterSet::from_addresses(&s(&["phase1.scale_factor", "phase1.caglioti_w"]), &[], &[("phase1.caglioti_w".into(), Some(lo), None)], m).unwrap();
    let mut mm = m.clone();
    mm.phases[0].w = w0 + 0.0008;
    let r = refine(&mut mm, l.data.as_ref().unwrap(), &set, Options { max_cycles: 20, eps: 0.1 }, |_| {});
    println!("Y2O3 with W ≥ {lo:.6}: W = {:.6} ({}; at its limit: {})", mm.phases[0].w, r.stop, r.params[1].at_limit);
    assert!((mm.phases[0].w - lo).abs() < 1e-9 && r.params[1].at_limit, "the user's limit is not held: {}", mm.phases[0].w);
}

#[test]
fn structural_edits() {
    let (_, l) = load("Y2O3", "y2o3_structure");
    let data = l.data.as_ref().unwrap();
    let m0 = l.model.clone();
    let y0 = m0.calculate(&data.x).ycalc;
    let mut m = m0.clone();
    // A spinel phase from a setting and a cell, with three atoms.
    let (g, _) = drs::spacegroup::SpaceGroup::from_symbol("F d -3 m:2").unwrap();
    let cell = drs::cell::Cell::new(8.083, 8.083, 8.083, 90.0, 90.0, 90.0).unwrap();
    let p = m.add_phase("spinel", g, cell).unwrap();
    m.add_atom(p, "Mg", "Mg", [0.125, 0.125, 0.125], 0.3, 0.0417).unwrap();
    m.add_atom(p, "Al", "Al", [0.5, 0.5, 0.5], 0.3, 0.0833).unwrap();
    m.add_atom(p, "O", "O", [0.262, 0.262, 0.262], 0.5, 0.1667).unwrap();
    assert!(m.add_atom(p, "O", "O", [0.1, 0.2, 0.3], 0.5, 1.0).is_err(), "a label given twice");
    m.phases[p].scale = 1e-4;
    let c = m.calculate(&data.x);
    let n = c.peaks.iter().filter(|k| k.phase == p).count();
    let added: f64 = c.ycalc.iter().zip(&y0).map(|(a, b)| a - b).sum();
    println!("spinel added: {n} peaks, {added:.1} counts");
    assert!(n > 0 && added > 0.0);
    // The new phase's values are fields too.
    assert!(m.fields().iter().any(|f| f.address == "phase2.O.atom_position.x" && f.fixed_by.is_none()));
    assert!(m.fields().iter().any(|f| f.address == "phase2.O.atom_position.y" && f.fixed_by == Some("site symmetry")));
    m.remove_atom(p, 1).unwrap();
    m.remove_phase(p).unwrap();
    assert_eq!(m.calculate(&data.x).ycalc, y0, "removed, the pattern is the original's");
    assert!(m.remove_phase(5).is_err());
    // A phase from a CIF.
    let cif = latin1(&root("7-FullProf-Tutorial-SymmetryModes/LaMnO3/SimulatedAnnealing/Lmo_san.cif"));
    let q = m.add_phase_from_cif("LaMnO3", &cif).unwrap();
    println!("from the CIF: {} atoms", match &m.phases[q].intensities {
        drs::powder::Intensities::Structure(st, _) => st.atoms.len(),
        _ => 0,
    });
    m.remove_phase(q).unwrap();
    // Background, excluded regions, profile, asymmetry, preferred orientation.
    assert!(m.set_background(BackgroundModel::Points(vec![(20.0, 1.0), (10.0, 2.0)])).is_err());
    assert!(m.set_background_terms(3).is_err(), "Y2O3's background is by points");
    m.set_background(BackgroundModel::Polynomial { origin: 40.0, coefficients: vec![10.0, 1.0] }).unwrap();
    m.set_background_terms(6).unwrap();
    assert_eq!(m.all_fields().iter().filter(|f| matches!(f, Field::Background(_))).count(), 6);
    let ne = m.excluded.len();
    m.add_excluded(100.0, 101.0).unwrap();
    assert!(m.add_excluded(3.0, 2.0).is_err());
    m.remove_excluded(m.excluded.iter().position(|e| e.0 == 100.0).unwrap()).unwrap();
    assert_eq!(m.excluded.len(), ne);
    assert!(m.set_asymmetry(0, Asymmetry::AxialDivergence { s_l: 0.02, h_l: 0.02 }).is_err(), "Y2O3 is pseudo-Voigt");
    m.set_shape(0, Shape::ThompsonCoxHastings).unwrap();
    m.set_asymmetry(0, Asymmetry::AxialDivergence { s_l: 0.02, h_l: 0.02 }).unwrap();
    assert!(m.fields().iter().any(|f| f.address == "phase1.axial_divergence_s_l"));
    m.set_asymmetry(0, Asymmetry::None).unwrap();
    assert!(m.set_preferred_orientation(0, [0.0, 0.0, 0.0], 1.0, 0.0).is_err());
    m.set_preferred_orientation_function(true);
    m.set_preferred_orientation(0, [0.0, 0.0, 1.0], 0.9, 0.1).unwrap();
    assert!(m.set_preferred_orientation(0, [0.0, 0.0, 1.0], -0.5, 0.0).is_err());
    m.calculate(&data.x);
}
