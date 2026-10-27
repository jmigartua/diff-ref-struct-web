//! Sequential and parametric refinement (P54, L317) and a series' outputs (P59), against FullProf's own sequential
//! refinement of LaMnO₃ on D20 (`manuals/fullprof/1-WinPLOTR/2_FP_sequential_refinements/LaMnO3_D20`: 71 patterns
//! from 68.9 to 168.1 K refined in turn, each from the values the one before ended at, `la0.log`; each refined state
//! in `la0_NNN.pcr`, the values and σ's in `la0.seq`).
//!
//! - **The chain** (`fullprof_chain`): the 16 patterns above the Néel point (la0175–la0190, 146.1–168.1 K), from
//!   FullProf's refined la0_174.pcr. The run has a magnetic phase (Jbt = 1) the core does not model; above the Néel
//!   point FullProf refines its moment to 0.3 μB, its sign changing from one pattern to the next (noise), so the core
//!   leaves it out there. Checked per pattern: the temperature from the data file's header is FullProf's PARVAL; every
//!   refined value is FullProf's within its σ, but one (frozen below); the esd's are FullProf's to 2.6 %, but the
//!   cell's (twice FullProf's: FullProf halves the cell's esd's, L053, L113) and the two end background points' (frozen);
//!   the cell volume is FullProf's within FullProf's σ; χ² is below FullProf's, by 0.27–0.71 % (frozen; 0.22 % of it
//!   at FullProf's own values, la0175: the 0.3 μB moment the core leaves out). la0186 is compared on its values only:
//!   FullProf's fit of it ended at χ² 4.32 after its magnetic moment's shift jumped to 55σ in cycle 16 (la0.log), where
//!   cycles 4 to 15 held 3.25, the core's χ².
//!   Frozen for Josu: the last background point (146.92°, flat beyond it to the data's end at 150.3°) is the one value
//!   beyond FullProf's σ (up to 1.5σ) and its √C 1/√2 of FullProf's to 1 % at every pattern; the core's would be
//!   FullProf's with the data beyond about 148.0° left out of its derivative (scanned on la0175), which no FullProf
//!   setting names. Not settled. The esd's are compared as √C, each over the root of its own χ².
//! - **The series' outputs** (`series_outputs`): the CSV and the multi-block CIF of the chain, every point read back;
//!   each point's `_cell_volume` with the table's su, the covariance's (Publish's P43).
//! - **Parametric fits** (`parametric_fits`): FullProf's own volumes and c of the 71 patterns fitted to three
//!   equations of T, against SciPy's curve_fit (`tests/data/make_parametric_reference.py`).

use drs::jacobian::Field;
use drs::sequential::{run, Point, Start};
use std::path::{Path, PathBuf};

mod common;
use common::latin1;

fn dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/1-WinPLOTR/2_FP_sequential_refinements/LaMnO3_D20")
}

/// FullProf's la0_174.pcr, its magnetic phase left out, as a control file; the 16 patterns above the Néel point.
fn chain(start: Start) -> (drs::control::Control, Vec<Point>) {
    let d = dir();
    let mut p = drs::pcr::read(&latin1(&d.join("la0_174.pcr"))).unwrap();
    p.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
    let c = drs::control::Control::from_pcr(&p, &d, "la0_174", &d).unwrap();
    let inputs: Vec<_> = (175..=190).map(|k| drs::sequential::read(&c, &d, Path::new(&format!("la0{k}.dat"))).unwrap()).collect();
    let points = run(&c, &inputs, start, |_, _| true).unwrap();
    (c, points)
}

#[test]
fn fullprof_chain() {
    let seq = drs::seq::read(&latin1(&dir().join("la0.seq"))).unwrap();
    let (c, points) = chain(Start::Previous);
    assert_eq!(points.len(), 16);
    let labels: Vec<Vec<String>> = vec![match &c.model.phases[0].intensities {
        drs::powder::Intensities::Structure(st, _) => st.atoms.iter().map(|a| a.label.clone()).collect(),
        _ => Vec::new(),
    }];
    let mut fails = Vec::new();
    let (mut lo, mut hi) = (f64::MAX, f64::MIN);
    let mut ranges = [(f64::MAX, f64::MIN); 4];
    for pt in &points {
        let b = seq.blocks.iter().rev().find(|b| b.name == pt.label).unwrap();
        // The temperature: the data file's header, FullProf's PARVAL.
        if pt.parameter != b.parameter {
            fails.push(format!("{}: T {:?} against PARVAL {:?}", pt.label, pt.parameter, b.parameter));
        }
        let mut compared = 0;
        for q in &b.parameters {
            let Some(f) = common::fullprof_field(&q.name, &labels) else { continue };
            let Some(j) = pt.set.free.iter().position(|fp| fp.fields[0].0 == f) else {
                // The magnetic moment, in the phase left out.
                assert!(q.name.ends_with("_ph2"), "{}: {} not refined", pt.label, q.name);
                continue;
            };
            compared += 1;
            let r = &pt.refinement.params[j];
            let dev = (r.value - q.value).abs() / q.sigma;
            let last_point = f == Field::Background(29);
            if dev > if last_point { 1.5 } else { 1.0 } {
                fails.push(format!("{}: {} {} against FullProf's {} ± {}: {dev:.2}σ", pt.label, q.name, r.value, q.value, q.sigma));
            }
            // The esd's: each class's ratio to FullProf's, its range over the chain frozen below.
            let class = match f {
                Field::Cell(..) => 0,
                Field::Background(29) => 1,
                Field::Background(2) => 2,
                _ => 3,
            };
            // √C, each esd over the root of its own χ² (FullProf's esd's are √(C χ²) too, L053): la0186's χ² 4.33
            // would otherwise raise FullProf's by 15 %.
            let ratio = r.esd / q.sigma * (b.chi2.unwrap() / pt.refinement.stats.chi2).sqrt();
            ranges[class].0 = ranges[class].0.min(ratio);
            ranges[class].1 = ranges[class].1.max(ratio);
        }
        assert_eq!(compared, 45, "{}", pt.label);
        // The cell volume.
        let (v, _) = drs::sequential::cell_volume(pt, 0);
        let (fv, fs) = b.phases[0].volume.unwrap();
        if (v - fv).abs() > fs {
            fails.push(format!("{}: V {v:.5} against FullProf's {fv} ± {fs}", pt.label));
        }
        // χ², but la0186's (FullProf's fit failed in its magnetic phase).
        let ratio = pt.refinement.stats.chi2 / b.chi2.unwrap();
        println!("{} T {:7.3}: {} cycles, χ² {:.4} against FullProf's {:.4} ({:+.2} %); V {v:.4} against {fv}", pt.label, pt.parameter.unwrap(), pt.refinement.cycles, pt.refinement.stats.chi2, b.chi2.unwrap(), 100.0 * (ratio - 1.0));
        if pt.label != "la0186" {
            lo = lo.min(ratio);
            hi = hi.max(ratio);
        }
    }
    // √C against FullProf's over the chain: the cell's twice FullProf's (FullProf halves them, L053, L113); the last
    // background point's 0.707–0.713, 1/√2 to 1 % at every pattern (its derivative, for Josu); the second point's
    // (12.40°, next to the region excluded below 10°) 0.95–0.96; every other value's FullProf's to 2.6 %.
    println!("esd ratios: cell {:?}, last point {:?}, point 2 {:?}, the rest {:?}", ranges[0], ranges[1], ranges[2], ranges[3]);
    let frozen = [(1.966, 2.016), (0.707, 0.714), (0.945, 0.960), (0.974, 1.019)];
    for (k, (r, f)) in ranges.iter().zip(&frozen).enumerate() {
        if r.0 < f.0 || r.1 > f.1 {
            fails.push(format!("esd ratio class {k} from {:.4} to {:.4}, frozen {:?}", r.0, r.1, f));
        }
    }
    if lo < 0.9925 || hi > 0.9975 {
        fails.push(format!("χ² against FullProf's from {lo:.5} to {hi:.5}, frozen 0.9925–0.9975"));
    }
    // Each pattern starts where the one before ended: la0176's first cycle starts at la0175's model on la0176.
    let d176 = &points[1].data;
    let start = drs::lsq::statistics(&points[0].model, d176, &points[0].model.calculate(&d176.x).ycalc, points[1].set.free.len());
    assert!((start.chi2 - points[1].refinement.history[0].chi2_start).abs() < 1e-9 * start.chi2, "{} {}", start.chi2, points[1].refinement.history[0].chi2_start);
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}

#[test]
fn each_from_the_control_file() {
    // From the control file's values each time: la0175 the same as in the chain, and a later one from la0_174's values,
    // ending where the chain ends to the convergence's precision.
    let (_, each) = chain(Start::Each);
    let (_, prev) = chain(Start::Previous);
    assert_eq!(each[0].refinement.stats.chi2, prev[0].refinement.stats.chi2);
    for (a, b) in each.iter().zip(&prev).skip(1) {
        let d = (a.refinement.stats.chi2 - b.refinement.stats.chi2).abs() / b.refinement.stats.chi2;
        assert!(d < 2e-3, "{}: {} against {}", a.label, a.refinement.stats.chi2, b.refinement.stats.chi2);
    }
}

#[test]
fn series_outputs() {
    let (c, points) = chain(Start::Previous);
    // The CSV: one row per point; each value and su as the refinement gives them.
    let csv = drs::sequential::series_csv(&points);
    let rows: Vec<&str> = csv.lines().collect();
    assert_eq!(rows.len(), 17);
    let header: Vec<&str> = rows[0].split(',').collect();
    assert_eq!(&header[..9], &["point", "label", "parameter", "cycles", "stop", "chi2", "rp", "rwp", "rexp"]);
    let col = |name: &str| header.iter().position(|h| *h == name).unwrap_or_else(|| panic!("{name}"));
    for (k, pt) in points.iter().enumerate() {
        let f: Vec<&str> = rows[k + 1].split(',').collect();
        assert_eq!(f[1], pt.label);
        assert_eq!(f[col("chi2")].parse::<f64>().unwrap(), pt.refinement.stats.chi2);
        let j = pt.set.free.iter().position(|p| p.name == "phase1.cell_a").unwrap();
        assert_eq!(f[col("phase1.cell_a")].parse::<f64>().unwrap(), pt.refinement.params[j].value);
        assert_eq!(f[col("phase1.cell_a su")].parse::<f64>().unwrap(), pt.refinement.params[j].esd);
        assert_eq!(f[col("phase1.cell_volume")].parse::<f64>().unwrap(), drs::sequential::cell_volume(pt, 0).0);
    }
    // The parametric fit reads it back.
    let v = drs::parametric::column(&csv, "phase1.cell_volume", "parameter").unwrap();
    assert_eq!(v.len(), 16);
    // The CIF: an index block, then each point's results file, read back one by one.
    let cif = drs::sequential::series_cif(&c, &points, "la0").unwrap();
    // P43: each point's `_cell_volume` carries the table's su, the covariance's (to its printed digits).
    for (pt, text) in points.iter().zip(drs::sequential::split_series(&cif)) {
        let blocks = drs::cif::parse(&text);
        let v = blocks.iter().find_map(|b| b.get("_cell_volume")).unwrap();
        let (val, su) = drs::cif::value_su(v).unwrap();
        let (tv, tsu) = drs::sequential::cell_volume(pt, 0);
        let digits = 10f64.powi(-(v.split('.').nth(1).map_or(0, |d| d.split('(').next().unwrap().len()) as i32));
        // The su is written rounded up (`cif_write::with_su`): within one unit of its last digit above the table's.
        assert!((val - tv).abs() <= 0.5 * digits && su >= tsu - 1e-12 && su - tsu < digits, "{}: _cell_volume {v} against the table's {tv} ({tsu})", pt.label);
    }
    let blocks: Vec<&str> = cif.lines().filter(|l| l.starts_with("data_")).collect();
    let mut unique = blocks.clone();
    unique.sort();
    unique.dedup();
    assert_eq!(unique.len(), blocks.len(), "block names unique");
    let parts = drs::sequential::split_series(&cif);
    assert_eq!(parts.len(), 16);
    for (pt, text) in points.iter().zip(&parts) {
        let r = drs::results::read(text).unwrap();
        assert!((r.chi2 - pt.refinement.stats.chi2).abs() < 1e-6 * pt.refinement.stats.chi2, "{}", pt.label);
        assert!(r.free.iter().any(|q| q.0 == "phase1.cell_a"));
    }
    println!("{} lines of CSV, {} blocks in the CIF", rows.len(), blocks.len());
}

#[test]
fn parametric_fits() {
    // FullProf's values of the 71 patterns (the last block of each name), as the reference script reads them.
    let seq = drs::seq::read(&latin1(&dir().join("la0.seq"))).unwrap();
    let mut last: std::collections::BTreeMap<String, &drs::seq::SeqBlock> = Default::default();
    for b in &seq.blocks {
        last.insert(b.name.clone(), b);
    }
    let volume: Vec<(f64, f64, f64)> = last.values().map(|b| { let (v, s) = b.phases[0].volume.unwrap(); (b.parameter.unwrap(), v, s) }).collect();
    let c: Vec<(f64, f64, f64)> = last.values().map(|b| { let q = b.parameters.iter().find(|q| q.name == "Cell_C_ph1_pat1").unwrap(); (b.parameter.unwrap(), q.value, q.sigma) }).collect();
    let text = std::fs::read_to_string(Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/parametric_reference.txt")).unwrap();
    // Per fit: its name, points, χ² and each coefficient's value and su.
    type Reference = (String, usize, f64, Vec<(f64, f64)>);
    let mut refs: Vec<Reference> = Vec::new();
    for l in text.lines().filter(|l| !l.starts_with('#')) {
        let t: Vec<&str> = l.split_whitespace().collect();
        if t[0] == "fit" {
            refs.push((t[1].to_string(), t[3].parse().unwrap(), t[5].parse().unwrap(), Vec::new()));
        } else {
            refs.last_mut().unwrap().3.push((t[1].parse().unwrap(), t[2].parse().unwrap()));
        }
    }
    let s = |pairs: &[(&str, f64)]| pairs.iter().map(|(n, v)| (n.to_string(), *v)).collect::<Vec<_>>();
    let cases = [
        ("volume", "V0*(1 + a*(T - 150) + b*(T - 150)^2)", s(&[("V0", 245.0), ("a", 1e-5), ("b", 1e-8)]), &volume),
        ("cell_c", "c0 + c1*T", s(&[("c0", 7.6), ("c1", 1e-4)]), &c),
        ("einstein", "V0 + A/(exp(th/T) - 1)", s(&[("V0", 245.0), ("A", 1.0), ("th", 300.0)]), &volume),
    ];
    for (name, eq, start, points) in cases {
        let r = refs.iter().find(|x| x.0 == name).unwrap();
        let e = drs::parametric::Equation::parse(eq, "T").unwrap();
        let f = drs::parametric::fit(&e, points, &start).unwrap();
        assert_eq!(f.points, r.1);
        println!("{name}: {} cycles, χ² {:.6} (SciPy {:.6})", f.cycles, f.chi2, r.2);
        assert!((f.chi2 - r.2).abs() < 1e-9 * r.2, "{name}: χ² {} against SciPy's {}", f.chi2, r.2);
        for ((n, v, su), (rv, rsu)) in f.coefficients.iter().zip(&r.3) {
            println!("  {n:3} {v:.10e} ± {su:.6e}   SciPy {rv:.10e} ± {rsu:.6e}");
            assert!((v - rv).abs() < 1e-4 * rsu, "{name}: {n} {v} against SciPy's {rv} ± {rsu}");
            assert!((su - rsu).abs() < 1e-4 * rsu, "{name}: {n}'s su {su} against SciPy's {rsu}");
        }
    }
    // The parser: precedence, unary minus, functions, errors.
    let e = drs::parametric::Equation::parse("-a*T^2 + 2^b - sqrt(c)/2", "T").unwrap();
    assert_eq!(e.coefficients, ["a", "b", "c"]);
    assert_eq!(e.value(3.0, &[1.0, 3.0, 16.0]), -9.0 + 8.0 - 2.0);
    assert!(drs::parametric::Equation::parse("a*(T", "T").is_err());
    assert!(drs::parametric::Equation::parse("T^2", "T").is_err(), "no coefficient");
    assert!(drs::parametric::Equation::parse("foo(T)*a", "T").is_err());
}

/// App request 42: a coefficient that fits to 0 (a volume that does not change with T) is fitted, with its su, where
/// central differences with a step scaled by the coefficient stopped the fit as singular; the derivatives are the
/// parsed equation's, exact (checked here against central differences of every function at a generic point).
#[test]
fn parametric_coefficient_at_zero() {
    let e = drs::parametric::Equation::parse("V0 + a*T", "T").unwrap();
    // Three equal volumes (app 0.43.0's self-check): a = 0, its su 0 (no scatter).
    let flat = [(100.0, 245.1, 0.01), (150.0, 245.1, 0.01), (200.0, 245.1, 0.01)];
    let f = drs::parametric::fit(&e, &flat, &[("V0".into(), 245.1), ("a".into(), 0.0)]).unwrap();
    assert!(f.coefficients[1].1.abs() < 1e-12 && f.coefficients[1].2.is_finite(), "{:?}", f.coefficients);
    // A slope near 0 with scatter: the straight line's closed form (weighted least squares, equal su's).
    let pts = [(100.0, 245.100, 0.01), (125.0, 245.112, 0.01), (150.0, 245.095, 0.01), (175.0, 245.104, 0.01), (200.0, 245.101, 0.01)];
    let f = drs::parametric::fit(&e, &pts, &[("V0".into(), 245.1), ("a".into(), 0.0)]).unwrap();
    let n = pts.len() as f64;
    let (sx, sy) = (pts.iter().map(|p| p.0).sum::<f64>(), pts.iter().map(|p| p.1).sum::<f64>());
    let sxx: f64 = pts.iter().map(|p| (p.0 - sx / n).powi(2)).sum();
    let slope = pts.iter().map(|p| (p.0 - sx / n) * (p.1 - sy / n)).sum::<f64>() / sxx;
    let resid: f64 = pts.iter().map(|p| (p.1 - (sy / n + slope * (p.0 - sx / n))).powi(2)).sum();
    let su = (resid / (n - 2.0) / sxx).sqrt();
    println!("a = {:.6e} ± {:.6e}, closed form {slope:.6e} ± {su:.6e}", f.coefficients[1].1, f.coefficients[1].2);
    assert!((f.coefficients[1].1 - slope).abs() < 1e-9 * su && (f.coefficients[1].2 - su).abs() < 1e-6 * su);
    // The exact derivatives against central differences.
    let g = drs::parametric::Equation::parse("a*exp(b*T) + ln(c*T) - log10(c) + sqrt(a*T) + sin(b) * cos(c) / tan(a) + abs(b - 1) + sinh(b) * cosh(a) - tanh(c) + c^b + T^a", "T").unwrap();
    let c = [0.7, 0.3, 1.9];
    for k in 0..3 {
        let h = 1e-6;
        let (mut up, mut dn) = (c, c);
        up[k] += h;
        dn[k] -= h;
        let fd = (g.value(2.5, &up) - g.value(2.5, &dn)) / (2.0 * h);
        let d = g.derivative(2.5, &c, k);
        assert!((d - fd).abs() < 1e-7 * d.abs().max(1.0), "coefficient {k}: {d} against {fd}");
    }
}
