//! The data operations of the Data module (R21; L066), against independent references and exact identities.
//!
//! - `against_references`: on four real files (the AutoBkg tutorial's synchrotron scan, PbSO4 X-ray and garnet neutron
//!   data of GSAS-II, Y2O3 of FullProf), read by the core's readers, compared with `tests/data/data_ops_references.txt`
//!   (by `make_data_ops_references.py`): SciPy's Savitzky–Golay filter and the σ from its coefficients, SciPy's
//!   `find_peaks` and `peak_widths`, pybaselines' arPLS as GSAS-II calls it, GSAS-II's peak candidates, and numpy's
//!   rebinning. Relative agreement 1e-9 (savgol, widths, rebinning), 1e-8 of the largest count (arPLS: banded
//!   solvers in a different order), indices exactly.
//! - `identities`: conversions (d·Q = 2π, 2θ → d → 2θ and 2θ → Q → 2θ, a change of wavelength and back), σ
//!   propagation (normalising, adding and subtracting a scan from itself, merging a scan with itself), counts
//!   conserved by rebinning, masking and deleting.
//! - `non_destructive`: the operations change neither the scan given nor the file.
//! - `operations_in_the_control_file`: every operation written as `[[pattern.operation]]` and read back the same.

use drs::control::syntax;
use drs::data::Pattern;
use drs::data_ops::{self, Normalization, Operation, RebinMode};
use drs::formats::{self, Axis, Scan, Wavelength};

fn root(rel: &str) -> String {
    format!("{}/../{rel}", env!("CARGO_MANIFEST_DIR"))
}

fn no_file(_: &str, _: Option<&str>) -> Result<Scan, String> {
    Err("no file".into())
}

fn close(a: f64, b: f64, rel: f64) -> bool {
    (a - b).abs() <= rel * a.abs().max(b.abs()).max(1e-300)
}

#[test]
fn against_references() {
    let text = std::fs::read_to_string(root("core/tests/data/data_ops_references.txt")).unwrap();
    let mut scan: Option<Scan> = None;
    let mut name = String::new();
    let mut smoothed: std::collections::BTreeMap<(usize, usize), Pattern> = Default::default();
    let (mut peaks, mut expected_peaks): (Vec<data_ops::FoundPeak>, usize) = (Vec::new(), 0);
    let mut base: Option<data_ops::AutoBackground> = None;
    let mut counts = std::collections::BTreeMap::<&str, usize>::new();
    let mut worst_base = 0f64;
    let mut seen_peaks = 0;
    let finish = |peaks: &Vec<data_ops::FoundPeak>, seen: usize, expected: usize, name: &str| {
        assert_eq!(seen, expected, "{name}: SciPy finds {expected} peaks");
        assert_eq!(peaks.len(), expected, "{name}: the core finds {} peaks, SciPy {expected}", peaks.len());
    };
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = line.split('\t').collect();
        let num = |i: usize| f[i].parse::<f64>().unwrap();
        let int = |i: usize| f[i].parse::<usize>().unwrap();
        match f[0] {
            "file" => {
                if scan.is_some() {
                    finish(&peaks, seen_peaks, expected_peaks, &name);
                }
                name = f[1].to_string();
                let mut s = formats::read_path(std::path::Path::new(&root(f[1])), Some(f[2])).unwrap().scans.remove(0);
                // GSAS-II keeps NCHAN − 1 points of a GSAS STD bank, the registry NCHAN (L062, accepted): the
                // references are on GSAS-II's points.
                if f[2] == "gsas" && s.pattern.len() == int(3) + 1 {
                    for v in [&mut s.pattern.x, &mut s.pattern.y, &mut s.pattern.sigma] {
                        v.pop();
                    }
                }
                // The same data as GSAS-II's reader: points and Σy.
                assert_eq!(s.pattern.len(), int(3), "{name}");
                assert!(close(s.pattern.y.iter().sum::<f64>(), num(4), 1e-12), "{name}: Σy");
                smoothed.clear();
                scan = Some(s);
                seen_peaks = 0;
            }
            "savgol" => {
                let s = scan.as_ref().unwrap();
                let (pts, order, i) = (int(1), int(2), int(3));
                let sm = smoothed.entry((pts, order)).or_insert_with(|| data_ops::smooth(&s.pattern, pts, order).unwrap());
                assert!(close(sm.y[i], num(4), 1e-9) || (sm.y[i] - num(4)).abs() < 1e-9, "{name}: savgol {pts},{order} at {i}: {} vs {}", sm.y[i], num(4));
                if f[5] != "-" {
                    assert!(close(sm.sigma[i], num(5), 1e-9), "{name}: savgol σ {pts},{order} at {i}: {} vs {}", sm.sigma[i], num(5));
                }
                *counts.entry("savgol").or_default() += 1;
            }
            "peaks" => {
                let s = scan.as_ref().unwrap();
                peaks = data_ops::peak_search(&s.pattern, num(1), 0.0);
                expected_peaks = int(2);
            }
            "peak" => {
                let p = &peaks.get(seen_peaks).unwrap_or_else(|| panic!("{name}: peak {} not found", f[1]));
                assert_eq!(p.index, int(1), "{name}: peak {seen_peaks}");
                assert!(close(p.prominence, num(2), 1e-12), "{name}: prominence at {}", p.index);
                assert!(close(p.fwhm, num(3), 1e-9), "{name}: width at {}: {} vs {}", p.index, p.fwhm, num(3));
                seen_peaks += 1;
                *counts.entry("peak").or_default() += 1;
            }
            "arpls" => {
                let s = scan.as_ref().unwrap();
                let b = data_ops::auto_background(&s.pattern, None, 10).unwrap();
                assert_eq!(b.log_lambda, num(1), "{name}: GSAS-II's default log λ");
                assert_eq!(10f64.powf(b.log_lambda).trunc(), num(2), "{name}: λ");
                assert_eq!(b.iterations, int(3), "{name}: iterations");
                base = Some(b);
            }
            "base" => {
                let s = scan.as_ref().unwrap();
                let b = base.as_ref().unwrap();
                let ymax = s.pattern.y.iter().copied().fold(0.0, f64::max);
                let d = (b.curve[int(1)] - num(2)).abs() / ymax;
                worst_base = worst_base.max(d);
                assert!(d < 1e-8, "{name}: arPLS at {}: {} vs {}", f[1], b.curve[int(1)], num(2));
                *counts.entry("base").or_default() += 1;
            }
            "gsas2_candidates" => {
                let s = scan.as_ref().unwrap();
                let got = data_ops::gsas2_peak_candidates(&s.pattern, None);
                let want: Vec<usize> = f[2].split(' ').filter(|v| !v.is_empty()).map(|v| v.parse().unwrap()).collect();
                assert_eq!(want.len(), int(1));
                assert_eq!(got, want, "{name}: GSAS-II's candidates");
                *counts.entry("gsas2_candidates").or_default() += want.len();
            }
            "rebin" => {
                let s = scan.as_ref().unwrap();
                let mut r = s.clone();
                Operation::RebinFactor { factor: int(1), mode: RebinMode::Sum }.apply(&mut r, &no_file).unwrap();
                assert!(close(r.pattern.y[int(2)], num(3), 1e-12), "{name}: rebin {} at {}", f[1], f[2]);
                *counts.entry("rebin").or_default() += 1;
            }
            other => panic!("unknown line {other}"),
        }
    }
    finish(&peaks, seen_peaks, expected_peaks, &name);
    println!("checked: {counts:?}; arPLS worst difference {worst_base:.2e} of the largest count");
    assert!(counts.len() == 5 && counts.values().all(|&n| n > 0));
}

fn synthetic() -> Scan {
    // 2θ from 10 to 100 in steps of 0.25 (exact in binary), two peaks on a sloping background, σ = √y.
    let x: Vec<f64> = (0..=360).map(|i| 10.0 + 0.25 * i as f64).collect();
    let y: Vec<f64> = x.iter().map(|&t| 100.0 + 0.5 * t + 900.0 * (-((t - 35.0) / 0.4f64).powi(2)).exp() + 400.0 * (-((t - 61.0) / 0.5f64).powi(2)).exp()).collect();
    let sigma = y.iter().map(|v| v.sqrt()).collect();
    let mut s = Scan::new("synthetic", Axis::TwoTheta, Pattern { x, y, sigma });
    s.wavelength = Some(Wavelength::single(1.5406));
    s.info.monitor = Some(25000.0);
    s.info.counting_time = Some(4.0);
    s
}

fn apply(s: &Scan, ops: &[Operation]) -> Scan {
    data_ops::apply_all(s, ops, &no_file).unwrap()
}

#[test]
fn identities() {
    let s = synthetic();
    let n = s.pattern.len();
    let lambda = 1.5406;
    // d · Q = 2π, point by point (Q from d and d from 2θ).
    let d = apply(&s, &[Operation::Convert { to: Axis::D, wavelength: None }]);
    let q = apply(&s, &[Operation::Convert { to: Axis::Q, wavelength: None }]);
    assert_eq!(d.axis, Axis::D);
    assert_eq!((d.pattern.len(), q.pattern.len()), (n, n));
    for i in 0..n {
        // d decreases as 2θ increases: the points are sorted again, so d's point i is 2θ's point n − 1 − i.
        let dq = d.pattern.x[n - 1 - i] * q.pattern.x[i];
        assert!(close(dq, 2.0 * std::f64::consts::PI, 1e-13), "d·Q = {dq}");
        assert_eq!(d.pattern.y[n - 1 - i], s.pattern.y[i]);
        let th = (s.pattern.x[i] / 2.0).to_radians();
        assert!(close(d.pattern.x[n - 1 - i], lambda / (2.0 * th.sin()), 1e-14));
    }
    // There and back: 2θ → d → 2θ, 2θ → Q → 2θ, d → Q → d.
    for to in [Axis::D, Axis::Q] {
        let back = apply(&s, &[Operation::Convert { to, wavelength: None }, Operation::Convert { to: Axis::TwoTheta, wavelength: None }]);
        for i in 0..n {
            assert!(close(back.pattern.x[i], s.pattern.x[i], 1e-13), "{to:?}: {} vs {}", back.pattern.x[i], s.pattern.x[i]);
            assert_eq!((back.pattern.y[i], back.pattern.sigma[i]), (s.pattern.y[i], s.pattern.sigma[i]));
        }
    }
    let dqd = apply(&d, &[Operation::Convert { to: Axis::Q, wavelength: None }, Operation::Convert { to: Axis::D, wavelength: None }]);
    for i in 0..n {
        assert!(close(dqd.pattern.x[i], d.pattern.x[i], 1e-15));
    }
    // A change of wavelength keeps d: sin θ′/λ′ = sin θ/λ, and back to λ gives the scan again.
    let w = apply(&s, &[Operation::Wavelength { to: 0.7093, from: None }]);
    let wd = apply(&w, &[Operation::Convert { to: Axis::D, wavelength: None }]);
    for i in 0..n {
        assert!(close(wd.pattern.x[i], d.pattern.x[i], 1e-13));
    }
    let wb = apply(&w, &[Operation::Wavelength { to: lambda, from: None }]);
    for i in 0..n {
        assert!(close(wb.pattern.x[i], s.pattern.x[i], 1e-13));
    }
    // To a longer wavelength, the points past 2θ = 180° go.
    let long = apply(&s, &[Operation::Wavelength { to: 2.5, from: None }]);
    let kept = s.pattern.x.iter().filter(|&&t| 2.5 / lambda * (t / 2.0f64).to_radians().sin() <= 1.0).count();
    assert_eq!(long.pattern.len(), kept);
    assert!(kept < n);

    // Normalising: y and σ by the same factor.
    for (op, f) in [
        (Normalization::Monitor(100000.0), 4.0),
        (Normalization::Time(1.0), 0.25),
        (Normalization::Maximum(1.0), 1.0 / s.pattern.y.iter().copied().fold(0.0, f64::max)),
        (Normalization::Factor(3.0), 3.0),
    ] {
        let r = apply(&s, &[Operation::Normalize(op.clone())]);
        for i in 0..n {
            assert!(close(r.pattern.y[i], f * s.pattern.y[i], 1e-15) && close(r.pattern.sigma[i], f * s.pattern.sigma[i], 1e-15), "{op:?}");
        }
    }
    // Adding a scan to itself: 2y, σ√2; subtracting: 0, σ√2; merging with itself: y, σ/√2.
    let me = s.clone();
    let load = move |_: &str, _: Option<&str>| Ok(me.clone());
    let file = |op: &str| match op {
        "add" => Operation::Add { file: "x".into(), reader: None, scale: 1.0 },
        "subtract" => Operation::Subtract { file: "x".into(), reader: None, scale: 1.0 },
        _ => Operation::Merge { file: "x".into(), reader: None, scale: 1.0 },
    };
    let add = data_ops::apply_all(&s, &[file("add")], &load).unwrap();
    let sub = data_ops::apply_all(&s, &[file("subtract")], &load).unwrap();
    let mer = data_ops::apply_all(&s, &[file("merge")], &load).unwrap();
    assert_eq!((add.pattern.len(), sub.pattern.len(), mer.pattern.len()), (n, n, n));
    for i in 0..n {
        let (y, e) = (s.pattern.y[i], s.pattern.sigma[i]);
        assert_eq!(add.pattern.y[i], 2.0 * y);
        assert_eq!(sub.pattern.y[i], 0.0);
        assert!(close(add.pattern.sigma[i], e * 2f64.sqrt(), 1e-15) && close(sub.pattern.sigma[i], e * 2f64.sqrt(), 1e-15));
        assert!(close(mer.pattern.y[i], y, 1e-15) && close(mer.pattern.sigma[i], e / 2f64.sqrt(), 1e-15));
    }
    // Adding a scan on a grid shifted by half a step: linear interpolation, σ² = (1−t)²σ_a² + t²σ_b² with t = ½.
    let mut half = s.clone();
    for v in half.pattern.x.iter_mut() {
        *v += 0.125;
    }
    let load_half = move |_: &str, _: Option<&str>| Ok(half.clone());
    let ah = data_ops::apply_all(&s, &[file("add")], &load_half).unwrap();
    assert_eq!(ah.pattern.len(), n - 1, "the first point is outside the shifted scan");
    for k in 0..ah.pattern.len() {
        let i = k + 1;
        let (ya, yb) = (s.pattern.y[i - 1], s.pattern.y[i]);
        let (ea, eb) = (s.pattern.sigma[i - 1], s.pattern.sigma[i]);
        assert!(close(ah.pattern.y[k], s.pattern.y[i] + 0.5 * (ya + yb), 1e-14));
        assert!(close(ah.pattern.sigma[k], (s.pattern.sigma[i].powi(2) + 0.25 * (ea * ea + eb * eb)).sqrt(), 1e-14));
    }

    // Rebinning: counts conserved (sum), σ² summed; a new step of twice the old equals a factor of 2 here.
    let r3 = apply(&s, &[Operation::RebinFactor { factor: 3, mode: RebinMode::Sum }]);
    assert_eq!(r3.pattern.len(), n / 3);
    assert!(close(r3.pattern.y.iter().sum::<f64>(), s.pattern.y[..n / 3 * 3].iter().sum::<f64>(), 1e-15));
    let v3: f64 = s.pattern.sigma[..3].iter().map(|e| e * e).sum();
    assert!(close(r3.pattern.sigma[0], v3.sqrt(), 1e-15));
    let rs = apply(&s, &[Operation::RebinStep { step: 0.5, mode: RebinMode::Sum }]);
    // 361 points cover 90.25°: 180 bins of 0.5° and a quarter-degree left over.
    assert_eq!(rs.pattern.len(), 180);
    let covered: f64 = s.pattern.y[..n - 1].iter().sum();
    assert!(close(rs.pattern.y.iter().sum::<f64>(), covered, 1e-14), "{} vs {covered}", rs.pattern.y.iter().sum::<f64>());
    let rf = apply(&s, &[Operation::RebinFactor { factor: 2, mode: RebinMode::Sum }]);
    for j in 0..rs.pattern.len() {
        // The bins of 0.5° start half a step before the first point: each holds two whole points.
        assert!(close(rs.pattern.y[j], rf.pattern.y[j], 1e-14) && close(rs.pattern.sigma[j], rf.pattern.sigma[j], 1e-14), "bin {j}");
    }
    let ra = apply(&s, &[Operation::RebinFactor { factor: 2, mode: RebinMode::Average }]);
    for j in 0..ra.pattern.len() {
        assert!(close(ra.pattern.y[j], rf.pattern.y[j] / 2.0, 1e-15) && close(ra.pattern.x[j], s.pattern.x[2 * j] + 0.125, 1e-15));
    }

    // Masking: no point lost, σ = ∞ inside the range; deleting: the points go.
    let m = apply(&s, &[Operation::Mask { from: 20.0, to: 30.0 }]);
    let inside = s.pattern.x.iter().filter(|&&t| (20.0..=30.0).contains(&t)).count();
    assert_eq!(inside, 41);
    assert_eq!(m.pattern.len(), n);
    assert_eq!(m.pattern.sigma.iter().filter(|e| e.is_infinite()).count(), inside);
    assert_eq!(m.pattern.y, s.pattern.y);
    let del = apply(&s, &[Operation::Delete { from: 20.0, to: 30.0 }, Operation::DeletePoints { x: vec![10.0, 100.0] }]);
    assert_eq!(del.pattern.len(), n - inside - 2);
    let mp = apply(&s, &[Operation::MaskPoints { x: vec![35.0, 61.0] }]);
    assert_eq!(mp.pattern.sigma.iter().filter(|e| e.is_infinite()).count(), 2);

    // Smoothing a polynomial of the window's order leaves it unchanged, ends included.
    let mut cubic = s.clone();
    cubic.pattern.y = s.pattern.x.iter().map(|&t| 3.0 - 0.2 * t + 0.01 * t * t - 1e-4 * t * t * t).collect();
    let sm = apply(&cubic, &[Operation::Smooth { points: 9, order: 3 }]);
    for i in 0..n {
        assert!((sm.pattern.y[i] - cubic.pattern.y[i]).abs() < 1e-9 * cubic.pattern.y[i].abs().max(1.0), "{i}");
    }
    // Peaks of the synthetic scan: at 35° and 61°, with their widths (Gaussian FWHM = 2√ln2 · w).
    let pk = data_ops::peak_search(&s.pattern, 50.0, 3.0);
    assert_eq!(pk.iter().map(|p| p.x).collect::<Vec<_>>(), [35.0, 61.0]);
    // The background under the peaks: arPLS follows the slope away from them.
    let b = data_ops::auto_background(&s.pattern, None, 10).unwrap();
    let at = |t: f64| b.curve[s.pattern.x.iter().position(|&v| v == t).unwrap()];
    for t in [15.0, 50.0, 80.0] {
        assert!((at(t) - (100.0 + 0.5 * t)).abs() < 3.0, "background at {t}: {}", at(t));
    }
}

#[test]
fn non_destructive() {
    let path = root("manuals/fullprof/runs/Y2O3/y2o3.dat");
    let before = std::fs::read(&path).unwrap();
    let s = formats::read_path(std::path::Path::new(&path), Some("free_steps")).unwrap().scans.remove(0);
    let copy = s.clone();
    let ops = [
        Operation::Mask { from: 20.0, to: 21.0 },
        Operation::Smooth { points: 5, order: 2 },
        Operation::RebinFactor { factor: 2, mode: RebinMode::Sum },
        Operation::Normalize(Normalization::Maximum(1000.0)),
        Operation::Convert { to: Axis::Q, wavelength: Some(1.5406) },
    ];
    let out = data_ops::apply_all(&s, &ops, &no_file).unwrap();
    assert_eq!(s, copy);
    assert_ne!(out.pattern, s.pattern);
    assert_eq!(std::fs::read(&path).unwrap(), before);
}

#[test]
fn operations_in_the_control_file() {
    let ops = vec![
        Operation::Mask { from: 20.0, to: 21.5 },
        Operation::MaskPoints { x: vec![33.12, 40.0] },
        Operation::Delete { from: 150.0, to: 160.0 },
        Operation::DeletePoints { x: vec![10.02] },
        Operation::Smooth { points: 7, order: 2 },
        Operation::RebinFactor { factor: 3, mode: RebinMode::Average },
        Operation::RebinStep { step: 0.05, mode: RebinMode::Sum },
        Operation::Normalize(Normalization::Monitor(1e5)),
        Operation::Normalize(Normalization::Time(1.0)),
        Operation::Normalize(Normalization::Maximum(100.0)),
        Operation::Normalize(Normalization::Factor(0.5)),
        Operation::Add { file: "../empty can.dat".into(), reader: Some("xys".into()), scale: 1.0 },
        Operation::Subtract { file: "empty.xye".into(), reader: None, scale: 0.97 },
        Operation::Merge { file: "detector 2.xye".into(), reader: None, scale: 1.0 },
        Operation::Convert { to: Axis::D, wavelength: Some(1.5406) },
        Operation::Convert { to: Axis::Q, wavelength: None },
        Operation::Wavelength { to: 0.7093, from: None },
    ];
    let mut text = String::new();
    for op in &ops {
        text.push_str("[[pattern.operation]]\n");
        for (k, v) in op.to_pairs() {
            text.push_str(&format!("{k} = {v}\n"));
        }
    }
    let doc = syntax::parse(&text);
    assert!(!doc.has_errors(), "{:?}", doc.diagnostics);
    let back: Vec<Operation> = doc.tables.iter().filter(|t| t.array).map(|t| Operation::from_table(t).unwrap()).collect();
    assert_eq!(back, ops);

    // In a whole control file: written, read back, and applied when the file is opened.
    let dir = std::env::temp_dir().join(format!("drs-data-ops-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let runs = std::path::PathBuf::from(root("manuals/fullprof/runs/Y2O3"));
    let pcr_path = std::fs::read_dir(&runs).unwrap().map(|e| e.unwrap().path()).find(|p| p.extension().is_some_and(|e| e == "pcr") && !p.to_string_lossy().contains("_calc")).unwrap();
    let pcr = drs::pcr::read(&String::from_utf8_lossy(&std::fs::read(&pcr_path).unwrap())).unwrap();
    let stem = pcr_path.file_stem().unwrap().to_string_lossy().to_string();
    let mut c = drs::control::Control::from_pcr(&pcr, &runs, &stem, &dir).unwrap();
    let applied = vec![Operation::Mask { from: 20.0, to: 21.0 }, Operation::RebinFactor { factor: 2, mode: RebinMode::Sum }, Operation::Normalize(Normalization::Factor(0.5))];
    c.data.as_mut().unwrap().operations = applied.clone();
    let file = dir.join("y2o3.drc");
    std::fs::write(&file, c.to_text()).unwrap();
    let o = drs::control::open(&file).unwrap();
    assert_eq!(o.control.data.as_ref().unwrap().operations, applied);
    let raw = formats::read_path(o.data_path.as_ref().unwrap(), c.data.as_ref().unwrap().reader.as_deref()).unwrap().scans.remove(0);
    assert_eq!(o.data.unwrap(), data_ops::apply_all(&raw, &applied, &no_file).unwrap().pattern);
    std::fs::remove_dir_all(&dir).ok();
}

/// App request 34: the wavelength operation on a scan that states none (the Y₂O₃ tutorial's data, λ1 = 1.5406 Å in
/// its .pcr), from the λ the user gives; the same as the scan given that λ; a disagreeing source refused; in the
/// control file as `from`.
#[test]
fn app34_wavelength_from_a_given_source() {
    let s = formats::read_path(std::path::Path::new(&root("manuals/fullprof/runs/Y2O3/y2o3.dat")), None).unwrap().scans.remove(0);
    assert_eq!(s.wavelength, None, "the tutorial's data states no wavelength");
    let op = Operation::Wavelength { to: 0.7093, from: None };
    assert!(data_ops::apply_all(&s, std::slice::from_ref(&op), &no_file).unwrap_err().contains("states no wavelength"));
    let given = apply(&s, &[Operation::Wavelength { to: 0.7093, from: Some(1.5406) }]);
    let mut stated = s.clone();
    stated.wavelength = Some(Wavelength::single(1.5406));
    assert_eq!(given.pattern, apply(&stated, &[op]).pattern);
    assert_eq!(given.wavelength, Some(Wavelength::single(0.7093)));
    // A source equal to the scan's own is taken; a different one refused, not silently preferred.
    assert_eq!(apply(&stated, &[Operation::Wavelength { to: 0.7093, from: Some(1.5406) }]).pattern, given.pattern);
    let e = data_ops::apply_all(&stated, &[Operation::Wavelength { to: 0.7093, from: Some(1.54) }], &no_file).unwrap_err();
    assert!(e.contains("disagrees"), "{e}");
    assert!(data_ops::apply_all(&s, &[Operation::Wavelength { to: 0.7093, from: Some(0.0) }], &no_file).is_err());
    // The control file: written with `from`, read back the same; without it, as before.
    for op in [Operation::Wavelength { to: 0.7093, from: Some(1.5406) }, Operation::Wavelength { to: 0.7093, from: None }] {
        let text = format!("[[pattern.operation]]\n{}\n", op.to_pairs().iter().map(|(k, v)| format!("{k} = {v}")).collect::<Vec<_>>().join("\n"));
        let doc = syntax::parse(&text);
        let t = doc.tables.into_iter().find(|t| t.header.is_some()).unwrap();
        assert_eq!(Operation::from_table(&t).unwrap(), op, "{text}");
    }
}
