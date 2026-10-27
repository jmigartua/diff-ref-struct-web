//! Background points (R26, APP_WORKFLOWS W4; L078).
//! - `against_references`: the natural cubic spline through Y2O3's 62 points against SciPy's `CubicSpline(bc_type =
//!   "natural")`, and the Chebyshev and polynomial fits against numpy's `chebfit` and `polyfit`
//!   (`tests/data/background_references.txt`, by `make_background_references.py`).
//! - `fullprof_bac_files`: FullProf's own `.bac` files of the runs, in its two formats (manual p. 8), read; the first
//!   format written at the data's points, as FullProf needs it (L081).
//! - `files_and_the_control_file`: the points in our file or a `.bac`, named by the control file, read back to the same
//!   model, pattern and text as when they are written in it; a spline kept through the control file and through the
//!   `.pcr` export (Nba < −5); a `.bac` written at the data's points; a `.pcr` with Nba = 1 and its `.bac` gives the
//!   background of the same points written in it.
//! - `fullprof_interpolation`: the core's background against FullProf's `.prf` on the twelve files of
//!   `runs/background_runs.py` FullProf ran (Nba = 17, 6, 62, 1 with a `.bac`; −6, −8, −17, −31, −62; points ending
//!   early, a point past the data, Thmax before the data's end), within the precision FullProf prints it to, by a
//!   bound derived in `offset_interval` (L081); the `.prf`'s offset a quarter of the least background, within the
//!   single precision FullProf computes it in (L087); χ² against FullProf's, each difference frozen (L087). Every
//!   variant is reported before the test fails.
//! - `fullprof_bac_pairs`: why FullProf's χ² is 19.9 with the pairs format (L079): it reads the numbers in order as
//!   the intensities of the first format (L081).

use drs::background::{self, Function, Interpolation};
use drs::control::Control;
use drs::powder::{BackgroundModel, Model};
use std::path::{Path, PathBuf};

fn root(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)
}

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap().iter().map(|&b| b as char).collect()
}

fn y2o3_points() -> Vec<(f64, f64)> {
    let c = Control::from_text(&std::fs::read_to_string(root("tests/data/broken/base.drc")).unwrap()).unwrap();
    match c.model.background {
        BackgroundModel::Points(p) => p,
        _ => panic!("Y2O3's background is points"),
    }
}

#[test]
fn against_references() {
    let pts = y2o3_points();
    assert_eq!(pts.len(), 62);
    let text = std::fs::read_to_string(root("tests/data/background_references.txt")).unwrap();
    let (mut n_spline, mut worst) = (0, 0f64);
    for l in text.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = l.split_whitespace().collect();
        let num = |i: usize| f[i].parse::<f64>().unwrap();
        match f[0] {
            "spline" => {
                let v = background::interpolate(&pts, Interpolation::Spline, &[num(1)])[0];
                let d = (v - num(2)).abs() / num(2).abs().max(1.0);
                worst = worst.max(d);
                assert!(d < 1e-12, "spline at {}: {v} vs {}", num(1), num(2));
                n_spline += 1;
            }
            "chebyshev" => {
                let terms = num(1) as usize;
                let (c, _) = background::fit(&pts, Function::Chebyshev { from: num(2), to: num(3), terms }).unwrap();
                for (k, v) in c.iter().enumerate() {
                    assert!((v - num(4 + k)).abs() < 1e-9 * num(4 + k).abs().max(1.0), "chebyshev {terms}: term {k} {v} vs {}", num(4 + k));
                }
            }
            "polynomial" => {
                let terms = num(1) as usize;
                let (c, _) = background::fit(&pts, Function::Polynomial { origin: num(2), terms }).unwrap();
                for (k, v) in c.iter().enumerate() {
                    assert!((v - num(3 + k)).abs() < 1e-8 * num(3 + k).abs().max(1.0), "polynomial: term {k} {v} vs {}", num(3 + k));
                }
            }
            other => panic!("{other}"),
        }
    }
    println!("spline at {n_spline} points within {worst:.1e} of SciPy; Chebyshev (6, 12 terms) and polynomial fits as numpy's");
    assert!(n_spline > 300);
}

#[test]
fn fullprof_bac_files() {
    // Format 1 (start step end, then intensities): NdNi_HTO.bac; format 2 (a comment, then pairs): LaB6rx.bac.
    let nd = background::read_path(&root("../manuals/fullprof/runs/NdNi-HTO/NdNi_HTO.bac")).unwrap();
    let head: Vec<f64> = latin1(&root("../manuals/fullprof/runs/NdNi-HTO/NdNi_HTO.bac")).lines().next().unwrap().split_whitespace().take(3).map(|v| v.parse().unwrap()).collect();
    assert_eq!(nd.points.len(), ((head[2] - head[0]) / head[1]).round() as usize + 1);
    assert_eq!(nd.points[0], (head[0], 164.0));
    let lab = background::read_path(&root("../manuals/fullprof/runs/LaB6/LaB6rx.bac")).unwrap();
    assert_eq!(lab.points.len(), 7245);
    assert_eq!(lab.points[0], (10.0383, 25.6106));
    // Written and read back, to the bit: format 2 as it is; format 1 at the data's points, where it holds the
    // interpolated background (FullProf reads only this one, with a value at every data point; L081).
    let pts = y2o3_points();
    assert_eq!(background::read(&background::write_bac_pairs(&pts, "Y2O3")).unwrap().points, pts);
    let x: Vec<f64> = (0..13091).map(|i| 15.0 + 0.01 * i as f64).collect();
    let at = background::read(&background::write_bac_at(&pts, Interpolation::Linear, &x, "Y2O3").unwrap()).unwrap().points;
    assert_eq!(at.len(), x.len());
    let want = background::interpolate(&pts, Interpolation::Linear, &x);
    assert!(at.iter().zip(&x).zip(&want).all(|((p, x), w)| (p.0 - x).abs() < 1e-9 && p.1 == *w));
    let ours = background::read(&background::write_ours(&pts, Interpolation::Spline)).unwrap();
    assert_eq!((ours.points, ours.interpolation), (pts, Some(Interpolation::Spline)));
}

#[test]
fn files_and_the_control_file() {
    let tmp = std::env::temp_dir().join(format!("drs_background_{}", std::process::id()));
    std::fs::create_dir_all(&tmp).unwrap();
    let run = root("../manuals/fullprof/runs/Y2O3");
    let p = drs::pcr::read(&latin1(&run.join("y2o3_structure.pcr"))).unwrap();
    let inline = Control::from_pcr(&p, &run, "y2o3_structure", &tmp).unwrap();
    let data = Model::from_pcr(&p, &run, "y2o3_structure").unwrap().data.unwrap();
    let y0 = inline.model.calculate(&data.x).ycalc;
    let close = |a: &[f64], b: &[f64]| a.len() == b.len() && a.iter().zip(b).all(|(u, v)| (u - v).abs() <= 1e-9 * v.abs().max(1.0));
    for file in ["y2o3.drb", "y2o3.bac"] {
        let mut c = inline.clone();
        c.background_file = Some(file.to_string());
        let bac = file.ends_with(".bac");
        if bac {
            // A .bac holds the background at every data point: without the data it is not written.
            assert!(c.background_file_text(None).unwrap().is_err());
            // Its points are the data's, not the 62: the refinement's terms on the 62 no longer apply.
            c.refinement.parameters.retain(|t| !t.iter().any(|x| matches!(x.0, drs::jacobian::Field::Background(_))));
        }
        let (name, text) = c.background_file_text(Some(&data.x)).unwrap().unwrap();
        std::fs::write(tmp.join(&name), text).unwrap();
        let drc = tmp.join("y2o3.drc");
        let t = c.to_text();
        assert!(t.contains(&format!("file = \"{file}\"")) && !t.contains("background_points = ["), "{file}: the points stay out of the control file");
        std::fs::write(&drc, &t).unwrap();
        let o = drs::control::open(&drc).unwrap();
        if bac {
            assert!(close(&o.control.model.calculate(&data.x).background, &inline.model.calculate(&data.x).background), "{file}: the background at the data's points");
            assert!(close(&o.control.model.calculate(&data.x).ycalc, &y0), "{file}: the pattern");
        } else {
            assert_eq!(o.control.model.background, inline.model.background, "{file}");
            assert_eq!(o.control.model.calculate(&data.x).ycalc, y0, "{file}: the pattern");
            assert_eq!(o.control.refinement.build(&o.control.model).unwrap(), inline.refinement.build(&inline.model).unwrap());
        }
        assert_eq!(o.control.to_text(), t, "{file}: the text written again");
        assert!(drs::check::check_path(&drc, &drs::check::Options::default()).is_empty(), "{file}: checks");
    }
    // A spline: through the control file, and through the .pcr export (Nba = −62) and back.
    let mut s = inline.clone();
    let BackgroundModel::Points(pts) = s.model.background.clone() else { panic!() };
    s.model.background = BackgroundModel::Spline(pts.clone());
    let back = Control::from_text(&s.to_text()).unwrap();
    assert_eq!(back.model.background, s.model.background);
    let e = drs::control::pcr_export::export(&s, "spl", Some(&data)).unwrap();
    let p2 = drs::pcr::read(&e.pcr).unwrap();
    assert_eq!(p2.flag("Nba"), -62);
    let m2 = Model::from_pcr_with_data(&p2, &tmp, "spl", Some(&run.join("y2o3.dat"))).unwrap().model;
    // FullProf adds a point of its own at the data's last point with the last value (L081): the export says so, and
    // the .pcr reads back with it; written again, that point is FullProf's and is dropped (Nba = −62 again).
    assert!(e.notes.iter().any(|n| n.contains("FullProf adds a point of its own")), "{:?}", e.notes);
    let mut with_end = pts.clone();
    with_end.push((*data.x.last().unwrap(), pts[pts.len() - 1].1));
    assert_eq!(m2.background, BackgroundModel::Spline(with_end));
    let mut s2 = s.clone();
    s2.model.background = m2.background.clone();
    let e2 = drs::control::pcr_export::export(&s2, "spl", Some(&data)).unwrap();
    assert_eq!(drs::pcr::read(&e2.pcr).unwrap().flag("Nba"), -62);
    assert!(!e2.notes.iter().any(|n| n.contains("FullProf adds a point of its own")), "{:?}", e2.notes);
    let ys = s.model.calculate(&data.x).ycalc;
    assert_ne!(ys, y0, "the spline is not the linear interpolation");
    // Nba = 1: the .bac of background_runs.py holds the same 62 points; its pattern is that of the points inline.
    let bdir = root("../manuals/fullprof/runs/background");
    let pb = drs::pcr::read(&latin1(&bdir.join("y2o3_bac62.pcr"))).unwrap();
    let lb = Model::from_pcr(&pb, &bdir, "y2o3_bac62").unwrap();
    assert!(close(&lb.model.calculate(&data.x).background, &inline.model.calculate(&data.x).background));
    let cb = Control::from_pcr(&pb, &bdir, "y2o3_bac62", &bdir).unwrap();
    assert_eq!(cb.background_file.as_deref(), Some("y2o3_bac62.bac"));
    let _ = std::fs::remove_dir_all(&tmp);
}

/// The rows of a `.prf`: (2θ, Backg, the Backg field as printed), the first n rows, n being the second number of the
/// second line (the reflection rows follow them).
fn prf_background(text: &str) -> Vec<(f64, f64, String)> {
    let n: usize = text.lines().nth(1).unwrap().split_whitespace().nth(1).unwrap().parse().unwrap();
    let start = text.lines().position(|l| l.trim_start().starts_with("2Theta")).unwrap() + 1;
    let rows: Vec<(f64, f64, String)> = text
        .lines()
        .skip(start)
        .take(n)
        .map(|l| {
            let f: Vec<&str> = l.split('\t').map(str::trim).collect();
            (f[0].parse().unwrap(), f[4].parse().unwrap(), f[4].to_string())
        })
        .collect();
    assert_eq!(rows.len(), n);
    rows
}

/// Whether one offset c puts the core's background within the reference's precision of FullProf's at every point.
///
/// FullProf's `.prf` prints Backg to a fixed number of decimals (0.1 count for Y2O3), less a constant (L079, L081).
/// The bound is derived from that precision, not chosen:
/// - FullProf's value B_i is printed rounded: the printed P_i is within h = ½·10^(−decimals) of B_i (0.05 here);
/// - FullProf computes in single precision: B_i and B_i − c carry each at most half an ulp, u_i = (|B_i| + |B_i −
///   c|)·2^(−24) (taken with the core's value for B_i; 10⁻⁶ here);
/// - so the core is right when there is one c with |ours_i − c − P_i| ≤ h + u_i at every i, that is when the interval
///   [max_i(d_i − h − u_i), min_i(d_i + h + u_i)] of c, d_i = ours_i − P_i, is not empty.
///
/// The offset need not be estimated: an estimate from rounded values (the median of d, as 0.19.0 did) is itself off
/// by up to h, which is the 0.0003 by which the linear background failed against a bound of h (L079). Asking only that
/// some c exists has no such error. The interval is returned, to be compared with FullProf's constant.
fn offset_interval(ours: &[f64], prf: &[(f64, f64, String)]) -> (f64, f64, f64) {
    let decimals = prf[0].2.split('.').nth(1).map_or(0, str::len) as i32;
    assert!(prf.iter().all(|r| r.2.split('.').nth(1).map_or(0, str::len) as i32 == decimals), "one precision in the column");
    let h = 0.5 * 10f64.powi(-decimals);
    let (mut lo, mut hi) = (f64::NEG_INFINITY, f64::INFINITY);
    let mut spread = (f64::INFINITY, f64::NEG_INFINITY);
    for (o, r) in ours.iter().zip(prf) {
        let d = o - r.1;
        let u = (o.abs() + r.1.abs()) * 2f64.powi(-24);
        lo = lo.max(d - h - u);
        hi = hi.min(d + h + u);
        spread = (spread.0.min(d), spread.1.max(d));
    }
    (lo, hi, spread.1 - spread.0)
}

/// The variants of `runs/background_runs.py` whose background the core predicts (L082, L087), all run by FullProf
/// (L085): asserted.
const PREDICTED: [&str; 12] = [
    "y2o3_lin16", "y2o3_spl16", "y2o3_spl6", "y2o3_lin6", "y2o3_spl8", "y2o3_spl31", "y2o3_spl62", "y2o3_spl6early", "y2o3_lin62", "y2o3_spl6beyond", "y2o3_spl6thmax", "y2o3_bacfirst",
];
/// Read by FullProf as something else (its Backg is the pairs read as intensities, L082): reported only; its χ² is
/// checked in `fullprof_bac_pairs`.
const ASKED: [&str; 1] = ["y2o3_bacpairs"];

/// |χ² − FullProf's| per variant, measured and frozen (L087). The background is FullProf's to the width of the offset
/// interval (10⁻⁵); the rest is the Bragg part's difference that `fullprof_profiles` measures on the same Y2O3 run
/// (Ycal rms 6.8·10⁻⁵ relative, the same in every variant within the print), which moves χ² by an amount that depends on
/// each variant's residuals. Each is also within that test's bound for χ², 5·10⁻⁴ + 1.5·10⁻³ χ².
const FROZEN_CHI2: [(&str, f64); 12] = [
    ("y2o3_lin16", 0.00006),
    ("y2o3_spl16", 0.00003),
    ("y2o3_spl6", 0.00105),
    ("y2o3_lin6", 0.00093),
    ("y2o3_spl8", 0.00056),
    ("y2o3_spl31", 0.00022),
    ("y2o3_spl62", 0.00018),
    ("y2o3_spl6early", 0.00022),
    ("y2o3_lin62", 0.00017),
    ("y2o3_spl6beyond", 0.00083),
    ("y2o3_spl6thmax", 0.00105),
    ("y2o3_bacfirst", 0.00017),
];

#[test]
fn fullprof_interpolation() {
    // FullProf 8.20's background (runs/background_runs.py, run on the Mac; L079, L081). The spline variants are read
    // with FullProf's own point at the data's last point (powder::Model::from_pcr, L081).
    let dir = root("../manuals/fullprof/runs/background");
    let mut waiting = Vec::new();
    let mut fails: Vec<String> = Vec::new();
    for v in PREDICTED.iter().chain(&ASKED) {
        let prf = dir.join(format!("{v}.prf"));
        if !prf.exists() {
            waiting.push(*v);
            continue;
        }
        let text = latin1(&dir.join(format!("{v}.pcr")));
        let p = drs::pcr::read(&text).unwrap();
        let l = Model::from_pcr(&p, &dir, v).unwrap();
        let data = l.data.as_ref().unwrap();
        let prf_text = latin1(&prf);
        let rows = prf_background(&prf_text);
        let x: Vec<f64> = rows.iter().map(|r| r.0).collect();
        let ours = l.model.calculate(&x).background;
        let (lo, hi, spread) = offset_interval(&ours, &rows);
        // FullProf's constant: a quarter of the background's least value (L082), within the interval. FullProf's least
        // value is its own single-precision background at that point (|B|·2⁻²⁴ from ours, as each point above), and
        // the quarter is rounded to single precision once more (c·2⁻²⁴): the quarter may lie that much outside (L087).
        let least = ours.iter().copied().fold(f64::INFINITY, f64::min);
        let quarter = least / 4.0;
        let delta = least.abs() * 2f64.powi(-24) / 4.0 + quarter.abs() * 2f64.powi(-24);
        // χ² as the .prf's first line gives it, and the core's with the .pcr's number of refined parameters.
        let refined: usize = text.lines().find(|l| l.contains("!Number of refined")).and_then(|l| l.split_whitespace().next()?.parse().ok()).unwrap_or(0);
        let chi2 = drs::lsq::statistics(&l.model, data, &l.model.calculate(&data.x).ycalc, refined).chi2;
        let fp_chi2 = prf_text.lines().next().unwrap().split("Chi2:").nth(1).and_then(|t| t.split_whitespace().next()).unwrap_or("?");
        println!("{v}: {} points; ours − FullProf's spans {spread:.5}; one offset fits in [{lo:.6}, {hi:.6}]; a quarter of the least background {quarter:.6}; χ² {chi2:.6} (FullProf {fp_chi2})", rows.len());
        if let (Some((_, f)), Ok(fp)) = (FROZEN_CHI2.iter().find(|f| f.0 == *v), fp_chi2.parse::<f64>()) {
            println!("    χ² − FullProf's {:+.6} (frozen at {f})", chi2 - fp);
            if (chi2 - fp).abs() > *f || *f > 5e-4 + 1.5e-3 * fp {
                fails.push(format!("{v}: χ² {chi2} against FullProf's {fp}: the frozen difference {f} grew"));
            }
        }
        if PREDICTED.contains(v) {
            if lo > hi {
                fails.push(format!("{v}: no offset puts the background within FullProf's printed precision (spread {spread:.4})"));
            }
            let outside = (lo - quarter).max(quarter - hi).max(0.0);
            println!("    the quarter outside the interval by {outside:.2e} (single precision allows {delta:.2e})");
            if outside > delta {
                fails.push(format!("{v}: the offset is not a quarter of the least background ({quarter:.7} outside [{lo:.7}, {hi:.7}] by {outside:.2e} > {delta:.2e})"));
            }
        }
    }
    assert!(fails.is_empty(), "{}", fails.join("\n"));
    if !waiting.is_empty() {
        println!("waiting for FullProf on the Mac (runs/background_runs.py run): {}", waiting.join(", "));
    }
}

#[test]
fn fullprof_bac_pairs() {
    // Josu's runs (L079): the 62 points' linear background at the 13091 data points, in a .bac with Nba = 1, gave
    // χ² 19.9 in the pairs format (format 2) and 1.80 in the first. FullProf 8.20 reads format 2 as format 1 (the
    // manual's footnote, p. 8: "This last format is not available at present"): the numbers x₁ y₁ x₂ y₂ … in order,
    // as the intensities at the successive data points. The core's χ² under that reading is FullProf's, within the
    // 0.05 to which 19.9 is given.
    let bdir = root("../manuals/fullprof/runs/background");
    let text = latin1(&bdir.join("y2o3_bac62.pcr"));
    let p = drs::pcr::read(&text).unwrap();
    let l = Model::from_pcr(&p, &bdir, "y2o3_bac62").unwrap();
    let data = l.data.unwrap();
    let refined: usize = text.lines().find(|l| l.contains("!Number of refined")).unwrap().split_whitespace().next().unwrap().parse().unwrap();
    let pts = y2o3_points();
    let b = background::interpolate(&pts, Interpolation::Linear, &data.x);
    let chi2 = |values: &[f64]| {
        let mut m = l.model.clone();
        m.background = BackgroundModel::Points(data.x.iter().zip(values).map(|(x, v)| (*x, *v)).collect());
        drs::lsq::statistics(&m, &data, &m.calculate(&data.x).ycalc, refined).chi2
    };
    let pairs: Vec<f64> = data.x.iter().zip(&b).flat_map(|(x, v)| [*x, *v]).collect();
    let (as_values, as_pairs) = (chi2(&b), chi2(&pairs[..b.len()]));
    println!("the background at every point: χ² {as_values:.3}; its pairs read as format 1: {as_pairs:.3} (FullProf 19.9)");
    assert!((as_pairs - 19.9).abs() <= 0.05, "{as_pairs}");
    // The first format's 1.80 is not explained: the core gives 1.19 with the same values (L081); the variants
    // bac62 and bacpairs of background_runs.py, written by the core, are for FullProf to run again.
    assert!((as_values - 1.19).abs() < 0.01, "{as_values}");
}
