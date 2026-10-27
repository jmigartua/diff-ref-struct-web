//! The reader registry against FullProf and GSAS-II reading the same files (L062).
//!
//! - `fullprof_reads`: every data file of manuals/fullprof/ that a single-pattern `.prf` was drawn from
//!   (`tests/data/fullprof_reads.txt`, by `make_fullprof_reads.py`) is read by the registry's format for the `.pcr`'s
//!   `Ins`; its points must be the `.prf`'s (2θ to 10⁻⁴, Yobs to the `.prf`'s printing of FullProf's single-precision
//!   value), and, where the `.out` written with the `.prf` is there, Σ w·y² over the non-excluded points must be
//!   FullProf's SumwYobsSQ to its four printed digits: that is FullProf's σ. w is 1/σ², or, for a `.pcr` with
//!   Iwg = 1, FullProf's maximum-likelihood weight 1/(Ycal·σ²/y) with the `.prf`'s Ycal. Five files differ, frozen
//!   below: FullProf keeps or leaves out the point at the end of an excluded region by a rule the runs do not settle.
//! - `gsas2_reads`: every file of `tests/data/gsas2_reads.txt` (by `make_gsas2_reads.py`, GSAS-II's own importers
//!   run on the files of the GSAS-II tutorials and tests, FullProf's files GSAS-II reads, and synthetic files for the
//!   formats with no real example) gives GSAS-II's number of points, sums of x, y and w = 1/σ², first and last
//!   points, and wavelengths. GSAS-II's TOF x is the file's divided by 100 (the registry keeps µs), and GSAS-II
//!   converts .qye to 2θ and multiplies y by 100: the comparison undoes these. Known departures of GSAS-II's reading,
//!   counted and frozen: it keeps NCHAN − 1 points of a GSAS STD bank (`j < Nch`); it reads the blank fields of a short
//!   line as zero-weight points (its test for blanks does not see the end of the line); it reads an ESD bank's values
//!   past NCHAN; and it stops a free-format file at x ≥ Tf, which with the step's rounding is one point more than
//!   FullProf's int((Tf − Ti)/step + 1.5). The registry follows the BANK record's NCHAN and FullProf's count. A
//!   free-format file's σ is FullProf's rule, not GSAS-II's, and is not compared here.
//! - `powder_cif_against_fullprof`: FullProf's own powder CIFs (`_prof.cif`, Rpa = −1) hold, for each point not
//!   excluded, its number in the data file, 2θ corrected for the zero shift, Yobs and σ: read by `powder_cif`, they
//!   must be the data file's points as the registry reads them (y and σ to the CIF's printing, x less one constant).
//! - `recognition`: every example file is identified as its format, or is ambiguous with its format among the
//!   candidates; the ambiguous ones are counted and frozen.
//! - `formats_without_example`: the formats with no real example file, listed and frozen.

mod common;

use common::latin1;
use drs::formats::{self, Axis, Probe};
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

fn root(rel: &str) -> String {
    format!("{}/../{rel}", env!("CARGO_MANIFEST_DIR"))
}

fn index(name: &str) -> String {
    std::fs::read_to_string(format!("{}/tests/data/{name}", env!("CARGO_MANIFEST_DIR"))).unwrap()
}

fn read_as(path: &str, id: &str) -> formats::Dataset {
    let bytes = std::fs::read(root(path)).unwrap_or_else(|e| panic!("{path}: {e}"));
    let f = formats::by_id(id).unwrap_or_else(|| panic!("no format {id}"));
    (f.read)(&Probe::new(path, &bytes)).unwrap_or_else(|e| panic!("{path} as {id}: {e}"))
}

/// Is the file identified as `id`, or ambiguous with `id` among the leading candidates? (true, false): identified;
/// (false, true): ambiguous with it; (false, false): neither.
fn recognised(path: &str, id: &str) -> (bool, bool, String) {
    let bytes = std::fs::read(root(path)).unwrap();
    let ident = formats::identify(&Probe::new(path, &bytes));
    let found = match ident.format() {
        Some(c) => c.id.to_string(),
        None => format!("ambiguous: {}", ident.ambiguous().iter().map(|c| c.id).collect::<Vec<_>>().join(", ")),
    };
    (ident.format().is_some_and(|c| c.id == id), ident.ambiguous().iter().any(|c| c.id == id), found)
}

/// A .prf point: 2θ, Yobs, Ycal, and half a unit of the last digit printed for Yobs.
struct PrfPoint {
    x: f64,
    yobs: f64,
    ycal: f64,
    half_unit: f64,
}

/// The .prf's points and excluded regions (their number is the last item of the header's third line).
fn prf(path: &str, start: usize, n: usize) -> (Vec<PrfPoint>, Vec<(f64, f64)>) {
    let text = latin1(Path::new(&root(path)));
    let lines: Vec<&str> = text.lines().collect();
    let counts: Vec<usize> = lines[2].split_whitespace().filter_map(|s| s.parse().ok()).collect();
    let nexcl = counts.last().copied().unwrap_or(0);
    let excl = (0..nexcl)
        .map(|k| {
            let v: Vec<f64> = lines[3 + k].split_whitespace().map(|s| s.parse().unwrap()).collect();
            (v[0], v[1])
        })
        .collect();
    let pts = lines[start - 1..start - 1 + n]
        .iter()
        .map(|l| {
            let s: Vec<&str> = l.split('\t').take(3).map(str::trim).collect();
            let decimals = s[1].split_once('.').map_or(0, |(_, d)| d.len());
            PrfPoint { x: s[0].parse().unwrap(), yobs: s[1].parse().unwrap(), ycal: s[2].parse().unwrap(), half_unit: 0.5 * 10f64.powi(-(decimals as i32)) }
        })
        .collect();
    (pts, excl)
}

#[test]
fn fullprof_reads() {
    let mut fails = Vec::new();
    let (mut files, mut points, mut sigma_checked) = (0, 0usize, 0);
    let mut frozen: Vec<(String, f64)> = Vec::new();
    let mut wavelengths = 0;
    for line in index("fullprof_reads.txt").lines().filter(|l| l.starts_with("READ")) {
        let f: Vec<&str> = line.split_whitespace().collect();
        let (id, data, prf_path, start, n) = (f[1], f[3], f[4], f[5].parse().unwrap(), f[6].parse().unwrap());
        let d = read_as(data, id);
        let p = &d.scans[0].pattern;
        let (pts, excl) = prf(prf_path, start, n);
        let Some(i0) = p.x.iter().position(|&x| (x - pts[0].x).abs() < 2e-3) else {
            fails.push(format!("{prf_path}: its first 2θ {} is not a point of {data}", pts[0].x));
            continue;
        };
        // FullProf prints Yobs divided by a power of ten when the counts are large.
        let fac = 10f64.powi((p.y[i0] / pts[0].yobs).log10().round() as i32);
        let mut bad = 0;
        for (k, q) in pts.iter().enumerate() {
            let i = i0 + k;
            // FullProf holds the intensities in single precision: 32914.99875 is 32915.0 there, printed 3292. (×10).
            if i >= p.len() || (p.x[i] - q.x).abs() > 1e-4 + 1e-9 || (p.y[i] as f32 as f64 - q.yobs * fac).abs() > q.half_unit * fac * 1.0001 + 1e-6 {
                if bad == 0 {
                    fails.push(format!("{prf_path}: point {k} is ({}, {}) in the .prf, ({:?}, {:?}) read from {data}", q.x, q.yobs * fac, p.x.get(i), p.y.get(i)));
                }
                bad += 1;
            }
        }
        points += pts.len();
        files += 1;
        // A wavelength found in the file's header (PSI's "lambda=", "WL:") must be the .pcr's λ1, which the .prf's
        // second line prints.
        if let Some(w) = d.scans[0].wavelength {
            let l1: f64 = latin1(Path::new(&root(prf_path))).lines().nth(1).unwrap().split_whitespace().nth(2).unwrap().parse().unwrap();
            if (w.lambda1 - l1).abs() > 5e-4 {
                fails.push(format!("{data}: the header's wavelength {} is not the .pcr's {l1}", w.lambda1));
            }
            wavelengths += 1;
        }
        if f[7] != "-" && bad == 0 {
            let want: f64 = f[7].parse().unwrap();
            let likelihood = f[9] == "1";
            // w = 1/σ², or with Iwg = 1 FullProf's maximum likelihood w = 1/(Ycal·σ²/y) (the core's Weighting).
            let terms: Vec<f64> = pts
                .iter()
                .enumerate()
                // Excluded regions closed, a ≤ x ≤ b, as the core's (see `frozen` below for the runs where FullProf's
                // end points differ).
                .filter(|(_, q)| !excl.iter().any(|&(a, b)| q.x >= a && q.x <= b))
                .map(|(k, q)| {
                    let i = i0 + k;
                    let (y, s2) = (p.y[i], p.sigma[i] * p.sigma[i]);
                    let var = if likelihood { q.ycal * fac * if y > 0.0 && s2 > 0.0 { s2 / y } else { 1.0 } } else { s2 };
                    if var.is_finite() && var > 0.0 {
                        y * y / var
                    } else {
                        0.0
                    }
                })
                .collect();
            let got: f64 = terms.iter().sum();
            // Four significant digits: half a unit of the fourth.
            let unit = 10f64.powf(want.log10().floor() - 3.0);
            if (got - want).abs() > 0.5 * unit * 1.0001 {
                frozen.push((prf_path.to_string(), (got - want).abs() / want));
                println!("  {prf_path}: Σ w·y² = {got:.5e} read as {id}, FullProf's SumwYobsSQ {want:.4e}");
            }
            sigma_checked += 1;
        }
    }
    println!("{files} files, {points} points as FullProf read them; σ checked against FullProf's SumwYobsSQ on {sigma_checked}; {wavelengths} wavelengths read from headers");
    assert!(fails.is_empty(), "{} differences:\n{}", fails.len(), fails.join("\n"));
    assert!(files >= 50 && sigma_checked >= 40);
    // Frozen (L062): five .prf files whose Σ w·y² is not FullProf's to half a unit of its fourth digit.
    // - Egypt_qa_Rietveld (3 files, 1.8·10⁻⁴): FullProf's sum keeps the point at 6.70, the end of the excluded region
    //   5.0–6.7 (6.9361e6 with it, 6.9347e6 without), while in K2NaAgCN, NdNi and Pd3MnD it leaves out the points at
    //   the ends of their regions; the rule for an end point is not in the manual and these runs do not settle it.
    // - citric-acid-sls_1 of the profile-matching tutorial (2 files, 2.3·10⁻⁵): SumYobs agrees to its four digits,
    //   SumwYobsSQ is 6.43265e8 for FullProf's 6.432e8; not explained (the tutorial's own .out, not a run of ours).
    let worst = frozen.iter().map(|f| f.1).fold(0.0, f64::max);
    println!("{} files frozen, worst relative difference {worst:.2e}", frozen.len());
    assert!(frozen.len() <= 5 && worst <= 1.9e-4, "frozen σ differences grew: {frozen:?}");
}

/// What GSAS-II's reader returns: x (TOF / 100; Q files in 2θ), y, w.
fn as_gsas2(id: &str, s: &formats::Scan, file: &str) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let p = &s.pattern;
    let w: Vec<f64> = p.sigma.iter().map(|&e| if e.is_finite() && e > 0.0 { 1.0 / (e * e) } else { 0.0 }).collect();
    let time_map = latin1(Path::new(&root(file))).contains("TIME_MAP");
    match (id, s.axis) {
        ("gsas", Axis::TimeOfFlight) if !time_map => (p.x.iter().map(|x| x / 100.0).collect(), p.y.clone(), w),
        ("qye", Axis::Q) => {
            let l = s.wavelength.unwrap().lambda1;
            (p.x.iter().map(|q| 2.0 * (l * q / (4.0 * std::f64::consts::PI)).asin().to_degrees()).collect(), p.y.iter().map(|y| 100.0 * y).collect(), w.iter().map(|w| w / 1e4).collect())
        }
        _ => (p.x.clone(), p.y.clone(), w),
    }
}

fn close(a: f64, b: f64, rel: f64) -> bool {
    (a - b).abs() <= rel * a.abs().max(b.abs()) + 1e-12
}

#[test]
fn gsas2_reads() {
    let mut fails = Vec::new();
    let mut departures: BTreeMap<&str, usize> = BTreeMap::new();
    let (mut banks, mut points, mut failed) = (0, 0, Vec::new());
    for line in index("gsas2_reads.txt").lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = line.split('\t').collect();
        let (file, id) = (f[0], f[2]);
        if f[3] == "FAILED" {
            failed.push(format!("{file}: {}", f[4]));
            continue;
        }
        let bank: usize = f[3].parse().unwrap();
        let n: usize = f[4].parse().unwrap();
        let sums: Vec<f64> = f[5..8].iter().map(|s| s.parse().unwrap()).collect();
        let first: Vec<f64> = f[8].split(' ').map(|s| s.parse().unwrap()).collect();
        let last: Vec<f64> = f[9].split(' ').map(|s| s.parse().unwrap()).collect();
        let d = read_as(file, id);
        let Some(s) = d.scans.get(bank - 1) else {
            fails.push(format!("{file}: bank {bank} not read ({} scans)", d.scans.len()));
            continue;
        };
        let (x, y, w) = as_gsas2(id, s, file);
        let text = latin1(Path::new(&root(file)));
        let record = text.lines().filter(|l| l.starts_with("BANK")).nth(bank - 1).unwrap_or("");
        // The common points, and why the counts differ.
        let prefix: Option<Vec<f64>> = f.get(11).filter(|p| **p != "-").map(|p| p.split(' ').map(|v| v.parse().unwrap()).collect());
        let mut sums = sums;
        let m = if x.len() == n {
            n
        } else if id == "gsas" && x.len() == n + 1 && !record.contains("ESD") && !record.contains("FXY") && !record.contains("ALT") {
            *departures.entry("GSAS STD bank: GSAS-II keeps NCHAN − 1 points").or_default() += 1;
            n
        } else if let Some(p) = prefix.as_ref().filter(|p| p[0] as usize == x.len()) {
            // GSAS-II reads the values of a bank past its NCHAN; the registry stops at NCHAN.
            let key = if p[2] == sums[1] { "GSAS bank: GSAS-II adds zero-weight points past NCHAN (blank fields)" } else { "GSAS bank: GSAS-II reads values past NCHAN" };
            *departures.entry(key).or_default() += 1;
            sums = p[1..].to_vec();
            x.len()
        } else if id == "gsas" && x.len() < n {
            *departures.entry("GSAS bank: GSAS-II adds zero-weight points for the blank fields of a short line").or_default() += 1;
            x.len()
        } else if id == "free_steps" && x.len() + 1 == n {
            *departures.entry("free format: GSAS-II reads to x ≥ Tf, one point more").or_default() += 1;
            x.len()
        } else {
            fails.push(format!("{file} bank {bank}: {} points read, GSAS-II reads {n}", x.len()));
            continue;
        };
        let with_prefix = prefix.as_ref().is_some_and(|p| p[0] as usize == m && m < n);
        let sum = |v: &[f64]| v[..m].iter().sum::<f64>();
        let compare_w = id != "free_steps";
        let full = m == n || with_prefix;
        let mut checks = vec![("Σy", sum(&y), sums[1]), ("first x", x[0], first[0]), ("first y", y[0], first[1])];
        if full {
            checks.push(("Σx", sum(&x), sums[0]));
        }
        if m == n {
            checks.push(("last x", x[m - 1], last[0]));
            checks.push(("last y", y[m - 1], last[1]));
        }
        if compare_w {
            checks.push(("Σw", sum(&w), sums[2]));
            checks.push(("first w", w[0], first[2]));
            if m == n {
                checks.push(("last w", w[m - 1], last[2]));
            }
        }
        for (what, got, want) in checks {
            if !close(got, want, 1e-9) {
                fails.push(format!("{file} bank {bank} ({id}): {what} = {got}, GSAS-II {want}"));
            }
        }
        if f[10] != "-" {
            let want: Vec<f64> = f[10].split(' ').map(|s| s.parse().unwrap()).collect();
            let digits = f[10].split(' ').next().unwrap().split_once('.').map_or(0, |(_, d)| d.len());
            let tol = 0.5 * 10f64.powi(-(digits as i32)) + 1e-12;
            match s.wavelength {
                None => fails.push(format!("{file}: no wavelength read, GSAS-II reads {}", f[10])),
                Some(l) => {
                    let got = [Some(l.lambda1), l.lambda2, l.ratio];
                    for (k, &wv) in want.iter().enumerate() {
                        if got[k].is_none_or(|g| (g - wv).abs() > tol) {
                            fails.push(format!("{file}: wavelength item {} is {:?}, GSAS-II {wv}", k + 1, got[k]));
                        }
                    }
                }
            }
        }
        banks += 1;
        points += m;
    }
    // Every scan the registry reads is one GSAS-II reads: the same number of banks, ranges or scans per file.
    let reference = index("gsas2_reads.txt");
    let mut per_file: BTreeMap<&str, (usize, &str)> = BTreeMap::new();
    for line in reference.lines().filter(|l| !l.starts_with('#') && !l.contains("\tFAILED\t")) {
        let f: Vec<&str> = line.split('\t').collect();
        let e = per_file.entry(f[0]).or_insert((0, f[2]));
        e.0 = e.0.max(f[3].parse().unwrap());
    }
    for (file, (n, id)) in &per_file {
        let got = read_as(file, id).scans.len();
        if got != *n {
            fails.push(format!("{file}: {got} scans read, GSAS-II reads {n}"));
        }
    }
    println!("{banks} banks, {points} points as GSAS-II read them");
    for (k, v) in &departures {
        println!("  departure: {k}: {v} banks");
    }
    for f in &failed {
        println!("  GSAS-II did not read {f}");
    }
    assert!(fails.is_empty(), "{} differences:\n{}", fails.len(), fails.join("\n"));
    // Frozen: the departures, and GSAS-II's failures: FullProf's profile CIFs are not valid CIF for PyCifRW
    // ("_diffrn_measurement_device_type?", a tag and its value run together), so GSAS-II reads none of the three.
    let want: BTreeMap<&str, usize> = [
        ("GSAS STD bank: GSAS-II keeps NCHAN − 1 points", 5),
        ("GSAS bank: GSAS-II adds zero-weight points for the blank fields of a short line", 1),
        ("GSAS bank: GSAS-II adds zero-weight points past NCHAN (blank fields)", 2),
        ("GSAS bank: GSAS-II reads values past NCHAN", 1),
        ("free format: GSAS-II reads to x ≥ Tf, one point more", 1),
    ]
    .into_iter()
    .collect();
    assert_eq!(departures, want);
    assert_eq!(failed.len(), 3);
}

#[test]
fn powder_cif_against_fullprof() {
    let cases = [
        ("manuals/fullprof/9-SAnn-OrganicCompounds/PPH3/Rietveld/PPH3-sls_prof.cif", "manuals/fullprof/9-SAnn-OrganicCompounds/PPH3/Rietveld/PPH3-sls.dat", "xys", 1.001),
        ("manuals/fullprof/9-SAnn-OrganicCompounds/CitricAcid_example/Rietveld/citric-acid-sls_1_prof.cif", "manuals/fullprof/9-SAnn-OrganicCompounds/CitricAcid_example/Rietveld/citric-acid-sls_1.dat", "xys", 1.001),
        ("manuals/fullprof/7-FullProf-Tutorial-SymmetryModes/LaMnO3/LaMnO3_prof.cif", "manuals/fullprof/7-FullProf-Tutorial-SymmetryModes/LaMnO3/LaMnO3.dat", "ill_d2b", 1.2272),
    ];
    let mut fails = Vec::new();
    let mut points = 0;
    for (cif, data, id, lambda) in cases {
        let c = read_as(cif, "powder_cif");
        let s = &c.scans[0];
        let d = read_as(data, id);
        let p = &d.scans[0].pattern;
        let blocks = drs::cif::parse(&latin1(Path::new(&root(cif))));
        let ids: Vec<usize> = blocks.iter().find_map(|b| b.column("_pd_proc_point_id")).unwrap().iter().map(|v| v.parse().unwrap()).collect();
        assert_eq!(ids.len(), s.pattern.len(), "{cif}");
        let shift = s.pattern.x[0] - p.x[ids[0] - 1];
        for (k, &i) in ids.iter().enumerate() {
            // y in single precision, as FullProf holds it (8955.49972 is 8955.5, printed 8956).
            let (x, y, e) = (p.x[i - 1], p.y[i - 1] as f32 as f64, p.sigma[i - 1]);
            // The CIF prints y and σ as "58(4)": σ rounded to one or two digits, y to the same place.
            let place = 10f64.powi((s.pattern.sigma[k].log10().floor() as i32 - 1).max(0));
            if (s.pattern.x[k] - x - shift).abs() > 1e-4 || (s.pattern.y[k] - y).abs() > 0.5 * place + 1e-6 || (s.pattern.sigma[k] - e.max(1.0)).abs() > 0.5 * place * 1.0001 + 0.5 {
                fails.push(format!("{cif} point {i}: ({}, {}, {}) in the CIF, ({}, {y}, {e}) in {data}", s.pattern.x[k], s.pattern.y[k], s.pattern.sigma[k], x + shift));
                break;
            }
        }
        if s.wavelength.map(|w| w.lambda1) != Some(lambda) {
            fails.push(format!("{cif}: wavelength {:?}, the CIF says {lambda}", s.wavelength));
        }
        points += ids.len();
    }
    println!("{points} points of FullProf's powder CIFs are their data files' points");
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}

#[test]
fn recognition() {
    let mut expected: Vec<(String, String)> = Vec::new();
    for line in index("fullprof_reads.txt").lines().filter(|l| l.starts_with("READ") || l.starts_with("FILE")) {
        let f: Vec<&str> = line.split_whitespace().collect();
        expected.push((f[3].to_string(), f[1].to_string()));
    }
    for line in index("gsas2_reads.txt").lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = line.split('\t').collect();
        expected.push((f[0].to_string(), f[2].to_string()));
    }
    expected.sort();
    expected.dedup();
    let (mut identified, mut ambiguous, mut fails) = (0, BTreeMap::<String, usize>::new(), Vec::new());
    for (file, id) in &expected {
        let (is, among, found) = recognised(file, id);
        if is {
            identified += 1;
        } else if among {
            *ambiguous.entry(format!("{id} ({found})")).or_default() += 1;
        } else {
            fails.push(format!("{file}: expected {id}, found {found}"));
        }
    }
    println!("{} files: {identified} identified as their format; ambiguous with it:", expected.len());
    for (k, v) in &ambiguous {
        println!("  {v} × {k}");
    }
    assert!(fails.is_empty(), "{}", fails.join("\n"));
    let total: usize = ambiguous.values().sum();
    // Frozen: the files left ambiguous are columns of numbers, which FullProf's Ins = 10, GSAS-II's xye and the table
    // all read (the .pcr's Ins, the extension or the user decides).
    assert_eq!(total, 17, "ambiguous files changed: {ambiguous:?}");
    assert!(ambiguous.keys().all(|k| k.starts_with("xys") || k.starts_with("xye")), "{ambiguous:?}");
}

#[test]
fn formats_without_example() {
    let mut with: BTreeSet<String> = BTreeSet::new();
    for line in index("fullprof_reads.txt").lines().filter(|l| !l.starts_with('#')) {
        with.insert(line.split_whitespace().nth(1).unwrap().to_string());
    }
    for line in index("gsas2_reads.txt").lines().filter(|l| !l.starts_with('#') && !l.contains("/synthetic/")) {
        with.insert(line.split('\t').nth(2).unwrap().to_string());
    }
    // Example files without a reference reading: Ins = 14 (FP_tests, no FullProf output).
    with.insert("isis_tof_banks".into());
    let without: Vec<&str> = formats::FORMATS.iter().map(|f| f.id).filter(|id| !with.contains(*id)).collect();
    let synthetic: BTreeSet<String> = index("gsas2_reads.txt").lines().filter(|l| l.contains("/synthetic/")).map(|l| l.split('\t').nth(2).unwrap().to_string()).collect();
    println!("formats with no real example file ({}):", without.len());
    for id in &without {
        let f = formats::by_id(id).unwrap();
        let note = if synthetic.contains(*id) {
            "synthetic file checked against GSAS-II"
        } else if f.not_read.is_some() {
            "recognised, not read"
        } else {
            "unit tests only"
        };
        println!("  {id:<18} {}  [{note}]", f.name);
    }
    assert_eq!(
        without,
        [
            "ill_d1b_old", "brookhaven", "dbws", "two_axis", "ill_d4_d20l", "socabim_uxd", "variable_time", "table", "bruker_raw_v1", "bruker_raw_v2", "bruker_raw_v4", "bruker_brml",
            "panalytical_xrdml", "rigaku_ras", "rigaku_rasx", "rigaku_txt", "qye", "fit2d_chi", "nexus_nxazint", "midas_zarr"
        ]
    );
}
