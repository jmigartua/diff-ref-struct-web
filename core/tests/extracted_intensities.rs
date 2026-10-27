//! The file of extracted intensities (P35): FullProf's `CODFILn_cltr.int` (manual, output files: Jbt = 2, More = 1,
//! Jvi = 11), written by `lebail::int_file` with FullProf's rule for the clusters of overlapped reflections (manual,
//! note on Jvi = 11: Rmub, Rmuc).
//!
//! - **The rule, against FullProf's four files whose run's `.hkl` is beside them** (LaMnO3's Le Bail fit, PPH3 and
//!   citric acid twice): the core writes the file from FullProf's own `.hkl` (2θ, FWHM, Iobs and its σ) and the
//!   result is compared line by line with FullProf's `.int`. LaMnO3's 414 lines agree. In the others two kinds of
//!   line differ, counted in `RULE`:
//!   - FullProf writes a near-zero observation as G = 0.0001 with σ = 1.0 (G up to 0.0003 in PPH3; scxm1's file has
//!     0.00076 with its own σ); the rule is not in the manual and the core keeps the value and its σ;
//!   - a cluster boundary decided at the `.hkl`'s print (2θ and FWHM to 10⁻⁴°, G to 10⁻³): the decision flips within
//!     half a unit of the last printed digit.
//! - **The core's own extraction** of LBF_lmopr (`fullprof_le_bail`), at FullProf's final values and `.hkl`
//!   intensities: the clusters against FullProf's `.int`, each cluster's G within `CORE` of FullProf's.
//! - `drs results --exports` writes `Exports/<control> <phase> intensities.int` for a Le Bail phase.

use drs::lebail::{clusters, int_file, Extracted, CLUSTER_DEFAULTS};
use std::path::{Path, PathBuf};

mod common;
use common::latin1;

fn manuals(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof").join(rel)
}

/// FullProf's runs with a `_cltr.int` and the `.hkl` it was written from: (folder, .hkl, .int, .pcr).
const RUNS: [(&str, &str, &str, &str); 4] = [
    ("7-FullProf-Tutorial-SymmetryModes/LaMnO3", "LBF_lmopr1.hkl", "LBF_lmopr1_cltr.int", "LBF_lmopr.pcr"),
    ("9-SAnn-OrganicCompounds/PPH3/Prof_Match/P21c", "PPH3-sls_11.hkl", "PPH3-sls_11_cltr.int", "PPH3-sls_1.pcr"),
    ("9-SAnn-OrganicCompounds/CitricAcid_example/PM", "citric-acid-sls_11.hkl", "citric-acid-sls_11_cltr.int", "citric-acid-sls_1.pcr"),
    ("9-SAnn-OrganicCompounds/CitricAcid_example/PM/P21c", "citric-acid-sls_11.hkl", "citric-acid-sls_11_cltr.int", "citric-acid-sls_1.pcr"),
];

/// A FullProf `.hkl` of profile matching (h k l, multiplicity, Iobs, σ, 2θ, FWHM per line, after two header lines).
fn hkl(text: &str) -> Vec<Extracted> {
    text.lines()
        .skip(2)
        .filter_map(|l| {
            let v: Vec<f64> = l.split_whitespace().map_while(|s| s.parse().ok()).collect();
            (v.len() >= 8).then(|| Extracted { hkl: [v[0] as i32, v[1] as i32, v[2] as i32], g: v[4], sigma: v[5], two_theta: v[6], fwhm: v[7] })
        })
        .collect()
}

/// The lines of an `.int`: h k l, G, σ, 2θ.
fn int_lines(text: &str) -> Vec<([i32; 3], f64, f64, f64)> {
    text.lines()
        .skip(3)
        .filter_map(|l| {
            let v: Vec<f64> = l.split_whitespace().filter_map(|s| s.parse().ok()).collect();
            (v.len() >= 9).then(|| ([v[0] as i32, v[1] as i32, v[2] as i32], v[3], v[4], v[8]))
        })
        .collect()
}

/// Rmub and Rmuc of the `.pcr`'s first phase (zeros: FullProf's defaults).
fn rule(pcr: &Path) -> (f64, f64) {
    let p = drs::pcr::read(&latin1(pcr)).unwrap();
    let f = &p.phases[0].flags;
    let (b, c) = (f.get("RMub").copied().unwrap_or(0.0), f.get("RMuc").copied().unwrap_or(0.0));
    if b == 0.0 && c == 0.0 {
        CLUSTER_DEFAULTS
    } else {
        (b, c)
    }
}

/// Per run: (lines, the same, near-zero lines FullProf floors, cluster boundaries decided at the `.hkl`'s print), as
/// measured in L256.
const RULE: [(usize, usize, usize, usize); 4] = [(414, 414, 0, 0), (479, 463, 14, 1), (243, 239, 4, 0), (243, 238, 5, 0)];

/// FullProf's clusters, read from its `.int`: a cluster ends at a line whose G is not −1.
fn fullprof_clusters(lines: &[([i32; 3], f64, f64, f64)]) -> Vec<std::ops::Range<usize>> {
    let mut out = Vec::new();
    let mut start = 0;
    for (k, l) in lines.iter().enumerate() {
        if l.1 != -1.0 {
            out.push(start..k + 1);
            start = k + 1;
        }
    }
    out
}

#[test]
fn fullprof_cluster_rule() {
    for ((dir, h, i, p), counts) in RUNS.iter().zip(RULE) {
        let d = manuals(dir);
        let rows = hkl(&latin1(&d.join(h)));
        let (rmub, rmuc) = rule(&d.join(p));
        let ours = int_lines(&int_file("test", 1.0, &rows, (rmub, rmuc)));
        let theirs = int_lines(&latin1(&d.join(i)));
        assert_eq!(ours.len(), theirs.len(), "{dir}");
        for (k, (o, t)) in ours.iter().zip(&theirs).enumerate() {
            assert_eq!(o.0, t.0, "{dir}: line {k}");
        }
        // The partitions: where they differ, the decision for the pair at the boundary must flip within the .hkl's
        // print (2θ and FWHM ± 0.00005°, G ± 0.0005).
        let (oc, tc) = (clusters(&rows, rmub, rmuc), fullprof_clusters(&theirs));
        let ends = |c: &[std::ops::Range<usize>]| c.iter().map(|r| r.end).collect::<std::collections::BTreeSet<_>>();
        let (oe, te) = (ends(&oc), ends(&tc));
        let mut ties = 0;
        let mut apart = std::collections::BTreeSet::new();
        for &e in oe.symmetric_difference(&te) {
            let (a, b) = (&rows[e - 1], &rows[e]);
            let (dt, half) = (b.two_theta - a.two_theta, 0.5 * (a.fwhm + b.fwhm));
            let near = (dt - half * rmub).abs() <= 1e-4 * (1.0 + rmub) || ((dt - half).abs() <= 1e-4 && b.g < 1.0);
            assert!(near, "{dir}: the cluster boundary after {:?} differs and is not decided at the print: ΔT {dt}, ½(H + H) {half}", a.hkl);
            ties += 1;
            // The lines of the clusters on both sides are left out of the line comparison.
            for c in oc.iter().chain(&tc).filter(|c| c.contains(&(e - 1)) || c.contains(&e)) {
                apart.extend(c.clone());
            }
        }
        let (mut same, mut floor) = (0, 0);
        for c in &oc {
            // Members printed 0.001 in the .hkl may be anything from 0 to 0.0015 (FullProf floors its print).
            let tiny = rows[c.clone()].iter().filter(|r| r.g <= 0.001).count() as f64;
            let n = c.len() as f64;
            for k in c.clone().filter(|k| !apart.contains(k)) {
                let (o, t) = (&ours[k], &theirs[k]);
                // The .hkl's print (half a unit of 10⁻³ per member), its floor, and FullProf's single precision.
                let g_ok = (o.1 - t.1).abs() <= 0.0005 * n + 0.001 * tiny + 2.4e-7 * t.1.abs() * n;
                let s_ok = (o.2 - t.2).abs() <= 0.0005 * n.sqrt() + 1e-6 || t.2 < 0.002;
                if g_ok && s_ok {
                    same += 1;
                } else if t.2 == 1.0 && t.1 <= 0.0003 + 1e-9 {
                    floor += 1;
                } else {
                    panic!("{dir} {:?}: ours G {} σ {}, FullProf's {} {}", o.0, o.1, o.2, t.1, t.2);
                }
            }
        }
        println!("{dir}: {} lines, {same} the same, {floor} near-zero lines floored by FullProf (G 0.0001, σ 1.0), {ties} cluster boundaries decided at the print ({} lines)", ours.len(), apart.len());
        assert_eq!((ours.len(), same, floor, ties), counts, "{dir}");
    }
}

/// The core's extraction of LBF_lmopr at FullProf's final values (its `.pcr`) and intensities (its `.hkl`), with
/// FullProf's reflection limit (S(1/d)max of its `.out`), against FullProf's `.int`, as measured in L256: (FullProf's
/// clusters, the core's that are FullProf's, the largest |ΔG| of a cluster in FullProf's σ, the largest relative
/// difference of σ). G: every cluster within 0.26 of FullProf's σ, the largest at 57.53° (G 0.043, −2.9 %); the first
/// reflection, 0 1 1 at 15.71°, at the data's low edge and weaker than its σ (FullProf 1.18 ± 1.43, the core 0.84 ±
/// 0.96), is −28 % (0.23 σ) and has the largest σ difference, −33 %. The other σ's: −0.9 … +1.5 %, median +1.3 %, but
/// one at 125.15° (+9.8 %); FullProf's σ formula is not in its manual (as in `reflection_diagnostics`, L100).
const CORE: (usize, usize, f64, f64) = (122, 122, 0.2585, 0.327);

#[test]
fn core_extraction() {
    let d = manuals("7-FullProf-Tutorial-SymmetryModes/LaMnO3");
    let out = latin1(&d.join("LBF_lmopr.out"));
    let p = drs::pcr::read(&latin1(&d.join("LBF_lmopr.pcr"))).unwrap();
    let l = drs::powder::Model::from_pcr(&p, &d, "LBF_lmopr").unwrap();
    let data = l.data.unwrap();
    let mut m = l.model;
    m.phases[0].reflection_limit = out.lines().filter_map(|l| l.split_once("S(1/d)max:")).map(|(_, r)| r.split_whitespace().next().unwrap().parse().unwrap()).next_back();
    let rows = drs::lebail::extracted(&m, &data, 0).unwrap();
    let text = int_file("Phase No.:   1 LeBail fit of LaMnO3", m.lambda[0], &rows, rule(&d.join("LBF_lmopr.pcr")));
    // The header as FullProf writes it: the format, then λ to five decimals with Itypdata = 0 and Ipow = 2.
    let fp = latin1(&d.join("LBF_lmopr1_cltr.int"));
    for (o, t) in text.lines().zip(fp.lines()).skip(1).take(2) {
        assert_eq!(o.trim_end(), t.trim_end());
    }
    let (ours, theirs) = (int_lines(&text), int_lines(&fp));
    let (oc, tc) = (fullprof_clusters(&ours), fullprof_clusters(&theirs));
    let key = |lines: &[([i32; 3], f64, f64, f64)], c: &std::ops::Range<usize>| lines[c.clone()].iter().map(|l| l.0).collect::<Vec<_>>();
    let (mut matched, mut dg, mut ds) = (0, 0.0f64, 0.0f64);
    for c in &tc {
        if let Some(o) = oc.iter().find(|o| key(&ours, o) == key(&theirs, c)) {
            matched += 1;
            let (g, s) = (ours[o.end - 1], theirs[c.end - 1]);
            dg = dg.max((g.1 - s.1).abs() / s.2);
            ds = ds.max((g.2 / s.2 - 1.0).abs());
        }
    }
    println!("LBF_lmopr: {} reflections (FullProf {}), {} clusters (FullProf {}), {matched} of FullProf's found; G within {dg:.3} of FullProf's σ, σ within {ds:.4}", ours.len(), theirs.len(), oc.len(), tc.len());
    assert_eq!((ours.len(), tc.len(), oc.len(), matched), (theirs.len(), CORE.0, CORE.0, CORE.1));
    assert!(dg <= CORE.2 && ds <= CORE.3, "G {dg}, σ {ds}");
}

/// `drs results --exports` writes a Le Bail phase's file, in FullProf's format, named by the file dictionary.
#[test]
fn intensities_exported() {
    let src = manuals("7-FullProf-Tutorial-SymmetryModes/LaMnO3");
    let d = Path::new(env!("CARGO_MANIFEST_DIR")).join("target/test-intensities");
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    for f in ["LBF_lmopr.pcr", "LBF_lmopr1.hkl", "LaMnO3.dat"] {
        std::fs::copy(src.join(f), d.join(f)).unwrap();
    }
    let o = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(["results", d.join("LBF_lmopr.pcr").to_str().unwrap(), "--exports", "--time", "2026-10-15 12.00.00"]).output().unwrap();
    let out = String::from_utf8_lossy(&o.stdout);
    assert!(o.status.success(), "{out}{}", String::from_utf8_lossy(&o.stderr));
    let p = drs::pcr::read(&latin1(&d.join("LBF_lmopr.pcr"))).unwrap();
    let rel = drs::formats::files::path("export_intensities", &[("control", "LBF_lmopr"), ("phase", p.phases[0].name.split_whitespace().collect::<Vec<_>>().join(" ").as_str())]).unwrap();
    assert!(out.contains(&rel), "{rel} not written:\n{out}");
    let text = std::fs::read_to_string(d.join(&rel)).unwrap();
    assert_eq!(text.lines().nth(1), Some("(3i4,2f16.5,i4,3f14.4)"));
    let lines = int_lines(&text);
    assert!(lines.len() > 300 && lines.iter().any(|l| l.1 == -1.0), "{} lines", lines.len());
    assert_eq!(drs::formats::files::kind_of(&rel).map(|k| k.id), Some("export_intensities"));
}

/// The rule on its own: the two conditions, the default Rmub and Rmuc, and a σ of 0 written as FullProf's near-zero line.
#[test]
fn cluster_conditions() {
    let r = |t: f64, g: f64| Extracted { hkl: [0, 0, 0], two_theta: t, fwhm: 0.2, g, sigma: g.sqrt() };
    // ΔT < ½(H + H)·Rmub: one cluster; a step further, two.
    assert_eq!(clusters(&[r(10.0, 100.0), r(10.19, 100.0)], 1.0, 0.2), vec![0..2]);
    assert_eq!(clusters(&[r(10.0, 100.0), r(10.21, 100.0)], 1.0, 0.2), vec![0..1, 1..2]);
    // With Rmub = 0.5, ΔT = 0.15 is joined only when the next is weaker than Rmuc of the cluster so far.
    assert_eq!(clusters(&[r(10.0, 100.0), r(10.15, 10.0)], 0.5, 0.2), vec![0..2]);
    assert_eq!(clusters(&[r(10.0, 100.0), r(10.15, 30.0)], 0.5, 0.2), vec![0..1, 1..2]);
    assert_eq!(CLUSTER_DEFAULTS, (1.0, 0.2));
    let text = int_file("t", 1.5406, &[r(10.0, 4.0), r(10.1, 9.0), Extracted { hkl: [1, 1, 1], two_theta: 20.0, fwhm: 0.2, g: 0.0, sigma: 0.0 }], CLUSTER_DEFAULTS);
    let l = int_lines(&text);
    assert_eq!((l[0].1, l[0].2), (-1.0, 2.0));
    assert_eq!((l[1].1, l[1].2), (13.0, 3.60555));
    assert_eq!((l[2].1, l[2].2), (0.0001, 1.0));
    assert_eq!(text.lines().nth(2), Some("  1.54060   0   2    0.0000"));
}
