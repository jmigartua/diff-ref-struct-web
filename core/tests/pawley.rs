//! The Pawley fit (P36, `pawley`; Pawley 1981): a phase's reflection intensities refined by least squares with the
//! cell and profile. FullProf has no Pawley fit, so no FullProf run of one can be the reference; the references are:
//!
//! - **Synthetic patterns** (`recovered_from_its_own_pattern`): the spinel's pattern calculated from known
//!   intensities (the tied ones in proportion to their multiplicities), the Pawley fit started from Le Bail's
//!   partition of equal intensities, with the cell and W moved: every intensity, the cell and W come back to 10⁻⁶ of
//!   what made the pattern; and each intensity's derivative column is the pattern's change with it (the pattern is
//!   linear in them).
//! - **The coincident reflections** (`coincident_reflections`): the groups the averaged metric finds are those of
//!   GSAS-II's rule (`GetPawleyConstr`: equal h² + k² + l² in m-3 and m-3m; equal h² + k² and |l| in 4/m and 4/mmm;
//!   equal h² + hk + k² and |l| hexagonal and trigonal) on the cubic, tetragonal and hexagonal phases of FullProf's
//!   runs, and none in an orthorhombic one.
//! - **FullProf's Le Bail fits** (`against_fullprof_le_bail`), the same data and profile decomposed by Le Bail's
//!   partition: LaMnO₃ (LBF_lmopr, neutrons, Pnma, 414 reflections), the spinel (runs/MgAl2O4-LeBail, X-rays, Fd-3m,
//!   58 reflections, 14 groups tied) and Co-YBCO (runs/Co-YBCO-LeBail, neutrons, P4/mmm, 83), each from FullProf's
//!   last values and intensities, refining what FullProf varied (its 8 parameters for LaMnO₃, none for the others)
//!   and the intensities. Checked: the weighted sum of squares below Le Bail's (a least-squares minimum), and every
//!   refined parameter of LaMnO₃ FullProf's within the Pawley fit's esd. Frozen, with their evidence: the parameters
//!   beyond FullProf's own esd (Le Bail's covariance leaves the intensities out; the Pawley fit's esd's are 3.5 to 3.9
//!   times FullProf's for the cell); the isolated reflections beyond 2σ of Le Bail's intensity (the two estimators
//!   part where the profile misfits the peaks: Le Bail's partition is the background-subtracted area under the window,
//!   Pawley's the weighted least-squares value, lower in 47 of LaMnO₃'s 56); and ten Le Bail partitions at the Pawley
//!   solution against FullProf's intensities, which they give back to 0.04σ for LaMnO₃ and 0.17σ for the spinel (the
//!   median), but not for Co-YBCO, whose FullProf fit varied nothing and did not settle (L315, listed for Josu).
//! - **The results file** as Publish's P39 reads a Pawley fit, and **the control file**: each intensity located.

use drs::jacobian::Field;
use drs::lsq::{refine, Options};
use drs::powder::{Intensities, Model};
use drs::refinable::{ParameterSet, Refinement};
use std::path::Path;

mod common;
use common::latin1;

fn root(rel: &str) -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof").join(rel)
}

fn spinel() -> (drs::pcr::Pcr, Model, drs::data::Pattern) {
    let dir = root("4-Analysis-Spinel-MgAl2O4/pcr-files");
    let p = drs::pcr::read(&latin1(&dir.join("MgAl2O4s_LeBail.pcr"))).unwrap();
    let l = Model::from_pcr_with_data(&p, &dir, "MgAl2O4s_LeBail", Some(&root("4-Analysis-Spinel-MgAl2O4/Data-files/MgAl2O4s.dat"))).unwrap();
    (p, l.model, l.data.unwrap())
}

fn list(m: &Model, p: usize) -> &Vec<([i32; 3], usize, f64)> {
    let Intensities::List(l) = &m.phases[p].intensities else { panic!("phase {} has a structure", p + 1) };
    l
}

/// The weighted sum of squares of a model.
fn s(m: &Model, data: &drs::data::Pattern) -> f64 {
    let st = drs::lsq::statistics(m, data, &m.calculate(&data.x).ycalc, 0);
    st.chi2 * st.n as f64
}

#[test]
fn recovered_from_its_own_pattern() {
    let (_, m0, data) = spinel();
    let mut truth = m0.clone();
    let mut rf = Refinement::default();
    drs::pawley::setup(&mut truth, &data, &mut rf, 0, 0).unwrap();
    // Known intensities: 1000 · (1 + k mod 7) per multiplicity's share, the coincident ones with equal |F|².
    let groups = drs::pawley::coincident(&truth, 0);
    let Intensities::List(l) = &mut truth.phases[0].intensities else { panic!() };
    for (k, r) in l.iter_mut().enumerate() {
        r.2 = 1000.0 * (1 + k % 7) as f64;
    }
    for g in &groups {
        for &k in g {
            l[k].2 = 1000.0 * (1 + g[0] % 7) as f64 * l[k].1 as f64 / l[g[0]].1 as f64;
        }
    }
    let y = truth.calculate(&data.x).ycalc;
    let synthetic = drs::data::Pattern { y: y.clone(), sigma: y.iter().map(|v| v.max(1.0).sqrt()).collect(), ..data.clone() };

    // The derivative columns: the pattern is linear in each intensity.
    let fields: Vec<Field> = (0..5).map(|k| Field::Intensity(0, k)).collect();
    let (_, jac) = truth.derivatives(&data.x, &fields);
    for (j, &f) in fields.iter().enumerate() {
        let mut moved = truth.clone();
        moved.set(f, truth.get(f) + 1.0);
        let y1 = moved.calculate(&data.x).ycalc;
        let worst = (0..y.len()).map(|i| (y1[i] - y[i] - jac[j][i]).abs()).fold(0.0, f64::max);
        let top = jac[j].iter().fold(0.0f64, |a, b| a.max(b.abs()));
        assert!(top > 0.0 && worst < 1e-9 * top, "{f:?}: derivative off by {worst} of {top}");
    }

    // From equal intensities, the cell and W moved.
    let mut m = m0.clone();
    m.set_field(Field::Cell(0, 0), m0.phases[0].cell.a + 0.002).unwrap();
    m.phases[0].w *= 1.2;
    let mut rf = Refinement { parameters: vec![vec![(Field::Cell(0, 0), 1.0)], vec![(Field::W(0), 1.0)]], ..Default::default() };
    let setup = drs::pawley::setup(&mut m, &synthetic, &mut rf, 0, drs::pawley::START_ITERATIONS).unwrap();
    // A pattern without noise: χ² goes to zero and every esd with it, so the guard on shifts of many esd's is off.
    rf.steering.guards.jump = f64::INFINITY;
    assert_eq!((setup.reflections, setup.ties.len()), (61, 15));
    let set = rf.build(&m).unwrap();
    let r = refine(&mut m, &synthetic, &set, Options { max_cycles: 30, eps: 1e-4 }, |_| {});
    println!("{} parameters, {} cycles: χ² {:.3e}, {}", set.free.len(), r.cycles, r.stats.chi2, r.stop);
    let (a, w) = (m.phases[0].cell.a, m.phases[0].w);
    assert!((a - truth.phases[0].cell.a).abs() < 1e-6 * a && (w - truth.phases[0].w).abs() < 1e-6 * w, "a {a}, W {w}");
    let (got, want) = (list(&m, 0), list(&truth, 0));
    let seen: Vec<usize> = set.free.iter().flat_map(|p| p.fields.iter().filter_map(|f| if let Field::Intensity(_, k) = f.0 { Some(k) } else { None })).collect();
    for &k in &seen {
        assert!((got[k].2 - want[k].2).abs() < 1e-6 * want[k].2, "{:?}: {} against {}", got[k].0, got[k].2, want[k].2);
    }
    println!("{} intensities recovered, the cell and W", seen.len());
}

#[test]
fn coincident_reflections() {
    // GSAS-II's rule (GSASIIstrIO.GetPawleyConstr) by the Laue class's sums of indices; in the classes it has no sum
    // for (−1, 2/m, mmm) it ties reflections whose d differ by less than 10⁻⁴, which the core does not (they part as
    // the cell is refined): there no two reflections coincide for every cell.
    fn gsas2(laue: &str, a: [i32; 3], b: [i32; 3]) -> bool {
        let sq = |h: [i32; 3]| h[0] * h[0] + h[1] * h[1];
        let hex = |h: [i32; 3]| h[0] * h[0] + h[1] * h[1] + h[0] * h[1];
        let all = |h: [i32; 3]| h[0] * h[0] + h[1] * h[1] + h[2] * h[2];
        let cross = |h: [i32; 3]| h[0] * h[1] + h[0] * h[2] + h[1] * h[2];
        match laue {
            "m-3" | "m-3m" => all(a) == all(b),
            "4/m" | "4/mmm" => sq(a) == sq(b) && a[2].abs() == b[2].abs(),
            "-3R" | "-3mR" => all(a) == all(b) && cross(a) == cross(b),
            "-3" | "-3m" | "6/m" | "6/mmm" => hex(a) == hex(b) && a[2].abs() == b[2].abs(),
            _ => false,
        }
    }
    let base = spinel().1;
    let phases: [(&str, &str, [f64; 6]); 14] = [
        ("F d -3 m", "m-3m", [8.083, 8.083, 8.083, 90.0, 90.0, 90.0]),
        ("P a -3", "m-3", [5.6, 5.6, 5.6, 90.0, 90.0, 90.0]),
        ("P 4/m", "4/m", [5.1, 5.1, 7.3, 90.0, 90.0, 90.0]),
        ("I 41/a", "4/m", [5.2, 5.2, 11.4, 90.0, 90.0, 90.0]),
        ("P 4/m m m", "4/mmm", [3.87, 3.87, 11.66, 90.0, 90.0, 90.0]),
        ("P 6/m", "6/m", [6.2, 6.2, 4.1, 90.0, 90.0, 120.0]),
        ("P 63/m m c", "6/mmm", [9.0, 9.0, 5.2, 90.0, 90.0, 120.0]),
        ("P -3", "-3", [6.0, 6.0, 7.2, 90.0, 90.0, 120.0]),
        ("R -3 c", "-3m", [4.76, 4.76, 13.0, 90.0, 90.0, 120.0]),
        ("R -3:R", "-3R", [5.5, 5.5, 5.5, 55.3, 55.3, 55.3]),
        ("R -3 c:R", "-3mR", [5.13, 5.13, 5.13, 55.3, 55.3, 55.3]),
        ("P n m a", "mmm", [5.75, 7.69, 5.54, 90.0, 90.0, 90.0]),
        ("P 1 21/c 1", "2/m", [5.1, 6.3, 7.2, 90.0, 101.0, 90.0]),
        ("P -1", "-1", [5.1, 6.3, 7.2, 81.0, 95.0, 103.0]),
    ];
    let mut checked = 0;
    for (symbol, laue, c) in phases {
        let (group, _) = drs::spacegroup::SpaceGroup::from_symbol(symbol).unwrap();
        let cell = drs::cell::Cell { a: c[0], b: c[1], c: c[2], alpha: c[3], beta: c[4], gamma: c[5] };
        let mut m = base.clone();
        m.phases[0].intensities = Intensities::List(drs::start::le_bail_list(&group, &cell, 0.7093, 100.0));
        m.phases[0].group = group;
        m.phases[0].cell = cell;
        m.phases[0].extract = false;
        let l = list(&m, 0);
        let groups = drs::pawley::coincident(&m, 0);
        let group_of = |k: usize| groups.iter().position(|g| g.contains(&k));
        let mut accidental = 0;
        for i in 0..l.len() {
            for j in i + 1..l.len() {
                let ours = group_of(i).is_some() && group_of(i) == group_of(j);
                assert_eq!(ours, gsas2(laue, l[i].0, l[j].0), "{symbol}: {:?} and {:?}", l[i].0, l[j].0);
                let (di, dj) = (cell.d_spacing(l[i].0.map(|v| v as f64)).unwrap(), cell.d_spacing(l[j].0.map(|v| v as f64)).unwrap());
                if !ours && (di - dj).abs() / di < 1e-4 {
                    accidental += 1;
                }
                checked += 1;
            }
        }
        println!("{symbol:12} ({laue:5}): {:3} reflections, {:2} groups tied; {accidental} pairs within 10⁻⁴ in d at this cell, not tied", l.len(), groups.len());
        assert!(matches!(laue, "mmm" | "2/m" | "-1") || !groups.is_empty(), "{symbol}: nothing tied");
    }
    assert!(checked > 10000, "{checked}");
}

/// What a run is held to: the parameters beyond FullProf's esd (all within the Pawley fit's); the isolated reflections
/// beyond 2σ of Le Bail's intensity, at most; and Le Bail's partitions at the Pawley solution against FullProf's
/// intensities, median and largest in σ, at most.
struct Frozen {
    beyond: &'static [&'static str],
    over: usize,
    partition: (f64, f64),
}

#[test]
fn against_fullprof_le_bail() {
    let mut fails = Vec::new();
    let runs: [(&str, &str, Option<&str>, Frozen); 3] = [
        ("7-FullProf-Tutorial-SymmetryModes/LaMnO3", "LBF_lmopr", None, Frozen { beyond: &["phase1.cell_b", "phase1.cell_c"], over: 2, partition: (0.05, 0.74) }),
        ("runs/MgAl2O4-LeBail", "MgAl2O4s_LeBail_conv", Some("MgAl2O4s.dat"), Frozen { beyond: &[], over: 4, partition: (0.18, 1.50) }),
        // FullProf's Le Bail fit of Co-YBCO with the profile it starts from (nothing varied, χ² 25.7, its last cycle
        // swinging between 24.5 and 25.7): the core's partitions from FullProf's intensities do not stay at them
        // (a P35 matter, listed for Josu), so that run holds the Pawley fit to its sum of squares only.
        ("runs/Co-YBCO-LeBail", "Cocuc_LeBail_conv", Some("Cocuc.dat"), Frozen { beyond: &[], over: 19, partition: (0.97, 6.4) }),
    ];
    for (rel, code, data_file, frozen) in runs {
        let dir = root(rel);
        let out = latin1(&dir.join(format!("{code}.out")));
        let p = drs::pcr::read(&latin1(&dir.join(format!("{code}.pcr")))).unwrap();
        let loaded = Model::from_pcr_with_data(&p, &dir, code, data_file.map(|f| dir.join(f)).as_deref()).unwrap();
        let data = loaded.data.as_ref().unwrap();
        let mut m = loaded.model.clone();
        let lim: Vec<f64> = out.lines().filter_map(|l| l.split_once("S(1/d)max:")).map(|(_, r)| r.split_whitespace().next().unwrap().parse().unwrap()).collect();
        m.phases[0].reflection_limit = lim.last().copied();
        let le_bail = s(&m, data);
        // What FullProf varied ("Number of Least-Squares parameters varied"): none in the two runs of L055's files.
        let varied: usize = out.lines().find_map(|l| l.split_once("Number of Least-Squares parameters varied:")).and_then(|(_, r)| r.trim().parse().ok()).unwrap();
        let mut rf = if varied > 0 { Refinement::from_pcr(&p, &m).unwrap() } else { Refinement::default() };
        let setup = drs::pawley::setup(&mut m, data, &mut rf, 0, 0).unwrap();
        let start = m.clone();
        let set = rf.build(&m).unwrap();
        let r = refine(&mut m, data, &set, Options { max_cycles: 30, eps: 0.1 }, |_| {});
        let pawley = s(&m, data);
        println!("{code}: {} reflections, {} parameters ({varied} of FullProf's), {} groups tied, {} cycles ({}): S {pawley:.1} against Le Bail's {le_bail:.1}", setup.reflections, set.free.len(), setup.ties.len(), r.cycles, r.stop);
        if !r.converged || pawley >= le_bail {
            fails.push(format!("{code}: S {pawley:.1} against Le Bail's {le_bail:.1}, {}", r.stop));
        }
        let fpp = common::fullprof_params(&out);
        let mut beyond = Vec::new();
        for (fp, q) in set.free.iter().zip(&r.params).filter(|(f, _)| !matches!(f.fields[0].0, Field::Intensity(..))) {
            let f = fp.fields[0].0;
            let Some((n, &(v, e))) = fpp.iter().find(|(n, _)| common::fullprof_field(n, &[vec![]]) == Some(f)) else {
                fails.push(format!("{code}: {} has no FullProf value", fp.name));
                continue;
            };
            let (ours, theirs) = ((q.value - v).abs() / q.esd, (q.value - v).abs() / e);
            println!("  {:26} {:12.6} ({:.1e}) FullProf {n:18} {v:12.6} ({e:.1e}): {ours:.2} of ours, {theirs:.2} of FullProf's esd", fp.name, q.value, q.esd);
            if ours > 1.0 {
                fails.push(format!("{code}: {} {:.6} against FullProf's {v}: {ours:.2} of the Pawley esd", fp.name, q.value));
            }
            if theirs > 1.0 {
                beyond.push(fp.name.clone());
            }
        }
        if beyond != frozen.beyond {
            fails.push(format!("{code}: beyond FullProf's esd {beyond:?}, frozen {:?}", frozen.beyond));
        }
        let c = isolated(&start, &m, &set, &r, data, 0);
        println!("{code}: {} isolated reflections: |ΔI|/σ median {:.3}, beyond 2σ {:?}; Pawley below Le Bail in {}; Le Bail's partitions at the Pawley solution against FullProf's: median {:.3}σ, largest {:.3}σ", c.count, c.median, c.over, c.lower, c.partition.0, c.partition.1);
        if c.over.len() > frozen.over || c.partition.0 > frozen.partition.0 || c.partition.1 > frozen.partition.1 {
            fails.push(format!("{code}: beyond 2σ {:?} (frozen {}), partitions {:.3}σ, {:.3}σ (frozen {:?})", c.over, frozen.over, c.partition.0, c.partition.1, frozen.partition));
        }
    }
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}

/// The Pawley fit against Le Bail's on the reflections isolated by FullProf's rule (clusters of one, Rmub 1.0, Rmuc
/// 0.2), in the Pawley fit's σ.
struct Isolated {
    count: usize,
    median: f64,
    over: Vec<[i32; 3]>,
    /// How many the Pawley fit puts below Le Bail.
    lower: usize,
    /// Ten Le Bail partitions at the Pawley solution's cell and profile against FullProf's intensities: median and
    /// largest.
    partition: (f64, f64),
}

fn isolated(start: &Model, end: &Model, set: &ParameterSet, r: &drs::lsq::Refinement, data: &drs::data::Pattern, p: usize) -> Isolated {
    let rows = drs::lebail::extracted(start, data, p).unwrap();
    let (l0, l1) = (list(start, p), list(end, p));
    let mut lb = end.clone();
    lb.phases[p].extract = true;
    for _ in 0..10 {
        lb.le_bail(data);
    }
    let l2 = list(&lb, p);
    let mut devs: Vec<(f64, [i32; 3], bool, f64)> = Vec::new();
    for c in drs::lebail::clusters(&rows, 1.0, 0.2).into_iter().filter(|c| c.len() == 1) {
        let h = rows[c.start].hkl;
        let k = l0.iter().position(|e| e.0 == h).unwrap();
        let Some(j) = set.free.iter().position(|q| q.fields[0].0 == Field::Intensity(p, k)) else { continue };
        let e = r.params[j].esd;
        devs.push(((l1[k].2 - l0[k].2).abs() / e, h, l1[k].2 < l0[k].2, (l2[k].2 - l0[k].2).abs() / e));
    }
    devs.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut part: Vec<f64> = devs.iter().map(|d| d.3).collect();
    part.sort_by(f64::total_cmp);
    let n = devs.len();
    Isolated { count: n, median: devs[n / 2].0, over: devs.iter().filter(|d| d.0 > 2.0).map(|d| d.1).collect(), lower: devs.iter().filter(|d| d.2).count(), partition: (part[n / 2], part[n - 1]) }
}

/// The results file as Publish's P39 reads a Pawley fit (drs-publish 0.20.0's T10, T11, F9, D9): the intensities among
/// the refined and the free parameters, named `phase1.reflection_intensity.h_k_l`; in the reflection loop each
/// intensity measured is the one refined and its σ the su from the covariance (the tied ones in proportion, one
/// parameter moving them); the ties in `_drs_pawley_tie_*`, read back as the control file holds them.
#[test]
fn results_file_as_publish_reads_it() {
    let dir = root("runs/MgAl2O4-LeBail");
    let p = drs::pcr::read(&latin1(&dir.join("MgAl2O4s_LeBail_conv.pcr"))).unwrap();
    let data = Model::from_pcr_with_data(&p, &dir, "MgAl2O4s_LeBail_conv", Some(&dir.join("MgAl2O4s.dat"))).unwrap().data.unwrap();
    let mut c = drs::control::Control::from_pcr(&p, &dir, "MgAl2O4s_LeBail_conv", &dir).unwrap();
    c.refinement.parameters.push(vec![(Field::Cell(0, 0), 1.0)]);
    drs::pawley::setup(&mut c.model, &data, &mut c.refinement, 0, 0).unwrap();
    let set = c.refinement.build(&c.model).unwrap();
    let r = refine(&mut c.model, &data, &set, Options { max_cycles: 10, eps: 0.1 }, |_| {});
    let (text, structures) = drs::results::write(&c, &data, Some(&set), Some(&r), "spinel pawley").unwrap();
    assert!(structures.is_empty(), "a phase without a structure has no structure CIF");
    let back = drs::results::read(&text).unwrap();
    // 1. The names.
    let named = |v: &[(String, f64, f64)]| v.iter().filter(|q| q.0.starts_with("phase1.") && q.0.contains("intensit")).count();
    assert_eq!(named(&back.free), set.free.len() - 1, "every intensity parameter among the free ones");
    assert!(named(&back.refined) >= set.free.len() - 1, "and among the refined values");
    assert!(back.free.iter().any(|q| q.0 == "phase1.reflection_intensity.1_1_1"));
    // 2. Each reflection's intensity and σ.
    let l = list(&c.model, 0);
    let mut seen = 0;
    for row in back.reflections.iter().filter(|x| x.0 == 1) {
        let k = l.iter().position(|e| e.0 == row.2).unwrap();
        let f = Field::Intensity(0, k);
        let Some((j, q)) = set.free.iter().enumerate().find(|(_, q)| q.fields.iter().any(|t| t.0 == f)) else { continue };
        let mult = q.fields.iter().find(|t| t.0 == f).unwrap().1;
        let (meas, sigma) = (row.4[5], row.4[6]);
        let rel = |a: f64, b: f64| (a - b).abs() <= 1e-6 * b.abs().max(1e-12);
        assert!(rel(meas, l[k].2), "{:?}: I meas {meas} against the refined {}", row.2, l[k].2);
        assert!(rel(sigma, mult * r.params[j].esd), "{:?}: σ {sigma} against {} × {}", row.2, mult, r.params[j].esd);
        seen += 1;
    }
    assert_eq!(seen, l.len().min(back.reflections.iter().filter(|x| x.0 == 1).count()));
    // 3. The ties.
    let ties = drs::pawley::ties(&c.model, &c.refinement);
    assert_eq!(ties.len(), 14);
    assert_eq!(back.pawley_ties, ties);
    assert!(back.pawley_ties.iter().any(|t| t.1 == vec![([3, 3, 3], 1.0), ([5, 1, 1], 3.0)]), "{:?}", back.pawley_ties);
    println!("{} reflections with their su's, {} ties, written and read back", seen, back.pawley_ties.len());
}

/// Each intensity of a Pawley phase in the control file (`control::locate`, the Editor's place of a value): its bytes
/// replaced by the value plus 0.125, the file read again moves that intensity and no other; and the file reads back
/// with the same parameters and ties.
#[test]
fn intensities_in_the_control_file() {
    let dir = root("runs/MgAl2O4-LeBail");
    let p = drs::pcr::read(&latin1(&dir.join("MgAl2O4s_LeBail_conv.pcr"))).unwrap();
    let data = Model::from_pcr_with_data(&p, &dir, "MgAl2O4s_LeBail_conv", Some(&dir.join("MgAl2O4s.dat"))).unwrap().data.unwrap();
    let mut c = drs::control::Control::from_pcr(&p, &dir, "MgAl2O4s_LeBail_conv", &dir).unwrap();
    drs::pawley::setup(&mut c.model, &data, &mut c.refinement, 0, 0).unwrap();
    let text = c.to_text();
    let back = drs::control::Control::from_text(&text).unwrap();
    assert_eq!(back.refinement.parameters, c.refinement.parameters);
    assert!(!back.refinement.steering.guards.singular);
    let n = list(&c.model, 0).len();
    for k in 0..n {
        let f = Field::Intensity(0, k);
        let a = drs::refinable::address(f, &c.model);
        let (_, r) = drs::control::locate::locate(&text, &a).unwrap_or_else(|| panic!("{a} not located"));
        let v: f64 = text[r.clone()].parse().unwrap();
        let moved = format!("{}{}{}", &text[..r.start], v + 0.125, &text[r.end..]);
        let m2 = drs::control::Control::from_text(&moved).unwrap().model;
        for j in 0..n {
            let d = m2.get(Field::Intensity(0, j)) - c.model.get(Field::Intensity(0, j));
            assert!((d - if j == k { 0.125 } else { 0.0 }).abs() < 1e-9, "{a}: intensity {j} moved by {d}");
        }
    }
    println!("{n} intensities located");
}
