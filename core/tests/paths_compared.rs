//! Explore's E11, E11a, E11b and E11c (L184): two refinements' paths compared, each run's ellipse at its end for a
//! pair, the final correlations and the run's stop kept in the results file. On the five Y₂O₃ runs of
//! `explore/compare/cases/` (README there: A and B reach the same minimum from two starts, C fixes the zero shift, D is
//! A cut at 4 cycles, E diverges), each run by a session here from its control file:
//! - A vs A: everything common, every difference exactly 0, `a_scaled` equal to A's `in_esd().from_end` to the bit;
//! - A vs B: 72 common, 11 steps with B held from step 8, every end within `ENDS_AGREE`;
//! - A vs C: 71 common, `pattern1.zero_shift` only in A, a, the scale and W apart at the ends;
//! - A vs D: D's 5 states A's first 5 to the bit, the difference 0 up to step 4 and not after;
//! - swap: B against A, B the reference, the difference negated;
//! - A vs E: E's stop `diverging-jump` in the summary;
//! - E11a: the pair's ellipses from each run's final correlations;
//! - E11b, E11c: the results file keeps the 2556 correlations and the stop, and a comparison from two results files
//!   read back gives the comparison of the two sessions;
//! - `drs compare`, with `--json` and `--pair`.

use drs::lsq::{self, Options, PathsCompared, Refinement, Run, Session, Side, Trajectory};
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

fn root(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)
}

/// The cases copied beside `y2o3.dat` (the README's instruction; nothing written into `manuals/` or `explore/`).
fn cases() -> &'static PathBuf {
    static DIR: OnceLock<PathBuf> = OnceLock::new();
    DIR.get_or_init(|| {
        let d = root("target/test-paths-compared");
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        for f in ["y2o3_a.drc", "y2o3_b.drc", "y2o3_c.drc", "y2o3_e.drc"] {
            std::fs::copy(root("../explore/compare/cases").join(f), d.join(f)).unwrap();
        }
        std::fs::copy(root("../manuals/fullprof/runs/Y2O3/y2o3.dat"), d.join("y2o3.dat")).unwrap();
        d
    })
}

struct Done {
    path: Trajectory,
    result: Refinement,
    control: drs::control::Control,
    data: drs::data::Pattern,
}

fn refine(file: &str, cycles: Option<usize>) -> Done {
    let o = drs::control::open(&cases().join(file)).unwrap();
    let set = o.control.refinement.build(&o.control.model).unwrap();
    let data = o.data.unwrap();
    let opts = Options { max_cycles: cycles.or(o.control.refinement.cycles).unwrap_or(10), eps: o.control.refinement.eps.unwrap_or(0.1) };
    let mut s = Session::new(o.control.model.clone(), data.clone(), set, opts);
    s.run(|_| true);
    let mut control = o.control;
    control.model = s.model.clone();
    Done { path: s.trajectory(), result: s.finish(), control, data }
}

/// The five runs, once for all the tests.
fn runs() -> &'static [Done; 5] {
    static RUNS: OnceLock<[Done; 5]> = OnceLock::new();
    RUNS.get_or_init(|| [refine("y2o3_a.drc", None), refine("y2o3_b.drc", None), refine("y2o3_c.drc", None), refine("y2o3_a.drc", Some(4)), refine("y2o3_e.drc", None)])
}

struct Ends {
    esd: Vec<f64>,
    code: String,
}

fn ends(d: &Done) -> Ends {
    Ends { esd: d.result.params.iter().map(|p| p.esd).collect(), code: d.result.reason.code().to_string() }
}

fn run<'a>(d: &'a Done, e: &'a Ends) -> Run<'a> {
    Run { path: &d.path, end_esd: Some(&e.esd), end_correlations: Some(&d.result.correlations), stop: Some((e.code.as_str(), d.result.stop.as_str())) }
}

fn compared(x: usize, y: usize, reference: Side) -> PathsCompared {
    let r = runs();
    let (ex, ey) = (ends(&r[x]), ends(&r[y]));
    lsq::compare(&run(&r[x], &ex), &run(&r[y], &ey), reference)
}

#[test]
fn the_runs_are_the_readmes() {
    // The README's table (explore/compare/cases/README.md): cycles, stop, χ² at the start and the end, a at the end. Since
    // the asymmetry is taken at the Bragg angle (L367; the cases' zero is refined with P1, P2) the χ²'s move in their
    // fourth or fifth figure, E's end from 2.0011 to 1.9928: the numbers here are the core's now, the README's to follow.
    let r = runs();
    let expect = [(10, "converged", 82.776, 1.1985, 10.60497), (7, "converged", 29.334, 1.1985, 10.60497), (8, "converged", 5.432, 2.4334, 10.60104), (4, "cycles", 82.776, 5.5113, f64::NAN), (5, "diverging-jump", 105.119, 1.9928, 10.60497)];
    for (k, (d, (n, code, c0, c1, a))) in r.iter().zip(expect).enumerate() {
        let a_end = d.result.params.iter().find(|p| p.name == "phase1.cell_a").unwrap().value;
        println!("run {k}: {} cycles, {}, χ² {} → {}, a {a_end}, {} free", d.path.cycles.len(), d.result.reason.code(), d.path.cycles[0].chi2_start, d.result.stats.chi2, d.path.set.free.len());
        assert_eq!((d.path.cycles.len(), d.result.reason.code()), (n, code), "run {k}");
        assert!((d.path.cycles[0].chi2_start - c0).abs() < 1e-3 && (d.path.cycles.last().unwrap().stats.chi2 - c1).abs() < 1e-4, "run {k}");
        assert!(a.is_nan() || (a_end - a).abs() < 1e-5, "run {k}: a {a_end}");
    }
}

#[test]
fn a_against_a() {
    let c = compared(0, 0, Side::A);
    assert_eq!((c.common.len(), c.only_a.len(), c.only_b.len(), c.fixed_differences.len()), (72, 0, 0, 0));
    assert!(c.difference.iter().flatten().all(|d| *d == Some(0.0)));
    assert!(c.end.iter().all(|e| e.delta == 0.0 && e.agree));
    let r = &runs()[0];
    let own = r.path.in_esd(Some(&ends(r).esd));
    assert_eq!(c.a_scaled, own.from_end);
    assert_eq!(c.b_scaled, own.from_end);
    assert_eq!(c.units[c.common.iter().position(|n| n == "phase1.cell_a").unwrap()], "Å");
}

#[test]
fn a_against_b() {
    let c = compared(0, 1, Side::A);
    assert_eq!((c.common.len(), c.only_a.len(), c.only_b.len(), c.fixed_differences.len()), (72, 0, 0, 0));
    assert_eq!((c.states_a.len(), c.states_b.len(), c.steps.len()), (11, 8, 11));
    assert_eq!(c.steps[10], (10, 7));
    assert_eq!(c.held_b, (0..11).map(|k| k >= 8).collect::<Vec<_>>());
    assert!(c.held_a.iter().all(|h| !h));
    assert_eq!((c.summary.agree, c.summary.disagree), (72, 0));
    let e = |n: &str| c.end.iter().find(|e| e.address == n).unwrap().clone();
    println!("a: {:?}\nzero: {:?}\nlargest {:?}", e("phase1.cell_a"), e("pattern1.zero_shift"), c.summary.largest);
    // The README: B − A for a about 1.3·10⁻⁷ Å, for 2θ₀ about 1.5·10⁻⁶°; 2.1·10⁻⁷ Å and 2.2·10⁻⁶° since L367.
    assert!((e("phase1.cell_a").delta - 2.1e-7).abs() < 0.1e-7 && (e("pattern1.zero_shift").delta - 2.2e-6).abs() < 0.1e-6);
    assert_eq!(c.summary.stop_b.as_ref().unwrap().0, "converged");
}

#[test]
fn a_against_c() {
    let c = compared(0, 2, Side::A);
    assert_eq!((c.common.len(), c.only_b.len(), c.fixed_differences.len()), (71, 0, 0));
    assert_eq!(c.only_a.iter().map(|o| o.0.as_str()).collect::<Vec<_>>(), ["pattern1.zero_shift"]);
    assert_eq!(c.only_a[0].1, 0.05);
    let e = |n: &str| c.end.iter().find(|e| e.address == n).unwrap().clone();
    for n in ["phase1.cell_a", "phase1.scale_factor", "phase1.caglioti_w"] {
        assert!(!e(n).agree, "{n}: {:?}", e(n));
    }
    assert!(e("phase1.cell_a").combined.unwrap().abs() > 50.0);
    println!("{} of 71 ends disagree; largest {:?}", c.summary.disagree, c.summary.largest);
}

#[test]
fn a_against_d() {
    let c = compared(0, 3, Side::A);
    assert_eq!(c.states_b.len(), 5);
    assert_eq!(&c.values_a[..5], &c.values_b[..]);
    for (k, row) in c.difference.iter().enumerate() {
        assert_eq!(row.iter().all(|d| *d == Some(0.0)), k <= 4, "step {k}");
    }
    assert_eq!(c.held_b, (0..11).map(|k| k >= 5).collect::<Vec<_>>());
    assert_eq!(c.summary.stop_b.as_ref().unwrap().0, "cycles");
}

#[test]
fn swapped() {
    let (ab, ba) = (compared(0, 1, Side::A), compared(1, 0, Side::B));
    // B's order is A's here; the reference the same run, so the same scale and the difference negated exactly.
    assert_eq!(ba.common, compared(1, 1, Side::A).common);
    assert_eq!(ab.scale, ba.scale);
    for (x, y) in ab.difference.iter().zip(&ba.difference) {
        assert!(x.iter().zip(y).all(|(p, q)| *p == q.map(|v| -v)));
    }
    for (x, y) in ab.end.iter().zip(&ba.end) {
        assert_eq!(x.combined, y.combined.map(|v| -v));
    }
}

#[test]
fn a_against_e() {
    let c = compared(0, 4, Side::A);
    assert_eq!(c.summary.stop_b.as_ref().unwrap().0, "diverging-jump");
    assert_eq!(c.states_b.len(), 6);
    println!("E: {}; ends {} agree, {} disagree", c.summary.stop_b.as_ref().unwrap().1, c.summary.agree, c.summary.disagree);
}

/// E11a: each run's (σa, σb, r) at its end and its ellipse, from its final correlations.
#[test]
fn e11a_the_pair() {
    let c = compared(0, 2, Side::A);
    let p = c.pair("phase1.caglioti_v", "phase1.caglioti_w").unwrap();
    for (d, sig, r, el, path, centre) in [(&runs()[0], p.sigmas_a, p.correlation_a, p.ellipse_a, &p.path_a, p.centre_a), (&runs()[2], p.sigmas_b, p.correlation_b, p.ellipse_b, &p.path_b, p.centre_b)] {
        let j = |n: &str| d.path.set.free.iter().position(|q| q.name == n).unwrap();
        let (jv, jw) = (j("phase1.caglioti_v"), j("phase1.caglioti_w"));
        assert_eq!(sig, (d.result.params[jv].esd, d.result.params[jw].esd));
        assert_eq!(r, Some(d.result.correlations[jv][jw]));
        assert_eq!(el, lsq::covariance_ellipse(sig.0, sig.1, r.unwrap()));
        assert_eq!(path.len(), d.path.cycles.len() + 1);
        assert_eq!(centre, (d.result.params[jv].value, d.result.params[jw].value));
    }
    println!("V–W: r_A {:?}, r_C {:?}", p.correlation_a, p.correlation_b);
    assert!(c.pair("phase1.caglioti_v", "pattern1.zero_shift").is_err(), "the zero shift is free in A only");
}

/// E11b, E11c: the results file keeps the correlations and the stop; two results files read back compare as the runs.
#[test]
fn e11bc_the_results_file() {
    let r = runs();
    let write = |d: &Done, name: &str| -> PathBuf {
        let set = &d.path.set;
        let (text, _) = drs::results::write(&d.control, &d.data, Some(set), Some(&d.result), name).unwrap();
        let f = cases().join(format!("{name} results.cif"));
        std::fs::write(&f, text).unwrap();
        f
    };
    let (fa, fc) = (write(&r[0], "y2o3_a 2026-10-13 12.00.00"), write(&r[2], "y2o3_c 2026-10-13 12.00.00"));
    let back = drs::results::read(&std::fs::read_to_string(&fa).unwrap()).unwrap();
    assert_eq!(back.correlations.len(), 72 * 71 / 2);
    assert_eq!((back.stop_code.as_deref(), back.stop_text.as_deref()), (Some("converged"), Some(r[0].result.stop.as_str())));
    let names: Vec<String> = r[0].path.set.free.iter().map(|p| p.name.clone()).collect();
    assert_eq!(back.correlation_matrix(&names).unwrap(), r[0].result.correlations, "to the bit");
    assert_eq!(back.free, r[0].result.params.iter().map(|p| (p.name.clone(), p.value, p.esd)).collect::<Vec<_>>());
    // The free parameters' esd's are kept beside the fields' su's: the second asymmetry parameter also moves P1 (× 3),
    // so the su of the field P1 is not the esd of the parameter named after it.
    let su = back.refined.iter().find(|f| f.0 == "phase1.berar_baldinozzi.1").unwrap().2;
    let esd = back.free.iter().find(|f| f.0 == "phase1.berar_baldinozzi.1").unwrap().2;
    println!("P1: the field's su {su}, the parameter's esd {esd}");
    assert!((su - esd).abs() > 0.1 * esd);
    // A file written before (its stop as text only; core 0.25.0's, kept when the committed one was written again, L367)
    // still reads.
    let old = drs::results::read(&std::fs::read_to_string(root("tests/data/results/a_modes_results_0.25.0.cif")).unwrap()).unwrap();
    assert_eq!((old.stop_code, old.stop_text.as_deref(), old.correlations.len()), (None, Some("1 cycles"), 0));
    // `drs compare` on the two results files, against the sessions' comparison. A run read back has its start rebuilt
    // from its end by taking the cycles' shifts back (`Trajectory::from_end`), and subtraction does not undo addition
    // exactly: the start rebuilt may differ from the session's own in the last bits, by an amount the platform's libm
    // decides through the run's numbers (here 0.05 and 0.05000000000000001 for the zero shift on macOS, L187; equal by
    // chance on Linux). So the sessions are read back the same way, from their ends and cycles: one computation on
    // both sides, which the file holds to the bit (the control file's floats and the cycles' loop are exact).
    let o = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(["compare", fa.to_str().unwrap(), fc.to_str().unwrap(), "--json"]).output().unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let files: PathsCompared = serde_json::from_slice(&o.stdout).unwrap();
    let read_back = |d: &Done| Done { path: Trajectory::from_end(&d.control.model, &d.path.set, d.path.cycles.clone()), result: d.result.clone(), control: d.control.clone(), data: d.data.clone() };
    let (ra, rc) = (read_back(&r[0]), read_back(&r[2]));
    let (ea, ec) = (ends(&ra), ends(&rc));
    let sessions = lsq::compare(&run(&ra, &ea), &run(&rc, &ec), Side::A);
    let rebuilt = |d: &Done| d.path.set.free.iter().filter(|p| p.value_in(&d.path.start) != p.value_in(&read_back(d).path.start)).count();
    println!("starts rebuilt from the end that differ from the session's own in some bit: A {}, C {}", rebuilt(&r[0]), rebuilt(&r[2]));
    assert_eq!((&files.common, &files.only_a, &files.states_a, &files.states_b), (&sessions.common, &sessions.only_a, &sessions.states_a, &sessions.states_b));
    assert_eq!(files.values_a.last(), sessions.values_a.last());
    assert_eq!((files.summary.agree, &files.summary.stop_a), (sessions.summary.agree, &sessions.summary.stop_a));
    for (x, y) in files.end.iter().zip(&sessions.end) {
        assert_eq!((x.delta, x.combined), (y.delta, y.combined), "{}", x.address);
    }
    assert_eq!((files.end_esd_a, files.end_correlations_a), (sessions.end_esd_a, sessions.end_correlations_a));
    // The command line on two control files, with a pair.
    let o = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(["compare", cases().join("y2o3_a.drc").to_str().unwrap(), cases().join("y2o3_b.drc").to_str().unwrap(), "--pair", "phase1.caglioti_v", "phase1.caglioti_w"]).output().unwrap();
    let out = String::from_utf8_lossy(&o.stdout);
    assert!(o.status.success() && out.contains("72 in common, 0 only in a, 0 only in b, 0 fixed values changed; ends: 72 agree") && out.contains("pair phase1.caglioti_v – phase1.caglioti_w"), "{out}");
}
