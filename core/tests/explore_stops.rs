//! Explore's E15 (L225): W11's stops and pauses kept in the run's record, on the Y₂O₃ runs of
//! `explore/compare/cases/` (copied beside `y2o3.dat`, as their README says):
//! - E15a: a run of 6 cycles stepped 2 at a time pauses after cycles 2 and 4 (by `step`), a value set in the second
//!   pause makes it changed with λ restarted, a caller's Stop becomes a pause (`Session::pause`), and the session file
//!   and the results file give the pauses back;
//! - E15b: case E stops diverging-jump on phase1.cell_a at cycle 5, kept in the trajectory and read back from its
//!   results file with its facts;
//! - E15c: the steering's notes as data: case A's first cycles shorten the pseudo-Voigt η's steps to their maximum
//!   shifts, as the words say, and a refused step carries its λ and χ².

use drs::lsq::{Change, Divergence, NoteKind, Options, PauseWhy, Session, Stop};
use std::path::{Path, PathBuf};

fn root(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)
}

/// Each test's own copy of the cases (L233: tests run in parallel, and a file copied while another test reads it is
/// read half written).
fn cases(test: &str) -> PathBuf {
    let d = root("target/test-explore-stops").join(test);
    let _ = std::fs::create_dir_all(&d);
    for f in ["y2o3_a.drc", "y2o3_e.drc"] {
        std::fs::copy(root("../explore/compare/cases").join(f), d.join(f)).unwrap();
    }
    std::fs::copy(root("../manuals/fullprof/runs/Y2O3/y2o3.dat"), d.join("y2o3.dat")).unwrap();
    d
}

fn session(test: &str, file: &str, cycles: usize) -> (Session, drs::control::Control) {
    let o = drs::control::open(&cases(test).join(file)).unwrap();
    let set = o.control.refinement.build(&o.control.model).unwrap();
    (Session::new(o.control.model.clone(), o.data.unwrap(), set, Options { max_cycles: cycles, eps: 0.1 }), o.control)
}

fn results_back(s: &Session, control: &drs::control::Control) -> drs::results::Results {
    let r = s.finish();
    let mut c = control.clone();
    c.model = s.model.clone();
    let (text, _) = drs::results::write(&c, &s.data, Some(&s.set), Some(&r), "run 2026-10-14 12.00.00").unwrap();
    drs::results::read(&text).unwrap()
}

#[test]
fn e15a_the_pauses_in_the_record() {
    let (mut s, control) = session("e15a", "y2o3_a.drc", 6);
    let lambda0 = s.lambda;
    s.step(2, |_| true);
    s.step(2, |_| true);
    s.set_value("pattern1.zero_shift", 0.05).unwrap();
    s.step(2, |_| true);
    assert_eq!(s.cycle, 6);
    // Live, the third step's end is a pause too; the finished run's record ends at cycle 6 instead.
    let t = s.trajectory();
    assert_eq!(t.pauses.iter().map(|p| (p.after, p.why)).collect::<Vec<_>>(), vec![(2, PauseWhy::Step), (4, PauseWhy::Step), (6, PauseWhy::Step)]);
    let r = s.finish();
    assert_eq!(r.pauses.len(), 2);
    let (p2, p4) = (&r.pauses[0], &r.pauses[1]);
    assert!(p2.changes.is_empty() && p2.restarted.is_none());
    assert_eq!(p4.changes, vec![Change::Value { address: "pattern1.zero_shift".into(), value: 0.05 }]);
    assert_eq!(p4.restarted, Some(lambda0), "λ restarted");
    assert!(p4.lambda.is_finite());
    // The session file and the results file give them back.
    let text = s.to_text();
    assert!(text.contains("pause 2 step") && text.contains("change 4 value pattern1.zero_shift 0.05"), "{text}");
    let back = Session::resume(&text, s.start.clone(), s.data.clone(), s.set.clone(), s.opts).unwrap();
    assert_eq!(back.pauses, s.pauses);
    let rb = results_back(&s, &control);
    assert_eq!(rb.pauses, r.pauses);
    // A caller's Stop (the user's Pause) becomes a pause: the run can go on.
    let (mut c, _) = session("e15a caller", "y2o3_a.drc", 6);
    c.step(5, |cy| cy.cycle < 3);
    assert_eq!(c.stop, Some(Stop::User));
    c.pause();
    assert!(c.stop.is_none());
    assert_eq!(c.pauses.iter().map(|p| (p.after, p.why)).collect::<Vec<_>>(), vec![(3, PauseWhy::Caller)]);
    c.step(1, |_| true);
    assert_eq!(c.cycle, 4);
}

#[test]
fn e15b_the_stops_facts() {
    let (mut s, control) = session("e15b", "y2o3_e.drc", 10);
    s.run(|_| true);
    let st = s.trajectory().stop.unwrap();
    let Stop::Diverging(Divergence::Jump { parameter, ratio, cycle }) = &st.stop else { panic!("{st:?}") };
    assert_eq!((parameter.as_str(), *cycle, st.after), ("phase1.cell_a", 5, 5));
    assert!(*ratio > 100.0);
    let rb = results_back(&s, &control);
    assert_eq!(rb.stop.as_ref(), Some(&st), "read back with its facts");
    assert_eq!(rb.stop_code.as_deref(), Some("diverging-jump"));
    // Every stop's facts go and come back.
    for stop in [
        Stop::Converged { eps: 0.1 },
        Stop::Minimum,
        Stop::Cycles { n: 7, le_bail: true },
        Stop::Diverging(Divergence::Rising { cycles: 3, from: 1.5, to: 2.25 }),
        Stop::Diverging(Divergence::Impossible { what: "phase1.cell_a".into(), value: -1.0 }),
        Stop::Singular { combination: vec![("phase1.caglioti_u".into(), 0.7), ("phase1.caglioti_w".into(), -0.7)], ratio: 1e-17, why: String::new(), tip: String::new() },
        Stop::User,
    ] {
        let f = stop.facts();
        assert_eq!(Stop::from_facts(stop.code(), |k| f.iter().find(|x| x.0 == k).map(|x| x.1.clone())), Some(stop.clone()), "{stop:?}");
    }
}

#[test]
fn e15c_the_steerings_notes_as_data() {
    let (mut s, _) = session("e15c", "y2o3_a.drc", 3);
    let mut cycles = Vec::new();
    s.run(|c| {
        cycles.push(c.clone());
        true
    });
    // Cycles 1 and 2: the η's steps shortened to their maximum shifts, as the words say. Since the asymmetry is taken at
    // the Bragg angle (L367; this run's zero is 0.05° and its P1, P2 refined), η0's first step is shortened too.
    let max = |k: usize| cycles[k].steering.iter().filter(|n| n.kind == NoteKind::MaximumShift).map(|n| n.address.clone().unwrap()).collect::<Vec<_>>();
    assert_eq!(max(0), ["phase1.pv_eta0", "phase1.pv_eta_slope"]);
    assert_eq!(max(1), ["phase1.pv_eta0", "phase1.pv_eta_slope"]);
    for c in &cycles {
        for n in c.steering.iter().filter(|n| n.kind == NoteKind::MaximumShift) {
            let a = n.address.as_ref().unwrap();
            assert!(c.notes.iter().any(|w| w.starts_with(a.as_str()) && w.contains("maximum shift")), "{a}");
            assert!(n.value.unwrap() > 0.0);
        }
        // A refused step: its λ and χ², as `tries` has them.
        let refused: Vec<(f64, f64)> = c.steering.iter().filter(|n| n.kind == NoteKind::Refused).map(|n| (n.lambda.unwrap(), n.chi2.unwrap())).collect();
        assert_eq!(refused, c.tries.iter().map(|t| (t.lambda, t.chi2)).collect::<Vec<_>>());
    }
    assert_eq!(cycles[2].steering.iter().filter(|n| n.kind == NoteKind::Refused).count(), 3);
    assert_eq!(NoteKind::MaximumShift.word(), "maximum-shift");
    // The results file gives them back, cycle by cycle.
    let (_, control) = session("e15c", "y2o3_a.drc", 3);
    let rb = results_back(&s, &control);
    assert_eq!(rb.cycles.iter().map(|c| c.steering.clone()).collect::<Vec<_>>(), cycles.iter().map(|c| c.steering.clone()).collect::<Vec<_>>());
}

/// The command-line path: `drs refine --trace` writes the notes as data, the stop's facts and its cycle, and a
/// session stepped by `--session --cycles` keeps its pauses.
#[test]
fn e15_on_the_command_line() {
    let d = cases("e15 cli");
    let drs = |args: &[&str]| {
        let o = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(args).output().unwrap();
        assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
        String::from_utf8_lossy(&o.stdout).to_string()
    };
    let (e, trace) = (d.join("y2o3_e.drc"), d.join("e15 e.json"));
    drs(&["refine", e.to_str().unwrap(), "--cycles", "10", "--trace", trace.to_str().unwrap()]);
    let json = std::fs::read_to_string(&trace).unwrap();
    assert!(json.contains("\"code\": \"diverging-jump\"") && json.contains("\"after\": 5, \"parameter\": \"phase1.cell_a\""), "{}", &json[..600.min(json.len())]);
    assert!(json.contains("\"cycle\": \"5\""));
    let (a, trace, sess) = (d.join("y2o3_a.drc"), d.join("e15 a.json"), d.join("e15 a.session"));
    let _ = std::fs::remove_file(&sess);
    let step = || drs(&["refine", a.to_str().unwrap(), "--session", sess.to_str().unwrap(), "--cycles", "2", "--trace", trace.to_str().unwrap()]);
    step();
    let json = std::fs::read_to_string(&trace).unwrap();
    assert!(json.contains("{\"kind\": \"maximum-shift\", \"address\": \"phase1.pv_eta_slope\""), "{}", &json[..600.min(json.len())]);
    // The session file keeps the cycles' numbers, not their notes: the second step's trace has the pauses.
    step();
    let json = std::fs::read_to_string(&trace).unwrap();
    assert!(json.contains("\"pauses\": [{\"after\": 2, \"why\": \"step\""), "{}", &json[json.len().saturating_sub(600)..]);
    assert!(std::fs::read_to_string(&sess).unwrap().contains("pause 2 step"));
}
