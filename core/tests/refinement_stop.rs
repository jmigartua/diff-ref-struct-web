//! A refinement that can be stopped (the app's request 6, L070): `lsq::refine_until` stopped after a cycle leaves the
//! model, the agreement factors, the esd's and the correlations exactly where a refinement of that many cycles leaves
//! them, with `stop` = "stopped by the user"; a Stop button's `AtomicBool`, set from another thread, stops it between
//! cycles.

use drs::lsq::{refine, refine_until, Options, STOPPED};
use drs::powder::Model;
use drs::refinable::ParameterSet;
use std::path::Path;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;

fn load() -> (drs::pcr::Pcr, drs::powder::Loaded) {
    // FullProf's starting file of the Y2O3 run.
    let d = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs/Y2O3");
    let text: String = std::fs::read(d.join("y2o3_structure.pcr.in")).unwrap().iter().map(|&b| b as char).collect();
    let p = drs::pcr::read(&text).unwrap();
    let l = Model::from_pcr(&p, &d, "y2o3_structure").unwrap();
    (p, l)
}

#[test]
fn stopped_after_a_cycle() {
    let (p, l) = load();
    let data = l.data.as_ref().unwrap();
    let set = ParameterSet::from_pcr(&p, &l.model).unwrap();
    let mut full = l.model.clone();
    let mut seen = 0;
    let rf = refine(&mut full, data, &set, Options { max_cycles: 6, eps: 0.0 }, |_| seen += 1);
    assert!(rf.cycles >= 3, "the run needs several cycles ({})", rf.cycles);
    for k in 1..rf.cycles {
        let (mut a, mut b) = (l.model.clone(), l.model.clone());
        let ra = refine(&mut a, data, &set, Options { max_cycles: k, eps: 0.0 }, |_| {});
        let mut n = 0;
        let rb = refine_until(&mut b, data, &set, Options { max_cycles: 6, eps: 0.0 }, |c| {
            n += 1;
            assert_eq!(c.cycle, n);
            n < k
        });
        assert_eq!(rb.stop, STOPPED);
        assert!(!rb.converged);
        assert_eq!((rb.cycles, n), (k, k));
        let bits = |m: &Model| m.all_fields().into_iter().map(|f| m.get(f).to_bits()).collect::<Vec<_>>();
        assert_eq!(bits(&a), bits(&b), "after {k} cycles the models differ");
        assert_eq!(ra.stats.chi2.to_bits(), rb.stats.chi2.to_bits());
        assert_eq!(ra.params.iter().map(|q| q.esd.to_bits()).collect::<Vec<_>>(), rb.params.iter().map(|q| q.esd.to_bits()).collect::<Vec<_>>());
        assert_eq!(ra.correlations, rb.correlations);
    }
    println!("Y2O3 from FullProf's start: {} cycles run (eps 0); stopped after each of 1..{} the same as that many cycles", rf.cycles, rf.cycles - 1);
}

#[test]
fn stop_button_from_another_thread() {
    let (p, l) = load();
    let set = ParameterSet::from_pcr(&p, &l.model).unwrap();
    let stop = Arc::new(AtomicBool::new(false));
    let (tx, rx) = std::sync::mpsc::channel::<usize>();
    let flag = stop.clone();
    let worker = std::thread::spawn(move || {
        let mut m = l.model.clone();
        let data = l.data.unwrap();
        refine_until(&mut m, &data, &set, Options { max_cycles: 50, eps: 0.0 }, |c| {
            tx.send(c.cycle).unwrap();
            // The window's thread sets the flag after the first cycle; the engine reads it between cycles.
            while c.cycle == 1 && !flag.load(Ordering::SeqCst) {
                std::thread::yield_now();
            }
            !flag.load(Ordering::SeqCst)
        })
    });
    assert_eq!(rx.recv().unwrap(), 1);
    stop.store(true, Ordering::SeqCst);
    let r = worker.join().unwrap();
    assert_eq!((r.cycles, r.stop.as_str()), (1, STOPPED));
}
