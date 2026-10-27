//! FullProf's step on the cell (L111–L113): FullProf moves a cell parameter by a fraction of the Gauss–Newton step,
//! the other parameters by the whole of it. Replayed in the core as a damping of the cell's parameters alone
//! (`steering`), the fraction reproduces FullProf's cell cycle after cycle:
//! - ½ with the pseudo-Voigt (Y2O3, Npr = 5): near the minimum (`shift_control/n3`, a, b, c 3 esd off) FullProf
//!   recovers 1.44 esd where the Gauss–Newton step recovers 2.95, and every other parameter moves as the core's,
//!   within 0.1 esd; far from it (`m10`, a 130 esd off) the core's half step lands on FullProf's a;
//! - ¼ with Thompson–Cox–Hastings (Npr = 7): LaB6 over 9 cycles and rutile–anatase over 12, to FullProf's print.
//!   There FullProf also prints the cell's esd at half of √(C_kk χ²) (L053, `fullprof_refinements`): it behaves as
//!   if the cell's derivative were doubled, which halves the step again;
//! - Egypt's five TCH phases, refined together, are not reproduced by one fraction for all their cells: counted.
//!
//! Not A = 1/a²: a Gauss–Newton step is the same to first order in a and in A, and the codeword's multiplier is 1;
//! the factor is FullProf's own and not in its manual. Counted and frozen here, listed for Josu (L113).

use drs::jacobian::Field;
use drs::lsq::{Options, Session};
use drs::powder::Model;
use drs::refinable::ParameterSet;
use std::path::Path;

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display())).iter().map(|&b| b as char).collect()
}

fn runs() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs")
}

fn load(d: &Path, file: &str, code: &str, data: Option<&Path>) -> (drs::pcr::Pcr, drs::powder::Loaded) {
    let mut p = drs::pcr::read(&latin1(&d.join(file))).unwrap();
    p.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
    let l = Model::from_pcr_with_data(&p, d, code, data).unwrap_or_else(|e| panic!("{file}: {e}"));
    (p, l)
}

/// A set with the cell's parameters damped by `f`.
fn cell_damped(set: &ParameterSet, f: f64) -> ParameterSet {
    let mut s = set.clone();
    for q in s.free.iter_mut().filter(|q| matches!(q.fields[0].0, Field::Cell(..))) {
        q.damping = f;
    }
    s
}

/// One cycle from `m`: the model after it.
fn one_cycle(m: &Model, data: &drs::data::Pattern, set: &ParameterSet) -> Model {
    let mut s = Session::new(m.clone(), data.clone(), set.clone(), Options { max_cycles: 1, eps: 0.0 });
    s.step(1, |_| true);
    s.model
}

/// FullProf's parameters in a .out by number: name, start, end, esd.
fn fullprof_parameters(out: &str) -> std::collections::BTreeMap<u32, (String, f64, f64, f64)> {
    let mut st = std::collections::BTreeMap::new();
    let mut en = std::collections::BTreeMap::new();
    for l in out.lines() {
        if let Some((a, r)) = l.split_once("-> Symbolic Name:") {
            let n: u32 = a.split_whitespace().last().unwrap().parse().unwrap();
            let w: Vec<&str> = r.split_whitespace().collect();
            st.entry(n).or_insert((w[0].to_string(), w[1].parse::<f64>().unwrap()));
        } else if l.contains("Parameter number") && l.contains("( +/-") {
            let (a, b) = l.split_once(':').unwrap();
            let n: u32 = a.split_whitespace().last().unwrap().parse().unwrap();
            let (v, e) = b.split_once("( +/-").unwrap();
            en.insert(n, (v.split_whitespace().last().unwrap().parse::<f64>().unwrap(), e.trim().trim_end_matches(')').trim().parse::<f64>().unwrap()));
        }
    }
    st.into_iter().map(|(n, (name, s))| (n, (name, s, en[&n].0, en[&n].1))).collect()
}

/// FullProf's cell values per cycle and phase.
fn fullprof_cells(out: &str) -> Vec<Vec<[f64; 6]>> {
    let l: Vec<&str> = out.lines().collect();
    let mut cyc: Vec<Vec<[f64; 6]>> = Vec::new();
    for i in 0..l.len() {
        if l[i].contains("=> CYCLE No.:") {
            cyc.push(Vec::new());
        }
        if l[i].contains("=> Cell parameters") && !cyc.is_empty() {
            let c: [f64; 6] = std::array::from_fn(|k| l[i + 1 + k].split_whitespace().next().unwrap().parse().unwrap());
            cyc.last_mut().unwrap().push(c);
        }
    }
    cyc
}

/// The frozen sizes (L113): n3's fraction and the other parameters' largest difference (esd's of FullProf); the
/// replays' largest deviation from FullProf's cell (Å); Egypt's cell values better fitted by ½ and by ¼.
const N3_FRACTION: (f64, f64) = (0.47, 0.51);
const N3_OTHERS: f64 = 0.089;
const M10_HALF: f64 = 5e-6;
const REPLAY: &[(&str, f64)] = &[("LaB6", 2.8e-6), ("rutile-anatase", 2.7e-5)];
const EGYPT: (usize, usize) = (10, 5);

#[test]
fn near_the_minimum_half_a_step() {
    let d = runs().join("Y2O3");
    let (p, l) = load(&d, "y2o3_structure.pcr.in", "y2o3_structure", Some(&d.join("y2o3.dat")));
    let data = l.data.clone().unwrap();
    let fp = fullprof_parameters(&latin1(&runs().join("shift_control/n3.out")));
    // n3's start: the tutorial's file (FullProf's minimum: its first cycle moves nothing by 10⁻⁶) with a, b, c at
    // 10.604829. The committed n3.pcr.in holds FullProf's values after the cycle; n3.out its start.
    let mut moved = l.model.clone();
    for k in 0..3 {
        moved.set(Field::Cell(0, k), 10.604829);
    }
    let set = ParameterSet::from_pcr(&p, &l.model).unwrap();
    let numbers = drs::refinable::pcr_parameter_numbers(&p, &l.model).unwrap();
    // Each program's step from n3 less its step from the unmoved file (the core's minimum is not exactly FullProf's).
    let rest = one_cycle(&l.model, &data, &set);
    let full = one_cycle(&moved, &data, &set);
    let half = one_cycle(&moved, &data, &cell_damped(&set, 0.5));
    let mut others: f64 = 0.0;
    let (mut frac, mut cell_half) = (0.0, 0.0);
    let mut line = Vec::new();
    for (n, (name, s0, s1, esd)) in &fp {
        if name.starts_with("Bck") || name.starts_with("Scale") || name.starts_with("Asym") {
            continue;
        }
        let f = numbers[n][0];
        let fp_step = (s1 - s0) / esd;
        let core = ((full.get(f) - moved.get(f)) - (rest.get(f) - l.model.get(f))) / esd;
        if f == Field::Cell(0, 0) {
            frac = fp_step / core;
            cell_half = ((half.get(f) - moved.get(f)) - 0.5 * (rest.get(f) - l.model.get(f))) / esd - fp_step;
        } else {
            others = others.max((core - fp_step).abs());
        }
        line.push(format!("{name} {fp_step:+.2}/{core:+.2}"));
    }
    println!("n3, FullProf/core in esd: {}; the cell's fraction {frac:.3}, the half step off by {cell_half:+.3} esd, the others within {others:.3} esd", line.join(", "));
    assert!((N3_FRACTION.0..=N3_FRACTION.1).contains(&frac), "{frac}");
    assert!(cell_half.abs() < 0.1 && others <= N3_OTHERS, "{cell_half} {others}");
    // Far from the minimum (m10, a 130 esd off): the core's half step lands on FullProf's a.
    let ds = runs().join("shift_control");
    let (pm, lm) = load(&ds, "m10.pcr.in", "m10", Some(&d.join("y2o3.dat")));
    let sm = ParameterSet::from_pcr(&pm, &lm.model).unwrap();
    let a = one_cycle(&lm.model, lm.data.as_ref().unwrap(), &cell_damped(&sm, 0.5)).get(Field::Cell(0, 0));
    println!("m10: the core's half step gives a = {a:.6}, FullProf 10.601567");
    assert!((a - 10.601567).abs() < M10_HALF);
}

#[test]
fn replayed_cycle_by_cycle() {
    let mut egypt = (0, 0);
    for (run, code, smooth) in [("LaB6", "LaB6rx", false), ("rutile-anatase", "rutana_tch", true), ("Egypt", "Egypt_qa_Rietveld", false)] {
        let d = runs().join(run);
        let (_, conv) = load(&d, &format!("{code}_conv.pcr"), code, None);
        let (p, start) = load(&d, &format!("{code}.pcr.in"), code, conv.data_path.as_deref());
        let out = latin1(&d.join(format!("{code}_conv.out")));
        let fp = fullprof_cells(&out);
        // As `fullprof_refinements` takes these runs.
        let mut m = start.model.clone();
        let lim: Vec<f64> = out.lines().filter_map(|l| l.split_once("S(1/d)max:")).map(|(_, r)| r.split_whitespace().next().unwrap().parse().unwrap()).collect();
        if lim.len() >= m.phases.len() {
            for (ph, l) in m.phases.iter_mut().zip(&lim) {
                ph.reflection_limit = Some(*l);
            }
        }
        if smooth {
            m.convention.peak_window = false;
        }
        let set = ParameterSet::from_pcr(&p, &m).unwrap();
        let np = m.phases.len();
        // The deviation of the core's cell from FullProf's at every cycle, per phase and parameter, by fraction.
        let dev: Vec<Vec<Vec<f64>>> = [1.0, 0.25, 0.5]
            .iter()
            .map(|&f| {
                let mut s = Session::new(m.clone(), start.data.clone().unwrap(), cell_damped(&set, f), Options { max_cycles: fp.len(), eps: 0.0 });
                let mut out = vec![Vec::new(); np * 6];
                for c in &fp {
                    s.step(1, |_| true);
                    for ph in 0..np {
                        for k in 0..6 {
                            out[ph * 6 + k].push(s.model.get(Field::Cell(ph, k)) - c[ph][k]);
                        }
                    }
                }
                out
            })
            .collect();
        let moved = |i: usize| fp.iter().any(|c| c[i / 6][i % 6] != fp[0][i / 6][i % 6]) || dev[0][i].iter().any(|x| x.abs() > 1e-5);
        let max = |f: usize, i: usize| dev[f][i].iter().fold(0.0f64, |a, x| a.max(x.abs()));
        let rms = |f: usize, i: usize| (dev[f][i].iter().map(|x| x * x).sum::<f64>() / dev[f][i].len() as f64).sqrt();
        if run == "Egypt" {
            for i in (0..np * 6).filter(|&i| moved(i)) {
                if rms(2, i) < rms(1, i) {
                    egypt.0 += 1;
                } else {
                    egypt.1 += 1;
                }
            }
            continue;
        }
        let worst = (0..np * 6).filter(|&i| moved(i)).map(|i| max(1, i)).fold(0.0, f64::max);
        let whole = (0..np * 6).filter(|&i| moved(i)).map(|i| max(0, i)).fold(0.0, f64::max);
        println!("{run}: {} cycles, the core with the cell's step at ¼ within {worst:.1e} Å of FullProf's cell at every cycle (the whole step: {whole:.1e})", fp.len());
        let bound = REPLAY.iter().find(|r| r.0 == run).unwrap().1;
        assert!(worst <= bound && whole > 10.0 * worst, "{run}: {worst} {whole}");
    }
    println!("Egypt: of its cell values, {} follow FullProf better with ½, {} with ¼", egypt.0, egypt.1);
    assert!(egypt == EGYPT, "{egypt:?}");
}
