//! The refinement steered and run as a session (W11–W12, L107–L110):
//! - `shift_control_against_fullprof`: FullProf 8.20's one cycle of Y2O3 with the cell's codeword 241.00 (m10), its
//!   multiplier 0.5 (m05) and R_pr = 0.5 (r05) (`manuals/fullprof/runs/shift_control/`): for every refined parameter,
//!   FullProf's shift in m05 and r05 over its shift in m10 is the core's (0.5 where the core damps it, 1 elsewhere),
//!   which checks the import (a lone multiplier is damping, the relaxation factors are the classes') and the classes
//!   (`steering::class_of`) against FullProf on every class present; the core's esd's are not changed by the damping,
//!   and FullProf's printed esd of the multiplied parameter is the core's times the multiplier (`lsq::fullprof_esd`,
//!   FullProf's eq. 3.7).
//! - `multipliers_by_meaning`: what the import makes of a lone multiplier, of shared ones, and of R_at … R_gl.
//! - `damping_never_changes_an_esd`, `max_shift_shortens_the_step`, `bounds_stop_reflect_periodic`: the steering.
//! - `links_with_offset`: Co-YBCO's mixed sites (Cu1/Co1 as p and s − p, equal B's) written as links give the set and
//!   the refinement FullProf's shared codewords give; one value across phases.
//! - `variables_are_mode_amplitudes`: case A's 21 mode amplitudes refined as new variables moving the coordinates
//!   give the refinement of the amplitudes themselves.
//! - `divergence_guards`: each guard stops the run with its reason as data.
//! - `session_pause_resume`: a session stepped, saved, resumed and run on gives the straight run.
//! - `plan_of_steps`: a plan run step by step reaches the minimum of the single run; written and read back.
//! - `restraint_weights_and_contributions`, `control_file_guardrails`, `command_line`.

use drs::control::Control;
use drs::jacobian::Field;
use drs::lsq::{fullprof_esd, Divergence, Options, Session, Stop};
use drs::powder::Model;
use drs::refinable::{field_of, ParameterSet, Refinement};
use drs::steering::{Behaviour, Class, Link, PlanStep, Variable};
use std::path::{Path, PathBuf};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof")
}

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display())).iter().map(|&b| b as char).collect()
}

fn pcr(p: &Path) -> drs::pcr::Pcr {
    let mut p = drs::pcr::read(&latin1(p)).unwrap();
    p.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
    p
}

/// A run of `shift_control`: its .pcr, model and data (Y2O3's data file).
fn shift_run(code: &str) -> (drs::pcr::Pcr, drs::powder::Loaded) {
    let d = root().join("runs/shift_control");
    let p = pcr(&d.join(format!("{code}.pcr.in")));
    let l = Model::from_pcr_with_data(&p, &d, code, Some(&root().join("runs/Y2O3/y2o3.dat"))).unwrap();
    (p, l)
}

/// Y2O3's converged run as a control file, with its data. FullProf's peak windows (a hard cut at Wdt·H, kept on
/// import where no exact translation exists) make χ² piecewise in the values: the default model's smooth window
/// (L041) is taken here, so that χ² is smooth and the minimum unique.
fn y2o3() -> (Control, drs::data::Pattern) {
    let d = root().join("runs/Y2O3");
    let p = pcr(&d.join("y2o3_structure_conv.pcr"));
    let mut c = Control::from_pcr(&p, &d, "y2o3_structure", &d).unwrap();
    c.model.convention.peak_window = false;
    c.model.convention.window_edge = false;
    let data = Model::from_pcr(&p, &d, "y2o3_structure").unwrap().data.unwrap();
    (c, data)
}

/// FullProf's parameters in a `.out`: symbolic name → (start, value after the cycles, esd).
fn fullprof_parameters(out: &str) -> Vec<(String, f64, f64, f64)> {
    let mut start = Vec::new();
    let mut end = std::collections::BTreeMap::new();
    for l in out.lines() {
        if let Some((_, r)) = l.split_once("-> Symbolic Name:").filter(|(_, r)| !start.iter().any(|(n, _): &(String, f64)| Some(n.as_str()) == r.split_whitespace().next())) {
            let w: Vec<&str> = r.split_whitespace().collect();
            start.push((w[0].to_string(), w[1].parse::<f64>().unwrap()));
        } else if l.contains("Parameter number") && l.contains("( +/-") {
            let (a, b) = l.split_once(':').unwrap().1.split_once("( +/-").unwrap();
            let w: Vec<&str> = a.split_whitespace().collect();
            end.insert(w[0].to_string(), (w[1].parse::<f64>().unwrap(), b.trim().trim_end_matches(')').trim().parse::<f64>().unwrap()));
        }
    }
    start.into_iter().map(|(n, s)| (n.clone(), s, end[&n].0, end[&n].1)).collect()
}

/// The first cycle's frozen difference with FullProf (L109, re-read in L113): parameters other than the cell whose
/// shift is more than 10 % off, and the largest factor between the two shifts. The cell's factor 2 is FullProf's half
/// step on the cell (`fullprof_cell_step`).
const FROZEN_COUNT: usize = 5;
const FROZEN_FACTOR: f64 = 4.29;

#[test]
fn shift_control_against_fullprof() {
    let d = root().join("runs/shift_control");
    let mut fails = Vec::new();
    // FullProf: each parameter's shift and esd in the three runs.
    let fp: Vec<Vec<(String, f64, f64, f64)>> = ["m10", "m05", "r05"].iter().map(|c| fullprof_parameters(&latin1(&d.join(format!("{c}.out"))))).collect();
    // The core: one cycle of each run, from the same start.
    let mut ours = Vec::new();
    for code in ["m10", "m05", "r05"] {
        let (p, l) = shift_run(code);
        let set = ParameterSet::from_pcr(&p, &l.model).unwrap();
        let mut s = Session::new(l.model.clone(), l.data.clone().unwrap(), set, Options { max_cycles: 1, eps: p.eps });
        let mut cyc = None;
        s.step(1, |c| {
            cyc = Some(c.clone());
            true
        });
        ours.push((p, s.model.clone(), s.set.clone(), cyc.unwrap()));
    }
    let numbers = drs::refinable::pcr_parameter_numbers(&ours[0].0, &ours[0].1).unwrap();
    assert_eq!(numbers.len(), fp[0].len());
    let (mut n_half, mut n_whole) = (0, 0);
    for (k, (name, start, end, esd)) in fp[0].iter().enumerate() {
        let shift = end - start;
        // FullProf's parameter k + 1 (its list is in order) and the core's parameter that moves its first value.
        let f = numbers[&(k as u32 + 1)][0];
        let j = ours[0].2.free.iter().position(|q| q.fields.iter().any(|x| x.0 == f)).unwrap();
        let cls = drs::steering::class_of(f);
        for (run, idx) in [("m05", 1), ("r05", 2)] {
            let fp_ratio = (fp[idx][k].2 - start) / shift;
            let core_ratio = ours[idx].3.params[j].shift / ours[0].3.params[j].shift;
            let expect = match run {
                "m05" if f == Field::Cell(0, 0) => 0.5,
                "r05" if cls == Class::Profile => 0.5,
                _ => 1.0,
            };
            // FullProf's ratio to its print (6 figures of the value), the core's exactly.
            let tol = 1e-6 * end.abs().max(1e-6) / shift.abs() + 1e-4;
            if (fp_ratio - expect).abs() > tol.max(5e-4) {
                fails.push(format!("{name}: FullProf's {run} shift is {fp_ratio:.4} of m10's, the core expects {expect}"));
            }
            if (core_ratio - expect).abs() > 1e-9 {
                fails.push(format!("{name} ({}): the core's {run} shift is {core_ratio} of m10's, {expect} expected", ours[0].2.free[j].name));
            }
            if expect == 0.5 && run == "r05" {
                n_half += 1;
            } else if run == "r05" {
                n_whole += 1;
            }
            // The computed (undamped) shift is the same in the three runs: convergence is judged on it.
            let (c0, c1) = (ours[0].3.params[j].computed, ours[idx].3.params[j].computed);
            if run == "m05" && (c0 - c1).abs() > 1e-9 * c0.abs().max(1e-300) {
                fails.push(format!("{name}: computed shift {c1} in m05, {c0} in m10"));
            }
        }
        // The esd of cell a: FullProf prints a·σ for the multiplied parameter (m05: 0.49989 of m10's; its σ(a) from the
        // reciprocal metric at the new a, (10.600783/10.601567)³ = 0.99978 of m10's in r05).
        if f == Field::Cell(0, 0) {
            let fp_m05 = fp[1][k].3 / esd;
            let fp_r05 = fp[2][k].3 / esd;
            let cubed = (fp[2][k].2 / end).powi(3);
            if (fp_r05 - cubed).abs() > 2e-5 || (fp_m05 - 0.5 * cubed).abs() > 2e-5 {
                fails.push(format!("FullProf's esd of a: m05/m10 {fp_m05:.5}, r05/m10 {fp_r05:.5}, (a ratio)³ {cubed:.5}"));
            }
            let e: Vec<f64> = ours.iter().map(|o| o.3.params[j].esd).collect();
            let fpe = fullprof_esd(&ours[1].2, j, e[1]);
            println!("cell a after one cycle: FullProf {end:.6} ({esd:.4e}), {:.6} ({:.4e}), {:.6} ({:.4e}); the core {:.6} ({:.4e}), {:.6} ({:.4e}, FullProf's reading {:.4e}), {:.6} ({:.4e})", fp[1][k].2, fp[1][k].3, fp[2][k].2, fp[2][k].3, ours[0].1.get(f), e[0], ours[1].1.get(f), e[1], fpe, ours[2].1.get(f), e[2]);
            // The core's esd σ = √((A⁻¹)_jj χ²_ν) is not changed by the damping: (A⁻¹)_jj is the start's in the three
            // runs, χ²_ν the values reached (a half step leaves it higher). FullProf takes χ²_ν at the start of the
            // cycle (7.54 printed in the three .out): the core's σ at that χ² is FullProf's to its χ² over the points
            // with Bragg contributions (`esd_profile`).
            let unit: Vec<f64> = ours.iter().map(|o| o.3.params[j].esd / o.3.stats.chi2.sqrt()).collect();
            if (unit[1] / unit[0] - 1.0).abs() > 1e-9 || (unit[2] / unit[0] - 1.0).abs() > 1e-9 {
                fails.push(format!("the core's esd of a over √χ²: {unit:?}"));
            }
            let at_start = unit[0] * ours[0].3.chi2_start.sqrt();
            println!("the core's esd of a at the start's χ² ({:.4}): {at_start:.4e}, FullProf's {esd:.4e}", ours[0].3.chi2_start);
            if (at_start / esd - 1.0).abs() > 0.01 {
                fails.push(format!("the core's esd of a at the start's χ²: {at_start:.4e} against FullProf's {esd:.4e}"));
            }
            if (fpe / e[1] - 0.5).abs() > 1e-12 {
                fails.push(format!("fullprof_esd: {fpe} for {}", e[1]));
            }
        }
    }
    println!("{} parameters: {n_half} damped by R_pr = 0.5 in r05 (the profile class: scale, U, V, W, η, X, asymmetry, cell), {n_whole} not (zero, background)", fp[0].len());
    assert!(n_half == 9 && n_whole == 63, "{n_half} {n_whole}");
    // The first cycle itself (m10) is not FullProf's (L109). This start has a 0.005 Å (≈ 130 esd) off the minimum of
    // FullProf's own run (whose first cycle moves nothing by more than 10⁻⁶): at 2θ = 140° the peaks are misplaced by
    // more than their width, far outside the range where a Gauss–Newton step is linear, and the two programs' steps
    // differ parameter by parameter, though both reach the same minimum (`fullprof_refinements`). Counted and frozen
    // here on the parameters translated as they are (the scale and the asymmetry are translated on import, ×2 and
    // P1 + 3·P2; the background points move less than their esd in both programs), and listed for Josu with the run
    // that would separate it (the converged file with a moved by a few esd, one cycle).
    // Re-read with that run (n3, L111–L113): the cell's factor 2.00 is FullProf's: it moves a cell parameter by half a
    // Gauss–Newton step near the minimum and far from it (`fullprof_cell_step`), and is left out here. The six others
    // are the far start's alone: near the minimum they move as the core's, within 0.09 esd.
    let (_, start) = shift_run("m10");
    let mut ratios = Vec::new();
    for (k, (name, s0, end, esd)) in fp[0].iter().enumerate() {
        let f = numbers[&(k as u32 + 1)][0];
        if matches!(f, Field::Scale(_) | Field::Asym(..) | Field::Cell(..)) {
            continue;
        }
        if matches!(f, Field::Background(_)) {
            assert!((end - s0).abs() < 1.2 * esd, "{name}: FullProf moves it {} esd", (end - s0) / esd);
            continue;
        }
        ratios.push((name.clone(), (ours[0].1.get(f) - start.model.get(f)) / (end - s0)));
    }
    let off = ratios.iter().filter(|r| (r.1 - 1.0).abs() > 0.1).count();
    let worst = ratios.iter().map(|r| if r.1 > 0.0 { r.1.max(1.0 / r.1) } else { f64::INFINITY }).fold(0.0, f64::max);
    let list: Vec<String> = ratios.iter().map(|r| format!("{} {:.3}", r.0, r.1)).collect();
    println!("the first cycle, the core's shift over FullProf's: {off} of {} beyond 10 %, the largest factor {worst:.2}: {}", ratios.len(), list.join(", "));
    assert!(off <= FROZEN_COUNT && worst <= FROZEN_FACTOR, "the first cycle's difference with FullProf grew: {off} parameters beyond 10 %, factor {worst}");
    // In the linear range the core's step is the Gauss–Newton step: from the minimum with a moved by 3 esd, one cycle
    // brings it back within a tenth of its esd.
    let (c, data) = y2o3();
    let set = c.refinement.build(&c.model).unwrap();
    let mut s = Session::new(c.model.clone(), data.clone(), set.clone(), Options { max_cycles: 30, eps: 1e-4 });
    s.run(|_| true);
    let (amin, sa) = (s.model.get(Field::Cell(0, 0)), s.finish().params.iter().find(|q| q.name == "phase1.cell_a").unwrap().esd);
    let mut moved = s.model.clone();
    for k in 0..3 {
        moved.set(Field::Cell(0, k), amin - 3.0 * sa);
    }
    let mut s1 = Session::new(moved, data, set, Options { max_cycles: 1, eps: 1e-4 });
    s1.step(1, |_| true);
    let back = (s1.model.get(Field::Cell(0, 0)) - amin) / sa;
    println!("a moved 3 esd ({sa:.2e}) from the minimum: one cycle brings it to {back:+.3} esd");
    assert!(back.abs() < 0.1);
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}

#[test]
fn multipliers_by_meaning() {
    let (p05, l05) = shift_run("m05");
    let r = Refinement::from_pcr(&p05, &l05.model).unwrap();
    assert_eq!(r.steering.damping, vec![("phase1.cell_a".to_string(), 0.5)]);
    // a, b and c carry the codeword 240.50; b and c only as the cubic symmetry ties them to a.
    assert!(r.parameters.iter().any(|t| t == &vec![(Field::Cell(0, 0), 1.0), (Field::Cell(0, 1), 1.0), (Field::Cell(0, 2), 1.0)]), "the lone multiplier is not a coefficient: {:?}", r.parameters);
    let (p, l) = shift_run("r05");
    let set = ParameterSet::from_pcr(&p, &l.model).unwrap();
    assert_eq!(set.steering.class_damping, [1.0, 1.0, 0.5, 1.0]);
    for q in &set.free {
        let want = if drs::steering::class_of(q.fields[0].0) == Class::Profile { 0.5 } else { 1.0 };
        assert_eq!(q.damping, want, "{}", q.name);
    }
    // Shared multipliers are the coefficients of a link: Co-YBCO's Cu1 and Co1 occupancies (71.00 and −71.00).
    let d = root().join("runs/Co-YBCO");
    let pc = pcr(&d.join("Cocus-asym.pcr.in"));
    let lc = Model::from_pcr(&pc, &d, "Cocus-asym").unwrap();
    let set = ParameterSet::from_pcr(&pc, &lc.model).unwrap();
    let occ = set.free.iter().find(|q| q.fields.iter().any(|x| matches!(x.0, Field::AtomOcc(0, 2)))).unwrap();
    assert_eq!(occ.fields, vec![(Field::AtomOcc(0, 2), 1.0), (Field::AtomOcc(0, 3), -1.0)]);
    assert_eq!(occ.damping, 1.0);
    // The parameters as the app shows them.
    let desc = set.describe(&lc.model);
    let o = desc.iter().find(|q| q.name == occ.name).unwrap();
    assert_eq!((o.class, o.fields.len(), o.fields[0].max_shift), ("atomic", 2, Some(0.2)));
    assert!(o.fields[0].min == Some(0.0) && o.fields[0].max.is_some(), "{o:?}");
    println!("m05: {:?}; r05: {:?}; Co-YBCO's mixed site: {:?}", r.steering.damping, [1.0, 1.0, 0.5, 1.0], o);
}

#[test]
fn damping_never_changes_an_esd() {
    let (p, mut l) = shift_run("m10");
    let data = l.data.clone().unwrap();
    // The smooth peak window (as `y2o3`): with FullProf's cut χ² is piecewise, and two converged runs can stop on
    // different pieces (0.07 esd apart here).
    l.model.convention.peak_window = false;
    l.model.convention.window_edge = false;
    let plain = ParameterSet::from_pcr(&p, &l.model).unwrap();
    let mut damped = plain.clone();
    damped.steering.class_damping = [0.3, 0.3, 0.3, 0.3];
    for q in damped.free.iter_mut() {
        q.damping = 0.3;
    }
    let opts = Options { max_cycles: 60, eps: 0.01 };
    let mut a = Session::new(l.model.clone(), data.clone(), plain.clone(), opts);
    let mut b = Session::new(l.model.clone(), data.clone(), damped.clone(), opts);
    // At the same values the esd's are the same whatever the damping.
    let (ra, rb) = (a.finish(), b.finish());
    assert!(ra.params.iter().zip(&rb.params).all(|(x, y)| x.esd == y.esd));
    // Both reach the same minimum, with the same esd's; the damped run takes more cycles.
    let (sa, sb) = (a.run(|_| true).clone(), b.run(|_| true).clone());
    let (ra, rb) = (a.finish(), b.finish());
    let mut worst: (f64, f64) = (0.0, 0.0);
    for (x, y) in ra.params.iter().zip(&rb.params) {
        worst.0 = worst.0.max((x.value - y.value).abs() / x.esd);
        worst.1 = worst.1.max((x.esd / y.esd - 1.0).abs());
    }
    println!("damping 1: {} cycles ({}), damping 0.3: {} cycles ({}); values within {:.3} esd, esd's within {:.1e}", ra.cycles, sa.text(), rb.cycles, sb.text(), worst.0, worst.1);
    assert!(sa.converged() && sb.converged() && rb.cycles > ra.cycles);
    assert!(worst.0 < 0.05 && worst.1 < 1e-3, "{worst:?}");
}

#[test]
fn max_shift_shortens_the_step() {
    let (p, l) = shift_run("m10");
    let mut set = ParameterSet::from_pcr(&p, &l.model).unwrap();
    set.steering.max_shift.push((Field::Cell(0, 0), 0.0005));
    let mut s = Session::new(l.model.clone(), l.data.clone().unwrap(), set, Options { max_cycles: 1, eps: 0.01 });
    let mut cyc = None;
    s.step(1, |c| {
        cyc = Some(c.clone());
        true
    });
    let c = cyc.unwrap();
    let a = c.params.iter().find(|q| q.name == "phase1.cell_a").unwrap();
    println!("cell a: computed {:.6}, applied {:.6}; notes {:?}", a.computed, a.shift, c.notes);
    assert!(a.computed > 0.0005 && (a.shift - 0.0005).abs() < 1e-15);
    assert!(c.notes.iter().any(|n| n.contains("phase1.cell_a") && n.contains("maximum shift")));
    // The dictionary's defaults (L108) did not bind: every other parameter moved by its computed shift.
    for q in c.params.iter().filter(|q| q.name != "phase1.cell_a") {
        assert!((q.shift - q.computed).abs() <= 1e-12 * q.computed.abs().max(1e-300), "{}: {} {}", q.name, q.shift, q.computed);
    }
}

#[test]
fn bounds_stop_reflect_periodic() {
    let (_, l) = shift_run("m10");
    let data = l.data.clone().unwrap();
    let one = |lo: f64, hi: f64, b: Behaviour| -> (f64, Vec<String>) {
        let mut set = ParameterSet::from_addresses(&["phase1.cell_a".to_string()], &[], &[("phase1.cell_a".to_string(), Some(lo), Some(hi))], &l.model).unwrap();
        set.behaviour.push((Field::Cell(0, 0), b));
        set.steering.marquardt = false;
        let mut s = Session::new(l.model.clone(), data.clone(), set, Options { max_cycles: 1, eps: 0.01 });
        let mut notes = Vec::new();
        s.step(1, |c| {
            notes = c.notes.clone();
            true
        });
        (s.model.get(Field::Cell(0, 0)), notes)
    };
    let v = 10.6;
    let (free, _) = one(0.0, 100.0, Behaviour::Stop);
    let d = free - v;
    let hi = v + 0.6 * d;
    let (stopped, _) = one(0.0, hi, Behaviour::Stop);
    let (reflected, nr) = one(0.0, hi, Behaviour::Reflect);
    let (wrapped, np) = one(v - 0.1 * d, hi, Behaviour::Periodic);
    println!("a: free {free:.7}, stopped {stopped:.7}, reflected {reflected:.7} {nr:?}, periodic {wrapped:.7} {np:?}");
    assert!(d > 0.0);
    assert!((stopped - hi).abs() < 1e-12);
    assert!((reflected - (v + 0.2 * d)).abs() < 1e-12 && nr.iter().any(|n| n.contains("reflected")));
    assert!((wrapped - (v + 0.3 * d)).abs() < 1e-12 && np.iter().any(|n| n.contains("wrapped")));
    // The physical bounds from the dictionary: 0 ≤ Occ ≤ the site's m/M (Y2O3: Y1 on 8b, Y2 on 24d, O on 48e of I a -3).
    let m = &l.model;
    let lim: Vec<(Option<f64>, Option<f64>)> = (0..3).map(|a| drs::refinable::limits(Field::AtomOcc(0, a), m)).collect();
    println!("occupancy limits: {lim:?}");
    for (a, want) in [(0, 8.0 / 48.0), (1, 24.0 / 48.0), (2, 1.0)] {
        assert!(lim[a].0 == Some(0.0) && (lim[a].1.unwrap() - want).abs() < 1e-4, "{a}: {:?}", lim[a]);
    }
    assert_eq!(drs::refinable::limits(Field::AtomB(0, 0), m), (Some(0.0), None));
}

/// The fields of a parameter set with their coefficients, by parameter, for comparison.
fn form(set: &ParameterSet) -> Vec<Vec<(Field, f64)>> {
    set.free.iter().map(|q| {
        let mut f = q.fields.clone();
        f.sort_by_key(|a| a.0);
        f
    }).collect()
}

#[test]
fn links_with_offset() {
    let d = root().join("runs/Co-YBCO");
    let p = pcr(&d.join("Cocus-asym_conv.pcr"));
    let c = Control::from_pcr(&p, &d, "Cocus-asym", &d).unwrap();
    let data = Model::from_pcr(&p, &d, "Cocus-asym").unwrap().data.unwrap();
    let m = &c.model;
    // FullProf's way: Co1's occupancy on Cu1's codeword with −1, both B's on one codeword.
    let shared = c.refinement.build(m).unwrap();
    // As links: Co1 occupancy = −1 × Cu1 occupancy + their sum; Co1 B = Cu1 B (and Co2/Cu2 likewise).
    let mut linked = c.refinement.clone();
    let follows = [(3, 2), (5, 4)];
    for terms in linked.parameters.iter_mut() {
        terms.retain(|t| !follows.iter().any(|&(co, _)| t.0 == Field::AtomOcc(0, co) || t.0 == Field::AtomB(0, co)));
    }
    linked.parameters.retain(|t| !t.is_empty());
    for (co, cu) in follows {
        let s = m.get(Field::AtomOcc(0, co)) + m.get(Field::AtomOcc(0, cu));
        linked.links.push(Link { value: Field::AtomOcc(0, co), follows: Field::AtomOcc(0, cu), coefficient: -1.0, offset: s });
        linked.links.push(Link { value: Field::AtomB(0, co), follows: Field::AtomB(0, cu), coefficient: 1.0, offset: 0.0 });
    }
    let lset = linked.build(m).unwrap();
    assert_eq!(form(&lset), form(&shared), "the links give FullProf's shared codewords");
    // The same refinement, the link holding p + (s − p) = s.
    let opts = Options { max_cycles: 3, eps: 1e-6 };
    let mut a = Session::new(m.clone(), data.clone(), shared, opts);
    let mut b = Session::new(m.clone(), data.clone(), lset, opts);
    a.run(|_| true);
    b.run(|_| true);
    let worst = m.all_fields().iter().map(|&f| (a.model.get(f) - b.model.get(f)).abs() / a.model.get(f).abs().max(1e-6)).fold(0.0, f64::max);
    let sum = b.model.get(Field::AtomOcc(0, 2)) + b.model.get(Field::AtomOcc(0, 3));
    let s0 = m.get(Field::AtomOcc(0, 2)) + m.get(Field::AtomOcc(0, 3));
    println!("Co-YBCO, 3 cycles: shared codewords and links agree within {worst:.1e}; Cu1 + Co1 occupancy {sum:.6} (start {s0:.6})");
    assert!(worst < 1e-12 && (sum - s0).abs() < 1e-12);
    // A link sets its value at the start: Co1's B made 0.5 above Cu1's.
    let mut off = linked.clone();
    off.links[1].offset = 0.5;
    let s = Session::new(m.clone(), data.clone(), off.build(m).unwrap(), opts);
    assert!((s.model.get(Field::AtomB(0, 3)) - s.model.get(Field::AtomB(0, 2)) - 0.5).abs() < 1e-12 && s.notes.iter().any(|n| n.contains("by its link")));
    // A linked value cannot also be refined, nor follow itself.
    let mut bad = c.refinement.clone();
    bad.links.push(Link { value: Field::AtomOcc(0, 3), follows: Field::AtomOcc(0, 2), coefficient: -1.0, offset: 0.0 });
    assert!(bad.build(m).unwrap_err().contains("also refined"));
    // One value across phases: Egypt's phase 2 W follows phase 1's.
    let d = root().join("runs/Egypt");
    let p = pcr(&d.join("Egypt_qa_Rietveld_conv.pcr"));
    let c = Control::from_pcr(&p, &d, "Egypt_qa_Rietveld", &d).unwrap();
    let mut r = c.refinement.clone();
    r.parameters.retain(|t| !t.iter().any(|x| x.0 == Field::Bov(1)));
    r.links.push(Link { value: Field::Bov(1), follows: Field::Bov(0), coefficient: 1.0, offset: 0.0 });
    let set = r.build(&c.model).unwrap();
    let b = set.free.iter().find(|q| q.fields.iter().any(|x| x.0 == Field::Bov(0))).unwrap();
    assert!(b.fields.contains(&(Field::Bov(1), 1.0)), "{:?}", b.fields);
}

#[test]
fn variables_are_mode_amplitudes() {
    let d = root().join("runs/modes/A_modes");
    let p = pcr(&d.join("a_modes.pcr"));
    let l = Model::from_pcr(&p, &d, "a_modes").unwrap();
    let data = l.data.clone().unwrap();
    let mut mm = l.model.clone();
    let rm = Refinement::from_pcr(&p, &mm).unwrap();
    // From a start away from the minimum: every amplitude at 0.9 of FullProf's.
    let ms = mm.phases[0].modes.clone().unwrap();
    for k in 0..ms.modes.len() {
        mm.set(Field::Mode(0, k), 0.9 * ms.modes[k].amplitude);
    }
    // The same structure without its modes: the coordinates where the amplitudes put them, and one variable per mode
    // moving them as the mode does (∂x/∂Q = c·V), valued at the amplitude.
    let mut mc = mm.clone();
    let ms = mm.phases[0].modes.clone().unwrap();
    mc.phases[0].modes = None;
    let mut rc = rm.clone();
    let refined: Vec<usize> = rm.parameters.iter().filter_map(|t| match t[0].0 {
        Field::Mode(_, k) => Some(k),
        _ => None,
    }).collect();
    rc.parameters.retain(|t| !t.iter().any(|x| matches!(x.0, Field::Mode(..))));
    for (k, mode) in ms.modes.iter().enumerate().filter(|(k, _)| refined.contains(k)) {
        let terms = ms.derivative(k).into_iter().map(|(a, j, c)| (Field::AtomX(0, a, j), c)).collect();
        rc.variables.push(Variable { name: mode.name.clone(), value: mode.amplitude, terms });
    }
    let (sm, sc) = (rm.build(&mm).unwrap(), rc.build(&mc).unwrap());
    assert_eq!(sm.free.len(), sc.free.len());
    let opts = Options { max_cycles: 4, eps: 1e-9 };
    let (mut a, mut b) = (Session::new(mm, data.clone(), sm, opts), Session::new(mc, data, sc, opts));
    let (mut ca, mut cb) = (Vec::new(), Vec::new());
    a.run(|c| {
        ca.push(c.stats.chi2);
        true
    });
    b.run(|c| {
        cb.push(c.stats.chi2);
        true
    });
    let (ra, rb) = (a.finish(), b.finish());
    let mut worst: (f64, f64) = (0.0, 0.0);
    for (k, mode) in ms.modes.iter().enumerate().filter(|(k, _)| refined.contains(k)) {
        let qa = ra.params.iter().find(|q| q.name.ends_with(&format!("{}.mode_amplitude", mode.name))).unwrap();
        let qb = rb.params.iter().find(|q| q.name == mode.name).unwrap();
        worst.0 = worst.0.max((qa.value - qb.value).abs() / qa.esd);
        worst.1 = worst.1.max((qa.esd / qb.esd - 1.0).abs());
        if k == 0 {
            println!("{}: amplitude {:.6} ({:.6}) as a mode, {:.6} ({:.6}) as a variable", mode.name, qa.value, qa.esd, qb.value, qb.esd);
        }
    }
    println!("case A, {} amplitudes of {}, 4 cycles: χ² {ca:?} as modes, {cb:?} as variables; amplitudes within {:.1e} esd, esd's within {:.1e}", refined.len(), ms.modes.len(), worst.0, worst.1);
    assert!(ca.iter().zip(&cb).all(|(x, y)| (x / y - 1.0).abs() < 1e-9));
    assert!(worst.0 < 1e-6 && worst.1 < 1e-6, "{worst:?}");
}

#[test]
fn divergence_guards() {
    let (c, data) = y2o3();
    let m = &c.model;
    let scale = |addrs: &[&str]| ParameterSet::from_addresses(&addrs.iter().map(|s| s.to_string()).collect::<Vec<_>>(), &[], &[], m).unwrap();
    let opts = Options { max_cycles: 20, eps: 1e-9 };
    // χ² rising: the scale over-relaxed (damping 2.5, every step applied) overshoots by 1.5 times its error each cycle.
    let mut set = scale(&["phase1.scale_factor"]);
    set.steering.marquardt = false;
    set.steering.guards.jump = 0.0;
    set.free[0].damping = 2.5;
    let mut start = m.clone();
    start.set(Field::Scale(0), 1.02 * m.get(Field::Scale(0)));
    let mut s = Session::new(start.clone(), data.clone(), set.clone(), opts);
    let stop = s.run(|_| true).clone();
    println!("rising: {} after {} cycles", stop.text(), s.cycle);
    assert!(matches!(stop, Stop::Diverging(Divergence::Rising { cycles: 3, .. })) && stop.code() == "diverging-rising" && s.cycle == 3);
    // A jump of many esd's (from cycle 1 here): the scale 50 % off.
    let mut set = scale(&["phase1.scale_factor"]);
    set.steering.guards.jump_from = 1;
    let mut start = m.clone();
    start.set(Field::Scale(0), 1.5 * m.get(Field::Scale(0)));
    let mut s = Session::new(start, data.clone(), set, opts);
    let stop = s.run(|_| true).clone();
    println!("jump: {}", stop.text());
    assert!(matches!(&stop, Stop::Diverging(Divergence::Jump { parameter, cycle: 1, ratio }) if parameter == "phase1.scale_factor" && *ratio > 100.0));
    // An impossible value set during a pause: W so negative that H² < 0; the session restarts its damping.
    let mut set = scale(&["phase1.caglioti_w"]);
    set.steering.lambda = 1e-3;
    let mut s = Session::new(m.clone(), data.clone(), set, opts);
    assert!(s.step(1, |_| true).is_none());
    s.set_value("phase1.caglioti_w", -1.0).unwrap();
    assert_eq!(s.lambda, 1e-3);
    assert!(s.notes.iter().any(|n| n.contains("λ restarted")));
    let stop = s.step(1, |_| true).unwrap().clone();
    println!("impossible: {}", stop.text());
    assert!(matches!(&stop, Stop::Diverging(Divergence::Impossible { what, value }) if what.contains("Gaussian H²") && *value < 0.0) && s.cycle == 1);
    // A singular matrix, named as R24 names it: two parameters that move the same two background points.
    let b0 = field_of("pattern1.background_points.10", m).unwrap();
    let b1 = field_of("pattern1.background_points.11", m).unwrap();
    let r = Refinement { parameters: vec![vec![(b0, 1.0), (b1, 1.0)], vec![(b1, 1.0), (b0, 1.0)], vec![(Field::Scale(0), 1.0)]], ..Default::default() };
    let set = r.build(m).unwrap();
    let mut s = Session::new(m.clone(), data.clone(), set.clone(), opts);
    let stop = s.run(|_| true).clone();
    println!("singular: {}", stop.text());
    match &stop {
        Stop::Singular { combination, ratio, .. } => {
            assert!(combination.len() == 2 && (combination[0].1.abs() - 0.5f64.sqrt()).abs() < 1e-6 && combination[0].1 * combination[1].1 < 0.0 && *ratio < 1e-12, "{combination:?}");
            assert!(stop.text().contains(&format!("{:+.2}·{} {:+.2}·{}", combination[0].1, combination[0].0, combination[1].1, combination[1].0)));
        }
        other => panic!("{other:?}"),
    }
    assert_eq!(s.cycle, 0);
    // With the guard off the direction is left out and the run goes on.
    let mut set = set;
    set.steering.guards.singular = false;
    let mut s = Session::new(m.clone(), data, set, opts);
    let stop = s.run(|_| true).clone();
    assert!(stop.converged() && s.history[0].singular == 1, "{stop:?}");
}

#[test]
fn session_pause_resume() {
    let (p, l) = shift_run("m10");
    let data = l.data.clone().unwrap();
    let set = ParameterSet::from_pcr(&p, &l.model).unwrap();
    let opts = Options { max_cycles: 4, eps: 1e-12 };
    let mut straight = Session::new(l.model.clone(), data.clone(), set.clone(), opts);
    let stop = straight.run(|_| true).clone();
    assert_eq!(stop, Stop::Cycles { n: 4, le_bail: false });
    // Step 2, saved as text, resumed on a fresh model, run on.
    let mut first = Session::new(l.model.clone(), data.clone(), set.clone(), opts);
    assert!(first.step(2, |_| true).is_none());
    let text = first.to_text();
    let mut resumed = Session::resume(&text, l.model.clone(), data.clone(), set.clone(), opts).unwrap();
    assert_eq!((resumed.cycle, resumed.lambda, resumed.history.len()), (2, first.lambda, 2));
    resumed.run(|_| true);
    let fields = set.fields();
    let worst = fields.iter().map(|&f| (straight.model.get(f) - resumed.model.get(f)).abs() / straight.model.get(f).abs().max(1e-12)).fold(0.0, f64::max);
    let chi = |s: &Session| s.history.iter().map(|c| c.stats.chi2).collect::<Vec<_>>();
    println!("4 cycles straight and 2 + 2 resumed from text: values within {worst:.1e}; χ² {:?} and {:?}", chi(&straight), chi(&resumed));
    assert!(worst < 1e-12);
    assert!(chi(&straight).iter().zip(&chi(&resumed)).all(|(a, b)| (a / b - 1.0).abs() < 1e-12));
    assert_eq!(resumed.history.len(), 4);
    let (h, g) = (&straight.history[3], &resumed.history[3]);
    assert!(h.params.iter().zip(&g.params).all(|(x, y)| x.name == y.name && (x.value - y.value).abs() <= 1e-12 * x.value.abs().max(1e-12) && (x.esd / y.esd - 1.0).abs() < 1e-9));
    // Each parameter's path over the cycles (the app's request 25).
    let paths = resumed.paths();
    assert!(paths.len() == set.free.len() && paths.iter().all(|(_, p)| p.len() == 4 && p.iter().map(|x| x.0).eq(1..=4)));
    // Another parameter set is refused.
    let mut other = set.clone();
    other.free.pop();
    assert!(Session::resume(&text, l.model.clone(), data, other, opts).unwrap_err().contains("free parameters"));
}

#[test]
fn plan_of_steps() {
    let (mut c, data) = y2o3();
    // From a start away from the minimum: the cell 0.05 % off, U doubled.
    for (f, k) in [(Field::Cell(0, 0), 1.0005), (Field::U(0), 2.0)] {
        let v = c.model.get(f);
        c.model.set(f, k * v);
    }
    c.refinement.plan = vec![
        PlanStep { name: "scale and global".into(), free: vec!["phase1.scale_factor".into(), "global".into()], damping: None, cycles: Some(10), eps: Some(0.01) },
        PlanStep { name: "profile".into(), free: vec!["profile".into()], damping: Some(0.8), cycles: Some(20), eps: Some(0.01) },
        PlanStep { name: "all".into(), free: vec!["all".into()], damping: None, cycles: Some(30), eps: Some(0.01) },
    ];
    // Written and read back.
    let back = Control::from_text(&c.to_text()).unwrap();
    assert_eq!(back.refinement, c.refinement);
    let mut m = c.model.clone();
    for k in 0..3 {
        let rk = c.refinement.for_step(k, &m).unwrap();
        let set = rk.build(&m).unwrap();
        for q in &set.free {
            assert!(c.refinement.plan[k].frees(&q.name, q.fields[0].0), "step {k}: {}", q.name);
        }
        if k == 1 {
            assert!(set.free.iter().all(|q| drs::steering::class_of(q.fields[0].0) == Class::Profile && (q.damping - 0.8).abs() < 1e-15));
        }
        let mut s = Session::new(m.clone(), data.clone(), set, Options { max_cycles: rk.cycles.unwrap(), eps: rk.eps.unwrap() });
        let stop = s.run(|_| true).clone();
        println!("step {} ({}): {} parameters, {} cycles, {}; χ² {:.5}", k + 1, c.refinement.plan[k].name, s.set.free.len(), s.cycle, stop.text(), s.history.last().unwrap().stats.chi2);
        assert!(stop.converged());
        m = s.model.clone();
    }
    // The single run from the same start reaches the same minimum.
    let set = c.refinement.build(&c.model).unwrap();
    let mut s = Session::new(c.model.clone(), data.clone(), set.clone(), Options { max_cycles: 60, eps: 0.01 });
    s.run(|_| true);
    let r = s.finish();
    let worst = set.free.iter().zip(&r.params).map(|(q, x)| (m.get(q.fields[0].0) - s.model.get(q.fields[0].0)).abs() / x.esd.max(1e-300)).fold(0.0, f64::max);
    println!("the single run: {} cycles, χ² {:.5}; the plan's values within {worst:.3} esd of it", r.cycles, r.stats.chi2);
    assert!(worst < 0.1);
    // A step that frees nothing is an error with its name.
    let mut bad = c.refinement.clone();
    bad.plan.push(PlanStep { name: "nothing".into(), free: vec!["phase7.cell_a".into()], damping: None, cycles: None, eps: None });
    assert!(bad.for_step(3, &m).unwrap_err().contains("'nothing'"));
}

#[test]
fn restraint_weights_and_contributions() {
    let path = root().join("FP_tests/lamn_3t2.pcr");
    let p = pcr(&path);
    let l = Model::from_pcr(&p, path.parent().unwrap(), "lamn_3t2").unwrap();
    let mut m = l.model.clone();
    let set = ParameterSet::from_pcr(&p, &m).unwrap();
    let r = drs::lsq::refine(&mut m, l.data.as_ref().unwrap(), &set, Options { max_cycles: 3, eps: p.eps }, |_| {});
    let dof = (r.stats.n - r.stats.p) as f64;
    assert!(!r.restraints.is_empty());
    let mut total = 0.0;
    for q in &r.restraints {
        let pull = ((q.target - q.value) / q.sigma).powi(2);
        // weight = χ²_ν/σ² (σ as pulled: FullProf's angle restraints with π/180 of their weight).
        let k = q.weight * q.sigma * q.sigma / r.stats.chi2;
        assert!((k - 1.0).abs() < 1e-12 || (k - std::f64::consts::PI / 180.0).abs() < 1e-12, "{}: {k}", q.label);
        assert!((q.contribution * dof - k * r.stats.chi2 * pull).abs() <= 1e-9 * (q.contribution * dof).max(1e-300), "{}", q.label);
        total += q.contribution;
    }
    assert!((drs::lsq::restraints_part(&r.restraints) - total).abs() <= 1e-15 * total.max(1e-300));
    println!("lamn_3t2: {} restraints add {total:.4} to χ² = {:.4}; the largest {:?}", r.restraints.len(), r.stats.chi2, r.restraints.iter().max_by(|a, b| a.contribution.total_cmp(&b.contribution)).map(|q| (&q.label, q.weight, q.contribution)));
}

#[test]
fn control_file_guardrails() {
    let (c, _) = y2o3();
    let text = c.to_text();
    let with = |extra_refine: &str, tables: &str| -> Vec<drs::control::syntax::Diagnostic> {
        let t = text.replacen("[refine]\n", &format!("[refine]\n{extra_refine}\n"), 1).replacen("[convention]", &format!("{tables}\n[convention]"), 1);
        drs::control::schema::read_text(&t).1
    };
    let has = |d: &[drs::control::syntax::Diagnostic], warning: bool, part: &str| d.iter().any(|x| x.warning == warning && x.message.contains(part));
    let d = with("damping = { atomic = 0.0, profile = 1.5 }", "");
    assert!(has(&d, false, "damping is above 0") && has(&d, true, "over-relaxation"), "{d:?}");
    let d = with("damping = { atomik = 1.0 }", "");
    assert!(d.iter().any(|x| x.code == "unknown-name" && x.replacement.as_deref() == Some("atomic")), "{d:?}");
    let d = with("limits = [ { value = \"pattern1.zero_shift\", behaviour = \"periodic\" } ]", "");
    assert!(has(&d, false, "periodic but has not two bounds"), "{d:?}");
    let d = with("lambda = -1\ndivergence = { rising = 2, jumps = 5 }", "");
    assert!(has(&d, false, "λ's start") && d.iter().any(|x| x.replacement.as_deref() == Some("jump")), "{d:?}");
    let d = with("", "[[refine.step]]\nname = \"empty\"\nfree = [\"phase3.cell_a\"]\n");
    assert!(has(&d, false, "frees no parameter"), "{d:?}");
    let d = with("steering = [ { value = \"phase1.cell_a\", max_shift = 0.0 } ]", "[[refine.link]]\nvalue = \"phase1.Y.atom_biso\"\nfollows = \"phase1.Y.atom_biso\"\n");
    assert!(has(&d, false, "maximum shift is above 0") && has(&d, false, "linked to itself"), "{d:?}");
    // A file with every new table reads back to the same refinement.
    let mut c2 = c.clone();
    let m = &c2.model;
    let by = field_of("phase1.Y:2.atom_biso", m).unwrap();
    let bo = field_of("phase1.O.atom_biso", m).unwrap();
    c2.refinement.parameters.retain(|t| t[0].0 != by);
    c2.refinement.links.push(Link { value: by, follows: field_of("phase1.Y.atom_biso", m).unwrap(), coefficient: 1.0, offset: 0.0 });
    c2.refinement.parameters.retain(|t| t[0].0 != bo);
    c2.refinement.variables.push(Variable { name: "B_O".into(), value: 0.0, terms: vec![(bo, 1.0)] });
    c2.refinement.steering.lambda = 1e-3;
    c2.refinement.steering.class_damping = [0.8, 1.0, 0.9, 1.0];
    c2.refinement.steering.damping.push(("phase1.cell_a".into(), 0.5));
    c2.refinement.steering.max_shift.push((Field::Cell(0, 0), 0.01));
    c2.refinement.steering.guards.rising = 5;
    c2.refinement.limits.push((Field::Zero, Some(-0.1), Some(0.1)));
    c2.refinement.behaviour.push((Field::Zero, Behaviour::Reflect));
    let t2 = c2.to_text();
    let back = Control::from_text(&t2).unwrap_or_else(|e| panic!("{e:?}\n{t2}"));
    assert_eq!(back.refinement, c2.refinement);
    assert_eq!(back.to_text(), t2);
    let set = back.refinement.build(&back.model).unwrap();
    assert!(set.free.iter().any(|q| q.name == "B_O" && q.value == Some(0.0)));
    println!("{}", t2.lines().skip_while(|l| !l.starts_with("[refine]")).take_while(|l| !l.starts_with("[convention]")).filter(|l| !l.starts_with("parameters") && !l.starts_with("  [")).collect::<Vec<_>>().join("\n"));
}

#[test]
fn command_line() {
    let (mut c, _) = y2o3();
    let tmp = std::env::temp_dir().join(format!("drs_steering_{}", std::process::id()));
    std::fs::create_dir_all(&tmp).unwrap();
    std::fs::copy(root().join("runs/Y2O3/y2o3.dat"), tmp.join("y2o3.dat")).unwrap();
    c.data.as_mut().unwrap().file = "y2o3.dat".into();
    c.refinement.plan = vec![
        PlanStep { name: "scale".into(), free: vec!["phase1.scale_factor".into()], damping: None, cycles: Some(3), eps: Some(0.1) },
        PlanStep { name: "all".into(), free: vec!["all".into()], damping: None, cycles: Some(5), eps: Some(0.1) },
    ];
    let drc = tmp.join("y2o3.drc");
    std::fs::write(&drc, c.to_text()).unwrap();
    let run = |args: &[&str]| -> String {
        let o = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).arg("refine").arg(&drc).args(args).output().unwrap();
        let s = format!("{}{}", String::from_utf8_lossy(&o.stdout), String::from_utf8_lossy(&o.stderr));
        assert!(o.status.success(), "{s}");
        s
    };
    let plan = run(&["--plan"]);
    assert!(plan.contains("step 1 of 2: scale (1 parameters") && plan.contains("step 2 of 2: all"), "{plan}");
    let one = run(&["--step", "2", "--parameters"]);
    assert!(one.contains("step 2 of 2") && !one.contains("step 1 of 2") && one.contains("shift"), "{one}");
    let sess = tmp.join("y2o3.session");
    let a = run(&["--session", sess.to_str().unwrap(), "--cycles", "1"]);
    let b = run(&["--session", sess.to_str().unwrap(), "--cycles", "1"]);
    let saved = std::fs::read_to_string(&sess).unwrap();
    println!("{}", b.lines().take(6).collect::<Vec<_>>().join("\n"));
    assert!(a.contains("Paused after cycle 1") || a.contains("every shift below"), "{a}");
    assert!(b.contains("resumed from") && saved.contains("\ncycle 2\n"), "{b}\n{saved}");
    let _ = std::fs::remove_dir_all(&tmp);
}
