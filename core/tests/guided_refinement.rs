//! Guided refinement (P55, `guide`, L302): the documented order of refinement (FullProf's manual, p. 14; McCusker et
//! al. 1999, §7) as a plan of steps for a model.
//! - **Y₂O₃** (cubic, X-ray, pseudo-Voigt, background points): the steps as the manual and McCusker give them, those
//!   that add nothing left out (its sites are full, its points kept fixed, no betas, no microstructure term).
//! - **Eight FullProf runs**: each plan cumulative, each step moving every value it names (`Refinement::for_step` of
//!   `guided`, built through the symmetry), the steps each run's model calls for: Pd₃MnD's five partly occupied
//!   atoms, Co-YBCO's S/L and H/L (Thompson–Cox–Hastings), NdNi-HTO's betas and, over its resolution file, its U, X, Y
//!   in the microstructure, prba10k's needles, NdNi-LTOa's quartic strain, case A's mode amplitudes.
//! - **The control file**: the guided refinement written and read back, plan and parameters the same.
//! - **The plan run** on Y₂O₃ from the tutorial's starting file: every step without divergence, the last to every
//!   shift below 0.1 esd, at a χ² below that of the file's own codewords refined to convergence (it frees more).
//! - **The plan run as one** (`lsq::run_plan`, P37): the same steps give the same χ² to the last bit; the results file
//!   (`drs results --plan`) records every cycle of every step, numbered on over the steps, each with the number and
//!   name of the step that ran it, and the parameters of each cycle are those its step frees.

use drs::control::Control;
use drs::guide::{documented_plan, guided};
use drs::powder::Model;
use drs::refinable::Refinement;
use std::path::Path;

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap().iter().map(|&b| b as char).collect()
}

fn dir(run: &str) -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs").join(run)
}

fn load(run: &str, code: &str) -> (drs::pcr::Pcr, Model) {
    let mut p = drs::pcr::read(&latin1(&dir(run).join(format!("{code}.pcr")))).unwrap();
    p.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
    let m = Model::from_pcr(&p, &dir(run), code).unwrap_or_else(|e| panic!("{run}: {e}")).model;
    (p, m)
}

#[test]
fn y2o3_documented_order() {
    let (_, m) = load("Y2O3", "y2o3_structure");
    let steps = documented_plan(&m);
    let got: Vec<(String, Vec<String>, Vec<String>)> = steps.iter().map(|s| (s.step.name.clone(), s.adds.clone(), s.drops.clone())).collect();
    let v = |x: &[&str]| x.iter().map(|s| s.to_string()).collect::<Vec<_>>();
    let want = vec![
        ("1 scale factor", v(&["phase1.scale_factor"]), v(&[])),
        ("2 zero point, first background term and cell", v(&["pattern1.zero_shift", "phase1.cell_a"]), v(&[])),
        ("3 positions of the heavier atoms", v(&["phase1.Y:2.atom_position.x"]), v(&[])),
        ("4 all positions and the overall B", v(&["phase1.overall_b", "phase1.O.atom_position.x", "phase1.O.atom_position.y", "phase1.O.atom_position.z"]), v(&[])),
        (
            "5 peak shape and asymmetry",
            v(&["phase1.caglioti_u", "phase1.caglioti_v", "phase1.caglioti_w", "phase1.pv_eta_slope", "phase1.pv_eta0", "phase1.berar_baldinozzi.1", "phase1.berar_baldinozzi.2", "phase1.berar_baldinozzi.3", "phase1.berar_baldinozzi.4"]),
            v(&[]),
        ),
        ("6 each atom's isotropic B for the overall B", v(&["phase1.Y.atom_biso", "phase1.Y:2.atom_biso", "phase1.O.atom_biso"]), v(&["phase1.overall_b"])),
        ("7 Sycos and Sysin", v(&["pattern1.sample_displacement", "pattern1.sample_transparency"]), v(&[])),
        ("8 all together, to convergence", v(&[]), v(&[])),
    ];
    let want: Vec<(String, Vec<String>, Vec<String>)> = want.into_iter().map(|(n, a, d)| (n.to_string(), a, d)).collect();
    assert_eq!(got, want);
    assert_eq!(steps.iter().map(|s| s.source).collect::<Vec<_>>(), ["FullProf manual, p. 14", "FullProf manual, p. 14", "McCusker et al. 1999, §7", "FullProf manual, p. 14", "FullProf manual, p. 14", "FullProf manual, p. 14", "FullProf manual, p. 14", "McCusker et al. 1999, §7"]);
    assert!(steps[..7].iter().all(|s| s.step.cycles == Some(5) && s.step.eps.is_none()));
    assert_eq!((steps[7].step.cycles, steps[7].step.eps), (None, Some(0.1)));
}

#[test]
fn plans_of_eight_runs() {
    let runs = [("Y2O3", "y2o3_structure"), ("NdNi-LTOa", "NdNi_LTOa"), ("NdNi-HTO", "NdNi_HTO"), ("prba10k", "prba10k"), ("Pd3MnD", "pdfinal"), ("modes/A_modes", "a_modes"), ("Co-YBCO", "Cocus-asym"), ("LaB6", "LaB6rx")];
    for (run, code) in runs {
        let (p, m) = load(run, code);
        let steps = documented_plan(&m);
        // Cumulative: each step is the one before, less what it drops, and what it adds.
        for w in steps.windows(2) {
            let mut expect: Vec<String> = w[0].step.free.iter().filter(|a| !w[1].drops.contains(a)).cloned().collect();
            expect.extend(w[1].adds.iter().cloned());
            if !w[1].step.name.contains("all together") {
                assert_eq!(w[1].step.free, expect, "{run}: {}", w[1].step.name);
            }
        }
        // Every value a step names is freed by it.
        let g = guided(&m, &Refinement::from_pcr(&p, &m).unwrap());
        for (k, s) in g.plan.iter().enumerate() {
            let set = g.for_step(k, &m).unwrap().build(&m).unwrap_or_else(|e| panic!("{run} step {}: {e}", k + 1));
            // Moved by a free parameter (one of the file's may move several values: NdNi-LTOa ties O2's B to O1's).
            let moved = set.fields();
            let missing: Vec<&String> = s.free.iter().filter(|a| !moved.contains(&drs::refinable::field_of(a, &m).unwrap())).collect();
            assert!(missing.is_empty(), "{run} step {}: {missing:?} not freed", s.name);
        }
        let step = |w: &str| steps.iter().find(|s| s.step.name.contains(w));
        let adds = |w: &str| step(w).map_or(Vec::new(), |s| s.adds.clone());
        println!("{run:14} {} steps: {}", steps.len(), steps.iter().map(|s| s.step.name.split_once(' ').unwrap().1.split(' ').next().unwrap()).collect::<Vec<_>>().join(", "));
        match run {
            "Pd3MnD" => assert_eq!(adds("occupancies"), ["phase1.MN.atom_occupancy", "phase1.PD.atom_occupancy", "phase1.MN:3.atom_occupancy", "phase1.PD:4.atom_occupancy", "phase1.D.atom_occupancy"]),
            "Co-YBCO" => assert!(adds("peak shape").iter().any(|a| a == "phase1.axial_divergence_s_l") && adds("peak shape").iter().any(|a| a == "phase1.axial_divergence_d_l")),
            "NdNi-HTO" => {
                assert!(!adds("betas").is_empty() && adds("betas").iter().all(|a| a.contains(".atom_betas.")), "{:?}", adds("betas"));
                // Over a resolution file, U, X and Y are the sample's: microstructure, not peak shape.
                assert!(adds("microstructure").iter().any(|a| a == "phase1.caglioti_u") && !adds("peak shape").iter().any(|a| a == "phase1.caglioti_u"));
            }
            "prba10k" => assert!(adds("microstructure").iter().any(|a| a == "phase1.shape_size")),
            "NdNi-LTOa" => assert!(adds("microstructure").iter().any(|a| a.starts_with("phase1.quartic_strain"))),
            "modes/A_modes" => assert!(adds("heavier").iter().all(|a| a.ends_with(".mode_amplitude")) && !adds("heavier").is_empty()),
            "Y2O3" => assert!(step("occupancies").is_none() && step("rest of the background").is_none() && step("betas").is_none()),
            _ => {}
        }
        assert!(steps.last().unwrap().step.name.contains("all together"));
    }
}

#[test]
fn plan_in_the_control_file() {
    let d = dir("Pd3MnD");
    let mut p = drs::pcr::read(&latin1(&d.join("pdfinal.pcr"))).unwrap();
    p.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
    let mut c = Control::from_pcr(&p, &d, "pdfinal", &d).unwrap();
    c.refinement = guided(&c.model, &c.refinement);
    let text = c.to_text();
    assert_eq!(text.matches("[[refine.step]]").count(), c.refinement.plan.len());
    let back = Control::from_text_in(&text, &d).unwrap();
    assert_eq!(back.refinement.plan, c.refinement.plan);
    assert_eq!(back.refinement.parameters.len(), c.refinement.parameters.len());
}

#[test]
fn plan_run_on_y2o3() {
    let d = dir("Y2O3");
    let p = drs::pcr::read(&latin1(&d.join("y2o3_structure.pcr.in"))).unwrap();
    let l = Model::from_pcr_with_data(&p, &d, "y2o3_structure", None).unwrap();
    let (m0, data) = (l.model, l.data.unwrap());
    let base = Refinement::from_pcr(&p, &m0).unwrap();
    let mut plain = m0.clone();
    let r0 = drs::lsq::refine(&mut plain, &data, &base.build(&m0).unwrap(), drs::lsq::Options { max_cycles: 50, eps: 0.1 }, |_| {});
    assert!(r0.converged, "{}", r0.stop);
    let g = guided(&m0, &base);
    let mut m = m0.clone();
    let mut last = None;
    for k in 0..g.plan.len() {
        let rk = g.for_step(k, &m).unwrap();
        let set = rk.build(&m).unwrap();
        let r = drs::lsq::refine(&mut m, &data, &set, drs::lsq::Options { max_cycles: rk.cycles.unwrap_or(50), eps: rk.eps.unwrap_or(0.1) }, |_| {});
        assert!(!matches!(r.reason, drs::lsq::Stop::Diverging(_) | drs::lsq::Stop::Singular { .. }), "step {}: {}", g.plan[k].name, r.stop);
        println!("{:42} {:3} parameters  χ² {:.4}  {}", g.plan[k].name, r.params.len(), r.stats.chi2, r.stop);
        last = Some(r);
    }
    let r = last.unwrap();
    println!("the file's own codewords: {} parameters, χ² {:.4}", r0.params.len(), r0.stats.chi2);
    assert!(r.converged && r.stop.contains("0.1"), "{}", r.stop);
    assert!(r.stats.chi2 < r0.stats.chi2, "{} against {}", r.stats.chi2, r0.stats.chi2);
}

#[test]
fn plan_run_into_the_results_file() {
    let d = dir("Y2O3");
    let p = drs::pcr::read(&latin1(&d.join("y2o3_structure.pcr.in"))).unwrap();
    let data = Model::from_pcr_with_data(&p, &d, "y2o3_structure", None).unwrap().data.unwrap();
    let out = Path::new(env!("CARGO_MANIFEST_DIR")).join("target/test-plan-results");
    let mut c = Control::from_pcr(&p, &d, "y2o3_structure", &out).unwrap();
    c.refinement = guided(&c.model, &c.refinement);
    c.refinement.cycles = Some(50);
    let g = c.refinement.clone();

    // Step by step, as plan_run_on_y2o3 runs it.
    let mut m = c.model.clone();
    let mut by_hand = Vec::new();
    for k in 0..g.plan.len() {
        let rk = g.for_step(k, &m).unwrap();
        let set = rk.build(&m).unwrap();
        let r = drs::lsq::refine(&mut m, &data, &set, drs::lsq::Options::stated(rk.cycles, rk.eps), |_| {});
        by_hand.push((r.history.len(), r.stats.chi2));
    }

    // As one run.
    let steps: Vec<usize> = (0..g.plan.len()).collect();
    let mut seen = Vec::new();
    let mut model = c.model.clone();
    let (r, session) = drs::lsq::run_plan(&mut model, &data, &g, &steps, |s| seen.push((s.index, s.stop.is_some())), |_| true).unwrap();
    assert_eq!(seen, steps.iter().flat_map(|&k| [(k, false), (k, true)]).collect::<Vec<_>>(), "each step announced as it starts and ends");
    assert_eq!(r.stats.chi2, by_hand.last().unwrap().1, "the same χ² as the steps run one by one");
    assert_eq!(model.get(drs::jacobian::Field::Cell(0, 0)), m.get(drs::jacobian::Field::Cell(0, 0)));
    let total: usize = by_hand.iter().map(|s| s.0).sum();
    assert_eq!((r.cycles, r.history.len(), r.cycle_steps.len()), (total, total, total));
    assert!(r.converged, "{}", r.stop);

    // Into the results file and back.
    let (text, _) = drs::results::write(&Control { model: model.clone(), ..c.clone() }, &data, Some(&session.set), Some(&r), "y2o3 plan").unwrap();
    let back = drs::results::read(&text).unwrap();
    assert_eq!(back.cycles.iter().map(|c| c.cycle).collect::<Vec<_>>(), (1..=total).collect::<Vec<_>>(), "numbered on over the steps");
    assert_eq!(back.cycle_steps, r.cycle_steps);
    let mut i = 0;
    for (k, (n, _)) in by_hand.iter().enumerate() {
        for cy in &back.cycles[i..i + n] {
            assert_eq!(back.cycle_steps[cy.cycle - 1], (k + 1, g.plan[k].name.clone()), "cycle {}", cy.cycle);
            let free = g.for_step(k, &c.model).unwrap().build(&c.model).unwrap();
            let names: Vec<&str> = cy.params.iter().map(|q| q.name.as_str()).collect();
            assert_eq!(names, free.free.iter().map(|f| f.name.as_str()).collect::<Vec<_>>(), "cycle {}: the parameters of step {}", cy.cycle, k + 1);
        }
        println!("step {} ({}): cycles {}–{}", k + 1, g.plan[k].name, i + 1, i + n);
        i += n;
    }
    // A run without a plan writes no step.
    let plain = text.lines().filter(|l| l.contains("_drs_cycle_step")).count();
    assert_eq!(plain, 2);
    let set = c.refinement.build(&c.model).unwrap();
    let mut m1 = c.model.clone();
    let r1 = drs::lsq::refine(&mut m1, &data, &set, drs::lsq::Options { max_cycles: 2, eps: 0.1 }, |_| {});
    let (t1, _) = drs::results::write(&c, &data, Some(&set), Some(&r1), "y2o3").unwrap();
    assert!(!t1.contains("_drs_cycle_step") && drs::results::read(&t1).unwrap().cycle_steps.is_empty());
}
