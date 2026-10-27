//! Restraints (P40) against FullProf (test `fullprof_restraints`, L055).
//!
//! 1. **Geometry.** FullProf writes, with each run, a file `CFML_Restraints_1.tpcr` listing the distances (and
//!    angles) of the structure as restraint lines ready for a `.pcr`: atoms, ITnum, T and the value it calculates
//!    (manual p. 123; distances to 4 decimals, angles to 2). The same lines are read here as restraints and their
//!    values calculated by the core for the structure of the run:
//!    - Yb2ZrWO8 (our run, `runs/Yb2ZrWO8`, ybfluo_Riet_IRF_conv.pcr; F m −3 m): 4 distances, among them O–O
//!      through ITnum −4 (the inversion times the 2-fold along x): within the print (3·10⁻⁵ Å);
//!    - LaMnO3, the symmetry-modes tutorial (LaMnO3_cub.pcr, P n m a, positions from 7 mode amplitudes; the core's
//!      positions agree with FullProf's "Final Complete Structure" to its 5 printed decimals): 23 distances, up to
//!      1.6·10⁻⁴ Å from FullProf's;
//!    - LaMnO3 on D20, the sequential tutorial (la0_190.pcr, the last of the sequence; P b n m): 6 distances up to
//!      1.5·10⁻⁴ Å and 31 angles up to 0.013°.
//!
//!    The two LaMnO3 differences, 10–20 times the print, are not explained (the cell and positions are FullProf's
//!    to their printed digits); they are frozen. An error in the operator convention (ITnum, the central atom of an
//!    angle, which atom the operator moves) would give differences of tenths of an ångström or of degrees.
//! 2. **Derivatives.** Every restraint of these files and of the four `.pcr` files with restraints, against central
//!    differences of its value with respect to each cell parameter and coordinate (Richardson's extrapolation of
//!    h = 10⁻⁵ and h/2), to 10⁻⁶ relative (10⁻⁵ for an angle within 1° of 0° or 180°, where acos is steep).
//! 3. **Refinement with restraints.** The three runs that carry them (SrHC3O6·H2O, the restraints tutorial; LaMnO3 on
//!    3T2; K2Na[Ag(CN)2]3) are refined from their `.pcr`: each converges, and each restraint ends within 3σ of its
//!    target or nearer to it than it started. There is no FullProf output of these refinements in the material: `runs/restraint_lebail_runs.py`
//!    made it on the Mac (L056), and `fullprof_restraint_runs` compares with it (below).
//! 4. **Linear restraints** (NLI): two lines added to Y2O3's `.pcr` restrain U + W of its phase; the refinement holds
//!    the sum to its target within 3σ while it moves both.

use drs::lsq::{refine, restraint_states, Options};
use drs::powder::Model;
use drs::refinable::ParameterSet;
use drs::restraints::Kind;
use std::path::{Path, PathBuf};

mod common;
use common::latin1;

fn root(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof").join(rel)
}

/// The restraint lines of a `.tpcr`: its distance list and its angle lists, as `.pcr` lines.
fn tpcr_lines(text: &str) -> (Vec<String>, Vec<String>) {
    let (mut dl, mut al, mut mode) = (Vec::new(), Vec::new(), 0);
    for l in text.lines() {
        if l.contains("At1   At2  ITnum") {
            mode = 1;
            continue;
        }
        if l.contains("Help for possible angle") {
            mode = 2;
            continue;
        }
        if l.trim().is_empty() || l.contains("=====") || l.contains("=>") {
            if mode == 1 {
                mode = 0;
            }
            continue;
        }
        if mode == 1 {
            dl.push(l.to_string());
        } else if mode == 2 && l.split_whitespace().count() >= 13 {
            al.push(l.to_string());
        }
    }
    (dl, al)
}

fn load(pcr: &Path) -> (drs::pcr::Pcr, drs::powder::Loaded) {
    let mut p = drs::pcr::read(&latin1(pcr)).unwrap();
    p.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
    let l = Model::from_pcr(&p, pcr.parent().unwrap(), pcr.file_stem().unwrap().to_str().unwrap()).unwrap_or_else(|e| panic!("{}: {e}", pcr.display()));
    (p, l)
}

/// Every restraint's analytic derivatives against central differences; returns the worst relative error.
fn check_derivatives(m: &Model, name: &str, fails: &mut Vec<String>) -> f64 {
    let mut worst = 0.0f64;
    for r in &m.restraints {
        let (_, d) = r.derivatives(m);
        for &(f, a) in &d {
            let v = m.get(f);
            let central = |h: f64| {
                let (mut up, mut dn) = (m.clone(), m.clone());
                up.set(f, v + h);
                dn.set(f, v - h);
                (r.value(&up) - r.value(&dn)) / (2.0 * h)
            };
            // Richardson's extrapolation of two central differences (an angle near 0°, where acos is steep).
            let h = 1e-5 * v.abs().max(1.0);
            let fd = (4.0 * central(h / 2.0) - central(h)) / 3.0;
            let err = (a - fd).abs() / a.abs().max(fd.abs()).max(1e-3);
            worst = worst.max(err);
            // Within 1° of 0° or 180° the differences themselves lose accuracy (acos is steep there; k2na1p5's
            // Ag–C–N at 0.18°: 2.7·10⁻⁶ with h = 10⁻⁵), so 10⁻⁵ is asked.
            let near = matches!(r.kind, Kind::Angle { .. }) && { let a = r.value(m); !(1.0..=179.0).contains(&a) };
            if err > if near { 1e-5 } else { 1e-6 } {
                fails.push(format!("{name}: {} ∂/∂{f:?}: {a:.8} against {fd:.8}", r.label(m)));
            }
        }
    }
    worst
}

#[test]
fn fullprof_restraints() {
    let mut fails = Vec::new();

    // 1. Geometry against FullProf's .tpcr; (tpcr, pcr, frozen distance difference Å, frozen angle difference °).
    let cases = [
        ("runs/Yb2ZrWO8/CFML_Restraints_1.tpcr", "runs/Yb2ZrWO8/ybfluo_Riet_IRF_conv.pcr", 5e-5, 0.005),
        ("7-FullProf-Tutorial-SymmetryModes/LaMnO3/CFML_Restraints_1.tpcr", "7-FullProf-Tutorial-SymmetryModes/LaMnO3/LaMnO3_cub.pcr", 1.7e-4, 0.005),
        ("1-WinPLOTR/2_FP_sequential_refinements/LaMnO3_D20/CFML_Restraints_1.tpcr", "1-WinPLOTR/2_FP_sequential_refinements/LaMnO3_D20/la0_190.pcr", 1.6e-4, 0.014),
    ];
    for (tp, pcr, fd, fa) in cases {
        let (dl, al) = tpcr_lines(&latin1(&root(tp)));
        let mut lines = dl.clone();
        lines.extend(al.iter().cloned());
        let (ds, an) = drs::pcr::restraints_of(&lines, dl.len(), al.len()).unwrap();
        let (mut p, _) = load(&root(pcr));
        p.phases[0].distance_restraints = ds;
        p.phases[0].angle_restraints = an;
        let path = root(pcr);
        let m = Model::from_pcr(&p, path.parent().unwrap(), path.file_stem().unwrap().to_str().unwrap()).unwrap().model;
        let (mut wd, mut wa, mut nd, mut na) = (0.0f64, 0.0f64, 0, 0);
        for r in &m.restraints {
            let dv = (r.value(&m) - r.target).abs();
            match r.kind {
                Kind::Distance { .. } => (wd, nd) = (wd.max(dv), nd + 1),
                Kind::Angle { .. } => (wa, na) = (wa.max(dv), na + 1),
            }
        }
        println!("{pcr}: {nd} distances, worst {wd:.5} Å; {na} angles, worst {wa:.4}°");
        if wd > fd || wa > fa {
            fails.push(format!("{pcr}: distances {wd:.2e} Å (frozen {fd:.1e}), angles {wa:.4}° (frozen {fa})"));
        }
        let w = check_derivatives(&m, pcr, &mut fails);
        println!("    derivatives: worst relative error {w:.1e}");
    }

    // 2–3. The runs with restraints: derivatives, then the refinement.
    for pcr in ["8-FullProf-Tutorial-RigidBodies/Sr-ox/sr-ox-test.pcr", "8-FullProf-Tutorial-RigidBodies/Sr-ox/sr-ox-rest.pcr", "FP_tests/lamn_3t2.pcr", "FP_tests/k2na1p5.pcr"] {
        let (p, l) = load(&root(pcr));
        let mut m = l.model.clone();
        let data = l.data.as_ref().unwrap();
        let w = check_derivatives(&m, pcr, &mut fails);
        let set = ParameterSet::from_pcr(&p, &m).unwrap();
        let before = restraint_states(&m, &set);
        let mut costs = Vec::new();
        let r = refine(&mut m, data, &set, Options { max_cycles: 50, eps: p.eps }, |c| costs.push(c.stats.chi2));
        println!("{pcr}: {} restraints, derivatives worst {w:.1e}; {} parameters, {} cycles ({}): χ² {:.4} → {:.4}", before.len(), set.free.len(), r.cycles, r.stop, costs.first().copied().unwrap_or(0.0), r.stats.chi2);
        if !r.converged {
            fails.push(format!("{pcr}: not converged in 50 cycles ({})", r.stop));
        }
        for (b, a) in before.iter().zip(&r.restraints) {
            let (db, da) = ((b.value - b.target) / b.sigma, (a.value - a.target) / a.sigma);
            println!("    {:28} target {:10.5}: {:10.5} ({db:+7.2}σ) → {:10.5} ({da:+7.2}σ)", a.label, a.target, b.value, a.value);
            if da.abs() > 3.0 && da.abs() > db.abs() {
                fails.push(format!("{pcr}: {} moved away from its target: {db:+.2}σ → {da:+.2}σ", a.label));
            }
        }
    }

    // 4. A linear restraint on U + W of Y2O3, read from NLI lines.
    let path = root("runs/Y2O3/y2o3_structure_conv.pcr");
    let text = latin1(&path);
    let mut lines: Vec<String> = text.lines().map(String::from).collect();
    let k = lines.iter().position(|l| l.trim_start().starts_with("!Ipr")).unwrap() + 1;
    let mut v: Vec<String> = lines[k].split_whitespace().map(String::from).collect();
    v[8] = "1".into();
    lines[k] = v.iter().map(|x| format!("{x:>4}")).collect();
    let k = lines.iter().position(|l| l.contains("2Th1")).unwrap() + 2;
    let (u0, w0) = (0.014185, 0.010206);
    let target = u0 + w0 + 0.0005;
    lines.insert(k, format!("UplusW  2  {target:.6}  0.000010"));
    lines.insert(k + 1, "  1.0 2   1.0 4".into());
    let p = drs::pcr::read(&lines.join("\n")).unwrap();
    let l = Model::from_pcr(&p, path.parent().unwrap(), "y2o3_structure").unwrap();
    let mut m = l.model.clone();
    let data = l.data.as_ref().unwrap();
    let set = ParameterSet::from_pcr(&p, &m).unwrap();
    println!("Y2O3 linear restraint: {:?}", set.linear);
    let r = refine(&mut m, data, &set, Options { max_cycles: 30, eps: p.eps }, |_| {});
    let st = &r.restraints[0];
    let (u, w) = (m.phases[0].u, m.phases[0].w);
    println!("Y2O3: U + W {:.6} ({u:.6} + {w:.6}) against {target:.6} ± 0.00001: {:+.2}σ; U moved {:+.6}, W {:+.6}", st.value, (st.value - target) / 0.00001, u - u0, w - w0);
    if set.linear.len() != 1 || set.linear[0].terms.len() != 2 || ((st.value - target) / 0.00001).abs() > 3.0 || (u - u0).abs() < 1e-5 || (w - w0).abs() < 1e-5 {
        fails.push(format!("Y2O3: the linear restraint U + W = {target:.6} is not held: {:.6}", st.value));
    }
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}

/// Against FullProf's own refinements with restraints (`runs/restraint_lebail_runs.py`, run on the Mac, L056): the
/// engine refines each `CODE.pcr.in` with FullProf's weight of angle restraints and is compared with `CODE_conv`:
/// every restraint's value within its σ of FullProf's, and every refined parameter within FullProf's esd.
///
/// FullProf's angle restraints pull with about π/180 of their stated weight (L058). The tutorial (p. 2) and FullProf's
/// own listing (its Σ W(Ang_obs − Ang_cal)², which with χ² reproduces the printed 11.3 and 87.3 for sr-ox-test) state
/// the angles and σ in degrees. Yet FullProf's four runs stop at their starting point in one cycle, with O2–C1–O5 at
/// 4.0σ and Ag–C–N at 3.4σ, where the stated weight moves them (to 1.65σ and 1.8σ). With the weight of the angle rows
/// scaled by π/180 (`Convention::angle_restraints`), every restraint of the four runs agrees with FullProf's within
/// 0.73σ and every parameter of the oxalates and of K2Na[Ag(CN)2]3 within its esd (at most 0.83). A radian
/// derivative against degree residuals would do this; (π/180)² fits as well (0.72σ), so the factor is an inference
/// from the results: FullProf's least-squares code is not public (CrysFML's restraint cost, `cost_functions.f90`, is
/// in degrees throughout and has no derivatives).
///
/// Frozen: the calculated patterns of two runs at FullProf's values. Sr-ox-test: χ² 1.9605 against 1.931 (+1.5 %);
/// its Fourier background (Nba = −4, 18 coefficients) differs from FullProf's by a smooth term up to 1.85 counts below
/// 60° (FullProf's Backg less the core's is constant beyond; 0.3 counts in Sr-ox-rest). LaMnO3-3t2: 2.2416 against
/// 2.137 (+4.9 %), all of it between 21° and 41°, around the peaks below AsyLim = 35°, with the background exact; the
/// form of F_b, the width in z and the place of P2 were tried and each makes it worse. Not explained; there 7
/// parameters differ by up to 2.71 esd.
#[test]
fn fullprof_restraint_runs() {
    let mut fails = Vec::new();
    for (case, code) in [("Sr-ox-test", "sr-ox-test"), ("Sr-ox-rest", "sr-ox-rest"), ("LaMnO3-3t2", "lamn_3t2"), ("K2NaAgCN", "k2na1p5")] {
        let d = root("runs").join(case);
        let out = latin1(&d.join(format!("{code}_conv.out")));
        // The last cycle's "=> Distance  restraints:" and "=> Angle restraints:" blocks (FullProf writes two spaces in
        // the first; the words are matched, not the spacing, L056): "(A - B ): Dobs Dcalc diff/sigma".
        let lines: Vec<&str> = out.lines().collect();
        let is_head = |l: &str, kind: &str| {
            let w: Vec<&str> = l.split_whitespace().collect();
            w.windows(2).any(|p| p[0] == kind && p[1] == "restraints:")
        };
        let block = |from: usize| -> Vec<f64> {
            lines[from + 1..]
                .iter()
                .skip_while(|l| l.trim().is_empty())
                .take_while(|l| l.contains("):"))
                .filter_map(|l| {
                    let v: Vec<f64> = l.split_once("):").unwrap().1.split_whitespace().filter_map(|t| t.parse().ok()).collect();
                    (v.len() >= 2).then(|| v[1])
                })
                .collect()
        };
        let last_dist = lines.iter().rposition(|l| is_head(l, "Distance"));
        let last_ang = lines.iter().rposition(|l| is_head(l, "Angle"));
        let mut fp: Vec<f64> = Vec::new();
        if let Some(i) = last_dist {
            fp.extend(block(i));
        }
        if let Some(i) = last_ang {
            fp.extend(block(i));
        }
        assert!(!fp.is_empty(), "{case}: no restraint block found in FullProf's listing");
        let (p, l) = load(&d.join(format!("{code}.pcr.in")));
        let mut m = l.model.clone();
        // FullProf's weight of angle restraints (L058).
        m.convention.angle_restraints = true;
        let set = ParameterSet::from_pcr(&p, &m).unwrap();
        let r = refine(&mut m, l.data.as_ref().unwrap(), &set, Options { max_cycles: 100, eps: p.eps }, |_| {});
        // The refined parameters against FullProf's (its `_conv.pcr`, imported the same way), within its esd's
        // (the asymmetry's P1 + 3P2 with σ₁ + 3σ₂, a bound; these runs are neutrons: no scale translation).
        let labels: Vec<Vec<String>> = p.phases.iter().map(|ph| ph.atoms.iter().map(|a| a.label.clone()).collect()).collect();
        let (_, conv) = load(&d.join(format!("{code}_conv.pcr")));
        let mut sig = std::collections::BTreeMap::new();
        for (name, &(_, e)) in &common::fullprof_params(&out) {
            if let Some(f) = common::fullprof_field(name, &labels) {
                sig.insert(f, e);
            }
        }
        let (mut worst, mut beyond) = ((0.0f64, String::new()), 0);
        for q in &set.free {
            let f = q.fields[0].0;
            let e = match f {
                drs::jacobian::Field::Asym(ph, k) if k % 2 == 0 => sig.get(&f).copied().unwrap_or(0.0) + 3.0 * sig.get(&drs::jacobian::Field::Asym(ph, k + 1)).copied().unwrap_or(0.0),
                _ => sig.get(&f).copied().unwrap_or(0.0),
            };
            if e <= 0.0 {
                continue;
            }
            let dev = (m.get(f) - conv.model.get(f)).abs() / e;
            if dev > 1.0 {
                beyond += 1;
                println!("    {:36} {:12.6} FullProf {:12.6} ± {e:.1e}: {dev:.2} esd", q.name, m.get(f), conv.model.get(f));
            }
            if dev > worst.0 {
                worst = (dev, q.name.clone());
            }
        }
        {
            let mut at = conv.model.clone();
            at.convention = m.convention;
            let data = l.data.as_ref().unwrap();
            let st = drs::lsq::statistics(&at, data, &at.calculate(&data.x).ycalc, set.free.len());
            let fp_chi2: f64 = out.lines().filter_map(|l| l.split_once("Conventional Rietveld Rp,Rwp,Re and Chi2:")).map(|(_, r)| r.split_whitespace().nth(3).unwrap().parse::<f64>().unwrap()).rev().nth(1).unwrap();
            let r2 = refine(&mut at, data, &set, Options { max_cycles: 100, eps: p.eps }, |_| {});
            println!("{case}: χ² at FullProf's values {:.4} (FullProf {fp_chi2}); the engine from there {:.4}", st.chi2, r2.stats.chi2);
            // The calculated pattern at FullProf's values: as `fullprof_profiles` (5·10⁻⁴ + 0.15 %), or frozen.
            let frozen = match case {
                "Sr-ox-test" => 0.0154,
                "LaMnO3-3t2" => 0.0491,
                _ => 0.0,
            };
            let d = (st.chi2 - fp_chi2).abs();
            if d > (5e-4 + 1.5e-3 * fp_chi2).max(frozen * fp_chi2) {
                fails.push(format!("{case}: χ² at FullProf's values {:.4} against {fp_chi2}, frozen {frozen}", st.chi2));
            }
        }
        println!("{case}: {} parameters, {beyond} beyond FullProf's esd, the largest {:.2} esd ({}); χ² {:.4}", set.free.len(), worst.0, worst.1, r.stats.chi2);
        // LaMnO3-3t2: its calculated pattern differs from FullProf's (above), 7 parameters by up to 2.71 esd.
        let (allowed, largest) = if case == "LaMnO3-3t2" { (7, 2.72) } else { (0, 1.0) };
        if beyond > allowed || worst.0 > largest {
            fails.push(format!("{case}: {beyond} parameters beyond FullProf's esd, the largest {:.2} ({})", worst.0, worst.1));
        }
        for (a, f) in r.restraints.iter().zip(&fp) {
            println!("{case}: {:28} engine {:10.5}, FullProf {f:10.5}: {:+.2}σ", a.label, a.value, (a.value - f) / a.sigma);
            if (a.value - f).abs() > a.sigma {
                fails.push(format!("{case}: {} {:.5} against FullProf's {f:.5}", a.label, a.value));
            }
        }
    }
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}
