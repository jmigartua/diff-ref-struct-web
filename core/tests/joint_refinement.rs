//! Joint refinement (P06, L322) against FullProf's three joint runs of `manuals/fullprof/runs/` (joint_runs.py, L316,
//! run on the Mac in L320): pbb (PbSO₄, D1A neutrons + Cu Kα X-rays, one phase), scor_rxn (Sr₂Co₂O₅ and CoO, D1A + Cu
//! Kα, both phases in both patterns), prbaco (PrBaCo₂O₅₊δ, BM16 synchrotron + D2B neutrons, background points).
//!
//! - `split_and_agreement`: each `.pcr` split into its patterns (`pcr::split_patterns`); at the values of CODE.pcr,
//!   which CODE_calc.out reports, each pattern's Σw(y − y_c)² over all its points and over those with Bragg
//!   contributions against the "Residual" FullProf prints for that pattern (4 digits), and the global χ², the patterns'
//!   Bragg χ² over FullProf's N − P + C weighted by W_PAT, against FullProf's "Global user-weigthed Chi2 (Bragg
//!   contrib.)" (3 digits). The patterns whose calculation is frozen in `fullprof_profiles` are frozen here by their
//!   relative difference.
//! - `refinement_against_fullprof`: pbb and scor_rxn refined by the joint engine from CODE.pcr.in, where FullProf's own
//!   run converged in its first cycle (its shifts are 0 to the print: pbb.out, lines 464–520): the engine converges; its
//!   minimum of Q = Σ W_k S_k is not above Q at FullProf's values; each parameter against FullProf's value in FullProf's
//!   esd's, and the engine's esd's by FullProf's rule (`Esd::FullProf`) against FullProf's (q = σ/σ_FullProf); what
//!   departs is counted and frozen with its reason.
//! - `prbaco_where_fullprof_diverged`: FullProf's refinement of prbaco rose from χ² 11.4 to 13.4 and ended its five
//!   cycles at 9.32 without saving them ("final(Chi2) > 1.2*initial(Chi2)": prbaco.run.log); the engine lowers the
//!   objective every cycle and is below FullProf's after the same five.
//! - `one_pattern_is_the_single_engine`: a joint refinement of one pattern (a single-pattern `.pcr`) reaches the
//!   single-pattern engine's values, esd's and χ².
//! - `results_file_as_publish_reads_it` (P41): the joint results file of scor_rxn (two phases in two patterns) read
//!   back (`joint::read_results`): a block per pattern, `|pattern k`, with its weight and agreement as the refinement
//!   left them; every reflection with its pattern and its joint phase; the phases' blocks and contributions; the
//!   overall block's joint χ², R factors and parameters.
//! - `models_per_pattern_without_the_data` (Publish's P50): each pattern's model, phases and free parameters read from
//!   the joint results file alone, in an empty folder, as the refinement left them; from a file written before core
//!   0.61.0, from its overall block's `.pcr` (in its folder) or its joint control file's patterns' files.

use drs::joint::{Esd, Joint};
use drs::lsq::Options;
use std::path::Path;

fn root(rel: &str) -> String {
    format!("{}/../{rel}", env!("CARGO_MANIFEST_DIR"))
}

fn latin1(path: &str) -> String {
    std::fs::read(path).unwrap_or_else(|e| panic!("{path}: {e}")).iter().map(|&b| b as char).collect()
}

fn joint(code: &str, file: &str) -> Joint {
    let dir = root(&format!("manuals/fullprof/runs/{code}"));
    Joint::from_pcr(&latin1(&format!("{dir}/{file}")), Path::new(&dir), code, true).unwrap_or_else(|e| panic!("{code}/{file}: {e}"))
}

/// FullProf's statistics of pattern k in the first cycle of an `.out`: (Residual over all the points, N − P + C), and
/// the same over the points with Bragg contributions.
fn fullprof_residuals(out: &str, k: usize) -> ((f64, usize), (f64, usize)) {
    let block = |head: &str| -> (f64, usize) {
        let start = out.find(&format!("{head}{k:>3}")).unwrap_or_else(|| panic!("no '{head}' for pattern {k}"));
        let lines: Vec<&str> = out[start..].lines().take(16).collect();
        let npc = lines.iter().find_map(|l| l.split_once("N-P+C:")).map(|(_, r)| r.split_whitespace().next().unwrap().parse().unwrap()).unwrap();
        let i = lines.iter().position(|l| l.contains("SumYdif")).unwrap();
        (lines[i + 1].split_whitespace().nth(4).unwrap().parse().unwrap(), npc)
    };
    (block("RELIABILITY FACTORS WITH ALL NON-EXCLUDED POINTS FOR PATTERN:"), block("RELIABILITY FACTORS FOR POINTS WITH BRAGG CONTRIBUTIONS FOR PATTERN:"))
}

/// The patterns whose calculation `fullprof_profiles` freezes (L322), by their relative difference in Σw(y − y_c)²
/// (all points, Bragg points): scor_rxn's X-ray pattern (every strong point 0.19% above FullProf's); prbaco's
/// synchrotron pattern (its window edges with Wdt 35, its strongest peaks within 1.2·10⁻⁴); prbaco's neutron pattern,
/// whose Σ here, the reflections not cut at the `.out`'s S(1/d)max, is 0.13% below FullProf's, where its `.prf`'s Ycal
/// (the cut applied, as `fullprof_profiles` applies it) gives 0.9% above. And prbaco's global χ², 11.329 against the
/// 11.4 printed, from its synchrotron pattern's 0.26%.
const RESIDUAL_FROZEN: &[(&str, usize, f64, f64)] = &[("scor_rxn", 2, 5.64e-4, 5.60e-4), ("prbaco", 1, 2.90e-3, 2.62e-3), ("prbaco", 2, 1.34e-3, 2.7e-4)];
const GLOBAL_FROZEN: &[(&str, f64)] = &[("prbaco", 0.0713)];

#[test]
fn split_and_agreement() {
    let mut fails = Vec::new();
    for (code, n, contributions, fp_global) in [("pbb", 2, vec![vec![true, true]], 3.17), ("scor_rxn", 2, vec![vec![true, true], vec![true, true]], 3.01), ("prbaco", 2, vec![vec![true, true]], 11.4)] {
        let m = drs::pcr::split_patterns(&latin1(&root(&format!("manuals/fullprof/runs/{code}/{code}.pcr")))).unwrap();
        assert_eq!((m.texts.len(), m.weights.clone(), m.refined.clone(), m.contributions.clone()), (n, vec![0.5, 0.5], vec![true, true], contributions), "{code}: the split");
        let j = joint(code, &format!("{code}.pcr"));
        let out = latin1(&root(&format!("manuals/fullprof/runs/{code}/{code}_calc.out")));
        let mut global = 0.0;
        for (k, p) in j.patterns.iter().enumerate() {
            let c = p.model.calculate(&p.data.x);
            let w = drs::lsq::weights(&p.model, &p.data, &c.ycalc);
            let bragg = drs::lsq::bragg_points(&p.data, &c.peaks);
            let (mut s, mut sb) = (0.0, 0.0);
            for i in 0..w.len() {
                let r = w[i] * (p.data.y[i] - c.ycalc[i]).powi(2);
                s += r;
                if bragg[i] {
                    sb += r;
                }
            }
            let ((fs, _), (fsb, npc_b)) = fullprof_residuals(&out, k + 1);
            global += p.weight * sb / npc_b as f64;
            let (ds, dsb) = ((s - fs) / fs, (sb - fsb) / fsb);
            println!("{code} pattern {} ({}): Σ {s:.1} against FullProf's {fs} ({ds:+.2e}); Bragg points {sb:.1} against {fsb} ({dsb:+.2e})", k + 1, p.label);
            let (fa, fb) = RESIDUAL_FROZEN.iter().find(|f| f.0 == code && f.1 == k + 1).map_or((5e-4, 5e-4), |f| (f.2.max(5e-4), f.3.max(5e-4)));
            if ds.abs() > fa || dsb.abs() > fb {
                fails.push(format!("{code} pattern {}: Σ {ds:+.4e}, Bragg {dsb:+.4e}", k + 1));
            }
        }
        println!("{code}: global χ² (Bragg) {global:.4} against FullProf's {fp_global}");
        let tol = GLOBAL_FROZEN.iter().find(|f| f.0 == code).map_or(5e-3 * fp_global, |f| f.1);
        if (global - fp_global).abs() > tol {
            fails.push(format!("{code}: global χ² {global:.4} against {fp_global}"));
        }
    }
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}

/// The parameters more than 0.25 of FullProf's esd from FullProf's value, frozen (count, largest), and the esd's with
/// q = σ/σ_FullProf outside 0.98–1.02, the cells of the TCH phases apart (count, largest |q − 1|), with their reasons
/// (L322):
/// - pbb: the engine's minimum is lower than FullProf's end (Q 13067.2 against 13134.9), where both calculate the
///   pattern alike (pbb_calc_1/2.prf exact, Σ to 4 digits); FullProf stopped where its own sum is not at its minimum,
///   its shifts 0 from the start (believed: FullProf's derivatives of the X-ray pattern's widths, U V W of pattern 2
///   moving 2.2–3.8 esd); the esd's of the X-ray pattern's own parameters 3–13% off FullProf's.
/// - scor_rxn: the calculation of its X-ray pattern differs (frozen in `fullprof_profiles`), and phase 1's scale is one
///   parameter for both patterns (codeword 11 in both: the same shift for scales of 9.2·10⁻⁴ and 5.8·10⁻⁷), which takes
///   up none of it: the minima differ (Q 7761.3 against 8050.1 at FullProf's values; global χ² 2.897 against 3.018),
///   the zeros (1.4 and 4.2 esd), the cells (1.5–2.9) and the X-ray widths (W 3.6) the most; the esd's of the X-ray
///   pattern's own parameters 3–14% off FullProf's.
///
/// The esd's compared are scaled to one global χ² (FullProf's rule carries that of the minimum); scor_rxn's S_400 of
/// strain model 8 against FullProf's Saa as 2·Saa·σ(Saa).
const PARAMS_FROZEN: &[(&str, usize, f64, usize, f64)] = &[("pbb", 12, 3.823, 7, 0.1310), ("scor_rxn", 17, 4.211, 14, 0.1436)];

/// FullProf's esd's by parameter number, and its values in CODE.pcr (the joint refinement built from it).
fn fullprof_esd(code: &str) -> std::collections::BTreeMap<u32, f64> {
    let out = latin1(&root(&format!("manuals/fullprof/runs/{code}/{code}.out")));
    out.lines()
        .filter(|l| l.contains("+/-") && l.contains("Parameter number"))
        .map(|l| {
            let n = l.split("Parameter number").nth(1).unwrap().split(':').next().unwrap().trim().parse().unwrap();
            let e = l.split("+/-").nth(1).unwrap().trim().trim_end_matches(')').trim().parse().unwrap();
            (n, e)
        })
        .collect()
}

#[test]
fn refinement_against_fullprof() {
    let mut fails = Vec::new();
    for code in ["pbb", "scor_rxn"] {
        let fp = joint(code, &format!("{code}.pcr"));
        let fp_values = fp.values();
        let fp_q = fp.objective();
        let esd = fullprof_esd(code);
        let mut j = joint(code, &format!("{code}.pcr.in"));
        let r = j.refine(Options { max_cycles: 30, eps: 0.01 }, Esd::FullProf, |_| true);
        println!("{code}: {} ({} cycles); Q {:.3} against {fp_q:.3} at FullProf's values; global χ² {:.4}", r.stop, r.cycles, j.objective(), r.statistics.chi2);
        if !r.converged {
            fails.push(format!("{code}: {}", r.stop));
        }
        if j.objective() > fp_q * (1.0 + 1e-6) {
            fails.push(format!("{code}: the engine's Q {} is above Q at FullProf's values {fp_q}", j.objective()));
        }
        // The esd's by FullProf's rule carry the global χ² of the minimum they belong to: compared at one χ².
        let scaling = (fp.statistics().chi2 / r.statistics.chi2).sqrt();
        let (mut beyond, mut worst, mut off, mut worst_q) = (0, 0.0f64, 0, 0.0f64);
        for (i, (name, v, e)) in r.params.iter().enumerate() {
            let n = j.parameters[i].number.expect("a FullProf number");
            // FullProf's strain model 8 refines Saa, the engine S_400 = Saa² (`powder::hexagonal_strain`).
            let fe = if name.ends_with("quartic_strain.400") { esd[&n] * 2.0 * fp_values[i].sqrt() } else { esd[&n] };
            let (dev, q) = ((v - fp_values[i]) / fe, e / fe * scaling);
            let tch_cell = name.contains("cell_");
            println!("    {name:44} {v:14.7} ({e:.2e}) FullProf {:14.7} ({fe:.2e}): {dev:+.3} esd, q {q:.4}", fp_values[i]);
            if dev.abs() > 0.25 {
                beyond += 1;
                worst = worst.max(dev.abs());
            }
            let q_bad = if tch_cell { !(1.93..=2.02).contains(&q) } else { (q - 1.0).abs() > 0.02 };
            if q_bad {
                off += 1;
                worst_q = worst_q.max(if tch_cell { (q / 2.0 - 1.0).abs() } else { (q - 1.0).abs() });
            }
        }
        let f = PARAMS_FROZEN.iter().find(|f| f.0 == code).unwrap();
        println!("{code}: {beyond} parameters beyond 0.25 esd (largest {worst:.3}); {off} esd's off FullProf's (largest |q − 1| {worst_q:.4})");
        if beyond > f.1 || worst > f.2 || off > f.3 || worst_q > f.4 {
            fails.push(format!("{code}: {beyond} beyond 0.25 esd up to {worst:.3}, {off} esd's off up to {worst_q:.4}; frozen {} up to {}, {} up to {}", f.1, f.2, f.3, f.4));
        }
    }
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}

#[test]
fn prbaco_where_fullprof_diverged() {
    let mut j = joint("prbaco", "prbaco.pcr");
    let start = j.objective();
    let mut qs = vec![start];
    let r = j.refine(Options { max_cycles: 5, eps: 0.01 }, Esd::Default, |c| {
        qs.push(c.objective);
        true
    });
    println!("prbaco: Q {:?}; global χ² {:.4} after {} cycles (FullProf's 11.4, 13.4, 11.3, 10.1, 9.32)", qs, r.statistics.chi2, r.cycles);
    assert!(qs.windows(2).all(|w| w[1] <= w[0]), "the objective rises");
    assert!(r.statistics.chi2 < 9.32, "global χ² {} after five cycles, FullProf's 9.32", r.statistics.chi2);
}

#[test]
fn one_pattern_is_the_single_engine() {
    let dir = root("manuals/fullprof/runs/Y2O3");
    let text = latin1(&format!("{dir}/y2o3_structure.pcr"));
    let mut j = Joint::from_pcr(&text, Path::new(&dir), "y2o3_structure", false).unwrap();
    let opts = Options { max_cycles: 3, eps: 1e-9 };
    let r = j.refine(opts, Esd::Default, |_| true);
    let p = drs::pcr::read(&text).unwrap();
    let mut l = drs::powder::Model::from_pcr(&p, Path::new(&dir), "y2o3_structure").unwrap();
    let set = drs::refinable::ParameterSet::from_pcr(&p, &l.model).unwrap();
    let s = drs::lsq::refine(&mut l.model, l.data.as_ref().unwrap(), &set, opts, |_| {});
    assert_eq!(r.params.len(), s.params.len());
    for (a, b) in r.params.iter().zip(&s.params) {
        assert_eq!(a.0, b.name);
        assert!((a.1 - b.value).abs() <= 1e-9 * b.value.abs().max(1e-6), "{}: {} against {}", a.0, a.1, b.value);
        assert!((a.2 - b.esd).abs() <= 1e-6 * b.esd, "{}: esd {} against {}", a.0, a.2, b.esd);
    }
    assert!((r.statistics.patterns[0].all.chi2 - s.stats.chi2).abs() < 1e-9 * s.stats.chi2);
}

#[test]
fn results_file_as_publish_reads_it() {
    let mut j = joint("scor_rxn", "scor_rxn.pcr");
    let r = j.refine(Options { max_cycles: 1, eps: 0.01 }, Esd::Default, |_| true);
    let name = "scor_rxn 2026-10-18 12.00.00";
    let text = j.results(Some(&r), name, &latin1(&root("manuals/fullprof/runs/scor_rxn/scor_rxn.pcr"))).unwrap();
    let back = drs::joint::read_results(&text).unwrap();
    let close = |a: f64, b: f64| (a - b).abs() <= 1e-6 * b.abs().max(1e-9);
    assert!(close(back.chi2, r.statistics.chi2_all) && close(back.chi2_bragg, r.statistics.chi2) && close(back.rwp, r.statistics.rwp) && close(back.objective, j.objective()));
    assert_eq!(back.parameters, j.parameters.len());
    assert_eq!(back.refined.len(), j.parameters.len());
    for (a, b) in back.refined.iter().zip(&r.params) {
        assert!(a.0 == b.0 && close(a.1, b.1) && close(a.2, b.2), "{a:?} against {b:?}");
    }
    assert_eq!(back.phase_blocks, vec![(1, format!("{name}|phase 1")), (2, format!("{name}|phase 2"))]);
    assert_eq!(back.patterns.len(), 2);
    for (k, p) in back.patterns.iter().enumerate() {
        let s = &r.statistics.patterns[k];
        assert_eq!((p.id, p.block.clone(), p.weight, p.refined, p.phases.clone()), (k + 1, format!("{name}|pattern {}", k + 1), 0.5, true, vec![1, 2]));
        assert!(close(p.rp, s.all.rp) && close(p.rwp, s.all.rwp) && close(p.rexp, s.all.rexp) && close(p.chi2, s.all.chi2) && close(p.chi2_bragg, s.bragg.chi2), "pattern {}", k + 1);
        assert_eq!(p.points, j.patterns[k].data.x.len());
        assert!(!p.reflections.is_empty() && p.reflections.iter().all(|&(ph, pat)| pat == k + 1 && (ph == 1 || ph == 2)), "pattern {}: reflections", k + 1);
        assert!(p.reflections.iter().any(|r| r.0 == 2), "pattern {}: CoO's reflections", k + 1);
    }
    assert_eq!(back.patterns[0].label, "scn");
}

#[test]
fn models_per_pattern_without_the_data() {
    // P50: each pattern's model from the joint results file alone, in a folder that holds nothing.
    let mut j = joint("scor_rxn", "scor_rxn.pcr");
    let r = j.refine(Options { max_cycles: 1, eps: 0.01 }, Esd::Default, |_| true);
    let pcr = latin1(&root("manuals/fullprof/runs/scor_rxn/scor_rxn.pcr"));
    let text = j.results(Some(&r), "scor_rxn 2026-10-19 12.00.00", &pcr).unwrap();
    let back = drs::joint::read_results(&text).unwrap();
    assert_eq!(back.name, "scor_rxn 2026-10-19 12.00.00");
    let empty = Path::new(env!("CARGO_MANIFEST_DIR")).join("target/test-joint-models-empty");
    let _ = std::fs::remove_dir_all(&empty);
    std::fs::create_dir_all(&empty).unwrap();
    let cs = back.pattern_controls(&empty).unwrap();
    assert_eq!(cs.len(), 2);
    for (k, c) in cs.iter().enumerate() {
        assert_eq!(c.model.to_text(), j.patterns[k].model.to_text(), "pattern {}: the refined model", k + 1);
        assert_eq!(c.model.phases.len(), back.patterns[k].phases.len(), "pattern {}: its phases, named by the joint ids", k + 1);
        assert_eq!(c.refinement, j.patterns[k].refinement, "pattern {}: its free parameters", k + 1);
    }
    // A file written before 0.61.0: no pattern control files, the models from the overall block's .pcr in its folder.
    let mut old = back.clone();
    for p in &mut old.patterns {
        p.control_file.clear();
    }
    let dir = root("manuals/fullprof/runs/scor_rxn");
    let fresh = Joint::from_pcr(&pcr, Path::new(&dir), "scor_rxn", false).unwrap();
    let cs = old.pattern_controls(Path::new(&dir)).unwrap();
    assert_eq!(cs.len(), 2);
    for (k, c) in cs.iter().enumerate() {
        assert_eq!(c.model.to_text(), fresh.patterns[k].model.to_text(), "pattern {}: the .pcr's model", k + 1);
    }
    // ... or from a joint control file's patterns' control files.
    let saved = Path::new(env!("CARGO_MANIFEST_DIR")).join("target/test-joint-models-saved");
    let _ = std::fs::remove_dir_all(&saved);
    std::fs::create_dir_all(&saved).unwrap();
    j.save(&saved.join("scor_rxn.drc")).unwrap();
    old.control_file = std::fs::read_to_string(saved.join("scor_rxn.drc")).unwrap();
    let cs = old.pattern_controls(&saved).unwrap();
    for (k, c) in cs.iter().enumerate() {
        assert_eq!(c.model.to_text(), j.patterns[k].model.to_text(), "pattern {}: from its control file", k + 1);
    }
    // Neither: refused.
    old.control_file.clear();
    assert!(old.pattern_controls(&empty).is_err());
    // `drs joint <results.cif>`: each pattern's control file, or pattern k's.
    std::fs::write(empty.join("scor_rxn results.cif"), &text).unwrap();
    let run = |extra: &[&str]| {
        let o = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).arg("joint").arg(empty.join("scor_rxn results.cif")).args(extra).output().unwrap();
        assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
        String::from_utf8(o.stdout).unwrap()
    };
    let all = run(&[]);
    assert_eq!(all.matches("# pattern ").count(), 2);
    let second = run(&["--pattern", "2"]);
    assert!(second.starts_with("# pattern 2 (scx), the joint phases 1, 2\n"), "{second}");
    assert!(second.contains(&back.pattern_controls(&empty).unwrap()[1].to_text()));
}
