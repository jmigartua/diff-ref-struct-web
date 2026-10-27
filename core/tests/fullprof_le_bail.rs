//! Le Bail fits against FullProf's (test `fullprof_le_bail`, P35, L055).
//!
//! The intensities of a Jbt = 2 phase are extracted at every cycle (`Model::le_bail`, the Rietveld formula, manual
//! eq. 3.15) while the engine refines the other parameters; all cycles are run, as FullProf does with Dum = 1.
//!
//! - **LBF_lmopr** (LaMnO3, neutrons, TCH; FullProf's tutorial output `LBF_lmopr.out`): the engine starts from the
//!   values FullProf's `.out` lists before its first cycle and runs its 15 cycles; the intensities start from
//!   `LBF_lmopr1.hkl` (Irf = 2), which FullProf wrote at the end of that run (the ones it started from are not in the
//!   material). Checked: N − P is FullProf's N − P + C (2461: the intensities are not counted), every refined
//!   parameter ends within FullProf's esd, and the engine's χ² after 15 cycles is not above FullProf's. Frozen: at
//!   FullProf's final values and `.hkl` intensities the core's χ² is 1.4139 against FullProf's 1.420 (−0.43 %; one
//!   more extraction gives 1.4080: the `.hkl` is believed to hold the intensities one extraction later than those of
//!   FullProf's last χ², not checked); after 15 cycles 1.3630 against 1.420 (−4.0 %), the engine starting from
//!   intensities closer to the end and FullProf's χ² still falling by 0.006 per cycle.
//! - **MgAl2O4s_LeBail** (spinel, X-rays; no FullProf output, only its `.pcr`, refined by FullProf, whose header gives
//!   the "global user-weighted χ² (Bragg contrib.)" of its last run, 1.848): no parameter refined, ten extractions
//!   from equal intensities. FullProf's global χ² is S over all points divided by N − P + C over the points with Bragg
//!   contributions (Yb2ZrWO8_conv.out: 5357/3114 = 1.720, printed 1.72; Egypt: 28130/6723 = 4.184, printed 4.18).
//!   The core: 1.8516 against 1.848 (+0.20 %, frozen; FullProf's starting intensities and Bragg points not known).
//! - **Cocuc_LeBail** (Co-YBCO, neutrons, pseudo-Voigt): its header (26.69) belongs to an earlier run than its values
//!   (round starting values, U = 0.086, a = 3.869); the engine's nine cycles are reported, not compared. FullProf runs
//!   of these two files are listed for Josu (L055).

use drs::lsq::{refine, statistics, Options};
use drs::powder::Model;
use drs::refinable::ParameterSet;
use std::path::Path;

mod common;
use common::{fullprof_field, fullprof_params, latin1};

fn root(rel: &str) -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof").join(rel)
}

/// FullProf's starting values: the "Parameter number … -> Symbolic Name: NAME value" lines before its first cycle.
fn fullprof_start(out: &str) -> Vec<(String, f64)> {
    let head = &out[..out.find("CYCLE No.").unwrap_or(out.len())];
    head.lines()
        .filter_map(|l| l.split_once("-> Symbolic Name:"))
        .filter_map(|(_, r)| {
            let mut t = r.split_whitespace();
            Some((t.next()?.to_string(), t.next()?.parse().ok()?))
        })
        .collect()
}

/// FullProf's global user-weighted χ² (Bragg contrib.): S over all the points used, over N − P of the points that
/// some reflection's window covers.
fn global_bragg_chi2(m: &Model, data: &drs::data::Pattern, p: usize) -> (f64, usize) {
    let c = m.calculate(&data.x);
    let ex = m.excluded_points(&data.x);
    let mut bragg = vec![false; data.x.len()];
    for pk in &c.peaks {
        m.peak_profile(&data.x, &ex, pk, |i, _| bragg[i] = true);
    }
    let st = statistics(m, data, &c.ycalc, p);
    let used = |i: usize| !ex[i] && data.sigma[i] > 0.0 && !m.excluded.iter().any(|&(a, b)| data.x[i] >= a && data.x[i] <= b);
    let nb = (0..data.x.len()).filter(|&i| bragg[i] && used(i)).count();
    (st.chi2 * (st.n - st.p) as f64 / (nb - p) as f64, nb)
}

#[test]
fn fullprof_le_bail() {
    let mut fails = Vec::new();

    // LBF_lmopr.
    let dir = root("7-FullProf-Tutorial-SymmetryModes/LaMnO3");
    let out = latin1(&dir.join("LBF_lmopr.out"));
    let p = drs::pcr::read(&latin1(&dir.join("LBF_lmopr.pcr"))).unwrap();
    let loaded = Model::from_pcr(&p, &dir, "LBF_lmopr").unwrap();
    let data = loaded.data.as_ref().unwrap();
    let mut m = loaded.model.clone();
    let lim: Vec<f64> = out.lines().filter_map(|l| l.split_once("S(1/d)max:")).map(|(_, r)| r.split_whitespace().next().unwrap().parse().unwrap()).collect();
    m.phases[0].reflection_limit = lim.last().copied();
    let fp_chi2: Vec<f64> = out
        .lines()
        .filter_map(|l| l.split_once("Conventional Rietveld Rp,Rwp,Re and Chi2:"))
        .map(|(_, r)| r.split_whitespace().nth(3).unwrap().parse().unwrap())
        .step_by(2)
        .collect();
    let fp_last = *fp_chi2.last().unwrap();
    let at = statistics(&m, data, &m.calculate(&data.x).ycalc, 8);
    println!("LBF_lmopr at FullProf's last values and intensities: χ² {:.4}, FullProf {fp_last}", at.chi2);
    if (at.chi2 - fp_last).abs() / fp_last > 0.0044 {
        fails.push(format!("LBF_lmopr: χ² at FullProf's values {:.4} against {fp_last}, frozen 0.43 %", at.chi2));
    }
    for (name, v) in fullprof_start(&out) {
        let f = fullprof_field(&name, &[vec![]]).unwrap_or_else(|| panic!("{name}"));
        m.set(f, v);
    }
    let set = ParameterSet::from_pcr(&p, &m).unwrap();
    let r = refine(&mut m, data, &set, Options { max_cycles: fp_chi2.len(), eps: p.eps }, |c| {
        println!("LBF_lmopr cycle {:2}: χ² {:.4}, FullProf {}", c.cycle, c.stats.chi2, fp_chi2[c.cycle - 1]);
    });
    let dof: usize = out.lines().filter_map(|l| l.split_once("N-P+C:")).map(|(_, r)| r.trim().parse().unwrap()).next().unwrap();
    if r.stats.n - r.stats.p != dof {
        fails.push(format!("LBF_lmopr: N − P {} against FullProf's {dof}", r.stats.n - r.stats.p));
    }
    let d = (r.stats.chi2 - fp_last) / fp_last;
    if !(-0.041..=0.0).contains(&d) {
        fails.push(format!("LBF_lmopr: χ² after {} cycles {:.4} against FullProf's {fp_last}: {d:+.4}, frozen −0.040", r.cycles, r.stats.chi2));
    }
    let fpp = fullprof_params(&out);
    for (fp, ours) in set.free.iter().zip(&r.params) {
        let f = fp.fields[0].0;
        let Some((n, &(v, e))) = fpp.iter().find(|(n, _)| fullprof_field(n, &[vec![]]) == Some(f)) else {
            fails.push(format!("LBF_lmopr: {} has no FullProf value", fp.name));
            continue;
        };
        let dev = (ours.value - v).abs() / e;
        println!("  {:28} {:12.6} ({:.1e}) FullProf {n} {v:12.6} ({e:.1e}): {dev:.2} esd", fp.name, ours.value, ours.esd);
        if dev > 1.0 {
            fails.push(format!("LBF_lmopr: {} {:.6} against FullProf's {v} ± {e}: {dev:.2} esd", fp.name, ours.value));
        }
    }

    // MgAl2O4s_LeBail and Cocuc_LeBail.
    for (rel, code, data_rel, frozen) in [
        ("4-Analysis-Spinel-MgAl2O4/pcr-files", "MgAl2O4s_LeBail", "4-Analysis-Spinel-MgAl2O4/Data-files/MgAl2O4s.dat", Some(0.0020)),
        ("3-Tutorial_Co_doped_YBCO/pcr-files", "Cocuc_LeBail", "3-Tutorial_Co_doped_YBCO/Cocuc.dat", None),
    ] {
        let dir = root(rel);
        let text = latin1(&dir.join(format!("{code}.pcr")));
        let p = drs::pcr::read(&text).unwrap();
        let loaded = Model::from_pcr_with_data(&p, &dir, code, Some(&root(data_rel))).unwrap();
        let data = loaded.data.as_ref().unwrap();
        let mut m = loaded.model.clone();
        let set = ParameterSet::from_pcr(&p, &m).unwrap();
        let r = refine(&mut m, data, &set, Options { max_cycles: p.ncy as usize, eps: p.eps }, |_| {});
        let header: f64 = text.lines().nth(1).and_then(|l| l.split('=').nth(1)).and_then(|v| v.trim().parse().ok()).unwrap();
        let (g, nb) = global_bragg_chi2(&m, data, set.free.len());
        println!("{code}: {} cycles, {} parameters: χ² {:.4}; global (Bragg contrib.) {g:.4} over {nb} − P, FullProf's .pcr header {header}", r.cycles, set.free.len(), r.stats.chi2);
        for q in &r.params {
            println!("  {:36} {:12.6} ({:.1e})", q.name, q.value, q.esd);
        }
        if let Some(f) = frozen {
            if (g - header).abs() / header > f {
                fails.push(format!("{code}: global χ² {g:.4} against FullProf's {header}, frozen {f}"));
            }
        }
    }
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}
