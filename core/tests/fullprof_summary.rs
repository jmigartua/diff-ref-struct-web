//! The summary of a run (P57, `listing::summary`) against FullProf's `.sum` of the same calculation, on the 84 runs
//! of `fullprof_profiles` (each has its `.sum`; the 9 patterns of its joint runs left out, L322, L341): the core's summary at the `.pcr`'s values, with the `.pcr`'s refined
//! parameters counted in P as FullProf counts them.
//!
//! - **All the points** (FullProf's first block): N − P + C, Rp, Rwp, Rexp, χ², the conventional Rp, Rwp, Rexp, the
//!   deviance D, the Durbin–Watson d, its bound Q_D (eq. 3.20) and the goodness of fit in σ's, each within FullProf's
//!   print (three significant figures for the R's, χ² and D) and 0.2 %; the GoF within χ²'s tolerance times
//!   √((N − P + C)/2). Agreeing in all: 93 runs (L367). The 24 others are named in `ALL_DIFFER`, each with its known cause.
//! - **D*:** FullProf's "Dev*" is not defined in its manual; D/(N − P + C) gives its print in `DEV_STAR` runs, D/N in
//!   a few (Y₂O₃), neither in the rest: the core gives D/(N − P + C).
//! - **The points with Bragg contributions** (FullProf's second block) differ where FullProf's peak windows do (L042):
//!   reported, not compared. FullProf's d there divides the numerator of all the points by S over the Bragg points
//!   (`BRAGG_D` runs within 0.2 %): the core keeps the paper's d over the block's own points.
//! - **The Bérar factor** ("Scor", with the core's Ycal) within 0.05 % in `SCOR` runs (against FullProf's own Ycal it
//!   agrees in 60 runs: `berar_factor`); **R_Bragg and R_F** within 2 % in `R_PHASE` of the phases' values: FullProf's
//!   I_calc and I_obs are the intensities within its window, the core's the whole (`reflection_diagnostics`, L100).
//! - The summary and the listing as text (`summary_text`, `listing_text`) carry these values, and `drs results`
//!   writes them.
use drs::powder::Model;
use std::path::Path;

mod common;
use common::latin1;

fn root(rel: &str) -> String {
    format!("{}/../{rel}", env!("CARGO_MANIFEST_DIR"))
}

/// FullProf's numbers after `key` on the last line that holds it.
fn after(text: &str, key: &str) -> Vec<f64> {
    text.lines().rfind(|l| l.contains(key)).map_or(Vec::new(), |l| {
        l.split_once(key).unwrap().1.split(|c: char| c.is_whitespace() || c == '(' || c == ')' || c == ':').filter_map(|w| w.parse().ok()).collect()
    })
}

/// The two blocks of a `.sum`: [N−P+C, Rp, Rwp, Rexp, χ², cRp, cRwp, cRexp, D, D*, d, Q_D, GoF σ].
fn blocks(sum: &str) -> Vec<Vec<f64>> {
    let ls: Vec<&str> = sum.lines().collect();
    let starts: Vec<usize> = ls.iter().enumerate().filter(|(_, l)| l.contains("RELIABILITY FACTORS")).map(|(i, _)| i).collect();
    let take = |i: usize| -> Vec<f64> {
        let block = ls[i..(i + 14).min(ls.len())].join("\n");
        let mut v = after(&block, "N-P+C:");
        let rs: Vec<&str> = block.lines().filter(|l| l.contains("=> Rp:")).collect();
        for (k, r) in rs.iter().take(2).enumerate() {
            let n: Vec<f64> = r.replace("Rp:", " ").replace("Rwp:", " ").replace("Rexp:", " ").replace("Chi2:", " ").split_whitespace().filter_map(|w| w.parse().ok()).collect();
            v.extend(n.iter().take(4 - k));
        }
        v.extend(after(&block, "Deviance:").iter().take(1));
        v.extend(after(&block, "Dev*  :").iter().take(1));
        v.extend(after(&block, "DW-Stat.:").iter().take(1));
        v.extend(after(&block, "DW-exp:").iter().take(1));
        v.extend(after(&block, "N-sigma of the GoF:").iter().take(1));
        v
    };
    // The last two blocks: all the points, then those with Bragg contributions.
    starts.iter().rev().take(2).rev().map(|&i| take(i)).collect()
}

/// The runs whose all-points block differs from FullProf's .sum, with the cause.
const ALL_DIFFER: [(&str, &str); 24] = [
    ("srsclrrt.prf", "P counted differently from FullProf (the mode course's case A; N − P + C 2622 against 2649)"),
    ("a_coords_calc.prf", "P counted differently from FullProf (case A by coordinates)"),
    ("b_parent_1_calc.prf", "P counted differently from FullProf (case B)"),
    ("b_parent_2_calc.prf", "P counted differently from FullProf (case B)"),
    ("e_calc.prf", "FullProf counts the points of the second excluded region in N (L098), frozen in fullprof_profiles"),
    ("LaMnO3ref.prf", "inconsistent: its .prf is not of its .pcr (fullprof_profiles)"),
    ("YVO3a_calc.prf", "frozen in fullprof_profiles (L039: FullProf's window with X = 0.645)"),
    ("LBF_lmopr.prf", "d: FullProf's Le Bail intensities are one extraction from the .hkl's (fullprof_le_bail)"),
    ("citric-acid-sls_1.prf", "Rp and cRp just past the print (2.20 against 2.21)"),
    ("irf1_narrow_x0.645_wdt20/ybfluo_Riet_IRF.prf", "a window experiment (L042): χ² and D"),
    ("yb_K2_x0.645/ybfluo_Riet_IRF.prf", "a window experiment (L042): d"),
    ("yb_all_x0.645/ybfluo_Riet_IRF.prf", "a window experiment (L042): d"),
    ("yb_irf1_x0.645/ybfluo_Riet_IRF.prf", "a window experiment (L042): d"),
    ("pdfinal_calc.prf", "d (1.1064 against 1.11): the window edge at point 166 frozen in fullprof_profiles (L267)"),
    // P05's TOF runs (L341): each a difference frozen in fullprof_profiles.
    ("arg_si_calc.prf", "TOF: FullProf's window rule not found, 586 window-edge points frozen in fullprof_profiles (L341)"),
    ("arg_si_plain.prf", "TOF: the window-edge points, as arg_si_calc"),
    ("arg_si_plain90.prf", "TOF: d alone, the window-edge points, as arg_si_calc"),
    ("alum_calc.prf", "TOF: its χ² frozen in fullprof_profiles (0.7 %, the window edges)"),
    ("CeO2_PEARL_calc.prf", "TOF: FullProf's σ of a RALF bank 100 times larger again, its χ² frozen in fullprof_profiles"),
    ("lamn_pol_calc.prf", "TOF: Iwg = 1, FullProf's χ² not reproduced, frozen in fullprof_profiles"),
    ("lamn_pol_irf.prf", "TOF: Iwg = 1, as lamn_pol_calc"),
    ("nac-osiris_calc.prf", "TOF: Iwg = 1, FullProf's χ² not reproduced, frozen in fullprof_profiles"),
    // P30 (L360): the anvil cell's texture, FullProf's up to 3.6 % below the core's, frozen in fullprof_profiles.
    ("yb_nor3.prf", "Nor = 3: Wright's cone integral not reproduced, the profile frozen in fullprof_profiles and texture"),
    // P41's TLS (L367).
    ("urea_test_calc.prf", "P: its .pcr's 'Number of refined parameters' is 0 under codewords naming 17, and FullProf varies none (N − P + C 4734 against 4751); χ² and D: the TLS centre, frozen in fullprof_profiles"),
];

/// Counts as measured in L267 (L256's, with prba10k and Pd3MnD): runs where D/(N − P + C) gives FullProf's Dev*, where FullProf's Bragg-block d is the
/// all-points numerator over S of the Bragg points, where Scor agrees, and the phases' R_Bragg and R_F that agree, of
/// all of them. L341 added P05's 14 TOF runs (gabriela's three patterns left out as a joint run's): 6 agree in the
/// all-points block (arg_si_abs1, abs3, ext10, ext100, invd and alum_noq), D* in 4 more, the Bragg-block d in 5, Scor
/// in 1, R_Bragg and R_F in 2 of 30 (FullProf's TOF windows, within which it sums I_calc and I_obs, not found).
/// L353 added L346's 16 compared runs (14 of background_terms, ureaz_calc, PPH3-sls_calc): all 16 agree in the
/// all-points block (cRp of y2o3_poly12_b7, b9 and more2_two once the core divides by Σ(y − b) as it is, negative
/// where the background lies above the counts), D* in 3, the Bragg-block d in all 16, Scor in 14, R_Bragg and R_F in
/// 11 of 36.
/// L354 added L351's citric-acid-sls_1_calc and citric_core: both agree in the all-points block, D* and the Bragg-block
/// d; Scor in neither; R_Bragg and R_F in 2 of 4.
/// L360 added texture/'s four runs: yb_nor2, yb_nor2_g2 and yb_nor2_one agree in the all-points block, the Bragg-block d
/// and Scor, all four in D*; R_Bragg and R_F in 3 of 8.
/// L367 added rigid/'s eleven TLS runs: ten agree in the all-points block (urea_test_calc named in `ALL_DIFFER`), all
/// eleven in D*, ten in the Bragg-block d, seven in Scor; R_Bragg and R_F in 12 of 22.
const DEV_STAR: usize = 75;
const BRAGG_D: usize = 84;
const SCOR: usize = 50;
const R_PHASE: (usize, usize) = (94, 252);

/// Half a unit of FullProf's last printed digit: `sf` significant figures, or `dec` decimals.
fn print_unit(v: f64, sf: Option<i32>, dec: i32) -> f64 {
    match sf {
        Some(n) if v != 0.0 => 0.5 * 10f64.powi(v.abs().log10().floor() as i32 - (n - 1)),
        _ => 0.5 * 10f64.powi(-dec),
    }
}

fn short(prf: &str) -> String {
    let mut parts = prf.rsplit('/');
    let file = parts.next().unwrap();
    match (parts.next(), parts.next()) {
        (Some(dir), Some("window_experiments")) => format!("{dir}/{file}"),
        _ => file.to_string(),
    }
}

#[test]
fn fullprof_summary() {
    let index = std::fs::read_to_string(root("core/tests/data/fullprof_profiles.txt")).unwrap();
    let (mut runs, mut differ, mut dev_star, mut bragg_d, mut scor, mut r_ok, mut r_all) = (0, Vec::new(), 0, 0, 0, 0, 0);
    let mut joint = 0;
    let mut fails = Vec::new();
    for line in index.lines().filter(|l| l.starts_with("COMPARE")) {
        let f: Vec<&str> = line.split_whitespace().collect();
        let (prf, pcr) = (f[1], f[2]);
        // A joint run's pattern (CODE.pcr#k, L322): its .sum reports every pattern with N − P + C of FullProf's own
        // P per pattern, which the core does not reproduce (`joint`); left out, counted.
        if pcr.contains('#') {
            joint += 1;
            continue;
        }
        let sum = latin1(Path::new(&root(&prf.replace(".prf", ".sum"))));
        let p = drs::pcr::read(&latin1(Path::new(&root(pcr)))).unwrap();
        let path = Path::new(pcr);
        let dir = root(path.parent().unwrap().to_str().unwrap());
        let l = Model::from_pcr_fullprof(&p, Path::new(&dir), path.file_stem().unwrap().to_str().unwrap().trim_end_matches(".pcr")).unwrap();
        let data = l.data.unwrap();
        let mut m = l.model;
        // As fullprof_profiles calculates it: FullProf's reflection limits and step from its .out.
        let out = std::fs::read(root(&prf.replace(".prf", ".out"))).map(|b| b.iter().map(|&c| c as char).collect::<String>()).unwrap_or_default();
        let lim: Vec<f64> = out.lines().filter_map(|l| l.split_once("S(1/d)max:")).filter_map(|(_, r)| r.split_whitespace().next()?.parse().ok()).collect();
        if lim.len() == m.phases.len() {
            for (ph, v) in m.phases.iter_mut().zip(&lim) {
                ph.reflection_limit = Some(*v);
            }
        }
        m.fullprof_step = out.lines().find_map(|l| l.split_once("Step:")).and_then(|(_, r)| r.split_whitespace().next()?.parse().ok());
        let set = drs::refinable::ParameterSet::from_pcr(&p, &m).unwrap();
        let s = drs::listing::summary(&m, &data, &set, None);
        runs += 1;
        let fp = blocks(&sum);
        let (a, t) = (&s.all, &s.all.stats);
        let ours = [(t.n + t.c - t.p) as f64, t.rp, t.rwp, t.rexp, t.chi2, t.rp_background_free, t.rwp_background_free, t.rexp_background_free, a.deviance, a.deviance_per_dof, t.durbin_watson, a.durbin_watson_bound, a.gof_sigmas];
        let theirs = &fp[0];
        assert_eq!(theirs.len(), 13, "{prf}: the .sum's first block");
        // (significant figures, decimals) of FullProf's print for each item.
        let print = [(None, 0), (Some(3), 0), (Some(3), 0), (None, 2), (Some(3), 0), (Some(3), 0), (Some(3), 0), (None, 2), (Some(3), 0), (Some(4), 0), (None, 4), (None, 4), (None, 3)];
        let chi2_tol = print_unit(theirs[4], Some(3), 0) + 2e-3 * theirs[4];
        let mut bad = Vec::new();
        for k in (0..13).filter(|&k| k != 9) {
            let tol = if k == 12 { print_unit(theirs[k], None, 3) + chi2_tol * (theirs[0] / 2.0).sqrt() } else { print_unit(theirs[k], print[k].0, print[k].1) + 2e-3 * theirs[k].abs() };
            if (ours[k] - theirs[k]).abs() > tol {
                bad.push(format!("{} ({:.4} against {})", ["N−P+C", "Rp", "Rwp", "Rexp", "χ²", "cRp", "cRwp", "cRexp", "D", "D*", "d", "Q_D", "GoF"][k], ours[k], theirs[k]));
            }
        }
        let name = short(prf);
        if !bad.is_empty() {
            println!("{name}: {} differ from FullProf's .sum", bad.join(", "));
            differ.push(name.clone());
        }
        if (ours[9] - theirs[9]).abs() <= print_unit(theirs[9], Some(4), 0) + 2e-3 * theirs[9] {
            dev_star += 1;
        }
        // FullProf's d over the Bragg points: the all-points numerator over S of the Bragg points.
        let s_all = t.chi2 * ((t.n + t.c).saturating_sub(t.p)) as f64;
        let b = &s.bragg.stats;
        let s_b = b.chi2 * ((b.n + b.c).saturating_sub(b.p)) as f64;
        if fp.len() > 1 && fp[1].len() > 10 && (t.durbin_watson * s_all / s_b - fp[1][10]).abs() <= 5e-4 + 2e-3 * fp[1][10] {
            bragg_d += 1;
        }
        if let Some(&sc) = after(&sum, "Scor:").first() {
            if (s.berar_factor / sc - 1.0).abs() <= 5e-4 + 5e-5 / sc {
                scor += 1;
            }
        }
        let rb: Vec<f64> = sum.lines().filter_map(|l| l.split_once("Bragg R-factor:")).filter_map(|(_, r)| r.split_whitespace().next()?.parse().ok()).collect();
        let rf: Vec<f64> = sum.lines().filter_map(|l| l.split_once("Rf-factor=")).filter_map(|(_, r)| r.split_whitespace().next()?.parse().ok()).collect();
        for (k, ph) in s.phases.iter().enumerate() {
            for (o, t) in [(ph.r_bragg, rb.get(k)), (ph.r_f, rf.get(k))] {
                if let Some(&t) = t {
                    r_all += 1;
                    if (o - t).abs() <= 0.005 + 0.02 * t {
                        r_ok += 1;
                    }
                }
            }
        }
        // The text carries the values.
        let text = drs::listing::summary_text(&name, &s);
        for v in [format!("χ² {:.4}", t.chi2), format!("Bérar factor {:.4}", s.berar_factor), format!("Durbin–Watson d {:.4}", t.durbin_watson)] {
            if !text.contains(&v) {
                fails.push(format!("{name}: the summary's text lacks {v}"));
            }
        }
    }
    println!("{runs} runs: all-points block as FullProf's in {}; D* = D/(N − P + C) in {dev_star}; FullProf's Bragg-block d explained in {bragg_d}; Scor in {scor}; R_Bragg and R_F {r_ok} of {r_all}", runs - differ.len());
    let named: Vec<String> = ALL_DIFFER.iter().map(|d| d.0.to_string()).collect();
    let mut sorted = differ.clone();
    sorted.sort();
    let mut expected = named.clone();
    expected.sort();
    assert_eq!(sorted, expected, "the runs whose all-points block differs");
    assert!(fails.is_empty(), "{}", fails.join("\n"));
    assert_eq!((runs, joint, dev_star, bragg_d, scor, r_ok, r_all), (117, 9, DEV_STAR, BRAGG_D, SCOR, R_PHASE.0, R_PHASE.1), "the counts");
}

/// `drs results` writes the summary and the listing beside the results file, named by the file dictionary; the
/// listing holds each cycle and ends with the summary.
#[test]
fn summary_and_listing_written() {
    let src = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs/rutile-anatase");
    let d = Path::new(env!("CARGO_MANIFEST_DIR")).join("target/test-listing");
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    for f in ["rutana_tch.pcr", "rutana.dat", "rx_cukalpha.irf"] {
        std::fs::copy(src.join(f), d.join(f)).unwrap();
    }
    let time = "2026-10-15 12.00.00";
    let o = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(["results", d.join("rutana_tch.pcr").to_str().unwrap(), "--cycles", "2", "--time", time]).output().unwrap();
    let out = String::from_utf8_lossy(&o.stdout);
    assert!(o.status.success(), "{out}{}", String::from_utf8_lossy(&o.stderr));
    let file = |id: &str| drs::formats::files::path(id, &[("control", "rutana_tch"), ("time", time)]).unwrap();
    for id in ["run_summary", "run_listing"] {
        assert!(out.contains(&file(id)), "{id} not written:\n{out}");
        assert_eq!(drs::formats::files::kind_of(&file(id)).map(|k| k.id), Some(id));
    }
    let summary = std::fs::read_to_string(d.join(file("run_summary"))).unwrap();
    let listing = std::fs::read_to_string(d.join(file("run_listing"))).unwrap();
    for part in ["Refinement: 2 cycles", "Agreement, all the points used", "Agreement, the points with Bragg contributions", "Bérar factor", "Rutile", "Anatase", "(P 4_2/m n m)", "weight fraction"] {
        assert!(summary.contains(part), "the summary lacks {part}:\n{summary}");
    }
    assert!(listing.contains("Cycle 1:") && listing.contains("Cycle 2:") && listing.ends_with(&summary), "{listing}");
    // The parameters at each cycle (value, shift, esd and shift/esd) and once more in the summary.
    assert_eq!(listing.lines().filter(|l| l.trim_start().starts_with("phase1.scale_factor")).count(), 3, "{listing}");
}
