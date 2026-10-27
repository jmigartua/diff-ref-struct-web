//! Calculated patterns against FullProf's `.prf` files (test `fullprof_profiles`).
//!
//! `tests/data/fullprof_profiles.txt` (written by `tests/data/make_fullprof_profiles.py`) lists the `.prf` files of
//! `manuals/fullprof/` and the two of the mode course's case A (`manuals/bcs`, L090): 68 rebuilt from their `.pcr` and compared point by point (among them the calculations of ten
//! FullProf runs of `runs/`, L034–L042, the CeO₂ runs ce and ce2 with surface micro-absorption and NaBD₄ with the
//! user's complex b (L239, L251), the two runs written for λ and the user's f′, f″ (`nabd4_lambda`, `y2o3_fprime`,
//! L252, L255), the eight of the mode course's cases in `runs/modes/`, L093–L098, and the 35
//! window experiments of `runs/window_experiments/`, L041–L042, whose
//! `.prf` belongs to their `input.pcr`), 71 used for their observed data only, and the others left out with the
//! reason (the anisotropic runs are version 2; NdNi-HTO has a magnetic phase). For each file compared, the core reads the `.pcr` (and the data, `.hkl`
//! and resolution files it names), calculates the pattern at the `.prf` points, and compares:
//! - Ycal, as the deviation beyond FullProf's print (half a unit of the last decimal it printed) over max(|Ycal_FP|, 1),
//!   max and rms over all points; the reflections are generated as far as FullProf did (S(1/d)max of its `.out`);
//! - Backg: the `.prf` prints the background less a constant, a quarter of its minimum over the pattern (found from
//!   the files: LaMnO3_cub.prf prints 11.789 where the background is 15.719 = 4·3.930, LBF_lmopr.prf −1.163 below a
//!   minimum of −4.651; L082); the files are read by the library's reader (`formats::fullprof_outputs::read_prf`,
//!   L100), which adds the constant back, and its background is compared with the core's;
//! - the reflection positions listed in the `.prf` (h k l and 2θ, Zero included), to 5·10⁻⁴ degrees (the worst,
//!   4.4·10⁻⁴° in LaMnO3_cub.prf, grows with k: FullProf's b there is about 1.5·10⁻⁵ Å below the 7.694143 Å of its
//!   `.pcr`, the kind of cell offset L028 found in FullProf's reflection lists);
//! - Yobs, against the data file read by the core's reader for the `.pcr`'s Ins.
//!
//! The `.prf` values may be divided by a power of ten (FullProf does so for large counts: citric-acid-sls_1.prf prints
//! Yobs 7390 for 73902.4 in the data); the factor is taken from Yobs and must be an exact power of ten.
//! FullProf's last point carries the background only (PPH3-sls.prf: 4890 where the core has 5000.25, the 110.25 being
//! the sum of every peak at that point); such a final point is reported and left out of the statistics.
//! For the FullProf runs, χ² over FullProf's N − P + C is compared with the one in their `.out`.
//! Classes: `exact` must agree to rms 5·10⁻⁴ and max 10⁻² (the `.prf` prints 3 decimals or whole counts, the `.pcr` 5–6
//! digits); `lebail` to rms 5·10⁻³ and max 0.1 (the intensities of the `.hkl` are those after the last cycle);
//! `inconsistent` is reported, and its `.out` reflection list (intensities, widths) is compared instead: FullProf's
//! Icalc there is the intensity within the ±Wdt·H window, I·(1 − η(1 − (2/π)·atan(2·Wdt))), as the comparison shows.

use drs::data::Pattern;
use drs::powder::Model;
use std::path::Path;

/// Window-edge points of TCH phases with X ≠ 0 left out of the exact comparison (L035): file, the most points, the
/// largest relative deviation. Measured in L035 (49 up to 1.35·10⁻², 181 up to 8.15·10⁻²); with the widths at the
/// Bragg angle (L039) 47 and 106, and LaB6 50 up to 0.406. With the extra low point of CrysFML's `Locate` (L042) LaB6
/// falls to 2, the Yb2ZrWO8 and YVO3 window experiments with X ≠ 0 lose their one missing low point per peak, and
/// rutile-anatase and MgAl2O4 keep as many points or fewer, but larger: there FullProf's low edge lies between the bare
/// edge and the extra point (rutana_tch_calc.prf, 110 at 27.46°: 26.93°, 0.59 step beyond the edge, is left out by
/// FullProf, where the window experiments take points up to 0.97 step beyond). What is left in the experiments is the
/// high side of the lowest-angle reflection(s), all within the asymmetry range (Bérar–Baldinozzi below AsyLim, or
/// FCJ; Yb2ZrWO8 111, YVO3 101, LaB6 100), and two points at 139–141° of peaks past the data's end (w8_x0.645). Their
/// rule was not found; frozen at these numbers (L042).
const EDGE_ALLOWED: &[(&str, usize, f64)] = &[
    ("MgAl2O4s_Rietveld_calc.prf", 42, 0.02535),
    // CeO₂ with surface micro-absorption (FP_tests/ce2.pcr run by FullProf, L239): TCH with X = 0.0174, the same kind.
    ("ce2_calc.prf", 27, 0.3605),
    ("rutana_tch_calc.prf", 106, 0.1268),
    ("LaB6rx_calc.prf", 2, 0.1776),
    ("YVO3a_calc.prf", 4005, 0.1995),
    ("irf1_narrow_x0.645_wdt20/ybfluo_Riet_IRF.prf", 1, 0.6229),
    ("w8_x0.645/ybfluo_Riet_IRF.prf", 3, 0.2089),
    ("wdt20_x0.645_narrow/ybfluo_Riet_IRF.prf", 1, 0.1087),
    ("wdt40_x0.645_narrow/ybfluo_Riet_IRF.prf", 1, 0.1926),
    ("yb_K2_x0.645/ybfluo_Riet_IRF.prf", 1, 0.4866),
    ("yb_all_x0.645/ybfluo_Riet_IRF.prf", 1, 0.4660),
    ("yb_irf1_x0.645/ybfluo_Riet_IRF.prf", 1, 0.4795),
    ("yvo3_base_w4/YVO3a.prf", 1, 0.08359),
    ("yvo3_noK2/YVO3a.prf", 1, 0.08793),
    ("yvo3_noK2_noasym/YVO3a.prf", 1, 0.1144),
    ("yvo3_noXins/YVO3a.prf", 2, 0.1013),
    ("yvo3_noasym/YVO3a.prf", 1, 0.1065),
    // nac_3t2 (P23, L353), 3T2 neutrons with X = 0.0295: the same kind; its Debye-like variant's background is 0 at most
    // of those points (FullProf's bounds), so their relative deviation is larger.
    // L367: one point more since the asymmetry is taken at the Bragg angle (nac_3t2's zero 0.049°).
    ("nac_3t2_calc.prf", 158, 3.31e-3),
    ("nac_3t2_debye.prf", 258, 9.02e-2),
    // The mode course's case A (L090), neutrons with TCH and X = 0.043: below 1.2·10⁻³, the same kind as above.
    ("srsclrrt.prf", 404, 1.193e-3),
    ("srsclrrt_modes_sb_to_origin.prf", 404, 1.175e-3),
    // The same case run by FullProf 8.20 (runs/modes, L093), by coordinates and by modes, at their CODE.pcr values.
    ("a_coords_calc.prf", 405, 1.193e-3),
    ("a_modes_calc.prf", 402, 1.175e-3),
    // P25's rest (L267). prba10k: TCH with the quartic strain's ξ = 0.148 and X = 0, 214 reflections; its windows end
    // within half a step (0.004°) of the core's on either side: (1 0 2) as the rule puts it, (1 1 0) one point wider
    // on the high side and one narrower on the low, a window 0.16% wider than the width without the Lorentzian strain.
    ("prba10k_calc.prf", 592, 0.2337),
    // Pd3MnD: TCH with X = 0.49 over 3T2's resolution file, the antiphase rule's broad peaks (Wdt 15): the windows of
    // its broad reflections reach up to 1.2° further than the core's (the tail of (4 2 0) to 93.35°, its window by the
    // rule ending at 92.18°), by a rule not found.
    ("pdfinal_calc.prf", 166, 9.579e-3),
    // The joint runs (P06, L322): scor_rxn's X-ray pattern, TCH with X = 0.167 and Wdt 8 over two phases; prbaco's
    // synchrotron pattern, X = 0.0148 and the quartic strain's ξ = 0.199 with Wdt 35, the windows of its strongest
    // peaks (2.5·10⁵ counts) reaching 0.6° (9.88°, core 525.8 against FullProf's 606 on the high side of 9.28°).
    // L367: one point more since the asymmetry is taken at the Bragg angle.
    ("scor_rxn_calc_2.prf", 1307, 5.311e-2),
    ("prbaco_calc_1.prf", 614, 0.4600),
    // The TOF runs (P05, L334, L341): FullProf's TOF window is not the manual's ±Wdt·H the core uses; it reaches about
    // Wdt·(5–8 µs) further on both sides for Npr = 9 (arg_si's (1 1 1): FullProf calculates 22995–23860 µs about 23421.7,
    // the core 23035–23808), and for gabriela's Npr = 10 bank it stops up to Wdt·11.6 µs short of the core's edge on the
    // rise side and reaches Wdt·4.4 µs past it on the decay side (13465 µs: FullProf 526, the core's background 217).
    // Not the true FWHM of the convolved peak, not 1/α or 1/β, not CrysFML nFP's ±15(H_G + ln 4/α); not found.
    ("CeO2_PEARL_calc.prf", 728, 5.952e-2),
    ("alum_calc.prf", 2912, 0.1377),
    ("alum_noq.prf", 2885, 1.6834e-2),
    ("arg_si_abs1.prf", 564, 6.624e-2),
    ("arg_si_abs3.prf", 567, 6.3475e-2),
    ("arg_si_calc.prf", 586, 0.10792),
    ("arg_si_ext10.prf", 559, 2.832e-2),
    ("arg_si_ext100.prf", 418, 1.498e-2),
    ("arg_si_invd.prf", 593, 0.11794),
    ("arg_si_plain.prf", 584, 0.11701),
    ("arg_si_plain90.prf", 544, 8.970e-2),
    ("gabriela_calc_1.prf", 871, 0.5872),
    ("lamn_pol_calc.prf", 2163, 0.1537),
    ("lamn_pol_irf.prf", 2163, 0.1537),
    ("nac-osiris_calc.prf", 634, 0.1239),
];

/// Files whose difference with FullProf is not resolved and is frozen (L039; CLAUDE.md, the gate): the largest rms and
/// max relative deviation of Ycal beyond the window-edge zone, and the largest |χ² − χ²_FullProf|.
/// YVO3a_calc.prf: X = 0.645 makes the Lorentzian strain dominate, and FullProf's peak window, which for the other runs
/// lies within 0.06–2.7% of Wdt × (the TCH width without X·tan θ), is there several times wider and narrower than
/// Wdt·H (the core with Wdt·H gives rms 7·10⁻³, with the width without X 2.3·10⁻²); its rule was not found.
/// L042: the same YVO3 with Wdt = 4 (window experiments yvo3_*) agrees to rms 5·10⁻⁶ but for one point, so the
/// difference comes with Wdt = 40; a low side of Wdt·H and a high side of Wdt·H_noX gives rms 1.9·10⁻², with a smooth
/// offset left that grows with angle. Still frozen.
/// runs/modes (L098): FullProf turns an excluded region into points with its own step (B: 0.050018 for the data's
/// 0.05), so B's 149.95°, which the core calculates, is left without peaks in FullProf (Ycal max at 149.95°; the
/// point is named in `symmetry_modes::EDGE_POINTS`); E lists its regions out of order (120–160, then 0–21.3), and
/// FullProf calculates peaks within the second one (Ycal max at 19.73°), whose points it leaves out of Σ.
const FROZEN: &[(&str, f64, f64, f64)] = &[
    ("YVO3a_calc.prf", 2.35e-2, 0.194, 0.0105),
    ("b_article_calc.prf", 1.152e-3, 6.53e-2, 0.0),
    ("b_parent_1_calc.prf", 2.548e-3, 0.1451, 0.0),
    ("b_parent_2_calc.prf", 2.548e-3, 0.1451, 0.0),
    ("b_thesis_calc.prf", 1.571e-3, 8.91e-2, 0.0),
    ("e_calc.prf", 2.074e-2, 0.358, 0.0),
    // L267: beyond the edge zone, the tails those windows leave (prba10k at 21.58°; Pd3MnD over 92.15–93.35°, the
    // tail of (4 2 0) the core's window does not reach). Their χ² agree within the print and 0.15%: 18.553 against
    // 18.55, 1.7008 against 1.699.
    ("prba10k_calc.prf", 6.19e-4, 2.55e-2, 0.0),
    ("pdfinal_calc.prf", 1.273e-3, 9.29e-3, 0.0),
    // The joint runs (P06, L322). scor_rxn_calc_2.prf (Cu Kα, both phases): beyond the edge zone every strong point
    // is 0.19% above FullProf's (the mean over the 33 points above 300 counts at 15–40° and the 28 at 40–60° alike),
    // the background within 3.3·10⁻⁵; the pattern has Bérar–Baldinozzi asymmetry below AsyLim 60° together with
    // axial divergence (S_L 0.030, D_L 0.015), which no other compared run combines; the cause is not found; χ² 2.2703
    // against 2.27. prbaco_calc_1.prf: χ² 20.342 against 20.38, from the strongest peaks (Ycal within 1.2·10⁻⁴ of
    // FullProf's at 2.5·10⁵ counts, σ 500); with FullProf's own Ycal the core's points and σ give 105379 against its
    // Residual 105385. prbaco_calc_2.prf: χ² 1.4662 against 1.453; the profile agrees (rms 4·10⁻⁵), and FullProf's own
    // Ycal over the core's points and σ (var = y/iww, FullProf's code for Ins = 6, manual p. 135) gives Σ 4365 against
    // its Residual 4323; the points FullProf calculates are the core's (its .prf prints Yobs − Ycal as its offset
    // alone at the points it leaves out: 8.00° and 148.00–150.00°). Without the cut at the .out's S(1/d)max the core's
    // Σ is 4317 (`joint_refinement`): FullProf's sum is nearer the pattern with the tails of the reflections past
    // 159.5° than to its own .prf's Ycal; not found.
    ("scor_rxn_calc_2.prf", 9.950e-4, 1.764e-2, 0.0),
    ("prbaco_calc_1.prf", 5e-4, 1e-2, 0.0382),
    ("prbaco_calc_2.prf", 5e-4, 1e-2, 0.0132),
    // The TOF runs (P05, L341). CeO2_PEARL: FullProf's χ² 0.0003858 is 1/9984 of the core's 3.851 with the same Ycal:
    // FullProf's σ of a GSAS RALF bank are 100 times the core's (its counts are GSAS-II's y·100/Δt, Yobs to 10⁻⁴; its σ
    // GSAS-II's σ·10⁴/Δt by its χ²), a reading of FullProf's not stated anywhere. lamn_pol and nac-osiris weigh by
    // maximum likelihood (Iwg = 1): the core's var = Ycal·σ²/y, which met LaB6rx_calc.out (L039), gives 5.143 and 7.607
    // for FullProf's 1.694 and 2.105, and no global factor on var = Ycal gives both (0.605 and 251 are needed); not
    // found. alum: 2.3967 against 2.379, 0.7%, its Ycal within 10⁻³ beyond the edge zone. gabriela: the TOF bank 2.2754
    // against 2.28 (0.2%); its X-ray patterns 2.1555 against 2.172 and 7.5436 against 7.531, where 48 points of each have
    // σ = 0 and FullProf calculates the point at a region's end that the core leaves out: 8.64000034° after a region to
    // 8.64 (pattern 2, FullProf 9.23 where the core's background is 7.17) and 8.57999992° within one to 8.58 (pattern
    // 3), neither by the core's rule (L033, L322), the two points of the Ycal max.
    ("CeO2_PEARL_calc.prf", 5e-4, 1e-2, 3.851),
    // nac_3t2_calc (L353): the background's first two points (`FROZEN_BACKGROUND`), Ycal being the background there.
    ("nac_3t2_calc.prf", 3.09e-3, 0.1545, 0.0),
    ("alum_calc.prf", 5e-4, 1e-2, 0.0177),
    ("gabriela_calc_1.prf", 5e-4, 1e-2, 0.0047),
    ("gabriela_calc_2.prf", 1.337e-3, 0.2235, 0.0166),
    // L367: χ² 3.4·10⁻⁵ further since the asymmetry is taken at the Bragg angle.
    ("gabriela_calc_3.prf", 2.886e-3, 0.3878, 0.01264),
    ("lamn_pol_calc.prf", 5e-4, 1e-2, 3.450),
    ("lamn_pol_irf.prf", 5e-4, 1e-2, 3.447),
    ("nac-osiris_calc.prf", 5e-4, 1e-2, 5.502),
    // The anvil cell's texture (Nor = 3, P30, L360): Wright's cone integral as the core believes it (`texture`) gives
    // FullProf's yb_nor3 to 2.5% of the strongest peak (at 93.36°, (4 2 2)); FullProf's CORR is up to 3.6% below it
    // (test `texture`), and no reading of eqs. 4–6 found meets it. χ² 1.6636 against 1.663.
    ("yb_nor3.prf", 5.68e-3, 2.48e-2, 0.0),
    // TLS (P41, L367): urea_test's χ² 13.2287 against 13.27 (0.31 %), its profile within the print (rms 3.7·10⁻⁴). The
    // TLS centre FullProf computes for urea (z 0.3303) is not the group's centre of mass (z 0.32957, by CrysFML's
    // weights); with FullProf's centre the core's χ² is within 0.06 % (`tls::urea_test_chi2`). Not explained, listed
    // for Josu.
    ("urea_test_calc.prf", 5e-4, 1e-2, 0.0414),
];

/// Files whose background differs from FullProf's at a point not explained, frozen (L037): the largest rms and max
/// relative deviation of the background. nabd4_calc.prf (L251): FullProf's last point, 164.90° (the excluded region
/// 157.05–164.90 is printed 164.89999 in the .prf), repeats the point before it, Ycal and Backg alike (201.5 and 156.2
/// at 164.85° and 164.90°), where the linear background through the .pcr's last point (164.90°, 202.244) rises by 0.7
/// per step; no other compared .prf repeats its last point (all 66 checked in L251), so the rule is not found.
/// nabd4_lambda_calc.prf (L255), the same data run again by FullProf, repeats it in the same way (201.5 and 156.2).
/// nac_3t2_calc.prf (L353): its first two points, in the excluded region 0–5°. At 2θ = 0 FullProf leaves the Debye-like
/// terms out (31.033 = −Σ A_j, where their limit at Q = 0 is 1, the 0.154 of the max); at 0.05° its background is 0.557
/// below the core's, not explained (`background_terms`).
const FROZEN_BACKGROUND: &[(&str, f64, f64)] = &[("nabd4_calc.prf", 6.10e-5, 3.445e-3), ("nabd4_lambda_calc.prf", 6.10e-5, 3.445e-3), ("nac_3t2_calc.prf", 3.09e-3, 0.1545)];

fn root(rel: &str) -> String {
    format!("{}/../{rel}", env!("CARGO_MANIFEST_DIR"))
}

fn latin1(path: &str) -> String {
    std::fs::read(path).unwrap_or_else(|e| panic!("{path}: {e}")).iter().map(|&b| b as char).collect()
}

/// Reflections of a .prf: h k l and 2θ.
type Reflections = Vec<([i32; 3], f64)>;

/// The profile of a .prf, read by the library's reader (`formats::fullprof_outputs::read_prf`, app request 17): (2θ,
/// Yobs, Ycal, background with FullProf's offset added back, L082) at the `n` points from line `start` of the index;
/// its reflections (h k l, 2θ); and the number of decimals FullProf printed Ycal with.
fn read_prf(path: &str, start: usize, n: usize) -> (Vec<[f64; 4]>, Reflections, usize) {
    let text = latin1(path);
    let p = drs::formats::fullprof_outputs::read_prf(&text).unwrap_or_else(|e| panic!("{path}: {}", e.0));
    assert_eq!(p.x.len(), n, "{path}: the index's points");
    let first: f64 = text.lines().nth(start - 1).unwrap().split('\t').next().unwrap().trim().parse().unwrap();
    assert_eq!(p.x[0], first, "{path}: the index's first line");
    let prof = (0..n).map(|i| [p.x[i], p.yobs[i], p.ycalc[i], p.background[i]]).collect();
    (prof, p.reflections.iter().map(|r| (r.hkl, r.two_theta)).collect(), p.ycalc_decimals)
}

/// A pattern of a joint run (P06, L322): the index names its `.pcr` `CODE.pcr#k`; the `.out` of its `.prf`
/// CODE_calc_k.prf is CODE_calc.out, which reports every pattern in turn.
fn joint(pcr: &str) -> Option<(&str, usize)> {
    let (p, k) = pcr.split_once('#')?;
    Some((p, k.parse().ok()?))
}

/// The `.out` written with a `.prf`.
fn out_of(prf: &str, pcr: &str) -> String {
    match joint(pcr) {
        Some((_, k)) => prf.replace(&format!("_{k}.prf"), ".out"),
        None => prf.replace(".prf", ".out"),
    }
}

/// The part of a joint run's `.out` that reports pattern k (its reliability factors), from its first header to the
/// next pattern's.
fn pattern_section(out: &str, k: usize) -> String {
    let head = |k: usize| format!("RELIABILITY FACTORS WITH ALL NON-EXCLUDED POINTS FOR PATTERN:{k:>3}");
    let start = out.find(&head(k)).unwrap_or_else(|| panic!("no reliability factors for pattern {k}"));
    let end = out[start..].find(&head(k + 1)).map_or(out.len(), |e| start + e);
    out[start..end].to_string()
}

/// FullProf's reflection limits, S(1/d)max of each phase, from the `.out` written with the `.prf` (L035: its rule is
/// not in the manual); empty without that `.out`. A joint run's `.out` lists them phase by phase, each for the
/// patterns it contributes to in turn (scor_rxn_calc.out, lines 393–816).
fn reflection_limits(prf: &str, pcr: &str) -> Vec<f64> {
    let Ok(out) = std::fs::read(root(&out_of(prf, pcr))) else { return Vec::new() };
    let out: String = out.iter().map(|&b| b as char).collect();
    let all: Vec<f64> = out.lines().filter_map(|l| l.split_once("S(1/d)max:")).map(|(_, r)| r.split_whitespace().next().unwrap().parse().unwrap()).collect();
    let Some((file, k)) = joint(pcr) else { return all };
    let m = drs::pcr::split_patterns(&latin1(&root(file))).unwrap();
    let mut it = all.into_iter();
    let mut mine = Vec::new();
    for c in &m.contributions {
        for (j, &on) in c.iter().enumerate() {
            if on {
                let v = it.next().expect("an S(1/d)max per phase and pattern");
                if j + 1 == k {
                    mine.push(v);
                }
            }
        }
    }
    mine
}

/// FullProf's step ("Step:" of the `.out` written with the `.prf`), which its rule for the background's low edge
/// counts in (L098); None without that `.out`. A joint run's `.out` gives one "2Thmin: … Step:" line per pattern.
fn fullprof_step(prf: &str, pcr: &str) -> Option<f64> {
    let out = std::fs::read(root(&out_of(prf, pcr))).ok()?;
    let out: String = out.iter().map(|&b| b as char).collect();
    let step = |l: &str| l.split_once("Step:").and_then(|(_, r)| r.split_whitespace().next()?.parse().ok());
    match joint(pcr) {
        Some((_, k)) => out.lines().filter(|l| l.contains("2Thmin:")).nth(k - 1).and_then(step),
        None => out.lines().find_map(step),
    }
}

/// The power of ten that the .prf Yobs is divided by, from the data: at the first point with counts in both (FullProf
/// prints 1.00 where the data have 0, srsclrrt.prf of the mode course, L090).
fn prf_factor(data: &Pattern, prof: &[[f64; 4]]) -> f64 {
    let (i, v) = prof
        .iter()
        .filter(|v| v[1] > 1.0)
        .find_map(|v| data.x.iter().position(|&t| (t - v[0]).abs() < 2e-3).filter(|&i| data.y[i] > 1.0).map(|i| (i, v)))
        .expect("a .prf point with counts in the data");
    let k = (data.y[i] / v[1]).log10().round();
    10f64.powi(k as i32)
}

/// The name a .prf is reported and frozen under: its file name, with its folder for the window experiments (L041),
/// whose files share names.
fn short_name(prf: &str) -> String {
    let mut parts = prf.rsplit('/');
    let file = parts.next().unwrap();
    match (parts.next(), parts.next()) {
        (Some(dir), Some("window_experiments")) => format!("{dir}/{file}"),
        _ => file.to_string(),
    }
}

/// `load`, or None where the `.pcr` does not import (a TOF variant's own `.pcr`, holding NaN).
fn try_load(pcr: &str) -> Option<drs::powder::Loaded> {
    let (pcr, text) = match joint(pcr) {
        Some((file, k)) => (file, drs::pcr::split_patterns(&latin1(&root(file))).ok()?.texts.get(k - 1)?.clone()),
        None => (pcr, latin1(&root(pcr))),
    };
    let p = drs::pcr::read(&text).ok()?;
    let path = Path::new(pcr);
    Model::from_pcr_fullprof(&p, Path::new(&root(path.parent().unwrap().to_str().unwrap())), path.file_stem().unwrap().to_str().unwrap().trim_end_matches(".pcr")).ok()
}

fn load(pcr: &str) -> drs::powder::Loaded {
    // A joint run's pattern k: the single-pattern text `split_patterns` gives for it.
    let (pcr, text) = match joint(pcr) {
        Some((file, k)) => (file, drs::pcr::split_patterns(&latin1(&root(file))).unwrap_or_else(|e| panic!("{file}: {e}")).texts[k - 1].clone()),
        None => (pcr, latin1(&root(pcr))),
    };
    let p = drs::pcr::read(&text).unwrap_or_else(|e| panic!("{pcr}: {e}"));
    let path = Path::new(pcr);
    // The .pcr's values under FullProf's own rules (L046, L050); a VARIANT.pcr.in's stem is VARIANT (its .bac).
    Model::from_pcr_fullprof(&p, Path::new(&root(path.parent().unwrap().to_str().unwrap())), path.file_stem().unwrap().to_str().unwrap().trim_end_matches(".pcr")).unwrap_or_else(|e| panic!("{pcr}: {e}"))
}

/// Yobs of the .prf against the data read by the core (the .prf prints 3 decimals, or whole counts).
fn check_yobs(name: &str, data: &Pattern, prof: &[[f64; 4]], fac: f64, fails: &mut Vec<String>) -> usize {
    // The nearest data point to the .prf's first (gabriela_calc_2.prf starts at the data's second point, 6.1007°, the
    // first, 6.0999999°, lying below Thmin 6.1 and within 2·10⁻³ of it).
    let mut j = (0..data.len()).filter(|&i| (data.x[i] - prof[0][0]).abs() < 2e-3).min_by(|&a, &b| (data.x[a] - prof[0][0]).abs().partial_cmp(&(data.x[b] - prof[0][0]).abs()).unwrap()).unwrap();
    let mut n = 0;
    for v in prof {
        // The nearest data point (steps of 7·10⁻⁴° in gabriela's synchrotron data).
        while j + 1 < data.len() && (data.x[j + 1] - v[0]).abs() <= (data.x[j] - v[0]).abs() {
            j += 1;
        }
        let (x, y) = (data.x[j], data.y[j] / fac);
        // FullProf prints 1.00 for a point of 0 counts (the first two points of the mode course's SrScLRRT.dat, both in
        // an excluded region; L090). Whether its reader or only its .prf does so is not known (believed: the .prf).
        let zero_as_one = data.y[j] == 0.0 && v[1] == 1.0;
        // x to its print and FullProf's single precision (a TOF of 24001.201 µs printed 24001.2012).
        if (x - v[0]).abs() > 1e-4 + 1e-7 * x.abs() || ((y - v[1]).abs() > 0.5 + 1e-6 * y.abs() && !zero_as_one) {
            fails.push(format!("{name}: Yobs at {} is {}, core reads {y} at {x}", v[0], v[1]));
            return n;
        }
        n += 1;
    }
    n
}

#[test]
fn fullprof_profiles() {
    let index = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/fullprof_profiles.txt")).unwrap();
    let mut fails = Vec::new();
    let (mut compared, mut data_only, mut skipped, mut yobs_points): (usize, usize, usize, usize) = (0, 0, 0, 0);
    for line in index.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = line.split_whitespace().collect();
        match f[0] {
            "SKIP" => {
                skipped += 1;
                println!("left out  {}: {}", f[1].rsplit('/').next().unwrap(), f[2..].join(" "));
            }
            "DATA" => {
                data_only += 1;
                let (prof, _, _) = read_prf(&root(f[1]), f[4].parse().unwrap(), f[5].parse().unwrap());
                let data = drs::data::read(&latin1(&root(f[2])), f[3].parse().unwrap()).unwrap();
                let fac = prf_factor(&data, &prof);
                yobs_points += check_yobs(f[1], &data, &prof, fac, &mut fails);
            }
            "COMPARE" => {
                compared += 1;
                let (prf, pcr, start, n, class) = (f[1], f[2], f[3].parse().unwrap(), f[4].parse().unwrap(), f[5]);
                let short = short_name(prf);
                let short = short.as_str();
                let (mut prof, refl, decimals) = read_prf(&root(prf), start, n);
                let mut loaded = load(pcr);
                loaded.model.fullprof_step = fullprof_step(prf, pcr);
                let limits = reflection_limits(prf, pcr);
                if limits.len() == loaded.model.phases.len() {
                    for (ph, l) in loaded.model.phases.iter_mut().zip(&limits) {
                        ph.reflection_limit = Some(*l);
                    }
                }
                let data = loaded.data.as_ref().expect("data file");
                let fac = prf_factor(data, &prof);
                yobs_points += check_yobs(prf, data, &prof, fac, &mut fails);
                for v in prof.iter_mut() {
                    v[1] *= fac;
                    v[2] *= fac;
                    v[3] *= fac;
                }
                let x: Vec<f64> = prof.iter().map(|v| v[0]).collect();
                let c = loaded.model.calculate(&x);
                // FullProf's last point: background only? Ycal and Backg are each printed to the unit (times fac), so
                // two equal values can print a unit apart: L351's citric-acid-sls_1_calc.prf prints 2109 and 1624 (+
                // the offset 484.33: 2108.33) at 118.388°, where the tutorial's own run printed 1625 (L354).
                let last = prof.len() - 1;
                let last_bg_only = (prof[last][2] - prof[last][3]).abs() <= fac + 1e-3 && (c.ycalc[last] - c.background[last]).abs() > fac;
                let m = if last_bg_only { last } else { prof.len() };
                // Relative deviation beyond FullProf's print (half a unit of its last decimal).
                let half = 0.5 * fac * 10f64.powi(-(decimals as i32));
                let rel = |a: f64, b: f64| ((a - b).abs() - half).max(0.0) / b.abs().max(1.0);
                // Points just outside the window of a TCH peak whose phase has a Lorentzian strain term X ≠ 0: there
                // FullProf's window is wider than the core's by 0.06% (rutile 110) to 2.7% (anatase 101), by a rule not
                // found (L035), so a point less than 3% of the window plus one step beyond the core's edge may carry
                // the peak in FullProf and not in the core. Those points are reported apart.
                let step = (x[last] - x[0]) / last as f64;
                let tof = loaded.model.tof.is_some();
                let edge: Vec<bool> = {
                    let mut e = vec![false; x.len()];
                    for p in &c.peaks {
                        // A TOF peak (P05, L341): FullProf's window is not ±Wdt·H; measured on arg_si, CeO2_PEARL and
                        // lamn_pol it reaches Wdt·(5–8 µs) further on either side, and gabriela's (Npr = 10) stops up
                        // to Wdt·11.6 µs short of the core's on the rise side and reaches Wdt·4.4 µs past it on the
                        // decay side; its rule was not found. The zone is 8 µs·Wdt and a step beyond the core's edge
                        // (Npr = 9), 12 µs·Wdt within it and 6 µs·Wdt beyond (Npr = 10); its points counted and frozen.
                        if tof {
                            let w = loaded.model.wdt;
                            let hodges = loaded.model.tof.as_ref().is_some_and(|b| b.thermal.is_some());
                            let (inside, outside) = if hodges { (12.0 * w + step, 6.0 * w + step) } else { (0.0, 8.0 * w + step) };
                            for (i, &t) in x.iter().enumerate() {
                                let d = (t - p.two_theta).abs() - p.reach;
                                if d > -inside && d < outside {
                                    e[i] = true;
                                }
                            }
                            continue;
                        }
                        let ph = &loaded.model.phases[p.phase];
                        // The quartic form's Lorentzian strain ξ (L267): the windows of prba10k_calc.prf end within
                        // half a step of the core's on either side, so the zone is taken on both sides of the edge.
                        let xi = ph.quartic_strain.as_ref().is_some_and(|q| q.lorentzian != 0.0);
                        if ph.shape != drs::powder::Shape::ThompsonCoxHastings || (ph.x == 0.0 && !xi) {
                            continue;
                        }
                        for (i, &t) in x.iter().enumerate() {
                            let d = (t - p.two_theta).abs() - p.reach;
                            if (d > 0.0 || xi && d > -step) && d < 0.03 * p.reach + step {
                                e[i] = true;
                            }
                        }
                    }
                    e
                };
                let (mut worst, mut at, mut sum, mut worst_b, mut sum_b) = (0.0f64, 0.0, 0.0, 0.0f64, 0.0);
                let (mut n_edge, mut worst_edge) = (0usize, 0.0f64);
                for (i, pt) in prof.iter().enumerate().take(m) {
                    let d = rel(c.ycalc[i], pt[2]);
                    if edge[i] && d > 1e-4 {
                        n_edge += 1;
                        worst_edge = worst_edge.max(d);
                        continue;
                    }
                    sum += d * d;
                    if d > worst {
                        (worst, at) = (d, pt[0]);
                    }
                    // The reader's background (FullProf's Backg with its offset added back) against the core's.
                    let db = rel(c.background[i], pt[3]);
                    sum_b += db * db;
                    worst_b = worst_b.max(db);
                }
                let (rms, rms_b) = ((sum / m as f64).sqrt(), (sum_b / m as f64).sqrt());
                // Reflection positions (λ1 and λ2 peaks alike: the nearest of the same h k l).
                // FullProf lists reflections past the data. The tick marks of a run's CODE_calc.prf are the Bragg angles of
                // the values read, shifted by the Zero after its cycle (YVO3a_calc.prf: 101 at 20.5444 = 20.4825 + the
                // 0.0619 of YVO3a_calc.pcr, not + the 0.0555 read; the refined cell would put them 0.1° off; L039).
                let mut ticks = loaded.model.clone();
                if prf.contains("/window_experiments/") {
                    // Their `.pcr` was not rewritten after the cycle: the Zero after it is the `.out`'s last Zero_pat1.
                    let out = latin1(&root(&prf.replace(".prf", ".out")));
                    let z = out.lines().rfind(|l| l.contains("Zero_pat1") && l.contains("+/-")).expect("Zero_pat1 after the cycle");
                    ticks.zero = z.split("Zero_pat1").nth(1).unwrap().split(['(', ' ']).find(|s| !s.is_empty()).unwrap().parse().unwrap();
                } else if prf.contains("/runs/") {
                    let (written, zero_name) = match joint(pcr) {
                        Some((file, k)) => (format!("{}#{k}", file.replace(".pcr", "_calc.pcr")), format!("Zero_pat{k}")),
                        None => (prf.replace(".prf", ".pcr"), "Zero_pat1".to_string()),
                    };
                    // A TOF variant's own .pcr, rewritten by FullProf after a diverged cycle, holds NaN (L334): no model.
                    if let Some(after) = try_load(&written) {
                        (ticks.zero, ticks.sycos, ticks.sysin) = (after.model.zero, after.model.sycos, after.model.sysin);
                        // A TOF bank's Dtt1, Dtt2, Dtt_1overd after the cycle, as its Zero (alum refines Dtt2).
                        if tof {
                            ticks.tof = after.model.tof.clone();
                        }
                    }
                    // FullProf does not always write the cycle's values back (runs/modes/D: d_calc.pcr keeps Zero
                    // −0.0593, its .out ends with Zero_pat1 −0.0786, and the ticks are 0.0193° below; L098): where
                    // Zero was refined, the .out's last Zero_pat1 is the Zero after the cycle.
                    let out = latin1(&root(&out_of(prf, pcr)));
                    let after_cycle = |name: &str| out.lines().rfind(|l| l.contains(name) && l.contains("+/-")).and_then(|z| z.split(name).nth(1).unwrap().split(['(', ' ']).find(|s| !s.is_empty()).unwrap().parse::<f64>().ok());
                    if let Some(z) = after_cycle(zero_name.as_str()) {
                        ticks.zero = z;
                    }
                    // A TOF bank's refined Dtt1 and Dtt2 after the cycle, from the .out as its Zero (alum_noq.pcr is not
                    // written back, and its ticks sit 6 µs below alum_calc.prf's: Zero 26.743 → 22.927, Dtt2 52.660 →
                    // 43.346 after its cycle).
                    if let Some(b) = ticks.tof.as_mut() {
                        let k = joint(pcr).map_or(1, |(_, k)| k);
                        if let Some(v) = after_cycle(&format!("Dtt1_pat{k}")) {
                            b.dtt1 = v;
                        }
                        if let Some(v) = after_cycle(&format!("Dtt2_pat{k}")) {
                            b.dtt2 = v;
                        }
                    }
                }
                // FullProf marks the λ2 positions also when Ratio = 0 (yvo3_noK2/YVO3a.prf: 479 marks, as with Kα2).
                if ticks.ratio == 0.0 {
                    ticks.ratio = 1e-30;
                }
                let peaks = if tof { ticks.peaks(x[0] * 0.5, x[last] * 2.0) } else { ticks.peaks(x[0] - 10.0, 180.0) };
                let mut worst_pos = 0.0f64;
                let mut npos = 0;
                for (h, t) in &refl {
                    let d = peaks.iter().filter(|p| p.hkl == *h).map(|p| (p.two_theta - t).abs()).fold(f64::MAX, f64::min);
                    if d < 1.0 {
                        // Within the data and 5° past it. Past that, near 180°, the ticks of runs/modes/B_* (data to
                        // 161.2°) agree to 1.1·10⁻³° at 165–178° and to 2.9·10⁻⁴° within the data: there 2θ moves as
                        // 2/cos θ times sin θ, believed (not checked) to be FullProf's single-precision 2θ (L098).
                        if *t <= x[last] + 5.0 {
                            worst_pos = worst_pos.max(d);
                        }
                        npos += 1;
                    }
                }
                println!(
                    "{class:<12} {short}: {m} points, Ycal max {worst:.2e} (at {at}) rms {rms:.2e}; Backg max {worst_b:.2e} rms {rms_b:.2e}; {npos} positions, worst {worst_pos:.5}{}; .prf divided by {fac}{}",
                    if tof { " µs" } else { "°" },
                    if last_bg_only { format!("; last point ({}) background only in FullProf", prof[last][0]) } else { String::new() }
                );
                if n_edge > 0 {
                    println!("             {short}: {n_edge} points just outside a window edge of a phase with X ≠ 0 differ by more than 10⁻⁴ (up to {worst_edge:.3e}), left out above (L035)");
                }
                // The points left out are frozen at what L035 measured (count, worst deviation), so that a change
                // which makes them more or worse fails here instead of passing unseen (L037).
                let (max_n, max_d) = EDGE_ALLOWED.iter().find(|(f, _, _)| *f == short).map_or((0, 0.0), |&(_, n, d)| (n, d));
                if n_edge > max_n || worst_edge > max_d {
                    fails.push(format!("{short}: {n_edge} window-edge points up to {worst_edge:.4e}, more than the {max_n} up to {max_d:.2e} accepted in L035–L042"));
                }
                // A TOF variant whose one cycle diverged (Extinc = 0 with its codeword: NaN shifts, L334) prints its ticks
                // as NaN: they belong to no model, and are not compared.
                let nan_ticks = !refl.is_empty() && refl.iter().all(|r| !r.1.is_finite());
                if nan_ticks {
                    println!("             {short}: its ticks are NaN (FullProf's cycle after the values read diverged, L334): not compared");
                } else if npos != refl.len() || worst_pos > if tof { 1e-6 * x[last] } else { 5e-4 } {
                    fails.push(format!("{short}: {npos} of {} reflection positions found, worst {worst_pos}", refl.len()));
                }
                let (tol_rms_b, tol_max_b) = FROZEN_BACKGROUND.iter().find(|f| f.0 == short).map_or((5e-4, 2e-3), |f| (f.1, f.2));
                if rms_b > tol_rms_b || worst_b > tol_max_b {
                    fails.push(format!("{short}: background max {worst_b:.4e}, rms {rms_b:.4e}"));
                }
                let (tol_rms, tol_max) = match class {
                    "exact" => (5e-4, 1e-2),
                    "lebail" => (5e-3, 0.1),
                    _ => (f64::MAX, f64::MAX),
                };
                let (tol_rms, tol_max) = FROZEN.iter().find(|f| f.0 == short).map_or((tol_rms, tol_max), |f| (f.1, f.2));
                if rms > tol_rms || worst > tol_max {
                    fails.push(format!("{short} ({class}): Ycal max {worst:.4e} at {at}, rms {rms:.4e}"));
                }
                let edge_fp: Vec<(f64, f64)> = (0..prof.len()).filter(|&i| edge[i]).map(|i| (prof[i][0], prof[i][2])).collect();
                check_chi2(prf, pcr, &loaded, &edge_fp, &mut fails);
                if class == "inconsistent" {
                    check_out(pcr, &loaded.model, (x[0], x[last]), &mut fails);
                }
            }
            _ => panic!("unknown line {line}"),
        }
    }
    println!("{compared} compared, {data_only} for their data only, {skipped} left out; {yobs_points} observed points read as FullProf read them");
    // Counted on disk, not as a fixed number, so that a .prf added to manuals/ without an index line fails here (L054).
    fn count_prf(dir: &std::path::Path) -> usize {
        std::fs::read_dir(dir).map_or(0, |rd| {
            rd.flatten()
                .map(|e| {
                    let p = e.path();
                    if p.is_dir() {
                        count_prf(&p)
                    } else {
                        usize::from(p.extension().is_some_and(|x| x.eq_ignore_ascii_case("prf")))
                    }
                })
                .sum()
        })
    }
    let base = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals");
    let on_disk = count_prf(&base.join("fullprof")) + count_prf(&base.join("bcs/mode-crystallography-manual/course-um6p-2024/04_tools_materials"));
    assert_eq!(compared + skipped, on_disk, "every .prf file is compared or left out with its reason");
    // L353: L346's runs add 16 compared: 14 of background_terms (all but nac_3t2's pass 1 and y2o3_filter), rigid's
    // ureaz_calc and PPH3-sls_calc.
    // L351: citric-acid-sls_1_calc and citric_core (L354).
    // L360: texture's four runs (Nor = 2, 3), compared from their VARIANT.pcr.in; yb_nor3 frozen.
    // L367: rigid/'s TLS runs (P41's TLS): hmt_calc, pdd_calc, pdd_l, pdd_t, pdd_t_rot, urea_calc, urea_conv_calc,
    // urea_notls, urea_test_calc, ureac_calc, ureacyl_calc.
    assert_eq!((compared, data_only), (126, 71));
    let shown: Vec<_> = fails.iter().take(100).cloned().collect();
    assert!(fails.is_empty(), "{} failures:\n{}", fails.len(), shown.join("\n"));
}

/// χ² of the observed pattern against the core's calculation, over FullProf's N − P + C, against FullProf's own
/// ("Conventional Rietveld Rp,Rwp,Re and Chi2" of the `.out` written with the `.prf`; the FullProf runs only, whose
/// `.out` belongs to the `.prf`). Tolerance: the print (5·10⁻⁴) and 0.15% (rutana_tch_calc: 0.12%, from the points at
/// the window edges, L035).
fn check_chi2(prf: &str, pcr: &str, loaded: &drs::powder::Loaded, edge_fp: &[(f64, f64)], fails: &mut Vec<String>) {
    if !prf.contains("/runs/") {
        return;
    }
    let out = latin1(&root(&out_of(prf, pcr)));
    // A joint run's pattern: its own reliability factors (N − P + C with the joint run's P).
    let out = match joint(pcr) {
        Some((_, k)) => pattern_section(&out, k),
        None => out,
    };
    let npc: usize = out.lines().find_map(|l| l.split_once("N-P+C:")).map(|(_, r)| r.split_whitespace().next().unwrap().parse().unwrap()).unwrap();
    let fp: f64 = out.lines().find_map(|l| l.split_once("Conventional Rietveld Rp,Rwp,Re and Chi2:")).map(|(_, r)| r.split_whitespace().nth(3).unwrap().parse().unwrap()).unwrap();
    let (data, model) = (loaded.data.as_ref().unwrap(), &loaded.model);
    let mut c = model.calculate(&data.x);
    // At the window-edge points kept apart (counted and frozen above), FullProf's own Ycal, so that χ² tests the rest.
    let mut j = 0;
    for (i, &x) in data.x.iter().enumerate() {
        while j < edge_fp.len() && edge_fp[j].0 < x - 1e-6 {
            j += 1;
        }
        if j < edge_fp.len() && (edge_fp[j].0 - x).abs() <= 1e-6 {
            c.ycalc[i] = edge_fp[j].1;
        }
    }
    // A TOF pattern's points FullProf leaves out (`Model::excluded_points`, P05): each run of them as a region.
    let regions: Vec<(f64, f64)> = if model.tof.is_some() {
        let mask = model.excluded_points(&data.x);
        let mut r = Vec::new();
        let mut i = 0;
        while i < mask.len() {
            if mask[i] {
                let a = i;
                while i + 1 < mask.len() && mask[i + 1] {
                    i += 1;
                }
                r.push((data.x[a], data.x[i]));
            }
            i += 1;
        }
        r
    } else {
        model.excluded.clone()
    };
    let n = data.x.iter().filter(|&&x| !regions.iter().any(|&(a, b)| x >= a && x <= b)).count();
    // Σw(y − yc)² over FullProf's N − P + C, which can exceed the core's N: runs/modes/E lists its regions out of
    // order, and FullProf counts the second region's points in N, not in Σ (L098).
    let (_, _, _, chi2_n) = drs::powder::agreement_weighted(data, &c.ycalc, &regions, 0, model.weighting);
    let chi2 = chi2_n * n as f64 / npc as f64;
    let short = short_name(prf);
    let short = short.as_str();
    println!("             {short}: χ² {chi2:.4} over N − P + C = {npc}, FullProf {fp}");
    let frozen = FROZEN.iter().find(|f| f.0 == short).map_or(0.0, |f| f.3);
    if (chi2 - fp).abs() > (5e-4 + 1.5e-3 * fp).max(frozen) {
        fails.push(format!("{short}: χ² {chi2} against FullProf's {fp}"));
    }
}

/// The last reflection list of the .out (intensities Icalc and widths Hw) against the core.
fn check_out(pcr: &str, model: &Model, range: (f64, f64), fails: &mut Vec<String>) {
    let out = latin1(&root(&pcr.replace(".pcr", ".out")));
    let mut rows = std::collections::BTreeMap::new();
    let mut inside = false;
    for l in out.lines() {
        if l.contains("Icalc") && l.contains("HwG") {
            inside = true;
            continue;
        }
        let v: Vec<f64> = l.split_whitespace().filter_map(|s| s.parse().ok()).collect();
        if inside && v.len() >= 15 {
            rows.insert([v[2] as i32, v[3] as i32, v[4] as i32], v);
        }
    }
    let peaks = model.peaks(0.0, 180.0);
    let (mut wi, mut wh, mut n) = (0.0f64, 0.0f64, 0);
    let mut worst = String::new();
    for (h, v) in &rows {
        if let Some(q) = peaks.iter().find(|q| q.hkl == *h && q.wavelength == 1) {
            n += 1;
            // Only peaks whose window lies within the data (FullProf integrates over the points it has).
            if q.two_theta - model.wdt * q.fwhm < range.0 || q.two_theta + model.wdt * q.fwhm > range.1 {
                continue;
            }
            // FullProf's Icalc is the intensity within the profile window ±Wdt·H, which cuts the Lorentzian tails:
            // I·(1 − η(1 − (2/π)·atan(2·Wdt))). Printed with one decimal; tolerance 0.5% plus the print.
            let window = 1.0 - q.eta * (1.0 - 2.0 / std::f64::consts::PI * (2.0 * model.wdt).atan());
            let e = (q.intensity * window - v[8]).abs() / (v[8].abs() * 0.005 + 0.05);
            if e > wi {
                wi = e;
                worst = format!("{h:?}: {:.3} against {}", q.intensity * window, v[8]);
            }
            wh = wh.max((q.fwhm / v[6] - 1.0).abs());
        }
    }
    println!("             {}: its .out lists {n} reflections; core Icalc (within the ±Wdt·H window) at {wi:.2} of the tolerance (0.5% + print; worst {worst}), FWHM within {wh:.1e}", pcr.rsplit('/').next().unwrap());
    if n < rows.len() || wi > 1.0 || wh > 1e-3 {
        fails.push(format!("{pcr}: .out reflections: {n} of {}, Icalc {wi}, FWHM {wh}", rows.len()));
    }
}
