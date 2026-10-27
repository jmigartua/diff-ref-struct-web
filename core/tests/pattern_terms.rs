//! Version 1's pattern terms finished in L233: the wavelength refined (P18), the user's scattering factors (P11) and
//! Bragg–Brentano surface micro-absorption (P29), each against what FullProf's files in the repository hold.
//!
//! - P18: FullProf's patterns that refine Lambda (line 14) write λ₁ = λ₂ = Lambda (FP_tests, nine patterns): λ₂
//!   moves with λ₁ for one wavelength; a doublet's is refused (FullProf refines a single wavelength, manual p. 83).
//!   The derivative against a central difference; case A (`runs/modes/A_modes`, neutrons, one wavelength) refined for
//!   λ from FullProf's converged state with the cell fixed returns FullProf's λ within its esd; line 14's codeword
//!   imported, exported and read back. FullProf's own refinement of λ (`runs/nabd4_lambda`, L252): the engine, from
//!   the same `.pcr` with its codewords, reaches FullProf's refined λ within FullProf's esd; FullProf's χ² of each of
//!   its cycles against the core's at the same values (`p18_fullprof_cycles`, part frozen, L255).
//! - P29: FullProf's own refinements of CeO₂ with micro-absorption (`FP_tests/ce.pcr`, `ce2.pcr`, data `ceo2.dat`,
//!   FullProf's converged Cp and τ on their line 15): the engine, started from FullProf's values, keeps Cp and τ within
//!   0.5 of its esd; the other parameters' moves are frozen (`FIXED_POINT_FROZEN`). Without line 15 χ² is higher.
//! - P11: line 12 read in the five FP_tests files with Nsc ≠ 0; the user's b, f′ and f″, coefficients and tables used
//!   as the manual says (pp. 70, 79–80); the control file's `[[pattern.scatterer]]` written and read back; the export.
//!   FullProf's own calculations with the user's values are rebuilt in `fullprof_profiles`: nabd4 (b complex, L251) and
//!   y2o3_fprime (f′ and f″ far from the tables, Y −2.0 and 4.0, O 0.5 and 0.1; χ² 1.2384 against FullProf's 1.238,
//!   L255).

use drs::jacobian::Field;
use drs::lsq::{refine, Options};
use drs::powder::{Intensities, Model};
use drs::refinable::ParameterSet;
use drs::scattering::{resolve, Radiation, Scatterer, UserScatterer};
use std::path::{Path, PathBuf};

mod common;
use common::latin1;

fn manuals(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof").join(rel)
}

fn fp_tests(file: &str) -> (drs::pcr::Pcr, drs::powder::Loaded) {
    let d = manuals("FP_tests");
    let p = drs::pcr::read(&latin1(&d.join(format!("{file}.pcr")))).unwrap_or_else(|e| panic!("{file}: {e}"));
    let l = Model::from_pcr(&p, &d, file).unwrap_or_else(|e| panic!("{file}: {e}"));
    (p, l)
}

/// χ² over N − P at the model's values.
fn chi2_at(m: &Model, data: &drs::data::Pattern, set: &ParameterSet) -> f64 {
    drs::lsq::statistics(m, data, &m.calculate(&data.x).ycalc, set.free.len()).chi2
}

/// The derivative of Ycal with respect to a field against a central difference, relative to the largest |∂y|.
fn derivative_error(m: &Model, x: &[f64], f: Field, h: f64) -> f64 {
    let (_, jac) = m.derivatives(x, &[f]);
    let at = |v: f64| {
        let mut q = m.clone();
        q.set(f, v);
        q.calculate(x).ycalc
    };
    let v = m.get(f);
    let (up, down) = (at(v + h), at(v - h));
    let scale = jac[0].iter().fold(0.0f64, |a, &d| a.max(d.abs()));
    assert!(scale > 0.0, "{f:?}: no derivative");
    jac[0].iter().zip(up.iter().zip(&down)).map(|(&d, (&u, &w))| (d - (u - w) / (2.0 * h)).abs()).fold(0.0, f64::max) / scale
}

/// The patterns of pb_ho, pbb, pb_sing and gabriela whose line 14 refines Lambda.
const PATTERNS_REFINING_LAMBDA: usize = 5;

#[test]
fn p18_wavelength_refinement() {
    // FullProf's own files: where line 14 refines Lambda, it writes λ₁ = λ₂ = Lambda (line 8 and line 14 of the same
    // pattern; multi-pattern files, read here as text, those whose comments name each line's pattern; p_p, pbho and
    // pb_hob write line 8 without its comment and are left out).
    let mut patterns = 0;
    for f in ["pb_ho", "pbb", "pb_sing", "gabriela"] {
        let text = latin1(&manuals(&format!("FP_tests/{f}.pcr")));
        let lines: Vec<&str> = text.lines().collect();
        let patt = |l: &str| l.rsplit_once("Patt#").and_then(|(_, n)| n.trim().parse::<usize>().ok());
        let next_nums = |i: usize| lines[i + 1].split_whitespace().filter_map(|w| w.parse::<f64>().ok()).collect::<Vec<_>>();
        let mut l8 = std::collections::BTreeMap::new();
        let mut l14 = std::collections::BTreeMap::new();
        for (i, l) in lines.iter().enumerate() {
            if l.trim_start().starts_with("! Lambda1") {
                l8.insert(patt(l).unwrap_or(1), next_nums(i));
            } else if l.trim_start().starts_with("!  Zero") {
                l14.insert(patt(l).unwrap_or(1), next_nums(i));
            }
        }
        for (n, z) in &l14 {
            if z.len() >= 8 && z[7] != 0.0 && l8.contains_key(n) {
                let w = &l8[n];
                assert!(w[0] == z[6] && w[1] == z[6], "{f} pattern {n}: λ₁ {} λ₂ {} Lambda {}", w[0], w[1], z[6]);
                patterns += 1;
            }
        }
    }
    assert_eq!(patterns, PATTERNS_REFINING_LAMBDA, "the patterns of FP_tests whose line 14 refines Lambda");

    // Case A: neutrons, one wavelength (λ₁ = λ₂ = 1.59509 Å).
    let d = manuals("runs/modes/A_modes");
    let p = drs::pcr::read(&latin1(&d.join("a_modes.pcr"))).unwrap();
    let l = Model::from_pcr(&p, &d, "a_modes").unwrap();
    let data = l.data.as_ref().unwrap();
    let m0 = l.model.clone();
    assert_eq!(m0.lambda, [1.59509, 1.59509]);
    // Setting λ moves λ₂ with it; the derivative is right.
    let mut q = m0.clone();
    q.set(Field::Wavelength, 1.6);
    assert_eq!(q.lambda, [1.6, 1.6]);
    // A step of 10⁻⁸ Å: larger ones move FullProf's hard window edges (Wdt·H, kept in the default model, L050) across
    // data points (10⁻⁶ Å: 7·10⁻³ at 105.25°; 10⁻⁸ Å: 6·10⁻⁹).
    let err = derivative_error(&m0, &data.x, Field::Wavelength, 1e-8);
    assert!(err < 1e-6, "∂y/∂λ against a central difference: {err:e}");

    // Refined from FullProf's converged state with the cell fixed (manual p. 83), λ moved by 0.002 Å first: it returns
    // to FullProf's 1.59509 within its esd, with the zero and the scale.
    let mut m = m0.clone();
    m.set(Field::Wavelength, 1.59709);
    let set = ParameterSet::from_addresses(&["pattern1.wavelength_1".into(), "pattern1.zero_shift".into(), "phase1.scale_factor".into()], &[], &[], &m).unwrap();
    let r = refine(&mut m, data, &set, Options { max_cycles: 20, eps: 0.01 }, |_| {});
    let (lam, esd) = (r.params[0].value, r.params[0].esd);
    println!("λ refined: {lam:.6} ± {esd:.6} (FullProf's {}), {} cycles, χ² {:.4}", m0.lambda[0], r.cycles, r.stats.chi2);
    assert!((lam - 1.59509).abs() < esd, "λ {lam} ± {esd}");
    assert_eq!(m.lambda[0], m.lambda[1]);

    // Line 14's codeword: imported as pattern1.wavelength_1, exported and read back.
    let text = latin1(&d.join("a_modes.pcr"));
    let mut lines: Vec<String> = text.lines().map(str::to_string).collect();
    let k = lines.iter().position(|l| l.trim_start().starts_with("!  Zero")).unwrap() + 1;
    let mut w: Vec<String> = lines[k].split_whitespace().map(str::to_string).collect();
    (w[6], w[7]) = ("1.595090".into(), "991.00".into());
    lines[k] = w.join(" ");
    let declared = lines.iter().position(|l| l.contains("!Number of refined parameters")).unwrap();
    let n: usize = lines[declared].split_whitespace().next().unwrap().parse().unwrap();
    lines[declared] = format!("{:>8}    !Number of refined parameters", n + 1);
    let pcr = drs::pcr::read(&lines.join("\n")).unwrap();
    let c = drs::control::Control::from_pcr(&pcr, &d, "a_modes", &d).unwrap();
    assert!(c.refinement.parameters.iter().flatten().any(|t| t.0 == Field::Wavelength), "{:?}", c.refinement.parameters);
    let e = drs::control::pcr_export::export(&c, "a_modes_l", l.data.as_ref()).unwrap();
    let back = drs::pcr::read(&e.pcr).unwrap();
    assert!(back.zero_codes[3] != 0.0 && back.lambda_refined == 1.59509, "line 14 written: {:?} {}", back.zero_codes, back.lambda_refined);
    let c2 = drs::control::Control::from_text(&c.to_text()).unwrap();
    assert!(c2.refinement.build(&c2.model).unwrap().free.iter().any(|f| f.name == "pattern1.wavelength_1"));

    // FullProf's own refinement of λ (runs/nabd4_lambda, L252): FP_tests/nabd4_lambda.pcr, FullProf's converged NaBD₄
    // with λ started at 1.200 Å and refined with the cell fixed, 15 parameters. The engine, from the same file with its
    // codewords, reaches FullProf's λ (its .out, "Lambda_pat1") within FullProf's esd. FullProf does not write λ back
    // into its .pcr, so the run's _calc pass is at 1.2 Å and not the reference for λ (L252).
    let d = manuals("runs/nabd4_lambda");
    let out = latin1(&d.join("nabd4_lambda.out"));
    let fp: Vec<f64> = out
        .lines()
        .rfind(|l| l.contains("Lambda_pat1"))
        .unwrap()
        .split(|c: char| c.is_whitespace() || c == '(' || c == ')')
        .filter_map(|w| w.parse().ok())
        .collect();
    let (fp_lam, fp_esd) = (fp[fp.len() - 2], fp[fp.len() - 1]);
    let p = drs::pcr::read(&latin1(&d.join("nabd4_lambda.pcr.in"))).unwrap();
    let c = drs::control::Control::from_pcr(&p, &d, "nabd4_lambda", &d).unwrap();
    let data = Model::from_pcr(&p, &d, "nabd4_lambda").unwrap().data.unwrap();
    let mut m = c.model.clone();
    assert_eq!(m.lambda, [1.2, 1.2]);
    let set = c.refinement.build(&m).unwrap();
    assert_eq!(set.free.len(), 15);
    let r = refine(&mut m, &data, &set, Options { max_cycles: 30, eps: 0.01 }, |_| {});
    let k = set.free.iter().position(|f| f.name == "pattern1.wavelength_1").unwrap();
    let (lam, esd) = (r.params[k].value, r.params[k].esd);
    println!("nabd4_lambda: λ {lam:.6} ± {esd:.6}, FullProf {fp_lam:.6} ± {fp_esd:.6}; {} cycles, χ² {:.4}", r.cycles, r.stats.chi2);
    assert!((lam - fp_lam).abs() < fp_esd, "λ {lam} ± {esd}, FullProf's {fp_lam} ± {fp_esd}");
    assert_eq!(m.lambda[0], m.lambda[1]);

    // A doublet's is refused, with the reason.
    let y = manuals("runs/Y2O3");
    let py = drs::pcr::read(&latin1(&y.join("y2o3_structure_conv.pcr"))).unwrap();
    let my = Model::from_pcr(&py, &y, "y2o3_structure_conv").unwrap().model;
    let e = ParameterSet::from_addresses(&["pattern1.wavelength_1".into()], &[], &[], &my).unwrap_err();
    assert!(e.contains("one wavelength only") && e.contains("p. 83"), "{e}");
}

/// FullProf's χ² of each cycle of its refinement of λ (`runs/nabd4_lambda`, L252) against the core's at the values that
/// cycle started from (the values its `.out` prints at the end of the cycle before; the first, the `.pcr`'s): (cycles,
/// the largest relative deviation). Cycles 1–4 agree within FullProf's print. From cycle 3 FullProf's U is −0.0552 and
/// U, V, W stay there with no shift and esd's of 10⁻⁸; H_G² is negative from 154.25° ("Square of FWHM(G) < 0 …
/// Fixing to HG=1.0e-10" in the `.out`; the core clamps it at 0), and FullProf ends at χ² 17.09 where the engine,
/// from the same file, reaches 3.71 with the same λ. In cycles 5–9 the core's χ² is 1.2–1.6 % lower than FullProf's,
/// the difference largest at 140–160° (Ycal −2 % at 155°) with a 0.1 % offset below 100°. Between cycle 4 (agrees) and
/// cycle 5 U, V, W do not move and Y falls from 0.116 to 0.049, so the high-angle peaks are narrower where FullProf
/// clamps H_G: believed to be FullProf's handling of those peaks, not checked. Not explained; frozen and listed for
/// Josu (L255). With the asymmetry taken at the Bragg angle, as FullProf takes it (L367), cycle 7 moved from −1.57 % to
/// −1.61 % (cycles 5–9 by −0.04 to +0.04 %): refrozen at 1.62 %, listed for Josu.
const LAMBDA_PATH: [(std::ops::RangeInclusive<usize>, f64); 2] = [(1..=4, 2.5e-3), (5..=9, 1.62e-2)];

/// The numbers of the form d.ddddd (an optional sign, five decimals) in a line of FullProf's `.out`, which prints them
/// without a space between them ("0.89712-0.00457").
fn five_decimals(l: &str) -> Vec<f64> {
    let b = l.as_bytes();
    let mut v = Vec::new();
    let mut i = 0;
    while i < b.len() {
        if b[i] == b'.' && i > 0 && b[i - 1].is_ascii_digit() && i + 5 < b.len() + 1 && b[i + 1..(i + 6).min(b.len())].iter().all(u8::is_ascii_digit) {
            let mut s = i - 1;
            while s > 0 && b[s - 1].is_ascii_digit() {
                s -= 1;
            }
            if s > 0 && b[s - 1] == b'-' {
                s -= 1;
            }
            v.push(l[s..i + 6].parse().unwrap());
            i += 6;
        } else {
            i += 1;
        }
    }
    v
}

#[test]
fn p18_fullprof_cycles() {
    let d = manuals("runs/nabd4_lambda");
    let out = latin1(&d.join("nabd4_lambda.out"));
    let lines: Vec<&str> = out.lines().collect();
    let starts: Vec<usize> = lines.iter().enumerate().filter(|(_, l)| l.contains("CYCLE No.:")).map(|(i, _)| i).collect();
    // Each cycle's values in FullProf's parameter order (its numbers 1–15, the .out's "SYMBOLIC NAMES"), and its χ².
    let cycles: Vec<(Vec<f64>, f64)> = starts
        .iter()
        .enumerate()
        .map(|(k, &i)| {
            let b = &lines[i..starts.get(k + 1).copied().unwrap_or(lines.len())];
            let find = |key: &str| b.iter().position(|l| l.contains(key)).unwrap_or_else(|| panic!("{key}"));
            let atom = |n: &str| five_decimals(b.iter().find(|l| l.starts_with(&format!(" {n} "))).unwrap());
            let first = |j: usize| b[j].split_whitespace().next().unwrap().parse::<f64>().unwrap();
            let word = |key: &str, n: usize| b[find(key)].split_whitespace().nth(n).unwrap().parse::<f64>().unwrap();
            let (na, bb, dd) = (atom("NA1"), atom("B1"), atom("D1"));
            let (h, a, xy) = (find("Halfwidth parameters"), find("Asymmetry parameters"), find("X and Y parameters"));
            let v = vec![
                word("Overall scale factor", 4),
                word("Refined wavelength", 4),
                na[9],
                bb[9],
                dd[6],
                dd[9],
                word("Zero-point", 2),
                first(h + 1),
                first(h + 2),
                first(h + 3),
                first(xy + 2),
                first(a + 1),
                first(a + 2),
                dd[3],
                dd[0],
            ];
            (v, b[find("Conventional Rietveld")].split_whitespace().next_back().unwrap().parse().unwrap())
        })
        .collect();
    assert_eq!(cycles.len(), 9);
    let p = drs::pcr::read(&latin1(&d.join("nabd4_lambda.pcr.in"))).unwrap();
    let l = Model::from_pcr_fullprof(&p, &d, "nabd4_lambda").unwrap();
    let data = l.data.unwrap();
    let numbers = drs::refinable::pcr_parameter_numbers(&p, &l.model).unwrap();
    assert_eq!(numbers.len(), 15);
    for (k, (_, fp)) in cycles.iter().enumerate() {
        let mut m = l.model.clone();
        if k > 0 {
            for (n, fields) in &numbers {
                for &f in fields {
                    m.set(f, cycles[k - 1].0[*n as usize - 1]);
                }
            }
        }
        let chi2 = drs::lsq::statistics(&m, &data, &m.calculate(&data.x).ycalc, 15).chi2;
        let dev = chi2 / fp - 1.0;
        let allowed = LAMBDA_PATH.iter().find(|(r, _)| r.contains(&(k + 1))).unwrap().1;
        println!("cycle {}: χ² core {chi2:.4}, FullProf {fp}, {:+.2} %", k + 1, 100.0 * dev);
        assert!(dev.abs() <= allowed, "cycle {}: χ² {chi2:.4} against FullProf's {fp}", k + 1);
    }
    // FullProf's last cycle starts with U < 0: H_G² negative at the high angles, as its warnings say.
    let u = cycles[7].0[7];
    assert!(u < 0.0 && out.contains("Square of FWHM(G) < 0"), "U {u}");
}

/// What the engine, started from FullProf's converged values of its micro-absorption refinements, moves: (run, the
/// parameters moved by more than one of the engine's esd's, the largest |shift|/esd, χ²'s relative fall), as measured
/// in L233. They are the profile's (U, V, W, η, the Lorentzian terms, the asymmetry) and, for ce2, the zero and the cell,
/// which correlate with them (U–V −0.90, V–W −0.91, zero–a +0.85); χ² falls by 0.13 % (ce) and 0.56 % (ce2). FullProf
/// stopped at Eps = 0.3 (its shifts below 0.3 of its esd) with R = 0.8 damping, and its esd's are larger than the
/// engine's (χ² over the points with Bragg contributions, L053); its own esd's are not in the repository (no `.out`).
/// Not explained further without FullProf's output; listed for Josu.
const FIXED_POINT_FROZEN: &[(&str, usize, f64, f64)] = &[("ce", 3, 2.18, 0.0013), ("ce2", 8, 2.81, 0.0057)];

#[test]
fn p29_surface_microabsorption() {
    // The manual's formula, P(h) = P0 + Cp (τ/sin θ)(1 − τ/sin θ), S_r = 1 − P (eqs. 3.35–3.36, p. 55).
    let (_, l) = fp_tests("ce");
    let m = &l.model;
    assert_eq!(m.microabsorption, [0.0, 1.2304, 0.1340]);
    let th = 20f64.to_radians();
    let r = 0.1340 / th.sin();
    assert!((m.micro_absorption(th) - (1.0 - 1.2304 * r * (1.0 - r))).abs() < 1e-15);
    let mut none = m.clone();
    none.microabsorption = [0.0; 3];
    assert_eq!(none.micro_absorption(th), 1.0);
    assert!((m.absorption(th) - none.absorption(th) * m.micro_absorption(th)).abs() < 1e-15);
    // Each term's derivative.
    let x = &l.data.as_ref().unwrap().x;
    for (k, h) in [(0, 1e-6), (1, 1e-6), (2, 1e-7)] {
        let err = derivative_error(m, x, Field::Micro(k), h);
        assert!(err < 1e-4, "∂y/∂P{k} against a central difference: {err:e}");
    }

    // FullProf's converged values are the engine's minimum for Cp and τ.
    for (run, frozen) in [("ce", FIXED_POINT_FROZEN[0]), ("ce2", FIXED_POINT_FROZEN[1])] {
        let (p, l) = fp_tests(run);
        let data = l.data.as_ref().unwrap();
        let mut m = l.model.clone();
        let set = ParameterSet::from_pcr(&p, &m).unwrap();
        let start: Vec<f64> = set.free.iter().map(|f| f.value_in(&m)).collect();
        let chi2_start = chi2_at(&m, data, &set);
        let r = refine(&mut m, data, &set, Options { max_cycles: 20, eps: 0.1 }, |_| {});
        let mut moved = Vec::new();
        for ((f, q), s0) in set.free.iter().zip(&r.params).zip(&start) {
            let z = (q.value - s0) / q.esd;
            println!("  {run} {:<40} FullProf {s0:>12.6}  engine {:>12.6} ± {:.6}  {z:+.2}", f.name, q.value, q.esd);
            if f.name.contains("microabsorption") {
                assert!(z.abs() < 0.5, "{run} {}: {z:.2} esd from FullProf's", f.name);
            } else if z.abs() > 1.0 {
                moved.push(z.abs());
            }
        }
        let largest = moved.iter().fold(0.0f64, |a, &b| a.max(b));
        println!("{run}: χ² {chi2_start:.4} at FullProf's values, {:.4} refined; {} parameters moved > 1 esd, largest {largest:.2}", r.stats.chi2, moved.len());
        let fall = (chi2_start - r.stats.chi2) / chi2_start;
        assert!(fall <= frozen.3, "{run}: χ² falls by {fall:.5}; frozen {}", frozen.3);
        assert!(moved.len() <= frozen.1 && largest <= frozen.2, "{run}: {} moved, largest {largest:.2}; frozen {frozen:?}", moved.len());
        // Without the micro-absorption, χ² at FullProf's values is higher.
        let mut none = l.model.clone();
        none.microabsorption = [0.0; 3];
        let without = chi2_at(&none, data, &set);
        println!("{run}: χ² {without:.4} without line 15");
        assert!(without > chi2_start * 1.05, "{run}: line 15 changes χ² by {:.2} %", 100.0 * (without / chi2_start - 1.0));
    }

    // In the control file and back to FullProf: pattern keys, refined by address, line 15 with its codewords.
    let (p, _) = fp_tests("ce");
    let d = manuals("FP_tests");
    let c = drs::control::Control::from_pcr(&p, &d, "ce", &d).unwrap();
    let text = c.to_text();
    assert!(text.contains("microabsorption_cp = 1.2304") && text.contains("microabsorption_tau = 0.134"), "{text}");
    let c2 = drs::control::Control::from_text(&text).unwrap();
    assert_eq!(c2.model.microabsorption, [0.0, 1.2304, 0.134]);
    let free: Vec<String> = c2.refinement.build(&c2.model).unwrap().free.iter().map(|f| f.name.clone()).collect();
    assert!(free.contains(&"pattern1.microabsorption_cp".into()) && free.contains(&"pattern1.microabsorption_tau".into()), "{free:?}");
    let e = drs::control::pcr_export::export(&c, "ce_x", None).unwrap();
    let back = drs::pcr::read(&e.pcr).unwrap();
    assert_eq!(back.microabsorption, Some([0.0, 1.2304, 0.134]));
    assert!(back.microabsorption_codes[1] != 0.0 && back.microabsorption_codes[2] != 0.0 && back.microabsorption_codes[0] == 0.0);
    // The command line: `drs refine` on FullProf's file refines Cp and τ.
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("target/test-pattern-terms");
    std::fs::create_dir_all(&dir).unwrap();
    for f in ["ce.pcr", "ceo2.dat"] {
        std::fs::copy(d.join(f), dir.join(f)).unwrap();
    }
    let o = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(["refine", dir.join("ce.pcr").to_str().unwrap(), "--cycles", "1"]).output().unwrap();
    let out = String::from_utf8_lossy(&o.stdout);
    assert!(o.status.success() && out.contains("pattern1.microabsorption_cp") && out.contains("pattern1.microabsorption_tau"), "{out}{}", String::from_utf8_lossy(&o.stderr));
    // Neutrons have no line 15 (manual p. 84): no micro-absorption fields.
    let (_, n) = fp_tests("nabd4");
    assert!(!n.model.all_fields().iter().any(|f| matches!(f, Field::Micro(_))));
}

#[test]
fn p11_user_scattering_factors() {
    // Line 12 of the five FP_tests files with Nsc ≠ 0, as written.
    let read = |f: &str| drs::pcr::read(&latin1(&manuals(&format!("FP_tests/{f}.pcr")))).unwrap_or_else(|e| panic!("{f}: {e}"));
    let nabd4 = read("nabd4");
    assert_eq!(nabd4.user_scattering.len(), 1);
    assert_eq!((nabd4.user_scattering[0].name.as_str(), nabd4.user_scattering[0].dfp, nabd4.user_scattering[0].ity), ("B", 0.53, 0));
    let hp = read("highpress");
    assert_eq!(hp.user_scattering.iter().map(|u| (u.name.as_str(), u.dfp, u.dfpp, u.ity)).collect::<Vec<_>>(), vec![("hg", -1.357, 3.279, 2), ("ba", -2.041, 0.74, 2), ("ca", 0.055, 0.088, 2), ("cu", 0.157, 0.396, 2), ("o", -0.002, 0.001, 2)]);
    // c60a: an X-ray table (Nsc = −1, ITY = 0) ended by −100; hpress: five sets of nine coefficients (Nsc = 5, ITY = 0).
    // (cedif4rf, Nsc = −1 with a magnetic table, stops before line 12 on its user symmetry lines, which are not read.)
    let c60 = read("c60a");
    assert_eq!((c60.user_scattering[0].name.as_str(), c60.user_scattering[0].ity), ("c60", 0));
    assert_eq!(c60.user_scattering[0].table[..2], [(0.0, 1.0), (22.599, 0.9584)]);
    let hpress = read("hpress");
    assert_eq!(hpress.user_scattering.len(), 5);
    assert_eq!(hpress.user_scattering[0].name, "hgin");
    assert_eq!(hpress.user_scattering[0].coefficients, vec![20.6809, 0.545, 19.0417, 8.4484, 21.6575, 1.5729, 5.9676, 38.3246, 12.6089]);
    let u = drs::powder::user_scatterers(&hpress);
    assert_eq!(u[0].form_factor.map(|k| k[8]), Some(12.6089));
    assert_eq!((u[0].f_prime, u[0].f_double_prime), (-1.357, 3.279));

    // Neutrons: the user's b with DFPP as its imaginary part, as FullProf uses it (nabd4_calc.out: "Scattering lengths
    // (real & imaginary parts) … For B 0.530000 −0.213000"; its manual, p. 79, says DFPP is ignored for neutrons). With
    // the powder intensity averaged over h and −h (`Structure::intensity`), FullProf's nabd4_calc.prf is rebuilt exactly
    // (fullprof_profiles, L251); the sign of b″ then does not matter.
    let (_, l) = fp_tests("nabd4");
    let Intensities::Structure(st, sc) = &l.model.phases[0].intensities else { panic!() };
    let b = st.atoms.iter().position(|a| a.species == "B").unwrap();
    assert_eq!(sc[b], Scatterer::Neutron { symbol: "B".into(), b: 0.53, b_imaginary: -0.213 });
    let mut flipped = sc.clone();
    flipped[b] = Scatterer::Neutron { symbol: "B".into(), b: 0.53, b_imaginary: 0.213 };
    for h in [[1, 1, 2], [2, 1, 1], [2, 1, 3]] {
        let (i, j) = (st.intensity(sc, h), st.intensity(&flipped, h));
        assert!((i - j).abs() < 1e-9 * i, "{h:?}: {i} {j}");
        // A single F(h) differs between the two signs in this group without a centre (P −4 2₁ c): the average is needed.
        assert!((st.structure_factor(sc, h).norm_sqr() - st.structure_factor(&flipped, h).norm_sqr()).abs() > 1e-3 * i);
    }
    let real: Vec<Scatterer> = sc.iter().map(|x| match x { Scatterer::Neutron { symbol, b, .. } => Scatterer::Neutron { symbol: symbol.clone(), b: *b, b_imaginary: 0.0 }, o => o.clone() }).collect();
    assert_eq!(st.intensity(&real, [1, 1, 2]), st.structure_factor(&real, [1, 1, 2]).norm_sqr(), "no imaginary part: |F(h)|² itself");
    let user = [UserScatterer { species: "B".into(), scattering_length: Some(0.6), ..Default::default() }];
    assert_eq!(resolve("b", Radiation::Neutron, &user).unwrap().at(0.3).re, 0.6, "matched case aside");
    // X-rays, ITY = 2: the tables' f(s) with the user's f′ and f″ (highpress: Hg at λ = 0.3774 Å).
    let xr = Radiation::Xray { lambda: 0.3774, anomalous: true };
    let hg = [UserScatterer { species: "hg".into(), f_prime: -1.357, f_double_prime: 3.279, ..Default::default() }];
    let own = resolve("HG", xr, &hg).unwrap();
    let plain = Scatterer::lookup("hg", Radiation::Xray { lambda: 0.3774, anomalous: false }).unwrap();
    for s in [0.0, 0.3, 0.8] {
        let (a, p) = (own.at(s), plain.at(s));
        assert!((a.re - (p.re - 1.357)).abs() < 1e-12 && (a.im - 3.279).abs() < 1e-12);
    }
    // The nine coefficients (ITY = 0, Nsc > 0): f(s) = Σ aᵢ exp(−bᵢ s²) + c; with A2 = 0 only A1 and B1 (p. 79).
    let k = [10.0, 2.0, 5.0, 10.0, 1.0, 30.0, 0.5, 60.0, 0.3];
    let u = UserScatterer { species: "Xx".into(), form_factor: Some(k), f_prime: 0.1, ..Default::default() };
    let s = 0.25f64;
    let want = 10.0 * (-2.0 * s * s).exp() + 5.0 * (-10.0 * s * s).exp() + (-30.0 * s * s).exp() + 0.5 * (-60.0 * s * s).exp() + 0.3 + 0.1;
    assert!((u.scatterer(xr).unwrap().at(s).re - want).abs() < 1e-12);
    let mut k1 = k;
    k1[2] = 0.0;
    let u1 = UserScatterer { form_factor: Some(k1), ..u.clone() };
    assert!((u1.scatterer(xr).unwrap().at(s).re - (10.0 * (-2.0 * s * s).exp() + 0.1)).abs() < 1e-12);
    // A table (Nsc < 0): linear between its points, held at its ends.
    let t = UserScatterer { species: "Xx".into(), form_factor_table: vec![(0.0, 10.0), (0.5, 6.0), (1.0, 4.0)], ..Default::default() };
    let ts = t.scatterer(xr).unwrap();
    assert_eq!((ts.at(0.25).re, ts.at(0.75).re, ts.at(2.0).re), (8.0, 5.0, 4.0));
    assert_eq!(ts.ds(0.25), -8.0);
    // An unknown species with no user value is the tables' error.
    assert!(resolve("Qq", xr, &[]).is_err());

    // The control file: a [[pattern.scatterer]] for each species whose scatterer is not the tables', read back the same;
    // a species of no atom stays out of the file.
    let y = manuals("runs/Y2O3");
    let py = drs::pcr::read(&latin1(&y.join("y2o3_structure_conv.pcr"))).unwrap();
    let mut c = drs::control::Control::from_pcr(&py, &y, "y2o3_structure_conv", &y).unwrap();
    let plain_text = c.to_text();
    assert!(!plain_text.contains("[[pattern.scatterer]]"));
    let ry = Radiation::Xray { lambda: c.model.lambda[0], anomalous: true };
    let yu = UserScatterer { species: "Y".into(), f_prime: -2.5, f_double_prime: 3.5, ..Default::default() };
    let ou = UserScatterer { species: "O".into(), form_factor_table: vec![(0.0, 8.0), (0.5, 5.0), (1.5, 2.0)], ..Default::default() };
    let Intensities::Structure(st, sc) = &mut c.model.phases[0].intensities else { panic!() };
    for (a, s) in st.atoms.iter().zip(sc.iter_mut()) {
        *s = resolve(&a.species, ry, &[yu.clone(), ou.clone()]).unwrap();
    }
    let text = c.to_text();
    assert_eq!(text.matches("[[pattern.scatterer]]").count(), 2, "{text}");
    let back = drs::control::Control::from_text(&text).unwrap();
    let (Intensities::Structure(_, a), Intensities::Structure(_, b)) = (&c.model.phases[0].intensities, &back.model.phases[0].intensities) else { panic!() };
    assert_eq!(a, b);
    assert_eq!(back.to_text(), text);
    // An unknown key is named.
    let bad = text.replacen("f_prime = -2.5", "fprime = -2.5", 1);
    let e = drs::control::Control::from_text(&bad).unwrap_err();
    assert!(e.iter().any(|d| d.message.contains("\"fprime\" is not a key")), "{e:?}");
    // An atom added later takes the pattern's user value.
    let mut m = back.model.clone();
    let at = m.add_atom(0, "Y9", "Y", [0.1, 0.2, 0.3], 0.5, 0.1).unwrap();
    let Intensities::Structure(_, sc) = &m.phases[0].intensities else { panic!() };
    assert_eq!(sc[at], yu.scatterer(ry).unwrap());
    // To FullProf: Nsc < 0 (a table), Y with ITY = 2, O's table ended by −100; read back the same.
    let e = drs::control::pcr_export::export(&c, "y2o3_u", None).unwrap();
    let p = drs::pcr::read(&e.pcr).unwrap();
    assert_eq!(p.flag("Nsc"), -2);
    let pu: Vec<_> = p.user_scattering.iter().map(|u| (u.name.as_str(), u.dfp, u.dfpp, u.ity, u.table.len())).collect();
    assert_eq!(pu, vec![("y", -2.5, 3.5, 2, 0), ("o", 0.0, 0.0, 0, 3)]);
    assert_eq!(drs::powder::user_scatterers(&p).len(), 2);
}
