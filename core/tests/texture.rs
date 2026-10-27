//! Preferred orientation along several axes (P30, L360) against FullProf 8.20's runs of `runs/version2_runs.py`
//! (`manuals/fullprof/runs/texture`, run on the Mac: L346): NCY = 0 variants of Yb2ZrWO8-PO-March (Fm-3m), whose
//! `.prf` belong to their `VARIANT.pcr.in` (`fullprof_profiles` compares the four profiles).
//!
//! - `against_fullprof_corr`: FullProf's CORR (the `.out`'s last column, the intensity's corrections) of each
//!   reflection over yb_nor2_one's (Nor = 2 with the one axis [1 1 0], r = 0.85), against the core's preferred
//!   orientation over its one-axis value. yb_nor2 (110 r 0.85 f 0.6, 001 r 1.2 f 0.4) and yb_nor2_g2 (the same with
//!   Pref2 = 0.3) agree at all 16 reflections within the print (CORR's 6 decimals, 2·10⁻⁶): Σ_j f_j D_j averaged over the
//!   equivalents, G2 not entering (yb_nor2_g2's CORR is yb_nor2's). yb_nor3 (Nor = 3, 20 steps) differs at all 16, the
//!   core's up to 3.75% above FullProf's (FullProf's 3.6% below, at (3 3 3)): frozen in `NOR3`, listed for Josu (L360). Wright et al. (1996) eqs. 4–6, which
//!   the News names, are not in `manuals/papers`; the core's cone integral (`texture`) is believed, and no reading
//!   tried (ψ = 90° − θ, θ, 2θ; φ over π, 2π, π/2; weights 1, sin φ; the equivalents averaged or not; the beam
//!   perpendicular to the force, integrated over the ring; per-axis factors) gives FullProf's values, which fall
//!   outside what the model reaches at any ψ for 13 of the 16 reflections.
//! - `derivatives`: the Jacobian's columns of each axis's r and fraction against central differences of the pattern.
//! - `refinement`: from both ratios and the scale moved, the engine recovers yb_nor2's from the core's own pattern.
//! - `control_file_and_export`: the control file and the `.pcr` export of yb_nor2 and yb_nor3, read back to the same
//!   texture and pattern; each value located in the text; the results file's numeric items (P48.1–2).

use drs::jacobian::Field;
use drs::lsq::{refine, Options};
use drs::powder::Model;
use drs::refinable::Refinement;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

fn dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs/texture")
}

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display())).iter().map(|&b| b as char).collect()
}

fn load(run: &str) -> drs::powder::Loaded {
    let p = drs::pcr::read(&latin1(&dir().join(format!("{run}.pcr.in")))).unwrap_or_else(|e| panic!("{run}: {e}"));
    Model::from_pcr_fullprof(&p, &dir(), run).unwrap_or_else(|e| panic!("{run}: {e}"))
}

/// The reflection list of a run's `.out`: h k l → (2θ, CORR).
fn corr(run: &str) -> BTreeMap<[i32; 3], (f64, f64)> {
    let text = latin1(&dir().join(format!("{run}.out")));
    let lines: Vec<&str> = text.lines().collect();
    let head = lines.iter().rposition(|l| l.trim_start().starts_with("No.  Code") && l.contains("CORR")).unwrap();
    let mut out = BTreeMap::new();
    for l in &lines[head + 1..] {
        let t: Vec<&str> = l.split_whitespace().collect();
        if t.len() < 15 {
            if out.is_empty() {
                continue;
            }
            break;
        }
        out.insert([t[2].parse().unwrap(), t[3].parse().unwrap(), t[4].parse().unwrap()], (t[7].parse().unwrap(), t[t.len() - 1].parse().unwrap()));
    }
    out
}

/// yb_nor3 (Nor = 3) against FullProf's CORR, frozen (L360): the reflections off by more than the print, and the
/// largest relative difference. It fails if either grows.
const NOR3: (usize, f64) = (16, 0.03752);

#[test]
fn against_fullprof_corr() {
    let one = load("yb_nor2_one").model;
    let fp_one = corr("yb_nor2_one");
    assert_eq!(fp_one.len(), 16);
    let po = |m: &Model, h: [i32; 3], tt: f64| m.preferred_orientation_at(&m.phases[0], h.map(|v| v as f64), (tt / 2.0).to_radians());
    for (run, frozen) in [("yb_nor2", None), ("yb_nor2_g2", None), ("yb_nor3", Some(NOR3))] {
        let m = load(run).model;
        let t = m.phases[0].texture.as_ref().unwrap();
        assert_eq!((t.axes.len(), t.steps), (2, if run == "yb_nor3" { 20 } else { 0 }), "{run}");
        let fp = corr(run);
        let (mut off, mut worst) = (0, 0.0f64);
        for (h, &(tt, c)) in &fp {
            let want = c / fp_one[h].1;
            let got = po(&m, *h, tt) / po(&one, *h, tt);
            let d = (got / want - 1.0).abs();
            if d > 2e-6 {
                off += 1;
            }
            worst = worst.max(d);
        }
        println!("{run}: {} reflections, {off} off by more than the print, worst {worst:.3e}", fp.len());
        match frozen {
            None => assert_eq!(off, 0, "{run}: worst {worst:.3e}"),
            Some((n, w)) => assert!(off <= n && worst <= w, "{run}: {off} off, worst {worst:.5}, frozen at {n}, {w}"),
        }
    }
    // FullProf's G2 does not enter Nor = 2: yb_nor2_g2's CORR is yb_nor2's.
    assert_eq!(corr("yb_nor2"), corr("yb_nor2_g2"));
}

#[test]
fn derivatives() {
    for run in ["yb_nor2", "yb_nor3"] {
        let l = load(run);
        let (m, data) = (l.model, l.data.unwrap());
        let n = m.phases[0].texture.as_ref().unwrap().len();
        assert_eq!(n, 3);
        let fields: Vec<Field> = (0..n).map(|k| Field::Texture(0, k)).collect();
        assert!(fields.iter().all(|f| m.all_fields().contains(f)), "{run}: the texture's values among the model's fields");
        let (_, jac) = m.derivatives(&data.x, &fields);
        for (j, &f) in fields.iter().enumerate() {
            let h = 1e-6;
            let (mut a, mut b) = (m.clone(), m.clone());
            a.set(f, m.get(f) + h);
            b.set(f, m.get(f) - h);
            let (ya, yb) = (a.calculate(&data.x).ycalc, b.calculate(&data.x).ycalc);
            let top = jac[j].iter().fold(0.0f64, |s, v| s.max(v.abs()));
            let worst = (0..data.x.len()).map(|i| ((ya[i] - yb[i]) / (2.0 * h) - jac[j][i]).abs()).fold(0.0, f64::max);
            assert!(top > 0.0 && worst < 1e-5 * top, "{run} {f:?}: off by {worst:.3e} of {top:.3e}");
        }
    }
}

#[test]
fn refinement() {
    let l = load("yb_nor2");
    let (truth, data) = (l.model, l.data.unwrap());
    let y = truth.calculate(&data.x).ycalc;
    let synthetic = drs::data::Pattern { y: y.clone(), sigma: y.iter().map(|v| v.max(1.0).sqrt()).collect(), ..data.clone() };
    // The two ratios and the scale, the fractions fixed: with all three values free, the 16 reflections of Fm-3m do
    // not fix them (from r₁ 0.9, f₁ 0.5, r₂ 1.1 the engine reaches χ² 1.6·10⁻⁷ at r₁ 0.764, f₁ 0.159, r₂ 1.126).
    let mut m = truth.clone();
    for (k, v) in [(0, 0.95), (2, 1.05)] {
        m.set(Field::Texture(0, k), v);
    }
    m.phases[0].scale *= 1.03;
    let mut rf = Refinement { parameters: vec![vec![(Field::Texture(0, 0), 1.0)], vec![(Field::Texture(0, 2), 1.0)], vec![(Field::Scale(0), 1.0)]], ..Default::default() };
    rf.steering.guards.jump = f64::INFINITY;
    let set = rf.build(&m).unwrap();
    let r = refine(&mut m, &synthetic, &set, Options { max_cycles: 30, eps: 1e-5 }, |_| {});
    let got = [m.get(Field::Texture(0, 0)), m.get(Field::Texture(0, 2)), m.phases[0].scale / truth.phases[0].scale];
    println!("{} cycles, χ² {:.3e}: r₁ {:.7}, r₂ {:.7}, scale × {:.7}", r.cycles, r.stats.chi2, got[0], got[1], got[2]);
    for (g, w) in got.iter().zip([0.85, 1.2, 1.0]) {
        assert!((g - w).abs() < 1e-6, "{got:?}");
    }
}

#[test]
fn control_file_and_export() {
    let d = dir();
    let tmp = Path::new(env!("CARGO_MANIFEST_DIR")).join("target/test-texture");
    let _ = std::fs::remove_dir_all(&tmp);
    std::fs::create_dir_all(&tmp).unwrap();
    for (run, nor) in [("yb_nor2", 2), ("yb_nor3", 3)] {
        let p = drs::pcr::read(&latin1(&d.join(format!("{run}.pcr.in")))).unwrap();
        let l = Model::from_pcr(&p, &d, run).unwrap();
        let data = l.data.unwrap();
        // Two values free, for the codes of the export and the su's of the results.
        let mut c = drs::control::Control::from_pcr(&p, &d, run, &tmp).unwrap();
        c.refinement.parameters.push(vec![(Field::Texture(0, 0), 1.0)]);
        c.refinement.parameters.push(vec![(Field::Texture(0, 1), 1.0)]);
        let text = c.to_text();
        let back = drs::control::Control::from_text_in(&text, &tmp).unwrap_or_else(|e| panic!("{run}: {e:?}"));
        assert_eq!(back.model.phases[0].texture, c.model.phases[0].texture, "{run}: the texture read back");
        assert_eq!(back.model.calculate(&data.x).ycalc, l.model.calculate(&data.x).ycalc, "{run}: the pattern read back");
        assert_eq!(back.refinement, c.refinement, "{run}: the free parameters read back");
        // Each value where the file writes it: three of two axes.
        let located: Vec<_> = drs::control::locate::all(&text, Some(&tmp)).into_iter().filter(|(a, ..)| a.contains("po_axis")).collect();
        assert_eq!(located.len(), 3, "{run}");
        for (a, _, r) in &located {
            let f = drs::refinable::field_of(a, &c.model).unwrap();
            assert!(matches!(f, Field::Texture(0, _)), "{a}");
            assert_eq!(text[r.clone()].parse::<f64>().unwrap(), c.model.get(f), "{a}");
        }
        // The .pcr export read again: Nor, the same texture and pattern, the same free parameters.
        let e = drs::control::pcr_export::export(&c, run, Some(&data)).unwrap();
        std::fs::write(tmp.join(format!("{run}.pcr")), &e.pcr).unwrap();
        for (name, t) in &e.files {
            std::fs::write(tmp.join(name), t).unwrap();
        }
        let p2 = drs::pcr::read(&e.pcr).unwrap_or_else(|er| panic!("{run}: {er}\n{}", e.pcr));
        assert_eq!(p2.flag("Nor"), nor, "{run}");
        let m2 = Model::from_pcr(&p2, &tmp, run).unwrap_or_else(|er| panic!("{run}: {er}")).model;
        assert_eq!(m2.phases[0].texture, c.model.phases[0].texture, "{run}: the texture through the .pcr");
        let (y1, y2) = (c.model.calculate(&data.x).ycalc, m2.calculate(&data.x).ycalc);
        let worst = y1.iter().zip(&y2).map(|(a, b)| (a - b).abs() / a.abs().max(1.0)).fold(0.0, f64::max);
        assert!(worst < 1e-9, "{run}: the pattern through the .pcr, {worst:.2e}");
        let free = |rf: &Refinement, m: &Model| rf.build(m).unwrap().free.iter().map(|f| f.fields[0].0).collect::<Vec<_>>();
        // The .pcr's codewords number them in another order (the fields' in the file).
        let sorted = |mut v: Vec<Field>| {
            v.sort();
            v
        };
        assert_eq!(sorted(free(&Refinement::from_pcr(&p2, &m2).unwrap(), &m2)), sorted(free(&c.refinement, &c.model)), "{run}: the free parameters through the .pcr");
        // The results file (P48.1–2): the function, the geometry, per axis h k l, r and fraction; after one cycle the
        // su's, the last fraction's that of the first (f₂ = 1 − f₁).
        let mut m = c.model.clone();
        let set = c.refinement.build(&m).unwrap();
        let r = refine(&mut m, &data, &set, Options { max_cycles: 1, eps: 1e-3 }, |_| {});
        let mut cr = c.clone();
        cr.model = m;
        let (cif, _) = drs::results::write(&cr, &data, Some(&set), Some(&r), run).unwrap();
        assert!(cif.contains("_drs_pref_orient_function          'march-dollase axes'"), "{run}");
        assert!(cif.contains(&format!("_drs_pref_orient_geometry          {}", if nor == 3 { "'anvil cell'" } else { "reflection" })), "{run}");
        let rows: Vec<Vec<&str>> = cif.lines().skip_while(|l| !l.starts_with("_drs_pref_orient_axis_fraction")).skip(1).take(2).map(|l| l.split_whitespace().collect()).collect();
        assert_eq!((rows[0][1..4].to_vec(), rows[1][1..4].to_vec()), (vec!["1", "1", "0"], vec!["0", "0", "1"]), "{run}: {rows:?}");
        let su = |s: &str| s.split_once('(').map(|x| x.1.trim_end_matches(')').to_string());
        assert!(su(rows[0][4]).is_some() && su(rows[0][5]).is_some() && su(rows[0][5]) == su(rows[1][5]), "{run}: {rows:?}");
        assert!(su(rows[1][4]).is_none(), "{run}: r₂ is fixed: {rows:?}");
        println!("{run}: control file {} kB, 3 values located, .pcr Nor {nor}, results rows {rows:?}", text.len() / 1024);
    }
}

#[test]
fn march_dollase_of_a_phase() {
    // P48.4: a phase's axes for Publish's F10; one March–Dollase axis with G2 as two axes, the second of r = 1.
    let one = drs::texture::of_phase(&load("yb_nor2_one").model.phases[0], true).unwrap();
    assert_eq!(one.axes.len(), 1);
    assert!((drs::texture::march_dollase(0.85, 1.0) - 0.85f64.powi(-3)).abs() < 1e-12);
    let mut ph = load("yb_nor2").model.phases[0].clone();
    ph.texture = None;
    ph.pref = [0.85, 0.3];
    let t = drs::texture::of_phase(&ph, true).unwrap();
    for c in [0.0, 0.3, 0.8, 1.0] {
        let want = 0.3 + 0.7 * drs::texture::march_dollase(0.85, c);
        assert!((t.factor(&[c, c], 0.0) - want).abs() < 1e-12);
    }
    assert!(drs::texture::of_phase(&ph, false).is_none());
}

