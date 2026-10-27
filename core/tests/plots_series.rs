//! Plots' requests PL5, PL11 and PL17 (R32, L170 step 4, L177), and two findings of core 0.29.0:
//! - PL5: windows on Tb₂Pd₂Sn's four magnetic peaks against T (the mock-up F5): the integral above the automatic
//!   background with σ, the centroid, the strongest peak, the background at 16.5°, 30° and 45°; the integral against
//!   the mock-up harness's sum; σ against finite differences for the line and for given points; a synthetic peak;
//! - PL11: LaMnO₃'s `la0.seq` (71 blocks) as a table, every block against its data file's header and its `.sum`;
//! - PL17: the projection with hidden parts on the mock-up F2's geometry, each visible piece checked against a
//!   brute-force test of every curve in front; a tent and a line by hand; a floor that lets a curve show beneath;
//! - (a) the D20 reader gives no wavelength: FullProf's Ins = 3 names the header field holding 1.870 `tr2`, and D1B's
//!   hocu.dat has 0.000 there with λ = 2.524 Å;
//! - (b) σ of no weight is `null` in JSON and reads back as no weight;
//! - the command line, `drs plots windows`, `drs plots seq`, `drs plots project`, with `--json`.

use drs::data::Pattern;
use drs::plots::{self, Displaced, Displacement, Window, WindowBackground, YScale};
use std::path::{Path, PathBuf};

fn root(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)
}

fn tb() -> PathBuf {
    root("../manuals/fullprof/1-WinPLOTR/1_Tb2pd2sn")
}

fn d20() -> PathBuf {
    root("../manuals/fullprof/1-WinPLOTR/2_FP_sequential_refinements/LaMnO3_D20")
}

fn drs(args: &[&str]) -> std::process::Output {
    std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(args).output().unwrap()
}

/// The mock-up F5's windows (plots/mockups/make_mockups.py, WINDOWS), on the peaks at 19.1°, 26.7°, 32.6° and 42.2°.
fn f5_windows() -> Vec<Window> {
    [("m1", 18.5, 19.8), ("m2", 26.0, 27.4), ("m3", 31.9, 33.3), ("m4", 41.5, 42.9)].map(|(n, a, b)| Window { name: n.into(), from: a, to: b }).to_vec()
}

#[test]
fn pl5_windows_against_t() {
    let set = plots::open_buf(&tb().join("Tb_vs_T.buf"), None).unwrap();
    let ws = f5_windows();
    let w = plots::windows(&set, &ws, &WindowBackground::Auto { log_lambda: None }, &[16.5, 30.0, 45.0]).unwrap();
    assert_eq!(w.members.len(), 20);
    // Against the mock-up harness's sum (plots/mockups/dump/src/main.rs): Σ (y − b)·h over the points in [from, to],
    // σ² = Σ (σ h)². Here each point counts over the part of its interval inside the window, so the two end points
    // (on the window's limits) count half: integral = sum − h/2 (net₁ + net_n), σ² = σ²_sum − ¾ h² (σ₁² + σ_n²).
    for (mw, m) in w.members.iter().zip(&set.members) {
        let p = &m.pattern;
        let b = drs::data_ops::auto_background(p, None, drs::data_ops::DEFAULT_MAX_ITER).unwrap().curve;
        let h = p.mean_step().unwrap();
        for v in &mw.windows {
            let idx: Vec<usize> = (0..p.len()).filter(|&i| p.x[i] >= v.from - 1e-9 && p.x[i] <= v.to + 1e-9).collect();
            let (f, l) = (idx[0], idx[idx.len() - 1]);
            let sum: f64 = idx.iter().map(|&i| (p.y[i] - b[i]) * h).sum();
            let var: f64 = idx.iter().map(|&i| (p.sigma[i] * h).powi(2)).sum();
            let expect = sum - h / 2.0 * ((p.y[f] - b[f]) + (p.y[l] - b[l]));
            assert!((v.integral - expect).abs() < 1e-9 * expect.abs().max(1.0), "{} {}: {} {expect}", m.label, v.name, v.integral);
            let expect_s = (var - 0.75 * h * h * (p.sigma[f].powi(2) + p.sigma[l].powi(2))).sqrt();
            assert!((v.sigma - expect_s).abs() < 1e-9 * expect_s, "{} {}", m.label, v.name);
            assert_eq!((v.points, v.unweighted), (idx.len(), 0));
        }
        // The background at the chosen x: 16.5°, 30° and 45° are points, where it is the curve's own value.
        for (k, x) in [16.5, 30.0, 45.0].iter().enumerate() {
            let i = (0..p.len()).min_by(|&i, &j| (p.x[i] - x).abs().total_cmp(&(p.x[j] - x).abs())).unwrap();
            assert!((mw.background_at[k].unwrap() - b[i]).abs() < 1e-6, "{} {x}", m.label);
        }
    }
    // F5 against T: the members sorted 4 … 30 K, then BT.
    let t = |k: usize| set.members[k].variable.as_ref().map(|v| v.value);
    let value = |k: usize, name: &str| w.members[k].windows.iter().find(|v| v.name == name).unwrap().clone();
    for k in 0..20 {
        let row: Vec<String> = ["m1", "m2", "m3", "m4"].iter().map(|n| format!("{:9.1} ± {:4.1}", value(k, n).integral, value(k, n).sigma)).collect();
        println!("{:10} T {:?}: {}  b {:?}", set.members[k].label, t(k), row.join("  "), w.members[k].background_at);
    }
    assert_eq!((t(0), t(18), t(19)), (Some(4.0), Some(30.0), None));
    for n in ["m1", "m2", "m3", "m4"] {
        let (cold, warm) = (value(0, n), value(18, n));
        assert!(cold.integral > 100.0 * cold.sigma, "{n} at 4 K: {} ± {}", cold.integral, cold.sigma);
        assert!(warm.integral < 0.05 * cold.integral, "{n}: 30 K {} against 4 K {}", warm.integral, cold.integral);
        // The order parameter: below half its 4 K value by 22 K (k = 10).
        assert!(value(10, n).integral < 0.5 * cold.integral, "{n} at 22 K");
        // The centroid on the peak the search finds.
        assert!((cold.centroid.unwrap() - cold.peak.as_ref().unwrap().x).abs() < 0.05, "{n}");
    }
    // At 30 K the peaks at 19.1° and 32.6° are gone (below 2σ, open in F5); the windows at 26.7° and 42.2° keep the
    // nuclear peaks at 26.4° and 42.1° (mock-up README: "about 120 counts·°" at 26.7°).
    for n in ["m1", "m3"] {
        assert!(value(18, n).integral < 2.0 * value(18, n).sigma, "{n}");
    }
    for (n, x) in [("m2", 26.4), ("m4", 42.1)] {
        let v = value(18, n);
        assert!(v.integral > 4.0 * v.sigma && (v.peak.as_ref().unwrap().x - x).abs() < 0.05, "{n}: {v:?}");
    }
    // The background doubles from 4 K to 30 K at 16.5° (the mock-up's lower panel).
    assert!(w.members[18].background_at[0].unwrap() > 2.0 * w.members[0].background_at[0].unwrap());
}

/// σ of the line and of given points against finite differences: the integral and the centroid are linear (the
/// centroid to first order) in every y, the line's ends included.
#[test]
fn pl5_sigma_by_finite_differences() {
    let set = plots::open_buf(&tb().join("Tb_vs_T.buf"), None).unwrap();
    let p = &set.members[0].pattern;
    let w = Window { name: "m1".into(), from: 18.53, to: 19.77 };
    let pts = vec![(16.0, 1000.0), (25.0, 1100.0), (46.0, 1050.0)];
    for line in [true, false] {
        let b = drs::background::interpolate(&pts, drs::background::Interpolation::Spline, &p.x);
        let v0 = plots::window_value(p, &b, line, &w, &[]).unwrap();
        let (mut s_i, mut s_c) = (0.0, 0.0);
        for j in 0..p.len() {
            let mut q = p.clone();
            q.y[j] += 1e-3;
            let v = plots::window_value(&q, &b, line, &w, &[]).unwrap();
            s_i += ((v.integral - v0.integral) / 1e-3 * p.sigma[j]).powi(2);
            s_c += ((v.centroid.unwrap() - v0.centroid.unwrap()) / 1e-3 * p.sigma[j]).powi(2);
        }
        println!("line {line}: {} ± {} (finite differences {}), centroid {} ± {} ({})", v0.integral, v0.sigma, s_i.sqrt(), v0.centroid.unwrap(), v0.centroid_sigma.unwrap(), s_c.sqrt());
        assert!((v0.sigma - s_i.sqrt()).abs() < 1e-6 * v0.sigma, "line {line}");
        assert!((v0.centroid_sigma.unwrap() - s_c.sqrt()).abs() < 1e-4 * s_c.sqrt(), "line {line}");
    }
    // A Gaussian of area 100 at 10 on a sloping background, the line between the window's ends: the area and the
    // centre back; a point of no weight left out and counted.
    let x: Vec<f64> = (0..=400).map(|k| k as f64 * 0.05).collect();
    let g = |x: f64| 100.0 / (0.3 * (2.0 * std::f64::consts::PI).sqrt()) * (-(x - 10.0).powi(2) / (2.0 * 0.3 * 0.3)).exp();
    let mut p = Pattern { y: x.iter().map(|&x| 50.0 + 3.0 * x + g(x)).collect(), sigma: vec![1.0; x.len()], x };
    let w = Window { name: "g".into(), from: 7.0, to: 13.0 };
    let v = plots::window_value(&p, &[], true, &w, &[]).unwrap();
    assert!((v.integral - 100.0).abs() < 1e-6 && (v.centroid.unwrap() - 10.0).abs() < 1e-9, "{v:?}");
    assert_eq!((v.points, v.unweighted), (121, 0));
    p.sigma[200] = f64::INFINITY;
    let v = plots::window_value(&p, &[], true, &w, &[]).unwrap();
    assert_eq!((v.points, v.unweighted), (121, 1));
    // A window outside the data, for the line: refused, not extrapolated.
    assert!(plots::window_value(&p, &[], true, &Window { name: "out".into(), from: 19.0, to: 21.0 }, &[]).is_err());
}

/// The layout of a `.sum`'s phase lines: "Bragg R-factor:   6.60       Vol:  245.013( 0.006)  Fract(%):  100.00( 1.85)"
/// and "Rf-factor=  5.06", the first phase's; GLOBAL_CHI_2 is "Global user-weigthed Chi2 (Bragg contrib.):   3.86" (the
/// first "Chi2:" of the R-factors is another number: la0_122 prints 3.80 there and 3.81 here, the .seq 3.80951).
fn sum_numbers(text: &str) -> (f64, f64, f64, f64, f64, f64, f64) {
    let after = |key: &str| -> f64 {
        let i = text.find(key).unwrap() + key.len();
        text[i..].split(|c: char| c.is_whitespace() || c == '(' || c == ')').find(|s| !s.is_empty()).unwrap().parse().unwrap()
    };
    let paren = |key: &str| -> f64 {
        let i = text.find(key).unwrap() + key.len();
        let r = &text[i..];
        r[r.find('(').unwrap() + 1..r.find(')').unwrap()].trim().parse().unwrap()
    };
    (after("Global user-weigthed Chi2 (Bragg contrib.):"), after("Bragg R-factor:"), after("Rf-factor="), after("Vol:"), paren("Vol:"), after("Fract(%):"), paren("Fract(%):"))
}

/// Within half a unit of the last decimal printed.
fn printed(v: f64, s: f64, decimals: i32) -> bool {
    (v - s).abs() <= 0.5 * 10f64.powi(-decimals) + 1e-9
}

#[test]
fn pl11_the_seq_as_a_table() {
    let text = std::fs::read_to_string(d20().join("la0.seq")).unwrap();
    let seq = drs::seq::read(&text).unwrap();
    assert_eq!(seq.blocks.len(), 71);
    assert!(seq.unread.is_empty() && seq.inconsistent().is_empty() && seq.repeated().is_empty(), "{:?}", seq.unread);
    assert_eq!(seq.names().len(), 46);
    let b = &seq.blocks[0];
    assert_eq!((b.name.as_str(), b.numpar, b.patterns, b.phase_count, b.index, b.parameter, b.numor, b.chi2), ("la0120", Some(46), Some(1), Some(2), Some(0), Some(68.861), Some(73430), Some(3.86369)));
    assert_eq!(b.phases[0], drs::seq::SeqPhase { phase: 1, r_bragg: Some(6.6040678), r_f: Some(5.0635004), volume: Some((245.01285, 0.56278920e-2)), weight_fraction: Some((100.0, 1.8546307)) });
    assert_eq!(b.parameters[0], drs::seq::SeqParameter { index: 1, name: "Scale_ph1_pat1".into(), value: 2.2853527, sigma: 0.29970573e-1 });
    // Every block: PARVAL against its data file's header (the sample temperature and the numor, FullProf's Ins = 3),
    // and its numbers against its .sum to the digits the .sum prints.
    for (k, b) in seq.blocks.iter().enumerate() {
        assert_eq!(b.index, Some(k as i64));
        let scan = drs::formats::read_path(&d20().join(format!("{}.dat", b.name)), None).unwrap().scans.remove(0);
        assert_eq!(b.parameter, scan.info.temperature, "{}", b.name);
        assert_eq!(b.numor.map(|n| n.to_string()), scan.info.other.iter().find(|o| o.0 == "numor").map(|o| o.1.clone()), "{}", b.name);
        let sum = std::fs::read_to_string(d20().join(format!("la0_{}.sum", &b.name[3..]))).unwrap();
        let (chi2, rb, rf, vol, svol, wf, swf) = sum_numbers(&sum);
        let ph = &b.phases[0];
        assert!(printed(b.chi2.unwrap(), chi2, 2) && printed(ph.r_bragg.unwrap(), rb, 2) && printed(ph.r_f.unwrap(), rf, 2), "{}", b.name);
        assert!(printed(ph.volume.unwrap().0, vol, 3) && printed(ph.volume.unwrap().1, svol, 3), "{}", b.name);
        assert!(printed(ph.weight_fraction.unwrap().0, wf, 2) && printed(ph.weight_fraction.unwrap().1, swf, 2), "{}", b.name);
    }
    let a = seq.column("Cell_A_ph1_pat1");
    assert!(a.iter().all(|v| v.is_some()));
    println!("a: {:?} at {} K … {:?} at {} K", a[0], seq.blocks[0].parameter.unwrap(), a[70], seq.blocks[70].parameter.unwrap());
    // A run appended to an earlier one (news, 14 February 2008): its files found twice.
    let twice = drs::seq::read(&format!("{text}{text}")).unwrap();
    assert_eq!((twice.blocks.len(), twice.repeated().len()), (142, 71));
    // A line the reader does not know is kept, not dropped.
    let odd = drs::seq::read(&text.replacen("GLOBAL_CHI_2", "GLOBAL_CHI_3", 1)).unwrap();
    assert_eq!(odd.unread.len(), 1);
    // The command line.
    let o = drs(&["plots", "seq", d20().join("la0.seq").to_str().unwrap(), "--json"]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let back: drs::seq::Seq = serde_json::from_slice(&o.stdout).unwrap();
    assert_eq!(back, seq);
    let o = drs(&["plots", "seq", d20().join("la0.seq").to_str().unwrap(), "--parameter", "Cell_A_ph1_pat1"]);
    assert!(String::from_utf8_lossy(&o.stdout).lines().nth(1).unwrap().contains("Cell_A_ph1_pat1 5.5421996 ± 0.00007368825"));
}

/// Whether (x, y) lies in the curtain of a curve: within its x range, at or below it, at or above its floor.
fn in_curtain(c: &Displaced, floor: Option<f64>, x: f64, y: f64) -> bool {
    let n = c.x.len();
    if x < c.x[0] || x > c.x[n - 1] {
        return false;
    }
    let k = c.x.partition_point(|&v| v <= x).clamp(1, n - 1) - 1;
    let top = c.y[k] + (c.y[k + 1] - c.y[k]) * (x - c.x[k]) / (c.x[k + 1] - c.x[k]);
    y <= top && floor.is_none_or(|f| y >= f + c.shift.1)
}

#[test]
fn pl17_the_projections_hidden_parts() {
    // The mock-up F2 (plots/mockups/make_mockups.py, fig_projection): the 20 files in the buffer's order, BT nearest,
    // each 0.22° and 2200 counts beyond the one before, each hiding down to 400 counts below the lowest count of all.
    let set = plots::open_buf(&tb().join("Tb_vs_T.buf"), None).unwrap();
    let pos = |k: usize| if k == 19 { 0 } else { k + 1 }; // the set is sorted by T, BT last; the buffer has BT first
    let overrides: Vec<(usize, (f64, f64))> = (0..20).map(|k| (k, (pos(k) as f64 * 0.22, pos(k) as f64 * 2200.0))).collect();
    let d = plots::displace(&set, Displacement::None, YScale::Linear, &overrides, None);
    let order: Vec<usize> = std::iter::once(19).chain(0..19).collect();
    let floor = set.members.iter().flat_map(|m| m.pattern.y.iter().copied()).fold(f64::INFINITY, f64::min) - 400.0;
    let pr = plots::project(&d, &order, Some(floor)).unwrap();
    // The nearest curve whole, in one piece.
    let front = d[19].as_ref().unwrap();
    assert_eq!(pr.curves[0].pieces, vec![front.x.iter().copied().zip(front.y.iter().copied()).collect::<Vec<_>>()]);
    // Every other curve, sampled twenty times a segment: visible exactly where no curve in front covers it.
    let mut checked = 0;
    for (r, c) in pr.curves.iter().enumerate() {
        let me = d[c.member].as_ref().unwrap();
        for piece in &c.pieces {
            for &(x, y) in piece {
                let k = me.x.partition_point(|&v| v <= x).clamp(1, me.x.len() - 1) - 1;
                let on = me.y[k] + (me.y[k + 1] - me.y[k]) * (x - me.x[k]) / (me.x[k + 1] - me.x[k]);
                assert!((y - on).abs() < 1e-9 * on.abs().max(1.0), "a piece's point off its curve");
            }
        }
        for i in 0..me.x.len() - 1 {
            for s in 0..20 {
                let x = me.x[i] + (s as f64 + 0.5) / 20.0 * (me.x[i + 1] - me.x[i]);
                let y = me.y[i] + (s as f64 + 0.5) / 20.0 * (me.y[i + 1] - me.y[i]);
                if c.pieces.iter().any(|p| (x - p[0].0).abs() < 1e-7 || (x - p[p.len() - 1].0).abs() < 1e-7) {
                    continue;
                }
                let hidden = order[..r].iter().any(|&f| in_curtain(d[f].as_ref().unwrap(), Some(floor), x, y));
                let shown = c.pieces.iter().any(|p| x > p[0].0 && x < p[p.len() - 1].0);
                assert_eq!(shown, !hidden, "{} at {x}", set.members[c.member].label);
                checked += 1;
            }
        }
    }
    let shares: Vec<String> = pr.curves.iter().map(|c| format!("{} {:.3} ({})", set.members[c.member].label, c.visible, c.pieces.len())).collect();
    println!("{checked} samples; visible: {}", shares.join(", "));
    assert!(pr.curves[1..].iter().all(|c| c.visible < 1.0 && c.visible > 0.5));
    // The oblique axis: each file's last 2θ at the foot of its curtain.
    assert_eq!(pr.axis[1], (0, (46.0 + 0.22, floor + 2200.0)));
    // Here no curve dips below the floor of one in front, so the floor changes nothing.
    assert_eq!(plots::project(&d, &order, None).unwrap().curves, pr.curves);
    // By hand: a tent (0, 0)–(1, 2)–(2, 0) in front of the line y = 1 from −1 to 3 hides it from 0.5 to 1.5.
    let c = |x: Vec<f64>, y: Vec<f64>, dy: f64| Some(Displaced { x, y, shift: (0.0, dy), baseline: dy, background: None });
    let two = [c(vec![0.0, 1.0, 2.0], vec![0.0, 2.0, 0.0], 0.0), c(vec![-1.0, 3.0], vec![1.0, 1.0], 0.0)];
    let p = plots::project(&two, &[0, 1], None).unwrap();
    assert_eq!(p.curves[1].pieces, vec![vec![(-1.0, 1.0), (0.5, 1.0)], vec![(1.5, 1.0), (3.0, 1.0)]]);
    assert!((p.curves[1].visible - 0.75).abs() < 1e-12);
    // A curve below the front one's floor shows; without a floor it is hidden.
    let low = [c(vec![0.0, 10.0], vec![0.0, 0.0], 0.0), c(vec![0.0, 10.0], vec![-5.0, -5.0], -5.0)];
    assert_eq!(plots::project(&low, &[0, 1], Some(-1.0)).unwrap().curves[1].visible, 1.0);
    assert!(plots::project(&low, &[0, 1], None).unwrap().curves[1].pieces.is_empty());
    // Refused: a member twice, a member not placed.
    assert!(plots::project(&two, &[0, 0], None).is_err());
    assert!(plots::project(&[None, two[1].clone()], &[0, 1], None).is_err());
    // The command line, in the set's order (sorted by T, BT last).
    let o = drs(&["plots", "project", tb().join("Tb_vs_T.buf").to_str().unwrap(), "--step", "2200", "--dx", "0.22", "--floor", &floor.to_string(), "--json"]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let back: plots::Projection = serde_json::from_slice(&o.stdout).unwrap();
    let mine = plots::project(&plots::displace(&set, Displacement::Step { dx: 0.22, dy: 2200.0 }, YScale::Linear, &[], None), &(0..20).collect::<Vec<_>>(), Some(floor)).unwrap();
    assert_eq!(back.curves.len(), 20);
    for (a, b) in back.curves.iter().zip(&mine.curves) {
        assert_eq!((a.member, a.pieces.len()), (b.member, b.pieces.len()));
        assert!((a.visible - b.visible).abs() < 1e-12);
    }
}

#[test]
fn pl5_the_command_line() {
    let buf = tb().join("Tb_vs_T.buf");
    let o = drs(&["plots", "windows", buf.to_str().unwrap(), "--window", "m1:18.5:19.8", "--window", "m3:31.9:33.3", "--at", "16.5,30", "--json"]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let back: plots::Windows = serde_json::from_slice(&o.stdout).unwrap();
    let set = plots::open_buf(&buf, None).unwrap();
    let mine = plots::windows(&set, &[f5_windows()[0].clone(), f5_windows()[2].clone()], &WindowBackground::Auto { log_lambda: None }, &[16.5, 30.0]).unwrap();
    assert_eq!(back, mine);
    // Given points, from a background file; the line.
    let dir = root("target/test-plots-series");
    std::fs::create_dir_all(&dir).unwrap();
    let points = dir.join("tb.drb");
    std::fs::write(&points, "format = \"diff-ref-struct background 1\"\ninterpolation = \"linear\"\npoints = [[16.0, 1000.0], [46.0, 1000.0]]\n").unwrap();
    let o = drs(&["plots", "windows", buf.to_str().unwrap(), "--window", "m1:18.5:19.8", "--background", points.to_str().unwrap(), "--at", "30"]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(String::from_utf8_lossy(&o.stdout).contains("b(30) 1000.000"));
    let o = drs(&["plots", "windows", buf.to_str().unwrap(), "--window", "m1:18.5:19.8", "--background", "line"]);
    assert!(o.status.success() && String::from_utf8_lossy(&o.stdout).starts_with("# windows above the line"));
    assert!(!drs(&["plots", "windows", buf.to_str().unwrap()]).status.success(), "no window given");
}

/// (a) FullProf reads Ins = 3's fourth line as xm, ti, x0, ome, chi, phi, tr1, tr2, step, tset, tcons, tsamp (manual
/// p. 133, `manuals/fullprof/FullProf_Manual.pdf`): no wavelength. The 1.870 of la0120.dat is in `tr2`; D1B's hocu.dat
/// has 0.000 there while its .pcr gives λ = 2.524 Å. The title's "1.87AA" is free text. So the reader states none.
#[test]
fn finding_a_the_d20_reader_states_no_wavelength() {
    let field = |path: &Path| -> f64 { std::fs::read_to_string(path).unwrap().lines().nth(3).unwrap().split_whitespace().nth(7).unwrap().parse().unwrap() };
    let fp = root("../manuals/fullprof");
    assert_eq!(field(&d20().join("la0120.dat")), 1.870);
    assert_eq!(field(&fp.join("FP_tests/hocu.dat")), 0.0);
    let pcr = std::fs::read_to_string(fp.join("FP_tests/hocu.pcr")).unwrap();
    let lambda: f64 = pcr.lines().skip_while(|l| !l.starts_with("! Lambda1")).nth(1).unwrap().split_whitespace().next().unwrap().parse().unwrap();
    assert_eq!(lambda, 2.524);
    for f in [d20().join("la0120.dat"), fp.join("FP_tests/hocu.dat")] {
        let d = drs::formats::read_path(&f, Some("ill_d1b_d20")).unwrap();
        assert_eq!(d.scans[0].wavelength, None, "{}", f.display());
    }
}

/// (b) A point of no weight (σ infinite) is `null` in JSON and reads back as no weight, through the core's own types.
#[test]
fn finding_b_no_weight_round_trips_through_json() {
    let p = Pattern { x: vec![1.0, 2.0, 3.0], y: vec![5.0, 0.0, 7.0], sigma: vec![2.0, f64::INFINITY, 3.0] };
    let text = serde_json::to_string(&p).unwrap();
    assert_eq!(text, r#"{"x":[1.0,2.0,3.0],"y":[5.0,0.0,7.0],"sigma":[2.0,null,3.0]}"#);
    assert_eq!(serde_json::from_str::<Pattern>(&text).unwrap(), p);
    // A set holding an xye file with a zero count (no weight, GSAS-II's rule), through `drs plots open --json`.
    let dir = root("target/test-plots-series/nullsigma");
    std::fs::create_dir_all(&dir).unwrap();
    let f = dir.join("zero.xye");
    std::fs::write(&f, "# a zero count\n10.0 4.0 2.0\n10.1 0.0 1.0\n10.2 9.0 3.0\n10.3 16.0 4.0\n").unwrap();
    let o = drs(&["plots", "open", f.to_str().unwrap(), "--json"]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(String::from_utf8_lossy(&o.stdout).contains("null"));
    let set: plots::PlotSet = serde_json::from_slice(&o.stdout).unwrap();
    let direct = drs::formats::read_path(&f, None).unwrap().scans.remove(0).pattern;
    assert!(direct.sigma[1].is_infinite());
    assert_eq!(set.members[0].pattern, direct);
}
