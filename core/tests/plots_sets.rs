//! Plots' PL12, PL18, PL19 and PL20 and Explore's E13a and E13b (L192):
//! - PL12: the plot set file `Plots/{set}.drv` read, written and read back to the same set; its members opened relative
//!   to the project's folder, each typed value with its source, the order; what is wrong said with its place; `drs plots
//!   save`, `drs plots open <set.drv>`, `drs plots export --csv`;
//! - PL18: the map's rows for the members apart, less the reference in a difference map;
//! - PL19: a step proposed for a stack, with its rule;
//! - PL20: members drawn after Data's operations;
//! - E13a: `Grid::map` against the stand-in set the app built (`app/src-tauri/src/modules/explore.rs`, `map_of`), to the
//!   bit, on a 5 × 5 grid of numbers and on a 5 × 5 χ² grid of Y₂O₃;
//! - E13b: `Ellipse::outline`, the window's parametrisation (`draw.js`, `ellipseOutline`), on its ellipse;
//! - PL23 (L202): a value read from a temperatures file kept in the set file with that source, through Save and Open,
//!   in the core and on the command line;
//! - PL24 (L202): a `.drv` known by its content wherever it lies and whatever its name, `drs data info` included;
//! - PL37 (L405, L406, L407): the map's clip (`map_range`, `difference_range`) and the regions saved by name kept by
//!   the set file, its refusals, and `plots::map_region`'s block of the Tb₂Pd₂Sn map, a region in 2θ and the same in d.
//! - PL41 (W4): the y interval kept by the set file (`y_interval = "own"`, absent shared) and `plots::extent`, the
//!   set's counts on a scale, equal to the map's colour range on the log scale.
//! - PL42 (L434): the x interval of the analysis kept by the set file (`x_interval = [a, b]`) and applied on opening,
//!   every member cut to it: on the ESRF run's names a map spanning exactly 2–20° and a window outside refused, a
//!   member that does not reach it refused by name; the Tb₂Pd₂Sn set cut to 18–34° giving the whole map's block of
//!   those columns (`map_region`), `drs plots map --x-interval` the same.

use drs::data::Pattern;
use drs::data_ops::{Normalization, Operation};
use drs::formats::Axis;
use drs::lsq::{covariance_ellipse, Grid};
use drs::plots::{self, Found, Order, SetFile, YScale};
use std::path::{Path, PathBuf};

fn root(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)
}

fn tb() -> PathBuf {
    root("../manuals/fullprof/1-WinPLOTR/1_Tb2pd2sn")
}

fn drs(args: &[&str]) -> std::process::Output {
    std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(args).output().unwrap()
}

/// A project with the Tb₂Pd₂Sn files in Data/ and a set file in Plots/ (plots/DESIGN.md §2.2's example).
fn project(name: &str) -> PathBuf {
    let d = root("target/test-plots-sets").join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(d.join("Data")).unwrap();
    std::fs::create_dir_all(d.join("Plots")).unwrap();
    for k in 1..=20 {
        let f = format!("tb_t{k:02}.dat");
        std::fs::copy(tb().join(&f), d.join("Data").join(&f)).unwrap();
    }
    d
}

const EXAMPLE: &str = r#"format = "drs plot set 1"
title = "Tb2Pd2.05Sn0.95, G41 (LLB), 1.5–30 K"
variable = { name = "T", unit = "K" }
x = "2theta"                 # 2theta | d | q | stl | tof | energy
y_scale = "sqrt"
order = "variable"
displacement = { rule = "step", dy = 3000.0 }
reference = "Data/tb_t20.dat"
windows = [ { name = "m1", from = 18.5, to = 19.8 }, { name = "m3", from = 31.9, to = 33.3 } ]

[[member]]
file = "Data/tb_t20.dat"

[[member]]
file = "Data/tb_t01.dat"
variable = { value = 1.5, source = "WinPLOTR demo 4 gives 1.5–30 K; the file says BT" }

[[member]]
file = "Data/tb_t02.dat"     # variable read from the title "T.=.4.K"
opacity = 0.3

[[member]]
file = "Data/tb_t12.dat"
colour = "accent"
visible = false
"#;

#[test]
fn pl12_the_plot_set_file() {
    let d = project("pl12");
    let f = SetFile::read(EXAMPLE).unwrap();
    assert_eq!((f.y_scale, f.order, f.windows.len(), f.members.len()), (YScale::Sqrt, Order::Variable, 2, 4));
    assert_eq!(f.displacement, plots::Displacement::Step { dx: 0.0, dy: 3000.0 });
    assert_eq!((f.members[2].opacity, f.members[3].visible, f.members[3].colour.as_deref()), (Some(0.3), Some(false), Some("accent")));
    // Written and read back: the same file.
    assert_eq!(SetFile::read(&f.to_text()).unwrap(), f);
    // Opened from Plots/: the files relative to the folder above; BT's typed value with its source; by T, so BT first.
    let path = d.join("Plots/tb.drv");
    std::fs::write(&path, EXAMPLE).unwrap();
    assert_eq!(plots::set_base(&path), d.join("Plots").parent().unwrap());
    let set = f.open(&plots::set_base(&path), None, &|_| Vec::new()).unwrap();
    let labels: Vec<&str> = set.members.iter().map(|m| m.label.as_str()).collect();
    assert_eq!(labels, ["tb_t01.dat", "tb_t02.dat", "tb_t12.dat", "tb_t20.dat"]);
    assert_eq!(set.members[0].variable.as_ref().map(|v| (v.value, v.found.clone())), Some((1.5, Found::Typed("WinPLOTR demo 4 gives 1.5–30 K; the file says BT".into()))));
    assert_eq!(set.members[1].variable.as_ref().map(|v| v.value), Some(4.0));
    assert_eq!(f.reference_in(&set), Some(3));
    let given = SetFile { order: Order::Given, ..f.clone() }.open(&plots::set_base(&path), None, &|_| Vec::new()).unwrap();
    assert_eq!(given.members[0].label, "tb_t20.dat");
    // What is wrong, with its place.
    let wrong = |text: &str| SetFile::read(text).unwrap_err().iter().map(|x| x.message.clone()).collect::<Vec<_>>().join("; ");
    assert!(wrong("title = \"x\"\n[[member]]\nfile = \"a.dat\"\n").contains("format"));
    assert!(wrong(&EXAMPLE.replace("opacity = 0.3", "opacity = 3.0")).contains("opacity runs from 0 to 1"));
    assert!(wrong(&EXAMPLE.replace(", source = \"WinPLOTR demo 4 gives 1.5–30 K; the file says BT\"", "")).contains("source"), "a typed value without its source");
    assert!(wrong(&EXAMPLE.replace("colour = \"accent\"", "color = \"accent\"")).contains("\"color\" is not a key of [[member]]"));
    assert!(wrong(&EXAMPLE.replace("[[member]]", "[[members]]")).contains("not a table of a plot set"));
    let e = SetFile::read(&EXAMPLE.replace("opacity = 0.3", "opacity = 3.0")).unwrap_err();
    assert_eq!(&EXAMPLE.replace("opacity = 0.3", "opacity = 3.0")[e[0].span.start..e[0].span.end], "3.0");
    // x is the view's quantity (PL10, PL27): the members open as their files give them, the set file keeps "d".
    let in_d = SetFile::read(&EXAMPLE.replace("x = \"2theta\"", "x = \"d\"")).unwrap();
    assert_eq!(in_d.x, "d");
    assert_eq!(in_d.open(&d, None, &|_| Vec::new()).unwrap().members.len(), 4);
    // A set written from a set: the typed value kept, read back the same members.
    let back = SetFile::read(&SetFile::of(&set).to_text()).unwrap().open(&plots::set_base(&path), None, &|_| Vec::new()).unwrap();
    assert_eq!(back.members, set.members);
}

#[test]
fn pl12_the_command_line() {
    let d = project("cli");
    // A buffer saved as a set file in Plots/: its files named from the project's folder.
    std::fs::copy(tb().join("Tb_vs_T.buf"), d.join("Data/Tb_vs_T.buf")).unwrap();
    let out = d.join("Plots/Tb.drv");
    let o = drs(&["plots", "save", d.join("Data/Tb_vs_T.buf").to_str().unwrap(), "--out", out.to_str().unwrap()]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let text = std::fs::read_to_string(&out).unwrap();
    assert!(text.starts_with("format = \"drs plot set 1\"") && text.contains("file = \"Data/tb_t01.dat\""), "{text}");
    assert!(!drs(&["plots", "save", d.join("Data/Tb_vs_T.buf").to_str().unwrap(), "--out", out.to_str().unwrap()]).status.success(), "never over a file");
    // Opened as the buffer opens.
    let json = |args: &[&str]| -> plots::PlotSet { serde_json::from_slice(&drs(args).stdout).unwrap() };
    let a = json(&["plots", "open", out.to_str().unwrap(), "--json"]);
    let b = json(&["plots", "open", d.join("Data/Tb_vs_T.buf").to_str().unwrap(), "--json"]);
    assert_eq!(a.members.len(), 20);
    for (x, y) in a.members.iter().zip(&b.members) {
        assert_eq!((&x.label, &x.variable, &x.pattern), (&y.label, &y.variable, &y.pattern));
    }
    // The CSV: 2θ, then y and σ of each member; each column said in a # line; every value as the member's.
    let o = drs(&["plots", "export", out.to_str().unwrap(), "--csv"]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let csv = String::from_utf8_lossy(&o.stdout).to_string();
    let comments: Vec<&str> = csv.lines().filter(|l| l.starts_with('#')).collect();
    assert_eq!(comments.len(), 21);
    assert!(comments[1].starts_with("# y1: Data/tb_t02.dat") && comments[1].contains("T = 4 K (the title's \"T.=.4.K\")"), "{}", comments[1]);
    let rows: Vec<Vec<&str>> = csv.lines().filter(|l| !l.starts_with('#')).map(|l| l.split(',').collect()).collect();
    assert_eq!((rows.len(), rows[0].len(), rows[0][0], rows[0][1]), (302, 41, "2theta_deg", "y1"));
    for (i, r) in rows[1..].iter().enumerate() {
        for (k, m) in a.members.iter().enumerate() {
            assert_eq!((r[1 + 2 * k].parse::<f64>().unwrap(), r[2 + 2 * k].parse::<f64>().unwrap()), (m.pattern.y[i], m.pattern.sigma[i]));
        }
    }
    // The file's windows, reference and scale serve the commands that take them.
    std::fs::write(d.join("Plots/ex.drv"), EXAMPLE).unwrap();
    let o = drs(&["plots", "windows", d.join("Plots/ex.drv").to_str().unwrap()]);
    assert!(o.status.success() && String::from_utf8_lossy(&o.stdout).contains("m3 "), "{}", String::from_utf8_lossy(&o.stderr));
    let m: plots::Map = serde_json::from_slice(&drs(&["plots", "map", d.join("Plots/ex.drv").to_str().unwrap(), "--json"]).stdout).unwrap();
    assert_eq!((m.reference, m.scale), (Some(3), YScale::Sqrt));
}

/// PL18: the rows apart, on the grid and the scale, less the reference.
#[test]
fn pl18_the_rows_apart() {
    let set = plots::open_buf(&tb().join("Tb_vs_T.buf"), None).unwrap();
    let grid = set.members[0].pattern.x.clone();
    let bt = 19;
    assert_eq!(set.members[bt].label, "tb_t01.dat");
    let m = plots::map(&set, &grid, YScale::Sqrt, None).unwrap();
    assert_eq!((m.apart.clone(), m.apart_z.len()), (vec![bt], 1));
    assert!(m.apart_z[0].iter().zip(&set.members[bt].pattern.y).all(|(z, y)| *z == Some(YScale::Sqrt.apply(*y))));
    let k30 = 18;
    let dm = plots::map(&set, &grid, YScale::Linear, Some(k30)).unwrap();
    let (bty, r) = (&set.members[bt].pattern.y, &set.members[k30].pattern.y);
    assert!(dm.apart_z[0].iter().enumerate().all(|(j, z)| *z == Some(bty[j] - r[j])));
    let largest = dm.z.iter().chain(&dm.apart_z).flatten().flatten().fold(0.0_f64, |a, v| a.max(v.abs()));
    assert_eq!(dm.range, (-largest, largest));
}

/// PL19: the step proposed for a stack, with its rule.
#[test]
fn pl19_a_proposed_step() {
    let set = plots::open_buf(&tb().join("Tb_vs_T.buf"), None).unwrap();
    let p = plots::proposed_step(&set, YScale::Linear).unwrap();
    // By hand: each member's 90th less 10th percentile (nearest rank), their median, then 1, 2 or 5 × 10^k above it.
    let mut spreads: Vec<f64> = set
        .members
        .iter()
        .map(|m| {
            let mut v = m.pattern.y.clone();
            v.sort_by(f64::total_cmp);
            let q = |f: f64| v[(f * (v.len() - 1) as f64).round() as usize];
            q(0.9) - q(0.1)
        })
        .collect();
    assert_eq!(p.spreads, spreads.iter().map(|s| Some(*s)).collect::<Vec<_>>());
    spreads.sort_by(f64::total_cmp);
    let median = (spreads[9] + spreads[10]) / 2.0;
    println!("Tb: median spread {median}, step {} ({})", p.dy, p.rule);
    assert!(p.dy >= median && p.dy < 2.5 * median && [1.0, 2.0, 5.0].iter().any(|m| (p.dy / m / 10f64.powf(p.dy.log10().floor()) - 1.0).abs() < 1e-12));
    assert!(p.rule.contains("median over 20 members"));
    let s = plots::proposed_step(&set, YScale::Sqrt).unwrap();
    assert!(s.dy < p.dy && s.scale == YScale::Sqrt);
    // Each member's bulk stays clear of the next one's: the step at least each member's 10th to 90th percentile for
    // half of them, by the rule.
    assert!(p.spreads.iter().flatten().filter(|&&x| x <= p.dy).count() >= 10);
}

/// PL20: a member drawn after Data's operations, as `apply_all` gives them.
#[test]
fn pl20_members_after_datas_operations() {
    let ops = |f: &str| -> Vec<Operation> {
        match f {
            "tb_t02.dat" => vec![Operation::Delete { from: 16.0, to: 18.0 }, Operation::Normalize(Normalization::Maximum(1.0))],
            _ => Vec::new(),
        }
    };
    let plain = plots::open_buf(&tb().join("Tb_vs_T.buf"), None).unwrap();
    let set = plots::open_buf_with(&tb().join("Tb_vs_T.buf"), None, &ops).unwrap();
    let k = set.members.iter().position(|m| m.file == "tb_t02.dat").unwrap();
    let scan = drs::formats::read_path(&tb().join("tb_t02.dat"), None).unwrap().scans.remove(0);
    let expect = drs::data_ops::apply_all(&scan, &ops("tb_t02.dat"), &|_, _| Err("none".into())).unwrap();
    assert_eq!(set.members[k].pattern, expect.pattern);
    assert_eq!((set.members[k].facts.points, set.members[k].facts.range.0), (280, 18.1));
    assert_eq!(set.members[k].pattern.max_count(), Some(1.0));
    assert_eq!(set.members[k].variable, plain.members[k].variable, "the header's temperature as read");
    for (a, b) in set.members.iter().zip(&plain.members).filter(|(a, _)| a.file != "tb_t02.dat") {
        assert_eq!(a.pattern, b.pattern);
    }
    // An operation that names another file reads it relative to the set's folder.
    let sub = |f: &str| -> Vec<Operation> { if f == "tb_t02.dat" { vec![Operation::Subtract { file: "tb_t20.dat".into(), reader: None, scale: 1.0 }] } else { Vec::new() } };
    let d = plots::open_buf_with(&tb().join("Tb_vs_T.buf"), None, &sub).unwrap();
    let (a, b) = (&plain.members[k].pattern, &plain.members.iter().find(|m| m.file == "tb_t20.dat").unwrap().pattern);
    assert!(d.members[k].pattern.y.iter().enumerate().all(|(i, y)| *y == a.y[i] - b.y[i]));
}

/// The app's stand-in (`explore.rs`, `map_of`, app 0.17.0): one member per value of b, its scan χ² along a.
fn stand_in(g: &Grid, scale: YScale) -> plots::Map {
    let (first, last) = (g.xa[0], g.xa[g.xa.len() - 1]);
    let members = g
        .xb
        .iter()
        .enumerate()
        .map(|(j, &vb)| plots::Member {
            file: g.b.clone(),
            scan: j,
            scans: g.xb.len(),
            label: g.b.clone(),
            comment: None,
            facts: plots::Facts { format: "chi2_grid".into(), axis: Axis::TwoTheta, points: g.xa.len(), range: (first, last), step: None, wavelength: None, title: None, temperature: None, monitor: None, counting_time: None },
            variable: Some(plots::Variable { value: vb, found: plots::Found::Typed("the χ² grid's second parameter".into()) }),
            candidates: Vec::new(),
            pattern: Pattern { x: g.xa.clone(), y: g.chi2.iter().map(|row| row[j]).collect(), sigma: vec![1.0; g.xa.len()] },
            others: Vec::new(),
            modified: None,
        })
        .collect();
    let set = plots::PlotSet { title: None, variable_name: g.b.clone(), variable_unit: String::new(), members, x_interval: None, x_interval_in: None };
    plots::map(&set, &g.xa, scale, None).unwrap()
}

#[test]
fn e13a_the_grid_as_a_map() {
    // A 5 × 5 grid of numbers with a point not computed and one ≤ 0.
    let mut chi2: Vec<Vec<f64>> = (0..5).map(|i| (0..5).map(|j| 1.0 + (i as f64 - 2.0).powi(2) + 0.5 * (j as f64 - 2.0).powi(2)).collect()).collect();
    chi2[0][4] = f64::NAN;
    chi2[4][0] = 0.0;
    let g = Grid { a: "phase1.caglioti_v".into(), b: "phase1.caglioti_w".into(), xa: (0..5).map(|i| -0.4 + 0.01 * i as f64).collect(), xb: (0..5).map(|j| 0.3 + 0.02 * j as f64).collect(), chi2, order: Vec::new(), computed: 24, dof: 0 };
    for scale in [YScale::Log, YScale::Linear, YScale::Sqrt] {
        assert_eq!(g.map(scale).unwrap(), stand_in(&g, scale), "{scale:?}");
    }
    let m = g.map(YScale::Log).unwrap();
    // Not computed (NaN) and ≤ 0 on log: no value; the range over the rest.
    assert_eq!((m.z[0][4], m.z[4][0]), (None, None));
    let computed: Vec<f64> = g.chi2.iter().flatten().copied().filter(|v| *v > 0.0).collect();
    assert_eq!(m.range, (computed.iter().copied().fold(f64::INFINITY, f64::min).log10(), computed.iter().copied().fold(f64::NEG_INFINITY, f64::max).log10()));
    assert_eq!((m.x_edges[0], m.rows[0].edges.0), (g.xa[0] - (g.xa[1] - g.xa[0]) / 2.0, g.xb[0] - (g.xb[1] - g.xb[0]) / 2.0));
    // Y₂O₃'s χ² on a 5 × 5 grid of V and W about their values, as Explore computes it.
    let d = root("../manuals/fullprof/runs/Y2O3");
    let text: String = std::fs::read(d.join("y2o3_structure.pcr")).unwrap().iter().map(|&b| b as char).collect();
    let p = drs::pcr::read(&text).unwrap();
    let loaded = drs::powder::Model::from_pcr(&p, &d, "y2o3_structure").unwrap();
    let set = drs::refinable::ParameterSet::from_pcr(&p, &loaded.model).unwrap();
    let j = |n: &str| set.free.iter().position(|q| q.name == n).unwrap();
    let grid = drs::lsq::chi2_grid(&loaded.model, loaded.data.as_ref().unwrap(), &set, (j("phase1.caglioti_v"), j("phase1.caglioti_w")), (0.01, 0.01), 5, |_, _| true);
    for scale in [YScale::Log, YScale::Linear] {
        assert_eq!(grid.map(scale).unwrap(), stand_in(&grid, scale), "{scale:?}, to the bit");
    }
    let lin = grid.map(YScale::Linear).unwrap();
    assert!(lin.z.iter().enumerate().all(|(r, row)| row.iter().enumerate().all(|(i, v)| *v == Some(grid.chi2[i][lin.rows[r].member]))));
    println!("Y2O3 V–W: χ² from {} to {}", lin.range.0, lin.range.1);
}

#[test]
fn e13b_the_ellipses_outline() {
    let (sa, sb, r) = (0.0064, 0.0131, -0.871);
    let e = covariance_ellipse(sa, sb, r).unwrap();
    let centre = (0.17, -0.39);
    for k in [1.0, 2.0] {
        let o = e.outline(centre, k, 90);
        assert_eq!(o.len(), 91);
        assert_eq!(o[0], o[90]);
        // The window's parametrisation (draw.js, ellipseOutline), the same operations in the same order.
        for (i, p) in o.iter().take(90).enumerate() {
            let t = 2.0 * std::f64::consts::PI * i as f64 / 90.0;
            let (u, v) = (k * e.l1.sqrt() * t.cos(), k * e.l2.sqrt() * t.sin());
            let q = (centre.0 + u * e.angle.cos() - v * e.angle.sin(), centre.1 + u * e.angle.sin() + v * e.angle.cos());
            assert!((p.0 - q.0).abs() <= 1e-15 && (p.1 - q.1).abs() <= 1e-15, "{i}");
            // On the k σ contour: (x − c)ᵀ C⁻¹ (x − c) = k².
            let (dx, dy) = (p.0 - centre.0, p.1 - centre.1);
            let det = sa * sa * sb * sb * (1.0 - r * r);
            let m = (sb * sb * dx * dx - 2.0 * r * sa * sb * dx * dy + sa * sa * dy * dy) / det;
            assert!((m - k * k).abs() < 1e-9 * k * k, "{i}: {m}");
        }
    }
}

#[test]
fn pl23_a_temperatures_files_value_kept() {
    // Plots' self-check (plots/REQUESTS.md PL23): tb_t03.dat takes 6.5 K from Tb.tem, the set is saved and opened again.
    let d = project("pl23");
    std::fs::write(d.join("Data/Tb.tem"), "tb_t03.dat TSAMP 6.5\n").unwrap();
    let t = plots::read_temperatures("tb_t03.dat TSAMP 6.5\n");
    let mut set = plots::open(&d, &[("Data/tb_t02.dat".into(), None), ("Data/tb_t03.dat".into(), None)], None, None).unwrap();
    let t03 = |s: &plots::PlotSet| s.members.iter().find(|m| m.file == "Data/tb_t03.dat").unwrap().clone();
    assert_eq!(t03(&set).variable.as_ref().map(|v| v.value), Some(6.0), "the title's T.=.6.K");
    assert_eq!(set.variable_from_file("Tb.tem", &t, false), 1);
    let f = SetFile::of(&set);
    let text = f.to_text();
    assert!(text.contains("variable = { value = 6.5, temperatures = \"Tb.tem\" }"), "{text}");
    let back = SetFile::read(&text).unwrap();
    assert_eq!(back, f);
    let opened = back.open(&d, None, &|_| Vec::new()).unwrap();
    let m = t03(&opened);
    assert_eq!(m.variable.as_ref().map(|v| (v.value, v.found.clone())), Some((6.5, Found::TemperatureFile("Tb.tem".into()))), "not back at the title's 6 K");
    assert_eq!(m.variable.as_ref().unwrap().found.words(), "the temperatures file Tb.tem");
    assert!(m.candidates.iter().any(|c| matches!(c.found, Found::Title(_)) && c.value == 6.0), "the title's value still a candidate");
    assert!(m.candidates.iter().any(|c| c.found == Found::TemperatureFile("Tb.tem".into()) && c.value == 6.5));
    // tb_t02.dat's value is its title's: not written, read again from the file.
    assert_eq!(f.members.iter().find(|x| x.file == "Data/tb_t02.dat").unwrap().value, None);
    // A value has one source, never none.
    let wrong = |t: &str| SetFile::read(t).unwrap_err().iter().map(|x| x.message.clone()).collect::<Vec<_>>().join("; ");
    assert!(wrong(&text.replace("temperatures = \"Tb.tem\"", "temperatures = \"Tb.tem\", source = \"typed\"")).contains("not both"));
    assert!(wrong(&text.replace(", temperatures = \"Tb.tem\"", "")).contains("never a value without one"));
    // The command line: --temperatures on a set file, saved, opened again.
    std::fs::write(d.join("Plots/two.drv"), "format = \"drs plot set 1\"\n\n[[member]]\nfile = \"Data/tb_t02.dat\"\n\n[[member]]\nfile = \"Data/tb_t03.dat\"\n").unwrap();
    let out = d.join("Plots/two again.drv");
    let o = drs(&["plots", "save", d.join("Plots/two.drv").to_str().unwrap(), "--temperatures", d.join("Data/Tb.tem").to_str().unwrap(), "--out", out.to_str().unwrap()]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(std::fs::read_to_string(&out).unwrap().contains("variable = { value = 6.5, temperatures = \"Tb.tem\" }"));
    let o = drs(&["plots", "open", out.to_str().unwrap(), "--json"]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let s: plots::PlotSet = serde_json::from_slice(&o.stdout).unwrap();
    assert_eq!(t03(&s).variable.map(|v| (v.value, v.found)), Some((6.5, Found::TemperatureFile("Tb.tem".into()))));
}

#[test]
fn pl24_a_set_file_by_its_content() {
    use drs::formats::project::{describe, describe_path, Description};
    use drs::formats::Probe;
    let d = project("pl24");
    // Outside Plots/, under another name and another extension: its format line names it.
    let text = EXAMPLE.to_string();
    let away = d.join("from a colleague.txt");
    std::fs::write(&away, &text).unwrap();
    let bytes = std::fs::read(&away).unwrap();
    let name = away.file_name().unwrap().to_string_lossy().into_owned();
    let p = Probe::new(&name, &bytes);
    assert!(matches!(describe(&p), Description::Project(ref k) if k.id == "plot_set"));
    assert!(matches!(describe_path(&away, &p), Description::Project(ref k) if k.id == "plot_set"));
    // A .drv by its extension too, when its first lines do not say (a set file being written).
    let p = Probe::new("set.drv", b"title = \"no format line yet\"\n");
    assert!(matches!(describe(&p), Description::Project(ref k) if k.id == "plot_set"));
    // The control file and the background points are still theirs.
    let drc = std::fs::read(root("../explore/compare/cases/y2o3_a.drc")).unwrap();
    assert!(matches!(describe(&Probe::new("y2o3_a.drc", &drc)), Description::Project(ref k) if k.id == "control"));
    // drs data info: no longer "unknown".
    let o = drs(&["data", "info", away.to_str().unwrap()]);
    let out = String::from_utf8_lossy(&o.stdout).to_string() + &String::from_utf8_lossy(&o.stderr);
    assert!(out.contains("plot set") && !out.contains("unknown"), "{out}");
}

#[test]
fn pl27_each_views_steps_and_scale() {
    let d = project("pl27");
    // The stack's step proposed on the linear scale, the projection's steps on the log scale, the curves on √, the map
    // on log, the members the same area over 20–30°, shown in d.
    let text = EXAMPLE.replace("displacement = { rule = \"step\", dy = 3000.0 }", "displacement = { rule = \"step\", dy = 3000.0, scale = \"linear\" }\nprojection = { rule = \"step\", dx = 0.2, dy = 0.15, scale = \"log\" }\nmap_scale = \"log\"\nnormalisation = { by = \"area\", from = 20.0, until = 30.0, to = 1.0 }").replace("x = \"2theta\"", "x = \"d\"");
    let f = SetFile::read(&text).unwrap();
    assert_eq!((f.displacement, f.displacement_scale), (plots::Displacement::Step { dx: 0.0, dy: 3000.0 }, Some(YScale::Linear)));
    assert_eq!((f.projection, f.projection_scale), (plots::Displacement::Step { dx: 0.2, dy: 0.15 }, Some(YScale::Log)));
    assert_eq!((f.y_scale, f.map_scale), (YScale::Sqrt, Some(YScale::Log)));
    assert_eq!(f.normalisation, plots::Normalisation::Area { from: 20.0, until: 30.0, to: 1.0 });
    assert_eq!(f.x, "d");
    // Written and read back, the same; the defaults write none of the new keys.
    assert_eq!(SetFile::read(&f.to_text()).unwrap(), f);
    let plain = SetFile::read(EXAMPLE).unwrap().to_text();
    assert!(!plain.contains("projection") && !plain.contains("map_scale") && !plain.contains("normalisation"), "{plain}");
    // Every normalisation goes and comes back.
    for n in [
        plots::Normalisation::Monitor { to: 1000.0 },
        plots::Normalisation::Time { to: 1.0 },
        plots::Normalisation::Maximum { to: 100.0 },
        plots::Normalisation::MaximumIn { from: 18.0, until: 20.0, to: 1.0 },
        plots::Normalisation::Peak { window: plots::Window { name: "m1".into(), from: 18.5, to: 19.8 }, background: plots::WindowBackground::Line, to: 1.0 },
        plots::Normalisation::Peak { window: plots::Window { name: "m1".into(), from: 18.5, to: 19.8 }, background: plots::WindowBackground::Auto { log_lambda: Some(5.0) }, to: 1.0 },
        plots::Normalisation::Peak { window: plots::Window { name: "m1".into(), from: 18.5, to: 19.8 }, background: plots::WindowBackground::Points { points: vec![(18.0, 900.0), (20.0, 950.0)], interpolation: drs::background::Interpolation::Linear }, to: 1.0 },
        plots::Normalisation::Reference { member: 2, from: 20.0, until: 30.0 },
    ] {
        let g = SetFile { normalisation: n.clone(), ..f.clone() };
        assert_eq!(SetFile::read(&g.to_text()).unwrap().normalisation, n);
    }
    let wrong = |t: &str| SetFile::read(t).unwrap_err().iter().map(|x| x.message.clone()).collect::<Vec<_>>().join("; ");
    assert!(wrong(&text.replace("map_scale = \"log\"", "map_scale = \"ln\"")).contains("map_scale is"));
    assert!(wrong(&text.replace("by = \"area\", from = 20.0, ", "by = \"area\", ")).contains("needs from"));
    // The command line: the set file's normalisation is the view's when --by is not given, and save writes --by's.
    let drv = d.join("Plots/pl27.drv");
    std::fs::write(&drv, &text).unwrap();
    let o = drs(&["plots", "normalise", drv.to_str().unwrap(), "--json"]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let n: Vec<plots::Normalised> = serde_json::from_slice(&o.stdout).unwrap();
    let opened = f.open(&d, None, &|_| Vec::new()).unwrap();
    assert_eq!(n, plots::normalise(&opened, &f.normalisation).unwrap());
    let out = d.join("Plots/pl27 saved.drv");
    let o = drs(&["plots", "save", drv.to_str().unwrap(), "--by", "maximum", "--to", "100", "--out", out.to_str().unwrap()]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let saved = SetFile::read(&std::fs::read_to_string(&out).unwrap()).unwrap();
    assert_eq!(saved.normalisation, plots::Normalisation::Maximum { to: 100.0 });
    assert_eq!((saved.projection, saved.map_scale), (f.projection, f.map_scale), "the rest kept");
}

/// PL37: the clip and the regions, kept, refused when wrong, and a region's block of the map.
#[test]
fn pl37_the_maps_clip_and_its_regions() {
    use plots::{MapX, Range, Region, XQuantity};
    let d = project("pl37");
    let regions = "map_range = { from = 1200.0, to = 20000.0 }\ndifference_range = { from = -5000.0, to = 5000.0 }\n\n[[regions]]\nname = \"m1 transition\"\nx = [18.95, 19.25]\nx_quantity = \"2theta\"\nvariable = [20.0, 26.0]\nmap_range = { from = 2000.0, to = 70000.0 }\nwindow = \"m1\"\n\n[[regions]]\nname = \"m3\"\nx = [31.9, 33.3]\nvariable = [1.5, 30.0]\n";
    let text = EXAMPLE.replace("\n[[member]]\nfile = \"Data/tb_t20.dat\"", &format!("{regions}\n[[member]]\nfile = \"Data/tb_t20.dat\""));
    let f = SetFile::read(&text).unwrap();
    assert_eq!((f.map_range, f.difference_range), (Some(Range { from: 1200.0, to: 20000.0 }), Some(Range { from: -5000.0, to: 5000.0 })));
    assert_eq!(f.regions.len(), 2);
    assert_eq!(f.regions[0], Region { name: "m1 transition".into(), x: [18.95, 19.25], x_quantity: Some(XQuantity::TwoTheta), variable: [20.0, 26.0], map_range: Some(Range { from: 2000.0, to: 70000.0 }), difference_range: None, window: Some("m1".into()) });
    assert_eq!((f.regions[1].x_quantity, f.regions[1].map_range, &f.regions[1].window), (None, None, &None));
    assert_eq!(f.members.len(), 4, "the members after the regions all read");
    // Written and read back: the same.
    assert_eq!(SetFile::read(&f.to_text()).unwrap(), f);
    // A .drv written before PL37 reads as before and writes none of the new keys.
    let old = SetFile::read(EXAMPLE).unwrap();
    assert_eq!((old.map_range, old.difference_range, old.regions.len()), (None, None, 0));
    let plain = old.to_text();
    assert!(!plain.contains("map_range") && !plain.contains("difference_range") && !plain.contains("regions"), "{plain}");
    // The refusals, each in plain words.
    let wrong = |t: &str| SetFile::read(t).unwrap_err().iter().map(|x| x.message.clone()).collect::<Vec<_>>().join("; ");
    let e = wrong(&text.replace("map_range = { from = 1200.0, to = 20000.0 }", "map_range = { from = 20000.0, to = 1200.0 }"));
    println!("{e}");
    assert!(e.contains("map_range: from = 20000 is not below to = 1200"), "{e}");
    let e = wrong(&text.replace("difference_range = { from = -5000.0, to = 5000.0 }", "difference_range = { from = 5000.0, to = 5000.0 }"));
    assert!(e.contains("difference_range: from = 5000 is not below to = 5000"), "{e}");
    let e = wrong(&text.replace("x = [18.95, 19.25]", "x = [19.25, 18.95]"));
    println!("{e}");
    assert!(e.contains("region \"m1 transition\": x = [19.25, 18.95] runs backwards"), "{e}");
    let e = wrong(&text.replace("variable = [20.0, 26.0]", "variable = [26.0, 20.0]"));
    println!("{e}");
    assert!(e.contains("region \"m1 transition\": variable = [26, 20] runs backwards"), "{e}");
    let e = wrong(&text.replace("map_range = { from = 2000.0, to = 70000.0 }", "map_range = { from = 70000.0, to = 2000.0 }"));
    assert!(e.contains("region \"m1 transition\": map_range: from = 70000 is not below to = 2000"), "{e}");
    let e = wrong(&text.replace("name = \"m3\"", "name = \"m1 transition\""));
    println!("{e}");
    assert!(e.contains("two regions are named \"m1 transition\""), "{e}");
    assert!(wrong(&text.replace("x_quantity = \"2theta\"", "x_quantity = \"tth\"")).contains("x_quantity is"));
    assert!(wrong(&text.replace("window = \"m1\"", "colour = \"m1\"")).contains("\"colour\" is not a key of [[regions]]"));
    // The map of the set as drs plots map draws it: rows 4–30 K, columns 16–46° by 0.1°; the region around the peak at
    // 19.1° from 20 to 26 K is rows 8–14 (20, 21, … 26 K) and columns 30–32 (19.0, 19.1, 19.2°).
    let set = f.open(&d, None, &|_| Vec::new()).unwrap();
    let all = plots::open_buf(&tb().join("Tb_vs_T.buf"), None).unwrap();
    let grid = plots::common_grid(&all, plots::GridRange::Overlap, plots::GridStep::Finest).unwrap();
    let map20 = plots::map(&all, &grid, YScale::Linear, None).unwrap();
    let x = MapX::of(&all);
    assert_eq!(x.quantity, Some(XQuantity::TwoTheta));
    let b = plots::map_region(&map20, &f.regions[0], &x).unwrap();
    assert_eq!(b.rows, (8..=14).collect::<Vec<_>>());
    assert_eq!(b.rows.iter().map(|&i| map20.rows[i].value).collect::<Vec<_>>(), [20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0]);
    assert_eq!(b.columns, [30, 31, 32]);
    assert_eq!(b.x, (18.95, 19.25));
    // Its range is the block's cells as the matrix holds them: nothing recomputed.
    let cells: Vec<f64> = b.rows.iter().flat_map(|&i| b.columns.iter().map(move |&j| (i, j))).filter_map(|(i, j)| map20.z[i][j]).collect();
    assert_eq!(cells.len(), 21);
    let (lo, hi) = cells.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(a, z), &v| (a.min(v), z.max(v)));
    assert_eq!(b.range, (lo, hi));
    println!("block rows {:?} columns {:?} x {:?} range {:?}; the map's {:?}", b.rows, b.columns, b.x, b.range, map20.range);
    assert!(b.range.1 <= map20.range.1 && b.range.0 >= map20.range.0 && b.range.1 > 10.0 * b.range.0, "the peak inside");
    // The same region in d: refused while the files state no λ, the members named; with G4.1's λ, the same block.
    let lambda = 2.4266_f64;
    let in_d = |t: f64| lambda / (2.0 * (t / 2.0).to_radians().sin());
    let r_d = Region { x: [in_d(19.25), in_d(18.95)], x_quantity: Some(XQuantity::D), ..f.regions[0].clone() };
    let e = plots::map_region(&map20, &r_d, &x).unwrap_err();
    println!("{e}");
    assert!(e.starts_with("region \"m1 transition\" is in d and the map in 2theta: converting it needs each member's wavelength, and these state none: tb_t10.dat, tb_t11.dat"), "{e}");
    assert!(e.contains("tb_t16.dat") && !e.contains("tb_t17.dat") && !e.contains("tb_t09.dat"), "only the members of the region's rows: {e}");
    let with_l = MapX { lambda: vec![Some(lambda); all.members.len()], ..x.clone() };
    let bd = plots::map_region(&map20, &r_d, &with_l).unwrap();
    assert_eq!((&bd.rows, &bd.columns, bd.range), (&b.rows, &b.columns, b.range));
    assert!((bd.x.0 - 18.95).abs() < 1e-9 && (bd.x.1 - 19.25).abs() < 1e-9, "{:?}", bd.x);
    // In Q too, through the same λ.
    let in_q = |t: f64| 4.0 * std::f64::consts::PI * (t / 2.0).to_radians().sin() / lambda;
    let r_q = Region { x: [in_q(18.95), in_q(19.25)], x_quantity: Some(XQuantity::Q), ..f.regions[0].clone() };
    assert_eq!(plots::map_region(&map20, &r_q, &with_l).unwrap().columns, b.columns);
    // The set-level check, from the files: a region in d on a set in 2θ whose files state no λ is refused.
    let fd = SetFile { regions: vec![r_d.clone()], ..f.clone() };
    let e = fd.check_regions(&set).unwrap_err();
    println!("{e}");
    assert!(e.contains("region \"m1 transition\" is in d and the set in 2theta") && e.contains("tb_t20.dat"), "{e}");
    assert!(f.check_regions(&set).is_ok(), "regions in the set's own x need no λ");
    // A region outside the map: refused, said.
    let away = Region { variable: [40.0, 50.0], ..f.regions[0].clone() };
    assert!(plots::map_region(&map20, &away, &x).unwrap_err().contains("no row of the map has its value of the variable from 40 to 50"));
    let away = Region { x: [50.0, 60.0], ..f.regions[0].clone() };
    assert!(plots::map_region(&map20, &away, &x).unwrap_err().contains("no column of the map lies from 50 to 60"));
    // A warning, not a refusal, for a range wholly outside the matrix's values.
    let w = plots::range_warning("the map's colour range", &Range { from: 1.0e6, to: 2.0e6 }, &map20).unwrap();
    println!("{w}");
    assert!(w.contains("lies wholly outside the map's values"), "{w}");
    assert!(plots::range_warning("the map's colour range", &Range { from: 1200.0, to: 1.0e6 }, &map20).is_none(), "partly inside is no warning");
    assert!(SetFile::read(&text.replace("from = 1200.0, to = 20000.0", "from = 1.0e6, to = 2.0e6")).is_ok(), "read, not refused");
    // The command line: describe says the clip and the regions; region prints the block.
    let drv = d.join("Plots/pl37.drv");
    std::fs::write(&drv, &text).unwrap();
    let o = drs(&["plots", "describe", drv.to_str().unwrap()]);
    let out = String::from_utf8_lossy(&o.stdout).to_string();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    println!("{out}");
    assert!(out.contains("the map's colour range clipped to 1200 to 20000 on the sqrt scale") && out.contains("difference map's colour range clipped to -5000 to 5000"), "{out}");
    assert!(out.contains("region \"m1 transition\": 2theta 18.95 to 19.25, T 20 to 26 K, its colour range on the map 2000 to 70000 (sqrt scale), drawn from the window m1"), "{out}");
    assert!(out.contains("region \"m3\": 2theta (the set's x) 31.9 to 33.3"), "{out}");
    let back: SetFile = serde_json::from_slice(&drs(&["plots", "describe", drv.to_str().unwrap(), "--json"]).stdout).unwrap();
    assert_eq!(back, f);
    let o = drs(&["plots", "region", drv.to_str().unwrap(), "--name", "m3", "--json"]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let m3: plots::MapBlock = serde_json::from_slice(&o.stdout).unwrap();
    assert_eq!(m3.rows.len(), 4, "the set file's four members, 1.5 K typed: every row");
    let o = drs(&["plots", "region", drv.to_str().unwrap(), "--name", "none"]);
    assert!(!o.status.success() && String::from_utf8_lossy(&o.stderr).contains("no region \"none\" in the set file (\"m1 transition\", \"m3\")"));
}

/// P65 (L416, L417, L419): a region keeps a colour range per view, `map_range` on the map and `difference_range` on the
/// difference map, each on its own map's scale; a file with one or none reads as before.
#[test]
fn p65_a_regions_range_per_view() {
    use plots::{Range, Region, XQuantity};
    let regions = "map_scale = \"log\"\n\n[[regions]]\nname = \"m1 transition\"\nx = [18.95, 19.25]\nx_quantity = \"2theta\"\nvariable = [20.0, 26.0]\nmap_range = { from = 3.434, to = 4.661 }\ndifference_range = { from = -40.0, to = 60.0 }\nwindow = \"m1\"\n\n[[regions]]\nname = \"only on the difference\"\nx = [31.9, 33.3]\nvariable = [1.5, 30.0]\ndifference_range = { from = -10.0, to = 10.0 }\n";
    let text = EXAMPLE.replace("\n[[member]]\nfile = \"Data/tb_t20.dat\"", &format!("{regions}\n[[member]]\nfile = \"Data/tb_t20.dat\""));
    let f = SetFile::read(&text).unwrap();
    assert_eq!(
        f.regions[0],
        Region { name: "m1 transition".into(), x: [18.95, 19.25], x_quantity: Some(XQuantity::TwoTheta), variable: [20.0, 26.0], map_range: Some(Range { from: 3.434, to: 4.661 }), difference_range: Some(Range { from: -40.0, to: 60.0 }), window: Some("m1".into()) }
    );
    assert_eq!((f.regions[1].map_range, f.regions[1].difference_range), (None, Some(Range { from: -10.0, to: 10.0 })));
    // Written and read back equal, both ranges in the text, each once.
    let back = f.to_text();
    assert_eq!(SetFile::read(&back).unwrap(), f);
    assert!(back.contains("map_range = { from = 3.434, to = 4.661 }\ndifference_range = { from = -40.0, to = 60.0 }\nwindow = \"m1\"\n"), "{back}");
    // Through serde, as the app and drs plots describe --json carry it; a region without the key takes none.
    let v = serde_json::to_value(&f.regions[0]).unwrap();
    assert_eq!(v["difference_range"], serde_json::json!({ "from": -40.0, "to": 60.0 }));
    let mut old = v.clone();
    old.as_object_mut().unwrap().remove("difference_range");
    assert_eq!(serde_json::from_value::<Region>(old).unwrap().difference_range, None);
    // A file written before P65 reads as before: the region's one range is the map's.
    let before = SetFile::read(&EXAMPLE.replace("\n[[member]]\nfile = \"Data/tb_t20.dat\"", "\n[[regions]]\nname = \"r\"\nx = [18.0, 19.0]\nvariable = [2.0, 4.0]\nmap_range = { from = 1.0, to = 2.0 }\n\n[[member]]\nfile = \"Data/tb_t20.dat\"")).unwrap();
    assert_eq!((before.regions[0].map_range, before.regions[0].difference_range), (Some(Range { from: 1.0, to: 2.0 }), None));
    assert!(!before.to_text().contains("difference_range"));
    // Refused as the map's range is: from not below to.
    let wrong = SetFile::read(&text.replace("difference_range = { from = -40.0, to = 60.0 }", "difference_range = { from = 60.0, to = -40.0 }")).unwrap_err();
    assert!(wrong[0].message.starts_with("region \"m1 transition\": difference_range: from = 60 is not below to = -40"), "{}", wrong[0].message);
    // Described with the map each range is on and that map's scale (the difference map on √ for a log map).
    let words = f.words().join("\n");
    assert!(words.contains("region \"m1 transition\": 2theta 18.95 to 19.25, T 20 to 26 K, its colour range on the map 3.434 to 4.661 (log scale), its colour range on the difference map -40 to 60 (sqrt scale), drawn from the window m1"), "{words}");
    assert!(words.contains("region \"only on the difference\": 2theta (the set's x) 31.9 to 33.3, T 1.5 to 30 K, its colour range on the difference map -10 to 10 (sqrt scale)"), "{words}");
    // drs plots describe says it; drs plots region gives the range of the map it is on and names the other.
    let d = project("p65");
    let path = d.join("Plots/p65.drv");
    std::fs::write(&path, &text).unwrap();
    let o = drs(&["plots", "describe", path.to_str().unwrap()]);
    assert!(o.status.success() && String::from_utf8_lossy(&o.stdout).contains("its colour range on the difference map -40 to 60 (sqrt scale)"), "{}", String::from_utf8_lossy(&o.stdout));
    let on_map = String::from_utf8_lossy(&drs(&["plots", "region", path.to_str().unwrap(), "--name", "m1 transition"]).stdout).to_string();
    assert!(on_map.contains("its colour range on the map 3.434 to 4.661 (log scale)\nits colour range on the difference map -40 to 60 (not this map's: shown with --reference)\n"), "{on_map}");
    let on_diff = drs(&["plots", "region", path.to_str().unwrap(), "--name", "m1 transition", "--reference", "4"]);
    let out = String::from_utf8_lossy(&on_diff.stdout);
    assert!(on_diff.status.success() && out.contains("its colour range on the difference map -40 to 60 (sqrt scale)\nits colour range on the map 3.434 to 4.661 (not this map's: shown with no --reference)\n"), "{out}{}", String::from_utf8_lossy(&on_diff.stderr));
}

#[test]
fn pl41_the_y_interval_and_the_sets_extent() {
    use plots::YInterval;
    // Kept by the set file: "own" written and read back equal, described; a .drv without it reads as shared and writes
    // none.
    let own = SetFile::read(&EXAMPLE.replace("\n[[member]]\nfile = \"Data/tb_t20.dat\"", "y_interval = \"own\"\n\n[[member]]\nfile = \"Data/tb_t20.dat\"")).unwrap();
    assert_eq!(own.y_interval, YInterval::Own);
    let text = own.to_text();
    assert!(text.contains("y_interval = \"own\"\n"), "{text}");
    assert_eq!(SetFile::read(&text).unwrap(), own);
    assert!(own.words().iter().any(|w| w == "the curves' y interval: each member's own"));
    let old = SetFile::read(EXAMPLE).unwrap();
    assert_eq!(old.y_interval, YInterval::Shared);
    assert!(!old.to_text().contains("y_interval"));
    assert!(old.words().iter().any(|w| w == "the curves' y interval: the whole set's"));
    let e = SetFile::read(&text.replace("y_interval = \"own\"", "y_interval = \"mine\"")).unwrap_err();
    assert!(e[0].message.contains("y_interval is \"shared\" (the whole set's) or \"own\" (each member's)"), "{}", e[0].message);
    assert_eq!(serde_json::to_value(own.y_interval).unwrap(), serde_json::json!("own"));
    // The set's extent on the log scale is the colour range of its map there, as drs plots map gives it.
    let files: Vec<(String, Option<String>)> = (1..=20).map(|k| (format!("tb_t{k:02}.dat"), None)).collect();
    let set = plots::open(&tb(), &files, None, None).unwrap();
    let grid = plots::common_grid(&set, plots::GridRange::Overlap, plots::GridStep::Finest).unwrap();
    let map = plots::map(&set, &grid, YScale::Log, None).unwrap();
    let ext = plots::extent(&set, YScale::Log, None).unwrap();
    assert_eq!(ext, map.range);
    let o = drs(&["plots", "map", tb().join("Tb_vs_T.buf").to_str().unwrap(), "--scale", "log"]);
    let out = String::from_utf8_lossy(&o.stdout);
    assert!(o.status.success() && out.contains(&format!("colour range {} to {}\n", ext.0, ext.1)), "{out}{}", String::from_utf8_lossy(&o.stderr));
    // Of the members named: one member's, within the set's; on the linear scale its counts.
    let one = plots::extent(&set, YScale::Linear, Some(&[3])).unwrap();
    assert_eq!(one, drs::data::extent(&set.members[3].pattern.y).unwrap());
    let lin = plots::extent(&set, YScale::Linear, None).unwrap();
    assert!(lin.0 <= one.0 && one.1 <= lin.1);
    assert_eq!(plots::extent(&set, YScale::Linear, Some(&[])), None, "no member, no extent");
}

/// The ESRF run's names (`plots_esrf`), each file written here from 0 to 27° in 0.01° steps.
fn esrf_like(name: &str, from: f64) -> plots::OpenedFile {
    let text: String = (0..=((27.0 - from) / 0.01).round() as usize).map(|i| format!("{:.2} {} 10.0\n", from + 0.01 * i as f64, 100 + i % 50)).collect();
    plots::OpenedFile { name: name.into(), text, ..Default::default() }
}

const ESRF_CUT: &str = r#"format = "drs plot set 1"
title = "three ESRF names, cut"
variable = { name = "T", unit = "°C" }
x_interval = [2.0, 20.0]

[[member]]
file = "scan_001_T0900C_t00000s_long.xye"

[[member]]
file = "scan_002_T0925C_t00810s_long.xye"

[[member]]
file = "scan_003_T0950C_t01622s_long.xye"
"#;

#[test]
fn pl42_the_x_interval_of_the_analysis() {
    use plots::{MapX, Region, Window, WindowBackground};
    let files: Vec<plots::OpenedFile> = ["scan_001_T0900C_t00000s_long.xye", "scan_002_T0925C_t00810s_long.xye", "scan_003_T0950C_t01622s_long.xye"].iter().map(|n| esrf_like(n, 0.0)).collect();
    // Read, written, described; a .drv without the field reads as before and writes none.
    let f = SetFile::read(ESRF_CUT).unwrap();
    assert_eq!(f.x_interval, Some([2.0, 20.0]));
    assert!(f.to_text().contains("x_interval = [2.0, 20.0]\n"), "{}", f.to_text());
    assert_eq!(SetFile::read(&f.to_text()).unwrap(), f);
    assert!(f.words().iter().any(|w| w == "the analysis cut to 2theta 2 to 20 on opening: every member's points outside it left out, the files untouched"), "{:?}", f.words());
    let old = SetFile::read(EXAMPLE).unwrap();
    assert_eq!(old.x_interval, None);
    assert!(!old.to_text().contains("x_interval"));
    assert!(SetFile::read(&ESRF_CUT.replace("[2.0, 20.0]", "[20.0, 2.0]")).unwrap_err()[0].message.contains("x_interval = [20, 2] runs backwards"));
    // Opened: every member cut to 2–20°, the map's columns spanning exactly that, in °C from the names (W8).
    let set = f.open_from(&files, Some("xye"), &|_| Vec::new()).unwrap();
    assert_eq!(set.x_interval, Some([2.0, 20.0]));
    assert!(set.members.iter().all(|m| m.facts.range == (2.0, 20.0) && m.facts.points == 1801), "{:?}", set.members[0].facts.range);
    let grid = plots::common_grid(&set, plots::GridRange::Overlap, plots::GridStep::Finest).unwrap();
    let map = plots::map(&set, &grid, YScale::Linear, None).unwrap();
    assert_eq!((map.x[0], *map.x.last().unwrap(), map.x.len()), (2.0, 20.0, 1801));
    assert_eq!(map.rows.iter().map(|r| r.value).collect::<Vec<_>>(), [900.0, 925.0, 950.0]);
    // A window outside the interval refused as such; one inside read.
    let beamstop = Window { name: "beamstop".into(), from: 0.5, to: 1.5 };
    let e = plots::windows(&set, &[beamstop], &WindowBackground::Line, &[]).unwrap_err();
    assert_eq!(e, "window \"beamstop\", 0.5 to 1.5, lies outside the analysis's x interval 2 to 20");
    assert!(plots::windows(&set, &[Window { name: "in".into(), from: 5.0, to: 6.0 }], &WindowBackground::Line, &[]).is_ok());
    // Kept by the set file written from the set, and by the proposals the host asks with it.
    assert_eq!(SetFile::of(&set).x_interval, Some([2.0, 20.0]));
    let p = plots::propose(&files, Some([2.0, 20.0]));
    assert_eq!(p.proposals.iter().map(|x| (x.kind.as_str(), x.set.x_interval)).collect::<Vec<_>>(), [("by_temperature", Some([2.0, 20.0])), ("against_time", Some([2.0, 20.0]))]);
    // A member whose range does not reach the interval refused by name; the files themselves untouched.
    let mut late = files.clone();
    late[1] = esrf_like("scan_002_T0925C_t00810s_long.xye", 21.0);
    assert_eq!(f.open_from(&late, Some("xye"), &|_| Vec::new()).unwrap_err(), "scan_002_T0925C_t00810s_long.xye: its range, 21 to 27, does not reach the x interval 2 to 20");
    assert!(files[0].text.starts_with("0.00 "));
    // Tb₂Pd₂Sn cut to 18–34°: the map of the cut set is the whole map's block of those columns (map_region).
    let names: Vec<(String, Option<String>)> = (1..=20).map(|k| (format!("tb_t{k:02}.dat"), None)).collect();
    let whole = plots::open(&tb(), &names, None, None).unwrap();
    let mut cut = whole.clone();
    cut.cut([18.0, 34.0], None).unwrap();
    let wmap = plots::map(&whole, &plots::common_grid(&whole, plots::GridRange::Overlap, plots::GridStep::Finest).unwrap(), YScale::Linear, None).unwrap();
    let cmap = plots::map(&cut, &plots::common_grid(&cut, plots::GridRange::Overlap, plots::GridStep::Finest).unwrap(), YScale::Linear, None).unwrap();
    let region = Region { name: "18–34°".into(), x: [18.0, 34.0], x_quantity: None, variable: [0.0, 100.0], map_range: None, difference_range: None, window: None };
    let block = plots::map_region(&wmap, &region, &MapX::of(&whole)).unwrap();
    assert_eq!(cmap.x, block.columns.iter().map(|&j| wmap.x[j]).collect::<Vec<_>>());
    assert_eq!(block.rows.len(), cmap.rows.len());
    for (r, &w) in block.rows.iter().enumerate() {
        assert_eq!(cmap.z[r], block.columns.iter().map(|&j| wmap.z[w][j]).collect::<Vec<_>>(), "row {r}");
    }
    let finite = cmap.z.iter().flatten().flatten().copied();
    assert_eq!(finite.fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), v| (a.min(v), b.max(v))), block.range);
    // drs plots map --x-interval gives the same columns; describe says the cut.
    let o = drs(&["plots", "map", tb().join("Tb_vs_T.buf").to_str().unwrap(), "--x-interval", "18:34", "--json"]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let m: plots::Map = serde_json::from_slice(&o.stdout).unwrap();
    assert_eq!(m.x, cmap.x);
    let d = project("pl42");
    let path = d.join("Plots/pl42.drv");
    std::fs::write(&path, EXAMPLE).unwrap();
    let o = drs(&["plots", "describe", path.to_str().unwrap(), "--x-interval", "18:34"]);
    assert!(String::from_utf8_lossy(&o.stdout).contains("the analysis cut to 2theta 18 to 34 on opening"), "{}{}", String::from_utf8_lossy(&o.stdout), String::from_utf8_lossy(&o.stderr));
}
