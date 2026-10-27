//! TLS displacements of rigid groups (P41, Publish's P47.2; L367) against FullProf 8.20's runs of `runs/version2_runs.py`
//! and the tutorials' (`manuals/fullprof/runs/rigid`, run on the Mac: L346, L351).
//!
//! - `betas_against_fullprof_ort`: every atom's β from T, L, S (Schomaker & Trueblood 1968, `rigid_body::tls_betas`)
//!   against the BETA lines FullProf writes to its `.ort` (FullProf's manual p. 39), as imported from each run's input.
//!   hmt's run and pdd's five (centre at the group's origin, P6 = −5.0): within the print, 10⁻⁵. Urea's six (P6 = −5.01,
//!   the centre FullProf computes): off by up to 5.45·10⁻⁵ (ureac, ureacyl) with the centre at the main group's centre
//!   of mass by CrysFML's atomic weights (z 0.32957), within the print only for z in 0.33032–0.33036. FullProf prints
//!   its centre of mass at z 0.326, which matches no weighting found (masses with or without the hydrogens, atomic
//!   numbers, neutron scattering lengths); 0.3303 is the group's origin z 0.32811 less (0.326 − 0.32811), as if the
//!   offset were added where it is taken away, but the centre of mass that would give it, 0.3259, matches no weighting
//!   either. Frozen in `UREA`, listed for Josu.
//! - `urea_test_chi2`: urea_test's χ² (13.239 against FullProf's 13.27, 0.23 %; its profile within rms 3.7·10⁻⁴) comes
//!   mostly from that centre: with z = 0.33034 the core gives 13.278, 0.06 % from FullProf's.
//! - `derivatives`: every TLS value's, the group values' and the cell's Jacobian columns against differences of the
//!   pattern (urea, hmt, pdd), within 10⁻⁴ of the column; urea's satellites follow their parent
//!   (`rigid_body::placement_derivative`).
//! - `control_file_pcr_and_results`: the TLS values through the control file (each where `locate` finds it), the
//!   `.pcr` export (Furth = 21, its codewords) and the results file (P47.2: T, L, S in FullProf's printed units with su's
//!   and the parameter each moves, s22 tied to s11 as urea's codewords tie them; the atoms' U with su's).

use drs::control::Control;
use drs::jacobian::Field;
use drs::lsq::{refine, Options};
use drs::powder::Model;
use drs::refinable::Refinement;
use std::path::{Path, PathBuf};

fn dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs/rigid")
}

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display())).iter().map(|&b| b as char).collect()
}

/// A run's input: VARIANT.pcr.in, or the base run's .pcr.
fn load(run: &str, file: &str) -> drs::powder::Loaded {
    let p = drs::pcr::read(&latin1(&dir().join(file))).unwrap_or_else(|e| panic!("{file}: {e}"));
    Model::from_pcr_fullprof(&p, &dir(), run).unwrap_or_else(|e| panic!("{run}: {e}"))
}

/// The BETA lines of a `.ort`, in its order (the phase's atoms').
fn ort_betas(file: &str) -> Vec<[f64; 6]> {
    let text = latin1(&dir().join(file));
    text.lines()
        .filter(|l| l.starts_with("BETA"))
        .map(|l| {
            // "BETA      C  1  0.00932 0.00932 0.00516 0.00033-0.00000 0.00000": the six numbers, glued signs apart.
            let rest = &l[l.find(|c: char| c.is_ascii_digit()).unwrap()..];
            let nums: Vec<f64> = rest.replace('-', " -").split_whitespace().skip(1).map(|t| t.parse().unwrap()).collect();
            std::array::from_fn(|i| nums[nums.len() - 6 + i])
        })
        .collect()
}

/// Urea's runs against their `.ort` (L367): the largest |Δβ| frozen.
const UREA: f64 = 5.46e-5;

fn runs() -> Vec<(&'static str, &'static str, &'static str)> {
    vec![
        ("hmt", "hmt.pcr.in", "hmt1.ort"),
        ("pdd", "pdd.pcr.in", "pdd1.ort"),
        ("pdd_l", "pdd_l.pcr.in", "pdd_l1.ort"),
        ("pdd_t", "pdd_t.pcr.in", "pdd_t1.ort"),
        ("pdd_t_rot", "pdd_t_rot.pcr.in", "pdd_t_rot1.ort"),
        ("pdd_calc", "pdd_calc.pcr", "pdd_calc1.ort"),
        ("urea", "urea.pcr.in", "urea1.ort"),
        ("urea_conv", "urea_conv.pcr.in", "urea_conv1.ort"),
        ("ureac", "ureac.pcr.in", "ureac1.ort"),
        ("ureacyl", "ureacyl.pcr.in", "ureacyl1.ort"),
        ("urea_notls", "urea_notls.pcr.in", "urea_notls1.ort"),
        ("urea_calc", "urea_calc.pcr", "urea_calc1.ort"),
    ]
}

#[test]
fn betas_against_fullprof_ort() {
    for (run, file, ort) in runs() {
        let m = load(run, file).model;
        let ph = &m.phases[0];
        let drs::powder::Intensities::Structure(st, _) = &ph.intensities else { panic!() };
        let want = ort_betas(ort);
        assert_eq!(want.len(), st.atoms.len(), "{run}");
        let worst = st.atoms.iter().zip(&want).map(|(a, w)| a.beta.unwrap().iter().zip(w).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max)).fold(0.0, f64::max);
        let tls = ph.rigid.iter().find_map(|g| g.displacements.as_ref()).unwrap();
        println!("{run:12}: {} atoms, centre less origin {:?}, worst |Δβ| {worst:.3e}", st.atoms.len(), tls.centre);
        let allowed = if run.starts_with("urea") { UREA } else { 5.5e-6 };
        assert!(worst <= allowed, "{run}: |Δβ| {worst:.2e} beyond {allowed:.1e}");
    }
}

#[test]
fn urea_test_chi2() {
    let l = load("urea_test", "urea_test.pcr");
    let (mut m, data) = (l.model, l.data.unwrap());
    let chi2 = |m: &Model| {
        let y = m.calculate(&data.x).ycalc;
        drs::powder::agreement(&data, &y, &m.excluded, 0).3
    };
    let ours = chi2(&m);
    // The centre FullProf's β imply (`betas_against_fullprof_ort`).
    let g = m.phases[0].rigid.iter().position(|g| g.displacements.is_some()).unwrap();
    let origin = m.phases[0].rigid[g].origin;
    m.phases[0].rigid[g].displacements.as_mut().unwrap().centre = [0.0, 0.0, 0.33034 - origin[2]];
    let k = drs::rigid_body::values_of(&m.phases[0].rigid[g]) - 1;
    let v = m.get(Field::Rigid(0, g, k));
    m.set(Field::Rigid(0, g, k), v);
    let theirs = chi2(&m);
    println!("urea_test χ²: {ours:.4} with the core's centre, {theirs:.4} with z = 0.33034; FullProf 13.27");
    assert!((ours - 13.27).abs() > 0.03 && (theirs - 13.27).abs() < 0.008, "{ours} {theirs}");
}

#[test]
fn derivatives() {
    for (run, file) in [("urea", "urea.pcr.in"), ("hmt", "hmt.pcr.in"), ("pdd", "pdd.pcr.in")] {
        let l = load(run, file);
        let (m, data) = (l.model, l.data.unwrap());
        let g = m.phases[0].rigid.iter().position(|g| g.displacements.is_some()).unwrap();
        let n = drs::rigid_body::values_of(&m.phases[0].rigid[g]);
        let mut fields: Vec<Field> = (0..n).map(|k| Field::Rigid(0, g, k)).filter(|&f| m.all_fields().contains(&f) && m.fixed_by(f).is_none()).collect();
        fields.extend((0..3).map(|k| Field::Cell(0, k)));
        let (_, jac) = m.derivatives(&data.x, &fields);
        let mut checked = 0;
        for (j, &f) in fields.iter().enumerate() {
            let h = 1e-6 * m.get(f).abs().max(1.0);
            let (mut a, mut b) = (m.clone(), m.clone());
            a.set(f, m.get(f) + h);
            b.set(f, m.get(f) - h);
            let (ya, yb) = (a.calculate(&data.x).ycalc, b.calculate(&data.x).ycalc);
            let num: Vec<f64> = (0..data.x.len()).map(|i| (ya[i] - yb[i]) / (2.0 * h)).collect();
            let top = num.iter().fold(0.0f64, |s, v| s.max(v.abs()));
            // A column whose differences are rounding (10⁴ ε of the pattern over the step) is zero by symmetry.
            let floor = 1e-12 * ya.iter().fold(0.0f64, |s, v| s.max(v.abs())) / h;
            if top < floor {
                assert!(jac[j].iter().all(|v| v.abs() < floor), "{run} {f:?}: rounding {floor:.1e}, but the Jacobian has {:.1e}", jac[j].iter().fold(0.0f64, |s, v| s.max(v.abs())));
                continue;
            }
            let worst = (0..data.x.len()).map(|i| (num[i] - jac[j][i]).abs()).fold(0.0, f64::max);
            assert!(worst < 1e-4 * top, "{run} {f:?}: off by {worst:.3e} of {top:.3e}");
            checked += 1;
        }
        println!("{run}: {checked} columns against differences");
        assert!(checked >= 10, "{run}: {checked}");
    }
}

#[test]
fn control_file_pcr_and_results() {
    for (run, file) in [("urea", "urea.pcr.in"), ("hmt", "hmt.pcr.in"), ("pdd", "pdd.pcr.in")] {
        let p = drs::pcr::read(&latin1(&dir().join(file))).unwrap();
        let c = Control::from_pcr(&p, &dir(), run, &dir()).unwrap_or_else(|e| panic!("{run}: {e}"));
        let data = load(run, file).data.unwrap();
        let g = c.model.phases[0].rigid.iter().position(|g| g.displacements.is_some()).unwrap();
        let n = 6 + 3 * c.model.phases[0].rigid[g].members.len();
        let values: Vec<f64> = (0..21).map(|i| c.model.get(Field::Rigid(0, g, n + i))).collect();
        // The control file, and each value where the file writes it.
        let text = c.to_text();
        let back = Control::from_text_in(&text, &dir()).unwrap_or_else(|e| panic!("{run}: {e:?}"));
        assert_eq!(back.model.to_text(), c.model.to_text(), "{run}: the model read back");
        assert_eq!(back.model.calculate(&data.x).ycalc, c.model.calculate(&data.x).ycalc, "{run}: the pattern read back");
        assert_eq!(back.refinement, c.refinement, "{run}: the free parameters read back");
        let located: Vec<_> = drs::control::locate::all(&text, Some(&dir())).into_iter().filter(|(a, ..)| a.contains("rigid_tls")).collect();
        assert_eq!(located.len(), 21, "{run}: 21 TLS values located");
        for (a, _, r) in &located {
            let f = drs::refinable::field_of(a, &c.model).unwrap_or_else(|e| panic!("{a}: {e}"));
            assert_eq!(text[r.clone()].parse::<f64>().unwrap(), c.model.get(f), "{a}");
        }
        // The .pcr export: Furth = 21, the same values, β and free parameters.
        let e = drs::control::pcr_export::export(&c, run, Some(&data)).unwrap_or_else(|er| panic!("{run}: {er}"));
        let p2 = drs::pcr::read(&e.pcr).unwrap_or_else(|er| panic!("{run}: {er}\n{}", e.pcr));
        assert_eq!(p2.phases[0].flag("Furth"), 21, "{run}");
        let m2 = Model::from_pcr(&p2, &dir(), run).unwrap_or_else(|er| panic!("{run}: {er}")).model;
        assert_eq!((0..21).map(|i| m2.get(Field::Rigid(0, g, n + i))).collect::<Vec<_>>(), values, "{run}: through the .pcr");
        let betas = |m: &Model| match &m.phases[0].intensities {
            drs::powder::Intensities::Structure(st, _) => st.atoms.iter().map(|a| a.beta).collect::<Vec<_>>(),
            _ => panic!(),
        };
        let worst = betas(&m2).iter().zip(betas(&c.model)).map(|(a, b)| a.unwrap().iter().zip(b.unwrap()).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max)).fold(0.0, f64::max);
        assert!(worst < 1e-9, "{run}: β through the .pcr, {worst:.1e}");
        let free = |rf: &Refinement, m: &Model| rf.build(m).unwrap().free.iter().map(|f| f.fields.clone()).collect::<Vec<_>>();
        assert_eq!(free(&Refinement::from_pcr(&p2, &m2).unwrap(), &m2), free(&c.refinement, &c.model), "{run}: the free parameters through the .pcr");
        // The results file after a cycle.
        let mut cr = c.clone();
        let set = c.refinement.build(&cr.model).unwrap();
        let r = refine(&mut cr.model, &data, &set, Options { max_cycles: 1, eps: 1e-3 }, |_| {});
        let (cif, _) = drs::results::write(&cr, &data, Some(&set), Some(&r), run).unwrap();
        let block = |item: &str| -> Vec<Vec<String>> {
            let i = cif.find(&format!("\n{item}\n")).unwrap_or_else(|| panic!("{run}: {item}"));
            let rest = &cif[i + 1..];
            let body = rest.lines().skip_while(|l| l.starts_with("_drs_")).take_while(|l| !l.starts_with("loop_") && !l.starts_with('_') && !l.is_empty());
            body.map(|l| l.split_whitespace().map(String::from).collect()).collect()
        };
        let rows = block("_drs_rigid_tls_value_group_id");
        assert_eq!(rows.len(), 21, "{run}");
        let row = |t: &str, c: &str| rows.iter().find(|r| r[1] == t && r[2] == c).unwrap().clone();
        // L in deg², as FullProf's .out prints it.
        let printed = row("L", "11")[3].split('(').next().unwrap().to_string();
        let l11: f64 = printed.parse().unwrap();
        let last = 10f64.powi(-(printed.split('.').nth(1).map_or(0, |d| d.len()) as i32));
        let want = cr.model.get(Field::Rigid(0, g, n + 6)) * (180.0 / std::f64::consts::PI).powi(2);
        assert!((l11 - want).abs() <= 0.5 * last + 1e-12, "{run}: L11 {printed} against {want}");
        if run == "urea" {
            let (s11, s22) = (row("S", "11"), row("S", "22"));
            assert!(s11[4] != "." && s11[4] == s22[4], "{run}: s22 moved by s11's parameter: {s11:?} {s22:?}");
            assert_eq!(s22[5].parse::<f64>().unwrap() / s11[5].parse::<f64>().unwrap(), -1.0, "{run}: s22 = −s11");
        }
        let atoms = block("_drs_rigid_tls_atom_label");
        let members = cr.model.phases[0].rigid.iter().enumerate().map(|(h, _)| drs::rigid_body::tls_members(&cr.model.phases[0].rigid, h).len()).sum::<usize>();
        assert_eq!(atoms.len(), members, "{run}: every TLS atom's U");
        assert!(atoms.iter().all(|r| r[1..].iter().any(|u| u.contains('('))), "{run}: U's with su's: {atoms:?}");
        println!("{run}: {} TLS values, {} atoms' U, worst |Δβ| through the .pcr {worst:.1e}", rows.len(), atoms.len());
    }
}
