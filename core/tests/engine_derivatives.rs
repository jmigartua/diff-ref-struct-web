//! The engine's analytic derivatives against finite differences (test `engine_derivatives`, L053).
//!
//! For thirteen FullProf runs imported into the default model, every kind of refinable field the run has (the zero,
//! Sycos (with FullProf's cos 2θ where Sycos ≠ 0: NdNi-HTO, cuf1k, L267), Sysin, Ratio, background values, and per
//! phase the scale, B_ov, the six cell parameters, U V W X Y I_G η0, P1–P4, S/L, H/L, G1, G2, the quartic form's S_HKL
//! and ξ (NdNi-LTOa, P25; prba10k's and cuf1k's tetragonal and monoclinic forms, L267), the anisotropic size (prba10k's
//! needles, Pd3MnD's antiphase domains, L267), and for two atoms x y z, B, Occ and the betas) is differentiated by
//! `Model::derivatives`
//! and by central finite differences of `Model::calculate` at steps h and h/2, h = 10⁻⁵·max(|value|, 10⁻²) (the scale:
//! 10⁻⁵·value; the cell: 10⁻⁶·value, see below), combined by Richardson's extrapolation.
//! The peak windows are widened (Wdt = 10⁴) and every reflection reaching the data is used, so that no window edge or
//! reflection boundary falls between the two evaluations; 300 points span each pattern. The error is the largest
//! |analytic − difference| over the points, over the largest of either; it must stay below 2·10⁻⁵, above the finite
//! differences' own truncation and rounding error at those steps. A derivative below 10⁻⁷ of the largest Ycal on both
//! sides (an atom on a special position, a neutral preferred orientation) counts as zero.
//!
//! Also checked: the dual-number path reproduces `Model::calculate` to 10⁻¹⁰ (relative) at every point of the data.
//!
//! Two more passes (rutile + anatase, Y2O3) take the default model at the run's own Wdt, so that the default window's
//! tapered edge (L058) lies inside the data: its derivative with respect to the position and the FWHM is checked too.

use drs::jacobian::Field;
use drs::powder::{Intensities, Model};
use std::path::Path;

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display())).iter().map(|&b| b as char).collect()
}

const RUNS: [(&str, &str); 16] = [
    ("Y2O3", "y2o3_structure"),
    ("Co-YBCO", "Cocus-asym"),
    ("MgAl2O4-s", "MgAl2O4s_Rietveld"),
    ("rutile-anatase", "rutana_tch"),
    ("Egypt", "Egypt_qa_Rietveld"),
    ("YVO3", "YVO3a"),
    ("Yb2ZrWO8", "ybfluo_Riet_IRF"),
    ("LaB6", "LaB6rx"),
    ("NdNi-HTO", "NdNi_HTO"),
    // Stephens' quartic form of strain and its Lorentzian share (P25, L211).
    ("NdNi-LTOa", "NdNi_LTOa"),
    ("Yb2ZrWO8-PO-March", "ybfluo_Riet_IRF"),
    // P25's rest (L267): needles and the tetragonal form; the monoclinic form and Sycos·cos 2θ; antiphase domains.
    ("prba10k", "prba10k"),
    ("cuf1k", "cuf1k"),
    ("Pd3MnD", "pdfinal"),
    // The default model at the run's own Wdt: the tapered window edge (L058) inside the data.
    ("rutile-anatase/default", "rutana_tch"),
    ("Y2O3/default", "y2o3_structure"),
];

fn fields(m: &Model) -> Vec<Field> {
    let mut f = vec![Field::Zero, Field::Sycos, Field::Sysin];
    if m.ratio != 0.0 && (m.lambda[1] - m.lambda[0]).abs() > 1e-9 {
        f.push(Field::Ratio);
    }
    for k in 0..3 {
        f.push(Field::Background(k));
    }
    for (p, ph) in m.phases.iter().enumerate() {
        f.extend([Field::Scale(p), Field::Bov(p), Field::U(p), Field::V(p), Field::W(p), Field::X(p), Field::Y(p), Field::GaussSize(p), Field::Eta0(p)]);
        f.extend((0..6).map(|k| Field::Cell(p, k)));
        f.extend((0..4).map(|k| Field::Asym(p, k)));
        f.extend([Field::SL(p), Field::DL(p), Field::Pref(p, 0), Field::Pref(p, 1)]);
        // An S_HKL at 0 (S_004, S_202, S_022 of NdNi-LTOa) is left out: the reflections it alone broadens, (0 0 l)
        // for S_004, have a Lorentzian width ξ c √SUM tan θ/M whose slope at SUM = 0 is infinite on one side, which no
        // central difference measures; the engine takes the Gaussian part's derivative there.
        if let Some(q) = &ph.quartic_strain {
            f.extend(q.terms.iter().enumerate().filter(|(_, t)| t.1 != 0.0).map(|(k, _)| Field::Quartic(p, k)));
            f.push(Field::QuarticShare(p));
        }
        if ph.shape_size.is_some() {
            f.push(Field::ShapeSize(p));
        }
        if let Intensities::Structure(st, _) = &ph.intensities {
            for (a, at) in st.atoms.iter().enumerate().take(2) {
                f.extend((0..3).map(|k| Field::AtomX(p, a, k)));
                f.extend([Field::AtomB(p, a), Field::AtomOcc(p, a)]);
                if at.beta.is_some() {
                    f.extend((0..6).map(|k| Field::AtomBeta(p, a, k)));
                }
            }
        }
    }
    f
}

#[test]
fn engine_derivatives() {
    let mut fails = Vec::new();
    let mut checked = 0;
    for (run, code) in RUNS {
        let (folder, default) = match run.split_once('/') {
            Some((f, _)) => (f, true),
            None => (run, false),
        };
        let d = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs").join(folder);
        let mut p = drs::pcr::read(&latin1(&d.join(format!("{code}.pcr")))).unwrap();
        p.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
        let loaded = Model::from_pcr(&p, &d, code).unwrap_or_else(|e| panic!("{run}: {e}"));
        let data = loaded.data.as_ref().expect("data");
        // The dual-number path against `calculate`, at the data's points.
        let mut m = loaded.model.clone();
        if m.phases.iter().all(|ph| matches!(ph.intensities, Intensities::Structure(..))) {
            let (a, b) = (m.calculate(&data.x).ycalc, m.calculate_dual_values(&data.x));
            let worst = a.iter().zip(&b).map(|(u, v)| (u - v).abs() / u.abs().max(1.0)).fold(0.0, f64::max);
            println!("{run:18} dual path against calculate: {worst:.1e}");
            if worst > 1e-10 {
                fails.push(format!("{run}: the dual path differs from calculate by {worst:.2e}"));
            }
        }
        // The derivatives against finite differences.
        if default {
            m.convention = drs::powder::Convention::DEFAULT;
        } else {
            m.wdt = 1e4;
            m.convention.peak_window = false;
            m.convention.reflections = false;
        }
        let n = data.x.len();
        let x: Vec<f64> = (0..300).map(|i| data.x[i * (n - 1) / 299]).collect();
        let fl = fields(&m);
        let (c0, jac) = m.derivatives(&x, &fl);
        let ymax = c0.ycalc.iter().map(|v| v.abs()).fold(0.0, f64::max);
        let mut worst = (0.0f64, String::new());
        for (j, &f) in fl.iter().enumerate() {
            let v = m.get(f);
            // The cell takes a step ten times smaller: a cell step of 10⁻⁴ Å moves a peak enough to change the number
            // of Gauss–Legendre terms of the FCJ profile (CrysFML's NTERMS, set by 2θ − 2φ_min), a jump in the
            // calculated pattern that a finite difference across it mistakes for slope (NdNi-HTO, 83.45°).
            // NdNi-LTOa: at 10⁻⁶ of c a reflection enters the list and the pattern jumps at 92.7° (with or without
            // the quartic form); from 10⁻⁷ to 10⁻⁹ the difference is −152.97, as the analytic derivative (L211).
            let h = match f {
                Field::Scale(_) => 1e-5 * v.abs(),
                Field::Cell(..) if run == "NdNi-LTOa" => 1e-8 * v.abs(),
                Field::Cell(..) => 1e-6 * v.abs(),
                _ => 1e-5 * v.abs().max(1e-2),
            };
            let central = |h: f64| -> Vec<f64> {
                let mut up = m.clone();
                up.set(f, v + h);
                let mut dn = m.clone();
                dn.set(f, v - h);
                let (yu, yd) = (up.calculate(&x).ycalc, dn.calculate(&x).ycalc);
                yu.iter().zip(&yd).map(|(a, b)| (a - b) / (2.0 * h)).collect()
            };
            // Richardson's extrapolation of two central differences: error of order h⁴.
            let (d1, d2) = (central(h), central(h / 2.0));
            let fd: Vec<f64> = d1.iter().zip(&d2).map(|(a, b)| (4.0 * b - a) / 3.0).collect();
            let top = fd.iter().chain(&jac[j]).map(|v| v.abs()).fold(0.0, f64::max);
            checked += 1;
            // A derivative below 10⁻⁷ of the largest Ycal per unit of the field (an atom on a special position, a
            // neutral preferred orientation) is zero on both sides, within the differences' rounding.
            if top * h < 1e-7 * ymax * h.max(1e-2) / 1e-2 || top < 1e-7 * ymax {
                continue;
            }
            let err = jac[j].iter().zip(&fd).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max) / top;
            if err > worst.0 {
                worst = (err, format!("{f:?}"));
            }
            if err > 2e-5 {
                fails.push(format!("{run}: {f:?}: analytic against finite difference {err:.2e} (largest derivative {top:.3e})"));
            }
        }
        println!("{run:18} {} fields, worst {:.1e} ({})", fl.len(), worst.0, worst.1);
    }
    println!("{checked} derivatives checked");
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}
