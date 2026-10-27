//! `Field::describe` (the app's request 39, L276, L292): every field of the models of eight FullProf runs labelled with
//! its dictionary id and name, phase, atom and component. The kinds of `Field` are matched here without a wildcard, so
//! a kind the core adds and this test does not describe stops it compiling; and every kind must occur in the runs.
//!
//! The runs: Y₂O₃ (two wavelengths, background points, two atoms labelled Y, X-ray Bragg–Brentano micro-absorption),
//! NdNi-LTOa (the quartic strain, betas, pseudo-Voigt), prba10k (needles, Chebyshev background), Pd3MnD (antiphase
//! size), the mode case A (mode amplitudes), Co-YBCO (the axial divergence), LaB6 (the Bérar–Baldinozzi P1–P4) and
//! FullProf's Le Bail fit of the spinel made a Pawley fit (the reflection intensities, P36), and arg_si (a TOF bank and
//! phase, P05), citric acid (a rigid Z-matrix group, P41), nac_3t2 (Debye-like background terms, P23) and yb_nor2 (a
//! texture of two axes, P30).

use drs::jacobian::Field;
use drs::powder::{BackgroundModel, Intensities, Model};
use std::path::Path;

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap().iter().map(|&b| b as char).collect()
}

fn model(run: &str, code: &str) -> Model {
    let d = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs").join(run);
    let mut p = drs::pcr::read(&latin1(&d.join(format!("{code}.pcr")))).unwrap();
    p.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
    // FullProf's Le Bail fit of the spinel made a Pawley fit (P36): its intensities are fields.
    if run == "MgAl2O4-LeBail" {
        let l = Model::from_pcr_with_data(&p, &d, code, Some(&d.join("MgAl2O4s.dat"))).unwrap();
        let mut m = l.model;
        drs::pawley::setup(&mut m, l.data.as_ref().unwrap(), &mut drs::refinable::Refinement::default(), 0, 0).unwrap();
        return m;
    }
    Model::from_pcr(&p, &d, code).unwrap_or_else(|e| panic!("{run}: {e}")).model
}

/// Each kind of field, numbered; no wildcard.
fn kind(f: Field) -> usize {
    use Field::*;
    match f {
        Zero => 0,
        Sycos => 1,
        Sysin => 2,
        Ratio => 3,
        Background(_) => 4,
        Scale(_) => 5,
        Bov(_) => 6,
        Cell(..) => 7,
        U(_) => 8,
        V(_) => 9,
        W(_) => 10,
        X(_) => 11,
        Y(_) => 12,
        GaussSize(_) => 13,
        Eta0(_) => 14,
        Asym(..) => 15,
        SL(_) => 16,
        DL(_) => 17,
        Pref(..) => 18,
        AtomX(..) => 19,
        AtomB(..) => 20,
        AtomOcc(..) => 21,
        AtomBeta(..) => 22,
        Mode(..) => 23,
        Quartic(..) => 24,
        QuarticShare(_) => 25,
        ShapeSize(_) => 26,
        Wavelength => 27,
        Micro(_) => 28,
        Intensity(..) => 29,
        Tof(_) => 30,
        TofProfile(..) => 31,
        Rigid(..) => 32,
        BackgroundTerm(_) => 33,
        Texture(..) => 34,
    }
}
const KINDS: usize = 35;

#[test]
fn every_field_described() {
    let runs = [("Y2O3", "y2o3_structure"), ("NdNi-LTOa", "NdNi_LTOa"), ("prba10k", "prba10k"), ("Pd3MnD", "pdfinal"), ("modes/A_modes", "a_modes"), ("Co-YBCO", "Cocus-asym"), ("LaB6", "LaB6rx"), ("NdNi-HTO", "NdNi_HTO"), ("MgAl2O4-LeBail", "MgAl2O4s_LeBail_conv"), ("arg_si", "arg_si"), ("../9-SAnn-OrganicCompounds/CitricAcid_example/Rietveld", "citric-acid-sls_1"), ("background_terms", "nac_3t2"), ("texture", "yb_nor2")];
    let mut seen = [0usize; KINDS];
    let mut shown = Vec::new();
    for (run, code) in runs {
        let m = model(run, code);
        for fi in m.fields() {
            let f = fi.field;
            let l = f.describe(&m);
            seen[kind(f)] += 1;
            assert_eq!(l.id, fi.id, "{run} {f:?}");
            assert!(drs::parameter_table::PARAMETERS.iter().any(|p| p.id == l.id && p.name == l.name), "{run} {f:?}: {} is not in the dictionary", l.id);
            assert_eq!(l.phase, f.phase(), "{run} {f:?}");
            // The atom: index and label as the structure has them, and the label in the address.
            let atom = match f {
                Field::AtomX(p, a, _) | Field::AtomB(p, a) | Field::AtomOcc(p, a) | Field::AtomBeta(p, a, _) => match &m.phases[p].intensities {
                    Intensities::Structure(st, _) => Some((a, st.atoms[a].label.clone())),
                    Intensities::List(_) => None,
                },
                _ => None,
            };
            assert_eq!(l.atom, atom, "{run} {f:?}");
            if let Some((_, label)) = &atom {
                assert!(fi.address.contains(&format!(".{label}.")) || fi.address.contains(&format!(".{label}:")), "{run}: {} for atom {label}", fi.address);
            }
            // The component.
            use Field::*;
            let want: Option<String> = match f {
                Background(k) => Some(match &m.background {
                    BackgroundModel::Points(p) | BackgroundModel::Spline(p) => format!("at {:.2}°", p[k].0),
                    _ => format!("order {k}"),
                }),
                Asym(_, k) => Some(["P1", "P2", "P3", "P4"][k].into()),
                AtomX(_, _, k) => Some(["x", "y", "z"][k].into()),
                AtomBeta(_, _, k) => Some(["β₁₁", "β₂₂", "β₃₃", "β₁₂", "β₁₃", "β₂₃"][k].into()),
                Mode(p, k) => Some(m.phases[p].modes.as_ref().unwrap().modes[k].name.clone()),
                // A group's origin or orientation by the group, an internal coordinate by its atom (P41).
                Rigid(p, g, k) => {
                    let gr = &m.phases[p].rigid[g];
                    let scope = if k < 6 { format!("group {}", g + 1) } else { gr.members[(k - 6) / 3].label.clone() };
                    Some(format!("{scope} {}", drs::rigid_body::value_component(gr, k)))
                }
                Quartic(p, k) => {
                    let t = m.phases[p].quartic_strain.as_ref().unwrap().terms[k].0;
                    Some(format!("S{}{}{}", t[0], t[1], t[2]))
                }
                // An axis of a texture by its number (P30).
                Texture(_, k) => Some(format!("axis {}", k / 2 + 1)),
                // A further term of the background by its number (P23).
                BackgroundTerm(k) => Some(format!("term {}", m.background_terms.id(k).1)),
                Intensity(p, k) => match &m.phases[p].intensities {
                    Intensities::List(l) => Some(format!("{} {} {}", l[k].0[0], l[k].0[1], l[k].0[2])),
                    Intensities::Structure(..) => None,
                },
                Zero | Sycos | Sysin | Ratio | Scale(_) | Bov(_) | Cell(..) | U(_) | V(_) | W(_) | X(_) | Y(_) | GaussSize(_) | Eta0(_) | SL(_) | DL(_) | Pref(..) | AtomB(..) | AtomOcc(..) | QuarticShare(_) | ShapeSize(_) | Wavelength | Micro(_) | Tof(_) | TofProfile(..) => None,
            };
            assert_eq!(l.component, want, "{run} {f:?}");
            if let Some(c) = &l.component {
                assert!(!c.is_empty(), "{run} {f:?}: an empty component");
            }
            if shown.len() < 40 && (l.component.is_some() || l.atom.is_some()) && !shown.iter().any(|s: &(usize, String)| s.0 == kind(f)) {
                shown.push((kind(f), format!("{run:14} {:30} {:22} {:?} {:?} {:?}", fi.address, l.name, l.phase, l.atom, l.component)));
            }
        }
    }
    for s in &shown {
        println!("{}", s.1);
    }
    let missing: Vec<usize> = (0..KINDS).filter(|&k| seen[k] == 0).collect();
    println!("{} fields described; per kind {seen:?}", seen.iter().sum::<usize>());
    assert!(missing.is_empty(), "kinds not met in the runs: {missing:?}");
}
