//! What waits for FullProf's runs of it (L345) is refused, its line or its term named, rather than misread or
//! calculated without a reference (L037):
//! - P30 (L360): several preferred-orientation directions per phase (Nor = 2) and the multi-axial March–Dollase
//!   model (Nor = 3), read since (`texture`); a file with Nor = 2, 3 without their lines is refused at the missing
//!   line, and a Nor beyond FullProf's four functions (0–3) is refused;
//! - P23 (L353): the Fourier filtering of the background (Nba = −2), which FullProf's run did not reproduce; the
//!   Debye-like terms and the 12-term polynomial, calculated since, are read from the material's two files
//!   (`FP_tests/nac_3t2.pcr`, `FP_tests/ndni_htr.pcr` with its magnetic phase set aside);
//! - P41 (L367): a rigid body's TLS displacements (Furth ≠ 0), calculated since against FullProf's rigid/ runs (`tls`):
//!   `FP_tests/urea.pcr` is read, its group's atoms given their β by T, L and S.

use std::path::{Path, PathBuf};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof")
}

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display())).iter().map(|&b| b as char).collect()
}

/// Y₂O₃'s file with the line after the comment `after` changed by `edit`.
fn y2o3_with(after: &str, edit: impl Fn(&str) -> String) -> String {
    let text = latin1(&root().join("runs/Y2O3/y2o3_structure.pcr"));
    let mut lines: Vec<String> = text.lines().map(String::from).collect();
    let i = lines.iter().position(|l| l.starts_with(after)).unwrap() + 1;
    lines[i] = edit(&lines[i]);
    lines.join("\n")
}

/// The n-th field of a line set to `v`, the rest as written.
fn field(line: &str, n: usize, v: &str) -> String {
    let mut f: Vec<&str> = line.split_whitespace().collect();
    f[n] = v;
    format!("  {}", f.join("  "))
}

#[test]
fn several_preferred_orientations_without_their_lines() {
    // Y₂O₃'s phase has no axes after line 30-1: the reader names the line it expected (P30, L360).
    for nor in ["2", "3"] {
        let e = drs::pcr::read(&y2o3_with("!Job Npr", |l| field(l, 6, nor))).unwrap_err();
        assert!(e.message.contains("preferred-orientation direction"), "Nor = {nor}: {}", e.message);
    }
    let e = drs::pcr::read(&y2o3_with("!Job Npr", |l| field(l, 6, "4"))).unwrap_err();
    assert!(e.message.contains("Nor = 4"), "{}", e.message);
    assert!(drs::pcr::read(&y2o3_with("!Job Npr", |l| field(l, 6, "1"))).is_ok());
}

#[test]
fn background_terms() {
    // The two files of the material with Nba = −1 and −3 (L345): read with their terms (P23, L353).
    for (file, nba, debye, inverse) in [("FP_tests/nac_3t2.pcr", -1, 6, 0), ("FP_tests/ndni_htr.pcr", -3, 0, 3)] {
        let path = root().join(file);
        let mut p = drs::pcr::read(&latin1(&path)).unwrap();
        assert_eq!(p.flag("Nba"), nba);
        // ndni_htr's second phase is magnetic (Jbt = 1), refused before the background: its nuclear phase alone.
        p.phases.retain(|ph| ph.flag("Jbt") != 1);
        let m = drs::powder::Model::from_pcr(&p, path.parent().unwrap(), path.file_stem().unwrap().to_str().unwrap()).unwrap_or_else(|e| panic!("{file}: {e}")).model;
        assert_eq!((m.background_terms.debye.len(), m.background_terms.inverse.len()), (debye, inverse), "{file}");
    }
    // Nba = −2, the Fourier filtering, stays refused: `background_terms::fourier_filter_refused`.
}

#[test]
fn tls_read() {
    let path = root().join("FP_tests/urea.pcr");
    let p = drs::pcr::read(&latin1(&path)).unwrap();
    let m = drs::powder::Model::from_pcr(&p, path.parent().unwrap(), "urea").unwrap_or_else(|e| panic!("{e}")).model;
    let g = &m.phases[0].rigid[0];
    assert!(g.tls() && g.displacements.is_some());
    let drs::powder::Intensities::Structure(st, _) = &m.phases[0].intensities else { panic!() };
    let members = drs::rigid_body::tls_members(&m.phases[0].rigid, 0);
    assert!(members.iter().all(|&(h, k)| st.atoms[m.phases[0].rigid[h].members[k].index].beta.is_some()), "every TLS atom has its β");
}
