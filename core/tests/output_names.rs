//! The parameter a line of an output names (the app's request 36): FullProf's symbolic names and the CIF's tags.
//!
//! - **FullProf:** every `.out` of `manuals/` with its `.pcr` beside it, whose model the core builds: each line
//!   "Parameter number n -> Symbolic Name: NAME" of its first list, NAME's field (`output_names::fullprof_field`) is
//!   one of the fields the parameter number n carries in the model built from the `.pcr` (`pcr_parameter_numbers`,
//!   the core's reading of the codewords): 3527 of 3563 names, none mapped to another number's field (the joint runs'
//!   each in its pattern's model, `Joint::fullprof_field`). The 36 others
//!   name what the core does not refine (λ₂'s own U, V, W; a reflection's own broadening, DHG2 and DHLo; Or1) or
//!   have a number the `.pcr` does not carry; a rigid group's ThRGB, PhRGB, ChRGB are its θ, φ, χ (P41, L347). In 15 files a name's number is not among the `.pcr`'s codewords: `NOT_OF_THE_PCR`, with the reason.
//! - **CIF:** the results file the core writes for two runs: every value of `_cell_`, the `_atom_site_` loop and
//!   the `_atom_site_aniso_` loop written with a standard uncertainty names a refined parameter of the run.

use drs::output_names::{address_of_cif, fullprof_field, fullprof_fields};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

mod common;
use common::latin1;

fn outs(d: &Path, v: &mut Vec<PathBuf>) {
    for e in std::fs::read_dir(d).unwrap() {
        let p = e.unwrap().path();
        if p.is_dir() {
            outs(&p, v);
        } else if p.extension().is_some_and(|x| x == "out") {
            v.push(p);
        }
    }
}

/// Counts as measured in L267 (with prba10k, cuf1k, lac and Pd3MnD built): (`.out` files with a parameter list, their models built, names compared, names that
/// map to their number's field). L322 added the six `.out` of the joint runs pbb, scor_rxn and prbaco, each name in its
/// pattern's model (`Joint::fullprof_field`): 276 names, all mapped. L341 added P05's TOF runs (23 `.out` with a list,
/// 21 built: arg_si_plain's .pcr holds NaN, arg_si_abs2b's ABSCOR2 is refused): 515 names, all mapped (Dtt1, Dtt2,
/// Dtt_1overD, Sig-n, Gam-n, Alpha0/1/Q, Beta0/1/Q, ExtTOF, TOF_Abs1; Dtt_1overD is Zerot with Hodges' terms).
/// L353 added L346's runs (40 `.out` with a list, 19 built: TLS and Nor ≥ 2 wait for their steps, y2o3_filter is
/// refused): 299 names, 294 mapped (the Debye-like amplitudes Bck_12 … Bck_16 among them); the rest in nac_3t2's three
/// `.out`, whose bet23_Ca1 the `.pcr` does not code (`NOT_OF_THE_PCR`).
/// L354 added L351's citric runs (4 `.out` with a list, 4 built; citric_free's too): 68 names, 64 mapped.
/// L360 builds texture/'s four (Nor = 2, 3, P30): 96 names, all mapped.
/// L367 builds rigid/'s TLS runs (16 more): 328 names, all mapped, with FullProf 8.20's t11 … s32 (T11 … in pdd's),
/// X_RGB, Y_RGB, Z_RGB (a group's origin by its first atom), Dm_In, TheIn, PhiIn (an atom's internal coordinates) and
/// Or1 (Pref1); with them 7 names of earlier runs map too (38 left unmapped: DHG2, DHLo, U-, V-, W-Lam2).
const COUNTS: (usize, usize, usize, usize) = (205, 181, 4354, 4316);

/// `.out` files whose names map to fields but whose parameter numbers the `.pcr` beside them does not carry, with the
/// reason.
const NOT_OF_THE_PCR: [(&str, &str); 15] = [
    ("fullprof/9-SAnn-OrganicCompounds/CitricAcid_example/PM/citric-acid-sls_1.out", "the .out of another version of the .pcr: its numbers 1–9 are none of the .pcr's codewords"),
    ("fullprof/runs/LaB6/LaB6rx.out", "Lambda_rat, FullProf's number 4, is the .pcr's VARY ratio, which numbers no codeword (L071)"),
    ("fullprof/runs/LaB6/LaB6rx_calc.out", "as LaB6rx.out"),
    ("fullprof/runs/LaB6/LaB6rx_conv.out", "as LaB6rx.out"),
    ("fullprof/runs/YVO3/YVO3a.out", "the .out's run refined U as parameter 32 (it ends at −0.091 with an esd of 3·10⁻¹⁰); the .pcr beside it does not refine U"),
    ("fullprof/runs/YVO3/YVO3a_conv.out", "as YVO3a.out"),
    ("fullprof/runs/background_terms/nac_3t2.out", "bet23_Ca1, FullProf's number 51, FullProf freed itself: the .pcr's codeword for Ca1's β23 is 0 (L353)"),
    ("fullprof/runs/background_terms/nac_3t2_calc.out", "as nac_3t2.out"),
    ("fullprof/runs/background_terms/nac_3t2_debye.out", "as nac_3t2.out"),
    ("fullprof/runs/window_experiments/yvo3_X0/YVO3a.out", "as YVO3a.out"),
    ("fullprof/runs/window_experiments/yvo3_base_w4/YVO3a.out", "as YVO3a.out"),
    ("fullprof/runs/window_experiments/yvo3_noK2/YVO3a.out", "as YVO3a.out"),
    ("fullprof/runs/window_experiments/yvo3_noK2_noasym/YVO3a.out", "as YVO3a.out"),
    ("fullprof/runs/window_experiments/yvo3_noXins/YVO3a.out", "as YVO3a.out"),
    ("fullprof/runs/window_experiments/yvo3_noasym/YVO3a.out", "as YVO3a.out"),
];

/// `.out` files where a name's atom label is shared by several atoms, so that the name is mapped by its number among
/// `fullprof_fields`: Pd₃MnD's antisite pairs (two MN and two PD atoms; Biso_MN_ph1 is parameters 4 and 6, L267).
const SHARED_LABELS: [&str; 2] = ["fullprof/runs/Pd3MnD/pdfinal.out", "fullprof/runs/Pd3MnD/pdfinal_calc.out"];

#[test]
fn fullprof_symbolic_names() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals");
    let mut files = Vec::new();
    outs(&root, &mut files);
    files.sort();
    let (mut with_list, mut built, mut compared, mut mapped) = (0, 0, 0, 0);
    let mut unmapped: BTreeMap<String, usize> = BTreeMap::new();
    let mut wrong = Vec::new();
    let mut apart = std::collections::BTreeSet::new();
    let mut shared = std::collections::BTreeSet::new();
    for out in files {
        let text = latin1(&out);
        // The first list: "Parameter number n -> Symbolic Name: NAME value".
        let list: Vec<(u32, String)> = text
            .lines()
            .filter_map(|l| {
                let (a, b) = l.split_once("-> Symbolic Name:")?;
                let n = a.split("Parameter number").nth(1)?.split_whitespace().next()?.parse().ok()?;
                Some((n, b.split_whitespace().next()?.to_string()))
            })
            .collect();
        if list.is_empty() {
            continue;
        }
        with_list += 1;
        let pcr = out.with_extension("pcr");
        // A joint run (P06, L322): each name in its pattern's model (`Joint::fullprof_field`), against the numbers that
        // pattern's codewords carry.
        let pcr_text = std::fs::read(&pcr).map(|b| b.iter().map(|&c| c as char).collect::<String>()).unwrap_or_default();
        if drs::pcr::is_multipattern(&pcr_text) {
            let stem = pcr.file_stem().unwrap().to_str().unwrap();
            let Ok(j) = drs::joint::Joint::from_pcr(&pcr_text, out.parent().unwrap(), stem, false) else { continue };
            let pcrs = drs::pcr::read_all(&pcr_text).unwrap();
            let numbers: Vec<_> = pcrs.iter().zip(&j.patterns).map(|(p, m)| drs::refinable::pcr_parameter_numbers(p, &m.model).unwrap()).collect();
            built += 1;
            for (n, name) in &list {
                compared += 1;
                match j.fullprof_field(name) {
                    Some((k, f)) => match numbers[k].get(n) {
                        Some(fields) if fields.contains(&f) => mapped += 1,
                        Some(fields) => wrong.push(format!("{}: {n} {name} → pattern {} {f:?}, the codewords' {fields:?}", out.strip_prefix(&root).unwrap().display(), k + 1)),
                        None => {
                            apart.insert(out.strip_prefix(&root).unwrap().display().to_string());
                        }
                    },
                    None => {
                        let kind = name.split("_ph").next().unwrap().split("_pat").next().unwrap();
                        let kind = kind.split_once('_').map_or(kind, |(k, _)| k);
                        *unmapped.entry(kind.to_string()).or_default() += 1;
                    }
                }
            }
            continue;
        }
        let Ok(p) = std::fs::read(&pcr).map(|b| b.iter().map(|&c| c as char).collect::<String>()).map_err(|_| ()).and_then(|t| drs::pcr::read(&t).map_err(|_| ())) else { continue };
        let stem = out.file_stem().unwrap().to_str().unwrap();
        let Ok(l) = drs::powder::Model::from_pcr(&p, out.parent().unwrap(), stem) else { continue };
        let Ok(numbers) = drs::refinable::pcr_parameter_numbers(&p, &l.model) else { continue };
        built += 1;
        for (n, name) in &list {
            compared += 1;
            match fullprof_field(name, &l.model) {
                Some(f) => match numbers.get(n) {
                    Some(fields) if fields.contains(&f) => mapped += 1,
                    // A label shared by several atoms: the name names each of them (`fullprof_fields`).
                    Some(fields) if fullprof_fields(name, &l.model).iter().any(|g| fields.contains(g)) => {
                        mapped += 1;
                        shared.insert(out.strip_prefix(&root).unwrap().display().to_string());
                    }
                    Some(fields) => wrong.push(format!("{}: {n} {name} → {f:?}, the codewords' {fields:?}", out.strip_prefix(&root).unwrap().display())),
                    None => {
                        apart.insert(out.strip_prefix(&root).unwrap().display().to_string());
                    }
                },
                None => {
                    let kind = name.split("_ph").next().unwrap().split("_pat").next().unwrap();
                    let kind = kind.split_once('_').map_or(kind, |(k, _)| k);
                    *unmapped.entry(kind.to_string()).or_default() += 1;
                }
            }
        }
    }
    println!("{with_list} .out with a parameter list, {built} models built: {compared} names, {mapped} mapped to their number's field; unmapped by kind: {unmapped:?}");
    for w in wrong.iter().take(20) {
        println!("  {w}");
    }
    println!("names whose number the .pcr's codewords do not carry, in: {apart:?}");
    assert!(wrong.is_empty(), "{} names map to another field than their number's", wrong.len());
    assert_eq!(apart.iter().cloned().collect::<Vec<_>>(), NOT_OF_THE_PCR.iter().map(|x| x.0.to_string()).collect::<Vec<_>>());
    assert_eq!(shared.iter().cloned().collect::<Vec<_>>(), SHARED_LABELS.to_vec(), "the files where a name's label is shared");
    assert_eq!((with_list, built, compared, mapped), COUNTS);
}

/// The core's own results file: every cell, atom-site and aniso value with an su names a refined parameter.
#[test]
fn cif_tags() {
    let mut total = 0;
    for (dir, code) in [("runs/Y2O3", "y2o3_structure"), ("runs/rutile-anatase", "rutana_tch")] {
        let d = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof").join(dir);
        let p = drs::pcr::read(&latin1(&d.join(format!("{code}.pcr")))).unwrap();
        let c = drs::control::Control::from_pcr(&p, &d, code, &d).unwrap();
        let data = drs::powder::Model::from_pcr(&p, &d, code).unwrap().data.unwrap();
        let mut m = c.model.clone();
        let set = c.refinement.build(&m).unwrap();
        let r = drs::lsq::refine(&mut m, &data, &set, drs::lsq::Options { max_cycles: 1, eps: 0.1 }, |_| {});
        let mut cc = c.clone();
        cc.model = m.clone();
        let (text, _) = drs::results::write(&cc, &data, Some(&set), Some(&r), "run").unwrap();
        let refined: Vec<String> = set.free.iter().flat_map(|f| f.fields.iter().map(|(fl, _)| drs::refinable::address(*fl, &m))).collect();
        // Walk the phase blocks: tags with an su, and the atom loops' rows.
        let (mut phase, mut checked) = (None::<usize>, 0);
        let mut lines = text.lines().peekable();
        while let Some(l) = lines.next() {
            // A phase's block: its id ends in "|phase <n>".
            if l.starts_with("_pd_block_id") {
                phase = l.split("|phase ").nth(1).and_then(|r| r.trim_end_matches('\'').trim().parse::<usize>().ok()).map(|n| n - 1);
            }
            let Some(k) = phase else { continue };
            // The cell's six parameters (the volume is derived from them).
            if (l.starts_with("_cell_length_") || l.starts_with("_cell_angle_")) && l.contains('(') {
                let tag = l.split_whitespace().next().unwrap();
                let a = address_of_cif(tag, k, None, &m).unwrap_or_else(|| panic!("{code}: {tag}"));
                assert!(refined.contains(&a), "{code}: {tag} → {a}, not refined");
                checked += 1;
            }
            if l.trim() == "loop_" && lines.peek().is_some_and(|n| n.starts_with("_atom_site_label") || n.starts_with("_atom_site_aniso_label")) {
                let mut cols = Vec::new();
                while lines.peek().is_some_and(|n| n.starts_with('_')) {
                    cols.push(lines.next().unwrap().trim().to_string());
                }
                while lines.peek().is_some_and(|n| !n.trim().is_empty() && !n.starts_with('_') && n.trim() != "loop_") {
                    let row: Vec<&str> = lines.next().unwrap().split_whitespace().collect();
                    for (col, v) in cols.iter().zip(&row).skip(1) {
                        if v.contains('(') {
                            let a = address_of_cif(col, k, Some(row[0]), &m).unwrap_or_else(|| panic!("{code}: {col} of {}", row[0]));
                            assert!(refined.contains(&a), "{code}: {col} of {} → {a}, not refined", row[0]);
                            checked += 1;
                        }
                    }
                }
            }
        }
        println!("{code}: {checked} values with an su, each a refined parameter");
        assert!(checked > 0, "{code}: {checked}");
        total += checked;
    }
    // Y₂O₃'s .pcr refines its cell, not its atoms; rutile + anatase refines both.
    assert!(total > 5, "{total}");
    // A tag of no parameter, and an atom the model lacks.
    let d = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs/Y2O3");
    let p = drs::pcr::read(&latin1(&d.join("y2o3_structure.pcr"))).unwrap();
    let m = drs::powder::Model::from_pcr(&p, &d, "y2o3_structure").unwrap().model;
    assert_eq!(address_of_cif("_cell_volume", 0, None, &m), None);
    assert_eq!(address_of_cif("_atom_site_fract_x", 0, Some("Q9"), &m), None);
    assert_eq!(address_of_cif("_CELL_LENGTH_A", 0, None, &m).as_deref(), Some("phase1.cell_a"));
}
