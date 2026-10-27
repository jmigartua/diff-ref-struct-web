//! Readers of FullProf's output shared by the FullProf comparison tests.
#![allow(dead_code)]

use drs::jacobian::Field;
use std::collections::BTreeMap;
use std::path::Path;

pub fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display())).iter().map(|&b| b as char).collect()
}

/// FullProf's refined parameters in its last cycle: name → (value, esd).
pub fn fullprof_params(out: &str) -> BTreeMap<String, (f64, f64)> {
    let mut m = BTreeMap::new();
    for l in out.lines().filter(|l| l.contains("Parameter number") && l.contains("+/-")) {
        let Some((_, r)) = l.split_once(':') else { continue };
        let mut t = r.split_whitespace();
        let name = t.next().unwrap().to_string();
        let rest: String = t.collect::<Vec<_>>().join(" ");
        let rest = rest.replace(['(', ')'], " ").replace("+/-", " ");
        let v: Vec<f64> = rest.split_whitespace().filter_map(|s| s.parse().ok()).collect();
        if v.len() >= 2 {
            m.insert(name, (v[0], v[1]));
        }
    }
    m
}

/// A FullProf parameter name as a field of the model (FullProf's space), given the atoms' labels per phase.
pub fn fullprof_field(name: &str, labels: &[Vec<String>]) -> Option<Field> {
    let n = name.trim_end_matches("_pat1");
    let (base, phase) = match n.rfind("_ph") {
        Some(i) => (&n[..i], n[i + 3..].parse::<usize>().ok()? - 1),
        None => (n, 0),
    };
    let atom = |label: &str| labels.get(phase)?.iter().position(|l| l == label);
    Some(match base {
        "Zero" => Field::Zero,
        "SyCos" => Field::Sycos,
        "SySin" => Field::Sysin,
        "Lambda_rat" => Field::Ratio,
        "Scale" => Field::Scale(phase),
        "Bover" => Field::Bov(phase),
        "U-Cagl" => Field::U(phase),
        "V-Cagl" => Field::V(phase),
        "W-Cagl" => Field::W(phase),
        "X-tan" => Field::X(phase),
        "Y-cos" => Field::Y(phase),
        "G-Size" => Field::GaussSize(phase),
        "EtaPV" => Field::Eta0(phase),
        "SL_asym" => Field::SL(phase),
        "DL_asym" => Field::DL(phase),
        "Pref1" => Field::Pref(phase, 0),
        "Pref2" => Field::Pref(phase, 1),
        b if b.starts_with("Asym") => Field::Asym(phase, b[4..].parse::<usize>().ok()? - 1),
        b if b.starts_with("Cell_") => Field::Cell(phase, "ABCDEF".find(&b[5..])?),
        b if b.starts_with("Bck_") => Field::Background(b[4..].parse().ok()?),
        b => {
            let (kind, label) = b.split_once('_')?;
            match kind {
                "X" => Field::AtomX(phase, atom(label)?, 0),
                "Y" => Field::AtomX(phase, atom(label)?, 1),
                "Z" => Field::AtomX(phase, atom(label)?, 2),
                "Biso" => Field::AtomB(phase, atom(label)?),
                "Occ" => Field::AtomOcc(phase, atom(label)?),
                k if k.starts_with("bet") => Field::AtomBeta(phase, atom(label)?, ["11", "22", "33", "12", "13", "23"].iter().position(|s| *s == &k[3..])?),
                _ => return None,
            }
        }
    })
}

