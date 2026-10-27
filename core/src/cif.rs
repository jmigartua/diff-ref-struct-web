//! A small CIF reader: enough of a crystal structure for structure factors (cell, space group, atom sites,
//! anisotropic displacements). Not a full CIF parser: one data block is read (the first with atom sites), a text field
//! is one value holding its lines, and numbers lose their standard uncertainties ("0.1234(5)" is 0.1234; `value_su`
//! keeps them).
//!
//! - Space group: the symmetry operators (`_space_group_symop_operation_xyz` or `_symmetry_equiv_pos_as_xyz`) when
//!   present, else the Hall symbol, else the Hermann–Mauguin symbol, as FullProf would read it (`from_fullprof_symbol`).
//! - Occupancy: `_atom_site_occupancy` is the chemical occupancy, converted to FullProf's Occ = occupancy × m/M.
//! - Displacements: `_atom_site_B_iso_or_equiv`, or `_atom_site_U_iso_or_equiv` (B = 8π²U); anisotropic U_ij or B_ij
//!   become FullProf's betas (β_ij = 2π² a*_i a*_j U_ij = a*_i a*_j B_ij / 4), with Biso then 0.
//! - Species: `_atom_site_type_symbol`, written as the tables write ions ("O2-" → "O-2", "Ti4+" → "Ti+4"), or the
//!   letters that start the label.

use crate::cell::Cell;
use crate::spacegroup::SpaceGroup;
use crate::structure_factor::{beta_from_u, biso_from_uiso, occupancy_from_chemical, Atom, Structure};
use crate::symop::SymOp;
use std::collections::HashMap;

/// The data items of a block: single values, and loops as (tags, rows).
#[derive(Debug, Default)]
pub struct Block {
    pub name: String,
    pub items: HashMap<String, String>,
    pub loops: Vec<(Vec<String>, Vec<Vec<String>>)>,
}

impl Block {
    /// A single value, or the first value of a looped tag.
    pub fn get(&self, tag: &str) -> Option<&str> {
        self.items.get(tag).map(String::as_str).or_else(|| self.column(tag).and_then(|c| c.first().copied()))
    }

    /// All values of a looped tag.
    pub fn column(&self, tag: &str) -> Option<Vec<&str>> {
        self.loops.iter().find_map(|(tags, rows)| {
            let i = tags.iter().position(|t| t == tag)?;
            Some(rows.iter().map(|r| r[i].as_str()).collect())
        })
    }
}

/// The tokens of a CIF text: quoted strings unquoted, comments dropped, a text field (between lines that start with
/// ';') one token holding its lines.
fn tokens(text: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut field: Option<Vec<&str>> = None;
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix(';') {
            match field.take() {
                Some(lines) => out.push(lines.join("\n")),
                None => field = Some(if rest.trim().is_empty() { Vec::new() } else { vec![rest] }),
            }
            continue;
        }
        if let Some(lines) = field.as_mut() {
            lines.push(line);
            continue;
        }
        let c: Vec<char> = line.chars().collect();
        let mut i = 0;
        while i < c.len() {
            if c[i].is_whitespace() {
                i += 1;
            } else if c[i] == '#' {
                break;
            } else if c[i] == '\'' || c[i] == '"' {
                let q = c[i];
                let mut j = i + 1;
                while j < c.len() && !(c[j] == q && (j + 1 == c.len() || c[j + 1].is_whitespace())) {
                    j += 1;
                }
                out.push(c[i + 1..j.min(c.len())].iter().collect());
                i = j + 1;
            } else {
                let mut j = i;
                while j < c.len() && !c[j].is_whitespace() {
                    j += 1;
                }
                out.push(c[i..j].iter().collect());
                i = j;
            }
        }
    }
    out
}

/// The data blocks of a CIF text.
pub fn parse(text: &str) -> Vec<Block> {
    let t = tokens(text);
    let mut blocks: Vec<Block> = Vec::new();
    let mut i = 0;
    let is_tag = |s: &str| s.starts_with('_');
    let is_key = |s: &str| {
        let l = s.to_ascii_lowercase();
        s.starts_with('_') || l == "loop_" || l.starts_with("data_") || l.starts_with("save_") || l == "global_"
    };
    while i < t.len() {
        let low = t[i].to_ascii_lowercase();
        if low.starts_with("data_") {
            blocks.push(Block { name: t[i][5..].to_string(), ..Default::default() });
            i += 1;
        } else if low == "loop_" {
            i += 1;
            let mut tags = Vec::new();
            while i < t.len() && is_tag(&t[i]) {
                tags.push(t[i].to_ascii_lowercase());
                i += 1;
            }
            let mut values = Vec::new();
            while i < t.len() && !is_key(&t[i]) {
                values.push(t[i].clone());
                i += 1;
            }
            if !tags.is_empty() {
                let rows = values.chunks(tags.len()).filter(|r| r.len() == tags.len()).map(|r| r.to_vec()).collect();
                if blocks.is_empty() {
                    blocks.push(Block::default());
                }
                blocks.last_mut().unwrap().loops.push((tags, rows));
            }
        } else if is_tag(&t[i]) {
            let tag = t[i].to_ascii_lowercase();
            let value = t.get(i + 1).filter(|v| !is_key(v)).cloned();
            if blocks.is_empty() {
                blocks.push(Block::default());
            }
            i += if value.is_some() { 2 } else { 1 };
            blocks.last_mut().unwrap().items.insert(tag, value.unwrap_or_else(|| "?".into()));
        } else {
            i += 1;
        }
    }
    blocks
}

/// A CIF number and its standard uncertainty ("0.1234(5)" is 0.1234 and 0.0005; "12(3)" 12 and 3); the uncertainty
/// is 0 when none is written. None for "?", "." and text.
pub fn value_su(s: &str) -> Option<(f64, f64)> {
    let v = number(s)?;
    let Some((a, b)) = s.split_once('(') else { return Some((v, 0.0)) };
    let digits: f64 = b.trim_end_matches(')').parse().ok()?;
    let decimals = a.split_once('.').map_or(0, |(_, d)| d.len());
    Some((v, digits * 10f64.powi(-(decimals as i32))))
}

/// A CIF number without its standard uncertainty; None for "?" and ".".
pub fn number(s: &str) -> Option<f64> {
    let s = s.split('(').next().unwrap_or(s);
    if s == "?" || s == "." {
        return None;
    }
    s.parse().ok()
}

/// A CIF type symbol as the tables write it: "O2-" → "O-2", "Ti4+" → "Ti+4", "Fe" → "Fe".
pub fn species_of(symbol: &str) -> String {
    let el: String = symbol.chars().take_while(|c| c.is_ascii_alphabetic()).collect();
    let rest = &symbol[el.len()..];
    let digits: String = rest.chars().filter(|c| c.is_ascii_digit()).collect();
    match rest.chars().find(|c| *c == '+' || *c == '-') {
        Some(sign) => format!("{el}{sign}{}", if digits.is_empty() { "1" } else { &digits }),
        None => el,
    }
}

/// The structure of the first block with atom sites.
pub fn structure(text: &str) -> Result<Structure, String> {
    let blocks = parse(text);
    let b = blocks
        .iter()
        .find(|b| b.column("_atom_site_fract_x").is_some())
        .ok_or("no block with atom sites (_atom_site_fract_x) in the CIF")?;
    let num = |tag: &str| -> Result<f64, String> {
        b.get(tag).and_then(number).ok_or(format!("{tag} missing or not a number"))
    };
    let cell = Cell::new(
        num("_cell_length_a")?,
        num("_cell_length_b")?,
        num("_cell_length_c")?,
        num("_cell_angle_alpha")?,
        num("_cell_angle_beta")?,
        num("_cell_angle_gamma")?,
    )
    .map_err(|e| e.to_string())?;
    let group = {
        let ops = b.column("_space_group_symop_operation_xyz").or_else(|| b.column("_symmetry_equiv_pos_as_xyz"));
        let hall = b.get("_space_group_name_hall").or_else(|| b.get("_symmetry_space_group_name_hall")).filter(|s| *s != "?");
        let hm = b.get("_space_group_name_h-m_alt").or_else(|| b.get("_symmetry_space_group_name_h-m")).filter(|s| *s != "?");
        if let Some(ops) = ops {
            let ops: Vec<SymOp> = ops.iter().map(|s| SymOp::parse_xyz(s)).collect::<Result<_, _>>()?;
            SpaceGroup::from_operators(&ops)?
        } else if let Some(h) = hall {
            SpaceGroup::from_hall(h)?
        } else if let Some(h) = hm {
            SpaceGroup::from_fullprof_symbol(h)?
        } else {
            return Err("no symmetry operators, Hall or Hermann–Mauguin symbol in the CIF".into());
        }
    };
    let col = |tag: &str| b.column(tag);
    let labels = col("_atom_site_label").ok_or("_atom_site_label missing")?;
    let xs = [col("_atom_site_fract_x"), col("_atom_site_fract_y"), col("_atom_site_fract_z")];
    let types = col("_atom_site_type_symbol");
    let occs = col("_atom_site_occupancy");
    let bs = col("_atom_site_b_iso_or_equiv");
    let us = col("_atom_site_u_iso_or_equiv");
    // Anisotropic displacements by label.
    let mut aniso: HashMap<String, [f64; 6]> = HashMap::new();
    if let Some(al) = col("_atom_site_aniso_label") {
        let names = ["11", "22", "33", "12", "13", "23"];
        for (kind, factor) in [("u", 1.0), ("b", 1.0 / (8.0 * std::f64::consts::PI * std::f64::consts::PI))] {
            let cols: Vec<Option<Vec<&str>>> = names.iter().map(|n| col(&format!("_atom_site_aniso_{kind}_{n}"))).collect();
            if cols.iter().all(Option::is_some) {
                for (i, l) in al.iter().enumerate() {
                    let u: Vec<f64> = cols.iter().map(|c| number(c.as_ref().unwrap()[i]).unwrap_or(0.0) * factor).collect();
                    aniso.insert(l.to_string(), beta_from_u(&cell, [u[0], u[1], u[2], u[3], u[4], u[5]]));
                }
            }
        }
    }
    let mut atoms = Vec::new();
    for (i, label) in labels.iter().enumerate() {
        let mut x = [0.0; 3];
        for k in 0..3 {
            x[k] = xs[k].as_ref().and_then(|c| number(c[i])).ok_or(format!("atom {label}: coordinate {k} missing"))?;
        }
        let species = match types.as_ref().map(|t| t[i]).filter(|t| *t != "?" && *t != ".") {
            Some(t) => species_of(t),
            None => label.chars().take_while(|c| c.is_ascii_alphabetic()).collect(),
        };
        let chemical = occs.as_ref().and_then(|o| number(o[i])).unwrap_or(1.0);
        let beta = aniso.get(*label).copied();
        let biso = if beta.is_some() {
            0.0
        } else {
            bs.as_ref()
                .and_then(|c| number(c[i]))
                .or_else(|| us.as_ref().and_then(|c| number(c[i])).map(biso_from_uiso))
                .unwrap_or(0.0)
        };
        let occ = occupancy_from_chemical(&group, x, chemical).map_err(|e| format!("atom {label}: {e}"))?;
        atoms.push(Atom { label: label.to_string(), species, x, biso, occ, beta, charge: None });
    }
    Ok(Structure { group, cell, atoms })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tokens_and_loops() {
        let text = "data_x\n_cell_length_a 5.1(2)\n_symmetry_space_group_name_H-M 'P n m a'\n;\ntext field\n;\nloop_\n_atom_site_label\n_atom_site_fract_x\nO1 0.25 # comment\nLa1 -0.05(2)\n";
        let b = &parse(text)[0];
        assert_eq!(b.get("_cell_length_a").and_then(number), Some(5.1));
        assert_eq!(b.get("_symmetry_space_group_name_h-m"), Some("P n m a"));
        assert_eq!(b.column("_atom_site_label").unwrap(), ["O1", "La1"]);
        assert_eq!(number(b.column("_atom_site_fract_x").unwrap()[1]), Some(-0.05));
    }

    #[test]
    fn ion_symbols() {
        assert_eq!(species_of("O2-"), "O-2");
        assert_eq!(species_of("Ti4+"), "Ti+4");
        assert_eq!(species_of("Cl-"), "Cl-1");
        assert_eq!(species_of("Pb"), "Pb");
    }
}
