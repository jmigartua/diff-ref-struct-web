//! FullProf's structure file, the CFL (CrysFML's "crystallographic file"): what FullProf writes for Bond_Str
//! (`CODFIL_n.cfl`, L081) and Bond_Str reads (its manual, `manuals/fullprof/docs/bond_str.pdf`, pp. 3–4). Keywords,
//! one per line, case insensitive, "!" a comment:
//! - `Title` text; `Cell` a b c α β γ; `SpGR` (or `Spg`) the space-group symbol;
//! - `Atom` label species x y z Biso Occ [2·spin, charge], Occ in FullProf's units (occ = m/M for a full site);
//! - Bond_Str's instructions: `DMAX` d_distances d_angles, `TOL` (%), and others kept as they are.
//!
//! A number may carry its standard deviation in parentheses, in units of its last digit ("0.0774(6)" is 0.0774 ±
//! 0.0006); they are kept in `Uncertainties`.

use crate::bonds::Uncertainties;
use crate::cell::Cell;
use crate::spacegroup::SpaceGroup;
use crate::structure_factor::{Atom, Structure};

/// A CFL file read.
#[derive(Clone, Debug)]
pub struct Cfl {
    pub title: String,
    pub structure: Structure,
    pub uncertainties: Uncertainties,
    /// `DMAX`: the largest distance listed and the largest bond length for angles (Å); 0 when not given.
    pub dmax: (f64, f64),
    /// `TOL`: the tolerance on the sum of ionic radii, in %; None when not given.
    pub tol: Option<f64>,
}

/// A number and its standard deviation ("0.0774(6)" → (0.0774, 0.0006); "90.0" → (90.0, 0.0)).
pub fn number(s: &str) -> Option<(f64, f64)> {
    let s = s.trim();
    let (v, e) = match s.find('(') {
        Some(i) => (&s[..i], s[i + 1..].trim_end_matches(')')),
        None => (s, ""),
    };
    let value: f64 = v.parse().ok()?;
    if e.is_empty() {
        return Some((value, 0.0));
    }
    let digits = v.split(['e', 'E']).next()?.split('.').nth(1).map_or(0, str::len) as i32;
    let esd: f64 = e.parse().ok()?;
    Some((value, esd * 10f64.powi(-digits)))
}

/// Read a CFL text.
pub fn read(text: &str) -> Result<Cfl, String> {
    let (mut title, mut cell, mut group, mut atoms, mut ux, mut dmax, mut tol) = (String::new(), None, None, Vec::new(), Vec::new(), (0.0, 0.0), None);
    for (n, line) in text.lines().enumerate() {
        let line = line.split('!').next().unwrap_or("").trim();
        let mut it = line.split_whitespace();
        let Some(key) = it.next() else { continue };
        let rest: Vec<&str> = it.collect();
        let err = |m: &str| format!("line {}: {m}", n + 1);
        let nums = |k: usize| -> Result<Vec<(f64, f64)>, String> { rest.iter().take(k).map(|t| number(t).ok_or_else(|| err(&format!("{t} is not a number")))).collect() };
        match key.to_lowercase().as_str() {
            "title" => title = rest.join(" "),
            "cell" => {
                let v = nums(6)?;
                if v.len() < 6 {
                    return Err(err("Cell needs a b c alpha beta gamma"));
                }
                let c = Cell::new(v[0].0, v[1].0, v[2].0, v[3].0, v[4].0, v[5].0).map_err(|e| err(&format!("{e:?}")))?;
                cell = Some((c, std::array::from_fn::<f64, 6, _>(|k| v[k].1)));
            }
            "spgr" | "spg" | "space_group" => group = Some(SpaceGroup::from_fullprof_symbol(&rest.join(" ")).map_err(|e| err(&e))?),
            "atom" => {
                if rest.len() < 7 {
                    return Err(err("Atom needs label species x y z Biso Occ"));
                }
                let v: Vec<(f64, f64)> = rest[2..7].iter().map(|t| number(t).ok_or_else(|| err(&format!("{t} is not a number")))).collect::<Result<_, _>>()?;
                atoms.push(Atom { label: rest[0].to_string(), species: rest[1].to_string(), x: [v[0].0, v[1].0, v[2].0], biso: v[3].0, occ: v[4].0, beta: None, charge: None });
                ux.push([v[0].1, v[1].1, v[2].1]);
            }
            "dmax" => {
                let v = nums(2)?;
                dmax = (v.first().map_or(0.0, |x| x.0), v.get(1).map_or(0.0, |x| x.0));
            }
            "tol" => tol = nums(1)?.first().map(|x| x.0),
            _ => {}
        }
    }
    let (cell, ucell) = cell.ok_or("no Cell")?;
    let group = group.ok_or("no SpGR")?;
    Ok(Cfl { title, structure: Structure { group, cell, atoms }, uncertainties: Uncertainties { cell: ucell, x: ux }, dmax, tol })
}

#[cfg(test)]
mod tests {
    use super::number;

    #[test]
    fn numbers_with_deviations() {
        assert_eq!(number("0.0774(6)"), Some((0.0774, 0.0006000000000000001)));
        assert_eq!(number("5.54018(7)").map(|v| (v.0, (v.1 * 1e5).round())), Some((5.54018, 7.0)));
        assert_eq!(number("0.355(15)").map(|v| (v.1 * 1e3).round()), Some(15.0));
        assert_eq!(number("90.0000"), Some((90.0, 0.0)));
        assert_eq!(number("12(3)"), Some((12.0, 3.0)));
    }
}
