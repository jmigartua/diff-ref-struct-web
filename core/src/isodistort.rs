//! ISODISTORT's CIF of a distorted structure in terms of its symmetry modes (P13; Campbell et al. 2006, J. Appl.
//! Cryst. 39, 607–614; L263): the mode list with labels and amplitudes, the modes' normalisation factors, the
//! coordinate deltas and their formulas, and the sparse matrix that relates deltas and modes.
//!
//! The tags read are those GSAS-II's importer of ISODISTORT CIFs reads (GSASIIobj.py and G2phase_CIF.py, changeset
//! 5067): `_iso_displacivemode_ID`, `_label`, `_value`; `_iso_displacivemodenorm_ID`, `_value`;
//! `_iso_deltacoordinate_ID`, `_label`, `_value`; `_iso_coordinate_label`, `_formula`;
//! `_iso_displacivemodematrix_row`, `_col`, `_value`. The matrix gives the coordinate deltas from the modes, a row per
//! delta and a column per mode: deltas (lattice units) = matrix·(amplitudes × norms), as ISODISTORT's file says above
//! the loop and GSAS-II reads it (its `Mode2VarMatrix`, inverted for the modes). Checked on ISODISTORT 6.12.2's CIF of
//! case A (Sr₂ScSbO₆, Fm-3m → P2₁/n; `manuals/bcs/isodistort/`, L278): the deltas from the modes to the print, and
//! the irreps' and modes' amplitudes as the core's `mode_decomposition` computes them from the same structures (test
//! `isodistort`). L263 read the matrix the other way round, on a fixture. The mode labels are split as ISODISTORT
//! writes them, parent[k]irrep(direction)[atom:Wyckoff:kind]point-irrep(component),
//! e.g. "Fm-3m[0,0,0]GM4+(a,a,0)[O1:e:dsp]T1u(a)".

use crate::cif;

/// One mode of the file.
#[derive(Clone, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct IsoMode {
    pub id: usize,
    pub label: String,
    /// The parts of the label (empty when the label does not have them).
    pub parent: String,
    pub k: String,
    pub irrep: String,
    pub direction: String,
    pub atom: String,
    pub wyckoff: String,
    /// "dsp" (displacive), "occ", "mag" …
    pub kind: String,
    pub point_irrep: String,
    pub value: f64,
    /// `_iso_displacivemodenorm_value` of the mode, when the file gives it.
    pub norm: Option<f64>,
}

/// What an ISODISTORT CIF holds about its modes.
#[derive(Clone, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct IsoDistort {
    pub modes: Vec<IsoMode>,
    /// (ID, label, value) of each coordinate delta ("Sr1_dx").
    pub deltas: Vec<(usize, String, f64)>,
    /// (label, formula): "Sr1_x" = "1/4 + Sr1_dx".
    pub formulas: Vec<(String, String)>,
    /// (row, column, value) of the mode matrix.
    pub matrix: Vec<(usize, usize, f64)>,
}

/// The parts of an ISODISTORT mode label: parent, k, irrep, direction, atom, Wyckoff letter, kind, point irrep.
pub fn split_label(label: &str) -> Option<[String; 8]> {
    let (parent, rest) = label.split_once('[')?;
    let (k, rest) = rest.split_once(']')?;
    let (irrep, rest) = rest.split_once('(')?;
    let (direction, rest) = rest.split_once(')')?;
    let rest = rest.strip_prefix('[')?;
    let (site, rest) = rest.split_once(']')?;
    let mut sp = site.split(':');
    let (atom, wyckoff, kind) = (sp.next()?.to_string(), sp.next().unwrap_or("").to_string(), sp.next().unwrap_or("").to_string());
    let point = rest.split('(').next().unwrap_or("").to_string();
    Some([parent.to_string(), format!("[{k}]"), irrep.to_string(), format!("({direction})"), atom, wyckoff, kind, point])
}

/// Read the modes of an ISODISTORT CIF (the first block with `_iso_displacivemode_label`).
pub fn isodistort_cif(text: &str) -> Result<IsoDistort, String> {
    let blocks = cif::parse(text);
    let b = blocks.iter().find(|b| b.column("_iso_displacivemode_label").is_some()).ok_or("no _iso_displacivemode_label loop: not an ISODISTORT mode CIF")?;
    let col = |t: &str| b.column(t);
    let ids = |t: &str| -> Result<Vec<usize>, String> { col(t).map_or(Ok(Vec::new()), |c| c.iter().map(|s| s.parse::<usize>().map_err(|_| format!("{t}: {s} is not an ID"))).collect()) };
    let nums = |t: &str| -> Result<Vec<f64>, String> { col(t).map_or(Ok(Vec::new()), |c| c.iter().map(|s| cif::number(s).ok_or(format!("{t}: {s} is not a number"))).collect()) };
    let labels = col("_iso_displacivemode_label").unwrap();
    let mode_ids = ids("_iso_displacivemode_id")?;
    let values = nums("_iso_displacivemode_value")?;
    let norm_ids = ids("_iso_displacivemodenorm_id")?;
    let norms = nums("_iso_displacivemodenorm_value")?;
    let mut out = IsoDistort::default();
    for (i, l) in labels.iter().enumerate() {
        let id = mode_ids.get(i).copied().unwrap_or(i + 1);
        let mut m = IsoMode { id, label: l.to_string(), value: values.get(i).copied().unwrap_or(0.0), ..Default::default() };
        if let Some([p, k, ir, d, a, w, kind, pt]) = split_label(l) {
            (m.parent, m.k, m.irrep, m.direction, m.atom, m.wyckoff, m.kind, m.point_irrep) = (p, k, ir, d, a, w, kind, pt);
        }
        m.norm = norm_ids.iter().position(|&n| n == id).and_then(|j| norms.get(j).copied()).or_else(|| if norm_ids.is_empty() { norms.get(i).copied() } else { None });
        out.modes.push(m);
    }
    if let Some(dl) = col("_iso_deltacoordinate_label") {
        let did = ids("_iso_deltacoordinate_id")?;
        let dv = nums("_iso_deltacoordinate_value")?;
        out.deltas = dl.iter().enumerate().map(|(i, l)| (did.get(i).copied().unwrap_or(i + 1), l.to_string(), dv.get(i).copied().unwrap_or(0.0))).collect();
    }
    if let (Some(cl), Some(cf)) = (col("_iso_coordinate_label"), col("_iso_coordinate_formula")) {
        out.formulas = cl.iter().zip(cf.iter()).map(|(l, f)| (l.to_string(), f.to_string())).collect();
    }
    let (rows, cols, vals) = (ids("_iso_displacivemodematrix_row")?, ids("_iso_displacivemodematrix_col")?, nums("_iso_displacivemodematrix_value")?);
    if rows.len() != cols.len() || rows.len() != vals.len() {
        return Err("the mode matrix's row, column and value columns differ in length".into());
    }
    out.matrix = rows.into_iter().zip(cols).zip(vals).map(|((r, c), v)| (r, c, v)).collect();
    Ok(out)
}

impl IsoDistort {
    /// The irreps of the file in order of appearance, with their modes' indices.
    pub fn irreps(&self) -> Vec<(String, Vec<usize>)> {
        let mut v: Vec<(String, Vec<usize>)> = Vec::new();
        for (i, m) in self.modes.iter().enumerate() {
            let key = format!("{}{}", m.k, m.irrep);
            match v.iter_mut().find(|(k, _)| *k == key) {
                Some((_, l)) => l.push(i),
                None => v.push((key, vec![i])),
            }
        }
        v
    }

    /// The coordinate deltas from the modes: δ_r = Σ_c M_rc A_c N_c (row ID, δ), rows in order of appearance.
    pub fn deltas_from_modes(&self) -> Vec<(usize, f64)> {
        let mut rows: Vec<(usize, f64)> = Vec::new();
        for &(r, c, v) in &self.matrix {
            let a = self.modes.iter().find(|m| m.id == c).map_or(0.0, |m| m.value * m.norm.unwrap_or(1.0));
            match rows.iter_mut().find(|x| x.0 == r) {
                Some(x) => x.1 += v * a,
                None => rows.push((r, v * a)),
            }
        }
        rows
    }

    /// The modes' amplitudes from the deltas: M x = δ solved (Gaussian elimination with partial pivoting), A = x/N;
    /// (mode ID, amplitude) in the modes' order. An error when the matrix is not square in the modes and deltas or is
    /// singular.
    pub fn modes_from_deltas(&self) -> Result<Vec<(usize, f64)>, String> {
        let n = self.modes.len();
        if self.deltas.len() != n {
            return Err(format!("{} deltas for {n} modes: the matrix is not square", self.deltas.len()));
        }
        let col = |id: usize| self.modes.iter().position(|m| m.id == id);
        let row = |id: usize| self.deltas.iter().position(|d| d.0 == id);
        let mut a = vec![vec![0.0; n + 1]; n];
        for &(r, c, v) in &self.matrix {
            let (Some(i), Some(j)) = (row(r), col(c)) else { return Err(format!("the matrix names row {r} or column {c}, not in the file")) };
            a[i][j] = v;
        }
        for (i, d) in self.deltas.iter().enumerate() {
            a[i][n] = d.2;
        }
        for k in 0..n {
            let p = (k..n).max_by(|&x, &y| a[x][k].abs().total_cmp(&a[y][k].abs())).unwrap();
            if a[p][k].abs() < 1e-12 {
                return Err("the mode matrix is singular".into());
            }
            a.swap(k, p);
            let pivot = a[k].clone();
            for (i, row) in a.iter_mut().enumerate() {
                if i != k {
                    let f = row[k] / pivot[k];
                    row.iter_mut().zip(&pivot).skip(k).for_each(|(x, p)| *x -= f * p);
                }
            }
        }
        Ok(self.modes.iter().enumerate().map(|(j, m)| (m.id, a[j][n] / a[j][j] / m.norm.unwrap_or(1.0))).collect())
    }
}
