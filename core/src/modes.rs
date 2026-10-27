//! Refinement by symmetry modes (R29, P13; L090): a distorted phase described as a reference structure plus the
//! amplitudes of its symmetry-adapted distortion modes.
//!
//! The model (FullProf's notes of 29 August 2008 – 16 January 2009, `manuals/fullprof/7-FullProf-Tutorial-SymmetryModes/
//! FullProf-Tutorial-Symmetry-Modes.pdf`, pp. 3–5; the course's `90-fullprof-pcr-anatomy.qmd`): the reference
//! structure is the parent (high-symmetry) structure written in the cell and the space group of the distorted phase
//! (AMPLIMODES's "Transformed high symmetry structure in the subgroup basis"); atom j of the distorted phase's
//! asymmetric unit is at
//!
//!   r_j = r⁰_j + Σ_m Q_m · c_m · V_mj,
//!
//! Q_m the amplitude of mode m (Å), V_mj its polarisation vector on atom j in fractions of the phase's cell, c_m a
//! coefficient (1 in every file of the material). The vectors are fixed by group theory; the amplitudes are what is
//! refined, in place of the coordinates.
//!
//! The normalisation (AMPLIMODES; Orayech et al. 2015, J. Appl. Cryst. 48, text of `manuals/bcs/
//! mode-crystallography-manual/articles/nb5132.md`, line 255): each mode's displacements, summed in square over the
//! atoms of the distorted structure's primitive cell, make 1 Å² in the metric of the reference (parent) cell
//! transformed to the distorted phase's basis. So a mode's amplitude is the root-mean-square size of its distortion
//! over the primitive cell, in Å, and amplitudes scale with the reference metric: the same structure described from a
//! parent of edge a₁ or a₂ has amplitudes in the ratio a₂/a₁ (the course's case B). `gram` computes the sums, `project`
//! inverts them: the amplitudes of a structure given by its coordinates.
//!
//! The modes of one irreducible representation (irrep) together make its distortion: amplitude A = √(Σ_m Q_m²) over
//! the irrep's modes and polarisation vector (Q_m / A), normalised (Orayech et al. 2015, line 262: "the global
//! amplitude Q is separated from the internal form of the distortion, given by the normalized mode"; AMPLIMODES's
//! "Summary of Amplitudes"; FullProf's "Global Amplitude of Representation" in its `.out`). FullProf's standard
//! deviations of these amplitudes and of the positions neglect the correlations between amplitudes (its note of 29
//! August 2008); `irrep_esd` and `position_esd` take a covariance and give FullProf's with a diagonal one.
//!
//! `read_amplimodes` reads AMPLIMODES's text output (Bilbao Crystallographic Server), as pasted in the course's
//! `04_tools_materials/0_1…0_4`: each analysis's parent and distorted structures, transformation, reference structure,
//! atom mappings, summary of amplitudes and the block for FullProf.

use crate::cell::Cell;
use crate::mat3::{self, Mat3};
use crate::site::Site;
use crate::spacegroup::SpaceGroup;

/// The polarisation vector of a mode on one atom.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ModeVector {
    /// The atom of the phase (its index in the structure).
    pub atom: usize,
    /// The vector, in fractions of the phase's cell edges.
    pub v: [f64; 3],
    /// FullProf's coefficient c_m (1 unless the user ties modes with known ratios).
    pub coeff: f64,
}

/// One symmetry-adapted mode.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Mode {
    /// Its name, as AMPLIMODES writes it: "A3_GM4+" (number, irrep).
    pub name: String,
    /// The irreducible representation of the parent's space group it transforms by ("GM4+", "X5+").
    pub irrep: String,
    /// The amplitude Q (Å).
    pub amplitude: f64,
    pub vectors: Vec<ModeVector>,
}

/// The modes of a phase and the reference structure they displace.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Modes {
    /// The reference position of each atom of the phase (fractions of the phase's cell).
    pub reference: Vec<[f64; 3]>,
    pub modes: Vec<Mode>,
    /// The reference (parent) cell in the phase's basis, the metric of the modes' normalisation, when it is known
    /// (AMPLIMODES's reference structure; a `.pcr` does not give it).
    pub reference_cell: Option<Cell>,
    /// FullProf's box |Q| ≤ value on every amplitude (`Max_Amplitude`, Å).
    pub max_amplitude: Option<f64>,
    /// FullProf's FullProf Studio output integers (the A_MODES line after the count), kept for the export.
    pub fst: Vec<i32>,
    /// The parent's space group and the transformation to this phase's basis, when known (an AMPLIMODES analysis,
    /// `Analysis::mode_parent`, or the control file; a `.pcr` does not give it): the reflections' parent indices and
    /// k-vectors (`ModeParent::reflection`).
    #[cfg_attr(feature = "serde", serde(default))]
    pub parent: Option<ModeParent>,
    /// The sign reference (E14b, Josu's answer 6 to Explore's §8): the amplitudes by mode name of the analysis the
    /// modes come from, against which each state's equivalent set is chosen and named (`reference_signs`).
    #[cfg_attr(feature = "serde", serde(default))]
    pub sign_reference: Option<SignReference>,
}

/// Where a mode phase's sign reference comes from and its amplitudes by mode name (E14b): the control file's
/// `mode_reference = { file, analysis }`. The file is an AMPLIMODES output (analysis k's block for FullProf, whose
/// amplitudes are in the setting of that analysis) or FullProf's `.sum` (its "Amplitudes of symmetry modes", in the
/// phase's own basis).
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SignReference {
    pub file: String,
    pub analysis: Option<usize>,
    pub amplitudes: Vec<(String, f64)>,
}

impl SignReference {
    /// The reference read from its file's text, for phase `phase` (from 0) of a `.sum`; with an AMPLIMODES analysis,
    /// also the parent it states (`Analysis::mode_parent`), which a `.pcr` never gives.
    pub fn read(file: &str, text: &str, analysis: Option<usize>, phase: usize) -> Result<(SignReference, Option<ModeParent>), String> {
        let analyses = read_amplimodes(text);
        if !analyses.is_empty() {
            let k = analysis.unwrap_or(1);
            let a = analyses.get(k.wrapping_sub(1)).ok_or(format!("{file}: no analysis {k} ({} in the file)", analyses.len()))?;
            let (_, _, amps) = a.fullprof.as_ref().ok_or(format!("{file}: analysis {k} has no block for FullProf, so no amplitudes by mode name"))?;
            return Ok((SignReference { file: file.to_string(), analysis: Some(k), amplitudes: amps.clone() }, a.mode_parent()));
        }
        let blocks = read_fullprof_amplitudes(text);
        if !blocks.is_empty() {
            let b = blocks.get(phase).ok_or(format!("{file}: no amplitudes for phase {}", phase + 1))?;
            return Ok((SignReference { file: file.to_string(), analysis, amplitudes: b.modes.iter().map(|m| (m.0.clone(), m.1)).collect() }, None));
        }
        Err(format!("{file}: neither an AMPLIMODES output nor FullProf's .sum with the amplitudes of symmetry modes"))
    }
}

/// The amplitudes FullProf reports for a mode phase (E14c): its `.sum` block "Amplitudes of symmetry modes" (name,
/// value, sigma) and its lines "Global Amplitude of Representation" (irrep, value, sigma as printed in parentheses;
/// none where it prints no sigma).
#[derive(Clone, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct FullProfAmplitudes {
    /// FullProf's phase number (from 1), from the "=> Phase No." line before the block.
    pub phase: usize,
    pub modes: Vec<(String, f64, f64)>,
    pub irreps: Vec<(String, f64, Option<f64>)>,
}

/// A number printed with its sigma in parentheses, "0.283(6)": (0.283, Some(0.006)); "0.00000": (0, None).
fn with_sigma(s: &str) -> Option<(f64, Option<f64>)> {
    match s.split_once('(') {
        Some((v, e)) => {
            let e = e.trim_end_matches(')');
            let decimals = v.split_once('.').map_or(0, |(_, d)| d.len());
            Some((v.parse().ok()?, Some(e.parse::<f64>().ok()? * 10f64.powi(-(decimals as i32)))))
        }
        None => Some((s.parse().ok()?, None)),
    }
}

/// FullProf's reported amplitudes from its `.sum` (or `.out`), one entry per mode phase in the order the file first
/// names it: the block "Amplitudes of symmetry modes" after its "=> Phase No." line and the global amplitudes of
/// representation that follow (E14c). A phase printed twice (`a_modes_calc.sum` repeats its results) keeps the last.
pub fn read_fullprof_amplitudes(text: &str) -> Vec<FullProfAmplitudes> {
    let mut out: Vec<FullProfAmplitudes> = Vec::new();
    let (mut in_modes, mut phase, mut at) = (false, 1usize, None::<usize>);
    for l in text.lines() {
        let t = l.trim();
        if let Some(rest) = t.strip_prefix("=> Phase No.") {
            phase = rest.split_whitespace().next().and_then(|v| v.parse().ok()).unwrap_or(phase);
        }
        if t.contains("Amplitudes of symmetry modes") {
            let k = match out.iter().position(|b| b.phase == phase) {
                Some(k) => {
                    out[k] = FullProfAmplitudes { phase, ..Default::default() };
                    k
                }
                None => {
                    out.push(FullProfAmplitudes { phase, ..Default::default() });
                    out.len() - 1
                }
            };
            at = Some(k);
            in_modes = true;
            continue;
        }
        if in_modes {
            let f: Vec<&str> = t.split_whitespace().collect();
            if f.len() == 3 && f[0] != "Name" {
                if let (Ok(v), Ok(e)) = (f[1].parse::<f64>(), f[2].parse::<f64>()) {
                    out[at.unwrap()].modes.push((f[0].to_string(), v, e));
                    continue;
                }
            }
            if f.first() == Some(&"Name") {
                continue;
            }
            in_modes = false;
        }
        if let Some(rest) = t.strip_prefix("=> Global Amplitude of Representation") {
            if let (Some((name, v)), Some(last)) = (rest.split_once(':'), at.and_then(|k| out.get_mut(k))) {
                if let Some((value, sigma)) = with_sigma(v.trim()) {
                    last.irreps.push((name.trim().to_string(), value, sigma));
                }
            }
        }
    }
    out
}

/// The parent of a mode phase: its space group (a symbol or number `SpaceGroup::from_symbol` reads) and the
/// transformation (P, p) from the parent's basis to the phase's, in ITA's convention (the columns of P are the phase's
/// basis vectors in the parent's: AMPLIMODES's "Transformation matrix").
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ModeParent {
    pub group: String,
    pub basis: Mat3,
    pub origin: [f64; 3],
}

/// A reflection of a mode phase in its parent's terms.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ParentReflection {
    /// The indices in the parent's cell, h·P⁻¹ (fractional for a superstructure reflection).
    pub hkl: [f64; 3],
    /// The k-vector it belongs to, in the parent's reciprocal basis: the parent indices less the nearest vector of the
    /// parent's reciprocal lattice (integer and allowed by its centring); 0 (Γ) for a reflection of the parent's lattice.
    pub k: [f64; 3],
}

impl ModeParent {
    /// The parent indices and k-vector of the phase's reflection h (h·P⁻¹, h a row vector, since reciprocal indices
    /// transform as h_phase = h_parent·P). Of the lattice vectors G of the parent's reciprocal lattice (integers with
    /// G·t integral for its centring translations t), the one nearest h·P⁻¹ in index length; ties go to the first in the
    /// order searched, so k is one member of its star.
    pub fn reflection(&self, hkl: [i32; 3]) -> Result<ParentReflection, String> {
        let (g, _) = SpaceGroup::from_symbol(&self.group)?;
        let inv = mat3::inverse(&self.basis).ok_or("the transformation is singular")?;
        let hp: [f64; 3] = [0, 1, 2].map(|j| (0..3).map(|i| hkl[i] as f64 * inv[i][j]).sum());
        let identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        let centring: Vec<[i32; 3]> = g.ops.iter().filter(|op| op.r == identity && op.t != [0, 0, 0]).map(|op| op.t).collect();
        let base = hp.map(|v| v.round() as i32);
        let mut best: Option<([i32; 3], f64)> = None;
        for dh in -1..=1 {
            for dk in -1..=1 {
                for dl in -1..=1 {
                    let gv = [base[0] + dh, base[1] + dk, base[2] + dl];
                    if centring.iter().any(|t| (0..3).map(|i| gv[i] * t[i]).sum::<i32>().rem_euclid(crate::symop::DEN) != 0) {
                        continue;
                    }
                    let d2: f64 = (0..3).map(|i| (hp[i] - gv[i] as f64).powi(2)).sum();
                    if best.is_none_or(|b| d2 < b.1 - 1e-9) {
                        best = Some((gv, d2));
                    }
                }
            }
        }
        let (gv, _) = best.ok_or("no vector of the parent's reciprocal lattice near the reflection")?;
        let clean = |v: f64| if v.abs() < 1e-9 { 0.0 } else { v };
        Ok(ParentReflection { hkl: hp.map(clean), k: [0, 1, 2].map(|i| clean(hp[i] - gv[i] as f64)) })
    }
}

/// An irrep's distortion: its modes, amplitude and polarisation vector.
#[derive(Clone, Debug, PartialEq)]
pub struct Irrep {
    pub name: String,
    /// The modes of this irrep (indices into `Modes::modes`), in order.
    pub modes: Vec<usize>,
    /// A = √(Σ Q²) over the irrep's modes (Å).
    pub amplitude: f64,
    /// Q_m / A for each of the irrep's modes (all zero when A = 0).
    pub polarisation: Vec<f64>,
}

/// The irrep of a mode named as AMPLIMODES names them ("A3_GM4+", "Q2_GM4-"): the text after the first '_'.
pub fn irrep_of_name(name: &str) -> Option<&str> {
    name.split_once('_').map(|(_, i)| i).filter(|i| !i.is_empty())
}

impl Modes {
    /// The modes of a FullProf phase (Jbt = 6): the atoms' positions as the reference, each V_MODES line on the atom
    /// with its label, the irrep of a mode with no vector taken from its name.
    pub fn from_pcr(m: &crate::pcr::Modes, atoms: &[crate::pcr::PcrAtom]) -> Result<Modes, String> {
        let mut modes: Vec<Mode> = m.amplitudes.iter().map(|(n, a)| Mode { name: n.clone(), irrep: irrep_of_name(n).unwrap_or("").to_string(), amplitude: *a, vectors: Vec::new() }).collect();
        for (k, label, irrep, v, coeff) in &m.vectors {
            let atom = atoms.iter().position(|a| &a.label == label).ok_or(format!("V_MODES: mode {k} names atom {label}, which is not in the phase"))?;
            let mode = &mut modes[k - 1];
            if mode.vectors.is_empty() {
                mode.irrep = irrep.clone();
            } else if &mode.irrep != irrep {
                return Err(format!("V_MODES: mode {k} is given as {} and as {irrep}", mode.irrep));
            }
            mode.vectors.push(ModeVector { atom, v: *v, coeff: *coeff });
        }
        Ok(Modes { reference: atoms.iter().map(|a| a.x).collect(), modes, reference_cell: None, max_amplitude: m.max_amplitude, fst: m.fst.clone(), parent: None, sign_reference: None })
    }

    /// The displacement Σ_m Q_m c_m V_mj of atom j (fractions of the cell).
    pub fn displacement(&self, j: usize) -> [f64; 3] {
        let mut u = [0.0; 3];
        for m in &self.modes {
            for v in m.vectors.iter().filter(|v| v.atom == j) {
                for (uk, vk) in u.iter_mut().zip(v.v) {
                    *uk += m.amplitude * v.coeff * vk;
                }
            }
        }
        u
    }

    /// The displacement of every atom for the given amplitudes, one per mode (E6: a morph r⁰ + t·u at a chosen t, or
    /// another solution's structure): Σ_m Q_m c_m V_mj, fractions of the cell.
    pub fn displacement_of(&self, amplitudes: &[f64]) -> Vec<[f64; 3]> {
        let mut u = vec![[0.0; 3]; self.reference.len()];
        for (m, q) in self.modes.iter().zip(amplitudes) {
            for v in &m.vectors {
                for (uk, vk) in u[v.atom].iter_mut().zip(v.v) {
                    *uk += q * v.coeff * vk;
                }
            }
        }
        u
    }

    /// The displacement of every atom by one irrep's modes alone, at their amplitudes (E6: the arrows of an irrep).
    /// The irreps' displacements sum to the whole.
    pub fn irrep_displacement(&self, irrep: &str) -> Vec<[f64; 3]> {
        self.displacement_of(&self.modes.iter().map(|m| if m.irrep == irrep { m.amplitude } else { 0.0 }).collect::<Vec<_>>())
    }

    /// The position of every atom: reference plus displacement.
    pub fn positions(&self) -> Vec<[f64; 3]> {
        (0..self.reference.len()).map(|j| {
            let u = self.displacement(j);
            [0, 1, 2].map(|k| self.reference[j][k] + u[k])
        }).collect()
    }

    /// ∂x_jk/∂Q_m = c_m V_mjk: the (atom, coordinate, coefficient) a mode moves.
    pub fn derivative(&self, m: usize) -> Vec<(usize, usize, f64)> {
        self.modes[m].vectors.iter().flat_map(|v| (0..3).filter(move |&k| v.v[k] != 0.0).map(move |k| (v.atom, k, v.coeff * v.v[k]))).collect()
    }

    /// The irreps, in the order their first mode appears.
    pub fn irreps(&self) -> Vec<Irrep> {
        let mut out: Vec<Irrep> = Vec::new();
        for (i, m) in self.modes.iter().enumerate() {
            match out.iter_mut().find(|r| r.name == m.irrep) {
                Some(r) => r.modes.push(i),
                None => out.push(Irrep { name: m.irrep.clone(), modes: vec![i], amplitude: 0.0, polarisation: Vec::new() }),
            }
        }
        for r in &mut out {
            r.amplitude = r.modes.iter().map(|&i| self.modes[i].amplitude.powi(2)).sum::<f64>().sqrt();
            r.polarisation = r.modes.iter().map(|&i| if r.amplitude > 0.0 { self.modes[i].amplitude / r.amplitude } else { 0.0 }).collect();
        }
        out
    }

    /// The standard deviation of an irrep's amplitude from the covariance of the modes' amplitudes, cov(m, n):
    /// σ_A² = Σ (Q_m/A)(Q_n/A) cov(m, n). With cov diagonal (σ_m²), FullProf's rule.
    pub fn irrep_esd(&self, r: &Irrep, cov: impl Fn(usize, usize) -> f64) -> f64 {
        let mut s = 0.0;
        for (a, &m) in r.polarisation.iter().zip(&r.modes) {
            for (b, &n) in r.polarisation.iter().zip(&r.modes) {
                s += a * b * cov(m, n);
            }
        }
        s.max(0.0).sqrt()
    }

    /// The standard deviations of atom j's coordinates from the covariance of the amplitudes:
    /// σ²(x_k) = Σ_mn c_m V_mjk c_n V_njk cov(m, n). With cov diagonal, FullProf's rule.
    pub fn position_esd(&self, j: usize, cov: impl Fn(usize, usize) -> f64) -> [f64; 3] {
        let d: Vec<[f64; 3]> = self.modes.iter().map(|m| {
            let mut g = [0.0; 3];
            for v in m.vectors.iter().filter(|v| v.atom == j) {
                for (gk, vk) in g.iter_mut().zip(v.v) {
                    *gk += v.coeff * vk;
                }
            }
            g
        }).collect();
        [0, 1, 2].map(|k| {
            let mut s = 0.0;
            for (m, gm) in d.iter().enumerate() {
                for (n, gn) in d.iter().enumerate() {
                    if gm[k] != 0.0 && gn[k] != 0.0 {
                        s += gm[k] * gn[k] * cov(m, n);
                    }
                }
            }
            s.max(0.0).sqrt()
        })
    }

    /// The number of images of each reference atom in the primitive cell of `group`: its multiplicity in the
    /// conventional cell over the number of centring translations.
    pub fn primitive_multiplicities(&self, group: &SpaceGroup) -> Result<Vec<f64>, String> {
        let n = group.centring.len().max(1) as f64;
        self.reference.iter().map(|&x| Site::new(group, x).map(|s| s.multiplicity as f64 / n)).collect()
    }

    /// The Gram matrix of the modes over the primitive cell of the distorted phase, in `metric`:
    /// G_mn = Σ_j mult_j (c_m V_mj)ᵀ g (c_n V_nj), g the metric tensor. AMPLIMODES's modes are orthonormal in the
    /// reference metric: G = 1 (the identity of the test `symmetry_modes`). The sum over the images of atom j is
    /// mult_j times atom j's term because the distorted group's rotations keep the reference metric.
    pub fn gram(&self, group: &SpaceGroup, metric: &Cell) -> Result<Vec<Vec<f64>>, String> {
        let mult = self.primitive_multiplicities(group)?;
        let g = metric.metric();
        let n = self.modes.len();
        let mut out = vec![vec![0.0; n]; n];
        for (a, ma) in self.modes.iter().enumerate() {
            for (b, mb) in self.modes.iter().enumerate() {
                let mut s = 0.0;
                for va in &ma.vectors {
                    for vb in mb.vectors.iter().filter(|v| v.atom == va.atom) {
                        let x = va.v.map(|c| c * va.coeff);
                        let y = vb.v.map(|c| c * vb.coeff);
                        let gy = mat3::mul_vec(&g, &y);
                        s += mult[va.atom] * (0..3).map(|k| x[k] * gy[k]).sum::<f64>();
                    }
                }
                out[a][b] = s;
            }
        }
        Ok(out)
    }

    /// The amplitudes of the displacements `u` (one per atom, fractions of the cell) in these modes: the least-squares
    /// solution of Σ_n G_mn Q_n = Σ_j mult_j (c_m V_mj)ᵀ g u_j, with the Gram matrix of `gram`. For orthonormal modes,
    /// Q_m = Σ_j mult_j (c_m V_mj)ᵀ g u_j: the projection AMPLIMODES makes of a distorted structure. Modes that move
    /// nothing are given amplitude 0. Also returns the residual: the root-mean-square size (Å) over the primitive
    /// cell of the part of u the modes do not describe (0 when u is a sum of these modes).
    pub fn project(&self, group: &SpaceGroup, metric: &Cell, u: &[[f64; 3]]) -> Result<(Vec<f64>, f64), String> {
        if u.len() != self.reference.len() {
            return Err(format!("{} displacements for {} atoms", u.len(), self.reference.len()));
        }
        let mult = self.primitive_multiplicities(group)?;
        let g = metric.metric();
        let gram = self.gram(group, metric)?;
        let n = self.modes.len();
        let mut b = vec![0.0; n];
        for (m, mode) in self.modes.iter().enumerate() {
            for v in &mode.vectors {
                let x = v.v.map(|c| c * v.coeff);
                let gu = mat3::mul_vec(&g, &u[v.atom]);
                b[m] += mult[v.atom] * (0..3).map(|k| x[k] * gu[k]).sum::<f64>();
            }
        }
        let live: Vec<usize> = (0..n).filter(|&m| gram[m][m] > 1e-12).collect();
        let a: Vec<Vec<f64>> = live.iter().map(|&i| live.iter().map(|&j| gram[i][j]).collect()).collect();
        let rhs: Vec<f64> = live.iter().map(|&i| b[i]).collect();
        let sol = solve(a, rhs).ok_or("the modes are linearly dependent")?;
        let mut q = vec![0.0; n];
        for (k, &i) in live.iter().enumerate() {
            q[i] = sol[k];
        }
        // The residual: |u|² − Qᵀ G Q over the primitive cell.
        let total: f64 = u.iter().enumerate().map(|(j, x)| mult[j] * mat3::quad(&g, x)).sum();
        let described: f64 = (0..n).map(|i| (0..n).map(|j| q[i] * gram[i][j] * q[j]).sum::<f64>()).sum();
        Ok((q, (total - described).max(0.0).sqrt()))
    }
}

/// An amplitude set that describes the same structure as another (E7): the parent's operation that relates them
/// (x,y,z in the parent's basis), the lattice translation of the parent added to it (in the parent's basis: an origin
/// shift of the phase when it is not a translation of the phase's lattice), and the amplitudes.
#[derive(Clone, Debug, PartialEq)]
pub struct Equivalent {
    pub operation: String,
    pub translation: [i32; 3],
    pub amplitudes: Vec<f64>,
}

impl Modes {
    /// Every amplitude set that describes the structure of `amplitudes` up to the parent's symmetry (E7): for each
    /// operation (W, w) of the parent's group and each translation T of its lattice, the structure moved by
    /// x ↦ P⁻¹(W(Px + p) + w + T − p) in the phase's basis ((P, p) the parent's transformation) is compared with the
    /// reference (each atom's nearest image of the same element within 1 Å, `displacements`) and projected on the
    /// modes (`project`); when the modes describe it whole (a residual below 10⁻³ Å) and the whole cell rebuilt from
    /// the amplitudes found is the whole moved cell (every atom of the cell, all the phase's operations applied, within
    /// 10⁻³ Å of one of the same element: P21, L218), its amplitudes are an equivalent set. A parent operation the
    /// phase's group does not hold changes signs (a domain of the distortion: the same structure, another setting of
    /// it); the identity gives the amplitudes back. Distinct sets only, the identity's first. `species` names each
    /// atom's element; `metric` is the reference cell (`reference_cell`, or the phase's).
    ///
    /// The whole cell is needed: `displacements` matches each atom of the asymmetric unit with the nearest image of a
    /// moved atom under the phase's own operations, which imposes the phase's symmetry on the moved structure. A parent
    /// operation outside the normaliser of the phase's group moves the structure to a domain that is not symmetric
    /// under that group in this setting; its asymmetric unit, symmetrised again, is another structure that the modes
    /// describe with no residual. Those sets move amplitude between irreps (case A: 20 of 24 before the check), which
    /// no symmetry operation does: it maps each irrep's subspace onto itself.
    pub fn equivalent_amplitudes(&self, species: &[String], group: &SpaceGroup, metric: &Cell, amplitudes: &[f64]) -> Result<Vec<Equivalent>, String> {
        let parent = self.parent.as_ref().ok_or("the phase's parent is not known (an AMPLIMODES analysis or the control file's mode_parent gives it)")?;
        let (pg, _) = SpaceGroup::from_symbol(&parent.group)?;
        let p = parent.basis;
        let inv = mat3::inverse(&p).ok_or("the transformation is singular")?;
        let u = self.displacement_of(amplitudes);
        let x: Vec<[f64; 3]> = (0..self.reference.len()).map(|j| [0, 1, 2].map(|k| self.reference[j][k] + u[j][k])).collect();
        let reference: Vec<(String, [f64; 3])> = species.iter().cloned().zip(self.reference.iter().copied()).collect();
        // The whole cell of a structure: every atom under every operation of the phase's group (centrings included).
        let whole = |pos: &[[f64; 3]]| -> Vec<(String, [f64; 3])> {
            let mut out = Vec::new();
            for (sp, xc) in species.iter().zip(pos) {
                for h in &group.ops {
                    out.push((sp.clone(), crate::site::SiteOp { r: h.r, t: h.t }.apply(*xc)));
                }
            }
            out
        };
        let cell = whole(&x);
        let element = |s: &str| s.chars().take_while(|c| c.is_ascii_alphabetic()).collect::<String>().to_ascii_lowercase();
        let gm = metric.metric();
        let same_cell = |a: &[(String, [f64; 3])], b: &[(String, [f64; 3])]| {
            a.iter().all(|(sa, xa)| {
                b.iter().any(|(sb, xb)| {
                    element(sa) == element(sb) && {
                        let d = [0, 1, 2].map(|k| {
                            let d = xa[k] - xb[k];
                            d - d.round()
                        });
                        mat3::quad(&gm, &d).sqrt() <= 1e-3
                    }
                })
            })
        };
        // The parent's lattice translations that are distinct in the phase's cell (P⁻¹T modulo 1).
        let n = p.iter().flatten().map(|v| v.abs().ceil() as i32).max().unwrap_or(1).max(1);
        // The shortest T for each (the first of each length in the order searched).
        let mut all: Vec<[i32; 3]> = Vec::new();
        for a in -n..=n {
            for b in -n..=n {
                for c in -n..=n {
                    all.push([a, b, c]);
                }
            }
        }
        all.sort_by_key(|t| t.iter().map(|v| v.abs()).sum::<i32>());
        let mut shifts: Vec<([i32; 3], [f64; 3])> = Vec::new();
        for t in all {
            let s = mat3::mul_vec(&inv, &t.map(f64::from)).map(|v| v - v.floor());
            let s = s.map(|v| if (v - 1.0).abs() < 1e-9 { 0.0 } else { v });
            if !shifts.iter().any(|(_, o)| (0..3).all(|k| (o[k] - s[k]).abs() < 1e-9)) {
                shifts.push((t, s));
            }
        }
        let mut out: Vec<Equivalent> = Vec::new();
        for op in &pg.ops {
            let w = op.t.map(|v| v as f64 / crate::symop::DEN as f64);
            for (t, _) in &shifts {
                let moved: Vec<(String, [f64; 3])> = species.iter().zip(&x).map(|(sp, xc)| {
                    let xp = mat3::mul_vec(&p, xc);
                    let xp = [0, 1, 2].map(|k| xp[k] + parent.origin[k]);
                    let y: [f64; 3] = std::array::from_fn(|i| (0..3).map(|j| op.r[i][j] as f64 * xp[j]).sum::<f64>() + w[i] + t[i] as f64 - parent.origin[i]);
                    (sp.clone(), mat3::mul_vec(&inv, &y))
                }).collect();
                let Ok(d) = displacements(&reference, &moved, group, metric, [0.0; 3], 1.0) else { continue };
                let Ok((q, residual)) = self.project(group, metric, &d) else { continue };
                if residual > 1e-3 || out.iter().any(|e| e.amplitudes.iter().zip(&q).all(|(a, b)| (a - b).abs() < 1e-6)) {
                    continue;
                }
                // The whole moved cell against the whole cell of the amplitudes found (P21).
                let moved_cell: Vec<(String, [f64; 3])> = cell.iter().map(|(sp, xc)| {
                    let xp = mat3::mul_vec(&p, xc);
                    let xp = [0, 1, 2].map(|k| xp[k] + parent.origin[k]);
                    let y: [f64; 3] = std::array::from_fn(|i| (0..3).map(|j| op.r[i][j] as f64 * xp[j]).sum::<f64>() + w[i] + t[i] as f64 - parent.origin[i]);
                    (sp.clone(), mat3::mul_vec(&inv, &y))
                }).collect();
                let uq = self.displacement_of(&q);
                let xq: Vec<[f64; 3]> = (0..self.reference.len()).map(|j| [0, 1, 2].map(|k| self.reference[j][k] + uq[j][k])).collect();
                if !same_cell(&moved_cell, &whole(&xq)) {
                    continue;
                }
                out.push(Equivalent { operation: op.to_xyz(), translation: *t, amplitudes: q });
            }
        }
        Ok(out)
    }

    /// Whether two amplitude sets describe the same structure up to the parent's symmetry (E7), within `tol` (Å) on
    /// every amplitude: the operation that relates them, or None.
    pub fn same_structure(&self, species: &[String], group: &SpaceGroup, metric: &Cell, a: &[f64], b: &[f64], tol: f64) -> Result<Option<Equivalent>, String> {
        Ok(self.equivalent_amplitudes(species, group, metric, a)?.into_iter().find(|e| e.amplitudes.iter().zip(b).all(|(x, y)| (x - y).abs() <= tol)))
    }

    /// Of the sets equivalent to `amplitudes`, the one nearest `reference` (E7: AMPLIMODES's signs as the convention,
    /// its analysis's amplitudes by mode name), in the sum of squared differences.
    pub fn nearest_to(&self, species: &[String], group: &SpaceGroup, metric: &Cell, amplitudes: &[f64], reference: &[f64]) -> Result<Equivalent, String> {
        let d = |e: &Equivalent| e.amplitudes.iter().zip(reference).map(|(a, b)| (a - b).powi(2)).sum::<f64>();
        self.equivalent_amplitudes(species, group, metric, amplitudes)?.into_iter().min_by(|a, b| d(a).total_cmp(&d(b))).ok_or("no equivalent set (not even the identity's)".into())
    }

    /// The sign reference's amplitudes in the order of the modes, by name; an error naming a mode it does not give.
    pub fn reference_amplitudes(&self) -> Option<Result<Vec<f64>, String>> {
        let r = self.sign_reference.as_ref()?;
        Some(self.modes.iter().map(|m| r.amplitudes.iter().find(|a| a.0 == m.name).map(|a| a.1).ok_or(format!("{} gives no amplitude for {}", r.file, m.name))).collect())
    }

    /// For a state's amplitudes, the equivalent set nearest the sign reference and the parent's operation that brings
    /// the state to it (E14b): the identity when the state is already in the reference's setting. None without a
    /// sign reference.
    pub fn reference_signs(&self, species: &[String], group: &SpaceGroup, metric: &Cell, amplitudes: &[f64]) -> Option<Result<Equivalent, String>> {
        let reference = match self.reference_amplitudes()? {
            Ok(r) => r,
            Err(e) => return Some(Err(e)),
        };
        Some(self.nearest_to(species, group, metric, amplitudes, &reference))
    }
}

impl Equivalent {
    /// Whether this is the state itself: the identity and no translation.
    pub fn is_identity(&self) -> bool {
        self.operation == "x,y,z" && self.translation == [0, 0, 0]
    }

    /// In words, for the view (E14b): "the reference's setting", or "brought to the reference's setting by 2 about
    /// [0 0 1] (-x,-y,z)", the rotation's kind and axis (or a mirror's normal) in the parent's basis, and the
    /// translation when there is one.
    pub fn words(&self) -> String {
        if self.is_identity() {
            return "the reference's setting".into();
        }
        let t = if self.translation == [0, 0, 0] { String::new() } else { format!(" + ({}, {}, {})", self.translation[0], self.translation[1], self.translation[2]) };
        let kind = crate::symop::SymOp::parse_xyz(&self.operation).ok().and_then(|op| {
            let rt = op.rot_type()?;
            // The axis: a non-zero column of R + I for a proper rotation (of R − I for an improper one, its normal).
            let sign = if op.det() > 0 { 1 } else { -1 };
            let m: [[i32; 3]; 3] = std::array::from_fn(|i| std::array::from_fn(|j| op.r[i][j] * sign + if i == j { 1 } else { 0 }));
            let axis = (0..3).map(|j| [m[0][j], m[1][j], m[2][j]]).find(|c| c.iter().any(|&v| v != 0)).map(|c| {
                let g = c.iter().fold(0i32, |g, &v| gcd(g, v.abs())).max(1);
                c.map(|v| v / g)
            });
            Some(match (rt, axis) {
                (crate::symop::RotType::One, _) => "the translation".to_string(),
                (crate::symop::RotType::Inversion, _) => "the inversion".to_string(),
                (crate::symop::RotType::Mirror, Some(a)) => format!("the mirror normal to [{} {} {}]", a[0], a[1], a[2]),
                (r, Some(a)) => format!("{} about [{} {} {}]", r.symbol(), a[0], a[1], a[2]),
                (r, None) => r.symbol().to_string(),
            })
        });
        format!("brought to the reference's setting by {} ({}){t}", kind.unwrap_or_else(|| "the operation".into()), self.operation)
    }
}

fn gcd(a: i32, b: i32) -> i32 {
    if b == 0 { a } else { gcd(b, a % b) }
}

/// Gaussian elimination with partial pivoting; None when singular.
fn solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for c in 0..n {
        let p = (c..n).max_by(|&i, &j| a[i][c].abs().total_cmp(&a[j][c].abs()))?;
        if a[p][c].abs() < 1e-12 {
            return None;
        }
        a.swap(c, p);
        b.swap(c, p);
        let pivot = a[c].clone();
        for r in c + 1..n {
            let f = a[r][c] / pivot[c];
            for (x, y) in a[r][c..].iter_mut().zip(&pivot[c..]) {
                *x -= f * y;
            }
            b[r] -= f * b[c];
        }
    }
    let mut x = vec![0.0; n];
    for r in (0..n).rev() {
        x[r] = (b[r] - (r + 1..n).map(|k| a[r][k] * x[k]).sum::<f64>()) / a[r][r];
    }
    Some(x)
}

/// The cell of basis P (columns: the new basis vectors in the old basis, ITA's convention, as AMPLIMODES prints its
/// "Transformation matrix") over a cell: metric Pᵀ g P.
/// A morph between a reference structure and the distorted one (Explore's structure view, R29; L148): each image's
/// Cartesian position r = r⁰ + t Σ u over the irreps chosen (all when none is named), the sum of the chosen
/// displacements per image (the arrows) and its length in Å.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Morph {
    pub positions: Vec<[f64; 3]>,
    pub total: Vec<[f64; 3]>,
    pub lengths: Vec<f64>,
}

pub fn morph(reference: &[[f64; 3]], displacements: &[(String, Vec<[f64; 3]>)], t: f64, chosen: Option<&[String]>) -> Morph {
    let used: Vec<&Vec<[f64; 3]>> = displacements.iter().filter(|(n, _)| chosen.is_none_or(|c| c.contains(n))).map(|(_, u)| u).collect();
    let total: Vec<[f64; 3]> = (0..reference.len()).map(|i| std::array::from_fn(|k| used.iter().map(|u| u.get(i).map_or(0.0, |v| v[k])).sum())).collect();
    let positions = reference.iter().zip(&total).map(|(r, u)| std::array::from_fn(|k| r[k] + t * u[k])).collect();
    let lengths = total.iter().map(|u| (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt()).collect();
    Morph { positions, total, lengths }
}

pub fn transform_cell(cell: &Cell, p: &Mat3) -> Result<Cell, String> {
    let g = cell.metric();
    let pt = [0, 1, 2].map(|i| [0, 1, 2].map(|j| p[j][i]));
    let h = mat3::mul(&pt, &mat3::mul(&g, p));
    let (a, b, c) = (h[0][0].sqrt(), h[1][1].sqrt(), h[2][2].sqrt());
    let ang = |x: f64, l: f64, m: f64| (x / (l * m)).clamp(-1.0, 1.0).acos().to_degrees();
    Cell::new(a, b, c, ang(h[1][2], b, c), ang(h[0][2], a, c), ang(h[0][1], a, b)).map_err(|e| e.to_string())
}

/// The displacements of a distorted structure from the reference: for each reference atom, the image (by the
/// group's operators and lattice translations) of a distorted atom of the same element nearest to it, after the
/// distorted structure is moved by `shift` (fractions; an origin change such as AMPLIMODES's "Sb to the origin").
/// Elements are compared without charge or case ("SR+2" is "Sr"). Refused when an atom has no partner within
/// `max_distance` Å in `metric`.
pub fn displacements(reference: &[(String, [f64; 3])], distorted: &[(String, [f64; 3])], group: &SpaceGroup, metric: &Cell, shift: [f64; 3], max_distance: f64) -> Result<Vec<[f64; 3]>, String> {
    let element = |s: &str| s.chars().take_while(|c| c.is_ascii_alphabetic()).collect::<String>().to_ascii_lowercase();
    let g = metric.metric();
    let mut out = Vec::new();
    for (species, r) in reference {
        let mut best: Option<(f64, [f64; 3])> = None;
        for (s, x) in distorted.iter().filter(|(s, _)| element(s) == element(species)) {
            let x = [0, 1, 2].map(|k| x[k] + shift[k]);
            for op in &group.ops {
                let y = crate::site::SiteOp { r: op.r, t: op.t }.apply(x);
                let d = [0, 1, 2].map(|k| {
                    let d = y[k] - r[k];
                    d - d.round()
                });
                let dist = mat3::quad(&g, &d).sqrt();
                if best.is_none_or(|(b, _)| dist < b) {
                    best = Some((dist, d));
                }
            }
            let _ = s;
        }
        match best {
            Some((dist, d)) if dist <= max_distance => out.push(d),
            Some((dist, _)) => return Err(format!("no {species} of the distorted structure within {max_distance} Å of the reference's at {r:?} (nearest {dist:.3} Å)")),
            None => return Err(format!("no {species} in the distorted structure")),
        }
    }
    Ok(out)
}

// AMPLIMODES's text output --------------------------------------------------------------------------------------------

/// A structure as AMPLIMODES prints it: space-group number, cell, atoms (label, Wyckoff letter, position).
#[derive(Clone, Debug, PartialEq)]
pub struct AmpliStructure {
    pub group: u32,
    pub cell: [f64; 6],
    pub atoms: Vec<(String, String, [f64; 3])>,
}

/// One line of AMPLIMODES's "Summary of Amplitudes".
#[derive(Clone, Debug, PartialEq)]
pub struct AmpliIrrep {
    pub k: String,
    pub irrep: String,
    pub direction: String,
    pub isotropy: String,
    pub dimension: usize,
    pub amplitude: f64,
}

/// AMPLIMODES's block for FullProf: atoms (label, species, position), V_MODES lines (mode, atom, irrep, vector,
/// coefficient), amplitudes (name, value).
pub type FullProfBlock = (Vec<(String, String, [f64; 3])>, Vec<(usize, String, String, [f64; 3], f64)>, Vec<(String, f64)>);

/// One analysis of AMPLIMODES (a block "Symmetry mode analysis").
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Analysis {
    pub high: Option<AmpliStructure>,
    pub low: Option<AmpliStructure>,
    /// The setting named for the low-symmetry structure ("P 1 21/n 1"), when given.
    pub low_setting: Option<String>,
    /// The transformation (P, p) from the parent's basis to the distorted phase's.
    pub transformation: Option<(Mat3, [f64; 3])>,
    pub reference: Option<AmpliStructure>,
    /// Atom mappings: reference label, its position, the distorted atom's label and position (the distorted
    /// phase's basis).
    pub mappings: Vec<(String, [f64; 3], String, [f64; 3])>,
    pub summary: Vec<AmpliIrrep>,
    pub global: Option<f64>,
    /// The block for FullProf.
    pub fullprof: Option<FullProfBlock>,
}

fn fraction(s: &str) -> Option<f64> {
    let s = s.trim();
    match s.split_once('/') {
        Some((a, b)) => Some(a.trim().parse::<f64>().ok()? / b.trim().parse::<f64>().ok()?),
        None => s.parse().ok(),
    }
}

fn floats(s: &str) -> Vec<f64> {
    s.split(|c: char| c.is_whitespace() || c == ',' || c == '(' || c == ')').filter(|t| !t.is_empty()).map_while(fraction).collect()
}

/// A structure starting at line i (group number, cell, count, atoms); the lines read.
fn structure(ls: &[&str], mut i: usize) -> Option<(AmpliStructure, usize)> {
    while i < ls.len() && (ls[i].trim().is_empty() || ls[i].starts_with("Setting used") || ls[i].starts_with("Transformation Matrix:")) {
        i += 1;
    }
    let group = ls.get(i)?.trim().parse::<u32>().ok()?;
    let c = floats(ls.get(i + 1)?);
    let n = ls.get(i + 2)?.trim().parse::<usize>().ok()?;
    if c.len() < 6 {
        return None;
    }
    let mut atoms = Vec::new();
    for l in ls.get(i + 3..i + 3 + n)? {
        let f: Vec<&str> = l.split('\t').map(str::trim).filter(|t| !t.is_empty()).collect();
        let f = if f.len() < 5 { l.split_whitespace().collect() } else { f };
        // element, number, [wyckoff], x, y, z
        let xyz: Vec<f64> = f.iter().rev().take(3).rev().filter_map(|t| fraction(t)).collect();
        if f.len() < 5 || xyz.len() != 3 {
            return None;
        }
        let wyckoff = if f.len() >= 6 { f[2].to_string() } else { String::new() };
        atoms.push((format!("{}{}", f[0], f[1]), wyckoff, [xyz[0], xyz[1], xyz[2]]));
    }
    Some((AmpliStructure { group, cell: [c[0], c[1], c[2], c[3], c[4], c[5]], atoms }, i + 3 + n))
}

/// The analyses of an AMPLIMODES text output, in order.
pub fn read_amplimodes(text: &str) -> Vec<Analysis> {
    let ls: Vec<&str> = text.lines().map(|l| l.trim_end()).collect();
    let starts: Vec<usize> = ls.iter().enumerate().filter(|(_, l)| l.trim() == "Symmetry mode analysis").map(|(i, _)| i).collect();
    let mut out = Vec::new();
    for (k, &s) in starts.iter().enumerate() {
        let end = starts.get(k + 1).copied().unwrap_or(ls.len());
        let mut a = Analysis::default();
        let mut i = s + 1;
        while i < end {
            let t = ls[i].trim();
            if t == "High symmetry structure" {
                if let Some((st, j)) = structure(&ls[..end], i + 1) {
                    a.high = Some(st);
                    i = j;
                    continue;
                }
            } else if t == "Low symmetry structure" {
                a.low_setting = ls.get(i + 1).and_then(|l| l.strip_prefix("Setting used:")).map(|s| s.trim().to_string());
                if let Some((st, j)) = structure(&ls[..end], i + 1) {
                    a.low = Some(st);
                    i = j;
                    continue;
                }
            } else if t == "Transformation matrix" && i + 3 < end {
                let rows: Vec<Vec<f64>> = (1..4).map(|r| floats(&ls[i + r].replace(['[', ']'], " "))).collect();
                if rows.iter().all(|r| r.len() == 4) {
                    let p = [0, 1, 2].map(|r| [rows[r][0], rows[r][1], rows[r][2]]);
                    a.transformation = Some((p, [rows[0][3], rows[1][3], rows[2][3]]));
                    i += 4;
                    continue;
                }
            } else if t == "Reference Structure" {
                if let Some((st, j)) = structure(&ls[..end], i + 1) {
                    a.reference = Some(st);
                    i = j;
                    continue;
                }
            } else if t == "Atom Mappings" {
                i += 2;
                // WP, site, reference atom, (x,y,z), distorted atom, (x,y,z): separated by tabs or (in 0_3) spaces.
                while i < end && !ls[i].trim().is_empty() {
                    let f: Vec<&str> = ls[i].split_whitespace().collect();
                    if f.len() >= 6 {
                        let (r, d) = (floats(f[3]), floats(f[5]));
                        if r.len() == 3 && d.len() == 3 {
                            a.mappings.push((f[2].to_string(), [r[0], r[1], r[2]], f[4].to_string(), [d[0], d[1], d[2]]));
                        }
                    }
                    i += 1;
                }
                continue;
            } else if t.starts_with('(') && a.summary.len() < 64 {
                // k-vector, irrep, direction, isotropy subgroup ("Fm-3m (225)"), dimension, amplitude.
                let f: Vec<&str> = t.split_whitespace().collect();
                let n = f.len();
                if n >= 6 {
                    if let (Ok(dimension), Ok(amplitude)) = (f[n - 2].parse(), f[n - 1].parse()) {
                        a.summary.push(AmpliIrrep { k: f[0].into(), irrep: f[1].into(), direction: f[2].into(), isotropy: f[3..n - 2].join(" "), dimension, amplitude });
                    }
                }
            } else if let Some(g) = t.strip_prefix("Global distortion:") {
                a.global = floats(g.trim_end_matches('Å')).first().copied();
            } else if t.starts_with("You can copy and paste") && a.fullprof.is_none() {
                let mut atoms = Vec::new();
                let mut vectors = Vec::new();
                let mut amps = Vec::new();
                let mut j = i + 1;
                // Atom lines (label, species, x, y, z, B, Occ, …) each followed by a codeword line.
                while j < end && !ls[j].trim_start().starts_with("V_MODES") {
                    let f: Vec<&str> = ls[j].split_whitespace().collect();
                    if f.len() >= 7 && f[0].parse::<f64>().is_err() && !ls[j].trim_start().starts_with('!') {
                        let x: Vec<f64> = f[2..5].iter().filter_map(|t| t.parse().ok()).collect();
                        if x.len() == 3 {
                            atoms.push((f[0].to_string(), f[1].to_string(), [x[0], x[1], x[2]]));
                        }
                    }
                    j += 1;
                }
                let nv = ls.get(j).map(|l| floats(&l.trim_start()[7..]).first().copied().unwrap_or(0.0) as usize).unwrap_or(0);
                j += 1;
                while j < end && vectors.len() < nv {
                    let f: Vec<&str> = ls[j].split_whitespace().collect();
                    if f.len() >= 7 && !ls[j].trim_start().starts_with('!') {
                        let v: Vec<f64> = f[3..7].iter().filter_map(|t| t.parse().ok()).collect();
                        if let (Ok(m), 4) = (f[0].parse::<usize>(), v.len()) {
                            vectors.push((m, f[1].to_string(), f[2].to_string(), [v[0], v[1], v[2]], v[3]));
                        }
                    }
                    j += 1;
                }
                while j < end && !ls[j].trim_start().starts_with("A_MODES") {
                    j += 1;
                }
                let na = ls.get(j).map(|l| floats(&l.trim_start()[7..]).first().copied().unwrap_or(0.0) as usize).unwrap_or(0);
                j += 1;
                while j < end && amps.len() < na {
                    let f: Vec<&str> = ls[j].split_whitespace().collect();
                    if f.len() >= 2 && !f[0].starts_with("Max_Amplitude") {
                        if let Ok(v) = f[1].parse::<f64>() {
                            amps.push((f[0].to_string(), v));
                        }
                    }
                    j += 1;
                }
                a.fullprof = Some((atoms, vectors, amps));
                i = j;
                continue;
            }
            i += 1;
        }
        // A heading repeated with nothing under it (0_2, 0_3: "Symmetry mode analysis" twice) is no analysis.
        if a.high.is_some() || !a.summary.is_empty() {
            out.push(a);
        }
    }
    out
}

impl Analysis {
    /// The modes of the analysis's block for FullProf: its atoms as the reference, each V_MODES line on its atom, the
    /// amplitudes of A_MODES in order (P13: they name the modes computed from the structures).
    pub fn modes(&self) -> Option<Modes> {
        let (atoms, lines, amps) = self.fullprof.as_ref()?;
        let mut modes: Vec<Mode> = Vec::new();
        for (k, label, irrep, v, c) in lines {
            let atom = atoms.iter().position(|a| a.0 == *label)?;
            while modes.len() < *k {
                let i = modes.len();
                modes.push(Mode { name: amps.get(i).map_or(format!("A{}", i + 1), |a| a.0.clone()), irrep: String::new(), amplitude: amps.get(i).map_or(0.0, |a| a.1), vectors: Vec::new() });
            }
            let m = &mut modes[k - 1];
            m.irrep = irrep.clone();
            m.vectors.push(ModeVector { atom, v: *v, coeff: *c });
        }
        Some(Modes { reference: atoms.iter().map(|a| a.2).collect(), modes, reference_cell: self.reference_cell(), max_amplitude: None, fst: Vec::new(), parent: self.mode_parent(), sign_reference: None })
    }

    /// The parent of the analysis: the high-symmetry structure's space group (its number) and the transformation.
    pub fn mode_parent(&self) -> Option<ModeParent> {
        let (h, (p, o)) = (self.high.as_ref()?, self.transformation.as_ref()?);
        Some(ModeParent { group: h.group.to_string(), basis: *p, origin: *o })
    }

    /// The reference cell in the distorted phase's basis: the parent's cell transformed by P.
    pub fn reference_cell(&self) -> Option<Cell> {
        let (h, (p, _)) = (self.high.as_ref()?, self.transformation.as_ref()?);
        let c = h.cell;
        transform_cell(&Cell::new(c[0], c[1], c[2], c[3], c[4], c[5]).ok()?, p).ok()
    }
}

#[cfg(test)]
mod morph_tests {
    /// L148: r⁰ + tΣu over the irreps chosen; the sum of their displacements.
    #[test]
    fn a_morph_between_reference_and_distorted() {
        let r0 = [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]];
        let u = vec![("GM1+".to_string(), vec![[0.1, 0.0, 0.0], [0.0, 0.2, 0.0]]), ("X5+".to_string(), vec![[0.0, 0.0, 0.4], [0.0, 0.0, 0.0]])];
        let m = super::morph(&r0, &u, 0.5, None);
        assert_eq!(m.total, vec![[0.1, 0.0, 0.4], [0.0, 0.2, 0.0]]);
        assert_eq!(m.positions, vec![[0.05, 0.0, 0.2], [1.0, 1.1, 1.0]]);
        assert!((m.lengths[0] - 0.17f64.sqrt()).abs() < 1e-12 && (m.lengths[1] - 0.2).abs() < 1e-12);
        let m = super::morph(&r0, &u, 1.0, Some(&["X5+".to_string()]));
        assert_eq!(m.positions[0], [0.0, 0.0, 0.4]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cubic_to_monoclinic_metric() {
        // AMPLIMODES's P2₁/c standard setting of Fm-3m (course file 0_1): 8.087481 5.718712 9.905101, β 144.735611.
        let c = Cell::new(8.087481, 8.087481, 8.087481, 90.0, 90.0, 90.0).unwrap();
        let p = [[0.0, 0.5, -0.5], [1.0, 0.0, -1.0], [0.0, -0.5, -0.5]];
        let t = transform_cell(&c, &p).unwrap();
        assert!((t.a - 8.087481).abs() < 1e-6 && (t.b - 5.718712).abs() < 1e-6 && (t.c - 9.905101).abs() < 1e-6, "{t:?}");
        assert!((t.beta - 144.735611).abs() < 1e-5 && (t.alpha - 90.0).abs() < 1e-9, "{t:?}");
    }

    #[test]
    fn irreps_of_amplitudes() {
        let m = |n: &str, a: f64| Mode { name: n.into(), irrep: irrep_of_name(n).unwrap().into(), amplitude: a, vectors: vec![] };
        let ms = Modes { reference: vec![], modes: vec![m("A1_GM1+", -0.3), m("A10_X5+", 0.110324), m("A11_X5+", -0.051054)], reference_cell: None, max_amplitude: None, fst: vec![], parent: None, sign_reference: None };
        let r = ms.irreps();
        assert_eq!(r.iter().map(|r| r.name.as_str()).collect::<Vec<_>>(), ["GM1+", "X5+"]);
        // FullProf's case A: X5+ 0.122(6) from 0.110324(4948) and −0.051054(8631), correlations neglected.
        assert!((r[1].amplitude - 0.12156).abs() < 1e-5);
        let s = [0.0, 0.004948, 0.008631];
        let e = ms.irrep_esd(&r[1], |a, b| if a == b { s[a] * s[a] } else { 0.0 });
        assert!((e - 0.0058).abs() < 1e-4, "{e}");
    }
}
