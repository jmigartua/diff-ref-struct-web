//! Reading FullProf control files (`.pcr`), the parts needed to rebuild the model of a constant-wavelength powder
//! pattern: the job flags, the wavelengths and instrument line, the angular range, the background, the excluded
//! regions, the zero and line shifts, and for each phase its atoms, symmetry modes, scale, profile, cell, preferred
//! orientation and asymmetry. Codewords are kept as they are written (`codes`), not interpreted.
//!
//! The format is the single-pattern ("COMM") format of the FullProf manual, appendix "CODFIL.PCR" (pp. 65–126): line 4
//! (Job … Aut), line 6 (resolution file), line 7 (Ipr … Ana), line 8 (Lambda1 … Rpolarz), line 9 (NCY … Sent0),
//! lines 10–11 (background points, excluded regions), 13 (number of refined parameters), 14 (Zero, Sycos, Sysin,
//! Lambda, More), 15 (micro-absorption), 17 (background coefficients), and per phase lines 18–33. The flags of lines 4,
//! 7, 19 and 19-1 are read by the names of the comment line FullProf writes above them, so that files of different
//! FullProf versions (with more or fewer columns) are read alike; without such a comment the manual's order is used.
//! Phases are delimited by FullProf's comment "Data for PHASE number"; lines of a phase after line 33 (hkl shifts,
//! size and strain models, propagation vectors) are kept unread in `extra`.
//!
//! Line 12, the user's scattering factors (Nsc ≠ 0, manual pp. 70 and 79), is read as written (`PcrScatterer`).
//!
//! The multi-pattern format (an NPATT line, FullProf_Manual.pdf pp. 65–66, 88) is split into one single-pattern text per
//! pattern (`split_patterns`, P06), each read as above.
//!
//! TOF patterns (Uni = 1, P05): line 8 is Bkpos, Wdt, Iabscor; line 14 Zero, Dtt1, Dtt2 [, Dtt_1overd], 2ThetaBank; line
//! 16 (Npr = 10) Hodges' Zerot, Dtt1t, Dtt2t, x-cross, Width; and per phase lines 26 (Scale, Extinc, …), 27 and 27-2
//! (Sig-2 … and Gam-2 …), 30 (Pref1, Pref2, alph0, beta0, alph1, beta1 [, alphQ, betaQ]), 31 (Npr = 10: alph0t …) and
//! 33 (ABSCOR1, ABSCOR2) (FullProf_Manual.pdf pp. 75–76, 82–85, 108–112; `PcrTof`).
//!
//! Not read: energy-dispersive patterns (Uni = 2).

use std::collections::BTreeMap;
use std::fmt;

#[derive(Clone, Debug, PartialEq)]
pub struct PcrError {
    pub line: usize,
    pub message: String,
}

impl fmt::Display for PcrError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "line {}: {}", self.line, self.message)
    }
}

/// The background of a pattern.
#[derive(Clone, Debug, PartialEq)]
pub enum Background {
    /// Nba ≥ 2: points (2θ, intensity) with linear interpolation (Nba ≤ −4 written as points: cubic splines).
    Points(Vec<(f64, f64)>),
    /// Nba = 0: polynomial Σ B_m (T/Bkpos − 1)^m, 6 coefficients (manual eq. 3.4).
    Polynomial(Vec<f64>),
    /// Nba = −5: Chebyshev polynomials, up to 24 coefficients (FullProf news, 17 June 2011).
    Chebyshev(Vec<f64>),
    /// Nba = −4: Fourier cosine series (manual eq. 3.4′).
    Fourier(Vec<f64>),
    /// Other values of Nba: coefficients as read, not modelled.
    Other(i32, Vec<f64>),
}

/// An atom of a phase as the `.pcr` file gives it.
#[derive(Clone, Debug, PartialEq)]
pub struct PcrAtom {
    pub label: String,
    pub species: String,
    pub x: [f64; 3],
    pub biso: f64,
    pub occ: f64,
    /// N_t: 0 isotropic, 2 anisotropic (betas read), 4 form-factor parameters.
    pub n_t: i32,
    pub beta: Option<[f64; 6]>,
    /// Codewords of x, y, z, Biso, Occ (line 25-2), and of the betas (line 25-4).
    pub codes: Vec<f64>,
    pub beta_codes: Option<Vec<f64>>,
    pub line: usize,
    /// A rigid-body phase's (|Jbt| = 4, P41): the numbers of line 25-1 past Occ (P6, THETA, PHI, Spc; an atom of a
    /// Z-matrix group its connectivity N1 N2 N3) and those of line 25-3 (the internal coordinates; a group's first
    /// atom also X0 Y0 Z0, CHI, P16, DEG, KIND). Empty in the other phases.
    pub rigid: (Vec<f64>, Vec<f64>),
    /// Their codewords: those of P9–P15 (line 25-4); P6–P8's are `codes` 5–7.
    pub rigid_codes: Vec<f64>,
}

/// Symmetry modes (FullProf's V_MODES / A_MODES block, Jbt = 6; FullProf's notes of 29 August 2008 – 16 January 2009
/// in `7-FullProf-Tutorial-SymmetryModes/FullProf-Tutorial-Symmetry-Modes.pdf`, pp. 3–5): the position of each atom
/// is its reference position plus Σ_m Q_m·coeff_m·V_m, the vectors in fractions of the phase's cell.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Modes {
    /// One line of V_MODES: (mode number, atom label, irrep, polarisation vector, coefficient).
    pub vectors: Vec<(usize, String, String, [f64; 3], f64)>,
    /// Amplitudes by mode number (1-based order of the A_MODES list): name, value.
    pub amplitudes: Vec<(String, f64)>,
    /// The codeword of each amplitude.
    pub codes: Vec<f64>,
    /// The integers after the number of amplitudes on the A_MODES line: what FullProf writes to FullProf Studio files
    /// (note of 11 January 2009: 2 = the modes of each irrep together as arrows, −2 as structures, ±N all modes one by
    /// one, −3 none, or 0/1 per mode, or Nrm and a list). Not used in the calculation.
    pub fst: Vec<i32>,
    /// `Max_Amplitude`: the box |Q| ≤ value on every amplitude (note of 16 January 2009), in Å.
    pub max_amplitude: Option<f64>,
}

/// A soft distance constraint (restraint) of a phase, line 45 (manual p. 123): the distance from atom 1 to atom 2
/// transformed by the operator (ITnum, T).
#[derive(Clone, Debug, PartialEq)]
pub struct PcrDistance {
    pub atoms: [String; 2],
    pub itnum: i32,
    pub t: [f64; 3],
    pub value: f64,
    pub sigma: f64,
}

/// A soft angle constraint, line 45-1 (manual p. 124): the angle at atom 2 transformed by (ITnum1, T) between atom 1
/// and atom 3 transformed by (ITnum2, t).
#[derive(Clone, Debug, PartialEq)]
pub struct PcrAngle {
    pub atoms: [String; 3],
    pub itnum: [i32; 2],
    pub t: [[f64; 3]; 2],
    pub value: f64,
    pub sigma: f64,
}

/// A general linear restraint (NLI > 0, line 7): Σ c·p_n ≈ value ± sigma over FullProf's parameter numbers n.
#[derive(Clone, Debug, PartialEq)]
pub struct PcrLinear {
    pub name: String,
    pub value: f64,
    pub sigma: f64,
    pub terms: Vec<(f64, u32)>,
}

/// A phase's profile in a TOF pattern (Uni = 1; FullProf_Manual.pdf pp. 108–112), as written. Its extinction is
/// `PcrPhase::shape1` (line 26's second value, Extinc).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PcrTof {
    /// Line 27: Sig-2, Sig-1, Sig-0 (variances, µs²/Å⁴, µs²/Å², µs²), then Sigma-Q (the older files' Xt, unused),
    /// Iso-GStrain (Yt), Iso-GSize (Z1) and Ani-LSize (Z0).
    pub sigma: [f64; 3],
    pub sigma_q: f64,
    pub gauss_strain: f64,
    pub gauss_size: f64,
    pub lorentz_size_ani: f64,
    /// Line 27-2: Gam-2, Gam-1, Gam-0 (µs/Å², µs/Å, µs), Iso-LorStrain, Iso-LorSize.
    pub gamma: [f64; 3],
    pub lorentz_strain: f64,
    pub lorentz_size: f64,
    /// Line 30: alph0, alph1, beta0, beta1, alphQ, betaQ (the last two 0 in the older files that do not write them).
    pub alpha: [f64; 2],
    pub beta: [f64; 2],
    pub alpha_q: f64,
    pub beta_q: f64,
    /// Line 31 (Npr = 10): alph0t, beta0t, alph1t, beta1t.
    pub thermal: Option<[f64; 4]>,
    /// Line 33: ABSCOR1, ABSCOR2.
    pub absorption: [f64; 2],
    /// Codewords: line 27-1 (seven), 27-3 (five), 30-1 (six or eight), 31-1 (four), line 33's two.
    pub sigma_codes: Vec<f64>,
    pub gamma_codes: Vec<f64>,
    pub shape_codes: Vec<f64>,
    pub thermal_codes: Vec<f64>,
    pub absorption_codes: [f64; 2],
}

/// A phase of the `.pcr` file.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PcrPhase {
    pub name: String,
    /// Line 19 flags by name (Nat, Dis, Ang, Jbt, Irf, Isy, Str, Furth, Nvk, Npr, More) and line 19-1 (Jvi, Jdi, …).
    pub flags: BTreeMap<String, f64>,
    /// Preferred-orientation direction Pr1 Pr2 Pr3.
    pub pref_direction: [f64; 3],
    pub atz: f64,
    pub symbol: String,
    /// User symmetry operators (SYMM lines, Isy ≠ 0).
    pub user_ops: Vec<String>,
    pub atoms: Vec<PcrAtom>,
    pub modes: Option<Modes>,
    pub scale: f64,
    pub shape1: f64,
    pub bov: f64,
    pub strain: [f64; 3],
    pub strain_model: i32,
    pub u: f64,
    pub v: f64,
    pub w: f64,
    pub x: f64,
    pub y: f64,
    pub gauss_size: f64,
    pub lorentz_size: f64,
    pub size_model: i32,
    /// The platelet–needle vector of LINE 36 (Size-Model ±1).
    pub size_vector: Option<[f64; 3]>,
    pub cell: [f64; 6],
    pub pref: [f64; 2],
    pub asym: [f64; 4],
    pub s_l: f64,
    pub d_l: f64,
    /// U, V, W of the second wavelength (Ratio < 0).
    pub uvw2: Option<[f64; 3]>,
    pub shp: [f64; 2],
    /// The quartic form of strain (|Str| = 1, Strain-Model ≠ 0): the S_HKL of LINE 37 in the order read, their
    /// codewords, and LINE 38's Lorentzian share ξ with its codeword (manual pp. 119–121; P25).
    pub quartic: Vec<f64>,
    pub quartic_codes: Vec<f64>,
    pub xi: Option<(f64, f64)>,
    /// Codeword lines, as written (in the order read).
    pub codes: Vec<Vec<f64>>,
    /// Codewords by line (L053): Scale, Shape1, Bov, Str1–3 (26-1); U, V, W, X, Y, GauSiz, LorSiz (27-1); the cell
    /// (29-1); Pref1, Pref2, Asy1–4, S_L, D_L (30-1).
    pub scale_codes: Vec<f64>,
    pub width_codes: Vec<f64>,
    pub cell_codes: Vec<f64>,
    pub pref_codes: Vec<f64>,
    /// Preferred orientation along several axes (Nor = 2, 3; P30): the axes after line 30-1, and their codewords in
    /// `Field::Texture`'s order (each axis's r and fraction, the last axis's r alone).
    pub texture: Option<crate::texture::MultiAxial>,
    pub texture_codes: Vec<f64>,
    /// Soft distance and angle constraints (Dis and Ang of line 19).
    pub distance_restraints: Vec<PcrDistance>,
    pub angle_restraints: Vec<PcrAngle>,
    /// Lines of the phase block after line 33, not read.
    pub extra: Vec<String>,
    pub line: usize,
    /// Its profile in a TOF pattern (Uni = 1), in place of lines 27 and 30's U … and Asy1 … (P05).
    pub tof: Option<PcrTof>,
}

impl PcrPhase {
    pub fn flag(&self, name: &str) -> i32 {
        self.flags.get(name).copied().unwrap_or(0.0) as i32
    }
}

/// One of line 12's user scattering factors (Nsc ≠ 0; manual p. 79): NAM, DFP, DFPP, ITY, and line 12b when it is
/// given: the 9 coefficients of an X-ray form factor (ITY = 0, Nsc > 0) or the 7 of a magnetic one (neutrons, ITY = 1,
/// Nsc > 0), or a table of sin θ/λ and f ended by −100 (Nsc < 0, X-rays with ITY = 0 or neutrons with ITY = 1).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PcrScatterer {
    pub name: String,
    /// Δf′ (X-rays) or the scattering length b (neutrons).
    pub dfp: f64,
    /// Δf″ (X-rays; ignored for neutrons).
    pub dfpp: f64,
    pub ity: i32,
    pub coefficients: Vec<f64>,
    pub table: Vec<(f64, f64)>,
}

/// A FullProf control file, single pattern.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Pcr {
    pub title: String,
    /// Line 4 (Job, Npr, Nph, Nba, Nex, Nsc, Nor, Dum, Iwg, Ilo, Ias, Res, Ste, Nre, Cry, Uni, Cor, Opt, Aut) and line 7
    /// (Ipr, Ppl, Ioc, Mat, Pcr, Ls1, Ls2, Ls3, NLI/Syo, Prf, Ins, Rpa, Sym, Hkl, Fou, Sho, Ana), by name.
    pub flags: BTreeMap<String, i32>,
    /// The data file named in FullProf's comment "Files => DAT-file: …", if any.
    pub data_file: Option<String>,
    pub resolution_file: Option<String>,
    pub lambda1: f64,
    pub lambda2: f64,
    pub ratio: f64,
    pub bkpos: f64,
    pub wdt: f64,
    pub cthm: f64,
    pub mur: f64,
    pub asym_limit: f64,
    pub rpolarz: f64,
    pub ncy: i32,
    /// Line 9's Eps: FullProf stops when every shift is below Eps times its esd.
    pub eps: f64,
    /// Line 9's relaxation factors R_at, R_an, R_pr, R_gl: each class's shifts are multiplied by its factor (manual
    /// LINE 9, p. 77: atomic, anisotropic, profile, global). None for a file built with `Default` (no damping).
    pub relax: Option<[f64; 4]>,
    pub thmin: f64,
    pub step: f64,
    pub thmax: f64,
    pub psd: f64,
    pub sent0: f64,
    pub background: Background,
    pub excluded: Vec<(f64, f64)>,
    /// Line 12: the user's scattering factors (Nsc ≠ 0).
    pub user_scattering: Vec<PcrScatterer>,
    /// Line 13: the number of refined parameters the file declares (FullProf's MAXS).
    pub declared_parameters: usize,
    pub zero: f64,
    pub sycos: f64,
    pub sysin: f64,
    pub lambda_refined: f64,
    pub microabsorption: Option<[f64; 3]>,
    /// The codewords of line 15's P0, Cp, τ.
    pub microabsorption_codes: [f64; 3],
    pub phases: Vec<PcrPhase>,
    /// Codeword lines of the global part, as written.
    pub codes: Vec<Vec<f64>>,
    /// Codewords of Zero, Sycos, Sysin, Lambda (line 14) and of each background value (L053).
    pub zero_codes: Vec<f64>,
    pub background_codes: Vec<f64>,
    /// MORE = 2, 3 (P23): the coefficients of the background profiles read from {data}_n.bac, and their codewords.
    pub profiles: Vec<f64>,
    pub profile_codes: Vec<f64>,
    /// Command lines (VARY, FIX, COMMANDS blocks) of the whole file, as written; not interpreted.
    pub commands: Vec<String>,
    /// General linear restraints (NLI > 0).
    pub linear_restraints: Vec<PcrLinear>,
    /// A TOF pattern (Uni = 1, P05): line 8's Iabscor (1 flat plate, 2 cylinder, 3 exponential); line 14's Dtt1, Dtt2,
    /// Dtt_1overd (0 when the file does not write it) and 2ThetaBank (degrees), Zero being `zero` and the four
    /// codewords `zero_codes`; line 16 (Npr = 10): Zerot, Dtt1t, Dtt2t, x-cross, Width and their codewords.
    pub iabscor: i32,
    pub dtt: [f64; 3],
    pub two_theta_bank: f64,
    pub thermal: Option<[f64; 5]>,
    pub thermal_codes: Vec<f64>,
}

impl Default for Background {
    fn default() -> Self {
        Background::Points(Vec::new())
    }
}

impl Pcr {
    pub fn flag(&self, name: &str) -> i32 {
        self.flags.get(name).copied().unwrap_or(0)
    }
}

const LINE4: [&str; 19] = ["Job", "Npr", "Nph", "Nba", "Nex", "Nsc", "Nor", "Dum", "Iwg", "Ilo", "Ias", "Res", "Ste", "Nre", "Cry", "Uni", "Cor", "Opt", "Aut"];
const LINE7: [&str; 17] = ["Ipr", "Ppl", "Ioc", "Mat", "Pcr", "Ls1", "Ls2", "Ls3", "NLI", "Prf", "Ins", "Rpa", "Sym", "Hkl", "Fou", "Sho", "Ana"];
const LINE19: [&str; 15] = ["Nat", "Dis", "Ang", "Pr1", "Pr2", "Pr3", "Jbt", "Irf", "Isy", "Str", "Furth", "ATZ", "Nvk", "Npr", "More"];
const LINE19_1: [&str; 10] = ["Jvi", "Jdi", "Hel", "Sol", "Mom", "Ter", "Brind", "RMua", "RMub", "RMuc"];

/// The numbers of a fixed-format field list, where FullProf may write them without a space ("1.5000-0.5000").
fn fused_numbers(s: &str) -> Vec<f64> {
    let mut spaced = String::new();
    let mut prev = ' ';
    for ch in s.chars() {
        if ch == '-' && (prev.is_ascii_digit() || prev == '.') {
            spaced.push(' ');
        }
        spaced.push(ch);
        prev = ch;
    }
    spaced.split_whitespace().map_while(|t| t.parse().ok()).collect()
}

/// The soft distance and angle constraints of a phase from the lines after its last profile line: the first Dis lines
/// that begin with two atom labels, then Ang lines with three (manual pp. 123–124).
pub fn restraints_of(extra: &[String], dis: usize, ang: usize) -> Result<(Vec<PcrDistance>, Vec<PcrAngle>), String> {
    let lines: Vec<&String> = extra.iter().filter(|l| !is_comment(l) && l.split_whitespace().next().is_some_and(|t| t.parse::<f64>().is_err())).collect();
    if lines.len() < dis + ang {
        return Err(format!("{} distance and {} angle restraints announced, {} lines found", dis, ang, lines.len()));
    }
    let mut d = Vec::new();
    for l in &lines[..dis] {
        let t: Vec<&str> = l.split_whitespace().collect();
        let v = fused_numbers(&t[2..].join(" "));
        if v.len() < 6 {
            return Err(format!("distance restraint \"{l}\": 6 numbers expected"));
        }
        d.push(PcrDistance { atoms: [t[0].into(), t[1].into()], itnum: v[0] as i32, t: [v[1], v[2], v[3]], value: v[4], sigma: v[5] });
    }
    let mut a = Vec::new();
    for l in &lines[dis..dis + ang] {
        let t: Vec<&str> = l.split_whitespace().collect();
        let v = fused_numbers(&t[3..].join(" "));
        if v.len() < 10 {
            return Err(format!("angle restraint \"{l}\": 10 numbers expected"));
        }
        a.push(PcrAngle {
            atoms: [t[0].into(), t[1].into(), t[2].into()],
            itnum: [v[0] as i32, v[1] as i32],
            t: [[v[2], v[3], v[4]], [v[5], v[6], v[7]]],
            value: v[8],
            sigma: v[9],
        });
    }
    Ok((d, a))
}

/// A line of the file: number (1-based), text, and the comment line just above it (if any).
#[derive(Clone, Debug)]
struct Line<'a> {
    no: usize,
    text: &'a str,
    comment: Option<&'a str>,
}

fn is_comment(l: &str) -> bool {
    l.trim_start().starts_with('!')
}

/// The numbers at the start of a line (stopping at the first item that is not a number, '!' or '#').
fn nums(s: &str) -> Vec<f64> {
    let s = s.split(['!', '#']).next().unwrap_or("");
    s.split_whitespace().map_while(|t| t.parse::<f64>().ok()).collect()
}

/// The column names of a FullProf comment line ("!Job Npr Nph …"), with "->Patt#…" and "/Codes" removed.
fn names(comment: &str) -> Vec<String> {
    let c = comment.trim_start().trim_start_matches('!');
    let c = c.split("->").next().unwrap_or(c);
    c.split_whitespace().filter(|w| !w.starts_with('/') && !w.starts_with('<')).map(|w| w.to_string()).collect()
}

struct Reader<'a> {
    lines: Vec<Line<'a>>,
    pos: usize,
    /// Command lines (VARY, FIX, and COMMANDS … END COMMANDS blocks), kept out of the line sequence.
    commands: Vec<String>,
}

/// A FullProf command line: VARY or FIX (global, news of 10 June 2011) or a line of a COMMANDS block.
fn is_command(t: &str) -> bool {
    let w = t.split_whitespace().next().unwrap_or("");
    matches!(w, "VARY" | "FIX" | "GVARY" | "GFIX" | "COMMANDS")
}

impl<'a> Reader<'a> {
    fn new(text: &'a [&'a str], first_no: usize) -> Reader<'a> {
        let mut lines = Vec::new();
        let mut commands = Vec::new();
        let mut in_commands = false;
        let mut last_comment: Option<&str> = None;
        for (i, &t) in text.iter().enumerate() {
            if t.trim().is_empty() {
                continue;
            }
            let up = t.trim().to_uppercase();
            if in_commands || is_command(t.trim()) {
                in_commands = (in_commands || up.starts_with("COMMANDS")) && !up.starts_with("END COMMANDS") && !up.starts_with("END_COMMANDS");
                commands.push(t.trim().to_string());
                continue;
            }
            if is_comment(t) {
                if t.trim().len() > 1 {
                    last_comment = Some(t);
                }
                continue;
            }
            lines.push(Line { no: first_no + i, text: t, comment: last_comment });
            last_comment = None;
        }
        Reader { lines, pos: 0, commands }
    }
    fn peek(&self) -> Option<&Line<'a>> {
        self.lines.get(self.pos)
    }
    fn next(&mut self, what: &str) -> Result<Line<'a>, PcrError> {
        let l = self.lines.get(self.pos).cloned().ok_or(PcrError { line: self.lines.last().map_or(0, |l| l.no), message: format!("file ends before {what}") })?;
        self.pos += 1;
        Ok(l)
    }
    fn numbers(&mut self, what: &str, at_least: usize) -> Result<(Vec<f64>, Line<'a>), PcrError> {
        let l = self.next(what)?;
        let v = nums(l.text);
        if v.len() < at_least {
            return Err(PcrError { line: l.no, message: format!("{what}: {at_least} numbers expected, found {}", v.len()) });
        }
        Ok((v, l))
    }
}

/// Flags of a line by the names of its comment (when the comment names as many columns), else by `order`.
fn by_names(l: &Line, v: &[f64], order: &[&str]) -> BTreeMap<String, f64> {
    let named = l.comment.map(names).filter(|n| n.len() >= v.len().min(order.len()) && n.iter().zip(order).filter(|(a, b)| a.eq_ignore_ascii_case(b)).count() >= 3);
    let keys: Vec<String> = match named {
        Some(n) => n,
        None => order.iter().map(|s| s.to_string()).collect(),
    };
    keys.into_iter().zip(v.iter().copied()).collect()
}

/// Read a `.pcr` file.
pub fn read(text: &str) -> Result<Pcr, PcrError> {
    let all: Vec<&str> = text.lines().collect();
    let title = all.first().copied().unwrap_or("");
    if is_multipattern(text) {
        return Err(PcrError { line: 1, message: "a multi-pattern file (NPATT): its patterns are read one by one through `split_patterns` (joint refinement, P06)".into() });
    }
    if !title.starts_with("COMM") {
        return Err(PcrError { line: 1, message: "the first line is not COMM (FullProf's single-pattern format) and no NPATT line follows".into() });
    }
    let mut p = Pcr { title: title[4..].trim().to_string(), ..Default::default() };
    for l in &all {
        if let Some(i) = l.find("DAT-file:") {
            let f = l[i + 9..].split(',').next().unwrap_or("").trim();
            if !f.is_empty() {
                p.data_file = Some(f.to_string());
            }
        }
    }
    // Split into the global part and the phase blocks at FullProf's "Data for PHASE number" comments.
    let starts: Vec<usize> = all.iter().enumerate().filter(|(_, l)| is_comment(l) && l.contains("Data for PHASE")).map(|(i, _)| i).collect();
    let global_end = starts.first().copied().unwrap_or(all.len());
    let mut r = Reader::new(&all[1..global_end], 2);

    let (v, l) = r.numbers("line 4 (Job …)", 15)?;
    for (k, x) in by_names(&l, &v, &LINE4) {
        p.flags.insert(k, x as i32);
    }
    if p.flag("Uni") != 0 && p.flag("Uni") != 1 {
        return Err(PcrError { line: l.no, message: "energy-dispersive patterns (Uni = 2) are not read yet".into() });
    }
    // Several preferred-orientation directions per phase (Nor = 2) and the multi-axial March–Dollase model of anvil
    // cells (Nor = 3), FullProf's News of 12 October 2003: their lines follow line 30 (`read_axes`, P30, L360).
    if p.flag("Nor") > 3 || p.flag("Nor") < 0 {
        return Err(PcrError { line: l.no, message: format!("Nor = {}: FullProf's preferred-orientation functions are 0 to 3", p.flag("Nor")) });
    }
    let tof = p.flag("Uni") == 1;
    if p.flag("Res") != 0 {
        let l = r.next("resolution file name")?;
        p.resolution_file = Some(l.text.split_whitespace().next().unwrap_or("").to_string());
    }
    let (v, l) = r.numbers("line 7 (Ipr …)", 10)?;
    for (k, x) in by_names(&l, &v, &LINE7) {
        p.flags.insert(k, x as i32);
    }
    // Ipr < 0: the name of a .spr file follows (written by FullProf with its own comment).
    if p.flag("Ipr") < 0 && r.peek().is_some_and(|l| nums(l.text).is_empty()) {
        r.next("spr file")?;
    }
    if p.flag("Cry") == 0 && tof {
        let (v, _) = r.numbers("line 8 (Bkpos Wdt Iabscor)", 3)?;
        (p.bkpos, p.wdt, p.iabscor) = (v[0], v[1], v[2] as i32);
    } else if p.flag("Cry") == 0 {
        let (v, _) = r.numbers("line 8 (Lambda1 …)", 9)?;
        (p.lambda1, p.lambda2, p.ratio, p.bkpos, p.wdt, p.cthm, p.mur, p.asym_limit, p.rpolarz) = (v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]);
    }
    let (v, _) = r.numbers("line 9 (NCY …)", 6)?;
    p.ncy = v[0] as i32;
    p.eps = v[1];
    p.relax = Some([v[2], v[3], v[4], v[5]]);
    if v.len() >= 11 {
        (p.thmin, p.step, p.thmax, p.psd, p.sent0) = (v[6], v[7], v[8], v[9], v[10]);
    } else if v.len() >= 9 && tof {
        // A TOF pattern's line 9 ends at TOF-max (no PSD, Sent0; p. 76).
        (p.thmin, p.step, p.thmax) = (v[6], v[7], v[8]);
    }
    let nba = p.flag("Nba");
    if p.flag("Cry") == 0 && !(-5..2).contains(&nba) {
        let mut pts = Vec::new();
        for _ in 0..nba.unsigned_abs() {
            let (v, l) = r.numbers("background point", 2)?;
            pts.push((v[0], v[1]));
            p.background_codes.push(v.get(2).copied().unwrap_or(0.0));
            p.codes.push(v[2..].to_vec());
            let _ = l;
        }
        p.background = Background::Points(pts);
    }
    for _ in 0..p.flag("Nex").max(0) {
        let (v, _) = r.numbers("excluded region", 2)?;
        p.excluded.push((v[0], v[1]));
    }
    let nsc = p.flag("Nsc");
    let xray = matches!(p.flag("Job"), 0 | 2);
    for _ in 0..nsc.unsigned_abs() {
        let l = r.next("line 12 (a user scattering factor)")?;
        let mut t = l.text.split_whitespace();
        let name = t.next().unwrap_or("").to_string();
        let v: Vec<f64> = t.map_while(|w| w.parse().ok()).collect();
        if v.len() < 3 {
            return Err(PcrError { line: l.no, message: "line 12: NAM, DFP, DFPP and ITY expected".into() });
        }
        let mut sc = PcrScatterer { name, dfp: v[0], dfpp: v[1], ity: v[2] as i32, ..Default::default() };
        match (nsc > 0, xray, sc.ity) {
            (true, true, 0) => sc.coefficients = r.numbers("line 12b (A1 B1 … C)", 9)?.0[..9].to_vec(),
            (true, false, 1) => sc.coefficients = r.numbers("line 12b (A a B b C c D)", 7)?.0[..7].to_vec(),
            (false, true, 0) | (false, false, 1) => loop {
                let (v, _) = r.numbers("line 12b (sin θ/λ, f)", 1)?;
                if v[0] == -100.0 {
                    break;
                }
                if v.len() < 2 {
                    return Err(PcrError { line: l.no, message: "line 12b: sin θ/λ and f expected, or −100".into() });
                }
                sc.table.push((v[0], v[1]));
            },
            _ => {}
        }
        p.user_scattering.push(sc);
    }
    let (v, _) = r.numbers("line 13 (number of refined parameters)", 1)?;
    p.declared_parameters = v[0].max(0.0) as usize;
    if p.flag("Cry") == 0 && tof {
        // Line 14: Zero, Dtt1, Dtt2, [Dtt_1overd,] each with its codeword, and 2ThetaBank (p. 82; the newer files write
        // Dtt_1overd, FullProf's .out: "TOF=TOF+dtt_1overD/dsp").
        let (v, l) = r.numbers("line 14 (Zero Code Dtt1 Code Dtt2 Code … 2ThetaBank)", 7)?;
        let n = if v.len() >= 9 { 4 } else { 3 };
        p.zero = v[0];
        p.dtt = [v[2], v[4], if n == 4 { v[6] } else { 0.0 }];
        p.two_theta_bank = v[2 * n];
        p.zero_codes = (0..n).map(|k| v[2 * k + 1]).collect();
        if n == 3 {
            p.zero_codes.push(0.0);
        }
        p.codes.push(p.zero_codes.clone());
        let _ = l;
        // Line 16 (Npr = 10): Hodges' thermal TOF(d) and the crossover (p. 84).
        if p.flag("Npr") == 10 {
            let (v, _) = r.numbers("line 16 (Zerot Dtt1t Dtt2t x-cross Width)", 5)?;
            p.thermal = Some([v[0], v[1], v[2], v[3], v[4]]);
            p.thermal_codes = r.numbers("line 16-1 (codewords)", 5)?.0;
            p.codes.push(p.thermal_codes.clone());
        }
    }
    if p.flag("Cry") == 0 && !tof {
        let (v, _) = r.numbers("line 14 (Zero …)", 6)?;
        (p.zero, p.sycos, p.sysin) = (v[0], v[2], v[4]);
        p.lambda_refined = v.get(6).copied().unwrap_or(0.0);
        p.codes.push(v[1..].iter().step_by(2).copied().collect());
        p.zero_codes = v[1..].iter().step_by(2).copied().take(4).collect();
        // Line 15 is read for 2θ X-ray patterns (Uni = 0, Job = 0 or 2) in Bragg–Brentano geometry (Ilo = 0) with
        // More ≠ 0 (manual p. 84), but MORE = 2, which reads the background profiles alone (P23, L353).
        let more = v.get(8).copied().unwrap_or(0.0);
        if more != 0.0 && more != 2.0 && p.flag("Ilo") == 0 && p.flag("Uni") == 0 && matches!(p.flag("Job"), 0 | 2) {
            let (m, _) = r.numbers("line 15 (micro-absorption)", 6)?;
            p.microabsorption = Some([m[0], m[2], m[4]]);
            p.microabsorption_codes = [m[1], m[3], m[5]];
            p.codes.push(vec![m[1], m[3], m[5]]);
        }
        // MORE = 2 and 3: background profiles read from files (FullProf's News of 28 July 2003, P23, L353), after line
        // 14 (and line 15 for MORE = 3): their number and coefficients, then the coefficients' codewords. FullProf
        // reads profile n from {data}_n.bac (`runs/background_terms/y2o3_more2_two.out`: "Background profile #2 read
        // from file: y2o3_2.bac").
        if more == 2.0 || more == 3.0 {
            let (c, l) = r.numbers("the number and coefficients of the background profiles", 1)?;
            let n = c[0].max(0.0) as usize;
            if c.len() < n + 1 {
                return Err(PcrError { line: l.no, message: format!("{n} background profiles: {n} coefficients expected after their number, found {}", c.len() - 1) });
            }
            p.profiles = c[1..=n].to_vec();
            let (k, _) = r.numbers("the codewords of the background profiles", n)?;
            p.profile_codes = k[..n].to_vec();
            p.codes.push(p.profile_codes.clone());
        }
    }
    if p.flag("Cry") == 0 {
        // Line 17: background coefficients in pairs of lines (values, codes) for −5 ≤ Nba ≤ 1.
        if (-5..=1).contains(&nba) {
            let mut coef = Vec::new();
            if nba == -2 {
                let (v, _) = r.numbers("FWINDOW", 1)?;
                coef.push(v[0]);
            } else {
                while let Some(l) = r.peek() {
                    let v = nums(l.text);
                    if v.is_empty() {
                        break;
                    }
                    r.pos += 1;
                    coef.extend(v);
                    let (c, _) = r.numbers("background codewords", 1)?;
                    p.background_codes.extend(&c);
                    p.codes.push(c);
                }
            }
            p.background = match nba {
                0 => Background::Polynomial(coef),
                -5 => Background::Chebyshev(coef),
                -4 => Background::Fourier(coef),
                n => Background::Other(n, coef),
            };
        }
    }
    p.commands.append(&mut r.commands);
    // Phases.
    for (k, &s) in starts.iter().enumerate() {
        let end = starts.get(k + 1).copied().unwrap_or_else(|| {
            all.iter().enumerate().skip(s).find(|(_, l)| is_comment(l) && (l.contains("2Th1/TOF1") || l.contains("2Th1"))).map_or(all.len(), |(i, _)| i)
        });
        let (ph, mut cmds) = read_phase(&all[s..end], s + 1, &p)?;
        p.phases.push(ph);
        p.commands.append(&mut cmds);
    }
    // General linear restraints (manual, NLI on line 7): after the plotting range of the pattern and the Nre lines of
    // parameters kept in a box, NLI pairs of lines: name, number of terms, value, sigma; then (coefficient, parameter
    // number) pairs.
    let nli = p.flag("NLI");
    if nli > 0 {
        let last = starts.last().copied().unwrap_or(global_end);
        let plot = all.iter().enumerate().skip(last).find(|(_, l)| is_comment(l) && l.contains("2Th1")).map_or(all.len(), |(i, _)| i);
        let rest: Vec<(usize, &str)> = all.iter().enumerate().skip(plot).filter(|(_, l)| !is_comment(l) && !l.trim().is_empty()).map(|(i, l)| (i + 1, *l)).collect();
        let skip = 1 + p.flag("Nre").max(0) as usize;
        for k in 0..nli as usize {
            let (Some(&(no, head)), Some(&(_, terms))) = (rest.get(skip + 2 * k), rest.get(skip + 2 * k + 1)) else {
                return Err(PcrError { line: all.len(), message: format!("NLI = {nli}: linear restraint {} missing", k + 1) });
            };
            let t: Vec<&str> = head.split_whitespace().collect();
            let v: Vec<f64> = t.iter().skip(1).map_while(|x| x.parse().ok()).collect();
            if t.is_empty() || v.len() < 3 {
                return Err(PcrError { line: no, message: "linear restraint: name, number of terms, value, sigma expected".into() });
            }
            let c = nums(terms);
            let n = v[0] as usize;
            if c.len() < 2 * n {
                return Err(PcrError { line: no + 1, message: format!("linear restraint {}: {n} (coefficient, parameter) pairs expected", t[0]) });
            }
            let terms = (0..n).map(|i| (c[2 * i], c[2 * i + 1] as u32)).collect();
            p.linear_restraints.push(PcrLinear { name: t[0].to_string(), value: v[1], sigma: v[2], terms });
        }
    }
    Ok(p)
}

fn read_phase(block: &[&str], first_no: usize, p: &Pcr) -> Result<(PcrPhase, Vec<String>), PcrError> {
    let mut r = Reader::new(block, first_no);
    let name_line = r.next("phase name")?;
    let mut ph = PcrPhase { name: name_line.text.trim().to_string(), line: name_line.no, ..Default::default() };
    let (v, l) = r.numbers("line 19 (Nat …)", 12)?;
    ph.flags = by_names(&l, &v, &LINE19);
    ph.pref_direction = [ph.flags.get("Pr1").copied().unwrap_or(0.0), ph.flags.get("Pr2").copied().unwrap_or(0.0), ph.flags.get("Pr3").copied().unwrap_or(1.0)];
    ph.atz = ph.flags.get("ATZ").copied().unwrap_or(0.0);
    if ph.flag("More") != 0 {
        let (v, l) = r.numbers("line 19-1 (Jvi …)", 6)?;
        for (k, x) in by_names(&l, &v, &LINE19_1) {
            ph.flags.insert(k, x);
        }
    }
    let jdi = ph.flag("Jdi");
    if jdi == 3 || jdi == 4 {
        let l = r.next("line 20 (distances)")?;
        if l.text.to_uppercase().contains("BVS") {
            let (v, _) = r.numbers("line 21 (BVS)", 2)?;
            r.next("cations")?;
            if v[1] > 0.0 {
                r.next("anions")?;
            }
        }
    }
    // Line 22: the space-group symbol in the first 20 columns.
    let l = r.next("space-group symbol")?;
    let s: String = l.text.chars().take(20).collect();
    ph.symbol = s.split("<--").next().unwrap_or("").trim().to_string();
    // A magnetic space group (FullProf's "<--Magnetic Space Group Symbol", Jbt = −6 with magnetic symmetry modes:
    // FP_tests/Ho2BaCuO5/TOF/test-tof.pcr, P2_1'/c): magnetic scattering is not calculated, and the phase is refused
    // here, where its lines first depart from a nuclear phase's (P46).
    if l.text.contains("Magnetic Space Group") || ph.flag("Jbt") == -6 {
        return Err(PcrError { line: l.no, message: format!("phase '{}': a magnetic phase in the magnetic space group {} (Jbt = {}): magnetic scattering is not calculated by the core", ph.name, l.text.split("number").next().unwrap_or("").split("<--").next().unwrap_or("").trim(), ph.flag("Jbt")) });
    }
    if ph.flag("Isy") != 0 {
        let (v, _) = r.numbers("line 23 (Nsym …)", 3)?;
        let nsym = v[0] as usize;
        while let Some(l) = r.peek() {
            let t = l.text.trim_start().to_uppercase();
            if t.starts_with("SYMM") {
                ph.user_ops.push(l.text.trim_start()[4..].trim().to_string());
            } else if !(t.starts_with("MSYM") || t.starts_with("DSYM") || t.starts_with("BAS")) {
                break;
            }
            r.pos += 1;
        }
        if ph.user_ops.len() != nsym {
            return Err(PcrError { line: l.no, message: format!("{} SYMM lines for Nsym = {nsym}", ph.user_ops.len()) });
        }
    }
    let jbt = ph.flag("Jbt");
    for _ in 0..ph.flag("Nat").max(0) {
        let l = r.next("atom")?;
        let f: Vec<&str> = l.text.split('#').next().unwrap_or("").split_whitespace().collect();
        let magnetic = jbt.abs() == 1 || jbt.abs() == 10 || jbt.abs() == 15;
        let skip = if magnetic { 2 } else { 0 }; // Mag, Vek before the coordinates
        let nv: Vec<f64> = f.iter().skip(2 + skip).map_while(|t| t.parse().ok()).collect();
        if f.len() < 2 || (nv.len() < 5 && !magnetic) {
            return Err(PcrError { line: l.no, message: "atom line: label, type, x, y, z, Biso, Occ expected".into() });
        }
        let n_t = if jbt == 0 || jbt == 2 || jbt == 6 || jbt == 3 { nv.get(7).copied().unwrap_or(0.0) as i32 } else { 0 };
        // Magnetic atom lines are fixed-width and their numbers may touch ("-0.00554-0.01836"); they are kept
        // but not used (magnetic phases are not calculated yet).
        let g = |i: usize| nv.get(i).copied().unwrap_or(0.0);
        let mut atom = PcrAtom { label: f[0].to_string(), species: f[1].to_string(), x: [g(0), g(1), g(2)], biso: g(3), occ: g(4), n_t, beta: None, codes: Vec::new(), beta_codes: None, line: l.no, rigid: (Vec::new(), Vec::new()), rigid_codes: Vec::new() };
        let (c, _) = r.numbers("atom codewords", 1)?;
        atom.codes = c.clone();
        ph.codes.push(c);
        if jbt.abs() == 4 {
            let (v, _) = r.numbers("atom line 25-3", 1)?;
            atom.rigid_codes = r.numbers("atom line 25-4", 1)?.0;
            atom.rigid = (nv.get(5..).unwrap_or(&[]).to_vec(), v);
        } else if magnetic {
            r.numbers("atom line 25-3", 1)?;
            r.numbers("atom line 25-4", 1)?;
        } else if n_t == 2 {
            let (b, _) = r.numbers("betas", 6)?;
            atom.beta = Some([b[0], b[1], b[2], b[3], b[4], b[5]]);
            atom.beta_codes = Some(r.numbers("beta codewords", 1)?.0);
        } else if n_t == 4 {
            for _ in 0..4 {
                r.numbers("form-factor parameters", 1)?;
            }
        }
        ph.atoms.push(atom);
    }
    // Symmetry modes.
    if r.peek().is_some_and(|l| l.text.trim_start().starts_with("V_MODES")) {
        let l = r.next("V_MODES")?;
        let n = nums(&l.text.trim_start()[7..]).first().copied().unwrap_or(0.0) as usize;
        let mut modes = Modes::default();
        for _ in 0..n {
            let l = r.next("mode vector")?;
            let f: Vec<&str> = l.text.split_whitespace().collect();
            if f.len() < 7 {
                return Err(PcrError { line: l.no, message: "mode vector: Nm Atm Irrep Vx Vy Vz Coeff".into() });
            }
            let num = |s: &str| s.parse::<f64>().map_err(|_| PcrError { line: l.no, message: format!("'{s}' is not a number") });
            modes.vectors.push((num(f[0])? as usize, f[1].to_string(), f[2].to_string(), [num(f[3])?, num(f[4])?, num(f[5])?], num(f[6])?));
        }
        let l = r.next("A_MODES")?;
        if !l.text.trim_start().starts_with("A_MODES") {
            return Err(PcrError { line: l.no, message: "A_MODES expected after the mode vectors".into() });
        }
        let head = nums(&l.text.trim_start()[7..]);
        let na = head.first().copied().unwrap_or(0.0) as usize;
        modes.fst = head.iter().skip(1).map(|&v| v as i32).collect();
        if let Some(l) = r.peek().filter(|l| l.text.trim_start().starts_with("Max_Amplitude")) {
            modes.max_amplitude = nums(&l.text.trim_start()[13..]).first().copied();
            r.pos += 1;
        }
        for _ in 0..na {
            let l = r.next("mode amplitude")?;
            let f: Vec<&str> = l.text.split_whitespace().collect();
            let a = f.get(1).and_then(|s| s.parse::<f64>().ok()).ok_or(PcrError { line: l.no, message: "mode amplitude: name, value, code".into() })?;
            modes.amplitudes.push((f[0].to_string(), a));
            modes.codes.push(f.get(2).and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0));
        }
        if let Some(&(m, ..)) = modes.vectors.iter().find(|v| v.0 == 0 || v.0 > na) {
            return Err(PcrError { line: l.no, message: format!("V_MODES names mode {m}, A_MODES has {na}") });
        }
        ph.modes = Some(modes);
    }
    // Integrated intensities (Cry ≠ 0 or Irf = 4): lines 26–30 of the single-crystal format (manual p. 106).
    if p.flag("Cry") != 0 || ph.flag("Irf") == 4 {
        let (v, _) = r.numbers("line 26 (Sc1 …)", 1)?;
        ph.scale = v[0];
        ph.codes.push(r.numbers("line 26-1 (codewords)", 1)?.0);
        r.numbers("line 27 (Ext1 …)", 1)?;
        ph.codes.push(r.numbers("line 27-1 (codewords)", 1)?.0);
        let (v, _) = r.numbers("line 29 (cell)", 6)?;
        ph.cell = [v[0], v[1], v[2], v[3], v[4], v[5]];
        ph.codes.push(r.numbers("line 29-1 (codewords)", 1)?.0);
        r.numbers("line 30 (Par1 …)", 1)?;
        ph.codes.push(r.numbers("line 30-1 (codewords)", 1)?.0);
        while let Some(l) = r.peek() {
            ph.extra.push(l.text.to_string());
            r.pos += 1;
        }
        return Ok((ph, r.commands));
    }
    // Lines 26–33.
    let (v, _) = r.numbers("line 26 (Scale …)", 7)?;
    (ph.scale, ph.shape1, ph.bov, ph.strain, ph.strain_model) = (v[0], v[1], v[2], [v[3], v[4], v[5]], v[6] as i32);
    ph.scale_codes = r.numbers("line 26-1 (codewords)", 1)?.0;
    ph.codes.push(ph.scale_codes.clone());
    if p.flag("Uni") == 1 {
        return read_tof_phase(r, ph, p);
    }
    let (v, _) = r.numbers("line 27 (U V W …)", 5)?;
    (ph.u, ph.v, ph.w, ph.x, ph.y) = (v[0], v[1], v[2], v[3], v[4]);
    ph.gauss_size = v.get(5).copied().unwrap_or(0.0);
    ph.lorentz_size = v.get(6).copied().unwrap_or(0.0);
    ph.size_model = v.get(7).copied().unwrap_or(0.0) as i32;
    ph.width_codes = r.numbers("line 27-1 (codewords)", 1)?.0;
    ph.codes.push(ph.width_codes.clone());
    let (v, _) = r.numbers("line 29 (cell)", 6)?;
    ph.cell = [v[0], v[1], v[2], v[3], v[4], v[5]];
    ph.cell_codes = r.numbers("line 29-1 (codewords)", 1)?.0;
    ph.codes.push(ph.cell_codes.clone());
    let (v, _) = r.numbers("line 30 (Pref1 …)", 6)?;
    (ph.pref, ph.asym) = ([v[0], v[1]], [v[2], v[3], v[4], v[5]]);
    if v.len() >= 8 {
        (ph.s_l, ph.d_l) = (v[6], v[7]);
    }
    ph.pref_codes = r.numbers("line 30-1 (codewords)", 1)?.0;
    ph.codes.push(ph.pref_codes.clone());
    read_axes(&mut r, p, &mut ph)?;
    if p.ratio < 0.0 {
        let (v, _) = r.numbers("line 32 (U2 V2 W2)", 3)?;
        ph.uvw2 = Some([v[0], v[1], v[2]]);
        r.numbers("line 32-1 (codewords)", 1)?;
    }
    let npr = if ph.flag("Npr") != 0 { ph.flag("Npr") } else { p.flag("Npr") };
    if npr == 4 || npr > 7 {
        let (v, _) = r.numbers("line 33 (Shp1 Shp2)", 4)?;
        ph.shp = [v[0], v[2]];
    }
    // LINE 36: the platelet–needle vector of Size-Model ±1 (manual p. 118).
    if ph.size_model.abs() == 1 && ph.flag("Sol") == 0 {
        let (v, _) = r.numbers("line 36 (platelet–needle vector)", 3)?;
        ph.size_vector = Some([v[0], v[1], v[2]]);
    }
    // LINE 37–38: the quartic form of strain, read where nothing comes before it but LINE 36 (no shifts).
    let size_lines = ph.size_model == 0 || ph.size_model.abs() == 1 || ph.size_model == -14;
    if ph.flag("Str").abs() == 1 && ph.strain_model != 0 && size_lines && ph.flag("Sol") == 0 {
        if let Some(n) = crate::microstructure::fullprof_strain_terms(ph.strain_model) {
            while ph.quartic.len() < n {
                let (v, _) = r.numbers("line 37 (S_HKL)", 1)?;
                let k = (n - ph.quartic.len()).min(v.len());
                ph.quartic.extend_from_slice(&v[..k]);
                ph.quartic_codes.extend(r.numbers("line 37-1 (codewords)", 1)?.0.into_iter().take(k));
            }
            if npr >= 7 {
                let (v, _) = r.numbers("line 38 (XI, code)", 2)?;
                ph.xi = Some((v[0], v[1]));
            }
        }
    }
    while let Some(l) = r.peek() {
        ph.extra.push(l.text.to_string());
        r.pos += 1;
    }
    let (dis, ang) = (ph.flag("Dis").max(0) as usize, ph.flag("Ang").max(0) as usize);
    if dis + ang > 0 {
        (ph.distance_restraints, ph.angle_restraints) = restraints_of(&ph.extra, dis, ang).map_err(|message| PcrError { line: ph.line, message })?;
    }
    Ok((ph, r.commands))
}

/// Nor = 2, 3 (FullProf's News of 12 October 2003): NPOD (and, for Nor = 3, the integration steps NOrstep), then NPOD
/// lines "h k l r code fraction code" (yb_nor2.pcr: "!  Pref-hkl  MD-r  Code  Fraction  Code"). Up to five axes; the
/// last fraction is FullProf's 1 minus the others, whatever the line says.
fn read_axes(r: &mut Reader<'_>, p: &Pcr, ph: &mut PcrPhase) -> Result<(), PcrError> {
    let nor = p.flag("Nor");
    if nor < 2 {
        return Ok(());
    }
    let (v, l) = r.numbers("the number of preferred-orientation directions (Nor = 2, 3)", 1)?;
    let n = v[0] as usize;
    if !(1..=5).contains(&n) || v[0].fract() != 0.0 {
        return Err(PcrError { line: l.no, message: format!("{} preferred-orientation directions: FullProf takes 1 to 5", v[0]) });
    }
    // NOrstep: 20 when it is not given (the News).
    let steps = if nor == 3 { v.get(1).copied().filter(|&x| x >= 1.0).unwrap_or(20.0) as usize } else { 0 };
    let mut t = crate::texture::MultiAxial { axes: Vec::new(), steps };
    let mut codes = Vec::new();
    for j in 0..n {
        let (v, l) = r.numbers("a preferred-orientation direction (h k l r code fraction code)", 7)?;
        if v[3] <= 0.0 {
            return Err(PcrError { line: l.no, message: format!("direction {}: the March–Dollase r must be above 0, not {}", j + 1, v[3]) });
        }
        t.axes.push(crate::texture::Axis { direction: [v[0], v[1], v[2]], r: v[3], fraction: v[5] });
        codes.push(v[4]);
        if j + 1 < n {
            codes.push(v[6]);
        }
    }
    t.close();
    ph.codes.push(codes.clone());
    (ph.texture, ph.texture_codes) = (Some(t), codes);
    Ok(())
}

/// Lines 27–33 of a phase in a TOF pattern (FullProf_Manual.pdf pp. 108–112), after line 26.
fn read_tof_phase(mut r: Reader<'_>, mut ph: PcrPhase, p: &Pcr) -> Result<(PcrPhase, Vec<String>), PcrError> {
    let mut t = PcrTof::default();
    let (v, _) = r.numbers("line 27 (Sig-2 Sig-1 Sig-0 …)", 3)?;
    let g = |v: &[f64], i: usize| v.get(i).copied().unwrap_or(0.0);
    (t.sigma, t.sigma_q, t.gauss_strain, t.gauss_size, t.lorentz_size_ani) = ([v[0], v[1], v[2]], g(&v, 3), g(&v, 4), g(&v, 5), g(&v, 6));
    ph.size_model = g(&v, 7) as i32;
    t.sigma_codes = r.numbers("line 27-1 (codewords)", 1)?.0;
    let (v, _) = r.numbers("line 27-2 (Gam-2 Gam-1 Gam-0 …)", 3)?;
    (t.gamma, t.lorentz_strain, t.lorentz_size) = ([v[0], v[1], v[2]], g(&v, 3), g(&v, 4));
    t.gamma_codes = r.numbers("line 27-3 (codewords)", 1)?.0;
    let (v, _) = r.numbers("line 29 (cell)", 6)?;
    ph.cell = [v[0], v[1], v[2], v[3], v[4], v[5]];
    ph.cell_codes = r.numbers("line 29-1 (codewords)", 1)?.0;
    ph.codes.push(ph.cell_codes.clone());
    let (v, _) = r.numbers("line 30 (Pref1 Pref2 alph0 beta0 alph1 beta1 …)", 6)?;
    ph.pref = [v[0], v[1]];
    (t.alpha, t.beta, t.alpha_q, t.beta_q) = ([v[2], v[4]], [v[3], v[5]], g(&v, 6), g(&v, 7));
    t.shape_codes = r.numbers("line 30-1 (codewords)", 1)?.0;
    ph.pref_codes = t.shape_codes.iter().take(2).copied().collect();
    let npr = if ph.flag("Npr") != 0 { ph.flag("Npr") } else { p.flag("Npr") };
    if npr == 10 {
        let (v, _) = r.numbers("line 31 (alph0t beta0t alph1t beta1t)", 4)?;
        t.thermal = Some([v[0], v[1], v[2], v[3]]);
        t.thermal_codes = r.numbers("line 31-1 (codewords)", 4)?.0;
    }
    read_axes(&mut r, p, &mut ph)?;
    // Line 33: "ABSCOR1 code ABSCOR2 code  ABS: ABSCOR1 ABSCOR2".
    if r.peek().is_some_and(|l| nums(l.text).len() >= 4) {
        let (v, _) = r.numbers("line 33 (ABSCOR1 code ABSCOR2 code)", 4)?;
        (t.absorption, t.absorption_codes) = ([v[0], v[2]], [v[1], v[3]]);
    }
    for c in [&t.sigma_codes, &t.gamma_codes, &t.shape_codes, &t.thermal_codes] {
        ph.codes.push(c.clone());
    }
    ph.codes.push(t.absorption_codes.to_vec());
    ph.tof = Some(t);
    while let Some(l) = r.peek() {
        ph.extra.push(l.text.to_string());
        r.pos += 1;
    }
    let (dis, ang) = (ph.flag("Dis").max(0) as usize, ph.flag("Ang").max(0) as usize);
    if dis + ang > 0 {
        (ph.distance_restraints, ph.angle_restraints) = restraints_of(&ph.extra, dis, ang).map_err(|message| PcrError { line: ph.line, message })?;
    }
    Ok((ph, r.commands))
}

/// An instrumental resolution file (`.irf`, Res ≠ 0 in the `.pcr`; manual p. 132).
#[derive(Clone, Debug, PartialEq)]
pub enum Resolution {
    /// Res = 1, 2, 3: U, V, W, X, Y, Z for each wavelength (the keyword format also gives the wavelengths, CTHM and
    /// the FCJ asymmetry S_L, D_L, kept here as written).
    Function { res: i32, uvwxyz: Vec<[f64; 6]>, keywords: BTreeMap<String, Vec<f64>> },
    /// Res = 4: a table of 2θ, H_G, H_L, linearly interpolated.
    Table(Vec<[f64; 3]>),
}

impl Resolution {
    /// Read a resolution file for the `.pcr` value of Res.
    pub fn read(text: &str, res: i32) -> Result<Resolution, String> {
        let lines: Vec<&str> = text.lines().filter(|l| !l.trim().is_empty()).collect();
        if res == 4 {
            let n = lines.get(1).and_then(|l| nums(l).first().copied()).ok_or("resolution table: number of points on line 2")? as usize;
            let mut t = Vec::new();
            for l in lines.iter().skip(2).take(n) {
                let v = nums(l);
                if v.len() < 3 {
                    return Err(format!("resolution table: '{l}'"));
                }
                t.push([v[0], v[1], v[2]]);
            }
            if t.len() != n {
                return Err(format!("resolution table: {} points for {n}", t.len()));
            }
            return Ok(Resolution::Table(t));
        }
        let mut keywords = BTreeMap::new();
        let mut uvwxyz = Vec::new();
        for (i, l) in lines.iter().enumerate() {
            if i == 0 || is_comment(l) {
                continue;
            }
            let t = l.trim_start();
            let word: String = t.chars().take_while(|c| c.is_ascii_alphabetic()).collect();
            if !word.is_empty() {
                keywords.insert(word.to_uppercase(), nums(&t[word.len()..]));
            } else {
                let v = nums(l);
                if v.len() >= 6 {
                    uvwxyz.push([v[0], v[1], v[2], v[3], v[4], v[5]]);
                }
            }
        }
        if uvwxyz.is_empty() {
            return Err("resolution file: no U V W X Y Z line".into());
        }
        Ok(Resolution::Function { res, uvwxyz, keywords })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names_of_a_comment() {
        assert_eq!(names("!Job Npr Nph Nba"), ["Job", "Npr", "Nph", "Nba"]);
        assert_eq!(names("!  Zero    Code -> Patt#1"), ["Zero", "Code"]);
        assert_eq!(nums("  20    !Number of refined"), [20.0]);
        assert_eq!(nums("La1 LA 0.1"), Vec::<f64>::new());
    }
}

/// A multi-pattern `.pcr` (FullProf_Manual.pdf pp. 65–66, 88; P06, L322) split into one single-pattern file per pattern,
/// each as FullProf's single-pattern format writes it, so that `read` reads it: the global lines (Nph Dum Ias Nre Cry Opt
/// Aut, Mat Pcr NLI Rpa Sym Sho, NCY …) joined with the pattern's own (its Job line, Ipr line, wavelengths, range,
/// background points, excluded regions, scattering factors, zero line and background coefficients), and the phases
/// that contribute to the pattern, each with its own Irf, Npr, Pr1–Pr3, Brind and Rmua–Rmuc and its own profile block
/// (scale, profile, cell, asymmetry), the atoms shared. The codewords are kept as written: a number that appears in
/// several patterns is one parameter of the joint refinement (an atom's coordinate, a cell shared by its codes).
#[derive(Clone, Debug, PartialEq)]
pub struct MultiPcr {
    pub title: String,
    /// The patterns' flags (1 refined, 0 excluded) and weights (W_PAT, normalised to sum 1).
    pub refined: Vec<bool>,
    pub weights: Vec<f64>,
    /// The data files, as written.
    pub files: Vec<String>,
    /// The single-pattern text of each pattern.
    pub texts: Vec<String>,
    /// The phases' names, and per phase whether it contributes to each pattern.
    pub phases: Vec<String>,
    pub contributions: Vec<Vec<bool>>,
}

/// Whether a `.pcr` is in the multi-pattern format: an NPATT line among its first lines (FullProf keeps the format of a
/// file whose first line is COMM, but writes NPATT when it has several patterns, as its refined pbb.pcr shows).
pub fn is_multipattern(text: &str) -> bool {
    text.lines().take(6).any(|l| l.trim_start().starts_with("NPATT"))
}

/// The pattern a FullProf comment or tagged line names ("-> Patt#: 2", "for Pattern#  2", "Patt# 1"), if any.
fn pattern_tag(line: &str) -> Option<usize> {
    let i = line.find("Patt").or_else(|| line.find("patt"))?;
    let rest = &line[i..];
    let j = rest.find('#')?;
    let digits: String = rest[j + 1..].trim_start_matches([':', ' ']).chars().take_while(|c| c.is_ascii_digit()).collect();
    digits.parse::<usize>().ok().filter(|&k| k >= 1).map(|k| k - 1)
}

fn names_of(comment: &str) -> Vec<String> {
    comment.trim_start().trim_start_matches('!').split_whitespace().take_while(|w| !w.starts_with("->") && !w.starts_with("<-")).map(str::to_string).collect()
}

fn values_of(line: &str) -> Vec<f64> {
    nums(line.split('!').next().unwrap_or(""))
}

fn number(v: f64) -> String {
    if v.fract() == 0.0 && v.abs() < 1e9 { format!("{}", v as i64) } else { format!("{v}") }
}

/// Split a multi-pattern `.pcr` into single-pattern ones.
pub fn split_patterns(text: &str) -> Result<MultiPcr, PcrError> {
    let all: Vec<&str> = text.lines().collect();
    let err = |line: usize, m: &str| PcrError { line, message: m.to_string() };
    let title = all.first().copied().unwrap_or("").trim_start_matches("COMM").trim().to_string();
    let (np_line, npatt_l) = all.iter().enumerate().find(|(_, l)| l.trim_start().starts_with("NPATT")).ok_or(err(1, "no NPATT line"))?;
    let t: Vec<&str> = npatt_l.split_whitespace().collect();
    let n: usize = t.get(1).and_then(|v| v.parse().ok()).ok_or(err(np_line + 1, "NPATT: the number of patterns expected"))?;
    let refined: Vec<bool> = (0..n).map(|k| t.get(2 + k).and_then(|v| v.parse::<i32>().ok()).is_none_or(|f| f != 0)).collect();
    let mut weights: Vec<f64> = all.iter().find(|l| l.trim_start().starts_with("W_PAT")).map(|l| l.split_whitespace().skip(1).filter_map(|v| v.parse().ok()).collect()).unwrap_or_default();
    if weights.len() != n {
        weights = vec![1.0; n];
    }
    let total: f64 = weights.iter().sum();
    let weights: Vec<f64> = weights.iter().map(|w| w / total).collect();
    // A phase's header comment is written once per pattern ("Data for PHASE number: 1 ==> Current R_Bragg for
    // Pattern# 1", then "# 2"): a start follows another only past a line that is not a comment.
    let mut phase_starts: Vec<usize> = Vec::new();
    for (i, _) in all.iter().enumerate().filter(|(_, l)| is_comment(l) && l.contains("Data for PHASE")) {
        if phase_starts.last().is_none_or(|&s| all[s..i].iter().any(|x| !is_comment(x) && !x.trim().is_empty())) {
            phase_starts.push(i);
        }
    }
    let global_end = phase_starts.first().copied().unwrap_or(all.len());
    let trailer = all.iter().enumerate().skip(phase_starts.last().copied().unwrap_or(global_end)).find(|(_, l)| is_comment(l) && l.contains("2Th1")).map_or(all.len(), |(i, _)| i);

    // The global part: each value line given to a pattern (its own tag, or its comment's) or to all.
    struct Group {
        comment: String,
        lines: Vec<String>,
        after_count: bool,
    }
    let mut global: BTreeMap<String, f64> = BTreeMap::new();
    let mut per: Vec<Vec<Group>> = (0..n).map(|_| Vec::new()).collect();
    let mut job: Vec<BTreeMap<String, f64>> = vec![BTreeMap::new(); n];
    let mut ipr: Vec<BTreeMap<String, f64>> = vec![BTreeMap::new(); n];
    let mut files = Vec::new();
    let (mut current, mut comment) = (None::<usize>, String::new());
    let mut count_line = None;
    let mut i = 1;
    while i < global_end {
        let l = all[i];
        let s = l.trim();
        if s.starts_with("NPATT") || s.starts_with("W_PAT") {
            i += 1;
            continue;
        }
        if is_comment(l) {
            let body = s.trim_start_matches('!').trim();
            if body.contains("File names of data") {
                for f in all.iter().skip(i + 1).filter(|x| !is_comment(x)).take(n) {
                    files.push(f.trim().to_string());
                }
                i += 1 + n;
                continue;
            }
            if !body.is_empty() {
                comment = l.to_string();
                // Line 16 of a TOF pattern (Npr = 10) belongs to the pattern of the line 14 above it, whose comment
                // names it ("Zerot/Code … " has no tag of its own: gabriela.pcr).
                if !(body.starts_with("Zerot") && current.is_some()) {
                    current = pattern_tag(l);
                }
            }
            i += 1;
            continue;
        }
        if s.contains("Number of refined") {
            count_line = Some(l.to_string());
            i += 1;
            continue;
        }
        let tagged = s.contains('!').then(|| pattern_tag(l)).flatten();
        match tagged.or(current) {
            Some(k) if k < n => {
                let names = names_of(&comment);
                if names.first().is_some_and(|w| w == "Job") {
                    job[k] = names.into_iter().zip(values_of(l)).collect();
                } else if names.first().is_some_and(|w| w == "Ipr") {
                    ipr[k] = names.into_iter().zip(values_of(l)).collect();
                } else {
                    match per[k].last_mut().filter(|g| g.comment == comment && !g.lines.is_empty() && tagged.is_none()) {
                        Some(g) => g.lines.push(l.to_string()),
                        None => per[k].push(Group { comment: comment.clone(), lines: vec![l.to_string()], after_count: count_line.is_some() }),
                    }
                }
            }
            Some(k) => return Err(err(i + 1, &format!("pattern {} named, but NPATT = {n}", k + 1))),
            None => {
                for (name, v) in names_of(&comment).into_iter().zip(values_of(l)) {
                    global.insert(name, v);
                }
            }
        }
        i += 1;
    }
    if files.len() != n {
        return Err(err(1, &format!("{} data files named for NPATT = {n}", files.len())));
    }

    // The phases: per phase its name, global lines, contributions, each pattern's lines, the shared body and each
    // pattern's profile block.
    struct Phase {
        name: String,
        nat: BTreeMap<String, f64>,
        jvi: BTreeMap<String, f64>,
        contributes: Vec<bool>,
        own: BTreeMap<usize, BTreeMap<String, f64>>,
        body: Vec<String>,
        blocks: BTreeMap<usize, Vec<String>>,
    }
    let mut phases = Vec::new();
    for (q, &s) in phase_starts.iter().enumerate() {
        let end = phase_starts.get(q + 1).copied().unwrap_or(trailer);
        let block = &all[s..end];
        let mut j = block.iter().position(|l| !is_comment(l)).ok_or(err(s + 1, "a phase without its name"))?;
        let name = block[j].trim().to_string();
        j += 1;
        let mut ph = Phase { name, nat: BTreeMap::new(), jvi: BTreeMap::new(), contributes: vec![true; n], own: BTreeMap::new(), body: Vec::new(), blocks: BTreeMap::new() };
        let mut last_comment = String::new();
        let mut contributing = Vec::new();
        let mut own_seen = 0usize;
        let mut in_block: Option<usize> = None;
        while j < block.len() {
            let l = block[j];
            let line_no = s + j + 1;
            if is_comment(l) && l.contains("Profile Parameters for Pattern") {
                let k = pattern_tag(l).ok_or(err(line_no, "a profile block without its pattern"))?;
                in_block = Some(k);
                ph.blocks.entry(k).or_default().push(l.to_string());
                j += 1;
                continue;
            }
            if let Some(k) = in_block {
                ph.blocks.entry(k).or_default().push(l.to_string());
                j += 1;
                continue;
            }
            if is_comment(l) {
                if !l.trim().trim_start_matches('!').trim().is_empty() {
                    last_comment = l.to_string();
                }
                if !ph.nat.is_empty() && own_seen >= 2 * contributing.len() && !l.contains("Contributions") {
                    ph.body.push(l.to_string());
                }
                j += 1;
                continue;
            }
            let names = names_of(&last_comment);
            if ph.nat.is_empty() && names.first().is_some_and(|w| w == "Nat") {
                ph.nat = names.into_iter().zip(values_of(l)).collect();
            } else if ph.jvi.is_empty() && names.first().is_some_and(|w| w == "Jvi") && contributing.is_empty() {
                ph.jvi = names.into_iter().zip(values_of(l)).collect();
            } else if last_comment.contains("Contributions") && contributing.is_empty() {
                let f = values_of(l);
                ph.contributes = (0..n).map(|k| f.get(k).is_some_and(|&x| x != 0.0)).collect();
                contributing = (0..n).filter(|&k| ph.contributes[k]).collect();
            } else if own_seen < 2 * contributing.len() {
                let k = pattern_tag(&last_comment).unwrap_or(contributing[own_seen / 2]);
                ph.own.entry(k).or_default().extend(names.into_iter().zip(values_of(l)));
                own_seen += 1;
            } else {
                ph.body.push(l.to_string());
            }
            j += 1;
        }
        if ph.nat.is_empty() {
            return Err(err(s + 1, "a phase without its Nat line"));
        }
        phases.push(ph);
    }

    // Each pattern's single-pattern text.
    let g = |name: &str, default: f64| global.get(name).copied().unwrap_or(default);
    let mut texts = Vec::new();
    for k in 0..n {
        let jb = |name: &str| job[k].get(name).copied().unwrap_or(0.0);
        let ip = |name: &str| ipr[k].get(name).copied().unwrap_or(0.0);
        if jb("Res") != 0.0 {
            return Err(err(1, &format!("pattern {}: a resolution file (Res ≠ 0) in a multi-pattern file is not read yet", k + 1)));
        }
        let nph = phases.iter().filter(|p| p.contributes[k]).count();
        let mut o = format!("COMM {title}\n! Files => DAT-file: {},  PCR-file: pattern {} of a multi-pattern file\n", files[k], k + 1);
        let line4 = [jb("Job"), jb("Npr"), nph as f64, jb("Nba"), jb("Nex"), jb("Nsc"), jb("Nor"), g("Dum", 0.0), jb("Iwg"), jb("Ilo"), g("Ias", 0.0), jb("Res"), jb("Ste"), g("Nre", 0.0), g("Cry", 0.0), jb("Uni"), jb("Cor"), g("Opt", 0.0), g("Aut", 0.0)];
        o.push_str(&format!("!{}\n{}\n", LINE4.join(" "), line4.iter().map(|v| number(*v)).collect::<Vec<_>>().join(" ")));
        let line7 = [ip("Ipr"), ip("Ppl"), ip("Ioc"), g("Mat", 0.0), g("Pcr", 1.0), ip("Ls1"), ip("Ls2"), ip("Ls3"), g("NLI", 0.0), ip("Prf"), ip("Ins"), g("Rpa", 0.0), g("Sym", 0.0), ip("Hkl"), ip("Fou"), g("Sho", 0.0), ip("Ana")];
        o.push_str(&format!("!{}\n{}\n", LINE7.join(" "), line7.iter().map(|v| number(*v)).collect::<Vec<_>>().join(" ")));
        let find = |pred: &dyn Fn(&str) -> bool| per[k].iter().find(|gr| pred(&gr.comment));
        // A TOF pattern (Uni = 1) has Bkpos, Wdt, Iabscor for its line 8 and TOF-min, Step, TOF-max for its range.
        let lambda = find(&|c| c.contains("Lambda1") || c.contains("Bkpos")).ok_or(err(1, &format!("pattern {}: no wavelength line", k + 1)))?;
        o.push_str(&format!("{}\n{}\n", lambda.comment, lambda.lines[0].split('!').next().unwrap_or("")));
        let range = find(&|c| c.contains("Thmin") || c.contains("TOF-min")).ok_or(err(1, &format!("pattern {}: no Thmin line", k + 1)))?;
        let ncy = ["NCY", "Eps", "R_at", "R_an", "R_pr", "R_gl"].map(|x| number(g(x, if x == "NCY" { 1.0 } else if x == "Eps" { 0.1 } else { 1.0 })));
        o.push_str(&format!("!NCY  Eps  R_at  R_an  R_pr  R_gl     Thmin       Step       Thmax    PSD    Sent0\n{} {}\n", ncy.join(" "), range.lines[0].split('!').next().unwrap_or("").trim()));
        let skip = |c: &str| c.contains("Lambda1") || c.contains("Thmin") || c.contains("Bkpos") || c.contains("TOF-min");
        for gr in per[k].iter().filter(|gr| !gr.after_count && !skip(&gr.comment)) {
            o.push_str(&format!("{}\n{}\n", gr.comment, gr.lines.join("\n")));
        }
        o.push_str(&format!("!\n{}\n", count_line.clone().unwrap_or_else(|| "      0    !Number of refined parameters".into())));
        for gr in per[k].iter().filter(|gr| gr.after_count) {
            o.push_str(&format!("{}\n{}\n", gr.comment, gr.lines.join("\n")));
        }
        for (q, ph) in phases.iter().enumerate().filter(|(_, p)| p.contributes[k]) {
            let own = ph.own.get(&k).cloned().unwrap_or_default();
            let nv = |m: &BTreeMap<String, f64>, x: &str, d: f64| m.get(x).copied().unwrap_or(d);
            let nat = [nv(&ph.nat, "Nat", 0.0), nv(&ph.nat, "Dis", 0.0), nv(&ph.nat, "Ang", 0.0), nv(&own, "Pr1", 0.0), nv(&own, "Pr2", 0.0), nv(&own, "Pr3", 1.0), nv(&ph.nat, "Jbt", 0.0), nv(&own, "Irf", 0.0), nv(&ph.nat, "Isy", 0.0), nv(&ph.nat, "Str", 0.0), nv(&ph.nat, "Furth", 0.0), nv(&ph.nat, "ATZ", 0.0), nv(&ph.nat, "Nvk", 0.0), nv(&own, "Npr", 0.0), 1.0];
            let jvi = [nv(&ph.jvi, "Jvi", 0.0), nv(&ph.jvi, "Jdi", 0.0), nv(&ph.jvi, "Hel", 0.0), nv(&ph.jvi, "Sol", 0.0), nv(&ph.jvi, "Mom", 0.0), nv(&ph.jvi, "Ter", 0.0), nv(&own, "Brind.", nv(&own, "Brind", 1.0)), nv(&own, "Rmua", 0.0), nv(&own, "Rmub", 0.0), nv(&own, "Rmuc", 0.0)];
            o.push_str(&format!("!-------------------------------------------------------------------------------\n!  Data for PHASE number: {:3}\n!-------------------------------------------------------------------------------\n{}\n!\n", q + 1, ph.name));
            o.push_str(&format!("!{}\n{}\n", LINE19.join(" "), nat.iter().map(|v| number(*v)).collect::<Vec<_>>().join(" ")));
            o.push_str(&format!("!{}\n{}\n", LINE19_1.join(" "), jvi.iter().map(|v| number(*v)).collect::<Vec<_>>().join(" ")));
            for l in &ph.body {
                o.push_str(l);
                o.push('\n');
            }
            for l in ph.blocks.get(&k).into_iter().flatten() {
                o.push_str(l);
                o.push('\n');
            }
        }
        for l in &all[trailer..] {
            o.push_str(l);
            o.push('\n');
        }
        texts.push(o);
    }
    Ok(MultiPcr { title, refined, weights, files, texts, phases: phases.iter().map(|p| p.name.clone()).collect(), contributions: phases.iter().map(|p| p.contributes.clone()).collect() })
}

/// Every pattern of a `.pcr` read: the one of a single-pattern file, or each of a multi-pattern file's
/// (`split_patterns`) in turn.
pub fn read_all(text: &str) -> Result<Vec<Pcr>, PcrError> {
    if is_multipattern(text) {
        split_patterns(text)?.texts.iter().map(|t| read(t)).collect()
    } else {
        Ok(vec![read(text)?])
    }
}

impl Pcr {
    /// Leave out the phases the core does not model, magnetic ones (|Jbt| = 1: magnetic structure factors are not
    /// calculated yet), and say which (P43: `drs sequential --leave-out-unmodelled` on FullProf's la0_174.pcr, as the
    /// core's own tests leave it out).
    pub fn leave_out_unmodelled(&mut self) -> Vec<String> {
        let mut said = Vec::new();
        for (k, ph) in self.phases.iter().enumerate() {
            if ph.flag("Jbt").abs() == 1 {
                said.push(format!("phase {} ({}): magnetic (Jbt = {}), whose structure factors the core does not calculate yet", k + 1, ph.name.trim(), ph.flag("Jbt")));
            }
        }
        self.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
        said
    }
}
