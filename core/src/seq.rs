//! FullProf's `.seq`, the results of a sequential refinement read as a table (Plots' PL11, R32).
//!
//! FullProf writes a `.seq` in a sequential refinement when Rpa = −2: "the summary of effectively refined parameters
//! together with their symbolic names" (FullProf news, 10 February 2004, `manuals/fullprof/news/Fullprof_News_2004.htm`);
//! a later run with the same code appends to it unless Rpa = −3 (news, 14 February 2008, `Fullprof_News_2008.htm`), so
//! a file may hold the same data file twice. WinPLOTR plots its parameters against the file order or the external
//! parameter (WinPLOTR deck, `manuals/fullprof/1-WinPLOTR/FPS_dec2015_WinPLOTR_TR.pdf` pp. 19–23). Neither document
//! gives the layout; it is read from the example `1-WinPLOTR/2_FP_sequential_refinements/LaMnO3_D20/la0.seq` (71
//! blocks) and checked against its run's data files and `.sum` files (test `plots`):
//!
//! ```text
//! NEW_REFINEMENT       la0120          the data file's code
//! NUMPAR    46                         the number of refined parameters listed below
//! N_PATT     1
//! NPHASE     2
//! PARVAL     0     68.8610    73430    the file's index, the external parameter and the numor
//! GLOBAL_CHI_2       3.86369
//!      1         R_Bragg(%):     6.6040678         5.0635004     phase 1: R_Bragg and R_F
//!      1   Unit_Cell_Volume:     245.01285        0.56278920E-02   the volume and its σ
//!      1    Weight_Fraction:     100.00000         1.8546307      the weight fraction (%) and its σ
//!      1                                    Scale_ph1_pat1     2.2853527        0.29970573E-01
//! ```
//!
//! PARVAL's second and third numbers are the data file's sample temperature and numor (`la0120.dat`: 68.861 K,
//! numor 73430; all 71 blocks agree). R_Bragg's second number is the `.sum`'s Rf-factor for a nuclear phase (la0_120:
//! 6.60 and 5.06); for the magnetic phase the `.sum` gives only the magnetic R-factor (4.47, the first number) and the
//! second is not identified. Lines the reader does not know are kept in `unread`, never dropped.

/// One phase's line group in a block.
#[derive(Clone, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SeqPhase {
    pub phase: usize,
    pub r_bragg: Option<f64>,
    /// R_Bragg's second number: the `.sum`'s Rf-factor for a nuclear phase.
    pub r_f: Option<f64>,
    pub volume: Option<(f64, f64)>,
    pub weight_fraction: Option<(f64, f64)>,
}

/// One refined parameter: FullProf's index, its symbolic name, value and σ.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SeqParameter {
    pub index: usize,
    pub name: String,
    pub value: f64,
    pub sigma: f64,
}

/// One refinement of the sequence.
#[derive(Clone, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SeqBlock {
    /// NEW_REFINEMENT: the data file's code.
    pub name: String,
    pub numpar: Option<usize>,
    pub patterns: Option<usize>,
    pub phase_count: Option<usize>,
    /// PARVAL: the file's index in the sequence, the external parameter, the numor.
    pub index: Option<i64>,
    pub parameter: Option<f64>,
    pub numor: Option<i64>,
    /// GLOBAL_CHI_2.
    pub chi2: Option<f64>,
    pub phases: Vec<SeqPhase>,
    pub parameters: Vec<SeqParameter>,
}

/// A `.seq` file: its blocks in the file's order, and the lines not understood (line number from 1, text).
#[derive(Clone, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Seq {
    pub blocks: Vec<SeqBlock>,
    pub unread: Vec<(usize, String)>,
}

fn number(s: &str) -> Option<f64> {
    s.replace(['D', 'd'], "E").parse().ok()
}

/// Reads a `.seq` (see the module).
pub fn read(text: &str) -> Result<Seq, String> {
    let mut seq = Seq::default();
    for (n, raw) in text.lines().enumerate() {
        let line = raw.trim();
        if line.is_empty() {
            continue;
        }
        let w: Vec<&str> = line.split_whitespace().collect();
        if w[0] == "NEW_REFINEMENT" {
            seq.blocks.push(SeqBlock { name: w.get(1).unwrap_or(&"").to_string(), ..Default::default() });
            continue;
        }
        let Some(b) = seq.blocks.last_mut() else {
            seq.unread.push((n + 1, raw.to_string()));
            continue;
        };
        let int = |k: usize| w.get(k).and_then(|v| v.parse::<i64>().ok());
        let real = |k: usize| w.get(k).and_then(|v| number(v));
        let ok = match w[0] {
            "NUMPAR" => int(1).map(|v| b.numpar = Some(v as usize)).is_some(),
            "N_PATT" => int(1).map(|v| b.patterns = Some(v as usize)).is_some(),
            "NPHASE" => int(1).map(|v| b.phase_count = Some(v as usize)).is_some(),
            "PARVAL" => {
                (b.index, b.parameter, b.numor) = (int(1), real(2), int(3));
                b.index.is_some() && b.parameter.is_some()
            }
            "GLOBAL_CHI_2" => real(1).map(|v| b.chi2 = Some(v)).is_some(),
            _ => match (w[0].parse::<usize>(), w.len()) {
                (Ok(k), 4) if w[1].ends_with(':') => {
                    let (v, e) = (real(2), real(3));
                    let ph = match b.phases.iter().position(|p| p.phase == k) {
                        Some(i) => &mut b.phases[i],
                        None => {
                            b.phases.push(SeqPhase { phase: k, ..Default::default() });
                            b.phases.last_mut().unwrap()
                        }
                    };
                    match (w[1], v, e) {
                        ("R_Bragg(%):", Some(v), Some(e)) => (ph.r_bragg, ph.r_f) = (Some(v), Some(e)),
                        ("Unit_Cell_Volume:", Some(v), Some(e)) => ph.volume = Some((v, e)),
                        ("Weight_Fraction:", Some(v), Some(e)) => ph.weight_fraction = Some((v, e)),
                        _ => {}
                    }
                    matches!(w[1], "R_Bragg(%):" | "Unit_Cell_Volume:" | "Weight_Fraction:") && v.is_some() && e.is_some()
                }
                (Ok(k), 4) => match (real(2), real(3)) {
                    (Some(value), Some(sigma)) => {
                        b.parameters.push(SeqParameter { index: k, name: w[1].to_string(), value, sigma });
                        true
                    }
                    _ => false,
                },
                _ => false,
            },
        };
        if !ok {
            seq.unread.push((n + 1, raw.to_string()));
        }
    }
    if seq.blocks.is_empty() {
        return Err("no NEW_REFINEMENT block: not a FullProf .seq".into());
    }
    Ok(seq)
}

impl Seq {
    /// Every parameter's name, in the order of first appearance.
    pub fn names(&self) -> Vec<String> {
        let mut out: Vec<String> = Vec::new();
        for p in self.blocks.iter().flat_map(|b| &b.parameters) {
            if !out.contains(&p.name) {
                out.push(p.name.clone());
            }
        }
        out
    }

    /// A parameter's value and σ in each block (None where the block does not refine it).
    pub fn column(&self, name: &str) -> Vec<Option<(f64, f64)>> {
        self.blocks.iter().map(|b| b.parameters.iter().find(|p| p.name == name).map(|p| (p.value, p.sigma))).collect()
    }

    /// The blocks whose NUMPAR differs from the parameters listed (none in a whole file).
    pub fn inconsistent(&self) -> Vec<usize> {
        (0..self.blocks.len()).filter(|&k| self.blocks[k].numpar.is_some_and(|n| n != self.blocks[k].parameters.len())).collect()
    }

    /// The data files refined more than once (a run appended to an earlier one, news of 14 February 2008).
    pub fn repeated(&self) -> Vec<String> {
        let mut out = Vec::new();
        for (k, b) in self.blocks.iter().enumerate() {
            if self.blocks[..k].iter().any(|c| c.name == b.name) && !out.contains(&b.name) {
                out.push(b.name.clone());
            }
        }
        out
    }
}
