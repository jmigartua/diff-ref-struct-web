//! FullProf's outputs of a refinement, read for the app's views of a FullProf run (app request 17, L100): the
//! profile file `CODFIL.prf` and the reflection file `CODFIL.hkl`.
//!
//! **`.prf`** (Prf = 1 or 3, "suitable for WinPLOTR"; manual p. 9 and the flag Prf, which do not print the layout:
//! read from FullProf 8.20's files in `manuals/fullprof/`). Line 1: the title, then `Chi2:`, `CELL:` and `SPGR:` of
//! the first phase. Line 2: the number of phases, the number of points, λ1, λ2, Zero and three more numbers. Line 3:
//! the number of reflections of each phase, then counts that end with the number of excluded regions, whose (low,
//! high) follow, one per line. Then a header line (`2Theta Yobs Ycal Yobs-Ycal Backg Posr (hkl) K`, tab separated),
//! one line per point (2θ, Yobs, Ycal, Yobs − Ycal shifted down for the plot, Backg), and one line per reflection
//! (2θ, four empty columns, the tick's vertical place, `(h k l)`, and `0 w : p`, w the wavelength and p the phase).
//! An older layout (the tutorials' files of WinPLOTR and the mode course) has a column `Bragg` after Backg: the
//! reflections ride in columns 6–9 of the first points' lines, with `0 w` alone, the phase given by the tick's place
//! (the first row, 0, is phase 1).
//!
//! Backg is printed less a quarter of the background's least value over the points (L082), so its least printed value
//! is ¾ of the least background, and the background is Backg + (least Backg)/3. The reader returns the background so
//! restored; it agrees with the core's to FullProf's print (`fullprof_profiles`). Yobs and Ycal may be the data's
//! divided by a power of ten (the test's `prf_factor`): the reader returns them as printed.
//!
//! A `.prf` of integrated intensities (`Sinthet/lamb Gobs Gcal`, simulated annealing on Cry ≠ 0) has no profile and
//! is refused.
//!
//! **`.hkl`** (`CODFIL.hkl`, written with Hkl ≠ 0, and `CODFILn.hkl`, which profile matching writes and reads with
//! Irf = 2; manual pp. 10–11 and the flags Hkl and Irf). A file holds one block per phase (and pattern), headed by a line
//! `Pattern# n Phase No.: k …` (or `Phase no.: k`), then `!` lines (space group, cell), then a line naming the
//! columns, then one line of numbers per reflection. The columns differ with Hkl: ±3 gives `h k l Mult Freal Fimag 2T
//! Intensity Code`, 2 gives `h k l Mult sinT/lamb 2T Fwhm F2 sf2`, and FullProf writes `h k k` for `h k l` in some
//! headers (read as h k l). The blocks of a `CODFILn.hkl` (Irf = 2) name no columns: their first five are the
//! manual's "h, k, l, Mult, Intensity", the last two 2θ and the FWHM (they equal the `2T` and `Fwhm` of the same run's
//! Hkl = 2 file, PPH3-sls), and the sixth is not named by the manual (`c6`). A file with no block heading (an input
//! list) is one block.

use crate::data::DataError;

/// A reflection of a `.prf`.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PrfReflection {
    /// 1, 2, …
    pub phase: usize,
    pub hkl: [i32; 3],
    /// 1 or 2.
    pub wavelength: usize,
    /// 2θ (degrees), shifted by Zero as FullProf draws it.
    pub two_theta: f64,
}

/// A FullProf `.prf` file.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Prf {
    pub title: String,
    /// χ² of the title line, if printed.
    pub chi2: Option<f64>,
    pub phases: usize,
    pub wavelengths: [f64; 2],
    /// Zero, SyCos and SySin of line 2 (degrees): the ticks are at 2θ_B + Zero + SyCos·cos θ + SySin·sin 2θ.
    pub zero: f64,
    pub sycos: f64,
    pub sysin: f64,
    /// The propagation vectors of line 3's count, if any (one line each before the excluded regions).
    pub k_vectors: Vec<[f64; 3]>,
    /// The excluded regions of the header (low, high).
    pub excluded: Vec<(f64, f64)>,
    pub x: Vec<f64>,
    pub yobs: Vec<f64>,
    pub ycalc: Vec<f64>,
    /// The background, FullProf's offset added back (Backg + least Backg / 3; L082).
    pub background: Vec<f64>,
    /// What FullProf took off the background in its Backg column (a quarter of the least background).
    pub background_offset: f64,
    /// Decimals of the Ycal column (its print is half a unit of the last).
    pub ycalc_decimals: usize,
    pub reflections: Vec<PrfReflection>,
}

/// A reflection line of a `.prf` as read: 2θ, Posr, h k l, wavelength, the phase if written.
type Tick = (f64, i64, [i32; 3], usize, Option<usize>);

/// Read a FullProf `.prf`.
pub fn read_prf(text: &str) -> Result<Prf, DataError> {
    let err = |s: String| DataError(format!(".prf: {s}"));
    let lines: Vec<&str> = text.lines().collect();
    if lines.iter().any(|l| l.contains("Sinthet/lamb")) {
        return Err(err("integrated intensities (Sinthet/lamb, Gobs, Gcal): no profile".into()));
    }
    let head = lines.iter().position(|l| l.contains("Yobs") && l.contains("Ycal")).ok_or_else(|| err("no column header (2Theta Yobs Ycal …)".into()))?;
    if head < 3 {
        return Err(err("fewer than three lines before the column header".into()));
    }
    let title_line = lines[0];
    let title = title_line.split("Chi2:").next().unwrap_or("").trim().to_string();
    let chi2 = title_line.split_once("Chi2:").and_then(|(_, r)| r.split_whitespace().next()?.parse().ok());
    let nums = |l: &str| -> Vec<f64> { l.split_whitespace().map_while(|s| s.parse().ok()).collect() };
    let l2 = nums(lines[1]);
    if l2.len() < 5 {
        return Err(err(format!("line 2 is not the phases, points, λ1, λ2 and Zero: {:?}", lines[1].trim())));
    }
    let (phases, npts) = (l2[0] as usize, l2[1] as usize);
    // Line 3: reflections per phase, then counts ending with the propagation vectors and the excluded regions, whose
    // lines follow in that order (NdNi_HTO.prf: "209 281 0 1 2", one k = (1, 0, 0), two regions).
    let l3 = nums(lines[2]);
    let header: Vec<Vec<f64>> = lines[3..head].iter().map(|l| nums(l)).filter(|v| !v.is_empty()).collect();
    let nex = l3.last().map_or(0, |&v| v as usize).min(header.len());
    let nk = if l3.len() >= 2 { (l3[l3.len() - 2] as usize).min(header.len() - nex) } else { 0 };
    let k_vectors = header[header.len() - nex - nk..header.len() - nex].iter().filter(|v| v.len() >= 3).map(|v| [v[0], v[1], v[2]]).collect();
    let excluded: Vec<(f64, f64)> = header[header.len() - nex..].iter().filter(|v| v.len() >= 2).map(|v| (v[0], v[1])).collect();
    let (mut x, mut yobs, mut ycalc, mut backg) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let mut decimals = 0;
    // Reflections as (2θ, Posr, hkl, wavelength, phase if written).
    let mut ticks: Vec<Tick> = Vec::new();
    let bragg = lines[head].contains("Bragg");
    for (k, l) in lines[head + 1..].iter().enumerate() {
        let f: Vec<&str> = l.split('\t').map(|s| s.trim()).collect();
        let point = f.len() >= 5 && !f[1].is_empty();
        if point {
            let v: Vec<f64> = f[..5].iter().map(|s| s.parse()).collect::<Result<_, _>>().map_err(|_| err(format!("line {}: not five numbers", head + k + 2)))?;
            decimals = decimals.max(f[2].split_once('.').map_or(0, |(_, d)| d.len()));
            x.push(v[0]);
            yobs.push(v[1]);
            ycalc.push(v[2]);
            backg.push(v[4]);
        }
        // A reflection: on a line of its own (2θ first, four empty columns), or, where the header has a Bragg column,
        // in columns 6–9 of a point's line.
        let (t, rest) = match (bragg, point) {
            (true, _) if f.len() >= 9 && !f[5].is_empty() => (f[5], &f[6..]),
            (false, false) if f.len() >= 8 && l.contains('(') => (f[0], &f[5..]),
            _ => continue,
        };
        let two_theta: f64 = t.parse().map_err(|_| err(format!("line {}: reflection without 2θ", head + k + 2)))?;
        let posr: i64 = rest[0].parse().unwrap_or(0);
        // (3I3): "(  0  0-10)" is 0, 0, −10.
        let inner = rest[1].trim().trim_start_matches('(').trim_end_matches(')');
        let h: Vec<i32> = if inner.len() % 3 == 0 && inner.is_ascii() {
            (0..inner.len() / 3).filter_map(|i| inner[3 * i..3 * i + 3].trim().parse().ok()).collect()
        } else {
            inner.split_whitespace().filter_map(|s| s.parse().ok()).collect()
        };
        // K: "0  w :  p", or "0  w" (the phase then by its tick row, Posr).
        let (a, p) = match rest[2].split_once(':') {
            Some((a, p)) => (a, p.split_whitespace().next().and_then(|s| s.parse().ok())),
            None => (rest[2], None),
        };
        let w = a.split_whitespace().last().and_then(|s| s.parse().ok()).unwrap_or(1);
        if h.len() == 3 {
            ticks.push((two_theta, posr, [h[0], h[1], h[2]], w, p));
        }
    }
    // Without a written phase, the tick rows from the top (Posr 0, then −275, …) are phases 1, 2, … (la0190.prf:
    // 159 reflections at 0 and 270 at −275, as its line 3 counts them).
    let mut rows: Vec<i64> = ticks.iter().map(|r| r.1).collect();
    rows.sort_unstable_by(|a, b| b.cmp(a));
    rows.dedup();
    let reflections = ticks
        .into_iter()
        .map(|(two_theta, posr, hkl, wavelength, p)| PrfReflection { phase: p.unwrap_or_else(|| rows.iter().position(|&r| r == posr).unwrap() + 1), hkl, wavelength, two_theta })
        .collect();
    if x.len() != npts {
        return Err(err(format!("line 2 announces {npts} points, {} read", x.len())));
    }
    let least = backg.iter().copied().fold(f64::INFINITY, f64::min);
    let background_offset = if least.is_finite() { least / 3.0 } else { 0.0 };
    let background = backg.iter().map(|b| b + background_offset).collect();
    Ok(Prf { title, chi2, phases, wavelengths: [l2[2], l2[3]], zero: l2[4], sycos: l2.get(5).copied().unwrap_or(0.0), sysin: l2.get(6).copied().unwrap_or(0.0), k_vectors, excluded, x, yobs, ycalc, background, background_offset, ycalc_decimals: decimals, reflections })
}

/// A block of a FullProf `.hkl`: one phase (of one pattern).
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct HklBlock {
    pub pattern: Option<usize>,
    pub phase: Option<usize>,
    /// The heading line, as written.
    pub heading: String,
    /// The column names, `h`, `k`, `l` first.
    pub columns: Vec<String>,
    pub rows: Vec<Vec<f64>>,
}

impl HklBlock {
    /// The values of a column by name (as FullProf prints it: `2T`, `Fwhm`, `F2`, `Intensity`, `Mult`, …).
    pub fn column(&self, name: &str) -> Option<Vec<f64>> {
        let i = self.columns.iter().position(|c| c.eq_ignore_ascii_case(name))?;
        Some(self.rows.iter().map(|r| r.get(i).copied().unwrap_or(f64::NAN)).collect())
    }

    /// The h k l of each row.
    pub fn hkl(&self) -> Vec<[i32; 3]> {
        self.rows.iter().map(|r| [r[0] as i32, r[1] as i32, r[2] as i32]).collect()
    }
}

/// Read a FullProf `.hkl` (`CODFIL.hkl` or `CODFILn.hkl`).
pub fn read_hkl(text: &str) -> Result<Vec<HklBlock>, DataError> {
    let mut blocks: Vec<HklBlock> = Vec::new();
    let new_block = |heading: &str| {
        let after = |key: &str| heading.split_once(key).and_then(|(_, r)| r.trim_start_matches(['.', ':', ' ']).split_whitespace().next()?.trim_end_matches(':').parse().ok());
        HklBlock { pattern: after("Pattern#"), phase: after("Phase No").or_else(|| after("Phase no")), heading: heading.trim().to_string(), columns: Vec::new(), rows: Vec::new() }
    };
    for l in text.lines() {
        let t = l.trim();
        if t.is_empty() {
            continue;
        }
        let bare = t.trim_start_matches(['>', '!']).trim_start().to_lowercase();
        if t.contains("Pattern#") || bare.starts_with("phase no") {
            blocks.push(new_block(t));
            continue;
        }
        if t.starts_with('!') {
            continue;
        }
        let words: Vec<&str> = t.trim_start_matches('>').split_whitespace().collect();
        if words.len() >= 3 && words[0] == "h" && words[1] == "k" {
            if blocks.is_empty() {
                blocks.push(new_block(""));
            }
            let b = blocks.last_mut().unwrap();
            b.columns = ["h", "k", "l"].iter().map(|s| s.to_string()).chain(words[3..].iter().map(|s| s.to_string())).collect();
            continue;
        }
        let v: Vec<f64> = words.iter().map_while(|s| s.parse().ok()).collect();
        if v.len() != words.len() || v.len() < 3 {
            // A count line ("518 0 0.00 SPGr: …") or text.
            continue;
        }
        if blocks.is_empty() {
            blocks.push(new_block(""));
        }
        blocks.last_mut().unwrap().rows.push(v);
    }
    for b in &mut blocks {
        if b.columns.is_empty() {
            let width = b.rows.iter().map(|r| r.len()).max().unwrap_or(0);
            let mut c: Vec<String> = ["h", "k", "l", "Mult", "Intensity"].iter().map(|s| s.to_string()).collect();
            if width == 8 && b.pattern.is_some() {
                // CODFILn.hkl of profile matching (Irf = 2): …, c6, 2θ, FWHM.
                c.extend(["c6", "2T", "Fwhm"].iter().map(|s| s.to_string()));
            }
            while c.len() < width {
                c.push(format!("c{}", c.len() + 1));
            }
            c.truncate(width.max(3));
            b.columns = c;
        }
    }
    // A list may be empty (la0177_2.hkl: a magnetic phase's heading and column names only).
    if blocks.is_empty() {
        return Err(DataError(".hkl: neither a block heading, column names nor a reflection".into()));
    }
    Ok(blocks)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn prf_reads_profile_reflections_and_background() {
        let t = " T   Chi2:  1.5  CELL: 1 1 1 90 90 90 SPGR: P 1\n  1      3     1.54060     1.54440     0.00250     0.00000     0.00000    0\n   2     0     1\n       0.00000      14.45000\n 2Theta\tYobs\tYcal\tYobs-Ycal\tBackg\tPosr\t(hkl)\tK\n     15.0000\t    25.0\t    19.2\t -2806.2\t    15.0\n     15.0100\t    22.0\t    19.2\t -2809.2\t    16.0\n     15.0200\t    22.0\t    19.2\t -2809.2\t    17.0\n     16.7086\t        \t        \t        \t        \t       0\t(  2  0  0)\t  0  1 :  1\n     16.75\t        \t        \t        \t        \t    -307\t(  2 -1  0)\t  0  2 :  2\n";
        let p = read_prf(t).unwrap();
        assert_eq!((p.phases, p.x.len(), p.excluded.clone(), p.chi2), (1, 3, vec![(0.0, 14.45)], Some(1.5)));
        assert_eq!(p.background, vec![20.0, 21.0, 22.0]);
        assert_eq!(p.reflections[1], PrfReflection { phase: 2, hkl: [2, -1, 0], wavelength: 2, two_theta: 16.75 });
        assert_eq!(p.ycalc_decimals, 1);
    }

    #[test]
    fn hkl_blocks_by_heading() {
        let t = " Pattern#  1 Phase No.:  1 phase 1:   2 reflections\n!  SPGR: I a -3;  CELL: 10 10 10 90 90 90\n   h   k   l    Mult    Freal       Fimag        2T            Intensity  Code\n   2   0   0       6    -31.6749     -1.3657     16.7061          3.2129    0\n   2   1   1      24    130.0799      6.5167     20.4974        143.1689    0\n";
        let b = read_hkl(t).unwrap();
        assert_eq!((b.len(), b[0].phase, b[0].pattern), (1, Some(1), Some(1)));
        assert_eq!(b[0].column("2T").unwrap(), vec![16.7061, 20.4974]);
        assert_eq!(b[0].hkl()[1], [2, 1, 1]);
        let le_bail = " Pattern#  1 Phase No:   1  phase 1:            Lambda:  1.001000 CELL:    11.4305\n     518       0    0.00  SPGr: P 2/m  <-- The number\n   1   1   0    4        704.811          5.770       6.3125       0.0234\n";
        let b = read_hkl(le_bail).unwrap();
        assert_eq!(b[0].columns, ["h", "k", "l", "Mult", "Intensity", "c6", "2T", "Fwhm"]);
        assert_eq!(b[0].column("Intensity").unwrap(), vec![704.811]);
    }
}

/// A refined parameter as FullProf's `.out` lists it after its last cycle ("Parameter number n : name value( +/- esd
/// )", manual p. 9): FullProf's own values, in its rules.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct OutParameter {
    pub number: usize,
    pub name: String,
    pub value: f64,
    pub esd: f64,
}

/// One of the `.out`'s statistics blocks of its last cycle: the R factors "not corrected for background" (%), χ² as
/// "R-Factors:" prints it (two decimals), the expected Rwp, N − P + C, and the conventional (background-free) Rp,
/// Rwp, Re and χ² (three decimals), "Values obtained using Ynet, but true sigma(y)".
#[derive(Clone, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct OutStatistics {
    pub rp: f64,
    pub rwp: f64,
    pub rexp: f64,
    pub chi2: f64,
    pub dof: usize,
    pub durbin_watson: f64,
    pub rp_conventional: f64,
    pub rwp_conventional: f64,
    pub rexp_conventional: f64,
    pub chi2_conventional: f64,
}

/// What a comparison reads of a FullProf `.out` (the app's request 49): the program's banner, the refined parameters
/// with their esd's after the last cycle, the two statistics blocks of the last cycle (over all the points not
/// excluded, and over those with Bragg contributions, by which FullProf scales its esd's: manual eqs. 3.7–3.8), R_Bragg
/// and R_F per phase, the cycles run and whether FullProf's convergence test was met. Layout read from FullProf 8.20's
/// files in `manuals/fullprof/runs/` (`*_conv.out`).
#[derive(Clone, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Out {
    /// "FullProf.2k (Version 8.20 - Feb2025-ILL JRC)".
    pub banner: String,
    pub parameters: Vec<OutParameter>,
    pub all_points: Option<OutStatistics>,
    pub bragg_points: Option<OutStatistics>,
    /// Per phase: R_Bragg and R_F (%).
    pub phases: Vec<(f64, f64)>,
    pub cycles: usize,
    pub converged: bool,
}

fn statistics_block(lines: &[&str]) -> OutStatistics {
    let nums = |s: &str| -> Vec<f64> { s.split_whitespace().map_while(|t| t.parse().ok()).collect() };
    let after = |l: &str, key: &str| -> Vec<f64> { l.split_once(key).map_or(Vec::new(), |(_, r)| nums(r)) };
    let mut s = OutStatistics::default();
    for l in lines {
        if l.contains("=> R-Factors:") {
            let v = after(l, "R-Factors:");
            (s.rp, s.rwp) = (v.first().copied().unwrap_or(f64::NAN), v.get(1).copied().unwrap_or(f64::NAN));
            s.chi2 = after(l, "Chi2:").first().copied().unwrap_or(f64::NAN);
            s.durbin_watson = after(l, "DW-Stat.:").first().copied().unwrap_or(f64::NAN);
        } else if l.contains("=> Expected :") {
            s.rexp = after(l, "Expected :").first().copied().unwrap_or(f64::NAN);
        } else if l.contains("=> N-P+C:") {
            s.dof = after(l, "N-P+C:").first().map_or(0, |&v| v as usize);
        } else if l.contains("Conventional Rietveld Rp,Rwp,Re and Chi2:") {
            let v = after(l, "Chi2:");
            if v.len() >= 4 {
                (s.rp_conventional, s.rwp_conventional, s.rexp_conventional, s.chi2_conventional) = (v[0], v[1], v[2], v[3]);
            }
        }
    }
    s
}

/// Read a FullProf `.out` (see `Out`).
pub fn read_out(text: &str) -> Result<Out, DataError> {
    let lines: Vec<&str> = text.lines().collect();
    let banner = lines
        .iter()
        .find_map(|l| l.split_once("PROGRAM ").filter(|_| l.contains("FullProf")).map(|(_, r)| r.trim().trim_end_matches('*').trim().to_string()))
        .ok_or_else(|| DataError(".out: no 'PROGRAM FullProf' banner".into()))?;
    let mut o = Out { banner, ..Default::default() };
    // The parameters: the last listing of each number.
    let mut by_number: std::collections::BTreeMap<usize, OutParameter> = std::collections::BTreeMap::new();
    for l in lines.iter().filter(|l| l.contains("Parameter number") && l.contains("+/-")) {
        let Some((head, rest)) = l.split_once(':') else { continue };
        let Some(number) = head.split_whitespace().last().and_then(|t| t.parse().ok()) else { continue };
        let mut t = rest.split_whitespace();
        let Some(name) = t.next() else { continue };
        let v: Vec<f64> = t.collect::<Vec<_>>().join(" ").replace(['(', ')'], " ").replace("+/-", " ").split_whitespace().filter_map(|s| s.parse().ok()).collect();
        if v.len() >= 2 {
            by_number.insert(number, OutParameter { number, name: name.to_string(), value: v[0], esd: v[1] });
        }
    }
    o.parameters = by_number.into_values().collect();
    let all = lines.iter().rposition(|l| l.contains("RELIABILITY FACTORS WITH ALL NON-EXCLUDED POINTS"));
    let bragg = lines.iter().rposition(|l| l.contains("RELIABILITY FACTORS FOR POINTS WITH BRAGG"));
    if let Some(a) = all {
        let end = bragg.filter(|&b| b > a).unwrap_or(lines.len());
        o.all_points = Some(statistics_block(&lines[a..end]));
    }
    if let Some(b) = bragg {
        let end = lines[b..].iter().position(|l| l.contains("Global user-weigthed") || l.contains("=> Phase:")).map_or(lines.len(), |k| b + k);
        o.bragg_points = Some(statistics_block(&lines[b..end]));
    }
    let first_num = |l: &str, key: &str| -> Option<f64> { l.split_once(key)?.1.split_whitespace().next()?.parse().ok() };
    // R_Bragg and R_F: the summary of the last cycle ("Bragg R-factor: … Vol: … Fract(%)", "Rf-factor= …") when
    // printed, otherwise the per-phase lines after the last statistics block.
    let summary: Vec<f64> = lines.iter().filter(|l| l.contains("Vol:")).filter_map(|l| first_num(l, "Bragg R-factor:")).collect();
    if !summary.is_empty() {
        let rf: Vec<f64> = lines.iter().filter_map(|l| first_num(l, "Rf-factor=")).collect();
        let n = summary.len();
        let rf = &rf[rf.len().saturating_sub(n)..];
        o.phases = summary.iter().enumerate().map(|(k, &r)| (r, rf.get(k).copied().unwrap_or(f64::NAN))).collect();
    } else if let Some(b) = bragg {
        let rb: Vec<f64> = lines[b..].iter().filter_map(|l| first_num(l, "Bragg R-factor:")).collect();
        let rf: Vec<f64> = lines[b..].iter().filter_map(|l| first_num(l, "RF-factor     :").or_else(|| first_num(l, "RF-factor :"))).collect();
        o.phases = rb.iter().enumerate().map(|(k, &r)| (r, rf.get(k).copied().unwrap_or(f64::NAN))).collect();
    }
    o.cycles = lines.iter().filter_map(|l| first_num(l, "CYCLE No.:")).next_back().map_or(0, |v| v as usize);
    o.converged = lines.iter().any(|l| l.contains("Convergence reached"));
    Ok(o)
}
