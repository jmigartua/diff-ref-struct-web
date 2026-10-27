//! The calculated powder pattern for constant-wavelength data, from the papers in `manuals/papers/` and FullProf's
//! manual; where FullProf's own calculation departs from them, `Convention` says, rule by rule, whether FullProf's rule
//! is followed (L043, L046, L050). The list below names FullProf's rules where they differ; `Convention` has the table.
//!
//! y_c(i) = Σ_φ S_φ Σ_h I_φh Ω(2θ_i − 2θ_φh) + b_i (FullProf manual eq. 3.1, p. 41), with:
//! - positions 2θ_h = 2 asin(λ/2d_h) + Zero + Sycos·cosθ + Sysin·sin2θ (degrees; manual pp. 56–57 and 83), for λ1 and,
//!   when Ratio ≠ 0 and λ2 ≠ λ1, for λ2 with intensity × Ratio;
//! - integrated intensities I_h = m_h · Lp(θ_h) · F_h² · exp(−2 B_ov s²) · P_h (eq. 3.2, p. 41; the overall
//!   temperature factor applied to F² as exp(−B_ov |h|²/2) = exp(−2 B_ov s²), p. 58), with F_h² the structure factor
//!   of `structure_factor` (FullProf's printed F², L030), or, for profile matching with intensities from a file
//!   (Jbt = 2, Irf = 2), the intensities of that file, taken as the total of both wavelengths and without the scale
//!   factor (as FullProf's profile-matching files show; not stated in the manual);
//! - Lorentz–polarisation for Debye–Scherrer or Bragg–Brentano geometry: Lp = 1/(2 sin²θ cosθ) for neutrons and
//!   ½(1 + Cthm·cos²2θ)/(2 sin²θ cosθ) for X-rays, eq. 3.26 (p. 53) with K = 0.5; FullProf computes twice that: its files
//!   agree only with this (PPH3-sls.prf, Cthm = 0, to the print; y2o3_structure.pcr, Cthm = 0.78); synchrotron
//!   (Ilo = 3): eq. 3.26 with K = Rpolarz, Cthm = 0 taken as 1 (checked on FullProf's joint run prbaco, L322);
//!   transmission flat plate (Ilo = 2), eqs. 3.27–3.28;
//! - widths: H² = (U tan θ + V) tan θ + W + I_G/cos²θ (eq. 3.20) for the pseudo-Voigt (Npr = 5, η = η0 + X·2θ, eq.
//!   3.21), and H_G² as H² with H_L = X tan θ + Y/cosθ for Thompson–Cox–Hastings (Npr = 7, eqs. 3.24–3.25), combined
//!   by eqs. 3.16–3.17; with a resolution file (Res ≠ 0) its instrumental H_G², H_L are added to the sample ones; all
//!   at the Bragg angle, without Zero and the displacement terms (L039);
//! - asymmetry below AsyLim: Finger–Cox–Jephcoat with S/L, H/L for Npr = 7 (FullProf's reading of its pair with
//!   `Convention::axial_divergence`), and the Bérar–Baldinozzi multiplier (the paper's; FullProf's F_b with
//!   `Convention::berar_baldinozzi`: `profile::berar_baldinozzi`) when Asy1–Asy4 are not zero, on top of FCJ when a
//!   phase has both;
//! - absorption (μR of line 8): the transmission of a cylinder, eq. 3.33 (p. 55), or of a flat plate in transmission
//!   (Ilo = 2), eq. 3.34;
//! - preferred orientation: the exponential of Dollase (1986), Table 1, no. 5 (Nor = 0; FullProf's eq. 3.31 has the opposite sign of G1) or March–Dollase (Nor = 1, eq. 3.32), averaged over the equivalent reflections;
//! - each peak computed within ±Wdt FWHM of its position (p. 75), the FCJ profile too; FullProf counts the width
//!   without the phase's X·tanθ for TCH (`Model::window_width`, L035) and takes one more low point when X ≠ 0 (L042);
//!   no peak inside the excluded regions, where FullProf's calculated profile is the background alone;
//! - every reflection whose window reaches the data; FullProf's: from the first point up to 5° past the last, or to
//!   `Phase::reflection_limit` when it is set;
//! - background: linear interpolation between points (constant beyond the end points), the polynomial of eq. 3.4,
//!   Chebyshev polynomials over the data range (Nba = −5), or the Fourier cosine series of eq. 3.4′ (Nba = −4).

use crate::cell::Cell;
use crate::data::Pattern;
use crate::pcr::{Background, Pcr, PcrPhase, Resolution};
use crate::profile::{berar_baldinozzi, fullprof_axial, gaussian, lorentzian, pseudo_voigt, taper, tch, Fcj};
use crate::reflections::generate;
use crate::scattering::Radiation;
use crate::spacegroup::SpaceGroup;
use crate::structure_factor::{Atom, Structure};
use std::path::{Path, PathBuf};

/// The peak shape of a phase.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Shape {
    Gaussian,
    Lorentzian,
    PseudoVoigt,
    ThompsonCoxHastings,
}

impl Shape {
    pub fn from_npr(npr: i32) -> Result<Shape, String> {
        match npr {
            0 => Ok(Shape::Gaussian),
            1 => Ok(Shape::Lorentzian),
            5 => Ok(Shape::PseudoVoigt),
            7 => Ok(Shape::ThompsonCoxHastings),
            // The pseudo-Voigt with axial divergence (Npr = 12, manual p. 69): Npr 5's η0 and X, with Finger, Cox and
            // Jephcoat's S_L, D_L (line 33) below AsyLim (gabriela.pcr's synchrotron patterns, P05).
            12 => Ok(Shape::PseudoVoigt),
            // The TOF functions (P05): the Thompson–Cox–Hastings pseudo-Voigt that the exponentials are convolved
            // with (`tof`).
            9 | 10 => Ok(Shape::ThompsonCoxHastings),
            13 => Err("profile function Npr = 13 (Ikeda–Carpenter ⊗ pseudo-Voigt, TOF): no FullProf run in the repository uses it, not calculated (L341)".into()),
            n => Err(format!("profile function Npr = {n} is not calculated yet (0, 1, 5, 7, 12 and the TOF 9, 10 are)")),
        }
    }
}

/// Where the intensities of a phase come from.
#[derive(Clone, Debug)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Intensities {
    /// From the structure (Rietveld).
    Structure(Box<Structure>, Vec<crate::scattering::Scatterer>),
    /// From a list (h k l, multiplicity, intensity), as profile matching reads it (Jbt = 2, Irf = 2).
    List(Vec<([i32; 3], usize, f64)>),
}

/// A phase of the pattern model.
#[derive(Clone, Debug)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Phase {
    pub name: String,
    pub group: SpaceGroup,
    pub cell: Cell,
    pub intensities: Intensities,
    pub scale: f64,
    pub bov: f64,
    pub shape: Shape,
    pub eta0: f64,
    pub u: f64,
    pub v: f64,
    pub w: f64,
    pub x: f64,
    pub y: f64,
    pub gauss_size: f64,
    pub asym: [f64; 4],
    pub s_l: f64,
    pub d_l: f64,
    pub pref: [f64; 2],
    pub pref_direction: [f64; 3],
    /// Largest 1/d (Å⁻¹) of the reflections of this phase, when it is set. FullProf also generates reflections past
    /// Thmax, as far as each phase's own (sample) broadening asks, and prints the limit in its `.out` ("S(1/d)max");
    /// its rule is not in the manual (L035). Unset, every reflection whose 2θ lies within 5° of the pattern is used.
    pub reflection_limit: Option<f64>,
    /// The intensities are extracted from the data at every cycle (Le Bail; FullProf's profile matching with
    /// constant scale, Jbt = 2): `Model::le_bail`.
    pub extract: bool,
    /// Quantitative analysis: FullProf's ATZ = Z·M_w·f²/t (manual p. 91) and the Brindley coefficient t (p. 93).
    pub atz: f64,
    pub brindley: f64,
    /// The weight fraction of this phase weighed into the sample as an internal standard (P33): the amorphous content
    /// follows from it (Orlhac et al. 2001, eqs. 3–4; `qpa::amorphous_content`). At most one phase carries it.
    #[cfg_attr(feature = "serde", serde(default))]
    pub internal_standard: Option<f64>,
    /// The relative range about |q| within which a valence sum is taken as plausible (`bonds::VALENCE_RANGE`, 0.10,
    /// when None; the control file's `valence_sum_range`, Explore's E18c).
    #[cfg_attr(feature = "serde", serde(default))]
    pub valence_sum_range: Option<f64>,
    /// Symmetry modes (R29, L090; FullProf's Jbt = 6): the atoms' positions are their reference positions plus the
    /// modes' displacements, Σ Q·c·V, and the amplitudes Q are refined in place of the coordinates. The structure's
    /// atoms hold the positions the modes give (`Modes::positions`), kept so by `Model::set`.
    #[cfg_attr(feature = "serde", serde(default))]
    pub modes: Option<crate::modes::Modes>,
    /// Rigid groups (P41, FullProf's Jbt = 4): each group's atoms placed by its origin, orientation and internal
    /// coordinates (`rigid_body`), which are refined in place of their coordinates. The structure's atoms hold the
    /// positions the groups give, kept so by `Model::set`. Empty: none.
    #[cfg_attr(feature = "serde", serde(default))]
    pub rigid: Vec<crate::rigid_body::Group>,
    /// U, V, W of the second wavelength's peaks, when they have their own (FullProf's Ratio < 0, line 32: a λ/2
    /// contamination, D1B in the mode course's case D; L098). Not refined.
    #[cfg_attr(feature = "serde", serde(default))]
    pub uvw2: Option<[f64; 3]>,
    /// The parent structure of a phase refined by its coordinates (Publish P17): its group and the transformation, as
    /// a mode phase holds it in `modes.parent`, so that its cell can be reduced to the parent's (`geometry::parent_cell`).
    #[cfg_attr(feature = "serde", serde(default))]
    pub parent: Option<crate::modes::ModeParent>,
    /// Anisotropic microstrain, Stephens' quartic form in hkl with its Lorentzian share (P25, L211): widths that
    /// depend on the reflection, added to the phase's U and X terms (`Model::peak_width`).
    #[cfg_attr(feature = "serde", serde(default))]
    pub quartic_strain: Option<crate::microstructure::QuarticStrain>,
    /// Anisotropic Lorentzian size by shape (platelets, needles) or by the antiphase rule (P25, L267).
    #[cfg_attr(feature = "serde", serde(default))]
    pub shape_size: Option<crate::microstructure::ShapeSize>,
    /// The phase's peak shape in a TOF pattern (P05; `Model::tof`): back-to-back exponentials convolved with the
    /// Thompson–Cox–Hastings pseudo-Voigt (`shape`), in place of U … Y and the asymmetry.
    #[cfg_attr(feature = "serde", serde(default))]
    pub tof: Option<crate::tof::Profile>,
    /// Preferred orientation along several axes (P30, FullProf's Nor = 2, 3): when set, it replaces `pref` and
    /// `pref_direction` (`texture`).
    #[cfg_attr(feature = "serde", serde(default))]
    pub texture: Option<crate::texture::MultiAxial>,
}

impl Phase {
    /// U, V, W for the peaks of wavelength 1 or 2.
    pub fn uvw(&self, wavelength: usize) -> [f64; 3] {
        match (wavelength, self.uvw2) {
            (2, Some(u)) => u,
            _ => [self.u, self.v, self.w],
        }
    }
}

/// The background model.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum BackgroundModel {
    Points(Vec<(f64, f64)>),
    Polynomial { origin: f64, coefficients: Vec<f64> },
    Chebyshev { from: f64, to: f64, coefficients: Vec<f64> },
    /// Points (2θ, intensity) joined by a natural cubic spline (R26, L078; FullProf's cubic splines for Nba < −5,
    /// with the point FullProf adds at the data's last point written in: L081).
    Spline(Vec<(f64, f64)>),
    /// Fourier cosine series Σ B_j cos(j·T/C), C = 180/π for 2θ (manual eq. 3.4′, p. 45).
    Fourier { coefficients: Vec<f64> },
}

impl BackgroundModel {
    pub fn at(&self, t: f64) -> f64 {
        match self {
            BackgroundModel::Points(p) => crate::background::linear_at(p, t),
            BackgroundModel::Spline(p) => crate::background::interpolate(p, crate::background::Interpolation::Spline, &[t])[0],
            // All coefficients 0 is no background, whatever the origin (FullProf's MORE runs with Bkpos = 0, whose
            // background is their profiles alone: runs/background_terms/y2o3_more2, L353).
            BackgroundModel::Polynomial { coefficients, .. } if coefficients.iter().all(|&c| c == 0.0) => 0.0,
            BackgroundModel::Polynomial { origin, coefficients } => {
                let u = t / origin - 1.0;
                coefficients.iter().rev().fold(0.0, |acc, c| acc * u + c)
            }
            BackgroundModel::Fourier { coefficients } => coefficients.iter().enumerate().map(|(j, b)| b * (j as f64 * t.to_radians()).cos()).sum(),
            BackgroundModel::Chebyshev { from, to, coefficients } => {
                let x = (2.0 * (t - from) / (to - from) - 1.0).clamp(-1.0, 1.0);
                let (mut t0, mut t1) = (1.0, x);
                let mut s = 0.0;
                for (j, c) in coefficients.iter().enumerate() {
                    let tj = match j {
                        0 => 1.0,
                        1 => x,
                        _ => {
                            let t2 = 2.0 * x * t1 - t0;
                            t0 = t1;
                            t1 = t2;
                            t2
                        }
                    };
                    s += c * tj;
                }
                s
            }
        }
    }
}

/// The instrumental resolution, added to the sample widths.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum InstrumentalWidths {
    None,
    Function { res: i32, uvwxyz: [f64; 6] },
    Table(Vec<[f64; 3]>),
}

impl InstrumentalWidths {
    /// (H_G², H_L) of the instrument at 2θ (degrees).
    pub fn at(&self, two_theta: f64) -> (f64, f64) {
        let th = (two_theta / 2.0).to_radians();
        match self {
            InstrumentalWidths::None => (0.0, 0.0),
            InstrumentalWidths::Function { res, uvwxyz: [u, v, w, x, y, z] } => match res {
                1 => ((u * th.tan() + v) * th.tan() + w, x * th.tan() + y / th.cos() + z),
                2 => ((u * th.tan() + v) * th.tan() + w, (x * two_theta + y) * two_theta + z),
                _ => ((u * two_theta + v) * two_theta + w, (x * two_theta + y) * two_theta + z),
            },
            InstrumentalWidths::Table(t) => {
                let i = t.iter().position(|r| r[0] >= two_theta).unwrap_or(t.len() - 1).max(1);
                let (a, b) = (t[i - 1], t[i]);
                let f = if b[0] > a[0] { (two_theta - a[0]) / (b[0] - a[0]) } else { 0.0 };
                let hg = a[1] + f * (b[1] - a[1]);
                (hg * hg, a[2] + f * (b[2] - a[2]))
            }
        }
    }
}

/// Whose rule the model follows, rule by rule, where FullProf departs from the papers or its manual (L043, L046,
/// L050). `false` is the default (the papers in `manuals/papers/`, or FullProf's manual where no paper sets the rule),
/// `true` is FullProf's rule.
///
/// | Field | Default | FullProf | Import of a `.pcr` |
/// |---|---|---|---|
/// | `lorentz_polarisation` | eq. 3.26 (p. 53) with K = ½ | twice that (L033) | exact: scale × 2 |
/// | `berar_baldinozzi` | paper's F_b (eq. 13) | F_b = (4u² − 3)F_a (L035) | exact: P1 + 3·P2, P3 + 3·P4 |
/// | `axial_divergence` | S/L, H/L as given (FCJ 1994) | D_L < 0.0002: (H/L + S/L, H/L − S/L) (CrysFML) | exact: `profile::fullprof_axial` |
/// | `preferred_orientation` | exp(−G1 α²) (Dollase 1986, Table 1, no. 5) | exp(+G1 α²) (manual eq. 3.31) | exact: G1 → −G1 |
/// | `excluded_points` | the points within the regions | also a point whose step reaches into one (L033) | exact on the data's grid |
/// | `peak_window` | ±Wdt·H, H the FWHM (manual p. 75) | width without X tan θ, one more low point if X ≠ 0 (L035, L042) | none: kept |
/// | `reflections` | every one whose window reaches the data | from the first point to `reflection_limit`, or 5° past the last (L035) | none: kept |
/// | `microstructure` | β* with each reflection's λ, no strain for U < 0 (microstructure paper, eq. 14) | λ1 and √|U| (L039) | output only |
/// | `window_edge` | the profile unchanged within ±Wdt·H, then tapered to zero by ±1.25·Wdt·H (`profile::taper`, L058) | cut at the edge | none: kept |
/// | `angle_restraints` | weight χ²_ν/σ², σ in degrees (restraints tutorial, p. 2) | π/180 of it (inferred from FullProf's refinements, L058) | none: the default |
/// | `background_edge` | linear background points: constant beyond the first and last | a data point less than one step (FullProf's: `fullprof_step`) below the first lies on the line through the first two (inferred from FullProf's `.prf`, L098) | none: kept |
/// | `displacement` | Sycos·cos θ (manual eq. 3.37, p. 56) | Sycos·cos 2θ (FullProf's `.prf` ticks and `.sub` profiles, L267) | exact when Sycos = 0; otherwise kept |
/// | `background_bounds` | the background function as it is | a polynomial background (with its further terms) 0 where it falls below 0 or above the largest observed count (`Model::background_ceiling`; FullProf's `.prf`, L353) | exact when it stays within them at the data's points; otherwise kept |
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Convention {
    pub lorentz_polarisation: bool,
    pub berar_baldinozzi: bool,
    pub axial_divergence: bool,
    pub preferred_orientation: bool,
    pub excluded_points: bool,
    pub peak_window: bool,
    pub reflections: bool,
    pub microstructure: bool,
    pub angle_restraints: bool,
    pub window_edge: bool,
    #[cfg_attr(feature = "serde", serde(default))]
    pub background_edge: bool,
    #[cfg_attr(feature = "serde", serde(default))]
    pub displacement: bool,
    #[cfg_attr(feature = "serde", serde(default))]
    pub background_bounds: bool,
}

impl Convention {
    /// The papers' rules everywhere.
    pub const DEFAULT: Convention = Convention {
        lorentz_polarisation: false,
        berar_baldinozzi: false,
        axial_divergence: false,
        preferred_orientation: false,
        excluded_points: false,
        peak_window: false,
        reflections: false,
        microstructure: false,
        angle_restraints: false,
        window_edge: false,
        background_edge: false,
        displacement: false,
        background_bounds: false,
    };
    /// FullProf's rules everywhere: a `.pcr` calculated as FullProf calculates it (`Model::from_pcr_fullprof`).
    pub const FULLPROF: Convention = Convention {
        lorentz_polarisation: true,
        berar_baldinozzi: true,
        axial_divergence: true,
        preferred_orientation: true,
        excluded_points: true,
        peak_window: true,
        reflections: true,
        microstructure: true,
        angle_restraints: true,
        window_edge: true,
        background_edge: true,
        displacement: true,
        background_bounds: true,
    };

    /// The names of the rules taken from FullProf.
    pub fn fullprof_rules(&self) -> Vec<&'static str> {
        [
            (self.lorentz_polarisation, "Lorentz–polarisation"),
            (self.berar_baldinozzi, "Bérar–Baldinozzi asymmetry"),
            (self.axial_divergence, "axial divergence (S_L, D_L)"),
            (self.preferred_orientation, "exponential preferred orientation"),
            (self.excluded_points, "excluded points"),
            (self.peak_window, "peak window"),
            (self.reflections, "reflections past the data"),
            (self.microstructure, "microstructure output"),
            (self.angle_restraints, "weight of angle restraints"),
            (self.window_edge, "peak window edge"),
            (self.background_edge, "linear background below its first point"),
            (self.displacement, "sample displacement with cos 2θ"),
            (self.background_bounds, "background 0 outside 0 and the largest count"),
        ]
        .iter()
        .filter(|r| r.0)
        .map(|r| r.1)
        .collect()
    }
}

/// The model of a constant-wavelength pattern.
#[derive(Clone, Debug)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Model {
    /// The papers' rules or FullProf's, rule by rule (`Convention`).
    pub convention: Convention,
    pub neutrons: bool,
    pub lambda: [f64; 2],
    pub ratio: f64,
    pub zero: f64,
    pub sycos: f64,
    pub sysin: f64,
    pub cthm: f64,
    pub rpolarz: f64,
    pub ilo: i32,
    /// μR (cylinder, Debye–Scherrer) or μt (flat plate in transmission, Ilo = 2).
    pub mur: f64,
    pub wdt: f64,
    pub asym_limit: f64,
    /// Nor = 1: modified March function; Nor = 0: exponential.
    pub march: bool,
    pub background: BackgroundModel,
    pub instrumental: InstrumentalWidths,
    pub phases: Vec<Phase>,
    /// Excluded regions (low, high): no peak is added there, as in FullProf, whose calculated profile there is the
    /// background alone.
    pub excluded: Vec<(f64, f64)>,
    /// The weights of the agreement factors (Iwg).
    pub weighting: Weighting,
    /// Distance and angle restraints (P40), used by the engine (`restraints`).
    pub restraints: Vec<crate::restraints::Restraint>,
    /// FullProf's step, as its `.out` prints it ("Step:"), when it is set: its rule for the linear background's low
    /// edge (`Convention::background_edge`) counts in it. For data of x, y, σ (Ins = 10) it is not the data's (D:
    /// 0.100078 for 0.1) and changes between runs of the same data (L098). Unset, the data's first step.
    #[cfg_attr(feature = "serde", serde(default))]
    pub fullprof_step: Option<f64>,
    /// Bragg–Brentano surface micro-absorption P0, Cp, τ (FullProf's eq. 3.36, line 15; P29, L233): all zero, none.
    #[cfg_attr(feature = "serde", serde(default))]
    pub microabsorption: [f64; 3],
    /// A TOF pattern's bank (P05): the pattern's x is then the time of flight (µs), `zero` its zero point, and the
    /// wavelength, Lp, U … Y and the asymmetry are not used.
    #[cfg_attr(feature = "serde", serde(default))]
    pub tof: Option<crate::tof::Bank>,
    /// Further terms of the background, added to `background` (P23, L353): Debye-like terms, inverse powers of 2θ and
    /// profiles read from files (`background_terms`). Empty: none.
    #[cfg_attr(feature = "serde", serde(default))]
    pub background_terms: crate::background_terms::BackgroundTerms,
    /// The largest count of the data FullProf read, above which its background is 0 (`Convention::background_bounds`,
    /// L353); none: no ceiling.
    #[cfg_attr(feature = "serde", serde(default))]
    pub background_ceiling: Option<f64>,
}

/// A peak of the calculated pattern.
#[derive(Clone, Debug, PartialEq)]
pub struct Peak {
    pub phase: usize,
    pub hkl: [i32; 3],
    pub multiplicity: usize,
    /// 1 or 2: the wavelength.
    pub wavelength: usize,
    pub two_theta: f64,
    pub fwhm: f64,
    pub eta: f64,
    /// Integrated intensity including the scale factor.
    pub intensity: f64,
    /// Half-width of the window over which the peak is calculated (degrees): Wdt times `window_width`.
    pub reach: f64,
    /// A TOF peak's rise and decay α, β (µs⁻¹; P05): its `two_theta` is then its time of flight (µs), and `fwhm`, `eta`
    /// the pseudo-Voigt's that the exponentials are convolved with.
    pub tof: Option<[f64; 2]>,
}

/// The calculated pattern at the points asked.
#[derive(Clone, Debug, PartialEq)]
pub struct Calculated {
    pub x: Vec<f64>,
    pub ycalc: Vec<f64>,
    pub background: Vec<f64>,
    pub peaks: Vec<Peak>,
    /// Each phase's profile alone (no background), point by point: `ycalc` is the background plus their sum (the
    /// app's request 7, L070). FullProf writes the same curve, with the background, in its `.sub` files.
    pub per_phase: Vec<Vec<f64>>,
}

impl Model {
    /// The calculated pattern checked for values that are not numbers (TU6): Ok when every Ycal is finite; otherwise
    /// the first point that is not, and the values of the model outside their limits (`refinable::limits`, the
    /// parameter dictionary's), which are where such a pattern comes from (a March–Dollase G₁ = 0 gives NaN).
    pub fn finite(&self, c: &Calculated) -> Result<(), String> {
        let Some(i) = c.ycalc.iter().position(|v| !v.is_finite()) else { return Ok(()) };
        let outside: Vec<String> = self
            .all_fields()
            .into_iter()
            .filter_map(|f| {
                let (v, lim) = (self.get(f), crate::refinable::limits(f, self));
                (crate::refinable::LimitState::of(v, lim) == crate::refinable::LimitState::Outside || !v.is_finite()).then(|| {
                    let b = |x: Option<f64>| x.map_or("none".to_string(), |x| x.to_string());
                    format!("{} = {v} (limits {}, {})", crate::refinable::address(f, self), b(lim.0), b(lim.1))
                })
            })
            .collect();
        Err(format!(
            "the calculated pattern is not a number at point {} (2θ = {}): {}",
            i + 1,
            c.x.get(i).copied().unwrap_or(f64::NAN),
            if outside.is_empty() { "no value of the model is outside its limits".to_string() } else { format!("outside their limits: {}", outside.join("; ")) }
        ))
    }

    /// Lorentz–polarisation factor at 2θ (degrees).
    pub fn lorentz_polarisation(&self, two_theta: f64) -> f64 {
        let th = (two_theta / 2.0).to_radians();
        let c2 = (2.0 * th).cos();
        match self.ilo {
            -1 => 1.0,
            2 => {
                let k = self.rpolarz;
                let (c1, c2m) = (self.cthm.sqrt(), self.cthm);
                let p = if self.neutrons { 1.0 } else { k * (1.0 + c1 * c2 * c2) / (1.0 + c1) + (1.0 - k) * (1.0 + c2m * c2 * c2) / (1.0 + c2m) };
                p / (2.0 * th).sin()
            }
            _ => {
                // X-rays: eq. 3.26 (p. 53) with K = ½ for characteristic radiation, ½(1 + Cthm·cos²2θ). FullProf
                // computes twice that (found from the files: PPH3-sls.prf with Cthm = 0 is reproduced exactly only so,
                // y2o3_structure.pcr with Cthm = 0.78 gives χ² 1.30 against FullProf's 1.20, and 14.4 with eq. 3.26
                // as printed; L033): a factor the scale absorbs. Ilo = 3 (synchrotron, K = Rpolarz): eq. 3.26 as
                // printed, 1 − K + K Cthm cos²2θ, in both conventions, but for Cthm = 0 (no monochromator), which
                // FullProf takes as cos²2θ_M = 1: prbaco_calc_1.prf (BM16, K = 0.1, Cthm = 0) is met so to rms
                // 6·10⁻³ where 1 − K alone leaves every peak at 0.90 of FullProf's (P06, L322). Cthm ≠ 0 with
                // Ilo = 3 is not checked.
                let pol = if self.neutrons {
                    1.0
                } else if self.ilo == 3 {
                    let m = if self.cthm == 0.0 { 1.0 } else { self.cthm };
                    1.0 - self.rpolarz + self.rpolarz * m * c2 * c2
                } else if self.convention.lorentz_polarisation {
                    1.0 + self.cthm * c2 * c2
                } else {
                    0.5 * (1.0 + self.cthm * c2 * c2)
                };
                pol / (2.0 * th.sin().powi(2) * th.cos())
            }
        }
    }

    /// The peaks of every phase between 2θ limits (degrees).
    pub fn peaks(&self, tmin: f64, tmax: f64) -> Vec<Peak> {
        if let Some(bank) = &self.tof {
            return self.tof_peaks(bank, tmin, tmax);
        }
        let mut out = Vec::new();
        for (ip, ph) in self.phases.iter().enumerate() {
            for (iw, &lambda) in self.lambda.iter().enumerate() {
                // λ2 contributes, with weight Ratio, only when it differs from λ1 (Cocus-asym.pcr has Ratio = 1 with
                // λ2 = λ1, and FullProf's intensities there are those of one wavelength).
                let two = self.ratio != 0.0 && (self.lambda[1] - self.lambda[0]).abs() > 1e-9;
                if iw == 1 && !two {
                    continue;
                }
                let weight = if iw == 0 { 1.0 } else { self.ratio };
                let s_max = ((tmax.min(179.0) / 2.0).to_radians().sin() + 0.02) / lambda;
                let dmin = 0.5 / s_max.min(1.0 / lambda);
                let list: Vec<([i32; 3], usize, f64)> = match &ph.intensities {
                    Intensities::Structure(..) => generate(&ph.group, &ph.cell, dmin, 1e9, None).iter().map(|r| (r.hkl, r.multiplicity, f64::NAN)).collect(),
                    Intensities::List(l) => l.clone(),
                };
                for (hkl, mult, given) in list {
                    let h = hkl.map(|v| v as f64);
                    let Some(d) = ph.cell.d_spacing(h) else { continue };
                    if ph.reflection_limit.is_some_and(|l| 1.0 / d > l) {
                        continue;
                    }
                    let sin_th = lambda / (2.0 * d);
                    if sin_th >= 1.0 {
                        continue;
                    }
                    let th = sin_th.asin();
                    let two_theta = 2.0 * th.to_degrees() + self.zero + self.sycos * self.displacement_factor(th) + self.sysin * (2.0 * th).sin();
                    if two_theta < tmin || two_theta > tmax {
                        continue;
                    }
                    // Widths at the Bragg angle, without Zero and the displacement terms: FullProf's `.out` widths (Hw,
                    // HwG, HwL, ETA) are reproduced to 2·10⁻⁵ this way and differ by up to 2.6·10⁻³ at the shifted angle
                    // (YVO3a_calc.out, Zero = 0.0555°, with a large Lorentzian strain; L039).
                    let bragg = 2.0 * th.to_degrees();
                    let (fwhm, eta) = self.peak_width(ph, hkl, bragg, iw + 1);
                    // The default's window reaches past ±Wdt·H, where its taper brings the profile to zero (L058).
                    let reach = if self.convention.peak_window {
                        self.wdt * self.peak_window_width(ph, hkl, bragg, iw + 1)
                    } else if self.convention.window_edge {
                        self.wdt * fwhm
                    } else {
                        self.wdt * fwhm / crate::profile::TAPER_START
                    };
                    let intensity = match &ph.intensities {
                        Intensities::Structure(st, sc) => {
                            let s = 1.0 / (2.0 * d);
                            let f2 = st.intensity(sc, hkl);
                            ph.scale * weight * mult as f64 * self.lorentz_polarisation(2.0 * th.to_degrees()) * f2 * (-2.0 * ph.bov * s * s).exp() * self.preferred_orientation_at(ph, h, th) * self.absorption(th)
                        }
                        // Profile matching (Jbt = 2): the intensity of the file is the total of both wavelengths, without
                        // the scale factor (found from FullProf's files: PPH3-sls_1.prf, Ratio = 1, λ2 = λ1, scale 10⁻⁴).
                        Intensities::List(_) => given * weight / if two { 1.0 + self.ratio } else { 1.0 },
                    };
                    out.push(Peak { phase: ip, hkl, multiplicity: mult, wavelength: iw + 1, two_theta, fwhm, eta, intensity, reach, tof: None });
                }
            }
        }
        out
    }

    /// The peaks of a TOF pattern whose time of flight lies between tmin and tmax (µs; P05, `tof`): each reflection's
    /// position, α, β, and the TCH pseudo-Voigt's FWHM and η from its d-spacing, and its intensity scale·m·F²·d⁴ sin θ_B
    /// ·A·E·exp(−2B_ov s²)·P (`tof` for the evidence). The window is ±Wdt·H (the manual's rule, p. 75), for FullProf's
    /// convention too: FullProf's own TOF window is wider by a few microseconds per Wdt, by no rule found (L341).
    fn tof_peaks(&self, bank: &crate::tof::Bank, tmin: f64, tmax: f64) -> Vec<Peak> {
        use crate::tof::{extinction, gauss_fwhm, Parameters};
        let mut out = Vec::new();
        let th = (bank.two_theta / 2.0).to_radians();
        let (sin_th, s2) = (th.sin(), th.sin().powi(2));
        for (ip, ph) in self.phases.iter().enumerate() {
            let Some(prof) = &ph.tof else { continue };
            let par = Parameters::of(self.zero, bank, prof);
            // The smallest d whose position reaches tmin, by bisection (the position grows with d), less a margin for
            // the windows.
            let (mut lo, mut hi) = (1e-3_f64, 1e3_f64);
            for _ in 0..200 {
                let m = (lo * hi).sqrt();
                if par.at(m).position < tmin {
                    lo = m;
                } else {
                    hi = m;
                }
            }
            let dmin = (lo * 0.9).max(0.05);
            let list: Vec<([i32; 3], usize, f64)> = match &ph.intensities {
                Intensities::Structure(..) => generate(&ph.group, &ph.cell, dmin, 1e9, None).iter().map(|r| (r.hkl, r.multiplicity, f64::NAN)).collect(),
                Intensities::List(l) => l.clone(),
            };
            let volume = ph.cell.volume();
            for (hkl, mult, given) in list {
                let h = hkl.map(|v| v as f64);
                let Some(d) = ph.cell.d_spacing(h) else { continue };
                let shape = par.at(d);
                if shape.position < tmin || shape.position > tmax {
                    continue;
                }
                let (fwhm, eta) = tch(gauss_fwhm(shape.sigma2), shape.lorentz);
                let reach = if self.convention.window_edge { self.wdt * fwhm } else { self.wdt * fwhm / crate::profile::TAPER_START };
                let lambda = 2.0 * d * sin_th;
                let intensity = match &ph.intensities {
                    Intensities::Structure(st, sc) => {
                        let s = 1.0 / (2.0 * d);
                        let f2 = st.intensity(sc, hkl);
                        let a = bank.absorption.at(prof.absorption * lambda, s2);
                        let e = extinction(prof.extinction, f2, lambda, volume, s2);
                        ph.scale * mult as f64 * f2 * d.powi(4) * sin_th * a * e * (-2.0 * ph.bov * s * s).exp() * self.preferred_orientation(ph, h)
                    }
                    Intensities::List(_) => given,
                };
                out.push(Peak { phase: ip, hkl, multiplicity: mult, wavelength: 1, two_theta: shape.position, fwhm, eta, intensity, reach, tof: Some([shape.alpha, shape.beta]) });
            }
        }
        out
    }

    /// The width FullProf's peak window is counted in (`Convention::peak_window`; Wdt of them on each side, manual p. 75,
    /// which says the FWHM, as the default convention takes it): the FWHM, except that
    /// for the TCH profile it leaves out the phase's Lorentzian strain term X·tan θ. Not in the manual; measured on the
    /// FullProf runs (L035): rutile in `rutana_tch_calc.prf`, where the windows of 11 rutile reflections are 0.525° to
    /// 1.244° wide each side and agree with this width to 0.003°, and MgAl2O4 (`MgAl2O4s_Rietveld_calc.prf`, whose only
    /// Lorentzian term is X), where they are Wdt·H_G; the Egyptian make-up and Co-YBCO (X = 0) use the whole FWHM.
    pub fn window_width(&self, ph: &Phase, two_theta: f64) -> f64 {
        self.window_width_of(ph, two_theta, 1)
    }

    /// `window_width` for the peaks of wavelength 1 or 2 (`Phase::uvw2`).
    pub fn window_width_of(&self, ph: &Phase, two_theta: f64, wavelength: usize) -> f64 {
        self.window_width_with(ph, two_theta, wavelength, (0.0, 0.0))
    }

    /// `window_width_of` for reflection hkl: its anisotropic strain's Gaussian part included and its Lorentzian part
    /// left out with X tan θ, as the Lorentzian strain; its anisotropic size included, as the size term Y/cos θ (P25;
    /// checked on prba10k_calc.prf and pdfinal_calc.prf, L267).
    pub fn peak_window_width(&self, ph: &Phase, hkl: [i32; 3], two_theta: f64, wavelength: usize) -> f64 {
        // The Lorentzian strain (ξ D_ST tan θ) left out, as X tan θ is; the anisotropic size kept, as Y/cos θ is
        // (L267: prba10k_calc.prf's windows, with the extra low point below; pdfinal_calc.prf).
        let h = hkl.map(|v| v as f64);
        let strain_g = ph.quartic_strain.as_ref().map_or(0.0, |q| q.widths(h, ph.cell.inv_d_squared(h), (two_theta / 2.0).to_radians().tan()).0);
        let size = ph.shape_size.map_or(0.0, |s| s.width(h, &ph.cell, two_theta));
        self.window_width_with(ph, two_theta, wavelength, (strain_g, size))
    }

    fn window_width_with(&self, ph: &Phase, two_theta: f64, wavelength: usize, extra: (f64, f64)) -> f64 {
        match ph.shape {
            Shape::ThompsonCoxHastings => {
                let th = (two_theta / 2.0).to_radians();
                let (t, c) = (th.tan(), th.cos());
                let (ig2, il) = self.instrumental.at(two_theta);
                let [u, v, w] = ph.uvw(wavelength);
                let hg2 = (u * t + v) * t + w + ph.gauss_size / (c * c) + ig2 + extra.0;
                tch(hg2.max(0.0).sqrt(), ph.y / c + il + extra.1).0
            }
            _ => self.width_with(ph, two_theta, wavelength, extra).0,
        }
    }

    /// The anisotropic strain's (H_G², H_L) for reflection hkl of a phase at 2θ (degrees, the Bragg angle; P25): zero
    /// without a quartic form.
    pub fn strain_widths(&self, ph: &Phase, hkl: [i32; 3], two_theta: f64) -> (f64, f64) {
        let h = hkl.map(|v| v as f64);
        let size = ph.shape_size.map_or(0.0, |s| s.width(h, &ph.cell, two_theta));
        let Some(q) = &ph.quartic_strain else { return (0.0, size) };
        let inv_d2 = ph.cell.inv_d_squared(h);
        let (g, l) = q.widths(h, inv_d2, (two_theta / 2.0).to_radians().tan());
        (g, l + size)
    }

    /// Gaussian and Lorentzian FWHM (H_G, H_L) of the TCH profile at 2θ (degrees; the Bragg angle, without Zero):
    /// the phase's terms (microstructure paper, eq. 19, p. 4) and the resolution file's added.
    pub fn width_components(&self, ph: &Phase, two_theta: f64) -> (f64, f64) {
        self.width_components_with(ph, two_theta, (0.0, 0.0))
    }

    /// `width_components` for reflection hkl, its anisotropic strain included (P25).
    pub fn peak_width_components(&self, ph: &Phase, hkl: [i32; 3], two_theta: f64) -> (f64, f64) {
        self.width_components_with(ph, two_theta, self.strain_widths(ph, hkl, two_theta))
    }

    fn width_components_with(&self, ph: &Phase, two_theta: f64, extra: (f64, f64)) -> (f64, f64) {
        let th = (two_theta / 2.0).to_radians();
        let (t, c) = (th.tan(), th.cos());
        let (ig2, il) = self.instrumental.at(two_theta);
        (((ph.u * t + ph.v) * t + ph.w + ph.gauss_size / (c * c) + ig2 + extra.0).max(0.0).sqrt(), ph.x * t + ph.y / c + il + extra.1)
    }

    /// Gaussian FWHM at 2θ.
    pub fn gauss_width(&self, ph: &Phase, two_theta: f64) -> f64 {
        self.width_components(ph, two_theta).0
    }

    /// FWHM and η at 2θ.
    pub fn width(&self, ph: &Phase, two_theta: f64) -> (f64, f64) {
        self.width_of(ph, two_theta, 1)
    }

    /// `width` for the peaks of wavelength 1 or 2: the second wavelength's U, V, W when the phase has its own
    /// (`Phase::uvw2`, FullProf's Ratio < 0).
    pub fn width_of(&self, ph: &Phase, two_theta: f64, wavelength: usize) -> (f64, f64) {
        self.width_with(ph, two_theta, wavelength, (0.0, 0.0))
    }

    /// FWHM and η of reflection hkl at its Bragg angle 2θ, for wavelength 1 or 2: `width_of` with the phase's
    /// anisotropic strain (P25). The profile's peaks are drawn with it.
    pub fn peak_width(&self, ph: &Phase, hkl: [i32; 3], two_theta: f64, wavelength: usize) -> (f64, f64) {
        self.width_with(ph, two_theta, wavelength, self.strain_widths(ph, hkl, two_theta))
    }

    fn width_with(&self, ph: &Phase, two_theta: f64, wavelength: usize, extra: (f64, f64)) -> (f64, f64) {
        let th = (two_theta / 2.0).to_radians();
        let (t, c) = (th.tan(), th.cos());
        let (ig2, il) = self.instrumental.at(two_theta);
        let [u, v, w] = ph.uvw(wavelength);
        let hg2 = (u * t + v) * t + w + ph.gauss_size / (c * c) + ig2 + extra.0;
        match ph.shape {
            Shape::ThompsonCoxHastings => tch(hg2.max(0.0).sqrt(), ph.x * t + ph.y / c + il + extra.1),
            Shape::PseudoVoigt => (hg2.max(0.0).sqrt(), ph.eta0 + ph.x * two_theta),
            Shape::Gaussian => (hg2.max(0.0).sqrt(), 0.0),
            Shape::Lorentzian => (hg2.max(0.0).sqrt(), 1.0),
        }
    }

    /// The Gaussian H_G² of a phase's peak at 2θ, before it is clamped at 0 (eq. 3.20): negative is a width no peak
    /// can have (FullProf stops on it), which the engine's divergence guard reports (W11).
    pub fn gaussian_h2(&self, ph: &Phase, two_theta: f64, wavelength: usize) -> f64 {
        let th = (two_theta / 2.0).to_radians();
        let (t, c) = (th.tan(), th.cos());
        let (ig2, _) = self.instrumental.at(two_theta);
        let [u, v, w] = ph.uvw(wavelength);
        (u * t + v) * t + w + ph.gauss_size / (c * c) + ig2
    }

    /// `gaussian_h2` for reflection hkl, its anisotropic strain's Gaussian part included (P25).
    pub fn peak_gaussian_h2(&self, ph: &Phase, hkl: [i32; 3], two_theta: f64, wavelength: usize) -> f64 {
        self.gaussian_h2(ph, two_theta, wavelength) + self.strain_widths(ph, hkl, two_theta).0
    }

    /// Absorption transmission factor at θ (eq. 3.33 for a cylinder, eq. 3.34 for a flat plate in transmission), times
    /// the surface micro-absorption (`micro_absorption`).
    pub fn absorption(&self, th: f64) -> f64 {
        self.transmission(th) * self.micro_absorption(th)
    }

    /// Bragg–Brentano surface micro-absorption at θ (FullProf's eqs. 3.35–3.36, manual p. 55; line 15, p. 84):
    /// S_r = 1 − P(h), P(h) = P0 + Cp (τ/sin θ)(1 − τ/sin θ). The manual cites Pitschke et al. and Suortti for its
    /// limits (its refs. [25], [26], not in the repository). 1 when P0, Cp and τ are zero.
    pub fn micro_absorption(&self, th: f64) -> f64 {
        let [p0, cp, tau] = self.microabsorption;
        if p0 == 0.0 && cp == 0.0 {
            return 1.0;
        }
        let r = tau / th.sin();
        1.0 - (p0 + cp * r * (1.0 - r))
    }

    fn transmission(&self, th: f64) -> f64 {
        if self.mur == 0.0 {
            return 1.0;
        }
        let s2 = th.sin().powi(2);
        if self.ilo == 2 {
            (-self.mur / th.cos()).exp() / th.cos()
        } else {
            (-(1.7133 - 0.0368 * s2) * self.mur + (0.0927 + 0.375 * s2) * self.mur * self.mur).exp()
        }
    }

    /// The preferred-orientation factor of reflection h, averaged over the reflections equivalent to h in the phase's
    /// Laue class, as Dollase (1986), J. Appl. Cryst. 19, 267–272, sums P over the m symmetry-equivalent planes (eq. 4,
    /// p. 268): the angle α to the preferred direction differs among them. FullProf's CORR is that average
    /// (Yb2ZrWO8-PO-March, G1 = 0.85 about [1 1 0]: P(2 2 0) 1.010306 in its `.out`, the average of the 12 equivalents
    /// 1.010305, h alone 1.628; L042). March–Dollase (Nor = 1): (r² cos²α + sin²α/r)^(−3/2), Dollase's Table 1, no. 6,
    /// the axial case of March (1932), Z. Kristallogr. 81, 285–297, p. 291 (density ∝ the cube of the deformation
    /// ellipsoid's radius); FullProf's eq. 3.32 with G1 = r. Exponential (Nor = 0): b + (1 − b) exp(−Gα²), Dollase's
    /// no. 5, with G = G1 and b = G2; FullProf's eq. 3.31 writes exp(+G1 α²) (`Convention::preferred_orientation`).
    pub fn preferred_orientation(&self, ph: &Phase, h: [f64; 3]) -> f64 {
        self.preferred_orientation_at(ph, h, std::f64::consts::FRAC_PI_2)
    }

    /// `preferred_orientation` at Bragg angle θ (radians), which an anvil cell's texture needs (P30; θ = 90° is the
    /// flat sample's ψ = 0). A phase's texture along several axes (`texture::MultiAxial`) takes the place of G1, G2:
    /// Σ_j f_j D_j averaged over the equivalent reflections.
    pub fn preferred_orientation_at(&self, ph: &Phase, h: [f64; 3], theta: f64) -> f64 {
        if let Some(t) = &ph.texture {
            let psi = t.psi(theta);
            let cos = texture_cosines(ph, h);
            return cos.iter().map(|c| t.factor(c, psi)).sum::<f64>() / cos.len() as f64;
        }
        let (g1, g2) = (ph.pref[0], ph.pref[1]);
        let nor_march = self.march;
        if (!nor_march && g1 == 0.0) || (nor_march && g1 == 1.0) {
            return 1.0;
        }
        // Angle between h and the preferred direction, in the reciprocal metric.
        let g = ph.cell.reciprocal_metric();
        let dot = |a: [f64; 3], b: [f64; 3]| (0..3).map(|i| (0..3).map(|j| a[i] * g[i][j] * b[j]).sum::<f64>()).sum::<f64>();
        let p = ph.pref_direction;
        let one = |h: [f64; 3]| {
            let cosa = (dot(h, p) / (dot(h, h) * dot(p, p)).sqrt()).clamp(-1.0, 1.0);
            let alpha = cosa.abs().acos();
            if nor_march {
                g2 + (1.0 - g2) * ((g1 * cosa).powi(2) + alpha.sin().powi(2) / g1).powf(-1.5)
            } else {
                // Dollase (1986), Table 1, no. 5 (p. 268): b + (1 − b) exp(−Gα²); FullProf's eq. 3.31 has exp(+G1 α²).
                let g = if self.convention.preferred_orientation { -g1 } else { g1 };
                g2 + (1.0 - g2) * (-g * alpha * alpha).exp()
            }
        };
        let eq = crate::reflections::equivalents(&ph.group, h.map(|v| v.round() as i32));
        eq.iter().map(|k| one(k.map(|v| v as f64))).sum::<f64>() / eq.len() as f64
    }

    /// Under FullProf's rule for linear background points (`Convention::background_edge`), the first two points
    /// (x₀, y₀), (x₁, y₁) when the data point t lies less than one step below x₀: there FullProf reads the background on
    /// the line through them, where elsewhere below x₀ (and above the last point) it holds the end point's value. The
    /// step is FullProf's (`fullprof_step`), else the data's first; a point exactly one step below is not on the line.
    /// Inferred from FullProf's `.prf` (L098): B's 12.60° (x₀ 12.626°) and E's 15.4667° (x₀ 15.4931°) on the line;
    /// A's 11.15° (x₀ 11.2°, step 0.05°) and LaMnO3's 4.40° (x₀ 4.45°, step 0.05°) constant; D's 0.79° (x₀ 0.89°)
    /// on the line, FullProf's step there being 0.100078°.
    /// With a file's b1–b4 (`background_terms::FileScale`) the test is made where the points are read, (1 + b3)·t + b4,
    /// the step stretched with it; the segment is then the points' own, to be read there (`background_x`).
    pub fn background_edge_segment(&self, x: &[f64], t: f64) -> Option<[(f64, f64); 2]> {
        let BackgroundModel::Points(p) = &self.background else { return None };
        if !self.convention.background_edge || p.len() < 2 || x.len() < 2 {
            return None;
        }
        let stretch = self.background_terms.file.as_ref().map_or(1.0, |f| 1.0 + f.stretch);
        let step = self.fullprof_step.unwrap_or((x[1] - x[0]).abs()) * stretch;
        let u = self.background_x(t);
        (u < p[0].0 && p[0].0 - u < step * (1.0 - 1e-6)).then_some([p[0], p[1]])
    }

    /// The background's function alone at the data's x, before a file's scale and offset: its points (or series) read
    /// at `background_x`, on FullProf's line below the first point where its rule applies.
    pub fn background_function_at(&self, x: &[f64], t: f64) -> f64 {
        let u = self.background_x(t);
        match self.background_edge_segment(x, t) {
            Some([(x0, y0), (x1, y1)]) => y0 + (y1 - y0) * (u - x0) / (x1 - x0),
            None => self.background.at(u),
        }
    }

    /// The derivative of the background at t with respect to a file's b_i (i = 0 … 3: scale, offset, stretch, shift;
    /// eq. 3.6). The stretch and shift move where the points are read: the function's slope there (the points joined
    /// linearly, its slope a segment's; at a point itself, the mean of the two), times the scale, times t for the
    /// stretch.
    pub fn background_file_derivative(&self, i: usize, x: &[f64], t: f64) -> f64 {
        let Some(f) = &self.background_terms.file else { return 0.0 };
        match i {
            0 => self.background_function_at(x, t),
            1 => 1.0,
            _ => {
                let u = self.background_x(t);
                let h = 1e-6 * u.abs().max(1.0);
                let at = |v: f64| match self.background_edge_segment(x, t) {
                    Some([(x0, y0), (x1, y1)]) => y0 + (y1 - y0) * (v - x0) / (x1 - x0),
                    None => self.background.at(v),
                };
                let slope = (at(u + h) - at(u - h)) / (2.0 * h);
                f.scale * slope * if i == 2 { t } else { 1.0 }
            }
        }
    }

    /// The background at the points x: its function and its further terms (`background_terms`), under FullProf's rule
    /// 0 where that falls below 0 or above `background_ceiling` (`Convention::background_bounds`).
    pub fn background_at(&self, x: &[f64]) -> Vec<f64> {
        let mut b = self.background_unbounded(x);
        if self.background_bounded() {
            for v in &mut b {
                if self.background_out_of_bounds(*v) {
                    *v = 0.0;
                }
            }
        }
        b
    }

    /// Whether FullProf's rule sets a background value b to 0 (`Convention::background_bounds`): below 0, or above the
    /// largest observed count. FullProf's runs (L353): nac_3t2_debye.prf prints 0 at the 1 504 points where eq. 3.4's
    /// sum is negative, and its Ycal there is the Bragg part alone; y2o3_poly12_b9.prf prints 0 from 144.95°, where B9's
    /// term first passes Y2O3's largest count, 11 248 (11 249.22 there, 11 240.66 at 144.94°). Where it applies:
    /// `background_bounded`.
    pub fn background_out_of_bounds(&self, b: f64) -> bool {
        b < 0.0 || self.background_ceiling.is_some_and(|c| b > c)
    }

    /// Whether FullProf's bounds apply: its rule, on a polynomial background (with its further terms). Not on points:
    /// FullProf's background goes below 0 between and beyond its points (LBF_lmopr.prf, −0.37 to −1.6 from 124.75°;
    /// runs/modes/D's d_calc.prf, −167.6 at 0.79° on the line below the first point, L353). Its Chebyshev and Fourier
    /// series are not known to be bounded: none of FullProf's runs takes them out of bounds.
    pub fn background_bounded(&self) -> bool {
        self.convention.background_bounds && matches!(self.background, BackgroundModel::Polynomial { .. })
    }

    /// The background at the points x without FullProf's bounds.
    pub fn background_unbounded(&self, x: &[f64]) -> Vec<f64> {
        let terms = !self.background_terms.is_empty();
        let file = self.background_terms.file.as_ref();
        x.iter()
            .map(|&t| {
                let b = self.background_function_at(x, t);
                // A file's b1 and b2 (eq. 3.6; b3, b4 in `background_x`).
                let b = file.map_or(b, |f| f.scale * b + f.offset);
                if terms {
                    b + self.background_terms.at(t, self.background_q(t))
                } else {
                    b
                }
            })
            .collect()
    }

    /// The Q (Å⁻¹) of the Debye-like background terms at x: 4π sin θ/λ₁ with 2θ = x, as FullProf (manual eq. 3.4, no
    /// zero shift: `nac_3t2_calc.prf`, L353); 2π/d at a TOF pattern's x (believed: no FullProf reference).
    pub fn background_q(&self, x: f64) -> f64 {
        match &self.tof {
            Some(_) => self.d_spacing_at(x).map_or(0.0, |d| 2.0 * std::f64::consts::PI / d),
            None => 4.0 * std::f64::consts::PI * (x / 2.0).to_radians().sin() / self.lambda[0],
        }
    }

    /// The background split into its parts at the points x (P49.2): the function (with the inverse powers of 2θ,
    /// which complete FullProf's 12-term polynomial), each Debye-like term, each profile from a file.
    pub fn background_parts(&self, x: &[f64]) -> Vec<(String, Vec<f64>)> {
        let t = &self.background_terms;
        let mut smooth = self.background_at(x);
        let mut parts = Vec::new();
        for j in 0..t.debye.len() {
            parts.push((format!("debye_{}", j + 1), x.iter().map(|&v| t.debye_at(j, self.background_q(v))).collect::<Vec<f64>>()));
        }
        for j in 0..t.profiles.len() {
            parts.push((format!("file_{}", j + 1), x.iter().map(|&v| t.profile_at(j, v)).collect::<Vec<f64>>()));
        }
        for (_, p) in &parts {
            for (s, v) in smooth.iter_mut().zip(p) {
                *s -= v;
            }
        }
        parts.insert(0, ("smooth".to_string(), smooth));
        parts
    }

    /// The abscissa the background is read at: the TOF scaled to 180/T_max for the Fourier series, whose argument is
    /// jπT/T_max (`tof::Bank::tmax`); x itself otherwise.
    /// The d-spacing (Å) at x of the pattern (P46): for a TOF pattern the inverse of T(d) with the bank's refined
    /// constants (`tof::Bank::d_spacing`); for a constant-wavelength one λ₁/(2 sin θ) with 2θ = x − Zero (the sample's
    /// displacement and transparency left out). None where there is none.
    pub fn d_spacing_at(&self, x: f64) -> Option<f64> {
        match &self.tof {
            Some(b) => b.d_spacing(self.zero, x),
            None => {
                let s = ((x - self.zero) / 2.0).to_radians().sin();
                (s > 0.0 && self.lambda[0] > 0.0).then(|| self.lambda[0] / (2.0 * s))
            }
        }
    }

    /// With a file's b1–b4 (FullProf's Nba = 1, eq. 3.6), (1 + b3)·t + b4, where its points are read.
    pub fn background_x(&self, t: f64) -> f64 {
        match (&self.tof, &self.background) {
            (Some(b), BackgroundModel::Fourier { .. }) if b.tmax > 0.0 => t * 180.0 / b.tmax,
            _ => self.background_terms.file.as_ref().map_or(t, |f| f.at(t)),
        }
    }

    /// What Sycos multiplies at the Bragg angle θ (radians): cos θ (the manual's eq. 3.37, p. 56), or FullProf's cos 2θ
    /// (`Convention::displacement`, L267).
    pub fn displacement_factor(&self, th: f64) -> f64 {
        if self.convention.displacement {
            (2.0 * th).cos()
        } else {
            th.cos()
        }
    }

    /// The calculated pattern at the points x (2θ, degrees).
    pub fn calculate(&self, x: &[f64]) -> Calculated {
        let background = self.background_at(x);
        let mut ycalc = background.clone();
        let mut per_phase = vec![vec![0.0; x.len()]; self.phases.len()];
        let peaks = self.select_peaks(x);
        let excluded = self.excluded_points(x);
        for pk in &peaks {
            let ph = &mut per_phase[pk.phase];
            self.peak_profile(x, &excluded, pk, |i, v| {
                ycalc[i] += v;
                ph[i] += v;
            });
        }
        Calculated { x: x.to_vec(), ycalc, background, peaks, per_phase }
    }

    /// The peaks `calculate` uses for the points x.
    pub fn select_peaks(&self, x: &[f64]) -> Vec<Peak> {
        let (lo, hi) = (x.iter().copied().fold(f64::MAX, f64::min), x.iter().copied().fold(f64::MIN, f64::max));
        if self.tof.is_some() {
            // Every TOF reflection whose window reaches the data (FullProf's .prf ticks run below the data too), within
            // half the first time of flight and twice the last: the reflections are generated only as far beyond the
            // data as the widest window there reaches (the largest reach over 32 d-spacings on each side), which spares
            // the generation of the many short d-spacings no window brings in.
            let (a, b) = self.tof_margins(lo, hi);
            let mut p = self.peaks((lo - a).max(lo * 0.5), (hi + b).min(hi * 2.0));
            p.retain(|k| k.two_theta + k.reach >= lo && k.two_theta - k.reach <= hi);
            return p;
        }
        if self.convention.reflections {
            // FullProf generates no reflection below the first point (its `.out` prints S(1/d)min at Thmin), and past
            // the last one as far as `reflection_limit` says, or 5° when it is unset (L035).
            let top = if self.phases.iter().any(|p| p.reflection_limit.is_some()) { 180.0 } else { hi + 5.0 };
            self.peaks(lo, top)
        } else {
            // Every reflection whose window reaches the data, on either side.
            let mut p = self.peaks(0.0, 180.0);
            p.retain(|k| k.two_theta + k.reach >= lo && k.two_theta - k.reach <= hi);
            p
        }
    }

    /// How far below `lo` and above `hi` (µs) a TOF reflection's window may reach the data from: the largest reach of
    /// any phase's peak over 32 d-spacings from the d of lo/2 to that of lo, and from hi to 2hi.
    fn tof_margins(&self, lo: f64, hi: f64) -> (f64, f64) {
        let Some(bank) = &self.tof else { return (0.0, 0.0) };
        let mut m = (0.0f64, 0.0f64);
        for ph in &self.phases {
            let Some(prof) = &ph.tof else { continue };
            let par = crate::tof::Parameters::of(self.zero, bank, prof);
            let reach = |d: f64| {
                let s = par.at(d);
                let h = tch(crate::tof::gauss_fwhm(s.sigma2), s.lorentz).0;
                if self.convention.window_edge { self.wdt * h } else { self.wdt * h / crate::profile::TAPER_START }
            };
            let widest = |t0: f64, t1: f64| match (bank.d_spacing(self.zero, t0), bank.d_spacing(self.zero, t1)) {
                (Some(d0), Some(d1)) => (0..32).map(|k| reach(d0 + (d1 - d0) * k as f64 / 31.0)).fold(0.0, f64::max),
                // Without d-spacings there, the whole of the old range.
                _ => f64::INFINITY,
            };
            m = (m.0.max(widest(lo * 0.5, lo)), m.1.max(widest(hi, hi * 2.0)));
        }
        m
    }

    /// The points left out of the calculation: the default excludes the points within the regions; FullProf also
    /// those a step away (`fullprof_excluded`).
    pub fn excluded_points(&self, x: &[f64]) -> Vec<bool> {
        // A TOF pattern: FullProf leaves out a point whose step forward reaches into a region, but not one whose step
        // back does (it calculates CeO2_PEARL_calc.prf's 2000.34 µs after a region ending at 2000, and alum_calc.prf's
        // first point, 24001.2, after one ending at 24000, and leaves out its last, 41697.15, before one from 41700).
        if self.convention.excluded_points && self.tof.is_some() {
            return (0..x.len())
                .map(|i| {
                    let step = if i + 1 < x.len() { x[i + 1] - x[i] } else if i > 0 { x[i] - x[i - 1] } else { 0.0 };
                    let tol = 1e-6 * step.abs();
                    self.excluded.iter().any(|&(a, b)| x[i] <= b && x[i] + step > a + tol)
                })
                .collect();
        }
        if self.convention.excluded_points {
            self.fullprof_excluded(x)
        } else {
            x.iter().map(|&t| self.excluded.iter().any(|&(a, b)| t >= a && t <= b)).collect()
        }
    }

    /// A constant-wavelength peak's Bragg angle 2θ (degrees) from its reflection's d and its wavelength: its position
    /// without the zero and the sample's displacement and transparency.
    pub fn bragg_two_theta(&self, pk: &Peak) -> f64 {
        let ph = &self.phases[pk.phase];
        let h = pk.hkl.map(|v| v as f64);
        let g = ph.cell.reciprocal_metric();
        let inv_d2: f64 = (0..3).map(|i| (0..3).map(|j| h[i] * g[i][j] * h[j]).sum::<f64>()).sum();
        let lambda = self.lambda[pk.wavelength.clamp(1, 2) - 1];
        2.0 * (lambda * inv_d2.sqrt() / 2.0).clamp(-1.0, 1.0).asin().to_degrees()
    }

    /// One peak's contribution I·Ω at each point of its window, handed to `sink(point, value)`.
    pub fn peak_profile(&self, x: &[f64], excluded: &[bool], pk: &Peak, mut sink: impl FnMut(usize, f64)) {
        let conv = self.convention;
        if let Some([alpha, beta]) = pk.tof {
            let (low, high) = (pk.two_theta - pk.reach, pk.two_theta + pk.reach);
            for i in x.partition_point(|&t| t < low)..x.len() {
                if x[i] > high {
                    break;
                }
                if excluded[i] {
                    continue;
                }
                let dx = x[i] - pk.two_theta;
                let omega = crate::tof::back_to_back(dx, alpha, beta, pk.fwhm, pk.eta).0;
                sink(i, pk.intensity * if conv.window_edge { omega } else { omega * taper(dx, pk.reach) });
            }
            return;
        }
        {
            let ph = &self.phases[pk.phase];
            // The asymmetry at the reflection's Bragg angle, the zero and the sample's shifts left out: below AsyLim, and
            // Bérar & Baldinozzi's 1/tan θ, 1/tan 2θ there (P41's pdd runs, Zero −0.20°: with the shifted position the
            // lowest peak, at 10.6°, is off by 1.6·10⁻³ of the strongest, with the Bragg angle by 10⁻⁵; L367).
            let bragg = self.bragg_two_theta(pk);
            let asym = bragg < self.asym_limit;
            // Axial divergence for TCH (Npr = 7) and the pseudo-Voigt (Npr = 12).
            let fcj = (asym && matches!(ph.shape, Shape::ThompsonCoxHastings | Shape::PseudoVoigt) && (ph.s_l != 0.0 || ph.d_l != 0.0)).then(|| {
                let (s_l, h_l) = if conv.axial_divergence { fullprof_axial(ph.s_l, ph.d_l) } else { (ph.s_l, ph.d_l) };
                Fcj::new(pk.eta, pk.fwhm, s_l, h_l, pk.two_theta)
            });
            let berar_too = asym && ph.asym.iter().any(|&a| a != 0.0);
            let berar = berar_too && fcj.is_none();
            // The window is ±reach about 2θ_h, also for the axial-divergence profile. With `Convention::peak_window`: (rutana_tch_calc.prf: the windows
            // of the rutile and anatase reflections, which have FCJ asymmetry, are symmetric about 2θ_h; L035). Its
            // points run from the last point at or below 2θ_h − reach to the last at or below 2θ_h + reach, so the
            // low side takes one point beyond the window: CrysFML's `Locate` (V(j) ≤ x < V(j+1)) as its profile range
            // is set (deps-pycrysfml 2258e77, CFML/Testing/nFP/Src/nFP_Reflections.f90 lines 205–208;
            // CFML/Src/CFML_Maths/Math_Locate.f90 lines 64–68), and as FullProf's windows show (L042).
            // The extra low point is found in the windows of TCH phases with a Lorentzian strain, X ≠ 0 (LaB6rx_calc.prf:
            // window-edge points 50 → 2) or ξ ≠ 0 (prba10k_calc.prf, L267), whose windows leave that strain out of the
            // width; in those without (PPH3, citric acid, Y2O3, Yb2ZrWO8) it is not, and their profiles agree to 10⁻⁵
            // without it. Why FullProf's range differs between them is not known (L042).
            let (low, high) = (pk.two_theta - pk.reach, pk.two_theta + pk.reach);
            // A Lorentzian strain of the quartic form (ξ ≠ 0) takes the extra low point as X ≠ 0 does (prba10k_calc.prf:
            // X = 0, ξ = 0.148, its windows end one point lower than ±Wdt H; L267).
            let lorentzian_strain = ph.x != 0.0 || ph.quartic_strain.as_ref().is_some_and(|q| q.lorentzian != 0.0);
            let extra = conv.peak_window && ph.shape == Shape::ThompsonCoxHastings && lorentzian_strain;
            let start = if extra { x.partition_point(|&t| t <= low).saturating_sub(1) } else { x.partition_point(|&t| t < low) };
            for i in start..x.len() {
                let dx = x[i] - pk.two_theta;
                if x[i] > high {
                    break;
                }
                if excluded[i] {
                    continue;
                }
                let omega = match &fcj {
                    // With both, the Bérar–Baldinozzi factor multiplies the axial-divergence profile, as FullProf does
                    // (window experiment yb_irf1_x0, Asy1 = 0.00901 with the file's S_L: 10⁻⁵ with the product,
                    // 2.9·10⁻² with FCJ alone; L042).
                    Some(f) if berar_too => f.at(x[i]) * berar_baldinozzi(dx / pk.fwhm, bragg, ph.asym, conv.berar_baldinozzi),
                    Some(f) => f.at(x[i]),
                    None => {
                        let base = match ph.shape {
                            Shape::Gaussian => gaussian(dx, pk.fwhm),
                            Shape::Lorentzian => lorentzian(dx, pk.fwhm),
                            _ => pseudo_voigt(dx, pk.eta, pk.fwhm),
                        };
                        if berar {
                            base * berar_baldinozzi(dx / pk.fwhm, bragg, ph.asym, conv.berar_baldinozzi)
                        } else {
                            base
                        }
                    }
                };
                let omega = if conv.window_edge { omega } else { omega * taper(dx, pk.reach) };
                sink(i, pk.intensity * omega);
            }
        }
    }
}

/// A model and its observed pattern, read from a FullProf `.pcr` file and the files it names.
#[derive(Clone, Debug)]
pub struct Loaded {
    pub model: Model,
    pub data: Option<Pattern>,
    pub data_path: Option<PathBuf>,
}

/// A phase's anisotropic size from its `.pcr` (P25, L267): LorSiz with Size-Model 1 (platelets), −1 (needles), both
/// with the vector of LINE 36, or −14 (the antiphase domains of FullProf's Pd3MnD run). The others are refused with
/// their reason: the user's selective rules (−2 to −9) and the particular ones (2 to 13; P26), the spherical harmonics
/// (15 to 22) and the quadratic form (Str = 2, 3): no FullProf run with its output in the repository uses them.
fn size_of(ph: &crate::pcr::PcrPhase) -> Result<Option<crate::microstructure::ShapeSize>, String> {
    use crate::microstructure::{ShapeSize, SizeShape};
    if ph.flag("Sol") != 0 {
        return Err(format!("hkl-dependent peak shifts (Sol = {}, LINE 34–35) are P26, not calculated (L267)", ph.flag("Sol")));
    }
    let shape = match ph.size_model {
        0 if ph.lorentz_size == 0.0 => return Ok(None),
        0 => return Err(format!("LorSiz {} without a size model (Size-Model 0)", ph.lorentz_size)),
        1 | -1 => {
            let t = ph.size_vector.ok_or("Size-Model ±1 without its platelet–needle vector (LINE 36)")?;
            if t == [0.0; 3] {
                return Err("Size-Model ±1 with the vector 0 0 0 (LINE 36)".into());
            }
            if ph.size_model == 1 { SizeShape::Platelet(t) } else { SizeShape::Needle(t) }
        }
        -14 => SizeShape::Antiphase,
        m @ -9..=-2 => return Err(format!("selective size broadening by the user's rules (Size-Model {m}) is P26, not calculated (L267)")),
        m @ 15..=22 => return Err(format!("anisotropic size by spherical harmonics (Size-Model {m}): no FullProf run with its output in the repository uses it (L267)")),
        m => return Err(format!("Size-Model {m}: a selective rule of the manual's list (pp. 23–24) not calculated; only −14, the antiphase rule of Pd3MnD, is (P26, L267)")),
    };
    Ok(Some(ShapeSize { shape, value: ph.lorentz_size }))
}

/// A phase's quartic form of strain from its `.pcr` (P25, L211, L267): |Str| = 1 with Strain-Model 2, −2, 3, 4 or 5
/// (the Laue classes 1 2/m 1, 1 1 2/m, mmm, 4/m and 4/mmm), with its Lorentzian share ξ. The rest is refused with the
/// reason it is not calculated:
/// - FullProf's strain models with Str = 0 (manual pp. 28–30) but model 8 (`hexagonal_strain`): described in words,
///   not by equations (NdNi-LTOb's model 6: "fluctuations along a, b and c with corr(a, b) = −1, corr(a, c) = 1,
///   corr(b, c) = 0");
/// - the generalised strain of Str = 1 with Strain-Model 0, or Str = 3 (STR4–STR15, ref. 17 of the manual, not in the
///   repository): no run;
/// - FullProf's form for −1 (the manual lists S_301 twice, p. 119) and its trigonal, hexagonal and cubic models (6–14):
///   the core builds those forms from the rotations (`quartic_terms`), but no FullProf run checks FullProf's reading.
fn quartic_of(ph: &crate::pcr::PcrPhase, group: &SpaceGroup) -> Result<Option<crate::microstructure::QuarticStrain>, String> {
    use crate::microstructure::{fullprof_strain_model, quartic_terms, QuarticStrain};
    let st = ph.flag("Str");
    if st == 0 && ph.strain_model == 8 {
        return hexagonal_strain(ph, group).map(Some);
    }
    if st == 0 && ph.strain_model != 0 {
        return Err(format!("FullProf's strain model {} (Str = 0) is described in words, not by an equation (manual pp. 28–30): not calculated (L211)", ph.strain_model));
    }
    if st == 0 {
        if ph.strain.iter().any(|&v| v != 0.0) {
            return Err("Str1–Str3 without a strain model".into());
        }
        return Ok(None);
    }
    if st.abs() != 1 || ph.strain_model == 0 {
        return Err(format!("Str = {st} with Strain-Model {} (the generalised strain or size of the manual's ref. 17, not in the repository) is not calculated (L211)", ph.strain_model));
    }
    let Some(order) = fullprof_strain_model(ph.strain_model) else {
        return Err(if ph.strain_model == 1 {
            "the quartic form of strain for −1 (Strain-Model 1) is not read: the manual lists S_301 twice and no S_310 (p. 119), so the order of its 15 terms is not known (L211)".to_string()
        } else {
            format!("the quartic form of strain with Strain-Model {}: no FullProf run checks this Laue class's terms; Strain-Models 2, −2, 3, 4 and 5 are read (L211, L267)", ph.strain_model)
        });
    };
    let allowed = quartic_terms(group)?;
    // FullProf's terms are the class's, but for 4/m's S_310, which FullProf has not (taken as 0; L267).
    let extra: Vec<&[u8; 3]> = allowed.iter().filter(|t| !order.contains(t)).collect();
    if order.iter().any(|t| !allowed.contains(t)) || extra.iter().any(|t| **t != [3, 1, 0]) {
        return Err(format!("Strain-Model {} is not the Laue class of the phase's group ({})", ph.strain_model, group.laue));
    }
    if ph.quartic.len() != order.len() {
        return Err(format!("{} S_HKL read for Strain-Model {}, {} expected", ph.quartic.len(), ph.strain_model, order.len()));
    }
    // The terms in the group's order (`quartic_terms`), FullProf's values put in their places.
    let terms = allowed.iter().map(|t| (*t, order.iter().position(|o| o == t).map_or(0.0, |i| ph.quartic[i]))).collect();
    Ok(Some(QuarticStrain::for_group(group, terms, ph.xi.map_or(0.0, |x| x.0))))
}

/// FullProf's strain model 8, "general anisotropic strain of hexagonal symmetry" (Str = 0; manual p. 29): with
/// 1/d² = A m + C l², m = h² + k² + hk, Var(1/d²) = Saa² m² + Scc² l⁴ + 2 Saa Scc Cac m l² and FWHM² = Var (P/(1/d²))²,
/// P = (180/π) √(8 ln 2) 10⁻³ tan θ, Str1–Str3 = Saa, Scc, Cac. It is Stephens' quartic form of the trigonal and
/// hexagonal classes, whose invariants are m² (S_400), l⁴ (S_004) and m l² (S_112), with S_400 = Saa², S_004 = Scc²,
/// S_112 = 2 Saa Scc Cac, Gaussian (ξ = 0); the other terms of −3 and −3m are 0. FullProf's own run decides the factor:
/// with the manual's 10⁻³ in P (S_400 = 100 Saa², σ² in 10⁻⁶) scor_rxn_calc_1.prf is missed by up to 290%, with 10⁻⁴
/// (σ² in 10⁻⁸, as the quartic form, p. 22) it is met to rms 7·10⁻⁵ (FullProf's joint run scor_rxn, P06, L322).
fn hexagonal_strain(ph: &crate::pcr::PcrPhase, group: &SpaceGroup) -> Result<crate::microstructure::QuarticStrain, String> {
    if !matches!(group.laue, "-3" | "-3m" | "6/m" | "6/mmm") {
        return Err(format!("FullProf's strain model 8 is hexagonal; the phase's Laue class is {}", group.laue));
    }
    let terms = crate::microstructure::quartic_terms(group)?;
    let [saa, scc, cac] = ph.strain;
    let value = |t: &[u8; 3]| match t {
        [4, 0, 0] => saa * saa,
        [0, 0, 4] => scc * scc,
        [1, 1, 2] => 2.0 * saa * scc * cac,
        _ => 0.0,
    };
    let terms = terms.iter().map(|t| (*t, value(t))).collect();
    Ok(crate::microstructure::QuarticStrain::for_group(group, terms, ph.xi.map_or(0.0, |x| x.0)))
}

/// A phase's TOF profile from its `.pcr` (P05, lines 26–33 of a TOF pattern): what FullProf's TOF runs check, the rest
/// refused with its reason.
fn tof_profile_of(ph: &crate::pcr::PcrPhase, npr: i32) -> Result<crate::tof::Profile, String> {
    let t = ph.tof.as_ref().ok_or("a TOF pattern's phase without its TOF profile lines")?;
    for (v, name) in [(t.gauss_strain, "Iso-GStrain"), (t.gauss_size, "Iso-GSize"), (t.lorentz_size_ani, "Ani-LSize"), (t.lorentz_strain, "Iso-LorStrain"), (t.lorentz_size, "Iso-LorSize")] {
        if v != 0.0 {
            return Err(format!("{name} = {v}: TOF size and strain terms are not calculated (no FullProf TOF run in the repository uses them; L341)"));
        }
    }
    if ph.flag("Str") != 0 || ph.strain_model != 0 || ph.size_model != 0 || ph.strain.iter().any(|&s| s != 0.0) {
        return Err("anisotropic strain or size in a TOF pattern is not calculated (no FullProf TOF run uses it; L341)".into());
    }
    if t.absorption[1] != 0.0 {
        return Err(format!("ABSCOR2 = {}: the second TOF absorption parameter is not calculated (FullProf's two runs with it stopped, L334)", t.absorption[1]));
    }
    let thermal = match (npr, t.thermal) {
        (10, Some(th)) => Some(th),
        (10, None) => return Err("Npr = 10 without line 31 (alph0t beta0t alph1t beta1t)".into()),
        _ => None,
    };
    let [s2, s1, s0] = t.sigma;
    let [g2, g1, g0] = t.gamma;
    Ok(crate::tof::Profile {
        sigma: [s0, s1, s2],
        sigma_q: t.sigma_q,
        gamma: [g0, g1, g2],
        alpha: [t.alpha[0], t.alpha[1], t.alpha_q],
        beta: [t.beta[0], t.beta[1], t.beta_q],
        thermal,
        extinction: ph.shape1,
        absorption: t.absorption[0],
    })
}

/// A data file read as FullProf reads it for a TOF pattern (P05): the registry's reader for `ins`, and for a GSAS file
/// (Ins = 12) the operations that turn GSAS-II's reading into FullProf's, returned so that a control file keeps them:
/// a CONST bank's x is the TOF in µs (GSAS-II divides it by 100 and multiplies back once the instrument says TOF), its
/// ALT counts as written (×1000: arg_si_calc.prf's Yobs), and a RALF bank's counts per 100 µs of its bin, y·100/Δt
/// (CeO2_PEARL_calc.prf's Yobs, to 10⁻⁴).
pub fn read_tof_data(text: &str, ins: i32) -> Result<(Pattern, Vec<crate::data_ops::Operation>), String> {
    use crate::data_ops::{Normalization, Operation};
    let f = crate::formats::by_ins(ins).ok_or(format!("data format Ins={ins} is not in FullProf's list"))?;
    let probe = crate::formats::Probe::new("data", text.as_bytes());
    let d = (f.read)(&probe).map_err(|e| e.to_string())?;
    let mut scan = d.scans.into_iter().next().ok_or(format!("Ins={ins}: no pattern in the file"))?;
    let mut ops = Vec::new();
    if f.id == "gsas" {
        let ralf = scan.axis == crate::formats::Axis::TimeOfFlight;
        crate::tof::time_axis(&mut scan);
        ops.push(Operation::Normalize(if ralf { Normalization::Step(100.0) } else { Normalization::Factor(1000.0) }));
    }
    let scan = crate::data_ops::apply_all(&scan, &ops, &|_: &str, _: Option<&str>| Err("no file".to_string()))?;
    Ok((scan.pattern, ops))
}

/// A file named in a `.pcr`, looked for as FullProf users lay out tutorials: in the folder of the `.pcr`, its parent,
/// and the parent's other folders (e.g. `pcr-files/` next to `Data-files/`), then the grandparent and its folders. A
/// Windows path keeps only its last part; a name found in no place as written is looked for in any case.
/// cos α between each reflection equivalent to h and each axis of the phase's texture, in the reciprocal metric:
/// [equivalent][axis] (P30). Empty without a texture.
pub fn texture_cosines(ph: &Phase, h: [f64; 3]) -> Vec<Vec<f64>> {
    let Some(t) = &ph.texture else { return Vec::new() };
    let g = ph.cell.reciprocal_metric();
    let dot = |a: [f64; 3], b: [f64; 3]| (0..3).map(|i| (0..3).map(|j| a[i] * g[i][j] * b[j]).sum::<f64>()).sum::<f64>();
    crate::reflections::equivalents(&ph.group, h.map(|v| v.round() as i32))
        .iter()
        .map(|k| {
            let k = k.map(|v| v as f64);
            t.axes.iter().map(|a| (dot(k, a.direction) / (dot(k, k) * dot(a.direction, a.direction)).sqrt()).clamp(-1.0, 1.0)).collect()
        })
        .collect()
}

pub fn find_file(dir: &Path, name: &str) -> Option<PathBuf> {
    let base = name.rsplit(['\\', '/']).next().unwrap_or(name);
    // The name in the folder itself first; then the path as named, from the folder ("../data/x.dat", as an export
    // names it), before the wider search: without it a fresh build of `background_terms` found its data only where
    // another test had left a copy (L398).
    if name.contains('/') && !name.contains('\\') && !dir.join(base).is_file() {
        let p = dir.join(name);
        if p.is_file() {
            return Some(p);
        }
    }
    let mut places = vec![dir.to_path_buf()];
    // The folder, its parent and the parent's folders, then the grandparent and its folders (the FullProf runs of
    // runs/window_experiments/ name the data of runs/Yb2ZrWO8/ and runs/YVO3/, L041).
    let mut level = dir;
    for _ in 0..2 {
        let Some(parent) = level.parent() else { break };
        places.push(parent.to_path_buf());
        if let Ok(rd) = std::fs::read_dir(parent) {
            let mut subs: Vec<PathBuf> = rd.filter_map(|e| e.ok().map(|e| e.path())).filter(|p| p.is_dir() && p != level).collect();
            subs.sort();
            places.extend(subs);
        }
        level = parent;
    }
    // FullProf runs on file systems that ignore case (macOS, Windows): case A of the mode course names
    // "srsclrrt.dat" for SrScLRRT.dat (L090). The exact name first, then any case.
    places.iter().map(|p| p.join(base)).find(|p| p.is_file()).or_else(|| {
        places.iter().find_map(|p| std::fs::read_dir(p).ok()?.filter_map(|e| e.ok()).map(|e| e.path()).find(|f| f.is_file() && f.file_name().is_some_and(|n| n.to_string_lossy().eq_ignore_ascii_case(base))))
    })
}

fn read_latin1(path: &Path) -> Result<String, String> {
    let b = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    Ok(String::from_utf8(b.clone()).unwrap_or_else(|_| b.iter().map(|&c| c as char).collect()))
}

/// The atoms of a phase, with symmetry modes applied (x = x0 + Σ A·coefficient·V).
pub fn atoms_of(ph: &PcrPhase) -> Vec<Atom> {
    let mut atoms: Vec<Atom> = ph.atoms.iter().map(|a| Atom { label: a.label.clone(), species: a.species.clone(), x: a.x, biso: a.biso, occ: a.occ, beta: a.beta, charge: None }).collect();
    if let Some(m) = &ph.modes {
        for (mode, label, _, v, coef) in &m.vectors {
            let Some(&(_, amp)) = m.amplitudes.get(mode - 1) else { continue };
            if let Some(a) = atoms.iter_mut().find(|a| &a.label == label) {
                for (xk, vk) in a.x.iter_mut().zip(v) {
                    *xk += amp * coef * vk;
                }
            }
        }
    }
    atoms
}

impl Model {
    /// Import a `.pcr` file into the default model (L047, L050): its values, refined under FullProf's rules, are
    /// translated wherever a rule differs and an exact translation exists (`translate_from_fullprof`); the rules with
    /// none (the peak window, the reflections past the data) stay FullProf's, and `convention` says which. `dir` is the
    /// folder the file is in (for the data, resolution and hkl files); `stem` its name without extension (CODFIL).
    pub fn from_pcr(p: &Pcr, dir: &Path, stem: &str) -> Result<Loaded, String> {
        Model::from_pcr_with_data(p, dir, stem, None)
    }

    /// As `from_pcr`, with the data file given (when the `.pcr` names one that is not at hand).
    pub fn from_pcr_with_data(p: &Pcr, dir: &Path, stem: &str, data_file: Option<&Path>) -> Result<Loaded, String> {
        let mut l = Model::from_pcr_fullprof_with_data(p, dir, stem, data_file)?;
        let x = l.data.as_ref().map(|d| d.x.clone());
        l.model.translate_from_fullprof(x.as_deref());
        Ok(l)
    }

    /// A `.pcr` file as FullProf calculates it: its values as written, under `Convention::FULLPROF` (the FullProf
    /// tests, `drs calc --fullprof`).
    pub fn from_pcr_fullprof(p: &Pcr, dir: &Path, stem: &str) -> Result<Loaded, String> {
        Model::from_pcr_fullprof_with_data(p, dir, stem, None)
    }

    /// Translate a model's values from FullProf's rules into the default ones, rule by rule, where the translation is
    /// exact; the rules without one keep FullProf's (L050):
    /// - Lorentz–polarisation: FullProf's X-ray Lp (Ilo = 0, 1, 4) is twice eq. 3.26, so the scale of each structure
    ///   phase is doubled (profile-matching intensities carry no Lp);
    /// - Bérar–Baldinozzi: FullProf's F_b is the paper's + 3·F_a, so P1 + 3·P2 → P1 and P3 + 3·P4 → P3;
    /// - axial divergence: FullProf's (S_L, D_L) become the paper's (S/L, H/L) (`profile::fullprof_axial`);
    /// - exponential preferred orientation: FullProf's exp(G1 α²) is Dollase's exp(−Gα²) with G = −G1;
    /// - excluded regions: given the data's points, each run of points FullProf leaves out becomes a region from its
    ///   first to its last point, which the default rule leaves out exactly; without the points the rule stays;
    /// - the linear background below its first point keeps FullProf's rule (`Convention::background_edge`): a
    ///   translation would add background points and renumber the refined ones;
    /// - microstructure: an output convention, set to the default (no effect on the pattern);
    /// - angle restraints: their σ is the user's statement, so the default weighs it as stated (FullProf's pull is π/180
    ///   of that, L058);
    /// - sample displacement: FullProf's Sycos·cos 2θ has no exact translation into Sycos·cos θ, so the rule stays where
    ///   Sycos ≠ 0 (L267).
    pub fn translate_from_fullprof(&mut self, x: Option<&[f64]>) {
        let c = self.convention;
        if c.lorentz_polarisation {
            if !self.neutrons && !matches!(self.ilo, -1 | 2 | 3) {
                for ph in self.phases.iter_mut().filter(|ph| matches!(ph.intensities, Intensities::Structure(..))) {
                    ph.scale *= 2.0;
                }
            }
            self.convention.lorentz_polarisation = false;
        }
        if c.berar_baldinozzi {
            for ph in &mut self.phases {
                ph.asym[0] += 3.0 * ph.asym[1];
                ph.asym[2] += 3.0 * ph.asym[3];
            }
            self.convention.berar_baldinozzi = false;
        }
        if c.axial_divergence {
            for ph in self.phases.iter_mut().filter(|ph| ph.s_l != 0.0 || ph.d_l != 0.0) {
                (ph.s_l, ph.d_l) = fullprof_axial(ph.s_l, ph.d_l);
            }
            self.convention.axial_divergence = false;
        }
        if c.preferred_orientation {
            if !self.march {
                for ph in &mut self.phases {
                    ph.pref[0] = -ph.pref[0];
                }
            }
            self.convention.preferred_orientation = false;
        }
        if let (true, Some(x), None) = (c.excluded_points, x, &self.tof) {
            let mask = self.fullprof_excluded(x);
            let mut regions = Vec::new();
            let mut i = 0;
            while i < x.len() {
                if mask[i] {
                    let first = i;
                    while i + 1 < x.len() && mask[i + 1] {
                        i += 1;
                    }
                    regions.push((x[first], x[i]));
                }
                i += 1;
            }
            self.excluded = regions;
            self.convention.excluded_points = false;
        }
        self.convention.microstructure = false;
        self.convention.angle_restraints = false;
        if self.sycos == 0.0 {
            self.convention.displacement = false;
        }
        // FullProf's bounds of the background change nothing where it stays within them at the data's points.
        if let (true, Some(x)) = (c.background_bounds, x) {
            if !self.background_bounded() || !self.background_unbounded(x).iter().any(|&b| self.background_out_of_bounds(b)) {
                self.convention.background_bounds = false;
                self.background_ceiling = None;
            }
        }
    }

    /// The points FullProf leaves out: a point is excluded when a step from it, to the next point or from the previous
    /// one, reaches into an excluded region (its `.prf`: citric-acid-sls_1.prf leaves out 38.958° for a region starting
    /// at 38.96°, and Cocus-asym_calc.prf, step 0.05°, gives the background alone at 5.019° after a region ending at
    /// 5.00°; L033).
    pub fn fullprof_excluded(&self, x: &[f64]) -> Vec<bool> {
        (0..x.len())
            .map(|i| {
                let step = if i + 1 < x.len() { x[i + 1] - x[i] } else if i > 0 { x[i] - x[i - 1] } else { 0.0 };
                let back = if i > 0 { x[i] - x[i - 1] } else { step };
                // A step that ends on the region's edge does not reach into it, to the data's rounding:
                // prbaco_calc_2.prf calculates 147.95° below a region from 148.00°, step 0.05 (L322).
                let tol = 1e-6 * step.abs().max(back.abs());
                self.excluded.iter().any(|&(a, b)| x[i] - back < b - tol && x[i] + step > a + tol)
            })
            .collect()
    }

    /// As `from_pcr_fullprof`, with the data file given.
    pub fn from_pcr_fullprof_with_data(p: &Pcr, dir: &Path, stem: &str, data_file: Option<&Path>) -> Result<Loaded, String> {
        if p.flag("Cry") != 0 {
            return Err("integrated-intensity jobs (Cry ≠ 0) have no profile to calculate".into());
        }
        let job = p.flag("Job");
        let tof = p.flag("Uni") == 1;
        let neutrons = matches!(job, 1 | 3 | -1 | -3);
        if tof && !matches!(job, -1 | -3) {
            return Err(format!("Job = {job} with Uni = 1: a TOF pattern is Job −1 (neutrons) or −3 (calculated)"));
        }
        if !tof && !matches!(job, 0..=3) {
            return Err(format!("Job = {job} is not a constant-wavelength X-ray or neutron pattern"));
        }
        let mut phases = Vec::new();
        for (k, ph) in p.phases.iter().enumerate() {
            let jbt = ph.flag("Jbt");
            let npr = if ph.flag("Npr") != 0 { ph.flag("Npr") } else { p.flag("Npr") };
            let shape = Shape::from_npr(npr)?;
            if tof != matches!(npr, 9 | 10 | 13) {
                return Err(format!("phase {}: Npr = {npr} in a {} pattern", k + 1, if tof { "TOF" } else { "constant-wavelength" }));
            }
            let tof_profile = if tof { Some(tof_profile_of(ph, npr).map_err(|e| format!("phase {}: {e}", k + 1))?) } else { None };
            let group = SpaceGroup::from_fullprof_symbol(&ph.symbol).map_err(|e| format!("phase {}: {e}", k + 1))?;
            let c = ph.cell;
            let cell = Cell::new(c[0], c[1], c[2], c[3], c[4], c[5]).map_err(|e| format!("phase {}: {e}", k + 1))?;
            let quartic_strain = if tof { None } else { quartic_of(ph, &group).map_err(|e| format!("phase {}: {e}", k + 1))? };
            let shape_size = if tof { None } else { size_of(ph).map_err(|e| format!("phase {}: {e}", k + 1))? };
            let mut modes = None;
            let mut rigid = Vec::new();
            let intensities = match (jbt, ph.flag("Irf")) {
                (0 | 4 | 6, 0) => {
                    if let (6, Some(m)) = (jbt, &ph.modes) {
                        modes = Some(crate::modes::Modes::from_pcr(m, &ph.atoms).map_err(|e| format!("phase {}: {e}", k + 1))?);
                    }
                    let mut st = Structure { group: group.clone(), cell, atoms: atoms_of(ph) };
                    // Rigid groups (P41): the atoms where the groups place them.
                    if jbt.abs() == 4 {
                        // In degrees whatever the file's DEG, so that each value has one unit.
                        rigid = crate::rigid_body::groups(ph).map_err(|e| format!("phase {}: {e}", k + 1))?.into_iter().map(crate::rigid_body::Group::in_degrees).collect();
                        let mut x: Vec<[f64; 3]> = st.atoms.iter().map(|a| a.x).collect();
                        crate::rigid_body::place(&mut rigid, &cell, &mut x);
                        for (a, x) in st.atoms.iter_mut().zip(&x) {
                            a.x = *x;
                        }
                        // TLS (P6 = 5.0x with Furth ≠ 0; P41, L367): T, L, S and the centre (`rigid_body::read_tls`).
                        let species: Vec<String> = st.atoms.iter().map(|a| a.species.clone()).collect();
                        crate::rigid_body::read_tls(ph, &mut rigid, &species, &x, &cell).map_err(|e| format!("phase {}: {e}", k + 1))?;
                        for (a, b) in crate::rigid_body::betas(&rigid, &cell, &x) {
                            st.atoms[a].beta = Some(b);
                        }
                    }
                    let radiation = if neutrons { Radiation::Neutron } else { Radiation::Xray { lambda: p.lambda1, anomalous: true } };
                    let user = user_scatterers(p);
                    let sc = st.atoms.iter().map(|a| crate::scattering::resolve(&a.species, radiation, &user)).collect::<Result<Vec<_>, _>>()?;
                    Intensities::Structure(Box::new(st), sc)
                }
                // Profile matching: the intensities of CODFILn.hkl (read with Irf = 2; with Irf = 0 FullProf iterates them
                // and writes the same file at the end, which is read if present).
                // Without that file (Irf = 0), every reflection the group allows up to 5° past Thmax, each with the same
                // intensity, 1: the starting point of a Le Bail extraction (`Model::le_bail`).
                (2, irf @ (0 | 2)) => match find_file(dir, &format!("{stem}{}.hkl", k + 1)) {
                    Some(path) => {
                        let text = read_latin1(&path)?;
                        let mut list = Vec::new();
                        for l in text.lines().skip(2) {
                            let v: Vec<f64> = l.split_whitespace().map_while(|s| s.parse().ok()).collect();
                            if v.len() >= 5 {
                                list.push(([v[0] as i32, v[1] as i32, v[2] as i32], v[3] as usize, v[4]));
                            }
                        }
                        Intensities::List(list)
                    }
                    None if irf == 0 => {
                        let lambda = p.lambda1.min(if p.lambda2 > 0.0 { p.lambda2 } else { p.lambda1 });
                        let top = (p.thmax + 5.0).min(179.0);
                        let dmin = lambda / (2.0 * (top / 2.0).to_radians().sin());
                        Intensities::List(generate(&group, &cell, dmin, 1e9, None).iter().map(|r| (r.hkl, r.multiplicity, 1.0)).collect())
                    }
                    None => return Err(format!("phase {}: profile matching with Irf = 2 needs the intensities of {stem}{}.hkl", k + 1, k + 1)),
                },
                (j, i) => return Err(format!("phase {}: Jbt = {j}, Irf = {i} is not calculated yet (Rietveld: Jbt 0, 4, 6 with Irf 0; profile matching: Jbt 2 with Irf 2)", k + 1)),
            };
            phases.push(Phase {
                name: ph.name.clone(),
                group,
                cell,
                intensities,
                scale: ph.scale,
                bov: ph.bov,
                shape,
                // A TOF phase's Shape1 is its extinction (`tof::Profile`).
                eta0: if ph.tof.is_some() { 0.0 } else { ph.shape1 },
                u: ph.u,
                v: ph.v,
                w: ph.w,
                x: ph.x,
                y: ph.y,
                gauss_size: ph.gauss_size,
                asym: ph.asym,
                // Npr = 12: S_L, D_L of line 33 ("Additional asymmetry parameters (S_L, D_L)", manual p. 112).
                s_l: if npr == 12 { ph.shp[0] } else { ph.s_l },
                d_l: if npr == 12 { ph.shp[1] } else { ph.d_l },
                pref: ph.pref,
                pref_direction: ph.pref_direction,
                reflection_limit: None,
                extract: jbt == 2,
                atz: ph.atz,
                brindley: ph.flags.get("Brind").copied().filter(|&b| b > 0.0).unwrap_or(1.0),
                internal_standard: None,
                valence_sum_range: None,
                modes,
                rigid,
                // FullProf's Ratio < 0: the second wavelength has its own U, V, W (manual p. 75, line 32).
                uvw2: if p.ratio < 0.0 { ph.uvw2 } else { None },
                parent: None,
                quartic_strain,
                shape_size,
                tof: tof_profile,
                texture: ph.texture.clone(),
            });
        }
        // The data file: the name FullProf noted in its comment, else CODFIL.dat.
        let data_path = match data_file {
            Some(f) => Some(f.to_path_buf()),
            None => {
                let name = p.data_file.clone().unwrap_or_else(|| stem.to_string());
                [name.clone(), format!("{name}.dat"), format!("{stem}.dat")].iter().find_map(|n| find_file(dir, n))
            }
        };
        let data = match &data_path {
            Some(path) if tof => Some(read_tof_data(&read_latin1(path)?, p.flag("Ins")).map_err(|e| format!("{}: {e}", path.display()))?.0),
            Some(path) => Some(crate::data::read(&read_latin1(path)?, p.flag("Ins")).map_err(|e| format!("{}: {e}", path.display()))?),
            None => None,
        };
        // Ste > 1: every Ste-th point kept, from the first (manual p. 64; prba10k_calc.prf, L267).
        let data = match (data, p.flag("Ste")) {
            (Some(d), ste) if ste > 1 => Some(crate::data_ops::rebin_factor(&d, ste as usize, crate::data_ops::RebinMode::Keep)?),
            (d, _) => d,
        };
        let (mut debye, mut inverse, mut bac_file) = (Vec::new(), Vec::new(), None);
        let background = match &p.background {
            // Nba < 0 written as points: cubic splines (manual p. 69: Nba < 0 with |Nba| > 4; −4 and −5 are read as
            // the Fourier and Chebyshev series). FullProf's spline is the natural one through the points and one more
            // of its own, at the last data point with the last point's value (L081: runs/background, Nba = −6 and −17,
            // within the .prf's print of 0.1 only so): the import writes that point, so that the default's natural
            // spline through the points is FullProf's curve. At the data's last point, not Thmax (spl6thmax), and not
            // when the last point lies past the data (spl6beyond): L087, eleven variants within FullProf's print.
            Background::Points(pts) if p.flag("Nba") < 0 => {
                let mut q = pts.clone();
                let end = data.as_ref().and_then(|d| d.x.last().copied()).unwrap_or(p.thmax);
                if let Some(&(xl, yl)) = q.last() {
                    if end > xl {
                        q.push((end, yl));
                    }
                }
                BackgroundModel::Spline(q)
            }
            Background::Points(pts) => BackgroundModel::Points(pts.clone()),
            Background::Polynomial(c) => BackgroundModel::Polynomial { origin: p.bkpos, coefficients: c.clone() },
            Background::Chebyshev(c) => {
                // A TOF pattern's over the .pcr's TOF-min to TOF-max (CeO2_PEARL_calc.out: "TOFmin: 1500.25 TOFmax:
                // 19380.78", the data reaching 19400.19; its .prf's background is met so to the print, to 2.7·10⁻³ over
                // the data's range); a constant-wavelength pattern's over the data's range.
                let (from, to) = match &data {
                    _ if tof && p.thmax > p.thmin => (p.thmin, p.thmax),
                    Some(d) if !d.is_empty() => (d.x[0], d.x[d.len() - 1]),
                    _ => (p.thmin, p.thmax),
                };
                BackgroundModel::Chebyshev { from, to, coefficients: c.clone() }
            }
            Background::Fourier(c) => BackgroundModel::Fourier { coefficients: c.clone() },
            // Nba = 1: the background of CODFIL.bac, b1·B[(1 + b3)2θ + b4] + b2 (manual eq. 3.6, p. 46; b1 = 0 read as
            // 1), its points joined linearly (believed: the manual does not say how it interpolates the file). The
            // file's points are kept as read and b1–b4 as four values of the background (`background_terms::FileScale`,
            // the app's request 47): B is read at (1 + b3)2θ + b4 (`runs/background_terms/y2o3_bac_b3b4.prf`, b3 = 0.01,
            // b4 = 0.5: within 0.06; y2o3_bac_b1b2, P23, L353).
            Background::Other(1, c) => {
                let b = |k: usize| c.get(k).copied().unwrap_or(0.0);
                let b1 = if b(0) == 0.0 { 1.0 } else { b(0) };
                if 1.0 + b(2) <= 0.0 {
                    return Err(format!("Nba = 1: b3 = {} stretches the .bac background to nothing (1 + b3 ≤ 0)", b(2)));
                }
                let path = find_file(dir, &format!("{stem}.bac")).ok_or(format!("background file {stem}.bac not found (Nba = 1)"))?;
                let pf = crate::background::read_path(&path)?;
                bac_file = Some(crate::background_terms::FileScale { file: path.file_name().map_or(String::new(), |n| n.to_string_lossy().into_owned()), scale: b1, offset: b(1), stretch: b(2), shift: b(3) });
                BackgroundModel::Points(pf.points)
            }
            // Nba = −1: the polynomial of six coefficients and Debye-like terms, Σ A_j sin(Q r_j)/(Q r_j) (manual eq.
            // 3.4, p. 45): the reader gives the six, then the amplitudes, then the distances (P23, L353).
            Background::Other(-1, c) if c.len() >= 6 && (c.len() - 6) % 2 == 0 => {
                let n = (c.len() - 6) / 2;
                debye = (0..n).map(|j| crate::background_terms::DebyeTerm { amplitude: c[6 + j], distance: c[6 + n + j] }).collect();
                BackgroundModel::Polynomial { origin: p.bkpos, coefficients: c[..6].to_vec() }
            }
            // Nba = −3: nine coefficients of the polynomial, then three of the inverse powers of 2θ (P23, L353).
            Background::Other(-3, c) => {
                inverse = c.iter().skip(9).copied().collect();
                BackgroundModel::Polynomial { origin: p.bkpos, coefficients: c.iter().take(9).copied().collect() }
            }
            // Nba = −2: the Fourier filtering of the background (manual eq. 3.7) is not reproduced by the runs of P23
            // (`runs/background_terms/y2o3_filter`, L353): refused.
            crate::pcr::Background::Other(n, _) => {
                let what = match n {
                    -1 => "the polynomial with Debye-like terms, written with too few values,",
                    -2 => "the Fourier filtering of the background",
                    _ => "this background",
                };
                return Err(format!("Nba = {n}: {what} is not calculated (P23)"));
            }
        };
        // MORE = 2, 3: profiles read from {data}_n.bac, each times its coefficient (P23, L353).
        let mut profiles = Vec::new();
        if !p.profiles.is_empty() {
            let data_stem = data_path.as_ref().and_then(|d| d.file_stem()).map_or(stem.to_string(), |s| s.to_string_lossy().to_string());
            for (k, &c) in p.profiles.iter().enumerate() {
                let name = format!("{data_stem}_{}.bac", k + 1);
                let path = find_file(dir, &name).ok_or(format!("background profile {name} not found (MORE = 2, 3 on line 14)"))?;
                let pf = crate::background::read_path(&path)?;
                profiles.push(crate::background_terms::Profile { file: path.to_string_lossy().to_string(), coefficient: c, points: pf.points });
            }
        }
        let background_terms = crate::background_terms::BackgroundTerms { debye, inverse, profiles, file: bac_file };
        // Ratio's sign is FullProf's flag for the second wavelength's own U, V, W (`Phase::uvw2`); the intensity ratio
        // is its size.
        let (mut lambda, mut ratio, mut cthm) = ([p.lambda1, p.lambda2], p.ratio.abs(), p.cthm);
        // A TOF bank (P05): line 14, line 16 (Npr = 10), Iabscor, and its resolution file's values (Res = 5, 6).
        let bank = if tof {
            let npr = p.flag("Npr");
            let thermal = match (npr, p.thermal) {
                (10, Some([zero, dtt1, dtt2, cross, width])) => Some(crate::tof::Thermal { zero, dtt1, dtt2, cross, width }),
                (10, None) => return Err("Npr = 10 without line 16 (Zerot Dtt1t Dtt2t x-cross Width)".into()),
                _ => None,
            };
            // Iabscor 0 (or an unknown value) with every ABSCOR1 = 0: no absorption, whichever function.
            let none = phases.iter().all(|ph| ph.tof.as_ref().is_none_or(|t| t.absorption == 0.0));
            let absorption = match crate::tof::Absorption::from_iabscor(p.iabscor) {
                Ok(a) => a,
                Err(_) if none => crate::tof::Absorption::Cylinder,
                Err(e) => return Err(e),
            };
            let mut two_theta = p.two_theta_bank;
            let instrument = match &p.resolution_file {
                None => None,
                Some(f) => {
                    let path = find_file(dir, f).ok_or(format!("resolution file {f} not found"))?;
                    let (i, twoth) = crate::tof::Instrument::read(&read_latin1(&path)?).map_err(|e| format!("{}: {e}", path.display()))?;
                    // TWOTH replaces the .pcr's 2ThetaBank (docs/TOF_FullProf.pdf).
                    if let Some(t) = twoth {
                        two_theta = t;
                    }
                    Some(i)
                }
            };
            Some(crate::tof::Bank { dtt1: p.dtt[0], dtt2: p.dtt[1], dtt_1overd: p.dtt[2], two_theta, thermal, absorption, instrument, tmax: p.thmax })
        } else {
            None
        };
        let instrumental = match &p.resolution_file {
            _ if tof => InstrumentalWidths::None,
            None => InstrumentalWidths::None,
            Some(f) => {
                let path = find_file(dir, f).ok_or(format!("resolution file {f} not found"))?;
                match Resolution::read(&read_latin1(&path)?, p.flag("Res"))? {
                    Resolution::Function { res, uvwxyz, keywords } => {
                        // The file's wavelengths, ratio and monochromator replace those of the `.pcr` (window
                        // experiment yb_irf1_x0, L042: its `.pcr` gives 1.5406, 1.5406, 1, Cthm 0.7899; FullProf's `.out`
                        // prints the 1.54056, 1.54440, 0.4374 and 0.8000 of LaB6a.irf's WAVE and CTHM lines).
                        if let Some(w) = keywords.get("WAVE").filter(|w| w.len() >= 3) {
                            (lambda, ratio) = ([w[0], w[1]], w[2]);
                        }
                        if let Some(c) = keywords.get("CTHM").and_then(|c| c.first()) {
                            cthm = *c;
                        }
                        // The axial divergence of the instrument (ASYM S_L D_L) is added to the phase's (window
                        // experiment yb_all_x0, L042: phase S_L = 0.0637 and the file's 0.0637 agree with FullProf's
                        // profile as 0.1274, to 10⁻⁵; with either alone they differ by up to 53%).
                        if let Some(a) = keywords.get("ASYM").filter(|a| a.len() >= 2) {
                            for ph in phases.iter_mut() {
                                ph.s_l += a[0];
                                ph.d_l += a[1];
                            }
                        }
                        InstrumentalWidths::Function { res, uvwxyz: uvwxyz[0] }
                    }
                    Resolution::Table(t) => InstrumentalWidths::Table(t),
                }
            }
        };
        let model = Model {
            convention: Convention::FULLPROF,
            neutrons,
            lambda,
            ratio,
            zero: p.zero,
            sycos: p.sycos,
            sysin: p.sysin,
            cthm,
            rpolarz: p.rpolarz,
            ilo: p.flag("Ilo"),
            mur: p.mur,
            wdt: p.wdt,
            asym_limit: p.asym_limit,
            march: p.flag("Nor") == 1,
            background,
            instrumental,
            phases,
            excluded: p.excluded.clone(),
            weighting: Weighting::from_iwg(p.flag("Iwg")),
            restraints: crate::restraints::from_pcr(p)?,
            fullprof_step: None,
            microabsorption: p.microabsorption.unwrap_or_default(),
            tof: bank,
            background_terms,
            // FullProf's ceiling of the background: the largest count of the data it read (L353).
            background_ceiling: data.as_ref().map(|d| d.y.iter().copied().fold(f64::NEG_INFINITY, f64::max)).filter(|v| v.is_finite()),
        };
        Ok(Loaded { model, data, data_path })
    }
}

/// The user's scattering factors of a `.pcr`'s line 12 (P11; manual pp. 70, 79–80), by FullProf's rules: for neutrons
/// (Job 1, 3) ITY = 0 gives b in DFP, ITY = 1 a magnetic form factor (left out: the core calculates no magnetic phase);
/// for X-rays DFP and DFPP are f′ and f″, with f(s) from the nine coefficients (ITY = 0, Nsc > 0), from a table
/// (ITY = 0, Nsc < 0), or from the tables (ITY = 2).
pub fn user_scatterers(p: &Pcr) -> Vec<crate::scattering::UserScatterer> {
    let neutrons = matches!(p.flag("Job"), 1 | 3 | -1 | -3);
    p.user_scattering
        .iter()
        .filter(|u| !(neutrons && u.ity == 1))
        .map(|u| {
            let mut s = crate::scattering::UserScatterer { species: u.name.clone(), ..Default::default() };
            if neutrons {
                // DFPP is b's imaginary part: FullProf uses it (nabd4_calc.out, L251), though its manual says it is ignored.
                (s.scattering_length, s.scattering_length_imaginary) = (Some(u.dfp), u.dfpp);
            } else {
                (s.f_prime, s.f_double_prime) = (u.dfp, u.dfpp);
                if u.ity == 0 && u.coefficients.len() == 9 {
                    let k = &u.coefficients;
                    s.form_factor = Some([k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8]]);
                }
                if u.ity == 0 {
                    s.form_factor_table = u.table.clone();
                }
            }
            s
        })
        .collect()
}

/// Agreement between an observed and a calculated pattern over the points not excluded (manual eqs. 3.9–3.12, p. 43):
/// (Rp, Rwp, Rexp, χ²) with weights 1/σ², n the number of points used and `p` refined parameters.
pub fn agreement(obs: &Pattern, ycalc: &[f64], excluded: &[(f64, f64)], p: usize) -> (f64, f64, f64, f64) {
    agreement_weighted(obs, ycalc, excluded, p, Weighting::Variance)
}

/// The weights of the observations (Iwg of the `.pcr`, manual LINE 4 and p. 43).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Weighting {
    /// Iwg = 0: w = 1/σ², the variances of the data.
    Variance,
    /// Iwg = 1: maximum likelihood, w = 1/σ²_calc with Poisson variances from the calculated profile, scaled by the
    /// data's σ²/y where the data are normalised (p. 43). Reproduces FullProf's χ² of LaB6rx_calc.out (12.94; with
    /// w = 1/σ² the core gives 15.35; L039).
    MaximumLikelihood,
    /// Iwg = 2: unit weights.
    Unit,
}

impl Weighting {
    pub fn from_iwg(iwg: i32) -> Weighting {
        match iwg {
            1 => Weighting::MaximumLikelihood,
            2 => Weighting::Unit,
            _ => Weighting::Variance,
        }
    }
}

/// Agreement factors (Rp, Rwp, Rexp, χ²) with the given weights (manual eqs. 3.9–3.12, p. 43).
pub fn agreement_weighted(obs: &Pattern, ycalc: &[f64], excluded: &[(f64, f64)], p: usize, weighting: Weighting) -> (f64, f64, f64, f64) {
    let (mut sa, mut sy, mut swd, mut swy, mut n) = (0.0, 0.0, 0.0, 0.0, 0usize);
    for (i, &x) in obs.x.iter().enumerate() {
        if excluded.iter().any(|&(a, b)| x >= a && x <= b) {
            continue;
        }
        let (y, s2) = (obs.y[i], obs.sigma[i] * obs.sigma[i]);
        let var = match weighting {
            Weighting::Variance => s2,
            Weighting::MaximumLikelihood => ycalc[i] * if y > 0.0 && s2 > 0.0 { s2 / y } else { 1.0 },
            Weighting::Unit => 1.0,
        };
        if var <= 0.0 {
            continue;
        }
        let w = 1.0 / var;
        let d = y - ycalc[i];
        sa += d.abs();
        sy += y.abs();
        swd += w * d * d;
        swy += w * y * y;
        n += 1;
    }
    let dof = (n.saturating_sub(p)).max(1) as f64;
    (100.0 * sa / sy, 100.0 * (swd / swy).sqrt(), 100.0 * (dof / swy).sqrt(), swd / dof)
}

/// A pseudo-random generator for simulations (SplitMix64, Steele, Lea & Flood 2014), reproducible from its seed.
pub struct Random(u64);

impl Random {
    pub fn new(seed: u64) -> Random {
        Random(seed)
    }
    /// Uniform in (0, 1).
    pub fn uniform(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        ((z >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    }
    /// A Poisson deviate of mean λ: multiplication of uniforms (Knuth) below 30, the transformed rejection with
    /// squeeze of Hörmann (Insurance: Mathematics and Economics 12, 39, 1993; "PTRS") above.
    pub fn poisson(&mut self, lambda: f64) -> f64 {
        // A mean that is not a finite positive number draws 0: the rejection loop below never accepts a NaN (TU6).
        if !lambda.is_finite() || lambda <= 0.0 {
            return 0.0;
        }
        if lambda < 30.0 {
            let l = (-lambda).exp();
            let (mut k, mut p) = (0.0, 1.0);
            loop {
                p *= self.uniform();
                if p <= l {
                    return k;
                }
                k += 1.0;
            }
        }
        let slam = lambda.sqrt();
        let loglam = lambda.ln();
        let b = 0.931 + 2.53 * slam;
        let a = -0.059 + 0.02483 * b;
        let invalpha = 1.1239 + 1.1328 / (b - 3.4);
        let vr = 0.9277 - 3.6224 / (b - 2.0);
        loop {
            let u = self.uniform() - 0.5;
            let v = self.uniform();
            let us = 0.5 - u.abs();
            let k = ((2.0 * a / us + b) * u + lambda + 0.43).floor();
            if us >= 0.07 && v <= vr {
                return k;
            }
            if k < 0.0 || (us < 0.013 && v > us) {
                continue;
            }
            if (v * invalpha / (a / (us * us) + b)).ln() <= -lambda + k * loglam - ln_gamma(k + 1.0) {
                return k;
            }
        }
    }
}

/// ln Γ(x) for x > 0 (Lanczos approximation, g = 7, n = 9; Numerical Recipes 3rd ed. §6.1 coefficients).
fn ln_gamma(x: f64) -> f64 {
    const C: [f64; 9] = [0.999_999_999_999_809_9, 676.520_368_121_885_1, -1_259.139_216_722_402_8, 771.323_428_777_653_1, -176.615_029_162_140_6, 12.507_343_278_686_905, -0.138_571_095_265_720_12, 9.984_369_578_019_572e-6, 1.505_632_735_149_311_6e-7];
    let x = x - 1.0;
    let t = x + 7.5;
    let s: f64 = C[0] + (1..9).map(|i| C[i] / (x + i as f64)).sum::<f64>();
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + s.ln()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn poisson_mean_and_variance() {
        let mut r = Random::new(7);
        for lambda in [3.5, 250.0] {
            let n = 200_000;
            let v: Vec<f64> = (0..n).map(|_| r.poisson(lambda)).collect();
            let m = v.iter().sum::<f64>() / n as f64;
            let var = v.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / n as f64;
            assert!((m - lambda).abs() < 4.0 * (lambda / n as f64).sqrt() + 1e-9, "mean {m} for {lambda}");
            assert!((var / lambda - 1.0).abs() < 0.02, "variance {var} for {lambda}");
        }
        assert!((ln_gamma(10.0) - 362880f64.ln()).abs() < 1e-12);
        // A mean that is not a number draws 0 and returns at once (TU6).
        assert_eq!((r.poisson(f64::NAN), r.poisson(f64::INFINITY), r.poisson(-1.0)), (0.0, 0.0, 0.0));
    }

    #[test]
    fn backgrounds() {
        let b = BackgroundModel::Points(vec![(10.0, 5.0), (20.0, 15.0)]);
        assert_eq!((b.at(0.0), b.at(15.0), b.at(30.0)), (5.0, 10.0, 15.0));
        // Σ B_m (T/Bkpos − 1)^m with Bkpos = 40: at T = 80, u = 1.
        let p = BackgroundModel::Polynomial { origin: 40.0, coefficients: vec![1.0, 2.0, 3.0] };
        assert_eq!(p.at(80.0), 6.0);
        // Chebyshev: T0..T3 at x = −1 alternate in sign.
        let c = BackgroundModel::Chebyshev { from: 10.0, to: 20.0, coefficients: vec![1.0, 1.0, 1.0, 1.0] };
        assert_eq!(c.at(10.0), 0.0);
        assert_eq!(c.at(20.0), 4.0);
    }

    #[test]
    fn resolution_table_interpolates() {
        let r = InstrumentalWidths::Table(vec![[20.0, 0.1, 0.01], [40.0, 0.3, 0.03]]);
        let (hg2, hl) = r.at(30.0);
        assert!((hg2.sqrt() - 0.2).abs() < 1e-12 && (hl - 0.02).abs() < 1e-12);
    }
}
