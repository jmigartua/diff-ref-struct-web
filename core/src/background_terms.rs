//! Further terms of the background (P23, L353), added to the background function (`Model::background`):
//!
//! - Debye-like terms Σ_j A_j sin(Q r_j)/(Q r_j), Q = 4π sin θ/λ₁ at the data's 2θ (FullProf's Nba = −1, manual eq.
//!   3.4, p. 45: the line pairs "Additional background coefficients", the amplitudes, then the distances r_j);
//! - inverse powers of 2θ, Σ_{j=1…3} c_j (2θ)^−j, 2θ in degrees (FullProf's Nba = −3, the last three of its twelve
//!   coefficients after nine of the polynomial in 2θ/Bkpos − 1; the News of 28 July 2003 gives the formula as an image
//!   that is lost: the form is read from FullProf's runs);
//! - profiles read from files, each times its coefficient (FullProf's MORE = 2 and 3 on line 14, "Additional
//!   background profiles", profile n from `{data}_n.bac`), joined linearly;
//! - a file background's scale, offset, stretch and shift, b1·B[(1 + b3)2θ + b4] + b2 (FullProf's Nba = 1, manual eq.
//!   3.6, p. 46): the function's points are the file's as read, and these four values move them (`FileScale`; the
//!   app's request 47, Publish's P57). Not added: they transform the function (`Model::background_unbounded`).
//!
//! Against FullProf's runs (`runs/background_terms`, L346): its background within its print at every point but two of
//! nac_3t2's (test `background_terms`, L353).
//!
//! The values are numbered for `Field::BackgroundTerm(k)`: each Debye-like term's amplitude and distance (2j, 2j + 1),
//! then the inverse powers' coefficients, then the profiles' coefficients, then a file's b1, b2, b3, b4.

/// One Debye-like term: amplitude A (counts) and distance r (Å).
#[derive(Clone, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct DebyeTerm {
    pub amplitude: f64,
    pub distance: f64,
}

/// A background profile read from a file, times its coefficient.
#[derive(Clone, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Profile {
    /// The file as the control file names it (relative to it), or as the `.pcr` import found it.
    pub file: String,
    pub coefficient: f64,
    /// The file's points (x, counts), joined linearly; outside them the end point's value.
    pub points: Vec<(f64, f64)>,
}

/// A file background's b1–b4 (FullProf's Nba = 1, eq. 3.6): the function's points (the file's, as read) are read at
/// (1 + stretch)·2θ + shift, times scale, plus offset. FullProf reads b1 = 0 as 1; the import stores 1.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct FileScale {
    /// The file the points were read from, as the `.pcr` import found it or the control file names it.
    pub file: String,
    pub scale: f64,
    pub offset: f64,
    pub stretch: f64,
    pub shift: f64,
}

impl Default for FileScale {
    fn default() -> Self {
        FileScale { file: String::new(), scale: 1.0, offset: 0.0, stretch: 0.0, shift: 0.0 }
    }
}

impl FileScale {
    /// b1–b4 in FullProf's order.
    pub fn values(&self) -> [f64; 4] {
        [self.scale, self.offset, self.stretch, self.shift]
    }

    /// Whether the four leave the points as they are (b1 = 1, the others 0).
    pub fn is_identity(&self) -> bool {
        self.values() == [1.0, 0.0, 0.0, 0.0]
    }

    /// Where the function is read for the data's x: (1 + b3)·x + b4.
    pub fn at(&self, x: f64) -> f64 {
        (1.0 + self.stretch) * x + self.shift
    }
}

/// The further terms of a pattern's background. Empty: none.
#[derive(Clone, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct BackgroundTerms {
    pub debye: Vec<DebyeTerm>,
    /// c_1, c_2, c_3 of Σ c_j (2θ)^−j (fewer when the last ones are not given).
    pub inverse: Vec<f64>,
    pub profiles: Vec<Profile>,
    /// A file background's b1–b4 (Nba = 1; the app's request 47).
    #[cfg_attr(feature = "serde", serde(default))]
    pub file: Option<FileScale>,
}

/// sin(x)/x and its derivative, by their series near 0.
fn sinc(x: f64) -> (f64, f64) {
    if x.abs() < 1e-4 {
        (1.0 - x * x / 6.0, -x / 3.0)
    } else {
        let (s, c) = x.sin_cos();
        (s / x, (c - s / x) / x)
    }
}

impl BackgroundTerms {
    pub fn is_empty(&self) -> bool {
        self.debye.is_empty() && self.inverse.is_empty() && self.profiles.is_empty() && self.file.is_none()
    }

    /// The number of values (`Field::BackgroundTerm`).
    pub fn len(&self) -> usize {
        2 * self.debye.len() + self.inverse.len() + self.profiles.len() + if self.file.is_some() { 4 } else { 0 }
    }

    /// What value k is: (part, term 0, 1, …, which value of a Debye-like term: 0 amplitude, 1 distance).
    pub fn part(&self, k: usize) -> Option<(Part, usize, usize)> {
        let nd = 2 * self.debye.len();
        let ni = nd + self.inverse.len();
        if k < nd {
            Some((Part::Debye, k / 2, k % 2))
        } else if k < ni {
            Some((Part::Inverse, k - nd, 0))
        } else if k < ni + self.profiles.len() {
            Some((Part::Profile, k - ni, 0))
        } else if self.file.is_some() && k < ni + self.profiles.len() + 4 {
            Some((Part::File, 0, k - ni - self.profiles.len()))
        } else {
            None
        }
    }

    pub fn value(&self, k: usize) -> f64 {
        match self.part(k) {
            Some((Part::Debye, j, 0)) => self.debye[j].amplitude,
            Some((Part::Debye, j, _)) => self.debye[j].distance,
            Some((Part::Inverse, j, _)) => self.inverse[j],
            Some((Part::Profile, j, _)) => self.profiles[j].coefficient,
            Some((Part::File, _, i)) => self.file.as_ref().map_or(0.0, |f| f.values()[i]),
            None => 0.0,
        }
    }

    pub fn set(&mut self, k: usize, v: f64) {
        match self.part(k) {
            Some((Part::Debye, j, 0)) => self.debye[j].amplitude = v,
            Some((Part::Debye, j, _)) => self.debye[j].distance = v,
            Some((Part::Inverse, j, _)) => self.inverse[j] = v,
            Some((Part::Profile, j, _)) => self.profiles[j].coefficient = v,
            Some((Part::File, _, i)) => {
                if let Some(f) = &mut self.file {
                    *[&mut f.scale, &mut f.offset, &mut f.stretch, &mut f.shift][i] = v;
                }
            }
            None => {}
        }
    }

    /// The dictionary id and the component of value k (`refinable::dictionary_id`, `address`).
    pub fn id(&self, k: usize) -> (&'static str, String) {
        match self.part(k) {
            Some((Part::Debye, j, 0)) => ("background_debye_amplitude", (j + 1).to_string()),
            Some((Part::Debye, j, _)) => ("background_debye_distance", (j + 1).to_string()),
            Some((Part::Inverse, j, _)) => ("background_inverse_power", (j + 1).to_string()),
            Some((Part::Profile, j, _)) => ("background_profile", (j + 1).to_string()),
            Some((Part::File, _, i)) => (["background_file_scale", "background_file_offset", "background_file_stretch", "background_file_shift"][i], "1".to_string()),
            None => ("background_profile", (k + 1).to_string()),
        }
    }

    /// Debye-like term j at Q (Å⁻¹).
    pub fn debye_at(&self, j: usize, q: f64) -> f64 {
        let d = &self.debye[j];
        d.amplitude * sinc(q * d.distance).0
    }

    /// The inverse powers at x (2θ in degrees); 0 at x ≤ 0.
    pub fn inverse_at(&self, x: f64) -> f64 {
        if x <= 0.0 {
            return 0.0;
        }
        self.inverse.iter().enumerate().map(|(j, c)| c * x.powi(-(j as i32 + 1))).sum()
    }

    /// Profile j at x.
    pub fn profile_at(&self, j: usize, x: f64) -> f64 {
        let p = &self.profiles[j];
        p.coefficient * crate::background::linear_at(&p.points, x)
    }

    /// All the terms at x, Q.
    pub fn at(&self, x: f64, q: f64) -> f64 {
        (0..self.debye.len()).map(|j| self.debye_at(j, q)).sum::<f64>() + self.inverse_at(x) + (0..self.profiles.len()).map(|j| self.profile_at(j, x)).sum::<f64>()
    }

    /// The derivative of the terms at x, Q with respect to value k (none for a file's b1–b4, which transform the
    /// function: `Model::background_file_derivative`).
    pub fn derivative(&self, k: usize, x: f64, q: f64) -> f64 {
        match self.part(k) {
            Some((Part::Debye, j, 0)) => sinc(q * self.debye[j].distance).0,
            Some((Part::Debye, j, _)) => self.debye[j].amplitude * q * sinc(q * self.debye[j].distance).1,
            Some((Part::Inverse, j, _)) if x > 0.0 => x.powi(-(j as i32 + 1)),
            Some((Part::Profile, j, _)) => crate::background::linear_at(&self.profiles[j].points, x),
            _ => 0.0,
        }
    }
}

/// The parts of the terms.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Part {
    Debye,
    Inverse,
    Profile,
    File,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn derivatives_by_differences() {
        let t = BackgroundTerms {
            debye: vec![DebyeTerm { amplitude: -37.0, distance: 1.75 }, DebyeTerm { amplitude: 33.0, distance: 2.5 }],
            inverse: vec![5.0, 40.0, 300.0],
            profiles: vec![Profile { file: "a.bac".into(), coefficient: 0.5, points: vec![(10.0, 2.0), (20.0, 4.0)] }],
            file: None,
        };
        for (x, q) in [(12.0, 1.07), (35.0, 2.9), (1.0, 1e-6)] {
            for k in 0..t.len() {
                let h = 1e-6 * t.value(k).abs().max(1.0);
                let (mut a, mut b) = (t.clone(), t.clone());
                a.set(k, t.value(k) + h);
                b.set(k, t.value(k) - h);
                let num = (a.at(x, q) - b.at(x, q)) / (2.0 * h);
                assert!((num - t.derivative(k, x, q)).abs() < 1e-6 * num.abs().max(1.0), "k {k} at {x}: {num} {}", t.derivative(k, x, q));
            }
        }
    }
}
