//! Powder diffraction data: the observed pattern, with standard deviations and any step (constant or variable).
//!
//! The readers are in the registry of `formats` (L062), which knows every format FullProf and GSAS-II read; `read`
//! here is the `.pcr` path, FullProf's `Ins` translated to the registry's format and the first scan of the file
//! taken, as FullProf takes it.

use std::fmt;

/// An observed powder pattern: scattering variable (2θ in degrees for constant wavelength), intensities and their
/// standard deviations.
#[derive(Clone, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Pattern {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    /// In JSON a point with no weight (σ infinite) is `null`, and `null` reads back as no weight (L177): JSON holds no
    /// infinity, and the app reads `null` as no weight already.
    #[cfg_attr(feature = "serde", serde(serialize_with = "no_weight::serialize", deserialize_with = "no_weight::deserialize"))]
    pub sigma: Vec<f64>,
}

/// σ in JSON: a number, or `null` for no weight (an infinite σ), both ways.
#[cfg(feature = "serde")]
mod no_weight {
    use serde::{Deserialize, Deserializer, Serialize, Serializer};

    pub fn serialize<S: Serializer>(v: &[f64], s: S) -> Result<S::Ok, S::Error> {
        v.iter().map(|x| x.is_finite().then_some(*x)).collect::<Vec<Option<f64>>>().serialize(s)
    }

    pub fn deserialize<'de, D: Deserializer<'de>>(d: D) -> Result<Vec<f64>, D::Error> {
        Ok(Vec::<Option<f64>>::deserialize(d)?.into_iter().map(|x| x.unwrap_or(f64::INFINITY)).collect())
    }
}

impl Pattern {
    pub fn len(&self) -> usize {
        self.x.len()
    }
    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }
    /// The mean step between points, (x_last − x_first)/(n − 1) (the outline's and the project's "Step", L145); none
    /// below two points.
    pub fn mean_step(&self) -> Option<f64> {
        (self.len() > 1).then(|| (self.x[self.len() - 1] - self.x[0]) / (self.len() - 1) as f64)
    }
    /// The largest intensity (the project's "Counts up to", Data's "the highest at"; L145).
    pub fn max_count(&self) -> Option<f64> {
        self.y.iter().copied().filter(|v| v.is_finite()).fold(None, |a: Option<f64>, v| Some(a.map_or(v, |m| m.max(v))))
    }
    /// Whether σ is the counting statistics at every point, σ² = max(y, 0) within 10⁻⁶ relative (the outline's "σ:
    /// √counts" against "from the file"; L145).
    pub fn sigma_is_counting(&self) -> bool {
        !self.is_empty() && self.y.iter().zip(&self.sigma).all(|(y, s)| (s * s - y.max(0.0)).abs() < 1e-6 * (1.0 + y.abs()))
    }
}

/// The points a pattern is calculated at when there is no data: from, from + step, … up to `to` (FullProf's Thmin,
/// Step, Thmax; L145). Empty for a step that is not positive.
pub fn grid(from: f64, step: f64, to: f64) -> Vec<f64> {
    if !(step > 0.0) || to < from {
        return Vec::new();
    }
    (0..=((to - from) / step + 1e-9) as usize).map(|k| from + k as f64 * step).collect()
}

/// The prominence proposed for a peak search: a fiftieth of the largest count, at least 1 (Data's Find Peaks dialog;
/// L147).
pub fn proposed_prominence(p: &Pattern) -> Option<f64> {
    p.max_count().map(|m| (m / 50.0).round().max(1.0))
}

/// The smallest and the largest of some abscissae (a selection's range; L147).
pub fn extent(x: &[f64]) -> Option<(f64, f64)> {
    let mut it = x.iter().copied().filter(|v| v.is_finite());
    let first = it.next()?;
    Some(it.fold((first, first), |(lo, hi), v| (lo.min(v), hi.max(v))))
}

/// The point a line of the data file's text is about: its first number read as an abscissa, and the point nearest it
/// if that is within half a step (a line selected in Open as Text selects its point, W6; L147).
pub fn point_of_line(x: &[f64], line: &str) -> Option<usize> {
    let word = line.split_whitespace().next()?;
    let v: f64 = word.replace(['d', 'D'], "e").parse().ok()?;
    let i = crate::diagnostics::nearest_point(x, v)?;
    let step = if x.len() > 1 { (x[1] - x[0]).abs() } else { 1.0 };
    ((x[i] - v).abs() <= step / 2.0).then_some(i)
}

#[cfg(test)]
mod pattern_tests {
    use super::*;

    /// The facts of a pattern the app shows (L145): the mean step, the largest count, whether σ is √counts; the grid.
    #[test]
    fn a_patterns_facts() {
        let p = Pattern { x: vec![10.0, 10.02, 10.04, 10.06], y: vec![4.0, 9.0, 100.0, 0.0], sigma: vec![2.0, 3.0, 10.0, 0.0] };
        assert!((p.mean_step().unwrap() - 0.02).abs() < 1e-12);
        assert_eq!(p.max_count(), Some(100.0));
        assert!(p.sigma_is_counting());
        let q = Pattern { sigma: vec![2.0, 3.0, 11.0, 0.0], ..p.clone() };
        assert!(!q.sigma_is_counting());
        assert_eq!(Pattern::default().mean_step(), None);
        assert_eq!(grid(15.0, 0.5, 16.0), vec![15.0, 15.5, 16.0]);
        assert!(grid(15.0, 0.0, 16.0).is_empty());
        assert_eq!(proposed_prominence(&p), Some(2.0));
        assert_eq!(extent(&[3.0, 1.0, 2.0]), Some((1.0, 3.0)));
        assert_eq!(extent(&[]), None);
        assert_eq!(point_of_line(&p.x, "  10.04  100  10"), Some(2));
        assert_eq!(point_of_line(&p.x, "1.0040D+01 100"), Some(2), "FORTRAN's D exponent");
        assert_eq!(point_of_line(&p.x, "10.5 1"), None, "beyond half a step");
        assert_eq!(point_of_line(&p.x, "# a comment"), None);
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct DataError(pub String);

impl fmt::Display for DataError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str(&self.0)
    }
}

/// Read a data file in FullProf format `ins` (the `Ins` flag of the `.pcr` file): the first scan of the file, read
/// by the registry's format for that `Ins`.
pub fn read(text: &str, ins: i32) -> Result<Pattern, DataError> {
    let f = crate::formats::by_ins(ins).ok_or_else(|| DataError(format!("data format Ins={ins} is not in FullProf's list")))?;
    if let Some(why) = f.not_read {
        return Err(DataError(format!("data format Ins={ins} ({}): {why}", f.id)));
    }
    let probe = crate::formats::Probe::new("data", text.as_bytes());
    let d = (f.read)(&probe)?;
    d.scans.into_iter().next().map(|s| s.pattern).ok_or(DataError(format!("Ins={ins}: no pattern in the file")))
}
