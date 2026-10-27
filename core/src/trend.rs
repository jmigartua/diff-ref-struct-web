//! A parameter across the versions of a control file (APP_WORKFLOWS W7, "the same parameter across the versions as a
//! trend"; L141). Each version is a control file kept by the app in `Versions/`, with, when its refinement was
//! accepted, the results file the core wrote for it. The value is read from the control file (the model it holds), the
//! standard uncertainty from the results file (`results::read`, the refined parameters with their su). The command
//! line: `drs trend <address> <control[:results]>…`.

use crate::control::{schema, syntax};
use crate::refinable::field_of;
use std::path::{Path, PathBuf};

/// One version's point of the trend.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TrendPoint {
    /// The value in the version's control file.
    pub value: f64,
    /// Its su from the version's results file; none without one, or when the parameter was not refined there.
    pub su: Option<f64>,
    /// Whether the parameter was refined in that version's run (it is among the results file's refined ones).
    pub refined: bool,
}

/// The value of the value at `address` in a control file's text. `base` is the folder its relative files are read
/// from (a background-points file): the original control file's folder, since a version is a copy kept elsewhere.
pub fn value_in(control_text: &str, base: Option<&Path>, address: &str) -> Result<f64, String> {
    let (c, diags) = schema::read_with(&syntax::parse(control_text), base);
    let c = c.ok_or_else(|| {
        let first = diags.iter().find(|d| !d.warning).map_or(String::new(), |d| format!(": {}", d.message));
        format!("the control file cannot be read{first}")
    })?;
    let f = field_of(address, &c.model)?;
    Ok(c.model.get(f))
}

/// The su of the parameter `address` in a results file's text, and whether it was refined there.
pub fn su_in(results_text: &str, address: &str) -> Result<Option<f64>, String> {
    let r = crate::results::read(results_text)?;
    Ok(r.refined.iter().find(|(n, _, _)| n == address).map(|(_, _, s)| *s))
}

/// The trend of `address` over the versions, oldest first: each version's control file and, if any, its results file.
/// A version that cannot be read, or whose model has no such value (a phase added later), is an error at its place.
pub fn trend(versions: &[(PathBuf, Option<PathBuf>)], base: Option<&Path>, address: &str) -> Vec<Result<TrendPoint, String>> {
    versions
        .iter()
        .map(|(control, results)| {
            let text = std::fs::read_to_string(control).map_err(|e| format!("{}: {e}", control.display()))?;
            let value = value_in(&text, base, address).map_err(|e| format!("{}: {e}", control.display()))?;
            let su = match results {
                Some(r) if r.exists() => {
                    let t = std::fs::read_to_string(r).map_err(|e| format!("{}: {e}", r.display()))?;
                    su_in(&t, address).map_err(|e| format!("{}: {e}", r.display()))?
                }
                _ => None,
            };
            Ok(TrendPoint { value, su, refined: su.is_some() })
        })
        .collect()
}
