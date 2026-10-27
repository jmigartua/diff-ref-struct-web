//! The user's notes on a refinement's reflections (APP_WORKFLOWS W8, L142): "impurity", "not in the model". A note is
//! scientific content, kept in the visible file the file dictionary names `{control} notes.toml` beside the control
//! file (`reflection_notes`), and carried into the results file (`results::with_notes`). A note is on a reflection
//! (phase and h k l) or on a place of the pattern (2θ). The command line: `drs notes <control>`.
//!
//! ```toml
//! # Notes on the reflections of y2o3_structure.drc (diff-ref-struct, W8)
//! [[note]]
//! phase = 1
//! hkl = [2, 2, 2]
//! text = "a shoulder on the low-angle side"
//!
//! [[note]]
//! two_theta = 31.25
//! text = "a peak no phase has: an impurity?"
//! ```

use crate::control::syntax::{self, Value};

#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Note {
    /// The phase, from 1, of the reflection the note is on.
    pub phase: Option<usize>,
    pub hkl: Option<[i32; 3]>,
    /// The place of the pattern the note is on, when it is on no reflection.
    pub two_theta: Option<f64>,
    pub text: String,
}

impl Note {
    /// Whether the note is on this reflection.
    pub fn is_on(&self, phase: usize, hkl: [i32; 3]) -> bool {
        self.phase == Some(phase) && self.hkl == Some(hkl)
    }
}

/// The notes of a notes file's text, in its order. Each note needs a text and a place: phase with hkl, or two_theta.
pub fn read(text: &str) -> Result<Vec<Note>, String> {
    let doc = syntax::parse(text);
    if let Some(d) = doc.diagnostics.iter().find(|d| !d.warning) {
        return Err(d.to_string());
    }
    let mut out = Vec::new();
    for t in doc.tables.iter().filter(|t| t.array && t.path == ["note"]) {
        let at = t.header.map_or(String::new(), |s| format!("line {}: ", s.line));
        let int = |k: &str| match t.get(k).map(|e| &e.value.value) {
            Some(Value::Integer(i)) => Ok(Some(*i)),
            None => Ok(None),
            Some(v) => Err(format!("{at}{k} is {}, not an integer", v.kind())),
        };
        let phase = int("phase")?.map(|p| p.max(0) as usize);
        let hkl = match t.get("hkl").map(|e| &e.value.value) {
            Some(Value::Array(a)) if a.len() == 3 => {
                let mut h = [0; 3];
                for (j, it) in a.iter().enumerate() {
                    match it.value {
                        Value::Integer(i) => h[j] = i as i32,
                        _ => return Err(format!("{at}hkl holds three integers")),
                    }
                }
                Some(h)
            }
            None => None,
            Some(_) => return Err(format!("{at}hkl holds three integers")),
        };
        let two_theta = match t.get("two_theta").map(|e| &e.value.value) {
            Some(Value::Float(x)) => Some(*x),
            Some(Value::Integer(i)) => Some(*i as f64),
            None => None,
            Some(v) => return Err(format!("{at}two_theta is {}, not a number", v.kind())),
        };
        let text = match t.get("text").map(|e| &e.value.value) {
            Some(Value::String(s)) => s.clone(),
            _ => return Err(format!("{at}a note needs its text (text = \"…\")")),
        };
        if !(phase.is_some() && hkl.is_some()) && two_theta.is_none() {
            return Err(format!("{at}a note is on a reflection (phase and hkl) or at a 2θ (two_theta)"));
        }
        out.push(Note { phase, hkl, two_theta, text });
    }
    Ok(out)
}

/// The notes file's text for `notes`, headed by the control file's name.
pub fn write(notes: &[Note], control: &str) -> String {
    let mut s = format!("# Notes on the reflections of {control} (diff-ref-struct, W8): one [[note]] each, on a reflection\n# (phase, hkl) or at a 2θ (two_theta).\n");
    for n in notes {
        s.push_str("\n[[note]]\n");
        if let Some(p) = n.phase {
            s.push_str(&format!("phase = {p}\n"));
        }
        if let Some([h, k, l]) = n.hkl {
            s.push_str(&format!("hkl = [{h}, {k}, {l}]\n"));
        }
        if let Some(t) = n.two_theta {
            s.push_str(&format!("two_theta = {}\n", syntax::float(t)));
        }
        s.push_str(&format!("text = {}\n", syntax::string(&n.text)));
    }
    s
}

/// The notes with the note on (phase, hkl) set to `text`, or removed when `text` is empty.
pub fn set_on_reflection(notes: &[Note], phase: usize, hkl: [i32; 3], text: &str) -> Vec<Note> {
    let mut out: Vec<Note> = notes.iter().filter(|n| !n.is_on(phase, hkl)).cloned().collect();
    if !text.trim().is_empty() {
        out.push(Note { phase: Some(phase), hkl: Some(hkl), two_theta: None, text: text.trim().to_string() });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A notes file written and read back is the same notes; a note on a reflection is replaced and removed; a note
    /// with no place, or no text, is an error that says why.
    #[test]
    fn notes_written_and_read_back() {
        let n = vec![
            Note { phase: Some(1), hkl: Some([2, 2, 2]), two_theta: None, text: "a shoulder on the \"low\" side".into() },
            Note { phase: None, hkl: None, two_theta: Some(31.25), text: "an impurity?".into() },
        ];
        let t = write(&n, "y2o3.drc");
        assert!(t.starts_with("# Notes on the reflections of y2o3.drc"));
        assert_eq!(read(&t).unwrap(), n);
        let m = set_on_reflection(&n, 1, [2, 2, 2], "not in the model");
        assert_eq!(m.len(), 2);
        assert_eq!(m.iter().find(|x| x.is_on(1, [2, 2, 2])).unwrap().text, "not in the model");
        assert_eq!(set_on_reflection(&m, 1, [2, 2, 2], "  ").len(), 1, "an empty text removes the note");
        assert!(read("[[note]]\ntext = \"x\"\n").unwrap_err().contains("on a reflection"));
        assert!(read("[[note]]\nphase = 1\nhkl = [1, 1, 0]\n").unwrap_err().contains("text"));
        assert_eq!(read("# nothing yet\n").unwrap(), vec![]);
    }
}
