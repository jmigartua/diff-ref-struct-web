//! The model written and read as text (the app's request 1, L058), behind the crate's `serde` feature: JSON of the
//! model as the core holds it (phases, atoms, background, instrument, rules, restraints; the space groups as their
//! Hall symbols and operators). The data are not part of it. The readable control file of P47 comes later.

use crate::powder::Model;

impl Model {
    /// The model as text (pretty-printed JSON).
    pub fn to_text(&self) -> String {
        serde_json::to_string_pretty(self).expect("a model is always written")
    }

    /// A model from its text.
    pub fn from_text(text: &str) -> Result<Model, String> {
        serde_json::from_str(text).map_err(|e| format!("model text: {e}"))
    }
}
