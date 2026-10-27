//! Where a control file writes a value (the app's request 40, L299): for an address of the file's model
//! (`phase1.cell_a`, `phase1.Y:2.atom_position.x`, `pattern1.zero_shift`), the line and the bytes of the value's text,
//! from the reader's own syntax tree. The keys are the parameter dictionary's ids (`refinable::dictionary_id`): a
//! pattern's value in `[pattern]`, its background in `[pattern.background]`, a phase's in its `[[phase]]`, an atom's in
//! the phase's `[[phase.atom]]` in order, a mode's amplitude in the phase's `[[phase.mode]]` in order; a component is
//! an element of the value's array (x of `atom_position`, P2 of `berar_baldinozzi`, the height of a background point,
//! the value of a `quartic_strain` row).
//!
//! With the control file's folder (`locate_in`, `all`, `address_at`; the app's request 41, L305) the file is read as
//! `schema::read_with` reads it, so that a file whose background points are kept in a file of their own has its
//! model; `all` gives every value's place from one reading, and `address_at` the value written at a byte.

use super::{schema, syntax};
use crate::jacobian::Field;
use crate::powder::Model;
use std::ops::Range;
use std::path::Path;

/// The line (from 1) and the bytes of the text where `text` writes the value `address` names; None when the text does
/// not read as a control file, the address names no value of its model, or the file does not write it (a value left
/// at its default: the micro-absorption at zero, background points kept in a file of their own). The file is read
/// without its folder: see `locate_in`.
pub fn locate(text: &str, address: &str) -> Option<(usize, Range<usize>)> {
    locate_in(text, address, None)
}

/// `locate`, the file read with its folder (`schema::read_with`): the files it names, such as background points kept in
/// a file of their own, are read from there.
pub fn locate_in(text: &str, address: &str, dir: Option<&Path>) -> Option<(usize, Range<usize>)> {
    let places = Places::read(text, dir)?;
    let f = crate::refinable::field_of(address, &places.model).ok()?;
    places.of(f)
}

/// Every value the file writes, from one reading: (address, line, bytes), in the order of `Model::fields`.
pub fn all(text: &str, dir: Option<&Path>) -> Vec<(String, usize, Range<usize>)> {
    let Some(places) = Places::read(text, dir) else { return Vec::new() };
    places.model.fields().into_iter().filter_map(|fi| places.of(fi.field).map(|(line, r)| (fi.address, line, r))).collect()
}

/// The address of the value written at byte `byte` of the text (its first byte up to just after its last), from one
/// reading; None when no value of the model is written there.
pub fn address_at(text: &str, byte: usize, dir: Option<&Path>) -> Option<String> {
    all(text, dir).into_iter().find(|(_, _, r)| r.start <= byte && byte <= r.end).map(|x| x.0)
}

/// A control file read once: its model and its syntax tree, the tables of each phase in the file's order.
struct Places {
    model: Model,
    doc: syntax::Document,
}

impl Places {
    fn read(text: &str, dir: Option<&Path>) -> Option<Places> {
        let doc = syntax::parse(text);
        let (control, _) = schema::read_with(&doc, dir);
        Some(Places { model: control?.model, doc })
    }

    /// The line and bytes where the file writes field `f`.
    fn of(&self, f: Field) -> Option<(usize, Range<usize>)> {
        let (m, doc) = (&self.model, &self.doc);
        let key = crate::refinable::dictionary_id(f, m);
        // The tables of each phase, in the file's order: the phase's own, its atoms' and its modes'.
        let (mut phases, mut atoms, mut modes): (Vec<&syntax::Table>, Vec<Vec<&syntax::Table>>, Vec<Vec<&syntax::Table>>) = (Vec::new(), Vec::new(), Vec::new());
        let mut groups: Vec<Vec<&syntax::Table>> = Vec::new();
        for t in &doc.tables {
            match (t.path.iter().map(String::as_str).collect::<Vec<_>>().as_slice(), t.array) {
                (["phase"], true) => {
                    phases.push(t);
                    atoms.push(Vec::new());
                    modes.push(Vec::new());
                    groups.push(Vec::new());
                }
                (["phase", "rigid_group"], true) => groups.last_mut()?.push(t),
                (["phase", "atom"], true) => atoms.last_mut()?.push(t),
                (["phase", "mode"], true) => modes.last_mut()?.push(t),
                _ => {}
            }
        }
        let named = |path: &[&str]| doc.tables.iter().find(|t| !t.array && t.path.iter().map(String::as_str).eq(path.iter().copied()));
        // The table, and the element of the value's array (with the element within a row) where the field is one.
        let (table, element): (&syntax::Table, Option<(usize, Option<usize>)>) = match f {
            Field::Zero | Field::Sycos | Field::Sysin | Field::Ratio | Field::Wavelength | Field::Micro(_) => (named(&["pattern"])?, None),
            Field::Background(k) => {
                let row = matches!(m.background, crate::powder::BackgroundModel::Points(_) | crate::powder::BackgroundModel::Spline(_));
                (named(&["pattern", "background"])?, Some((k, row.then_some(1))))
            }
            // A further term of the background (P23): the element of its key's array.
            Field::BackgroundTerm(k) => {
                let part = m.background_terms.part(k)?;
                // A file's b1–b4 (the app's request 47): each its own key, a number.
                let element = (part.0 != crate::background_terms::Part::File).then_some((part.1, None));
                (named(&["pattern", "background"])?, element)
            }
            Field::AtomX(p, a, k) | Field::AtomBeta(p, a, k) => (*atoms.get(p)?.get(a)?, Some((k, None))),
            Field::AtomB(p, a) | Field::AtomOcc(p, a) => (*atoms.get(p)?.get(a)?, None),
            Field::Mode(p, k) => (*modes.get(p)?.get(k)?, None),
            // A rigid group's origin and orientation in its table; an internal coordinate in its atom's, the element of
            // its rigid_length or rigid_angle (P41).
            Field::Rigid(p, g, k) if k < 6 => (*groups.get(p)?.get(g)?, Some((k % 3, None))),
            // A TLS value (L367): the element of rigid_tls_t, rigid_tls_l (six each) or rigid_tls_s (nine) in the
            // group's table.
            Field::Rigid(p, g, k) if m.phases.get(p)?.rigid.get(g).and_then(|gr| crate::rigid_body::tls_index(gr, k)).is_some() => {
                let i = crate::rigid_body::tls_index(&m.phases[p].rigid[g], k)?;
                (*groups.get(p)?.get(g)?, Some((if i < 12 { i % 6 } else { i - 12 }, None)))
            }
            Field::Rigid(p, g, k) => {
                let gr = m.phases.get(p)?.rigid.get(g)?;
                let (len, ang) = gr.coordinates.split();
                let c = (k - 6) % 3;
                let i = len.iter().position(|&x| x == c).or_else(|| ang.iter().position(|&x| x == c))?;
                (*atoms.get(p)?.get(gr.members.get((k - 6) / 3)?.index)?, Some((i, None)))
            }
            Field::Asym(p, k) => (*phases.get(p)?, Some((k, None))),
            // An axis's r or fraction (P30): the element of po_axis_ratio or po_axis_fraction.
            Field::Texture(p, k) => (*phases.get(p)?, Some((k / 2, None))),
            Field::Quartic(p, k) => (*phases.get(p)?, Some((k, Some(3)))),
            // A Pawley fit's intensity (P36): the fifth element of its row, [h, k, l, multiplicity, intensity].
            Field::Intensity(p, k) => (*phases.get(p)?, Some((k, Some(4)))),
            _ => (*phases.get(f.phase()?)?, None),
        };
        let mut item = &table.get(key)?.value;
        if let Some((i, inner)) = element {
            let syntax::Value::Array(a) = &item.value else { return None };
            item = a.get(i)?;
            if let Some(j) = inner {
                let syntax::Value::Array(b) = &item.value else { return None };
                item = b.get(j)?;
            }
        }
        Some((item.span.line, item.span.start..item.span.end))
    }
}
