//! A control file written again after its model changed, as the person wrote it (L066): its comments, its order and the
//! way each value is written are kept; a value is rewritten only where it changed (compared as values, so "1.50" stays
//! when the value is 1.5); a key the canonical file has and this one lacks is added at the end of its table; a key it
//! has and the canonical file lacks is removed with its line; a table the canonical file has and this one lacks is
//! added after the last table of its kind (or at the end); a table it has and the canonical file lacks is removed.
//!
//! Tables are matched by their path and their place among the tables of that path (the second `[[phase]]`, its third
//! `[[phase.atom]]` or `[[phase.mode]]`).

use super::schema;
use super::syntax::{self, Document, Table};
use super::Control;

/// A table's identity: its path and, for an element of an array of tables, its place among those of its parent.
fn keys(doc: &Document) -> Vec<String> {
    let mut out = Vec::new();
    let mut phase = 0usize;
    let mut counts: std::collections::BTreeMap<String, usize> = std::collections::BTreeMap::new();
    for t in &doc.tables {
        let path = t.path.join(".");
        if path == "phase" && t.array {
            phase += 1;
        }
        let scope = if path == "phase.atom" || path == "phase.mode" { format!("phase{phase}.") } else { String::new() };
        let k = format!("{scope}{path}");
        let n = counts.entry(k.clone()).or_default();
        *n += 1;
        out.push(if t.array { format!("{k}#{n}") } else { k });
    }
    out
}

/// The text of a table as the canonical document has it (its header and its entries).
fn table_text(doc: &Document, t: &Table) -> String {
    let start = t.header_lines.0;
    let mut s = doc.source[start..t.end()].to_string();
    if !s.ends_with('\n') {
        s.push('\n');
    }
    s
}

/// Rewrite `old` (a control file's text) for the control `c`.
pub fn update(old: &str, c: &Control) -> String {
    let new_text = schema::write(c);
    let new = syntax::parse(&new_text);
    let doc = syntax::parse(old);
    if doc.has_errors() {
        // A file that does not parse is not merged into: the canonical text replaces it.
        return new_text;
    }
    let (ok, nk) = (keys(&doc), keys(&new));
    // Edits as (start, end, replacement), applied from the end.
    let mut edits: Vec<(usize, usize, String)> = Vec::new();
    for (i, t) in doc.tables.iter().enumerate() {
        match nk.iter().position(|k| *k == ok[i]) {
            None => {
                // A table the model no longer has: its header and entries go (comments before it stay).
                if t.header.is_some() {
                    edits.push((t.header_lines.0, t.end(), String::new()));
                }
            }
            Some(j) => {
                let nt = &new.tables[j];
                for e in &t.entries {
                    match nt.get(&e.key) {
                        None => edits.push((e.lines.0, e.lines.1, String::new())),
                        Some(ne) if !ne.value.value.same(&e.value.value) => {
                            edits.push((e.value.span.start, e.value.span.end, new_text[ne.value.span.start..ne.value.span.end].to_string()));
                        }
                        Some(_) => {}
                    }
                }
                let added: String = nt.entries.iter().filter(|ne| t.get(&ne.key).is_none()).map(|ne| new_text[ne.lines.0..ne.lines.1].to_string()).collect();
                if !added.is_empty() {
                    let at = t.end();
                    let lead = if at > 0 && !doc.source[..at].ends_with('\n') { "\n" } else { "" };
                    edits.push((at, at, format!("{lead}{added}")));
                }
            }
        }
    }
    // Tables to add: after the last old table that precedes them in the canonical order, or at the end.
    for (j, nt) in new.tables.iter().enumerate() {
        if ok.contains(&nk[j]) {
            continue;
        }
        let anchor = (0..j).rev().find_map(|q| ok.iter().position(|k| *k == nk[q]));
        let at = anchor.map_or(doc.source.len(), |i| doc.tables[i].end());
        let lead = if at > 0 && !doc.source[..at].ends_with('\n') { "\n" } else { "" };
        edits.push((at, at, format!("{lead}\n{}", table_text(&new, nt))));
    }
    // Insertions at the same place keep their canonical order; deletions and replacements do not overlap.
    edits.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    let mut out = String::with_capacity(old.len() + 256);
    let mut pos = 0;
    for (s, e, r) in edits {
        if s < pos {
            continue;
        }
        out.push_str(&doc.source[pos..s]);
        out.push_str(&r);
        pos = e;
    }
    out.push_str(&doc.source[pos..]);
    out
}
