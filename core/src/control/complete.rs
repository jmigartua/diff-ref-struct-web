//! Completion for the control file's editor (P49): the names valid at a place of a `.drc`, each with what the
//! parameter dictionary says of it.
//!
//! At a line that opens a table (`[` typed): the tables a control file may hold (`schema::TABLES`), less the single
//! tables already written. At a key: the keys of the table the line is in, less those it already holds, each with
//! its dictionary name, unit and meaning when the key is a parameter's id (`parameters::lookup`), or the operation's
//! setting in `[[pattern.operation]]`. At a value: the closed choices the reader accepts (`schema::choices`), the
//! scattering species of the tables for `species`, the 530 settings' symbols for `space_group`, and for the
//! refinement's addresses (`parameters`, a step's `free`, a link's `value` and `follows`) every address of the
//! model the file describes that can be set by itself (`Model::fields`), when the rest of the file reads.
//!
//! Checked in `completion`: every key offered is one the reader accepts in that table and every key it accepts is
//! offered; every table offered is one it accepts; every address offered is one `refinable::field_of` reads back.

use super::schema::{choices, TABLES};
use super::syntax;

/// What a completion is.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum CompletionKind {
    Table,
    Key,
    Value,
    Address,
}

/// One name valid at the place.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Completion {
    pub label: String,
    pub kind: CompletionKind,
    /// A short line: the parameter's name and unit, the table's role.
    pub detail: String,
    /// The dictionary's meaning, when there is one.
    pub documentation: String,
    /// The text to insert in place of the word being typed.
    pub insert: String,
}

/// A table's header (`[a.b]` or `[[a.b]]`, a comment after it allowed): its path and whether it is an array's.
fn header(line: &str) -> Option<(Vec<String>, bool)> {
    let t = code(line).trim();
    let (inner, array) = match t.strip_prefix("[[").and_then(|r| r.strip_suffix("]]")) {
        Some(i) => (i, true),
        None => (t.strip_prefix('[')?.strip_suffix(']')?, false),
    };
    let inner = inner.trim();
    let ok = inner.starts_with(|c: char| c.is_ascii_alphabetic()) && inner.chars().all(|c| c.is_ascii_alphanumeric() || matches!(c, '_' | '-' | '.' | ' '));
    ok.then(|| (inner.split('.').map(|s| s.trim().to_string()).collect(), array))
}

/// A line without its comment (a `#` outside a string).
fn code(line: &str) -> &str {
    let mut quoted = false;
    for (i, c) in line.char_indices() {
        match c {
            '"' => quoted = !quoted,
            '#' if !quoted => return &line[..i],
            _ => {}
        }
    }
    line
}

/// The brackets a piece of a line opens less those it closes, strings and comments aside.
fn depth(piece: &str) -> i32 {
    let mut quoted = false;
    let mut d = 0;
    for c in code(piece).chars() {
        match c {
            '"' => quoted = !quoted,
            '[' if !quoted => d += 1,
            ']' if !quoted => d -= 1,
            _ => {}
        }
    }
    d
}

/// The table a line is in: the last header above it, with its line; the root when none.
fn table_at(lines: &[&str], line: usize) -> (Vec<String>, bool, usize) {
    (0..line).rev().find_map(|i| header(lines[i]).map(|(p, a)| (p, a, i))).unwrap_or((Vec::new(), false, 0))
}

/// The statements of a table's lines `from..to`: (key, first line, last line), an array's value over several lines
/// held together.
fn statements(lines: &[&str], from: usize, to: usize) -> Vec<(String, usize, usize)> {
    let mut v: Vec<(String, usize, usize)> = Vec::new();
    let mut open = 0;
    for (i, l) in lines.iter().enumerate().take(to).skip(from) {
        if open > 0 {
            open += depth(l);
            if let Some(last) = v.last_mut() {
                last.2 = i;
            }
            continue;
        }
        if let Some((k, rest)) = code(l).split_once('=') {
            v.push((k.trim().to_string(), i, i));
            open = depth(rest).max(0);
        }
    }
    v
}

fn key_completion(key: &str) -> Completion {
    let (detail, documentation) = match crate::parameters::lookup(key) {
        Some(p) => (format!("{} ({})", p.name, p.unit), p.meaning.to_string()),
        None => (String::new(), String::new()),
    };
    Completion { label: key.to_string(), kind: CompletionKind::Key, detail, documentation, insert: format!("{key} = ") }
}

/// The scattering species of the tables, as a control file writes them ("Fe+3", "O", "He").
fn species() -> Vec<String> {
    let mut v: Vec<String> = crate::scattering_table::XRAY_FORM
        .iter()
        .map(|f| {
            let mut c = f.0.chars();
            c.next().map(|h| h.to_uppercase().collect::<String>() + c.as_str()).unwrap_or_default()
        })
        .collect();
    for n in crate::scattering_table::NEUTRON_B.iter() {
        let mut c = n.0.chars();
        let s = c.next().map(|h| h.to_string() + &c.as_str().to_lowercase()).unwrap_or_default();
        if !v.contains(&s) {
            v.push(s);
        }
    }
    v
}

/// The completions at a place of a control file's text: `line` and `column` from 1, the column counted in characters.
pub fn complete(text: &str, line: usize, column: usize) -> Vec<Completion> {
    let lines: Vec<&str> = text.split('\n').collect();
    let Some(current) = lines.get(line.saturating_sub(1)) else { return Vec::new() };
    let prefix: String = current.chars().take(column.saturating_sub(1)).collect();
    let trimmed = prefix.trim_start();
    let (path, array, top) = table_at(&lines, line - 1);
    // The table's statements (from its header to the next), and the one the place is in, if any: a key's line, or a
    // line inside an array value opened above.
    let end = (line..lines.len()).find(|&i| header(lines[i]).is_some()).unwrap_or(lines.len());
    let first = if header(lines.get(top).copied().unwrap_or("")).is_some() { top + 1 } else { 0 };
    let stmts = statements(&lines, first, end);
    let at = stmts.iter().find(|s| s.1 < line && line - 1 <= s.2).cloned();
    // A table's header, at a "[" outside a value.
    if trimmed.starts_with('[') && at.is_none() {
        let typed = trimmed.trim_start_matches('[');
        let doc = syntax::parse(text);
        return TABLES
            .iter()
            .filter(|(path, _, _)| !path.is_empty())
            .filter(|(path, array, _)| *array || !doc.tables.iter().any(|t| !t.array && t.path == *path))
            .map(|(path, array, _)| {
                let name = path.join(".");
                let label = if *array { format!("[[{name}]]") } else { format!("[{name}]") };
                Completion { label: label.clone(), kind: CompletionKind::Table, detail: if *array { "a table that may be repeated".into() } else { "a table".into() }, documentation: String::new(), insert: label }
            })
            .filter(|c| c.label.trim_start_matches('[').starts_with(typed))
            .collect();
    }
    let here: Vec<&str> = path.iter().map(String::as_str).collect();
    let Some(&(tpath, _, keys)) = TABLES.iter().find(|(p, a, _)| *p == here.as_slice() && (*a == array || p.is_empty())) else { return Vec::new() };
    let held: Vec<String> = stmts.iter().filter(|s| at.as_ref() != Some(*s)).map(|s| s.0.clone()).collect();
    // The text before the place within its statement's value, when the place is in a value.
    let value_before = match &at {
        Some((k, f, _)) if *f < line - 1 => Some((k.clone(), format!("{}\n{}", lines[*f..line - 1].iter().map(|l| code(l)).collect::<Vec<_>>().join("\n"), prefix))),
        Some((k, _, _)) => prefix.split_once('=').map(|(_, v)| (k.clone(), v.to_string())),
        None => None,
    };
    let (statement_first, statement_last) = at.as_ref().map_or((line - 1, line - 1), |s| (s.1, s.2));
    match value_before {
        // A key.
        None => {
            let typed = trimmed.trim();
            let mut out: Vec<Completion> = keys.iter().filter(|k| !held.iter().any(|h| h == *k) && k.starts_with(typed)).map(|k| key_completion(k)).collect();
            if tpath == ["pattern", "operation"] {
                // The settings of the operation the table names.
                let op = (first..end).find_map(|i| lines[i].split_once('=').filter(|(k, _)| k.trim() == "op").map(|(_, v)| v.trim().trim_matches('"').to_string()));
                if let Some(entry) = op.and_then(|o| crate::operation_table::OPERATIONS.iter().find(|e| e.id == o)) {
                    for st in entry.settings.iter().filter(|st| !held.iter().any(|h| h == st.key) && st.key.starts_with(typed)) {
                        out.push(Completion { label: st.key.to_string(), kind: CompletionKind::Key, detail: format!("{} ({})", st.name, st.unit), documentation: st.meaning.to_string(), insert: format!("{} = ", st.key) });
                    }
                }
            }
            out
        }
        // A value.
        Some((key, value)) => {
            let key = key.trim();
            let typed = value.rsplit([',', '[', '\n']).next().unwrap_or("").trim_start().trim_start_matches('"').to_string();
            let mut out = Vec::new();
            if let Some((_, _, vals)) = choices().into_iter().find(|(p, k, _)| *p == tpath && *k == key) {
                for v in vals.into_iter().filter(|v| v.starts_with(&typed)) {
                    let insert = if v == "true" || v == "false" { v.to_string() } else { format!("\"{v}\"") };
                    out.push(Completion { label: v.to_string(), kind: CompletionKind::Value, detail: String::new(), documentation: String::new(), insert });
                }
            } else if key == "species" {
                for s in species().into_iter().filter(|s| s.starts_with(&typed)) {
                    out.push(Completion { label: s.clone(), kind: CompletionKind::Value, detail: "a scattering species of the tables".into(), documentation: String::new(), insert: format!("\"{s}\"") });
                }
            } else if key == "space_group" {
                // Each setting by a symbol that reads back as that setting: the short symbol, else the full one, else the
                // short one with the setting's choice ("P n -3 n:1").
                let norm = |x: &str| x.replace('_', "").to_uppercase();
                for st in crate::sgtable::SETTINGS.iter().filter(|st| norm(st.hm).starts_with(&norm(&typed)) || norm(st.hm_full).starts_with(&norm(&typed))) {
                    let reads = |sym: &str| crate::spacegroup::SpaceGroup::from_fullprof_symbol(sym).is_ok_and(|g| g.hall == st.hall);
                    let with_choice = format!("{}:{}", st.hm, st.choice);
                    let symbol = [st.hm, st.hm_full, &with_choice].into_iter().find(|s| reads(s)).unwrap_or(st.hm);
                    let detail = format!("No. {}{}", st.number, if st.choice.is_empty() { String::new() } else { format!(", setting {}", st.choice) });
                    out.push(Completion { label: symbol.to_string(), kind: CompletionKind::Value, detail, documentation: format!("Hall symbol {}", st.hall), insert: format!("\"{symbol}\"") });
                }
            } else if matches!((tpath, key), (["refine"], "parameters") | (["refine", "step"], "free") | (["refine", "link"], "value" | "follows")) {
                // The model's addresses, from the file without this line; those that can be set by themselves.
                let without: String = lines.iter().enumerate().filter(|(i, _)| *i < statement_first || *i > statement_last).map(|(_, l)| *l).collect::<Vec<_>>().join("\n");
                // When the rest does not read without the statement (a table just begun), without the whole table.
                let without_table: String = lines.iter().enumerate().filter(|(i, _)| top == 0 && first == 0 || *i < top || *i >= end).map(|(_, l)| *l).collect::<Vec<_>>().join("\n");
                if let Some(c) = super::Control::from_text(&without).ok().or_else(|| super::Control::from_text(&without_table).ok()) {
                    for f in c.model.fields().into_iter().filter(|f| f.fixed_by.is_none() && f.address.starts_with(&typed)) {
                        let (detail, documentation) = crate::parameters::lookup(f.id).map_or((String::new(), String::new()), |p| (format!("{} ({})", p.name, p.unit), p.meaning.to_string()));
                        out.push(Completion { label: f.address.clone(), kind: CompletionKind::Address, detail, documentation, insert: format!("\"{}\"", f.address) });
                    }
                }
            }
            out
        }
    }
}
