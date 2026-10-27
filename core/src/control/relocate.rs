//! A control file moved to another folder (the app's request 45, File ▸ Move To…, L337): every path it holds rewritten
//! so that it reaches the same files from the new folder, the rest of its text (comments, order, the way each value is
//! written) kept as it was.
//!
//! The paths a control file holds (`paths`), which are the core's to know: the data file (`[pattern] data`), the
//! points of the background (`[pattern.background] file`) and its profiles (`profile_files`, P23), the files of the
//! data operations (each setting whose unit is "path" in the operations' dictionary: `file` of add, subtract and
//! merge), a mode phase's sign reference (`[[phase]] mode_reference = { file = … }`); in a joint control file, each
//! pattern's control file (`[[joint.pattern]] control`).

use super::syntax::{self, Document, Item, Span, Value};
use std::path::{Path, PathBuf};

/// A path a control file names: where its value is in the text, the table and key that hold it, and the path as
/// written.
#[derive(Clone, Debug, PartialEq)]
pub struct PathEntry {
    pub span: Span,
    /// "pattern", "pattern.background", "pattern.operation", "phase", "joint.pattern".
    pub table: String,
    /// "data", "file", "profile_files", "control", or "mode_reference.file".
    pub key: String,
    pub path: String,
}

/// The settings of the data operations that are paths (the operations' dictionary: unit "path").
fn operation_path_keys(op: &str) -> Vec<&'static str> {
    crate::operation_table::OPERATIONS.iter().find(|o| o.id == op).map_or(Vec::new(), |o| o.settings.iter().filter(|s| s.unit == "path").map(|s| s.key).collect())
}

fn string_at(it: &Item) -> Option<(String, Span)> {
    match &it.value {
        Value::String(s) => Some((s.clone(), it.span)),
        _ => None,
    }
}

/// Every path a control file's document names, in the order of the text.
pub fn paths(doc: &Document) -> Vec<PathEntry> {
    let mut out = Vec::new();
    for t in &doc.tables {
        let table = t.path.join(".");
        let keys: Vec<&str> = match table.as_str() {
            "pattern" => vec!["data"],
            "pattern.background" => vec!["file"],
            "joint.pattern" => vec!["control"],
            "pattern.operation" => match t.get("op").and_then(|e| string_at(&e.value)) {
                Some((op, _)) => operation_path_keys(&op),
                None => Vec::new(),
            },
            _ => Vec::new(),
        };
        for e in &t.entries {
            if keys.contains(&e.key.as_str()) {
                if let Some((path, span)) = string_at(&e.value) {
                    out.push(PathEntry { span, table: table.clone(), key: e.key.clone(), path });
                }
            }
            // The background's profiles (P23): an array of files.
            if table == "pattern.background" && e.key == "profile_files" {
                if let Value::Array(a) = &e.value.value {
                    for (path, span) in a.iter().filter_map(string_at) {
                        out.push(PathEntry { span, table: table.clone(), key: e.key.clone(), path });
                    }
                }
            }
            // A mode phase's sign reference: an inline table with its file.
            if table == "phase" && e.key == "mode_reference" {
                if let Value::Table(inner) = &e.value.value {
                    for x in inner.iter().filter(|x| x.key == "file") {
                        if let Some((path, span)) = string_at(&x.value) {
                            out.push(PathEntry { span, table: table.clone(), key: "mode_reference.file".into(), path });
                        }
                    }
                }
            }
        }
    }
    out.sort_by_key(|p| p.span.start);
    out
}

/// The file a path names, from the folder the control file is read in.
pub fn resolve(path: &str, folder: &Path) -> PathBuf {
    folder.join(path)
}

/// The text of a control file written in folder `from`, as it reads from folder `to`: each path it holds rewritten
/// relative to `to` (absolute when the two share nothing below the root, as on another volume), its comments, order
/// and formatting kept. Refused, the file named, when a path names a file that is not there, and when the text does
/// not parse.
pub fn relocate(text: &str, from: &Path, to: &Path) -> Result<String, String> {
    if !to.is_dir() {
        return Err(format!("{} is not a folder", to.display()));
    }
    let doc = syntax::parse(text);
    if doc.has_errors() {
        let d = doc.diagnostics.iter().find(|d| !d.warning).map_or(String::new(), |d| format!(" (line {}: {})", d.span.line, d.message));
        return Err(format!("the control file does not read{d}: its paths are not rewritten"));
    }
    let mut edits = Vec::new();
    for p in paths(&doc) {
        let file = resolve(&p.path, from);
        if !file.exists() {
            return Err(format!("[{}] {} = \"{}\": {} is not there", p.table, p.key, p.path, file.display()));
        }
        let new = super::relative(&file, to);
        edits.push((p.span, syntax::string(&new)));
    }
    let mut out = String::with_capacity(text.len() + 64);
    let mut pos = 0;
    for (span, s) in edits {
        out.push_str(&text[pos..span.start]);
        out.push_str(&s);
        pos = span.end;
    }
    out.push_str(&text[pos..]);
    Ok(out)
}
