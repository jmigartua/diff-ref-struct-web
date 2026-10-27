//! The file dictionary in use (R27, L091): the name of each kind of file of a project started in the app, filled from
//! its pattern (`path`), and the kind of a file from its place and name in the project (`kind_of`). The dictionary is
//! `ecosystem/data/files.yml`, generated into `file_table`. A project that comes from FullProf keeps FullProf's names,
//! which `project::describe` names.

pub use super::file_table::{FileKind, FILES};

/// The placeholders a pattern may use ({stem}: the control file's name without extension, or the data file's before a
/// control file exists, L116; {name}: a publication set's, the app's request 37).
pub const PLACEHOLDERS: [&str; 12] = ["project", "control", "stem", "phase", "data", "original", "instrument", "time", "n", "set", "operation", "name"];

/// The kind with this id.
pub fn kind(id: &str) -> Option<&'static FileKind> {
    FILES.iter().find(|k| k.id == id)
}

/// The placeholders of a pattern, in order.
pub fn placeholders(pattern: &str) -> Vec<&str> {
    pattern.split('{').skip(1).filter_map(|s| s.split_once('}').map(|(p, _)| p)).collect()
}

/// The path, relative to the project's folder, of a file of kind `id` with its placeholders filled from `values`
/// (name → value). Refused for an unknown kind, a placeholder without a value, or a value holding a path separator.
pub fn path(id: &str, values: &[(&str, &str)]) -> Result<String, String> {
    let k = kind(id).ok_or(format!("{id}: not a kind of the file dictionary"))?;
    let mut name = k.pattern.to_string();
    for p in placeholders(k.pattern) {
        let v = values.iter().find(|(n, _)| *n == p).map(|(_, v)| *v).ok_or(format!("{id}: the name needs {{{p}}}"))?;
        if v.contains(['/', '\\']) {
            return Err(format!("{id}: {{{p}}} = \"{v}\" is not a name"));
        }
        name = name.replacen(&format!("{{{p}}}"), v, 1);
    }
    Ok(if k.folder.is_empty() { name } else { format!("{}/{name}", k.folder) })
}

/// A time as the dictionary writes it in names, YYYY-MM-DD HH.MM.SS (the versions of L080), from its parts.
pub fn time_text(year: i32, month: u32, day: u32, hour: u32, minute: u32, second: u32) -> String {
    format!("{year:04}-{month:02}-{day:02} {hour:02}.{minute:02}.{second:02}")
}

/// The time now, as `time_text` writes it, in UTC (the core keeps no time-zone table; the app passes its local time).
pub fn now_text() -> String {
    let secs = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map_or(0, |d| d.as_secs()) as i64;
    let (days, rem) = (secs.div_euclid(86400), secs.rem_euclid(86400));
    // Days since 1970-01-01 to a civil date (H. Hinnant's algorithm).
    let z = days + 719468;
    let era = z.div_euclid(146097);
    let doe = z - era * 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let m = if mp < 10 { mp + 3 } else { mp - 9 } as u32;
    let y = (yoe + era * 400 + i64::from(m <= 2)) as i32;
    time_text(y, m, d, (rem / 3600) as u32, (rem % 3600 / 60) as u32, (rem % 60) as u32)
}

/// Whether `name` fits `pattern`: its literal parts in order, a placeholder standing for any non-empty text.
fn fits(pattern: &str, name: &str) -> bool {
    let mut parts: Vec<&str> = Vec::new();
    let mut rest = pattern;
    while let Some(i) = rest.find('{') {
        parts.push(&rest[..i]);
        rest = rest[i..].split_once('}').map_or("", |(_, r)| r);
    }
    parts.push(rest);
    let (first, last) = (parts[0], parts[parts.len() - 1]);
    if parts.len() == 1 {
        return name == first;
    }
    if !name.starts_with(first) || !name.ends_with(last) || name.len() < first.len() + last.len() {
        return false;
    }
    let mut s = &name[first.len()..name.len() - last.len()];
    for p in &parts[1..parts.len() - 1] {
        // Each placeholder takes at least one character before the next literal part.
        match s.get(1..).and_then(|t| t.find(p)) {
            Some(j) => s = &s[j + 1 + p.len()..],
            None => return false,
        }
    }
    !s.is_empty()
}

/// The kind of a file from its path relative to the project's folder: the kinds whose folder holds it and whose
/// pattern it fits, the most specific (the longest literal pattern) first.
pub fn kind_of(rel: &str) -> Option<&'static FileKind> {
    let rel = rel.replace('\\', "/");
    let (folder, name) = rel.rsplit_once('/').unwrap_or(("", &rel));
    // The most literal text first: "{data} processed.xys" before "{original}", which fits everything in Data/.
    let literal = |k: &FileKind| placeholders(k.pattern).iter().fold(k.pattern.len(), |n, p| n - p.len() - 2);
    FILES.iter().filter(|k| k.folder == folder && fits(k.pattern, name)).max_by_key(|k| literal(k))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn patterns_filled_and_recognised() {
        let p = path("version", &[("control", "y2o3_structure"), ("time", &time_text(2026, 10, 10, 16, 42, 7))]).unwrap();
        assert_eq!(p, "Versions/y2o3_structure 2026-10-10 16.42.07.drc");
        assert_eq!(kind_of(&p).unwrap().id, "version");
        assert_eq!(kind_of("Versions/y2o3_structure refinement log.md").unwrap().id, "refinement_log");
        assert_eq!(kind_of("Exports/FullProf/y2o3_structure1.hkl").unwrap().id, "fullprof_reflections_export");
        assert_eq!(kind_of("Data/y2o3.dat").unwrap().id, "data");
        assert_eq!(kind_of("Data/y2o3 processed.xys").unwrap().id, "processed_data");
        assert!(path("version", &[("control", "a/b"), ("time", "t")]).is_err());
        assert!(path("structure", &[]).is_err());
    }
}
