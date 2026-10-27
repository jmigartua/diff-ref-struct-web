//! Completion for the control file's editor (P49, `control::complete`), on Y₂O₃'s control file (its `.pcr`
//! converted) and on tables added to it.
//!
//! - **Keys:** at an empty line of every table the reader knows, the keys offered are the table's keys less those it
//!   holds (and an operation's settings in `[[pattern.operation]]`); each, written there as `key = 0`, is one the
//!   reader accepts (no "is not a key of").
//! - **Tables:** at a `[`, the tables offered are the reader's, less the single tables already written; each,
//!   written as a header, is one the reader accepts.
//! - **Addresses:** inside `[refine] parameters` (a line of the array) and a link's `follows`, every address offered
//!   is one `refinable::field_of` reads back, with the dictionary's name and unit; a prefix narrows them.
//! - **Values:** the closed choices, each accepted by the reader; the species; the space groups.

use drs::control::complete::{complete, CompletionKind};
use drs::control::schema::{read_text, TABLES};
use std::path::Path;

mod common;
use common::latin1;

fn y2o3() -> String {
    let d = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs/Y2O3");
    let p = drs::pcr::read(&latin1(&d.join("y2o3_structure.pcr"))).unwrap();
    drs::control::Control::from_pcr(&p, &d, "y2o3_structure", &d).unwrap().to_text()
}

/// The text with `insert` lines put after line index `at` (0-based), and the 1-based number of the first.
fn with(text: &str, at: usize, insert: &[&str]) -> (String, usize) {
    let mut lines: Vec<&str> = text.split('\n').collect();
    for (k, l) in insert.iter().enumerate() {
        lines.insert(at + 1 + k, l);
    }
    (lines.join("\n"), at + 2)
}

/// A header line: `[name]` or `[[name]]` (an array value's line starts with `["` or a number).
fn is_head(l: &str) -> bool {
    let t = l.trim_start();
    t.starts_with("[[") || t.strip_prefix('[').is_some_and(|r| r.starts_with(|c: char| c.is_ascii_alphabetic()))
}

fn not_a_key(text: &str) -> Vec<String> {
    read_text(text).1.into_iter().filter(|d| d.message.contains("is not a key of") || d.message.contains("is not a table of")).map(|d| d.to_string()).collect()
}

#[test]
fn keys_offered_are_read() {
    let base = y2o3();
    assert!(not_a_key(&base).is_empty());
    let lines: Vec<&str> = base.split('\n').collect();
    let op = drs::operation_table::OPERATIONS.iter().find(|o| !o.settings.is_empty()).unwrap();
    let operation_line = format!("op = \"{}\"", op.id);
    let mut total = 0;
    for (path, array, keys) in TABLES {
        let name = path.join(".");
        let head = if *array { format!("[[{name}]]") } else { format!("[{name}]") };
        // The table's header in the file, or the table added at the end; an operation's table names one.
        let (text, line) = match lines.iter().position(|l| l.split('#').next().unwrap().trim() == head) {
            Some(i) => with(&base, i, &[""]),
            None if path.is_empty() => with(&base, 0, &[""]),
            None => {
                let mut add = vec![head.as_str()];
                if *path == ["pattern", "operation"] {
                    add.push(&operation_line);
                }
                add.push("");
                let (t, l) = with(&base, lines.len() - 1, &add);
                (t, l + add.len() - 1)
            }
        };
        let offered: Vec<String> = complete(&text, line, 1).into_iter().inspect(|c| assert_eq!(c.kind, CompletionKind::Key)).map(|c| c.label).collect();
        // The table's keys less those it holds.
        let tl: Vec<&str> = text.split('\n').collect();
        let start = (0..line - 1).rev().find(|&i| is_head(tl[i])).map_or(0, |i| i + 1);
        let stop = (line..tl.len()).find(|&i| is_head(tl[i])).unwrap_or(tl.len());
        let held: Vec<&str> = tl[start..stop].iter().filter_map(|l| l.split_once('=').map(|(k, _)| k.trim())).collect();
        let mut expected: Vec<String> = keys.iter().filter(|k| !held.contains(k)).map(|k| k.to_string()).collect();
        if *path == ["pattern", "operation"] {
            expected.extend(op.settings.iter().map(|s| s.key.to_string()));
        }
        assert_eq!(offered, expected, "{head}");
        for k in &offered {
            let (t, _) = with(&text, line - 2, &[&format!("{k} = 0")]);
            assert!(not_a_key(&t).is_empty(), "{head} {k}: {:?}", not_a_key(&t));
        }
        // A typed beginning narrows them.
        if let Some(k) = offered.first() {
            let (t, l) = with(&text, line - 2, &[&k[..2]]);
            assert!(complete(&t, l, 3).iter().all(|c| c.label.starts_with(&k[..2])));
        }
        println!("{head:24} {} keys offered", offered.len());
        total += offered.len();
    }
    assert!(total > 50, "{total}");
}

/// Each table alone (under its parent), so that every key it takes is offered and read.
#[test]
fn every_key_offered_is_read() {
    let mut total = 0;
    for (path, array, keys) in TABLES {
        let name = path.join(".");
        let mut add: Vec<String> = vec!["format = \"diff-ref-struct 1\"".into()];
        if path.len() > 1 {
            add.push(if path[0] == "phase" { "[[phase]]".into() } else { format!("[{}]", path[0]) });
        }
        if !path.is_empty() {
            add.push(if *array { format!("[[{name}]]") } else { format!("[{name}]") });
        }
        if *path == ["pattern", "operation"] {
            add.push("op = \"mask\"".into());
        }
        let text = add.join("\n") + "\n";
        let line = add.len() + 1;
        let offered: Vec<String> = complete(&text, line, 1).into_iter().map(|c| c.label).collect();
        let mut expected: Vec<String> = keys.iter().filter(|k| !(path.is_empty() && **k == "format") && !(*path == ["pattern", "operation"] && **k == "op")).map(|k| k.to_string()).collect();
        if *path == ["pattern", "operation"] {
            expected.extend(drs::operation_table::OPERATIONS[0].settings.iter().map(|s| s.key.to_string()));
        }
        assert_eq!(offered, expected, "[{name}]");
        for k in &offered {
            let t = format!("{text}{k} = 0\n");
            assert!(not_a_key(&t).is_empty(), "[{name}] {k}: {:?}", not_a_key(&t));
        }
        total += offered.len();
    }
    println!("{total} keys offered over {} tables, each read", TABLES.len());
    assert!(total > 140, "{total}");
}

#[test]
fn tables_offered_are_read() {
    let base = y2o3();
    let n = base.split('\n').count();
    let (text, line) = with(&base, n - 1, &["["]);
    let offered = complete(&text, line, 2);
    assert!(offered.iter().all(|c| c.kind == CompletionKind::Table));
    let labels: Vec<&str> = offered.iter().map(|c| c.label.as_str()).collect();
    println!("{labels:?}");
    // The single tables the file holds are not offered again; the others and every array table are.
    for held in ["[pattern]", "[pattern.background]", "[pattern.fullprof]", "[refine]", "[convention]"] {
        assert!(!labels.contains(&held), "{held}");
    }
    assert_eq!(labels.len(), TABLES.len() - 1 - 5);
    for c in &offered {
        let mut add = vec![c.insert.as_str()];
        if c.label == "[[pattern.operation]]" {
            add.push("op = \"mask\"");
        }
        let (t, _) = with(&base, n - 1, &add);
        assert!(not_a_key(&t).is_empty(), "{}: {:?}", c.label, not_a_key(&t));
    }
    // A typed beginning narrows them; a line of an array value is not a header.
    assert_eq!(complete(&text.replace("\n[", "\n[[phase.a"), line, 10).iter().map(|c| c.label.as_str()).collect::<Vec<_>>(), ["[[phase.atom]]"]);
}

#[test]
fn addresses_offered_are_read() {
    let base = y2o3();
    let c = drs::control::Control::from_text(&base).unwrap();
    let lines: Vec<&str> = base.split('\n').collect();
    // A line inside [refine] parameters = [ … ], after its "[\"".
    let i = lines.iter().position(|l| l.trim_start().starts_with("[\"pattern1.zero_shift\"")).unwrap();
    let col = lines[i].find('"').unwrap() + 2;
    let offered = complete(&base, i + 1, col);
    let settable: Vec<String> = c.model.fields().into_iter().filter(|f| f.fixed_by.is_none()).map(|f| f.address).collect();
    println!("{} addresses offered, {} settable fields", offered.len(), settable.len());
    assert_eq!(offered.iter().map(|c| c.label.clone()).collect::<Vec<_>>(), settable);
    for o in &offered {
        assert_eq!(o.kind, CompletionKind::Address);
        drs::refinable::field_of(&o.label, &c.model).unwrap_or_else(|e| panic!("{e}"));
    }
    let zero = offered.iter().find(|o| o.label == "pattern1.zero_shift").unwrap();
    assert_eq!(zero.detail, "Zero shift (deg 2θ)");
    assert!(!zero.documentation.is_empty());
    // A prefix narrows them, in the array and in a link.
    let (t, l) = with(&base, i - 1, &["  [\"phase1.cell_"]);
    let cell: Vec<String> = complete(&t, l, 17).into_iter().map(|c| c.label).collect();
    assert_eq!(cell, ["phase1.cell_a", "phase1.cell_b", "phase1.cell_c", "phase1.cell_alpha", "phase1.cell_beta", "phase1.cell_gamma"].iter().filter(|a| settable.contains(&a.to_string())).map(|a| a.to_string()).collect::<Vec<_>>());
    assert!(!cell.is_empty());
    let n = lines.len();
    // A link just begun: the rest of the file reads without its table.
    let (t, l) = with(&base, n - 1, &["[[refine.link]]", "follows = \"phase1.Y:2.atom_"]);
    let link: Vec<String> = complete(&t, l + 1, 32).into_iter().map(|c| c.label).collect();
    println!("link: {link:?}");
    assert_eq!(link, ["phase1.Y:2.atom_position.x", "phase1.Y:2.atom_biso", "phase1.Y:2.atom_occupancy"]);
}

/// The settings whose offered symbol reads back as the setting, as measured in L260.
const SG_READ_BACK: usize = 530;

#[test]
fn values_offered() {
    let base = y2o3();
    let lines: Vec<&str> = base.split('\n').collect();
    // Each closed choice the editor offers for the geometry, the weighting and the profile is one the reader takes.
    for key in ["geometry", "weighting", "profile", "radiation", "preferred_orientation_function"] {
        let i = lines.iter().position(|l| l.starts_with(&format!("{key} ="))).unwrap();
        let col = key.len() + 4;
        let offered = complete(&base, i + 1, col);
        assert!(offered.len() >= 2, "{key}: {offered:?}");
        for o in &offered {
            let mut l2 = lines.clone();
            let line = format!("{key} = {}", o.insert);
            l2[i] = &line;
            // The value is read (a choice may leave other keys of the file without use: not the value's fault).
            let d: Vec<String> = read_text(&l2.join("\n")).1.into_iter().filter(|d| d.span.line == i + 1).map(|d| d.to_string()).collect();
            assert!(d.is_empty(), "{key} = {}: {d:?}", o.insert);
        }
        println!("{key}: {:?}", offered.iter().map(|o| o.label.as_str()).collect::<Vec<_>>());
    }
    // The species, from the scattering tables.
    let i = lines.iter().position(|l| l.starts_with("species =")).unwrap();
    let all: Vec<String> = complete(&base, i + 1, 12).into_iter().map(|c| c.label).collect();
    for s in ["Y", "O", "Y+3", "O-2", "Fe", "Fe+3", "D"] {
        assert!(all.contains(&s.to_string()), "{s}");
    }
    let (t, l) = with(&base, i - 1, &["species = \"Fe+"]);
    assert!(complete(&t, l, 15).iter().all(|c| c.label.starts_with("Fe+")));
    // The space groups: the 530 settings, each by a symbol that reads back as that setting when the Hall symbol is
    // not given (the short symbol, else the full one, else the short one with its choice); "21" matches "2_1".
    let i = lines.iter().position(|l| l.starts_with("space_group =")).unwrap();
    let (t, l) = with(&base, i - 1, &["space_group = \""]);
    let sg = complete(&t, l, 16);
    assert_eq!(sg.len(), 530);
    let back = sg.iter().filter(|c| drs::spacegroup::SpaceGroup::from_fullprof_symbol(&c.label).is_ok_and(|g| c.documentation == format!("Hall symbol {}", g.hall))).count();
    println!("space groups: {} offered, {back} read back as their setting", sg.len());
    assert_eq!(back, SG_READ_BACK);
    let (t, l) = with(&base, i - 1, &["space_group = \"P 21/"]);
    let sg = complete(&t, l, 21);
    assert!(sg.iter().any(|c| c.label == "P 2_1/c" && c.detail == "No. 14, setting b1"), "{sg:?}");
}
