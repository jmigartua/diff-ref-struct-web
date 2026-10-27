//! The package's commands against the app's (W1, PL41; L426): every Tauri command of
//! `app/src-tauri/src/modules/plots.rs` is served here by its own name, or by its text form (`_from`, W5) when the app's
//! takes paths, or is listed below as not in the browser with why. The test fails naming a command of the app that is
//! none of the three (a command added to the app and not ported, as `plots_extent` after Plots' step 20), and a command
//! here that the app does not have and that is not one of the browser's own.

use std::collections::BTreeSet;

/// The app's commands that take paths, and the commands here that take the files' names and text instead.
const TEXT_FORMS: [(&str, &str); 5] = [
    ("plots_open_files", "plots_open_from"),
    ("plots_open_buffer", "plots_open_buffer_from"),
    ("plots_add", "plots_add_from"),
    ("plots_open_set", "plots_open_set_from"),
    ("plots_temperatures", "plots_temperatures_from"),
];

/// The app's commands not in the browser, and why.
const NOT_HERE: [(&str, &str); 10] = [
    ("plots_export_csv", "writes a CSV to a path; no text form asked yet"),
    ("plots_reload", "reads the members' files again from disk"),
    ("plots_write_text", "the host writes (downloads) the text itself"),
    ("plots_figure_widths", "Publish's presets: drs-publish is not in the package"),
    ("plots_seq_open", "FullProf's .seq by its path"),
    ("plots_seq_column", "a .seq opened by plots_seq_open"),
    ("plots_seq_close", "a .seq opened by plots_seq_open"),
    ("plots_results_open", "results files by their paths"),
    ("plots_results_set_variable", "a series opened by plots_results_open"),
    ("plots_results_close", "a series opened by plots_results_open"),
];

/// Commands served here that the app is to gain, with the request that asks for them: the app's binding is the Plots
/// stream's to write (accepted here whether or not the app has it yet).
const AHEAD: [(&str, &str); 1] = [("plots_patterns", "W9 (a): the members' scans out of every SetView, asked for by the views that draw them; the app's #[tauri::command] in app/src-tauri/src/modules/plots.rs is Plots'")];

/// The browser's own commands, which the app has no reason to have.
const OWN: [&str; 4] = ["plots_names_in", "plots_save_text", "plots_propose", "describe"];

fn app_commands() -> BTreeSet<String> {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../app/src-tauri/src/modules/plots.rs");
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let lines: Vec<&str> = text.lines().collect();
    let mut out = BTreeSet::new();
    for (i, l) in lines.iter().enumerate() {
        if l.trim() == "#[tauri::command]" {
            // the next line that declares a fn (an #[allow] may come between)
            let decl = lines[i + 1..].iter().find(|x| x.contains("fn ")).expect("a command's fn");
            let name = decl.split("fn ").nth(1).unwrap().split(['(', '<']).next().unwrap().trim();
            out.insert(name.to_string());
        }
    }
    out
}

#[test]
fn the_packages_commands_are_the_apps() {
    let app = app_commands();
    assert!(app.len() >= 30, "{} commands read from the app's modules/plots.rs: {app:?}", app.len());
    let here: BTreeSet<String> = drs_wasm::COMMANDS.iter().map(|c| c.to_string()).chain(std::iter::once("describe".to_string())).collect();
    let mut wrong: Vec<String> = Vec::new();
    for c in &app {
        let served = here.contains(c) || TEXT_FORMS.iter().any(|(a, w)| a == c && here.contains(*w)) || NOT_HERE.iter().any(|(a, _)| a == c);
        if !served {
            wrong.push(format!("the app's {c} is neither served by the package, nor given a text form, nor listed as not in the browser"));
        }
    }
    for c in &here {
        let known = app.contains(c) || TEXT_FORMS.iter().any(|(_, w)| w == c) || OWN.contains(&c.as_str()) || AHEAD.iter().any(|(a, _)| a == c);
        if !known {
            wrong.push(format!("the package's {c} is not a command of the app, nor a text form of one, nor one of the browser's own"));
        }
    }
    for a in TEXT_FORMS.iter().map(|x| x.0).chain(NOT_HERE.iter().map(|x| x.0)) {
        if !app.contains(a) {
            wrong.push(format!("{a} is listed here and the app no longer has it"));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
    // Each command the dispatcher names is a command call() answers (not "no command").
    for c in &here {
        let e = drs_wasm::call(c, "{}").err().unwrap_or_default();
        assert!(!e.starts_with("no command"), "{c}: {e}");
    }
}
