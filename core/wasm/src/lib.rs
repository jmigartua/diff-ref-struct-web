//! The diff-ref-struct core in the browser (W1, I10, R34; L424): the Plots module's commands compiled to WebAssembly,
//! with the app's Tauri commands' names, arguments and JSON (`host`, a port of `app/src-tauri/src/modules/plots.rs`),
//! so that the module's `invoke(cmd, args)` is served by `dispatch(cmd, argsJson)`. The core computes every number as
//! in the app and the command line; this crate only carries JSON in and out.
//!
//! Each command is exported by its name, taking its arguments as JSON (the keys as the module passes them, camelCase as
//! Tauri takes them) and returning its result as JSON; `dispatch` takes the name and the arguments and returns a
//! promise of the parsed result, rejected with the core's words.

use serde::de::DeserializeOwned;
use serde::Deserialize;
use wasm_bindgen::prelude::*;

pub mod host;

use drs::plots::{OpenedFile, Range, Window};
use host::{NormChoice, RegionView, Settings, XChoice};

thread_local! {
    static SETS: host::Sets = host::Sets::default();
}

fn args<T: DeserializeOwned>(json: &str) -> Result<T, String> {
    serde_json::from_str(if json.trim().is_empty() { "{}" } else { json }).map_err(|e| format!("the command's arguments: {e}"))
}

fn out<T: serde::Serialize>(v: T) -> Result<String, String> {
    serde_json::to_string(&v).map_err(|e| e.to_string())
}

#[derive(Deserialize)]
struct Id {
    id: u32,
}

#[derive(Deserialize)]
struct Files {
    files: Vec<OpenedFile>,
    #[serde(default)]
    reader: Option<String>,
    #[serde(default)]
    title: Option<String>,
    #[serde(default)]
    rebin: Option<f64>,
}

#[derive(Deserialize)]
struct SetOther {
    id: u32,
    member: usize,
    name: String,
    value: f64,
    note: String,
}

#[derive(Deserialize)]
struct Patterns {
    id: u32,
    #[serde(default)]
    members: Option<Vec<usize>>,
}

#[derive(Deserialize)]
struct NamesIn {
    path: String,
    text: String,
}

#[derive(Deserialize)]
struct Extent {
    id: u32,
    scale: String,
    #[serde(default)]
    members: Option<Vec<usize>>,
}

#[derive(Deserialize)]
struct Propose {
    files: Vec<OpenedFile>,
}

#[derive(Deserialize)]
struct IdFiles {
    id: u32,
    files: Vec<OpenedFile>,
}

/// A buffer by its text and, as W5 sends it, its path (`root` taken and not needed: no Data operations here).
#[derive(Deserialize)]
struct BufferFrom {
    #[serde(default)]
    path: Option<String>,
    text: String,
    files: Vec<OpenedFile>,
}

/// A set file by its name (`name`, or W5's `path`) and text.
#[derive(Deserialize)]
struct SetFrom {
    #[serde(alias = "path")]
    name: String,
    text: String,
    files: Vec<OpenedFile>,
}

#[derive(Deserialize)]
struct Save {
    id: u32,
    root: String,
    name: String,
    settings: Settings,
}

#[derive(Deserialize)]
struct Scale {
    id: u32,
    scale: String,
}

#[derive(Deserialize)]
struct TextOf {
    id: u32,
    name: String,
    text: String,
}

#[derive(Deserialize)]
struct SetVariable {
    id: u32,
    member: usize,
    value: f64,
    note: String,
}

#[derive(Deserialize)]
struct Displace {
    id: u32,
    rule: String,
    dx: f64,
    dy: f64,
    scale: String,
}

#[derive(Deserialize)]
struct MapArgs {
    id: u32,
    scale: String,
    reference: Option<usize>,
}

#[derive(Deserialize)]
struct RegionArgs {
    id: u32,
    scale: String,
    reference: Option<usize>,
    region: RegionView,
}

#[derive(Deserialize)]
struct RangeArgs {
    id: u32,
    scale: String,
    reference: Option<usize>,
    range: Range,
}

#[derive(Deserialize)]
struct Derive {
    id: u32,
    members: Vec<usize>,
    kind: String,
    #[serde(default)]
    root: String,
    name: String,
}

#[derive(Deserialize)]
struct Project {
    id: u32,
    rule: String,
    dx: f64,
    dy: f64,
    scale: String,
    floor: Option<f64>,
    order: Option<Vec<usize>>,
}

#[derive(Deserialize)]
struct Difference {
    id: u32,
    a: usize,
    b: usize,
}

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct ReadVariable {
    id: u32,
    source: String,
    replace_typed: bool,
}

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct Windows {
    id: u32,
    windows: Vec<Window>,
    background: String,
    log_lambda: Option<f64>,
    at: Vec<f64>,
}

#[derive(Deserialize)]
struct ViewRule {
    id: u32,
    norm: Option<NormChoice>,
    x: Option<XChoice>,
}

#[derive(Deserialize)]
struct Describe {
    name: String,
    text: String,
}

/// What a file is (the core's `formats::project::describe`, by its name and content): `kind` "project", "data" or
/// "unknown", the kind's or the format's `id` and `name`, a project file's `role` and `source`.
fn describe_file(name: &str, text: &str) -> serde_json::Value {
    use drs::formats::project::{describe, Description, Role};
    match describe(&drs::formats::Probe::new(name, text.as_bytes())) {
        Description::Project(k) => serde_json::json!({ "kind": "project", "id": k.id, "name": k.name, "role": match k.role { Role::Control => "control", Role::Input => "input", Role::Output => "output", Role::Other => "other", Role::Unknown => "unknown" }, "source": k.source }),
        Description::Data { id, name } => serde_json::json!({ "kind": "data", "id": id, "name": name, "role": null, "source": null }),
        Description::Unknown => serde_json::json!({ "kind": "unknown", "id": "unknown", "name": null, "role": null, "source": null }),
    }
}

/// The commands, by their names: the app's (`plots_*`), with `_from` for those that take files by name and content
/// in place of paths, and `describe`.
pub const COMMANDS: [&str; 29] = [
    "plots_open_from",
    "plots_names_in",
    "plots_extent",
    "plots_propose",
    "plots_patterns",
    "plots_save_text",
    "plots_new",
    "plots_open_buffer_from",
    "plots_add_from",
    "plots_get",
    "plots_close",
    "plots_save",
    "plots_open_set_from",
    "plots_proposed_step",
    "plots_temperatures_from",
    "plots_set_variable",
    "plots_set_other",
    "plots_displace",
    "plots_map",
    "plots_map_region",
    "plots_range_warning",
    "plots_derive",
    "plots_project",
    "plots_difference",
    "plots_reverse",
    "plots_sort",
    "plots_read_variable",
    "plots_windows",
    "plots_view_rule",
];

/// One command by its name, its arguments and its result as JSON.
pub fn call(cmd: &str, json: &str) -> Result<String, String> {
    SETS.with(|s| match cmd {
        "plots_open_from" => {
            let a: Files = args(json)?;
            out(host::open_from(s, a.files, a.reader.as_deref(), a.title, a.rebin)?)
        }
        "plots_names_in" => {
            let a: NamesIn = args(json)?;
            out(serde_json::json!({ "files": host::names_in(&a.path, &a.text)? }))
        }
        "plots_extent" => {
            let a: Extent = args(json)?;
            out(host::extent(s, a.id, &a.scale, a.members)?)
        }
        "plots_propose" => out(host::propose(args::<Propose>(json)?.files)),
        "plots_set_other" => {
            let a: SetOther = args(json)?;
            out(host::set_other(s, a.id, a.member, &a.name, a.value, &a.note)?)
        }
        "plots_patterns" => {
            let a: Patterns = args(json)?;
            out(host::patterns(s, a.id, a.members)?)
        }
        "plots_new" => out(host::new_set(s)?),
        "plots_open_buffer_from" => {
            let a: BufferFrom = args(json)?;
            out(host::open_buffer_from(s, a.path.as_deref(), &a.text, a.files)?)
        }
        "plots_add_from" => {
            let a: IdFiles = args(json)?;
            out(host::add_from(s, a.id, a.files)?)
        }
        "plots_get" => out(host::get(s, args::<Id>(json)?.id)?),
        "plots_close" => {
            s.close(args::<Id>(json)?.id);
            out(())
        }
        // plots_save_text: W5's name for the same, { path, text } read from the same answer
        "plots_save" | "plots_save_text" => {
            let a: Save = args(json)?;
            out(host::save(s, a.id, &a.root, &a.name, a.settings)?)
        }
        "plots_open_set_from" => {
            let a: SetFrom = args(json)?;
            out(host::open_set_from(s, &a.name, &a.text, a.files)?)
        }
        "plots_proposed_step" => {
            let a: Scale = args(json)?;
            out(host::proposed(s, a.id, &a.scale)?)
        }
        "plots_temperatures_from" => {
            let a: TextOf = args(json)?;
            out(host::temperatures_from(s, a.id, &a.name, &a.text)?)
        }
        "plots_set_variable" => {
            let a: SetVariable = args(json)?;
            out(host::set_variable(s, a.id, a.member, a.value, a.note)?)
        }
        "plots_displace" => {
            let a: Displace = args(json)?;
            out(host::displace(s, a.id, &a.rule, a.dx, a.dy, &a.scale)?)
        }
        "plots_map" => {
            let a: MapArgs = args(json)?;
            out(host::map(s, a.id, &a.scale, a.reference)?)
        }
        "plots_map_region" => {
            let a: RegionArgs = args(json)?;
            out(host::map_region(s, a.id, &a.scale, a.reference, a.region)?)
        }
        "plots_range_warning" => {
            let a: RangeArgs = args(json)?;
            out(host::range_warning(s, a.id, &a.scale, a.reference, a.range)?)
        }
        "plots_derive" => {
            let a: Derive = args(json)?;
            out(host::derive(s, a.id, a.members, &a.kind, &a.root, &a.name)?)
        }
        "plots_project" => {
            let a: Project = args(json)?;
            out(host::project(s, a.id, &a.rule, a.dx, a.dy, &a.scale, a.floor, a.order)?)
        }
        "plots_difference" => {
            let a: Difference = args(json)?;
            out(host::difference(s, a.id, a.a, a.b)?)
        }
        "plots_reverse" => out(host::reverse(s, args::<Id>(json)?.id)?),
        "plots_sort" => out(host::sort(s, args::<Id>(json)?.id)?),
        "plots_read_variable" => {
            let a: ReadVariable = args(json)?;
            out(host::read_variable(s, a.id, &a.source, a.replace_typed)?)
        }
        "plots_windows" => {
            let a: Windows = args(json)?;
            out(host::windows(s, a.id, a.windows, &a.background, a.log_lambda, a.at)?)
        }
        "plots_view_rule" => {
            let a: ViewRule = args(json)?;
            out(host::set_view_rule(s, a.id, a.norm, a.x)?)
        }
        "describe" => {
            let a: Describe = args(json)?;
            out(describe_file(&a.name, &a.text))
        }
        o => Err(format!("no command \"{o}\" in the browser's core (drs-wasm {}): {}, describe", env!("CARGO_PKG_VERSION"), COMMANDS.join(", "))),
    })
}

/// `invoke(cmd, args)` as the app's shell serves it: the command by its name, its arguments as JSON; a promise of the
/// parsed result, rejected with the core's words.
#[wasm_bindgen]
pub fn dispatch(cmd: &str, args_json: &str) -> js_sys::Promise {
    match call(cmd, args_json) {
        Ok(j) => js_sys::Promise::resolve(&js_sys::JSON::parse(&j).unwrap_or_else(|_| JsValue::from_str(&j))),
        Err(e) => js_sys::Promise::reject(&JsValue::from_str(&e)),
    }
}

/// The core's version, as `drs --version` gives it.
#[wasm_bindgen]
pub fn version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}

/// The command names `dispatch` serves.
#[wasm_bindgen]
pub fn commands() -> Vec<String> {
    COMMANDS.iter().map(|c| c.to_string()).chain(std::iter::once("describe".to_string())).collect()
}

macro_rules! export {
    ($($name:ident),* $(,)?) => {
        $(
            #[doc = concat!("The command `", stringify!($name), "`: its arguments as JSON in, its result as JSON out.")]
            #[wasm_bindgen]
            pub fn $name(args_json: &str) -> Result<String, JsValue> {
                call(stringify!($name), args_json).map_err(|e| JsValue::from_str(&e))
            }
        )*
    };
}

export!(
    plots_open_from,
    plots_names_in,
    plots_extent,
    plots_propose,
    plots_patterns,
    plots_save_text,
    plots_new,
    plots_open_buffer_from,
    plots_add_from,
    plots_get,
    plots_close,
    plots_save,
    plots_open_set_from,
    plots_proposed_step,
    plots_temperatures_from,
    plots_set_variable,
    plots_set_other,
    plots_displace,
    plots_map,
    plots_map_region,
    plots_range_warning,
    plots_derive,
    plots_project,
    plots_difference,
    plots_reverse,
    plots_sort,
    plots_read_variable,
    plots_windows,
    plots_view_rule,
    describe,
);
