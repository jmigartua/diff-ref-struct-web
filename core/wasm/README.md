# drs-wasm: the core in the browser

The web stream's W1 (I10, R34; L423, L424): the diff-ref-struct core compiled to WebAssembly, with the Plots module's
commands as the app serves them, so that the module runs in a browser unchanged behind its one `invoke(cmd, args)`.
No number is computed here: `src/host.rs` is a port of the app's `app/src-tauri/src/modules/plots.rs` (the sets by id,
the window's words mapped to the core's types), and the core computes as in the app and on the command line.

## Use

```js
import init, { dispatch, commands, version } from "./pkg/drs_wasm.js";
await init();                                   // or initSync({ module: bytes })
const set = await dispatch("plots_open_from", JSON.stringify({ files: [{ name: "tb_t02.dat", text }] }));
const map = await dispatch("plots_map", JSON.stringify({ id: set.id, scale: "linear", reference: null }));
```

`dispatch(cmd, argsJson)` returns a promise of the parsed result, rejected with the core's words. Each command is also
exported by its name, JSON in and JSON out. The names, arguments (camelCase, as Tauri takes them) and results are the
Tauri commands', except where a browser has no path:

| the app's command | here | what differs |
|---|---|---|
| `plots_open_files { paths, root }` | `plots_open_from { files, reader?, title?, root?, rebin?, xInterval? }` | files by `{ name, text, comment?, variable?, modified? }`; names relative to the first file's folder; `rebin` a step in x, every member rebinned on opening (W9); `xInterval: [a, b]` every member cut to it (PL42) |
| `plots_open_buffer { path, root }` | `plots_open_buffer_from { path?, text, files, root?, xInterval? }` | the buffer's text, its files by the path the buffer's folder gives them or by name |
| `plots_add { id, paths }` | `plots_add_from { id, files }` | files by name and text |
| `plots_open_set { path, root }` | `plots_open_set_from { name or path, text, files, root?, xInterval? }` | the `.drv`'s text, its members' files by the path joined to the folder above `Plots/` or by name |
| `plots_temperatures { id, path }` | `plots_temperatures_from { id, name, text }` | the temperatures file's text |
| `plots_save { … }` → the file's path | `plots_save { … }` and `plots_save_text { … }` → `{ relative, path, text }` | the file returned for the host to save, `path` `Plots/<name>.drv` (W5) |
| — | `plots_names_in { path, text, root? }` → `{ files }` | the files a `.buf` or a `.drv` names, joined as the app joins paths: the host reads them (W5) |
| — | `plots_propose { files, xInterval? }` → `{ proposals, classes, left_out }` | the series the files offer, each a set file to open as it is, with its `rebin` and the `reader` its members share (`plots::propose`, W2, W9) |
| `plots_rebin { id, step }` | the same | the set rebinned from its measured scans, kept the first time (no files to read again); a set opened rebinned refuses another step |
| `plots_derive { … }` → `DerivedView` | the same, with `text` | the derived file returned for the host to save |
| `plots_new { root }` | `plots_new {}` | no project |
| — | `describe { name, text }` → `{ kind, id, name, role, source }` | the core's `formats::project::describe` |

The same as the app's: `plots_patterns` (W9: `[{ member, x, y, sigma }]`, the members shown), `plots_get`, `plots_close`, `plots_proposed_step`, `plots_set_variable`, `plots_set_other`, `plots_displace`,
`plots_extent` (PL41), `plots_map`, `plots_map_region`, `plots_range_warning`, `plots_project`, `plots_difference`, `plots_reverse`,
`plots_sort`, `plots_read_variable`, `plots_windows`, `plots_view_rule` (a control file's λ refused: type it). Not
here: Data's operations (they live in the project's `.drs/data/`), FullProf's `.seq`, the results series, the CSV
export, `plots_reload`, `plots_write_text`, `plots_figure_widths`. `tests/commands.rs` (`cargo test`) reads the app's
`#[tauri::command]`s and fails, naming it, on a command here that is not the app's or the app's that is not here
(a command ahead of the app's binding accepted when declared),
unless declared with its reason.

## Build and test

You need Rust with the `wasm32-unknown-unknown` target, `wasm-pack` (0.15) and node (22). From `core/wasm`:

```sh
npm test        # wasm-pack build --target web --release, then node --test tests/*.test.mjs
```

The test loads the package in node, gives it the Tb₂Pd₂Sn files by name and text, and checks the map equal to
`drs plots map --json`'s to the bit (the command line built from `core/` when it is not); it prints the `.wasm`'s size.
Text files only: a binary format (Bruker RAW, BRML) passed as text is not read as its bytes.
