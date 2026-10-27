# web: requests to the core, the Plots module and the app's shell

What the web stream needs and does not edit itself (L423). Numbered W1, W2 …; numbers never change. Each says what,
why and the evidence; the merging session passes them on, and a request is closed by naming the version that provides
it.

## To the core

### W1. The core in the browser: a set from names and text, and a WebAssembly build (I10 step 1; filed by the merging session, L423)

- **What**: (a) `plots::open_from(files: &[OpenedFile { name, text, variable: Option<(f64, String)> }], reader, title)
  -> Result<PlotSet, String>`: the set from file names and contents, the readers' text functions (`data::read(text,
  ins)`, the registry's identification on the text and the name) with no path touched; `plots::open` becomes a wrapper
  that reads the paths and calls it. The same for a `.drv` text whose members' contents are given by name
  (`SetFile::open_from`). (b) A crate `core/wasm` (`drs-wasm`, `wasm-bindgen`), the plots subset exported for the
  browser host as the app's commands are (`plots_open_from`, `plots_map`, `plots_displace`, `plots_project`,
  `plots_windows`, `plots_derive`, `plots_map_region`, `plots_range_warning`, `plots_page_data`, `describe`, the set
  file read and written as text), each taking and returning the JSON the Tauri commands take and return, so that the
  module's `invoke(cmd, args)` is served by `dispatch(cmd, args)` with the same names and shapes. Built with
  `wasm-pack build --target web`; the core's own crate unchanged in what it computes (the `std::fs` uses stay behind
  the path wrappers). (c) Test: in node, the wasm package loaded, the Tb₂Pd₂Sn files given by name and text, the map
  equal to `drs plots map`'s to the bit (`core/wasm/tests/`), and a size line (the `.wasm` in kB) in the test's
  report.
- **Why**: the one architecture that keeps L139 in a browser (I10's assessment): the same crate, no second
  implementation; the files read in the browser and never uploaded.

### W2. The series recognised: `plots::propose(files)` (I10 step 3; filed by the merging session, L423)

- **What**: from the files' facts (the registry's: format, points, step, λ, title, temperature where the header or the
  name gives it with its source, the elapsed time from the name or the modification time, the exposure class by the
  esd's median, PL2's sources), the structure of a run and the sets it offers: one scan per temperature (a ramp: the
  temperatures distinct, ordered), holds (runs of equal temperature of three scans or more, against time since the
  run's first), the whole run against time, each class apart when two are found; each proposal a `SetFile` the host
  opens as it is, with a title that says what it is and how many scans, and the variables with their sources; the
  files left out named with why.
- **Why**: "the app offers the possible map plots from the files" (Josu, I10): what F07 and the timeline S01–S12 of
  the team's ESRF run were found by hand. Tested on the 101 names and facts of that run (the three sets of L021 found:
  25 by temperature, 100 against time, the 1270 °C hold of 39) and on Tb₂Pd₂Sn's buffer (one set, by temperature).

### W3. A second variable per member, for the map's right-hand axis (I10 step 0; filed by the merging session, L423)

- **What**: a member of a set may carry, beside the set's variable, named values with their sources
  (`variables = { T = { value = 1270, source = "…" } }` in the member's table of the `.drv`; the registry's and
  PL2's sources fill `T` and `t` when found); `plots::map` returns each row's other variables, so that the shell
  labels the right-hand axis with them; `propose` (W2) fills them. Absent by default; a `.drv` written so far reads
  as before.
- **Why**: Josu, I10: the hundred-scan map shows hours, should it be temperatures? Many scans share a temperature
  (the holds), so time is the axis and the temperature belongs on a second axis at the right, each row labelled
  (the team's ESRF project asked it from the start).

### W5. The path commands as text, and the package's three exports (I10 step 2; filed by the web stream, 25 October 2026)

- **What**: beside W1's list, in the package's `dispatch(cmd, args_json) -> result_json` (wasm-pack `--target web`,
  `--out-name drs_wasm`; exports `default` (init), `dispatch`, and `version` if it can):
  - `plots_open_from { files: [{ name, text }], reader, title, root }` → the SetView `plots_open_files` returns (W1);
  - `plots_add_from { id, files: [{ name, text }] }` → SetView, as `plots_add`;
  - `plots_names_in { path, text, root }` → `{ files: [path …] }`: the files a `.buf` or a `.drv` names, joined to the
    folder the core resolves them from (`set_base`, the buffer's folder) as the app joins real paths;
  - `plots_open_buffer_from { path, text, files: [{ name, text }], root }` → SetView, as `plots_open_buffer`;
  - `plots_open_set_from { path, text, files, root }` → `[SetView, Settings]`, as `plots_open_set` (a member the files
    lack named as missing, as the app names a file not on disk);
  - `plots_save_text { id, root, name, settings }` → `{ path: "Plots/<name>.drv", text }`: what `plots_save` would
    write, not written.
  Names are the host's virtual paths (`/local/<folder>/<path>`), strings the core joins and relativises as it does real
  ones. No Data operations (`.drs/data/`) in the browser: the members as measured.
- **Why**: the module asks with paths (`plots_open_files`, `plots_add`, `plots_open_set`, `plots_open_buffer`,
  `plots_save`); the host reads the files the core says it needs and hands their text over, so that the host never
  parses a set file or a buffer (L139). `app/ui/host-web/dispatch.js` sends exactly these; `web/tests/host.test.mjs`
  pins the shapes.

### W8. The ESRF names read: a temperature in °C and the elapsed time (I10 step 3; filed by the web stream, 25 October 2026)

- **What**: the variable from a name such as `scan_001_T0900C_t00000s_long.xye`: `T0900C` a temperature of
  900 °C, kept in °C or given in K with its conversion said in the source (as the headers' °C are,
  `core/src/formats/mod.rs:479`), and `t00000s` an elapsed time in s (W3's `t`); the `long`/`short` part left to W2's
  exposure classes. Tested on the 101 names of the team's ESRF run.
- **Why**: drs-wasm 0.79.0 opens the 101 files as one set in the browser (web/tests/package.test.mjs) with every
  member drawn apart, no `T` found: `temperature_in_name` (`core/src/plots.rs:266`) takes `_300K`, `T300`, `1p5K` and
  no `C` suffix; so neither the map by temperature nor W2's proposals (one scan per temperature, the 1270 °C hold of 39
  against time) can come from these names until the core reads them.

### W9. A large set in the browser: its patterns not in every SetView, its overlay drawn in reasonable time (I10 step 3; filed by the web stream, 25 October 2026; for the core and Plots)

- **What**: (a) the core: the members' x and y out of the SetView the set commands return (asked by the views that draw
  them, or a decimated overlay from the core); (b) Plots: the Overlay (and its linked zoom) of 101 members × 27 000
  points drawn without blocking the page for minutes.
- **Why**, measured on the ESRF `renamed/` folder (101 `.xye`, 112 MB), drs-wasm 0.79.0:
  - in node, `plots_open_from` takes 30 s and returns 71 MB of JSON; `plots_map` 2.6 s (log, 101 MB of JSON) and 1.1 s
    (linear, 57 MB); `describe` 0.18 s a file; `plots_difference` 0.07 s;
  - in Chrome, through the page: read in 1.4 s, described in 66 s, opened in 20 s more; the host's panels redraw the
    101 members in 124 ms; switching to the Overlay then blocks the page's thread for more than ten minutes (the
    renderer at 100 % and 3.4 GB; nothing logged to the console), and the Map came only about 16 minutes after the set
    opened. Tb₂Pd₂Sn (20 × 301 points) is instant. The time is the module's drawing, not the core's nor the host's;
    not profiled further.

### W10. A proposal opened as its set file: its `.drv` text from `plots_propose` (I10 step 4; filed by the web stream, 25 October 2026)

- **What**: each proposal of `plots_propose` with its set file as text (`SetFile::to_text`) beside the structure, or
  `plots_open_set_from` taking the structure (`set`); the host then opens a proposal as a set file (Open Plot Set),
  with the values it types and their sources.
- **Why**: drs-wasm 0.80.0 gives `set` as the `SetFile` structure, and no command opens one; the host opens a proposal
  by its members' files (`plots_open_from`), so what the proposal types is lost: on the ESRF run every time of the
  "against time" proposal (typed from the files' modification times, `elapsed_in_name` not taking `t00000s`, which
  sits inside a word). The page says so in its status line when it happens. The host writing `.drv` text itself
  would be a second writer of the core's format.
- **Measured on core 0.81.0** (26 October 2026, the package rebuilt at 0.81.0 and run in node on the 101 ESRF files,
  while writing `web/manual/`): `plots_propose` types the time from the names now (`t00000s`), so the proposals
  against time are right, but opening one from its files makes each member's **variable the temperature** in its name
  and the time only its second variable: the hold's 39 members open at T = 1270 °C, every map row with edges
  [1270, 1270], so the hold's map collapses to one line (the 100-scan proposal opens by T as well). The
  by-temperature proposal is unaffected. The set files typing the time (the project's `hold-1270.drv`,
  `all-scans.drv`) open right through `plots_open_set_from` (t in min and h, map rows apart). The manuals tell the
  team to use them meanwhile (`web/manual/user-manual.md` §7, §20; the team's use case §4).

### W11. The package's answers in the app's shapes, and a test of the shapes (filed by the web stream, 26 October 2026)

- **What**: (a) `plots_patterns` in drs-wasm answering each pattern as the app's bridge does, `{ member, x, y, sigma }`
  (`app/src-tauri/src/modules/plots.rs:190`), not `{ index, key, x, y, sigma }` (`core/wasm/src/host.rs:112`);
  (b) a test that compares, command by command, the JSON keys of the package's answers with the bridge's structures,
  beside the existing test of the command names (L427), so that a shape that drifts fails before it is published.
- **Why**: found while taking the guide's screenshots (`web/site/tools/shots.mjs`): on the public page (web a1c308d,
  core 0.81.0) the Overlay draws no curve for the Tb₂Pd₂Sn buffer, locally and published alike; the module files
  the patterns by `p.member` (`app/ui/modules/plots/index.js:613`), every pattern fell under one key. The Stack and
  the maps do not read `plots_patterns` and are unaffected. The web host now renames the field in the dispatcher
  (`SHAPES` in `app/ui/host-web/dispatch.js`, a rename, no number touched; tested), which retires when the core
  answers so itself.

### W13. The run's identifiers out of the core's sources and tests (filed by the web stream, 26 October 2026)

- **Met** 27 October 2026 (L439), done on main at Josu's word (the core stream merges main before its 0.82.0): the names
  `scan_<seq>_…`, the scan's number read before its temperature, the folder from `DRS_ESRF_RENAMED` only; `plots_esrf`
  2 of 2 and `plots_propose` 5 of 5, the package's tests 9 of 9, with the folder.

- **What**: the ESRF run's identifiers (its proposal number, its sample, its project's folder name) replaced by
  neutral ones in `core/src/plots.rs` (a comment), `core/tests/plots_esrf.rs`, `core/tests/plots_propose.rs` and
  `core/wasm/tests/map.test.mjs`: the 101 names kept in their structure (`scan_001_T0900C_t00000s_long.xye` …) so that
  the tests test the same rules, the private folder read from an environment variable only (as
  `web/tests/package.test.mjs` now does with `DRS_ESRF_RENAMED`).
- **Why**: the source of the page goes public under the Apache License 2.0 (Josu, 26 October 2026); the public
  selection (`web/tools/source.mjs`) is guarded against the team's unpublished data, and these four files are the only
  ones it stops on (`python3 web/site/tools/release.py`: "source" fails). The team's data stay unpublished.

### W14. To the merging session: the log by tool, the decisions of 26 October, Pages (filed by the web stream)

- **What**: (a) a field `tools: [web]` on the log entries about the page, and a `public:` one-line summary on an entry
  that names a team's data, so that the guide's Log (`web/site/tools/build.py`, which today finds the entries by their
  words and leaves out six that name the run) shows every step; (b) the decisions taken by Josu on 26 October 2026
  logged: the page open source (code Apache-2.0, documents CC BY 4.0, CITATION.cff, a DOI per release), the source
  published as a subset by an Action, requests as the public repository's issues with an e-mail path, the vocabulary
  appendix named Conventions, semantic versions from diff-ref-struct Web 1.0.0, the example computed by the core, the
  AI-assistance statement (its wording in `web/public/NOTICE`), English only for now; (c) after merging, the public
  repository's Pages switched to the branch `gh-pages` (Settings ▸ Pages), then `publish-source.yml` run by hand once
  W13 is met, then Zenodo linked to the repository for the DOIs.
- **Why**: the web stream edits only `web/`, `app/ui/host-web/` and its Actions (L423).

## To the Plots module

### W4. Step 0's controls (I10; filed by the merging session, L423)

- **What**: (a) a toggle of the y interval, shared by the whole set (the set's extent, the core's) or adapted to the
  member shown, for the row strip and the curves (Overlay, Stack); (b) the right-hand axis of the Map and the
  Difference Map labelled with a member's second variable when the core gives one (W3); (c) the members table richer:
  an icon by the file's kind (the core's `describe`), a checkbox per member and all / none / by exposure class / by a
  range of the variable, sort by variable, time or name, the member's facts on hover; dim and solo as they are.
- **Why**: Josu's asks on the present page that need no web host (I10 step 0); the web host takes the module as it is,
  so each lands in both.

### W7. A host's way to open a set, and Save's words from the host (I10 step 2; filed by the web stream, 25 October 2026)

- **What**: (a) `plots.adopt(set, title, path?, settings?)` exported (the module's own `adopt`), so that a host that
  opened a set itself shows it in a tab; today the web host hands the paths to Open CSV as a Set…'s open panel. (b) The
  module's checks of `IN_APP` read as "a core answers" (W6's `hasCore`), so that a host without Tauri runs the module.
  (c) Save's status line taking the host's word for where the file went (the web host downloads it: `plots_save`
  answering `{ path, note }`, the note said), rather than "Saved: Plots/x.drv".
- **Why**: the module is unchanged in the web host only through the Tauri stand-in and the dialog detour
  (web/DESIGN.md, step 2); these three let it be hosted plainly.

### W12. The inspector's panels' ? opening their manual sections (filed by the web stream, 26 October 2026)

- **What**: the module's panels (members, regions, windows, the member) each with a ? that calls the host's `help(id)`
  with a stable id (`help-members`, `help-regions`, `help-windows`, `help-member`); the host opens the manual at that
  anchor (`<a id="help-…">` in `web/manual/user-manual.md`). The view's ? does so from 1.0.0 (the host maps the view
  shown to `help-overlay` … `help-series`).
- **Why**: Josu's ask (26 October 2026): each panel's ? opening its own section of the manual.

## To the app's shell

### W6. The host taken from the page, not from Tauri (I10 step 2; filed by the web stream, 25 October 2026)

- **What**: `core.js` takes its `invoke`, `listen` and open panel from a host object the page sets
  (`window.__DRS_HOST__ = { invoke, listen, open, ask }`), Tauri's being the app's host; `IN_APP` stays "inside the
  Mac app", and a new `hasCore` says a core answers `invoke`. table.js and plot.js ask the clipboard through the host.
- **Why**: the web host now sets `window.__TAURI__` itself (`app/ui/host-web/tauri-shim.js`) so that `IN_APP` holds and
  the module runs; a page posing as Tauri may take a path meant for the app only. With W6 the stand-in goes.
