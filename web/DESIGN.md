# web — design

Grows with the steps (README.md). The assessment and the architecture decided are I10's page; what is settled here:

- **One module, three hosts**: the Tauri app, the static page `drs plots show` writes (I08), and the browser host of
  this stream. The host provides what the app's shell provides (`invoke`, the window's state, the panels) and nothing
  the module computes.
- **The dispatcher**: `app/ui/host-web/dispatch.js` maps `invoke(cmd, args)` onto the wasm package's exports (W1),
  command names and JSON shapes those of `app/src-tauri/src/modules/plots.rs`, so that the module is unchanged.
- **Files**: picked (`showDirectoryPicker`, Chromium; `<input webkitdirectory>` elsewhere) or dropped; read in the
  browser; names and text to `plots_open_from`; handles kept for live mode. Nothing uploaded, ever.
- **Layout**: the inspector on the left (the folder's tree, the members table), `<` to shrink it, a draggable
  separator; the views at the right as in the app.

## Step 2: the host as built

1. **The module unchanged.** The page imports the module and the shell's service files from `app/ui/`. `core.js`
   decides `IN_APP` from `window.__TAURI__` when it is evaluated, so `tauri-shim.js` is the first import of `host.js`:
   it gives `core.invoke` (the dispatcher), `dialog.open` (the files read in the page: the paths the host hands it, else
   the tree's selection, filtered by the dialog's extensions), `dialog.ask` (`confirm`) and `event.listen` (nothing to
   hear). app.js is not loaded: it pulls in every module and the native window's commands. The host does for the module
   what app.js does: `state.tabs`, `state.activeTab` on the bus's `tab`, `plots.close` on `tab-closed`,
   `window.__renderPanels`, a file tab rendered by `plots.view`. A shell that takes its host from the page (W6) retires
   the stand-in.
2. **Files as virtual paths.** Each file read has a path under `/local/<folder>/`, and that root is
   `state.project.root`, so the module's Save (which needs a project) and the core's relative names work as in the app.
   The dispatcher reads the texts of the paths a command names and sends the text form: `plots_open_files` →
   `plots_open_from`, `plots_add` → `plots_add_from`, `plots_open_set`/`plots_open_buffer` → `plots_names_in` (the core
   says which files the set file or buffer names; the host never parses one), then `…_from` with those texts;
   `plots_save` → its `{ relative, text }`, the text downloaded (step 3: the package's real forms, below). A file is
   read only when a command needs it.
3. **Opening files.** The tree's selection (or a folder's "set", or files dropped loose) goes to the module's own Open
   CSV as a Set… ("a CSV of several columns, or several data files": `plots_open_files`), with the paths handed to its
   open panel; a `.drv` or a `.buf` alone goes to `plots.open` as the app's openFile sends it. The extension only
   routes, as the module's own dialogs filter; what a file is and its variable are the core's. A module entry for a
   host (W7) would replace the detour through a dialog command.
4. **The package is announced, never probed.** `wasm/manifest.json` says `"package": null` on the stub, so the page
   never asks for a file that is not there (a 404 is a console error). `build-wasm.sh` writes wasm-pack's output beside
   it and names `drs_wasm.js` in it.
5. **Layout.** The inspector is a grid column `var(--insp-w)`, clamped to 220–640 px while leaving the views 360 px
   (`separator.js`); `<` makes it a 36 px rail; a double-click on the separator does the same; ←/→ move it by 16 px
   (64 with ⇧); kept as `drs-web.inspector` in localStorage. The views are the shell's dockview with only the file tabs
   (the app's side zones are not added; their edge stubs are hidden). Under 640 px the inspector stacks above the views.

### The package the host uses (W1, core 0.79.0; L424)

`app/ui/host-web/wasm/drs_wasm.js` and `drs_wasm_bg.wasm`, from `web/tools/build-wasm.sh` (`wasm-pack build core/wasm
--target web --out-name drs_wasm`), contract in `core/wasm/README.md`:

- `default(input?)`: wasm-pack's initialiser, which loads the `.wasm` beside the module;
- `dispatch(cmd: string, args_json: string) -> Promise<value>`: the parsed result, rejected with the core's words (a
  string), as Tauri rejects; the sets kept by id inside the package, as `plots.rs` keeps them. (W5 asked for JSON text;
  the dispatcher takes either);
- `commands() -> string[]`, `version() -> string`: what it serves, and the core's version in the bar.

The text forms as the package has them: `plots_open_from { files }`, `plots_add_from { id, files }`,
`plots_open_buffer_from { text, files }`, `plots_open_set_from { name, text, files }`, `plots_save` → `{ relative,
text }`, `describe { name, text }` → `{ kind, id, name, role, source }`. `plots_names_in` (W5) and `plots_propose`
(W2) are not in 0.79.0.

## Step 3: the files routed, the proposals offered, the tests honest

1. **What a file is, the core's.** Every file of a folder is described once by the core (`describe`, by its name and
   text; the dispatcher keeps the answer by the file's size and date). The tree shows each file's kind (`xye`,
   `free_steps`, `buffer_file`, `?` for what no reader recognises) and lists the folder by **folders** or by **kinds**
   (a group per kind, data first, the project's files, then the unrecognised), the choice kept (`drs-web.tree.by`).
   A file over 64 MB is not read (`TEXT_LIMIT`).
2. **Routing.** Files opened (a double-click, Open as Plot Set, a folder's "set", a drop) go by what the core says:
   each plot set file or WinPLOTR buffer among them opens as itself, through the module's `plots.open` (as the app's
   openFile, and the tab already open brought forward); otherwise the data files open as one set (Open CSV as a Set…,
   `plots_open_files` → `plots_open_from`); what is neither is left out and named in the status line with what the
   core says it is. So a buffer dropped with its files opens as the buffer (L425's refusal). On the stub, before any
   `describe`, the extensions `.drv` and `.buf` route as the module's dialogs filter.
3. **The files a buffer or a set file names.** With `plots_names_in` (W5) the core says which; until then the host
   hands over every data file of the folder (by `describe`: an unrecognised file given to `plots_open_set_from` is
   refused, `winplotr.mlf` in the Tb folder), named as the core resolves them: a buffer's from its own folder, a set
   file's from its base (`setBase`, the core's `plots::set_base`: the project's folder when the file is in Plots/),
   with `../` where needed (`relativeTo`). The core takes its own by name and names a missing one.
4. **The proposals.** When a folder is opened the host asks `plots_propose` with the folder's data files, if the
   package serves it, and lists what comes back under the files, "Sets in this folder": each proposal's title from the
   core, its number of members, its variables' sources on hover, Open to open it (a proposal's set file is held in the
   page under `Plots/` and opened by the module's Open Plot Set, never written; or its files as one set). The answer
   is read as `{ proposals: [{ title, text | files, members?, variables? }], left_out: [{ name, why }] }` (an array
   alone taken too) until W2 settles it. Without the command the panel says the core has no proposer yet.
5. **A command the package lacks** (`plots_extent` in 0.79.0, which the module's y toggle asks) is refused by the
   dispatcher at once in the core's own first clause, "no command "plots_extent" in the browser's core (drs-wasm
   0.79.0)", without the list of every command the core appends; the module shows it in the status line. Nothing is
   logged to the console; the dispatcher counts every call (`__drsWebState().calls`) so that a check sees how often.
6. **Same files dropped again.** A drop opens what was dropped, whether new to the folder or the same files again
   (`files.add` returns where each landed).
7. **The manifest.** The committed `wasm/manifest.json` is the stub's; `build-wasm.sh` writes the built one and marks
   it skip-worktree, so a local build leaves `git status` clean; `--stub` puts the committed one back. The tests read
   the manifest as it is: the stub's, or a package whose files are there (`web/tests/package.test.mjs` runs the real
   package in node then, skipped otherwise).

## Steps 4 and 5: live mode, offline, the address, the public site

1. **The page is the host** (W6, Plots' step 21): `window.__DRS_HOST__` answers `invoke` (the dispatcher), `open` (the
   files read here: the paths the host hands the next panel, else the tree's selection), `save` (a path under the
   folder, the dispatcher downloading what is written there), `ask` and `listen`. `tauri-shim.js` is gone.
2. **Live mode** reads the folder's handle again (`fromDirectoryHandle`) on a timer and diffs the listing by path and
   stamp (`diffListing`). New files are described; if any, the core proposes again (with the files' modification
   times, `OpenedFile::modified`), and each live set (a tab opened from a proposal, kept by kind and class) gets the
   files the new proposal of its kind and class has and it lacks (`growth`), through the module's Add Files… on its
   tab (brought forward for it, the tab in use put back). Nothing is ever removed from a set: a file that vanishes is
   said. The re-proposal reads every data file again: the ESRF run's 101 take 19 s in node, so live mode proposes only
   when files are added. A background tab's timers are slowed by the browser (once a minute after five minutes).
3. **Offline.** `sw.js` installs onto an empty cache the list `precache.json` (every file the page loads, from its
   links and imports, `import()` of literals included: the Surface's three.js) and the package; it serves the
   published site cache-first (its `VERSION` stamped by `stage.mjs`, so that a new version installs anew and the old
   cache is deleted) and a local checkout ("dev") network-first, keeping the listed files fresh. A request not listed
   passes to the network untouched and is never kept: the user's files are never fetched anyway. A node test fails when
   the list is not the page's.
4. **The address.** `hash.js` encodes where the set came from (each file a `file=` parameter, a set file or buffer by
   its path, or a proposal by kind and class) and the view, the scale, the zoom, the clip of the view shown and the
   region, from the tab's state as the module keeps it (`view`, `mapScale`, `mapZoom`, `clip`, `selregion`), written
   by `replaceState` once a second when it changed. On load, when the folder is opened, the set is opened as it was and
   the state set in the tab; the zoom once the core's map is in (its first drawing fits the whole map). A file missing
   is named.
5. **The public site.** `stage.mjs` lays out exactly the files the page loads, the host at the root; the Action
   (`publish-web.yml`) builds and pushes it as `publish-site.yml` does the project site. The page says what it is and
   that nothing is uploaded (the start page and the line under the inspector) and its version (the start page, the
   offline line). Licences: dockview-core and three.js, MIT, in `LICENCES.md`; the project's own code carries no
   licence statement yet (Josu's to give).

## Step 6: the guide

1. **Where.** At `/guide/` of the same address, not a repository of its own and not at the root: the app keeps its
   address (links already sent, installed copies, the service worker's scope), and one Action publishes both, so the
   guide never describes a version that is not the one beside it. Set aside: the app moved to `/app/` with the site at
   the root (every shared link and installed copy broken for a front page).
2. **One source per manual.** `web/manual/*.md` is read on GitHub, rendered by the site and typeset as PDF; the
   generator drops the H1 and the hand-made contents (the site's outline and LaTeX's contents say them), rewrites the
   figure paths and the links between manuals. Team-only parts are marked in the source and left out of the public
   edition; a team's manual is refused by the public build.
3. **PDFs built locally, committed.** TeX and the typefaces are on the Mac; the Action renders HTML only and warns when
   a PDF is older than its manual (a hash of the source each PDF was built from, `downloads/pdfs.json`). Set aside:
   TeX Live in the Action (minutes on every push, fonts to fetch there too).
4. **Screenshots from the app, by the address.** Each scene is a view link (hash.js) on the public example folder, so
   the pictures are what a user's link shows, regenerated by one command when the app changes.
5. **The page computes nothing; the site neither.** The start page's check asks the browser what it supports and
   says it; no number about data is computed or shown that the core did not give (the screenshots' numbers are the
   core's, the recipes' numbers come from runs of the core named with them).

## Step 7: the open-tool kit in 1.0.0

1. **Records are data, everything else renders them.** Releases, roadmap and conventions are YAML in `site/data/`; the
   log is the project's own. `build.py` writes the pages, the manual's appendices, the CHANGELOG and the feeds;
   `stage.mjs` takes the page's version from the releases record; `release.py` checks that every place agrees.
2. **Private by construction.** A team's manual is marked in its own first line, never named in a script; the guard's
   words live in `private/`, outside the public source; the guide, its PDFs and the source selection are all checked
   against them, and the log leaves out (and counts) an entry that names a team's data.
3. **Two branches, one repository.** The public repository's `main` holds the source (by `publish-source.yml`), its
   `gh-pages` the built page and guide (by `publish-web.yml`); the source Action refuses to run while Pages still
   serves `main`, so the page is never replaced by its source.
4. **The example is computed.** No permission to ask, nothing to cite but the core; reproducible to the byte.
5. **Screenshots are a test.** The check reads the view's main canvas (its largest): healthy views give 155 to 1047
   colours in a sample, the blank Overlay of 0.4.0 gave 67, so the threshold is 120; a WebGL view by its capture's
   size.
