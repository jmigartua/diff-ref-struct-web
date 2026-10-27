# web

The analysis in the browser (I10, R34, L423): the core compiled to WebAssembly and run in the browser, the Plots
module and the Map viewer as they are behind the single `invoke` they use, the user's files read in the browser and
never uploaded, the series recognised by the core and its maps offered; a static site, installable, usable offline,
with a live mode on a folder (Chromium) for the instrument's case. Public as **diff-ref-struct Web** at
<https://jmigartua.github.io/diff-ref-struct-web/>, the repository `diff-ref-struct-web` holding the built app only,
pushed by an Action from this repository as the project site is (`publish-site.yml`).

Worked in its own stream: worktree `../20261009_diff-ref-str-web`, branch `web`. It owns this folder and
`app/ui/host-web/` (the browser host: the page, the dispatcher onto the wasm core, the inspector, the file access, the
service worker) and the Action `publish-web.yml`. It edits nothing else: what it needs from the core, the Plots module
or the app's shell goes into `REQUESTS.md` here (W1, W2 …); the merging session passes it on and keeps the streams in
step. The module's guard holds over the host (PL14): the browser computes nothing, the core does, in wasm (L139).

The assessment, the five steps and the asks answered are I10's page on the project site
(`ecosystem/data/ideas.yml`); the design grows here in `DESIGN.md` as the steps are built.

## Manuals

- `manual/user-manual.md`: the page for its users, every control, for web a1c308d with core 0.81.0 (26 October 2026);
  the limits of that version said plainly (§20).
- a team's worked use case (an ESRF ID22 in-situ run) in `manual/`, for its experimental team
  (unpublished data: not for the public site without the team's word). It gives a one-line conversion of the
  project's `*-map.drv` sets to `Data/` with `rebin = 0.01`, checked in the core 0.81.0 package in node (the three
  sets open with their clips and regions, the maps equal to the desktop's: 2699 columns, 0–28.26 √counts for the
  hold).

## Steps (from I10's proposal, 25 October 2026)

1. The core in the browser (the core stream, W1): `plots::open_from`, a `wasm` build, the Tb₂Pd₂Sn map to the bit.
2. The web host: the app's shell and the Plots module behind the wasm dispatcher; the inspector on the left with `<`
   and a draggable separator; a folder picked or dropped; the Map with the clip, cursors, zoom and regions on local
   files. Hosted locally first.
3. The series detector (the core, W2) offered by the host; the chosen proposal opened.
4. Live mode and offline: the folder read again on a timer (Chromium), the service worker, the `.drv` down and up,
   the view in the hash.
5. Public: `diff-ref-struct-web` built and pushed by `publish-web.yml`; a page for the instrument's users.

## Step 2: the browser host (25 October 2026), hosted locally

`app/ui/host-web/` is the page. It loads the shell's service files and the Plots module from `app/ui/` as they are
(relative imports `../core.js`, `../layout.js`, `../registry.js`, `../modules/plots/index.js`; no copy, no bundler) and
serves their single `invoke` with `dispatch.js`. Run it:

    web/tools/serve.sh            # app/ui/ on http://localhost:8765/, the page at /host-web/
    node --test "web/tests/*.test.mjs"   # the pure parts (a glob, for Node 22 too)
    web/tools/build-wasm.sh       # once core/wasm exists (W1): the package into host-web/wasm/, named in its manifest
    web/tools/build-wasm.sh --stub

`serve.sh` needs `app/ui/vendor/` (dockview, three), which `app/tools/vendor.js` generates (`cd app && npm ci && node
tools/vendor.js`; it is gitignored, as in the app).

| file | what it does |
|---|---|
| `index.html` | the page: the bar (Open Folder…, Open Files…, the Plots menu, the core's state, the theme), the inspector at the left, the separator, the views (dockview), the status line, the shell's three dialogs |
| `host.css` | the page's own layout on `../style.css` (the kit's variables, light and dark) |
| `host.js` | the shell for the module: the tabs, the panels it draws into, the tree, the files into the module, the menu, the width and the rail |
| `tauri-shim.js` | `window.__TAURI__` stood in for, so that `core.js`'s `IN_APP` and `invoke` hold: its `invoke` the dispatcher's, its open panel the files read here (until W6) |
| `dispatch.js` | `invoke(cmd, args)` onto the package's `dispatch`; the path commands in their text forms (W5); the stub's refusal |
| `files.js` | the folder picked or dropped, its handles kept, each file's text read when the core needs it, downloads |
| `tree.js`, `separator.js` | the tree from a file list; the width clamped and kept (pure, tested in node) |
| `wasm/manifest.json` | `"package": null` on the stub; `build-wasm.sh` names the package there |

On the stub every command is refused with one sentence: *The core is not in this page yet: "plots_open_files" waits
for its WebAssembly build (request W1, web/REQUESTS.md). Nothing is computed in the browser, so nothing is answered
without it.* That sentence appears in the status line when a set is opened. The clipboard and a text downloaded (the
viewer page) are the host's own and answer anyway.

Checked in Chrome (Interceptor, the isolated profile): the page loads on the stub with no console error or warning,
no uncaught exception and no failed request (an error recorder installed before the page's first script); a folder
dropped (Chromium's directory handle) becomes the tree, hidden files left out, `tb_2K` before `tb_10K`, its handle
kept; a ⇧-click range and Open as Plot Set reach the module, which shows the W1 refusal; a `.drv` dropped alone goes to
Open Plot Set (`plots_open_set`, refused the same way); `<` makes a 36 px rail and back; the separator's arrow keys stop
at 640 px and the width is kept in localStorage; the Plots menu lists the module's 55 commands. The drop was a
synthetic event carrying a directory handle, not a drag from the Finder, which waits for a check by hand.
Not yet possible without W1: the Map itself (its clip, cursors, zoom and regions are the module's and come with the
core's answers), Save Set's download and Open Set's members.

## Step 3: the files routed, the proposals offered, the tests honest (25 October 2026)

On core 0.79.0 (drs-wasm, L424), the package built with `web/tools/build-wasm.sh` (1156 kB; the manifest it writes is
skip-worktree, so `git status` stays clean; `--stub` returns to the committed one). The dispatcher now speaks the
package's real contract (`core/wasm/README.md`): the merging session's await kept; `plots_save` answered as `{ relative,
text }` (it had been mapped to a `plots_save_text` that does not exist); `commands()` read, so that a command the
package lacks is refused at once in the core's words.

- **What a file is** comes from the core's `describe`, asked once per file: each row of the tree carries its kind, and
  the tree lists the folder by folders or by kinds.
- **Routing**: a `.buf` or a `.drv` among the files opened (picked, double-clicked, dropped with its files or alone)
  opens as itself, through `plots_open_buffer_from` / `plots_open_set_from`, with the folder's data files named as the
  core resolves them (until `plots_names_in`, W5, says which); data files alone open as a set; anything else is left out
  and named. The same set file opened twice brings its tab forward.
- **Proposals**: "Sets in this folder" under the files asks `plots_propose` when a folder is opened; 0.79.0 has none,
  and the panel says so (W2, core 0.80.0).
- **The y toggle** (`plots_extent`, not in 0.79.0): the status line says `no command "plots_extent" in the browser's
  core (drs-wasm 0.79.0): the app's core has it, this build not yet`; nothing in the console.
- **Tests** (`node --test web/tests/`): `host.test.mjs` (23) reads the manifest as it is, the stub's or a built
  package's; `package.test.mjs` (5) runs the real package in node through the host's dispatcher (the Tb folder's kinds,
  both buffers with 20 members, Save and the `.drv` opened again, the refusal of a missing command, the ESRF folder),
  skipped when no package is built.

Checked in Chrome (Interceptor, isolated profile; the repository and the ESRF folder served from one scratch root, the
page in a frame with an error recorder put in before its first script, files dropped through a DataTransfer or a
directory handle in eval): no console error or warning throughout.
- The Tb₂Pd₂Sn folder: 23 files described (2 `buffer_file`, 20 `free_steps`, 1 unrecognised, `winplotr.mlf`);
  `Tb_vs_T.buf` double-clicked opens with 20 members and its map (19 rows, one member apart); `data_buffer_file.buf`
  dropped loose with its 20 files on a fresh page opens as the buffer (L425's refusal gone), dropped again brings its
  tab forward; ten `.dat` alone make a 10-member set; two `.dat` with the `.mlf` a 2-member set, the `.mlf` named left
  out; the `.mlf` alone "Nothing to open as a plot set".
- The ESRF `renamed/` folder (101 `.xye`, 112 MB, read only): read in 1.4 s, all 101 described `xye` in 66 s; the
  proposals panel says the core has no proposer yet; the 101 opened as one set in 20 s more, every member drawn apart
  (no `T`: the core reads `T300`, `_300K`, not `T0900C`; W8); then the module's Overlay of 101 × 27 000 points blocked
  the page for more than ten minutes (W9; in node the core opens them in 30 s and maps them in 2.6 s).
- The y toggle's refusal in the status line, `plots_extent` asked once per map and toggle.

Filed: W8 (the core: `T0900C` in °C and `t00000s` from a name), W9 (the core and Plots: a large set's patterns out of
every SetView, and its overlay drawn in reasonable time). W6 and W7 stay with Plots and the shell.

## Steps 4 and 5: live mode, offline, the address, and the public site (25 October 2026)

On core 0.80.0 (drs-wasm 1251 kB). The page is now the modules' host as Plots' step 21 defines it
(`window.__DRS_HOST__ = { invoke, listen, open, save, ask }`, W6): the Tauri stand-in is gone and `IN_APP` is false in
the browser. Save answers `{ path, note }` (W7): the module says "Downloaded as <name>.drv … nothing is written into
the folder". A `.buf` or `.drv` goes with the files `plots_names_in` names (W5), by their paths.

- **Live mode** (Chromium, a folder opened by Open Folder… or dropped): the folder read again every 2 to 60 s (Live,
  every, Now; kept); new files appear in the tree and are described; when files are added the core proposes again,
  and a set opened from a proposal grows with the files the core now proposes for the same kind and class, through the
  module's Add Files… (`plots_add_from`), its map asked again; a file that vanishes is said in the live line and the
  status line and stays in its set. `live.js` holds the diff, the growth and the line.
- **Offline**: `sw.js` keeps the files the page loads (`precache.json`, 33 files, from `web/tools/precache.mjs`, which
  walks the page's links and imports) and the package (`wasm/manifest.json`), 36 files, 4.4 MB, and never anything
  else; `manifest.webmanifest` and `icon.svg` make the page installable. The line under the inspector says what is
  kept and the version.
- **The address**: the set shown (its files, its set file or buffer, or the core's proposal by kind and class), the
  view, the colour scale, the map's zoom, the clip and the region kept in the hash (`hash.js`), and brought back when
  the folder is opened again with those files; otherwise the page says which are missing.
- **The public site** (`web/tools/stage.mjs`): the files the page loads, the worker and the package, the host at the
  root (`../core.js` → `./core.js` in the host's own files only), the version stamped in `sw.js` and `version.json`,
  `LICENCES.md` with the two third-party libraries the page ships (dockview-core 4.13.1 and three.js 0.186.1, both MIT,
  their texts in `web/licences/`), a README. 36 files, 4.5 MB. `.github/workflows/publish-web.yml` builds the package,
  generates the vendored libraries, runs the tests and the offline list's check, stages and pushes to
  `diff-ref-struct-web` with the deploy key `WEB_DEPLOY_KEY`, on a push to main that touches `app/ui/`, `core/`, `web/`
  or itself.

Tests: `node --test web/tests/` (host 23, package 5, steps 4 and 5: 10).

Checked in Chrome (Interceptor, isolated profile), the repository served from a scratch root:
- live: a scratch folder (`live/`, never `renamed/`) holding three ESRF files copied with their times; the core
  proposed "3 scans against time, 27.033 min", opened as a live set read every 2 s; the 4th and 5th files copied one
  by one each grew the set by one (4, then 5 members); with the Map shown the 6th grew it to 6 and the map was asked
  again; a file removed was said ("vanished … the set keeps it as read") and the set kept 6;
- offline: the server stopped, the page loaded again from the worker, the core 0.80.0 running; the worker keeps 36
  files, 4.4 MB, none of them a data file;
- the address: `#…set=buf:Tb_vs_T.buf&view=map&scale=log&zoom=20,30,5,20&clip=3.2,4.5` loaded, the folder dropped,
  the buffer opened on the Map zoomed to 2θ 20–30°, T 5–20 K, its colour range 3.2 to 4.5;
- the staged site served at its root: the Tb buffer opened through the proposals and the routing, no console error.
- Not yet good: an ESRF set of raw scans (27 000 points each) blocks the page while the Map or the Overlay is drawn,
  for 94 s with 5 members (W9; the core's rebin on opening, 0.81.0, passed by the host when a proposal says it).

## Step 6: the guide, the tool's own site (26 October 2026)

`web/site/` is a Quarto site in the family kit (theme.css, site.css, fonts.html, docs-close.html copied unchanged from
`ecosystem/site`; this site's rules in `project.css`, the accent the app's green #0a7f6f), published by
`publish-web.yml` beside the app at **/guide/** (the app stays at the root, so that its links, its installed copies and
the address's views hold). The app's start page links to it and its ? opens the manual there.

| page | what |
|---|---|
| Overview | what the page does, the seven views with screenshots, why the numbers can be trusted, privacy, who it is for |
| Start | a check of what this browser supports (feature detection only), then the public Tb₂Pd₂Sn folder (FullProf's WinPLOTR tutorial) from opening to a pinned map, every number the core's |
| Manual | `web/manual/user-manual.md` as a page, and as PDF; the quick reference (`web/manual/quick-reference.md`), page and one-page PDF |
| Recipes | a ramp, a hold (with a tested script that writes a `.drv` typing each scan's time), two patterns compared, live at the instrument, sharing |
| How it works | the core in wasm, the checks, the files never sent, offline, the publication; the formats (the manual's §21, one source) |
| Versions, About | `data/versions.yml`; credits, licences, citing |

One source for each manual: `tools/build.py` writes the site's copy (`_gen/`) and `tools/pdf.sh` the PDF's from the
same Markdown; a part between `<!-- team-only -->` marks is left out of the public edition. The PDFs are typeset by
LuaLaTeX in the site's typefaces (Source Sans 3, Source Serif 4, IBM Plex Mono, SIL OFL, fetched into `.fonts/`) and
committed in `downloads/`; `build.py` warns when one is older than its manual (`--check` fails). The ESRF use case is
a team's edition (unpublished data): `pdf.sh team` writes it as PDF and one self-contained HTML page into
`web/manual/team/`, never into the site; the Action stops if a trace of it reaches the guide.

Screenshots are taken from the app, never drawn: `web/site/tools/shots.mjs` drives Chrome over the DevTools protocol
(Node's WebSocket; `--experimental-websocket` on Node 20) on the locally served page, the folder through the page's
own input, each scene's view set through the address, as a user's link would.

    web/tools/serve.sh 8790 &
    node --experimental-websocket web/site/tools/shots.mjs --port 8790     # assets/figures/: views at 3×, windows at 2×, light and dark (WebP)
    web/site/tools/pdf.sh                                                  # downloads/*.pdf, ../manual/team/*
    cd web/site && quarto render                                           # _site/, as the Action renders it

Found while taking them, and fixed in the host: the **Overlay drew no curves** on the public page (web a1c308d): the
package's `plots_patterns` names each pattern's member `index` where the module reads `member`
(`app/src-tauri/src/modules/plots.rs:190` against `core/wasm/src/host.rs:112`). The dispatcher now renames it (`SHAPES`,
no number touched, tested), W11 asks the core to answer so and to test the shapes. And **Sum, Average and Difference of
Members** are downloaded (the core gave their text; nothing saved it). Tests: `node --test web/tests/*.test.mjs` 40 of
40 with the package built at core 0.81.0.

The landing's pictures (26 October 2026, later): each view cropped to its plot areas and taken at 3× in WebP, in the
light and the dark theme (the page's `<picture>` takes the one of the site's theme), the Map as a wide feature card,
every screenshot opening at full size in a lightbox (`assets/includes/lightbox.html`). The whole-window PNG is kept
only for the PDF's figure. The open-tool kit, the standard proposed for every public tool of the family (log,
releases in the manuals, roadmap, Conventions, requests, licence and citation, the skill), is `PROPOSAL-open-tool.md`.

## Step 7: 1.0.0, the open-tool kit applied (26 October 2026)

Josu's decisions (all of `PROPOSAL-open-tool.md`'s recommendations, English only for now; to be logged, W14): open
source (code Apache-2.0, documents CC BY 4.0), the source published as a subset by an Action, requests as the public
repository's issues with an e-mail path, the appendix named Conventions, semantic versions from 1.0.0, the example
computed by the core, the AI-assistance statement (proposed wording in `public/NOTICE` and the guide's About).

| what | where |
|---|---|
| the records | `site/data/releases.yml` (semantic versions; the page's version, CITATION's, the manuals' Appendix B, `CHANGELOG.md`, a feed), `roadmap.yml` (Appendix C), `conventions.yml` with the project's file dictionary (Appendix A) |
| the guide's new pages | Releases, Log (the project's log, entries about the page, none naming a team's data), Roadmap, Requests (the public issues, read at render), Conventions, Privacy and accessibility; the known issues on the home page; Atom feeds |
| the page | Help ▾ (the guide, the manual at the view shown, the manual of this version, what's new, Report a problem…, Ask for a feature…, e-mail, Cite this view, About); the view's ? at its manual section; Try the example; the version `diff-ref-struct Web 1.0.0 (build …)` |
| offline | the guide kept by the page's worker too (`guide/offline.json` from `site/tools/offline.py`), network first |
| the example | `examples/tools/make_example.py` → `app/ui/host-web/examples/ramp/`: 21 scans computed by `drs simulate` (periclase 300–1300 K, corundum above 1000 K), byte-reproducible, CC BY 4.0 |
| open source | `public/` (README, LICENSE, LICENSE-docs, NOTICE, CITATION.cff, .zenodo.json, CONTRIBUTING, CODE_OF_CONDUCT, the issue forms), `tools/source.mjs` (the public selection, guarded by `private/guard.txt`), `.github/workflows/publish-source.yml` (by hand or a tag `web-v*`, to main; the page and guide to gh-pages by `publish-web.yml`) |
| the release gate | `site/tools/release.py [--shots]`: versions, notes, PDFs, the guide's links and guard, the source's guard, licences, tests, every view drawn (`shots.mjs --check`, which fails the 0.4.0 Overlay: 67 colours against ≥ 120) |

The gate today: every line OK but **source**, which stops on four of the core's files that name the team's run (W13).
Until W13 is met, the source is not published; the page and the guide are, by `publish-web.yml`, once Pages serves
`gh-pages` (W14).
