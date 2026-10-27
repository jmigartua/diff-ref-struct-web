# Changelog — diff-ref-struct

Release notes. Generated from `ecosystem/data/releases.yml` by `ecosystem/tools/build.py`; the same text is on the project site. Do not edit by hand.

## web 1.1.0 — 2026-10-27 — diff-ref-struct Web 1.1.0: the 2θ interval of the analysis

On core 0.82.0 (L442): a set file's x_interval cuts every member to an interval of 2θ as the set opens, every view, window and map working on the cut set. Archived as doi:10.5281/zenodo.23683238. Its own release notes: https://jmigartua.github.io/diff-ref-struct-web/guide/releases/

### Added

- A set file's `x_interval = [from, to]` honoured by the page (core 0.82.0, PL42).

### Changed

- The page's core 0.82.0; the patterns in the app's words (W11 met).

### Known limitations

- No control in the page for the interval yet (Plots' step 25).

## web 1.0.0 — 2026-10-27 — diff-ref-struct Web 1.0.0: the guide, open source, the record

The page in the browser gets its own site (/guide/) with the manual in HTML and PDF, a licence and a citation, a place for requests, and a record (releases, log, roadmap, conventions) rendered in the guide and the manuals; its versions become semantic (L436, L437). Its own release notes: https://jmigartua.github.io/diff-ref-struct-web/guide/releases/

### Added

- The guide at /guide/: start, manual and quick reference (HTML and PDF), recipes, how it works, conventions, releases, log, roadmap, requests, about.
- Help ▾ in the page: the manual at the view shown, what's new, report a problem, ask for a feature, cite this view; Try the example (computed by the core); the guide offline.
- Open source: Apache-2.0 for the code, CC BY 4.0 for the documents, CITATION.cff, issue forms.

### Fixed

- The Overlay draws its curves in the browser (the core's patterns read under their name, W11).
- Sum, Average and Difference of Members download their file in the browser.

### Known limitations

- Proposals against time open by temperature (W10); the source is published once the core's files name no team's data (W13).

## app 0.54.0 — 2026-10-25 — Plots' steps 20 to 22 in the Mac app: the y interval shared or own, the richer members table, the right-hand axis, a large set drawn without blocking, on core 0.81.0 (L430)

What the web app's round brought to the shared Plots module, in the desktop app too: the controls Josu asked for on the present page (I10 step 0), the host taken from the page, and a hundred raw scans opened and drawn in seconds instead of a blocked window.

### Added

- Plots, step 20: a toggle of the y interval, shared by the whole set or adapted to the member shown (the row strip, Overlay, Stack), kept in the set file (`y_interval`, core 0.80.0); the members table with an icon by the file's kind, a checkbox per member and all / none / by exposure class / by a range of the variable, three orders, the member's facts on hover (W4, L425).
- Plots, step 20 and 21: the Map's and the Difference Map's right-hand axis labelled with a member's second variable with its source (core 0.79.0, W3).
- Plots, step 22: a set's patterns asked when a view draws them and no longer carried in every answer; the Overlay of a hundred members of 27 000 points drawn in slices across frames; the core's rebin to 0.01° offered when a large set is opened (W9, L427–L428).

### Changed

- The shell takes its host (invoke, listen, the panels, ask) from a host object the page sets, Tauri's being the app's; `hasCore` says a core answers (W6, step 21). The app behaves as before: the self-checks say so.

### Known limitations

- The exposure classes a run's files make are told apart by the counts' scale σ²/y (core 0.81.0, L428), Josu's to accept.

## app 0.53.1 — 2026-10-25 — A region's own colour range per view (P65, L420)



### Fixed

- A region saved on the difference map keeps its colour range in its own field of the set file (`difference_range`, core 0.78.0) and applies it to that map only; one saved on the map likewise; Save no longer keeps one in the tab (Plots' step 19).

## drs-publish 0.36.1 — 2026-10-25 — A region's own colour range per view (P65 closed, L420)



### Changed

- F3's difference map takes a region's own `difference_range` as its clip, and the map the region's `map_range`, each left aside on the other without a note (core 0.78.0 gives `Region` both).

## app 0.53.0 — 2026-10-24 — The interactive map of a series (I08, R33): the clip on the colour bar, cursors, the zoom in the T–2θ plane, regions, and the map's page without the app (L405–L417)

Plots' Map and Difference Map read for phase transitions: the colour range clipped on the bar itself so that the weak details come up, two sections through the map under the pointer, a zoom in the temperature against 2θ plane with the clip following it, regions saved by name in the set file, and the same map in a page the command line writes without the app. On core 0.77.0 and drs-publish 0.36.0.

### Added

- Plots, step 16: two marks on the colour bar set the range shown (Set Colour Range… types it), values beyond saturated, per view and on the scale shown, a double-click back to the core's range; the same clip on the Surface's height; a row cursor (that member's pattern below the map) and a column cursor (that 2θ's intensity against T beside it), pinned by a click (L406).
- Plots, step 17: a rectangle zooms both axes, the wheel over an axis zooms that axis, a double-click shows the whole; Range from the View clips to the core's block of the rectangle; Save Region… names a region (a 2θ window times a range of T, in the quantity shown), listed beside the windows, saved in the set file with the clip (core 0.76.0, PL37) (L408).
- Plots, step 18: the Map viewer as one module hosted by the app and by the page `drs plots show <set.drv> --browser` writes beside the set (core 0.77.0, PL38): the same map in a browser on a machine without the app, nothing fetched (L412, L414).
- `open -a diff-ref-struct <set.drv> --args --view map|difference|surface` opens the set in Plots on that view, which `drs plots show` uses when the app is installed (PL39, L417).
- Publish's F3 map and difference map take the set's clip and a named region, the caption saying the range shown and the saturation (drs-publish 0.36.0, PL40, L416).

### Known limitations

- A region's own colour range does not say which map it was set on (Publish's P65): Publish takes it for the map only; the core's and the module's fix follows.
- The strips of the cursors need a plot area over 260 px tall (the row) and 520 px wide (the column); a short window shows the mark and the readout only.

## drs-publish 0.36.0 — 2026-10-24 — F3's map takes the plot set's clip and a named region (PL40, L416)



### Added

- The map and the difference map drawn over the set file's map_range or difference_range, every value beyond it in the end colours, the colour bar's axis the clip's range (I08 step 2, R33).
- A figure's region = "…": the map is the core's block of a region the set saves (plots::map_region), on its bounds, with its own colour range when it has one; a region the core refuses refuses the figure in the core's words.
- The caption says the region, the range shown, what is saturated and the core's range.

### Changed

- F3's map is drawn on the set's map_scale (else y_scale), the difference map on the core's difference_scale of it.

### Known limitations

- A region's own map_range is taken for the map only, since the set file does not say which map it was set on (publish/REQUESTS.md P65).

## drs-publish 0.35.0 — 2026-10-24 — The checks before export as DESIGN.md §5.4 writes them (L415)



### Added

- Check 2: T1–T6 read back from their LaTeX row by row, each number against the results file's to the digit it is written to.
- Check 3: every phase's block; checkCIF's powder tests POWD_002–004 as errors; no su of (0).
- Check 4: each figure's height, its texts and lines read from its PGF, each raster's resolution and pixel width.
- Check 5: pixel sizes within a kind, one-line labels, an inset over the data measured on the drawing.

### Fixed

- A set named without its folder (check x.drp inside Publish/) finds its sources again.

## 0.76.0 — 2026-10-23 — The plot set file keeps the map's clip and its named regions (Plots' PL37, I08, R33; L407)



### Added

- `SetFile::map_range` and `difference_range` (`plots::Range { from, to }`, on the map's scale), read, written and described; absent by default.
- `SetFile::regions` (`[[regions]]`): `plots::Region { name, x, x_quantity, variable, map_range, window }`, `x_quantity` the set's `x` when absent.
- `plots::map_region(&Map, &Region, &MapX) -> Result<MapBlock, String>`: the rows and columns of a map inside a region and the smallest and largest value of its cells, the region's bounds converted to the map's x quantity with each row's λ; `plots::MapX::of` (a set as its files give it) and `MapX::of_view` (a set shown through `x_view`).
- `plots::range_warning`: a colour range wholly outside a map's values said; `SetFile::check_regions` (a region the set cannot be converted to, the members without λ named); `SetFile::x_quantity`, `SetFile::words`.
- `drs plots describe <set.drv> [--json]` and `drs plots region <set.drv> --name <n> [--scale …] [--reference k] [--json]`; `drs plots map` says when the set file's clip lies wholly outside the map's values.
- Test `pl37_the_maps_clip_and_its_regions` in `plots_sets`.

### Changed

- The file dictionary's `plot_set` names the clip and the regions; a `[[member]]` and a `[[regions]]` are the set file's two tables.

### Known limitations

- `map_region` takes the map's x as a third argument (`MapX`), not in PL37's two-argument form: a `Map` does not say which quantity its columns are in nor each row's λ (the app's map in d is built from `x_view`'s patterns while the members' facts stay in 2θ), so the caller states it; Plots' step 17 builds it with `MapX::of_view`.
- The app's Tauri commands that read and save a set do not carry the new fields yet: this release edits only `core/`, the binding is Plots' step 17.

## 0.77.0 — 2026-10-24 — The Map viewer's page without the app: plots::page_data and drs plots show (Plots' PL38, L413)

L413: the page's data in one JSON in the shape the app's Map viewer reads, and a self-contained page holding it and the viewer, written beside the set, or the app opened onto the set when it is installed.

### Added

- `plots::page_data` and `PageData` (with `PageMap`, `PageCells`, `PageApart`, `PageMember`, `PageRegion`, `PageBlocks`): the map, the difference map, the clip, the regions with their blocks or refusal, the members with their temperatures' sources and patterns, the core's version and the date.
- `plots::page_html`, `difference_scale`, and the viewer's assets embedded (`PAGE_SHELL`, `VIEWER_BUNDLE`, `VIEWER_SOURCES`).
- `drs plots show <set.drv>`: the app when installed (macOS), else `<set>.map.html` beside the set, opened in the browser; `--view`, `--browser`, `--out`, `--no-open`, `--no-patterns`.
- `MapBlock::counted_from` (0): where a block's rows and columns count from.
- Test `plots_page` (6 cases).

### Changed

- `drs plots region` prints rows and columns counted from 1 and their indices from 0, and the scale as its word.

### Known limitations

- The page's shell opens on the map: it does not yet read `view`, so `--view difference` reaches the page only as the field.
- That the bundle agrees with the app's sources is checked by the bundler (`--check`), not by the core's test, which checks the asset against its own SOURCES.txt.
- The app taking `--view` waits for PL39; an app before it opens Plots on the view the file says.

## 0.78.0 — 2026-10-25 — A region's colour range per view (Publish's P65, L419)

L419: a saved region keeps its own colour range on the map and, apart, on the difference map, each on its map's scale, so that each is applied to its own view only.

### Added

- `Region::difference_range` beside `map_range`: read, written, refused when reversed, carried in `page_data`'s regions.
- Test `p65_a_regions_range_per_view` (plots_sets); the page's JSON holds both ranges (plots_page).

### Changed

- `drs plots describe` names the map each range is on and its scale.
- `drs plots region --reference` draws the difference map on its own scale (√ for a log map) and prints the region's range for the map shown, naming the other.

### Known limitations

- The app's `plots::Region` literal (`app/src-tauri/src/modules/plots.rs:655`) needs the new field before the app builds on 0.78.0: the Plots stream's half of P65.

## 0.79.0 — 2026-10-25 — The core in the browser: a set from names and text, the crate drs-wasm, a member's other variables (W1, W3; L424)

L424: the Plots module's commands compiled to WebAssembly with the app's names and JSON, files given by name and content; each member may carry other variables than the set's, with their sources, for the map's right-hand axis.

### Added

- `plots::open_from`, `open_from_with` and `OpenedFile`: a set from names and text, nothing read from disk; `SetFile::open_from` for a `.drv` whose members' files are given by name.
- `Member::others` (`Other`), `SetMember::variables` (`variables = { T = { value, source } }`), `MapRow::others` (`RowOther`): W3.
- The crate `core/wasm` (`drs-wasm`): `dispatch(cmd, argsJson)`, the Plots module's commands by name (the app's, with `_from` where files come by name and text), `describe`; `wasm-pack build --target web`; the `.wasm` 1155 kB.
- Tests `plots_open_from` (3 cases) and `core/wasm/tests/map.test.mjs` (5 cases, `npm test`).

### Changed

- `plots::open_with` reads its files' bytes and shares `open_from`'s reader; what it gives is unchanged (the suite).

### Known limitations

- Text files only in the browser: a binary format passed as text is not its bytes.
- Not ported to the browser: Data's operations, a control file's λ, `.seq`, the results series, the CSV export, reload; `page_data` (it reads the clock, which panics in wasm32).
- The `.wasm` is built without wasm-opt; `t` (the elapsed time) has no source yet (W2).

## 0.80.0 — 2026-10-25 — The series a run's files offer, the path commands as text in the browser, the y interval and the set's extent (W2, W3, W5, PL41; L426)

L426: from a run's files, the sets they offer (by temperature, the holds, against time, per exposure class), each a set file to open as it is; t from the elapsed seconds in the names; in the browser, the files a buffer or a set file names and those commands by path; the curves' y interval kept by the set file and the set's extent.

### Added

- `plots::propose`, `Proposal`, `Proposals`, `ExposureClass`, `SAME_TEMPERATURE`, `SAME_WRITE`; `drs plots propose` (W2).
- `plots::elapsed_in_name` and `PlotSet::fill_elapsed`: t since the first member when every name gives seconds (W3); `OpenedFile::modified`, `Member::modified`; `Other::unit`, `RowOther::unit`.
- `SetFile::y_interval` (`YInterval`: shared, own) and `plots::extent` (PL41).
- drs-wasm: `plots_names_in`, `plots_save_text`, `plots_propose`, `plots_extent`; `plots_open_buffer_from` and `plots_open_set_from` by `{ path, text, files, root }` (W5).
- Tests `plots_propose` (5 cases), `pl41_the_y_interval_and_the_sets_extent` (plots_sets), `core/wasm/tests/commands.rs` (the package's commands against the app's); `map.test.mjs` from 5 cases to 8.

### Changed

- A set opened by path carries each file's modification time (`Member::modified`); one from names and text carries the host's, if given. `RowOther`'s JSON has `unit`.
- drs-wasm's `plots_open_from` names the members relative to the first file's folder, as the app does, and takes `reader` and `title`; `MemberView` carries `modified`; `plots_save` returns `path` beside `relative`.

### Known limitations

- W2's test on the ESRF run's 101 files waits for their names (`file_index.csv`): not given to this session, the files never copied into the repository; `elapsed_in_name`'s pattern is unchecked against them.
- A modification time is taken only by `propose` and only when no two files' fall within a second; never on opening.
- The `.wasm` is built without wasm-opt: 1254 kB (1155 kB at 0.79.0).

## 0.81.0 — 2026-10-25 — The ESRF run's three sets found, the exposure classes by the counts' scale, a large set rebinned on opening (W2, W8, W9; L428)

L428: `drs plots propose` on the ESRF run's 101 files offers the three sets of its L021 as they were built by hand; a name's °C read and kept; the exposure classes told apart by σ²/y; a large set rebinned to 0.01° on opening, and the members' patterns asked for apart instead of carried by every SetView.

### Added

- W8: `T0900C`, `900C`, `900°C` read as °C and `t00000s` as elapsed seconds; `Found::unit`, `Variable::in_unit`, `PlotSet::temperatures_in`: a set in °C when its sources say °C, else in K with the °C converted.
- `PlotSet::rebin`, `SetFile::rebin` (`rebin = 0.01`), `data_ops::rebin_step` public; `CLASS_STEP`, `LARGE_SCAN`, `LARGE_STEP`; `Proposal::reader`.
- drs-wasm: `plots_patterns { id, members }`, `plots_open_from`'s `rebin`, `plots_set_other` and `MemberView.others` (Plots' step 21).
- Test `plots_esrf` (2 cases: the 101 names; with the folder, the three sets); `map.test.mjs` 9 cases (the ESRF one with the folder).

### Changed

- `propose`: classes by the median of σ²/y above 2° (`ExposureClass::median_scale`, was `median_esd`); one set by temperature over all classes (the longest exposure's last scan at each), one against time over all classes, the longest hold of each class against the time since its own first scan; time in min up to 5 h.
- drs-wasm's SetView carries no member's x, y, σ (W9): `plots_patterns` gives them.
- The table recogniser tries no delimiter absent from the data and no decimal comma without a comma: the same answer, half the time.
- `elapsed_in_name` reads a lower-case t before the seconds (`hold_t12.5sec` is now 12.5 s; L426's test had guessed it was not read).

### Known limitations

- The exposure classes by σ²/y in place of the σ/y the order asked for: Josu's to accept (L428).
- The app's `plots_patterns` binding and the module's drawing from it are Plots' (W9): until then the package's members carry no pattern.
- `plots_propose` on the 101 ESRF files takes 10 s in node; `plots_open_from` 9.7 s, 2.4 s with the proposal's reader.

## 0.82.0 — 2026-10-26 — The x interval of the analysis, set once (PL42; L435)

L435: a set file may cut every member to an interval of x on opening, so that every view, window, derived file, export, page and proposal works on the cut set; `drs plots … --x-interval a:b` the same.

### Added

- `SetFile::x_interval` (`x_interval = [from, to]`), `PlotSet::cut`, `PlotSet::interval_of`, `PlotSet::x_interval`; `propose(files, x_interval)`.
- `--x-interval a:b` on every `drs plots` command that opens a set.
- drs-wasm: `xInterval` on `plots_open_from`, `plots_open_buffer_from`, `plots_open_set_from`, `plots_propose`; `plots_rebin` (the app's, step 22).
- Test `pl42_the_x_interval_of_the_analysis` (plots_sets); `map.test.mjs` 10 cases.

### Changed

- `windows` refuses a window outside a cut set's interval in those words.
- drs-wasm's `plots_patterns` answers in the app's `PatternView` (`member`, the members shown); `Settings` carries `y_interval`.

### Known limitations

- `open`, `open_from` and `open_buf` keep their signatures: their callers cut (the app's three opening bindings take the argument in Plots' step 25).
- Publish draws the cut set already; its notes do not yet say the cut.

## 0.75.0 — 2026-10-21 — The cell's box for any phase (Explore's E22, L401)



### Added

- `Cell::edges`: the cell's twelve edges in Cartesian Å, for any phase; `geometry::scene` uses it.
- Test `explore_e22`.

## 0.74.0 — 2026-10-21 — Per-format reading options and an order of preference among the readers (the app's request 51, L400)

Each reader declares the options it takes, reads with them, and a control file keeps them with its data; an order of preference decides among readers that read a file equally well.

### Added

- Module `formats::options`: `FormatOption`, `OptionKind`, `OptionValue`, `Options`, `validate`; `Format::options` for twenty readers.
- `formats::read_with`, `read_path_with`, `identify_with`; `Identification::preferred` and `by_preference`.
- `[pattern] reader_options = { … }` in the control file (`DataSource::reader_options`), applied by `control::open` (`data_scan_with`) and the checks.
- `drs data formats` lists each reader's options; `drs data identify --prefer`; `drs data read --option key=value --prefer`.
- Test `data_format_options`.

## 0.73.0 — 2026-10-21 — A refinement compared with GSAS-II (the app's request 50, L399)

The core prepares a GSAS-II run of a control file through GSASIIscriptable, reads what GSAS-II found in the core's terms and compares the two refinements, as it does with FullProf since 0.72.0.

### Added

- `compare::gsas2_prepare` and `gsas2_read` (module `compare::gsas2`, feature `serde`), with the fixed script `compare/gsas2.py`.
- `drs compare-with gsas2 prepare|read`.
- Test `gsas2_comparison` and the reference run `manuals/gsas2/runs/y2o3_tch`: a simulated TCH pattern fitted by both programs, χ² 1.1313 against 1.1333, 11 of 12 parameters within GSAS-II's esd.

### Known limitations

- Not carried to GSAS-II: the core's background (GSAS-II's Chebyshev polynomial instead), Bérar–Baldinozzi's asymmetry, preferred orientation, μR, a resolution file's widths; refused: a phase without atoms, time of flight.
- The scales are not compared: GSAS-II's phase fraction is on its own normalisation.
- GSAS-II's background-free R's and its R on F² are its own definitions, compared as printed.

## 0.72.0 — 2026-10-21 — A refinement compared with FullProf (the app's request 49, L398)

The core prepares a FullProf run of a control file, reads FullProf's outputs and compares the two refinements: χ², the R factors, every parameter with both esd's and Δ/σ, both calculated patterns.

### Added

- Module `compare`: `fullprof_prepare` (the run's folder, the command and a manifest of what was carried, approximated or refused), `fullprof_read` (the `Comparison`), `write` and `read` (`comparison.json`), `Run::save`, `Run::load`.
- `fullprof_outputs::read_out`: a FullProf `.out`'s banner, refined parameters with esd's, statistics of the last cycle, R_Bragg and R_F per phase, cycles and convergence.
- `pcr_export::Export::fixed`: the addresses of the free parameters the export writes fixed.
- `drs compare-with fullprof prepare|read`.
- Test `fullprof_comparison`: on FullProf's stored Y2O3 run 72 of 72 parameters within FullProf's esd; with `FULLPROF` set, FullProf run on the prepared folder.

### Fixed

- A `.pcr`'s files are found by the path it names when its folder does not hold them (`powder::find_file`): `background_terms` failed on a fresh build, where no other test had left a copy of its data.

### Known limitations

- R_F differs from FullProf's (Y2O3: 1.626 against 1.49) for a reason not found; frozen in the test.
- FullProf's halved esd's of a Thompson–Cox–Hastings phase's cell parameters are compared as printed, so those rows can read 'differs' at up to twice the rule; each such row says so.
- The run's folder name is the app's; it is not in the file dictionary yet.

## drs-publish 0.34.0 — 2026-10-20 — A series of refinements sent by Plots (L391)



### Added

- A series of refinements sent by Plots becomes a set with F5 or F4, each point's temperature with its source in the caption (PL34).

## app 0.52.1 — 2026-10-20 — A project's first data file never an export (L388)



### Fixed

- A project no longer takes one of the app's own exports as its first data file.

## drs-publish 0.33.0 — 2026-10-20 — TLS read; texture against ρ; relative paths; version 2 complete in Publish (L387)



### Added

- T16: each TLS group's T, L and S with su's (P47).

### Changed

- F10 against ρ (P59); relative paths in D2 (P61).

## app 0.52.0 — 2026-10-20 — External Programs checked; the Comparison view's frame (L382)



### Added

- Settings ▸ External Programs: FullProf and GSAS-II found and checked by running them.
- Refine's Comparison view, waiting for the core to run and compare (requests 49, 50).

## app 0.51.1 — 2026-10-20 — A background file left alone when unchanged (L379)



### Changed

- A background file is written only when its points changed, through the core; FullProf's own .bac stays byte for byte (request 48).

## app 0.51.0 — 2026-10-20 — TLS in Refine; the changed-on-disk question's buttons (L375)



### Added

- A rigid group's TLS: T, L and S shown and refined as 21 values named by the core, with their centre; each atom's displacement marked by its group's TLS (P41).

### Fixed

- Every question's buttons have a tooltip; the rigid item's internal coordinates no longer run into the TLS values.

## app 0.50.1 — 2026-10-20 — b1–b4 in Refine; every FullProf file named (L370)



### Added

- The scaled background's b1–b4 shown and refined (request 47).

### Fixed

- A background file with b1–b4 is no longer changed by each save.
- FullProf's empty .ort files named; the files the core cannot name fall from 15 to 10.

## drs-publish 0.32.0 — 2026-10-20 — The background file and b1–b4; a joint refinement's tables per pattern (L369)



### Added

- T13 gives a background file and its b1–b4 (P57); T13, T14, T17 of a joint refinement per pattern (P50).

## app 0.50.0 — 2026-10-19 — Texture along several axes in Refine (L366)



### Added

- Each texture axis's direction, March ratio and fraction shown and refined (P30).

## drs-publish 0.31.0 — 2026-10-19 — Texture along several axes read (L365)



### Added

- T14 along several axes; F10 the axes' March–Dollase pole densities (P48).

## app 0.49.0 — 2026-10-19 — A TOF refinement runs; further background terms in Refine (L362)



### Added

- A time-of-flight bank refined from Refine (request 46).
- Debye-like terms, the 12-term polynomial's inverse powers and backgrounds from files, shown and refined (P23).

### Known limitations

- Not installable until the core names FullProf's .ort files (request 47, L359).

## drs-publish 0.30.0 — 2026-10-19 — The background's parts read (L361)



### Added

- T13's further background terms with su's; F11 the background and its parts (P49).

## app 0.48.0 — 2026-10-19 — Rigid groups in Refine (L356)



### Added

- A phase's rigid groups shown and refined in place of their atoms' coordinates: origin, orientation, internal lengths and angles (P41).

## drs-publish 0.29.0 — 2026-10-19 — Rigid groups read; T15 (L355)



### Added

- T15: each rigid group's placement and its atoms' internal coordinates; T2's atoms in a group marked G.

## app 0.47.0 — 2026-10-19 — Time of flight in Refine; File ▸ Move To… (L350)



### Added

- A time-of-flight pattern in Refine against TOF, d or Q, its bank and back-to-back profile values refinable, the ticks in µs (P05).
- File ▸ Move To… moves a control file and rewrites its paths through the core (request 45).

### Fixed

- The self-check's Editor probe no longer auto-saved (L344).

### Known limitations

- A TOF refinement's Run is stopped by the core's checks until request 46.

## drs-publish 0.28.0 — 2026-10-19 — Time of flight read (L349)



### Added

- F1 against a bank's time of flight or d and Q; T1's bank; T17's TOF values with su's; D1 and D2 for TOF (P46).

## drs-publish 0.27.0 — 2026-10-19 — Each item's state in the Publish sidebar (L348)



### Added

- The sidebar says of each item built, stale (with its reason), failing or not built.

## 0.62.0 — 2026-10-19 — Rigid groups refined (P41) and written for Publish (P47)

L347: a rigid group's origin, orientation and internal coordinates refined in place of its atoms' coordinates with exact derivatives, against FullProf's refinements of PPH₃ and citric acid; groups kept in the control file, the `.pcr` export and the results file.

### Added

- `Phase::rigid`; `Field::Rigid`; `rigid_body::place`, `Group::value`, `set_value`, `derivative`, `cell_derivative`, `in_degrees`, `values_of`, `value_id`, `value_component`, `Coordinates::split`, `name`; `PcrAtom::rigid_codes`.
- Dictionary ids `rigid_origin`, `rigid_orientation`, `rigid_length`, `rigid_angle`.
- The control file's `[[phase.rigid_group]]` and an atom's `rigid_group`, `rigid_length`, `rigid_angle`, `rigid_connect`; the `.pcr` export of Jbt = 4 phases.
- The results file's `_drs_rigid_group_*` and `_drs_rigid_atom_*` loops and `_atom_site_refinement_flags_posn`.
- Test `rigid_refinement`; runs citric_core and citric_free named in `version2_runs.py`.

### Changed

- A Jbt = 4 phase's atoms are placed by their groups, and its codewords for θ, φ and line 25-4 are read (they were dropped).
- FullProf's ThRGB, PhRGB, ChRGB name a group's θ, φ, χ.

### Known limitations

- Citric acid's θ ends 5.2 esd from FullProf's (not explained, frozen); PPH₃'s angles 1.4–1.8 esd (frozen).
- A satellite's following of a turned parent is unchecked; TLS still waits for FullProf's runs.

## 0.63.0 — 2026-10-19 — Further background terms (P23) and their parts for Publish (P49)

L353: Debye-like terms, the 12-term polynomial with inverse powers, Nba = 1's b3 and b4 and profiles from files, calculated, refined and kept in the control file, the `.pcr` export and the results file, against FullProf's 16 runs of L346; FullProf's bounds on its background found and kept in its layer; the Fourier filter refused.

### Added

- Module `background_terms` (`BackgroundTerms`, `DebyeTerm`, `Profile`); `Model::background_terms`, `background_q`, `background_parts`, `background_unbounded`, `background_out_of_bounds`; `Field::BackgroundTerm`.
- `Convention::background_bounds`, `Model::background_ceiling`, `background_bounded`: FullProf's polynomial background 0 below 0 and above the data's largest count.
- Dictionary ids `background_debye_amplitude`, `background_debye_distance`, `background_inverse_power`, `background_profile`.
- The control file's `background_debye_amplitude`, `background_debye_distance`, `background_inverse_power`, `profile_files`, `background_profile`, `background_ceiling` and `[convention] background_bounds`; the `.pcr` export of Nba = −1, −3 and MORE = 2.
- The results file's `_drs_bkg_calc_smooth`, `_drs_bkg_calc_debye_n`, `_drs_bkg_calc_file_n`; `Results::background_parts`.
- `Pcr::profiles`, `profile_codes`. Test `background_terms`.

### Changed

- Nba = −1, −3, MORE = 2, 3 and Nba = 1 with b3, b4 are read and calculated (they were refused); MORE = 2 no longer reads a micro-absorption line.
- FullProf's Bck_12 … Bck_17 name the Debye-like amplitudes.
- L346's 42 runs counted in `fullprof_profiles`, `berar_factor`, `fullprof_outputs`, `output_names`, `control_file`.

### Fixed

- A polynomial background whose coefficients are all 0 is 0 whatever its origin (Bkpos = 0 gave NaN).
- The background-free Rp divides by Σ(y − b) as it is (manual p. 44), negative where the background lies above the counts; it was held above 10⁻³⁰⁰.

### Known limitations

- The Fourier filter (Nba = −2) is refused: not reproduced (L353).
- nac_3t2: 2θ = 0.05° 0.557 from FullProf's background, and the residual at FullProf's printed final values −0.64 % (frozen, for Josu).

## 0.64.0 — 2026-10-19 — A time-of-flight control file checked as it is calculated (the app's request 46)

L354: the checks Refine runs before a Run read a TOF bank's data on its time axis, need no wavelength, compare its range and expect its zero and excluded regions in µs; an xys TOF file is no longer read ×100.

### Added

- `control::data_scan`: a control file's data, a GSAS bank alone put on the µs axis. Test `fullprof_tof::tof_checks`.

### Changed

- The checks of a TOF control file: no wavelength; the data's range in µs; `zero_shift`, `excluded_regions` expected in µs; FullProf's TOF-min alone compared with the data.

### Fixed

- `control::open` scaled every TOF pattern's x by 100: alum's and lamn_pol's xys files (in µs) were read 100 times too far.

## 0.65.0 — 2026-10-19 — Texture along several axes (P30) and its numbers for Publish (P48)

L360: FullProf's Nor = 2 (several March–Dollase axes, each with its fraction) read, calculated, refined and kept in the control file, the `.pcr` export and the results file, against FullProf's texture runs within their print; Nor = 3 (the anvil cell) in the believed form of Wright et al. (1996), FullProf's up to 3.6 % below it, frozen.

### Added

- Module `texture` (`MultiAxial`, `Axis`, `march_dollase`, `of_phase`); `Phase::texture`; `Model::preferred_orientation_at`; `powder::texture_cosines`; `Field::Texture`.
- Dictionary ids `po_axis_ratio`, `po_axis_fraction`, `po_axis_steps`; the control file's `po_axis_direction`, `po_axis_ratio`, `po_axis_fraction`, `po_axis_steps`; the `.pcr` export of Nor = 2, 3.
- The results file's `_drs_pref_orient_function`, `_geometry`, `_steps`, `_index_h/k/l`, `_G1`, `_G2` and the loop `_drs_pref_orient_axis_*` with su's (P48.1–2).
- `drs texture <file>`: each phase's axes and March–Dollase's density at 0–90°.
- `PcrPhase::texture`, `texture_codes`. Test `texture`.

### Changed

- Nor = 2, 3 are read (they were refused); a Nor beyond 0–3 is refused by name.
- texture/'s four runs compared in `fullprof_profiles` (115); `field_labels` meets `Field::Texture` in yb_nor2.

### Known limitations

- Nor = 3: FullProf's yb_nor3 is up to 3.6 % below the core's cone integral at all 16 reflections; frozen in `texture` and `fullprof_profiles`, for Josu (L360).
- Spherical harmonics (P48.3) wait for a run that uses them.

## 0.66.0 — 2026-10-19 — FullProf's .ort named; a .bac background's b1–b4 as values (the app's request 47, Publish's P57)

L364: FullProf's ORTEP files named by the manual (p. 39); Nba = 1's file kept as read with b1–b4 as four refinable values of the background, in the control file, the .pcr export and the results file, with the background's further terms in its text.

### Added

- File kind `fullprof_ortep` (CODFILn.ort, an output).
- `background_terms::FileScale` (`BackgroundTerms::file`), `Part::File`; `Model::background_function_at`, `background_file_derivative`.
- Dictionary ids `background_file_scale`, `background_file_offset`, `background_file_stretch`, `background_file_shift`; the control file's keys of those names.
- The results file's `_drs_bkg_file`, `_drs_bkg_file_scale`, `_offset`, `_stretch`, `_shift` with su's (P57.2). Test `background_terms::file_scale_values`.

### Changed

- Nba = 1: the .bac's points kept as read and b1–b4 applied where the background is read (they moved the points before); its codewords read and written (they were refused).
- A control file from a .pcr with Nba = 1 names the .bac whatever b1–b4 are.
- `_pd_proc_ls_background_function` names the further terms (P57.1).
- `file_kinds.txt` regenerated over every run of the repository.

## 0.67.0 — 2026-10-19 — TLS displacements of rigid groups (P41, Publish's P47.2); the asymmetry at the Bragg angle

L367: a rigid group's T, L, S give its atoms' β as FullProf's do (hmt and pdd within the print of their .ort; urea off by up to 5.45·10⁻⁵ through a centre not found, frozen for Josu), refined, kept in the control file, the .pcr export and the results file; Bérar–Baldinozzi's asymmetry taken at the Bragg angle, as FullProf takes it.

### Added

- `rigid_body::Tls`, `Group::displacements`, `TLS_NAMES`, `tls_betas`, `betas`, `beta_derivative`, `read_tls`, `tls_of`, `centre_of_mass`, `tls_members`, `tls_centre`, `placement_derivative`; `scattering_table::ATOMIC_WEIGHT` (CrysFML's AtWe).
- Dictionary ids `rigid_tls_t`, `rigid_tls_l`, `rigid_tls_s`; the control file's keys of those names and `tls_centre`; `Model::bragg_two_theta`; `geometry::tls_u`.
- The results file's `_drs_rigid_tls_group_id`, `_centre_x/y/z`, the loops `_drs_rigid_tls_value_*` (T Å², L deg², S Å·deg with su's, the parameter and multiplier) and `_drs_rigid_tls_atom_*` (U with su's) (P47.2).
- `drs rigid` prints T, L, S, the centre and β. Tests `tls`, `waiting_for_runs::tls_read`.

### Changed

- P6 = 5.0x with Furth ≠ 0 is read and calculated (it was refused); the .pcr export writes Furth = 21 and the TLS lines.
- The Bérar–Baldinozzi asymmetry and AsyLim at the Bragg angle (without the zero and the sample's shifts), in the calculation and the Jacobian.
- FullProf's output names t11 … s32, X/Y/Z_RGB, Dm_In, TheIn, PhiIn and Or1 mapped to their fields.
- A TLS atom's β is fixed by "TLS"; the CIF's anisotropic U of TLS atoms carries su's through the covariance.

### Fixed

- The Jacobian of a rigid group's values and of the cell moves its satellites with it (urea's hydrogens stayed put).

### Known limitations

- Urea's TLS centre: FullProf's β fit z 0.33032–0.33036, its printed centre of mass is 0.326, CrysFML's weights give 0.32957; frozen (5.45·10⁻⁵ in β; urea_test's χ² 0.31 %), for Josu.
- FullProf's .pcr holds one TLS group per phase: the export refuses several.

## 0.68.0 — 2026-10-20 — A .bac with b1–b4 written as held (the app's request 48); D against ρ (P59); relative paths across the root (P61)

L374: a background file with b1–b4 is saved as its points are held and left alone when unchanged; March–Dollase's density against the angle for Publish's F10; paths that share only "/" written relative.

### Added

- `Control::background_file_to_write` (`None` when the file already holds the points).
- `texture::march_dollase_at`, `texture::march_dollase_curve`.
- Tests `background_terms::bac_with_b1_b4_written_as_held`, `paths_relative_across_the_root`; `texture::tests::the_curve_against_the_angle`.

### Changed

- `Control::background_file_text` writes a file's points as held when the background holds its b1–b4 (it wrote the background with b1–b4 at the data's points).
- `control::relative`: absolute only without a common root; the control file's paths and `_drs_bkg_file` relative to the control file across "/".
- `drs background store` leaves an unchanged file alone; `drs texture` samples D by `march_dollase_curve`.

## 0.69.0 — 2026-10-20 — Explore's E18–E21: the cell's su's and plausible ranges, the Williamson–Hall axes, a state at given values, the χ² map's contours

L376: what the Explore stream needs to restart, each with its test and command-line path.

### Added

- `PhaseState::cell_su` (E18a); `Contact::plausible` (E18b); `ValenceSum::expected`, `plausible`, `with_range`, `bonds::VALENCE_RANGE`, `VALENCE_RANGE_KEY`, `Phase::valence_sum_range` and the control file's `valence_sum_range` (E18c).
- `Breadths::hkl`, `d_star`, `beta_g_star`, `beta_l_star`, `beta_tch_star` (E19).
- `lsq::point`, `PointState` (E20a); `lsq::Contour`, `Grid::contours`, `contours_at_sigma`, `minimum`, `Grid::dof` (E21).
- `drs state` prints the cell's su's, the contacts' plausible distances and the valence sums' ranges, and `--at name=value,…`; `drs refine --map` prints the 1σ and 2σ contours. Test `explore_e18_e21`; `fullprof_microstructure::williamson_hall_axes`.

### Changed

- `chi2_grid`'s nodes are `lsq::point`'s: the links enforced as a step enforces them.
- `drs mic` prints d* from `Breadths::d_star`.

## 0.70.0 — 2026-10-20 — The Python bindings (P53): drs.run and drs.open, a wheel tested against the command line

L380, L381: the core's scripting interface in Python through PyO3 (Q1, L021). `drs.run` runs any command of the command line, whose binary the wheel carries; `drs.open` gives a project whose values, pattern, refinement, state, checks, export and results come from the core in process.

### Added

- Crate `core/python` (`drs-python`; package `diff-ref-struct`, module `drs`): `run`, `binary`, `open`/`Project`, `read_data`, `DrsError`, `__version__`; PyO3 0.23, `abi3-py39`.
- `core/python/build_wheel.sh` (the binary into the package, then maturin); `core/python/tests/test_drs.py` (9 cases against the command line); `core/python/README.md`.

### Known limitations

- Wheels built and tested here for Linux x86-64 on Python 3.13; the macOS wheel is built on the Mac with the same script; Python 3.9–3.12 not run.

## 0.71.0 — 2026-10-20 — A WinPLOTR header's temperature with its unit (Plots' PL32)

L383: the temperature line of a header that WinPLOTR wrote, read in K when its unit is known (the line's, or °C for a Bruker .UXD source, L243), never as kelvin by default; case E's t7.dat is 179.99 °C = 453.14 K.

### Added

- `formats::winplotr_temperature`; the registry gives its value to whichever reader takes the file's columns (xys, xye, table), its line as written beside it, or why it was not read.
- Test `winplotr_header`: case E through the three readers, a plot set (its source the file's header), a sequential refinement and `drs data read`; the units known and unknown.

### Fixed

- The xys reader no longer takes "File" of WinPLOTR's first line as the instrument.

## 0.61.0 — 2026-10-19 — A joint refinement's model per pattern from its results file (Publish's P50); rigid bodies placed as FullProf places them (P41)

L345: each pattern block of a joint results file holds its control file, read back without the data; module `rigid_body` reads and places FullProf's rigid groups against its files and prints; TLS, Nor = 2, 3, MORE = 2, 3 and Nba = −1, −2, −3 refused by name until FullProf's runs of them (`version2_runs.py`).

### Added

- `Joint::pattern_control`; `PatternResults::control_file`, `JointResults::control_file`, `name`, `pattern_controls`; `drs joint <results.cif> [--pattern k]`.
- Module `rigid_body`: `Group`, `Member`, `Satellite`, `Coordinates`, `groups`, `positions`, `rotation`, `satellite_axis`; `PcrAtom::rigid`; `drs rigid`.
- Tests `rigid_bodies`, `waiting_for_runs`, `joint_refinement::models_per_pattern_without_the_data`; `manuals/fullprof/runs/version2_runs.py`.

### Changed

- A joint results file's pattern blocks carry `_drs_control_file`.
- Nor ≥ 2 and MORE = 2, 3 are refused by the `.pcr` reader instead of misread; Nba = −1, −2, −3 and TLS are refused with the term named.

### Known limitations

- TLS, several preferred-orientation directions, the Debye-like, filtered, inverse-power and file-profile backgrounds wait for FullProf's runs (L345).
- A rigid group's values are not refined yet: the powder model takes the positions FullProf wrote.
- Four rigid-body files' written positions and two atoms of citric acid's Cartesian print depart from FullProf's own values: frozen and listed in L345.

## 0.60.0 — 2026-10-19 — A control file's paths rewritten when it is moved (the app's request 45)

L342: `control::relocate` and `drs move` rewrite every path a control file holds so that it reaches the same files from another folder, the rest of its text kept.

### Added

- `control::relocate`, `control::relocate::paths` and `PathEntry`; `drs move <file.drc> <folder>`; `Operation::file_mut`.
- Test `control_relocate`.

### Changed

- A joint refinement's patterns' control files carry the data operations' files relative to their folder.

### Known limitations

- A multi-bank GSAS file is not yet read as one pattern per bank (request 44).

## 0.59.0 — 2026-10-19 — Time-of-flight neutron patterns (P05), against FullProf's runs; Publish's P46

L341: a TOF bank and its phases calculated as FullProf does (Npr 9 and Hodges' Npr 10, the back-to-back profile, absorption, extinction), refined with analytic derivatives, kept in the control file and the `.pcr` export, and written to the results file as Publish reads it.

### Added

- Module `tof`: `Bank`, `Thermal`, `Profile`, `Instrument`, `Absorption`, `Parameters`, `back_to_back` with its partials, `extinction`, `time_axis`; `Bank::d_spacing`, `equation`; `Model::d_spacing_at`.
- `Field::Tof`, `Field::TofProfile`, their derivatives and 29 dictionary ids; the `.pcr`'s TOF codewords and symbolic names; the control file's TOF keys and `[pattern.resolution] function = "tof"`; the `.pcr` export of TOF patterns.
- The results file's `_pd_meas_time_of_flight`, `_pd_meas_2theta_fixed`, `_pd_calibration_conversion_eqn`, `_drs_tof_parameter_*`, `_drs_refln_tof`; `Results::axis`, `two_theta_fixed`, `conversion`, `tof_parameters`, `abscissa_in`; `PatternResults::axis`.
- Tests `fullprof_tof` (derivatives, refinement, results); FullProf's TOF runs in `fullprof_profiles` and the other tests that count FullProf's runs.

### Changed

- e^z E₁(z) by the continued fraction for Re z ≥ 0 from |z| = 2 (to 3·10⁻¹⁶, where the series lost up to 10⁻⁸).
- A TOF pattern's reflections generated only as far as the windows reach the data.
- `VARY back_n` of a background term its codeword refines is accepted (it changes nothing).

### Fixed

- The deposition CIF of a TOF pattern no longer generates reflections without bound (no wavelength).

### Known limitations

- FullProf's TOF peak window, the σ of a RALF bank, Iwg = 1 in TOF and the esd of Extinc are not reproduced: frozen and listed in L341.
- ABSCOR2, the TOF size and strain terms and Npr 13 are refused; magnetic phases (Ho2BaCuO5's four banks) are refused.
- The flight path is not written: no file gives it.

## app 0.46.0 — 2026-10-19 — One Open… for folders and files; Export as CIF; Duplicate; Report an Issue; thin-line icons; smart zoom (L340)



### Added

- File ▸ Open… (⌘O) opens a folder or a file; File ▸ Export ▸ CIF… writes Publish's CIF; File ▸ Duplicate copies a control file beside it.
- Help ▸ Report an Issue… opens Mail to Josu with the versions and the diagnostics report.
- Settings ▸ General ▸ At launch, open; the Editor's font family; two-finger double-tap fits the plot.

### Changed

- Every icon in thin lines (L337).

### Fixed

- Edit ▸ Emoji & Symbols and Start Dictation no longer repeat; Settings ▸ Shortcuts' Change works after a relaunch.

### Known limitations

- File ▸ Move To… waits for the core to rewrite a control file's paths (request 45).

## drs-publish 0.26.0 — 2026-10-19 — Distances, angles and bond-valence sums against temperature or composition (L338)



### Added

- F7 draws a series' bond distances, angles, mean distances and bond-valence sums from the core's results files, with su's.
- T4 takes ⟨d⟩ and the bond-valence sums with their su's from the results file.

## app 0.45.1 — 2026-10-19 — The self-check's Copy checks read back the app's own write (L336)



### Fixed

- The self-check's Copy lines can no longer be failed by another program copying meanwhile (L330).

## drs-publish 0.25.0 — 2026-10-19 — The set's changes on the shared Undo; instrument and profile parameters (L335)



### Added

- Each change to a Publish set is one named step of Edit ▸ Undo and Redo.
- T17: every instrument and profile parameter the core lists, with its symbol, unit and su.

## app 0.45.0 — 2026-10-19 — Undo as a shared service; a joint refinement saved as our own control files (L333)



### Added

- Edit ▸ Undo and Redo name the change, one stack per document, used by Refine, the Editor, Data, a joint refinement and Explore.
- A multi-pattern .pcr imported into our own joint control file and one per pattern; its weights, flags and contributions set through the core (request 43).

## drs-publish 0.24.0 — 2026-10-19 — The rest of version 2 designed; background and preferred-orientation tables (L332)



### Added

- T13, the background: five models, each point or term with its su; T14, the preferred orientation: March–Dollase or the exponential, the axis, G1 and G2 with su's.
- Designed, to be filled when the core gives them: time of flight (F1, T1, T17, D2), rigid bodies (T15) and TLS (T16), texture (F10), background contributions (F11).

## app 0.44.0 — 2026-10-18 — Joint refinement in Refine; the shared table everywhere (L328)



### Added

- A multi-pattern .pcr opens in Refine as a joint refinement: each pattern's weight, flags, plot and agreement, the joint χ², run and accepted as the core's results file (P06).
- Refine's and Data's tables are the shared table: sort, choose columns, find, copy; the columns chosen are remembered.

### Changed

- The parametric fit uses the core's exact derivatives (request 42).

### Known limitations

- A joint refinement's settings are not yet written back: the .pcr stays as it is until the core gives a control file for several patterns (request 43).
- Installed with one self-check line failing (L330): the copied plot's clipboard read back as text, likely the Mac's shared clipboard during the run.

## drs-publish 0.23.0 — 2026-10-18 — A joint refinement read; the core's volume su; F1 of every pattern (L327)



### Added

- F1, T1 and D2 read a joint refinement a pattern at a time; F1 of every pattern, a panel each.

### Changed

- The CIF's cell-volume su is the core's one (P43).

## app 0.43.0 — 2026-10-18 — The Pawley fit and a series of patterns in Refine (L324)



### Added

- Refine ▸ Pawley Fit for Phase…: reflection intensities refined with the cell and profile, the coincident reflections tied (P36).
- Refine ▸ Refine a Series…, the Series view with su bars, a parametric fit to an equation you write, Series as CSV and CIF (P54, P59).

### Fixed

- A relaunched tab group comes back whole (L321).

### Known limitations

- Plots' own check of its saved set fails on a relaunch: Plots' code.

## drs-publish 0.22.0 — 2026-10-18 — The core's Pawley fit and sequential series read; parametric fits (L323)



### Added

- T10, T11, F9, D9 read a Pawley fit, its tied reflections counted and marked.
- F4, F5, F7, T7, D8 read the core's sequential series; P54's parametric fits as F7's line and T12.

## 0.58.0 — 2026-10-18 — A joint refinement's control file and settings (app request 43)

L329: a joint refinement saved as a joint control file naming one control file per pattern, read back exactly; its weights, flags and phases' contributions set through the core.

### Added

- `Joint::save`, `control_files`, `open`, `set_weight`, `set_weights`, `set_refined`, `set_contribution`, `objective_of_pattern`; `joint::is_joint_control`; `Member`'s title, refinement, data source, FullProf record and background file.
- `Field::with_phase`, `Field::is_structural`.
- The file dictionary's `pattern_control`; `drs joint --save`; `drs joint`, `drs refine` and `drs results` on a joint control file.
- Test `joint_control`.

## 0.57.1 — 2026-10-18 — A parametric fit's coefficient at 0 (app 42); one cell-volume su (P43)

L326: exact derivatives in parametric fits; the cell volume's su through the covariance everywhere the core gives it.

### Added

- `Equation::derivative`; `ParameterSet::cell_volume`; `Pcr::leave_out_unmodelled`; `drs sequential --leave-out-unmodelled`; `PhaseLine::volume_su`.

### Fixed

- A parametric fit whose coefficient fits to 0 no longer stops as singular.
- `_cell_volume`'s su in the results file and the deposition CIF is the covariance's, as the series table's.

## app 0.42.1 — 2026-10-18 — What kept 0.41.1 and 0.42.0 from being installed (L319)



### Fixed

- Edit ▸ Spelling and Grammar's actions reach a window behind another.
- View ▸ Bottom Panel on the Right after a relaunch; ⧉ floating a group of tabs.
- Stop during Run the Plan to the End ends the plan, between steps too.
- A relaunch that hung building its windows after a full-screen tab group: a tab group comes back without full screen.

### Known limitations

- Plots' own check of its saved set fails on a relaunch: Plots' code, left to its stream.
- Installed on Josu's Mac with one FAIL (L321): after a relaunch, a macOS tab group of three windows came back as two.

## drs-publish 0.21.0 — 2026-10-18 — A series of refinements; the guided run by its steps (L318)



### Added

- F7 from refinements: refined parameters against the series' variable, with su's; T7 of parameters, a column per point.
- T9 by the plan steps the core records; example y2o3-guided-tfg.

## 0.57.0 — 2026-10-18 — Joint refinement (P06) and its results file (P41)

L322: several patterns refined together, from FullProf's multi-pattern files, against FullProf's joint runs pbb, scor_rxn and prbaco; a results file with a block per pattern, as Publish reads it.

### Added

- Module `joint` (`Joint`, `refine`, `Esd`, `statistics`, `results`, `read_results`, `fullprof_field`); `drs joint`, and `drs refine`/`drs results` on a multi-pattern `.pcr`.
- `pcr::split_patterns`, `pcr::read_all`, `pcr::is_multipattern`; `refinable::pcr_numbered_set`.
- FullProf's strain model 8 (hexagonal), read, refined, written back and named.
- Test `joint_refinement`; the joint runs in `fullprof_profiles`, `berar_factor` and `output_names`.

### Fixed

- Synchrotron polarisation (Ilo = 3) with Cthm = 0: cos²2θ_M = 1, as FullProf.
- A point one step before an excluded region's edge is kept, as FullProf keeps it, to the data's rounding.

### Known limitations

- pbb's and scor_rxn's minima below FullProf's ends, prbaco's neutron Σ 0.13% from FullProf's: frozen with their reasons (L322), for Josu.
- No Le Bail phases, restraints, steering or plans in a joint refinement yet; no control file of our own for several patterns.

## 0.56.0 — 2026-10-18 — Sequential and parametric refinement (P54); a series' table and CIF (P59)

L317: a series of patterns refined in turn, as FullProf's sequential mode or each from the control file, against FullProf's sequential D20 run of LaMnO₃; refined values fitted to equations of temperature, against SciPy.

### Added

- Module `sequential` (`read`, `run`, `Start`, `cell_volume`, `series_csv`, `series_cif`, `split_series`); `drs sequential`.
- Module `parametric` (`Equation`, `fit`, `column`); `drs parametric`.
- The file dictionary's series_table and series_cif.

### Known limitations

- The last background point of the D20 run: 1/√2 of FullProf's √C and up to 1.5σ from its value (L317), for Josu.
- A magnetic phase is not modelled: the D20 run is compared above its Néel point only.

## 0.55.0 — 2026-10-18 — The Pawley fit (P36)

L315: a phase's reflection intensities refined by least squares with the cell and profile, the coincident reflections tied with equal |F|²; against FullProf's Le Bail fits and the spinel's own pattern; the results file as Publish's P39 reads it.

### Added

- `Field::Intensity` (`phase1.reflection_intensity.h_k_l`); module `pawley` (`coincident`, `parameters`, `setup`, `ties`); `drs pawley <file> --to <file.drc>`.
- The results file: a Pawley phase's I meas and σ the refined intensity and its su; `_drs_pawley_tie_*`, read back as `Results::pawley_ties`.
- `lsq::householder_ql` and `JACOBI_LIMIT`: Householder + QL above 120 parameters.

### Changed

- The dictionary's `reflection_intensity`: refined by a Pawley fit, ≥ 0; reference Pawley (1981).
- A Pawley phase's intensities are fields of the model (`Model::fields`); a Le Bail phase's are not.

### Known limitations

- Accidental overlaps are not tied (GSAS-II ties d within 10⁻⁴ in −1, 2/m, mmm): left out of each cycle.
- The core's Le Bail partitions from FullProf's final intensities of Cocuc_LeBail_conv move them (median 0.97σ, largest 6.4σ): for Josu (L315).

## app 0.42.0 — 2026-10-16 — State restoration: tab groups, panel windows, full screen and split sizes (L314)



### Added

- macOS tab groups rejoined at launch with the same tab in front; panels in windows of their own back at their frames; full screen back on its tab; a split's groups at their kept sizes.

### Known limitations

- L312's six FAILs not yet addressed; a relaunch self-check run once ended early through the quit path, cause not found; full screen's restoration not yet confirmed in the app.

## drs-publish 0.20.0 — 2026-10-16 — A whole-pattern decomposition's tables and figure (L313)



### Added

- T10 the fit, T11 the extracted intensities with their clusters, F9 the intensities as sticks, D9 FullProf's .int; from the core's Le Bail fit, ready for Pawley's (P36).

## 0.54.0 — 2026-10-16 — A plan run recorded as one, each cycle with its step (Publish's P37)

L310: a control file's plan of steps run as one refinement, its cycles numbered on over the steps and each recorded with the step that ran it, so that Publish reads the strategy from the core.

### Added

- `lsq::run_plan` and `PlanStepRun`; `lsq::Refinement::cycle_steps`, each cycle's step number and name.
- The results file's `_drs_cycle_step` and `_drs_cycle_step_name`, read back as `Results::cycle_steps`.
- `drs results <file.drc> --plan`.

### Changed

- `drs refine --plan` runs through `run_plan`: its cycle count is the run's, no longer the last step's.

## app 0.41.1 — 2026-10-16 — 0.40.0's three FAILs; request 41 used (L309)



### Changed

- The Editor finds a parameter's line in a control file whose background points are in a file of their own (request 41).

### Fixed

- View ▸ Bottom Panel on the Right after a relaunch: only the panel and its Log move, to the left of the inspector.
- The diagnostics report leaves every path out of its log lines.

### Known limitations

- Explore's own services check reads 0 of 0 parameters on a relaunch: Explore's code, left to its stream.

## drs-publish 0.19.0 — 2026-10-16 — The refinement's strategy as run, and a series' CIF and CSV (L308)



### Added

- T9: the refinement's strategy as run, its stages from the recorded cycles, named by the control file's plan when it matches; in the supplementary material too (strategy.csv).
- D8: a series' multi-block CIF of every point, and a CSV of every point's R factors, χ² and refined values with su's.

## app 0.41.0 — 2026-10-16 — P55 in Refine: the core's guided refinement as the plan (L306)



### Added

- Refine ▸ Guided Plan…: the documented order of refinement from the core, each step's source and the values it adds and drops by name, written as the control file's plan.
- Run the Plan to the End: each step run and accepted, a version each, stopping where the core's plan stops.

### Changed

- The app's own demonstration plan replaced by the core's.

### Known limitations

- 0.40.0's three FAILs (L307): the bottom panel on the right, the inspector's ? with Refine in front, the diagnostics report copied whole.

## 0.53.0 — 2026-10-16 — locate with the control file's folder; every value's place from one reading (the app's request 41)

L305: a control file read with its folder, so that one keeping its background points apart has its model; every value's place, and the value at a byte, from one reading.

### Added

- `control::locate::locate_in(text, address, dir)`, `all(text, dir)`, `address_at(text, byte, dir)`.
- `drs locate <file.drc> [<address> | --at <byte>]`, read with the file's folder; without an address, every value.

### Fixed

- A control file with its background points in a file of their own, and refined, was located nowhere (no model without its folder; L303).

## drs-publish 0.18.0 — 2026-10-16 — The command line's commands for scripts (L304)



### Added

- drs-publish new --control <file.drc>: a set with the default items of the newest results file.
- check --json; preview <set.drp> <item>: one item made alone, the same bytes as in the build; figure|table|file without a set.

## app 0.40.0 — 2026-10-16 — A parameter's line in the control file; Refine's problems marked by their value; the audit's small conventions (L303)



### Added

- Jump to Parameter and the linked selection between Refine and the Editor, from the core's control::locate (request 40).
- Refine's problems marked on the outline and the inspector by their address (§9).
- The bottom panel docks on the right; the modules' order changes by dragging; Paste and Match Style; the standard menu items in ⌘K and the context-menu commands in a menu; ? opens the page for the panel in front; Increase contrast; Advanced's diagnostics; Plots' marker size.

### Known limitations

- A control file whose background points are in a file of their own (as drs new writes it) is not located until the core reads it with its folder (request 41).

## 0.52.0 — 2026-10-16 — P55: guided refinement, the documented order of refinement as a plan of steps

L302: FullProf's manual's sequence of refinement (p. 14), completed by McCusker et al. 1999 §7, offered for a model as a plan of steps, each naming the values it frees; written to the control file and run.

### Added

- `guide::documented_plan(&Model) -> Vec<GuidedStep>` (the step, its source, what it adds and drops) and `guide::guided(&Model, &Refinement) -> Refinement`.
- `drs plan <file.pcr|file.drc> [--to <file.drc>]`; the plan runs with `drs refine <file.drc> --plan`.
- Test `guided_refinement`.

### Known limitations

- Preferred orientation, the wavelength and the micro-absorption are not in the manual's list: they are refined only in the last step, and only when the file refines them already.
- The overall B is kept at its refined value once the atoms' B's are free (equivalent to the transfer the manual describes), not moved into them.

## app 0.39.1 — 2026-10-16 — Parameters named by the core; the document types read by their keys (L300)



### Changed

- Each parameter's phase, atom and component come from the core's Field::describe: the app keeps no list of the core's fields (L276, request 39).
- The document types are read from Info.plist by their keys, whatever the keys' order (L296).

### Known limitations

- APP_CONVENTIONS.md: 26 items unmet or partly met, listed with their sizes in app/README.md.

## 0.51.0 — 2026-10-16 — control::locate: where a control file writes a parameter (the app's request 40)

L299: for a control file's text and an address of its model, the line and the bytes of the value's text, from the reader's syntax tree; for the Editor's linked selection.

### Added

- `control::locate::locate(text, address) -> Option<(usize, Range<usize>)>`; `drs locate <file.drc> <address>`.
- Test `control_locate`: 620 values of eight runs, each edited where located and read back.

## drs-publish 0.17.0 — 2026-10-16 — Provenance in each file of a build, read back

The next largest item DESIGN.md planned that needs no core work (§5.3, §27): every file a build writes says where it came from, and drs-publish reads it back.

### Added

- Each file carries its item's header lines (the crate's version, the set, the item, the preset; each source's name and SHA-256; the core that wrote it) where its format keeps text: the PDF's XMP packet, comments in .tex, .pgf, .cif, .md, .html, SVG and EPS, a PNG's iTXt, a TIFF's ImageDescription, a .docx's core properties; CSV, .xye, .vesta and the Bilbao .txt keep theirs in the manifest.
- drs-publish provenance <file>: the file's provenance, the manifest's sources, and whether the file is unchanged (its SHA-256); shown in the app's inspector (publish_provenance).

## app 0.39.0 — 2026-10-16 — Explore and Plots on the shell's shared services; the parameters' labels from the core; on core 0.50.0 (L292–L296)



### Changed

- Explore draws its plots, 3D view, time slider, tables and linked selection with the shell's shared services, and has its own Settings pane (E8, E8a).
- Plots draws with the shared plot: measured margins, a figure's text in points (PL13, PL30).
- On core 0.50.0: Field::describe gives each refinable value's phase, atom and component for display (request 39).

## app 0.38.1 — 2026-10-16 — The two FAILs of 0.38.0's installed check (L290)



### Fixed

- Refine's phase name field named for VoiceOver and explained by a tooltip, whatever its value.
- A printed figure's width capped from its own proportions as well as its height, so that a landscape print stays on one page; the figure never broken across pages.

### Known limitations

- Field::describe (request 39) and an address's line in a control file (request 40) wait for the core.

## drs-publish 0.16.0 — 2026-10-15 — F1's options; the decimal marker follows the language

Josu's L284 and L285: a decimal point everywhere, Spanish included; Basque keeps the comma. Then the next largest item DESIGN.md planned that needs no core work: F1's options (§2.1, §26).

### Added

- F1's y_scale (sqrt, log: plots::YScale), x_axis (d, q: data_ops::to_axis with λ₁), difference_in = sigma ((Yobs − Ycalc)/σ in a strip, diagnostics::cumulative's weighted difference), show_background, cumulative (Cumulative::gsas2) and magnify_from = [2θ, n]; in the inspector; the caption names each, in the set's language.
- Check 4: every figure and table in its language's decimal marker, the comma in Basque, the point elsewhere.
- Three figures in the example y2o3-tfg (rietveld-sqrt, rietveld-d, rietveld-log).

### Changed

- The decimal marker follows the language (L285): Basque's comma, the point in Spanish and every other; a decimal key may only say the language's own; Choose Language… offers the language alone.

## 0.50.0 — 2026-10-15 — Field::describe: each refinable value labelled by the core (the app's request 39)

L292: a field's dictionary id and name, phase, atom and component, so that the app labels every parameter without a list of the core's fields.

### Added

- `Field::describe(&Model) -> refinable::FieldLabel` (id, name, phase, atom (index, label), component).
- Test `field_labels`: every field of eight FullProf runs, all 29 kinds, matched without a wildcard.

## app 0.38.0 — 2026-10-15 — The shell's services for Explore and Plots; project windows as macOS tabs; the document icons' family



### Added

- plot.js for Explore (E8, E8a) and Plots (PL13, PL30): numbered paths, filled patches, contour lines over a map, a symmetric-log scale, hatched spans, areas filled below a curve, linked axes, hundreds of series drawn per pixel column, text in points with measured margins, log ticks at 2 and 5, a figure's factor in its unit.
- Shared services: the table's column of names, a time slider (Refine's run on it first), linked selection of a parameter between views, a module's own Settings pane, a 3D canvas on three.js.
- Project windows join macOS tab groups (Window ▸ Merge All Windows, Move Tab to New Window, Show Tab Bar, Show All Tabs); the Settings and the panels' windows never; the page moves under the tab bar (L281).
- The document icons in one family with the app's icon: braces, stacked curves, a page with a figure, a pattern with ticks (L281).

### Fixed

- A split kept: a Refine document saved from the Editor is read again in place, its tab staying where it was (L283).
- The self-check's print to a file never asks to replace an earlier one (L283).

### Known limitations

- Field::describe (request 39) and an address's line in the control file (request 40) are asked of the core; the modules move onto the services in their own steps.

## drs-publish 0.15.0 — 2026-10-15 — EHU Sans as installed on the Mac; a thesis in Spanish or Basque, with the decimal comma

Josu installed the EHU fonts (L274): a TFG is set in EHU Sans Light as the Faculty's model has it. Then the largest item DESIGN.md still planned that needs no core work: the thesis's language (§5.1, §25).

### Added

- language = "en" | "es" | "eu" in a set: every label, legend, caption, table head, row label and note in that language; strings.yml holds the 169 texts' Spanish and Basque, written into src/strings_table.rs by tools/make_strings_table.py (--check, --add); a text untranslated stays English and the build says so.
- decimal = "comma": the comma in the figures' maths and through siunitx's marker, semicolons between numbers in a list, the range's word in the language; the CIF, the CSV and the data files keep the point.
- The example sr2scsbo6-eu-tfg.drp (case A, a TFG in Basque with the comma); the app's Choose Language… and the inspector's Language row.

### Changed

- A face not in TeX Live is looked for installed on the Mac first (~/Library/Fonts, /Library/Fonts), loaded by name (EHUSans-Light upright, Bold, Italic, BoldItalic), the project's Fonts/ kept as the fallback; the manifest says where it came from.
- The notes on the amplitude set reported (T5, T6, F6) are written as whole sentences, so that they translate.

### Known limitations

- The translations are a first version, for Josu to read. The core's own words quoted in a note (the operation to the reference's setting) stay English.
- EHU Sans has no Light Italic: the maths' italics are Regular-weight Italic beside Light text, as in the Faculty's model.

## app 0.37.0 — 2026-10-15 — Panels in macOS windows of their own; a parameter's phase from the core (L276)



### Added

- The sidebar, the inspector and the bottom panel detached into macOS windows of their own, following their window's file and selection, one Refine document shared through the core, docked back by closing them (APP_CONVENTIONS §1, §3).

### Changed

- A parameter's label takes its phase from the core's Field::phase(): the app no longer lists the core's fields, so a new one does not break the build (L276); the component's name waits for Field::describe (request 39).

### Fixed

- P45's numbers kept while the document's pattern is the same; the inspector shows a phase's parameters' impact.
- A split made in the window is kept; the settings a window writes are not read back over its own (Default Set); the Plots tests build the core's drs once.

### Known limitations

- Native window tabbing stays Josu's decision; the tabs-kept-at-opening fix of this release was not run in the app before the merge.

## 0.49.0 — 2026-10-15 — P13 closed: the ISODISTORT reader checked on ISODISTORT's own CIF; the mode matrix read the right way round

L278: ISODISTORT 6.12.2's mode CIF of case A read, and its 7 irreps and 12 modes found by the core's own decomposition of the same structures, each amplitude to ISODISTORT's five decimals. The mode matrix gives the deltas from the modes.

### Added

- `IsoDistort::deltas_from_modes` (δ = M·(A × N)) and `modes_from_deltas`.
- `drs isodistort` prints each irrep's amplitude and each delta beside the one the modes give.
- Test `isodistort`: ISODISTORT's CIF of case A against the file's deltas and the core's decomposition.

### Changed

- `IsoDistort::matrix_times_deltas` removed: it applied the matrix to the deltas, the reverse of ISODISTORT's definition (L263's reading).

### Known limitations

- ISODISTORT's strain modes (`_iso_strainmode_*`) are not read.
- The core's irrep names stay computed (GM+[3]a); ISODISTORT's labels are matched in the test, not carried onto the decomposition.

## app 0.36.0 — 2026-10-15 — P45 in Refine: each parameter's impact and the normal probability plot; the Publish set and Fit Peak's Kα doublet from the core (L277)



### Added

- Refine: each parameter's impact on χ² and on the calculated pattern, in the inspector and ranked under Refine ▸ Parameter Impact; the normal probability plot of the weighted residuals (core 0.48.0, P45).
- Fit Peak fits the Kα doublet from the scan's wavelengths (core 0.47.0, request 38).

### Changed

- A project's publication set is found through the core's file dictionary (request 37).

## drs-publish 0.14.0 — 2026-10-15 — Plots dropped into Publish; T5 and F6 from the modes the core computes; EHU Sans's files found whatever their case

App 0.34.0's drag and drop (L265) and core 0.46.0's mode decomposition (L263) taken in.

### Added

- module::drop_plot: a Refine plot a Rietveld plot of its newest results file, a Data plot a pattern figure of a one-file plot set the core writes; the app's view and sidebar drop targets (dropTarget).
- source::computed: the modes of a phase refined by coordinates with no AMPLIMODES file, against the parent the set names ([[phase]] parent, parent_basis, parent_origin); T5 and F6 from them, the irreps named as the core names them.
- The example sr2scsbo6-computed-tfg.drp and its Fm-3m parent CIF.

### Changed

- EHU Sans's files found whatever their case (EHUSans-Light.otf for EHUSans-light.otf).

### Fixed

- A project face's files met in another folder than the one luaotfload cached them from were looked for there: its stale cache entries are removed before the build (Tex::forget_other_fonts).

### Known limitations

- EHU Sans's files not yet in a project: no TFG set in it so far.

## 0.48.0 — 2026-10-15 — P45: each parameter's impact on χ² and on the pattern; the normal probability plot of the weighted residuals

L272: refinement diagnostics as GSAS-II defines them. Every parameter, refined or not, is weighed by what refining it alone would do: its shift in its own esd's, the fall of χ², and where it moves the calculated pattern. The normal probability plot uses the exact expected normal order statistics.

### Added

- `impact`: `impacts`, `impact_of`, `Impact` (slope ∂χ²/∂p, shift, esds, fall of χ², the largest change and its 2θ, the change for one esd, the ranges, the change point by point), `normal_probability`, `NormalProbability`, `central_line`, `expected_normal_order_statistics`, `ln_normal_cdf`, `normal_quantile`.
- `drs impact <file> [--cycles n] [--parameter a] [--top k] [--json]` and `drs normal <file> [--cycles n] [--points] [--json]`.
- Test `refinement_diagnostics`; `tests/data/expected_normal_order_statistics.py` (mpmath).

### Known limitations

- Toby (2024) and Abrahams & Keve (1971) are not in `manuals/papers/`: the definitions are GSAS-II's code's.
- The restraints' terms are not in the χ² the impact weighs.

## app 0.35.0 — 2026-10-15 — Landscape printing on one page; Spelling and Grammar; Customize Toolbar…; the tabs' split groups restored



### Added

- Edit ▸ Spelling and Grammar with the Mac's own items and selectors (Show Spelling and Grammar ⌘:, Check Document Now ⌘;, Check Spelling While Typing, grammar, autocorrection), answered by WebKit for the Editor's text (APP_CONVENTIONS §3).
- View ▸ Customize Toolbar…: the visible module's toolbar chosen and ordered, Default Set, kept in the settings for every window (§1, §3).
- The tabs' split groups (Split Right, Split Down) kept with the project and split again when its window reopens (P46).

### Fixed

- Print… in landscape gives one page: a figure is kept within the page's height less its margins, read from Page Setup (L266).
- A split made in the window is kept at once, from dockview's own events.

### Known limitations

- The split groups come back at equal sizes; requests 37 and 38 wait for the core.

## drs-publish 0.13.0 — 2026-10-15 — Word manuscripts for the journals that take Word (L257)

Acta C, IUCrData and Powder Diffraction: the tables and the text parts as .docx beside the LaTeX, the figures in each journal's own formats.

### Added

- docx.rs: a WordprocessingML writer of the crate's own (stored ZIP, CRC-32, no dates): paragraphs, runs from the tables' LaTeX, booktabs-ruled tables, inline PNG figures.
- presets.yml word: the IUCr's (A4, figures embedded, arabic) and Powder Diffraction's (letter, double-spaced, Times New Roman 12 pt, Symbol-font Greek, figures apart, roman), cited; each table's .docx and {set} manuscript.docx.
- Powder Diffraction's TIFF in the ICDD's colour modes (CMYK, 1-bit) by Ghostscript, its DateTime fixed for a reproducible build; check 4 for the Word files; the app's inspector names the manuscript.

### Changed

- A unit after a number: no space before ° and % (the HTML, Markdown and CSV exports too).

### Known limitations

- The Symbol-font Greek of Powder Diffraction's Word files is unchecked in Word itself (QuickLook does not render it).

## 0.47.0 — 2026-10-15 — P25's rest against FullProf: needles, antiphase domains, the quartic form of every Laue class; FullProf's Sycos·cos 2θ; the app's requests 37 and 38

L266: the anisotropic size (platelets, needles, Cu₃Au's antiphase domains) and the quartic strain of any Laue class, checked reflection by reflection on FullProf's runs of prba10k, lac, cuf1k and Pd3MnD; FullProf's Ste read; FullProf found to shift by Sycos·cos 2θ, not its manual's cos θ, a compatibility rule; the Publish module's set file in the file dictionary and the Kα doublet from a scan's wavelengths.

### Added

- `microstructure::ShapeSize` and `SizeShape` (Platelet, Needle, Antiphase): FullProf's Size-Models 1, −1 and −14, LorSiz's width k_s λ LorSiz g/cos θ; `Phase::shape_size`, the field `ShapeSize` (dictionary `shape_size`), read from and written to the `.pcr` (LINE 27, LINE 36) and the control file (`shape_size`, `shape_size_kind`, `shape_size_direction`).
- `microstructure::quartic_forms` and `QuarticStrain::for_group`: the quartic form of any Laue class from its rotations' orbit sums; FullProf's Strain-Models 4 and 5 (4/m, 4/mmm) read and written.
- `Convention::displacement` (the control file's `displacement`): FullProf's Sycos·cos 2θ; `Model::displacement_factor`.
- The data operation `rebin` with mode `keep` (every n-th point), FullProf's Ste, on import and back on export.
- `output_names::fullprof_fields`: every atom a FullProf name can mean when atoms share a label; `StrainN`, `LorStr` and `Ani-LSiz` named.
- `peakfit::doublet_of`, `PeakRequest::with_wavelength`, `drs peak … --doublet scan` (the app's request 38).
- The file kind `publication_set`, `Publish/{name}.drp` (the app's request 37), and the placeholder `{name}`.
- Tests `fullprof_anisotropic_size` (325 reflections of four runs; the files and the engine) and `fullprof_displacement`.

### Changed

- The peak window under FullProf's rules includes the anisotropic size and leaves out the Lorentzian strain ξ·D_ST tan θ, with the extra low point when X ≠ 0 or ξ ≠ 0 (L266).
- `quartic_terms` no longer refuses the tetragonal, trigonal, hexagonal and cubic classes.
- prba10k and Pd3MnD are compared in `fullprof_profiles` (70 files); their edge points, prba10k's Bérar factor and Pd3MnD's Durbin–Watson d are frozen; the counting tests count the new runs' files.

### Fixed

- A `.pcr` with Sol ≠ 0 had its shift line read as the platelets' vector; such phases are now refused (P26).
- The `.hkl` ticks of NdNi-HTO and cuf1k (Sycos ≠ 0), left unresolved in L098, agree with FullProf's `.prf`.

### Known limitations

- Size by the quadratic form (Str = −1, 2, 3), spherical harmonics (Size-Models 15–22) and the selective rules other than −14 are not calculated; hkl-dependent peak shifts (Sol) are P26.
- The trigonal, hexagonal and cubic quartic forms follow from the rotations; no FullProf run checks them yet.
- NdNi-HTO's nuclear profile against FullProf's `.sub`: 203 of 2990 points beyond its print, by up to 1.95 counts, at peak-window edges; frozen.

## app 0.34.0 — 2026-10-15 — Drag and drop between modules; every control named for VoiceOver and explained; the installed check's three FAILs



### Added

- Drag and drop between modules, a shared service (APP_CONVENTIONS §7, §11): files from the sidebar onto a module's icon, Data (drawn over), Refine (a data file starts a refinement from its control file; a CIF adds its phase) and the Editor; a plot by its grip into a module that takes plots; what a drop will do said in words over the target.
- Every control named for VoiceOver and explained by a tooltip (§5, §10), a shared service: the dialogs' controls with their own words, the module templates' controls written, plots described in words; audited on five pages in a headless browser and in each module of the running app.

### Fixed

- Copy as Image goes to the pasteboard through AppKit, so a menu chosen while the window is behind still copies (L264).
- Print: the plot takes the page's width at its proportions (a viewBox from its own size), with no trailing blank page; the Editor's print takes its own tab.
- The self-check fits a peak the core finds in the open file, and runs Finder, Print and the drops before its long waits (L264).

### Known limitations

- Publish's target for a plot (§7) is the Publish stream's, with the service's two calls; requests 37 and 38 wait for the core.

## 0.46.0 — 2026-10-15 — Version 1's step D: symmetry modes computed from the structures; ISODISTORT's mode CIF read

L263: P13's modes computed from a parent and a distorted structure, as AMPLIMODES computes them, the irreps' k, parity, dimension, modes and amplitudes from the class sums of the parent's group without tables, checked on the course's ten AMPLIMODES analyses; ISODISTORT's mode CIF read to the tags GSAS-II reads.

### Added

- `mode_decomposition`: `decompose`, `Decomposition` (reference structure, displacements, irreps with star, parity, dimension, reality, modes and amplitudes per orbit), `name_by_modes`, `to_modes`, `fullprof_block`, `field`, `outside`, `from_amplimodes`, `from_cif`, `parse_basis`; `drs decompose`.
- `modes::Analysis::modes`: the modes of AMPLIMODES's block for FullProf.
- `isodistort`: `isodistort_cif`, `IsoDistort`, `IsoMode`, `split_label`; `drs isodistort`.
- Tests `mode_decomposition` (56 irreps of 10 AMPLIMODES analyses; 12 modes in their computed subspaces; a supercell case) and `isodistort` (a fixture).

### Known limitations

- AMPLIMODES adjusts the input before it decomposes it (Sr's x from −0.00003 to 0): from the input, 2 of 56 amplitudes differ in the 4th decimal, frozen; from its atom mappings, none.
- An irrep's label number (GM4+), isotropy subgroup and order-parameter direction are not computed: they need tables of irreps. The name is k, parity and dimension, or AMPLIMODES's when its modes are given.
- The ISODISTORT reader is checked on a fixture only: no file written by ISODISTORT is available yet.
- A polar subgroup's uniform translation is taken off by the mean displacement; unchecked against AMPLIMODES.

## drs-publish 0.12.0 — 2026-10-15 — Each preset in its journal's or university's typeface (L258); the set file watched in the app (L253)

Josu's L258: every character still typeset by LaTeX, in the face the preset's instructions name; free TeX faces.

### Added

- latex::FACES: TeX Gyre Termes with Termes Math (Times), TeX Gyre Heros (Helvetica) and Arimo (Arial) with New Computer Modern Sans Math and the face's letters and digits, EHU Sans from the project's Fonts/EHU Sans/, Latin Modern.
- presets.yml typeface: face and alternatives, cited; a set's typeface = "…" among them (module::set_typeface); the manifest records the face; check 8 says when a set is not in its preset's face.
- The app: Choose Typeface…, the inspector's Typeface; diskChanged and reloadTab for a set file changed on disk (§7).

### Changed

- The check of L103 follows the set's face (check::fonts(files, face); drs-publish check --face).
- The IUCr's left margin 11 mm, ACS's 11.5 mm (wider digits); a legend above the axes given the room of all its rows.

### Known limitations

- EHU Sans's files are needed from Josu; until they are in Fonts/EHU Sans/ a TFG is set in Latin Modern, said.

## app 0.33.0 — 2026-10-15 — Version 1's step C in the window: completion, single-peak fitting, the symmetry elements, the outputs' names; the app's documents from the Finder; Print and Page Setup



### Added

- P49: completion in the Editor for a control file, the core's names (control::complete): tables, keys with their meaning, values, the model's addresses.
- P51: Data ▸ Fit Peak… (⌥⌘K): one peak fitted by the core on the selection's range, its position, intensity, FWHM and η each with its esd, drawn over the points.
- P09: Refine ▸ Space Group… lists each operator's symmetry element (symbol, sense, location, screw or glide part), the core's (request 35).
- Editor ▸ Jump to Parameter also finds FullProf's symbolic names and a CIF's tags, by the core (request 36).
- The app's documents registered with the Finder: .drc, .drv and .drp owned, a results file (CIF) as an alternative, folders and files taken on the Dock icon; each opens in its project's window; plain placeholder icons (APP_CONVENTIONS §7).
- File ▸ Print… (⌘P) of the active plot or the Editor's file, and Page Setup… (⇧⌘P), through WebKit's print operation and the Mac's sheets (§3).

### Fixed

- The Settings window's self-check puts its change back when the main window says so, not after a time of its own (Publish's FAIL, L257).

### Known limitations

- Fit Peak holds one wavelength's peak: a Kα₂ peak needs λ₂/λ₁ from the core (request 38).
- Page Setup's choices last until the app quits; the icons are placeholders until Josu decides their design; the extensions .drv and .drp are Josu's to confirm.
- The Publish set's place in the file dictionary is asked of the core (request 37).

## 0.45.0 — 2026-10-15 — Version 1's step C: symmetry elements, single-peak fitting, completion for the editor, the names of the outputs

L260: P09's symmetry elements of every operator of the 530 settings, against cctbx (the app's request 35); P51's single-peak fitting with calibrated esd's, against FullProf on Y₂O₃; P49's completion source for the control file's editor, every name it offers read back by the reader; FullProf's symbolic names and the CIF's tags as the model's parameters (the app's request 36).

### Added

- `symmetry_elements`: `element`, `element_in`, `SymmetryElement`, `ElementKind` and `SpaceGroup::symmetry_elements` (kind, symbol, sense, axis or plane, intrinsic translation, location as ITA prints it); `drs sg` prints each operator's element.
- `peakfit`: `fit`, `PeakRequest`, `PeakShape`, `PeakFit`, `Estimate` (pseudo-Voigt, Gaussian or Lorentzian over a linear background, a Kα₂ doublet when asked, each parameter with its esd); `drs peak`.
- `control::complete`: `complete(text, line, column)`, `Completion`, `CompletionKind` (tables, keys with the dictionary's name, unit and meaning, choices, species, space groups, the model's addresses); `schema::TABLES` and `schema::choices`; `drs complete`.
- `output_names`: `fullprof_field`, `address_of_fullprof`, `cif_field`, `address_of_cif`; `drs name`.
- Tests `symmetry_elements` (7388 operators against cctbx, none differs), `peak_fit`, `completion`, `output_names`; `tests/data/make_symmetry_elements.py` and `symmetry_elements_530.txt`.

### Known limitations

- R 3 c and R −3 c in rhombohedral axes: the core names their glide n where the full symbol keeps the hexagonal c (2 of 529 letters, frozen).
- Single-peak fitting holds η in [0, 1]; FullProf's peaks reach η 1.0076. Against Y₂O₃'s whole-pattern profile, the bounds are frozen at 4.2 esd in position, 19.5 % in width and 8.8 % in intensity.
- 42 of FullProf's 2702 symbolic names in `manuals/` name what the core does not refine (λ₂'s own U, V, W; rigid bodies; a reflection's own broadening).

## app 0.32.0 — 2026-10-15 — Version 1's step B in the core: amorphous content, extracted intensities, the summary and listing; Publish's SVG and tables as Markdown, HTML and CSV (L256–L258)



### Added

- On core 0.44.0: amorphous content with an internal standard (P33), the file of extracted intensities in FullProf's format (P35), the summary and listing of a run (P57).
- Publish: SVG beside every PDF; tables as Markdown, HTML and CSV; every preset checked against its journal's own instructions.

## drs-publish 0.11.0 — 2026-10-15 — SVG; the tables as Markdown, HTML and CSV; the presets checked against their journals

What DESIGN.md still held that needs no new core work (L230): §2.4's SVG and data exports, §4's presets out of draft, each value compared with its target's instructions and cited (DESIGN.md §21).

### Added

- An SVG beside every figure's and table's PDF (pdftocairo -svg: glyphs as Latin Modern outlines; its size in points); the check of L103 reads SVGs too.
- Each table as Markdown, HTML (MathJax with mhchem) and CSV (export.rs).
- A figure's sub- and superscripts no smaller than the preset's minimum (\DeclareMathSizes in its own document); a text below the minimum said at the build.
- Presets: status checked and the date; no_colour_in_captions (check 7 follows the preset).

### Changed

- Every preset checked on 2026-10-15: IUCr (minimum 6.5 pt, column 88 mm, edge 0.35 pt, 600 dpi, no PDF), ACS (177.8 mm, letter, edge 0.5 pt, article), RSC (171 mm, 600 dpi, inset), Elsevier (190 mm, 10 pt body, legend 7 pt, EPS first), APS (2 mm characters: 8.5 pt minimum; letter, 179 mm), Powder Diffraction (ended with Cambridge in 2025; measured widths), the UPV/EHU's (the official EHU typeface noted for the TFG; the PhD note corrected); each journal's own deposition instead of the IUCr's inherited strings.
- F8: no text below the minimum (every second column number when two digits do not fit); a mode named by its symbol.

### Known limitations

- For Josu: the EHU typeface and figure text; Word manuscripts at Acta C, IUCrData and Powder Diffraction; Powder Diffraction's publisher; APS's 1 mm data points against a Rietveld plot's thousands.

## 0.44.0 — 2026-10-15 — Version 1's step B: amorphous content, the file of extracted intensities, the summary and listing

L255: FullProf's runs of the core's two files taken in, the refined λ reached and the user's f′, f″ rebuilt. L256: P33's amorphous content with an internal standard (Orlhac et al. 2001), P35's file of extracted intensities in FullProf's format with its cluster rule, P57's summary and listing, each against FullProf's or the paper's numbers.

### Added

- `qpa::amorphous_content`, `with_standard` and `brindley_tau`; the phase's `internal_standard` and the control file's `internal_standard_weight_fraction`; `drs refine` prints the amorphous content.
- `lebail::Extracted`, `clusters`, `extracted`, `int_file` and `CLUSTER_DEFAULTS`: FullProf's CODFILn_cltr.int; `drs results --exports` writes `Exports/<control> <phase> intensities.int`.
- `listing::summary`, `summary_text` and `listing_text`; `drs results` writes `Results/<control> <time> summary.txt` and `listing.txt`; the file dictionary's `export_intensities`, `run_summary` and `run_listing`.
- `lsq::bragg_points` and `lsq::statistics_over`.
- Tests `amorphous_content`, `extracted_intensities`, `fullprof_summary`, `p18_fullprof_cycles`.

### Changed

- `fullprof_profiles` compares nabd4_lambda_calc and y2o3_fprime_calc (68); `prf_files` reads 200; `berar_factor` counts 131; `p18_wavelength_refinement` checks the refined λ against FullProf's.

### Known limitations

- FullProf's refinement of λ (nabd4_lambda) ends in a worse minimum with U < 0; its χ² from cycle 5 is 1.2–1.6 % above the core's at the same values, frozen (L255).
- The summary gives no density: the core has no atomic masses yet.

## app 0.31.0 — 2026-10-15 — A file changed on disk noticed in every module; a context menu on every surface



### Added

- A file open in a tab and changed outside the app is read again when it has no unsaved changes, or asked about (read the file on disk, or keep my changes); the module's badge says so meanwhile (APP_CONVENTIONS §7, §2).
- A context menu on the plot, a tab, Refine's outline and inspector rows, a table and the activity bar, from the registry's commands and the module's own actions at that place (§6); a module can be hidden from the activity bar and shown again.

### Known limitations

- Explore, Plots and Publish have no hook of their own yet: their tabs are read again by closing and opening them, and their own saves may be taken for a change.

## 0.43.0 — 2026-10-15 — FullProf's references for P29 and P11 taken in: b's imaginary part and the Friedel average

L251: FullProf's runs of ce, ce2, nabd4 and highpress (L239). Surface micro-absorption and the user's neutron scattering lengths are reproduced from FullProf's own calculations, with two rules found there.

### Added

- `Scatterer::Neutron::b_imaginary`, `UserScatterer::scattering_length_imaginary` and the control file's `scattering_length_imaginary`: an absorbing nucleus's b″, which FullProf takes from DFPP (its manual says it ignores it).
- `Structure::intensity` and `intensity_derivatives`: the powder |F|², the average over h and −h, which differ when a scatterer has an imaginary part in a group without a centre.
- Two files for FullProf to run, in `manuals/fullprof/FP_tests/`: `nabd4_lambda.pcr` (λ refined with the cell fixed) and `y2o3_fprime.pcr` (the user's f′ and f″), with `y2o3.dat`.

### Changed

- The calculated pattern, its derivatives and the reflection rows take the averaged |F|²: unchanged where no scatterer has an imaginary part.
- `fullprof_profiles` compares ce_calc, ce2_calc and nabd4_calc (66 in all) and leaves out highpress_calc with its reason; `prf_files` reads 196; `berar_factor` counts 127.

### Fixed

- A user's b with DFPP ≠ 0 (nabd4): the pattern was 7·10⁻³ off FullProf's, now 1.3·10⁻⁴.

## drs-publish 0.10.0 — 2026-10-15 — The checks against the three papers; F2, F8, T7 and T8

McCusker et al. 1999, Young et al. 1982 and Toby 2006 (L238) read into the checks and T1's words; the catalogue's last figures and tables; case E at 453 K (L243).

### Added

- Check 9, what a Rietveld paper reports (McCusker §§10–11, Young's guidelines): the plot, the high angles magnified for X-rays, T1, the background-free R's, the instrument and temperature, the deposited data, the largest structural correlation, the method's name.
- Check 6: χ² of 2 or more and R_F above 5 % to be explained (Toby 2006, p. 69).
- F2, models of one data file compared ([[source.models]]); F8, the correlation matrix and the parameters' paths (Trajectory::from_end, in_esd).
- T7, a series with a column per point (T2 or T5); T8, AMPLIMODES's input (the parent, the distorted structure, (P, p)).
- T1: each agreement factor defined in its note; (Δ/σ)_max of the last cycle.

### Changed

- T1 names R_I and R_F, not R_Bragg, which may mean either (Toby 2006, p. 68); s.u. for esd in the checks and T1.
- The example results written again by core 0.41.0 (correlations, cycles); case E of the La2CoMnO6 series at 453 K, 179.99 °C in its data file's header (L243), not the course's 500 K.

## app 0.30.0 — 2026-10-15 — The Settings window with panes; a project's own settings; with atoms and phases added, the jump to a parameter, Explore's microstructure, the series sent to Publish, on core 0.42.0 (L240–L247)



### Added

- The Settings window (⌘,), one for the app: General, Appearance, Editor, Plots, Data formats, Refinement, External programs, Advanced; each setting explained in a line, every change applied at once in every window (APP_CONVENTIONS §4).
- New settings: the Editor's font size, tab width, line numbers, wrapping and saving by itself; Refine's plot axis and the intensity scale; a version before each run; the notification when a run ends; the recent projects kept; where FullProf and GSAS-II are; the Log's length; Reset All.
- File ▸ Project Settings…: a project's own cycles, Pause every, plot axis and scale, Yobs as a line and tab width, said in the status bar.
- Refine: Add Atom…, Add Phase… from a symbol and a cell, Add Phase from CIF… (P48); Editor ▸ Jump to Parameter (⌃⌘J) from a line of the core's outputs (P49).
- Explore: the Microstructure tab (size and strain per cycle, each reflection's breadths); the wavelength and micro-absorption in the traces and the pair view.
- Plots: Send Series to Publish… with each temperature's source.
- The Bérar factor printed with the refinement and kept in the results file (core 0.42.0, P43).

### Changed

- Case E's temperature is 453 K (179.99 °C), Josu's word (L243).

### Fixed

- A window holding an old copy of the settings no longer undoes a change made in another, and a setting put back to its default reaches every window; Refine's Pause every field keeps its value.

### Known limitations

- Some §4 choices are said as information, not offered: the update service, the minimap, an app-wide convergence, the publication style (Publish's), the readers' order (the core's).

## 0.42.0 — 2026-10-15 — The Bérar factor (P43)

L243: Bérar & Lelann's factor for serially correlated residuals, as FullProf reads the paper, checked against FullProf's printed factor in the 119 outputs of the repository that carry it.

### Added

- `lsq::berar_factor` and `Refinement::berar_factor`: √(S″/S) of Bérar & Lelann (1991), eqs. 10–12 with the appendix's 2π, on the points with Bragg contributions.
- The results file's `_drs_berar_factor` (read back as `Results::berar_factor`); `drs refine` prints it.
- Test `berar_factor`: 58 outputs compared, 34 within 0.05 %, 17 within the rounding of the .prf's Ycal, 7 frozen; the other 61 left out by reason, counted.

### Known limitations

- The default follows FullProf's reading of the paper, whose equations and appendix differ (L243); the equations alone give 8–12 % less.
- LBF_lmopr (+0.09 %) and sr-ox-rest_conv (+0.88 %) differ from FullProf's factor without a found reason.

## app 0.28.0 — 2026-10-15 — Atoms and phases added in the app; from an output's line to its parameter



### Added

- Refine ▸ Add Atom…, Add Phase… (a symbol and a cell) and Add Phase from CIF…, through the core's edit (P48).
- Editor ▸ Jump to Parameter (⌃⌘J): the parameter a line of the core's outputs names, selected in Refine's outline and inspector (P49).

### Fixed

- The self-check's P46 line in a fresh home: the reopening is said not checked, and the recording is checked (L239).

### Known limitations

- FullProf's .out names and the CIF's tags wait for the core (request 36); removing an atom or a phase is not offered; completion in the Editor and single-peak fitting wait for the core.

## app 0.27.0 — 2026-10-15 — Version 1's gaps closing: windows and tabs reopened, Refine's plot in d or Q, the space-group view; Explore's Cell and Bonds per cycle; Publish's VESTA, Bilbao and checks; on core 0.41.0 (L232–L237)



### Added

- Windows and their tabs reopened where they were at launch (P46); Refine's plot in 2θ, d or Q, on a linear, √ or log scale (P51); Refine ▸ Space Group… with the operators and the Wyckoff positions (P09).
- Explore: the Cell (a, b, c, angles, V with esd's) and the Bonds (distances, angles, bond-valence sums with esd's) per cycle.
- Publish: VESTA files, the Bilbao server's format and the data (D5–D7); the checks before export listed in the Problems panel.
- Plots: Figure Size… at the journals' widths; Open Results Series… (refinements against T).

### Changed

- On core 0.41.0: the wavelength refined (P18), scattering factors given by the user (P11), Bragg–Brentano surface micro-absorption (P29).

## drs-publish 0.9.0 — 2026-10-14 — VESTA, the Bilbao server's format, the data; the checks before export

What needed no new core work (L230): D5–D7 and the checks 1–7 of DESIGN.md §5.4, shown in the app's Problems panel.

### Added

- D5: the mode phase in VESTA, a .vesta per irrep with its displacements as vectors (Å along the axes; FullProf 8.20 writes ten times as much).
- D6: the phase in the Bilbao server's format, its non-standard setting said; D7: the profile as .xye and CSV, the CIF's loops as CSV.
- checks::run and drs-publish check <set.drp>: one version, the tables equal to the CIF, the CIF as checkCIF sees it, the preset, ThesisFigures, the refinement as reported, the captions.

### Fixed

- Every figure was 1.2 mm wider than its width class: the PGF's last line ended in a space after the picture (found by check 4).

## drs-publish 0.8.0 — 2026-10-14 — The set edited in the app; the CIF, pdCIF and supplementary material; E14b

The module's second step (L229).

### Added

- [[file]] items: D1 the CIF for deposition, D2 the pdCIF, D4 the supplementary material.
- The set edited through its text: items added, removed, moved; their options (module::add_item, remove_item, move_item, options, set_option).
- T5, T6 and F6 report the amplitude set nearest the sign reference and name the operation (E14b, core 0.38.0).

### Known limitations

- The core reads one AMPLIMODES run as three analyses (P24); PL28 b waits on the core's map scale (P25).

## drs-publish 0.7.0 — 2026-10-14 — The Publish module in the app, first step

L207: the app's Publish module calls drs-publish as a library.

### Added

- module.rs: a set's items and state, the preset written into a set, a set made from a plot set (PL28 a; drs-publish new --plots), an item's page as pixels with its size in mm, the TeX found.
- build::build_with: progress per step and Stop.

### Fixed

- poppler's tools found where Homebrew and MacTeX put them, for an app opened from the Finder.

## drs-publish 0.6.0 — 2026-10-13 — The mirror solutions said in T5, T6 and F6

After Explore's Modes (L203): which amplitude set the tables and the figure show, and how many describe the structure.

### Added

- The notes and the caption name the set shown and the number of equivalent sets (Modes::equivalent_amplitudes).

### Known limitations

- 20 of case A's 24 sets changed an irrep's amplitude (P21, met by core 0.39.0).

## drs-publish 0.5.0 — 2026-10-13 — F3 and the series view from a plot set

PL15 on core 0.33.0: figures drawn from a plot set file (Plots/{name}.drv).

### Added

- F3: a stack, a projection with hidden parts, a map or a difference map, of the whole pattern or of chosen windows (plots::displace, project, map, windows).
- The series view: each window's integral, centroid, peak or background against the variable, with σ.
- Example: the Tb₂Pd₂Sn buffer saved as a .drv (tb2pd2sn-tfg).

## drs-publish 0.4.0 — 2026-10-13 — A phase refined by coordinates in F6; B and U_eq from the results file

On core 0.28.0 (P19, P20): F6 through geometry::decompose and Scene::reference_fract; T2's B and U_eq/B_eq as the results file writes them.

### Added

- F6 for a phase refined by coordinates, described by another file's mode phase (modes, shift).
- Example: case A refined by coordinates (sr2scsbo6-coords-tfg).

## drs-publish 0.3.0 — 2026-10-12 — Every calculation is the core's

After L139 (L156's merge): F1's difference, T1's range and step, T2's B and F6's structure taken from the core.

### Added

- tests/nophysics.rs: the guard against π, roots, powers, trigonometry and symmetry operations outside the PGF writer.

### Changed

- F1: diagnostics::difference; T1: data::extent and Pattern::mean_step; T2: the refined B, no 8π²U; F6: geometry::scene and modes::morph, only projected by the crate.

## drs-publish 0.2.0 — 2026-10-11 — Bonds, the modes' tables and figures, the cell against the variable

The catalogue's next items on core 0.25.0 (L156): T4–T6 and F4–F6.

### Added

- T4: the CIF's distances and angles with their su's, the linkage angles between named elements, ⟨d⟩ and bond-valence sums from the core with the set's charges.
- T5: the irreps' amplitudes with AMPLIMODES's k, direction, isotropy subgroup and number of modes; T6: each mode's signed amplitude.
- F4: the cell reduced to the parent's against temperature or composition, the phases as bands; F5: the amplitudes as bars or against a series; F6: the displacements on the reference structure, projected, with polyhedra.

### Known limitations

- No su's on ⟨d⟩ and the bond-valence sums (P18); no parent cell for a phase refined by coordinates (P17).

## 0.41.0 — 2026-10-14 — The wavelength refined, the user's scattering factors, surface micro-absorption (P18, P11, P29)

L233: version 1's step A. The wavelength is a refinable value; a .pcr's line 12 and the control file give the user's scattering factors; Bragg–Brentano surface micro-absorption is calculated and refined, FullProf's CeO₂ refinements kept at their Cp and τ. The Bérar factor (P43) waits for its paper.

### Added

- `Field::Wavelength` (`pattern1.wavelength_1`): λ₁ refined, λ₂ with it for one wavelength; a doublet's refused; line 14's Lambda codeword imported and exported.
- `Field::Micro` and `Model::microabsorption` (P0, Cp, τ; FullProf's eqs. 3.35–3.36): `Model::micro_absorption`, refinable, line 15 read with its codewords and exported, control-file keys `microabsorption_p0`, `_cp`, `_tau`; three dictionary entries.
- Line 12 read (`pcr::PcrScatterer`, `powder::user_scatterers`); `scattering::UserScatterer`, `scattering::resolve`, `Scatterer::XrayTable`; the control file's `[[pattern.scatterer]]`; Nsc and line 12 in the export; `Model::user_scatterers`.
- Test `pattern_terms` (P18, P29, P11).

### Changed

- A `.pcr` with Nsc ≠ 0 or with surface micro-absorption is read (they were refused).
- Line 15 is read only for X-ray 2θ patterns in Bragg–Brentano geometry, as the manual says.
- `explore_stops`: each test reads its own copy of the cases.

### Known limitations

- The refined λ and the user's scattering factors are not compared with FullProf's numbers: no FullProf output in the repository has them (runs asked of the merging session, L233).
- FullProf's ce and ce2, started from FullProf's values, move 3 and 8 profile parameters by up to 2.8 of the engine's esd (χ² −0.13 %, −0.56 %): frozen.
- The Bérar factor (P43) is not built: Bérar & Lelann (1991) is not in the repository.
- The app does not compile on 0.41.0 at `app/src-tauri/src/modules/refine.rs:741` (the match over `Field` needs `Wavelength` and `Micro`).

## app 0.26.0 — 2026-10-14 — Windows and tabs reopened at launch; Refine's plot in d or Q, on √ or log; a phase's space group



### Added

- State restoration (P46): the windows open at the quit come back where they were, each project with its tabs; Settings ▸ Reopen windows and tabs at launch.
- Refine's plot in 2θ, d or Q, the conversion the core's (to_axis), and on a linear, √ or log intensity scale (P51).
- Refine ▸ Space Group… (P09): symbols, lattice and centring, the operators with their kind, the Wyckoff positions with the atoms on each, all from the core.

### Fixed

- The Problems panel opens for a problem not shown before, even when the count of problems stays the same, and when the file with it comes to the front (L208).

### Known limitations

- The tabs come back in the main area, not in the groups they were split into; the symmetry elements wait for the core (request 35); single-peak fitting in Data waits for the core (L230, part B).

## app 0.25.0 — 2026-10-14 — A past cycle refined again; the publication set edited, with its CIF, pdCIF and supplementary material; Send to Publish built in LaTeX; on core 0.40.0 (L225–L230)



### Added

- Explore: a past cycle opened and refined again from there (W11).
- Publish: the publication set edited in the app (items added, removed, moved, their options); the CIF for deposition, the pdCIF and the supplementary material built from it; the mode tables and figures show the set nearest the AMPLIMODES reference and name the operation.
- Plots: Send to Publish builds the set in LaTeX, compared with Plots' own export.

### Changed

- On core 0.40.0: a run's pauses, its stop's facts and the steering's notes kept as data (E15); the plot set file keeps each view's steps and scale (PL27).
- With app 0.24.0: auto-save with the versions kept, numbers with a point, the system's accent colour, Rename… of the control file only, Delete to the Trash, Quit asks every window, ⌘C on a plot copies it as an image; Export ▸ Plot's PDF at its size (PL31).

## app 0.24.0 — 2026-10-14 — Files save by themselves, the versions kept; the system's accent; Rename; the Trash; Quit asks



### Added

- Auto-save (L224): a file saves two seconds after its last change and when its window closes, goes to the background or the app quits; File ▸ Revert To lists Last Opened and the file's versions (snapshots and accepted refinements), what a revert replaces kept as a version first.
- File ▸ Rename… renames a control file, the control file only (L224).
- A hidden panel's stub edge brings it back; arrow keys in the sidebar's list and Refine's outline; ⌘C on a plot copies it as an image in every module.

### Changed

- Buttons, selections and highlights take the system's accent colour, the app's teal when multicolour (L224).
- A file the app removes goes to the Trash (macOS 15 and later).
- Quitting asks every window, as closing a window does, about what could not be saved.

### Fixed

- Export ▸ Plot as PDF was 4/3 too large: the plot's CSS pixels are now taken at 96 to the inch (PL31).
- The app compiles on core 0.40.0 (E15's run record in Explore's Rust side).

### Known limitations

- The files named by a renamed control file's old stem (versions, the refinement log, notes) keep that name.

## 0.40.0 — 2026-10-14 — A run's pauses, stop and steering kept as data (E15); the set file's views (PL27)

L225: Explore's E15 and Plots' PL27. A run's record keeps where it paused and what changed, why it stopped with the facts, and what the steering did, as data written to the results file and read back. The plot set file keeps each view's steps and scale, the map's scale and the normalisation.

### Added

- `lsq::Pause` (after, why, λ, changes, λ restarted), `Session::pause()`, `Trajectory::pauses`; the session text and the results file (`_drs_pause_*`) keep them.
- `lsq::Stopped` (the stop and its cycle) in `Trajectory::stop`; `Stop::facts` and `Stop::from_facts`; the results file's `_drs_stop_after` and `_drs_stop_{fact}`.
- `lsq::SteeringNote` and `NoteKind` (maximum-shift, at-limit, reflected, wrapped, refused) in `Cycle::steering`; the results file's `_drs_steering_*` loop.
- `drs refine --trace` writes the steering's notes as data, the stop's facts with its cycle, and the pauses.
- `plots::SetFile`: the stack's step with its scale, `projection` with its steps and scale, `map_scale`, `normalisation`; `drs plots normalise` takes the file's normalisation without `--by`, and `drs plots save --by` writes it.
- Tests `explore_stops` (4) and `pl27_each_views_steps_and_scale`.

### Changed

- `SetFile::open` accepts an x quantity other than 2θ: the view's, the members read in 2θ.
- `lsq::Cycle` has the field `steering`; `Trajectory` has `pauses` and `stop` (`Trajectory::new` builds one without them).

### Known limitations

- The app does not compile on 0.40.0 at two struct literals in `app/src-tauri/src/modules/explore.rs` (lines 28 and 69), the Explore stream's to update (L225).
- A session resumed from its file keeps the earlier cycles' numbers, not their notes.
- PL28 is Publish's (drs-publish); the core gives PL27's data.

## app 0.23.0 — 2026-10-14 — The Publish module in the app; Plots' Figure Look and Send to Publish; Explore's Modes against the AMPLIMODES reference; on core 0.39.0 (L211–L222)



### Added

- Publish, the sixth module (L207): the publication set's items, an item at its printed size, a build with progress and Stop, the preset chosen; a set made from a plot set.
- Plots: Send to Publish… hands a plot set over; Figure Look for all five views (on white, open axes, scale bars, colour bars' scales), figures at 8.5, 12 and 17 cm.
- Explore's Modes tab: esd bars at every cycle, each cycle named against the AMPLIMODES or FullProf reference (the mirror marked), FullProf's reported amplitudes beside the core's.
- Explore: the run's stop and steering notes marked on the traces (W11).

### Changed

- On core 0.39.0: anisotropic strain by Stephens' quartic form (P25 begun); the equivalent amplitude sets only of the same structure (P21: case A has 4, not 24); files written as the view shows the set.
- With app 0.22.1 and 0.21.x: the quartic terms named S<hkl>; the Problems panel shown on a new problem; nine small APP_CONVENTIONS items met; Settings ▸ Shortcuts…; a .pcr checked as it is typed.

## app 0.22.1 — 2026-10-13 — The Problems panel shown; Close Window ⇧⌘W; APP_CONVENTIONS.md gone through



### Added

- Window ▸ Bring All to Front; Edit ▸ Use Selection for Find (⌘E) in the Editor; the settings in ⌘K.
- A window closed with unsaved files names them and asks first.

### Changed

- The Problems panel is shown when the Editor's check finds a new problem (the default layout hides the bottom panel, L208).
- Edit ▸ Find (⌘F) also finds in the table shown; a click on the status bar shows the Log; the toolbar's visibility is kept with the layout.
- A refinement runs at macOS's utility priority, as Explore's grids (L150).

### Fixed

- Close Window had ⌘W, Close Tab's shortcut; it is ⇧⌘W.
- The app compiles on core 0.38.0's quartic-strain parameters (S400 …, P25).

### Known limitations

- APP_CONVENTIONS.md items not yet met are listed in app/README.md (0.22.1), small and large, with what each needs.
- Not yet confirmed in the built app: its self-check stalled while five streams compiled at once, and is to run again.

## 0.39.0 — 2026-10-13 — Equivalent amplitude sets only of the same structure (P21)

L216: Publish's P21. `Modes::equivalent_amplitudes` offered, for each parent operation outside the normaliser of the phase's group, another structure as the same. Case A now has 4 sets, not 24.

### Added

- Test `p21_equivalent_sets_keep_each_irreps_amplitude`.

### Fixed

- `Modes::equivalent_amplitudes` (and `same_structure`, `nearest_to`, `reference_signs`, `drs modes --equivalent`) checks the whole moved cell against the whole cell of the amplitudes found: no set moves amplitude between irreps.

### Known limitations

- Publish's test `the_mirror_solutions_are_counted_and_frozen` expects 24 sets; with this version there are 4 (the Publish stream's to update).

## 0.38.0 — 2026-10-13 — Files written as the view shows the set; the core's CSV read back in d and Q; the modes' esd's at every cycle, the sign reference, FullProf's reported amplitudes

L212: Plots' PL25 and PL26, Explore's E14a–c.

### Added

- `PlotSet::as_shown`, `Normalisation::words`, `csv_with`, `csv_long_with`, `derived_text_with` (PL25); `drs plots export` and `derive` with `--by …`.
- `state::summary_at` with a covariance given directly, `geometry::Propagator::from_covariance` (E14a).
- `modes::SignReference`, `Modes::sign_reference`, `reference_amplitudes`, `reference_signs`, `Equivalent::words`, `is_identity`; the control file's `mode_reference = { file, analysis }` (E14b).
- `modes::FullProfAmplitudes`, `read_fullprof_amplitudes`; `drs modes` prints FullProf's reported amplitudes from the `.sum` beside a `.pcr`, and a state against its sign reference (E14c).
- Tests `pl25_the_view_written_as_shown`, `pl26_the_cores_csv_read_back_in_d_and_q`, `e14a_esds_at_every_cycle`, `e14b_the_sign_reference`, `e14c_fullprof_reported_amplitudes`.

### Changed

- `modes::Modes` has a `sign_reference` field (a struct literal needs it).

### Fixed

- The table reader did not know the core's own x columns `d_angstrom` and `q_inverse_angstrom`: a CSV in d was read in 2θ, a long one as two columns (PL26).

### Known limitations

- Core 0.37.0's `Field::Quartic` and `Field::QuarticShare` are not covered by the app's `match` at `app/src-tauri/src/modules/refine.rs:742`: the app does not compile on 0.37.0 or later until the app stream adds the one-line arm (L212).

## 0.37.0 — 2026-10-13 — Anisotropic strain by Stephens' quartic form, against FullProf's run of Nd₂NiO₄

L211: P25 begun, the gate rule kept. The quartic form of strain for the Laue classes −1, 2/m and mmm, checked reflection by reflection on the one FullProf run in the repository that uses it; the rest of P25 has no FullProf reference.

### Added

- `microstructure::QuarticStrain` (S_HKL and the Lorentzian share ξ), `quartic_terms` (the terms a group allows), `Phase::quartic_strain`; `Model::peak_width`, `peak_width_components`, `peak_window_width`, `peak_gaussian_h2`, `strain_widths`: widths that depend on the reflection.
- `Field::Quartic`, `Field::QuarticShare` with analytic derivatives; `phase1.quartic_strain.400`, `phase1.quartic_strain_lorentzian`; ξ held within 0–1.
- The `.pcr`'s LINE 37–38 read (Str = 1, Strain-Model 2, −2, 3) and written; the control file's `quartic_strain` and `quartic_strain_lorentzian`; the dictionary's `quartic_strain` and `quartic_strain_lorentzian`.
- Test `fullprof_anisotropic_strain` (3): NdNi-LTOa's 366 reflections of the `.out` and 366 rows of the `.mic`; NdNi-LTOa in `engine_derivatives`.

### Changed

- `drs mic` leaves out magnetic phases and says so; the peaks of `Model::peaks` take their reflection's anisotropic strain.
- A `.pcr` with anisotropic size or strain is refused with the reason for its model (no FullProf run with output, words without an equation, Stephens' relations not in the repository, selective broadening is P26), where it said "version 2 (P25)".

### Fixed

- Pd3MnD's size model −14 was called spherical harmonics in the profile index: it is selective broadening (P26).

### Known limitations

- Not calculated, for want of a FullProf reference with its output in the repository: the quartic form of the tetragonal, trigonal, hexagonal and cubic classes (Stephens' relations), −1 from a `.pcr`, size by the quadratic form, spherical harmonics, platelets and needles, FullProf's strain models with Str = 0.
- The peak window with the quartic form's Lorentzian part is not checked on a FullProf profile: the run's second phase is magnetic.

## app 0.21.1 — 2026-10-13 — Run Next Plan Step ⌃⌥⌘R; ⇧⌘D on File ▸ Open Data File… only



### Fixed

- Two shortcuts were on two commands in the menu bar at once (found by app 0.21.0): Run Next Plan Step is now ⌃⌥⌘R, One Cycle keeps ⌥⇧⌘R; ⇧⌘D is File ▸ Open Data File…'s only, and Data's menu item runs the same command without a shortcut of its own (L209).

## app 0.21.0 — 2026-10-13 — A .pcr checked in the Editor; shortcuts that can be changed



### Added

- Settings ▸ Shortcuts…: a command's shortcut changed, removed or given back; each checked against the system's own and the app's other commands; kept in the settings and taken by the menus at once (APP_CONVENTIONS §6).

### Changed

- A FullProf .pcr in the Editor is checked by the core as you type (check::check_pcr, as drs check), its problems underlined and listed as a control file's (request 13).
- Help ▸ Keyboard Shortcuts points to Settings ▸ Shortcuts; the Editor's manual says how a data file's lines are linked to its points (the core's point_lines).

### Known limitations

- Two shortcuts are on two commands in the menu bar at the same time: ⇧⌘D (File and Data ▸ Open Data File…) and ⌥⇧⌘R (Refine ▸ One Cycle and Run Next Plan Step); either can be changed in Settings ▸ Shortcuts.

## app 0.20.0 — 2026-10-13 — Explore's Modes tab; Plots normalises and converts as views, and writes and reads sets as CSV; the core's requests used across the app (L201–L206)



### Added

- Explore: the Modes tab, each irrep's and each mode's amplitude per cycle with the end's esd's, and the amplitude sets equivalent to the end's (the mirror solutions) when the parent is known.
- Plots: Normalise… (monitor, time, maximum, the same area over a range, a peak's integral, a reference member) and Convert X… (2θ, d, Q, sin θ/λ, 1/d) as views, said in each; Export ▸ Plot Set as CSV…; Open CSV as a Set….
- Data: Export the Processed Pattern…; a point selected in the plot selects its line in the Editor, and the other way.
- Refine: Export for FullProf into Exports/FullProf/.

### Changed

- Find Peaks and Automatic Background take their words from the core; the sidebar names files by the core's file dictionary.
- A temperature read from a temperatures file is kept by the plot set file; a .drv is known by its content (core 0.36.0).

### Known limitations

- A CSV in d or Q reads back as 2θ (PL26); a CSV does not record the view's normalisation (PL25).

## 0.36.0 — 2026-10-13 — A temperatures file's value kept by the plot set file; a plot set known by its content

L202: Plots' PL23 and PL24.

### Added

- A set file member's `variable = { value = …, temperatures = "…" }`: a value read from a FullProf temperatures file, kept with that source through Save and Open (PL23).
- `formats::project::PLOT_SET`: a plot set known by its `format` line wherever it lies, and by the extension `.drv` (PL24).
- Tests `pl23_a_temperatures_files_value_kept`, `pl24_a_set_file_by_its_content` (`plots_sets`).

### Changed

- `plots::SetMember::value` is `Option<(f64, Found)>` (typed or a temperatures file), where it was `Option<(f64, String)>`.

### Fixed

- `drs plots save` from a `.drv` lost the values chosen with `--temperatures` and `--set`.

## app 0.19.1 — 2026-10-13 — The processed pattern exported, the line of a point both ways, the analyses in the core's words, exports for FullProf



### Added

- Data ▸ Export the Processed Pattern… writes the pattern after the operations by the core (formats::write, request 12), proposed as Data/<data> processed.xys.
- Refine ▸ Export for FullProf writes the .pcr and the files it names by the core (control::pcr_export) into Exports/FullProf/ (request 32).

### Changed

- A point selected in the plot selects its line in the Editor and the other way, by the core's table of each point's line (formats::point_lines, request 15), FullProf's free format included.
- Find Peaks and Automatic Background are worded from the core's dictionary (data_ops::settings, request 16).
- The sidebar names files through the core's describe_in and describe_path (request 31); the app's own rule for a results file is gone.

### Fixed

- After Save Project As the data open and Start a Refinement… is enabled (confirmed in the built app; the fix of 0.18.1).

### Known limitations

- A .pcr in the Editor is coloured, not yet checked by the core's check_pcr (request 13).

## 0.35.0 — 2026-10-13 — A CSV of several columns read as a set; a set in d, Q, sin θ/λ and 1/d as a view; a results file named a result by its place and its first line

L199: Plots' PL8 and PL10; the app's requests 31 and 32 (met by 0.27.0) and the command line's last gap.

### Added

- The table reader (PL8): several y columns named give one scan each, σ paired by name; a long table (a member column) one scan per member; empty cells; PL7's # lines read back as each member's file, format and variable with its source; the columns not read said (`TableLayout::series`, `member`, `unread`).
- `plots::Found::CsvLine`: a member's variable from its CSV line.
- `plots::x_view`, `XQuantity`, `LambdaFrom`, `XMember`, `XView` (PL10); `drs plots view`.
- `formats::project::describe_path`; `drs data info` names a file by its place without `--in`.
- Test `plots_tables` (5).

### Changed

- A table read with several y columns named gives several scans, where it gave the first.

### Fixed

- A results file was named a control file by its content (it holds its control file); its first line now names it (`project::RESULTS`).
- The 41-column series CSV dropped 19 of its 20 members silently (STUDY.md §10.3).

## app 0.19.0 — 2026-10-13 — Plot sets saved and opened; members after Data's operations; Explore draws the core's ellipses and map; on core 0.35.0 (L196–L200)



### Added

- Plots: Save Plot Set (⌘S), Save Plot Set As…, Open Plot Set… and a double-click on a .drv, the set kept as Plots/{name}.drv with its order, rule, reference, windows and each typed value with its source.
- Plots: each member drawn after the operations Data keeps for its file, and said so; Reload Members.
- Plots: the stack opens with the core's proposed step, its rule said; Read Variable From… a FullProf temperatures file; a series point below 2σ drawn open.

### Changed

- Explore draws the χ² map and every ellipse from the core's points (Grid::map, Ellipse::outline); the 3D view's listeners go when the view goes.
- With app 0.18.1's fixes: a plan step that frees nothing no longer stops every run; a Refine plot beside another tab acts on its own document; Save Project As opens its data.

### Known limitations

- A temperature read from a temperatures file is not kept by the plot set file (PL23); a .drv is known by its extension, not yet by its content (PL24).

## app 0.18.1 — 2026-10-13 — Runs start after a plan is written; excluded regions drawn in a split view kept; Start a Refinement… after Save Project As



### Changed

- The tour's player drags and clicks on the active tab's plot.

### Fixed

- A plan step naming parameters not yet refined made the whole control file an error, and every run was refused: Add Plan Step now lists the parameters it names (addresses, scopes, prefixes) among the refined ones, a link takes its value out of them, and a change that would leave a step freeing nothing is refused.
- An excluded region drawn, an edge moved or a background point picked on a Refine plot in a split group acted on the active document, silently: each now acts on the plot's own document and a refusal is said.
- After Save Project As the new project's data did not open, and Start a Refinement… stayed disabled.

### Known limitations

- One self-check of the built app: the plan and the excluded region pass; the Save Project As scene was not reached (the check's script waited for a dialog the data window rightly does not open).

## 0.34.0 — 2026-10-13 — Normalisation as a view and the constant counting time of multi-scan files; a set as CSV in 2θ, d or Q; FullProf's temperatures file; a window's value below 2σ flagged

L195 (L170, step 5): Plots' PL6, PL7, PL21, PL22.

### Added

- `plots::Normalisation`, `normalise`, `Normalised`, `PlotSet::normalised` (PL6): monitor, time, maximum (in a range), the same area over a range, a peak's integral, a reference member's area; a view, no file changed.
- `data_ops::constant_time` and `RATE` (PL6): the scans of a multi-scan UXD or XRDML brought to one counting time per point, overlaps keeping the first; `drs data constant-time`.
- `plots::csv` in any x quantity and `csv_long`, `plots::converted` (PL7); `drs plots export --csv --long --x d|q --wavelength`.
- `formats::write`: xye in d and Q, and qye; the core's xye reader takes the axis back from the header (PL7).
- `plots::read_temperatures`, `TemperatureFile`, `PlotSet::variable_from_file`, `Found::TemperatureFile` (PL21); `--temperatures` for every `drs plots` command.
- `plots::WindowValue::significant`, `significance`, `WINDOW_SIGNIFICANCE` = 2 (PL22).
- `drs plots normalise`.
- Test `plots_views` (5).

### Changed

- The UXD reader marks `_CPS` blocks and the XRDML reader intensities in cps as rates (`Info::other`, ('intensity', 'counts per second')).

### Known limitations

- The constant counting time is given by the caller: WinPLOTR's own choice is in a note (S12) not in the repository.
- A temperatures file naming a file twice with different values gives it none (FullProf's own example does so; which FullProf takes is not stated).
- A normalisation's factor is taken as exact.

## app 0.18.0 — 2026-10-13 — Two refinements compared in Explore; Plots' linked overlay and the projection's apart member; on core 0.33.0 (L188, L189, L194)



### Added

- Explore's compare view (E12 over the core's E11): the run against a version chosen in the sidebar, the two paths, their difference in σ of the reference, each end's agreement, and the pair's two ellipses, every number the core's.
- Plots' linked overlay (F6): the whole range on the left, a window on the right, one cursor on both, the core's difference of two members under them; Reverse Order, Sort by Variable, Read Variable From… the file's name.
- The tour's handles on Explore's selects and sliders.

### Changed

- Plots' projection keeps a member without a temperature apart, its peaks where they were measured.
- Built on core 0.33.0: the plot set file, the map's rows apart, a proposed stack step and members after Data's operations are in the core, for the module's next step.

### Known limitations

- A series point below 2σ is not yet drawn open: the core's flag (PL22) comes in its next step.

## 0.33.0 — 2026-10-13 — The plot set file; the map's rows apart, a proposed stack step, members after Data's operations; the χ² grid as a map and the ellipse's outline; the E11 tests one computation on every platform

L192: L187's two E11 tests; Plots' PL12, PL18–PL20; Explore's E13a, E13b.

### Added

- `plots::SetFile` (PL12): `Plots/{set}.drv` read, written and opened (`read`, `to_text`, `open`, `of`, `reference_in`), `set_base`, `relative_to`; the file dictionary's `plot_set`.
- `drs plots save`, `drs plots export --csv` (`plots::csv`, the wide CSV in 2θ), `drs plots step`; every `drs plots` command opens a `.drv`.
- `plots::Map::apart_z` (PL18): the rows of the members apart, less the reference in a difference map.
- `plots::proposed_step` and `ProposedStep` (PL19), with its rule.
- `plots::open_with`, `open_buf_with` and `Operations` (PL20): members after Data's operations.
- `lsq::Grid::map` (E13a) and `lsq::Ellipse::outline` (E13b).
- Test `plots_sets` (7).

### Changed

- A map's colour range covers the rows apart too (Tb₂Pd₂Sn: 905 to 68658, the difference map ±66204).
- `lsq::compare`'s combined σ is √(σA² + σB²) by products, a sum and a root, the same bits in either order on every platform (it was `hypot`, which macOS's C library does not keep symmetric).

### Fixed

- `paths_compared`'s `swapped` and `e11bc_the_results_file` failed on macOS by the last bit (L187): the first by `hypot`'s order, the second by comparing a start rebuilt from the end with the session's own; both now one computation.

### Known limitations

- A run read back from its results file has its start rebuilt from its end, which may differ from the run's own start in the last bits (`Trajectory::from_end`).
- A set file's x other than 2θ is refused until the conversion of a set (PL10); the CSV is the wide form in 2θ (d, Q and the long form are PL7's).

## app 0.17.1 — 2026-10-13 — TU10's items: plan step gate, Reset Layout keeps tabs, the tour on its own settings



### Changed

- A tour launch keeps its settings in a folder of its own; Help ▸ Tour plays on the default settings and layout.
- The player reads a step's answer: for file dialogs in any window, the zone action, problems:fix:N as the N-th problem's Fix; it returns to its window after a window: new scene and writes a timings file.
- Convert to Another Wavelength on a scan without λ offers the λ of a control file (request 34, core 0.32.0).

### Fixed

- Run Next Plan Step was disabled until something was free; the plan's free list now also takes spaces.
- Reset Layout closed every tab; Move Tab to New Window threw on an Explore tab (now refused with the reason).
- Save Project As with the data referenced lost the reference at the next layout save; the project reopened with no file.

### Known limitations

- Take 17 in the built app: 15 of 27 scenes, 116 of 207 commands; no run starts in 3.2 to 3.6 (not yet traced); 6.1b's Start a Refinement still disabled.

## 0.32.0 — 2026-10-13 — Two refinements' paths compared; the results file keeps the correlations and the stop; the wavelength operation from a given λ

L184: Explore's E11, E11a–c and the app's request 34.

### Added

- `lsq::compare` and `PathsCompared` (E11): two runs' paths by parameter address, in common and apart, the fixed values changed between the starts, both paths on the reference's scale, the steps with the shorter run held, the difference at each step, the ends compared within `ENDS_AGREE` = 3 combined σ, a summary with both stops; `Side`, `Run`, `EndDifference`, `CompareSummary`.
- `PathsCompared::pair` and `PairCompared` (E11a): each run's path, end, σ's, correlation and covariance ellipse for two parameters.
- The results file's free parameters with value and esd (`_drs_free_*`) and their correlations (`_drs_correlation_*`, E11b), and the stop as `_drs_stop_code` and `_drs_stop_text` (E11c); `Results::free`, `correlations`, `correlation_matrix`, `stop_code`, `stop_text`.
- `drs compare <a> <b> [--reference a|b] [--pair x y] [--json]`, on results files or control files.
- `Operation::Wavelength`'s `from` (app 34): the λ of a scan that states none; in the control file and the operations' dictionary.
- Tests `paths_compared` (9) and `app34_wavelength_from_a_given_source`.

### Changed

- The results file's `_drs_stop` is now `_drs_stop_code` and `_drs_stop_text`; an earlier file's `_drs_stop` reads as the text.
- `Operation::Wavelength` has a second field, `from`: code that builds it adds `from: None` (the app's data module).

### Known limitations

- A field's su in the results file is not the esd of the parameter named after it when another parameter also moves that field; the free parameters' esd's are now kept beside them.

## app 0.17.0 — 2026-10-13 — Plots: the projection with hidden parts, the series of windows against temperature, FullProf's .seq (L183)



### Added

- Projection (F2): the members displaced along an oblique temperature axis, each drawn only where no curve in front hides it, as the core works it out (plots::project).
- Series (F5): Add Window…, Remove the Windows, Background at…; for each member and window the core's integral above the background with σ, its centroid and strongest peak, drawn against the temperature with σ bars; the windows shown as bands in the overlay and the stack.
- FullProf's .seq opened in a Plots tab: every refined parameter, χ², each phase's R_Bragg, volume and weight fraction against the external parameter or the file order.
- The manual's sections and ? topics for the projection, the series and the sequence.

### Fixed

- The stack's y axis named one member's file; it reads "counts + each member's step".
- Export ▸ Plot's figure listed every member when the legend was hidden; a linear axis forced to zero flattened a cell parameter's curve.

## app 0.16.1 — 2026-10-13 — Two questions in a row kept again; Revert leaves the file saved; the tour follows into a new window



### Changed

- The tour plays each step in the window in front, the one the tour itself opened or brought forward, through a mailbox on the Rust side (scene 6.1); a window it opens runs no self-check.
- The tour's checks read the clipboard's types from the system's pasteboard; a target not found says which handles of its kind are shown; problems:fix:N is the N-th Fix button shown.

### Fixed

- Smooth, Mask a Range and Delete a Range kept no operation: the dialog asked twice in a row was cancelled by the first question's late close event (found by the tour, scene 1.7).
- After Revert To ▸ Last Saved the Editor still marked the file unsaved, and Snapshot asked to save first (scene 2.1).

### Known limitations

- The whole tour in the built app: 3 of 26 scenes with tutorial/script.yml as it is (50 of 197 commands), 8 of 26 with its mismatches corrected (73); TU7's corrections are not merged yet.

## app 0.16.0 — 2026-10-13 — Plots, the fifth module: several data files at once, overlaid, stacked, as a map and a difference map (L180)



### Added

- The Plots module (R32, L169) on core 0.30.0: a plot set opened from several files or WinPLOTR's .buf, each file with its temperature and where it was read; Overlay (dimming, Solo Member, Dim Others, Restore All, an opacity per file), Stack (the core's displacements, Set Displacement…, labels at the baselines), Map and Difference Map (the core's rows, edges and colour range).
- Sum, Average and Difference of Members write a new file, Plots/{set} {operation}.xys, through the core, never overwriting.
- A file that does not state its temperature is drawn apart until one is typed, with its source; with no monitor every view says "no monitor: raw counts".
- The plot component: an opacity per series, series displaced with their values, labels at the baselines, an image layer with a colour bar.
- The Plots manual page and its ? topics; the no-physics guard reads every file of the module and its bridge.

### Known limitations

- The plot set lives in memory until the core writes its file (PL12): Save Plot Set is off.
- The projection and the series view come with the core's PL17, PL5 and PL11 (core 0.31.0) in the module's next step.

## app 0.15.0 — 2026-10-13 — The tour of a refinement: the tutorial's script played in the app (TU1), on stable handles (TU5)



### Added

- Help ▸ Tour of a Refinement and DRS_APP_TOUR=<script.yml>: the tutorial's script played in the window, scene by scene (commands, toolbar fields, actions, waits), each scene's checks over a read-only API, a caption bar with the narrator's line, a ring on the control in use, Pause, Step and Take Over (Escape); a report of every scene, step, check and command, shown and written as JSON.
- TU5's handles: data-tour on the plots (with plot.pixelOf), legend keys, the sidebar's files, the dividers, the app's dialogs' buttons, Data's selects and background controls, Refine's outline, inspector rows and Reflections controls.
- The tour's checks read by an interpreter of a JavaScript subset: the window's CSP refuses eval, and keeps doing so.
- The manual's section on the tour, opened by ? while a tour plays.

### Fixed

- Close Tab on the text beside a plot left no tab active: the tab used before becomes the active one.
- Undo after a click on a plot undid a text field's typing: a click on a plot takes the focus from a text field.
- Convert to Another Wavelength on a scan with no wavelength broke every later operation: refused with the reason (request 34).

### Known limitations

- The whole tour in the built app: 5 of 26 scenes pass with the script as it is (50 of 176 commands), 8 of 26 with its own mismatches corrected (62); the rest are the script's (TU7) or Explore's missing handles. Smooth, Mask a Range and Delete a Range in scene 1.7 leave no operation, not yet traced.

## 0.31.0 — 2026-10-13 — Plots (c) and (d): windows against T, FullProf's .seq as a table, the projection's hidden parts

L177 (L170, step 4): Plots' PL5, PL11 and PL17; the two findings of 0.29.0 closed.

### Added

- `plots::Window`, `WindowBackground`, `WindowValue`, `MemberWindows`, `Windows`, `windows`, `window_value` (PL5): per window and member the integral above the automatic background, given points or the line between the window's ends, with σ; the centroid with σ; the strongest peak; the background at chosen x.
- Module `seq` (PL11): FullProf's `.seq` read as a table (`Seq`, `SeqBlock`, `SeqPhase`, `SeqParameter`; `names`, `column`, `inconsistent`, `repeated`), lines not understood kept.
- `plots::Projected`, `Projection`, `project` (PL17): the visible pieces of each curve of a projection, front to back, and the oblique axis.
- `drs plots windows`, `drs plots seq`, `drs plots project`, with `--json`.
- Serde on `data_ops::FoundPeak` and `background::Interpolation`.
- Test `plots_series` (7).

### Changed

- A point of no weight is `null` in a pattern's JSON both ways: `null` now reads back as an infinite σ.

### Fixed

- A pattern's JSON with a point of no weight did not read back into the core (known in 0.29.0).

### Known limitations

- The D20 reader states no wavelength, by decision (L177): FullProf's Ins = 3 names the field holding 1.870 `tr2` and reads no wavelength.
- `plots::Resampled` writes σ as `null` both outside a member's range and for a point of no weight.
- The `.seq`'s second R_Bragg number of a magnetic phase is not identified.

## app 0.14.0 — 2026-10-13 — The manual, Help and What's New working; the whole manual written (L174)



### Added

- The manual in a window of its own: contents, sections, search, links to the site, maths in KaTeX, offline. The app's manual (the window, the project and its files, Settings, the keyboard and the workflows W1–W14 step by step), one manual per module (Data, Refine, Editor, Explore: every view, panel, dialog and field with the core's function behind it), the site's Files and Data operations pages made by the site's own writers, and the parameter dictionary; every command's table drawn from the command registry.
- ? on every panel's header and Help ▸ Manual for What Is on Screen: the module's manual at the section for its view. ⌘K finds the manual's sections.
- What's New opens by itself after an update, the releases new since the version last used marked; every release newest first by date.
- The self-check opens every viewer, manual page and Help item, and fails on a body with no height or no text.

### Changed

- A version's data read with the core's control::open_in (request 33, core 0.28.0).

### Fixed

- What's New and Keyboard Shortcuts opened empty (L174): their body was flex: 1, which WebKit keeps at zero height in a dialog with only a maximum height; flex: 1 1 auto gives it its height. Checked in the running app: 818 × 0 px before, 818 × 588 px after.

## app 0.13.1 — 2026-10-12 — The run's record is the document's: the core's Trajectory kept through Accept and read back from a version's results file (E9)



### Added

- Request 33: control::open_in(path, base).

### Changed

- The document keeps the core's Trajectory (lsq::Session::trajectory): set when a run begins, grown each cycle, kept with the pending result and after Accept, dropped on Discard; a version opened from Versions/ reads its run back from its results file's cycles (Trajectory::from_end). Explore's commands take the record from the document; no shift travels through the window (explore/REQUESTS.md E9's shell part). Their shifts argument is optional and ignored when the document holds a record.

### Fixed

- A version opened from Versions/ failed when its data file is named relative to the project's folder: the document opens it with its files resolved from the parent of Versions/ (request 33 asks the core for control::open_in).

## 0.30.0 — 2026-10-13 — Plots (b): a common grid and the map of a set, the sum, average and difference of members

L175 (L170, step 3): Plots' PL4 and PL9.

### Added

- `plots::GridRange`, `GridStep`, `same_grid`, `common_grid`, `Resampling`, `Resampled`, `resample`, `on_grid` (PL4 a): members on one grid, as measured, interpolated or bin-averaged, σ propagated.
- `plots::Map`, `MapRow`, `map` (PL4 b): rows at the variable's values with their edges, members without a value apart, x edges, the colour range on the linear, √ or log scale; the difference map less a reference member on a symmetric range.
- `plots::Derivation`, `Derived`, `derive`, `derived_text` (PL9): the sum, plain and σ-weighted average and difference of members with σ, written as XYDATA (.xys) naming its members and their temperatures.
- `data_ops::interpolate_at`, public.
- File dictionary kind `derived_data`: `Plots/{set} {operation}.xys`.
- `drs plots grid`, `drs plots map`, `drs plots derive`.
- Test `plots`: `pl4_grid_and_map`, `pl9_derived_patterns`.

### Known limitations

- Whether FullProf skips the '!' lines a derived file holds after its six header lines is not stated in its manual (not checked).

## 0.29.0 — 2026-10-13 — Plots (a): several files as one set, WinPLOTR's buffers, each member's temperature with its source, the displacements of a stack

L172 (L170, step 2): Plots' PL1–PL3.

### Added

- Module `plots`: `PlotSet`, `Member`, `Facts`, `open`, `Buffer`, `read_buf`, `open_buf` (PL1); `Variable`, `Found`, `temperature_in`, `temperature_in_name`, `temperature_candidates`, `PlotSet::set_variable`, `sort_by_variable`, `without_variable` (PL2); `YScale`, `Displacement`, `Displaced`, `displace` (PL3).
- `drs plots open` and `drs plots displace`, with `--json`.
- Serde on `data::Pattern` and `formats::Scan`.
- Test `plots` (4).

### Known limitations

- The D20 reader does not give the wavelength its title states (1.87 Å).
- A σ of "no weight" (infinite) is written as null in JSON and does not read back into the core.

## 0.28.0 — 2026-10-13 — B with its su in the results file; a control file opened from another folder; a phase refined by coordinates described by modes in the 3D scene

L171 (L170, step 1): Publish's P20 and P19, the app's 33.

### Added

- `structure_factor::b_from_u`, `u_from_b`: B ↔ U with their su's (P20).
- The results file: each atom's B_iso (B_eq for an anisotropic one) with its su, `_drs_atom_site_B_iso_or_equiv` (P20).
- `control::open_in(path, base)`: a control file whose files are named from another folder (app 33); `--base <folder>` for every command that reads a control file.
- `geometry::decompose`: a phase refined by its coordinates described by another file's modes in the same setting, with the residual (P19); `Scene::reference_fract`; `drs structure --modes <file> [--amplimodes <output>]`.
- Test `small_requests` (3).

### Changed

- `geometry::Scene` gained `reference_fract`: code that builds a `Scene` by hand adds it.

## 0.27.0 — 2026-10-12 — The open requests to the core: U_eq, the parent cell and the valence sums in the results file, a file named by its place, FullProf's exports in Exports/FullProf; a pattern that is not a number refused

L168 (L153, step 2): the requests of the app, Explore, Publish and the tutorial checked against the core; those already met named, the rest built.

### Added

- `Model::finite`: a calculated pattern with a value that is not a number refused, with its first point and the values outside their limits (TU6).
- `geometry::u_equivalent`, `u_eq`: U_eq of an anisotropic atom with its su (P15); `geometry::Valence`, `valence`: each charged atom's coordination, mean distance and bond-valence sum with su's through the covariance (P18).
- `Phase::parent`: the parent of a phase refined by its coordinates, the control file's `mode_parent_*` read and written without modes (P17).
- The results file: `_drs_parent_cell_*` for every structure phase with a parent, the `_drs_valence_*` loop (P17, P18).
- `formats::project::describe_in`: a file named by its place in the project before its content (app 31); `drs data info --in <project>`.
- `drs modes --equivalent [--reference <file>]` (E7); `--amplimodes` gives a mode phase its parent.
- Tests `open_requests` (6) and `simulate_not_a_number`.

### Changed

- The CIF and the results file write U_eq for an anisotropic atom, where they wrote 0 (P15).
- FullProf's exports are in `Exports/FullProf/` (the file dictionary; Josu, L123; app 32): `drs convert` from a control file writes there.
- `drs calc` and `drs simulate` stop on a pattern that is not a number.
- `powder::Phase` gained `parent`: code that builds a `Phase` by hand adds it.

### Fixed

- `drs simulate` never ended on a model whose pattern holds NaN (the Poisson sampler's rejection loop; TU6).

### Known limitations

- Publish's P19 and P20 are not filed yet.

## app 0.13.0 — 2026-10-12 — Start a Refinement from the data and the structure (L162)

The app stream's first step of L160, on core 0.26.0's `drs new`.

### Added

- Data ▸ Start a Refinement… from the structure: one CIF per phase, or a formula, a space group and a cell (Le Bail); the wavelength from the file, typed, or an anode's Kα lines; neutrons. Before opening, what the core wrote and every default with its source; Open or Discard.

### Known limitations

- npm install must precede npm run build (three.js).

## 0.26.0 — 2026-10-12 — A refinement started from a data file and a structure: drs new

L155 (L153, step 1): the core writes a control file from a data file and a CIF, or a formula, a space group and a cell, with the automatic background in its own file, the widths from the data's peaks, the scale by least squares, the scale and the background free and every default stated at the head of the file with its source. On Y2O3 the file refines to FullProf's minimum with FullProf's background points; with its own, the differences are frozen.

### Added

- Module `start`: `start` (a control file from a data file, `Known` phases and an `Instrument`), `Start` (the control, its text, the background file, `stated` defaults, the data's range, peaks and median FWHM) and `Start::write` (refuses to overwrite); `cell_from` (a cell from its free values in the crystal system's order), `formula_elements`, `le_bail_list`, `median_fwhm`, `least_squares_scales`; `KALPHA_RATIO`, `DEFAULT_PEAK_WINDOW`, `DEFAULT_ASYMMETRY_LIMIT`.
- `drs new <data> (--cif <file> | --phase <formula> <space group> <cell>)… [--wavelength | --anode] [--neutron] [--polarisation] [--window] [--asymmetry-limit] [--profile] [--out] [--as]`.
- Test `start` (4): Y2O3 from its data and a CIF against FullProf's run; the defaults; a Le Bail phase from a cell; the command line.

### Changed

- The parameter dictionary's `monochromator_polarisation`: without a monochromator Cthm = cos² 0 = 1 (FullProf's manual, eq. 3.26), not 0, which leaves out the polarisation's angle dependence.

### Known limitations

- With its automatic background (101 points, all refined) Y2O3 refines to χ² 1.1673 against FullProf's 1.199, with the scale, η0, X, W and the zero 1.5–8.4 FullProf esd's from FullProf's values; with FullProf's 62 points every parameter is within 1.04 esd. Frozen in `start`, for Josu to accept.
- The control file holds no plan: the steps after the scale and background are the user's to write.

## app 0.12.4 — 2026-10-12 — A Modules menu: every module with its own actions (L151)

Installed in /Applications as the one copy to open.

### Added

- Modules in the menu bar: Data, Refine, Editor and Explore, each with Show (⌘1 …) and all its actions, whichever module is active; Explore's actions are at hand from Refine.

## app 0.12.3 — 2026-10-12 — The Mac stays usable while the app works (L150)

The two loads found when Josu's Mac slowed with the app open, cut. Installed in /Applications as the one copy to open.

### Changed

- Explore's χ² map uses half the logical processors (8 of 16 here) at the background's priority (macOS QoS utility).
- Explore's 3D view asks for the low-power GPU, keeps no copy of each frame, and lets its GPU context go when another view is shown.

### Known limitations

- Measured idle and during a refinement: about 70–90 MB for the app and 100–130 MB for its window, no CPU when nothing runs. If the Mac still slows, the view open at the time is to be measured.

## app 0.12.2 — 2026-10-12 — The activity bar labelled with thin line icons; crossing dividers dragged together (L149)

Josu's three requests of L149, written into APP_CONVENTIONS §1 and §2. Installed in /Applications as the one copy to open.

### Added

- Where a vertical and a horizontal divider meet, one drag moves both (the move cursor shows there).

### Changed

- Every button of the activity bar has its name under its icon, the modules' included; the bar is 60 px wide.
- Thin line icons in the bar: the modules' at stroke 1.1, the bottom glyphs replaced by line icons of the same weight.

## app 0.12.1 — 2026-10-12 — Every calculation the app showed is the core's (L140 done, L145–L148)

The 35 calculations L140 found in the app (its window and its Rust side, borderline ones included) are moved into the core and tested there; a guard of six tests fails if any is written in the app again. Installed in /Applications as the one copy to open; the running app's self-check passed on it.

### Changed

- The app's Rust side: strong correlations, who is this peak, the data's step, largest count and σ, the calculation grid, the run's defaults, limit states, the map's half-widths, the widths' sampling, a background function's defaults: the core's (L145).
- Refine: Yobs − Ycal, the cumulative χ² on counts, the shares in per cent, seconds per cycle, the strong cells of the heat map, the restraints' deviation: the core's; the sparkline labels the first and last χ² (L146).
- Data: the largest count and the proposed prominence, the selection's range, the point a text line is about: the core's (L147).
- Explore: a profile on its background, the H_G² < 0 range, the default reflection and pair, the χ² colour bar, the structure's morph and arrows: the core's (L148).

### Fixed

- Explore's χ² map no longer shows a time per point guessed from the window's clock (L148).

## app 0.12.0 — 2026-10-12 — Every calculation in the core (L139); the five features still missing from L125: trend, notes, Dock, Help, split and windows

L139–L144. Josu's rule that the app only shows is written into APP_CONVENTIONS §0, and Explore's path in esd units moved into the core; item 7 of L125's second group is done. Installed in /Applications as the one copy to open.

### Added

- Refine's bottom panel: a free parameter across the versions, as a trend, values from each version's control file and su's from its results file (the core's `trend`, `drs trend`; W7, L141).
- Notes on reflections: a note column in Reflections and Refine ▸ Note on the Reflection…, kept in `<control> notes.toml` and carried into the results file (the core's `notes`, `drs notes`; W8, L142).
- A refinement's progress in the Dock icon and a notification when it ends with the app in the background (§9, L143).
- Help at the bottom of the activity bar: search (⌘K), keyboard shortcuts, the manual (§2, L143).
- View ▸ Split Right and Split Down; Window ▸ Move Tab to New Window, and a tab dragged out of the window opens in a window of its own (§1, L144).

### Changed

- Explore's traces, its table's shift/esd and distance columns, the χ² map's span over the path and the cycles' largest shift are the core's (`lsq::Trajectory::in_esd`, `pair_span`, `Cycle::largest_shift`; L139).

### Fixed

- Explore drew every trace starting at 0 after a run: the start was taken from the document's refined values; it is now the core's start (L139).

### Known limitations

- L140 lists 26 calculations still in the app, to move into the core one at a time.
- Not seen yet in the running app: the trend's drawing, the notification on screen, the real drag-out gesture (L144).

## app 0.11.4 — 2026-10-12 — The second group of the audit (L125): export, names, layout, physics from the core

Items 1–6 of L125's second group ("short of what was asked"), L133–L138, each guarded by a test. Installed in /Applications as the one copy to open.

### Added

- File ▸ Export ▸ Plot… writes the plot in view as PDF, SVG or PNG (300 dpi) in Data, Refine and Explore, on white whatever the window's appearance (L134).

### Changed

- One tick rule for every plot: round counts on √ and log axes, Explore's cycle axis thinned (L133).
- An unnamed phase is proposed by the core's formula: Y₂O₃, not "YO" (L135).
- Refine names parameters by the dictionary ("Cell a", "Background at 15.48°", "Atom position x (Y)"), never by their address, in the bottom panel, the Free parameters tab, the steering and the core's notes; the address on hover (L136).
- Explore's controls are a bar of their own under the view tabs, wrapping in a narrow window (L137).
- Data's d-spacings and Explore's covariance ellipse are the core's (`data_ops::to_axis`, `lsq::covariance_ellipse`) (L138).

### Fixed

- Ycal and the difference strip were missing from exported figures drawn from a dark window (L134).
- Data showed no d-spacings with a wavelength taken from a control file (L138).
- Explore's list of runs ran the cycle count into the title (L137).

### Known limitations

- Item 7 of the group is still to do: W7's trend across versions, W8's notes on reflections, tabs dragged out and split views (§1), Help in the left bar (§2), the Dock's progress and a notification when a run ends (§9).

## app 0.11.3 — 2026-10-12 — The fixes of the audit (L125): Data's operations, the File menu, shortcuts, the cycles, the warnings

Fixes L126–L131, made on main and each guarded by a test that fails on the code before it. Installed in /Applications as the one copy to open.

### Fixed

- Data ▸ Mask a Range…, Delete a Range…, Smooth…, Convert to Another Wavelength… and the file operations wait for their answers and act (L126); a test guards every dialog question of the app.
- File ▸ Organise as a Project… and Undo Organise in the File menu; Not now hides the strip for the session only (L127); a test that every File, Edit and View command is in its menu.
- No shortcut does two things: One Cycle ⌥⇧⌘R, Open as Text ⌃⌘T, Select Points ⇧⌘L, Pick Background Points ⌥⌘B, Fix Selected ⌥⇧⌘F, Discard ⇧⌘⌫; the macOS shortcuts keep their meaning (L128), tested.
- The control file's Cycles decide a run; Settings' value is a default for a file that states none (L129), tested.
- A Continue does not ask again about warnings already accepted in the run (L130), tested.
- "Cycles at most"; a converged run says "Converged after k of N cycles" and that Convergence 0 runs them all (L131), tested.

## app 0.11.2 — 2026-10-11 — Josu's try of 0.11.1 fixed (the status bar's error, results as results, ours and FullProf's apart, tabs named by module); W14, FullProf's files apart with Organise as a Project and Undo



### Added

- W14 (L123): Organise as a Project, offered in the sidebar and the File menu when FullProf's files are at a folder's top. The dialog lists every step first; FullProf's files move into FullProf/ as they are, the working .pcr is imported as {project}.drc, the data into Data/, the .irf into Instrument/, a .bac the .pcr reads into {project} background.drb; a reference folder is never organised; File ▸ Undo Organise puts everything back from .drs/organise.json. The sidebar shows FullProf/ as one closed group, FullProf (original files).
- Requests 31 (the dictionary's role in describe()) and 32 (the FullProf exports under Exports/FullProf/).

### Changed

- Ours and FullProf's told apart (L122, W5): a .drc and a .pcr, a .drb and a .bac, our results file and a .prf have different icons; every FullProf file carries an F on its icon, said in words in the tooltip.
- A file's tab carries its module's icon and name (L122): Refine and the Editor on one file are told apart.

### Fixed

- The status bar's script error with Explore's tab open (a tab without a file; L122): the status bar's words are a pure function with a node test.
- A results file (pdCIF) listed under Data (L122): the file dictionary's role by place and name comes before the content's format (request 31 asks the core's describe() to do the same).

### Known limitations

- The exports for FullProf still go where the dictionary puts them (FullProf/) until request 32 lands; in an organised project the sidebar keys FullProf/'s files by content, so an export written there shows as FullProf's.

## app 0.11.1 — 2026-10-11 — 0.11.0's plan steps, variables and links asked in the app's own dialogs; one Run with checkpoints; §8's tick and table rules



### Added

- askForm: several fields in one of the app's dialogs (text, number, several lines), used by the plan's steps, the variables and the links.

### Changed

- W11 as changed by Josu (L119): one Run, continuous by default, with 'Pause every [k]' beside it (k = 5 pauses at 5, 10, 15…; at each pause Continue, Accept or Stop); every cycle recorded, paused or not; One Cycle (⌥⌘R) stays; 'Step n' is gone; Cycles stays the maximum.
- The shared plot draws as many ticks as its labels allow (measured), takes a common offset or factor of ten into the axis label (×10⁻⁵, + 10.6048) and uses the true minus (§8, L119); the shared table's header is fixed and its first column stays in view when it scrolls sideways.
- The sidebar's files as a code editor shows them (W5, L120): an icon per kind of file keyed by the core's kinds (data, control, version, results, CIF, background points, resolution, FullProf's outputs, notes, exports, unknown), a colour per role (the group's bar and its files' icons; Okabe–Ito, light and dark), the state on the name (a dot for unsaved changes, dimmed when stale, red when missing, a lock on a version, a mark on an open file), and the same in words below the name and on hover (§10).

### Fixed

- A plan's step, a new variable and a link were asked with the browser's prompt(), which Tauri's WebKit view does not show: they did nothing in the window (L118). A form dialog of the app's own asks for all their fields at once, each with a line of help.

## drs-publish 0.1.0 — 2026-10-11 — The Publish crate's first build: the Rietveld plot and three tables from a results file, every character in LaTeX

The first build of R30's catalogue (publish/DESIGN.md): F1 and T1–T3 from a results file, read from a publication set (.drp) and written by `drs-publish build`, every figure and table typeset by LuaLaTeX in Latin Modern (L103).

### Added

- `drs-publish build <set.drp> [--out <folder>]`: every item of a publication set rebuilt from its results file into the set's folder (Publish/<name>/): figures as PGF, PDF and the preset's formats; tables as LaTeX and PDF; a preview document; a manifest with each file's SHA-256.
- F1, the Rietveld plot (observed, calculated, difference, ticks, an inset; a monochrome variant); T1 crystal data and refinement, T2 atomic parameters, T3 distances (preview).
- Presets as data (`presets.yml`): IUCr, ACS, RSC, Elsevier, APS, UPV/EHU TFG, TFM and PhD, slides.
- The check of L103 as a test: every PDF and EPS must list Latin Modern faces only.

### Known limitations

- Josu's §10 decisions recorded: no Word tables (the LaTeX source for journals that want an editable table); a table in a document not set in Latin Modern switches its own text and needs Latin Modern Math; HTML and Markdown are data exports.
- The other items of the catalogue come in later builds.

## app 0.10.1 — 2026-10-11 — Every question in the app's own dialog: Overlay a Version…, Snapshot…, Rebin…, Normalise… and the others act again



### Fixed

- window.prompt and confirm show nothing in Tauri's WebKit view, so Overlay a Version…, the Snapshot… note (Refine and Editor), Rebin…, Normalise…, the wavelength for a conversion and the fitted function were silently cancelled (Josu, 0.10.0). The app's own dialogs ask instead: a text, a number, or a choice among options; Overlay a Version… lists the versions with their time, trigger, χ² and note.

## app 0.11.0 — 2026-10-11 — W11 and W12: the refinement paused and resumed with the engine's state, every stop the core's reason; the steering, the restraints' shares and the plan run step by step

L117, on core 0.25.0 (lsq::Session, steering and plan of core 0.24.0): a run is the core's session, run without the documents' lock, paused and resumed to the bit; W12's interface shows and edits what the core builds from the control file, and the plan runs one step at a time, each accepted step a version.

### Added

- Pause (⌥⌘.) and Resume (⌃⌘R): the session kept in the document while paused, resumed with its λ and cycle count; Step continues a paused run; Stop ends a paused run where it is.
- During a pause: freeing, fixing and the steering changed, the session taking them up on Resume with its damping restarted (the core's note).
- The stop reasons from the core's codes in plain words, with a singular matrix's combination, why and what to do; per cycle λ, the core's time, the refused steps and the steering's notes.
- Steering in the inspector: λ's start, Marquardt, the classes' damping, each free parameter's state, damping, bounds, behaviour, maximum shift and link; links and variables added, edited, removed; every change written to the control file by the core.
- Restraints in the inspector: target, value, Δ/σ, weight and share of χ², and their part; the part in every cycle and result.
- The plan: steps added, edited, moved, removed; Run a step or the next; each accepted step a version naming it.
- Tests: pause and resume against the uninterrupted run; a change during a pause; the steering kept in the control file; the plan step by step; the restraints' shares (32 tests).

### Changed

- A run no longer holds the documents' lock: Stop, the view and every command reach the document while it runs.
- The background file of a data file alone is the file dictionary's "{stem} background.drb", named by the core (request 29); the forbidden reflections are the core's diagnostics::forbidden_in (request 30).
- The window's title-bar style is set on macOS only, so that the tests build on Linux.

### Fixed

- Stop during a run: until 0.10.0 it waited on the lock the run held and could not reach it.

### Known limitations

- The steering's fields are typed (no sliders); dialogs for plan steps, variables and links are the system's prompts.
- A run's progress is not in the Dock, and no notification is sent when a long run ends while the app is in the background.

## 0.25.0 — 2026-10-11 — What Publish and Explore asked of the core: a refined state's geometry with its su's, one su rule, the trajectory of a run, χ² maps, a state's bundle, the 3D scene, the signs of mode amplitudes

L113–L116: FullProf's half step on the cell traced and frozen (a factor of its own, ½ with the pseudo-Voigt, ¼ with TCH); the Publish stream's P2–P11 (the results file and the CIF complete, positions and bonds of the refined state with propagated su's, charges, one su rule rounded up, a phase's name, the parent cell); the Explore stream's E1–E7 (a run's trajectory kept, persisted and rebuilt to the bit, χ² on a grid, a state's bundle, the scene of the 3D view, the eigen-spectrum, an irrep's displacement, equivalent amplitude sets); Josu's decisions on the su rule, the background file's name and Co-YBCO.

### Added

- Module `geometry`: `Propagator` (every su of a derived quantity through the covariance), `structure`, `unique_labels`, `geometry` (positions reduced to the cell, bonds once, angles, with su's), `parent_cell`, `contacts_with_su`; `scene` with `Image`, `Face`, `Polyhedron`, `SceneRule` and `hull` (E4).
- Module `state`: `summary` (a state's bundle: cell, volume, strain, widths with the unclamped H_G², reflections, contacts, valence sums, modes and irreps with su's, breadths, restraints), `reflection_profile`, `plausible` (bond-valence and ionic-radius distances) (E3).
- `lsq`: per cycle the covariance, the spectrum of the scaled normal matrix and the time; `Trajectory` with `state_at` and `from_end`, `state_at`, `correlations_of`; `Session::start` and `trajectory`; `Refinement::history` and `shift_su_max`; `chi2_grid`, `Grid`, `coarse_to_fine` (E1, E2); background-free Rp, Rwp, Rexp in `Statistics` (P2).
- `check::directions` and `Direction` (E5); `Modes::displacement_of`, `irrep_displacement` (E6), `equivalent_amplitudes`, `same_structure`, `nearest_to` and `Equivalent` (E7).
- The results file: R_F per phase, the background-free R's, `_refine_ls_shift/su_max`, `_refine_special_details`, the mode amplitudes and irreps with su's, the reference structure, the parent cell, every cycle's loops (`_drs_cycle_*`, read back into `Results::cycles`); `results::trace_json`.
- The CIF: `_pd_block_id` and the diffractogram's, `_pd_phase_name`, `_computing_*`, `_refine_ls_R_factor_all`, `_refine_ls_R_Fsqd_factor`, `_refine_ls_number_constraints`, the weighting scheme, `_refine_ls_shift/su_max`, `_refine_special_details`, the profile, background and preferred-orientation functions; `cif_write::phase_name`, `deposition_in`.
- `Atom::charge` and `charge = q` in the control file; `bonds::atom_ion`, `with_charges` (P5).
- `drs refine --trace <file.json>` and `--map <a> <b> [--span] [--grid]`; `drs state`; `drs structure --cartesian`; `drs bonds` on a .pcr or control file at its refined state (`--cycles`, `--charges`); `drs modes --amplimodes` prints the summary per irrep.
- The file dictionary's placeholder {stem} (L116).
- Tests `fullprof_cell_step` (2), `publish_requests` (6), `explore_requests` (5).

### Changed

- Every su is printed by one rule, `cif_write::with_su`, rounded up (IUCr Rec. 12, Josu's decision): 2.2865 ± 0.0032 is 2.287(4); `drs bonds` on a CFL and the weight fractions included (FullProf's print under `--fullprof`).
- The CIF's positions are reduced to the cell with su's from the covariance (a mode phase's through its amplitudes), its bonds listed once and its labels unique; values fixed by symmetry carry no su.
- A phase's CIF and data block are named by `phase_name` (Y2O3, Sr2ScSbO6), not the .pcr's line.
- `bonds` takes angles by atan2 in Cartesian axes: 180° at an inversion centre is exact, without su.
- The background points' file is "{stem} background.drb" (the app's request 29).
- `Cycle` and `Refinement` gained fields (covariance, spectrum, seconds; history): code that builds them by hand adds them.
- The committed results file of case A and GSAS-II's reading of it, regenerated.

### Fixed

- Positions written outside the cell in the CIF (O11 −0.77349).

### Known limitations

- FullProf moves the cell by ½ (pseudo-Voigt) or ¼ (TCH) of the Gauss–Newton step; the core takes the whole step. Reproduced by damping the cell, counted and frozen (Egypt not reproduced by one fraction), for Josu (L113).
- Case A's Bond_Str standard deviations: 40 of 132 above Bond_Str's formula by up to 0.32 of a last digit; unexplained, frozen, for Josu (L114).
- A state rebuilt from a run does not rebuild a Le Bail phase's extracted intensities.

## app 0.10.0 — 2026-10-11 — W7's diagnostics and W3: the difference in σ, the cumulative χ², who is this peak, the correlation heat map, the parameters' paths, two versions overlaid; a refinement started from a previous control file

On core 0.23.0: every diagnostic of W7 drawn from the core's numbers, and W3's start from a previous control file through the core's rebase, with what does not fit shown before the copy opens. Pause and Resume (W11) and W12's interface wait for core 0.24.0.

### Added

- Legend keys for (Yobs − Ycal)/σ in the strip and the cumulative χ² in the frame (diagnostics::cumulative), off until clicked; the legend's choices kept across redraws.
- 'Who is this peak?': a click on the pattern opens a card with the background's and each phase's share at the nearest point (diagnostics::shares_at) and the reflections reaching it.
- After a run: the correlation matrix as a heat map with the pair and its r under the pointer; each parameter's path over the cycles with its esd band and last shift (cycle events carry every parameter's value, esd and shift).
- Overlay a Version…: another version's Ycal calculated by the core at the pattern's points, drawn over the plot; from the Refine menu or a version's context menu; Clear Overlays.
- Start a Refinement… (W3) in Data: a .drc or .pcr copied into the project and pointed at the data by drs::rebase; what does not fit listed with suggestions before Open or Discard; the copy opened in Refine and the Editor.

### Known limitations

- W3's proposal of control files found on the disk (Spotlight), a phase from a CIF and a control file from nothing come later.
- Pause and Resume with the engine's state, the stop reasons as data, links, bounds, damping and the plan wait for core 0.24.0 (REQUESTS.md 23–28).

## 0.24.0 — 2026-10-11 — The refinement as a session, steered from the control file, with a plan of steps

L107–L110: W11's engine (Step n, pause and resume with λ and the history, a session a script resumes, every stop and divergence as data) and W12's steering, all in the control file with the reader's checks (damping per parameter and class, maximum shifts from the dictionary, λ's start, bounds that stop, reflect or wrap, links with an offset, new variables, a plan of steps); a .pcr's multipliers imported by their meaning and its relaxation factors as the classes' damping, checked against FullProf's shift_control runs parameter by parameter.

### Added

- `lsq::Session`: `step(n)`, `run`, `set_value` and `set_parameters` during a pause (λ restarted), `finish`, `to_text` and `resume`; `refine` and `refine_until` run a session to its end.
- `lsq::Stop` and `Divergence` with codes and texts: converged, minimum, cycles, diverging (χ² rising k cycles, a jump of many esd's, an impossible value: `lsq::impossible`), singular (the combination named as `check` names it, `check::null_direction`), user.
- Per cycle: each parameter's shift/esd and computed shift; λ; the steps refused (`Try`); the restraints with weight and contribution to χ² (`restraint_states_at`); notes.
- Module `steering`: `Class` and `class_of` (FullProf's relaxation classes), `Behaviour` (stop, reflect, periodic), `Link`, `Variable`, `Guards`, `Steering`, `PlanStep` and its selectors.
- `refinable`: `Free::damping` and `value`; `ParameterSet::behaviour`, `links`, `steering`, `enforce_links`, `describe` (each free parameter's class, damping, coefficients, bounds, behaviour, maximum shift and link, for the app); `Refinement::for_step`; the occupancy's upper bound, a full site.
- The control file: `lambda`, `marquardt`, `damping`, `divergence`, `steering` and `behaviour` in `[refine]`; `[[refine.link]]`, `[[refine.variable]]`, `[[refine.step]]`; each checked as it is read.
- The parameter dictionary: each refinable entry's FullProf class (with its source: the manual, measured, or believed) and its largest shift per cycle (measured, L108), on the Parameters page and in `parameter_table`.
- `pcr::Pcr::relax` (R_at, R_an, R_pr, R_gl); `lsq::fullprof_esd` (FullProf's a·σ, for comparisons); `Model::gaussian_h2`.
- The app's requests 23–28 (W11, W12), with `Session::paths` (25) and `lsq::restraints_part` (27); 30: `diagnostics::forbidden_in`, used by `drs reflections --forbidden`.
- `drs refine --session <file>` (Step n with --cycles, resumed from the file), `--plan`, `--step <k>`, `--parameters`; the restraints printed with weight and contribution.
- Test `steering` (13): shift_control against FullProf, the import, damping and esd's, maximum shift, bounds, links, variables as mode amplitudes, the guards, a resumed session, the plan, restraint weights, the reader's checks, the command line.

### Changed

- A .pcr's lone multiplier is the damping of its parameter (a change of variable before, with no effect); the relaxation factors are read and damp the classes, so mode case A (0.8) and LaMnO3 3T2 (0.85) now move as FullProf moves them. The .pcr export writes both back.
- A singular normal matrix stops the run with the combination the data cannot fix (left out of the step before; `divergence = { singular = false }` keeps the old way).
- Convergence is judged on the computed shift, before damping and the maximum shift.
- A lone value's coefficient is ±1 after the import's translations (the same set however the file wrote it).
- The committed results file of case A and GSAS-II's reading of it, regenerated (the cell 10⁻⁵ Å away, a tenth of its su).

### Fixed

- A value the import's translations give a row without a codeword of its own was left out of its parameter: H/L when D_L < 0.0002 and S_L is refined (Co-YBCO's Cocus-asym), P1 when only P2 is coded (L108).

### Known limitations

- The first cycle of FullProf's shift_control runs (a 130 esd off) is not FullProf's: 6 of 7 comparable shifts differ by more than 10 %, up to a factor 4.29; frozen, listed for Josu with the FullProf run that would separate it (L109). In the linear range the core's step returns a to −0.04 esd.
- One cell for two patterns waits for multi-pattern refinement.
- FullProf's LINE 47 bounds (Nre, IBOUND) are skipped on import: the four .pcr files of the material that use them do not import for other reasons.
- A saved session does not hold a Le Bail phase's extracted intensities: they are extracted again at the next cycle.
- The app still runs `refine_until`: Step, Pause and Resume, the stop reasons and the steering are the core's side of requests 23–28. Request 29 (a background file beside a data file alone) is left for Josu's naming (L110).
- Co-YBCO's frozen difference between predicted and refined correlations re-measured with the corrected axial pair: 0.03323, from 0.03169 (L108); for Josu to accept.

## app 0.9.2 — 2026-10-11 — The Rietveld plot laid out as WinPLOTR; the panels always recoverable

The difference curve gets a strip of its own with its own scale, and the ticks a band, so a zoom never squashes the pattern (Josu, 0.9.1); three toggles in the toolbar show or hide the panels, and a hidden panel always comes back.

### Added

- Panel toggles (Project, bottom panel, Inspector) at the right of the toolbar, lit when the panel is shown.

### Fixed

- The difference was drawn in the frame's scale with an offset over the whole pattern: zoomed in, the pattern became a sliver. It is now a strip under the tick band, scaled over the points in view, with a zero line; the ticks sit in rows of a band under the frame.
- A hidden Inspector could not be brought back in every case: View ▸ Show … and the new toolbar toggles add the panel back beside the main group, or beside any group, or alone.

## app 0.9.1 — 2026-10-11 — Work on a Copy made safe: no copy into itself, no write into a reference folder, the projects folder checked

A serious bug of 0.8.1 found by Josu: a reference run copied into itself 85 levels deep, and files written into a reference folder. One guard now refuses every write into a reference folder; the copy lists its source before writing and refuses a destination inside the source or a reference folder; the projects folder is never inside either and is reset at launch when it is.

### Added

- Tests on temporary copies with a temporary git repository standing for manuals/: a copy into the source or its subfolder refused with nothing written; a projects folder inside a reference folder refused; a reference folder opened as it is, then refined, accepted, saved, snapshotted, logged, its results and exports and a background file written, the editor and the data operations: the folder listed the same before and after; New Project inside a reference folder refused.

### Fixed

- copy_to lists the source before writing, leaves symbolic links out, refuses a destination equal to or inside the source (canonical paths) or inside a reference folder, and removes a copy that comes out deeper or larger than its listing.
- projectsDir is never set inside a reference folder or a project; New Project proposes ~/Documents/diff-ref-struct and refuses, with the reason and a choice, a place inside either; a stored folder that breaks the rule is reset at launch with a message.
- One guard (project::guard_write) used by every write into a project — state, imported .drc, save, versions, refinement log, results, exports, background files, the editor, data operations, new projects, the old snapshots' migration — refuses a reference folder; a data file opened on its own inside one is read-only.

## app 0.9.0 — 2026-10-11 — The background by hand, the Reflections view, excluded regions on the plot, the cycles in view

On core 0.23.0 (APP_WORKFLOWS.md W4, W7, W8, W11): background points picked on the plot and kept where the user chooses; a Reflections view on the core's partitioned intensities, forbidden reflections on demand, FullProf's .prf and .hkl in the same view; excluded regions dragged on the plot; cycles and convergence in the toolbar, Step n, every stop explained and logged, only accepted results as versions, the results file written.

### Added

- Background points picked on the plot (W4) in Data and Refine: click adds, drag moves, Delete removes; a live list in the sidebar; the curve through them linear or cubic spline (the core's); kept in a .drb file beside the data (Data) or named by the control file (Refine, the file dictionary's name) or in the control file itself; synced both ways with the file opened as text; compared with the automatic background and a function fitted by the core.
- The Reflections view (W8): h k l, multiplicity, d, 2θ, FWHM, η, |F|², I calc, I obs partitioned by the model, σ, (I obs − I calc)/σ from diagnostics::reflections; filter by phase and by the plot's view; row ↔ tick; forbidden reflections on demand with their conditions; FullProf's .prf and .hkl opened read-only in the same view (request 17, core 0.23.0); Export ▸ Reflections and Profile as CSV through the core.
- The results file (pdCIF, R28) and one structure CIF per phase written by the core into Results/ when a result is accepted, with the version's time.
- Excluded regions drawn, added by a drag, moved by their edges and removed on the plot (W7), written to the control file.
- Cycles and Convergence as toolbar fields and inspector rows synced with the control file's [refine]; Step n; 'cycle k of N', the χ² sparkline and the time per cycle; every stop in plain words and in the refinement log; only an accepted result becomes a version (W11).
- The plot component: handles (draggable points), pick and region modes, editable bands, tick clicks, a view callback.

### Changed

- The snapshot before every run is gone: runs, steps and stops add nothing to Versions/ (W11).

### Known limitations

- Pause and Resume with the engine's state, the divergence reasons as data, each parameter's path, links, bounds, damping, the restraints' shares and the refinement plan wait for the core (REQUESTS.md 23–28).
- The background file of a data file without a control file is named <data> background.drb, a name the file dictionary does not give yet (REQUESTS.md 29); the forbidden reflections' d-limits repeat the command line's two lines (REQUESTS.md 30).
- The cumulative χ², the weighted difference, 'who is this peak', the correlation heat map and two versions overlaid come with 0.10.0 (the core has them since 0.23.0).

## 0.23.0 — 2026-10-11 — The mode cases agree with FullProf; FullProf's .prf and .hkl read; the reflections view's numbers; the results file in pdCIF

L098–L101: the 1–6 % between the core and FullProf on the mode cases B, D and E found and fixed (the reflection limits, D's λ/2 widths, FullProf's background below its first point), with single points at excluded-region edges explained and frozen; the sixteen new .prf files indexed; public readers of FullProf's .prf and .hkl (app request 17); the core's side of W7–W8 (app requests 19–21); the results file written in pdCIF as Josu chose, read back and read by GSAS-II, with one structure CIF per phase; the file dictionary's names decided.

### Added

- `formats::fullprof_outputs`: `read_prf` (both layouts, the background with FullProf's offset added back, propagation vectors, Zero, SyCos, SySin) and `read_hkl` (blocks by the column names FullProf prints; the Irf = 2 lists); `drs fullprof <file.prf|file.hkl>`; test `fullprof_outputs` over the 188 `.prf` and 168 `.hkl` of `manuals/`.
- `diagnostics`: per reflection d, 2θ, FWHM, η, |F|², I calc, I obs partitioned by the model, its σ and (I obs − I calc)/σ (`reflections`); each phase's and the background's share of Ycal at a point (`shares`, `shares_at`); the cumulative χ² and √w(y − yc) (`cumulative`, with `fraction` and GSAS-II's scaling); test `reflection_diagnostics` against five FullProf `.out`.
- `reflections::forbidden` and `absence_condition`: the absent reflections in a range with the condition each breaks ("h0l: l = 2n").
- `modes::ModeParent` (parent group, P, p) with `reflection`: a mode phase's reflection in the parent's cell and its k-vector; `Analysis::mode_parent`; the control file's `mode_parent_group`, `mode_parent_basis`, `mode_parent_origin`; test `symmetry_modes::parent_reflections`.
- `results`: the results file (pdCIF: overall, phase and pattern blocks; DDL1 names, the rest under `_drs_`) with `write`, `read`, and the CSV exports `reflections_csv`, `profile_csv`; `drs results [--cycles n] [--project folder] [--time t] [--exports]`; test `results_file` with GSAS-II's reading (`tests/data/make_gsas2_results.py`).
- `Phase::uvw2` (U, V, W of λ2 for FullProf's Ratio < 0, the control file's `lambda2_uvw`); `Convention::background_edge` and `Model::fullprof_step`.
- `drs reflections <file> [--forbidden]`; `drs calc --cumulative --shares`; `formats::files::now_text`; `cif::value_su`; `lsq::weights` public.
- File kinds `reflection_notes`, `export_reflections`, `export_profile`.

### Changed

- Σw(y − yc)² and χ² leave out the points the calculation leaves without peaks (`Model::excluded_points`, FullProf's rule under its convention), as FullProf's `.out` counts them (L098).
- The file dictionary: every kind decided (L099); `results` is "{control} {time} results.cif"; the {results} placeholder removed.
- The CIF reader keeps text fields as values.
- `symmetry_modes`: the mode cases compared through FullProf 8.20's runs (`mode_runs_against_fullprof`), the B parents' check kept in `cases_b_d_e_patterns`; `fullprof_profiles` reads the `.prf` with the library's reader, compares 63 files, takes FullProf's step and the Zero after the cycle from the `.out`.
- A `.bac` written at the data's points holds the model's background there.

### Fixed

- FullProf's Ratio < 0 (D): |Ratio| as the λ2 weight and its own U, V, W (negative λ2 peaks before).

### Known limitations

- Frozen with their reasons (L098): B's two and E's one excluded-region edge points, E's out-of-order region (peaks inside it in FullProf); the window-edge points of A; B's ticks past the data near 180°; NdNi-HTO's .hkl against its .prf ticks.
- I obs against FullProf: rutile-anatase (TCH with axial divergence) lower in the core for 9 strong reflections, up to 0.13 of I calc, not resolved; FullProf's σ of I obs differs by 2–6 % (its formula is not in the manual).
- FullProf's Sent0 is not implemented.
- The time in results names is UTC from the command line; the app passes its local time.

## app 0.8.1 — 2026-10-11 — The audit's fixes: tabs by file, FullProf's copies apart, the plot first, Work on a Copy

The small fixes of the audit of 0.8.0 (L093; APP_WORKFLOWS.md W7, W9 and the W5 fixes), on core 0.22.0. The reflections proposal (DESIGN.md §10) stands as amended in W8; requests 17 and 19–22 ask the core what the Reflections view and the diagnostics need.

### Added

- Work on a Copy (W9): a folder inside a git repository or read-only is offered as a copy into the projects folder, or opened as it is with nothing written (no .drs/, no imported .drc, no version, no save, no data operation kept).
- 'Use λ from <control file>' in the Data inspector when the data file states no wavelength and a control file of the project does (W7).
- The inspector's summary of the pattern (data, instrument, background, agreement, phases) when nothing is selected (W7).
- Settings: Draw Yobs as a line (markers by default, W7).

### Changed

- Tabs are named by their file, the control file's title is a subtitle; a phase without a name gets one proposed from its elements or its space group, marked proposed (W7; REQUESTS.md 22).
- FullProf's copies (_calc, _conv) and kept originals (.pcr.in) are a closed group 'FullProf copies', not working control files; their inputs are not checked (W5 fixes).
- The bottom panel stays folded until a refinement runs; the plot takes at least half the window (W7).
- The rail's bottom icons are labelled (W7).

### Fixed

- No red 'missing' for a data file named without its extension (y2o3 → y2o3.dat, as FullProf reads it).

### Known limitations

- The composition formula for a proposed phase name (Y₂O₃) waits for the core (REQUESTS.md 22); the app proposes the elements.
- W7's excluded regions on the plot, the diagnostics and 'who is this peak', and W8's Reflections view come with 0.9.0 and 0.10.0.

## 0.22.0 — 2026-10-11 — Refinement by symmetry modes; the file dictionary; the results file compared

L090–L092: a phase described by a reference structure and the amplitudes of its symmetry modes, refined and checked against FullProf's mode refinement and AMPLIMODES's analyses (R29, step 1); every kind of file of a project named in a dictionary (R27); pdCIF, NeXus/HDF5 and a format of our own compared for the results file (R28), for Josu to choose.

### Added

- `modes`: a phase's reference positions and modes (`Modes`, `Mode`, `ModeVector`); positions, displacements and their derivatives; irreps with amplitude √(ΣQ²), polarisation vector and esd from the covariance; the Gram matrix over the primitive cell in a reference metric (`gram`), projection of a structure onto the modes (`project`), `transform_cell`, `displacements`; `read_amplimodes` for AMPLIMODES's text output.
- `Phase::modes`; `Field::Mode`, refined with the chain rule through the coordinates; FullProf's Max_Amplitude as a limit; `ParameterSet::field_covariance`.
- Dictionary entries `mode_amplitude` and `mode_vector`; the control file's `[[phase.mode]]`, `mode_reference_cell`, `mode_max_amplitude`, `mode_fst`; the `.pcr` export of a mode phase (Jbt = 6, V_MODES, A_MODES).
- `drs modes` (a phase's modes, irreps and positions, the normalisation, a structure projected onto the modes; AMPLIMODES's analyses); `drs refine` reports each irrep's amplitude, esd and polarisation vector.
- `ecosystem/data/files.yml`: 18 kinds of file of a project, generated into `formats::file_table` and the page "Files"; `formats::files` (`path`, `kind_of`, `time_text`); `drs files`.
- Tests `symmetry_modes` and `file_dictionary`; case A of the mode course in `fullprof_profiles`; `manuals/fullprof/runs/modes_runs.py` for FullProf on the Mac.
- The page "The results file: the choice" (R28).

### Changed

- The `.pcr` reader keeps each mode vector's irrep, the amplitudes' codewords, the A_MODES line's FullProf Studio integers and Max_Amplitude; the amplitudes are refined from their codewords.
- `drs convert` without `--to` writes the dictionary's names (`{stem}.drc`; `FullProf/{stem}.pcr`).
- A file a `.pcr` names is also found in any case (FullProf ran on file systems that ignore case).

### Fixed

- `fullprof_profiles` takes the `.prf`'s scale from a point with counts in both files (a first point of 0 counts gave a scale of 0).

### Known limitations

- Modes are read from FullProf and AMPLIMODES, not yet from ISODISTORT, nor computed from the structures (R29, later steps).
- Run with FullProf on the Mac (L093): case B's two parents agree as predicted; the core's Σw(y − yc)² differs from FullProf's by 1–6 % on cases B, D and E (not explained); sixteen new .prf files not yet in fullprof_profiles' index.
- Case A's profiles: 404 window-edge points (TCH, X ≠ 0) up to 1.2·10⁻³, frozen as in L035; FullProf prints Yobs 1 at two points of 0 counts (excluded).
- The file names marked proposed, and the results file's format, are Josu's to accept or choose (L091, L092).

## app 0.8.0 — 2026-10-10 — The first contact: each project in its own window, the data file first, the files by role

Redirected by Josu after trying 0.7.0 (L075; APP_WORKFLOWS.md W1, W2, W5, W6): a data file opens on its own in a window of its own; projects never close each other; the sidebar lists the files by role, each named by the core; versions live in a visible, locked Versions folder with a refinement log; any file opens as text.

### Added

- Windows (W1): each project in its own window of the one process; Open Project… (a folder pick) never closes a project; the Window menu by project name; the project pop-up in the title bar with each window's state (refining, unsaved changes, errors); ⌘K lists the open projects; ⌘` cycles; each window's layout restored; the welcome page when nothing is open.
- The data file first (W2): File ▸ Open Data File… (⇧⌘D, toolbar) with filters from the core's registry by readable name, 'All data files the app reads' first; the file opens as an untitled project in its own window, not moved; Save Project As… asks where and whether the data are referenced or copied; the first look: header, pattern with σ; the peaks and the automatic background are suggested in a strip above the plot and explained in their dialogs (prominence, significance; smoothness), never done unasked (Josu, first try).
- The sidebar by role (W5): Data, Control, Versions, Inputs (missing in red), Outputs (stale dimmed), Other; a switch to the folder as on disk; hover information for every file (file_info); Quick Look on Space; Open With ▸ and Open as Text in the context menu.
- Versions in Versions/ (W5): '<control name> YYYY-MM-DD HH.MM.SS.drc', each file exactly as it was, locked; metadata in .drs/versions/; hover shows what changed since the version before; 'Versions/<control name> refinement log.md' with one entry per step; the .drs/history snapshots of 0.7.0 migrated.
- The look after Josu's try: a wider grey rail with the tools, grey dividers 6 px wide between the parts, every panel detachable into a floating window inside the app window (⧉ on its header) and docked back (⊡); the window drags by its title bar; a parameter's symbol marked in colour in its formula on its card and in the manual.
- Any text file in the Editor (W6); a data file opened as text beside its plot, the selections kept together (point ↔ line by the line's first number).
- Kind, name and role of every file from the core (formats::project::describe, core 0.20.0; request 14 closed): a test opens every run of manuals/fullprof/runs and every FullProf tutorial as a project (on copies), checks the roles against FullProf's manual where it leaves no doubt, and counts the files the core does not know, frozen at 15 of 2698 (527 with core 0.17.0). The data file is credited to the working control file, the one FullProf's outputs belong to, never to a _calc or _conv copy.

### Fixed

- New Project… is a folder pick and a name, never inside the open project (L075); File ▸ Open Project… (⌘O) is a folder pick.
- A second window (a project opened while one is open) drew nothing: the capability covered the main window only; now every window (w*). A page's errors go to the status line and the launch log.
- .fou, .pcr.in, .inp and the other files that were named 'Table of columns' or nothing: named by the core or 'unknown', never guessed (L075).
- A .pcr's data file named without extension is read as <name>.dat, as FullProf does.

### Known limitations

- 15 of 2698 files in the runs and tutorials are unknown to the core (.mlf, .ort, one .cif, .spr, .py, .zip); the sidebar says 'unknown'.
- macOS's window tab bar was tried and set aside: with the toolbar in the title bar it covered it (Josu, first try); the windows are plain, switched by ⌘`, the Window menu, the pop-up and ⌘K.
- Point ↔ line sync works for one-point-per-line formats only; FullProf's free format waits for the line of each point from the core (REQUESTS.md 15).
- File associations (double-click in Finder) and panels detached into separate macOS windows (they float inside the app window today) are still to come.

## 0.21.0 — 2026-10-10 — The data operations' settings in a dictionary; a data file written; the line of each point; FullProf's background settled on eleven variants

L087–L088: every background variant FullProf ran is within its printed precision, and the χ² differences are explained and frozen. Three of the app's requests are answered (12, 15, 16).

### Added

- `ecosystem/data/operations.yml`: the ten data operations and two analyses, with 30 settings (name, unit, meaning, range, default, proposed value, raising and lowering, source); `operation_table.rs` and the page "Data operations: the settings", both written by `build.py`; `data_ops::settings`, `keys_of` and the `DEFAULT_*` constants.
- `drs data settings [<operation>]`; `--help` after `drs data op <operation>`, `peaks` or `background` prints that operation's settings.
- `formats::write` (xys, xye), `format_of`, `NO_WEIGHT`; `drs data read|op|apply --out f.xye|f.xys` writes that format.
- `formats::point_lines`: the line of each point of a scan in its file (one point per row, or FullProf's free format).
- Tests `operation_dictionary` and `data_write`.

### Changed

- A data operation refuses a key that is not its own (another operation's key was accepted and ignored).
- `formats::read` takes σ ≥ 10³⁰ as no weight; the `.pcr` export writes its processed data through `formats::write`.
- `background_points` reports every variant: the quarter of the least background within single precision (2⁻²⁴, L087); twelve χ² differences with FullProf frozen (at most 0.00105).
- `fullprof_profiles` indexes the ten new background `.prf` files.

### Known limitations

- χ² differs from FullProf's by up to 0.00104 (0.08 %) on the background variants: the Bragg part of Ycal is about 10⁻⁴ off, the same in every variant (frozen; for Josu to accept).
- `point_lines` gives no line for the ILL formats (a header to each block), GSAS RALF (x in its own units) or a TOF XYDATA file.
- xye cannot hold y ≤ 0: GSAS-II's reader reads it as y = 0 with no weight; xys keeps it.
- `app/REQUESTS.md` is the app stream's: requests 12–16 are closed there when it merges main.

## 0.20.0 — 2026-10-10 — FullProf's background read; drs check with its data; every file of a project named; distances, angles, bond valences and a CIF for deposition

L081–L084: the background rules found in FullProf's own outputs, a `.pcr` checked with its data, the files of a FullProf project by what they are (W5), and the first outputs of phase 2 (P12, the CIF), each against FullProf's.

### Added

- `bonds`: distances and angles with standard deviations, bond-valence sums with CrysFML's table, distortion and instability indices; Bond_Str's rules kept as `bond_str_*` (occupancies normalised to the first atom, no cell term, ends independent); `drs bonds`.
- `cfl`: FullProf's CFL structure file with the standard deviations it prints; `bv_table` (generated by `make_bv_table.py` from CrysFML).
- `cif_write::deposition` and `drs cif`: a CIF of a refinement for deposition, with the engine's uncertainties.
- `formats::project`: 73 kinds of a FullProf project's files with role and source; `describe` (kind, data, or unknown); `drs data info` uses it.
- `drs check <file.pcr> [<data>]`; `check::check_path_with_data`, `check_pcr_with_data`; the hint `data-elsewhere`.
- `background::write_bac` (FullProf's first format), `write_bac_at`, `bac_text`, `write_bac_pairs`; `drs background convert --data`.
- Tests `bond_str`, `cif_deposition`, `file_kinds`; `background_points::fullprof_bac_pairs`; `check_diagnostics::pcr_data`; new variants in `runs/background_runs.py`.

### Changed

- A spline read from a `.pcr` (Nba < −5) has FullProf's own point at the data's last point; the export drops it again.
- A `.bac` is written in the first format at the data's points (FullProf reads no other); `Control::background_file_text` takes the data's points; the `.pcr` export writes the points in `<stem>.bac` (Nba = 1) when the control file keeps them in a `.bac` or they exceed 9999; `drs check` warns of a `.bac` in pairs (`bac-format`).
- `background_points::fullprof_interpolation` bounds the background by FullProf's printed precision (L081) and checks the offset, a quarter of the least background.
- A FullProf Occ is read as atoms per primitive cell when the documented reading puts an atom above 1 and this one does not (bond valence, CIF, the occupancy check).

### Fixed

- `drs data info` named FullProf's logs listings, and missed most `.prf` and `.hkl` files.
- The model's text reads floats exactly (serde_json `float_roundtrip`): one value of 13091 came back one ulp off.
- `fullprof_metric` and `fullprof_profiles` account for FullProf's background runs (L079).

### Known limitations

- Run with FullProf on the Mac (L085): `bacfirst` gives χ² 1.19 and `bacpairs` 19.9, as predicted (the 1.80 of L079 came from a hand-written file); the offset rule fails for Nba = −31 by less than the fourth decimal; ten new `.prf` files are not yet in `fullprof_profiles`' index.
- Nba = −6: χ² 0.0010 below FullProf's with the background within its printed precision (frozen).
- LaMnO3ref: Rp and Rwp 0.22 % above FullProf's CIF at its `.pcr`'s values (frozen); LaMnO3's cell uncertainties up to 2.25 times FullProf's, which refined symmetry modes (frozen).
- Bond valence needs charged species ("O-2"); a neutron `.pcr` gives elements only.
- The app does not use `describe` yet (request 14, W5).

## 0.19.0 — 2026-10-10 — Background points in a file or inline, splines, a function from points; a control file rebased onto new data

L078: the core's part of R25 and R26 (APP_WORKFLOWS W3, W4).

### Added

- `background`: our `.drb` file and FullProf's `.bac` (both formats) read and written; linear and natural cubic-spline interpolation; Chebyshev and polynomial fits to points.
- `BackgroundModel::Spline`; `[pattern.background]` `interpolation` and `file` (the points in a file the control file names, or in it); `Control::from_text_in`, `background_file_text`; `schema::read_with`.
- A `.pcr` with Nba = 1 (`<stem>.bac`, b1 and b2) imported; Nba ≤ −6 read as splines, a spline of six points or more exported as Nba = −n.
- `rebase`: a control pointed at new data, what does not fit as diagnostics with suggestions; `drs rebase`.
- `drs background read | convert | eval | fit [--into] | store`.
- Tests `background_points` and `rebase`; `manuals/fullprof/runs/background_runs.py` and `runs/background/` for FullProf on the Mac.

### Known limitations

- Which Nba FullProf reads as splines, and whether its spline is the natural one, waits for `background_runs.py` on the Mac.
- Run with FullProf on the Mac (L079): Nba = 17 is linear (agreement to the `.prf`'s rounding, the test's bound to be set from it), Nba = −17 a natural spline (within 0.079), Nba = −6 not explained; a `.bac` needs a value at every data point; the export holds on all nine runs.
- Nba = 1 with b3 or b4 ≠ 0, and refining b1–b4, are refused.
- The `.drb` format and extension are for Josu to accept.

## 0.18.0 — 2026-10-10 — The core's checks (R24), and LaB6's export numbered without a hole

L076: `drs check` grows into the core's checks of a control file or a `.pcr`, from the syntax to the normal matrix before a run, each problem said where, what, why and what to do, with correlated parameters explained from the parameter dictionary. LaB6's export no longer leaves a hole in FullProf's parameter numbers (L071). L077 records what FullProf and GSAS-II check.

### Added

- `check`: `check_path`, `check_control`, `check_pcr`; `Diagnostic` (severity, level, stable code, span or addresses, what, why, suggestion, replacement, dictionary entries) as text or JSON; `predicted_correlations`.
- `drs check <file.drc|file.pcr> [--json] [--correlation c] [--no-run]`.
- The dictionary's `correlations` (with, why, tip, source) on 12 entries, in `parameter_table::Correlation` and on the Parameters page.
- `lsq::normal_matrix`; `refinable::pcr_parameter_numbers`; `Pcr::declared_parameters`; `Field::phase`.
- Test `check_diagnostics` and `tests/data/broken/` (`make_broken.py`): 18 broken files with their `.pcr` exports for FullProf.

### Changed

- The control file's diagnostics carry a code and, for an unknown name, the replacement.

### Fixed

- The `.pcr` export numbers the coded parameters 1…N without the ratio of `VARY ratio`, and declares N (LaB6, "Hole in the matrix", L071); `control_file` checks it on 62 runs.

### Known limitations

- FullProf has not yet been run on the new LaB6 export nor on the broken files' exports (cloud).
- Predicted and refined correlations differ by up to 0.032 on Co-YBCO, whose cell moves while it refines (frozen).
- The app does not show the checks yet.

## app 0.7.0 — 2026-10-10 — The Editor and versions, with the guardrails of R24; Refine on the control file

The third step of Josu's plan (L067) on core 0.17.0: the control file (.drc) is the document; the Editor checks it as it is typed with the core's checks; versions in .drs/history; Run stops on errors and asks on warnings.

### Added

- The Editor module (CodeMirror 6, vendored): a .drc checked by the core as typed, each problem underlined where it is with a hover (what, why, the suggestion) and a quick fix when mechanical; a Problems panel in the bottom panel; the number of errors on the Editor's icon; folding and an outline from the core's parser; Find and Replace; a .pcr coloured by FullProf's lines; any text file of the project opened.
- Versions (R22): Snapshot… with a note into .drs/history/<file>/NNN_date_time.ext with a sidecar (trigger, note, χ², Rwp, cycles, converged); listed in the sidebar; Diff side by side, Restore (what it replaces is kept), Remove; File ▸ Revert To ▸ Last Saved and Browse All Versions…; a snapshot before every refinement, its sidecar filled with the result.
- Guardrails (R24): Refine ▸ Run first runs the core's checks of the document as it is; errors stop the run and are listed in the Refinement panel with Open in Editor; warnings ask. Refine ▸ Check on demand.
- Refine ▸ Stop (⌘.) ends a run after the cycle in progress (core 0.17.0, request 6).
- Synced both ways: a .drc saved in the Editor is re-read by Refine; a save in Refine refreshes the Editor.

### Changed

- Save, open and import go through drs::control only: a .pcr is imported by writing its .drc beside it (an existing .drc is opened instead); a .drc opens through control::open; Save rewrites the .drc keeping the user's comments (control::update). The 0.5.0 .drc of the model's JSON is gone (L072).
- Every change goes through the core: set_field (the symmetry's ties follow), the edit module, the refinement's parameters as the control file writes them, the set built by the core; the app's interim copies of the ties are deleted (requests 2 and 3).

### Known limitations

- Only a .drc is checked; a .pcr's checks, the refinement's holes (two parameters moving one value are accepted silently today), the physics checks and marks at the inspector's fields wait for the core's check module (REQUESTS.md 13).
- Refine still draws with its own plot (0.8.0).

## app 0.6.0 — 2026-10-10 — The Data module: every format the core reads, the shared plot and table, operations kept apart from the file

The second step of Josu's plan (L067): the Data module on the core's reader registry (core 0.15.0) and data operations (core 0.16.0). No reader and no operation is written in the app; the measured file is never written.

### Added

- Any file the core reads opens in Data: a double-click in the project, Open File…, Data ▸ Open Data File… (⇧⌘D), or a drop on the window; the format recognised by content; a bad file gives a message naming it.
- The shared plot component (ui/plot.js, APP_CONVENTIONS §8): x in 2θ, d or Q, y linear, √ or log, box zoom, wheel and pinch zoom, two-finger pan, double-click to fit, a crosshair readout, a legend whose entries toggle, Copy Plot as Image (⌥⌘C); other open files overlaid; a scan pop-up for files with several scans.
- The shared table component (ui/table.js): sortable and resizable columns, a choice of columns, Find, multiple selection, Copy as tab-separated text; its selection synced with the plot (Select Points, ⇧⌘S).
- Operations, the core's data_ops, kept as an ordered list in .drs/data/<file>.ops.json and applied on every read: mask (⌘⌫) and delete the selected points or a range, smooth, rebin, normalise, add, subtract, merge another file, convert to another wavelength; listed in the sidebar, each removable; Edit ▸ Undo and Redo step through them; Edit ▸ Delete removes the selected points.
- Find Peaks (⇧⌘K) and Automatic Background (⇧⌘B), the core's analyses, drawn over the pattern and listed in the inspector.
- The inspector: the file, the scan (points, masked, axis, range, wavelength), what the file says of itself, the selection, the peaks, the background.

### Changed

- Undo, Redo and Delete in the Edit menu go to the visible module (Refine's document or Data's operations).

### Known limitations

- Export of the processed pattern waits for the core's writer (REQUESTS.md 12).
- The operations' parameters are asked in plain prompts; proper sheets come with the Settings window of a later release. Refine still draws with its own plot until 0.8.0.

## 0.17.0 — 2026-10-10 — Control files end .drc; the .pcr export fixed; a refinement that can be stopped; each phase's curve

L070: the control file's syntax is accepted with the extension `.drc`. The four faults FullProf found in the `.pcr` export (L069) are fixed. The engine stops between cycles when asked, and the calculation gives each phase's curve (the app's requests 6 and 7).

### Added

- `[pattern.fullprof]` in the control file: an imported `.pcr`'s Thmin, Step, Thmax and excluded regions as written (`control::FullProfRecord`), for the export.
- `lsq::refine_until` and `lsq::STOPPED`: a refinement stopped by its caller after a cycle; Ctrl-C in `drs refine`.
- `Calculated::per_phase`: each phase's curve without the background; `drs calc --phases`.
- Tests `refinement_stop` and `phase_curves`; in `control_file`, the exported `.pcr`'s resolution line, range, excluded regions and commands against the original on 62 runs.

### Changed

- Control files end `.drc` (Josu, L070); the content is valid TOML.

### Fixed

- The `.pcr` export: the resolution file's name alone on its line; Thmin, Step, Thmax as the `.pcr` gave them; the excluded regions as written; VARY lines straight after the wavelength line (L069).

### Known limitations

- FullProf has not yet been run on the corrected exports (cloud); `export_runs.py` is to be run again on the Mac.
- Run with FullProf on the Mac (L071): eight runs give FullProf's numbers exactly; LaB6 stops with a hole in the matrix (the export codes the parameters 1–10 and 12–19, counting `VARY ratio` among them).
- Regions changed after the import are written by the rule of L033, which FullProf's handling at the data's ends may not follow.
- NdNi-HTO's nuclear phase differs from FullProf's `.sub` and `.prf` (peaks 0.004–0.015° higher); the run is not imported (magnetic phase) and the difference is not explained.

## app 0.5.0 — 2026-10-10 — The frame: an umbrella of modules over a project folder

The first step of the second design (L060, APP_CONVENTIONS.md): the window of the conventions with two modules registered, Data (empty) and Refine (0.4.0's view moved in as it was). The 0.4.0 shell is deleted.

### Added

- A project is a folder (R20): File ▸ Open… lists it in the sidebar with every file's type from the core's registry of formats (31, core 0.15.0); a double-click opens a file in its module; Open File…, Import ▸ FullProf .pcr…, Open Recent, New Project…, Show in Finder; the layout kept in .drs/project.json.
- The activity bar: Data and Refine (⌘1, ⌘2, badge while refining); light/dark/automatic, Settings, What's New (the release notes, bundled), Manual (the dictionary from the core) at the bottom.
- The dock on dockview-core: sidebar | main | inspector, the bottom panel with the module's panel and the Log; draggable dividers, double-click to collapse, minimum and maximum sizes, sizes remembered per project, View ▸ Hide Sidebar ⌃⌘S / Inspector ⌥⌘I / Bottom Panel ⇧⌘Y / Toolbar ⌥⌘T, Window ▸ Reset Layout; file tabs reordered by dragging.
- The menu bar of APP_CONVENTIONS §3 complete, items enabled only when they can act, later items disabled with their release on the card; one menu per visible module; the Window and Help menus handed to macOS (the window list, the Help search field).
- The module contract (icon, menu, ⌘K commands, manual entry, sidebar section, inspector, bottom panel, views); the dictionary read from the core (parameters::dictionary), no longer from the YAML at build.
- Refine: File ▸ Save writes the model as the core's text beside the .pcr, as <name>.drc (a placeholder extension until the core's control-file syntax is accepted).
- View ▸ Zoom In/Out/Actual Size (the views' text), Appearance ▸ Light/Dark/Automatic; context menus on the project's files; Help ▸ Keyboard Shortcuts.

### Changed

- The Refine commands are the module's (refine_*); the plot resets with a double-click or Refine ▸ Fit Plot (⇧⌘0), since ⌘0 is Actual Size.
- The window's toolbar sits in the title bar (titleBarStyle Overlay).

### Known limitations

- Data is an empty view until app 0.6.0; the Editor and versions until 0.7.0; Refine's interim copies of the core's ties until 0.8.0.
- Panels do not yet float or detach, tabs do not yet split or drag out, no native window tabbing (0.8.0).
- Settings is a dialog, not the standard window with panes; manuals are simple viewers; no file associations, no watcher for files changed outside the app, no auto-save.

## 0.16.0 — 2026-10-10 — Our own control file, and the data operations of the Data module

L066: the control file of diff-ref-struct is a strict subset of TOML, for Josu to accept. It is read and written by the core with exact positions, comments kept, and a lossless round trip on every imported run; it is exported to `.pcr` where FullProf can express the model. The Data module's operations are non-destructive and checked against SciPy, pybaselines and GSAS-II's code.

### Added

- `control::syntax`: the parser of the TOML subset: every key and value with its line, column and span; comments and order kept; plain messages; recovery after an error; an outline of the tables for folding.
- `control::schema`: the control file under the parameter dictionary's names with units: pattern (data by relative path, reader id, instrument, corrections, excluded regions, background, resolution, operations), phases, atoms, restraints, refinement (parameters as their terms, equalities, limits, linear restraints, cycles, convergence), the convention switches; format line `diff-ref-struct 1`; canonical writer; unknown keys with a suggestion.
- `control::update`: a file written again after a change keeps its comments, order and notation; only changed values are rewritten.
- `control::pcr_export`: a `.pcr` with its `.irf`, `CODFILn.hkl` and processed `.xys`, the import's translations inverted exactly; notes for what FullProf cannot express.
- `refinable::Refinement`: a refinement as data (`build`, `from_pcr`), reproducing `ParameterSet::from_pcr`.
- `data_ops`: mask, delete, smoothing, rebinning, normalising, adding, subtracting, merging, conversion between 2θ, d and Q, change of wavelength; peak search, GSAS-II's peak candidates, automatic background (arPLS).
- `formats::project`: the kinds of a project's other files.
- `SpaceGroup::fullprof_symbol`: the symbol FullProf reads for a group.
- `drs convert <file.pcr> --to <file.toml>` and back; `drs check <file.toml>`; `drs calc` and `drs refine` read control files; `drs data info | op | apply | peaks | background`.
- Tests `control_file` (62 runs to the bit, both ways) and `data_operations` (four real files against SciPy, pybaselines and GSAS-II; identities); script `make_data_ops_references.py`.

### Known limitations

- The syntax is for Josu to accept.
- FullProf itself was not run on the exported `.pcr` files (cloud); they were read back by the core.
- Run with FullProf on the Mac (L069), the export is refused for six of the nine runs and differs for one: the resolution file's line, the excluded regions, the range and the place of `VARY ratio`.
- GSAS-II's rejection of peak picks closer than two FWHM is not reproduced (it needs the instrument's profile).
- Conversions from time of flight or energy need the instrument's constants and are not done.

## 0.15.0 — 2026-10-10 — The reader registry: every data format of FullProf's list, a general table, and GSAS-II's

L062: one registry of 31 data formats for `drs` and the app. Each format has a readable name, and files are recognised by content and then by extension. The readers are checked against FullProf's and GSAS-II's reading of the same files.

### Added

- `formats`: the registry (`FORMATS`, `identify`, `read`, `read_path`, `by_id`, `by_ins`); `Dataset` and `Scan` with axis (2θ, TOF, Q, d, E), wavelength and `Info` (title, instrument, date, sample, temperature, monitor, counting time, step).
- FullProf's list: Ins 0–12 and Ins 14 (ISIS TOF banks, from the 2003 News), readable ids (`free_steps`, `ill_d1a_d2b_old`, `ill_d1b_old`, `ill_d1b_d20`, `brookhaven`, `two_axis`, `ill_d2b`, `ill_d4_d20l`, `psi_dmc`, `socabim_uxd`, `xys`, `variable_time`, `gsas`, `isis_tof_banks`); Ins −4 (DBWS) recognised.
- `table`: the general table reader (delimiter, decimal mark, header lines, columns and quantity recognised or given as `TableOptions`).
- GSAS-II's formats: Bruker RAW 2, 3, 4 (1 recognised) and BRML, PANalytical XRDML, Rigaku RAS, RASX and TXT, xye, qye, Fit2D/pyFAI chi, powder CIF; NeXus/HDF5 and MIDAS zarr recognised (L063).
- Zip, DEFLATE and a small XML reader in the core: no new dependency.
- `drs data formats`, `drs data identify <file>`, `drs data read <file> [--as <id>] [--scan n] [--out]` with the table options.
- Test `data_formats` (FullProf's reads, GSAS-II's reads, FullProf's powder CIFs, recognition, formats without example); scripts `make_fullprof_reads.py`, `make_gsas2_reads.py`, `make_synthetic_formats.py`; `manuals/gsas2/` with the GSAS-II tutorials' data files.

### Changed

- `data::read(text, Ins)` reads through the registry; `data.rs` keeps the pattern type.
- GSAS files: every bank, NCHAN points per STD and ESD bank, TOF in µs, no weight (σ = ∞) for y ≤ 0, as GSAS-II.
- The FullProf readers check the structure of real files (integer I fields, written counters, ILL end marks), which `identify` uses.

### Known limitations

- Five FullProf runs' Σ w·y² differ from FullProf's by up to 1.8·10⁻⁴; frozen: the point at an excluded region's end, and one unexplained (citric acid).
- Ins 9 (UXD) and Ins 14 are read as believed: no layout in the manual, and no FullProf output for them.
- 20 formats have no real example file; 10 of them are checked on synthetic files against GSAS-II.
- NeXus/HDF5 and MIDAS zarr are recognised, not read.

## 0.14.0 — 2026-10-10 — A smooth peak window, FullProf's angle-restraint weight, and the model editable from the app

L058: the default model's peaks reach zero smoothly; FullProf's restrained refinements reproduced with its angle restraints at π/180 of their stated weight; and the app's requests 1–5: the model enumerated, edited, written and read as text, the parameter set built from addresses, the whole dictionary in the core with display names.

### Added

- `profile::taper` and `Convention::window_edge`: the default window keeps each peak to ±Wdt·H and tapers it to zero by ±1.25·Wdt·H (FullProf's cut stays in the compatibility layer, kept on import).
- `Convention::angle_restraints`: FullProf's angle restraints at π/180 of their stated weight (inferred from its runs; the default and the import keep the stated weight).
- `edit`: `Model::fields` (`FieldInfo`: address, id, value, limits, fixed_by), `Model::all_fields`, `Model::fixed_by`, `Model::set_field`; `add_phase`, `add_phase_from_cif`, `remove_phase`, `add_atom`, `remove_atom`, `set_background`, `set_background_terms`, `add_excluded`, `remove_excluded`, `set_shape`, `set_asymmetry`, `set_preferred_orientation_function`, `set_preferred_orientation` (app request 3).
- `refinable::field_of` and `ParameterSet::from_addresses` (free addresses, additive equalities, user limits); `ParameterSet::user_limits` and `limits_of` (app request 2).
- `text`: `Model::to_text` and `Model::from_text`, behind the `serde` feature (on by default); `Serialize`/`Deserialize` on the model's types, the space groups as Hall symbols and operators (app request 1).
- `parameter_table` carries every field of every entry and `GROUPS`; `parameters::dictionary()`, `parameters::groups()` (app request 4); a `name` per entry in `parameters.yml` (app request 5).
- Tests `model_editing` and `model_text`; `fullprof_restraint_runs` runs by default; two default-model passes in `engine_derivatives`; the tapered window in `esd_profile`.

### Changed

- The default model's peak window (the taper); `default_convention`'s table re-measured.
- `from_pcr`'s symmetry and columns are shared with `from_addresses`; the engine takes a set's user limits.
- `make_fullprof_profiles.py` gives the `CODE_conv.prf` files their own reason.

### Known limitations

- The weight of FullProf's angle restraints is inferred (π/180 and (π/180)² fit equally); its code is not public.
- Sr-ox-test's Fourier background and LaMnO3-3t2's peaks at 21–41° differ from FullProf's (χ² +1.5 % and +4.9 % at its values); frozen, not explained.
- `Model::set` stays the engine's unchecked setter; the app's checked `set` is `Model::set_field`.
- The app's requests 6 (stopping a refinement), 7 (per-phase curves) and 8 (data without a .pcr) are open.

## app 0.2.0 — 2026-10-10 — Refine, with the live plot

The app refines a `.pcr` with the core's engine and redraws the plot after every cycle (L057).

### Added

- **Refine**: the `.pcr`'s coded parameters refined off the window's thread, Ycal and the difference redrawn after every cycle.
- A refinement panel: the cycles (Rwp, χ², largest shift/esd), the statistics, R_Bragg and R_F per phase, the parameters as value(esd) with limits marked, the correlations above 0.8.

### Known limitations

- The refined values are not saved: the `.pcr` is not changed, and there is no document yet.
- Parameters are listed by their dictionary ids; labels and help from the dictionary come with the app's own work (branch `app`).

## 0.13.0 — 2026-10-10 — Restraints, Le Bail fits and quantitative phase analysis

Phase 2, step 5 with restraints (L055): distance, angle and linear restraints in the engine, Le Bail extraction of intensities, weight fractions with esd's, and a test that checks the esd's against the χ² surface itself.

### Added

- `restraints`: distance and angle restraints through FullProf's ITnum operators (manual p. 123), values and exact derivatives; `Model::restraints`; read from a `.pcr`'s Dis and Ang lines.
- `refinable::LinearRestraint`: Σ c·p ≈ v ± σ on named parameters; `ParameterSet::linear`, read from NLI lines.
- The engine adds restraints as observations weighted by χ²_ν/σ², as FullProf's; `Refinement::restraints` and `lsq::restraint_states`.
- `lebail`: `Model::le_bail`, one extraction by the Rietveld formula (eq. 3.15) for every Jbt = 2 phase; the engine extracts at every cycle and runs all its cycles; Irf = 0 without an `.hkl` generates the reflections.
- `qpa`: `weight_fractions` and `weight_fractions_esd` (Hill & Howard, manual p. 31), esd's from the covariance; `Phase::atz`, `Phase::brindley`; dictionary entries `qpa_atz`, `brindley_coefficient`.
- `Refinement::covariance`.
- `drs refine` prints the restraints (target, value, Δ/σ) and the weight fractions, and refines Le Bail phases.
- Tests `fullprof_restraints` (FullProf's `.tpcr` distances and angles, derivatives, three restrained refinements, a linear restraint), `fullprof_le_bail` (LBF_lmopr, MgAl2O4s_LeBail), `esd_profile` (the esd's against the profile of S), and the QPA part of `fullprof_refinements`.
- `manuals/fullprof/runs/restraint_lebail_runs.py`: FullProf references for restraints and Le Bail fits, to run on the Mac (test `fullprof_restraint_runs`, ignored until then).

### Changed

- Profile matching with Irf = 0 no longer needs an `.hkl` file.

### Known limitations

- No FullProf output of a refinement with restraints in the material yet: the restrained refinements are checked for convergence and their restraints, not against FullProf.
- LaMnO3's FullProf `.tpcr` distances differ from the core's by up to 1.6·10⁻⁴ Å, and D20's angles by 0.0135°: not explained, frozen.
- Le Bail: LBF_lmopr's χ² at FullProf's values is 0.43 % below FullProf's, and MgAl2O4s_LeBail's global χ² 0.20 % above; frozen.
- Weight-fraction esd's are smaller than FullProf's, which treats the sum as independent of each scale.
- Brindley coefficients other than 1 are not checked (none in the material).
- Rigid bodies, TLS and moment restraints are not read.

## 0.12.0 — 2026-10-09 — The least-squares engine

Phase 2, step 4 (L053): the parameters of a `.pcr`'s codewords, named by the dictionary, with constraints and limits, refined by Levenberg–Marquardt with analytic derivatives. Starting from the tutorials' files, it reproduces five of FullProf's converged refinements: χ² to the tolerance of `fullprof_profiles`, every parameter within FullProf's esd. The differences on the other four are explained where possible and frozen.

### Added

- `dual`: forward-mode dual numbers; the profile functions (`gaussian`, `lorentzian`, `pseudo_voigt`, `tch`, `berar_baldinozzi`, `Fcj`) are generic over them.
- `jacobian`: the refinable fields (`Field`), `Model::get`/`set`, and `Model::derivatives` (Ycal and its exact derivatives at every point); `Structure::structure_factor_derivatives` and `Scatterer::ds`.
- `refinable`: `ParameterSet::from_pcr`, the codewords as free parameters named by the dictionary, shared codewords as constraints with their multipliers, the crystal system's cell constraints, the import's translations; `refinable::limits` (P39).
- `lsq`: `refine` (Levenberg–Marquardt, the normal matrix scaled and solved by eigen-decomposition, singular directions left out, limits as an active set, convergence on shift/esd < Eps), each cycle an event (`Cycle`); `statistics` (Rp, Rwp, Rexp, χ², Durbin–Watson, R_Bragg and R_F per phase); covariance, esd's and correlations.
- `drs refine <file.pcr> [--cycles N] [--data <file>]`: R-factors per cycle, the parameters with their esd's, limits reached, R_Bragg and R_F, correlations above 0.8.
- The parameter dictionary has `limits` for every refinable entry, on the Parameters page and in `parameter_table`.
- Test `engine_derivatives`: every kind of field of ten FullProf runs against finite differences (to 2·10⁻⁵); the dual path against `calculate` (10⁻¹⁰).
- Test `fullprof_refinements`: the nine converged runs of L052 refined from their starting files and compared with FullProf's: χ², N − P, each parameter, each esd.
- Test `parameter_dictionary::engine_names_and_limits`: the engine's names and limits against the dictionary.

### Changed

- `Model::calculate` is built from `select_peaks`, `excluded_points` and `peak_profile`, which the derivatives reuse.
- The `.pcr` reader keeps the codewords of the atoms, betas, scale, widths, cell, preferred orientation, zero and background, and Eps.

### Known limitations

- esd's are scaled by χ² over all points used; FullProf's by χ² over the points with Bragg contributions, so FullProf's are larger by √(χ²_Bragg/χ²_all).
- With the TCH profile, FullProf prints half the esd of the cell parameters and I_G that the engine finds; why is not known (L053).
- FullProf's peak window depends on how many cycles it has run (rutile + anatase, LaB6): its later cycles are not reproduced in general; LaB6's differences are frozen.
- Co-YBCO: FullProf stops 2.7 % of χ² above the minimum the engine reaches.
- No restraints (P40), Le Bail, quantitative phase analysis, magnetic phases, wavelength or anisotropic size and strain refinement yet.
- A value read outside its limits is not moved onto them; it may only move towards them.

## app 0.1.1 — 2026-10-09 — A .pcr opened through the translated import

The app opens a `.pcr` through the core's import into the default model (L050, L051).

### Changed

- A `.pcr`'s values are translated into the default model on opening, instead of being calculated under FullProf's rules; the header says which of FullProf's rules are kept (the peak window, the reflections past the data).

## 0.11.0 — 2026-10-09 — The papers decide the defaults; a .pcr is translated into the default model

Phase 2 (L050): the six papers in `manuals/papers/` settle the default rules; `Convention` goes rule by rule; and importing a `.pcr` translates its values exactly where FullProf's rules differ, keeping FullProf's only where no exact translation exists.

### Added

- `Model::from_pcr_fullprof` (the `.pcr`'s values under FullProf's rules) and `Model::translate_from_fullprof`; `profile::fullprof_axial`; `Model::fullprof_excluded`.
- Test `default_convention::pcr_import`: 47 FullProf runs imported, χ² identical to FullProf's rules; the effect of the two kept rules frozen per run.

### Changed

- `Convention` is a struct with one switch per rule, `Convention::DEFAULT` and `Convention::FULLPROF`, and `fullprof_rules()`.
- `Model::from_pcr` translates: X-ray scale × 2, Bérar–Baldinozzi P1 + 3·P2 and P3 + 3·P4, the axial pair, G1 → −G1 for the exponential, excluded regions on the data's grid; FullProf's peak window and reflection range kept.
- Default Bérar–Baldinozzi: the paper's sign (positive on the low-angle side, as FullProf; the manual's eq. 3.29 prints the opposite) and F_b = (8u³ − 12u)e^{−u²} (eq. 13).
- Default FCJ takes S/L, H/L as given (`Fcj::new`); FullProf's small-D_L reading is a convention.
- Default exponential preferred orientation is Dollase's exp(−G1 α²).
- `drs calc`, `drs simulate`, `drs mic` import by translation; `--fullprof` gives the values as written under FullProf's rules. The app opens a .pcr through the translated import.

### Known limitations

- The peak window and the reflections past the data keep FullProf's rules on import: no exact translation (up to 5.5% of χ² on YVO3 if switched).
- NdNi-LTOa, NdNi-LTOb and Pd3MnD are not imported: anisotropic size and strain are version 2.
- The app's change could not be built in the cloud session (no GTK/WebKit).

## app 0.1.0 — 2026-10-09 — The walking skeleton of the Mac app

The first build of the app (L049): a Tauri 2 window around the core, linked as a library.

### Added

- **Open .pcr…**: a FullProf control file and its data are read through the core and the pattern calculated, under the FullProf convention (L047).
- The Rietveld plot: Yobs, Ycal, the difference, a row of reflection markers per phase, excluded regions shaded; zoom, pan and a readout under the cursor; Rp, Rwp and χ².

### Known limitations

- No refinement yet: the core's least-squares engine is the next step of phase 2.
- No menus, documents, editor or live plot yet (phase 3).
- The app is not signed: macOS asks for confirmation the first time it is opened.

## 0.10.0 — 2026-10-09 — The default model from FullProf's manual, FullProf's rules behind a switch, and the parameter dictionary

Phase 2 (L046): where FullProf computes something other than its manual or the microstructure paper prints, the core's default now follows the documents and `Convention::FullProf` keeps FullProf's rules, which the FullProf tests use. The parameter dictionary names every parameter the core reads.

### Added

- `powder::Convention` (Default, FullProf) and `--fullprof` on `drs calc`, `drs simulate` and `drs mic`.
- The parameter dictionary in the core: `parameter_table` (generated from `ecosystem/data/parameters.yml`, 45 entries) and `parameters::FIELDS`, mapping every field of `Model`, `Phase` and `Atom`.
- Test `parameter_dictionary`: every model field mapped or named as no parameter; the table matches the YAML.
- Test `default_convention`: each default rule checked, and the default's distance from FullProf's convention on ten runs frozen (largest 4.5·10⁻² of a peak, Y2O3 at 29°).

### Changed

- Default X-ray Lorentz–polarisation is eq. 3.26 with K = ½, half FullProf's: a default model's scale is twice FullProf's.
- Default Bérar–Baldinozzi asymmetry is eq. 3.29 as printed (z, F_b = 2(2z² − 3)F_a); `profile::berar_baldinozzi` takes a `fullprof` flag.
- Default peak window is ±Wdt·FWHM, with no extra low point; default reflections are every one whose window reaches the data; default excluded points are those within the regions.
- Default microstructure takes β* with each reflection's own λ, and gives no strain for U < 0.

### Known limitations

- The original papers (Bérar & Baldinozzi 1993 above all) could not be read: the defaults follow FullProf's manual and are provisional until they are (L046).
- The FullProf differences frozen in 0.9.0 (YVO3 at Wdt = 40, the rutile-anatase and MgAl2O4 window edges, 15 points of 12 experiments) are unchanged.

## 0.9.0 — 2026-10-09 — The peak window's edge points, the parametric resolution file and March–Dollase as FullProf does them

Phase 2, step 3 continued (L042): the 35 FullProf window experiments and a preferred-orientation run in the profile test. The core now has CrysFML's low edge, the wavelengths, monochromator and axial divergence a parametric resolution file gives, Bérar–Baldinozzi on top of FCJ, and March–Dollase averaged over equivalent reflections.

### Added

- `fullprof_profiles` compares 53 files: the 35 window experiments and Yb2ZrWO8-PO-March. `fullprof_corrections` covers 11 runs, 1707 reflections, preferred orientation included.
- Files named in a `.pcr` are also looked for in the grandparent folder and its folders.

### Changed

- The window of a TCH peak whose phase has X ≠ 0 starts one point below 2θ − Wdt·H, as CrysFML's `Locate` does. LaB6's edge points fall from 50 to 2.
- A parametric resolution file's WAVE and CTHM replace the `.pcr`'s wavelengths, Ratio and Cthm. Its ASYM S_L, D_L is added to the phase's instead of only filling a phase that has none.
- With both FCJ and Bérar–Baldinozzi asymmetry, the Bérar–Baldinozzi factor multiplies the FCJ profile.
- The preferred-orientation factor is averaged over the reflections equivalent to h.

### Known limitations

- YVO3 at Wdt = 40 still differs (rms 2.34·10⁻²). Rutile-anatase, MgAl2O4 and 15 points of 12 window experiments differ at window edges. All are frozen in `fullprof_profiles` and listed for Josu (L042).
- The exponential preferred-orientation function has no FullProf reference.

## 0.8.0 — 2026-10-09 — Isotropic size and strain, and the corrections, against FullProf's runs

Phase 2, step 3 (L039): the apparent size and strain FullProf writes in its `.mic` files, widths at the Bragg angle, Lp, absorption and the weights of Iwg checked reflection by reflection on FullProf's runs, and four runs of tutorial 6 in the profile test.

### Added

- `microstructure` and `drs mic`: integral breadths, apparent size and maximum strain per reflection, as FullProf's `.mic` gives them.
- `Weighting` (Iwg): maximum likelihood and unit weights in the agreement factors.
- Test `fullprof_microstructure`: the `.out` widths of 9 runs (2710 reflections) and every `.mic` row of YVO3, Yb2ZrWO8 and NdNi-HTO (688).
- Test `fullprof_corrections`: m·Lp·A·P against FullProf's CORR, 1691 reflections of 10 runs, absorption included.
- `fullprof_profiles` compares Yb2ZrWO8, LaB6 and YVO3 (17 files); `fullprof_spacegroups` covers 37 `.out` files.

### Changed

- Widths are evaluated at the Bragg angle, without Zero and the displacement terms (FullProf's `.out` widths to 2·10⁻⁵).
- Surface micro-absorption is refused with a message instead of being ignored.

### Known limitations

- The peak window when X ≠ 0 and the YVO3 profile differ from FullProf by rules not found; frozen in `fullprof_profiles` and listed for Josu (L039).
- Preferred orientation is not checked: no FullProf run uses it. Magnetic phases and anisotropic size and strain are not calculated.

## 0.7.0 — 2026-10-09 — FullProf's own structure factors, and the pattern of its reference runs

Phase 1 closed against FullProf's own structure factors, and the calculated pattern brought to FullProf's on the five reference runs (L035). The pattern changes where FullProf and its manual differ, so this is a minor release.

### Added

- Test `fullprof_reference_sf`: F, real and imaginary, of the 11 phases of the five FullProf runs (1992 reflections), within FullProf's print and single precision.
- `Phase::reflection_limit` (largest 1/d of a phase's reflections) and `Model::window_width`.
- `fullprof_profiles` compares the five runs (`exact`), counting deviations beyond FullProf's print, and their χ² against FullProf's.

### Changed

- Bérar–Baldinozzi asymmetry as FullProf computes it: z of the opposite sign, F_b = (4z² − 3)·F_a (Y2O3 from rms 1.7·10⁻² to 7·10⁻⁵).
- The peak window: ±Wdt widths about 2θ_h, also for the FCJ profile, the width for TCH leaving out the phase's X·tan θ.
- Reflections are generated from the first point up, no longer 5° below it.
- A point is excluded when a step to or from it reaches into an excluded region.

### Known limitations

- FullProf's reflection limit past Thmax and the last 0.06–2.7% of its window when X ≠ 0 follow rules not found: the tests take the first from FullProf's `.out` and report the points of the second apart (L035).
- Cocus-firstpeak gives χ² 12.2 against its header's 3.63, unexplained (no FullProf run of it).

## 0.6.0 — 2026-10-09 — The calculated powder pattern, as FullProf computes it

Phase 2, steps 1 and 2 (L033): the core reads FullProf `.pcr` files and the data formats of the tutorials, and calculates the pattern from them, point by point as FullProf does.

### Added

- `pcr`: single-pattern FullProf `.pcr` files (job, flags, wavelengths, background, excluded regions, phases, atoms with betas and symmetry modes, profile, asymmetry, preferred orientation, absorption), codewords kept as written; resolution files. All 151 tutorial `.pcr` files read.
- `data`: FullProf data formats Ins 0, 1, 3, 5, 6, 8, 10 (ported from the manual's Fortran) and GSAS Ins 12, with σ and variable steps.
- `profile`: pseudo-Voigt, TCH pseudo-Voigt, Bérar–Baldinozzi asymmetry, and Finger–Cox–Jephcoat axial divergence ported from CrysFML's `Prof_Val` (commit 8ed8670) with its Gauss–Legendre tables (`fcj_table`).
- `powder`: the calculated pattern: positions with Zero, Sycos, Sysin; Caglioti widths and the resolution file; Lorentz–polarisation by radiation; overall B; preferred orientation; absorption; backgrounds (points, polynomial, Chebyshev, Fourier); several phases; excluded regions; agreement factors; Poisson simulation.
- `drs calc <file.pcr> [--data f] [--out f]` (2θ, Ycal, Backg, Yobs, σ) and `drs simulate` (XYDATA with Poisson noise).
- Test `fullprof_profiles`: 9 `.prf` files rebuilt from their `.pcr`, Ycal to rms 1–4·10⁻⁴ (7 files), 2.7·10⁻³ (Le Bail), Backg to 5·10⁻⁴, positions to 5·10⁻⁴°; the observed data of 80 `.prf` files (301429 points) read as FullProf read them; the 76 files left out listed with their reason.

### Known limitations

- Size and strain broadening, TLS, Cry ≠ 0, magnetic phases, multi-pattern and time-of-flight files are refused with a message.
- No `.xrdml`, `.raw` or `.brml` reader: the material holds no example. Ins 5, 8 and 12 are not checked against FullProf.
- LaMnO3ref.prf was not drawn from the state of its `.pcr` (rms 1.6·10⁻²); Cocus-firstpeak, Egypt and Sr-oxalate give a χ² well above their `.pcr` header, unexplained (L033).

## 0.5.0 — 2026-10-09 — Structure factors with FullProf's conventions

The last item of phase 1 (L030): structure factors for X-rays, with anomalous terms, and for neutrons, with FullProf's occupancies and displacement parameters. They are checked against every F² FullProf printed in the material (12026 reflections), against gemmi for the four reference structures, and by symmetry properties in all 530 settings.

### Added

- Structure factors summed over all the operators of the group, with FullProf's Occ (chemical occupancy × m/M), Biso and betas; conversions between Occ and chemical occupancy, and from U to betas and B. |F|² is FullProf's printed F²: f′ and f″ included, the overall B left out.
- Scattering tables generated from CrysFML (FullProf's library): 214 X-ray form factors, ions included; f′ and f″ for 98 elements at the Kα1 of Cr, Fe, Cu, Mo and Ag, taken at the anode nearest to λ as FullProf takes them; 108 neutron scattering lengths.
- A small CIF reader: cell; symmetry operators, Hall or H-M symbol; atom sites; U or B; anisotropic U_ij or B_ij; ionic type symbols.
- `drs sf <file.cif | symbol + cell + --atom …> --dmin <d> [--xray <λ or anode> | --neutron]`: h k l, multiplicity, d, |F|, F².
- `SpaceGroup::from_fullprof_symbol`: a symbol read as FullProf reads it, e.g. "F d -3 m" as origin choice 2.
- Test `fullprof_structure_factors`: PPH3, citric acid (X-rays) and LaMnO3 (neutrons, special positions). Of 12026 printed F², all agree within the precision of the printed atoms, at most 0.37 of that bound. Controls without f′, f″ or with chemical occupancies fail.
- Test `gemmi_structure_factors`: 10 phases of the four reference structures (Y2O3, MgAl2O4, Co-doped YBCO, the QPA phases), 1435 reflections, X-rays and neutrons. With gemmi's coefficients the agreement is 6·10⁻⁸; table by table, only O²⁻ differs.
- Test `structure_factor_properties` (530 settings, 182320 reflections): absences, Friedel and Bijvoet pairs, real F with −1 at the origin, equivalent reflections, and the occupancy convention against the structure written out in P 1. Test `cif_structures`.

### Known limitations

- The material holds no FullProf structure factors for the four reference structures. They were checked against gemmi; FullProf's own (Hkl = −3) can be produced by running FullProf on their tutorial files.
- The O²⁻ X-ray form factor of CrysFML (two Gaussians) differs from gemmi's (four); the sources of each were not checked.
- f′ and f″ come from five tabulated anodes only, as in FullProf; values at other energies (synchrotron) cannot be given yet.
- The CIF reader takes the first block with atom sites, skips text fields and ignores standard uncertainties.

## 0.4.0 — 2026-10-09 — Reflections as FullProf lists them

The fourth item of phase 1 (L028): the unique reflections of a phase between two d-spacings or in a 2θ range, with systematic absences derived from the operators, FullProf's choice of representative, multiplicities, d and 2θ, checked against 14337 reflections printed by FullProf and against direct structure factors in all 530 settings.

### Added

- Systematic absences from the operators: h absent when an operator {R|t} has h·R = h and h·t not an integer (centring, screw axes, glide planes); reflections forbidden only by special positions are kept.
- Unique reflections chosen as FullProf (CrysFML Hkl_Uni) chooses them: index box and region of the Laue class, by space-group number and unique axis; multiplicities with Friedel mates; sorted by decreasing d.
- `drs hkl <symbol> <a b c α β γ> --dmin <d>` (or `--lambda <λ> --2theta <min> <max>`): h k l, multiplicity, d, 2θ.
- Test `fullprof_reflections`: 9 phases and 14337 reflections of FullProf's `.out` files, row for row: h k l, multiplicity, d and 2θ to the printed decimals.
- Test `reflection_properties`: in all 530 settings, absences against a directly computed structure factor (197980 absent, 506920 allowed), multiplicities against a count in a sphere, one FullProf representative per set in every standard setting.

### Known limitations

- FullProf's choice of representative is checked against FullProf output only for the Laue classes mmm and 2/m; for the others it rests on the port of CrysFML's rules and the property tests.
- In 42 non-standard settings (monoclinic with unique axis a, R groups on rhombohedral axes) FullProf's rule does not single out one reflection per set; a stated fallback is used and such reflections are marked '*'.
- One pair in FullProf's PPH3 profile-matching list is in an order that contradicts FullProf's own printed d; not explained.
- No structure factors yet: the last item of phase 1.

## 0.3.0 — 2026-10-09 — Sites, Wyckoff positions and special-position constraints

The third item of phase 1 (L026): for any position in a space group, its site symmetry, multiplicity, orbit and Wyckoff position, and the constraints its symmetry puts on the coordinates and the anisotropic displacement parameters, checked against spglib for every Wyckoff position of the 230 standard settings and against FullProf.

### Added

- Site-symmetry group of a position (exact operators), site point group, multiplicity and orbit in the cell.
- Wyckoff letter and oriented site-symmetry symbol for the 530 settings, from spglib's Wyckoff table.
- Constraints on x, y, z and on the betas (and so on the U_ij of a CIF) as exact linear relations over free parameters, the earliest variables free as in FullProf (x, 1/4, z; x, 2x, z; β22 = β11, β12 = β11/2).
- FullProf's space-group symbols (CrysFML's table): older names such as "B m a b" or "F m 3 m" are read, and a group given by any Hall symbol is matched to its table setting.
- `drs site <symbol> <x> <y> <z>`: site operators, multiplicity, Wyckoff position, equivalent positions, constraints.
- Test `spglib_wyckoff`: 3462 positions, every Wyckoff position of the 230 standard settings; multiplicity, letter, site symbol, site point group, free coordinates, free betas (by the character formula), and multiplicity × site order = group order.
- Test `fullprof_sites`: FullProf's multiplicities for the 96 atoms of its `.out` files; the betas and codes of the 69 anisotropic atoms of its `.pcr` files obey the core's constraints.
- Tests `crysfml_symbols` (612 FullProf symbols) and `u_ij_follow_beta_relations` (530 settings).

### Known limitations

- Wyckoff letters of non-standard settings come from the table and were not checked against an outside analysis; the 230 standard settings were.
- A position counts as special within 0.001 of the symmetry element (`--tol` to change it); FullProf's own tolerance was not checked.
- 29 cells of CrysFML outside the 530 settings (A 1, I -1, F 2/m, C 4 2 21, F 4/m m m and others) work as groups but have no Wyckoff letters.
- No reflections or structure factors yet: the next items of phase 1.

## 0.2.1 — 2026-10-09 — Command-line fixes after the review of 0.2.0

Found when 0.2.0 was reviewed and merged (L025). No change to what the core computes; all tests pass as before.

### Fixed

- `drs sg --help` and `drs cell --help` (or `-h`) print the usage instead of reading `--help` as a symbol or a number.
- `drs` no longer panics when its output goes to a pipe that closes early (`drs sg --list | head`).
- The eight clippy warnings of 0.2.0 (style only): the settings table is a `static`, and three expressions are simplified.

## 0.2.0 — 2026-10-09 — Space groups from Hall symbols

The second item of phase 1 (L024): space groups generated from Hall symbols, as general operators, with the 530 settings of the 230 types, every one checked against spglib and the FullProf ones against FullProf.

### Added

- Symmetry operators {R|t} with exact translations (twelfths): product, kind of rotation, reading and writing "x,y,z" triplets.
- A Hall-symbol parser: generators from the lattice symbol, matrix symbols and change of origin; the group closed under multiplication; centring vectors, inversion (and whether at the origin), point group, Laue class, crystal system.
- The table of the 530 settings (Hall number, number, choice, Hall and Hermann–Mauguin symbols) from spglib 2.7.0's table of settings; lookup by number, Hermann–Mauguin symbol (with :choice) or Hall symbol.
- `drs sg <symbol>` lists a group's operators (one per rotation, or every operator with --all); `drs sg --list` lists the 530 settings.
- Test `spglib_530`: all 530 settings, 7388 operators, point groups, crystal systems, centring and centrosymmetry against spglib, no disagreement.
- Test `fullprof_spacegroups`: all 13 FullProf `.out` files (P n m a, P 21/c, P 21/a, P 2/m): symbols, setting, system, Laue class, point group, multiplicity, centrosymmetry and operators, no disagreement.

### Known limitations

- Change-of-origin vectors are read; general change-of-basis matrices in Hall symbols are not (none of the 530 settings uses one).
- The International Tables, Vol. A, were not compared directly; spglib and FullProf were.
- No Wyckoff positions, site symmetry or reflection conditions yet: the next items of phase 1.
- `drs` stops with a "broken pipe" message when its output is cut short (e.g. piped to `head`).

## 0.1.0 — 2026-10-09 — The Rust core and the unit cell

The first build of the core (L023): the crate `drs` with its command line, and the unit cell, checked against FullProf's own output.

### Added

- The `drs` crate in `core/`: a library for the app and the `drs` command line (`drs cell …`, `drs --version`).
- The unit cell: direct and reciprocal parameters, metric tensors, volume, d(hkl), 2θ for a wavelength, FullProf's Cartesian frame and the Busing–Levy B matrix.
- Test `fullprof_metric`: the metric block of all 13 FullProf `.out` files in `manuals/fullprof/` recomputed from the printed cell, with no disagreement.

### Known limitations

- No space groups, reflections or structure factors yet: the next items of phase 1.
- Cell parameters carry no standard uncertainties yet.
