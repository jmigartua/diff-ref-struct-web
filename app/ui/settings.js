/* diff-ref-struct · the Settings window's model (APP_CONVENTIONS §4, app 0.30.0): its panes, every setting with its
   default and a one-line explanation, which a project may override, and how a value is found (the project's, else the
   app's, else the default). Rows of kind "info" say what the app does where §4 names a choice it does not offer yet, so
   that nothing pretends to be a setting. No module imports: the node tests read it as it is. */

/** Each pane: id, title, glyph, rows. A row: key (in the settings), label, help (one line), kind (bool, select, number,
    path, info, button, formats), default, options ([value, words]), min/max/step/unit for a number, overridable (a
    project may set its own), action (a button's). */
export const PANES = [
  { id: "general", title: "General", glyph: "⚙", rows: [
    { key: "launchOpens", label: "At launch, open", kind: "select", default: "session", options: [["session", "the windows of the last session"], ["recent", "the most recent project or data file"], ["welcome", "the welcome page"]], help: "What the app shows when it starts with nothing asked of it: the windows open when it last quit, each where it was with its tabs; the newest of File ▸ Open Recent; or the welcome page. A file opened from the Finder opens as asked whatever is chosen." },
    { key: "recentCount", label: "Recent projects kept", kind: "number", default: 10, min: 1, max: 30, step: 1, help: "How many projects File ▸ Open Recent lists, the newest first." },
    { key: "projectsDir", label: "Projects folder", kind: "path", default: "", help: "Where New Project and Work on a Copy put a project; empty is ~/Documents/diff-ref-struct." },
    { key: "_shortcuts", label: "Keyboard shortcuts", kind: "button", action: "shortcuts", button: "Shortcuts…", help: "Change a command's shortcut, or give it none; each is checked against the system's and the app's." },
    { key: "_updates", label: "Updates", kind: "info", text: "None, decided (L337): each version is installed by its builders; What's New opens by itself after each new version.", help: "" },
  ] },
  { id: "appearance", title: "Appearance", glyph: "◐", rows: [
    { key: "theme", label: "Appearance", kind: "select", default: "system", options: [["system", "Automatic (the system's)"], ["light", "Light"], ["dark", "Dark"]], help: "Light, dark, or the system's choice; the activity bar's button cycles the same three." },
    { key: "_accent", label: "Accent colour", kind: "info", text: "accent", help: "Buttons, selections and highlights take the accent of System Settings ▸ Appearance; the app's teal when it is multicolour." },
    { key: "zoom", label: "Size of the views", kind: "select", default: 1, options: [[0.8, "80 %"], [0.9, "90 %"], [1, "100 %"], [1.1, "110 %"], [1.25, "125 %"], [1.5, "150 %"]], help: "The text and controls of every window; View ▸ Zoom In and Zoom Out change it too." },
    { key: "tips", label: "Hover cards", kind: "bool", default: true, help: "Explain controls and parameters when the pointer rests on them, from the parameter dictionary." },
    { key: "_palette", label: "Plot colours", kind: "info", text: "One colour-blind-safe palette, with the same meaning in light and dark (Yobs, Ycal, the difference, each phase).", help: "" },
  ] },
  { id: "editor", title: "Editor", glyph: "✎", rows: [
    { key: "editorFont", label: "Font", kind: "select", default: "system", options: [["system", "the system's monospaced (SF Mono)"], ["Menlo", "Menlo"], ["Monaco", "Monaco"], ["SF Mono", "SF Mono"], ["Courier New", "Courier New"], ["Andale Mono", "Andale Mono"], ["PT Mono", "PT Mono"]], help: "The Editor's typeface: monospaced, so that FullProf's columns stay aligned; every one listed comes with macOS." },
    { key: "editorFontSize", label: "Font size", kind: "number", default: 12.5, min: 10, max: 20, step: 0.5, unit: "px", help: "The size of the text in the Editor, in a monospaced font." },
    { key: "tabWidth", label: "Tab width", kind: "select", default: 4, options: [[2, "2 spaces"], [4, "4 spaces"], [8, "8 spaces"]], overridable: true, help: "How wide a tab character is drawn." },
    { key: "lineNumbers", label: "Line numbers", kind: "bool", default: true, help: "The number of each line in the gutter, as the core's checks name them." },
    { key: "wrap", label: "Wrap long lines", kind: "bool", default: false, help: "Long lines continue on the next line instead of scrolling sideways." },
    { key: "autosave", label: "Save files by themselves", kind: "bool", default: true, help: "Two seconds after the last change, and when a window closes or the app quits; versions are kept as before (L224)." },
    { key: "_minimap", label: "Minimap", kind: "info", text: "None, decided (L337): control files are short; the outline and Next Problem (F8) take its place.", help: "" },
  ] },
  { id: "plots", title: "Plots", glyph: "⌇", rows: [
    { key: "refineXAxis", label: "Refine's plot against", kind: "select", default: "2theta", options: [["2theta", "2θ"], ["d", "d"], ["q", "Q"]], overridable: true, help: "The x of a Refine document's plot when it opens; d and Q are the core's conversion at λ1." },
    { key: "yScale", label: "Intensity scale", kind: "select", default: "linear", options: [["linear", "linear"], ["sqrt", "√"], ["log", "log"]], overridable: true, help: "The intensity axis of Data's and Refine's plots when they open." },
    { key: "yobsLine", label: "Yobs as a line", kind: "bool", default: false, overridable: true, help: "Yobs drawn as a line; markers by default, as in WinPLOTR and in journals." },
    { key: "markerSize", label: "Marker size", kind: "select", default: 1, options: [[0.75, "small"], [1, "normal"], [1.5, "large"], [2, "larger"]], overridable: true, help: "The size of the points every plot draws on screen (Yobs's markers, a data file's points); a figure's are set in Publish." },
    { key: "_publication", label: "Publication style", kind: "info", text: "Publish's presets (journals, thesis), set in the Publish module.", help: "" },
  ] },
  { id: "formats", title: "Data formats", glyph: "▤", rows: [
    { key: "_formats", label: "The core's readers", kind: "formats", help: "A file is recognised by its content first, then its extension; the core decides which reader reads it." },
    { key: "_formatOrder", label: "Order of the readers", kind: "info", text: "The core's, for now: when two readers read a file equally well, a preference set here will decide (app request 51).", help: "" },
    { key: "_formatOptions", label: "Each format's options", kind: "info", text: "To come with the core's declarations (app request 51): a GSAS file's bank, a column file's columns, a wavelength where the file states none, a time-of-flight unit.", help: "" },
  ] },
  { id: "refinement", title: "Refinement", glyph: "↻", rows: [
    { key: "cycles", label: "Default cycles", kind: "number", default: "", min: 1, max: 1000, step: 1, overridable: true, help: "Used when a control file states no cycles; a control file's own Cycles always decide." },
    { key: "pauseEvery", label: "Pause every", kind: "number", default: "", min: 1, max: 1000, step: 1, unit: "cycles", overridable: true, help: "A checkpoint every k cycles, to look before going on (W11); empty runs to the end." },
    { key: "snapshotBeforeRun", label: "A version before each run", kind: "bool", default: false, help: "The control file kept in Versions/ before every run, besides the version each accepted result makes." },
    { key: "notifyRunEnd", label: "Notify when a run ends", kind: "bool", default: true, help: "A system notification when a refinement ends or pauses while the app is in the background." },
    { key: "_convergence", label: "Convergence", kind: "info", text: "The control file's own (Refine's toolbar); the core has no app-wide default to override it.", help: "" },
  ] },
  { id: "external", title: "External programs", glyph: "⇄", rows: [
    { key: "fullprofPath", label: "FullProf", kind: "program", program: "fullprof", default: "", file: true, help: "FullProf's fp2k program, for Refine ▸ Compare with ▸ FullProf; Check runs it as a comparison will and shows its banner." },
    { key: "gsas2Path", label: "GSAS-II", kind: "program", program: "gsas2", default: "", file: false, help: "GSAS-II's installation folder (the one holding GSAS-II/GSASII), for Refine ▸ Compare with ▸ GSAS-II; Check imports GSASIIscriptable in its Python." },
    { key: "gsas2Python", label: "GSAS-II's Python", kind: "path", default: "", file: true, help: "Empty: the installation's bin/python; or another Python in which GSAS-II can be imported (pip, conda)." },
    { key: "_compare", label: "Comparisons", kind: "info", text: "Refine ▸ Compare with ▸ FullProf or GSAS-II runs both from the same control file and shows them side by side; it waits for the core to write their inputs and read their outputs (app requests 49, 50).", help: "" },
  ] },
  { id: "advanced", title: "Advanced", glyph: "⚑", rows: [
    { key: "logLines", label: "Log lines shown", kind: "number", default: 200, min: 50, max: 5000, step: 50, help: "How many of the latest messages the Log keeps on screen." },
    { key: "_folder", label: "Settings file", kind: "button", action: "folder", button: "Show in Finder", help: "The app's settings, recent projects and last session, beside each other." },
    { key: "_diagnostics", label: "Diagnostics", kind: "button", action: "diagnostics", button: "Diagnostics Report…", help: "The versions, macOS and the Mac, the windows, the settings (no paths, no names) and the recent log, to copy into a report of a problem." },
    { key: "_reset", label: "Reset all settings", kind: "button", action: "reset", button: "Reset All…", help: "Every setting, shortcut and toolbar back to its default; projects, their files and their own settings are not touched." },
  ] },
];

/** The settings a row writes (no "_" keys: those are information and actions). */
export const PREF_KEYS = [];
/** Every setting's default. */
export const DEFAULTS = {};
/** The settings a project may set its own value of (§4: "a project can override a few, and says so"). */
export const OVERRIDABLE = [];
function index() {
  const rows = PANES.flatMap((p) => p.rows);
  PREF_KEYS.length = 0; PREF_KEYS.push(...rows.filter((r) => !r.key.startsWith("_")).map((r) => r.key));
  for (const k of Object.keys(DEFAULTS)) delete DEFAULTS[k];
  Object.assign(DEFAULTS, Object.fromEntries(rows.filter((r) => !r.key.startsWith("_")).map((r) => [r.key, r.default])));
  OVERRIDABLE.length = 0; OVERRIDABLE.push(...rows.filter((r) => r.overridable).map((r) => r.key));
}
index();
/** A module's own pane (E8's module slot, app 0.38.0, §4, §11): `{ id, title, glyph, rows }` with rows as above, its
    keys named by the module (`explore.trail`); shown in the Settings window after the app's own panes, its settings
    kept, merged and overridden as the app's are. A pane already registered is replaced. */
export function registerPane(pane) {
  if (!pane || !pane.id || !Array.isArray(pane.rows)) return;
  const k = PANES.findIndex((p) => p.id === pane.id);
  const p = { glyph: "◇", ...pane, module: true };
  if (k >= 0) PANES.splice(k, 1, p); else PANES.push(p);
  index();
}

/** A setting's value in force: the project's own, else the app's, else the default. */
export function effective(key, app = {}, project = null) {
  if (project && Object.prototype.hasOwnProperty.call(project, key) && OVERRIDABLE.includes(key)) return project[key];
  const v = app[key];
  return v === undefined || v === null ? DEFAULTS[key] : v;
}

/** The settings a window writes: its own state, with the preferences taken from what is on disk, except those this
    window has just changed itself (`own`), so that a window holding an old copy never undoes the Settings window. */
export function mergeForSave(mine = {}, disk = {}, own = []) {
  const out = { ...mine };
  for (const k of PREF_KEYS) {
    if (own.includes(k)) continue;
    // a preference absent on disk is at its default: an old copy must not bring it back (0.30.0's self-check)
    if (Object.prototype.hasOwnProperty.call(disk, k)) out[k] = disk[k]; else delete out[k];
  }
  return out;
}

/** What the Settings window writes: what is on disk, with every preference as this window holds it, a preference it
    no longer holds (put back to its default) removed, and the other keys it changed itself (`own`: the shortcuts). */
export function settingsWindowMerge(mine = {}, disk = {}, own = []) {
  const out = { ...disk };
  for (const k of [...PREF_KEYS, ...own]) {
    if (Object.prototype.hasOwnProperty.call(mine, k)) out[k] = mine[k]; else delete out[k];
  }
  return out;
}

/** What opens at launch (§4, app 0.46.0): "session", "recent" or "welcome". A settings file of an earlier version
    keeps its choice: Reopen windows and tabs at launch switched off is the welcome page. */
export function launchOpens(settings = {}) {
  const v = settings.launchOpens;
  if (v === "session" || v === "recent" || v === "welcome") return v;
  return settings.restoreWindows === false ? "welcome" : "session";
}

/** The Editor's font family as CSS (§4, app 0.46.0): the chosen family, then the system's monospaced faces after it. */
export function editorFontFamily(choice) {
  const fallback = "ui-monospace, Menlo, monospace";
  if (!choice || choice === "system") return fallback;
  return `"${String(choice).replace(/"/g, "")}", ${fallback}`;
}
