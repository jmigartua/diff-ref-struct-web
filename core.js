/* diff-ref-struct · the window's shared state and small helpers (the services the modules may use, DESIGN.md §4).
   No physics: numbers come from the core through the app's commands, names from the parameter dictionary (the
   core's table, asked once at start). */

import { effective } from "./settings.js";

export const $ = (id) => document.getElementById(id);
export const esc = (s) => String(s ?? "").replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/"/g, "&quot;");
export const el = (tag, cls, html) => { const e = document.createElement(tag); if (cls) e.className = cls; if (html != null) e.innerHTML = html; return e; };

/** Inside the Mac app (Tauri). A page that hosts the modules otherwise (the web host, W6) sets `window.__DRS_HOST__`
    and leaves `__TAURI__` alone. */
export const IN_APP = !!window.__TAURI__;
/** The host (W6, PL13q): what answers the modules' `invoke`, `listen`, the open and save panels and `ask`. The page's
    own object when it sets one (`window.__DRS_HOST__ = { invoke, listen, open, save, openAny, ask }`, each optional
    but `invoke`), else Tauri's inside the app, else none. Read at each call, so that a host set after this file is
    evaluated is seen. */
const tauriHost = () => (window.__TAURI__ ? {
  invoke: (cmd, args) => window.__TAURI__.core.invoke(cmd, args),
  listen: (ev, fn) => window.__TAURI__.event.listen(ev, (e) => fn(e.payload)),
  open: (opts) => window.__TAURI__.dialog.open(opts),
  save: (opts) => window.__TAURI__.dialog.save(opts),
  openAny: (title) => (window.__TAURI__.dialog && window.__TAURI__.dialog.openAny ? window.__TAURI__.dialog.openAny() : window.__TAURI__.core.invoke("open_panel", { title })),
  ask: (message, title, ok) => window.__TAURI__.dialog.ask(message, { title, kind: "warning", okLabel: ok, cancelLabel: "Cancel" }),
} : null);
export const host = () => window.__DRS_HOST__ || tauriHost();
/** Whether a core answers `invoke` (the app's through Tauri, or the page's host): what a module asks before asking the
    core; `IN_APP` is for what belongs to the Mac app alone. */
export const hasCore = () => { const h = host(); return !!(h && typeof h.invoke === "function"); };
export const APP = "diff-ref-struct";
export const SITE = "https://jmigartua.github.io/diff-ref-struct-site/";
export const ISSUES = "https://github.com/jmigartua/diff-ref-struct/issues";
/** Pages of the project site the Help menu and the cards link to: [path, title]. */
export const DOCS = [
  ["manual/parameters/index.html", "Parameters: the theory"],
  ["design/index.html", "Design: the decisions"],
  ["roadmap/index.html", "Roadmap"],
  ["releases/index.html", "Release notes"],
  ["log/index.html", "Project log"],
  ["index.html", "Project site"],
];

/** The window's state. */
/** A setting's value in force (APP_CONVENTIONS §4): the open project's own (File ▸ Project Settings…), else the app's
    (the Settings window), else its default (settings.js). */
export const setting = (key) => effective(key, state.settings, state.projectState && state.projectState.settings);

export const state = {
  project: null,          // { root, name, files } — the folder on disk (R20)
  projectState: {},       // .drs/project.json: layout and per-module state
  tabs: new Map(),        // tab id → { id, module, path, title, docId? }
  activeTab: null,
  activeModule: "data",   // the module whose sidebar section, menu and panes are shown
  docs: new Map(),        // Refine: doc id → { id, view, refine, sel, plot, el }
  settings: { theme: "system", tips: true, cycles: "", eps: "", zoom: 1 },
  recent: [],             // recent projects (folders)
  versions: ["?", "?"],
  dictionary: { groups: [], parameters: [] },
  byId: new Map(),
  ext: "drc",             // our control file's extension, from the Rust side (DESIGN.md §9, answer 1)
  log: [],                // the status messages, newest last (C9.3)
};

/** A tiny event bus: "project", "tabs", "module", "doc", "select", "refine", "layout", "theme". */
const handlers = new Map();
export const bus = {
  on(ev, fn) { if (!handlers.has(ev)) handlers.set(ev, []); handlers.get(ev).push(fn); return () => { const l = handlers.get(ev); const i = l.indexOf(fn); if (i >= 0) l.splice(i, 1); }; },
  emit(ev, arg) { for (const fn of [...(handlers.get(ev) || [])]) fn(arg); },
};

export const activeTab = () => (state.activeTab == null ? null : state.tabs.get(state.activeTab) || null);
/** The Refine document of the active tab, if the active tab is a Refine document. */
export const activeDoc = () => { const t = activeTab(); return t && t.module === "refine" && t.docId != null ? state.docs.get(t.docId) || null : null; };
export const active = activeDoc;

/** Linked selection (E8, app 0.38.0): a parameter chosen in one view, selected in the others. `address` is the core's
    (`phase1.cell_a`), `docId` the Refine document it belongs to (null: the active one), `from` the view that chose it,
    which does not answer itself. Refine selects it in its outline and inspector; the Editor's line waits for the core's
    place of an address in the control file (app/REQUESTS.md 40). */
export function selectParameter(address, docId = null, from = "") { if (address) bus.emit("parameter", { address, docId, from }); }

export function status(msg, error = false) {
  const s = $("status"); if (s) { s.textContent = msg; s.classList.toggle("error", error); }
  state.log.push({ t: new Date(), msg, error });
  if (state.log.length > 500) state.log.shift();
  bus.emit("log");
}

// ------------------------------------------------------------------ the app's commands
export async function invoke(cmd, args) {
  const h = host();
  if (!h || typeof h.invoke !== "function") throw new Error("no core answers: not inside the app, and the page set no host");
  return h.invoke(cmd, args);
}
export async function listen(ev, fn) {
  const h = host();
  if (!h || typeof h.listen !== "function") return () => {};
  return h.listen(ev, fn);
}
/** False while the open project is a reference folder opened as it is (W9): nothing is written into it. */
export const writable = () => !(state.project && state.project.readonly);
/** The host's open panel (W6): `opts` as Tauri's dialog takes them ({ multiple, directory, title, defaultPath,
    filters }); null when the host has none. */
export async function pickFiles(opts) {
  const h = host();
  if (!h || typeof h.open !== "function") return null;
  return h.open(opts);
}
/** The host's save panel (W6): the path chosen, or null when the host has none or the user cancelled. */
export async function pickSavePath(opts) {
  const h = host();
  if (!h || typeof h.save !== "function") return null;
  return h.save(opts);
}
export async function pickFolder(title = "Open a project folder", defaultPath = undefined) {
  return pickFiles({ directory: true, multiple: false, title, defaultPath });
}
/** File ▸ Open… (§3, app 0.46.0): AppKit's open panel, a folder or a file in one (the tour answers it as it answers
    the other dialogs). */
export async function pickAny(title = "Open a project's folder, a control file or a data file") {
  const h = host();
  if (!h) return null;
  if (typeof h.openAny === "function") return h.openAny(title);
  return pickFiles({ multiple: false, title });
}
export async function pickFile(filters) {
  return pickFiles({ multiple: false, filters });
}
export async function ask(message, title, ok = "OK") {
  const h = host();
  if (!h || typeof h.ask !== "function") return window.confirm(message);
  return h.ask(message, title, ok);
}
/* The dialogs below are reused question after question. A dialog's "close" event comes a moment after close(): when a
   command asks twice in a row (Smooth: the points, then the order), the first question's close arrived once the
   dialog was open for the second, and answered it "cancelled" (found by the tour, scene 1.7: Smooth, Mask a Range and
   Delete a Range never kept an operation). A close while the dialog is open again is ignored. */
/** A text asked in the app's own dialog (never window.prompt, which Tauri's WebKit view does not show; app 0.10.1).
    opts: { label ("Name"), ok ("OK"), placeholder }. Null when cancelled. */
export function askText(title, text, value = "", opts = {}) {
  return new Promise((resolve) => {
    const d = $("dlg-name"); if (!d) { resolve(null); return; }
    $("dlg-name-title").textContent = title; $("dlg-name-text").textContent = text || "";
    $("dlg-name-label").textContent = opts.label || "Value"; $("dlg-name-ok").textContent = opts.ok || "OK";
    const inp = $("dlg-name-input"); inp.value = value == null ? "" : String(value); inp.placeholder = opts.placeholder || "";
    const done = (v) => { d.onclose = null; d.close(); resolve(v); };
    $("dlg-name-ok").onclick = () => done(inp.value);
    $("dlg-name-cancel").onclick = () => done(null);
    d.onclose = () => { if (!d.open) resolve(null); };   // a close while open again is the previous question's (below)
    inp.onkeydown = (e) => { if (e.key === "Enter") { e.preventDefault(); done(inp.value); } };
    d.showModal(); inp.focus(); inp.select();
  });
}
/** A number asked in the app's own dialog: null when cancelled, NaN when not a number. */
export async function askNumber(title, text, value, opts = {}) {
  const v = await askText(title, text, value, opts); if (v == null) return null;
  const n = parseFloat(String(v).trim().replace(",", ".")); return isFinite(n) ? n : NaN;
}
/** One of several choices, in the app's own dialog: options [{ value, label }]; the value chosen, or null. */
export function pickOne(title, text, options, selected = null, ok = "OK") {
  return new Promise((resolve) => {
    const d = $("dlg-pick"); if (!d) { resolve(null); return; }
    $("dlg-pick-title").textContent = title; $("dlg-pick-text").textContent = text || ""; $("dlg-pick-ok").textContent = ok;
    const sel = $("dlg-pick-select"); sel.innerHTML = options.map((o) => `<option value="${esc(String(o.value))}" ${String(o.value) === String(selected) ? "selected" : ""}>${esc(o.label)}</option>`).join("");
    const done = (v) => { d.onclose = null; d.close(); resolve(v); };
    $("dlg-pick-ok").onclick = () => done(sel.value);
    $("dlg-pick-cancel").onclick = () => done(null);
    d.onclose = () => { if (!d.open) resolve(null); };   // a close while open again is the previous question's (below)
    sel.onkeydown = (e) => { if (e.key === "Enter") { e.preventDefault(); done(sel.value); } };
    d.showModal(); sel.focus();
  });
}
/** Several fields asked at once in the app's own dialog (app 0.11.1): fields [{ id, label, value, placeholder, help,
    kind: "text" | "number" | "textarea" }]; the values by id, or null when cancelled. A number field returns a
    number, null when blank, NaN when not a number. */
export function askForm(title, text, fields, ok = "OK") {
  return new Promise((resolve) => {
    const d = $("dlg-form"); if (!d) { resolve(null); return; }
    $("dlg-form-title").textContent = title; $("dlg-form-text").textContent = text || ""; $("dlg-form-ok").textContent = ok;
    const box = $("dlg-form-fields"); box.innerHTML = "";
    for (const f of fields) {
      const row = el("label", "row"); row.innerHTML = `<span>${esc(f.label)}</span>`;
      const inp = f.kind === "textarea" ? el("textarea") : el("input");
      if (f.kind !== "textarea") inp.type = "text";
      inp.id = "dlg-form-" + f.id; inp.value = f.value == null ? "" : String(f.value); inp.placeholder = f.placeholder || ""; inp.spellcheck = false;
      if (f.kind === "textarea") { inp.rows = 3; inp.cols = 44; } else inp.size = f.kind === "number" ? 10 : 34;
      row.appendChild(inp); box.appendChild(row);
      if (f.help) box.appendChild(el("p", "help", esc(f.help)));
    }
    const values = () => { const out = {}; for (const f of fields) { const v = $("dlg-form-" + f.id).value; out[f.id] = f.kind === "number" ? (v.trim() === "" ? null : (isFinite(parseFloat(v.trim().replace(",", "."))) ? parseFloat(v.trim().replace(",", ".")) : NaN)) : v; } return out; };
    const done = (v) => { d.onclose = null; d.close(); resolve(v); };
    $("dlg-form-ok").onclick = () => done(values());
    $("dlg-form-cancel").onclick = () => done(null);
    d.onclose = () => { if (!d.open) resolve(null); };   // a close while open again is the previous question's (below)
    box.querySelectorAll("input").forEach((inp) => { inp.onkeydown = (e) => { if (e.key === "Enter") { e.preventDefault(); done(values()); } }; });
    d.showModal(); const first = box.querySelector("input, textarea"); if (first) { first.focus(); first.select && first.select(); }
  });
}
/** A question whose body is a list to read before answering (Refine ▸ Guided Plan…'s steps): the form dialog with HTML in
    place of fields. True for OK, false for Cancel. */
export function askHtml(title, text, html, ok = "OK") {
  return new Promise((resolve) => {
    const d = $("dlg-form"); if (!d) { resolve(false); return; }
    $("dlg-form-title").textContent = title; $("dlg-form-text").textContent = text || ""; $("dlg-form-ok").textContent = ok;
    $("dlg-form-fields").innerHTML = html;
    const done = (v) => { d.onclose = null; d.close(); resolve(v); };
    $("dlg-form-ok").onclick = () => done(true);
    $("dlg-form-cancel").onclick = () => done(false);
    d.onclose = () => { if (!d.open) resolve(false); };
    d.showModal(); $("dlg-form-ok").focus();
  });
}
export async function openExternal(url) {
  if (IN_APP && window.__TAURI__.opener) return window.__TAURI__.opener.openUrl(url);
  window.open(url, "_blank");
}
export async function revealInFinder(path) {
  if (IN_APP && window.__TAURI__.opener) return window.__TAURI__.opener.revealItemInDir(path);
}
export const docUrl = (path) => SITE + path;

// ------------------------------------------------------------------ numbers
/** value(esd) in the crystallographers' notation: the esd's first significant digits in brackets, two of them when
    they start with 1 (so 0.000082 → (82), 0.0123 → (12), 0.03 → (3)). */
export function withEsd(v, e) {
  if (!(e > 0) || !isFinite(e)) return fmt(v);
  const lead = Math.floor(Math.log10(e));
  let digits = Math.max(0, -lead + (e / Math.pow(10, lead) < 2 ? 1 : 0));
  digits = Math.min(digits, 10);
  return `${v.toFixed(digits)}(${Math.round(e * Math.pow(10, digits))})`;
}
/** A value without esd: plain when it reads well, exponential when it would not. */
export function fmt(v) {
  if (v == null || !isFinite(v)) return "";
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 1e6 || a < 1e-4) return v.toExponential(4).replace(/\.?0+e/, "e");
  return String(parseFloat(v.toPrecision(7)));
}
export function fmtSize(n) {
  if (n < 1024) return `${n} B`;
  if (n < 1024 * 1024) return `${(n / 1024).toFixed(0)} kB`;
  return `${(n / 1024 / 1024).toFixed(1)} MB`;
}

// ------------------------------------------------------------------ the dictionary (from the core)
export function setDictionary(d) {
  state.dictionary = { groups: (d.groups || []).map((g) => (Array.isArray(g) ? { key: g[0], title: g[1] } : g)), parameters: d.parameters || [] };
  state.byId = new Map(state.dictionary.parameters.map((p) => [p.id, p]));
}
export const entry = (id) => state.byId.get(id) || null;
export const groupTitle = (key) => (state.dictionary.groups.find((g) => g.key === key) || {}).title || key;
export const nameOf = (id) => { const e = entry(id); return e ? e.name : id; };

/** Renders LaTeX into an element with KaTeX when it is loaded (vendored, no network), else as text. With `symbol`
    (the parameter's own, as the dictionary writes it), its occurrences in the formula are marked (class `hl`), so the
    eye finds the parameter in the equation (Josu, 0.8.0). */
export function tex(target, src, display = false, symbol = null) {
  if (!src) { target.textContent = ""; return; }
  const marked = symbol ? markSymbol(src, symbol) : src;
  if (window.katex) { try { window.katex.render(marked, target, { throwOnError: false, displayMode: display, trust: true, strict: false }); return; } catch (e) { /* fall through */ } }
  target.textContent = src;
}
/** "U_i, V_i, W_i" → each token wrapped in \htmlClass{hl}{…} where it stands alone in the equation ("W" not in "W_i"). */
export function markSymbol(eq, symbol) {
  let out = eq;
  for (let tok of String(symbol).split(",")) {
    tok = tok.trim();
    if (tok.includes("=")) tok = tok.split("=")[0].trim();        // "S_c = -2s/R": the symbol is S_c
    if (!tok) continue;
    const re = new RegExp("(?<![A-Za-z\\\\])" + tok.replace(/[.*+?^${}()|[\]\\]/g, "\\$&") + "(?![A-Za-z0-9_{])", "g");
    out = out.replace(re, (m) => `\\htmlClass{hl}{${m}}`);
  }
  return out;
}
