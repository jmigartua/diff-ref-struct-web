/* diff-ref-struct Web · the browser host (I10 steps 2 to 5, R34; web/DESIGN.md): the Plots module and the shell's service
   files loaded as they are from app/ui/ (relative imports, no copy, no bundler), behind the dispatcher onto the core
   in wasm (dispatch.js). The page holds what the app's window holds for the module: the tabs (dockview, layout.js),
   the panels the module draws into (its members table, its member), the status line, the dialogs; and the files: a
   folder picked or files dropped, read in the browser, never uploaded. The browser computes nothing (L139). */
import { $, el, esc, state, bus, status, fmtSize } from "../core.js";
import { register, run, cmd, enabled } from "../registry.js";
import { initDock, ensureStart, openTab, activateTab, relayout } from "../layout.js";
import { plots } from "../modules/plots/index.js";
import { createDispatcher, loadPackage } from "./dispatch.js";
import { createFiles, fromDirectoryHandle, fromFileList, fromDrop } from "./files.js";
import { buildTree, groupByKind, visibleRows, filesUnder, clickSelect } from "./tree.js";
import { clampWidth, dragged, keyed, kept, keep, RAIL } from "./separator.js";
import { diffListing, nothing, growth, liveLine, INTERVALS, DEFAULT_INTERVAL } from "./live.js";
import { encode, decode, missingFor } from "./hash.js";

const files = createFiles();
let dispatcher = createDispatcher({ pkg: null, fs: files });
const core = { ready: false, version: null, why: "", checked: false };
let folderHandle = null;   // the folder's handle where the browser gives one: live mode reads it again (step 4)

/** The page is the modules' host (W6, Plots' step 21): core.js asks it for `invoke` (the core in wasm, dispatch.js) and
    the open and save panels (the files read here; a download). `IN_APP` stays false: this is not the Mac app. */
window.__DRS_HOST__ = {
  invoke: (cmd, args) => dispatcher.invoke(cmd, args || {}),
  listen: async () => () => {},
  open: (opts) => pick(opts || {}),
  // a save panel: the path the module proposes, under the folder shown; the dispatcher downloads what is written there
  save: async (opts) => { const d = (opts && opts.defaultPath) || "download"; return d.startsWith("/") ? d : files.root ? files.absOf(d) : "/local/" + d; },
  ask: async (message) => window.confirm(message),
};

// ------------------------------------------------------------------ the open panel, answered from the files read here
/** The paths the host hands the next open panel (Open as Plot Set, a drop), else the tree's selection. */
let pending = null;
async function pick(opts) {
  let paths = pending; pending = null;
  if (!paths) {
    const exts = (opts.filters || []).flatMap((f) => f.extensions || []).filter((x) => x !== "*").map((x) => x.toLowerCase());
    const every = (opts.filters || []).some((f) => (f.extensions || []).includes("*"));
    paths = tree.sel.map(files.absOf).filter((p) => every || !exts.length || exts.includes(p.split(".").pop().toLowerCase()));
    if (!paths.length) {
      status(`${opts.title ? opts.title + ": " : ""}choose ${exts.length && !every ? `a .${exts.join(" or .")} file` : "the files"} in the folder's tree first (click, ⌘-click, ⇧-click), then the command again. Files are picked with Open Folder… or dropped on the page.`, true);
      return null;
    }
  }
  return opts.multiple ? paths : paths[0];
}

// ------------------------------------------------------------------ the core
async function startCore() {
  try {
    const got = await loadPackage({
      base: import.meta.url,
      getJson: async (url) => { const r = await fetch(url, { cache: "no-cache" }); if (!r.ok) throw `HTTP ${r.status}`; return r.json(); },
      load: (url) => import(url),
    });
    if (got.pkg) { dispatcher = createDispatcher({ pkg: got.pkg, fs: files }); core.ready = true; core.version = got.version; }
    core.why = got.why || "";
  } catch (e) { core.why = `the package did not load: ${e}`; }
  core.checked = true;
  paintCore();
}
function paintCore() {
  const c = $("core-chip");
  c.textContent = core.ready ? `core ${core.version || ""} in WebAssembly` : "core: waiting for W1";
  c.title = core.ready ? "The core runs in this page as WebAssembly: every number is its own." : `The core is not in this page yet (request W1: the core compiled to WebAssembly). The page opens and lists your files; every view waits for the core. (${core.why})`;
  c.classList.toggle("off", !core.ready);
  document.querySelectorAll(".wh-core-note").forEach((n) => { n.hidden = core.ready; });
}

// ------------------------------------------------------------------ the folder's tree
const tree = { root: buildTree([]), open: new Set(), shut: new Set(), sel: [], anchor: null, by: "folders" };
try { tree.by = localStorage.getItem("drs-web.tree.by") === "kinds" ? "kinds" : "folders"; } catch (e) { /* none kept */ }
/** What each file is, the core's `describe` (step 3), by its path in the folder: { kind, id, name }. */
const kinds = new Map();
const kindOf = (rel) => kinds.get(rel) || null;
function setFiles() {
  tree.open = new Set(buildTree(files.list()).children.filter((c) => c.dir).map((c) => c.path));   // the first level open
  tree.sel = tree.sel.filter((p) => files.has(files.absOf(p)));
  for (const k of [...kinds.keys()]) if (!files.has(files.absOf(k))) kinds.delete(k);
  state.project = files.root ? { root: files.root, name: files.name, readonly: false, files: [] } : null;
  rebuild();
}
/** The tree in the way chosen: by folders, or by the core's kinds (a group per kind, open unless closed). */
function rebuild() {
  tree.root = tree.by === "kinds" ? groupByKind(files.list(), kindOf, files.name || "") : buildTree(files.list(), files.name || "");
  renderTree();
}
const openSet = () => (tree.by === "kinds" ? new Set(tree.root.children.map((g) => g.path).filter((p) => !tree.shut.has(p))) : tree.open);
/** Every file described by the core, one after the other, the tree redrawn as they come (a new folder starts again). */
let describing = 0;
async function describeAll() {
  const gen = ++describing;
  if (!dispatcher.ready) return;
  const list = files.list(); let n = 0;
  for (const f of list) {
    if (gen !== describing) return;
    if (!kinds.has(f.path)) {
      try { kinds.set(f.path, await dispatcher.kindOf(files.absOf(f.path))); }
      catch (e) { kinds.set(f.path, { kind: "unknown", id: "unread", name: `not described: ${e}` }); }
    }
    n += 1;
    if (n % 10 === 0 || n === list.length) { if (tree.by === "kinds") rebuild(); else renderTree(); }
  }
}
function renderTree() {
  const box = $("tree"); box.innerHTML = "";
  $("folder-name").textContent = files.name || "";
  $("tree-count").textContent = tree.root.files ? `${tree.root.files} file${tree.root.files === 1 ? "" : "s"} · ${fmtSize(tree.root.size)}` : "";
  $("btn-openset").disabled = !tree.sel.length;
  $("btn-openset").textContent = tree.sel.length > 1 ? `Open ${tree.sel.length} as Plot Set` : "Open as Plot Set";
  if (!tree.root.files) { box.appendChild(el("div", "empty small", "No folder yet: <b>Open Folder…</b>, or drop a folder or files anywhere on the page.")); return; }
  const open = openSet(), byKind = tree.by === "kinds";
  const rows = visibleRows(tree.root, open);
  const order = rows.filter((r) => !r.node.dir).map((r) => r.node.path);
  for (const { node, depth } of rows) {
    const r = el("div", "wh-row" + (node.dir ? " dir" : "") + (!node.dir && tree.sel.includes(node.path) ? " sel" : ""));
    r.setAttribute("role", "treeitem"); r.setAttribute("aria-level", String(depth + 1));
    if (node.dir) r.setAttribute("aria-expanded", String(open.has(node.path))); else r.setAttribute("aria-selected", String(tree.sel.includes(node.path)));
    r.style.paddingLeft = `${8 + depth * 14}px`;
    const k = node.dir ? null : kindOf(node.path);
    const what = k ? (k.kind === "unknown" ? "not recognised by the core's readers" : `${k.name} (${k.kind === "data" ? "data" : "a project's file"}, the core's describe)`) : core.ready ? "being described by the core" : "described when the core is in the page";
    const setOk = !byKind || node.kind === "data" || node.kind === "pending";
    r.title = node.dir ? (byKind ? `${node.name}: ${node.files} file${node.files === 1 ? "" : "s"}, as the core's describe says.` : `${node.path}/: ${node.files} file${node.files === 1 ? "" : "s"}. Double-click opens or closes it; "set" opens its data files as one plot set.`) : `${node.path} · ${fmtSize(node.size)} · ${what}. Click to choose (⌘ for several, ⇧ for a range); double-click opens it in Plots.`;
    r.innerHTML = node.dir
      ? `<span class="wh-caret">${open.has(node.path) ? "▾" : "▸"}</span><span class="wh-name">${esc(node.name)}</span><span class="sub">${node.files}</span>${setOk ? `<button class="small wh-dirset" title="Open the data files of ${esc(node.name)} as one plot set" aria-label="Open ${esc(node.name)} as a plot set">set</button>` : ""}`
      : `<span class="wh-caret"></span><span class="wh-name">${esc(node.name)}</span>${k && !byKind ? `<span class="wh-kind${k.kind === "unknown" ? " none" : ""}">${esc(k.kind === "unknown" ? "?" : k.id)}</span>` : ""}<span class="sub">${fmtSize(node.size)}</span>`;
    r.onclick = (e) => {
      if (node.dir) {
        if (e.target.classList.contains("wh-dirset")) { openPaths(filesUnder(node)); return; }
        const set = byKind ? tree.shut : tree.open, add = byKind ? open.has(node.path) : !open.has(node.path);
        if (add) set.add(node.path); else set.delete(node.path);
        renderTree(); return;
      }
      Object.assign(tree, clickSelect(tree.sel, tree.anchor, node.path, order, { toggle: e.metaKey || e.ctrlKey, range: e.shiftKey }));
      renderTree();
    };
    if (!node.dir) r.ondblclick = () => openPaths([node.path]);
    box.appendChild(r);
  }
}

// ------------------------------------------------------------------ files into the module
/** What a file is for routing: the core's `describe` (asked now if the tree has not yet), else, on the stub, its
    extension as the module's own dialogs filter (.drv, .buf). */
async function routeKind(rel) {
  let k = kindOf(rel);
  if (!k && dispatcher.ready) { try { k = await dispatcher.kindOf(files.absOf(rel)); kinds.set(rel, k); } catch (e) { k = null; } }
  if (k) return k;
  return /\.drv$/i.test(rel) ? { kind: "project", id: "plot_set" } : /\.buf$/i.test(rel) ? { kind: "project", id: "buffer_file" } : null;
}
const SET_FILES = ["plot_set", "buffer_file"];
/** Files of the folder opened in Plots (step 3): each plot set file or WinPLOTR buffer among them as itself (its
    members the folder's files it names, taken by the core: `plots_open_set_from`, `plots_open_buffer_from`); else the
    data files as one plot set (`plots_open_from`). What is neither is left out and named with what the core says it is. */
async function openPaths(rels) {
  if (!rels.length) return;
  const ks = []; for (const r of rels) ks.push([r, await routeKind(r)]);
  const sets = ks.filter(([, k]) => k && SET_FILES.includes(k.id));
  if (sets.length) { for (const [r, k] of sets) await openSetFile(r, k.id); return; }
  const data = ks.filter(([, k]) => !k || k.kind === "data").map(([r]) => r);
  const left = ks.filter(([, k]) => k && k.kind !== "data").map(([r, k]) => `${r} (${k.kind === "unknown" ? "no reader recognises it" : k.name})`);
  if (!data.length) { status(`Nothing to open as a plot set: ${left.join("; ")}.`, true); return; }
  pending = data.map(files.absOf);
  const tab = await opened(() => run("plots.opencsv"));   // Open CSV as a Set…: "a CSV of several columns (or several data files)" → plots_open_files
  pending = null;
  if (tab) sources.set(tab.id, { by: "files", files: data });
  if (left.length) { const s = $("status"); status(`${s.textContent} Left out, not data: ${left.join("; ")}.`, s.classList.contains("error")); }
}
/** Where each tab's set came from (tab id → the hash's `set`): what the address keeps and live mode grows. */
const sources = new Map();
/** The tab a module command opened (the one new in state.tabs after it), or null. */
async function opened(fn) {
  const before = new Set(state.tabs.keys());
  await fn();
  const id = [...state.tabs.keys()].find((k) => !before.has(k));
  return id ? state.tabs.get(id) : null;
}
/** A set file or a buffer of the folder opened by the module, as the app's openFile does. */
async function openSetFile(rel, id) {
  // the file already open: its tab brought forward, as the app's openFile does
  const ex = [...state.tabs.values()].find((t) => t.module === "plots" && t.path === rel);
  if (ex) { activateTab(ex.id); status(`${rel} is open already: its tab brought forward.`); return; }
  try {
    const tab = await plots.open({ name: rel.split("/").pop(), path: rel, absolute: files.absOf(rel), kind: id, dir: false });
    if (tab) sources.set(tab.id, { by: id === "plot_set" ? "drv" : "buf", path: rel });
    if (tab) { state.tabs.set(tab.id, tab); openTab(tab); renderPanels(); status(`${rel} opened in Plots: ${tab.plots && tab.plots.set ? tab.plots.set.members.length : 0} members, its files taken from the folder by the core.`); }
  } catch (e) { status(String(e), true); }
}

// ------------------------------------------------------------------ the sets the core finds in a folder (W2, step 3)
const proposals = { state: "none", list: [], left: [], classes: [], why: "" };
/** Each proposal's files as the folder's paths: the core names its set's members relative to the first file's folder
    (`plots_propose`, core 0.80.0), a file of several scans once. */
function withPaths(list, base) {
  const join = (b, rel) => { const out = b.split("/"); for (const part of rel.split("/")) { if (part === "..") out.pop(); else if (part && part !== ".") out.push(part); } return out.join("/"); };
  return list.map((p) => ({ ...p, paths: [...new Set(((p.set && p.set.members) || []).map((m) => (m.file.startsWith("/") ? m.file : join(base, m.file))))] }));
}
async function askProposals() {
  const gen = describing;
  Object.assign(proposals, { list: [], left: [], classes: [], why: "" });
  if (!files.root) proposals.state = "none";
  else if (!dispatcher.ready) proposals.state = "nocore";
  else if (!dispatcher.has("plots_propose")) proposals.state = "noproposer";
  else {
    proposals.state = "asking"; renderProposals();
    try {
      const paths = files.paths(), base = paths.length ? paths[0].replace(/\/[^/]*$/, "") : files.root;
      const r = await dispatcher.propose(paths);
      if (gen !== describing) return;
      proposals.list = withPaths(r.proposals || [], base);
      proposals.left = r.left_out || [];
      proposals.classes = r.classes || [];
      proposals.state = "done";
    } catch (e) { proposals.state = "error"; proposals.why = String(e); }
  }
  renderProposals();
}
function renderProposals() {
  const box = $("proposals"); if (!box) return;
  box.innerHTML = `<div class="sec-head"><span>Sets in this folder</span><span class="sub">${proposals.state === "done" ? `${proposals.list.length} found by the core` : ""}</span></div>`;
  const say = (t) => box.appendChild(el("div", "empty small", t));
  if (proposals.state === "none") return say("Open a folder: the core proposes the plot sets its files make (one scan per temperature, the holds, the whole run against time).");
  if (proposals.state === "nocore") return say("The core is not in this page: no proposals without it (W1).");
  if (proposals.state === "noproposer") return say(`The core in this page (drs-wasm ${esc(dispatcher.version || "")}) has no proposer yet: the sets a folder makes (one scan per temperature, the holds, the run against time) come with its series detector, W2 in core 0.80.0. Meanwhile choose the files, or a folder's "set".`);
  if (proposals.state === "asking") return say("The core reads the folder's files and proposes its sets…");
  if (proposals.state === "error") return say(`The core proposed nothing: ${esc(proposals.why)}`);
  if (!proposals.list.length) say("The core finds no set in these files.");
  proposals.list.forEach((p, k) => {
    const title = p.title || p.name || `Proposal ${k + 1}`;
    const n = p.paths ? p.paths.length : null;
    const live = [...lives.values()].some((l) => l.kind === p.kind && (l.class ?? null) === (p.class ?? null));
    const row = el("div", "wh-prop", `<button class="small">Open</button><span class="wh-name">${esc(title)}</span>${live ? `<span class="wh-kind" title="Opened, and grown by live mode with the new scans the core proposes for it">live</span>` : ""}${n != null ? `<span class="sub">${n}</span>` : ""}`);
    row.title = [`${title} (the core's ${p.kind}${p.class != null ? `, exposure class ${p.class + 1}` : ""})`, p.left_out && p.left_out.length ? `left out of it: ${p.left_out.map((x) => `${x[0]} (${x[1]})`).join("; ")}` : ""].filter(Boolean).join("\n");
    row.querySelector("button").onclick = () => openProposal(p, title);
    box.appendChild(row);
  });
  if (proposals.classes.length > 1) say(`${proposals.classes.length} exposure classes, told apart by the core by the median of their esd's.`);
  if (proposals.left.length) say(`Left out of every set: ${proposals.left.map((x) => (Array.isArray(x) ? `${esc(files.relOf(x[0]))} (${esc(x[1])})` : esc(String(x.name || x)))).join("; ")}.`);
}
/** A proposal opened: its set file as text when the core gives one (W10: held in the page under Plots/, never written,
    opened by the module's Open Plot Set), else its files as one set (`plots_open_from`, with the proposal's `rebin`
    when it says one, core 0.81.0). Recorded, so that live mode grows it with the new scans the core proposes for it. */
async function openProposal(p, title) {
  const text = p.text || (p.set && p.set.text) || null;
  let tab = null;
  if (text) {
    const rel = `Plots/${title.replace(/[\/:]/g, "-").slice(0, 80)}.drv`; files.put(rel, text); kinds.delete(rel);
    tab = await opened(() => openSetFile(rel, "plot_set"));
  } else if (p.paths && p.paths.length) {
    const rebin = p.rebin ?? (p.set && p.set.rebin);
    dispatcher.nextOpening(rebin != null ? { rebin } : null);
    pending = p.paths.slice();
    tab = await opened(() => run("plots.opencsv"));
    pending = null;
    const typed = ((p.set && p.set.members) || []).filter((m) => m.value != null).length;
    if (tab && typed) status(`${title}: opened from its ${p.paths.length} files. The ${typed} values its set file types (each with its source) are not carried: the core gives a proposal as a structure the page cannot open as a set file yet (W10); the members take what the core reads in their files.`, true);
  } else { status(`${title}: the core's proposal names no files to open.`, true); return; }
  if (!tab) return;
  sources.set(tab.id, { by: "proposal", kind: p.kind, class: p.class ?? null });
  lives.set(tab.id, { kind: p.kind, class: p.class ?? null, title, paths: new Set(p.paths || []) });
  renderProposals(); paintLive();
}
/** A folder shown anew: its files described, its sets proposed, the address's view brought back when its files are
    here, live mode offered when the browser gave the folder's handle. */
async function folderOpened() {
  lives.clear(); paintLive(); startLive();
  await describeAll(); await askProposals();
  await restoreView();
}

async function pickFolder() {
  try {
    if (window.showDirectoryPicker) {
      const h = await window.showDirectoryPicker({ id: "drs-web", mode: "read" });
      status(`Reading ${h.name}/ …`);
      const entries = await fromDirectoryHandle(h);
      folderHandle = h; files.replace(h.name, entries); setFiles(); folderOpened();
      status(`${h.name}/: ${entries.length} files listed, read here only when the core needs them; nothing leaves this machine.`);
    } else $("in-folder").click();
  } catch (e) { if (e && e.name !== "AbortError") status(`The folder could not be read: ${e}`, true); }
}
function takeList(input, asFolder) {
  const { folder, entries } = fromFileList(input.files);
  input.value = "";
  if (!entries.length) return;
  if (asFolder) { folderHandle = null; files.replace(folder || "folder", entries); setFiles(); folderOpened(); status(`${folder}/: ${entries.length} files listed; nothing leaves this machine.`); return; }
  return addAndOpen(entries);
}
async function dropped(dt) {
  const d = await fromDrop(dt);
  if (!d.entries.length) { status("Nothing readable was dropped (a file or a folder, from the Finder).", true); return; }
  if (d.folder) { folderHandle = d.handle || null; files.replace(d.folder, d.entries); setFiles(); folderOpened(); status(`${d.folder}/: ${d.entries.length} files listed; nothing leaves this machine. Choose files in the tree, or a folder's "set", to open them in Plots.`); return; }
  await addAndOpen(d.entries);
}
/** Loose files join the folder shown and open in Plots (a set file or a buffer as itself, the rest as one set). */
async function addAndOpen(entries) {
  const landed = files.add(entries); setFiles(); describeAll();
  if (proposals.state === "none") askProposals();
  await openPaths(landed);   // what was dropped, whether new to the folder or the same files again
}

// ------------------------------------------------------------------ live mode (step 4): the folder read again by its handle
/** The sets opened from a proposal, grown with the new scans the core proposes for them: tab id → { kind, class, title,
    paths }. */
const lives = new Map();
const live = { on: true, every: DEFAULT_INTERVAL, timer: null, busy: false, line: "", vanished: new Set() };
try { const v = JSON.parse(localStorage.getItem("drs-web.live") || "null"); if (v) { live.on = v.on !== false; if (INTERVALS.includes(v.every)) live.every = v.every; } } catch (e) { /* none kept */ }
function startLive() { clearTimeout(live.timer); live.vanished = new Set(); if (folderHandle) live.timer = setTimeout(tick, live.every * 1000); }
/** One reading: the folder walked again, what changed said, new files described, the live sets grown. Never removes a
    member: a scan that vanishes is said, and stays in its set. */
async function tick() {
  live.timer = null;
  if (!folderHandle || live.busy) return;
  if (live.on) {
    live.busy = true;
    try {
      const entries = await fromDirectoryHandle(folderHandle);
      const next = new Map(entries.map((e) => [e.rel, `${e.size}:${e.lastModified}`]));
      const d = diffListing(files.stamps(), next);
      if (!nothing(d)) {
        files.sync(entries); setFiles();
        for (const p of d.vanished) live.vanished.add(p);
        for (const p of d.added) live.vanished.delete(p);
        await describeAll();
        let grown = 0;
        if (d.added.length && dispatcher.has("plots_propose")) { await askProposals(); grown = await grow(); }
        live.line = liveLine(d, grown) || live.line;
        const inSets = d.vanished.filter((p) => [...lives.values()].some((l) => l.paths.has(files.absOf(p))));
        if (inSets.length) status(`Live: ${inSets.join(", ")} vanished from the folder; the set keeps ${inSets.length === 1 ? "it" : "them"} as read (nothing is removed from a set by itself).`, true);
      }
    } catch (e) { live.line = `the folder could not be read again: ${e}`; }
    live.busy = false;
  }
  paintLive();
  live.timer = setTimeout(tick, live.every * 1000);
}
/** The live sets grown: for each, the files the core now proposes for its kind and class that it lacks, added by the
    module's own Add Files… (`plots_add` → `plots_add_from`), the tab brought forward for it and the one in use back. */
async function grow() {
  let n = 0;
  for (const [id, l] of lives) {
    const tab = state.tabs.get(id); if (!tab) { lives.delete(id); continue; }
    const add = growth(proposals.list, l, l.paths); if (!add.length) continue;
    const back = state.activeTab;
    if (back !== id) activateTab(id);
    pending = add; await run("plots.add"); pending = null;
    for (const p of add) l.paths.add(p);
    n += add.length;
    if (back != null && back !== id && state.tabs.has(back)) activateTab(back);
  }
  return n;
}
function paintLive() {
  const box = $("live"); if (!box) return;
  box.hidden = !folderHandle;
  if (!folderHandle) return;
  const sel = $("live-every");
  if (!sel.options.length) sel.innerHTML = INTERVALS.map((s) => `<option value="${s}">${s} s</option>`).join("");
  sel.value = String(live.every);
  $("live-on").checked = live.on;
  const sets = lives.size ? `${lives.size} live set${lives.size === 1 ? "" : "s"}` : "no live set (open a proposal to grow one)";
  $("live-line").textContent = `${live.on ? `Reading the folder every ${live.every} s` : "Paused"}; ${sets}.${live.vanished.size ? ` Vanished: ${[...live.vanished].join(", ")}.` : ""}${live.line ? ` ${live.line}` : ""}`;
}
function installLive() {
  const keepLive = () => { try { localStorage.setItem("drs-web.live", JSON.stringify({ on: live.on, every: live.every })); } catch (e) { /* */ } };
  $("live-on").onchange = (e) => { live.on = e.target.checked; keepLive(); paintLive(); if (live.on) { clearTimeout(live.timer); tick(); } };
  $("live-every").onchange = (e) => { live.every = Number(e.target.value) || DEFAULT_INTERVAL; keepLive(); startLive(); paintLive(); };
  $("live-now").onclick = () => { clearTimeout(live.timer); tick(); };
}

// ------------------------------------------------------------------ the view in the address (step 4)
/** The active tab's view as the hash keeps it: where its set came from, the view, the colour scale, the map's zoom,
    the clip of the view shown, the region chosen. Read from the tab's state as the module keeps it; written with
    replaceState (no history entry per change). */
function viewNow() {
  const t = state.activeTab != null ? state.tabs.get(state.activeTab) : null;
  const src = t && sources.get(t.id); if (!src || !t.plots) return null;
  const ps = t.plots, key = ps.view === "difference" ? "difference" : "map";
  const clip = ps.clip && ps.clip[key];
  const region = ps.selregion != null && ps.regions && ps.regions[ps.selregion] ? ps.regions[ps.selregion].name : null;
  return { folder: files.name, set: src, view: ps.view, scale: ps.mapScale, zoom: ps.mapZoom || null, clip: clip ? { from: clip.from, to: clip.to } : null, region };
}
let lastHash = null;
let restoring = false;
function keepView() {
  if (restoring) return;   // the address is the view being brought back until it is drawn
  const h = encode(viewNow());
  if (h && h !== lastHash && h !== location.hash) { history.replaceState(null, "", h); }
  lastHash = h || lastHash;
}
/** The address's view brought back once its files are here: the set opened as it was (its files, its set file or
    buffer, or the same proposal), then its view, scale, zoom, clip and region set in the tab and drawn. */
let wanted = decode(location.hash);
async function restoreView() {
  const w = wanted; if (!w || !files.root) return;
  const miss = missingFor(w.set, (rel) => files.has(files.absOf(rel)));
  if (miss.length) { status(`The address keeps a view of ${w.folder ? `${w.folder}/` : "a folder"}'s files; ${miss.length} of them are not in the folder open (${miss.slice(0, 3).join(", ")}${miss.length > 3 ? " …" : ""}): open that folder to bring the view back.`, true); return; }
  wanted = null;
  let tab = null;
  if (w.set.by === "files") { pending = w.set.files.map(files.absOf); tab = await opened(() => run("plots.opencsv")); pending = null; if (tab) sources.set(tab.id, w.set); }
  else if (w.set.by === "drv" || w.set.by === "buf") tab = await opened(() => openSetFile(w.set.path, w.set.by === "drv" ? "plot_set" : "buffer_file"));
  else if (w.set.by === "proposal") {
    const p = proposals.list.find((q) => q.kind === w.set.kind && (q.class ?? null) === (w.set.class ?? null));
    if (!p) { status(`The address names the core's proposal "${w.set.kind}", which the core does not make of this folder.`, true); return; }
    await openProposal(p, p.title); tab = state.activeTab != null ? state.tabs.get(state.activeTab) : null;
  }
  if (!tab || !tab.plots) return;
  const ps = tab.plots;
  restoring = true;
  try {
    if (w.view) ps.view = w.view;
    if (w.scale) { ps.mapScale = w.scale; ps.map = null; }
    if (w.clip) { const key = ps.view === "difference" ? "difference" : "map"; ps.clip = { ...(ps.clip || {}), [key]: { from: w.clip.from, to: w.clip.to, scale: key === "difference" && ps.mapScale === "log" ? "sqrt" : ps.mapScale } }; }
    bus.emit("theme");   // the module draws the active tab again (its own listener), with the state set here
    renderPanels();
    // the zoom once the core's map is in the tab: its first drawing fits the whole map, which clears a zoom set before
    if (w.zoom && (ps.view === "map" || ps.view === "difference")) {
      for (let k = 0; k < 100 && !ps.map; k++) await new Promise((r) => setTimeout(r, 100));
      await new Promise((r) => setTimeout(r, 200));
      ps.mapZoom = w.zoom; bus.emit("theme");
    }
    if (w.region) { const row = [...document.querySelectorAll("#members .child")].find((r) => r.firstElementChild && r.firstElementChild.textContent === w.region); if (row) row.click(); }
  } finally { restoring = false; }
  status(`The view in the address brought back: ${tab.title}, ${ps.view}${w.zoom ? ", zoomed" : ""}${w.clip ? ", its colour range clipped" : ""}${w.region ? `, region ${w.region}` : ""}.`);
}

// ------------------------------------------------------------------ offline (step 4): the service worker
const offline = { version: null, line: "" };
async function startOffline() {
  try { offline.version = await (await fetch("version.json", { cache: "no-cache" })).json(); } catch (e) { offline.version = null; }
  if (!("serviceWorker" in navigator)) { offline.line = "Offline: not in this browser (no service worker)."; return paintOffline(); }
  try {
    await navigator.serviceWorker.register("sw.js");
    navigator.serviceWorker.addEventListener("message", (e) => { if (e.data && e.data.kept) { const k = e.data.kept; offline.line = `Offline: ${k.files} files, ${fmtSize(k.bytes)} kept by this browser (web ${k.version}${core.version ? `, core ${core.version}` : ""})${k.guide ? `, and the guide's ${k.guide} files` : ""}: the page and its manual open again without the network.`; paintOffline(); } });
    const askKept = () => { if (navigator.serviceWorker.controller) navigator.serviceWorker.controller.postMessage({ ask: "kept" }); };   // a message to the worker, not a dialog of the app's (named apart from the dialog helpers)
    navigator.serviceWorker.addEventListener("controllerchange", askKept);
    offline.line = navigator.serviceWorker.controller ? "Offline: counting what this browser keeps…" : "Offline: being kept now (the first visit installs it; it is ready from the next).";
    paintOffline(); askKept();
  } catch (e) { offline.line = `Offline: not kept (${e.message || e}).`; paintOffline(); }
}
function paintOffline() { const n = $("offline-line"); if (n) n.textContent = offline.line; }

// ------------------------------------------------------------------ the panels the module draws into
function renderPanels() {
  const sc = $("insp-scroll"), top = sc.scrollTop;
  const m = $("members"); m.innerHTML = ""; plots.sidebar(m);
  const i = $("member"); i.innerHTML = ""; plots.inspector(i);
  sc.scrollTop = top;
  renderMenu();
}
window.__renderPanels = renderPanels;
bus.on("tab", (tabId) => { state.activeTab = tabId; state.activeModule = "plots"; renderPanels(); });
bus.on("tab-closed", (tabId) => {
  const t = state.tabs.get(tabId); if (!t) return;
  if (plots.close) plots.close(t);
  state.tabs.delete(tabId);
  if (state.activeTab === tabId) state.activeTab = null;
  renderPanels();
});

/** The guide: the page's own site (web/site/, published at /guide/ beside the page), its manual, recipes and quick
    reference. Absolute, so that a local checkout (no guide served) links to the published one. */
const GUIDE = "https://jmigartua.github.io/diff-ref-struct-web/guide/";
/** The public repository (its issues are the place for requests) and the address for those without an account. */
const REPO = "jmigartua/diff-ref-struct-web";
const CONTACT = "josu.igartua@ehu.eus";
/** The manual's section for each view (the stable anchors `<a id="help-…">` of web/manual/user-manual.md). */
const VIEW_HELP = { overlay: "help-overlay", stack: "help-stack", projection: "help-projection", map: "help-map", difference: "help-difference", surface: "help-surface", series: "help-series" };
/** The version words: diff-ref-struct Web 1.0.0 (version.json, stamped by web/tools/stage.mjs from the releases
    record), the core's, the build; "dev" in a local checkout. */
function versionWords() {
  const v = offline.version || {};
  return { web: v.version || (v.web && v.web !== "dev" ? v.web : "dev"), build: v.build || v.web || null, core: core.version || v.core || null };
}
function openGuide(path) { const w = window.open(GUIDE + path, "_blank"); if (w) w.opener = null; return !!w; }
/** The manual at the section of the view shown (the active tab's view), or at its start. */
function manualHere() {
  const t = state.activeTab != null ? state.tabs.get(state.activeTab) : null;
  const id = t && t.plots ? VIEW_HELP[t.plots.view] : null;
  const ok = openGuide(`manual/index.html${id ? "#" + id : ""}`);
  status(ok ? `The manual opened in a new tab${id ? `, at the ${t.plots.view} view` : ""}.` : `The manual: ${GUIDE}manual/ (the browser kept the tab from opening).`, !ok);
}
/** What a request carries without the user typing it: the version, the browser, the view's address (no data). */
function requestFields() {
  const v = versionWords();
  return { version: `diff-ref-struct Web ${v.web}${v.build ? ` (build ${v.build})` : ""}, core ${v.core || "?"}`, browser: navigator.userAgent, address: location.href };
}
function issueUrl(template) {
  const q = new URLSearchParams({ template: `${template}.yml`, ...requestFields() });
  return `https://github.com/${REPO}/issues/new?${q}`;
}
function mailUrl(kind) {
  const f = requestFields();
  const body = `${kind === "bug" ? "What happened, and what you expected:" : "What you would like to do, with which data, and what you do today instead:"}\n\n\n---\nVersion: ${f.version}\nBrowser: ${f.browser}\nThe view's address: ${f.address}\n(Do not attach unpublished data: a description and the address are enough.)`;
  return `mailto:${CONTACT}?subject=${encodeURIComponent(`diff-ref-struct Web: ${kind === "bug" ? "a problem" : "a feature"}`)}&body=${encodeURIComponent(body)}`;
}
/** One line to cite the view shown: the tool's and the core's versions, the date, the address. */
async function citeView() {
  const v = versionWords(), day = new Date().toISOString().slice(0, 10);
  const line = `J. M. Igartua, diff-ref-struct Web ${v.web} (core ${v.core || "?"}), https://jmigartua.github.io/diff-ref-struct-web/; view of ${day}: ${location.href}`;
  try { await dispatcher.invoke("clipboard_text", { text: line }); status(`Copied: ${line}`); } catch (e) { status(`The citation (copy it by hand): ${line}`); }
}
const HELP_MENU = [
  ["The guide", () => openGuide("index.html")],
  ["Manual for this view", () => manualHere()],
  ["Manual of this version (PDF)", () => { const v = versionWords().web; openGuide(v && v !== "dev" ? `downloads/user-manual-${v}.pdf` : "downloads/user-manual.pdf"); }],
  ["What's new in this version", () => openGuide("releases/index.html")],
  ["-"],
  ["Report a problem…", () => { const w = window.open(issueUrl("bug"), "_blank"); if (w) w.opener = null; status("The bug form opened in a new tab, with the version, the browser and this view's address filled in (no data)."); }],
  ["Ask for a feature…", () => { const w = window.open(issueUrl("feature"), "_blank"); if (w) w.opener = null; status("The feature form opened in a new tab, with the version filled in."); }],
  ["Write instead (e-mail)…", () => { location.href = mailUrl("bug"); }],
  ["-"],
  ["Cite this view", () => citeView()],
  ["About, licence and credits", () => openGuide("about/index.html")],
];
function renderHelp() {
  const box = $("menu-help"); box.innerHTML = "";
  for (const [label, act] of HELP_MENU) {
    if (label === "-") { box.appendChild(el("div", "wh-msep")); continue; }
    const b = el("button", "wh-mitem", `<span class="wh-check"></span>${esc(label)}`);
    b.setAttribute("role", "menuitem");
    b.onclick = () => { closeHelp(); act(); };
    box.appendChild(b);
  }
}
function openHelp() { const m = $("menu-help"); m.hidden = false; $("btn-help").setAttribute("aria-expanded", "true"); renderHelp(); const f = m.querySelector("button"); if (f) f.focus(); }
function closeHelp() { $("menu-help").hidden = true; $("btn-help").setAttribute("aria-expanded", "false"); }

/** Try the example: a heating run computed by the core (web/examples/tools/make_example.py), the page's own files
    beside it (examples/ramp/), opened as a folder would be, its proposal opened on the Map. Fetched from the page's
    address only; nothing of the user's is read. */
async function tryExample() {
  try {
    status("The example: reading its 21 scans from the page…");
    const m = await (await fetch("examples/ramp/manifest.json", { cache: "no-cache" })).json();
    const t0 = Date.parse("2026-10-26T09:00:00Z");
    const entries = m.files.map((f, k) => ({
      rel: f.name, size: f.size, lastModified: t0 + k * 1800 * 1000, handle: null,
      getFile: async () => { const r = await fetch(`examples/ramp/${f.name}`); if (!r.ok) throw new Error(`HTTP ${r.status}`); const text = await r.text(); return { text: async () => text }; },
    }));
    folderHandle = null; files.replace("example-ramp", entries); setFiles();
    await folderOpened();
    const p = proposals.list.find((q) => q.kind === "by_temperature");
    if (p) {
      await openProposal(p, p.title);
      const tab = state.activeTab != null ? state.tabs.get(state.activeTab) : null;
      if (tab && tab.plots) { tab.plots.view = "map"; tab.plots.mapScale = "sqrt"; tab.plots.map = null; bus.emit("theme"); renderPanels(); }
    }
    status(`The example: ${m.title}. ${m.note} Computed by the core (${m.generated_by || "drs"}); the guide's Start page walks through it.`);
  } catch (e) { status(`The example could not be read: ${e}`, true); }
}

function paintVersion() {
  const n = $("start-version"); if (!n) return;
  const v = offline.version;
  const w = versionWords();
  n.innerHTML = `diff-ref-struct Web ${esc(w.web)}${w.build && w.build !== w.web ? ` (build ${esc(w.build)})` : ""}${v && v.built ? `, built ${esc(v.built)}` : ""}${core.version ? ` · core ${esc(core.version)} in WebAssembly` : ""} · no account, no tracking, no upload · <a href="${GUIDE}releases/index.html" target="_blank" rel="noopener">what's new</a> · open source (Apache-2.0)`;
}
function renderStart(box) {
  box.innerHTML = `<div class="wh-start">
    <h1>diff-ref-struct Web</h1>
    <p class="lede">Several data files at once as a plot set, the map with its colour range, cursors, zoom and regions: the app's Plots module in the browser, the core doing every number in WebAssembly.</p>
    <p><button class="primary" data-act="folder">Open Folder…</button> <button data-act="files">Open Files…</button> <button data-act="example" title="A heating run computed by the core (21 scans, 300 to 1300 K, a second phase appearing above 1000 K), opened as if its folder had been dropped">Try the example</button></p>
    <p class="sub">Or drop a folder or files anywhere on this page. They are read in this browser and never leave this machine: nothing is uploaded.</p>
    <p class="sub wh-core-note"${core.ready ? " hidden" : ""}>The core is not in this page yet (request W1, its WebAssembly build): your files open and are listed; every view waits for the core and says so.</p>
    <p class="sub">In Chromium (Chrome, Edge, Arc), a folder opened with Open Folder… is read again every few seconds: new scans appear, and a set opened from the core's proposals grows with them. The page keeps itself for offline use and can be installed as an app; the view shown is kept in the address.</p>
    <p class="sub">New here? <a href="${GUIDE}start/index.html" target="_blank" rel="noopener">Start in five minutes</a> with public example data, then <a href="${GUIDE}recipes/index.html" target="_blank" rel="noopener">the recipes</a> (a ramp, a hold, two patterns compared, a run followed live) and <a href="${GUIDE}manual/index.html" target="_blank" rel="noopener">the manual</a> (also as PDF).</p>
    <p class="sub wh-version" id="start-version"></p>
  </div>`;
  paintVersion();
  box.querySelector('[data-act="folder"]').onclick = pickFolder;
  box.querySelector('[data-act="files"]').onclick = () => $("in-files").click();
  box.querySelector('[data-act="example"]').onclick = tryExample;
}

// ------------------------------------------------------------------ the Plots menu (the module's own, as the app's menu bar lists it)
function renderMenu() {
  const box = $("menu-plots"); if (box.hidden) return;
  box.innerHTML = "";
  for (const id of plots.menu) {
    if (id === "-") { box.appendChild(el("div", "wh-msep")); continue; }
    const c = cmd(id); if (!c) continue;
    const b = el("button", "wh-mitem" + (c.check && c.check() ? " on" : ""), `<span class="wh-check">${c.check && c.check() ? "✓" : ""}</span>${esc(c.title)}`);
    b.disabled = !enabled(c); b.title = c.disabled && !enabled(c) ? c.disabled : c.what || "";
    b.setAttribute("role", "menuitem");
    b.onclick = () => { closeMenu(); run(id); };
    box.appendChild(b);
  }
}
function openMenu() { const m = $("menu-plots"); m.hidden = false; $("btn-plots").setAttribute("aria-expanded", "true"); renderMenu(); const f = m.querySelector("button:not(:disabled)"); if (f) f.focus(); }
function closeMenu() { $("menu-plots").hidden = true; $("btn-plots").setAttribute("aria-expanded", "false"); }

// ------------------------------------------------------------------ the inspector's width and rail
const shape = kept(localStorage);
function applyShape() {
  const insp = $("insp"), w = shape.rail ? RAIL : clampWidth(shape.width, window.innerWidth);
  document.documentElement.style.setProperty("--insp-w", w + "px");
  insp.classList.toggle("rail", shape.rail);
  const b = $("btn-rail"); b.textContent = shape.rail ? ">" : "<";
  b.title = shape.rail ? "Show the inspector (the folder's tree and the members)" : "Shrink the inspector to a rail";
  b.setAttribute("aria-label", b.title); b.setAttribute("aria-expanded", String(!shape.rail));
  $("sep").setAttribute("aria-valuenow", String(w));
}
function installSeparator() {
  const sep = $("sep");
  sep.addEventListener("pointerdown", (e) => {
    if (e.button !== 0) return;
    if (shape.rail) { shape.rail = false; }
    const w0 = clampWidth(shape.width, window.innerWidth), x0 = e.clientX;
    sep.setPointerCapture(e.pointerId); document.body.classList.add("wh-dragging");
    const move = (ev) => { shape.width = dragged(w0, x0, ev.clientX, window.innerWidth); applyShape(); };
    const up = () => { sep.removeEventListener("pointermove", move); sep.removeEventListener("pointerup", up); document.body.classList.remove("wh-dragging"); keep(localStorage, shape); };
    sep.addEventListener("pointermove", move); sep.addEventListener("pointerup", up);
    e.preventDefault();
  });
  sep.addEventListener("dblclick", () => { shape.rail = !shape.rail; applyShape(); keep(localStorage, shape); });
  sep.addEventListener("keydown", (e) => { const w = keyed(clampWidth(shape.width, window.innerWidth), e.key, e.shiftKey, window.innerWidth); if (w !== shape.width && !shape.rail) { shape.width = w; applyShape(); keep(localStorage, shape); e.preventDefault(); } });
  $("btn-rail").onclick = () => { shape.rail = !shape.rail; applyShape(); keep(localStorage, shape); };
  window.addEventListener("resize", applyShape);
}

// ------------------------------------------------------------------ the theme (the kit's variables, style.css)
const THEMES = ["system", "light", "dark"];
function applyTheme(t) {
  if (t === "system") document.documentElement.removeAttribute("data-theme"); else document.documentElement.setAttribute("data-theme", t);
  state.settings.theme = t; $("btn-theme").textContent = { system: "◐ auto", light: "☀ light", dark: "☾ dark" }[t];
  $("btn-theme").title = `Appearance: ${t === "system" ? "the system's" : t}. Click for the next.`;
  bus.emit("theme");
}

// ------------------------------------------------------------------ start
async function boot() {
  let theme = "system"; try { theme = localStorage.getItem("drs-web.theme") || "system"; } catch (e) { /* */ }
  applyTheme(THEMES.includes(theme) ? theme : "system");
  $("btn-theme").onclick = () => { const t = THEMES[(THEMES.indexOf(state.settings.theme) + 1) % THEMES.length]; applyTheme(t); try { localStorage.setItem("drs-web.theme", t); } catch (e) { /* */ } };
  state.activeModule = "plots";
  register(plots.commands);
  if (plots.init) await plots.init();
  initDock($("dock"), {
    file: (tabId, box) => { const t = state.tabs.get(tabId); return t ? plots.view(t, box) : null; },
    start: renderStart,
    tabModule: (tabId) => (state.tabs.get(tabId) ? { id: plots.id, name: plots.name, icon: plots.icon } : null),
    help: () => manualHere(),
  }, () => {});
  ensureStart(); relayout();   // the dock laid out once its first panel is in (the app does so after its default layout)
  applyShape(); installSeparator();
  $("btn-folder").onclick = pickFolder;
  $("btn-files").onclick = () => $("in-files").click();
  $("btn-openset").onclick = () => openPaths(tree.sel);
  $("in-folder").onchange = (e) => takeList(e.target, true);
  $("in-files").onchange = (e) => takeList(e.target, false);
  $("btn-plots").onclick = (e) => { e.stopPropagation(); closeHelp(); if ($("menu-plots").hidden) openMenu(); else closeMenu(); };
  $("btn-help").onclick = (e) => { e.stopPropagation(); closeMenu(); if ($("menu-help").hidden) openHelp(); else closeHelp(); };
  document.addEventListener("click", (e) => { if (!e.target.closest("#menu-plots")) closeMenu(); if (!e.target.closest("#menu-help")) closeHelp(); });
  document.addEventListener("keydown", (e) => { if (e.key === "Escape") { closeMenu(); closeHelp(); } });
  // a drop anywhere on the page: a folder becomes the folder shown, files join it and open in Plots
  const hasFiles = (e) => e.dataTransfer && [...e.dataTransfer.types].includes("Files");
  document.addEventListener("dragover", (e) => { if (!hasFiles(e)) return; e.preventDefault(); e.dataTransfer.dropEffect = "copy"; document.body.classList.add("wh-dropping"); });
  document.addEventListener("dragleave", (e) => { if (e.relatedTarget == null) document.body.classList.remove("wh-dropping"); });
  document.addEventListener("drop", (e) => { if (!hasFiles(e)) return; e.preventDefault(); document.body.classList.remove("wh-dropping"); dropped(e.dataTransfer).catch((err) => status(`The drop could not be read: ${err}`, true)); });
  // the tree by folders or by the core's kinds (step 3), the choice kept
  const paintBy = () => document.querySelectorAll("#tree-by button").forEach((x) => { x.classList.toggle("on", x.dataset.by === tree.by); x.setAttribute("aria-pressed", String(x.dataset.by === tree.by)); });
  document.querySelectorAll("#tree-by button").forEach((x) => { x.onclick = () => { tree.by = x.dataset.by; try { localStorage.setItem("drs-web.tree.by", tree.by); } catch (e) { /* */ } paintBy(); rebuild(); }; });
  paintBy();
  installLive();
  setInterval(keepView, 1000);
  window.addEventListener("hashchange", () => { if (location.hash !== lastHash) { wanted = decode(location.hash); restoreView(); } });
  window.__drsWebState = () => ({ live: { on: live.on, every: live.every, sets: lives.size, line: live.line, vanished: [...live.vanished] }, hash: location.hash, offline: offline.line, core: core.ready, coreChecked: core.checked, version: core.version, files: tree.root.files, folder: files.name, selected: tree.sel.length, tabs: state.tabs.size, rail: shape.rail, width: shape.width, handle: !!folderHandle, by: tree.by, described: kinds.size, kinds: [...kinds.values()].reduce((o, k) => ((o[k.kind === "data" || k.kind === "project" ? k.id : k.kind] = (o[k.kind === "data" || k.kind === "project" ? k.id : k.kind] || 0) + 1), o), {}), proposals: proposals.state, calls: { ...dispatcher.calls } });
  renderTree(); renderPanels(); renderProposals(); paintCore();
  await startCore();
  startOffline().then(paintVersion);
  paintVersion();
  status(core.ready ? `The core ${core.version || ""} runs in this page. Open a folder or drop files.` : "The page is ready; the core waits for its WebAssembly build (W1). Open a folder or drop files: they are listed here and never uploaded.");
}
boot().catch((e) => { console.error("diff-ref-struct Web did not start", e); status(`The page did not start: ${e}`, true); });
