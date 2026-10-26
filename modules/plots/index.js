/* diff-ref-struct · the Plots module (R32, L154; plots/DESIGN.md, accepted by Josu, L169). Several data files as one
   plot set, the counterpart of WinPLOTR's buffers, shown as an overlay, a stack, a map, a projection with hidden parts,
   the map as a surface in 3D (surface3d.js, step 15) and a series of window quantities against the temperature
   (DESIGN.md §10, steps 1 and 2); and FullProf's `.seq`, a
   sequential refinement's parameters against its external parameter (PL11).
   - The set is the core's (`plots::open`, `open_buf`), kept by the module's Rust side (plots.rs) while it is shown:
     the core gives no plot-set file yet (PL12), so it lives in memory and Save stays off.
   - Each member's temperature is the core's, with where it was found; a file that states none is drawn apart until
     the user types one, which keeps its source (Josu, L169). With no monitor the counts are raw, and the views say
     "no monitor: raw counts" (L169).
   - Every number is the core's (L139): the displaced curves and their baselines (`plots::displace`), the map's rows,
     edges, values and colour range (`plots::map`), the sums, averages and differences (`plots::derive`, written to
     Plots/), the projection's visible pieces and oblique axis (`plots::project`, PL17), the windows' integrals,
     centroids and peaks with σ (`plots::windows`, PL5), a `.seq`'s blocks and columns (`seq::read`, `Seq::column`,
     PL11). The window turns them into pixels and colours: a member's colour from its place in the set's order, its
     opacity (dimming), the plot's scale.
   - Normalise… and Convert X… are views (PL6, PL10): the Rust side keeps the rule with the set and every view draws
     the members as the core gives them, normalised (`plots::normalise`) and in the x quantity chosen (`plots::x_view`),
     each view saying so; the files and Data's operations untouched. Export ▸ Plot Set as CSV… writes the core's CSV
     (`plots::csv`, `csv_long`, in 2θ, d or Q by `plots::converted`, PL7); a CSV of several columns opens as a set (PL8).
   - Step 16 (I08, R33): the map's colour range clipped on its colour bar (two marks dragged, or typed), per view and on
     the scale shown, the same pair the surface's height; and a row and a column cursor, the row's member pattern below
     the map and the column's values against the temperature beside it, every number the core's. The clip lives in the
     tab until the set file keeps it (PL37). */
// W6 (step 21): the core is asked through core.js's host (the app's Tauri, or the page's own); `hasCore` says one
// answers, `IN_APP` is for what is the Mac app's alone (the app's own files for the Dev page); the open and save
// panels are the host's (`pickFiles`, `pickSavePath`)
import { $, el, esc, state, bus, activeTab, status, invoke, IN_APP, hasCore, pickFiles, pickSavePath, askForm, askText, pickOne, fmt, writable } from "../../core.js";
import { createPlot } from "../../plot.js";
import { barValue } from "../../mapclip.js";   // the shell's place on the colour bar, for the self-check's drag
// step 18 (I08 step 5, L409): the Map view's viewer, factored out so that the core's page hosts the same map; the
// module is its app host (the core's commands behind it) and draws nothing of the map itself
import { createMapViewer, clipOf, diffScaleOf, zWordOf, kindIcon, kindFamily, sortMembers, exposureClasses, membersInRange } from "../../viewer/map.js";
import { VIEWER_SOURCES, bundleViewer } from "../../viewer/bundle.js";
import { openTab, activateTab } from "../../layout.js";
import { createTable } from "../../table.js";
import { inLightPalette, sideBySide } from "../../svgexport.js";
import { run as runCommand } from "../../registry.js";

// ------------------------------------------------------------------ the set's colours (display only)
// A member's colour from its place in the set's order (the core's, by the variable): dark violet for the first (the
// coldest) to orange for the last, magma's middle part (the mock-ups' ramp, plots/mockups); grey for a member drawn
// apart (no value of the variable).
// In dark mode the ramp is magma's lighter part, so that the coldest members stay visible on the dark background
// (APP_CONVENTIONS §10: every colour has its light and its dark value); a figure exported is drawn in the light one.
const MAGMA = { light: ["#2c115f", "#51127c", "#721f81", "#932b80", "#b73779", "#d8456c", "#f1605d", "#fc8961"], dark: ["#8c2981", "#a3307e", "#bc3978", "#d3436e", "#e95462", "#f8765c", "#fe9f6d", "#fec98d"] };
const isDark = () => { const t = document.documentElement.getAttribute("data-theme"); return t === "dark" || (t !== "light" && window.matchMedia("(prefers-color-scheme: dark)").matches); };
function rampColour(t) {
  const stops = MAGMA[isDark() ? "dark" : "light"];
  const k = Math.max(0, Math.min(1, t)) * (stops.length - 1);
  const i = Math.min(stops.length - 2, Math.floor(k)), f = k - i;
  const a = stops[i].match(/\w\w/g).map((h) => parseInt(h, 16)), b = stops[i + 1].match(/\w\w/g).map((h) => parseInt(h, 16));
  return `rgb(${a.map((v, q) => Math.round(v + (b[q] - v) * f)).join(",")})`;
}
/** A window's colour and marker by its place in the set (ThesisFigures C7b, C7c), Okabe & Ito's colours and the shapes
    in Publish's order (publish/src/patterns.rs WINDOW_COLOURS, WINDOW_SHAPES), so that the screen and the paper agree. */
const WINDOW_COLOURS = ["#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442"];
const WINDOW_SHAPES = ["circle", "square", "diamond", "triangle"];
const DIMMED = 0.22;   // a dimmed member's opacity (Dim Others, Solo Member)
const VIEWS = [["overlay", "Overlay"], ["stack", "Stack"], ["projection", "Projection"], ["map", "Map"], ["difference", "Difference Map"], ["surface", "Surface"], ["series", "Series"]];
/** The series view's quantities of a window (PL5), each with where the core gives it and its σ when it gives one. */
const QUANTITIES = [
  ["integral", "integral above the background", (w) => w.integral, (w) => w.sigma, "counts·°"],
  ["centroid", "centroid", (w) => w.centroid, (w) => w.centroid_sigma, "°"],
  ["peak_x", "strongest peak: position", (w) => (w.peak ? w.peak.x : null), () => null, "°"],
  ["peak_height", "strongest peak: height", (w) => (w.peak ? w.peak.height : null), () => null, "counts"],
  ["peak_fwhm", "strongest peak: width", (w) => (w.peak ? w.peak.fwhm : null), () => null, "°"],
];

// ------------------------------------------------------------------ the tab's state
const tabOf = (id) => "plots:" + id;
const activePlots = () => { const t = activeTab(); return t && t.module === "plots" ? t : null; };
/** The tab's own state: the set as the core gave it, the view, the members chosen (in the order clicked), each
    member's opacity, the displacement rule, the y scale, the reference of the difference map, the core's answers. */
const stateOf = (tab) => (tab.plots ||= { set: null, view: "overlay", sel: [], opacity: {}, rule: { kind: "step", dx: 0, dy: null }, prule: { kind: "step", dx: null, dy: null }, scale: "linear", mapScale: "log", reference: null, disp: null, map: null, proj: null, wins: [], wbg: "auto", at: [], quantity: "integral", winres: null, seq: null, param: null, xby: "parameter", seqcol: null, linked: true, zoomwin: null, diff: null, nc: null, xc: null, xg: {}, look: false, pbar: null, sbar: null, figmm: 160, figname: "UPV/EHU TFG, full width", figfont: { label_pt: 11, tick_pt: 10 }, surfOf: "map", clip: { map: null, difference: null }, cursors: true, pin: null, cell: null, regions: [], selregion: null, rangeFromView: false, clipBefore: { map: null, difference: null }, mapZoom: null, regionsRefused: null, yInterval: "shared", sortBy: "set", extent: null });
/** A member's pattern from a map of patterns by member index (step 22), or null. */
const patternOf = (pats, m) => (pats ? pats.get(m.index) || null : null);
/** The word of the x quantity the map is shown in (PL37's `x_quantity`): the view's, else the files' own (none for TOF
    or energy, or members in different quantities). */
function xWordOf(ps) {
  if (ps.set && ps.set.x_quantity) return ps.set.x_quantity;
  const ax = members(ps)[0] ? members(ps)[0].axis : "";
  const w = { "2θ (°)": "2theta", "d (Å)": "d", "Q (Å⁻¹)": "q" }[ax] || null;
  return w && members(ps).every((m) => m.axis === ax) ? w : null;
}
const members = (ps) => (ps.set ? ps.set.members : []);
const memberByKey = (ps, key) => members(ps).find((m) => m.key === key) || null;
const placed = (ps) => members(ps).filter((m) => m.variable);
// in a figure's look the ramp stops short of its palest end, so that the warmest curve and its label read on white
const colourOf = (ps, m) => { if (!m.variable) return "--muted"; const p = placed(ps); const k = p.indexOf(m); return rampColour((p.length > 1 ? k / (p.length - 1) : 0) * (ps.look ? 0.8 : 1)); };
const opacityOf = (ps, m) => (ps.opacity[m.key] == null ? 1 : ps.opacity[m.key]);
/** The curve under the pointer drawn a little thicker (DESIGN.md §4: hovering a curve highlights it; the readout names
    it), never in a figure. */
const widthOf = (ps, m) => (m.key === ps.hoverKey && !ps.look ? 2 : 1);
/** A member's value of the variable as text, with its source, or what to do when it has none (L169). */
function variableText(ps, m, long = false) {
  const s = ps.set;
  if (!m.variable) return long ? `no ${s.variable_name}: the file states none; drawn apart until one is typed (Set ${s.variable_name === "T" ? "Temperature" : "Value"}…)` : `no ${s.variable_name}`;
  const v = `${s.variable_name} = ${fmt(m.variable.value)} ${s.variable_unit}`;
  return long ? `${v}, from ${m.variable.source}` : v;
}
/** "no monitor: raw counts" (L169) when the members state no monitor; how many when some do. */
function monitorNote(ps) {
  const ms = members(ps); if (!ms.length || (ps.set && ps.set.normalisation)) return "";
  const without = ms.filter((m) => m.monitor == null).length;
  if (without === ms.length) return "no monitor: raw counts";
  return without ? `${without} of ${ms.length} members have no monitor: raw counts` : "";
}

// ------------------------------------------------------------------ the view's normalisation and x (PL6, PL10)
/** The x quantity shown (PL10): the core's label ("d (Å)"), and its name without the unit. */
const xLabelOf = (ps) => (ps.set && ps.set.x_label) || "2θ (°)";
const xNameOf = (ps) => xLabelOf(ps).replace(/ \([^)]*\)$/, "");
/** Whether the view is in the files' own x, where the windows (ranges of the files' x) are drawn as bands. */
const inFilesX = (ps) => !ps.set || !ps.set.x_quantity || members(ps).every((m) => m.axis === ps.set.x_label);
/** What the views' y is: the counts, or the counts normalised (PL6). */
const yWord = (ps) => (ps.set && ps.set.normalisation ? "normalised counts" : "counts");
/** The normalisations the core makes (PL6), in the order offered. */
const NORMS = [
  ["none", "none: the counts as read"],
  ["monitor", "the monitor: each member's counts × the value / its file's monitor"],
  ["time", "the counting time: each member's counts × the value / its file's counting time"],
  ["maximum", "the maximum: each member's largest y made the value (in a range, when one is given)"],
  ["area", "the same area: each member's area over a range made the value"],
  ["peak", "a peak's integral: each member's integral over a range, above its background, made the value"],
  ["reference", "a reference member: each member's area over a range made the chosen member's"],
];
/** The x quantities the core shows (PL10; WinPLOTR's 2THETA, D, Q, STL, S). */
const XQS = [["", "the files' own x"], ["2theta", "2θ (°)"], ["d", "d (Å)"], ["q", "Q = 4π sin θ/λ (Å⁻¹)"], ["stl", "sin θ/λ (Å⁻¹)"], ["s", "s = 1/d = 2 sin θ/λ (Å⁻¹)"]];
/** The normalisation shown, in words, from the window's own choice as the core echoes it. */
function normText(ps) {
  const c = ps.set && ps.set.normalisation; if (!c) return "";
  const to = c.to == null ? 1 : c.to, ax = members(ps)[0] ? members(ps)[0].axis.replace(/ \([^)]*\)$/, "") : "2θ";
  const r = c.from != null && c.until != null ? `${ax} ${fmt(c.from)} to ${fmt(c.until)}` : "";
  const ref = ps.set.norm_reference != null && ps.set.members[ps.set.norm_reference] ? ps.set.members[ps.set.norm_reference].label : "the reference";
  return { monitor: `each member's counts × ${fmt(to)} / its monitor`, time: `each member's counts × ${fmt(to)} / its counting time`, maximum: `each member's largest y${r ? ` over ${r}` : ""} made ${fmt(to)}`, area: `the same area, ${fmt(to)}, for every member over ${r}`, peak: `each member's integral over ${r} above its ${c.background === "line" ? "line between the ends" : "automatic background"} made ${fmt(to)}`, reference: `each member's area over ${r} made that of ${ref}` }[c.by] || c.by;
}
/** What the view does to the set, said in every view (PL6, PL10): its normalisation, its x, the members it leaves out. */
function viewText(ps) {
  const s = ps.set; if (!s) return "";
  const out = s.members.filter((m) => m.shown.left_out);
  const parts = [];
  if (s.view_error) parts.push(`The view as chosen could not be made (${s.view_error}): the members drawn as read.`);
  if (s.normalisation) parts.push(`Normalised by the core, a view (the files untouched): ${normText(ps)}.`);
  if (s.x_quantity) parts.push(`x in ${s.x_label}, converted by the core, a view${s.x_fallback ? `; λ ${fmt(s.x_fallback[0])} Å from ${s.x_fallback[1]} for the members whose files state none` : ""}.`);
  if (s.without_lambda.length) parts.push(`Without a wavelength, so not shown in ${s.x_label}: ${s.without_lambda.map((k) => s.members[k].label).join(", ")} (Convert X… gives one).`);
  const other = out.filter((m) => !s.without_lambda.includes(m.index));
  if (other.length) parts.push(`Left out of the view: ${other.map((m) => `${m.label} (${m.shown.left_out})`).join("; ")}.`);
  return parts.join(" ");
}

function setTitle(set, fallback) { return set.title || fallback; }
/** The project's folder, where Data keeps each file's operations (PL20) and where Plots/ is; null outside a project. */
const rootOf = () => (state.project && state.project.root) || null;
/** A tab for a set the core opened: the reference the warmest placed member until one is set; the settings a set file
    holds (PL12) applied when it came from one. */
function tabFor(set, title, path = null, settings = null, drv = null) {
  const id = tabOf(set.id);
  const tab = { id, module: "plots", path, title, setId: set.id, drv };
  const ps = stateOf(tab); ps.set = set;
  ps.timing = { t0: performance.now(), opened: performance.now(), firstFrame: null, whole: null, patterns: null };   // step 22: the clock of a large set's drawing (ms since the window's start)
  ps.reference = placed(ps).length ? placed(ps)[placed(ps).length - 1].key : null;   // the last placed member (the warmest, by the core's order)
  if (settings) applySettings(ps, settings);
  // the view's normalisation and x quantity the set file keeps (PL27), as the Rust side opened it
  if (set.normalisation) ps.nc = { ...set.normalisation, member: set.norm_reference };
  if (set.x_quantity) ps.xc = { quantity: set.x_quantity, control: null, typed: null };
  return tab;
}
/** The set file's settings (PL12) in the tab's state: the order, the stack's rule, the y scale, the reference, the
    windows, each member's opacity. */
function applySettings(ps, st) {
  ps.order = st.order || "variable";
  // PL27: the stack's step and the projection's steps apart, each with its scale. A file written before PL27 holds one
  // displacement, whose dx is the projection's (as Publish reads it, publish/DESIGN.md §15); a stack displacement that
  // only copies the projection's steps (written for Publish 0.8.0) is no stack step of its own
  const legacy = !st.prule || st.prule === "none";
  if (st.rule && st.rule !== "none") {
    const copy = !legacy && st.prule === st.rule && st.pdx === st.dx && st.pdy === st.dy;
    if (!copy) ps.rule = { kind: st.rule, dx: legacy ? 0 : st.dx || 0, dy: st.dy, proposed: null, scale: st.stack_scale || st.scale };
    if (legacy && st.dx) ps.prule = { kind: st.rule, dx: st.dx, dy: st.dy, scale: st.scale };
  }
  if (!legacy) ps.prule = { kind: st.prule, dx: st.pdx, dy: st.pdy, scale: st.pscale || st.scale };
  if (st.map_scale) ps.mapScale = st.map_scale;
  if (st.scale) ps.scale = st.scale;
  if (st.reference != null && ps.set.members[st.reference]) ps.reference = ps.set.members[st.reference].key;
  ps.wins = (st.windows || []).map((w) => ({ name: w.name, from: w.from, to: w.to }));
  ps.opacity = {}; (st.opacity || []).forEach((o, i) => { if (o != null && ps.set.members[i]) ps.opacity[ps.set.members[i].key] = o; });
  // PL37 (step 17): the map's and the difference map's clip, on the scale the file's map is drawn on, and the regions
  ps.clip = { map: st.map_range ? { from: st.map_range.from, to: st.map_range.to, scale: ps.mapScale } : null, difference: st.difference_range ? { from: st.difference_range.from, to: st.difference_range.to, scale: diffScaleOf(ps.mapScale) } : null };
  ps.regions = (st.regions || []).map((r) => ({ name: r.name, x: [r.x[0], r.x[1]], x_quantity: r.x_quantity || null, variable: [r.variable[0], r.variable[1]], map_range: r.map_range ? { from: r.map_range.from, to: r.map_range.to } : null, difference_range: r.difference_range ? { from: r.difference_range.from, to: r.difference_range.to } : null, window: r.window || null }));
  ps.regionsRefused = st.regions_refused || null; ps.selregion = null; ps.mapZoom = null; ps.rangeFromView = false;
  ps.yInterval = st.y_interval === "own" ? "own" : "shared";   // PL41 met by core 0.80.0: the set file keeps it
}
async function adopt(set, title, path = null, settings = null, drv = null) {
  const tab = tabFor(set, title, path, settings, drv);
  const ps = stateOf(tab);
  state.tabs.set(tab.id, tab);
  openTab(tab);
  if (window.__renderPanels) window.__renderPanels();
  const ops = set.members.filter((m) => m.operations.length).length;
  status(`Plots: ${members(ps).length} member${members(ps).length === 1 ? "" : "s"}${set.apart.length ? `, ${set.apart.length} without a ${set.variable_name} (drawn apart)` : ""}${ops ? `, ${ops} drawn after Data's operations` : ""}. ${drv ? `From ${drv}.` : "Not saved yet: Save (⌘S) writes it as Plots/{name}.drv."}`);
  return tab;
}
let untitled = 0;
async function newSet() {
  if (!hasCore()) return;
  try { const set = await invoke("plots_new", { root: rootOf() }); untitled += 1; const tab = await adopt(set, `Plot set ${untitled}`); await addFiles(tab); } catch (e) { status(String(e), true); }
}
async function openBuffer() {
  if (!hasCore()) return;
  const path = await pickFiles({ multiple: false, title: "Open a WinPLOTR buffer (.buf)", filters: [{ name: "WinPLOTR buffer", extensions: ["buf"] }] });
  if (!path) return;
  try { const set = await invoke("plots_open_buffer", { path, root: rootOf() }); await adopt(set, setTitle(set, path.split("/").pop()), path); } catch (e) { status(String(e), true); }
}
/** Open Plot Set… (PL12): a `.drv` read by the core, its members after Data's operations, its settings applied. */
async function openSetFile(path = null) {
  if (!hasCore()) return;
  if (!path) path = await pickFiles({ multiple: false, title: "Open a plot set (.drv)", defaultPath: rootOf() ? `${rootOf()}/Plots` : undefined, filters: [{ name: "Plot set", extensions: ["drv"] }] });
  if (!path) return;
  try { const [set, st] = await invoke("plots_open_set", { path, root: rootOf() }); await adopt(set, path.split("/").pop().replace(/\.drv$/i, ""), path, st, path.split("/").pop()); } catch (e) { status(String(e), true); }
}
/** The settings a set file keeps (PL12), from the tab's state. */
function settingsOf(ps) {
  const ms = members(ps);
  // PL27 (core 0.40.0): the stack's step, set by hand, with the scale it was set on, and the projection's steps with
  // theirs, apart; the map's colour scale. The core's proposed step is not written, so that it is proposed again.
  // Publish 0.8.0 still reads the projection's Δx from the stack's displacement (PL28): with no stack step of its own,
  // the projection's steps are written there too, and the y scale written is the one the steps were set on.
  const stack = ps.rule.dy != null && !ps.rule.proposed ? ps.rule : null;
  const proj = ps.prule.dy != null && ps.prule.dx ? ps.prule : null;
  const d = stack || proj;
  return { order: ps.order || "variable", rule: d ? d.kind : "none", dx: d ? d.dx || 0 : 0, dy: d ? d.dy : 0, scale: (d && d.scale) || ps.scale, stack_scale: d ? d.scale || ps.scale : null, prule: proj ? proj.kind : "none", pdx: proj ? proj.dx : 0, pdy: proj ? proj.dy : 0, pscale: proj ? proj.scale || ps.scale : null, map_scale: ps.mapScale, reference: ps.reference && ms.some((m) => m.key === ps.reference) ? ms.findIndex((m) => m.key === ps.reference) : null, windows: ps.wins, opacity: ms.map((m) => (ps.opacity[m.key] == null ? null : ps.opacity[m.key])),
    // PL37 (step 17): each clip written as { from, to } only while its scale is the set's map scale (the difference map's
    // on the scale it is drawn on for that map scale); the regions as the window keeps them
    y_interval: ps.yInterval === "own" ? "own" : "shared",   // PL41 met by core 0.80.0
    map_range: (() => { const c = clipOf(ps, "map", ps.mapScale); return c ? { from: c.from, to: c.to } : null; })(),
    difference_range: (() => { const c = clipOf(ps, "difference", diffScaleOf(ps.mapScale)); return c ? { from: c.from, to: c.to } : null; })(),
    regions: ps.regions.map((r) => ({ name: r.name, x: r.x, x_quantity: r.x_quantity || null, variable: r.variable, map_range: r.map_range || null, difference_range: r.difference_range || null, window: r.window || null })) };
}
/** Save (PL12): the set as Plots/{name}.drv through the core; the name asked the first time (Save As… asks again). */
async function saveSet(tab = activePlots(), askName = false) {
  if (!tab || !hasCore() || !stateOf(tab).set) return;
  const root = rootOf();
  if (!root || !writable()) { status("A plot set is saved in the project's Plots/ folder: open a project folder of your own (not a reference folder) first.", true); return; }
  let name = tab.saveName;
  if (!name || askName) {
    const v = await askText("Save Plot Set", "The set is written as Plots/{name}.drv in the project: its members (each file relative to the project's folder), each typed temperature with its source, the order, the stack's rule, the y scale, the reference, the windows and each member's opacity.", (tab.saveName || tab.title || "plot set").replace(/\.(buf|drv)$/i, ""), { label: "Name", ok: "Save" });
    if (v == null) return;
    name = v.trim().replace(/\//g, "-"); if (!name) return;
  }
  try {
    // W7: the host says where the file went ({ path, note }: the app's path in the project, the web host's download)
    const saved = await invoke("plots_save", { id: tab.setId, root, name, settings: settingsOf(stateOf(tab)) });
    const rel = typeof saved === "string" ? saved : saved.path, hostNote = typeof saved === "string" ? "" : saved.note || "";
    tab.saveName = name; tab.drv = rel; tab.dirty = false;
    const others = await reloadOthers(rel, tab);
    // step 17 (PL37 met): the clip and the regions go with the file
    const sv = settingsOf(stateOf(tab)); const kept = [sv.map_range ? "the map's colour range" : "", sv.difference_range ? "the difference map's colour range" : "", sv.regions.length ? `${sv.regions.length} region${sv.regions.length === 1 ? "" : "s"}` : ""].filter(Boolean);
    // P65 (step 19, core 0.78.0): a region's own range on each map goes with the file in its own field
    status(`${hostNote ? `${hostNote} (${rel})` : `Saved: ${rel}`}${kept.length ? `, with ${kept.join(", ")}` : ""}.${others ? ` Its other tab read again from the file.` : ""}`); bus.emit("file-saved", `${root}/${rel}`);
    if (window.__renderPanels) window.__renderPanels();
  } catch (e) { status(String(e), true); }
}
/** Send to Publish… (R30, PL15, PL29): the set saved as Plots/{name}.drv, with the steps and the scale Publish's views
    need, and handed to the Publish module (window.__sendTo, L220), which makes its publication set of the view shown
    and opens it; Build there typesets the figure. */
async function sendToPublish(tab = activePlots()) {
  if (!tab || !hasCore()) return null;
  const ps = stateOf(tab);
  await saveSet(tab);
  if (!tab.drv) return null;
  const st = settingsOf(ps);
  // the view Publish draws: the one shown, the stack from the overlay, and the map (or the difference map) from the
  // surface, whose numbers it draws on paper (DESIGN.md §8.1: the map and the projection carry them in print)
  const shown = ps.view === "surface" ? (ps.surfOf === "difference" ? "difference" : "map") : ps.view;
  const view = { stack: "F3, stack", projection: "F3, projection", map: "F3, map", difference: "F3, difference map", series: "the series view" }[shown] || "F3";
  const pview = ["stack", "projection", "map", "difference", "series"].includes(shown) ? shown : "stack";
  const lacks = [
    !st.dx ? "its projection needs the projection's Δx (Set Displacement… in the projection, then Save)" : "",
    st.scale === "log" ? "its difference map is refused on a set saved on the log scale (PL27, PL28)" : "",
    ps.set.normalisation || ps.set.x_quantity ? "a normalisation or an x quantity shown is not kept by the set file: Publish draws the counts as read, in 2θ (PL27)" : "",
  ].filter(Boolean);
  // the hand-over Publish defines (PL29, L220): it makes the publication set of this set and view, or opens the one
  // made before, and shows it; Build there typesets the figure with LuaLaTeX in Latin Modern
  let publishTab = null;
  if (window.__sendTo) publishTab = await window.__sendTo("publish", { drv: tab.drv, view: pview });
  status(`${tab.drv} saved and sent to Publish (${view})${publishTab ? `: ${publishTab.path} opened; Build (⌘B) there makes the figure in LaTeX` : window.__sendTo ? ": Publish made nothing (its status says why)" : ": Publish is not in this app"}.${lacks.length ? ` Publish: ${lacks.join("; ")}.` : ""}`, !!window.__sendTo && !publishTab);
  return { drv: tab.drv, settings: st, lacks, publishTab, view: pview };
}
/** Set Scale Bar… (step 7): the figure's scale bar in the y units shown, or the core's round step when left empty. */
async function setScaleBar(tab = activePlots()) {
  if (!tab) return;
  const ps = stateOf(tab);
  const v = await askForm("Set Scale Bar", `The scale bar a figure of the stack or the projection draws in place of the y axis, whose numbers hold for one member only; in ${yWord(ps)}, on the linear scale. Empty: the core's round step (plots::proposed_step, 1, 2 or 5 × 10ᵏ).`, [{ id: "v", label: `Length (${yWord(ps)})`, value: ps.sbar == null ? "" : ps.sbar, kind: "number" }], "Set");
  if (!v) return;
  if (v.v != null && !(isFinite(v.v) && v.v > 0)) { status("The scale bar's length is a positive number, or empty for the core's step.", true); return; }
  ps.sbar = v.v; ps.look = true; renderActive();
}
/** Figure Size… (steps 7 and 9): the width a figure is exported at, a journal's or a thesis's figure width as Publish's
    presets give it (drs-publish's presets, through the Rust side), each class (single, onehalf, double, full) with its
    preset's name. */
let widthsCache = null;
async function setFigureSize(tab = activePlots()) {
  if (!tab || !hasCore()) return;
  const ps = stateOf(tab);
  try { widthsCache = widthsCache || await invoke("plots_figure_widths"); } catch (e) { status(String(e), true); return; }
  const opts = widthsCache.map((w) => ({ value: `${w.preset}:${w.class}`, label: `${w.name.replace(/ \(.*$/, "")}: ${w.class} column, ${fmt(w.mm)} mm (text ${fmt(w.label_pt)}/${fmt(w.tick_pt)} pt there)` }));
  const cur = widthsCache.find((w) => w.mm === ps.figmm);
  const pick = await pickOne("Figure Size", "The width a figure is exported at in its figure look (Export ▸ Plot): a journal's or a thesis's figure width, from Publish's presets; its height in the window's proportion; its text at the preset's sizes, the margins measured from them.", opts, cur ? `${cur.preset}:${cur.class}` : opts[0] && opts[0].value, "Set");
  if (!pick) return;
  const w = widthsCache.find((x) => `${x.preset}:${x.class}` === pick);
  ps.figmm = w.mm; ps.figname = `${w.name.replace(/ \(.*$/, "")}, ${w.class}`; ps.figfont = { label_pt: w.label_pt, tick_pt: w.tick_pt }; ps.look = true; renderActive();
}
/** A set file written while another tab shows it (a tab the launch reopened from the last quit, P46, or one opened
    before): that tab read again from the file by the core, so that it shows what was written and its own Save does not
    put its older settings back. Returns how many tabs were read again. */
async function reloadOthers(rel, except) {
  const root = rootOf(); if (!root || !hasCore()) return 0;
  let n = 0;
  for (const x of [...state.tabs.values()]) {
    if (x.module !== "plots" || x === except || x.setId == null || (x.drv !== rel && x.path !== rel)) continue;
    try {
      const [set, st] = await invoke("plots_open_set", { path: `${root}/${rel}`, root });
      const ps = stateOf(x); const old = x.setId;
      x.setId = set.id; ps.set = set; applySettings(ps, st);
      ps.nc = set.normalisation ? { ...set.normalisation, member: set.norm_reference } : null;
      ps.xc = set.x_quantity ? { quantity: set.x_quantity, control: null, typed: null } : null;
      setSet(x, set);
      invoke("plots_close", { id: old }).catch(() => {});
      n++;
    } catch (e) { status(`${rel}: not read again in its other tab: ${e}`, true); }
  }
  return n;
}
/** Reload Members (PL20): the members read again by the core, after Data's operations as they are now. */
async function reloadMembers(tab = activePlots()) {
  if (!tab || !hasCore()) return;
  try { const s = await invoke("plots_reload", { id: tab.setId }); setSet(tab, s); const n = s.members.filter((m) => m.operations.length).length; status(`The members read again: ${n} drawn after Data's operations.`); } catch (e) { status(String(e), true); }
}
/** Read Variable From… a temperatures file (PL21): FullProf's "name TSAMP value" lines, read by the core. */
async function readTemperatures(tab) {
  const path = await pickFiles({ multiple: false, title: "FullProf's temperatures file (a data file's name and TSAMP with its temperature per line)", defaultPath: rootOf() || undefined });
  if (!path) return;
  try {
    const r = await invoke("plots_temperatures", { id: tab.setId, path });
    setSet(tab, r.set);
    status(`${path.split("/").pop()}: ${r.taken} member${r.taken === 1 ? "" : "s"} took their temperature from it${r.conflicts.length ? `; named twice with two values, so given nothing: ${r.conflicts.join(", ")}` : ""}${r.unread ? `; ${r.unread} line(s) not "name TSAMP value", kept unread` : ""}. Typed values are kept.`, !r.taken);
  } catch (e) { status(String(e), true); }
}
async function addFiles(tab = activePlots()) {
  if (!hasCore() || !tab) return;
  const paths = await pickFiles({ multiple: true, title: "Add data files to the plot set" });
  if (!paths || !paths.length) return;
  try { const set = await invoke("plots_add", { id: tab.setId, paths: [].concat(paths) }); setSet(tab, set); status(`${[].concat(paths).length} file(s) added: ${set.members.length} members.`); } catch (e) { status(String(e), true); }
}
/** A FullProf `.seq` (PL11) in a Plots tab of its own: its parameters against its external parameter (PARVAL). */
async function openSeqPath(path, title = path.split("/").pop(), rel = null) {
  const v = await invoke("plots_seq_open", { path });
  const id = "plots:seq:" + v.id;
  const tab = { id, module: "plots", path: rel, title, seqId: v.id };
  const ps = stateOf(tab); ps.seq = v; ps.view = "series"; ps.param = v.names.includes("Cell_A_ph1_pat1") ? "Cell_A_ph1_pat1" : v.names[0] || "chi2";
  return tab;
}
async function openSeq() {
  if (!hasCore()) return;
  const path = await pickFiles({ multiple: false, title: "Open FullProf's sequential results (.seq)", filters: [{ name: "FullProf sequential results", extensions: ["seq"] }] });
  if (!path) return;
  try { const tab = await openSeqPath(path); state.tabs.set(tab.id, tab); openTab(tab); if (window.__renderPanels) window.__renderPanels(); status(`${tab.title}: ${stateOf(tab).seq.blocks.length} refinements, ${stateOf(tab).seq.names.length} parameters.`); } catch (e) { status(String(e), true); }
}
/** A series of refinements (step 9): results files read by the core (results::read), each with its value of the
    variable from its refinement's data file as the core reads it, in a Plots tab of their own. */
async function openResultsPaths(paths, title = null) {
  const v = await invoke("plots_results_open", { paths });
  const id = "plots:results:" + v.id;
  const tab = { id, module: "plots", path: null, title: title || `${v.points.length} refinements`, rsId: v.id };
  const ps = stateOf(tab); ps.rs = v; ps.view = "series"; ps.xby = "variable";
  ps.param = v.names.find((n) => /cell_a$/.test(n)) || "chi2";
  return tab;
}
async function openResults() {
  if (!hasCore()) return;
  const paths = await pickFiles({ multiple: true, title: "Open the results files of a series (refinements against the temperature)", defaultPath: rootOf() ? `${rootOf()}/Results` : undefined, filters: [{ name: "Results file", extensions: ["cif"] }] });
  if (!paths || !paths.length) return;
  try {
    const tab = await openResultsPaths([].concat(paths)); state.tabs.set(tab.id, tab); openTab(tab); if (window.__renderPanels) window.__renderPanels();
    const v = stateOf(tab).rs; const out = v.points.filter((p) => !p.variable);
    status(`${v.points.length} refinements, ${v.names.length} refined parameters${out.length ? `; no temperature for ${out.map((p) => p.label).join(", ")}: double-click one in the sidebar to type it, with its source` : ""}.`, out.length > 0);
  } catch (e) { status(String(e), true); }
}
/** A refinement's temperature typed with its source (L169), for one whose data file says none. */
async function setResultsVariable(tab, k) {
  const ps = stateOf(tab); const p = ps.rs.points[k];
  const v = await askForm(`Set Temperature of ${p.label}`, `${p.variable ? `Now ${fmt(p.variable.value)} K, from ${p.variable.source}.` : `The core reads no temperature in its data file${p.data_file ? ` (${p.data_file.split("/").pop()})` : ""}${p.no_data ? `: ${p.no_data}` : ""}.`} A typed value keeps your note of where it came from (L169).`, [
    { id: "value", label: "T (K)", value: p.variable ? p.variable.value : "", kind: "number" },
    { id: "note", label: "Its source", value: "", placeholder: "e.g. the data file's header; the experiment's log" },
  ], "Set");
  if (!v) return;
  if (v.value == null || !isFinite(v.value)) { status("T is a number.", true); return; }
  if (!String(v.note || "").trim()) { status("Say where the value comes from: it is kept with it (L169).", true); return; }
  try { ps.rs = await invoke("plots_results_set_variable", { id: tab.rsId, k, value: v.value, note: String(v.note).trim() }); renderActive(); } catch (e) { status(String(e), true); }
}
/** Send Series to Publish… (step 10): the refinements of a results series handed to Publish as the series of its F5
    (the modes' amplitudes against the temperature) or F4 (the cell against it): each results file relative to the
    project with its temperature and that value's source. The hand-over is Publish's (PL34): it takes a series when it
    says so (`window.__sendTo.takes` holds "series"); until then nothing is sent, and the status says why. */
const publishTakesSeries = () => !!(window.__sendTo && Array.isArray(window.__sendTo.takes) && window.__sendTo.takes.includes("series"));
async function sendSeriesToPublish(tab = activePlots(), figure = null) {
  if (!tab || !hasCore() || tab.rsId == null) return null;
  const ps = stateOf(tab); const v = ps.rs; const root = rootOf();
  const out = v.points.filter((p) => !p.variable);
  if (out.length) { status(`Type the temperature of ${out.map((p) => p.label).join(", ")} first (a double-click in the sidebar): Publish's series needs a value for each refinement.`, true); return null; }
  if (!root || !v.points.every((p) => p.path.startsWith(root + "/"))) { status("Publish's series names its results files in the project: open the series from the project's Results/.", true); return null; }
  if (!figure) figure = await pickOne("Send Series to Publish", "What Publish draws of the series, in LaTeX: F5, the modes' amplitudes of each irrep against the temperature, or F4, the cell reduced to the parent's against it.", [{ value: "mode-amplitudes", label: "F5: the modes' amplitudes against T" }, { value: "cell", label: "F4: the cell against T" }], "mode-amplitudes", "Send");
  if (!figure) return null;
  const series = v.points.map((p) => ({ results: p.path.slice(root.length + 1), temperature: p.variable.value, source: p.variable.source }));
  if (!publishTakesSeries()) { status(`Publish does not take a series of refinements yet (request PL34): its F5 and F4 need a set naming the ${series.length} results files with their temperatures (${series.map((x) => `${fmt(x.temperature)} K`).join(", ")}), which Publish makes when it does.`, true); return { sent: false, series, figure }; }
  const t = await window.__sendTo("publish", { series, figure, title: tab.title });
  status(t ? `The series sent to Publish: ${t.path} opened; Build (⌘B) there typesets ${figure === "cell" ? "F4" : "F5"} in LaTeX.` : "Publish made nothing of the series (its status says why).", !t);
  return { sent: !!t, series, figure, publishTab: t };
}
/** The results series' quantities: χ², R_wp, R_p, then every refined parameter (results::read). */
function resultsQuantities(v) {
  return [["chi2", "χ²"], ["rwp", "R_wp"], ["rp", "R_p"], ...v.names.map((n) => [n, n])];
}
/** A refinement's value of a quantity, with its su where the results file gives one; null where it does not refine it. */
function resultsValue(p, q) {
  if (q === "chi2" || q === "rwp" || q === "rp") return [p[q], null];
  const r = p.refined.find((x) => x[0] === q); return r ? [r[1], r[2]] : null;
}

/** The `.seq` tab's quantities: every refined parameter (`Seq::names`), χ², and each phase's R_Bragg, volume and
    weight fraction (the block's own lines). */
function seqQuantities(v) {
  const phases = [...new Set(v.blocks.flatMap((b) => b.phases.map((p) => p.phase)))].sort((a, b) => a - b);
  return [["chi2", "χ² (GLOBAL_CHI_2)"], ...phases.flatMap((k) => [[`rb:${k}`, `R_Bragg, phase ${k} (%)`], [`vol:${k}`, `cell volume, phase ${k} (Å³)`], [`wf:${k}`, `weight fraction, phase ${k} (%)`]]), ...v.names.map((n) => [n, n])];
}
/** The set again from the core (after a change): the answers that depended on it are dropped. */
function setSet(tab, set) {
  const ps = stateOf(tab); ps.set = set; ps.disp = null; ps.map = null; ps.proj = null; ps.winres = null; ps.diff = null;
  ps.pin = null; ps.cell = null;   // step 16: the cursors' cell is of the map that is gone
  ps.sel = ps.sel.filter((k) => memberByKey(ps, k));
  if (ps.reference && !memberByKey(ps, ps.reference)) ps.reference = null;
  renderTab(tab); if (window.__renderPanels) window.__renderPanels(); renderBottom();
}

// ------------------------------------------------------------------ the view
const views = new Map();   // tab id → { redraw, plot, figure }
function renderTab(tab) { const v = views.get(tab.id); if (v) v.redraw(); }
function renderActive() { const t = activePlots(); if (t) renderTab(t); if (window.__renderPanels) window.__renderPanels(); renderBottom(); }

function viewOf(tab, box) {
  const ps = stateOf(tab);
  const isSeq = tab.seqId != null || tab.rsId != null;
  box.innerHTML = `<div class="ctabs">${VIEWS.filter(([k]) => !isSeq || k === "series").map(([k, t]) => `<span class="ctab" data-tab="${k}">${t}</span>`).join("")}
      <span class="axis pl-series" hidden><label>quantity <select class="pl-q">${QUANTITIES.map(([k, t]) => `<option value="${k}">${t}</option>`).join("")}<option value="background">background at x</option></select></label><label>background <select class="pl-wbg"><option value="auto">automatic</option><option value="line">line between the ends</option></select></label><button class="small pl-addwin" title="Add Window…: a named 2θ range">+ window</button></span>
      <span class="axis pl-seqctl" hidden><label>y <select class="pl-param"></select></label><label>x <select class="pl-xby"><option value="parameter">PARVAL's external parameter</option><option value="index">file order (PARVAL's index)</option></select></label></span>
      <span class="axis pl-scale"><label>y <select class="pl-y"><option value="linear">linear</option><option value="sqrt">√</option><option value="log">log</option></select></label></span>
      <span class="axis pl-linkctl" hidden><button class="small pl-link" title="Linked Zoom: the whole range on the left, a window of it on the right, zoom and cursor linked (F6)">Linked zoom</button></span>
      <span class="axis pl-rule" hidden><button class="small pl-disp" title="Set Displacement…: the step between members, or per unit of the variable">Δ</button><span class="pl-rule-text sub"></span></span>
      <span class="axis pl-winctl" hidden><button class="small pl-draw" title="Draw Window: drag across the plot to draw a window; drag its edges to move them; click a band to choose its window">Draw window</button></span>
      <span class="axis pl-lookctl" hidden><button class="small pl-look" title="Figure Look: the view as a publication figure would show it, on white, open axes, text in ink, a scale bar for a stack or a projection; Export ▸ Plot writes it at the figure's width (Figure Size…)">Figure look</button></span>
      <span class="axis pl-surfctl" hidden><label>of <select class="pl-sof"><option value="map">the map</option><option value="difference">the difference map</option></select></label></span>
      <span class="axis pl-mapscale" hidden><label>colour <select class="pl-ms"><option value="linear">linear</option><option value="sqrt">√</option><option value="log">log</option></select></label></span>
      <span class="axis pl-clipctl" hidden><button class="small pl-clip" title="Set Colour Range…: the range shown on the colour scale, typed; the two marks on the colour bar drag it; values beyond it are drawn saturated at the end colours; a double-click on the bar returns to the core's range">Range…</button><span class="pl-clip-text sub"></span></span>
      <span class="axis pl-curctl" hidden><button class="small pl-cur" title="Row and Column Cursors: the row under the pointer, its member's pattern below the map; the column, its values against the temperature beside it; a click pins them, a click again releases">Cursors</button><button class="small pl-follow" title="Range from the View: the colour range follows the zoom, taken from the core's block of the rectangle shown (plots::map_region); off puts the clip back as it was">Range from view</button><button class="small pl-saveregion" title="Save Region…: the rectangle shown saved by name in the set file, with its colour range when Range from the View is on">Save region…</button></span>
      <span class="axis pl-yintctl" hidden><button class="small pl-yint" title="The y interval: shared by the whole set (its extent, the core's) or adapted to the member chosen (Overlay, Stack) and to the row under the pointer (the strip under the map)">y: shared</button></span>
      <span class="subtitle pl-title"></span></div>
    <div class="pl-areas"><div class="pl-mapcol"><div class="plotarea pl-area"><canvas class="plot"></canvas><div class="readout"></div><div class="legend pl-legend" hidden></div></div>
        <div class="plotarea pl-rowstrip" hidden><canvas class="plot pl-rcanvas"></canvas><div class="readout pl-rread"></div></div></div>
      <div class="plotarea pl-colstrip" hidden><canvas class="plot pl-ccanvas"></canvas><div class="readout pl-cread"></div></div>
      <div class="plotarea pl-surf" hidden><div class="readout pl-sread"></div></div>
      <div class="plotarea pl-zoom" hidden><canvas class="plot pl-zcanvas"></canvas><div class="readout pl-zread"></div></div></div>
    <div class="pl-note sub"></div>`;
  const q = (s) => box.querySelector(s);
  const note = q(".pl-note"), ySel = q(".pl-y"), msSel = q(".pl-ms");
  q(".pl-clip").onclick = () => setColourRange(tab);
  q(".pl-cur").onclick = () => toggleCursors(tab);
  q(".pl-follow").onclick = () => toggleFollow(tab);
  q(".pl-saveregion").onclick = () => saveRegion(tab);
  q(".pl-yint").onclick = () => toggleYInterval(tab);
  box.querySelectorAll(".ctab").forEach((t) => { t.onclick = () => { ps.view = t.dataset.tab; redraw(); if (window.__renderPanels) window.__renderPanels(); }; });
  ySel.onchange = () => { ps.scale = ySel.value; ps.disp = null; redraw(); };
  msSel.onchange = () => { ps.mapScale = msSel.value; ps.map = null; redraw(); };
  q(".pl-sof").onchange = () => { ps.surfOf = q(".pl-sof").value; redraw(); };
  q(".pl-disp").onclick = () => setDisplacement(tab);
  q(".pl-draw").onclick = () => drawWindow(tab);
  q(".plot").addEventListener("click", (e) => { if (e.altKey && !isSeq && ps.hoverKey) soloToggle(tab, ps.hoverKey); });
  q(".pl-look").onclick = () => { ps.look = !ps.look; redraw(); if (window.__renderPanels) window.__renderPanels(); };
  q(".pl-q").onchange = () => { ps.quantity = q(".pl-q").value; redraw(); };
  q(".pl-wbg").onchange = () => { ps.wbg = q(".pl-wbg").value; ps.winres = null; redraw(); };
  q(".pl-addwin").onclick = () => addWindow(tab);
  q(".pl-link").onclick = () => { ps.linked = !ps.linked; redraw(); if (window.__renderPanels) window.__renderPanels(); };
  q(".pl-param").onchange = () => { ps.param = q(".pl-param").value; ps.seqcol = null; redraw(); };
  q(".pl-xby").onchange = () => { ps.xby = q(".pl-xby").value; redraw(); };
  // the plot's options stay this view's: the right margin follows the view (the stack's labels, the map's colour bar)
  const popts = { xLabel: "2θ (°)", yLabel: "counts", nearestInY: true, right: 16 };
  const legendEl = q(".pl-legend");
  const zoomEl = q(".pl-zoom");
  // F6: the overlay's window, a second plot linked to the first: a drag on the left chooses the window, a zoom or pan
  // on the right moves it, the pointer on either draws a cursor on both
  let zplot = null;
  // the Surface view (DESIGN.md §8.1 item 6): the core's map as a surface on the shared 3D canvas (surface3d.js), made
  // when the view is shown and let go when another one is, so that no GPU context is held by a view not shown (L150)
  const surfEl = q(".pl-surf");
  let surf = null, surfPending = null, surfFor = null;
  function ensureSurface() {
    if (surf) return Promise.resolve(surf);
    if (surfPending) return surfPending;
    return (surfPending = import("../../surface3d.js").then((m) => {
      surfPending = null;
      surf = m.createSurface(surfEl, { label: `${tab.title}: the map as a surface, the temperature in depth; the members table holds its numbers`, readout: q(".pl-sread") });
      surfFor = null; redraw(); return surf;
    }));
  }
  function dropSurface() { if (surf) { surf.dispose(); surf = null; surfFor = null; } }
  const linkPointer = (x) => { plot.setCursor(x); if (zplot) zplot.setCursor(x); };
  popts.onPointer = linkPointer;
  // the windows on the plot (DESIGN.md §4): drawn by a drag while Draw Window is on, their edges dragged, a band clicked
  // to choose its window; the series view follows at once (the core's windows again for the new ranges)
  popts.onBandCreate = (a, b) => {
    if (ps.drawWin) { ps.drawWin = false; addWindowRange(tab, a, b); return; }
    if (ps.view === "overlay" && ps.linked) { ps.zoomwin = [a, b]; redraw(); }
  };
  popts.onBandChange = (id, i, a, b) => {
    if (id === "zoomwin") { ps.zoomwin = [a, b]; redraw(); return; }
    if (id === "windows" && ps.wins[i]) moveWindow(tab, i, a, b);
  };
  // dimming by the pointer (DESIGN.md §4): the member under it highlighted, and ⌥-click solos it (again: all restored)
  const keyOfSeries = (id) => (id == null ? null : String(id).replace(/:\d+$/, ""));
  popts.onHover = (n) => {
    if (isSeq || !ps.set) return;
    const key = n && n.s ? keyOfSeries(n.s.id) : null;
    const k2 = key && memberByKey(ps, key) ? key : null;
    if (k2 !== ps.hoverKey && (ps.view === "overlay" || ps.view === "stack" || ps.view === "projection")) { ps.hoverKey = k2; redraw(); }
  };
  popts.onBandSelect = (id, i) => { if (id === "windows" && ps.wins[i]) { ps.selwin = i; if (window.__renderPanels) window.__renderPanels(); } };
  const plot = createPlot(q(".plot"), q(".readout"), legendEl, popts);
  /** The map shown now (the Map's or the Difference Map's, as the core gave it), or null. */
  function mapShown() { const diff = ps.view === "difference"; const got = ps.map && ps.map.key.startsWith(diff ? "d:" : "m:") ? ps.map : null; return got ? { got, diff } : null; }
  // step 18 (I08 step 5, L409): the Map view is the viewer's (app/ui/viewer/map.js), hosted here on the shared plot
  // with the tab's state as its state; steps 16 and 17 (the clip and its marks, the cursors' strips, the zoom kept,
  // range from the view, the regions) live there, drawn once for the app and the core's page. The host answers what
  // needs the core and draws again on every change the viewer tells.
  let panelsTimer = null;
  const host = {
    mapRegion: (region, c) => (hasCore() ? invoke("plots_map_region", { id: tab.setId, scale: c.scale, reference: c.reference, region }) : Promise.reject(new Error("the app's core is not here"))),
    rangeWarning: (range, c) => (hasCore() ? invoke("plots_range_warning", { id: tab.setId, scale: c.scale, reference: c.reference, range }) : Promise.resolve(null)),
    status: (text, bad) => status(text, !!bad),
    onChange: () => { renderTab(tab); clearTimeout(panelsTimer); panelsTimer = setTimeout(() => { if (window.__renderPanels) window.__renderPanels(); }, 120); },
    ask: (question, value) => askText("Plots", question, value, { label: "Name", ok: "OK" }),
  };
  const viewer = createMapViewer({ plot, plotOpts: popts, rowEl: q(".pl-rowstrip"), rowCanvas: q(".pl-rcanvas"), rowReadout: q(".pl-rread"), colEl: q(".pl-colstrip"), colCanvas: q(".pl-ccanvas"), colReadout: q(".pl-cread"), note, clipText: q(".pl-clip-text") }, host, ps);
  /** The viewer's data from the tab (PL38's shape; the core's `page_data` writes the same for the page): the maps as
      the core gave them, the members with their patterns as the set holds them, the regions (their blocks asked of the
      core through the host, so none precomputed here). */
  function viewerData() {
    const set = ps.set;
    const mapOf = (key) => { const got = ps.map && ps.map.key.startsWith(key === "difference" ? "d:" : "m:") ? ps.map : null; return got ? { map: got.v.map, apart: got.v.apart, scale: got.scale, reference: key === "difference" ? got.ref : null } : null; };
    return {
      title: tab.title, path: tab.drv || tab.path || null, core: null, written: null,
      variable: { name: set.variable_name, unit: set.variable_unit }, x: { label: xLabelOf(ps), quantity: xWordOf(ps) }, y_word: yWord(ps),
      members: set.members.map((m) => ({ key: m.key, label: m.label, variable: m.variable ? { value: m.variable.value, source: m.variable.source } : null, pattern: m.shown.left_out ? null : patternOf(fetchPatterns(), m), kind: m.format, kind_name: m.format_name, points: m.points, counting_time: m.counting_time, time: m.modified, others: m.others || [], facts: memberFacts(ps, m) })),
      maps: { map: mapOf("map"), difference: mapOf("difference") },
      clip: { map: null, difference: null },
      regions: ps.regions.map((r) => ({ ...r, blocks: { map: null, difference: null }, refused: null })),
      regions_refused: ps.regionsRefused, same_grid: set.same_grid,
      y_extent: ps.extent && ps.map && ps.extent.key.startsWith(ps.map.scale + ":") ? ps.extent.v : null,
    };
  }
  const zopts = { xLabel: "2θ (°)", yLabel: "counts", nearestInY: true, onPointer: linkPointer,
    onView: (x0, x1) => { if (ps.view !== "overlay" || !ps.linked) return; const w = ps.zoomwin; if (w && w[0] === x0 && w[1] === x1) return; ps.zoomwin = [x0, x1]; drawLeftBand(); } };
  function ensureZoom() {
    if (zplot) return zplot;
    zplot = createPlot(q(".pl-zcanvas"), q(".pl-zread"), null, zopts);
    return zplot;
  }
  /** The set's y extent on a scale, the core's (`plots_extent`, PL41 asks `plots::extent`), asked once per scale and
      view rule: the shared y interval of the row strip under the map (step 20, W4). */
  function fetchExtent(scale) {
    const key = `${scale}:${ps.set && ps.set.normalisation ? JSON.stringify(ps.set.normalisation) : ""}:${ps.set ? ps.set.x_quantity || "" : ""}`;
    if (ps.extent && ps.extent.key === key) return ps.extent.v;
    if (asked === "extent:" + key || !hasCore()) return null;
    asked = "extent:" + key;
    invoke("plots_extent", { id: tab.setId, scale, members: null }).then((v) => { asked = null; ps.extent = { key, v: v ? { from: v[0], to: v[1] } : null }; redraw(); }).catch((e) => { asked = null; status(String(e), true); });
    return null;
  }
  let leftSeries = [];
  function drawLeftBand() {
    const band = ps.zoomwin ? [{ id: "zoomwin", kind: "bands", x: [], ranges: [ps.zoomwin.slice()], editable: true, color: "--accent" }] : [];
    plot.setSeries([...band, ...leftSeries], true);
  }
  let asked = null;   // the request in flight, so that one view asks once

  async function fetchDisp() {
    const r = ps.rule; const key = `${r.kind}:${r.dx}:${r.dy}:${ps.scale}`;
    if (ps.disp && ps.disp.key === key) return ps.disp.v;
    if (asked === key || !hasCore()) return null;
    asked = key;
    try { const v = await invoke("plots_displace", { id: tab.setId, rule: r.kind, dx: r.dx || 0, dy: r.dy || 0, scale: ps.scale }); ps.disp = { key, v }; asked = null; redraw(); return v; } catch (e) { asked = null; status(String(e), true); return null; }
  }
  /** The members' patterns as the view shows them (W9 a, step 22): asked of the core once per view rule
      (`plots_patterns`, every member; a set view carries no x and y since core 0.81.0 and the app's bridge of this
      step), kept by the tab as a map by member index; a set view that still carries them (an older core) is read as
      it is. Null while the core answers (the view draws again when they come). */
  function fetchPatterns() {
    const set = ps.set; if (!set) return null;
    if (set.members.length && set.members.every((m) => Array.isArray(m.x))) return new Map(set.members.map((m, i) => [i, { x: m.x, y: m.y, sigma: m.sigma || null }]));
    const key = `${set.id}:${set.normalisation ? JSON.stringify(set.normalisation) : ""}:${set.x_quantity || ""}:${set.members.map((m) => m.key).join("|")}`;
    if (ps.patterns && ps.patterns.key === key) return ps.patterns.v;
    if (asked === "patterns:" + key || !hasCore()) return null;
    asked = "patterns:" + key;
    const t0 = performance.now();
    invoke("plots_patterns", { id: tab.setId, members: null }).then((list) => {
      asked = null;
      ps.patterns = { key, v: new Map(list.map((p) => [p.member, { x: p.x, y: p.y, sigma: p.sigma }])) };
      if (ps.timing) ps.timing.patterns = performance.now() - t0;
      redraw();
    }).catch((e) => { asked = null; status(String(e), true); });
    return null;
  }
  async function fetchMap(diff) {
    const ref = diff ? members(ps).findIndex((m) => m.key === ps.reference) : -1;
    const scale = diff && ps.mapScale === "log" ? "sqrt" : ps.mapScale;   // a difference has no log: √ keeps its sign and shows the small ones
    const key = `${diff ? "d" : "m"}:${scale}:${ref}`;
    if (ps.map && ps.map.key === key) return ps.map;
    if (asked === key || !hasCore()) return null;
    asked = key;
    try { const v = await invoke("plots_map", { id: tab.setId, scale, reference: diff ? ref : null }); ps.map = { key, v, scale, ref }; asked = null; redraw(); return ps.map; } catch (e) { asked = null; status(String(e), true); return null; }
  }

  async function fetchProj() {
    const r = ps.prule; const key = `${r.kind}:${r.dx}:${r.dy}:${ps.scale}`;
    if (ps.proj && ps.proj.key === key) return ps.proj.v;
    if (asked === key || !hasCore()) return null;
    asked = key;
    try { const v = await invoke("plots_project", { id: tab.setId, rule: r.kind, dx: r.dx || 0, dy: r.dy || 0, scale: ps.scale, floor: null, order: null }); ps.proj = { key, v }; asked = null; redraw(); return v; } catch (e) { asked = null; status(String(e), true); return null; }
  }
  async function fetchWins() {
    const key = JSON.stringify([ps.wins, ps.wbg, ps.at]);
    if (ps.winres && ps.winres.key === key) return ps.winres.v;
    if (asked === key || !hasCore()) return null;
    asked = key;
    try { const v = await invoke("plots_windows", { id: tab.setId, windows: ps.wins, background: ps.wbg, logLambda: null, at: ps.at }); ps.winres = { key, v }; asked = null; redraw(); return v; } catch (e) { asked = null; status(String(e), true); return null; }
  }
  async function fetchDiff() {
    const pair = diffPair(ps); if (!pair) return null;
    const key = pair.join(":");
    if (ps.diff && ps.diff.key === key) return ps.diff.v;
    if (asked === key || !hasCore()) return null;
    asked = key;
    try { const v = await invoke("plots_difference", { id: tab.setId, a: pair[0], b: pair[1] }); ps.diff = { key, v }; asked = null; redraw(); return v; } catch (e) { asked = null; status(String(e), true); return null; }
  }
  async function fetchSeqCol() {
    const name = ps.param, key = "seq:" + name;
    if (ps.seqcol && ps.seqcol.key === key) return ps.seqcol.v;
    if (asked === key || !hasCore()) return null;
    asked = key;
    try { const v = await invoke("plots_seq_column", { id: tab.seqId, name }); ps.seqcol = { key, v }; asked = null; redraw(); return v; } catch (e) { asked = null; status(String(e), true); return null; }
  }

  /** A figure's scale bar (steps 6 and 7), for the stack and the projection, whose y numbers hold for one member only:
      the length set by Set Scale Bar…, else the core's round step on the scale shown (plots::proposed_step: 1, 2 or
      5 × 10ᵏ); none on √ or log, where a length is not the same everywhere. */
  function scaleBarOf(set, yw) {
    if (ps.scale !== "linear") return null;
    if (ps.sbar != null) return { value: ps.sbar, label: `${fmt(ps.sbar)} ${yw}`, own: true };
    const key = `${set.id}:${ps.scale}:${set.normalisation ? JSON.stringify(set.normalisation) : ""}`;
    if ((!ps.pbar || ps.pbar.key !== key) && hasCore() && asked !== "pbar") {
      asked = "pbar"; invoke("plots_proposed_step", { id: tab.setId, scale: ps.scale }).then((p) => { asked = null; ps.pbar = { key, v: p }; redraw(); }).catch((e) => { asked = null; status(String(e), true); });
    }
    const p = ps.pbar && ps.pbar.key === key ? ps.pbar.v : null;
    return p ? { value: p.dy, label: `${fmt(p.dy)} ${yw}`, own: false } : null;
  }
  /** What the figure's look does, said in the note of the view shown. */
  const lookText = (bar) => `Figure look: on white, open axes, text in ink${bar ? `, a scale bar of ${bar.label} (${bar.own ? "set by Set Scale Bar…" : "the core's round step, plots::proposed_step; Set Scale Bar… gives another"}) in place of the y axis` : ps.view === "stack" || ps.view === "projection" ? (ps.scale === "linear" ? "" : ", the y axis kept (a scale bar is drawn on the linear scale only)") : ""}; Export ▸ Plot writes it ${fmt(ps.figmm)} mm wide, ${ps.figname}, its text at ${fmt(ps.figfont.label_pt)} and ${fmt(ps.figfont.tick_pt)} pt (Figure Size…). A quick look in the app's typeface: a paper's figure is Publish's, in LaTeX. `;

  /** The results series (step 9): a quantity of every refinement against its temperature or the file order. */
  function redrawResults() {
    const v = ps.rs;
    box.querySelectorAll(".ctab").forEach((t) => t.classList.toggle("active", t.dataset.tab === "series"));
    q(".pl-scale").hidden = true; q(".pl-rule").hidden = true; q(".pl-mapscale").hidden = true; q(".pl-series").hidden = true; q(".pl-seqctl").hidden = false; q(".pl-lookctl").hidden = false; legendEl.hidden = true;
    q(".pl-look").classList.toggle("on", !!ps.look);
    const opts = resultsQuantities(v).map(([k, t]) => `<option value="${esc(k)}">${esc(t)}</option>`).join("");
    if (q(".pl-param").innerHTML !== opts) q(".pl-param").innerHTML = opts;
    const xo = `<option value="variable">T (K), from each refinement's data file</option><option value="index">file order</option>`;
    if (q(".pl-xby").innerHTML !== xo) q(".pl-xby").innerHTML = xo;
    q(".pl-param").value = ps.param; q(".pl-xby").value = ps.xby === "index" ? "index" : "variable";
    q(".pl-title").textContent = `${tab.title} · ${v.points.length} refinements`;
    popts.right = 16; popts.yFromZero = false;
    const label = (resultsQuantities(v).find(([k]) => k === ps.param) || [0, ps.param])[1];
    const xs = [], ys = [], sg = [], apart = [];
    v.points.forEach((p, k) => {
      const y = resultsValue(p, ps.param); if (!y) return;
      const x = ps.xby === "index" ? k + 1 : p.variable ? p.variable.value : null;
      if (x == null) { apart.push(p); return; }
      xs.push(x); ys.push(y[0]); sg.push(y[1]);
    });
    plot.setAxes({ xLabel: ps.xby === "index" ? "file order" : "T (K)", yLabel: label, yScale: "linear", look: ps.look ? "figure" : null, xSpan: null, scaleBar: null });
    plot.setSeries([{ id: "results", kind: "both", points: true, marker: "circle", label, x: xs, y: ys, sigma: sg.some((e) => e != null) ? sg.map((e) => (e == null ? 0 : e)) : null, color: "#0072B2", width: 1 }]);
    const notRefined = v.points.filter((p) => !resultsValue(p, ps.param)).map((p) => p.label);
    note.textContent = `${label} in each of the ${v.points.length} refinements (their results files, read by the core: results::read), with its su where the file gives one, against ${ps.xby === "index" ? "the order the files were opened in" : "each refinement's temperature, read by the core from its data file (its header, title or name) or typed with its source"}; ${xs.length} point(s).${apart.length ? ` No temperature, so not drawn: ${apart.map((p) => `${p.label} (${p.no_data || "its data file says none the core reads"})`).join("; ")}; a double-click in the sidebar types one.` : ""}${notRefined.length ? ` Not refined in: ${notRefined.join(", ")}.` : ""}`;
  }
  /** The `.seq` tab (PL11): a quantity of every refinement against PARVAL's external parameter or the file order. */
  function redrawSeq() {
    const v = ps.seq;
    box.querySelectorAll(".ctab").forEach((t) => t.classList.toggle("active", t.dataset.tab === "series"));
    q(".pl-scale").hidden = true; q(".pl-rule").hidden = true; q(".pl-mapscale").hidden = true; q(".pl-series").hidden = true; q(".pl-seqctl").hidden = false; q(".pl-lookctl").hidden = true; legendEl.hidden = true;
    const opts = seqQuantities(v).map(([k, t]) => `<option value="${esc(k)}">${esc(t)}</option>`).join("");
    if (q(".pl-param").innerHTML !== opts) q(".pl-param").innerHTML = opts;
    q(".pl-param").value = ps.param; q(".pl-xby").value = ps.xby;
    q(".pl-title").textContent = `${tab.title} · ${v.blocks.length} refinements`;
    popts.right = 16; popts.yFromZero = false;
    let pts = null;   // [value, σ] per block, as the core gives them
    if (ps.param === "chi2") pts = v.blocks.map((b) => (b.chi2 == null ? null : [b.chi2, null]));
    else if (/^(rb|vol|wf):/.test(ps.param)) {
      const [what, k] = ps.param.split(":"); const ph = (b) => b.phases.find((p) => p.phase === Number(k));
      pts = v.blocks.map((b) => { const p = ph(b); if (!p) return null; return what === "rb" ? (p.r_bragg == null ? null : [p.r_bragg, null]) : what === "vol" ? p.volume : p.weight_fraction; });
    } else { const c = fetchSeqCol(); void c; pts = ps.seqcol && ps.seqcol.key === "seq:" + ps.param ? ps.seqcol.v : null; }
    if (!pts) { note.textContent = hasCore() ? "The core reads the column…" : "The .seq needs the app (the core reads it)."; return; }
    const xs = [], ys = [], sg = [];
    v.blocks.forEach((b, i) => { const x = ps.xby === "index" ? b.index : b.parameter; const p = pts[i]; if (x == null || !p) return; xs.push(x); ys.push(p[0]); sg.push(p[1]); });
    const label = (seqQuantities(v).find(([k]) => k === ps.param) || [0, ps.param])[1];
    plot.setAxes({ xLabel: ps.xby === "index" ? "file order (PARVAL's index)" : "PARVAL's external parameter", yLabel: label, yScale: "linear" });
    plot.setSeries([{ id: "seq", kind: "both", points: true, label, x: xs, y: ys, sigma: sg.some((e) => e != null) ? sg.map((e) => (e == null ? 0 : e)) : null, color: "--phase-0", width: 1 }]);
    const warn = [v.repeated.length ? ` Refined more than once (a run appended to an earlier one): ${v.repeated.join(", ")}.` : "", v.inconsistent.length ? ` ${v.inconsistent.length} block(s) list fewer or more parameters than their NUMPAR.` : "", v.unread ? ` ${v.unread} line(s) the core did not read, kept.` : ""].join("");
    note.textContent = `${label} in each of the ${v.blocks.length} refinements of ${tab.title} (FullProf's sequential results, read by the core: seq::read), with its σ where the file gives one, against ${ps.xby === "index" ? "the file order (PARVAL's first number)" : "PARVAL's second number, the external parameter of each refinement (in FullProf's LaMnO₃ example, the data file's sample temperature)"}; ${xs.length} point(s).${warn}`;
  }

  /** The chosen window's band marked on the plot. */
  function markWin() { if (ps.selwin != null && ps.wins[ps.selwin]) plot.setSelection("windows", [ps.selwin]); else if (plot.state.sel.has("windows")) plot.clearSelection(); }
  function redraw() {
    if (isSeq) { if (ps.rs) redrawResults(); else if (ps.seq) redrawSeq(); return; }
    if (!ps.set) return;
    redrawSet();
    const vt = viewText(ps); if (vt && ps.set.members.length) note.textContent = `${vt} ${note.textContent}`;
  }
  function redrawSet() {
    const set = ps.set;
    // the figure's look (steps 6 and 7): every view as a publication figure would show it; the stack and the projection
    // add their x span and scale bar below
    const figLook = ps.look;
    const surface = ps.view === "surface";
    q(".pl-mapcol").style.display = surface ? "none" : ""; surfEl.hidden = !surface; q(".pl-surfctl").hidden = !surface; q(".pl-sof").value = ps.surfOf;
    if (!surface) dropSurface();
    q(".pl-lookctl").hidden = surface; q(".pl-look").classList.toggle("on", figLook);
    const mapView = ps.view === "map" || ps.view === "difference";
    if (!mapView) viewer.hide();   // step 18: the viewer's marks and strips go with its view; its zoom stays in the state
    plot.setAxes({ look: figLook ? "figure" : null, xSpan: null, scaleBar: null, keepView: mapView });
    box.querySelectorAll(".ctab").forEach((t) => t.classList.toggle("active", t.dataset.tab === ps.view));
    const isMap = ps.view === "map" || ps.view === "difference" || surface;
    const mapDiff = ps.view === "difference" || (surface && ps.surfOf === "difference");
    // step 16: the clip's control in the map views and the surface; the cursors' controls in the map views only
    q(".pl-clipctl").hidden = !isMap; q(".pl-curctl").hidden = !mapView; q(".pl-cur").classList.toggle("on", !!ps.cursors); q(".pl-follow").classList.toggle("on", !!ps.rangeFromView);
    // step 20 (W4): the y interval shared by the set or adapted to the member chosen: Overlay and Stack fit their y to
    // the chosen members' series (the shared plot's own fit, screen layout), the strip under the map takes the set's
    // extent from the core or fits its row
    const curves = ps.view === "overlay" || ps.view === "stack";
    q(".pl-yintctl").hidden = !(curves || mapView); q(".pl-yint").textContent = `y: ${ps.yInterval === "own" ? "own" : "shared"}`; q(".pl-yint").classList.toggle("on", ps.yInterval === "own");
    popts.fitTo = curves && ps.yInterval === "own" && ps.sel.length ? ps.sel.slice() : null; zopts.fitTo = popts.fitTo;
    q(".pl-scale").hidden = isMap; q(".pl-rule").hidden = ps.view !== "stack" && ps.view !== "projection"; q(".pl-mapscale").hidden = !isMap;
    q(".pl-series").hidden = ps.view !== "series"; q(".pl-seqctl").hidden = true; legendEl.hidden = ps.view !== "series";
    const linked = ps.view === "overlay" && ps.linked && !ps.look;   // a figure is one plot
    q(".pl-linkctl").hidden = ps.view !== "overlay"; q(".pl-link").classList.toggle("on", linked); zoomEl.hidden = !linked;
    plot.setMode(linked || ps.drawWin ? "region" : "zoom"); if (!linked) plot.setCursor(null);
    q(".pl-winctl").hidden = !(ps.view === "overlay" || ps.view === "stack") || !inFilesX(ps) || ps.look; q(".pl-draw").classList.toggle("on", !!ps.drawWin);
    q(".pl-q").value = ps.quantity; q(".pl-wbg").value = ps.wbg;
    [...q(".pl-q").options].forEach((o) => { o.disabled = o.value === "background" && !ps.at.length; });
    ySel.value = ps.scale; msSel.value = ps.mapScale;
    [...msSel.options].forEach((o) => { o.disabled = mapDiff && o.value === "log"; });
    if (mapDiff && msSel.value === "log") msSel.value = "sqrt";
    q(".pl-title").textContent = `${tab.title} · ${set.members.length} member${set.members.length === 1 ? "" : "s"}${set.normalisation ? " · normalised" : ""}${set.x_quantity ? ` · in ${xNameOf(ps)}` : ""}`;
    const xl = xLabelOf(ps), yw = yWord(ps);
    const mon = monitorNote(ps);
    const apart = set.apart.map((k) => set.members[k].label);
    const apartText = apart.length ? ` Drawn apart, no ${set.variable_name} in the file: ${apart.join(", ")} (Set ${set.variable_name === "T" ? "Temperature" : "Value"}… types one, with its source).` : "";
    if (!set.members.length) { plot.setSeries([]); note.textContent = "No member yet: Plots ▸ Add Files… or Open Buffer…. In memory: the set is not saved (the plot-set file waits for the core, PL12)."; return; }
    popts.yFromZero = ps.view !== "series";
    // the right margin holds what each view writes there: the shared plot measures it (PL30 d, app 0.38.0) from the
    // stack's baseline labels, the projection's axis labels (point labels, PL35) and the map's colour bar
    popts.right = ps.view === "overlay" || ps.view === "series" ? 16 : "auto";
    // editable where they can be drawn: the overlay and the stack, in the files' x
    const windowBands = ps.wins.length && inFilesX(ps) ? [{ id: "windows", kind: "bands", x: [], ranges: ps.wins.map((w) => [w.from, w.to]), color: "--accent", editable: (ps.view === "overlay" || ps.view === "stack") && !ps.look }] : [];
    // step 22 (W9 b): the overlay and the stack of many members drawn in slices across frames by the shared plot, the
    // status saying how far; the tab's clock takes the first frame and the whole drawing (a large set's timing)
    const many = set.members.length > 30 && (ps.view === "overlay" || ps.view === "stack");
    popts.slice = many ? { per: 10, onProgress: (done, total) => { if (ps.timing && ps.timing.firstFrame == null) ps.timing.firstFrame = performance.now() - ps.timing.t0; if (done < total) status(`${ps.view === "stack" ? "Stack" : "Overlay"}: drawing ${done} of ${total} members…`); else { if (ps.timing && ps.timing.whole == null) ps.timing.whole = performance.now() - ps.timing.t0; status(`${ps.view === "stack" ? "Stack" : "Overlay"}: ${total} members drawn${ps.timing && ps.timing.firstFrame != null ? ` (first frame ${Math.round(ps.timing.firstFrame)} ms, whole ${Math.round(ps.timing.whole)} ms after the set)` : ""}.`); } } } : null;
    zopts.slice = popts.slice;
    if (ps.view === "overlay") {
      plot.setAxes({ xLabel: xl, yLabel: yw, yScale: ps.scale });
      const pats = fetchPatterns();
      if (!pats) { note.textContent = hasCore() ? "The core gives the members' patterns…" : "The overlay needs a core (the members' patterns)."; return; }
      const curves = set.members.filter((m) => !m.shown.left_out && patternOf(pats, m)).map((m) => ({ id: m.key, kind: "line", label: m.label, x: patternOf(pats, m).x, y: patternOf(pats, m).y, color: colourOf(ps, m), opacity: opacityOf(ps, m), width: widthOf(ps, m), legend: false }));
      if (!linked) {
        const t0 = performance.now();
        plot.setSeries([...windowBands, ...curves], true); markWin();
        // the overlay's first frame (the axes and the first slice, drawn at once) from the set's arrival; the whole
        // drawing when the last slice is done (onProgress), at once when not sliced
        if (ps.timing && ps.timing.firstFrame == null) { ps.timing.firstFrame = performance.now() - ps.timing.t0; if (!popts.slice) ps.timing.whole = performance.now() - ps.timing.t0; }
        void t0;
        note.textContent = `Every member as measured, in the core's order by ${set.variable_name}; chosen members bright, the others dimmed (Solo Member, Dim Others, the opacity bars in the sidebar).${mon ? " " + mon + "." : ""}${apartText}`;
      } else {
        // F6: the whole range on the left, the window on the right, the core's difference of two members under both
        const pair = diffPair(ps); fetchDiff();
        const dv = ps.diff && pair && ps.diff.key === pair.join(":") ? ps.diff.v : null;
        const strip = dv ? [{ id: "difference", kind: "line", strip: true, label: `${set.members[pair[0]].label} − ${set.members[pair[1]].label}`, legend: false, x: dv.x, y: dv.y, color: "--ink", width: 0.8 }] : [];
        leftSeries = [...windowBands, ...curves, ...strip];
        drawLeftBand();
        const z = ensureZoom();
        z.setAxes({ xLabel: xl, yLabel: yw, yScale: ps.scale });
        z.setSeries([...windowBands, ...curves, ...strip], true);
        const zv = z.getView();
        if (ps.zoomwin && (zv.x0 !== ps.zoomwin[0] || zv.x1 !== ps.zoomwin[1])) z.setView({ x0: ps.zoomwin[0], x1: ps.zoomwin[1], y0: null, y1: null });
        else if (!ps.zoomwin) z.fitAll();
        note.textContent = `Left, the whole range; right, the window (drag across the left plot to choose it; zoom or pan the right one and the band on the left follows); the pointer on either draws the cursor on both. Under both, ${dv ? `the core's difference ${set.members[pair[0]].label} − ${set.members[pair[1]].label} (plots::derive, σ propagated)` : pair ? "the core's difference…" : "the difference of two members, when two are chosen or a reference is set"}: the first member chosen (or the first of the set) less the reference (Set Reference), or less the second chosen. Chosen members bright, the others dimmed.${mon ? " " + mon + "." : ""}${apartText}`;
      }
    } else if (ps.view === "stack") {
      if (ps.rule.dy == null || (ps.rule.proposed && ps.rule.proposed.scale !== ps.scale)) {
        // PL19: the core proposes the step on the scale shown, and says its rule; Set Displacement… overrides it
        const key = "step:" + ps.scale;
        if (hasCore() && asked !== key) { asked = key; invoke("plots_proposed_step", { id: tab.setId, scale: ps.scale }).then((p) => { asked = null; if (p) { ps.rule = { kind: "step", dx: 0, dy: p.dy, proposed: { scale: ps.scale, rule: p.rule } }; ps.disp = null; } else ps.rule = { ...ps.rule, proposed: null }; redraw(); }).catch((e) => { asked = null; status(String(e), true); }); }
      }
      q(".pl-rule-text").textContent = ps.rule.dy == null ? "no step set" : ps.rule.proposed ? `Δy ${fmt(ps.rule.dy)} per member, the core's proposal` : ps.rule.kind === "variable" ? `Δy ${fmt(ps.rule.dy)} per ${set.variable_unit}${ps.rule.dx ? `, Δx ${fmt(ps.rule.dx)}` : ""}` : `Δy ${fmt(ps.rule.dy)} per member${ps.rule.dx ? `, Δx ${fmt(ps.rule.dx)}` : ""}`;
      if (ps.rule.dy == null) { plot.setSeries([]); note.textContent = "The stack raises each member by a step: the core proposes one from the data (PL19), or Set Displacement… (Δ) gives it, in the y units shown (counts on linear, √counts on √, decades on log), per member or per unit of the variable."; return; }
      const d = fetchDisp(); void d;
      const v = ps.disp && ps.disp.v;
      if (!v) { note.textContent = hasCore() ? "The core displaces the members…" : "The stack needs the app (the core displaces the members)."; return; }
      const sbar = figLook ? scaleBarOf(set, yw) : null;
      plot.setAxes(figLook ? { xLabel: xl, yLabel: sbar ? "" : `${yw}, bottom member`, yScale: ps.scale, look: "figure", xSpan: null, scaleBar: sbar } : { xLabel: xl, yLabel: `${yw} + each member's step`, yScale: ps.scale });
      const series = [...windowBands];
      const pats = fetchPatterns();   // the readout's counts (`values`): the member's pattern as the set holds it
      v.forEach((dv) => {
        const m = set.members[dv.member]; if (!dv.placed) return;
        const p = patternOf(pats, m);
        series.push({ id: m.key, kind: "line", label: m.label, legend: false, x: dv.x, y: dv.y, values: p ? p.y : dv.y, scaled: true, color: colourOf(ps, m), opacity: opacityOf(ps, m), width: widthOf(ps, m), baselineLabel: { text: m.variable ? `${fmt(m.variable.value)} ${set.variable_unit}` : ps.look ? m.label : `${m.label} (no ${set.variable_name})`, at: dv.baseline, priority: (m.variable ? 0 : 2) + (ps.sel.includes(m.key) ? 1 : 0) } });
      });
      plot.setSeries(series, true); markWin();
      const unplaced = v.filter((dv) => !dv.placed).map((dv) => set.members[dv.member].label);
      note.textContent = `${figLook ? lookText(sbar) : ""}${ps.rule.proposed ? `The step is the core's proposal (plots::proposed_step): ${ps.rule.proposed.rule}; Set Displacement… (Δ) gives another. ` : ""}Each member raised by the core's displacement (plots::displace), its label at its baseline; the axis's numbers hold for the member at the bottom only.${unplaced.length ? ` Not placed by a step per ${set.variable_unit}: ${unplaced.join(", ")} (no ${set.variable_name}).` : ""}${mon ? " " + mon + "." : ""}${apartText}`;
    } else if (ps.view === "projection") {
      const r = ps.prule;
      q(".pl-rule-text").textContent = r.dy == null ? "no steps set" : `Δx ${fmt(r.dx)}°, Δy ${fmt(r.dy)} per ${r.kind === "variable" ? set.variable_unit : "member"}`;
      if (r.dy == null) { plot.setSeries([]); note.textContent = "The projection moves each member behind the one before by Δx and Δy, and the nearer hide the farther (WinPLOTR's pseudo-3D): Set Displacement… (Δ) gives the two steps, in degrees and in the y units shown."; return; }
      const pr0 = fetchProj(); void pr0;
      const pr = ps.proj && ps.proj.v && ps.proj.v.projection;
      if (!pr) { note.textContent = hasCore() ? "The core works out the hidden parts…" : "The projection needs the app (the core works out the hidden parts)."; return; }
      // the figure's look (step 6): the x axis over the front member's range, a scale bar of the core's round step
      // (plots::proposed_step on the scale shown) in place of the y axis, whose numbers hold for the front member only
      let span = null, bar = null;
      if (figLook) {
        const front = pr.curves.find((c) => c.member === pr.order[0]);
        const pcs = front ? front.pieces.filter((pc) => pc.length) : [];
        span = pcs.length ? [pcs[0][0][0], pcs[pcs.length - 1][pcs[pcs.length - 1].length - 1][0]] : null;
        bar = scaleBarOf(set, yw);
      }
      plot.setAxes(figLook ? { xLabel: `${xl}, front member`, yLabel: bar ? "" : `${yw}, front member`, yScale: ps.scale, look: "figure", xSpan: span, scaleBar: bar } : { xLabel: `${xl} of the front member, each one behind moved by Δx`, yLabel: `${yw} + each member's step`, yScale: ps.scale });
      const series = [];
      for (const c of pr.curves) {
        const m = set.members[c.member];
        c.pieces.forEach((pc, i) => series.push({ id: `${m.key}:${i}`, kind: "line", label: m.label, legend: false, x: pc.map((p) => p[0]), y: pc.map((p) => p[1]), scaled: true, color: colourOf(ps, m), opacity: opacityOf(ps, m), width: widthOf(ps, m) }));
      }
      // the oblique axis joins the members placed by their temperature; a member without one stands apart (L183), at its
      // measured 2θ, labelled with its name at its own end
      const onAxis = pr.axis.filter((a) => set.members[a[0]].variable);
      const axisLabels = pr.axis.map((a) => { const m = set.members[a[0]]; return { x: a[1][0], y: a[1][1], text: m.variable ? (figLook ? fmt(m.variable.value) : `${fmt(m.variable.value)} ${set.variable_unit}`) : `${m.label} (no ${set.variable_name}${figLook ? "" : ", apart"})`, priority: m.variable ? 0 : 2, tick: figLook && !!m.variable }; });
      // the axis's title over its far end (the figure's look): the variable and its unit, as an axis is labelled
      const far = onAxis.length ? onAxis[onAxis.length - 1] : null;
      if (figLook && far) axisLabels.push({ x: far[1][0], y: far[1][1], text: `${set.variable_name} (${set.variable_unit})`, priority: 5, above: true });
      series.push({ id: "axis", kind: "line", x: onAxis.map((a) => a[1][0]), y: onAxis.map((a) => a[1][1]), scaled: true, color: "--ink", width: 0.8, pointLabels: axisLabels });
      plot.setSeries(series, true);
      const hidden = pr.curves.filter((c) => c.visible < 0.999).length;
      const apartP = ps.proj.v.apart.map((k) => set.members[k].label);
      note.textContent = `${figLook ? `${lookText(bar)}The x axis over the front member's range, the ${set.variable_name} on the oblique axis with its ticks. ` : ""}Each member moved behind the one before (plots::displace) and each curve hiding what lies under it (plots::project): ${hidden} of ${pr.curves.length} curves partly hidden; the oblique axis at the curves' ends carries the ${set.variable_name}. The axis's numbers hold for the front member only.${apartP.length ? ` Apart, at the 2θ where ${apartP.length === 1 ? "it was" : "they were"} measured and off the axis: ${apartP.join(", ")} (no ${set.variable_name}; ${ps.prule.kind === "variable" ? "a rule by temperature does not place it" : "raised by its step, not moved in 2θ"}).` : ""}${mon ? " " + mon + "." : ""}`;
    } else if (ps.view === "series") {
      if (!ps.wins.length) { plot.setSeries([]); note.textContent = "The series view follows windows of 2θ across the set: Add Window… (+ window) names one, from and to; the core gives, for every member, the integral above the background with σ, the centroid with σ and the strongest peak (plots::windows). The windows show as bands in the overlay and the stack."; return; }
      const w0 = fetchWins(); void w0;
      const W = ps.winres && ps.winres.v;
      if (!W) { note.textContent = hasCore() ? "The core integrates the windows…" : "The series needs the app (the core integrates the windows)."; return; }
      const xMembers = W.members.filter((mw) => set.members[mw.member].variable);
      const xs = xMembers.map((mw) => set.members[mw.member].variable.value);
      const series = [];
      if (ps.quantity === "background") {
        ps.at.forEach((x, j) => series.push({ id: "at" + j, kind: "both", points: true, label: `background at ${fmt(x)}°`, x: xs, y: xMembers.map((mw) => (mw.background_at[j] == null ? NaN : mw.background_at[j])), width: 1 }));
      } else {
        const [, qt, get, sig, unit] = QUANTITIES.find(([k]) => k === ps.quantity);
        ps.wins.forEach((w, j) => {
          const ys = xMembers.map((mw) => { const v = get(mw.windows[j]); return v == null ? NaN : v; });
          const ss = xMembers.map((mw) => sig(mw.windows[j]));
          // F5: a point the core finds below 2σ is drawn open (its flag `significant`, PL22); filled where the core says nothing
          const open = xMembers.map((mw) => mw.windows[j].significant === false);
          series.push({ id: "w" + j, kind: "both", points: true, open, label: `${w.name} (${fmt(w.from)}–${fmt(w.to)}°)`, x: xs, y: ys, sigma: ss.some((e) => e != null) ? ss.map((e) => (e == null ? 0 : e)) : null, width: 1, color: WINDOW_COLOURS[j % WINDOW_COLOURS.length], marker: WINDOW_SHAPES[j % WINDOW_SHAPES.length] });
        });
        plot.setAxes({ xLabel: `${set.variable_name} (${set.variable_unit})`, yLabel: `${qt} (${set.normalisation ? unit.replace("counts", "normalised counts") : unit})`, yScale: ps.scale });
      }
      if (ps.quantity === "background") plot.setAxes({ xLabel: `${set.variable_name} (${set.variable_unit})`, yLabel: `the ${ps.wbg === "auto" ? "automatic " : ""}background (${yw})`, yScale: ps.scale });
      plot.setSeries(series);
      const apartVals = W.members.filter((mw) => !set.members[mw.member].variable).map((mw) => `${set.members[mw.member].label}: ${mw.windows.map((wv) => `${wv.name} ${fmt(wv.integral)} ± ${fmt(wv.sigma)}`).join(", ")}`);
      note.textContent = `Each window's ${ps.quantity === "background" ? "background" : (QUANTITIES.find(([k]) => k === ps.quantity) || [])[1]} for every member against its ${set.variable_name} (plots::windows), σ as bars where the core gives one${W.members.some((mw) => mw.windows.some((wv) => wv.significant != null)) ? ", a point below 2σ open (the core's)" : "; a point below 2σ will be drawn open when the core flags it (PL22)"}; the background ${ps.wbg === "auto" ? "is each member's automatic background (GSAS-II's arPLS)" : "is the line between the member's values at the window's ends"}.${inFilesX(ps) ? "" : ` The windows are ranges of the files' x (${members(ps)[0] ? members(ps)[0].axis : "2θ (°)"}), not of ${xl}: the series is taken there.`}${apartVals.length ? ` Drawn apart (no ${set.variable_name}), their integrals: ${apartVals.join("; ")}.` : ""}${mon ? " " + mon + "." : ""}`;
    } else if (surface) {
      // the same numbers as the map (plots::map): x across, each member's row at its temperature in depth, the value up
      // and in colour; the members without a temperature are left out of it, and the note says which
      const diff = ps.surfOf === "difference";
      if (diff && !ps.reference) { dropSurface(); note.textContent = "The surface of the difference map subtracts a reference member from every other: choose one in the sidebar and Set Reference, or show the surface of the map."; return; }
      fetchMap(diff);
      const got = ps.map && ps.map.key.startsWith(diff ? "d:" : "m:") ? ps.map : null;
      if (!got) { note.textContent = hasCore() ? "The core computes the map…" : "The surface needs the app (the core computes its map)."; return; }
      const mv = got.v, M = mv.map;
      const refM = diff ? set.members[got.ref] : null;
      const zWord = zWordOf(got.scale, yw);
      // step 16: the same pair that clips the map's colours clips the surface's height
      const clip = clipOf(ps, diff ? "difference" : "map", got.scale), clipKey = clip ? `${clip.from}|${clip.to}` : "";
      q(".pl-clip-text").textContent = clip ? `${fmt(clip.from)} to ${fmt(clip.to)}` : "the core's range";
      if (!surf) { ensureSurface(); note.textContent = "The surface is being made…"; return; }
      if (surfFor !== mv || surf.xl !== xl || surf.clipKey !== clipKey) {
        surfFor = mv; surf.xl = xl; surf.clipKey = clipKey;
        try { surf.set({
          map: M, scale: got.scale, ramp: diff ? "diverging" : "sequential", clip: clip ? { from: clip.from, to: clip.to } : null,
          labels: { x: xl, v: `${set.variable_name} (${set.variable_unit})`, z: `${diff ? "less the reference, " : ""}${yw}${got.scale === "log" ? ", log scale" : got.scale === "sqrt" ? (diff ? ", signed √ scale" : ", √ scale") : ""}` },
          cellText: (j, r) => { const row = M.rows[r], z = M.z[r][j]; return `${xNameOf(ps)} ${M.x[j].toFixed(4)}  ${set.variable_name} ${fmt(row.value)} ${set.variable_unit} (${set.members[row.member].label})  ${z == null ? "no value" : `${diff ? "less the reference, " : ""}${zWord} ${z.toFixed(3)}`}`; },
        }); } catch (e) { surfFor = null; note.textContent = `The surface could not be drawn: ${e.message || e}.`; status(`The surface could not be drawn: ${e.message || e}`, true); return; }
      }
      const st3 = surf.stats();
      note.textContent = `${diff ? `Every member less ${refM ? refM.label : "the reference"} (the core's difference map), red above zero and blue below. ` : ""}The map as a surface: ${xNameOf(ps)} across, each member's row at its ${set.variable_name} in depth (so that unequal steps show as unequal gaps), the value up and in colour on ${clip ? `the range shown, ${fmt(clip.from)} to ${fmt(clip.to)} ${zWord}, the clip set on the map's colour bar (values beyond it saturated at the floor and the top; the core's range is ${M.range[0].toFixed(3)} to ${M.range[1].toFixed(3)}; Range… or a double-click on the bar in the ${diff ? "difference map" : "map"} returns to it)` : `the core's range, ${M.range[0].toFixed(3)} to ${M.range[1].toFixed(3)} ${zWord}`} (plots::map); each member's pattern drawn over it as a line. ${st3 ? `${st3.rows} rows × ${st3.columns} columns. ` : ""}Drag to turn, the wheel to zoom, a double-click to start again; the pointer gives a cell's value. Export ▸ Plot writes the view as an image; a paper's figure of the same numbers is the map or the projection, through Publish.${mv.apart.length ? ` Not on the surface, no ${set.variable_name} in the file: ${mv.apart.map((a) => set.members[a.member].label).join(", ")} (Set ${set.variable_name === "T" ? "Temperature" : "Value"}… types one).` : ""}${mon ? " " + mon + "." : ""}`;
    } else {
      // step 18: the Map and the Difference Map are the viewer's (steps 16 and 17 inside it); the module fetches the
      // core's map and hands it over with the set's words, and adds what is the app's (the monitor's note)
      const diff = ps.view === "difference";
      if (diff && !ps.reference) { viewer.hide(); plot.setSeries([]); note.textContent = "The difference map subtracts a reference member from every other: choose one in the sidebar and Set Reference."; return; }
      const mp = fetchMap(diff); void mp;
      const got = ps.map && ps.map.key.startsWith(diff ? "d:" : "m:") ? ps.map : null;
      if (!got) { note.textContent = hasCore() ? "The core computes the map…" : "The map needs the app (the core computes it)."; return; }
      if (!diff && ps.cursors) fetchExtent(got.scale);   // the strip's shared y interval, the core's extent on the map's scale
      viewer.show(viewerData(), diff ? "difference" : "map", { look: figLook });
      if (mon) note.textContent += ` ${mon}.`;
    }
  }
  /** The view shown, drawn once more into SVG for File ▸ Export ▸ Plot (L134). */
  function figure() {
    const stem = (tab.title || "plot set").replace(/\.[^.]+$/, "");
    if (isSeq) {
      if (ps.rs && ps.look) { const wpx = Math.round(ps.figmm / 25.4 * 96); return { svg: plot.toSVG({ w: wpx, h: Math.round(wpx * 0.62), ...ps.figfont }), name: `${stem} ${ps.param} figure` }; }
      return ps.seq || ps.rs ? { svg: plot.toSVG(), name: `${stem} ${ps.param}` } : null;
    }
    const set = ps.set; if (!set || !set.members.length) return null;
    if (ps.view === "surface") {
      // a WebGL view is an image: its PNG in an SVG at the canvas's size, so that Export ▸ Plot and Print take it as
      // they take a plot (in the light palette: the surface is drawn again in its colours first)
      if (!surf) return null;
      let im; try { im = surf.image(); } catch (e) { status(`The surface could not be drawn: ${e.message || e}`, true); return null; }
      if (!im.w || !im.h) return null;
      return { svg: `<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="${im.w}" height="${im.h}" viewBox="0 0 ${im.w} ${im.h}"><image x="0" y="0" width="${im.w}" height="${im.h}" xlink:href="${im.url}"/></svg>`, name: `${stem} surface` };
    }
    if (ps.look) {
      // the figure's width (Figure Size…, in mm at 96 px per inch), its text at the preset's sizes in points (PL30 a, the
      // shared plot's `toSVG`), its height in the window's proportion, kept between 0.62 and 1 of the width
      const cv = q(".pl-area canvas"); const wpx = Math.round(ps.figmm / 25.4 * 96);
      const r = cv && cv.clientWidth ? cv.clientHeight / cv.clientWidth : 0.75;
      const hpx = Math.round(wpx * (r < 0.62 ? 0.62 : r > 1 ? 1 : r));
      return { svg: plot.toSVG({ w: wpx, h: hpx, ...ps.figfont }), name: `${stem} ${ps.view} figure` };
    }
    if (ps.view === "overlay" && ps.linked && zplot) {
      const a = { svg: plot.toSVG(), w: q(".pl-area canvas").clientWidth, h: q(".pl-area canvas").clientHeight }, b = { svg: zplot.toSVG(), w: q(".pl-zcanvas").clientWidth, h: q(".pl-zcanvas").clientHeight };
      return { svg: sideBySide([a, b]), name: `${stem} overlay linked` };
    }
    return { svg: plot.toSVG(), name: `${stem} ${ps.view}` };
  }
  views.set(tab.id, { redraw, plot, figure, viewer, viewerData, zoom: () => zplot, surface: () => surf, note: () => note.textContent, hover: (id) => popts.onHover(id == null ? null : { s: { id } }), band: { create: (a, b) => popts.onBandCreate(a, b), change: (i, a, b) => popts.onBandChange("windows", i, a, b), select: (i) => popts.onBandSelect("windows", i) }, cell: (c) => viewer.cell(c), cellClick: (c) => viewer.cellClick(c), strips: () => viewer.strips(), mapShown });
  redraw();
  return () => { dropSurface(); views.delete(tab.id); };
}

// ------------------------------------------------------------------ the commands' actions
/** The view's rule (PL6, PL10) set on the Rust side, which checks it by making the view once; the answers that
    depended on the old one dropped. The members' λs typed one by one are passed by their index now. */
async function applyViewRule(tab, nc, xc, done = null) {
  const ps = stateOf(tab);
  const x = xc ? { ...xc, given: Object.entries(ps.xg).map(([key, [l, n]]) => { const m = memberByKey(ps, key); return m ? [m.index, l, n] : null; }).filter(Boolean) } : null;
  try {
    const s = await invoke("plots_view_rule", { id: tab.setId, norm: nc, x });
    const xChanged = (ps.set.x_quantity || null) !== (s.x_quantity || null);
    ps.nc = nc; ps.xc = xc;
    if (ps.rule.proposed) ps.rule = { kind: "step", dx: 0, dy: null, proposed: null };   // the core proposes the step again, on the y shown
    if (xChanged) ps.zoomwin = null;
    setSet(tab, s);
    const out = s.members.filter((m) => m.shown.left_out).length;
    status(done || `${s.normalisation ? `Normalised: ${normText(ps)}` : "The counts as read"}; ${s.x_quantity ? `x in ${s.x_label}` : "x as in the files"}${out ? `; ${out} of ${s.members.length} members left out of the view (the reason in the note and the members' list)` : ""}. A view: the files untouched.`, out > 0 && out === s.members.length);
    return s;
  } catch (e) { status(String(e), true); return null; }
}
/** Normalise… (PL6): a rule from the core's, its value and its range asked; a view, said in every view. */
async function normaliseView(tab = activePlots()) {
  if (!tab || !hasCore()) return;
  const ps = stateOf(tab); const cur = ps.nc;
  const by = await pickOne("Normalise", "The members as the views show them, normalised by the core (plots::normalise): a view, the files and Data's operations untouched, and every view says it. A member the rule cannot normalise (no monitor, no counting time, no positive area) is left out of the view, with the reason.", NORMS.map(([value, label]) => ({ value, label })), cur ? cur.by : "area", "Next");
  if (!by) return;
  if (by === "none") { await applyViewRule(tab, null, ps.xc, "Normalisation removed: the counts as read."); return; }
  const nc = { by, to: null, from: null, until: null, member: null, background: null };
  if (by === "reference") {
    const c = chosen(ps);
    if (!c.length) { status("Choose the reference member in the sidebar first, then Normalise… ▸ a reference member.", true); return; }
    nc.member = c[0].index;
  }
  const ax = members(ps)[0] ? members(ps)[0].axis : "2θ (°)";
  const same = cur && cur.by === by ? cur : null;
  const w = ps.wins[0];
  const fields = [];
  if (by !== "reference") fields.push({ id: "to", label: "Normalised to", value: same && same.to != null ? same.to : 1, kind: "number" });
  if (by === "maximum" || by === "area" || by === "peak" || by === "reference") {
    const help = by === "maximum" ? "Empty for the largest y over the whole pattern." : by === "peak" && w ? `The first window, ${w.name}, by default.` : "In the files' x, as the windows are.";
    fields.push({ id: "from", label: `From, ${ax}`, value: same && same.from != null ? same.from : by === "peak" && w ? w.from : "", kind: "number", help });
    fields.push({ id: "until", label: `To, ${ax}`, value: same && same.until != null ? same.until : by === "peak" && w ? w.to : "", kind: "number" });
  }
  if (by === "peak") fields.push({ id: "background", label: "Background (auto or line)", value: (same && same.background) || "auto", help: "auto: each member's automatic background (GSAS-II's arPLS); line: the line between the member's values at the range's ends." });
  if (fields.length) {
    const v = await askForm("Normalise", `${NORMS.find(([k]) => k === by)[1]}${by === "reference" ? `: ${chosen(ps)[0].label}` : ""}.`, fields, "Normalise");
    if (!v) return;
    if ([v.to, v.from, v.until].some((x) => x != null && !isFinite(x))) { status("The value and the range are numbers.", true); return; }
    nc.to = v.to == null ? null : v.to; nc.from = v.from == null ? null : v.from; nc.until = v.until == null ? null : v.until;
    if (by === "peak") { const b = String(v.background || "auto").trim().toLowerCase(); if (b !== "auto" && b !== "line") { status(`The background is "auto" or "line", not "${v.background}".`, true); return; } nc.background = b; }
  }
  await applyViewRule(tab, nc, ps.xc);
}
/** Convert X… (PL10): the members in another x quantity, by the core; the members without a λ listed, and a λ for
    them asked: a control file's λ1, or one typed. */
async function convertX(tab = activePlots()) {
  if (!tab || !hasCore()) return;
  const ps = stateOf(tab);
  const q = await pickOne("Convert X", "The members shown in another x quantity, converted by the core (plots::x_view) with each member's λ: the one set for it (Set Wavelength…), else its file's, else a control file's λ1 or one typed for all; a view, the files untouched. A member that needs a λ and has none is left out of the view and listed.", XQS.map(([value, label]) => ({ value, label })), ps.set.x_quantity || "", "Show");
  if (q == null) return;
  if (q === "") { await applyViewRule(tab, ps.nc, null, "x as in the files."); return; }
  const prev = ps.xc || {};
  const xc = { quantity: q, control: prev.control || null, typed: prev.typed == null ? null : prev.typed };
  const s = await applyViewRule(tab, ps.nc, xc); if (!s || !s.without_lambda.length) return;
  // the core lists the members without a λ: one for them, from a control file of the project or typed
  const lams = (state.project && state.project.wavelengths) || [];
  const names = s.without_lambda.map((k) => s.members[k].label);
  const pick = await pickOne("Wavelength", `${names.length} of ${s.members.length} members state no wavelength, so they cannot be shown in ${s.x_label}: ${names.slice(0, 6).join(", ")}${names.length > 6 ? " …" : ""}. Their λ:`, [...lams.map(([rel, l]) => ({ value: "c:" + rel, label: `λ1 = ${fmt(l)} Å, of the control file ${rel}` })), { value: "typed", label: "a λ typed now" }, { value: "none", label: "none: leave them out, listed" }], lams.length ? "c:" + lams[0][0] : "typed", "Show");
  if (!pick || pick === "none") return;
  if (pick.startsWith("c:")) { await applyViewRule(tab, ps.nc, { ...xc, control: `${state.project.root}/${pick.slice(2)}`, typed: null }); return; }
  const v = await askForm("Wavelength", `The λ for the members whose files state none, in ${s.x_label}. A member's own λ is set with Set Wavelength….`, [{ id: "l", label: "λ (Å)", value: prev.typed == null ? "" : prev.typed, kind: "number", placeholder: "1.5406" }], "Show");
  if (!v || v.l == null || !isFinite(v.l)) return;
  await applyViewRule(tab, ps.nc, { ...xc, control: null, typed: v.l });
}
/** Set Wavelength… (PL10): a λ typed for the chosen member, with its note, used by Convert X… before its file's. */
async function setMemberLambda(tab = activePlots()) {
  const ps = stateOf(tab); const c = chosen(ps); if (!c.length || !hasCore()) return;
  const m = c[0]; const had = ps.xg[m.key];
  const v = await askForm(`Wavelength of ${m.label}`, `${m.wavelength == null ? "Its file states none." : `Its file states ${fmt(m.wavelength)} Å.`} The λ typed here is used by Convert X… before the file's, with your note of where it comes from; empty removes it.`, [
    { id: "l", label: "λ (Å)", value: had ? had[0] : "", kind: "number" },
    { id: "note", label: "Its source", value: had ? had[1] : "", placeholder: "e.g. the instrument's log" },
  ], "Set");
  if (!v) return;
  if (v.l == null) delete ps.xg[m.key];
  else { if (!isFinite(v.l)) { status("λ is a number.", true); return; } if (!String(v.note || "").trim()) { status("Say where the λ comes from: it is kept with it.", true); return; } ps.xg[m.key] = [v.l, String(v.note).trim()]; }
  if (ps.xc) await applyViewRule(tab, ps.nc, ps.xc);
  else { status(`${m.label}: λ ${v.l == null ? "removed" : `${fmt(v.l)} Å`}, used when Convert X… shows another quantity.`); renderActive(); }
}
/** Export ▸ Plot Set as CSV… (PL7): the core's CSV, wide or long, in 2θ, d or Q, where the user chooses. */
async function exportCsv(tab = activePlots()) {
  if (!tab || !hasCore()) return;
  const ps = stateOf(tab); const set = ps.set;
  const last = ps.csv || {};
  const v = await askForm("Export the Plot Set as CSV", `The core's CSV (plots::csv, csv_long): wide, x then y and σ of each member on the common grid, or long, one row per member and point; what each column is (its file, scan, format, ${set.variable_name} and its source) in # lines above the header, which Open CSV as a Set… reads back. In d or Q each member is converted with its own λ, or the one given here when its file states none. The members as the view shows them${set.normalisation ? `: normalised (${normText(ps)}), the rule and each member's factor in # lines (PlotSet::as_shown, PL25), the files as read` : ", after Data's operations"}.`, [
    { id: "shape", label: "Shape (wide or long)", value: last.shape || "wide" },
    { id: "x", label: "x (2theta, d or q)", value: last.x || (set.x_quantity === "d" || set.x_quantity === "q" ? set.x_quantity : "2theta") },
    { id: "range", label: "Wide: the grid (overlap or union)", value: last.range || "overlap", help: "overlap: where every member has points; union: every member's points, a cell left empty where one has none." },
    { id: "lambda", label: "λ (Å) where a file states none", value: last.lambda != null ? last.lambda : set.x_fallback ? set.x_fallback[0] : "", kind: "number" },
  ], "Choose File…");
  if (!v) return;
  const shape = String(v.shape).trim().toLowerCase(), x = String(v.x).trim().toLowerCase().replace("2θ", "2theta"), range = String(v.range).trim().toLowerCase();
  if (!["wide", "long"].includes(shape) || !["2theta", "d", "q"].includes(x) || !["overlap", "union"].includes(range)) { status("Shape: wide or long; x: 2theta, d or q; grid: overlap or union.", true); return; }
  if (v.lambda != null && !isFinite(v.lambda)) { status("λ is a number.", true); return; }
  ps.csv = { shape, x, range, lambda: v.lambda };
  const stem = (tab.title || "plot set").replace(/\.[^.]+$/, "");
  const name = `${stem}${x === "2theta" ? "" : ` in ${x === "q" ? "Q" : x}`}${shape === "long" ? " long" : ""}.csv`;
  const path = await pickSavePath({ title: "Export the plot set as CSV", defaultPath: rootOf() ? `${rootOf()}/Exports/${name}` : name, filters: [{ name: "CSV", extensions: ["csv"] }] });
  if (!path) return;
  try {
    const r = await invoke("plots_export_csv", { id: tab.setId, path, long: shape === "long", x, wavelength: v.lambda == null ? null : v.lambda, range });
    status(`${r.path.split("/").pop()}: ${r.members} members, ${r.lines} lines, the core's CSV${set.normalisation ? ", normalised as the view shows it, the rule and the factors in its # lines" : ""}.`);
    bus.emit("file-saved", r.path);
  } catch (e) { status(String(e), true); }
}
/** Open CSV as a Set… (PL8): a CSV of several columns (or several data files) read by the core as a new set, one
    member per y column, each with the file and the temperature its # lines state. */
async function openFilesAsSet() {
  if (!hasCore()) return;
  const paths = await pickFiles({ multiple: true, title: "Open a CSV of several columns (or several data files) as a plot set", defaultPath: rootOf() || undefined, filters: [{ name: "Tables", extensions: ["csv", "tsv", "txt"] }, { name: "Every file", extensions: ["*"] }] });
  if (!paths || !paths.length) return;
  const list = [].concat(paths);
  try { const set = await invoke("plots_open_files", { paths: list, root: rootOf() }); await adopt(set, setTitle(set, list[0].split("/").pop().replace(/\.[^.]+$/, ""))); } catch (e) { status(String(e), true); }
}
async function setDisplacement(tab = activePlots()) {
  if (!tab) return;
  const ps = stateOf(tab); const set = ps.set;
  const proj = ps.view === "projection", r = proj ? ps.prule : ps.rule;   // the projection keeps steps of its own
  const v = await askForm(proj ? "Set the Projection's Steps" : "Set Displacement", `Each member is ${proj ? "moved behind the one before" : "raised"} by the core (plots::displace): by a step per member in the set's order, or in proportion to its ${set.variable_name} from the first member's. Δy in the y units shown (${ps.scale === "linear" ? "counts" : ps.scale === "sqrt" ? "√counts" : "decades"}), Δx in degrees.`, [
    { id: "kind", label: "Rule (step or variable)", value: r.kind, help: `step: the k-th member raised by k·Δy; variable: by (${set.variable_name} − ${set.variable_name}₀)·Δy, so that unequal steps in ${set.variable_name} show as unequal gaps (a member without one is not placed).` },
    { id: "dy", label: "Δy", value: r.dy == null ? "" : r.dy, kind: "number" },
    { id: "dx", label: proj ? "Δx (°)" : "Δx (°, 0 for a plain stack)", value: r.dx || (proj ? "" : 0), kind: "number" },
  ], "Set");
  if (!v) return;
  const kind = String(v.kind).trim().toLowerCase();
  if (kind !== "step" && kind !== "variable") { status(`The rule is "step" or "variable", not "${v.kind}".`, true); return; }
  if (v.dy == null || !isFinite(v.dy) || !isFinite(v.dx == null ? 0 : v.dx)) { status("Δy and Δx are numbers.", true); return; }
  if (proj) { ps.prule = { kind, dy: v.dy, dx: v.dx || 0, scale: ps.scale }; ps.proj = null; } else { ps.rule = { kind, dy: v.dy, dx: v.dx || 0, proposed: null, scale: ps.scale }; ps.disp = null; ps.view = "stack"; }
  renderTab(tab); if (window.__renderPanels) window.__renderPanels();
}
/** Add Window…: a named 2θ range whose integral, centroid and peak the core gives for every member (PL5). */
async function addWindow(tab = activePlots()) {
  if (!tab) return;
  const ps = stateOf(tab);
  const v = await askForm("Add Window", "A named range of 2θ. For every member the core gives the integral above the background with σ, the centroid with σ and the strongest peak in it (plots::windows); the series view draws them against the temperature.", [
    { id: "name", label: "Name", value: `w${ps.wins.length + 1}` },
    { id: "from", label: "From (°)", value: "", kind: "number" },
    { id: "to", label: "To (°)", value: "", kind: "number" },
  ], "Add");
  if (!v) return;
  const name = String(v.name || "").trim();
  if (!name || ps.wins.some((w) => w.name === name)) { status(`A window needs a name of its own${name ? `: "${name}" is taken` : ""}.`, true); return; }
  if (v.from == null || v.to == null || !isFinite(v.from) || !isFinite(v.to) || !(v.to > v.from)) { status("From and to are numbers, to above from.", true); return; }
  ps.wins = [...ps.wins, { name, from: v.from, to: v.to }]; ps.winres = null; ps.view = "series";
  renderActive();
}
/** Draw Window: the next drag across the overlay or the stack draws a window (then it is off again). */
function drawWindow(tab = activePlots()) {
  if (!tab) return;
  const ps = stateOf(tab);
  ps.drawWin = !ps.drawWin;
  if (ps.drawWin && ps.view !== "overlay" && ps.view !== "stack") ps.view = "overlay";
  if (ps.drawWin) ps.look = false;
  status(ps.drawWin ? "Drag across the plot to draw a window; its edges can then be dragged, and a click on a band chooses it." : "Draw Window off.");
  renderActive();
}
/** A window drawn on the plot: named w1, w2 … (the first name not taken), chosen. */
function addWindowRange(tab, a, b) {
  const ps = stateOf(tab);
  let k = ps.wins.length + 1; while (ps.wins.some((w) => w.name === `w${k}`)) k++;
  ps.wins = [...ps.wins, { name: `w${k}`, from: a, to: b }]; ps.winres = null; ps.selwin = ps.wins.length - 1;
  status(`Window w${k}: ${fmt(a)} to ${fmt(b)} (in the files' x). Its edges can be dragged; Rename Window… names it.`);
  renderActive();
}
/** A window's edges dragged on the plot. */
function moveWindow(tab, i, a, b) {
  const ps = stateOf(tab);
  ps.wins = ps.wins.map((w, k) => (k === i ? { ...w, from: a, to: b } : w)); ps.winres = null; ps.selwin = i;
  status(`Window ${ps.wins[i].name}: ${fmt(a)} to ${fmt(b)}; the series view follows.`);
  renderActive();
}
/** Rename Window…: the chosen window's name (asked, or given), unique in the set. */
async function renameWindow(tab = activePlots(), i = null, name = null) {
  if (!tab) return false;
  const ps = stateOf(tab); if (i == null) i = ps.selwin; const w = ps.wins[i]; if (!w) return false;
  if (name == null) { name = await askText("Rename Window", `The window from ${fmt(w.from)} to ${fmt(w.to)}: its name in the series view, the legend and the set file.`, w.name, { label: "Name", ok: "Rename" }); if (name == null) return false; }
  name = String(name).trim();
  if (!name || ps.wins.some((x, k) => k !== i && x.name === name)) { status(`A window needs a name of its own${name ? `: "${name}" is taken` : ""}.`, true); return false; }
  ps.wins = ps.wins.map((x, k) => (k === i ? { ...x, name } : x)); ps.winres = null;
  renderActive(); return true;
}
/** Remove Window: the chosen window only. */
function removeWindow(tab = activePlots(), i = null) {
  if (!tab) return;
  const ps = stateOf(tab); if (i == null) i = ps.selwin; const w = ps.wins[i]; if (!w) return;
  ps.wins = ps.wins.filter((_, k) => k !== i); ps.winres = null; ps.selwin = null;
  if (ps.quantity === "background" && !ps.at.length && !ps.wins.length) ps.quantity = "integral";
  status(`Window ${w.name} removed.`); renderActive();
}
function removeWindows(tab = activePlots()) { const ps = stateOf(tab); ps.wins = []; ps.winres = null; ps.selwin = null; if (ps.quantity === "background" && !ps.at.length) ps.quantity = "integral"; renderActive(); }
async function backgroundAt(tab = activePlots()) {
  if (!tab) return;
  const ps = stateOf(tab);
  const v = await askForm("Background at", "The background the windows are taken above, read by the core at these 2θ (comma-separated), for every member: the series view's \"background at x\".", [{ id: "at", label: "2θ (°)", value: ps.at.join(", "), placeholder: "16.5, 30, 45" }], "Show");
  if (!v) return;
  const at = String(v.at).split(",").map((t) => t.trim()).filter(Boolean).map((t) => parseFloat(t.replace(",", ".")));
  if (at.some((x) => !isFinite(x))) { status("2θ values, separated by commas.", true); return; }
  ps.at = at; ps.winres = null; if (at.length) { ps.quantity = "background"; ps.view = "series"; }
  renderActive();
}
// ------------------------------------------------------------------ step 16: the clip and the cursors (the viewer's, step 18)
/** The tab's viewer (step 18), or null before its view is made. */
const viewerOf = (tab) => { const v = tab && views.get(tab.id); return v && v.viewer ? v.viewer : null; };
/** The clip of a view set (a drag of the bar's marks, `quiet`, or typed), the viewer's; the surface follows the state. */
function setClip(tab, key, from, to, scale, quiet = false) { const v = viewerOf(tab); if (v) v.setClip(key, from, to, scale, quiet); }
/** Back to the core's range (a double-click on the bar, Reset Colour Range, or both fields empty). */
function resetClip(tab, key) { const v = viewerOf(tab); if (v) v.resetClip(key); }
/** The view whose clip a command means: the difference map's in that view (or its surface), else the map's. */
const clipKeyOf = (ps) => (ps.view === "difference" || (ps.view === "surface" && ps.surfOf === "difference") ? "difference" : "map");
/** Set Colour Range…: the clip typed, on the scale the map is drawn on, inside the core's range (I08 step 1). */
async function setColourRange(tab = activePlots()) {
  if (!tab) return;
  const ps = stateOf(tab); const key = clipKeyOf(ps);
  const got = ps.map && ps.map.key.startsWith(key === "difference" ? "d:" : "m:") ? ps.map : null;
  if (!got) { status("Show the map, the difference map or the surface first: the range is on its colour scale.", true); return; }
  const M = got.v.map, cur = clipOf(ps, key, got.scale), what = key === "difference" ? "difference map" : "map";
  const v = await askForm("Set Colour Range", `The range the ${what} and its surface show, on the ${got.scale} colour scale they are drawn on, inside the core's range ${M.range[0].toFixed(3)} to ${M.range[1].toFixed(3)}: values below it are drawn in the first colour, above it in the last (saturated), so that the weak details come up as the strong peaks are clipped. A field left empty keeps the core's end; both empty, the core's range. The two marks on the colour bar drag the same pair.`, [
    { id: "from", label: "From", value: cur ? cur.from : "", kind: "number" },
    { id: "to", label: "To", value: cur ? cur.to : "", kind: "number" },
  ], "Set");
  if (!v) return;
  // the pair checked and set by the viewer (inside the core's range, in order; empty fields the core's ends); in the
  // surface view the viewer holds the map's data from its last showing, so the pair is kept in the state directly
  const vw = viewerOf(tab);
  if (vw && vw.data() && vw.data().maps[key]) { vw.typedClip(key, v.from, v.to); return; }
  if (v.from == null && v.to == null) { if (cur) { ps.clip[key] = null; renderTab(tab); } else status("The core's range, as before."); return; }
  const from = v.from == null ? M.range[0] : v.from, to = v.to == null ? M.range[1] : v.to;
  if (!isFinite(from) || !isFinite(to)) { status("From and to are numbers.", true); return; }
  if (!(to > from)) { status("To lies above from.", true); return; }
  if (to <= M.range[0] || from >= M.range[1]) { status(`The pair ${fmt(from)} to ${fmt(to)} lies wholly outside the core's range ${M.range[0].toFixed(3)} to ${M.range[1].toFixed(3)}: nothing of the ${what} would be seen.`, true); return; }
  ps.clip[key] = { from, to, scale: got.scale }; renderTab(tab);
  status(`${what === "difference map" ? "The difference map" : "The map"} shows ${fmt(from)} to ${fmt(to)} on its ${got.scale} colour scale; values beyond are saturated.`);
}
// ------------------------------------------------------------------ step 17: the zoom's range and the regions (PL37)
/** The map a tab shows now and what its core call needs: the view's key, the scale and the reference, or null. */
function mapCall(tab) {
  const ps = stateOf(tab); const vw = views.get(tab.id); const shown = vw && vw.mapShown ? vw.mapShown() : null;
  if (!shown) return null;
  return { key: shown.diff ? "difference" : "map", scale: shown.got.scale, reference: shown.diff ? shown.got.ref : null, M: shown.got.v.map, diff: shown.diff };
}
// The zoom's range and the regions are the viewer's (step 18): the module's words below hand each command to the
// tab's viewer, after the view the command needs is shown, and mark the tab dirty where the set file changes.
/** "Range from the view": the viewer asks the core's block of the rectangle shown through the host. */
function followRange(tab) { const v = viewerOf(tab); if (v) return v.followRange(); }
/** Range from the View on or off. */
function toggleFollow(tab = activePlots()) {
  if (!tab) return;
  if (!mapCall(tab)) { status("Show the map or the difference map first: the range follows its zoom.", true); return; }
  viewerOf(tab).toggleFollow();
}
/** A region added to the tab after the core accepts it on the map shown. */
async function addRegion(tab, name, rect, clip = null) {
  const v = viewerOf(tab); if (!v) return null;
  const r = await v.addRegion(name, rect, clip); if (r) tab.dirty = true; return r;
}
/** Save Region…: the rectangle shown named (the viewer asks the name), with the clip when Range from the View is on. */
async function saveRegion(tab = activePlots()) {
  if (!tab) return null;
  if (!mapCall(tab)) { status("Show the map or the difference map first: a region is the rectangle it shows.", true); return null; }
  const r = await viewerOf(tab).saveRegion(); if (r) tab.dirty = true; return r;
}
/** A region chosen: the map shown if it is not, then the viewer zooms to the core's block and applies its range. */
async function chooseRegion(tab, i) {
  const ps = stateOf(tab); if (!ps.regions[i]) return false;
  if (ps.view !== "map" && ps.view !== "difference") { ps.view = "map"; renderActive(); }
  const v = viewerOf(tab); return v ? v.chooseRegion(i) : false;
}
async function renameRegion(tab = activePlots(), i = null, name = null) {
  const v = viewerOf(tab); if (!v) return false;
  const ok = await v.renameRegion(i, name); if (ok) tab.dirty = true; return ok;
}
function removeRegion(tab = activePlots(), i = null) { const v = viewerOf(tab); if (v) { v.removeRegion(i); tab.dirty = true; } }
// ------------------------------------------------------------------ step 22 (W9 c): a large set rebinned on opening
/** How large a set is before the rebin is offered: members × points, the first member's points standing for all. */
const LARGE_POINTS = 2000000;
const isLarge = (set) => !!(set && set.members.length && set.members.length * (set.members[0].points || 0) > LARGE_POINTS);
/** A large set just opened: the core's rebin offered ("101 scans of 27 000 points: open rebinned to 0.01°, or as
    measured"), the choice kept in the tab (`ps.rebin`) and said in the note. The step is the one the ESRF case was
    worked on (0.01°, W9); the core's own proposal comes with 0.81.0, and the field is read when it is there
    (`set.rebin_proposal`). Nothing is asked twice: a set answered keeps its choice. */
async function offerRebin(tab) {
  const ps = stateOf(tab); const set = ps.set;
  if (!set || !isLarge(set) || ps.rebin !== undefined || !hasCore()) return;
  const n = set.members.length, pts = set.members[0].points, step = set.rebin_proposal != null ? set.rebin_proposal : 0.01;
  const v = await pickOne("A Large Set", `${n} scans of ${pts} points: open rebinned to ${fmt(step)}° (the core's rebin, each new bin the average of its points; ${set.rebin_proposal != null ? "the core's proposal" : "the step the ESRF case was worked on; the core's proposal comes with 0.81.0"}), or as measured? Rebinned, the views draw in a moment; as measured, the overlay of every scan is drawn in slices and takes longer. The choice is kept with this tab (Reload Members keeps it) and said in the note.`, [{ value: "rebin", label: `rebinned to ${fmt(step)}°` }, { value: "raw", label: "as measured" }], "rebin", "Open");
  if (v == null) { ps.rebin = null; return; }
  if (v === "raw") { ps.rebin = null; renderActive(); return; }
  try {
    status(`Rebinning ${n} scans to ${fmt(step)}° (the core)…`);
    const s = await invoke("plots_rebin", { id: tab.setId, step });
    ps.rebin = step; ps.patterns = null; setSet(tab, s);
    status(`${n} scans rebinned to ${fmt(step)}° by the core: ${s.members[0] ? s.members[0].points : "?"} points each.`);
  } catch (e) { ps.rebin = null; status(`Not rebinned: ${e}`, true); }
}
/** What the note says of the rebin, in every view. */
const rebinText = (ps) => (ps.rebin ? ` Rebinned to ${fmt(ps.rebin)}° by the core on opening (each new bin the average of its points); as measured on Reload Members after Rebin… ▸ as measured.` : "");
/** Rebin…: the choice made again (the core's rebin to a step, or as measured). */
async function rebinSet(tab = activePlots()) {
  if (!tab || !hasCore()) return;
  const ps = stateOf(tab); const set = ps.set; if (!set) return;
  const v = await askForm("Rebin", `The members read again by the core rebinned to a step (each new bin the average of its points), or as measured (empty). Now: ${ps.rebin ? `rebinned to ${fmt(ps.rebin)}°` : "as measured"}.`, [{ id: "step", label: "Step (°)", value: ps.rebin || "", kind: "number" }], "Rebin");
  if (!v) return;
  const step = v.step == null ? null : v.step;
  if (step != null && !(isFinite(step) && step > 0)) { status("The step is a positive number, or empty for the scans as measured.", true); return; }
  try { const s = await invoke("plots_rebin", { id: tab.setId, step }); ps.rebin = step; ps.patterns = null; setSet(tab, s); status(step ? `${set.members.length} members rebinned to ${fmt(step)}° by the core.` : "The members as measured again."); } catch (e) { status(String(e), true); }
}

// ------------------------------------------------------------------ step 20 (W4): the y interval, the members' list
/** The y interval shared by the whole set or adapted to the member chosen (Overlay, Stack) and to the row under the
    pointer (the strip under the map); kept in the tab (the set file has no field for it: PL41). */
function toggleYInterval(tab = activePlots()) {
  if (!tab) return;
  const ps = stateOf(tab); ps.yInterval = ps.yInterval === "own" ? "shared" : "own";
  status(ps.yInterval === "shared" ? "The y interval shared by the whole set: the curves fitted to every member, the strip under the map on the set's extent (the core's). Kept in the tab (PL41: the set file has no field for it yet)." : `The y interval adapted: Overlay and Stack fitted to the member${chosen(ps).length === 1 ? "" : "s"} chosen${chosen(ps).length ? ` (${chosen(ps).map((m) => m.label).join(", ")})` : " (choose one in the sidebar)"}, the strip under the map to the row under the pointer.`);
  renderActive();
}
/** A member's facts, the words shown on hover: its file, variable with its source, kind, points, counting time, the
    file's time, Data's operations, and what the view did to it. */
function memberFacts(ps, m) {
  const when = m.modified == null ? "" : `\nfile time ${new Date(m.modified * 1000).toISOString().replace("T", " ").slice(0, 16)}`;
  const others = (m.others || []).map((o) => `\n${o.name} = ${fmt(o.value)}, from ${o.source}`).join("");
  return `${m.file}${m.scans > 1 ? ` (scan ${m.scan + 1} of ${m.scans})` : ""}\n${variableText(ps, m, true)}${others}\n${m.format_name} (${kindFamily(m.format, m.format_name).words}), ${m.points} points, ${fmt(m.range[0])} to ${fmt(m.range[1])}${m.wavelength == null ? "" : `, λ ${fmt(m.wavelength)} Å`}${m.counting_time == null ? "\ncounting time not stated" : `\ncounting time ${fmt(m.counting_time)} s (its exposure class)`}${m.monitor == null ? "\nno monitor: raw counts" : `\nmonitor ${fmt(m.monitor)}`}${when}${m.operations.length ? `\ndrawn after Data's operations: ${m.operations.join("; ")}` : "\nas measured (no operation in Data)"}${m.shown.left_out ? `\nleft out of the view: ${m.shown.left_out}` : ""}${m.shown.factor != null ? `\nnormalised: × ${fmt(m.shown.factor)}` : ""}${m.shown.lambda != null ? `\nλ ${fmt(m.shown.lambda)} Å from ${m.shown.lambda_from}` : ""}`;
}
/** One member chosen or not by its checkbox (its key in the tab's choice, in the order clicked). */
function pickOneMember(tab, key, on) {
  const ps = stateOf(tab); const sel = ps.sel.filter((k) => k !== key); ps.sel = on ? [...sel, key] : sel;
  if (window.__renderPanels) window.__renderPanels(); renderBottomSelection(); renderTab(tab);
}
/** Every member chosen, or none. */
function pickAll(tab = activePlots(), on = true) {
  if (!tab) return; const ps = stateOf(tab);
  ps.sel = on ? members(ps).map((m) => m.key) : [];
  status(on ? `Every member chosen (${ps.sel.length}).` : "No member chosen.");
  if (window.__renderPanels) window.__renderPanels(); renderBottomSelection(); renderTab(tab);
}
/** The members of one exposure class: the counting time the files state (W2's class by the esd's median is the
    core's, later). */
async function pickClass(tab = activePlots(), which = null) {
  if (!tab) return; const ps = stateOf(tab);
  const classes = exposureClasses(members(ps).map((m) => ({ counting_time: m.counting_time })));
  if (!classes.length) return;
  if (which == null) {
    const v = await pickOne("Choose by Exposure Class", "The members whose files state the same counting time (the exposure class as the files give it; the core's class by the esd's median comes with W2).", classes.map((c, i) => ({ value: String(i), label: `${c.label}: ${c.members.length} member${c.members.length === 1 ? "" : "s"}` })), "0", "Choose");
    if (v == null) return; which = parseInt(v, 10);
  }
  const c = classes[which]; if (!c) return;
  ps.sel = c.members.map((k) => members(ps)[k].key);
  status(`${ps.sel.length} member${ps.sel.length === 1 ? "" : "s"} chosen: counting time ${c.label}.`);
  if (window.__renderPanels) window.__renderPanels(); renderBottomSelection(); renderTab(tab);
}
/** The members whose variable lies in a range, typed. */
async function pickRange(tab = activePlots(), from = null, to = null) {
  if (!tab) return; const ps = stateOf(tab); const set = ps.set;
  if (from == null || to == null) {
    const v = await askForm("Choose by a Range of the Variable", `The members whose ${set.variable_name} lies from a value to another (ends included); a member without one is never chosen.`, [{ id: "from", label: `From (${set.variable_unit})`, value: "", kind: "number" }, { id: "to", label: `To (${set.variable_unit})`, value: "", kind: "number" }], "Choose");
    if (!v || v.from == null || v.to == null || !isFinite(v.from) || !isFinite(v.to)) { if (v) status("From and to are numbers.", true); return; }
    from = v.from; to = v.to;
  }
  const idx = membersInRange(members(ps).map((m) => ({ variable: m.variable })), from, to);
  ps.sel = idx.map((k) => members(ps)[k].key);
  status(`${ps.sel.length} member${ps.sel.length === 1 ? "" : "s"} chosen: ${set.variable_name} from ${fmt(from)} to ${fmt(to)} ${set.variable_unit}.`);
  if (window.__renderPanels) window.__renderPanels(); renderBottomSelection(); renderTab(tab);
}
/** The members' list order (the sidebar's and the page's): the set's own, by the variable, by the file's time, by name. */
async function setListOrder(tab = activePlots(), by = null) {
  if (!tab) return; const ps = stateOf(tab);
  if (by == null) { by = await pickOne("List Order", "The order of the members' list in the sidebar (the set's own order, the core's, is not changed: Sort by Variable and Reverse Order do that).", [{ value: "set", label: "the set's order" }, { value: "variable", label: `by ${ps.set.variable_name} (those without one last)` }, { value: "time", label: "by the file's time (those without one last)" }, { value: "name", label: "by name" }], ps.sortBy || "set", "Order"); if (by == null) return; }
  ps.sortBy = by; if (window.__renderPanels) window.__renderPanels();
}
/** Row and Column Cursors: on or off, in the map and the difference map (shown if another view is). */
function toggleCursors(tab = activePlots()) {
  if (!tab) return;
  const ps = stateOf(tab);
  if (!ps.cursors && ps.view !== "map" && ps.view !== "difference") { ps.view = "map"; renderActive(); }
  const v = viewerOf(tab); if (v) v.toggleCursors();
}
/** Export ▸ Viewer Page… (Dev; step 18): the viewer's page written from this tab's data, as the core's `drs plots show`
    will write it (PL38): the page shell, the JSON of the set's maps (the map and, with a reference, the difference map,
    each saved region's block from the core), and the viewer bundled from the same sources the node bundler reads, all
    fetched from the app itself; nothing from the network. Returns { path, html, bytes }. */
async function exportViewerPage(tab = activePlots(), path = null) {
  if (!tab || !hasCore() || !stateOf(tab).set) return null;
  const ps = stateOf(tab); const set = ps.set; const root = rootOf();
  const stem = (tab.title || "plot set").replace(/\.[^.]+$/, "");
  if (!path) path = await pickSavePath({ title: "Write the Map viewer's page (the set's data and the viewer inside, nothing fetched)", defaultPath: root ? `${root}/Exports/${stem} map viewer.html` : `${stem} map viewer.html`, filters: [{ name: "HTML page", extensions: ["html"] }] });
  if (!path) return null;
  try {
    const ref = ps.reference ? members(ps).findIndex((m) => m.key === ps.reference) : -1;
    const mapV = await invoke("plots_map", { id: tab.setId, scale: ps.mapScale, reference: null });
    const diffV = ref >= 0 ? await invoke("plots_map", { id: tab.setId, scale: diffScaleOf(ps.mapScale), reference: ref }).catch(() => null) : null;
    const maps = { map: { map: mapV.map, apart: mapV.apart, scale: ps.mapScale, reference: null }, difference: diffV ? { map: diffV.map, apart: diffV.apart, scale: diffScaleOf(ps.mapScale), reference: ref } : null };
    const regions = [];
    for (const r of ps.regions) {
      const blocks = { map: null, difference: null }; let refused = null;
      for (const key of ["map", "difference"]) {
        if (!maps[key]) continue;
        try { blocks[key] = await invoke("plots_map_region", { id: tab.setId, scale: maps[key].scale, reference: maps[key].reference, region: r }); } catch (e) { refused = String(e); }
      }
      regions.push({ ...r, blocks, refused });
    }
    const mapClip = clipOf(ps, "map", ps.mapScale), diffClip = clipOf(ps, "difference", diffScaleOf(ps.mapScale));
    // the members' patterns for the page (W9 a): asked of the core here, the set view carrying none
    const pageList = await invoke("plots_patterns", { id: tab.setId, members: null });
    const pagePats = new Map(pageList.map((p) => [p.member, { x: p.x, y: p.y }]));
    const data = {
      title: tab.title, path: tab.drv || tab.path || null, core: null, written: `written by the app on ${new Date().toISOString().slice(0, 10)}`,
      variable: { name: set.variable_name, unit: set.variable_unit }, x: { label: xLabelOf(ps), quantity: xWordOf(ps) }, y_word: yWord(ps),
      members: set.members.map((m) => ({ label: m.label, variable: m.variable ? { value: m.variable.value, source: m.variable.source } : null, pattern: m.shown.left_out ? null : patternOf(pagePats, m) })),
      maps, clip: { map: mapClip ? { from: mapClip.from, to: mapClip.to } : null, difference: diffClip ? { from: diffClip.from, to: diffClip.to } : null },
      regions, regions_refused: ps.regionsRefused, same_grid: set.same_grid,
    };
    // the sources as the node bundler reads them, from the app's own files; the hash as SOURCES.txt names it
    const sha = async (text) => [...new Uint8Array(await crypto.subtle.digest("SHA-256", new TextEncoder().encode(text)))].map((b) => b.toString(16).padStart(2, "0")).join("");
    const sources = [];
    for (const name of VIEWER_SOURCES) { const text = await (await fetch(name)).text(); sources.push({ name, text, sha256: await sha(text) }); }
    const shell = await (await fetch("viewer/page.html")).text();
    const safe = (t) => t.replace(/<\//g, "<\\/");
    const html = shell.replace("<!--DRS_PAGE_DATA-->", safe(JSON.stringify(data))).replace("<!--DRS_VIEWER_BUNDLE-->", safe(bundleViewer(sources)));
    const bytes = await invoke("plots_write_text", { path, text: html });
    status(`${path.split("/").pop()}: the Map viewer's page, ${Math.round(bytes / 1024)} kB, the set's data and the viewer inside (${regions.length} region(s)${diffV ? ", the difference map too" : ""}); open it in a browser. The core's own page is drs plots show (PL38).`);
    bus.emit("file-saved", path);
    return { path, html, bytes };
  } catch (e) { status(`The viewer's page could not be written: ${e}`, true); return null; }
}
function chosen(ps) { return ps.sel.map((k) => memberByKey(ps, k)).filter(Boolean); }
/** The two members of the overlay's difference (F6): the first chosen (or the set's first) less the reference (or the
    second chosen); their indices in the set, or null. */
function diffPair(ps) {
  const ms = members(ps); const c = chosen(ps);
  const a = c.length ? c[0] : ms[0];
  const b = c.length > 1 ? c[1] : memberByKey(ps, ps.reference);
  if (!a || !b || a.key === b.key) return null;
  return [ms.indexOf(a), ms.indexOf(b)];
}
async function reorder(cmd, tab = activePlots()) {
  if (!tab || !hasCore()) return;
  try { const s = await invoke(cmd, { id: tab.setId }); stateOf(tab).order = cmd === "plots_reverse" ? "given" : "variable"; setSet(tab, s); status(cmd === "plots_reverse" ? "The members in the reverse order: the stack, the projection and the colours follow." : "The members sorted by temperature (the core's), those without one last."); } catch (e) { status(String(e), true); }
}
/** Read Variable From… (L183): each member's temperature from one of the sources the core found in its file. */
async function readVariable(tab = activePlots()) {
  if (!tab || !hasCore()) return;
  const ps = stateOf(tab); const set = ps.set;
  const choice = await pickOne(`Read ${set.variable_name} From`, `Each member's ${set.variable_name} is taken from the source chosen, where the core finds it in the file; a member whose file says nothing there keeps its value, and a value you typed is kept.`, [
    { value: "name", label: "the file's name (\"_300K\", \"_T300\")" },
    { value: "first", label: "the file itself, as opened (its header, the buffer's comment, its title, its name)" },
    { value: "tempfile", label: "a temperatures file (FullProf's: a file name and TSAMP per line)" },
  ], "name", "Read");
  if (!choice) return;
  if (choice === "tempfile") { await readTemperatures(tab); return; }
  try {
    const [s, n] = await invoke("plots_read_variable", { id: tab.setId, source: choice, replaceTyped: false });
    setSet(tab, s);
    status(n ? `${n} member${n === 1 ? "" : "s"} took ${set.variable_name} from ${choice === "name" ? "the file's name" : "the file"}; the others kept theirs.` : `No file says its ${set.variable_name} ${choice === "name" ? "in its name" : "itself"}: every member kept its value.`, !n);
  } catch (e) { status(String(e), true); }
}
function solo(tab = activePlots()) {
  const ps = stateOf(tab); const c = chosen(ps); if (!c.length) return;
  for (const m of members(ps)) ps.opacity[m.key] = m.key === c[0].key ? 1 : DIMMED;
  ps.soloKey = c[0].key;
  renderActive();
}
/** ⌥-click on a member (its curve or its row, DESIGN.md §4): that member alone at full opacity, the others dimmed; the
    same ⌥-click again brings every member back. */
function soloToggle(tab, key) {
  const ps = stateOf(tab); const m = memberByKey(ps, key); if (!m) return;
  if (ps.soloKey === key) { ps.opacity = {}; ps.soloKey = null; status(`Every member back at full opacity.`); }
  else { for (const x of members(ps)) ps.opacity[x.key] = x.key === key ? 1 : DIMMED; ps.soloKey = key; ps.sel = [key]; status(`${m.label} alone, the others dimmed; ⌥-click it again to bring them back.`); }
  renderActive();
}
function dimOthers(tab = activePlots()) {
  const ps = stateOf(tab); const keep = new Set(ps.sel);
  for (const m of members(ps)) ps.opacity[m.key] = keep.has(m.key) ? 1 : DIMMED;
  renderActive();
}
function restoreAll(tab = activePlots()) { const ps = stateOf(tab); ps.opacity = {}; ps.soloKey = null; renderActive(); }
function setReference(tab = activePlots()) {
  const ps = stateOf(tab); const c = chosen(ps); if (!c.length) return;
  ps.reference = c[0].key; ps.map = null; ps.view = "difference";
  status(`Reference of the difference map: ${c[0].label}.`); renderActive();
}
async function setVariable(tab = activePlots()) {
  const ps = stateOf(tab); const c = chosen(ps); if (!c.length || !hasCore()) return;
  const m = c[0], set = ps.set;
  const v = await askForm(`Set ${set.variable_name} of ${m.label}`, `${m.variable ? `Now ${variableText(ps, m, true)}.` : `The file states no ${set.variable_name}${m.title ? ` (its title: "${m.title}")` : ""}.`} A typed value keeps your note of where it came from (L169).`, [
    { id: "value", label: `${set.variable_name} (${set.variable_unit})`, value: m.variable ? m.variable.value : "", kind: "number" },
    { id: "note", label: "Its source", value: m.variable && m.variable.typed ? m.variable.source.replace(/^typed: /, "") : "", placeholder: "e.g. the experiment's log; WinPLOTR's demo page", help: "Where the value comes from, kept with it and written into what is made from the set." },
  ], "Set");
  if (!v) return;
  if (v.value == null || !isFinite(v.value)) { status(`${set.variable_name} is a number.`, true); return; }
  if (!String(v.note || "").trim()) { status("Say where the value comes from: it is kept with it (L169).", true); return; }
  try { const s = await invoke("plots_set_variable", { id: tab.setId, member: m.index, value: v.value, note: String(v.note).trim() }); setSet(tab, s); status(`${m.label}: ${set.variable_name} = ${fmt(v.value)} ${set.variable_unit}, typed (${String(v.note).trim()}).`); } catch (e) { status(String(e), true); }
}
async function derive(kind, tab = activePlots()) {
  const ps = stateOf(tab); const c = chosen(ps);
  if (!state.project || !state.project.root || !writable()) { status("A derived pattern is written into the project's Plots/ folder: open a project folder of your own (not a reference folder) first.", true); return; }
  const name = (tab.title || "plot set").replace(/\.[^.]+$/, "");
  try {
    const d = await invoke("plots_derive", { id: tab.setId, members: c.map((m) => m.index), kind, root: state.project.root, name });
    status(`The ${d.derivation} of ${c.length} member${c.length === 1 ? "" : "s"} written to ${d.relative}: ${d.points} points${d.dropped ? `, ${d.dropped} grid points left out where a member has none` : ""}${ps.set.normalisation ? "; normalised as the view shows it, the rule and the factors in its ! lines (PL25)" : ""}.`);
    bus.emit("file-saved", d.path);
  } catch (e) { status(String(e), true); }
}
async function openInData(tab = activePlots()) {
  const ps = stateOf(tab); const c = chosen(ps); if (!c.length || !hasCore()) return;
  const m = c[0], root = state.project && state.project.root;
  if (root && m.path.startsWith(root + "/") && window.__openPath) { await window.__openPath(m.path.slice(root.length + 1), "data"); return; }
  try { await invoke("window_open", { kind: "data", path: m.path }); } catch (e) { status(String(e), true); }
}
function select(tab, key, add) {
  const ps = stateOf(tab);
  if (add) { const i = ps.sel.indexOf(key); if (i >= 0) ps.sel.splice(i, 1); else ps.sel.push(key); }
  else ps.sel = [key];
  if (window.__renderPanels) window.__renderPanels(); renderBottomSelection();
}

// ------------------------------------------------------------------ the bottom panel: the members as a list
let table = null;
function renderBottom() {
  const box = $("pl-table"); if (!box) return;
  const t = activePlots(); const ps = t ? stateOf(t) : null;
  if (ps && ps.seq) { renderSeqTable(box, ps); return; }
  if (ps && ps.rs) { renderResultsTable(box, ps); return; }
  if (box.dataset.kind === "seq" || box.dataset.kind === "results") { box.dataset.kind = ""; table = null; }
  if (!ps || !ps.set || !ps.set.members.length) { table = null; box.innerHTML = `<div class="empty">The members of the plot set, one row each: the file, its format, points and range, its wavelength and monitor, its ${ps && ps.set ? ps.set.variable_name : "T"} and where it was found.</div>`; return; }
  const set = ps.set;
  const rows = set.members.map((m) => ({ m, key: m.key, label: m.label, value: m.variable ? m.variable.value : null, source: m.variable ? m.variable.source : "none: drawn apart", format: m.format_name, points: m.points, from: m.range[0], to: m.range[1], step: m.step, lambda: m.wavelength, monitor: m.monitor, ops: m.operations.length ? m.operations.join("; ") : "none: as measured", time: m.counting_time, title: m.title || "", factor: m.shown.factor, shownl: m.shown.lambda, inview: m.shown.left_out ? `left out: ${m.shown.left_out}` : m.shown.lambda != null ? `λ ${fmt(m.shown.lambda)} Å from ${m.shown.lambda_from}` : "shown" }));
  if (!table || !box.querySelector(".dtable")) {
    box.innerHTML = "";
    table = createTable(box, [
      { id: "label", title: "member", get: (r) => r.label, num: false, width: 160 },
      { id: "value", title: `${set.variable_name} (${set.variable_unit})`, get: (r) => r.value, fmt: (v) => (v == null ? "—" : fmt(v)) },
      { id: "source", title: "from", get: (r) => r.source, num: false, width: 220 },
      { id: "format", title: "format", get: (r) => r.format, num: false, width: 200 },
      { id: "points", title: "points", get: (r) => r.points, fmt: (v) => String(v) },
      { id: "from", title: "from (°)", get: (r) => r.from, fmt: (v) => fmt(v) },
      { id: "to", title: "to (°)", get: (r) => r.to, fmt: (v) => fmt(v) },
      { id: "step", title: "step", get: (r) => r.step, fmt: (v) => (v == null ? "" : fmt(v)) },
      { id: "lambda", title: "λ (Å)", get: (r) => r.lambda, fmt: (v) => (v == null ? "not stated" : fmt(v)) },
      { id: "monitor", title: "monitor", get: (r) => r.monitor, fmt: (v) => (v == null ? "none: raw counts" : fmt(v)) },
      { id: "ops", title: "Data's operations", get: (r) => r.ops, num: false, width: 220 },
      { id: "factor", title: "factor (view)", get: (r) => r.factor, fmt: (v) => (v == null ? "" : fmt(v)) },
      { id: "inview", title: "in the view", get: (r) => r.inview, num: false, width: 240 },
      { id: "time", title: "time (s)", get: (r) => r.time, fmt: (v) => (v == null ? "" : fmt(v)), shown: false },
      { id: "title", title: "title", get: (r) => r.title, num: false, width: 280, shown: false },
    ], { onSelect: (idx) => { const t2 = activePlots(); if (!t2 || !idx.length) return; stateOf(t2).sel = idx.map((i) => rows[i].key); if (window.__renderPanels) window.__renderPanels(); } });
  }
  table.setRows(rows);
  renderBottomSelection(rows);
}
/** A results series' refinements as rows: the file, its temperature and where it came from, χ², R_wp, the data file. */
function renderResultsTable(box, ps) {
  const rows = ps.rs.points.map((p) => ({ label: p.label, t: p.variable ? p.variable.value : null, from: p.variable ? p.variable.source : p.no_data || "none: not drawn", chi2: p.chi2, rwp: p.rwp, data: p.data_file ? p.data_file.split("/").pop() : "" }));
  if (box.dataset.kind !== "results" || !table) {
    box.innerHTML = ""; box.dataset.kind = "results";
    table = createTable(box, [
      { id: "label", title: "refinement (results file)", get: (r) => r.label, num: false, width: 240 },
      { id: "t", title: "T (K)", get: (r) => r.t, fmt: (v) => (v == null ? "—" : fmt(v)) },
      { id: "from", title: "from", get: (r) => r.from, num: false, width: 240 },
      { id: "chi2", title: "χ²", get: (r) => r.chi2, fmt: (v) => fmt(v) },
      { id: "rwp", title: "R_wp", get: (r) => r.rwp, fmt: (v) => fmt(v) },
      { id: "data", title: "data file", get: (r) => r.data, num: false, width: 200 },
    ], {});
  }
  table.setRows(rows);
}
/** A `.seq`'s refinements as rows: the data file's code, PARVAL's index, external parameter and numor, χ². */
function renderSeqTable(box, ps) {
  const rows = ps.seq.blocks.map((b) => ({ name: b.name, index: b.index, parameter: b.parameter, numor: b.numor, chi2: b.chi2, n: b.parameters.length }));
  if (box.dataset.kind !== "seq" || !table) {
    box.innerHTML = ""; box.dataset.kind = "seq";
    table = createTable(box, [
      { id: "name", title: "refinement (data file)", get: (r) => r.name, num: false, width: 160 },
      { id: "index", title: "PARVAL index", get: (r) => r.index, fmt: (v) => (v == null ? "" : String(v)) },
      { id: "parameter", title: "PARVAL parameter", get: (r) => r.parameter, fmt: (v) => (v == null ? "" : fmt(v)) },
      { id: "numor", title: "numor", get: (r) => r.numor, fmt: (v) => (v == null ? "" : String(v)) },
      { id: "chi2", title: "χ²", get: (r) => r.chi2, fmt: (v) => (v == null ? "" : fmt(v)) },
      { id: "n", title: "parameters", get: (r) => r.n, fmt: (v) => String(v) },
    ], {});
  }
  table.setRows(rows);
}
function renderBottomSelection(rows) {
  const t = activePlots(); if (!t || !table || stateOf(t).seq) return;
  const ps = stateOf(t); const keys = (rows || ps.set.members.map((m) => ({ key: m.key }))).map((r) => r.key);
  table.setSelection(ps.sel.map((k) => keys.indexOf(k)).filter((i) => i >= 0), false);
}

// ------------------------------------------------------------------ the module
const hasSet = () => { const t = activePlots(); return !!t && !!stateOf(t).set; };
const hasMembers = () => { const t = activePlots(); return !!t && members(stateOf(t)).length > 0; };
const nChosen = () => { const t = activePlots(); return t ? chosen(stateOf(t)).length : 0; };
const hasSeq = () => { const t = activePlots(); return !!t && !!stateOf(t).seq; };
const viewCmd = ([k, title]) => ({ id: "plots.view." + k, group: "plots", title: `Show ${title}`, what: { overlay: "Every member as measured, dimmed at will.", stack: "Each member raised by the core's displacement, labelled at its baseline.", projection: "WinPLOTR's pseudo-3D: each member moved behind the one before, the nearer hiding the farther (plots::project), an oblique axis for the temperature.", map: "The set as a map: one row per member at its temperature, colour from the core's range.", difference: "Every member less the reference member, on a diverging colour scale symmetric about zero.", surface: "The map as a surface in 3D: each member's row at its temperature in depth, the value up and in colour (plots::map); turned by dragging.", series: "Windows' integrals, centroids and peaks against the temperature (plots::windows); for a .seq, its refined parameters against its external parameter." }[k], when: k === "series" ? () => hasSet() || hasSeq() : hasSet, check: () => { const t = activePlots(); return !!t && stateOf(t).view === k; }, run: () => { const t = activePlots(); stateOf(t).view = k; renderActive(); } });

export const plots = {
  id: "plots",
  name: "Plots",
  icon: `<svg viewBox="0 0 20 20" width="20" height="20" fill="none" stroke="currentColor" stroke-width="1.1" stroke-linecap="round" stroke-linejoin="round"><path d="M2.5 15.5 L6 15.5 L7.5 11 L9 15.5 L17.5 15.5"/><path d="M2.5 11.5 L6.5 11.5 L8 7.5 L9.5 11.5 L17.5 11.5" opacity=".7"/><path d="M2.5 7.5 L7 7.5 L8.5 3.5 L10 7.5 L17.5 7.5" opacity=".45"/></svg>`,
  empty: "Plots shows several data files at once (R32): Plots ▸ Open Buffer… opens a WinPLOTR .buf as a set, New Plot Set… starts one from files of any format. Each member's temperature comes from its file, with its source; the views are an overlay, a stack and a map.",
  /** A WinPLOTR buffer opens in Plots before any other module (the Editor opens every file). */
  // a plot set by its content (PL24: the core's describe knows `format = "drs plot set 1"`), by its place in Plots/, or
  // by its extension
  claims: (file) => !file.dir && (file.kind === "buffer_file" || file.kind === "fullprof_sequential" || file.kind === "plot_set" || file.dict_kind === "plot_set" || /\.drv$/i.test(file.name || "")),
  opens: (file) => !file.dir && (file.kind === "buffer_file" || file.kind === "fullprof_sequential" || file.kind === "plot_set" || file.dict_kind === "plot_set" || /\.drv$/i.test(file.name || "")),
  /** A set a host opened itself, shown in a tab (W7): the set as the core gave it, its title, the file's path and the
      set file's settings when it came from one; the tab is made, opened and made active. */
  adopt: (set, title, path = null, settings = null) => adopt(set, title, path, settings, path && /\.drv$/i.test(path) ? path.split("/").pop() : null),
  /** File ▸ Save (⌘S) saves the plot set shown (PL12). */
  canSave: () => { const t = activePlots(); return !!t && !!stateOf(t).set && members(stateOf(t)).length > 0 && !!rootOf() && writable(); },
  save: () => saveSet(),
  commands: [
    { id: "plots.new", group: "plots", title: "New Plot Set…", what: "A plot set of data files of any format chosen in the dialog, each read by the core's registry, its temperature from its file with the source (R32, PL1–PL2). Kept in memory: the plot-set file waits for the core (PL12).", toolbar: "New Set", run: () => newSet() },
    { id: "plots.add", group: "plots", title: "Add Files…", what: "Data files added to the plot set shown; the core sorts the members again by temperature.", when: hasSet, run: () => addFiles() },
    { id: "plots.openbuf", group: "plots", title: "Open Buffer…", what: "A WinPLOTR buffer (.buf) opened as a plot set: its files relative to it, the text after a name (\"T = 5 K\") read as the temperature, its TITLE line as the set's title (plots::open_buf).", toolbar: "Buffer", run: () => openBuffer() },
    { id: "plots.opencsv", group: "plots", title: "Open CSV as a Set…", what: "A CSV of several columns (the core's own wide or long CSV, or any table) read by the core as a new plot set (PL8): one member per y column, σ paired by name, each member's file and temperature from the # lines above the header, with their source.", run: () => openFilesAsSet() },
    ...VIEWS.map(viewCmd),
    { id: "plots.surface.copy", group: "plots", title: "Copy Surface as Image", what: "The surface as it is turned now onto the clipboard as a PNG (the shared 3D canvas's Copy).", when: () => { const t = activePlots(); const v = t && views.get(t.id); return !!(v && v.surface()); }, run: async () => { const t = activePlots(); const s = t && views.get(t.id) && views.get(t.id).surface(); if (!s) return; try { await s.copyImage(); status("The surface is on the clipboard as an image (PNG)."); } catch (e) { status("The image could not be copied: " + e, true); } } },
    { id: "plots.normalise", group: "plots", title: "Normalise…", what: "The members normalised as a view, by the core (plots::normalise, PL6): by the monitor, the counting time, the maximum (in a range), the same area over a range, a peak's integral above its background, or a reference member's area; every view says it, the files untouched.", sep: true, when: hasMembers, check: () => { const t = activePlots(); return !!t && !!stateOf(t).set && !!stateOf(t).set.normalisation; }, run: () => normaliseView() },
    { id: "plots.convertx", group: "plots", title: "Convert X…", what: "The members shown in 2θ, d, Q, sin θ/λ or s = 1/d, converted by the core with each member's λ (plots::x_view, PL10); the members without one listed, and a control file's λ1 or a typed one offered for them; a view, the files untouched.", when: hasMembers, check: () => { const t = activePlots(); return !!t && !!stateOf(t).set && !!stateOf(t).set.x_quantity; }, run: () => convertX() },
    { id: "plots.lambda", group: "plots", title: "Set Wavelength…", what: "A λ typed for the chosen member, with its source, which Convert X… uses before its file's.", when: () => nChosen() > 0, run: () => setMemberLambda() },
    { id: "plots.export.csv", group: "plots", title: "Plot Set as CSV…", what: "The plot set as CSV, the core's text (plots::csv, csv_long; PL7): wide on the common grid or long, in 2θ, d or Q, what each column is in # lines; written into the project's Exports/ by default.", when: hasMembers, run: () => exportCsv() },
    { id: "plots.window.add", group: "plots", title: "Add Window…", what: "A named 2θ range: for every member the core gives the integral above the background with σ, the centroid with σ and the strongest peak (plots::windows); drawn against the temperature in the series view, as a band in the overlay and the stack.", sep: true, when: hasMembers, run: () => addWindow() },
    { id: "plots.window.draw", group: "plots", title: "Draw Window", what: "The next drag across the overlay or the stack draws a window; its edges are then dragged on the plot, and a click on a band chooses its window. The series view follows each change.", when: hasMembers, check: () => { const t = activePlots(); return !!t && !!stateOf(t).drawWin; }, run: () => drawWindow() },
    { id: "plots.window.rename", group: "plots", title: "Rename Window…", what: "The chosen window's name (a click on its band or its row in the sidebar chooses it).", when: () => { const t = activePlots(); return !!t && stateOf(t).selwin != null && !!stateOf(t).wins[stateOf(t).selwin]; }, run: () => renameWindow() },
    { id: "plots.window.remove", group: "plots", title: "Remove Window", what: "The chosen window removed, the others kept.", when: () => { const t = activePlots(); return !!t && stateOf(t).selwin != null && !!stateOf(t).wins[stateOf(t).selwin]; }, run: () => removeWindow() },
    { id: "plots.window.clear", group: "plots", title: "Remove the Windows", what: "Every window of the set removed.", when: () => { const t = activePlots(); return !!t && stateOf(t).wins.length > 0; }, run: () => removeWindows() },
    { id: "plots.window.at", group: "plots", title: "Background at…", what: "The background the windows are taken above, read by the core at chosen 2θ for every member, against the temperature.", when: hasMembers, run: () => backgroundAt() },
    { id: "plots.results", group: "plots", title: "Open Results Series…", what: "The results files of a series of refinements (Results/*.cif) read by the core (results::read): every refined parameter with its su, χ² and R_wp against each refinement's temperature, read by the core from the refinement's data file or typed with its source, or against the file order.", run: () => openResults() },
    { id: "plots.results.publish", group: "plots", title: "Send Series to Publish…", what: "The refinements of the results series shown handed to Publish as its F5 (the modes' amplitudes against T) or F4 (the cell against T), each results file with its temperature and that value's source; Publish typesets the figure in LaTeX.", when: () => { const t = activePlots(); return !!t && t.rsId != null && publishTakesSeries(); }, disabled: "Open a results series (Open Results Series…); Publish takes a series of refinements with its request PL34.", run: () => sendSeriesToPublish() },
    { id: "plots.seq", group: "plots", title: "Open Sequence (.seq)…", what: "FullProf's sequential results (a .seq) read by the core (seq::read): every refined parameter, χ² and each phase's R_Bragg, volume and weight fraction against the external parameter (PARVAL) or the file order.", run: () => openSeq() },
    { id: "plots.clip", group: "plots", title: "Set Colour Range…", what: "The range the map, the difference map and the surface show on their colour scale, typed (the two marks on the colour bar drag the same pair): values below it drawn in the first colour, above it in the last, saturated, so that the weak details come up as the strong peaks are clipped (I08, R33); per view, on the scale shown. Kept with the tab until the set file keeps it (PL37).", when: () => { const t = activePlots(); return !!t && ["map", "difference", "surface"].includes(stateOf(t).view) && !!stateOf(t).map; }, disabled: "Show the map, the difference map or the surface first.", run: () => setColourRange() },
    { id: "plots.clip.reset", group: "plots", title: "Reset Colour Range", what: "The map's (or the difference map's) colour range back to the core's; a double-click on the colour bar does the same.", when: () => { const t = activePlots(); return !!t && !!stateOf(t).clip[clipKeyOf(stateOf(t))]; }, run: () => { const t = activePlots(); if (t) resetClip(t, clipKeyOf(stateOf(t))); } },
    { id: "plots.cursors", group: "plots", title: "Row and Column Cursors", what: "In the map and the difference map, the row under the pointer (that member's pattern, as the set holds it, drawn below the map) and the column (that 2θ's values against the temperature, beside it): the two sections through the map that show a transition as a kink (I08). A click on the map pins them, a click again releases; the readout gives the cell as before.", when: hasSet, check: () => { const t = activePlots(); return !!t && !!stateOf(t).cursors; }, run: () => toggleCursors() },
    { id: "plots.follow", group: "plots", title: "Range from the View", what: "The map's colour range follows its zoom: the range of the rectangle shown, taken from the core's block of it (plots::map_region), never read off the matrix in the window (I08 step 4, PL37); off puts the colour range back as it was.", when: () => { const t = activePlots(); return !!t && ["map", "difference"].includes(stateOf(t).view) && !!stateOf(t).map; }, check: () => { const t = activePlots(); return !!t && !!stateOf(t).rangeFromView; }, run: () => toggleFollow() },
    { id: "plots.region.save", group: "plots", title: "Save Region…", what: "The rectangle the map shows (its zoom, or the whole) saved by name in the set file (PL37), in the quantity the map is shown in, with its colour range when Range from the View is on; choosing it in the sidebar zooms to it again.", when: () => { const t = activePlots(); return !!t && ["map", "difference"].includes(stateOf(t).view) && !!stateOf(t).map; }, run: () => saveRegion() },
    { id: "plots.region.rename", group: "plots", title: "Rename Region…", what: "The chosen region's name (a click on its row in the sidebar chooses it).", when: () => { const t = activePlots(); return !!t && stateOf(t).selregion != null && !!stateOf(t).regions[stateOf(t).selregion]; }, run: () => renameRegion() },
    { id: "plots.region.remove", group: "plots", title: "Remove Region", what: "The chosen region removed, the others kept; Save writes the set file without it.", when: () => { const t = activePlots(); return !!t && stateOf(t).selregion != null && !!stateOf(t).regions[stateOf(t).selregion]; }, run: () => removeRegion() },
    { id: "plots.yinterval", group: "plots", title: "Y Interval Adapted", what: "The y interval adapted to the member chosen (Overlay and Stack fitted to its curve) and to the row under the pointer (the strip under the map), or shared by the whole set (the curves fitted to every member, the strip on the set's extent, the core's) when off (W4, I10 step 0). Kept in the tab: the set file has no field for it yet (PL41).", when: hasMembers, check: () => { const t = activePlots(); return !!t && stateOf(t).yInterval === "own"; }, run: () => toggleYInterval() },
    { id: "plots.pick.all", group: "plots", title: "Choose All Members", what: "Every member chosen (the checkboxes in the sidebar).", sep: true, when: hasMembers, run: () => pickAll(activePlots(), true) },
    { id: "plots.pick.none", group: "plots", title: "Choose No Member", what: "No member chosen.", when: hasMembers, run: () => pickAll(activePlots(), false) },
    { id: "plots.pick.class", group: "plots", title: "Choose by Exposure Class…", what: "The members whose files state one counting time, chosen together (the exposure class as the files give it; the core's class by the esd's median comes with W2).", when: hasMembers, run: () => pickClass() },
    { id: "plots.pick.range", group: "plots", title: "Choose by a Range of the Variable…", what: "The members whose temperature lies from a value to another, chosen together.", when: hasMembers, run: () => pickRange() },
    { id: "plots.list.order", group: "plots", title: "List Order…", what: "The members' list in the sidebar ordered as the set is, by the variable, by the file's time or by name; the set's own order is the core's (Sort by Variable, Reverse Order).", when: hasMembers, run: () => setListOrder() },
    { id: "plots.viewerpage", group: "plots", title: "Viewer Page… (Dev)", what: "The Map viewer's page written from this set, as drs plots show writes it (PL38, I08 step 5): one HTML file with the set's maps, members, colour range and regions as JSON and the viewer bundled inside, nothing fetched; for looking at the page without the app.", when: hasMembers, run: () => exportViewerPage() },
    { id: "plots.linked", group: "plots", title: "Linked Zoom", what: "The overlay as two plots: the whole range on the left and a window of it on the right, the zoom and the cursor linked, the core's difference of two members under both (F6).", when: hasSet, check: () => { const t = activePlots(); return !!t && stateOf(t).linked; }, run: () => { const t = activePlots(); const p2 = stateOf(t); p2.linked = !p2.linked; p2.view = "overlay"; renderActive(); } },
    { id: "plots.reverse", group: "plots", title: "Reverse Order", what: "The members in the reverse order: the front of the projection, the bottom of the stack and the colours follow.", sep: true, when: hasMembers, run: () => reorder("plots_reverse") },
    { id: "plots.sort", group: "plots", title: "Sort by Variable", what: "The members sorted by temperature (the core's PlotSet::sort_by_variable), those without one last.", when: hasMembers, run: () => reorder("plots_sort") },
    { id: "plots.readvar", group: "plots", title: "Read Variable From…", what: "Each member's temperature from the file's name, or from the file itself as opened, where the core finds it; a temperatures file waits for the core's reader (PL21).", when: hasMembers, run: () => readVariable() },
    { id: "plots.solo", group: "plots", title: "Solo Member", what: "The chosen member at full opacity, every other dimmed.", sep: true, when: () => nChosen() > 0, run: () => solo() },
    { id: "plots.dimothers", group: "plots", title: "Dim Others", what: "The chosen members at full opacity, the others dimmed.", when: () => nChosen() > 0, run: () => dimOthers() },
    { id: "plots.restore", group: "plots", title: "Restore All", what: "Every member at full opacity.", when: hasMembers, run: () => restoreAll() },
    { id: "plots.displacement", group: "plots", title: "Set Displacement…", what: "The stack's rule: a step per member, or in proportion to the temperature, in the y units shown; the core displaces (plots::displace).", sep: true, when: hasMembers, run: () => setDisplacement() },
    { id: "plots.reference", group: "plots", title: "Set Reference", what: "The chosen member is subtracted from every other in the difference map (plots::map with a reference).", when: () => nChosen() > 0, run: () => setReference() },
    { id: "plots.variable", group: "plots", title: "Set Temperature…", what: "The chosen member's temperature typed by you, with a note of where it comes from (L169): a file that states none is drawn apart until then.", when: () => nChosen() > 0, run: () => setVariable() },
    { id: "plots.sum", group: "plots", title: "Sum of Members", what: "The chosen members summed with σ on the common grid, written as a new file Plots/{set} sum.xys that names its members (plots::derive).", sep: true, when: () => nChosen() > 0, run: () => derive("sum") },
    { id: "plots.average", group: "plots", title: "Average of Members", what: "The chosen members averaged with σ on the common grid, written as Plots/{set} average.xys (plots::derive).", when: () => nChosen() > 0, run: () => derive("average") },
    { id: "plots.difference", group: "plots", title: "Difference of Members", what: "The first member chosen less the second, with σ, written as Plots/{set} difference.xys (plots::derive).", when: () => nChosen() === 2, disabled: "Choose two members: the first chosen less the second.", run: () => derive("difference") },
    { id: "plots.opendata", group: "plots", title: "Open in Data", what: "The chosen member's file in Data: its operations, peaks and background.", sep: true, when: () => nChosen() > 0, run: () => openInData() },
    { id: "plots.save", group: "plots", title: "Save Plot Set", what: "The plot set as Plots/{name}.drv in the project, through the core (PL12): the members, each typed temperature with its source, the order, the stack's rule, the y scale, the reference, the windows, each member's opacity, the map's colour range and the regions (PL37). File ▸ Save (⌘S) does the same.", when: () => plots.canSave(), disabled: "Open a plot set with members in a project of your own.", run: () => saveSet() },
    { id: "plots.saveas", group: "plots", title: "Save Plot Set As…", what: "The plot set under another name in Plots/.", when: () => plots.canSave(), run: () => saveSet(activePlots(), true) },
    { id: "plots.openset", group: "plots", title: "Open Plot Set…", what: "A plot set file (.drv) read by the core: its members after Data's operations, the typed temperatures with their sources, the order, the rule, the reference, the windows. A double-click on a .drv in the project does the same.", run: () => openSetFile() },
    { id: "plots.publish", group: "plots", title: "Send to Publish…", what: "The plot set saved as Plots/{name}.drv, with the projection's steps and the scale they were set on, and handed to Publish, which makes its publication set of the view shown (F3's stack, projection, map or difference map, or the series view) and opens it; Build there typesets the figure with LuaLaTeX in Latin Modern.", when: () => plots.canSave(), disabled: "Open a plot set with members in a project of your own.", run: () => sendToPublish() },
    { id: "plots.figurelook", group: "plots", title: "Figure Look", what: "Every view as a publication figure would show it: on white, open axes, text in ink; the stack and the projection with a scale bar in place of the y axis, the projection's x axis over the front member's range and its temperature axis with ticks and title; a map's colour bar saying its scale; Export ▸ Plot writes it at the figure's width.", when: hasSet, check: () => { const t = activePlots(); return !!t && !!stateOf(t).look; }, run: () => { const t = activePlots(); const p2 = stateOf(t); p2.look = !p2.look; renderActive(); } },
    { id: "plots.scalebar", group: "plots", title: "Set Scale Bar…", what: "The length of the figure's scale bar in the stack and the projection, in the y units shown; empty for the core's round step (plots::proposed_step).", when: hasSet, run: () => setScaleBar() },
    { id: "plots.figsize", group: "plots", title: "Figure Size…", what: "The width a figure is exported at in its figure look, and its text's sizes: a journal's or a thesis's figure width and sizes, from Publish's presets.", when: hasSet, run: () => setFigureSize() },
    { id: "plots.reload", group: "plots", title: "Reload Members", what: "The members read again by the core, after the operations Data keeps for their files now (PL20).", when: hasMembers, run: () => reloadMembers() },
  ],
  menu: ["plots.new", "plots.openset", "plots.openbuf", "plots.opencsv", "plots.add", "plots.seq", "plots.results", "plots.results.publish", "-", "plots.save", "plots.saveas", "plots.reload", "plots.export.csv", "plots.publish", "-", "plots.view.overlay", "plots.view.stack", "plots.view.projection", "plots.view.map", "plots.view.difference", "plots.view.surface", "plots.surface.copy", "plots.view.series", "plots.linked", "plots.clip", "plots.clip.reset", "plots.cursors", "plots.follow", "plots.region.save", "plots.region.rename", "plots.region.remove", "plots.yinterval", "plots.viewerpage", "-", "plots.pick.all", "plots.pick.none", "plots.pick.class", "plots.pick.range", "plots.list.order", "-", "plots.figurelook", "plots.scalebar", "plots.figsize", "-", "plots.normalise", "plots.convertx", "plots.lambda", "-", "plots.reverse", "plots.sort", "plots.readvar", "-", "plots.window.add", "plots.window.draw", "plots.window.rename", "plots.window.remove", "plots.window.clear", "plots.window.at", "-", "plots.solo", "plots.dimothers", "plots.restore", "-", "plots.displacement", "plots.reference", "plots.variable", "-", "plots.sum", "plots.average", "plots.difference", "-", "plots.opendata"],
  manual: { title: "Plots", html: "<p>Several data files at once as a plot set: an overlay with dimming, a stack, a projection with hidden parts, a map and a difference map, and windows' quantities against the temperature, normalised and in 2θ, d, Q, sin θ/λ or 1/d as views; exported as CSV; FullProf's .seq against its external parameter. The manual's Plots page says more.</p>" },
  /** The figure File ▸ Export ▸ Plot writes (L134): the view shown in the active Plots tab. */
  figure: () => { const t = activeTab(); const v = t && views.get(t.id); const f = v && v.figure ? v.figure() : null; return f ? { ...f, after: () => v.redraw() } : null; },
  /** The view on screen, for Help ▸ Manual for What Is on Screen. */
  helpView: (tab) => (tab && tab.plots ? (tab.plots.seq ? "sequence" : tab.plots.rs ? "series" : tab.plots.view) : null),

  async open(file) {
    if (!hasCore()) return null;
    if (file.kind === "fullprof_sequential") return openSeqPath(file.absolute, file.name, file.path);
    if (file.kind === "plot_set" || file.dict_kind === "plot_set" || /\.drv$/i.test(file.name || "")) {
      const [set, st] = await invoke("plots_open_set", { path: file.absolute, root: rootOf() });
      const tab = tabFor(set, file.name.replace(/\.drv$/i, ""), file.path, st, file.path);
      tab.saveName = file.name.replace(/\.drv$/i, "");
      return tab;
    }
    const set = await invoke("plots_open_buffer", { path: file.absolute, root: rootOf() });
    return tabFor(set, file.name, file.path);
  },
  view(tab, box) { return viewOf(tab, box); },
  sidebar(box) {
    const t = activePlots();
    const tabs = [...state.tabs.values()].filter((x) => x.module === "plots");
    if (!tabs.length) { box.appendChild(el("div", "empty", esc(plots.empty))); return; }
    const list = el("div", "sourcelist");
    list.appendChild(el("div", "sec-head", "<span>Plot sets</span><span class=\"sub\">in memory, not saved</span>"));
    for (const x of tabs) {
      const n = x.plots && x.plots.set ? x.plots.set.members.length : 0;
      const what = x.plots && x.plots.seq ? `FullProf sequence, ${x.plots.seq.blocks.length} refinements` : x.plots && x.plots.rs ? `results series, ${x.plots.rs.points.length} refinements` : `${n} member${n === 1 ? "" : "s"}`;
      const row = el("div", "child" + (t && x.id === t.id ? " sel" : ""), `<span>${esc(x.title)}<br><span class="sub">${esc(what)}</span></span>`);
      row.style.padding = "4px 12px"; row.style.cursor = "default"; row.onclick = () => activateTab(x.id); list.appendChild(row);
    }
    box.appendChild(list);
    if (!t) return;
    const ps = stateOf(t);
    if (ps.rs) {
      box.appendChild(el("div", "sec-head", `<span>Refinements</span><span class="sub">T (K)</span>`));
      const rl = el("div", "sourcelist");
      ps.rs.points.forEach((p, k) => {
        const row = el("div", "child pl-member", `<span class="pl-lab">${esc(p.label)}</span><span class="pl-var sub">${esc(p.variable ? `${fmt(p.variable.value)}${p.variable.typed ? " typed" : ""}` : "none")}</span>`);
        row.title = `${p.data_file || "no data file"}\n${p.variable ? `T = ${fmt(p.variable.value)} K, from ${p.variable.source}` : `no temperature: ${p.no_data || "its data file says none the core reads"} (double-click to type one)`}\nχ² ${fmt(p.chi2)}, R_wp ${fmt(p.rwp)}`;
        row.ondblclick = () => setResultsVariable(t, k); rl.appendChild(row);
      });
      box.appendChild(rl);
      box.appendChild(el("div", "sec-head", `<span>Quantities</span><span class="sub">${ps.rs.names.length} parameters</span>`));
      const ul = el("div", "sourcelist");
      for (const [k, label] of resultsQuantities(ps.rs)) {
        const row = el("div", "child" + (k === ps.param ? " sel" : ""), `<span>${esc(label)}</span>`); row.style.padding = "2px 12px"; row.style.cursor = "default";
        row.onclick = () => { ps.param = k; renderActive(); }; ul.appendChild(row);
      }
      box.appendChild(ul); return;
    }
    if (ps.seq) {
      box.appendChild(el("div", "sec-head", `<span>Quantities</span><span class="sub">${ps.seq.names.length} parameters</span>`));
      const ul = el("div", "sourcelist");
      for (const [k, label] of seqQuantities(ps.seq)) {
        const row = el("div", "child" + (k === ps.param ? " sel" : ""), `<span>${esc(label)}</span>`); row.style.padding = "2px 12px"; row.style.cursor = "default";
        row.onclick = () => { ps.param = k; ps.seqcol = null; renderActive(); }; ul.appendChild(row);
      }
      box.appendChild(ul); return;
    }
    const set = ps.set; if (!set) return;
    if (ps.wins.length) {
      box.appendChild(el("div", "sec-head", `<span>Windows</span><span class="sub">2θ (°)</span>`));
      const wl = el("div", "sourcelist");
      ps.wins.forEach((w, k) => {
        const row = el("div", "child" + (k === ps.selwin ? " sel" : ""), `<span>${esc(w.name)}</span><span class="sub">${fmt(w.from)}–${fmt(w.to)} <button class="small pl-winx" title="Remove Window ${esc(w.name)}" aria-label="Remove Window ${esc(w.name)}">×</button></span>`);
        row.style.padding = "2px 12px"; row.style.display = "flex"; row.style.justifyContent = "space-between"; row.style.cursor = "default";
        row.title = `${w.name}: ${fmt(w.from)} to ${fmt(w.to)}. Click to choose it (its band is marked), double-click to rename it; its edges are dragged on the plot.`;
        row.onclick = (e) => { if (e.target.classList.contains("pl-winx")) { removeWindow(t, k); return; } ps.selwin = k; renderActive(); };
        row.ondblclick = () => renameWindow(t, k);
        wl.appendChild(row);
      });
      box.appendChild(wl);
    }
    // step 17 (PL37): the regions saved by name, beside the windows; a click zooms to one, a double-click renames it
    if (ps.regions.length) {
      box.appendChild(el("div", "sec-head", `<span>Regions</span><span class="sub">${esc(xNameOf(ps))} × ${esc(set.variable_name)}</span>`));
      const rl = el("div", "sourcelist");
      ps.regions.forEach((r, k) => {
        const row = el("div", "child" + (k === ps.selregion ? " sel" : ""), `<span>${esc(r.name)}</span><span class="sub">${fmt(r.x[0])}–${fmt(r.x[1])}${r.x_quantity && r.x_quantity !== xWordOf(ps) ? ` ${esc(r.x_quantity)}` : ""}, ${fmt(r.variable[0])}–${fmt(r.variable[1])} <button class="small pl-regx" title="Remove Region ${esc(r.name)}" aria-label="Remove Region ${esc(r.name)}">×</button></span>`);
        row.style.padding = "2px 12px"; row.style.display = "flex"; row.style.justifyContent = "space-between"; row.style.cursor = "default";
        // P65: a region's own colour range per view, the map's and the difference map's apart
        const ranges = [r.map_range ? `on the map ${fmt(r.map_range.from)} to ${fmt(r.map_range.to)}` : "", r.difference_range ? `on the difference map ${fmt(r.difference_range.from)} to ${fmt(r.difference_range.to)}` : ""].filter(Boolean);
        row.title = `${r.name}: ${r.x_quantity || xWordOf(ps) || "x"} ${fmt(r.x[0])} to ${fmt(r.x[1])}, ${set.variable_name} ${fmt(r.variable[0])} to ${fmt(r.variable[1])} ${set.variable_unit}${ranges.length ? `; its colour range ${ranges.join(", ")}` : ""}${r.window ? `; drawn from the window ${r.window}` : ""}. Click to zoom the map to it (its colour range on the view shown applied), double-click to rename it.`;
        row.onclick = (e) => { if (e.target.classList.contains("pl-regx")) { removeRegion(t, k); return; } chooseRegion(t, k); };
        row.ondblclick = () => renameRegion(t, k);
        rl.appendChild(row);
      });
      box.appendChild(rl);
    }
    box.appendChild(el("div", "sec-head", `<span>Members</span><span class="sub">${esc(set.variable_name)} (${esc(set.variable_unit)}) · opacity</span>`));
    // step 20 (W4): the list's order, and the members chosen all, none, by exposure class or by a range of the variable
    const ctl = el("div", "pl-memctl sub", `<label>order <select class="pl-sort" aria-label="The members' order in the list"><option value="set">the set's</option><option value="variable">by ${esc(set.variable_name)}</option><option value="time">by the file's time</option><option value="name">by name</option></select></label> <button class="small pl-pick-all" title="Every member chosen">all</button><button class="small pl-pick-none" title="No member chosen">none</button><button class="small pl-pick-class" title="The members of one exposure class (the counting time the files state)">class…</button><button class="small pl-pick-range" title="The members whose ${esc(set.variable_name)} lies in a range">range…</button>`);
    ctl.querySelector(".pl-sort").value = ps.sortBy || "set";
    ctl.querySelector(".pl-sort").onchange = (e) => { ps.sortBy = e.target.value; if (window.__renderPanels) window.__renderPanels(); };
    ctl.querySelector(".pl-pick-all").onclick = () => pickAll(t, true);
    ctl.querySelector(".pl-pick-none").onclick = () => pickAll(t, false);
    ctl.querySelector(".pl-pick-class").onclick = () => pickClass(t);
    ctl.querySelector(".pl-pick-range").onclick = () => pickRange(t);
    box.appendChild(ctl);
    const ul = el("div", "sourcelist pl-members");
    for (const k of sortMembers(set.members.map((m) => ({ label: m.label, variable: m.variable, time: m.modified })), ps.sortBy)) {
      const m = set.members[k];
      const on = ps.sel.includes(m.key);
      const row = el("div", "child pl-member" + (on ? " sel" : ""));
      row.title = memberFacts(ps, m);
      // the icon by the file's kind (the core's format, kindIcon's family), a checkbox per member (chosen), the label,
      // the variable, the opacity
      row.innerHTML = `<input type="checkbox" class="pl-pick" ${on ? "checked" : ""} aria-label="${esc(m.label)} chosen">${kindIcon(m.format, m.format_name)}<span class="pl-sw" style="background:${m.variable ? colourOf(ps, m) : "var(--muted)"}"></span><span class="pl-lab">${esc(m.label)}${m.operations.length ? ` <span class="pl-ops" title="drawn after Data's operations">ops ${m.operations.length}</span>` : ""}</span><span class="pl-var sub">${esc(m.shown.left_out ? "out" : m.variable ? `${fmt(m.variable.value)}${m.variable.typed ? " typed" : ""}` : "apart")}</span><input type="range" class="pl-op" min="0" max="1" step="0.05" value="${opacityOf(ps, m)}" aria-label="Opacity of ${esc(m.label)}">`;
      row.onclick = (e) => { if (e.target.classList.contains("pl-op") || e.target.classList.contains("pl-pick")) return; if (e.altKey) { soloToggle(t, m.key); return; } select(t, m.key, e.metaKey || e.shiftKey); };
      row.ondblclick = () => { if (!m.variable) { select(t, m.key, false); setVariable(t); } };
      row.querySelector(".pl-pick").onchange = (e) => { pickOneMember(t, m.key, e.target.checked); };
      const op = row.querySelector(".pl-op"); op.oninput = () => { ps.opacity[m.key] = parseFloat(op.value); renderTab(t); };
      ul.appendChild(row);
    }
    box.appendChild(ul);
    const mon = monitorNote(ps);
    if (mon) box.appendChild(el("div", "empty small", esc(mon)));
    const vt = viewText(ps);
    if (vt) box.appendChild(el("div", "empty small", esc(vt)));
    const withOps = set.members.filter((m) => m.operations.length);
    box.appendChild(el("div", "empty small", esc(withOps.length ? `${withOps.length} member${withOps.length === 1 ? "" : "s"} drawn after the operations Data keeps for ${withOps.length === 1 ? "its file" : "their files"} (marked "ops"); the others as measured. Reload Members reads them again.` : "Every member as measured: Data keeps no operation for these files.")));
    if (set.apart.length) box.appendChild(el("div", "empty small", esc(`Drawn apart (no ${set.variable_name} in the file): ${set.apart.map((k) => set.members[k].label).join(", ")}. Double-click one, or Set ${set.variable_name === "T" ? "Temperature" : "Value"}…, to type one with its source.`)));
  },
  inspector(box) {
    const inner = el("div", "inspector"); box.appendChild(inner);
    const t = activePlots(); const ps = t ? stateOf(t) : null;
    const c = ps ? chosen(ps) : [];
    if (ps && ps.rs) {
      const v = ps.rs;
      inner.appendChild(el("div", "insp-head", `<div class="kind">results series</div><div class="title">${esc(t.title)}</div><div class="sub">${v.points.length} refinements · ${v.names.length} parameters</div>`));
      inner.appendChild(el("div", "insp-sec", `<h3>Read by the core</h3><div>results::read: each refinement's refined parameters with their su, χ² and the R factors; its temperature from the data file its control file names (read by plots::open: the header, the title, the name), or typed with its source. A refinement whose data file says none is not drawn until one is typed (L169).</div>`));
      return;
    }
    if (ps && ps.seq) {
      const v = ps.seq;
      inner.appendChild(el("div", "insp-head", `<div class="kind">FullProf sequence</div><div class="title">${esc(t.title)}</div><div class="sub">${v.blocks.length} refinements · ${v.names.length} parameters</div>`));
      inner.appendChild(el("div", "insp-sec", `<h3>Read by the core</h3><div>seq::read: the blocks of NEW_REFINEMENT, each with PARVAL's index, external parameter and numor, χ², each phase's R_Bragg, volume and weight fraction, and every refined parameter with σ.${v.repeated.length ? ` Refined more than once: ${esc(v.repeated.join(", "))}.` : ""}${v.unread ? ` ${v.unread} line(s) not read, kept.` : ""}</div>`));
      return;
    }
    if (!ps || !ps.set) { inner.appendChild(el("div", "empty", "Open a plot set: Plots ▸ Open Buffer… or New Plot Set…; or a FullProf sequence: Open Sequence (.seq)….")); return; }
    if (!c.length) {
      const set = ps.set;
      inner.appendChild(el("div", "insp-head", `<div class="kind">plot set</div><div class="title">${esc(t.title)}</div><div class="sub">${set.members.length} members · ${set.same_grid ? "one grid shared" : "different grids: the map and sums use the common grid"}</div>`));
      inner.appendChild(el("div", "insp-sec", `<h3>File</h3><div>${t.drv ? `<code>${esc(t.drv)}</code>: Save (⌘S) writes it again.` : "Not saved yet: Save (⌘S) writes the set as <code>Plots/{name}.drv</code> in the project (the core's set file, PL12)."}</div>`));
      const mon = monitorNote(ps); if (mon) inner.appendChild(el("div", "insp-sec", `<h3>Counts</h3><div>${esc(mon)}. Normalise… shows them by a monitor, a counting time, a maximum, the same area, a peak's integral or a reference member, as a view (the core's, PL6).</div>`));
      const vt = viewText(ps); if (vt) inner.appendChild(el("div", "insp-sec", `<h3>The view</h3><div>${esc(vt)}</div>`));
      inner.appendChild(el("div", "empty", "Click a member in the sidebar (⌘-click for several) to see what its file says."));
      return;
    }
    const m = c[0], set = ps.set;
    inner.appendChild(el("div", "insp-head", `<div class="kind">member${c.length > 1 ? ` (1 of ${c.length} chosen)` : ""}</div><div class="title">${esc(m.label)}</div><div class="sub">${esc(m.file)}${m.scans > 1 ? ` · scan ${m.scan + 1} of ${m.scans}` : ""}</div>`));
    inner.appendChild(el("div", "insp-sec", `<h3>${esc(set.variable_name)}</h3><div>${esc(variableText(ps, m, true))}</div>${m.candidates.length > 1 ? `<div class="sub">also: ${m.candidates.slice(1).map((v) => esc(`${fmt(v.value)} ${set.variable_unit} from ${v.source}`)).join("; ")}</div>` : ""}`));
    const facts = [["format", m.format_name], ["axis", m.axis], ["points", String(m.points)], ["range", `${fmt(m.range[0])} to ${fmt(m.range[1])}`], ["step", m.step == null ? "" : fmt(m.step)], ["wavelength", m.wavelength == null ? "not stated" : `${fmt(m.wavelength)} Å`], ["monitor", m.monitor == null ? "none: raw counts" : fmt(m.monitor)], ["counting time", m.counting_time == null ? "not stated" : `${fmt(m.counting_time)} s`], ["title", m.title || ""], ["buffer comment", m.comment || ""]];
    inner.appendChild(el("div", "insp-sec", `<h3>Data's operations</h3>${m.operations.length ? `<div>Drawn after the operations Data keeps for this file, in order:</div><div class="insp-children">${m.operations.map((o) => `<div class="child"><span>${esc(o)}</span></div>`).join("")}</div>` : "<div>None: the file as measured.</div>"}`));
    if (set.normalisation || set.x_quantity) inner.appendChild(el("div", "insp-sec", `<h3>In the view</h3><div>${esc(m.shown.left_out ? `Left out: ${m.shown.left_out}.` : [m.shown.factor != null ? `normalised: × ${fmt(m.shown.factor)} (the core's factor)` : "", m.shown.lambda != null ? `in ${set.x_label} with λ ${fmt(m.shown.lambda)} Å from ${m.shown.lambda_from}` : set.x_quantity ? `in ${set.x_label}` : ""].filter(Boolean).join("; ") + ".")}</div>`));
    inner.appendChild(el("div", "insp-sec", `<h3>What the file says</h3><div class="insp-children">${facts.filter((f) => f[1]).map(([k, v]) => `<div class="child"><span>${esc(k)}</span><span>${esc(v)}</span></div>`).join("")}</div>`));
    const acts = el("div", "insp-sec"); acts.innerHTML = `<h3>Actions</h3>`;
    // each action with its tooltip, the command's own words (APP_CONVENTIONS §10: no control without a name)
    const actTip = (id) => { const c = plots.commands.find((x) => x.id === id); return c ? c.what : ""; };
    for (const [label, fn, id] of [[`Set ${set.variable_name === "T" ? "Temperature" : "Value"}…`, () => setVariable(t), "plots.variable"], ["Open in Data", () => openInData(t), "plots.opendata"], ["Solo Member", () => solo(t), "plots.solo"], ["Set Reference", () => setReference(t), "plots.reference"], ["Set Wavelength…", () => setMemberLambda(t), "plots.lambda"]]) { const b = el("button", "small", esc(label)); b.style.margin = "2px 4px 2px 0"; b.title = actTip(id); b.onclick = fn; acts.appendChild(b); }
    inner.appendChild(acts);
  },
  bottom(box) { box.innerHTML = `<section class="pl-bottom"><div id="pl-table"></div></section>`; table = null; renderBottom(); },
  close(tab) { views.delete(tab.id); if (hasCore() && tab.setId != null) invoke("plots_close", { id: tab.setId }).catch(() => {}); if (hasCore() && tab.seqId != null) invoke("plots_seq_close", { id: tab.seqId }).catch(() => {}); if (hasCore() && tab.rsId != null) invoke("plots_results_close", { id: tab.rsId }).catch(() => {}); },

  /** The self-check (DRS_APP_SELFTEST, L144) when a buffer is the startup file: the set, each view, the core's numbers. */
  async selftest(note, sleep) {
    // the startup buffer's tab, found by its file: on a relaunch the app first reopens the tabs of the last quit (P46),
    // and the active one may be a set file an earlier run saved (L366, L370, L375: the check then ran on that tab)
    const t0 = activePlots();
    const t = [...state.tabs.values()].find((x) => x.module === "plots" && x.setId != null && /\.buf$/i.test(x.path || "")) || t0;
    if (!t) { await note("plots: no Plots tab"); return; }
    const ps = stateOf(t); const set = ps.set;
    const wait = async (f, ms = 20000) => { const t0 = Date.now(); while (!f()) { if (Date.now() - t0 > ms) return false; await sleep(100); } return true; };
    // step 22 (W9): a large set opened at start (the ESRF run, 101 scans of 27 000 points): the times measured and
    // nothing else checked. The times are the tab's (`ps.timing`: the set's arrival, the overlay's first frame, its
    // whole drawing), written by the module as it draws.
    if (stateOf(t).set && stateOf(t).set.members.length > 50) {
      const sp = stateOf(t); sp.view = "overlay"; sp.linked = false; renderActive();
      await wait(() => !!(sp.timing && sp.timing.whole != null), 1200000);
      const tm = sp.timing || {};
      const pts = sp.set.members[0] ? sp.set.members[0].points : 0;
      await note(`plots: large set ${t.title}: ${sp.set.members.length} members of ${pts} points (the first's); the set arrived ${tm.opened != null ? `${Math.round(tm.opened)} ms` : "?"} after the launch's start; the overlay's first frame ${tm.firstFrame != null ? `${Math.round(tm.firstFrame)} ms` : "?"} after the set, whole ${tm.whole != null ? `${Math.round(tm.whole)} ms` : "not finished"} after the set${tm.patterns != null ? `; the patterns asked of the core in ${Math.round(tm.patterns)} ms` : ""}: ${tm.whole != null ? "OK" : "FAIL"}`);
      await note("done"); return;
    }
    if (t !== t0) { activateTab(t.id); await sleep(500); await note(`plots: the startup buffer's tab ${t.title} made active (the launch had made ${t0 ? t0.title : "none"} active)`); }
    await note(`plots: opened ${t.title}: ${set.members.length} members, ${set.apart.length} apart (${set.apart.map((k) => set.members[k].label).join(", ")}), first ${variableText(ps, set.members[0], true)}; ${monitorNote(ps) || "monitors stated"}; same grid ${set.same_grid}`);
    for (const [k] of VIEWS) {
      if (k === "stack") { ps.rule = { kind: "step", dx: 0, dy: null, proposed: null }; ps.disp = null; }   // the core proposes the step (PL19)
      if (k === "projection") { ps.prule = { kind: "step", dx: 0.22, dy: 2200, scale: "linear" }; ps.proj = null; }
      if (k === "series") { ps.wins = [{ name: "m1", from: 18.5, to: 19.8 }, { name: "m2", from: 26.0, to: 27.4 }, { name: "m3", from: 31.9, to: 33.3 }, { name: "m4", from: 41.5, to: 42.9 }]; ps.winres = null; ps.quantity = "integral"; ps.scale = "log"; }
      if (k === "overlay") { ps.linked = true; ps.zoomwin = [17.6, 21.0]; }
      if (k === "surface") ps.surfOf = "map";
      ps.view = k; renderActive();
      const surfaceOf = () => { const vv = views.get(t.id); const sf = vv && vv.surface(); return sf && sf.stats() ? sf : null; };
      const ok = await wait(() => (k === "overlay" ? !!(ps.diff && ps.diff.v) : k === "stack" ? !!(ps.disp && ps.disp.v) : k === "projection" ? !!(ps.proj && ps.proj.v) : k === "series" ? !!(ps.winres && ps.winres.v) : k === "surface" ? !!surfaceOf() : !!(ps.map && ps.map.key.startsWith(k === "difference" ? "d:" : "m:"))));
      await sleep(300);
      const v = views.get(t.id); const series = k === "surface" ? (surfaceOf() ? [surfaceOf()] : []) : v ? v.plot.state.series : [];
      const cv = t && document.querySelector(k === "surface" ? ".pl-surf canvas" : ".pl-area canvas"); const r = cv ? cv.getBoundingClientRect() : { width: 0, height: 0 };
      let what = "";
      if (k === "overlay") { const z = v && v.zoom(); if (v && z) v.plot.setCursor(19.1); const zv = z ? z.getView() : null; const d = ps.diff && ps.diff.v; what = `${series.length} series on the left (the window's band and the difference among them); the right plot's window ${zv ? `${fmt(zv.x0)} to ${fmt(zv.x1)}°` : "none"}; the core's difference ${d ? `${set.members[d.a].label} − ${set.members[d.b].label}, ${d.x.length} points` : "none"}; a cursor set at 19.1° on the left: ${v && v.plot.state.cursor === 19.1 ? "drawn" : "not drawn"}`; }
      if (k === "stack") what = ps.disp && ps.disp.v ? `${series.length} curves, the core's proposed step ${fmt(ps.rule.dy)} (${ps.rule.proposed ? ps.rule.proposed.rule : "not the core's"}), baselines ${ps.disp.v.filter((d) => d.placed).slice(0, 3).map((d) => d.baseline).join(", ")}…` : "no displacement";
      if (k === "projection") { const P = ps.proj && ps.proj.v; const pr = P && P.projection; const ap = P && P.apart.length ? pr.curves.find((c) => c.member === P.apart[0]) : null; const apX = ap ? ap.pieces.flat().map((q) => q[0]) : []; what = pr ? `${pr.curves.length} curves, ${pr.curves.filter((c) => c.visible < 0.999).length} partly hidden, ${pr.curves.flatMap((c) => c.pieces).length} visible pieces, ${pr.axis.length} axis points; apart: ${P.apart.map((m) => set.members[m].label).join(", ") || "none"}${ap ? `, its 2θ from ${fmt(apX[0])} to ${fmt(apX[apX.length - 1])}° as measured` : ""}` : "no projection"; }
      if (k === "series") { const W = ps.winres && ps.winres.v; const first = W && W.members.find((mw) => set.members[mw.member].variable); const last = W && [...W.members].reverse().find((mw) => set.members[mw.member].variable); what = W ? `${W.members.length} members × ${ps.wins.length} windows; m1 at ${fmt(set.members[first.member].variable.value)} K ${fmt(first.windows[0].integral)} ± ${fmt(first.windows[0].sigma)}, at ${fmt(set.members[last.member].variable.value)} K ${fmt(last.windows[0].integral)} ± ${fmt(last.windows[0].sigma)}` : "no windows"; }
      if (k === "map" || k === "difference") { const M = ps.map && ps.map.v.map; what = M ? `${M.rows.length} rows × ${M.x.length} columns, range ${M.range[0].toFixed(3)} to ${M.range[1].toFixed(3)}, ${ps.map.v.apart.length} row(s) apart` : "no map"; }
      if (k === "surface") {
        // the core's map as a surface (DESIGN.md §8.1 item 6): its rows and columns all drawn, the readout of a cell (the
        // 22 K row at 19.1°, the strongest magnetic peak, as the map's cut in F3) in the core's own numbers
        const sf = surfaceOf(), S = sf && sf.stats(), M = ps.map && ps.map.v.map;
        const ri = M ? M.rows.findIndex((row) => row.value === 22) : -1, cj = M ? M.x.findIndex((x) => x >= 19.1) : -1;
        if (sf && ri >= 0 && cj >= 0) sf.showCell({ i: ri, j: cj });
        const ro = (document.querySelector(".pl-sread") || {}).textContent || "";
        const zc = M && ri >= 0 && cj >= 0 ? M.z[ri][cj] : null;
        const same = !!(S && M && S.rows === M.rows.length && S.columns === M.x.length && zc != null && ro.includes("T 22 K") && ro.includes(zc.toFixed(3)));
        what = S ? `${S.rows} rows × ${S.columns} columns as the core's map, ${S.vertices} vertices, ${S.triangles} triangles, ticks ${S.ticks.join("/")}; ${ps.map.v.apart.length} member(s) left out (no T); the readout at 22 K, ${M && cj >= 0 ? M.x[cj].toFixed(2) : "?"}°: "${ro}" ${same ? "(the core's value)" : "(not the core's: FAIL)"}` : "no surface";
        if (!same) series.length = 0;
        if (sf) sf.showCell(null);
      }
      await note(`plots: ${k}: ${what}; canvas ${Math.round(r.width)}×${Math.round(r.height)} px: ${ok && r.height > 40 && series.length ? "OK" : "FAIL"}`);
      // the view as a figure, by Export ▸ Plot's own path, into the project's Exports/ (the self-check's copy), to be looked at
      if (state.project && state.project.root && writable()) {
        const fig = inLightPalette(() => plots.figure());   // as File ▸ Export ▸ Plot draws it
        if (fig && fig.after) fig.after();
        if (fig) { try { const out = await invoke("export_figure", { svg: fig.svg, path: `${state.project.root}/Exports/selfcheck plots ${k}.png`, dpi: 150 }); await note(`plots: ${k} exported: ${out}`); } catch (e) { await note(`plots: ${k} not exported: ${e}`); } }
      }
    }
    // step 15: the surface let go once another view is shown (L150: no GPU context for a view not shown)
    { const v = views.get(t.id); const gone = !!v && !v.surface() && !document.querySelector(".pl-surf canvas"); await note(`plots: the surface after the ${ps.view} is shown: ${gone ? "let go, its WebGL canvas removed" : "still held"}: ${gone ? "OK" : "FAIL"}`); }
    // step 13: windows edited on the plot (DESIGN.md §4), through the plot's own band callbacks: m1's edges dragged, the
    // series following; a window drawn, chosen, renamed and removed; m1 put back
    {
      const v = views.get(t.id); const was = { ...ps.wins[0] };
      ps.view = "stack"; ps.scale = "linear"; ps.look = false; renderActive(); await sleep(400);
      const bands = (v ? v.plot.state.series : []).find((x) => x.id === "windows");
      const before = ps.winres && ps.winres.v ? ps.winres.v.members[0].windows[0].integral : null;
      v.band.change(0, 18.3, 20.0);
      ps.view = "series"; renderActive(); await wait(() => !!(ps.winres && ps.winres.v)); await sleep(300);
      const after = ps.winres && ps.winres.v ? ps.winres.v.members[0].windows[0].integral : null;
      const n0 = ps.wins.length;
      ps.view = "overlay"; ps.linked = false; drawWindow(t); await sleep(300); const drawing = ps.drawWin && v.plot.state.mode === "region";
      v.band.create(36.5, 38.0); await sleep(300);
      const drawn = ps.wins[ps.wins.length - 1];
      const renamed = await renameWindow(t, ps.wins.length - 1, "drawn");
      v.band.select(0); const chose = ps.selwin === 0;
      removeWindow(t, ps.wins.length - 1); await sleep(200);
      ps.wins = ps.wins.map((w, k) => (k === 0 ? was : w)); ps.winres = null; ps.selwin = null; ps.linked = true; renderActive(); await sleep(300);
      const ok = !!(bands && bands.editable && before != null && after != null && after !== before && drawing && drawn && drawn.from === 36.5 && drawn.to === 38 && renamed && chose && ps.wins.length === n0 && ps.wins[0].from === was.from);
      await note(`plots: windows on the plot: the stack's bands ${bands && bands.editable ? "editable" : "not editable"}; m1 dragged to 18.3–20°, its integral at the first member ${before != null ? fmt(before) : "?"} → ${after != null ? fmt(after) : "?"} (the series followed); Draw Window ${drawing ? "on (region mode)" : "not on"}, a window drawn ${drawn ? `${drawn.name} ${fmt(drawn.from)}–${fmt(drawn.to)}` : "none"}, renamed "drawn" ${renamed ? "yes" : "no"}, a band chosen ${chose ? "yes" : "no"}, removed, ${ps.wins.length} windows as before: ${ok ? "OK" : "FAIL"}`);
    }
    // step 14: dimming by the pointer (DESIGN.md §4): a curve hovered drawn thicker, ⌥-click solos it, again restores
    {
      const v = views.get(t.id); ps.view = "overlay"; ps.linked = false; ps.look = false; ps.opacity = {}; ps.soloKey = null; renderActive(); await sleep(400);
      const key = ps.set.members[3].key;
      v.hover(key); await sleep(300);
      const w = ((v.plot.state.series || []).find((x) => x.id === key) || {}).width;
      soloToggle(t, key); await sleep(300);
      const dimmed = ps.set.members.filter((m) => m.key !== key && opacityOf(ps, m) < 1).length, alone = opacityOf(ps, ps.set.members[3]) === 1;
      soloToggle(t, key); await sleep(300);
      const back = ps.set.members.every((m) => opacityOf(ps, m) === 1);
      v.hover(null); ps.linked = true; renderActive(); await sleep(300);
      await note(`plots: dimming by the pointer: ${ps.set.members[3].label} hovered drawn ${w}px wide; ⌥-click: it alone, ${dimmed} others dimmed; ⌥-click again: all back ${back ? "yes" : "no"}: ${w === 2 && alone && dimmed === ps.set.members.length - 1 && back ? "OK" : "FAIL"}`);
    }
    // the order and the variable (L183): reversed, sorted back by the core, read from the names (the Tb files' hold none)
    const first = () => stateOf(t).set.members[0].label;
    await reorder("plots_reverse", t); await sleep(300); const rev = first();
    await reorder("plots_sort", t); await sleep(300); const back = first();
    const [, nn] = await invoke("plots_read_variable", { id: t.setId, source: "name", replaceTyped: false });
    await note(`plots: order: reversed, the first member ${rev}; sorted by the core, ${back}; read from the file names: ${nn} value(s) (the Tb files' names hold an index, "t01"): ${rev !== back && nn === 0 ? "OK" : "FAIL"}`);
    // PL20: a member drawn after the operations Data keeps for its file, and saying so. A fresh copy keeps none (L200):
    // the check then gives tb_t20.dat a factor of 2 through Data's own command and reads the members again
    let withOps = stateOf(t).set.members.filter((m) => m.operations.length), made = "";
    if (!withOps.length && state.project && state.project.root && writable()) {
      const m20 = stateOf(t).set.members.find((m) => m.file.endsWith("tb_t20.dat"));
      if (m20 && m20.path.startsWith(state.project.root + "/")) {
        try { await invoke("data_set_ops", { root: state.project.root, path: m20.path.slice(state.project.root.length + 1), ops: [{ op: "normalize", kind: "factor", value: 2.0 }] }); made = " (a factor of 2 given to tb_t20.dat by Data's command for the check)"; } catch (e) { made = ` (Data's command refused: ${e})`; }
        await reloadMembers(t); await sleep(300);
        withOps = stateOf(t).set.members.filter((m) => m.operations.length);
      }
    }
    await note(`plots: Data's operations${made}: ${withOps.length ? withOps.map((m) => `${m.label} drawn after ${m.operations.join("; ")}`).join(", ") : "none kept for these files"}; ${stateOf(t).set.members.length - withOps.length} member(s) as measured: ${withOps.length ? "OK" : writable() ? "FAIL" : "not checked (a reference folder: nothing written)"}`);
    // PL21: FullProf's temperatures file beside the buffer (Tb.tem: tb_t02.dat TSAMP 4.0, tb_t03.dat TSAMP 6.5)
    const tem = state.project.files.find((x) => !x.dir && /\.tem$/i.test(x.name));
    if (tem) {
      const r = await invoke("plots_temperatures", { id: t.setId, path: `${state.project.root}/${tem.path}` }); setSet(t, r.set); await sleep(300);
      const m3 = r.set.members.find((m) => m.file.endsWith("tb_t03.dat"));
      await note(`plots: temperatures file ${tem.name}: ${r.taken} member(s) took a value; tb_t03.dat ${m3 && m3.variable ? `${fmt(m3.variable.value)} K from ${m3.variable.source}` : "none"}: ${r.taken === 2 && m3 && m3.variable && m3.variable.value === 6.5 ? "OK" : "FAIL"}`);
    } else await note("plots: no temperatures file in the project: not checked");
    // PL12: BT typed with its source, the set saved as Plots/selfcheck.drv, opened again by the file's double-click path
    if (writable()) {
      const sp = stateOf(t); const bt = sp.set.apart.length ? sp.set.apart[0] : null;
      if (bt != null) setSet(t, await invoke("plots_set_variable", { id: t.setId, member: bt, value: 1.5, note: "WinPLOTR's demo page (self-check)" }));
      sp.wins = sp.wins.length ? sp.wins : [{ name: "m1", from: 18.5, to: 19.8 }];
      const rel = (await invoke("plots_save", { id: t.setId, root: state.project.root, name: "selfcheck", settings: settingsOf(sp) })).path;
      const reread = await reloadOthers(rel, t);
      bus.emit("file-saved", rel); await sleep(1500);
      // the file as written, read back by the core (what a later launch opens), whatever tab shows it
      let back = null, backSet = null;
      try { [backSet, back] = await invoke("plots_open_set", { path: `${state.project.root}/${rel}`, root: state.project.root }); await invoke("plots_close", { id: backSet.id }); } catch (e) { await note(`plots: ${rel} not read back: ${e}: FAIL`); }
      if (window.__openPath) await window.__openPath(rel, "plots"); await sleep(1200);
      const o = activePlots(); const os = o ? stateOf(o) : null;
      // tb_t03.dat's value and its source as the set was saved: 6.5 K from Tb.tem when the copy holds a temperatures file
      // (PL23 keeps it), the title's 6 K when it holds none (the merging session's copy, L384)
      const t03 = sp.set.members.find((m) => m.file.endsWith("tb_t03.dat"));
      const sameT03 = () => { const m = backSet.members.find((x) => x.file.endsWith("tb_t03.dat")); return !!(m && t03 && m.variable && t03.variable && m.variable.value === t03.variable.value && m.variable.source === t03.variable.source); };
      const fileOk = !!(back && backSet && backSet.members.length === sp.set.members.length && backSet.members.some((m) => m.file.endsWith("tb_t01.dat") && m.variable && m.variable.typed) && back.windows.length === sp.wins.length && back.prule === "step" && back.pdx === 0.22 && back.pdy === 2200 && back.pscale === "linear" && back.map_scale === sp.mapScale && sameT03());
      await note(`plots: ${rel} read back by the core: ${back ? `${backSet.members.length} members, ${back.windows.length} window(s), the projection's ${back.prule} Δx ${fmt(back.pdx)} Δy ${fmt(back.pdy)} on the ${back.pscale} scale, the map's ${back.map_scale}, the y scale ${back.scale}; tb_t03.dat ${(() => { const m = backSet.members.find((x) => x.file.endsWith("tb_t03.dat")); return m && m.variable ? `${fmt(m.variable.value)} K from ${m.variable.source} (PL23)` : "no T"; })()}` : "not read"}${reread ? `; ${reread} other tab of it (reopened at launch) read again` : ""}: ${fileOk ? "OK" : "FAIL"}`);
      const b2 = os && os.set ? os.set.members.find((m) => m.file.endsWith("tb_t01.dat")) : null;
      const t3 = os && os.set ? os.set.members.find((m) => m.file.endsWith("tb_t03.dat")) : null;
      const okSet = !!(o && o.drv && os.set && os.set.members.length === sp.set.members.length && b2 && b2.variable && b2.variable.typed && os.wins.length === sp.wins.length && os.prule.dx === 0.22 && os.prule.dy === 2200 && (os.prule.scale || os.scale) === "linear" && os.mapScale === sp.mapScale && o !== t);
      await note(`plots: saved ${rel} and opened it again: ${os && os.set ? `${os.set.members.length} members, tb_t01.dat ${b2 && b2.variable ? `${fmt(b2.variable.value)} K, ${b2.variable.source}` : "no T"}, ${os.wins.length} window(s), the projection's steps Δx ${fmt(os.prule.dx)}°, Δy ${fmt(os.prule.dy)} on the ${os.prule.scale || os.scale} scale (PL27: apart from the stack's, the core's proposed stack step not written), the map's colour scale ${os.mapScale}, the reference ${os.reference ? memberByKey(os, os.reference).label : "none"}; tb_t03.dat ${t3 && t3.variable ? `${fmt(t3.variable.value)} K from ${t3.variable.source}` : "no T"} (kept by the set file with its source, PL23)` : "not opened"}: ${okSet ? "OK" : "FAIL"}`);
      if (o && o.id !== t.id) { await sleep(300); activateTab(t.id); }
    }
    // step 5 (PL6, PL10, PL7, PL8): Normalise… and Convert X… as views, the CSV written and read back as a set
    if (activePlots() !== t) { activateTab(t.id); await sleep(400); }
    const exportView = async (name) => {
      if (!(state.project && state.project.root && writable())) return;
      const fig = inLightPalette(() => plots.figure()); if (fig && fig.after) fig.after();
      if (fig) { try { const out = await invoke("export_figure", { svg: fig.svg, path: `${state.project.root}/Exports/selfcheck plots ${name}.png`, dpi: 150 }); await note(`plots: ${name} exported: ${out}`); } catch (e) { await note(`plots: ${name} not exported: ${e}`); } }
    };
    {
      const sp = stateOf(t); const vw = views.get(t.id);
      const before = sp.set.members.find((m) => m.file.endsWith("tb_t10.dat"));
      // the members' patterns as the view shows them are asked of the core (plots_patterns, step 22), not read off the set view
      const patternsNow = async () => new Map((await invoke("plots_patterns", { id: t.setId, members: null })).map((p) => [p.member, p]));
      const p0 = await patternsNow(); const y0 = before && p0.get(before.index) ? p0.get(before.index).y[100] : null;
      await applyViewRule(t, { by: "area", to: 1, from: 20, until: 30, member: null, background: null }, null); await sleep(300);
      const p1 = await patternsNow();
      sp.view = "overlay"; sp.linked = false; renderActive(); await sleep(600);
      let s = sp.set; const fs = s.members.map((m) => m.shown.factor).filter((f) => f != null);
      const m10 = s.members.find((m) => m.file.endsWith("tb_t10.dat"));
      const said = vw ? vw.note() : "";
      const series = vw ? vw.plot.state.series.filter((x) => x.kind !== "bands").length : 0;
      await note(`plots: Normalise… by the same area over 20–30°: ${fs.length} of ${s.members.length} members normalised, factors ${fs.length ? `${fmt(fs[0])} (${s.members[0].label}) … ${fmt(fs[fs.length - 1])} (${s.members[s.members.length - 1].label})` : "none"}; tb_t10.dat y[100] ${fmt(y0)} shown as ${m10 && p1.get(m10.index) ? fmt(p1.get(m10.index).y[100]) : "?"} (× ${m10 ? fmt(m10.shown.factor) : "?"}); ${series} curves; the note: "${said.slice(0, 110)}…": ${fs.length === s.members.length && /Normalised by the core/.test(said) && series === s.members.length ? "OK" : "FAIL"}`);
      await exportView("normalised overlay");
      await applyViewRule(t, { by: "monitor", to: 1, member: null, background: null, from: null, until: null }, null); await sleep(300);
      s = sp.set; const out = s.members.filter((m) => m.shown.left_out);
      await note(`plots: Normalise… by the monitor: ${out.length} of ${s.members.length} members left out, the core's reason "${out.length ? out[0].shown.left_out : ""}": ${out.length === s.members.length ? "OK" : "FAIL"}`);
      // PL10: in d, the Tb files state no λ: every member listed; with λ typed, the map in d
      await applyViewRule(t, null, { quantity: "d", control: null, typed: null }); await sleep(300);
      const without = sp.set.without_lambda.length;
      await applyViewRule(t, { by: "area", to: 1, from: 20, until: 30, member: null, background: null }, { quantity: "d", control: null, typed: 1.5406 }); await sleep(300);
      s = sp.set;
      sp.view = "map"; sp.map = null; renderActive();
      const okMap = await wait(() => !!(sp.map && sp.map.key.startsWith("m:"))); await sleep(400);
      const M = okMap ? sp.map.v.map : null;
      await note(`plots: Convert X… to d: ${without} of ${s.members.length} members listed without a λ; with λ 1.5406 Å typed, ${s.members.filter((m) => m.shown.lambda != null).length} converted (${s.members[0].shown.lambda_from}); the map in ${s.x_label}, ${M ? `${M.rows.length} rows × ${M.x.length} columns, d from ${fmt(M.x[0])} to ${fmt(M.x[M.x.length - 1])} Å` : "none"}; the note says "${(vw ? vw.note() : "").slice(0, 90)}…": ${without === s.members.length && M && M.rows.length && M.x[0] < M.x[M.x.length - 1] && /x in d/.test(vw ? vw.note() : "") ? "OK" : "FAIL"}`);
      await exportView("normalised map in d");
      // PL7 and PL8: the CSV in d, wide and long, written by the core and read back as sets
      if (state.project && state.project.root && writable()) {
        const back = [];
        // wide in d, long in 2θ: the core's reader takes a CSV in d back as 2θ (PL26)
        for (const [long, x] of [[false, "d"], [true, "2theta"]]) {
          const path = `${state.project.root}/Exports/selfcheck plots${x === "d" ? " in d" : ""}${long ? " long" : ""}.csv`;
          try {
            const r = await invoke("plots_export_csv", { id: t.setId, path, long, x, wavelength: 1.5406, range: "overlap" });
            const b = await invoke("plots_open_files", { paths: [path], root: state.project.root });
            const t02 = b.members.find((m) => m.label === "tb_t02.dat");
            back.push(`${long ? "long" : "wide"} in ${x}: ${r.members} members, ${r.lines} lines; read back as a set of ${b.members.length} members, ${b.apart.length} apart, tb_t02.dat ${t02 && t02.variable ? `${fmt(t02.variable.value)} K from ${t02.variable.source}` : "no T"}, x ${b.members[0] ? b.members[0].axis : "?"}`);
            await invoke("plots_close", { id: b.id });
            if (b.members.length !== r.members) back.push("FAIL");
          } catch (e) { back.push(`${long ? "long" : "wide"}: ${e} FAIL`); }
        }
        await note(`plots: Export ▸ Plot Set as CSV (wide in d, long in 2θ), ${back.join("; ")}: ${back.some((x) => /FAIL/.test(x)) ? "FAIL" : "OK"}`);
      } else await note("plots: the CSV is not checked (no project of one's own)");
      await applyViewRule(t, null, null); sp.view = "overlay"; renderActive(); await sleep(300);
    }
    // step 16 (I08, R33): the clip on the colour bar and the row and column cursors, on the Tb₂Pd₂Sn map
    if (activePlots() !== t) { activateTab(t.id); await sleep(400); }
    {
      const sp = stateOf(t); const vw = views.get(t.id);
      const surfStats = () => { const sf = vw && vw.surface(); return sf && sf.stats() ? sf.stats() : null; };
      const image = () => (vw ? vw.plot.state.series.find((x) => x.kind === "image") : null);
      sp.view = "map"; sp.look = false; sp.cursors = true; sp.pin = null; sp.cell = null; sp.clip = { map: null, difference: null }; sp.map = null; renderActive();
      const okM = await wait(() => !!(sp.map && sp.map.key.startsWith("m:"))); await sleep(400);
      const M = okM ? sp.map.v.map : null;
      let clipOk = false, what = "no map";
      if (M) {
        // the two marks dragged to 0.4 and 0.8 of the bar's height (the shell's place on it), through the bar's own callback
        const from = barValue(0.4, M.range[0], M.range[1]), to = barValue(0.8, M.range[0], M.range[1]);
        const img0 = image(); if (img0 && img0.bar.clip) img0.bar.clip.onChange(from, to); await sleep(400);
        const img = image();
        const clipped = !!(img && img.bar.clip && img.bar.clip.from === from && img.bar.clip.to === to && sp.clip.map && sp.clip.map.scale === sp.map.scale);
        const said = /Colour range shown: .* saturated/.test(vw.note());
        await exportView("map clipped");
        // the surface's height on the same pair
        sp.view = "surface"; sp.surfOf = "map"; renderActive(); await wait(() => !!surfStats()); await sleep(400);
        const S = surfStats(); const surfSame = !!(S && S.range && S.range[0] === from && S.range[1] === to);
        // the difference map keeps a clip of its own (none set here)
        sp.view = "difference"; renderActive(); await wait(() => !!(sp.map && sp.map.key.startsWith("d:"))); await sleep(400);
        const imgD = image(); const ownD = !!(imgD && imgD.bar.clip && imgD.bar.clip.from == null && imgD.bar.clip.to == null);
        // a double-click on the bar: the core's range again
        sp.view = "map"; renderActive(); await wait(() => !!(sp.map && sp.map.key.startsWith("m:"))); await sleep(400);
        const img2 = image(); if (img2 && img2.bar.clip) img2.bar.clip.onReset(); await sleep(300);
        const img3 = image(); const reset = !!(img3 && img3.bar.clip && img3.bar.clip.from == null && !sp.clip.map);
        clipOk = clipped && said && surfSame && ownD && reset;
        what = `${fmt(from)} to ${fmt(to)} of the core's ${M.range[0].toFixed(3)} to ${M.range[1].toFixed(3)} on the ${sp.map.scale} scale: the map ${clipped ? "clipped" : "not clipped"}, the note ${said ? "says the saturation" : "does not say it"}, the surface's height ${surfSame ? "on the same pair" : `not (${S && S.range ? S.range.map((x) => fmt(x)).join(" to ") : "no range"})`}, the difference map's own clip ${ownD ? "none" : "not its own"}, a double-click on the bar ${reset ? "returns to the core's range" : "does not"}`;
      }
      await note(`plots: the clip on the colour bar (I08 step 1): ${what}: ${clipOk ? "OK" : "FAIL"}`);
      // the cursors: the cell at 22 K, 19.1° (the strongest magnetic peak), its two strips, a pin and its release
      let curOk = false; what = "no map";
      if (M) {
        const ri = M.rows.findIndex((row) => row.value === 22), cj = M.x.findIndex((x) => x >= 19.1);
        if (ri >= 0 && cj >= 0) {
          vw.cell({ id: "map", j: cj, r: ri }); await sleep(400);
          const st2 = vw.strips(); const rs = st2.row && st2.row.state.series[0], cs = st2.col && st2.col.state.series[0];
          const m22 = sp.set.members[M.rows[ri].member];
          // the strip's pattern is the one the viewer's data holds for the member (asked of the core once, step 22)
          const patOf = (i) => { const d = vw.viewerData(); return d.members[i] ? d.members[i].pattern : null; };
          const p22 = patOf(M.rows[ri].member);
          const rowOk = !!(rs && p22 && rs.x === p22.x && rs.y === p22.y && rs.label.includes(m22.label));
          // step 23: the column strip's path has the variable up (y) and the core's values across (x), the rows holding a value only
          const k22 = cs ? cs.y.indexOf(M.rows[ri].value) : -1;
          const colOk = !!(cs && k22 >= 0 && cs.x[k22] === M.z[ri][cj] && cs.x.length === M.z.filter((zr) => zr[cj] != null).length);
          const cc = vw.plot.state.cellCursor; const marked = !!(cc && cc.j === cj && cc.r === ri);
          const rcv = document.querySelector(".pl-rcanvas"), ccv = document.querySelector(".pl-ccanvas"), areaEl = document.querySelector(".pl-areas");
          const rr = rcv ? rcv.getBoundingClientRect() : { height: 0 }, cr = ccv ? ccv.getBoundingClientRect() : { width: 0 }, ar = areaEl ? areaEl.getBoundingClientRect() : { width: 0, height: 0 };
          // the strips need room (the module's container queries): the row strip a plot area over 260 px tall, the
          // column strip one over 520 px wide; a relaunch reopens the window with a short area (step 16's relaunch, 100 px)
          const roomRow = ar.height > 260, roomCol = ar.width > 520;
          vw.cellClick({ id: "map", j: cj, r: ri }); await sleep(300); const pinned = !!(sp.pin && sp.pin.r === ri && sp.pin.j === cj);
          vw.cell({ id: "map", j: cj + 5, r: 0 }); await sleep(300);
          const held = !!(st2.row.state.series[0] && p22 && st2.row.state.series[0].x === p22.x);
          await exportView("map cursors");
          vw.cellClick({ id: "map", j: cj + 5, r: 0 }); await sleep(400);
          const m0 = sp.set.members[M.rows[0].member];
          const p0r = patOf(M.rows[0].member); void m0;
          const released = !sp.pin && !!(st2.row.state.series[0] && p0r && st2.row.state.series[0].x === p0r.x);
          curOk = rowOk && colOk && marked && pinned && held && released && (rr.height > 40 || !roomRow) && (cr.width > 40 || !roomCol);
          what = `over ${fmt(M.rows[ri].value)} K, ${M.x[cj].toFixed(2)}°: the row strip ${rowOk ? `draws ${m22.label}'s pattern as the set holds it, ${rs.x.length} points` : "does not draw the member's pattern"} (${roomRow ? `${Math.round(rr.height)} px high` : `not shown: the plot area is ${Math.round(ar.height)} px tall, under the 260 px it needs`}); the column strip ${colOk ? `${cs.x.length} values against T, at 22 K ${cs.x[k22].toFixed(3)} (the core's)` : "does not draw the core's column"} (${roomCol ? `${Math.round(cr.width)} px wide` : `not shown: the plot area is ${Math.round(ar.width)} px wide, under the 520 px it needs`}); the cell marked on the map ${marked ? "yes" : "no"}; a click pins ${pinned ? "yes" : "no"}, the strips held while the pointer moves ${held ? "yes" : "no"}, a click again releases ${released ? "yes" : "no"}`;
        } else what = `no row at 22 K or no column at 19.1° (${ri}, ${cj})`;
      }
      await note(`plots: the row and column cursors (I08 step 3): ${what}: ${curOk ? "OK" : "FAIL"}`);
      // step 23 (L431): the strips' frames are the map's — the row strip's left and right and its x view, the column
      // strip's top and bottom and its y view — at the whole map and in a zoom; each strip's own label carries no scale
      // word (the plot says it once). A strip without room (the container queries) has no frame to compare: said so.
      let frOk = false; what = "no map";
      if (M) {
        const ri = M.rows.findIndex((row) => row.value === 22), cj = M.x.findIndex((x) => x >= 19.1);
        // the check's window is short (the relaunch's ~100 px, the first run's 243 px): the plot area is given room for
        // the check (a minimum of 960 × 640 px, overflowing the window) so that both strips are shown, and given back after
        const areaEl = document.querySelector(".pl-areas"); const was = areaEl ? { h: areaEl.style.minHeight, w: areaEl.style.minWidth } : null;
        if (areaEl) { areaEl.style.minHeight = "640px"; areaEl.style.minWidth = "960px"; await sleep(500); }
        vw.cellClick({ id: "map", j: cj, r: ri }); await sleep(400);   // pinned, so that the strips stay while the view changes
        const st3 = vw.strips(); const ar = areaEl ? areaEl.getBoundingClientRect() : { width: 0, height: 0 };
        const roomRow = ar.height > 260, roomCol = ar.width > 520;
        const same = () => {
          const g = vw.plot.state.geom, rg = st3.row ? st3.row.state.geom : null, cg = st3.col ? st3.col.state.geom : null;
          const mv = vw.plot.getView(), rv = st3.row ? st3.row.getView() : null, cv = st3.col ? st3.col.getView() : null;
          const row = !!(g && rg && rg.left === g.left && rg.right === g.right && rv.x0 === mv.x0 && rv.x1 === mv.x1);
          const col = !!(g && cg && cg.top === g.top && cg.bottom === g.bottom && cv.y0 === g.y0 && cv.y1 === g.y1);
          return { g, rg, cg, row, col, mv };
        };
        const whole = same();
        vw.plot.setView({ x0: 18.95, x1: 19.25, y0: 20, y1: 26 }); await sleep(400);
        const zoomed = same(); await exportView("map strips framed");
        vw.plot.fitAll(); await sleep(300);
        const rl = st3.row ? st3.row.state.yLabel : "", cl = st3.col ? st3.col.state.xLabel : "";
        const once = !/\(√\)|\(log\)/.test(rl) && !/\(√\)|\(log\)/.test(cl);   // the plot adds " (√)" or " (log)" to its y label itself, so the strip's own carries none
        const rowOk = !roomRow || (whole.row && zoomed.row), colOk = !roomCol || (whole.col && zoomed.col);
        frOk = rowOk && colOk && once && !!whole.g;
        const edges = (q) => (q ? `${q.left}–${q.w - q.right} px across, ${q.top}–${q.bottom} px down` : "none");
        what = `the plot area given ${Math.round(ar.width)} × ${Math.round(ar.height)} px for the check; the map's frame ${edges(whole.g)}; the row strip ${roomRow ? (rowOk ? `the same left and right (${whole.rg.left}–${whole.rg.w - whole.rg.right} px) and the same x view, whole and zoomed to 18.95–19.25°` : `not the map's (${edges(whole.rg)}, zoomed ${edges(zoomed.rg)}; x ${whole.rg ? `${st3.row.getView().x0}–${st3.row.getView().x1}` : "?"})`) : `not shown (the plot area is ${Math.round(ar.height)} px tall, under the 260 px it needs)`}; the column strip ${roomCol ? (colOk ? `the same top and bottom (${whole.cg.top}–${whole.cg.bottom} px) and the same y view, whole (${fmt(whole.g.y0)}–${fmt(whole.g.y1)} K) and zoomed to 20–26 K` : `not the map's (${edges(whole.cg)}, zoomed ${edges(zoomed.cg)})`) : `not shown (the plot area is ${Math.round(ar.width)} px wide, under the 520 px it needs)`}; the labels "${rl}" and "${cl}" ${once ? "say the unit once" : "repeat the scale"}`;
        sp.pin = null; vw.cellClick({ id: "map", j: cj, r: ri }); await sleep(200); sp.pin = null;
        if (areaEl) { areaEl.style.minHeight = was.h; areaEl.style.minWidth = was.w; await sleep(300); }
      }
      await note(`plots: the strips' frames are the map's (L431; step 23): ${what}: ${frOk ? "OK" : "FAIL"}`);
      sp.pin = null; sp.view = "overlay"; renderActive(); await sleep(300);
    }
    // step 17 (I08 step 4, R33; PL37 met by core 0.76.0): the zoom in the T–2θ plane kept across redraws, the wheel on
    // each axis, "range from the view" against drs plots region's numbers (L407), a region saved, the clip and the
    // region written to the set file, read back by the core, opened again and the region chosen; a double-click back
    if (activePlots() !== t) { activateTab(t.id); await sleep(400); }
    {
      const sp = stateOf(t); const vw = views.get(t.id); const image = () => vw.plot.state.series.find((x) => x.kind === "image");
      sp.view = "map"; sp.look = false; sp.cursors = false; sp.pin = null; sp.rangeFromView = false; sp.clip = { map: null, difference: null }; sp.clipBefore = { map: null, difference: null }; sp.regions = []; sp.mapZoom = null; sp.mapScale = "linear"; sp.map = null; renderActive();
      const okM = await wait(() => !!(sp.map && sp.map.key.startsWith("m:") && sp.map.scale === "linear")); await sleep(400);
      const M = okM ? sp.map.v.map : null;
      let what = "no map on the linear scale", zoomOk = false, fileOk = false;
      if (M) {
        // a rectangle dragged (the plot's box zoom ends in the view set): 18.95–19.25° by 20–26 K, L407's region
        vw.plot.setView({ x0: 18.95, x1: 19.25, y0: 20, y1: 26 }); await sleep(300);
        const z1 = vw.plot.getView(); const zoomed = !!(sp.mapZoom && z1.x0 === 18.95 && z1.y0 === 20 && z1.y1 === 26);
        // the zoom kept across a redraw for the clip (step 16's drag through the bar's callback)
        const lo = barValue(0.3, M.range[0], M.range[1]), hi = barValue(0.9, M.range[0], M.range[1]);
        image().bar.clip.onChange(lo, hi); await sleep(300);
        const z2 = vw.plot.getView(); const kept = z2.x0 === 18.95 && z2.x1 === 19.25 && z2.y0 === 20 && z2.y1 === 26 && !!sp.clip.map;
        // the wheel over each axis (a wheel event at the axis): the y axis's margin zooms T alone, the x axis's band x alone
        const cv = document.querySelector(".pl-area canvas"); const g = vw.plot.state.geom; const cr = cv.getBoundingClientRect();
        const wheel = (x, y) => cv.dispatchEvent(new WheelEvent("wheel", { clientX: cr.left + x, clientY: cr.top + y, deltaY: -240, bubbles: true, cancelable: true }));
        wheel(g.left - 12, (g.top + g.bottom) / 2); await sleep(250); const z3 = vw.plot.getView();
        const yZoomed = z3.y1 - z3.y0 < z2.y1 - z2.y0 - 1e-9 && z3.x0 === z2.x0 && z3.x1 === z2.x1;
        wheel((g.left + g.w - g.right) / 2, g.bottom + 12); await sleep(250); const z4 = vw.plot.getView();
        const xZoomed = z4.x1 - z4.x0 < z3.x1 - z3.x0 - 1e-9 && z4.y0 === z3.y0 && z4.y1 === z3.y1;
        // the rectangle again, then Range from the View: the clip from the core's block, as drs plots region gives it for
        // L407's region (rows 8–14, columns 30–32 of the map, its cells 2718 to 45849 counts)
        vw.plot.setView({ x0: 18.95, x1: 19.25, y0: 20, y1: 26 }); await sleep(250);
        toggleFollow(t); const gotB = await wait(() => !!(sp.clip.map && sp.clip.map.block)); await sleep(300);
        const B = gotB ? sp.clip.map.block : null, c = sp.clip.map;
        // the block is the rectangle's cells: every row of the map within 20–26 K and no other, three columns within
        // 18.95–19.25°, and the cells' range as drs plots region gives it on the buffer alone (L407: 2718 to 45849; its
        // indices are 8–14 and 30–32 there, the check's set holding the typed members too, so its indices differ)
        const inT = (i) => M.rows[i].value >= 20 && M.rows[i].value <= 26;
        const followOk = !!(B && B.rows.length === 7 && M.rows.every((row, i) => inT(i) === B.rows.includes(i)) && B.columns.length === 3 && B.columns.every((j) => M.x[j] >= 18.95 && M.x[j] <= 19.25) && B.range[0] > 2717 && B.range[0] < 2719 && B.range[1] > 45848 && B.range[1] < 45850 && c.from === B.range[0] && c.to === B.range[1] && /from the view/.test(vw.note()));
        await exportView("map region range from the view");
        // a region saved with that range; the switch off puts the dragged clip back
        const region = await addRegion(t, "m1 transition", { x: [18.95, 19.25], variable: [20, 26] }, c ? { from: c.from, to: c.to } : null);
        toggleFollow(t); await sleep(250);
        const back0 = !sp.rangeFromView && !!sp.clip.map && sp.clip.map.from === lo && sp.clip.map.to === hi;
        // a double-click: the whole map again, the clip kept
        vw.plot.fitAll(); await sleep(250); const whole = sp.mapZoom == null && !!sp.clip.map;
        zoomOk = zoomed && kept && yZoomed && xZoomed && followOk && !!region && back0 && whole;
        what = `a rectangle 18.95–19.25° × 20–26 K: zoomed ${zoomed ? "yes" : "no"}, kept across the clip's redraw ${kept ? "yes" : "no"}; the wheel over the y axis zooms T alone ${yZoomed ? "yes" : "no"} (${fmt(z3.y0)}–${fmt(z3.y1)} K), over the x axis x alone ${xZoomed ? "yes" : "no"} (${fmt(z4.x0)}–${fmt(z4.x1)}°); Range from the View: ${B ? `the core's block rows ${B.rows[0] + 1}–${B.rows[B.rows.length - 1] + 1}, columns ${B.columns[0] + 1}–${B.columns[B.columns.length - 1] + 1}, cells ${fmt(B.range[0])} to ${fmt(B.range[1])} (the rows within 20–26 K, the columns within 18.95–19.25°; drs plots region on the buffer alone, L407: 2718 to 45849), the clip ${c ? `${fmt(c.from)} to ${fmt(c.to)}` : "none"}` : "no block"} ${followOk ? "as the core's" : "FAIL"}; region "m1 transition" ${region ? "accepted by the core" : "refused"}; the switch off puts the clip back ${back0 ? "yes" : "no"}; a double-click shows the whole map ${whole ? "yes" : "no"}`;
        // the clip and the region written by the core, read back, the set opened again and the region chosen there
        if (state.project && state.project.root && writable()) {
          const root = state.project.root;
          t.saveName = "selfcheck regions"; await saveSet(t); await sleep(400);
          let back = null, err = "";
          try { const [bs, bst] = await invoke("plots_open_set", { path: `${root}/${t.drv}`, root }); back = bst; await invoke("plots_close", { id: bs.id }); } catch (e) { err = String(e); }
          const read = !!(back && back.map_range && Math.abs(back.map_range.from - sp.clip.map.from) < 1e-9 && Math.abs(back.map_range.to - sp.clip.map.to) < 1e-9 && back.map_scale === "linear" && back.regions.length === 1 && back.regions[0].name === "m1 transition" && back.regions[0].x_quantity === "2theta" && back.regions[0].map_range && !back.regions_refused);
          if (window.__openPath) await window.__openPath(t.drv, "plots"); await sleep(1200);
          const o = activePlots(); const os = o ? stateOf(o) : null;
          const reopened = !!(o && o !== t && os.clip.map && os.clip.map.from === sp.clip.map.from && os.regions.length === 1 && os.mapScale === "linear");
          let chosenOk = false, cw = "";
          if (reopened) {
            os.view = "map"; renderActive(); await wait(() => !!(os.map && os.map.key.startsWith("m:"))); await sleep(400);
            const ok2 = await chooseRegion(o, 0); await sleep(500);
            const ov = views.get(o.id).plot.getView();
            // the region's own colour range (the view's when it was saved), not the set's clip, is what choosing it applies
            const rm = region && region.map_range;
            chosenOk = ok2 && !!os.mapZoom && Math.abs(ov.x0 - 18.95) < 1e-9 && Math.abs(ov.x1 - 19.25) < 1e-9 && ov.y0 === 20 && ov.y1 === 26 && !!os.clip.map && !!rm && os.clip.map.from === rm.from && os.clip.map.to === rm.to && /Colour range shown/.test(views.get(o.id).note());
            cw = `zoomed to ${fmt(ov.x0)}–${fmt(ov.x1)}°, ${fmt(ov.y0)}–${fmt(ov.y1)} K, its own colour range ${os.clip.map ? `${fmt(os.clip.map.from)} to ${fmt(os.clip.map.to)}` : "none"} applied${rm ? ` (saved ${fmt(rm.from)} to ${fmt(rm.to)})` : ""}`;
            views.get(o.id).plot.fitAll(); await sleep(200);
          }
          fileOk = read && reopened && chosenOk;
          await note(`plots: the clip and the region in the set file (PL37): ${t.drv} written with the map's colour range ${fmt(sp.clip.map.from)} to ${fmt(sp.clip.map.to)} on the linear scale and 1 region; read back by the core ${back ? `${back.map_range ? `map_range ${fmt(back.map_range.from)} to ${fmt(back.map_range.to)}` : "no map_range"}, ${back.regions.length} region(s)${back.regions_refused ? `, refused: ${back.regions_refused}` : ""}` : `not read: ${err}`} ${read ? "as saved" : "FAIL"}; opened again ${reopened ? "with its clip and region restored" : "without them"}; the region chosen there ${chosenOk ? cw : `not as saved (${cw || "not chosen"})`}: ${fileOk ? "OK" : "FAIL"}`);
          if (o && o.id !== t.id) { activateTab(t.id); await sleep(300); }
        } else { fileOk = true; await note("plots: the clip and the region in the set file: not checked (no project of one's own)"); }
      }
      await note(`plots: the zoom in the T–2θ plane and range from the view (I08 step 4): ${what}: ${zoomOk ? "OK" : "FAIL"}`);
      sp.regions = []; sp.selregion = null; sp.clip = { map: null, difference: null }; sp.mapZoom = null; sp.rangeFromView = false; sp.mapScale = "log"; sp.map = null; sp.cursors = true; sp.view = "overlay"; renderActive(); await sleep(300);
    }
    // step 18 (I08 step 5, L409): the Map view through the factored viewer (steps 16 and 17 above ran through it), and
    // the viewer's page written from this set as drs plots show will write it (PL38): the JSON and the bundle inside,
    // nothing fetched from the network (the SVG namespace is the one http the sources name)
    if (activePlots() !== t) { activateTab(t.id); await sleep(400); }
    {
      const sp = stateOf(t); const vw = views.get(t.id);
      sp.view = "map"; renderActive(); await wait(() => !!(sp.map && sp.map.key.startsWith("m:"))); await sleep(300);
      const through = !!(vw && vw.viewer && vw.viewer.which() === "map" && vw.plot.state.series.some((x) => x.kind === "image" && x.id === "map"));
      let pageOk = false, pw = "not written (no project of one's own)";
      if (state.project && state.project.root && writable()) {
        sp.regions = [{ name: "m1 transition", x: [18.95, 19.25], x_quantity: "2theta", variable: [20, 26], map_range: null, window: null }];
        const p = await exportViewerPage(t, `${state.project.root}/Exports/selfcheck viewer.html`); sp.regions = []; renderActive();
        if (p) {
          const h = p.html;
          const hasData = h.includes('id="drs-page-data"') && /"maps":\{"map":\{"map":\{/.test(h) && h.includes('"regions":[{"name":"m1 transition"') && /"blocks":\{"map":\{"rows":\[/.test(h) && /"members":\[\{"label":"tb_t01\.dat"/.test(h);
          const hasViewer = h.includes("diff-ref-struct Map viewer, one classic script") && h.includes("function createMapViewer(") && !/<!--DRS_/.test(h);
          const noHttp = !/\bhttps?:\/\/(?!www\.w3\.org\/)/.test(h);
          pageOk = hasData && hasViewer && noHttp && p.bytes > 100000;
          pw = `${p.path.split("/").pop()}, ${Math.round(p.bytes / 1024)} kB: the JSON ${hasData ? "with the maps, the members and the region's block" : "missing or short"}, the viewer ${hasViewer ? "bundled, the markers replaced" : "missing"}, network references ${noHttp ? "none" : "FOUND"}`;
        } else pw = "not written (the status says why)";
      } else pageOk = true;
      await note(`plots: the Map view through the factored viewer (step 18): ${through ? "the viewer draws the map on the shared plot" : "not the viewer's"}; the viewer's page ${pw}: ${through && pageOk ? "OK" : "FAIL"}`);
      sp.view = "overlay"; renderActive(); await sleep(200);
    }
    // step 19 (P65): a region's own colour range per view: one region saved on the map carries the map's clip in
    // map_range, one saved on the difference map carries the difference's in difference_range; chosen on either view,
    // each applies the range of the view shown only (a range on log₁₀ counts means nothing on the signed √ difference)
    if (activePlots() !== t) { activateTab(t.id); await sleep(400); }
    {
      const sp = stateOf(t); const vw = views.get(t.id); const image = () => vw.plot.state.series.find((x) => x.kind === "image");
      sp.regions = []; sp.clip = { map: null, difference: null }; sp.mapZoom = null; sp.rangeFromView = false; sp.mapScale = "log"; sp.map = null; sp.cursors = false;
      sp.view = "map"; renderActive(); const okM = await wait(() => !!(sp.map && sp.map.key.startsWith("m:"))); await sleep(300);
      let ok = false, what = "no map";
      if (okM) {
        const M = sp.map.v.map, mLo = barValue(0.4, M.range[0], M.range[1]), mHi = barValue(0.8, M.range[0], M.range[1]);
        image().bar.clip.onChange(mLo, mHi); await sleep(250);
        const rm = await addRegion(t, "p65 map", { x: [18.95, 19.25], variable: [20, 26] }, { from: mLo, to: mHi });
        sp.view = "difference"; renderActive(); const okD = await wait(() => !!(sp.map && sp.map.key.startsWith("d:"))); await sleep(300);
        const D = okD ? sp.map.v.map : null, dLo = D ? barValue(0.3, D.range[0], D.range[1]) : 0, dHi = D ? barValue(0.7, D.range[0], D.range[1]) : 1;
        if (D) image().bar.clip.onChange(dLo, dHi); await sleep(250);
        const rd = D ? await addRegion(t, "p65 difference", { x: [26.0, 27.4], variable: [4, 12] }, { from: dLo, to: dHi }) : null;
        const saved = !!(rm && rm.map_range && !rm.difference_range && rm.map_range.from === mLo && rd && rd.difference_range && !rd.map_range && rd.difference_range.from === dLo);
        // on the difference map: the map's region brings no range, the difference's its own
        sp.clip.difference = null; renderTab(t); await sleep(200);
        await chooseRegion(t, 0); await sleep(300); const noneD = !sp.clip.difference;
        await chooseRegion(t, 1); await sleep(300); const ownD = !!(sp.clip.difference && sp.clip.difference.from === dLo && sp.clip.difference.to === dHi);
        // on the map: the difference's region brings no range, the map's its own
        sp.view = "map"; sp.clip.map = null; renderActive(); await wait(() => !!(sp.map && sp.map.key.startsWith("m:"))); await sleep(300);
        await chooseRegion(t, 1); await sleep(300); const noneM = !sp.clip.map;
        await chooseRegion(t, 0); await sleep(300); const ownM = !!(sp.clip.map && sp.clip.map.from === mLo && sp.clip.map.to === mHi);
        const sv = settingsOf(sp);
        const written = !!(sv.regions[0] && sv.regions[0].map_range && !sv.regions[0].difference_range && sv.regions[1] && sv.regions[1].difference_range && !sv.regions[1].map_range);
        ok = saved && noneD && ownD && noneM && ownM && written;
        what = `"p65 map" saved on the map with ${fmt(mLo)} to ${fmt(mHi)} (log₁₀ counts), "p65 difference" on the difference map with ${fmt(dLo)} to ${fmt(dHi)} (signed √): ${saved ? "each in its view's field only" : "not as saved"}; on the difference map the map's region applies nothing ${noneD ? "yes" : "no"} and the difference's its own ${ownD ? "yes" : "no"}; on the map the difference's applies nothing ${noneM ? "yes" : "no"} and the map's its own ${ownM ? "yes" : "no"}; the settings carry each in its field ${written ? "yes" : "no"} (the difference's waits for the core's Region, P65)`;
      }
      await note(`plots: a region's own colour range per view (P65, step 19): ${what}: ${ok ? "OK" : "FAIL"}`);
      sp.regions = []; sp.clip = { map: null, difference: null }; sp.mapZoom = null; sp.cursors = true; sp.view = "overlay"; renderActive(); await sleep(200);
    }
    // step 20 (W4, I10 step 0): the y interval shared by the set or adapted (the strip under the map on the core's
    // extent or the row's own; Overlay and Stack fitted to every member or to the one chosen); the members chosen all,
    // none, by exposure class, by a range of the variable; the list's order; the icons and the facts; the right-hand
    // axis, absent until the core gives each row's second variable (W3, core 0.79.0)
    if (activePlots() !== t) { activateTab(t.id); await sleep(400); }
    {
      const sp = stateOf(t); const vw = views.get(t.id); const n = sp.set.members.length;
      sp.view = "map"; sp.cursors = true; sp.pin = null; sp.cell = null; sp.yInterval = "shared"; sp.mapScale = "log"; sp.map = null; sp.extent = null; sp.sel = []; renderActive();
      const okM = await wait(() => !!(sp.map && sp.map.key.startsWith("m:")) && !!(sp.extent && sp.extent.v)); await sleep(300);
      let what = "no map or no extent from the core", ok = false;
      if (okM) {
        const M = sp.map.v.map, ri = M.rows.findIndex((r) => r.value === 22), cj = M.x.findIndex((x) => x >= 19.1);
        const img = vw.plot.state.series.find((x) => x.kind === "image");
        const rightAxis = img && img.rightLabels ? `${img.rightLabels.length} labels (${img.rightTitle})` : "absent (no variables per row from the core: W3, core 0.79.0)";
        // the strip under the map: shared, the set's extent on the map's scale (plots_extent); own, the row's fit
        vw.cell({ id: "map", j: cj, r: ri }); await sleep(300);
        const ext = sp.extent.v, strips = vw.strips(), rv1 = strips.row ? strips.row.getView() : null;
        const sharedOk = !!(ext && rv1 && rv1.y0 === ext.from && rv1.y1 === ext.to);
        toggleYInterval(t); await sleep(300); vw.cell({ id: "map", j: cj + 1, r: ri }); await sleep(300);
        const rv2 = strips.row ? strips.row.getView() : null; const ownOk = !!(rv2 && rv2.y0 == null);
        // Overlay and Stack: fitted to every member, then to the 22 K member chosen (a lower top than the set's)
        const key22 = sp.set.members[M.rows[ri].member].key;
        sp.view = "overlay"; sp.linked = false; sp.look = false; sp.yInterval = "shared"; sp.sel = []; renderActive(); await sleep(400);
        const yAll = vw.plot.state.geom ? vw.plot.state.geom.y1 : null;
        sp.sel = [key22]; sp.yInterval = "own"; renderActive(); await sleep(400);
        const yOne = vw.plot.state.geom ? vw.plot.state.geom.y1 : null;
        const ovOk = yAll != null && yOne != null && yOne < yAll;
        sp.view = "stack"; sp.scale = "linear"; sp.yInterval = "shared"; sp.sel = []; sp.rule = { kind: "step", dx: 0, dy: null, proposed: null }; sp.disp = null; renderActive(); await wait(() => !!(sp.disp && sp.disp.v)); await sleep(400);
        const sAll = vw.plot.state.geom ? vw.plot.state.geom.y1 : null;
        sp.sel = [key22]; sp.yInterval = "own"; renderActive(); await sleep(400);
        const sOne = vw.plot.state.geom ? vw.plot.state.geom.y1 : null;
        const stOk = sAll != null && sOne != null && sOne < sAll;
        sp.yInterval = "shared";
        // the members chosen: all, none, by exposure class, by a range of the variable (20–26 K: seven members)
        pickAll(t, true); await sleep(150); const all = sp.sel.length === n;
        pickAll(t, false); await sleep(150); const none = sp.sel.length === 0;
        const classes = exposureClasses(sp.set.members.map((m) => ({ counting_time: m.counting_time })));
        await pickClass(t, 0); await sleep(150); const cls = classes.length >= 1 && sp.sel.length === classes[0].members.length;
        await pickRange(t, 20, 26); await sleep(150); const rng = sp.sel.length === 7;
        // the list's order: by name first tb_t01.dat; by the file's time a permutation with every file timed; by the
        // variable never decreasing; the icons by the file's kind and the facts on hover in the sidebar
        const facts = sp.set.members.map((m) => ({ label: m.label, variable: m.variable, time: m.modified }));
        const byName = sortMembers(facts, "name"), byTime = sortMembers(facts, "time"), byVar = sortMembers(facts, "variable");
        await setListOrder(t, "name"); await sleep(300);
        const firstRow = document.querySelector(".pl-members .pl-member .pl-lab");
        const nameOk = !!firstRow && firstRow.textContent.trim().startsWith(sp.set.members[byName[0]].label) && sp.set.members[byName[0]].label === "tb_t01.dat";
        const timeOk = byTime.length === n && new Set(byTime).size === n && sp.set.members.every((m) => m.modified != null);
        const varOk = byVar.every((k, i) => i === 0 || !sp.set.members[k].variable || !sp.set.members[byVar[i - 1]].variable || sp.set.members[k].variable.value >= sp.set.members[byVar[i - 1]].variable.value);
        const icons = document.querySelectorAll(".pl-members .pl-member .kind").length, boxes = document.querySelectorAll(".pl-members .pl-pick").length;
        const title = (document.querySelector(".pl-members .pl-member") || {}).title || "";
        const factsOk = /points/.test(title) && /counting time/.test(title) && /file time/.test(title);
        await setListOrder(t, "set"); sp.sel = []; if (window.__renderPanels) window.__renderPanels();
        ok = sharedOk && ownOk && ovOk && stOk && all && none && cls && rng && nameOk && timeOk && varOk && icons === n && boxes === n && factsOk;
        what = `the strip under the map on the set's extent ${ext ? `${fmt(ext.from)} to ${fmt(ext.to)} (log₁₀ counts, the core's)` : "none"} ${sharedOk ? "yes" : "no"}, then fitted to the row ${ownOk ? "yes" : "no"}; Overlay fitted to every member (top ${fmt(yAll)}) then to the 22 K member (top ${fmt(yOne)}) ${ovOk ? "yes" : "no"}; Stack likewise (${fmt(sAll)} → ${fmt(sOne)}) ${stOk ? "yes" : "no"}; chosen all ${all ? n : "no"}, none ${none ? "yes" : "no"}, by exposure class (${classes.map((c) => `${c.label}: ${c.members.length}`).join(", ")}) ${cls ? "yes" : "no"}, by 20–26 K ${rng ? "7" : sp.sel.length}; the list by name (${sp.set.members[byName[0]].label} first) ${nameOk ? "yes" : "no"}, by the file's time ${timeOk ? "yes" : "no"}, by ${sp.set.variable_name} ${varOk ? "yes" : "no"}; ${icons} icons and ${boxes} checkboxes for ${n} members, the facts on hover ${factsOk ? "yes" : "no"}; the right-hand axis ${rightAxis}`;
      }
      await note(`plots: step 0's controls (W4, I10): ${what}: ${ok ? "OK" : "FAIL"}`);
      sp.view = "overlay"; sp.linked = true; renderActive(); await sleep(200);
    }
    // step 21 (W3 lit, W6, W7): a second variable typed on each placed member with its source (the member's index as
    // minutes, "t"), kept by the set file, given back on the map's rows (core 0.79.0: `others`) and drawn as the
    // right-hand axis, its source said in the readout; the core asked through core.js's host (the app's Tauri);
    // Save's words from the host
    if (activePlots() !== t) { activateTab(t.id); await sleep(400); }
    {
      const sp = stateOf(t); const vw = views.get(t.id);
      let what = "", ok = false;
      try {
        const placed = sp.set.members.map((m, i) => (m.variable ? i : -1)).filter((i) => i >= 0);
        let s = sp.set;
        for (const i of placed) s = await invoke("plots_set_other", { id: t.setId, member: i, name: "t", value: i, note: "the self-check's clock (the member's index, in minutes)" });
        setSet(t, s); await sleep(200);
        const withOthers = s.members.filter((m) => m.others && m.others.length).length;
        sp.view = "map"; sp.cursors = true; sp.pin = null; sp.cell = null; sp.map = null; sp.mapScale = "log"; renderActive();
        await wait(() => !!(sp.map && sp.map.key.startsWith("m:"))); await sleep(400);
        const M = sp.map.v.map, img = vw.plot.state.series.find((x) => x.kind === "image");
        const rowsWith = M.rows.filter((r) => Array.isArray(r.others) && r.others.length).length;
        const drawn = !!(img && img.rightLabels && img.rightLabels.length === rowsWith && img.rightTitle === "t");
        await exportView("map right axis");   // the map with its right-hand axis, to be looked at
        const ri = M.rows.findIndex((r) => r.value === 22), cj = M.x.findIndex((x) => x >= 19.1);
        const readout = img && img.cellText ? img.cellText(cj, ri) : "";
        const said = /\bt 1[0-9]\b|\bt \d+\b/.test(readout) && /self-check's clock/.test(readout);
        const noted = /right-hand axis labels each row with its t/.test(vw.note());
        // the set file keeps them (variables = { t = { value, source } }), read back by the core
        let fileOk = false, back = "not saved (no project of one's own)";
        if (state.project && state.project.root && writable()) {
          t.saveName = "selfcheck others"; await saveSet(t); await sleep(400);
          try {
            const [bs] = await invoke("plots_open_set", { path: `${state.project.root}/${t.drv}`, root: state.project.root });
            const n2 = bs.members.filter((m) => m.others && m.others.length && m.others[0].name === "t").length;
            fileOk = n2 === placed.length; back = `${n2} of ${placed.length} members with t read back from ${t.drv}`;
            await invoke("plots_close", { id: bs.id });
          } catch (e) { back = `not read back: ${e}`; }
        } else fileOk = true;
        // W6: the core is asked through the host; W7: Save answered where the file went
        const viaHost = hasCore();
        ok = withOthers === placed.length && rowsWith === placed.length && drawn && said && noted && fileOk && viaHost;
        what = `t typed on ${withOthers} of ${placed.length} placed members (the self-check's clock); the map's rows carry others on ${rowsWith} rows; the right-hand axis ${drawn ? `drawn: ${img.rightLabels.length} labels, titled "${img.rightTitle}"` : "not drawn"}; the readout at 22 K ${said ? "says t with its source" : `does not: "${readout.slice(0, 120)}"`}; the note ${noted ? "names the axis" : "is silent"}; ${back}; the core asked through core.js's host (${IN_APP ? "the app's" : "the page's"}) ${viaHost ? "yes" : "no"}`;
      } catch (e) { what = `failed: ${e}`; }
      await note(`plots: the right-hand axis lit, the host and Save's words (W3, W6, W7; step 21): ${what}: ${ok ? "OK" : "FAIL"}`);
      sp.view = "overlay"; renderActive(); await sleep(200);
    }
    // step 22 (W9): the patterns asked of the core and not carried by the set view; the overlay of this small set
    // drawn whole (slices are for more than thirty members); the core's rebin to a step and back as measured; the y
    // interval kept by the set file (PL41)
    if (activePlots() !== t) { activateTab(t.id); await sleep(400); }
    {
      const sp = stateOf(t); const vw = views.get(t.id);
      let what = "", ok = false;
      try {
        sp.view = "overlay"; sp.linked = false; renderActive(); await wait(() => !!(sp.patterns && sp.patterns.v)); await sleep(300);
        const noXY = sp.set.members.every((m) => !Array.isArray(m.x) && !Array.isArray(m.y));
        const asked = !!(sp.patterns && sp.patterns.v && sp.patterns.v.size === sp.set.members.length);
        const curves = vw.plot.state.series.filter((x) => x.kind === "line").length;
        const first = sp.set.members[0], pf = sp.patterns.v.get(first.index);
        const drawn = curves === sp.set.members.length && !!pf && vw.plot.state.series.some((x) => x.id === first.key && x.x === pf.x);
        const whole = !vw.plot.state.drawToken || (vw.plot.state.series.length <= 30 && !(vw.plot.state.pendingDraw));
        // the core's rebin: the members read again rebinned to 0.5° (the Tb files' own step is 0.1°: fewer points), then as measured again
        const n0 = first.points;
        const r1 = await invoke("plots_rebin", { id: t.setId, step: 0.5 }); sp.rebin = 0.5; sp.patterns = null; setSet(t, r1); await sleep(300);
        const n1 = stateOf(t).set.members[0].points;
        const r2 = await invoke("plots_rebin", { id: t.setId, step: null }); sp.rebin = null; sp.patterns = null; setSet(t, r2); await sleep(300);
        const n2 = stateOf(t).set.members[0].points;
        const rebinOk = n1 < n0 && n2 === n0;
        // PL41: the y interval kept by the set file
        let yiOk = false, yiw = "not checked (no project of one's own)";
        if (state.project && state.project.root && writable()) {
          sp.yInterval = "own"; t.saveName = "selfcheck yinterval"; await saveSet(t); await sleep(300);
          const [bs, bst] = await invoke("plots_open_set", { path: `${state.project.root}/${t.drv}`, root: state.project.root }); await invoke("plots_close", { id: bs.id });
          yiOk = bst.y_interval === "own"; yiw = `y_interval ${bst.y_interval} read back from ${t.drv}`;
          sp.yInterval = "shared";
        } else yiOk = true;
        ok = noXY && asked && drawn && whole && rebinOk && yiOk;
        what = `the set view carries no x or y ${noXY ? "yes" : "no"}; the patterns asked of the core for ${sp.patterns ? sp.patterns.v.size : 0} members in ${sp.timing && sp.timing.patterns != null ? `${Math.round(sp.timing.patterns)} ms` : "?"}; the overlay's ${curves} curves drawn from them ${drawn ? "yes" : "no"}, in one frame (${sp.set.members.length} members, slices from 31) ${whole ? "yes" : "no"}; rebinned to 0.5° by the core: ${n0} → ${n1} points, as measured again ${n2} ${rebinOk ? "yes" : "no"}; ${yiw} ${yiOk ? "yes" : "no"}`;
      } catch (e) { what = `failed: ${e}`; }
      await note(`plots: a large set's patterns, slices and rebin (W9; step 22): ${what}: ${ok ? "OK" : "FAIL"}`);
      sp.rebin = null; sp.view = "overlay"; renderActive(); await sleep(200);
    }
    // steps 6 and 7: each of the five views in the figure's look, exported by Export ▸ Plot's path at 17 cm; Send to Publish
    if (activePlots() !== t) { activateTab(t.id); await sleep(400); }
    {
      const sp = stateOf(t); const vw = views.get(t.id);
      sp.prule = { kind: "step", dx: 0.22, dy: 2200, scale: "linear" }; sp.look = true; sp.figmm = 160; sp.figname = "UPV/EHU TFG, full"; sp.sbar = null;
      const ready = { stack: () => !!(sp.disp && sp.disp.v && sp.pbar && sp.pbar.v), projection: () => !!(sp.proj && sp.proj.v && sp.pbar && sp.pbar.v), map: () => !!(sp.map && sp.map.key.startsWith("m:")), difference: () => !!(sp.map && sp.map.key.startsWith("d:")), series: () => !!(sp.winres && sp.winres.v) };
      for (const k of ["stack", "projection", "map", "difference", "series"]) {
        sp.scale = k === "series" ? "log" : "linear"; sp.disp = null; sp.proj = null; sp.map = null; sp.winres = null; sp.view = k; renderActive();
        const ok = await wait(ready[k]); renderActive(); await sleep(600);
        const pst = vw ? vw.plot.state : {};
        const img = (pst.series || []).find((x) => x.kind === "image");
        const what = k === "stack" || k === "projection" ? `scale bar ${pst.scaleBar ? pst.scaleBar.label : "none"}${k === "projection" && pst.xSpan ? `, x axis ${fmt(pst.xSpan[0])} to ${fmt(pst.xSpan[1])}°` : ""}` : img ? `colour bar "${img.bar.label}"` : `y "${pst.yLabel}"`;
        const good = ok && pst.look === "figure" && (k === "stack" || k === "projection" ? !!pst.scaleBar : k === "difference" ? /signed √ scale/.test(img ? img.bar.label : "") : k === "map" ? /log scale/.test(img ? img.bar.label : "") : true);
        await note(`plots: ${k} in the figure's look: ${what}: ${good ? "OK" : "FAIL"}`);
        if (state.project && state.project.root && writable()) {
          for (const ext of k === "projection" ? ["png", "pdf"] : ["png"]) {
            const fig = inLightPalette(() => plots.figure()); if (fig && fig.after) fig.after();
            if (fig) { try { const out = await invoke("export_figure", { svg: fig.svg, path: `${state.project.root}/Exports/selfcheck plots ${k} figure.${ext}`, dpi: 300 }); await note(`plots: ${k} figure exported (17 cm): ${out}`); } catch (e) { await note(`plots: ${k} figure not exported: ${e} FAIL`); } }
          }
        }
      }
      // step 9: the stack and the projection at a journal's widths, from Publish's presets (IUCr: single and double)
      if (state.project && state.project.root && writable()) {
        const widths = await invoke("plots_figure_widths");
        const iucr = widths.filter((w) => w.preset === "iucr" && (w.class === "single" || w.class === "double"));
        for (const w of iucr) for (const k of ["stack", "projection"]) {
          sp.scale = "linear"; sp.view = k; sp.figmm = w.mm; sp.figname = `IUCr, ${w.class}`; renderActive(); await wait(ready[k]); renderActive(); await sleep(500);
          const fig = inLightPalette(() => plots.figure()); if (fig && fig.after) fig.after();
          if (fig) { try { const out = await invoke("export_figure", { svg: fig.svg, path: `${state.project.root}/Exports/selfcheck plots ${k} iucr ${w.class}.png`, dpi: 300 }); await note(`plots: ${k} at IUCr's ${w.class} column, ${fmt(w.mm)} mm: ${out}: OK`); } catch (e) { await note(`plots: ${k} at IUCr's ${w.class} column not exported: ${e} FAIL`); } }
        }
        await note(`plots: Publish's presets give ${widths.length} figure widths (${[...new Set(widths.map((w) => w.preset))].join(", ")}): ${widths.length > 5 && iucr.length === 2 ? "OK" : "FAIL"}`);
        sp.figmm = 160; sp.figname = "UPV/EHU TFG, full";
      }
      sp.scale = "linear"; sp.view = "projection"; renderActive(); await sleep(300);
      if (state.project && state.project.root && writable()) {
        t.saveName = "selfcheck publish";
        const sent = await sendToPublish(t); await sleep(600);
        let back = null;
        if (sent) { try { const [bs, bst] = await invoke("plots_open_set", { path: `${state.project.root}/${sent.drv}`, root: state.project.root }); back = bst; await invoke("plots_close", { id: bs.id }); } catch (e) { back = null; } }
        await note(`plots: Send to Publish: ${sent ? `${sent.drv} written, displacement ${back ? `${back.rule} Δx ${fmt(back.dx)} Δy ${fmt(back.dy)} on the ${back.scale} scale` : "not read back"}; handed to Publish: ${sent.publishTab ? `${sent.publishTab.path} opened (${sent.view})` : "nothing opened"}` : "not saved"}: ${sent && back && back.dx === 0.22 && back.dy === 2200 && back.scale === "linear" && sent.publishTab ? "OK" : "FAIL"}`);
        // step 8: the whole path to a paper's figure: Publish builds the projection in LaTeX (its own Build command, in its
        // tab), and Plots exports its own figure of the same view at Publish's width (16 cm, upvehu-tfg), to compare them
        if (sent && sent.publishTab) {
          const t0 = Date.now(); await runCommand("publish.build");
          let built = null;
          // up to 30 minutes: LuaLaTeX's first run in a fresh home folder builds its font cache first (12 minutes and
          // more under nice, step 14's first run; the launcher now keeps that cache between runs)
          while (Date.now() - t0 < 1800000) { try { const v = await invoke("publish_open", { path: sent.publishTab.setPath }); const it = v.items.find((i) => i.pdf); if (it) { built = it; break; } } catch (e) { /* not yet */ } await sleep(2000); }
          await note(`plots: Publish built the projection: ${built ? `${built.id}, ${built.pdf}, in ${Math.round((Date.now() - t0) / 1000)} s` : "no PDF in 30 minutes"}: ${built ? "OK" : "FAIL"}`);
          activateTab(t.id); await sleep(600);
          sp.figmm = 160; sp.view = "projection"; sp.look = true; renderActive(); await sleep(800);
          for (const ext of ["pdf", "png"]) {
            const fig = inLightPalette(() => plots.figure()); if (fig && fig.after) fig.after();
            if (fig) { try { const out = await invoke("export_figure", { svg: fig.svg, path: `${state.project.root}/Exports/selfcheck plots projection 16 cm.${ext}`, dpi: 300 }); await note(`plots: Plots' own projection figure at Publish's width, 16 cm: ${out}`); } catch (e) { await note(`plots: Plots' own projection figure not exported: ${e} FAIL`); } }
          }
        }
      }
      sp.look = false; renderActive(); await sleep(200);
    }
    const f = plots.figure(); await note(`plots: Export ▸ Plot's figure of the ${ps.view}: ${f ? Math.round(f.svg.length / 1024) + " kB of SVG" : "none"}: ${f && f.svg.length > 1000 ? "OK" : "FAIL"}`);
    // a FullProf .seq beside the buffer (the self-check's copy holds LaMnO3's la0.seq): its tab and its series
    const seqFile = state.project && state.project.files ? state.project.files.find((x) => !x.dir && x.kind === "fullprof_sequential") : null;
    if (seqFile && window.__openPath) {
      await window.__openPath(seqFile.path, "plots"); await sleep(800);
      const st = activePlots(); const sp = st ? stateOf(st) : null;
      const ok = !!(sp && sp.seq) && await wait(() => !!(sp.seqcol && sp.seqcol.v));
      await sleep(300);
      const col = ok ? sp.seqcol.v.filter(Boolean) : [];
      const cv = document.querySelector(".pl-area canvas"); const r = cv ? cv.getBoundingClientRect() : { width: 0, height: 0 };
      await note(`plots: sequence ${seqFile.name}: ${sp && sp.seq ? `${sp.seq.blocks.length} refinements, ${sp.seq.names.length} parameters; ${sp.param} in ${col.length} of them, ${col.length ? `${fmt(col[0][0])} ± ${fmt(col[0][1])} to ${fmt(col[col.length - 1][0])} ± ${fmt(col[col.length - 1][1])}` : ""}` : "not opened"}; canvas ${Math.round(r.width)}×${Math.round(r.height)} px: ${ok && col.length && r.height > 40 ? "OK" : "FAIL"}`);
      if (ok && state.project.root && writable()) {
        const fig = inLightPalette(() => plots.figure()); if (fig && fig.after) fig.after();
        if (fig) { try { const out = await invoke("export_figure", { svg: fig.svg, path: `${state.project.root}/Exports/selfcheck plots sequence.png`, dpi: 150 }); await note(`plots: sequence exported: ${out}`); } catch (e) { await note(`plots: sequence not exported: ${e}`); } }
      }
    } else await note("plots: no .seq in the project: the sequence is not checked");
    // step 9: a series of refinements, its results files in the project's Results/ (the self-check's copy holds
    // La2CoMnO6's cases D and E with their data files)
    const rfiles = state.project && state.project.files ? state.project.files.filter((x) => !x.dir && /^Results\/.*\.cif$/i.test(x.path)) : [];
    if (rfiles.length > 1) {
      const rt = await openResultsPaths(rfiles.map((x) => `${state.project.root}/${x.path}`), "La2CoMnO6 series"); state.tabs.set(rt.id, rt); openTab(rt); await sleep(800);
      const rp = stateOf(rt); const v = rp.rs;
      const out = v.points.filter((p) => !p.variable);
      await note(`plots: results series: ${v.points.length} refinements, ${v.names.length} parameters; ${v.points.map((p) => `${p.label.slice(0, 14)}… T ${p.variable ? `${fmt(p.variable.value)} K from ${p.variable.source}` : `none (${p.no_data || "its data file says none the core reads"})`}`).join("; ")}; ${rp.param} ${v.points.map((p) => { const y = resultsValue(p, rp.param); return y ? `${fmt(y[0])} ± ${fmt(y[1])}` : "—"; }).join(", ")}: ${v.points.length === rfiles.length && v.points.some((p) => p.variable) ? "OK" : "FAIL"}`);
      // case E's temperature is read by the core from t7.dat's WinPLOTR header since core 0.71.0 (PL32: 179.99 °C,
      // 453.14 K, Josu's word L243); a refinement still without one is typed as Josu gave it
      const tE = v.points.find((p) => p.data_file && p.data_file.endsWith("t7.dat"));
      await note(`plots: case E's temperature from its data file (PL32): ${tE && tE.variable ? `${fmt(tE.variable.value)} K from ${tE.variable.source}` : "none"}: ${tE && tE.variable && Math.abs(tE.variable.value - 453.14) < 1e-6 && tE.variable.source === "the file's header" ? "OK" : "FAIL"}`);
      for (const p of out) rp.rs = await invoke("plots_results_set_variable", { id: rt.rsId, k: v.points.indexOf(p), value: 453.14, note: "Josu (L243): t7.dat's TEMPERATURE 179.99 is °C" });
      rp.look = true; rp.figmm = 88.5; renderActive(); await sleep(600);
      const rv = views.get(rt.id); const n = rv ? (rv.plot.state.series[0] || { x: [] }).x.length : 0;
      await note(`plots: results series with the typed temperature: ${n} points of ${rp.param} against T: ${n === v.points.length ? "OK" : "FAIL"}`);
      if (state.project.root && writable()) {
        const fig = inLightPalette(() => plots.figure()); if (fig && fig.after) fig.after();
        if (fig) { try { const o = await invoke("export_figure", { svg: fig.svg, path: `${state.project.root}/Exports/selfcheck plots results ${rp.param}.png`, dpi: 300 }); await note(`plots: results series exported: ${o}`); } catch (e) { await note(`plots: results series not exported: ${e} FAIL`); } }
      }
      // step 10: the series to Publish's F5 (PL34): sent when Publish takes a series, else said
      const sent = await sendSeriesToPublish(rt, "mode-amplitudes");
      if (sent && sent.sent) {
        const t0 = Date.now(); await runCommand("publish.build");
        let built = null;
        while (Date.now() - t0 < 600000) { try { const pv = await invoke("publish_open", { path: sent.publishTab.setPath }); const it = pv.items.find((i) => i.pdf); if (it) { built = it; break; } } catch (e) { /* not yet */ } await sleep(2000); }
        await note(`plots: Send Series to Publish: F5 of ${sent.series.length} refinements built: ${built ? `${built.pdf}: OK` : "no PDF in 10 minutes: FAIL"}`);
      } else await note(`plots: Send Series to Publish: ${sent ? `Publish takes no series yet (PL34); the series it would get: ${sent.series.map((x) => `${x.results.split("/").pop().slice(0, 14)}… ${fmt(x.temperature)} K from ${x.source.slice(0, 40)}`).join("; ")}: not checked` : "nothing to send: FAIL"}`);
    } else await note("plots: no results files in the project's Results/: the results series is not checked");
    await note("done");
  },
  async init() {
    bus.on("theme", () => renderActive());
    bus.on("settings", () => renderActive());
    bus.on("tab", () => renderBottom());
    window.__drsPlotsState = () => { const t = activePlots(); return t ? { view: stateOf(t).view, sel: stateOf(t).sel, members: members(stateOf(t)).length } : null; };
    // the module's styles, until they move into the shell's style sheet
    const st = document.createElement("style");
    st.textContent = `.pl-areas{flex:1;min-height:0;display:flex;gap:6px}.pl-area,.pl-zoom{position:relative;flex:1;min-height:0;min-width:0}.pl-zoom{flex:0.7}.pl-zoom canvas{width:100%;height:100%}.pl-link.on,.pl-cur.on,.pl-follow.on,.pl-yint.on{background:var(--sel);border-color:var(--accent)}.pl-pick{margin:0 2px 0 0;flex:none}.pl-member .kind{flex:none;color:var(--muted)}.pl-memctl{padding:2px 12px 4px;display:flex;flex-wrap:wrap;gap:4px;align-items:center}.pl-memctl select{font-size:11px}.pl-mapcol{flex:1;min-width:0;min-height:0;display:flex;flex-direction:column;gap:6px}.pl-rowstrip{position:relative;flex:0 0 30%;min-height:0}.pl-colstrip{position:relative;flex:0 0 28%;min-width:0}.pl-rowstrip canvas,.pl-colstrip canvas{width:100%;height:100%}.pl-areas{container-type:size}@container (max-height:260px){.pl-rowstrip{display:none}}@container (max-width:520px){.pl-colstrip{display:none}}.pl-ops{font-size:10px;color:var(--accent);border:1px solid var(--accent);border-radius:4px;padding:0 3px;margin-left:4px}.pl-area canvas{width:100%;height:100%}.pl-note{flex:none;padding:4px 12px 8px;font-size:11px;line-height:1.4}.pl-member{display:flex;align-items:center;gap:6px;padding:2px 12px;cursor:default}.pl-member.sel{background:var(--sel)}.pl-sw{width:12px;height:3px;flex:none}.pl-lab{flex:1;min-width:0;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}.pl-var{flex:none;width:52px;text-align:right}.pl-op{width:56px;flex:none}.pl-bottom{height:100%;display:flex;flex-direction:column;overflow:hidden}.pl-bottom .dtable{height:100%;display:flex;flex-direction:column}.pl-bottom .dt td:first-child,.pl-bottom .dt th:first-child{position:sticky;left:0;background:var(--panel);text-align:left;z-index:1}`;
    document.head.appendChild(st);
  },
};
