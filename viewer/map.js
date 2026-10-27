/* diff-ref-struct · the Map viewer (Plots' step 18, I08 step 5, R33; L409): the Map view's part of the Plots module,
   factored out so that two hosts draw the same map and nothing is drawn twice: the app (the Plots module, through
   its shared plot and the core's commands) and the page the core writes for a machine without the app (`drs plots
   show`, PL38), where the data is one JSON and the core cannot be asked. The viewer draws the core's numbers from
   its data and computes none of them (L139, PL14): the image series with its clip and marks, the row and column
   cursors' strips, the zoom in the T–2θ plane kept across redraws, "range from the view" (the core's block of the
   rectangle, through the host), the regions and the jump to one. It imports only the four shell files it leans on:
   plot.js, mapclip.js, ticks.js, svgexport.js (through plot.js). Everything it needs beyond its data it asks of the
   host (below); every change of its state is told to the host, which draws again (`show`): the viewer never redraws
   itself, so the app's note and panels stay in step with it.

   The data, one object (the core's PL38 `page_data` writes it; the app host builds it from `plots_map` and the set):
     title, path                 the set's title and the file's path (words; the page shows them)
     core, written               the core's version and the date the page was written (null in the app)
     variable: { name, unit }    "T", "K"
     x: { label, quantity }      the label the map is drawn in ("2θ (°)", "d (Å)") and its word ("2theta", "d", "q",
                                 "stl", "s"; null for TOF, energy or mixed members)
     y_word                      "counts" or "normalised counts"
     members: [{ key, label, variable: { value, source } | null, pattern: { x, y } | null, kind, kind_name,
                 points, counting_time, time, facts }]
                                 by their index in the set, as the map's rows name them; the pattern as the set holds
                                 it (the row cursor's strip; the map's own row when absent, as in the page); step 20
                                 (W4): `key` the host's name for the member (its index as text on the page), `kind`
                                 the file's format id and `kind_name` its name (the core's describe; the icon), `points`,
                                 `counting_time` (the exposure class), `time` (the file's modification time, seconds),
                                 `facts` (the words shown on hover)
     y_extent: { from, to } | null
                                 the set's y extent on the map's scale (the core's), the shared y interval of the row
                                 strip (step 20, W4); null: the strip fits the row shown
     maps: { map, difference }   each { map: { x, x_edges, rows: [{ member, value, edges: [lo, hi], others: [{ name,
                                 value, source }] }], z: [[…]], range: [lo, hi] }, apart: [{ member, z: […] | null }],
                                 scale: "linear" | "sqrt" | "log", reference: index | null } as `plots_map` returns it,
                                 or null when not made; a row's `others` (W3, core 0.79.0: the member's other variables
                                 with their sources) labels the right-hand axis with the first of them and the readout
                                 says them all; none, no right axis
     clip: { map, difference }   each { from, to } | null, the set file's (PL37), on that map's scale
     regions: [{ name, x: [a, b], x_quantity, variable: [a, b], map_range: { from, to } | null,
                 difference_range: { from, to } | null, window,
                 blocks: { map: block | null, difference: block | null }, refused: string | null }]
                                 a region's own colour range per view (P65, step 19): `map_range` on the map's scale,
                                 `difference_range` on the difference map's; the one of the view shown is written when
                                 the region is saved and applied when it is chosen
                                 a block as `plots::map_region` gives it: { rows: [i…], columns: [j…], x: [lo, hi],
                                 range: [lo, hi] }; precomputed by the core for the page, asked of the core in the app
     regions_refused             the core's word on the file's regions at opening, or null
     same_grid                   whether the members share one grid

   The host: { mapRegion(region, { key, scale, reference }) → Promise<block> (the core's `map_region`; the page answers
   with the saved blocks or throws, saying so), rangeWarning(range, { key, scale, reference }) → Promise<string | null>
   (the core's `range_warning`; the page answers null), status(text, bad), onChange() (the state changed: draw again),
   ask → Promise<string | null>, given a question and a default value (a name typed; the page's prompt) }.

   The state, the host's object (the Plots tab's own state in the app), read and written here: clip { map, difference }
   (each { from, to, scale, block?, warning? } | null), clipBefore, cursors, pin, cell, mapZoom, rangeFromView,
   regions, selregion, regionsRefused; step 20 (W4): yInterval ("shared" | "own"), sel (the keys of the members
   chosen), sortBy ("set" | "variable" | "time" | "name", the members' list order). */
import { createPlot } from "../plot.js";
import { keptZoom, followSwitch, stripFrames } from "../mapclip.js";

/** A number for the eye (as the shell's `fmt`): 7 significant digits at most, exponent form far from one. */
export function fmtNumber(v) {
  if (v == null || !isFinite(v)) return "";
  if (v === 0) return "0";
  const a = v < 0 ? -v : v;
  if (a >= 1e6 || a < 1e-4) return v.toExponential(4).replace(/\.?0+e/, "e");
  return String(parseFloat(v.toPrecision(7)));
}
/** The clip of a view in force: kept only while the map is drawn on the scale it was set on. */
export function clipOf(state, key, scale) { const c = state.clip && state.clip[key]; return c && c.scale === scale ? c : null; }
/** The scale the difference map is drawn on for the map's colour scale: a difference has no log, √ keeps its sign. */
export const diffScaleOf = (mapScale) => (mapScale === "log" ? "sqrt" : mapScale);
/** The word for the value's scale: "log₁₀ counts", "√counts", "counts". */
export const zWordOf = (scale, yw) => (scale === "log" ? `log₁₀ ${yw}` : scale === "sqrt" ? `√${yw}` : yw);
/** The field of a region that holds its colour range on a view (P65): `map_range` on the map, `difference_range` on
    the difference map. */
export const rangeFieldOf = (which) => (which === "difference" ? "difference_range" : "map_range");
/** A region's own colour range for the view shown (P65, step 19): the one saved on that view, never the other's (a
    range on the map's log₁₀ scale means nothing on the difference's signed √); null when it has none there. */
export function regionRangeFor(region, which) {
  const r = region ? region[rangeFieldOf(which)] : null;
  return r && r.from != null && r.to != null ? { from: r.from, to: r.to } : null;
}

// ---------------------------------------------------------------- the members' list (step 20, W4), pure parts
/** The family of a file's kind, from the core's format id and name (the icon's letter and its words): columns of text
    (xy, xye, xys, csv), WinPLOTR and FullProf, Bruker, PANalytical, GSAS, Rigaku, DBWS, HDF/NeXus, zarr, CIF; "file"
    when none of them. */
export function kindFamily(kind, kindName = "") {
  const k = `${kind || ""} ${kindName || ""}`.toLowerCase();
  const of = (re, letter, words) => (re.test(k) ? { letter, words } : null);
  return of(/winplotr|fullprof|\bbuf\b|\bdat\b|\bprf\b/, "F", "WinPLOTR / FullProf") || of(/bruker|\braw\b|\buxd\b/, "B", "Bruker") || of(/panalytical|xrdml|philips/, "P", "PANalytical") || of(/gsas|fxye|\bgsa\b/, "G", "GSAS") || of(/rigaku|\bras\b|\brasx\b/, "R", "Rigaku") || of(/dbws/, "D", "DBWS") || of(/nexus|hdf|h5|nxazint/, "H", "HDF5 / NeXus") || of(/zarr|midas/, "Z", "zarr") || of(/\bcif\b/, "C", "CIF") || of(/\bxye?\b|\bxys\b|csv|tsv|column|text|table|\bxy\b/, "≡", "columns of text") || { letter: "·", words: "a data file" };
}
/** The icon of a file's kind: a small square with the family's letter (SVG text, the markup's own). */
export function kindIcon(kind, kindName = "") {
  const f = kindFamily(kind, kindName);
  const t = f.letter.replace(/[<>&]/g, "");
  return `<svg class="kind" viewBox="0 0 16 16" width="14" height="14" aria-label="${f.words}" role="img"><rect x="1.5" y="1.5" width="13" height="13" rx="3" fill="none" stroke="currentColor" stroke-width="1.1"/><text x="8" y="11.6" text-anchor="middle" font-size="9.5" font-family="-apple-system, sans-serif" fill="currentColor">${t}</text></svg>`;
}
/** The members' list order: the set's own, by the variable (those without one last), by the file's time (those
    without one last), or by name; a stable order, no number changed. */
export function sortMembers(members, by) {
  const idx = members.map((m, i) => i);
  if (by === "set" || !by) return idx;
  const key = by === "variable" ? (m) => (m.variable ? m.variable.value : null) : by === "time" ? (m) => (m.time == null ? null : m.time) : (m) => String(m.label || "");
  return idx.sort((a, b) => {
    const ka = key(members[a]), kb = key(members[b]);
    if (ka == null && kb == null) return a - b;
    if (ka == null) return 1;
    if (kb == null) return -1;
    if (typeof ka === "string") return ka.localeCompare(kb, undefined, { numeric: true }) || a - b;
    return ka < kb ? -1 : ka > kb ? 1 : a - b;
  });
}
/** The exposure classes of the members: the distinct counting times the files state (W2's class by the esd's
    median is the core's, later), each with its members' indices; "not stated" for those without one. */
export function exposureClasses(members) {
  const out = [];
  members.forEach((m, i) => {
    const t = m.counting_time == null ? null : m.counting_time;
    let c = out.find((x) => x.time === t);
    if (!c) { c = { time: t, label: t == null ? "not stated" : `${fmtNumber(t)} s`, members: [] }; out.push(c); }
    c.members.push(i);
  });
  return out.sort((a, b) => (a.time == null ? 1 : b.time == null ? -1 : a.time - b.time));
}
/** The members whose variable lies within [from, to] (ends included); those without one never. */
export function membersInRange(members, from, to) {
  return members.map((m, i) => (m.variable && m.variable.value >= from && m.variable.value <= to ? i : -1)).filter((i) => i >= 0);
}
/** A map's row's other variables (W3, core 0.79.0): `others`, a list of { name, value, source } as `plots_map` and
    `page_data` write it; a `variables` object ({ name: { value, source } } or { name: value }) is read too. Each as
    { name, value, source }. */
export function rowVariables(row) {
  if (!row) return [];
  if (Array.isArray(row.others)) return row.others.filter((o) => o && typeof o.value === "number" && o.name).map((o) => ({ name: String(o.name), value: o.value, source: o.source || null }));
  const v = row.variables;
  if (!v || typeof v !== "object") return [];
  return Object.keys(v).map((name) => { const x = v[name]; if (typeof x === "number") return { name, value: x, source: null }; if (x && typeof x === "object" && typeof x.value === "number") return { name, value: x.value, source: x.source || null }; return null; }).filter(Boolean);
}
/** The first other variable of a map's row (W3), the right-hand axis's label: { name, value, source } or null. */
export function rowVariable(row) {
  const all = rowVariables(row);
  return all.length ? all[0] : null;
}

/** els: { plot, plotOpts } (the host's plot and its options object: the app's shared plot), or { canvas, readout } (a
    plot made here, the page's); rowEl, rowCanvas, rowReadout, colEl, colCanvas, colReadout (the cursors' strips);
    note (an element the viewer writes its words into); clipText (the clip in short, optional); regions (an element
    for the regions' list, optional: the app lists them in its sidebar). */
export function createMapViewer(els, host, state) {
  const popts = els.plotOpts || { xLabel: "2θ (°)", yLabel: "T (K)", nearestInY: true, right: "auto" };
  const plot = els.plot || createPlot(els.canvas, els.readout, null, popts);
  let data = null, which = null, look = false, active = false;
  let rplot = null, cplot = null, followTimer = null;
  const current = () => (data && which ? data.maps[which] : null);
  const diff = () => which === "difference";
  const stripsOn = () => active && !!state.cursors && !look;
  const ctx = () => { const m = current(); return m ? { key: which, scale: m.scale, reference: diff() ? m.reference : null } : null; };
  const xWord = () => (data && data.x ? data.x.quantity : null);
  const xName = () => (data && data.x && data.x.label ? data.x.label.replace(/ \([^)]*\)$/, "") : "x");

  // ---------------------------------------------------------------- the cursors' strips (step 16)
  function ensureStrips() {
    if (!rplot) {
      // step 23 (L431): the strips take the map's frame edges at each draw (`stripFrames`), so the row strip shares
      // the map's vertical lines and the column strip, drawn with the variable up (a path: the values across), its rows
      rplot = createPlot(els.rowCanvas, els.rowReadout, null, { xLabel: "2θ (°)", yLabel: "counts", right: 16 });
      cplot = createPlot(els.colCanvas, els.colReadout, null, { xLabel: "counts", yLabel: "T (K)", right: 16, yFromZero: false });
    }
    return [rplot, cplot];
  }
  /** The cursors' cell in words: its x and its member. */
  function cellWords(c) {
    const m = current(); if (!m || !c) return "";
    const M = m.map;
    const mIdx = typeof c.r === "string" ? (m.apart[parseInt(c.r.slice(6))] || {}).member : M.rows[c.r] ? M.rows[c.r].member : null;
    const mem = mIdx != null ? data.members[mIdx] : null;
    return `${xName()} ${M.x[c.j] == null ? "?" : M.x[c.j].toFixed(3)}${mem ? `, ${mem.label}` : ""}`;
  }
  function redrawStrips() {
    if (!stripsOn()) { if (rplot) { rplot.setSeries([]); cplot.setSeries([]); } return; }
    const [rp, cp] = ensureStrips();
    const m = current(), c = state.pin || state.cell, xl = data.x.label, yw = data.y_word, vn = data.variable.name, vu = data.variable.unit;
    const frames = stripFrames(plot.state.geom);
    if (!m || !c) { rp.setAxes({ xLabel: xl, yLabel: `${yw}: the row under the pointer`, frame: frames.row }); cp.setAxes({ xLabel: "the column under the pointer", yLabel: `${vn} (${vu})`, frame: frames.col }); rp.setSeries([]); cp.setSeries([]); return; }
    const M = m.map, scale = m.scale, zWord = zWordOf(scale, yw);
    const apartK = typeof c.r === "string" ? parseInt(c.r.slice(6)) : -1;
    const row = apartK >= 0 ? null : M.rows[c.r];
    const mIdx = row ? row.member : m.apart[apartK] ? m.apart[apartK].member : null;
    const mem = mIdx == null ? null : data.members[mIdx];
    if (!mem || M.x[c.j] == null) { rp.setSeries([]); cp.setSeries([]); return; }
    const T = mem.variable ? `${vn} = ${fmtNumber(mem.variable.value)} ${vu}` : `no ${vn}`;
    const colour = memberColour(mIdx);
    // the row cursor: the Map draws the member's pattern as the set holds it (when the data carries it; the page draws
    // the map's own row); the Difference Map the core's row of it (the member less the reference), what that map shows
    // the y label says the unit once: the plot itself adds the scale's word to its y label, so none is added here (L431)
    if (diff() || !mem.pattern) {
      const z = apartK >= 0 ? m.apart[apartK].z || null : M.z[c.r];
      const zs = z ? z.map((v) => (v == null ? NaN : v)) : null;
      rp.setAxes({ xLabel: xl, yLabel: `${diff() ? "less the reference, " : ""}${yw}`, yScale: scale, frame: frames.row });
      rp.setSeries(zs ? [{ id: "row", kind: "line", label: `${mem.label}${diff() ? " less the reference" : ""}, ${T}`, x: M.x, y: zs, scaled: true, values: zs, color: colour, width: 1 }] : [], true);
    } else {
      rp.setAxes({ xLabel: xl, yLabel: yw, yScale: scale, frame: frames.row });
      rp.setSeries([{ id: `m${mIdx}`, kind: "line", label: `${mem.label}, ${T}`, x: mem.pattern.x, y: mem.pattern.y, color: colour, width: 1 }], true);
    }
    // the strip's x is the map's (its frame's left and right the map's, step 23); its y the set's extent (shared, the
    // core's: the map's range for the difference's rows, `y_extent` on the map's scale for a member's pattern) or the
    // row's own (step 20, W4)
    const v = plot.getView();
    const shared = state.yInterval === "shared" ? (diff() || !mem.pattern ? { from: M.range[0], to: M.range[1] } : data.y_extent || null) : null;
    rp.setView({ x0: v.x0, x1: v.x1, y0: shared ? shared.from : null, y1: shared ? shared.to : null });
    // the column cursor: the matrix's column at the cell's x, the variable up (a path, the values across on the map's
    // scale), the frame's top and bottom and the y view the map's (step 23), so that each row sits on its line
    const held = M.rows.map((r, i) => i).filter((i) => M.z[i][c.j] != null);
    const xs = held.map((i) => M.rows[i].value), ys = held.map((i) => M.z[i][c.j]);
    cp.setAxes({ xLabel: `${diff() ? "less the reference, " : ""}${zWord}`, yLabel: `${vn} (${vu})`, yScale: "linear", frame: frames.col });
    cp.setSeries([{ id: "column", kind: "both", path: true, points: true, radius: 2.5, label: `${xName()} ${M.x[c.j].toFixed(3)}`, x: ys, y: xs, color: "--ink", width: 1 }]);
    const g = plot.state.geom;
    if (g) { cp.fitAll(); const cv = cp.getView(); cp.setView({ x0: cv.x0, x1: cv.x1, y0: g.y0, y1: g.y1 }); }
  }
  // step 23: the map's frame moves with its canvas (its bottom with the height, its right with the bar), so the strips
  // take it again after each resize of the map; the plot's own observer, registered first, has redrawn the map by then
  if (typeof ResizeObserver !== "undefined" && plot.canvas) new ResizeObserver(() => { if (stripsOn() && rplot) redrawStrips(); }).observe(plot.canvas);
  /** A member's colour by its place among the members with a value (magma's middle part, as the module's), grey apart. */
  const MAGMA = ["#2c115f", "#51127c", "#721f81", "#932b80", "#b73779", "#d8456c", "#f1605d", "#fc8961"];
  function memberColour(k) {
    const mem = data.members[k]; if (!mem || !mem.variable) return "--muted";
    const placed = data.members.filter((x) => x.variable); const i = placed.indexOf(mem);
    const t = placed.length > 1 ? i / (placed.length - 1) : 0;
    const q = t * (MAGMA.length - 1), a = Math.min(MAGMA.length - 2, Math.floor(q)), f = q - a;
    const ca = MAGMA[a].match(/\w\w/g).map((h) => parseInt(h, 16)), cb = MAGMA[a + 1].match(/\w\w/g).map((h) => parseInt(h, 16));
    return `rgb(${ca.map((v, i2) => Math.round(v + (cb[i2] - v) * f)).join(",")})`;
  }
  popts.onCell = (c) => {
    if (!stripsOn() || !c) return;   // the pointer off the map: the cursors stay where it left them
    const same = state.cell && state.cell.j === c.j && state.cell.r === c.r;
    state.cell = { j: c.j, r: c.r };
    if (state.pin || same) return;
    plot.setCellCursor(state.cell); redrawStrips();
  };
  popts.onCellClick = (c) => {
    if (!stripsOn() || !c) return;
    if (state.pin) { state.pin = null; state.cell = { j: c.j, r: c.r }; plot.setCellCursor(state.cell); redrawStrips(); host.status("The cursors released: they follow the pointer again."); }
    else { state.pin = { j: c.j, r: c.r }; plot.setCellCursor(state.pin); redrawStrips(); host.status(`The cursors pinned at ${cellWords(state.pin)}; a click on the map releases them.`); }
  };

  // ---------------------------------------------------------------- the zoom kept and range from the view (step 17)
  plot.onViewChange((v) => {
    if (rplot && stripsOn()) redrawStrips();   // step 23: the row strip's x view and the column strip's y view follow the map's (its frame is set by then)
    if (!active) return;
    const z = keptZoom(v, plot.state.fit);
    const same = (z == null && state.mapZoom == null) || (z && state.mapZoom && z.x0 === state.mapZoom.x0 && z.x1 === state.mapZoom.x1 && z.y0 === state.mapZoom.y0 && z.y1 === state.mapZoom.y1);
    state.mapZoom = z;
    if (!same && state.rangeFromView) { clearTimeout(followTimer); followTimer = setTimeout(() => followRange(), 150); }
  });
  /** The rectangle shown as a region for the core: the zoom, or the whole map; x in the quantity the map is shown in. */
  function regionOfView(name = "the view") {
    const c = ctx(); if (!c) return null;
    const z = state.mapZoom, g = plot.state.geom, v = plot.getView();
    const x = z ? [z.x0, z.x1] : [v.x0, v.x1];
    const y = z && z.y0 != null ? [z.y0, z.y1] : g ? [g.y0, g.y1] : null;
    if (!y) return null;
    return { name, x, x_quantity: xWord(), variable: y, map_range: null, window: null };
  }
  /** "Range from the view": the clip set to the core's range of the rectangle shown (the host's block); nothing read
      off the matrix here. */
  async function followRange() {
    const c = ctx(), r = regionOfView(); if (!c || !r || !state.rangeFromView) return;
    try {
      const b = await host.mapRegion(r, c);
      if (!state.rangeFromView) return;
      state.clip[c.key] = { from: b.range[0], to: b.range[1], scale: c.scale, block: b, warning: null };
      host.onChange();
    } catch (e) { host.status(`Range from the view: ${e.message || e}`, true); }
  }
  function toggleFollow() {
    const c = ctx(); if (!c) { host.status("Show the map or the difference map first: the range follows its zoom.", true); return; }
    state.rangeFromView = !state.rangeFromView;
    const s = followSwitch(state.clip[c.key], state.clipBefore[c.key], state.rangeFromView);
    state.clip[c.key] = s.clip; state.clipBefore[c.key] = s.before;
    if (state.rangeFromView) { host.status("The colour range follows the zoom: the core's range of the rectangle shown (plots::map_region)."); followRange(); }
    else host.status(`The colour range as it was${state.clip[c.key] ? ` (${fmtNumber(state.clip[c.key].from)} to ${fmtNumber(state.clip[c.key].to)})` : " (the core's)"}.`);
    host.onChange();
  }

  // ---------------------------------------------------------------- the clip (step 16)
  /** The clip of a view set (a drag of the bar's marks, `quiet`, or typed): kept with the scale it was set on. */
  function setClip(key, from, to, scale, quiet = false) {
    state.clip[key] = quiet ? { from, to, scale, warning: null } : { from, to, scale };
    host.onChange();
    if (!quiet) host.status(`${key === "difference" ? "The difference map" : "The map"} shows ${fmtNumber(from)} to ${fmtNumber(to)} on its ${scale} colour scale; values beyond are saturated.`);
  }
  function resetClip(key) {
    if (!state.clip[key]) return;
    state.clip[key] = null; host.onChange();
    host.status(`${key === "difference" ? "The difference map" : "The map"} on the core's colour range again.`);
  }
  /** A pair typed: inside the core's range, in order; a field left empty keeps the core's end, both empty the core's range. */
  function typedClip(key, from, to) {
    const m = data.maps[key]; if (!m) return false;
    const M = m.map, what = key === "difference" ? "difference map" : "map";
    if (from == null && to == null) { if (state.clip[key]) resetClip(key); else host.status("The core's range, as before."); return true; }
    const a = from == null ? M.range[0] : from, b = to == null ? M.range[1] : to;
    if (!isFinite(a) || !isFinite(b)) { host.status("From and to are numbers.", true); return false; }
    if (!(b > a)) { host.status("To lies above from.", true); return false; }
    if (b <= M.range[0] || a >= M.range[1]) { host.status(`The pair ${fmtNumber(a)} to ${fmtNumber(b)} lies wholly outside the core's range ${M.range[0].toFixed(3)} to ${M.range[1].toFixed(3)}: nothing of the ${what} would be seen.`, true); return false; }
    setClip(key, a, b, m.scale);
    return true;
  }

  // ---------------------------------------------------------------- the regions (step 17)
  /** A region added after the core accepts it on the map shown; `clip` (when given) is written into the field of the
      view shown (P65): `map_range` on the map, `difference_range` on the difference map. */
  async function addRegion(name, rect, clip = null) {
    name = String(name || "").trim();
    if (!name || state.regions.some((r) => r.name === name)) { host.status(`A region needs a name of its own${name ? `: "${name}" is taken` : ""}.`, true); return null; }
    const c = ctx(); if (!c) { host.status("Show the map or the difference map first: a region is the rectangle it shows.", true); return null; }
    const region = { name, x: [rect.x[0], rect.x[1]], x_quantity: rect.x_quantity === undefined ? xWord() : rect.x_quantity, variable: [rect.variable[0], rect.variable[1]], map_range: null, difference_range: null, window: rect.window || null };
    if (clip) region[rangeFieldOf(c.key)] = { from: clip.from, to: clip.to };
    let b;
    try { b = await host.mapRegion(region, c); } catch (e) { host.status(`Region "${name}" not saved: ${e.message || e}`, true); return null; }
    state.regions = [...state.regions, region]; state.selregion = state.regions.length - 1;
    const own = regionRangeFor(region, c.key);
    host.status(`Region "${name}": ${xName()} ${fmtNumber(region.x[0])} to ${fmtNumber(region.x[1])}, ${data.variable.name} ${fmtNumber(region.variable[0])} to ${fmtNumber(region.variable[1])} ${data.variable.unit}; the core's block: rows ${b.rows[0] + 1}–${b.rows[b.rows.length - 1] + 1}, columns ${b.columns[0] + 1}–${b.columns[b.columns.length - 1] + 1}, its cells ${fmtNumber(b.range[0])} to ${fmtNumber(b.range[1])}${own ? `; its colour range on the ${c.key === "difference" ? "difference map" : "map"} ${fmtNumber(own.from)} to ${fmtNumber(own.to)}` : ""}. Save (⌘S) writes it to the set file.`);
    host.onChange();
    return region;
  }
  /** Save Region…: the rectangle shown named (the host asks the name), with the clip when Range from the View is on. */
  async function saveRegion() {
    const c = ctx(), r = regionOfView();
    if (!c || !r) { host.status("Show the map or the difference map first: a region is the rectangle it shows.", true); return null; }
    let k = state.regions.length + 1; while (state.regions.some((x) => x.name === `r${k}`)) k++;
    const clip = state.rangeFromView ? clipOf(state, c.key, c.scale) : null;
    const name = await host.ask(`Save the rectangle shown, ${xName()} ${fmtNumber(r.x[0])} to ${fmtNumber(r.x[1])} and ${data.variable.name} ${fmtNumber(r.variable[0])} to ${fmtNumber(r.variable[1])} ${data.variable.unit}, as a region by name${r.x_quantity ? ` (in ${r.x_quantity})` : ""}${clip ? `, with its colour range ${fmtNumber(clip.from)} to ${fmtNumber(clip.to)} (Range from the View is on)` : ""}. Choosing it zooms to it again${clip ? " and applies its colour range" : ""}.`, `r${k}`);
    if (name == null) return null;
    return addRegion(name, r, clip);
  }
  /** A region chosen: the map zoomed to the core's block of it (its bounds in the map's quantity), its own colour
      range on the view shown applied (P65: the map's or the difference map's, never the other's; the block's range
      when it has none there and Range from the View is on); refused in the core's words. */
  async function chooseRegion(i) {
    const r = state.regions[i]; if (!r) return false;
    const c = ctx();
    if (!c) { state.selregion = i; host.status(`Region "${r.name}": the map is being computed; choose it again.`, true); host.onChange(); return false; }
    try {
      const b = await host.mapRegion(r, c);
      state.selregion = i;
      state.mapZoom = { x0: b.x[0], x1: b.x[1], y0: r.variable[0], y1: r.variable[1] };
      const own = regionRangeFor(r, c.key);
      if (own) state.clip[c.key] = { from: own.from, to: own.to, scale: c.scale };
      else if (state.rangeFromView) state.clip[c.key] = { from: b.range[0], to: b.range[1], scale: c.scale, block: b, warning: null };
      const other = regionRangeFor(r, c.key === "difference" ? "map" : "difference");
      host.status(`Region "${r.name}": rows ${b.rows[0] + 1}–${b.rows[b.rows.length - 1] + 1}, columns ${b.columns[0] + 1}–${b.columns[b.columns.length - 1] + 1} of the map (${xName()} ${fmtNumber(b.x[0])} to ${fmtNumber(b.x[1])}), its cells ${fmtNumber(b.range[0])} to ${fmtNumber(b.range[1])}${own ? `; its colour range on the ${c.key === "difference" ? "difference map" : "map"} ${fmtNumber(own.from)} to ${fmtNumber(own.to)} applied` : other ? `; its colour range is the ${c.key === "difference" ? "map's" : "difference map's"}, not applied here` : ""}.`);
      host.onChange();
      return true;
    } catch (e) { host.status(`Region "${r.name}": ${e.message || e}`, true); return false; }
  }
  async function renameRegion(i, name = null) {
    if (i == null) i = state.selregion; const r = state.regions[i]; if (!r) return false;
    if (name == null) { name = await host.ask(`The region ${xName()} ${fmtNumber(r.x[0])} to ${fmtNumber(r.x[1])}, ${data.variable.name} ${fmtNumber(r.variable[0])} to ${fmtNumber(r.variable[1])}: its name.`, r.name); if (name == null) return false; }
    name = String(name).trim();
    if (!name || state.regions.some((x, k) => k !== i && x.name === name)) { host.status(`A region needs a name of its own${name ? `: "${name}" is taken` : ""}.`, true); return false; }
    state.regions = state.regions.map((x, k) => (k === i ? { ...x, name } : x));
    host.onChange(); return true;
  }
  function removeRegion(i = null) {
    if (i == null) i = state.selregion; const r = state.regions[i]; if (!r) return;
    state.regions = state.regions.filter((_, k) => k !== i); state.selregion = null;
    host.status(`Region "${r.name}" removed; Save (⌘S) writes the set file without it.`); host.onChange();
  }
  /** The regions' list (the page's; the app lists them in its sidebar): a click zooms to one, × removes it. */
  function renderRegions() {
    const box = els.regions; if (!box) return;
    box.innerHTML = "";
    if (!state.regions.length) { box.textContent = "No region saved: Save region… names the rectangle shown."; return; }
    state.regions.forEach((r, k) => {
      const row = document.createElement("div"); row.className = "region" + (k === state.selregion ? " sel" : "");
      row.innerHTML = `<span class="name"></span><span class="sub"></span><button class="small x" title="Remove region" aria-label="Remove region">×</button>`;
      row.querySelector(".name").textContent = r.name;
      const ranges = [r.map_range ? `map ${fmtNumber(r.map_range.from)} to ${fmtNumber(r.map_range.to)}` : "", r.difference_range ? `difference ${fmtNumber(r.difference_range.from)} to ${fmtNumber(r.difference_range.to)}` : ""].filter(Boolean);
      row.querySelector(".sub").textContent = `${fmtNumber(r.x[0])}–${fmtNumber(r.x[1])}${r.x_quantity && r.x_quantity !== xWord() ? ` ${r.x_quantity}` : ""}, ${fmtNumber(r.variable[0])}–${fmtNumber(r.variable[1])}${ranges.length ? ` · ${ranges.join("; ")}` : ""}${r.refused ? " · refused" : ""}`;
      row.title = r.refused ? `${r.name}: ${r.refused}` : `${r.name}: click to zoom the map to it (its colour range on the view shown applied), double-click to rename it.`;
      row.onclick = (e) => { if (e.target.classList.contains("x")) { removeRegion(k); return; } chooseRegion(k); };
      row.ondblclick = () => renameRegion(k);
      box.appendChild(row);
    });
  }

  // ---------------------------------------------------------------- the drawing
  /** The map drawn: `which` "map" or "difference"; `opts.look` the figure's look (no cursors, no marks). The note is
      written; the zoom and the clip in force put back. */
  function show(d, w, opts = {}) {
    data = d; which = w; look = !!opts.look; active = true;
    popts.axisWheel = true;
    const m = current();
    const strips = stripsOn();
    if (els.rowEl) els.rowEl.hidden = !strips; if (els.colEl) els.colEl.hidden = !strips;
    if (!strips) plot.setCellCursor(null);
    if (!m) { plot.setSeries([]); if (els.note) els.note.textContent = diff() ? "The difference map subtracts a reference member from every other: choose one and Set Reference." : "The map is being computed…"; renderRegions(); return; }
    const M = m.map, yw = data.y_word, vn = data.variable.name, vu = data.variable.unit, zWord = zWordOf(m.scale, yw);
    const refM = diff() && m.reference != null ? data.members[m.reference] : null;
    const clip = clipOf(state, which, m.scale);
    if (els.clipText) els.clipText.textContent = clip ? `${fmtNumber(clip.from)} to ${fmtNumber(clip.to)}` : "the core's range";
    // the core's word on a typed or restored clip (range_warning), asked once per clip; a dragged pair is not asked
    if (clip && clip.warning === undefined) {
      clip.warning = null;
      Promise.resolve(host.rangeWarning({ from: clip.from, to: clip.to }, ctx())).then((wd) => { if (wd) { clip.warning = wd; host.onChange(); } }).catch(() => {});
    }
    plot.setAxes({ xLabel: data.x.label, yLabel: `${vn} (${vu})`, yScale: "linear", look: look ? "figure" : null, xSpan: null, scaleBar: null, keepView: true });
    plot.setSeries([{
      id: "map", kind: "image", xEdges: M.x_edges,
      rows: M.rows.map((r, i) => ({ lo: r.edges[0], hi: r.edges[1], z: M.z[i] })),
      range: M.range, ramp: diff() ? "diverging" : "sequential",
      bar: { label: `${diff() ? `${refM ? `less ${refM.variable ? `${fmtNumber(refM.variable.value)} ${vu} (${refM.label})` : refM.label}` : "less the reference"}, ${yw}` : yw}${m.scale === "log" ? ", log scale" : m.scale === "sqrt" ? (diff() ? ", signed √ scale" : ", √ scale") : ""}`, scale: m.scale,
        clip: { from: clip ? clip.from : null, to: clip ? clip.to : null, onChange: (from, to) => setClip(which, from, to, m.scale, true), onReset: () => resetClip(which) } },
      apart: m.apart.map((a) => ({ label: data.members[a.member] ? data.members[a.member].label : `member ${a.member + 1}`, z: a.z, note: a.z ? "" : `no ${vn}: its row in the difference map comes with the core's PL18` })),
      // the right-hand axis (W3, step 20): each row's second variable as the core gives it; none, no axis
      ...rightAxisOf(M),
      cellText: (j, r) => {
        const x = M.x[j];
        if (typeof r === "string") { const a = m.apart[parseInt(r.slice(6))]; const z = a.z ? a.z[j] : null; return `${xName()} ${x.toFixed(4)}  ${data.members[a.member].label} (no ${vn})  ${z == null ? "" : `${zWord} ${z.toFixed(3)}`}`; }
        const row = M.rows[r], z = M.z[r][j];
        // the row's other variables with their sources (W3), said in the readout as the pointer hovers the row
        const others = rowVariables(row).map((o) => `${o.name} ${fmtNumber(o.value)}${o.source ? ` (${o.source})` : ""}`).join(", ");
        return `${xName()} ${x.toFixed(4)}  ${vn} ${fmtNumber(row.value)} ${vu} (${data.members[row.member].label})${others ? `  ${others}` : ""}  ${z == null ? "no value" : `${diff() ? "less the reference, " : ""}${zWord} ${z.toFixed(3)}`}`;
      },
    }], true);
    // the zoom (step 17): the rectangle kept put back after the redraw, the whole map when none
    if (state.mapZoom) { const zv = plot.getView(); if (zv.x0 !== state.mapZoom.x0 || zv.x1 !== state.mapZoom.x1 || zv.y0 !== state.mapZoom.y0 || zv.y1 !== state.mapZoom.y1) plot.setView(state.mapZoom); }
    else if (!plot.state.fit) plot.fitAll();
    // the cursors (step 16): the cell marked on the map and the two strips drawn for it
    plot.setCellCursor(strips ? state.pin || state.cell : null); redrawStrips();
    if (els.note) els.note.textContent = noteText(m, clip, strips);
    renderRegions(); renderMembers();
  }
  /** The right-hand axis's labels from the rows' second variables (W3): one per row at its centre, the axis titled
      with the variable's name; nothing when no row carries one (core 0.79.0 gives them). */
  function rightAxisOf(M) {
    const labels = [], names = new Set();
    M.rows.forEach((r) => { const rv = rowVariable(r); if (rv) { labels.push({ at: r.value, text: fmtNumber(rv.value) }); names.add(rv.name); } });
    return labels.length ? { rightLabels: labels, rightTitle: [...names].join(", ") } : {};
  }
  /** The members' list (the page's; the app lists them in its sidebar) with the same parts: an icon by the file's
      kind, a checkbox per member (the keys chosen in `state.sel`), the facts on hover, the order chosen. */
  function renderMembers() {
    const box = els.members; if (!box || !data) return;
    const vn = data.variable.name, vu = data.variable.unit;
    const sel = new Set(state.sel || []);
    box.innerHTML = "";
    const order = sortMembers(data.members, state.sortBy || "set");
    for (const k of order) {
      const mem = data.members[k];
      const row = document.createElement("div"); row.className = "region member" + (sel.has(mem.key) ? " sel" : "");
      row.innerHTML = `<input type="checkbox" class="pick" aria-label="${mem.label.replace(/"/g, "&quot;")} chosen">${kindIcon(mem.kind, mem.kind_name)}<span class="name"></span><span class="sub"></span>`;
      row.querySelector(".pick").checked = sel.has(mem.key);
      row.querySelector(".name").textContent = mem.label;
      row.querySelector(".sub").textContent = mem.variable ? `${fmtNumber(mem.variable.value)} ${vu}` : `no ${vn}`;
      row.title = mem.facts || (mem.variable ? `${vn} = ${mem.variable.value} ${vu}, from ${mem.variable.source}` : `no ${vn}: drawn apart`);
      row.querySelector(".pick").onchange = (e) => { pick(mem.key, e.target.checked); };
      box.appendChild(row);
    }
  }
  /** A member chosen or not (its key in `state.sel`); the host told. */
  function pick(key, on) {
    const sel = (state.sel || []).filter((k) => k !== key);
    state.sel = on ? [...sel, key] : sel;
    host.onChange();
  }
  /** All, none, an exposure class (by its index in `exposureClasses`), or a range of the variable chosen. */
  function pickAll(on) { state.sel = on ? data.members.map((m) => m.key) : []; host.onChange(); }
  function pickClass(i) { const c = exposureClasses(data.members)[i]; if (!c) return; state.sel = c.members.map((k) => data.members[k].key); host.onChange(); }
  function pickRange(from, to) { state.sel = membersInRange(data.members, from, to).map((k) => data.members[k].key); host.onChange(); }
  function setSort(by) { state.sortBy = by; host.onChange(); }
  function toggleYInterval() {
    state.yInterval = state.yInterval === "own" ? "shared" : "own";
    host.status(state.yInterval === "shared" ? "The y interval shared by the whole set: the row strip on the set's extent (the core's)." : "The y interval adapted to the row shown.");
    host.onChange();
  }
  function noteText(m, clip, strips) {
    const M = m.map, yw = data.y_word, vn = data.variable.name, vu = data.variable.unit, zWord = zWordOf(m.scale, yw);
    const refM = diff() && m.reference != null ? data.members[m.reference] : null;
    const clipText = clip
      ? ` Colour range shown: ${fmtNumber(clip.from)} to ${fmtNumber(clip.to)} ${zWord}, ${state.rangeFromView && clip.block ? "from the view (the core's block of the rectangle shown, plots::map_region)" : "set on the colour bar (Range… types it)"}: values below ${fmtNumber(clip.from)} are drawn in the first colour and values above ${fmtNumber(clip.to)} in the last, saturated; the core's range stays ${M.range[0].toFixed(3)} to ${M.range[1].toFixed(3)}, and a double-click on the bar returns to it. Saved with the set file (PL37).${clip.warning ? ` The core says: ${clip.warning}.` : ""}`
      : ` Drag the two marks on the colour bar (or Range…) to clip the range shown: values beyond the marks are drawn saturated at the end colours, so that the weak details come up as the strong peaks are clipped.`;
    const curText = strips ? ` Cursors: the row under the pointer (its member's pattern, below, its y interval ${state.yInterval === "shared" ? "shared by the whole set" : "its own"}) and the column (its values against ${vn}, beside)${state.pin ? `, pinned at ${cellWords(state.pin)}: a click on the map releases them` : "; a click on the map pins them"}.` : "";
    const rightText = M.rows.some((r) => rowVariable(r)) ? ` The right-hand axis labels each row with its ${[...new Set(M.rows.map((r) => (rowVariable(r) || {}).name).filter(Boolean))].join(", ")} (the member's second variable, the core's).` : "";
    const z = state.mapZoom;
    const zoomText = z ? ` Zoomed to ${xName()} ${fmtNumber(z.x0)} to ${fmtNumber(z.x1)}, ${vn} ${fmtNumber(z.y0)} to ${fmtNumber(z.y1)} ${vu} (a rectangle dragged zooms both axes, the wheel over an axis zooms that axis, a double-click shows the whole map; Save Region… names it).` : ` Drag a rectangle to zoom both axes, the wheel over an axis to zoom that axis.`;
    const regText = `${state.regions.length ? ` Regions saved: ${state.regions.map((r) => r.name).join(", ")} (choosing one zooms to it and applies its colour range).` : ""}${state.regionsRefused ? ` The core refuses a region of the file on this set: ${state.regionsRefused}.` : ""}`;
    const apart = m.apart.map((a) => (data.members[a.member] || {}).label).filter(Boolean);
    const apartText = apart.length ? ` Drawn apart, no ${vn} in the file: ${apart.join(", ")}.` : "";
    return `${diff() ? `Every member less ${refM ? refM.label : "the reference"} (the core's difference map, a range symmetric about zero), red above and blue below. ` : ""}One row per member at its ${vn}, from the midpoint to its neighbours (plots::map); colour range from the core: ${M.range[0].toFixed(3)} to ${M.range[1].toFixed(3)} ${zWord}.${clipText}${zoomText}${regText}${curText}${rightText} Grid: ${data.same_grid ? "the members' own points (they share them)" : "the common grid (the overlap, at the finest step)"}.${apartText}`;
  }
  /** Another view is shown: the cursors' marks and strips go, the zoom stays in the state for the next show. */
  function hide() { active = false; popts.axisWheel = false; plot.setCellCursor(null); if (rplot) { rplot.setSeries([]); cplot.setSeries([]); } if (els.rowEl) els.rowEl.hidden = true; if (els.colEl) els.colEl.hidden = true; }
  function toggleCursors() {
    state.cursors = !state.cursors; state.pin = null;
    host.status(state.cursors ? "Cursors on: the row under the pointer below the map, the column beside it; a click pins them." : "Cursors off.");
    host.onChange();
  }

  return {
    plot, show, hide, setClip, resetClip, typedClip, toggleCursors, toggleFollow, followRange, regionOfView, addRegion, saveRegion, chooseRegion, renameRegion, removeRegion,
    pick, pickAll, pickClass, pickRange, setSort, toggleYInterval,
    cell: (c) => popts.onCell(c), cellClick: (c) => popts.onCellClick(c), strips: () => ({ row: rplot, col: cplot }),
    current, which: () => which, data: () => data,
  };
}
