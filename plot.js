/* diff-ref-struct · the shared plot component (APP_CONVENTIONS §8, C8.1): one component for every module. Series
   (points, lines, marker rows, bands), axes with a label and a scale (linear, √, log on y), zoom with a box, the
   wheel or a pinch, pan with two fingers (or a drag with Space / the middle button), double-click to fit, a
   crosshair readout, a selection of points by a box (in select mode), Copy as an image. Plain canvas, no library;
   it draws what it is given and computes nothing but screen coordinates. The Data module was its first user;
   Refine's Rietveld plot is on it since app 0.8.0 (reflection ticks as marker rows at a height, the difference offset). */
import { axisTicks, sqrtTicks, logTicks, symlogTicks, symlog, symlogInv, fmtTick } from "./ticks.js";
import { svgRecorder, legendItems, drawLegend, svgText } from "./svgexport.js";
import { clipPlace, barValue, clipMoved, cellAt, columnAt, blockOf } from "./mapclip.js";   // PL13l, PL13m: a map's clip, the cell under the pointer, the cells in view; pure

// the theme's colours, read once a drawing (PL13: hundreds of series asked the style sheet each time)
let cssCache = null;
const css = (name) => {
  if (cssCache && cssCache.has(name)) return cssCache.get(name);
  const v = getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  if (cssCache) cssCache.set(name, v);
  return v;
};
/** A size in points as the canvas's pixels (1/96 in, the SVG's unit): 11 pt → 14.7 px (PL30). */
const ptPx = (pt) => (pt * 96) / 72;
const colourOf = (c, fallback) => (c && c.startsWith && c.startsWith("--") ? css(c) || fallback : c || fallback);
const PALETTE = ["--phase-0", "--phase-1", "--phase-2", "--phase-3", "--phase-4", "--diff", "--obs"];
// PL13 (Plots, L169): the image layer's colour ramps, a value's place in the core's range turned into a colour (display
// only): viridis (Matplotlib's, sampled) for intensities, RdBu (ColorBrewer's, reversed: red positive) for differences.
// Exported for the 3D surface (surface3d.js, PL13k), so that a map and its surface share their colours.
const RAMPS = {
  sequential: ["#440154", "#482475", "#414487", "#355f8d", "#2a788e", "#21918c", "#22a884", "#44bf70", "#7ad151", "#bddf26", "#fde725"],
  diverging: ["#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f"],
};
export function rampColour(name, t) {
  const stops = RAMPS[name] || RAMPS.sequential;
  const k = Math.max(0, Math.min(1, t)) * (stops.length - 1);
  const i = Math.min(stops.length - 2, Math.floor(k)), f = k - i;
  const a = stops[i].match(/\w\w/g).map((h) => parseInt(h, 16)), b = stops[i + 1].match(/\w\w/g).map((h) => parseInt(h, 16));
  return `rgb(${a.map((v, q) => Math.round(v + (b[q] - v) * f)).join(",")})`;
}

/** opts: { xLabel, yLabel, yScale: "linear" | "sqrt" | "log", onSelect(indices by series id), onHover(x, y),
    onPick(x, y) (a click in "pick" mode: a background point added, W4), onHandleMove(id, i, x, y) (a handle dragged), onHandleDrag(id, i, x, y) (while it is dragged; E10l),
    onHandleSelect(id, indices), onBandCreate(a, b) (a drag in "region" mode, W7), onBandChange(id, i, a, b) (an edge
    of an editable band dragged), onBandSelect(id, i), onMarker(id, i) (a tick clicked, W8), onView(x0, x1) }
    Series kinds: "both" (points and line), "points", "line", "markers" (ticks: `row` in the tick band under the
    frame, or `at` a height in data units), "vlines", "bands" (`ranges`, `editable`), "handles" (draggable points,
    drawn as squares). A series with `strip: true` is drawn in its own strip under the tick band, with its own
    scale over the visible points (the difference curve, as WinPLOTR lays it out), never in the frame's scale.
    PL13 (Plots, L169), additive: a series may give `opacity` (0–1, dimming), `scaled: true` (its y are already on the
    plot's y scale: the core's displaced curves, never transformed again; `values`, the y the readout shows), `baselineLabel: { text, at }` (a label in the
    right margin at the height `at`, on the plot's scale; `opts.right` widens the margin). Kind "image": a map drawn
    as cells, `xEdges`, `rows: [{ lo, hi, z: [] }]` (z on the scale the core chose, null for no value), `range: [lo,
    hi]` (the core's colour range), `ramp` ("sequential" | "diverging"), `bar: { label, scale: "linear" | "sqrt" |
    "log" }` (the colour bar at the right, its ticks labelled as the axis rule labels a y axis on that scale), `apart:
    [{ label, z | null, note }]` (rows drawn in a band under the frame, the members without a value of the variable),
    and `cellText(j, row)` (the readout over a cell; row is an index of `rows`, or "apart:k"). `pointLabels: [{ x, y, text, priority }]`:
    labels at a series' own points (y on its scale), thinned with the baseline labels. `legend: false`: a series kept
    out of the legend (it keeps its label for the readout). `opts.yFromZero: false`: a linear y axis spans the data,
    not down to zero (a cell parameter against the temperature). Step 3: `open: [bool]` per point (drawn hollow: below
    2σ by the core's flag, F5), `setCursor(x)` (a dashed cursor set from another plot) and `opts.onPointer(x | null)`
    (the pointer's x, null on leaving), for cursors linked between plots (F6). `opts.nearestInY`: the readout
    names the curve nearest the pointer in x and y. E10c (Explore), additive: a series with `path: true` is a path in
    the plane (x in any order: a pair's path, an ellipse), drawn whole and clipped to the frame, with `dash` and, for
    its points, `radius`. PL13f (Plots, step 6), additive: `setAxes({ look: "figure", xSpan: [a, b], scaleBar: { value,
    label } })`: a figure's look, open axes without a grid or a frame, the x axis drawn only over [a, b] (the front
    member's range), and a scale bar of `value` y units (the core's) in place of the y axis; a point label may give
    `tick: true` (a short tick at its point, the label beside it) or `above: true` (the label over its point: an
    axis's title). Without them a plot draws as before. PL13g (Plots, step 7), additive, in a figure's look only: the
    text in ink rather than grey, points drawn larger, a log or √ y axis said as ", log scale" or ", √ scale", the scale
    bar's label running up from its foot, a map kept in its frame; `toSVG({ w, h })` draws the figure at a given size
    (a printed width) on white. A series may give `marker` ("circle", "square", "diamond", "triangle"): its points in
    that shape, hollow where `open`, in every look.
    App 0.38.0, for Explore (E8, E8a) and Plots (PL13, PL30), additive: a path series may give `fill: true` (closed and
    filled, `fillOpacity`: the core's ellipse outlines as patches) and `numbered: true` or `numbers: []` (its points
    numbered); an image may give `contours: [{ x: [], y: [], label, dash }]` (lines the caller gives, drawn over it,
    never computed here); `setAxes({ yScale: "symlog" })` (linear within ±1, logarithmic beyond); a bands series may
    give `hatch: true`; a line may give `fillBelow: true` (filled down to the frame in `fillColor`, the background by
    default: a stack's hidden parts); `opts.font = { label_pt, tick_pt }` and `toSVG({ w, h, label_pt, tick_pt })` (text
    in points) with the margins then measured from what is drawn, as in a figure's look; `opts.right = "auto"`;
    `onViewChange(fn)` and `linkAxes(plots, { x, y })` (axes linked between plots). PL35 (Plots, step 12): the measured
    right margin counts the point labels that reach past the frame. E10m (Explore), additive:
    `opts.xIntegers` (read at each drawing): an x axis of whole numbers, labelled at whole numbers only. A line of more points than pixels is
    drawn as each pixel column's span (what the pixels can show), so that hundreds of series stay quick.
    PL13l (Plots, step 16; I08, R33), additive: an image's `bar.clip: { from, to, onChange(from, to), onReset() }`: the
    range shown is [from, to] (null: the core's end), on the scale the values are on; a value beyond it takes the end
    colour; two marks on the bar at the bounds are dragged (`onChange` on every move, the pair kept in order and inside
    the core's range), and a double-click on the bar calls `onReset`; the bar's axis stays the core's range, its colours
    showing the saturation. `opts.onCell({ id, j, r } | null)` on the pointer's move over an image's cell (`r` an index
    of `rows`, or "apart:k") and `opts.onCellClick(cell)` on a click on one; `setCellCursor({ j, r } | null)` marks a
    column and a row of the image. Without them a map draws as before.
    PL13m (Plots, step 17; I08 step 4), additive: `setAxes({ …, keepView: true })` keeps the y zoom (a redraw of a
    zoomed map for its clip or its cursors); `opts.axisWheel`: the wheel over the y axis's margin zooms y about the
    pointer, over the x axis's band zooms x, inside the frame as before; an image draws only the cells in view (the
    column span by `blockOf`).
    PL13p (Plots, step 20; W4), additive: an image's `rightLabels: [{ at, text }]` (a label at the height `at`, data
    units, in the right margin before the colour bar: a row's second variable, W3) with `rightTitle`, the labels
    thinned as the baseline labels are and the margin measured for them; `opts.fitTo: [ids]` (read at each fit): the
    y range fitted to those series alone, the others drawn as they fall (a y interval adapted to the member chosen). */
/** Settings ▸ Plots ▸ Marker size (§4): the scale of the points the plots draw on screen, which the app keeps in the
    page's `--marker-scale`; each plot draws again when it changes (the "drs-markers" event). */
export function markerScale() {
  const v = parseFloat(getComputedStyle(document.documentElement).getPropertyValue("--marker-scale"));
  return isFinite(v) && v > 0 ? v : 1;
}

export function createPlot(canvas, readout, legendEl, opts = {}) {
  let ctx = canvas.getContext("2d");
  const viewListeners = [];   // PL13: told of every change of the view (linked axes)
  let exporting = null;   // { w, h } while the plot draws into the SVG recorder (toSVG)
  const st = { series: [], x0: 0, x1: 1, y0: null, y1: null, drag: null, mode: "zoom", sel: new Map(), geom: null, xLabel: opts.xLabel || "x", yLabel: opts.yLabel || "y", yScale: opts.yScale || "linear", fit: true, look: null, xSpan: null, scaleBar: null };
  const T = { sqrt: (y) => (y >= 0 ? Math.sqrt(y) : -Math.sqrt(-y)), log: (y) => (y > 0 ? Math.log10(y) : NaN), linear: (y) => y, symlog };
  const Tinv = { sqrt: (v) => (v >= 0 ? v * v : -v * v), log: (v) => Math.pow(10, v), linear: (v) => v, symlog: symlogInv };
  // PL30: the text's sizes in pixels (11 px unless given in points), and whether the margins are measured
  let fontPt = opts.font || null;
  const fontOf = (px) => `${px}px -apple-system, sans-serif`;
  const ty = (y) => T[st.yScale](y);

  function layout() {
    if (exporting) return exporting;
    const dpr = window.devicePixelRatio || 1;
    const w = canvas.clientWidth, h = canvas.clientHeight;
    canvas.width = Math.round(w * dpr); canvas.height = Math.round(h * dpr);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    return { w, h };
  }
  /** A series' y on the plot's scale: as given when the series is `scaled` (PL13), else transformed. */
  const fy = (s) => (s.scaled ? (v) => (v == null ? NaN : v) : ty);
  function xRange() {
    let lo = Infinity, hi = -Infinity;
    for (const s of st.series) {
      if (s.hidden || s.kind === "bands") continue;   // bands (excluded regions, windows) carry `ranges`, no x: they never set the range
      if (s.kind === "image") { const e = s.xEdges; if (e && e.length) { lo = Math.min(lo, e[0], e[e.length - 1]); hi = Math.max(hi, e[0], e[e.length - 1]); } continue; }
      if (s.x.length) { lo = Math.min(lo, s.x[0], s.x[s.x.length - 1]); hi = Math.max(hi, s.x[0], s.x[s.x.length - 1]); }
    }
    return isFinite(lo) ? [lo, hi] : [0, 1];
  }
  function yRange(x0, x1) {
    let lo = Infinity, hi = -Infinity, ilo = Infinity, ihi = -Infinity;
    // PL13p: fitted to the series named in `opts.fitTo` when they are drawn, else to every series
    const fitTo = opts.fitTo && opts.fitTo.length && st.series.some((s) => !s.hidden && opts.fitTo.includes(s.id)) ? opts.fitTo : null;
    for (const s of st.series) {
      if (s.hidden || (fitTo && !fitTo.includes(s.id))) continue;
      if (s.kind === "image") { for (const r of s.rows) { ilo = Math.min(ilo, r.lo); ihi = Math.max(ihi, r.hi); } continue; }
      if (s.kind === "markers" && s.at != null) { const v = ty(s.at); if (isFinite(v)) { lo = Math.min(lo, v); hi = Math.max(hi, v); } continue; }
      if (!s.y || s.strip) continue;
      const f = fy(s);
      for (let i = 0; i < s.x.length; i++) { if (s.x[i] < x0 || s.x[i] > x1) continue; const v = f(s.y[i]); if (!isFinite(v)) continue; if (v < lo) lo = v; if (v > hi) hi = v; }
    }
    if (isFinite(ilo) && !isFinite(lo)) return [ilo, ihi];   // an image alone: its rows' edges, no padding
    if (!isFinite(lo)) return [0, 1];
    if (st.yScale === "linear" && opts.yFromZero !== false) lo = Math.min(lo, 0);   // PL13: off for quantities whose change is the point
    const pad = (hi - lo) * 0.05 || 1;
    return [lo - (st.yScale === "linear" ? 0 : pad), hi + pad];
  }

  function draw() {
    cssCache = new Map();
    try { drawNow(); } finally { cssCache = null; }
  }
  function drawNow() {
    const { w, h } = layout();
    ctx.clearRect(0, 0, w, h);
    if (!st.series.length || w < 40 || h < 40) { st.geom = null; return; }
    const ink = css("--ink"), muted = st.look === "figure" ? ink : css("--muted"), line = css("--line");   // PL13g: a figure's text in ink
    // PL30: the text in points for a figure (`opts.font`, `toSVG`'s sizes), else 11 px; the margins measured from what
    // is drawn there when the text is in points or the look is a figure's, else the screen's fixed ones
    const tickPx = fontPt ? ptPx(fontPt.tick_pt || 9) : 11, labelPx = fontPt ? ptPx(fontPt.label_pt || 10) : 11;
    const measured = !!fontPt || st.look === "figure";
    const legendBottom = legendEl ? legendEl.offsetTop + legendEl.offsetHeight : 0;
    const top = measured ? Math.max(Math.round(tickPx / 2) + 4, legendBottom + 6) : Math.max(30, legendBottom + 10);
    const axisH = measured ? Math.round(tickPx + 8 + labelPx + 6) : 34;
    // the tick band (one row per marker series with `row`) and the strip (the difference) sit under the frame; an
    // image's rows apart (PL13) sit right under the frame, above the tick band
    const rowsN = st.series.filter((s) => !s.hidden && s.kind === "markers" && s.row != null).reduce((m, s) => Math.max(m, s.row + 1), 0);
    const bandH = rowsN ? 8 + 10 * rowsN : 0;
    const strips = st.series.filter((s) => !s.hidden && s.strip && s.y);
    const stripH = strips.length ? Math.max(50, Math.round((h - axisH - top) * 0.18)) : 0;
    const images = st.series.filter((s) => !s.hidden && s.kind === "image");
    const apartN = images.reduce((m, s) => Math.max(m, (s.apart || []).length), 0);
    const apartH = apartN ? 8 + 18 * apartN : 0;
    const bottom = h - axisH - stripH - bandH - apartH;
    const apartTop = bottom + 8, bandTop = bottom + apartH, stripTop = bandTop + bandH, stripBottom = stripTop + stripH;
    const [y0, y1] = st.y0 == null ? yRange(st.x0, st.x1) : [st.y0, st.y1];
    // the y ticks first: their widest label sets a measured left margin (PL30 b)
    const yt = st.yScale === "linear" ? axisTicks(y0, y1, bottom - top, null) : st.yScale === "sqrt" ? sqrtTicks(y0, y1, bottom - top) : st.yScale === "symlog" ? symlogTicks(y0, y1) : logTicks(y0, y1, bottom - top);
    ctx.font = fontOf(tickPx);
    let left = 64;
    if (measured) {
      const widest = yt.ticks.reduce((m, t) => Math.max(m, ctx.measureText(t.label).width), 0);
      left = st.scaleBar && st.yScale === "linear" ? Math.round(labelPx + 26) : Math.round(widest + 10 + labelPx + 8);
    }
    // PL30 d: a right margin measured from what is drawn there (the baseline and point labels, a colour bar)
    let right = opts.right && opts.right !== "auto" ? opts.right : 16;
    if (opts.right === "auto" || (measured && !opts.right)) {
      let r = 12;
      for (const s of st.series) {
        if (s.hidden) continue;
        if (s.baselineLabel) r = Math.max(r, ctx.measureText(s.baselineLabel.text).width + 12);
        // PL35 (Plots, step 12): a point label beyond the frame (the projection's temperatures at its oblique axis's
        // far end): at the fraction f of the x range its point lies, its right edge stays inside the canvas when
        // right ≥ w − left − (w − 4 − left − offset − width) / f; screen layout only
        if (s.pointLabels && st.x1 > st.x0) {
          for (const p of s.pointLabels) {
            const f = (p.x - st.x0) / (st.x1 - st.x0);
            if (!(f > 0) || f > 1.02) continue;   // a point outside the view draws no label
            const d = p.above ? -4 : p.tick ? 9 : 6;
            r = Math.max(r, w - left - (w - 4 - left - d - ctx.measureText(p.text).width) / f);
          }
          r = Math.min(r, w * 0.45);
        }
        if (s.kind === "image" && s.bar) {
          const [zl, zh] = s.range || [0, 1];
          const bt = s.bar.scale === "log" ? logTicks(zl, zh, bottom - top) : s.bar.scale === "sqrt" ? sqrtTicks(zl, zh, bottom - top) : axisTicks(zl, zh, bottom - top, null);
          r = Math.max(r, rightLabelsWidth(s) + 14 + 12 + 5 + bt.ticks.reduce((m, t) => Math.max(m, ctx.measureText(t.label).width), 0) + 6 + labelPx + 4);
        }
      }
      right = Math.round(r);
    }
    const sx = (x) => left + (x - st.x0) / (st.x1 - st.x0) * (w - left - right);
    const sy = (v) => bottom - (v - y0) / (y1 - y0) * (bottom - top);
    // the strip's own scale over the visible points, symmetric about zero
    let sa = 0;
    for (const s of strips) for (let i = 0; i < s.x.length; i++) { if (s.x[i] < st.x0 || s.x[i] > st.x1) continue; const v = Math.abs(s.y[i]); if (isFinite(v) && v > sa) sa = v; }
    sa = sa || 1;
    const sys = (v) => stripTop + stripH / 2 - v / sa * (stripH / 2 - 3);
    st.geom = { left, right, top, bottom, w, h, sx, sy, y0, y1, bandTop, stripTop, stripBottom, sys, axisTop: stripBottom, apartTop };
    // grid and axes
    ctx.font = fontOf(tickPx); ctx.fillStyle = muted; ctx.strokeStyle = line; ctx.lineWidth = 1;
    const xt = axisTicks(st.x0, st.x1, w - left - right, (l) => ctx.measureText(l).width);
    // E10m (Explore): an axis of whole numbers (cycles, steps) labelled at whole numbers only
    if (opts.xIntegers && xt.step < 1) xt.ticks = xt.ticks.filter((t) => Math.abs(t.value - Math.round(t.value)) < 1e-9).map((t) => ({ value: Math.round(t.value), label: String(Math.round(t.value)) }));
    const fig = st.look === "figure";
    // PL13f: a figure's x axis, a line with ticks only over its span (the front member's range), no grid
    const [xa, xb] = fig && st.xSpan ? [Math.max(st.x0, st.xSpan[0]), Math.min(st.x1, st.xSpan[1])] : [st.x0, st.x1];
    if (fig) { ctx.strokeStyle = ink; ctx.beginPath(); ctx.moveTo(Math.round(sx(xa)), Math.round(stripBottom) + 0.5); ctx.lineTo(Math.round(sx(xb)), Math.round(stripBottom) + 0.5); ctx.stroke(); }
    ctx.textAlign = "center";
    for (const t of xt.ticks) {
      if (fig && (t.value < xa - 1e-9 || t.value > xb + 1e-9)) continue;
      const X = Math.round(sx(t.value)) + 0.5;
      if (fig) { ctx.strokeStyle = ink; ctx.beginPath(); ctx.moveTo(X, stripBottom); ctx.lineTo(X, stripBottom + 4); ctx.stroke(); ctx.strokeStyle = line; }
      else { ctx.globalAlpha = 0.35; ctx.beginPath(); ctx.moveTo(X, top); ctx.lineTo(X, bottom); ctx.stroke(); ctx.globalAlpha = 1; }
      ctx.fillText(t.label, X, stripBottom + 5 + tickPx);
    }
    ctx.font = fontOf(labelPx);
    ctx.fillText(withFactor(st.xLabel, xt.suffix, fig), fig ? (sx(xa) + sx(xb)) / 2 : left + (w - left - right) / 2, h - 4);
    ctx.font = fontOf(tickPx);
    if (stripH) {
      // the strip: a zero line, its scale at the left
      ctx.strokeStyle = line; ctx.globalAlpha = 0.5; ctx.beginPath(); ctx.moveTo(left, Math.round(sys(0)) + 0.5); ctx.lineTo(w - right, Math.round(sys(0)) + 0.5); ctx.stroke(); ctx.globalAlpha = 1;
      ctx.textAlign = "right"; ctx.fillStyle = muted; ctx.fillText(fmtTick(sa, sa), left - 6, sys(sa) + 4); ctx.fillText(fmtTick(-sa, sa), left - 6, sys(-sa) + 4);
      ctx.strokeStyle = ink; ctx.globalAlpha = 0.4; ctx.strokeRect(left + 0.5, stripTop + 0.5, w - left - right - 1, stripH - 1); ctx.globalAlpha = 1;
      ctx.textAlign = "center";
    }
    ctx.textAlign = "right";
    // the y ticks (above): on the transformed scale, labelled with the true values (linear: the §8 rule; √, log and
    // symlog: the value)
    if (fig && st.scaleBar && st.yScale === "linear") {
      // PL13f: a scale bar of the core's value in place of the y axis, from the bottom of the frame up, its label beside it
      const Y0 = Math.round(bottom) + 0.5, Y1 = Math.round(sy(y0 + st.scaleBar.value)) + 0.5, X = left - 14.5;
      ctx.strokeStyle = ink; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.moveTo(X, Y0); ctx.lineTo(X, Y1); ctx.stroke(); ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(X - 4, Y0); ctx.lineTo(X + 4, Y0); ctx.moveTo(X - 4, Y1); ctx.lineTo(X + 4, Y1); ctx.stroke();
      ctx.save(); ctx.translate(X - 8, Y0); ctx.rotate(-Math.PI / 2); ctx.textAlign = "left"; ctx.fillStyle = ink; ctx.fillText(st.scaleBar.label, 0, 0); ctx.restore();   // PL13g: from the bar's foot up
      ctx.fillStyle = muted;
    } else {
      if (fig) { ctx.strokeStyle = ink; ctx.beginPath(); ctx.moveTo(left + 0.5, top); ctx.lineTo(left + 0.5, bottom); ctx.stroke(); ctx.strokeStyle = line; }
      for (const t of yt.ticks) {
        const Y = Math.round(sy(t.value)) + 0.5;
        if (fig) { ctx.strokeStyle = ink; ctx.beginPath(); ctx.moveTo(left - 4, Y); ctx.lineTo(left, Y); ctx.stroke(); ctx.strokeStyle = line; }
        else { ctx.globalAlpha = 0.25; ctx.beginPath(); ctx.moveTo(left, Y); ctx.lineTo(w - right, Y); ctx.stroke(); ctx.globalAlpha = 1; }
        ctx.fillText(t.label, left - 6, Y + tickPx * 0.36);
      }
      const said = fig ? (st.yScale === "sqrt" ? ", √ scale" : st.yScale === "log" ? ", log scale" : st.yScale === "symlog" ? ", symmetric log scale" : "") : (st.yScale === "sqrt" ? " (√)" : st.yScale === "log" ? " (log)" : st.yScale === "symlog" ? " (symlog)" : "");   // PL13g
      ctx.font = fontOf(labelPx);
      ctx.save(); ctx.translate(measured ? Math.round(labelPx * 0.8) + 2 : 12, top + (bottom - top) / 2); ctx.rotate(-Math.PI / 2); ctx.textAlign = "center"; ctx.fillText(withFactor(st.yLabel + said, yt.suffix, fig), 0, 0); ctx.restore();
      ctx.font = fontOf(tickPx);
    }
    // PL13: the image layer (a map) under everything else, its colour bar in the right margin, its rows apart under the
    // frame; colours from the core's range, nothing computed but each value's place in that range
    for (const s of images) {
      const [zlo, zhi] = s.range || [0, 1];
      // PL13l: the range shown is the clip's when the bar has one (null bounds: the core's ends); a value beyond it is
      // placed at 0 or 1, the end colours
      const clip = s.bar && s.bar.clip ? s.bar.clip : null;
      const cfrom = clip && clip.from != null ? clip.from : zlo, cto = clip && clip.to != null ? clip.to : zhi;
      const place = clip ? (z) => clipPlace(z, cfrom, cto) : (z) => (zhi > zlo ? (z - zlo) / (zhi - zlo) : 0.5);
      const e = s.xEdges;
      // PL13m: the column span in view, so that a zoomed map draws its cells only
      const span = blockOf({ x0: st.x0, x1: st.x1, y0: -Infinity, y1: Infinity }, e, []) || { j0: 0, j1: -1 };
      const cells = (z, Y0, Y1) => {
        for (let j = span.j0; j <= span.j1 && j < z.length; j++) {
          if (z[j] == null || !isFinite(z[j])) continue;
          const X0 = sx(Math.max(e[j], st.x0)), X1 = sx(Math.min(e[j + 1], st.x1));
          ctx.fillStyle = rampColour(s.ramp, place(z[j])); ctx.fillRect(X0, Y0, X1 - X0 + 0.5, Y1 - Y0);
        }
      };
      ctx.save(); ctx.beginPath(); ctx.rect(left, top, w - left - right, bottom - top); ctx.clip();
      for (const r of s.rows) { const Ya = sy(r.hi), Yb = sy(r.lo); if (Yb < top || Ya > bottom) continue; cells(r.z, Math.max(top, Ya), Math.min(bottom, Yb) + 0.5); }
      // E8 c: contour lines over the map, as the caller gives them (the core's), never computed here
      for (const c of s.contours || []) {
        ctx.strokeStyle = colourOf(c.color, ink); ctx.globalAlpha = 0.75; ctx.lineWidth = c.width || 1; ctx.setLineDash(c.dash || []);
        ctx.beginPath(); let pen = false;
        for (let i = 0; i < c.x.length; i++) { const v = c.y[i]; if (!isFinite(v) || !isFinite(c.x[i])) { pen = false; continue; } const X = sx(c.x[i]), Y = sy(v); pen ? ctx.lineTo(X, Y) : ctx.moveTo(X, Y); pen = true; }
        ctx.stroke(); ctx.setLineDash([]); ctx.globalAlpha = 1;
        if (c.label && c.x.length) { const m = c.x.length >> 1; ctx.fillStyle = ink; ctx.textAlign = "center"; ctx.fillText(c.label, sx(c.x[m]), sy(c.y[m]) - 3); }
      }
      ctx.restore();
      (s.apart || []).forEach((a, i) => {
        const Y0 = apartTop + 18 * i, Y1 = Y0 + 14;
        ctx.save(); ctx.beginPath(); ctx.rect(left, Y0, w - left - right, Y1 - Y0); ctx.clip();
        if (a.z) cells(a.z, Y0, Y1);
        else { ctx.fillStyle = muted; ctx.textAlign = "center"; ctx.fillText(a.note || "", left + (w - left - right) / 2, Y0 + 11); }
        ctx.restore();
        ctx.strokeStyle = ink; ctx.globalAlpha = 0.4; ctx.strokeRect(left + 0.5, Y0 + 0.5, w - left - right - 1, Y1 - Y0 - 1); ctx.globalAlpha = 1;
        ctx.fillStyle = muted; ctx.textAlign = "right"; ctx.fillText(a.label, left - 6, Y0 + 11);
      });
      // PL13p: the right-hand axis of a map, each row labelled with its second variable (W3), before the colour bar
      const rlw = rightLabelsWidth(s);
      if (rlw) {
        ctx.font = fontOf(tickPx); ctx.fillStyle = muted; ctx.textAlign = "left"; ctx.strokeStyle = ink;
        const taken = [];
        for (const l of s.rightLabels) {
          const Y = sy(l.at); if (!isFinite(Y) || Y < top - 2 || Y > bottom + 2) continue;
          if (taken.some((y) => Math.abs(y - Y) < tickPx + 2)) continue;
          taken.push(Y);
          ctx.beginPath(); ctx.moveTo(w - right, Math.round(Y) + 0.5); ctx.lineTo(w - right + 4, Math.round(Y) + 0.5); ctx.stroke();
          ctx.fillText(l.text, w - right + 6, Y + tickPx * 0.36);
        }
        if (s.rightTitle) { ctx.font = fontOf(labelPx); ctx.save(); ctx.translate(w - right + rlw - 2, top + (bottom - top) / 2); ctx.rotate(-Math.PI / 2); ctx.textAlign = "center"; ctx.fillText(s.rightTitle, 0, 0); ctx.restore(); ctx.font = fontOf(tickPx); }
      }
      if (s.bar) {
        const bx = w - right + 14 + rlw, bw = 12, n = 64;
        // the bar's axis is the core's range; with a clip its colours show the saturation (one colour beyond each bound)
        const onBar = (v) => bottom - clipPlace(v, zlo, zhi) * (bottom - top);
        for (let q = 0; q < n; q++) { ctx.fillStyle = rampColour(s.ramp, clip ? place(barValue((q + 0.5) / n, zlo, zhi)) : (q + 0.5) / n); ctx.fillRect(bx, bottom - (q + 1) * (bottom - top) / n, bw, (bottom - top) / n + 0.5); }
        ctx.strokeStyle = ink; ctx.globalAlpha = 0.5; ctx.strokeRect(bx + 0.5, top + 0.5, bw - 1, bottom - top - 1); ctx.globalAlpha = 1;
        const by = onBar;
        st.geom.bar = { id: s.id, x0: bx, x1: bx + bw, top, bottom, zlo, zhi, clip: clip ? [cfrom, cto] : null };
        // PL13l: the clip's two marks on the bar, dragged; not in a figure, not in an export
        if (clip && !exporting && !fig) {
          for (const v of [cfrom, cto]) {
            const Y = Math.round(onBar(v)) + 0.5;
            ctx.fillStyle = css("--accent"); ctx.strokeStyle = css("--bg"); ctx.lineWidth = 1;
            ctx.beginPath(); ctx.moveTo(bx - 8, Y - 5); ctx.lineTo(bx - 1, Y); ctx.lineTo(bx - 8, Y + 5); ctx.closePath(); ctx.fill(); ctx.stroke();
            ctx.strokeStyle = ink; ctx.beginPath(); ctx.moveTo(bx, Y); ctx.lineTo(bx + bw, Y); ctx.stroke();
          }
          ctx.lineWidth = 1;
        }
        let bt = s.bar.scale === "log" ? logTicks(zlo, zhi) : s.bar.scale === "sqrt" ? sqrtTicks(zlo, zhi, bottom - top) : axisTicks(zlo, zhi, bottom - top, null);
        // E10c: a log bar narrower than one 1-2-5 step (Explore's χ² about its minimum) labelled by the linear rule
        if (s.bar.scale === "log" && bt.ticks.length < 2) { const lin = axisTicks(Math.pow(10, zlo), Math.pow(10, zhi), bottom - top, null); bt = { ticks: lin.ticks.map((t) => ({ value: Math.log10(t.value), label: t.label })), suffix: lin.suffix }; }
        ctx.fillStyle = muted; ctx.textAlign = "left";
        for (const t of bt.ticks) { const Y = Math.round(by(t.value)) + 0.5; ctx.strokeStyle = ink; ctx.beginPath(); ctx.moveTo(bx + bw, Y); ctx.lineTo(bx + bw + 3, Y); ctx.stroke(); ctx.fillText(t.label, bx + bw + 5, Y + 4); }
        ctx.save(); ctx.translate(w - 6, top + (bottom - top) / 2); ctx.rotate(-Math.PI / 2); ctx.textAlign = "center"; ctx.fillText((s.bar.label || "") + (bt.suffix || ""), 0, 0); ctx.restore();
      }
      // PL13l: the cell cursor, a column and a row of the image marked (the module's row and column cursors); drawn in
      // an export too, as the view shows it (the clip's marks are not: they are handles)
      if (st.cellCursor) {
        const c = st.cellCursor, ac = css("--accent");
        ctx.save(); ctx.beginPath(); ctx.rect(left, top, w - left - right, h - top); ctx.clip();
        ctx.strokeStyle = ac; ctx.fillStyle = ac; ctx.lineWidth = 1;
        if (c.j != null && c.j >= 0 && c.j + 1 < e.length) {
          const X0 = sx(e[c.j]), X1 = sx(e[c.j + 1]);
          ctx.globalAlpha = 0.15; ctx.fillRect(X0, top, Math.max(1, X1 - X0), bottom - top); ctx.globalAlpha = 0.9;
          ctx.setLineDash([3, 3]); ctx.beginPath(); ctx.moveTo(Math.round(X0) + 0.5, top); ctx.lineTo(Math.round(X0) + 0.5, bottom); ctx.moveTo(Math.round(X1) + 0.5, top); ctx.lineTo(Math.round(X1) + 0.5, bottom); ctx.stroke(); ctx.setLineDash([]);
        }
        let Y0 = null, Y1 = null;
        if (typeof c.r === "number" && s.rows[c.r]) { Y0 = Math.max(top, sy(s.rows[c.r].hi)); Y1 = Math.min(bottom, sy(s.rows[c.r].lo)); }
        else if (typeof c.r === "string" && s.apart) { const a = parseInt(c.r.slice(6)); if (s.apart[a]) { Y0 = apartTop + 18 * a; Y1 = Y0 + 14; } }
        if (Y0 != null) {
          ctx.globalAlpha = 0.15; ctx.fillRect(left, Y0, w - left - right, Math.max(1, Y1 - Y0)); ctx.globalAlpha = 0.9;
          ctx.setLineDash([3, 3]); ctx.beginPath(); ctx.moveTo(left, Math.round(Y0) + 0.5); ctx.lineTo(w - right, Math.round(Y0) + 0.5); ctx.moveTo(left, Math.round(Y1) + 0.5); ctx.lineTo(w - right, Math.round(Y1) + 0.5); ctx.stroke(); ctx.setLineDash([]);
        }
        ctx.globalAlpha = 1; ctx.restore();
      }
    }
    ctx.textAlign = "right";
    // bands and regions first
    for (const s of st.series) {
      if (s.hidden || s.kind !== "bands") continue;
      const selB = st.sel.get(s.id);
      s.ranges.forEach(([a, b], i) => {
        if (b < st.x0 || a > st.x1) return;
        const X0 = sx(Math.max(a, st.x0)), X1 = sx(Math.min(b, st.x1));
        ctx.fillStyle = selB && selB.has(i) ? css("--accent") : colourOf(s.color, muted); ctx.globalAlpha = selB && selB.has(i) ? 0.22 : s.hatch ? 0.05 : 0.12;
        ctx.fillRect(X0, top, X1 - X0, bottom - top);
        if (s.hatch) {   // E8 e: a hatched span, its lines at 45°, 6 px apart, so that it reads without its colour
          ctx.save(); ctx.beginPath(); ctx.rect(X0, top, X1 - X0, bottom - top); ctx.clip();
          ctx.strokeStyle = colourOf(s.color, muted); ctx.globalAlpha = 0.45; ctx.lineWidth = 1; ctx.beginPath();
          const H = bottom - top;
          for (let xx = X0 - H; xx < X1; xx += 6) { ctx.moveTo(xx, bottom); ctx.lineTo(xx + H, top); }
          ctx.stroke(); ctx.restore();
        }
        ctx.globalAlpha = 1;
        if (s.editable) { ctx.strokeStyle = selB && selB.has(i) ? css("--accent") : colourOf(s.color, muted); ctx.lineWidth = 1; ctx.setLineDash([3, 3]); ctx.beginPath(); ctx.moveTo(Math.round(sx(a)) + 0.5, top); ctx.lineTo(Math.round(sx(a)) + 0.5, bottom); ctx.moveTo(Math.round(sx(b)) + 0.5, top); ctx.lineTo(Math.round(sx(b)) + 0.5, bottom); ctx.stroke(); ctx.setLineDash([]); }
      });
    }
    // the band being created
    if (st.drag && st.drag.band) { const d = st.drag; ctx.fillStyle = css("--accent"); ctx.globalAlpha = 0.15; ctx.fillRect(Math.min(d.x, d.cx), top, Math.abs(d.cx - d.x), bottom - top); ctx.globalAlpha = 1; }
    // series
    let k = 0;
    for (const s of st.series) {
      if (s.hidden || s.kind === "bands" || s.kind === "image") continue;
      const color = colourOf(s.color, null) || css(PALETTE[k++ % PALETTE.length]);
      const x = s.x, y = s.y, f = fy(s), op = s.opacity == null ? 1 : s.opacity;
      if (s.strip && y) {
        ctx.save(); ctx.beginPath(); ctx.rect(left, stripTop, w - left - right, stripH); ctx.clip();
        ctx.strokeStyle = color; ctx.lineWidth = s.width || 1; ctx.beginPath(); let pen = false;
        for (let i = 0; i < x.length; i++) { if (x[i] < st.x0 - 1e-9 || x[i] > st.x1 + 1e-9) { pen = false; continue; } const v = y[i]; if (!isFinite(v)) { pen = false; continue; } const X = sx(x[i]), Y = sys(v); pen ? ctx.lineTo(X, Y) : ctx.moveTo(X, Y); pen = true; }
        ctx.stroke(); ctx.restore();
        continue;
      }
      if (s.path) {   // E10c (Explore): a path in the plane, x in any order, drawn whole, clipped to the frame; `dash`, `radius`
        ctx.save(); ctx.beginPath(); ctx.rect(left, top, w - left - right, bottom - top); ctx.clip();
        ctx.strokeStyle = color; ctx.lineWidth = s.width || 1.2; ctx.globalAlpha = op; ctx.setLineDash(s.dash || []);
        if (s.kind !== "points") {
          ctx.beginPath(); let pen = false; for (let i = 0; i < x.length; i++) { const v = f(y[i]); if (!isFinite(v)) { pen = false; continue; } pen ? ctx.lineTo(sx(x[i]), sy(v)) : ctx.moveTo(sx(x[i]), sy(v)); pen = true; }
          // E8 b: a patch (the core's ellipse outline), closed and filled
          if (s.fill) { ctx.closePath(); ctx.fillStyle = color; ctx.globalAlpha = op * (s.fillOpacity == null ? 0.15 : s.fillOpacity); ctx.fill(); ctx.globalAlpha = op; }
          ctx.stroke();
        }
        ctx.setLineDash([]); ctx.fillStyle = color; ctx.strokeStyle = ink; ctx.lineWidth = 0.6;
        if (s.kind !== "line") for (let i = 0; i < x.length; i++) { const v = f(y[i]); if (!isFinite(v)) continue; ctx.beginPath(); ctx.arc(sx(x[i]), sy(v), s.radius || 3, 0, 2 * Math.PI); ctx.fill(); ctx.stroke(); }
        // E8 a: the points numbered (a path over the cycles), a number left out where it would touch the last one drawn
        if (s.numbered || s.numbers) {
          ctx.fillStyle = muted; ctx.textAlign = "left"; let last = null;
          for (let i = 0; i < x.length; i++) {
            const v = f(y[i]); if (!isFinite(v)) continue;
            const X = sx(x[i]), Y = sy(v);
            if (last && Math.hypot(X - last[0], Y - last[1]) < tickPx + 2) continue;
            ctx.fillText(String(s.numbers ? s.numbers[i] : i), X + (s.radius || 3) + 2, Y - 3); last = [X, Y];
          }
          ctx.textAlign = "right";
        }
        ctx.restore(); continue;
      }
      let i0 = 0; while (i0 < x.length && x[i0] < st.x0) i0++; i0 = Math.max(0, i0 - 1);
      let i1 = x.length - 1; while (i1 > 0 && x[i1] > st.x1) i1--; i1 = Math.min(x.length - 1, i1 + 1);
      if (s.kind === "handles") {
        const sel = st.sel.get(s.id);
        for (let i = 0; i < x.length; i++) {
          if (x[i] < st.x0 || x[i] > st.x1) continue;
          const v = ty(y[i]); if (!isFinite(v)) continue;
          const X = sx(x[i]), Y = sy(v), r = sel && sel.has(i) ? 5 : 4;
          ctx.fillStyle = sel && sel.has(i) ? css("--accent") : color; ctx.strokeStyle = css("--bg"); ctx.lineWidth = 1;
          ctx.fillRect(X - r, Y - r, 2 * r, 2 * r); ctx.strokeRect(X - r + 0.5, Y - r + 0.5, 2 * r - 1, 2 * r - 1);
        }
        if (!s.line) continue;
      }
      if (s.kind === "markers") {
        const selM = st.sel.get(s.id);
        ctx.strokeStyle = color; ctx.lineWidth = 1;
        ctx.lineWidth = s.width || 1; ctx.globalAlpha = s.alpha == null ? 1 : s.alpha;
        for (let i = i0; i <= i1; i++) {
          const X = Math.round(sx(x[i])) + 0.5;
          // at a height (data units) when the series gives one, else in rows above the axis
          const Yc = s.at != null ? sy(ty(s.at)) : (s.row != null ? bandTop + 8 + 10 * s.row : null);
          const hit = selM && selM.has(i);
          if (hit) { ctx.save(); ctx.strokeStyle = css("--accent"); ctx.lineWidth = 2.5; }
          ctx.beginPath();
          if (Yc != null) { ctx.moveTo(X, Yc - (hit ? 7 : 4)); ctx.lineTo(X, Yc + (hit ? 7 : 4)); } else { ctx.moveTo(X, bottom - 14); ctx.lineTo(X, bottom - 6); }
          ctx.stroke();
          if (hit) ctx.restore();
        }
        ctx.globalAlpha = 1;
        continue;
      }
      if (s.kind === "points" || s.points) {
        const stride = Math.max(1, Math.floor((i1 - i0) / (w * 2)));
        // Settings ▸ Plots ▸ Marker size (§4): the points on screen scaled; a figure's keep their printed size
        const ms = markerScale();
        const sel = st.sel.get(s.id);
        // σ as error bars when fewer than one point per 4 px (otherwise they would blacken the plot)
        const bars = s.sigma && !s.scaled && (i1 - i0) < (w - left - right) / 4;
        for (let i = i0; i <= i1; i += stride) {
          const v = f(y[i]); if (!isFinite(v)) continue;
          const masked = s.sigma && s.sigma[i] == null;
          const X = sx(x[i]), Y = sy(v);
          if (bars && !masked && s.sigma[i] > 0) {
            const ya = ty(y[i] - s.sigma[i]), yb = ty(y[i] + s.sigma[i]);
            if (isFinite(ya) && isFinite(yb)) { ctx.strokeStyle = color; ctx.globalAlpha = 0.5; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(X + 0.5, sy(ya)); ctx.lineTo(X + 0.5, sy(yb)); ctx.moveTo(X - 2, sy(ya) + 0.5); ctx.lineTo(X + 3, sy(ya) + 0.5); ctx.moveTo(X - 2, sy(yb) + 0.5); ctx.lineTo(X + 3, sy(yb) + 0.5); ctx.stroke(); ctx.globalAlpha = 1; }
          }
          if (sel && sel.has(i)) { ctx.fillStyle = css("--accent"); ctx.fillRect(X - 2.5, Y - 2.5, 5, 5); continue; }
          if (s.marker && !masked) {   // PL13g: a series' marker shape, hollow for an open point (C7c: one shape per window)
            const r = fig ? 3 : 2.5 * ms; ctx.beginPath();
            if (s.marker === "circle") ctx.arc(X, Y, r, 0, 2 * Math.PI);
            else if (s.marker === "diamond") { ctx.moveTo(X, Y - r - 1); ctx.lineTo(X + r + 1, Y); ctx.lineTo(X, Y + r + 1); ctx.lineTo(X - r - 1, Y); ctx.closePath(); }
            else if (s.marker === "triangle") { ctx.moveTo(X, Y - r - 1); ctx.lineTo(X + r + 1, Y + r); ctx.lineTo(X - r - 1, Y + r); ctx.closePath(); }
            else ctx.rect(X - r, Y - r, 2 * r, 2 * r);
            ctx.globalAlpha = op; ctx.lineWidth = 1; ctx.strokeStyle = color;
            if (s.open && s.open[i]) { ctx.fillStyle = fig ? "#ffffff" : css("--bg"); ctx.fill(); ctx.stroke(); } else { ctx.fillStyle = color; ctx.fill(); }
            ctx.globalAlpha = 1; continue;
          }
          if (masked) { ctx.strokeStyle = color; ctx.globalAlpha = 0.5; ctx.strokeRect(X - 2, Y - 2, 4, 4); ctx.globalAlpha = 1; }
          else if (s.open && s.open[i]) { ctx.strokeStyle = color; ctx.lineWidth = 1; ctx.fillStyle = css("--bg"); ctx.globalAlpha = op; ctx.fillRect(X - 2.5, Y - 2.5, 5, 5); ctx.strokeRect(X - 2.5, Y - 2.5, 5, 5); ctx.globalAlpha = 1; }   // PL13: an open point (F5: below 2σ, the core's flag)
          else if (fig) { ctx.fillStyle = color; ctx.globalAlpha = op; ctx.fillRect(X - 2, Y - 2, 4, 4); ctx.globalAlpha = 1; }   // PL13g: a figure's points larger
          else { ctx.fillStyle = color; ctx.globalAlpha = op; ctx.fillRect(X - ms, Y - ms, 2 * ms, 2 * ms); ctx.globalAlpha = 1; }
        }
      }
      if (s.kind === "line" || s.kind === "both" || (s.kind === "handles" && s.line)) {
        const cols = Math.max(1, Math.round(w - left - right));
        // PL13: filled down to the frame in the background (a stack's hidden parts, F2), before the line itself
        if (s.fillBelow) {
          const k2 = Math.max(1, Math.floor((i1 - i0) / (cols * 2)));
          ctx.save(); ctx.beginPath(); ctx.rect(left, top, w - left - right, bottom - top); ctx.clip();
          ctx.beginPath(); let firstX = null, lastX = null;
          for (let i = i0; i <= i1; i += k2) { const v = f(y[i]); if (!isFinite(v)) continue; const X = sx(x[i]), Y = sy(v); if (firstX == null) { ctx.moveTo(X, bottom); ctx.lineTo(X, Y); firstX = X; } else ctx.lineTo(X, Y); lastX = X; }
          if (firstX != null) { ctx.lineTo(lastX, bottom); ctx.closePath(); ctx.fillStyle = colourOf(s.fillColor || "--bg", "#ffffff"); ctx.globalAlpha = 1; ctx.fill(); }
          ctx.restore();
        }
        ctx.strokeStyle = color; ctx.lineWidth = s.width || 1.2; ctx.globalAlpha = op; ctx.beginPath(); let pen = false;
        if (i1 - i0 > 3 * cols) {
          // PL13: more points than pixels: each pixel column's first, lowest, highest and last point, in order, which
          // draws what the full line would show at this width and keeps hundreds of series quick
          let col = null, first = 0, lo = 0, hi = 0, lastY = 0;
          const flush = () => { if (col == null) return; const X = left + col + 0.5; pen ? ctx.lineTo(X, first) : ctx.moveTo(X, first); ctx.lineTo(X, lo); ctx.lineTo(X, hi); ctx.lineTo(X, lastY); pen = true; };
          for (let i = i0; i <= i1; i++) {
            const v = f(y[i]); if (!isFinite(v)) { flush(); col = null; pen = false; continue; }
            const Y = sy(v), c = Math.floor(sx(x[i]) - left);
            if (c !== col) { flush(); col = c; first = lo = hi = lastY = Y; } else { if (Y < lo) lo = Y; if (Y > hi) hi = Y; lastY = Y; }
          }
          flush();
        } else {
          for (let i = i0; i <= i1; i++) { const v = f(y[i]); if (!isFinite(v)) { pen = false; continue; } const X = sx(x[i]), Y = sy(v); pen ? ctx.lineTo(X, Y) : ctx.moveTo(X, Y); pen = true; }
        }
        ctx.stroke(); ctx.globalAlpha = 1;
      }
      if (s.kind === "vlines") {
        ctx.strokeStyle = color; ctx.lineWidth = 1; ctx.setLineDash([3, 3]);
        for (let i = 0; i < x.length; i++) { if (x[i] < st.x0 || x[i] > st.x1) continue; const X = Math.round(sx(x[i])) + 0.5; ctx.beginPath(); ctx.moveTo(X, top); ctx.lineTo(X, bottom); ctx.stroke(); }
        ctx.setLineDash([]);
      }
    }
    // PL13: labels at the series' baselines, in the right margin, in each series' colour and opacity; a label that
    // would touch one already drawn is left out, those with the higher `priority` and then the brightest drawn first
    ctx.textAlign = "left";
    k = 0;
    const labels = [];
    for (const s of st.series) {
      if (s.hidden || s.kind === "bands" || s.kind === "image") continue;
      const color = colourOf(s.color, null) || css(PALETTE[k++ % PALETTE.length]);
      if (!s.baselineLabel) continue;
      const Y = sy(s.scaled ? s.baselineLabel.at : ty(s.baselineLabel.at));
      if (!isFinite(Y) || Y < top - 2 || Y > bottom + 2) continue;
      labels.push({ Y, color, op: s.opacity == null ? 1 : s.opacity, text: s.baselineLabel.text, rank: s.baselineLabel.priority || 0 });
    }
    // PL13 (Plots, step 2): labels at a series' own points (`pointLabels: [{ x, y, text }]`, y on the series' scale), as
    // the projection's oblique axis carries its temperatures; thinned as the baseline labels are
    for (const s of st.series) {
      if (s.hidden || !s.pointLabels) continue;
      const f = fy(s), color = colourOf(s.color, css("--muted"));
      for (const p of s.pointLabels) {
        const X = sx(p.x), Y = sy(f(p.y));
        if (!isFinite(X) || !isFinite(Y) || Y < top - 2 || Y > bottom + 2 || X < left || X > w - 4) continue;
        if (p.tick) { ctx.strokeStyle = color; ctx.beginPath(); ctx.moveTo(Math.round(X) + 0.5, Math.round(Y) + 0.5); ctx.lineTo(Math.round(X) + 5.5, Math.round(Y) + 0.5); ctx.stroke(); }   // PL13f: the axis's tick
        labels.push(p.above ? { Y: Y - 14, X: X - 4, color, op: 1, text: p.text, rank: p.priority || 0 } : { Y, X: X + (p.tick ? 9 : 6), color, op: 1, text: p.text, rank: p.priority || 0 });
      }
    }
    const taken = [];
    for (const l of labels.sort((a, b) => b.rank - a.rank || b.op - a.op)) {
      if (taken.some((y) => Math.abs(y - l.Y) < 11)) continue;
      taken.push(l.Y);
      ctx.fillStyle = l.color; ctx.globalAlpha = Math.max(0.45, l.op);
      ctx.fillText(l.text, l.X == null ? w - right + 6 : l.X, l.Y + 4); ctx.globalAlpha = 1;
    }
    // PL13 (step 3): a cursor at x set from outside (another plot's pointer, F6), dashed across the frame and the strip
    if (st.cursor != null && st.cursor >= st.x0 && st.cursor <= st.x1) {
      const X = Math.round(sx(st.cursor)) + 0.5;
      ctx.strokeStyle = css("--accent"); ctx.lineWidth = 1; ctx.setLineDash([4, 3]); ctx.beginPath(); ctx.moveTo(X, top); ctx.lineTo(X, stripH ? stripBottom : bottom); ctx.stroke(); ctx.setLineDash([]);
    }
    // the box being dragged
    if (st.drag && st.drag.box) {
      const d = st.drag;
      ctx.strokeStyle = css("--accent"); ctx.lineWidth = 1; ctx.setLineDash([4, 3]);
      ctx.strokeRect(Math.min(d.x, d.cx) + 0.5, Math.min(d.y, d.cy) + 0.5, Math.abs(d.cx - d.x), Math.abs(d.cy - d.y)); ctx.setLineDash([]);
    }
    if (!fig || images.length) { ctx.strokeStyle = ink; ctx.globalAlpha = 0.6; ctx.strokeRect(left + 0.5, top + 0.5, w - left - right - 1, bottom - top - 1); ctx.globalAlpha = 1; }
    renderLegend();
    if (opts.onView && (st.lastView == null || st.lastView[0] !== st.x0 || st.lastView[1] !== st.x1)) { st.lastView = [st.x0, st.x1]; opts.onView(st.x0, st.x1); }
    // linked axes (PL13, F6): every listener told of a change of the view, the y range too
    const vk = `${st.x0}|${st.x1}|${st.y0}|${st.y1}`;
    if (vk !== st.lastViewKey) { st.lastViewKey = vk; for (const fn of viewListeners) { try { fn(getView()); } catch (e) { /* a listener's own */ } } }
  }
  /** PL13p: the width an image's right-hand labels take (their widest text, the tick, the title's line), 0 without. */
  function rightLabelsWidth(s) {
    if (!s.rightLabels || !s.rightLabels.length) return 0;
    const widest = s.rightLabels.reduce((m, l) => Math.max(m, ctx.measureText(l.text).width), 0);
    return Math.round(6 + widest + 6 + (s.rightTitle ? 14 : 0));
  }
  /** PL30 f: a common factor said in the unit for a figure ("intensity (10³ counts)"), as papers write it; on screen,
      and with an offset, the axis rule's own suffix after the label. */
  function withFactor(label, suffix, fig) {
    const m = /^ ×(10[⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+)$/.exec(suffix || "");
    if (!fig || !m) return label + (suffix || "");
    const u = /^(.*)\(([^()]*)\)\s*$/.exec(label);
    return u ? `${u[1]}(${m[1]} ${u[2]})` : `${label} (×${m[1]})`;
  }
  function renderLegend() {
    if (!legendEl) return;
    let k = 0;
    legendEl.innerHTML = st.series.filter((s) => s.label && s.kind !== "bands" && s.legend !== false).map((s) => {
      const color = colourOf(s.color, null) || css(PALETTE[k++ % PALETTE.length]);
      return `<span class="key${s.hidden ? " off" : ""}" data-id="${s.id}" data-tour="legend:${s.id}"><span class="sw" style="background:${color}"></span>${s.label}</span>`;
    }).join("");
    legendEl.querySelectorAll(".key").forEach((e) => { e.onclick = () => { const s = st.series.find((q) => q.id === e.dataset.id); if (s) { s.hidden = !s.hidden; draw(); } }; });
  }

  // ------------------------------------------------------------------ the API
  /** §10: what the plot shows, in words, for VoiceOver (the plot is an image to it): its series by name, against its
      axes. Each series is named in the legend too, so that none is told apart by its colour alone. */
  function describe() {
    const names = st.series.filter((s) => !s.hidden && s.label).map((s) => s.label);
    const scale = { sqrt: ", on a square-root scale", log: ", on a logarithmic scale", linear: "" }[st.yScale] || "";
    if (!canvas.setAttribute) return;   // a test's stand-in canvas
    canvas.setAttribute("role", "img");
    canvas.setAttribute("aria-label", `Plot of ${names.length ? names.join("; ") : "nothing yet"}: ${st.yLabel}${scale} against ${st.xLabel}`);
  }
  function setSeries(series, keepView = false) {
    st.series = series.map((s) => ({ kind: "both", ...s }));
    if (!keepView || st.fit) { [st.x0, st.x1] = xRange(); st.y0 = st.y1 = null; st.fit = true; }
    describe(); draw();
  }
  function setAxes({ xLabel, yLabel, yScale, look, xSpan, scaleBar, keepView }) { if (xLabel != null) st.xLabel = xLabel; if (yLabel != null) st.yLabel = yLabel; if (yScale) st.yScale = yScale; if (look !== undefined) st.look = look; if (xSpan !== undefined) st.xSpan = xSpan; if (scaleBar !== undefined) st.scaleBar = scaleBar; if (!keepView) st.y0 = st.y1 = null; describe(); draw(); }   // PL13m: keepView keeps a y zoom
  function fitAll() { [st.x0, st.x1] = xRange(); st.y0 = st.y1 = null; st.fit = true; draw(); }
  const modeCursor = (mode) => (mode === "select" ? "cell" : mode === "pick" ? "copy" : mode === "region" ? "col-resize" : "crosshair");
  function setMode(mode) { st.mode = mode; canvas.style.cursor = modeCursor(mode); }
  /** PL13l: a column and a row of the image marked ({ j, r }, r an index of `rows` or "apart:k"), or none (null). */
  function setCellCursor(c) { const k = c ? `${c.j}|${c.r}` : ""; if (k === st.cellCursorKey) return; st.cellCursorKey = k; st.cellCursor = c; draw(); }
  function setSelection(id, indices) { st.sel.set(id, new Set(indices)); draw(); }
  function clearSelection() { st.sel.clear(); draw(); }
  function getView() { return { x0: st.x0, x1: st.x1, y0: st.y0, y1: st.y1 }; }
  function setView(v) { if (v) { st.x0 = v.x0; st.x1 = v.x1; st.y0 = v.y0; st.y1 = v.y1; st.fit = false; draw(); } }
  /** The plot as it is on screen, drawn once more into SVG with its legend (File ▸ Export ▸ Plot, L134). */
  function toSVG(size = null) {
    const w = size && size.w ? size.w : canvas.clientWidth, h = size && size.h ? size.h : canvas.clientHeight;
    const keptFont = fontPt;
    if (size && (size.label_pt || size.tick_pt)) fontPt = { label_pt: size.label_pt || size.tick_pt, tick_pt: size.tick_pt || size.label_pt };   // PL30 a
    const bg = getComputedStyle(canvas).backgroundColor;
    const rec = svgRecorder(w, h, st.look === "figure" ? "#ffffff" : bg && bg !== "rgba(0, 0, 0, 0)" ? bg : css("--panel") || "#ffffff");   // PL13g: a figure on white
    const saved = ctx;
    ctx = rec; exporting = { w, h };
    try {
      draw();
      drawLegend(rec, legendItems(legendEl), 64, 18, css("--ink") || "#000");
    } finally { ctx = saved; exporting = null; fontPt = keptFont; draw(); }
    return svgText(rec);
  }
  /** The plot as drawn, a PNG on the pasteboard. In the app AppKit puts it there (`clipboard_png`): the web view's
      clipboard refuses a write no user's action started, as a menu chosen while the window is behind (0.33.0's
      installed self-check); in a browser, the page's clipboard. */
  async function copyImage() {
    return new Promise((resolve, reject) => canvas.toBlob(async (blob) => {
      try {
        // W6 (PL13q): the clipboard through the host (the page's, else Tauri's), the browser's when there is none
        const h = window.__DRS_HOST__ || (window.__TAURI__ ? { invoke: (c, a) => window.__TAURI__.core.invoke(c, a) } : null);
        if (h && typeof h.invoke === "function") {
          const png = await new Promise((ok, no) => { const r = new FileReader(); r.onload = () => ok(String(r.result).split(",")[1]); r.onerror = () => no(r.error); r.readAsDataURL(blob); });
          await h.invoke("clipboard_png", { png });
        }
        else await navigator.clipboard.write([new ClipboardItem({ "image/png": blob })]);
        resolve(true);
      } catch (e) { reject(e); }
    }, "image/png"));
  }
  function toX(px) { const g = st.geom; return g ? st.x0 + (px - g.left) / (g.w - g.left - g.right) * (st.x1 - st.x0) : 0; }
  /** The handle under the pointer (within 7 px), the nearest first. */
  function handleAt(px, py) {
    const g = st.geom; if (!g) return null;
    let best = null, bd = 49;
    for (const s of st.series) {
      if (s.hidden || s.kind !== "handles") continue;
      for (let i = 0; i < s.x.length; i++) { const v = ty(s.y[i]); if (!isFinite(v)) continue; const dx = g.sx(s.x[i]) - px, dy = g.sy(v) - py, d2 = dx * dx + dy * dy; if (d2 < bd) { bd = d2; best = { s, i }; } }
    }
    return best;
  }
  /** An editable band's edge under the pointer (within 5 px): { s, i, edge: 0 | 1 }, or the band itself: { s, i, edge: -1 }. */
  function bandAt(px) {
    const g = st.geom; if (!g) return null;
    for (const s of st.series) {
      if (s.hidden || s.kind !== "bands") continue;
      for (let i = 0; i < s.ranges.length; i++) {
        const [a, b] = s.ranges[i];
        if (s.editable && Math.abs(g.sx(a) - px) <= 5) return { s, i, edge: 0 };
        if (s.editable && Math.abs(g.sx(b) - px) <= 5) return { s, i, edge: 1 };
      }
    }
    for (const s of st.series) {
      if (s.hidden || s.kind !== "bands") continue;
      for (let i = 0; i < s.ranges.length; i++) { const [a, b] = s.ranges[i]; if (px >= g.sx(a) && px <= g.sx(b)) return { s, i, edge: -1 }; }
    }
    return null;
  }
  /** A marker tick under the pointer (within 4 px in x, and in y when the row sits at a height). */
  function markerAt(px, py) {
    const g = st.geom; if (!g) return null;
    let best = null, bd = 16;
    for (const s of st.series) {
      if (s.hidden || s.kind !== "markers") continue;
      if (s.at != null) { const Y = g.sy(ty(s.at)); if (Math.abs(Y - py) > 10) continue; }
      else if (s.row != null) { const Y = g.bandTop + 8 + 10 * s.row; if (Math.abs(Y - py) > 7) continue; }
      for (let i = 0; i < s.x.length; i++) { const d = Math.abs(g.sx(s.x[i]) - px); if (d * d < bd) { bd = d * d; best = { s, i }; } }
    }
    return best;
  }
  /** PL13l: the image's cell under the pointer, { s, j, r } (r an index of `rows`, or "apart:k"), the cell itself and
      never a blend (mapclip.js's `cellAt`), or null. */
  function imageCell(px, py) {
    const g = st.geom; if (!g || px < g.left || px > g.w - g.right) return null;
    for (const s of st.series) {
      if (s.hidden || s.kind !== "image") continue;
      const xv = toX(px), e = s.xEdges;
      if (py >= g.top && py <= g.bottom) { const c = cellAt(xv, toY(py), e, s.rows); if (c) return { s, j: c.j, r: c.r }; }
      const j = columnAt(xv, e), a = Math.floor((py - g.apartTop) / 18);   // a row apart: its band under the frame
      if (j >= 0 && s.apart && a >= 0 && a < s.apart.length && py - g.apartTop - 18 * a <= 14) return { s, j, r: "apart:" + a };
    }
    return null;
  }
  /** PL13: the readout over an image's cell (the series' own `cellText`), or null. */
  function imageCellAt(px, py) {
    const c = imageCell(px, py);
    return c && c.s.cellText ? c.s.cellText(c.j, c.r) : null;
  }
  /** PL13l: the clip's mark under the pointer (within 8 px of it, beside the bar): { s, mark: 0 | 1 }, or null. */
  function barMarkAt(px, py) {
    const g = st.geom; const b = g && g.bar; if (!b || !b.clip || st.look === "figure") return null;
    if (px < b.x0 - 12 || px > b.x1 + 4) return null;
    const s = st.series.find((q) => q.id === b.id && q.kind === "image"); if (!s || !s.bar || !s.bar.clip) return null;
    const ys = b.clip.map((v) => b.bottom - clipPlace(v, b.zlo, b.zhi) * (b.bottom - b.top));
    let best = null, bd = 8;
    ys.forEach((Y, k) => { const d = Math.abs(Y - py); if (d <= bd) { bd = d; best = k; } });
    return best == null ? null : { s, mark: best };
  }
  /** PL13l: whether the pointer is on the colour bar (a double-click there resets the clip). */
  function onBarAt(px, py) {
    const b = st.geom && st.geom.bar; if (!b || !b.clip) return null;
    if (px < b.x0 - 12 || px > b.x1 + 4 || py < b.top || py > b.bottom) return null;
    return st.series.find((q) => q.id === b.id && q.kind === "image") || null;
  }
  function toY(py) { const g = st.geom; return g ? g.y0 + (g.bottom - py) / (g.bottom - g.top) * (g.y1 - g.y0) : 0; }
  /** PL13 (`opts.nearestInY`): the curve nearest the pointer in x and y, for a stack or an overlay of many series. */
  function nearestXY(x, py) {
    const g = st.geom; let best = null, bd = Infinity;
    for (const s of st.series) {
      if (s.hidden || !s.y || s.kind !== "line" && s.kind !== "both" && s.kind !== "points") continue;
      let a = 0, b = s.x.length - 1; if (b < 0) continue;
      while (a < b) { const m = (a + b) >> 1; if (s.x[m] < x) a = m + 1; else b = m; }
      const i = a > 0 && Math.abs(s.x[a - 1] - x) < Math.abs(s.x[a] - x) ? a - 1 : a;
      const v = fy(s)(s.y[i]); if (!isFinite(v)) continue;
      const d = Math.abs(g.sy(v) - py); if (d < bd) { bd = d; best = { s, i }; }
    }
    return best;
  }
  function nearest(x) {
    // the nearest point of the first visible series with y, by x
    for (const s of st.series) {
      if (s.hidden || !s.y || s.kind === "markers" || s.kind === "vlines" || s.kind === "bands") continue;
      let a = 0, b = s.x.length - 1; while (a < b) { const m = (a + b) >> 1; if (s.x[m] < x) a = m + 1; else b = m; }
      const i = a > 0 && Math.abs(s.x[a - 1] - x) < Math.abs(s.x[a] - x) ? a - 1 : a;
      return { s, i };
    }
    return null;
  }

  // ------------------------------------------------------------------ the mouse and the trackpad
  canvas.addEventListener("wheel", (e) => {
    if (!st.geom) return; e.preventDefault();
    const [f0, f1] = xRange();
    // PL13m (Plots, step 17): the wheel over an axis zooms that axis alone about the pointer: the y axis's margin for
    // y (a map's temperature), the x axis's band for x; inside the frame the wheel works as before
    if (opts.axisWheel) {
      const g = st.geom;
      if (e.offsetX < g.left && e.offsetY >= g.top && e.offsetY <= g.bottom) {
        const [full0, full1] = yRange(f0, f1);
        const [ya, yb] = st.y0 == null ? yRange(st.x0, st.x1) : [st.y0, st.y1];
        const t = toY(e.offsetY), f = Math.exp(e.deltaY * 0.0015);
        let a = t - (t - ya) * f, b = t + (yb - t) * f;
        if (b - a < (full1 - full0) * 1e-4) return;
        st.y0 = Math.max(full0, a); st.y1 = Math.min(full1, b); st.fit = false; draw(); return;
      }
      if (e.offsetY > g.bottom && e.offsetX >= g.left && e.offsetX <= g.w - g.right) {
        const t = toX(e.offsetX), f = Math.exp(e.deltaY * 0.0015);
        let a = t - (t - st.x0) * f, b = t + (st.x1 - t) * f;
        if (b - a < (f1 - f0) * 1e-4) return;
        st.x0 = Math.max(f0, a); st.x1 = Math.min(f1, b); st.fit = false; draw(); return;
      }
    }
    if (e.ctrlKey || e.metaKey) {   // a pinch (WebKit sends it as a wheel with ctrlKey), or ⌘-wheel: zoom about the cursor
      const t = toX(e.offsetX), f = Math.exp(e.deltaY * 0.01);
      let a = t - (t - st.x0) * f, b = t + (st.x1 - t) * f;
      if (b - a < (f1 - f0) * 1e-4) return;
      st.x0 = Math.max(f0, a); st.x1 = Math.min(f1, b); st.fit = false; draw(); return;
    }
    // two-finger scroll: pan
    const span = st.x1 - st.x0;
    const dx = e.deltaX * span / (st.geom.w - st.geom.left - st.geom.right);
    let a = st.x0 + dx; a = Math.max(f0, Math.min(a, f1 - span));
    st.x0 = a; st.x1 = a + span;
    if (e.deltaY && st.y0 != null) { const dy = e.deltaY * (st.y1 - st.y0) / (st.geom.bottom - st.geom.top); st.y0 += dy; st.y1 += dy; }
    else if (e.deltaY && Math.abs(e.deltaY) > Math.abs(e.deltaX)) { const t = toX(e.offsetX), f = Math.exp(e.deltaY * 0.0015); let a2 = t - (t - st.x0) * f, b2 = t + (st.x1 - t) * f; if (b2 - a2 >= (f1 - f0) * 1e-4) { st.x0 = Math.max(f0, a2); st.x1 = Math.min(f1, b2); } }
    st.fit = false; draw();
  }, { passive: false });
  canvas.addEventListener("mousedown", (e) => {
    if (!st.geom) return;
    // PL13l: a clip's mark on the colour bar taken by the pointer
    const bm = e.button === 0 ? barMarkAt(e.offsetX, e.offsetY) : null;
    if (bm) { st.drag = { clipMark: bm.mark, x: e.offsetX, y: e.offsetY, cx: e.offsetX, cy: e.offsetY, moved: false }; return; }
    const h = handleAt(e.offsetX, e.offsetY);
    if (h && e.button === 0) { st.drag = { handle: h, x: e.offsetX, y: e.offsetY, cx: e.offsetX, cy: e.offsetY, moved: false }; return; }
    const bd = st.mode === "region" || st.mode === "zoom" ? bandAt(e.offsetX) : null;
    if (bd && bd.edge >= 0 && e.button === 0) { st.drag = { bandEdge: bd, x: e.offsetX, y: e.offsetY, cx: e.offsetX, cy: e.offsetY, moved: false, a: bd.s.ranges[bd.i][0], b: bd.s.ranges[bd.i][1] }; return; }
    const pan = e.button === 1 || e.shiftKey || spaceDown;
    st.drag = { x: e.offsetX, y: e.offsetY, cx: e.offsetX, cy: e.offsetY, x0: st.x0, x1: st.x1, box: !pan && st.mode !== "region" && st.mode !== "pick", band: !pan && st.mode === "region", pan, click: true };
  });
  window.addEventListener("mouseup", (e) => {
    const d = st.drag; if (!d) return; st.drag = null;
    if (d.clipMark != null) { draw(); return; }   // PL13l: the clip's mark let go; every move already told the caller
    if (d.handle) {
      const { s, i } = d.handle;
      if (d.moved) { if (opts.onHandleMove) opts.onHandleMove(s.id, i, s.x[i], s.y[i]); }
      else { const cur = st.sel.get(s.id); const next = e.metaKey && cur ? new Set(cur) : new Set(); if (next.has(i)) next.delete(i); else next.add(i); st.sel.set(s.id, next); if (opts.onHandleSelect) opts.onHandleSelect(s.id, [...next]); }
      draw(); return;
    }
    if (d.bandEdge) {
      const { s, i } = d.bandEdge; const [a, b] = s.ranges[i];
      if (d.moved && opts.onBandChange) opts.onBandChange(s.id, i, Math.min(a, b), Math.max(a, b));
      else { st.sel.set(s.id, new Set([i])); if (opts.onBandSelect) opts.onBandSelect(s.id, i); }
      draw(); return;
    }
    const still = Math.abs(d.cx - d.x) <= 3 && Math.abs(d.cy - d.y) <= 3;
    if (still && d.click) {
      // a click: a tick, a band, or a point picked
      const m = markerAt(d.x, d.y);
      if (m) { st.sel.set(m.s.id, new Set([m.i])); if (opts.onMarker) opts.onMarker(m.s.id, m.i); draw(); return; }
      if (st.mode === "pick") { if (opts.onPick) opts.onPick(toX(d.x), toY(d.y)); return; }
      // PL13l: a click on an image's cell (the module pins or releases its cursors there)
      if (opts.onCellClick && st.mode === "zoom") { const c = imageCell(d.x, d.y); if (c) { opts.onCellClick({ id: c.s.id, j: c.j, r: c.r }); return; } }
      const bd = bandAt(d.x);
      if (bd && bd.edge < 0 && (st.mode === "region" || bd.s.editable)) { st.sel.set(bd.s.id, new Set([bd.i])); if (opts.onBandSelect) opts.onBandSelect(bd.s.id, bd.i); draw(); return; }
      if (st.mode === "region") { for (const s of st.series) if (s.kind === "bands") st.sel.delete(s.id); draw(); return; }
      // a plain click on the pattern in zoom mode ("who is this peak", W7)
      if (st.mode === "zoom" && opts.onClick && d.y >= st.geom.top && d.y <= st.geom.bottom) opts.onClick(toX(d.x), toY(d.y));
      return;
    }
    if (d.band && Math.abs(d.cx - d.x) > 3) { const a = toX(Math.min(d.x, d.cx)), b = toX(Math.max(d.x, d.cx)); if (opts.onBandCreate) opts.onBandCreate(a, b); draw(); return; }
    if (d.box && Math.abs(d.cx - d.x) > 4 && Math.abs(d.cy - d.y) > 4) {
      const xa = toX(Math.min(d.x, d.cx)), xb = toX(Math.max(d.x, d.cx)), ya = toY(Math.max(d.y, d.cy)), yb = toY(Math.min(d.y, d.cy));
      if (st.mode === "select") {
        const out = {};
        for (const s of st.series) {
          if (s.hidden || !s.y || s.kind === "markers" || s.kind === "vlines" || s.kind === "bands" || s.noSelect) continue;
          const idx = [];
          for (let i = 0; i < s.x.length; i++) { const v = ty(s.y[i]); if (s.x[i] >= xa && s.x[i] <= xb && v >= ya && v <= yb) idx.push(i); }
          if (idx.length) out[s.id] = idx;
        }
        for (const id of Object.keys(out)) st.sel.set(id, new Set(e.metaKey && st.sel.has(id) ? [...st.sel.get(id), ...out[id]] : out[id]));
        if (!Object.keys(out).length && !e.metaKey) st.sel.clear();
        draw();
        if (opts.onSelect) opts.onSelect(Object.fromEntries([...st.sel].map(([k, v]) => [k, [...v]])));
      } else { st.x0 = xa; st.x1 = xb; st.y0 = ya; st.y1 = yb; st.fit = false; }
    }
    draw();
  });
  let spaceDown = false;
  window.addEventListener("keydown", (e) => { if (e.code === "Space" && document.activeElement === document.body) spaceDown = true; });
  window.addEventListener("keyup", (e) => { if (e.code === "Space") spaceDown = false; });
  canvas.addEventListener("mousemove", (e) => {
    if (!st.geom) return;
    const d = st.drag;
    if (d && d.clipMark != null) {
      // PL13l: the mark dragged along the bar: the value at its height on the core's range, the pair kept in order
      const b = st.geom.bar; const s = b && st.series.find((q) => q.id === b.id && q.kind === "image");
      if (!s || !s.bar || !s.bar.clip) return;
      d.moved = true;
      const f = (b.bottom - e.offsetY) / (b.bottom - b.top);
      const [from, to] = clipMoved(s.bar.clip, d.clipMark, barValue(f, b.zlo, b.zhi), b.zlo, b.zhi);
      if (s.bar.clip.onChange) s.bar.clip.onChange(from, to); else { s.bar.clip.from = from; s.bar.clip.to = to; draw(); }
      return;
    }
    if (!d) {   // PL13l: the hand on a clip's mark; the mode's cursor back once it leaves (a module's own cursor untouched)
      const onMark = !!barMarkAt(e.offsetX, e.offsetY);
      if (onMark && !st.onMark) canvas.style.cursor = "row-resize"; else if (!onMark && st.onMark) canvas.style.cursor = modeCursor(st.mode);
      st.onMark = onMark;
    }
    if (d && d.handle) {
      const { s, i } = d.handle;
      if (Math.abs(e.offsetX - d.x) > 2 || Math.abs(e.offsetY - d.y) > 2) d.moved = true;
      if (d.moved) { s.x[i] = toX(e.offsetX); s.y[i] = Tinv[st.yScale](toY(e.offsetY)); if (s.keepOrder) { /* the module re-sorts on move end */ } draw(); if (opts.onHandleDrag) opts.onHandleDrag(s.id, i, s.x[i], s.y[i]); }   // E10l: while dragging (Explore's dragged point, §2.7)
      return;
    }
    if (d && d.bandEdge) {
      const { s, i, edge } = d.bandEdge;
      if (Math.abs(e.offsetX - d.x) > 2) d.moved = true;
      if (d.moved) { const r = s.ranges[i]; r[edge] = toX(e.offsetX); draw(); }
      return;
    }
    if (d) {
      if (d.pan) {
        const [f0, f1] = xRange();
        const dx = (toX(d.x) - toX(e.offsetX)); const span = d.x1 - d.x0;
        let a = Math.min(Math.max(f0, d.x0 + dx), f1 - span); st.x0 = a; st.x1 = a + span; st.fit = false;
      } else { d.cx = e.offsetX; d.cy = e.offsetY; }
      draw();
    }
    const x = toX(e.offsetX), n = opts.nearestInY ? nearestXY(x, e.offsetY) : nearest(x);
    const ic = imageCell(e.offsetX, e.offsetY);
    if (opts.onCell) opts.onCell(ic ? { id: ic.s.id, j: ic.j, r: ic.r } : null);   // PL13l
    const cell = ic && ic.s.cellText ? ic.s.cellText(ic.j, ic.r) : null;
    if (readout && cell) { readout.textContent = cell; if (opts.onHover) opts.onHover(null); return; }
    if (readout) {
      const g = st.geom;
      readout.textContent = n ? `${st.xLabel.split(" ")[0]} ${n.s.x[n.i].toFixed(4)}  ${n.s.label ? n.s.label + " " : ""}${(n.s.values || n.s.y)[n.i].toFixed(2)}${n.s.sigma ? (n.s.sigma[n.i] == null ? "  masked" : `  σ ${n.s.sigma[n.i].toFixed(2)}`) : ""}` : `${x.toFixed(4)}, ${toY(e.offsetY).toFixed(2)}`;
      void g;
    }
    if (opts.onHover) opts.onHover(n);
    if (opts.onPointer) opts.onPointer(x);
  });
  canvas.addEventListener("mouseleave", () => { if (readout) readout.textContent = ""; if (opts.onPointer) opts.onPointer(null); if (opts.onCell) opts.onCell(null); });
  // PL13l: a double-click on the colour bar returns the clip to the core's range; elsewhere it fits the plot
  canvas.addEventListener("dblclick", (e) => { const s = onBarAt(e.offsetX, e.offsetY); if (s) { if (s.bar.clip.onReset) s.bar.clip.onReset(); return; } fitAll(); });
  new ResizeObserver(draw).observe(canvas);
  window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", draw);
  window.addEventListener("drs-markers", () => draw());
  window.matchMedia("(prefers-contrast: more)").addEventListener("change", draw);   // §10: Increase contrast's colours
  function selection(id) { return [...(st.sel.get(id) || [])]; }
  /** PL13 (step 3): a cursor drawn at x (null: none), set from another plot's pointer (linked cursors, F6). */
  function setCursor(x) { if (st.cursor === x) return; st.cursor = x; draw(); }
  /** PL13 (F6): `fn(view)` after every change of the view (a zoom, a pan, a fit), for axes linked between plots. */
  function onViewChange(fn) { viewListeners.push(fn); return () => { const i = viewListeners.indexOf(fn); if (i >= 0) viewListeners.splice(i, 1); }; }
  /** PL30 a: the text's sizes in points, or null for the screen's. */
  function setFont(f) { fontPt = f || null; draw(); }
  /** TU5: a point in data units (x, and y in counts, or the frame's middle) in the window's coordinates, for the tour
      player; null before the first drawing. Also the frame's top and bottom, for a height given as a fraction. */
  function pixelOf(x, y = null) {
    const g = st.geom; if (!g) return null;
    const r = canvas.getBoundingClientRect();
    return { x: r.left + g.sx(x), y: r.top + (y == null ? (g.top + g.bottom) / 2 : g.sy(ty(y))), top: r.top + g.top, bottom: r.top + g.bottom };
  }
  /** The data under a point of the window (a context menu's place): x and y in the axes' units, the pixels' mapping
      reversed, or null outside the axes. */
  function dataAt(clientX, clientY) {
    const g = st.geom; if (!g) return null;
    const r = canvas.getBoundingClientRect(); const px = clientX - r.left, py = clientY - r.top;
    if (px < g.left || px > g.w - g.right || py < g.top || py > g.bottom) return null;
    return { x: toX(px), y: Tinv[st.yScale](toY(py)) };
  }
  const api = { setSeries, setAxes, fitAll, setMode, setSelection, clearSelection, selection, getView, setView, copyImage, toSVG, draw, pixelOf, dataAt, state: st, setCursor, onViewChange, setFont, setCellCursor };
  if (!canvas.dataset.tour) canvas.dataset.tour = "plot";   // TU5: the tour's handle on every plot
  canvas.__drsPlot = api;
  return api;
}

/** PL13 (F6): the axes of several plots linked: a zoom, a pan or a fit in one is made in the others, for x and, when
    asked, y. Returns a function that unlinks them. */
export function linkAxes(plots, { x = true, y = false } = {}) {
  let busy = false;
  const offs = plots.map((p) => p.onViewChange((v) => {
    if (busy) return;
    busy = true;
    try {
      for (const q of plots) {
        if (q === p) continue;
        const cur = q.getView();
        q.setView({ x0: x ? v.x0 : cur.x0, x1: x ? v.x1 : cur.x1, y0: y ? v.y0 : cur.y0, y1: y ? v.y1 : cur.y1 });
      }
    } finally { busy = false; }
  }));
  return () => offs.forEach((off) => off());
}
