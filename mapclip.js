/* diff-ref-struct · the clip on a map's colour bar and the cell under a pointer (Plots' step 16, I08 step 1 and 3,
   R33; PL13l). Pure screen layout, so that node tests it: a value's place in the range shown, clamped to [0, 1] (a
   value beyond the clip takes the end colour, as plot.js's ramp and surfacegrid.js's height saturate it); the value a
   point of the bar stands for, on the scale the map is drawn on; a clip's bound moved along the bar, kept in order and
   inside the core's range; and the cell of a map under a point, taken whole, never interpolated. No value is made
   here: the bounds come from the user's hand on the bar (or typed) and every cell's value is the core's. */

/** A value's place in the range [from, to], 0 at `from`, 1 at `to`, clamped: below `from` it is 0 (the first colour),
    above `to` it is 1 (the last). A pair not in order places everything in the middle, as plot.js did for an empty
    range. */
export function clipPlace(z, from, to) {
  if (!(to > from)) return 0.5;
  const t = (z - from) / (to - from);
  return t < 0 ? 0 : t > 1 ? 1 : t;
}

/** The value at a fraction `f` of the bar (0 at its foot, the core's `lo`; 1 at its top, `hi`), clamped to the bar. */
export function barValue(f, lo, hi) {
  const t = f < 0 ? 0 : f > 1 ? 1 : f;
  return lo + t * (hi - lo);
}

/** A clip's pair with one bound moved to `v` on the bar: `which` 0 moves `from`, 1 moves `to`. The pair stays inside
    the core's [lo, hi] and in order, the moved bound stopping one thousandth of the range short of the other. A null
    bound stands for the core's own end. Returns [from, to]. */
export function clipMoved(clip, which, v, lo, hi) {
  const gap = (hi - lo) * 1e-3;
  let from = clip && clip.from != null ? clip.from : lo, to = clip && clip.to != null ? clip.to : hi;
  const w = v < lo ? lo : v > hi ? hi : v;
  if (which === 0) from = w > to - gap ? to - gap : w;
  else to = w < from + gap ? from + gap : w;
  return [from, to];
}

/** The drawn cells inside a rectangle of the plane (Plots' step 17, the zoom in the T–2θ plane; PL13m): the columns
    whose cell overlaps [x0, x1] (`xEdges` ascending) as the span j0…j1, and the rows whose band [lo, hi) overlaps
    [y0, y1] as i0…i1; null when no cell is inside. Screen layout: which of the core's cells the view shows, nothing
    about their values (the block's range is the core's, `plots::map_region`). */
export function blockOf(view, xEdges, rows) {
  let j0 = -1, j1 = -1;
  for (let q = 0; q + 1 < xEdges.length; q++) if (xEdges[q + 1] > view.x0 && xEdges[q] < view.x1) { if (j0 < 0) j0 = q; j1 = q; }
  let i0 = -1, i1 = -1;
  rows.forEach((row, i) => { if (row.hi > view.y0 && row.lo < view.y1) { if (i0 < 0) i0 = i; i1 = i; } });
  if (j0 < 0 || (rows.length && i0 < 0)) return null;
  return { j0, j1, i0, i1 };
}

/** The zoom a map keeps across its redraws (step 17): the plot's view while it is zoomed, none (the whole) once the
    plot is fitted again (a double-click), so that a redraw for a clip or a cursor puts the same rectangle back. */
export function keptZoom(view, fit) {
  return fit || !view ? null : { x0: view.x0, x1: view.x1, y0: view.y0, y1: view.y1 };
}

/** "Range from the view" turned on or off (step 17): on, the clip in force is remembered and the view's range (the
    core's block range, `from`–`to`) takes over; off, the remembered clip comes back. Returns { clip, before }. */
export function followSwitch(clip, before, on, range = null) {
  if (on) return { clip: range ? { from: range[0], to: range[1] } : clip, before: clip };
  return { clip: before === undefined ? clip : before, before: null };
}

/** The points a line keeps per pixel column (Plots' step 22, W9 b; PL13r): of the points `i0…i1` (in x order),
    `colOf(i)` their pixel column and `yOf(i)` their pixel height, at most four per column in drawing order, the
    column's first, its lowest, its highest and its last, so that the line drawn is what the full line would show at
    this width; a point with no height (NaN) breaks the line (a null entry). Screen layout: a choice of points, nothing
    summed or averaged; a series with fewer points than three per column is drawn whole by the caller. */
export function pickPerColumn(i0, i1, colOf, yOf) {
  const out = [];
  let col = null, first = -1, lo = -1, hi = -1, last = -1, loY = 0, hiY = 0;
  const flush = () => { if (col == null) return; out.push(first); if (lo !== first) out.push(lo); if (hi !== lo && hi !== first) out.push(hi); if (last !== hi && last !== lo && last !== first) out.push(last); };
  for (let i = i0; i <= i1; i++) {
    const Y = yOf(i);
    if (!isFinite(Y)) { flush(); col = null; out.push(null); continue; }
    const c = colOf(i);
    if (c !== col) { flush(); col = c; first = lo = hi = last = i; loY = hiY = Y; }
    else { if (Y < loY) { lo = i; loY = Y; } if (Y > hiY) { hi = i; hiY = Y; } last = i; }
  }
  flush();
  return out;
}

/** The slices a drawing of `n` series is cut into, `per` series each, for one frame apiece (step 22): [[from, to)…];
    one slice when `per` is not a positive number or `n` fits in it. */
export function sliceBounds(n, per) {
  if (!(per > 0) || n <= per) return n > 0 ? [[0, n]] : [];
  const out = [];
  for (let a = 0; a < n; a += per) out.push([a, Math.min(n, a + per)]);
  return out;
}

/** The column whose edges hold `xv` (`xEdges` ascending, one more than the columns), or −1. */
export function columnAt(xv, xEdges) {
  for (let q = 0; q + 1 < xEdges.length; q++) if (xv >= xEdges[q] && xv < xEdges[q + 1]) return q;
  return -1;
}

/** The cell of a map under a point in data units: the column `j` whose edges hold `xv` and the row `r` whose [lo, hi)
    holds `yv`. The cell itself, never a neighbour or a blend; null when the point lies between no edges or on no row. */
export function cellAt(xv, yv, xEdges, rows) {
  const j = columnAt(xv, xEdges);
  if (j < 0) return null;
  const r = rows.findIndex((row) => yv >= row.lo && yv < row.hi);
  return r < 0 ? null : { j, r };
}
