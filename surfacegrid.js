/* diff-ref-struct · the layout of a 3D surface (Plots' Surface view, plots/DESIGN.md §8.1 item 6; PL13k): a map as the
   core gives it (`plots::map`: columns at x, one row per member at its value of the variable, each value on the scale
   the core chose, the core's colour range) placed in a box of the scene, as plot.js places a map's cells in pixels.
   Screen layout only: where each value goes and its place in the colour range, the triangles between four values the
   core gave (a cell with no value leaves a hole), and the axes' ticks by the one tick rule (ticks.js). No value is made
   here. Pure, so that node tests it; surface3d.js draws it with three.js. */
import { axisTicks, sqrtTicks, logTicks } from "./ticks.js";

/** The box: x across (scene x, −w/2 … w/2), the variable in depth (scene z, +d/2 at the first row, nearest the eye,
    to −d/2 at the last), the value up (scene y, 0 at the range's low end to h at its high end). */
export const BOX = { w: 3, d: 2, h: 1.2 };

const place = (v, lo, hi) => (hi > lo ? (v - lo) / (hi - lo) : 0.5);
const clamp01 = (t) => (t < 0 ? 0 : t > 1 ? 1 : t);

/** map: { x: [nc], rows: [{ value }] (nr, in the core's order), z: [nr][nc] (null for no value), range: [lo, hi] },
    scale: "linear" | "sqrt" | "log" (the scale z is on, for the value axis's ticks). Returns { nr, nc, positions
    (Float32Array, 3 per vertex, vertex i·nc + j), fraction (Float32Array, the vertex's place in the range, 0–1, NaN
    with no value), index (Uint32Array, 3 per triangle), at: { x(v), v(v), z(v) } (a value's scene coordinate on each
    axis), ticks: { x, v, z } (each [{ value, label }], with `suffix`), box, range (the [lo, hi] the height spans) }.
    PL13l (Plots, step 16): `clip: { from, to }` (either null for the core's end) is the range the height and the colour
    span instead of the core's, as the map's colour bar clips them: a value beyond it sits at the floor or the top. */
export function surfaceGrid(map, scale = "linear", box = BOX, clip = null) {
  const nc = map.x.length, nr = map.rows.length;
  const x0 = map.x[0], x1 = map.x[nc - 1];
  const v0 = nr ? map.rows[0].value : 0, v1 = nr ? map.rows[nr - 1].value : 0;
  const lo = clip && clip.from != null ? clip.from : map.range[0], hi = clip && clip.to != null ? clip.to : map.range[1];
  const at = {
    x: (v) => (place(v, x0, x1) - 0.5) * box.w,
    v: (v) => (0.5 - place(v, v0, v1)) * box.d,
    z: (v) => clamp01(place(v, lo, hi)) * box.h,
  };
  const positions = new Float32Array(nr * nc * 3), fraction = new Float32Array(nr * nc);
  const xs = map.x.map(at.x);
  for (let i = 0; i < nr; i++) {
    const zr = map.z[i], vz = at.v(map.rows[i].value);
    for (let j = 0; j < nc; j++) {
      const k = i * nc + j, zv = zr[j];
      positions[3 * k] = xs[j]; positions[3 * k + 1] = zv == null ? 0 : at.z(zv); positions[3 * k + 2] = vz;
      fraction[k] = zv == null ? NaN : clamp01(place(zv, lo, hi));
    }
  }
  const tri = [];
  for (let i = 0; i + 1 < nr; i++) for (let j = 0; j + 1 < nc; j++) {
    const a = i * nc + j, b = a + 1, c = a + nc, d = c + 1;
    if (isNaN(fraction[a]) || isNaN(fraction[b]) || isNaN(fraction[c]) || isNaN(fraction[d])) continue;
    tri.push(a, c, b, b, c, d);
  }
  const zt = scale === "log" ? logTicks(lo, hi) : scale === "sqrt" ? sqrtTicks(lo, hi, 240) : axisTicks(lo, hi, 240, null);
  return {
    nr, nc, positions, fraction, index: Uint32Array.from(tri), at, box, range: [lo, hi],
    ticks: { x: axisTicks(x0, x1, 360, null), v: nr > 1 ? axisTicks(Math.min(v0, v1), Math.max(v0, v1), 240, null) : { ticks: [{ value: v0, label: String(v0) }], suffix: "" }, z: zt },
  };
}

/** The cell (row i, column j) of the vertex nearest a point of the surface: of a triangle's three vertices (a hit's
    `face.a`, `b`, `c`), the one nearest `p` ([x, y, z] in the scene), as the 2D map's readout takes the cell under the
    pointer. */
export function nearestCell(g, verts, p) {
  let best = -1, bd = Infinity;
  for (const k of verts) {
    const dx = g.positions[3 * k] - p[0], dy = g.positions[3 * k + 1] - p[1], dz = g.positions[3 * k + 2] - p[2];
    const dd = dx * dx + dy * dy + dz * dz;
    if (dd < bd) { bd = dd; best = k; }
  }
  return best < 0 ? null : { i: Math.floor(best / g.nc), j: best % g.nc };
}
