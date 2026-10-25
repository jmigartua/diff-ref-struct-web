/* diff-ref-struct · the shared 3D surface (Plots' Surface view, plots/DESIGN.md §8.1 item 6; PL13k): a map as the core
   gives it drawn as a surface on the shared 3D canvas (canvas3d.js: the renderer on the low-power GPU, the orbit, the
   zoom, the keyboard, the theme, Copy). The layout is surfacegrid.js's; here three.js objects are made of it: the
   surface coloured on the map's own ramps (plot.js's viridis, RdBu for a difference), each member's row drawn over it
   as a line, the box's three axes with their ticks and titles, a mark at the cell under the pointer and the readout's
   text for it, which the module gives (`cellText`, as the 2D map's). Nothing animates by itself (Reduce motion). It
   draws exactly the core's values; a WebGL canvas is an image to VoiceOver, so `label` names it and the module's
   members table holds the numbers. */
import { createCanvas3D } from "./canvas3d.js";
import { rampColour } from "./plot.js";
import { surfaceGrid, nearestCell } from "./surfacegrid.js";

const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
// the starting view: from the front, above and to the left, so that the coldest row is nearest and the peaks rise
const TILT = { x: 0.62, y: -0.55 };

/** container: the element to fill; opts: { label (VoiceOver's words), readout (an element for the cell's text) }. */
export function createSurface(container, opts = {}) {
  const c3 = createCanvas3D(container, { label: opts.label || "The map as a surface in three dimensions" });
  const { THREE } = c3;
  const live = new AbortController(), on = { signal: live.signal };
  const st = { data: null, grid: null, mesh: null, mark: null, holder: null, stats: null };

  /** A text as a sprite that faces the eye, `px` its height in the label canvas, `size` its height in the scene. */
  function label(text, colour, size, anchor = "center") {
    const font = getComputedStyle(document.body).fontFamily || "-apple-system, sans-serif";
    const px = 48, cv = document.createElement("canvas"), cx = cv.getContext("2d");
    cx.font = `${px}px ${font}`;
    const w = Math.ceil(cx.measureText(text).width) + 8;
    cv.width = w; cv.height = Math.ceil(px * 1.3);
    cx.font = `${px}px ${font}`; cx.fillStyle = colour; cx.textBaseline = "middle";
    cx.textAlign = anchor === "left" ? "left" : anchor === "right" ? "right" : "center";
    cx.fillText(text, anchor === "left" ? 4 : anchor === "right" ? w - 4 : w / 2, cv.height / 2);
    const tex = new THREE.CanvasTexture(cv); tex.colorSpace = THREE.SRGBColorSpace; tex.anisotropy = 4;
    const sp = new THREE.Sprite(new THREE.SpriteMaterial({ map: tex, depthTest: false, transparent: true }));
    sp.scale.set((size * cv.width) / cv.height, size, 1);
    sp.center.set(anchor === "left" ? 0 : anchor === "right" ? 1 : 0.5, 0.5);
    sp.renderOrder = 2;
    return sp;
  }
  // a label put at a point (three.js's `position` is not writable: it is set, never assigned)
  const placed = (sp, x, y, z) => { sp.position.set(x, y, z); return sp; };
  const segments = (pts, colour, opacity = 1) => {
    const g = new THREE.BufferGeometry(); g.setAttribute("position", new THREE.Float32BufferAttribute(pts, 3));
    return new THREE.LineSegments(g, new THREE.LineBasicMaterial({ color: colour, transparent: opacity < 1, opacity }));
  };

  /** Everything made again from the data, in the theme's colours now (a figure is drawn in the light ones). */
  function build() {
    if (st.holder) st.holder.traverse((o) => { if (o.material && o.material.map) o.material.map.dispose(); });   // the labels' textures
    c3.clear(); st.mesh = null; st.mark = null; st.holder = null;
    const d = st.data; if (!d) return;
    const g = (st.grid = surfaceGrid(d.map, d.scale, undefined, d.clip || null));   // PL13l: the height on the clip's range when one is set
    const { box, at } = g;
    const ink = css("--ink") || "#222222", muted = css("--muted") || "#777777";
    const holder = (st.holder = new THREE.Group()); holder.rotation.set(TILT.x, TILT.y, 0, "XYZ");
    c3.group.add(holder);
    // the surface: one vertex per value the core gave, its colour its place in the core's range on the map's ramp
    const geo = new THREE.BufferGeometry();
    geo.setAttribute("position", new THREE.BufferAttribute(g.positions, 3));
    const col = new Float32Array(g.fraction.length * 3), c = new THREE.Color();
    for (let k = 0; k < g.fraction.length; k++) { c.set(isNaN(g.fraction[k]) ? muted : rampColour(d.ramp, g.fraction[k])); col[3 * k] = c.r; col[3 * k + 1] = c.g; col[3 * k + 2] = c.b; }
    geo.setAttribute("color", new THREE.BufferAttribute(col, 3));
    geo.setIndex(new THREE.BufferAttribute(g.index, 1));
    geo.computeVertexNormals();
    st.mesh = new THREE.Mesh(geo, new THREE.MeshLambertMaterial({ vertexColors: true, side: THREE.DoubleSide, polygonOffset: true, polygonOffsetFactor: 1, polygonOffsetUnits: 1 }));
    holder.add(st.mesh);
    // each member's row over it, as a line where it has values: the patterns themselves, at their temperatures
    const rows = [];
    for (let i = 0; i < g.nr; i++) for (let j = 0; j + 1 < g.nc; j++) {
      const a = i * g.nc + j, b = a + 1;
      if (isNaN(g.fraction[a]) || isNaN(g.fraction[b])) continue;
      rows.push(g.positions[3 * a], g.positions[3 * a + 1], g.positions[3 * a + 2], g.positions[3 * b], g.positions[3 * b + 1], g.positions[3 * b + 2]);
    }
    holder.add(segments(rows, ink, 0.35));
    // the box: the base's outline, and the three axes along the edges the eye sees first
    const W = box.w / 2, D = box.d / 2, H = box.h;
    holder.add(segments([-W, 0, D, W, 0, D, W, 0, D, W, 0, -D, W, 0, -D, -W, 0, -D, -W, 0, -D, -W, 0, D], muted, 0.6));
    holder.add(segments([-W, 0, -D, -W, H, -D, -W, H, -D, W, H, -D, -W, H, -D, -W, H, D], muted, 0.25));
    const tk = 0.05, ts = 0.085, tickPts = [];
    for (const t of g.ticks.x.ticks) { const X = at.x(t.value); tickPts.push(X, 0, D, X, 0, D + tk); holder.add(placed(label(t.label, muted, ts), X, 0, D + 2.4 * tk)); }
    for (const t of g.ticks.v.ticks) { const Z = at.v(t.value); tickPts.push(W, 0, Z, W + tk, 0, Z); holder.add(placed(label(t.label, muted, ts, "left"), W + 1.6 * tk, 0, Z)); }
    for (const t of g.ticks.z.ticks) { const Y = at.z(t.value); tickPts.push(-W, Y, -D, -W - tk, Y, -D); holder.add(placed(label(t.label, muted, ts, "right"), -W - 1.6 * tk, Y, -D)); }
    holder.add(segments(tickPts, ink, 0.8));
    const L = d.labels || {};
    holder.add(placed(label((L.x || "x") + (g.ticks.x.suffix || ""), ink, 0.1), 0, 0, D + 6 * tk));
    holder.add(placed(label((L.v || "") + (g.ticks.v.suffix || ""), ink, 0.1, "left"), W + 1.6 * tk, 0, -D - 3 * tk));
    holder.add(placed(label((L.z || "") + (g.ticks.z.suffix || ""), ink, 0.1, "left"), -W, H + 3 * tk, -D));
    // the mark at the cell under the pointer: a line from the base up to the value
    st.mark = segments([0, 0, 0, 0, 1, 0], css("--accent") || "#0a84ff", 1); st.mark.visible = false; st.mark.renderOrder = 3; st.mark.material.depthTest = false;
    holder.add(st.mark);
    st.stats = { rows: g.nr, columns: g.nc, vertices: g.fraction.length, triangles: g.index.length / 3, ticks: [g.ticks.x.ticks.length, g.ticks.v.ticks.length, g.ticks.z.ticks.length], range: g.range };
    c3.frame([0, H / 2, 0], Math.hypot(W, D, H / 2) * 1.05);
  }

  // the cell under the pointer: the surface's triangle hit, the nearest of its three vertices, the module's words for it
  const ray = new THREE.Raycaster(), ndc = new THREE.Vector2();
  function cellAt(clientX, clientY) {
    if (!st.mesh) return null;
    const r = c3.canvas.getBoundingClientRect();
    ndc.set(((clientX - r.left) / r.width) * 2 - 1, -((clientY - r.top) / r.height) * 2 + 1);
    ray.setFromCamera(ndc, c3.camera);
    const hit = ray.intersectObject(st.mesh, false)[0];
    if (!hit || !hit.face) return null;
    const p = st.mesh.worldToLocal(hit.point.clone());
    return nearestCell(st.grid, [hit.face.a, hit.face.b, hit.face.c], [p.x, p.y, p.z]);
  }
  function showCell(cell) {
    if (!st.mark) return;
    if (!cell) { st.mark.visible = false; if (opts.readout) opts.readout.textContent = ""; c3.render(); return; }
    const k = cell.i * st.grid.nc + cell.j, P = st.grid.positions, a = st.mark.geometry.attributes.position;
    a.setXYZ(0, P[3 * k], 0, P[3 * k + 2]); a.setXYZ(1, P[3 * k], P[3 * k + 1], P[3 * k + 2]); a.needsUpdate = true;
    st.mark.visible = true;
    if (opts.readout && st.data && st.data.cellText) opts.readout.textContent = st.data.cellText(cell.j, cell.i);
    c3.render();
  }
  c3.canvas.addEventListener("mousemove", (e) => { if (e.buttons) return; showCell(cellAt(e.clientX, e.clientY)); }, on);
  c3.canvas.addEventListener("mouseleave", () => showCell(null), on);
  c3.onTheme(build);

  return {
    canvas: c3.canvas,
    /** data: { map (the core's: x, rows, z, range), scale ("linear" | "sqrt" | "log"), ramp ("sequential" |
        "diverging"), labels: { x, v, z }, cellText(j, i), clip ({ from, to } | null, PL13l: the range the height spans) }. */
    set(data) { st.data = data; build(); },
    /** The view now as a PNG (a data URL) at the canvas's size in CSS pixels, drawn in the theme's colours of the
        moment (File ▸ Export ▸ Plot draws in the light ones). */
    image() { build(); const w = c3.canvas.clientWidth, h = c3.canvas.clientHeight; c3.render(); return { url: c3.canvas.toDataURL("image/png"), w, h }; },
    /** What is drawn, counted (the self-check's). */
    stats: () => st.stats,
    /** The readout of a cell, as the pointer gives it (the self-check's). */
    showCell,
    reset: c3.reset, copyImage: c3.copyImage,
    dispose() { live.abort(); c3.dispose(); },
  };
}
