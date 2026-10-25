/* diff-ref-struct · the shared 3D canvas (Explore's E8, explore/DESIGN.md §4; app 0.38.0, §11): three.js (vendored by
   tools/vendor.js, L112) set up once for every module's 3D view. What every such view needs and none should build again:
   - the renderer on the integrated GPU (`low-power`, L150), its size following the container, nothing animating by
     itself (Reduce motion: the view changes only when asked);
   - an orthographic camera over a group the module fills (`group`), framed on a centre and a radius it gives;
   - orbit by dragging, zoom by the wheel, a double-click back to the start; and from the keyboard (the canvas takes
     the focus, a visible ring): ← → ↑ ↓ orbit, + − zoom, 0 back to the start (§6, §10);
   - the theme's colours (`colour("--ink")`) read again when the appearance changes, `onTheme(fn)` to redraw then;
   - the object under the pointer (`pick`), Copy as Image through the pasteboard as the plots do, `dispose()`.
   It draws exactly what the module puts in `group`, which comes from the core (atoms, bonds, faces, arrows): no
   physics here. A WebGL canvas is an image to VoiceOver: `label` names it, and the module gives the same data as a
   table beside it. */
import * as THREE from "./vendor/three/three.module.js";

const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

/** container: the element to fill. opts: { label (VoiceOver's words), onPick(object | null, event) (a click),
    onHover(object | null) }. */
export function createCanvas3D(container, opts = {}) {
  const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: false, powerPreference: "low-power", preserveDrawingBuffer: false });
  renderer.setPixelRatio(window.devicePixelRatio || 1);
  const canvas = renderer.domElement;
  Object.assign(canvas.style, { width: "100%", height: "100%", display: "block" });
  canvas.tabIndex = 0;
  canvas.setAttribute("role", "img");
  canvas.setAttribute("aria-label", opts.label || "A view in three dimensions");
  canvas.dataset.tip = "Drag to turn, the wheel to zoom, a double-click to start again; with the keyboard: arrows turn, + and − zoom, 0 starts again.";
  container.appendChild(canvas);
  const live = new AbortController(), on = { signal: live.signal };
  const scene = new THREE.Scene();
  const camera = new THREE.OrthographicCamera(-10, 10, 10, -10, -500, 500);
  scene.add(new THREE.AmbientLight(0xffffff, 1.6));
  const key = new THREE.DirectionalLight(0xffffff, 2.2); key.position.set(3, 5, 8); scene.add(key);
  const group = new THREE.Group(); scene.add(group);
  const st = { zoom: 1, rot: new THREE.Quaternion(), centre: new THREE.Vector3(), radius: 8, themed: [] };

  function render() {
    const w = container.clientWidth || 300, h = container.clientHeight || 200;
    renderer.setSize(w, h, false);
    scene.background = new THREE.Color(css("--bg") || "#ffffff");
    const aspect = w / h, half = (st.radius * 1.15) / st.zoom;
    camera.left = -half * aspect; camera.right = half * aspect; camera.top = half; camera.bottom = -half;
    camera.updateProjectionMatrix();
    group.quaternion.copy(st.rot);
    group.position.copy(st.centre.clone().applyQuaternion(st.rot).multiplyScalar(-1));
    camera.position.set(0, 0, 10); camera.lookAt(0, 0, 0);
    renderer.render(scene, camera);
  }
  /** The view framed on `centre` ([x, y, z]) with `radius`, as the module measures its scene (the core's box). */
  function frame(centre, radius) { st.centre = new THREE.Vector3(...centre); st.radius = radius > 0 ? radius : 1; render(); }
  function turn(dx, dy) { st.rot.premultiply(new THREE.Quaternion().setFromEuler(new THREE.Euler(dy, dx, 0, "XYZ"))); render(); }
  function reset() { st.zoom = 1; st.rot = new THREE.Quaternion(); render(); }
  /** Along a direction: "a", "b", "c" (the module's axes as x, y, z), or the start. */
  function look(axis) { const e = { a: new THREE.Euler(0, Math.PI / 2, 0), b: new THREE.Euler(-Math.PI / 2, 0, 0), c: new THREE.Euler(0, 0, 0) }[axis]; st.rot = e ? new THREE.Quaternion().setFromEuler(e) : new THREE.Quaternion(); render(); }

  // the mouse and the trackpad
  let drag = null;
  canvas.addEventListener("mousedown", (e) => { drag = { x: e.clientX, y: e.clientY, moved: false }; canvas.focus(); }, on);
  window.addEventListener("mousemove", (e) => {
    if (!drag) return;
    const dx = (e.clientX - drag.x) / 150, dy = (e.clientY - drag.y) / 150;
    if (Math.abs(e.clientX - drag.x) + Math.abs(e.clientY - drag.y) > 2) drag.moved = true;
    drag.x = e.clientX; drag.y = e.clientY; turn(dx, dy);
  }, on);
  window.addEventListener("mouseup", (e) => { if (drag && !drag.moved && opts.onPick) opts.onPick(pick(e.clientX, e.clientY), e); drag = null; }, on);
  canvas.addEventListener("wheel", (e) => { e.preventDefault(); st.zoom = Math.max(0.2, Math.min(8, st.zoom * Math.exp(-e.deltaY / 400))); render(); }, { passive: false, signal: live.signal });
  canvas.addEventListener("dblclick", reset, on);
  canvas.addEventListener("keydown", (e) => {
    const step = Math.PI / 24;
    const k = { ArrowLeft: () => turn(-step, 0), ArrowRight: () => turn(step, 0), ArrowUp: () => turn(0, -step), ArrowDown: () => turn(0, step), "+": () => { st.zoom = Math.min(8, st.zoom * 1.2); render(); }, "=": () => { st.zoom = Math.min(8, st.zoom * 1.2); render(); }, "-": () => { st.zoom = Math.max(0.2, st.zoom / 1.2); render(); }, 0: reset }[e.key];
    if (k) { e.preventDefault(); k(); }
  }, on);
  // the object under the pointer
  const ray = new THREE.Raycaster(), ndc = new THREE.Vector2();
  function pick(clientX, clientY, objects = group.children) {
    const r = canvas.getBoundingClientRect();
    ndc.set(((clientX - r.left) / r.width) * 2 - 1, -((clientY - r.top) / r.height) * 2 + 1);
    ray.setFromCamera(ndc, camera);
    const hit = ray.intersectObjects(objects, true)[0];
    return hit ? hit.object : null;
  }
  if (opts.onHover) canvas.addEventListener("mousemove", (e) => { if (!drag) opts.onHover(pick(e.clientX, e.clientY)); }, on);
  // the theme: the background at every drawing, the module's own colours through onTheme
  const themed = () => { for (const fn of st.themed) { try { fn(); } catch (e) { /* the module's */ } } render(); };
  window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", themed, on);
  const mo = new MutationObserver(themed); mo.observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme", "data-accent"] });
  const ro = new ResizeObserver(() => render()); ro.observe(container);

  /** The view as a PNG on the pasteboard (as the plots: AppKit's, so that a menu chosen while the window is behind
      still copies), rendered once more and read in the same step. */
  async function copyImage() {
    render();
    const data = canvas.toDataURL("image/png");
    // W6 (PL13q): through the host (the page's, else Tauri's), the browser's clipboard when there is none
    const h = window.__DRS_HOST__ || (window.__TAURI__ ? { invoke: (c, a) => window.__TAURI__.core.invoke(c, a) } : null);
    if (h && typeof h.invoke === "function") await h.invoke("clipboard_png", { png: data.split(",")[1] });
    else await navigator.clipboard.write([new ClipboardItem({ "image/png": await (await fetch(data)).blob() })]);
    return true;
  }
  return {
    THREE, scene, group, camera, renderer, canvas,
    render, frame, look, reset, pick, copyImage,
    /** A theme colour (`--ink`, `--phase-0` …) as three.js takes it. */
    colour: (name, fallback = "#888888") => new THREE.Color(css(name) || fallback),
    /** `fn()` when the appearance changes, before the view is drawn again (the module re-colours its objects). */
    onTheme: (fn) => { st.themed.push(fn); },
    /** Everything the module put in `group` taken away, its geometries and materials freed. */
    clear() { for (const o of [...group.children]) { group.remove(o); o.traverse((q) => { if (q.geometry) q.geometry.dispose(); if (q.material) (Array.isArray(q.material) ? q.material : [q.material]).forEach((m) => m.dispose()); }); } render(); },
    dispose() { live.abort(); ro.disconnect(); mo.disconnect(); renderer.dispose(); renderer.forceContextLoss(); if (canvas.parentNode) canvas.parentNode.removeChild(canvas); },
  };
}
