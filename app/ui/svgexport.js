/* diff-ref-struct · a plot drawn once more into SVG for File ▸ Export ▸ Plot (APP_CONVENTIONS §8, L134). The same code
   that draws a plot on screen draws it into canvas2svg's recorder (vendor/canvas2svg, MIT), so the figure is the plot.
   canvas2svg has no dashed lines and clears to white: setLineDash is recorded here as stroke-dasharray, and a cleared
   area takes the plot's own background. The core turns the SVG into PDF or PNG (src-tauri/src/figure.rs). */

let patched = false;
function patch(C2S) {
  if (patched) return;
  patched = true;
  const P = C2S.prototype;
  P.setLineDash = function (a) { this.__drsDash = a && a.length ? a.join(",") : null; };
  P.getLineDash = function () { return this.__drsDash ? this.__drsDash.split(",").map(Number) : []; };
  const stroke = P.stroke;
  P.stroke = function (...a) {
    stroke.apply(this, a);
    if (this.__drsDash && this.__currentElement && this.__currentElement.setAttribute) this.__currentElement.setAttribute("stroke-dasharray", this.__drsDash);
  };
  const strokeRect = P.strokeRect;
  P.strokeRect = function (...a) {
    strokeRect.apply(this, a);
    const last = this.__closestGroupOrSvg().lastChild;
    if (this.__drsDash && last && last.setAttribute) last.setAttribute("stroke-dasharray", this.__drsDash);
  };
  // canvas2svg moves the current path into the clip without writing its outline (an empty d): every clipped drawing
  // vanished (the Rietveld plot's difference strip). The outline is written first.
  const clip = P.clip;
  P.clip = function (...a) {
    if (this.__currentElement && this.__currentElement.nodeName === "path" && this.__applyCurrentDefaultPath) this.__applyCurrentDefaultPath();
    return clip.apply(this, a);
  };
  const clearRect = P.clearRect;
  P.clearRect = function (x, y, w, h) {
    clearRect.call(this, x, y, w, h);
    const last = this.__closestGroupOrSvg().lastChild;
    if (last && last.setAttribute && this.__drsBg) last.setAttribute("fill", this.__drsBg);
  };
}

/** Runs `fn` with the window's light palette, so that a figure is drawn on white for a document whatever the window's
    appearance; the window's own appearance is put back before anything is painted (one task, no flicker). */
export function inLightPalette(fn) {
  const root = document.documentElement, had = root.getAttribute("data-theme");
  root.setAttribute("data-theme", "light");
  try { return fn(); } finally { if (had == null) root.removeAttribute("data-theme"); else root.setAttribute("data-theme", had); }
}

/** A recorder of `w` × `h` CSS pixels whose cleared areas take `background`. */
export function svgRecorder(w, h, background) {
  const C2S = window.C2S;
  if (!C2S) throw new Error("The plot cannot be exported: the SVG recorder (vendor/canvas2svg) is not loaded.");
  patch(C2S);
  const r = new C2S(Math.max(1, Math.round(w)), Math.max(1, Math.round(h)));
  r.__drsBg = background;
  return r;
}

/** A stand-in for a <canvas> whose context is the recorder, for drawing functions that take a canvas. */
export function canvasFor(recorder, w, h) {
  return { clientWidth: w, clientHeight: h, width: w, height: h, getContext: () => recorder };
}

/** The keys of an HTML legend (`.key` with a `.sw` swatch), as [{ color, label }], those switched off left out. */
export function legendItems(el) {
  if (!el) return [];
  return [...el.querySelectorAll(".key")].filter((k) => !k.classList.contains("off")).map((k) => {
    const sw = k.querySelector(".sw");
    return { color: sw ? getComputedStyle(sw).backgroundColor : "", label: k.textContent.trim() };
  }).filter((k) => k.label);
}

/** Draws a legend row into the recorder from x, at baseline y: a short line in each colour and its label. */
export function drawLegend(ctx, items, x, y, ink, font = "11px -apple-system, sans-serif") {
  ctx.save();
  ctx.font = font; ctx.textAlign = "left"; ctx.lineWidth = 2;
  let at = x;
  for (const it of items) {
    ctx.strokeStyle = it.color || ink; ctx.beginPath(); ctx.moveTo(at, y - 4); ctx.lineTo(at + 14, y - 4); ctx.stroke();
    ctx.fillStyle = ink; ctx.fillText(it.label, at + 18, y);
    at += 18 + ctx.measureText(it.label).width + 16;
  }
  ctx.restore();
}

/** The SVG text of a recorder. */
export const svgText = (r) => r.getSerializedSvg(true);

/** Several SVGs side by side in one figure (Explore's peak-shape view has two panels). */
export function sideBySide(parts) {
  const w = parts.reduce((s, p) => s + p.w, 0), h = Math.max(...parts.map((p) => p.h));
  let x = 0;
  const inner = parts.map((p) => { const s = p.svg.replace(/^<svg /, `<svg x="${x}" y="0" `); x += p.w; return s; }).join("");
  return `<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="${w}" height="${h}" viewBox="0 0 ${w} ${h}">${inner}</svg>`;
}

/** SVGs stacked top to bottom (a legend row over a panel). */
export function stacked(parts) {
  const w = Math.max(...parts.map((p) => p.w)), h = parts.reduce((s, p) => s + p.h, 0);
  let y = 0;
  const inner = parts.map((p) => { const s = p.svg.replace(/^<svg /, `<svg x="0" y="${y}" `); y += p.h; return s; }).join("");
  return `<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="${w}" height="${h}" viewBox="0 0 ${w} ${h}">${inner}</svg>`;
}

/** A picture (a PNG data URL, the 3D view's) as an SVG of w × h. */
export function imageSVG(dataURL, w, h) {
  return `<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="${w}" height="${h}" viewBox="0 0 ${w} ${h}"><image x="0" y="0" width="${w}" height="${h}" xlink:href="${dataURL}" href="${dataURL}"/></svg>`;
}

/** A legend row alone, as an SVG of width w. */
export function legendSVG(items, w, background, ink) {
  if (!items.length) return null;
  const r = svgRecorder(w, 24, background);
  r.clearRect(0, 0, w, 24);
  drawLegend(r, items, 56, 16, ink);
  return { svg: svgText(r), w, h: 24 };
}
