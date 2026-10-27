/* diff-ref-struct Web · the inspector's width (web/DESIGN.md §4): the separator dragged, clamped, kept in
   localStorage; `<` shrinks the inspector to a rail and back. The arithmetic here is the layout's (pixels), not the
   data's. Pure apart from `keep`/`kept`, which are given the storage. */

export const RAIL = 36;                    // the rail's width, px
export const LIMITS = { min: 220, max: 640, views: 360, start: 320 };   // the inspector's least and most; the views' least

/** The inspector's width for a wanted `w` in a window `total` px wide: within [min, max], and never so wide that the
    views get less than `views` px (then the least of the two wins: min). A width that is not a number is `start`. */
export function clampWidth(w, total = Infinity, lim = LIMITS) {
  if (!Number.isFinite(w)) w = lim.start;
  const room = Number.isFinite(total) ? total - lim.views : lim.max;
  return Math.round(Math.max(lim.min, Math.min(w, lim.max, room)));
}

/** The width after the separator is dragged from `x0` to `x` (the pointer's x), the inspector at the left starting
    `w0` wide. */
export const dragged = (w0, x0, x, total, lim = LIMITS) => clampWidth(w0 + (x - x0), total, lim);

/** The width after an arrow key on the focused separator: 16 px a press, 64 with ⇧. */
export const keyed = (w0, key, shift, total, lim = LIMITS) => (key === "ArrowLeft" || key === "ArrowRight" ? clampWidth(w0 + (key === "ArrowRight" ? 1 : -1) * (shift ? 64 : 16), total, lim) : w0);

const KEY = "drs-web.inspector";
/** What the storage keeps, { width, rail }, or the defaults when it keeps nothing readable (a private window). */
export function kept(storage) {
  try { const v = JSON.parse(storage.getItem(KEY) || "null"); if (v && typeof v === "object") return { width: Number.isFinite(v.width) ? v.width : LIMITS.start, rail: !!v.rail }; } catch (e) { /* none kept */ }
  return { width: LIMITS.start, rail: false };
}
export function keep(storage, v) { try { storage.setItem(KEY, JSON.stringify({ width: v.width, rail: !!v.rail })); } catch (e) { /* storage refused */ } }
