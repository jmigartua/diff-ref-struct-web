/* diff-ref-struct · the one tick rule of every plot (APP_CONVENTIONS §8, L119, L133): as many ticks as the space allows
   without their labels touching, the part common to long numbers taken into the axis label (+ 10.6048, ×10⁻⁵) so that
   each label shows only the digits that change, the true minus, the same digits along an axis. The shared plot
   component (plot.js) and Explore's views (modules/explore) both use it. Pure: tested by node. */

export const MINUS = "−";
export function niceStep(span, target) {
  const raw = span / Math.max(1, target);
  const p = Math.pow(10, Math.floor(Math.log10(raw)));
  const m = raw / p;
  return (m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10) * p;
}
const decimals = (step) => Math.min(10, Math.max(0, -Math.floor(Math.log10(step) + 1e-9)));
const minus = (s) => String(s).replace(/-/g, MINUS);
const SUP = { "-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹" };
export const pow10 = (e) => "10" + String(e).split("").map((c) => SUP[c] || c).join("");

/** The ticks of a linear axis over [lo, hi] on `px` pixels. `measure(label)` gives a label's width in pixels; null for
    a vertical axis, where 26 px per label is the rule. Returns { ticks: [{ value, label }], suffix, step, offset,
    factor }: `suffix` is appended to the axis label (" ×10⁻⁵", " + 10.6048", or both). */
export function axisTicks(lo, hi, px, measure) {
  const span = hi - lo;
  if (!(span > 0) || !isFinite(span)) return { ticks: [], suffix: "", step: 0, offset: 0, factor: 1 };
  const mag = Math.max(Math.abs(lo), Math.abs(hi));
  let offset = 0;
  if (mag / span > 200) { const p = Math.pow(10, Math.floor(Math.log10(span))); offset = Math.floor(lo / p) * p; }
  const m2 = Math.max(Math.abs(lo - offset), Math.abs(hi - offset)) || 1;
  let factor = 1;
  if (m2 >= 1e5 || m2 < 1e-3) factor = Math.pow(10, Math.floor(Math.log10(m2)));
  const a = (lo - offset) / factor, b = (hi - offset) / factor;
  let target = measure ? Math.max(2, Math.floor(px / 70)) : Math.max(2, Math.min(8, Math.floor(px / 48)));
  let step = niceStep(b - a, target), ticks = [];
  for (let guard = 0; guard < 8; guard++) {
    step = niceStep(b - a, target);
    const vals = [];
    for (let t = Math.ceil(a / step - 1e-9) * step; t <= b + step * 1e-9; t += step) vals.push(Math.abs(t) < step * 1e-9 ? 0 : t);
    ticks = vals.map((v) => ({ value: v * factor + offset, label: minus(step >= 1 ? String(Math.round(v)) : v.toFixed(decimals(step))) }));
    if (!measure || ticks.length <= 2) break;
    const maxw = Math.max(...ticks.map((t) => measure(t.label))) + 14;
    if (maxw * ticks.length <= px) break;
    const next = Math.max(2, Math.floor(px / maxw));
    if (next >= target) target = Math.max(2, target - 1); else target = next;
    if (target <= 2 && niceStep(b - a, 2) === step) break;
  }
  const offsetDigits = Math.max(0, -Math.floor(Math.log10(span)));
  const suffix = (factor !== 1 ? ` ×${pow10(Math.round(Math.log10(factor)))}` : "") + (offset ? ` ${offset > 0 ? "+" : MINUS} ${Math.abs(offset).toFixed(Math.min(12, offsetDigits))}` : "");
  return { ticks, suffix, step, offset, factor };
}

/** The ticks of a √ axis: round values in the data's own units (0, 2000, 4000 counts…), placed on the √ scale and
    thinned where they would crowd (`minGap` pixels). `y0`, `y1` are the axis's ends in √ units; each tick's `value`
    is in √ units too, its label in counts. */
export function sqrtTicks(y0, y1, px, minGap = 26) {
  const T = (y) => (y >= 0 ? Math.sqrt(y) : -Math.sqrt(-y)), Tinv = (v) => (v >= 0 ? v * v : -v * v);
  const lin = axisTicks(Tinv(y0), Tinv(y1), px, null);
  const placed = lin.ticks.map((t) => ({ value: T(t.value), label: t.label })).filter((t) => t.value >= y0 - 1e-12 && t.value <= y1 + 1e-12);
  const kept = [];
  for (const t of placed.slice().reverse()) {           // from the top, where √ crowds the ticks
    const at = ((t.value - y0) / (y1 - y0)) * px;
    if (!kept.length || Math.abs(at - kept[kept.length - 1].at) >= minGap) kept.push({ ...t, at });
  }
  return { ticks: kept.reverse().map(({ value, label }) => ({ value, label })), suffix: lin.suffix };
}

/** The ticks of a log₁₀ axis: the decades (10, 100, 1000…), with 2× and 5× between them when fewer than three
    decades show, or when `px` (the axis's length) leaves 70 px a decade or more, so that a value between decades can be
    read (PL30 e). `y0`, `y1` are the axis's ends in log₁₀ units; each tick's `value` is in log₁₀ units, its label the
    number itself (10ⁿ when it is long). */
export function logTicks(y0, y1, px = 0) {
  const out = [];
  const lo = Math.floor(y0), hi = Math.ceil(y1);
  const few = hi - lo < 3 || (px > 0 && px / Math.max(1e-9, y1 - y0) >= 70);
  for (let e = lo; e <= hi; e++) for (const m of few ? [1, 2, 5] : [1]) {
    const v = m * Math.pow(10, e), lv = Math.log10(v);
    if (lv < y0 - 1e-12 || lv > y1 + 1e-12) continue;
    const label = e >= 5 || e <= -4 ? (m === 1 ? pow10(e) : `${m}×${pow10(e)}`) : minus(e >= 0 ? String(v) : v.toFixed(-e));
    out.push({ value: lv, label });
  }
  return { ticks: out, suffix: "" };
}

/** The ticks of a symmetric-log axis (E8 d): linear within ±1, logarithmic beyond, as `sym(t) = t` for |t| ≤ 1 and
    sign(t)·(1 + log₁₀|t|) beyond. `y0`, `y1` are the axis's ends on that scale; the ticks are 0, ±1, ±10, ±100 …,
    each `value` on the scale, its label the number, with the true minus. */
export const symlog = (t) => (Math.abs(t) <= 1 ? t : Math.sign(t) * (1 + Math.log10(Math.abs(t))));
export const symlogInv = (v) => (Math.abs(v) <= 1 ? v : Math.sign(v) * Math.pow(10, Math.abs(v) - 1));
export function symlogTicks(y0, y1) {
  const out = [];
  const big = Math.max(Math.abs(y0), Math.abs(y1));
  const top = big <= 1 ? 0 : Math.ceil(big - 1);
  const vals = [0, 1, -1];
  for (let e = 1; e <= top; e++) vals.push(Math.pow(10, e), -Math.pow(10, e));
  for (const v of vals.sort((a, b) => a - b)) {
    const s = symlog(v);
    if (s < y0 - 1e-12 || s > y1 + 1e-12) continue;
    const a = Math.abs(v), e = Math.round(Math.log10(a || 1));
    out.push({ value: s, label: v === 0 ? "0" : (v < 0 ? MINUS : "") + (e >= 5 ? pow10(e) : String(a)) });
  }
  return { ticks: out, suffix: "" };
}

/** For an axis of n equally spaced items (cycles): the stride k so that labelling every k-th one leaves `minGap` pixels
    between labels of width `width` (px per item = px / (n − 1)). The first and the last are always labelled. */
export function labelStride(n, px, width, minGap = 10) {
  if (n <= 1) return 1;
  const per = px / (n - 1);
  for (const k of [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]) if (per * k >= width + minGap) return k;
  return Math.ceil(n / 2);
}

/** One value at a step's precision (the strip's scale under a plot): exponent form only for very large or very small
    numbers, the true minus. */
export const fmtTick = (v, step) => minus(Math.abs(v) >= 1e5 || (Math.abs(v) < 1e-3 && v !== 0) ? v.toExponential(1) : step >= 1 ? String(Math.round(v)) : v.toFixed(Math.min(6, Math.max(0, -Math.floor(Math.log10(step))))));
