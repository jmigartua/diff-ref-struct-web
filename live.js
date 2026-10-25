/* diff-ref-struct Web · live mode (I10 step 4; web/DESIGN.md): a folder opened by its handle (Chromium) read again
   every few seconds. What changed between two readings, by each file's path in the folder and its stamp (size and
   date): the files that appeared, vanished or changed. Pure, so that node tests it. */

/** The difference between two listings (Map path → stamp): { added, vanished, changed }, each by name. */
export function diffListing(prev, next) {
  const by = (a, b) => a.localeCompare(b, undefined, { numeric: true, sensitivity: "base" });
  const added = [], vanished = [], changed = [];
  for (const [p, s] of next) { if (!prev.has(p)) added.push(p); else if (prev.get(p) !== s) changed.push(p); }
  for (const p of prev.keys()) if (!next.has(p)) vanished.push(p);
  return { added: added.sort(by), vanished: vanished.sort(by), changed: changed.sort(by) };
}
export const nothing = (d) => !d.added.length && !d.vanished.length && !d.changed.length;

/** The intervals offered, in s, and the one taken first. */
export const INTERVALS = [2, 5, 10, 30, 60];
export const DEFAULT_INTERVAL = 10;

/** The members to add to a set opened from a proposal: the files the core proposes now for the same kind and class,
    less those the set holds (paths). `proposals` as plots_propose gives them; `match` { kind, class }. */
export function growth(proposals, match, held) {
  const p = (proposals || []).find((q) => q.kind === match.kind && (q.class ?? null) === (match.class ?? null));
  if (!p) return [];
  return p.paths.filter((f) => !held.has(f));
}

/** What live mode says of a reading, in one line (null when nothing changed). */
export function liveLine(d, grown = 0, at = new Date()) {
  if (nothing(d) && !grown) return null;
  const t = at.toTimeString().slice(0, 8);
  const part = (n, w) => (n ? `${n} ${w}` : "");
  const parts = [part(d.added.length, d.added.length === 1 ? "new file" : "new files"), part(d.changed.length, "changed"), part(d.vanished.length, "vanished (still in the sets opened: nothing is removed from a set by itself)"), grown ? `${grown} scan${grown === 1 ? "" : "s"} added to the live set` : ""].filter(Boolean);
  return `${t}: ${parts.join(", ")}.`;
}
