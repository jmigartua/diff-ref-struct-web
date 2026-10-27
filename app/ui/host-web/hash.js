/* diff-ref-struct Web · the view in the address (I10 step 4; web/DESIGN.md): the set shown and how (map or difference,
   the colour scale, the zoom, the clip, the region) in the URL's hash, so that a link or a reload brings the view back
   when the same files are open again. Only names and the numbers the view already holds: no file's content, nothing
   computed. Pure, so that node tests it. */

const num = (v) => { const x = Number(v); return Number.isFinite(x) ? x : null; };
const nums = (s, n) => { const a = String(s || "").split(",").map(num); return a.length === n && a.every((x) => x !== null) ? a : null; };
const VIEWS = ["overlay", "stack", "projection", "map", "difference", "surface", "series"];
const SCALES = ["linear", "sqrt", "log"];

/** The hash of a view state: { folder, set: { by: "files" | "drv" | "buf" | "proposal", files?, path?, kind?, class? },
    view, scale, zoom: { x0, x1, y0, y1 }, clip: { from, to }, region }. Absent parts are left out. */
export function encode(st) {
  if (!st || !st.set) return "";
  const q = new URLSearchParams();
  if (st.folder) q.set("folder", st.folder);
  const s = st.set;
  if (s.by === "files") { q.set("set", "files"); for (const f of s.files) q.append("file", f); }   // one parameter per file: a name may hold any character
  else if (s.by === "drv" || s.by === "buf") q.set("set", `${s.by}:${s.path}`);
  else if (s.by === "proposal") q.set("set", `proposal:${s.kind}:${s.class ?? ""}`);
  else return "";
  if (st.view && VIEWS.includes(st.view)) q.set("view", st.view);
  if (st.scale && SCALES.includes(st.scale)) q.set("scale", st.scale);
  if (st.zoom) q.set("zoom", [st.zoom.x0, st.zoom.x1, st.zoom.y0, st.zoom.y1].map(String).join(","));
  if (st.clip) q.set("clip", `${st.clip.from},${st.clip.to}`);
  if (st.region) q.set("region", st.region);
  return "#" + q.toString();
}

/** The view state of a hash, or null when it names no set; a part that does not read is left out. */
export function decode(hash) {
  const q = new URLSearchParams(String(hash || "").replace(/^#/, ""));
  const raw = q.get("set"); if (!raw) return null;
  const i = raw.indexOf(":");
  const by = i < 0 ? raw : raw.slice(0, i), rest = i < 0 ? "" : raw.slice(i + 1);
  let set = null;
  if (by === "files") { const files = q.getAll("file").filter(Boolean); if (files.length) set = { by, files }; }
  else if ((by === "drv" || by === "buf") && rest) set = { by, path: rest };
  else if (by === "proposal") { const [kind, cls] = rest.split(":"); if (kind) set = { by, kind, class: cls === "" || cls == null ? null : num(cls) }; }
  if (!set) return null;
  const st = { folder: q.get("folder") || null, set, view: null, scale: null, zoom: null, clip: null, region: null };
  if (VIEWS.includes(q.get("view"))) st.view = q.get("view");
  if (SCALES.includes(q.get("scale"))) st.scale = q.get("scale");
  const z = nums(q.get("zoom"), 4); if (z && z[0] < z[1] && z[2] < z[3]) st.zoom = { x0: z[0], x1: z[1], y0: z[2], y1: z[3] };
  const c = nums(q.get("clip"), 2); if (c && c[0] < c[1]) st.clip = { from: c[0], to: c[1] };
  if (q.get("region")) st.region = q.get("region");
  return st;
}

/** Whether the files a hash's set needs are among the folder's (paths in the folder): which are missing. */
export function missingFor(set, has) {
  if (!set) return [];
  if (set.by === "files") return set.files.filter((f) => !has(f));
  if (set.by === "drv" || set.by === "buf") return has(set.path) ? [] : [set.path];
  return [];
}
