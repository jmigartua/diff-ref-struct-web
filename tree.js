/* diff-ref-struct Web · the folder's tree (web/DESIGN.md §3): the files read in the browser as a tree of folders and
   files, by their names only (what a file is, the core says). Pure, so that node tests it. */

/** Names in the order a person reads them: "scan_2" before "scan_10", case aside. */
export const byName = (a, b) => a.localeCompare(b, undefined, { numeric: true, sensitivity: "base" });

/** The tree of a file list [{ path: "a/b/c.xy", size }], paths relative to the folder picked: each node
    { name, path, dir, children (folders first, then files, each by name), files (how many files under it), size }.
    A path's empty parts and "." are dropped; the same path twice is kept once (the last). */
export function buildTree(list, rootName = "") {
  const root = { name: rootName, path: "", dir: true, children: [], files: 0, size: 0 };
  const dirs = new Map([["", root]]);
  const seen = new Map();
  for (const f of list) {
    const parts = String(f.path).split("/").filter((p) => p && p !== ".");
    if (!parts.length) continue;
    seen.set(parts.join("/"), { parts, size: f.size || 0 });
  }
  for (const { parts, size } of seen.values()) {
    let dir = root;
    for (let k = 0; k < parts.length - 1; k++) {
      const p = parts.slice(0, k + 1).join("/");
      let d = dirs.get(p);
      if (!d) { d = { name: parts[k], path: p, dir: true, children: [], files: 0, size: 0 }; dirs.set(p, d); dir.children.push(d); }
      dir = d;
    }
    dir.children.push({ name: parts[parts.length - 1], path: parts.join("/"), dir: false, size });
  }
  (function finish(d) {
    d.children.sort((a, b) => (a.dir !== b.dir ? (a.dir ? -1 : 1) : byName(a.name, b.name)));
    d.files = 0; d.size = 0;
    for (const c of d.children) { if (c.dir) finish(c); d.files += c.dir ? c.files : 1; d.size += c.size; }
  })(root);
  return root;
}

/** The files grouped by what the core says they are (`describe`, `kindOf(path)` → { kind, id, name } or null while it
    is being asked): one group per kind's name, the data formats first, then the project's files, then what no reader
    recognises; each file under its group by its path in the folder. A tree of the same shape as buildTree's. */
const KIND_ORDER = { data: 0, project: 1, unknown: 2 };
export function groupByKind(list, kindOf, rootName = "") {
  const root = { name: rootName, path: "", dir: true, children: [], files: 0, size: 0 };
  const groups = new Map();
  for (const f of list) {
    const k = kindOf(f.path);
    const label = k ? k.name || (k.kind === "unknown" ? "not recognised by a reader" : k.id) : "being described by the core…";
    const key = (k ? k.kind : "pending") + ":" + label;
    let g = groups.get(key);
    if (!g) { g = { name: label, path: "kind:" + key, dir: true, kind: k ? k.kind : "pending", id: k ? k.id : null, children: [], files: 0, size: 0 }; groups.set(key, g); }
    g.children.push({ name: f.path, path: f.path, dir: false, size: f.size || 0 });
    g.files += 1; g.size += f.size || 0;
  }
  root.children = [...groups.values()].sort((a, b) => (KIND_ORDER[a.kind] ?? 3) - (KIND_ORDER[b.kind] ?? 3) || byName(a.name, b.name));
  for (const g of root.children) { g.children.sort((a, b) => byName(a.name, b.name)); root.files += g.files; root.size += g.size; }
  return root;
}

/** Every file under a node, in the tree's order (a folder's Open as Plot Set takes them). */
export function filesUnder(node) {
  if (!node.dir) return [node.path];
  return node.children.flatMap(filesUnder);
}

/** The rows shown: the tree walked in order, a folder's children only when it is open (`open`: a Set of folder
    paths); each row { node, depth }. The root itself is not a row. */
export function visibleRows(root, open) {
  const out = [];
  (function walk(d, depth) {
    for (const c of d.children) { out.push({ node: c, depth }); if (c.dir && open.has(c.path)) walk(c, depth + 1); }
  })(root, 0);
  return out;
}

/** The selection after a click on `path` among the visible file rows `order`: alone, toggled (⌘), or the range from
    the anchor (⇧). Returns { sel: string[], anchor }. */
export function clickSelect(sel, anchor, path, order, { toggle = false, range = false } = {}) {
  if (range && anchor != null && order.includes(anchor)) {
    const a = order.indexOf(anchor), b = order.indexOf(path);
    const [i, j] = a < b ? [a, b] : [b, a];
    return { sel: order.slice(i, j + 1), anchor };
  }
  if (toggle) return { sel: sel.includes(path) ? sel.filter((p) => p !== path) : [...sel, path], anchor: path };
  return { sel: [path], anchor: path };
}
