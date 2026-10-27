/* diff-ref-struct Web · the user's files (web/DESIGN.md §3): a folder picked (`showDirectoryPicker`, Chromium;
   `<input webkitdirectory>` elsewhere) or files and folders dropped, read as text in the browser when the core needs
   them, never uploaded. Each file gets a path under a virtual root (`/local/<folder>/<path in it>`), which is what the
   module and the core see where the app gives real paths; the handles are kept for live mode (step 4). */

/** The virtual root of a folder by its name: the module's `state.project.root`. */
export const virtualRoot = (name) => "/local/" + (String(name || "files").replace(/[\/\\]/g, "-") || "files");
/** A hidden file or folder (".git", ".DS_Store", ".drs") on the way: left out of the tree. */
export const isHidden = (rel) => String(rel).split("/").some((p) => p.startsWith("."));

/** The host's files: entries by their virtual path, and the services the dispatcher asks of them. */
export function createFiles(doc = globalThis.document, nav = globalThis.navigator) {
  const entries = new Map();     // virtual path → { rel, name, size, lastModified, getFile, handle }
  let root = null, name = null;
  const files = {
    get root() { return root; },
    get name() { return name; },
    /** A new folder in place of the old one (`entries` [{ rel, size, lastModified, getFile, handle }]). */
    replace(folderName, list) { entries.clear(); name = folderName; root = virtualRoot(folderName); files.add(list); },
    /** Files joining the folder shown (dropped loose). The same file again (its name, size and date) takes its place;
        another file of a name taken gets " (2)", " (3)" … before its extension. Returns the paths in the folder each
        one landed at, in the order given (the hidden ones left out). */
    add(list) {
      const landed = [];
      if (!root) { name = "dropped files"; root = virtualRoot(name); }
      const same = (a, b) => a.size === b.size && a.lastModified === b.lastModified;
      for (const e of list) {
        if (isHidden(e.rel)) continue;
        let rel = e.rel, k = 1;
        while (entries.has(root + "/" + rel) && !same(entries.get(root + "/" + rel), e)) { k += 1; rel = e.rel.replace(/(\.[^./]*)?$/, (x) => ` (${k})${x}`); }
        entries.set(root + "/" + rel, { ...e, rel, name: rel.split("/").pop() });
        landed.push(rel);
      }
      return landed;
    },
    /** The files shown (a proposal's set file, held in the page, is not one). */
    list: () => [...entries.values()].filter((e) => !e.held).map((e) => ({ path: e.rel, size: e.size })),
    /** Every file's virtual path, in the order read (the dispatcher's: the files a buffer or a set file may name). */
    paths: () => [...entries.keys()].filter((p) => !entries.get(p).held),
    absOf: (rel) => root + "/" + rel,
    relOf: (p) => (root && p.startsWith(root + "/") ? p.slice(root.length + 1) : p),
    has: (p) => entries.has(p),
    entry: (p) => entries.get(p) || null,
    size: (p) => (entries.get(p) || {}).size || 0,
    /** A file's modification time in ms (the browser's File.lastModified), for the core's propose (W2). */
    modified: (p) => { const e = entries.get(p); return e && e.lastModified ? e.lastModified : null; },
    /** The files shown by their path in the folder and stamp, for live mode's diff (live.js). */
    stamps: () => new Map([...entries.values()].filter((e) => !e.held).map((e) => [e.rel, `${e.size}:${e.lastModified}`])),
    /** The folder read again (live mode): its files replaced by `list`, the held ones kept. */
    sync(list) {
      for (const [p, e] of [...entries]) if (!e.held) entries.delete(p);
      for (const e of list) if (!isHidden(e.rel)) entries.set(root + "/" + e.rel, { ...e, name: e.rel.split("/").pop() });
    },
    /** What tells a file read again from the one read before: its size and date. */
    stamp: (p) => { const e = entries.get(p); return e ? `${e.size}:${e.lastModified}` : ""; },
    /** A text held in the page under a path of the folder (a proposal's set file), not listed in the tree. */
    put(rel, text) {
      if (!root) { name = "dropped files"; root = virtualRoot(name); }
      const f = { size: text.length, lastModified: Date.now(), getFile: async () => ({ text: async () => text }), handle: null, held: true };
      entries.set(root + "/" + rel, { ...f, rel, name: rel.split("/").pop() });
      return root + "/" + rel;
    },
    async text(p) {
      const e = entries.get(p);
      if (!e) throw `${p}: not among the files read in this page (pick its folder or drop it first).`;
      return (await e.getFile()).text();
    },
    /** A text handed to the user as a download (the browser's Downloads): Save Set's .drv, the viewer page. */
    download(fileName, text, type = "text/plain") {
      const url = URL.createObjectURL(new Blob([text], { type }));
      const a = doc.createElement("a"); a.href = url; a.download = fileName; a.hidden = true;
      doc.body.appendChild(a); a.click(); a.remove();
      setTimeout(() => URL.revokeObjectURL(url), 10000);
    },
    clipboardText: (text) => nav.clipboard.writeText(text),
    clipboardPng: (png) => nav.clipboard.write([new ClipboardItem({ "image/png": new Blob([new Uint8Array(png)], { type: "image/png" }) })]),
  };
  return files;
}

/** Every file under a directory handle (the File System Access API), paths relative to it; hidden ones skipped. */
export async function fromDirectoryHandle(dir, prefix = "") {
  const out = [];
  for await (const h of dir.values()) {
    if (h.name.startsWith(".")) continue;
    const rel = prefix + h.name;
    if (h.kind === "directory") out.push(...(await fromDirectoryHandle(h, rel + "/")));
    else { const f = await h.getFile(); out.push({ rel, size: f.size, lastModified: f.lastModified, getFile: () => h.getFile(), handle: h }); }
  }
  return out;
}

/** The files of `<input webkitdirectory>` (or plain `multiple`): the folder's name and the files under it. */
export function fromFileList(list) {
  const fs = [...list];
  const first = fs.find((f) => f.webkitRelativePath);
  const folder = first ? first.webkitRelativePath.split("/")[0] : null;
  const entries = fs.map((f) => {
    const rel = f.webkitRelativePath ? f.webkitRelativePath.split("/").slice(1).join("/") : f.name;
    return { rel, size: f.size, lastModified: f.lastModified, getFile: async () => f, handle: null };
  });
  return { folder, entries };
}

/** What was dropped: one folder (it becomes the folder shown, its handle kept where the browser gives one), or loose
    files and folders (they join the folder shown). The items are taken in the drop's own task, as the browser
    requires; they are read after. */
export async function fromDrop(dt) {
  const items = [...(dt.items || [])].filter((i) => i.kind === "file");
  const got = items.map((i) => (i.getAsFileSystemHandle ? i.getAsFileSystemHandle() : Promise.resolve(entryOf(i))));
  const asFiles = items.map((i) => (i.getAsFile ? i.getAsFile() : null));   // the fallback, taken now as well
  const resolved = await Promise.all(got.map((p) => p.catch(() => null)));
  // an item the browser gives no handle or entry for (a drop from another page, a synthetic one): its file as it is
  const plain = items.length ? asFiles.filter((f, k) => f && !resolved[k]) : [...(dt.files || [])];
  const handles = resolved.filter(Boolean);
  if (plain.length && !handles.length) return { folder: null, entries: fromFileList(plain).entries };
  if (handles.length === 1 && isDir(handles[0])) return { folder: handles[0].name, entries: await walk(handles[0], ""), handle: handles[0].kind ? handles[0] : null };
  const entries = plain.length ? fromFileList(plain).entries : [];
  for (const h of handles) entries.push(...(isDir(h) ? await walk(h, h.name + "/") : await walk(h, "")));
  return { folder: null, entries };
}
function entryOf(item) { return item.webkitGetAsEntry ? item.webkitGetAsEntry() : null; }
const isDir = (h) => h.kind === "directory" || h.isDirectory === true;
/** A handle (File System Access) or an entry (the older FileSystemEntry API), walked alike: a folder's files with
    `prefix` before their paths in it, a file as `prefix` + its name. */
async function walk(h, prefix) {
  if (h.kind === "directory") return fromDirectoryHandle(h, prefix);
  if (h.kind === "file") return [await fileEntry(prefix + h.name, () => h.getFile(), h)];
  if (h.isFile) return [await fileEntry(prefix + h.name, () => new Promise((ok, bad) => h.file(ok, bad)), null)];
  const out = [], reader = h.createReader();
  for (;;) {
    const batch = await new Promise((ok, bad) => reader.readEntries(ok, bad));   // a folder's entries come in batches
    if (!batch.length) break;
    for (const c of batch) if (!c.name.startsWith(".")) out.push(...(await walk(c, c.isDirectory ? prefix + c.name + "/" : prefix)));
  }
  return out;
}
async function fileEntry(rel, getFile, handle) { const f = await getFile(); return { rel, size: f.size, lastModified: f.lastModified, getFile, handle }; }
