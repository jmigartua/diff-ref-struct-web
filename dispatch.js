/* diff-ref-struct Web · the dispatcher (I10, R34; web/DESIGN.md §2): the browser host's one `invoke(cmd, args)`, the
   same names and JSON as the app's Tauri commands (app/src-tauri/src/modules/plots.rs), served by the core compiled to
   WebAssembly (drs-wasm, core/wasm/README.md). The browser computes nothing (L139): this file passes the module's
   questions to the core and its answers back, and does the one thing the Tauri side does besides, reading and writing
   files, here as names and text read in the browser (nothing uploaded) and downloads.

   - With the package (app/ui/host-web/wasm/, announced by wasm/manifest.json): every command goes to its
     `dispatch(cmd, argsJson)`, a promise of the parsed value rejected with the core's words (core 0.79.0, L424, L425);
     the commands that name paths go as the package's text forms (PATH_COMMANDS): the files' names and contents in
     place of the paths, named as the core resolves them; what a file is, the core says (`describe`).
   - Without it (the stub): every command is refused with one plain sentence naming W1.
   Pure: no DOM, no fetch of its own (the loader and the files are given), so that node tests it. */

/** The refusal the stub gives every command: one sentence, the command named, the request that brings the core. */
export function refusal(cmd) {
  return `The core is not in this page yet: "${cmd}" waits for its WebAssembly build (request W1, web/REQUESTS.md). Nothing is computed in the browser, so nothing is answered without it.`;
}
/** A command the package does not serve, in the core's own words (drs-wasm's "no command … in the browser's core
    (drs-wasm X)") without the list of every command it serves, which the status line cannot hold. */
export function missing(cmd, version) {
  return `no command "${cmd}" in the browser's core (drs-wasm ${version || "?"}): the app's core has it, this build not yet`;
}

// ------------------------------------------------------------------ paths (strings: the host's virtual paths)
const dirOf = (p) => String(p).replace(/\/[^/]*$/, "");
const baseName = (p) => String(p || "").split("/").pop();
/** `path` relative to the folder `dir` ("a/b.xy", "../c.xy"), both absolute virtual paths. */
export function relativeTo(path, dir) {
  const a = String(path).split("/").filter(Boolean), b = String(dir).split("/").filter(Boolean);
  let k = 0; while (k < a.length - 1 && k < b.length && a[k] === b[k]) k++;
  return [...Array(b.length - k).fill(".."), ...a.slice(k)].join("/");
}
/** The folder a set file's members are named from: the core's `plots::set_base` (its folder, or the project's when
    the file is in Plots/), so that the members are handed over under the names the .drv keeps. */
export function setBase(path) {
  const dir = dirOf(path);
  return baseName(dir) === "Plots" ? dirOf(dir) : dir;
}

/** The commands that name files on disk in the app, and the package's text forms the host sends instead
    (core/wasm/README.md's table). Each: the command the core is asked, and its arguments from the module's with the
    host's files (`h`: the dispatcher's helpers). */
export const PATH_COMMANDS = {
  // data files as one set: W1's open_from, the files by their virtual paths (the core finds a variable in a name)
  plots_open_files: { to: "plots_open_from", args: async (a, h) => ({ files: await h.texts(a.paths), ...h.takeOpening() }) },
  plots_add: { to: "plots_add_from", args: async (a, h) => ({ id: a.id, files: await h.texts(a.paths) }) },
  // a WinPLOTR buffer: its text, and the files it may name under the names it writes them by (relative to its folder)
  plots_open_buffer: { to: "plots_open_buffer_from", args: async (a, h) => { const text = await h.fs.text(a.path); return { path: a.path, text, files: await h.named(a.path, text, dirOf(a.path), a.root), root: a.root ?? null }; } },
  // a plot set file: its text and its members' files under the names it keeps (relative to its base)
  plots_open_set: { to: "plots_open_set_from", args: async (a, h) => { const text = await h.fs.text(a.path); return { path: a.path, text, files: await h.named(a.path, text, setBase(a.path), a.root), root: a.root ?? null }; } },
  // Save: the core gives the .drv as `{ relative, path, text }` (written nowhere); the host downloads it and says so
  // (W7: the module shows the host's note beside the path)
  plots_save: { to: "plots_save", args: (a) => a, after: (r, h) => { const rel = r.path || r.relative, name = baseName(rel); h.fs.download(name, r.text); return { path: rel, note: `Downloaded as ${name} (your browser's Downloads): nothing is written into the folder; drop it on the page with the folder open to open it again` }; } },
  // the Map viewer's page (Viewer Page…, Dev): the module's text, downloaded
  plots_write_text: { host: (a, h) => { h.fs.download(baseName(a.path), a.text); return a.text.length; } },
};

/** What the host answers itself, with or without the core: the clipboard (plot.js, table.js ask the app's). */
export const HOST_COMMANDS = {
  clipboard_text: (a, h) => h.fs.clipboardText(a.text),
  clipboard_png: (a, h) => h.fs.clipboardPng(a.png),
};

/** The largest file the host hands the core as text (a folder's images or archives are not read). */
export const TEXT_LIMIT = 64 * 1024 * 1024;

/** The dispatcher. `pkg`: the wasm package's module (`dispatch`, `commands`, `version`), or null for the stub. `fs`:
    the host's files: { paths(), has(path), size(path), stamp(path), text(path), download(name, text),
    clipboardText(text), clipboardPng(bytes) }. */
export function createDispatcher({ pkg = null, fs = null } = {}) {
  const served = pkg && typeof pkg.commands === "function" ? new Set(pkg.commands()) : null;
  const version = pkg && typeof pkg.version === "function" ? pkg.version() : null;
  const calls = {};
  const has = (cmd) => !!pkg && (!served || served.has(cmd));
  const core = async (cmd, args) => {
    calls[cmd] = (calls[cmd] || 0) + 1;
    if (!pkg) throw refusal(cmd);
    if (!has(cmd)) throw missing(cmd, version);
    let out;
    // drs-wasm 0.79.0 (L424): `dispatch` returns a promise of the parsed value, rejected with the core's words; a
    // string is parsed too
    try { out = await pkg.dispatch(cmd, JSON.stringify(args ?? {})); }
    catch (e) { throw typeof e === "string" ? e : e && e.message ? e.message : String(e); }   // the core's own message, as Tauri gives it
    if (out === undefined || out === null || out === "") return null;
    return typeof out === "string" ? JSON.parse(out) : out;
  };

  // what each file is, the core's `describe`, asked once per file as read (its size and date)
  const kinds = new Map();
  let opening = null;
  // a file's modification time in s, as the core's propose reads it (OpenedFile::modified), where the host knows it
  const modified = (p) => { const t = fs.modified ? fs.modified(p) : null; return t ? { modified: t / 1000 } : {}; };
  async function kindOf(path) {
    const stamp = fs.stamp(path), got = kinds.get(path);
    if (got && got.stamp === stamp) return got.kind;
    let kind;
    if (fs.size(path) > TEXT_LIMIT) kind = { kind: "unknown", id: "large", name: `not read: larger than ${TEXT_LIMIT / 1048576} MB` };
    else kind = await core("describe", { name: baseName(path), text: await fs.text(path) });
    kinds.set(path, { stamp, kind });
    return kind;
  }
  const isData = async (p) => (await kindOf(p)).kind === "data";
  const h = {
    fs,
    async texts(paths) {
      const out = [];
      for (const p of [].concat(paths || [])) out.push({ name: p, text: await fs.text(p) });
      return out;
    },
    /** The files a buffer or a set file names: those the core says (`plots_names_in`, W5, core 0.80.0), by their
        paths, the ones the folder lacks left for the core to name; or, from a package without it (0.79.0), every
        data file of the folder (the core's `describe`) named relative to `dir`, among which the core takes its own. */
    async named(path, text, dir, root) {
      const out = [];
      if (has("plots_names_in")) {
        const { files } = await core("plots_names_in", { path, text, root: root ?? null });
        for (const p of files) if (fs.has(p)) out.push({ name: p, text: await fs.text(p), ...modified(p) });
        return out;
      }
      for (const p of fs.paths()) if (p !== path && (await isData(p))) out.push({ name: relativeTo(p, dir), text: await fs.text(p) });
      return out;
    },
    /** What the next opening of files carries beside them (a proposal's `rebin`, core 0.81.0), taken once. */
    takeOpening() { const o = opening; opening = null; return o || {}; },
  };

  async function invoke(cmd, args = {}) {
    if (HOST_COMMANDS[cmd]) return HOST_COMMANDS[cmd](args, h);
    const p = PATH_COMMANDS[cmd];
    if (p && p.host) return p.host(args, h);
    if (!p) return core(cmd, args);
    if (!pkg) { calls[cmd] = (calls[cmd] || 0) + 1; throw refusal(cmd); }
    const r = await core(p.to, await p.args(args, h));
    return p.after ? p.after(r, h) : r;
  }
  /** The sets the core finds in the files (W2's `plots_propose`), or null when the package has no proposer. */
  async function propose(paths) {
    if (!has("plots_propose")) return null;
    const files = []; for (const p of paths) if (await isData(p)) files.push({ name: p, text: await fs.text(p), ...modified(p) });
    return core("plots_propose", { files });
  }
  /** The arguments the next `plots_open_files` carries beside its files ({ rebin } from a proposal). */
  const nextOpening = (o) => { opening = o && Object.keys(o).length ? o : null; };
  return { invoke, ready: !!pkg, core, has, kindOf, propose, nextOpening, version, calls };
}

/** The package, when the manifest announces one: `{ package: "drs_wasm.js" }` in wasm/manifest.json, written by
    web/tools/build-wasm.sh beside wasm-pack's output (`--target web`: a default export that loads the .wasm, and
    `dispatch`). `load(url)` imports a module (the page's `import()`), `getJson(url)` reads the manifest. Null on the
    stub (the committed manifest says `"package": null`), so that no missing file is ever asked for. */
export async function loadPackage({ base, getJson, load }) {
  let m;
  try { m = await getJson(new URL("wasm/manifest.json", base).href); } catch (e) { return { pkg: null, why: `no wasm/manifest.json (${e})` }; }
  if (!m || !m.package) return { pkg: null, why: m && m.why ? m.why : "the manifest names no package" };
  const mod = await load(new URL("wasm/" + m.package, base).href);
  if (typeof mod.default === "function") await mod.default();
  if (typeof mod.dispatch !== "function") return { pkg: null, why: `wasm/${m.package} exports no dispatch(cmd, args)` };
  return { pkg: mod, why: null, manifest: m, version: typeof mod.version === "function" ? mod.version() : m.core || null };
}
