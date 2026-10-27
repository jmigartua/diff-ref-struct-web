// diff-ref-struct Web · the browser host's pure parts (step 2): the separator's clamp, the tree from a file list, the
// dispatcher's stub answers and its text forms of the path commands, the package's loader, the files' names.
// node --test web/tests/
import test from "node:test";
import assert from "node:assert/strict";
import { readFileSync, readdirSync, existsSync } from "node:fs";
import { clampWidth, dragged, keyed, kept, keep, LIMITS, RAIL } from "../../app/ui/host-web/separator.js";
import { buildTree, groupByKind, filesUnder, visibleRows, clickSelect } from "../../app/ui/host-web/tree.js";
import { createDispatcher, refusal, missing, loadPackage, PATH_COMMANDS, relativeTo, setBase } from "../../app/ui/host-web/dispatch.js";
import { virtualRoot, isHidden, createFiles, fromFileList } from "../../app/ui/host-web/files.js";

const HOST = new URL("../../app/ui/host-web/", import.meta.url);

// ------------------------------------------------------------------ the separator
test("separator: the width clamped to the inspector's limits and the views' least", () => {
  assert.equal(clampWidth(100, 2000), LIMITS.min);
  assert.equal(clampWidth(5000, 2000), LIMITS.max);
  assert.equal(clampWidth(400, 2000), 400);
  assert.equal(clampWidth(400.6, 2000), 401);
  // a narrow window: the views keep their 360 px, the inspector never under its least
  assert.equal(clampWidth(500, 700), 700 - LIMITS.views);
  assert.equal(clampWidth(500, 400), LIMITS.min);
  assert.equal(clampWidth(NaN, 2000), LIMITS.start);
  assert.equal(clampWidth(undefined), LIMITS.start);
  assert.ok(RAIL < LIMITS.min);
});
test("separator: a drag and the arrow keys move the width, clamped", () => {
  assert.equal(dragged(320, 320, 420, 2000), 420);
  assert.equal(dragged(320, 320, 0, 2000), LIMITS.min);
  assert.equal(dragged(320, 320, 1900, 2000), LIMITS.max);
  assert.equal(keyed(320, "ArrowRight", false, 2000), 336);
  assert.equal(keyed(320, "ArrowLeft", true, 2000), 256);
  assert.equal(keyed(320, "Enter", false, 2000), 320);
});
test("separator: the width and the rail kept in the storage, the defaults when it keeps nothing or refuses", () => {
  const m = new Map(); const st = { getItem: (k) => (m.has(k) ? m.get(k) : null), setItem: (k, v) => m.set(k, v) };
  assert.deepEqual(kept(st), { width: LIMITS.start, rail: false });
  keep(st, { width: 412, rail: true });
  assert.deepEqual(kept(st), { width: 412, rail: true });
  m.set("drs-web.inspector", "{not json");
  assert.deepEqual(kept(st), { width: LIMITS.start, rail: false });
  const refusing = { getItem: () => { throw new Error("denied"); }, setItem: () => { throw new Error("denied"); } };
  assert.deepEqual(kept(refusing), { width: LIMITS.start, rail: false });
  assert.doesNotThrow(() => keep(refusing, { width: 300, rail: false }));
});

// ------------------------------------------------------------------ the tree
const LIST = [
  { path: "run/scan_10.xye", size: 10 }, { path: "run/scan_2.xye", size: 20 }, { path: "run/sub/a.dat", size: 5 },
  { path: "notes.txt", size: 1 }, { path: "./Plots/set.drv", size: 3 }, { path: "run//scan_2.xye", size: 21 }, { path: "", size: 9 },
];
test("tree: folders first, names in reading order, counts and sizes, the same path once", () => {
  const t = buildTree(LIST, "beamtime");
  assert.equal(t.name, "beamtime");
  assert.deepEqual(t.children.map((c) => c.name), ["Plots", "run", "notes.txt"]);
  const run = t.children[1];
  assert.deepEqual(run.children.map((c) => c.name), ["sub", "scan_2.xye", "scan_10.xye"]);
  assert.equal(run.files, 3);
  assert.equal(run.children[1].size, 21);   // the last of the same path wins
  assert.equal(t.files, 5);
  assert.equal(t.size, 3 + 10 + 21 + 5 + 1);
  assert.equal(run.children[0].path, "run/sub");
});
test("tree: every file under a folder, the rows shown by the folders open, a click's selection", () => {
  const t = buildTree(LIST);
  const run = t.children.find((c) => c.name === "run");
  assert.deepEqual(filesUnder(run), ["run/sub/a.dat", "run/scan_2.xye", "run/scan_10.xye"]);
  assert.deepEqual(visibleRows(t, new Set()).map((r) => r.node.path), ["Plots", "run", "notes.txt"]);
  const rows = visibleRows(t, new Set(["run"]));
  assert.deepEqual(rows.map((r) => [r.node.path, r.depth]), [["Plots", 0], ["run", 0], ["run/sub", 1], ["run/scan_2.xye", 1], ["run/scan_10.xye", 1], ["notes.txt", 0]]);
  const order = ["a", "b", "c", "d"];
  assert.deepEqual(clickSelect([], null, "b", order), { sel: ["b"], anchor: "b" });
  assert.deepEqual(clickSelect(["b"], "b", "d", order, { range: true }), { sel: ["b", "c", "d"], anchor: "b" });
  assert.deepEqual(clickSelect(["b"], "b", "a", order, { range: true }).sel, ["a", "b"]);
  assert.deepEqual(clickSelect(["b"], "b", "c", order, { toggle: true }), { sel: ["b", "c"], anchor: "c" });
  assert.deepEqual(clickSelect(["b", "c"], "c", "b", order, { toggle: true }).sel, ["c"]);
});

test("tree: the files grouped by the core's kind, data first, what no reader knows last, the pending apart", () => {
  const k = { "a.xye": { kind: "data", id: "xye", name: "xye columns" }, "b.xye": { kind: "data", id: "xye", name: "xye columns" }, "x.buf": { kind: "project", id: "buffer_file", name: "buffer" }, "m.mlf": { kind: "unknown", id: "unknown", name: null } };
  const t = groupByKind([{ path: "m.mlf", size: 1 }, { path: "x.buf", size: 1 }, { path: "b.xye", size: 2 }, { path: "a.xye", size: 3 }, { path: "new.dat", size: 4 }], (p) => k[p] || null, "run");
  assert.deepEqual(t.children.map((g) => [g.name, g.files]), [["xye columns", 2], ["buffer", 1], ["not recognised by a reader", 1], ["being described by the core…", 1]]);
  assert.deepEqual(t.children[0].children.map((f) => f.path), ["a.xye", "b.xye"]);
  assert.equal(t.files, 5); assert.equal(t.size, 11);
});

// ------------------------------------------------------------------ the dispatcher
/** The host's files as the dispatcher asks them, from a map of path → text (all of them data unless `kinds` says). */
const fakeFs = (texts = {}) => {
  const downloads = [], clip = [];
  return {
    downloads, clip,
    paths: () => Object.keys(texts),
    has: (p) => p in texts,
    size: (p) => (texts[p] || "").length,
    stamp: (p) => `${(texts[p] || "").length}:1`,
    text: async (p) => { if (!(p in texts)) throw `${p}: not read`; return texts[p]; },
    download: (name, text) => downloads.push([name, text]),
    clipboardText: async (t) => { clip.push(t); },
    clipboardPng: async () => {},
  };
};
// the commands the Plots module and W1 name (app/src-tauri/src/modules/plots.rs)
const COMMANDS = ["plots_open_files", "plots_open_set", "plots_open_buffer", "plots_new", "plots_add", "plots_map", "plots_map_region", "plots_range_warning", "plots_displace", "plots_project", "plots_windows", "plots_derive", "plots_save", "plots_close", "plots_extent", "describe"];
test("dispatcher, stub: every command refused with one plain sentence naming W1 and the command", async () => {
  const d = createDispatcher({ pkg: null, fs: fakeFs({ "/local/x/a.xy": "1 2" }) });
  assert.equal(d.ready, false);
  for (const c of COMMANDS) {
    await assert.rejects(d.invoke(c, { id: 1, paths: ["/local/x/a.xy"], path: "/local/x/a.xy" }), (e) => {
      assert.equal(typeof e, "string", "rejected with a string, as Tauri rejects");
      assert.equal(e, refusal(c));
      assert.match(e, /W1/); assert.match(e, new RegExp(`"${c}"`));
      return true;
    });
  }
  assert.equal(await d.propose(["/local/x/a.xy"]), null);
});
test("dispatcher, stub: the host's own services answer without the core (clipboard, a text downloaded)", async () => {
  const fs = fakeFs();
  const d = createDispatcher({ pkg: null, fs });
  await d.invoke("clipboard_text", { text: "a\tb" });
  assert.deepEqual(fs.clip, ["a\tb"]);
  assert.equal(await d.invoke("plots_write_text", { path: "/local/x/Exports/set map.html", text: "<html>" }), 6);
  assert.deepEqual(fs.downloads, [["set map.html", "<html>"]]);
});
/** A package as drs-wasm 0.79.0 is (core/wasm/README.md): `dispatch` a promise of the parsed value, rejected with
    the core's words; `commands()`, `version()`. It records each call; `describe` says data unless `kinds` names a file. */
const fakePkg = (answers = {}, { kinds = {}, served = null } = {}) => {
  const calls = [];
  return {
    calls,
    version: () => "0.79.0",
    commands: () => served || ["plots_open_from", "plots_add_from", "plots_open_buffer_from", "plots_open_set_from", "plots_save", "plots_map", "describe"],
    dispatch: async (cmd, json) => {
      const a = JSON.parse(json);
      if (cmd === "describe") return kinds[a.name] || { kind: "data", id: "xye", name: "x, y, σ columns" };
      calls.push([cmd, a]);
      if (answers[cmd] instanceof Error) throw answers[cmd].message;
      return cmd in answers ? answers[cmd] : { cmd };
    },
  };
};
test("dispatcher, package: a set command goes as it is, its JSON both ways, the answer awaited (a value or JSON text)", async () => {
  const pkg = fakePkg({ plots_map: { map: { rows: 2 } } });
  const d = createDispatcher({ pkg, fs: fakeFs() });
  assert.equal(d.ready, true); assert.equal(d.version, "0.79.0");
  assert.deepEqual(await d.invoke("plots_map", { id: 3, scale: "log", reference: null }), { map: { rows: 2 } });
  assert.deepEqual(pkg.calls, [["plots_map", { id: 3, scale: "log", reference: null }]]);
  const text = createDispatcher({ pkg: { dispatch: async () => '{"a":1}' }, fs: fakeFs() });
  assert.deepEqual(await text.invoke("plots_get", { id: 1 }), { a: 1 });
});
test("dispatcher, package: the core's refusal reaches the module as its own message", async () => {
  const d = createDispatcher({ pkg: fakePkg({ plots_map: new Error("no plot set 9") }), fs: fakeFs() });
  await assert.rejects(d.invoke("plots_map", { id: 9 }), (e) => e === "no plot set 9");
});
test("dispatcher, package: a command the package does not serve refused at once in the core's words, without its list", async () => {
  const pkg = fakePkg();
  const d = createDispatcher({ pkg, fs: fakeFs() });
  assert.equal(d.has("plots_extent"), false);
  for (let k = 0; k < 3; k++) await assert.rejects(d.invoke("plots_extent", { id: 1, scale: "log", members: null }), (e) => e === missing("plots_extent", "0.79.0"));
  assert.match(missing("plots_extent", "0.79.0"), /^no command "plots_extent" in the browser's core \(drs-wasm 0\.79\.0\)/);
  assert.deepEqual(pkg.calls, [], "never sent to the package");
  assert.equal(d.calls.plots_extent, 3);
});
test("dispatcher, package: data files by path go as names and text (plots_open_files → plots_open_from, plots_add → plots_add_from)", async () => {
  const pkg = fakePkg();
  const d = createDispatcher({ pkg, fs: fakeFs({ "/local/r/a.xye": "A", "/local/r/b.xye": "B" }) });
  await d.invoke("plots_open_files", { paths: ["/local/r/a.xye", "/local/r/b.xye"], root: "/local/r" });
  await d.invoke("plots_add", { id: 4, paths: ["/local/r/b.xye"] });
  assert.deepEqual(pkg.calls, [
    ["plots_open_from", { files: [{ name: "/local/r/a.xye", text: "A" }, { name: "/local/r/b.xye", text: "B" }] }],
    ["plots_add_from", { id: 4, files: [{ name: "/local/r/b.xye", text: "B" }] }],
  ]);
  assert.deepEqual(Object.keys(PATH_COMMANDS).sort(), ["plots_add", "plots_derive", "plots_open_buffer", "plots_open_files", "plots_open_set", "plots_save", "plots_write_text"]);
});
test("dispatcher, package: the patterns in the bridge's words (member, not index: W11), every number as the core gave it", async () => {
  const pats = [{ index: 0, key: "a.xye#1", x: [1, 2], y: [3, 4], sigma: [null, 2] }, { index: 1, key: "b.xye#1", x: [1], y: [5], sigma: [1] }];
  const d = createDispatcher({ pkg: fakePkg({ plots_patterns: pats }, { served: ["plots_patterns"] }), fs: fakeFs() });
  const got = await d.invoke("plots_patterns", { id: 1, members: null });
  assert.deepEqual(got.map((p) => p.member), [0, 1]);
  assert.deepEqual(got.map((p) => [p.x, p.y, p.sigma]), pats.map((p) => [p.x, p.y, p.sigma]));
  // a core that already says member is passed as it is
  const as = [{ member: 3, x: [], y: [], sigma: [] }];
  assert.deepEqual(await createDispatcher({ pkg: fakePkg({ plots_patterns: as }, { served: ["plots_patterns"] }), fs: fakeFs() }).invoke("plots_patterns", { id: 1 }), as);
});
test("dispatcher, package: a derived file (Sum, Average, Difference) downloaded, the module told where it went", async () => {
  const fs = fakeFs();
  const pkg = fakePkg({ plots_derive: { path: "Plots/s sum.xys", relative: "Plots/s sum.xys", points: 301, dropped: 0, derivation: "sum", text: "! sum\n16 1 1\n" } }, { served: ["plots_derive"] });
  const r = await createDispatcher({ pkg, fs }).invoke("plots_derive", { id: 1, members: [0, 1], kind: "sum", root: "/local/Tb", name: "s" });
  assert.deepEqual(fs.downloads, [["s sum.xys", "! sum\n16 1 1\n"]]);
  assert.match(r.relative, /^s sum\.xys, downloaded to your browser's Downloads/);
  assert.equal(r.points, 301);
});
const FOLDER = { "/local/Tb/Tb_vs_T.buf": "BUF", "/local/Tb/tb_t01.dat": "1", "/local/Tb/sub/tb_t02.dat": "2", "/local/Tb/winplotr.mlf": "MLF", "/local/Tb/Plots/s.drv": "DRV" };
const KINDS = { "Tb_vs_T.buf": { kind: "project", id: "buffer_file", name: "list of file names (.buf)" }, "winplotr.mlf": { kind: "unknown", id: "unknown", name: null }, "s.drv": { kind: "project", id: "plot_set", name: "plot set" } };
test("dispatcher, package 0.79.0: a buffer gets the folder's data files (the core's describe) named from its folder; the core takes its own", async () => {
  const pkg = fakePkg({}, { kinds: KINDS });
  const d = createDispatcher({ pkg, fs: fakeFs(FOLDER) });
  await d.invoke("plots_open_buffer", { path: "/local/Tb/Tb_vs_T.buf", root: "/local/Tb" });
  assert.deepEqual(pkg.calls, [["plots_open_buffer_from", { path: "/local/Tb/Tb_vs_T.buf", text: "BUF", files: [{ name: "tb_t01.dat", text: "1" }, { name: "sub/tb_t02.dat", text: "2" }], root: "/local/Tb" }]]);
  assert.equal((await d.kindOf("/local/Tb/winplotr.mlf")).kind, "unknown");
});
test("dispatcher, package 0.79.0: a set file in Plots/ gets the data files named from the project's folder; Save downloads", async () => {
  const fs = fakeFs(FOLDER);
  const pkg = fakePkg({ plots_open_set_from: [{ id: 1 }, { order: "variable" }], plots_save: { relative: "Plots/my set.drv", text: "format = 1\n" } }, { kinds: KINDS });
  const d = createDispatcher({ pkg, fs });
  assert.deepEqual(await d.invoke("plots_open_set", { path: "/local/Tb/Plots/s.drv", root: "/local/Tb" }), [{ id: 1 }, { order: "variable" }]);
  assert.deepEqual(pkg.calls[0], ["plots_open_set_from", { path: "/local/Tb/Plots/s.drv", text: "DRV", files: [{ name: "tb_t01.dat", text: "1" }, { name: "sub/tb_t02.dat", text: "2" }], root: "/local/Tb" }]);
  const saved = await d.invoke("plots_save", { id: 2, root: "/local/Tb", name: "my set", settings: { order: "variable" } });
  assert.equal(saved.path, "Plots/my set.drv"); assert.match(saved.note, /^Downloaded as my set\.drv/);
  assert.deepEqual(pkg.calls[1], ["plots_save", { id: 2, root: "/local/Tb", name: "my set", settings: { order: "variable" } }]);
  assert.deepEqual(fs.downloads, [["my set.drv", "format = 1\n"]]);
});
test("dispatcher, package with plots_names_in (W5): only the files the core says the buffer names are read", async () => {
  const pkg = fakePkg({ plots_names_in: { files: ["/local/Tb/tb_t01.dat", "/local/Tb/gone.dat"] } }, { kinds: KINDS, served: ["plots_names_in", "plots_open_buffer_from", "describe"] });
  const d = createDispatcher({ pkg, fs: fakeFs(FOLDER) });
  await d.invoke("plots_open_buffer", { path: "/local/Tb/Tb_vs_T.buf", root: "/local/Tb" });
  assert.deepEqual(pkg.calls, [["plots_names_in", { path: "/local/Tb/Tb_vs_T.buf", text: "BUF", root: "/local/Tb" }], ["plots_open_buffer_from", { path: "/local/Tb/Tb_vs_T.buf", text: "BUF", files: [{ name: "/local/Tb/tb_t01.dat", text: "1" }], root: "/local/Tb" }]]);
});
test("dispatcher: the proposer asked with the data files only when the package has one (W2)", async () => {
  assert.equal(await createDispatcher({ pkg: fakePkg({}, { kinds: KINDS }), fs: fakeFs(FOLDER) }).propose(Object.keys(FOLDER)), null);
  const pkg = fakePkg({ plots_propose: { proposals: [{ title: "2 scans" }] } }, { kinds: KINDS, served: ["plots_propose", "describe"] });
  const r = await createDispatcher({ pkg, fs: fakeFs(FOLDER) }).propose(Object.keys(FOLDER));
  assert.deepEqual(r, { proposals: [{ title: "2 scans" }] });
  assert.deepEqual(pkg.calls, [["plots_propose", { files: [{ name: "/local/Tb/tb_t01.dat", text: "1" }, { name: "/local/Tb/sub/tb_t02.dat", text: "2" }] }]]);
});
test("paths: relative to a folder, and a set file's base as the core's set_base", () => {
  assert.equal(relativeTo("/local/Tb/sub/a.dat", "/local/Tb"), "sub/a.dat");
  assert.equal(relativeTo("/local/Tb/a.dat", "/local/Tb/sub"), "../a.dat");
  assert.equal(relativeTo("/local/Tb/a.dat", "/local/Tb"), "a.dat");
  assert.equal(setBase("/local/Tb/Plots/s.drv"), "/local/Tb");
  assert.equal(setBase("/local/Tb/elsewhere/s.drv"), "/local/Tb/elsewhere");
});

// ------------------------------------------------------------------ the package's loader
const MANIFEST = JSON.parse(readFileSync(new URL("wasm/manifest.json", HOST), "utf8"));
test("loader: the stub's manifest asks for no package (no missing file fetched)", async () => {
  const loads = [];
  const stub = { package: null, why: "the core's WebAssembly build is request W1" };
  const got = await loadPackage({ base: "http://h/host-web/", getJson: async (u) => { assert.equal(u, "http://h/host-web/wasm/manifest.json"); return stub; }, load: async (u) => { loads.push(u); } });
  assert.equal(got.pkg, null); assert.match(got.why, /W1/); assert.deepEqual(loads, []);
  const none = await loadPackage({ base: "http://h/host-web/", getJson: async () => { throw "HTTP 404"; }, load: async () => { throw new Error("not asked"); } });
  assert.equal(none.pkg, null); assert.match(none.why, /manifest/);
});
test("loader: the manifest in wasm/ as it is here: the stub's, or a built package whose files are there", () => {
  if (!MANIFEST.package) { assert.match(MANIFEST.why, /W1/); return; }
  assert.ok(existsSync(new URL("wasm/" + MANIFEST.package, HOST)), `wasm/${MANIFEST.package}`);
  assert.ok(existsSync(new URL("wasm/" + (MANIFEST.wasm || "drs_wasm_bg.wasm"), HOST)), "the .wasm beside it");
  assert.ok(MANIFEST.bytes > 0);
});
test("loader: a package named by the manifest is loaded, initialised and checked for dispatch", async () => {
  let inited = 0;
  const mod = { default: async () => { inited++; }, dispatch: async () => null, version: () => "0.79.0" };
  const got = await loadPackage({ base: "http://h/host-web/", getJson: async () => ({ package: "drs_wasm.js" }), load: async (u) => { assert.equal(u, "http://h/host-web/wasm/drs_wasm.js"); return mod; } });
  assert.equal(got.pkg, mod); assert.equal(inited, 1); assert.equal(got.version, "0.79.0");
  const bad = await loadPackage({ base: "http://h/host-web/", getJson: async () => ({ package: "drs_wasm.js" }), load: async () => ({ default: async () => {} }) });
  assert.equal(bad.pkg, null); assert.match(bad.why, /dispatch/);
});

// ------------------------------------------------------------------ the files
test("files: virtual paths under a root by the folder's name, hidden ones left out, a name taken twice", async () => {
  assert.equal(virtualRoot("beam time"), "/local/beam time");
  assert.equal(virtualRoot("a/b"), "/local/a-b");
  assert.equal(virtualRoot(""), "/local/files");
  assert.ok(isHidden(".git/config")); assert.ok(isHidden("a/.DS_Store")); assert.ok(!isHidden("a/b.xy"));
  const f = createFiles({}, {});
  const file = (text, lastModified = 1) => ({ size: text.length, lastModified, getFile: async () => ({ text: async () => text }) });
  f.replace("run", [{ rel: "a.xye", ...file("A") }, { rel: ".hidden", ...file("H") }, { rel: "sub/b.xye", ...file("B") }]);
  assert.equal(f.root, "/local/run");
  assert.deepEqual(f.list().map((x) => x.path), ["a.xye", "sub/b.xye"]);
  assert.equal(await f.text("/local/run/sub/b.xye"), "B");
  await assert.rejects(f.text("/local/run/none.xye"), /not among the files read/);
  assert.deepEqual(f.add([{ rel: "a.xye", ...file("A") }]), ["a.xye"]);                   // the same file again: its place
  assert.deepEqual(f.add([{ rel: "a.xye", ...file("other", 2) }, { rel: ".x", ...file("H") }]), ["a (2).xye"]);   // another of that name: a (2)
  assert.deepEqual(f.list().map((x) => x.path), ["a.xye", "sub/b.xye", "a (2).xye"]);
  const fresh = createFiles({}, {}); fresh.add([{ rel: "x.drv", ...file("D") }]);
  assert.equal(fresh.root, "/local/dropped files");
});
test("files: <input webkitdirectory>'s list, the folder's name taken off each path", () => {
  const { folder, entries } = fromFileList([{ name: "a.xye", webkitRelativePath: "run/a.xye", size: 3, lastModified: 1 }, { name: "b.xye", webkitRelativePath: "run/s/b.xye", size: 4, lastModified: 1 }]);
  assert.equal(folder, "run"); assert.deepEqual(entries.map((e) => e.rel), ["a.xye", "s/b.xye"]);
  const loose = fromFileList([{ name: "c.drv", webkitRelativePath: "", size: 1, lastModified: 1 }]);
  assert.equal(loose.folder, null); assert.deepEqual(loose.entries.map((e) => e.rel), ["c.drv"]);
});

// ------------------------------------------------------------------ the host loads the shell and the module as they are
test("host: the shell's files and the Plots module imported from app/ui/ (no copy); the page the modules' host (W6)", () => {
  const src = readFileSync(new URL("host.js", HOST), "utf8");
  const imports = [...src.matchAll(/^import\s+(?:[^"]*?from\s+)?"([^"]+)"/gm)].map((m) => m[1]);
  for (const p of ["../core.js", "../registry.js", "../layout.js", "../modules/plots/index.js"]) assert.ok(imports.includes(p), p);
  assert.match(src, /window\.__DRS_HOST__ = \{/); assert.doesNotMatch(src, /__TAURI__/);
  const own = readdirSync(HOST).filter((n) => n.endsWith(".js"));
  assert.deepEqual(own.sort(), ["dispatch.js", "files.js", "hash.js", "host.js", "live.js", "separator.js", "sw.js", "tree.js"]);
  const html = readFileSync(new URL("index.html", HOST), "utf8");
  assert.match(html, /never leave this machine/);
  assert.match(html, /<script type="module" src="host.js">/);
  assert.match(html, /id="proposals"/); assert.match(html, /data-by="kinds"/);
  assert.match(html, /rel="manifest" href="manifest.webmanifest"/); assert.match(html, /id="live"/); assert.match(html, /id="offline-line"/);
  const man = JSON.parse(readFileSync(new URL("manifest.webmanifest", HOST), "utf8"));
  assert.equal(man.start_url, "./"); assert.equal(man.display, "standalone"); assert.ok(man.icons.length);
  const sw = readFileSync(new URL("sw.js", HOST), "utf8");
  assert.match(sw, /^const VERSION = "dev";$/m, "stage.mjs stamps this line");
});
