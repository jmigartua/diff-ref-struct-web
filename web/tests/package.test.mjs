// diff-ref-struct Web · the host's dispatcher on the real package (step 3): drs-wasm as web/tools/build-wasm.sh puts it
// in app/ui/host-web/wasm/, loaded in node, served the files of a folder as the page serves them (virtual paths, text
// read when the core needs it). Skipped, saying why, when no package is built (the committed manifest is the stub's).
// node --test web/tests/
import test from "node:test";
import assert from "node:assert/strict";
import { readFileSync, readdirSync, statSync, existsSync } from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { createDispatcher, missing } from "../../app/ui/host-web/dispatch.js";

const here = path.dirname(fileURLToPath(import.meta.url));
const top = path.resolve(here, "../..");
const wasmDir = path.join(top, "app/ui/host-web/wasm");
const manifest = JSON.parse(readFileSync(path.join(wasmDir, "manifest.json"), "utf8"));
const built = !!manifest.package && existsSync(path.join(wasmDir, manifest.package));
const skip = built ? false : "no package in app/ui/host-web/wasm/ (web/tools/build-wasm.sh builds it)";
const TB = path.join(top, "manuals/fullprof/1-WinPLOTR/1_Tb2pd2sn");
// a team's ESRF run (101 .xye, read only, unpublished): only where DRS_ESRF_RENAMED names its folder
const ESRF = process.env.DRS_ESRF_RENAMED || "/nonexistent/private-esrf-run";

let pkg = null;
async function load() {
  if (pkg) return pkg;
  pkg = await import(path.join(wasmDir, manifest.package));
  pkg.initSync({ module: readFileSync(path.join(wasmDir, manifest.wasm || "drs_wasm_bg.wasm")) });
  return pkg;
}
/** A folder on disk as the page holds one: each file under /local/<name>/, its text read when asked; extra texts held. */
function folderFs(dir, name, held = {}) {
  const root = "/local/" + name, files = new Map();
  for (const f of readdirSync(dir)) if (!f.startsWith(".")) { const st = statSync(path.join(dir, f)); if (st.isFile()) files.set(`${root}/${f}`, { disk: path.join(dir, f), size: st.size, mtime: st.mtimeMs }); }
  const downloads = [];
  return {
    root, downloads, held,
    paths: () => [...files.keys()],
    has: (p) => files.has(p) || p in held,
    size: (p) => (files.get(p) || { size: (held[p] || "").length }).size,
    stamp: (p) => { const f = files.get(p); return f ? `${f.size}:${f.mtime}` : "held"; },
    text: async (p) => (p in held ? held[p] : readFileSync(files.get(p).disk, "utf8")),
    download: (n, t) => downloads.push([n, t]),
  };
}
const count = (xs) => xs.reduce((o, k) => ((o[k] = (o[k] || 0) + 1), o), {});

test("package: the Tb₂Pd₂Sn folder by the core's kinds", { skip }, async () => {
  const fs = folderFs(TB, "Tb2Pd2Sn"), d = createDispatcher({ pkg: await load(), fs });
  const ks = []; for (const p of fs.paths()) ks.push(await d.kindOf(p));
  assert.deepEqual(count(ks.map((k) => (k.kind === "unknown" ? "unknown" : k.id))), { buffer_file: 2, free_steps: 20, unknown: 1 });
});
test("package: each buffer opens with the twenty files it names, taken by the core from the folder's data files", { skip }, async () => {
  const fs = folderFs(TB, "Tb2Pd2Sn"), d = createDispatcher({ pkg: await load(), fs });
  for (const b of ["Tb_vs_T.buf", "data_buffer_file.buf"]) {
    const set = await d.invoke("plots_open_buffer", { path: `${fs.root}/${b}`, root: fs.root });
    assert.equal(set.members.length, 20, b);
    assert.equal(set.apart.length, 1, `${b}: one member without a T (tb_t01)`);
  }
});
test("package: Save downloads the core's .drv; held under Plots/ it opens again with its members and not the folder's others", { skip }, async () => {
  const fs = folderFs(TB, "Tb2Pd2Sn"), d = createDispatcher({ pkg: await load(), fs });
  const five = fs.paths().filter((p) => /tb_t0[1-5]\.dat$/.test(p));
  const set = await d.invoke("plots_open_files", { paths: five, root: fs.root });
  assert.equal(set.members.length, 5);
  const settings = { order: "variable", rule: "none", dx: 0, dy: 0, scale: "linear", prule: "none", pdx: 0, pdy: 0, map_scale: null, reference: null, windows: [], opacity: [], map_range: null, difference_range: null, regions: [] };
  const saved = await d.invoke("plots_save", { id: set.id, root: fs.root, name: "five", settings });
  assert.equal(saved.path, "Plots/five.drv"); assert.match(saved.note, /^Downloaded as five\.drv/);
  const [[name, text]] = fs.downloads;
  assert.equal(name, "five.drv"); assert.match(text, /file = "tb_t02\.dat"/);
  fs.held[`${fs.root}/Plots/five.drv`] = text;
  const [again] = await d.invoke("plots_open_set", { path: `${fs.root}/Plots/five.drv`, root: fs.root });
  assert.deepEqual(again.members.map((m) => m.label).sort(), set.members.map((m) => m.label).sort());
});
test("package: a command it does not serve (plots_extent, 0.79.0) refused at once in the core's words", { skip }, async () => {
  const p = await load(), d = createDispatcher({ pkg: p, fs: folderFs(TB, "Tb2Pd2Sn") });
  if (d.has("plots_extent")) return;   // ported by then (L425)
  await assert.rejects(d.invoke("plots_extent", { id: 1, scale: "log", members: null }), (e) => e === missing("plots_extent", p.version()));
});
test("package: the ESRF run (101 .xye) described as data, opened as one set, the proposer asked when there is one", { skip: skip || (!existsSync(ESRF) && "the ESRF folder is not on this machine") }, async () => {
  const fs = folderFs(ESRF, "renamed"), d = createDispatcher({ pkg: await load(), fs });
  const ks = []; for (const p of fs.paths()) ks.push(await d.kindOf(p));
  assert.equal(ks.length, 101);
  assert.deepEqual(count(ks.map((k) => k.id)), { xye: 101 });
  const set = await d.invoke("plots_open_files", { paths: fs.paths(), root: fs.root });
  assert.equal(set.members.length, 101);
  console.log(`# ESRF: ${set.members.length} members, ${set.apart.length} apart, variable ${set.variable_name} (${set.variable_unit}); first ${set.members[0].label} from ${set.members[0].variable ? set.members[0].variable.source : "none"}`);
  const r = await d.propose(fs.paths());
  if (!d.has("plots_propose")) assert.equal(r, null);
  else console.log(`# ESRF proposals: ${JSON.stringify(r).slice(0, 400)}`);
});
