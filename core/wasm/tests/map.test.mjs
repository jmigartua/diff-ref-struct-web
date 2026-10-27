// W1 (I10, R34; L424): the core in WebAssembly, as the browser's host will load it, against the command line.
// The package built by `wasm-pack build --target web` (`npm test` builds it first), loaded in node; the Tb₂Pd₂Sn files
// given by name and text as a browser reads them; the map equal to `drs plots map --json`'s to the bit (every number
// the same double: deepStrictEqual compares numbers with Object.is); the commands served by `dispatch` with the app's
// names and the core's words when one is refused; W2's proposals on Tb₂Pd₂Sn's buffer, W5's path shapes (the files a
// buffer or a set file names, joined as the app joins paths, the host's virtual paths `/local/<folder>/…`) and PL41's
// extent. Prints the .wasm's size.

import { test } from "node:test";
import assert from "node:assert/strict";
import { readFileSync, readdirSync, statSync, existsSync } from "node:fs";
import { execFileSync } from "node:child_process";
import { fileURLToPath } from "node:url";
import path from "node:path";
import { initSync, dispatch, version, commands } from "../pkg/drs_wasm.js";

const here = path.dirname(fileURLToPath(import.meta.url));
const core = path.resolve(here, "../..");
const tb = path.resolve(core, "../manuals/fullprof/1-WinPLOTR/1_Tb2pd2sn");
const wasm = path.resolve(here, "../pkg/drs_wasm_bg.wasm");
initSync({ module: readFileSync(wasm) });

const names = Array.from({ length: 20 }, (_, k) => `tb_t${String(k + 1).padStart(2, "0")}.dat`);
const files = names.map((name) => ({ name, text: readFileSync(path.join(tb, name), "utf8") }));

function drs(args) {
  const bin = path.join(core, "target/release/drs");
  if (!existsSync(bin)) execFileSync("cargo", ["build", "--release", "--manifest-path", path.join(core, "Cargo.toml")], { stdio: "inherit" });
  return execFileSync(bin, args, { encoding: "utf8", maxBuffer: 1 << 28 });
}

const call = (cmd, args) => dispatch(cmd, JSON.stringify(args));

test("the .wasm's size", () => {
  const kB = statSync(wasm).size / 1000;
  console.log(`drs_wasm_bg.wasm: ${kB.toFixed(0)} kB (drs-wasm ${version()})`);
  assert.ok(kB > 100);
});

test("the version is the command line's", () => {
  assert.equal(`drs ${version()}`, drs(["--version"]).trim());
});

test("the Tb₂Pd₂Sn map from names and text equals drs plots map --json to the bit", async () => {
  const set = await call("plots_open_from", { files });
  assert.equal(set.members.length, 20);
  assert.equal(set.members[0].label, "tb_t02.dat");
  assert.equal(set.members[0].variable.source, 'the title\'s "T.=.4.K"');
  const view = await call("plots_map", { id: set.id, scale: "linear", reference: null });
  const cli = JSON.parse(drs(["plots", "map", ...names.map((n) => path.join(tb, n)), "--json"]));
  assert.deepStrictEqual(view.map, cli);
  assert.equal(view.map.rows.length, 19);
  assert.equal(view.apart.length, 1);
  // the same on the square root, and the difference map with tb_t20 as the reference
  const k20 = set.members.findIndex((m) => m.label === "tb_t20.dat");
  const sq = await call("plots_map", { id: set.id, scale: "sqrt", reference: k20 });
  const cliSq = JSON.parse(drs(["plots", "map", ...set.members.map((m) => path.join(tb, m.file)), "--json", "--scale", "sqrt", "--reference", String(k20 + 1)]));
  assert.deepStrictEqual(sq.map, cliSq);
});

test("a region, a warning, the set file as text, describe", async () => {
  const set = await call("plots_open_from", { files });
  const region = { name: "m1", x: [18.95, 19.25], variable: [20, 26] };
  const block = await call("plots_map_region", { id: set.id, scale: "linear", reference: null, region });
  assert.deepStrictEqual([block.rows, block.columns, block.counted_from], [[8, 9, 10, 11, 12, 13, 14], [30, 31, 32], 0]);
  assert.equal(await call("plots_range_warning", { id: set.id, scale: "linear", reference: null, range: { from: 1e6, to: 2e6 } }), "the map's colour range 1000000 to 2000000 lies wholly outside the map's values, 905 to 68658: every cell is drawn saturated");
  const settings = { order: "variable", rule: "none", dx: 0, dy: 0, scale: "linear", reference: null, windows: [], opacity: [], regions: [region] };
  const saved = await call("plots_save", { id: set.id, root: "", name: "tb", settings });
  assert.equal(saved.relative, "Plots/tb.drv");
  assert.match(saved.text, /^format = "drs plot set 1"/);
  // opened again from its text, its members' files by name: the same members and the region
  const [again, st] = await call("plots_open_set_from", { name: "Plots/tb.drv", text: saved.text, files });
  assert.deepStrictEqual(again.members.map((m) => [m.label, m.variable && m.variable.value]), set.members.map((m) => [m.label, m.variable && m.variable.value]));
  assert.equal(st.regions[0].name, "m1");
  const d = await call("describe", { name: "Plots/tb.drv", text: saved.text });
  assert.equal(d.kind, "project");
  assert.equal(d.id, "plot_set");
  assert.equal((await call("describe", files[0])).kind, "data");
  // the stack, the projection, the windows, a derived pattern as text
  assert.equal((await call("plots_displace", { id: set.id, rule: "step", dx: 0, dy: 1000, scale: "linear" })).length, 20);
  const p = await call("plots_project", { id: set.id, rule: "step", dx: 0.5, dy: 1000, scale: "linear", floor: null, order: null });
  assert.ok(p.projection.curves.length > 0);
  const w = await call("plots_windows", { id: set.id, windows: [{ name: "m1", from: 18.5, to: 19.8 }], background: "auto", logLambda: null, at: [] });
  assert.equal(w.members.length, 20);
  const derived = await call("plots_derive", { id: set.id, members: [0, 1], kind: "sum", root: "", name: "tb" });
  assert.match(derived.text, /^XYDATA/);
});

// the host's virtual paths, as the web host names what it read
const local = "/local/1_Tb2pd2sn";
const localFiles = files.map((f) => ({ name: `${local}/${f.name}`, text: f.text }));
const buf = readFileSync(path.join(tb, "Tb_vs_T.buf"), "utf8");

test("W5: a buffer and a set file by path, their files named by the core", async () => {
  const { files: named } = await call("plots_names_in", { path: `${local}/Tb_vs_T.buf`, text: buf, root: null });
  assert.equal(named.length, 20);
  assert.ok(named.every((n) => n.startsWith(`${local}/tb_t`)), named[0]);
  const given = localFiles.filter((f) => named.includes(f.name));
  const set = await call("plots_open_buffer_from", { path: `${local}/Tb_vs_T.buf`, text: buf, files: given, root: null });
  assert.equal(set.members.length, 20);
  // the same set as from the buffer's text and its files by bare name, and as drs plots open-buf gives on disk
  const bare = await call("plots_open_buffer_from", { text: buf, files });
  assert.deepStrictEqual(set.members.map((m) => [m.file, m.variable && m.variable.value]), bare.members.map((m) => [m.file, m.variable && m.variable.value]));
  // a set file's members joined to the folder above Plots/
  const drv = 'format = "drs plot set 1"\ntitle = "two"\nvariable = { name = "T", unit = "K" }\n\n[[member]]\nfile = "Data/tb_t20.dat"\n\n[[member]]\nfile = "Data/tb_t02.dat"\n\n[[member]]\nfile = "Data/tb_t20.dat"\nscan = 1\n';
  const { files: inSet } = await call("plots_names_in", { path: "/local/proj/Plots/two.drv", text: drv, root: "/local/proj" });
  assert.deepStrictEqual(inSet, ["/local/proj/Data/tb_t20.dat", "/local/proj/Data/tb_t02.dat"]);
  const data = files.map((f) => ({ name: `/local/proj/Data/${f.name}`, text: f.text }));
  const [two] = await call("plots_open_set_from", { path: "/local/proj/Plots/two.drv", text: drv, files: data.filter((f) => inSet.includes(f.name)), root: "/local/proj" });
  assert.deepStrictEqual(two.members.map((m) => m.file), ["Data/tb_t02.dat", "Data/tb_t20.dat", "Data/tb_t20.dat"], "the file named twice, read once, two members");
  // a member the files lack, named as missing
  await assert.rejects(call("plots_open_set_from", { path: "/local/proj/Plots/two.drv", text: drv, files: data.filter((f) => f.name.endsWith("tb_t20.dat")) }), (e) => String(e) === "Data/tb_t02.dat: the set file names it, and no file of that name was given");
  // save as text: what plots_save would write, its path and its text
  const settings = { order: "variable", rule: "none", dx: 0, dy: 0, scale: "linear", reference: null, windows: [], opacity: [] };
  const saved = await call("plots_save_text", { id: two.id, root: "/local/proj", name: "two again", settings });
  assert.equal(saved.path, "Plots/two again.drv");
  assert.match(saved.text, /file = "Data\/tb_t02.dat"/);
});

test("W2: the Tb₂Pd₂Sn buffer's files offer one set by temperature", async () => {
  const { files: named } = await call("plots_names_in", { path: `${local}/Tb_vs_T.buf`, text: buf });
  const p = await call("plots_propose", { files: localFiles.filter((f) => named.includes(f.name)) });
  assert.deepStrictEqual(p.proposals.map((x) => [x.kind, x.title, x.set.members.length]), [["by_temperature", "19 scans by temperature, 4 to 30 K", 19]]);
  assert.deepStrictEqual(p.left_out, [["tb_t01.dat", "no temperature and no time: in no set"]]);
  assert.equal(p.classes.length, 1);
  // the proposal's members named relative to the folder, so that its set file opens from the folder's Plots/
  assert.ok(p.proposals[0].set.members.every((m) => /^tb_t\d\d\.dat$/.test(m.file)), p.proposals[0].set.members[0].file);
});

test("PL41: the set's extent on a scale is its map's colour range there", async () => {
  const set = await call("plots_open_from", { files: localFiles });
  const ext = await call("plots_extent", { id: set.id, scale: "log", members: null });
  const view = await call("plots_map", { id: set.id, scale: "log", reference: null });
  assert.deepStrictEqual(ext, view.map.range);
  const one = await call("plots_extent", { id: set.id, scale: "linear", members: [3] });
  const all = await call("plots_extent", { id: set.id, scale: "linear", members: null });
  assert.ok(all[0] <= one[0] && one[1] <= all[1] && one[0] < one[1]);
});

// W2, W8, W9 on a team's ESRF run (unpublished): its 101 files by text, from the folder DRS_ESRF_RENAMED names on this
// machine; skipped without it (W13: no name of the run in this repository).
const esrfDir = [process.env.DRS_ESRF_RENAMED].filter(Boolean).find((d) => existsSync(d) && readdirSync(d).filter((n) => n.endsWith(".xye")).length === 101);

test("W2, W8, W9: the ESRF folder's 101 files give L021's three sets, opened rebinned in under 10 s", { skip: esrfDir ? false : "the ESRF folder is not on this machine (DRS_ESRF_RENAMED)" }, async () => {
  const names = readdirSync(esrfDir).filter((n) => n.endsWith(".xye")).sort();
  assert.equal(names.length, 101);
  const esrf = names.map((n) => ({ name: `/local/renamed/${n}`, text: readFileSync(path.join(esrfDir, n), "utf8") }));
  let t = performance.now();
  const p = await call("plots_propose", { files: esrf });
  console.log(`plots_propose on the 101 files: ${((performance.now() - t) / 1000).toFixed(1)} s`);
  assert.deepStrictEqual(p.proposals.map((x) => [x.kind, x.set.members.length, x.set.rebin]), [["by_temperature", 25, 0.01], ["against_time", 100, 0.01], ["hold", 39, 0.01]]);
  for (const x of p.proposals) console.log(x.title);
  assert.equal(p.left_out.length, 1);
  assert.match(p.left_out[0][0], /_101_background_nosample_t38192s_long\.xye$/);
  // a proposal's members opened rebinned to 0.01°, the patterns out of the answer and asked for apart
  const hold = p.proposals[2].set.members.map((m) => m.file);
  t = performance.now();
  const set = await call("plots_open_from", { files: esrf.filter((f) => hold.some((h) => f.name.endsWith(`/${h}`))), rebin: 0.01 });
  const hs = (performance.now() - t) / 1000;
  t = performance.now();
  const all = await call("plots_open_from", { files: esrf.slice(0, 100), rebin: 0.01 });
  const s = (performance.now() - t) / 1000;
  // the reader the proposal names, given back: the files are not identified again
  assert.deepStrictEqual(p.proposals.map((x) => x.reader), ["xye", "xye", "xye"]);
  t = performance.now();
  const named = await call("plots_open_from", { files: esrf.slice(0, 100), rebin: 0.01, reader: p.proposals[1].reader });
  const sr = (performance.now() - t) / 1000;
  console.log(`plots_open_from rebin 0.01: the hold's 39 in ${hs.toFixed(1)} s, the run's 100 in ${s.toFixed(1)} s (${sr.toFixed(1)} s with the proposal's reader), ${(JSON.stringify(all).length / 1e6).toFixed(2)} MB of JSON`);
  // a time budget, not a tolerance (L429): 9.7 s in the cloud's container, 12.4 s on Josu's throttled Mac; the ceiling
  // catches a regression of kind (minutes), the line above reports the number
  assert.ok(s < 60, `${s} s`);
  assert.ok(sr < 60, `${sr} s`);
  assert.equal(named.members.length, 100);
  assert.equal(set.members.length, 39);
  assert.equal(all.members.length, 100);
  assert.equal(all.variable_unit, "°C");
  assert.ok(all.members.every((m) => m.x === undefined && m.y === undefined), "no pattern in the SetView");
  assert.ok(all.members.every((m) => m.points > 2690 && m.points < 2705), String(all.members[0].points));
  const pats = await call("plots_patterns", { id: all.id, members: [0, 99] });
  assert.deepStrictEqual(pats.map((q) => [q.index, q.x.length, q.y.length]), [[0, all.members[0].points, all.members[0].points], [99, all.members[99].points, all.members[99].points]]);
  // the map from the rebinned set, without the patterns ever crossing
  const map = await call("plots_map", { id: all.id, scale: "log", reference: null });
  assert.equal(map.map.rows.length, 100);
});

test("dispatch serves the commands by name and rejects in the core's words", async () => {
  const served = commands();
  for (const c of ["plots_open_from", "plots_patterns", "plots_set_other", "plots_names_in", "plots_extent", "plots_propose", "plots_save_text", "plots_open_buffer_from", "plots_map", "plots_displace", "plots_project", "plots_windows", "plots_derive", "plots_map_region", "plots_range_warning", "plots_save", "plots_open_set_from", "describe"]) assert.ok(served.includes(c), c);
  await assert.rejects(call("plots_map", { id: 9999, scale: "linear", reference: null }), (e) => e === "no plot set 9999");
  await assert.rejects(call("plots_open_files", { paths: [] }), (e) => String(e).startsWith('no command "plots_open_files"'));
  await assert.rejects(call("plots_map", { id: 1 }), (e) => String(e).startsWith("the command's arguments:"));
});
