// diff-ref-struct Web · steps 4 and 5's pure parts: the view in the address (hash.js), live mode's diff and growth
// (live.js), the offline list (web/tools/precache.mjs) and the public site's staging (web/tools/stage.mjs).
// node --test web/tests/
import test from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { encode, decode, missingFor } from "../../app/ui/host-web/hash.js";
import { diffListing, nothing, growth, liveLine, INTERVALS, DEFAULT_INTERVAL } from "../../app/ui/host-web/live.js";
import { references, closure, hostList } from "../tools/precache.mjs";
import { siteDest, rewriteHost, plan } from "../tools/stage.mjs";

const HOST = new URL("../../app/ui/host-web/", import.meta.url);

// ------------------------------------------------------------------ the hash
test("hash: a view written and read back, every part", () => {
  const st = { folder: "renamed", set: { by: "files", files: ["a 1.xye", "sub/b|c.xye"] }, view: "difference", scale: "sqrt", zoom: { x0: 18.5, x1: 21, y0: 900, y1: 1300 }, clip: { from: -2.5e3, to: 1e4 }, region: "the 110 peak" };
  const h = encode(st);
  assert.match(h, /^#folder=renamed&set=files&file=a\+1\.xye&file=sub%2Fb%7Cc\.xye&/);
  assert.deepEqual(decode(h), st);
  for (const set of [{ by: "drv", path: "Plots/my set.drv" }, { by: "buf", path: "Tb_vs_T.buf" }, { by: "proposal", kind: "against_time", class: null }, { by: "proposal", kind: "hold", class: 1 }]) {
    const back = decode(encode({ set, view: "map" }));
    assert.deepEqual(back.set, set); assert.equal(back.view, "map"); assert.equal(back.zoom, null);
  }
});
test("hash: what does not read is left out, a hash naming no set is none", () => {
  assert.equal(decode(""), null); assert.equal(decode("#view=map"), null); assert.equal(decode("#set=nonsense"), null); assert.equal(decode("#set=files"), null);
  const d = decode("#set=buf:x.buf&view=pie&scale=cubic&zoom=1,2,3&clip=5,1&region=");
  assert.deepEqual(d, { folder: null, set: { by: "buf", path: "x.buf" }, view: null, scale: null, zoom: null, clip: null, region: null });
  assert.equal(encode(null), ""); assert.equal(encode({ set: { by: "elsewhere" } }), "");
});
test("hash: the files a view needs that the folder lacks", () => {
  const has = (p) => ["a.xye", "Plots/s.drv"].includes(p);
  assert.deepEqual(missingFor({ by: "files", files: ["a.xye", "b.xye"] }, has), ["b.xye"]);
  assert.deepEqual(missingFor({ by: "drv", path: "Plots/s.drv" }, has), []);
  assert.deepEqual(missingFor({ by: "buf", path: "x.buf" }, has), ["x.buf"]);
  assert.deepEqual(missingFor({ by: "proposal", kind: "hold" }, has), []);
});

// ------------------------------------------------------------------ live mode
test("live: the difference between two readings, by name", () => {
  const prev = new Map([["s_2.xye", "1:1"], ["s_10.xye", "1:1"], ["gone.xye", "5:1"], ["edit.xye", "3:1"]]);
  const next = new Map([["s_2.xye", "1:1"], ["s_10.xye", "1:1"], ["s_11.xye", "2:2"], ["s_3.xye", "2:2"], ["edit.xye", "4:2"]]);
  assert.deepEqual(diffListing(prev, next), { added: ["s_3.xye", "s_11.xye"], vanished: ["gone.xye"], changed: ["edit.xye"] });
  assert.ok(nothing(diffListing(prev, prev)));
  assert.ok(INTERVALS.includes(DEFAULT_INTERVAL));
});
test("live: a set grows with the files the core now proposes for its kind and class, never loses one", () => {
  const props = [{ kind: "against_time", class: null, paths: ["/l/a", "/l/b", "/l/c"] }, { kind: "hold", class: 0, paths: ["/l/b", "/l/c"] }];
  assert.deepEqual(growth(props, { kind: "against_time", class: null }, new Set(["/l/a", "/l/x"])), ["/l/b", "/l/c"]);
  assert.deepEqual(growth(props, { kind: "hold", class: 0 }, new Set(["/l/b", "/l/c"])), []);
  assert.deepEqual(growth(props, { kind: "by_temperature", class: null }, new Set()), [], "a kind the core no longer proposes: nothing added");
  const at = new Date(2026, 9, 25, 14, 3, 7);
  assert.equal(liveLine({ added: [], vanished: [], changed: [] }, 0, at), null);
  assert.equal(liveLine({ added: ["a"], vanished: ["b"], changed: [] }, 1, at), "14:03:07: 1 new file, 1 vanished (still in the sets opened: nothing is removed from a set by itself), 1 scan added to the live set.");
});

// ------------------------------------------------------------------ the offline list
test("precache: references by kind of file, local ones only", () => {
  assert.deepEqual(references("a/index.html", `<link href="x.css"><script src="../y.js"></script><a href="https://h/z.js"><img src="data:a"><a href="#top">`), ["x.css", "../y.js"]);
  assert.deepEqual(references("m.js", `import a from "./a.js";\nimport "./side.js";\nexport { b } from "../b.js";\nconst c = await import("./c.js");\nconst d = import(name);`), ["./a.js", "./side.js", "../b.js", "./c.js"]);
  assert.deepEqual(references("m.webmanifest", `{"icons":[{"src":"icon.svg"}]}`), ["icon.svg"]);
});
test("precache: the closure walks imports once each and names what is missing", () => {
  const fs = { "w/index.html": `<script src="h.js"></script>`, "w/h.js": `import "../c.js";\nimport("./lazy.js");`, "c.js": `import "./w/h.js";\nimport "./gone.js";`, "w/lazy.js": "" };
  const r = closure(["w/index.html"], (p) => (p in fs ? fs[p] : null));
  assert.deepEqual(r.files, ["w/index.html", "w/h.js", "c.js", "w/lazy.js"]);
  assert.deepEqual(r.missing, ["gone.js"]);
});
test("precache: the committed list is the page's as it is now, every file there, the stand-in gone", () => {
  const l = hostList();
  assert.deepEqual(l.missing, []);
  const committed = JSON.parse(readFileSync(new URL("precache.json", HOST), "utf8"));
  assert.deepEqual(committed.files, l.files.map((f) => f.path), "node web/tools/precache.mjs");
  for (const f of ["index.html", "host.js", "../core.js", "../modules/plots/index.js", "../vendor/dockview/dockview-core.esm.js", "../surface3d.js", "../vendor/three/three.core.js", "manifest.webmanifest", "icon.svg", "wasm/manifest.json"]) assert.ok(committed.files.includes(f), f);
  assert.ok(!committed.files.some((f) => /tauri-shim|tests\/|tools\//.test(f)));
});

// ------------------------------------------------------------------ the staging
test("stage: the host at the site's root, the shell's files beside it, nothing outside app/ui/", () => {
  assert.equal(siteDest("host.js"), "host.js");
  assert.equal(siteDest("../core.js"), "core.js");
  assert.equal(siteDest("../modules/plots/index.js"), "modules/plots/index.js");
  assert.equal(siteDest("wasm/drs_wasm.js"), "wasm/drs_wasm.js");
  assert.equal(siteDest("../../core/src/lib.rs"), null);
  const p = plan(["index.html", "../core.js", "../vendor/three/three.module.js"]);
  assert.deepEqual(p.map((s) => [s.src, s.dest, s.host]), [["host-web/index.html", "index.html", true], ["core.js", "core.js", false], ["vendor/three/three.module.js", "vendor/three/three.module.js", false]]);
  assert.throws(() => plan(["core.js", "../core.js"]), /both/);
  assert.throws(() => plan(["../../x.js"]), /outside/);
});
test("stage: the host's references to ../ made the root's, its other text untouched", () => {
  const t = `import { a } from "../core.js";\nimport "../side.js";\nconst m = await import('../lazy.js');\n<link href="../style.css"><script src="../x.js">\nconst parent = "../not-a-reference"; // "../c.xy" in a comment`;
  assert.equal(rewriteHost(t), `import { a } from "./core.js";\nimport "./side.js";\nconst m = await import('./lazy.js');\n<link href="./style.css"><script src="./x.js">\nconst parent = "../not-a-reference"; // "../c.xy" in a comment`);
  const host = readFileSync(new URL("host.js", HOST), "utf8");
  assert.ok(!/\.\.\//.test(rewriteHost(host).split("\n").filter((l) => /^\s*import\b/.test(l)).join("\n")));
});
