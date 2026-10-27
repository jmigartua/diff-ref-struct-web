// diff-ref-struct Web · the public source (the open-tool kit, §6 of web/PROPOSAL-open-tool.md): the files the page is
// built from, taken from the project's private repository and laid out for the public repository's main branch.
//   - core/ (the engine and its WebAssembly crate), every file git tracks;
//   - app/ui/: the files the page loads (host-web/precache.json, the same list the service worker keeps), the browser
//     host app/ui/host-web/ whole, and what builds the vendored libraries (app/package.json, its lock, tools/vendor.js);
//   - web/: the page's tools, tests, the guide's sources and the public manuals;
//   - web/public/*: the repository's own files at its root (README, LICENSE, LICENSE-docs, NOTICE, CITATION.cff,
//     CONTRIBUTING, CODE_OF_CONDUCT, .zenodo.json, the issue forms).
// Left out: web/private/ (the guard's words), every team's edition of a manual (a manual marked
// `<!-- edition: team -->`, and web/manual/team/), anything git does not track. Then the guard: no file may hold a word
// of web/private/guard.txt; a file that does stops the run, named, and nothing is written.
//   node web/tools/source.mjs <out> [--audit]      --audit: list what would go, by folder, and what the guard finds
import { execFileSync } from "node:child_process";
import { readFileSync, writeFileSync, mkdirSync, rmSync, copyFileSync, existsSync, statSync } from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const top = path.resolve(here, "../..");
const tracked = (dir) => execFileSync("git", ["ls-files", "-z", dir], { cwd: top }).toString("utf8").split("\0").filter(Boolean);

/** The files of the public source, as paths in the project's repository → paths in the public one. */
export function selection() {
  const pick = new Map();
  for (const f of tracked("core")) pick.set(f, f);
  const list = JSON.parse(readFileSync(path.join(top, "app/ui/host-web/precache.json"), "utf8")).files;
  for (const rel of list) { const f = path.posix.normalize(path.posix.join("app/ui/host-web", rel)); if (existsSync(path.join(top, f))) pick.set(f, f); }
  for (const f of tracked("app/ui/host-web")) pick.set(f, f);
  for (const f of ["app/package.json", "app/package-lock.json", "app/tools/vendor.js"]) if (existsSync(path.join(top, f))) pick.set(f, f);
  for (const f of tracked("web")) {
    if (f.startsWith("web/private/") || f.startsWith("web/manual/team/") || f.startsWith("web/public/")) continue;
    if (f.startsWith("web/manual/") && f.endsWith(".md") && /<!--\s*edition:\s*team/.test(readFileSync(path.join(top, f), "utf8").slice(0, 2000))) continue;
    pick.set(f, f);
  }
  for (const f of tracked("web/public")) pick.set(f, f.slice("web/public/".length));
  return pick;
}

/** The guard's words (web/private/guard.txt), and the files of the selection that hold one. */
export function guard(pick) {
  const g = path.join(top, "web/private/guard.txt");
  const words = existsSync(g) ? readFileSync(g, "utf8").split("\n").map((w) => w.trim()).filter((w) => w && !w.startsWith("#")) : [];
  const hits = [];
  for (const f of pick.keys()) {
    const buf = readFileSync(path.join(top, f));
    if (buf.includes(0)) continue;   // a binary file (an image, the .wasm): not text to search
    const low = buf.toString("utf8").toLowerCase();
    const found = words.filter((w) => low.includes(w.toLowerCase()));
    if (found.length) hits.push([f, found]);
  }
  return { words: words.length, hits };
}

if (process.argv[1] && fileURLToPath(import.meta.url) === path.resolve(process.argv[1])) {
  const out = process.argv[2] && !process.argv[2].startsWith("--") ? path.resolve(process.argv[2]) : null;
  const audit = process.argv.includes("--audit");
  const pick = selection();
  const g = guard(pick);
  if (audit) {
    const by = new Map();
    for (const f of pick.keys()) { const k = f.split("/").slice(0, 2).join("/"); const n = by.get(k) || [0, 0]; n[0] += 1; n[1] += statSync(path.join(top, f)).size; by.set(k, n); }
    for (const [k, [n, b]] of [...by].sort()) console.log(`${k.padEnd(28)} ${String(n).padStart(5)} files ${(b / 1024).toFixed(0).padStart(7)} kB`);
    console.log(`guard: ${g.words} words, ${g.hits.length} files hold one`);
    for (const [f, w] of g.hits) console.log(`  ${f}: ${w.join(", ")}`);
  }
  if (g.hits.length) { console.error(`source: ${g.hits.length} files hold a guarded word; nothing written (web/private/guard.txt)`); process.exit(1); }
  if (!out) process.exit(0);
  rmSync(out, { recursive: true, force: true }); mkdirSync(out, { recursive: true });
  for (const [from, to] of pick) { mkdirSync(path.dirname(path.join(out, to)), { recursive: true }); copyFileSync(path.join(top, from), path.join(out, to)); }
  console.log(`source: ${pick.size} files → ${out}`);
}
