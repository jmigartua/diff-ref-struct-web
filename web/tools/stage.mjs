// diff-ref-struct Web · the public site staged (I10 step 5): the files the page loads (host-web/precache.json), the
// service worker and the core's package, laid out with the host at the root: app/ui/host-web/x → x, app/ui/x → x, the
// host's own references to "../" rewritten to "./". Nothing else of app/ui/ goes (no other module, no test, no tool).
// The version stamped in sw.js (so that a new one installs anew) and version.json; the third-party licences beside.
//   node web/tools/stage.mjs <out> [--web <version>]     (the package must be built: web/tools/build-wasm.sh)
import { readFileSync, writeFileSync, mkdirSync, rmSync, copyFileSync, existsSync, statSync, readdirSync } from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

/** Where a listed file goes in the site: a path relative to host-web/ ("../core.js", "host.js") → its path at the
    site's root ("core.js", "host.js"); null for one outside app/ui/. */
export function siteDest(rel) {
  const p = path.posix.normalize(path.posix.join("host-web", rel));   // relative to app/ui/
  if (p.startsWith("..")) return null;
  return p.startsWith("host-web/") ? p.slice("host-web/".length) : p;
}
/** A host file's text with its references to the shell's files ("../x") made the site root's ("./x"): imports, import()
    and the page's href/src. Other files are not touched (their relative imports hold as they are). */
export function rewriteHost(text) {
  return text.replace(/((?:\bfrom\s+|\bimport\s*\(\s*|\bimport\s+|\b(?:href|src)=))(["'])\.\.\//g, "$1$2./");
}
/** The staging plan from the precache list: { src (under app/ui/), dest (under the site), host (rewrite) }. A file of
    the host and one of the shell would collide at the root only by a name both have: refused. */
export function plan(list) {
  const out = [], seen = new Map();
  for (const rel of list) {
    const dest = siteDest(rel); if (!dest) throw new Error(`${rel}: outside app/ui/`);
    const src = path.posix.normalize(path.posix.join("host-web", rel));
    if (seen.has(dest) && seen.get(dest) !== src) throw new Error(`${dest}: both ${seen.get(dest)} and ${src}`);
    seen.set(dest, src);
    out.push({ src, dest, host: src.startsWith("host-web/") && /\.(js|html|json|webmanifest|css)$/.test(src) });
  }
  return out;
}

const here = path.dirname(fileURLToPath(import.meta.url));
const top = path.resolve(here, "../..");
const ui = path.join(top, "app/ui");

if (process.argv[1] && fileURLToPath(import.meta.url) === path.resolve(process.argv[1])) {
  const out = path.resolve(process.argv[2] || "");
  if (!process.argv[2]) { console.error("stage: node web/tools/stage.mjs <out> [--web <version>]"); process.exit(2); }
  const web = process.argv.includes("--web") ? process.argv[process.argv.indexOf("--web") + 1] : "local";
  // the tool's version: the first entry of the releases record (web/site/data/releases.yml), semantic from 1.0.0
  const releases = readFileSync(path.join(top, "web/site/data/releases.yml"), "utf8");
  const version = (releases.match(/^- version: "([^"]+)"/m) || [])[1] || "0.0.0";
  const manifest = JSON.parse(readFileSync(path.join(ui, "host-web/wasm/manifest.json"), "utf8"));
  if (!manifest.package) { console.error("stage: no package in app/ui/host-web/wasm/: web/tools/build-wasm.sh first"); process.exit(1); }
  const list = JSON.parse(readFileSync(path.join(ui, "host-web/precache.json"), "utf8")).files;
  const pkgFiles = [`wasm/${manifest.package}`, `wasm/${manifest.wasm || "drs_wasm_bg.wasm"}`];
  const steps = plan([...list, "sw.js", ...pkgFiles]);
  rmSync(out, { recursive: true, force: true }); mkdirSync(out, { recursive: true });
  const built = new Date().toISOString().slice(0, 10);
  const pkgVersion = (() => { try { return JSON.parse(readFileSync(path.join(ui, "host-web/wasm/package.json"), "utf8")).version; } catch (e) { return null; } })();
  let bytes = 0;
  for (const s of steps) {
    const from = path.join(ui, s.src), to = path.join(out, s.dest);
    if (!existsSync(from)) { console.error(`stage: ${s.src} is listed and missing`); process.exit(1); }
    mkdirSync(path.dirname(to), { recursive: true });
    if (s.dest === "sw.js") writeFileSync(to, readFileSync(from, "utf8").replace('const VERSION = "dev";', `const VERSION = ${JSON.stringify(`${version}+${web}`)};`));
    else if (s.dest === "version.json") writeFileSync(to, JSON.stringify({ version, build: web, web, built, core: pkgVersion }) + "\n");
    else if (s.dest === "wasm/manifest.json") writeFileSync(to, JSON.stringify({ package: manifest.package, wasm: manifest.wasm, bytes: manifest.bytes, core: pkgVersion }) + "\n");
    else if (s.host) writeFileSync(to, rewriteHost(readFileSync(from, "utf8")));
    else copyFileSync(from, to);
    bytes += statSync(to).size;
  }
  // the worker's list as the site lays it out
  writeFileSync(path.join(out, "precache.json"), JSON.stringify({ note: "The files the page loads, for sw.js (web/tools/precache.mjs, staged by web/tools/stage.mjs).", files: list.map(siteDest) }, null, 1) + "\n");
  // the third-party licences of what the site ships (web/licences/)
  const lic = [["dockview-core", "vendor/dockview/", "dockview-core.txt"], ["three.js", "vendor/three/", "three.txt"]].filter(([, dir]) => steps.some((s) => s.dest.startsWith(dir)));
  writeFileSync(path.join(out, "LICENCES.md"), `# Third-party code in diff-ref-struct Web\n\nThe site ships these libraries unchanged, each under its licence, reproduced here:\n\n${lic.map(([n, dir, f]) => `## ${n} (\`${dir}\`)\n\n\`\`\`\n${readFileSync(path.join(top, "web/licences", f), "utf8").trim()}\n\`\`\`\n`).join("\n")}\nThe rest is diff-ref-struct's own: the app's Plots module and shell (\`app/ui/\`), the browser host and the core compiled to WebAssembly (\`wasm/\`), from the private repository diff-ref-struct.\n`);
  writeFileSync(path.join(out, ".nojekyll"), "");
  // the example the core computed (Try the example), served beside the page: fetched by the page, not precached
  const ex = path.join(ui, "host-web/examples/ramp");
  if (existsSync(ex)) { mkdirSync(path.join(out, "examples/ramp"), { recursive: true }); for (const f of readdirSync(ex)) { copyFileSync(path.join(ex, f), path.join(out, "examples/ramp", f)); bytes += statSync(path.join(ex, f)).size; } }
  // the licences and the attribution travel with every copy of the built page (Apache-2.0 §4)
  for (const f of ["LICENSE", "LICENSE-docs", "NOTICE", "CITATION.cff"]) copyFileSync(path.join(top, "web/public", f), path.join(out, f));
  writeFileSync(path.join(out, "README.md"), `# diff-ref-struct Web ${version} (built page)\n\nThe built page and its guide (/guide/), published by the project's Action at each change; every run replaces this branch whole. The source is on the repository's main branch. Live at <https://jmigartua.github.io/diff-ref-struct-web/>.\n\ndiff-ref-struct Web ${version} · build ${web} · core ${pkgVersion} · ${built}. Code under the Apache License 2.0 (LICENSE), documents under CC BY 4.0 (LICENSE-docs); attribution in NOTICE; third-party licences in LICENCES.md; citation in CITATION.cff.\n`);
  console.log(`stage: ${steps.length} files and the example, ${(bytes / 1024).toFixed(0)} kB, diff-ref-struct Web ${version} (build ${web}), core ${pkgVersion} → ${out}`);
}
