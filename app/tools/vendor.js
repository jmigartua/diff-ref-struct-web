// Copies the vendored libraries from node_modules into ui/vendor/ (generated, not committed), so that the window
// loads them from the bundle: no network at run time (CSP: 'self'). KaTeX for the equations of the parameter cards;
// dockview-core for the dock (DESIGN.md §9, answer 6).
"use strict";
const fs = require("fs");
const path = require("path");

const ui = path.join(__dirname, "..", "ui", "vendor");
const nm = path.join(__dirname, "..", "node_modules");
const version = (pkg) => JSON.parse(fs.readFileSync(path.join(nm, pkg, "package.json"))).version;

// KaTeX
const katex = path.join(nm, "katex", "dist");
if (!fs.existsSync(katex)) { console.error("vendor.js: katex is not installed (npm install)"); process.exit(1); }
const kdest = path.join(ui, "katex");
fs.rmSync(kdest, { recursive: true, force: true });
fs.mkdirSync(path.join(kdest, "fonts"), { recursive: true });
for (const f of ["katex.min.js", "katex.min.css"]) fs.copyFileSync(path.join(katex, f), path.join(kdest, f));
for (const f of fs.readdirSync(path.join(katex, "fonts"))) if (f.endsWith(".woff2")) fs.copyFileSync(path.join(katex, "fonts", f), path.join(kdest, "fonts", f));
fs.writeFileSync(path.join(kdest, "VERSION"), version("katex") + "\n");
console.log(`ui/vendor/katex: KaTeX ${version("katex")}`);

// dockview-core
const dv = path.join(nm, "dockview-core", "dist");
if (!fs.existsSync(dv)) { console.error("vendor.js: dockview-core is not installed (npm install)"); process.exit(1); }
const ddest = path.join(ui, "dockview");
fs.rmSync(ddest, { recursive: true, force: true });
fs.mkdirSync(ddest, { recursive: true });
fs.copyFileSync(path.join(dv, "dockview-core.esm.js"), path.join(ddest, "dockview-core.esm.js"));
fs.copyFileSync(path.join(dv, "styles", "dockview.css"), path.join(ddest, "dockview.css"));
fs.writeFileSync(path.join(ddest, "VERSION"), version("dockview-core") + "\n");
console.log(`ui/vendor/dockview: dockview-core ${version("dockview-core")}`);

// js-yaml: the tour player reads tutorial/script.yml in the window (TU1)
const jy = path.join(nm, "js-yaml", "dist", "js-yaml.mjs");
if (!fs.existsSync(jy)) { console.error("vendor.js: js-yaml is not installed (npm install)"); process.exit(1); }
const jydest = path.join(ui, "js-yaml");
fs.rmSync(jydest, { recursive: true, force: true });
fs.mkdirSync(jydest, { recursive: true });
fs.copyFileSync(jy, path.join(jydest, "js-yaml.mjs"));
fs.writeFileSync(path.join(jydest, "VERSION"), version("js-yaml") + "\n");
console.log(`ui/vendor/js-yaml: js-yaml ${version("js-yaml")}`);

// The tutorial's script (TU1), bundled for Help ▸ Tour of a Refinement; the tutorial stream owns it
const tourSrc = path.join(__dirname, "..", "..", "tutorial", "script.yml");
const tourDest = path.join(__dirname, "..", "ui", "tour-script.yml");
if (fs.existsSync(tourSrc)) { fs.copyFileSync(tourSrc, tourDest); console.log("ui/tour-script.yml: the tutorial's script"); } else { fs.rmSync(tourDest, { force: true }); console.log("ui/tour-script.yml: no tutorial/script.yml; Help ▸ Tour says so"); }

// CodeMirror 6, bundled from ui/editor/cm.src.js into one ES module (esbuild)
const { buildSync } = require("esbuild");
const cmdest = path.join(ui, "codemirror");
fs.rmSync(cmdest, { recursive: true, force: true });
fs.mkdirSync(cmdest, { recursive: true });
buildSync({ entryPoints: [path.join(__dirname, "..", "ui", "editor", "cm.src.js")], bundle: true, format: "esm", outfile: path.join(cmdest, "codemirror.js"), minify: true, legalComments: "none", logLevel: "error" });
fs.writeFileSync(path.join(cmdest, "VERSION"), `codemirror ${version("codemirror")}, @codemirror/view ${version("@codemirror/view")}, @codemirror/lint ${version("@codemirror/lint")}, @codemirror/merge ${version("@codemirror/merge")}, legacy-modes ${version("@codemirror/legacy-modes")}\n`);
console.log(`ui/vendor/codemirror: ${fs.readFileSync(path.join(cmdest, "VERSION"), "utf8").trim()} (${Math.round(fs.statSync(path.join(cmdest, "codemirror.js")).size / 1024)} kB)`);

// three.js (the Explore module's 3D view, explore/DESIGN.md §4, decided by Josu on 11 October 2026): the ESM build and
// the core it imports, copied as they are; the module imports ui/vendor/three/three.module.js.
const three = path.join(nm, "three", "build");
if (!fs.existsSync(three)) { console.error("vendor.js: three is not installed (npm install)"); process.exit(1); }
const tdest = path.join(ui, "three");
fs.rmSync(tdest, { recursive: true, force: true });
fs.mkdirSync(tdest, { recursive: true });
for (const f of ["three.module.js", "three.core.js"]) fs.copyFileSync(path.join(three, f), path.join(tdest, f));
fs.writeFileSync(path.join(tdest, "VERSION"), version("three") + "\n");
console.log(`ui/vendor/three: three ${version("three")}`);

// canvas2svg (MIT): the recorder that turns a plot's drawing into SVG for File ▸ Export ▸ Plot (L134); a classic
// script that defines window.C2S, loaded by index.html before the app.
const c2s = path.join(nm, "canvas2svg", "canvas2svg.js");
if (!fs.existsSync(c2s)) { console.error("vendor.js: canvas2svg is not installed (npm install)"); process.exit(1); }
const cdest = path.join(ui, "canvas2svg");
fs.rmSync(cdest, { recursive: true, force: true });
fs.mkdirSync(cdest, { recursive: true });
fs.copyFileSync(c2s, path.join(cdest, "canvas2svg.js"));
fs.copyFileSync(path.join(nm, "canvas2svg", "LICENSE"), path.join(cdest, "LICENSE"));
fs.writeFileSync(path.join(cdest, "VERSION"), version("canvas2svg") + "\n");
console.log(`ui/vendor/canvas2svg: canvas2svg ${version("canvas2svg")}`);
