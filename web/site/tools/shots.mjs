// diff-ref-struct Web · the site's screenshots, taken from the app itself (not drawn): headless Chrome driven over the
// DevTools protocol (Node 22's WebSocket, no library), the page served locally by web/tools/serve.sh with the wasm
// package built (web/tools/build-wasm.sh), the public example folder opened through the page's own folder input: the
// Tb₂Pd₂Sn patterns of FullProf's WinPLOTR tutorial (manuals/fullprof/1-WinPLOTR/1_Tb2pd2sn, the ILL's FullProf
// documentation page). Each scene sets the view through the page's address (hash.js), as a user's link would.
//   web/tools/serve.sh 8790 &   then   node [--experimental-websocket on Node 20] web/site/tools/shots.mjs [--port 8790] [--only map]
// Writes web/site/assets/figures/: view-<scene>-<theme>.webp (the plot areas at 3×), app-<scene>-<theme>.webp (the
// window, 1440 × 900 CSS px at 2×), app-<scene>.png (light, for the PDF, scenes marked pdf). The app follows the emulated theme (its
// "auto"). Nothing of the page is changed for a shot.
import { spawn } from "node:child_process";
import { writeFileSync, mkdtempSync, rmSync, existsSync } from "node:fs";
import { tmpdir } from "node:os";
import path from "node:path";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const top = path.resolve(here, "../../..");
const OUT = path.join(here, "../assets/figures");
const FOLDER = path.join(top, "manuals/fullprof/1-WinPLOTR/1_Tb2pd2sn");
const arg = (k, d) => { const i = process.argv.indexOf(k); return i > 0 ? process.argv[i + 1] : d; };
const PORT = Number(arg("--port", 8790)), ONLY = arg("--only", null);
const BASE = arg("--base", `http://127.0.0.1:${PORT}/host-web/`), OUTDIR = arg("--out", null);
const CHECK = process.argv.includes("--check");   // the release check: draw every scene, take nothing, fail on a blank view
const checks = [];
const THEMES = arg("--theme", CHECK ? "light" : "light,dark").split(",");
let THEME = "light";
const CHROME = process.env.CHROME || ["/Applications/_____Google Chrome.app/Contents/MacOS/Google Chrome", "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"].find(existsSync);
const W = 1440, H = 900;
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

/** The scenes: the address's view (hash.js's parameters), then what is clicked, then the shot. */
const SCENES = [
  { name: "folder", hash: null, after: [] },
  { name: "example", folder: false, hash: null, after: [["eval", "document.querySelector('[data-act=example]').click()"], ["until", "window.__drsWebState().tabs >= 1 && !!document.querySelector('.pl-areas canvas')"], ["wait", 4000]] },
  { name: "map", pdf: true, hash: "set=buf:Tb_vs_T.buf&view=map&scale=log&clip=3.2,4.5", after: [["members"], ["pin", 0.506, 0.33]] },
  { name: "overlay", hash: "set=buf:Tb_vs_T.buf&view=overlay", after: [["members"], ["wait", 6000]] },
  { name: "stack", hash: "set=buf:Tb_vs_T.buf&view=stack", after: [] },
  { name: "difference", hash: "set=buf:Tb_vs_T.buf&view=difference&scale=sqrt", after: [["wait", 2000]] },
  { name: "projection", hash: "set=buf:Tb_vs_T.buf&view=projection", after: [["button", "Δ"], ["form", ["step", "2500", "0.6"]], ["wait", 2500]] },
  { name: "series", hash: "set=buf:Tb_vs_T.buf&view=series", after: [["button", "+ window"], ["form", ["19.2°", "18.6", "19.8"]], ["button", "+ window"], ["form", ["30.8°", "30.4", "31.2"]], ["wait", 2500]] },
  { name: "surface", hash: "set=buf:Tb_vs_T.buf&view=surface&scale=log&clip=3.2,4.5", after: [] },
];

async function cdp(wsUrl) {
  const ws = new WebSocket(wsUrl); let id = 0; const wait = new Map(), events = [];
  await new Promise((r, j) => { ws.onopen = r; ws.onerror = j; });
  ws.onmessage = (m) => { const d = JSON.parse(m.data); if (d.id && wait.has(d.id)) { const { r, j } = wait.get(d.id); wait.delete(d.id); d.error ? j(new Error(JSON.stringify(d.error))) : r(d.result); } else events.push(d); };
  const send = (method, params = {}) => new Promise((r, j) => { const k = ++id; wait.set(k, { r, j }); ws.send(JSON.stringify({ id: k, method, params })); });
  return { send, events, close: () => ws.close() };
}
async function evaluate(c, expr) {
  const r = await c.send("Runtime.evaluate", { expression: expr, awaitPromise: true, returnByValue: true });
  if (r.exceptionDetails) throw new Error(r.exceptionDetails.text + " " + JSON.stringify(r.exceptionDetails.exception?.description));
  return r.result.value;
}
async function until(c, expr, what, ms = 60000) {
  const t0 = Date.now();
  while (Date.now() - t0 < ms) { try { if (await evaluate(c, expr)) return; } catch (e) { /* page not ready */ } await sleep(250); }
  throw new Error(`timed out waiting for ${what}`);
}

async function shoot(c, scene) {
  const url = `${BASE}${scene.hash ? "#folder=1_Tb2pd2sn&" + scene.hash : ""}`;
  await c.send("Page.navigate", { url: "about:blank" }); await sleep(200);
  await c.send("Page.navigate", { url });
  await until(c, "window.__drsWebState && window.__drsWebState().core", "the core in the page");
  if (scene.folder === false) {   // a scene of the page alone (Try the example): no folder of ours
    for (const [kind, a] of scene.after) { if (kind === "eval") await evaluate(c, a); if (kind === "wait") await sleep(a); if (kind === "until") await until(c, a, "the example", 120000); }
    await sleep(1500);
    const save = (name, data) => { const f = path.join(OUTDIR || OUT, name); writeFileSync(f, Buffer.from(data, "base64")); console.log(`shots: ${path.relative(top, f)}`); };
    if (!CHECK) save(`app-${scene.name}-${THEME}.webp`, (await c.send("Page.captureScreenshot", { format: "webp", quality: 90 })).data);
    else { const n = await evaluate(c, "window.__drsWebState().tabs"); checks.push([scene.name, n >= 1, { tabs: n }]); console.log(`check: ${scene.name}: ${n >= 1 ? "OK" : "FAIL"} tabs ${n}`); }
    return;
  }
  // the folder through the page's own folder input (Open Folder…'s fallback, <input webkitdirectory>)
  const { root } = await c.send("DOM.getDocument", { depth: 0 });
  const { nodeId } = await c.send("DOM.querySelector", { nodeId: root.nodeId, selector: "#in-folder" });
  await c.send("DOM.setFileInputFiles", { nodeId, files: [FOLDER] });
  await until(c, "window.__drsWebState().described >= 23 && window.__drsWebState().proposals === 'done'", "the folder described and proposed");
  if (scene.hash) {
    await until(c, "window.__drsWebState().tabs >= 1", "the set opened from the address");
    await sleep(2500);
  } else {
    await evaluate(c, "document.querySelector('#tree-by [data-by=\"kinds\"]').click()");
    await sleep(600);
  }
  for (const [kind, a, b] of scene.after) {
    if (kind === "wait") await sleep(a);
    if (kind === "form") await evaluate(c, `(() => { const d = document.getElementById('dlg-form'); if (!d.open) throw new Error('no form open'); const ins = [...d.querySelectorAll('#dlg-form-fields input, #dlg-form-fields select')]; ${JSON.stringify(a)}.forEach((v, i) => { const el = ins[i]; if (!el) return; el.value = v; el.dispatchEvent(new Event('input', { bubbles: true })); el.dispatchEvent(new Event('change', { bubbles: true })); }); document.getElementById('dlg-form-ok').click(); })()`);
    if (kind === "eval") console.log("eval:", JSON.stringify(await evaluate(c, a)).slice(0, 1500));
    if (kind === "dbl") await evaluate(c, `(() => { const r = [...document.querySelectorAll('#tree .wh-row')].find((e) => e.querySelector('.wh-name')?.textContent === ${JSON.stringify(a)}); r.dispatchEvent(new MouseEvent('dblclick', { bubbles: true })); })()`);
    if (kind === "tab") await evaluate(c, `(() => { const x = [...document.querySelectorAll('.wh-views .ctab')].find((e) => e.textContent.trim() === ${JSON.stringify(a)} && e.offsetParent); if (!x) throw new Error('no tab ${a}'); x.click(); })()`);
    if (kind === "members") await evaluate(c, "document.querySelector('#members').scrollIntoView({ block: 'start' })");
    if (kind === "button") await evaluate(c, `(() => { const x = [...document.querySelectorAll('button')].find((e) => e.textContent.trim() === ${JSON.stringify(a)} && e.offsetParent); if (!x) throw new Error('no button ${a}'); x.click(); })()`);
    if (kind === "pin") {
      const r = await evaluate(c, "(() => { const cv = [...document.querySelectorAll('.wh-views canvas')].sort((p, q) => q.width * q.height - p.width * p.height)[0]; const b = cv.getBoundingClientRect(); return { x: b.left, y: b.top, w: b.width, h: b.height }; })()");
      const x = r.x + r.w * a, y = r.y + r.h * b;
      await c.send("Input.dispatchMouseEvent", { type: "mouseMoved", x, y });
      await sleep(300);
      await c.send("Input.dispatchMouseEvent", { type: "mousePressed", x, y, button: "left", clickCount: 1 });
      await c.send("Input.dispatchMouseEvent", { type: "mouseReleased", x, y, button: "left", clickCount: 1 });
    }
    await sleep(1200);
  }
  await sleep(800);
  if (CHECK) {   // the release check: the view's main 2D canvas (its largest) drawn: ≥ 120 colours in a sample (healthy views 155–1047, the blank Overlay of 0.4.0 67); a WebGL view by its capture
    const r = await evaluate(c, "(() => { const a = [...document.querySelectorAll('.pl-areas')].find((e) => e.offsetParent); if (!a) return { none: true }; let n = 0, main = null; for (const cv of a.querySelectorAll('canvas')) { const g = cv.getContext('2d'); if (!g || !cv.width || !cv.height) continue; n += 1; if (!main || cv.width * cv.height > main.width * main.height) main = cv; } if (!main) return { canvases: 0, colours: 0 }; const d = main.getContext('2d').getImageData(0, 0, main.width, main.height).data, seen = new Set(); for (let i = 0; i < d.length; i += 4 * 31) seen.add((d[i] << 16) | (d[i + 1] << 8) | d[i + 2]); return { canvases: n, colours: seen.size }; })()");
    let ok = r.none ? !scene.hash : (r.canvases === 0 || r.colours >= 120);
    if (!ok && !r.none) {   // a WebGL view (the surface) is not readable as 2D: its area's capture, compressed, says whether it is drawn
      const box = await evaluate(c, "(() => { const b = [...document.querySelectorAll('.pl-areas')].find((e) => e.offsetParent).getBoundingClientRect(); return { x: b.left, y: b.top, width: b.width, height: b.height }; })()");
      const png = (await c.send("Page.captureScreenshot", { format: "png", clip: { ...box, scale: 1 } })).data;
      r.png = Math.round(png.length * 0.75); ok = r.png > 60000;
    }
    checks.push([scene.name, ok, r]);
    console.log(`check: ${scene.name} ${THEME}: ${ok ? "OK" : "FAIL"} ${JSON.stringify(r)}`);
    return;
  }
  // the view alone (the plot areas, without toolbar and note) at 3× in WebP, for the site's cards and its lightbox;
  // the whole window at 2× in WebP; and, in the light theme, the whole window in PNG for the PDF (LaTeX takes no WebP)
  const area = scene.hash ? await evaluate(c, "(() => { const a = [...document.querySelectorAll('.pl-areas')].find((e) => e.offsetParent); if (!a) return null; const b = a.getBoundingClientRect(); return { x: b.left, y: b.top, width: b.width, height: b.height }; })()") : null;
  const save = (name, data) => { const f = path.join(OUTDIR || OUT, name); writeFileSync(f, Buffer.from(data, "base64")); console.log(`shots: ${path.relative(top, f)}`); };
  if (area) save(`view-${scene.name}-${THEME}.webp`, (await c.send("Page.captureScreenshot", { format: "webp", quality: 90, clip: { ...area, scale: 1.5 } })).data);
  save(`app-${scene.name}-${THEME}.webp`, (await c.send("Page.captureScreenshot", { format: "webp", quality: 90 })).data);
  if (THEME !== "light" || !scene.pdf) return;
  const shot = await c.send("Page.captureScreenshot", { format: "png", captureBeyondViewport: false });
  if (process.argv.includes("--debug")) for (const e of c.events) if (e.method === "Runtime.exceptionThrown" || (e.method === "Runtime.consoleAPICalled" && /error|warn/.test(e.params.type))) console.log("page:", JSON.stringify(e.params).slice(0, 600));
  const file = path.join(OUTDIR || OUT, `app-${scene.name}.png`);
  writeFileSync(file, Buffer.from(shot.data, "base64"));
  console.log(`shots: ${path.relative(top, file)}`);
}

const profile = mkdtempSync(path.join(tmpdir(), "drs-shots-"));
const chrome = spawn(CHROME, [...(process.env.HEADFUL ? [] : ["--headless=new"]), `--remote-debugging-port=${PORT + 1000}`, `--user-data-dir=${profile}`, `--window-size=${W},${H}`, "--hide-scrollbars", "--force-color-profile=srgb", "--no-first-run", "--no-default-browser-check", "about:blank"], { stdio: "ignore" });
try {
  let list = null;
  for (let k = 0; k < 60 && !list; k++) { try { list = await (await fetch(`http://127.0.0.1:${PORT + 1000}/json/list`)).json(); } catch (e) { await sleep(250); } }
  const page = list.find((t) => t.type === "page");
  const c = await cdp(page.webSocketDebuggerUrl);
  await c.send("Page.enable"); await c.send("Runtime.enable"); await c.send("DOM.enable");
  await c.send("Emulation.setDeviceMetricsOverride", { width: W, height: H, deviceScaleFactor: 2, mobile: false });
  for (THEME of THEMES) {
    await c.send("Emulation.setEmulatedMedia", { features: [{ name: "prefers-color-scheme", value: THEME }] });
    for (const s of SCENES) if (!ONLY || s.name === ONLY) await shoot(c, s);
  }
  c.close();
} finally { chrome.kill(); await sleep(300); rmSync(profile, { recursive: true, force: true }); }
if (CHECK) { const bad = checks.filter(([, ok]) => !ok); console.log(`check: ${checks.length - bad.length} of ${checks.length} views drawn`); process.exit(bad.length ? 1 : 0); }
