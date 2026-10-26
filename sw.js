/* diff-ref-struct Web · the service worker (I10 step 4; web/DESIGN.md): the page, the shell's files, the module and
   the core's package kept by the browser, so that the app opens again without the network. It keeps code only: the
   user's files are never fetched (they are read from the disk by the page) and never cached. The list is
   precache.json (web/tools/precache.mjs) and the package wasm/manifest.json names, and nothing else is ever put in the
   cache (a request not listed passes to the network untouched). VERSION is stamped by
   web/tools/stage.sh for the published site, so that a new version installs anew; "dev" (a local checkout) asks the
   network first and keeps what it gets for offline. */
const VERSION = "a1c308d";
const CACHE = "drs-web-" + VERSION;

async function listed() {
  const list = await (await fetch("precache.json", { cache: "no-cache" })).json();
  const urls = ["./", ...list.files];
  try {
    const m = await (await fetch("wasm/manifest.json", { cache: "no-cache" })).json();
    if (m.package) urls.push("wasm/" + m.package, "wasm/" + (m.wasm || "drs_wasm_bg.wasm"));
  } catch (e) { /* no package: the stub */ }
  return urls;
}
self.addEventListener("install", (e) => e.waitUntil((async () => {
  await caches.delete(CACHE);   // a worker installs onto an empty cache: what it keeps is its list, nothing older
  const c = await caches.open(CACHE);
  await c.addAll((await listed()).map((u) => new Request(u, { cache: "reload" })));
  await self.skipWaiting();
})()));
self.addEventListener("activate", (e) => e.waitUntil((async () => {
  for (const k of await caches.keys()) if (k.startsWith("drs-web-") && k !== CACHE) await caches.delete(k);
  await self.clients.claim();
})()));
self.addEventListener("fetch", (e) => {
  const r = e.request;
  if (r.method !== "GET" || new URL(r.url).origin !== self.location.origin) return;
  e.respondWith((async () => {
    const c = await caches.open(CACHE);
    const kept = () => c.match(r, { ignoreSearch: true }).then((m) => m || (r.mode === "navigate" ? c.match("./") : undefined));
    const listed = await c.match(r, { ignoreSearch: true });   // only the page's own files are kept: never anything else
    if (VERSION !== "dev" && listed) return listed;
    try {
      const got = await fetch(r);
      if (got.ok && listed) c.put(r, got.clone());   // dev: the listed file kept fresh
      return got;
    } catch (err) {
      const m = await kept(); if (m) return m;
      throw err;
    }
  })());
});
/** What the page's offline line says: the version, how many files and bytes are kept. */
self.addEventListener("message", (e) => {
  if (!e.data || e.data.ask !== "kept") return;
  e.waitUntil((async () => {
    const c = await caches.open(CACHE), keys = await c.keys();
    let bytes = 0; for (const k of keys) { const m = await c.match(k); if (m) bytes += (await m.clone().arrayBuffer()).byteLength; }
    e.source.postMessage({ kept: { version: VERSION, files: keys.length, bytes } });
  })());
});
