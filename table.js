/* diff-ref-struct · the shared table component (APP_CONVENTIONS §8, C8.2): sortable and resizable columns, a choice
   of columns, find, multiple selection (click, ⇧-click, ⌘-click), Copy as tab-separated text, and a selection synced
   with the plot where both show the same data. Rows are drawn only where visible (a scan has thousands). The first
   column stays in view when the table scrolls sideways, and a module with many variables lists them as rows, never
   as columns (§8, L119). */
import { el, esc, state } from "./core.js";

/** The columns a table's user chose to show or hide, by the table's key (§8: the choice is remembered, app 0.44.0):
    kept in the app's settings, `tableColumns[key][column id]` true or false; a column never chosen keeps its default. */
function chosenColumns(key) { const all = (state.settings && state.settings.tableColumns) || {}; return (key && all[key]) || {}; }
function rememberColumn(key, id, shown) {
  if (!key || !state.settings) return;
  const all = { ...(state.settings.tableColumns || {}) };
  all[key] = { ...(all[key] || {}), [id]: shown };
  state.settings.tableColumns = all;
  try { if (window.__saveSettings) window.__saveSettings(["tableColumns"]); } catch (e) { /* outside the app */ }
}

/** columns: [{ id, title, get(row, i) -> value, fmt?(value, row) -> text, width?, shown?, num?, tip? }]; rows: an array.
    opts: { onSelect(indices), rowClass?(i), labels?: true, key?, fitRows? } — `labels`: the first column names the rows (a module's
    variables, E8a), set left and in the monospaced face; it stays in view when the table scrolls sideways, its header
    too, whatever the table (§8). `key`: the table's name in the settings, under which the columns chosen are kept. `fitRows`: a short table (an inspector's,
    a sheet's) as tall as its rows, up to that many, instead of filling its container. */
export function createTable(container, columns, opts = {}) {
  const chosen = chosenColumns(opts.key);
  const st = { columns: columns.map((c) => ({ ...c, shown: c.id in chosen ? !!chosen[c.id] : c.shown !== false })), n: 0, order: null, sortBy: null, sortDir: 1, sel: new Set(), anchor: null, filter: "", rowH: 22 };
  container.classList.add("dtable");
  if (opts.labels) container.classList.add("dt-labels");
  container.innerHTML = `<div class="dt-head"><input class="dt-find" type="search" placeholder="Find…" aria-label="Find in the table"><span class="dt-count"></span><button class="dt-cols small" title="Choose the columns">Columns</button><button class="dt-copy small" title="Copy the selected rows (or all) as tab-separated text (⌘C)">Copy</button></div><div class="dt-scroll"><div class="dt-spacer"></div><table class="dt"><thead></thead><tbody></tbody></table></div>`;
  const scroll = container.querySelector(".dt-scroll"), spacer = container.querySelector(".dt-spacer"), thead = container.querySelector("thead"), tbody = container.querySelector("tbody"), table = container.querySelector("table");
  let rows = [];

  const visible = () => st.columns.filter((c) => c.shown);
  function value(c, i) { return c.get(rows[i], i); }
  function text(c, i) { const v = value(c, i); return c.fmt ? c.fmt(v, rows[i]) : v == null ? "" : String(v); }

  function computeOrder() {
    let idx = [...Array(rows.length).keys()];
    if (st.filter) { const f = st.filter.toLowerCase(); idx = idx.filter((i) => visible().some((c) => text(c, i).toLowerCase().includes(f))); }
    if (st.sortBy) { const c = st.columns.find((q) => q.id === st.sortBy); idx.sort((a, b) => { const va = value(c, a), vb = value(c, b); return (va == null ? 1 : vb == null ? -1 : va < vb ? -1 : va > vb ? 1 : 0) * st.sortDir; }); }
    st.order = idx;
    container.querySelector(".dt-count").textContent = `${idx.length}${idx.length !== rows.length ? ` of ${rows.length}` : ""} rows`;
  }
  function renderHead() {
    thead.innerHTML = "<tr>" + visible().map((c) => `<th data-id="${c.id}"${c.tip ? ` title="${esc(c.tip)}"` : ""} style="${c.width ? `width:${c.width}px` : ""}">${esc(c.title)}${st.sortBy === c.id ? (st.sortDir > 0 ? " ▲" : " ▼") : ""}<span class="dt-resize"></span></th>`).join("") + "</tr>";
    thead.querySelectorAll("th").forEach((th) => {
      th.onclick = (e) => { if (e.target.classList.contains("dt-resize")) return; const id = th.dataset.id; if (st.sortBy === id) { if (st.sortDir < 0) { st.sortBy = null; st.sortDir = 1; } else st.sortDir = -1; } else { st.sortBy = id; st.sortDir = 1; } computeOrder(); renderHead(); renderBody(); };
      const r = th.querySelector(".dt-resize");
      r.onmousedown = (e) => { e.preventDefault(); const c = st.columns.find((q) => q.id === th.dataset.id); const x0 = e.clientX, w0 = th.offsetWidth; const mv = (ev) => { c.width = Math.max(40, w0 + ev.clientX - x0); th.style.width = c.width + "px"; }; const up = () => { window.removeEventListener("mousemove", mv); window.removeEventListener("mouseup", up); }; window.addEventListener("mousemove", mv); window.addEventListener("mouseup", up); };
    });
  }
  function renderBody() {
    const n = st.order.length;
    spacer.style.height = n * st.rowH + "px";
    const first = Math.max(0, Math.floor(scroll.scrollTop / st.rowH) - 5), last = Math.min(n, Math.ceil((scroll.scrollTop + scroll.clientHeight) / st.rowH) + 5);
    table.style.transform = `translateY(${first * st.rowH}px)`;
    const cols = visible();
    let html = "";
    for (let k = first; k < last; k++) {
      const i = st.order[k];
      html += `<tr data-i="${i}" class="${st.sel.has(i) ? "selected" : ""}${opts.rowClass ? " " + opts.rowClass(i) : ""}">` + cols.map((c) => `<td class="${c.num !== false ? "num" : ""}">${esc(text(c, i))}</td>`).join("") + "</tr>";
    }
    tbody.innerHTML = html;
  }
  scroll.addEventListener("scroll", renderBody);
  tbody.addEventListener("click", (e) => {
    const tr = e.target.closest("tr"); if (!tr) return;
    const i = parseInt(tr.dataset.i);
    if (e.shiftKey && st.anchor != null) { const a = st.order.indexOf(st.anchor), b = st.order.indexOf(i); for (let k = Math.min(a, b); k <= Math.max(a, b); k++) st.sel.add(st.order[k]); }
    else if (e.metaKey) { if (st.sel.has(i)) st.sel.delete(i); else st.sel.add(i); st.anchor = i; }
    else { st.sel.clear(); st.sel.add(i); st.anchor = i; }
    renderBody();
    if (opts.onSelect) opts.onSelect([...st.sel]);
  });
  container.querySelector(".dt-find").oninput = (e) => { st.filter = e.target.value.trim(); computeOrder(); renderBody(); };
  container.querySelector(".dt-copy").onclick = () => copy();
  container.querySelector(".dt-cols").onclick = (e) => {
    const old = document.querySelector(".ctx"); if (old) old.remove();
    const m = el("div", "ctx");
    for (const c of st.columns) { const b = el("button", null, `${c.shown ? "✓ " : "　"}${esc(c.title)}`); b.onclick = () => { c.shown = !c.shown; rememberColumn(opts.key, c.id, c.shown); m.remove(); renderHead(); computeOrder(); renderBody(); }; m.appendChild(b); }
    const r = e.target.getBoundingClientRect(); m.style.left = r.left + "px"; m.style.top = r.bottom + 4 + "px"; document.body.appendChild(m);
    setTimeout(() => document.addEventListener("click", function h(ev) { if (!m.contains(ev.target)) { m.remove(); document.removeEventListener("click", h); } }), 0);
  };
  container.addEventListener("keydown", (e) => {
    if ((e.metaKey || e.ctrlKey) && e.key.toLowerCase() === "c") { e.preventDefault(); copy(); }
    if (e.key === "ArrowDown" || e.key === "ArrowUp") {
      e.preventDefault();
      const cur = st.anchor != null ? st.order.indexOf(st.anchor) : -1;
      const k = Math.max(0, Math.min(st.order.length - 1, cur + (e.key === "ArrowDown" ? 1 : -1)));
      const i = st.order[k]; if (i == null) return;
      if (!e.shiftKey) st.sel.clear(); st.sel.add(i); st.anchor = i; renderBody();
      const y = k * st.rowH; if (y < scroll.scrollTop || y > scroll.scrollTop + scroll.clientHeight - st.rowH) scroll.scrollTop = y - scroll.clientHeight / 2;
      if (opts.onSelect) opts.onSelect([...st.sel]);
    }
  });
  container.tabIndex = 0;

  async function copy() {
    const cols = visible();
    const idx = st.sel.size ? st.order.filter((i) => st.sel.has(i)) : st.order;
    const tsv = [cols.map((c) => c.title).join("\t"), ...idx.map((i) => cols.map((c) => text(c, i)).join("\t"))].join("\n");
    // in the app the pasteboard's own (a copy from a menu or ⌘K is no user's action to the web view's clipboard)
    // W6 (PL13q): through the host (the page's, else Tauri's), the browser's clipboard when there is none
    const h = window.__DRS_HOST__ || (window.__TAURI__ ? { invoke: (c, a) => window.__TAURI__.core.invoke(c, a) } : null);
    try { if (h && typeof h.invoke === "function") await h.invoke("clipboard_text", { text: tsv }); else await navigator.clipboard.writeText(tsv); } catch (e) { /* the clipboard may be refused outside the app */ }
    return idx.length;
  }
  function setRows(r) {
    rows = r; st.sel.clear(); st.anchor = null;
    if (opts.fitRows) container.style.height = `${Math.min(Math.max(rows.length, 1), opts.fitRows) * st.rowH + 64}px`;
    computeOrder(); renderHead(); renderBody();
  }
  function setSelection(indices, scrollTo = true) {
    st.sel = new Set(indices); st.anchor = indices.length ? indices[0] : null; renderBody();
    if (scrollTo && indices.length) { const k = st.order.indexOf(indices[0]); if (k >= 0) scroll.scrollTop = Math.max(0, k * st.rowH - scroll.clientHeight / 2); }
  }
  const selection = () => [...st.sel];
  new ResizeObserver(renderBody).observe(scroll);
  /** The columns shown, by id (the self-check's). */
  const shownColumns = () => visible().map((c) => c.id);
  return { setRows, setSelection, selection, copy, shownColumns, rowCount: () => rows.length, refresh: () => { computeOrder(); renderHead(); renderBody(); } };
}
