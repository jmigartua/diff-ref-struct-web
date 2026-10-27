/* diff-ref-struct · the dock (APP_CONVENTIONS §1, C1.1–C1.9): sidebar | main | inspector, the bottom panel under
   the main area, on dockview-core (DESIGN.md §9, answer 6). Dividers (grey, 4 px, the panels' corners rounded) drag and double-click to collapse,
   panels have minimum and maximum sizes, every panel detaches into a floating window inside the app window (the ⧉
   button on its header; C1.8) and docks back, the layout is saved (in the project's .drs/project.json when a project
   is open, else in the settings) and Window ▸ Reset Layout returns to the default. The main area's tabs are the open
   files. */
import { createDockview, DefaultTab } from "./vendor/dockview/dockview-core.esm.js";
import { $, el, state, bus } from "./core.js";

export let dock = null;
let dockRoot = null;          // the dock's container element (TU5's divider handles)

class ModuleTab extends DefaultTab {
  init(params) {
    super.init(params);
    const m = render && render.tabModule ? render.tabModule((params.params || {}).tab) : null;
    if (!m) return;
    const i = document.createElement("span"); i.className = "tab-mod"; i.innerHTML = m.icon; i.title = `In ${m.name}`;
    this.element.insertBefore(i, this.element.firstChild);
    this.element.dataset.module = m.id;
    const name = () => this.element.setAttribute("aria-label", `${params.api.title || ""} in ${m.name}`);
    name(); this.addDisposables(params.api.onDidTitleChange(name));
  }
}
let render = null;          // { sidebar(el), inspector(el), bottom(el), log(el), file(tab, el) -> dispose, start(el) }
let onSave = null;          // (json) => void
let suppressSave = false;
const hidden = new Set();   // zones hidden by the user: "sidebar" | "inspector" | "bottom"
const SIZES = { sidebar: { w: 250, min: 180, max: 520 }, inspector: { w: 380, min: 240, max: 640 }, bottom: { h: 230, min: 120, max: 600 }, bottomRight: { w: 420, min: 260, max: 900 } };
let bottomSide = "bottom";   // §1: the bottom panel docked under the main area, or on its right ("right"); kept with the layout

export function initDock(container, renderers, save) {
  render = renderers; onSave = save; dockRoot = container;
  bus.on("layout", () => renderStubs());
  dock = createDockview(container, {
    theme: { name: "drs", className: "dv-theme-drs", gap: 4 },
    disableFloatingGroups: false,
    floatingGroupBounds: "boundedWithinViewport",
    createRightHeaderActionComponent: (group) => headerActions(group),
    // a file's tab carries its module's icon, with the module's name in the tooltip and the accessible name, so that
    // two tabs of one file (Refine and the Editor) are told apart (L122)
    createTabComponent: (opts) => (opts.name === "file" ? new ModuleTab() : undefined),
    createComponent: (opts) => {
      const body = el("div", "pane pane-" + opts.name);
      let dispose = null;
      return {
        element: body,
        init(params) {
          const p = params.params || {};
          if (opts.name === "file") { dispose = render.file(p.tab, body, params.api) || null; return; }
          if (opts.name === "start") { render.start(body); return; }
          dispose = render[opts.name] ? render[opts.name](body, params.api) || null : null;
        },
        update(params) { if (opts.name === "file" && render.fileUpdate) render.fileUpdate(params.params.tab, body); },
        dispose() { if (dispose) dispose(); },
      };
    },
  });
  dock.onDidLayoutChange(() => { if (!settling) layoutDefault = false; clearTimeout(saveLayout._t); saveLayout._t = setTimeout(saveLayout, 300); clearTimeout(tagDividers._t); tagDividers._t = setTimeout(tagDividers, 120); });
  dock.onDidActivePanelChange((p) => { if (p && p.params && p.params.tab) bus.emit("tab", p.params.tab); });
  dock.onDidRemovePanel((p) => {
    if (p.params && p.params.tab && !resetting) bus.emit("tab-closed", p.params.tab);
    for (const z of Object.keys(ZONE_PANEL)) if (p.id === ZONE_PANEL[z] && !suppressSave) { hidden.add(z); setTimeout(() => bus.emit("layout"), 0); }
    if (!resetting && !dock.panels.some((q) => q.params && (q.params.tab || q.id === "start"))) setTimeout(ensureStart, 0);
  });
  installSashDoubleClick(container); installSashCorners(container);
  // the dock takes the container's size, now and whenever the window changes
  const fit = () => { try { dock.layout(container.clientWidth, container.clientHeight, true); } catch (e) { /* */ } };
  new ResizeObserver(fit).observe(container);
  window.addEventListener("resize", fit);
  relayout._fit = fit;
}

/** The ⧉ / ⊡ button at the right of every group's header: detach into a floating panel, dock back (C1.8). */
function headerActions(group) {
  const element = el("div", "dv-actions-drs");
  // ? : the manual of the module shown, at the section for its view (APP_CONVENTIONS §5, L174)
  const q = el("button", "dv-act dv-help", "?");
  q.title = "The manual for what is on screen"; q.setAttribute("aria-label", "Manual for what is on screen");
  q.onclick = (e) => { e.stopPropagation(); if (render && render.help) render.help(zoneOf() || null); };
  element.appendChild(q);
  const b = el("button", "dv-act", "⧉");
  element.appendChild(b);
  let sub = null;
  // a panel (the sidebar, the inspector, the bottom panel) goes to a macOS window of its own (0.37.0); a group of file
  // tabs floats inside the window, as before (its tabs move to windows of their own by Move Tab to New Window)
  const zoneOf = () => Object.keys(ZONE_PANEL).find((k) => group.panels.some((p) => p.id === ZONE_PANEL[k]));
  const refresh = () => { const floating = group.api.location.type === "floating"; b.textContent = floating ? "⊡" : "⧉"; b.title = floating ? "Dock this panel back into the window" : zoneOf() && render && render.detachWindow ? "Detach this panel into a window of its own; closing that window brings it back" : "Detach this panel: a window of its own inside the app window, moved by its header"; b.setAttribute("aria-label", b.title); };
  b.onclick = async (e) => {
    e.stopPropagation();
    if (group.api.location.type === "floating") { dockBack(group); return; }
    const z = zoneOf();
    if (z && render && render.detachWindow && (await render.detachWindow(z))) return;
    detach(group);
  };
  return {
    element,
    init() { refresh(); try { sub = group.api.onDidLocationChange(refresh); } catch (e) { /* */ } },
    dispose() { if (sub) sub.dispose(); },
  };
}
export function detach(group) {
  if (!dock) return;
  const z = Object.keys(ZONE_PANEL).find((k) => group.panels.some((p) => p.id === ZONE_PANEL[k]));
  // the layout's container (dockview's API object has no element: the ⧉ of a group of tabs threw here)
  const box = (dockRoot || document.body).getBoundingClientRect();
  const w = z === "bottom" ? 640 : z ? SIZES[z].w + 40 : 720, h = z === "bottom" ? 300 : 420;
  try { group.locked = false; } catch (e) { /* */ }
  dock.addFloatingGroup(group, { position: { left: Math.max(20, (box.width - w) / 2), top: 60 }, width: w, height: h });
  saveLayout(); bus.emit("layout");
}
export function dockBack(group) {
  if (!dock) return;
  const main = mainGroup();
  const z = Object.keys(ZONE_PANEL).find((k) => group.panels.some((p) => p.id === ZONE_PANEL[k]));
  const target = main && main !== group ? main : dock.groups.find((g) => g !== group && g.api.location.type === "grid");
  if (!target) return;
  const position = z === "sidebar" ? "left" : z === "inspector" ? "right" : z === "bottom" ? (bottomSide === "right" ? "right" : "bottom") : "center";
  // the sidebar and the inspector back at the layout's own edges, not beside the main group (with a split the inspector
  // came back between its two groups: 0.42.1's relaunch)
  if (z === "sidebar" || z === "inspector") {
    const g = dock.addGroup({ direction: z === "sidebar" ? "left" : "right" });
    for (const p of group.panels.slice()) { try { p.api.moveTo({ group: g }); } catch (e) { console.warn("dock back", e); } }
  } else for (const p of group.panels.slice()) { try { p.api.moveTo({ group: target, position }); } catch (e) { console.warn("dock back", e); } }
  if (z) { const p = zonePanel(z); if (p) decorate(p, z); }
  saveLayout(); bus.emit("layout");
}

export function saveLayout() {
  if (!dock || suppressSave || !onSave) return;
  // the toolbar's visibility with the panels' (APP_CONVENTIONS §1: sizes and visibility remembered per project)
  const tb = document.getElementById("toolbar");
  try { onSave({ dock: dock.toJSON(), hidden: [...hidden], toolbarHidden: !!(tb && tb.hidden), bottomSide }); } catch (e) { /* nothing */ }
}

const ZONE_PANEL = { sidebar: "sidebar", inspector: "inspector", bottom: "bottom" };

function mainGroup() {
  const p = dock.panels.find((q) => q.params && (q.params.tab || q.id === "start"));
  return p ? p.group : null;
}
function zonePanel(zone) { return dock.getPanel(ZONE_PANEL[zone]); }

/** The default layout, as a grid with explicit sizes: sidebar | (main over bottom) | inspector. */
export function buildDefault() {
  suppressSave = true;
  dock.clear(); hidden.clear(); bottomSide = "bottom";
  const box = dockRoot ? dockRoot.parentElement : null;
  const W = Math.max(700, (box && box.clientWidth) || dock.width || 1200), H = Math.max(400, (box && box.clientHeight) || dock.height || 800);
  const sb = SIZES.sidebar.w, ins = SIZES.inspector.w, bt = SIZES.bottom.h;
  const centreW = W - sb - ins;
  let gid = 0;
  const leaf = (views, size) => ({ type: "leaf", data: { views, activeView: views[0], id: String(++gid) }, size });
  // the bottom panel stays folded until something fills it (a refinement's run, the Data table on demand): the
  // plot takes at least half the window (W7)
  void bt;
  const grid = {
    root: { type: "branch", data: [
      leaf(["sidebar"], sb),
      leaf(["start"], centreW),
      leaf(["inspector"], ins),
    ], size: H },
    width: W, height: H, orientation: "HORIZONTAL",
  };
  const panels = {
    sidebar: { id: "sidebar", contentComponent: "sidebar", title: "Project" },
    start: { id: "start", contentComponent: "start", title: "Start" },
    inspector: { id: "inspector", contentComponent: "inspector", title: "Inspector" },
  };
  hidden.add("bottom");
  try { dock.fromJSON({ grid, panels, activeGroup: "2" }); }
  catch (e) {
    console.error("default layout refused, building incrementally", e);
    dock.clear();
    dock.addPanel({ id: "start", component: "start", title: "Start" });
    const main = mainGroup();
    addZone("sidebar", main); addZone("inspector", main); addZone("bottom", main);
  }
  for (const z of ["sidebar", "inspector", "bottom"]) { const p = zonePanel(z); if (p) decorate(p, z); }
  const s = dock.getPanel("start"); if (s) s.api.setActive();
  suppressSave = false;
  relayout();
  saveLayout();
}
function addZone(zone, main) {
  const sz = SIZES[zone];
  const common = { id: ZONE_PANEL[zone], component: zone, title: zone === "sidebar" ? "Project" : zone === "inspector" ? "Inspector" : "Bottom panel" };
  let p;
  if (zone === "bottom" && bottomSide === "right") { const to = rightOfMain(main); p = dock.addPanel({ ...common, position: to.referenceGroup ? to : { direction: "right" } }); }
  else if (zone === "bottom") p = dock.addPanel({ ...common, position: { referenceGroup: main, direction: "below" }, initialHeight: sz.h });
  else p = dock.addPanel({ ...common, position: { referenceGroup: main, direction: zone === "sidebar" ? "left" : "right" }, initialWidth: sz.w });
  decorate(p, zone);
  const size = () => { try { if (zone === "bottom" && bottomSide === "right") return; else if (zone === "bottom") p.api.setSize({ height: sz.h }); else p.api.setSize({ width: sz.w }); } catch (e) { /* */ } };
  size(); setTimeout(size, 0);
  return p;
}
function decorate(p, zone) {
  const sz = SIZES[zone];
  try { p.group.locked = p.group.api.location.type === "floating" ? false : "no-drop-target"; p.group.header.hidden = false; } catch (e) { /* older dockview */ }
  const free = Number.MAX_SAFE_INTEGER;
  try {
    if (zone === "bottom" && bottomSide === "right") p.api.setConstraints({ minimumWidth: SIZES.bottomRight.min, maximumWidth: SIZES.bottomRight.max, minimumHeight: 0, maximumHeight: free });
    else if (zone === "bottom") p.api.setConstraints({ minimumHeight: sz.min, maximumHeight: sz.max, minimumWidth: 0, maximumWidth: free });
    else p.api.setConstraints({ minimumWidth: sz.min, maximumWidth: sz.max });
  } catch (e) { /* */ }
  if (zone === "bottom") {
    try { if (!dock.getPanel("log")) dock.addPanel({ id: "log", component: "log", title: "Log", position: { referencePanel: p.id, direction: "within" }, inactive: true }); } catch (e) { /* */ }
    try { p.api.setActive(); } catch (e) { /* */ }
  }
}

/** Restores a saved layout, or builds the default. */
export function restore(saved) {
  if (saved && saved.dock) {
    suppressSave = true;
    try {
      dock.fromJSON(saved.dock);
      hidden.clear(); for (const z of saved.hidden || []) hidden.add(z);
      bottomSide = saved.bottomSide === "right" ? "right" : "bottom";
      { const tb = document.getElementById("toolbar"); if (tb) tb.hidden = !!saved.toolbarHidden; }
      for (const z of ["sidebar", "inspector", "bottom"]) { const p = zonePanel(z); if (p) decorate(p, z); }
      // file tabs do not survive a restart yet: their panels are dropped, the start page shown
      for (const p of dock.panels.slice()) if (p.params && p.params.tab) p.api.close();
      ensureStart();
      suppressSave = false;
      relayout(); renderStubs();
      return;
    } catch (e) { console.warn("saved layout unusable, rebuilding", e); suppressSave = false; }
  }
  buildDefault();
  renderStubs();
}

export function ensureStart() {
  if (!dock) return;
  if (dock.panels.some((q) => q.params && q.params.tab)) { const s = dock.getPanel("start"); if (s) s.api.close(); return; }
  if (!dock.getPanel("start")) {
    const ref = zonePanel("bottom") || zonePanel("sidebar") || zonePanel("inspector");
    const pos = ref ? { referenceGroup: ref.group, direction: ref.id === "bottom" ? "above" : ref.id === "sidebar" ? "right" : "left" } : undefined;
    dock.addPanel({ id: "start", component: "start", title: "Start", position: pos });
  }
}

/** Show or hide a zone (View ▸ Show/Hide …). */
/** The status bar's click (§9: "the full log a click away"): the bottom panel shown with its Log tab in front. */
export function showLog() {
  if (!dock) return;
  showZone("bottom");
  const p = dock.getPanel("log"); if (p) { try { p.api.setActive(); } catch (e) { /* */ } }
}
export function toggleZone(zone) {
  const p = zonePanel(zone);
  if (p) { hidden.add(zone); p.api.close(); }
  else { hidden.delete(zone); showZone(zone); }
  saveLayout(); bus.emit("layout");
}
/** The zone shown whatever happened to it (closed, detached and closed, its group gone): added back beside the main
    group, or beside any grid group, or alone (app 0.9.2: a hidden inspector could not be recovered). */
export function showZone(zone) {
  if (!dock) return;
  const p = zonePanel(zone);
  if (p) { try { p.api.setActive(); } catch (e) { /* */ } return; }
  hidden.delete(zone);
  const main = mainGroup() || dock.groups.find((g) => g.api.location.type === "grid") || null;
  if (main) addZone(zone, main);
  else { dock.addPanel({ id: ZONE_PANEL[zone], component: zone, title: zone === "sidebar" ? "Project" : zone === "inspector" ? "Inspector" : "Bottom panel" }); const q = zonePanel(zone); if (q) decorate(q, zone); }
  relayout();
}
export const zoneShown = (zone) => !!zonePanel(zone);
/** Where the bottom panel goes on the right: left of the inspector when the inspector is the layout's last column, else
    at the layout's right edge; the side panels' widths are put back after the move. */
function rightOfMain(main, own = null) {
  // left of the inspector when it is the layout's last column; else at the layout's right edge: right of every group of
  // tabs and every row of them, whatever a relaunch did with the inspector (0.41.1's install: at the left; 0.42.1's
  // relaunch: between the groups of a split)
  const ins = zonePanel("inspector");
  // the layout's own right edge: its container's (dockview's API object has no element of its own)
  const edge = dockRoot ? dockRoot.getBoundingClientRect().right : Infinity;
  const root = dockRoot ? dockRoot.getBoundingClientRect() : null;
  const ir = ins && ins.group ? ins.group.element.getBoundingClientRect() : null;
  // the inspector's column the layout's last and its whole height (0.42.1's relaunch: a column shared below it)
  if (ins && ins.group !== own && ins.group.api.location.type === "grid" && ir.right >= edge - 8 && (!root || ir.height >= root.height - 40)) return { referenceGroup: ins.group, direction: "left" };
  return { direction: "right" };
}
/** §1: where the bottom panel is docked, "bottom" (under the main area) or "right" (beside it). */
export const bottomDock = () => bottomSide;
/** View ▸ Bottom Panel on the Right / at the Bottom: the panel's group (its tabs, the Log among them) moved beside or
    under the main area, its size and limits turned with it; kept with the layout. A panel in a window of its own or
    floating keeps its place: the side applies when it docks back. */
export function setBottomDock(side) {
  bottomSide = side === "right" ? "right" : "bottom";
  // a bottom panel left floating (a relaunch restores a floating group as it was kept) is docked back first, then moved
  { const f = zonePanel("bottom"); if (f && f.group && f.group.api.location.type === "floating") dockBack(f.group); }
  const p = zonePanel("bottom");
  if (p && dock && p.group.api.location.type === "grid") {
    const main = mainGroup() || dock.groups.find((g) => g !== p.group && g.api.location.type === "grid");
    if (main && main !== p.group) {
      // only the bottom zone's own panels (the Log) go with it: a file tab that a relaunch left in its group stays
      // (0.40.0's relaunch run: the tabs went along, and the group they made was taken for the main area's)
      const old = p.group, active = old.activePanel, others = old.panels.filter((q) => q !== p && !(q.params && q.params.tab) && q.id !== "start");
      // the side panels' widths, put back after the move (a column taken away is shared out among the others)
      const kept = {}; for (const z of ["sidebar", "inspector"]) { const q = zonePanel(z); if (q && q.group.api.location.type === "grid") kept[z] = q.api.width; }
      // the panel moved beside its reference as Dock Back moves a panel, the Log and the others into its new group
      const to = bottomSide === "right" ? rightOfMain(main, p.group) : { referenceGroup: main, direction: "below" };
      try {
        if (to.referenceGroup) p.api.moveTo({ group: to.referenceGroup, position: to.direction === "below" ? "bottom" : to.direction });
        else { const g = dock.addGroup({ direction: to.direction }); p.api.moveTo({ group: g }); }   // the layout's edge
      } catch (e) { console.warn("bottom panel", e); }
      for (const q of others) { try { q.api.moveTo({ group: p.group, position: "center" }); } catch (e) { console.warn("bottom panel", e); } }
      decorate(p, "bottom");
      // under the main area, its usual height; on the right, half of the tabs' group it sits beside (a width given
      // here would be taken from the inspector)
      try { if (bottomSide !== "right") p.api.setSize({ height: SIZES.bottom.h }); } catch (e) { /* */ }
      for (const z of Object.keys(kept)) { const q = zonePanel(z); try { if (q && Math.abs(q.api.width - kept[z]) > 1) q.api.setSize({ width: kept[z] }); } catch (e) { /* */ } }
      if (active) { try { active.api.setActive(); } catch (e) { /* */ } }
    }
  }
  relayout(); saveLayout(); bus.emit("layout");
}
/** The bottom panel's tab carries the active module's name. */
export function setBottomTitle(title) { const p = zonePanel("bottom"); if (p) p.api.setTitle(title); }

/** A file tab in the main area. */
export function openTab(tab, opts = {}) {
  const ex = dock.getPanel("tab:" + tab.id);
  if (ex) { ex.api.setActive(); return ex; }
  const start = dock.getPanel("start");
  // beside: in a group to the right of the tab named (Open as Text beside its plot, W6)
  const ref = opts.beside ? dock.getPanel("tab:" + opts.beside) : null;
  const main = ref ? ref.group : mainGroup();
  const p = dock.addPanel({ id: "tab:" + tab.id, component: "file", tabComponent: "file", title: tab.title, params: { tab: tab.id }, position: main ? { referenceGroup: main, direction: opts.beside ? "right" : "within" } : undefined });
  if (start) start.api.close();
  return p;
}
/** A file tab dragged out of the window (§1) calls `fn(tabId)`: the drag of a tab that ends outside the window's frame
    (macOS gives the drop's place in screen points), which dockview leaves where it was. */
export function onTabDraggedOut(fn) {
  let dragged = null;
  document.addEventListener("dragstart", (e) => {
    const t = e.target && e.target.closest && e.target.closest(".dv-tab");
    const p = t && dock ? dock.panels.find((q) => q.params && q.params.tab && q.view && q.view.tab && q.view.tab.element && (q.view.tab.element === t || t.contains(q.view.tab.element) || q.view.tab.element.contains(t))) : null;
    dragged = p ? p.params.tab : null;
  }, true);
  document.addEventListener("dragend", (e) => {
    const id = dragged; dragged = null;
    if (!id || (e.screenX === 0 && e.screenY === 0)) return;
    const outside = e.screenX < window.screenX || e.screenX > window.screenX + window.outerWidth || e.screenY < window.screenY || e.screenY > window.screenY + window.outerHeight;
    if (outside) fn(id);
  }, true);
}

/** View ▸ Split Right / Split Down (APP_CONVENTIONS §1): the active file tab moved into a new group beside ("right")
    or below ("below") its own, so that two files show side by side or one above the other. Tabs are then dragged
    between the groups. Nothing to split when the tab is alone in its group. */
export function canSplit() { const p = dock && dock.activePanel; return !!(p && p.params && p.params.tab && p.group && p.group.panels.length > 1); }
export function splitActive(direction) {
  if (!canSplit()) return false;
  const p = dock.activePanel;
  const g = dock.addGroup({ referenceGroup: p.group, direction });
  p.api.moveTo({ group: g });
  p.api.setActive();
  return true;
}
/** The groups that hold file tabs (P46's rest, 0.35.0), in the order they appear: each with its tabs in order and,
    from the second on, the earlier group it was split from and how, "right" or "below", read from where the groups
    lie on screen. */
export function tabGroups() {
  const groups = [];
  // each tab by the group its header sits in on screen (the group view around it), else by the panel's group: the
  // installed app's check counted one group after Split Down made a second (0.37.0, the merging session's run)
  const holder = (p) => { const t = p.view && p.view.tab && p.view.tab.element; return (t && t.closest && t.closest(".dv-groupview")) || (p.group && p.group.element) || p.group; };
  for (const p of dock.panels) {
    if (!(p.params && p.params.tab)) continue;
    const key = holder(p);
    let g = groups.find((x) => x.key === key);
    if (!g) { g = { key, group: p.group, tabs: [] }; groups.push(g); }
    g.tabs.push(p.params.tab);
  }
  for (const g of groups) { const el = g.key && g.key.getBoundingClientRect ? g.key : g.group && g.group.element; g.rect = el ? el.getBoundingClientRect() : { left: 0, top: 0, right: 0, bottom: 0, width: 0, height: 0 }; }
  groups.sort((a, b) => a.rect.top - b.rect.top || a.rect.left - b.rect.left);
  // each group's size on screen (§1, app 0.42.0: the split sizes kept with the project's tabs)
  const size = (g) => ({ w: Math.round(g.rect.width || 0), h: Math.round(g.rect.height || 0) });
  return groups.map((g, k) => Object.assign(placeOf(g, k), size(g)));
  function placeOf(g, k) {
    if (k === 0) return { tabs: g.tabs };
    const r = g.rect;
    const over = (a0, a1, b0, b1) => Math.min(a1, b1) - Math.max(a0, b0) > 4;
    for (let j = k - 1; j >= 0; j--) {
      const q = groups[j].rect;
      if (r.left >= q.right - 4 && over(r.top, r.bottom, q.top, q.bottom)) return { tabs: g.tabs, from: j, direction: "right" };
      if (r.top >= q.bottom - 4 && over(r.left, r.right, q.left, q.right)) return { tabs: g.tabs, from: j, direction: "below" };
    }
    return { tabs: g.tabs, from: 0, direction: r.left > groups[0].rect.left + 4 ? "right" : "below" };
  }
}
/** §1 (app 0.42.0): the group of a file tab given the size it was kept with, in pixels (a split's sizes come back with
    the project's tabs, not equal). */
export function sizeTabGroup(tabId, w, h) {
  const p = dock && dock.getPanel("tab:" + tabId); if (!p || !p.group) return false;
  try { if (w > 0) p.group.api.setSize({ width: w }); if (h > 0) p.group.api.setSize({ height: h }); return true; } catch (e) { return false; }
}
/** `fn()` when the layout's sizes change (a divider dragged), for the sizes kept with the project's tabs. */
export function onLayoutResized(fn) { if (dock) dock.onDidLayoutChange(() => fn()); }
/** `fn()` after a group of tabs is added or removed or a tab moves between groups (the groups kept with the project,
    P46's rest: the layout's own save is debounced and may come late; 0.35.0's self-check). */
export function onGroupsChanged(fn) {
  if (!dock) return;
  for (const ev of ["onDidAddGroup", "onDidRemoveGroup", "onDidMovePanel"]) if (dock[ev]) dock[ev](() => setTimeout(fn, 0));
}
/** The tabs `tabIds` moved into a new group split off `direction` ("right", "below") of the group that holds
    `besideTabId` (the groups restored with the window, P46). */
export function splitTabs(tabIds, besideTabId, direction) {
  const ref = dock.getPanel("tab:" + besideTabId); if (!ref) return false;
  const panels = tabIds.map((id) => dock.getPanel("tab:" + id)).filter(Boolean);
  if (!panels.length) return false;
  const g = dock.addGroup({ referenceGroup: ref.group, direction });
  for (const p of panels) p.api.moveTo({ group: g });
  panels[0].api.setActive();
  return true;
}
/** The file tab a dock tab's header element shows (a context menu on it, §6), or null. */
export function tabIdAt(el) {
  const t = el && el.closest ? el.closest(".dv-tab") : null; if (!t || !dock) return null;
  const p = dock.panels.find((q) => q.params && q.params.tab != null && q.view && q.view.tab && q.view.tab.element && (q.view.tab.element === t || t.contains(q.view.tab.element) || q.view.tab.element.contains(t)));
  return p ? p.params.tab : null;
}
export function closeTab(tabId) { const p = dock.getPanel("tab:" + tabId); if (p) p.api.close(); }
export function setTabTitle(tabId, title) { const p = dock.getPanel("tab:" + tabId); if (p) p.api.setTitle(title); }
export function activateTab(tabId) { const p = dock.getPanel("tab:" + tabId); if (p) p.api.setActive(); }
export function tabOrder() { return dock.panels.filter((p) => p.params && p.params.tab).map((p) => p.params.tab); }
export function nextTab(delta = 1) {
  const ids = tabOrder(); if (ids.length < 2) return;
  const i = ids.indexOf(state.activeTab);
  activateTab(ids[(i + delta + ids.length) % ids.length]);
}

/** Where a vertical and a horizontal divider meet, one drag moves both (as VS Code and TRAE do; APP_CONVENTIONS §1,
    L149): over the crossing the pointer shows the move cursor, and a press there starts dockview's own drag of each
    divider (each follows the pointer from the document, so one movement drives both). The crossing is within 8 px of
    both dividers' centre lines, a T-junction included. */
function installSashCorners(container) {
  const NEAR = 8;
  const sashes = () => [...container.querySelectorAll(".dv-sash")].filter((e) => e.offsetParent !== null).map((e) => ({ e, r: e.getBoundingClientRect() })).filter((q) => q.r.width > 0 && q.r.height > 0);
  const crossing = (x, y) => {
    const all = sashes();
    const vert = all.filter((q) => q.r.height > q.r.width && Math.abs(x - (q.r.left + q.r.right) / 2) <= NEAR && y >= q.r.top - NEAR && y <= q.r.bottom + NEAR);
    const hori = all.filter((q) => q.r.width > q.r.height && Math.abs(y - (q.r.top + q.r.bottom) / 2) <= NEAR && x >= q.r.left - NEAR && x <= q.r.right + NEAR);
    return vert.length && hori.length ? [vert[0].e, hori[0].e] : null;
  };
  container.addEventListener("pointermove", (e) => {
    if (e.buttons) return;
    container.classList.toggle("sash-corner", !!crossing(e.clientX, e.clientY));
  });
  container.addEventListener("pointerleave", () => container.classList.remove("sash-corner"));
  container.addEventListener("pointerdown", (e) => {
    if (!e.isTrusted || e.button !== 0) return;   // the presses sent below come back here; they pass
    const both = crossing(e.clientX, e.clientY); if (!both) return;
    e.preventDefault(); e.stopPropagation();
    for (const s of both) s.dispatchEvent(new PointerEvent("pointerdown", { bubbles: true, cancelable: true, clientX: e.clientX, clientY: e.clientY, pointerId: e.pointerId, pointerType: e.pointerType, button: 0, buttons: 1, isPrimary: true }));
  }, true);
}

/** A double-click on a divider collapses the panel beside it (C1.3); the divider nearest the double-click decides. */
/** A hidden panel leaves a stub edge where it was (APP_CONVENTIONS §1, L224): a thin strip at the window's left (the
    sidebar), right (the inspector) or bottom (the bottom panel), which a click or a double-click brings the panel back
    from, as a double-click on its divider took it away. The View menu and the toolbar's toggles do the same. */
const STUB_TITLES = { sidebar: "Show the Sidebar (⌃⌘S)", inspector: "Show the Inspector (⌥⌘I)", bottom: "Show the Bottom Panel (⇧⌘Y)" };
export function renderStubs() {
  if (!dockRoot || !dock) return;
  if (getComputedStyle(dockRoot).position === "static") dockRoot.style.position = "relative";
  for (const z of ["sidebar", "inspector", "bottom"]) {
    let st = dockRoot.querySelector(`:scope > .zone-stub[data-zone="${z}"]`);
    const show = !zonePanel(z);
    if (show && !st) {
      st = document.createElement("div"); st.className = "zone-stub"; st.dataset.zone = z; st.dataset.tour = "stub:" + z;
      st.title = STUB_TITLES[z]; st.setAttribute("role", "button"); st.setAttribute("aria-label", STUB_TITLES[z]);
      st.onclick = () => { if (!zonePanel(z)) toggleZone(z); };
      dockRoot.appendChild(st);
    }
    if (!show && st) st.remove();
  }
}
function installSashDoubleClick(container) {
  container.addEventListener("dblclick", (e) => {
    const sash = e.target.closest(".dv-sash");
    if (!sash) return;
    const r = sash.getBoundingClientRect(), c = container.getBoundingClientRect();
    const vertical = r.height > r.width;
    if (vertical) {
      const side = zonePanel("sidebar"), insp = zonePanel("inspector");
      if (side && Math.abs(r.left - (side.group.element.getBoundingClientRect().right)) < 12) return toggleZone("sidebar");
      if (insp && Math.abs(r.right - (insp.group.element.getBoundingClientRect().left)) < 12) return toggleZone("inspector");
      if (r.left - c.left < c.width / 2 && side) return toggleZone("sidebar");
      if (insp) return toggleZone("inspector");
    } else {
      if (zonePanel("bottom")) return toggleZone("bottom");
    }
  });
}

/** Window ▸ Reset Layout: the default arrangement, the open tabs kept in the main area and the one in use still active
    (TU10: the reset cleared the dock, and every tab's removal was taken for a closing). */
let resetting = false;
export function resetLayout() {
  const tabs = dock ? dock.panels.filter((p) => p.params && p.params.tab != null).map((p) => ({ id: p.params.tab, title: p.title || p.api.title || "" })) : [];
  const active = dock && dock.activePanel && dock.activePanel.params ? dock.activePanel.params.tab : null;
  resetting = true;
  try { buildDefault(); for (const t of tabs) openTab(t); }
  finally { resetting = false; }
  if (active != null) activateTab(active);
  bus.emit("layout"); markDefault();
}
/** Whether the layout is the default one (the tour's check after Window ▸ Reset Layout): true from a reset until the
    user next changes the layout. */
let layoutDefault = false, settling = false;
function markDefault() { layoutDefault = true; settling = true; setTimeout(() => { settling = false; }, 800); }
export function layoutIsDefault() { return layoutDefault; }
/** TU5: each zone's divider tagged data-tour="divider:<zone>" (the sash beside the sidebar's, the inspector's or the
    bottom panel's group), for the tour player; done again after every change of the layout. */
export function tagDividers() {
  if (!dock) return;
  if (!dockRoot) return;
  const sashes = [...dockRoot.querySelectorAll(".dv-sash")].filter((e) => e.offsetParent !== null).map((e) => ({ e, r: e.getBoundingClientRect() })).filter((q) => q.r.width > 0 && q.r.height > 0);
  for (const q of sashes) if (q.e.dataset.tour && q.e.dataset.tour.startsWith("divider:")) delete q.e.dataset.tour;
  for (const [zone, pid] of Object.entries(ZONE_PANEL)) {
    const p = dock.getPanel(pid); if (!p || !p.group || !p.group.element) continue;
    if (p.group.api.location.type !== "grid") continue;
    const g = p.group.element.getBoundingClientRect();
    let best = null, bd = 14;
    for (const { e, r } of sashes) {
      const cx = (r.left + r.right) / 2, cy = (r.top + r.bottom) / 2; let dist = Infinity;
      if (zone === "sidebar" && r.height > r.width && cy > g.top && cy < g.bottom) dist = Math.abs(cx - g.right);
      if (zone === "inspector" && r.height > r.width && cy > g.top && cy < g.bottom) dist = Math.abs(cx - g.left);
      if (zone === "bottom" && r.width > r.height && cx > g.left && cx < g.right) dist = Math.abs(cy - g.top);
      if (dist < bd) { bd = dist; best = e; }
    }
    if (best) best.dataset.tour = "divider:" + zone;
  }
}
export function relayout() { if (relayout._fit) relayout._fit(); else if (dock) dock.layout(dock.width, dock.height, true); }
