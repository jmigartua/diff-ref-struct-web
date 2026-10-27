/* diff-ref-struct · the command registry (the difes shell's way). Every command of the app and of its modules is
   declared once and generated everywhere: the native menu bar, the toolbar, the ⌘K palette and the hover cards.
   A command cannot be in one surface and missing in another, because all surfaces are this list (C3.8). */
import { esc } from "./core.js";

export const COMMANDS = [];
const byId = new Map();

/** cmd: { id, title, group ("App" | "File" | "Edit" | "View" | "Window" | "Help" | a module id), shortcut?
    ("⇧⌘S"), what, affects?, run(arg), when?() -> bool, disabled? (why, when `when` is false by design), menu?:
    false (not in the menu bar), palette?: false, toolbar?: label, sep?: true (a separator before it), dynamic?(d)
    -> title, predefined?: a Tauri predefined item, submenu?: [ids] (a submenu of other commands), check?() -> bool,
    field?: { get() -> text, set(text), width?, placeholder? } (an input in the toolbar after its label, not a menu item) } */
export function register(cmds) {
  for (const c of cmds) { if (byId.has(c.id)) { Object.assign(byId.get(c.id), c); continue; } byId.set(c.id, c); COMMANDS.push(c); }
}
export const cmd = (id) => byId.get(id) || null;
export const enabled = (c) => !c.when || !!c.when();
export async function run(id, arg) {
  const c = byId.get(id);
  if (!c) { console.warn("no command", id); return; }
  if (!enabled(c)) return;
  return c.run(arg);
}

// ------------------------------------------------------------------ the generated surfaces
export function paletteItems(q, extra = []) {
  const items = [];
  for (const c of COMMANDS) if (enabled(c) && c.palette !== false && (!c.predefined || c.palette === true)) items.push({ group: c.group, glyph: "›", title: c.title, sub: c.what || "", kbd: c.shortcut || "", run: () => c.run() });
  items.push(...extra);
  const words = q.toLowerCase().split(/\s+/).filter(Boolean);
  return items.filter((it) => words.every((w) => (it.title + " " + it.sub + " " + it.group).toLowerCase().includes(w))).slice(0, 60);
}

export function cardHTML(id) {
  const c = byId.get(id); if (!c) return "";
  return `<div class="vc-h"><b>${esc(c.title)}</b>${c.shortcut ? `<kbd>${esc(c.shortcut)}</kbd>` : ""}</div>
    <div class="vc-what">${esc(c.what || "")}</div>
    ${c.affects ? `<div class="vc-aff">${esc(c.affects)}</div>` : ""}
    ${!enabled(c) && c.disabled ? `<div class="vc-aff">${esc(c.disabled)}</div>` : ""}
    <div class="vc-f">${c.menu === false ? "" : "in the menu bar · "}in ⌘K${c.group ? " · " + esc(c.group) : ""}</div>`;
}

/** The commands of a menu group, in registry order. */
export function menuItems(group) {
  return COMMANDS.filter((c) => c.group === group && c.menu !== false);
}

/** "⇧⌘S" → Tauri's accelerator "CmdOrCtrl+Shift+S". */
export function accelerator(sc) {
  if (!sc) return undefined;
  const mods = [];
  if (sc.includes("⌃")) mods.push("Ctrl");
  if (sc.includes("⌥")) mods.push("Alt");
  if (sc.includes("⇧")) mods.push("Shift");
  if (sc.includes("⌘")) mods.push("CmdOrCtrl");
  const key = sc.replace(/[⌘⇧⌥⌃]/g, "");
  const names = { ",": "Comma", ".": "Period", "↩": "Enter", "⌫": "Backspace", "⇥": "Tab", "␣": "Space", "+": "Plus", "−": "Minus", "-": "Minus", "=": "Equal", "←": "ArrowLeft", "→": "ArrowRight", "↑": "ArrowUp", "↓": "ArrowDown", "⌦": "Delete", "`": "Backquote" };
  return [...mods, names[key] || key.toUpperCase()].join("+");
}

/** Shortcut text → a predicate on a keyboard event (browser mode; in the app the native menus own the keys). */
export function shortcutMatches(sc, e) {
  if (!sc) return false;
  const want = { meta: sc.includes("⌘"), shift: sc.includes("⇧"), alt: sc.includes("⌥"), ctrl: sc.includes("⌃") };
  const key = sc.replace(/[⌘⇧⌥⌃]/g, "");
  const k = e.key.length === 1 ? e.key.toUpperCase() : e.key;
  const map = { "⇥": "Tab", "⌫": "Backspace", "↩": "Enter", "−": "-" };
  const want_key = (map[key] || key).toUpperCase();
  const keyOk = want_key === k.toUpperCase() || (want_key === "+" && (e.key === "+" || e.key === "="));
  return keyOk && e.metaKey === want.meta && e.ctrlKey === want.ctrl && e.shiftKey === want.shift && e.altKey === want.alt;
}
