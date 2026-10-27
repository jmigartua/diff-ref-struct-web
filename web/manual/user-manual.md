# diff-ref-struct Web — user manual

For **diff-ref-struct Web 1.1.0 with core 0.82.0**, at <https://jmigartua.github.io/diff-ref-struct-web/>. The page
shows the version it runs at the foot of its start page and in the line under the inspector, and its **Help ▾** menu
opens the manual of that version. The manual ends with three appendices made from the guide's records: **A, Conventions**
(the words, files, suffixes and names), **B, Release notes** (every version up to this one) and **C, Roadmap**.

The site's **Recipes** work through the analyses the page is made for (a temperature ramp, a hold at constant
temperature, two patterns compared, a run followed live), and its **Quick reference** fits the controls on one page.

---

## Contents

1. [What the page is](#1-what-the-page-is)
2. [What you need](#2-what-you-need)
3. [The page at a glance](#3-the-page-at-a-glance)
4. [Bringing files in](#4-bringing-files-in)
5. [The folder's tree](#5-the-folders-tree)
6. [What opens as what](#6-what-opens-as-what)
7. [Sets in this folder: the core's proposals](#7-sets-in-this-folder-the-cores-proposals)
8. [A plot set in its tab](#8-a-plot-set-in-its-tab)
9. [The members and the member](#9-the-members-and-the-member)
10. [The views](#10-the-views)
11. [The Map and the Difference Map in detail](#11-the-map-and-the-difference-map-in-detail)
12. [Normalising, converting x, windows and the Series view](#12-normalising-converting-x-windows-and-the-series-view)
13. [Saving and exporting](#13-saving-and-exporting)
14. [The Plots menu, item by item](#14-the-plots-menu-item-by-item)
15. [Live mode: a folder still being written](#15-live-mode-a-folder-still-being-written)
16. [The view in the address: links and reloads](#16-the-view-in-the-address-links-and-reloads)
17. [Offline and installed](#17-offline-and-installed)
18. [Mouse, keyboard, appearance](#18-mouse-keyboard-appearance)
19. [Messages you will meet, and what to do](#19-messages-you-will-meet-and-what-to-do)
20. [Limits of this version](#20-limits-of-this-version)
21. [The file formats the core reads](#21-the-file-formats-the-core-reads)
22. Appendices: A Conventions, B Release notes, C Roadmap (made from the guide's records, in the site and the PDF)

---

## 1. What the page is

diff-ref-struct Web is the **Plots module** of the diff-ref-struct Mac app, running in a browser page. It shows many
powder-diffraction patterns at once as a **plot set**. You can view them as:

- an overlay;
- a stack;
- a WinPLOTR-style projection;
- a colour **map** of intensity against 2θ and the set's variable (temperature or time);
- a difference map against a reference pattern;
- a 3D surface;
- a series of quantities taken from windows of 2θ.

Three things set it apart from a plotting web page:

- **The numbers are the core's.** The page contains the diff-ref-struct core, the same Rust engine the Mac app uses,
  compiled to WebAssembly. The core computes every number on screen: reading the files, the temperatures found in their
  names, the map's grid and colour range, the stack's steps, the projection, the differences, the windows' integrals. The
  page draws what the core returns and computes nothing itself. A map in the browser therefore equals the map in the
  app, number for number.
- **Your files never leave your machine.** The browser reads a folder or files where they are. Nothing is uploaded,
  there is no server behind the page, and you need no account. The page has no tracking. The line at the foot of the
  inspector says so as well.
- **It recognises series.** Open a folder of scans from a temperature ramp or a time-resolved run and the core proposes
  the plot sets it finds: one scan per temperature, each hold at constant temperature against time, and the whole run
  against time. Each proposal is one click away.

What it is **not**: it is not the refinement program. The Rietveld refinement, the Data module's operations,
Publish's figures and the project's files belong to the Mac app. The browser page is for looking at data,
comparing patterns, following a series and finding where things change.

## 2. What you need

| what | needed |
|---|---|
| **Browser** | Any current browser. **Chromium-based browsers** (Chrome, Edge, Arc, Brave) do the most: they open a folder in place, so that **live mode** works (§15). Firefox and Safari open a folder through the browser's folder dialog instead, which takes a copy of the list once; live mode is not offered there. |
| **Network** | Only for the first visit. After it the page keeps itself and opens offline (§17). |
| **Memory** | Ordinary sets (tens of files) need nothing special. A run of a hundred synchrotron scans of 27 000 points each (112 MB) works; allow a few hundred MB of memory for the tab. |
| **Files** | Text data files up to **64 MB** each (§21). Larger files are listed but not read. |

**Permissions.** When you open a folder, Chromium asks "Let site view files?". The page asks only to *read* and never
writes into the folder. Firefox and Safari may word their dialog as an "upload". It is not one: the files stay in the
browser.

<a id="help-page"></a>

## 3. The page at a glance

![The page with a plot set open: the Tb₂Pd₂Sn patterns of FullProf's WinPLOTR tutorial on the Map, colour on a log
scale clipped to 3.2–4.5, the cursors pinned at 2θ 31.6° (row strip below, column strip at the right). At the left,
the inspector with the plot set and its members.](../site/assets/figures/app-map.png)

From top to bottom: **the bar**; at the left **the inspector**, and the **separator** at its right edge; at the right
**the views**, one tab per plot set, each with its view tabs, toolbar, plot and note; at the foot **the status line**.

**The bar**, left to right:

- **Open Folder…**: chooses a folder of data files (§4).
- **Open Files…**: chooses one or more files: data files, a plot set file (`.drv`) or a WinPLOTR buffer (`.buf`).
- **Plots ▾**: the module's menu, the same list as the app's Plots menu (§14). Items that cannot be used now are
  greyed out, and hovering one says why.
- **The core's chip**: "core 0.81.0 in WebAssembly" once the core has loaded. If it reads "core: waiting…", the core
  did not load; reload the page.
- **Help ▾**: the guide; **Manual for this view** (the manual opened at the section of the view shown); **What's new in
  this version**; **Report a problem…** and **Ask for a feature…** (a form of the public repository's issues, filled in
  with the version, the browser and the view's address, never your data; or the same by e-mail); **Cite this view** (one
  line for a lab book or a paper's methods, copied to the clipboard).
- **Appearance**: "◐ auto" (the system's), "☀ light", "☾ dark". Click for the next one. The choice is kept.

**The inspector** (the column at the left), top to bottom:

- **Folder**: its name. `<` shrinks the inspector to a narrow rail and `>` brings it back.
- **Files**: the folder's tree (§5), with **Open as Plot Set** and, in Chromium, the **Live** controls (§15).
- **Sets in this folder**: the core's proposals (§7).
- **Plot sets, Windows, Regions, Members**: the panels of the plot set shown (§9).
- **The member**: the facts of the member chosen (fold it open).
- **The privacy line**, and under it the offline line (§17).

**The separator** between inspector and views: drag it to set the inspector's width (220–640 px), or double-click it
to shrink the inspector to a rail and back. With the keyboard, focus it and use ← → (16 px a step, 64 px with ⇧). The
width is kept for your next visit. On a narrow window (under 640 px) the inspector stacks above the views.

**The views** take the rest of the page: one tab per plot set, each with its view tabs and toolbar (§8). Before
anything is open, a start page shows the two Open buttons and the version.

**The status line** at the foot says what the last action did ("Plots: 39 members …") or why it could not. Messages in
red are refusals or problems. Read it whenever something seems not to happen.

## 4. Bringing files in

There are three ways, and every one of them reads the files in the browser only. A fourth needs no files of yours:
**Try the example** on the start page opens a synthetic heating run computed by the core (21 scans of a periclase sample
from 300 to 1300 K, with corundum appearing above about 1000 K), as if its folder had been dropped.

**Open Folder…** (bar or start page). Choose a folder. Its files appear in the tree, sub-folders included. Hidden files
and folders (names starting with `.`, such as `.DS_Store` or `.git`) are left out. The page lists the files at once and
reads a file's content only when the core needs it. Then:

1. The core **describes** every file, one after another: what it is (an `.xye` file, a FullProf free-format file, a
   buffer, a plot set file …, or "?" for something no reader recognises). The tree fills in each file's kind as the
   answers come. For a hundred large synchrotron files this takes about a minute.
2. The core **proposes** the plot sets the folder's data files make (§7).
3. In Chromium, **live mode** starts watching the folder (§15).

**Open Files…** (bar or start page). Choose one or more files. They join the folder already shown, or make a new one
if none is, and they **open in Plots at once** (§6 says as what).

**Drag and drop.** Drop anything from the Finder (or Explorer) anywhere on the page:

- **a folder** becomes the folder shown, exactly as with Open Folder…. In Chromium a dropped folder is also watched by
  live mode;
- **loose files** join the folder shown and open in Plots at once, as with Open Files…. Dropping the same files again
  opens them again.

A page holds **one folder at a time**. Opening another folder replaces the tree. Plot sets already open stay in their
tabs.

## 5. The folder's tree

Each row of the tree shows a file's name, its kind as the core says it (`xye`, `free_steps`, `buffer_file`, `plot_set`,
or `?`) and its size. Hover a row for the file's path, its size and the core's full description.

**folders | kinds** (the switch on the Files header) chooses how the tree is listed:

- **folders**: as on disk. The first level of sub-folders is open, and a double-click (or click) on a folder opens or
  closes it.
- **kinds**: one group per kind as the core describes the files: the data files first, then the project's files (plot
  sets, buffers), then what is not recognised. This makes it easy to see at a glance which files the core can read.

The choice is kept for your next visit. The header also counts the files and their total size.

**Choosing files.** A click chooses a file, **⌘-click** (Ctrl-click on Windows and Linux) adds or removes one, and
**⇧-click** chooses a range. **Open as Plot Set** then opens the chosen files; it says "Open 12 as Plot Set" when 12
are chosen.

**Opening directly.**

- A **double-click** on a file opens it in Plots.
- The **set** button on a folder's row opens all the data files under that folder as one plot set. This is the
  quickest way when a folder holds exactly one series.

## 6. What opens as what

Whatever you open (a double-click, Open as Plot Set, a folder's "set", files dropped or picked), the page asks the core
what each file is and routes it accordingly:

| What is among the files | What opens |
|---|---|
| a **plot set file** (`.drv`) | the set it describes, with its members, variables, view settings, clips and regions. Its member files are taken from the folder by the names the `.drv` keeps (see below). |
| a **WinPLOTR buffer** (`.buf`) | the buffer's set, its files taken from the folder by the names it lists. |
| **data files** only | one new plot set with those files as its members, in name order. The core reads each member's variable from the file if it can: a temperature in its header, or in its name (`_300K`, `T300`, `T0900C` …), or an elapsed time (`t16530s`). |
| something else (a `.pcr`, an image, an unrecognised file) | left out. The status line names it and says what the core says it is. |

**A `.drv` or `.buf` needs its files where it names them.** A plot set file written by the app names its members
relative to the project folder, for example `Data/scan_001.xye` for a `.drv` kept in `Plots/`. So open the **project
folder** (the one that holds `Plots/` and `Data/`), not the `Plots/` folder alone, and then open the `.drv` from the
tree. A member whose file is missing is named in the status line.

**The same set file opened twice** brings its tab forward instead of opening a second copy.

## 7. Sets in this folder: the core's proposals

When a folder is opened, the core reads all its data files and proposes the plot sets they make. The section **Sets in
this folder** lists them, each with an **Open** button, a title and the number of members. Hover a proposal for the
core's kind and exposure class and for the files it leaves out, with the reason.

The kinds of proposal:

- **By temperature**: one scan per distinct temperature, ordered. Where several scans share a temperature, the one of
  the longest exposure is taken. Examples: *"19 scans by temperature, 4 to 30 K"* (the Tb₂Pd₂Sn
  neutron patterns of FullProf's WinPLOTR tutorial); *"25 scans by temperature, 900 to 1300 °C (the longest exposure
  at each temperature, of 3 exposure classes), rebinned to 0.01°"* (an in-situ synchrotron run of 101 scans).
- **Against time**: the whole run, each scan at its time since the first. Example: *"100 scans against time, 10.036 h
  from …001… (3 exposure classes), rebinned to 0.01°"*.
- **A hold**: a run of three or more scans at the same temperature, against time since the hold's first scan. Other
  holds of the same class are named in the title. Example: *"the hold at 1270 °C: 39 scans against time, 166.433 min
  …"*.

Where the variable comes from:

- **The temperature** is taken from a file's header or its name.
- **The time** is taken from the name (`t00000s`) or, failing that, from the files' modification times.
- **Exposure classes** separate scans counted differently, such as long and short scans. The core finds them from the
  counts' scale, the median of σ²/y above 2°. Scans of one class differ in that scale by less than 10 %. When there is
  more than one class, the section says how many.

**Rebinning on opening.** For large scans (tens of thousands of points), the proposal says *"rebinned to 0.01°"*. The
core then averages the points onto a 0.01° grid as the set opens, so that the map and the overlay draw in about a
second instead of minutes. The files themselves are not touched.

Files that fit no proposal (a background scan with no temperature in its name, for example) are listed under "Left out
of every set" with the reason.

**Open** opens the proposal as a new plot set. A proposal opened this way is also remembered for **live mode**: when
new scans of the same kind and class appear in the folder, they are added to it (§15).

> **Important in this version: proposals against time.** A proposal opens from its files, and each member's variable is
> read again from the file, not from the proposal. Where a file's name gives both a temperature and a time (as in
> `…_T1270C_t16530s_…`), the members open **by temperature**, with the time only as a second variable. A
> **by-temperature** proposal is therefore right. A **hold** or a **run against time** opens with every scan of a hold
> at the same temperature, so its map collapses onto one line. Until this is fixed, open sets against time from a
> plot set file (`.drv`) that types each scan's time (the site's recipe *A hold: growth curves* shows how). The status line says when a
> proposal's values were not carried.

## 8. A plot set in its tab

Each plot set opens in a tab of its own in the views area. Tabs can be dragged to rearrange them, and closed with their
×. The **?** at the top right of the views opens this manual at the section of the view shown. Above the plot, the tab
shows:

- **The title strip**: the set's title, its number of members, and "· normalised" or "· in d" when a view rule is in
  force (§12).
- **The view tabs**: Overlay, Stack, Projection, Map, Difference Map, Surface, Series (§10).
- **The toolbar** of the view shown. Only the controls that apply to the view are visible:

| control | in | does |
|---|---|---|
| **y** linear / √ / log | Overlay, Stack, Projection, Series | the intensity axis's scale |
| **Linked zoom** | Overlay | a second, zoomed plot beside the first, with a difference strip (§10) |
| **Δ** and its rule | Stack, Projection | the steps between members (Set Displacement…) |
| **Draw window** | Overlay, Stack | the next drag across the plot draws a window of 2θ (§12) |
| **Figure look** | all but Surface | the view as a figure: white, open axes, ink text, a scale bar for a stack or a projection |
| **colour** linear / √ / log | Map, Difference Map, Surface | the colour scale |
| **Range…** and its clip | Map, Difference Map, Surface | type the colour range (§11) |
| **Cursors** | Map, Difference Map | the row and column strips (§11) |
| **Range from view** | Map, Difference Map | the colour range follows the zoom (§11) |
| **Save region…** | Map, Difference Map | name the rectangle shown (§11) |
| **y: shared / y: own** | Overlay, Stack, maps | one y range for the whole set, or each member's own |
| **of** the map / the difference map | Surface | which map the surface shows |
| **quantity**, **background**, **+ window** | Series | what the series follows (§12) |

Below the plot, a short note explains what the view shows and mentions anything left out or drawn apart.

## 9. The members and the member

When a plot set is shown, the inspector holds its panels:

**Plot sets** ("in memory, not saved"): one row per open set. Click a row to switch to that set.

**Windows** (2θ, °): the set's windows of 2θ (§12). Click a row to choose the window, double-click to rename it, and
× removes it.

**Regions**: the named rectangles of the map (§11).

**Members**: one row per member, with its variable's name and unit in the header. Each row shows:

- a **checkbox**: whether the member is *chosen*. Chosen members are what Solo, Dim, Set Reference, Difference of
  Members and the inspector act on;
- an **icon** for the file's kind (a letter, for example F for FullProf formats, ≡ for a column table);
- a **colour swatch**: its colour in the overlay, along a dark-to-bright ramp in the set's order. Grey means a member
  with no value of the variable;
- the **label** (the file's name), with a badge "ops N" if the Mac app's Data module keeps operations on the file;
- the **value of the variable**, marked "typed" if it was typed rather than read, "apart" if it has none, or "out" if
  it is left out of the view;
- an **opacity** slider (0 to 1). Turn a member down to see others behind it.

Mouse on the member rows:

- a **click** chooses one member; **⌘-click** or **⇧-click** adds or removes a member;
- **⌥-click** (Alt-click) *solos* the member: the others are dimmed until you ⌥-click again;
- a **double-click** on a member without a value opens **Set Temperature…**;
- **hovering** shows the member's facts: file, value and its source, format, points, range, wavelength, counting time
  and exposure class, monitor, the file's date, and the normalisation factor.

The control row above the members:

- **order**: the set's own order, by the variable, by the file's time, or by name;
- **all**, **none**: choose every member, or none;
- **class…**: choose the members of one exposure class, as the files state their counting time ("10 s: 4 members");
- **range…**: choose the members whose variable lies between two values, ends included.

Notes at the foot of the panel say, for example, "no monitor: raw counts", or which members are "Drawn apart (no T in
the file)".

**The member** (fold it open): with nothing chosen it describes the set: whether all members share one grid, the
counts, and how the view is set. With a member chosen it shows:

- the variable's value with its source, and other candidate values ("also: …");
- the Mac app's Data operations on the file;
- how the member is drawn;
- **What the file says**: format, axis, points, range, step, wavelength, monitor, counting time, title;
- the **Actions**: Set Temperature…, Solo Member, Set Reference, Set Wavelength… (and Open in Data, which belongs to
  the Mac app and does nothing in the browser).

**Every value keeps its source.** Set Temperature… asks for the value **and where it comes from** (a logbook, a
thermocouple file, "nominal set point"). The page refuses a value without a source: the source is saved with the value.
**Read Variable From…** reads the values again from the files' names or from the files themselves.

## 10. The views

All views share the same mouse controls (§18). In short: drag a box to zoom, double-click to see everything, ⇧-drag
to pan.

<a id="help-overlay"></a>

### Overlay

Every member drawn as a curve on the same axes, coloured along the set's order. Hover a curve to thicken it.
⌥-click a curve to solo it, and ⌥-click again to show all.

With **y: own**, the y range fits the chosen members. With **y: shared**, it is the core's extent of the whole set, so
that the scale does not jump as members are switched.

**Linked zoom** splits the overlay in two:

- The left plot shows the whole pattern. Drag across it to choose the 2θ window that the right plot shows.
- Zooming or panning the right plot moves the band on the left.
- Under both plots a **difference strip** shows the first chosen member less the reference member, or less the second
  chosen member. It is the quickest way to see what changed between two scans.

Large sets are drawn in slices. The status line counts them ("Overlay: drawing 60 of 100 members…") while the page
stays responsive.

<a id="help-stack"></a>

### Stack

Each member raised by a step above the one before and labelled at its baseline with its value ("1270 °C"), or with its
name and "(no T)". The core proposes the step from the data; **Δ** or **Set Displacement…** sets your own:

- **rule**: *step* (the same Δy for each member) or *variable* (Δy proportional to the variable, so that the spacing
  shows the temperature or time steps);
- **Δy**: the vertical step, in the y axis's units;
- **Δx** (°): a horizontal shift per member, which gives a slanted stack. 0 gives a plain stack.

The y label reads "counts + each member's step": the numbers on the y axis hold for the first member only. In Figure
look a scale bar replaces the y axis for that reason.

<a id="help-projection"></a>

### Projection

WinPLOTR's pseudo-3D view: each member drawn behind the one before, shifted by Δx and Δy, its hidden parts removed,
with an oblique axis carrying the variable. The view stays empty until you give it steps with **Δ** (Set the
Projection's Steps); it keeps steps of its own, separate from the Stack's.

### Map

Intensity as colour, with 2θ across and one row per member up the side at its value of the variable. Each row spans
halfway to its neighbours, so uneven steps in temperature or time show as rows of different heights. This is the view
for following peaks through a series: a phase appearing, a peak moving, a transition. §11 covers it in detail.

<a id="help-difference"></a>

### Difference Map

The map of each member **less a reference member**. Choose a member and **Set Reference** first (until you do, the
warmest member placed is the reference). Values above zero and below zero are told apart by the colour ramp. The log
scale is not offered, since differences can be negative: log falls back to signed √.

<a id="help-surface"></a>

### Surface

The map, or the difference map (**of** on the toolbar), as a 3D surface:

- drag to turn it;
- the wheel zooms;
- a double-click starts again;
- the pointer reads the cell under it.

The map's colour range also limits the height. **Copy Surface as Image** (menu) copies it to the clipboard as a PNG.

### Series

The quantity chosen for each **window** of 2θ, against the variable, with ±σ error bars. Points below 2σ are drawn
open. Each window has its own colour and marker. With no window yet, the view tells you to add one (§12).

<a id="help-map"></a>

## 11. The Map and the Difference Map in detail

**Reading the map.** The colour bar at the right gives the scale: "counts", "√counts" or "log₁₀ counts". In the
Difference Map it says what was subtracted ("less 1300 °C (name), counts, signed √ scale"). The readout follows the
pointer and shows:

- 2θ;
- the variable;
- the member's label;
- any second variable of that row, with its source;
- the cell's value ("no value" where a cell is empty).

**The right-hand axis.** When the members carry a second variable (for a set against time, the temperature of each
scan), the right-hand axis labels each row with it. You can then read the time on the left and the temperature on the
right.

**Colour scale** (**colour** on the toolbar): linear, √ or log. √ is often best for diffraction: it compresses the
strong peaks less than log and still shows the weak ones.

**The colour range: the clip.** The core gives the full range of the data. Strong peaks usually take the whole colour
ramp, leaving the weak features in one colour. To go down to the details, **clip** the range:

- **Drag** either of the two marks on the colour bar (the pointer changes to a double arrow). Values beyond the clip
  are drawn saturated. The two marks stay in order and inside the core's range.
- **Range…** (or Plots ▸ Set Colour Range…) to **type** the two values, *From* and *To*. An empty field keeps the
  core's end at that side, and both empty restores the core's range. If the values lie wholly outside the data, the page
  refuses them and says why. If they hide most of the map, the core's warning is added to the note under the map
  ("The core says: …").
- **Double-click the colour bar** (or Reset Colour Range) to return to the core's range.

The clip is kept per view and per scale. It is saved with the set (§13) and carried in the page's address (§16).
Changing the scale starts again from the core's range for that scale.

**Range from view.** Switched on, the colour range follows the zoom: whenever you zoom, the core gives the range of the
rectangle shown, and the colours spread over that rectangle alone. Switched off, the earlier clip comes back.

**Zoom.**

- **Drag a rectangle** to zoom both axes.
- **The wheel over the left margin** zooms the variable's axis only, and **over the bottom band** zooms 2θ only.
- **⌘-wheel** or a pinch zooms about the pointer.
- **A double-click** on the map shows the whole map.

The zoom is kept when you switch views and back.

**Cursors** (the toolbar's Cursors, or Row and Column Cursors):

- the **row strip** under the map shows the pattern of the row under the pointer, on the same 2θ axis as the map. In
  the Difference Map it shows the difference row;
- the **column strip** beside the map shows the intensity at the pointer's 2θ against the variable: a peak's growth
  through the series as a curve;
- **click the map** to **pin** the cursors at a point ("The cursors pinned at …"), then move freely; click again to
  release them.

Turning Cursors on switches to the Map if another view is shown. The strips' y range follows **y: shared / own** as in
the Overlay. Cursors are not drawn in Figure look.

**Regions: named rectangles.** A region keeps a place of the map worth coming back to, such as a peak's onset or a
window where a phase grows.

1. Zoom to the rectangle (and, if you like, switch Range from view on, or clip the colours).
2. **Save region…** and give it a name (r1, r2 … by default). The core computes the rectangle's block of rows and
   columns, and the status line names it. With Range from view on, the region keeps its colour range too.
3. The region is listed under **Regions** in the inspector. **Click** it to come back to it: the map zooms to the
   region and applies its own colour range.
4. **Double-click** a region (or Rename Region…) to rename it, and **×** (or Remove Region) to remove it.

Two regions cannot share a name. Regions are saved with the set (§13). A set file made in the Mac app brings its regions
with it.

**Members drawn apart.** A member with no value of the variable has no place on the map's axis. It is drawn as a
separate row and named in the note ("Drawn apart, no T in the file: …"). Give it a value with Set Temperature… to place
it.

<a id="help-series"></a>

## 12. Normalising, converting x, windows and the Series view

These change **the view**, never the files ("A view: the files untouched", says the status line).

**Normalise…** puts the members on a common scale. It offers:

- **none**: the counts as read;
- **the monitor**: each member divided by the monitor its file states;
- **the counting time**: each member divided by the counting time its file states;
- **the maximum**: the highest point of each member set to the same value;
- **the same area**: each member's area set to the same value, between *From* and *To* 2θ;
- **a peak's integral**: one peak's integral above its background, between *From* and *To*;
- **a reference member**: every member scaled to the reference (choose it with Set Reference first).

The fields are *Normalised to* (the value), *From, 2θ (°)*, *To, 2θ (°)* and *Background* (auto, or a line between the
ends). Members the rule cannot normalise (no monitor, for example) are left out of the view and listed. The title strip
then shows "· normalised".

**Convert X…** draws the set against another variable:

- the files' own x;
- 2θ (°);
- d (Å);
- Q = 4π sin θ/λ (Å⁻¹);
- sin θ/λ (Å⁻¹);
- s = 1/d (Å⁻¹).

Every conversion but the first needs each member's wavelength. If some files state none, a chooser offers *a λ typed
now* (with its source) or *none: leave them out*, listed. **Set Wavelength…** sets the λ of the chosen members, with
its source.

**Windows** are ranges of 2θ followed through the set: a peak, a group of peaks, a band. There are three ways to make
one:

- **+ window** on the Series toolbar, or **Add Window…**: a *Name*, *From (°)* and *To (°)*;
- **Draw window** (Overlay, Stack), then drag across the plot. The window is named w1, w2 …;
- the windows a `.drv` already carries.

To work with a window:

- drag a window's band edges on the plot to move them;
- click the band to choose the window;
- rename or remove it in the inspector or the menu.

**The Series view** plots, for every window, one **quantity** against the variable:

- **integral above the background** (counts·°), with σ;
- **centroid** (°), with σ;
- **strongest peak: position**, **height**, **width**;
- **background at x**: the background's value at the 2θ positions given with Background at….

The **background** under a window is either **automatic** (the arPLS baseline GSAS-II uses) or a **line between the
ends** of the window. The integral of a peak's window against temperature or time is the usual way to follow a phase's
growth.

## 13. Saving and exporting

**Save Plot Set** (Plots ▾) writes the set as a **plot set file**, `<name>.drv`. In the browser it is **downloaded**
to your Downloads folder: nothing is written into the folder you opened. The `.drv` keeps:

- the members (by their file names) and the variable, with every typed value's source;
- the order;
- the stack's and projection's steps;
- the y and colour scales;
- the reference member;
- the windows, the opacities, the clips and the regions;
- the y interval.

**The 2θ interval of the analysis.** A set file may say `x_interval = [2.0, 20.0]` among its first lines. Every member
is then cut to that interval of 2θ as the set opens: every view, window, map, colour range and difference is that of
the cut set, the files untouched, and a window outside the interval is refused. Use it to leave out a beamstop's first
degree, or a range of no interest, once for the whole analysis. In this version it is typed in the set file (a text
editor will do); a control in the page comes later (Appendix C).

It is a plain text (TOML) file, read by the Mac app and the command line `drs` as well as by this page.

To open it again later:

1. Move the downloaded `.drv` into the folder's `Plots/` sub-folder, or next to the data files. It must be placed where
   its member names resolve: §6.
2. Open the folder.
3. Double-click the `.drv`.

Or drop the `.drv` on the page while the folder is open.

**Save Plot Set As…**: the same, with a new name.

**Copying.** The bottom table's **Copy** (or ⌘C on selected rows) copies its rows as tab-separated text, ready to paste
into a spreadsheet. **Copy Surface as Image** copies the surface as a PNG.

**Viewer Page… (Dev)**: downloads a single HTML file, `<set> map viewer.html`, holding the map viewer and the set's
maps, members, clips and regions. It opens in any browser, offline, without this page or the files. It is marked Dev
because it is a developer's tool; it has not been checked in the browser host in this version.

**Not in the browser in this version** (they belong to the Mac app or need a later core build; §20):

- figure export (PDF, SVG, PNG);
- Plot Set as CSV…;
- Send to Publish….

A screenshot of the view in Figure look is the stop-gap for a figure.

## 14. The Plots menu, item by item

The **Plots ▾** menu lists the module's commands in the app's order. A greyed item explains itself on hover. Marks:
**✓** a toggle that is on; **(app)** an item that belongs to the Mac app and does nothing useful in the browser;
**(later)** an item this browser core does not serve yet: it answers in the status line *no command "…" in the
browser's core (drs-wasm 0.81.0): the app's core has it, this build not yet*.

**Opening**

| item | does |
|---|---|
| New Plot Set… | an empty set, then Add Files… for its members (choose them in the tree first) |
| Open Plot Set… | opens the `.drv` chosen in the tree |
| Open Buffer… | opens the WinPLOTR `.buf` chosen in the tree |
| Open CSV as a Set… | the chosen files as one set (what Open as Plot Set does) |
| Add Files… | adds the files chosen in the tree to the set shown |
| Open Sequence (.seq)… | FullProf's sequential results against their external parameter **(later)** |
| Open Results Series… | a series of refinements against T **(later)** |
| Send Series to Publish… | **(app)**, always greyed here |

The open panels of the browser are the tree: **choose the files in the tree first, then the command.** With nothing
suitable chosen, the status line says what to choose.

**Saving and output**

| item | does |
|---|---|
| Save Plot Set / Save Plot Set As… | downloads the `.drv` (§13) |
| Reload Members | reads the members' files again **(later)** |
| Plot Set as CSV… | wide or long CSV of the set **(later)** |
| Send to Publish… | **(app)**: saves (downloads), then says Publish is not in this app |

**Views and map tools**

| item | does |
|---|---|
| Show Overlay / Stack / Projection / Map / Difference Map / Surface / Series | switch the view (✓ on the one shown) |
| Copy Surface as Image | the surface as a PNG on the clipboard |
| Linked Zoom | the overlay's linked pair (✓) |
| Set Colour Range… / Reset Colour Range | type the clip / back to the core's range (map views and Surface) |
| Row and Column Cursors | the strips (✓) |
| Range from the View | the clip follows the zoom (✓) |
| Save Region… / Rename Region… / Remove Region | the map's regions |
| Y Interval Adapted | ✓ = "y: own", off = "y: shared" |
| Viewer Page… (Dev) | the self-contained map page, downloaded |

**Choosing members**: Choose All Members, Choose No Member, Choose by Exposure Class…, Choose by a Range of the
Variable…, List Order…

**Figure**

| item | does |
|---|---|
| Figure Look | the view as a figure (✓) |
| Set Scale Bar… | the scale bar's length in a stack or projection figure (empty: the core's round step) |
| Figure Size… | journal and thesis widths from Publish's presets **(later)** |

**Normalisation and x**: Normalise…, Convert X…, Set Wavelength… (§12).

**Order**: Reverse Order, Sort by Variable, Read Variable From… (the file's name or the file itself work; *a
temperatures file* is **(later)**).

**Windows**: Add Window…, Draw Window (✓), Rename Window…, Remove Window, Remove the Windows, Background at….

**Opacity**: Solo Member, Dim Others (others at 22 % opacity), Restore All.

**Per member**: Set Displacement…, Set Reference (also switches to the Difference Map), Set Temperature….

**Derived**: Sum of Members, Average of Members, Difference of Members (exactly two chosen: the first less the second).
The core computes the pattern with its σ; the browser downloads it as `<set> sum.xys` (average, difference) and the
status line says so.

**Open in Data**: **(app)**.

<a id="help-live"></a>

## 15. Live mode: a folder still being written

At a beamline or a diffractometer, scans keep arriving in a folder. In **Chromium**, with a folder opened by **Open
Folder…** or dropped, the page can watch it:

- **Live** (checkbox, on by default): read the folder again on a timer;
- **every**: 2, 5, 10, 30 or 60 s (10 s by default);
- **Now**: read it at once.

Your choice is kept for the next visit. At each reading:

- **new files** appear in the tree and are described by the core;
- if files were added, the core **proposes again**, and every set **opened from a proposal** grows with the new scans
  the core now proposes for the same kind and class. They are added to the set as Add Files… would, and the map is
  computed again. Such sets are marked **live** in the proposals list;
- a file that **vanishes** is named in the live line (and in the status line if it belongs to a set), but it **stays in
  its set**: nothing is ever removed from a set by itself.

The live line under the controls says what is happening: "Reading the folder every 10 s; 1 live set. …".

Notes:

- Only sets opened **from a proposal** grow. A set opened from chosen files stays as it was opened.
- Re-proposing reads every data file again. For a hundred large scans this takes tens of seconds, which is why it is
  done only when files have been added.
- Browsers slow the timers of a tab in the background (to once a minute after five minutes). Keep the tab in front for
  quick updates.
- Firefox and Safari do not give the page a lasting handle on a folder, so live mode is not offered there.

<a id="help-address"></a>

## 16. The view in the address: links and reloads

The page keeps the view shown in its address, after the `#`, updating it about once a second:

- the folder's name;
- where the set came from: its files, its `.drv` or `.buf`, or the core's proposal by kind and class;
- the view;
- the colour scale;
- the map's zoom;
- the clip;
- the region chosen.

For example:

```text
…/diff-ref-struct-web/#folder=1_Tb2pd2sn&set=buf%3ATb_vs_T.buf&view=map&scale=log&zoom=24%2C36%2C4%2C30&clip=3.2%2C4.5
```

The address holds **only names and numbers already on screen**, never a file's content. So:

- **Reload** the page, open the same folder again, and the view comes back as it was: the set opened, the view, the
  scale, the zoom, the clip and the region.
- **Send the link** to a colleague who has the same files: when they open their copy of the folder, they see the same
  view. If files are missing, the page names them and says to open that folder.
- **Bookmark** views you come back to.

The zoom is applied once the core's map has been drawn. You may see the whole map for a moment first.

**Cite this view** (Help ▾) copies one line with the tool's version, the core's, the date and this address, for a lab
book or a paper's methods: *J. M. Igartua, diff-ref-struct Web 1.0.0 (core 0.81.0), view of 26 October 2026:
https://…/#folder=…*. With the same files, the line brings the same view back.

## 17. Offline and installed

On the first visit the browser keeps the page's own files: the page, the module, the core in WebAssembly, about 36 files
and 4.4 MB. From the next visit on, the page opens **without the network**. The line under the inspector says what is
kept: "Offline: 36 files, 4.4 MB kept by this browser (web …, core 0.81.0)". It never keeps your data files.

**The guide offline too.** On the same first visit the browser also keeps the guide's pages and this manual's PDF
(not the screenshots, which it keeps once you have seen them), so the manual opens at a beamline without network.

**Install** the page as an app from the browser: in Chrome, the install icon in the address bar or ⋮ ▸ Cast, save and
share ▸ Install page as app; in Safari, File ▸ Add to Dock; in Edge, ⋯ ▸ Apps ▸ Install. It then opens in its own
window, from the Dock or the Start menu.

**Updates.** When a new version is published, the browser fetches it in the background and replaces the old one. A tab
already open may keep running the old version until it is reloaded. If the version shown (start page, offline line)
is older than expected, **close every tab of the page and open it again**.

## 18. Mouse, keyboard, appearance

**On any plot**

| action | does |
|---|---|
| drag a box | zoom both axes |
| ⇧-drag, middle-drag, or Space + drag | pan |
| two-finger scroll | pan (a mostly vertical scroll zooms x) |
| ⌘-wheel or pinch | zoom about the pointer |
| wheel over an axis's margin (maps) | zoom that axis only |
| double-click | fit the whole plot |
| double-click on a colour bar | the core's colour range again |
| drag a colour bar's mark | clip the colour range |
| click (map with Cursors on) | pin or release the cursors |
| ⌥-click a curve (Overlay) | solo it, and again to restore |
| drag on the plot with Draw window on | a new window |
| drag a window's edge | move it |

**In the tree and the members**: click, ⌘-click (Ctrl-click), ⇧-click to choose; double-click to open (tree) or to
set a value (member without one); ⌥-click to solo (members).

**In the bottom table**:

- click a header to sort (▲, ▼, off);
- drag a column's edge to resize it;
- ↑ ↓ move the selection;
- ⌘C copies the selected rows;
- **Find…** searches the table, and **Columns** shows or hides columns (time and title are hidden by default).

**Elsewhere**: Escape closes the Plots menu; ← → on the focused separator set the inspector's width.

**Keyboard shortcuts.** The Plots commands have no keyboard shortcuts in the browser: use the menu. In particular
**⌘S is the browser's "Save page"**, not Save Plot Set. The status line's "Save (⌘S)" is the Mac app's wording.

**Appearance**: the ◐ button cycles auto, light and dark.

**Accessibility**: the tree, the menu, the separator and the status line carry their roles for screen readers. The
status line is announced politely as it changes.

## 19. Messages you will meet, and what to do

| the status line says | meaning, and what to do |
|---|---|
| *Nothing to open as a plot set: x.pcr (…)* | the files chosen are not data files. Choose data files, a `.drv` or a `.buf`. |
| *… choose the files in the folder's tree first (click, ⌘-click, ⇧-click), then the command again* | a menu command needs files: choose them in the tree, then the command. |
| *Plots: N members, k without a T (drawn apart) …* | opened; k files gave no value of the variable. Set Temperature… places them. |
| *no command "…" in the browser's core (drs-wasm 0.81.0): the app's core has it, this build not yet* | the command needs the Mac app or a later web version (§20). |
| *The core says: …* (under a map) | the colour range typed hides much of the map; the core says how much. |
| *To lies above from.* / *The pair … lies wholly outside the core's range* | the clip typed is reversed, or outside the data. |
| *Show the map or the difference map first: a region is the rectangle it shows.* | switch to a map view before Save Region…. |
| *A region needs a name of its own …* / *A window needs a name of its own …* | the name is taken. |
| *Say where the value comes from: it is kept with it.* | a typed temperature or λ needs its source. |
| *Choose two members: the first chosen less the second.* | Difference of Members needs exactly two chosen. |
| *Downloaded as x.drv (your browser's Downloads): nothing is written into the folder …* | Save worked: the file is in Downloads. |
| *The address keeps a view of …'s files; n of them are not in the folder open …* | the link or reload names files of another folder: open that one. |
| *Live: … vanished from the folder; the set keeps it as read* | a scan disappeared; the set keeps the data it read. |
| *x.drv is open already: its tab brought forward.* | the set was open; nothing new opened. |
| *…: opened from its N files. The M values its set file types … are not carried* | a proposal opened from its files; the variable was read again from the files (§7). |

If the page seems stuck on a very large set, wait for the status line. A run of 101 synchrotron scans of 27 000 points (112 MB) takes
about a minute to describe and about 40 s to propose. If nothing changes for several minutes, reload: the address
brings the view back.

## 20. Limits of this version

What the browser does **not** do in diff-ref-struct Web 1.0.0 with core 0.81.0, and why. Each limit is on the
roadmap (Appendix C) with what is decided about it.

- **Proposals against time open by temperature** (§7). A hold's or a run's members take the temperature in their
  names as their variable, and the time becomes the second variable. The map of a hold collapses. Use a `.drv` that
  types the times.

- **No figure export** (PDF, SVG, PNG): the Mac app's Export ▸ Plot. Use Figure look and a screenshot, or Viewer
  Page… for an interactive page.
- **No CSV export** (Plot Set as CSV…), no Figure Size…, no Reload Members, no Open Sequence (.seq)…, no Open Results
  Series…, no *temperatures file* in Read Variable From…: the browser's core does not serve these yet, and the status
  line says so.
- **No Data operations, refinement, Publish or projects**: those are the Mac app's.
- **Text files only**: binary formats (Bruker RAW and BRML, Rigaku RASX, NeXus/HDF5, zarr) are not read by the
  browser in this version, even though the core knows them (§21).
- **The strips under the map** (Cursors) may sit slightly off the map's frame in some window sizes: a peak in the row
  strip a degree or so left of its column, the strip's y label saying its scale twice ("counts (log) (log)"). It was
  seen on the public page on 26 October 2026; in other captures the frames match. Read 2θ on the map itself. A fix is
  under way.
- **No control in the page for the 2θ interval**: type `x_interval = [from, to]` in the set file (§13).
- **Coming next**: that control, and the map's default colour scale linear in place of log.

## 21. The file formats the core reads

The core recognises a file by its **content first**, its extension second, and says which reader it chose (hover the
file in the tree). Text formats, read in the browser:

| kind (as the tree shows it) | format | extensions |
|---|---|---|
| `xye` | x, y, σ columns (Topas .xye, GSAS-II's rules) | .xye .xy .dat .txt |
| `xys` | x, y, σ columns with six header lines (FullProf XYDATA) | .dat .xys .xy |
| `qye` | Q, y, σ columns, wavelength in a comment | .qye |
| `table` | table of columns: CSV, TSV, blanks; delimiter, decimal mark, header and columns recognised | .csv .tsv .txt .dat .xy |
| `free_steps` | FullProf free format: Ti, step, Tf, then the intensities | .dat |
| `ill_d1a_d2b_old`, `ill_d1b_old`, `ill_d1b_d20`, `ill_d2b`, `ill_d4_d20l` | ILL instruments (D1A, D2B, D1B, D20, D4) | .dat |
| `psi_dmc` | PSI DMC, HRPT | .dat |
| `brookhaven` | Brookhaven synchrotron (y and σ on alternate lines) | .dat |
| `two_axis` | general two-axis instrument | .dat |
| `variable_time` | X-rays with variable counting time | .dat |
| `dbws` | DBWS synchrotron format | .dat |
| `gsas` | GSAS powder data (STD, ESD, FXY, FXYE, ALT) | .gsas .gss .gsa .gda .fxye .fxy .raw .dat .xra .cwn |
| `isis_tof_banks` | ISIS TOF multibank | .gss .dat |
| `socabim_uxd` | Socabim UXD (Bruker DIFFRAC-AT text) | .uxd |
| `panalytical_xrdml` | PANalytical XRDML | .xrdml .xml |
| `rigaku_ras`, `rigaku_txt` | Rigaku RAS, Rigaku text export | .ras .txt |
| `fit2d_chi` | Fit2D or pyFAI .chi (2θ or Q) | .chi .qchi |
| `powder_cif` | powder CIF (pdCIF) | .cif |

The project files it opens: plot sets (`.drv`) and WinPLOTR buffers (`.buf`).

Recognised by the core but **not read in the browser** (binary): Bruker RAW v1–v4, Bruker BRML, Rigaku RASX,
NeXus/HDF5, MIDAS zarr. Use the Mac app or convert them to `.xye`.
