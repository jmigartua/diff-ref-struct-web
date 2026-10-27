# diff-ref-struct Web — quick reference

For diff-ref-struct Web 1.1.0 with core 0.82.0; every control in full in [the manual](user-manual.md).

## Bring files in

- **Try the example** (start page): a heating run computed by the core.

- **Open Folder…** or drop a folder: the tree, the core's description of every file, the core's proposals.
- **Open Files…** or drop files: they open in Plots at once.
- A `.drv` or `.buf` opens as itself, with its files taken from the folder it names them from.
- Nothing is uploaded; text files up to 64 MB.

## The tree

- **folders | kinds**: as on disk, or by what the core says the files are.
- Click, ⌘-click, ⇧-click to choose; **Open as Plot Set**.
- Double-click a file to open it; **set** on a folder opens its data files.

## Sets in this folder

- The core proposes: **by temperature**, **a hold** against time, **the run** against time.
- Large scans are rebinned to 0.01° on opening; the files are untouched.
- Proposals against time open by temperature in this version: use a `.drv` that types the times.

## Views

- **Overlay**: ⌥-click a curve to solo it; **Linked zoom** gives a zoomed plot and a difference strip.
- **Stack**: Δ sets the step (*step* or *variable*), Δy and Δx.
- **Projection**: give it steps with Δ first.
- **Map**: colour linear, √ or log; clip, cursors, regions.
- **Difference Map**: Set Reference first; signed √ in place of log.
- **Surface**: drag to turn it, wheel to zoom, double-click to reset.
- **Series**: windows' integral, centroid or peak against the variable.

## On any plot

- Drag a box: zoom. Double-click: everything.
- ⇧-drag, middle-drag or Space-drag: pan.
- ⌘-wheel or pinch: zoom about the pointer.
- Wheel over an axis's margin (maps): that axis only.

## The map's colour

- Drag the bar's two marks: **clip**. Values beyond are saturated.
- **Range…**: type From and To; an empty field keeps the core's end.
- Double-click the bar: the core's range again.
- **Range from view**: the range follows the zoom.

## Cursors and regions

- **Cursors**: row strip below, column strip at the right.
- Click the map: pin the cursors; click again: release them.
- **Save region…**: names the rectangle shown; click it in Regions to come back; double-click renames; × removes.

## Members

- Checkbox: chosen. ⌘- or ⇧-click: add. ⌥-click: solo.
- Opacity slider per member.
- **order**, **all**, **none**, **class…**, **range…**.
- Every typed value keeps its source (Set Temperature…, Set Wavelength…).

## Windows and the Series

- **+ window** or **Add Window…**: name, from, to (°).
- **Draw window**, then drag across the Overlay or the Stack.
- Background: *automatic* (arPLS) or *line between the ends*.

## View rules (the files untouched)

- **Normalise…**: monitor, counting time, maximum, area, a peak, a reference.
- **Convert X…**: 2θ, d, Q, sin θ/λ, 1/d (λ typed with its source where a file has none).

## Keep and share

- **Help ▾ ▸ Cite this view**: one line with the versions and the address.
- **Help ▾ ▸ Report a problem… / Ask for a feature…**: the form, filled in for you.

- **Save Plot Set**: the `.drv` is downloaded; put it in the folder's `Plots/`.
- The address keeps the set, view, scale, zoom, clip and region: reload it or send it.
- ⌘S is the browser's Save Page, not Save Plot Set.
- `x_interval = [from, to]` in a set file: every member cut to that 2θ interval on opening.

## Live mode (Chrome, Edge, Arc)

- **Live**, every 2–60 s, **Now**.
- Sets opened from a proposal grow with the new scans; nothing is removed by itself.
- Keep the tab in front: background tabs are slowed.

## Not in the browser yet

- Figure export, CSV export, Figure Size….
- Sequences (.seq), results series, Data's operations.
- Binary formats (Bruker RAW and BRML, Rigaku RASX, HDF5).
