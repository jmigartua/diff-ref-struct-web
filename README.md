# diff-ref-struct Web

Series of powder diffraction patterns in the browser: overlays, stacks, WinPLOTR-style projections, maps of intensity
against 2θ and temperature or time, difference maps, surfaces and growth curves of chosen peaks. The numbers are
computed by the diff-ref-struct core, compiled to WebAssembly and run in your browser; your files are read there and
never uploaded.

- **The app**: <https://jmigartua.github.io/diff-ref-struct-web/>
- **The guide** (start, manual in HTML and PDF, recipes, releases, roadmap, conventions):
  <https://jmigartua.github.io/diff-ref-struct-web/guide/>
- **Requests**: the [bug form](../../issues/new?template=bug.yml) or the [feature form](../../issues/new?template=feature.yml);
  from the app, **Help ▾** fills in the version and the view for you. Never attach unpublished data.

## What is here

The source the page is built from, published from the project's repository at each release (by its Action; this
repository is not edited by hand):

| folder | what |
|---|---|
| `core/` | the engine (Rust): readers, plot sets, maps, differences, windows; `core/wasm/` its WebAssembly build |
| `app/ui/` | the shell's services and the Plots module, shared with the diff-ref-struct Mac app; `app/ui/host-web/` the browser page |
| `web/` | the page's tools and tests, the guide's sources (`web/site/`, Quarto), the manuals (`web/manual/`) |

Build and run locally: `web/tools/build-wasm.sh` (Rust with the `wasm32-unknown-unknown` target and `wasm-pack`), then
`web/tools/serve.sh` and open <http://localhost:8765/host-web/>. Tests: `node --test web/tests/*.test.mjs`.

## Licence and credit

The code is under the **Apache License 2.0** (`LICENSE`); the guide, the manuals, the screenshots and the example data
under **CC BY 4.0** (`LICENSE-docs`). Attribution, third-party software and how the tool was made: `NOTICE`. To cite
it: `CITATION.cff` (GitHub's "Cite this repository"), with the version you used.

Josu M. Igartua, Departamento de Física, UPV/EHU · josu.igartua@ehu.eus
