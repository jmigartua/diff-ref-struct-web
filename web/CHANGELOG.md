# Changelog of diff-ref-struct Web

*Generated from web/site/data/releases.yml by web/site/tools/build.py; edit the record, not this file.*

## 1.0.0 · The guide, open source, the record

*2026-10-27 · published · core 0.81.0 · build b4820d8*

The page gains its own site, a licence and a citation, a place for requests, and a record (releases, log, roadmap, conventions) that the site and the manuals render; versions become semantic.

**Added**

- The guide at /guide/: overview, a five-minute start with a check of what the browser supports, the manual and a one-page quick reference (HTML and PDF), five recipes, how it works, formats, conventions, releases, log, roadmap, requests, about.
- The manuals as PDF, typeset by LuaLaTeX in the site's typefaces; every edition with three appendices: Conventions, Release notes, Roadmap.
- The guide offline: after a first visit the guide's pages and the manual's PDFs open without network, as the page does.
- A Help menu in the page: the guide, the manual at the section of the view shown, what's new, Report a problem…, Ask for a feature…, Cite this view.
- Try the example: a synthetic heating run computed by the core (21 scans, 300 to 1300 K, a second phase appearing), opened as if a folder had been dropped.
- Open source: the code under the Apache License 2.0, the guide, manuals and screenshots under CC BY 4.0; CITATION.cff; issue forms for bugs and features.
- Semantic versions: diff-ref-struct Web 1.0.0, with the core's version and the build beside it.

**Fixed**

- The Overlay draws its curves: the browser's core named each pattern's member `index` where the page reads `member`; the page reads both (a rename, no number touched).
- Sum, Average and Difference of Members download their file (they were computed and lost).

**Known**

- Proposals against time open with each scan's temperature as its variable: open a hold or a run from a .drv that types the times (recipe: A hold).
- No figure or CSV export in the browser; text formats only.

*Closes W11; log L431, L433, L434; evidence: Web host tests 40 of 40 with the core 0.81.0 package; every guide screenshot taken from the app in both themes; the public PDFs and the guide checked for team data by the leak guard.; archived as [doi:10.5281/zenodo.23682727](https://doi.org/10.5281/zenodo.23682727).*

## 0.4.0 · The strips take the map's frame

*2026-10-26 · published · core 0.81.0 · build a1c308d*

The Plots module's step 23: the row strip under the map takes the map's frame, so that a peak sits under its column; the scale's word said once.

**Known**

- The Overlay draws no curves (the patterns' member named otherwise by the browser's core); the Stack and the maps are not affected.
- On the public page, in some window sizes, the strips may still sit off the map's frame.

*log L431, L432, L433.*

## 0.3.0 · Public at jmigartua.github.io/diff-ref-struct-web

*2026-10-25 · published · core 0.81.0 · build 24d8619*

The first publication by its Action, after the page's steps 4 and 5.

**Added**

- Live mode (Chromium): the folder read again every 2 to 60 s, new files described, sets opened from a proposal grown with the new scans.
- Offline: the page keeps its 36 files (4.4 MB) and opens without the network; installable as an app.
- The view in the address: the set, view, colour scale, zoom, clip and region kept in the link.
- Series recognised by the core: one scan per temperature, the holds, the run against time; exposure classes from the counts' scale σ²/y; large scans rebinned to 0.01° on opening.
- A hundred raw synchrotron scans opened and overlaid in seconds: the patterns asked when drawn, the overlay drawn in slices.

*log L427, L428, L429.*

## 0.2.0 · The core in the browser: files described, buffers and set files routed

*2026-10-25 · not published · core 0.79.0 · build bf451a3*

**Added**

- The core compiled to WebAssembly answers every command; each file described by the core (its kind in the tree, listed by folders or by kinds).
- A .drv or a .buf among the files opens as itself, its members taken from the folder.

*log L424, L425, L427.*

## 0.1.0 · The browser host

*2026-10-25 · not published · core — · build —*

**Added**

- The app's Plots module loaded as it is, in a page: the inspector, the folder's tree, the separator, the views.

*log L423.*
