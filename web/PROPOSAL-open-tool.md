# The open-tool kit: what a public tool of the family carries

A proposal, 26 October 2026, from the web stream. **Decided the same day by Josu: every recommendation, English only for now; applied in diff-ref-struct Web 1.0.0 (web/README.md, step 7).** It starts from diff-ref-struct Web's guide (step 6)
and generalises it into a standard that every tool we make public follows. The standard would become a skill working
beside the family design kit (the Quarto site skill) and `_SIBLING`.

**In one line:** a public tool is four things kept in step: **the app**, **its guide**, **its source**, and **its record**.
The record holds the releases, the log, the roadmap, the requests and the conventions. Each of these has one data file
as its source, and the guide, the manuals, the app and the repository all render that same file. A release is the
moment they are all regenerated together and frozen.

---

## 1. What I think

Making the web app public and open source changes three things that the guide does not yet answer:

- **who may reuse the work, and how to credit it**: a licence, a citation and a DOI;
- **where the source lives**: today the public repository holds only the built page, and the source sits in a private
  monorepo, next to things that must not go public;
- **how users talk back**: a place for bugs and feature requests, and a path from a request to a release note.

Your three asks fit here: a log, releases that also appear in the manuals whatever the version, and a vocabulary
appendix. Each one is a view of the record. Building them once, as data with generators, is what makes them reusable
by the next tool.

Today showed why this matters. The screenshots for the guide found a bug on the public page: the Overlay drew no
curves. A record that regenerates the screenshots at every release is also a test. That bug would have stopped the
release.

## 2. The shape of a public tool

| part | where | what it holds |
|---|---|---|
| **App** | `<tool>/` (the root of the address) | the tool itself; its version visible; links to the guide; *Report a problem* and *Ask for a feature* in its menu |
| **Guide** | `<tool>/guide/` | the site in the family kit (below) |
| **Source** | a public repository | the code the app is built from, with its licence, citation and contribution notes |
| **Record** | `data/*.yml` in the source | releases, log, roadmap, requests, conventions: the single sources everything renders |

The guide's sections. Those marked **new** are not in the guide yet:

1. **Overview**: the home.
2. **Start**: a browser check, then a public example from opening to a first result.
3. **Manual**: HTML and PDF, versioned, with three appendices: **Conventions**, **Release notes** and **Roadmap**.
4. **Recipes**: the analyses the tool is made for.
5. **How it works**: architecture, privacy, offline use, how the numbers are checked.
6. **Releases** (new): every version, newest first, with what it added, changed, fixed and still lacks. Each entry
   links to its manual edition and its log entries. An Atom feed.
7. **Log** (new): the dated record of the work on this tool: findings, decisions with their reasons, work done.
   Filtered from the project's log, with a feed.
8. **Roadmap** (new): decided, under way, considered, set aside, each with its reason.
9. **Requests** (new): the place for users. How to report a bug or ask for a feature, what is open, what was done
   and in which version.
10. **Conventions** (new): the vocabulary, the files and their suffixes, the naming of scans, units and signs,
    versions and ids. One page, and an appendix of every manual.
11. **About**: credits, licence, how to cite, privacy, accessibility.

## 3. One source per record, many renderings

This is the heart of the kit: each record is written once, as data, and rendered everywhere it is needed.

| source | rendered as |
|---|---|
| `releases.yml` | the Releases page; **Appendix B, Release notes**, in every manual edition (all releases up to that edition's version, so a PDF from version 1.2 says what 1.2 is); `CHANGELOG.md`; an Atom feed; the app's *What's new* on first start after an update; the GitHub release; the version in `CITATION.cff` |
| `log.yml` (the project's, entries tagged `tools: [web]`) | the Log page and feed; each release links the entries it closes |
| `roadmap.yml` (from requests, ideas and decisions) | the Roadmap page; **Appendix C** of the manual; the Versions page's *Coming next* |
| GitHub Issues (bug and feature forms) | the Requests page, built from the issues at render time; a request closed by a release is named in that release's notes |
| `conventions.yml` (with the file dictionary `ecosystem/data/files.yml`, which already generates the core's `file_table.rs`) | the Conventions page; **Appendix A** of every manual; the app's hover text for file kinds |

**Versions.** The tool takes a semantic version of its own (diff-ref-struct Web 1.0.0), shown beside the core's
(core 0.81.0) and the build (a1c308d). The manual's cover, its running head and its release-notes appendix carry that
version. Every release keeps its manual edition at `guide/v/<version>/manual.pdf`. The app's ? opens the manual of its
own version, falling back to the latest. That answers "whatever the version is": a user of an older copy reads the
manual that matches it.

## 4. Conventions: the name for the vocabulary appendix

I propose **Conventions**, as Appendix A of each manual and a page of the guide. "Glossary" covers only words. This
appendix also covers how files are named, what their suffixes mean, how scans should be named so that the core can read
them, and what the ids and versions mean. Alternatives set aside: *Lexicon* (words only), *Field guide* (it suits
recipes better).

Its sections, with this tool's entries:

| section | entries (diff-ref-struct Web) |
|---|---|
| **Words** | plot set, member, variable and its source, second variable, proposal, exposure class, rebin, clip, region, window, reference member, live set |
| **Files and suffixes** | `.drv` plot set, `.buf` WinPLOTR buffer, `.drp` publication set, `.drc` control file, `.xye` / `.xys` / `.qye` data columns, `-rebin0p01.xys` (a scan rebinned to 0.01°), `Plots/` `Data/` `Rebinned/` `Publish/` (the project's folders), `.drs/` (the Mac app's data operations, hidden) |
| **Naming scans** | `T0900C` (900 °C), `_300K` / `T300` (K), `t16530s` (seconds since the run's first scan), `_long` / `_short` (exposure), a three-digit sequence number first so that names sort in acquisition order |
| **Units and signs** | 2θ in degrees, d in Å, Q in Å⁻¹, √counts, log₁₀ counts; the decimal point (the decimal comma in Basque editions) |
| **Versions and ids** | web 1.2.0 · core 0.81.0 · build a1c308d; L… log entries, I… ideas, W… web requests, #… public issues |

The files section is generated from the file dictionary, so the manual, the core and the guide never disagree on what
a suffix means.

## 5. Requests: the place for users

**Recommended:** the public source repository's **GitHub Issues**, with two issue forms, plus a path for users without
a GitHub account.

- **Bug report form**: what happened, what was expected, and the version and browser, prefilled when the report starts
  from the app. It asks for the page's address at that moment, which carries the view and no data. It warns: *do not
  attach unpublished data; a description and the address are enough*.
- **Feature request form**: what you want to do, with which data, and what you do today instead.
- **Discussions** for questions and showing work. Not needed at first.
- **No account**: *Report a problem* in the app opens either a GitHub issue with the fields filled in, or an e-mail
  with the same fields, so that nobody is turned away.
- **The flow**: issue → labels (`bug`, `feature`, `question`; `confirmed`, `planned`, `done`) → a feature becomes an
  idea in the seed-place (I…, as the project's rule says) and enters the roadmap only on your decision → the release
  that does it names the issue, and the issue is closed with the version.
- **The Requests page** lists the open requests by label and the ones closed in each version. It is built from the
  issues at every render, so the guide shows what users asked for and what came of it.

Set aside:

- a form on a third-party service, because it sends users' words to a server we do not control and breaks the
  no-upload promise in spirit;
- e-mail only, because nothing is public and the next user cannot see that the bug is known.

## 6. Open source and attribution

Recommendations (each a decision for you):

- **Licence for the code: Apache-2.0** (or MIT). Both are permissive. Apache-2.0 adds an explicit patent grant and a
  NOTICE file that keeps attribution with every copy, which is the "proper attribution" you asked for. GPL-3.0 if you
  want derived works to stay open; it would make reuse inside closed tools impossible.
- **Licence for the guide, the manuals and the screenshots: CC BY 4.0.** Reuse is free with credit.
- **Citation**: a `CITATION.cff` in the repository, and a **Zenodo DOI per release** (Zenodo archives each GitHub
  release automatically). The guide's *How to cite* gives the DOI of the version the user ran.
- **Credits**: authors; third-party libraries with their licences (already in `LICENCES.md`); example data with its
  source; the typefaces (SIL OFL); and, your decision, a line saying the code and the documents were written with AI
  assistance under your direction.
- **What must not go public**:
  - FullProf's and GSAS-II's manuals and papers in `manuals/` (their publishers' copyright);
  - the ESRF team's data and manual;
  - anything a leak guard catches (a per-project list of private words, checked on every render and publication, as
    the guide's Action already does).
- **The public source repository.** I recommend that an Action publishes the subset the app is built from (`core/`,
  `app/ui/` with the shared modules, `web/`) into the public repository, with each release tagged. The monorepo stays
  private. Set aside:
  - making the whole monorepo public: it would need an audit of 400+ log entries and of `manuals/`;
  - copying the code by hand: the copy would drift from the source.

## 7. Your four ideas, as I would do them

1. **The guide available offline.** The service worker already keeps the app (36 files, 4.4 MB). It would also keep
   the guide's pages, the manual PDF and the quick reference, about 3 MB more without the screenshots. The manual then
   opens at a beamline without network. The screenshots would be fetched on demand and kept once seen. Small; no
   decision needed.
2. **Each panel's ? opening its own section.** Every manual section gets a stable id (`{#map-clip}`, `{#live-mode}`),
   listed in the conventions so that renaming a heading never breaks a link. The module's ? passes the panel and view
   on screen, and the host opens `manual/#<id>`. Needs a small hook in the Plots module (a request to its stream).
   The ids also make the manual's cross-references (§11) into links.
3. **Spanish and Basque editions.** Manuals and guide first, the app's interface later (it has no translation layer
   yet).
   - Basque goes through the EHULKU skill's two stages (calques flagged against the EHU rules), keeping the decimal
     comma.
   - Spanish takes the decimal point, as your rule L284 says.
   - Each edition is a folder (`guide/es/`, `guide/eu/`), each manual a PDF per language.
   - The conventions stay one file, with a column per language, so that terms are translated once.
   - Translations are regenerated when the English changes, and the release notes say which editions lag.
4. **A "try the example" button.** I would no longer use FullProf's tutorial files for this. I propose an example the
   core **computes**:
   - `drs simulate` already makes a pattern with counting noise from a control file;
   - a public structure, for example corundum or silicon from the COD (public domain), with a cell that grows with
     temperature and a second phase that appears above a set temperature, gives a ramp of 20 synthetic scans;
   - the series is written once by a script in the repository, versioned, and shipped with the page, about 100 kB.
   It has no licence question, and it shows the core at work. The button opens it as if a folder had been dropped, with
   the Start page's tour on top of it. The FullProf files stay as the second example, linked and not hosted. Hosting
   them would need the ILL's permission, and asking for it is your decision.

## 8. Further proposals

- **Screenshots as a release check.** `shots.mjs` runs at every release, in both themes. A view that comes out blank or
  changed beyond a threshold stops the release, as today's blank Overlay would have.
- **A Known issues box** on the home page and in the app's start page, from `releases.yml`'s `known` of the current
  version, so that users read a fault before they meet it.
- **The evidence in each release note**: the tests that passed (counts), the core version, the browsers checked.
  Consistent with the project's rule that a claim carries its evidence.
- **"Cite this view"**: the app's address plus the version and a date, copied as one line for a lab book or a paper's
  methods.
- **Privacy and accessibility statements**: short, true, and the same shape for every tool. The privacy statement
  lists the only things kept: local storage keys and the cache, never data.
- **A release gate**: one command that checks everything before a version goes out. It fails with a list of what
  is missing:
  - the version bumped in every place;
  - release notes written;
  - manual edition built;
  - PDFs current;
  - screenshots fresh;
  - links unbroken;
  - leak guard clean;
  - licence and citation present;
  - tests green.

## 9. The skill

Proposed name **`_OPENTOOL`** (the open-tool kit). It works with the Quarto skill (the family design kit) for the
pages and with `_SIBLING` for editions. Created through CreateSkill when you say so. What it does:

- **scaffold**: for a new public tool, the guide's skeleton (the sections of §2), the record's data files with their
  generators, the PDF pipeline (LuaLaTeX in the kit's typefaces), the screenshot tool, the leak guard, the issue forms,
  `CITATION.cff`, `LICENSE`, `NOTICE`, the Action;
- **release**: the gate of §8, then regenerating and freezing everything for a version, tagging it and asking Zenodo
  for the DOI;
- **check**: the gate alone, at any time;
- **triage**: read the open issues, propose labels, and file features as seed-place ideas for your decision.

## 10. The order I propose for diff-ref-struct Web

| phase | what | needs |
|---|---|---|
| A, now | Releases page and manual appendix from `releases.yml`; Log page from the project's log tagged `web`; Conventions page and appendix; Roadmap page and appendix; the guide offline; screenshots as a release check | nothing: in the web stream's folder |
| B | licence, `CITATION.cff`, the public source repository by Action, issue forms, Requests page, *Report a problem* in the app, semantic versions starting at 1.0.0, Zenodo | your decisions in §6 and §5 |
| C | per-panel ?, the computed example and its button, Spanish and Basque editions | a request to the Plots stream; the example script with the core; EHULKU |
| D | `_OPENTOOL` as a skill, written from what A–C settled | CreateSkill |

## 11. Decisions for you

1. **Licence**: Apache-2.0 for the code and CC BY 4.0 for documents (recommended), MIT, or GPL-3.0.
2. **Source**: a public subset published by an Action (recommended), or the whole monorepo after an audit.
3. **Requests**: GitHub Issues with forms and an e-mail fallback (recommended).
4. **The appendix's name**: *Conventions* (recommended).
5. **Versions**: start semantic versions at diff-ref-struct Web 1.0.0 with the next release, or keep the build hash.
6. **The example**: computed by the core (recommended), or FullProf's tutorial hosted with the ILL's permission.
7. **AI-assistance statement** in the credits: yes or no, and its wording.
