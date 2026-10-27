#!/usr/bin/env python3
"""diff-ref-struct Web · the release gate (the open-tool kit, web/PROPOSAL-open-tool.md §8): everything a version needs
before it goes out, checked at once, each line OK or FAIL with what to do. Exit 1 when any line fails.

  python3 web/site/tools/release.py            the checks that need no browser
  python3 web/site/tools/release.py --shots    also every view drawn (web/tools/serve.sh 8790 running; Chrome)

Checked:
  versions    the release record's version in CITATION.cff and in the manuals' first lines; its core is core/'s
  notes       the current release has a title, a summary and something added, changed or fixed
  pdfs        the public PDFs are current with their sources, and the manual's edition of this version is kept
  guide       the rendered guide (_site/) has no broken local link and no guarded word
  source      the public source's selection holds no guarded word (web/tools/source.mjs)
  licence     LICENSE, LICENSE-docs, NOTICE, CITATION.cff and the issue forms are in web/public/
  tests       the page's tests pass (node --test web/tests/*.test.mjs)
  shots       (--shots) every view of the example set drawn, none blank
"""
import os
import re
import subprocess
import sys
import urllib.parse
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
SITE = HERE.parent
WEB = SITE.parent
TOP = WEB.parent
results = []


def line(name, ok, say):
    results.append((name, ok, say))
    print(f"{'OK  ' if ok else 'FAIL'} {name:9} {say}")


def run(cmd, cwd=TOP):
    p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    return p.returncode, (p.stdout + p.stderr).strip()


rel = yaml.safe_load((SITE / "data/releases.yml").read_text(encoding="utf-8"))
cur = rel[0]
v = cur["version"]

# versions
cff = (WEB / "public/CITATION.cff").read_text(encoding="utf-8")
cff_v = re.search(r'^version:\s*"?([^"\n]+)"?', cff, re.M)
core_v = re.search(r'^version\s*=\s*"([^"]+)"', (TOP / "core/Cargo.toml").read_text(encoding="utf-8"), re.M)
man = [(f, re.search(r"diff-ref-struct Web ([0-9][0-9.]*) with core", (WEB / "manual" / f).read_text(encoding="utf-8")[:1500])) for f in ("user-manual.md", "quick-reference.md")]
bad = []
if not cff_v or cff_v.group(1) != v:
    bad.append(f"CITATION.cff says {cff_v.group(1) if cff_v else '?'}")
for f, m in man:
    if not m or m.group(1) != v:
        bad.append(f"{f} says {m.group(1) if m else '?'}")
if core_v and core_v.group(1) != cur["core"]:
    bad.append(f"the record's core {cur['core']}, core/Cargo.toml {core_v.group(1)}")
line("versions", not bad, f"{v} everywhere" if not bad else "; ".join(bad))

# notes
has = any(cur.get(k) for k in ("added", "changed", "fixed"))
line("notes", bool(cur.get("title") and cur.get("summary") and has), f"{v}: {cur.get('title', '(no title)')}" if has else f"{v} has no added, changed or fixed")

# pdfs
code, out = run([sys.executable, str(HERE / "build.py"), "--check"], cwd=SITE)
kept = (SITE / f"downloads/user-manual-{v}.pdf").exists()
line("pdfs", code == 0 and kept, "current, the edition of this version kept" if code == 0 and kept else (out.splitlines()[-1] if code else f"downloads/user-manual-{v}.pdf missing: web/site/tools/pdf.sh"))

# guide: links and guard
site = SITE / "_site"
if not site.exists():
    line("guide", False, "not rendered: cd web/site && quarto render")
else:
    broken = []
    for p in site.rglob("*.html"):
        if "assets/includes" in p.as_posix():
            continue
        for m in re.finditer(r'(?:href|src)="([^"#?]+)', p.read_text(encoding="utf-8", errors="replace")):
            u = m.group(1)
            if re.match(r"^[a-z]+:|^//|^data:", u):
                continue
            q = (p.parent / urllib.parse.unquote(u)).resolve()
            if q.is_dir():
                q = q / "index.html"
            if not q.exists():
                broken.append(f"{p.relative_to(site)} → {u}")
    words = [w.strip().lower() for w in (WEB / "private/guard.txt").read_text(encoding="utf-8").splitlines() if w.strip() and not w.startswith("#")] if (WEB / "private/guard.txt").exists() else []
    leaks = [p.relative_to(site).as_posix() for p in site.rglob("*") if p.is_file() and p.suffix in (".html", ".xml", ".json", ".md") and any(w in p.read_text(encoding="utf-8", errors="replace").lower() for w in words)]
    line("guide", not broken and not leaks, "no broken link, no guarded word" if not broken and not leaks else "; ".join((broken + [f"guarded word in {x}" for x in leaks])[:6]))

# source
code, out = run(["node", "web/tools/source.mjs", "--audit"])
hits = [l.strip() for l in out.splitlines() if l.startswith("  ")]
line("source", code == 0, "no guarded word in the public selection" if code == 0 else f"{len(hits)} files hold a guarded word: " + "; ".join(hits[:4]))

# licence
need = ["LICENSE", "LICENSE-docs", "NOTICE", "CITATION.cff", "README.md", "CONTRIBUTING.md", ".github/ISSUE_TEMPLATE/bug.yml", ".github/ISSUE_TEMPLATE/feature.yml"]
miss = [f for f in need if not (WEB / "public" / f).exists()]
line("licence", not miss, "licences, NOTICE, citation and the issue forms present" if not miss else "missing: " + ", ".join(miss))

# tests
tests = sorted(str(p.relative_to(TOP)) for p in (WEB / "tests").glob("*.test.mjs"))
code, out = run(["node", "--test", *tests])
m = re.search(r"# pass (\d+)\n# fail (\d+)", out)
line("tests", code == 0, f"{m.group(1)} pass, {m.group(2)} fail" if m else out.splitlines()[-1] if out else "no output")

# shots
if "--shots" in sys.argv:
    code, out = run(["node", "--experimental-websocket", "web/site/tools/shots.mjs", "--port", "8790", "--check"])
    last = [l for l in out.splitlines() if l.startswith("check:")]
    line("shots", code == 0, last[-1][7:] if last else out.splitlines()[-1] if out else "no output")

failed = [n for n, ok, _ in results if not ok]
print(f"\nrelease {v}: " + ("ready" if not failed else f"not ready ({', '.join(failed)})"))
sys.exit(1 if failed else 0)
