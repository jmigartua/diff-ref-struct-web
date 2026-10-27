#!/usr/bin/env python3
"""diff-ref-struct Web · the guide's offline list (Quarto's post-render): every page, style, script, feed and the
current manuals' PDFs of the rendered guide, written to _site/offline.json, which the page's service worker reads when it
installs (app/ui/host-web/sw.js) so that the guide opens without network. The screenshots are left out (they are kept
once seen) and so are the PDFs of earlier versions."""
import json
from pathlib import Path

SITE = Path(__file__).resolve().parent.parent / "_site"
KEEP = {".html", ".css", ".js", ".svg", ".xml", ".json"}
files = []
for p in sorted(SITE.rglob("*")):
    if not p.is_file():
        continue
    rel = p.relative_to(SITE).as_posix()
    if rel == "offline.json" or rel.startswith("assets/figures/") or rel.startswith("assets/includes/"):
        continue
    if p.suffix in KEEP or rel in ("downloads/user-manual.pdf", "downloads/quick-reference.pdf"):
        files.append(rel)
(SITE / "offline.json").write_text(json.dumps({"note": "The guide's files the page's service worker keeps for offline use (tools/offline.py).", "files": files}, indent=1) + "\n", encoding="utf-8")
print(f"offline: {len(files)} files listed for the service worker")
