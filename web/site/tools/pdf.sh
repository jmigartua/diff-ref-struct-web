#!/bin/sh
# diff-ref-struct Web · the manuals as PDF, typeset by LuaLaTeX through pandoc (Quarto's), from the same Markdown the site
# shows (tools/build.py --source: the H1 and the hand-made contents dropped, a team's parts left out of the public
# edition). Writes:
#   downloads/user-manual.pdf, downloads/quick-reference.pdf      public, published with the site; their sources'
#                                                                  hashes in downloads/pdfs.json (build.py checks them)
#   ../manual/team/<name>.pdf and .html                            each team's edition (unpublished data): never in the
#                                                                  site, never staged
# The typefaces are the site's (Source Sans 3, Source Serif 4, IBM Plex Mono; SIL OFL), fetched once into .fonts/
# (git-ignored) from their makers' releases. Run from anywhere: web/site/tools/pdf.sh [user-manual|quick-reference|team]
set -eu
here=$(cd "$(dirname "$0")" && pwd)
site=$(cd "$here/.." && pwd)
manual=$(cd "$site/../manual" && pwd)
fonts="$site/.fonts"
latex="$site/latex"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
only="${1:-all}"

fetch_fonts() {
  [ -f "$fonts/SourceSerif4-Regular.otf" ] && [ -f "$fonts/SourceSans3-Regular.otf" ] && [ -f "$fonts/IBMPlexMono-Regular.otf" ] && return
  mkdir -p "$fonts/zip"
  echo "pdf: fetching the site's typefaces (SIL OFL) into $fonts"
  curl -sSL -o "$fonts/zip/serif.zip" https://github.com/adobe-fonts/source-serif/releases/download/4.005R/source-serif-4.005_Desktop.zip
  curl -sSL -o "$fonts/zip/sans.zip" https://github.com/adobe-fonts/source-sans/releases/download/3.052R/OTF-source-sans-3.052R.zip
  curl -sSL -o "$fonts/zip/mono.zip" "https://github.com/IBM/plex/releases/download/%40ibm%2Fplex-mono%401.1.0/ibm-plex-mono.zip"
  for w in Regular It Semibold SemiboldIt Bold BoldIt; do
    unzip -jo -q "$fonts/zip/serif.zip" "source-serif-4.005_Desktop/OTF/SourceSerif4-$w.otf" -d "$fonts"
    unzip -jo -q "$fonts/zip/sans.zip" "OTF/SourceSans3-$w.otf" -d "$fonts"
  done
  for w in Regular Italic Medium; do unzip -jo -q "$fonts/zip/mono.zip" "ibm-plex-mono/fonts/complete/otf/IBMPlexMono-$w.otf" -d "$fonts"; done
}

# pdf <name> <edition> <out.pdf> <title> <short title> <subtitle> <audience> [quick]
pdf() {
  name=$1; edition=$2; out=$3; title=$4; short=$5; sub=$6; audience=$7; kind=${8:-manual}
  python3 "$here/build.py" --source "$name" --edition "$edition" > "$tmp/$name.md"
  version="diff-ref-struct Web $(python3 "$here/build.py" --version) with core $(python3 -c "import yaml;print(yaml.safe_load(open('$site/data/releases.yml'))[0]['core'])")"
  [ -n "$version" ] || version="web, core 0.81.0"
  {
    printf '\\newcommand{\\fontpath}{%s}\n' "$fonts"
    printf '\\newcommand{\\doctitle}{%s}\n\\newcommand{\\doctitleshort}{%s}\n' "$title" "$short"
    printf '\\newcommand{\\docsub}{%s}\n\\newcommand{\\docaudience}{%s}\n' "$sub" "$audience"
    printf '\\newcommand{\\docversion}{For %s · %s}\n' "$version" "$(date +%Y-%m-%d)"
    cat "$latex/manual.tex"
    [ "$kind" = quick ] && cat "$latex/quick.tex"
  } > "$tmp/head.tex"
  if [ "$kind" = quick ]; then
    printf '{\\drshead\\Large diff-ref-struct \\textcolor{drsgreen}{Web} · quick reference}\\hfill{\\small\\color{drsfaint}\\docversion\\ · jmigartua.github.io/diff-ref-struct-web/guide}\\par\\vspace{2pt}{\\color{drsline}\\rule{\\linewidth}{0.4pt}}\\par\\vspace{4pt}\n\\begin{multicols}{3}\n' > "$tmp/before.tex"
    printf '\\end{multicols}\n' > "$tmp/after.tex"
    set -- -V "geometry=a4paper,landscape,margin=1.0cm" -V classoption=fontsize=10.5pt -B "$tmp/before.tex" -A "$tmp/after.tex"
  else
    set -- -V "classoption=a4paper,11pt,DIV=12" -B "$latex/title.tex"
  fi
  quarto pandoc "$tmp/$name.md" -f markdown -o "$out" --pdf-engine=lualatex -V documentclass=scrartcl \
    --shift-heading-level-by=-1 --lua-filter "$latex/tablewidths.lua" --resource-path "$manual" \
    -H "$tmp/head.tex" "$@" -M title= -M lang=en-GB
  echo "pdf: $out ($(du -k "$out" | cut -f1) kB)"
}

record() {   # the source's hash beside the public PDF, which build.py compares with the manual as it is now
  python3 - "$site/downloads/pdfs.json" "$1" "$(python3 "$here/build.py" --source "$1" --edition public | shasum -a 256 | cut -d' ' -f1)" <<'EOF'
import json, sys, datetime
path, name, sha = sys.argv[1:4]
try: rec = json.load(open(path))
except FileNotFoundError: rec = {}
rec[name] = {"sha256": sha, "built": datetime.date.today().isoformat()}
json.dump(rec, open(path, "w"), indent=2, sort_keys=True); open(path, "a").write("\n")
EOF
}

fetch_fonts
mkdir -p "$site/downloads" "$manual/team"
if [ "$only" = all ] || [ "$only" = user-manual ]; then
  pdf user-manual public "$site/downloads/user-manual.pdf" "User manual" "User manual" \
    "Powder diffraction series in the browser: every control, what the core computes behind it, the limits of this version" \
    "For anyone who opens a folder of diffraction patterns in the page."
  record user-manual
  # the edition of this version kept beside the current one (the app's Help ▾ ▸ Manual of this version opens it)
  cp "$site/downloads/user-manual.pdf" "$site/downloads/user-manual-$(python3 "$here/build.py" --version).pdf"
fi
if [ "$only" = all ] || [ "$only" = quick-reference ]; then
  pdf quick-reference public "$site/downloads/quick-reference.pdf" "Quick reference" "Quick reference" "" "" quick
  record quick-reference
fi
if [ "$only" = all ] || [ "$only" = team ]; then
  # every team's edition (a manual whose first line says <!-- edition: team; title: …; short: …; sub: …; audience: … -->),
  # as PDF and as one self-contained HTML page in the site's look, into ../manual/team/ (never staged, never published)
  for src in "$manual"/*.md; do
    head -c 2000 "$src" | grep -q '<!-- edition: team' || continue
    name=$(basename "$src" .md)
    field() { python3 -c "import re,sys; m=re.search(r'<!--\s*edition:\s*team;(.*?)-->', open(sys.argv[1], encoding='utf-8').read(2000), re.S); d=dict(p.split(':',1) for p in [x.strip() for x in m.group(1).split(';')] if ':' in p); print(d.get(sys.argv[2],'').strip())" "$src" "$1"; }
    title=$(field title); short=$(field short); sub=$(field sub); audience=$(field audience)
    pdf "$name" team "$manual/team/$name.pdf" "$title" "$short" "$sub" "$audience"
    python3 "$here/build.py" --source "$name" --edition team > "$tmp/team.md"
    plain=$(printf '%s' "$title" | sed 's/\\//g')
    printf '<link href="https://fonts.googleapis.com/css2?family=Source+Serif+4:opsz,wght@8..60,400;8..60,600&family=Source+Sans+3:ital,wght@0,400;0,600;1,400&family=IBM+Plex+Mono:wght@400&display=swap" rel="stylesheet" />\n<style>body{background:var(--paper);color:var(--ink);font-family:var(--sans);margin:0}.page-article{padding-top:28px;padding-bottom:64px}</style>\n' > "$tmp/head.html"
    printf '<div class="topline"></div><article class="page-article"><p class="sub" style="font-family:var(--mono);font-size:12px;letter-spacing:.06em;text-transform:uppercase;color:var(--teal)">diff-ref-struct Web · a worked use case</p><h1>%s</h1>\n' "$plain" > "$tmp/before.html"
    printf '</article>\n' > "$tmp/after.html"
    quarto pandoc "$tmp/team.md" -f markdown -t html5 -s --embed-resources -o "$manual/team/$name.html" \
      --css "$site/assets/styles/theme.css" --css "$site/assets/styles/site.css" --css "$site/assets/styles/project.css" \
      -H "$tmp/head.html" -B "$tmp/before.html" -A "$tmp/after.html" --toc --toc-depth=2 \
      -M title="$plain · diff-ref-struct Web" -M pagetitle="$plain · diff-ref-struct Web" -V title= 2>/dev/null
    echo "pdf: $manual/team/$name.html"
  done
fi
