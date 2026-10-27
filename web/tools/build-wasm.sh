#!/bin/sh
# diff-ref-struct Web · the core's WebAssembly package (W1) built into the browser host's folder and named in its
# manifest, so that the page loads it (app/ui/host-web/dispatch.js, loadPackage). `--stub` puts the stub's manifest back
# (no package: every command refused naming W1). The package itself is not committed (wasm/.gitignore), and the
# manifest a local build writes is kept out of `git status` (skip-worktree) until --stub puts the committed one back.
set -eu
top=$(cd "$(dirname "$0")/../.." && pwd)
out="$top/app/ui/host-web/wasm"
stub() {
  printf '{\n  "package": null,\n  "why": "the core'"'"'s WebAssembly build is request W1 (web/REQUESTS.md); web/tools/build-wasm.sh writes the package here and names it in this file"\n}\n' > "$out/manifest.json"
  printf '# wasm-pack'"'"'s output (W1), built by web/tools/build-wasm.sh: not committed; the manifest says whether it is here\n*\n!.gitignore\n!manifest.json\n' > "$out/.gitignore"
}
if [ "${1:-}" = "--stub" ]; then
  find "$out" -maxdepth 1 -type f ! -name manifest.json ! -name .gitignore -delete
  git -C "$top" update-index --no-skip-worktree app/ui/host-web/wasm/manifest.json 2>/dev/null || true
  stub; echo "build-wasm: the stub's manifest (no package)"; exit 0
fi
crate="$top/core/wasm"
[ -f "$crate/Cargo.toml" ] || { echo "build-wasm: no core/wasm/Cargo.toml yet: the crate is request W1"; exit 1; }
command -v wasm-pack >/dev/null 2>&1 || { echo "build-wasm: wasm-pack is not installed (cargo install wasm-pack; rustup target add wasm32-unknown-unknown)"; exit 1; }
wasm-pack build "$crate" --target web --release --out-dir "$out" --out-name drs_wasm
stub   # the host's .gitignore back over wasm-pack's (its "*" would hide the manifest)
bytes=$(wc -c < "$out/drs_wasm_bg.wasm" | tr -d ' ')
printf '{\n  "package": "drs_wasm.js",\n  "wasm": "drs_wasm_bg.wasm",\n  "bytes": %s,\n  "built": "%s"\n}\n' "$bytes" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$out/manifest.json"
git -C "$top" update-index --skip-worktree app/ui/host-web/wasm/manifest.json 2>/dev/null || true
echo "build-wasm: app/ui/host-web/wasm/drs_wasm.js and drs_wasm_bg.wasm ($((bytes / 1024)) kB); the manifest names them (--stub returns to the stub)"
