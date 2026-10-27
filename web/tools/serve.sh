#!/bin/sh
# diff-ref-struct Web · the browser host served locally as a static folder (step 2): app/ui/ served, the page at
# /host-web/. The shell's vendored libraries (dockview, three) are generated into app/ui/vendor/ by app/tools/vendor.js
# (npm ci in app/ first); without them the page cannot load the dock.
set -eu
top=$(cd "$(dirname "$0")/../.." && pwd)
port="${1:-8765}"
[ -f "$top/app/ui/vendor/dockview/dockview-core.esm.js" ] || { echo "serve: app/ui/vendor/ is missing: cd app && npm ci && node tools/vendor.js"; exit 1; }
echo "serve: http://localhost:$port/host-web/"
cd "$top/app/ui" && exec python3 -m http.server "$port" --bind 127.0.0.1
