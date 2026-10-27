#!/bin/sh
# Build the Python wheel of diff-ref-struct (P53, L380): the drs command line from the core, put inside the package,
# then the native module and the package by maturin. Usage, from core/python (maturin installed, e.g. `pip install
# maturin`):  ./build_wheel.sh  → target/wheels/diff_ref_struct-<version>-cp39-abi3-<platform>.whl
set -eu
here=$(cd "$(dirname "$0")" && pwd)
cargo build --release --manifest-path "$here/../Cargo.toml" --bin drs
mkdir -p "$here/python/drs/bin"
cp "$here/../target/release/drs" "$here/python/drs/bin/drs"
cd "$here"
maturin build --release "$@"
