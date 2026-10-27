#!/usr/bin/env python3
"""Write core/src/wyckoff_table.rs from spglib's database/Wyckoff.csv (https://github.com/spglib/spglib, BSD-3).

    python3 make_wyckoff_table.py <path to Wyckoff.csv> <spglib commit> > ../../src/wyckoff_table.rs

For each of the 530 settings (spglib's Hall numbers), each Wyckoff position: letter, multiplicity (conventional
cell), oriented site-symmetry symbol and the first coordinate triplet. The core uses the letter and the symbol as
names; multiplicities, site groups and orbits it computes itself, and the test `spglib_wyckoff` checks them."""
import re, sys

path, commit = sys.argv[1], sys.argv[2]
rows, index, h = [], {}, None
for line in open(path):
    line = line.rstrip("\n")
    m = re.match(r"^(\d+):", line)
    if m:
        h = int(m.group(1)); index[h] = [len(rows), 0]; continue
    f = line.split(":")
    if len(f) > 5 and f[2]:
        rows.append((f[3], int(f[2]), f[4], f[5].strip("()")))
        index[h][1] += 1
assert sorted(index) == list(range(1, 531))
print(f"""//! The Wyckoff positions of the 530 space-group settings: letter, multiplicity, site-symmetry symbol, coordinates.
//!
//! Source: spglib's `database/Wyckoff.csv` (https://github.com/spglib/spglib, commit {commit[:7]}, BSD-3 licence),
//! written by `core/tests/data/make_wyckoff_table.py`. Do not edit. The letters and site-symmetry symbols are used as
//! names; multiplicities, site groups and orbits are computed by `site`, and the test `spglib_wyckoff` checks them.

/// One Wyckoff position of a setting, as listed in the table (most special first is the table's order, `a` last).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct WyckoffEntry {{
    /// Wyckoff letter.
    pub letter: char,
    /// Multiplicity in the conventional cell (centring included).
    pub multiplicity: u16,
    /// Oriented site-symmetry symbol, e.g. ".m.", "2.22", "-1".
    pub site: &'static str,
    /// First coordinate triplet, e.g. "x,1/4,z", "x,2x,z".
    pub coord: &'static str,
}}

/// For Hall number h (1–530), `INDEX[h - 1]` = (first row in `ENTRIES`, number of positions).
pub static INDEX: [(u16, u8); 530] = [""")
for k in range(1, 531):
    print(f"    ({index[k][0]}, {index[k][1]}),")
print(f"];\n\npub static ENTRIES: [WyckoffEntry; {len(rows)}] = [")
for letter, mult, site, coord in rows:
    print(f'    WyckoffEntry {{ letter: \'{letter}\', multiplicity: {mult}, site: "{site}", coord: "{coord}" }},')
print("""];

/// The Wyckoff positions of a setting, in the table's order (general position first).
pub fn positions(hall_number: u16) -> &'static [WyckoffEntry] {
    let (start, n) = INDEX[hall_number as usize - 1];
    &ENTRIES[start as usize..start as usize + n as usize]
}""")
