#!/usr/bin/env python3
"""Write core/src/sg_fullprof.rs, the space-group symbols FullProf accepts, with their Hall symbols, from CrysFML
(the crystallographic library of FullProf, J. Rodríguez-Carvajal, ILL; LGPL-3): the `spgr_info_type` table of
Src/CFML_Symmetry_Tables.f90.

    python3 make_fullprof_symbols.py <path to CFML_Symmetry_Tables.f90> <crysfml commit> > ../../src/sg_fullprof.rs

Used to read the symbols written in FullProf .pcr files that the main table (`sgtable`, from spglib) names
differently, such as the older "B m a b" (now "B m e b") or "F m 3 m" (now "F m -3 m")."""
import re, sys

path, commit = sys.argv[1], sys.argv[2]
rows = []
for line in open(path, encoding="latin-1"):
    # Fortran strings: a quote inside is written twice ("P 3 2""" is the Hall symbol P 3 2").
    m = re.search(r'spgr_info_type\(\s*(\d+)\s*,\s*"((?:[^"]|"")*)"\s*,\s*"((?:[^"]|"")*)"', line)
    if m:
        rows.append((int(m.group(1)), m.group(2).replace('""', '"').strip(), m.group(3).replace('""', '"').strip()))
print(f"""//! The space-group symbols of CrysFML, the crystallographic library of FullProf: number, Hermann–Mauguin symbol
//! as FullProf writes it (upper case), Hall symbol.
//!
//! Source: CrysFML (J. Rodríguez-Carvajal, ILL; LGPL-3), `Src/CFML_Symmetry_Tables.f90`, table `Spgr_Info`, from the
//! copy at https://github.com/easyscience/crysfml, commit {commit[:7]}; written by
//! `core/tests/data/make_fullprof_symbols.py`. Do not edit. Used only to read FullProf's names of space groups;
//! the test `crysfml_symbols` checks every entry against the main table.

/// (number, Hermann–Mauguin symbol in CrysFML's upper-case form, Hall symbol)
pub static FULLPROF_SYMBOLS: [(u16, &str, &str); {len(rows)}] = [""")
for n, hm, hall in rows:
    print(f'    ({n}, "{hm}", "{hall.replace(chr(34), chr(92) + chr(34))}"),')
print("];")
