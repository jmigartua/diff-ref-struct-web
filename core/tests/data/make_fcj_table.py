#!/usr/bin/env python3
"""Write core/src/fcj_table.rs, the Gauss-Legendre abscissas and weights that CrysFML's Finger-Cox-Jephcoat axial
divergence uses (Src/CFML_PowderProfiles_Finger.f90, subroutine Init_Prof_Val; CrysFML, J. Rodriguez-Carvajal, ILL,
LGPL-3), copied as written.

    python3 make_fcj_table.py <path to CFML_PowderProfiles_Finger.f90> <crysfml commit> > ../../src/fcj_table.rs
"""
import re, sys

path, commit = sys.argv[1], sys.argv[2]
text = open(path, encoding="latin-1").read()
start = text.index("Subroutine Init_Prof_Val()")
end = text.index("End Subroutine Init_Prof_Val")
body = re.sub(r"&\s*\n\s*", " ", text[start:end])
arrays = {"xp": [0.0] * 1884, "wp": [0.0] * 1884}
n = 0
for name, a, b, vals in re.findall(r"\b(xp|wp)\((\d+):(\d+)\)\s*=\s*\(/(.*?)/\)", body, re.S):
    v = [float(x.replace("D", "E")) for x in vals.replace("\n", " ").split(",") if x.strip()]
    a, b = int(a), int(b)
    assert len(v) == b - a + 1, (name, a, b, len(v))
    arrays[name][a:b + 1] = v
    n += len(v)
assert n == 2 * 1883, n
print(f"""//! Gauss–Legendre abscissas and weights of CrysFML's Finger–Cox–Jephcoat axial-divergence profile: the positive
//! halves of the rules with 6, 10, 20, 40, 60, 80, 100, 150, 200, 300, 400, 600, 800 and 1000 points, stored one after
//! another from index 1 (index 0 unused).
//!
//! Source: CrysFML (J. Rodríguez-Carvajal, ILL; LGPL-3), `Src/CFML_PowderProfiles_Finger.f90`, subroutine
//! `Init_Prof_Val`, from the copy at https://github.com/easyscience/crysfml, commit {commit[:7]}; written by
//! `core/tests/data/make_fcj_table.py`. Do not edit.
#![allow(clippy::approx_constant, clippy::excessive_precision)]

/// Abscissas (index 0 unused).
pub static XP: [f64; 1884] = [""")
for x in arrays["xp"]:
    print(f"    {x!r},")
print("];\n\n/// Weights (index 0 unused).\npub static WP: [f64; 1884] = [")
for x in arrays["wp"]:
    print(f"    {x!r},")
print("];")
