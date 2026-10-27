#!/usr/bin/env python3
"""Write fullprof_reflections.txt, the reference for the test `fullprof_reflections`, from the FullProf .out files
in manuals/fullprof/ (standard library only; run from this folder: python3 make_fullprof_reflections.py > fullprof_reflections.txt).

FullProf prints, after the last cycle, the reflection list of each phase: "No. Code H K L Mult ... 2theta ... d-hkl"
(two layouts, read by the header's column names). For each such table one PHASE line gives the Hall symbol of the
phase, the cell the list was computed with and the first wavelength ("=> Wavelengths:"),
then one REFL line per reflection, in FullProf's order, with file and line:

    PHASE "<Hall symbol>" <a> <b> <c> <alpha> <beta> <gamma> <lambda1> <file>
    REFL <h> <k> <l> <mult> <2theta> <d> <file>:<line>

Every value is copied as printed, except the cell, which is computed from printed values. FullProf prints a
"=> Cell parameters:" block after each cycle: per parameter the new value, the shift just applied and the standard
deviation. The list is computed in the last cycle with the cell before that cycle's shifts, i.e. value − shift of
the last block whose shifts are not all zero (a block with zero shifts reports no new cycle). Two cases decided it:
LBF_lmopr.out, last block a = 5.748007 with shift 0.000003, and d(200) = 2.874002 = (5.748007 − 0.000003)/2 in the list;
LaMnO3_cub.out, last blocks b = 7.694143 with shifts 0.000014 then 0, and d(020) = 3.847064 = (7.694143 − 0.000014)/2.
"""
import os, re, sys

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../manuals/fullprof"))
REPO = os.path.dirname(os.path.dirname(ROOT))


def main():
    files = sorted(os.path.join(d, f) for d, _, fs in os.walk(ROOT) for f in fs if f.lower().endswith(".out"))
    print("# fullprof_reflections.txt: FullProf's reflection lists, the reference for the core's reflection tests.")
    print("# Written by make_fullprof_reflections.py (core/tests/data) from manuals/fullprof/**/*.out; do not edit.")
    print("# PHASE \"<Hall symbol>\" a b c alpha beta gamma lambda1 <file>:<line of the cell block>;  REFL h k l mult 2theta d <file>:<line>")
    n_phases = n_refl = 0
    for path in files:
        rel = os.path.relpath(path, REPO)
        lines = open(path, encoding="latin-1").read().splitlines()
        heads = [i for i, l in enumerate(lines) if re.search(r"No\.\s+Code\s+H\s+K\s+L\s+Mult", l)]
        if not heads:
            continue
        hall = next(l.split(":", 1)[1].strip() for l in lines if "Hall Symbol:" in l)
        lam = next(l.split(":", 1)[1].split()[0] for l in lines if "=> Wavelengths:" in l)
        # the cell of the last cycle: value − shift in the last "=> Cell parameters:" block with a non-zero shift
        blocks = [i for i, l in enumerate(lines[: heads[0]]) if l.strip().startswith("=> Cell parameters:")]
        vals = {b: [lines[b + 1 + j].split()[:2] for j in range(6)] for b in blocks}
        moved = [b for b in blocks if any(float(v[1]) != 0.0 for v in vals[b])]
        cblock = moved[-1] if moved else blocks[-1]
        cell = [f"{float(v) - float(sh):.6f}" for v, sh in vals[cblock]]
        print(f'PHASE "{hall}" {" ".join(cell)} {lam} {rel}:{cblock + 1}')
        n_phases += 1
        for h0 in heads:
            cols = lines[h0].split()
            ih, imult = cols.index("H"), cols.index("Mult")
            itt = next(i for i, c in enumerate(cols) if c.startswith("2theta"))
            idd = cols.index("d-hkl")
            for j in range(h0 + 1, len(lines)):
                f = lines[j].split()
                if not f:
                    if j > h0 + 1:
                        break
                    continue
                if not f[0].isdigit() or len(f) != len(cols):
                    break
                print("REFL", f[ih], f[ih + 1], f[ih + 2], f[imult], f[itt], f[idd], f"{rel}:{j + 1}")
                n_refl += 1
    print(f"# phases: {n_phases}, reflections: {n_refl}")


main()
