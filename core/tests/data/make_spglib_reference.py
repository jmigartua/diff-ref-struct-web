#!/usr/bin/env python3
"""Write spglib_530.txt, the reference for the test `spglib_530` (core/tests/spglib_530.rs).

spglib (https://spglib.readthedocs.io, BSD-3) is used only here, in the test-data step: `pip install spglib`,
then `python3 make_spglib_reference.py > spglib_530.txt`. The core never depends on it.

For each of spglib's 530 Hall numbers, from `get_spacegroup_type` (spglib's table of settings) and
`get_symmetry_from_database` (spglib's stored operators, a table separate from the Hall symbols):

    # <hall_number> <number> <choice or -> | <hall symbol> | <international> | <international_full> | <point group>
    r11 r12 r13 r21 r22 r23 r31 r32 r33  t1 t2 t3

one line per operator, rotations as integers, translations as integers in twelfths (exact: spglib's are
multiples of 1/6), reduced to 0..11.
"""
import spglib

print(f"# spglib {spglib.__version__}: get_spacegroup_type and get_symmetry_from_database, Hall numbers 1-530")
print("# translations in twelfths; written by make_spglib_reference.py, do not edit")
for h in range(1, 531):
    t = spglib.get_spacegroup_type(h)
    s = spglib.get_symmetry_from_database(h)
    print(f"# {h} {t.number} {t.choice or '-'} | {t.hall_symbol} | {t.international} | "
          f"{t.international_full} | {t.pointgroup_international}")
    for r, tr in zip(s["rotations"], s["translations"]):
        tw = []
        for x in tr:
            v = round(float(x) * 12)
            assert abs(float(x) * 12 - v) < 1e-6, (h, tr)
            tw.append(v % 12)
        print(" ".join(str(int(v)) for v in r.flatten()), " ", " ".join(map(str, tw)))
