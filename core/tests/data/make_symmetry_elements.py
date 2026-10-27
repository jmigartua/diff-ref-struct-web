#!/usr/bin/env python3
"""Write symmetry_elements_530.txt, the reference for the test `symmetry_elements` (core/tests/symmetry_elements.rs).

cctbx (https://cctbx.github.io, the Computational Crystallography Toolbox; BSD-3-Clause-LBNL) is used only here, in
the test-data step: `pip install cctbx-base`, then
`python3 make_symmetry_elements.py spglib_530.txt > symmetry_elements_530.txt`. The core never depends on it.

For each of the 530 settings of spglib_530.txt, the group cctbx builds from its listed operators (`sgtbx.space_group`
with `expand_smx`; some of spglib's Hall spellings are not cctbx's), and for each of its operators {R|t}, centring
translations included (`all_ops`), cctbx's geometric reading of the operator
(`rot_mx_info`, `translation_part_info`, the analysis behind `show_geometrical_elements`):

    # <hall_number> | <hall symbol>
    r11 .. r33  t1 t2 t3 | type sense | ev1 ev2 ev3 | w1 w2 w3 | p1 p2 p3

R as integers and t in twelfths reduced to 0..11 (the operator as the core keeps it); the rotation type (1, 2, 3,
4, 6, negative when improper) and its sense (+1, -1, 0); the axis direction ev (integers; for a mirror the normal's
axis, for -1 zeros); the intrinsic part w (the screw or glide component) and a point p of the element (cctbx's
origin_shift), both in 1/24.
"""
import sys
from cctbx import sgtbx

def frac24(x, what):
    v = x * 24
    n = round(float(v))
    assert abs(float(v) - n) < 1e-9, (what, x)
    return n

halls = []
for line in open(sys.argv[1]):
    if line.startswith("# ") and "|" in line and line[2].isdigit():
        head, hall = line[2:].split("|")[0:2]
        halls.append((int(head.split()[0]), hall.strip(), []))
    elif halls and not line.startswith("#") and line.strip():
        v = [int(x) for x in line.split()]
        halls[-1][2].append(sgtbx.rt_mx(sgtbx.rot_mx(v[:9]), sgtbx.tr_vec(v[9:12], 12)))
assert len(halls) == 530, len(halls)
from importlib.metadata import version
cctbx_version = version("cctbx-base")
print(f"# cctbx-base {cctbx_version}: rot_mx_info and translation_part_info for the 530 settings of spglib_530.txt")
print("# R, t in twelfths | type sense | axis | intrinsic part in 1/24 | a point of the element in 1/24; written by make_symmetry_elements.py, do not edit")
for h, hall, ops in halls:
    g = sgtbx.space_group()
    for op in ops:
        g.expand_smx(op.new_denominators(g.r_den(), g.t_den()))
    assert g.order_z() == len(ops), (h, g.order_z(), len(ops))
    print(f"# {h} | {hall}")
    for op in g.all_ops():
        op = op.mod_positive()
        r = [v // op.r().den() for v in op.r().num()]
        den = op.t().den()
        t = [((v * 12) // den) % 12 for v in op.t().num()]
        assert all((v * 12) % den == 0 for v in op.t().num()), (h, op)
        ri = op.r().info()
        tp = sgtbx.translation_part_info(op)
        w = [frac24(x, "intrinsic") for x in tp.intrinsic_part().as_rational()]
        p = [frac24(x, "point") for x in tp.origin_shift().as_rational()]
        print(" ".join(map(str, r)), " ", " ".join(map(str, t)), "|", ri.type(), ri.sense(), "|",
              " ".join(map(str, ri.ev())), "|", " ".join(map(str, w)), "|", " ".join(map(str, p)))
