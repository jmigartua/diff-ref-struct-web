#!/usr/bin/env python3
"""Write core/src/bv_table.rs, the bond-valence tables of CrysFML, the crystallographic library of FullProf and of
its program Bond_Str (J. Rodríguez-Carvajal, J. González-Platas, N. A. Katcho, ILL; LGPL-3), from the source of the
PyPI package crysfml 0.8.0 (`pip download --no-deps --no-binary :all: crysfml==0.8.0`, then `repo/CFML`):

- `BVS_Table` (Src/CFML_Tables/Tab_Set_BVST.f90, 247 cations): for each cation, d0 and B0 of s = exp((d0 − d)/B0)
  with each of the 14 anions of `BVS_Anions`, and the reference of each pair (`REF_BVS`, 36 entries);
- `BVS_Anions` and `BVS_Anions_Rion` (Src/CFML_BVS_Tables.f90): the anions and their ionic radii;
- `Chem_Info` (Src/CFML_Tables/Tab_Set_ScatterT.f90, 108 elements): the oxidation states `Oxid` and ionic radii
  `Rion` that `Get_Ionic_Radius` returns for a cation.

    python3 make_bv_table.py <path to repo/CFML> > ../../src/bv_table.rs

The numbers are copied as written in the Fortran source; nothing is typed by hand."""
import json, re, sys

root = sys.argv[1]
def read(rel):
    t = open(f"{root}/{rel}", encoding="latin-1").read()
    # Join Fortran continuation lines ("&" at the end of a line, optionally "&" at the start of the next).
    return re.sub(r"&[ \t]*(![^\n]*)?\n[ \t]*&?", " ", t)

bvst = read("Src/CFML_Tables/Tab_Set_BVST.f90")
tables = read("Src/CFML_BVS_Tables.f90")
chem = read("Src/CFML_Tables/Tab_Set_ScatterT.f90")

num = r"\s*([-+0-9.]+)\s*"
arr = lambda n: r"\[" + ",".join([num] * n) + r"\]"
bvs = re.findall(r'BVS_Par_Type\(\s*"([^"]*)"\s*,\s*' + arr(14) + r"\s*,\s*" + arr(14) + r"\s*,\s*" + arr(14) + r"\s*\)", bvst)
anions = re.findall(r'"([^"]*)"', re.search(r"BVS_ANIONS\s*=\s*\[([^\]]*)\]", tables).group(1))
rion = re.findall(r"[-+0-9.]+", re.search(r"BVS_ANIONS_RION\s*=\s*\[([^\]]*)\]", tables).group(1))
refs = re.findall(r'"([^"]*)"', re.search(r"REF_BVS\s*=\s*\[(.*?)\]", tables, re.S).group(1))
info = re.findall(r'chem_info_type\(\s*"([^"]*)"\s*,\s*"[^"]*"\s*,\s*(\d+)\s*,' + ",".join([num] * 4) + r",\s*" + arr(5) + r"\s*,\s*" + arr(5), chem)
assert len(bvs) == 247 and len(anions) == 14 and len(rion) == 14 and len(refs) == 36 and len(info) == 108, (len(bvs), len(anions), len(rion), len(refs), len(info))

def f(s):
    """A Fortran real as a Rust f64 literal."""
    return s if "." in s else s + ".0"

out = []
w = out.append
w("//! Bond-valence tables of CrysFML, the crystallographic library of FullProf and of its program Bond_Str: d0 and B0")
w("//! of s = exp((d0 − d)/B0) for 247 cations with 14 anions, with their references; the anions' ionic radii; the")
w("//! cations' ionic radii by oxidation state (Bond_Str counts a bond when d ≤ (r_cation + r_anion)(1 + tol)).")
w("//!")
w("//! Source: CrysFML (J. Rodríguez-Carvajal, J. González-Platas, N. A. Katcho, ILL; LGPL-3), the PyPI package crysfml")
w("//! 0.8.0, `Src/CFML_Tables/Tab_Set_BVST.f90` (`BVS_Table`), `Src/CFML_BVS_Tables.f90` (`BVS_Anions`,")
w("//! `BVS_Anions_Rion`, `REF_BVS`) and `Src/CFML_Tables/Tab_Set_ScatterT.f90` (`Chem_Info`: `Oxid`, `Rion`); written by")
w("//! `core/tests/data/make_bv_table.py`. Do not edit. The test `bond_valence` checks them against Bond_Str's outputs.")
w("// Table values are data, not approximations of mathematical constants.")
w("#![allow(clippy::approx_constant)]")
w("")
w("/// The anions of the table, in its column order, as CrysFML writes them (\"O-2\", \"CL-1\").")
w(f"pub const ANIONS: [&str; 14] = [{', '.join(json.dumps(a.strip()) for a in anions)}];")
w("")
w("/// The anions' ionic radii (Å), in the order of `ANIONS`.")
w(f"pub const ANION_RADII: [f64; 14] = [{', '.join(f(r) for r in rion)}];")
w("")
w("/// The references of the parameters, by the number each pair carries (0: unknown).")
w("pub const REFERENCES: [&str; 36] = [")
for r in refs:
    w(f"    {json.dumps(r.strip(), ensure_ascii=False)},")
w("];")
w("")
w("/// A cation (\"LA+3\"), and for each anion of `ANIONS`: d0 (Å, 0 when the pair is not in the table), B0 (Å) and")
w("/// the reference number.")
w("pub type Cation = (&'static str, [f64; 14], [f64; 14], [u8; 14]);")
w("")
w(f"pub static CATIONS: [Cation; {len(bvs)}] = [")
for e in bvs:
    sym, vals = e[0].strip(), e[1:]
    d0, b0, rf = vals[:14], vals[14:28], vals[28:]
    w(f'    ("{sym}", [{", ".join(f(v) for v in d0)}], [{", ".join(f(v) for v in b0)}], [{", ".join(str(int(v)) for v in rf)}]),')
w("];")
w("")
w("/// An element (\"LA\", \"W\" in CrysFML's upper case), its oxidation states (0: none) and the ionic radius (Å) of each.")
w("pub type Ionic = (&'static str, [i8; 5], [f64; 5]);")
w("")
w(f"pub static IONIC_RADII: [Ionic; {len(info)}] = [")
for e in info:
    sym = e[0].strip()
    oxid, radii = e[6:11], e[11:16]
    w(f'    ("{sym.upper()}", [{", ".join(str(int(v)) for v in oxid)}], [{", ".join(f(v) for v in radii)}]),')
w("];")
print("\n".join(out))
