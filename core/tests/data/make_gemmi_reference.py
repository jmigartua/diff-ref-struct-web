#!/usr/bin/env python3
"""Write tests/data/gemmi_structure_factors.txt: structure factors computed by gemmi (an independent library, MPL-2.0,
https://github.com/project-gemmi/gemmi) for the reference structures of L020, read from the tutorial .pcr files in
manuals/fullprof/, for the test `gemmi_structure_factors`.

    pip install gemmi          # test-data step only; gemmi is not a dependency of the core
    python3 make_gemmi_reference.py <repository root> > gemmi_structure_factors.txt

Phases (file, phase number in the file):
- Y2O3, I a -3: 2-Tutorial-Y2O3/pcr-files/y2o3_structure.pcr, phase 1;
- MgAl2O4, F d -3 m: 4-Analysis-Spinel-MgAl2O4/pcr-files/MgAl2O4b_Rietveld.pcr, phase 1;
- Co-doped YBCO, P 4/m m m: 3-Tutorial_Co_doped_YBCO/pcr-files/Cocus-asym.pcr, phase 1;
- quantitative phase analysis: 5-Quantitative_Phase_Analysis/Anatase-Rutile/pcr-files/rutana_tch.pcr, phases 1-2
  (rutile, anatase), and EgyptianMakeup/pcr-files/Egypt_qa_Rietveld.pcr, phases 1-5 (laurionite, phosgenite,
  cerussite, galena, gypsum).

For each phase gemmi builds the group from FullProf's Hall symbol for the .pcr symbol (CrysFML's table, as copied in
core/src/sg_fullprof.rs; gemmi generates the operators itself), places the atoms with chemical occupancy
Occ·M/m (M the order of the group, m the number of distinct images of the atom under gemmi's operators: FullProf's Occ
is chemical occupancy × m/M, manual p. 102), lets gemmi convert them for its sum over all operators
(`change_occupancies_to_crystallographic`, dividing by the order of the site-symmetry group it finds), and B = Biso, and computes, for every unique reflection with d ≥ 0.9 Å:
- F for X-rays with gemmi's IT92 form factors, ions included (`IT92_set_ignore_charge(False)`), without anomalous terms;
- F for neutrons with gemmi's Neutron92 scattering lengths (fm).
gemmi's coefficients for every species are written too, so that the test can tell the table from the summation.
Lines:
    PHASE <name> "<.pcr symbol>" "<Hall>" a b c alpha beta gamma <file>
    ATOM <label> <species> x y z Biso Occ <line>
    XRAY <species> a1 a2 a3 a4 b1 b2 b3 b4 c        (gemmi IT92, as gemmi stores them: single precision)
    NEUTRON <species> b                              (gemmi Neutron92, fm)
    REF h k l Re(F_xray) Im(F_xray) Re(F_neutron) Im(F_neutron)   (electrons; fm)
"""
import math, re, sys
from pathlib import Path
import gemmi

root = Path(sys.argv[1])
FP = "manuals/fullprof"
PHASES = [
    ("Y2O3", "2-Tutorial-Y2O3/pcr-files/y2o3_structure.pcr", 1),
    ("MgAl2O4", "4-Analysis-Spinel-MgAl2O4/pcr-files/MgAl2O4b_Rietveld.pcr", 1),
    ("Co-YBCO", "3-Tutorial_Co_doped_YBCO/pcr-files/Cocus-asym.pcr", 1),
    ("rutile", "5-Quantitative_Phase_Analysis/Anatase-Rutile/pcr-files/rutana_tch.pcr", 1),
    ("anatase", "5-Quantitative_Phase_Analysis/Anatase-Rutile/pcr-files/rutana_tch.pcr", 2),
    ("laurionite", "5-Quantitative_Phase_Analysis/EgyptianMakeup/pcr-files/Egypt_qa_Rietveld.pcr", 1),
    ("phosgenite", "5-Quantitative_Phase_Analysis/EgyptianMakeup/pcr-files/Egypt_qa_Rietveld.pcr", 2),
    ("cerussite", "5-Quantitative_Phase_Analysis/EgyptianMakeup/pcr-files/Egypt_qa_Rietveld.pcr", 3),
    ("galena", "5-Quantitative_Phase_Analysis/EgyptianMakeup/pcr-files/Egypt_qa_Rietveld.pcr", 4),
    ("gypsum", "5-Quantitative_Phase_Analysis/EgyptianMakeup/pcr-files/Egypt_qa_Rietveld.pcr", 5),
]
DMIN = 0.9

# FullProf's Hall symbols, from CrysFML's table as generated into the core.
FULLPROF_HALL = {}
for m in re.finditer(r'\(\d+, "([^"]*)", "((?:[^"\\]|\\")*)"\)', (root / "core/src/sg_fullprof.rs").read_text()):
    FULLPROF_HALL.setdefault(m.group(1).replace(" ", ""), m.group(2).replace('\\"', '"'))


def pcr_phase(rel, number):
    """Symbol, cell and atoms of phase `number` of a .pcr file (atoms with N_t = 0: isotropic)."""
    L = (root / FP / rel).read_text(encoding="latin-1").splitlines()
    sg_lines = [i for i, l in enumerate(L) if "<--Space group symbol" in l]
    i = sg_lines[number - 1]
    symbol = L[i].split("<--")[0].strip()
    atoms = []
    j = i + 2
    while not L[j].startswith("!"):
        f = L[j].split()
        if len(f) >= 11 and re.match(r"[A-Za-z]", f[0]) and re.match(r"[A-Za-z]", f[1]):
            assert int(f[9]) == 0, f"{rel}:{j + 1}: anisotropic atom (N_t = {f[9]}) not handled"
            atoms.append((f[0], f[1], *f[2:7], j + 1))
        j += 1
    k = next(k for k in range(j, len(L)) if "#Cell Info" in L[k])
    cell = L[k + 1].split("#")[0].split()[:6]
    return symbol, cell, atoms


def element_charge(species):
    m = re.match(r"([A-Za-z]{1,2})([+-]\d+)?$", species)
    el = m.group(1).capitalize()
    return el, int(m.group(2)) if m.group(2) else 0


gemmi.IT92_set_ignore_charge(False)
print("# Structure factors by gemmi", gemmi.__version__, "for the reference structures of L020; written by make_gemmi_reference.py. Do not edit.")
for name, rel, number in PHASES:
    symbol, cell, atoms = pcr_phase(rel, number)
    hall = FULLPROF_HALL[symbol.replace(" ", "").upper()]
    ops = gemmi.symops_from_hall(hall)
    all_ops = list(ops)
    M = len(all_ops)
    st = gemmi.SmallStructure()
    st.cell = gemmi.UnitCell(*map(float, cell))
    st.spacegroup_hall = hall
    st.determine_and_set_spacegroup("H")
    assert st.spacegroup is not None and set(map(str, st.spacegroup.operations())) == set(map(str, ops)), name
    print(f'PHASE {name} "{symbol}" "{hall}" {" ".join(cell)} {FP}/{rel}')
    species = []
    for label, sp, x, y, z, biso, occ, line in atoms:
        print(f"ATOM {label} {sp} {x} {y} {z} {biso} {occ} {line}")
        xyz = [float(x), float(y), float(z)]
        images = []
        for op in all_ops:
            p = [v % 1.0 for v in op.apply_to_xyz(xyz)]
            if not any(all(min(abs(a - b), 1 - abs(a - b)) < 1e-3 for a, b in zip(p, q)) for q in images):
                images.append(p)
        el, charge = element_charge(sp)
        site = gemmi.SmallStructure.Site()
        site.label = label
        site.element = gemmi.Element(el)
        site.charge = charge
        site.fract = gemmi.Fractional(*xyz)
        site.occ = float(occ) * M / len(images)
        site.u_iso = float(biso) / (8 * math.pi ** 2)
        st.add_site(site)
        if sp not in species:
            species.append(sp)
    st.setup_cell_images()
    # gemmi sums every operator over the listed sites; its own conversion divides each chemical occupancy by the order of
    # the site-symmetry group, which it finds itself.
    st.change_occupancies_to_crystallographic()
    for sp in species:
        el, charge = element_charge(sp)
        e = gemmi.Element(el)
        it = gemmi.IT92_get_exact(e, charge) if charge else e.it92
        if it is None:
            print(f"XRAY {sp} none")
        else:
            c = it.get_coefs()
            print(f"XRAY {sp} " + " ".join(f"{v:.9g}" for v in c))
        print(f"NEUTRON {sp} {e.neutron92.get_coefs()[0]:.9g}")
    cx = gemmi.StructureFactorCalculatorX(st.cell)
    cn = gemmi.StructureFactorCalculatorN(st.cell)
    hkls = gemmi.make_miller_array(st.cell, st.spacegroup, DMIN, 0, True)
    for h in sorted((tuple(int(v) for v in r) for r in hkls), key=lambda h: (-st.cell.calculate_d(list(h)), h)):
        fx = cx.calculate_sf_from_small_structure(st, list(h))
        fn = cn.calculate_sf_from_small_structure(st, list(h))
        print(f"REF {h[0]} {h[1]} {h[2]} {fx.real:.10g} {fx.imag:.10g} {fn.real:.10g} {fn.imag:.10g}")
