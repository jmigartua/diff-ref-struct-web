#!/usr/bin/env python3
"""Write tests/data/fullprof_structure_factors.txt: structures and the F² FullProf printed for them, copied from the
FullProf files in manuals/fullprof/ with file and line, for the test `fullprof_structure_factors`.

    python3 make_fullprof_structure_factors.py <repository root> > fullprof_structure_factors.txt

Three phases, the files of manuals/fullprof/ that print F² together with the structure they come from:
- PPH3 (triphenylphosphine, P 21/c, 34 atoms, X-rays λ = 1.001 Å): atoms from the last "New parameters" block of
  PPH3-sls.out (the refinement had converged: shifts zero), species from its "Initial parameters" table, the X-ray
  coefficients from its "X-ray scattering coeff." table, the cell from the header of PPH3-sls.hkl, and every reflection
  of PPH3-sls.hkl (h k l, Mult, sinT/lamb, F2).
- Citric acid (P 21/c, 21 atoms, X-rays λ = 1.001 Å): the same files of CitricAcid_example/Rietveld.
- LaMnO3 (P n m a, 4 atoms on special and general positions, neutrons λ = 1.2272 Å): atoms from "Final Complete
  Structure with contribution of all modes" of Lmo_san.out (simulated annealing on symmetry-mode amplitudes), the cell
  from Lmo_san.pcr, and the reflection table of Lmo_san.out (F2cal(Nuc+Mag), no magnetic phase: Ical(mag) = 0).

Numbers are copied as printed (esd's in brackets removed), so the test can take each one's precision from its
decimals. Lines:
    PHASE <name> "<Hall>" a b c alpha beta gamma <xray|neutron> <lambda> <structure file> <reflection file>
    COEF <species> a1 b1 a2 b2 a3 b3 a4 b4 c f0 Z dfp dfpp <line in the structure file>
    ATOM <label> <species> x y z Biso Occ <line in the structure file>
    REF h k l mult sinT/lambda F2 <line in the reflection file>
"""
import re, sys
from pathlib import Path

root = Path(sys.argv[1])
FP = "manuals/fullprof"


def lines(rel):
    return (root / rel).read_text(encoding="latin-1").splitlines()


def nobracket(s):
    return re.sub(r"\(\d+\)", "", s)


def hall_of(out):
    for l in out:
        m = re.match(r"\s*=>\s+Hall Symbol:\s*(.+?)\s*$", l)
        if m:
            return m.group(1)
    raise SystemExit("no Hall symbol")


def rietveld(name, folder, stem):
    """A Rietveld .out with its .hkl (X-rays): PPH3 and citric acid."""
    outf, hklf = f"{FP}/{folder}/{stem}.out", f"{FP}/{folder}/{stem}.hkl"
    out, hkl = lines(outf), lines(hklf)
    lam = next(l.split()[2] for l in out if l.strip().startswith("=> Wavelengths:"))
    cell = re.search(r"CELL:\s+(.*)$", hkl[1]).group(1).split()
    print(f'PHASE {name} "{hall_of(out)}" {" ".join(cell)} xray {lam} {outf} {hklf}')
    # Species: the "Initial parameters" table, "A1     P    0.14126 ...".
    start = next(i for i, l in enumerate(out) if l.strip().startswith("=> Initial parameters"))
    species = {}
    for l in out[start:start + 400]:
        m = re.match(r"^ (\S+)\s+([A-Za-z]{1,2}[+-]?\d*)\s+-?\d+\.\d+\s", l)
        if m:
            species.setdefault(m.group(1), m.group(2))
        if l.strip().startswith("=> The normalized") or l.strip().startswith("=>-------> PROFILE"):
            break
    # X-ray coefficients printed by FullProf.
    i = next(i for i, l in enumerate(out) if "X-ray scattering coeff." in l)
    for j in range(i + 1, i + 30):
        f = out[j].split()
        if len(f) == 14 and re.match(r"[A-Za-z]{1,2}[+-]?\d*$", f[0]):
            print(f"COEF {' '.join(f)} {j + 1}")
        elif f and f[0].startswith("---"):
            break
    # Atoms: the last "New parameters" block.
    i = max(i for i, l in enumerate(out) if l.strip().startswith("=> New parameters"))
    n = 0
    for j in range(i, i + 400):
        m = re.match(r"^ (\S+)\s+(-?\d+\.\d+) \S+ \S+\s+(-?\d+\.\d+) \S+ \S+\s+(-?\d+\.\d+) \S+ \S+\s+(-?\d+\.\d+) \S+ \S+\s+(-?\d+\.\d+)", out[j])
        if m and m.group(1) in species:
            lab = m.group(1)
            print(f"ATOM {lab} {species[lab]} {m.group(2)} {m.group(3)} {m.group(4)} {m.group(5)} {m.group(6)} {j + 1}")
            n += 1
        if n and out[j].strip().startswith("Atom") and "r9" in out[j]:
            break
    assert n == len(species), (name, n, len(species))
    # Reflections: h k l Mult sinT/lamb 2T Fwhm F2 sf2.
    for j, l in enumerate(hkl[3:], 4):
        f = l.split()
        if len(f) == 9:
            print(f"REF {f[0]} {f[1]} {f[2]} {f[3]} {f[4]} {f[7]} {j}")


def lamno3():
    d = f"{FP}/7-FullProf-Tutorial-SymmetryModes/LaMnO3/SimulatedAnnealing"
    outf, pcrf = f"{d}/Lmo_san.out", f"{d}/Lmo_san.pcr"
    out, pcr = lines(outf), lines(pcrf)
    lam = next(l.split()[-1] for l in out if "Wavelength (Angstr.)" in l)
    k = next(i for i, l in enumerate(pcr) if "#Cell Info" in l)
    cell = pcr[k + 1].split("#")[0].split()
    print(f'PHASE LaMnO3 "{hall_of(out)}" {" ".join(cell)} neutron {lam} {outf} {outf}')
    i = next(i for i, l in enumerate(out) if "Final Complete Structure with contribution of all modes" in l)
    for j in range(i + 2, i + 20):
        f = out[j].split()
        if not f or f[0] != "Atom":
            break
        print(f"ATOM {f[1]} {f[2]} {nobracket(f[3])} {nobracket(f[4])} {nobracket(f[5])} {f[6]} {f[7]} {j + 1}")
    i = next(i for i, l in enumerate(out) if "F2cal(Nuc+Mag)" in l)
    for j in range(i + 1, len(out)):
        l = out[j]
        if l.strip().startswith("---") or l.strip().startswith("=>"):
            break
        m = re.match(r"^\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+0\s+(\d+)\s", l)
        if m:
            f = l.split()
            assert f[-4] == "0.0000", l  # Ical(mag): no magnetic contribution
            print(f"REF {f[0]} {f[1]} {f[2]} {f[4]} {f[-1]} {f[-3]} {j + 1}")


print("# Structures and F² printed by FullProf, from manuals/fullprof/; written by make_fullprof_structure_factors.py. Do not edit.")
rietveld("PPH3", "9-SAnn-OrganicCompounds/PPH3/Rietveld", "PPH3-sls")
rietveld("citric-acid", "9-SAnn-OrganicCompounds/CitricAcid_example/Rietveld", "citric-acid-sls_1")
lamno3()
