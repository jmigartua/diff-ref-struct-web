#!/usr/bin/env python3
"""What FullProf reads from each data file of manuals/fullprof/, for the test `data_formats` (L062).

For every .prf of a single pattern whose .pcr names the data file and its format (Ins), a line
    READ <format id> <Ins> <data file> <prf> <first profile line> <points> <SumwYobsSQ> <N-P+C> <Iwg>
The .prf gives the points FullProf read (2θ, Yobs); the .out written with it, when there is one, gives
SumwYobsSQ = Σ y²/σ² over the non-excluded points, in its last block "RELIABILITY FACTORS WITH ALL NON-EXCLUDED
POINTS" (4 significant digits), which tests σ ('-' without an .out), and the .pcr's weighting Iwg (0: w = 1/σ²; 1: maximum
likelihood, w = 1/Ycal·y/σ²). And for every data file of a .pcr without a
.prf, a line
    FILE <format id> <Ins> <data file>
(recognition and reading only).

    python3 make_fullprof_reads.py <repository root> > fullprof_reads.txt
"""
import re
import sys
from pathlib import Path

root = Path(sys.argv[1])
FP = root / "manuals/fullprof"
IDS = {0: "free_steps", 1: "ill_d1a_d2b_old", 2: "ill_d1b_old", 3: "ill_d1b_d20", 4: "brookhaven", 5: "two_axis", 6: "ill_d2b",
       7: "ill_d4_d20l", 8: "psi_dmc", 9: "socabim_uxd", 10: "xys", 11: "variable_time", 12: "gsas", 14: "isis_tof_banks"}


def rel(p):
    return str(p.relative_to(root))


def pcr_ins(pcr):
    lines = pcr.read_text(encoding="latin-1").splitlines()
    for i, l in enumerate(lines):
        if l.startswith("!Ipr") and "Ins" in l:
            head, vals = l[1:].split(), lines[i + 1].split()
            try:
                return int(vals[head.index("Ins")])
            except (ValueError, IndexError):
                return None
    return None


def pcr_iwg(pcr):
    lines = pcr.read_text(encoding="latin-1").splitlines()
    for i, l in enumerate(lines):
        if l.startswith("!Job") and "Iwg" in l:
            head, vals = l[1:].split(), lines[i + 1].split()
            try:
                return int(vals[head.index("Iwg")])
            except (ValueError, IndexError):
                return 0
    return 0


def data_of(pcr):
    """The data file of a single-pattern .pcr: FullProf's comment "DAT-file:", else CODFIL.dat."""
    text = pcr.read_text(encoding="latin-1")
    if re.search(r"^!\s*Files?\s*=>.*DAT-files?:.*,.*DAT", text, re.M) or "NPATT" in text:
        return None  # several patterns
    m = re.search(r"DAT-file:\s*([^,\s]+)", text)
    names = [m.group(1)] if m else []
    names.append(pcr.stem)
    for n in names:
        for cand in (pcr.parent / n, pcr.parent / (n + ".dat")):
            # Integrated intensities (.int, simulated annealing) are not powder patterns.
            if cand.is_file() and cand.suffix.lower() not in (".pcr", ".int"):
                return cand
    return None


def profile(prf):
    lines = prf.read_text(encoding="latin-1").splitlines()
    start = next((i + 2 for i, l in enumerate(lines) if l.startswith(" 2Theta") and "Yobs" in l), None)
    if start is None:
        return None
    n = 0
    for l in lines[start - 1:]:
        f = l.split("\t")
        if len(f) < 5 or not f[1].strip():
            break
        n += 1
    return start, n


def out_sums(out):
    if not out.exists():
        return "-", "-"
    text = out.read_text(encoding="latin-1")
    blocks = text.split("RELIABILITY FACTORS WITH ALL NON-EXCLUDED POINTS")
    if len(blocks) < 2 or "Patt#:  2" in text:
        return "-", "-"
    b = blocks[-1]
    npc = re.search(r"N-P\+C:\s*(\d+)", b)
    sums = re.search(r"SumwYobsSQ.*\n\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", b)
    if not npc or not sums:
        return "-", "-"
    return sums.group(4), npc.group(1)


print("# What FullProf reads from the data files of manuals/fullprof/, for the test data_formats; written by")
print("# make_fullprof_reads.py. Do not edit. READ <id> <Ins> <data> <prf> <first line> <points> <SumwYobsSQ> <N-P+C> <Iwg>;")
print("# FILE <id> <Ins> <data> (no .prf).")
seen = set()
for prf in sorted(FP.rglob("*.prf")):
    stem = prf.stem
    for suffix in ("_calc",):
        if stem.endswith(suffix) and prf.with_name(stem[: -len(suffix)] + ".pcr").exists():
            stem = stem[: -len(suffix)]
    pcr = prf.with_name(stem + ".pcr")
    if not pcr.exists():
        continue
    ins = pcr_ins(pcr)
    data = data_of(pcr)
    p = profile(prf)
    if ins not in IDS or data is None or p is None:
        continue
    total, npc = out_sums(prf.with_suffix(".out"))
    print(f"READ {IDS[ins]} {ins} {rel(data)} {rel(prf)} {p[0]} {p[1]} {total} {npc} {pcr_iwg(pcr)}")
    seen.add(data)
for pcr in sorted(FP.rglob("*.pcr")):
    ins = pcr_ins(pcr)
    data = data_of(pcr)
    if ins in IDS and data is not None and data not in seen:
        print(f"FILE {IDS[ins]} {ins} {rel(data)}")
        seen.add(data)
