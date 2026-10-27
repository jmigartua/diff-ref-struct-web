#!/usr/bin/env python3
"""Write fullprof_sites.txt, a reference of FullProf atom-site data for the core's tests.

Standard library only. Run from this folder:
    python3 make_fullprof_sites.py > fullprof_sites.txt

It reads the FullProf files under ../../../manuals/fullprof (paths relative to this script) and writes two kinds
of record, every value copied as a string exactly as FullProf printed it, each line ending with <file>:<line>
(the file relative to the repository root):

MULT  from the .out files: the "Initial parameters" atom table (header `Atom Ntyp ... Spc Mult`) gives, for each
      atom, its coordinates and the site multiplicity FullProf computed; the phase's Hall symbol is the
      `Hall Symbol:` line of the space-group block before the table (every .out here holds one phase).
          MULT "<hall>" <label> <x> <y> <z> <mult> <file>:<line of the atom>
      In rigid-body/Z-matrix tables (header `P1 P2 ... P8`) the .out states that P1, P2, P3 are x, y, z.

BETA  from the .pcr files (all, recursively): atoms with N_t = 2 (column located by the phase's `!Atom ...`
      comment header, `N_t` or `N_type`), whose codes line is followed by the line
      beta11 beta22 beta33 beta12 beta13 beta23 and its codes line (FullProf manual, LINE 25, appendix).
          BETA "<space group symbol>" <label> <x> <y> <z> <b11> <b22> <b33> <b12> <b13> <b23>
               <c11> <c22> <c33> <c12> <c13> <c23> <file>:<line of the betas line>
      Only phases with Isy = 0 (symmetry from the symbol on the line ending `<--Space group symbol`) and
      Jbt in {0, 2, 3, -3, 10} are read; Jbt = 10 phases use the Atom Typ Mag Vek X Y Z ... N_type layout.
      Skipped: magnetic phases (Jbt = 1, -1, 5, 7, -7, -10, ...), symmetry-mode phases (Jbt = 6, -6), rigid-body
      phases (Jbt = 4, -4), phases with symmetry operators (Isy != 0) or without a plain space-group line, and any
      phase whose atom block does not parse to exactly Nat atoms.
"""
import os
import re
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
ROOT = os.path.join(REPO, "manuals", "fullprof")

NUM = re.compile(r"^[-+]?(\d+\.?\d*|\.\d+)([eEdD][-+]?\d+)?$")


def isnum(s):
    return bool(NUM.match(s))


def rel(path):
    return os.path.relpath(path, REPO).replace(os.sep, "/")


def walk(ext):
    found = []
    for dp, dn, fn in os.walk(ROOT):
        dn.sort()
        for f in sorted(fn):
            if f.lower().endswith(ext):
                found.append(os.path.join(dp, f))
    return sorted(found)


def read_lines(path):
    with open(path, "rb") as fh:
        return fh.read().decode("latin-1").splitlines()


# ---------------------------------------------------------------- .out files: MULT records

def mult_records(path, out, stats):
    lines = read_lines(path)
    halls = [(i, l.split("Hall Symbol:", 1)[1].strip()) for i, l in enumerate(lines) if "Hall Symbol:" in l]
    phases = [i for i, l in enumerate(lines) if re.match(r"\s*=> Number of phases:", l)]
    if phases and int(lines[phases[0]].split(":")[1]) != 1:
        stats["skip"]["out file with more than one phase"] += 1
        return
    for h, line in enumerate(lines):
        t = line.split()
        if not (len(t) > 3 and t[0] == "Atom" and t[1] == "Ntyp" and t[-1] == "Mult"):
            continue
        before = [s for i, s in halls if i < h]
        if not before:
            stats["skip"]["atom table without a Hall symbol before it"] += 1
            continue
        hall = before[-1]
        ix, iy, iz = (t.index("X"), t.index("Y"), t.index("Z")) if "X" in t else (t.index("P1"), t.index("P2"), t.index("P3"))
        started = False
        for k in range(h + 1, len(lines)):
            s = lines[k].split()
            if not s:
                if started:
                    break
                continue
            if s[0].startswith("=>") or s[0].startswith("==>"):
                break
            if lines[k].startswith("  "):
                # indented continuation: B11.. / P9.. headers, Codes:, P9--> lines, SAnn warnings and limits
                continue
            # atom line (label in column 2): label, type, numbers ..., last = Mult (integer)
            if len(s) >= iz + 1 and all(isnum(x) for x in s[2:]) and re.match(r"^\d+$", s[-1]):
                started = True
                out.append('MULT "%s" %s %s %s %s %s %s:%d' % (hall, s[0], s[ix], s[iy], s[iz], s[-1], rel(path), k + 1))
            else:
                break


# ---------------------------------------------------------------- .pcr files: BETA records

def data(line):
    """Tokens of a data line, without trailing `#...` remarks."""
    return line.split("#", 1)[0].split()


def is_comment(line):
    s = line.strip()
    return s == "" or s.startswith("!")


def header_cols(line):
    cols = []
    for tok in line.strip().lstrip("!").split():
        if tok.startswith("/"):
            break
        cols.append(tok.split("/", 1)[0])
        if "/" in tok:
            break
    return cols


def next_data(lines, k, stop):
    while k < stop and is_comment(lines[k]):
        k += 1
    return k


def numbers(tokens, n):
    """First n tokens if they are numbers and anything after them is a remark (`<-...`)."""
    if len(tokens) < n or not all(isnum(x) for x in tokens[:n]):
        return None
    rest = tokens[n:]
    if rest and not rest[0].startswith("<-"):
        return None
    return tokens[:n]


def beta_records(path, out, stats):
    lines = read_lines(path)
    starts = [i for i, l in enumerate(lines) if re.match(r"\s*!\s*Nat\s+Dis\s", l)]
    for p, h in enumerate(starts):
        end = starts[p + 1] if p + 1 < len(starts) else len(lines)
        stats["phases"] += 1
        hdr = lines[h].strip().lstrip("!").split()
        vals = data(lines[h + 1])
        d = dict(zip(hdr, vals))
        try:
            nat, jbt, isy = int(d["Nat"]), int(d["Jbt"]), int(d["Isy"])
        except (KeyError, ValueError):
            stats["skip"]["phase line unreadable"] += 1
            continue
        if nat == 0:
            stats["no atoms"] += 1
            continue
        if jbt in (6, -6):
            stats["skip"]["symmetry modes (Jbt=%d)" % jbt] += 1
            continue
        if jbt in (4, -4):
            stats["skip"]["rigid body (Jbt=%d)" % jbt] += 1
            continue
        if jbt not in (0, 2, 3, -3, 10):
            stats["skip"]["magnetic (Jbt=%d)" % jbt] += 1
            continue
        if isy != 0:
            stats["skip"]["symmetry operators (Isy=%d)" % isy] += 1
            continue
        sg_k = sg = None
        atom_h = None
        for k in range(h + 2, end):
            if sg_k is None and "<--" in lines[k]:
                left, right = lines[k].split("<--", 1)
                if right.strip() == "Space group symbol":
                    sg_k, sg = k, left.strip()
            if re.match(r"\s*!\s*Atom\s", lines[k]):
                atom_h = k
                break
        if sg is None:
            stats["skip"]["no plain '<--Space group symbol' line"] += 1
            continue
        if atom_h is None or atom_h < sg_k:
            stats["skip"]["no '!Atom' header after the space group"] += 1
            continue
        cols = header_cols(lines[atom_h])
        try:
            ix, iy, iz = cols.index("X"), cols.index("Y"), cols.index("Z")
            nt = cols.index("N_t") if "N_t" in cols else cols.index("N_type")
        except ValueError:
            stats["skip"]["'!Atom' header without X Y Z N_t"] += 1
            continue

        def atom_line(k):
            t = data(lines[k])
            return (len(t) > nt and not isnum(t[0]) and not t[0].startswith("<") and not isnum(t[1])
                    and all(isnum(x) for x in t[2:nt + 1]) and re.match(r"^-?\d+$", t[nt]))

        found, recs, ok = 0, [], True
        k = next_data(lines, atom_h + 1, end)
        while found < nat:
            if k >= end or not atom_line(k):
                ok = False
                break
            t = data(lines[k])
            found += 1
            c = next_data(lines, k + 1, end)
            if c >= end or numbers(data(lines[c]), 5) is None:
                ok = False
                break
            if t[nt] == "2":
                b = next_data(lines, c + 1, end)
                cb = next_data(lines, b + 1, end)
                bv = numbers(data(lines[b]), 6) if b < end else None
                cv = numbers(data(lines[cb]), 6) if cb < end else None
                if bv is None or cv is None:
                    ok = False
                    break
                recs.append('BETA "%s" %s %s %s %s %s %s %s:%d' % (
                    sg, t[0], t[ix], t[iy], t[iz], " ".join(bv), " ".join(cv), rel(path), b + 1))
                k = next_data(lines, cb + 1, end)
            else:
                # other atom types carry their own extra lines: move to the next atom line
                k = next_data(lines, c + 1, end)
                if found < nat:
                    while k < end and not lines[k].lstrip().startswith("!-") and not atom_line(k):
                        k = next_data(lines, k + 1, end)
        if not ok:
            stats["skip"]["atom block did not parse to Nat atoms"] += 1
            continue
        stats["read"] += 1
        out.extend(recs)


def main():
    outs, pcrs = walk(".out"), walk(".pcr")
    mstats = {"skip": Counter()}
    pstats = {"skip": Counter(), "phases": 0, "no atoms": 0, "read": 0}
    mult, beta = [], []
    for f in outs:
        mult_records(f, mult, mstats)
    for f in pcrs:
        beta_records(f, beta, pstats)
    skip = pstats["skip"] + mstats["skip"]
    print("# fullprof_sites.txt: FullProf atom-site data, the reference for the core's site tests.")
    print("# Written by make_fullprof_sites.py (core/tests/data) from the files under manuals/fullprof; do not edit.")
    print("# Every value is a string copied from a FullProf file; each line ends with <file>:<line> (repository-relative).")
    print('# MULT "<Hall symbol>" <label> <x> <y> <z> <multiplicity>   (.out, Initial parameters table; Mult column)')
    print('# BETA "<space group symbol>" <label> <x> <y> <z> <b11> <b22> <b33> <b12> <b13> <b23> <c11> <c22> <c33> <c12> <c13> <c23>')
    print("#      (.pcr, atoms with N_t = 2; b = anisotropic betas, c = their codewords; line = the betas line)")
    for r in mult:
        print(r)
    for r in beta:
        print(r)
    reasons = ", ".join("%s %d" % (k, v) for k, v in sorted(skip.items()))
    print("# .pcr phases: %d, without atoms: %d, read: %d" % (pstats["phases"], pstats["no atoms"], pstats["read"]))
    print("# MULT records: %d, BETA records: %d, files read: %d, phases skipped: %d (%s)" % (
        len(mult), len(beta), len(outs) + len(pcrs), sum(skip.values()), reasons))


if __name__ == "__main__":
    main()
