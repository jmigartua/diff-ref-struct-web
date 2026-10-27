#!/usr/bin/env python3
"""Write tests/data/fullprof_profiles.txt: the index of the .prf files of manuals/fullprof/ (and of the mode course's
case A in manuals/bcs, L090) for the test
`fullprof_profiles`, each with the .pcr whose model rebuilds it, or the reason it is left out.

    python3 make_fullprof_profiles.py <repository root> > fullprof_profiles.txt

A .prf (written by FullProf with Prf = 1 or 3) holds, point by point, 2θ, Yobs, Ycal, Yobs − Ycal and Backg, and then
the reflection positions with their h k l. Lines:
    COMPARE <prf> <pcr> <first profile line> <points> <class>
        class: `exact` (converged Rietveld or profile matching: the .pcr holds the parameters of the .prf),
               `lebail` (Le Bail: the intensities FullProf wrote in its .hkl after the last cycle, while the .prf
                        was drawn with those of the cycle before; weak reflections differ),
               `inconsistent` (the .prf was not calculated from the state that the .pcr and .out record).
        For the FullProf runs of runs/ (L034) and runs/modes/ (modes_runs.py, L093) the reference is CODE_calc.prf,
        rebuilt from CODE.pcr (see below).
    DATA <prf> <data file> <Ins> <first profile line> <points>    (the observed pattern alone: data readers)
    SKIP <prf> <reason>
Each line ends with the .prf line where the profile begins, so every number can be found in the file.
"""
import re, sys
from pathlib import Path

root = Path(sys.argv[1])
FP = root / "manuals/fullprof"


def rel(p):
    return str(p.relative_to(root))


def profile_start(lines):
    for i, l in enumerate(lines):
        if (l.startswith(" 2Theta") or l.startswith(" T.O.F.")) and "Yobs" in l:
            return i + 2  # 1-based line of the first point
    return None


def npoints(lines, start):
    n = 0
    for l in lines[start - 1:]:
        f = l.split("\t")
        if len(f) < 5 or not f[1].strip():
            break
        n += 1
    return n


INCONSISTENT = {
    # Its profile corresponds to Zero ≈ 0.0207 and other intensities at high angle, while its .pcr, .out, .prf header
    # and the .prf's own reflection positions give Zero = 0.01671; the reflection intensities and widths of its .out
    # agree with the core (see the test).
    "7-FullProf-Tutorial-SymmetryModes/LaMnO3/LaMnO3ref.prf",
}
# FullProf runs with anisotropic size or strain (runs/README.md, L038, L216).
# The quartic form of NdNi-LTOa is calculated (P25, L211); its pattern is not, for its magnetic phase.
ANISOTROPIC = {
    "NdNi-LTOa": "magnetic phase (phase 2, Jbt = −1, NdNi_LTOa.pcr): magnetic structure factors are not calculated yet; the anisotropic strain of its nuclear phase (Str = 1, strain model 3) is checked reflection by reflection in fullprof_anisotropic_strain",
    "NdNi-LTOb": "FullProf's strain model 6 (Str = 0, Str1–Str2), described in words, not by an equation (manual pp. 28–30): not calculated (L211)",
    # The rest of P25 (L216, L267): Pd3MnD's antiphase domains and prba10k's needles are calculated and compared;
    # highpress's platelets are, but not its peak shifts; lac and cuf1k have a magnetic phase, their nuclear phase's
    # widths are checked in fullprof_anisotropic_size.
    "highpress": "hkl-dependent peak shifts (Sol = 1, ModS = 1, Shift1 −0.02439 along [0 0 1]: highpress.pcr, lines 81–84) are P26, not calculated (L267); its platelet size is (Size-Model 1), and its user scattering factors are checked in pattern_terms (L239)",
    "lac": "magnetic phase (phase 2, Jbt = 1, lac.pcr): magnetic structure factors are not calculated yet; the tetragonal strain of its nuclear phase is checked reflection by reflection in fullprof_anisotropic_size (L267)",
    "cuf1k": "magnetic phase (phase 2, Jbt = 1, cuf1k.pcr): magnetic structure factors are not calculated yet; the monoclinic strain of its nuclear phase is checked reflection by reflection in fullprof_anisotropic_size (L267)",
}
LEBAIL = {"7-FullProf-Tutorial-SymmetryModes/LaMnO3/LBF_lmopr.prf"}
# The TOF variants of tof_runs.py (P05, L331, L334): NCY = 0 from a refined base with one change; FullProf rewrote each
# variant's .pcr after its one cycle (NaN where the cycle diverged), so its .prf belongs to VARIANT.pcr.in, the input
# `tof_runs.py inputs` writes again (arg_si_abs2 and arg_si_abs2b stopped with "Singular matrix" and have no .prf).
TOF_VARIANTS = {"arg_si_plain", "arg_si_plain90", "arg_si_abs1", "arg_si_abs3", "arg_si_ext10", "arg_si_ext100", "arg_si_invd", "alum_noq", "lamn_pol_irf"}
# The runs of version2_runs.py (P41, P30, P23; L345, L346): its variants are NCY = 0 from a refined base, compared with
# VARIANT.pcr.in as the TOF variants; the runs the core does not calculate yet are left out with their reason.
V2_VARIANTS = {"urea_notls", "pdd_t", "pdd_l", "pdd_t_rot", "nac_3t2_debye", "y2o3_bac_b1b2", "y2o3_bac_b3b4", "y2o3_more2", "y2o3_more3",
               "y2o3_more2_two", "y2o3_poly12_b1", "y2o3_poly12_b2", "y2o3_poly12_b7", "y2o3_poly12_b9", "y2o3_poly12_b10", "y2o3_poly12_b11", "y2o3_poly12_b12",
               "yb_nor2_one", "yb_nor2", "yb_nor2_g2", "yb_nor3", "citric_core"}
V2_TLS = {"urea", "urea_conv", "ureac", "ureacyl", "urea_test", "urea_notls", "hmt", "pdd", "pdd_t", "pdd_l", "pdd_t_rot"}
V2_SKIP = {
    "y2o3_filter": "Nba = -2, the Fourier filtering of the background: not reproduced (background_terms, L353), not calculated",
}

# The mode-crystallography course's FullProf runs (manuals/bcs, R29, L090): case A, Sr2ScSbO6 refined by symmetry
# modes and by coordinates, both converged (their .pcr is the file FullProf ran; FullProf 8.20 reproduces them).
MODES = root / "manuals/bcs/mode-crystallography-manual/course-um6p-2024/04_tools_materials"

print("# Index of the .prf files of manuals/fullprof/ and of the mode course (manuals/bcs) for the test fullprof_profiles; written by make_fullprof_profiles.py. Do not edit.")
for prf in sorted(MODES.rglob("*.prf")) if MODES.is_dir() else []:
    lines = prf.read_text(encoding="latin-1").splitlines()
    start = profile_start(lines)
    print(f"COMPARE {rel(prf)} {rel(prf.with_suffix('.pcr'))} {start} {npoints(lines, start)} exact")
for prf in sorted(FP.rglob("*.prf")):
    name = str(prf.relative_to(FP))
    lines = prf.read_text(encoding="latin-1").splitlines()
    start = profile_start(lines)
    if start is None:
        print(f"SKIP {rel(prf)} integrated intensities (columns Sinthet/lamb, Gobs, Gcal: simulated annealing), no profile")
        continue
    n = npoints(lines, start)
    pcr = prf.with_suffix(".pcr")
    m = re.match(r"(la0)(\d{3})\.prf$", prf.name)
    if m:
        print(f"DATA {rel(prf)} {rel(prf.with_suffix('.dat'))} 3 {start} {n}")
        own = prf.with_name(f"la0_{m.group(2)}.pcr")
        line = next(i for i, l in enumerate(own.read_text(encoding="latin-1").splitlines(), 1) if l.startswith("Magnetic phase"))
        print(f"SKIP {rel(prf)} magnetic phase (Jbt = 1, line {line} of {own.name}): magnetic structure factors are not calculated yet")
        continue
    if prf.parent.parent.name == "window_experiments":
        # FullProf run on modified copies of a refined run, NCY = 0, to measure its peak windows (L041): the .prf
        # belongs to the values read, input.pcr (identical to the .pcr FullProf wrote back where nothing is refined).
        print(f"COMPARE {rel(prf)} {rel(prf.parent / 'input.pcr')} {start} {n} exact")
        continue
    if prf.parent.name == "background" and prf.parent.parent.name == "runs":
        # FullProf's background runs (background_runs.py, L079): their Backg and chi2 are compared in background_points
        # (L082), where the spline's rule was found.
        print(f"SKIP {rel(prf)} a background run (background_runs.py, L079): Backg and chi2 compared in background_points (L082)")
        continue
    if prf.parent.parent.name == "runs" and prf.stem.endswith("_conv") and prf.parent.name != "rigid":
        # Converged refinements from CODE.pcr.in (converge_runs.py, L052; restraint_lebail_runs.py, L056): no pass 2
        # was run on them; their patterns are compared through χ² at FullProf's values in the tests named.
        tests = {"Sr-ox-test": "fullprof_restraints", "Sr-ox-rest": "fullprof_restraints", "LaMnO3-3t2": "fullprof_restraints",
                 "K2NaAgCN": "fullprof_restraints", "MgAl2O4-LeBail": "none (Le Bail: fullprof_le_bail uses the tutorial's file)",
                 "Co-YBCO-LeBail": "none (Le Bail: fullprof_le_bail uses the tutorial's file)"}.get(prf.parent.name, "fullprof_refinements")
        print(f"SKIP {rel(prf)} a converged refinement (CODE_conv, L052, L056): drawn from the parameters of its cycle before the last, and no CODE_conv_calc pass exists; its pattern is checked through chi2 at FullProf's values in: {tests}")
        continue
    in_runs = prf.parent.parent.name == "runs" or prf.parent.parent.relative_to(FP).as_posix() == "runs/modes"
    joint = re.match(r"(.+?)(_calc)?_(\d+)\.prf$", prf.name)
    code_pcr = prf.parent / f"{joint.group(1)}.pcr" if joint else None
    if in_runs and joint and code_pcr.exists() and any(l.startswith("NPATT") for l in code_pcr.read_text(encoding="latin-1").splitlines()[:6]):
        # A joint run (joint_runs.py, P06, L320): one .prf per pattern; pass 2's CODE_calc_k.prf belongs to pattern k
        # of CODE.pcr, the values FullProf read (the index names it CODE.pcr#k).
        if not joint.group(2):
            print(f"SKIP {rel(prf)} pass 1 of a joint run (the refinement): drawn from the parameters of its cycle before the last; the reference is {joint.group(1)}_calc_{joint.group(3)}.prf")
        else:
            print(f"COMPARE {rel(prf)} {rel(code_pcr)}#{joint.group(3)} {start} {n} exact")
        continue
    base = prf.stem[: -len("_calc")] if prf.stem.endswith("_calc") else prf.stem
    if in_runs and prf.parent.name in ("rigid", "texture", "background_terms"):
        if prf.stem in V2_SKIP:
            print(f"SKIP {rel(prf)} {V2_SKIP[prf.stem]}")
            continue
    if in_runs and prf.stem in V2_VARIANTS:
        print(f"COMPARE {rel(prf)} {rel(prf.with_name(prf.stem + '.pcr.in'))} {start} {n} exact")
        continue
    if in_runs and prf.stem in TOF_VARIANTS:
        print(f"COMPARE {rel(prf)} {rel(prf.with_name(prf.stem + '.pcr.in'))} {start} {n} exact")
        continue
    if in_runs and not prf.stem.endswith("_calc"):
        print(f"SKIP {rel(prf)} pass 1 of a run (the refinement): drawn from the parameters of its cycle before the last; the reference is {prf.stem}_calc.prf")
        continue
    text = pcr.read_text(encoding="latin-1")
    if prf.name == "la0.prf":
        print(f"SKIP {rel(prf)} magnetic phase (Jbt = 1, line {next(i for i, l in enumerate(text.splitlines(), 1) if l.startswith('Magnetic phase'))} of la0.pcr): magnetic structure factors are not calculated yet")
        continue
    if "Simulated Annealing" in lines[0] and re.search(r"^\s*1\s+0\s+1\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+0\s+7\s+3", text, re.M):
        print(f"SKIP {rel(prf)} simulated annealing on integrated intensities (Cry = 3 in the .pcr): the profile is FullProf's cluster profile (.spr), not a model to rebuild")
        continue
    if prf.stem.endswith("_calc") and prf.parent.name in ANISOTROPIC:
        print(f"SKIP {rel(prf)} {ANISOTROPIC[prf.parent.name]}")
        continue
    if prf.stem.endswith("_calc") and prf.parent.name == "NdNi-HTO":
        print(f"SKIP {rel(prf)} magnetic phase (phase 2, Jbt = 1, Nvk = 1, NdNi_HTO.pcr): magnetic structure factors are not calculated yet; its nuclear phase is checked in fullprof_microstructure")
        continue
    if prf.stem.endswith("_calc"):
        # runs/: FullProf asked for NCY = 0 still makes one cycle ("Number of cycles: 1" in CODE_calc.out) and writes
        # the shifted values back into CODE_calc.pcr; its outputs belong to the values it read, those of CODE.pcr.
        pcr = prf.with_name(prf.stem[: -len("_calc")] + ".pcr")
    cls = "inconsistent" if name in INCONSISTENT else "lebail" if name in LEBAIL else "exact"
    print(f"COMPARE {rel(prf)} {rel(pcr)} {start} {n} {cls}")
