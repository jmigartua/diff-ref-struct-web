#!/usr/bin/env python3
"""Write core/tests/data/file_kinds.txt: what every file of manuals/fullprof/runs, the tutorials (1-…, 11-…) and
FP_tests is, for the test `file_kinds` (L081, W5).

The kinds and their sources come from a survey of the 2649 files (L081): FullProf's manual (July 2001; Input files
p. 8, Output files pp. 9–12), the FullProf News, the tutorials and the manuals of the suite's other programs. Each rule
below is the signature the survey found in every file of the kind, written here apart from the core's classifier
(`core/src/formats/project.rs`) so that the two can be compared. "data" is a file the core's reader registry reads.

    python3 make_file_kinds.py > file_kinds.txt          (from core/tests/data)
"""
import os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
FP = os.path.normpath(os.path.join(HERE, "../../../manuals/fullprof"))

# id → the source that says what the kind is (the same sources the classifier cites).
SOURCES = {
    "control": "L070 (core/src/control)",
    "fullprof_pcr": "manual p. 8, Input files: CODFIL.pcr",
    "fullprof_pcr_new": "manual p. 8, Input files: CODFIL.new; p. 72 (Pcr = 2)",
    "fullprof_pcr_original": "runs/README.md (runs/run_fullprof.py)",
    "fullprof_background": "manual p. 8, Input files: FILE.bac; p. 73 (Ppl = 2)",
    "fullprof_difference": "manual p. 73 (Ppl = 3)",
    "fullprof_reflections": "manual p. 11, Output files: CODFIL.hkl (Hkl)",
    "fullprof_reflections_hkl_minus2": "not documented (content; the .pcr has Hkl = −2)",
    "fullprof_phase_reflections": "manual p. 9, Input files: CODFILn.hkl (Irf); p. 90 (Jbt = 2, Irf 0 → 2)",
    "fullprof_magnetic_reflections": "manual p. 10, Output files: CODFILn.hkl; p. 90 (Irf ≤ 0, Isy)",
    "fullprof_bad_reflections": "not documented (content)",
    "shelx_hklf4": "8-FullProf-Tutorial-RigidBodies/tutorial-rgb.inf",
    "fullprof_integrated_intensities": "manual p. 9, Input files: CODFIL.int; p. 12, Output files: CODFILn.int; p. 138",
    "fullprof_resolution": "manual p. 9, Input files: MYRESOL.irf",
    "fullprof_listing": "manual p. 9, Output files: CODFIL.out",
    "fullprof_summary": "manual p. 10, Output files: CODFIL.sum",
    "fullprof_log": "manual p. 7 (fp2k CODFIL DATFIL LOGFIL)",
    "fullprof_run_log": "runs/run_fullprof.py",
    "fullprof_profile": "manual p. 9, Output files: CODFIL.prf",
    "fullprof_integrated_profile": "News 19 June 2020",
    "fullprof_symmetry": "manual p. 10, Output files: CODFIL.sym",
    "fullprof_fourier": "manual p. 10, Output files: CODFIL.fou (Fou)",
    "fullprof_fourier_input": "manual p. 10, Output files: CODFILn.inp (Fou = 3, 4)",
    "fullprof_shelx_template": "manual p. 10, Output files: CODFILn.ins (Fou = 2)",
    "fullprof_microstructure": "manual pp. 11–12, Output files: CODFILn.mic",
    "fullprof_distances": "manual p. 11, Output files: CODFILn.dis",
    "fullprof_phase_profile": "manual p. 12, Output files: CODFILn.sub",
    "fullprof_atoms": "manual p. 12, Output files: CODFILn.atm",
    "fullprof_studio": "News 14 July 2004",
    "fullprof_ortep": "manual p. 39 (Jbt = 4: x#.ortep, reduced co-ordinates and temperature parameters for ORTEP; believed its present name CODFILn.ort: its ATOM and BETA lines)",
    "fullprof_mcif": "News 19 May 2016",
    "fullprof_cif": "News 2002 (Rpa = −1)",
    "fullprof_profile_cif": "News 27 November 2019",
    "fullprof_simple_cif": "not documented (content)",
    "fullprof_distances_cif": "News 22 March 2008",
    "cif": "Hall, Allen & Brown, Acta Cryst. A47, 655 (1991)",
    "vesta": "News 18 May 2015",
    "vesta_grid": "News 22 October 2014",
    "fullprof_microstructure_binary": "News April 2002",
    "fullprof_strain": "News 14 May 2002",
    "fullprof_size": "News 14 May 2002",
    "fullprof_cfl": "Y2O3 tutorial p. 21; docs/bond_str.pdf p. 1",
    "crysfml_cfl": "docs/bond_str.pdf p. 1",
    "bond_valence_summary": "docs/bond_str.pdf p. 3",
    "fullprof_restraints": "News 20 November 2012",
    "mol_tpcr_output": "docs/FullProf-Restraints-RigidBodies.pdf p. 15",
    "fullprof_sa_profile": "News 7 July 2010",
    "fullprof_sequential": "News 10 February 2004",
    "fullprof_sa_solution": "News 27 July 2006 (.so1)",
    "fullprof_sa_log": "not documented (content)",
    "fullprof_sa_configurations": "not documented (content)",
    "background_points": "FAULTS manual p. 32",
    "winplotr_background_points": "not documented (content)",
    "buffer_file": "WinPLOTR manual p. 15",
    "peak_list": "FullProfAPP manual (Export Peaks)",
    "dicvol_input": "not documented (content: WinPLOTR's header)",
    "dicvol_output": "not documented (content: DICVOL's header)",
    "dicvol_condensed": "not documented (content: DICVOL's header)",
    "dicvol_cell": "not documented (content: DICVOL's header)",
    "dicvol_summary": "not documented (content: WDICVOL04 line)",
    "treor_output": "not documented (content: TREOR's header)",
    "check_group_output": "not documented (content: CheckGroup's banner)",
    "formal_charges_output": "not documented (content: Formal_Charges' banner)",
    "zmatrix_fhz": "docs/FullProf-Restraints-RigidBodies.pdf p. 5",
    "datared_input": "News 22 February 2022",
    "faults_input": "FAULTS manual p. 7",
    "magnetic_structure_factors": "News 12 February 2012",
    "get_log_results": "not documented (content: Get_Log's banner)",
    "text_information": "manual p. 8 (fp2k.inf)",
    "text": "the extension",
    "script": "the extension",
    "pdf": "the extension",
}


def kind(path):
    name = os.path.basename(path)
    low = name.lower()
    ext = low.rsplit(".", 1)[1] if "." in low else ""
    raw = open(path, "rb").read()
    if ext in ("ggrid", "pgrid"):
        return "vesta_grid"
    if ext == "bin" and (b"Strain Diagram" in raw[:200] or b"Size Diagram" in raw[:200]):
        return "fullprof_microstructure_binary"
    if ext == "pdf":
        return "pdf"
    if ext in ("zip", "pyc"):
        return "unknown"
    t = raw.decode("latin-1")
    lines = t.splitlines()
    l0 = lines[0].strip() if lines else ""
    l1 = lines[1].strip() if len(lines) > 1 else ""
    if low.endswith(".pcr.in"):
        return "fullprof_pcr_original"
    if low.endswith(".run.log") or low == "run.log":
        return "fullprof_run_log"
    if re.search(r'^\s*format\s*=\s*"diff-ref-struct ', t, re.M):
        return "control"
    if any(l.startswith("!Job Npr") for l in lines):
        return "fullprof_pcr_new" if ext == "new" else "fullprof_pcr"
    if ext == "log" and "L O G   F I L E" in t[:4000]:
        return "fullprof_log"
    if "PROGRAM FullProf" in t[:4000] and ext in ("out", "sum"):
        return {"out": "fullprof_listing", "sum": "fullprof_summary"}[ext]
    if ext == "dis":
        return "fullprof_distances"
    if ext == "bvs":
        return "bond_valence_summary"
    if ext == "prf":
        return "fullprof_integrated_profile" if l1.startswith("Sinthet/lamb") else "fullprof_profile"
    if ext == "hkl":
        if l0.startswith("Pattern#") and "Lambda:" in l0:
            return "fullprof_phase_reflections"
        if l0.startswith("Pattern#") or l0.startswith(">Phase no.:"):
            return "fullprof_reflections"
        if l0.startswith(">") and "F2" in l0:
            return "fullprof_reflections_hkl_minus2"
        if l0.startswith("! Phase No.:"):
            return "fullprof_magnetic_reflections"
        if l0.startswith("LIST OF REFLECTIONS WITH DELTA/SIGMA"):
            return "fullprof_bad_reflections"
        return "shelx_hklf4"
    if ext == "bac":
        if "difference pattern" in l0.lower():
            return "fullprof_difference"
        return "fullprof_background"
    simple = {
        "sym": "fullprof_symmetry", "fou": "fullprof_fourier", "mic": "fullprof_microstructure", "sub": "fullprof_phase_profile",
        "atm": "fullprof_atoms", "fst": "fullprof_studio", "ort": "fullprof_ortep", "mcif": "fullprof_mcif", "vesta": "vesta", "str": "fullprof_strain",
        "siz": "fullprof_size", "irf": "fullprof_resolution", "int": "fullprof_integrated_intensities", "spr": "fullprof_sa_profile",
        "seq": "fullprof_sequential", "so1": "fullprof_sa_solution", "so2": "fullprof_sa_solution", "anl": "fullprof_sa_log",
        "bgr": "background_points", "buf": "buffer_file", "aps": "peak_list", "dic": "dicvol_input", "ind": "dicvol_output",
        "ord": "dicvol_output", "ucp": "dicvol_condensed", "cry": "dicvol_cell", "cnd": "treor_output", "spg": "check_group_output",
        "fch": "formal_charges_output", "fhz": "zmatrix_fhz", "red": "datared_input", "flts": "faults_input",
        "msf": "magnetic_structure_factors", "inf": "text_information", "txt": "text", "md": "text", "py": "script",
        "sh": "script", "bat": "script", "drc": "control",
    }
    if ext == "inp":
        return "fullprof_fourier_input" if l0.startswith("! Template file to be used with the program FOURIER") else "unknown"
    if ext == "ins":
        return "fullprof_shelx_template" if l0.startswith("TITL") else "unknown"
    if ext == "cfg":
        return "fullprof_sa_configurations" if l0.startswith("! NT") else "unknown"
    if ext == "res":
        return "get_log_results" if "Program Get_Log" in t[:2000] else "unknown"
    if ext == "cfl":
        return "fullprof_cfl" if "CFL-file generated from FullProf" in t[:400] else "crysfml_cfl"
    if ext == "tpcr":
        if l0.startswith("FILE with lines for soft distance"):
            return "fullprof_restraints"
        return "mol_tpcr_output" if "PROGRAM MOL_TPCR" in t[:4000] else "unknown"
    if ext == "sum":
        return "dicvol_summary" if "WDICVOL04" in t else "unknown"
    if ext == "cif":
        if "Magnetic CIF file generated by FullProf" in t[:400]:
            return "fullprof_mcif"
        if "FullProf-generated Simple CIF" in l1:
            return "fullprof_simple_cif"
        if "FullProf-generated CIF" in l1:
            return "fullprof_cif"
        if any(l.startswith("# Powder diffraction pattern number") for l in lines[:5]):
            return "fullprof_profile_cif"
        if "UNIT CELL INFORMATION" in t[:3000] and "MOLECULAR GEOMETRY" in t:
            return "fullprof_distances_cif"
        if any(i in t for i in ("_pd_meas_counts_total", "_pd_meas_intensity_total", "_pd_proc_intensity_total", "_pd_proc_intensity_net", "_pd_calc_intensity_total")):
            return "data"
        return "cif" if "_cell_length_a" in t else "unknown"
    if ext == "xy" and l0.startswith("! Data file"):
        return "winplotr_background_points"
    if ext in ("dat", "xys", "xye", "xy", "xyd", "gss"):
        return "data"
    if ext in simple:
        return simple[ext]
    return "unknown"


def main():
    files = []
    for top in sorted(os.listdir(FP)):
        if top == "runs" or top == "FP_tests" or re.match(r"\d+-", top):
            for d, _, fs in os.walk(os.path.join(FP, top)):
                for f in fs:
                    files.append(os.path.relpath(os.path.join(d, f), FP))
    files.sort()
    used = set()
    rows = []
    for f in files:
        k = kind(os.path.join(FP, f))
        used.add(k)
        rows.append(f"{f}\t{k}")
    print("# What each file of manuals/fullprof (runs, tutorials, FP_tests) is (L081); written by make_file_kinds.py.")
    print("# Kinds and their sources: 'kind <id>\\t<source>'; then '<path>\\t<id>'.")
    for k in sorted(used):
        if k in SOURCES:
            print(f"kind {k}\t{SOURCES[k]}")
    print("\n".join(rows))


main()
