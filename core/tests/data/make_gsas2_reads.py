#!/usr/bin/env python3
"""What GSAS-II reads from each example data file, for the test `data_formats` (L062).

Runs GSAS-II's own importers (GSASII/imports/G2pwd_*.py) on the files below, through GSAS-II's scripting loader
(GSASIIscriptable.import_generic, every bank), and writes core/tests/data/gsas2_reads.txt: per bank, the number of
points, the sums of x, y and w = 1/σ², the first and last points, and the wavelengths GSAS-II takes from the file.

Usage (from the repository root; GSAS-II needs numpy and scipy):
    python3 core/tests/data/make_gsas2_reads.py <GSAS-II checkout>
The checkout used for the committed table: https://github.com/AdvancedPhotonSource/GSAS-II, commit fe944cd.
"""
import math
import os
import re
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))

# (file, GSAS-II importer module, class, the registry's format id)
FILES = [
    ('manuals/gsas2/testinp/PBSO4.CWN', 'G2pwd_fxye', 'GSAS_ReaderClass', 'gsas'),
    ('manuals/gsas2/testinp/PBSO4.XRA', 'G2pwd_fxye', 'GSAS_ReaderClass', 'gsas'),
    ('manuals/gsas2/tutorials/CWNeutron/data/garnet.raw', 'G2pwd_fxye', 'GSAS_ReaderClass', 'gsas'),
    ('manuals/gsas2/tutorials/MCsimanneal/data/Quinoli.gda', 'G2pwd_fxye', 'GSAS_ReaderClass', 'gsas'),
    ('manuals/gsas2/tutorials/MDtutorials/CWInstDemo/data/NIST660CBI.gsas', 'G2pwd_fxye', 'GSAS_ReaderClass', 'gsas'),
    ('manuals/gsas2/tutorials/MDtutorials/k_vec_tutorial/data/Er2Ge2O7_op2112_Ge113_2K.gss', 'G2pwd_fxye', 'GSAS_ReaderClass', 'gsas'),
    ('manuals/gsas2/tutorials/Magnetic-I/data/LaMnO3_50k.gsas', 'G2pwd_fxye', 'GSAS_ReaderClass', 'gsas'),
    ('manuals/gsas2/tutorials/Magnetic-IV/data/PrSrMnO.gda', 'G2pwd_fxye', 'GSAS_ReaderClass', 'gsas'),
    ('manuals/gsas2/tutorials/RMCProfile-I/data/sf6_190gsas.dat', 'G2pwd_fxye', 'GSAS_ReaderClass', 'gsas'),
    ('manuals/gsas2/tutorials/RMCProfile-II/data/gem05984.gss', 'G2pwd_fxye', 'GSAS_ReaderClass', 'gsas'),
    ('manuals/gsas2/tutorials/RMCProfile-III/data/105870-pe-GaPO4_ed.gsas', 'G2pwd_fxye', 'GSAS_ReaderClass', 'gsas'),
    ('manuals/gsas2/tutorials/TOF Calibration/data/PG3_17541.gsa', 'G2pwd_fxye', 'GSAS_ReaderClass', 'gsas'),
    ('manuals/gsas2/tutorials/TOF Sequential Single Peak Fit/data/435.gda', 'G2pwd_fxye', 'GSAS_ReaderClass', 'gsas'),
    ('manuals/gsas2/tutorials/TOF-CW Joint Refinement/data/11BM_NAC.fxye', 'G2pwd_fxye', 'GSAS_ReaderClass', 'gsas'),
    ('manuals/fullprof/FP_tests/CeO2_PEARL.dat', 'G2pwd_fxye', 'GSAS_ReaderClass', 'gsas'),
    ('manuals/fullprof/FP_tests/arg_si.dat', 'G2pwd_fxye', 'GSAS_ReaderClass', 'gsas'),
    ('manuals/fullprof/FP_tests/prbaco10k.dat', 'G2pwd_fxye', 'GSAS_ReaderClass', 'gsas'),
    ('manuals/fullprof/FP_tests/Dy2BaCuO5_1p5K/44329-1_10raw.gss', 'G2pwd_fxye', 'GSAS_ReaderClass', 'gsas'),
    ('manuals/gsas2/tutorials/MDtutorials/CWInstDemo/data/LaB6_Jan2018.raw', 'G2pwd_BrukerRAW', 'raw_ReaderClass', 'bruker_raw_v3'),
    ('manuals/gsas2/tutorials/StackingFaults-II/data/Keokuk kaolinite.RAW', 'G2pwd_BrukerRAW', 'raw_ReaderClass', 'bruker_raw_v3'),
    ('manuals/gsas2/tutorials/StackingFaults-III/data/KGA1-CONT.RAW', 'G2pwd_BrukerRAW', 'raw_ReaderClass', 'bruker_raw_v3'),
    ('manuals/gsas2/tutorials/AutoBkg/data/test_RampDown-00000.xye', 'G2pwd_xye', 'xye_ReaderClass', 'xye'),
    ('manuals/gsas2/tutorials/MCsimanneal/data/x3lactos.xye', 'G2pwd_xye', 'xye_ReaderClass', 'xye'),
    ('manuals/gsas2/tutorials/Magnetic-III/data/Ba6Co6.xye', 'G2pwd_xye', 'xye_ReaderClass', 'xye'),
    ('manuals/gsas2/tutorials/Magnetic-V/data/PG3_20254-2_10K.xye', 'G2pwd_xye', 'xye_ReaderClass', 'xye'),
    ('manuals/gsas2/tutorials/Magnetic-II/data/Cr2WO6_T4K.dat', 'G2pwd_xye', 'xye_ReaderClass', 'xye'),
    ('manuals/fullprof/FP_tests/ceo2_dl.dat', 'G2pwd_xye', 'xye_ReaderClass', 'xye'),
    ('manuals/fullprof/runs/Y2O3/y2o3.dat', 'G2pwd_FP', 'fp_ReaderClass', 'free_steps'),
    ('manuals/fullprof/FP_tests/b8p.dat', 'G2pwd_FP', 'fp_ReaderClass', 'free_steps'),
    ('manuals/fullprof/8-FullProf-Tutorial-RigidBodies/Urea/urea.dat', 'G2pwd_FP', 'fp_ReaderClass', 'free_steps'),
    ('manuals/fullprof/FP_tests/highpress.dat', 'G2pwd_FP', 'fp_ReaderClass', 'free_steps'),
    ('manuals/fullprof/7-FullProf-Tutorial-SymmetryModes/LaMnO3/LaMnO3_prof.cif', 'G2pwd_CIF', 'CIFpwdReader', 'powder_cif'),
    ('manuals/fullprof/9-SAnn-OrganicCompounds/PPH3/Rietveld/PPH3-sls_prof.cif', 'G2pwd_CIF', 'CIFpwdReader', 'powder_cif'),
    ('manuals/fullprof/9-SAnn-OrganicCompounds/CitricAcid_example/Rietveld/citric-acid-sls_prof.cif', 'G2pwd_CIF', 'CIFpwdReader', 'powder_cif'),
]

# Synthetic files (written by make_synthetic_formats.py) for the formats with no real example.
SYNTHETIC = [
    ('manuals/gsas2/synthetic/scan.xrdml', 'G2pwd_Panalytical', 'Panalytical_ReaderClass', 'panalytical_xrdml'),
    ('manuals/gsas2/synthetic/scan.brml', 'G2pwd_BrukerBRML', 'brml_ReaderClass', 'bruker_brml'),
    ('manuals/gsas2/synthetic/scan.ras', 'G2pwd_rigaku', 'Rigaku_rasReaderClass', 'rigaku_ras'),
    ('manuals/gsas2/synthetic/scan.rasx', 'G2pwd_rigaku', 'Rigaku_rasReaderClass', 'rigaku_rasx'),
    ('manuals/gsas2/synthetic/scan.txt', 'G2pwd_rigaku', 'Rigaku_txtReaderClass', 'rigaku_txt'),
    ('manuals/gsas2/synthetic/scan.csv', 'G2pwd_csv', 'csv_ReaderClass', 'table'),
    ('manuals/gsas2/synthetic/scan.chi', 'G2pwd_xye', 'xye_ReaderClass', 'fit2d_chi'),
    ('manuals/gsas2/synthetic/scan.qye', 'G2pwd_xye', 'xye_ReaderClass', 'qye'),
    ('manuals/gsas2/synthetic/scan_v4.raw', 'G2pwd_BrukerRAW', 'raw_ReaderClass', 'bruker_raw_v4'),
    ('manuals/gsas2/synthetic/scan_v2.raw', 'G2pwd_BrukerRAW', 'raw_ReaderClass', 'bruker_raw_v2'),
]


def wavelengths(rd):
    """The wavelengths GSAS-II takes from the file: instdict['wave'], the Wave of a Q file, or Ka1/Ka2/ratio comments."""
    w = rd.instdict.get('wave') if hasattr(rd, 'instdict') else None
    if w is None and getattr(rd, 'Wave', None):
        w = rd.Wave
    if w is not None:
        return ' '.join(repr(float(v)) for v in (w if isinstance(w, list) else [w]))
    found = {}
    for c in rd.comments if isinstance(rd.comments, list) else []:
        m = re.match(r'\s*(Ka1|Ka2|Ka2/Ka1)\s*=\s*([-0-9.eE+]+)', str(c))
        if m:
            found[m.group(1)] = m.group(2)
    if 'Ka1' in found:
        return ' '.join(found[k] for k in ('Ka1', 'Ka2', 'Ka2/Ka1') if k in found)
    return '-'


def main():
    g2 = os.path.abspath(sys.argv[1])
    sys.path.insert(0, g2)
    import importlib
    import types
    # GSAS-II's CIF importer asks the user, through its wx dialogs, which x, intensity and weight items to use when a
    # block has several; a stub answers with the first of each list (and every block), as a user taking the defaults.
    stub = types.ModuleType('GSASII.GSASIIctrlGUI')
    stub.MultipleChoicesSelector = lambda choices, labels, *a, **k: [0, 0, 0, 0][:len(choices)]
    stub.MultipleBlockSelector = lambda choices, *a, **k: list(range(len(choices)))
    sys.modules['GSASII.GSASIIctrlGUI'] = stub
    from GSASII import GSASIIscriptable as G2sc
    out = [
        '# What GSAS-II reads from each example file (GSASII commit fe944cd), for the test data_formats; written by',
        '# make_gsas2_reads.py. Do not edit. Tab-separated: file, importer, format id, bank (1-based), points,',
        '# sum x, sum y, sum w (w = 1/σ², 0 for no weight), first x y w, last x y w, wavelengths (or -), and for a GSAS bank',
        '# with more values than its NCHAN the first NCHAN points\' count and sums (or -); or FAILED and why.',
        '# x as GSAS-II\'s reader leaves it: 2θ in degrees, GSAS TOF divided by 100.',
    ]
    for path, mod, cls, fid in FILES + SYNTHETIC:
        full = os.path.join(ROOT, path)
        if not os.path.exists(full):
            print('missing', path, file=sys.stderr)
            continue
        m = importlib.import_module('GSASII.imports.' + mod)
        reader = getattr(m, cls)()
        cwd = os.getcwd()
        os.chdir(os.path.join(ROOT, 'manuals', 'gsas2'))  # BRML and RASX readers unzip into the working directory
        try:
            rds = G2sc.import_generic(full, [reader])
        except Exception as e:
            # A scripted import cannot answer GSAS-II's questions (which of several intensity items of a CIF).
            out.append('\t'.join([path, f'{mod}.{cls}', fid, 'FAILED', type(e).__name__ + ': ' + str(e).replace(ROOT + '/', '')]))
            continue
        finally:
            os.chdir(cwd)
        if cls == 'raw_ReaderClass':
            # G2pwd_BrukerRAW.py ends with self.repeat = False (its line 273), so the loader stops after the first
            # range of a version 2 or 3 file; its reader reads the next ranges when asked for them by number.
            b = len(rds)
            while b < 16:
                b += 1
                try:
                    ok = reader.Reader(full, buffer={}, blocknum=b)
                except Exception:
                    ok = False
                if not ok or not len(reader.powderdata[0]) or any(len(r.powderdata[0]) == len(reader.powderdata[0]) and r.powderdata[0][0] == reader.powderdata[0][0] for r in rds):
                    break
                import copy
                rds.append(copy.deepcopy(reader))
        nchan = []
        if fid == 'gsas':
            with open(full, encoding='latin-1') as fh:
                nchan = [int(l.split()[2]) for l in fh if l.startswith('BANK')]
        for k, rd in enumerate(rds):
            x, y, w = (list(map(float, a)) for a in rd.powderdata[:3])
            n = min(len(x), len(y), len(w))
            fmt = lambda i: f'{x[i]!r} {y[i]!r} {w[i]!r}'
            # A GSAS bank holding more values than its NCHAN: the sums over the first NCHAN points too.
            m = nchan[k] if k < len(nchan) else n
            prefix = f'{m} {math.fsum(x[:m])!r} {math.fsum(y[:m])!r} {math.fsum(w[:m])!r}' if m < n else '-'
            out.append('\t'.join([path, f'{mod}.{cls}', fid, str(k + 1), str(n), repr(math.fsum(x[:n])), repr(math.fsum(y[:n])),
                                  repr(math.fsum(w[:n])), fmt(0), fmt(n - 1), wavelengths(rd), prefix]))
    with open(os.path.join(ROOT, 'core', 'tests', 'data', 'gsas2_reads.txt'), 'w') as f:
        f.write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
