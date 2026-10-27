#!/usr/bin/env python3
"""What GSAS-II reads from the core's results file and structure CIF, for the test `results_file` (R28, L101).

Runs GSAS-II's own importers on core/tests/data/results/a_modes_results.cif (the results file of `drs results
manuals/fullprof/runs/modes/A_modes/a_modes.pcr --cycles 1 --time "2026-10-11 12.00.00"`) and on a_modes_phase.cif
(the structure CIF written beside it): the powder importer (imports/G2pwd_CIF.py) and the phase importer
(imports/G2phase_CIF.py), through GSAS-II's scripting loader (GSASIIscriptable.import_generic), as
make_gsas2_reads.py runs them. Writes core/tests/data/gsas2_results.txt:
    profile  <points> <sum x> <sum y> <sum w> <first x y w> <last x y w> <wavelength>
    phase    <file> <a b c alpha beta gamma> <space group> <atoms> <sum of x+y+z over the atoms>

Usage (from the repository root; GSAS-II needs numpy, scipy and PyCifRW, and its phase importer the compiled space-group
module pyspg, which a source checkout lacks: built from its sources/pyspg.for and sources/spsubs/*.for with
`python3 -m numpy.f2py -c pyspg.for spsubs/*.for -m pyspg --lower -I sources/INCLDS` and put on PYTHONPATH):
    python3 core/tests/data/make_gsas2_results.py <GSAS-II checkout>
The checkout used for the committed table: https://github.com/AdvancedPhotonSource/GSAS-II, commit 305a159.
"""
import math
import os
import sys
import types

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
DIR = os.path.join(ROOT, 'core', 'tests', 'data', 'results')


def main():
    g2 = os.path.abspath(sys.argv[1])
    sys.path.insert(0, g2)
    # The CIF importers ask, through wx dialogs, which items and blocks to use; a stub takes the first of each list and
    # every block, as a user taking the defaults (make_gsas2_reads.py).
    stub = types.ModuleType('GSASII.GSASIIctrlGUI')
    stub.MultipleChoicesSelector = lambda choices, labels, *a, **k: [0, 0, 0, 0][:len(choices)]
    stub.MultipleBlockSelector = lambda choices, *a, **k: list(range(len(choices)))
    sys.modules['GSASII.GSASIIctrlGUI'] = stub
    # GSAS-II imports its compiled pyspg as GSASII.pyspg when no binary path is set.
    try:
        import pyspg
        sys.modules['GSASII.pyspg'] = pyspg
    except ImportError:
        pass
    from GSASII import GSASIIscriptable as G2sc
    from GSASII.imports import G2pwd_CIF, G2phase_CIF
    out = ['# What GSAS-II (commit 305a159) reads from the core\'s results file and structure CIF, for the test results_file;',
           '# written by make_gsas2_results.py. Do not edit.']
    rds = G2sc.import_generic(os.path.join(DIR, 'a_modes_results.cif'), [G2pwd_CIF.CIFpwdReader()])
    for rd in rds:
        x, y, w = (list(map(float, a)) for a in rd.powderdata[:3])
        wave = rd.instdict.get('wave')
        out.append(' '.join(['profile', str(len(x)), repr(math.fsum(x)), repr(math.fsum(y)), repr(math.fsum(w)),
                             repr(x[0]), repr(y[0]), repr(w[0]), repr(x[-1]), repr(y[-1]), repr(w[-1]), repr(wave)]))
    for f in ('a_modes_results.cif', 'a_modes_phase.cif'):
        for rd in G2sc.import_generic(os.path.join(DIR, f), [G2phase_CIF.CIFPhaseReader()]):
            gen = rd.Phase['General']
            cell = gen['Cell'][1:7]
            atoms = rd.Phase['Atoms']
            cx = gen['AtomPtrs'][0]
            s = math.fsum(a[cx] + a[cx + 1] + a[cx + 2] for a in atoms)
            out.append(' '.join(['phase', f] + [repr(float(v)) for v in cell] + [gen['SGData']['SpGrp'].replace(' ', ''), str(len(atoms)), repr(s)]))
    with open(os.path.join(ROOT, 'core', 'tests', 'data', 'gsas2_results.txt'), 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    print('\n'.join(out))


if __name__ == '__main__':
    main()
