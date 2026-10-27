#!/usr/bin/env python3
"""Independent references for the data operations, for the test `data_operations` (L066).

On real files, read by GSAS-II's own importers (the reads the test `data_formats` checks against the core's readers),
this writes core/tests/data/data_ops_references.txt with:
- `savgol`: scipy.signal.savgol_filter(y, points, order, mode="interp") at sampled points, ends included, and the σ
  of the smoothed interior points from scipy.signal.savgol_coeffs, √(Σ c_k² σ_k²);
- `peaks`: scipy.signal.find_peaks(y, prominence=P) and peak_widths(y, peaks, rel_height=0.5): each peak's index,
  prominence and width in x (the interpolated positions mapped to x with numpy.interp);
- `arpls`: pybaselines.whittaker.arpls(y, lam=int(10**logLam), max_iter=10), as GSAS-II's autoBkgCalc calls it
  (GSASIIpwd.py), logLam GSAS-II's default for N points (GSASIIpwdGUI.py, min(10, (int(10 log10(N)^1.5) - 9.5)/10)):
  the iterations and the baseline at sampled points;
- `gsas2_candidates`: the peak candidates of GSAS-II's automatic peak search (GSASIIpwdGUI.OnAutoSearch, its lines
  805-818 run as written, with no background file): local maxima at least the standard deviation of y;
- `rebin`: numpy sums of each n consecutive points (reshape), at sampled bins.

Usage (from the repository root; needs numpy, scipy, pybaselines and a GSAS-II checkout):
    python3 core/tests/data/make_data_ops_references.py <GSAS-II checkout>
Committed table: GSAS-II commit fe944cd, numpy and scipy of /tmp venv (versions written in the header),
pybaselines 1.2.1.
"""
import importlib
import math
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))

# (file, GSAS-II importer module, class, the registry's format id, prominence for find_peaks)
FILES = [
    ('manuals/gsas2/tutorials/AutoBkg/data/test_RampDown-00000.xye', 'G2pwd_xye', 'xye_ReaderClass', 'xye', 200.0),
    ('manuals/gsas2/testinp/PBSO4.XRA', 'G2pwd_fxye', 'GSAS_ReaderClass', 'gsas', 100.0),
    ('manuals/gsas2/tutorials/CWNeutron/data/garnet.raw', 'G2pwd_fxye', 'GSAS_ReaderClass', 'gsas', 100.0),
    ('manuals/fullprof/runs/Y2O3/y2o3.dat', 'G2pwd_FP', 'fp_ReaderClass', 'free_steps', 50.0),
]
SMOOTHING = [(7, 2), (15, 4), (5, 1)]


def sample(n, every=40, ends=8):
    """Indices: the first and last `ends`, and every n/every in between."""
    s = set(range(min(ends, n))) | set(range(max(0, n - ends), n)) | set(range(0, n, max(1, n // every)))
    return sorted(s)


def main():
    g2 = os.path.abspath(sys.argv[1])
    sys.path.insert(0, g2)
    import numpy as np
    import numpy.ma as ma
    import scipy
    import scipy.signal as ss
    import pybaselines
    import pybaselines.whittaker
    from GSASII import GSASIIscriptable as G2sc
    out = [
        '# Independent references for the data operations (test data_operations), written by make_data_ops_references.py.',
        f'# Do not edit. numpy {np.__version__}, scipy {scipy.__version__}, pybaselines {pybaselines.__version__}, GSAS-II fe944cd.',
        '# Tab-separated. file <path> <format id> <points> <sum y>; then lines for that file:',
        '# savgol <points> <order> <i> <smoothed y> <σ or - at the ends>; peaks <prominence> <count>;',
        '# peak <i> <prominence> <width in x>; arpls <logLam> <lam> <iterations>; base <i> <baseline>;',
        '# gsas2_candidates <count> <indices, space-separated>; rebin <factor> <i> <sum>.',
    ]
    for path, mod, cls, fid, prom in FILES:
        full = os.path.join(ROOT, path)
        m = importlib.import_module('GSASII.imports.' + mod)
        rd = G2sc.import_generic(full, [getattr(m, cls)()])[0]
        x = np.array(rd.powderdata[0], dtype=float)
        y = np.array(rd.powderdata[1], dtype=float)
        w = np.array(rd.powderdata[2], dtype=float)
        n = min(len(x), len(y), len(w))
        x, y, w = x[:n], y[:n], w[:n]
        sig = np.where(w > 0, 1.0 / np.sqrt(np.where(w > 0, w, 1.0)), np.inf)
        out.append('\t'.join(['file', path, fid, str(n), repr(math.fsum(y))]))
        for pts, order in SMOOTHING:
            sm = ss.savgol_filter(y, pts, order, mode='interp')
            c = ss.savgol_coeffs(pts, order, use='dot')
            h = pts // 2
            for i in sample(n):
                s = '-'
                if h <= i < n - h:
                    # use='dot': smoothed y_i = Σ_k c_k y_(i-h+k).
                    s = repr(math.sqrt(math.fsum((c[k] * sig[i - h + k]) ** 2 for k in range(pts))))
                out.append('\t'.join(['savgol', str(pts), str(order), str(i), repr(float(sm[i])), s]))
        peaks, props = ss.find_peaks(y, prominence=prom)
        wl, _, lips, rips = ss.peak_widths(y, peaks, rel_height=0.5, prominence_data=(props['prominences'], props['left_bases'], props['right_bases']))
        idx = np.arange(n, dtype=float)
        out.append('\t'.join(['peaks', repr(prom), str(len(peaks))]))
        for k, p in enumerate(peaks):
            width = float(np.interp(rips[k], idx, x) - np.interp(lips[k], idx, x))
            out.append('\t'.join(['peak', str(int(p)), repr(float(props['prominences'][k])), repr(width)]))
        log_lam = min(10, float(int(10 * np.log10(n) ** 1.5) - 9.5) / 10.)
        lam = int(10 ** log_lam)
        base, params = pybaselines.whittaker.arpls(y, lam=lam, max_iter=10)
        out.append('\t'.join(['arpls', repr(log_lam), str(lam), str(len(params['tol_history']))]))
        for i in sample(n, every=50):
            out.append('\t'.join(['base', str(i), repr(float(base[i]))]))
        # GSASIIpwdGUI.OnAutoSearch, lines 805-818, the whole pattern, no background file (bxye = 0).
        y0 = y - 0.0
        ysig = 1.0 * np.std(y0)
        offset = [-1, 1]
        ymask = ma.array(y0, mask=(y0 < ysig))
        for off in offset:
            ymask = ma.array(ymask, mask=(ymask - np.roll(y0, off) <= 0.))
        indx = ymask.nonzero()
        cand = [int(i) for i in indx[0]]
        out.append('\t'.join(['gsas2_candidates', str(len(cand)), ' '.join(map(str, cand))]))
        for f in (2, 3):
            r = y[: n // f * f].reshape(-1, f).sum(axis=1)
            for i in sample(len(r), every=20, ends=3):
                out.append('\t'.join(['rebin', str(f), str(i), repr(float(r[i]))]))
    with open(os.path.join(ROOT, 'core', 'tests', 'data', 'data_ops_references.txt'), 'w') as fh:
        fh.write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
