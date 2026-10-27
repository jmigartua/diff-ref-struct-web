#!/usr/bin/env python3
"""SciPy's fits of FullProf's own sequential results, the reference of the core's parametric fits (test
`parametric_fits`, P54, L317).

Reads manuals/fullprof/1-WinPLOTR/2_FP_sequential_refinements/LaMnO3_D20/la0.seq (FullProf's .seq of the D20 run, 71
blocks: the external parameter PARVAL, the cell volume of phase 1 and the refined Cell_C_ph1_pat1, each with its
σ) and fits, with scipy.optimize.curve_fit (sigma = σ, absolute_sigma=False: the coefficients' su's scaled by χ²_ν),
three equations of T:
    volume   V0*(1 + a*(T - 150) + b*(T - 150)^2)          start V0 245, a 1e-5, b 1e-8
    cell_c   c0 + c1*T                                     start c0 7.6, c1 1e-4
    einstein V0 + A/(exp(th/T) - 1)                        start V0 245, A 1, th 300
Writes core/tests/data/parametric_reference.txt: per fit its name, the points, χ² and χ²_ν, and each coefficient's
value and su to 17 digits.

Usage (numpy and scipy):  python3 core/tests/data/make_parametric_reference.py
"""
import os

import numpy as np
import scipy
from scipy.optimize import curve_fit

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
SEQ = os.path.join(ROOT, 'manuals/fullprof/1-WinPLOTR/2_FP_sequential_refinements/LaMnO3_D20/la0.seq')
OUT = os.path.join(ROOT, 'core/tests/data/parametric_reference.txt')


def blocks():
    out, cur = [], None
    for line in open(SEQ, encoding='latin-1'):
        t = line.split()
        if not t:
            continue
        if t[0] == 'NEW_REFINEMENT':
            cur = {'name': t[1]}
            out.append(cur)
        elif t[0] == 'PARVAL':
            cur['T'] = float(t[2])
        elif len(t) >= 4 and t[0] == '1' and t[1] == 'Unit_Cell_Volume:':
            cur['V'] = (float(t[2]), float(t[3]))
        elif len(t) == 4 and t[1] == 'Cell_C_ph1_pat1':
            cur['c'] = (float(t[2]), float(t[3]))
    return out


def main():
    b = blocks()
    # A file may repeat a data file (an appended second run); the last block of each name is kept, as FullProf's
    # later run is the one its .pcr files hold.
    seen = {}
    for x in b:
        seen[x['name']] = x
    b = sorted(seen.values(), key=lambda x: x['name'])
    T = np.array([x['T'] for x in b])
    fits = [
        ('volume', lambda T, V0, a, bb: V0 * (1 + a * (T - 150) + bb * (T - 150) ** 2), [245.0, 1e-5, 1e-8], 'V'),
        ('cell_c', lambda T, c0, c1: c0 + c1 * T, [7.6, 1e-4], 'c'),
        ('einstein', lambda T, V0, A, th: V0 + A / (np.exp(th / T) - 1), [245.0, 1.0, 300.0], 'V'),
    ]
    lines = [f'# scipy {scipy.__version__}, numpy {np.__version__}: curve_fit(sigma, absolute_sigma=False) of {os.path.relpath(SEQ, ROOT)}']
    for name, f, p0, key in fits:
        y = np.array([x[key][0] for x in b])
        s = np.array([x[key][1] for x in b])
        p, cov = curve_fit(f, T, y, p0=p0, sigma=s, absolute_sigma=False, maxfev=100000)
        r = (y - f(T, *p)) / s
        chi2 = float(np.sum(r * r))
        lines.append(f'fit {name} points {len(T)} chi2 {chi2!r} chi2_reduced {chi2 / (len(T) - len(p))!r}')
        for k, v in enumerate(p):
            lines.append(f'  {k} {float(v)!r} {float(np.sqrt(cov[k][k]))!r}')
    open(OUT, 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
