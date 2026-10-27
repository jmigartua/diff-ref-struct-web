#!/usr/bin/env python3
"""Independent references for the background points (test `background_points`, L078).

On Y2O3's 62 background points (core/tests/data/broken/base.drc, read with Python's tomllib: the control file is
TOML), writes core/tests/data/background_references.txt with:
- `spline`: scipy.interpolate.CubicSpline(x, y, bc_type="natural") at 2θ from 15.48 to 144.97 in steps of 0.37;
- `chebyshev`: numpy.polynomial.chebyshev.chebfit of the points mapped to [−1, 1] over [from, to] = [15.48, 144.97],
  with 6 and 12 terms (least squares, unweighted);
- `polynomial`: numpy.polyfit in u = 2θ/origin − 1, origin 80, 6 terms, as FullProf's eq. 3.4 writes it.

Usage (from the repository root; numpy and scipy): python3 core/tests/data/make_background_references.py
"""
import os
import tomllib

import numpy as np
import scipy
from scipy.interpolate import CubicSpline

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))


def main():
    d = tomllib.load(open(os.path.join(ROOT, 'core', 'tests', 'data', 'broken', 'base.drc'), 'rb'))
    pts = np.array(d['pattern']['background']['background_points'], dtype=float)
    x, y = pts[:, 0], pts[:, 1]
    out = [f'# Background references (make_background_references.py), numpy {np.__version__}, scipy {scipy.__version__}. Do not edit.',
           '# spline <2θ> <value>; chebyshev <terms> <from> <to> <coefficients>; polynomial <terms> <origin> <coefficients>']
    s = CubicSpline(x, y, bc_type='natural')
    t = np.arange(x[0], x[-1], 0.37)
    for v in t:
        out.append(f'spline {float(v)!r} {float(s(v))!r}')
    lo, hi = float(x[0]), float(x[-1])
    for n in (6, 12):
        c = np.polynomial.chebyshev.chebfit(2 * (x - lo) / (hi - lo) - 1, y, n - 1)
        out.append(f'chebyshev {n} {lo!r} {hi!r} ' + ' '.join(repr(float(v)) for v in c))
    origin = 80.0
    c = np.polyfit(x / origin - 1, y, 5)[::-1]
    out.append(f'polynomial 6 {origin!r} ' + ' '.join(repr(float(v)) for v in c))
    open(os.path.join(ROOT, 'core', 'tests', 'data', 'background_references.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
