#!/usr/bin/env python3
"""The broken files of the test `check_diagnostics` (R24, L076): each one a run with one mistake, that the core's checks
must report as exactly one diagnostic, the one its name gives (`<code>.drc`, `<code>--<variant>.drc`, or `.pcr`).

All but one are derived from the Y2O3 run (manuals/fullprof/runs/Y2O3/y2o3_structure.pcr), which the checks pass with nothing
to report: its control file (`base.drc`, by `drs convert`) and its export (`base.pcr`). Each `.drc` is exported to the
`.pcr` next to it where the export can be made, for FullProf to say what it makes of the same mistake (run on the Mac;
FullProf finds the data, y2o3.dat, in manuals/fullprof/runs/Y2O3/). `singular.drc` comes from Yb2ZrWO8 (`base-tch.drc`),
whose profile is Thompson–Cox–Hastings.

Usage (from the repository root, after `cargo build --release` in core/):
    python3 core/tests/data/make_broken.py
"""
import os
import re
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
OUT = os.path.join(ROOT, 'core', 'tests', 'data', 'broken')
DRS = os.path.join(ROOT, 'core', 'target', 'release', 'drs')
RUN = os.path.join(ROOT, 'manuals', 'fullprof', 'runs', 'Y2O3', 'y2o3_structure.pcr')
# A run with a Thompson–Cox–Hastings profile (I_G exists only there), which the checks also pass clean.
RUN_TCH = os.path.join(ROOT, 'manuals', 'fullprof', 'runs', 'Yb2ZrWO8', 'ybfluo_Riet_IRF.pcr')


def drs(*args):
    r = subprocess.run([DRS, *args], capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def sub(text, old, new, count=1):
    """Replace exactly `count` occurrences of `old` (a regular expression)."""
    out, n = re.subn(old, new, text, count=count, flags=re.M)
    if n != count:
        sys.exit(f'pattern not found: {old!r}')
    return out


def add_parameter(text, entry):
    """A parameter added at the end of [refine] parameters."""
    return sub(text, r'^(parameters = \[\n(?:  .*\n)*?)\]', lambda m: m.group(1) + f'  {entry},\n]')


ATOM_GD = '''
[[phase.atom]]
label = "Gd"
species = "Gd"
atom_position = [0.25, 0.25, 0.25]      # fraction of the cell edge
atom_biso = 0.23928                     # Å²
atom_occupancy = 0.16666
'''

# (file name, what is broken, the edit of base.drc)
DRC = [
    ('syntax.drc', 'a number written with two decimal points',
     lambda t: sub(t, r'^zero_shift = 0\.0025', 'zero_shift = 0.0025.1')),
    ('unknown-name.drc', 'a key misspelt (zero_shft)',
     lambda t: sub(t, r'^zero_shift =', 'zero_shft =')),
    ('unit.drc', 'a value commented in another unit than the dictionary\'s (Å for a zero shift in deg 2θ)',
     lambda t: sub(t, r'^(zero_shift = 0\.0025 +)# deg 2θ', r'\1# Å')),
    ('range.drc', 'a negative scale factor',
     lambda t: sub(t, r'^scale_factor = 0\.', 'scale_factor = -0.')),
    ('data-missing.drc', 'the data file misnamed',
     lambda t: sub(t, r'y2o3\.dat"', 'y2o3.dta"')),
    ('wavelength.drc', 'a wavelength of zero',
     lambda t: sub(t, r'^wavelength_1 = 1\.5406', 'wavelength_1 = 0.0')),
    ('range-data.drc', 'Thmin and Thmax of another pattern (25–95° for data from 15° to 145.9°)',
     lambda t: sub(t, r'^range = \[15\.0, 0\.01, 145\.9\]', 'range = [25.0, 0.01, 95.9]')),
    ('occupancy.drc', 'a second atom (Gd) on the full site of Y at 1/4 1/4 1/4',
     lambda t: sub(t, r'^(atom_occupancy = 0\.16666\n)', lambda m: m.group(1) + ATOM_GD)),
    ('symmetry-fixed.drc', 'the x of Y at 1/4 1/4 1/4 (fixed by its site, 8b of I a -3) set free',
     lambda t: add_parameter(t, '["phase1.Y.atom_position.x", 1.0]')),
    ('no-reflections.drc', 'a cell of 0.5 Å, whose reflections all fall beyond 2θ = 180°',
     lambda t: sub(t, r'^(cell_[abc]) = 10\.604939', r'\1 = 0.5', count=3)),
    ('refinement.drc', 'the Caglioti U given to two parameters',
     lambda t: add_parameter(t, '["phase1.caglioti_u", 1.0]')),
    ('no-effect.drc', 'G2 refined with G1 = 1 (March–Dollase: P = 1 whatever G2)',
     lambda t: add_parameter(t, '["phase1.po_g2", 1.0]')),
    ('singular--overall-b.drc', 'the overall B refined with the B of every atom',
     lambda t: add_parameter(add_parameter(add_parameter(add_parameter(t, '["phase1.overall_b", 1.0]'),
                                                          '["phase1.Y.atom_biso", 1.0]'), '["phase1.Y:2.atom_biso", 1.0]'),
                             '["phase1.O.atom_biso", 1.0]')),
    ('correlation.drc', 'the sample displacement refined with the zero shift',
     lambda t: add_parameter(t, '["pattern1.sample_displacement", 1.0]')),
]


def flags_line(text):
    """The index of the line after '!Job Npr …' and its fields."""
    lines = text.split('\n')
    i = next(k for k, l in enumerate(lines) if l.startswith('!Job')) + 1
    return lines, i, lines[i].split()


def pcr_aut0_hole(text):
    """A hole in the parameter numbers: the codeword of parameter 3 renumbered N + 1, N + 1 declared (Aut = 0)."""
    n = int(re.search(r'^ *(\d+) +!Number of refined parameters', text, re.M).group(1))
    text = re.sub(r'^ *\d+( +!Number of refined parameters)', f'{n + 1:>6}\\1', text, count=1, flags=re.M)
    # The codeword 31.0 (parameter 3, multiplier 1) wherever it is written.
    out, k = re.subn(r'(?<![\d.])31\.0(?![\d])', f'{10 * (n + 1) + 1}.0', text)
    if k == 0:
        sys.exit('codeword 31.0 not found')
    return out


def pcr_count(text):
    """The declared count two above the numbers used (Aut = 0)."""
    n = int(re.search(r'^ *(\d+) +!Number of refined parameters', text, re.M).group(1))
    return re.sub(r'^ *\d+( +!Number of refined parameters)', f'{n + 2:>6}\\1', text, count=1, flags=re.M)


def pcr_irf(text):
    """Res = 1 with a resolution file that does not exist."""
    lines, i, f = flags_line(text)
    f[11] = '1'
    lines[i] = ''.join(f'{v:>4}' for v in f)
    lines.insert(i + 1, '!  Resolution file for Pattern#   1')
    lines.insert(i + 2, 'nothere.irf')
    return '\n'.join(lines)


PCR = [
    ('parameter-hole.pcr', 'no codeword uses parameter 3 (Aut = 0)', pcr_aut0_hole),
    ('parameter-count.pcr', 'two more parameters declared than the codewords number (Aut = 0)', pcr_count),
    ('irf-missing.pcr', 'a resolution file that is not there', pcr_irf),
]


def main():
    os.makedirs(OUT, exist_ok=True)
    for f in os.listdir(OUT):
        if f.endswith(('.drc', '.pcr', '.txt', '.irf')):
            os.remove(os.path.join(OUT, f))
    base = os.path.join(OUT, 'base.drc')
    code, msg = drs('convert', RUN, '--to', base)
    if code:
        sys.exit(msg)
    code, msg = drs('convert', base, '--to', os.path.join(OUT, 'base.pcr'))
    if code:
        sys.exit(msg)
    text = open(base, encoding='utf-8').read()
    index = ['# The broken files (make_broken.py): file, what is broken, and the .pcr export next to it (or why there is none).']
    for name, what, edit in DRC:
        path = os.path.join(OUT, name)
        open(path, 'w', encoding='utf-8').write(edit(text))
        pcr = path[:-4] + '.pcr'
        code, msg = drs('convert', path, '--to', pcr)
        index.append(f'{name}\t{what}\t' + (os.path.basename(pcr) if code == 0 else 'no export: ' + msg.strip().splitlines()[-1]))
    # From the TCH run: I_G refined with U and W.
    tch = os.path.join(OUT, 'base-tch.drc')
    code, msg = drs('convert', RUN_TCH, '--to', tch)
    if code:
        sys.exit(msg)
    path = os.path.join(OUT, 'singular.drc')
    t = open(tch, encoding='utf-8').read()
    open(path, 'w', encoding='utf-8').write(add_parameter(add_parameter(t, '["phase1.caglioti_u", 1.0]'), '["phase1.caglioti_w", 1.0]'))
    code, msg = drs('convert', path, '--to', path[:-4] + '.pcr')
    index.append('singular.drc\tthe Gaussian size I_G of Yb2ZrWO8 (TCH) refined with U and W (1/cos²θ = tan²θ + 1)\t' + ('singular.pcr' if code == 0 else 'no export: ' + msg.strip().splitlines()[-1]))
    pcr_base = open(os.path.join(OUT, 'base.pcr'), encoding='latin-1').read()
    for name, what, edit in PCR:
        open(os.path.join(OUT, name), 'w', encoding='latin-1').write(edit(pcr_base))
        index.append(f'{name}\t{what}\t(itself a .pcr)')
    open(os.path.join(OUT, 'index.txt'), 'w', encoding='utf-8').write('\n'.join(index) + '\n')
    print('\n'.join(index))


if __name__ == '__main__':
    main()
