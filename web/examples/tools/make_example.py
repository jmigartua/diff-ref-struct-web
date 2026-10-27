#!/usr/bin/env python3
"""Write the example data set of diff-ref-struct Web: a synthetic in-situ heating run, 300 K to 1300 K.

One FullProf control file (.pcr) per temperature is written into a temporary folder and `drs simulate` turns each into
a simulated observed pattern (the calculated pattern with Poisson noise, seed = the scan's number). The patterns and a
manifest go to app/ui/host-web/examples/ramp/ (the page's own files, served beside it). Re-running the script writes the same files byte for byte (same drs version).

The two phases (laboratory X-rays, Cu Ka1, 2theta 15-75 deg):
  periclase MgO, F m -3 m, a = 4.2112 A at 300 K, linear expansion 1.3e-5 /K;
  corundum alpha-Al2O3, R -3 c (hexagonal axes), a = 4.7589 A, c = 12.991 A at 300 K, expansion 8e-6 /K on a and c.
Corundum's weight fraction is zero up to 1000 K and grows along a logistic centred at 1100 K to one half at 1300 K;
periclase's falls in step. No real reaction is modelled: the numbers are made up to give a plausible-looking series.

Usage: python3 web/examples/tools/make_example.py [--pcr-dir <folder>]
  DRS=<path to drs> overrides the default binary; --pcr-dir keeps a copy of the .pcr files there.
Python 3 standard library only.
"""

import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

DRS_DEFAULT = "/Users/User/Desktop/2026/20261001/20261009_diff-ref-str/core/target/release/drs"
HERE = Path(__file__).resolve().parent
OUT = HERE.parents[2] / "app" / "ui" / "host-web" / "examples" / "ramp"   # served beside the page (Try the example)

# --- the run ---------------------------------------------------------------------------------------------------
TEMPERATURES = [300 + 50 * k for k in range(21)]  # K
LAMBDA = 1.5406  # Cu Ka1, A
THMIN, STEP, THMAX = 15.0, 0.05, 75.0

# --- the instrument (pseudo-Voigt, Npr = 5) ---------------------------------------------------------------------
U, V, W = 0.020, -0.010, 0.020  # FWHM^2 = U tan^2 + V tan + W (deg^2): about 0.14 deg at 15-45 deg
ETA0 = 0.50  # Lorentzian fraction (Shape1)
BKPOS = 45.0
BACKGROUND = [150.0, -20.0, 30.0, 0.0, 0.0, 0.0]  # polynomial in (2theta/Bkpos - 1): about 176 to 150 counts

# --- the phases -------------------------------------------------------------------------------------------------
MGO_A300, MGO_ALPHA = 4.2112, 1.3e-5
COR_A300, COR_C300, COR_ALPHA = 4.7589, 12.991, 8.0e-6
MGO_ZM = 4 * (24.305 + 15.999)  # Z * formula mass: 161.216
COR_ZM = 6 * (2 * 26.982 + 3 * 15.999)  # 611.778
S0 = 1.4e-2  # periclase's scale when it is the whole sample (sets the counts: a few thousand at the strongest peak)
T_MID, T_WIDTH, T_ONSET, T_END, W_END = 1100.0, 40.0, 1000.0, 1300.0, 0.5


def mgo_a(t):
    return MGO_A300 * (1 + MGO_ALPHA * (t - 300))


def cor_ac(t):
    f = 1 + COR_ALPHA * (t - 300)
    return COR_A300 * f, COR_C300 * f


def b_iso(b300, t):
    """Isotropic displacement, rising linearly by 0.15 A^2 from 300 K to 1300 K."""
    return b300 + 0.15 * (t - 300) / 1000


def corundum_fraction(t):
    """Weight fraction of corundum: 0 up to T_ONSET, a logistic about T_MID rescaled to W_END at T_END."""
    if t <= T_ONSET:
        return 0.0
    lg = lambda x: 1 / (1 + math.exp(-(x - T_MID) / T_WIDTH))
    g = (lg(t) - lg(T_ONSET)) / (lg(T_END) - lg(T_ONSET))
    return W_END * min(1.0, max(0.0, g))


def scales(t):
    """Scale factors giving the weight fractions by Hill & Howard's rule, W_i proportional to S_i Z M V."""
    w = corundum_fraction(t)
    v_mgo = mgo_a(t) ** 3
    a, c = cor_ac(t)
    v_cor = math.sqrt(3) / 2 * a * a * c
    ref = MGO_ZM * MGO_A300**3
    return S0 * (1 - w) * ref / (MGO_ZM * v_mgo), S0 * w * ref / (COR_ZM * v_cor), w


def pcr(n, t):
    s_mgo, s_cor, w = scales(t)
    a = mgo_a(t)
    ca, cc = cor_ac(t)
    zeros5 = "                  0.00      0.00      0.00      0.00      0.00"
    bg = "".join(f"{b:12.3f}" for b in BACKGROUND)
    return f"""COMM Synthetic heating run, scan {n:02d}, {t} K: periclase + corundum (diff-ref-struct Web example)
! Written by web/examples/tools/make_example.py; the counts are simulated, not measured.
! Corundum weight fraction {w:.4f}
!Job Npr Nph Nba Nex Nsc Nor Dum Iwg Ilo Ias Res Ste Nre Cry Uni Cor Opt Aut
   0   5   2   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
!
!Ipr Ppl Ioc Mat Pcr Ls1 Ls2 Ls3 NLI Prf Ins Rpa Sym Hkl Fou Sho Ana
   0   0   1   0   1   0   0   0   0   0  10   0   0   0   0   0   0
!
! Lambda1  Lambda2    Ratio    Bkpos    Wdt    Cthm     muR   AsyLim   Rpolarz  2nd-muR -> Patt# 1
 {LAMBDA:8.6f} {LAMBDA:8.6f}  0.00000 {BKPOS:8.3f}   8.000  0.0000  0.0000   45.00    0.0000  0.0000
!
!NCY  Eps  R_at  R_an  R_pr  R_gl     Thmin       Step       Thmax    PSD    Sent0
  1  0.10  1.00  1.00  1.00  1.00 {THMIN:11.4f} {STEP:10.6f} {THMAX:10.4f}   0.000   0.000
!
       0    !Number of refined parameters
!
!  Zero    Code    SyCos    Code   SySin    Code  Lambda     Code MORE ->Patt# 1
  0.00000    0.0  0.00000    0.0  0.00000    0.0 0.000000    0.00   0
!   Background coefficients/codes  for Pattern#  1  (Polynomial of 6th degree)
{bg}
        0.00        0.00        0.00        0.00        0.00        0.00
!-------------------------------------------------------------------------------
!  Data for PHASE number:   1
!-------------------------------------------------------------------------------
Periclase MgO
!
!Nat Dis Ang Pr1 Pr2 Pr3 Jbt Irf Isy Str Furth       ATZ    Nvk Npr More
   2   0   0 0.0 0.0 1.0   0   0   0   0   0 {MGO_ZM:14.3f}   0   5   0
!
F m -3 m                 <--Space group symbol
!Atom   Typ       X        Y        Z     Biso       Occ     In Fin N_t Spc /Codes
Mg1    Mg      0.00000  0.00000  0.00000 {b_iso(0.30, t):8.5f}  0.02083   0   0   0    0
{zeros5}
O1     O       0.50000  0.50000  0.50000 {b_iso(0.35, t):8.5f}  0.02083   0   0   0    0
{zeros5}
!-------> Profile Parameters for Pattern #   1  ----> Phase #   1
!  Scale          Shape1      Bov      Str1      Str2      Str3   Strain-Model
 {s_mgo:.7E}   {ETA0:7.5f}   0.00000   0.00000   0.00000   0.00000       0
       0.00000     0.000     0.000     0.000     0.000     0.000
!       U            V            W             X            Y         GauSiz      LorSiz    Size-Model
 {U:12.6f} {V:12.6f} {W:12.6f}     0.000000     0.000000     0.000000     0.000000       0
         0.00         0.00         0.00         0.00         0.00         0.00         0.00
!     a          b         c        alpha      beta       gamma      # Cell Info
 {a:10.6f} {a:10.6f} {a:10.6f}  90.000000  90.000000  90.000000
    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000
!  Pref1    Pref2      Asy1     Asy2     Asy3     Asy4
  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000
     0.00     0.00     0.00     0.00     0.00     0.00
!-------------------------------------------------------------------------------
!  Data for PHASE number:   2
!-------------------------------------------------------------------------------
Corundum Al2O3
!
!Nat Dis Ang Pr1 Pr2 Pr3 Jbt Irf Isy Str Furth       ATZ    Nvk Npr More
   2   0   0 0.0 0.0 1.0   0   0   0   0   0 {COR_ZM:14.3f}   0   5   0
!
R -3 c                   <--Space group symbol
!Atom   Typ       X        Y        Z     Biso       Occ     In Fin N_t Spc /Codes
Al1    Al      0.00000  0.00000  0.35230 {b_iso(0.30, t):8.5f}  0.33333   0   0   0    0
{zeros5}
O1     O       0.30640  0.00000  0.25000 {b_iso(0.35, t):8.5f}  0.50000   0   0   0    0
{zeros5}
!-------> Profile Parameters for Pattern #   1  ----> Phase #   2
!  Scale          Shape1      Bov      Str1      Str2      Str3   Strain-Model
 {s_cor:.7E}   {ETA0:7.5f}   0.00000   0.00000   0.00000   0.00000       0
       0.00000     0.000     0.000     0.000     0.000     0.000
!       U            V            W             X            Y         GauSiz      LorSiz    Size-Model
 {U:12.6f} {V:12.6f} {W:12.6f}     0.000000     0.000000     0.000000     0.000000       0
         0.00         0.00         0.00         0.00         0.00         0.00         0.00
!     a          b         c        alpha      beta       gamma      # Cell Info
 {ca:10.6f} {ca:10.6f} {cc:10.6f}  90.000000  90.000000 120.000000
    0.00000    0.00000    0.00000    0.00000    0.00000    0.00000
!  Pref1    Pref2      Asy1     Asy2     Asy3     Asy4
  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000
     0.00     0.00     0.00     0.00     0.00     0.00
!  2Th1/TOF1    2Th2/TOF2  Pattern to plot
 {THMIN:11.3f} {THMAX:11.3f}       1
"""


def name(n, t):
    return f"example_{n:02d}_{t}K"


def main():
    drs = os.environ.get("DRS", DRS_DEFAULT)
    keep = None
    if len(sys.argv) == 3 and sys.argv[1] == "--pcr-dir":
        keep = Path(sys.argv[2])
    elif len(sys.argv) != 1:
        sys.exit(__doc__)
    version = subprocess.run([drs, "--version"], capture_output=True, text=True, check=True).stdout.strip()
    OUT.mkdir(parents=True, exist_ok=True)
    for old in OUT.glob("example_*.xys"):
        old.unlink()
    files = []
    with tempfile.TemporaryDirectory() as tmp:
        for n, t in enumerate(TEMPERATURES, start=1):
            stem = name(n, t)
            Path(tmp, stem + ".pcr").write_text(pcr(n, t), encoding="ascii", newline="\n")
            # run in the temporary folder with a relative name: the file's title names the .pcr, never the folder
            r = subprocess.run([drs, "simulate", stem + ".pcr", "--seed", str(n), "--out", stem + ".xys"], cwd=tmp, capture_output=True, text=True)
            if r.returncode != 0:
                sys.exit(f"drs simulate {stem}.pcr failed:\n{r.stderr}{r.stdout}")
            shutil.copyfile(Path(tmp, stem + ".xys"), OUT / (stem + ".xys"))
            files.append({"name": stem + ".xys", "size": (OUT / (stem + ".xys")).stat().st_size, "temperature_K": t, "seed": n})
        if keep:
            keep.mkdir(parents=True, exist_ok=True)
            for p in sorted(Path(tmp).glob("*.pcr")):
                shutil.copyfile(p, keep / p.name)
    manifest = {
        "title": "Synthetic in-situ heating run: periclase MgO with corundum growing above 1000 K, 300 to 1300 K",
        "licence": "CC BY 4.0",
        "note": "The counts are simulated, not measured.",
        "files": files,
        "generated_by": version,
        "seed_rule": "drs simulate --seed n for scan n (1 to 21, in temperature order)",
        "command": "python3 web/examples/tools/make_example.py",
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    total = sum(f["size"] for f in files)
    print(f"{len(files)} files, {total} bytes, {version}, in {OUT}")


if __name__ == "__main__":
    main()
