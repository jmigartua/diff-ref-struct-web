# A refinement of diff-ref-struct's control file by GSAS-II, for a comparison (the app's request 50).
# Written by the core (compare::gsas2_prepare); run with GSAS-II's Python, the GSAS-II folder on PYTHONPATH:
#     python compare.py
# It reads compare_input.json (the data, the instrument file, each phase's CIF, the values the CIFs do not hold and
# the refinement flags, all chosen by the core), builds compare.gpx, refines the cycles asked for and writes
# compare_output.json: GSAS-II's version, every varied parameter with its esd, the cells and atoms with esd's, the
# residuals and the calculated pattern. It decides nothing of its own.
import json
import sys

import numpy as np
from GSASII import GSASIIscriptable as G2sc
from GSASII import GSASIIpath

cfg = json.load(open("compare_input.json"))
gpx = G2sc.G2Project(newgpx="compare.gpx")
hist = gpx.add_powder_histogram(cfg["data"], cfg["instprm"], fmthint="Topas xye")
sample = hist.SampleParameters
sample["Type"] = cfg["geometry"]
sample["Gonio. radius"] = cfg["goniometer_radius"]
for key, value in cfg["sample_values"].items():
    sample[key][0] = value
if cfg["excluded"]:
    hist.Excluded(cfg["excluded"])
phases = []
for p in cfg["phases"]:
    ph = gpx.add_phase(p["cif"], phasename=p["name"], histograms=[hist], fmthint="CIF")
    # No sample broadening of GSAS-II's own (its defaults: 1 µm, 1000 microstrain): the core's widths are the
    # instrument's U, V, W, X, Y.
    hap = ph.data["Histograms"][hist.name]
    hap["Size"][1][0] = 1.0e5
    hap["Mustrain"][1][0] = 0.0
    phases.append(ph)
# The background: GSAS-II's own function, refined (the core's is not carried).
b = cfg["background"]
hist.set_refinements({"Background": {"type": b["type"], "no. coeffs": b["terms"], "refine": True}})
hist.set_refinements({"Instrument Parameters": cfg["instrument_refine"], "Sample Parameters": cfg["sample_refine"]})
# The histogram's scale stays fixed: each phase's fraction carries the core's scale (both refined are one value).
hist.clear_refinements({"Sample Parameters": ["Scale"]})
for ph, p in zip(phases, cfg["phases"]):
    atoms = {a["label"]: a["refine"] for a in p["atoms"] if a["refine"]}
    ph.set_refinements({"Cell": p["refine_cell"], "Atoms": atoms})
    ph.set_HAP_refinements({"Scale": p["refine_scale"]}, histograms=[hist])
gpx.set_Controls("cycles", cfg["cycles"])
gpx.refine()

cov = gpx["Covariance"]["data"]
rvals = cov.get("Rvals", {})
vary = list(cov.get("varyList", []))
values = list(cov.get("variables", []))
sig = list(cov.get("sig", []))


def number(v):
    """A float, or None where GSAS-II has none (JSON holds no NaN)."""
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    return v if np.isfinite(v) else None


out = {
    "version": str(GSASIIpath.GetVersionNumber()),
    "python": sys.version.split()[0],
    "converged": bool(rvals.get("converged", False)),
    "chisq": number(rvals.get("chisq")),
    "gof": number(rvals.get("GOF")),
    "nobs": int(rvals.get("Nobs", 0)),
    "nvars": int(rvals.get("Nvars", 0)),
    "residuals": {k: number(v) for k, v in hist.residuals.items()},
    "varied": {n: [number(v), number(e)] for n, v, e in zip(vary, values, sig)},
    "phases": [],
    "x": [number(v) for v in hist.getdata("x")],
    "yobs": [number(v) for v in hist.getdata("yobs")],
    "ycalc": [number(v) for v in hist.getdata("ycalc")],
}
for ph in phases:
    cell, cellsig = ph.get_cell_and_esd()[:2]
    atoms = []
    for a in ph.atoms():
        atoms.append({"label": a.label, "xyz": [number(v) for v in a.coordinates], "uiso": number(a.uiso), "frac": number(a.occupancy)})
    out["phases"].append({"id": ph.id, "name": ph.name, "cell": {k: number(v) for k, v in cell.items()}, "cell_esd": {k: number(v) for k, v in cellsig.items()}, "atoms": atoms})
json.dump(out, open("compare_output.json", "w"), allow_nan=False)
print("compare_output.json written")
