#!/usr/bin/env python3
"""Write core/src/scattering_table.rs, the scattering tables of CrysFML, the crystallographic library of FullProf
(J. Rodríguez-Carvajal, ILL; LGPL-3), from Src/CFML_Scattering_Chemical_Tables.f90:

- `Xray_Form` (214 entries, atoms and ions): f(s) = Σ a_i exp(−b_i s²) + c, s = sinθ/λ, in electrons;
- `Anomalous_ScFac` (98 elements): f′ and f″ at the Kα1 lines of Cr, Fe, Cu, Mo, Ag, in electrons;
- `Chem_Info` (108 entries): `SctF`, the neutron (Fermi) coherent scattering length, in 10⁻¹² cm, and `AtWe`, the
  atomic weight (the centre of mass of a TLS rigid group, P41);
- `Xray_Wavelengths` (7 anodes): Kα1, Kα2, Kβ in Å.

    python3 make_scattering_tables.py <path to CFML_Scattering_Chemical_Tables.f90> <crysfml commit> > ../../src/scattering_table.rs

The numbers are copied as written in the Fortran source; nothing is typed by hand."""
import re, sys

path, commit = sys.argv[1], sys.argv[2]
text = open(path, encoding="latin-1").read()
# Join Fortran continuation lines ("&" at the end of a line, optionally "&" at the start of the next).
text = re.sub(r"&[ \t]*(![^\n]*)?\n[ \t]*&?", " ", text)
num = r"\s*([-+0-9.EeDd]+)\s*"
arr5 = r"\(/" + ",".join([num] * 5) + r"/\)"
arr4 = r"\(/" + ",".join([num] * 4) + r"/\)"
xray = re.findall(r'xray_form_type\(\s*"([^"]*)"\s*,\s*(\d+)\s*,\s*' + arr4 + r"\s*,\s*" + arr4 + r"\s*," + num + r"\)", text)
anom = re.findall(r'Anomalous_Sc_Type\(\s*"([^"]*)"\s*,\s*' + arr5 + r"\s*,\s*" + arr5 + r"\s*\)", text)
chem = re.findall(r'chem_info_type\(\s*"([^"]*)"\s*,\s*"([^"]*)"\s*,\s*(\d+)\s*,' + ",".join([num] * 4) + r",\s*\(/[^/]*/\)\s*,\s*\(/[^/]*/\)\s*," + ",".join([num] * 3) + r"\)", text)
waves = re.findall(r'Xray_Wavelength_type\(\s*"([^"]*)"\s*,\s*\(/' + num + "," + num + r"/\)\s*," + num + r"\)", text)
assert len(xray) == 214 and len(anom) == 98 and len(chem) == 108 and len(waves) == 7, (len(xray), len(anom), len(chem), len(waves))
def f(s):
    """A Fortran real as a Rust f64 literal."""
    s = s.replace("D", "E").replace("d", "e")
    return s if "." in s or "e" in s.lower() else s + ".0"

print(f"""//! Scattering tables of CrysFML, the crystallographic library of FullProf: X-ray form factors (four Gaussians and a
//! constant), the anomalous terms f′ and f″ at five anodes, neutron coherent scattering lengths, and the anode lines.
//!
//! Source: CrysFML (J. Rodríguez-Carvajal, ILL; LGPL-3), `Src/CFML_Scattering_Chemical_Tables.f90`, tables
//! `Xray_Form`, `Anomalous_ScFac`, `Chem_Info` (fields `SctF` and `AtWe`) and `Xray_Wavelengths`, from the copy at
//! https://github.com/easyscience/crysfml, commit {commit[:7]}; written by `core/tests/data/make_scattering_tables.py`.
//! Do not edit. The test `fullprof_structure_factors` checks the entries FullProf printed in its `.out` files.
// Table values such as 0.318 are data, not approximations of mathematical constants.
#![allow(clippy::approx_constant)]

/// An X-ray form factor: (symbol in CrysFML's lower case, e.g. "fe+3"; Z; a1..a4; b1..b4; c).
pub type XrayForm = (&'static str, u8, [f64; 4], [f64; 4], f64);

/// X-ray form factors, f(s) = Σ a_i exp(−b_i s²) + c, s = sinθ/λ in Å⁻¹, f in electrons.
pub static XRAY_FORM: [XrayForm; {len(xray)}] = [""")
for e in xray:
    a = ", ".join(f(x) for x in e[2:6]); b = ", ".join(f(x) for x in e[6:10])
    print(f'    ("{e[0].strip()}", {e[1]}, [{a}], [{b}], {f(e[10])}),')
print(f"""];

/// Anomalous scattering: (element in lower case; f′ at Cr, Fe, Cu, Mo, Ag Kα1; f″ at the same), in electrons.
pub static ANOMALOUS: [(&str, [f64; 5], [f64; 5]); {len(anom)}] = [""")
for e in anom:
    print(f'    ("{e[0].strip()}", [{", ".join(f(x) for x in e[1:6])}], [{", ".join(f(x) for x in e[6:11])}]),')
print(f"""];

/// Neutron coherent scattering lengths: (symbol in CrysFML's upper case, e.g. "FE"; Z; b in 10⁻¹² cm).
pub static NEUTRON_B: [(&str, u8, f64); {len(chem)}] = [""")
for e in chem:
    print(f'    ("{e[0].strip()}", {e[2]}, {f(e[7])}),')
print(f"""];

/// Anode lines: (element in upper case; Kα1, Kα2, Kβ in Å). The first five are the columns of `ANOMALOUS`.
pub static ANODES: [(&str, f64, f64, f64); {len(waves)}] = [""")
for e in waves:
    print(f'    ("{e[0].strip()}", {f(e[1])}, {f(e[2])}, {f(e[3])}),')
print("];")
print(f"""
/// Atomic weights: (symbol in CrysFML's upper case; AtWe, g/mol).
pub static ATOMIC_WEIGHT: [(&str, f64); {len(chem)}] = [""")
for e in chem:
    print(f'    ("{e[0].strip()}", {f(e[3])}),')
print("];")
