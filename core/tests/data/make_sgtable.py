#!/usr/bin/env python3
"""Write core/src/sgtable.rs, the table of the 530 space-group settings, from spglib's table of settings
(`get_spacegroup_type`, spglib's spg_database.c, BSD-3; https://spglib.readthedocs.io).

Only the names are taken: Hall number, space-group number, setting choice, Hall symbol, spglib's name string and
full Hermann–Mauguin symbol. The short symbol of a setting is spglib's, except for the monoclinic groups (3–15),
where spglib gives the type's symbol ("P 2_1/c" for all nine settings of no. 14): there it is the full symbol with
its 1s dropped ("P 1 2_1/a 1" -> "P 2_1/a"), the International Tables' rule for monoclinic short symbols. The operators are not: the core generates them from the Hall symbol, and the test
`spglib_530` compares them with spglib's separately stored operators.
Run: python3 make_sgtable.py > ../../src/sgtable.rs   (needs `pip install spglib`)."""
import spglib

def q(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'

print(f"""//! The 530 space-group settings: Hall number, number, setting choice, Hall symbol, Hermann–Mauguin symbols.
//!
//! Source: spglib {spglib.__version__}, `get_spacegroup_type` (the table of settings in spglib's `spg_database.c`,
//! BSD-3 licence, https://spglib.readthedocs.io), written by `core/tests/data/make_sgtable.py`. Do not edit.
//! Only names come from there; the operators are generated from the Hall symbols by `spacegroup`.

/// One setting of a space-group type.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Setting {{
    /// Serial number 1–530 (spglib's Hall number).
    pub hall_number: u16,
    /// Space-group type number 1–230 (International Tables, Vol. A).
    pub number: u16,
    /// Setting choice: unique axis and cell choice, origin choice 1 or 2, H or R axes; "" for the standard one.
    pub choice: &'static str,
    /// Hall symbol, from which the operators are generated.
    pub hall: &'static str,
    /// Short Hermann–Mauguin symbol of this setting, spaces between the parts, screw axes as 2_1.
    pub hm: &'static str,
    /// spglib's name string, alternatives separated by " = " (e.g. "P 2_1/c = P 1 2_1/a 1").
    pub international: &'static str,
    /// Full Hermann–Mauguin symbol.
    pub hm_full: &'static str,
}}

pub static SETTINGS: [Setting; 530] = [""")
for h in range(1, 531):
    t = spglib.get_spacegroup_type(h)
    if 3 <= t.number <= 15:
        parts = t.international_full.split()
        hm = " ".join([parts[0]] + [x for x in parts[1:] if x != "1"])
    else:
        hm = t.international.split(" = ")[0]
    print(f"    Setting {{ hall_number: {h}, number: {t.number}, choice: {q(t.choice)}, hall: {q(t.hall_symbol)}, "
          f"hm: {q(hm)}, international: {q(t.international)}, hm_full: {q(t.international_full)} }},")
print("];")
