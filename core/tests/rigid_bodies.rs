//! Rigid bodies (P41, Publish's P47): `rigid_body` against FullProf's placements of the groups of the 17 rigid-body
//! `.pcr` files under `manuals/fullprof` (Jbt = 4, none with −4) and its prints of citric acid's molecule.
//!
//! - `euler_matrix_as_fullprof_prints_it`: R(θ, φ, χ) = R_z(φ) R_y(θ) R_z(χ) against the Euler matrices FullProf prints
//!   for citric acid (`citric-acid-sls_1.out`, lines 413–416) and for urea's Z-matrix (`urea.tpcr`, Mol_tpcr's
//!   output), to their five decimals.
//! - `z_matrix_as_fullprof_converts_it`: citric acid's Z-matrix of 21 atoms in the molecule's Cartesian frame against
//!   the table FullProf prints from it (`.out`, lines 424–445; the internal values the `.out` prints, so nothing
//!   rounded comes in): 19 within the print's 5·10⁻⁶ Å and 10⁻⁶ for FullProf's own arithmetic (four reach 5.0–5.5·10⁻⁶;
//!   single-precision reals in STRMOD, believed, not checked); A19 and A20, the two oxygens placed through hydrogens at bond
//!   angles of 41° and 124°, at 1.6·10⁻⁵ and 2.8·10⁻⁵ Å: not explained, counted and frozen (L345).
//! - `groups_placed_as_fullprof_wrote_them`: each file's atoms placed from its groups' values (the spherical,
//!   Cartesian, cylindrical and Z-matrix kinds, degrees and radians, satellites) against the fractional x, y, z FullProf
//!   wrote beside them (P1–P3, "for information", manual p. 35), atom by atom within the bound their prints allow: the
//!   rounding of every number the placement reads (half a unit of its last printed decimal, each through its
//!   derivative, summed) and of P1–P3's five decimals. 13 files within it, all their atoms; the 4 that depart are
//!   frozen with what was found:
//!   - `urea_test.pcr` (4 atoms, the satellites' hydrogens): its satellites' χ = ∓0.2 rad, but its P1–P3 are those of
//!     χ = 0, every atom within the bound then: the χ were written after FullProf placed the atoms;
//!   - `naph_san.pcr` (6 atoms, up to 4.5 × the bound, 1.7·10⁻³ Å): its positions are not a rigid image of its
//!     internal values (the best orientation and origin fitted to them leave 3.9·10⁻⁴ Å, scratch fit);
//!   - `ureaz.pcr` (3 atoms, 1.4 ×): FullProf's own Ur3 and Ur4 are not mirror images (0.14466 0.64466 against
//!     −0.14468 0.35536) where their internal values are (the same distance and angle, the torsion 180°);
//!   - `PPH3-sls.pcr` and its backup (10 atoms, 1.4 ×, up to 5.2·10⁻⁴ Å): not explained.
//!
//!   FullProf's runs of these files with NCY = 0, writing their P1–P3 again from their own values, would settle them.
//! - `satellites_hang_from_their_parent`: a satellite's origin and the polar angles of its z axis from its parent's atoms
//!   (P16 = N1.N2), against the X0, θ and φ FullProf wrote in the urea files.

use drs::cell::Cell;
use drs::rigid_body::{groups, positions, rotation, satellite_axis};
use std::path::{Path, PathBuf};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof")
}

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display())).iter().map(|&b| b as char).collect()
}

/// The 17 files and, for those that depart, (atoms over the bound, the largest ratio to it).
const FILES: [(&str, Option<(usize, f64)>); 17] = [
    ("8-FullProf-Tutorial-RigidBodies/Naphthalene/naph_san.pcr", Some((6, 4.5))),
    ("8-FullProf-Tutorial-RigidBodies/Urea/urea.pcr", None),
    ("8-FullProf-Tutorial-RigidBodies/Urea/urea_conv.pcr", None),
    ("8-FullProf-Tutorial-RigidBodies/Urea/urea_or.pcr", None),
    ("8-FullProf-Tutorial-RigidBodies/Urea/urea_start.pcr", None),
    ("8-FullProf-Tutorial-RigidBodies/Urea/urea_test.pcr", Some((4, 1900.0))),
    ("8-FullProf-Tutorial-RigidBodies/Urea/ureac.pcr", None),
    ("8-FullProf-Tutorial-RigidBodies/Urea/ureacyl.pcr", None),
    ("8-FullProf-Tutorial-RigidBodies/Urea/ureaz.pcr", Some((3, 1.4))),
    ("9-SAnn-OrganicCompounds/CitricAcid_example/Rietveld/citric-acid-sls_1.pcr", None),
    ("9-SAnn-OrganicCompounds/CitricAcid_example/SAN/P21c/ZMatrix/citric-acid-sls_1_SAN.pcr", None),
    ("9-SAnn-OrganicCompounds/PPH3/Rietveld/PPH3-sls.pcr", Some((10, 1.4))),
    ("9-SAnn-OrganicCompounds/PPH3/Rietveld/PPH3-sls_backup.pcr", Some((10, 1.4))),
    ("9-SAnn-OrganicCompounds/PPH3/SAN/PPH3-sls_1_SAN.pcr", None),
    ("FP_tests/hmt.pcr", None),
    ("FP_tests/pdd.pcr", None),
    ("FP_tests/urea.pcr", None),
];

fn rigid_phases(file: &str) -> (String, Vec<drs::pcr::PcrPhase>) {
    let text = latin1(&root().join(file));
    let p = drs::pcr::read(&text).unwrap_or_else(|e| panic!("{file}: {e:?}"));
    let phases: Vec<_> = p.phases.into_iter().filter(|ph| ph.flag("Jbt").abs() == 4).collect();
    (text, phases)
}

fn cell_of(ph: &drs::pcr::PcrPhase) -> Cell {
    let c = ph.cell;
    Cell::new(c[0], c[1], c[2], c[3], c[4], c[5]).unwrap()
}

/// Half a unit of the last decimal of the n-th number of a line (0 past its numbers).
fn half_unit(line: &str, n: usize) -> f64 {
    let toks: Vec<&str> = line.split('#').next().unwrap().split_whitespace().filter(|t| t.parse::<f64>().is_ok()).collect();
    toks.get(n).map_or(0.0, |t| 0.5 * 10f64.powi(-(t.split('.').nth(1).map_or(0, |d| d.len()) as i32)))
}

/// A change of one of a group's values, by h.
type Nudge = Box<dyn Fn(&mut drs::rigid_body::Group, f64)>;

fn norm(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

#[test]
#[allow(clippy::approx_constant)] // FullProf's printed matrices, 0.70711 among them
fn euler_matrix_as_fullprof_prints_it() {
    // (φ, θ, χ) as printed, the matrix's rows.
    let citric = rotation(218.194f64.to_radians(), 25.219f64.to_radians(), 24.060f64.to_radians());
    let printed = [[-0.82295, -0.09919, -0.55939], [0.06306, 0.96261, -0.26346], [0.56461, -0.25209, -0.78592]];
    let urea = rotation(90.0f64.to_radians(), 135.00002f64.to_radians(), (-180.0f64).to_radians());
    let printed_urea = [[0.0, 0.70711, -0.70711], [0.0, 0.70711, 0.70711], [1.0, 0.0, 0.0]];
    for (r, p) in [(citric, printed), (urea, printed_urea)] {
        for i in 0..3 {
            for j in 0..3 {
                assert!((r[i][j] - p[i][j]).abs() <= 5.1e-6, "R[{i}][{j}] = {} against {}", r[i][j], p[i][j]);
            }
        }
    }
}

#[test]
fn z_matrix_as_fullprof_converts_it() {
    let dir = root().join("9-SAnn-OrganicCompounds/CitricAcid_example/Rietveld");
    let (_, phases) = rigid_phases("9-SAnn-OrganicCompounds/CitricAcid_example/Rietveld/citric-acid-sls_1.pcr");
    let mut g = groups(&phases[0]).unwrap().remove(0);
    let out = latin1(&dir.join("citric-acid-sls_1.out"));
    let lines: Vec<&str> = out.lines().collect();
    let z = lines.iter().position(|l| l.contains("Type of Molecular description: Z-MATRIX")).unwrap();
    let c = lines.iter().position(|l| l.contains("Type of Molecular description: CARTESIAN COORDINATES")).unwrap();
    let row = |l: &str| -> [f64; 3] {
        let t: Vec<f64> = l.split_whitespace().skip(2).take(3).map(|x| x.parse().unwrap()).collect();
        [t[0], t[1], t[2]]
    };
    for (k, m) in g.members.iter_mut().enumerate() {
        m.internal = row(lines[z + 5 + k]);
    }
    let mut over = Vec::new();
    for (k, x) in g.molecular().iter().enumerate() {
        let p = row(lines[c + 5 + k]);
        let d = (0..3).map(|i| (x[i] - p[i]).abs()).fold(0.0, f64::max);
        if d > 6e-6 {
            over.push((g.members[k].label.clone(), d));
        }
    }
    println!("over the print: {over:?}");
    assert_eq!(over.iter().map(|o| o.0.as_str()).collect::<Vec<_>>(), ["A19", "A20"], "frozen (L345)");
    assert!(over.iter().all(|o| o.1 < 3.0e-5), "frozen at 1.6e-5 and 2.8e-5 Å");
}

#[test]
fn groups_placed_as_fullprof_wrote_them() {
    let mut report = Vec::new();
    for (file, frozen) in FILES {
        let (text, phases) = rigid_phases(file);
        assert_eq!(phases.len(), 1, "{file}: one rigid-body phase");
        let ph = &phases[0];
        let lines: Vec<&str> = text.lines().collect();
        let cell = cell_of(ph);
        let m = cell.crystal_to_orthonormal();
        let gs = groups(ph).unwrap_or_else(|e| panic!("{file}: {e}"));
        assert_eq!(gs.iter().map(|g| g.members.len()).sum::<usize>(), ph.atoms.len(), "{file}: every atom in a group");
        let placed = positions(ph).unwrap();
        // The bound: each number the placement reads, through its derivative, half a unit of its last decimal.
        let mut bound: Vec<[f64; 3]> = ph.atoms.iter().map(|_| std::array::from_fn(|i| (0..3).map(|j| m[i][j].abs() * 0.5e-5).sum())).collect();
        for (gi, g) in gs.iter().enumerate() {
            let head = &ph.atoms[g.members[0].index];
            let (l1, l3) = (lines[head.line - 1], lines[head.line + 1]);
            let mut inputs: Vec<(Nudge, f64)> = vec![
                (Box::new(|g, h| g.theta += h), half_unit(l1, 6)),
                (Box::new(|g, h| g.phi += h), half_unit(l1, 7)),
                (Box::new(|g, h| g.origin[0] += h), half_unit(l3, 3)),
                (Box::new(|g, h| g.origin[1] += h), half_unit(l3, 4)),
                (Box::new(|g, h| g.origin[2] += h), half_unit(l3, 5)),
                (Box::new(|g, h| g.chi += h), half_unit(l3, 6)),
            ];
            for (k, mm) in g.members.iter().enumerate() {
                let l = lines[ph.atoms[mm.index].line + 1];
                for c in 0..3 {
                    inputs.push((Box::new(move |g, h| g.members[k].internal[c] += h), half_unit(l, c)));
                }
            }
            let base = g.cartesian(&cell);
            for (f, du) in inputs {
                let h = 1e-7;
                let mut g2 = gs[gi].clone();
                f(&mut g2, h);
                for (k, x) in g2.cartesian(&cell).iter().enumerate() {
                    let b = &mut bound[g.members[k].index];
                    for i in 0..3 {
                        b[i] += ((x[i] - base[k][i]) / h).abs() * du;
                    }
                }
            }
        }
        let mut over = 0;
        let (mut worst, mut largest) = (0.0f64, 0.0f64);
        for (k, a) in ph.atoms.iter().enumerate() {
            let d = drs::mat3::mul_vec(&m, &[placed[k][0] - a.x[0], placed[k][1] - a.x[1], placed[k][2] - a.x[2]]);
            let r = (0..3).map(|i| d[i].abs() / bound[k][i]).fold(0.0, f64::max);
            worst = worst.max(r);
            largest = largest.max(norm(d));
            over += usize::from(r > 1.0);
        }
        report.push(format!("{file}: {} atoms, {over} over, worst {worst:.3} × the bound, {largest:.2e} Å", ph.atoms.len()));
        match frozen {
            None => assert_eq!(over, 0, "{file}: {worst} × the bound"),
            Some((n, r)) => assert!(over <= n && worst <= r, "{file}: frozen at {n} atoms, {r} × the bound; now {over}, {worst}"),
        }
    }
    println!("{}", report.join("\n"));
    // urea_test: its P1–P3 are those of satellites with χ = 0.
    let (text, phases) = rigid_phases("8-FullProf-Tutorial-RigidBodies/Urea/urea_test.pcr");
    let ph = &phases[0];
    let mut gs = groups(ph).unwrap();
    let _ = text;
    let cell = cell_of(ph);
    let m = cell.crystal_to_orthonormal();
    for g in gs.iter_mut().filter(|g| g.satellite.is_some()) {
        assert!((g.chi.abs() - 0.2).abs() < 1e-12);
        g.chi = 0.0;
        for (mm, x) in g.members.iter().zip(g.fractional(&cell)) {
            let w = ph.atoms[mm.index].x;
            assert!(norm(drs::mat3::mul_vec(&m, &[x[0] - w[0], x[1] - w[1], x[2] - w[2]])) < 6e-5, "{}", mm.label);
        }
    }
}

#[test]
fn satellites_hang_from_their_parent() {
    let mut n = 0;
    for (file, _) in FILES.iter().filter(|f| f.0.contains("Urea/") || f.0.ends_with("FP_tests/urea.pcr")) {
        let (_, phases) = rigid_phases(file);
        let ph = &phases[0];
        let cell = cell_of(ph);
        let gs = groups(ph).unwrap();
        for g in (0..gs.len()).filter(|&g| gs[g].satellite.is_some()) {
            let (o, t, p) = satellite_axis(&gs, g, &cell).unwrap();
            let s = &gs[g];
            let tol = if s.degrees { 5e-3 } else { 5e-5 };
            assert!((0..3).all(|i| (o[i] - s.origin[i]).abs() <= 1e-5), "{file}: {} origin {o:?} against {:?}", s.members[0].label, s.origin);
            assert!((t - s.theta).abs() <= tol && (p - s.phi).abs() <= tol, "{file}: {} θ, φ = {t}, {p} against {}, {}", s.members[0].label, s.theta, s.phi);
            n += 1;
        }
    }
    assert_eq!(n, 16, "two satellites in each of the eight urea files with them (ureaz has none)");
    // `drs rigid`: the groups, the satellites' anchors and each atom against the file.
    let o = std::process::Command::new(env!("CARGO_BIN_EXE_drs")).arg("rigid").arg(root().join("FP_tests/urea.pcr")).output().unwrap();
    let out = String::from_utf8(o.stdout).unwrap();
    assert!(o.status.success() && out.contains("3 groups") && out.contains("satellite of group 1 at its atom 3 (from 1)"), "{out}");
    assert_eq!(out.lines().filter(|l| l.ends_with(" Å")).count(), 8);
}
