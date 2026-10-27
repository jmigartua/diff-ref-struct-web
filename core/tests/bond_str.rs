//! Distances, angles and bond-valence sums (P12; L081) against FullProf's own, written by Bond_Str in `CODFILn.dis`
//! (manual p. 11) from the CFL file FullProf writes after a refinement, which each `.dis` repeats at its top.
//!
//! The `.dis` files: `7-FullProf-Tutorial-SymmetryModes/LaMnO3` (two, Pnma), `runs/Yb2ZrWO8` and
//! `runs/Yb2ZrWO8-PO-March` (three, Fm-3m with three cations on one site), and the sequential refinement
//! `1-WinPLOTR/2_FP_sequential_refinements/LaMnO3_D20` (71, Pbnm, with angles); for the distances only, case A's
//! two (`runs/modes/A_modes`, P2₁/n, L114), whose CFL gives every atom FullProf's species "SR+2".
//!
//! Each number is compared within the precision Bond_Str prints it to: a value with k decimals is within h =
//! ½·10^(−k) of FullProf's, plus FullProf's single precision (8 roundings of 2^(−24) of the value, a bound on the
//! few operations of a distance or a sum); a standard deviation printed as e units of the last digit is within
//! (e ± ½)·10^(−k). Two of Bond_Str's numbers lose more to single precision, bounded in `angle_agrees` (an angle of a
//! flat triangle, 0.03° where it is 0) and `distortion_agrees` (a difference of two numbers near 1, times 10⁴).
//! - `distances`: every atom's neighbours, the same images at the same positions, each distance, and each standard
//!   deviation by Bond_Str's formula (`bond_str_distance_sigma`); the default (the first-order propagation, which
//!   moves both ends of an atom bonded to its own image together) differs from it at six distances, all of that kind;
//! - `angles`: every angle of every triangle, its value, and its standard deviation by Bond_Str's formula from the
//!   three distances (`bond_str_angle_sigma`), the distances' without the cell's standard deviations, which Bond_Str
//!   leaves out (`bond_str_sigma`: an Mn–Mn distance of c/2 is printed with none); the default standard deviation,
//!   the angle's own derivative with the cell's, differs at 2047 of 2160;
//! - `valence_sums`: with Bond_Str's occupancies (normalised to the first atom), each bond's valence, the
//!   coordination, the mean distance, the distortion, the sum with its standard deviation, the instability indices.

use drs::bonds::{self, Contact};
use drs::cfl;
use std::path::{Path, PathBuf};

fn root(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)
}

fn files() -> Vec<PathBuf> {
    let mut out = Vec::new();
    for d in ["../manuals/fullprof/7-FullProf-Tutorial-SymmetryModes/LaMnO3", "../manuals/fullprof/runs/Yb2ZrWO8", "../manuals/fullprof/runs/Yb2ZrWO8-PO-March", "../manuals/fullprof/1-WinPLOTR/2_FP_sequential_refinements/LaMnO3_D20"] {
        let mut v: Vec<PathBuf> = std::fs::read_dir(root(d)).unwrap().map(|e| e.unwrap().path()).filter(|p| p.extension().is_some_and(|e| e == "dis")).collect();
        v.sort();
        out.extend(v);
    }
    out
}

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap().iter().map(|&b| b as char).collect()
}

/// "1.9715(7)", "2.6986(  11)" → (value, esd, decimals, printed esd units).
fn printed(s: &str) -> (f64, f64, i32, f64) {
    let s: String = s.split_whitespace().collect();
    let v = s.split('(').next().unwrap();
    let k = v.split('.').nth(1).map_or(0, str::len) as i32;
    let e: f64 = s.split('(').nth(1).map_or(0.0, |t| t.trim_end_matches(')').parse().unwrap());
    (v.parse().unwrap_or_else(|_| panic!("not a number: {s}")), e * 10f64.powi(-k), k, e)
}

/// Within the printed precision of FullProf's value (single precision included).
fn agrees(ours: f64, p: (f64, f64, i32, f64)) -> bool {
    (ours - p.0).abs() <= 0.5 * 10f64.powi(-p.2) + 8.0 * ours.abs() * 2f64.powi(-24)
}

/// An angle within the printed precision, and within what single precision does to an angle FullProf takes as the arc
/// cosine of a cosine from three distances (some 8 roundings, δ = 8·2^(−24)): Δθ ≤ δ/|sin θ|, and √(2δ) where the
/// triangle is flat (θ near 0° or 180°, where the arc cosine's slope is infinite).
fn angle_agrees(ours: f64, p: (f64, f64, i32, f64)) -> bool {
    let delta = 8.0 * 2f64.powi(-24);
    let single = (delta / ours.to_radians().sin().abs().max(1e-12)).min((2.0 * delta).sqrt()).to_degrees();
    (ours - p.0).abs() <= 0.5 * 10f64.powi(-p.2) + single
}

/// Bond_Str's distortion 10⁴(⟨d²⟩/⟨d⟩² − 1) from n bonds, in single precision: the ratio, near 1, carries the
/// roundings of the n squares and n distances summed, the mean squared and the division (2n + 4 of 2^(−24)), times 10⁴.
fn distortion_agrees(ours: f64, n: usize, p: (f64, f64, i32, f64)) -> bool {
    (ours - p.0).abs() <= 0.5 * 10f64.powi(-p.2) + 1e4 * (2 * n + 4) as f64 * 2f64.powi(-24)
}

/// A standard deviation within half a unit of the printed digit.
fn esd_agrees(ours: f64, p: (f64, f64, i32, f64)) -> bool {
    (ours - p.1).abs() <= 0.5 * 10f64.powi(-p.2) + 1e-12
}

/// A printed distance: the neighbour, the distance as printed, the position.
type Row = (usize, (f64, f64, i32, f64), [f64; 3]);

struct Dis {
    cfl: cfl::Cfl,
    text: String,
}

fn read(path: &Path) -> Dis {
    let text = latin1(path);
    let start = text.find("Content of the FullProf generated CFL file").expect("the CFL in the .dis");
    let end = text.find("=> PARVAL").or_else(|| text.find("Metric information")).unwrap();
    let body: String = text[start..end].lines().skip(1).map(|l| l.trim()).collect::<Vec<_>>().join("\n");
    Dis { cfl: cfl::read(&body).unwrap_or_else(|e| panic!("{}: {e}", path.display())), text }
}

/// The atom index of a label.
fn index(d: &Dis, label: &str) -> usize {
    d.cfl.structure.atoms.iter().position(|a| a.label == label).unwrap_or_else(|| panic!("label {label}"))
}

fn near(a: [f64; 3], b: [f64; 3], tol: f64) -> bool {
    (0..3).all(|k| (a[k] - b[k]).abs() <= tol)
}

/// Counted and frozen (L114): case A's two .dis (the symmetry-mode refinement of Sr2ScSbO6), whose standard
/// deviations stand above Bond_Str's formula from the CFL's printed ones at 40 of its 132 distances (Sb1–O11_2 0.0016
/// where the formula gives 0.00152, O11–O11_3 0.003 for 0.0024), by at most 0.32 of a last digit past the print's half
/// unit; the distances agree, all 132. Neither the CFL's printed esd's nor the unrounded ones of the amplitudes taken
/// as independent (smaller: O11's x 0.000253 against "(3)") nor one factor on either reproduces them: unexplained,
/// listed for Josu. The default σ is off at the same 40.
const N_CASE_A_DEFAULT_OFF: usize = 40;
const N_FORMULA_OFF: usize = 40;
const FORMULA_EXCESS: f64 = 0.33;

#[test]
fn distances() {
    let (mut n_files, mut n_dist, mut n_sigma_off) = (0, 0, 0);
    // Case A's two .dis (L114): the mode refinement's CFL, whose species are all "SR+2" (FullProf's), so only its
    // distances are compared, not its valence sums.
    let case_a = root("../manuals/fullprof/runs/modes/A_modes");
    let mut paths = files();
    paths.extend(["a_modes_1.dis", "a_modes_calc_1.dis"].map(|f| case_a.join(f)));
    // Bond_Str's σ beyond its formula from the CFL's printed standard deviations: (count, largest excess in units of
    // the last digit beyond the print's half unit).
    let (mut n_formula_off, mut formula_excess, mut n_case_a): (usize, f64, usize) = (0, 0.0, 0);
    for path in paths {
        let d = read(&path);
        let st = &d.cfl.structure;
        let ours = bonds::contacts(st, Some(&d.cfl.uncertainties), d.cfl.dmax.0);
        // FullProf's rows: "  i  k  j  (A)-(B): d(s)  x y z  (t) op", in sections by atom.
        let mut theirs: Vec<Vec<Row>> = vec![Vec::new(); st.atoms.len()];
        let mut current = None;
        for l in d.text.lines() {
            if let Some(r) = l.split("to atom:").nth(1) {
                current = Some(index(&d, r.split_whitespace().next().unwrap()));
                continue;
            }
            if l.contains("BOND-VALENCE") || l.contains("Angles around atom") {
                current = None;
            }
            let Some(i) = current else { continue };
            let Some((head, tail)) = l.split_once("):") else { continue };
            if head.split_whitespace().take(3).filter(|t| t.parse::<usize>().is_ok()).count() != 3 {
                continue;
            }
            let k = index(&d, head.rsplit("-(").next().unwrap().trim());
            // The distance with its esd, then x y z.
            let tail = tail.replace("( ", "(").replace("( ", "(");
            let f: Vec<&str> = tail.split_whitespace().collect();
            let pos = [f[1], f[2], f[3]].map(|t| t.parse::<f64>().unwrap());
            theirs[i].push((k, printed(f[0]), pos));
        }
        for (i, list) in theirs.iter().enumerate() {
            let mine: &Vec<Contact> = &ours[i];
            assert_eq!(mine.len(), list.len(), "{}: atom {}: {} neighbours, FullProf {}", path.display(), st.atoms[i].label, mine.len(), list.len());
            for (k, p, pos) in list {
                let c = mine.iter().find(|c| c.atom == *k && near(c.position, *pos, 6e-5)).unwrap_or_else(|| panic!("{}: {}–{} at {pos:?} not found", path.display(), st.atoms[i].label, st.atoms[*k].label));
                assert!(agrees(c.distance, *p), "{}: {}–{}: {} against {:?}", path.display(), st.atoms[i].label, st.atoms[*k].label, c.distance, p);
                if p.3 > 0.0 && !esd_agrees(c.sigma, *p) {
                    // Explained: an atom bonded to its own image, whose two ends Bond_Str takes as independent.
                    let fp = bonds::bond_str_distance_sigma(st, &d.cfl.uncertainties, i, c);
                    println!("{}: {}–{} {:?}: σ {:.5}, Bond_Str's formula {:.5}, printed {:?}", path.file_name().unwrap().to_string_lossy(), st.atoms[i].label, st.atoms[*k].label, pos, c.sigma, fp, p);
                    assert!((c.atom == i && esd_agrees(fp, *p)) || path.starts_with(&case_a), "{}: an unexplained standard deviation", path.display());
                    n_sigma_off += 1;
                }
                if p.3 > 0.0 {
                    let fp = bonds::bond_str_distance_sigma(st, &d.cfl.uncertainties, i, c);
                    if !esd_agrees(fp, *p) {
                        // Counted and frozen (L114): case A only, `N_FORMULA_OFF`.
                        println!("{}: {}–{} {:?}: Bond_Str's formula {:.6}, printed {:?}", path.file_name().unwrap().to_string_lossy(), st.atoms[i].label, st.atoms[*k].label, pos, fp, p);
                        assert!(path.starts_with(&case_a), "{}: Bond_Str's σ", path.display());
                        n_formula_off += 1;
                        formula_excess = formula_excess.max(((fp - p.1).abs() - 0.5 * 10f64.powi(-p.2)) * 10f64.powi(p.2));
                    }
                }
                n_dist += 1;
                n_case_a += usize::from(path.starts_with(&case_a));
            }
        }
        n_files += 1;
    }
    println!("{n_files} .dis files, {n_dist} distances as Bond_Str's; {n_sigma_off} standard deviations off");
    assert!(n_files >= 75);
    // Counted and frozen (L081): the six O1–O1 distances of LaMnO3_cub, an atom to its own image; (L114) case A's 40.
    println!("Bond_Str's formula off at {n_formula_off} of case A's {n_case_a}, by {formula_excess:.3} of a last digit past the print's half unit");
    assert_eq!(n_sigma_off, 6 + N_CASE_A_DEFAULT_OFF);
    assert!(n_formula_off <= N_FORMULA_OFF && formula_excess <= FORMULA_EXCESS, "{n_formula_off} off by {formula_excess}");
}

/// "( O1   ) :  0.57740  0.01360 -0.25000  ( O2   ) :  0.22540  0.19130 -0.03930" → the two positions.
fn two_positions(l: &str) -> Option<([f64; 3], [f64; 3])> {
    let parts: Vec<&str> = l.split(") :").collect();
    if parts.len() != 3 {
        return None;
    }
    let a: Vec<f64> = parts[1].split_whitespace().take(3).filter_map(|t| t.parse().ok()).collect();
    let b: Vec<f64> = parts[2].split_whitespace().take(3).filter_map(|t| t.parse().ok()).collect();
    (a.len() == 3 && b.len() == 3).then(|| ([a[0], a[1], a[2]], [b[0], b[1], b[2]]))
}

#[test]
fn angles() {
    let (mut n_angles, mut n_default_off) = (0, 0);
    for path in files() {
        let d = read(&path);
        if d.cfl.dmax.1 <= 0.0 {
            continue;
        }
        let st = &d.cfl.structure;
        let u = &d.cfl.uncertainties;
        let contacts = bonds::contacts(st, Some(u), d.cfl.dmax.0);
        let ours = bonds::angles(st, Some(u), &contacts, d.cfl.dmax.1);
        let lines: Vec<&str> = d.text.lines().collect();
        let mut centre = None;
        let mut n_file = 0;
        for (n, l) in lines.iter().enumerate() {
            if let Some(r) = l.split("Angles around atom:").nth(1) {
                centre = Some(index(&d, r.split_whitespace().next().unwrap()));
                continue;
            }
            let (Some(i), Some(dl)) = (centre, l.split_once("d12 =")) else { continue };
            let ds: Vec<String> = dl.1.replace("d23 =", "|").replace("d13 =", "|").split('|').map(|t| t.trim().to_string()).collect();
            let (d1, d2, d12) = (printed(&ds[0]), printed(&ds[1]), printed(&ds[2]));
            let ang = |m: usize| printed(lines[n + m].rsplit("):").next().unwrap());
            let (a_centre, a_second, a_first) = (ang(1), ang(2), ang(3));
            let (p, q) = two_positions(lines[n + 4]).unwrap();
            let c = &contacts[i];
            // The pair in Bond_Str's order or the other (its operators come in another order than ours).
            let found = ours.iter().find_map(|a| {
                if a.centre != i {
                    return None;
                }
                let (x, y) = (c[a.first].position, c[a.second].position);
                if near(x, p, 6e-5) && near(y, q, 6e-5) {
                    Some((a, false))
                } else if near(x, q, 6e-5) && near(y, p, 6e-5) {
                    Some((a, true))
                } else {
                    None
                }
            });
            let (a, swapped) = found.unwrap_or_else(|| panic!("{}: angle at {} {p:?} {q:?} not found", path.display(), st.atoms[i].label));
            let (cj, ck, at_j, at_k) = if swapped { (&c[a.second], &c[a.first], a.at_second.0, a.at_first.0) } else { (&c[a.first], &c[a.second], a.at_first.0, a.at_second.0) };
            for (ours, theirs, what) in [(cj.distance, d1, "d12"), (ck.distance, d2, "d23"), (a.between.0, d12, "d13")] {
                assert!(agrees(ours, theirs), "{}: {what} at {}: {ours} against {theirs:?}", path.display(), st.atoms[i].label);
            }
            for (ours, theirs, what) in [(a.angle.0, a_centre, "angle"), (at_j, a_first, "at the first"), (at_k, a_second, "at the second")] {
                assert!(angle_agrees(ours, theirs), "{}: {what} at {}: {ours} against {theirs:?}", path.display(), st.atoms[i].label);
            }
            // Bond_Str's standard deviations, from its distances' (bond_str_angle_sigma).
            let s = |a: (usize, [f64; 3]), b: (usize, [f64; 3])| bonds::bond_str_sigma(st, u, a, b);
            let o = (i, st.atoms[i].x);
            let (s1, s2, s12) = (s(o, (cj.atom, cj.position)), s(o, (ck.atom, ck.position)), s((cj.atom, cj.position), (ck.atom, ck.position)));
            let (t1, t2, t12) = ((cj.distance, s1), (ck.distance, s2), (a.between.0, s12));
            for (sig, theirs, what) in [(s1, d1, "σ(d12)"), (s2, d2, "σ(d23)"), (s12, d12, "σ(d13)"), (bonds::bond_str_angle_sigma(t1, t2, t12), a_centre, "σ(angle)"), (bonds::bond_str_angle_sigma(t12, t1, t2), a_first, "σ(at the first)"), (bonds::bond_str_angle_sigma(t12, t2, t1), a_second, "σ(at the second)")] {
                assert!(esd_agrees(sig, theirs), "{}: {what} at {}: {sig} against {theirs:?}", path.display(), st.atoms[i].label);
            }
            if !esd_agrees(a.angle.1, a_centre) {
                n_default_off += 1;
            }
            n_angles += 1;
            n_file += 1;
        }
        assert_eq!(n_file, ours.len(), "{}: every angle", path.display());
    }
    println!("{n_angles} triangles as Bond_Str's, values and standard deviations; the default σ of the angle (its own derivative) differs from Bond_Str's at {n_default_off}");
    assert!(n_angles > 1000);
    // Counted and frozen (L081): Bond_Str leaves out the cell's standard deviations and takes an angle's from the
    // triangle's three distances as if independent; the default does neither.
    assert!(n_default_off <= 2047, "{n_default_off}");
}

/// The numbers after a marker on a line: "Sums:    2.810(  4)" → "2.810(4)".
fn after<'a>(l: &'a str, marker: &str) -> Option<&'a str> {
    l.split(marker).nth(1).map(str::trim)
}

#[test]
fn valence_sums() {
    let (mut n_atoms, mut n_bonds) = (0, 0);
    for path in files() {
        let d = read(&path);
        let st = &d.cfl.structure;
        let u = &d.cfl.uncertainties;
        // Bond_Str's figures use its distances' standard deviations.
        let mut contacts = bonds::contacts(st, Some(u), d.cfl.dmax.0);
        for (i, list) in contacts.iter_mut().enumerate() {
            for c in list.iter_mut() {
                c.sigma = bonds::bond_str_distance_sigma(st, u, i, c);
            }
        }
        let occ = bonds::bond_str_occupancies(st).unwrap();
        let tol = d.cfl.tol.filter(|t| *t > 0.1).unwrap_or(20.0) / 100.0;
        let (sums, gii) = bonds::valence_sums(st, &contacts, &occ, tol).unwrap();
        let name = path.file_name().unwrap().to_string_lossy().to_string();
        // The parameters and radii Bond_Str printed are the table's.
        for l in d.text.lines().filter(|l| l.contains(" d0=")) {
            let f: Vec<&str> = l.split_whitespace().collect();
            let (cat, an) = (f[2], f[6]);
            let (d0, b0, _) = bonds::parameters(cat, an).unwrap_or_else(|| panic!("{name}: {cat}–{an}"));
            assert_eq!((format!("{d0:.3}"), format!("{b0:.3}")), (after(l, "d0=").unwrap().split_whitespace().next().unwrap().to_string(), after(l, "B0=").unwrap().split_whitespace().next().unwrap().to_string()), "{name}: {l}");
        }
        for l in d.text.lines().filter(|l| l.contains("(Eff. radius):")) {
            for part in l.split("radius):").skip(1) {
                let (sp, r) = part.trim().split_once('(').unwrap();
                let (el, q) = bonds::ion(sp.trim()).unwrap();
                assert_eq!(format!("{:.3}", bonds::ionic_radius(&el, q)), r.split(')').next().unwrap().trim(), "{name}: radius of {sp}");
            }
        }
        let lines: Vec<&str> = d.text.lines().collect();
        for (n, l) in lines.iter().enumerate() {
            let Some(r) = after(l, "coordination of atom:") else { continue };
            let label = r.split_whitespace().next().unwrap();
            let i = index(&d, label);
            let p_occ = printed(after(l, "occupancy:").unwrap());
            assert!(agrees(occ[i], p_occ), "{name}: occupancy of {label}: {} against {p_occ:?}", occ[i]);
            let s = sums.iter().find(|s| s.atom == i).unwrap_or_else(|| panic!("{name}: no sum for {label}"));
            // The bonds, in order.
            let rows: Vec<&str> = lines[n + 1..].iter().take_while(|l| !l.contains("Coordination number")).filter(|l| l.contains(")-(") && l.contains(") :")).copied().collect();
            assert_eq!(rows.len(), s.bonds.len(), "{name}: bonds of {label}");
            // Each printed bond is one of ours with the same neighbour and distance (the order is the operators').
            let mut used = vec![false; s.bonds.len()];
            for row in &rows {
                let neighbour = row.split(")-(").nth(1).unwrap().split(')').next().unwrap().trim();
                let t = row.rsplit(") :").next().unwrap().replace("( ", "(").replace("( ", "(").replace("( ", "(");
                let f: Vec<&str> = t.split_whitespace().collect();
                let (pd, pv) = (printed(f[0]), printed(f[1]));
                let m = (0..s.bonds.len()).find(|&m| !used[m] && st.atoms[contacts[i][s.bonds[m].contact].atom].label == neighbour && agrees(contacts[i][s.bonds[m].contact].distance, pd)).unwrap_or_else(|| panic!("{name}: {row}"));
                used[m] = true;
                let b = &s.bonds[m];
                let c = &contacts[i][b.contact];
                assert!(agrees(c.distance, pd) && esd_agrees(c.sigma, pd), "{name}: {label} bond {} against {pd:?}", c.distance);
                assert!(agrees(b.valence, pv) && esd_agrees(b.sigma, pv), "{name}: {label} bond valence {} ({}) against {pv:?}", b.valence, b.sigma);
                n_bonds += 1;
            }
            let block: Vec<&str> = lines[n..].iter().take(rows.len() + 20).copied().collect();
            let find = |m: &str| block.iter().find_map(|l| after(l, m)).unwrap_or_else(|| panic!("{name}: {m}"));
            let cn: usize = find("Coordination number:").split_whitespace().next().unwrap().parse().unwrap();
            let eff = printed(find("Eff.Coor. number:").split_whitespace().next().unwrap());
            let mean = printed(&find("Average distance  :").split("Distortion").next().unwrap().replace(' ', ""));
            let dist = printed(find("Distortion:").split_whitespace().next().unwrap());
            let sum = printed(&find("Sums:").replace(' ', ""));
            assert_eq!(cn, s.coordination, "{name}: {label}");
            assert!(distortion_agrees(s.distortion, s.coordination, dist), "{name}: {label} distortion: {} against {dist:?}", s.distortion);
            for (ours, theirs, what) in [(s.effective_coordination, eff, "effective coordination"), (s.mean_distance.0, mean, "mean distance"), (s.sum.0, sum, "valence sum")] {
                assert!(agrees(ours, theirs), "{name}: {label} {what}: {ours} against {theirs:?}");
            }
            assert!(esd_agrees(s.mean_distance.1, mean) && esd_agrees(s.sum.1, sum), "{name}: {label} σ: {:?} {:?}", s.mean_distance, s.sum);
            n_atoms += 1;
        }
        let gl = |m: &str| d.text.lines().find(|l| l.contains(m) && (l.trim_end().ends_with("/100") || l.trim_end().ends_with('%'))).map(|l| printed(l.rsplit('=').next().unwrap().split_whitespace().next().unwrap()));
        for (ours, theirs, what) in [(gii.old, gl("Old Global Instability"), "GII"), (gii.a, gl("Normalized   GII(a)"), "GII(a)"), (gii.b, gl("Normalized   GII(b)"), "GII(b)"), (gii.c, gl("Normalized   GII(c)"), "GII(c)")] {
            let theirs = theirs.unwrap_or_else(|| panic!("{name}: {what}"));
            assert!(agrees(ours, theirs), "{name}: {what}: {ours} against {theirs:?}");
        }
    }
    println!("{n_atoms} atoms' bond-valence sums and {n_bonds} bonds as Bond_Str's, with the instability indices");
    assert!(n_atoms > 250);
}

#[test]
fn from_the_runs() {
    // The structure FullProf refined, imported from its .pcr, gives the distances Bond_Str listed from the CFL FullProf
    // wrote after that refinement; the sums with the chemical occupancies are Bond_Str's where the first atom's site
    // is full, and twice them for Yb2ZrWO8 (Bond_Str's Yb is 1, the chemistry's 0.5; L081).
    for (dir, stem, dis) in [
        ("../manuals/fullprof/runs/Yb2ZrWO8", "ybfluo_Riet_IRF", "ybfluo_Riet_IRF_1.dis"),
        ("../manuals/fullprof/7-FullProf-Tutorial-SymmetryModes/LaMnO3", "LaMnO3ref", "LaMnO3ref_1.dis"),
        ("../manuals/fullprof/7-FullProf-Tutorial-SymmetryModes/LaMnO3", "LaMnO3_cub", "LaMnO3_cub_1.dis"),
    ] {
        let d = read(&root(dir).join(dis));
        let p = drs::pcr::read(&latin1(&root(dir).join(format!("{stem}.pcr")))).unwrap();
        let m = drs::powder::Model::from_pcr(&p, &root(dir), stem).unwrap().model;
        let drs::powder::Intensities::Structure(st, _) = &m.phases[0].intensities else { panic!() };
        let mut st = (**st).clone();
        st.cell = m.phases[0].cell;
        // The charges as the CFL gives them (a neutron .pcr names the elements only).
        for (a, b) in st.atoms.iter_mut().zip(&d.cfl.structure.atoms) {
            a.species = b.species.clone();
        }
        let (ours, theirs) = (bonds::contacts(&st, None, d.cfl.dmax.0), bonds::contacts(&d.cfl.structure, None, d.cfl.dmax.0));
        for (a, b) in ours.iter().zip(&theirs) {
            let mut x: Vec<f64> = a.iter().map(|c| c.distance).collect();
            let mut y: Vec<f64> = b.iter().map(|c| c.distance).collect();
            x.sort_by(f64::total_cmp);
            y.sort_by(f64::total_cmp);
            assert_eq!(x.len(), y.len(), "{stem}");
            // The CFL rounds each coordinate to the last digit of its standard deviation (at most σ/2 away), the
            // others to 5 decimals: two ends, each moving the distance by at most that times the longest cell edge.
            let c = &d.cfl.structure.cell;
            let worst = d.cfl.uncertainties.x.iter().flatten().fold(1e-5f64, |a, s| a.max(s / 2.0));
            let bound = 2.0 * worst * c.a.max(c.b).max(c.c);
            assert!(x.iter().zip(&y).all(|(u, v)| (u - v).abs() <= bound), "{stem}: {x:?} {y:?} (bound {bound})");
        }
        let chem = bonds::chemical_occupancies(&st).unwrap();
        let bs = bonds::bond_str_occupancies(&d.cfl.structure).unwrap();
        let (s1, _) = bonds::valence_sums(&st, &ours, &chem, 0.2).unwrap();
        let (s2, _) = bonds::valence_sums(&d.cfl.structure, &theirs, &bs, 0.2).unwrap();
        let ratio = s2[0].sum.0 / s1[0].sum.0;
        println!("{stem}: first atom's sum {:.3} (chemical occupancies) against Bond_Str's {:.3}", s1[0].sum.0, s2[0].sum.0);
        let want = if stem.starts_with("ybfluo") { 2.0 } else { 1.0 };
        assert!((ratio - want).abs() < 1e-3, "{stem}: {ratio}");
    }
}
