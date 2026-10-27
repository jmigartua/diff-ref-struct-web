//! The readers of FullProf's `.prf` and `.hkl` (`formats::fullprof_outputs`, app request 17, L100) over every such file
//! in `manuals/` (FullProf's tutorials and runs, and the mode course's case A):
//! - every `.prf` reads, but the integrated-intensity ones (simulated annealing), which are refused; the points are
//!   those its line 2 announces; each phase has the reflections line 3 counts for it (both wavelengths); its
//!   excluded regions are those of the `.out` written with it, where there is one;
//! - the background, FullProf's offset added back, is compared with the core's in `fullprof_profiles`, which reads the
//!   `.prf` files with this reader;
//! - every `.hkl` reads into blocks with rows; where a block heading counts its reflections ("444 reflections"), the
//!   block has at most that many rows (FullProf: "the # of eff. reflections may be lower") and at least one; where a run's `.hkl` and `.prf` come from the same pass, the `2T` of each reflection
//!   of the `.hkl` is a tick of the `.prf` (λ1, shifted by the `.prf`'s Zero, SyCos·cos 2θ and SySin·sin 2θ: FullProf's
//!   SyCos multiplies cos 2θ, not its manual's cos θ, L267).

use drs::formats::fullprof_outputs::{read_hkl, read_prf};
use std::path::{Path, PathBuf};

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap().iter().map(|&b| b as char).collect()
}

fn files(dir: &Path, ext: &str, out: &mut Vec<PathBuf>) {
    for e in std::fs::read_dir(dir).into_iter().flatten().flatten() {
        let p = e.path();
        if p.is_dir() {
            files(&p, ext, out);
        } else if p.extension().is_some_and(|x| x.eq_ignore_ascii_case(ext)) {
            out.push(p);
        }
    }
}

fn all(ext: &str) -> Vec<PathBuf> {
    let base = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals");
    let mut v = Vec::new();
    files(&base.join("fullprof"), ext, &mut v);
    files(&base.join("bcs/mode-crystallography-manual/course-um6p-2024/04_tools_materials"), ext, &mut v);
    v.sort();
    v
}

#[test]
fn prf_files() {
    let (mut read, mut refused, mut regions_checked) = (0, 0, 0);
    let mut fails = Vec::new();
    for f in all("prf") {
        let text = latin1(&f);
        let name = f.display().to_string();
        match read_prf(&text) {
            Err(e) => {
                if e.0.contains("integrated intensities") {
                    refused += 1;
                } else {
                    fails.push(format!("{name}: {}", e.0));
                }
            }
            Ok(p) => {
                read += 1;
                // Line 3: the reflections of each phase.
                let line3: Vec<usize> = text.lines().nth(2).unwrap().split_whitespace().filter_map(|s| s.parse().ok()).collect();
                for ph in 1..=p.phases {
                    let n = p.reflections.iter().filter(|r| r.phase == ph).count();
                    if line3.get(ph - 1) != Some(&n) {
                        fails.push(format!("{name}: phase {ph} has {n} reflections, line 3 says {:?}", line3.get(ph - 1)));
                    }
                }
                if p.yobs.len() != p.x.len() || p.background.len() != p.x.len() || p.x.windows(2).any(|w| w[1] <= w[0]) {
                    fails.push(format!("{name}: columns of unequal length or 2θ not increasing"));
                }
                // The excluded regions against the `.out` written with the .prf.
                let out = f.with_extension("out");
                if out.is_file() {
                    let o = latin1(&out);
                    if let Some(i) = o.lines().position(|l| l.contains("Excluded regions for Pattern")) {
                        let regions: Vec<(f64, f64)> = o
                            .lines()
                            .skip(i + 2)
                            .map_while(|l| {
                                let v: Vec<f64> = l.split_whitespace().filter_map(|s| s.parse().ok()).collect();
                                (v.len() == 2).then(|| (v[0], v[1]))
                            })
                            .collect();
                        regions_checked += 1;
                        let same = regions.len() == p.excluded.len() && regions.iter().zip(&p.excluded).all(|(a, b)| (a.0 - b.0).abs() < 1e-3 && (a.1 - b.1).abs() < 1e-3);
                        if !same {
                            fails.push(format!("{name}: excluded regions {:?}, the .out's {regions:?}", p.excluded));
                        }
                    }
                }
            }
        }
    }
    println!("{read} .prf files read, {refused} of integrated intensities refused; the excluded regions of {regions_checked} against their .out");
    assert!(fails.is_empty(), "{}", fails.join("\n"));
    // L322: the joint runs pbb, scor_rxn and prbaco add 12 (one per pattern and pass). L334: P05's TOF runs add 25:
    // arg_si's two passes and seven variants (9), CeO2_PEARL (2), alum with alum_noq (3), lamn_pol with lamn_pol_irf
    // (3), nac-osiris (2), gabriela's three patterns in each pass (6).
    // L346: the runs of P41, P30 and P23 add 42 (rigid 22, texture 4, background_terms 16).
    // L351: the citric runs add 3 (citric-acid-sls_1's two passes, citric_core; citric_free wrote none).
    assert_eq!((read, refused), (288, 3));
}

#[test]
fn hkl_files() {
    let (mut read, mut counted, mut ticks) = (0, 0, 0);
    let mut fails = Vec::new();
    for f in all("hkl") {
        let name = f.display().to_string();
        let blocks = match read_hkl(&latin1(&f)) {
            Ok(b) => b,
            Err(e) => {
                fails.push(format!("{name}: {}", e.0));
                continue;
            }
        };
        read += 1;
        for b in &blocks {
            if let Some(n) = b.heading.split_once(" reflections").and_then(|(a, _)| a.split_whitespace().last()?.parse::<usize>().ok()) {
                counted += 1;
                // FullProf's heading: "(The # of eff. reflections may be lower)".
                if b.rows.len() > n || b.rows.is_empty() {
                    fails.push(format!("{name}: block {:?} announces {n} reflections, {} read", b.phase, b.rows.len()));
                }
            }
        }
        // A run's CODE_calc.hkl (Hkl = −3) and CODE_calc.prf: the 2T of the .hkl are ticks of the .prf (λ1).
        let prf = f.with_extension("prf");
        // NdNi_HTO_calc and cuf1k_calc, the two with SyCos ≠ 0, place their ticks with SyCos·cos 2θ (L267; with cos θ
        // they missed by 0.004°–0.164° and 0.0002°–0.0035°, which L098 had left unresolved).
        if f.to_string_lossy().contains("/runs/") && f.file_stem().unwrap().to_string_lossy().ends_with("_calc") && prf.is_file() {
            let Ok(p) = read_prf(&latin1(&prf)) else { continue };
            // A TOF pattern of Npr = 9 (P05, L341): the .hkl's column is Dtt1·d alone, and the tick Zero + Dtt1 d + Dtt2 d²
            // + Dtt_1overd/d, with the values after FullProf's cycle that the .prf's second line holds where a constant
            // wavelength's Zero, SyCos and SySin are (nac-osiris_calc.prf: −2.30209, 17591.13867, −0.67840; its .pcr's
            // zero, −2.279, is the one before), Dtt_1overd the .pcr's; Npr = 10 prints the tick itself.
            let pcr = drs::pcr::read_all(&latin1(&f.with_extension("pcr"))).ok().and_then(|v| v.into_iter().next());
            let tof = pcr.filter(|q| q.flag("Uni") == 1 && q.flag("Npr") == 9).map(|q| (p.zero, [p.sycos, p.sysin, q.dtt[2]]));
            for b in &blocks {
                let (Some(phase), Some(tt)) = (b.phase, b.column("2T")) else { continue };
                ticks += 1;
                for (h, t) in b.hkl().iter().zip(&tt) {
                    let th = (t / 2.0).to_radians();
                    let shifted = match tof {
                        Some((zero, [d1, d2, d3])) => {
                            let d = t / d1;
                            zero + t + d2 * d * d + d3 / d
                        }
                        None => t + p.zero + p.sycos * (2.0 * th).cos() + p.sysin * (2.0 * th).sin(),
                    };
                    // A TOF tick is printed in single precision (nac-osiris_calc.prf's (2 1 1) at 73607.4766 µs, where its
                    // step is 0.0078 µs): within 10⁻⁷ of it more.
                    let tol = if tof.is_some() { 2e-3 + 1e-7 * shifted.abs() } else { 2e-3 };
                    let hit = p.reflections.iter().any(|r| r.phase == phase && r.wavelength == 1 && r.hkl == *h && (r.two_theta - shifted).abs() < tol);
                    if !hit {
                        fails.push(format!("{name}: phase {phase} {h:?} at {t} is not a tick of {}", prf.display()));
                        break;
                    }
                }
            }
        }
    }
    println!("{read} .hkl files read; {counted} blocks against their announced count; {ticks} blocks against their .prf's ticks");
    assert!(fails.is_empty(), "{}", fails.join("\n"));
    assert_eq!(read, all("hkl").len());
}
