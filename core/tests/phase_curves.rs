//! Each phase's calculated curve (the app's request 7, L070).
//! - `phases_sum_to_the_pattern`: in every multi-phase run, the background plus the phases' curves is Ycalc.
//!
//! FullProf writes each phase's curve with the background in `.sub` files (Ipr = 3); the only ones in the material are
//! NdNi-HTO's, a run the core does not import (magnetic phase), whose nuclear phase the core calculates with its peaks
//! shifted from FullProf's (L070). The verified multi-phase runs have no `.sub` file yet.

use drs::powder::Model;
use std::path::Path;

fn latin1(p: &Path) -> String {
    std::fs::read(p).unwrap().iter().map(|&b| b as char).collect()
}

#[test]
fn phases_sum_to_the_pattern() {
    let base = Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs");
    let mut checked = 0;
    for run in std::fs::read_dir(&base).unwrap().map(|e| e.unwrap().path()).filter(|p| p.is_dir()) {
        for f in std::fs::read_dir(&run).unwrap().map(|e| e.unwrap().path()) {
            let name = f.file_name().unwrap().to_string_lossy().to_string();
            if !name.ends_with(".pcr") || name.contains("_calc") || name.contains("_conv") {
                continue;
            }
            let stem = name.trim_end_matches(".pcr");
            // A joint run's patterns each (P06, L322).
            // A file the reader refuses is left out (none since P30's step read runs/texture, L360).
            for mut p in drs::pcr::read_all(&latin1(&f)).unwrap_or_default() {
                p.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
                let Ok(l) = Model::from_pcr(&p, &run, stem) else { continue };
                let Some(data) = &l.data else { continue };
                let c = l.model.calculate(&data.x);
                assert_eq!(c.per_phase.len(), l.model.phases.len());
                for i in 0..data.x.len() {
                    let sum = c.background[i] + c.per_phase.iter().map(|v| v[i]).sum::<f64>();
                    assert!((sum - c.ycalc[i]).abs() <= 1e-12 * c.ycalc[i].abs().max(1.0), "{name} at {}: {sum} vs {}", data.x[i], c.ycalc[i]);
                }
                // A phase's curve is the pattern of that phase alone.
                for k in 0..l.model.phases.len() {
                    let mut alone = l.model.clone();
                    for (j, ph) in alone.phases.iter_mut().enumerate() {
                        if j != k {
                            ph.scale = 0.0;
                        }
                    }
                    let a = alone.calculate(&data.x);
                    for i in 0..data.x.len() {
                        let v = a.ycalc[i] - a.background[i];
                        assert!((v - c.per_phase[k][i]).abs() <= 1e-9 * a.ycalc[i].abs().max(1.0), "{name} phase {k} at {}", data.x[i]);
                    }
                }
                if l.model.phases.len() > 1 {
                    checked += 1;
                }
            }
        }
    }
    println!("{checked} multi-phase runs: the phases' curves sum to Ycalc, each the pattern of its phase alone");
    assert!(checked >= 2);
}
