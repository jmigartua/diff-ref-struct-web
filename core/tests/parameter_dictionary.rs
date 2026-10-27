//! The parameter dictionary against the core (test `parameter_dictionary`, L044, L046).
//!
//! Checks that:
//! - every field of the model's structures (`powder::Model`, `powder::Phase`, `structure_factor::Atom`, read from
//!   their source) is listed in `parameters::FIELDS`, with dictionary ids or the reason it holds none;
//! - every id there is in the dictionary, and every refinable or read parameter of the dictionary is held by a field;
//! - the dictionary's own `core` lists (in `parameter_table`, generated from `ecosystem/data/parameters.yml`) say the
//!   same as `FIELDS`;
//! - the generated table holds exactly the ids of `ecosystem/data/parameters.yml`, in its order (build.py was run),
//!   with each entry's name, group and kind, and the groups (L058: the whole dictionary is in the core).

use drs::parameters::{lookup, FIELDS, GROUPS, PARAMETERS};
use std::collections::BTreeSet;

fn root(rel: &str) -> String {
    format!("{}/{rel}", env!("CARGO_MANIFEST_DIR"))
}

/// The field names of `pub struct <name> { … }` in a source file.
fn fields(file: &str, name: &str) -> Vec<String> {
    let src = std::fs::read_to_string(root(file)).unwrap();
    let start = src.find(&format!("pub struct {name} {{")).unwrap_or_else(|| panic!("{name} not in {file}"));
    let body = &src[start..start + src[start..].find("\n}").unwrap()];
    body.lines().skip(1).filter_map(|l| l.trim().strip_prefix("pub ")).filter_map(|l| l.split_once(':').map(|(f, _)| f.trim().to_string())).collect()
}

#[test]
fn parameter_dictionary() {
    let listed: BTreeSet<&str> = FIELDS.iter().map(|f| f.field).collect();
    assert_eq!(listed.len(), FIELDS.len(), "a field listed twice in parameters::FIELDS");
    let mut n = 0;
    for (file, name) in [("src/powder.rs", "Model"), ("src/powder.rs", "Phase"), ("src/structure_factor.rs", "Atom")] {
        for f in fields(file, name) {
            let key = format!("{name}.{f}");
            assert!(listed.contains(key.as_str()), "{key} is not in parameters::FIELDS: give its dictionary ids or say why it holds none");
            n += 1;
        }
    }
    assert_eq!(n, FIELDS.len(), "parameters::FIELDS lists a field the structures no longer have");
    for f in FIELDS {
        assert!(f.ids.is_empty() != f.not_a_parameter.is_empty(), "{}: ids or a reason, not both or neither", f.field);
        for id in f.ids {
            let p = lookup(id).unwrap_or_else(|| panic!("{}: '{id}' is not in the dictionary", f.field));
            assert!(p.core.contains(&f.field), "{id}: the dictionary does not name {} among its core fields", f.field);
        }
    }
    for p in PARAMETERS {
        assert!(!p.core.is_empty(), "{}: no core field", p.id);
        for c in p.core {
            let f = FIELDS.iter().find(|f| f.field == *c).unwrap_or_else(|| panic!("{}: core field {c} is not in FIELDS", p.id));
            assert!(f.ids.contains(&p.id), "{}: {c} does not hold it in FIELDS", p.id);
        }
    }
    // The generated table holds the ids of the YAML file, in its order.
    let yml = std::fs::read_to_string(root("../ecosystem/data/parameters.yml")).unwrap();
    let ids: Vec<&str> = yml.lines().filter_map(|l| l.strip_prefix("  - id: ")).map(str::trim).collect();
    let table: Vec<&str> = PARAMETERS.iter().map(|p| p.id).collect();
    assert_eq!(ids, table, "core/src/parameter_table.rs is out of date: run ecosystem/tools/build.py");
    // Every entry's name, group and kind as the YAML gives them (L058), and the groups with their titles.
    let mut entry = "";
    for l in yml.lines() {
        if let Some(id) = l.strip_prefix("  - id: ") {
            entry = id.trim();
            continue;
        }
        let Some(p) = lookup(entry) else { continue };
        let field = |k: &str| l.strip_prefix(&format!("    {k}: ")).map(|v| v.trim().trim_matches('"').to_string());
        if let Some(v) = field("name") {
            assert_eq!(p.name, v, "{entry}: name");
        }
        if let Some(v) = field("group") {
            assert_eq!(p.group, v, "{entry}: group");
            assert!(GROUPS.iter().any(|g| g.0 == v), "{entry}: group {v} is not in GROUPS");
        }
        if let Some(v) = field("kind") {
            assert_eq!(p.kind, v, "{entry}: kind");
            assert_eq!(p.refinable, v == "parameter", "{entry}: refinable");
        }
    }
    for p in PARAMETERS {
        assert!(!p.name.is_empty() && !p.meaning.is_empty() && !p.equation.is_empty() && !p.range.is_empty(), "{}: an empty field", p.id);
    }
    let keys: Vec<&str> = yml.lines().filter_map(|l| l.trim().strip_prefix("- {key: ")).map(|r| r.split(',').next().unwrap()).collect();
    assert_eq!(keys, GROUPS.iter().map(|g| g.0).collect::<Vec<_>>(), "GROUPS");
    assert_eq!(drs::parameters::dictionary().len(), PARAMETERS.len());
    println!("{} fields of Model, Phase and Atom; {} dictionary entries, {} refinable", FIELDS.len(), PARAMETERS.len(), PARAMETERS.iter().filter(|p| p.refinable).count());
}

/// The engine's names and limits against the dictionary (L053): every field a refinement can move, on runs with each
/// background model and profile, is named by a refinable id of the dictionary, and the dictionary's `limits` say
/// "None." where `refinable::limits` never sets a bound, a bound where it sets one, and "none" for the cases left
/// unbounded by a conditional limit.
#[test]
fn engine_names_and_limits() {
    use drs::jacobian::Field;
    use drs::powder::{Intensities, Model};
    let runs = [("Y2O3", "y2o3_structure"), ("rutile-anatase", "rutana_tch"), ("LaB6", "LaB6rx"), ("Yb2ZrWO8", "ybfluo_Riet_IRF"), ("Yb2ZrWO8-PO-March", "ybfluo_Riet_IRF")];
    let mut fails = Vec::new();
    let mut seen = BTreeSet::new();
    for (run, code) in runs {
        let d = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs").join(run);
        let text: String = std::fs::read(d.join(format!("{code}.pcr"))).unwrap().iter().map(|&b| b as char).collect();
        let p = drs::pcr::read(&text).unwrap();
        let m = Model::from_pcr(&p, &d, code).unwrap().model;
        let mut fl = vec![Field::Zero, Field::Sycos, Field::Sysin, Field::Ratio, Field::Background(0)];
        for (k, ph) in m.phases.iter().enumerate() {
            fl.extend([Field::Scale(k), Field::Bov(k), Field::U(k), Field::V(k), Field::W(k), Field::X(k), Field::Y(k), Field::GaussSize(k), Field::Eta0(k)]);
            fl.extend((0..6).map(|c| Field::Cell(k, c)));
            fl.extend([Field::Asym(k, 0), Field::SL(k), Field::DL(k), Field::Pref(k, 0), Field::Pref(k, 1)]);
            if let Intensities::Structure(..) = ph.intensities {
                fl.extend([Field::AtomX(k, 0, 0), Field::AtomB(k, 0), Field::AtomOcc(k, 0), Field::AtomBeta(k, 0, 0)]);
            }
        }
        for f in fl {
            let id = drs::refinable::dictionary_id(f, &m);
            seen.insert(id);
            let Some(e) = PARAMETERS.iter().find(|e| e.id == id) else {
                fails.push(format!("{run}: {f:?} is named {id}, not in the dictionary"));
                continue;
            };
            if !e.refinable {
                fails.push(format!("{run}: {f:?} is named {id}, which the dictionary does not mark refinable"));
            }
            let bounded = drs::refinable::limits(f, &m) != (None, None);
            // A conditional limit ("≥ 0 for …; none otherwise") names the unbounded case.
            let ok = if bounded { e.limits != "None." } else { e.limits == "None." || e.limits.contains("none") };
            if !ok {
                fails.push(format!("{run}: {f:?} ({id}): the engine's limits {:?} against the dictionary's \"{}\"", drs::refinable::limits(f, &m), e.limits));
            }
        }
    }
    // The class and the largest shift per cycle (W12, L107–L108): every field's dictionary class is the engine's
    // (`steering::class_of`), and its max_shift a positive number or "none".
    for run in ["Y2O3/y2o3_structure_conv.pcr", "Egypt/Egypt_qa_Rietveld_conv.pcr", "NdNi-HTO/NdNi_HTO_conv.pcr", "modes/A_modes/a_modes.pcr"] {
        let d = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs");
        let path = d.join(run);
        let text: String = std::fs::read(&path).unwrap().iter().map(|&b| b as char).collect();
        let mut p = drs::pcr::read(&text).unwrap();
        p.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
        let stem = path.file_stem().unwrap().to_str().unwrap().to_string();
        let m = drs::powder::Model::from_pcr(&p, path.parent().unwrap(), &stem).unwrap().model;
        for f in m.all_fields() {
            let id = drs::refinable::dictionary_id(f, &m);
            let e = PARAMETERS.iter().find(|e| e.id == id).unwrap();
            if e.class != drs::steering::class_of(f).word() {
                fails.push(format!("{run}: {id}: class {} in the dictionary, {} in the engine", e.class, drs::steering::class_of(f).word()));
            }
        }
    }
    // A TOF pattern's fields (P05), with Hodges' terms (gabriela's first pattern, Npr = 10) and without (arg_si): named
    // by refinable ids, the dictionary's limits and class the engine's.
    let base = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../manuals/fullprof/runs");
    for (run, k) in [("arg_si/arg_si.pcr", 0), ("gabriela/gabriela.pcr", 0)] {
        let path = base.join(run);
        let text: String = std::fs::read(&path).unwrap().iter().map(|&b| b as char).collect();
        let p = drs::pcr::read_all(&text).unwrap().remove(k);
        let stem = path.file_stem().unwrap().to_str().unwrap().to_string();
        let m = drs::powder::Model::from_pcr(&p, path.parent().unwrap(), &stem).unwrap().model;
        let fields = m.all_fields();
        assert!(fields.iter().any(|f| matches!(f, Field::Tof(_))) && fields.iter().any(|f| matches!(f, Field::TofProfile(..))), "{run}: no TOF fields");
        for f in fields {
            let id = drs::refinable::dictionary_id(f, &m);
            seen.insert(id);
            let Some(e) = PARAMETERS.iter().find(|e| e.id == id) else {
                fails.push(format!("{run}: {f:?} is named {id}, not in the dictionary"));
                continue;
            };
            if !e.refinable {
                fails.push(format!("{run}: {f:?} is named {id}, which the dictionary does not mark refinable"));
            }
            let bounded = drs::refinable::limits(f, &m) != (None, None);
            if !(if bounded { e.limits != "None." } else { e.limits == "None." || e.limits.contains("none") }) {
                fails.push(format!("{run}: {f:?} ({id}): the engine's limits {:?} against the dictionary's \"{}\"", drs::refinable::limits(f, &m), e.limits));
            }
            if e.class != drs::steering::class_of(f).word() {
                fails.push(format!("{run}: {id}: class {} in the dictionary, {} in the engine", e.class, drs::steering::class_of(f).word()));
            }
        }
    }
    for e in PARAMETERS.iter().filter(|e| e.refinable) {
        if drs::steering::Class::from_word(e.class).is_none() || e.class_source.is_empty() {
            fails.push(format!("{}: class \"{}\" ({})", e.id, e.class, e.class_source));
        }
        if !(e.max_shift == "none" || e.max_shift.parse::<f64>().is_ok_and(|v| v > 0.0)) {
            fails.push(format!("{}: max_shift \"{}\"", e.id, e.max_shift));
        }
    }
    // Every refinable entry carries its limits.
    for e in PARAMETERS.iter().filter(|e| e.refinable && e.limits.is_empty()) {
        fails.push(format!("{}: no limits in the dictionary", e.id));
    }
    println!("{} dictionary ids named by the engine: {}", seen.len(), seen.iter().copied().collect::<Vec<_>>().join(", "));
    assert!(fails.is_empty(), "{}", fails.join("\n"));
}
