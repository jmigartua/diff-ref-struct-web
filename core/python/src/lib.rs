//! Python bindings to the diff-ref-struct core (P53, Q1/L021; L380): the native part of the `drs` package.
//!
//! No physics here: each function calls the core as the command line does, and gives its result to Python as JSON text,
//! which `drs/__init__.py` turns into dicts and lists. The command line itself travels in the package and is run by
//! `drs.run` (the same binary, so the same functions with the same output).

use drs::control::Control;
use drs::data::Pattern;
use drs::lsq::{Refinement, Statistics};
use drs::refinable::ParameterSet;
use pyo3::prelude::*;
use serde_json::{json, Value};
use std::path::{Path, PathBuf};

pyo3::create_exception!(_drs, DrsError, pyo3::exceptions::PyException);

fn err(e: impl std::fmt::Display) -> PyErr {
    DrsError::new_err(e.to_string())
}

/// JSON's numbers have no NaN or infinity: those become null.
fn num(v: f64) -> Value {
    if v.is_finite() {
        json!(v)
    } else {
        Value::Null
    }
}

fn nums(v: &[f64]) -> Value {
    Value::Array(v.iter().map(|&x| num(x)).collect())
}

fn stats(s: &Statistics) -> Value {
    json!({
        "n": s.n, "p": s.p, "c": s.c, "rp": num(s.rp), "rwp": num(s.rwp), "rexp": num(s.rexp), "chi2": num(s.chi2),
        "durbin_watson": num(s.durbin_watson),
        "phases": s.phases.iter().map(|&(rb, rf)| json!({"r_bragg": num(rb), "r_f": num(rf)})).collect::<Vec<_>>(),
    })
}

fn latin1(path: &Path) -> Result<String, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    Ok(String::from_utf8(bytes.clone()).unwrap_or_else(|_| bytes.iter().map(|&b| b as char).collect()))
}

/// The core's version (the package's).
#[pyfunction]
fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

/// A data file read by the core's readers (`formats::read_path`): its format and scans, each with its axis, x, y, σ.
#[pyfunction]
#[pyo3(signature = (path, format = None))]
fn read_data(path: &str, format: Option<&str>) -> PyResult<String> {
    let d = drs::formats::read_path(Path::new(path), format).map_err(err)?;
    let scans: Vec<Value> = d
        .scans
        .iter()
        .map(|s| json!({"label": s.label, "axis": format!("{:?}", s.axis), "x": nums(&s.pattern.x), "y": nums(&s.pattern.y), "sigma": nums(&s.pattern.sigma)}))
        .collect();
    Ok(json!({"format": d.format, "scans": scans}).to_string())
}

/// A refinement project: a control file (`.drc`) or a FullProf `.pcr` imported into the default model, with its data.
#[pyclass(unsendable)]
struct Project {
    control: Control,
    data: Option<Pattern>,
    /// The file's folder and stem (a `.pcr`'s export and results are named from them).
    dir: PathBuf,
    stem: String,
    refinement: Option<(ParameterSet, Refinement)>,
}

impl Project {
    fn data(&self) -> PyResult<&Pattern> {
        self.data.as_ref().ok_or_else(|| err("the project has no data file"))
    }

    fn parameter_set(&self) -> PyResult<ParameterSet> {
        self.control.refinement.build(&self.control.model).map_err(err)
    }
}

#[pymethods]
impl Project {
    /// Open a control file or a `.pcr` (its data found as `drs` finds them; `base`: the folder a control file's files
    /// are named from, as `--base`).
    #[new]
    #[pyo3(signature = (path, base = None))]
    fn new(path: &str, base: Option<&str>) -> PyResult<Self> {
        let p = Path::new(path);
        let dir = std::fs::canonicalize(p).ok().and_then(|q| q.parent().map(Path::to_path_buf)).unwrap_or_else(|| PathBuf::from("."));
        let stem = p.file_stem().and_then(|s| s.to_str()).unwrap_or("").to_string();
        if drs::control::is_control_file(p) {
            let o = match base {
                Some(b) => drs::control::open_in(p, Path::new(b)),
                None => drs::control::open(p),
            }
            .map_err(err)?;
            Ok(Project { control: o.control, data: o.data, dir, stem, refinement: None })
        } else {
            let pcr = drs::pcr::read(&latin1(p).map_err(err)?).map_err(|e| err(format!("{path}: {e}")))?;
            let control = Control::from_pcr(&pcr, &dir, &stem, &dir).map_err(err)?;
            let data = drs::powder::Model::from_pcr(&pcr, &dir, &stem).map_err(err)?.data;
            Ok(Project { control, data, dir, stem, refinement: None })
        }
    }

    fn title(&self) -> String {
        self.control.title.clone()
    }

    /// The control file's text (`Control::to_text`).
    fn control_text(&self) -> String {
        self.control.to_text()
    }

    /// Write the control file.
    fn save(&self, path: &str) -> PyResult<()> {
        std::fs::write(path, self.control.to_text()).map_err(|e| err(format!("{path}: {e}")))
    }

    /// Every value of the model with its address, dictionary id, value, limits and what fixes it (`Model::fields`).
    fn fields_json(&self) -> String {
        let f: Vec<Value> = self
            .control
            .model
            .fields()
            .into_iter()
            .map(|i| json!({"address": i.address, "id": i.id, "value": num(i.value), "lower": i.limits.0.map(num), "upper": i.limits.1.map(num), "fixed_by": i.fixed_by}))
            .collect();
        Value::Array(f).to_string()
    }

    /// A value by its address.
    fn get(&self, address: &str) -> PyResult<f64> {
        let f = drs::refinable::field_of(address, &self.control.model).map_err(err)?;
        Ok(self.control.model.get(f))
    }

    /// Set a value by its address, checked (`Model::set_field`: the values tied to it follow).
    fn set(&mut self, address: &str, value: f64) -> PyResult<()> {
        let f = drs::refinable::field_of(address, &self.control.model).map_err(err)?;
        self.control.model.set_field(f, value).map_err(err)
    }

    /// The free parameters of the refinement: name, value and the fields each moves with their coefficients.
    fn free_json(&self) -> PyResult<String> {
        let set = self.parameter_set()?;
        let m = &self.control.model;
        let f: Vec<Value> = set
            .free
            .iter()
            .map(|p| json!({"name": p.name, "value": num(p.value_in(m)), "fields": p.fields.iter().map(|&(f, c)| json!([drs::refinable::address(f, m), c])).collect::<Vec<_>>()}))
            .collect();
        Ok(Value::Array(f).to_string())
    }

    /// The pattern at the model's values: x, Yobs, σ, Ycal, the background, each phase's curve, and the statistics
    /// (χ² over N − P, P the refinement's free parameters).
    fn calculate_json(&self) -> PyResult<String> {
        let d = self.data()?;
        let c = self.control.model.calculate(&d.x);
        let np = self.parameter_set().map(|s| s.free.len()).unwrap_or(0);
        let s = drs::lsq::statistics(&self.control.model, d, &c.ycalc, np);
        Ok(json!({"x": nums(&d.x), "yobs": nums(&d.y), "sigma": nums(&d.sigma), "ycalc": nums(&c.ycalc), "background": nums(&c.background), "phases": c.per_phase.iter().map(|v| nums(v)).collect::<Vec<_>>(), "stats": stats(&s)}).to_string())
    }

    /// Refine: up to `cycles` (the file's, else 10) to convergence `eps` (the file's, else 0.1), the model moved to the
    /// values reached; the cycles run, whether it converged, why it stopped, the statistics and each free parameter's
    /// value and esd.
    #[pyo3(signature = (cycles = None, eps = None))]
    fn refine_json(&mut self, cycles: Option<usize>, eps: Option<f64>) -> PyResult<String> {
        let d = self.data()?.clone();
        let set = self.parameter_set()?;
        if set.free.is_empty() {
            return Err(err("no parameter is refined"));
        }
        let max_cycles = cycles.or(self.control.refinement.cycles).unwrap_or(10).max(1);
        let eps = eps.or(self.control.refinement.eps).unwrap_or(0.1);
        let r = drs::lsq::refine(&mut self.control.model, &d, &set, drs::lsq::Options { max_cycles, eps }, |_| {});
        let out = json!({
            "cycles": r.cycles, "converged": r.converged, "stop": r.stop, "stats": stats(&r.stats),
            "params": r.params.iter().map(|p| json!({"name": p.name, "value": num(p.value), "esd": num(p.esd)})).collect::<Vec<_>>(),
        });
        self.refinement = Some((set, r));
        Ok(out.to_string())
    }

    /// The statistics with free parameters at the values given, as a step sets them (`lsq::point`, E20a); the model
    /// is not changed.
    fn point_json(&self, values: Vec<(String, f64)>) -> PyResult<String> {
        let d = self.data()?;
        let set = self.parameter_set()?;
        let v: Vec<(&str, f64)> = values.iter().map(|(k, x)| (k.as_str(), *x)).collect();
        let p = drs::lsq::point(&self.control.model, d, &set, &v).map_err(err)?;
        Ok(json!({"ycalc": nums(&p.ycalc), "stats": stats(&p.statistics)}).to_string())
    }

    /// The checks of the control file as `drs check` makes them (`check::check_control`): severity, level, code,
    /// what, why, a suggestion and the line.
    fn check_json(&self) -> String {
        let d = drs::check::check_control(&self.control.to_text(), &self.dir, &drs::check::Options::default());
        let v: Vec<Value> = d
            .iter()
            .map(|x| json!({"severity": x.severity.name(), "level": x.level.name(), "code": x.code, "what": x.what, "why": x.why, "suggestion": x.suggestion, "addresses": x.addresses, "line": x.span.as_ref().map(|s| s.line)}))
            .collect();
        Value::Array(v).to_string()
    }

    /// The state's derived quantities (`state::summary`, E3, E18): per phase the cell with su's, the contacts and the
    /// valence sums; with the su's of the last refinement when there is one.
    fn state_json(&self) -> PyResult<String> {
        let d = self.data()?;
        let (set, r) = match &self.refinement {
            Some((s, r)) => (Some(s), Some(r)),
            None => (None, None),
        };
        let st = drs::state::summary(&self.control.model, d, set, r, &drs::state::Options::default());
        let phases: Vec<Value> = st
            .phases
            .iter()
            .map(|ph| {
                let contacts: Vec<Value> = ph
                    .contacts
                    .iter()
                    .enumerate()
                    .flat_map(|(i, cs)| cs.iter().map(move |c| (i, c)))
                    .map(|(i, c)| json!({"from": ph.labels[i], "to": ph.labels[c.atom], "distance": num(c.distance), "su": num(c.sigma), "plausible": c.plausible.as_ref().map(|(bv, ionic)| json!({"bond_valence": bv.iter().map(|(n, d)| json!([n, d.map(num)])).collect::<Vec<_>>(), "ionic_radii": num(*ionic)}))}))
                    .collect();
                let valence: Vec<Value> = ph.valence.as_ref().map_or(Vec::new(), |(v, _)| v.iter().map(|s| json!({"atom": ph.labels[s.atom], "sum": num(s.sum.0), "su": num(s.sum.1), "expected": num(s.expected), "plausible": [num(s.plausible.0), num(s.plausible.1)]})).collect());
                json!({"name": ph.name, "cell": ph.cell_su.iter().map(|&(v, e)| json!([num(v), num(e)])).collect::<Vec<_>>(), "contacts": contacts, "valence": valence})
            })
            .collect();
        Ok(json!({"stats": stats(&st.stats), "phases": phases}).to_string())
    }

    /// The `.pcr` FullProf would read (`pcr_export::export`), and the files to write beside it.
    fn pcr_json(&self) -> PyResult<String> {
        let e = drs::control::pcr_export::export(&self.control, &self.stem, self.data.as_ref()).map_err(err)?;
        Ok(json!({"pcr": e.pcr, "files": e.files, "notes": e.notes}).to_string())
    }

    /// The results file (R28, pdCIF) of the last refinement, or of the model's values before any (`results::write`),
    /// named `name`.
    fn results(&self, name: &str) -> PyResult<String> {
        let d = self.data()?;
        let (set, r) = match &self.refinement {
            Some((s, r)) => (Some(s), Some(r)),
            None => (None, None),
        };
        Ok(drs::results::write(&self.control, d, set, r, name).map_err(err)?.0)
    }
}

#[pymodule]
fn _drs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add_function(wrap_pyfunction!(read_data, m)?)?;
    m.add_class::<Project>()?;
    m.add("DrsError", m.py().get_type::<DrsError>())?;
    Ok(())
}
