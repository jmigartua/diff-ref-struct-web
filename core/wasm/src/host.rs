//! The Plots module's host commands for the browser (W1, I10, R34; L424): the app's Tauri commands
//! (`app/src-tauri/src/modules/plots.rs`, app 0.53.0) ported to run beside the core in WebAssembly, with the same names,
//! the same arguments and the same JSON, so that the module's `invoke(cmd, args)` is served by `dispatch`. Every number
//! is the core's (L139): this file keeps the sets by id, maps the window's words to the core's types, and names
//! things, as the app's does.
//!
//! What differs from the app's, because a browser has no paths: files come by name and content (`plots_open_from`,
//! `plots_open_buffer_from`, `plots_add_from`, `plots_open_set_from`, `plots_temperatures_from`); what the app writes
//! to the project is returned as text for the host to save (`plots_save`, `plots_derive`); Data's operations (kept in
//! the project's `.drs/data/`) and a control file's λ (read from the project) are not read, and a λ is typed instead.
//! The port is a copy of the app's code where nothing touches the disk, kept so to be compared line for line.

use drs::formats::Axis;
use drs::plots::{self, Derivation, Displacement, GridRange, GridStep, LambdaFrom, Normalisation, OpenedFile, PlotSet, XQuantity, YScale};
use serde::Serialize;
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

/// The sets the window shows, by id. A set's members name their files as given (relative to the folder the host read
/// them from).
#[derive(Default)]
pub struct Sets {
    inner: Mutex<(u32, HashMap<u32, Kept>)>,
}

struct Kept {
    set: PlotSet,
    base: PathBuf,
    /// How the window shows the set: a normalisation and an x quantity (PL6, PL10).
    rule: ViewRule,
    /// The step the members were rebinned to on opening (W9), kept for the set file.
    rebin: Option<f64>,
    /// The set as measured (cut to the analysis's interval, not rebinned), kept the first time `rebin` changes the
    /// step: the browser has no files to read again (the app reads them, `plots_rebin`).
    raw: Option<PlotSet>,
}

impl Sets {
    fn insert(&self, set: PlotSet, base: PathBuf) -> u32 {
        let mut g = self.inner.lock().unwrap();
        g.0 += 1;
        let id = g.0;
        g.1.insert(id, Kept { set, base, rule: ViewRule::default(), rebin: None, raw: None });
        id
    }
    fn with<T>(&self, id: u32, f: impl FnOnce(&mut Kept) -> Result<T, String>) -> Result<T, String> {
        let mut g = self.inner.lock().unwrap();
        let k = g.1.get_mut(&id).ok_or(format!("no plot set {id}"))?;
        f(k)
    }
    pub fn close(&self, id: u32) {
        self.inner.lock().unwrap().1.remove(&id);
    }
}

// ------------------------------------------------------------------------------------------------------- the views

/// A member's value of the variable and where it came from, in the core's words (`Found::words`).
#[derive(Serialize, Clone, Debug)]
pub struct VariableView {
    pub value: f64,
    pub source: String,
    pub typed: bool,
}

/// One member as the window lists it: what the reader says (`plots::Facts`), its value of the variable with its source,
/// and every source that gave one (so that a disagreement is seen).
#[derive(Serialize, Clone, Debug)]
pub struct MemberView {
    /// The member's identity across sorts: its file and scan.
    pub key: String,
    pub index: usize,
    pub file: String,
    /// The file on disk (the set's folder joined to the name), for Open in Data.
    pub path: String,
    pub label: String,
    pub scan: usize,
    pub scans: usize,
    pub comment: Option<String>,
    pub format: String,
    pub format_name: String,
    pub axis: String,
    pub points: usize,
    pub range: (f64, f64),
    pub step: Option<f64>,
    pub wavelength: Option<f64>,
    pub title: Option<String>,
    pub monitor: Option<f64>,
    pub counting_time: Option<f64>,
    /// The file's modification time, seconds since the epoch (W4's sort by time): the host's, given with the file.
    pub modified: Option<f64>,
    pub variable: Option<VariableView>,
    pub candidates: Vec<VariableView>,
    /// The operations Data keeps for the member's file, in Data's words (PL20): the member is drawn after them.
    pub operations: Vec<String>,
    /// The member's other variables with their sources (W3, core 0.79.0; Plots' step 21): the map's right-hand axis.
    pub others: Vec<OtherView>,
    /// The member's factor and λ in the view, or why the view leaves it out (PL6, PL10).
    pub shown: ShownMember,
}

/// A member's other variable (W3): its name, value and source in words.
#[derive(Serialize, Clone, Debug)]
pub struct OtherView {
    pub name: String,
    pub value: f64,
    pub source: String,
}

/// One member's scan as the view shows it (W9 (a): asked for by the views that draw it, not carried by every SetView;
/// the app's `PatternView`): the member by its index in the set, the reader's scan after Data's operations,
/// normalised and in the x quantity chosen (x, y, σ; None for a point with no weight).
#[derive(Serialize, Clone, Debug)]
pub struct PatternView {
    pub member: usize,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub sigma: Vec<Option<f64>>,
}

/// The members' scans as the view shows them (W9 (a), the app's `patterns`): those named by their index in the set,
/// else all; a member the view leaves out has none.
pub fn patterns(sets: &Sets, id: u32, members: Option<Vec<usize>>) -> Result<Vec<PatternView>, String> {
    sets.with(id, |k| {
        let s = shown(k, true)?;
        Ok(s.index
            .iter()
            .enumerate()
            .filter(|(_, &i)| members.as_ref().is_none_or(|ms| ms.contains(&i)))
            .map(|(q, &i)| {
                let p = &s.set.members[q].pattern;
                PatternView { member: i, x: p.x.clone(), y: p.y.clone(), sigma: p.sigma.iter().map(|&e| e.is_finite().then_some(e)).collect() }
            })
            .collect())
    })
}

/// The set rebinned to a step (W9 (c), the app's `rebin`; Plots' step 22), or as measured (none): from the set as
/// measured, kept the first time (the browser cannot read the files again); a set opened rebinned has no measured
/// scans kept, and changing its step is refused.
pub fn rebin(sets: &Sets, id: u32, step: Option<f64>) -> Result<SetView, String> {
    if let Some(s) = step {
        if !(s.is_finite() && s > 0.0) {
            return Err(format!("a rebin step of {s}: give a positive number"));
        }
    }
    sets.with(id, |k| {
        if k.raw.is_none() {
            if let Some(r) = k.rebin {
                return Err(format!("the set was opened rebinned to {r}: open its files again to choose another step"));
            }
            k.raw = Some(k.set.clone());
        }
        let mut set = k.raw.clone().unwrap();
        if let Some(s) = step {
            set.rebin(s)?;
        }
        k.set = set;
        k.rebin = step;
        Ok(view(id, k))
    })
}

/// A set as the window shows it.
#[derive(Serialize, Clone, Debug)]
pub struct SetView {
    pub id: u32,
    pub title: Option<String>,
    pub variable_name: String,
    pub variable_unit: String,
    pub members: Vec<MemberView>,
    /// Whether every member has the same points (`plots::same_grid`).
    pub same_grid: bool,
    /// The members without a value of the variable, drawn apart until one is typed (L169).
    pub apart: Vec<usize>,
    /// The folder the members' names are relative to.
    pub base: String,
    /// The normalisation shown (PL6), in the window's words, and its reference member's index (by reference); None for
    /// the counts as read.
    pub normalisation: Option<NormChoice>,
    pub norm_reference: Option<usize>,
    /// The x quantity shown (PL10): its word and its label (the core's `XQuantity`); None for the files' own x.
    pub x_quantity: Option<String>,
    pub x_label: String,
    /// The λ given for the members whose files state none, and where it comes from (a control file, typed).
    pub x_fallback: Option<(f64, String)>,
    /// The members the x quantity cannot show for want of a λ (`XView::without`).
    pub without_lambda: Vec<usize>,
    /// Why the view could not be made as chosen (the members are then drawn as read).
    pub view_error: Option<String>,
}

fn variable_view(v: &plots::Variable) -> VariableView {
    VariableView { value: v.value, source: v.found.words(), typed: matches!(v.found, plots::Found::Typed(_)) }
}

/// The operations Data keeps for a file: none in the browser (they live in the project's `.drs/data/`).
fn data_ops_words(_base: &Path, _file: &str) -> Vec<String> {
    Vec::new()
}

fn view(id: u32, k: &Kept) -> SetView {
    let s = &k.set;
    let (sh, view_error) = match shown(k, true) {
        Ok(v) => (Some(v), None),
        Err(e) => (None, Some(e)),
    };
    let members = s
        .members
        .iter()
        .enumerate()
        .map(|(i, m)| MemberView {
            key: format!("{}#{}", m.file, m.scan),
            index: i,
            file: m.file.clone(),
            path: k.base.join(&m.file).to_string_lossy().into_owned(),
            label: m.label.clone(),
            scan: m.scan,
            scans: m.scans,
            comment: m.comment.clone(),
            format: m.facts.format.clone(),
            format_name: drs::formats::by_id(&m.facts.format).map_or(m.facts.format.clone(), |f| f.name.to_string()),
            axis: m.facts.axis.label().to_string(),
            points: m.facts.points,
            range: m.facts.range,
            step: m.facts.step,
            wavelength: m.facts.wavelength,
            title: m.facts.title.clone(),
            monitor: m.facts.monitor,
            counting_time: m.facts.counting_time,
            modified: m.modified,
            variable: m.variable.as_ref().map(variable_view),
            candidates: m.candidates.iter().map(variable_view).collect(),
            operations: data_ops_words(&k.base, &m.file),
            others: m.others.iter().map(|o| OtherView { name: o.name.clone(), value: o.value, source: o.found.words() }).collect(),
            shown: sh.as_ref().map(|v| v.members[i].clone()).unwrap_or_default(),
        })
        .collect();
    let x = k.rule.x.as_ref();
    SetView {
        id,
        title: s.title.clone(),
        variable_name: s.variable_name.clone(),
        variable_unit: s.variable_unit.clone(),
        members,
        same_grid: plots::same_grid(s),
        apart: s.without_variable(),
        base: k.base.to_string_lossy().into_owned(),
        normalisation: k.rule.norm.clone(),
        norm_reference: k.rule.norm_ref.as_ref().and_then(|r| s.members.iter().position(|m| key_of(m) == *r)),
        x_quantity: x.map(|r| r.quantity.word().to_string()),
        x_label: x.map_or_else(|| s.members.first().map_or(Axis::TwoTheta, |m| m.facts.axis).label().to_string(), |r| r.quantity.label().to_string()),
        x_fallback: x.and_then(|r| r.fallback.as_ref().map(|f| (f.0, lambda_words(&f.1)))),
        without_lambda: sh.as_ref().map(|v| v.without.clone()).unwrap_or_default(),
        view_error,
    }
}

// ------------------------------------------------------------------------------- the view's rule (PL6, PL10)

/// A normalisation as the window chooses it (PL6): `by` none, monitor, time, maximum (over [from, until] when they are
/// given), area (the same area over [from, until]), peak (the integral over [from, until] above its `background`,
/// "auto" or "line") or reference (member `member`'s area over [from, until]); `to` the value each member is brought to
/// (1 when not given).
#[derive(serde::Deserialize, Serialize, Clone, Debug, Default, PartialEq)]
pub struct NormChoice {
    pub by: String,
    pub to: Option<f64>,
    pub from: Option<f64>,
    pub until: Option<f64>,
    pub member: Option<usize>,
    pub background: Option<String>,
}

/// An x quantity as the window chooses it (PL10): `quantity` 2theta, d, q, stl or s; a λ typed for some members (by
/// index, with the user's note); for the others whose files state none, a control file's λ1 (its path) or a λ typed.
#[derive(serde::Deserialize, Serialize, Clone, Debug, Default)]
pub struct XChoice {
    pub quantity: String,
    #[serde(default)]
    pub given: Vec<(usize, f64, String)>,
    pub control: Option<String>,
    pub typed: Option<f64>,
}

/// How the window shows a set (PL6, PL10), kept with it; members named by their file and scan, so that a reordering or
/// an added file leaves the rule as chosen.
#[derive(Clone, Debug, Default)]
struct ViewRule {
    norm: Option<NormChoice>,
    /// The reference member of a normalisation by reference.
    norm_ref: Option<String>,
    x: Option<XRule>,
}

#[derive(Clone, Debug)]
struct XRule {
    quantity: XQuantity,
    given: Vec<(String, f64, LambdaFrom)>,
    fallback: Option<(f64, LambdaFrom)>,
}

/// One member as the view shows it (PL6, PL10): its factor, its λ and where it comes from, or why the view leaves it
/// out (the core's `why_not`).
#[derive(Serialize, Clone, Debug, Default)]
pub struct ShownMember {
    pub factor: Option<f64>,
    pub lambda: Option<f64>,
    pub lambda_from: Option<String>,
    pub left_out: Option<String>,
}

/// The set as the view shows it: its members normalised and in the x quantity chosen, those the view leaves out
/// removed; `index` gives each shown member's index in the set, `members` what the view did to each member of the set.
struct Shown {
    set: PlotSet,
    index: Vec<usize>,
    members: Vec<ShownMember>,
    without: Vec<usize>,
}

impl Shown {
    /// A member of the set, by its place among the members shown.
    fn pos(&self, member: usize, what: &str) -> Result<usize, String> {
        self.index.iter().position(|&i| i == member).ok_or_else(|| format!("{what} is left out of the view: {}", self.members.get(member).and_then(|m| m.left_out.clone()).unwrap_or_else(|| "no such member".into())))
    }
    /// Each member of the set with its item among those shown, None for a member left out.
    fn spread<T>(&self, items: Vec<T>) -> Vec<Option<T>> {
        let mut out: Vec<Option<T>> = (0..self.members.len()).map(|_| None).collect();
        for (p, t) in items.into_iter().enumerate() {
            out[self.index[p]] = Some(t);
        }
        out
    }
}

fn key_of(m: &plots::Member) -> String {
    format!("{}#{}", m.file, m.scan)
}

fn lambda_words(f: &LambdaFrom) -> String {
    match f {
        LambdaFrom::File => "the file".into(),
        LambdaFrom::ControlFile(c) => format!("the control file {c}"),
        LambdaFrom::Typed(n) => format!("typed ({n})"),
    }
}

/// The window's normalisation as the core's (`plots::Normalisation`), the reference member found by its key.
fn normalisation(c: &NormChoice, set: &PlotSet, reference: Option<&str>) -> Result<Normalisation, String> {
    let to = c.to.unwrap_or(1.0);
    if !(to.is_finite() && to > 0.0) {
        return Err(format!("normalised to {to}: give a positive number"));
    }
    let range = |what: &str| match (c.from, c.until) {
        (Some(a), Some(b)) if a.is_finite() && b.is_finite() && b > a => Ok((a, b)),
        _ => Err(format!("{what} needs a range: from, and until above it")),
    };
    Ok(match c.by.as_str() {
        "none" => Normalisation::None,
        "monitor" => Normalisation::Monitor { to },
        "time" => Normalisation::Time { to },
        "maximum" if c.from.is_none() && c.until.is_none() => Normalisation::Maximum { to },
        "maximum" => {
            let (from, until) = range("the maximum in a range")?;
            Normalisation::MaximumIn { from, until, to }
        }
        "area" => {
            let (from, until) = range("the same area")?;
            Normalisation::Area { from, until, to }
        }
        "peak" => {
            let (from, until) = range("a peak's integral")?;
            let background = match c.background.as_deref() {
                None | Some("auto") => plots::WindowBackground::Auto { log_lambda: None },
                Some("line") => plots::WindowBackground::Line,
                Some(o) => return Err(format!("no background \"{o}\" (auto, line)")),
            };
            Normalisation::Peak { window: plots::Window { name: "peak".into(), from, to: until }, background, to }
        }
        "reference" => {
            let (from, until) = range("a reference member's area")?;
            let key = reference.ok_or("a normalisation by a reference member needs the member: choose one")?;
            let member = set.members.iter().position(|m| key_of(m) == key).ok_or(format!("the reference member {key} is no longer in the set"))?;
            Normalisation::Reference { member, from, until }
        }
        o => return Err(format!("no normalisation \"{o}\" (none, monitor, time, maximum, area, peak, reference)")),
    })
}

/// A set file's normalisation (PL27) as the window's choice, in its words, with the reference member by its key; none for
/// `Normalisation::None` or a background by points, which the window does not offer.
fn choice_of(n: &Normalisation, set: &PlotSet) -> (Option<NormChoice>, Option<String>) {
    let c = |by: &str, to: Option<f64>, from: Option<f64>, until: Option<f64>| NormChoice { by: by.into(), to, from, until, member: None, background: None };
    match n {
        Normalisation::None => (None, None),
        Normalisation::Monitor { to } => (Some(c("monitor", Some(*to), None, None)), None),
        Normalisation::Time { to } => (Some(c("time", Some(*to), None, None)), None),
        Normalisation::Maximum { to } => (Some(c("maximum", Some(*to), None, None)), None),
        Normalisation::MaximumIn { from, until, to } => (Some(c("maximum", Some(*to), Some(*from), Some(*until))), None),
        Normalisation::Area { from, until, to } => (Some(c("area", Some(*to), Some(*from), Some(*until))), None),
        Normalisation::Peak { window, background, to } => match background {
            plots::WindowBackground::Auto { .. } => (Some(NormChoice { background: Some("auto".into()), ..c("peak", Some(*to), Some(window.from), Some(window.to)) }), None),
            plots::WindowBackground::Line => (Some(NormChoice { background: Some("line".into()), ..c("peak", Some(*to), Some(window.from), Some(window.to)) }), None),
            plots::WindowBackground::Points { .. } => (None, None),
        },
        Normalisation::Reference { member, from, until } => match set.members.get(*member) {
            Some(m) => (Some(NormChoice { member: Some(*member), ..c("reference", None, Some(*from), Some(*until)) }), Some(key_of(m))),
            None => (None, None),
        },
    }
}

/// The window's x quantity as the view's rule: the λs typed for members by their keys, and the λ for the members whose
/// files state none, a control file's λ1 (read as the project reads it, `project::wavelengths`) or one typed.
fn x_rule(c: &XChoice, set: &PlotSet) -> Result<XRule, String> {
    let quantity = XQuantity::from_word(&c.quantity).ok_or(format!("no x quantity \"{}\" (2theta, d, q, stl, s)", c.quantity))?;
    let positive = |l: f64| if l.is_finite() && l > 0.0 { Ok(l) } else { Err(format!("a wavelength of {l} Å: give a positive number")) };
    let given = c
        .given
        .iter()
        .map(|(i, l, note)| {
            let m = set.members.get(*i).ok_or(format!("no member {}", i + 1))?;
            Ok((key_of(m), positive(*l)?, LambdaFrom::Typed(note.clone())))
        })
        .collect::<Result<Vec<_>, String>>()?;
    let fallback = match (&c.control, c.typed) {
        // a control file's λ1 is read from the project, which a browser does not have: typed instead
        (Some(path), _) => return Err(format!("{path}: a control file's λ is not read in the browser; type the λ")),
        (None, Some(l)) => Some((positive(l)?, LambdaFrom::Typed("for the members whose files state none".into()))),
        (None, None) => None,
    };
    Ok(XRule { quantity, given, fallback })
}

/// The set as the view shows it (PL6, PL10), the files untouched: each member normalised by the core
/// (`plots::normalise`), then, with `with_x`, shown in the x quantity chosen (`plots::x_view`); a member the core cannot
/// normalise or convert is left out, with its reason.
fn shown(k: &Kept, with_x: bool) -> Result<Shown, String> {
    let mut set = k.set.clone();
    let n = set.members.len();
    let mut info = vec![ShownMember::default(); n];
    let mut without = Vec::new();
    if let Some(c) = &k.rule.norm {
        let rule = normalisation(c, &k.set, k.rule.norm_ref.as_deref())?;
        for x in plots::normalise(&k.set, &rule)? {
            match x.pattern {
                Some(p) => {
                    set.members[x.member].pattern = p;
                    info[x.member].factor = x.factor;
                }
                None => info[x.member].left_out = Some(format!("not normalised: {}", x.why_not.unwrap_or_default())),
            }
        }
    }
    if let (true, Some(r)) = (with_x, &k.rule.x) {
        let given: Vec<(usize, f64, LambdaFrom)> = r.given.iter().filter_map(|(key, l, f)| set.members.iter().position(|m| key_of(m) == *key).map(|i| (i, *l, f.clone()))).collect();
        let v = plots::x_view(&set, r.quantity, &given, r.fallback.clone());
        without = v.without.clone();
        for x in v.members {
            let i = x.member;
            if info[i].left_out.is_some() {
                continue;
            }
            match x.pattern {
                Some(p) => {
                    set.members[i].pattern = p;
                    info[i].lambda = x.lambda;
                    info[i].lambda_from = x.from.as_ref().map(lambda_words);
                }
                None => info[i].left_out = Some(format!("not in {}: {}", r.quantity.label(), x.why_not.unwrap_or_default())),
            }
        }
    }
    let index: Vec<usize> = (0..n).filter(|&i| info[i].left_out.is_none()).collect();
    if index.len() < n {
        set.members = index.iter().map(|&i| set.members[i].clone()).collect();
    }
    Ok(Shown { set, index, members: info, without })
}

/// The set normalised as the view shows it, with the lines a file written from it carries to say so (PL25,
/// `PlotSet::as_shown`): the set as read and no line when the view normalises nothing; refused, naming the member, when a
/// member cannot be normalised by the rule.
fn as_shown(k: &Kept) -> Result<(PlotSet, Vec<String>), String> {
    match &k.rule.norm {
        Some(c) => k.set.as_shown(&normalisation(c, &k.set, k.rule.norm_ref.as_deref())?),
        None => Ok((k.set.clone(), Vec::new())),
    }
}

/// Normalise… and Convert X… (PL6, PL10): the view's rule set, checked by making the view once; the rule as it was
/// when the core refuses the new one (a reference member without a positive area, a control file without λ).
pub fn set_view_rule(sets: &Sets, id: u32, norm: Option<NormChoice>, x: Option<XChoice>) -> Result<SetView, String> {
    sets.with(id, |k| {
        let norm = norm.filter(|c| c.by != "none");
        let norm_ref = match norm.as_ref().filter(|c| c.by == "reference") {
            Some(c) => Some(c.member.and_then(|i| k.set.members.get(i)).map(key_of).ok_or("a normalisation by a reference member needs the member: choose one")?),
            None => None,
        };
        let x = x.map(|c| x_rule(&c, &k.set)).transpose()?;
        let old = std::mem::replace(&mut k.rule, ViewRule { norm, norm_ref, x });
        if let Err(e) = shown(k, true) {
            k.rule = old;
            return Err(e);
        }
        Ok(view(id, k))
    })
}


// ------------------------------------------------------------------------------------------------- opening a set

/// Files given by name and content opened as a new set (the app's `plots_open_files`, PL8): read by the core
/// (`plots::open_from`), every scan a member; a CSV of several columns gives one member per y column.
/// With `rebin`, every member rebinned to that step in its x on opening (W9, `PlotSet::rebin`: the counts averaged, σ
/// propagated, as `drs data op rebin step=… mode=average`), the step kept for the set file.
/// With `x_interval`, every member cut to it first (PL42, `PlotSet::cut`, in the members' own x), kept by the set.
pub fn open_from(sets: &Sets, files: Vec<OpenedFile>, reader: Option<&str>, title: Option<String>, rebin: Option<f64>, x_interval: Option<[f64; 2]>) -> Result<SetView, String> {
    if files.is_empty() {
        return Err("no file chosen".into());
    }
    let (base, files) = relative_to_first(files);
    let mut set = plots::open_from(&files, reader, title)?;
    if let Some(iv) = x_interval {
        set.cut(iv, None)?;
    }
    if let Some(step) = rebin {
        set.rebin(step)?;
    }
    let id = sets.insert(set, base);
    sets.with(id, |k| {
        k.rebin = rebin;
        Ok(view(id, k))
    })
}

/// The names relative to the first file's folder, as the app's `open_files` makes them of real paths (the host's
/// virtual paths `/local/<folder>/…`); a name in another folder kept whole. The folder, and the files.
fn relative_to_first(files: Vec<OpenedFile>) -> (PathBuf, Vec<OpenedFile>) {
    let base = files.first().and_then(|f| Path::new(&f.name).parent()).map(Path::to_path_buf).unwrap_or_default();
    let files = files.into_iter().map(|f| OpenedFile { name: Path::new(&f.name).strip_prefix(&base).map_or(f.name.clone(), |r| r.to_string_lossy().into_owned()), ..f }).collect();
    (base, files)
}

/// A file among those given, by the path the core joins (`base` and a name) or by the name itself.
fn given<'a>(files: &'a [OpenedFile], base: &Path, name: &str) -> Option<&'a OpenedFile> {
    let full = base.join(name).to_string_lossy().into_owned();
    files.iter().find(|f| f.name == full).or_else(|| files.iter().find(|f| f.name == name))
}

/// The files a buffer or a set file names (W5), joined to the folder the core reads them from as the app joins real
/// paths: a set file's members to the folder above `Plots/` (`plots::set_base`), a buffer's names to its own folder.
/// A set file is known by its content (`drs plot set 1`, PL24), else by its extension.
pub fn names_in(path: &str, text: &str) -> Result<Vec<String>, String> {
    let p = Path::new(path);
    if text.trim_start().starts_with("format") && text.contains(plots::SET_FORMAT) || path.to_lowercase().ends_with(".drv") {
        let f = plots::SetFile::read(text).map_err(|d| format!("{path}: {}", d.iter().map(|x| x.message.clone()).collect::<Vec<_>>().join("; ")))?;
        let base = plots::set_base(p);
        let mut out: Vec<String> = Vec::new();
        for m in &f.members {
            let full = base.join(&m.file).to_string_lossy().into_owned();
            if !out.contains(&full) {
                out.push(full);
            }
        }
        Ok(out)
    } else {
        let base = p.parent().map(Path::to_path_buf).unwrap_or_default();
        Ok(plots::read_buf(text, &|_| true).entries.iter().map(|(n, _)| base.join(n).to_string_lossy().into_owned()).collect())
    }
}

/// An empty set (New Plot Set…): files are added to it with `add_from`.
pub fn new_set(sets: &Sets) -> Result<SetView, String> {
    let set = plots::open_from(&[], None, None)?;
    let id = sets.insert(set, PathBuf::new());
    sets.with(id, |k| Ok(view(id, k)))
}

/// A WinPLOTR buffer given by its text, its files by name and content (the app's `plots_open_buffer`): its names, its
/// comments and its title read by the core (`plots::read_buf`), the files it names taken from `files`.
pub fn open_buffer_from(sets: &Sets, path: Option<&str>, text: &str, files: Vec<OpenedFile>, x_interval: Option<[f64; 2]>) -> Result<SetView, String> {
    // the buffer's folder, as the app's `open_buffer` takes it; the names relative to it
    let base = path.and_then(|p| Path::new(p).parent()).map(Path::to_path_buf).unwrap_or_default();
    let b = plots::read_buf(text, &|name| given(&files, &base, name).is_some());
    if b.entries.is_empty() {
        return Err(format!("{}: no file named", path.unwrap_or("the buffer")));
    }
    let opened: Vec<OpenedFile> = b
        .entries
        .iter()
        .map(|(name, comment)| given(&files, &base, name).map(|f| OpenedFile { name: name.clone(), comment: comment.clone(), ..f.clone() }).ok_or(format!("{name}: the buffer names it, and no file of that name was given")))
        .collect::<Result<_, _>>()?;
    let mut set = plots::open_from(&opened, None, b.title)?;
    // the x interval of the analysis (PL42)
    if let Some(iv) = x_interval {
        set.cut(iv, None)?;
    }
    let id = sets.insert(set, base);
    sets.with(id, |k| Ok(view(id, k)))
}

/// The set's y extent on a scale, as the view shows it (PL41, W4's shared y interval; the app's `plots_extent`): the
/// core's `plots::extent` over the members named (by their index in the set), else all shown.
pub fn extent(sets: &Sets, id: u32, scale: &str, members: Option<Vec<usize>>) -> Result<Option<(f64, f64)>, String> {
    let scale = scale_of(scale)?;
    sets.with(id, |k| {
        let s = shown(k, true)?;
        let at: Option<Vec<usize>> = members.map(|ms| (0..s.index.len()).filter(|&p| ms.contains(&s.index[p])).collect());
        Ok(plots::extent(&s.set, scale, at.as_deref()))
    })
}

/// The series the files offer (W2, `plots::propose`): each a set file the host opens as it is, its members' files and
/// the files left out named relative to the first file's folder (as `open_from` names members), the folder whose
/// `Plots/` the set file is saved in.
pub fn propose(files: Vec<OpenedFile>, x_interval: Option<[f64; 2]>) -> plots::Proposals {
    plots::propose(&relative_to_first(files).1, x_interval)
}

/// A member's other variable typed with its source (W3; Plots' step 21, the app's `set_other`): kept with the set
/// (`Member::others`, `Found::Typed`), written into the set file as `variables = { name = { value, source } }`, and
/// given back on the map's rows for the right-hand axis. A value already there under the name is replaced.
pub fn set_other(sets: &Sets, id: u32, member: usize, name: &str, value: f64, note: &str) -> Result<SetView, String> {
    let name = name.trim();
    if name.is_empty() {
        return Err("a variable needs a name".into());
    }
    if !value.is_finite() {
        return Err(format!("{name}: give a number"));
    }
    if note.trim().is_empty() {
        return Err(format!("{name}: say where the value comes from (L169)"));
    }
    sets.with(id, |k| {
        let m = k.set.members.get_mut(member).ok_or(format!("no member {}", member + 1))?;
        m.others.retain(|o| o.name != name);
        m.others.push(plots::Other { name: name.to_string(), value, found: plots::Found::Typed(note.trim().to_string()), unit: String::new() });
        Ok(view(id, k))
    })
}

/// Files added to a set (Add Files…, the app's `plots_add`), by name and content: their members joined to the set's,
/// the whole sorted by the variable again.
pub fn add_from(sets: &Sets, id: u32, files: Vec<OpenedFile>) -> Result<SetView, String> {
    let mut more = plots::open_from(&files, None, None)?;
    sets.with(id, |k| {
        // a set cut to the analysis's interval (PL42) takes the files added cut alike
        if let Some(iv) = k.set.x_interval {
            more.cut(iv, k.set.x_interval_in)?;
        }
        k.set.members.extend(more.members);
        k.set.sort_by_variable();
        Ok(view(id, k))
    })
}

pub fn get(sets: &Sets, id: u32) -> Result<SetView, String> {
    sets.with(id, |k| Ok(view(id, k)))
}

// --------------------------------------------------------------------------------------------- the set file (PL12)

/// What the window keeps of a set that the set file holds besides its members (PL12, PL27): the order, the stack's
/// rule and the scale it was set on, the projection's steps and their scale, the y scale, the map's colour scale, the
/// reference member, the windows, each member's opacity. (The view's normalisation and x quantity are the Rust side's
/// own rule, written and read with them.)
#[derive(serde::Deserialize, Serialize, Clone, Debug, Default)]
pub struct Settings {
    /// "variable", "name" or "given".
    pub order: String,
    /// The stack's displacement: "none", "step" or "variable".
    pub rule: String,
    pub dx: f64,
    pub dy: f64,
    pub scale: String,
    /// The scale the stack's step was set on (PL27); none: `scale`.
    #[serde(default)]
    pub stack_scale: Option<String>,
    /// The projection's steps (PL27): "none", "step" or "variable", Δx, Δy, and the scale they were set on.
    #[serde(default)]
    pub prule: String,
    #[serde(default)]
    pub pdx: f64,
    #[serde(default)]
    pub pdy: f64,
    #[serde(default)]
    pub pscale: Option<String>,
    /// The map's colour scale (PL27).
    #[serde(default)]
    pub map_scale: Option<String>,
    /// The reference member's index in the set.
    pub reference: Option<usize>,
    pub windows: Vec<plots::Window>,
    /// Each member's opacity, by its index in the set (1 when not given).
    pub opacity: Vec<Option<f64>>,
    /// The map's and the difference map's colour range as clipped (PL37, I08; Plots' step 17), on the map's scale;
    /// none: the core's range.
    #[serde(default)]
    pub map_range: Option<plots::Range>,
    #[serde(default)]
    pub difference_range: Option<plots::Range>,
    /// The regions of the variable–x plane saved by name (PL37), in the window's words.
    #[serde(default)]
    pub regions: Vec<RegionView>,
    /// What the core refused of the file's regions when the set was opened (`SetFile::check_regions`: a region in
    /// another quantity whose conversion needs a λ the members lack), in its words; the regions are kept as read.
    #[serde(default)]
    pub regions_refused: Option<String>,
    /// The y interval (PL41, core 0.80.0): "shared" or "own"; none written reads as shared.
    #[serde(default)]
    pub y_interval: Option<String>,
}

/// A region as the window names it (PL37, `plots::Region`): `x_quantity` as a word ("2theta", "d", "q", "stl", "s"), or
/// none for the set's x.
#[derive(serde::Deserialize, Serialize, Clone, Debug, PartialEq)]
pub struct RegionView {
    pub name: String,
    pub x: [f64; 2],
    #[serde(default)]
    pub x_quantity: Option<String>,
    pub variable: [f64; 2],
    /// The region's own colour range on the map (P65, Plots' step 19), and on the difference map, apart: the window
    /// writes the one of the view the region was saved on and applies it to that view only.
    #[serde(default)]
    pub map_range: Option<plots::Range>,
    #[serde(default)]
    pub difference_range: Option<plots::Range>,
    #[serde(default)]
    pub window: Option<String>,
}

impl RegionView {
    /// The core's region, with both of its own ranges (P65, core 0.78.0): the map's in `map_range`, the difference
    /// map's in `difference_range`, each applied to its view only.
    fn to_core(&self) -> Result<plots::Region, String> {
        let x_quantity = self.x_quantity.as_deref().map(|w| XQuantity::from_word(w).ok_or(format!("region \"{}\": no x quantity \"{w}\" (2theta, d, q, stl, s)", self.name))).transpose()?;
        Ok(plots::Region { name: self.name.clone(), x: self.x, x_quantity, variable: self.variable, map_range: self.map_range, difference_range: self.difference_range, window: self.window.clone() })
    }
    fn of_core(r: &plots::Region) -> RegionView {
        RegionView { name: r.name.clone(), x: r.x, x_quantity: r.x_quantity.map(|q| q.word().to_string()), variable: r.variable, map_range: r.map_range, difference_range: r.difference_range, window: r.window.clone() }
    }
}

/// A file the app would write to the project, for the browser's host to save: its name relative to the project's folder
/// and its text.
#[derive(Serialize, Clone, Debug)]
pub struct SavedText {
    pub relative: String,
    /// The same name, as W5's `plots_save_text` asks it.
    pub path: String,
    pub text: String,
}

/// The set as the project's set file, `Plots/{name}.drv` (the file dictionary's `plot_set`), through the core
/// (`SetFile::of`, `to_text`): each member's file relative to the project's folder (`plots::relative_to`), its typed
/// value with its source, its opacity; the set's order, rule, y scale, reference and windows. The browser has no project
/// to write to: the file's name in the project and its text are returned, for the host to save (the app writes it).
// The app's code kept as it is, `!(from < to)` refusing a NaN as well as a pair out of order.
#[allow(clippy::neg_cmp_op_on_partial_ord)]
pub fn save(sets: &Sets, id: u32, root: &str, name: &str, st: Settings) -> Result<SavedText, String> {
    let rel = drs::formats::files::path("plot_set", &[("set", name)])?;
    sets.with(id, |k| {
        let mut f = plots::SetFile::of(&k.set);
        f.rebin = k.rebin;
        for (i, m) in f.members.iter_mut().enumerate() {
            m.file = plots::relative_to(&k.base.join(&m.file), Path::new(root));
            m.opacity = st.opacity.get(i).copied().flatten().filter(|&o| o != 1.0);
        }
        f.order = match st.order.as_str() {
            "name" => plots::Order::Name,
            "given" => plots::Order::Given,
            _ => plots::Order::Variable,
        };
        let rule_of = |rule: &str, dx: f64, dy: f64| match rule {
            "step" => Displacement::Step { dx, dy },
            "variable" => Displacement::Variable { dx, dy },
            _ => Displacement::None,
        };
        f.displacement = rule_of(&st.rule, st.dx, st.dy);
        f.displacement_scale = st.stack_scale.as_deref().and_then(YScale::from_word);
        f.projection = rule_of(&st.prule, st.pdx, st.pdy);
        f.projection_scale = st.pscale.as_deref().and_then(YScale::from_word);
        f.map_scale = st.map_scale.as_deref().and_then(YScale::from_word);
        f.y_scale = YScale::from_word(&st.scale).unwrap_or(YScale::Linear);
        f.y_interval = if st.y_interval.as_deref() == Some("own") { plots::YInterval::Own } else { plots::YInterval::Shared };
        // the view's normalisation and x quantity (PL27), as the view shows the set
        if let Some(c) = &k.rule.norm {
            f.normalisation = normalisation(c, &k.set, k.rule.norm_ref.as_deref())?;
        }
        if let Some(x) = &k.rule.x {
            f.x = x.quantity.word().to_string();
        }
        f.reference = st.reference.and_then(|r| f.members.get(r)).map(|m| m.file.clone());
        f.windows = st.windows.clone();
        // the map's clip and the regions (PL37, step 17), as the window keeps them; the core refuses a pair not in order
        for (what, r) in [("map_range", st.map_range), ("difference_range", st.difference_range)] {
            if let Some(r) = r {
                if !(r.from < r.to) {
                    return Err(format!("{what}: from = {} is not below to = {}", r.from, r.to));
                }
            }
        }
        f.map_range = st.map_range;
        f.difference_range = st.difference_range;
        f.regions = st.regions.iter().map(RegionView::to_core).collect::<Result<Vec<_>, String>>()?;
        Ok(SavedText { relative: rel.clone(), path: rel.clone(), text: f.to_text() })
    })
}

/// A set file given by its name and text, its members' files by name and content (the app's `plots_open_set`), through
/// the core (`SetFile::read`, `open_from`): the typed values with their sources, the order; and the settings the file
/// holds, for the window.
/// With `x_interval`, the set cut to it in the file's x quantity after its own (PL42): the file's interval replaced.
pub fn open_set_from(sets: &Sets, name: &str, text: &str, files: Vec<OpenedFile>, x_interval: Option<[f64; 2]>) -> Result<(SetView, Settings), String> {
    let f = plots::SetFile::read(text).map_err(|d| format!("{name}: {}", d.iter().map(|x| x.message.clone()).collect::<Vec<_>>().join("; ")))?;
    // the folder above Plots/ (`plots::set_base`), as the app's `open_set_file`; each member's file found among those
    // given by the path joined to it, or by its own name; one not given is named by the core as missing
    let base = plots::set_base(Path::new(name));
    let named: Vec<OpenedFile> = f.members.iter().filter_map(|m| given(&files, &base, &m.file).map(|g| OpenedFile { name: m.file.clone(), ..g.clone() })).collect();
    let mut set = f.open_from(&named, None, &|_| Vec::new())?;
    if let Some(iv) = x_interval {
        set.cut(iv, f.x_quantity())?;
    }
    let reference = f.reference_in(&set);
    let opacity = set.members.iter().map(|m| f.members.iter().find(|x| x.file == m.file && x.scan.is_none_or(|s| s == m.scan + 1)).and_then(|x| if x.visible == Some(false) { Some(0.0) } else { x.opacity })).collect();
    let words = |d: Displacement| match d {
        Displacement::None => ("none", 0.0, 0.0),
        Displacement::Step { dx, dy } => ("step", dx, dy),
        Displacement::Variable { dx, dy } => ("variable", dx, dy),
    };
    let (rule, dx, dy) = words(f.displacement);
    let (prule, pdx, pdy) = words(f.projection);
    let order = match f.order {
        plots::Order::Variable => "variable",
        plots::Order::Name => "name",
        plots::Order::Given => "given",
    };
    let st = Settings {
        order: order.into(),
        rule: rule.into(),
        dx,
        dy,
        scale: f.y_scale.word().into(),
        stack_scale: f.displacement_scale.map(|s| s.word().to_string()),
        prule: prule.into(),
        pdx,
        pdy,
        pscale: f.projection_scale.map(|s| s.word().to_string()),
        map_scale: f.map_scale.map(|s| s.word().to_string()),
        reference,
        windows: f.windows.clone(),
        opacity,
        map_range: f.map_range,
        difference_range: f.difference_range,
        regions: f.regions.iter().map(RegionView::of_core).collect(),
        // a region the core refuses on this set (its conversion needs a λ the members lack) is said, the set still opened
        regions_refused: f.check_regions(&set).err(),
        y_interval: Some(f.y_interval.word().to_string()),
    };
    // the view's normalisation and x quantity the file keeps (PL27), as the window's rule; the members' typed λs are
    // not kept by the file
    let (norm, norm_ref) = choice_of(&f.normalisation, &set);
    let x = XQuantity::from_word(&f.x).filter(|q| *q != XQuantity::TwoTheta || set.members.iter().any(|m| m.facts.axis != Axis::TwoTheta)).map(|quantity| XRule { quantity, given: Vec::new(), fallback: None });
    let id = sets.insert(set, base);
    sets.with(id, |k| {
        k.rule = ViewRule { norm, norm_ref, x };
        k.rebin = f.rebin;
        Ok((view(id, k), st))
    })
}

// ------------------------------------------------------------------------------------------ the stack's step (PL19)

/// The core's proposed step for the stack on the scale shown (`plots::proposed_step`), with its rule in words.
pub fn proposed(sets: &Sets, id: u32, scale: &str) -> Result<Option<plots::ProposedStep>, String> {
    let scale = scale_of(scale)?;
    sets.with(id, |k| Ok(plots::proposed_step(&shown(k, true)?.set, scale)))
}

// -------------------------------------------------------------------------------- FullProf's temperatures file (PL21)

/// What reading a temperatures file did: the set, how many members took a value, the names given twice with
/// different values (which attribute nothing), the lines not read.
#[derive(Serialize, Clone, Debug)]
pub struct TemperaturesRead {
    pub set: SetView,
    pub taken: usize,
    pub conflicts: Vec<String>,
    pub unread: usize,
}

/// FullProf's temperatures file read by the core (`plots::read_temperatures`) and each member whose file it names
/// given that value with the file as its source (`PlotSet::variable_from_file`); typed values kept.
pub fn temperatures_from(sets: &Sets, id: u32, name: &str, text: &str) -> Result<TemperaturesRead, String> {
    let t = plots::read_temperatures(text);
    let name = Path::new(name).file_name().map_or(name.to_string(), |f| f.to_string_lossy().into_owned());
    sets.with(id, |k| {
        let taken = k.set.variable_from_file(&name, &t, false);
        Ok(TemperaturesRead { set: view(id, k), taken, conflicts: t.conflicts.iter().map(|c| c.0.clone()).collect(), unread: t.unread.len() })
    })
}

/// The user's value of the variable for a member, with their note of its source (L169: what a file does not say is
/// typed, and keeps where it came from); the core sorts the members again.
pub fn set_variable(sets: &Sets, id: u32, member: usize, value: f64, note: String) -> Result<SetView, String> {
    sets.with(id, |k| {
        k.set.set_variable(member, value, &note)?;
        Ok(view(id, k))
    })
}

// ------------------------------------------------------------------------------------------------- displacements

fn scale_of(s: &str) -> Result<YScale, String> {
    match s {
        "linear" => Ok(YScale::Linear),
        "sqrt" => Ok(YScale::Sqrt),
        "log" => Ok(YScale::Log),
        _ => Err(format!("no y scale \"{s}\" (linear, sqrt, log)")),
    }
}

/// One member displaced (`plots::Displaced`), with its index in the set; None for a member the rule cannot place (the
/// rule by the variable, for a member without a value).
#[derive(Serialize, Clone, Debug)]
pub struct DisplacedView {
    pub member: usize,
    pub placed: bool,
    pub x: Vec<f64>,
    /// y on the scale shown, the shift added.
    pub y: Vec<Option<f64>>,
    pub shift: (f64, f64),
    pub baseline: f64,
}

/// Each member displaced by the set's rule (`plots::displace`): `rule` "none", "step" (k·dx, k·dy) or "variable"
/// ((v − v₀)·dx, (v − v₀)·dy), on the y scale shown.
pub fn displace(sets: &Sets, id: u32, rule: &str, dx: f64, dy: f64, scale: &str) -> Result<Vec<DisplacedView>, String> {
    let rule = match rule {
        "none" => Displacement::None,
        "step" => Displacement::Step { dx, dy },
        "variable" => Displacement::Variable { dx, dy },
        _ => return Err(format!("no displacement rule \"{rule}\" (none, step, variable)")),
    };
    let scale = scale_of(scale)?;
    sets.with(id, |k| {
        let s = shown(k, true)?;
        let d = plots::displace(&s.set, rule, scale, &[], None);
        Ok(s.spread(d)
            .into_iter()
            .enumerate()
            .map(|(i, d)| match d.flatten() {
                Some(d) => DisplacedView { member: i, placed: true, x: d.x, y: d.y.into_iter().map(|v| v.is_finite().then_some(v)).collect(), shift: d.shift, baseline: d.baseline },
                None => DisplacedView { member: i, placed: false, x: Vec::new(), y: Vec::new(), shift: (0.0, 0.0), baseline: 0.0 },
            })
            .collect())
    })
}

// ------------------------------------------------------------------------------------------------------------ the map

/// A member without a value of the variable, as the map draws it apart (L169): its values on the map's grid and scale.
#[derive(Serialize, Clone, Debug)]
pub struct ApartRow {
    pub member: usize,
    /// None where the member has no point or the scale has no value; the whole row None when the core does not give
    /// it (a difference map: request PL18).
    pub z: Option<Vec<Option<f64>>>,
}

#[derive(Serialize, Clone, Debug)]
pub struct MapView {
    pub map: plots::Map,
    pub apart: Vec<ApartRow>,
    /// How each member was put on the grid (`plots::Resampling`), in words.
    pub resampling: Vec<String>,
}

/// The map of a set (`plots::map`) on the common grid (`plots::common_grid`, the overlap at the finest step: the
/// members' own points when they share them); with a reference member, the difference map. The members without a
/// value are drawn apart: for the plain map, their rows on the grid and scale (`plots::on_grid`, `YScale::apply`); for
/// the difference map the core does not give them yet (PL18).
pub fn map(sets: &Sets, id: u32, scale: &str, reference: Option<usize>) -> Result<MapView, String> {
    let scale = scale_of(scale)?;
    sets.with(id, |k| {
        let (m, s, res) = map_of(k, scale, reference)?;
        let apart = m.apart.iter().zip(&m.apart_z).map(|(&a, z)| ApartRow { member: a, z: Some(z.clone()) }).collect();
        let resampling = s.spread(res).iter().enumerate().map(|(i, r)| r.as_ref().map_or_else(|| s.members[i].left_out.clone().unwrap_or_default(), |r| format!("{:?}", r.how))).collect();
        Ok(MapView { map: m, apart, resampling })
    })
}

/// The map as the window draws it (the set as the view shows it, the overlap at the finest step, the scale, a reference
/// for the difference map), its rows' members by their index in the set; with the view and each member's resampling.
fn map_of(k: &Kept, scale: YScale, reference: Option<usize>) -> Result<(plots::Map, Shown, Vec<plots::Resampled>), String> {
    if k.set.members.is_empty() {
        return Err("the set has no member: Add Files… or Open Buffer…".into());
    }
    let s = shown(k, true)?;
    if s.set.members.is_empty() {
        return Err("the view leaves every member out: see why in the members' list".into());
    }
    let reference = reference.map(|r| s.pos(r, "the reference member")).transpose()?;
    let grid = plots::common_grid(&s.set, GridRange::Overlap, GridStep::Finest)?;
    let mut m = plots::map(&s.set, &grid, scale, reference)?;
    let res = plots::on_grid(&s.set, &grid);
    // the members by their index in the set; the rows of the members apart are the core's, less the reference in a
    // difference map (PL18)
    for r in m.rows.iter_mut() {
        r.member = s.index[r.member];
    }
    m.apart = m.apart.iter().map(|&a| s.index[a]).collect();
    m.reference = m.reference.map(|r| s.index[r]);
    Ok((m, s, res))
}

/// The x the window's map is drawn in (PL37, `plots::MapX`): the view's quantity with each member's λ as the view
/// converted it (`MapX::of_view`, the view made again on the set as read, which decides each member's λ as `shown`
/// does), else the files' (`MapX::of`); indexed by the members' indices in the set, as the map's rows name them.
fn map_x(k: &Kept) -> plots::MapX {
    match &k.rule.x {
        Some(r) => {
            let given: Vec<(usize, f64, LambdaFrom)> = r.given.iter().filter_map(|(key, l, f)| k.set.members.iter().position(|m| key_of(m) == *key).map(|i| (i, *l, f.clone()))).collect();
            let v = plots::x_view(&k.set, r.quantity, &given, r.fallback.clone());
            plots::MapX::of_view(&k.set, &v)
        }
        None => plots::MapX::of(&k.set),
    }
}

/// A region's block of the map as the window draws it (PL37, I08's zoom and "range from the view"; step 17):
/// `plots::map_region`, the rows and columns inside and the range of its cells, a reading of the matrix the core
/// computed; refused in the core's words (a λ the conversion needs and members lack, no row, no column, no value).
pub fn map_region(sets: &Sets, id: u32, scale: &str, reference: Option<usize>, region: RegionView) -> Result<plots::MapBlock, String> {
    let scale = scale_of(scale)?;
    let r = region.to_core()?;
    sets.with(id, |k| {
        let (m, _, _) = map_of(k, scale, reference)?;
        plots::map_region(&m, &r, &map_x(k))
    })
}

/// The core's word on a colour range against the map's values (PL37, `plots::range_warning`): a warning when it lies
/// wholly outside them, every cell then drawn saturated; none otherwise.
pub fn range_warning(sets: &Sets, id: u32, scale: &str, reference: Option<usize>, range: plots::Range) -> Result<Option<String>, String> {
    let scale = scale_of(scale)?;
    sets.with(id, |k| {
        let (m, _, _) = map_of(k, scale, reference)?;
        Ok(plots::range_warning(if reference.is_some() { "the difference map's colour range" } else { "the map's colour range" }, &range, &m))
    })
}

// ------------------------------------------------------------------------------------------------- derived patterns

/// A derived pattern written: where, and what the core says of it.
#[derive(Serialize, Clone, Debug)]
pub struct DerivedView {
    pub path: String,
    pub relative: String,
    pub points: usize,
    pub dropped: usize,
    pub derivation: String,
    /// The derived pattern's file, for the host to save (the browser has no project).
    pub text: String,
}

/// The sum, average or difference of members (`plots::derive` on the common grid), written as the file dictionary's
/// `derived_data` (`Plots/{set} {operation}.xys`, `plots::derived_text_with`) in the project. Nothing is overwritten.
/// The members as the view shows them (PL25): normalised by the view's rule, the file's `!` lines saying so and each
/// member's factor (`PlotSet::as_shown`); in the files' x.
pub fn derive(sets: &Sets, id: u32, members: Vec<usize>, kind: &str, root: &str, name: &str) -> Result<DerivedView, String> {
    let _ = root;
    let derivation = match kind {
        "sum" => Derivation::Sum,
        "average" => Derivation::Average,
        "weighted_average" => Derivation::WeightedAverage,
        "difference" => Derivation::Difference,
        _ => return Err(format!("no derivation \"{kind}\" (sum, average, weighted_average, difference)")),
    };
    sets.with(id, |k| {
        let (set, notes) = as_shown(k)?;
        let grid = plots::common_grid(&set, GridRange::Overlap, GridStep::Finest)?;
        let d = plots::derive(&set, &members, derivation, &grid)?;
        let rel = drs::formats::files::path("derived_data", &[("set", name), ("operation", derivation.word())])?;
        // the browser has no project to write to: the file's text is returned for the host to save (the app writes it)
        Ok(DerivedView { path: rel.clone(), relative: rel, points: d.pattern.len(), dropped: d.dropped, derivation: derivation.word().to_string(), text: plots::derived_text_with(&set, &d, name, &notes) })
    })
}

// ------------------------------------------------------------------------------------------- the projection (PL17)

/// The projection with its members apart (L183): the core's projection, and the members without a value of the
/// variable, which stand apart at their measured 2θ.
#[derive(Serialize, Clone, Debug)]
pub struct ProjectionView {
    pub projection: plots::Projection,
    pub apart: Vec<usize>,
}

/// The projection with hidden parts (`plots::displace`, then `plots::project`): the members placed by the rule, front
/// to back in `order` (the set's order of placed members when none is given), each curve hiding what lies under it down
/// to `floor` plus its shift (without end when none). The window draws the visible pieces as they are.
///
/// A member without a value of the variable has no place behind the others (L169, L183): with a step rule it keeps
/// the core's own step in y (its place at the back, as in the stack) but is not moved in 2θ, its peaks where they were
/// measured; the core applies that as the member's override (`plots::displace`'s `overrides`). The rule by the
/// variable places no such member at all.
#[allow(clippy::too_many_arguments)] // the window's command's own fields, passed through
pub fn project(sets: &Sets, id: u32, rule: &str, dx: f64, dy: f64, scale: &str, floor: Option<f64>, order: Option<Vec<usize>>) -> Result<ProjectionView, String> {
    let rule = match rule {
        "step" => Displacement::Step { dx, dy },
        "variable" => Displacement::Variable { dx, dy },
        _ => return Err(format!("no projection rule \"{rule}\" (step, variable)")),
    };
    let scale = scale_of(scale)?;
    sets.with(id, |k| {
        let s = shown(k, true)?;
        let apart = s.set.without_variable();
        let first = plots::displace(&s.set, rule, scale, &[], None);
        let overrides: Vec<(usize, (f64, f64))> = apart.iter().filter_map(|&m| first[m].as_ref().map(|d| (m, (0.0, d.shift.1)))).collect();
        let d = plots::displace(&s.set, rule, scale, &overrides, None);
        let order = match order {
            Some(o) => o.iter().map(|&m| s.pos(m, "a member of the order")).collect::<Result<Vec<_>, String>>()?,
            None => (0..d.len()).filter(|&i| d[i].is_some()).collect(),
        };
        // the members by their index in the set
        let mut p = plots::project(&d, &order, floor)?;
        p.order = p.order.iter().map(|&i| s.index[i]).collect();
        for c in p.curves.iter_mut() {
            c.member = s.index[c.member];
        }
        for a in p.axis.iter_mut() {
            a.0 = s.index[a.0];
        }
        Ok(ProjectionView { projection: p, apart: apart.iter().map(|&i| s.index[i]).collect() })
    })
}

// ------------------------------------------------------------------------------------- F6: the overlay's difference

/// One member less another on the common grid with σ (`plots::derive`, a difference, not written to a file): the strip
/// under the linked overlay (F6), of the members as the view shows them.
#[derive(Serialize, Clone, Debug)]
pub struct DifferenceView {
    pub a: usize,
    pub b: usize,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub sigma: Vec<f64>,
    pub dropped: usize,
}

pub fn difference(sets: &Sets, id: u32, a: usize, b: usize) -> Result<DifferenceView, String> {
    sets.with(id, |k| {
        let s = shown(k, true)?;
        let (pa, pb) = (s.pos(a, "the first member")?, s.pos(b, "the member subtracted")?);
        let grid = plots::common_grid(&s.set, GridRange::Overlap, GridStep::Finest)?;
        let d = plots::derive(&s.set, &[pa, pb], Derivation::Difference, &grid)?;
        Ok(DifferenceView { a, b, x: d.pattern.x, y: d.pattern.y, sigma: d.pattern.sigma, dropped: d.dropped })
    })
}

// ------------------------------------------------------------------------------------ the order and the variable

/// The members in the reverse of their order (Reverse Order): a reordering, nothing computed.
pub fn reverse(sets: &Sets, id: u32) -> Result<SetView, String> {
    sets.with(id, |k| {
        k.set.members.reverse();
        Ok(view(id, k))
    })
}

/// The members sorted by the variable (`PlotSet::sort_by_variable`: those with a value by it, then the others).
pub fn sort(sets: &Sets, id: u32) -> Result<SetView, String> {
    sets.with(id, |k| {
        k.set.sort_by_variable();
        Ok(view(id, k))
    })
}

/// Read Variable From… (L183): each member's value taken from one of the sources the core found for it
/// (`plots::temperature_candidates`): "name", the file's name (`Found::FileName`); "first", the first source in the
/// core's order (the header, the buffer's comment, the title, the name), as the set was opened. A member whose files say
/// nothing there keeps its value; a typed value is kept unless `typed` is false. The core sorts the members again.
/// Returns the set and how many members took a value from the source.
pub fn read_variable(sets: &Sets, id: u32, source: &str, replace_typed: bool) -> Result<(SetView, usize), String> {
    sets.with(id, |k| {
        let mut n = 0;
        for m in k.set.members.iter_mut() {
            if !replace_typed && matches!(m.variable.as_ref().map(|v| &v.found), Some(plots::Found::Typed(_))) {
                continue;
            }
            let found = match source {
                "name" => m.candidates.iter().find(|c| matches!(c.found, plots::Found::FileName(_))).cloned(),
                "first" => m.candidates.first().cloned(),
                _ => return Err(format!("no source \"{source}\" (name, first)")),
            };
            if let Some(v) = found {
                m.variable = Some(v);
                n += 1;
            }
        }
        k.set.sort_by_variable();
        Ok((view(id, k), n))
    })
}

// ---------------------------------------------------------------------------------------------------- windows (PL5)

/// The windows' numbers for every member (`plots::windows`): the integral above the background with σ, the centroid
/// with σ, the strongest peak; the background at `at`. `background` "auto" (each member's automatic background, at
/// `log_lambda` or GSAS-II's default) or "line" (the line between the member's values at the window's ends). The
/// members normalised as the view shows them (PL6), in the files' x.
pub fn windows(sets: &Sets, id: u32, list: Vec<plots::Window>, background: &str, log_lambda: Option<f64>, at: Vec<f64>) -> Result<plots::Windows, String> {
    let bg = match background {
        "auto" => plots::WindowBackground::Auto { log_lambda },
        "line" => plots::WindowBackground::Line,
        _ => return Err(format!("no window background \"{background}\" (auto, line)")),
    };
    sets.with(id, |k| {
        // the windows are in the files' x: the members normalised as shown, not converted
        let s = shown(k, false)?;
        let mut w = plots::windows(&s.set, &list, &bg, &at)?;
        for m in w.members.iter_mut() {
            m.member = s.index[m.member];
        }
        Ok(w)
    })
}
