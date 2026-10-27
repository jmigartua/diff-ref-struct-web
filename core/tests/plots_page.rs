//! Plots' PL38 (I08 step 5, R33; L409, L412, L413): the page without the app.
//! - `plots::page_data` on the Tb₂Pd₂Sn set with a clip and two regions: the map's rows and range, the difference map
//!   on the scale the module draws it on, the regions' blocks as `map_region` gives them (L407's rows 8–14 and columns
//!   30–32 for 18.95–19.25° at 20–26 K) or their refusal in the core's words, the members' temperatures with their
//!   sources, the scale words; the keys the viewer's header names (`app/ui/viewer/map.js`), no more and no fewer;
//! - a set file without a reference gives no difference map;
//! - `drs plots show --browser --no-open` writes `<set>.map.html` beside the set: the page's JSON and the viewer's
//!   bundle inside it, no `http` other than the two SVG namespaces, nothing left of the shell's markers; with
//!   `--no-patterns` the members' patterns left out;
//! - a set with one member refused in words;
//! - the viewer's assets present and self-consistent (`SOURCES.txt` against the bundle's header), the shell's markers;
//! - `drs plots region` counts from 1 for a person and says the indices from 0, as the blocks do.

use drs::plots::{self, MapX, PageData, Range, Region, SetFile, YScale};
use std::path::{Path, PathBuf};

fn root(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)
}

fn tb() -> PathBuf {
    root("../manuals/fullprof/1-WinPLOTR/1_Tb2pd2sn")
}

fn drs(args: &[&str]) -> std::process::Output {
    std::process::Command::new(env!("CARGO_BIN_EXE_drs")).args(args).output().unwrap()
}

/// The Tb₂Pd₂Sn set (20 files, 19 with T) and its set file: the map on log with a clip, the reference tb_t20, a region
/// around the magnetic peak at 19.1° and one drawn in d (refused: the files state no λ).
fn tb_set(reference: bool) -> (plots::PlotSet, SetFile) {
    let set = plots::open_buf(&tb().join("Tb_vs_T.buf"), None).unwrap();
    let mut f = SetFile::of(&set);
    f.map_scale = Some(YScale::Log);
    f.map_range = Some(Range { from: 3.7, to: 4.46 });
    if reference {
        f.reference = Some(set.members.iter().find(|m| m.label == "tb_t20.dat").unwrap().file.clone());
    }
    f.regions = vec![
        Region { name: "m1 at 20-26 K".into(), x: [18.95, 19.25], x_quantity: None, variable: [20.0, 26.0], map_range: Some(Range { from: 3.9, to: 4.4 }), difference_range: Some(Range { from: -40.0, to: 60.0 }), window: Some("m1".into()) },
        Region { name: "in d".into(), x: [4.25, 4.5], x_quantity: Some(plots::XQuantity::D), variable: [12.0, 26.0], map_range: None, difference_range: None, window: None },
    ];
    (set, f)
}

fn keys(v: &serde_json::Value) -> Vec<String> {
    let mut k: Vec<String> = v.as_object().unwrap().keys().cloned().collect();
    k.sort();
    k
}

#[test]
fn page_data_holds_the_map_and_the_regions() {
    let (set, f) = tb_set(true);
    let x = MapX::of(&set);
    let k20 = f.reference_in(&set).unwrap();
    let d = plots::page_data(&set, &f, &x, Some(k20), "Plots/tb.drv", true).unwrap();
    // The map as plots::map draws it on the file's map scale, on the overlap at the finest step.
    let grid = plots::common_grid(&set, plots::GridRange::Overlap, plots::GridStep::Finest).unwrap();
    let map = plots::map(&set, &grid, YScale::Log, None).unwrap();
    let pm = d.maps.map.as_ref().unwrap();
    assert_eq!((pm.scale.as_str(), pm.reference, pm.map.rows.len(), pm.map.x.len()), ("log", None, 19, 301));
    assert_eq!((&pm.map.rows, pm.map.range, &pm.map.z), (&map.rows, map.range, &map.z));
    assert_eq!(pm.apart.len(), 1);
    assert_eq!((set.members[pm.apart[0].member].label.as_str(), pm.apart[0].z.as_ref().map(|z| z.len())), ("tb_t01.dat", Some(301)));
    // The difference map on the square root, as the module draws it for a log map scale.
    let diff = plots::map(&set, &grid, YScale::Sqrt, Some(k20)).unwrap();
    let pd = d.maps.difference.as_ref().unwrap();
    assert_eq!((pd.scale.as_str(), pd.reference, pd.map.range), ("sqrt", Some(k20), diff.range));
    assert_eq!(plots::difference_scale(YScale::Linear), YScale::Linear);
    // The clip and the regions: m1's blocks as map_region gives them, rows 8–14 and columns 30–32 from 0.
    assert_eq!((d.clip.map, d.clip.difference), (Some(Range { from: 3.7, to: 4.46 }), None));
    let m1 = &d.regions[0];
    let b = m1.blocks.map.as_ref().unwrap();
    assert_eq!(b, &plots::map_region(&map, &f.regions[0], &x).unwrap());
    assert_eq!((b.rows.clone(), b.columns.clone(), b.counted_from), ((8..=14).collect::<Vec<_>>(), vec![30, 31, 32], 0));
    assert_eq!(m1.blocks.difference.as_ref().unwrap(), &plots::map_region(&diff, &f.regions[0], &x).unwrap());
    assert_eq!((m1.refused.as_deref(), m1.map_range, m1.window.as_deref(), m1.x_quantity.as_deref()), (None, Some(Range { from: 3.9, to: 4.4 }), Some("m1"), None));
    // P65: each of the region's ranges, the map's and the difference map's.
    assert_eq!((m1.map_range, m1.difference_range), (Some(Range { from: 3.9, to: 4.4 }), Some(Range { from: -40.0, to: 60.0 })));
    assert_eq!((d.regions[1].map_range, d.regions[1].difference_range), (None, None));
    let in_d = &d.regions[1];
    assert_eq!((in_d.blocks.map.as_ref(), in_d.blocks.difference.as_ref(), in_d.x_quantity.as_deref()), (None, None, Some("d")));
    assert!(in_d.refused.as_deref().unwrap().starts_with("region \"in d\" is in d and the map in 2theta: converting it needs each member's wavelength, and these state none: tb_t06.dat"), "{:?}", in_d.refused);
    assert!(d.regions_refused.as_deref().unwrap().starts_with("region \"in d\" is in d and the set in 2theta"));
    // The members by their index in the set, each temperature with its source in the core's words.
    assert_eq!(d.members.len(), 20);
    for (m, p) in set.members.iter().zip(&d.members) {
        assert_eq!(p.label, m.label);
        assert_eq!(p.variable.as_ref().map(|v| (v.value, v.source.clone())), m.variable.as_ref().map(|v| (v.value, v.found.words())));
        assert_eq!(p.pattern.as_ref().map(|q| (&q.x, &q.y)), Some((&m.pattern.x, &m.pattern.y)));
    }
    assert_eq!(d.members[0].variable.as_ref().map(|v| (v.value, v.source.as_str())), Some((4.0, "the title's \"T.=.4.K\"")));
    assert_eq!(d.members[19].variable, None);
    // The words: x as the app labels it, the core's version, today's date.
    assert_eq!((d.x.label.as_str(), d.x.quantity.as_deref(), d.y_word.as_str(), d.variable.name.as_str(), d.variable.unit.as_str()), ("2θ (°)", Some("2theta"), "counts", "T", "K"));
    assert_eq!((d.path.as_str(), d.core.as_str(), d.same_grid, d.view.as_str()), ("Plots/tb.drv", env!("CARGO_PKG_VERSION"), true, "map"));
    assert_eq!((d.written.len(), &d.written[4..5], &d.written[7..8]), (10, "-", "-"));
    // The keys the viewer's header names, field by field (and the view and the blocks' counted_from beyond it).
    let v = serde_json::to_value(&d).unwrap();
    assert_eq!(keys(&v), ["clip", "core", "maps", "members", "path", "regions", "regions_refused", "same_grid", "title", "variable", "view", "written", "x", "y_word"]);
    assert_eq!(keys(&v["maps"]["map"]), ["apart", "map", "reference", "scale"]);
    assert_eq!(keys(&v["maps"]["map"]["map"]), ["range", "rows", "x", "x_edges", "z"]);
    assert_eq!(keys(&v["maps"]["map"]["map"]["rows"][0]), ["edges", "member", "others", "value"], "W3: each row's other variables beside the contract's keys");
    assert_eq!(keys(&v["maps"]["map"]["apart"][0]), ["member", "z"]);
    assert_eq!(keys(&v["members"][0]), ["label", "pattern", "variable"]);
    assert_eq!(keys(&v["members"][0]["variable"]), ["source", "value"]);
    assert_eq!(keys(&v["regions"][0]), ["blocks", "difference_range", "map_range", "name", "refused", "variable", "window", "x", "x_quantity"]);
    assert_eq!((&v["regions"][0]["map_range"], &v["regions"][0]["difference_range"]), (&serde_json::json!({ "from": 3.9, "to": 4.4 }), &serde_json::json!({ "from": -40.0, "to": 60.0 })));
    assert_eq!(keys(&v["regions"][0]["blocks"]["map"]), ["columns", "counted_from", "range", "rows", "x"]);
    assert_eq!(v["maps"]["map"]["map"]["rows"][0]["edges"], serde_json::json!([3.0, 5.0]));
    assert!(v["maps"]["map"]["map"]["range"].is_array() && v["regions"][0]["blocks"]["map"]["x"].is_array());
    assert_eq!(serde_json::from_value::<PageData>(v).unwrap(), d);
    // Without patterns, the members as they are and no pattern.
    let lean = plots::page_data(&set, &f, &x, Some(k20), "Plots/tb.drv", false).unwrap();
    assert!(lean.members.iter().all(|m| m.pattern.is_none()));
    assert_eq!(lean.maps, d.maps);
}

#[test]
fn no_reference_no_difference_map() {
    let (set, mut f) = tb_set(false);
    f.map_scale = None;
    let d = plots::page_data(&set, &f, &MapX::of(&set), None, "tb.drv", true).unwrap();
    assert_eq!(d.maps.difference, None);
    assert_eq!(d.maps.map.as_ref().unwrap().scale, "linear");
    assert!(d.regions.iter().all(|r| r.blocks.difference.is_none()));
    let v = serde_json::to_value(&d).unwrap();
    assert!(v["maps"]["difference"].is_null() && v["clip"]["difference"].is_null());
}

#[test]
fn one_member_is_refused() {
    let set = plots::open(&tb(), &[("tb_t02.dat".to_string(), None)], None, None).unwrap();
    let f = SetFile::of(&set);
    let e = plots::page_data(&set, &f, &MapX::of(&set), None, "one.drv", true).unwrap_err();
    assert_eq!(e, "a map has one row per value of T, and this set has 1 member with one (tb_t02.dat): no map, so no page");
}

/// A project with the Tb files in Data/ and the set file in Plots/, its clip and regions as `tb_set`'s.
fn project(name: &str) -> PathBuf {
    let d = root("target/test-plots-page").join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(d.join("Data")).unwrap();
    std::fs::create_dir_all(d.join("Plots")).unwrap();
    for k in 1..=20 {
        let f = format!("tb_t{k:02}.dat");
        std::fs::copy(tb().join(&f), d.join("Data").join(&f)).unwrap();
    }
    std::fs::copy(tb().join("Tb_vs_T.buf"), d.join("Data/Tb_vs_T.buf")).unwrap();
    let out = d.join("Plots/tb.drv");
    assert!(drs(&["plots", "save", d.join("Data/Tb_vs_T.buf").to_str().unwrap(), "--out", out.to_str().unwrap()]).status.success());
    let mut f = SetFile::read(&std::fs::read_to_string(&out).unwrap()).unwrap();
    let (_, g) = tb_set(false);
    (f.map_scale, f.map_range, f.regions) = (g.map_scale, g.map_range, g.regions);
    f.reference = Some("Data/tb_t20.dat".into());
    std::fs::write(&out, f.to_text()).unwrap();
    d
}

fn page_json(html: &str) -> serde_json::Value {
    let open = "<script id=\"drs-page-data\" type=\"application/json\">";
    let a = html.find(open).expect("the page's data element") + open.len();
    let b = a + html[a..].find("</script>").unwrap();
    serde_json::from_str(&html[a..b]).unwrap()
}

#[test]
fn drs_plots_show_writes_one_page() {
    let d = project("show");
    let set = d.join("Plots/tb.drv");
    let o = drs(&["plots", "show", set.to_str().unwrap(), "--browser", "--no-open"]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let page = d.join("Plots/tb.map.html");
    let html = std::fs::read_to_string(&page).unwrap();
    assert!(String::from_utf8_lossy(&o.stdout).starts_with(&format!("{}: the page of {}", page.display(), set.display())));
    // The JSON: the set's map, regions and members, the file and the core's version and date on it.
    let v = page_json(&html);
    assert_eq!((v["path"].as_str(), v["core"].as_str()), (set.to_str(), Some(env!("CARGO_PKG_VERSION"))));
    assert_eq!(v["written"].as_str().unwrap().len(), 10);
    assert_eq!((v["maps"]["map"]["scale"].as_str(), v["maps"]["difference"]["scale"].as_str()), (Some("log"), Some("sqrt")));
    assert_eq!(v["regions"][0]["blocks"]["map"]["columns"], serde_json::json!([30, 31, 32]));
    // The page's JSON holds both of the region's ranges (P65).
    assert_eq!((&v["regions"][0]["map_range"], &v["regions"][0]["difference_range"]), (&serde_json::json!({ "from": 3.9, "to": 4.4 }), &serde_json::json!({ "from": -40.0, "to": 60.0 })));
    assert_eq!(v["members"].as_array().unwrap().len(), 20);
    assert!(v["members"][0]["pattern"]["x"].is_array());
    // The bundle whole, `</` escaped; the shell's markers gone; nothing fetched.
    assert!(html.contains(&plots::VIEWER_BUNDLE.replace("</", "<\\/")));
    assert!(!html.contains(plots::DATA_MARKER) && !html.contains(plots::BUNDLE_MARKER));
    assert_eq!(html.matches("</script>").count(), 3, "the data, the bundle and the shell's own script, each closed once");
    let mut http: Vec<&str> = html.match_indices("http").map(|(i, _)| &html[i..(i + 30).min(html.len())]).collect();
    http.retain(|h| !h.starts_with("http://www.w3.org/2000/svg") && !h.starts_with("http://www.w3.org/1999/xlink"));
    assert_eq!(http, Vec::<&str>::new(), "no http reference other than the SVG namespaces");
    // Again over its own page; without patterns, smaller and the same maps.
    let lean = d.join("Plots/lean.html");
    assert!(drs(&["plots", "show", set.to_str().unwrap(), "--browser", "--no-open", "--no-patterns", "--out", lean.to_str().unwrap()]).status.success());
    let w = page_json(&std::fs::read_to_string(&lean).unwrap());
    assert!(w["members"].as_array().unwrap().iter().all(|m| m["pattern"].is_null()));
    assert_eq!(w["maps"], v["maps"]);
    assert!(std::fs::metadata(&lean).unwrap().len() < std::fs::metadata(&page).unwrap().len());
    assert!(drs(&["plots", "show", set.to_str().unwrap(), "--browser", "--no-open"]).status.success(), "over its own page");
    // Never over another file.
    let other = d.join("Plots/notes.html");
    std::fs::write(&other, "<p>mine</p>").unwrap();
    let o = drs(&["plots", "show", set.to_str().unwrap(), "--browser", "--no-open", "--out", other.to_str().unwrap()]);
    assert!(!o.status.success() && String::from_utf8_lossy(&o.stderr).contains("is not a page of drs plots show: nothing written"));
    assert_eq!(std::fs::read_to_string(&other).unwrap(), "<p>mine</p>");
    // The difference view needs a reference.
    let mut f = SetFile::read(&std::fs::read_to_string(&set).unwrap()).unwrap();
    f.reference = None;
    let bare = d.join("Plots/bare.drv");
    std::fs::write(&bare, f.to_text()).unwrap();
    let o = drs(&["plots", "show", bare.to_str().unwrap(), "--browser", "--no-open", "--view", "difference"]);
    assert!(!o.status.success() && String::from_utf8_lossy(&o.stderr).contains("names no reference member"));
    // One member: refused in words.
    let mut one = f.clone();
    one.members.truncate(1);
    let path = d.join("Plots/one.drv");
    std::fs::write(&path, one.to_text()).unwrap();
    let o = drs(&["plots", "show", path.to_str().unwrap(), "--browser", "--no-open"]);
    assert!(!o.status.success() && String::from_utf8_lossy(&o.stderr).contains("a map has one row per value of T, and this set has 1 member with one"), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(!d.join("Plots/one.map.html").exists());
}

#[test]
fn the_viewers_assets_are_present_and_consistent() {
    // SOURCES.txt: one "name sha256" line per source; the bundle's header names the same sources and hashes, in order.
    let sources: Vec<(String, String)> = plots::VIEWER_SOURCES
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let (n, h) = l.split_once(' ').unwrap_or_else(|| panic!("core/assets/viewer/SOURCES.txt: \"{l}\" is not \"name sha256\""));
            assert!(h.len() == 64 && h.chars().all(|c| c.is_ascii_hexdigit()), "core/assets/viewer/SOURCES.txt: {n}'s hash \"{h}\" is not a sha256");
            (n.to_string(), h.to_string())
        })
        .collect();
    assert!(!sources.is_empty(), "core/assets/viewer/SOURCES.txt names no source: regenerate the assets with the app's bundler (npm run bundle-viewer)");
    let header = plots::VIEWER_BUNDLE.lines().next().unwrap_or("");
    let at = header.find("Sources (app/ui/): ").unwrap_or_else(|| panic!("core/assets/viewer/map.js: its first line names no sources (\"Sources (app/ui/): …\"); regenerate it with the app's bundler"));
    let named: Vec<(String, String)> = header[at + "Sources (app/ui/): ".len()..]
        .trim_end_matches('.')
        .split("; ")
        .map(|s| {
            let (n, h) = s.split_once(" sha256 ").unwrap_or_else(|| panic!("core/assets/viewer/map.js: \"{s}\" in its header is not \"name sha256 hash\""));
            (n.to_string(), h.to_string())
        })
        .collect();
    assert_eq!(named, sources, "the bundle's header and SOURCES.txt name different sources or hashes: the assets are stale or edited by hand; regenerate them with the app's bundler");
    assert!(header.contains("Generated by app/ui/tools/bundle-viewer.mjs"), "core/assets/viewer/map.js is not the bundler's");
    for (n, _) in &sources {
        assert!(plots::VIEWER_BUNDLE.contains(&format!("\n// ---- {n} ----\n")), "core/assets/viewer/map.js lacks the part of {n} that its header names");
    }
    assert!(plots::VIEWER_BUNDLE.contains("function createMapViewer("), "core/assets/viewer/map.js holds no createMapViewer, which the page's shell mounts");
    // The shell: each marker once, where the page's data and the bundle go.
    for m in [plots::DATA_MARKER, plots::BUNDLE_MARKER] {
        assert_eq!(plots::PAGE_SHELL.matches(m).count(), 1, "core/assets/viewer/page.html must hold the marker {m} once");
    }
    assert!(plots::PAGE_SHELL.contains("<script id=\"drs-page-data\" type=\"application/json\"><!--DRS_PAGE_DATA--></script>"));
}

#[test]
fn drs_plots_region_counts_from_one_and_says_the_indices() {
    let d = project("region");
    let set = d.join("Plots/tb.drv");
    let o = drs(&["plots", "region", set.to_str().unwrap(), "--name", "m1 at 20-26 K"]);
    let out = String::from_utf8_lossy(&o.stdout);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(out.starts_with("# region \"m1 at 20-26 K\" on the map, log scale, x in 2theta\n"), "{out}");
    assert!(out.contains("rows 9–15 of 19 (T = 20 to 26 K), columns 31–33 of 301 (2theta 19 to 19.2); counted from 1, the indices from 0 being rows 8–14, columns 30–32\n"), "{out}");
    let j: plots::MapBlock = serde_json::from_slice(&drs(&["plots", "region", set.to_str().unwrap(), "--name", "m1 at 20-26 K", "--json"]).stdout).unwrap();
    assert_eq!((j.rows[0], j.columns.clone(), j.counted_from), (8, vec![30, 31, 32], 0));
}
