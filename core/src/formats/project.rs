//! The files of a project by what they are (app request 9; L066; W5, L080, L081): every file of a FullProf run folder
//! named by its role, with the source that says what it is; unknown when no source does; never a "table" by content
//! when the extension is FullProf's. Data files are the registry's (`identify`): `describe` asks it when the file is
//! none of these.
//!
//! The order of the tests: the file's name where it is this project's own (`.pcr.in`, `.run.log`); content
//! signatures, each found in every file of its kind in `manuals/fullprof/runs`, the tutorials and FP_tests (2649
//! files, test `file_kinds`); then the extension. The sources: FullProf's manual (July 2001; Input files p. 8, Output
//! files pp. 9–12, and the pages of the flags), the FullProf News (`manuals/fullprof/news`), the tutorials and the
//! manuals of the suite's other programs; "not documented" where only the content says what a file is.

use super::Probe;

/// What a file is for in a project.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Role {
    /// A control file: a refinement is run from it.
    Control,
    /// Read by a run (or by another program of the suite).
    Input,
    /// Written by a run.
    Output,
    /// Documents, scripts, others' files.
    Other,
    /// Nothing says what it is.
    Unknown,
}

/// A kind of project file: its id, readable name, role and the source that says what it is.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Kind {
    pub id: &'static str,
    pub name: &'static str,
    pub role: Role,
    pub source: &'static str,
}

const fn k(id: &'static str, name: &'static str, role: Role, source: &'static str) -> Kind {
    Kind { id, name, role, source }
}

use Role::{Control as C, Input as I, Other as O, Output as W};

pub const CONTROL: Kind = k("control", "diff-ref-struct control file (.drc)", C, "L070 (core/src/control)");
pub const RESULTS: Kind = k("results", "Results of a refinement", W, "its first line, as core/src/results.rs writes it (R28, L092)");
pub const PLOT_SET: Kind = k("plot_set", "diff-ref-struct plot set (.drv)", I, "its line format = \"drs plot set 1\", as core/src/plots.rs writes it (PL12, PL24)");
pub const PUBLICATION_SET: Kind = k("publication_set", "diff-ref-struct publication set (.drp)", I, "its line format = \"drs-publish set 1\", as publish/src/set.rs writes it (the app's request 37)");
pub const OUR_BACKGROUND: Kind = k("drs_background", "diff-ref-struct background points (.drb)", I, "L078 (core/src/background.rs)");
pub const PCR: Kind = k("fullprof_pcr", "FullProf control file (.pcr)", C, "manual p. 8, Input files: CODFIL.pcr");
pub const PCR_NEW: Kind = k("fullprof_pcr_new", "FullProf control file written anew (.new)", W, "manual p. 8, Input files: CODFIL.new; p. 72 (Pcr = 2)");
pub const PCR_ORIGINAL: Kind = k("fullprof_pcr_original", "FullProf control file as the tutorial gave it (.pcr.in)", I, "runs/README.md (runs/run_fullprof.py)");
pub const BACKGROUND: Kind = k("fullprof_background", "FullProf background file (.bac)", I, "manual p. 8, Input files: FILE.bac; p. 73 (Ppl = 2)");
pub const DIFFERENCE: Kind = k("fullprof_difference", "FullProf difference pattern (.bac, Ppl = 3)", W, "manual p. 73 (Ppl = 3)");
pub const REFLECTIONS: Kind = k("fullprof_reflections", "FullProf reflection list (.hkl)", W, "manual p. 11, Output files: CODFIL.hkl (Hkl)");
pub const REFLECTIONS_HKL2: Kind = k("fullprof_reflections_hkl_minus2", "FullProf reflection list, Hkl = −2 (.hkl)", W, "not documented (content; the .pcr has Hkl = −2)");
pub const PHASE_REFLECTIONS: Kind = k("fullprof_phase_reflections", "FullProf reflections of a phase, profile matching (CODFILn.hkl)", I, "manual p. 9, Input files: CODFILn.hkl (Irf); p. 90 (Jbt = 2, Irf 0 → 2)");
pub const MAGNETIC_REFLECTIONS: Kind = k("fullprof_magnetic_reflections", "FullProf observed and calculated F² of a phase (CODFILn.hkl)", W, "manual p. 10, Output files: CODFILn.hkl; p. 90 (Irf ≤ 0, Isy)");
pub const BAD_REFLECTIONS: Kind = k("fullprof_bad_reflections", "FullProf reflections fitted worst (_bad_refl.hkl)", W, "not documented (content)");
pub const SHELX_HKL: Kind = k("shelx_hklf4", "single-crystal intensities, SHELX HKLF 4 (.hkl)", I, "8-FullProf-Tutorial-RigidBodies/tutorial-rgb.inf");
pub const INTEGRATED: Kind = k("fullprof_integrated_intensities", "FullProf integrated intensities (.int)", I, "manual p. 9, Input files: CODFIL.int; p. 12, Output files: CODFILn.int; p. 138");
pub const RESOLUTION: Kind = k("fullprof_resolution", "FullProf instrumental resolution file (.irf)", I, "manual p. 9, Input files: MYRESOL.irf");
pub const LISTING: Kind = k("fullprof_listing", "FullProf listing (.out)", W, "manual p. 9, Output files: CODFIL.out");
pub const SUMMARY: Kind = k("fullprof_summary", "FullProf summary (.sum)", W, "manual p. 10, Output files: CODFIL.sum");
pub const LOG: Kind = k("fullprof_log", "FullProf log of a run (.log)", W, "manual p. 7 (fp2k CODFIL DATFIL LOGFIL)");
pub const RUN_LOG: Kind = k("fullprof_run_log", "FullProf screen output kept by the runs scripts (.run.log)", W, "runs/run_fullprof.py");
pub const PROFILE: Kind = k("fullprof_profile", "FullProf observed and calculated profile (.prf)", W, "manual p. 9, Output files: CODFIL.prf");
pub const INTEGRATED_PROFILE: Kind = k("fullprof_integrated_profile", "FullProf observed and calculated integrated intensities (.prf)", W, "News 19 June 2020");
pub const SYMMETRY: Kind = k("fullprof_symmetry", "FullProf symmetry listing (.sym)", W, "manual p. 10, Output files: CODFIL.sym");
pub const FOURIER: Kind = k("fullprof_fourier", "FullProf structure factors for Fourier maps (.fou)", W, "manual p. 10, Output files: CODFIL.fou (Fou)");
pub const FOURIER_INPUT: Kind = k("fullprof_fourier_input", "FullProf input for (G)FOURIER (.inp)", W, "manual p. 10, Output files: CODFILn.inp (Fou = 3, 4)");
pub const SHELX_TEMPLATE: Kind = k("fullprof_shelx_template", "FullProf input for SHELXS (.ins)", W, "manual p. 10, Output files: CODFILn.ins (Fou = 2)");
pub const MICROSTRUCTURE: Kind = k("fullprof_microstructure", "FullProf size and strain file (.mic)", W, "manual pp. 11–12, Output files: CODFILn.mic");
pub const DISTANCES: Kind = k("fullprof_distances", "FullProf distances, angles and bond valences (.dis)", W, "manual p. 11, Output files: CODFILn.dis");
pub const PHASE_PROFILE: Kind = k("fullprof_phase_profile", "FullProf calculated profile of a phase (.sub)", W, "manual p. 12, Output files: CODFILn.sub");
pub const ATOMS: Kind = k("fullprof_atoms", "FullProf atoms in the cell (.atm)", W, "manual p. 12, Output files: CODFILn.atm");
pub const RPA: Kind = k("fullprof_parameters_appended", "FullProf refined parameters, appended run after run (.rpa)", W, "manual p. 9, Output files: CODFIL.rpa");
pub const SAV: Kind = k("fullprof_reflections_interval", "FullProf reflections between two angles (.sav)", W, "manual p. 11, Output files: CODFIL.sav");
pub const SIM: Kind = k("fullprof_simulated", "FullProf simulated pattern (.sim)", W, "manual p. 12, Output files: CODFIL.sim");
pub const SCH: Kind = k("fullprof_schakal", "FullProf input for SCHAKAL or STRUPLO (.sch)", W, "manual p. 12, Output files: CODFILn.sch");
pub const SHP: Kind = k("fullprof_shape", "FullProf numerical peak shape (.shp)", I, "manual p. 9, Input files: CODFIL.shp");
pub const COR: Kind = k("fullprof_corrections", "FullProf intensity corrections (.cor)", I, "manual p. 9, Input files: CODFIL.cor");
pub const HLP: Kind = k("fullprof_constraint_help", "FullProf distance and angle restraints to paste (DCONSTRn.hlp)", W, "manual p. 11, Output files: DCONSTRn.hlp");
pub const STUDIO: Kind = k("fullprof_studio", "FullProf Studio structure (.fst)", W, "News 14 July 2004");
pub const ORTEP: Kind = k("fullprof_ortep", "FullProf atoms for ORTEP (CODFILn.ort)", W, "manual p. 39 (Jbt = 4: x#.ortep, reduced co-ordinates and temperature parameters for ORTEP; believed its present name CODFILn.ort: its ATOM and BETA lines)");
pub const MCIF: Kind = k("fullprof_mcif", "magnetic CIF written by FullProf (.mcif)", W, "News 19 May 2016");
pub const FULLPROF_CIF: Kind = k("fullprof_cif", "CIF written by FullProf", W, "News 2002 (Rpa = −1)");
pub const PROFILE_CIF: Kind = k("fullprof_profile_cif", "powder pattern CIF written by FullProf (_prof.cif)", W, "News 27 November 2019");
pub const SIMPLE_CIF: Kind = k("fullprof_simple_cif", "simple CIF written by FullProf (_parent.cif)", W, "not documented (content)");
pub const DISTANCES_CIF: Kind = k("fullprof_distances_cif", "distances and angles as CIF, written by FullProf (_dis.cif)", W, "News 22 March 2008");
pub const STRUCTURE_CIF: Kind = k("cif", "crystal structure (CIF)", I, "Hall, Allen & Brown, Acta Cryst. A47, 655 (1991)");
pub const VESTA: Kind = k("vesta", "VESTA file (.vesta)", W, "News 18 May 2015");
pub const VESTA_GRID: Kind = k("vesta_grid", "VESTA volumetric grid, binary (.ggrid, .pgrid)", W, "News 22 October 2014");
pub const MIC_BINARY: Kind = k("fullprof_microstructure_binary", "FullProf size or strain diagram, binary (.bin)", W, "News April 2002");
pub const STRAIN: Kind = k("fullprof_strain", "FullProf strain along directions (.str)", W, "News 14 May 2002");
pub const SIZE: Kind = k("fullprof_size", "FullProf size along directions (.siz)", W, "News 14 May 2002");
pub const FULLPROF_CFL: Kind = k("fullprof_cfl", "structure for Bond_Str written by FullProf (.cfl)", W, "Y2O3 tutorial p. 21; docs/bond_str.pdf p. 1");
pub const CRYSFML_CFL: Kind = k("crysfml_cfl", "CrysFML input file (.cfl)", I, "docs/bond_str.pdf p. 1");
pub const BVS: Kind = k("bond_valence_summary", "bond-valence summary (_sum.bvs)", W, "docs/bond_str.pdf p. 3");
pub const RESTRAINTS: Kind = k("fullprof_restraints", "FullProf restraints to paste into a .pcr (.tpcr)", W, "News 20 November 2012");
pub const MOL_TPCR: Kind = k("mol_tpcr_output", "Mol_tpcr rigid bodies to paste into a .pcr (.tpcr)", W, "docs/FullProf-Restraints-RigidBodies.pdf p. 15");
pub const SA_PROFILE: Kind = k("fullprof_sa_profile", "FullProf profile for simulated annealing (.spr)", W, "News 7 July 2010");
pub const SEQUENTIAL: Kind = k("fullprof_sequential", "FullProf results of a sequential refinement (.seq)", W, "News 10 February 2004");
pub const SA_SOLUTION: Kind = k("fullprof_sa_solution", "FullProf simulated-annealing solution (.so1, .so2)", W, "News 27 July 2006 (.so1)");
pub const SA_LOG: Kind = k("fullprof_sa_log", "FullProf simulated-annealing temperature log (.anl)", W, "not documented (content)");
pub const SA_CONFIGURATIONS: Kind = k("fullprof_sa_configurations", "FullProf simulated-annealing configurations (.cfg)", W, "not documented (content)");
pub const BACKGROUND_POINTS: Kind = k("background_points", "background points (.bgr)", I, "FAULTS manual p. 32");
pub const WINPLOTR_BACKGROUND: Kind = k("winplotr_background_points", "background points exported as x y (.xy)", I, "not documented (content)");
pub const BUFFER: Kind = k("buffer_file", "list of file names (.buf)", I, "WinPLOTR manual p. 15");
pub const PEAKS: Kind = k("peak_list", "peak positions, intensities and background (.aps)", W, "FullProfAPP manual (Export Peaks)");
pub const DICVOL_INPUT: Kind = k("dicvol_input", "DICVOL input (.dic)", I, "not documented (content: WinPLOTR's header)");
pub const DICVOL_OUTPUT: Kind = k("dicvol_output", "DICVOL listing (.ind, .ord)", W, "not documented (content: DICVOL's header)");
pub const DICVOL_CONDENSED: Kind = k("dicvol_condensed", "DICVOL condensed results (.ucp)", W, "not documented (content: DICVOL's header)");
pub const DICVOL_CELL: Kind = k("dicvol_cell", "DICVOL cell (.cry)", W, "not documented (content: DICVOL's header)");
pub const DICVOL_SUMMARY: Kind = k("dicvol_summary", "DICVOL cell summary (_dicvol_cell.sum)", W, "not documented (content: WDICVOL04 line)");
pub const TREOR: Kind = k("treor_output", "TREOR listing (.cnd)", W, "not documented (content: TREOR's header)");
pub const CHECK_GROUP: Kind = k("check_group_output", "CheckGroup listing (.spg)", W, "not documented (content: CheckGroup's banner)");
pub const FORMAL_CHARGES: Kind = k("formal_charges_output", "Formal_Charges listing (.fch)", W, "not documented (content: Formal_Charges' banner)");
pub const ZMATRIX: Kind = k("zmatrix_fhz", "molecule as a Z-matrix (.fhz)", I, "docs/FullProf-Restraints-RigidBodies.pdf p. 5");
pub const DATARED: Kind = k("datared_input", "DataRed control file (.red)", I, "News 22 February 2022");
pub const FAULTS: Kind = k("faults_input", "FAULTS control file (.flts)", C, "FAULTS manual p. 7");
pub const MSF: Kind = k("magnetic_structure_factors", "magnetic structure factors for Fourier (.msf)", I, "News 12 February 2012");
pub const GET_LOG: Kind = k("get_log_results", "Get_Log summary of FullProf logs (.res)", W, "not documented (content: Get_Log's banner)");
pub const INFO: Kind = k("text_information", "text (.inf)", O, "manual p. 8 (fp2k.inf)");
pub const TEXT: Kind = k("text", "text", O, "the extension");
pub const SCRIPT: Kind = k("script", "script", O, "the extension");
pub const PDF: Kind = k("pdf", "document (PDF)", O, "the extension");
pub const UNKNOWN: Kind = k("unknown", "unknown", Role::Unknown, "nothing says what it is");

/// The items of a powder pattern's loop in a CIF (pdCIF).
const PATTERN_ITEMS: [&str; 5] = ["_pd_meas_counts_total", "_pd_meas_intensity_total", "_pd_proc_intensity_total", "_pd_proc_intensity_net", "_pd_calc_intensity_total"];

/// The extensions of intensity data (FullProf's manual p. 8: FILE.dat, "dat, uxd, acq"; Ins = 10: x y σ; GSAS).
const DATA_EXTENSIONS: [&str; 9] = ["dat", "uxd", "acq", "xy", "xys", "xye", "xyd", "gss", "gsa"];

/// A line among the first `n` that starts (after blanks) with `s`.
fn starts(head: &str, n: usize, s: &str) -> bool {
    head.lines().take(n).any(|l| l.trim_start().starts_with(s))
}

fn line(head: &str, i: usize) -> &str {
    head.lines().nth(i).unwrap_or("")
}

/// What a file is, when it is not data: by its name where it is this project's own, by content, then by extension.
/// None when nothing says (the registry's readers may still read it as data; `describe` then says "unknown").
pub fn kind(p: &Probe) -> Option<Kind> {
    let ext = p.extension();
    let base = p.name.rsplit(['/', '\\']).next().unwrap_or(p.name).to_lowercase();
    // Binary files: by their own header.
    if !p.is_text() {
        let b = p.bytes;
        if b.windows(8).any(|w| w == b"!Job Npr") {
            return Some(PCR);
        }
        return match ext.as_str() {
            "ggrid" | "pgrid" => Some(VESTA_GRID),
            "bin" if b.windows(13).take(64).any(|w| w == b"Strain Diagra" || w == b"-Size Diagram") => Some(MIC_BINARY),
            "pdf" => Some(PDF),
            _ => None,
        };
    }
    let text = p.text();
    let head: String = text.chars().take(8192).collect();
    let first = line(&head, 0).trim_start();
    let second = line(&head, 1).trim_start();
    // This project's own names.
    if base.ends_with(".pcr.in") {
        return Some(PCR_ORIGINAL);
    }
    if base.ends_with(".run.log") || base == "run.log" {
        return Some(RUN_LOG);
    }
    // Content. A results file first: it holds its control file inside (`_drs_control_file`), which would name it one.
    if first.starts_with("# The results of a refinement, written by diff-ref-struct") {
        return Some(RESULTS);
    }
    // A plot set (PL24): by its format line, wherever the file lies and whatever its name.
    if head.lines().any(|l| l.trim_start().starts_with("format") && l.contains(&format!("\"{}\"", crate::plots::SET_FORMAT))) {
        return Some(PLOT_SET);
    }
    // A publication set (the app's request 37): by its format line, as a plot set.
    if head.lines().any(|l| l.trim_start().starts_with("format") && l.contains("\"drs-publish set 1\"")) {
        return Some(PUBLICATION_SET);
    }
    if head.lines().any(|l| l.trim_start().starts_with("format") && l.contains("\"diff-ref-struct ")) {
        return Some(if head.contains("diff-ref-struct background") { OUR_BACKGROUND } else { CONTROL });
    }
    if text.lines().any(|l| l.starts_with("!Job Npr")) {
        return Some(if ext == "new" { PCR_NEW } else { PCR });
    }
    if head.contains("PROGRAM FullProf") {
        if head.contains("L O G   F I L E") {
            return Some(LOG);
        }
        if head.lines().filter(|l| !l.trim().is_empty()).take(3).all(|l| l.starts_with(" =>")) {
            return Some(RUN_LOG);
        }
        return Some(match ext.as_str() {
            "sum" => SUMMARY,
            "log" => LOG,
            _ => LISTING,
        });
    }
    if head.contains("====== PROGRAM BOND_STR ======") {
        return Some(DISTANCES);
    }
    if first.starts_with("Subroutine Calc_BVS") {
        return Some(BVS);
    }
    if first.contains("CELL:") && first.contains("SPGR:") && head.lines().take(12).any(|l| l.contains("Yobs") && l.contains("Ycal")) {
        return Some(PROFILE);
    }
    if second.starts_with("Sinthet/lamb") && second.contains("Gobs") {
        return Some(INTEGRATED_PROFILE);
    }
    // Reflections.
    if first.starts_with("Pattern#") && first.contains("Phase No") && first.contains("Lambda:") {
        return Some(PHASE_REFLECTIONS);
    }
    if (first.starts_with("Pattern#") && first.contains("Phase No.")) || first.starts_with(">Phase no.:") {
        return Some(REFLECTIONS);
    }
    if first.starts_with(">") && first.contains(" h ") && first.contains("F2") {
        return Some(REFLECTIONS_HKL2);
    }
    if first.starts_with("LIST OF REFLECTIONS WITH DELTA/SIGMA") {
        return Some(BAD_REFLECTIONS);
    }
    if ext == "hkl" && first.starts_with("! Phase No.:") {
        return Some(MAGNETIC_REFLECTIONS);
    }
    if ext == "fou" && first.starts_with("Phase No.:") {
        return Some(FOURIER);
    }
    if head.contains("SYMMETRY INFORMATION ON PHASE") {
        return Some(SYMMETRY);
    }
    if first.starts_with("! MICRO-STRUCTURAL ANALYSIS FROM FULLPROF") {
        return Some(MICROSTRUCTURE);
    }
    if first.starts_with("! Template file to be used with the program FOURIER") {
        return Some(FOURIER_INPUT);
    }
    if ext == "ins" && first.starts_with("TITL") {
        return Some(SHELX_TEMPLATE);
    }
    if first.starts_with("TITLE Atom coordinates of phase no:") {
        return Some(ATOMS);
    }
    if starts(&head, 400, "SPACEG") && starts(&head, 400, "CELL") {
        return Some(STUDIO);
    }
    if second.starts_with("# Magnetic CIF file generated by FullProf") || head.lines().take(3).any(|l| l.contains("Magnetic CIF file generated by FullProf")) {
        return Some(MCIF);
    }
    if second.contains("FullProf-generated Simple CIF") {
        return Some(SIMPLE_CIF);
    }
    if second.contains("FullProf-generated CIF") {
        return Some(FULLPROF_CIF);
    }
    if head.lines().take(5).any(|l| l.starts_with("# Powder diffraction pattern number")) {
        return Some(PROFILE_CIF);
    }
    if ext == "cif" && head.contains("UNIT CELL INFORMATION") && head.contains("MOLECULAR GEOMETRY") {
        return Some(DISTANCES_CIF);
    }
    if first.starts_with("#VESTA_FORMAT_VERSION") {
        return Some(VESTA);
    }
    if first.starts_with("FILE containing strain along any direction") {
        return Some(STRAIN);
    }
    if first.starts_with("FILE containing size") {
        return Some(SIZE);
    }
    if head.lines().take(3).any(|l| l.contains("CFL-file generated from FullProf")) {
        return Some(FULLPROF_CFL);
    }
    if first.starts_with("FILE with lines for soft distance and angle constraints") {
        return Some(RESTRAINTS);
    }
    if head.contains("PROGRAM MOL_TPCR") {
        return Some(MOL_TPCR);
    }
    if first.starts_with("! Npoi") {
        return Some(SA_PROFILE);
    }
    if head.lines().take(5).any(|l| l.trim_start().starts_with("NEW_REFINEMENT")) {
        return Some(SEQUENTIAL);
    }
    if head.lines().take(5).any(|l| l.trim_start().starts_with("!  Data for PHASE number:")) {
        return Some(SA_SOLUTION);
    }
    if first.starts_with("! Simulated Anneling Run") || first.starts_with("! Simulated Annealing Run") {
        return Some(SA_LOG);
    }
    if first.starts_with("! NT") && first.contains("R-optimum") {
        return Some(SA_CONFIGURATIONS);
    }
    if ext == "bac" {
        if first.contains("Difference pattern") {
            return Some(DIFFERENCE);
        }
        return Some(BACKGROUND);
    }
    if ext == "sub" && first.contains("Phase No") {
        return Some(PHASE_PROFILE);
    }
    if head.lines().take(3).any(|l| l.trim_start().starts_with('(')) && ext == "int" {
        return Some(INTEGRATED);
    }
    if first.starts_with("! Data file") && starts(&head, 3, "! Format") {
        return Some(if ext == "bgr" || starts(&head, 4, "! Number of background points") { BACKGROUND_POINTS } else { WINPLOTR_BACKGROUND });
    }
    if first.starts_with("# Input file for WDICVOL") {
        return Some(DICVOL_INPUT);
    }
    if head.contains("DICVOL CONDENSED RESULTS") {
        return Some(DICVOL_CONDENSED);
    }
    if head.contains("DICVOL input parameters") {
        return Some(DICVOL_OUTPUT);
    }
    if first.starts_with("Dicvol04 solution") {
        return Some(DICVOL_CELL);
    }
    if head.contains("WDICVOL04") && ext == "sum" {
        return Some(DICVOL_SUMMARY);
    }
    if head.contains("TREOR input parameters") {
        return Some(TREOR);
    }
    if head.contains("PROGRAM CHECK_GROUP") {
        return Some(CHECK_GROUP);
    }
    if head.contains("FORMAL CHARGES") {
        return Some(FORMAL_CHARGES);
    }
    if head.contains("Program Get_Log") {
        return Some(GET_LOG);
    }
    // A CIF with a powder pattern is data (the registry's); one with a cell and no pattern, a structure.
    if ext == "cif" {
        if PATTERN_ITEMS.iter().any(|t| text.contains(t)) {
            return None;
        }
        if text.contains("_cell_length_a") {
            return Some(STRUCTURE_CIF);
        }
    }
    // By extension.
    Some(match ext.as_str() {
        "drc" => CONTROL,
        "drb" => OUR_BACKGROUND,
        "drv" => PLOT_SET,
        "drp" => PUBLICATION_SET,
        "pcr" => PCR,
        "new" => PCR_NEW,
        "out" => LISTING,
        "sum" => SUMMARY,
        "prf" => PROFILE,
        "sym" => SYMMETRY,
        "mic" => MICROSTRUCTURE,
        "irf" => RESOLUTION,
        "hkl" if text.lines().take(3).all(|l| l.split_whitespace().count() >= 5 && l.split_whitespace().all(|t| t.parse::<f64>().is_ok())) => SHELX_HKL,
        "hkl" => REFLECTIONS,
        "int" => INTEGRATED,
        "fou" => FOURIER,
        "dis" => DISTANCES,
        "sub" => PHASE_PROFILE,
        "atm" => ATOMS,
        "rpa" => RPA,
        "sav" => SAV,
        "sim" => SIM,
        "sch" => SCH,
        "shp" => SHP,
        "cor" => COR,
        "hlp" => HLP,
        "fst" => STUDIO,
        // A rigid-body phase's atoms for ORTEP (app request 47); FullProf writes it empty for other phases.
        "ort" => ORTEP,
        "mcif" => MCIF,
        "vesta" => VESTA,
        "ggrid" | "pgrid" => VESTA_GRID,
        "str" => STRAIN,
        "siz" => SIZE,
        "cfl" => CRYSFML_CFL,
        "bvs" => BVS,
        "spr" => SA_PROFILE,
        "seq" => SEQUENTIAL,
        "so1" | "so2" | "so3" => SA_SOLUTION,
        "anl" => SA_LOG,
        "bgr" => BACKGROUND_POINTS,
        "buf" => BUFFER,
        "aps" => PEAKS,
        "dic" => DICVOL_INPUT,
        "ind" | "ord" => DICVOL_OUTPUT,
        "ucp" => DICVOL_CONDENSED,
        "cry" => DICVOL_CELL,
        "cnd" => TREOR,
        "fhz" => ZMATRIX,
        "red" => DATARED,
        "flts" => FAULTS,
        "msf" => MSF,
        "inf" => INFO,
        "txt" | "md" => TEXT,
        "py" | "sh" | "bat" => SCRIPT,
        "pdf" => PDF,
        // FullProf's own extensions that nothing documents: named unknown, never read as a table.
        "mlf" | "tpcr" | "cfg" | "inp" | "ins" | "res" | "pyc" => UNKNOWN,
        _ => return None,
    })
}

/// What a file is: a kind of project file, data (the registry's format), or unknown.
#[derive(Clone, Debug, PartialEq)]
pub enum Description {
    Project(Kind),
    /// A data file: the registry's format id and name.
    Data { id: &'static str, name: &'static str },
    Unknown,
}

impl Description {
    pub fn id(&self) -> &'static str {
        match self {
            Description::Project(k) => k.id,
            Description::Data { .. } => "data",
            Description::Unknown => "unknown",
        }
    }
}

/// What a file is, given its path relative to the project's folder (the app's request 31): the file dictionary's kind
/// when its place and name fit one, unless that kind is the measured data (`Data/{original}`, which any name fits:
/// the content says what the data are); else `describe`. A results file (`Results/… results.cif`, a powder CIF the
/// registry reads as a pattern) is so named a result, never data.
/// The kind is the dictionary's id and name, its role mapped onto `Role` (a version is a control file kept, the app's
/// state is other), its source "the file dictionary".
pub fn describe_in(rel: &str, p: &Probe) -> Description {
    match crate::formats::files::kind_of(rel) {
        Some(k) if k.role != "data" => {
            let role = match k.role {
                "control" | "version" => Role::Control,
                "input" => Role::Input,
                "output" => Role::Output,
                _ => Role::Other,
            };
            Description::Project(Kind { id: k.id, name: k.name, role, source: "the file dictionary (ecosystem/data/files.yml), by the file's place and name" })
        }
        _ => describe(p),
    }
}

/// What a file is when no project folder is given (the app's request 31 on the command line): `describe_in` with the
/// file's path taken from each folder above it, its own first, up to four levels; the first that the file dictionary
/// names by place and name (other than the measured data) decides. A results file under `Results/` is so a result
/// wherever the project lies. Else `describe`.
pub fn describe_path(path: &std::path::Path, p: &Probe) -> Description {
    let parts: Vec<String> = path.components().filter_map(|c| match c {
        std::path::Component::Normal(s) => Some(s.to_string_lossy().into_owned()),
        _ => None,
    }).collect();
    for k in 1..=parts.len().min(4) {
        let rel = parts[parts.len() - k..].join("/");
        if let Some(kd) = crate::formats::files::kind_of(&rel) {
            if kd.role != "data" {
                return describe_in(&rel, p);
            }
        }
    }
    describe(p)
}

/// What a file is: `kind`, else data when one of the registry's readers leads clearly, else unknown.
pub fn describe(p: &Probe) -> Description {
    if let Some(k) = kind(p) {
        return if k.id == "unknown" { Description::Unknown } else { Description::Project(k) };
    }
    match super::identify(p).format() {
        Some(c) => Description::Data { id: c.id, name: c.name },
        // FullProf's FILE.dat is intensity data whatever its layout (manual p. 8); the .pcr's Ins says which.
        None if DATA_EXTENSIONS.contains(&p.extension().as_str()) => Description::Data { id: "unrecognised", name: "intensity data (the layout is not recognised alone: the .pcr's Ins says it)" },
        None => Description::Unknown,
    }
}
