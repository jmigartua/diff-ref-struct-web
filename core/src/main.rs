//! `drs`, the command line of the diff-ref-struct core. Every core feature has a path from here.

use drs::cell::Cell;
use drs::mat3::Mat3;
use drs::sgtable::SETTINGS;
use drs::site::Site;
use drs::spacegroup::SpaceGroup;
use std::process::ExitCode;

const USAGE: &str = "\
drs — diff-ref-struct core

USAGE:
  drs cell <a> <b> <c> <alpha> <beta> <gamma> [--hkl <h> <k> <l>] [--lambda <Å>]
      Metric information of a unit cell (Å, degrees), as FullProf prints it:
      volume, reciprocal cell, metric tensors, Cartesian frame, Busing–Levy B.
      With --hkl: d(hkl), and 2θ when --lambda is given.
  drs sg <symbol> [--all]
      A space group and its symmetry operators. <symbol> is a number 1–230, a Hermann–Mauguin symbol
      (\"P 21/c\", \"Pnma\", \"F d -3 m:2\" for a setting choice) or a Hall symbol (\"-P 2ybc\").
      Prints one operator per rotation and the centring vectors; --all prints every operator. Each operator with its
      symmetry element as International Tables lists it (P09): the symbol with the sense of a 3-, 4- or 6-fold axis
      (2_1, 4_2-, m, c, n, d, -3+ …), the screw or glide component and the location (x,1/4,z; 0,y,1/4; x,x,x).
  drs sg --list
      The 530 settings: Hall number, number, choice, Hermann–Mauguin and Hall symbols.
  drs site <symbol> <x> <y> <z> [--tol <t>]
      A position in a space group: site-symmetry operators and point group, multiplicity, Wyckoff letter
      and site symbol, the equivalent positions in the cell, and the constraints the site puts on x, y, z
      and on the anisotropic displacement parameters (betas; the U_ij of a CIF obey the same relations).
      Coordinates as decimals or fractions (1/4); --tol is the tolerance in fractions (default 0.001).
  drs hkl <symbol> <a> <b> <c> <alpha> <beta> <gamma> --dmin <d> [--dmax <d>] [--lambda <Å>]
  drs hkl <symbol> <a> <b> <c> <alpha> <beta> <gamma> --lambda <Å> --2theta <min> <max>
      The reflections of a phase as FullProf lists them: one per set of equivalents of the Laue class,
      systematic absences from the operators, sorted by increasing 2θ (decreasing d): h k l, multiplicity,
      d, and 2θ when a wavelength is given. A '*' marks a representative not chosen by FullProf's rule
      (non-standard settings only).
  drs sf <file.cif> --dmin <d> [--xray <λ or anode> | --neutron] [--no-anomalous] [--dmax <d>]
  drs sf <symbol> <a> <b> <c> <alpha> <beta> <gamma> --atom <label> <species> <x> <y> <z> <Biso> [<occ>] ...
         --dmin <d> [--xray <λ or anode> | --neutron] [--no-anomalous] [--dmax <d>]
      Structure factors F(hkl) = Σ Occ·f(s)·exp(−B s²)·Σ_ops T·exp(2πi h·(Rx+t)), with FullProf's conventions, for
      the unique reflections with d ≥ dmin: h k l, multiplicity, d, |F|, F² (F² as FullProf prints it in its .hkl
      files: all operators, f′ and f″ included, no overall B). The structure comes from a CIF (cell, symmetry
      operators or Hall or H-M symbol, atom sites, U or B, anisotropic U_ij or B_ij) or from the command line: a
      space-group symbol read as FullProf reads it (\"F d -3 m\" is origin choice 2), six cell parameters, and one
      --atom per atom with its chemical occupancy (default 1; FullProf's Occ = occ·m/M is printed). Species are
      elements or ions as the tables write them (\"Fe\", \"O-2\", \"Ti+4\"). X-rays by default, Cu Kα1; --xray takes
      a wavelength in Å or an anode (Cr Fe Cu Mo Ag Co Ni, Kα1); f′, f″ are those of the anode nearest to λ.
      Neutrons: scattering lengths in 10⁻¹² cm, F in the same unit. Tables: CrysFML's (FullProf's).
  drs calc <file.pcr|file.drc> [--data <file>] [--out <file>] [--fullprof] [--phases] [--cumulative] [--shares]
      The calculated powder pattern of a FullProf control file, imported into the default model: the .pcr's values are
      translated wherever FullProf's rules differ from the papers' and an exact translation exists, and FullProf's peak
      window and reflections past the data are kept (L050). With --fullprof, the values as written under all of
      FullProf's rules, as FullProf calculates them. Reads the .pcr and the
      files it names (data, resolution, .hkl for profile matching), and writes 2θ, Ycal, Backg and, when the data
      file is found (or given with --data), Yobs and σ (to --out, or the terminal), with Rp, Rwp, Rexp and χ² on the excluded-free points.
      Constant-wavelength X-rays or neutrons; profiles Npr 0, 1, 5, 7; FCJ and Bérar–Baldinozzi asymmetry;
      backgrounds by points, polynomial or Chebyshev; resolution files (Res 1–4). --phases adds each phase's curve
      (without the background) as a column. With the data: --cumulative adds Σw(y − yc)² from the first point as a
      fraction of the whole, and √w(y − yc) ((y − yc)/σ under variance weights); --shares the background's and each
      phase's share of Ycal at each point (W7).
  drs reflections <file.pcr|file.drc> [--data <file>] [--fullprof] [--forbidden]
      Each reflection of the pattern (W8): phase, wavelength, h k l, multiplicity, 2θ, d, FWHM, η, |F|², I calc, I obs
      partitioned by the model as Rietveld's method shares the counts above the background (not a measurement where
      reflections overlap), its σ through the same partition, (I obs − I calc)/σ, and for a mode phase whose parent is
      known (mode_parent_group, mode_parent_basis) its indices in the parent's cell and its k-vector. --forbidden: the
      reflections each phase's space group forbids within the data's range, with the condition each breaks.
  drs simulate <file.pcr> [--data <file>] [--out <file>] [--seed <n>] [--from <2θ> --to <2θ> --step <Δ2θ>] [--fullprof]
      A simulated observed pattern: the calculated pattern of the .pcr with Poisson noise (the counts set by the
      scale factors), at the data points or at Thmin, Step, Thmax of the .pcr (or --from/--to/--step); writes
      2θ, Yobs and σ = √Yobs in the X, Y, σ format (Ins=10).
  drs mic <file.pcr> [--out <file>] [--fullprof]
      Size and strain per reflection, as FullProf writes them in its .mic file: for each TCH phase and each
      reflection between Thmin and Thmax (λ1 and λ2), s = 1/d, the observed and instrumental Gaussian and Lorentzian
      FWHM (HGo HGi HLo HLi), the sample's integral breadths (betaG betaL, Å⁻¹ × 1000), the apparent size (Å; 99999 =
      no size broadening) and maximum strain (× 10⁴), h k l, 2θ (Bragg angle), FWHM and η of the whole profile, and
      Beta_TCH. Isotropic size (Y, IG) and strain (U, X), and Stephens' quartic form of strain (Str = 1, Strain-Model
      2, −2, 3) with its Lorentzian share; the resolution file gives the instrument's widths. Magnetic phases are
      left out, and said. β* takes each reflection's own λ and a negative U gives no strain; --fullprof takes λ1 and
      √|U| as FullProf.
  drs refine <file.pcr|file.drc> [--cycles <N>] [--data <file>] [--session <file>] [--plan | --step <k>] [--parameters]
             [--trace <file.json>] [--map <a> <b> [--span <s>esd | <value>] [--grid <n>]]
      Least-squares refinement of a FullProf control file, imported into the default model as by `drs calc`. The
      refined parameters are those of the .pcr's codewords, named by the parameter dictionary (phase1.cell_a,
      phase1.O1.atom_biso, …); codewords shared by several values are equality constraints with their multipliers, a
      lone multiplier is the damping of its parameter's shifts and line 9's R_at, R_an, R_pr, R_gl the damping of
      each class (W12); the crystal system and the atoms' sites constrain the rest. Levenberg–Marquardt with the
      normal matrix solved by eigen-decomposition, analytic derivatives, limits held as bounds (or reflected, or
      periodic), each value's shift per cycle at most the dictionary's maximum shift; it stops when every computed
      shift is below the .pcr's Eps times its esd, or after N cycles (default: the .pcr's NCY), or when it diverges
      (χ² rising 3 cycles in a row, a shift of 100 esd's from cycle 3, an impossible value) or the matrix is singular,
      naming the combination the data cannot fix. Each cycle's notes say where a step was shortened, held, reflected
      or wrapped; --parameters prints each parameter's value, esd, shift and shift/esd at every cycle.
      --session <file>: the run as a session a script can resume (W11): with --cycles N it runs N cycles and pauses
      (Step N); the state (values, λ, cycle count, history) is saved in the file and the next call resumes from it.
      --plan runs the control file's plan ([[refine.step]]: what each step frees, its damping, cycles and
      convergence) step after step, each from the values the previous reached; --step <k> runs step k alone.
      Profile-matching phases (Jbt = 2) have their intensities extracted at every cycle (Le Bail), and then all the
      cycles are run.
      Prints Rp, Rwp, Rexp, χ² and the Durbin–Watson d per cycle, then each parameter with its esd, R_Bragg and R_F
      per phase, the restraints (distances, angles and linear combinations: target, value, Δ/σ; their weight is χ²_ν, as
      FullProf's), the weight fractions (ATZ > 0 in two phases or more; Hill & Howard, manual p. 31, esd's from the
      covariance), the amorphous content when a phase is an internal standard (internal_standard_weight_fraction in
      the .drc; Orlhac et al. 2001, eqs. 3–4), and the correlations above 0.8. esd = √(C_kk χ²), χ² over N − P (FullProf scales by the χ² over the
      points with Bragg contributions). Ctrl-C stops it after the cycle in progress, with the results at that cycle.
      --trace <file.json>: the run's record (E1): per cycle χ², R's, λ, the refused steps, each parameter's value,
      shift, esd and computed shift, the spectrum of the scaled normal matrix, the covariance and the time; the
      steering's notes as data, the stop's facts with its cycle and the pauses (E15).
      --map <a> <b>: χ² on a grid of the two refined parameters about the values reached (E2), each over ± its span
      (--span 3esd, the default, or a value in its units), --grid n points a side (21), the points computed coarse
      to fine on every thread; printed as a b χ².
  drs data formats
      The reader registry: every data format of FullProf's list (with its Ins, the .pcr translation), the general
      table, and the formats GSAS-II reads; for each, its id, name, usual extensions and where its layout comes from.
  drs data identify <file> [--prefer <id>,<id>,…]
      The formats a data file may be in, by content first, then extension: each reader's confidence (1 signature,
      0.8 structure, 0.5 loose columns; +0.1 for a usual extension). The file is identified when the best leads the
      next by 0.1, otherwise the candidates are listed; --prefer: among tied readers the first named reads it (the
      app's request 51), and the line says the preference chose it.
  drs data read <file> [--as <id>] [--scan <n>] [--out <file>] [--option <key>=<value> …] [--prefer <id>,…]
                [--delimiter <c>] [--decimal <c>] [--header <lines>] [--x <col>] [--y <col>] [--sigma <col>|none]
                [--axis 2theta|d|q|tof]
      Read a data file, as the format it is identified as or the one given with --as, and write each scan (or scan
      n) as x, y, σ after a header with what the file says: format, quantity on x, wavelength, title, instrument,
      date, sample, temperature, monitor, counting time, step. σ = inf marks a point with no weight. With --out
      <file.xye|file.xys>, the scan is written in that format instead, to be read back (σ 1e30 for no weight); data op
      and data apply do the same. The table
      options (columns from 1) choose what the general table reader would otherwise recognise, and imply --as table.
      --option: a reader's option as drs data formats lists it (scan=2, wavelength=1.5406, axis=tof, tof_unit=ms,
      x_column=1 …; the app's request 51); --prefer: the order of preference among tied readers.
  drs convert <file.pcr> [--to <file.drc>]
  drs convert <file.drc> [--to <file.pcr>]
      A FullProf control file imported into a control file of diff-ref-struct (the model as `drs calc` imports it,
      the refinement as `drs refine` reads it, the data file by its path relative to the new file and its reader), or
      a control file exported to a .pcr where FullProf can express it: values translated back to FullProf's rules,
      codewords from the refinement, the resolution in a .irf and the reflections of profile-matching phases in
      CODFILn.hkl next to it, processed data in <stem>.xys (Ins = 10). Prints what the .pcr cannot hold.
  drs rigid <file.pcr>
      The rigid groups of a rigid-body phase (Jbt = 4, P41): each group's kind (spherical, Cartesian, cylindrical,
      Z-matrix), origin, θ, φ, χ and satellite anchor, and each atom placed from them, against the x, y, z the file
      writes beside them (FullProf's placement), with the difference in Å.
  drs texture <file.pcr|file.drc>
      Each phase's preferred orientation as axes (P30): the function, each axis's direction, March's ratio r and
      fraction (FullProf's Nor = 2, 3; one March–Dollase axis with its untextured fraction G2 as a second axis of
      r = 1), the anvil cell's integration steps, and March–Dollase's density D(r, α) of each axis (multiples of a
      random distribution) at α = 0°, 15°, …, 90° from it.
  drs move <file.drc> <folder>
      A control file moved to another folder (the app's File ▸ Move To…, request 45): every path it holds (the data,
      the background points, the data operations' files, a mode phase's sign reference; a joint file's patterns'
      control files) rewritten to reach the same files from there, its comments, order and formatting kept; the file
      written in the folder, and the old one removed. Refused, the file named, when a path names a file that is not
      there, or the folder holds a file of the same name.
  drs check <file.drc|file.pcr> [<data>] [--json] [--correlation <c>] [--no-run]
      The core's checks (R24): the file (syntax, unknown names with the nearest known one, units, values outside the
      engine's limits), the model (data and resolution files, wavelength, range against the data, occupancies beyond a
      site, values fixed by a site's symmetry set free, phases without reflections), the refinement (holes and counts
      of a .pcr's parameter numbers, ties that do not build, parameters that do not change the pattern) and, before a
      run, the normal matrix at the starting values: its null and near-null directions named as the parameters they
      combine, and the correlations at or above c (default 0.9), explained from the parameter dictionary. Each as
      file:line:column (or the model's addresses): severity[code]: what, why, fix. --json writes them as JSON, as the
      app reads them. Exit status 1 when there are errors. A .pcr's data are the file given after it, as FullProf
      takes them on its command line (fp2k CODFIL DATFIL), else CODFIL.dat next to it; a control file names its own.
  drs data info <file> [--in <project folder>]
      What a file of a project is (W5, L081): a data file (its format, scans, points, range, wavelength), or one of
      some 70 kinds of FullProf's and its suite's files (control files, listings, profiles, reflection lists, Fourier
      files, .fst, .mcif, Bond_Str's, ...), with its role and the source that says what it is; unknown otherwise.
      With --in, a file of a project the app started is named first by its place and name (the file dictionary):
      a results file in Results/ is a result, not data.
  drs data op <file> <operation> [key=value ...] [--as <id>] [--out <file>]
      One data operation applied to the file's first scan, the file left unchanged; the keys are those of the control
      file's [[pattern.operation]]: mask from= to= | mask points=[x, ...]; delete (the same); smooth points= order=
      (Savitzky–Golay); rebin factor= | step= [mode=sum|average]; normalize by=monitor|time|maximum to= | by=factor
      factor=; add | subtract | merge file= [reader=] [scale=]; convert to=2theta|d|q [wavelength=];
      wavelength to=. Writes x, y, σ (σ propagated; inf = no weight); with --out <file.xye|file.xys>, that format
      (formats::write, app request 12). `drs data op <operation> --help`: its settings in the dictionary's words.
  drs data settings [<operation>]
      What each data operation and analysis does and what each of its settings means: unit, range, default, how a
      value is proposed, what raising and lowering it does, the definition it follows (app request 16; the manual page
      \"Data operations: the settings\" shows the same words); --help after op, peaks or background shows theirs.
  drs data apply <file.drc> [--out <file>]
      The data of a control file after its operations, in order.
  drs data peaks <file> [--prominence <counts>] [--significance <n>] [--as <id>] [--out <file>]
      The peaks: local maxima with SciPy's prominence ≥ the given counts (default 0) and ≥ n σ (default 3), their
      height and width at half the prominence.
  drs data background <file> [--log-lambda <l>] [--curve] [--as <id>] [--out <file>]
      The automatic background of GSAS-II: arPLS (Baek et al. 2015) with λ = 10^l (default GSAS-II's for the number of
      points), at most 10 iterations; points every N/100 (or the whole curve with --curve).
  drs background read <file.drb|file.bac|file.drc>
      Background points (R26): our file (.drb: format = \"diff-ref-struct background 1\", interpolation, points), FullProf's
      .bac in either of its formats (manual p. 8), or a control file's own points; printed as 2θ, intensity.
  drs background convert <in> --to <out.drb|out.bac> [--data <file>] [--spline]
      The points written in another format (the extension's). A .bac is FullProf's first format (start step end,
      then the background at every point), the only one FullProf reads: at the data's points (--data), or at the
      points themselves when they are equally spaced.
  drs background eval <points> --data <file> [--spline] [--out <file>]
      The background at the data's points: linear interpolation (FullProf's Nba ≥ 2) or a natural cubic spline
      (FullProf's cubic splines, Nba < −5; a .pcr's points are read with the point FullProf adds at the data's
      last point, L081).
  drs background fit <points> (--chebyshev <n> [--from <2θ> --to <2θ>] | --polynomial <n> --origin <2θ>) [--into <file.drc>]
      A function fitted to the points by least squares: the coefficients and the rms of the residuals; with --into,
      the control file's background becomes that function, as starting values (its comments kept).
  drs background store <file.drc> (--file <name.drb|name.bac> | --inline)
      Where a control file keeps its background points: in a file it names (written next to it), or in itself. A .bac
      is written at the control file's data points.
  drs new <data> (--cif <file.cif> | --phase <formula> <space group> <cell values>)… [--out <file.drc>]
          [--as <reader>] [--wavelength <λ1>[,<λ2>,<ratio>] | --anode <Cu|Mo|…>] [--neutron]
          [--polarisation <Cthm>] [--window <Wdt>] [--asymmetry-limit <2θ>] [--profile <shape>]
      A control file started from a data file and what is known of each phase (L153): a CIF (its structure), or
      a formula, a space group and the cell's free values (\"Y2O3\" \"I a -3\" 10.6; no atoms: Le Bail). The
      scale of each phase with atoms and every point of the automatic background (arPLS, kept in
      \"<stem> background.drb\") are free, the rest fixed; the widths start from the data's peaks. Every value not
      given and not read from the data is a default, listed at the head of the file with its source.
      --out: default <data stem>.drc beside the data. Nothing is overwritten.
  drs compare <a> <b> [--reference a|b] [--pair <address> <address>] [--json]
      Two refinements' paths compared (Explore's E11): each of <a>, <b> a results file (its control file, cycles,
      su's, correlations and stop inside) or a control file refined first. Prints the parameters in common and apart,
      the fixed values changed between the two starts, both runs' last χ² and stop, and the ends' differences sorted
      by |difference|/√(σA² + σB²) (agreeing within 3). --pair: each run's path, end, σ's and ellipse for two
      parameters. --json writes lsq::PathsCompared (an object of paths and pair with --pair).
  drs compare-with fullprof prepare <file.drc> <folder> [--cycles <n>]
  drs compare-with fullprof read <file.drc> <folder> [--json]
  drs compare-with gsas2 prepare <file.drc> <folder> [--cycles <n>]
  drs compare-with gsas2 read <file.drc> <folder> [--json]
      A refinement checked against FullProf (the app's request 49; compare::fullprof_prepare, fullprof_read).
      prepare: the run's folder (compare.pcr and the files it names, the data, run.json) and the command, printed with
      the manifest of what was carried over, approximated or refused. Run it there: FULLPROF=<fp2k's folder> fp2k
      compare <data>, nothing on its input. read: FullProf's outputs read, the control file refined (a copy) with the
      same free parameters and cycles, and the two compared: χ², the R factors, each parameter with both esd's and
      Δ/σ, the parameters not compared; written into the folder as comparison.json (--json: printed too).
      gsas2: the same against GSAS-II (request 50): compare.py, compare_input.json, the data as compare.xye, the
      instrument file and each phase's CIF; run there with GSAS-II's Python, PYTHONPATH=<the folder holding GSASII>.
  drs data constant-time <file.uxd|file.xrdml> --time <s> [--out <file.xys>]
      The scans of a multi-scan file brought to one counting time per point and joined (PL6, WinPLOTR's constant time
      normalisation): counts × time / the scan's counting time, σ alike; at an x an earlier scan has, the earlier kept.
  drs plots <command> … --x-interval a:b
      Every drs plots command that opens a set takes the x interval of the analysis (PL42): every member cut to a ≤ x
      ≤ b on opening (in the .drv's x quantity, else the files' own), the files untouched, as a set file's x_interval;
      a member that holds no point inside is refused by name; save writes it, propose gives it to its proposals.
  drs plots save <file.buf | file.drv | data files…> [--out <Plots/set.drv>] [--name <set>]
      The set as a plot set file (PL12): Plots/<set>.drv by default (the file dictionary), the control file's TOML
      subset, the files relative to the project's folder. Every drs plots command opens a .drv as it opens a .buf:
      the files relative to the folder above Plots/ (or --base), the typed values with their sources, the order.
  drs plots export <set> --csv [--long] [--x 2theta|d|q] [--wavelength λ] [--range …] [--step …] [--out <file.csv>]
      The set as CSV (PL7, PL12): wide, x then y and σ of each member on the common grid, or --long, one row per member
      and point on its own points; each member's file, scan, format and value of the variable with its source in #
      lines. In d or Q each member is converted with its own λ, or --wavelength when its file states none. With
      --by and the options of drs plots normalise, the members as the view shows them, normalised in the files' x
      and then converted, the rule and each member's factor written in # lines (PL25); likewise drs plots derive.
  drs plots normalise <set> --by monitor|time|maximum|area|peak|reference [--to v] [--within from:to]
          [--member k] [--background auto|line] [--json]
      Each member's factor for a normalisation as a view (PL6): no file and no operation changed. maximum (in a range
      with --within), the same area over a range, a peak's integral above its background, or member k's area over a
      range (--by reference).
  drs plots view <set> --x 2theta|d|q|stl|s [--lambda k=λ]… [--lambda-from <file.drc> | --wavelength λ] [--json]
      The set in another x quantity as a view (PL10): each member with its own λ (--lambda for member k, else its
      file's, else the control file's λ1 or the one typed); the members without one listed.
  Every drs plots command takes --temperatures <file>: FullProf's temperatures file (lines like data_23.dat TSAMP 1.5; PL21),
  each member its file names taking that value as its source; a name given twice with different values takes none.
  drs plots step <set> [--scale linear|sqrt|log] [--json]
      A step proposed for a stack (PL19): the median over the members of their 90th less 10th percentile on the scale,
      rounded up to 1, 2 or 5 × 10^k, with its rule.
  drs plots open <file.buf | data files…> [--as <reader>] [--set <k>=<value>:<note>]… [--json]
      Several data files as one set (R32, Plots; PL1–PL2): every scan a member with the reader's facts (format,
      axis, points, range, step, λ, title, temperature, monitor, time) and the set's variable, the temperature (K),
      with where it was found: the header, a WinPLOTR buffer's comment, the title (\"T.=.4.K\"), the file name
      (\"_300K\", \"_T300\"); none when nothing says it, never a guess. Members sorted by it, those without one last.
      A .buf names its files relative to it. --set gives member k (from 1, in the order shown) a typed value.
  drs plots displace <file.buf | data files…> (--step <dy> [--dx <dx>] | --per-unit <dy> [--dx <dx>]) [--scale linear|sqrt|log] [--json]
      Each member displaced as a stack or a projection (PL3): k·dx and k·dy for the k-th member, or in proportion to
      the variable from the first member's; on the y scale shown (in √ the shift is added to √y). Prints each
      member's shift and baseline; --json the displaced curves.
  drs plots grid <file.buf | data files…> [--range overlap|union] [--step finest|coarsest|<step>] [--json]
  drs plots map <file.buf | data files…> [--range …] [--step …] [--scale linear|sqrt|log] [--reference <k>] [--json]
      The members on one grid (PL4): their own points when they share them, else interpolated where the grid is finer
      and averaged over each bin where it is coarser, σ propagated. The map: one row per member at its value of the
      variable (edges at the midpoints), those without one apart, the colour range; --reference k (from 1) gives the
      difference map, each member less member k, on a range symmetric about zero.
  drs plots derive <file.buf | data files…> (--sum | --average | --weighted-average | --difference) --members <k,…>
          [--range …] [--step …] [--out <file.xys>] [--name <set>]
      The sum, average or difference of members (from 1, in the order shown) with σ (PL9), written as an XYDATA file
      (.xys, negative y kept) that names its members and their values of the variable; by default
      Plots/<set> <operation>.xys (the file dictionary).
  drs plots windows <file.buf | data files…> --window <name:from:to>… [--background auto|line|<points file>]
          [--log-lambda v] [--at x,…] [--json]
      For each window and member (PL5): the integral of y above the background with σ, the centroid with σ, the
      strongest peak of the peak search; the background at the x given by --at. auto: each member's automatic
      background; line: the line between the member's values at the window's ends; a .drb or .bac: its points.
  drs plots project <file.buf | data files…> (--step <dy> | --per-unit <dy>) [--dx <dx>] [--scale …]
          [--order <k,…>] [--floor <y>] [--json]
      The projection with hidden parts (PL17): the members displaced, in --order front to back (from 1; the set's
      order by default), each curve hiding what lies under it down to --floor plus its shift (without end if no
      floor); the visible pieces of each curve.
  drs plots describe <set.drv> [--json]
      What a plot set file says beyond its members (PL12, PL27, PL37): the variable, x and scales, the map's and the
      difference map's colour range as clipped (map_range, difference_range), the windows and the regions saved by
      name; a region the set cannot be shown in (no λ to convert it) is said as refused.
  drs plots region <set.drv> --name <region> [--scale linear|sqrt|log] [--reference <k>] [--json]
      A region's block of the map (PL37): the rows and columns inside it and the smallest and largest value of its
      cells, the map as drs plots map draws it, on the set's map_scale (else its y_scale), the difference map with
      --reference on its own scale (the square root for a log map); the region's own colour range for that map, and
      the other map's named (P65); the region's bounds
      converted to the map's x with each member's λ (refused, the members named, when the conversion needs one a
      file does not state).
  drs plots show <set.drv> [--view map|difference] [--browser] [--out <page.html>] [--no-open] [--no-patterns]
          [--reference <k>] [--wavelength λ]
      The set's map without the app's window (PL38, I08): when the app is installed (macOS) and --browser is not
      given, open -a diff-ref-struct <set.drv> --args --view <v> (the app takes --view from PL39). Otherwise one
      self-contained page beside the set, <set>.map.html (or --out): the app's Map viewer and the page's data (the map,
      the difference map with the file's reference, the clip, the regions with their blocks, the members' patterns
      unless --no-patterns) inline, nothing fetched; the set, the file, the core's version and the date on it; opened
      with open (xdg-open on Linux) unless --no-open. Never over a file that is not such a page.
  drs plots propose <data files… | file.buf> [--json]
      The series a run's files offer (W2): by temperature (one scan per temperature), the holds (three or more scans
      at one temperature, against the time) and the whole run against the time, per exposure class (the esd's
      median, a class where it more than doubles); the time from the elapsed seconds every name gives (_3600s),
      else from the files' modification times; each a plot set file (--json gives them), the files left out with why.
  drs plots seq <file.seq> [--parameter <name>] [--json]
      FullProf's .seq as a table (PL11): per block the data file, PARVAL's index, external parameter and numor, χ²,
      each phase's R_Bragg, volume and weight fraction with σ, and every refined parameter with its value and σ.
  drs rebase <file.drc|file.pcr> <data> --to <new.drc> [--as <reader>]
      A control file pointed at new data (R25, W3): written to new.drc with the data file relative to it, and what
      does not fit the new data reported as the checks report (drs check): another wavelength or axis, range,
      instrument or step, excluded regions and background points outside the data, a Chebyshev range of other
      angles; each with a suggestion. The old .pcr's range and regions are dropped from [pattern.fullprof].
  drs files [<path in a project>]
      The file dictionary (R27): every kind of file of a project started in the app, its folder and name pattern,
      role and status (decided or proposed); with a path relative to the project's folder, the kind it is. Without
      --to, drs convert writes the dictionary's name: <stem>.drc beside a .pcr, Exports/FullProf/<stem>.pcr from a .drc (L123).
  drs modes <file.pcr|file.drc> [--amplimodes <output.txt> [--analysis k]] [--project <file.pcr> [--shift x,y,z]]
            [--equivalent] [--reference <file.pcr|file.drc>]
  drs modes <amplimodes-output.txt>
      Refinement by symmetry modes (R29): a phase's modes (FullProf's Jbt = 6, or [[phase.mode]]), each with its
      irrep and amplitude, the irreps' amplitudes √(ΣQ²) and polarisation vectors, the positions they give. With the
      reference metric (an AMPLIMODES analysis's parent and transformation, or mode_reference_cell), each mode's
      Σ|cV|² over the primitive cell (1 Å² for AMPLIMODES's modes), and with --project another structure's
      coordinates (a coordinate refinement) decomposed into these modes. --equivalent lists the amplitude sets
      that describe the same structure up to the parent's symmetry (E7), and --reference the one nearest another
      file's amplitudes (by mode name), its signs as the convention. Given AMPLIMODES's text output itself, its
      analyses: parent, distorted structure, reference metric, summary of amplitudes.
  drs decompose <parent.cif> <distorted.cif> --basis <a',b',c'> [--origin x,y,z] [--fullprof]
  drs decompose --amplimodes <output.txt> [--analysis k] [--mapped] [--fullprof]
      The symmetry modes of a distorted structure computed from its parent (P13, L263), as AMPLIMODES computes them:
      the reference structure, the displacements, and for each irrep of the parent's group that the distorted group
      keeps, its k-vector star, parity, dimension, number of modes and amplitude (Å, normalised in the distorted
      phase's primitive cell), orbit by orbit. --basis gives the distorted phase's basis in the parent's (ITA:
      a-b,a+b,c or c,b,-a-c), --origin its origin. Irreps are named from k, parity and dimension (GM+[3]a); the
      number of a label such as GM4+ needs tables, and is taken from AMPLIMODES's modes when given. --amplimodes takes
      the structures and transformation of an AMPLIMODES analysis and prints its summary beside (--mapped: its atom
      mappings' positions, which AMPLIMODES may have adjusted, instead of the input). --fullprof prints the block for
      the .pcr (atoms, V_MODES, A_MODES).
  drs isodistort <file.cif>
      The modes of an ISODISTORT CIF (P13): each mode's label split (parent, k, irrep, direction, atom, Wyckoff,
      kind), amplitude and normalisation; each irrep's amplitude (the root of its modes' squares); each coordinate
      delta as the file gives it and as the modes give it through the mode matrix, δ = M·(A × N). Checked on
      ISODISTORT 6.12.2's CIF of case A against the core's own decomposition (L278).
  drs state <file.pcr|file.drc> [--cycles <n>] [--cycle <k>] [--step <°>] [--dmax <Å>] [--at <name>=<value>,…]
      A state's derived quantities in one bundle (E3): the agreement, each phase's cell and volume with su's (E18a),
      FWHM, η, the unclamped Gaussian H_G² (marked where negative) and H_L every --step degrees (10), the contacts up to
      --dmax (3 Å) by atom, operator and translation with su's and a cation–anion contact's plausible distances (E18b),
      the bond-valence sums when the atoms have charges, each against |q| and its plausible range (E18c), the irreps'
      amplitudes with su's, the restraints. At the file's values, or after --cycles n of refinement at its cycle k.
      --at: the agreement with those free parameters at those values, as a step sets them (E20a).
  drs structure <file.pcr|file.drc> --cartesian [--phase <p>] [--box x,y,z:x,y,z] [--within <Å>] [--cycles <n>] [--cycle <k>]
          [--modes <file.pcr|file.drc> [--amplimodes <output.txt> [--analysis k]]]
      A phase's scene for a 3D view (E4): every image in the box (the cell by default) in fractions and Cartesian Å,
      named by atom, operator and translation; the bonds (cation–anion within 1.2 times the ionic radii, or every pair
      within --within Å); the coordination polyhedra (faces of the hull, volume, Baur's distortion); the cell's edges;
      for a mode phase, each image's reference position (fractions and Cartesian) and each irrep's displacement on
      it. --modes <file> describes a phase refined by its coordinates by that file's mode phase in the same setting
      (P19; its reference metric from --amplimodes <output> [--analysis k] when given), with the residual.
  drs bonds <file.cfl|file.pcr|file.drc> [--dmax <Å>] [--angles <Å>] [--tol <%>] [--fullprof]
      Distances, angles and bond-valence sums of each structure phase (P12), as FullProf's Bond_Str lists them in
      CODFILn.dis: the neighbours of each atom within dmax (default the CFL's DMAX, else 3.0 Å), with their standard
      deviations when the file gives them (a CFL); the angles between neighbours within --angles; each cation–anion
      bond within (r_cation + r_anion)(1 + tol) (default 20 %), its valence exp((d0 − d)/B0) with CrysFML's table, the
      sums, the distortion and the global instability indices. --fullprof: Bond_Str's occupancies (normalised to the
      first atom) and standard deviations (without the cell's).
  drs cif <file.pcr|file.drc> [--cycles <n>] [--out <file.cif>]
      A CIF for deposition (P12): each structure phase with its formula and Z, symmetry, cell, experiment, agreement
      factors, atoms, bonds and angles, with the standard uncertainties of the engine's covariance after n cycles
      (default 0: at the file's values); \"?\" where the authors must say.
  drs notes <control>
      The user's notes on the reflections of a control file (W8): the notes file beside it, <control> notes.toml, one
      line per note (phase and h k l, or 2θ, and the text). drs results carries them into the results file.
  drs trend <address> <control[:results]>… [--base <folder>]
      A parameter across versions of a control file (W7): for each control file (oldest first) the value at the
      address (phase1.cell_a, …) and, given the version's results file after a colon, its su. --base: the folder
      the versions' relative files are read from (default: each file's own).
  drs joint <file.pcr|joint.drc> [--cycles <n>] [--esd default|fullprof] [--fullprof] [--results] [--project <folder>]
            [--time \"…\"] [--save <joint.drc>]
  drs joint <joint results.cif> [--pattern <k>]
      Each pattern's control file from a joint refinement's results file, its model, phases and free parameters as
      the refinement left them, read without the data (Publish's P50); from a file written before core 0.61.0, from
      its source .pcr or joint control file, with the files they name in the results file's folder.
      A joint refinement (P06) of a multi-pattern .pcr (NPATT, W_PAT): each pattern its own model and data, the
      structures shared, the parameters tied by FullProf's codeword numbers; Q = Σ W_k Σw(y − yc)² minimised (FullProf's
      objective). Prints each cycle, each pattern's and the joint agreement (global χ²: Σ W_k χ²_k over the points
      with Bragg contributions, FullProf's) and the parameters; --esd fullprof gives FullProf's esd's, √(C Σ W_k χ²_k),
      the default √(C Q/(N − P)); --fullprof keeps the values under FullProf's rules. --results (or --project) writes
      Results/<control> <time> results.cif with a block per pattern (P41). A joint control file (format = \"diff-ref-struct
      joint 1\": [[joint.pattern]] each pattern's control file, weight, refined flag and phases; [[joint.tie]] the
      parameters they share) is read as a multi-pattern .pcr; --save writes it, beside <joint> pattern <k>.drc, at the
      values reached (app request 43). drs refine and drs results on a multi-pattern .pcr or a joint control file run it.
  drs results <file.drc|file.pcr> [--cycles <n> | --plan] [--project <folder>] [--time \"YYYY-MM-DD HH.MM.SS\"] [--exports]
      The results of a refinement of n cycles (default 0: at the file's values) as the file dictionary names them
      (R28): Results/<control> <time> results.cif, a powder CIF with every phase, the profile (2θ, Yobs, weight,
      background, Ycal, each phase's curve), the reflections (with I obs partitioned by the model), the agreement
      factors, the refined parameters with their standard uncertainties and the control file; and one structure CIF
      per phase, Results/<control> <phase>.cif; the summary and the listing to read, Results/<control> <time>
      summary.txt and listing.txt (P57: FullProf's .sum and .out in the core's words). --exports adds Exports/<control> reflections.csv and profile.csv, and
      for each Le Bail phase its extracted intensities, Exports/<control> <phase> intensities.int (FullProf's
      CODFILn_cltr.int: overlapped reflections grouped by FullProf's rule, Rmub 1.0, Rmuc 0.2; P35). The
      time is UTC unless given; the folder is the file's unless --project names another. --plan runs the control
      file's plan of steps instead ([[refine.step]], as drs refine --plan), each step's cycles and convergence its
      own; the results file's cycles are numbered on over the steps, each with its step's number and name
      (_drs_cycle_step, _drs_cycle_step_name; P37).
  drs plan <file.pcr|file.drc> [--to <file.drc>]
      Guided refinement (P55): the documented order of refinement (FullProf's manual, p. 14; McCusker et al. 1999,
      §7) for the file's model, as a plan of steps, each with its source and the values it adds and drops: scale;
      zero, first background term, cell; the heavier atoms' positions, then all and the overall B; peak shape and
      asymmetry; occupancies of partly occupied sites; individual B's; the background; betas; Sycos and Sysin;
      microstructure; all together to a largest shift/esd of 0.10. --to writes the control file with the plan
      ([[refine.step]]), every value it frees added to the refined parameters; drs refine <file.drc> --plan runs it.
  drs pawley <file.pcr|file.drc> --to <file.drc> [--phase <n>] [--iterations <k>]
      A Pawley fit (P36; Pawley 1981): phase n's (default 1) reflection intensities refined by least squares with
      the cell and profile, in place of a structure or of Le Bail's partition. A phase with a structure takes the
      list of its reflections to 5° past the data; the intensities start from k Le Bail partitions (default 10) at
      the file's cell and profile; reflections that coincide for every cell the crystal system allows are tied with
      equal |F|² (as GSAS-II's Pawley fit ties them); the scale factor is no longer refined, and combinations of
      accidentally overlapped reflections the data cannot fix are left out of each cycle (the guard singular off).
      --to writes the control file: its intensities named phase<n>.reflection_intensity.<h>_<k>_<l>, the tied ones
      one parameter. drs refine runs it; drs results writes each intensity with its su from the covariance.
  drs sequential <file.drc|file.pcr> <data file>… [--each] [--leave-out-unmodelled] [--project <folder>] [--time \"…\"]
      A sequential refinement (P54): the data files refined in turn with the file's model and refinement (its plan,
      when it has one), each read with the file's reader and data operations and labelled by its name; each starts
      from the values the one before ended at, as FullProf's sequential mode, or with --each from the file's values.
      The external parameter is the temperature the data file's header gives. Writes, as the file dictionary names
      them (P59), Results/<control> <time> series.csv (per pattern: its label, parameter, cycles, stop, χ², R's,
      each phase's cell volume and su, every free parameter's value and su) and series.cif (an index block, then
      every pattern's results file). --leave-out-unmodelled (a .pcr) leaves out the phases the core does not model
      (magnetic, |Jbt| = 1) and says which.
  drs parametric <series.csv> <column> --equation <equation> [--start <name>=<value>,…] [--x <column>] [--variable <name>]
      A parametric fit (P54): a column of a series table (a refined value or a cell volume, its su from the column
      \"<column> su\") against the external parameter (--x, default parameter), fitted to an equation of the variable
      (default T: numbers, + − × / ^, parentheses, exp ln log10 sqrt sin cos tan abs sinh cosh tanh, pi; every other
      name a coefficient, starting at 1 unless --start gives it) by weighted least squares. Prints each coefficient
      with its su (√(C χ²_ν), as SciPy's curve_fit with absolute_sigma=False), χ² and χ²_ν.
  drs locate <file.drc> [<address> | --at <byte>]
      Where the control file writes the value of an address (phase1.cell_a, phase1.Y:2.atom_position.x,
      pattern1.zero_shift; the app's requests 40 and 41): its line and bytes, and the text there; nothing when the file
      does not write it (a value at its default). The file is read with its folder (background points in a file of
      their own). Without an address, every value it writes, from one reading; --at, the address of the value written
      at that byte.
  drs impact <file.drc|file.pcr> [--cycles <n>] [--parameter <address>] [--top <k>] [--json]
      Each parameter's impact (P45): what refining it alone would do to first order, refined or not, the largest
      fall of chi2 first. Its slope d(chi2)/dp, the shift and that shift in its own esd's, the fall of chi2, the
      largest change of the calculated pattern (counts, at 2theta), the change for one esd and the 2theta ranges it
      acts on. --cycles: refine n cycles first (default 0); --parameter: one parameter, with the change point by
      point; --top: the first k (default all).
  drs normal <file.drc|file.pcr> [--cycles <n>] [--points] [--json]
      The normal probability plot of the weighted residuals (P45; Abrahams & Keve 1971, GSAS-II's Error Analysis):
      their slope and intercept against the expected normal order statistics over the central 68 %, the number of
      points; --points lists each pair (expected, residual).
  drs peak <data file> <from> <to> [--gaussian | --lorentzian] [--doublet <λ1> <λ2> <ratio> | --doublet scan]
      Fits one peak between 2θ from and to (P51): a pseudo-Voigt (η refined, or held at 0 or 1) over a linear
      background, weighted least squares, each parameter with its esd (√(C χ²), as the engine's): the position, the
      integrated intensity, the FWHM, η and the background. --doublet adds the λ2 peak at the same width and shape
      with the given intensity ratio (Cu Kα: 1.54056 1.54439 0.5); --doublet scan takes the file's own λ1, λ2 and
      ratio (the app's request 38). The data file is read as `drs data` reads it.
  drs name <file.pcr|file.drc> <name> [--phase <k>] [--label <atom>]
      The parameter a name of an output names (the app's request 36), as an address of the file's model: FullProf's
      symbolic names of its .out and .sum (Scale_ph1_pat1, Cell_A_ph1_pat1, X_O1_ph1_pat1, A3_GM1_ph1 …), or a CIF
      tag of a phase's block (_cell_length_a; _atom_site_fract_x with --label; --phase k, default 1).
  drs complete <file.drc> <line> <column>
      The names valid at a place of a control file (P49, the editor's completion), line and column from 1: the
      tables at a [, the keys of the table the line is in (less those it holds) with the dictionary's name and
      unit, and at a value the closed choices, the species, the space-group symbols or the model's addresses.
  drs fullprof <file.prf|file.hkl>
      FullProf's outputs of a run, read (app request 17): a .prf's title, phases, wavelengths, Zero, SyCos, SySin,
      propagation vectors and excluded regions, then 2θ, Yobs, Ycal and the background (FullProf's Backg with the
      quarter of the least background it takes off added back, L082), then the reflections (phase, wavelength, h k l,
      2θ); a .hkl's blocks (pattern, phase, columns as FullProf names them, one line per reflection).
  drs --version

Every command that reads a control file takes --base <folder>: its data, operations' and background files are named
from that folder instead of the file's own (a version in Versions/ names them from the project's folder).
";

/// The folder `--base` names, when given.
static BASE: std::sync::OnceLock<std::path::PathBuf> = std::sync::OnceLock::new();

/// A control file opened with its files named from `--base` when given, else from its own folder.
fn open_control(path: &std::path::Path) -> Result<drs::control::Opened, String> {
    match BASE.get() {
        Some(b) => drs::control::open_in(path, b),
        None => drs::control::open(path),
    }
}

#[cfg(unix)]
extern "C" {
    fn signal(signum: i32, handler: usize) -> usize;
}

fn main() -> ExitCode {
    // Rust ignores SIGPIPE, so `drs sg --list | head` would panic on the closed pipe; restore the
    // Unix default (SIGPIPE = 13, SIG_DFL = 0 on macOS and Linux) so the program just stops.
    #[cfg(unix)]
    unsafe {
        signal(13, 0);
    }
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    // --base <folder>, for every command: a control file's data, operations' and background files are named from
    // that folder (a version in Versions/, named from the project's; the app's request 33).
    if let Some(i) = args.iter().position(|a| a == "--base") {
        if i + 1 < args.len() {
            let _ = BASE.set(std::path::PathBuf::from(args.remove(i + 1)));
            args.remove(i);
        }
    }
    let result = match args.first().map(String::as_str) {
        Some("cell") => cmd_cell(&args[1..]),
        Some("sg") => cmd_sg(&args[1..]),
        Some("site") => cmd_site(&args[1..]),
        Some("hkl") => cmd_hkl(&args[1..]),
        Some("sf") => cmd_sf(&args[1..]),
        Some("calc") => cmd_calc(&args[1..], false),
        Some("simulate") => cmd_calc(&args[1..], true),
        Some("mic") => cmd_mic(&args[1..]),
        Some("refine") => cmd_refine(&args[1..]),
        Some("data") => cmd_data(&args[1..]),
        Some("convert") => cmd_convert(&args[1..]),
        Some("move") => cmd_move(&args[1..]),
        Some("rigid") => cmd_rigid(&args[1..]),
        Some("texture") => cmd_texture(&args[1..]),
        Some("check") => cmd_check(&args[1..]),
        Some("background") => cmd_background(&args[1..]),
        Some("rebase") => cmd_rebase(&args[1..]),
        Some("new") => cmd_new(&args[1..]),
        Some("plots") => cmd_plots(&args[1..]),
        Some("compare") => cmd_compare(&args[1..]),
        Some("compare-with") => cmd_compare_with(&args[1..]),
        Some("bonds") => cmd_bonds(&args[1..]),
        Some("state") => cmd_state(&args[1..]),
        Some("structure") => cmd_structure(&args[1..]),
        Some("modes") => cmd_modes(&args[1..]),
        Some("decompose") => cmd_decompose(&args[1..]),
        Some("isodistort") => cmd_isodistort(&args[1..]),
        Some("files") => cmd_files(&args[1..]),
        Some("cif") => cmd_cif(&args[1..]),
        Some("fullprof") => cmd_fullprof(&args[1..]),
        Some("peak") => cmd_peak(&args[1..]),
        Some("impact") => cmd_impact(&args[1..]),
        Some("locate") => cmd_locate(&args[1..]),
        Some("plan") => cmd_plan(&args[1..]),
        Some("pawley") => cmd_pawley(&args[1..]),
        Some("sequential") => cmd_sequential(&args[1..]),
        Some("parametric") => cmd_parametric(&args[1..]),
        Some("normal") => cmd_normal(&args[1..]),
        Some("complete") => cmd_complete(&args[1..]),
        Some("name") => cmd_name(&args[1..]),
        Some("reflections") => cmd_reflections(&args[1..]),
        Some("results") => cmd_results(&args[1..]),
        Some("joint") => cmd_joint(&args[1..]),
        Some("trend") => cmd_trend(&args[1..]),
        Some("notes") => cmd_notes(&args[1..]),
        Some("--version") | Some("-V") => {
            println!("drs {}", env!("CARGO_PKG_VERSION"));
            Ok(())
        }
        Some("--help") | Some("-h") | None => {
            print!("{USAGE}");
            Ok(())
        }
        Some(other) => Err(format!("unknown command '{other}'\n\n{USAGE}")),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("drs: {e}");
            ExitCode::FAILURE
        }
    }
}

fn num(s: &str, what: &str) -> Result<f64, String> {
    s.parse::<f64>()
        .map_err(|_| format!("{what}: '{s}' is not a number"))
}

fn cmd_cell(args: &[String]) -> Result<(), String> {
    if wants_help(args) {
        print!("{USAGE}");
        return Ok(());
    }
    let mut pos = Vec::new();
    let (mut hkl, mut lambda) = (None, None);
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--hkl" => {
                let v = args
                    .get(i + 1..i + 4)
                    .ok_or("--hkl needs three indices h k l")?;
                hkl = Some([num(&v[0], "h")?, num(&v[1], "k")?, num(&v[2], "l")?]);
                i += 4;
            }
            "--lambda" => {
                lambda = Some(num(
                    args.get(i + 1).ok_or("--lambda needs a wavelength in Å")?,
                    "lambda",
                )?);
                i += 2;
            }
            s => {
                pos.push(num(s, "cell parameter")?);
                i += 1;
            }
        }
    }
    if pos.len() != 6 {
        return Err(format!(
            "cell needs six parameters a b c alpha beta gamma, got {}\n\n{USAGE}",
            pos.len()
        ));
    }
    let cell =
        Cell::new(pos[0], pos[1], pos[2], pos[3], pos[4], pos[5]).map_err(|e| e.to_string())?;
    let r = cell.reciprocal();

    println!("Direct cell parameters:");
    println!(
        "  a = {:10.4}   b = {:10.4}   c = {:10.4}",
        cell.a, cell.b, cell.c
    );
    println!(
        "  alpha = {:8.3}   beta = {:8.3}   gamma = {:8.3}",
        cell.alpha, cell.beta, cell.gamma
    );
    println!("  Direct cell volume = {:.4} Å³", cell.volume());
    println!("Reciprocal cell parameters:");
    println!(
        "  a* = {:10.6}   b* = {:10.6}   c* = {:10.6}",
        r.a, r.b, r.c
    );
    println!(
        "  alpha* = {:8.3}   beta* = {:8.3}   gamma* = {:8.3}",
        r.alpha, r.beta, r.gamma
    );
    println!("  Reciprocal cell volume = {:.8} Å⁻³", r.volume());
    pair(
        "Direct and reciprocal metric tensors",
        "G",
        &cell.metric(),
        4,
        "G*",
        &cell.reciprocal_metric(),
        6,
    );
    println!("Cartesian frame: x along a, y in the ab-plane, z = x × y");
    pair(
        "Crystal to orthonormal / inverse",
        "M",
        &cell.crystal_to_orthonormal(),
        4,
        "M⁻¹",
        &cell.orthonormal_to_crystal(),
        6,
    );
    let b = cell.busing_levy_b();
    let binv = drs::mat3::inverse(&b).expect("B of a valid cell is invertible");
    pair(
        "Busing–Levy B (h_cart = B·h) / inverse",
        "B",
        &b,
        6,
        "B⁻¹",
        &binv,
        4,
    );

    if let Some(h) = hkl {
        let d = cell.d_spacing(h).ok_or("(0 0 0) has no d-spacing")?;
        print!("hkl = ({} {} {}):  d = {:.5} Å", h[0], h[1], h[2], d);
        match lambda.map(|l| (l, cell.two_theta(h, l))) {
            Some((l, Some(tt))) => println!(",  2θ = {tt:.4}° for λ = {l} Å"),
            Some((l, None)) => println!(",  beyond the limiting sphere for λ = {l} Å"),
            None => println!(),
        }
    }
    Ok(())
}

fn pair(title: &str, ln: &str, l: &Mat3, ld: usize, rn: &str, r: &Mat3, rd: usize) {
    println!("{title}:");
    println!("  {ln:>36}  {rn:>36}");
    for i in 0..3 {
        let f = |m: &Mat3, d: usize| {
            m[i].iter()
                .map(|x| format!("{:12.*}", d, x + 0.0))
                .collect::<String>()
        };
        println!("  {}  {}", f(l, ld), f(r, rd));
    }
}

fn wants_help(args: &[String]) -> bool {
    args.iter().any(|a| a == "--help" || a == "-h")
}

fn cmd_sg(args: &[String]) -> Result<(), String> {
    if wants_help(args) {
        print!("{USAGE}");
        return Ok(());
    }
    if args.iter().any(|a| a == "--list") {
        println!("Hall  No. Choice Hermann–Mauguin        Hall symbol");
        for s in SETTINGS.iter() {
            println!(
                "{:>4} {:>4} {:<5} {:<22} {}",
                s.hall_number, s.number, s.choice, s.hm, s.hall
            );
        }
        return Ok(());
    }
    let all = args.iter().any(|a| a == "--all");
    let sym = args
        .iter()
        .filter(|a| *a != "--all")
        .cloned()
        .collect::<Vec<_>>()
        .join(" ");
    if sym.trim().is_empty() {
        return Err(format!("sg needs a symbol\n\n{USAGE}"));
    }
    let (g, others) = SpaceGroup::from_symbol(&sym)?;
    match &g.setting {
        Some(s) => {
            println!("Space group {}  ({}),  no. {}", s.hm, s.hm_full, s.number);
            println!(
                "  Hall symbol: {}   setting: {}   (table entry {} of 530)",
                s.hall,
                if s.choice.is_empty() {
                    "standard"
                } else {
                    s.choice
                },
                s.hall_number
            );
        }
        None => println!(
            "Space group from the Hall symbol {} (not one of the 530 table settings)",
            g.hall
        ),
    }
    println!(
        "  Crystal system: {}   Laue class: {}   Point group: {}",
        g.crystal_system, g.laue, g.point_group
    );
    println!(
        "  Lattice: {}   centrosymmetric: {}{}   general multiplicity: {}",
        g.lattice,
        if g.centrosymmetric { "yes" } else { "no" },
        if g.inversion_at_origin {
            " (−1 at the origin)"
        } else if g.centrosymmetric {
            " (−1 off the origin)"
        } else {
            ""
        },
        g.multiplicity()
    );
    let frac = |t: i32| -> String {
        if t == 0 {
            "0".into()
        } else {
            let g = gcd(t, drs::symop::DEN);
            format!("{}/{}", t / g, drs::symop::DEN / g)
        }
    };
    let centring: Vec<String> = g
        .centring
        .iter()
        .map(|c| format!("({},{},{})", frac(c[0]), frac(c[1]), frac(c[2])))
        .collect();
    println!("  Centring vectors: {}", centring.join("  "));
    let list = if all {
        g.ops.clone()
    } else {
        g.coset_representatives()
    };
    if all {
        println!("All {} operators:", list.len());
    } else {
        println!(
            "{} operators, one per rotation (add the centring vectors for all {}):",
            list.len(),
            g.multiplicity()
        );
    }
    // Each operator's symmetry element as International Tables lists it (P09): symbol, intrinsic part, location.
    for (i, op) in list.iter().enumerate() {
        let e = drs::symmetry_elements::element_in(op, &g.centring);
        let sense = match e.sense {
            1 => "+",
            -1 => "-",
            _ => "",
        };
        println!("  {:>3}  {:<28} {:<5} {:<14} {}", i + 1, op.to_xyz(), format!("{}{sense}", e.symbol), e.translation, e.location);
    }
    if !others.is_empty() {
        println!("Other settings with this symbol:");
        for s in others {
            println!(
                "  {} (choice {}): {}   e.g. drs sg \"{}:{}\"",
                s.hm, s.choice, s.hall, s.hm, s.choice
            );
        }
    }
    Ok(())
}

fn gcd(a: i32, b: i32) -> i32 {
    if b == 0 {
        a.abs()
    } else {
        gcd(b, a % b)
    }
}

/// A coordinate as a decimal or a fraction such as 1/4 or -2/3.
fn coordinate(s: &str) -> Result<f64, String> {
    match s.split_once('/') {
        Some((a, b)) => Ok(num(a, "coordinate")? / num(b, "coordinate")?),
        None => num(s, "coordinate"),
    }
}

fn cmd_site(args: &[String]) -> Result<(), String> {
    if wants_help(args) {
        print!("{USAGE}");
        return Ok(());
    }
    let mut tol = drs::site::TOL;
    let mut rest: Vec<String> = Vec::new();
    let mut i = 0;
    while i < args.len() {
        if args[i] == "--tol" {
            tol = num(args.get(i + 1).ok_or("--tol needs a value")?, "tolerance")?;
            i += 2;
        } else {
            rest.push(args[i].clone());
            i += 1;
        }
    }
    if rest.len() < 4 {
        return Err(format!("site needs a symbol and three coordinates\n\n{USAGE}"));
    }
    let k = rest.len() - 3;
    let x = [coordinate(&rest[k])?, coordinate(&rest[k + 1])?, coordinate(&rest[k + 2])?];
    let (g, _) = SpaceGroup::from_symbol(&rest[..k].join(" "))?;
    let s = Site::with_tolerance(&g, x, tol)?;
    match &g.setting {
        Some(t) => println!("Space group {} (no. {}, Hall {}), general multiplicity {}", t.hm, t.number, t.hall, g.multiplicity()),
        None => println!("Space group from the Hall symbol {}, general multiplicity {}", g.hall, g.multiplicity()),
    }
    let p = s.position;
    println!("Position ({:.5}, {:.5}, {:.5})", p[0], p[1], p[2]);
    match s.wyckoff {
        Some(w) => println!("  Wyckoff position {}{}   site symmetry {}   (letters and symbols: spglib's table)", w.multiplicity, w.letter, w.site),
        None => println!("  Wyckoff position: not in the table for this setting"),
    }
    println!("  Multiplicity {}   site point group {} (order {})", s.multiplicity, s.point_group, s.ops.len());
    println!("Site-symmetry operators:");
    for op in &s.ops {
        let sop = drs::symop::SymOp::new(op.r, op.t);
        let kind = sop.rot_type().map_or("?", |k| k.symbol());
        let mut xyz = String::new();
        for (i, row) in op.r.iter().enumerate() {
            let mut e = String::new();
            for (j, v) in ["x", "y", "z"].iter().enumerate() {
                match row[j] {
                    0 => {}
                    c => e.push_str(&format!("{}{}{v}", if c < 0 { "-" } else if e.is_empty() { "" } else { "+" }, if c.abs() > 1 { c.abs().to_string() } else { String::new() })),
                }
            }
            let t = drs::site::Rational::new(op.t[i] as i64, drs::symop::DEN as i64);
            if !t.is_zero() {
                e.push_str(&format!("{}{t}", if t.num > 0 { "+" } else { "" }));
            }
            xyz.push_str(&e);
            if i < 2 {
                xyz.push(',');
            }
        }
        println!("  {xyz:<28} {kind}");
    }
    println!("Equivalent positions ({}):", s.orbit.len());
    for (i, o) in s.orbit.iter().enumerate() {
        println!("  {:>3}  {:9.5} {:9.5} {:9.5}", i + 1, o[0], o[1], o[2]);
    }
    let c = &s.coordinates;
    println!("Coordinates: ({})   free: {}", (0..3).map(|i| c.expression(i)).collect::<Vec<_>>().join(", "),
        if c.free.is_empty() { "none".into() } else { c.free.iter().map(|&i| c.names[i]).collect::<Vec<_>>().join(" ") });
    let b = &s.adp;
    println!("Anisotropic displacements: {}   free: {}",
        (0..6).map(|i| format!("{}={}", b.names[i], b.expression(i))).collect::<Vec<_>>().join("  "),
        b.free.iter().map(|&i| b.names[i]).collect::<Vec<_>>().join(" "));
    Ok(())
}

fn cmd_hkl(args: &[String]) -> Result<(), String> {
    if wants_help(args) {
        print!("{USAGE}");
        return Ok(());
    }
    let (mut dmin, mut dmax, mut lambda, mut tt) = (None, None, None, None);
    let mut rest: Vec<String> = Vec::new();
    let mut i = 0;
    let val = |i: usize, what: &str| -> Result<f64, String> { num(args.get(i).ok_or(format!("{what} needs a value"))?, what) };
    while i < args.len() {
        match args[i].as_str() {
            "--dmin" => { dmin = Some(val(i + 1, "--dmin")?); i += 2; }
            "--dmax" => { dmax = Some(val(i + 1, "--dmax")?); i += 2; }
            "--lambda" => { lambda = Some(val(i + 1, "--lambda")?); i += 2; }
            "--2theta" => { tt = Some((val(i + 1, "--2theta")?, val(i + 2, "--2theta")?)); i += 3; }
            a => { rest.push(a.to_string()); i += 1; }
        }
    }
    if rest.len() < 7 {
        return Err(format!("hkl needs a symbol and six cell parameters\n\n{USAGE}"));
    }
    let k = rest.len() - 6;
    let c: Vec<f64> = rest[k..].iter().map(|x| num(x, "cell parameter")).collect::<Result<_, _>>()?;
    let cell = Cell::new(c[0], c[1], c[2], c[3], c[4], c[5]).map_err(|e| e.to_string())?;
    let (g, _) = SpaceGroup::from_symbol(&rest[..k].join(" "))?;
    // A 2θ range becomes a d range through Bragg's law.
    if let Some((t1, t2)) = tt {
        let w = lambda.ok_or("--2theta needs --lambda")?;
        let d_of = |t: f64| w / (2.0 * (t.min(179.999).to_radians() / 2.0).sin());
        dmin = Some(d_of(t1.max(t2)));
        dmax = Some(d_of(t1.min(t2).max(1e-6)));
    }
    let dmin = dmin.ok_or("give --dmin, or --lambda with --2theta")?;
    let dmax = dmax.unwrap_or(f64::MAX);
    let list = drs::reflections::generate(&g, &cell, dmin, dmax, lambda);
    let name = g.setting.map_or(g.hall.clone(), |s| format!("{} (no. {}, Laue class {})", s.hm, s.number, g.laue));
    println!("Reflections of {name}: {} unique, d from {dmin:.5} Å{}", list.len(),
        if dmax < 1e6 { format!(" to {dmax:.5} Å") } else { String::new() });
    println!("{:>4} {:>4} {:>4} {:>5} {:>11}{}", "h", "k", "l", "mult", "d (Å)", if lambda.is_some() { "    2θ (°)" } else { "" });
    for r in &list {
        let t = r.two_theta.map_or(String::new(), |t| format!(" {t:9.3}"));
        println!("{:>4} {:>4} {:>4} {:>5} {:>11.6}{t}{}", r.hkl[0], r.hkl[1], r.hkl[2], r.multiplicity, r.d,
            if r.fullprof_choice { "" } else { " *" });
    }
    Ok(())
}

fn cmd_sf(args: &[String]) -> Result<(), String> {
    use drs::scattering::{anode, Radiation, Scatterer};
    use drs::structure_factor::{chemical_from_occupancy, occupancy_from_chemical, Atom, Structure};
    if wants_help(args) {
        print!("{USAGE}");
        return Ok(());
    }
    let (mut dmin, mut dmax, mut neutron, mut anomalous) = (None, f64::MAX, false, true);
    let mut xray: Option<String> = None;
    let mut atoms_raw: Vec<Vec<String>> = Vec::new();
    let mut rest: Vec<String> = Vec::new();
    let val = |i: usize, what: &str| -> Result<f64, String> { num(args.get(i).ok_or(format!("{what} needs a value"))?, what) };
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--dmin" => { dmin = Some(val(i + 1, "--dmin")?); i += 2; }
            "--dmax" => { dmax = val(i + 1, "--dmax")?; i += 2; }
            "--neutron" => { neutron = true; i += 1; }
            "--no-anomalous" => { anomalous = false; i += 1; }
            "--xray" => { xray = Some(args.get(i + 1).ok_or("--xray needs a wavelength or an anode")?.clone()); i += 2; }
            "--atom" => {
                let mut j = i + 1;
                let mut a = Vec::new();
                while j < args.len() && !args[j].starts_with("--") && a.len() < 7 {
                    a.push(args[j].clone());
                    j += 1;
                }
                if a.len() < 6 {
                    return Err("--atom needs <label> <species> <x> <y> <z> <Biso> [<occ>]".into());
                }
                atoms_raw.push(a);
                i = j;
            }
            a => { rest.push(a.to_string()); i += 1; }
        }
    }
    let dmin = dmin.ok_or("give --dmin")?;
    if neutron && xray.is_some() {
        return Err("--xray and --neutron exclude each other".into());
    }
    let (radiation, rad_text) = if neutron {
        (Radiation::Neutron, "neutrons (scattering lengths, 10⁻¹² cm)".to_string())
    } else {
        let x = xray.unwrap_or_else(|| "Cu".into());
        let (lambda, name) = match anode(&x) {
            Some(a) => (a.1, format!("{} Kα1", a.0)),
            None => (num(&x, "--xray")?, format!("λ = {x} Å")),
        };
        let k = drs::scattering::nearest_anode(lambda);
        let title = |n: &str| n[..1].to_string() + &n[1..].to_ascii_lowercase();
        let name = anode(&x).map_or(name, |a| format!("{} Kα1", title(a.0)));
        let anom = if anomalous { format!(", f′ and f″ of {} Kα1", title(drs::scattering_table::ANODES[k].0)) } else { ", no anomalous terms".into() };
        (Radiation::Xray { lambda, anomalous }, format!("X-rays, {name} ({lambda} Å){anom}"))
    };
    let st = if rest.len() == 1 && std::path::Path::new(&rest[0]).is_file() {
        if !atoms_raw.is_empty() {
            return Err("--atom is for a structure given on the command line, not with a CIF".into());
        }
        // Latin-1 when not UTF-8 (ICSD files are): every byte a character.
        let bytes = std::fs::read(&rest[0]).map_err(|e| format!("{}: {e}", rest[0]))?;
        let text = String::from_utf8(bytes.clone()).unwrap_or_else(|_| bytes.iter().map(|&b| b as char).collect());
        drs::cif::structure(&text).map_err(|e| format!("{}: {e}", rest[0]))?
    } else {
        if rest.len() < 7 {
            return Err(format!("sf needs a CIF file, or a symbol, six cell parameters and --atom ...\n\n{USAGE}"));
        }
        let k = rest.len() - 6;
        let c: Vec<f64> = rest[k..].iter().map(|x| num(x, "cell parameter")).collect::<Result<_, _>>()?;
        let cell = Cell::new(c[0], c[1], c[2], c[3], c[4], c[5]).map_err(|e| e.to_string())?;
        let group = SpaceGroup::from_fullprof_symbol(&rest[..k].join(" "))?;
        if atoms_raw.is_empty() {
            return Err("no atoms: give one --atom <label> <species> <x> <y> <z> <Biso> [<occ>] per atom".into());
        }
        let mut atoms = Vec::new();
        for a in &atoms_raw {
            let x = [coordinate(&a[2])?, coordinate(&a[3])?, coordinate(&a[4])?];
            let chemical = if a.len() > 6 { num(&a[6], "occupancy")? } else { 1.0 };
            let occ = occupancy_from_chemical(&group, x, chemical).map_err(|e| format!("atom {}: {e}", a[0]))?;
            atoms.push(Atom { label: a[0].clone(), species: a[1].clone(), x, biso: num(&a[5], "Biso")?, occ, beta: None, charge: None });
        }
        Structure { group, cell, atoms }
    };
    let mut scatterers = Vec::new();
    for a in &st.atoms {
        let (s, fell_back) = Scatterer::lookup_or_element(&a.species, radiation)?;
        if fell_back {
            println!("note: {} '{}' is not in the X-ray table; its neutral atom {} is used", a.label, a.species, drs::scattering::element_of(&a.species));
        }
        scatterers.push(s);
    }
    let g = &st.group;
    let name = g.setting.map_or(format!("Hall {}", g.hall), |s| format!("{} (no. {}, Hall {})", s.hm, s.number, s.hall));
    let c = &st.cell;
    println!("Structure factors of {name}, M = {}", g.multiplicity());
    println!("Cell {:.5} {:.5} {:.5} {:.4} {:.4} {:.4}; {rad_text}", c.a, c.b, c.c, c.alpha, c.beta, c.gamma);
    println!("{:<8} {:<6} {:>9} {:>9} {:>9} {:>8} {:>5} {:>8} {:>9}", "atom", "type", "x", "y", "z", "Biso", "m", "occ", "Occ(FP)");
    for a in &st.atoms {
        let m = drs::site::Site::new(g, a.x)?.multiplicity;
        let chem = chemical_from_occupancy(g, a.x, a.occ)?;
        println!("{:<8} {:<6} {:>9.5} {:>9.5} {:>9.5} {:>8.4} {:>5} {:>8.4} {:>9.5}{}", a.label, a.species, a.x[0], a.x[1], a.x[2], a.biso, m, chem, a.occ,
            if a.beta.is_some() { "  (anisotropic)" } else { "" });
    }
    let list = drs::reflections::generate(g, c, dmin, dmax, None);
    println!("{} unique reflections with d ≥ {dmin} Å", list.len());
    println!("{:>4} {:>4} {:>4} {:>5} {:>10} {:>12} {:>14}", "h", "k", "l", "mult", "d (Å)", "|F|", "F²");
    for r in &list {
        let f = st.structure_factor(&scatterers, r.hkl);
        println!("{:>4} {:>4} {:>4} {:>5} {:>10.5} {:>12.4} {:>14.4}", r.hkl[0], r.hkl[1], r.hkl[2], r.multiplicity, r.d, f.abs(), f.norm_sqr());
    }
    Ok(())
}

fn cmd_calc(args: &[String], simulate: bool) -> Result<(), String> {
    use drs::powder::{agreement_weighted, Random};
    use std::io::Write;
    if wants_help(args) {
        print!("{USAGE}");
        return Ok(());
    }
    let (mut out, mut seed, mut from, mut to, mut step, mut data) = (None, 1u64, None, None, None, None::<String>);
    let (mut file, mut fullprof, mut per_phase, mut cumulative, mut shares) = (None, false, false, false, false);
    let mut i = 0;
    while i < args.len() {
        let next = |i: usize, what: &str| args.get(i + 1).cloned().ok_or(format!("{what} needs a value"));
        match args[i].as_str() {
            "--fullprof" => { fullprof = true; i += 1; }
            "--phases" => { per_phase = true; i += 1; }
            "--cumulative" => { cumulative = true; i += 1; }
            "--shares" => { shares = true; i += 1; }
            "--out" => { out = Some(next(i, "--out")?); i += 2; }
            "--data" => { data = Some(next(i, "--data")?); i += 2; }
            "--seed" => { seed = next(i, "--seed")?.parse().map_err(|_| "--seed: an integer")?; i += 2; }
            "--from" => { from = Some(num(&next(i, "--from")?, "--from")?); i += 2; }
            "--to" => { to = Some(num(&next(i, "--to")?, "--to")?); i += 2; }
            "--step" => { step = Some(num(&next(i, "--step")?, "--step")?); i += 2; }
            a => { file = Some(a.to_string()); i += 1; }
        }
    }
    let file = file.ok_or(format!("give a .pcr or control file\n\n{USAGE}"))?;
    let (loaded, grid) = load_model(&file, data.as_deref(), fullprof)?;
    let x: Vec<f64> = match (&loaded.data, from, to, step, grid) {
        (_, Some(a), Some(b), Some(s), _) if s > 0.0 => (0..=((b - a) / s + 1e-9) as usize).map(|k| a + k as f64 * s).collect(),
        (Some(d), _, _, _, _) => d.x.clone(),
        (_, _, _, _, Some((a, s, b))) => (0..=((b - a) / s + 1e-9) as usize).map(|k| a + k as f64 * s).collect(),
        _ => return Err("no data file found and no Thmin, Step, Thmax: give --from, --to, --step".into()),
    };
    let c = loaded.model.calculate(&x);
    loaded.model.finite(&c).map_err(|e| format!("{file}: {e}"))?;
    let mut w: Box<dyn Write> = match &out {
        Some(f) => Box::new(std::io::BufWriter::new(std::fs::File::create(f).map_err(|e| format!("{f}: {e}"))?)),
        None => Box::new(std::io::BufWriter::new(std::io::stdout())),
    };
    let werr = |e: std::io::Error| e.to_string();
    let m = &loaded.model;
    let phases: Vec<String> = m.phases.iter().map(|p| format!("{} ({})", p.name, p.group.setting.map_or(p.group.hall.clone(), |s| s.hm.to_string()))).collect();
    if simulate {
        let mut r = Random::new(seed);
        writeln!(w, "XYDATA\nTITLE: simulated by drs {} from {file} (Poisson noise, seed {seed})\n! phases: {}\n! 2theta  Yobs  sigma\n!\n!", env!("CARGO_PKG_VERSION"), phases.join(", ")).map_err(werr)?;
        for (i, &t) in x.iter().enumerate() {
            let y = r.poisson(c.ycalc[i]);
            writeln!(w, "{t:12.5} {y:14.2} {:12.3}", y.max(1.0).sqrt()).map_err(werr)?;
        }
        if out.is_some() {
            eprintln!("{} points simulated from {file} ({} peaks; seed {seed})", x.len(), c.peaks.len());
        }
        return Ok(());
    }
    writeln!(w, "# drs {} calc {file}: {} ({} peaks); {}", env!("CARGO_PKG_VERSION"), phases.join(", "), c.peaks.len(), if m.neutrons { "neutrons" } else { "X-rays" }).map_err(werr)?;
    let kept = m.convention.fullprof_rules();
    writeln!(w, "# {}", if fullprof { "FullProf's rules, the .pcr's values as written".to_string() } else if kept.is_empty() { "translated into the default model".to_string() } else { format!("translated into the default model; FullProf's rules kept (no exact translation): {}", kept.join(", ")) }).map_err(werr)?;
    let obs = loaded.data.as_ref().filter(|d| d.len() == x.len() && from.is_none());
    if let Some(d) = obs {
        let (rp, rwp, rexp, chi2) = agreement_weighted(d, &c.ycalc, &m.excluded, 0, m.weighting);
        let line = format!("Rp {rp:.2}  Rwp {rwp:.2}  Rexp {rexp:.2}  chi2 {chi2:.3} ({} points, excluded regions left out; no parameter counted)", d.len());
        writeln!(w, "# {line}").map_err(werr)?;
        if out.is_some() {
            eprintln!("{line}");
        }
        let cum = cumulative.then(|| drs::diagnostics::cumulative(m, d, &c.ycalc));
        let frac = cum.as_ref().map(|c| c.fraction());
        let shares_head: String = if shares { std::iter::once("   bg share".to_string()).chain((1..=m.phases.len()).map(|k| format!(" phase{k} share"))).collect() } else { String::new() };
        writeln!(w, "# 2theta        Ycal          Backg          Yobs         sigma{}{}{}", phase_head(per_phase, m.phases.len()), if cumulative { "  cum.chi2/total  w^1/2(y-yc)" } else { "" }, shares_head).map_err(werr)?;
        for (i, xi) in x.iter().enumerate() {
            let extra = match (&cum, &frac) {
                (Some(c), Some(f)) => format!(" {:14.6} {:12.4}", f[i], c.weighted_difference[i]),
                _ => String::new(),
            };
            let sh = if shares { let s = drs::diagnostics::shares(&c, i); std::iter::once(s.background).chain(s.phases).map(|v| format!(" {v:12.5}")).collect() } else { String::new() };
            writeln!(w, "{:10.5} {:14.4} {:14.4} {:14.4} {:12.4}{}{extra}{sh}", xi, c.ycalc[i], c.background[i], d.y[i], d.sigma[i], phase_cols(per_phase, &c, i)).map_err(werr)?;
        }
    } else {
        writeln!(w, "# 2theta        Ycal          Backg{}", phase_head(per_phase, m.phases.len())).map_err(werr)?;
        for (i, xi) in x.iter().enumerate() {
            writeln!(w, "{:10.5} {:14.4} {:14.4}{}", xi, c.ycalc[i], c.background[i], phase_cols(per_phase, &c, i)).map_err(werr)?;
        }
    }
    Ok(())
}

/// Thmin, Step, Thmax of a `.pcr`.
type Grid = (f64, f64, f64);

/// A .pcr or control file with its data (or `data` given): the model and, for a .pcr, Thmin, Step, Thmax.
fn load_model(file: &str, data: Option<&str>, fullprof: bool) -> Result<(drs::powder::Loaded, Option<Grid>), String> {
    use drs::powder::Model;
    let path = std::path::Path::new(file);
    if drs::control::is_control_file(path) {
        if fullprof {
            return Err("--fullprof reads a .pcr; a control file states its own rules in [convention]".into());
        }
        let o = open_control(path)?;
        let data = match data {
            Some(f) => Some(drs::formats::read_path(std::path::Path::new(f), None).map_err(|e| e.to_string())?.scans.remove(0).pattern),
            None => o.data,
        };
        Ok((drs::powder::Loaded { model: o.control.model, data, data_path: o.data_path }, None))
    } else {
        let bytes = std::fs::read(path).map_err(|e| format!("{file}: {e}"))?;
        let text: String = String::from_utf8(bytes.clone()).unwrap_or_else(|_| bytes.iter().map(|&b| b as char).collect());
        let pcr = drs::pcr::read(&text).map_err(|e| format!("{file}: {e}"))?;
        // The absolute folder, so that its parent can be searched for the data and resolution files.
        let dir = std::fs::canonicalize(path).ok().and_then(|p| p.parent().map(|d| d.to_path_buf())).unwrap_or_else(|| std::path::PathBuf::from("."));
        let dir = dir.as_path();
        let stem = path.file_stem().and_then(|s| s.to_str()).unwrap_or("");
        let data_file = data.map(std::path::Path::new);
        let l = if fullprof { Model::from_pcr_fullprof_with_data(&pcr, dir, stem, data_file)? } else { Model::from_pcr_with_data(&pcr, dir, stem, data_file)? };
        Ok((l, (pcr.step > 0.0).then_some((pcr.thmin, pcr.step, pcr.thmax))))
    }
}

/// `drs reflections`: each reflection with its width, I calc, I obs partitioned by the model, σ and residual; the
/// forbidden ones; a mode phase's parent indices and k-vectors (W8, app requests 19–20, L100).
fn cmd_reflections(args: &[String]) -> Result<(), String> {
    if wants_help(args) || args.is_empty() {
        print!("{USAGE}");
        return Ok(());
    }
    let (mut file, mut data, mut fullprof, mut forbidden) = (None, None::<String>, false, false);
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--fullprof" => (fullprof, i) = (true, i + 1),
            "--forbidden" => (forbidden, i) = (true, i + 1),
            "--data" => (data, i) = (Some(args.get(i + 1).cloned().ok_or("--data needs a file")?), i + 2),
            a => (file, i) = (Some(a.to_string()), i + 1),
        }
    }
    let file = file.ok_or(format!("give a .pcr or control file\n\n{USAGE}"))?;
    let (l, _) = load_model(&file, data.as_deref(), fullprof)?;
    let m = &l.model;
    let d = l.data.as_ref().ok_or("no data: give --data (I obs is partitioned from the observed pattern)")?;
    if forbidden {
        println!("# phase\th\tk\tl\tmult\td\t2theta\tcondition");
        for (ip, r) in drs::diagnostics::forbidden_in(m, &d.x) {
            println!("{}\t{}\t{}\t{}\t{}\t{:.5}\t{}\t{}", ip + 1, r.hkl[0], r.hkl[1], r.hkl[2], r.multiplicity, r.d, r.two_theta.map_or("-".into(), |t| format!("{t:.4}")), r.condition);
        }
        return Ok(());
    }
    println!("# phase\twavelength\th\tk\tl\tmult\t2theta\td\tFWHM\teta\tF2\tIcalc\tIobs\tsigma\t(Iobs-Icalc)/sigma\tparent hkl\tk");
    for r in drs::diagnostics::reflections(m, d) {
        let parent = m.phases[r.phase].modes.as_ref().and_then(|ms| ms.parent.as_ref()).and_then(|p| p.reflection(r.hkl).ok());
        let (ph, kv) = parent.map_or(("-".to_string(), "-".to_string()), |p| (format!("{:.3} {:.3} {:.3}", p.hkl[0], p.hkl[1], p.hkl[2]), format!("{:.3} {:.3} {:.3}", p.k[0], p.k[1], p.k[2])));
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.5}\t{:.5}\t{:.4}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.2}\t{ph}\t{kv}",
            r.phase + 1, r.wavelength, r.hkl[0], r.hkl[1], r.hkl[2], r.multiplicity, r.two_theta, r.d, r.fwhm, r.eta, r.f2.map_or("-".into(), |f| format!("{f:.3}")), r.icalc, r.iobs, r.sigma, r.residual
        );
    }
    Ok(())
}

/// `drs mic`: the microstructure quantities of every reflection of the TCH phases of a .pcr.
/// The header of the per-phase columns of `drs calc --phases`.
fn phase_head(on: bool, n: usize) -> String {
    if on { (1..=n).map(|k| format!("      phase{k:<6}")).collect() } else { String::new() }
}

/// Each phase's curve at point i (no background), for `drs calc --phases`.
fn phase_cols(on: bool, c: &drs::powder::Calculated, i: usize) -> String {
    if on { c.per_phase.iter().map(|v| format!(" {:14.4}", v[i])).collect() } else { String::new() }
}

/// Set by Ctrl-C during `drs refine`: the engine stops after the cycle it is in (`lsq::refine_until`).
static INTERRUPTED: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

#[cfg(unix)]
extern "C" fn on_interrupt(_: i32) {
    INTERRUPTED.store(true, std::sync::atomic::Ordering::SeqCst);
}

fn cmd_mic(args: &[String]) -> Result<(), String> {
    use drs::powder::Model;
    use std::io::Write;
    if wants_help(args) {
        print!("{USAGE}");
        return Ok(());
    }
    let (mut out, mut file, mut fullprof) = (None, None, false);
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--fullprof" => {
                fullprof = true;
                i += 1;
            }
            "--out" => {
                out = Some(args.get(i + 1).cloned().ok_or("--out needs a value")?);
                i += 2;
            }
            a => {
                file = Some(a.to_string());
                i += 1;
            }
        }
    }
    let file = file.ok_or(format!("give a .pcr file\n\n{USAGE}"))?;
    let path = std::path::Path::new(&file);
    let bytes = std::fs::read(path).map_err(|e| format!("{file}: {e}"))?;
    let text: String = String::from_utf8(bytes.clone()).unwrap_or_else(|_| bytes.iter().map(|&b| b as char).collect());
    let mut pcr = drs::pcr::read(&text).map_err(|e| format!("{file}: {e}"))?;
    // Magnetic phases (Jbt = ±1) are left out, and said: their intensities are not calculated (L211).
    for (k, ph) in pcr.phases.iter().enumerate() {
        if ph.flag("Jbt").abs() == 1 {
            eprintln!("{file}: phase {} is magnetic (Jbt = {}): left out, its intensities are not calculated", k + 1, ph.flag("Jbt"));
        }
    }
    pcr.phases.retain(|ph| ph.flag("Jbt").abs() != 1);
    let dir = std::fs::canonicalize(path).ok().and_then(|p| p.parent().map(|d| d.to_path_buf())).unwrap_or_else(|| std::path::PathBuf::from("."));
    let stem = path.file_stem().and_then(|s| s.to_str()).unwrap_or("");
    let m = if fullprof { Model::from_pcr_fullprof(&pcr, &dir, stem)?.model } else { Model::from_pcr(&pcr, &dir, stem)?.model };
    let mut w: Box<dyn Write> = match &out {
        Some(f) => Box::new(std::io::BufWriter::new(std::fs::File::create(f).map_err(|e| format!("{f}: {e}"))?)),
        None => Box::new(std::io::BufWriter::new(std::io::stdout())),
    };
    let werr = |e: std::io::Error| e.to_string();
    writeln!(w, "# drs {} mic {file}", env!("CARGO_PKG_VERSION")).map_err(werr)?;
    writeln!(w, "#  s=1/d     HGo       HGi       HLo       HLi       betaG     betaL    App-size  Max-strain    h    k    l    2theta     HpVo      eta   Beta_TCH").map_err(werr)?;
    let mut peaks = m.peaks(pcr.thmin, pcr.thmax);
    peaks.sort_by(|a, b| (a.phase, a.two_theta).partial_cmp(&(b.phase, b.two_theta)).unwrap());
    for p in &peaks {
        let Some(b) = m.breadths(p.phase, p.hkl, p.wavelength) else { continue };
        let size = b.size.filter(|&v| v < 99999.0).unwrap_or(99999.0);
        writeln!(
            w,
            "{:9.4} {:9.4} {:9.4} {:9.4} {:9.4} {:9.4} {:9.4} {:10.2} {:10.4} {:4} {:4} {:4} {:9.4} {:8.4} {:8.4} {:9.4}",
            b.d_star, b.hg_obs, b.hg_ins, b.hl_obs, b.hl_ins, b.beta_g, b.beta_l, size, b.strain.unwrap_or(0.0), p.hkl[0], p.hkl[1], p.hkl[2], b.two_theta, b.h_obs, b.eta_obs, b.beta_tch
        )
        .map_err(werr)?;
    }
    Ok(())
}

/// `drs refine`: least-squares refinement of a .pcr's coded parameters.
/// `drs files`: the file dictionary (R27, L091), or the kind of a file of a project.
fn cmd_files(args: &[String]) -> Result<(), String> {
    use drs::formats::files::{kind_of, FILES};
    if wants_help(args) {
        print!("{USAGE}");
        return Ok(());
    }
    if let Some(rel) = args.first() {
        let k = kind_of(rel).ok_or(format!("{rel}: no kind of the file dictionary has this place and name"))?;
        println!("{rel}: {} ({}), {}: {}", k.name, k.id, k.role, k.contents);
        return Ok(());
    }
    for k in FILES {
        let place = if k.folder.is_empty() { k.pattern.to_string() } else { format!("{}/{}", k.folder, k.pattern) };
        println!("{:42} {:28} {:8} {}", place, k.id, k.role, k.status);
    }
    Ok(())
}

/// `drs fullprof`: FullProf's `.prf` or `.hkl`, read (app request 17, L100).
fn cmd_fullprof(args: &[String]) -> Result<(), String> {
    use drs::formats::fullprof_outputs::{read_hkl, read_prf};
    if wants_help(args) || args.is_empty() {
        print!("{USAGE}");
        return Ok(());
    }
    let path = &args[0];
    let text: String = std::fs::read(path).map_err(|e| format!("{path}: {e}"))?.iter().map(|&b| b as char).collect();
    if path.to_lowercase().ends_with(".hkl") {
        for b in read_hkl(&text).map_err(|e| format!("{path}: {}", e.0))? {
            println!("# pattern {} phase {}: {} reflections — {}", b.pattern.map_or("-".into(), |v| v.to_string()), b.phase.map_or("-".into(), |v| v.to_string()), b.rows.len(), b.heading);
            println!("{}", b.columns.join("\t"));
            for r in &b.rows {
                println!("{}", r.iter().map(|v| v.to_string()).collect::<Vec<_>>().join("\t"));
            }
        }
        return Ok(());
    }
    let p = read_prf(&text).map_err(|e| format!("{path}: {}", e.0))?;
    println!("# {}{}", p.title, p.chi2.map_or(String::new(), |c| format!("  (chi2 {c})")));
    println!("# {} phase(s), {} points; lambda {} {}; Zero {} SyCos {} SySin {}", p.phases, p.x.len(), p.wavelengths[0], p.wavelengths[1], p.zero, p.sycos, p.sysin);
    for k in &p.k_vectors {
        println!("# propagation vector {} {} {}", k[0], k[1], k[2]);
    }
    for (a, b) in &p.excluded {
        println!("# excluded {a} {b}");
    }
    println!("# background = Backg + {} (L082)", p.background_offset);
    println!("2theta\tYobs\tYcal\tbackground");
    for i in 0..p.x.len() {
        println!("{}\t{}\t{}\t{}", p.x[i], p.yobs[i], p.ycalc[i], p.background[i]);
    }
    println!("phase\twavelength\th\tk\tl\t2theta");
    for r in &p.reflections {
        println!("{}\t{}\t{}\t{}\t{}\t{}", r.phase, r.wavelength, r.hkl[0], r.hkl[1], r.hkl[2], r.two_theta);
    }
    Ok(())
}

/// `drs modes`: a phase's symmetry modes (R29, L090), or the analyses of an AMPLIMODES output.
fn cmd_modes(args: &[String]) -> Result<(), String> {
    use drs::modes::{self, Modes};
    if wants_help(args) || args.is_empty() {
        print!("{USAGE}");
        return Ok(());
    }
    let (mut file, mut ampli, mut analysis, mut project, mut shift) = (None::<String>, None::<String>, 1usize, None::<String>, [0.0; 3]);
    let (mut equivalent, mut reference) = (false, None::<String>);
    let mut i = 0;
    while i < args.len() {
        let next = |i: usize, what: &str| args.get(i + 1).cloned().ok_or(format!("{what} needs a value"));
        match args[i].as_str() {
            "--amplimodes" => (ampli, i) = (Some(next(i, "--amplimodes")?), i + 2),
            "--analysis" => (analysis, i) = (next(i, "--analysis")?.parse().map_err(|_| "--analysis: a number from 1")?, i + 2),
            "--project" => (project, i) = (Some(next(i, "--project")?), i + 2),
            "--equivalent" => (equivalent, i) = (true, i + 1),
            "--reference" => (reference, i) = (Some(next(i, "--reference")?), i + 2),
            "--shift" => {
                let v: Vec<f64> = next(i, "--shift")?.split(',').filter_map(|t| t.trim().parse().ok()).collect();
                shift = <[f64; 3]>::try_from(v).map_err(|_| "--shift: x,y,z")?;
                i += 2;
            }
            a => (file, i) = (Some(a.to_string()), i + 1),
        }
    }
    let file = file.ok_or(format!("give a .pcr, a control file or an AMPLIMODES output\n\n{USAGE}"))?;
    let path = std::path::Path::new(&file);
    let read = |f: &str| -> Result<String, String> { std::fs::read(f).map(|b| String::from_utf8(b.clone()).unwrap_or_else(|_| b.iter().map(|&c| c as char).collect())).map_err(|e| format!("{f}: {e}")) };
    let text = read(&file)?;
    let analyses = modes::read_amplimodes(&text);
    if !analyses.is_empty() && !drs::control::is_control_file(path) && drs::pcr::read(&text).is_err() {
        for (k, a) in analyses.iter().enumerate() {
            println!("analysis {}: parent {}, distorted {}{}", k + 1, a.high.as_ref().map_or("?".into(), |h| format!("group {} cell {:?}", h.group, h.cell)), a.low.as_ref().map_or("?".into(), |l| format!("group {} cell {:?}", l.group, l.cell)), a.low_setting.as_ref().map_or(String::new(), |s| format!(" ({s})")));
            if let Some(c) = a.reference_cell() {
                println!("  reference metric: {:.6} {:.6} {:.6} {:.4} {:.4} {:.4}", c.a, c.b, c.c, c.alpha, c.beta, c.gamma);
            }
            if let Some(pa) = a.mode_parent() {
                // The keys of the phase's table in a control file: its reflections' parent indices and k-vectors (L100).
                let rows: Vec<String> = pa.basis.iter().map(|r| format!("[{}, {}, {}]", r[0], r[1], r[2])).collect();
                println!("  control file: mode_parent_group = \"{}\"; mode_parent_basis = [{}]; mode_parent_origin = [{}, {}, {}]", pa.group, rows.join(", "), pa.origin[0], pa.origin[1], pa.origin[2]);
            }
            for s in &a.summary {
                println!("  {:8} {:8} {:16} {:18} dim {}  {:.4} Å", s.k, s.irrep, s.direction, s.isotropy, s.dimension, s.amplitude);
            }
            if let Some(g) = a.global {
                println!("  global distortion {g} Å; {} atom mappings{}", a.mappings.len(), if a.fullprof.is_some() { "; a block for FullProf" } else { "" });
            }
        }
        return Ok(());
    }
    let model = if drs::control::is_control_file(path) {
        open_control(path)?.control.model
    } else {
        let pcr = drs::pcr::read(&text).map_err(|e| format!("{file}: {e}"))?;
        let dir = std::fs::canonicalize(path).ok().and_then(|p| p.parent().map(|d| d.to_path_buf())).unwrap_or_else(|| std::path::PathBuf::from("."));
        let stem = path.file_stem().and_then(|s| s.to_str()).unwrap_or("");
        drs::powder::Model::from_pcr(&pcr, &dir, stem)?.model
    };
    let mut ampli_parent = None;
    let metric = match &ampli {
        Some(f) => {
            let all = modes::read_amplimodes(&read(f)?);
            let a = all.get(analysis - 1).ok_or(format!("{f}: no analysis {analysis}"))?;
            // AMPLIMODES's summary per irrep (Publish P9): k-vector, irrep, direction, isotropy subgroup, dimension, amplitude.
            println!("{f}, analysis {analysis}: AMPLIMODES's summary");
            for s in &a.summary {
                println!("  {:8} {:8} {:16} {:18} dim {}  {:.4} Å", s.k, s.irrep, s.direction, s.isotropy, s.dimension, s.amplitude);
            }
            ampli_parent = a.mode_parent();
            Some(a.reference_cell().ok_or(format!("{f}: analysis {analysis} has no parent cell and transformation"))?)
        }
        None => None,
    };
    let mut any = false;
    for (p, ph) in model.phases.iter().enumerate() {
        let (Some(ms), drs::powder::Intensities::Structure(st, _)) = (&ph.modes, &ph.intensities) else { continue };
        any = true;
        let metric = metric.or(ms.reference_cell);
        println!("phase {} ({}): {} modes on {} atoms", p + 1, ph.name.trim(), ms.modes.len(), st.atoms.len());
        let gram = metric.map(|c| ms.gram(&ph.group, &c)).transpose()?;
        for (k, m) in ms.modes.iter().enumerate() {
            let atoms: Vec<&str> = m.vectors.iter().map(|v| st.atoms[v.atom].label.as_str()).collect();
            let norm = gram.as_ref().map_or(String::new(), |g| format!("  Σ|cV|² {:.5} Å²", g[k][k]));
            println!("  {:12} {:8} {:10.6} Å  on {}{norm}", m.name, m.irrep, m.amplitude, atoms.join(" "));
        }
        println!("  irreps: amplitude (Å) and polarisation vector");
        for ir in ms.irreps() {
            println!("  {:8} {:10.6}   {:?}", ir.name, ir.amplitude, ir.polarisation.iter().map(|c| format!("{c:.4}")).collect::<Vec<_>>());
        }
        println!("  positions: reference + displacement");
        for (j, a) in st.atoms.iter().enumerate() {
            let u = ms.displacement(j);
            let len = drs::mat3::quad(&ph.cell.metric(), &u).sqrt();
            println!("  {:8} {:10.6} {:10.6} {:10.6}   u {:9.6} {:9.6} {:9.6}  {len:.5} Å", a.label, a.x[0], a.x[1], a.x[2], u[0], u[1], u[2]);
        }
        if let Some(f) = &project {
            let metric = metric.ok_or("--project needs the reference metric: give --amplimodes <file> [--analysis k], or a control file with mode_reference_cell")?;
            let other = drs::pcr::read(&read(f)?).map_err(|e| format!("{f}: {e}"))?;
            let q = other.phases.get(p).ok_or(format!("{f} has no phase {}", p + 1))?;
            let reference: Vec<(String, [f64; 3])> = st.atoms.iter().zip(&ms.reference).map(|(a, x)| (a.species.clone(), *x)).collect();
            let distorted: Vec<(String, [f64; 3])> = q.atoms.iter().map(|a| (a.species.clone(), a.x)).collect();
            let u = modes::displacements(&reference, &distorted, &ph.group, &metric, shift, 1.0)?;
            let (amps, residual) = ms.project(&ph.group, &metric, &u)?;
            let mut pm: Modes = ms.clone();
            for (m, a) in pm.modes.iter_mut().zip(&amps) {
                m.amplitude = *a;
            }
            println!("  {f} projected onto the modes (shift {shift:?}): the part no mode describes {residual:.2e} Å");
            for m in &pm.modes {
                println!("  {:12} {:10.6} Å", m.name, m.amplitude);
            }
            for ir in pm.irreps() {
                println!("  {:8} {:10.6}", ir.name, ir.amplitude);
            }
        }
        // E7: the amplitude sets that describe the same structure up to the parent's symmetry, and the one nearest a
        // reference's amplitudes (another file's mode phase, by mode name) as the sign convention.
        if equivalent || reference.is_some() {
            let species: Vec<String> = st.atoms.iter().map(|a| a.species.clone()).collect();
            let m3 = metric.unwrap_or(ph.cell);
            // The parent: the phase's, else the AMPLIMODES analysis given.
            let mut ms = ms.clone();
            if ms.parent.is_none() {
                ms.parent = ampli_parent.clone();
            }
            let q: Vec<f64> = ms.modes.iter().map(|m| m.amplitude).collect();
            if equivalent {
                let eq = ms.equivalent_amplitudes(&species, &ph.group, &m3, &q)?;
                println!("  {} amplitude sets describe this structure (the parent's operation and translation):", eq.len());
                for e in &eq {
                    println!("  {:24} {:?}  {}", e.operation, e.translation, e.amplitudes.iter().map(|v| format!("{v:+.5}")).collect::<Vec<_>>().join(" "));
                }
            }
            if let Some(f) = &reference {
                let rp = std::path::Path::new(f);
                let other = if drs::control::is_control_file(rp) {
                    open_control(rp)?.control.model
                } else {
                    let pcr = drs::pcr::read(&read(f)?).map_err(|e| format!("{f}: {e}"))?;
                    let dir = std::fs::canonicalize(rp).ok().and_then(|p| p.parent().map(|d| d.to_path_buf())).unwrap_or_else(|| std::path::PathBuf::from("."));
                    drs::powder::Model::from_pcr(&pcr, &dir, rp.file_stem().and_then(|s| s.to_str()).unwrap_or(""))?.model
                };
                let oms = other.phases.get(p).and_then(|o| o.modes.as_ref()).ok_or(format!("{f}: phase {} has no modes", p + 1))?;
                let refq: Vec<f64> = ms.modes.iter().map(|m| oms.modes.iter().find(|o| o.name == m.name).map(|o| o.amplitude).ok_or(format!("{f}: no mode {}", m.name))).collect::<Result<_, _>>()?;
                let near = ms.nearest_to(&species, &ph.group, &m3, &q, &refq)?;
                println!("  nearest {f}'s amplitudes: {} {:?}", near.operation, near.translation);
                for (m, a) in ms.modes.iter().zip(&near.amplitudes) {
                    println!("  {:12} {a:+10.6} Å", m.name);
                }
            }
        }
        // E14b: the control file's sign reference, and the state against it.
        if let Some(sr) = &ms.sign_reference {
            let species: Vec<String> = st.atoms.iter().map(|a| a.species.clone()).collect();
            let q: Vec<f64> = ms.modes.iter().map(|m| m.amplitude).collect();
            let m3 = metric.unwrap_or(ph.cell);
            match ms.reference_signs(&species, &ph.group, &m3, &q) {
                Some(Ok(e)) => println!("  sign reference {}{}: {}", sr.file, sr.analysis.map_or(String::new(), |k| format!(", analysis {k}")), e.words()),
                Some(Err(e)) => println!("  sign reference {}: {e}", sr.file),
                None => {}
            }
        }
    }
    // E14c: FullProf's reported amplitudes, from the .sum beside a .pcr.
    if !drs::control::is_control_file(path) {
        let sum = path.with_extension("sum");
        if sum.is_file() {
            for (k, b) in modes::read_fullprof_amplitudes(&read(&sum.to_string_lossy())?).iter().enumerate() {
                println!("FullProf's reported amplitudes ({}, mode phase {}):", sum.display(), k + 1);
                for (n, v, e) in &b.modes {
                    println!("  {n:12} {v:+10.6} ({e})");
                }
                for (n, v, e) in &b.irreps {
                    println!("  {n:8} {v:.5}{}", e.map_or(String::new(), |e| format!(" ({e})")));
                }
            }
        }
    }
    if !any {
        return Err(format!("{file}: no phase is described by symmetry modes"));
    }
    Ok(())
}

fn cmd_refine(args: &[String]) -> Result<(), String> {
    use drs::lsq::{Options, Session};
    use drs::powder::Model;
    use drs::refinable::ParameterSet;
    if wants_help(args) {
        print!("{USAGE}");
        return Ok(());
    }
    if args.first().is_some_and(|f| is_multipattern_file(f)) {
        // A multi-pattern .pcr: the joint refinement (P06).
        return cmd_joint(args);
    }
    let (mut file, mut cycles, mut data, mut session_file) = (None, None, None::<String>, None::<String>);
    let (mut plan, mut only_step, mut each) = (false, None::<usize>, false);
    let (mut trace, mut map, mut span, mut grid) = (None::<String>, None::<(String, String)>, "3esd".to_string(), 21usize);
    let mut i = 0;
    while i < args.len() {
        let next = |i: usize, what: &str| args.get(i + 1).cloned().ok_or(format!("{what} needs a value"));
        match args[i].as_str() {
            "--trace" => {
                trace = Some(next(i, "--trace")?);
                i += 2;
            }
            "--map" => {
                let b = args.get(i + 2).cloned().ok_or("--map needs two parameters")?;
                map = Some((next(i, "--map")?, b));
                i += 3;
            }
            "--span" => {
                span = next(i, "--span")?;
                i += 2;
            }
            "--grid" => {
                grid = next(i, "--grid")?.parse::<usize>().ok().filter(|&n| n >= 2).ok_or("--grid: a number of points from 2")?;
                i += 2;
            }
            "--cycles" => {
                cycles = Some(next(i, "--cycles")?.parse::<usize>().map_err(|_| "--cycles: a positive integer")?);
                i += 2;
            }
            "--data" => {
                data = Some(next(i, "--data")?);
                i += 2;
            }
            "--session" => {
                session_file = Some(next(i, "--session")?);
                i += 2;
            }
            "--plan" => {
                plan = true;
                i += 1;
            }
            "--step" => {
                only_step = Some(next(i, "--step")?.parse::<usize>().ok().filter(|&k| k >= 1).ok_or("--step: a step number from 1")?);
                i += 2;
            }
            "--parameters" => {
                each = true;
                i += 1;
            }
            a => {
                file = Some(a.to_string());
                i += 1;
            }
        }
    }
    let file = file.ok_or(format!("give a .pcr or control file\n\n{USAGE}"))?;
    let path = std::path::Path::new(&file);
    let (loaded, set, ncy, eps, refinement) = if drs::control::is_control_file(path) {
        let o = open_control(path)?;
        let set = o.control.refinement.build(&o.control.model)?;
        let data = match data.as_deref() {
            Some(f) => Some(drs::formats::read_path(std::path::Path::new(f), None).map_err(|e| e.to_string())?.scans.remove(0).pattern),
            None => o.data,
        };
        let (n, e) = (o.control.refinement.cycles.unwrap_or(10), o.control.refinement.eps.unwrap_or(0.1));
        (drs::powder::Loaded { model: o.control.model, data, data_path: o.data_path }, set, n, e, Some(o.control.refinement))
    } else {
        let bytes = std::fs::read(path).map_err(|e| format!("{file}: {e}"))?;
        let text: String = String::from_utf8(bytes.clone()).unwrap_or_else(|_| bytes.iter().map(|&b| b as char).collect());
        let pcr = drs::pcr::read(&text).map_err(|e| format!("{file}: {e}"))?;
        let dir = std::fs::canonicalize(path).ok().and_then(|p| p.parent().map(|d| d.to_path_buf())).unwrap_or_else(|| std::path::PathBuf::from("."));
        let stem = path.file_stem().and_then(|s| s.to_str()).unwrap_or("");
        let loaded = Model::from_pcr_with_data(&pcr, &dir, stem, data.as_deref().map(std::path::Path::new))?;
        let set = ParameterSet::from_pcr(&pcr, &loaded.model)?;
        (loaded, set, pcr.ncy.max(1) as usize, pcr.eps, None)
    };
    let obs = loaded.data.as_ref().ok_or("no data file found: give --data")?;
    let mut model = loaded.model.clone();
    // The plan's steps to run (W12): all of them (--plan), or one (--step k); else the refinement as written.
    let steps: Vec<usize> = match (&refinement, plan, only_step) {
        (_, false, None) => Vec::new(),
        (Some(rf), _, _) if rf.plan.is_empty() => return Err(format!("{file} has no plan ([[refine.step]] tables)")),
        (Some(rf), _, Some(k)) if k > rf.plan.len() => return Err(format!("{file}'s plan has {} steps", rf.plan.len())),
        (Some(_), _, Some(k)) => vec![k - 1],
        (Some(rf), true, None) => (0..rf.plan.len()).collect(),
        (None, _, _) => return Err("a plan is written in a control file (.drc): a .pcr has none".into()),
    };
    if steps.is_empty() && set.free.is_empty() && !model.has_extracted() {
        return Err(format!("{file}: no parameter has a refinement code"));
    }
    // The phase names without the .pcr's keywords on the same line ("Name: Rutile  VARY xyz b").
    let phases: Vec<String> = model
        .phases
        .iter()
        .map(|p| {
            let n = p.name.trim().trim_start_matches("Name:").trim();
            n.split("  ").next().unwrap_or(n).split(" VARY").next().unwrap_or(n).split(" FIX").next().unwrap_or(n).trim().to_string()
        })
        .collect();
    println!("drs {} refine {file}: {}; {} points, {} parameters", env!("CARGO_PKG_VERSION"), phases.join(", "), obs.len(), set.free.len());
    let kept = model.convention.fullprof_rules();
    if !kept.is_empty() {
        println!("FullProf's rules kept (no exact translation): {}", kept.join(", "));
    }
    // Ctrl-C stops the refinement after the cycle in progress (SIGINT = 2), with the results at that cycle.
    #[cfg(unix)]
    unsafe {
        signal(2, on_interrupt as *const () as usize);
    }
    let on_cycle = |c: &drs::lsq::Cycle| {
        let big = c.params.iter().filter(|p| p.esd > 0.0).map(|p| ((p.computed / p.esd).abs(), p.name.as_str())).fold((0.0, ""), |a, b| if b.0 > a.0 { b } else { a });
        let s = &c.stats;
        println!("{:>5} {:8.2} {:8.2} {:8.2} {:10.4} {:8.4} {:9.0e}  {:.3} ({}){}", c.cycle, s.rp, s.rwp, s.rexp, s.chi2, s.durbin_watson, c.marquardt, big.0, big.1, if c.singular > 0 { format!("; {} singular directions left out", c.singular) } else { String::new() });
        for n in &c.notes {
            println!("      {n}");
        }
        if each {
            for p in &c.params {
                println!("      {:<43} {:16.6} {:12.6} shift {:12.6} ({:+.3} esd){}", p.name, p.value, p.esd, p.shift, p.shift_esd, if p.at_limit { "  at its limit" } else { "" });
            }
        }
        !INTERRUPTED.load(std::sync::atomic::Ordering::SeqCst)
    };
    let header = || println!("{:>5} {:>8} {:>8} {:>8} {:>10} {:>8} {:>9}  largest |computed shift|/esd", "cycle", "Rp", "Rwp", "Rexp", "chi2", "DW", "λ");
    let (r, trajectory) = if steps.is_empty() {
        let max_cycles = if session_file.is_some() { ncy } else { cycles.unwrap_or(ncy) }.max(1);
        let opts = Options { max_cycles, eps };
        let mut session = match session_file.as_deref().filter(|f| std::path::Path::new(f).exists()) {
            Some(f) => {
                let text = std::fs::read_to_string(f).map_err(|e| format!("{f}: {e}"))?;
                let s = Session::resume(&text, model.clone(), obs.clone(), set.clone(), opts)?;
                println!("resumed from {f} after cycle {} (λ = {})", s.cycle, s.lambda);
                s
            }
            None => Session::new(model.clone(), obs.clone(), set.clone(), opts),
        };
        for n in &session.notes {
            println!("{n}");
        }
        header();
        match (&session_file, cycles) {
            (Some(_), Some(n)) => {
                if session.step(n, on_cycle).is_none() {
                    println!("Paused after cycle {}", session.cycle);
                }
            }
            _ => {
                session.run(on_cycle);
            }
        }
        if let Some(f) = &session_file {
            std::fs::write(f, session.to_text()).map_err(|e| format!("{f}: {e}"))?;
            println!("session saved in {f}");
        }
        model = session.model.clone();
        (session.finish(), session.trajectory())
    } else {
        let rf = refinement.clone().unwrap();
        let on_step = |s: &drs::lsq::PlanStepRun| match s.stop {
            None => {
                println!("\nstep {} of {}: {} ({} parameters, up to {} cycles)", s.index + 1, rf.plan.len(), s.name, s.parameters, s.max_cycles);
                header();
            }
            Some(stop) => {
                println!("step {}: {}", s.index + 1, stop.text());
                if matches!(stop, drs::lsq::Stop::User | drs::lsq::Stop::Diverging(_) | drs::lsq::Stop::Singular { .. }) {
                    println!("the plan stops at step {}", s.index + 1);
                }
            }
        };
        let (r, session) = drs::lsq::run_plan(&mut model, obs, &rf, &steps, on_step, on_cycle)?;
        (r, session.trajectory())
    };
    let set = &trajectory.set;
    println!("{} after {} cycles: {}", if r.converged { "Stopped" } else { "Not converged" }, r.cycles, r.stop);
    if let drs::lsq::Stop::Singular { why, tip, .. } = &r.reason {
        println!("  why: {why}\n  what to do: {tip}");
    }
    let s = &r.stats;
    println!("N = {}, P = {}: Rp {:.2}  Rwp {:.2}  Rexp {:.2}  chi2 {:.4}  DW {:.4}", s.n, s.p, s.rp, s.rwp, s.rexp, s.chi2, s.durbin_watson);
    println!("  Bérar factor {:.4}: the esd's times it estimate the probable errors (local correlations; Bérar & Lelann 1991)", r.berar_factor);
    let fractions = drs::qpa::weight_fractions_esd(&model, set, &r.covariance);
    for (k, (rb, rf)) in s.phases.iter().enumerate() {
        let w = fractions[k].filter(|_| fractions.iter().flatten().count() > 1).map_or(String::new(), |(w, e)| format!("  weight fraction {} %", drs::cif_write::with_su(100.0 * w, 100.0 * e, 2)));
        println!("  phase {} ({}): R_Bragg {rb:.2}  R_F {rf:.2}{w}", k + 1, phases[k]);
    }
    // The amorphous content when a phase is an internal standard (P33; Orlhac et al. 2001, eqs. 3–4).
    match drs::qpa::amorphous_content(&model, set, &r.covariance) {
        Ok(Some(a)) => {
            let pc = |v: (f64, f64)| drs::cif_write::with_su(100.0 * v.0, 100.0 * v.1, 2);
            println!(
                "  amorphous content {} % of the sample, {} % of the mixture with {} % of phase {} ({}) as the internal standard",
                pc(a.amorphous_in_sample),
                pc(a.amorphous_in_mixture),
                100.0 * a.standard_fraction,
                a.standard + 1,
                phases[a.standard]
            );
            if a.amorphous_in_sample.0 < 0.0 {
                println!(
                    "    negative: the standard's weighed fraction {} % exceeds the {:.2} % the refinement gives it (a fully crystalline mixture); check the weighing, the standard's ATZ and the model",
                    100.0 * a.standard_fraction,
                    100.0 * drs::qpa::weight_fractions(&model)[a.standard].unwrap_or(0.0)
                );
            }
            for (k, w) in a.in_sample.iter().enumerate() {
                if let Some(w) = w {
                    println!("    phase {} ({}): {} % of the sample", k + 1, phases[k], pc(*w));
                }
            }
        }
        Ok(None) => {}
        Err(e) => println!("  amorphous content not calculated: {e}"),
    }
    println!("\n{:<45} {:>16} {:>12}", "parameter", "value", "esd");
    for (p, fp) in r.params.iter().zip(&set.free) {
        // A value imported outside its limits may only move towards them (it is not moved onto them).
        let (lo, hi) = set.limits_of(fp.fields[0].0, &model);
        let outside = lo.is_some_and(|b| p.value < b) || hi.is_some_and(|b| p.value > b);
        let note = if outside {
            format!("  outside its limits ({}, {})", lo.map_or("−∞".into(), |b| b.to_string()), hi.map_or("∞".into(), |b| b.to_string()))
        } else if p.at_limit {
            "  at its limit".into()
        } else {
            String::new()
        };
        println!("{:<45} {:16.6} {:12.6}{note}", p.name, p.value, p.esd);
    }
    // Symmetry modes (R29): each irrep's amplitude √(Σ Q²), its esd with the amplitudes' covariance, its polarisation.
    for (p, ph) in model.phases.iter().enumerate() {
        let Some(ms) = &ph.modes else { continue };
        println!("\nphase {} by symmetry modes: irrep amplitudes (Å) and polarisation vectors", p + 1);
        for ir in ms.irreps() {
            let e = ms.irrep_esd(&ir, |a, b| set.field_covariance(&r.covariance, drs::jacobian::Field::Mode(p, a), drs::jacobian::Field::Mode(p, b)));
            let pol: Vec<String> = ir.modes.iter().zip(&ir.polarisation).map(|(&k, c)| format!("{} {c:.4}", ms.modes[k].name)).collect();
            println!("  {:8} {:10.6} {:10.6}   ({})", ir.name, ir.amplitude, e, pol.join(", "));
        }
    }
    if !r.restraints.is_empty() {
        println!("\n{:<36} {:>12} {:>12} {:>10} {:>12} {:>10}", "restraint", "target", "value", "Δ/σ", "weight", "to χ²");
        for q in &r.restraints {
            println!("{:<36} {:12.5} {:12.5} {:10.3} {:12.4e} {:10.4}", q.label, q.target, q.value, (q.target - q.value) / q.sigma, q.weight, q.contribution);
        }
    }
    let mut pairs = Vec::new();
    for j in 0..r.params.len() {
        for k in j + 1..r.params.len() {
            if r.correlations[j][k].abs() > 0.8 {
                pairs.push((r.correlations[j][k], j, k));
            }
        }
    }
    pairs.sort_by(|a, b| b.0.abs().total_cmp(&a.0.abs()));
    println!("\nCorrelations above 0.8: {}", if pairs.is_empty() { "none" } else { "" });
    for (c, j, k) in pairs {
        println!("  {:+.3}  {}  {}", c, r.params[j].name, r.params[k].name);
    }
    // E1: the run's record.
    if let Some(f) = &trace {
        std::fs::write(f, drs::results::trace_json(&trajectory, &r.reason)).map_err(|e| format!("{f}: {e}"))?;
        println!("\nthe run's record ({} cycles) written in {f}", trajectory.cycles.len());
    }
    // E2: χ² on a grid of two parameters about the values reached.
    if let Some((na, nb)) = &map {
        let index = |n: &str| set.free.iter().position(|p| p.name == *n).ok_or(format!("--map: {n} is not a refined parameter (they are: {})", set.free.iter().map(|p| p.name.as_str()).collect::<Vec<_>>().join(", ")));
        let (a, b) = (index(na)?, index(nb)?);
        let half = |k: usize| -> Result<f64, String> {
            match span.strip_suffix("esd") {
                Some(x) => Ok(x.trim().parse::<f64>().map_err(|_| format!("--span {span}"))? * r.params[k].esd),
                None => span.trim().parse::<f64>().map_err(|_| format!("--span {span}: a number of esd's (\"3esd\") or a value")),
            }
        };
        let (sa, sb) = (half(a)?, half(b)?);
        if sa <= 0.0 || sb <= 0.0 {
            return Err(format!("--map: a span of 0 ({na} ± {sa}, {nb} ± {sb}): give --span as a value"));
        }
        let g = drs::lsq::chi2_grid(&model, obs, set, (a, b), (sa, sb), grid, |done, total| {
            eprint!("\r{done}/{total} points");
            !INTERRUPTED.load(std::sync::atomic::Ordering::SeqCst)
        });
        eprintln!();
        println!("\nchi2 on a {grid} x {grid} grid: {na} {} ± {sa}, {nb} {} ± {sb}", r.params[a].value, r.params[b].value);
        println!("{:>16} {:>16} {:>14}", na, nb, "chi2");
        for (i, x) in g.xa.iter().enumerate() {
            for (j, y) in g.xb.iter().enumerate() {
                println!("{x:16.8} {y:16.8} {:14.6}", g.chi2[i][j]);
            }
        }
        // The conditional cut's contours at the 1σ and 2σ ellipses' levels (E21).
        for c in g.contours_at_sigma(&[1.0, 2.0]) {
            let pieces = 1 + c.x.iter().filter(|v| v.is_nan()).count();
            println!("\ncontour chi2 = {:.6}: {} points in {pieces} piece(s)", c.level, c.x.iter().filter(|v| v.is_finite()).count());
            for (x, y) in c.x.iter().zip(&c.y).filter(|p| p.0.is_finite()) {
                println!("{x:16.8} {y:16.8}");
            }
        }
    }
    Ok(())
}

fn cmd_data(args: &[String]) -> Result<(), String> {
    use drs::formats::{self, table, Axis, Probe};
    use std::io::Write;
    if args.is_empty() || (wants_help(args) && !matches!(args[0].as_str(), "op" | "peaks" | "background")) {
        print!("{USAGE}");
        return Ok(());
    }
    if wants_help(args) {
        // `drs data peaks --help` and the others: the settings in the dictionary's words (app request 16, L088).
        let ids: Vec<&str> = match args[0].as_str() {
            "peaks" => vec!["peak_search"],
            "background" => vec!["auto_background"],
            _ => args.iter().map(String::as_str).filter(|a| drs::data_ops::settings(a).is_some()).collect(),
        };
        let usage = USAGE.lines().filter(|l| l.trim_start().starts_with(&format!("drs data {} ", args[0]))).collect::<Vec<_>>().join("\n");
        println!("{usage}\n");
        return print_settings(if ids.is_empty() { None } else { Some(&ids) });
    }
    match args[0].as_str() {
        "formats" => {
            println!("{:<18} {:>4}  {:<64} {:<26} source", "id", "Ins", "name", "extensions");
            for f in formats::FORMATS {
                let ins = f.ins.map_or(String::new(), |i| i.to_string());
                let status = f.not_read.map_or(String::new(), |w| format!(" [recognised, not read: {w}]"));
                println!("{:<18} {:>4}  {:<64} {:<26} {}{status}", f.id, ins, f.name, f.extensions.join(" "), f.source);
                // Its options (the app's request 51).
                for o in f.options {
                    use formats::options::OptionKind;
                    let kind = match o.kind {
                        OptionKind::Integer { min } => format!("integer from {min}"),
                        OptionKind::Choice(w) => w.join(" | "),
                        OptionKind::Number { unit } => format!("number ({unit})"),
                        OptionKind::Column { none } => format!("column from 1{}", if none { ", 0 none" } else { "" }),
                    };
                    println!("{:<24}{} = {kind}; default {}: {}", "", o.id, if o.default.is_empty() { "from the file" } else { o.default }, o.help);
                }
            }
            Ok(())
        }
        "identify" => {
            let file = args.get(1).ok_or("data identify: a file is needed")?;
            let bytes = std::fs::read(file).map_err(|e| format!("{file}: {e}"))?;
            let p = Probe::new(file, &bytes);
            let prefer = prefer_list(&args[2..])?;
            let prefer: Vec<&str> = prefer.iter().map(String::as_str).collect();
            let id = formats::identify_with(&p, &prefer);
            if id.candidates.is_empty() {
                println!("{file}: no reader recognises this file");
                return Ok(());
            }
            match id.format() {
                Some(c) if id.by_preference() => println!("{file}: {} — {} (tied with others: chosen by the order of preference)", c.id, c.name),
                Some(c) => println!("{file}: {} — {}", c.id, c.name),
                None => println!("{file}: ambiguous between {}", id.ambiguous().iter().map(|c| c.id).collect::<Vec<_>>().join(", ")),
            }
            println!("{:<18} {:>7} {:>9} {:>6}  name", "candidate", "content", "extension", "score");
            for c in &id.candidates {
                let reads = formats::by_id(c.id).is_some_and(|f| f.reads());
                println!("{:<18} {:>7.2} {:>9} {:>6.2}  {}{}", c.id, c.content, if c.extension { "yes" } else { "no" }, c.score, c.name, if reads { "" } else { " (recognised, not read)" });
            }
            Ok(())
        }
        "read" => {
            let file = args.get(1).ok_or("data read: a file is needed")?;
            let (mut fmt, mut scan, mut out): (Option<String>, Option<usize>, Option<String>) = (None, None, None);
            let mut o = table::TableOptions::default();
            let mut tabled = false;
            let (mut opts, mut prefer) = (formats::options::Options::new(), Vec::<String>::new());
            let mut i = 2;
            let col = |s: &str| s.parse::<usize>().ok().filter(|&c| c >= 1).map(|c| c - 1).ok_or(format!("column '{s}': a number from 1"));
            while i < args.len() {
                let v = args.get(i + 1).ok_or(format!("{} needs a value", args[i]))?;
                match args[i].as_str() {
                    "--as" => fmt = Some(v.clone()),
                    "--option" => {
                        let (k, x) = v.split_once('=').ok_or(format!("--option {v}: key=value"))?;
                        let x = x.trim();
                        let value = match (x.parse::<i64>(), x.parse::<f64>()) {
                            (Ok(n), _) => formats::options::OptionValue::Integer(n),
                            (_, Ok(f)) => formats::options::OptionValue::Number(f),
                            _ => formats::options::OptionValue::Text(x.to_string()),
                        };
                        opts.insert(k.trim().to_string(), value);
                    }
                    "--prefer" => prefer = v.split(',').map(|w| w.trim().to_string()).collect(),
                    "--scan" => scan = Some(v.parse::<usize>().map_err(|_| format!("--scan {v}"))?),
                    "--out" => out = Some(v.clone()),
                    "--delimiter" => o.delimiter = Some(if v == "tab" || v == "\\t" { '\t' } else if v == "space" { ' ' } else { v.chars().next().ok_or("--delimiter")? }),
                    "--decimal" => o.decimal = v.chars().next(),
                    "--header" => o.header_lines = Some(v.parse().map_err(|_| format!("--header {v}"))?),
                    "--x" => o.x = Some(col(v)?),
                    "--y" => o.y = Some(col(v)?),
                    "--sigma" => o.sigma = Some(if v == "none" { None } else { Some(col(v)?) }),
                    "--axis" => {
                        o.axis = Some(match v.to_lowercase().as_str() {
                            "2theta" | "tth" => Axis::TwoTheta,
                            "d" => Axis::D,
                            "q" => Axis::Q,
                            "tof" => Axis::TimeOfFlight,
                            _ => return Err(format!("--axis {v}: 2theta, d, q or tof")),
                        })
                    }
                    a => return Err(format!("data read: unknown option {a}")),
                }
                if !matches!(args[i].as_str(), "--as" | "--scan" | "--out" | "--option" | "--prefer") {
                    tabled = true;
                }
                i += 2;
            }
            let bytes = std::fs::read(file).map_err(|e| format!("{file}: {e}"))?;
            let d = if tabled {
                if fmt.as_deref().is_some_and(|f| f != "table") {
                    return Err("the table options read the file as a table (--as table)".into());
                }
                let text: String = match std::str::from_utf8(&bytes) {
                    Ok(t) => t.to_string(),
                    Err(_) => bytes.iter().map(|&b| b as char).collect(),
                };
                table::read_with(&text, &o).map_err(|e| format!("{file}: {e}"))?.0
            } else {
                let prefer: Vec<&str> = prefer.iter().map(String::as_str).collect();
                formats::read_with(&Probe::new(file, &bytes), fmt.as_deref(), &opts, &prefer).map_err(|e| format!("{file}: {e}"))?
            };
            // A file of a format the core writes (app request 12): the scan asked for, or the first.
            if let Some((f, id)) = out.as_deref().and_then(|f| formats::format_of(f).map(|id| (f, id))) {
                let k = scan.unwrap_or(1);
                let s = d.scans.get(k - 1).ok_or(format!("{file}: no scan {k}"))?;
                std::fs::write(f, formats::write(s, id)?).map_err(|e| format!("{f}: {e}"))?;
                println!("{f}: scan {k} written as {id} ({} points)", s.pattern.len());
                return Ok(());
            }
            let mut w: Box<dyn Write> = match &out {
                Some(f) => Box::new(std::fs::File::create(f).map_err(|e| format!("{f}: {e}"))?),
                None => Box::new(std::io::stdout().lock()),
            };
            let werr = |e: std::io::Error| e.to_string();
            let name = formats::by_id(d.format).map_or("", |f| f.name);
            for (k, s) in d.scans.iter().enumerate() {
                if scan.is_some_and(|n| n != k + 1) {
                    continue;
                }
                writeln!(w, "# format: {} ({name})", d.format).map_err(werr)?;
                writeln!(w, "# scan {} of {}: {}, {} points, x = {}", k + 1, d.scans.len(), s.label, s.pattern.len(), s.axis.label()).map_err(werr)?;
                if let Some(l) = s.wavelength {
                    writeln!(w, "# wavelength: {}{}{}", l.lambda1, l.lambda2.map_or(String::new(), |v| format!(" {v}")), l.ratio.map_or(String::new(), |v| format!(" ratio {v}"))).map_err(werr)?;
                }
                let inf = &s.info;
                for (key, v) in [("title", &inf.title), ("instrument", &inf.instrument), ("date", &inf.date), ("sample", &inf.sample)] {
                    if let Some(v) = v {
                        writeln!(w, "# {key}: {v}").map_err(werr)?;
                    }
                }
                for (key, v) in [("temperature", inf.temperature), ("monitor", inf.monitor), ("counting time", inf.counting_time), ("step", inf.step)] {
                    if let Some(v) = v {
                        writeln!(w, "# {key}: {v}").map_err(werr)?;
                    }
                }
                for (key, v) in &inf.other {
                    writeln!(w, "# {key}: {v}").map_err(werr)?;
                }
                for j in 0..s.pattern.len() {
                    writeln!(w, "{:14.6} {:16.6} {:14.6}", s.pattern.x[j], s.pattern.y[j], s.pattern.sigma[j]).map_err(werr)?;
                }
                writeln!(w).map_err(werr)?;
            }
            Ok(())
        }
        "info" => {
            let file = args.get(1).ok_or("data info: a file is needed")?;
            let bytes = std::fs::read(file).map_err(|e| format!("{file}: {e}"))?;
            let p = Probe::new(file, &bytes);
            // With --in <project folder>: by the file's place first (the file dictionary; the app's request 31).
            let within = args.iter().position(|a| a == "--in").and_then(|i| args.get(i + 1));
            let desc = match within {
                Some(dir) => {
                    let rel = std::path::Path::new(file).strip_prefix(dir).map_err(|_| format!("{file} is not in {dir}"))?;
                    formats::project::describe_in(&rel.to_string_lossy().replace('\\', "/"), &p)
                }
                None => formats::project::describe_path(std::path::Path::new(file), &p),
            };
            match desc {
                formats::project::Description::Project(k) => {
                    println!("{file}: {} ({}; {:?}), {} bytes; source: {}", k.name, k.id, k.role, bytes.len(), k.source);
                    return Ok(());
                }
                formats::project::Description::Unknown => {
                    println!("{file}: unknown ({} bytes): no source says what it is", bytes.len());
                    return Ok(());
                }
                formats::project::Description::Data { .. } => {}
            }
            let id = formats::identify(&p);
            let Some(c) = id.format() else {
                return match id.candidates.is_empty() {
                    true => Err(format!("{file}: no reader recognises this file")),
                    false => Err(format!("{file}: ambiguous between {} (drs data identify)", id.ambiguous().iter().map(|c| c.id).collect::<Vec<_>>().join(", "))),
                };
            };
            println!("{file}: data, {} ({})", c.name, c.id);
            match formats::read(&p, Some(c.id)) {
                Ok(d) => {
                    for (k, s) in d.scans.iter().enumerate() {
                        let pt = &s.pattern;
                        let (lo, hi) = if !pt.is_empty() { (pt.x[0], pt.x[pt.len() - 1]) } else { (0.0, 0.0) };
                        println!("scan {}: {}, {} points, x = {} from {lo} to {hi}{}", k + 1, s.label, pt.len(), s.axis.label(), s.wavelength.map_or(String::new(), |w| format!(", λ = {} Å", w.lambda1)));
                    }
                }
                Err(e) => println!("not read: {e}"),
            }
            Ok(())
        }
        "op" | "apply" | "peaks" | "background" | "constant-time" => data_process(args),
        "settings" => print_settings(args.get(1).map(|a| vec![a.as_str()]).as_deref()),
        other => Err(format!("data: unknown subcommand '{other}' (formats, identify, info, read, op, apply, peaks, background, settings)")),
    }
}

/// The dictionary of the operations and analyses (app request 16, L088): every setting in the core's words; all of
/// them, or those of `wanted`.
fn print_settings(wanted: Option<&[&str]>) -> Result<(), String> {
    use drs::operation_table::OPERATIONS;
    let ops: Vec<_> = OPERATIONS.iter().filter(|o| wanted.is_none_or(|w| w.contains(&o.id))).collect();
    if ops.is_empty() {
        return Err(format!("data settings: {} is not an operation or analysis ({})", wanted.unwrap_or_default().join(", "), OPERATIONS.iter().map(|o| o.id).collect::<Vec<_>>().join(", ")));
    }
    for o in ops {
        println!("{} ({}, {}): {}", o.id, o.name, o.kind, o.meaning);
        println!("  σ: {}", o.sigma);
        println!("  definition: {}", o.source);
        for st in o.settings {
            println!("  {} — {} [{}]: {}", st.key, st.name, st.unit, st.meaning);
            println!("      range: {}; default: {}; proposed: {}", st.range, st.default, st.proposed);
            if !st.raise.is_empty() {
                println!("      raising it: {}; lowering it: {}", st.raise, st.lower);
            }
            if !st.source.is_empty() {
                println!("      source: {}", st.source);
            }
        }
        println!();
    }
    Ok(())
}

/// x, y, σ of a scan after a header, as `drs data read` writes them; or the file of a format the core writes
/// (`formats::write`, app request 12) when `format` names one.
fn write_scan(w: &mut dyn std::io::Write, s: &drs::formats::Scan, header: &[String], format: Option<&str>) -> Result<(), String> {
    let werr = |e: std::io::Error| e.to_string();
    if let Some(id) = format {
        return w.write_all(drs::formats::write(s, id)?.as_bytes()).map_err(werr);
    }
    for h in header {
        writeln!(w, "# {h}").map_err(werr)?;
    }
    writeln!(w, "# {} points, x = {}{}", s.pattern.len(), s.axis.label(), s.wavelength.map_or(String::new(), |l| format!(", wavelength {}", l.lambda1))).map_err(werr)?;
    for j in 0..s.pattern.len() {
        writeln!(w, "{:14.6} {:16.6} {:14.6}", s.pattern.x[j], s.pattern.y[j], s.pattern.sigma[j]).map_err(werr)?;
    }
    Ok(())
}

/// `drs data op | apply | peaks | background`.
fn data_process(args: &[String]) -> Result<(), String> {
    use drs::control::syntax;
    use drs::data_ops::{self, Operation};
    use std::io::Write;
    let sub = args[0].as_str();
    let file = args.get(1).ok_or(format!("data {sub}: a file is needed"))?.clone();
    let (mut reader, mut out, mut rest): (Option<String>, Option<String>, Vec<String>) = (None, None, Vec::new());
    let mut i = 2;
    while i < args.len() {
        match args[i].as_str() {
            "--as" | "--out" => {
                let v = args.get(i + 1).ok_or(format!("{} needs a value", args[i]))?.clone();
                if args[i] == "--as" { reader = Some(v) } else { out = Some(v) }
                i += 2;
            }
            a => {
                rest.push(a.to_string());
                i += 1;
            }
        }
    }
    let dir = std::path::Path::new(&file).parent().map(|d| d.to_path_buf()).unwrap_or_default();
    let read = |f: &str, r: Option<&str>| -> Result<drs::formats::Scan, String> {
        let d = drs::formats::read_path(std::path::Path::new(f), r).map_err(|e| format!("{f}: {e}"))?;
        d.scans.into_iter().next().ok_or(format!("{f}: no scan"))
    };
    let mut w: Box<dyn Write> = match &out {
        Some(f) => Box::new(std::io::BufWriter::new(std::fs::File::create(f).map_err(|e| format!("{f}: {e}"))?)),
        None => Box::new(std::io::BufWriter::new(std::io::stdout())),
    };
    let werr = |e: std::io::Error| e.to_string();
    let fmt = out.as_deref().and_then(drs::formats::format_of);
    match sub {
        "op" => {
            // The operation as the control file writes it: op name, then key=value pairs.
            let name = rest.first().ok_or("data op: name the operation (mask, delete, smooth, rebin, normalize, add, subtract, merge, convert, wavelength)")?;
            let mut t = format!("[[pattern.operation]]\nop = {}\n", syntax::string(name));
            for kv in &rest[1..] {
                let (k, v) = kv.split_once('=').ok_or(format!("'{kv}': key=value"))?;
                let v = v.trim();
                let quoted = v.parse::<f64>().is_err() && !v.starts_with('[') && !v.starts_with('"');
                t.push_str(&format!("{} = {}\n", k.trim(), if quoted { syntax::string(v) } else { v.to_string() }));
            }
            let doc = syntax::parse(&t);
            if let Some(d) = doc.diagnostics.iter().find(|d| !d.warning) {
                return Err(format!("data op: {}", d.message));
            }
            let op = Operation::from_table(doc.tables.iter().find(|t| t.array).ok_or("data op: no operation")?).map_err(|d| format!("data op: {}", d.message))?;
            // Files named by the operation are relative to the data file's folder, as in a control file.
            let load = |f: &str, r: Option<&str>| read(&dir.join(f).to_string_lossy(), r);
            let scan = data_ops::apply_all(&read(&file, reader.as_deref())?, std::slice::from_ref(&op), &load)?;
            let pairs: Vec<String> = op.to_pairs().iter().map(|(k, v)| format!("{k} = {v}")).collect();
            write_scan(&mut w, &scan, &[format!("{file}, operation {{ {} }}", pairs.join(", "))], fmt)
        }
        "constant-time" => {
            let t = opt_value(args, "--time").ok_or("data constant-time: --time <seconds per point>")?;
            let t = num(&t, "--time")?;
            let d = drs::formats::read_path(std::path::Path::new(&file), reader.as_deref()).map_err(|e| format!("{file}: {e}"))?;
            let scan = data_ops::constant_time(&d.scans, t)?;
            write_scan(&mut w, &scan, &[format!("{file}: {} scans brought to {t} s per point (PL6)", d.scans.len())], fmt)
        }
        "apply" => {
            let o = open_control(std::path::Path::new(&file))?;
            let d = o.control.data.as_ref().ok_or(format!("{file}: no data file in [pattern]"))?;
            let scan = o.scan.as_ref().ok_or(format!("{file}: no data"))?;
            let ops: Vec<String> = d.operations.iter().map(|op| op.name().to_string()).collect();
            write_scan(&mut w, scan, &[format!("{} after {} operation(s): {}", d.file, ops.len(), ops.join(", "))], fmt)
        }
        "peaks" => {
            let (mut prom, mut sig) = (data_ops::DEFAULT_MIN_PROMINENCE, data_ops::DEFAULT_MIN_SIGNIFICANCE);
            let mut k = 0;
            while k < rest.len() {
                let v = rest.get(k + 1).ok_or(format!("{} needs a value", rest[k]))?;
                match rest[k].as_str() {
                    "--prominence" => prom = num(v, "--prominence")?,
                    "--significance" => sig = num(v, "--significance")?,
                    a => return Err(format!("data peaks: unknown option {a}")),
                }
                k += 2;
            }
            let scan = read(&file, reader.as_deref())?;
            let peaks = data_ops::peak_search(&scan.pattern, prom, sig);
            writeln!(w, "# {file}: {} peaks (prominence ≥ {prom} and ≥ {sig} σ)", peaks.len()).map_err(werr)?;
            writeln!(w, "# {:>12} {:>14} {:>14} {:>10}", "x", "height", "prominence", "fwhm").map_err(werr)?;
            for p in peaks {
                writeln!(w, "{:14.5} {:14.3} {:14.3} {:10.5}", p.x, p.height, p.prominence, p.fwhm).map_err(werr)?;
            }
            Ok(())
        }
        "background" => {
            let mut log_lambda = None;
            let mut curve = false;
            let mut k = 0;
            while k < rest.len() {
                match rest[k].as_str() {
                    "--log-lambda" => {
                        log_lambda = Some(num(rest.get(k + 1).ok_or("--log-lambda needs a value")?, "--log-lambda")?);
                        k += 2;
                    }
                    "--curve" => {
                        curve = true;
                        k += 1;
                    }
                    a => return Err(format!("data background: unknown option {a}")),
                }
            }
            let scan = read(&file, reader.as_deref())?;
            let b = data_ops::auto_background(&scan.pattern, log_lambda, data_ops::DEFAULT_MAX_ITER)?;
            writeln!(w, "# {file}: arPLS background, log10 λ = {}, {} iteration(s)", b.log_lambda, b.iterations).map_err(werr)?;
            if curve {
                for (x, y) in scan.pattern.x.iter().zip(&b.curve) {
                    writeln!(w, "{x:14.6} {y:16.6}").map_err(werr)?;
                }
            } else {
                for (x, y) in &b.points {
                    writeln!(w, "{x:14.6} {y:16.6}").map_err(werr)?;
                }
            }
            Ok(())
        }
        _ => unreachable!(),
    }
}

/// A .pcr file read as `drs calc` reads it: the .pcr, its absolute folder and its stem.
fn read_pcr_file(path: &std::path::Path) -> Result<(drs::pcr::Pcr, std::path::PathBuf, String), String> {
    let bytes = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let text: String = String::from_utf8(bytes.clone()).unwrap_or_else(|_| bytes.iter().map(|&b| b as char).collect());
    let pcr = drs::pcr::read(&text).map_err(|e| format!("{}: {e}", path.display()))?;
    let dir = std::fs::canonicalize(path).ok().and_then(|p| p.parent().map(|d| d.to_path_buf())).unwrap_or_else(|| std::path::PathBuf::from("."));
    let stem = path.file_stem().and_then(|s| s.to_str()).unwrap_or("").to_string();
    Ok((pcr, dir, stem))
}

/// The folder of a file about to be written (created if needed).
fn out_dir(out: &std::path::Path) -> Result<std::path::PathBuf, String> {
    let d = out.parent().filter(|d| !d.as_os_str().is_empty()).map(|d| d.to_path_buf()).unwrap_or_else(|| std::path::PathBuf::from("."));
    std::fs::create_dir_all(&d).map_err(|e| format!("{}: {e}", d.display()))?;
    Ok(d)
}

fn cmd_rigid(args: &[String]) -> Result<(), String> {
    use drs::rigid_body::{groups, positions, satellite_axis};
    if wants_help(args) || args.len() != 1 {
        print!("{USAGE}");
        return Ok(());
    }
    let text: String = std::fs::read(&args[0]).map_err(|e| format!("{}: {e}", args[0]))?.iter().map(|&b| b as char).collect();
    let p = drs::pcr::read(&text).map_err(|e| format!("{}: {e:?}", args[0]))?;
    let mut any = false;
    for (k, ph) in p.phases.iter().enumerate().filter(|(_, ph)| ph.flag("Jbt").abs() == 4) {
        any = true;
        let c = ph.cell;
        let cell = drs::cell::Cell::new(c[0], c[1], c[2], c[3], c[4], c[5]).map_err(|e| format!("phase {}: the cell: {e:?}", k + 1))?;
        let m = cell.crystal_to_orthonormal();
        let mut gs = groups(ph).map_err(|e| format!("phase {}: {e}", k + 1))?;
        let x = positions(ph).map_err(|e| format!("phase {}: {e}", k + 1))?;
        let species: Vec<String> = ph.atoms.iter().map(|a| a.species.clone()).collect();
        drs::rigid_body::read_tls(ph, &mut gs, &species, &x, &cell).map_err(|e| format!("phase {}: {e}", k + 1))?;
        let betas = drs::rigid_body::betas(&gs, &cell, &x);
        println!("phase {} {}: {} groups", k + 1, ph.name.trim(), gs.len());
        for (gi, g) in gs.iter().enumerate() {
            let unit = if g.degrees { "deg" } else { "rad" };
            println!("  group {} ({:?}{}, P6 = {}): origin {:.5} {:.5} {:.5}, θ {} φ {} χ {} {unit}", gi + 1, g.coordinates, if g.tls() { ", TLS" } else { "" }, g.option, g.origin[0], g.origin[1], g.origin[2], g.theta, g.phi, g.chi);
            if let (Some(s), Some((o, t, f))) = (g.satellite, satellite_axis(&gs, gi, &cell)) {
                println!("    satellite of group {} at its atom {} (from {}): origin {:.5} {:.5} {:.5}, θ {t:.5} φ {f:.5} {unit}", s.parent + 1, s.atom, s.from, o[0], o[1], o[2]);
            }
            for mm in &g.members {
                let (a, w) = (x[mm.index], ph.atoms[mm.index].x);
                let d = drs::mat3::mul_vec(&m, &[a[0] - w[0], a[1] - w[1], a[2] - w[2]]);
                println!("    {:<6} {:9.5} {:9.5} {:9.5}   file {:9.5} {:9.5} {:9.5}   {:.1e} Å", mm.label, a[0], a[1], a[2], w[0], w[1], w[2], (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt());
            }
            // TLS (L367): T, L, S in the units of FullProf's .out (Å², deg², Å·deg), its centre, and the β they give.
            if let (Some(t), Some(c)) = (&g.displacements, drs::rigid_body::tls_centre(g)) {
                let (d, f) = (180.0 / std::f64::consts::PI, |v: &[f64], s: f64| v.iter().map(|x| format!("{:.6}", x * s)).collect::<Vec<_>>().join(" "));
                println!("    TLS centre {:.5} {:.5} {:.5}", c[0], c[1], c[2]);
                println!("    T (Å²)   11 22 33 12 13 23: {}", f(&t.t, 1.0));
                println!("    L (deg²) 11 22 33 12 13 23: {}", f(&t.l, d * d));
                println!("    S (Å·deg) 11 22 33 12 13 23 21 31 32: {}", f(&t.s, d));
                for (a, b) in betas.iter().filter(|(a, _)| drs::rigid_body::tls_members(&gs, gi).iter().any(|&(h, mm)| gs[h].members[mm].index == *a)) {
                    println!("    β {:<6} {}", ph.atoms[*a].label, b.iter().map(|v| format!("{v:.5}")).collect::<Vec<_>>().join(" "));
                }
            }
        }
    }
    if !any {
        return Err(format!("{}: no rigid-body phase (Jbt = 4)", args[0]));
    }
    Ok(())
}

fn cmd_texture(args: &[String]) -> Result<(), String> {
    if wants_help(args) || args.len() != 1 {
        print!("{USAGE}");
        return Ok(());
    }
    let m = load_model(&args[0], None, false)?.0.model;
    for (k, ph) in m.phases.iter().enumerate() {
        let Some(t) = drs::texture::of_phase(ph, m.march) else {
            println!("phase {} {}: exponential, G1 {} G2 {} along [{} {} {}]", k + 1, ph.name.trim(), ph.pref[0], ph.pref[1], ph.pref_direction[0], ph.pref_direction[1], ph.pref_direction[2]);
            continue;
        };
        let how = if ph.texture.is_none() { "March–Dollase" } else if t.steps > 0 { "March–Dollase axes, an anvil cell" } else { "March–Dollase axes" };
        println!("phase {} {}: {how}{}", k + 1, ph.name.trim(), if t.steps > 0 { format!(", {} steps", t.steps) } else { String::new() });
        // D against ρ every 15° (`texture::march_dollase_curve`, Publish's P59).
        let angles: Vec<f64> = drs::texture::march_dollase_curve(1.0, 6).iter().map(|p| p.0).collect();
        println!("  axis  direction           r         fraction   D at {}°", angles.iter().map(|a| format!("{a:.0}")).collect::<Vec<_>>().join(", "));
        for (j, a) in t.axes.iter().enumerate() {
            let d: Vec<String> = drs::texture::march_dollase_curve(a.r, 6).iter().map(|p| format!("{:.4}", p.1)).collect();
            println!("  {:<4}  [{} {} {}]  {:>9.5}  {:>9.5}   {}", j + 1, a.direction[0], a.direction[1], a.direction[2], a.r, a.fraction, d.join(" "));
        }
    }
    Ok(())
}

fn cmd_move(args: &[String]) -> Result<(), String> {
    if wants_help(args) || args.len() != 2 {
        print!("{USAGE}");
        return Ok(());
    }
    let (file, folder) = (std::path::Path::new(&args[0]), std::path::Path::new(&args[1]));
    let from = file.parent().map(|p| if p.as_os_str().is_empty() { std::path::Path::new(".") } else { p }).unwrap_or(std::path::Path::new("."));
    let name = file.file_name().ok_or(format!("{}: no file name", file.display()))?;
    let target = folder.join(name);
    if target.exists() {
        return Err(format!("{} is there already", target.display()));
    }
    let text = std::fs::read_to_string(file).map_err(|e| format!("{}: {e}", file.display()))?;
    let moved = drs::control::relocate(&text, from, folder)?;
    std::fs::write(&target, moved).map_err(|e| format!("{}: {e}", target.display()))?;
    std::fs::remove_file(file).map_err(|e| format!("{}: {e}", file.display()))?;
    println!("{} moved to {}", file.display(), target.display());
    Ok(())
}

fn cmd_convert(args: &[String]) -> Result<(), String> {
    use drs::control::{self, Control};
    if wants_help(args) || args.is_empty() {
        print!("{USAGE}");
        return Ok(());
    }
    let (mut input, mut output) = (None, None);
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--to" => {
                output = Some(args.get(i + 1).ok_or("--to needs a file")?.clone());
                i += 2;
            }
            a => {
                input = Some(a.to_string());
                i += 1;
            }
        }
    }
    let input = input.ok_or("convert: give the file to convert")?;
    // Without --to, the file dictionary's name (R27, L091): the control file beside the .pcr, named after it; a
    // .pcr exported from a control file in the folder Exports/FullProf/ next to it (L123).
    let output = match output {
        Some(o) => o,
        None => {
            let p = std::path::Path::new(&input);
            let stem = p.file_stem().and_then(|s| s.to_str()).ok_or("convert: a file name")?;
            let dir = p.parent().unwrap_or(std::path::Path::new(""));
            let (kind, key) = if input.to_lowercase().ends_with(".pcr") { ("control", "project") } else { ("fullprof_pcr_export", "control") };
            dir.join(drs::formats::files::path(kind, &[(key, stem)])?).to_string_lossy().into_owned()
        }
    };
    let (inp, outp) = (std::path::Path::new(&input), std::path::Path::new(&output));
    let ext = |p: &std::path::Path| p.extension().map(|e| e.to_string_lossy().to_lowercase()).unwrap_or_default();
    match (ext(inp).as_str(), ext(outp).as_str()) {
        ("pcr", "drc") => {
            let (pcr, dir, stem) = read_pcr_file(inp)?;
            let base = out_dir(outp)?;
            let c = Control::from_pcr(&pcr, &dir, &stem, &base)?;
            let text = c.to_text();
            std::fs::write(outp, &text).map_err(|e| format!("{output}: {e}"))?;
            let free = c.refinement.build(&c.model).map(|s| s.free.len()).unwrap_or(0);
            println!("{output}: {} phase(s), {} refined parameter(s), data {}", c.model.phases.len(), free, c.data.as_ref().map_or("none".to_string(), |d| d.file.clone()));
            Ok(())
        }
        ("drc", "pcr") => {
            let o = open_control(inp)?;
            let base = out_dir(outp)?;
            let mut c = o.control;
            // The data file named from the .pcr's folder.
            if let (Some(d), Some(p)) = (c.data.as_mut(), o.data_path.as_ref()) {
                d.file = control::relative(p, &base);
            }
            let stem = outp.file_stem().and_then(|s| s.to_str()).ok_or("--to: a file name")?;
            let e = control::pcr_export::export(&c, stem, o.data.as_ref())?;
            std::fs::write(outp, &e.pcr).map_err(|er| format!("{output}: {er}"))?;
            println!("{output}");
            for (f, t) in &e.files {
                let p = base.join(f);
                std::fs::write(&p, t).map_err(|er| format!("{}: {er}", p.display()))?;
                println!("{}", p.display());
            }
            for n in &e.notes {
                println!("note: {n}");
            }
            Ok(())
        }
        _ => Err("convert: from .pcr to .drc, or from .drc to .pcr".into()),
    }
}

fn cmd_check(args: &[String]) -> Result<(), String> {
    use drs::check::{self, Options, Severity};
    if wants_help(args) || args.is_empty() {
        print!("{USAGE}");
        return Ok(());
    }
    let (mut files, mut json, mut opts) = (Vec::new(), false, Options::default());
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--json" => json = true,
            "--no-run" => opts.run = false,
            "--correlation" => {
                opts.correlation = num(args.get(i + 1).ok_or("--correlation needs a value")?, "--correlation")?;
                i += 1;
            }
            a => files.push(a.to_string()),
        }
        i += 1;
    }
    let file = files.first().cloned().ok_or("check: give a .drc or .pcr file")?;
    let data = files.get(1).map(std::path::Path::new);
    if files.len() > 2 || (data.is_some() && !file.to_lowercase().ends_with(".pcr")) {
        return Err("check: one file, and a .pcr's data file (a control file names its own)".into());
    }
    let diags = check::check_path_with_data(std::path::Path::new(&file), data, &opts);
    let count = |s: Severity| diags.iter().filter(|d| d.severity == s).count();
    if json {
        println!("[{}]", diags.iter().map(|d| d.to_json(&file)).collect::<Vec<_>>().join(",\n "));
    } else {
        for d in &diags {
            println!("{}", d.to_text(&file));
        }
        println!("{file}: {} error(s), {} warning(s), {} hint(s)", count(Severity::Error), count(Severity::Warning), count(Severity::Hint));
    }
    if check::has_errors(&diags) {
        return Err(format!("{file}: {} error(s)", count(Severity::Error)));
    }
    Ok(())
}

/// Background points from a points file (ours or a .bac) or a control file's own points.
/// Points and the interpolation the file states.
type Points = (Vec<(f64, f64)>, Option<drs::background::Interpolation>);

fn points_of(file: &str) -> Result<Points, String> {
    use drs::powder::BackgroundModel;
    let path = std::path::Path::new(file);
    if drs::control::is_control_file(path) {
        let o = open_control(path)?;
        return match o.control.model.background {
            BackgroundModel::Points(p) => Ok((p, Some(drs::background::Interpolation::Linear))),
            BackgroundModel::Spline(p) => Ok((p, Some(drs::background::Interpolation::Spline))),
            _ => Err(format!("{file}: its background is a function, not points")),
        };
    }
    let pf = drs::background::read_path(path)?;
    Ok((pf.points, pf.interpolation))
}

fn cmd_background(args: &[String]) -> Result<(), String> {
    use drs::background::{self, Function, Interpolation};
    use drs::powder::BackgroundModel;
    if wants_help(args) || args.len() < 2 {
        print!("{USAGE}");
        return Ok(());
    }
    let (sub, file) = (args[0].as_str(), args[1].clone());
    let mut opt = std::collections::BTreeMap::<String, String>::new();
    let mut flags = std::collections::BTreeSet::<String>::new();
    let mut i = 2;
    while i < args.len() {
        let a = args[i].clone();
        if matches!(a.as_str(), "--spline" | "--inline") {
            flags.insert(a);
            i += 1;
        } else {
            opt.insert(a.clone(), args.get(i + 1).ok_or(format!("{a} needs a value"))?.clone());
            i += 2;
        }
    }
    let how_of = |stated: Option<Interpolation>| if flags.contains("--spline") { Interpolation::Spline } else { stated.unwrap_or(Interpolation::Linear) };
    match sub {
        "read" => {
            let (p, how) = points_of(&file)?;
            println!("# {file}: {} points, {}", p.len(), how.map_or("interpolation not stated (FullProf takes it from Nba)", |h| h.name()));
            for (x, y) in p {
                println!("{x:12.4} {y:14.4}");
            }
            Ok(())
        }
        "convert" => {
            let to = opt.get("--to").ok_or("convert: --to <file>")?;
            let (p, how) = points_of(&file)?;
            let x = match opt.get("--data") {
                Some(d) => Some(drs::formats::read_path(std::path::Path::new(d), None).map_err(|e| format!("{d}: {e}"))?.scans.remove(0).pattern.x),
                None => None,
            };
            background::write_path(std::path::Path::new(to), &p, how_of(how), x.as_deref())?;
            println!("{to}: {} points", p.len());
            Ok(())
        }
        "eval" => {
            let data = opt.get("--data").ok_or("eval: --data <file>")?;
            let (p, how) = points_of(&file)?;
            let scan = drs::formats::read_path(std::path::Path::new(data), None).map_err(|e| format!("{data}: {e}"))?.scans.remove(0);
            let b = background::interpolate(&p, how_of(how), &scan.pattern.x);
            let mut out = format!("# {file} at the points of {data}: {}\n", how_of(how).name());
            for (x, v) in scan.pattern.x.iter().zip(&b) {
                out.push_str(&format!("{x:12.4} {v:14.4}\n"));
            }
            match opt.get("--out") {
                Some(o) => std::fs::write(o, out).map_err(|e| format!("{o}: {e}"))?,
                None => print!("{out}"),
            }
            Ok(())
        }
        "fit" => {
            let (p, _) = points_of(&file)?;
            let n = |k: &str| opt.get(k).map(|v| num(v, k)).transpose();
            let f = if let Some(t) = n("--chebyshev")? {
                Function::Chebyshev { from: n("--from")?.unwrap_or(p[0].0), to: n("--to")?.unwrap_or(p[p.len() - 1].0), terms: t as usize }
            } else if let Some(t) = n("--polynomial")? {
                Function::Polynomial { origin: n("--origin")?.ok_or("--polynomial needs --origin <2θ> (FullProf's Bkpos)")?, terms: t as usize }
            } else {
                return Err("fit: --chebyshev <n> or --polynomial <n> --origin <2θ>".into());
            };
            let (c, rms) = background::fit(&p, f)?;
            println!("# {file}: {f:?} fitted to {} points, rms of the residuals {rms:.4}", p.len());
            println!("{}", c.iter().map(|v| drs::control::syntax::float(*v)).collect::<Vec<_>>().join(" "));
            if let Some(drc) = opt.get("--into") {
                let text = std::fs::read_to_string(drc).map_err(|e| format!("{drc}: {e}"))?;
                let dir = std::path::Path::new(drc).parent().map(|d| d.to_path_buf()).unwrap_or_default();
                let mut ctl = drs::control::Control::from_text_in(&text, &dir).map_err(|d| format!("{drc}: {}", d.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("; ")))?;
                ctl.model.background = match f {
                    Function::Chebyshev { from, to, .. } => BackgroundModel::Chebyshev { from, to, coefficients: c },
                    Function::Polynomial { origin, .. } => BackgroundModel::Polynomial { origin, coefficients: c },
                };
                ctl.background_file = None;
                // The refinement's terms on the old points no longer apply.
                ctl.refinement.parameters.retain(|t| !t.iter().any(|x| matches!(x.0, drs::jacobian::Field::Background(_))));
                std::fs::write(drc, drs::control::update::update(&text, &ctl)).map_err(|e| format!("{drc}: {e}"))?;
                println!("{drc}: the background is now this function (its terms not refined yet)");
            }
            Ok(())
        }
        "store" => {
            let path = std::path::Path::new(&file);
            let dir = path.parent().map(|d| d.to_path_buf()).unwrap_or_default();
            let text = std::fs::read_to_string(path).map_err(|e| format!("{file}: {e}"))?;
            let mut ctl = drs::control::Control::from_text_in(&text, &dir).map_err(|d| format!("{file}: {}", d.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("; ")))?;
            if !matches!(ctl.model.background, BackgroundModel::Points(_) | BackgroundModel::Spline(_)) {
                return Err(format!("{file}: its background is a function, not points"));
            }
            if flags.contains("--inline") {
                ctl.background_file = None;
            } else {
                let to = opt.get("--file").ok_or("store: --file <name> or --inline")?;
                ctl.background_file = Some(to.clone());
                // A .bac holds the background at every data point: the control file's data are read for it.
                let x = open_control(path).ok().and_then(|o| o.data).map(|d| d.x);
                // A file that already holds these points is left as it is (the app's request 48).
                match ctl.background_file_to_write(&dir, x.as_deref()).unwrap()? {
                    Some((name, t)) => std::fs::write(dir.join(&name), t).map_err(|e| format!("{name}: {e}"))?,
                    None => println!("{to}: unchanged, left as it is"),
                }
            }
            // The canonical text: the points move between the file and the control file.
            std::fs::write(path, ctl.to_text()).map_err(|e| format!("{file}: {e}"))?;
            println!("{file}: background points {}", ctl.background_file.as_ref().map_or("written in the control file".to_string(), |f| format!("in {f}")));
            Ok(())
        }
        other => Err(format!("background: unknown subcommand '{other}' (read, convert, eval, fit, store)")),
    }
}

/// A structure for `drs bonds`: name, structure, standard deviations, the file's DMAX and TOL.
type BondsItem = (String, drs::structure_factor::Structure, Option<drs::bonds::Uncertainties>, (f64, f64), Option<f64>);

/// `drs bonds` on a `.pcr` or control file: the refined state's distances, angles and valence sums with su's (P4, P5).
fn bonds_refined(file: &str, opt: &std::collections::BTreeMap<String, f64>, charges: &[(String, i32)]) -> Result<(), String> {
    use drs::bonds;
    use drs::cif_write::with_su;
    let path = std::path::Path::new(file);
    let cycles = opt.get("--cycles").copied().unwrap_or(0.0) as usize;
    let (mut model, data, set) = if file.to_lowercase().ends_with(".pcr") {
        let text: String = std::fs::read(path).map_err(|e| format!("{file}: {e}"))?.iter().map(|&b| b as char).collect();
        let p = drs::pcr::read(&text).map_err(|e| format!("{file}: {e:?}"))?;
        let dir = path.parent().map(|d| d.to_path_buf()).unwrap_or_default();
        let stem = path.file_stem().map(|s| s.to_string_lossy().to_string()).unwrap_or_default();
        let l = drs::powder::Model::from_pcr(&p, &dir, &stem)?;
        let set = drs::refinable::ParameterSet::from_pcr(&p, &l.model)?;
        (l.model, l.data, set)
    } else {
        let o = open_control(path)?;
        let set = o.control.refinement.build(&o.control.model)?;
        (o.control.model, o.data, set)
    };
    // Without data the values are the file's and no su is propagated.
    let prop = match &data {
        Some(d) => {
            let r = drs::lsq::refine(&mut model, d, &set, drs::lsq::Options { max_cycles: cycles, eps: 0.1 }, |_| {});
            println!("# {file}: the state after {} cycles (Rwp {:.2} %); su's through the covariance", r.cycles, r.stats.rwp);
            Some(drs::geometry::Propagator::new(&model, &set, &r))
        }
        None => {
            println!("# {file}: no data, the file's values without su's");
            None
        }
    };
    for p in 0..model.phases.len() {
        let Some(st0) = drs::geometry::structure(&model, p) else { continue };
        let dmax = opt.get("--dmax").copied().unwrap_or(3.0);
        let dangl = opt.get("--angles").copied().unwrap_or(0.0);
        let tol = opt.get("--tol").copied().unwrap_or(20.0) / 100.0;
        let (mut st, contacts, angles) = drs::geometry::contacts_with_su(&model, p, prop.as_ref(), dmax, dangl)?;
        if !charges.is_empty() {
            st = bonds::with_charges(&st, charges);
        }
        let labels = drs::geometry::unique_labels(&st0);
        println!("# phase {} ({}): distances up to {dmax} Å", p + 1, drs::cif_write::phase_name(&model, p));
        for (k, list) in contacts.iter().enumerate() {
            println!("{}", labels[k]);
            for c in list {
                println!("  {:>6} {:>14}  {:8.4} {:8.4} {:8.4}  {}", labels[c.atom], with_su(c.distance, c.sigma, 4), c.position[0], c.position[1], c.position[2], st.group.ops[c.op].to_xyz());
            }
        }
        if dangl > 0.0 {
            println!("# angles between neighbours within {dangl} Å");
            for a in &angles {
                let c = &contacts[a.centre];
                println!("  {}-{}-{}  {}", labels[c[a.first].atom], labels[a.centre], labels[c[a.second].atom], with_su(a.angle.0, a.angle.1, 2));
            }
        }
        if st.atoms.iter().all(|a| bonds::atom_ion(a).is_some()) {
            let occ = bonds::bond_str_occupancies(&st)?;
            let (sums, gii) = bonds::valence_sums(&st, &contacts, &occ, tol)?;
            println!("# bond-valence sums (tol {:.0} %)", tol * 100.0);
            for v in &sums {
                println!("  {:6} {:>3}  valence {:+}  sum {}", labels[v.atom], v.coordination, v.valence, with_su(v.sum.0, v.sum.1, 3));
            }
            println!("  global instability index {:.4}", gii.a);
        } else {
            println!("# no valence sums: give the charges (--charges Y=3,O=-2, or charge = … in the control file)");
        }
    }
    Ok(())
}

fn cmd_bonds(args: &[String]) -> Result<(), String> {
    use drs::bonds;
    if wants_help(args) || args.is_empty() {
        print!("{USAGE}");
        return Ok(());
    }
    let file = &args[0];
    let mut opt = std::collections::BTreeMap::<String, f64>::new();
    let mut fullprof = false;
    let mut charges: Vec<(String, i32)> = Vec::new();
    let mut i = 1;
    while i < args.len() {
        if args[i] == "--fullprof" {
            fullprof = true;
            i += 1;
            continue;
        }
        if args[i] == "--charges" {
            // Publish P5: "Y=3,O=-2" (by element or label), for the valence sums.
            let v = args.get(i + 1).ok_or("--charges needs a list, e.g. Y=3,O=-2")?;
            for part in v.split(',') {
                let (k, q) = part.split_once('=').ok_or(format!("--charges: {part} is not label=charge"))?;
                charges.push((k.trim().to_string(), q.trim().parse::<i32>().map_err(|_| format!("--charges: {q} is not a whole number"))?));
            }
            i += 2;
            continue;
        }
        let v = args.get(i + 1).ok_or(format!("{} needs a value", args[i]))?;
        opt.insert(args[i].clone(), num(v, &args[i])?);
        i += 2;
    }
    let path = std::path::Path::new(file);
    let ext = path.extension().map(|e| e.to_string_lossy().to_lowercase()).unwrap_or_default();
    // A .pcr or control file: its refined state (Publish P4), --cycles n from the file's values (0: at them), every
    // su propagated through the covariance as the CIF's.
    if ext != "cfl" {
        return bonds_refined(file, &opt, &charges);
    }
    // (name, structure, standard deviations, the file's DMAX and TOL)
    let mut items: Vec<BondsItem> = Vec::new();
    let structures = |m: &drs::powder::Model| {
        m.phases.iter().filter_map(|ph| match &ph.intensities {
            drs::powder::Intensities::Structure(st, _) => Some((ph.name.clone(), (**st).clone(), None, (0.0, 0.0), None)),
            _ => None,
        }).collect::<Vec<_>>()
    };
    match ext.as_str() {
        "cfl" => {
            let c = drs::cfl::read(&std::fs::read_to_string(path).map_err(|e| format!("{file}: {e}"))?).map_err(|e| format!("{file}: {e}"))?;
            items.push((c.title.clone(), c.structure, Some(c.uncertainties), c.dmax, c.tol));
        }
        "pcr" => {
            let text: String = std::fs::read(path).map_err(|e| format!("{file}: {e}"))?.iter().map(|&b| b as char).collect();
            let p = drs::pcr::read(&text).map_err(|e| format!("{file}: {e:?}"))?;
            let dir = path.parent().map(|d| d.to_path_buf()).unwrap_or_default();
            let stem = path.file_stem().map(|s| s.to_string_lossy().to_string()).unwrap_or_default();
            items.extend(structures(&drs::powder::Model::from_pcr(&p, &dir, &stem)?.model));
        }
        _ => items.extend(structures(&open_control(path)?.control.model)),
    }
    if items.is_empty() {
        return Err(format!("{file}: no structure phase"));
    }
    for (name, st, esd, dm, tol) in items {
        let dmax = opt.get("--dmax").copied().unwrap_or(if dm.0 > 0.0 { dm.0 } else { 3.0 });
        let dangl = opt.get("--angles").copied().unwrap_or(dm.1);
        let tol = opt.get("--tol").copied().or(tol.filter(|t| *t > 0.1)).unwrap_or(20.0) / 100.0;
        let mut contacts = bonds::contacts(&st, esd.as_ref(), dmax);
        if fullprof {
            if let Some(u) = &esd {
                for (k, list) in contacts.iter_mut().enumerate() {
                    for c in list.iter_mut() {
                        c.sigma = bonds::bond_str_distance_sigma(&st, u, k, c);
                    }
                }
            }
        }
        // The one su rule (IUCr, rounded up: `cif_write::with_su`); --fullprof prints as Bond_Str does (fixed decimals,
        // the su rounded to the nearest unit of the last place), to compare with its .dis.
        let with_places = |v: f64, s: f64, places: usize| match (fullprof, s > 0.0) {
            (true, true) => format!("{v:.places$}({:.0})", s * 10f64.powi(places as i32)),
            (true, false) => format!("{v:.places$}"),
            (false, _) => drs::cif_write::with_su(v, s, places),
        };
        let with = |v: f64, s: f64| with_places(v, s, 4);
        println!("# {name}: distances up to {dmax} Å");
        for (k, list) in contacts.iter().enumerate() {
            println!("{}", st.atoms[k].label);
            for c in list {
                println!("  {:>6} {:>12}  {:8.4} {:8.4} {:8.4}  {}", st.atoms[c.atom].label, with(c.distance, c.sigma), c.position[0], c.position[1], c.position[2], st.group.ops[c.op].to_xyz());
            }
        }
        if dangl > 0.0 {
            println!("# angles between neighbours within {dangl} Å");
            for a in bonds::angles(&st, esd.as_ref(), &contacts, dangl) {
                let c = &contacts[a.centre];
                println!("  {}-{}-{}  {}", st.atoms[c[a.first].atom].label, st.atoms[a.centre].label, st.atoms[c[a.second].atom].label, with_places(a.angle.0, a.angle.1, 2));
            }
        }
        let occ = if fullprof { bonds::bond_str_occupancies(&st)? } else { bonds::chemical_occupancies(&st)? };
        match bonds::valence_sums(&st, &contacts, &occ, tol) {
            Ok((sums, gii)) => {
                println!("# bond-valence sums (tol {:.0} %)\n  atom  coord  eff.coord  <d>  distortion(1e-4)  valence  sum", tol * 100.0);
                for s in &sums {
                    println!("  {:>6} {:>3} {:6.2} {:>12} {:10.3} {:6.2} {:>12}", st.atoms[s.atom].label, s.coordination, s.effective_coordination, with(s.mean_distance.0, s.mean_distance.1), s.distortion, s.valence, with_places(s.sum.0, s.sum.1, 3));
                }
                println!("  GII {:.2}/100; normalised a {:.2}/100, b {:.2} %, c {:.2}/100", gii.old, gii.a, gii.b, gii.c);
            }
            Err(e) => println!("# no bond-valence sums: {e}"),
        }
    }
    Ok(())
}

/// A .pcr or control file's model, data and parameters, refined `--cycles n` (0 by default) and taken at `--cycle k`
/// of that run (the last by default), with the refinement at that state (its covariance, no further cycle): for
/// `drs state` and `drs structure` (E1, E3, E4).
type StateAt = (drs::powder::Model, drs::data::Pattern, drs::refinable::ParameterSet, drs::lsq::Refinement, usize);

fn state_at_cycle(file: &str, opt: &std::collections::BTreeMap<String, String>) -> Result<StateAt, String> {
    let path = std::path::Path::new(file);
    let (model, data, set) = if file.to_lowercase().ends_with(".pcr") {
        let text: String = std::fs::read(path).map_err(|e| format!("{file}: {e}"))?.iter().map(|&b| b as char).collect();
        let p = drs::pcr::read(&text).map_err(|e| format!("{file}: {e:?}"))?;
        let dir = path.parent().map(|d| d.to_path_buf()).unwrap_or_default();
        let stem = path.file_stem().map(|s| s.to_string_lossy().to_string()).unwrap_or_default();
        let l = drs::powder::Model::from_pcr(&p, &dir, &stem)?;
        let set = drs::refinable::ParameterSet::from_pcr(&p, &l.model)?;
        (l.model, l.data.ok_or(format!("{file}: no data file"))?, set)
    } else {
        let o = open_control(path)?;
        let set = o.control.refinement.build(&o.control.model)?;
        (o.control.model, o.data.ok_or(format!("{file}: no data"))?, set)
    };
    let cycles = opt.get("--cycles").map(|v| v.parse::<usize>().map_err(|_| format!("--cycles {v}"))).transpose()?.unwrap_or(0);
    let mut s = drs::lsq::Session::new(model, data.clone(), set.clone(), drs::lsq::Options { max_cycles: cycles, eps: 0.1 });
    if cycles > 0 {
        s.run(|_| true);
    }
    let t = s.trajectory();
    let k = opt.get("--cycle").map(|v| v.parse::<usize>().map_err(|_| format!("--cycle {v}"))).transpose()?.unwrap_or(t.cycles.len());
    if k > t.cycles.len() {
        return Err(format!("--cycle {k}: the run has {} cycles", t.cycles.len()));
    }
    let m = t.state_at(k);
    let r = drs::lsq::Session::new(m.clone(), data.clone(), set.clone(), drs::lsq::Options { max_cycles: 0, eps: 0.1 }).finish();
    Ok((m, data, set, r, k))
}

fn options_of(args: &[String], flags: &[&str]) -> Result<(String, std::collections::BTreeMap<String, String>), String> {
    let mut file = None;
    let mut opt = std::collections::BTreeMap::new();
    let mut i = 0;
    while i < args.len() {
        if args[i].starts_with("--") {
            if flags.contains(&args[i].as_str()) {
                opt.insert(args[i].clone(), String::new());
                i += 1;
            } else {
                opt.insert(args[i].clone(), args.get(i + 1).ok_or(format!("{} needs a value", args[i]))?.clone());
                i += 2;
            }
        } else {
            file = Some(args[i].clone());
            i += 1;
        }
    }
    Ok((file.ok_or(format!("give a .pcr or control file\n\n{USAGE}"))?, opt))
}

/// `drs state`: a state's derived quantities in one bundle (E3).
fn cmd_state(args: &[String]) -> Result<(), String> {
    use drs::cif_write::with_su;
    if wants_help(args) || args.is_empty() {
        print!("{USAGE}");
        return Ok(());
    }
    let (file, opt) = options_of(args, &[])?;
    let (m, data, set, r, k) = state_at_cycle(&file, &opt)?;
    let step = opt.get("--step").map(|v| num(v, "--step")).transpose()?.unwrap_or(10.0);
    let dmax = opt.get("--dmax").map(|v| num(v, "--dmax")).transpose()?.unwrap_or(3.0);
    let o = drs::state::Options { widths: Some((data.x[0], *data.x.last().unwrap(), step)), dmax, dangl: 0.0, reference_cells: Vec::new() };
    let st = drs::state::summary(&m, &data, Some(&set), Some(&r), &o);
    let s = &st.stats;
    println!("# {file} at cycle {k}: Rp {:.3}  Rwp {:.3}  Rexp {:.3}  chi2 {:.4}; {} reflections", s.rp, s.rwp, s.rexp, s.chi2, st.reflections.len());
    for (p, ph) in st.phases.iter().enumerate() {
        // The cell with its su's (E18a).
        let cs: Vec<String> = ph.cell_su.iter().map(|&(v, e)| with_su(v, e, 6)).collect();
        println!("\nphase {} ({}): a {} b {} c {} alpha {} beta {} gamma {}; V {} A^3", p + 1, ph.name, cs[0], cs[1], cs[2], cs[3], cs[4], cs[5], cs[6]);
        println!("  {:>8} {:>10} {:>8} {:>12} {:>10}", "2theta", "FWHM", "eta", "H_G^2", "H_L");
        for w in &ph.widths {
            println!("  {:8.2} {:10.5} {:8.4} {:12.6} {:10.5}{}", w.two_theta, w.fwhm, w.eta, w.h_g2, w.h_l, if w.h_g2 < 0.0 { "  H_G^2 < 0" } else { "" });
        }
        for (i, cs) in ph.contacts.iter().enumerate() {
            for c in cs {
                // A cation–anion contact's plausible distances (E18b): bond valence at the cation's coordination, ionic radii.
                let plausible = c.plausible.as_ref().map_or(String::new(), |(bv, ionic)| {
                    let b: Vec<String> = bv.iter().map(|(n, d)| format!("bond valence at {n}: {}", d.map_or("none".to_string(), |d| format!("{d:.4}")))).collect();
                    format!("  plausible: {}, ionic radii {ionic:.4}", b.join(", "))
                });
                println!("  {:>6} {:>6} {:>14}  op {} t {:?}{plausible}", ph.labels[i], ph.labels[c.atom], with_su(c.distance, c.sigma, 4), c.op + 1, c.translation);
            }
        }
        if let Some((sums, gii)) = &ph.valence {
            for v in sums {
                println!("  valence sum {:>6} {}  against {} (plausible {:.3}–{:.3}, {})", ph.labels[v.atom], with_su(v.sum.0, v.sum.1, 3), v.expected, v.plausible.0, v.plausible.1, drs::bonds::VALENCE_RANGE_KEY);
            }
            println!("  global instability index {:.4}", gii.a);
        }
        for (n, a, e) in &ph.irreps {
            println!("  irrep {n:8} {}", with_su(*a, *e, 4));
        }
    }
    for q in &st.restraints {
        println!("restraint {:<30} target {:.5} value {:.5}", q.label, q.target, q.value);
    }
    // --at name=value,…: the state with those free parameters at those values (E20a, `lsq::point`).
    if let Some(at) = opt.get("--at") {
        let pairs: Vec<(String, f64)> = at
            .split(',')
            .map(|kv| kv.split_once('=').ok_or(format!("--at {kv}: name=value")).and_then(|(k, v)| Ok((k.trim().to_string(), num(v.trim(), "--at")?))))
            .collect::<Result<_, _>>()?;
        let values: Vec<(&str, f64)> = pairs.iter().map(|(k, v)| (k.as_str(), *v)).collect();
        let pt = drs::lsq::point(&m, &data, &set, &values)?;
        let s = &pt.statistics;
        println!("\nat {at}: Rp {:.3}  Rwp {:.3}  chi2 {:.4} (at cycle {k}: {:.4})", s.rp, s.rwp, s.chi2, st.stats.chi2);
    }
    Ok(())
}

/// `drs structure --cartesian`: a phase's scene for a 3D view (E4).
fn cmd_structure(args: &[String]) -> Result<(), String> {
    if wants_help(args) || args.is_empty() {
        print!("{USAGE}");
        return Ok(());
    }
    let (file, opt) = options_of(args, &["--cartesian"])?;
    let (m, _, _, _, k) = state_at_cycle(&file, &opt)?;
    let p = opt.get("--phase").map(|v| v.parse::<usize>().map_err(|_| format!("--phase {v}"))).transpose()?.unwrap_or(1).max(1) - 1;
    let three = |v: &str| -> Result<[f64; 3], String> {
        let f: Vec<f64> = v.split(',').map(|x| x.trim().parse::<f64>().map_err(|_| format!("--box: {v}"))).collect::<Result<_, _>>()?;
        (f.len() == 3).then(|| [f[0], f[1], f[2]]).ok_or(format!("--box: {v} is not x,y,z"))
    };
    let (lo, hi) = match opt.get("--box") {
        Some(b) => {
            let (a, c) = b.split_once(':').ok_or("--box lo:hi, e.g. 0,0,0:1,1,1")?;
            (three(a)?, three(c)?)
        }
        None => ([0.0; 3], [1.0; 3]),
    };
    let rule = match opt.get("--within") {
        Some(d) => drs::geometry::SceneRule::Within(num(d, "--within")?),
        None => drs::geometry::SceneRule::Ionic(1.2),
    };
    // P19: a phase refined by its coordinates, described by another file's mode phase in the same setting.
    let (m, residual) = match opt.get("--modes") {
        Some(f) => {
            let fp = std::path::Path::new(f);
            let source = if drs::control::is_control_file(fp) {
                open_control(fp)?.control.model
            } else {
                let text: String = std::fs::read(fp).map_err(|e| format!("{f}: {e}"))?.iter().map(|&b| b as char).collect();
                let pcr = drs::pcr::read(&text).map_err(|e| format!("{f}: {e}"))?;
                let dir = std::fs::canonicalize(fp).ok().and_then(|x| x.parent().map(|d| d.to_path_buf())).unwrap_or_else(|| std::path::PathBuf::from("."));
                drs::powder::Model::from_pcr(&pcr, &dir, fp.file_stem().and_then(|s| s.to_str()).unwrap_or(""))?.model
            };
            let metric = match opt.get("--amplimodes") {
                Some(a) => {
                    let text: String = std::fs::read(a).map_err(|e| format!("{a}: {e}"))?.iter().map(|&b| b as char).collect();
                    let k = opt.get("--analysis").map(|v| v.parse::<usize>().map_err(|_| format!("--analysis {v}"))).transpose()?.unwrap_or(1).max(1);
                    drs::modes::read_amplimodes(&text).get(k - 1).and_then(|x| x.reference_cell())
                }
                None => None,
            };
            let (dm, r) = drs::geometry::decompose(&m, p, &source, p, metric, [0.0; 3])?;
            (dm, Some(r))
        }
        None => (m, None),
    };
    let sc = drs::geometry::scene(&m, p, lo, hi, &rule)?;
    if let Some(r) = residual {
        println!("# the coordinates described by {}'s modes: the part they do not describe {r:.2e} Å", opt["--modes"]);
    }
    println!("# {file} at cycle {k}, phase {}: {} images, {} bonds, {} polyhedra; Cartesian axes Å (a along x, b in the xy plane)", p + 1, sc.images.len(), sc.bonds.len(), sc.polyhedra.len());
    for (i, im) in sc.images.iter().enumerate() {
        println!("atom {i:4} {:>6} op {:2} t {:?}  fract {:9.5} {:9.5} {:9.5}  cart {:10.5} {:10.5} {:10.5}", im.label, im.op + 1, im.translation, im.fract[0], im.fract[1], im.fract[2], im.cart[0], im.cart[1], im.cart[2]);
    }
    for (i, j, d) in &sc.bonds {
        println!("bond {i:4} {j:4} {d:.4}");
    }
    for ph in &sc.polyhedra {
        println!("polyhedron {:4} ({}): {} vertices, {} faces, V {:.4} A^3, distortion {:.5}", ph.centre, sc.images[ph.centre].label, ph.vertices.len(), ph.faces.len(), ph.volume, ph.distortion);
    }
    for (a, b) in &sc.edges {
        println!("edge {:10.5} {:10.5} {:10.5}  {:10.5} {:10.5} {:10.5}", a[0], a[1], a[2], b[0], b[1], b[2]);
    }
    for (i, r) in sc.reference_fract.iter().enumerate() {
        println!("reference {i:4} fract {:9.5} {:9.5} {:9.5}  cart {:10.5} {:10.5} {:10.5}", r[0], r[1], r[2], sc.reference[i][0], sc.reference[i][1], sc.reference[i][2]);
    }
    for (name, u) in &sc.displacements {
        for (i, v) in u.iter().enumerate() {
            if v.iter().any(|x| x.abs() > 1e-9) {
                println!("displacement {name:8} {i:4} {:10.5} {:10.5} {:10.5}", v[0], v[1], v[2]);
            }
        }
    }
    Ok(())
}

fn cmd_cif(args: &[String]) -> Result<(), String> {
    if wants_help(args) || args.is_empty() {
        print!("{USAGE}");
        return Ok(());
    }
    let file = &args[0];
    let mut opt = std::collections::BTreeMap::<String, String>::new();
    let mut i = 1;
    while i < args.len() {
        opt.insert(args[i].clone(), args.get(i + 1).ok_or(format!("{} needs a value", args[i]))?.clone());
        i += 2;
    }
    let cycles = opt.get("--cycles").map(|v| num(v, "--cycles")).transpose()?.unwrap_or(0.0) as usize;
    let path = std::path::Path::new(file);
    let (mut model, data, set, title) = if file.to_lowercase().ends_with(".pcr") {
        let text: String = std::fs::read(path).map_err(|e| format!("{file}: {e}"))?.iter().map(|&b| b as char).collect();
        let p = drs::pcr::read(&text).map_err(|e| format!("{file}: {e:?}"))?;
        let dir = path.parent().map(|d| d.to_path_buf()).unwrap_or_default();
        let stem = path.file_stem().map(|s| s.to_string_lossy().to_string()).unwrap_or_default();
        let l = drs::powder::Model::from_pcr(&p, &dir, &stem)?;
        let set = drs::refinable::ParameterSet::from_pcr(&p, &l.model)?;
        (l.model, l.data.ok_or(format!("{file}: no data file"))?, set, p.title.clone())
    } else {
        let o = open_control(path)?;
        let set = o.control.refinement.build(&o.control.model)?;
        (o.control.model, o.data.ok_or(format!("{file}: no data"))?, set, o.control.title)
    };
    let r = drs::lsq::refine(&mut model, &data, &set, drs::lsq::Options { max_cycles: cycles, eps: 0.1 }, |_| {});
    let text = drs::cif_write::deposition(&model, &data, Some(&set), Some(&r), &title)?;
    match opt.get("--out") {
        Some(o) => {
            std::fs::write(o, text).map_err(|e| format!("{o}: {e}"))?;
            println!("{o}: written ({} cycles, Rwp {:.2} %)", r.cycles, r.stats.rwp);
        }
        None => print!("{text}"),
    }
    Ok(())
}

/// `drs results`: the results file of a refinement and its phases' structure CIFs, named and placed by the file
/// dictionary (R28, L101); with --exports the CSV of the reflections and of the profile.
fn cmd_notes(args: &[String]) -> Result<(), String> {
    if wants_help(args) || args.is_empty() {
        print!("{USAGE}");
        return Ok(());
    }
    let path = std::path::Path::new(&args[0]);
    let stem = path.file_stem().map(|s| s.to_string_lossy().to_string()).unwrap_or_default();
    let file = path.parent().unwrap_or(std::path::Path::new(".")).join(drs::formats::files::path("reflection_notes", &[("control", &stem)])?);
    let notes = match std::fs::read_to_string(&file) {
        Ok(t) => drs::notes::read(&t).map_err(|e| format!("{}: {e}", file.display()))?,
        Err(_) => Vec::new(),
    };
    println!("# {} note{} in {}", notes.len(), if notes.len() == 1 { "" } else { "s" }, file.display());
    for n in &notes {
        let place = match (n.phase, n.hkl, n.two_theta) {
            (Some(p), Some([h, k, l]), _) => format!("phase {p} ({h} {k} {l})"),
            (_, _, Some(t)) => format!("2θ {t}°"),
            _ => String::new(),
        };
        println!("{place}\t{}", n.text);
    }
    Ok(())
}

fn cmd_new(args: &[String]) -> Result<(), String> {
    use drs::start::{Instrument, Known};
    if wants_help(args) || args.is_empty() {
        print!("{USAGE}");
        return Ok(());
    }
    let data = std::path::PathBuf::from(&args[0]);
    let (mut known, mut inst, mut out, mut reader) = (Vec::new(), Instrument::default(), None::<std::path::PathBuf>, None::<String>);
    let mut i = 1;
    let next = |i: usize, k: &str| args.get(i + 1).cloned().ok_or(format!("{k} needs a value"));
    while i < args.len() {
        let k = args[i].as_str();
        match k {
            "--cif" => {
                let f = next(i, k)?;
                let text = std::fs::read_to_string(&f).map_err(|e| format!("{f}: {e}"))?;
                let name = std::path::Path::new(&f).file_stem().map_or(f.clone(), |s| s.to_string_lossy().into_owned());
                known.push(Known::Cif { name, text });
                i += 2;
            }
            "--phase" => {
                let (formula, group) = (next(i, k)?, args.get(i + 2).cloned().ok_or("--phase <formula> <space group> <cell values>")?);
                let mut cell = Vec::new();
                let mut j = i + 3;
                while let Some(v) = args.get(j).and_then(|a| a.parse::<f64>().ok()) {
                    cell.push(v);
                    j += 1;
                }
                if cell.is_empty() {
                    return Err("--phase <formula> <space group> <cell values>: the cell's values are missing".into());
                }
                known.push(Known::Cell { formula, group, cell });
                i = j;
            }
            "--wavelength" => {
                let v: Vec<f64> = next(i, k)?.split(',').map(|x| num(x.trim(), k)).collect::<Result<_, _>>()?;
                inst.wavelength = Some(match v.as_slice() {
                    [l] => (*l, None),
                    [l1, l2, r] => (*l1, Some((*l2, *r))),
                    _ => return Err("--wavelength <λ1> or <λ1>,<λ2>,<ratio>".into()),
                });
                i += 2;
            }
            "--anode" => (inst.anode, i) = (Some(next(i, k)?), i + 2),
            "--neutron" => (inst.neutrons, i) = (true, i + 1),
            "--polarisation" => (inst.monochromator_polarisation, i) = (Some(num(&next(i, k)?, k)?), i + 2),
            "--window" => (inst.peak_window, i) = (Some(num(&next(i, k)?, k)?), i + 2),
            "--asymmetry-limit" => (inst.asymmetry_limit, i) = (Some(num(&next(i, k)?, k)?), i + 2),
            "--profile" => {
                let p = next(i, k)?;
                inst.profile = Some(match p.as_str() {
                    "gaussian" => drs::powder::Shape::Gaussian,
                    "lorentzian" => drs::powder::Shape::Lorentzian,
                    "pseudo-voigt" => drs::powder::Shape::PseudoVoigt,
                    "thompson-cox-hastings" => drs::powder::Shape::ThompsonCoxHastings,
                    _ => return Err(format!("--profile {p}: gaussian, lorentzian, pseudo-voigt or thompson-cox-hastings")),
                });
                i += 2;
            }
            "--out" => (out, i) = (Some(next(i, k)?.into()), i + 2),
            "--as" => (reader, i) = (Some(next(i, k)?), i + 2),
            _ => return Err(format!("new: unknown option {k}")),
        }
    }
    let out = out.unwrap_or_else(|| data.with_extension("drc"));
    let s = drs::start::start(&data, reader.as_deref(), &known, &inst, &out)?;
    s.write(&out)?;
    println!("{}: {} phase(s), {} free parameters, data {} ({:.3}–{:.3}°), background {} points in {}", out.display(), s.control.model.phases.len(), s.control.refinement.parameters.len(), data.display(), s.range.0, s.range.1, match &s.control.model.background {
        drs::powder::BackgroundModel::Points(p) => p.len(),
        _ => 0,
    }, s.background_file.0);
    println!("{} peaks found, median FWHM {:.4}°", s.peaks, s.fwhm);
    for d in &s.stated {
        println!("  {} = {}: {}", d.key, d.value, d.source);
    }
    Ok(())
}

/// The value after an option, when given.
fn opt_value(args: &[String], key: &str) -> Option<String> {
    args.iter().position(|a| a == key).and_then(|i| args.get(i + 1)).cloned()
}

/// A run for `drs compare`: its path, its end esd's and correlations (by its free parameters) and its stop (code,
/// text); from a results file read back, or a control file refined here.
struct ComparedRun {
    path: drs::lsq::Trajectory,
    end_esd: Vec<f64>,
    correlations: Option<Vec<Vec<f64>>>,
    stop: Option<(String, String)>,
}

fn compared_run(file: &str) -> Result<ComparedRun, String> {
    let path = std::path::Path::new(file);
    if file.to_lowercase().ends_with(".cif") {
        let text = std::fs::read_to_string(path).map_err(|e| format!("{file}: {e}"))?;
        let r = drs::results::read(&text)?;
        let dir = path.parent().unwrap_or(std::path::Path::new("."));
        let control = drs::control::Control::from_text_in(&r.control_file, dir).map_err(|d| format!("{file}: its control file: {}", d.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("; ")))?;
        let set = control.refinement.build(&control.model)?;
        let names: Vec<String> = set.free.iter().map(|p| p.name.clone()).collect();
        // Each free parameter's esd; a file without them (before core 0.32.0) gives its fields' su's.
        let end_esd = names.iter().map(|n| r.free.iter().chain(&r.refined).find(|f| &f.0 == n).map_or(0.0, |f| f.2)).collect();
        let correlations = r.correlation_matrix(&names);
        let stop = r.stop_code.clone().map(|c| (c, r.stop_text.clone().unwrap_or_default()));
        let path = drs::lsq::Trajectory::from_end(&control.model, &set, r.cycles);
        return Ok(ComparedRun { path, end_esd, correlations, stop });
    }
    let o = open_control(path)?;
    let set = o.control.refinement.build(&o.control.model)?;
    let data = o.data.ok_or(format!("{file}: no data"))?;
    let opts = drs::lsq::Options { max_cycles: o.control.refinement.cycles.unwrap_or(10).max(1), eps: o.control.refinement.eps.unwrap_or(0.1) };
    let mut s = drs::lsq::Session::new(o.control.model, data, set, opts);
    s.run(|_| true);
    let r = s.finish();
    Ok(ComparedRun { path: s.trajectory(), end_esd: r.params.iter().map(|p| p.esd).collect(), correlations: Some(r.correlations), stop: Some((r.reason.code().to_string(), r.stop.clone())) })
}

/// `--prefer a,b,…` among a command's arguments: the order of preference among the readers (request 51).
fn prefer_list(args: &[String]) -> Result<Vec<String>, String> {
    match args.iter().position(|a| a == "--prefer") {
        Some(i) => Ok(args.get(i + 1).ok_or("--prefer needs reader ids, a,b,…")?.split(',').map(|w| w.trim().to_string()).collect()),
        None => Ok(Vec::new()),
    }
}

/// `drs compare-with`: a refinement compared with FullProf or GSAS-II (the app's requests 49, 50).
#[cfg(not(feature = "serde"))]
fn cmd_compare_with(_args: &[String]) -> Result<(), String> {
    Err("compare-with: this drs is built without serde, which the comparison runs' files need".into())
}

/// `drs compare-with`: a refinement compared with FullProf or GSAS-II (the app's requests 49, 50).
#[cfg(feature = "serde")]
fn cmd_compare_with(args: &[String]) -> Result<(), String> {
    use drs::compare;
    if wants_help(args) || args.len() < 4 {
        print!("{USAGE}");
        return Ok(());
    }
    let (program, action, file, dir) = (args[0].as_str(), args[1].as_str(), std::path::Path::new(&args[2]), std::path::Path::new(&args[3]));
    let (mut cycles, mut json) = (None, false);
    let mut i = 4;
    while i < args.len() {
        match args[i].as_str() {
            "--cycles" => {
                cycles = Some(args.get(i + 1).and_then(|v| v.parse::<usize>().ok()).ok_or("--cycles: a positive integer")?);
                i += 2;
            }
            "--json" => (json, i) = (true, i + 1),
            a => return Err(format!("compare-with: {a}?")),
        }
    }
    let opened = open_control(file)?;
    let print_manifest = |m: &[compare::ManifestItem]| {
        for x in m.iter().filter(|x| x.carried != compare::Carried::Carried) {
            println!("  {:?}{}: {}", x.carried, if x.name.is_empty() { String::new() } else { format!(" {}", x.name) }, x.why);
        }
    };
    match (program, action) {
        ("fullprof", "prepare") => {
            let run = compare::fullprof_prepare(&opened, dir, cycles)?;
            println!("{}: {} {} ({} cycles), {}", run.dir.display(), run.program, run.args.join(" "), run.cycles, run.env.iter().map(|(k, v)| format!("{k} = {v}")).collect::<Vec<_>>().join(", "));
            println!("{} parameter(s) carried", run.manifest.iter().filter(|x| x.carried == compare::Carried::Carried).count());
            print_manifest(&run.manifest);
            Ok(())
        }
        ("gsas2", "prepare") => {
            let run = compare::gsas2_prepare(&opened, dir, cycles)?;
            println!("{}: {} {} ({} cycles), {}", run.dir.display(), run.program, run.args.join(" "), run.cycles, run.env.iter().map(|(k, v)| format!("{k} = {v}")).collect::<Vec<_>>().join(", "));
            println!("{} value(s) carried", run.manifest.iter().filter(|x| x.carried == compare::Carried::Carried).count());
            print_manifest(&run.manifest);
            Ok(())
        }
        ("fullprof" | "gsas2", "read") => {
            let run = compare::Run::load(dir)?;
            let c = if program == "fullprof" { compare::fullprof_read(dir, &run, &opened)? } else { compare::gsas2_read(dir, &run, &opened)? };
            let path = compare::write(dir, &c)?;
            if json {
                println!("{}", serde_json::to_string(&c).map_err(|e| e.to_string())?);
                return Ok(());
            }
            print_comparison(&c);
            println!("written: {}", path.display());
            Ok(())
        }
        _ => Err(format!("compare-with {program} {action}: fullprof|gsas2 prepare|read")),
    }
}

/// A comparison as `drs compare-with … read` prints it.
#[cfg(feature = "serde")]
fn print_comparison(c: &drs::compare::Comparison) {
    let conv = |p: &drs::compare::Program| format!("{} cycles, {}", p.cycles, if p.converged { "converged" } else { "not converged" });
    println!("{} ({}; {}) against {} ({}; {})", c.core.name, c.core.version, conv(&c.core), c.other.name, c.other.version, conv(&c.other));
    println!("\n{:28} {:>12} {:>12} {:>12}", "agreement", "this app", c.other.name, "difference");
    for m in &c.agreement {
        let o = m.other.map_or("—".to_string(), |v| format!("{v:.4}"));
        let d = m.difference.map_or("—".to_string(), |v| format!("{v:+.4}"));
        println!("{:28} {:>12.4} {:>12} {:>12}", m.name, m.core, o, d);
    }
    println!("\n{:36} {:>24} {:>24} {:>8} {:>7}  verdict ({})", "parameter", "this app", c.other.name, "Δ/σ", "σ ratio", c.rule);
    for p in &c.parameters {
        println!(
            "{:36} {:>24} {:>24} {:>8.2} {:>7.3}  {}{}",
            p.address,
            drs::cif_write::with_su(p.core_value, p.core_esd, 1),
            drs::cif_write::with_su(p.other_value, p.other_esd, 1),
            p.delta_sigma,
            p.esd_ratio,
            if p.agrees { "agrees" } else { "differs" },
            p.note.as_ref().map_or(String::new(), |n| format!(" — {n}"))
        );
    }
    let agree = c.parameters.iter().filter(|p| p.agrees).count();
    println!("{agree} of {} parameters agree", c.parameters.len());
    if !c.not_compared.is_empty() {
        println!("\nnot compared:");
        for n in &c.not_compared {
            println!("  {} (refined by {}): {}", n.name, n.refined_by, n.why);
        }
    }
    for n in &c.notes {
        println!("note: {n}");
    }
}

/// `drs compare`: two refinements' paths compared (E11).
fn cmd_compare(args: &[String]) -> Result<(), String> {
    use drs::lsq::{Run, Side};
    if wants_help(args) || args.len() < 2 {
        print!("{USAGE}");
        return Ok(());
    }
    let (mut files, mut reference, mut pair, mut json) = (Vec::new(), Side::A, None::<(String, String)>, false);
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--json" => (json, i) = (true, i + 1),
            "--reference" => {
                reference = match args.get(i + 1).map(String::as_str) {
                    Some("a") | Some("A") => Side::A,
                    Some("b") | Some("B") => Side::B,
                    o => return Err(format!("--reference {}: a or b", o.unwrap_or(""))),
                };
                i += 2;
            }
            "--pair" => {
                let (a, b) = (args.get(i + 1).cloned(), args.get(i + 2).cloned());
                pair = Some((a.ok_or("--pair needs two addresses")?, b.ok_or("--pair needs two addresses")?));
                i += 3;
            }
            f => (files, i) = ({ files.push(f.to_string()); files }, i + 1),
        }
    }
    let [fa, fb] = &files[..] else { return Err("compare: give two files, a and b".into()) };
    let (ra, rb) = (compared_run(fa)?, compared_run(fb)?);
    fn run(r: &ComparedRun) -> Run<'_> {
        Run { path: &r.path, end_esd: Some(&r.end_esd), end_correlations: r.correlations.as_deref(), stop: r.stop.as_ref().map(|(c, t)| (c.as_str(), t.as_str())) }
    }
    let c = drs::lsq::compare(&run(&ra), &run(&rb), reference);
    let p = pair.as_ref().map(|(a, b)| c.pair(a, b)).transpose()?;
    if json {
        #[cfg(feature = "serde")]
        {
            let text = match &p {
                Some(p) => format!("{{\"paths\": {}, \"pair\": {}}}", serde_json::to_string(&c).map_err(|e| e.to_string())?, serde_json::to_string(p).map_err(|e| e.to_string())?),
                None => serde_json::to_string(&c).map_err(|e| e.to_string())?,
            };
            println!("{text}");
            return Ok(());
        }
        #[cfg(not(feature = "serde"))]
        return Err("--json needs the core built with serde".into());
    }
    let s = &c.summary;
    println!("# a: {fa}\n# b: {fb}\n# reference {:?}", c.reference);
    println!("{} in common, {} only in a, {} only in b, {} fixed values changed; ends: {} agree, {} disagree (within {} combined σ)", s.common, s.only_a, s.only_b, s.fixed_differences, s.agree, s.disagree, drs::lsq::ENDS_AGREE);
    let stop = |x: &Option<(String, String)>| x.as_ref().map_or("not known".to_string(), |(c, t)| format!("{c} ({t})"));
    println!("a: {} states, last χ² {}, stop {}", c.states_a.len(), s.last_chi2_a.map_or("-".into(), |v| v.to_string()), stop(&s.stop_a));
    println!("b: {} states, last χ² {}, stop {}", c.states_b.len(), s.last_chi2_b.map_or("-".into(), |v| v.to_string()), stop(&s.stop_b));
    for (side, list) in [("a", &c.only_a), ("b", &c.only_b)] {
        for (n, v0, v1) in list.iter() {
            println!("only in {side}: {n} from {v0} to {v1}");
        }
    }
    for (n, va, vb) in &c.fixed_differences {
        println!("fixed, changed: {n} a {va} b {vb}");
    }
    let mut order: Vec<usize> = (0..c.end.len()).collect();
    order.sort_by(|&x, &y| c.end[y].combined.map_or(0.0, f64::abs).total_cmp(&c.end[x].combined.map_or(0.0, f64::abs)));
    println!("\n{:<45} {:>14} {:>12} {:>10}  unit", "parameter", "b − a", "/σ ref", "/σ comb");
    for k in order {
        let e = &c.end[k];
        let f = |v: Option<f64>| v.map_or("-".to_string(), |v| format!("{v:.3}"));
        println!("{:<45} {:>14.6e} {:>12} {:>10}{} {}", e.address, e.delta, f(e.in_scale), f(e.combined), if e.agree { " " } else { "*" }, c.units[k]);
    }
    if let Some(p) = p {
        println!("\npair {} – {}", p.a, p.b);
        for (side, centre, sig, r, el) in [("a", p.centre_a, p.sigmas_a, p.correlation_a, p.ellipse_a), ("b", p.centre_b, p.sigmas_b, p.correlation_b, p.ellipse_b)] {
            println!("{side}: end ({}, {}), σ ({}, {}), r {}, ellipse {:?}", centre.0, centre.1, sig.0, sig.1, r.map_or("-".into(), |r| format!("{r:.4}")), el);
        }
    }
    Ok(())
}

fn crate_range(p: &drs::data::Pattern) -> String {
    drs::data::extent(&p.x).map_or("-".into(), |(a, b)| format!("{a:.5} to {b:.5}"))
}

/// `drs plots seq`: FullProf's `.seq` as a table (PL11).
fn plots_seq(path: &str, parameter: Option<String>, json: bool) -> Result<(), String> {
    let text = std::fs::read(path).map_err(|e| format!("{path}: {e}"))?;
    let seq = drs::seq::read(&String::from_utf8_lossy(&text))?;
    if json {
        return print_json(&seq);
    }
    let names = seq.names();
    println!("# {path}: {} blocks, {} parameters refined in some block{}", seq.blocks.len(), names.len(), if seq.unread.is_empty() { String::new() } else { format!(", {} lines not read", seq.unread.len()) });
    for b in &seq.blocks {
        let phases: Vec<String> = b
            .phases
            .iter()
            .map(|p| format!("ph{} R_Bragg {} V {} wf {}", p.phase, p.r_bragg.map_or("-".into(), |v| format!("{v}")), p.volume.map_or("-".into(), |v| format!("{}({})", v.0, v.1)), p.weight_fraction.map_or("-".into(), |v| format!("{}({})", v.0, v.1))))
            .collect();
        let value = parameter.as_ref().map(|n| b.parameters.iter().find(|p| &p.name == n).map_or(format!("  {n} -"), |p| format!("  {n} {} ± {}", p.value, p.sigma)));
        println!(
            "{:10} {:4} {:10} {:7} chi2 {:9} {}  {} parameters{}",
            b.name,
            b.index.map_or("-".into(), |v| v.to_string()),
            b.parameter.map_or("-".into(), |v| v.to_string()),
            b.numor.map_or("-".into(), |v| v.to_string()),
            b.chi2.map_or("-".into(), |v| v.to_string()),
            phases.join("; "),
            b.parameters.len(),
            value.unwrap_or_default()
        );
    }
    for (n, l) in &seq.unread {
        println!("not read, line {n}: {l}");
    }
    Ok(())
}

#[cfg(feature = "serde")]
fn print_json<T: serde::Serialize>(v: &T) -> Result<(), String> {
    println!("{}", serde_json::to_string(v).map_err(|e| e.to_string())?);
    Ok(())
}

#[cfg(not(feature = "serde"))]
fn print_json<T>(_: &T) -> Result<(), String> {
    Err("--json needs the core built with serde".into())
}

/// `drs plots`: several files as one set (R32; PL1–PL5, PL9, PL11, PL17).
/// The series a run's files offer (W2): each file read alone, with its modification time as a browser's host gives it,
/// a buffer's files with its comments; the proposals' titles, the classes and the files left out, or the whole as JSON.
fn plots_propose(files: &[String], json: bool, x_interval: Option<[f64; 2]>) -> Result<(), String> {
    use drs::plots::{self, OpenedFile};
    let opened = |path: &std::path::Path, name: String, comment: Option<String>| -> Result<OpenedFile, String> {
        let bytes = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
        let modified = std::fs::metadata(path).and_then(|m| m.modified()).ok().and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok()).map(|d| d.as_secs_f64());
        Ok(OpenedFile { name, text: String::from_utf8_lossy(&bytes).into_owned(), comment, modified, ..Default::default() })
    };
    let given: Vec<OpenedFile> = if files.len() == 1 && files[0].to_lowercase().ends_with(".buf") {
        let path = std::path::Path::new(&files[0]);
        let base = path.parent().unwrap_or(std::path::Path::new(""));
        let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", files[0]))?;
        plots::read_buf(&text, &|n| base.join(n).is_file()).entries.into_iter().map(|(n, c)| opened(&base.join(&n), n, c)).collect::<Result<_, _>>()?
    } else {
        files.iter().map(|f| opened(std::path::Path::new(f), f.clone(), None)).collect::<Result<_, _>>()?
    };
    let p = plots::propose(&given, x_interval);
    if json {
        return print_json(&p);
    }
    let classes: Vec<String> = p.classes.iter().map(|c| format!("{} scans, σ²/y {:.3} to {:.3}", c.scans, c.median_scale.0, c.median_scale.1)).collect();
    println!("# {} set(s) offered; exposure classes: {}", p.proposals.len(), classes.join("; "));
    for (k, x) in p.proposals.iter().enumerate() {
        println!("{}. {}: {}", k + 1, x.kind, x.title);
        // the first few scans it leaves out, with why; --json names every one
        for (f, why) in x.left_out.iter().take(3) {
            println!("   left out {f}: {why}");
        }
        if x.left_out.len() > 3 {
            println!("   … and {} more left out (--json names each with why)", x.left_out.len() - 3);
        }
    }
    for (f, why) in &p.left_out {
        println!("in no set {f}: {why}");
    }
    Ok(())
}

fn cmd_plots(args: &[String]) -> Result<(), String> {
    use drs::plots::{self, Displacement, YScale};
    if wants_help(args) || args.len() < 2 {
        print!("{USAGE}");
        return Ok(());
    }
    let sub = args[0].as_str();
    let (mut files, mut reader, mut json, mut sets) = (Vec::<String>::new(), None::<String>, false, Vec::<(usize, f64, String)>::new());
    let (mut rule, mut scale) = (Displacement::None, YScale::Linear);
    let (mut dx, mut step, mut per_unit) = (0.0, None::<f64>, None::<f64>);
    let mut i = 1;
    while i < args.len() {
        let a = args[i].as_str();
        let next = || args.get(i + 1).cloned().ok_or(format!("{a} needs a value"));
        match a {
            "--as" => (reader, i) = (Some(next()?), i + 2),
            "--json" => (json, i) = (true, i + 1),
            "--dx" => (dx, i) = (num(&next()?, a)?, i + 2),
            "--step" => (step, i) = (num(&next()?, a).ok(), i + 2),
            "--per-unit" => (per_unit, i) = (Some(num(&next()?, a)?), i + 2),
            "--scale" => {
                scale = match next()?.as_str() {
                    "linear" => YScale::Linear,
                    "sqrt" => YScale::Sqrt,
                    "log" => YScale::Log,
                    o => return Err(format!("--scale {o}: linear, sqrt or log")),
                };
                i += 2;
            }
            "--set" => {
                let v = next()?;
                let (k, rest) = v.split_once('=').ok_or("--set <k>=<value>:<note>")?;
                let (val, note) = rest.split_once(':').unwrap_or((rest, "typed on the command line"));
                sets.push((k.parse::<usize>().map_err(|_| format!("--set {v}"))?, num(val, "--set")?, note.to_string()));
                i += 2;
            }
            "--range" | "--reference" | "--members" | "--out" | "--name" | "--window" | "--background" | "--log-lambda" | "--at" | "--order" | "--floor" | "--parameter" | "--x-interval" => i += 2,
            "--csv" | "--long" | "--browser" | "--no-open" | "--no-patterns" => i += 1,
            "--x" | "--wavelength" | "--temperatures" | "--by" | "--to" | "--within" | "--member" | "--lambda-from" | "--lambda" | "--view" => i += 2,
            "--sum" | "--average" | "--weighted-average" | "--difference" => i += 1,
            _ => (files, i) = ({ files.push(a.to_string()); files }, i + 1),
        }
    }
    if files.is_empty() {
        return Err("plots: give a .buf or data files".into());
    }
    if sub == "seq" {
        return plots_seq(&files[0], opt_value(args, "--parameter"), json);
    }
    // The x interval of the analysis (PL42): every member cut to it on opening, as a set file's `x_interval`.
    let x_interval = match opt_value(args, "--x-interval") {
        Some(v) => {
            let (a, b) = v.split_once(':').ok_or("--x-interval a:b")?;
            Some([num(a, "--x-interval")?, num(b, "--x-interval")?])
        }
        None => None,
    };
    if sub == "propose" {
        return plots_propose(&files, json, x_interval);
    }
    let mut set_file = None;
    let mut set = if files.len() == 1 && files[0].to_lowercase().ends_with(".buf") {
        plots::open_buf(std::path::Path::new(&files[0]), reader.as_deref())?
    } else if files.len() == 1 && files[0].to_lowercase().ends_with(".drv") {
        let path = std::path::Path::new(&files[0]);
        let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", files[0]))?;
        let f = plots::SetFile::read(&text).map_err(|d| format!("{}: {}", files[0], d.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("; ")))?;
        let base = BASE.get().cloned().unwrap_or_else(|| plots::set_base(path));
        let s = f.open(&base, reader.as_deref(), &|_| Vec::new())?;
        if !args.iter().any(|a| a == "--scale") {
            scale = f.y_scale;
        }
        set_file = Some(f);
        s
    } else {
        plots::open(std::path::Path::new(""), &files.iter().map(|f| (f.clone(), None)).collect::<Vec<_>>(), reader.as_deref(), None)?
    };
    if let Some(iv) = x_interval {
        // in the set file's x quantity, else in the members' own
        set.cut(iv, set_file.as_ref().and_then(|f: &plots::SetFile| f.x_quantity()))?;
        if let Some(f) = set_file.as_mut() {
            f.x_interval = set.x_interval;
        }
    }
    // FullProf's temperatures file as the variable's source (PL21): each member its file names takes its value.
    if let Some(f) = opt_value(args, "--temperatures") {
        let t = plots::read_temperatures(&std::fs::read_to_string(&f).map_err(|e| format!("{f}: {e}"))?);
        for (name, values) in &t.conflicts {
            eprintln!("{f}: {name} is given {} values ({}): it takes none", values.len(), values.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(", "));
        }
        for (n, l) in &t.unread {
            eprintln!("{f}, line {n}: not \"name TSAMP value\": {l}");
        }
        let name = std::path::Path::new(&f).file_name().map_or(f.clone(), |x| x.to_string_lossy().into_owned());
        let n = set.variable_from_file(&name, &t, true);
        eprintln!("{f}: {n} of {} members took their temperature from it", set.members.len());
    }
    for (k, v, note) in &sets {
        set.set_variable(k.checked_sub(1).ok_or("--set: members count from 1")?, *v, note)?;
    }
    match sub {
        "open" => {
            if json {
                #[cfg(feature = "serde")]
                {
                    println!("{}", serde_json::to_string(&set).map_err(|e| e.to_string())?);
                    return Ok(());
                }
                #[cfg(not(feature = "serde"))]
                return Err("--json needs the core built with serde".into());
            }
            println!("# {} members{}; variable {} ({})", set.members.len(), set.title.as_ref().map_or(String::new(), |t| format!(", \"{t}\"")), set.variable_name, set.variable_unit);
            for (k, m) in set.members.iter().enumerate() {
                let f = &m.facts;
                println!(
                    "{:3} {:24} {:14} {:5} points {:9.4}–{:9.4} step {} λ {}  {}",
                    k + 1,
                    m.label,
                    f.format,
                    f.points,
                    f.range.0,
                    f.range.1,
                    f.step.map_or("-".into(), |s| format!("{s:.5}")),
                    f.wavelength.map_or("-".into(), |l| format!("{l}")),
                    match &m.variable {
                        Some(v) => format!("{} = {} {} ({})", set.variable_name, v.value, set.variable_unit, v.found.words()),
                        None => format!("{}: nothing says it", set.variable_name),
                    }
                );
                for c in m.candidates.iter().skip(1) {
                    if Some(c.value) != m.variable.as_ref().map(|v| v.value) {
                        println!("    also {} from {}", c.value, c.found.words());
                    }
                }
            }
            Ok(())
        }
        "displace" | "project" => {
            rule = match (step, per_unit) {
                (Some(dy), None) => Displacement::Step { dx, dy },
                (None, Some(dy)) => Displacement::Variable { dx, dy },
                (None, None) if dx != 0.0 => Displacement::Step { dx, dy: 0.0 },
                (None, None) => set_file.as_ref().map_or(rule, |f| f.displacement),
                _ => return Err("give --step or --per-unit, not both".into()),
            };
            let d = plots::displace(&set, rule, scale, &[], None);
            if sub == "project" {
                let order: Vec<usize> = match opt_value(args, "--order") {
                    Some(o) => o.split(',').map(|v| v.trim().parse::<usize>().ok().and_then(|k| k.checked_sub(1)).ok_or(format!("--order: {v}"))).collect::<Result<_, _>>()?,
                    None => (0..d.len()).filter(|&k| d[k].is_some()).collect(),
                };
                let floor = opt_value(args, "--floor").map(|v| num(&v, "--floor")).transpose()?;
                let pr = plots::project(&d, &order, floor)?;
                if json {
                    return print_json(&pr);
                }
                println!("# the projection, front to back, {:?} on the {:?} scale{}", rule, scale, floor.map_or(String::new(), |f| format!(", each curve hiding down to {f} plus its shift")));
                for c in &pr.curves {
                    println!("{:24} {:5.1} % visible in {} piece{}", set.members[c.member].label, 100.0 * c.visible, c.pieces.len(), if c.pieces.len() == 1 { "" } else { "s" });
                }
                return Ok(());
            }
            if json {
                #[cfg(feature = "serde")]
                {
                    println!("{}", serde_json::to_string(&d).map_err(|e| e.to_string())?);
                    return Ok(());
                }
                #[cfg(not(feature = "serde"))]
                return Err("--json needs the core built with serde".into());
            }
            println!("# {:?} on the {:?} scale", rule, scale);
            for (m, o) in set.members.iter().zip(&d) {
                match o {
                    Some(o) => println!("{:24} shift x {:+.4} y {:+.4}  baseline {:.4}", m.label, o.shift.0, o.shift.1, o.baseline),
                    None => println!("{:24} not placed: no value of {}", m.label, set.variable_name),
                }
            }
            Ok(())
        }
        "save" => {
            let mut f = set_file.clone().unwrap_or_else(|| plots::SetFile::of(&set));
            // The view's normalisation (PL27), with the options of drs plots normalise.
            if let Some(rule) = norm_rule(args)? {
                f.normalisation = rule;
            }
            // The values chosen on this command (--temperatures, --set) kept with their sources (PL23).
            let now = plots::SetFile::of(&set);
            for m in f.members.iter_mut() {
                let mine: Vec<&plots::SetMember> = now.members.iter().filter(|x| x.file == m.file && (m.scan.is_none() || x.scan == m.scan)).collect();
                if let [x] = mine[..] {
                    m.value = x.value.clone();
                }
            }
            let name = opt_value(args, "--name").or_else(|| set.title.clone()).unwrap_or_else(|| "plot set".into());
            let out = match opt_value(args, "--out") {
                Some(o) => o,
                None => drs::formats::files::path("plot_set", &[("set", &name)])?,
            };
            // The files named from the new file's project folder (the folder above Plots/, or --base).
            let from = if files.len() == 1 && (files[0].to_lowercase().ends_with(".buf") || files[0].to_lowercase().ends_with(".drv")) {
                let p = std::path::Path::new(&files[0]);
                if files[0].to_lowercase().ends_with(".drv") { BASE.get().cloned().unwrap_or_else(|| plots::set_base(p)) } else { p.parent().map(|d| d.to_path_buf()).unwrap_or_default() }
            } else {
                std::path::PathBuf::new()
            };
            let to = BASE.get().cloned().unwrap_or_else(|| plots::set_base(std::path::Path::new(&out)));
            for m in f.members.iter_mut() {
                m.file = plots::relative_to(&from.join(&m.file), &to);
            }
            if let Some(r) = f.reference.as_mut() {
                *r = plots::relative_to(&from.join(&*r), &to);
            }
            if !out.to_lowercase().ends_with(".drv") {
                return Err(format!("{out}: a plot set is a .drv file"));
            }
            let path = std::path::Path::new(&out);
            if path.exists() {
                return Err(format!("{out} exists: nothing written"));
            }
            if let Some(dir) = path.parent().filter(|d| !d.as_os_str().is_empty()) {
                std::fs::create_dir_all(dir).map_err(|e| format!("{}: {e}", dir.display()))?;
            }
            std::fs::write(path, f.to_text()).map_err(|e| format!("{out}: {e}"))?;
            println!("{out}: {} members (files relative to the project's folder, the folder above Plots/)", f.members.len());
            Ok(())
        }
        "step" => {
            let p = plots::proposed_step(&set, scale).ok_or("no member has two values on this scale")?;
            if json {
                return print_json(&p);
            }
            println!("dy {}  ({})", p.dy, p.rule);
            Ok(())
        }
        "export" => {
            if !args.iter().any(|a| a == "--csv") {
                return Err("export: give --csv".into());
            }
            let to = match opt_value(args, "--x").as_deref() {
                None | Some("2theta") => drs::formats::Axis::TwoTheta,
                Some("d") => drs::formats::Axis::D,
                Some("q") => drs::formats::Axis::Q,
                Some(o) => return Err(format!("--x {o}: 2theta, d or q")),
            };
            let lambda = opt_value(args, "--wavelength").map(|v| num(&v, "--wavelength")).transpose()?;
            // As the view shows it (PL25): normalised in the files' x, then converted; the file says how.
            let (set, notes) = match norm_rule(args)? {
                Some(rule) => set.as_shown(&rule)?,
                None => (set, Vec::new()),
            };
            let set = plots::converted(&set, to, lambda)?;
            let range = match opt_value(args, "--range").as_deref() {
                None | Some("overlap") => plots::GridRange::Overlap,
                Some("union") => plots::GridRange::Union,
                Some(o) => return Err(format!("--range {o}: overlap or union")),
            };
            let gstep = match opt_value(args, "--step").as_deref() {
                None | Some("finest") => plots::GridStep::Finest,
                Some("coarsest") => plots::GridStep::Coarsest,
                Some(v) => plots::GridStep::Given(num(v, "--step")?),
            };
            let text = if args.iter().any(|a| a == "--long") { plots::csv_long_with(&set, &notes) } else { plots::csv_with(&set, &plots::common_grid(&set, range, gstep)?, &notes) };
            match opt_value(args, "--out") {
                Some(o) => {
                    if std::path::Path::new(&o).exists() {
                        return Err(format!("{o} exists: nothing written"));
                    }
                    std::fs::write(&o, text).map_err(|e| format!("{o}: {e}"))?;
                    println!("{o}: {} members", set.members.len());
                }
                None => print!("{text}"),
            }
            Ok(())
        }
        "view" => {
            let q = opt_value(args, "--x").ok_or("view: --x 2theta|d|q|stl|s")?;
            let quantity = plots::XQuantity::from_word(&q).ok_or(format!("--x {q}: 2theta, d, q, stl or s"))?;
            let given: Vec<(usize, f64, plots::LambdaFrom)> = args
                .windows(2)
                .filter(|w| w[0] == "--lambda")
                .map(|w| {
                    let (k, v) = w[1].split_once('=').ok_or(format!("--lambda {}: k=λ", w[1]))?;
                    let k = k.parse::<usize>().ok().and_then(|k| k.checked_sub(1)).ok_or(format!("--lambda {}: members count from 1", w[1]))?;
                    Ok((k, num(v, "--lambda")?, plots::LambdaFrom::Typed("typed on the command line".into())))
                })
                .collect::<Result<_, String>>()?;
            let fallback = match (opt_value(args, "--lambda-from"), opt_value(args, "--wavelength")) {
                (Some(f), _) => {
                    let text = std::fs::read_to_string(&f).map_err(|e| format!("{f}: {e}"))?;
                    let c = drs::control::Control::from_text(&text).map_err(|d| format!("{f}: {}", d.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("; ")))?;
                    Some((c.model.lambda[0], plots::LambdaFrom::ControlFile(std::path::Path::new(&f).file_name().map_or(f.clone(), |x| x.to_string_lossy().into_owned()))))
                }
                (None, Some(v)) => Some((num(&v, "--wavelength")?, plots::LambdaFrom::Typed("typed on the command line".into()))),
                (None, None) => None,
            };
            let v = plots::x_view(&set, quantity, &given, fallback);
            if json {
                return print_json(&v);
            }
            println!("# in {} (a view: no file and no operation changed)", quantity.label());
            for x in &v.members {
                let m = &set.members[x.member];
                match &x.pattern {
                    Some(p) => println!("{:24} {:5} points {}{}", m.label, p.len(), crate_range(p), x.lambda.map_or(String::new(), |l| format!("  λ {l} Å from {:?}", x.from.clone().unwrap()))),
                    None => println!("{:24} not shown: {}", m.label, x.why_not.clone().unwrap_or_default()),
                }
            }
            if !v.without.is_empty() {
                println!("without a wavelength: {}", v.without.iter().map(|&k| set.members[k].label.clone()).collect::<Vec<_>>().join(", "));
            }
            Ok(())
        }
        "normalise" => {
            // Without --by, the set file's own normalisation (PL27).
            let rule = match norm_rule(args)? {
                Some(r) => r,
                None => set_file.as_ref().map(|f| f.normalisation.clone()).filter(|n| *n != plots::Normalisation::None).ok_or("normalise: --by monitor|time|maximum|area|peak|reference")?,
            };
            let n = plots::normalise(&set, &rule)?;
            if json {
                return print_json(&n);
            }
            println!("# {rule:?} (a view: no file and no operation changed)");
            for x in &n {
                let m = &set.members[x.member];
                match x.factor {
                    Some(f) => println!("{:24} × {f}", m.label),
                    None => println!("{:24} not normalised: {}", m.label, x.why_not.clone().unwrap_or_default()),
                }
            }
            Ok(())
        }
        "windows" => {
            let mut list: Vec<plots::Window> = args
                .windows(2)
                .filter(|w| w[0] == "--window")
                .map(|w| {
                    let p: Vec<&str> = w[1].split(':').collect();
                    match p[..] {
                        [name, a, b] => Ok(plots::Window { name: name.to_string(), from: num(a, "--window")?, to: num(b, "--window")? }),
                        _ => Err(format!("--window {}: name:from:to", w[1])),
                    }
                })
                .collect::<Result<_, String>>()?;
            if list.is_empty() {
                list = set_file.as_ref().map(|f| f.windows.clone()).unwrap_or_default();
            }
            if list.is_empty() {
                return Err("give one --window name:from:to or more".into());
            }
            let log_lambda = opt_value(args, "--log-lambda").map(|v| num(&v, "--log-lambda")).transpose()?;
            let background = match opt_value(args, "--background").as_deref() {
                None | Some("auto") => plots::WindowBackground::Auto { log_lambda },
                Some("line") => plots::WindowBackground::Line,
                Some(f) => {
                    let pf = drs::background::read(&std::fs::read_to_string(f).map_err(|e| format!("{f}: {e}"))?)?;
                    plots::WindowBackground::Points { points: pf.points, interpolation: pf.interpolation.unwrap_or(drs::background::Interpolation::Linear) }
                }
            };
            let at: Vec<f64> = match opt_value(args, "--at") {
                Some(v) => v.split(',').map(|x| num(x.trim(), "--at")).collect::<Result<_, _>>()?,
                None => Vec::new(),
            };
            let w = plots::windows(&set, &list, &background, &at)?;
            if json {
                return print_json(&w);
            }
            println!("# windows above {}", match &background {
                plots::WindowBackground::Auto { .. } => "each member's automatic background".to_string(),
                plots::WindowBackground::Line => "the line between each window's ends".to_string(),
                plots::WindowBackground::Points { points, interpolation } => format!("{} given points ({})", points.len(), interpolation.name()),
            });
            for mw in &w.members {
                let m = &set.members[mw.member];
                let v = m.variable.as_ref().map_or("-".to_string(), |v| format!("{} = {}", set.variable_name, v.value));
                let cells: Vec<String> = mw
                    .windows
                    .iter()
                    .map(|x| format!("{} {:.6} ± {:.6}{}{}", x.name, x.integral, x.sigma, x.centroid.map_or(String::new(), |c| format!(" at {c:.4}")), x.peak.as_ref().map_or(String::new(), |q| format!(" peak {:.3}", q.x))))
                    .collect();
                let bg: Vec<String> = w.at.iter().zip(&mw.background_at).map(|(x, b)| format!("b({x}) {}", b.map_or("-".into(), |b| format!("{b:.3}")))).collect();
                println!("{:16} {:10} {}  {}", m.label, v, cells.join("; "), bg.join(" "));
            }
            Ok(())
        }
        "describe" => {
            // What a set file says (PL12, PL27, PL37): its views, the map's clip and its regions, each checked.
            let f = set_file.as_ref().ok_or("describe: give a plot set file (.drv)")?;
            if json {
                return print_json(f);
            }
            for l in f.words() {
                println!("{l}");
            }
            if let Err(e) = f.check_regions(&set) {
                println!("refused: {e}");
            }
            Ok(())
        }
        "show" => {
            // The set's map in the app when it is installed, else in a page of its own (PL38, I08 step 5).
            let f = set_file.as_ref().ok_or("show: give a plot set file (.drv)")?;
            let path = std::path::Path::new(&files[0]);
            let view = opt_value(args, "--view").unwrap_or_else(|| "map".into());
            if view != "map" && view != "difference" {
                return Err(format!("--view {view}: map or difference"));
            }
            let reference = opt_value(args, "--reference").map(|v| v.parse::<usize>().map_err(|_| format!("--reference {v}")).and_then(|k| k.checked_sub(1).ok_or("--reference counts from 1".to_string()))).transpose()?;
            let reference = reference.or_else(|| f.reference_in(&set));
            if view == "difference" && reference.is_none() {
                return Err("--view difference: the set file names no reference member (reference = …), so there is no difference map; give --reference <k>".into());
            }
            if !args.iter().any(|a| a == "--browser") && app_installed() {
                let o = std::process::Command::new("open").args(["-a", "diff-ref-struct"]).arg(path).args(["--args", "--view", &view]).status().map_err(|e| format!("open -a diff-ref-struct: {e}"))?;
                if !o.success() {
                    return Err(format!("open -a diff-ref-struct {}: {o}", files[0]));
                }
                println!("{}: opened in diff-ref-struct on the {view} view (the app takes --view from PL39; an app before it opens Plots on the view the file says)", files[0]);
                return Ok(());
            }
            // The set as the view shows it: normalised as the file says, in the file's x.
            let set = if f.normalisation == plots::Normalisation::None { set } else { set.as_shown(&f.normalisation)?.0 };
            let as_read = plots::MapX::of(&set);
            let (set, mx) = match f.x_quantity() {
                Some(q) if as_read.quantity != Some(q) => {
                    let fallback = opt_value(args, "--wavelength").map(|v| num(&v, "--wavelength")).transpose()?.map(|l| (l, plots::LambdaFrom::Typed("--wavelength".into())));
                    let v = plots::x_view(&set, q, &[], fallback);
                    if !v.without.is_empty() {
                        return Err(format!("the set file shows the set in {}, and these members state no wavelength to convert them: {}; give --wavelength", q.word(), v.without.iter().map(|&k| set.members[k].label.as_str()).collect::<Vec<_>>().join(", ")));
                    }
                    let mx = plots::MapX::of_view(&set, &v);
                    let mut shown = set.clone();
                    for x in &v.members {
                        match &x.pattern {
                            Some(p) => shown.members[x.member].pattern = p.clone(),
                            None => return Err(format!("{}: {}", set.members[x.member].label, x.why_not.clone().unwrap_or_else(|| "not shown in this quantity".into()))),
                        }
                    }
                    (shown, mx)
                }
                _ => (set, as_read),
            };
            let mut data = plots::page_data(&set, f, &mx, reference, &files[0], !args.iter().any(|a| a == "--no-patterns"))?;
            data.view = view.clone();
            let html = plots::page_html(&data)?;
            let out = match opt_value(args, "--out") {
                Some(o) => std::path::PathBuf::from(o),
                None => path.with_extension("map.html"),
            };
            // Over a page this command wrote, never over another file.
            if out.exists() && !std::fs::read_to_string(&out).is_ok_and(|t| t.contains("id=\"drs-page-data\"")) {
                return Err(format!("{} exists and is not a page of drs plots show: nothing written", out.display()));
            }
            std::fs::write(&out, &html).map_err(|e| format!("{}: {e}", out.display()))?;
            println!("{}: the page of {} ({} members, {} kB, core {}, {}), nothing fetched", out.display(), files[0], set.members.len(), html.len().div_ceil(1000), data.core, data.written);
            if !args.iter().any(|a| a == "--no-open") {
                let opener = if cfg!(target_os = "macos") { "open" } else { "xdg-open" };
                match std::process::Command::new(opener).arg(&out).status() {
                    Ok(st) if st.success() => {}
                    Ok(st) => eprintln!("{opener} {}: {st}; open the page by hand", out.display()),
                    Err(e) => eprintln!("{opener} {}: {e}; open the page by hand", out.display()),
                }
            }
            Ok(())
        }
        "region" => {
            // A region's block of the map (PL37): its rows and columns and the range of its cells, the map as drs plots
            // map draws it (the members as their files give them, the overlap at the finest step, the map's scale).
            let f = set_file.as_ref().ok_or("region: give a plot set file (.drv)")?;
            let name = opt_value(args, "--name").ok_or("region: --name <region>")?;
            let r = f.regions.iter().find(|r| r.name == name).ok_or_else(|| format!("no region \"{name}\" in the set file ({})", if f.regions.is_empty() { "it has none".to_string() } else { f.regions.iter().map(|r| format!("\"{}\"", r.name)).collect::<Vec<_>>().join(", ") }))?;
            let reference = opt_value(args, "--reference").map(|v| v.parse::<usize>().map_err(|_| format!("--reference {v}")).and_then(|k| k.checked_sub(1).ok_or("--reference counts from 1".to_string()))).transpose()?;
            // The map on the set's map scale, the difference map on its own (the square root for a log map, P65).
            let map_scale = f.map_scale.unwrap_or(f.y_scale);
            let scale = if args.iter().any(|a| a == "--scale") { scale } else if reference.is_some() { plots::difference_scale(map_scale) } else { map_scale };
            let grid = plots::common_grid(&set, plots::GridRange::Overlap, plots::GridStep::Finest)?;
            let mp = plots::map(&set, &grid, scale, reference)?;
            let mx = plots::MapX::of(&set);
            let b = plots::map_region(&mp, r, &mx)?;
            if json {
                return print_json(&b);
            }
            let q = mx.quantity.map_or("x", |q| q.word());
            println!("# region \"{}\" on the {}map, {} scale, x in {q}", r.name, if reference.is_some() { "difference " } else { "" }, mp.scale.word());
            // Counted from 1 here, as the page counts them; the indices (--json, PL38's page data) count from 0.
            let (r0, r1, c0, c1) = (b.rows[0], b.rows[b.rows.len() - 1], b.columns[0], b.columns[b.columns.len() - 1]);
            println!("rows {}–{} of {} ({} = {} to {} {}), columns {}–{} of {} ({q} {} to {}); counted from 1, the indices from 0 being rows {r0}–{r1}, columns {c0}–{c1}", r0 + 1, r1 + 1, mp.rows.len(), set.variable_name, mp.rows[r0].value, mp.rows[r1].value, set.variable_unit, c0 + 1, c1 + 1, mp.x.len(), mp.x[c0], mp.x[c1]);
            println!("range {} to {} (the block's cells)", b.range.0, b.range.1);
            // The region's own range for the map shown, and the other's said (P65).
            let (own, other) = if reference.is_some() { ((&r.difference_range, "the difference map"), (&r.map_range, "the map")) } else { ((&r.map_range, "the map"), (&r.difference_range, "the difference map")) };
            if let Some(m) = own.0 {
                println!("its colour range on {} {} to {} ({} scale)", own.1, m.from, m.to, mp.scale.word());
                if let Some(w) = plots::range_warning("its colour range", m, &mp) {
                    eprintln!("{w}");
                }
            }
            if let Some(m) = other.0 {
                println!("its colour range on {} {} to {} (not this map's: shown with {})", other.1, m.from, m.to, if reference.is_some() { "no --reference" } else { "--reference" });
            }
            Ok(())
        }
        "grid" | "map" | "derive" => {
            let range = match opt_value(args, "--range").as_deref() {
                None | Some("overlap") => plots::GridRange::Overlap,
                Some("union") => plots::GridRange::Union,
                Some(o) => return Err(format!("--range {o}: overlap or union")),
            };
            let gstep = match opt_value(args, "--step").as_deref() {
                None | Some("finest") => plots::GridStep::Finest,
                Some("coarsest") => plots::GridStep::Coarsest,
                Some(v) => plots::GridStep::Given(num(v, "--step")?),
            };
            let grid = plots::common_grid(&set, range, gstep)?;
            match sub {
                "grid" => {
                    let res = plots::on_grid(&set, &grid);
                    if json {
                        #[cfg(feature = "serde")]
                        {
                            println!("{}", serde_json::to_string(&(&grid, &res)).map_err(|e| e.to_string())?);
                            return Ok(());
                        }
                    }
                    println!("# {} points from {} to {}; the members share one grid: {}", grid.len(), grid[0], grid[grid.len() - 1], if plots::same_grid(&set) { "yes" } else { "no" });
                    for (m, r) in set.members.iter().zip(&res) {
                        println!("{:24} {:?}, {} of {} points", m.label, r.how, r.y.iter().filter(|v| v.is_some()).count(), grid.len());
                    }
                    Ok(())
                }
                "map" => {
                    let reference = opt_value(args, "--reference").map(|v| v.parse::<usize>().map_err(|_| format!("--reference {v}")).and_then(|k| k.checked_sub(1).ok_or("--reference counts from 1".to_string()))).transpose()?;
                    let reference = reference.or_else(|| set_file.as_ref().and_then(|f| f.reference_in(&set)));
                    let mp = plots::map(&set, &grid, scale, reference)?;
                    // The set file's clip (PL37), said when it lies wholly outside the map's values.
                    let clip = set_file.as_ref().and_then(|f| if reference.is_some() { f.difference_range.map(|r| ("the difference map's colour range", r)) } else { f.map_range.map(|r| ("the map's colour range", r)) });
                    if let Some(w) = clip.and_then(|(what, r)| plots::range_warning(what, &r, &mp)) {
                        eprintln!("{w}");
                    }
                    if json {
                        #[cfg(feature = "serde")]
                        {
                            println!("{}", serde_json::to_string(&mp).map_err(|e| e.to_string())?);
                            return Ok(());
                        }
                    }
                    println!("# map: {} rows by {} columns on the {:?} scale; colour range {} to {}{}", mp.rows.len(), mp.x.len(), mp.scale, mp.range.0, mp.range.1, mp.reference.map_or(String::new(), |r| format!("; less {}", set.members[r].label)));
                    for r in &mp.rows {
                        println!("{:24} {} = {} {}, row {} to {}", set.members[r.member].label, set.variable_name, r.value, set.variable_unit, r.edges.0, r.edges.1);
                    }
                    for k in &mp.apart {
                        println!("{:24} apart: no value of {}", set.members[*k].label, set.variable_name);
                    }
                    Ok(())
                }
                _ => {
                    let flags = [("--sum", plots::Derivation::Sum), ("--average", plots::Derivation::Average), ("--weighted-average", plots::Derivation::WeightedAverage), ("--difference", plots::Derivation::Difference)];
                    let chosen: Vec<plots::Derivation> = flags.iter().filter(|f| args.iter().any(|a| a == f.0)).map(|f| f.1).collect();
                    let [derivation] = chosen[..] else { return Err("give one of --sum, --average, --weighted-average, --difference".into()) };
                    let members: Vec<usize> = opt_value(args, "--members").ok_or("--members <k,…>")?.split(',').map(|v| v.trim().parse::<usize>().ok().and_then(|k| k.checked_sub(1)).ok_or(format!("--members: {v}"))).collect::<Result<_, _>>()?;
                    // As the view shows it (PL25), with lines saying how.
                    let (set, notes) = match norm_rule(args)? {
                        Some(rule) => set.as_shown(&rule)?,
                        None => (set.clone(), Vec::new()),
                    };
                    let d = plots::derive(&set, &members, derivation, &grid)?;
                    let name = opt_value(args, "--name").or_else(|| set.title.clone()).unwrap_or_else(|| "plot set".into());
                    let out = match opt_value(args, "--out") {
                        Some(o) => o,
                        None => drs::formats::files::path("derived_data", &[("set", &name), ("operation", derivation.word())])?,
                    };
                    if !out.to_lowercase().ends_with(".xys") {
                        return Err(format!("{out}: a derived pattern is written as XYDATA, a .xys file"));
                    }
                    let path = std::path::Path::new(&out);
                    if path.exists() {
                        return Err(format!("{out} exists: nothing written"));
                    }
                    if let Some(dir) = path.parent().filter(|d| !d.as_os_str().is_empty()) {
                        std::fs::create_dir_all(dir).map_err(|e| format!("{}: {e}", dir.display()))?;
                    }
                    std::fs::write(path, plots::derived_text_with(&set, &d, &name, &notes)).map_err(|e| format!("{out}: {e}"))?;
                    println!("{out}: the {} of {} members, {} points{}", derivation.word(), members.len(), d.pattern.len(), if d.dropped > 0 { format!(" ({} left out where a member has none)", d.dropped) } else { String::new() });
                    Ok(())
                }
            }
        }
        o => Err(format!("plots: unknown subcommand '{o}' (open, save, describe, region, show, export, step, normalise, view, displace, project, windows, grid, map, derive, propose, seq)")),
    }
}

/// Whether the Mac app is installed (PL38): Spotlight knows its bundle identifier, or it is in /Applications. Never on
/// another system.
fn app_installed() -> bool {
    if !cfg!(target_os = "macos") {
        return false;
    }
    let found = std::process::Command::new("mdfind").arg("kMDItemCFBundleIdentifier == 'eus.ehu.diffrefstruct'").output().is_ok_and(|o| o.status.success() && !String::from_utf8_lossy(&o.stdout).trim().is_empty());
    found || std::path::Path::new("/Applications/diff-ref-struct.app").exists()
}

fn cmd_trend(args: &[String]) -> Result<(), String> {
    if wants_help(args) || args.len() < 2 {
        print!("{USAGE}");
        return Ok(());
    }
    let address = &args[0];
    let (mut base, mut items) = (None::<std::path::PathBuf>, Vec::new());
    let mut i = 1;
    while i < args.len() {
        if args[i] == "--base" {
            base = Some(args.get(i + 1).ok_or("--base needs a folder")?.into());
            i += 2;
            continue;
        }
        let (c, r) = match args[i].split_once(':') {
            Some((c, r)) => (std::path::PathBuf::from(c), Some(std::path::PathBuf::from(r))),
            None => (std::path::PathBuf::from(&args[i]), None),
        };
        items.push((c, r));
        i += 1;
    }
    println!("# {address} across {} version{}", items.len(), if items.len() == 1 { "" } else { "s" });
    for ((c, _), p) in items.iter().zip(drs::trend::trend(&items, base.as_deref(), address)) {
        match p {
            Ok(p) => println!("{}\t{}\t{}", c.display(), p.value, p.su.map_or("-".to_string(), |s| s.to_string())),
            Err(e) => println!("{}\terror: {e}", c.display()),
        }
    }
    Ok(())
}

/// A run for the diagnostics: the model, its data and parameter set from a control file or a `.pcr`, refined n cycles.
fn diagnosed_run(file: &str, cycles: usize) -> Result<(drs::powder::Model, drs::data::Pattern, drs::refinable::ParameterSet), String> {
    let path = std::path::Path::new(file);
    let (mut model, data, set) = if file.to_lowercase().ends_with(".pcr") {
        let text: String = std::fs::read(path).map_err(|e| format!("{file}: {e}"))?.iter().map(|&b| b as char).collect();
        let p = drs::pcr::read(&text).map_err(|e| format!("{file}: {e:?}"))?;
        let dir = path.parent().map(|d| d.to_path_buf()).unwrap_or_default();
        let stem = path.file_stem().map(|s| s.to_string_lossy().to_string()).unwrap_or_default();
        let l = drs::powder::Model::from_pcr(&p, &dir, &stem)?;
        let set = drs::refinable::ParameterSet::from_pcr(&p, &l.model)?;
        (l.model, l.data.ok_or(format!("{file}: no data file"))?, set)
    } else {
        let o = open_control(path)?;
        let set = o.control.refinement.build(&o.control.model)?;
        (o.control.model, o.data.ok_or(format!("{file}: no data"))?, set)
    };
    if cycles > 0 {
        drs::lsq::refine(&mut model, &data, &set, drs::lsq::Options { max_cycles: cycles, eps: 0.1 }, |_| {});
    }
    Ok((model, data, set))
}

fn cmd_plan(args: &[String]) -> Result<(), String> {
    if wants_help(args) || args.is_empty() {
        print!("{USAGE}");
        return Ok(());
    }
    let file = &args[0];
    let to = opt_value(args, "--to");
    let path = std::path::Path::new(file);
    let mut control = if file.to_lowercase().ends_with(".pcr") {
        let text: String = std::fs::read(path).map_err(|e| format!("{file}: {e}"))?.iter().map(|&b| b as char).collect();
        let p = drs::pcr::read(&text).map_err(|e| format!("{file}: {e:?}"))?;
        let dir = path.parent().map(|d| d.to_path_buf()).unwrap_or_default();
        let stem = path.file_stem().map(|s| s.to_string_lossy().to_string()).unwrap_or_default();
        let out_dir = to.as_ref().and_then(|t| std::path::Path::new(t).parent().map(|d| d.to_path_buf())).unwrap_or_else(|| dir.clone());
        drs::control::Control::from_pcr(&p, &dir, &stem, &out_dir)?
    } else {
        open_control(path)?.control
    };
    let steps = drs::guide::documented_plan(&control.model);
    for s in &steps {
        println!("{} ({}; {} values, {} cycles{})", s.step.name, s.source, s.step.free.len(), s.step.cycles.map_or("as the refinement's".to_string(), |c| c.to_string()), s.step.eps.map_or(String::new(), |e| format!(", every shift below {e} esd")));
        if !s.adds.is_empty() {
            println!("  adds  {}", s.adds.join(", "));
        }
        if !s.drops.is_empty() {
            println!("  drops {}", s.drops.join(", "));
        }
    }
    if let Some(t) = to {
        control.refinement = drs::guide::guided(&control.model, &control.refinement);
        std::fs::write(&t, control.to_text()).map_err(|e| format!("{t}: {e}"))?;
        println!("written {t}: {} refined parameters, a plan of {} steps", control.refinement.parameters.len(), control.refinement.plan.len());
    }
    Ok(())
}

fn cmd_pawley(args: &[String]) -> Result<(), String> {
    if wants_help(args) || args.is_empty() {
        print!("{USAGE}");
        return Ok(());
    }
    let file = &args[0];
    let to = opt_value(args, "--to").ok_or("pawley: --to <file.drc> names the control file to write")?;
    let phase = opt_value(args, "--phase").map_or(Ok(1), |v| v.parse::<usize>().map_err(|_| format!("--phase {v}: a phase number")))?;
    let iterations = opt_value(args, "--iterations").map_or(Ok(drs::pawley::START_ITERATIONS), |v| v.parse::<usize>().map_err(|_| format!("--iterations {v}: a number")))?;
    let path = std::path::Path::new(file);
    let (mut control, data) = if file.to_lowercase().ends_with(".pcr") {
        let text: String = std::fs::read(path).map_err(|e| format!("{file}: {e}"))?.iter().map(|&b| b as char).collect();
        let p = drs::pcr::read(&text).map_err(|e| format!("{file}: {e:?}"))?;
        let dir = path.parent().map(|d| d.to_path_buf()).unwrap_or_default();
        let stem = path.file_stem().map(|s| s.to_string_lossy().to_string()).unwrap_or_default();
        let out_dir = std::path::Path::new(&to).parent().map(|d| d.to_path_buf()).unwrap_or_else(|| dir.clone());
        let data = drs::powder::Model::from_pcr(&p, &dir, &stem)?.data.ok_or(format!("{file}: no data file"))?;
        (drs::control::Control::from_pcr(&p, &dir, &stem, &out_dir)?, data)
    } else {
        let o = open_control(path)?;
        (o.control, o.data.ok_or(format!("{file}: no data"))?)
    };
    if phase == 0 || phase > control.model.phases.len() {
        return Err(format!("--phase {phase}: the file has {} phases", control.model.phases.len()));
    }
    let s = drs::pawley::setup(&mut control.model, &data, &mut control.refinement, phase - 1, iterations)?;
    println!("phase {phase} ({}): {} reflections, {} intensities refined, from {iterations} Le Bail partitions", control.model.phases[phase - 1].name, s.reflections, s.parameters);
    for t in &s.ties {
        println!("  tied, equal |F|²: {}", t.iter().map(|h| format!("{} {} {}", h[0], h[1], h[2])).collect::<Vec<_>>().join(" = "));
    }
    std::fs::write(&to, control.to_text()).map_err(|e| format!("{to}: {e}"))?;
    println!("written {to}: {} refined parameters", control.refinement.parameters.len());
    Ok(())
}

fn cmd_sequential(args: &[String]) -> Result<(), String> {
    use drs::formats::files;
    if wants_help(args) || args.len() < 2 {
        print!("{USAGE}");
        return Ok(());
    }
    let file = &args[0];
    let (mut each, mut project, mut time, mut data) = (false, None::<String>, None::<String>, Vec::new());
    let mut leave_out = false;
    let mut i = 1;
    while i < args.len() {
        let next = |i: usize| args.get(i + 1).cloned().ok_or(format!("{} needs a value", args[i]));
        match args[i].as_str() {
            "--each" => (each, i) = (true, i + 1),
            "--leave-out-unmodelled" => (leave_out, i) = (true, i + 1),
            "--project" => (project, i) = (Some(next(i)?), i + 2),
            "--time" => (time, i) = (Some(next(i)?), i + 2),
            a if a.starts_with("--") => return Err(format!("sequential: unknown option {a}")),
            a => (data, i) = ({ data.push(a.to_string()); data }, i + 1),
        }
    }
    let path = std::path::Path::new(file);
    let base = path.parent().map(|d| d.to_path_buf()).unwrap_or_default();
    let stem = path.file_stem().map(|s| s.to_string_lossy().to_string()).unwrap_or_default();
    let folder = project.map(std::path::PathBuf::from).unwrap_or_else(|| base.clone());
    let control = if file.to_lowercase().ends_with(".pcr") {
        let text: String = std::fs::read(path).map_err(|e| format!("{file}: {e}"))?.iter().map(|&b| b as char).collect();
        let mut p = drs::pcr::read(&text).map_err(|e| format!("{file}: {e:?}"))?;
        if leave_out {
            for s in p.leave_out_unmodelled() {
                println!("left out: {s}");
            }
        }
        drs::control::Control::from_pcr(&p, &base, &stem, &folder)?
    } else {
        open_control(path)?.control
    };
    let cwd = std::env::current_dir().map_err(|e| e.to_string())?;
    let inputs: Vec<drs::sequential::Input> = data.iter().map(|d| drs::sequential::read(&control, &base, &cwd.join(d))).collect::<Result<_, _>>()?;
    println!("{} patterns, each from {}", inputs.len(), if each { "the file's values" } else { "the values the one before ended at" });
    let start = if each { drs::sequential::Start::Each } else { drs::sequential::Start::Previous };
    let points = drs::sequential::run(&control, &inputs, start, |_, _| !INTERRUPTED.load(std::sync::atomic::Ordering::SeqCst))?;
    for (k, p) in points.iter().enumerate() {
        let s = &p.refinement.stats;
        println!("{:4} {:20} {:>10} {:3} cycles  χ² {:9.4}  Rwp {:6.2}  {}", k + 1, p.label, p.parameter.map_or("".into(), |t| format!("{t:.3}")), p.refinement.cycles, s.chi2, s.rwp, p.refinement.stop);
    }
    let time = time.unwrap_or_else(files::now_text);
    let mut written = Vec::new();
    for (kind, body) in [("series_table", drs::sequential::series_csv(&points)), ("series_cif", drs::sequential::series_cif(&control, &points, &format!("{stem} {time}"))?)] {
        let rel = files::path(kind, &[("control", &stem), ("time", &time)])?;
        let full = folder.join(&rel);
        if let Some(d) = full.parent() {
            std::fs::create_dir_all(d).map_err(|e| format!("{}: {e}", d.display()))?;
        }
        std::fs::write(&full, body).map_err(|e| format!("{}: {e}", full.display()))?;
        written.push(rel);
    }
    println!("written in {}:", if folder.as_os_str().is_empty() { std::path::Path::new(".") } else { folder.as_path() }.display());
    for w in written {
        println!("  {w}");
    }
    Ok(())
}

fn cmd_parametric(args: &[String]) -> Result<(), String> {
    if wants_help(args) || args.len() < 2 {
        print!("{USAGE}");
        return Ok(());
    }
    let text = std::fs::read_to_string(&args[0]).map_err(|e| format!("{}: {e}", args[0]))?;
    let column = &args[1];
    let equation = opt_value(args, "--equation").ok_or("parametric: --equation names the equation to fit")?;
    let variable = opt_value(args, "--variable").unwrap_or_else(|| "T".into());
    let x = opt_value(args, "--x").unwrap_or_else(|| "parameter".into());
    let mut start = Vec::new();
    if let Some(s) = opt_value(args, "--start") {
        for part in s.split(',').filter(|p| !p.trim().is_empty()) {
            let (n, v) = part.split_once('=').ok_or(format!("--start {part}: name=value"))?;
            start.push((n.trim().to_string(), num(v.trim(), "--start")?));
        }
    }
    let points = drs::parametric::column(&text, column, &x)?;
    let eq = drs::parametric::Equation::parse(&equation, &variable)?;
    let f = drs::parametric::fit(&eq, &points, &start)?;
    println!("{column} = {equation}: {} points, {} cycles, χ² {:.6}, χ²_ν {:.6}", f.points, f.cycles, f.chi2, f.chi2_reduced);
    for (n, v, su) in &f.coefficients {
        println!("  {n:12} {v:16.9e} ± {su:.4e}");
    }
    Ok(())
}

fn cmd_locate(args: &[String]) -> Result<(), String> {
    if wants_help(args) || args.is_empty() {
        print!("{USAGE}");
        return Ok(());
    }
    let text = std::fs::read_to_string(&args[0]).map_err(|e| format!("{}: {e}", args[0]))?;
    let dir = std::path::Path::new(&args[0]).parent().map(|d| d.to_path_buf());
    let dir = dir.as_deref().filter(|d| !d.as_os_str().is_empty()).or(Some(std::path::Path::new(".")));
    use drs::control::locate;
    match (opt_value(args, "--at"), args.get(1).filter(|a| !a.starts_with("--"))) {
        (Some(b), _) => {
            let byte = b.parse::<usize>().map_err(|_| format!("--at {b}: a byte offset"))?;
            println!("{}", locate::address_at(&text, byte, dir).unwrap_or_else(|| format!("no value written at byte {byte}")));
        }
        (None, Some(a)) => match locate::locate_in(&text, a, dir) {
            Some((line, r)) => println!("line {line}, bytes {}..{}: {}", r.start, r.end, &text[r.clone()]),
            None => println!("{a}: not written in {}", args[0]),
        },
        (None, None) => {
            for (a, line, r) in locate::all(&text, dir) {
                println!("{a}\tline {line}\tbytes {}..{}\t{}", r.start, r.end, &text[r.clone()]);
            }
        }
    }
    Ok(())
}

fn cmd_impact(args: &[String]) -> Result<(), String> {
    if wants_help(args) || args.is_empty() {
        print!("{USAGE}");
        return Ok(());
    }
    let (mut cycles, mut one, mut top, mut json) = (0usize, None::<String>, usize::MAX, false);
    let mut i = 1;
    while i < args.len() {
        let next = |i: usize| args.get(i + 1).cloned().ok_or(format!("{} needs a value", args[i]));
        match args[i].as_str() {
            "--cycles" => (cycles, i) = (num(&next(i)?, "--cycles")? as usize, i + 2),
            "--parameter" => (one, i) = (Some(next(i)?), i + 2),
            "--top" => (top, i) = (num(&next(i)?, "--top")? as usize, i + 2),
            "--json" => (json, i) = (true, i + 1),
            a => return Err(format!("impact: unknown option {a}")),
        }
    }
    let (model, data, set) = diagnosed_run(&args[0], cycles)?;
    let list = match &one {
        Some(name) => vec![drs::impact::impact_of(&model, &data, &set, name)?],
        None => drs::impact::impacts(&model, &data, &set)?,
    };
    let list: Vec<_> = list.into_iter().take(top).collect();
    if json {
        return print_json(&list);
    }
    println!("# parameter  refined  value  d(chi2)/dp  shift  esds  chi2_fall  largest_change  at_2theta  per_esd  ranges");
    for p in &list {
        let ranges = p.regions.iter().map(|r| format!("{:.2}-{:.2}", r.0, r.1)).collect::<Vec<_>>();
        let shown = if ranges.len() > 4 { format!("{} … ({} ranges)", ranges[..4].join(" "), ranges.len()) } else { ranges.join(" ") };
        println!("{} {} {} {:.4e} {:+.4e} {:+.3} {:.4e} {:.4} {:.3} {:.4} {}", p.name, if p.refined { "yes" } else { "no" }, p.value, p.slope, p.shift, p.esds, p.gain, p.largest, p.at, p.per_esd, shown);
    }
    if let Some(p) = list.first().filter(|_| one.is_some()) {
        println!("# 2theta  change of Ycal for the shift");
        for (x, c) in data.x.iter().zip(&p.change) {
            println!("{x} {c}");
        }
    }
    Ok(())
}

fn cmd_normal(args: &[String]) -> Result<(), String> {
    if wants_help(args) || args.is_empty() {
        print!("{USAGE}");
        return Ok(());
    }
    let (mut cycles, mut points, mut json) = (0usize, false, false);
    let mut i = 1;
    while i < args.len() {
        let next = |i: usize| args.get(i + 1).cloned().ok_or(format!("{} needs a value", args[i]));
        match args[i].as_str() {
            "--cycles" => (cycles, i) = (num(&next(i)?, "--cycles")? as usize, i + 2),
            "--points" => (points, i) = (true, i + 1),
            "--json" => (json, i) = (true, i + 1),
            a => return Err(format!("normal: unknown option {a}")),
        }
    }
    let (model, data, _) = diagnosed_run(&args[0], cycles)?;
    let yc = model.calculate(&data.x).ycalc;
    let np = drs::impact::normal_probability(&model, &data, &yc);
    if json {
        return print_json(&np);
    }
    println!("{} points; over the central {} (expected value within ±1): slope {:.4}, intercept {:.4}", np.residuals.len(), np.central, np.slope, np.intercept);
    println!("slope 1 and intercept 0: normal residuals of the stated sigma; the slope is the goodness of fit of the well-fitted points, the intercept a bias");
    if points {
        println!("# expected  residual");
        for (e, r) in np.expected.iter().zip(&np.residuals) {
            println!("{e} {r}");
        }
    }
    Ok(())
}

fn is_multipattern_file(file: &str) -> bool {
    let lower = file.to_lowercase();
    let text = || std::fs::read(file).map(|b| b.iter().map(|&c| c as char).collect::<String>());
    (lower.ends_with(".pcr") && text().is_ok_and(|t| drs::pcr::is_multipattern(&t))) || (lower.ends_with(".drc") && text().is_ok_and(|t| drs::joint::is_joint_control(&t)))
}

fn cmd_joint(args: &[String]) -> Result<(), String> {
    use drs::formats::files;
    use drs::joint::{Esd, Joint};
    if wants_help(args) || args.is_empty() {
        print!("{USAGE}");
        return Ok(());
    }
    let file = &args[0];
    if file.to_lowercase().ends_with(".cif") {
        // A joint results file: each pattern's control file, its model read without the data (Publish's P50).
        let path = std::path::Path::new(file);
        let text = std::fs::read_to_string(path).map_err(|e| format!("{file}: {e}"))?;
        let r = drs::joint::read_results(&text).map_err(|e| format!("{file}: {e}"))?;
        let only = opt_value(args, "--pattern").map(|k| num(&k, "--pattern").map(|k| k as usize)).transpose()?;
        let cs = r.pattern_controls(path.parent().unwrap_or(std::path::Path::new(".")))?;
        for (k, c) in cs.iter().enumerate().filter(|(k, _)| only.is_none_or(|o| o == k + 1)) {
            let p = &r.patterns[k];
            let phases: Vec<String> = p.phases.iter().map(|g| g.to_string()).collect();
            println!("# pattern {} ({}), the joint phases {}", p.id, p.label, phases.join(", "));
            print!("{}", c.to_text());
        }
        return Ok(());
    }
    let (mut cycles, mut project, mut time, mut esd, mut fullprof, mut write) = (None::<usize>, None::<String>, None::<String>, Esd::Default, false, false);
    let mut save = None::<String>;
    let mut i = 1;
    while i < args.len() {
        let next = |i: usize| args.get(i + 1).cloned().ok_or(format!("{} needs a value", args[i]));
        match args[i].as_str() {
            "--cycles" => (cycles, i) = (Some(num(&next(i)?, "--cycles")? as usize), i + 2),
            "--project" => (project, write, i) = (Some(next(i)?), true, i + 2),
            "--time" => (time, i) = (Some(next(i)?), i + 2),
            "--results" => (write, i) = (true, i + 1),
            "--fullprof" => (fullprof, i) = (true, i + 1),
            "--save" => (save, i) = (Some(next(i)?), i + 2),
            "--esd" => {
                esd = match next(i)?.as_str() {
                    "default" => Esd::Default,
                    "fullprof" => Esd::FullProf,
                    e => return Err(format!("--esd {e}: default or fullprof")),
                };
                i += 2;
            }
            a => return Err(format!("joint: unknown option {a}")),
        }
    }
    let path = std::path::Path::new(file);
    let stem = path.file_stem().map(|s| s.to_string_lossy().to_string()).unwrap_or_default();
    let dir = path.parent().map(|d| d.to_path_buf()).unwrap_or_default();
    let text: String = std::fs::read(path).map_err(|e| format!("{file}: {e}"))?.iter().map(|&b| b as char).collect();
    let mut j = if drs::joint::is_joint_control(&text) { Joint::open(path)? } else { Joint::from_pcr(&text, &dir, &stem, fullprof)? };
    println!("{}: {} patterns, {} phases, {} parameters ({} shared)", file, j.patterns.len(), j.phases.len(), j.parameters.len(), j.parameters.iter().filter(|q| q.members.iter().map(|m| m.0).collect::<std::collections::BTreeSet<_>>().len() > 1).count());
    for (k, p) in j.patterns.iter().enumerate() {
        println!("  pattern {}: {} ({}, λ {}), weight {}, {}; phases {:?}", k + 1, p.file, if p.model.neutrons { "neutrons" } else { "X-rays" }, p.model.lambda[0], p.weight, if p.refined { "refined" } else { "calculated only" }, p.phases.iter().map(|g| g + 1).collect::<Vec<_>>());
    }
    let opts = drs::lsq::Options::stated(cycles.or(j.cycles), j.eps);
    let r = j.refine(opts, esd, |c| {
        println!("cycle {:>3}: global χ² {:.4} (from {:.4}), Q {:.3}, λ {:.0e}", c.cycle, c.chi2, c.chi2_start, c.objective, c.marquardt);
        true
    });
    println!("{}", r.stop);
    for (k, p) in r.statistics.patterns.iter().enumerate() {
        println!("  pattern {}: Rp {:.2} %, Rwp {:.2} %, Rexp {:.2} %, χ² {:.4} (Bragg points {:.4})", k + 1, p.all.rp, p.all.rwp, p.all.rexp, p.all.chi2, p.bragg.chi2);
    }
    println!("  joint: Rp {:.2} %, Rwp {:.2} %, Rexp {:.2} %, χ² {:.4}, global χ² (Bragg points) {:.4}", r.statistics.rp, r.statistics.rwp, r.statistics.rexp, r.statistics.chi2_all, r.statistics.chi2);
    for (n, v, e) in &r.params {
        println!("  {n:<44} {v:>14.7} ({e:.2e})");
    }
    if write {
        let folder = project.map(std::path::PathBuf::from).unwrap_or_else(|| dir.clone());
        let time = time.unwrap_or_else(files::now_text);
        let body = j.results(Some(&r), &format!("{stem} {time}"), &text)?;
        let rel = files::path("results", &[("control", &stem), ("time", &time)])?;
        let full = folder.join(&rel);
        if let Some(d) = full.parent() {
            std::fs::create_dir_all(d).map_err(|e| format!("{}: {e}", d.display()))?;
        }
        std::fs::write(&full, body).map_err(|e| format!("{}: {e}", full.display()))?;
        println!("written in {}: {rel}", folder.display());
    }
    if let Some(to) = save {
        for p in j.save(std::path::Path::new(&to))? {
            println!("written: {}", p.display());
        }
    }
    Ok(())
}

fn cmd_results(args: &[String]) -> Result<(), String> {
    use drs::formats::files;
    if wants_help(args) || args.is_empty() {
        print!("{USAGE}");
        return Ok(());
    }
    let file = &args[0];
    if is_multipattern_file(file) {
        // A multi-pattern .pcr: its joint refinement's results file (P06, P41).
        let mut a = args.to_vec();
        a.push("--results".into());
        return cmd_joint(&a);
    }
    let (mut cycles, mut project, mut time, mut exports, mut plan) = (0usize, None::<String>, None::<String>, false, false);
    let mut i = 1;
    while i < args.len() {
        let next = |i: usize| args.get(i + 1).cloned().ok_or(format!("{} needs a value", args[i]));
        match args[i].as_str() {
            "--cycles" => (cycles, i) = (num(&next(i)?, "--cycles")? as usize, i + 2),
            "--project" => (project, i) = (Some(next(i)?), i + 2),
            "--time" => (time, i) = (Some(next(i)?), i + 2),
            "--exports" => (exports, i) = (true, i + 1),
            "--plan" => (plan, i) = (true, i + 1),
            a => return Err(format!("results: unknown option {a}")),
        }
    }
    let path = std::path::Path::new(file);
    let stem = path.file_stem().map(|s| s.to_string_lossy().to_string()).unwrap_or_default();
    let folder = project.map(std::path::PathBuf::from).unwrap_or_else(|| path.parent().map(|d| d.to_path_buf()).unwrap_or_default());
    let (mut control, data) = if file.to_lowercase().ends_with(".pcr") {
        let text: String = std::fs::read(path).map_err(|e| format!("{file}: {e}"))?.iter().map(|&b| b as char).collect();
        let p = drs::pcr::read(&text).map_err(|e| format!("{file}: {e:?}"))?;
        let dir = path.parent().map(|d| d.to_path_buf()).unwrap_or_default();
        let c = drs::control::Control::from_pcr(&p, &dir, &stem, &folder)?;
        let data = drs::powder::Model::from_pcr(&p, &dir, &stem)?.data.ok_or(format!("{file}: no data file"))?;
        (c, data)
    } else {
        let o = open_control(path)?;
        (o.control, o.data.ok_or(format!("{file}: no data"))?)
    };
    // The refinement as written (n cycles), or its plan run step by step (P37): each cycle then records its step.
    let (r, set) = if plan {
        let rf = control.refinement.clone();
        if rf.plan.is_empty() {
            return Err(format!("{file} has no plan ([[refine.step]] tables)"));
        }
        let steps: Vec<usize> = (0..rf.plan.len()).collect();
        let (r, session) = drs::lsq::run_plan(&mut control.model, &data, &rf, &steps, |_| {}, |_| true)?;
        (r, session.set)
    } else {
        let set = control.refinement.build(&control.model)?;
        let r = drs::lsq::refine(&mut control.model, &data, &set, drs::lsq::Options { max_cycles: cycles, eps: 0.1 }, |_| {});
        (r, set)
    };
    let record = &r.history;
    let time = time.unwrap_or_else(files::now_text);
    let (text, structures) = drs::results::write(&control, &data, Some(&set), Some(&r), &format!("{stem} {time}"))?;
    // the user's notes on the reflections (W8), from the notes file beside the control file
    let notes_file = path.parent().unwrap_or(std::path::Path::new(".")).join(files::path("reflection_notes", &[("control", &stem)])?);
    let notes = match std::fs::read_to_string(&notes_file) {
        Ok(t) => drs::notes::read(&t).map_err(|e| format!("{}: {e}", notes_file.display()))?,
        Err(_) => Vec::new(),
    };
    let text = drs::results::with_notes(&text, &notes);
    let mut written = Vec::new();
    let mut put = |rel: String, body: &str| -> Result<(), String> {
        let full = folder.join(&rel);
        if let Some(d) = full.parent() {
            std::fs::create_dir_all(d).map_err(|e| format!("{}: {e}", d.display()))?;
        }
        std::fs::write(&full, body).map_err(|e| format!("{}: {e}", full.display()))?;
        written.push(rel);
        Ok(())
    };
    put(files::path("results", &[("control", &stem), ("time", &time)])?, &text)?;
    // The summary and the listing (P57), to read beside the record.
    let summary = drs::listing::summary(&control.model, &data, &set, (r.cycles > 0).then_some(&r));
    let title = format!("{} ({}), run {time}", path.display(), control.title.trim());
    put(files::path("run_summary", &[("control", &stem), ("time", &time)])?, &drs::listing::summary_text(&title, &summary))?;
    put(files::path("run_listing", &[("control", &stem), ("time", &time)])?, &drs::listing::listing_text(&title, record, &summary))?;
    for (phase, cif) in &structures {
        let name: String = phase.split_whitespace().collect::<Vec<_>>().join(" ").chars().map(|c| if c == '/' || c == '\\' { '_' } else { c }).collect();
        put(files::path("deposition_cif", &[("control", &stem), ("phase", &name)])?, cif)?;
    }
    if exports {
        put(files::path("export_reflections", &[("control", &stem)])?, &drs::results::reflections_csv(&drs::diagnostics::reflections(&control.model, &data)))?;
        put(files::path("export_profile", &[("control", &stem)])?, &drs::results::profile_csv(&control.model, &data))?;
        // A Le Bail phase's extracted intensities, in FullProf's CODFILn_cltr.int (P35).
        for (k, ph) in control.model.phases.iter().enumerate().filter(|(_, p)| p.extract) {
            let rows = drs::lebail::extracted(&control.model, &data, k)?;
            let name: String = ph.name.split_whitespace().collect::<Vec<_>>().join(" ").chars().map(|c| if c == '/' || c == '\\' { '_' } else { c }).collect();
            let text = drs::lebail::int_file(&format!("Phase No.: {:3} {}", k + 1, ph.name), control.model.lambda[0], &rows, drs::lebail::CLUSTER_DEFAULTS);
            put(files::path("export_intensities", &[("control", &stem), ("phase", &name)])?, &text)?;
        }
    }
    println!("{} cycles, Rwp {:.2} %, chi2 {:.4}; written in {}:", r.cycles, r.stats.rwp, r.stats.chi2, folder.display());
    for w in written {
        println!("  {w}");
    }
    Ok(())
}

fn cmd_rebase(args: &[String]) -> Result<(), String> {
    use drs::check::Severity;
    if wants_help(args) || args.len() < 2 {
        print!("{USAGE}");
        return Ok(());
    }
    let (file, data) = (std::path::Path::new(&args[0]), std::path::Path::new(&args[1]));
    let (mut to, mut reader) = (None, None);
    let mut i = 2;
    while i < args.len() {
        let v = args.get(i + 1).ok_or(format!("{} needs a value", args[i]))?.clone();
        match args[i].as_str() {
            "--to" => to = Some(v),
            "--as" => reader = Some(v),
            a => return Err(format!("rebase: unknown option {a}")),
        }
        i += 2;
    }
    let to = to.ok_or("rebase: --to <new.drc>")?;
    let r = drs::rebase::rebase_path(file, data, reader.as_deref(), std::path::Path::new(&to))?;
    for d in &r.diagnostics {
        println!("{}", d.to_text(&to));
    }
    let count = |s: Severity| r.diagnostics.iter().filter(|d| d.severity == s).count();
    println!("{to}: written; {} error(s), {} warning(s), {} hint(s)", count(Severity::Error), count(Severity::Warning), count(Severity::Hint));
    if drs::check::has_errors(&r.diagnostics) {
        return Err(format!("{to}: {} error(s)", count(Severity::Error)));
    }
    Ok(())
}

/// The normalisation `--by …` asks for (PL6, PL25), none without `--by`.
fn norm_rule(args: &[String]) -> Result<Option<drs::plots::Normalisation>, String> {
    use drs::plots;
    let range = |what: &str| -> Result<(f64, f64), String> {
        let v = opt_value(args, "--within").ok_or(format!("--by {what} needs --within <from:to>"))?;
        let (a, b) = v.split_once(':').ok_or("--within <from:to>")?;
        Ok((num(a, "--within")?, num(b, "--within")?))
    };
    let to = opt_value(args, "--to").map(|v| num(&v, "--to")).transpose()?.unwrap_or(1.0);
    let Some(by) = opt_value(args, "--by") else { return Ok(None) };
    let rule = match by.as_str() {
        "monitor" => plots::Normalisation::Monitor { to },
        "time" => plots::Normalisation::Time { to },
        "maximum" => match opt_value(args, "--within") {
            Some(_) => {
                let (from, until) = range("maximum")?;
                plots::Normalisation::MaximumIn { from, until, to }
            }
            None => plots::Normalisation::Maximum { to },
        },
        "area" => {
            let (from, until) = range("area")?;
            plots::Normalisation::Area { from, until, to }
        }
        "peak" => {
            let (from, until) = range("peak")?;
            let background = match opt_value(args, "--background").as_deref() {
                None | Some("auto") => plots::WindowBackground::Auto { log_lambda: None },
                Some("line") => plots::WindowBackground::Line,
                Some(o) => return Err(format!("--background {o}: auto or line")),
            };
            plots::Normalisation::Peak { window: plots::Window { name: "peak".into(), from, to: until }, background, to }
        }
        "reference" => {
            let (from, until) = range("reference")?;
            let k = opt_value(args, "--member").ok_or("--by reference needs --member <k>")?;
            let member = k.parse::<usize>().ok().and_then(|k| k.checked_sub(1)).ok_or(format!("--member {k}: from 1"))?;
            plots::Normalisation::Reference { member, from, until }
        }
        o => return Err(format!("--by {o}: monitor, time, maximum, area, peak or reference")),
    };
    Ok(Some(rule))
}

fn cmd_peak(args: &[String]) -> Result<(), String> {
    use drs::peakfit::{fit, PeakRequest, PeakShape};
    if wants_help(args) || args.len() < 3 {
        print!("{USAGE}");
        return Ok(());
    }
    let (from, to) = (num(&args[1], "from")?, num(&args[2], "to")?);
    let shape = if args.iter().any(|a| a == "--gaussian") {
        PeakShape::Gaussian
    } else if args.iter().any(|a| a == "--lorentzian") {
        PeakShape::Lorentzian
    } else {
        PeakShape::PseudoVoigt
    };
    let d = drs::formats::read_path(std::path::Path::new(&args[0]), None).map_err(|e| format!("{}: {e}", args[0]))?;
    let scan = d.scans.first().ok_or("no scan in the file")?;
    let doublet = match args.iter().position(|a| a == "--doublet") {
        // The scan's own λ₁, λ₂ and ratio (the app's request 38).
        Some(i) if args.get(i + 1).map(String::as_str) == Some("scan") => {
            let w = scan.wavelength.ok_or("--doublet scan: the file states no wavelength")?;
            Some(drs::peakfit::doublet_of(&w).ok_or("--doublet scan: the file states no second wavelength with its ratio")?)
        }
        Some(i) => {
            let v = |k: usize| args.get(i + k).ok_or("--doublet <λ1> <λ2> <ratio>".to_string()).and_then(|a| num(a, "--doublet"));
            Some((v(2)? / v(1)?, v(3)?))
        }
        None => None,
    };
    let pat = &scan.pattern;
    let r = fit(&pat.x, &pat.y, &pat.sigma, &PeakRequest { from, to, shape, doublet })?;
    let show = |e: drs::peakfit::Estimate, p: usize| match e.esd {
        Some(s) => drs::cif_write::with_su(e.value, s, p),
        None => format!("{} (held)", e.value),
    };
    println!("{} points between {from} and {to}: χ² {:.4}, {} iterations{}", r.points, r.chi2, r.iterations, if r.converged { "" } else { ", not converged" });
    println!("  position   {}°", show(r.position, 6));
    println!("  intensity  {}", show(r.intensity, 3));
    println!("  FWHM       {}°", show(r.fwhm, 6));
    println!("  η          {}", show(r.eta, 4));
    println!("  background {} + {} (2θ − {:.4})", show(r.background[0], 3), show(r.background[1], 3), r.centre);
    Ok(())
}

fn cmd_name(args: &[String]) -> Result<(), String> {
    if wants_help(args) || args.len() < 2 {
        print!("{USAGE}");
        return Ok(());
    }
    let path = std::path::Path::new(&args[0]);
    let model = if args[0].to_lowercase().ends_with(".pcr") {
        let text: String = std::fs::read(path).map_err(|e| format!("{}: {e}", args[0]))?.iter().map(|&b| b as char).collect();
        let p = drs::pcr::read(&text).map_err(|e| format!("{}: {e:?}", args[0]))?;
        let stem = path.file_stem().map(|s| s.to_string_lossy().to_string()).unwrap_or_default();
        drs::powder::Model::from_pcr(&p, path.parent().unwrap_or(std::path::Path::new(".")), &stem)?.model
    } else {
        open_control(path)?.control.model
    };
    let name = &args[1];
    let phase = opt_value(args, "--phase").map(|v| v.parse::<usize>().map_err(|_| format!("--phase {v}"))).transpose()?.unwrap_or(1).saturating_sub(1);
    let label = opt_value(args, "--label");
    let address = if name.starts_with('_') { drs::output_names::address_of_cif(name, phase, label.as_deref(), &model) } else { drs::output_names::address_of_fullprof(name, &model) };
    match address {
        Some(a) => {
            println!("{a}");
            Ok(())
        }
        None => Err(format!("{name}: no parameter of the model")),
    }
}

fn cmd_complete(args: &[String]) -> Result<(), String> {
    if wants_help(args) || args.len() < 3 {
        print!("{USAGE}");
        return Ok(());
    }
    let text = std::fs::read_to_string(&args[0]).map_err(|e| format!("{}: {e}", args[0]))?;
    let line: usize = args[1].parse().map_err(|_| format!("line {}", args[1]))?;
    let column: usize = args[2].parse().map_err(|_| format!("column {}", args[2]))?;
    for c in drs::control::complete::complete(&text, line, column) {
        println!("{:<28} {:<8} {}", c.label, format!("{:?}", c.kind).to_lowercase(), c.detail);
    }
    Ok(())
}

fn cmd_decompose(args: &[String]) -> Result<(), String> {
    use drs::mode_decomposition as md;
    if wants_help(args) || args.is_empty() {
        print!("{USAGE}");
        return Ok(());
    }
    let read = |f: &str| -> Result<String, String> { std::fs::read(f).map(|b| b.iter().map(|&c| c as char).collect()).map_err(|e| format!("{f}: {e}")) };
    let analysis: usize = opt_value(args, "--analysis").map_or(Ok(1), |v| v.parse().map_err(|_| format!("--analysis {v}")))?;
    let mut summary = Vec::new();
    let mut amplimodes = None;
    let (parent, distorted, p, origin) = if let Some(f) = opt_value(args, "--amplimodes") {
        let all = drs::modes::read_amplimodes(&read(&f)?);
        let a = all.get(analysis.wrapping_sub(1)).ok_or(format!("{f}: no analysis {analysis} ({} in the file)", all.len()))?.clone();
        let (parent, mut distorted, p, o) = md::from_amplimodes(&a)?;
        if args.iter().any(|x| x == "--mapped") {
            distorted.atoms = a.mappings.iter().map(|m| md::Atom { label: m.2.clone(), species: m.2.chars().take_while(|c| c.is_ascii_alphabetic()).collect(), position: m.3 }).collect();
        }
        summary = a.summary.clone();
        amplimodes = Some((a, all));
        (parent, distorted, p, o)
    } else {
        let files: Vec<&String> = args.iter().enumerate().filter(|(i, a)| !a.starts_with("--") && (*i == 0 || !matches!(args[i - 1].as_str(), "--basis" | "--origin" | "--analysis"))).map(|(_, a)| a).collect();
        let [pf, df] = files[..] else { return Err(format!("give the parent's and the distorted structure's CIF\n\n{USAGE}")) };
        let basis = opt_value(args, "--basis").ok_or("--basis: the distorted phase's basis in the parent's, as \"a-b,a+b,c\"")?;
        let origin = match opt_value(args, "--origin") {
            Some(v) => {
                let o: Vec<f64> = v.split(',').map(|t| drs::cif::number(t.trim()).or_else(|| t.trim().split_once('/').and_then(|(a, b)| Some(a.parse::<f64>().ok()? / b.parse::<f64>().ok()?)))).collect::<Option<_>>().ok_or("--origin: x,y,z")?;
                <[f64; 3]>::try_from(o).map_err(|_| "--origin: x,y,z")?
            }
            None => [0.0; 3],
        };
        (md::from_cif(&read(pf)?)?, md::from_cif(&read(df)?)?, md::parse_basis(&basis)?, origin)
    };
    let mut d = md::decompose(&parent, &distorted, &p, origin)?;
    // The irreps named after AMPLIMODES's modes, when the analysis has its block for FullProf.
    // (the analysis's own, else that of another analysis of the file with the same parent and transformation).
    let named = amplimodes.as_ref().and_then(|(a, all): &(drs::modes::Analysis, Vec<drs::modes::Analysis>)| a.modes().or_else(|| all.iter().find(|b| b.fullprof.is_some() && b.transformation == a.transformation && b.high.as_ref().map(|h| h.group) == a.high.as_ref().map(|h| h.group)).and_then(|b| b.modes())));
    if let Some(m) = named {
        let n = d.name_by_modes(&m);
        println!("{n} irreps named after the analysis's modes");
    }
    let c = &d.reference_cell;
    println!("reference cell {:.6} {:.6} {:.6} {:.4} {:.4} {:.4}; supercell {}×{}×{} of the parent's, |G/nZ³| {}, |H| {}, {} primitive cells of the distorted phase", c.a, c.b, c.c, c.alpha, c.beta, c.gamma, d.supercell, d.supercell, d.supercell, d.order, d.order_h, d.primitive_cells);
    println!("reference structure:");
    for ((l, sp, x, _), u) in d.reference.iter().zip(&d.displacements) {
        println!("  {l:8} {sp:3} {:9.6} {:9.6} {:9.6}   u {:9.6} {:9.6} {:9.6}", x[0], x[1], x[2], u[0], u[1], u[2]);
    }
    println!("global distortion {:.4} Å, largest displacement {:.4} Å, not kept by H {:.1e} Å, translation removed {:?} Å", d.global, d.largest, d.residual, d.translation.map(|v| (v * 1e6).round() / 1e6));
    for r in &d.irreps {
        let star: Vec<String> = r.star.iter().map(|k| format!("({},{},{})", k[0], k[1], k[2])).collect();
        println!("  {:14} k {:28} parity {:2} dim {} {:?}  modes {}  {:.4} Å", r.name, star.join(" "), r.parity.map_or("".into(), |p| if p > 0 { "+".to_string() } else { "-".to_string() }), r.dimension, r.reality, r.modes, r.amplitude);
        for o in &r.orbits {
            println!("      {:8} modes {}  {:.4} Å  ({})", parent.atoms[o.orbit].label, o.modes.len(), o.amplitude, o.modes.iter().map(|m| format!("{:.6}", m.amplitude)).collect::<Vec<_>>().join(", "));
        }
    }
    if !summary.is_empty() {
        println!("AMPLIMODES:");
        for s in &summary {
            println!("  {:8} {:8} {:16} modes {}  {:.4} Å", s.k, s.irrep, s.direction, s.dimension, s.amplitude);
        }
    }
    if args.iter().any(|x| x == "--fullprof") {
        print!("{}", d.fullprof_block());
    }
    Ok(())
}

fn cmd_isodistort(args: &[String]) -> Result<(), String> {
    if wants_help(args) || args.is_empty() {
        print!("{USAGE}");
        return Ok(());
    }
    let text: String = std::fs::read(&args[0]).map_err(|e| format!("{}: {e}", args[0]))?.iter().map(|&c| c as char).collect();
    let iso = drs::isodistort::isodistort_cif(&text)?;
    println!("{} modes, {} deltas, {} matrix entries", iso.modes.len(), iso.deltas.len(), iso.matrix.len());
    for m in &iso.modes {
        println!("  {:>3} {:10} {:8} {:14} {:6} {:2} {:4} {:10.5}{}", m.id, m.k, m.irrep, m.direction, m.atom, m.wyckoff, m.kind, m.value, m.norm.map_or(String::new(), |n| format!("  norm {n:.5}")));
    }
    for (irrep, ms) in iso.irreps() {
        let amp = ms.iter().map(|&i| iso.modes[i].value.powi(2)).sum::<f64>().sqrt();
        println!("  {irrep}: modes {}; amplitude {amp:.5}", ms.iter().map(|&i| iso.modes[i].id.to_string()).collect::<Vec<_>>().join(" "));
    }
    let from_modes = iso.deltas_from_modes();
    for (id, label, v) in &iso.deltas {
        let m = from_modes.iter().find(|x| x.0 == *id).map_or(0.0, |x| x.1);
        println!("  {label:10} {v:10.5}  from the modes {m:10.5}");
    }
    Ok(())
}
