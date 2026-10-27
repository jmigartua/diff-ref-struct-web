# diff-ref-struct from Python

The core's scripting interface (P53; Q1, L021; built in L380–L381): Python bindings to the Rust core through PyO3. No
physics is written in Python; every number comes from the core.

## Two ways in

- `drs.run(args)` runs any command of the `drs` command line. It uses the binary built from the same core, which
  travels inside the package, and returns its exit code, standard output and standard error. Every command and option
  of `drs --help` works as on the command line, `--base` and Ctrl-C included.
- `drs.open(path)` opens a control file (`.drc`) or a FullProf `.pcr` with its data as a `Project`. Its methods call the
  core in the Python process:

| method | what the core does | the command line's equivalent |
|---|---|---|
| `fields()`, `get(address)`, `set(address, value)` | the model's values by their addresses, checked; tied values follow | `drs locate`, the control file |
| `free()` | the refinement's free parameters and the fields each moves | `drs refine` (its header) |
| `calculate()` | x, Yobs, σ, Ycal, background, each phase's curve, the agreement | `drs calc` |
| `refine(cycles, eps)` | the refinement; the project's model moves to the values reached | `drs refine --cycles` |
| `point({name: value})` | the agreement at given values of free parameters (E20a) | `drs state --at` |
| `state()` | the cell with su's, contacts, valence sums (E3, E18) | `drs state` |
| `check()` | the checks of the control file | `drs check` |
| `to_pcr()` | the `.pcr` FullProf reads | `drs convert` |
| `results(name)` | the results file (pdCIF) | `drs results` |
| `save(path)`, `control_text()` | the control file | `drs convert` |

`drs.read_data(path)` reads a data file with the core's readers, as `drs data read` does.

Results are plain dicts and lists of floats, which `numpy.asarray` takes as they are. Errors raise `drs.DrsError`.

```python
import drs
p = drs.open("y2o3_structure.pcr")
r = p.refine(cycles=3)
print(r["stats"]["chi2"], p.get("phase1.cell_a"))
print(drs.run(["sg", "P 21/c"]).stdout)
```

## Build and test

You need Rust and maturin (`pip install maturin`). From `core/python`:

```sh
./build_wheel.sh            # → target/wheels/diff_ref_struct-<version>-cp39-abi3-<platform>.whl
pip install target/wheels/diff_ref_struct-*.whl
python -m unittest discover -s tests -v
```

- The wheel uses Python's stable ABI from 3.9, so one wheel serves every later Python on its platform. On the Mac, the
  same script builds the macOS wheel.
- The tests compare the bindings with the command line on the Y₂O₃ tutorial:
  - `drs.run` gives the binary's output to the byte;
  - `calculate` gives `drs calc`'s Ycal and χ²;
  - `refine` gives `drs refine`'s values, esd's and cycles.
- The package's version is the core's: `drs.__version__` equals `drs --version`.
