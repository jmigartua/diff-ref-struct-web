"""diff-ref-struct from Python (P53, Q1/L021; L380): the core's scripting interface.

Two ways in, both the core's own code; nothing is calculated in Python:

* ``drs.run(["refine", "y2o3.pcr", "--cycles", "3"])`` runs any command of the ``drs`` command line, the very binary
  built with this package, and returns its exit code and output. Every command and option of ``drs --help`` works.
* ``drs.open(path)`` gives a :class:`Project` (a control file ``.drc`` or a FullProf ``.pcr`` with its data) whose
  methods call the core in this process: values by their addresses, the calculated pattern, a refinement, the state,
  the checks, the ``.pcr`` export and the results file. Results are plain dicts and lists (``numpy.asarray`` takes
  them as they are).

>>> import drs
>>> p = drs.open("y2o3_structure.pcr")
>>> p.get("phase1.cell_a")
>>> r = p.refine(cycles=3)
>>> r["stats"]["chi2"], p.get("phase1.cell_a")
"""

import json as _json
import os as _os
import subprocess as _subprocess
from typing import Dict, List, NamedTuple, Optional, Sequence, Tuple

from . import _drs
from ._drs import DrsError

__version__ = _drs.version()

__all__ = ["DrsError", "Project", "Run", "binary", "open", "read_data", "run", "__version__"]


def binary() -> str:
    """The path of the ``drs`` command line inside this package."""
    path = _os.path.join(_os.path.dirname(__file__), "bin", "drs.exe" if _os.name == "nt" else "drs")
    if not _os.path.exists(path):
        raise DrsError(f"the drs command line is not in the package ({path}): build the wheel with build_wheel.sh")
    return path


class Run(NamedTuple):
    """What a command printed: its exit code, standard output and standard error."""

    code: int
    stdout: str
    stderr: str


def run(args: Sequence[str], check: bool = True, cwd: Optional[str] = None) -> Run:
    """Run ``drs`` with ``args`` (the words after ``drs`` on the command line) and return its output.

    With ``check`` (the default) a failing command raises :class:`DrsError` with its message.
    """
    p = _subprocess.run([binary(), *[str(a) for a in args]], capture_output=True, text=True, cwd=cwd)
    out = Run(p.returncode, p.stdout, p.stderr)
    if check and p.returncode != 0:
        raise DrsError(p.stderr.strip() or f"drs {' '.join(map(str, args))}: exit code {p.returncode}")
    return out


def read_data(path: str, format: Optional[str] = None) -> Dict:
    """A data file read by the core's readers: ``{"format": …, "scans": [{"label", "axis", "x", "y", "sigma"}, …]}``.

    ``format``: a reader's id (``drs data formats``), else the core identifies the file.
    """
    return _json.loads(_drs.read_data(str(path), format))


class Project:
    """A refinement project: a control file (``.drc``) or a ``.pcr`` imported into the core's default model.

    ``base``: the folder a control file's data and other files are named from (``--base``).
    """

    def __init__(self, path: str, base: Optional[str] = None):
        self._p = _drs.Project(str(path), None if base is None else str(base))

    @property
    def title(self) -> str:
        return self._p.title()

    def control_text(self) -> str:
        """The control file's text, as the core writes it."""
        return self._p.control_text()

    def save(self, path: str) -> None:
        """Write the control file."""
        self._p.save(str(path))

    def fields(self) -> List[Dict]:
        """Every value of the model: ``address``, dictionary ``id``, ``value``, ``lower``, ``upper``, ``fixed_by``."""
        return _json.loads(self._p.fields_json())

    def get(self, address: str) -> float:
        """A value by its address (``phase1.cell_a``, ``phase1.O1.atom_position.x``, …)."""
        return self._p.get(address)

    def set(self, address: str, value: float) -> None:
        """Set a value by its address; the values tied to it follow (a cubic cell's b and c, a site's coordinates)."""
        self._p.set(address, float(value))

    def free(self) -> List[Dict]:
        """The refinement's free parameters: ``name``, ``value`` and the ``fields`` each moves, ``[address, coefficient]``."""
        return _json.loads(self._p.free_json())

    def calculate(self) -> Dict:
        """The pattern at the model's values: ``x``, ``yobs``, ``sigma``, ``ycalc``, ``background``, ``phases`` (each
        phase's curve) and ``stats`` (``rp``, ``rwp``, ``rexp``, ``chi2``, ``durbin_watson``, ``n``, ``p``, ``c``)."""
        return _json.loads(self._p.calculate_json())

    def refine(self, cycles: Optional[int] = None, eps: Optional[float] = None) -> Dict:
        """Refine (up to the file's cycles, else 10; convergence the file's, else 0.1); the project's model moves to the
        values reached. Returns ``cycles``, ``converged``, ``stop``, ``stats`` and ``params`` (``name``, ``value``,
        ``esd``)."""
        return _json.loads(self._p.refine_json(cycles, eps))

    def point(self, values: Dict[str, float]) -> Dict:
        """The statistics with free parameters at the values given, as a refinement step sets them; the model is left
        as it is. Returns ``ycalc`` and ``stats``."""
        pairs: List[Tuple[str, float]] = [(str(k), float(v)) for k, v in values.items()]
        return _json.loads(self._p.point_json(pairs))

    def check(self) -> List[Dict]:
        """The checks ``drs check`` makes of the control file: ``severity``, ``level``, ``code``, ``what``, ``why``,
        ``suggestion``, ``addresses``, ``line``."""
        return _json.loads(self._p.check_json())

    def state(self) -> Dict:
        """The state (``drs state``): ``stats`` and per phase its ``cell`` (a, b, c, α, β, γ, V as [value, su]),
        ``contacts`` and ``valence`` sums, with the su's of the last refinement."""
        return _json.loads(self._p.state_json())

    def to_pcr(self) -> Dict:
        """The ``.pcr`` FullProf reads (``pcr``), the ``files`` to write beside it and the export's ``notes``."""
        return _json.loads(self._p.pcr_json())

    def results(self, name: str = "drs") -> str:
        """The results file (pdCIF) of the last refinement, or of the model's values before any."""
        return self._p.results(name)


def open(path: str, base: Optional[str] = None) -> Project:  # noqa: A001 (as Python's own open, for files)
    """Open a control file or a ``.pcr`` as a :class:`Project`."""
    return Project(path, base)
