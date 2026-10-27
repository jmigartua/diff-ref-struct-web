"""The Python bindings against the command line (P53, L381), on the Y2O3 tutorial.

Run with the wheel installed:  python -m unittest discover -s core/python/tests -v
(DRS_RUNS names the folder of FullProf's runs, manuals/fullprof/runs, when the tests are not run from the repository.)
"""

import os
import re
import subprocess
import tempfile
import unittest

import drs

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.environ.get("DRS_RUNS", os.path.join(HERE, "..", "..", "..", "manuals", "fullprof", "runs"))
Y2O3 = os.path.join(RUNS, "Y2O3", "y2o3_structure.pcr")


def cli(*args):
    """The command line's own output, the binary in the package run directly."""
    return subprocess.run([drs.binary(), *args], capture_output=True, text=True, check=True).stdout


class Package(unittest.TestCase):
    def test_version_is_the_command_lines(self):
        self.assertEqual(drs.run(["--version"]).stdout.strip(), f"drs {drs.__version__}")

    def test_run_is_the_command_line(self):
        # Every command through drs.run gives the binary's output to the byte.
        for args in (["cell", "5.4", "5.4", "5.4", "90", "90", "90", "--hkl", "1", "1", "1", "--lambda", "1.5406"], ["sg", "P 21/c"], ["calc", Y2O3]):
            self.assertEqual(drs.run(args).stdout, cli(*args), args)
        help_text = drs.run(["--help"]).stdout
        commands = sorted(set(re.findall(r"^  drs ([a-z]+)", help_text, re.M)))
        self.assertGreaterEqual(len(commands), 40, commands)

    def test_a_failing_command_raises(self):
        with self.assertRaises(drs.DrsError) as e:
            drs.run(["refine", "no-such-file.pcr"])
        self.assertIn("no-such-file", str(e.exception))
        self.assertNotEqual(drs.run(["unknown-command"], check=False).code, 0)


class Project(unittest.TestCase):
    def test_values_by_address(self):
        p = drs.open(Y2O3)
        fields = {f["address"]: f for f in p.fields()}
        a = fields["phase1.cell_a"]
        self.assertEqual(a["id"], "cell_a")
        self.assertEqual(p.get("phase1.cell_a"), a["value"])
        # A cubic cell: b and c follow a, and are fixed by the crystal system.
        p.set("phase1.cell_a", 10.61)
        self.assertEqual((p.get("phase1.cell_b"), p.get("phase1.cell_c")), (10.61, 10.61))
        self.assertEqual(fields["phase1.cell_b"]["fixed_by"], "crystal system")
        with self.assertRaises(drs.DrsError):
            p.set("phase1.cell_b", 10.0)
        with self.assertRaises(drs.DrsError):
            p.get("phase1.no_such_value")

    def test_calculate_is_drs_calc(self):
        # Ycal and the agreement against `drs calc`'s printed table and header.
        c = drs.open(Y2O3).calculate()
        out = cli("calc", Y2O3)
        chi2 = float(re.search(r"chi2 ([0-9.]+)", out).group(1))
        rows = [l.split() for l in out.splitlines() if l and not l.startswith("#")]
        self.assertEqual(len(rows), len(c["x"]))
        for i in (0, 1000, len(rows) - 1):
            self.assertAlmostEqual(float(rows[i][0]), c["x"][i], places=5)
            self.assertAlmostEqual(float(rows[i][1]), c["ycalc"][i], places=4)
            self.assertAlmostEqual(float(rows[i][2]), c["background"][i], places=4)
        # `drs calc` counts no parameter; with the refinement's free parameters χ² is over N − P.
        s = c["stats"]
        self.assertAlmostEqual(s["chi2"] * (s["n"] - s["p"] + s["c"]) / s["n"], chi2, places=3)
        self.assertEqual(len(c["phases"]), 1)
        self.assertAlmostEqual(c["background"][500] + c["phases"][0][500], c["ycalc"][500], places=9)

    def test_refine_is_drs_refine(self):
        p = drs.open(Y2O3)
        r = p.refine(cycles=3)
        out = cli("refine", Y2O3, "--cycles", "3")
        for prm in r["params"]:
            m = re.search(rf"^{re.escape(prm['name'])}\s+(\S+)\s+(\S+)$", out, re.M)
            self.assertIsNotNone(m, prm["name"])
            self.assertAlmostEqual(float(m.group(1)), prm["value"], places=5, msg=prm["name"])
            self.assertAlmostEqual(float(m.group(2)), prm["esd"], places=5, msg=prm["name"])
        self.assertEqual(p.get("phase1.cell_a"), next(x["value"] for x in r["params"] if x["name"] == "phase1.cell_a"))
        # The cycles run, as the command line's table counts them (it converges before the third).
        cycles = [int(m.group(1)) for m in re.finditer(r"^\s+(\d+)\s+\d+\.\d\d\s+\d+\.\d\d\s+\d+\.\d\d\s", out, re.M)]
        self.assertEqual(r["cycles"], max(cycles))
        self.assertLessEqual(r["cycles"], 3)
        # The state carries the refinement's su's: a's is its esd.
        st = p.state()
        esd = next(x["esd"] for x in r["params"] if x["name"] == "phase1.cell_a")
        self.assertAlmostEqual(st["phases"][0]["cell"][0][1], esd, delta=1e-6 * esd)
        # A point leaves the model as it is.
        a = p.get("phase1.cell_a")
        q = p.point({"phase1.cell_a": a + 1e-3})
        self.assertGreater(q["stats"]["chi2"], r["stats"]["chi2"])
        self.assertEqual(p.get("phase1.cell_a"), a)
        # The results file names the refinement's values.
        cif = p.results("y2o3 test")
        self.assertIn("_pd_proc_ls_prof_wR_factor", cif)
        self.assertIn("_drs_free_name", cif)

    def test_control_file_round_trip_and_export(self):
        # Saved elsewhere and opened again with its files named from the run's folder (--base): the same text and
        # pattern.
        p = drs.open(Y2O3)
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "y2o3.drc")
            p.save(path)
            q = drs.open(path, base=os.path.dirname(Y2O3))
            self.assertEqual(q.control_text(), p.control_text())
            self.assertEqual(q.calculate()["ycalc"], p.calculate()["ycalc"])
        e = p.to_pcr()
        self.assertIn("Space group symbol", e["pcr"])
        self.assertEqual([x for x in p.check() if x["severity"] == "error"], [])

    def test_read_data(self):
        d = drs.read_data(os.path.join(RUNS, "Y2O3", "y2o3.dat"))
        c = drs.open(Y2O3).calculate()
        self.assertEqual(len(d["scans"]), 1)
        self.assertEqual(len(d["scans"][0]["x"]), len(c["x"]))
        self.assertTrue(all(isinstance(v, float) for v in d["scans"][0]["y"][:10]))

    def test_a_project_without_values_fails_clearly(self):
        with self.assertRaises(drs.DrsError):
            drs.open(os.path.join(RUNS, "no-such-run.pcr"))


if __name__ == "__main__":
    unittest.main()
