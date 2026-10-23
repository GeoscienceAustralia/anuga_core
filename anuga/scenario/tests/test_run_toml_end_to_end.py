"""
End-to-end smoke test for the ``anuga_toml_run`` runner.

Unlike the other scenario tests (which exercise the ``setup_*`` modules in
isolation against synthetic domains), this drives the *whole* runner the way a
user does: it lays down the tiny ``simple`` dam-break scenario in a temporary
directory and invokes ``anuga_toml_run`` on it as a subprocess, then validates
the SWW it produces.

To run regardless of the working directory, it is deliberately self-contained:

* **Runner** — preferred is the installed ``anuga_toml_run`` console command
  (``shutil.which``); falling back to ``scripts/anuga_toml_run.py`` in a source
  checkout. The test only skips if neither can be found.
* **Inputs** — the shipped ``examples/run_toml/simple/`` files are used when a
  checkout is locatable (so the committed example is smoke-tested); otherwise
  the equivalent inputs are written inline, so the test still runs from any
  directory against an installed ANUGA.

The scenario is tiny (two CSV polygons + a short TOML) and runs ~30 s of model
time on ~500 triangles, so it takes a few seconds. Marked slow because it builds
a mesh, evolves, and writes GeoTIFFs.
"""

import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from glob import glob
from pathlib import Path

import pytest

try:
    import numpy as np
    from anuga.file.netcdf import NetCDFFile
    HAS_MODULE = True
    SKIP_REASON = ''
except ImportError as _e:  # pragma: no cover - import guard
    HAS_MODULE = False
    SKIP_REASON = str(_e)

# The anuga_toml_run runner imports the scenario boundary/mesh setup, which pulls
# in the geodata interface (fiona/rasterio/shapely) via spatialInputUtil. When
# that stack is unavailable (e.g. a broken conda geodata install in CI) the
# runner subprocess fails; skip rather than report a spurious failure, matching
# the other geodata-dependent tests.
try:
    import fiona  # noqa: F401
    import rasterio  # noqa: F401
    import shapely  # noqa: F401
    HAS_GEODATA = True
except ImportError:
    HAS_GEODATA = False


# ---------------------------------------------------------------------------
# Inline copy of examples/run_toml/simple/ (used when no checkout is found)
# ---------------------------------------------------------------------------

_DAM_BREAK_TOML = """\
[project]
scenario               = "dam_break"
output_base_directory  = "OUTPUT/"
yieldstep              = 2.0
finaltime              = 30.0
projection_information = -56
flow_algorithm         = "DE1"

[mesh]
bounding_polygon = "bounding_polygon.csv"
default_res      = 20.0
[[mesh.boundary_tags]]
tag   = "south"
edges = [0]
[[mesh.boundary_tags]]
tag   = "east"
edges = [1]
[[mesh.boundary_tags]]
tag   = "north"
edges = [2]
[[mesh.boundary_tags]]
tag   = "west"
edges = [3]

[boundary_conditions]
[[boundary_conditions.boundaries]]
tag  = "south"
type = "Reflective"
[[boundary_conditions.boundaries]]
tag  = "east"
type = "Reflective"
[[boundary_conditions.boundaries]]
tag  = "north"
type = "Reflective"
[[boundary_conditions.boundaries]]
tag  = "west"
type = "Reflective"

[initial_conditions]
[[initial_conditions.elevation]]
polygon = "All"
value   = 0.0
[[initial_conditions.friction]]
polygon = "All"
value   = 0.03
[[initial_conditions.stage]]
polygon = "dam.csv"
value   = 4.0
[[initial_conditions.stage]]
polygon = "All"
value   = 0.0
[[initial_conditions.xmomentum]]
polygon = "All"
value   = 0.0
[[initial_conditions.ymomentum]]
polygon = "All"
value   = 0.0
"""

_BOUNDING_POLYGON_CSV = "0,0\n100,0\n100,100\n0,100\n"
_DAM_CSV = "0,0\n40,0\n40,100\n0,100\n"


# ---------------------------------------------------------------------------
# Locate the runner (installed command or source script) and example data
# ---------------------------------------------------------------------------

def _search_roots():
    """Candidate checkout roots, walking up from cwd, this file, and an optional
    ``ANUGA_SOURCE_ROOT``. cwd is the reliable anchor because an installed
    (copied) ANUGA imports from site-packages, so ``__file__`` does not reach
    ``scripts/``/``examples/``."""
    starts = [Path.cwd().resolve(), Path(__file__).resolve()]
    env_root = os.environ.get('ANUGA_SOURCE_ROOT')
    if env_root:
        starts.insert(0, Path(env_root).resolve() / '_')  # parents include env_root
    seen = set()
    for start in starts:
        for parent in start.parents:
            if parent not in seen:
                seen.add(parent)
                yield parent


def _find_runner():
    """Path to the runner: the installed ``anuga_toml_run`` console command if on
    PATH, else ``scripts/anuga_toml_run.py`` from a checkout, else None."""
    installed = shutil.which('anuga_toml_run')
    if installed:
        return installed
    for root in _search_roots():
        candidate = root / 'scripts' / 'anuga_toml_run.py'
        if candidate.is_file():
            return str(candidate)
    return None


def _find_example_dir():
    """The shipped ``examples/run_toml/simple/`` directory if a checkout is
    locatable, else None (inputs are then generated inline)."""
    for root in _search_roots():
        toml = root / 'examples' / 'run_toml' / 'simple' / 'dam_break.toml'
        if toml.is_file():
            return str(toml.parent)
    return None


RUNNER = _find_runner()
EXAMPLE_DIR = _find_example_dir()


def _runner_cmd(runner):
    """A ``.py`` script is run through the interpreter; an installed console
    command is executed directly."""
    if runner.endswith('.py'):
        return [sys.executable, runner]
    return [runner]


def _stage_inputs(work):
    """Populate the work dir with the dam-break scenario: copy the shipped
    example files when available, otherwise write the inline equivalents."""
    if EXAMPLE_DIR is not None:
        for name in os.listdir(EXAMPLE_DIR):
            src = os.path.join(EXAMPLE_DIR, name)
            if os.path.isfile(src):
                shutil.copy(src, work)
        return
    (Path(work) / 'dam_break.toml').write_text(_DAM_BREAK_TOML)
    (Path(work) / 'bounding_polygon.csv').write_text(_BOUNDING_POLYGON_CSV)
    (Path(work) / 'dam.csv').write_text(_DAM_CSV)


@unittest.skipUnless(HAS_MODULE, SKIP_REASON)
@unittest.skipUnless(RUNNER is not None, 'anuga_toml_run runner not found')
@unittest.skipUnless(HAS_GEODATA,
                     'requires geodata interface (fiona, rasterio, shapely)')
class TestRunTomlEndToEnd(unittest.TestCase):

    @pytest.mark.slow
    def test_simple_dam_break_produces_valid_sww(self):
        with tempfile.TemporaryDirectory() as work:
            _stage_inputs(work)

            # Run the runner the way a user would, serially.
            env = dict(os.environ, OMP_NUM_THREADS='1')
            proc = subprocess.run(
                _runner_cmd(RUNNER) + ['dam_break.toml'],
                cwd=work, env=env,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, timeout=600)

            tail = '\n'.join(proc.stdout.splitlines()[-25:])
            self.assertEqual(
                proc.returncode, 0,
                'anuga_toml_run exited %d. Output tail:\n%s'
                % (proc.returncode, tail))

            # The runner creates OUTPUT/RUN_<timestamp>_dam_break/dam_break.sww
            matches = glob(os.path.join(
                work, 'OUTPUT', 'RUN_*_dam_break', 'dam_break.sww'))
            self.assertEqual(
                len(matches), 1,
                'expected exactly one SWW, found %r. Output tail:\n%s'
                % (matches, tail))
            sww = matches[0]
            self.assertGreater(os.path.getsize(sww), 0, 'SWW is empty')

            self._check_sww(sww)

    def _check_sww(self, sww):
        fid = NetCDFFile(sww, 'r')
        try:
            # --- time vector: finaltime 30 s, yieldstep 2 s -> 16 frames ---
            time = fid.variables['time'][:]
            self.assertEqual(len(time), 16)
            self.assertAlmostEqual(float(time[0]), 0.0, places=6)
            self.assertAlmostEqual(float(time[-1]), 30.0, places=6)

            # --- stage: finite, and the 4 m reservoir column is present ---
            stage = fid.variables['stage'][:]
            self.assertTrue(np.all(np.isfinite(stage)), 'non-finite stage')
            self.assertGreater(float(np.max(stage)), 3.5)
            self.assertLess(float(np.max(stage)), 4.5)

            # --- CRS metadata: projection_information = -56 must propagate ---
            # to the SWW as UTM zone 56 (southern) / EPSG:32756. This guards
            # the runner's set_epsg() step (commit c4ef89f5).
            self.assertEqual(int(fid.zone), 56)
            self.assertEqual(int(fid.epsg), 32756)
        finally:
            fid.close()


if __name__ == '__main__':
    unittest.main()
