"""
Option C: Validate anuga.scenario setup_* functions using the real Towradgi mesh.

This script is intentionally a *smoke test*, not a quantitative calibration check.
It loads the existing Towradgi mesh (towradgi.tsh), applies scenario setup functions,
and evolves for a short time to confirm that the pipeline runs end-to-end without
errors on a realistic, large-scale mesh (~250,000 triangles).

Skip conditions
---------------
Skipped automatically (pytest.skip) when:
  - The ``DEM_bridges/towradgi.tsh`` mesh file is not present (data not downloaded).
  - The ``anuga.scenario`` module cannot be imported.

Running the test
----------------
From the repository root::

    pytest validation_tests/case_studies/towradgi/validate_scenario.py -v

Or directly from the case study directory::

    cd validation_tests/case_studies/towradgi
    python -m pytest validate_scenario.py -v

Data download
-------------
If the mesh file is absent, download Towradgi data first::

    cd validation_tests/case_studies/towradgi
    python data_download.py
"""

import os
import sys
import shutil
import tempfile

import pytest

# ---------------------------------------------------------------------------
# Skip guards
# ---------------------------------------------------------------------------

_CASE_DIR = os.path.dirname(os.path.abspath(__file__))
_MESH_FILE = os.path.join(_CASE_DIR, 'DEM_bridges', 'towradgi.tsh')

_DATA_AVAILABLE = os.path.exists(_MESH_FILE)
_DATA_SKIP = pytest.mark.skipif(
    not _DATA_AVAILABLE,
    reason='Towradgi mesh not found. Run data_download.py to get it.')

try:
    from anuga.scenario.setup_initial_conditions import setup_initial_conditions
    from anuga.scenario.setup_boundary_conditions import setup_boundary_conditions
    from anuga.scenario.setup_rainfall import setup_rainfall
    from anuga.scenario.setup_riverwalls import setup_riverwalls
    import anuga
    _HAS_MODULE = True
    _MODULE_SKIP_REASON = ''
except ImportError as _e:
    _HAS_MODULE = False
    _MODULE_SKIP_REASON = str(_e)

_MODULE_SKIP = pytest.mark.skipif(not _HAS_MODULE, reason=_MODULE_SKIP_REASON)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_rain_csv(path, times, rates_mm_hr):
    with open(path, 'w') as f:
        f.write('time,rain_mm_hr\n')
        for t, r in zip(times, rates_mm_hr):
            f.write('%g,%g\n' % (t, r))


def _make_project(rain_csv=None, output_dir=None, rw_csv_files=None,
                  elevation=0.0, friction=0.04):
    """Minimal project-like object for scenario setup functions."""
    from unittest.mock import MagicMock
    import anuga.utilities.spatialInputUtil as su

    p = MagicMock()

    for qty, val in [('elevation', elevation), ('friction', friction),
                     ('stage', 0.0), ('xmomentum', 0.0), ('ymomentum', 0.0)]:
        setattr(p, f'{qty}_data',       [['All', val]])
        setattr(p, f'{qty}_clip_range', [[float('-inf'), float('inf')]])
        setattr(p, f'{qty}_mean',       None)
        setattr(p, f'{qty}_additions',  [])

    if rain_csv is not None:
        p.rain_data = [[rain_csv, 0.0, 'linear', 'All', 1.0]]
    else:
        p.rain_data = []

    # Riverwalls: read real Towradgi riverwall CSV files if provided
    if rw_csv_files:
        p.riverwalls, p.riverwall_par = su.readListOfRiverWalls(rw_csv_files)
    else:
        p.riverwalls = {}
        p.riverwall_par = {}

    p.output_dir = output_dir or tempfile.mkdtemp()
    p.spatial_text_output_dir = 'SPATIAL_TEXT'
    return p


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

@_DATA_SKIP
@_MODULE_SKIP
class TestTowradgiScenario:

    @pytest.fixture(autouse=True)
    def _setup(self):
        self.tmpdir = tempfile.mkdtemp()
        yield
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _load_domain(self):
        """Load the Towradgi mesh and configure domain basics."""
        domain = anuga.create_domain_from_file(_MESH_FILE)
        domain.set_name('validate_scenario')
        domain.set_datadir(self.tmpdir)
        domain.set_flow_algorithm('DE0')
        return domain

    def test_setup_initial_conditions_on_real_mesh(self):
        """setup_initial_conditions runs without error on ~250k-triangle mesh."""
        domain = self._load_domain()
        project = _make_project(elevation=0.0, friction=0.04)
        setup_initial_conditions(domain, project)
        import numpy as np
        elev = domain.get_quantity('elevation').get_values(location='centroids')
        assert len(elev) == domain.get_number_of_triangles()

    def test_setup_boundary_conditions_reflective(self):
        """All boundary tags get a Reflective boundary without error."""
        domain = self._load_domain()
        setup_initial_conditions(domain, _make_project(elevation=0.0))
        all_tags = list(set(domain.boundary.values()))
        project = _make_project()
        project.boundary_data = [[t, 'Reflective'] for t in all_tags]
        project.boundary_tags = all_tags
        setup_boundary_conditions(domain, project)

    def test_setup_rainfall_creates_operators(self):
        """setup_rainfall attaches Rate_operator(s) to the real mesh."""
        domain = self._load_domain()
        setup_initial_conditions(domain, _make_project(elevation=0.0))
        rain_csv = os.path.join(self.tmpdir, 'rain.csv')
        _write_rain_csv(rain_csv, [0, 3600], [20.0, 20.0])
        project = _make_project(rain_csv=rain_csv, output_dir=self.tmpdir)
        setup_rainfall(domain, project)
        from anuga.operators.rate_operators import Rate_operator
        ops = [op for op in domain.fractional_step_operators
               if isinstance(op, Rate_operator)]
        assert len(ops) == 1

    def test_setup_riverwalls_with_real_csvs(self):
        """setup_riverwalls builds riverwall data from real Towradgi CSV files."""
        rw_dir = os.path.join(_CASE_DIR, 'Model', 'Riverwalls')
        rw_files = sorted([
            os.path.join(rw_dir, f)
            for f in os.listdir(rw_dir) if f.endswith('.csv')
        ])
        if not rw_files:
            pytest.skip('No riverwall CSV files found in Model/Riverwalls/')
        domain = self._load_domain()
        setup_initial_conditions(domain, _make_project(elevation=0.0))
        project = _make_project(
            rw_csv_files=rw_files, output_dir=self.tmpdir)
        os.makedirs(os.path.join(self.tmpdir, 'SPATIAL_TEXT'), exist_ok=True)
        setup_riverwalls(domain, project)

    def test_short_evolve_with_rainfall(self):
        """Entire setup pipeline + 30s evolve runs without error."""
        domain = self._load_domain()
        setup_initial_conditions(domain, _make_project(elevation=0.0,
                                                       friction=0.04))
        all_tags = list(set(domain.boundary.values()))
        project_bc = _make_project()
        project_bc.boundary_data = [[t, 'Reflective'] for t in all_tags]
        project_bc.boundary_tags = all_tags
        setup_boundary_conditions(domain, project_bc)
        rain_csv = os.path.join(self.tmpdir, 'rain.csv')
        _write_rain_csv(rain_csv, [0, 100], [10.0, 10.0])
        project_rain = _make_project(rain_csv=rain_csv,
                                     output_dir=self.tmpdir)
        setup_rainfall(domain, project_rain)
        for _ in domain.evolve(yieldstep=10.0, finaltime=30.0):
            pass


if __name__ == '__main__':
    sys.exit(pytest.main([__file__, '-v']))
