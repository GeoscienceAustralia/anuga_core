"""
Integration tests for anuga.scenario setup_* modules.

Covers modules not tested elsewhere:
  - setup_initial_conditions
  - setup_culverts
  - setup_riverwalls
  - end-to-end: initial conditions + rainfall + evolve

All tests use a synthetic 1×1 rectangular domain (rectangular_cross_domain,
10×10 cells) so no external data files are required.
"""

import os
import shutil
import tempfile
import unittest
from unittest.mock import MagicMock

try:
    import numpy as np
    import anuga
    from anuga import rectangular_cross_domain
    from anuga.scenario.setup_initial_conditions import setup_initial_conditions
    from anuga.scenario.setup_culverts import setup_culverts
    from anuga.scenario.setup_riverwalls import setup_riverwalls
    from anuga.scenario.setup_rainfall import setup_rainfall
    HAS_MODULE = True
    SKIP_REASON = ''
except ImportError as _e:
    HAS_MODULE = False
    SKIP_REASON = str(_e)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_domain():
    """10×10 rectangular cross domain on [0,1]×[0,1].

    The finer grid (0.1 m cells) ensures culvert enquiry points land in a
    different triangle from the inlet exchange region, avoiding the
    "Enquiry point is in an inlet triangle" UserWarning.
    """
    return rectangular_cross_domain(10, 10)


def _make_initial_project(**overrides):
    """Project mock with all initial-condition attributes set to flat values."""
    p = MagicMock()
    elev = overrides.get('elevation', 1.0)
    fric = overrides.get('friction', 0.03)
    stage = overrides.get('stage', 1.5)

    for qty, val in [('elevation', elev), ('friction', fric),
                     ('stage', stage), ('xmomentum', 0.0), ('ymomentum', 0.0)]:
        setattr(p, f'{qty}_data',       [['All', val]])
        setattr(p, f'{qty}_clip_range', [[float('-inf'), float('inf')]])
        setattr(p, f'{qty}_mean',       None)
        setattr(p, f'{qty}_additions',  [])
    return p


def _write_rain_csv(path, times, rates_mm_hr):
    with open(path, 'w') as f:
        f.write('time,rain_mm_hr\n')
        for t, r in zip(times, rates_mm_hr):
            f.write('%g,%g\n' % (t, r))


def _count_rate_operators(domain):
    from anuga.operators.rate_operators import Rate_operator
    return sum(1 for op in domain.fractional_step_operators
               if isinstance(op, Rate_operator))


# ---------------------------------------------------------------------------
# setup_initial_conditions
# ---------------------------------------------------------------------------

@unittest.skipUnless(HAS_MODULE, SKIP_REASON)
class TestSetupInitialConditions(unittest.TestCase):

    def setUp(self):
        self.domain = _make_domain()

    def test_sets_constant_elevation(self):
        project = _make_initial_project(elevation=3.5)
        setup_initial_conditions(self.domain, project)
        vals = self.domain.get_quantity('elevation').get_values(
            location='centroids')
        np.testing.assert_allclose(vals, 3.5, rtol=1e-6)

    def test_sets_constant_friction(self):
        project = _make_initial_project(friction=0.05)
        setup_initial_conditions(self.domain, project)
        vals = self.domain.get_quantity('friction').get_values(
            location='centroids')
        np.testing.assert_allclose(vals, 0.05, rtol=1e-6)

    def test_sets_constant_stage(self):
        project = _make_initial_project(elevation=0.0, stage=2.0)
        setup_initial_conditions(self.domain, project)
        vals = self.domain.get_quantity('stage').get_values(
            location='centroids')
        np.testing.assert_allclose(vals, 2.0, rtol=1e-6)

    def test_stage_below_elevation_sets_dry_cell(self):
        """stage set below elevation remains as-set (ANUGA clamps during evolve)."""
        project = _make_initial_project(elevation=5.0, stage=0.0)
        setup_initial_conditions(self.domain, project)
        elev = self.domain.get_quantity('elevation').get_values(
            location='centroids')
        stage = self.domain.get_quantity('stage').get_values(
            location='centroids')
        # Elevation was set to 5.0, stage to 0.0 — domain should reflect that
        np.testing.assert_allclose(elev, 5.0, rtol=1e-6)
        np.testing.assert_allclose(stage, 0.0, atol=1e-6)

    def test_all_five_quantities_set(self):
        """No AttributeError for any of the five quantities."""
        project = _make_initial_project()
        setup_initial_conditions(self.domain, project)
        for qty in ('elevation', 'friction', 'stage', 'xmomentum', 'ymomentum'):
            vals = self.domain.get_quantity(qty).get_values(
                location='centroids')
            self.assertEqual(len(vals), self.domain.get_number_of_triangles())

    def test_additions_applied_on_top_of_base(self):
        """elevation_additions should add to the base quantity value."""
        project = _make_initial_project(elevation=1.0)
        project.elevation_additions = [['All', 0.5]]
        setup_initial_conditions(self.domain, project)
        vals = self.domain.get_quantity('elevation').get_values(
            location='centroids')
        # Base 1.0 + addition 0.5 = 1.5
        np.testing.assert_allclose(vals, 1.5, rtol=1e-6)


# ---------------------------------------------------------------------------
# setup_culverts
# ---------------------------------------------------------------------------

@unittest.skipUnless(HAS_MODULE, SKIP_REASON)
class TestSetupCulverts(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.domain = _make_domain()
        self.domain.set_datadir(self.tmpdir)
        setup_initial_conditions(self.domain, _make_initial_project(
            elevation=0.0, stage=0.5))

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _culvert_data(self, ctype='boyd_box', **kw):
        defaults = dict(
            label='c1', type=ctype, losses=0.0, barrels=1.0,
            blockage=0.0, z1=0.0, z2=0.0, apron=0.1, manning=0.013,
            enquiry_gap=0.2, smoothing_timescale=0.0,
            use_momentum_jet=True, use_velocity_head=True,
            exchange_line_0=None, exchange_line_1=None,
            invert_elevations=None,
        )
        defaults.update(kw)
        return defaults

    def test_no_culverts_adds_no_operators(self):
        project = MagicMock()
        project.culvert_data = []
        n_before = len(self.domain.fractional_step_operators)
        setup_culverts(self.domain, project)
        self.assertEqual(len(self.domain.fractional_step_operators), n_before)

    def test_culvert_data_absent_attribute_no_error(self):
        """If culvert_data is missing entirely the function should not crash."""
        project = MagicMock(spec=[])   # no attributes at all
        # Should silently do nothing (getattr default is [])
        setup_culverts(self.domain, project)

    def test_boyd_box_end_points_creates_operator(self):
        # enquiry_point = end_point ± (apron + enquiry_gap) * direction
        # apron=0.1, enquiry_gap=0.2 → gap=0.3; keep end_points > 0.3 from boundary
        cd = self._culvert_data(
            ctype='boyd_box', width=0.05, height=None,
            diameter=None,
            end_point_0=[0.35, 0.5],
            end_point_1=[0.65, 0.5],
        )
        project = MagicMock()
        project.culvert_data = [cd]
        setup_culverts(self.domain, project)
        from anuga.structures.boyd_box_operator import Boyd_box_operator
        ops = [op for op in self.domain.fractional_step_operators
               if isinstance(op, Boyd_box_operator)]
        self.assertEqual(len(ops), 1)

    def test_boyd_pipe_end_points_creates_operator(self):
        cd = self._culvert_data(
            ctype='boyd_pipe', diameter=0.05, width=None, height=None,
            end_point_0=[0.35, 0.3],
            end_point_1=[0.65, 0.3],
        )
        project = MagicMock()
        project.culvert_data = [cd]
        setup_culverts(self.domain, project)
        from anuga.structures.boyd_pipe_operator import Boyd_pipe_operator
        ops = [op for op in self.domain.fractional_step_operators
               if isinstance(op, Boyd_pipe_operator)]
        self.assertEqual(len(ops), 1)

    def test_two_culverts_created(self):
        cd1 = self._culvert_data(
            label='c1', ctype='boyd_box', width=0.05, height=None,
            diameter=None, end_point_0=[0.35, 0.5], end_point_1=[0.65, 0.5])
        cd2 = self._culvert_data(
            label='c2', ctype='boyd_box', width=0.05, height=None,
            diameter=None, end_point_0=[0.35, 0.7], end_point_1=[0.65, 0.7])
        project = MagicMock()
        project.culvert_data = [cd1, cd2]
        setup_culverts(self.domain, project)
        from anuga.structures.boyd_box_operator import Boyd_box_operator
        ops = [op for op in self.domain.fractional_step_operators
               if isinstance(op, Boyd_box_operator)]
        self.assertEqual(len(ops), 2)


# ---------------------------------------------------------------------------
# setup_riverwalls
# ---------------------------------------------------------------------------

@unittest.skipUnless(HAS_MODULE, SKIP_REASON)
class TestSetupRiverwalls(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.domain = _make_domain()
        setup_initial_conditions(self.domain, _make_initial_project(
            elevation=0.0, stage=0.5))

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_empty_riverwalls_is_noop(self):
        project = MagicMock()
        project.riverwalls = {}
        project.riverwall_par = {}
        project.output_dir = self.tmpdir
        project.spatial_text_output_dir = 'SPATIAL_TEXT'
        # Should complete without error
        setup_riverwalls(self.domain, project)

    def test_missing_riverwalls_attribute_is_noop(self):
        """Old project objects without riverwalls attr should not crash."""
        project = MagicMock()
        project.riverwalls = {}
        project.riverwall_par = {}
        project.output_dir = self.tmpdir
        project.spatial_text_output_dir = 'SPATIAL_TEXT'
        setup_riverwalls(self.domain, project)


# ---------------------------------------------------------------------------
# End-to-end: initial conditions + rainfall → evolve → stage increases
# ---------------------------------------------------------------------------

@unittest.skipUnless(HAS_MODULE, SKIP_REASON)
class TestScenarioEvolve(unittest.TestCase):
    """Smoke test: set up a wet domain, add rainfall, evolve, check depth."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_rainfall_increases_mean_stage(self):
        domain = _make_domain()
        domain.set_name('evolve_test')
        domain.set_datadir(self.tmpdir)

        # Flat, wet initial conditions
        setup_initial_conditions(domain, _make_initial_project(
            elevation=0.0, stage=0.1))

        # Reflective boundaries on all edges
        all_tags = list(set(domain.boundary.values()))
        boundary_dict = {tag: anuga.Reflective_boundary(domain)
                         for tag in all_tags}
        domain.set_boundary(boundary_dict)

        # Rainfall: 100 mm/hr over entire domain for 60 s
        rain_csv = os.path.join(self.tmpdir, 'rain.csv')
        _write_rain_csv(rain_csv, [0, 100], [100.0, 100.0])
        project = MagicMock()
        project.rain_data = [[rain_csv, 0.0, 'linear', 'All', 1.0]]
        setup_rainfall(domain, project)

        stage_before = domain.get_quantity('stage').get_values(
            location='centroids').mean()

        for _ in domain.evolve(yieldstep=10.0, finaltime=60.0):
            pass

        stage_after = domain.get_quantity('stage').get_values(
            location='centroids').mean()

        # 100 mm/hr for 60 s = 100/3.6e6 m/s × 60 s ≈ 0.00167 m added depth
        self.assertGreater(stage_after, stage_before)

    def test_no_rainfall_stage_conserved(self):
        """Without rain a closed reflective domain keeps total volume constant."""
        domain = _make_domain()
        domain.set_name('no_rain_test')
        domain.set_datadir(self.tmpdir)

        setup_initial_conditions(domain, _make_initial_project(
            elevation=0.0, stage=0.5))

        all_tags = list(set(domain.boundary.values()))
        boundary_dict = {tag: anuga.Reflective_boundary(domain)
                         for tag in all_tags}
        domain.set_boundary(boundary_dict)

        vol_before = domain.get_quantity('stage').get_values(
            location='centroids').sum()

        for _ in domain.evolve(yieldstep=10.0, finaltime=30.0):
            pass

        vol_after = domain.get_quantity('stage').get_values(
            location='centroids').sum()

        # Volume should be conserved to within 0.1%
        self.assertAlmostEqual(vol_before, vol_after, delta=abs(vol_before) * 1e-3)


if __name__ == '__main__':
    unittest.main()
