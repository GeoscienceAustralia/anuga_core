"""
A negative starttime (a burn-in before the forcing data begin at t=0).

Inlet, stage-boundary and rainfall CSV series start at t=0. Before that the
inlet and boundary hold their first value and the rain is zero, so a domain
can evolve from a negative starttime without the interpolators raising.
"""
import os
import shutil
import tempfile
import unittest
from unittest.mock import MagicMock

try:
    import numpy as np
    from anuga import rectangular_cross_domain
    from anuga.scenario.timeseries import make_timeseries_function
    from anuga.scenario.setup_inlets import setup_inlets
    from anuga.scenario.setup_rainfall import setup_rainfall
    from anuga.scenario.setup_boundary_conditions import \
        setup_boundary_conditions
    HAS_MODULE = True
    SKIP_REASON = ''
except ImportError as _e:
    HAS_MODULE = False
    SKIP_REASON = str(_e)


def _write(path, text):
    with open(path, 'w') as f:
        f.write(text)
    return path


@unittest.skipUnless(HAS_MODULE, SKIP_REASON)
class TestTimeseriesFunction(unittest.TestCase):

    def test_holds_first_value_before_start(self):
        f = make_timeseries_function([0.0, 10.0], [2.0, 4.0])
        self.assertEqual(f(-100.0), 2.0)
        self.assertAlmostEqual(float(f(5.0)), 3.0)

    def test_fixed_value_before_start(self):
        f = make_timeseries_function([0.0, 10.0], [2.0, 4.0], before=0.0)
        self.assertEqual(f(-1.0), 0.0)

    def test_after_end_still_raises(self):
        f = make_timeseries_function([0.0, 10.0], [2.0, 4.0])
        with self.assertRaises(ValueError):
            f(11.0)


@unittest.skipUnless(HAS_MODULE, SKIP_REASON)
class TestNegativeStarttime(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_burn_in_runs_inlet_only(self):
        """From t=-20 s to 0 only the 1 m^3/s inlet adds water; the rain
        (36 mm/hr) starts at t=0."""
        domain = rectangular_cross_domain(10, 5, len1=100.0, len2=50.0)
        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('stage', 0.0)
        domain.set_starttime(-20.0)

        inlet_csv = _write(os.path.join(self.tmp, 'inlet.csv'),
                           'time,discharge\n0.0,1.0\n1000.0,1.0\n')
        line_csv = _write(os.path.join(self.tmp, 'line.csv'),
                          '50.0,10.0\n50.0,40.0\n')
        stage_csv = _write(os.path.join(self.tmp, 'stage.csv'),
                           'time,stage\n0.0,0.0\n1000.0,0.0\n')
        rain_csv = _write(os.path.join(self.tmp, 'rain.csv'),
                          'time,rate\n0.0,36.0\n1000.0,36.0\n')

        project = MagicMock()
        project.inlet_data = [['in', line_csv, inlet_csv, 0.0]]
        project.rain_data = [[rain_csv, 0.0, 'linear', 'All']]
        project.boundary_data = [['left', 'Stage', stage_csv, 0.0],
                                 ['right', 'Reflective'],
                                 ['bottom', 'Reflective'],
                                 ['top', 'Reflective']]
        project.boundary_tags = ['left', 'right', 'bottom', 'top']

        setup_inlets(domain, project)
        setup_rainfall(domain, project)
        setup_boundary_conditions(domain, project)

        volumes = {}
        for t in domain.evolve(yieldstep=10.0, finaltime=10.0):
            volumes[round(t)] = domain.get_water_volume()

        self.assertEqual(sorted(volumes), [-20, -10, 0, 10])
        # Burn-in: inlet only, 1 m^3/s for 20 s
        self.assertAlmostEqual(volumes[0], 20.0, delta=0.5)
        # Then inlet plus rain: 10 s * (1 + 36e-3/3600 * 5000 m^2)
        self.assertAlmostEqual(volumes[10] - volumes[0], 10.5, delta=0.5)
        self.assertTrue(np.isfinite(volumes[10]))


if __name__ == '__main__':
    unittest.main()
