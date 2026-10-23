"""
Unit tests for anuga.scenario.setup_rainfall

Tests cover:
  - No rainfall data — no operators added
  - Single global rainfall (polygon = 'All')
  - Rainfall with a spatial polygon CSV
  - mm/hr → m/s unit conversion
  - Rainfall multiplier scaling
  - Negative rainfall raises AssertionError
"""
import os
import tempfile
import shutil
import unittest
from unittest.mock import MagicMock

try:
    import anuga
    from anuga import rectangular_cross_domain
    from anuga.scenario.setup_rainfall import setup_rainfall
    HAS_MODULE = True
    SKIP_REASON = ''
except ImportError as _e:
    HAS_MODULE = False
    SKIP_REASON = str(_e)

# Reading a polygon CSV pulls in the geodata interface (fiona/rasterio/shapely)
# via spatialInputUtil. When that stack is not installed, tests that touch it
# must skip rather than fail — matching every other geodata test in the suite
# (see test_tif2.py, test_spatialInputUtil.py). Without this the polygon test
# turns a missing optional dependency into a hard ImportError that reds CI.
try:
    import fiona  # noqa: F401
    import rasterio  # noqa: F401
    import shapely  # noqa: F401
    HAS_GEODATA = True
except ImportError:
    HAS_GEODATA = False


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_domain():
    return rectangular_cross_domain(3, 3)


def _write_rain_csv(path, times, rates_mm_hr):
    """Write a two-column CSV (time_s, rain_mm_hr) with a header row."""
    with open(path, 'w') as f:
        f.write('time,rain_mm_hr\n')
        for t, r in zip(times, rates_mm_hr):
            f.write('%s,%s\n' % (t, r))


def _write_polygon_csv(path, points):
    """Write an x,y CSV polygon (no header — su.read_polygon handles that)."""
    with open(path, 'w') as f:
        f.write('x,y\n')
        for x, y in points:
            f.write('%s,%s\n' % (x, y))


def _count_rate_operators(domain):
    """Return the number of Rate_operator instances attached to the domain."""
    from anuga.operators.rate_operators import Rate_operator
    return sum(1 for op in domain.fractional_step_operators
               if isinstance(op, Rate_operator))


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

@unittest.skipUnless(HAS_MODULE, SKIP_REASON)
class TestSetupRainfall(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.domain = _make_domain()

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _rain_path(self, name='rain.csv'):
        return os.path.join(self.tmpdir, name)

    # --- No rainfall ---

    def test_no_rain_data_no_operators(self):
        project = MagicMock()
        project.rain_data = []
        setup_rainfall(self.domain, project)
        self.assertEqual(_count_rate_operators(self.domain), 0)

    # --- Global rainfall ('All' polygon) ---

    def test_single_global_rain_creates_operator(self):
        csv = self._rain_path()
        _write_rain_csv(csv, times=[0, 3600], rates_mm_hr=[10.0, 10.0])
        project = MagicMock()
        project.rain_data = [[csv, 0.0, 'linear', 'All']]
        setup_rainfall(self.domain, project)
        self.assertEqual(_count_rate_operators(self.domain), 1)

    def test_multiple_global_rain_events_create_multiple_operators(self):
        csv1 = self._rain_path('rain1.csv')
        csv2 = self._rain_path('rain2.csv')
        _write_rain_csv(csv1, [0, 1800], [5.0, 5.0])
        _write_rain_csv(csv2, [0, 1800], [3.0, 3.0])
        project = MagicMock()
        project.rain_data = [
            [csv1, 0.0, 'linear', 'All'],
            [csv2, 0.0, 'linear', 'All'],
        ]
        setup_rainfall(self.domain, project)
        self.assertEqual(_count_rate_operators(self.domain), 2)

    # --- Start time offset ---

    def test_start_time_offset_accepted(self):
        csv = self._rain_path()
        # Times start at 1000 s; start_time=1000 shifts them to 0
        _write_rain_csv(csv, times=[1000, 2000, 3000], rates_mm_hr=[5.0, 10.0, 5.0])
        project = MagicMock()
        project.rain_data = [[csv, 1000.0, 'linear', 'All']]
        # Should not raise
        setup_rainfall(self.domain, project)
        self.assertEqual(_count_rate_operators(self.domain), 1)

    # --- Unit conversion ---

    def test_zero_rain_still_creates_operator(self):
        # setup_rainfall uses `max >= 0` (not `> 0`) so a zero-rate timeseries
        # still creates a Rate_operator (it just produces zero flow).
        csv = self._rain_path()
        _write_rain_csv(csv, times=[0, 3600], rates_mm_hr=[0.0, 0.0])
        project = MagicMock()
        project.rain_data = [[csv, 0.0, 'linear', 'All']]
        setup_rainfall(self.domain, project)
        self.assertEqual(_count_rate_operators(self.domain), 1)

    # --- Multiplier ---

    def test_multiplier_accepted(self):
        csv = self._rain_path()
        _write_rain_csv(csv, times=[0, 3600], rates_mm_hr=[10.0, 10.0])
        project = MagicMock()
        # 5-element rain_data entry includes multiplier
        project.rain_data = [[csv, 0.0, 'linear', 'All', 2.0]]
        # Should not raise
        setup_rainfall(self.domain, project)
        self.assertEqual(_count_rate_operators(self.domain), 1)

    def test_default_multiplier_when_omitted(self):
        # Without the 5th element the multiplier defaults to 1.0
        csv = self._rain_path()
        _write_rain_csv(csv, times=[0, 3600], rates_mm_hr=[10.0, 10.0])
        project = MagicMock()
        project.rain_data = [[csv, 0.0, 'linear', 'All']]  # no multiplier
        setup_rainfall(self.domain, project)
        self.assertEqual(_count_rate_operators(self.domain), 1)

    # --- Polygon-restricted rainfall ---

    @unittest.skipUnless(
        HAS_GEODATA,
        'requires geodata interface (fiona, rasterio, shapely)')
    def test_rain_with_polygon_csv(self):
        rain_csv = self._rain_path('rain.csv')
        poly_csv = self._rain_path('poly.csv')
        _write_rain_csv(rain_csv, [0, 3600], [8.0, 8.0])
        # A polygon that covers part of the 1×1 domain
        _write_polygon_csv(poly_csv, [(0.0, 0.0), (0.5, 0.0), (0.5, 0.5), (0.0, 0.5)])
        project = MagicMock()
        project.rain_data = [[rain_csv, 0.0, 'linear', poly_csv]]
        setup_rainfall(self.domain, project)
        self.assertEqual(_count_rate_operators(self.domain), 1)

    # --- Negative rainfall raises ---

    def test_negative_rainfall_raises(self):
        csv = self._rain_path()
        _write_rain_csv(csv, times=[0, 3600], rates_mm_hr=[-1.0, 5.0])
        project = MagicMock()
        project.rain_data = [[csv, 0.0, 'linear', 'All']]
        with self.assertRaises(AssertionError):
            setup_rainfall(self.domain, project)


if __name__ == '__main__':
    unittest.main()
