"""Tests for anuga/rain: Arr_hub_rain, Single_pattern, Raster_rate_operator,
ARR_rate_operator, Arr_grd, Arr_ifd_rain."""

import io
import os
import tempfile
import unittest
import warnings
import zipfile

import numpy as np
import pytest

from anuga import Domain, Reflective_boundary, rectangular_cross_domain
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.rain import ARR_rate_operator, Arr_grd, Arr_ifd_rain
from anuga.rain.raster_rate_operator import Raster_rate_operator


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_domain():
    """Minimal 4-triangle domain."""
    a = [0.0, 0.0]; b = [0.0, 2.0]; c = [2.0, 0.0]
    d = [0.0, 4.0]; e = [2.0, 2.0]; f = [4.0, 0.0]
    points   = [a, b, c, d, e, f]
    vertices = [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]]
    domain = Domain(points, vertices)
    domain.set_quantity('elevation', 0.0)
    domain.set_quantity('stage', 1.0)
    domain.set_quantity('friction', 0.0)
    Br = Reflective_boundary(domain)
    domain.set_boundary({'exterior': Br})
    return domain


def _make_hub_file(tmp_dir):
    """Write a minimal ARR hub text file and return its path."""
    content = (
        '[INPUTDATA]\n'
        'Latitude,-33.035717\n'
        'Longitude,151.265069\n'
        '[RIVREG]\n'
        'Division,South East Coast (NSW)\n'
        'RivRegNum,10\n'
        'RivRegName,Hunter River\n'
        '[RIVREG_META]\n'
        'TimeAccessed,17 September 2019 02:43PM\n'
        'Version,2016_v1\n'
        '[LONGARF]\n'
        'header\n'
        'a,0.1\n'
        'b,0.2\n'
        'c,0.3\n'
        'd,0.4\n'
        'e,0.5\n'
        'f,0.6\n'
        'g,0.7\n'
        'h,0.8\n'
        'i,0.9\n'
        '[LOSSES]\n'
        'header\n'
        'IL,5.0\n'
        'CL,2.0\n'
        '[TP]\n'
        'Tpat_code,EC\n'
        'Tpatlabel,East Coast\n'
        '[ATP]\n'
        'ATpat_code,ATEC\n'
        'ATpatlabel,AT East Coast\n'
    )
    path = os.path.join(tmp_dir, 'test_hub.txt')
    with open(path, 'w') as f:
        f.write(content)
    return path


def _make_pattern_zip(tmp_dir, tpat_code='EC', n_patterns=2,
                      ev_dur=60, tstep=10, ev_frq='1%AEP'):
    """
    Create a synthetic patterns zip with AllStats and Increments CSVs.

    Each pattern row has header + n_patterns data rows.
    With ev_dur=60, tstep=10 there are 6 time steps; percentages sum to 100.
    """
    tstps = ev_dur // tstep
    percentages = [100.0 / tstps] * tstps  # uniform distribution

    def make_increments_csv():
        lines = ['Ev_ID,Ev_dur,Tstep,Zone,Ev_Frq,' +
                 ','.join(f'T{i+1}' for i in range(tstps))]
        for i in range(1, n_patterns + 1):
            row = [f'ID{i}', str(ev_dur), str(tstep), 'EC', ev_frq]
            row += [f'{p:.4f}' for p in percentages]
            lines.append(','.join(row))
        return '\n'.join(lines).encode('utf-8')

    def make_allstats_csv():
        lines = ['Ev_ID,Stat1,Stat2']
        for i in range(1, n_patterns + 1):
            lines.append(f'ID{i},val1,val2')
        return '\n'.join(lines).encode('utf-8')

    zip_path = os.path.join(tmp_dir, 'patterns.zip')
    with zipfile.ZipFile(zip_path, 'w') as zf:
        zf.writestr(tpat_code + '_Increments.csv', make_increments_csv())
        zf.writestr(tpat_code + '_AllStats.csv',   make_allstats_csv())
    return zip_path


class _SimpleRasterData:
    """Minimal Raster_time_slice_data substitute for testing."""

    def __init__(self, n_slices=3, grid_size=5, time_step=600.0):
        nx = ny = grid_size
        self.x = np.linspace(0.0, 4.0, nx)
        self.y = np.linspace(0.0, 4.0, ny)
        self.time_step = time_step
        self.times = np.arange(n_slices) * time_step
        # Uniform depth (m per time_step) that increases with slice index
        self.data_slices = np.array(
            [(i + 1) * 0.001 * np.ones((ny, nx)) for i in range(n_slices)]
        )

    def extract_data_at_locations(self, locations):
        """Return uniform depth (same for every location) for each slice."""
        locations = np.asarray(locations)
        n_loc = len(locations)
        result = np.zeros((len(self.times), n_loc))
        for i, ds in enumerate(self.data_slices):
            result[i, :] = np.mean(ds)
        return result


# ---------------------------------------------------------------------------
# Tests: Arr_hub_rain
# ---------------------------------------------------------------------------

class TestArrHubRain(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.hub_path = _make_hub_file(self.tmp)

    def test_loads_location(self):
        from anuga.rain.arr_hub_rain import Arr_hub_rain
        hub = Arr_hub_rain(self.hub_path)
        self.assertEqual(hub.Loc_Lat, '-33.035717')
        self.assertEqual(hub.Loc_Lon, '151.265069')

    def test_loads_river_region(self):
        from anuga.rain.arr_hub_rain import Arr_hub_rain
        hub = Arr_hub_rain(self.hub_path)
        self.assertEqual(hub.RivName, 'Hunter River')
        self.assertEqual(hub.RivNum,  '10')

    def test_loads_version(self):
        from anuga.rain.arr_hub_rain import Arr_hub_rain
        hub = Arr_hub_rain(self.hub_path)
        self.assertEqual(hub.Version, '2016_v1')

    def test_loads_losses(self):
        from anuga.rain.arr_hub_rain import Arr_hub_rain
        hub = Arr_hub_rain(self.hub_path)
        self.assertEqual(hub.ARR_IL, '5.0')
        self.assertEqual(hub.ARR_CL, '2.0')

    def test_loads_tpat_code(self):
        from anuga.rain.arr_hub_rain import Arr_hub_rain
        hub = Arr_hub_rain(self.hub_path)
        self.assertEqual(hub.Tpat_code, 'EC')

    def test_missing_file_raises(self):
        from anuga.rain.arr_hub_rain import Arr_hub_rain
        with self.assertRaises(Exception):
            Arr_hub_rain(os.path.join(self.tmp, 'nonexistent.txt'))

    def test_accessible_from_package(self):
        import anuga.rain
        self.assertTrue(hasattr(anuga.rain, 'Arr_hub_rain'))


# ---------------------------------------------------------------------------
# Tests: ARR_point_rainfall_patterns and Single_pattern
# ---------------------------------------------------------------------------

class TestArrPatterns(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.zip_path = _make_pattern_zip(
            self.tmp, tpat_code='EC', n_patterns=5,
            ev_dur=60, tstep=10, ev_frq='1%AEP')

    def test_loads_patterns(self):
        from anuga.rain.arr_hub_rain import ARR_point_rainfall_patterns
        prp = ARR_point_rainfall_patterns(self.zip_path, 'EC')
        self.assertIsInstance(prp.linesInc, list)
        self.assertIsInstance(prp.linesAStat, list)

    def test_increments_labels(self):
        from anuga.rain.arr_hub_rain import ARR_point_rainfall_patterns
        prp = ARR_point_rainfall_patterns(self.zip_path, 'EC')
        self.assertIn('Ev_ID', prp.INCS_Labels[0])

    def test_missing_code_raises(self):
        from anuga.rain.arr_hub_rain import ARR_point_rainfall_patterns
        with self.assertRaises(IOError):
            ARR_point_rainfall_patterns(self.zip_path, 'BADCODE')

    def test_missing_zip_raises(self):
        from anuga.rain.arr_hub_rain import ARR_point_rainfall_patterns
        with self.assertRaises(IOError):
            ARR_point_rainfall_patterns('/no/such/file.zip', 'EC')

    def test_accessible_from_package(self):
        import anuga.rain
        self.assertTrue(hasattr(anuga.rain, 'ARR_point_rainfall_patterns'))


class TestSinglePattern(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.zip_path = _make_pattern_zip(
            self.tmp, tpat_code='EC', n_patterns=5,
            ev_dur=60, tstep=10, ev_frq='1%AEP')
        from anuga.rain.arr_hub_rain import ARR_point_rainfall_patterns
        self.prp = ARR_point_rainfall_patterns(self.zip_path, 'EC')

    def test_index_range(self):
        from anuga.rain.arr_hub_rain import Single_pattern
        with self.assertRaises(AssertionError):
            Single_pattern(self.prp, index=0)
        with self.assertRaises(AssertionError):
            Single_pattern(self.prp, index=721)

    def test_tplot_starts_zero(self):
        from anuga.rain.arr_hub_rain import Single_pattern
        sp = Single_pattern(self.prp, index=1, Ev_dep=100.0)
        self.assertEqual(sp.Tplot[0], 0.0)

    def test_tplot_length(self):
        from anuga.rain.arr_hub_rain import Single_pattern
        sp = Single_pattern(self.prp, index=1, Ev_dep=100.0)
        # 60-min event with 10-min steps → 6 steps + t=0
        self.assertEqual(len(sp.Tplot), 7)

    def test_rplot_sums_to_ev_dep(self):
        from anuga.rain.arr_hub_rain import Single_pattern
        Ev_dep = 80.0
        sp = Single_pattern(self.prp, index=1, Ev_dep=Ev_dep)
        self.assertAlmostEqual(float(np.sum(sp.Rplot)), Ev_dep, places=2)

    def test_rplot_starts_zero(self):
        from anuga.rain.arr_hub_rain import Single_pattern
        sp = Single_pattern(self.prp, index=1, Ev_dep=100.0)
        self.assertEqual(sp.Rplot[0], 0.0)

    def test_ev_dur_and_tstep(self):
        from anuga.rain.arr_hub_rain import Single_pattern
        sp = Single_pattern(self.prp, index=1, Ev_dep=100.0)
        self.assertEqual(sp.Ev_dur, 60)
        self.assertEqual(sp.Tstep, 10)

    def test_accessible_from_package(self):
        import anuga.rain
        self.assertTrue(hasattr(anuga.rain, 'Single_pattern'))


# ---------------------------------------------------------------------------
# Tests: Raster_rate_operator
# ---------------------------------------------------------------------------

class TestRasterRateOperator(unittest.TestCase):

    def setUp(self):
        self.domain = _make_domain()
        self.n = self.domain.number_of_elements  # 4

    def _make_raster(self, n_slices=3, time_step=600.0):
        return _SimpleRasterData(n_slices=n_slices, time_step=time_step)

    def test_init_no_error(self):
        from anuga.rain.raster_rate_operator import Raster_rate_operator
        raster = self._make_raster()
        op = Raster_rate_operator(self.domain, raster)
        self.assertIsNotNone(op)

    def test_rates_cache_shape(self):
        from anuga.rain.raster_rate_operator import Raster_rate_operator
        raster = self._make_raster(n_slices=4, time_step=300.0)
        op = Raster_rate_operator(self.domain, raster)
        self.assertEqual(op._rates_cache.shape, (4, self.n))

    def test_rates_in_ms(self):
        from anuga.rain.raster_rate_operator import Raster_rate_operator
        raster = self._make_raster(n_slices=3, time_step=600.0)
        op = Raster_rate_operator(self.domain, raster)
        # Slice 0: depth = 0.001 m/slice → rate = 0.001/600 m/s
        expected = 0.001 / 600.0
        np.testing.assert_allclose(op._rates_cache[0], expected, rtol=1e-6)

    def test_call_sets_rate(self):
        from anuga.rain.raster_rate_operator import Raster_rate_operator
        raster = self._make_raster(n_slices=3, time_step=600.0)
        op = Raster_rate_operator(self.domain, raster)
        # At t=0 we are in slice 0
        op()
        np.testing.assert_allclose(op.rate, op._rates_cache[0], rtol=1e-6)

    def test_call_advances_slice(self):
        from anuga.rain.raster_rate_operator import Raster_rate_operator
        raster = self._make_raster(n_slices=3, time_step=600.0)
        op = Raster_rate_operator(self.domain, raster)
        op()  # t=0, slice 0
        # Force domain time to 650 s → raster[1] (times=[0,600,1200])
        self.domain.set_time(650.0)
        op()
        np.testing.assert_allclose(op.rate, op._rates_cache[1], rtol=1e-6)

    def test_no_redundant_set_rate(self):
        """set_rate must only be called when the time slice changes."""
        from anuga.rain.raster_rate_operator import Raster_rate_operator
        raster = self._make_raster(n_slices=3, time_step=600.0)
        op = Raster_rate_operator(self.domain, raster)
        op()
        last_idx = op._last_slice_idx
        # Advance time by 1 second — still in same slice
        self.domain.set_time(1.0)
        op()
        self.assertEqual(op._last_slice_idx, last_idx)

    def test_time_offset(self):
        from anuga.rain.raster_rate_operator import Raster_rate_operator
        raster = self._make_raster(n_slices=3, time_step=600.0)
        # With offset=600 and t=0, effective raster time = 600 → slice index 1
        op = Raster_rate_operator(self.domain, raster, time_offset=600.0)
        op()
        self.assertEqual(op._last_slice_idx, 1)

    def test_accessible_from_package(self):
        import anuga.rain
        self.assertTrue(hasattr(anuga.rain, 'Raster_rate_operator'))

    def test_adds_water_to_stage(self):
        from anuga.rain.raster_rate_operator import Raster_rate_operator
        raster = self._make_raster(n_slices=3, time_step=600.0)
        op = Raster_rate_operator(self.domain, raster)
        stage_before = self.domain.quantities['stage'].centroid_values.copy()
        for _ in self.domain.evolve(yieldstep=1.0, finaltime=1.0):
            pass
        stage_after = self.domain.quantities['stage'].centroid_values
        self.assertTrue(np.all(stage_after >= stage_before))


# ---------------------------------------------------------------------------
# Tests: ARR_rate_operator
# ---------------------------------------------------------------------------

class TestARRRateOperator(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.zip_path = _make_pattern_zip(
            self.tmp, tpat_code='EC', n_patterns=5,
            ev_dur=60, tstep=10, ev_frq='1%AEP')
        from anuga.rain.arr_hub_rain import ARR_point_rainfall_patterns, Single_pattern
        prp = ARR_point_rainfall_patterns(self.zip_path, 'EC')
        self.pattern = Single_pattern(prp, index=1, Ev_dep=100.0)
        self.domain = _make_domain()
        self.n = self.domain.number_of_elements

    def test_init_no_error(self):
        from anuga.rain.raster_rate_operator import ARR_rate_operator
        depth_array = np.full(self.n, 50.0)
        op = ARR_rate_operator(self.domain, depth_array, self.pattern)
        self.assertIsNotNone(op)

    def test_rates_cache_shape(self):
        from anuga.rain.raster_rate_operator import ARR_rate_operator
        depth_array = np.full(self.n, 50.0)
        op = ARR_rate_operator(self.domain, depth_array, self.pattern)
        n_steps = len(self.pattern.Rplot)
        self.assertEqual(op._rates_cache.shape, (n_steps, self.n))

    def test_first_rate_is_zero(self):
        """Rplot[0] = 0 so the first rate slice must be all zeros."""
        from anuga.rain.raster_rate_operator import ARR_rate_operator
        depth_array = np.full(self.n, 50.0)
        op = ARR_rate_operator(self.domain, depth_array, self.pattern)
        np.testing.assert_array_equal(op._rates_cache[0], 0.0)

    def test_rates_positive_after_first(self):
        from anuga.rain.raster_rate_operator import ARR_rate_operator
        depth_array = np.full(self.n, 50.0)
        op = ARR_rate_operator(self.domain, depth_array, self.pattern)
        self.assertTrue(np.all(op._rates_cache[1:] >= 0.0))
        self.assertTrue(np.any(op._rates_cache[1:] > 0.0))

    def test_rate_sum_equals_total_depth(self):
        """Integrated rate × Tstep × 1000 should equal depth_array per centroid."""
        from anuga.rain.raster_rate_operator import ARR_rate_operator
        D = 80.0
        depth_array = np.full(self.n, D)
        op = ARR_rate_operator(self.domain, depth_array, self.pattern)
        Tstep_sec = self.pattern.Tstep * 60.0
        total_mm = np.sum(op._rates_cache, axis=0) * Tstep_sec * 1000.0
        np.testing.assert_allclose(total_mm, D, rtol=1e-5)

    def test_wrong_depth_shape_raises(self):
        from anuga.rain.raster_rate_operator import ARR_rate_operator
        with self.assertRaises(AssertionError):
            ARR_rate_operator(self.domain, np.zeros(self.n + 1), self.pattern)

    def test_call_sets_rate(self):
        """At t=0 the operator should apply Rplot[1] (the first non-zero step)."""
        from anuga.rain.raster_rate_operator import ARR_rate_operator
        depth_array = np.full(self.n, 50.0)
        op = ARR_rate_operator(self.domain, depth_array, self.pattern)
        op()
        np.testing.assert_allclose(op.rate, op._rates_cache[1], rtol=1e-6)

    def test_call_advances_step(self):
        """After Tstep+1 seconds the operator should switch to Rplot[2]."""
        from anuga.rain.raster_rate_operator import ARR_rate_operator
        depth_array = np.full(self.n, 50.0)
        op = ARR_rate_operator(self.domain, depth_array, self.pattern)
        op()  # t=0, idx=1
        Tstep_sec = self.pattern.Tstep * 60.0
        self.domain.set_time(Tstep_sec + 1.0)
        op()
        self.assertEqual(op._last_slice_idx, 2)
        np.testing.assert_allclose(op.rate, op._rates_cache[2], rtol=1e-6)

    def test_spatial_variation(self):
        """Depth variation → proportional rate variation."""
        from anuga.rain.raster_rate_operator import ARR_rate_operator
        depth_array = np.array([10.0, 20.0, 30.0, 40.0])
        op = ARR_rate_operator(self.domain, depth_array, self.pattern)
        # At step 1 (all percentages equal), rates should scale with depth
        Tstep_sec = self.pattern.Tstep * 60.0
        step1_rates = op._rates_cache[1]
        ratios = step1_rates / step1_rates[0]
        np.testing.assert_allclose(ratios, depth_array / depth_array[0], rtol=1e-6)

    def test_accessible_from_package(self):
        import anuga.rain
        self.assertTrue(hasattr(anuga.rain, 'ARR_rate_operator'))


# ---------------------------------------------------------------------------
# Tests: the time span of the rain data
# ---------------------------------------------------------------------------

class _Pattern:
    """10 mm over three 10-minute steps."""
    Tstep = 10
    Ev_dur = 30
    Ev_dep = 10.0
    Tplot = np.array([0.0, 10.0, 20.0, 30.0])
    Rplot = np.array([0.0, 2.0, 3.0, 5.0])


class _Raster:
    """2 mm then 3 mm, in two 10-minute slices."""
    time_step = 600.0
    times = np.array([0.0, 600.0])

    def extract_data_at_locations(self, locations):
        self.locations = np.array(locations)
        return np.tile([[0.002], [0.003]], (1, len(locations)))


def _flat_domain(mode, starttime=0.0):
    domain = rectangular_cross_domain(4, 4, len1=100.0, len2=100.0)
    domain.set_name('rain_span')
    domain.set_store(False)
    domain.set_quantity('elevation', 0.0)
    domain.set_quantity('friction', 0.03)
    domain.set_quantity('stage', 0.0)
    domain.set_boundary({tag: Reflective_boundary(domain)
                         for tag in domain.get_boundary_tags()})
    domain.set_starttime(starttime)
    domain.set_compute_mode(mode)
    return domain


def _delivered_mm(domain, finaltime, yieldstep=60.0):
    for _ in domain.evolve(yieldstep=yieldstep, finaltime=finaltime):
        pass
    return domain.quantities['stage'].get_integral() / 1.0e4 * 1000.0


@pytest.mark.parametrize('mode', ['legacy', 'unified'])
@pytest.mark.parametrize('starttime, finaltime', [(0.0, 1800.0), (0.0, 3600.0), (-1800.0, 3600.0)])
def test_arr_delivers_the_storm_depth_and_no_more(mode, starttime, finaltime):
    domain = _flat_domain(mode, starttime)
    ARR_rate_operator(domain, np.full(len(domain), 10.0), _Pattern)
    assert _delivered_mm(domain, finaltime) == pytest.approx(10.0, abs=1e-9)


def test_rain_operators_report_their_time_step():
    domain = _flat_domain('legacy')
    assert ARR_rate_operator(domain, np.full(len(domain), 10.0), _Pattern).time_step == 600.0
    assert Raster_rate_operator(domain, _Raster()).time_step == 600.0


@pytest.mark.parametrize('mode', ['legacy', 'unified'])
@pytest.mark.parametrize('divisor', [1, 3])
def test_yieldstep_from_time_step_delivers_the_exact_depth(mode, divisor):
    domain = _flat_domain(mode)
    rain = ARR_rate_operator(domain, np.full(len(domain), 10.0), _Pattern)
    delivered = _delivered_mm(domain, 1800.0, yieldstep=rain.time_step / divisor)
    assert delivered == pytest.approx(10.0, abs=1e-9)


@pytest.mark.parametrize('mode', ['legacy', 'unified'])
def test_yieldstep_off_the_time_steps_misses_the_depth(mode):
    """Why the docs say to use time_step: 1.5 x time_step skips a change of rate."""
    domain = _flat_domain(mode)
    rain = ARR_rate_operator(domain, np.full(len(domain), 10.0), _Pattern)
    delivered = _delivered_mm(domain, 1800.0, yieldstep=1.5 * rain.time_step)
    assert abs(delivered - 10.0) > 1e-3


@pytest.mark.parametrize('mode', ['legacy', 'unified'])
@pytest.mark.parametrize('starttime', [0.0, -1800.0])
def test_raster_delivers_its_slices_and_no_more(mode, starttime):
    domain = _flat_domain(mode, starttime)
    Raster_rate_operator(domain, _Raster())
    assert _delivered_mm(domain, 3600.0) == pytest.approx(5.0, abs=1e-9)


@pytest.mark.parametrize('make', ['arr', 'raster'])
def test_default_rate_applies_outside_the_data(make):
    domain = _flat_domain('legacy')
    if make == 'arr':
        op = ARR_rate_operator(domain, np.full(len(domain), 10.0), _Pattern,
                               default_rate=lambda t: 1.0e-6 * t)
        end = 1800.0
    else:
        op = Raster_rate_operator(domain, _Raster(), default_rate=lambda t: 1.0e-6 * t)
        end = 1200.0

    domain.set_time(end)
    op()
    np.testing.assert_allclose(op.rate, 1.0e-6 * end)

    # A default that varies with time is followed, and stays a centroid array.
    domain.set_time(end + 100.0)
    op()
    np.testing.assert_allclose(op.rate, 1.0e-6 * (end + 100.0))
    assert op.rate_type == 'centroid_array'

    domain.set_time(0.0)
    op()
    np.testing.assert_allclose(op.rate, op._rates_cache[1 if make == 'arr' else 0])


def test_raster_is_read_at_absolute_coordinates():
    domain = _flat_domain('legacy')
    domain.set_georeference(Geo_reference(zone=56, xllcorner=300000.0, yllcorner=6000000.0))
    raster = _Raster()
    Raster_rate_operator(domain, raster)
    np.testing.assert_allclose(raster.locations,
                               domain.get_centroid_coordinates(absolute=True))


# ---------------------------------------------------------------------------
# Tests: ARR_rate_operator inputs
# ---------------------------------------------------------------------------

def test_ev_dep_has_no_effect_on_the_rates():
    class Scaled(_Pattern):
        Ev_dep = 100.0
        Rplot = _Pattern.Rplot * 10.0

    domain = _flat_domain('legacy')
    depth = np.full(len(domain), 10.0)
    np.testing.assert_allclose(ARR_rate_operator(domain, depth, _Pattern)._rates_cache,
                               ARR_rate_operator(domain, depth, Scaled)._rates_cache)


@pytest.mark.parametrize('bad', [np.nan, -9999.0, np.inf])
def test_arr_rejects_bad_depths_only_where_it_acts(bad):
    domain = _flat_domain('legacy')
    depth = np.full(len(domain), 10.0)
    depth[0] = bad
    with pytest.raises(ValueError, match='depth_array'):
        ARR_rate_operator(domain, depth, _Pattern)
    ARR_rate_operator(domain, depth, _Pattern, indices=[1, 2, 3])


def test_ifd_duration_warns_on_a_mismatch():
    domain = _flat_domain('legacy')
    depth = np.full(len(domain), 10.0)
    with pytest.warns(UserWarning, match='60 min.*30 min'):
        ARR_rate_operator(domain, depth, _Pattern, ifd_duration=60)
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        ARR_rate_operator(domain, depth, _Pattern, ifd_duration=30)
        ARR_rate_operator(domain, depth, _Pattern, ifd_duration='30')
        ARR_rate_operator(domain, depth, _Pattern)


# ---------------------------------------------------------------------------
# Tests: Arr_grd and Arr_ifd_rain
# ---------------------------------------------------------------------------

XLL, YLL, CELL, NCOLS, NROWS = 150.0, -34.0, 0.025, 4, 3


def _grid_values():
    """Row 0 is the northernmost; each value identifies its cell."""
    return [[100.0 * (row + 1) + col for col in range(NCOLS)] for row in range(NROWS)]


@pytest.fixture
def make_grd(tmp_path):
    def make(values=None, name='catchment_depth_60min_01000aep.txt.asc'):
        values = _grid_values() if values is None else values
        text = (f'ncols {NCOLS}\nnrows {NROWS}\nxllcorner {XLL}\nyllcorner {YLL}\n'
                f'cellsize {CELL}\nNODATA_value -9999\n')
        text += ''.join(' '.join(str(v) for v in row) + '\n' for row in values)
        with zipfile.ZipFile(tmp_path / 'IFD.zip', 'w') as archive:
            archive.writestr(name, text)
        return Arr_grd(60, 1000, name, str(tmp_path / 'IFD.zip'), str(tmp_path))
    return make


@pytest.mark.parametrize('fx', [0.0, 0.1, 0.49, 0.51, 0.9])
@pytest.mark.parametrize('fy', [0.0, 0.1, 0.49, 0.51, 0.9])
def test_grid_lookup_returns_the_cell_holding_the_point(make_grd, fx, fy):
    grd = make_grd()
    values = _grid_values()
    for from_south in range(NROWS):
        for col in range(NCOLS):
            lon = XLL + (col + fx) * CELL
            lat = YLL + (from_south + fy) * CELL
            expected = values[NROWS - 1 - from_south][col]
            assert grd.get_rain_at_point(lon, lat) == expected
            assert grd.get_rain_at_points([lon], [lat])[0] == expected


def test_grid_lookup_on_the_outer_edges(make_grd):
    grd = make_grd()
    east, north = XLL + NCOLS * CELL, YLL + NROWS * CELL
    # (Lon, Lat, row from the north, column)
    for lon, lat, row, col in ((XLL, YLL, NROWS - 1, 0),
                               (east, YLL, NROWS - 1, NCOLS - 1),
                               (XLL, north, 0, 0),
                               (east, north, 0, NCOLS - 1)):
        assert grd.get_rain_at_point(lon, lat) == _grid_values()[row][col]


@pytest.mark.parametrize('lon, lat', [(XLL - 0.01, YLL + 0.01), (XLL + NCOLS * CELL + 0.01, YLL + 0.01),
                                      (XLL + 0.01, YLL - 0.01), (XLL + 0.01, YLL + NROWS * CELL + 0.01)])
def test_grid_lookup_rejects_points_outside(make_grd, lon, lat):
    grd = make_grd()
    with pytest.raises(ValueError, match='outside the IFD grid'):
        grd.get_rain_at_point(lon, lat)
    with pytest.raises(ValueError, match='1 of 2 points'):
        grd.get_rain_at_points([XLL + 0.01, lon], [YLL + 0.01, lat])


def test_nodata_is_nan_with_a_warning(make_grd):
    values = _grid_values()
    values[0][1] = -9999
    grd = make_grd(values)
    assert grd.NODATA_value == -9999.0
    assert np.isnan(grd.IFD_Data[0, 1])
    assert np.nanmin(grd.IFD_Data) > 0.0

    north = YLL + (NROWS - 0.5) * CELL
    with pytest.warns(UserWarning, match='1 of 2 points fall on NODATA'):
        out = grd.get_rain_at_points([XLL + 1.5 * CELL, XLL + 0.5 * CELL], [north, north])
    assert np.isnan(out[0]) and out[1] == values[0][0]


def test_plot_at_cell_centres_with_nodata(make_grd):
    plt = pytest.importorskip('matplotlib.pyplot')
    values = _grid_values()
    values[0][1] = -9999
    grd = make_grd(values)
    ax = grd.plot()
    # Plotted at the cell centres, the mesh's cell edges are the grid's edges.
    assert ax.collections[0].get_datalim(ax.transData).x0 == pytest.approx(XLL)
    plt.close('all')


def test_arr_ifd_rain_opens_a_grid_by_duration_and_aep(tmp_path, make_grd):
    make_grd()  # writes tmp_path/IFD.zip
    ifd = Arr_ifd_rain(str(tmp_path), Lat=-33.97, Lon=150.01)
    grd = ifd.open_grd(Dur=60, Frq=1000)
    assert (int(grd.dur), int(grd.frq)) == (60, 1000)
    assert grd.get_rain_at_point(XLL + 0.5 * CELL, YLL + 0.5 * CELL) == _grid_values()[-1][0]


def _ifd_dir(tmp_path, zips):
    """Write zip files of grids, each grid's values set to its own AEP label."""
    for zip_name, grids in zips.items():
        with zipfile.ZipFile(tmp_path / zip_name, 'w') as archive:
            archive.writestr('IFD2016_projection _epsg4283.prj', 'Projection GEOGRAPHIC')
            for dur, label in grids:
                text = (f'ncols {NCOLS}\r\nnrows {NROWS}\r\nxllcorner {XLL}\r\n'
                        f'yllcorner {YLL}\r\ncellsize {CELL}\r\nNODATA_value  -9999\r\n')
                text += (' '.join([str(int(label))] * NCOLS) + ' \r\n') * NROWS
                archive.writestr(f'catchment_depth_{dur}min_{label}aep.txt.asc', text)
    return str(tmp_path)


AEP_GRIDS = [(60, label) for label in ('00050', '00100', '01000', '39347', '63212')]


@pytest.mark.parametrize('aep, label', [(0.05, 50), (0.1, 100), (1.0, 1000), (1, 1000),
                                        (39.347, 39347), (39.35, 39347), (63.2, 63212)])
def test_open_grd_by_aep_percent(tmp_path, aep, label):
    ifd = Arr_ifd_rain(_ifd_dir(tmp_path, {'Rare_Extent.zip': AEP_GRIDS}), -33.97, 150.01)
    grd = ifd.open_grd(Dur=60, aep_percent=aep)
    assert int(grd.frq) == label
    assert grd.IFD_Data[0, 0] == label


def test_open_grd_frq_is_the_file_label_not_a_return_period(tmp_path):
    """Frq=100 is the 1 in 1000 grid; 1% AEP is 1000."""
    ifd = Arr_ifd_rain(_ifd_dir(tmp_path, {'Rare_Extent.zip': AEP_GRIDS}), -33.97, 150.01)
    assert ifd.open_grd(Dur=60, Frq=100).IFD_Data[0, 0] == ifd.open_grd(Dur=60, aep_percent=0.1).IFD_Data[0, 0]
    assert ifd.open_grd(Dur=60, Frq=1000).IFD_Data[0, 0] == ifd.open_grd(Dur=60, aep_percent=1.0).IFD_Data[0, 0]


@pytest.mark.parametrize('kwargs', [dict(Dur=60), dict(Dur=60, Frq=1000, aep_percent=1.0)])
def test_open_grd_needs_exactly_one_probability(tmp_path, kwargs):
    ifd = Arr_ifd_rain(_ifd_dir(tmp_path, {'x.zip': AEP_GRIDS}), -33.97, 150.01)
    with pytest.raises(ValueError, match='exactly one'):
        ifd.open_grd(**kwargs)


@pytest.mark.parametrize('kwargs', [dict(Dur=30, aep_percent=1.0), dict(Dur=60, aep_percent=2.0),
                                    dict(Dur=60, Frq=2000), dict(Dur=60, aep_percent=0.11)])
def test_open_grd_raises_and_lists_what_is_available(tmp_path, kwargs):
    ifd = Arr_ifd_rain(_ifd_dir(tmp_path, {'x.zip': AEP_GRIDS}), -33.97, 150.01)
    with pytest.raises(ValueError, match=r'No IFD grid.*Durations \(min\): \[60\].*0\.05'):
        ifd.open_grd(**kwargs)


def test_arr_ifd_rain_reads_zips_under_any_name(tmp_path):
    (tmp_path / 'notes.zip').write_text('not a zip')
    with zipfile.ZipFile(tmp_path / 'other.zip', 'w') as archive:
        archive.writestr('readme.txt', 'nothing to see')
    _ifd_dir(tmp_path, {'IFD_Extent.zip': [(60, '01000')],
                        'Rare_Extent.zip': [(60, '01000'), (60, '00050')],
                        'Very_frequent_Extent.zip': [(60, '99999')]})
    with pytest.warns(UserWarning, match='notes.zip'):
        ifd = Arr_ifd_rain(str(tmp_path), -33.97, 150.01)
    assert len(ifd.file_list) == 4  # .prj and readme.txt skipped
    assert ifd.open_grd(Dur=60, aep_percent=1.0).IFD_Data[0, 0] == 1000
    assert ifd.open_grd(Dur=60, aep_percent=99.999).IFD_Data[0, 0] == 99999


def test_arr_ifd_rain_warns_when_it_finds_no_grids(tmp_path):
    with pytest.warns(UserWarning, match='No IFD grids'):
        Arr_ifd_rain(str(tmp_path), -33.97, 150.01)


def test_a_depth_of_zero_warns(make_grd):
    """The Bureau's grids hold 0, not NODATA, over the ocean."""
    values = _grid_values()
    values[0][1] = 0
    grd = make_grd(values)
    north = YLL + (NROWS - 0.5) * CELL
    with pytest.warns(UserWarning, match='1 of 2 points read a depth of 0'):
        out = grd.get_rain_at_points([XLL + 1.5 * CELL, XLL + 0.5 * CELL], [north, north])
    assert out[0] == 0.0 and out[1] == values[0][0]


def test_documented_centroid_lookup(make_grd):
    """The recipe in docs/source/setup_anuga_script/rainfall.rst."""
    pytest.importorskip('pyproj')
    from anuga.coordinate_transforms.redfearn import epsg_to_ll, ll_to_epsg

    easting, northing = ll_to_epsg(YLL + 0.03, XLL + 0.05, 32756)
    domain = _flat_domain('legacy')
    domain.set_georeference(Geo_reference(epsg=32756, xllcorner=float(easting),
                                          yllcorner=float(northing)))
    xy = domain.get_centroid_coordinates(absolute=True)
    lats, lons = epsg_to_ll(xy[:, 0], xy[:, 1], domain.get_epsg())
    assert make_grd().get_rain_at_points(lons, lats).shape == (len(domain),)


if __name__ == '__main__':
    unittest.main()
