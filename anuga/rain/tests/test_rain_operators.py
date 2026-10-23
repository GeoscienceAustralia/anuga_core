"""Tests for anuga/rain: Arr_hub_rain, Single_pattern, Raster_rate_operator,
ARR_rate_operator."""

import io
import os
import tempfile
import unittest
import zipfile

import numpy as np

from anuga import Domain, Reflective_boundary


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


if __name__ == '__main__':
    unittest.main()
