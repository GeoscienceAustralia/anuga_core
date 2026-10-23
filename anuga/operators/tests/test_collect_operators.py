"""
Tests for Collect_max_stage_operator and Collect_max_quantities_operator.
"""

import unittest
import numpy as num
import anuga
from anuga import Reflective_boundary


def make_domain():
    """4-triangle domain, 1m water over 0.5m elevation."""
    a = [0.0, 0.0]; b = [0.0, 2.0]; c = [2.0, 0.0]
    d = [0.0, 4.0]; e = [2.0, 2.0]; f = [4.0, 0.0]
    points = [a, b, c, d, e, f]
    vertices = [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]]
    domain = anuga.Domain(points, vertices)
    domain.set_quantity('elevation', 0.5)
    domain.set_quantity('stage', 1.0)
    domain.set_quantity('friction', 0.0)
    domain.set_boundary({'exterior': Reflective_boundary(domain)})
    return domain


class Test_collect_max_stage_operator(unittest.TestCase):

    def setUp(self):
        self.domain = make_domain()

    def tearDown(self):
        import os
        try:
            os.remove('domain.sww')
        except OSError:
            pass

    def test_construction(self):
        from anuga.operators.collect_max_stage_operator import Collect_max_stage_operator
        operator = Collect_max_stage_operator(self.domain)
        self.assertIsNotNone(operator)

    def test_max_stage_quantity_initialised(self):
        """max_stage quantity should exist and be initialised to -1e100."""
        from anuga.operators.collect_max_stage_operator import Collect_max_stage_operator
        operator = Collect_max_stage_operator(self.domain)
        self.assertIn('max_stage', self.domain.quantities)
        max_stage_vals = self.domain.quantities['max_stage'].centroid_values
        self.assertTrue(num.all(max_stage_vals <= -1.0e+99),
                        "max_stage should be initialised to -1e100")

    def test_call_updates_max_stage(self):
        """After calling the operator, max_stage should be >= stage."""
        from anuga.operators.collect_max_stage_operator import Collect_max_stage_operator
        operator = Collect_max_stage_operator(self.domain)
        self.domain.timestep = 1.0
        operator()
        max_stage_vals = self.domain.quantities['max_stage'].centroid_values
        stage_vals = self.domain.quantities['stage'].centroid_values
        self.assertTrue(num.all(max_stage_vals >= stage_vals),
                        "max_stage should be >= stage after operator call")

    def test_parallel_safe(self):
        from anuga.operators.collect_max_stage_operator import Collect_max_stage_operator
        operator = Collect_max_stage_operator(self.domain)
        self.assertTrue(operator.parallel_safe())


class Test_collect_max_quantities_operator(unittest.TestCase):

    def setUp(self):
        self.domain = make_domain()

    def tearDown(self):
        import os
        try:
            os.remove('domain.sww')
        except OSError:
            pass

    def test_construction_defaults(self):
        from anuga.operators.collect_max_quantities_operator import Collect_max_quantities_operator
        operator = Collect_max_quantities_operator(self.domain)
        self.assertIsNotNone(operator)

    def test_update_frequency_and_start_time_stored(self):
        from anuga.operators.collect_max_quantities_operator import Collect_max_quantities_operator
        operator = Collect_max_quantities_operator(
            self.domain, update_frequency=3, collection_start_time=5.0)
        self.assertEqual(operator.update_frequency, 3)
        self.assertEqual(operator.collection_start_time, 5.0)

    def test_call_updates_max_stage_when_time_exceeds_start(self):
        """After setting domain time > collection_start_time and calling operator,
        max_stage should be updated (>= -1e30)."""
        from anuga.operators.collect_max_quantities_operator import Collect_max_quantities_operator
        operator = Collect_max_quantities_operator(
            self.domain, update_frequency=1, collection_start_time=0.0)
        # domain.get_time() uses starttime + relative_time
        self.domain.relative_time = 1.0
        operator()
        # max_stage should now have been updated from -max_float to actual stage values
        self.assertTrue(num.all(operator.max_stage >= -1.0e30),
                        "max_stage should be updated after operator call")

    def test_call_no_update_before_start_time(self):
        """Before collection_start_time, max_stage should remain at initialised value."""
        from anuga.operators.collect_max_quantities_operator import Collect_max_quantities_operator
        from anuga.config import max_float
        operator = Collect_max_quantities_operator(
            self.domain, update_frequency=1, collection_start_time=10.0)
        # domain.time defaults to 0.0, which is <= collection_start_time
        operator()
        # max_stage should remain at the initialised -max_float
        self.assertTrue(num.all(operator.max_stage < -1.0e30),
                        "max_stage should not be updated before collection_start_time")

    def test_parallel_safe(self):
        from anuga.operators.collect_max_quantities_operator import Collect_max_quantities_operator
        operator = Collect_max_quantities_operator(self.domain)
        self.assertTrue(operator.parallel_safe())


class Test_collect_max_quantities_operator_extra(unittest.TestCase):
    """Tests for uncovered Collect_max_quantities_operator methods."""

    def setUp(self):
        import anuga
        self.domain = anuga.rectangular_cross_domain(2, 2)
        self.domain.set_quantity('elevation', 0.0)
        self.domain.set_quantity('stage', 1.0)

    def test_velocity_zero_height_not_none(self):
        """velocity_zero_height kwarg stored when not None (line 80)."""
        from anuga.operators.collect_max_quantities_operator import Collect_max_quantities_operator
        op = Collect_max_quantities_operator(self.domain, velocity_zero_height=0.01)
        self.assertAlmostEqual(op.velocity_zero_height, 0.01)

    def test_statistics_quantities(self):
        """statistics returns a string (lines 119-120)."""
        from anuga.operators.collect_max_quantities_operator import Collect_max_quantities_operator
        op = Collect_max_quantities_operator(self.domain)
        msg = op.statistics()
        self.assertIsInstance(msg, str)

    def test_timestepping_statistics_quantities(self):
        """timestepping_statistics returns a string (lines 124, 126-127)."""
        from anuga.operators.collect_max_quantities_operator import Collect_max_quantities_operator
        op = Collect_max_quantities_operator(self.domain)
        msg = op.timestepping_statistics()
        self.assertIsInstance(msg, str)


class Test_collect_max_stage_operator_extra(unittest.TestCase):
    """Tests for uncovered Collect_max_stage_operator methods."""

    def setUp(self):
        import anuga
        self.domain = anuga.rectangular_cross_domain(2, 2)
        self.domain.set_quantity('elevation', 0.0)
        self.domain.set_quantity('stage', 1.0)

    def test_statistics(self):
        """statistics returns a string (lines 63-64)."""
        from anuga.operators.collect_max_stage_operator import Collect_max_stage_operator
        op = Collect_max_stage_operator(self.domain)
        msg = op.statistics()
        self.assertIsInstance(msg, str)

    def test_timestepping_statistics(self):
        """timestepping_statistics returns a string (lines 68, 70-71)."""
        from anuga.operators.collect_max_stage_operator import Collect_max_stage_operator
        op = Collect_max_stage_operator(self.domain)
        msg = op.timestepping_statistics()
        self.assertIsInstance(msg, str)

    def test_save_centroid_data_to_csv(self):
        """save_centroid_data_to_csv delegates to max_stage (line 76)."""
        import tempfile
        import os
        from anuga.operators.collect_max_stage_operator import Collect_max_stage_operator
        op = Collect_max_stage_operator(self.domain)
        op()  # run once to populate max_stage
        orig = os.getcwd()
        with tempfile.TemporaryDirectory() as tmpdir:
            os.chdir(tmpdir)
            try:
                op.save_centroid_data_to_csv()
            finally:
                os.chdir(orig)


class Test_sww_merge_parallel_max_quantities(unittest.TestCase):
    """Verify that sww_merge_parallel correctly assembles the four max centroid
    quantities (max_stage_c, max_depth_c, max_speed_c, max_uh_c) written by
    Collect_max_quantities_operator / set_collect_max_quantities from
    per-partition parallel SWW files into the merged global file.

    Files are created in non-smooth parallel format (3*n_vols == n_points),
    which is what ANUGA parallel runs produce and what
    _sww_merge_parallel_non_smooth handles.
    """

    def setUp(self):
        import tempfile
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    # ------------------------------------------------------------------

    def _write_partition_sww(self, path, n_global_tris,
                             tri_l2g, tri_full_flag,
                             times,
                             max_stage_c, max_depth_c,
                             max_speed_c, max_uh_c):
        """Write a minimal non-smooth parallel SWW file with max centroid
        quantities.  Stage/elevation/momentum data are all zeros; only the
        max quantity values are checked in the assertions."""
        import numpy as np
        from anuga.file.sww import Write_sww
        from anuga.file.netcdf import NetCDFFile
        from anuga.config import netcdf_mode_w, netcdf_int
        from anuga.coordinate_transforms.geo_reference import Geo_reference

        n_local = len(tri_l2g)
        n_pts   = 3 * n_local          # non-smooth: 3 vertex copies per tri

        fid = NetCDFFile(path, netcdf_mode_w)
        sww = Write_sww(
            static_quantities=['elevation'],
            dynamic_quantities=['stage', 'xmomentum', 'ymomentum'],
            static_c_quantities=['max_stage_c', 'max_depth_c',
                                 'max_speed_c', 'max_uh_c'],
            dynamic_c_quantities=[]
        )
        t = num.array(times, dtype=float)
        sww.store_header(fid, t, n_local, n_pts, smoothing=False, order=1)

        # Trivial geometry: unit triangles stacked at the same position
        x   = num.tile([0., 1., 0.], n_local).astype(num.float32)
        y   = num.tile([0., 0., 1.], n_local).astype(num.float32)
        pts = num.column_stack([x, y])
        vols = num.arange(n_pts).reshape(-1, 3)
        geo  = Geo_reference()
        sww.store_triangulation(fid, pts, vols, points_georeference=geo)

        # Parallel metadata
        fid.number_of_global_triangles = int(n_global_tris)
        fid.number_of_global_nodes = int(3 * n_global_tris)
        fid.createVariable('tri_l2g',      netcdf_int, ('number_of_volumes',))
        fid.createVariable('tri_full_flag', netcdf_int, ('number_of_volumes',))
        fid.variables['tri_l2g'][:]       = num.array(tri_l2g,      dtype=num.int32)
        fid.variables['tri_full_flag'][:] = num.array(tri_full_flag, dtype=num.int32)

        sww.store_static_quantities(fid, elevation=num.zeros(n_pts, num.float32))
        sww.store_static_quantities_centroid(
            fid,
            max_stage_c=num.array(max_stage_c, dtype=num.float32),
            max_depth_c=num.array(max_depth_c, dtype=num.float32),
            max_speed_c=num.array(max_speed_c, dtype=num.float32),
            max_uh_c   =num.array(max_uh_c,    dtype=num.float32),
        )

        n_steps = len(times)
        zeros   = num.zeros((n_steps, n_pts), dtype=num.float32)
        for q in ('stage', 'xmomentum', 'ymomentum'):
            for i in range(n_steps):
                fid.variables[q][i] = zeros[i]

        fid.close()

    # ------------------------------------------------------------------

    def test_max_quantities_survive_merge(self):
        """max_stage/depth/speed/uh are correctly assembled into the global
        file when the two partitions each own half the triangles."""
        import os
        from anuga.utilities.sww_merge import sww_merge_parallel
        from anuga.file.netcdf import NetCDFFile
        from anuga.config import netcdf_mode_r

        base = os.path.join(self.tmpdir, 'run')
        times = [0.0, 1.0, 2.0]

        # Partition 0: global triangles 0, 1
        ms0  = [1.1, 2.2];  md0 = [0.5, 0.8]
        msp0 = [0.3, 0.6];  mu0 = [0.4, 0.7]
        self._write_partition_sww(
            base + '_P2_0.sww', 4,
            tri_l2g=[0, 1], tri_full_flag=[1, 1],
            times=times,
            max_stage_c=ms0, max_depth_c=md0,
            max_speed_c=msp0, max_uh_c=mu0)

        # Partition 1: global triangles 2, 3
        ms1  = [3.3, 4.4];  md1 = [1.2, 1.5]
        msp1 = [0.9, 1.1];  mu1 = [1.0, 1.2]
        self._write_partition_sww(
            base + '_P2_1.sww', 4,
            tri_l2g=[2, 3], tri_full_flag=[1, 1],
            times=times,
            max_stage_c=ms1, max_depth_c=md1,
            max_speed_c=msp1, max_uh_c=mu1)

        sww_merge_parallel(base, 2)
        self.assertTrue(os.path.exists(base + '.sww'))

        fid = NetCDFFile(base + '.sww', netcdf_mode_r)
        got_ms  = num.array(fid.variables['max_stage_c'][:], dtype=num.float32)
        got_md  = num.array(fid.variables['max_depth_c'][:], dtype=num.float32)
        got_msp = num.array(fid.variables['max_speed_c'][:], dtype=num.float32)
        got_mu  = num.array(fid.variables['max_uh_c'][:],    dtype=num.float32)
        fid.close()

        num.testing.assert_allclose(got_ms,  num.array(ms0  + ms1,  dtype=num.float32), rtol=1e-5)
        num.testing.assert_allclose(got_md,  num.array(md0  + md1,  dtype=num.float32), rtol=1e-5)
        num.testing.assert_allclose(got_msp, num.array(msp0 + msp1, dtype=num.float32), rtol=1e-5)
        num.testing.assert_allclose(got_mu,  num.array(mu0  + mu1,  dtype=num.float32), rtol=1e-5)

    def test_merge_ignores_ghost_triangle_values(self):
        """Ghost triangles (tri_full_flag=0) must not overwrite the correct
        max values contributed by the owning partition."""
        import os
        from anuga.utilities.sww_merge import sww_merge_parallel
        from anuga.file.netcdf import NetCDFFile
        from anuga.config import netcdf_mode_r

        SENTINEL = -999.0   # placed in ghost slots; must never appear in output

        base  = os.path.join(self.tmpdir, 'run_ghost')
        times = [0.0, 1.0, 2.0]

        # Partition 0: full tris 0,1 + ghost copies of 2,3
        self._write_partition_sww(
            base + '_P2_0.sww', 4,
            tri_l2g=[0, 1, 2, 3], tri_full_flag=[1, 1, 0, 0],
            times=times,
            max_stage_c=[1.1, 2.2, SENTINEL, SENTINEL],
            max_depth_c=[0.5, 0.8, SENTINEL, SENTINEL],
            max_speed_c=[0.3, 0.6, SENTINEL, SENTINEL],
            max_uh_c   =[0.4, 0.7, SENTINEL, SENTINEL])

        # Partition 1: full tris 2,3 + ghost copies of 0,1
        self._write_partition_sww(
            base + '_P2_1.sww', 4,
            tri_l2g=[0, 1, 2, 3], tri_full_flag=[0, 0, 1, 1],
            times=times,
            max_stage_c=[SENTINEL, SENTINEL, 3.3, 4.4],
            max_depth_c=[SENTINEL, SENTINEL, 1.2, 1.5],
            max_speed_c=[SENTINEL, SENTINEL, 0.9, 1.1],
            max_uh_c   =[SENTINEL, SENTINEL, 1.0, 1.2])

        sww_merge_parallel(base, 2)

        fid = NetCDFFile(base + '.sww', netcdf_mode_r)
        got_ms = num.array(fid.variables['max_stage_c'][:], dtype=num.float32)
        got_mu = num.array(fid.variables['max_uh_c'][:],    dtype=num.float32)
        fid.close()

        num.testing.assert_allclose(
            got_ms, num.array([1.1, 2.2, 3.3, 4.4], dtype=num.float32), rtol=1e-5)
        num.testing.assert_allclose(
            got_mu, num.array([0.4, 0.7, 1.0, 1.2], dtype=num.float32), rtol=1e-5)
        self.assertFalse(num.any(got_ms < -100.),
                         "Sentinel ghost values must not appear in merged output")


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromModule(__import__('__main__'))
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
