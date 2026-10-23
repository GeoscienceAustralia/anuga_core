"""Tests for shallow_water_domain flow algorithm setters."""
import unittest
import anuga
import numpy as num


def _make_domain():
    domain = anuga.rectangular_cross_domain(3, 3, len1=1.0, len2=1.0)
    domain.set_quantity('elevation', 0.0)
    domain.set_quantity('stage', 0.5)
    return domain


class Test_set_flow_algorithm(unittest.TestCase):

    def test_DE0_sets_algorithm(self):
        domain = _make_domain()
        domain.set_flow_algorithm('DE0')
        self.assertEqual(domain.get_flow_algorithm(), 'DE0')

    def test_DE1_sets_algorithm(self):
        domain = _make_domain()
        domain.set_flow_algorithm('DE1')
        self.assertEqual(domain.get_flow_algorithm(), 'DE1')

    def test_DE2_sets_algorithm(self):
        domain = _make_domain()
        domain.set_flow_algorithm('DE2')
        self.assertEqual(domain.get_flow_algorithm(), 'DE2')

    def test_DE0_7_sets_algorithm(self):
        domain = _make_domain()
        domain.set_flow_algorithm('DE0_7')
        self.assertEqual(domain.get_flow_algorithm(), 'DE0_7')

    def test_DE1_7_sets_algorithm(self):
        domain = _make_domain()
        domain.set_flow_algorithm('DE1_7')
        self.assertEqual(domain.get_flow_algorithm(), 'DE1_7')

    def test_unknown_algorithm_raises(self):
        domain = _make_domain()
        with self.assertRaises(Exception):
            domain.set_flow_algorithm('INVALID')

    def test_dot_replaced_by_underscore(self):
        """Dots in algorithm name are normalised to underscores."""
        domain = _make_domain()
        domain.set_flow_algorithm('DE0.7')
        self.assertEqual(domain.get_flow_algorithm(), 'DE0_7')

    def test_DE1_dot_7(self):
        domain = _make_domain()
        domain.set_flow_algorithm('DE1.7')
        self.assertEqual(domain.get_flow_algorithm(), 'DE1_7')


class Test_timestep_output_limits(unittest.TestCase):

    def test_update_timestep_clamps_stale_yield_deadline_to_zero(self):
        domain = _make_domain()
        domain.relative_time = 1000.0
        domain.relative_yieldtime = num.nextafter(domain.relative_time, -num.inf)
        domain.relative_finaltime = None
        domain.flux_timestep = 10.0

        domain.update_timestep(yieldstep=1.0, finaltime=None)

        self.assertEqual(domain.timestep, 0.0)
        self.assertEqual(domain._get_max_timestep_to_output_times(1.0, None), 0.0)

    def test_update_timestep_clamps_stale_final_deadline_to_zero(self):
        domain = _make_domain()
        domain.relative_time = 1000.0
        domain.relative_yieldtime = domain.relative_time + 10.0
        domain.relative_finaltime = num.nextafter(domain.relative_time, -num.inf)
        domain.flux_timestep = 10.0

        domain.update_timestep(yieldstep=10.0, finaltime=domain.relative_finaltime)

        self.assertEqual(domain.timestep, 0.0)
        self.assertEqual(domain._get_max_timestep_to_output_times(10.0, domain.relative_finaltime), 0.0)


class Test_DE0_defaults(unittest.TestCase):

    def setUp(self):
        self.domain = _make_domain()
        self.domain.set_flow_algorithm('DE0')

    def test_cfl(self):
        self.assertAlmostEqual(self.domain.CFL, 0.9)

    def test_timestepping_euler(self):
        self.assertEqual(self.domain.timestepping_method, 'euler')

    def test_discontinuous_elevation(self):
        self.assertTrue(self.domain.get_using_discontinuous_elevation())

    def test_compute_fluxes_method(self):
        self.assertEqual(self.domain.get_compute_fluxes_method(), 'DE')

    def test_beta_w(self):
        self.assertAlmostEqual(self.domain.beta_w, 0.5)

    def test_default_order(self):
        self.assertEqual(self.domain.default_order, 2)

    def test_minimum_allowed_height(self):
        self.assertAlmostEqual(self.domain.minimum_allowed_height, 1.0e-12)

    def test_edge_coordinates_set(self):
        self.assertIsNotNone(self.domain.edge_coordinates)


class Test_DE1_defaults(unittest.TestCase):

    def setUp(self):
        self.domain = _make_domain()
        self.domain.set_flow_algorithm('DE1')

    def test_cfl(self):
        self.assertAlmostEqual(self.domain.CFL, 1.0)

    def test_timestepping_rk2(self):
        self.assertEqual(self.domain.timestepping_method, 'rk2')

    def test_discontinuous_elevation(self):
        self.assertTrue(self.domain.get_using_discontinuous_elevation())

    def test_compute_fluxes_method(self):
        self.assertEqual(self.domain.get_compute_fluxes_method(), 'DE')

    def test_beta_w(self):
        self.assertAlmostEqual(self.domain.beta_w, 1.0)

    def test_minimum_allowed_height(self):
        self.assertAlmostEqual(self.domain.minimum_allowed_height, 1.0e-5)

    def test_edge_coordinates_set(self):
        self.assertIsNotNone(self.domain.edge_coordinates)


class Test_DE2_defaults(unittest.TestCase):

    def setUp(self):
        self.domain = _make_domain()
        self.domain.set_flow_algorithm('DE2')

    def test_cfl(self):
        self.assertAlmostEqual(self.domain.CFL, 1.0)

    def test_timestepping_rk3(self):
        self.assertEqual(self.domain.timestepping_method, 'rk3')

    def test_beta_w(self):
        self.assertAlmostEqual(self.domain.beta_w, 1.0)

    def test_minimum_allowed_height(self):
        self.assertAlmostEqual(self.domain.minimum_allowed_height, 1.0e-5)

    def test_compute_fluxes_method(self):
        self.assertEqual(self.domain.get_compute_fluxes_method(), 'DE')


class Test_DE0_7_defaults(unittest.TestCase):

    def setUp(self):
        self.domain = _make_domain()
        self.domain.set_flow_algorithm('DE0_7')

    def test_cfl(self):
        self.assertAlmostEqual(self.domain.CFL, 0.9)

    def test_timestepping_euler(self):
        self.assertEqual(self.domain.timestepping_method, 'euler')

    def test_beta_w(self):
        self.assertAlmostEqual(self.domain.beta_w, 0.7)

    def test_beta_w_dry(self):
        self.assertAlmostEqual(self.domain.beta_w_dry, 0.1)

    def test_minimum_allowed_height(self):
        self.assertAlmostEqual(self.domain.minimum_allowed_height, 1.0e-12)

    def test_compute_fluxes_method(self):
        self.assertEqual(self.domain.get_compute_fluxes_method(), 'DE')


class Test_DE1_7_defaults(unittest.TestCase):

    def setUp(self):
        self.domain = _make_domain()
        self.domain.set_flow_algorithm('DE1_7')

    def test_cfl(self):
        self.assertAlmostEqual(self.domain.CFL, 1.0)

    def test_timestepping_rk2(self):
        self.assertEqual(self.domain.timestepping_method, 'rk2')

    def test_beta_w(self):
        self.assertAlmostEqual(self.domain.beta_w, 0.75)

    def test_beta_w_dry(self):
        self.assertAlmostEqual(self.domain.beta_w_dry, 0.1)

    def test_minimum_allowed_height(self):
        self.assertAlmostEqual(self.domain.minimum_allowed_height, 1.0e-12)

    def test_compute_fluxes_method(self):
        self.assertEqual(self.domain.get_compute_fluxes_method(), 'DE')

    def test_edge_coordinates_set(self):
        self.assertIsNotNone(self.domain.edge_coordinates)


class Test_algorithm_verbose(unittest.TestCase):
    """Test that verbose=True doesn't crash the setters."""

    def test_DE0_verbose(self):
        domain = anuga.rectangular_cross_domain(2, 2)
        domain.verbose = True
        domain.set_flow_algorithm('DE0')  # should not raise

    def test_DE1_verbose(self):
        domain = anuga.rectangular_cross_domain(2, 2)
        domain.verbose = True
        domain.set_flow_algorithm('DE1')  # should not raise

    def test_DE2_verbose(self):
        domain = anuga.rectangular_cross_domain(2, 2)
        domain.verbose = True
        domain.set_flow_algorithm('DE2')  # should not raise

    def test_DE0_7_verbose(self):
        domain = anuga.rectangular_cross_domain(2, 2)
        domain.verbose = True
        domain.set_flow_algorithm('DE0_7')  # should not raise

    def test_DE1_7_verbose(self):
        domain = anuga.rectangular_cross_domain(2, 2)
        domain.verbose = True
        domain.set_flow_algorithm('DE1_7')  # should not raise


class Test_domain_misc_methods(unittest.TestCase):
    """Tests for misc. shallow_water Domain methods that are easy to cover."""

    def setUp(self):
        self.domain = _make_domain()

    def test_print_algorithm_parameters(self):
        self.domain.print_algorithm_parameters()  # should not raise

    def test_get_store(self):
        self.domain.set_store(True)
        self.assertTrue(self.domain.get_store())
        self.domain.set_store(False)
        self.assertFalse(self.domain.get_store())

    def test_get_store_centroids(self):
        self.domain.set_store_centroids(True)
        self.assertTrue(self.domain.get_store_centroids())

    def test_get_starttime_as_datetime(self):
        self.domain.set_starttime(0.0)
        dt = self.domain.get_starttime(datetime=True)
        self.assertIsNotNone(dt)

    def test_set_timezone_string(self):
        self.domain.set_timezone('UTC')
        tz = self.domain.get_timezone()
        self.assertIsNotNone(tz)

    def test_set_timezone_none(self):
        self.domain.set_timezone(None)
        tz = self.domain.get_timezone()
        self.assertIsNotNone(tz)

    def test_set_timezone_unknown_raises(self):
        with self.assertRaises(Exception):
            self.domain.set_timezone(12345)

    def test_set_checkpointing(self):
        import tempfile
        import os
        with tempfile.TemporaryDirectory() as tmpdir:
            ckdir = os.path.join(tmpdir, 'checkpoints')
            self.domain.set_checkpointing(checkpoint=True,
                                          checkpoint_dir=ckdir,
                                          checkpoint_step=5)
            self.assertTrue(os.path.exists(ckdir))

    def test_set_checkpointing_false(self):
        self.domain.set_checkpointing(checkpoint=False)  # should not raise

    def test_set_checkpointing_with_time(self):
        import tempfile
        import os
        with tempfile.TemporaryDirectory() as tmpdir:
            ckdir = os.path.join(tmpdir, 'ck_t')
            self.domain.set_checkpointing(checkpoint=True,
                                          checkpoint_dir=ckdir,
                                          checkpoint_time=60.0)
            self.assertEqual(self.domain.checkpoint_time, 60.0)
            self.assertEqual(self.domain.checkpoint_step, 0)

    def test_set_compute_fluxes_method_invalid(self):
        with self.assertRaises(Exception):
            self.domain.set_compute_fluxes_method('INVALID')

    def test_set_extrapolate_velocity_false(self):
        self.domain.set_extrapolate_velocity(False)
        self.assertFalse(self.domain.extrapolate_velocity_second_order)

    def test_set_use_optimise_dry_cells_false(self):
        self.domain.set_use_optimise_dry_cells(False)
        self.assertEqual(self.domain.optimise_dry_cells, 0)

    def test_set_use_kinematic_viscosity_true(self):
        self.domain.set_use_kinematic_viscosity(True)
        self.assertIsNotNone(self.domain.kv_operator)

    def test_set_use_kinematic_viscosity_toggle(self):
        self.domain.set_use_kinematic_viscosity(True)
        self.assertIsNotNone(self.domain.kv_operator)
        self.domain.set_use_kinematic_viscosity(False)
        self.assertIsNone(self.domain.kv_operator)

    def test_set_betas(self):
        self.domain.set_betas(0.8, 0.1, 0.7, 0.05, 0.6, 0.0)
        self.assertAlmostEqual(self.domain.beta_w, 0.8)
        self.assertAlmostEqual(self.domain.beta_w_dry, 0.1)
        self.assertAlmostEqual(self.domain.beta_uh, 0.7)
        self.assertAlmostEqual(self.domain.beta_vh, 0.6)

    def test_get_minimum_storable_height(self):
        self.domain.set_minimum_storable_height(0.01)
        self.assertAlmostEqual(self.domain.get_minimum_storable_height(), 0.01)

    def test_set_store_vertices_uniquely(self):
        # set_store_vertices_uniquely(True) means NOT smooth
        self.domain.set_store_vertices_uniquely(True)
        self.assertFalse(self.domain.smooth)

    def test_set_store_vertices_uniquely_with_reduction(self):
        import numpy as np
        self.domain.set_store_vertices_uniquely(True, reduction=min)

    def test_get_global_wet_element_count(self):
        count = self.domain.get_global_wet_element_count()
        self.assertIsInstance(count, (int, float))

    def test_get_global_max_stage(self):
        self.domain.set_quantity('stage', 0.5)
        max_s = self.domain.get_global_max_stage()
        self.assertGreaterEqual(max_s, 0.5)

    def test_get_global_max_speed(self):
        max_speed = self.domain.get_global_max_speed()
        self.assertIsNotNone(max_speed)

    def test_get_water_volume(self):
        v = self.domain.get_water_volume()
        self.assertIsNotNone(v)

    def test_get_fractional_step_volume_integral(self):
        v = self.domain.get_fractional_step_volume_integral()
        self.assertIsNotNone(v)

    def test_get_boundary_flux_integral_raises_for_non_DE(self):
        """get_boundary_flux_integral raises for non-DE flux methods."""
        self.domain.set_compute_fluxes_method('original')
        with self.assertRaises(Exception):
            self.domain.get_boundary_flux_integral()


if __name__ == '__main__':
    unittest.main()
