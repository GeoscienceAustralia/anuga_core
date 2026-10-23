"""Tests for the ADER-2 (Cauchy-Kovalewski) timestepping scheme."""

import unittest
import numpy as np
import anuga
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
from anuga.shallow_water.shallow_water_domain import Domain


def _make_domain(nx=10, ny=10, length=100.0):
    """Small rectangular domain for fast unit tests."""
    points, vertices, boundary = rectangular_cross(nx, ny, len1=length, len2=length)
    domain = Domain(points, vertices, boundary)
    domain.set_store(False)
    domain.set_default_order(2)
    Br = anuga.Reflective_boundary(domain)
    domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
    return domain


class TestAder2Setup(unittest.TestCase):

    def test_set_timestepping_method_string(self):
        """domain.set_timestepping_method('ader2') is accepted."""
        domain = _make_domain()
        domain.set_timestepping_method('ader2')
        self.assertEqual(domain.get_timestepping_method(), 'ader2')

    def test_timestep_fluxcalls_is_one(self):
        """ader2 uses timestep_fluxcalls == 1."""
        domain = _make_domain()
        domain.set_timestepping_method('ader2')
        self.assertEqual(domain.timestep_fluxcalls, 1)

    def test_invalid_method_still_raises(self):
        """Unrecognised timestepping method raises."""
        domain = _make_domain()
        with self.assertRaises(Exception):
            domain.set_timestepping_method('bogus')


class TestAder2WellBalance(unittest.TestCase):

    def test_still_water_flat_bed(self):
        """Still water over flat bed stays still (well-balanced)."""
        domain = _make_domain()
        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('stage', 1.0)
        domain.set_timestepping_method('ader2')

        for t in domain.evolve(yieldstep=2.0, finaltime=4.0):
            pass

        stage = domain.quantities['stage'].centroid_values
        xmom  = domain.quantities['xmomentum'].centroid_values
        ymom  = domain.quantities['ymomentum'].centroid_values
        self.assertTrue(np.allclose(stage, 1.0, atol=1e-14),
                        f'stage drifted: max dev = {np.max(np.abs(stage-1.0)):.2e}')
        self.assertTrue(np.allclose(xmom, 0.0, atol=1e-14))
        self.assertTrue(np.allclose(ymom, 0.0, atol=1e-14))

    def test_still_water_sloped_bed(self):
        """Still water over sloped bed stays still (well-balanced)."""
        domain = _make_domain()
        domain.set_quantity('elevation', lambda x, y: 0.01 * x)
        domain.set_quantity('stage', 1.0)  # uniform stage = sloped free surface
        domain.set_timestepping_method('ader2')

        for t in domain.evolve(yieldstep=2.0, finaltime=4.0):
            pass

        stage = domain.quantities['stage'].centroid_values
        xmom  = domain.quantities['xmomentum'].centroid_values
        self.assertTrue(np.allclose(stage, 1.0, atol=1e-10),
                        f'stage drifted: max dev = {np.max(np.abs(stage-1.0)):.2e}')
        self.assertTrue(np.allclose(xmom, 0.0, atol=1e-10))


class TestAder2MassConservation(unittest.TestCase):

    def test_dam_break_mass_conserved(self):
        """Total volume is conserved under ADER-2 with reflective walls."""
        domain = _make_domain(nx=15, ny=15)

        def stage(x, y):
            return np.where(x > 50, 0.5, 1.0)

        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('stage', stage)
        domain.set_timestepping_method('ader2')

        vol0 = domain.compute_total_volume()
        for t in domain.evolve(yieldstep=5.0, finaltime=5.0):
            pass
        vol1 = domain.compute_total_volume()

        self.assertAlmostEqual(vol0, vol1, delta=vol0 * 1e-12)

    def test_mass_same_as_rk2(self):
        """ADER-2 and RK2 conserve the same total mass."""
        def run(method):
            domain = _make_domain(nx=12, ny=12)
            domain.set_quantity('elevation', 0.0)
            domain.set_quantity('stage', lambda x, y: np.where(x > 50, 0.5, 1.0))
            domain.set_timestepping_method(method)
            for t in domain.evolve(yieldstep=5.0, finaltime=5.0):
                pass
            return domain.compute_total_volume()

        self.assertAlmostEqual(run('ader2'), run('rk2'), delta=1e-8)


class TestAder2Consistency(unittest.TestCase):

    def test_dam_break_stage_in_range(self):
        """Stage stays within [h_right, h_left] during a dam break."""
        domain = _make_domain(nx=15, ny=15)
        h_left, h_right = 1.0, 0.5

        def stage(x, y):
            return np.where(x > 50, h_right, h_left)

        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('stage', stage)
        domain.set_timestepping_method('ader2')

        for t in domain.evolve(yieldstep=5.0, finaltime=5.0):
            pass

        stage_cv = domain.quantities['stage'].centroid_values
        # Allow small reconstruction overshoot (second-order scheme)
        self.assertTrue(np.all(stage_cv >= h_right - 1e-6))
        self.assertTrue(np.all(stage_cv <= h_left + 1e-6))

    def test_ader2_close_to_rk2(self):
        """ADER-2 and RK2 give similar final states (max diff < 0.05)."""
        def run(method):
            domain = _make_domain(nx=15, ny=15)
            domain.set_quantity('elevation', 0.0)
            domain.set_quantity('stage', lambda x, y: np.where(x > 50, 0.5, 1.0))
            domain.set_timestepping_method(method)
            for t in domain.evolve(yieldstep=5.0, finaltime=5.0):
                pass
            return domain.quantities['stage'].centroid_values.copy()

        s_rk2   = run('rk2')
        s_ader2 = run('ader2')
        max_diff = np.max(np.abs(s_rk2 - s_ader2))
        self.assertLess(max_diff, 0.05,
                        f'ADER-2 and RK2 diverged: max diff = {max_diff:.4f}')


class TestAder2NonNegative(unittest.TestCase):

    def test_no_negative_depths(self):
        """Depths remain non-negative under ADER-2."""
        domain = _make_domain(nx=15, ny=15)

        def stage(x, y):
            return np.where(x > 50, 0.5, 1.0)

        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('stage', stage)
        domain.set_timestepping_method('ader2')

        for t in domain.evolve(yieldstep=5.0, finaltime=5.0):
            pass

        stage_cv = domain.quantities['stage'].centroid_values
        elev_cv  = domain.quantities['elevation'].centroid_values
        depths   = stage_cv - elev_cv
        self.assertTrue(np.all(depths >= -1e-10),
                        f'Negative depth detected: min depth = {depths.min():.2e}')


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromModule(__import__(__name__))
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
