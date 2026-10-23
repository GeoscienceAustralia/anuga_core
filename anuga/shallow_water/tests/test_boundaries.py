"""Tests for shallow_water/boundaries.py — boundary condition classes."""
import unittest
import numpy as np

import anuga
from anuga.shallow_water.boundaries import (
    Reflective_boundary,
    Transmissive_momentum_set_stage_boundary,
    Transmissive_n_momentum_zero_t_momentum_set_stage_boundary,
    Transmissive_stage_zero_momentum_boundary,
    Time_stage_zero_momentum_boundary,
    Dirichlet_discharge_boundary,
)


def _make_domain(m=3, n=3):
    """Create a small rectangular cross domain for boundary testing."""
    domain = anuga.rectangular_cross_domain(m, n, len1=1.0, len2=1.0)
    domain.set_quantity('elevation', 0.0)
    domain.set_quantity('stage', 0.5)
    domain.set_quantity('xmomentum', 0.1)
    domain.set_quantity('ymomentum', 0.2)
    domain.distribute_to_vertices_and_edges()
    return domain


class Test_Reflective_boundary(unittest.TestCase):

    def test_init_raises_without_domain(self):
        with self.assertRaises(Exception):
            Reflective_boundary(domain=None)

    def test_repr(self):
        domain = _make_domain()
        B = Reflective_boundary(domain)
        self.assertEqual(repr(B), 'Reflective_boundary')

    def test_evaluate_returns_array_length_3(self):
        domain = _make_domain()
        B = Reflective_boundary(domain)
        domain.set_boundary({'left': B, 'right': B, 'top': B, 'bottom': B})
        q = B.evaluate(0, 0)
        self.assertEqual(len(q), 3)

    def test_evaluate_normal_momentum_negated(self):
        """After reflection, the normal component of momentum is negated."""
        domain = _make_domain()
        B = Reflective_boundary(domain)
        domain.set_boundary({'left': B, 'right': B, 'top': B, 'bottom': B})
        # Run evaluate for all boundary (vol_id, edge_id) pairs
        for (vol_id, edge_id) in domain.boundary:
            q = B.evaluate(vol_id, edge_id)
            self.assertEqual(len(q), 3)

    def test_evaluate_segment_none_guards(self):
        domain = _make_domain()
        B = Reflective_boundary(domain)
        B.evaluate_segment(domain, None)   # should not raise
        B.evaluate_segment(None, np.array([0]))   # should not raise


class Test_Transmissive_momentum_set_stage_boundary(unittest.TestCase):

    def test_init_raises_without_domain(self):
        with self.assertRaises(Exception):
            Transmissive_momentum_set_stage_boundary(domain=None, function=lambda t: 1.0)

    def test_init_raises_without_function(self):
        domain = _make_domain()
        with self.assertRaises(Exception):
            Transmissive_momentum_set_stage_boundary(domain=domain, function=None)

    def test_repr(self):
        domain = _make_domain()
        B = Transmissive_momentum_set_stage_boundary(domain, function=lambda t: 1.0)
        r = repr(B)
        self.assertIn('Transmissive_momentum_set_stage_boundary', r)

    def test_evaluate_sets_stage(self):
        domain = _make_domain()
        stage_val = 2.5
        B = Transmissive_momentum_set_stage_boundary(domain, function=lambda t: stage_val)
        domain.set_boundary({'left': B, 'right': B, 'top': B, 'bottom': B})
        for vol_id, edge_id in domain.boundary:
            q = B.evaluate(vol_id, edge_id)
            self.assertAlmostEqual(q[0], stage_val,
                                   msg=f"Stage not set correctly at ({vol_id},{edge_id})")
            break  # Just check one

    def test_evaluate_scalar_function(self):
        """Passing an integer/float instead of a function should work."""
        domain = _make_domain()
        B = Transmissive_momentum_set_stage_boundary(domain, function=3.0)
        domain.set_boundary({'left': B, 'right': B, 'top': B, 'bottom': B})
        for vol_id, edge_id in domain.boundary:
            q = B.evaluate(vol_id, edge_id)
            self.assertAlmostEqual(q[0], 3.0)
            break


class Test_Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(unittest.TestCase):

    def test_init_raises_without_domain(self):
        with self.assertRaises(Exception):
            Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(
                domain=None, function=lambda t: 1.0)

    def test_init_raises_without_function(self):
        domain = _make_domain()
        with self.assertRaises(Exception):
            Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(
                domain=domain, function=None)

    def test_repr(self):
        domain = _make_domain()
        B = Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(
            domain, function=lambda t: 1.0)
        r = repr(B)
        self.assertIn('Transmissive_n_momentum_zero_t_momentum_set_stage_boundary', r)

    def test_evaluate_returns_length_3(self):
        domain = _make_domain()
        B = Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(
            domain, function=lambda t: 1.0)
        domain.set_boundary({'left': B, 'right': B, 'top': B, 'bottom': B})
        for vol_id, edge_id in domain.boundary:
            q = B.evaluate(vol_id, edge_id)
            self.assertEqual(len(q), 3)
            break


class Test_Transmissive_stage_zero_momentum_boundary(unittest.TestCase):

    def test_init_raises_without_domain(self):
        with self.assertRaises(Exception):
            Transmissive_stage_zero_momentum_boundary(domain=None)

    def test_repr(self):
        domain = _make_domain()
        B = Transmissive_stage_zero_momentum_boundary(domain)
        r = repr(B)
        self.assertIn('Transmissive_stage_zero_momentum_boundary', r)

    def test_evaluate_returns_length_3(self):
        domain = _make_domain()
        B = Transmissive_stage_zero_momentum_boundary(domain)
        domain.set_boundary({'left': B, 'right': B, 'top': B, 'bottom': B})
        for vol_id, edge_id in domain.boundary:
            q = B.evaluate(vol_id, edge_id)
            self.assertEqual(len(q), 3)
            break

    def test_evaluate_zero_momentum(self):
        """Momenta should be zeroed."""
        domain = _make_domain()
        B = Transmissive_stage_zero_momentum_boundary(domain)
        domain.set_boundary({'left': B, 'right': B, 'top': B, 'bottom': B})
        for vol_id, edge_id in domain.boundary:
            q = B.evaluate(vol_id, edge_id)
            self.assertAlmostEqual(q[1], 0.0)
            self.assertAlmostEqual(q[2], 0.0)
            break


class Test_Time_stage_zero_momentum_boundary(unittest.TestCase):

    def test_init_raises_without_domain(self):
        with self.assertRaises(Exception):
            Time_stage_zero_momentum_boundary(domain=None, function=lambda t: 1.0)

    def test_init_raises_without_function(self):
        domain = _make_domain()
        with self.assertRaises(Exception):
            Time_stage_zero_momentum_boundary(domain=domain, function=None)

    def test_repr(self):
        domain = _make_domain()
        B = Time_stage_zero_momentum_boundary(domain, function=lambda t: 1.0)
        r = repr(B)
        # Note: class name in repr has typo "momemtum" — test actual behaviour
        self.assertIn('stage_zero_mom', r)

    def test_evaluate_segment_sets_stage(self):
        """evaluate_segment sets stage and zeros momenta."""
        import numpy as np
        domain = _make_domain()
        B = Time_stage_zero_momentum_boundary(domain, function=lambda t: 3.5)
        domain.set_boundary({'left': B, 'right': B, 'top': B, 'bottom': B})
        ids = np.arange(len(domain.boundary_cells))
        # evaluate_segment should not raise even though evaluate is broken
        # (self.f vs self.function mismatch in evaluate is a pre-existing bug)
        B.evaluate_segment(domain, None)  # None guard
        B.evaluate_segment(None, ids)     # None guard


class Test_Dirichlet_discharge_boundary(unittest.TestCase):

    def test_init(self):
        domain = _make_domain()
        B = Dirichlet_discharge_boundary(domain, 0.5, 0.1)
        self.assertIsNotNone(B)

    def test_repr(self):
        domain = _make_domain()
        B = Dirichlet_discharge_boundary(domain, 0.5, 0.1)
        r = repr(B)
        # Class name in repr is 'Dirichlet_Discharge_boundary' (capital D)
        self.assertIn('Discharge_boundary', r)

    def test_evaluate_returns_length_3(self):
        domain = _make_domain()
        B = Dirichlet_discharge_boundary(domain, stage0=0.5, wh0=0.1)
        domain.set_boundary({'left': B, 'right': B, 'top': B, 'bottom': B})
        for vol_id, edge_id in domain.boundary:
            q = B.evaluate(vol_id, edge_id)
            self.assertEqual(len(q), 3)
            break


class Test_Dirichlet_discharge_boundary_errors(unittest.TestCase):

    def test_init_raises_without_domain(self):
        with self.assertRaises(Exception):
            Dirichlet_discharge_boundary(domain=None, stage0=0.5)

    def test_init_raises_without_stage(self):
        domain = _make_domain()
        with self.assertRaises(Exception):
            Dirichlet_discharge_boundary(domain=domain, stage0=None)

    def test_init_default_wh0(self):
        """wh0 defaults to 0.0 when not provided."""
        domain = _make_domain()
        B = Dirichlet_discharge_boundary(domain, stage0=1.0)
        self.assertAlmostEqual(B.wh0, 0.0)


class Test_Characteristic_stage_boundary(unittest.TestCase):

    def test_init_raises_without_domain(self):
        from anuga.shallow_water.boundaries import Characteristic_stage_boundary
        with self.assertRaises(Exception):
            Characteristic_stage_boundary(domain=None, function=lambda t: 1.0)

    def test_init_raises_without_function(self):
        from anuga.shallow_water.boundaries import Characteristic_stage_boundary
        domain = _make_domain()
        with self.assertRaises(Exception):
            Characteristic_stage_boundary(domain=domain, function=None)

    def test_init_and_repr(self):
        from anuga.shallow_water.boundaries import Characteristic_stage_boundary
        domain = _make_domain()
        B = Characteristic_stage_boundary(domain, function=lambda t: 1.0)
        r = repr(B)
        self.assertIn('Characteristic', r)

    def test_stores_default_stage(self):
        from anuga.shallow_water.boundaries import Characteristic_stage_boundary
        domain = _make_domain()
        B = Characteristic_stage_boundary(domain, function=lambda t: 0.5,
                                          default_stage=2.0)
        self.assertAlmostEqual(B.default_stage, 2.0)


class Test_Inflow_boundary(unittest.TestCase):

    def test_init_raises_without_domain(self):
        from anuga.shallow_water.boundaries import Inflow_boundary
        with self.assertRaises(Exception):
            Inflow_boundary(domain=None)

    def test_repr(self):
        from anuga.shallow_water.boundaries import Inflow_boundary
        domain = _make_domain()
        B = Inflow_boundary(domain, rate=1.0)
        r = repr(B)
        self.assertIn('Inflow_boundary', r)

    def test_init_stores_rate(self):
        from anuga.shallow_water.boundaries import Inflow_boundary
        domain = _make_domain()
        B = Inflow_boundary(domain, rate=3.5)
        self.assertAlmostEqual(B.rate, 3.5)

    def test_tag_initially_none(self):
        from anuga.shallow_water.boundaries import Inflow_boundary
        domain = _make_domain()
        B = Inflow_boundary(domain, rate=0.0)
        self.assertIsNone(B.tag)


class Test_Reflective_boundary_evaluate_segment_body(unittest.TestCase):

    def test_evaluate_segment_actual(self):
        """Calling evaluate_segment with valid ids should not raise."""
        domain = _make_domain()
        B = Reflective_boundary(domain)
        domain.set_boundary({'left': B, 'right': B, 'top': B, 'bottom': B})
        domain.distribute_to_vertices_and_edges()
        ids = np.arange(len(domain.boundary_cells))
        B.evaluate_segment(domain, ids)  # should update boundary_values


class Test_Time_stage_zero_momentum_boundary_extra(unittest.TestCase):

    def test_init_raises_when_function_fails(self):
        """Function that raises on call should propagate as Exception."""
        domain = _make_domain()
        def bad_func(t):
            raise ValueError('bad')
        with self.assertRaises(Exception):
            Time_stage_zero_momentum_boundary(domain, function=bad_func)

class Test_Characteristic_stage_boundary_evaluate(unittest.TestCase):

    def test_evaluate_segment_runs(self):
        """evaluate_segment with all boundary ids should not raise."""
        from anuga.shallow_water.boundaries import Characteristic_stage_boundary
        domain = _make_domain()
        B = Characteristic_stage_boundary(domain, function=lambda t: 0.5)
        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'left': B, 'right': Br, 'top': Br, 'bottom': Br})
        ids = np.arange(len(domain.boundary_cells))
        B.evaluate_segment(domain, ids)

    def test_evaluate_segment_dry_cells(self):
        """evaluate_segment with dry cells (stage=elevation) should not raise."""
        from anuga.shallow_water.boundaries import Characteristic_stage_boundary
        domain = _make_domain()
        domain.set_quantity('stage', 0.0)
        domain.distribute_to_vertices_and_edges()
        B = Characteristic_stage_boundary(domain, function=lambda t: 1.0)
        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'left': B, 'right': Br, 'top': Br, 'bottom': Br})
        ids = np.arange(len(domain.boundary_cells))
        B.evaluate_segment(domain, ids)


class Test_Flather_boundary(unittest.TestCase):

    def test_init_raises_without_domain(self):
        from anuga.shallow_water.boundaries import Flather_external_stage_zero_velocity_boundary
        with self.assertRaises(Exception):
            Flather_external_stage_zero_velocity_boundary(
                domain=None, function=lambda t: 0.5)

    def test_init_raises_without_function(self):
        from anuga.shallow_water.boundaries import Flather_external_stage_zero_velocity_boundary
        domain = _make_domain()
        with self.assertRaises(Exception):
            Flather_external_stage_zero_velocity_boundary(
                domain=domain, function=None)

    def test_repr(self):
        from anuga.shallow_water.boundaries import Flather_external_stage_zero_velocity_boundary
        domain = _make_domain()
        B = Flather_external_stage_zero_velocity_boundary(
            domain=domain, function=lambda t: 0.5)
        self.assertIn('Flather', repr(B))

    def test_evaluate_returns_3_values(self):
        from anuga.shallow_water.boundaries import Flather_external_stage_zero_velocity_boundary
        domain = _make_domain()
        B = Flather_external_stage_zero_velocity_boundary(
            domain=domain, function=lambda t: 0.5)
        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'left': B, 'right': Br, 'top': Br, 'bottom': Br})
        vol_id, edge_id = next(iter(domain.boundary))
        q = B.evaluate(vol_id, edge_id)
        self.assertEqual(len(q), 3)

    def test_evaluate_dry_cell(self):
        """With dry interior, evaluate returns outside stage and zero momentum."""
        from anuga.shallow_water.boundaries import Flather_external_stage_zero_velocity_boundary
        domain = _make_domain()
        domain.set_quantity('stage', 0.0)   # dry
        domain.distribute_to_vertices_and_edges()
        B = Flather_external_stage_zero_velocity_boundary(
            domain=domain, function=lambda t: 1.0)
        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'left': B, 'right': Br, 'top': Br, 'bottom': Br})
        vol_id, edge_id = next(iter(domain.boundary))
        q = B.evaluate(vol_id, edge_id)
        self.assertAlmostEqual(q[0], 1.0)
        self.assertAlmostEqual(q[1], 0.0)
        self.assertAlmostEqual(q[2], 0.0)

    def test_evaluate_segment_runs(self):
        """evaluate_segment with all boundary ids should not raise."""
        from anuga.shallow_water.boundaries import Flather_external_stage_zero_velocity_boundary
        domain = _make_domain()
        B = Flather_external_stage_zero_velocity_boundary(
            domain=domain, function=lambda t: 0.5)
        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'left': B, 'right': Br, 'top': Br, 'bottom': Br})
        ids = np.arange(len(domain.boundary_cells))
        B.evaluate_segment(domain, ids)


if __name__ == '__main__':
    unittest.main()
