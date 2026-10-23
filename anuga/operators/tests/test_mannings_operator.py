"""Tests for Mannings_operator."""

import unittest
import numpy as np
import anuga
from anuga.operators.mannings_operator import Mannings_operator


def _make_domain(friction=0.03):
    domain = anuga.rectangular_cross_domain(4, 4, len1=4.0, len2=4.0)
    domain.set_flow_algorithm('DE0')
    domain.set_quantity('elevation', 0.0)
    domain.set_quantity('stage', 1.0)
    domain.set_quantity('xmomentum', 0.5)
    domain.set_quantity('ymomentum', 0.0)
    domain.set_quantity('friction', friction)
    domain.set_boundary({'left': anuga.Reflective_boundary(domain),
                         'right': anuga.Reflective_boundary(domain),
                         'top': anuga.Reflective_boundary(domain),
                         'bottom': anuga.Reflective_boundary(domain)})
    return domain


class Test_Mannings_operator(unittest.TestCase):

    def test_init(self):
        domain = _make_domain()
        op = Mannings_operator(domain)
        self.assertIsNotNone(op)
        self.assertEqual(op.exp_gamma_max, 0.0)
        self.assertEqual(op.exp_gamma_min, 1.0)

    def test_init_verbose(self):
        domain = _make_domain()
        op = Mannings_operator(domain, verbose=True)
        self.assertIsNotNone(op)

    def test_call_reduces_momentum(self):
        """Manning friction should reduce momentum magnitudes."""
        domain = _make_domain(friction=0.05)
        op = Mannings_operator(domain)
        xmom_before = domain.quantities['xmomentum'].centroid_values.copy()
        domain.timestep = 0.1
        op()
        xmom_after = domain.quantities['xmomentum'].centroid_values
        # All x-momentum should be reduced in magnitude (friction damps)
        np.testing.assert_array_less(np.abs(xmom_after), np.abs(xmom_before) + 1e-12)

    def test_call_zero_friction(self):
        """Zero friction should leave momentum unchanged."""
        domain = _make_domain(friction=0.0)
        op = Mannings_operator(domain)
        xmom_before = domain.quantities['xmomentum'].centroid_values.copy()
        domain.timestep = 0.1
        op()
        xmom_after = domain.quantities['xmomentum'].centroid_values
        np.testing.assert_allclose(xmom_after, xmom_before, atol=1e-12)

    def test_parallel_safe(self):
        domain = _make_domain()
        op = Mannings_operator(domain)
        self.assertTrue(op.parallel_safe())

    def test_statistics(self):
        domain = _make_domain()
        op = Mannings_operator(domain)
        msg = op.statistics()
        self.assertIn('Manning', msg)

    def test_timestepping_statistics(self):
        domain = _make_domain(friction=0.03)
        op = Mannings_operator(domain)
        domain.timestep = 0.1
        op()
        msg = op.timestepping_statistics()
        self.assertIn('Manning', msg)
        # Stats should reset after call
        self.assertEqual(op.exp_gamma_max, 0.0)
        self.assertEqual(op.exp_gamma_min, 1.0)

    def test_dry_cells_handled(self):
        """Cells with zero depth should not cause divide-by-zero."""
        domain = _make_domain()
        domain.set_quantity('stage', -1.0)   # dry everywhere
        domain.set_quantity('elevation', 0.0)
        op = Mannings_operator(domain)
        domain.timestep = 0.1
        op()   # should not raise


if __name__ == '__main__':
    unittest.main()
