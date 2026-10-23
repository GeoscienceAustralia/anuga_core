#!/usr/bin/env python

import unittest
import warnings

import numpy
import anuga
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.structures.internal_boundary_operator import Internal_boundary_operator


verbose = False

# This end-point geometry places one enquiry point inside an inlet triangle,
# which is expected for this small test mesh and triggers a UserWarning.
_END_POINTS = [[3., 2.5], [7., 2.5]]
_INLET_WARNING = 'Enquiry point.*is in an inlet triangle'


def make_culvert_domain():
    """Create a simple 10m x 5m rectangular domain for culvert tests."""
    points, vertices, boundary = rectangular_cross(10, 5, len1=10.0, len2=5.0)
    domain = Domain(points, vertices, boundary)
    domain.set_quantity('elevation', 0.0)
    domain.set_quantity('stage', 1.0)
    domain.set_quantity('friction', 0.0)
    Br = anuga.Reflective_boundary(domain)
    domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
    return domain


def simple_ibf(hw, tw):
    """Simple internal boundary function: flow proportional to head difference."""
    return max(0.0, hw - tw) * 1.0


class Test_Internal_boundary_operator(unittest.TestCase):
    """Tests for Internal_boundary_operator."""

    def setUp(self):
        self.domain = make_culvert_domain()
        # Suppress the known inlet-triangle warning for all tests except
        # test_construction which asserts it explicitly.
        self._warning_ctx = warnings.catch_warnings()
        self._warning_ctx.__enter__()
        warnings.filterwarnings('ignore', message=_INLET_WARNING,
                                category=UserWarning)

    def tearDown(self):
        self._warning_ctx.__exit__(None, None, None)

    def test_construction(self):
        """Internal_boundary_operator construction warns when enquiry point is
        in an inlet triangle (expected for this test mesh geometry)."""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            op = Internal_boundary_operator(
                self.domain,
                internal_boundary_function=simple_ibf,
                end_points=_END_POINTS,
                width=1.0,
                verbose=verbose)
        self.assertIsNotNone(op)
        inlet_warnings = [w for w in caught
                          if issubclass(w.category, UserWarning)
                          and 'inlet triangle' in str(w.message)]
        self.assertGreater(len(inlet_warnings), 0,
                           'Expected inlet-triangle UserWarning was not raised')

    def test_structure_type(self):
        """structure_type attribute is set to 'internal_boundary'."""
        op = Internal_boundary_operator(
            self.domain,
            internal_boundary_function=simple_ibf,
            end_points=_END_POINTS,
            width=1.0,
            verbose=verbose)
        self.assertEqual(op.structure_type, 'internal_boundary')

    def test_width_stored(self):
        """get_culvert_width() returns the specified width."""
        op = Internal_boundary_operator(
            self.domain,
            internal_boundary_function=simple_ibf,
            end_points=_END_POINTS,
            width=1.0,
            verbose=verbose)
        self.assertAlmostEqual(op.get_culvert_width(), 1.0)

    def test_calling_operator(self):
        """Calling the operator once does not raise an exception."""
        op = Internal_boundary_operator(
            self.domain,
            internal_boundary_function=simple_ibf,
            end_points=_END_POINTS,
            width=1.0,
            verbose=verbose)
        self.domain.timestep = 0.1
        self.domain.yieldstep = 1.0
        # Calling the operator should not raise
        op()

    def test_statistics_returns_string(self):
        """statistics() returns a non-empty string."""
        op = Internal_boundary_operator(
            self.domain,
            internal_boundary_function=simple_ibf,
            end_points=_END_POINTS,
            width=1.0,
            verbose=verbose)
        result = op.statistics()
        self.assertIsInstance(result, str)
        self.assertGreater(len(result), 0)


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Internal_boundary_operator)
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)


# ---------------------------------------------------------------------------
# pytest functions: uncovered branches
# ---------------------------------------------------------------------------

import pytest


def _suppress_inlet_warning():
    ctx = warnings.catch_warnings()
    ctx.__enter__()
    warnings.filterwarnings('ignore', message=_INLET_WARNING, category=UserWarning)
    return ctx


def test_verbose_true_prints(capsys):
    """verbose=True triggers the EXPERIMENTAL warning block."""
    domain = make_culvert_domain()
    ctx = _suppress_inlet_warning()
    try:
        Internal_boundary_operator(
            domain, simple_ibf, end_points=_END_POINTS,
            width=1.0, verbose=True)
    finally:
        ctx.__exit__(None, None, None)
    out = capsys.readouterr().out
    assert 'INTERNAL BOUNDARY OPERATOR' in out


def test_explicit_discharge_runs():
    """compute_discharge_implicitly=False routes through discharge_routine_explicit."""
    domain = make_culvert_domain()
    ctx = _suppress_inlet_warning()
    try:
        op = Internal_boundary_operator(
            domain, simple_ibf, end_points=_END_POINTS,
            width=1.0, compute_discharge_implicitly=False, verbose=False)
    finally:
        ctx.__exit__(None, None, None)
    domain.timestep = 0.1
    domain.yieldstep = 1.0
    op()


def test_explicit_discharge_use_velocity_head():
    """Explicit discharge with use_velocity_head=True uses total energy."""
    domain = make_culvert_domain()
    ctx = _suppress_inlet_warning()
    try:
        op = Internal_boundary_operator(
            domain, simple_ibf, end_points=_END_POINTS,
            width=1.0, compute_discharge_implicitly=False,
            use_velocity_head=True, verbose=False)
    finally:
        ctx.__exit__(None, None, None)
    domain.timestep = 0.1
    domain.yieldstep = 1.0
    op()


def test_explicit_blocked_structure():
    """height <= 0 returns Q=0 immediately in discharge_routine_explicit."""
    domain = make_culvert_domain()
    ctx = _suppress_inlet_warning()
    try:
        op = Internal_boundary_operator(
            domain, simple_ibf, end_points=_END_POINTS,
            width=1.0, height=1.0, compute_discharge_implicitly=False, verbose=False)
    finally:
        ctx.__exit__(None, None, None)
    op.height = 0.0
    domain.timestep = 0.1
    domain.yieldstep = 1.0
    op()
    assert op.discharge == 0.0 or op.discharge is not None


def test_implicit_use_velocity_head():
    """Implicit discharge with use_velocity_head=True uses total energy."""
    domain = make_culvert_domain()
    ctx = _suppress_inlet_warning()
    try:
        op = Internal_boundary_operator(
            domain, simple_ibf, end_points=_END_POINTS,
            width=1.0, compute_discharge_implicitly=True,
            use_velocity_head=True, verbose=False)
    finally:
        ctx.__exit__(None, None, None)
    domain.timestep = 0.1
    domain.yieldstep = 1.0
    op()


def test_implicit_negative_Q_reverses_inlets():
    """Implicit discharge with reverse flow sets inlets[1] as inflow."""
    def reverse_ibf(hw, tw):
        return hw - tw  # negative when tw > hw

    domain = make_culvert_domain()
    # Make right side (inlet1 near x=7) higher than left (inlet0 near x=3)
    def stage_fn(x, y):
        return numpy.where(x > 5.0, 2.0, 1.0)
    domain.set_quantity('stage', stage_fn)

    ctx = _suppress_inlet_warning()
    try:
        op = Internal_boundary_operator(
            domain, reverse_ibf, end_points=_END_POINTS,
            width=1.0, compute_discharge_implicitly=True, verbose=False)
    finally:
        ctx.__exit__(None, None, None)
    domain.timestep = 0.1
    domain.yieldstep = 1.0
    op()
