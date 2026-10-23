"""Tests for Wind_stress_operator and Barometric_pressure_operator."""

import pytest
import numpy as np
import anuga
from anuga import Domain, Reflective_boundary
from anuga.config import rho_w, rho_a, eta_w
from math import pi


def _make_domain(stage=1.0):
    """6-point, 4-triangle flat domain at given depth, timestep=1 s."""
    points = [[0, 0], [0, 2], [2, 0], [0, 4], [2, 2], [4, 0]]
    vertices = [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]]
    domain = Domain(points, vertices)
    domain.set_quantity('elevation', 0.0)
    domain.set_quantity('stage', stage)
    domain.set_boundary({'exterior': Reflective_boundary(domain)})
    domain.timestep = 1.0
    return domain


# ---------------------------------------------------------------------------
# Barometric_pressure_operator
# ---------------------------------------------------------------------------

class TestBarometricPressureOperator:

    def test_constant_pressure_no_update(self):
        """Uniform pressure has zero gradient — momentum must not change."""
        from anuga.operators.barometric_pressure import Barometric_pressure_operator
        domain = _make_domain()
        P = Barometric_pressure_operator(domain, pressure=101325.0)
        xmom0 = domain.quantities['xmomentum'].centroid_values.copy()
        ymom0 = domain.quantities['ymomentum'].centroid_values.copy()
        P()
        np.testing.assert_allclose(
            domain.quantities['xmomentum'].centroid_values, xmom0, atol=1e-14)
        np.testing.assert_allclose(
            domain.quantities['ymomentum'].centroid_values, ymom0, atol=1e-14)

    def test_callable_pressure_updates_xmom(self):
        """p(t,x,y)=x gives dp/dx=1, dp/dy=0: only xmom increases."""
        from anuga.operators.barometric_pressure import Barometric_pressure_operator
        domain = _make_domain(stage=1.0)   # depth h = 1.0

        def pressure(t, x, y):
            return x   # gradient (1, 0) over every triangle

        P = Barometric_pressure_operator(domain, pressure=pressure)
        P()

        expected = 1.0 * 1.0 / rho_w   # dt * h * px / rho_w
        np.testing.assert_allclose(
            domain.quantities['xmomentum'].centroid_values,
            expected, rtol=1e-10)
        np.testing.assert_allclose(
            domain.quantities['ymomentum'].centroid_values,
            0.0, atol=1e-14)

    def test_use_coordinates_false_needs_file_function(self):
        """use_coordinates=False requires a single-quantity file_function."""
        from anuga.operators.barometric_pressure import Barometric_pressure_operator
        domain = _make_domain()
        with pytest.raises(ValueError):
            Barometric_pressure_operator(domain, use_coordinates=False)
        with pytest.raises(ValueError):
            Barometric_pressure_operator(domain, pressure=lambda t, x, y: x,
                                         use_coordinates=False)

    def test_parallel_safe(self):
        from anuga.operators.barometric_pressure import Barometric_pressure_operator
        domain = _make_domain()
        assert Barometric_pressure_operator(domain).parallel_safe() is True


# ---------------------------------------------------------------------------
# Wind_stress_operator
# ---------------------------------------------------------------------------

class TestWindStressOperator:

    def test_east_wind_updates_xmom_only(self):
        """phi=0 (east): xmom increases, ymom stays zero."""
        from anuga.operators.wind_stress_operator import Wind_stress_operator
        domain = _make_domain()
        speed = 10.0
        W = Wind_stress_operator(domain, speed=speed, phi=0.0)
        W()

        const = eta_w * rho_a / rho_w
        # u=speed, v=0 → S=const*speed → xmom += dt*S*speed
        expected_xmom = 1.0 * (const * speed) * speed
        np.testing.assert_allclose(
            domain.quantities['xmomentum'].centroid_values,
            expected_xmom, rtol=1e-10)
        np.testing.assert_allclose(
            domain.quantities['ymomentum'].centroid_values,
            0.0, atol=1e-12)

    def test_north_wind_updates_ymom_only(self):
        """phi=90 (north): ymom increases, xmom stays zero."""
        from anuga.operators.wind_stress_operator import Wind_stress_operator
        domain = _make_domain()
        speed = 10.0
        W = Wind_stress_operator(domain, speed=speed, phi=90.0)
        W()

        const = eta_w * rho_a / rho_w
        expected_ymom = 1.0 * (const * speed) * speed
        np.testing.assert_allclose(
            domain.quantities['xmomentum'].centroid_values,
            0.0, atol=1e-12)
        np.testing.assert_allclose(
            domain.quantities['ymomentum'].centroid_values,
            expected_ymom, rtol=1e-10)

    def test_callable_speed_phi_matches_constant(self):
        """Callable speed/phi produce the same update as scalar equivalents."""
        from anuga.operators.wind_stress_operator import Wind_stress_operator
        domain1 = _make_domain()
        domain2 = _make_domain()

        W1 = Wind_stress_operator(domain1, speed=5.0, phi=45.0)
        W1()

        W2 = Wind_stress_operator(
            domain2,
            speed=lambda t, x, y: np.full(len(x), 5.0),
            phi=lambda t, x, y: np.full(len(x), 45.0),
        )
        W2()

        np.testing.assert_allclose(
            domain1.quantities['xmomentum'].centroid_values,
            domain2.quantities['xmomentum'].centroid_values, rtol=1e-10)
        np.testing.assert_allclose(
            domain1.quantities['ymomentum'].centroid_values,
            domain2.quantities['ymomentum'].centroid_values, rtol=1e-10)

    def test_zero_speed_no_update(self):
        """Zero wind speed: no momentum change."""
        from anuga.operators.wind_stress_operator import Wind_stress_operator
        domain = _make_domain()
        W = Wind_stress_operator(domain, speed=0.0, phi=0.0)
        W()
        np.testing.assert_allclose(
            domain.quantities['xmomentum'].centroid_values, 0.0, atol=1e-14)
        np.testing.assert_allclose(
            domain.quantities['ymomentum'].centroid_values, 0.0, atol=1e-14)

    def test_parallel_safe(self):
        from anuga.operators.wind_stress_operator import Wind_stress_operator
        domain = _make_domain()
        assert Wind_stress_operator(domain).parallel_safe() is True
