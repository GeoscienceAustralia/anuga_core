"""Tests for LOW_FROUDE_TANGENTIAL_UPWIND: the edge-tangential momentum is
upwinded with the mass flux instead of being diffused at the gravity-wave
speed."""

import numpy as num
import pytest

import anuga
from anuga.config import g, epsilon
from anuga.config import (LOW_FROUDE_OFF, LOW_FROUDE_1, LOW_FROUDE_2,
                          LOW_FROUDE_TANGENTIAL_UPWIND as TU)
from anuga.shallow_water.sw_domain_openmp_ext import flux_function_central as flux_function


def _flux(normal, ql, qr, low_froude):
    """Edge flux between two flat-bed states of depth ql[0], qr[0] (bed at 0)."""
    ql = num.array(ql, float)
    qr = num.array(qr, float)
    hl, hr = ql[0], qr[0]
    edgeflux = num.zeros(3)
    max_speed, pressure_flux = flux_function(num.array(normal, float), ql, qr, hl, hr, hl, hr,
                                             edgeflux, epsilon, 0.0, g, 0.0, hl, hr, low_froude)
    return edgeflux, max_speed, pressure_flux


def test_shear_without_normal_flow_exchanges_no_momentum():
    # Edge normal along x; the two sides differ only in the tangential (y)
    # velocity and nothing flows through the edge.
    ql, qr = [1.0, 0.0, 1.0], [1.0, 0.0, 0.0]
    for mode in (LOW_FROUDE_OFF, LOW_FROUDE_1, LOW_FROUDE_2):
        f0, s0, p0 = _flux([1.0, 0.0], ql, qr, mode)
        f1, s1, p1 = _flux([1.0, 0.0], ql, qr, mode | TU)
        assert abs(f0[2]) > 1.0e-3            # central upwind diffuses the shear
        assert num.allclose(f1, [0.0, 0.0, 0.0], atol=1e-14)
        assert num.allclose(f0[:2], f1[:2])
        assert s0 == s1 and p0 == p1


@pytest.mark.parametrize('mode', [LOW_FROUDE_OFF, LOW_FROUDE_1, LOW_FROUDE_2])
@pytest.mark.parametrize('normal', [(1.0, 0.0), (0.6, 0.8), (-0.8, 0.6)])
@pytest.mark.parametrize('sign', [1.0, -1.0])
def test_tangential_flux_is_mass_flux_times_upwind_velocity(mode, normal, sign):
    n = num.array(normal)
    t = num.array([-n[1], n[0]])
    hl, hr = 1.2, 0.9
    unl, unr = sign * 0.8, sign * 0.5       # normal velocities, flow in direction sign
    vtl, vtr = 0.7, -0.4                    # tangential velocities
    ql = [hl, *(hl * (unl * n + vtl * t))]
    qr = [hr, *(hr * (unr * n + vtr * t))]
    f0, _, _ = _flux(normal, ql, qr, mode)
    f1, _, _ = _flux(normal, ql, qr, mode | TU)

    # Mass and normal momentum fluxes are unchanged by the flag
    assert num.allclose(f1[0], f0[0])
    assert num.allclose(f1[1:] @ n, f0[1:] @ n)
    # Tangential momentum flux = mass flux * upwind tangential velocity
    assert f1[0] * sign > 0
    upwind = vtl if f1[0] > 0 else vtr
    assert num.allclose(f1[1:] @ t, f1[0] * upwind)


def test_domain_api():
    domain = anuga.rectangular_cross_domain(2, 2)
    assert not domain.get_tangential_upwind()
    domain.set_low_froude(LOW_FROUDE_2)
    domain.set_tangential_upwind()
    assert domain.get_tangential_upwind()
    assert domain.low_froude == LOW_FROUDE_2 | TU
    domain.set_tangential_upwind(False)
    assert domain.low_froude == LOW_FROUDE_2
    domain.set_low_froude(LOW_FROUDE_1 | TU)
    assert domain.get_tangential_upwind()
    with pytest.raises(AssertionError):
        domain.set_low_froude(3)


def _shear_layer(tangential_upwind, compute_mode):
    """Frictionless flat-bed shear layer: u = 1 m/s for y > 10, still water
    below. The exact solution keeps the step; return the velocity difference
    across it after 10 s, away from the ends."""
    domain = anuga.rectangular_cross_domain(60, 20, 60.0, 20.0)
    domain.set_compute_mode(compute_mode)
    domain.set_store(False)
    domain.set_quantity('elevation', 0.0)
    domain.set_quantity('friction', 0.0)
    domain.set_quantity('stage', 1.0)
    domain.set_quantity('xmomentum', lambda x, y: num.where(y > 10.0, 1.0, 0.0))
    domain.set_tangential_upwind(tangential_upwind)
    br = anuga.Reflective_boundary(domain)
    bt = anuga.Transmissive_boundary(domain)
    domain.set_boundary({'left': bt, 'right': bt, 'top': br, 'bottom': br})
    for _ in domain.evolve(yieldstep=10.0, finaltime=10.0):
        pass
    x, y = domain.centroid_coordinates.T
    h = domain.quantities['stage'].centroid_values - domain.quantities['elevation'].centroid_values
    u = domain.quantities['xmomentum'].centroid_values / h
    mid = (x > 20.0) & (x < 40.0)
    return u[mid & (y > 10.0) & (y < 12.0)].mean() - u[mid & (y > 8.0) & (y < 10.0)].mean()


@pytest.mark.parametrize('compute_mode', ['legacy', 'unified'])
def test_shear_layer_is_kept(compute_mode):
    jump_central = _shear_layer(False, compute_mode)
    jump_upwind = _shear_layer(True, compute_mode)
    # Central upwind smears the step across several cells; the contact-wave
    # upwinding keeps most of it.
    assert jump_central < 0.6
    assert jump_upwind > 0.75
