"""The Rouse near-bed concentration ratio d*(Z) -- [S-4], PHYSICS_SPEC 4.3.

Deposition is driven by the concentration AT THE BED, but the transported
quantity is depth-averaged; d* = c_b/c bridges them. The kernel evaluates a
fitted polynomial rather than the quadrature, so the fit is checked against the
quadrature it stands in for -- and the COEFFICIENTS are read out of the kernel
source, so this tests what is compiled rather than a reimplementation of the
same idea.

The mode 1 / mode 2 comparison is in test_sediment_gpu.py.
"""
import os
import re

import numpy as np
import pytest
from scipy.integrate import quad

import anuga
from anuga import Reflective_boundary, rectangular_cross_domain

LEN = 500.0

_KERNEL_SRC = os.path.join(os.path.dirname(anuga.__file__),
                           'shallow_water', 'gpu', 'core_kernels.c')
# Only present in a source checkout; a wheel install has no .c files.
#
# This has to abort the module at COLLECTION time, not mark its tests skipped:
# the coefficients are read at import (below), so a pytestmark skipif still
# lets that import run and the whole session dies on a collection error --
# every test in the run, not just this file. An editable install never shows
# it, because there the .c file is right there in the source tree.
if not os.path.exists(_KERNEL_SRC):
    pytest.skip('kernel source not available (needed to read the fitted '
                'coefficients)', allow_module_level=True)


def _fit_from_kernel():
    """Read the fit's coefficients out of the kernel source.

    Returns None if the source is there but does not contain the fit -- a
    checkout whose kernel predates it, say. Anything raising here would raise
    during COLLECTION and take the whole session down with it, so the caller
    turns a None into a module-level skip instead.
    """
    try:
        src = open(_KERNEL_SRC).read()
        rows = re.findall(r'\{([-+0-9.e]+), ([-+0-9.e]+), ([-+0-9.e]+), ([-+0-9.e]+)\},',
                          src)
        coeffs = np.array([[float(v) for v in r] for r in rows[:7]])
        limits = {k: float(re.search(r'#define ANUGA_ROUSE_%s\s+([0-9.e+-]+)' % k,
                                     src).group(1))
                  for k in ('Z_LO', 'Z_HI', 'AH_LO', 'AH_HI')}
    except (OSError, AttributeError, ValueError, IndexError):
        return None
    if coeffs.shape != (7, 4):
        return None
    return coeffs, limits


_FIT = _fit_from_kernel()
if _FIT is None:
    pytest.skip('the Rouse fit could not be read from the kernel source',
                allow_module_level=True)
COEFFS, LIMITS = _FIT


def dstar_fit(Z, a_h):
    """The kernel's fit, evaluated here from its own coefficients."""
    Z = min(max(Z, LIMITS['Z_LO']), LIMITS['Z_HI'])
    a_h = min(max(a_h, LIMITS['AH_LO']), LIMITS['AH_HI'])
    L = np.log(a_h)
    P = 0.0
    for i in range(6, -1, -1):
        P = P * Z + (COEFFS[i][0] + COEFFS[i][1] * L
                     + COEFFS[i][2] * L ** 2 + COEFFS[i][3] * L ** 3)
    return max(np.exp(-Z * L + P), 1.0)


def dstar_quadrature(Z, a_h, z0_h=1e-4, h=1.0):
    """[S-4] by quadrature, with the corrected (h-z) Rouse factor.

    The published DL09 print has (z-a) here; that is a typo in the paper, not
    in the spec.
    """
    a, z0 = a_h * h, z0_h * h
    num = quad(lambda z: np.log(z / z0), a, h, limit=200)[0]
    den = quad(lambda z: (((h - z) / (h - a)) * (a / z)) ** Z * np.log(z / z0),
               a, h, limit=200)[0]
    return num / den


def test_the_fit_reproduces_the_quadrature_over_its_range():
    errs = [abs(dstar_fit(Z, ah) / dstar_quadrature(Z, ah) - 1.0)
            for Z in np.geomspace(LIMITS['Z_LO'], LIMITS['Z_HI'], 17)
            for ah in np.geomspace(LIMITS['AH_LO'], LIMITS['AH_HI'], 9)]
    assert max(errs) < 0.015, 'max relative error %.3f%%' % (100 * max(errs))


def test_the_fit_holds_in_the_extended_low_ah_band():
    """The range was extended down to a/h = 1e-3 to reach anugaSed's regime."""
    errs = [abs(dstar_fit(Z, ah) / dstar_quadrature(Z, ah) - 1.0)
            for Z in np.geomspace(LIMITS['Z_LO'], LIMITS['Z_HI'], 13)
            for ah in np.geomspace(1e-3, 0.01, 5)]
    assert max(errs) < 0.015
    assert LIMITS['AH_LO'] <= 1e-3 + 1e-12


def test_d_star_is_never_below_one():
    """DL09: the near-bed concentration always exceeds the depth average."""
    assert all(dstar_fit(Z, ah) >= 1.0
               for Z in np.geomspace(0.001, 10, 40)
               for ah in (0.01, 0.05, 0.15))


def test_d_star_increases_with_the_rouse_number():
    assert all(dstar_fit(z2, 0.05) > dstar_fit(z1, 0.05)
               for z1, z2 in zip([0.02, 0.1, 0.5, 1.0], [0.1, 0.5, 1.0, 2.0]))


def test_d_star_matches_DL09_figure_4_for_small_Z():
    assert 1.0 <= dstar_fit(0.09, 0.05) <= 3.0


def test_the_well_mixed_limit_is_recovered():
    assert abs(dstar_fit(LIMITS['Z_LO'], 0.05) - 1.0) < 0.05


def test_out_of_range_inputs_are_clamped_not_extrapolated():
    """An 8th-degree polynomial extrapolates catastrophically; the fit is
    clamped at its range edges instead."""
    assert dstar_fit(50.0, 0.05) == dstar_fit(LIMITS['Z_HI'], 0.05)
    assert dstar_fit(1.0, 1e-6) == dstar_fit(1.0, LIMITS['AH_LO'])
    assert np.isfinite(dstar_fit(1e6, 1e-9))
    assert dstar_fit(1e6, 1e-9) >= 1.0


def _lake(d_star_mode=0, dt=1.0):
    d = rectangular_cross_domain(10, 10, len1=LEN, len2=LEN)
    d.set_flow_algorithm('DE0')
    d.set_low_froude(0)
    d.store = False
    d.set_quantity('elevation', lambda x, y: -0.01 * x)
    d.set_quantity('stage', lambda x, y: -0.01 * x + 1.0)
    d.set_quantity('friction', 0.03)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    d.evolve_max_timestep = dt
    d.sediment_d_star_mode = d_star_mode
    return d


def test_constant_mode_keeps_the_well_mixed_semantics():
    d = _lake(d_star_mode=0)
    d.add_sediment_fraction(name='s', diameter=1e-4, d_star=1.0, initial_concentration=0.02)
    d.evolve_to_end(finaltime=30.0)
    assert np.isfinite(d.get_tracer('s')).all()


def test_rouse_deposits_at_least_as_fast_as_the_well_mixed_limit():
    """d* >= 1 concentrates sediment near the bed, so deposition cannot be
    slower than the d* = 1 case. Direction, not magnitude -- the magnitude is
    what the fit tests above cover."""
    const = _lake(d_star_mode=0)
    const.add_sediment_fraction(name='s', diameter=1e-4, d_star=1.0, initial_concentration=0.02)
    const.evolve_to_end(finaltime=30.0)

    rouse = _lake(d_star_mode=1)
    rouse.add_sediment_fraction(name='s', diameter=1e-4, initial_concentration=0.02)
    rouse.evolve_to_end(finaltime=30.0)

    assert (rouse.get_tracer('s').mean()
            <= const.get_tracer('s').mean() * (1.0 + 1e-9))
