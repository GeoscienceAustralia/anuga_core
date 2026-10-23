"""Morphological acceleration factor M (`morphological_factor`).

Every step's bed change from the suspended exchange [G-4] and from bedload
[G-5] is multiplied by M; the erodible base [L-5] is respected; the water
column is untouched. Checked in a still lake (deposition only, no advection)
and on a sloping channel with bedload.
"""
import numpy as np
import pytest

import anuga
from anuga import Dirichlet_boundary, Reflective_boundary, rectangular_cross_domain

LEN = 100.0


def still_lake(M=None, mode='legacy', c0=0.02, base_depth=None):
    d = rectangular_cross_domain(6, 6, len1=LEN, len2=LEN)
    d.set_flow_algorithm('DE0')
    d.set_low_froude(0)
    d.set_compute_mode(mode)
    d.store = False
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', 1.0)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    kw = {} if M is None else dict(morphological_factor=M)
    d.initialize_sediment_operator(porosity=0.3, bed_evolution=True, **kw)
    d.add_sediment_fraction('sand', diameter=2.0e-4, tau_c_star=0.0,
                            initial_concentration=c0)
    if base_depth is not None:
        d.set_erodible_base(depth=base_depth)
    d.evolve_max_timestep = 0.5
    return d


def one_step(d, dt=0.5):
    """Bed change and change of the conserved sediment mass h c over one
    step. The mass, not the concentration: the bed rises M times more, so
    the depth and with it c = m/h differ between runs while m does not."""
    z0 = d.quantities['elevation'].centroid_values.copy()
    m0 = d.tracer_conserved_values[0].copy()
    for _ in d.evolve(yieldstep=dt, finaltime=dt):
        pass
    return (d.quantities['elevation'].centroid_values - z0,
            d.tracer_conserved_values[0] - m0)


@pytest.mark.parametrize('mode', ['legacy', 'unified'])
def test_the_bed_change_scales_by_M_and_the_water_column_does_not(mode):
    dz1, dc1 = one_step(still_lake(mode=mode))
    dz5, dc5 = one_step(still_lake(M=5.0, mode=mode))
    k = 20
    assert dz1[k] > 0.0                      # deposition raises the bed
    assert dz5[k] == pytest.approx(5.0 * dz1[k], rel=1e-12)
    assert dc5[k] == pytest.approx(dc1[k], rel=1e-12)


def test_M_of_one_is_the_default_and_changes_nothing():
    d = still_lake()
    assert d.sediment_morphological_factor == 1.0
    dz_default, _ = one_step(still_lake())
    dz_one, _ = one_step(still_lake(M=1.0))
    assert np.array_equal(dz_default, dz_one)


def test_the_erodible_base_is_respected_under_M():
    """Entrainment from a thin erodible layer: with M large the bed would be
    dug far below the base in one step unless the cap is scaled."""
    d = rectangular_cross_domain(10, 4, len1=LEN, len2=40.0)
    d.set_flow_algorithm('DE1')
    d.store = False
    d.set_quantity('elevation', 0.0)
    d.set_quantity('friction', 0.0)
    d.set_quantity('stage', 1.0)
    d.set_quantity('xmomentum', 1.5)
    Bd = Dirichlet_boundary([1.0, 1.5, 0.0])
    Br = Reflective_boundary(d)
    d.set_boundary({'left': Bd, 'right': Bd, 'top': Br, 'bottom': Br})
    d.initialize_sediment_operator(porosity=0.3, bed_evolution=True,
                                   morphological_factor=50.0)
    d.set_sediment_friction('larsen_lamb', k_s=0.05)
    d.add_sediment_fraction('sand', diameter=2.0e-4, initial_concentration=0.0)
    d.set_tracer_boundary('sand', 'left', 0.0)
    base = 0.002
    d.set_erodible_base(depth=base)
    for _ in d.evolve(yieldstep=5.0, finaltime=20.0):
        pass
    z = d.quantities['elevation'].centroid_values
    assert z.min() >= -base - 1e-12, z.min()
    # and it was actually eroded to the base somewhere
    assert z.min() < -0.5 * base


def test_bedload_bed_change_scales_by_M():
    def run(M):
        d = rectangular_cross_domain(20, 2, len1=LEN, len2=10.0)
        d.set_flow_algorithm('DE1')
        d.store = False
        # a hump on the bed so that the bedload divergence is non-zero
        d.set_quantity('elevation', lambda x, y: 0.05 * np.exp(-((x - 50.0) / 8.0) ** 2))
        d.set_quantity('friction', 0.0)
        d.set_quantity('stage', 1.0)
        d.set_quantity('xmomentum', 1.5)
        Bd = Dirichlet_boundary([1.0, 1.5, 0.0])
        Br = Reflective_boundary(d)
        d.set_boundary({'left': Bd, 'right': Bd, 'top': Br, 'bottom': Br})
        d.initialize_sediment_operator(porosity=0.3, bed_evolution=True,
                                       morphological_factor=M)
        d.set_sediment_friction('larsen_lamb', k_s=0.05)
        d.set_deposition(tau_d=0.0, law='threshold')      # no suspended exchange
        d.add_sediment_fraction('sand', diameter=1.0e-3, tau_c_star=0.0,
                                initial_concentration=0.0)
        d.set_bedload('wong_parker_eq24', open_boundaries=['left', 'right'])
        z0 = d.quantities['elevation'].centroid_values.copy()
        # one step, below the CFL limit (~0.02 s here), so that the bed
        # change is a single application of the kernel
        d.evolve_max_timestep = 0.01
        for _ in d.evolve(yieldstep=0.01, finaltime=0.01):
            pass
        return d.quantities['elevation'].centroid_values - z0
    dz1, dz4 = run(1.0), run(4.0)
    assert np.abs(dz1).max() > 0.0
    # rk2: the second stage sees a state the first has moved by O(dt), so
    # the ratio is 4 to about 1e-9, not to roundoff
    assert np.allclose(dz4, 4.0 * dz1, rtol=1e-6, atol=1e-18)


def test_bad_values():
    d = still_lake()
    with pytest.raises(ValueError):
        d.set_sediment_parameters(morphological_factor=0.0)
    with pytest.raises(ValueError):
        d.set_sediment_parameters(morphological_factor=-3.0)
    d.set_sediment_parameters(morphological_factor=10.0)
    assert 'morphological M    : 10' in d.sediment_summary()
