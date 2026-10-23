"""Exner bed evolution [G-4]: LM15 Example 2, deposition in quiescent water.

Liu, Mohammadian, Kurganov & Infante Sedano (2015), JCP 300, 202-218, Example
2, which they describe as "designed to verify that the added source terms will
not affect the well-balanced property in quiescent water with uniform sediment
deposition".

That is the question a fixed-bed well-balanced scheme does not answer: it does
not stay well-balanced once z is time-varying. Here the bed AGGRADES under a
depositing plume while the free surface must stay flat and the velocity exactly
zero.

LM15's setup: domain [0,2]x[0,1], eta = 1, u = v = Z = 0,
c(x,y,0) = 0.7 exp(-5(x-0.9)^2 - 50(y-0.5)^2), porosity 0.28, d = 0.01,
rho_s = 2400, rho_f = 1000.

Two deliberate deviations. Reflective rather than zero-order-extrapolation
boundaries, since a lake at rest is preserved exactly by reflection and
boundary noise cannot then be mistaken for a well-balancedness failure -- which
makes the test stricter. And c_max/c_pack are raised for this test only: LM15's
peak c = 0.7 exceeds both defaults ([L-2] 0.30, [L-4] 0.65) and clipping it
would change the experiment. This test is about well-balancedness, not the
limiters.

The mode 1 / mode 2 comparison is in test_sediment_gpu.py.
"""
import numpy as np
import pytest

from anuga import Reflective_boundary, rectangular_cross_domain

LAMBDA = 0.28
D_G = 0.01
RHO_S, RHO_W = 2400.0, 1000.0
T_END = 20.0


def build():
    d = rectangular_cross_domain(42, 21, len1=2.0, len2=1.0)
    d.set_flow_algorithm('DE0')
    d.set_low_froude(0)
    d.store = False
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', 1.0)
    d.set_quantity('xmomentum', 0.0)
    d.set_quantity('ymomentum', 0.0)
    d.set_quantity('friction', 0.0)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    d.sediment_c_max = 0.80
    d.sediment_c_pack = 0.80
    d.sediment_porosity = LAMBDA
    x, y = d.centroid_coordinates[:, 0], d.centroid_coordinates[:, 1]
    c0 = 0.7 * np.exp(-5.0 * (x - 0.9) ** 2 - 50.0 * (y - 0.5) ** 2)
    d.add_sediment_fraction(name='sand', diameter=D_G, rho_s=RHO_S, tau_c_star=0.0,
                     initial_concentration=c0)
    return d, c0


@pytest.fixture(scope='module')
def settled():
    d, c0 = build()
    z0 = d.quantities['elevation'].centroid_values.copy()
    w0 = d.quantities['stage'].centroid_values.copy()
    m0 = float((d.tracer_conserved_values[0] * d.areas).sum())
    d.evolve_to_end(finaltime=T_END)
    return d, c0, z0, w0, m0


def test_the_free_surface_stays_flat(settled):
    d, _, _, w0, _ = settled
    w1 = d.quantities['stage'].centroid_values
    assert np.allclose(w1, w0, rtol=0, atol=1e-10)


def test_no_spurious_momentum_is_generated_by_the_moving_bed(settled):
    d = settled[0]
    u = d.quantities['xmomentum'].centroid_values
    v = d.quantities['ymomentum'].centroid_values
    assert np.abs(u).max() < 1e-10
    assert np.abs(v).max() < 1e-10


def test_the_bed_rises_where_the_plume_was(settled):
    d, _, z0, _, _ = settled
    dz = d.quantities['elevation'].centroid_values - z0
    assert dz.max() > 0.0
    com = float((d.centroid_coordinates[:, 0] * dz).sum() / max(dz.sum(), 1e-30))
    assert abs(com - 0.9) < 0.15, 'bed rise centred at x = %.3f, plume at 0.9' % com


def test_the_per_cell_budget_is_exact(settled):
    """The water is at rest, so there is no lateral transport and each cell
    keeps its own sediment: (1-lambda) dz = the mass that cell lost.

    Much stronger than "the bed barely moved where c was small" -- the Gaussian
    is never exactly zero, so that weaker form passes on any roughly-right
    answer.
    """
    d, c0, z0, w0, _ = settled
    dz = d.quantities['elevation'].centroid_values - z0
    m0_cell = c0 * (w0 - z0)
    m1_cell = d.tracer_conserved_values[0]
    assert np.allclose((1.0 - LAMBDA) * dz, m0_cell - m1_cell,
                       rtol=0, atol=1e-12)
    # Close to complete deposition, but not exactly: deposition decays
    # exponentially, so a little is still suspended at t = 20 s.
    dz_full = m0_cell / (1.0 - LAMBDA)
    assert np.abs(dz - dz_full).max() < 1e-4 * dz_full.max()


def test_suspended_loss_equals_bed_gain(settled):
    d, _, z0, _, m0 = settled
    dz = d.quantities['elevation'].centroid_values - z0
    m1 = float((d.tracer_conserved_values[0] * d.areas).sum())
    bed_gain = float(((1.0 - LAMBDA) * dz * d.areas).sum())
    assert abs(bed_gain - (m0 - m1)) <= 1e-8 * max(abs(m0 - m1), 1.0)


def test_depth_falls_by_exactly_the_bed_rise(settled):
    """Stage is held, so h = w - z falls by dz. This is the quiescent-water
    behaviour LM15 Example 2 requires: the pore space in the new bed is filled
    from the water column."""
    d, _, z0, w0, _ = settled
    z1 = d.quantities['elevation'].centroid_values
    h1 = np.maximum(d.quantities['stage'].centroid_values - z1, 0.0)
    assert np.allclose(h1, (w0 - z0) - (z1 - z0), rtol=0, atol=1e-12)
