"""Tracer mass budget: mass(t) - mass(0) == boundary flux integral.

The tracer counterpart of the water balance. `get_tracer_mass` integrates the
conserved m = h*c; the flux kernel records what crosses each domain-boundary
edge, per substep; `tracer_flux_integral_operator` time-integrates it with the
timestepping method's own weights.
"""

import numpy as num
import pytest

import anuga


def _closed(n=12):
    d = anuga.rectangular_cross_domain(n, n)
    d.set_quantity('elevation', lambda x, y: -0.1 * x)
    d.set_quantity('stage', 0.5)
    b = anuga.Reflective_boundary(d)
    d.set_boundary({'left': b, 'right': b, 'top': b, 'bottom': b})
    return d


def _draining(n=16):
    d = anuga.rectangular_cross_domain(n, n)
    d.set_quantity('elevation', lambda x, y: -0.3 * x)
    d.set_quantity('stage', lambda x, y: num.maximum(0.5 - 0.3 * x, 0.05))
    Bt = anuga.Transmissive_boundary(d)
    Br = anuga.Reflective_boundary(d)
    d.set_boundary({'left': Br, 'right': Bt, 'top': Br, 'bottom': Br})
    return d


def test_mass_is_the_integral_of_h_times_c():
    d = _closed()
    d.add_tracer('c', initial_value=0.3)
    # h = 0.5 - (-0.1x) so the mass is sum(h*c*area); compare against the
    # conserved variable directly rather than restating the formula.
    expected = float(num.sum(d.tracer_conserved_values[0] * d.areas))
    assert num.isclose(d.get_tracer_mass('c'), expected, rtol=1e-12)


def test_a_closed_domain_conserves_tracer_mass():
    d = _closed()
    d.add_tracer('c', initial_value=0.3)
    m0 = d.get_tracer_mass('c')
    for _ in d.evolve(yieldstep=0.5, finaltime=2.0):
        pass
    change, flux, disc = d.check_tracer_conservation('c')
    assert abs(change) < 1e-12 * max(m0, 1.0), 'mass moved on a closed domain'
    assert abs(flux) < 1e-12 * max(m0, 1.0), 'flux across a closed boundary'
    assert abs(disc) < 1e-12 * max(m0, 1.0)


def test_the_budget_balances_when_tracer_leaves():
    """The real test: tracer leaves, and the flux integral accounts for it."""
    d = _draining()
    d.add_tracer('c', initial_value=0.4)
    m0 = d.get_tracer_mass('c')
    for _ in d.evolve(yieldstep=0.25, finaltime=3.0):
        pass
    change, flux, disc = d.check_tracer_conservation('c')

    assert change < -0.1 * m0, 'the test is vacuous unless tracer actually left'
    assert flux < 0.0, 'net flux should be outward'
    assert abs(disc) < 1e-10 * abs(change), \
        'mass change and boundary flux integral disagree: %g vs %g' % (change, flux)


def test_two_tracers_are_accounted_separately():
    d = _draining()
    d.add_tracer('a', initial_value=0.4)
    d.add_tracer('b', initial_value=0.1)
    for _ in d.evolve(yieldstep=0.5, finaltime=2.0):
        pass
    ca, fa, da = d.check_tracer_conservation('a')
    cb, fb, db = d.check_tracer_conservation('b')
    assert abs(da) < 1e-10 * abs(ca)
    assert abs(db) < 1e-10 * abs(cb)
    # b started at a quarter of a's concentration, so it should have lost
    # about a quarter as much -- not the same amount, which would mean the
    # rows are aliased.
    assert not num.isclose(ca, cb, rtol=1e-3)


def test_the_baseline_is_taken_after_the_initial_value():
    # _tracer_initial_mass must reflect the seeded field, not an empty domain.
    d = _closed()
    d.add_tracer('c', initial_value=0.25)
    change, _, _ = d.check_tracer_conservation('c')
    assert abs(change) < 1e-14, 'baseline was taken before the tracer was seeded'


def test_a_second_tracer_does_not_reset_the_first_accounting():
    d = _draining()
    d.add_tracer('a', initial_value=0.4)
    for _ in d.evolve(yieldstep=0.5, finaltime=1.0):
        pass
    integral_before = d.get_tracer_boundary_flux_integral('a')
    assert integral_before != 0.0
    d.add_tracer('b', initial_value=0.2)      # reallocates every tracer array
    assert d.get_tracer_boundary_flux_integral('a') == integral_before


def test_set_tracer_before_the_run_rebases_the_baseline():
    """add_tracer(name) then set_tracer(name, field) is the usual idiom."""
    d = _closed()
    d.add_tracer('c')                      # seeds zero
    x = d.centroid_coordinates[:, 0]
    d.set_tracer('c', num.where(x < 0.5, 1.0, 0.0))
    assert d.get_tracer_mass('c') > 0.0
    change, _, _ = d.check_tracer_conservation('c')
    assert abs(change) < 1e-14, 'baseline did not follow set_tracer'


def test_set_tracer_mid_run_is_not_silently_absorbed():
    """A mid-run set_tracer really does break the budget; it must show."""
    d = _closed()
    d.add_tracer('c', initial_value=0.2)
    for _ in d.evolve(yieldstep=0.5, finaltime=1.0):
        pass
    before = d.get_tracer_mass('c')
    d.set_tracer('c', 0.6)                 # a deliberate intervention
    added = d.get_tracer_mass('c') - before
    assert added > 0.0

    _, _, disc = d.check_tracer_conservation('c')
    assert abs(disc - added) < 1e-10 * added, \
        'mid-run set_tracer was rebased away instead of showing as a discrepancy'


def _fed(n=16):
    """Water enters on the left; the domain itself starts clean."""
    d = anuga.rectangular_cross_domain(n, n)
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', 1.0)
    d.set_quantity('xmomentum', 0.5)
    Bi = anuga.Dirichlet_boundary([1.0, 0.5, 0.0])
    Bo = anuga.Transmissive_boundary(d)
    Br = anuga.Reflective_boundary(d)
    d.set_boundary({'left': Bi, 'right': Bo, 'top': Br, 'bottom': Br})
    return d


def test_tracer_arriving_on_an_inflow_is_counted_as_it_enters():
    """The other sign of the budget: set_tracer_boundary feeding a clean domain.

    Pairs the two halves -- a prescribed inflow concentration is only applied
    where water flows in, and what it brings has to appear in the integral.
    """
    d = _fed()
    d.add_tracer('c', initial_value=0.0)
    d.set_tracer_boundary('c', 'left', 1.0)
    assert d.get_tracer_mass('c') == 0.0

    for _ in d.evolve(yieldstep=0.25, finaltime=1.0):
        pass
    change, flux, disc = d.check_tracer_conservation('c')

    assert change > 0.0, 'no tracer entered through the inflow'
    assert flux > 0.0, 'net flux should be inward'
    assert abs(disc) < 1e-10 * abs(change), \
        'inflow budget does not balance: change %g vs flux %g' % (change, flux)


def test_an_unset_inflow_boundary_brings_nothing():
    """The documented default: an unset boundary carries c = 0."""
    d = _fed()
    d.add_tracer('c', initial_value=0.0)      # deliberately no set_tracer_boundary
    for _ in d.evolve(yieldstep=0.25, finaltime=1.0):
        pass
    change, flux, _ = d.check_tracer_conservation('c')
    assert abs(change) < 1e-14, 'clean water carried tracer in'
    assert abs(flux) < 1e-14
