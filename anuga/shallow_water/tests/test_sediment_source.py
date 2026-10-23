"""The suspended sediment source term: settling, deposition, entrainment, limiters.

The exchange E - D of [G-3], applied as a fractional step. Covered here:

  settling      [S-1] Ferguson & Church, against the spec's verified value
  deposition    [D-1], against the analytic m0 exp(-v_s t/h) on a fixed bed
  entrainment   [E-1] Shields, including that it is off below threshold
  limiters      [L-1] positivity, [L-2] the concentration ceiling, [L-4] packing

The mode 1 / mode 2 comparisons are in test_sediment_gpu.py.
"""
import numpy as np
import warnings

import pytest

from anuga import Reflective_boundary, rectangular_cross_domain
from anuga import Sediment_transport_operator

LEN = 500.0
DEPTH, DIAM = 1.0, 1.0e-4


def still(nxy=10, depth=1.0, dt=None):
    """A flat lake at rest: no advection, so only the source term acts."""
    d = rectangular_cross_domain(nxy, nxy, len1=LEN, len2=LEN)
    d.set_flow_algorithm('DE0')
    d.set_low_froude(0)
    d.store = False
    d.set_quantity('elevation', 0.0)
    d.set_quantity('stage', depth)
    d.set_quantity('xmomentum', 0.0)
    d.set_quantity('ymomentum', 0.0)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    if dt is not None:
        d.evolve_max_timestep = dt
    return d


def channel(nxy=10, depth=1.0, slope=0.01, n_manning=0.03, dt=1.0):
    d = rectangular_cross_domain(nxy, nxy, len1=LEN, len2=LEN)
    d.set_flow_algorithm('DE0')
    d.set_low_froude(0)
    d.store = False
    d.set_quantity('elevation', lambda x, y: -slope * x)
    d.set_quantity('stage', lambda x, y: -slope * x + depth)
    d.set_quantity('friction', n_manning)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    d.evolve_max_timestep = dt
    return d


def tilted(d_star_mode=1, floor=0.01):
    d = rectangular_cross_domain(10, 10, len1=LEN, len2=LEN)
    d.set_flow_algorithm('DE0')
    d.set_low_froude(0)
    d.store = False
    d.set_quantity('elevation', lambda x, y: -0.01 * x)
    d.set_quantity('stage', lambda x, y: -0.01 * x + 1.0)
    d.set_quantity('friction', 0.03)
    d.set_boundary({t: Reflective_boundary(d) for t in d.get_boundary_tags()})
    d.evolve_max_timestep = 1.0
    d.sediment_d_star_mode = d_star_mode
    d.sediment_a_h_floor = floor
    return d


# ---------------------------------------------------------------- settling

def test_settling_velocity_matches_the_published_value():
    """[S-1] Ferguson & Church for 0.045 mm quartz; P13 give 1.75e-3 m/s."""
    v_s = still().settling_velocity(4.5e-5)
    assert abs(v_s - 1.75e-3) / 1.75e-3 < 0.02, 'v_s = %.6e m/s' % v_s


def test_settling_velocity_increases_with_grain_size():
    d = still()
    assert (d.settling_velocity(1e-3) > d.settling_velocity(1e-4)
            > d.settling_velocity(1e-5))


def test_natural_grain_constants_fall_slower_than_smooth_spheres():
    d = still()
    assert d.settling_velocity(1e-4, C1=1.0, C2=1.1) != d.settling_velocity(1e-4)


# ---------------------------------------------------------------- deposition

def test_deposition_follows_the_analytic_decay():
    """On a FIXED bed and still water, [D-1] reduces to dm/dt = -v_s m / h,
    whose solution is m0 exp(-v_s t/h). The bed is held so that h stays
    constant, which is what the analytic solution assumes."""
    d = still(depth=DEPTH, dt=1.0)
    d.sediment_bed_evolution = False
    d.add_sediment_fraction(name='sand', diameter=DIAM, initial_concentration=0.05)
    v_s = d.sediment_settling_velocity[0]
    d.evolve_to_end(finaltime=60.0)
    exact = 0.05 * np.exp(-v_s * 60.0 / DEPTH)
    assert abs(d.get_tracer('sand').mean() - exact) / exact < 0.02


def test_the_source_integration_is_first_order_in_dt():
    errs = {}
    for dt in (4.0, 1.0):
        d = still(depth=DEPTH, dt=dt)
        d.sediment_bed_evolution = False
        d.add_sediment_fraction(name='sand', diameter=DIAM, initial_concentration=0.05)
        v_s = d.sediment_settling_velocity[0]
        d.evolve_to_end(finaltime=60.0)
        exact = 0.05 * np.exp(-v_s * 60.0 / DEPTH)
        errs[dt] = abs(d.get_tracer('sand').mean() - exact) / exact
    assert errs[1.0] < errs[4.0], errs


def test_d_star_zero_disables_deposition_entirely():
    """The sediment class then behaves as a plain tracer, bit for bit."""
    base = still(depth=DEPTH, dt=1.0)
    base.add_tracer('plain', beta=1.0)
    base.set_tracer('plain', 0.05)
    base.evolve_to_end(finaltime=20.0)

    zero = still(depth=DEPTH, dt=1.0)
    zero.add_sediment_fraction(name='zero', diameter=1.0e-4, d_star=0.0,
                        initial_concentration=0.05)
    zero.evolve_to_end(finaltime=20.0)

    assert np.allclose(zero.get_tracer('zero'), base.get_tracer('plain'),
                       rtol=0, atol=1e-15)


# ---------------------------------------------------------------- limiters

def test_an_aggressive_settler_never_drives_mass_negative():
    """[L-1]: the source may remove at most the sediment PRESENT. A 5 mm grain
    settles many times the depth per step, so an unguarded scheme goes
    negative -- and a negative m flips deposition's sign and starts creating
    sediment."""
    d = still(depth=DEPTH, dt=1.0)
    d.add_sediment_fraction(name='fast', diameter=5.0e-3, initial_concentration=0.05)
    d.evolve_to_end(finaltime=30.0)
    assert d.tracer_conserved_values[0].min() >= 0.0
    assert d.get_tracer('fast').max() < 1e-6, 'it should deposit essentially all'


def test_concentration_stays_under_c_max():
    d = still(depth=DEPTH, dt=1.0)
    d.sediment_c_max = 0.10
    d.add_sediment_fraction(name='capped', diameter=1.0e-4, initial_concentration=0.05)
    d.evolve_to_end(finaltime=10.0)
    assert d.get_tracer('capped').max() <= 0.10 + 1e-12


# ---------------------------------------------------------------- API

def test_a_second_grain_size_registers_at_the_next_index():
    """The constructor returns the OPERATOR now, so check the order directly."""
    d = still()
    d.add_sediment_fraction(name='a', diameter=1e-4)
    d.add_sediment_fraction(name='b', diameter=5e-4)

    assert d.get_sediment_names() == ['a', 'b'], 'registration order lost'
    assert d.get_tracer_index('a') == 0
    assert d.get_tracer_index('b') == 1
    assert (d.sediment_settling_velocity[1]
            > d.sediment_settling_velocity[0]), 'per-grain-size v_s must differ'


def test_however_many_grain_sizes_there_is_one_operator():
    """One fractional step, however many grain sizes.

    Two operators in the list would apply the bed exchange twice per timestep.
    """
    d = still()
    d.add_sediment_fraction(name='a', diameter=1e-4)
    d.add_sediment_fraction(name='b', diameter=5e-4)

    n = sum(isinstance(op, Sediment_transport_operator)
            for op in d.fractional_step_operators)
    assert n == 1, 'the domain has %d sediment operators, expected 1' % n
    assert d.n_sediment_classes == 2


def test_initialize_returns_the_one_operator_however_often_it_is_called():
    d = still()
    first = d.initialize_sediment_operator()
    d.add_sediment_fraction(name='a', diameter=1e-4)
    second = d.initialize_sediment_operator(porosity=0.28)

    assert first is second, 'a second call built a separate operator'
    assert d.sediment_porosity == 0.28, 'the second call was ignored'
    n = sum(isinstance(op, Sediment_transport_operator)
            for op in d.fractional_step_operators)
    assert n == 1


def test_the_two_calls_do_not_share_a_parameter():
    """The split is what stops a parameter being silently dropped.

    Under the old single entry point, `diameter=` without `name=` added
    nothing and raised nothing, and a domain-wide `rho_w=` on a second call
    was discarded. Neither can be expressed now: a per-grain-size parameter
    is not in initialize_sediment_operator's signature, and a domain-wide one
    is not in add_sediment_fraction's, so both are a TypeError rather than a silent
    no-op.
    """
    d = still()
    with pytest.raises(TypeError):
        d.initialize_sediment_operator(diameter=5.0e-5)
    # A domain-wide name would otherwise be swallowed by **settling_kwargs and
    # resurface as a TypeError naming settling_velocity, which tells the caller
    # nothing about what they did wrong.
    with pytest.raises(TypeError, match='property of the run'):
        d.add_sediment_fraction(name='sand', diameter=2.0e-4, rho_w=1025.0)
    with pytest.raises(TypeError, match='properties of the run'):
        d.add_sediment_fraction(name='sand', diameter=2.0e-4, rho_w=1025.0, c_max=0.2)
    assert d.get_sediment_names() == [], 'nothing should have been added'


def test_a_grain_size_needs_a_diameter():
    d = still()
    with pytest.raises(TypeError):
        d.add_sediment_fraction(name='sand')


def test_initialize_alone_registers_no_grain_size():
    """Domain-wide setup is legal before any grain size exists."""
    d = still()
    op = d.initialize_sediment_operator(porosity=0.28)
    assert d.get_sediment_names() == []
    assert op in d.fractional_step_operators
    assert d.sediment_porosity == 0.28


def test_add_sediment_fraction_creates_the_operator_if_initialize_was_skipped():
    d = still()
    d.add_sediment_fraction(name='sand', diameter=2.0e-4)
    ops = [o for o in d.fractional_step_operators
           if isinstance(o, Sediment_transport_operator)]
    assert len(ops) == 1
    assert d.get_sediment_names() == ['sand']


def test_the_two_calls_may_come_in_either_order():
    """initialize after add still applies, and R follows the new rho_w."""
    d = still()
    d.add_sediment_fraction(name='sand', diameter=2.0e-4, rho_s=2650.0)
    assert abs(d.sediment_R[0] - 1.65) < 1e-12

    d.initialize_sediment_operator(rho_w=1250.0)
    assert abs(d.sediment_R[0] - (2650.0 / 1250.0 - 1.0)) < 1e-12, (
        'R kept the density it was registered with')


def test_set_bed_material_does_not_reset_the_water_density():
    """It used to take rho_w=1000.0 as a DEFAULT, so every call clobbered it.

    Choosing an erosion law says nothing about the density of the water, and
    the reset was silent: R and v_s quietly reverted to their fresh-water
    values partway through a setup.
    """
    d = still()
    d.initialize_sediment_operator(rho_w=1025.0)
    d.add_sediment_fraction(name='sand', diameter=2.0e-4, rho_s=2650.0)
    R = d.sediment_R[0]

    d.set_bed_material('cohesive', tau_crit=0.5)

    assert d.sediment_rho_w == 1025.0
    assert d.sediment_R[0] == R


def test_evolving_with_no_grain_size_warns_rather_than_doing_nothing():
    """A bare operator is legal -- it is how operator order is controlled --
    but evolving that way transports nothing.

    Bedload does not rescue it: core_apply_bedload returns immediately when
    n_sediment_classes is 0, because the diameter and R that set the Shields
    stress live on a grain size. Without the warning the run completes with the
    bed untouched and nothing said.
    """
    d = still()
    d.initialize_sediment_operator()
    d.set_bedload('wong_parker_eq24')
    z0 = d.quantities['elevation'].centroid_values.copy()

    with pytest.warns(UserWarning, match='no grain size is registered'):
        d.evolve_to_end(finaltime=2.0)

    assert np.abs(d.quantities['elevation'].centroid_values - z0).max() == 0.0


def test_the_no_grain_size_warning_is_raised_once_not_per_timestep():
    d = still()
    d.initialize_sediment_operator()
    with pytest.warns(UserWarning) as record:
        d.evolve_to_end(finaltime=5.0)
    n = sum('no grain size is registered' in str(w.message) for w in record)
    assert n == 1, 'warned %d times; it must not fire every timestep' % n


def test_no_warning_once_a_grain_size_is_registered_after_the_operator():
    """The operator-ordering pattern must stay silent."""
    d = still()
    d.initialize_sediment_operator()
    d.add_sediment_fraction('sand', diameter=2.0e-4)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        d.evolve_to_end(finaltime=2.0)
    assert not [x for x in w if 'no grain size' in str(x.message)]


def test_a_non_positive_diameter_is_rejected():
    with pytest.raises(ValueError):
        still().add_sediment_fraction(name='c', diameter=0.0)


def test_mixing_add_tracer_and_add_sediment_fraction_is_rejected():
    """Grain size s must occupy tracer slot s, so interleaving is refused
    rather than silently breaking that correspondence."""
    d = still()
    d.add_tracer('plain')
    with pytest.raises(ValueError):
        d.add_sediment_fraction(name='s', diameter=1e-4)


# ---------------------------------------------------------------- entrainment

def test_no_entrainment_below_the_critical_shields_stress():
    d = still(depth=DEPTH, dt=1.0)
    d.add_sediment_fraction(name='sand', diameter=DIAM, tau_c_star=0.04,
                     initial_concentration=0.0)
    d.evolve_to_end(finaltime=20.0)
    assert float(np.abs(d.get_tracer('sand')).max()) == 0.0


def test_a_flowing_channel_entrains_from_a_clean_bed():
    d = channel()
    d.add_sediment_fraction(name='sand', diameter=DIAM, tau_c_star=0.04,
                     initial_concentration=0.0)
    d.evolve_to_end(finaltime=60.0)
    assert d.get_tracer('sand').max() > 0.0


def test_tau_c_star_zero_disables_entrainment():
    d = channel()
    d.add_sediment_fraction(name='sand', diameter=DIAM, tau_c_star=0.0,
                     initial_concentration=0.0)
    d.evolve_to_end(finaltime=60.0)
    assert float(np.abs(d.get_tracer('sand')).max()) == 0.0


@pytest.mark.slow
def test_violent_flow_stays_bounded_near_c_max():
    """[L-2] caps the SOURCE, not the state, so advective inflow can carry a
    cell slightly over c_max. The tolerance admits that; what it must not do is
    run away, which the longer run checks."""
    d = channel(depth=3.0, slope=0.05, n_manning=0.05)
    d.add_sediment_fraction(name='sand', diameter=DIAM, tau_c_star=0.04,
                     initial_concentration=0.0)
    d.evolve_to_end(finaltime=60.0)
    c = d.get_tracer('sand')
    assert np.isfinite(c).all()
    assert c.max() / d.sediment_c_max - 1.0 <= 0.02

    d.evolve_to_end(finaltime=300.0)
    c_long = d.get_tracer('sand')
    assert np.isfinite(c_long).all()
    assert c_long.max() / d.sediment_c_max - 1.0 <= 0.02


# ---------------------------------------------------------------- packing [L-4]

def test_deposition_only_never_creates_mass():
    d = tilted()
    d.add_sediment_fraction(name='s', diameter=1e-4, tau_c_star=0.0, initial_concentration=0.02)
    m0 = float((d.tracer_conserved_values[0] * d.areas).sum())
    d.evolve_to_end(finaltime=15.0)
    m1 = float((d.tracer_conserved_values[0] * d.areas).sum())
    assert m1 <= m0 * (1 + 1e-9), 'deposition may only remove'
    assert m1 >= -0.01 * m0, 'and must not overshoot into negative mass'


def test_a_near_still_start_does_not_deposit_everything_at_once():
    """[L-4]. Equilibrium Rouse d* at vanishing shear is unbounded, and without
    the packing bound this removed the ENTIRE suspended load in under a second
    against a physical timescale h/v_s of about 125 s. The packing-limited
    prediction is about 26%."""
    d = tilted()
    d.add_sediment_fraction(name='s', diameter=1e-4, tau_c_star=0.0, initial_concentration=0.02)
    m0 = float((d.tracer_conserved_values[0] * d.areas).sum())
    d.evolve_to_end(finaltime=1.0)
    m1 = float((d.tracer_conserved_values[0] * d.areas).sum())
    assert 1.0 - m1 / m0 < 0.6


def test_c_pack_is_exposed_and_binds():
    rates = {}
    for c_pack in (0.65, 0.05):
        d = tilted()
        d.sediment_c_pack = c_pack
        d.add_sediment_fraction(name='s', diameter=1e-4, tau_c_star=0.0,
                         initial_concentration=0.02)
        m0 = float((d.tracer_conserved_values[0] * d.areas).sum())
        d.evolve_to_end(finaltime=1.0)
        m1 = float((d.tracer_conserved_values[0] * d.areas).sum())
        rates[c_pack] = 1.0 - m1 / m0
    assert rates[0.05] < rates[0.65], rates
