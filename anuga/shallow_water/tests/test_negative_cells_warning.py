"""Tests for the volume-based "possible loss of conservation" warning emitted
by Domain.update_conserved_quantities().

The warning fires when clamping negative-depth cells to zero depth adds a large
enough fraction of the total water volume. Warning on *volume* rather than cell
count means a handful of cells clamped by a near-zero depth (normal in
wetting/drying) stays silent, while a genuinely large loss still warns.
"""

import unittest
import warnings

import numpy as num
import pytest

import anuga


def make_domain(threshold=None):
    """Flat 1 m of water over a 0 m bed on a small mesh."""
    domain = anuga.rectangular_cross_domain(4, 4)
    # White-box tests: they call update_conserved_quantities() directly (outside
    # evolve, so no gpu_interface is set up) and inspect host centroid_values.
    # Pin to legacy (mode 1) so the CPU openmp kernel runs regardless of the
    # ANUGA_DEFAULT_COMPUTE_MODE the suite is exercised under.
    domain.set_compute_mode('legacy')
    domain.set_quantity('elevation', 0.0)
    domain.set_quantity('stage', 1.0)
    Br = anuga.Reflective_boundary(domain)
    domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
    domain.timestep = 0.0   # so update_conserved_quantities makes no flux change
    if threshold is not None:
        domain.set_negative_volume_warning_fraction(threshold)
    return domain


def _force_negative_and_update(domain, deficit, n_cells=3):
    """Set n_cells to a negative depth of `deficit`, then run the update, and
    return any warnings raised."""
    stage_c = domain.quantities['stage'].centroid_values
    bed_c = domain.quantities['elevation'].centroid_values
    stage_c[:n_cells] = bed_c[:n_cells] - deficit
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        domain.update_conserved_quantities()
    return w


class Test_negative_cells_warning(unittest.TestCase):

    def test_default_threshold_from_config(self):
        from anuga.config import negative_volume_warning_fraction
        domain = make_domain()
        self.assertEqual(domain.get_negative_volume_warning_fraction(),
                         negative_volume_warning_fraction)

    def test_setter_getter(self):
        domain = make_domain()
        domain.set_negative_volume_warning_fraction(0.02)
        self.assertEqual(domain.get_negative_volume_warning_fraction(), 0.02)

    def test_tiny_deficit_is_silent(self):
        """A near-zero depth clamp adds negligible volume: no warning by default."""
        domain = make_domain()
        w = _force_negative_and_update(domain, deficit=1.0e-7)
        self.assertEqual(len(w), 0)

    def test_large_deficit_warns(self):
        """A large volume added by clamping warns, and the message reports it."""
        domain = make_domain()
        w = _force_negative_and_update(domain, deficit=0.5)
        self.assertEqual(len(w), 1)
        msg = str(w[-1].message)
        self.assertIn('loss of conservation', msg)
        self.assertIn('m^3', msg)

    def test_threshold_zero_warns_on_any_added_volume(self):
        """threshold=0.0 restores warn-on-any-loss behaviour."""
        domain = make_domain(threshold=0.0)
        w = _force_negative_and_update(domain, deficit=1.0e-7)
        self.assertEqual(len(w), 1)

    def test_no_negative_cells_no_warning(self):
        """Nothing clamped -> no warning regardless of threshold."""
        domain = make_domain(threshold=0.0)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            domain.update_conserved_quantities()
        self.assertEqual(len(w), 0)

    def test_noise_floor_suppresses_femto_losses(self):
        """A nearly-dry domain can clamp deficits that are a large fraction of an
        essentially-zero total volume; the absolute noise floor suppresses those
        even with the relative threshold set to 0 (warn-on-any)."""
        domain = make_domain(threshold=0.0)   # relative test always passes
        stage_c = domain.quantities['stage'].centroid_values
        bed_c = domain.quantities['elevation'].centroid_values
        areas = domain.areas
        # Femto-scale water everywhere, a few cells pushed slightly negative.
        stage_c[:] = bed_c + 1.0e-12
        stage_c[:3] = bed_c[:3] - 1.0e-12
        added = float(num.sum(1.0e-12 * areas[:3]))
        self.assertLess(added, 1.0e-9)   # far below the noise floor
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            domain.update_conserved_quantities()
        self.assertEqual(len(w), 0)

    def test_kernel_returns_count_and_volume(self):
        """The openmp kernel reports both the clamped-cell count and the exact
        added volume, and still clamps the cells to the bed."""
        from anuga.shallow_water.sw_domain_openmp_ext import \
            update_conserved_quantities
        domain = make_domain()
        stage_c = domain.quantities['stage'].centroid_values
        bed_c = domain.quantities['elevation'].centroid_values
        areas = domain.areas
        stage_c[:3] = bed_c[:3] - 0.5
        expected_volume = float(num.sum(0.5 * areas[:3]))

        count, volume = update_conserved_quantities(
            domain, 0.0, update_domain_c_struct=True)

        self.assertEqual(count, 3)
        self.assertTrue(num.isclose(volume, expected_volume))
        # cells were clamped up to the bed (zero depth)
        self.assertTrue(num.allclose(stage_c[:3], bed_c[:3]))


def rain_on_a_dry_valley(algorithm, mode, n=40):
    """Rain on a dry, tilted V valley with walls all round (issue #530)."""
    domain = anuga.rectangular_cross_domain(n, n, len1=400.0, len2=400.0)
    domain.set_compute_mode(mode)
    domain.set_flow_algorithm(algorithm)
    domain.set_store(False)
    domain.set_quantity('elevation',
                        lambda x, y: 0.02*num.abs(x - 200.0) + 0.01*y)
    domain.set_quantity('friction', 0.03)
    domain.set_quantity('stage', expression='elevation')

    operator = anuga.Rate_operator(domain, rate=50.0e-3/3600.0)

    Br = anuga.Reflective_boundary(domain)
    domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
    return domain, operator


@pytest.mark.parametrize('algorithm', ['DE1', 'DE2'])
def test_runge_kutta_average_comes_before_negative_cells_are_fixed(algorithm):
    """Rain on a dry slope must not create water (issue #530).

    The first step on a dry domain runs to the first yield, and leaves a thin
    film at rest. The timestep of the next step comes from that film, the
    water accelerates down the slope within the step, and the later
    Runge-Kutta stages drive some cells negative. Averaging those stages with
    the start of the step brings the depths back above zero, so nothing is
    lost, but only if the negative depths are not set to zero first.

    The Euler start-up is turned off: it would hide the stages under test.
    """

    domain, operator = rain_on_a_dry_valley(algorithm, 'legacy')
    domain.set_euler_startup_steps(0)

    with warnings.catch_warnings():
        warnings.simplefilter('error')    # no "negative cells" warning
        for t in domain.evolve(yieldstep=600.0, finaltime=660.0):
            pass

    rained = operator.cumulative_influx
    assert rained > 1000.0
    assert abs(domain.get_water_volume() - rained) < 1.0e-9*rained


@pytest.mark.parametrize('mode', ['legacy', 'unified'])
@pytest.mark.parametrize('algorithm', ['DE1', 'DE2'])
def test_euler_startup_keeps_rain_on_a_dry_slope_conservative(algorithm, mode):
    """With a short yield step the averaging above is not enough: the
    Runge-Kutta stages reuse a timestep taken from a film at rest, and cells
    stay negative. Starting the run with Euler steps avoids it (issue #530).
    """

    def gain(startup_steps):
        domain, operator = rain_on_a_dry_valley(algorithm, mode)
        if startup_steps is not None:
            domain.set_euler_startup_steps(startup_steps)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            for t in domain.evolve(yieldstep=50.0, finaltime=100.0):
                pass
        assert domain.get_timestepping_method() in ('rk2', 'rk3')
        return domain.get_water_volume() - operator.cumulative_influx

    assert abs(gain(None)) < 1.0e-8       # the default start-up
    assert gain(0) > 1.0                  # what it prevents


def test_euler_startup_steps_setting():
    domain = anuga.rectangular_cross_domain(4, 4)
    assert domain.get_euler_startup_steps() == 10

    domain.set_euler_startup_steps(3)
    assert domain.get_euler_startup_steps() == 3

    with pytest.raises(ValueError):
        domain.set_euler_startup_steps(-1)


@pytest.mark.parametrize('mode', ['legacy', 'unified'])
def test_euler_startup_counts_steps_over_the_life_of_the_domain(mode):
    """Three start-up steps, then Runge-Kutta, across two calls of evolve;
    the method the user set is the one seen at every yield."""

    domain = anuga.rectangular_cross_domain(6, 6)
    domain.set_compute_mode(mode)
    domain.set_flow_algorithm('DE1')
    domain.set_store(False)
    domain.set_quantity('elevation', 0.0)
    domain.set_quantity('stage', lambda x, y: 1.0 + 0.1*x)
    Br = anuga.Reflective_boundary(domain)
    domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
    domain.set_euler_startup_steps(3)

    used = []
    def count(name):
        original = getattr(domain, name)

        def counted(yieldstep, finaltime):
            used.append(name)
            return original(yieldstep, finaltime)
        setattr(domain, name, counted)

    count('evolve_one_euler_step')
    count('evolve_one_rk2_step')

    for t in domain.evolve(yieldstep=0.02, finaltime=0.04):
        assert domain.get_timestepping_method() == 'rk2'
    for t in domain.evolve(yieldstep=0.02, finaltime=0.2):
        assert domain.get_timestepping_method() == 'rk2'

    assert len(used) > 6
    assert used[:3] == ['evolve_one_euler_step']*3
    assert set(used[3:]) == {'evolve_one_rk2_step'}


if __name__ == '__main__':
    unittest.main()
