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


if __name__ == '__main__':
    unittest.main()
