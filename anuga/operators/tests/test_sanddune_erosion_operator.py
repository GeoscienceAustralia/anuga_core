"""Tests for Sanddune_erosion_operator.

Regression coverage for the two crashes fixed on branch
``fix/sanddune-erosion-base-scalar``:

1. A default scalar ``base`` was indexed as an array (``self.base[ind]``),
   raising ``TypeError: 'float' object is not subscriptable``.
2. Whole-domain application (``indices=None``) passed ``None`` straight into
   NumPy fancy indexing (``neighbours[None]``), which inserts a new axis
   instead of selecting everything and corrupted the neighbour arrays,
   eventually surfacing as a broadcasting ``ValueError``.

Both are exercised together by applying the operator to the whole domain with
the default scalar ``base`` and running it through a timestep.
"""

import unittest

import numpy as num

import anuga
from anuga import Reflective_boundary
from anuga.operators.sanddune_erosion_operator import Sanddune_erosion_operator


def make_domain():
    """4-triangle domain with moving water over a 0.5 m bed."""
    a = [0.0, 0.0]; b = [0.0, 2.0]; c = [2.0, 0.0]
    d = [0.0, 4.0]; e = [2.0, 2.0]; f = [4.0, 0.0]
    points = [a, b, c, d, e, f]
    vertices = [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]]
    domain = anuga.Domain(points, vertices)
    domain.set_quantity('elevation', 0.5)
    domain.set_quantity('stage', 1.0)
    domain.set_quantity('friction', 0.0)
    domain.set_quantity('xmomentum', 2.0)
    domain.set_quantity('ymomentum', 3.0)
    domain.set_boundary({'exterior': Reflective_boundary(domain)})
    return domain


class Test_sanddune_erosion_operator(unittest.TestCase):

    def setUp(self):
        self.domain = make_domain()

    def test_base_does_not_lift_a_bed_that_starts_below_it(self):
        """#305: `base` is a limit on erosion, not a level to restore to.

        Applied as a bare maximum it lifts a bed that already sits below it.
        """
        d = anuga.rectangular_cross_domain(20, 10, 20.0, 10.0)
        d.set_quantity('elevation', -0.5)          # everywhere BELOW base
        d.set_quantity('stage', 1.0)
        d.set_quantity('friction', 0.03)
        d.set_flow_algorithm('DE1')
        d.set_datadir('.')
        d.set_name('test_305_sanddune')
        d.set_quantities_to_be_stored(None)
        R = Reflective_boundary(d)
        d.set_boundary({'left': R, 'right': R, 'top': R, 'bottom': R})
        Sanddune_erosion_operator(
            d, base=0.0,
            polygon=[[2.0, 2.0], [18.0, 2.0], [18.0, 8.0], [2.0, 8.0]])

        z0 = d.quantities['elevation'].centroid_values.copy()
        for t in d.evolve(yieldstep=2.0, finaltime=4.0):
            pass
        dz = d.quantities['elevation'].centroid_values - z0

        assert dz.max() <= 1.0e-10, (
            'raised a bed already below base by %g m' % dz.max())

    def test_whole_domain_scalar_base_runs(self):
        """indices=None + default scalar base: previously crashed, now runs.

        Covers both fixed bugs at once (scalar-base indexing and the
        None-as-newaxis fancy-index corruption).
        """
        n = self.domain.number_of_elements
        stage_c = self.domain.quantities['stage'].centroid_values
        elev_c = self.domain.quantities['elevation'].centroid_values
        height_before = num.sum(stage_c - elev_c)

        operator = Sanddune_erosion_operator(self.domain)   # indices=None, base=0.0

        self.domain.timestep = 1.0
        operator()   # must not raise

        # Water column (height) is preserved: stage is rebuilt as elev + depth.
        height_after = num.sum(
            self.domain.quantities['stage'].centroid_values
            - self.domain.quantities['elevation'].centroid_values)
        self.assertTrue(num.allclose(height_before, height_after))

        # Erosion never cuts below the base level (0.0 here) anywhere.
        self.assertTrue(
            num.all(self.domain.quantities['elevation'].centroid_values
                    >= 0.0 - 1.0e-10))
        # Shapes untouched.
        self.assertEqual(
            self.domain.quantities['elevation'].centroid_values.shape, (n,))

    def test_scalar_base_expanded_to_full_array(self):
        """A scalar base is broadcast to a per-element array of the right length."""
        n = self.domain.number_of_elements
        operator = Sanddune_erosion_operator(self.domain, base=0.25)
        self.assertEqual(num.shape(operator.base), (n,))
        self.assertTrue(num.allclose(operator.base, 0.25))

    def test_per_element_base_array_accepted(self):
        """A correctly sized per-element base array is accepted as-is."""
        n = self.domain.number_of_elements
        base = num.linspace(0.0, 0.2, n)
        operator = Sanddune_erosion_operator(self.domain, base=base)
        self.assertTrue(num.allclose(operator.base, base))

    def test_wrong_length_base_array_raises(self):
        """A base array whose length != number_of_elements is rejected."""
        n = self.domain.number_of_elements
        with self.assertRaises(ValueError):
            Sanddune_erosion_operator(self.domain, base=num.zeros(n + 1))

    def test_base_level_is_respected(self):
        """Erosion cannot lower elevation below a raised base level."""
        base_level = 0.5   # equal to the bed, so no erosion is permitted
        operator = Sanddune_erosion_operator(self.domain, base=base_level)
        self.domain.timestep = 5.0
        operator()
        self.assertTrue(
            num.all(self.domain.quantities['elevation'].centroid_values
                    >= base_level - 1.0e-10))

    def test_subset_indices_with_full_domain_base_array(self):
        """Primary usage (see anuga-clinic notebook3): a full-domain base array
        applied to a subset of triangles selected by indices.

        Also a regression for the ``self.indices != []`` guard, which raised a
        broadcast ValueError under NumPy 2.x once Region turned a non-empty
        indices list into an ndarray.
        """
        # Notebook pattern: base is the initial (full-domain) elevation array.
        base = self.domain.get_quantity('elevation').get_values(
            location='centroids').copy()
        indices = [0, 1]
        operator = Sanddune_erosion_operator(
            self.domain, base=base, indices=indices, Ra=45)

        elev_before = self.domain.quantities['elevation'].centroid_values.copy()
        self.domain.timestep = 1.0
        operator()   # must not raise
        elev_after = self.domain.quantities['elevation'].centroid_values

        # Triangles outside the erosion indices are left untouched.
        outside = [i for i in range(self.domain.number_of_elements)
                   if i not in indices]
        self.assertTrue(num.allclose(elev_before[outside], elev_after[outside]))

    def test_empty_indices_is_a_noop(self):
        """indices=[] means "apply nowhere": nothing changes and nothing raises."""
        operator = Sanddune_erosion_operator(self.domain, indices=[])
        elev_before = self.domain.quantities['elevation'].centroid_values.copy()
        stage_before = self.domain.quantities['stage'].centroid_values.copy()
        self.domain.timestep = 1.0
        operator()
        self.assertTrue(num.allclose(
            elev_before, self.domain.quantities['elevation'].centroid_values))
        self.assertTrue(num.allclose(
            stage_before, self.domain.quantities['stage'].centroid_values))


if __name__ == '__main__':
    unittest.main()
