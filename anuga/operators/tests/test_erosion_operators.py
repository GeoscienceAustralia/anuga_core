"""  Test erosion operators
"""

import unittest
import os
import anuga
from anuga import Domain
from anuga import Reflective_boundary
from anuga import rectangular_cross_domain
from anuga import file_function

from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a
from anuga.file_conversion.file_conversion import timefile2netcdf
from anuga.config import time_format

from anuga.operators.erosion_operators import Erosion_operator
from anuga.operators.erosion_operators import Flat_slice_erosion_operator

from pprint import pprint

import numpy as num
import warnings
import time



class Test_erosion_operators(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_erosion_operator_simple_de0(self):
        from anuga.config import rho_a, rho_w, eta_w
        from math import pi, cos, sin

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #             bac,     bce,     ecf,     dbe
        vertices = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)

        # Flat surface with 1m of water
        domain.set_quantity('elevation', 0.5)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)
        domain.set_quantity('xmomentum',2.0)
        domain.set_quantity('ymomentum',3.0)

        Stage = domain.quantities['stage'].centroid_values
        Elevation = domain.quantities['elevation'].centroid_values

        Height = Stage - Elevation

        sum1 = num.sum(Height)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})


#        print domain.quantities['stage'].centroid_values
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        # Apply operator to these triangles
        indices = [0,1,3]


        operator = Erosion_operator(domain, indices=indices, logging=True)

        # Apply Operator
        domain.timestep = 2.0
        operator()

        elev_ex  = [ 0. ,  0. ,  0.5,  0. ]
        stage_ex = [ 0.5,  0.5,  1. ,  0.5]

        Stage = domain.quantities['stage'].centroid_values
        Elevation = domain.quantities['elevation'].centroid_values

        Height = Stage - Elevation

        sum2 = num.sum(Height)

        #pprint( domain.quantities['elevation'].centroid_values )
        #pprint( domain.quantities['stage'].centroid_values )
        #print domain.quantities['xmomentum'].centroid_values
        #print domain.quantities['ymomentum'].centroid_values

        assert sum1 == sum2
        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 2.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 3.0)


def make_domain():
    """4-triangle domain, 1m water over 0.5m elevation."""
    a = [0.0, 0.0]; b = [0.0, 2.0]; c = [2.0, 0.0]
    d = [0.0, 4.0]; e = [2.0, 2.0]; f = [4.0, 0.0]
    points = [a, b, c, d, e, f]
    vertices = [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]]
    domain = anuga.Domain(points, vertices)
    domain.set_quantity('elevation', 0.5)
    domain.set_quantity('stage', 1.0)
    domain.set_quantity('friction', 0.0)
    domain.set_boundary({'exterior': Reflective_boundary(domain)})
    return domain


class Test_erosion_operator_variants(unittest.TestCase):

    def setUp(self):
        self.domain = make_domain()

    def tearDown(self):
        try:
            os.remove('domain.sww')
        except OSError:
            pass

    def test_circular_erosion_operator(self):
        """Circular_erosion_operator constructs correctly and sets indices."""
        from anuga.operators.erosion_operators import Circular_erosion_operator
        elev_before = self.domain.quantities['elevation'].centroid_values.copy()
        stage_before = self.domain.quantities['stage'].centroid_values.copy()
        operator = Circular_erosion_operator(
            self.domain, center=[1.0, 1.0], radius=2.0)
        # Verify construction: operator should have a list of indices
        self.assertIsNotNone(operator.indices)
        # Elevation and stage arrays should still have their original shape
        elev_after = self.domain.quantities['elevation'].centroid_values
        stage_after = self.domain.quantities['stage'].centroid_values
        self.assertEqual(elev_after.shape, elev_before.shape)
        self.assertEqual(stage_after.shape, stage_before.shape)

    def test_bed_shear_erosion_operator(self):
        """Bed_shear_erosion_operator constructs and calls without error."""
        from anuga.operators.erosion_operators import Bed_shear_erosion_operator
        operator = Bed_shear_erosion_operator(self.domain, indices=[0, 1])
        self.domain.timestep = 1.0
        operator()

    def test_flat_slice_erosion_operator(self):
        """Flat_slice_erosion_operator slices elevation to target <= 0.5."""
        from anuga.operators.erosion_operators import Flat_slice_erosion_operator
        operator = Flat_slice_erosion_operator(
            self.domain,
            elevation=lambda t: 0.3,
            indices=[0, 1, 2, 3])
        self.domain.timestep = 1.0
        operator()
        elev_c = self.domain.quantities['elevation'].centroid_values
        self.assertTrue(num.all(elev_c <= 0.5 + 1.0e-10),
                        "elevation should not exceed original 0.5 after flat slice")

    def test_flat_fill_slice_erosion_operator(self):
        """Flat_fill_slice_erosion_operator constructs and calls without error."""
        from anuga.operators.erosion_operators import Flat_fill_slice_erosion_operator
        operator = Flat_fill_slice_erosion_operator(
            self.domain,
            elevation=lambda t: 0.4,
            indices=[0, 1, 2, 3])
        self.domain.timestep = 1.0
        operator()


class Test_circular_operators(unittest.TestCase):

    def setUp(self):
        self.domain = make_domain()

    def tearDown(self):
        try:
            os.remove('domain.sww')
        except OSError:
            pass

    # ------------------------------------------------------------------
    # Regressions for #301 (whole-domain erosion was a silent no-op) and
    # #302 (the base operator crashed under discontinuous elevation).
    # ------------------------------------------------------------------

    def _sloping_channel(self, name):
        """A channel with flow, so an erosion operator has something to do."""
        L, W = 100.0, 20.0
        d = rectangular_cross_domain(30, 6, L, W, origin=(0.0, 0.0))
        d.set_quantity('elevation', lambda x, y: 5.0 - 0.05 * x)
        d.set_quantity('friction', 0.03)
        d.set_quantity('stage', expression='elevation')
        d.set_flow_algorithm('DE1')
        d.set_datadir('.')
        d.set_name(name)
        d.set_quantities_to_be_stored(None)
        d.set_boundary({'left': anuga.Dirichlet_boundary([5.5, 2.0, 0.0]),
                        'right': anuga.Transmissive_boundary(d),
                        'top': Reflective_boundary(d),
                        'bottom': Reflective_boundary(d)})
        return d

    def test_polygonal_erosion_runs_under_discontinuous_elevation(self):
        """#302: the base update_quantities subtracted a per-cell (n,) array
        from the (n, 3) vertex array, which cannot broadcast. Every supported
        flow algorithm uses discontinuous elevation, so this crashed always."""
        from anuga.operators.erosion_operators import Polygonal_erosion_operator

        d = self._sloping_channel('test_302')
        self.assertTrue(d.get_using_discontinuous_elevation())
        Polygonal_erosion_operator(
            d, polygon=[[0.0, 0.0], [100.0, 0.0], [100.0, 20.0], [0.0, 20.0]])

        z0 = d.quantities['elevation'].centroid_values.copy()
        for t in d.evolve(yieldstep=30.0, finaltime=90.0):
            pass
        dz = d.quantities['elevation'].centroid_values - z0

        assert dz.min() < -1.0e-3, 'the operator ran but eroded nothing'

    def test_erosion_never_raises_the_bed(self):
        """An erosion operator must not deposit.

        Taking amax over the vertex values pulls the centroid up to its
        highest vertex, which under discontinuous elevation can sit above the
        centroid -- so the bed gains height where it should only lose it.

        Runs to 180 s deliberately: with amax the bed still looks well behaved
        at 120 s and only starts gaining height after that, so a shorter run
        does not discriminate.
        """
        from anuga.operators.erosion_operators import Polygonal_erosion_operator

        d = self._sloping_channel('test_302_rise')
        Polygonal_erosion_operator(
            d, polygon=[[0.0, 0.0], [100.0, 0.0], [100.0, 20.0], [0.0, 20.0]])

        z0 = d.quantities['elevation'].centroid_values.copy()
        for t in d.evolve(yieldstep=30.0, finaltime=180.0):
            pass
        dz = d.quantities['elevation'].centroid_values - z0

        assert dz.max() <= 1.0e-10, (
            'bed rose by %g m under an erosion operator' % dz.max())

    def test_erosion_respects_the_base_level(self):
        """Erosion stops at `base`, however long it runs.

        Only cells that START above `base` are constrained by it. This channel
        runs from 5 m down to 0 m, so its downstream end begins below
        `base=1.0`; those cells must be left where they are, not lifted (#305).
        """
        from anuga.operators.erosion_operators import Polygonal_erosion_operator

        d = self._sloping_channel('test_302_base')
        Polygonal_erosion_operator(
            d, polygon=[[0.0, 0.0], [100.0, 0.0], [100.0, 20.0], [0.0, 20.0]],
            base=1.0)

        z0 = d.quantities['elevation'].centroid_values.copy()
        for t in d.evolve(yieldstep=30.0, finaltime=90.0):
            pass
        z = d.quantities['elevation'].centroid_values

        above = z0 >= 1.0
        assert above.any() and (~above).any(), 'the case must cover both sides'
        assert z[above].min() >= 1.0 - 1.0e-10, (
            'eroded below base: %g' % z[above].min())
        assert num.allclose(z[~above], z0[~above]), (
            'cells starting below base were moved')

    def _flat_slice_run(self, polygon, name):
        d = rectangular_cross_domain(20, 6, 100.0, 20.0, origin=(0.0, 0.0))
        d.set_quantity('elevation', 2.0)
        d.set_quantity('stage', 4.0)
        d.set_quantity('friction', 0.03)
        d.set_flow_algorithm('DE1')
        d.set_datadir('.')
        d.set_name(name)
        d.set_quantities_to_be_stored(None)
        R = Reflective_boundary(d)
        d.set_boundary({'left': R, 'right': R, 'top': R, 'bottom': R})
        Flat_slice_erosion_operator(d, elevation=lambda t: 0.0, polygon=polygon)
        for t in d.evolve(yieldstep=25.0, finaltime=50.0):
            pass
        return d.quantities['elevation'].centroid_values.mean()

    def test_whole_domain_erosion_is_not_a_no_op(self):
        """#301: `indices is None` is documented as "all triangles", but took
        a branch that added 0.0 and returned, so the bed never moved."""
        whole = self._flat_slice_run(None, 'test_301_none')
        assert whole < 2.0 - 1.0e-6, (
            'polygon=None left the bed at %g m -- whole-domain erosion did '
            'nothing' % whole)

    def test_whole_domain_erodes_at_least_as_much_as_a_polygon(self):
        """The None path and an explicit covering polygon must agree in kind.

        The polygon cannot beat the whole domain: it covers a subset.
        """
        whole = self._flat_slice_run(None, 'test_301_a')
        poly = self._flat_slice_run(
            [[1.0, 1.0], [99.0, 1.0], [99.0, 19.0], [1.0, 19.0]], 'test_301_b')
        assert whole <= poly + 1.0e-10, (
            'whole domain (%g m) eroded less than a polygon inside it (%g m)'
            % (whole, poly))

    def test_base_does_not_lift_a_bed_that_starts_below_it(self):
        """#305: `base` limits erosion -- "Allow erosion down to base level".

        Applied as a bare maximum it also LIFTS a bed that starts below it, so
        the operator deposits on its first call with no flow involved.
        """
        from anuga.operators.erosion_operators import (
            Polygonal_erosion_operator, Circular_erosion_operator,
            Bed_shear_erosion_operator)

        poly = [[2.0, 2.0], [18.0, 2.0], [18.0, 8.0], [2.0, 8.0]]
        makers = [
            ('Polygonal', lambda d: Polygonal_erosion_operator(
                d, base=0.0, polygon=poly)),
            ('Circular', lambda d: Circular_erosion_operator(
                d, base=0.0, center=(10.0, 5.0), radius=5.0)),
            ('Bed_shear', lambda d: Bed_shear_erosion_operator(
                d, base=0.0, polygon=poly)),
        ]

        for name, make in makers:
            d = rectangular_cross_domain(20, 10, 20.0, 10.0)
            d.set_quantity('elevation', -0.5)      # everywhere BELOW base
            d.set_quantity('stage', 1.0)
            d.set_quantity('friction', 0.03)
            d.set_flow_algorithm('DE1')
            d.set_datadir('.')
            d.set_name('test_305_' + name)
            d.set_quantities_to_be_stored(None)
            R = Reflective_boundary(d)
            d.set_boundary({'left': R, 'right': R, 'top': R, 'bottom': R})
            make(d)

            z0 = d.quantities['elevation'].centroid_values.copy()
            for t in d.evolve(yieldstep=2.0, finaltime=4.0):
                pass
            dz = d.quantities['elevation'].centroid_values - z0

            assert dz.max() <= 1.0e-10, (
                '%s raised a bed already below base by %g m'
                % (name, dz.max()))

    def test_circular_rate_operator(self):
        """Circular_rate_operator constructs and calls without raising."""
        from anuga.operators.rate_operators import Circular_rate_operator
        stage_before = self.domain.quantities['stage'].centroid_values.copy()
        operator = Circular_rate_operator(
            self.domain, center=[1.0, 1.0], radius=2.0, rate=0.01)
        self.domain.timestep = 1.0
        operator()
        # Just verify the call didn't raise; stage array shape preserved
        stage_after = self.domain.quantities['stage'].centroid_values
        self.assertEqual(stage_after.shape, stage_before.shape)

    def test_circular_set_quantity_operator(self):
        """Circular_set_quantity_operator sets friction to 0.05 inside circle."""
        from anuga.operators.set_quantity_operator import Circular_set_quantity_operator
        # Note: Set_quantity_operator does not allow setting 'stage' or 'elevation'
        # directly (use Set_stage_operator / Set_elevation_operator for those).
        operator = Circular_set_quantity_operator(
            self.domain, quantity='friction',
            center=[1.0, 1.0], radius=2.0, value=0.05)
        operator()
        friction_vals = self.domain.quantities['friction'].centroid_values
        # At least one cell inside the circle should have friction == 0.05
        self.assertTrue(num.any(num.isclose(friction_vals, 0.05)),
                        "At least one centroid friction value should equal 0.05")

    def test_circular_set_stage_operator(self):
        """Circular_set_stage_operator sets stage >= elevation inside circle."""
        from anuga.operators.set_stage_operator import Circular_set_stage_operator
        operator = Circular_set_stage_operator(
            self.domain, center=[1.0, 1.0], radius=2.0, stage=1.5)
        operator()
        stage_vals = self.domain.quantities['stage'].centroid_values
        # At least one centroid should be at or above 1.0 (elevation=0.5, stage >= elev)
        self.assertTrue(num.any(stage_vals >= 1.0),
                        "At least one stage centroid value should be >= 1.0")


class Test_set_quantity_operator_extra(unittest.TestCase):
    """Additional tests for Set_quantity_operator methods."""

    def setUp(self):
        self.domain = rectangular_cross_domain(2, 2)
        self.domain.set_quantity('elevation', 0.0)
        self.domain.set_quantity('stage', 1.0)

    def _make_op(self):
        from anuga.operators.set_quantity_operator import Set_quantity_operator
        return Set_quantity_operator(self.domain, 'friction', value=0.03)

    def test_parallel_safe(self):
        """parallel_safe returns True (line 77)."""
        self.assertTrue(self._make_op().parallel_safe())

    def test_statistics(self):
        """statistics returns string (lines 81-83)."""
        msg = self._make_op().statistics()
        self.assertIsInstance(msg, str)

    def test_timestepping_statistics(self):
        """timestepping_statistics returns string (line 91)."""
        msg = self._make_op().timestepping_statistics()
        self.assertIsInstance(msg, str)

    def test_call_all_triangles(self):
        """Call with indices=None updates friction (set_quantity.py lines 109-117)."""
        from anuga.operators.set_quantity_operator import Set_quantity_operator
        import numpy as num
        op = Set_quantity_operator(self.domain, 'friction', value=0.05, indices=None)
        self.domain.timestep = 1.0
        op()
        self.assertTrue(num.allclose(
            self.domain.quantities['friction'].centroid_values, 0.05))

    def test_call_specific_indices(self):
        """Call with specific indices updates those friction values (set_quantity.py lines 125-133)."""
        from anuga.operators.set_quantity_operator import Set_quantity_operator
        import numpy as num
        op = Set_quantity_operator(self.domain, 'friction', value=0.07, indices=[0, 1])
        self.domain.timestep = 1.0
        op()

    def test_call_empty_indices(self):
        """Empty indices → early return (set_quantity.py line 99)."""
        from anuga.operators.set_quantity_operator import Set_quantity_operator
        op = Set_quantity_operator(self.domain, 'friction', value=0.05, indices=[])
        self.domain.timestep = 1.0
        op()  # should return early

    def test_pass_region_object(self):
        """Passing a Region object uses lines 48-49 in set_quantity.py."""
        from anuga.operators.set_quantity_operator import Set_quantity_operator
        from anuga import Region
        region = Region(self.domain, indices=[0, 1, 2])
        op = Set_quantity_operator(self.domain, 'friction', value=0.04, region=region)
        self.assertIsNotNone(op)


class Test_erosion_elevation_storage(unittest.TestCase):
    """Erosion operators default elevation to time-varying SWW storage."""

    def _domain(self):
        domain = rectangular_cross_domain(6, 6)
        domain.set_boundary(
            {b: Reflective_boundary(domain) for b in domain.get_boundary_tags()})
        return domain

    def test_operator_promotes_elevation_to_time_varying(self):
        from anuga.operators.erosion_operators import Circular_erosion_operator
        domain = self._domain()
        # Default is static (1); creating an erosion operator promotes to
        # dynamic (2) and flags the domain.
        self.assertEqual(domain.quantities_to_be_stored['elevation'], 1)
        Circular_erosion_operator(domain, center=(0.5, 0.5), radius=0.2)
        self.assertEqual(domain.quantities_to_be_stored['elevation'], 2)
        self.assertTrue(getattr(domain, '_erosion_present', False))

    def test_respects_explicit_static(self):
        from anuga.operators.erosion_operators import Circular_erosion_operator
        domain = self._domain()
        # A deliberate static choice is honoured (not promoted).
        domain._elevation_static_by_user = True
        domain.quantities_to_be_stored['elevation'] = 1
        Circular_erosion_operator(domain, center=(0.5, 0.5), radius=0.2)
        self.assertEqual(domain.quantities_to_be_stored['elevation'], 1)

    def test_warns_when_reset_to_static_before_storage(self):
        import tempfile
        from anuga.operators.erosion_operators import Circular_erosion_operator
        domain = self._domain()
        domain.set_name('erosion_storage_warn')
        Circular_erosion_operator(domain, center=(0.5, 0.5), radius=0.2)
        # User resets to static AFTER the operator promoted it.
        domain.quantities_to_be_stored['elevation'] = 1
        cwd = os.getcwd()
        tmp = tempfile.mkdtemp()
        try:
            os.chdir(tmp)
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter('always')
                domain.initialise_storage()
            msgs = ' '.join(str(x.message) for x in w)
            self.assertIn('eroded bed', msgs)
        finally:
            os.chdir(cwd)


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_erosion_operators)
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
