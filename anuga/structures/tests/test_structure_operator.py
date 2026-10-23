#!/usr/bin/env python

import unittest
import warnings

import numpy
import anuga
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.structures.inlet import level_stages_to_average
from anuga.structures.inlet_enquiry import Inlet_enquiry


verbose = False

# This end-point geometry places enquiry points inside inlet triangles on
# this small test mesh — expected behaviour that raises a UserWarning.
_INLET_WARNING = 'Enquiry point.*is in an inlet triangle'


def make_domain():
    """Create a simple 10m x 5m rectangular domain for testing."""
    points, vertices, boundary = rectangular_cross(10, 5, len1=10.0, len2=5.0)
    domain = Domain(points, vertices, boundary)
    domain.set_quantity('elevation', 0.0)
    domain.set_quantity('stage', 1.0)
    domain.set_quantity('friction', 0.0)
    Br = anuga.Reflective_boundary(domain)
    domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
    return domain


class Test_Structure_operator(unittest.TestCase):
    """Tests for the Structure_operator base class."""

    def setUp(self):
        self.domain = make_domain()
        self._warning_ctx = warnings.catch_warnings()
        self._warning_ctx.__enter__()
        warnings.filterwarnings('ignore', message=_INLET_WARNING,
                                category=UserWarning)

    def tearDown(self):
        self._warning_ctx.__exit__(None, None, None)

    def _make_operator(self):
        """Helper: create a Structure_operator with all required parameters.

        The culvert runs along the x-axis so that auto-computed enquiry points
        stay inside the 10m x 5m domain.
        """
        return anuga.Structure_operator(
            self.domain,
            end_points=[[3., 2.5], [7., 2.5]],
            width=1.0,
            manning=0.013,
            enquiry_gap=0.0,
            verbose=verbose)

    def test_construction(self):
        """Structure_operator construction warns when enquiry points are in
        inlet triangles (expected for this test mesh geometry)."""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            op = anuga.Structure_operator(
                self.domain,
                end_points=[[3., 2.5], [7., 2.5]],
                width=1.0,
                manning=0.013,
                enquiry_gap=0.0,
                verbose=verbose)
        self.assertIsNotNone(op)
        inlet_warnings = [w for w in caught
                          if issubclass(w.category, UserWarning)
                          and 'inlet triangle' in str(w.message)]
        self.assertGreater(len(inlet_warnings), 0,
                           'Expected inlet-triangle UserWarning was not raised')

    def test_get_culvert_length(self):
        """get_culvert_length returns a positive value."""
        op = self._make_operator()
        length = op.get_culvert_length()
        self.assertGreater(length, 0.0)

    def test_get_culvert_width(self):
        """get_culvert_width returns the specified width."""
        op = self._make_operator()
        self.assertAlmostEqual(op.get_culvert_width(), 1.0)

    def test_repr_returns_string(self):
        """str() on the operator returns a string (via __repr__ or __str__)."""
        op = self._make_operator()
        result = str(op)
        self.assertIsInstance(result, str)

    def test_statistics_returns_string(self):
        """statistics() returns a non-empty string."""
        op = self._make_operator()
        result = op.statistics()
        self.assertIsInstance(result, str)
        self.assertGreater(len(result), 0)

    def test_discharge_routine_raises(self):
        """Base class discharge_routine raises (NotImplementedError or similar)."""
        op = self._make_operator()
        with self.assertRaises(Exception):
            op.discharge_routine()

    def test_constructor_momentum_jet_and_zero_outflow_raises(self):
        """use_momentum_jet=True with zero_outflow_momentum=True raises Exception."""
        with self.assertRaises(Exception):
            anuga.Structure_operator(
                self.domain,
                end_points=[[3., 2.5], [7., 2.5]],
                width=1.0,
                use_momentum_jet=True,
                zero_outflow_momentum=True,
                verbose=verbose)

    def test_constructor_no_geometry_raises(self):
        """Omitting both exchange_lines and end_points raises Exception."""
        with self.assertRaises(Exception):
            anuga.Structure_operator(
                self.domain,
                width=1.0,
                verbose=verbose)

    def test_description_parameter(self):
        """Non-None description is stored verbatim."""
        op = anuga.Structure_operator(
            self.domain,
            end_points=[[3., 2.5], [7., 2.5]],
            width=1.0,
            enquiry_gap=0.0,
            description='my culvert',
            verbose=verbose)
        self.assertEqual(op.description, 'my culvert')

    def test_setters(self):
        """set_culvert_* setters store the given values."""
        op = self._make_operator()
        op.set_culvert_height(2.0)
        self.assertEqual(op.culvert_height, 2.0)
        op.set_culvert_width(1.5)
        self.assertEqual(op.culvert_width, 1.5)
        op.set_culvert_z1(0.5)
        self.assertEqual(op.culvert_z1, 0.5)
        op.set_culvert_z2(0.3)
        self.assertEqual(op.culvert_z2, 0.3)
        op.set_culvert_blockage(0.1)
        self.assertEqual(op.culvert_blockage, 0.1)
        op.set_culvert_barrels(2)
        self.assertEqual(op.culvert_barrels, 2)

    def test_enquiry_points_provided(self):
        """Providing explicit enquiry_points skips auto-computation (lines 425-429)."""
        op = anuga.Structure_operator(
            self.domain,
            end_points=[[3., 2.5], [7., 2.5]],
            width=1.0,
            enquiry_points=[[2.0, 2.5], [8.0, 2.5]],
            verbose=verbose)
        self.assertIsNotNone(op)

    def test_skew_culvert_4point_exchange_lines(self):
        """4-point exchange lines trigger the n_exchange_0==4 branch."""
        el0 = numpy.array([[2., 2.], [2., 3.], [2., 2.5], [3., 2.5]])
        el1 = numpy.array([[7., 2.], [7., 3.], [7., 2.5], [8., 2.5]])
        op = anuga.Structure_operator(
            self.domain,
            exchange_lines=[el0, el1],
            width=1.0,
            enquiry_gap=0.0,
            verbose=verbose)
        self.assertIsNotNone(op)
        self.assertGreater(op.culvert_length, 0.0)

    def test_print_statistics(self):
        """print_statistics() writes the statistics string to stdout."""
        import io
        import sys
        op = self._make_operator()
        captured = io.StringIO()
        sys.stdout = captured
        try:
            op.print_statistics()
        finally:
            sys.stdout = sys.__stdout__
        self.assertGreater(len(captured.getvalue()), 0)

    def test_timestepping_statistics(self):
        """timestepping_statistics() returns a comma-separated string."""
        op = self._make_operator()
        op.case = 'test'
        result = op.timestepping_statistics()
        self.assertIsInstance(result, str)
        self.assertIn(',', result)

    def test_print_timestepping_statistics(self):
        """print_timestepping_statistics() writes to stdout."""
        import io
        import sys
        op = self._make_operator()
        op.case = 'test'
        captured = io.StringIO()
        sys.stdout = captured
        try:
            op.print_timestepping_statistics()
        finally:
            sys.stdout = sys.__stdout__
        self.assertGreater(len(captured.getvalue()), 0)

    def test_get_culvert_apron(self):
        """get_culvert_apron returns the apron value (defaults to width)."""
        op = self._make_operator()
        apron = op.get_culvert_apron()
        self.assertAlmostEqual(apron, 1.0)

    def test_get_master_proc(self):
        """get_master_proc returns 0 for serial domains."""
        op = self._make_operator()
        self.assertEqual(op.get_master_proc(), 0)

    def test_enquiry_getters(self):
        """All 15 enquiry getter pairs return two-element lists."""
        op = self._make_operator()
        self.assertEqual(len(op.get_enquiry_stages()), 2)
        self.assertEqual(len(op.get_enquiry_depths()), 2)
        self.assertEqual(len(op.get_enquiry_positions()), 2)
        self.assertEqual(len(op.get_enquiry_xmoms()), 2)
        self.assertEqual(len(op.get_enquiry_ymoms()), 2)
        self.assertEqual(len(op.get_enquiry_elevations()), 2)
        self.assertEqual(len(op.get_enquiry_water_depths()), 2)
        self.assertEqual(len(op.get_enquiry_invert_elevations()), 2)
        self.assertEqual(len(op.get_enquiry_velocitys()), 2)
        self.assertEqual(len(op.get_enquiry_xvelocitys()), 2)
        self.assertEqual(len(op.get_enquiry_yvelocitys()), 2)
        self.assertEqual(len(op.get_enquiry_speeds()), 2)
        self.assertEqual(len(op.get_enquiry_velocity_heads()), 2)
        self.assertEqual(len(op.get_enquiry_total_energys()), 2)
        self.assertEqual(len(op.get_enquiry_specific_energys()), 2)

    def test_statistics_non_constant_elevation(self):
        """statistics() warns when inlet triangles have non-uniform elevation."""
        points, vertices, boundary = rectangular_cross(10, 5, len1=10.0, len2=5.0)
        domain = Domain(points, vertices, boundary)
        domain.set_quantity('elevation', lambda x, y: x * 0.1)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0.0)
        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', message=_INLET_WARNING,
                                    category=UserWarning)
            op = anuga.Structure_operator(
                domain,
                end_points=[[3., 2.5], [7., 2.5]],
                width=1.0,
                enquiry_gap=0.0,
                verbose=verbose)
        result = op.statistics()
        self.assertIn('non-constant', result)


class Test_Inlet_enquiry(unittest.TestCase):
    """Tests for the Inlet_enquiry class."""

    def setUp(self):
        self.domain = make_domain()

    def tearDown(self):
        pass

    def test_construction(self):
        """Inlet_enquiry can be constructed without error."""
        region = [[2.5, 0.], [2.5, 2.5]]
        enquiry_pt = [1.5, 2.5]
        inlet = Inlet_enquiry(
            self.domain,
            region=region,
            enquiry_pt=enquiry_pt,
            verbose=verbose)
        self.assertIsNotNone(inlet)

    def test_enquiry_pt_stored(self):
        """enquiry_pt attribute is stored correctly."""
        region = [[2.5, 0.], [2.5, 2.5]]
        enquiry_pt = [1.5, 2.5]
        inlet = Inlet_enquiry(
            self.domain,
            region=region,
            enquiry_pt=enquiry_pt,
            verbose=verbose)
        self.assertTrue(numpy.allclose(inlet.enquiry_pt, enquiry_pt))

    def test_enquiry_index_set(self):
        """enquiry_index is set and is a valid triangle index (>= 0)."""
        region = [[2.5, 0.], [2.5, 2.5]]
        enquiry_pt = [1.5, 2.5]
        inlet = Inlet_enquiry(
            self.domain,
            region=region,
            enquiry_pt=enquiry_pt,
            verbose=verbose)
        self.assertGreaterEqual(inlet.enquiry_index, 0)
        num_triangles = len(self.domain)
        self.assertLess(inlet.enquiry_index, num_triangles)


class Test_level_stages_to_average(unittest.TestCase):
    """The inlet write-back kernel behind Inlet.set_average_depth() (issue #229).

    Fills raise the lowest stages to a common level; drawdowns lower the highest
    stages to a common level, clamping at the bed — the "water finds its level"
    model the structure operators assume, applied to STAGE rather than depth.
    """

    @staticmethod
    def _volume(stages, beds, areas):
        return float(numpy.dot(numpy.maximum(stages - beds, 0.0), areas))

    def test_level_surface_stays_level(self):
        """A level lake over a sloping bed must stay level through a transfer."""
        beds = numpy.array([0.0, 1.0, 2.0, 3.0])
        stages = numpy.full(4, 5.0)
        areas = numpy.ones(4)
        avg = float(numpy.dot(stages - beds, areas) / areas.sum())

        out = level_stages_to_average(stages, beds, areas, avg, avg - 0.5)
        numpy.testing.assert_allclose(out, 4.5, rtol=0, atol=1e-12,
                                      err_msg='drawdown should stay level')

        out = level_stages_to_average(stages, beds, areas, avg, avg + 0.25)
        numpy.testing.assert_allclose(out, 5.25, rtol=0, atol=1e-12,
                                      err_msg='fill should stay level')

    def test_zero_change_is_exactly_a_no_op(self):
        """No transfer must be bit-exact — the well-balanced case."""
        beds = numpy.array([0.0, 1.0, 2.0, 3.0])
        stages = numpy.array([5.0, 4.5, 3.5, 3.2])
        areas = numpy.array([1.0, 2.0, 1.0, 2.0])
        avg = self._volume(stages, beds, areas) / areas.sum()

        out = level_stages_to_average(stages, beds, areas, avg, avg)
        numpy.testing.assert_array_equal(out, stages)

    def test_volume_change_is_exact(self):
        """The volume moved is exactly what the caller asked for."""
        beds = numpy.array([0.0, 1.0, 2.0, 3.0])
        stages = numpy.array([5.0, 4.5, 3.5, 3.2])
        areas = numpy.array([1.0, 2.0, 1.0, 2.0])
        total_area = float(areas.sum())
        avg = self._volume(stages, beds, areas) / total_area

        for delta in (-0.15, -0.02, 0.3, 2.0):
            out = level_stages_to_average(stages, beds, areas, avg, avg + delta)
            moved = (self._volume(out, beds, areas)
                     - self._volume(stages, beds, areas))
            self.assertAlmostEqual(moved, delta * total_area, places=9)

    def test_fill_raises_only_the_lowest(self):
        """Water finds its level: a high cell is never pulled down."""
        beds = numpy.zeros(4)
        stages = numpy.array([5.0, 4.0, 3.0, 2.0])
        areas = numpy.ones(4)
        avg = stages.mean()

        out = level_stages_to_average(stages, beds, areas, avg, avg + 0.5)
        self.assertEqual(out[0], 5.0)          # highest untouched
        self.assertAlmostEqual(out[2], out[3], places=12)  # lowest at one level
        self.assertGreaterEqual(out[1], 4.0)

    def test_drawdown_lowers_only_the_highest(self):
        """A low cell is never lifted by a removal."""
        beds = numpy.zeros(4)
        stages = numpy.array([5.0, 4.0, 3.0, 2.0])
        areas = numpy.ones(4)
        avg = stages.mean()

        out = level_stages_to_average(stages, beds, areas, avg, avg - 0.5)
        self.assertEqual(out[3], 2.0)          # lowest untouched
        self.assertAlmostEqual(out[0], out[1], places=9)  # highest at one level
        numpy.testing.assert_allclose(out, [3.5, 3.5, 3.0, 2.0], atol=1e-9)

    def test_wet_dry_clamp(self):
        """Cells that reach their bed dry out; the deep cells supply the rest."""
        beds = numpy.array([4.9, 4.8, 0.0, 0.0])
        stages = numpy.full(4, 5.0)
        areas = numpy.ones(4)
        avg = self._volume(stages, beds, areas) / areas.sum()

        out = level_stages_to_average(stages, beds, areas, avg, avg - 2.0)
        self.assertAlmostEqual(out[0], 4.9, places=12)   # clamped at bed
        self.assertAlmostEqual(out[1], 4.8, places=12)
        self.assertAlmostEqual(out[2], out[3], places=9)
        moved = (self._volume(out, beds, areas)
                 - self._volume(stages, beds, areas))
        self.assertAlmostEqual(moved, -2.0 * areas.sum(), places=9,
                               msg='clamping must not change the volume moved')

    def test_draining_everything_leaves_it_dry(self):
        """Asking for more water than the inlet holds empties it, no further."""
        beds = numpy.array([0.0, 0.5])
        stages = numpy.array([1.0, 2.0])
        areas = numpy.ones(2)
        avg = self._volume(stages, beds, areas) / areas.sum()

        out = level_stages_to_average(stages, beds, areas, avg, avg - 10.0)
        numpy.testing.assert_array_equal(out, beds)

    def test_flat_bed_fill_matches_old_uniform_depth(self):
        """On a flat bed from level water, leveling == the old set_depths()."""
        beds = numpy.zeros(4)
        stages = numpy.full(4, 0.4)
        areas = numpy.ones(4)

        out = level_stages_to_average(stages, beds, areas, 0.4, 0.7)
        numpy.testing.assert_allclose(out, 0.7, rtol=0, atol=1e-12)
        out = level_stages_to_average(stages, beds, areas, 0.4, 0.1)
        numpy.testing.assert_allclose(out, 0.1, rtol=0, atol=1e-12)


class Test_Structure_well_balanced(unittest.TestCase):
    """A structure must not disturb a lake at rest (issue #229).

    The inlet write-back used to set a uniform DEPTH, which on a sloping bed
    tilts the water surface onto the bed — a lake at rest picked up an error of
    about half the bed elevation range across the inlet, every timestep, with no
    flow through the structure at all.
    """

    def _lake_at_rest(self, slope_denominator, with_culvert=True):
        domain = anuga.rectangular_cross_domain(30, 15, len1=200.0, len2=50.0)
        domain.set_flow_algorithm('DE0')
        domain.set_name('well_balanced')
        domain.store = False
        domain.set_quantity('elevation',
                            lambda x, y: -5.0 + x / slope_denominator)
        domain.set_quantity('stage', 1.0)       # level surface: nothing to drive
        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        if with_culvert:
            anuga.Boyd_box_operator(
                domain, end_points=[[60.0, 25.0], [140.0, 25.0]],
                # Enquiry points well clear of the inlet regions: an enquiry
                # cell just outside a wide inlet is strongly coupled to what the
                # operator writes, which makes the run amplify roundoff and
                # measures something other than well-balancedness.
                enquiry_points=[[40.0, 25.0], [160.0, 25.0]],
                losses=1.5, width=20.0, height=3.0, apron=5.0,
                use_momentum_jet=False, use_velocity_head=False,
                manning=0.013, verbose=False)

        for _ in domain.evolve(yieldstep=0.5, finaltime=1.0):
            pass

        stage = domain.quantities['stage'].centroid_values
        return float(numpy.abs(stage - 1.0).max())

    def test_lake_at_rest_on_a_sloping_bed(self):
        """The steeper the bed, the worse the old behaviour was; now: nothing."""
        for slope_denominator in (50.0, 200.0):
            deviation = self._lake_at_rest(slope_denominator)
            # Was 6.8e-02 (1/50) and 1.8e-02 (1/200) with the uniform-depth
            # write; both are now at roundoff level over this interval.
            self.assertLess(
                deviation, 1e-5,
                'a culvert passing no flow disturbed a lake at rest by %.3e m '
                'on a 1/%g bed' % (deviation, slope_denominator))

    def test_matches_a_domain_without_the_structure(self):
        """The structure should be as quiet as not having one at all."""
        self.assertLess(self._lake_at_rest(50.0, with_culvert=False), 1e-12)


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Structure_operator)
    suite.addTests(unittest.TestLoader().loadTestsFromTestCase(Test_Inlet_enquiry))
    suite.addTests(unittest.TestLoader().loadTestsFromTestCase(Test_level_stages_to_average))
    suite.addTests(unittest.TestLoader().loadTestsFromTestCase(Test_Structure_well_balanced))
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
