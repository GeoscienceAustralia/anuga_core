import unittest
import os

import numpy
from anuga.structures import internal_boundary_functions
from anuga.structures.internal_boundary_functions import hecras_internal_boundary_function
from anuga.utilities.system_tools import get_pathname_from_package


class Test_internal_boundary_functions(unittest.TestCase):

    def setUp(self):

        path = get_pathname_from_package('anuga.structures')
        self.input_file = filename1 = os.path.join(path, 'tests', 'data', 'hecras_bridge_table.csv')
        self.hb = hecras_internal_boundary_function(self.input_file, verbose=False)

        return

    def tearDown(self):
        pass

        return

    def test_basic_properties(self):

        assert self.hb.name == self.input_file
        assert len(self.hb.hw_max_given_tw) == len(self.hb.nonfree_flow_tw)
        assert len(self.hb.nonfree_flow_curves) == len(self.hb.nonfree_flow_tw)

        return

    def test_function_values(self):

        # Stationary states #
        assert numpy.allclose(self.hb(-3., -3.), 0.)
        assert numpy.allclose(self.hb(-2.5, -2.5), 0.)
        assert numpy.allclose(self.hb(0., 0.), 0.)
        assert numpy.allclose(self.hb(1., 1.), 0.)
        assert numpy.allclose(self.hb(2.8, 2.8), 0.)
        assert numpy.allclose(self.hb(-2.146, -2.146), 0.)

        # Known values (looked up manually in table) #

        # Top of one nonfree-flow curve
        assert numpy.allclose(self.hb(-2.759, -2.947), 9.425)

        # Slightly above previous (i.e. free flow curve)
        assert numpy.allclose(self.hb(-2.747, -2.947), 9.736)

        # Value on most extreme nonfree flow curve
        assert numpy.allclose(self.hb(2.909, 2.894), 5.118)

        # Value on least-extreme nonfree flow curve
        assert numpy.allclose(self.hb(-3.26, -3.266), 0.27)

        # Value below least-extreme nonfree flow curve, but hw still valid on
        # free flow curve
        assert numpy.allclose(self.hb(-2.747, -3.4), 9.736)

        # Near top of free flow curve
        assert numpy.allclose(self.hb(2.468, -1.4), 82.89)
        assert numpy.allclose(self.hb(2.468, -2.0), 82.89)
        assert numpy.allclose(self.hb(2.468, -5.0), 82.89)

        # Logical ordering
        assert self.hb(2.4, 2.0) > self.hb(2.2, 2.0)
        assert self.hb(1.5, 1.0) > self.hb(1.5, 1.01)

        # Things which should throw errors
        def tw_too_big():
            return self.hb(9.00, 2.0)

        def hw_too_big():
            return self.hb(10.00, 0.)

        def hw_too_small():
            return self.hb(-4.0, -4.)

        self.assertRaises(Exception, lambda: tw_too_big())
        self.assertRaises(Exception, lambda: hw_too_big())
        self.assertRaises(Exception, lambda: hw_too_small())


        # Check case where sign reversal is allowed
        self.hb.allow_sign_reversal=True
        Q1 = self.hb( -2.75, -2.95)
        Q2 = self.hb( -2.95, -2.75)
        assert numpy.allclose(Q1, -Q2)

        return


if __name__ == "__main__":
    suite = unittest.makeSuite(Test_internal_boundary_functions, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
