"""  Test set operators - w_uh_vh elevation erosion.
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

from anuga.operators.set_w_uh_vh_operator import *

import numpy as num
import warnings
import time



class Test_set_w_uh_vh_operators(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


    def test_set_w_uh_vh_operator_simple(self):
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

        #Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('xmomentum', 7.0)
        domain.set_quantity('ymomentum', 8.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})


#        print domain.quantities['w_uh_vh'].centroid_values
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        # Apply operator to these triangles
        indices = [0,1,3]

        w_uh_vh = [3.0, 4.0, 5.0]


        operator = Set_w_uh_vh_operator(domain, w_uh_vh=w_uh_vh, indices=indices)

        # Apply Operator
        domain.timestep = 2.0
        operator()

        stage_ex = [ 3.,  3.,   1.,  3.]
        xmom_ex = [ 4.,  4.,   7.,  4.]
        ymom_ex = [ 5.,  5.,   8.,  5.]


        #print domain.quantities['stage'].centroid_values
        #print domain.quantities['xmomentum'].centroid_values
        #print domain.quantities['ymomentum'].centroid_values

        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, xmom_ex)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, ymom_ex)

    def test_set_w_uh_vh_operator_simple_time(self):
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

        #Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('xmomentum', 7.0)
        domain.set_quantity('ymomentum', 8.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})


#        print domain.quantities['w_uh_vh'].centroid_values
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        # Apply operator to these triangles
        indices = [0,1,3]

        w_uh_vh = lambda t : [3.0, 4.0, 5.0]


        operator = Set_w_uh_vh_operator(domain, w_uh_vh=w_uh_vh, indices=indices)

        # Apply Operator
        domain.timestep = 2.0
        operator()

        stage_ex = [ 3.,  3.,   1.,  3.]
        xmom_ex = [ 4.,  4.,   7.,  4.]
        ymom_ex = [ 5.,  5.,   8.,  5.]


        #print domain.quantities['stage'].centroid_values
        #print domain.quantities['xmomentum'].centroid_values
        #print domain.quantities['ymomentum'].centroid_values

        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, xmom_ex)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, ymom_ex)


    def test_set_w_uh_vh_operator_time(self):
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

        #Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('xmomentum', 7.0)
        domain.set_quantity('ymomentum', 8.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})


#        print domain.quantities['w_uh_vh'].centroid_values
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        # Apply operator to these triangles
        indices = [0,1,3]

        w_uh_vh = lambda t : [t, t+1, t+2]


        operator = Set_w_uh_vh_operator(domain, w_uh_vh=w_uh_vh, indices=indices)

        # Apply Operator
        domain.timestep = 2.0
        domain.set_time(1.0)
        operator()

        t = domain.get_time()
        stage_ex = [ t,  t,   1.,  t]
        xmom_ex = [ t+1,  t+1,   7.,  t+1]
        ymom_ex = [ t+2,  t+2,   8.,  t+2]


        #print domain.quantities['stage'].centroid_values
        #print domain.quantities['xmomentum'].centroid_values
        #print domain.quantities['ymomentum'].centroid_values

        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, xmom_ex)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, ymom_ex)




class Test_set_w_uh_vh_operator_extra(unittest.TestCase):
    """Tests for uncovered methods in Set_w_uh_vh_operator."""

    def setUp(self):
        self.domain = rectangular_cross_domain(2, 2)
        self.domain.set_quantity('elevation', 0.0)
        self.domain.set_quantity('stage', 1.0)
        Br = Reflective_boundary(self.domain)
        self.domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

    def tearDown(self):
        try:
            import os
            os.remove('domain.sww')
        except OSError:
            pass

    def _make_op(self, w_uh_vh=None, indices=None):
        from anuga.operators.set_w_uh_vh_operator import Set_w_uh_vh_operator
        return Set_w_uh_vh_operator(
            self.domain, w_uh_vh=w_uh_vh, indices=indices)

    def test_parallel_safe(self):
        op = self._make_op(w_uh_vh=[1.0, 0.0, 0.0])
        self.assertTrue(op.parallel_safe())

    def test_statistics(self):
        op = self._make_op(w_uh_vh=[1.0, 0.0, 0.0])
        msg = op.statistics()
        self.assertIsInstance(msg, str)

    def test_timestepping_statistics(self):
        op = self._make_op(w_uh_vh=[1.0, 0.0, 0.0])
        msg = op.timestepping_statistics()
        self.assertIsInstance(msg, str)

    def test_call_none_w_uh_vh(self):
        """w_uh_vh=None triggers early return (line 85)."""
        op = self._make_op(w_uh_vh=None)
        self.domain.timestep = 1.0
        op()  # should return early without raising

    def test_call_all_triangles(self):
        """Call with indices=None updates all stage/xmom/ymom (lines 92-94)."""
        op = self._make_op(w_uh_vh=[1.5, 0.1, 0.2], indices=None)
        self.domain.timestep = 1.0
        op()
        import numpy as num
        self.assertTrue(num.allclose(
            self.domain.quantities['stage'].centroid_values, 1.5))

    def test_call_empty_indices(self):
        """Empty indices list → early return (line 80)."""
        op = self._make_op(w_uh_vh=[1.5, 0.1, 0.2], indices=[])
        self.domain.timestep = 1.0
        op()  # should return without updating

    def test_set_w_uh_vh_scalar(self):
        """set_w_uh_vh with scalar value (line 109)."""
        op = self._make_op(w_uh_vh=[1.0, 0.0, 0.0])
        op.set_w_uh_vh(2.0)  # scalar type → float cast


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_set_w_uh_vh_operators)
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
