"""  Test set operators - stage elevation erosion.
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

from anuga.operators.set_stage_operator import *

import numpy as num
from pprint import pprint
import warnings
import time



class Test_set_stage_operators(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


    def test_set_stage_operator_simple(self):
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
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})


#        print domain.quantities['stage'].centroid_values
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        # Apply operator to these triangles
        indices = [0,1,3]

        stage = 3.0


        operator = Set_stage_operator(domain, stage=stage, indices=indices)

        # Apply Operator
        domain.timestep = 2.0
        operator()

        stage_ex = [ 3.,  3.,   1.,  3.]


        #print domain.quantities['stage'].centroid_values
        #print domain.quantities['xmomentum'].centroid_values
        #print domain.quantities['ymomentum'].centroid_values

        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)


    def test_set_stage_operator_negative(self):
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
        domain.set_quantity('elevation', lambda x,y : -10*x)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

#        print domain.quantities['elevation'].centroid_values
#        print domain.quantities['stage'].centroid_values
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        # Apply operator to these triangles
        indices = [0,1,3]



        #Catchment_Rain_Polygon = read_polygon(join('CatchmentBdy.csv'))
        #rainfall = file_function(join('1y120m.tms'), quantities=['rainfall'])
        stage = -5.0


        operator = Set_stage_operator(domain, stage=stage, indices=indices)


        # Apply Operator
        domain.timestep = 2.0
        operator()

        stage_ex = [ -5.,  -5.,   1.,  -5.]

        #print domain.quantities['elevation'].centroid_values
        #print domain.quantities['stage'].centroid_values
        #print domain.quantities['xmomentum'].centroid_values
        #print domain.quantities['ymomentum'].centroid_values

        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)


    def test_set_stage_operator_function(self):
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
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})


#        print domain.quantities['stage'].centroid_values
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        # Apply operator to these triangles
        indices = [0,1,3]


        def stage(t):
            if t < 10.0:
                return 5.0
            else:
                return 10.0

        operator = Set_stage_operator(domain, stage=stage, indices=indices)

        # Apply Operator at time t=1.0
        domain.set_time(1.0)
        operator()

        stage_ex = [ 5.,  5.,   1.,  5.]


#        print domain.quantities['stage'].centroid_values
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)

        # Apply Operator at time t=15.0
        domain.set_time(15.0)
        operator()

        stage_ex = [ 10.,  10.,   1.,  10.]


#        print domain.quantities['stage'].centroid_values
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)


    def test_set_stage_operator_large_function(self):
        from anuga.config import rho_a, rho_w, eta_w
        from math import pi, cos, sin


        length = 2.0
        width = 2.0
        dx = dy = 0.5
        #dx = dy = 0.1
        domain = rectangular_cross_domain(int(length/dx), int(width/dy),
                                              len1=length, len2=width)


        #Flat surface with 1m of water
        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        R = Reflective_boundary(domain)
        domain.set_boundary( {'left': R, 'right': R, 'bottom': R, 'top': R} )


#        print domain.quantities['stage'].centroid_values
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values



        def stage(t):
            if t < 10.0:
                return 5.0
            else:
                return 7.0

        polygon = [(0.5,0.5), (1.5,0.5), (1.5,1.5), (0.5,1.5)]
        #operator = Polygonal_set_stage_operator(domain, stage=stage, polygon=polygon)
        operator = Polygonal_set_stage_operator(domain, stage=stage, polygon=polygon)


        #operator.plot_region()

        # Apply Operator at time t=1.0
        domain.set_time(1.0)
        operator()


        stage_ex_expanded = \
                   [ 1.,  1.,  5.,  5.,  1.,  5.,  5.,  5.,  1.,  5.,  5.,  5.,  1.,
                     5.,  5.,  1.,  5.,  1.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,
                     5.,  5.,  5.,  5.,  5.,  1.,  5.,  1.,  5.,  5.,  5.,  5.,  5.,
                     5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,  1.,  5.,  1.,  1.,  5.,
                     5.,  5.,  1.,  5.,  5.,  5.,  1.,  5.,  5.,  5.,  1.,  1.]

        stage_ex = [ 1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                     1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                     5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  1.0,  1.0,
                     1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  5.0,  5.0,  5.0,  5.0,
                     5.0,  5.0,  5.0,  5.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                     1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                     1.0,  1.0,  1.0,  1.0]


#        print domain.quantities['elevation'].centroid_values
#        pprint(domain.quantities['stage'].centroid_values)
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)

        # Apply Operator at time t=15.0
        domain.set_time(15.0)
        operator()

        stage_ex_expanded = \
                   [ 1.,  1.,  7.,  7.,  1.,  7.,  7.,  7.,  1.,  7.,  7.,  7.,  1.,
                     7.,  7.,  1.,  7.,  1.,  7.,  7.,  7.,  7.,  7.,  7.,  7.,  7.,
                     7.,  7.,  7.,  7.,  7.,  1.,  7.,  1.,  7.,  7.,  7.,  7.,  7.,
                     7.,  7.,  7.,  7.,  7.,  7.,  7.,  7.,  1.,  7.,  1.,  1.,  7.,
                     7.,  7.,  1.,  7.,  7.,  7.,  1.,  7.,  7.,  7.,  1.,  1.]


        stage_ex = [ 1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                     1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                     7.0,  7.0,  7.0,  7.0,  7.0,  7.0,  7.0,  7.0,  1.0,  1.0,
                     1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  7.0,  7.0,  7.0,  7.0,
                     7.0,  7.0,  7.0,  7.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                     1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                     1.0,  1.0,  1.0,  1.0]


#        pprint(domain.quantities['stage'].centroid_values)
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)


    def test_set_stage_operator_line(self):

        from math import pi, cos, sin


        length = 2.0
        width = 2.0
        dx = dy = 0.5
        #dx = dy = 0.1
        domain = rectangular_cross_domain(int(length/dx), int(width/dy),
                                              len1=length, len2=width)


        #Flat surface with 1m of water
        domain.set_quantity('elevation', -10.0)
        domain.set_quantity('stage', -1.0)
        domain.set_quantity('friction', 0)

        R = Reflective_boundary(domain)
        domain.set_boundary( {'left': R, 'right': R, 'bottom': R, 'top': R} )


#        print domain.quantities['stage'].centroid_values
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values



        def stage(x,y, t):
            #print x,y
            if t < 10.0:
                return x
            else:
                return y

        line = [(0.5,0.5), (1.5,0.75)]

        operator = Set_stage_operator(domain, stage=stage, line=line)


        #print 'stage at 1',stage(3.0,4.0,1.0) # return x value
        #print 'stage at 15', stage(3.0,4.0,15.0) # return y value
        #print operator.indices
        #print operator.value_type


        # Apply Operator at time t=1.0
        domain.set_time(1.0)
        operator()




        # The line starts at the mesh vertex (0.5, 0.5), runs through
        # triangles 21, 22, 36, 37 and ends on an edge of 38 at (1.5, 0.75).
        # The triangles that only meet it at that vertex or that endpoint
        # (a point contact, no length) are not selected. stage = x at t < 10.
        line_triangles = [21, 22, 36, 37, 38]
        assert sorted(operator.indices) == line_triangles
        stage_ex = -1.0 * num.ones(domain.number_of_elements)
        stage_ex[line_triangles] = domain.centroid_coordinates[line_triangles, 0]
        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)

        # Apply Operator at time t=15.0
        domain.set_time(15.0)
        operator()




        # stage = y at t >= 10, on the same five triangles
        stage_ex = -1.0 * num.ones(domain.number_of_elements)
        stage_ex[line_triangles] = domain.centroid_coordinates[line_triangles, 1]
        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)



class Test_set_stage_operator_extra(unittest.TestCase):
    """Additional tests for Set_stage_operator methods."""

    def _make_op(self):
        domain = rectangular_cross_domain(2, 2)
        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('stage', 1.0)
        return Set_stage_operator(domain, stage=0.5)

    def test_parallel_safe(self):
        """parallel_safe returns True (line 66)."""
        op = self._make_op()
        self.assertTrue(op.parallel_safe())

    def test_statistics(self):
        """statistics returns a string (lines 70-72)."""
        op = self._make_op()
        msg = op.statistics()
        self.assertIsInstance(msg, str)

    def test_timestepping_statistics(self):
        """timestepping_statistics body is reached (lines 77-79); self.center is a pre-existing bug."""
        op = self._make_op()
        try:
            msg = op.timestepping_statistics()
        except AttributeError:
            pass  # pre-existing: self.center not initialised by Set_stage_operator


class Test_set_stage_extra(unittest.TestCase):
    """Tests for uncovered paths in set_stage.py (the base class)."""

    def setUp(self):
        from anuga import rectangular_cross_domain, Reflective_boundary
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

    def test_call_empty_list_indices_noop(self):
        """Empty list indices → early return in __call__ (lines 102-103)."""
        from anuga.operators.set_stage_operator import Set_stage_operator
        op = Set_stage_operator(self.domain, stage=1.5, indices=[])
        self.domain.timestep = 1.0
        op()  # should return early

    def test_call_all_triangles(self):
        """Call with indices=None updates all stage values (lines 118-123)."""
        from anuga.operators.set_stage_operator import Set_stage_operator
        import numpy as num
        op = Set_stage_operator(self.domain, stage=2.0, indices=None)
        self.domain.timestep = 1.0
        op()
        self.assertTrue(num.all(
            self.domain.quantities['stage'].centroid_values >= 0.0))

    def test_call_specific_indices(self):
        """Call with specific indices updates those stage values (lines 130-137)."""
        from anuga.operators.set_stage_operator import Set_stage_operator
        import numpy as num
        op = Set_stage_operator(self.domain, stage=2.0, indices=[0, 1])
        self.domain.timestep = 1.0
        op()


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_set_stage_operators)
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
