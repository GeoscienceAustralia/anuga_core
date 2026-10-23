"""  Test set_quantity
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

from anuga.operators.set_quantity import *
from anuga.operators.set_stage import *

import numpy as num
from pprint import pprint
import warnings
import time



class Test_set_quantity(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


    def test_set_quantity_simple(self):
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


        update_stage = Set_quantity(domain, "stage", value=stage, indices=indices, test_stage=False)

        # Apply Operator
        update_stage()

        stage_ex = [ 3.,  3.,   1.,  3.]


        #print domain.quantities['stage'].centroid_values
        #print domain.quantities['xmomentum'].centroid_values
        #print domain.quantities['ymomentum'].centroid_values

        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)


    def test_set_quantity_negative(self):
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
        domain.set_quantity('elevation', lambda x,y : -2*x)
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


        try:
            update_stage = Set_quantity(domain, 'stage', value=stage, indices=indices)
        except AssertionError:
            pass
        except Exception:
            self.fail('Unexpected exception thrown:', e)
        else:
            self.fail('ExpectedException not thrown')


        update_stage = Set_quantity(domain, 'stage', value=stage, indices=indices, test_stage=False)
        # Apply Operator
        domain.timestep = 2.0
        update_stage()

        stage_ex = [ -5.,  -5.,   1.,  -5.]

        #print domain.quantities['elevation'].centroid_values
        #print domain.quantities['stage'].centroid_values
        #print domain.quantities['xmomentum'].centroid_values
        #print domain.quantities['ymomentum'].centroid_values

        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)

        update_stage = Set_stage(domain, stage=stage, indices=indices)
        domain.timestep = 2.0
        update_stage()

        stage_ex = [-1.33333333, -2.66666667,  1.,         -1.33333333]


        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)


    def test_set_quantity_function(self):
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

        operator = Set_quantity(domain, 'stage', value=stage, indices=indices, test_stage=False)

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


    def test_set_quantity_large_mesh_function(self):
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
        operator = operator = Set_quantity(domain, 'stage', value=stage, polygon=polygon, test_stage=False)



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

        operator = Set_quantity(domain, 'stage', value=stage, line=line, test_stage=False)


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



class Test_set_quantity_extra(unittest.TestCase):
    """Cover previously uncovered paths in set_quantity.py."""

    def _make_domain(self):
        points = [[0.0, 0.0], [0.0, 2.0], [2.0, 0.0], [0.0, 4.0], [2.0, 2.0], [4.0, 0.0]]
        vertices = [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]]
        domain = Domain(points, vertices)
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)
        domain.set_boundary({'exterior': Reflective_boundary(domain)})
        return domain

    def test_update_valueerror_no_indices(self):
        """Lines 112-117: ValueError during update with indices=None is caught."""
        domain = self._make_domain()

        def bad_value(x, y):
            raise ValueError('out of range')

        sq = Set_quantity(domain, 'stage', value=bad_value, test_stage=False)
        # Should not raise — ValueError is caught and logged
        sq()

    def test_update_valueerror_with_indices(self):
        """Lines 131-132: ValueError during update with indices set is caught."""
        domain = self._make_domain()

        def bad_value(x, y):
            raise ValueError('out of range')

        sq = Set_quantity(domain, 'stage', value=bad_value, indices=[0, 1], test_stage=False)
        # Should not raise — ValueError is caught and logged
        sq()

    def test_get_value_explicit_t(self):
        """Lines 156->160: get_value with explicit t skips domain.get_time()."""
        domain = self._make_domain()
        sq = Set_quantity(domain, 'stage', value=lambda t: t * 2.0, test_stage=False)
        result = sq.get_value(t=5.0)
        self.assertAlmostEqual(result, 10.0)

    def test_get_value_xy_type(self):
        """Line 163: get_value with value_type == 'x,y'."""
        domain = self._make_domain()
        import numpy as np
        sq = Set_quantity(domain, 'stage', value=lambda x, y: x + y, test_stage=False)
        x = np.array([1.0, 2.0])
        y = np.array([3.0, 4.0])
        result = sq.get_value(x=x, y=y, t=0.0)
        self.assertTrue(num.allclose(result, [4.0, 6.0]))


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_set_quantity)
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
