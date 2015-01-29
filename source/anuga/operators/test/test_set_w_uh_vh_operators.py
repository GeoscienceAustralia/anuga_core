"""  Test set operators - w_uh_vh elevation erosion.
"""

import unittest, os
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
        domain.time = 1.0
        operator()

        t = domain.time
        stage_ex = [ t,  t,   1.,  t]
        xmom_ex = [ t+1,  t+1,   7.,  t+1]
        ymom_ex = [ t+2,  t+2,   8.,  t+2]


        #print domain.quantities['stage'].centroid_values
        #print domain.quantities['xmomentum'].centroid_values
        #print domain.quantities['ymomentum'].centroid_values

        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, xmom_ex)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, ymom_ex)




if __name__ == "__main__":
    suite = unittest.makeSuite(Test_set_w_uh_vh_operators, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
