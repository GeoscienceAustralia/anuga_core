"""  Test set operators - stage elevation erosion.
"""

import unittest, os
import anuga
from anuga import Domain
from anuga import Reflective_boundary
from anuga import rectangular_cross_domain
from anuga import file_function
from anuga import Region


from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a
from anuga.file_conversion.file_conversion import timefile2netcdf
from anuga.config import time_format

from anuga.operators.set_friction_operators import *

import numpy as np
from pprint import pprint
import warnings
import time



class Test_set_friction_operators(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


    def test_set_friction_operator_float(self):
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

        region = Region(domain, indices=indices)

        friction_float = 3.0


        operator = Set_depth_friction_operator(domain, friction=friction_float, region=region)
        
        # Apply Operator
        domain.timestep = 2.0
        operator()



        #print domain.quantities['stage'].centroid_values
        #print domain.quantities['xmomentum'].centroid_values
        #print domain.quantities['ymomentum'].centroid_values

        assert np.allclose(domain.quantities['friction'].centroid_values[indices], friction_float)
        assert np.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert np.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)

 


    def test_set_friction_operator_function(self):
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


        def friction(h):
            if h < 10.0:
                return 5.0
            else:
                return 10.0

        operator = Set_depth_friction_operator(domain, friction=friction, indices=indices)

        # Apply Operator at time t=1.0
        domain.set_time(1.0)
        operator()

        friction_ex = [ 5.,  5.,   0.0,  5.]



        assert np.allclose(domain.quantities['friction'].centroid_values, friction_ex)
        assert np.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert np.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)

        # Apply Operator at time t=15.0
        domain.set_quantity('stage', 15.0)
        operator()

        friction_ex = [ 10.,  10.,   0.,  10.]



        assert np.allclose(domain.quantities['friction'].centroid_values, friction_ex)
        assert np.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert np.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)


    def test_friction_none(self):
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        vertices = [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]]

        domain = Domain(points, vertices)

        # Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        # Apply operator to these triangles
        indices = [0, 1, 3]

        region = Region(domain, indices=indices)

        # Test that friction=None raises ValueError
        with self.assertRaises(ValueError):
            Set_depth_friction_operator(domain, friction=None, region=region)

    def test_friction_string(self):
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        vertices = [[1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4]]

        domain = Domain(points, vertices)

        # Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        # Apply operator to these triangles
        indices = [0, 1, 3]

        region = Region(domain, indices=indices)

        # Test that friction as a string raises ValueError
        with self.assertRaises(ValueError):
            Set_depth_friction_operator(domain, friction="invalid_string", region=region)

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_set_friction_operators, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
