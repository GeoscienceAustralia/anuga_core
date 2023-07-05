"""  Test erosion operators
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

from anuga.operators.erosion_operators import Erosion_operator

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

            
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_erosion_operators, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
