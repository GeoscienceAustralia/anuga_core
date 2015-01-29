"""  Test set operators - stage elevation erosion.
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

from anuga.operators.set_elevation_operator import *

import numpy as num
import warnings
import time



class Test_set_elevation_operator(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


    def test_set_elevation_operator_simple(self):
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

        stage_c = domain.quantities['stage'].centroid_values
        elev_c = domain.quantities['elevation'].centroid_values

        height_c = stage_c - elev_c

        integral0 = num.sum(height_c)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})


#        print domain.quantities['stage'].centroid_values
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        # Apply operator to these triangles
        indices = [0,1,3]

        elev = 3.0


        operator = Set_elevation_operator(domain, elevation=elev, indices=indices)
        
        # Apply Operator
        domain.timestep = 2.0
        operator()


        height_c = stage_c - elev_c

        integral1 = num.sum(height_c)

        assert integral0 == integral1

        stage_ex = [ 3.66666667,  3.33333333,  2.33333333,  3.66666667]

        elev_ex = [ 2.66666667,  2.33333333,  1.33333333,  2.66666667]
        

#        print domain.quantities['elevation'].centroid_values
#        print domain.quantities['stage'].centroid_values
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        assert num.allclose(domain.quantities['elevation'].centroid_values, elev_ex)
        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)

 
    def test_set_elevation_operator_negative(self):
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


        stage_c = domain.quantities['stage'].centroid_values
        elev_c = domain.quantities['elevation'].centroid_values

        height_c = stage_c - elev_c

        integral0 = num.sum(height_c)


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
        elev = -5.0


        operator = Set_elevation_operator(domain, elevation=elev, indices=indices)


        # Apply Operator
        domain.timestep = 2.0
        operator()
        
        height_c = stage_c - elev_c
        integral1 = num.sum(height_c)
        assert integral0 == integral1

        elev_ex  = [-4.88888889, -4.77777778, -5.77777778, -4.88888889]
        stage_ex = [-2.55555556, -1.11111111,  0.55555556, -2.55555556]


#        print domain.quantities['elevation'].centroid_values
#        print domain.quantities['stage'].centroid_values
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        assert num.allclose(domain.quantities['elevation'].centroid_values, elev_ex)
        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)


    def test_set_elevation_operator_small_function(self):
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
        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})


#        print domain.quantities['stage'].centroid_values
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        # Apply operator to these triangles
        indices = [0,1,3]


        def elev(t):
            if t < 10.0:
                return 5.0
            else:
                return 7.0

        operator = Set_elevation_operator(domain, elevation=elev, indices=indices)

        # Apply Operator at time t=1.0
        domain.set_time(1.0)
        operator()


        elev_ex  = [ 4.44444444,  3.88888889, 2.22222222,  4.44444444]
        stage_ex = [ 5.44444444,  4.88888889,  3.22222222,  5.44444444]

#        print domain.quantities['elevation'].centroid_values
#        print domain.quantities['stage'].centroid_values
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        assert num.allclose(domain.quantities['elevation'].centroid_values, elev_ex)
        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)

        # Apply Operator at time t=15.0
        domain.set_time(15.0)
        operator()


        elev_ex =  [ 6.59259259,  6.18518519,  3.85185185,  6.59259259]
        stage_ex = [ 7.59259259,  7.18518519,  4.85185185,  7.59259259]

#        print domain.quantities['elevation'].centroid_values
#        print domain.quantities['stage'].centroid_values
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        assert num.allclose(domain.quantities['elevation'].centroid_values, elev_ex)
        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)



    def test_set_polygonal_elevation_operator_large_function(self):
        from anuga.config import rho_a, rho_w, eta_w
        from math import pi, cos, sin


        length = 2.0
        width = 2.0
        dx = dy = 0.5
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

        # Apply operator to these triangles
        indices = [0,1,3]


        def elev(t):
            if t < 10.0:
                return 5.0
            else:
                return 7.0

        polygon = [(0.5,0.5), (1.5,0.5), (1.5,1.5), (0.5,1.5)]
        operator = Polygonal_set_elevation_operator(domain, elevation=elev, polygon=polygon)

        # Apply Operator at time t=1.0
        domain.set_time(1.0)
        operator()




        elev_ex = [ 0.        ,  0.        ,  0.41666667,  0.41666667,  0.        ,
        0.41666667,  1.25      ,  0.83333333,  0.        ,  0.83333333,
        1.25      ,  0.41666667,  0.        ,  0.41666667,  0.41666667,
        0.        ,  0.41666667,  0.        ,  0.83333333,  1.25      ,
        2.91666667,  2.91666667,  4.16666667,  4.16666667,  2.91666667,
        4.16666667,  4.16666667,  2.91666667,  0.41666667,  1.25      ,
        0.83333333,  0.        ,  0.83333333,  0.        ,  0.41666667,
        1.25      ,  4.16666667,  2.91666667,  2.91666667,  4.16666667,
        4.16666667,  4.16666667,  2.91666667,  2.91666667,  0.83333333,
        1.25      ,  0.41666667,  0.        ,  0.41666667,  0.        ,
        0.        ,  0.41666667,  1.25      ,  0.41666667,  0.        ,
        0.83333333,  1.25      ,  0.83333333,  0.        ,  0.41666667,
        0.41666667,  0.41666667,  0.        ,  0.        ]

        stage_ex = [ 1.        ,  1.        ,  1.41666667,  1.41666667,  1.        ,
        1.41666667,  2.25      ,  1.83333333,  1.        ,  1.83333333,
        2.25      ,  1.41666667,  1.        ,  1.41666667,  1.41666667,
        1.        ,  1.41666667,  1.        ,  1.83333333,  2.25      ,
        3.91666667,  3.91666667,  5.16666667,  5.16666667,  3.91666667,
        5.16666667,  5.16666667,  3.91666667,  1.41666667,  2.25      ,
        1.83333333,  1.        ,  1.83333333,  1.        ,  1.41666667,
        2.25      ,  5.16666667,  3.91666667,  3.91666667,  5.16666667,
        5.16666667,  5.16666667,  3.91666667,  3.91666667,  1.83333333,
        2.25      ,  1.41666667,  1.        ,  1.41666667,  1.        ,
        1.        ,  1.41666667,  2.25      ,  1.41666667,  1.        ,
        1.83333333,  2.25      ,  1.83333333,  1.        ,  1.41666667,
        1.41666667,  1.41666667,  1.        ,  1.        ]





#        from pprint import pprint
#        pprint (domain.quantities['elevation'].centroid_values)
#        pprint (domain.quantities['stage'].centroid_values)
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        assert num.allclose(domain.quantities['elevation'].centroid_values, elev_ex)
        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)

        # Apply Operator at time t=15.0
        domain.set_time(15.0)
        operator()




        elev_ex = [ 0.        ,  0.        ,  0.89583333,  0.89583333,  0.        ,
        0.89583333,  2.47916667,  1.58333333,  0.        ,  1.58333333,
        2.47916667,  0.89583333,  0.        ,  0.89583333,  0.89583333,
        0.        ,  0.89583333,  0.        ,  1.58333333,  2.47916667,
        4.8125    ,  4.8125    ,  6.25      ,  6.25      ,  4.8125    ,
        6.25      ,  6.25      ,  4.8125    ,  0.89583333,  2.47916667,
        1.58333333,  0.        ,  1.58333333,  0.        ,  0.89583333,
        2.47916667,  6.25      ,  4.8125    ,  4.8125    ,  6.25      ,
        6.25      ,  6.25      ,  4.8125    ,  4.8125    ,  1.58333333,
        2.47916667,  0.89583333,  0.        ,  0.89583333,  0.        ,
        0.        ,  0.89583333,  2.47916667,  0.89583333,  0.        ,
        1.58333333,  2.47916667,  1.58333333,  0.        ,  0.89583333,
        0.89583333,  0.89583333,  0.        ,  0.        ]

        stage_ex = [ 1.        ,  1.        ,  1.89583333,  1.89583333,  1.        ,
        1.89583333,  3.47916667,  2.58333333,  1.        ,  2.58333333,
        3.47916667,  1.89583333,  1.        ,  1.89583333,  1.89583333,
        1.        ,  1.89583333,  1.        ,  2.58333333,  3.47916667,
        5.8125    ,  5.8125    ,  7.25      ,  7.25      ,  5.8125    ,
        7.25      ,  7.25      ,  5.8125    ,  1.89583333,  3.47916667,
        2.58333333,  1.        ,  2.58333333,  1.        ,  1.89583333,
        3.47916667,  7.25      ,  5.8125    ,  5.8125    ,  7.25      ,
        7.25      ,  7.25      ,  5.8125    ,  5.8125    ,  2.58333333,
        3.47916667,  1.89583333,  1.        ,  1.89583333,  1.        ,
        1.        ,  1.89583333,  3.47916667,  1.89583333,  1.        ,
        2.58333333,  3.47916667,  2.58333333,  1.        ,  1.89583333,
        1.89583333,  1.89583333,  1.        ,  1.        ]

#        pprint (domain.quantities['elevation'].centroid_values)
#        pprint (domain.quantities['stage'].centroid_values)
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        assert num.allclose(domain.quantities['elevation'].centroid_values, elev_ex)
        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)


    def test_set_elevation_operator_large_function(self):
        from anuga.config import rho_a, rho_w, eta_w
        from math import pi, cos, sin


        length = 2.0
        width = 2.0
        dx = dy = 0.5
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

        # Apply operator to these triangles
        indices = [0,1,3]


        def elev(t):
            if t < 10.0:
                return 5.0
            else:
                return 7.0

        polygon = [(0.5,0.5), (1.5,0.5), (1.5,1.5), (0.5,1.5)]
        operator = Set_elevation_operator(domain, elevation=elev, polygon=polygon)

        # Apply Operator at time t=1.0
        domain.set_time(1.0)
        operator()

        stage_ex = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]


        elev_ex = [ 0.        ,  0.        ,  0.41666667,  0.41666667,  0.        ,
        0.41666667,  1.25      ,  0.83333333,  0.        ,  0.83333333,
        1.25      ,  0.41666667,  0.        ,  0.41666667,  0.41666667,
        0.        ,  0.41666667,  0.        ,  0.83333333,  1.25      ,
        2.91666667,  2.91666667,  4.16666667,  4.16666667,  2.91666667,
        4.16666667,  4.16666667,  2.91666667,  0.41666667,  1.25      ,
        0.83333333,  0.        ,  0.83333333,  0.        ,  0.41666667,
        1.25      ,  4.16666667,  2.91666667,  2.91666667,  4.16666667,
        4.16666667,  4.16666667,  2.91666667,  2.91666667,  0.83333333,
        1.25      ,  0.41666667,  0.        ,  0.41666667,  0.        ,
        0.        ,  0.41666667,  1.25      ,  0.41666667,  0.        ,
        0.83333333,  1.25      ,  0.83333333,  0.        ,  0.41666667,
        0.41666667,  0.41666667,  0.        ,  0.        ]

        stage_ex = [ 1.        ,  1.        ,  1.41666667,  1.41666667,  1.        ,
        1.41666667,  2.25      ,  1.83333333,  1.        ,  1.83333333,
        2.25      ,  1.41666667,  1.        ,  1.41666667,  1.41666667,
        1.        ,  1.41666667,  1.        ,  1.83333333,  2.25      ,
        3.91666667,  3.91666667,  5.16666667,  5.16666667,  3.91666667,
        5.16666667,  5.16666667,  3.91666667,  1.41666667,  2.25      ,
        1.83333333,  1.        ,  1.83333333,  1.        ,  1.41666667,
        2.25      ,  5.16666667,  3.91666667,  3.91666667,  5.16666667,
        5.16666667,  5.16666667,  3.91666667,  3.91666667,  1.83333333,
        2.25      ,  1.41666667,  1.        ,  1.41666667,  1.        ,
        1.        ,  1.41666667,  2.25      ,  1.41666667,  1.        ,
        1.83333333,  2.25      ,  1.83333333,  1.        ,  1.41666667,
        1.41666667,  1.41666667,  1.        ,  1.        ]





#        from pprint import pprint
#        pprint (domain.quantities['elevation'].centroid_values)
#        pprint (domain.quantities['stage'].centroid_values)
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        assert num.allclose(domain.quantities['elevation'].centroid_values, elev_ex)
        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)

        # Apply Operator at time t=15.0
        domain.set_time(15.0)
        operator()

        stage_ex = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]


        elev_ex = [ 0.        ,  0.        ,  0.89583333,  0.89583333,  0.        ,
        0.89583333,  2.47916667,  1.58333333,  0.        ,  1.58333333,
        2.47916667,  0.89583333,  0.        ,  0.89583333,  0.89583333,
        0.        ,  0.89583333,  0.        ,  1.58333333,  2.47916667,
        4.8125    ,  4.8125    ,  6.25      ,  6.25      ,  4.8125    ,
        6.25      ,  6.25      ,  4.8125    ,  0.89583333,  2.47916667,
        1.58333333,  0.        ,  1.58333333,  0.        ,  0.89583333,
        2.47916667,  6.25      ,  4.8125    ,  4.8125    ,  6.25      ,
        6.25      ,  6.25      ,  4.8125    ,  4.8125    ,  1.58333333,
        2.47916667,  0.89583333,  0.        ,  0.89583333,  0.        ,
        0.        ,  0.89583333,  2.47916667,  0.89583333,  0.        ,
        1.58333333,  2.47916667,  1.58333333,  0.        ,  0.89583333,
        0.89583333,  0.89583333,  0.        ,  0.        ]

        stage_ex = [ 1.        ,  1.        ,  1.89583333,  1.89583333,  1.        ,
        1.89583333,  3.47916667,  2.58333333,  1.        ,  2.58333333,
        3.47916667,  1.89583333,  1.        ,  1.89583333,  1.89583333,
        1.        ,  1.89583333,  1.        ,  2.58333333,  3.47916667,
        5.8125    ,  5.8125    ,  7.25      ,  7.25      ,  5.8125    ,
        7.25      ,  7.25      ,  5.8125    ,  1.89583333,  3.47916667,
        2.58333333,  1.        ,  2.58333333,  1.        ,  1.89583333,
        3.47916667,  7.25      ,  5.8125    ,  5.8125    ,  7.25      ,
        7.25      ,  7.25      ,  5.8125    ,  5.8125    ,  2.58333333,
        3.47916667,  1.89583333,  1.        ,  1.89583333,  1.        ,
        1.        ,  1.89583333,  3.47916667,  1.89583333,  1.        ,
        2.58333333,  3.47916667,  2.58333333,  1.        ,  1.89583333,
        1.89583333,  1.89583333,  1.        ,  1.        ]

#        pprint (domain.quantities['elevation'].centroid_values)
#        pprint (domain.quantities['stage'].centroid_values)
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        assert num.allclose(domain.quantities['elevation'].centroid_values, elev_ex)
        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)


    def test_set_circular_elevation_operator_large_function(self):
        from anuga.config import rho_a, rho_w, eta_w
        from math import pi, cos, sin


        length = 2.0
        width = 2.0
        dx = dy = 0.5
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

        # Apply operator to these triangles


        def elev(t):
            if t < 10.0:
                return 5.0
            else:
                return 7.0

        operator = Circular_set_elevation_operator(domain, elevation=elev, center=[1.0,1.0], radius=1.0)

        # Apply Operator at time t=1.0
        domain.set_time(1.0)
        operator()




        elev_ex = [ 2.08333333,  2.08333333,  3.75      ,  3.75      ,  4.58333333,
        4.58333333,  5.        ,  5.        ,  4.58333333,  5.        ,
        5.        ,  4.58333333,  2.08333333,  3.75      ,  3.75      ,
        2.08333333,  4.58333333,  4.58333333,  5.        ,  5.        ,
        5.        ,  5.        ,  5.        ,  5.        ,  5.        ,
        5.        ,  5.        ,  5.        ,  4.58333333,  5.        ,
        5.        ,  4.58333333,  5.        ,  4.58333333,  4.58333333,
        5.        ,  5.        ,  5.        ,  5.        ,  5.        ,
        5.        ,  5.        ,  5.        ,  5.        ,  5.        ,
        5.        ,  4.58333333,  4.58333333,  3.75      ,  2.08333333,
        2.08333333,  3.75      ,  5.        ,  4.58333333,  4.58333333,
        5.        ,  5.        ,  5.        ,  4.58333333,  4.58333333,
        3.75      ,  3.75      ,  2.08333333,  2.08333333]




        stage_ex = [ 3.08333333,  3.08333333,  4.75      ,  4.75      ,  5.58333333,
        5.58333333,  6.        ,  6.        ,  5.58333333,  6.        ,
        6.        ,  5.58333333,  3.08333333,  4.75      ,  4.75      ,
        3.08333333,  5.58333333,  5.58333333,  6.        ,  6.        ,
        6.        ,  6.        ,  6.        ,  6.        ,  6.        ,
        6.        ,  6.        ,  6.        ,  5.58333333,  6.        ,
        6.        ,  5.58333333,  6.        ,  5.58333333,  5.58333333,
        6.        ,  6.        ,  6.        ,  6.        ,  6.        ,
        6.        ,  6.        ,  6.        ,  6.        ,  6.        ,
        6.        ,  5.58333333,  5.58333333,  4.75      ,  3.08333333,
        3.08333333,  4.75      ,  6.        ,  5.58333333,  5.58333333,
        6.        ,  6.        ,  6.        ,  5.58333333,  5.58333333,
        4.75      ,  4.75      ,  3.08333333,  3.08333333]





#        from pprint import pprint
#        pprint (domain.quantities['elevation'].centroid_values)
#        pprint (domain.quantities['stage'].centroid_values)
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        assert num.allclose(domain.quantities['elevation'].centroid_values, elev_ex)
        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)

        # Apply Operator at time t=15.0
        domain.set_time(15.0)
        operator()




        elev_ex = [ 3.64583333,  3.64583333,  5.97916667,  5.97916667,  6.72916667,
        6.72916667,  7.        ,  7.        ,  6.72916667,  7.        ,
        7.        ,  6.72916667,  3.64583333,  5.97916667,  5.97916667,
        3.64583333,  6.72916667,  6.72916667,  7.        ,  7.        ,
        7.        ,  7.        ,  7.        ,  7.        ,  7.        ,
        7.        ,  7.        ,  7.        ,  6.72916667,  7.        ,
        7.        ,  6.72916667,  7.        ,  6.72916667,  6.72916667,
        7.        ,  7.        ,  7.        ,  7.        ,  7.        ,
        7.        ,  7.        ,  7.        ,  7.        ,  7.        ,
        7.        ,  6.72916667,  6.72916667,  5.97916667,  3.64583333,
        3.64583333,  5.97916667,  7.        ,  6.72916667,  6.72916667,
        7.        ,  7.        ,  7.        ,  6.72916667,  6.72916667,
        5.97916667,  5.97916667,  3.64583333,  3.64583333]


        stage_ex = [ 4.64583333,  4.64583333,  6.97916667,  6.97916667,  7.72916667,
        7.72916667,  8.        ,  8.        ,  7.72916667,  8.        ,
        8.        ,  7.72916667,  4.64583333,  6.97916667,  6.97916667,
        4.64583333,  7.72916667,  7.72916667,  8.        ,  8.        ,
        8.        ,  8.        ,  8.        ,  8.        ,  8.        ,
        8.        ,  8.        ,  8.        ,  7.72916667,  8.        ,
        8.        ,  7.72916667,  8.        ,  7.72916667,  7.72916667,
        8.        ,  8.        ,  8.        ,  8.        ,  8.        ,
        8.        ,  8.        ,  8.        ,  8.        ,  8.        ,
        8.        ,  7.72916667,  7.72916667,  6.97916667,  4.64583333,
        4.64583333,  6.97916667,  8.        ,  7.72916667,  7.72916667,
        8.        ,  8.        ,  8.        ,  7.72916667,  7.72916667,
        6.97916667,  6.97916667,  4.64583333,  4.64583333]


#        from pprint import pprint
#        pprint (domain.quantities['elevation'].centroid_values)
#        pprint (domain.quantities['stage'].centroid_values)
#        pprint (domain.quantities['xmomentum'].centroid_values)
#        pprint (domain.quantities['ymomentum'].centroid_values)

        assert num.allclose(domain.quantities['elevation'].centroid_values, elev_ex)
        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)

    def test_set_elevation_operator_center_radius(self):
        from math import pi, cos, sin


        length = 2.0
        width = 2.0
        dx = dy = 0.5
        domain = rectangular_cross_domain(int(length/dx), int(width/dy),
                                              len1=length, len2=width)


        #Flat surface with 1m of water
        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        R = Reflective_boundary(domain)
        domain.set_boundary( {'left': R, 'right': R, 'bottom': R, 'top': R} )

        from pprint import pprint
        #pprint(domain.quantities['stage'].centroid_values)
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        # Apply operator to these triangles


        def elev(t):
            if t < 10.0:
                return 5.0
            else:
                return 7.0

        operator = Set_elevation_operator(domain, elevation=elev, center=[1.0,1.0], radius=1.0)

        # Apply Operator at time t=1.0
        domain.set_time(1.0)
        operator()


        #pprint(domain.quantities['elevation'].centroid_values)

        elev_ex = [ 2.08333333,  2.08333333,  3.75      ,  3.75      ,  4.58333333,
        4.58333333,  5.        ,  5.        ,  4.58333333,  5.        ,
        5.        ,  4.58333333,  2.08333333,  3.75      ,  3.75      ,
        2.08333333,  4.58333333,  4.58333333,  5.        ,  5.        ,
        5.        ,  5.        ,  5.        ,  5.        ,  5.        ,
        5.        ,  5.        ,  5.        ,  4.58333333,  5.        ,
        5.        ,  4.58333333,  5.        ,  4.58333333,  4.58333333,
        5.        ,  5.        ,  5.        ,  5.        ,  5.        ,
        5.        ,  5.        ,  5.        ,  5.        ,  5.        ,
        5.        ,  4.58333333,  4.58333333,  3.75      ,  2.08333333,
        2.08333333,  3.75      ,  5.        ,  4.58333333,  4.58333333,
        5.        ,  5.        ,  5.        ,  4.58333333,  4.58333333,
        3.75      ,  3.75      ,  2.08333333,  2.08333333]




        #pprint(domain.quantities['stage'].centroid_values)
        
        stage_ex = [ 3.08333333,  3.08333333,  4.75      ,  4.75      ,  5.58333333,
        5.58333333,  6.        ,  6.        ,  5.58333333,  6.        ,
        6.        ,  5.58333333,  3.08333333,  4.75      ,  4.75      ,
        3.08333333,  5.58333333,  5.58333333,  6.        ,  6.        ,
        6.        ,  6.        ,  6.        ,  6.        ,  6.        ,
        6.        ,  6.        ,  6.        ,  5.58333333,  6.        ,
        6.        ,  5.58333333,  6.        ,  5.58333333,  5.58333333,
        6.        ,  6.        ,  6.        ,  6.        ,  6.        ,
        6.        ,  6.        ,  6.        ,  6.        ,  6.        ,
        6.        ,  5.58333333,  5.58333333,  4.75      ,  3.08333333,
        3.08333333,  4.75      ,  6.        ,  5.58333333,  5.58333333,
        6.        ,  6.        ,  6.        ,  5.58333333,  5.58333333,
        4.75      ,  4.75      ,  3.08333333,  3.08333333]



#        from pprint import pprint
#        pprint (domain.quantities['elevation'].centroid_values)
#        pprint (domain.quantities['stage'].centroid_values)
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        assert num.allclose(domain.quantities['elevation'].centroid_values, elev_ex)
        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)

        # Apply Operator at time t=15.0
        domain.set_time(15.0)
        operator()




        elev_ex = [ 3.64583333,  3.64583333,  5.97916667,  5.97916667,  6.72916667,
        6.72916667,  7.        ,  7.        ,  6.72916667,  7.        ,
        7.        ,  6.72916667,  3.64583333,  5.97916667,  5.97916667,
        3.64583333,  6.72916667,  6.72916667,  7.        ,  7.        ,
        7.        ,  7.        ,  7.        ,  7.        ,  7.        ,
        7.        ,  7.        ,  7.        ,  6.72916667,  7.        ,
        7.        ,  6.72916667,  7.        ,  6.72916667,  6.72916667,
        7.        ,  7.        ,  7.        ,  7.        ,  7.        ,
        7.        ,  7.        ,  7.        ,  7.        ,  7.        ,
        7.        ,  6.72916667,  6.72916667,  5.97916667,  3.64583333,
        3.64583333,  5.97916667,  7.        ,  6.72916667,  6.72916667,
        7.        ,  7.        ,  7.        ,  6.72916667,  6.72916667,
        5.97916667,  5.97916667,  3.64583333,  3.64583333]


        stage_ex = [ 4.64583333,  4.64583333,  6.97916667,  6.97916667,  7.72916667,
        7.72916667,  8.        ,  8.        ,  7.72916667,  8.        ,
        8.        ,  7.72916667,  4.64583333,  6.97916667,  6.97916667,
        4.64583333,  7.72916667,  7.72916667,  8.        ,  8.        ,
        8.        ,  8.        ,  8.        ,  8.        ,  8.        ,
        8.        ,  8.        ,  8.        ,  7.72916667,  8.        ,
        8.        ,  7.72916667,  8.        ,  7.72916667,  7.72916667,
        8.        ,  8.        ,  8.        ,  8.        ,  8.        ,
        8.        ,  8.        ,  8.        ,  8.        ,  8.        ,
        8.        ,  7.72916667,  7.72916667,  6.97916667,  4.64583333,
        4.64583333,  6.97916667,  8.        ,  7.72916667,  7.72916667,
        8.        ,  8.        ,  8.        ,  7.72916667,  7.72916667,
        6.97916667,  6.97916667,  4.64583333,  4.64583333]


#        from pprint import pprint
#        pprint (domain.quantities['elevation'].centroid_values)
#        pprint (domain.quantities['stage'].centroid_values)
#        pprint (domain.quantities['xmomentum'].centroid_values)
#        pprint (domain.quantities['ymomentum'].centroid_values)

        assert num.allclose(domain.quantities['elevation'].centroid_values, elev_ex)
        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)

        operator = Set_elevation(domain, elevation=0.0)


        #print operator.value_type
        
        operator()

        #from pprint import pprint
        #pprint (domain.quantities['elevation'].centroid_values)
#        pprint (domain.quantities['stage'].centroid_values)
#        pprint (domain.quantities['xmomentum'].centroid_values)
#        pprint (domain.quantities['ymomentum'].centroid_values)

        assert num.allclose(domain.quantities['elevation'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['stage'].centroid_values, 1.0)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)


        operator = Set_elevation(domain, elevation=lambda t: t, indices = [0,1,3])

        operator()



        elev_ex = [ 11.25 ,  10.   ,   5.625,   6.875,   2.5  ,   3.125,   0.625,
         0.   ,   0.   ,   0.   ,   0.   ,   0.   ,   0.   ,   0.   ,
         0.   ,   0.   ,   1.875,   1.25 ,   0.   ,   0.625,   0.625,
         0.625,   0.   ,   0.   ,   0.   ,   0.   ,   0.   ,   0.   ,
         0.   ,   0.   ,   0.   ,   0.   ,   0.   ,   0.   ,   0.   ,
         0.   ,   0.   ,   0.   ,   0.   ,   0.   ,   0.   ,   0.   ,
         0.   ,   0.   ,   0.   ,   0.   ,   0.   ,   0.   ,   0.   ,
         0.   ,   0.   ,   0.   ,   0.   ,   0.   ,   0.   ,   0.   ,
         0.   ,   0.   ,   0.   ,   0.   ,   0.   ,   0.   ,   0.   ,   0.   ]

        stage_ex = [ 12.25 ,  11.   ,   6.625,   7.875,   3.5  ,   4.125,   1.625,
         1.   ,   1.   ,   1.   ,   1.   ,   1.   ,   1.   ,   1.   ,
         1.   ,   1.   ,   2.875,   2.25 ,   1.   ,   1.625,   1.625,
         1.625,   1.   ,   1.   ,   1.   ,   1.   ,   1.   ,   1.   ,
         1.   ,   1.   ,   1.   ,   1.   ,   1.   ,   1.   ,   1.   ,
         1.   ,   1.   ,   1.   ,   1.   ,   1.   ,   1.   ,   1.   ,
         1.   ,   1.   ,   1.   ,   1.   ,   1.   ,   1.   ,   1.   ,
         1.   ,   1.   ,   1.   ,   1.   ,   1.   ,   1.   ,   1.   ,
         1.   ,   1.   ,   1.   ,   1.   ,   1.   ,   1.   ,   1.   ,   1.   ]



#        from pprint import pprint
#        pprint (domain.quantities['elevation'].centroid_values)
#        pprint (domain.quantities['stage'].centroid_values)
#        pprint (domain.quantities['xmomentum'].centroid_values)
#        pprint (domain.quantities['ymomentum'].centroid_values)

        assert num.allclose(domain.quantities['elevation'].centroid_values, elev_ex)
        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)

    def test_set_elevation_operator_center_radius_de1(self):
        from math import pi, cos, sin


        length = 2.0
        width = 2.0
        dx = dy = 0.5
        domain = rectangular_cross_domain(int(length/dx), int(width/dy),
                                              len1=length, len2=width)


        #Flat surface with 1m of water
        domain.set_flow_algorithm('DE1')
        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        R = Reflective_boundary(domain)
        domain.set_boundary( {'left': R, 'right': R, 'bottom': R, 'top': R} )

        from pprint import pprint
        #pprint(domain.quantities['stage'].centroid_values)
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        # Apply operator to these triangles


        def elev(t):
            if t < 10.0:
                return 5.0
            else:
                return 7.0

        operator = Set_elevation_operator(domain, elevation=elev, center=[1.0,1.0], radius=1.0)

        # Apply Operator at time t=1.0
        domain.set_time(1.0)
        operator()


        #pprint(domain.quantities['elevation'].centroid_values)

        elev_ex = [ 0.,  0.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,  0.,
        5.,  5.,  0.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,
        5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,
        5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,  0.,  0.,  5.,
        5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,  0.,  0.]


        #pprint(domain.quantities['stage'].centroid_values)
        
        stage_ex = [ 1.,  1.,  6.,  6.,  6.,  6.,  6.,  6.,  6.,  6.,  6.,  6.,  1.,
        6.,  6.,  1.,  6.,  6.,  6.,  6.,  6.,  6.,  6.,  6.,  6.,  6.,
        6.,  6.,  6.,  6.,  6.,  6.,  6.,  6.,  6.,  6.,  6.,  6.,  6.,
        6.,  6.,  6.,  6.,  6.,  6.,  6.,  6.,  6.,  6.,  1.,  1.,  6.,
        6.,  6.,  6.,  6.,  6.,  6.,  6.,  6.,  6.,  6.,  1.,  1.]



#        from pprint import pprint
#        pprint (domain.quantities['elevation'].centroid_values)
#        pprint (domain.quantities['stage'].centroid_values)
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        assert num.allclose(domain.quantities['elevation'].centroid_values, elev_ex)
        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)

        # Apply Operator at time t=15.0
        domain.set_time(15.0)
        operator()


        #pprint(domain.quantities['elevation'].centroid_values)

        elev_ex = [ 0.,  0.,  7.,  7.,  7.,  7.,  7.,  7.,  7.,  7.,  7.,  7.,  0.,
        7.,  7.,  0.,  7.,  7.,  7.,  7.,  7.,  7.,  7.,  7.,  7.,  7.,
        7.,  7.,  7.,  7.,  7.,  7.,  7.,  7.,  7.,  7.,  7.,  7.,  7.,
        7.,  7.,  7.,  7.,  7.,  7.,  7.,  7.,  7.,  7.,  0.,  0.,  7.,
        7.,  7.,  7.,  7.,  7.,  7.,  7.,  7.,  7.,  7.,  0.,  0.]


        #pprint(domain.quantities['stage'].centroid_values)

        stage_ex = [ 1.,  1.,  8.,  8.,  8.,  8.,  8.,  8.,  8.,  8.,  8.,  8.,  1.,
        8.,  8.,  1.,  8.,  8.,  8.,  8.,  8.,  8.,  8.,  8.,  8.,  8.,
        8.,  8.,  8.,  8.,  8.,  8.,  8.,  8.,  8.,  8.,  8.,  8.,  8.,
        8.,  8.,  8.,  8.,  8.,  8.,  8.,  8.,  8.,  8.,  1.,  1.,  8.,
        8.,  8.,  8.,  8.,  8.,  8.,  8.,  8.,  8.,  8.,  1.,  1.]
        

#        from pprint import pprint
#        pprint (domain.quantities['elevation'].centroid_values)
#        pprint (domain.quantities['stage'].centroid_values)
#        pprint (domain.quantities['xmomentum'].centroid_values)
#        pprint (domain.quantities['ymomentum'].centroid_values)

        assert num.allclose(domain.quantities['elevation'].centroid_values, elev_ex)
        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)

        operator = Set_elevation(domain, elevation=0.0)


        #print operator.value_type
        
        operator()

        #from pprint import pprint
#        pprint (domain.quantities['elevation'].centroid_values)
#        pprint (domain.quantities['stage'].centroid_values)
#        pprint (domain.quantities['xmomentum'].centroid_values)
#        pprint (domain.quantities['ymomentum'].centroid_values)

        assert num.allclose(domain.quantities['elevation'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['stage'].centroid_values, 1.0)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)


        operator = Set_elevation(domain, elevation=lambda t: t, indices = [0,1,3])

        operator()

        #pprint (domain.quantities['elevation'].centroid_values)
        #pprint (domain.quantities['stage'].centroid_values)

        elev_ex = [ 15.,  15.,   0.,  15.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.]



        stage_ex = [ 16.,  16.,   1.,  16.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,
         1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,
         1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,
         1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,
         1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,
         1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.]
        

#        from pprint import pprint
#        pprint (domain.quantities['elevation'].centroid_values)
#        pprint (domain.quantities['stage'].centroid_values)
#        pprint (domain.quantities['xmomentum'].centroid_values)
#        pprint (domain.quantities['ymomentum'].centroid_values)

        assert num.allclose(domain.quantities['elevation'].centroid_values, elev_ex)
        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)

            
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_set_elevation_operator, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
