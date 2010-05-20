#!/usr/bin/env python

import unittest, os
import os.path
from math import pi, sqrt

from swb_domain import *

from anuga.config import g

import numpy as num

# Get gateway to C implementation of flux function for direct testing
from swb_domain_ext import flux_function_c
from swb_domain_ext import gravity_c
from swb_domain_ext import compute_fluxes_c


class Test_swb_domain_ext(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_compute_fluxes0(self):
        # Do a full triangle and check that fluxes cancel out for
        # the constant stage case

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
        domain.set_quantity('stage', [[2,2,2], [2,2,2], [2,2,2], [2,2,2]])
        domain.check_integrity()
        domain.distribute_to_vertices_and_edges()
        
        assert num.allclose(domain.neighbours,
                            [[-1,1,-1], [2,3,0], [-1,-1,1],[1,-1,-1]])
        assert num.allclose(domain.neighbour_edges,
                            [[-1,2,-1], [2,0,1], [-1,-1,0],[1,-1,-1]])

        edgeflux = num.zeros(3, num.float)
        edgeflux0 = num.zeros(3, num.float)
        edgeflux1 = num.zeros(3, num.float)
        edgeflux2 = num.zeros(3, num.float)

        # Flux across right edge of volume 1
        normal = domain.get_normal(1, 0)
        ql = domain.get_evolved_quantities(vol_id=1, edge=0)
        qr = domain.get_evolved_quantities(vol_id=2, edge=2)

        #print ql
        #print qr
        local_max_speed = flux_function_c(normal, ql, qr, edgeflux0, g)

        # Check that flux seen from other triangles is inverse
        (ql, qr) = (qr, ql)
        normal = domain.get_normal(2, 2)
        local_max_speed = flux_function_c(normal, ql, qr, edgeflux, g)

        assert num.allclose(edgeflux0 + edgeflux, 0.)

        # Flux across upper edge of volume 1
        normal = domain.get_normal(1, 1)
        ql = domain.get_evolved_quantities(vol_id=1, edge=1)
        qr = domain.get_evolved_quantities(vol_id=3, edge=0)
        #print ql
        #print qr
        local_max_speed = flux_function_c(normal, ql, qr, edgeflux1, g)

        # Check that flux seen from other triangles is inverse
        (ql, qr) = (qr, ql)
        normal = domain.get_normal(3, 0)
        local_max_speed = flux_function_c(normal, ql, qr, edgeflux, g)

        assert num.allclose(edgeflux1 + edgeflux, 0.)

        # Flux across lower left hypotenuse of volume 1
        normal = domain.get_normal(1, 2)
        ql = domain.get_evolved_quantities(vol_id=1, edge=2)
        qr = domain.get_evolved_quantities(vol_id=0, edge=1)
        #print ql
        #print qr
        local_max_speed = flux_function_c(normal, ql, qr, edgeflux2, g)

        # Check that flux seen from other triangles is inverse
        (ql, qr) = (qr, ql)
        normal = domain.get_normal(0, 1)
        max_speed = flux_function_c(normal, ql, qr, edgeflux, g)
        assert num.allclose(edgeflux2 + edgeflux, 0.)

        # Scale by edgelengths, add up anc check that total flux is zero
        e0 = domain.edgelengths[1, 0]
        e1 = domain.edgelengths[1, 1]
        e2 = domain.edgelengths[1, 2]

        assert num.allclose(e0*edgeflux0 + e1*edgeflux1 + e2*edgeflux2, 0.)

        # Now check that compute_flux yields zeros as well
        w  = domain.quantities['stage']
        uh = domain.quantities['xmomentum']
        vh = domain.quantities['ymomentum']
        h  = domain.quantities['height']
        z  = domain.quantities['elevation']
        u  = domain.quantities['xvelocity']
        v  = domain.quantities['yvelocity']
        timestep = compute_fluxes_c(2.0,domain,w,uh,vh,h,z,u,v)

        assert num.allclose(timestep,  0.1064794275)


        for name in ['stage', 'xmomentum', 'ymomentum']:
            assert num.allclose(domain.quantities[name].explicit_update[1], 0.0)

    def test_compute_fluxes1(self):
        #Use values from previous version
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
        val0 = 2. + 2.0/3
        val1 = 4. + 4.0/3
        val2 = 8. + 2.0/3
        val3 = 2. + 8.0/3

        domain.set_quantity('stage', [[val0, val0, val0], [val1, val1, val1],
                                      [val2, val2, val2], [val3, val3, val3]])
        domain.check_integrity()
        domain.set_default_order(1)
        domain.distribute_to_vertices_and_edges()


        w  = domain.quantities['stage']
        uh = domain.quantities['xmomentum']
        vh = domain.quantities['ymomentum']
        h  = domain.quantities['height']
        z  = domain.quantities['elevation']
        u  = domain.quantities['xvelocity']
        v  = domain.quantities['yvelocity']


        edgeflux = num.zeros(3, num.float)
        edgeflux0 = num.zeros(3, num.float)
        edgeflux1 = num.zeros(3, num.float)
        edgeflux2 = num.zeros(3, num.float)


        # Flux across right edge of volume 1
        normal = domain.get_normal(1, 0)    # Get normal 0 of triangle 1
        assert num.allclose(normal, [1, 0])

  
        ql = domain.get_evolved_quantities(vol_id=1, edge=0)
        #print ql, val1
        assert num.allclose(ql, [val1, 0, 0, val1, 0, 0, 0])

        qr = domain.get_evolved_quantities(vol_id=2, edge=2)
        assert num.allclose(qr, [val2, 0, 0, val2, 0, 0, 0])

        local_max_speed = flux_function_c(normal, ql, qr, edgeflux0, g)

        #print edgeflux0
        #print local_max_speed
        
        # Flux across edge in the east direction (as per normal vector)
        assert num.allclose(edgeflux0, [-15.3598804, 253.71111111, 0.])
        assert num.allclose(local_max_speed, 9.21592824046)

        #Flux across edge in the west direction (opposite sign for xmomentum)
        normal_opposite = domain.get_normal(2, 2)   # Get normal 2 of triangle 2
        assert num.allclose(normal_opposite, [-1, 0])

        max_speed = flux_function_c(normal_opposite, ql, qr, edgeflux, g)
        assert num.allclose(edgeflux, [-15.3598804, -253.71111111, 0.])

        #Flux across upper edge of volume 1
        normal = domain.get_normal(1, 1)
        ql = domain.get_evolved_quantities(vol_id=1, edge=1)
        qr = domain.get_evolved_quantities(vol_id=3, edge=0)
        max_speed = flux_function_c(normal, ql, qr, edgeflux1, g)

        assert num.allclose(edgeflux1, [2.4098563, 0., 123.04444444])
        assert num.allclose(max_speed, 7.22956891292)

        #Flux across lower left hypotenuse of volume 1
        normal = domain.get_normal(1, 2)
        ql = domain.get_evolved_quantities(vol_id=1, edge=2)
        qr = domain.get_evolved_quantities(vol_id=0, edge=1)
        max_speed = flux_function_c(normal, ql, qr, edgeflux2, g)

        assert num.allclose(edgeflux2, [9.63942522, -61.59685738, -61.59685738])
        assert num.allclose(max_speed, 7.22956891292)

        #Scale, add up and check that compute_fluxes is correct for vol 1
        e0 = domain.edgelengths[1, 0]
        e1 = domain.edgelengths[1, 1]
        e2 = domain.edgelengths[1, 2]

        total_flux = -(e0*edgeflux0 +
                       e1*edgeflux1 +
                       e2*edgeflux2) / domain.areas[1]

        assert num.allclose(total_flux, [-0.68218178, -166.6, -35.93333333])


        timestep = compute_fluxes_c(2.0,domain,w,uh,vh,h,z,u,v)
        assert num.allclose(timestep, 0.0511510624314)

        for i, name in enumerate(['stage', 'xmomentum', 'ymomentum']):
            assert num.allclose(total_flux[i],
                                domain.quantities[name].explicit_update[1])


        assert num.allclose(domain.quantities['stage'].explicit_update,
                            [0., -0.68218178, -111.77316251, -35.68522449])
        assert num.allclose(domain.quantities['xmomentum'].explicit_update,
                            [-69.68888889, -166.6, 69.68888889, 0])
        assert num.allclose(domain.quantities['ymomentum'].explicit_update,
                            [-69.68888889, -35.93333333, 0., 69.68888889])

    def test_compute_fluxes2(self):
        #Random values, incl momentum
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
        domain.set_default_order(1)

        w  = domain.quantities['stage']
        uh = domain.quantities['xmomentum']
        vh = domain.quantities['ymomentum']
        h  = domain.quantities['height']
        z  = domain.quantities['elevation']
        u  = domain.quantities['xvelocity']
        v  = domain.quantities['yvelocity']

        zl = 10    # Assume flat bed
        val0 = zl + 2. + 2.0/3
        val1 = zl + 4. + 4.0/3
        val2 = zl + 8. + 2.0/3
        val3 = zl + 2. + 8.0/3


        edgeflux = num.zeros(3, num.float)
        edgeflux0 = num.zeros(3, num.float)
        edgeflux1 = num.zeros(3, num.float)
        edgeflux2 = num.zeros(3, num.float)

        domain.set_quantity('elevation', zl*num.ones((4, 3), num.int)) #array default#

        domain.set_quantity('stage', [[val0, val0-1, val0-2],
                                      [val1, val1+1, val1],
                                      [val2, val2-2, val2],
                                      [val3-0.5, val3, val3]])

        domain.set_quantity('xmomentum',
                            [[1,2,3], [3,4,5], [1,-1,0], [0,-2,2]])

        domain.set_quantity('ymomentum',
                            [[1,-1,0], [0,-3,2], [0,1,0], [-1,2,2]])

        domain.check_integrity()
        domain.distribute_to_vertices_and_edges()

        
        # Flux across right edge of volume 1
        normal = domain.get_normal(1, 0)
        ql = domain.get_evolved_quantities(vol_id=1, edge=0)
        qr = domain.get_evolved_quantities(vol_id=2, edge=2)
        max_speed = flux_function_c(normal, ql, qr, edgeflux0, g)

        # Flux across upper edge of volume 1
        normal = domain.get_normal(1, 1)
        ql = domain.get_evolved_quantities(vol_id=1, edge=1)
        qr = domain.get_evolved_quantities(vol_id=3, edge=0)
        max_speed = flux_function_c(normal, ql, qr, edgeflux1, g)

        # Flux across lower left hypotenuse of volume 1
        normal = domain.get_normal(1, 2)
        ql = domain.get_evolved_quantities(vol_id=1, edge=2)
        qr = domain.get_evolved_quantities(vol_id=0, edge=1)
        max_speed = flux_function_c(normal, ql, qr, edgeflux2, g)

        # Scale, add up and check that compute_fluxes is correct for vol 1
        e0 = domain.edgelengths[1, 0]
        e1 = domain.edgelengths[1, 1]
        e2 = domain.edgelengths[1, 2]

        total_flux = -(e0*edgeflux0 +
                       e1*edgeflux1 +
                       e2*edgeflux2) / domain.areas[1]

        timestep = compute_fluxes_c(2.0,domain,w,uh,vh,h,z,u,v)
        assert num.allclose(timestep, 0.0529903533177)

        for i, name in enumerate(['stage', 'xmomentum', 'ymomentum']):
            assert num.allclose(total_flux[i],
                                domain.quantities[name].explicit_update[1])

    # FIXME (Ole): Need test like this for fluxes in very shallow water.
    def test_compute_fluxes3(self):
        #Random values, incl momentum
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
        domain.set_default_order(1)


        w  = domain.quantities['stage']
        uh = domain.quantities['xmomentum']
        vh = domain.quantities['ymomentum']
        h  = domain.quantities['height']
        z  = domain.quantities['elevation']
        u  = domain.quantities['xvelocity']
        v  = domain.quantities['yvelocity']
        
        val0 = 2.+2.0/3
        val1 = 4.+4.0/3
        val2 = 8.+2.0/3
        val3 = 2.+8.0/3

        zl = -3.75    # Assume constant bed (must be less than stage)
        domain.set_quantity('elevation', zl*num.ones((4, 3), num.int)) #array default#

        edgeflux = num.zeros(3, num.float)
        edgeflux0 = num.zeros(3, num.float)
        edgeflux1 = num.zeros(3, num.float)
        edgeflux2 = num.zeros(3, num.float)


        domain.set_quantity('stage', [[val0, val0-1, val0-2],
                                      [val1, val1+1, val1],
                                      [val2, val2-2, val2],
                                      [val3-0.5, val3, val3]])

        domain.set_quantity('xmomentum',
                            [[1,2,3], [3,4,5], [1,-1,0], [0,-2,2]])

        domain.set_quantity('ymomentum',
                            [[1,-1,0], [0,-3,2], [0,1,0], [-1,2,2]])

        domain.check_integrity()
        domain.distribute_to_vertices_and_edges()        

        # Flux across right edge of volume 1
        normal = domain.get_normal(1, 0)
        ql = domain.get_evolved_quantities(vol_id=1, edge=0)
        qr = domain.get_evolved_quantities(vol_id=2, edge=2)
        max_speed = flux_function_c(normal, ql, qr, edgeflux0, g)

        #print max_speed, edgeflux0
        assert num.allclose(edgeflux0, [-10.51926294, 577.81462335, -3.64772873], rtol=1.0e-2)

        # Flux across upper edge of volume 1
        normal = domain.get_normal(1, 1)
        ql = domain.get_evolved_quantities(vol_id=1, edge=1)
        qr = domain.get_evolved_quantities(vol_id=3, edge=0)
        max_speed = flux_function_c(normal, ql, qr, edgeflux1, g)

        #print max_speed, edgeflux1
        assert num.allclose(edgeflux1, [ 5.93945967, 19.14204647,  377.48851296], rtol=1.0e-2)
        
        # Flux across lower left hypotenuse of volume 1
        normal = domain.get_normal(1, 2)
        ql = domain.get_evolved_quantities(vol_id=1, edge=2)
        qr = domain.get_evolved_quantities(vol_id=0, edge=1)
        max_speed = flux_function_c(normal, ql, qr, edgeflux2, g)

        #print max_speed, edgeflux2
        assert num.allclose(edgeflux2, [17.21047971, -192.90578267, -203.05771337], rtol=1.0e-2)

        
        # Scale, add up and check that compute_fluxes is correct for vol 1
        e0 = domain.edgelengths[1, 0]
        e1 = domain.edgelengths[1, 1]
        e2 = domain.edgelengths[1, 2]

        total_flux = -(e0*edgeflux0 +
                       e1*edgeflux1 +
                       e2*edgeflux2) / domain.areas[1]

        timestep = compute_fluxes_c(100.0,domain,w,uh,vh,h,z,u,v)

        assert num.allclose(timestep, 0.0438142267244, rtol=1.0e-2)

        for i, name in enumerate(['stage', 'xmomentum', 'ymomentum']):
            assert num.allclose(total_flux[i],
                                domain.quantities[name].explicit_update[1])



    def test_compute_fluxes4(self):
        #Random values, incl momentum
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
        domain.set_default_order(1)


        w  = domain.quantities['stage']
        uh = domain.quantities['xmomentum']
        vh = domain.quantities['ymomentum']
        h  = domain.quantities['height']
        z  = domain.quantities['elevation']
        u  = domain.quantities['xvelocity']
        v  = domain.quantities['yvelocity']
        
        val0 = 2.+2.0/3
        val1 = 4.+4.0/3
        val2 = 8.+2.0/3
        val3 = 2.+8.0/3

        zl = -3.75    # Assume constant bed (must be less than stage)
        domain.set_quantity('elevation', [0.0, -2.0, 0.0 ,0.0], location='centroids') 

        edgeflux = num.zeros(3, num.float)
        edgeflux0 = num.zeros(3, num.float)
        edgeflux1 = num.zeros(3, num.float)
        edgeflux2 = num.zeros(3, num.float)


        domain.set_quantity('stage', -1.0)

        domain.set_quantity('xmomentum',0.0)

        domain.set_quantity('ymomentum',0.0)

        domain.check_integrity()
        domain.distribute_to_vertices_and_edges()        


        print w.edge_values
        print u.edge_values
        print v.edge_values
        print h.edge_values
        print z.edge_values


        domain.compute_fluxes()

        for i, name in enumerate(['stage', 'xmomentum', 'ymomentum']):
            print domain.quantities[name].explicit_update[1]
            #assert num.allclose(total_flux[i],
            #                    domain.quantities[name].explicit_update[1])


#################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_swb_domain_ext, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
