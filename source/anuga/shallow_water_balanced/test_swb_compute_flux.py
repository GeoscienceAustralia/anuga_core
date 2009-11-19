#!/usr/bin/env python

import unittest, os
import os.path
from math import pi, sqrt
import tempfile

from anuga.config import g, epsilon
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a
from anuga.utilities.numerical_tools import mean
from anuga.utilities.polygon import is_inside_polygon
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.abstract_2d_finite_volumes.quantity import Quantity
from anuga.geospatial_data.geospatial_data import Geospatial_data
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross

from anuga.utilities.system_tools import get_pathname_from_package
from swb_domain import *

import numpy as num

# Get gateway to C implementation of flux function for direct testing
from shallow_water_ext import flux_function_central as flux_function




class Test_swb_clean(unittest.TestCase):
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

        assert num.allclose(domain.neighbours,
                            [[-1,1,-1], [2,3,0], [-1,-1,1],[1,-1,-1]])
        assert num.allclose(domain.neighbour_edges,
                            [[-1,2,-1], [2,0,1], [-1,-1,0],[1,-1,-1]])

        zl = zr = 0.     # Assume flat bed

        edgeflux = num.zeros(3, num.float)
        edgeflux0 = num.zeros(3, num.float)
        edgeflux1 = num.zeros(3, num.float)
        edgeflux2 = num.zeros(3, num.float)
        H0 = 0.0

        # Flux across right edge of volume 1
        normal = domain.get_normal(1, 0)
        ql = domain.get_conserved_quantities(vol_id=1, edge=0)
        qr = domain.get_conserved_quantities(vol_id=2, edge=2)
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux0,
                                  epsilon, g, H0)

        # Check that flux seen from other triangles is inverse
        (ql, qr) = (qr, ql)
        normal = domain.get_normal(2, 2)
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux,
                                  epsilon, g, H0)

        assert num.allclose(edgeflux0 + edgeflux, 0.)

        # Flux across upper edge of volume 1
        normal = domain.get_normal(1, 1)
        ql = domain.get_conserved_quantities(vol_id=1, edge=1)
        qr = domain.get_conserved_quantities(vol_id=3, edge=0)
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux1,
                                  epsilon, g, H0)

        # Check that flux seen from other triangles is inverse
        (ql, qr) = (qr, ql)
        normal = domain.get_normal(3, 0)
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux,
                                  epsilon, g, H0)

        assert num.allclose(edgeflux1 + edgeflux, 0.)

        # Flux across lower left hypotenuse of volume 1
        normal = domain.get_normal(1, 2)
        ql = domain.get_conserved_quantities(vol_id=1, edge=2)
        qr = domain.get_conserved_quantities(vol_id=0, edge=1)
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux2,
                                  epsilon, g, H0)

        # Check that flux seen from other triangles is inverse
        (ql, qr) = (qr, ql)
        normal = domain.get_normal(0, 1)
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux,
                                  epsilon, g, H0)
        assert num.allclose(edgeflux2 + edgeflux, 0.)

        # Scale by edgelengths, add up anc check that total flux is zero
        e0 = domain.edgelengths[1, 0]
        e1 = domain.edgelengths[1, 1]
        e2 = domain.edgelengths[1, 2]

        assert num.allclose(e0*edgeflux0 + e1*edgeflux1 + e2*edgeflux2, 0.)

        # Now check that compute_flux yields zeros as well
        domain.compute_fluxes()

        for name in ['stage', 'xmomentum', 'ymomentum']:
            assert num.allclose(domain.quantities[name].explicit_update[1], 0)

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

        zl = zr = 0.    # Assume flat bed

        edgeflux = num.zeros(3, num.float)
        edgeflux0 = num.zeros(3, num.float)
        edgeflux1 = num.zeros(3, num.float)
        edgeflux2 = num.zeros(3, num.float)
        H0 = 0.0

        # Flux across right edge of volume 1
        normal = domain.get_normal(1, 0)    # Get normal 0 of triangle 1
        assert num.allclose(normal, [1, 0])

        ql = domain.get_conserved_quantities(vol_id=1, edge=0)
        assert num.allclose(ql, [val1, 0, 0])

        qr = domain.get_conserved_quantities(vol_id=2, edge=2)
        assert num.allclose(qr, [val2, 0, 0])

        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux0,
                                  epsilon, g, H0)

        # Flux across edge in the east direction (as per normal vector)
        assert num.allclose(edgeflux0, [-15.3598804, 253.71111111, 0.])
        assert num.allclose(max_speed, 9.21592824046)

        #Flux across edge in the west direction (opposite sign for xmomentum)
        normal_opposite = domain.get_normal(2, 2)   # Get normal 2 of triangle 2
        assert num.allclose(normal_opposite, [-1, 0])

        max_speed = flux_function(normal_opposite, ql, qr, zl, zr, edgeflux,
                                  epsilon, g, H0)
        assert num.allclose(edgeflux, [-15.3598804, -253.71111111, 0.])

        #Flux across upper edge of volume 1
        normal = domain.get_normal(1, 1)
        ql = domain.get_conserved_quantities(vol_id=1, edge=1)
        qr = domain.get_conserved_quantities(vol_id=3, edge=0)
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux1,
                                  epsilon, g, H0)

        assert num.allclose(edgeflux1, [2.4098563, 0., 123.04444444])
        assert num.allclose(max_speed, 7.22956891292)

        #Flux across lower left hypotenuse of volume 1
        normal = domain.get_normal(1, 2)
        ql = domain.get_conserved_quantities(vol_id=1, edge=2)
        qr = domain.get_conserved_quantities(vol_id=0, edge=1)
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux2,
                                  epsilon, g, H0)

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

        domain.compute_fluxes()

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
        val0 = 2. + 2.0/3
        val1 = 4. + 4.0/3
        val2 = 8. + 2.0/3
        val3 = 2. + 8.0/3

        zl = zr = 0    # Assume flat zero bed
        edgeflux = num.zeros(3, num.float)
        edgeflux0 = num.zeros(3, num.float)
        edgeflux1 = num.zeros(3, num.float)
        edgeflux2 = num.zeros(3, num.float)
        H0 = 0.0

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

        # Flux across right edge of volume 1
        normal = domain.get_normal(1, 0)
        ql = domain.get_conserved_quantities(vol_id=1, edge=0)
        qr = domain.get_conserved_quantities(vol_id=2, edge=2)
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux0,
                                  epsilon, g, H0)

        # Flux across upper edge of volume 1
        normal = domain.get_normal(1, 1)
        ql = domain.get_conserved_quantities(vol_id=1, edge=1)
        qr = domain.get_conserved_quantities(vol_id=3, edge=0)
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux1,
                                  epsilon, g, H0)

        # Flux across lower left hypotenuse of volume 1
        normal = domain.get_normal(1, 2)
        ql = domain.get_conserved_quantities(vol_id=1, edge=2)
        qr = domain.get_conserved_quantities(vol_id=0, edge=1)
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux2,
                                  epsilon, g, H0)

        # Scale, add up and check that compute_fluxes is correct for vol 1
        e0 = domain.edgelengths[1, 0]
        e1 = domain.edgelengths[1, 1]
        e2 = domain.edgelengths[1, 2]

        total_flux = -(e0*edgeflux0 +
                       e1*edgeflux1 +
                       e2*edgeflux2) / domain.areas[1]

        domain.compute_fluxes()

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

        val0 = 2.+2.0/3
        val1 = 4.+4.0/3
        val2 = 8.+2.0/3
        val3 = 2.+8.0/3

        zl = zr = -3.75    # Assume constant bed (must be less than stage)
        domain.set_quantity('elevation', zl*num.ones((4, 3), num.int)) #array default#

        edgeflux = num.zeros(3, num.float)
        edgeflux0 = num.zeros(3, num.float)
        edgeflux1 = num.zeros(3, num.float)
        edgeflux2 = num.zeros(3, num.float)
        H0 = 0.0

        domain.set_quantity('stage', [[val0, val0-1, val0-2],
                                      [val1, val1+1, val1],
                                      [val2, val2-2, val2],
                                      [val3-0.5, val3, val3]])

        domain.set_quantity('xmomentum',
                            [[1,2,3], [3,4,5], [1,-1,0], [0,-2,2]])

        domain.set_quantity('ymomentum',
                            [[1,-1,0], [0,-3,2], [0,1,0], [-1,2,2]])

        domain.check_integrity()

        # Flux across right edge of volume 1
        normal = domain.get_normal(1, 0)
        ql = domain.get_conserved_quantities(vol_id=1, edge=0)
        qr = domain.get_conserved_quantities(vol_id=2, edge=2)
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux0,
                                  epsilon, g, H0)

        # Flux across upper edge of volume 1
        normal = domain.get_normal(1, 1)
        ql = domain.get_conserved_quantities(vol_id=1, edge=1)
        qr = domain.get_conserved_quantities(vol_id=3, edge=0)
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux1,
                                  epsilon, g, H0)

        # Flux across lower left hypotenuse of volume 1
        normal = domain.get_normal(1, 2)
        ql = domain.get_conserved_quantities(vol_id=1, edge=2)
        qr = domain.get_conserved_quantities(vol_id=0, edge=1)
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux2,
                                  epsilon, g, H0)

        # Scale, add up and check that compute_fluxes is correct for vol 1
        e0 = domain.edgelengths[1, 0]
        e1 = domain.edgelengths[1, 1]
        e2 = domain.edgelengths[1, 2]

        total_flux = -(e0*edgeflux0 +
                       e1*edgeflux1 +
                       e2*edgeflux2) / domain.areas[1]

        domain.compute_fluxes()

        for i, name in enumerate(['stage', 'xmomentum', 'ymomentum']):
            assert num.allclose(total_flux[i],
                                domain.quantities[name].explicit_update[1])

    def test_flux_optimisation(self):
        """test_flux_optimisation

        Test that fluxes are correctly computed using
        dry and still cell exclusions
        """

        from anuga.config import g
        import copy

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

        #Set up for a gradient of (3,0) at mid triangle (bce)
        def slope(x, y):
            return 3*x

        h = 0.1
        def stage(x, y):
            return slope(x, y) + h

        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage)

        # Allow slope limiters to work
        # (FIXME (Ole): Shouldn't this be automatic in ANUGA?)
        domain.distribute_to_vertices_and_edges()

        initial_stage = copy.copy(domain.quantities['stage'].vertex_values)

        domain.set_boundary({'exterior': Reflective_boundary(domain)})

        #  Check that update arrays are initialised to zero
        assert num.allclose(domain.get_quantity('stage').explicit_update, 0)
        assert num.allclose(domain.get_quantity('xmomentum').explicit_update, 0)
        assert num.allclose(domain.get_quantity('ymomentum').explicit_update, 0)

        # Get true values
        domain.optimise_dry_cells = False
        domain.compute_fluxes()
        stage_ref = copy.copy(domain.get_quantity('stage').explicit_update)
        xmom_ref = copy.copy(domain.get_quantity('xmomentum').explicit_update)
        ymom_ref = copy.copy(domain.get_quantity('ymomentum').explicit_update)

        # Try with flux optimisation
        domain.optimise_dry_cells = True
        domain.compute_fluxes()

        assert num.allclose(stage_ref,
                            domain.get_quantity('stage').explicit_update)
        assert num.allclose(xmom_ref,
                            domain.get_quantity('xmomentum').explicit_update)
        assert num.allclose(ymom_ref,
                            domain.get_quantity('ymomentum').explicit_update)


#################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_swb_clean, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
