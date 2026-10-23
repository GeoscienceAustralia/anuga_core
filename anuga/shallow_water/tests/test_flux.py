#!/usr/bin/env python
"""Tests for flux functions and compute_fluxes — direct C interface tests."""
import unittest
from math import pi, sqrt
import numpy as num

from anuga.config import g, epsilon
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross, rectangular
from anuga.abstract_2d_finite_volumes.quantity import Quantity
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.shallow_water.boundaries import Reflective_boundary
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions import Dirichlet_boundary

from anuga.shallow_water.sw_domain_openmp_ext import flux_function_central as flux_function
from anuga.shallow_water.sw_domain_openmp_ext import rotate

class Test_Shallow_Water(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


    def test_rotate(self):
        normal = num.array([0.0, -1.0])

        q = num.array([1.0, 2.0, 3.0])

        r = rotate(q, normal, direction = 1)
        assert r[0] == 1
        assert r[1] == -3
        assert r[2] == 2

        w = rotate(r, normal, direction = -1)
        assert num.allclose(w, q)

        # Check error check
        try:
            rotate(r, num.array([1, 1, 1]))
        except Exception:
            pass
        else:
            raise Exception('Should have raised an exception')

    # Individual flux tests


    def test_flux_zero_case(self):
        ql = num.zeros(3, float)
        qr = num.zeros(3, float)
        normal = num.zeros(2, float)
        edgeflux = num.zeros(3, float)
        zl = zr = ze = 0.
        H0 = 1.0e-3 # As suggested in the manual

        h = hl = hr = hle = hre = hc = hc_n = 0.0
        low_froude = 1

        max_speed, pressure_flux = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux, epsilon, ze, g, H0, hc, hc_n, low_froude)

        assert num.allclose(edgeflux, [0, 0, 0])
        assert num.allclose(pressure_flux, 0.5*g*h**2)
        assert max_speed == 0.


    def test_flux_constants(self):
        w = 2.0

        normal = num.array([1.,0])
        ql = num.array([w, 0, 0])
        qr = num.array([w, 0, 0])
        edgeflux = num.zeros(3, float)
        zr = zl = ze = 0
        hl = hr = h = w - (zl+zr) / 2
        H0 = 0.0

        hle = hre = hc = h
        hc_n = 0

        for low_froude in [1, 2]:
            max_speed, pressure_flux = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux, epsilon, ze, g, H0, hc, hc_n, low_froude)
            assert max_speed == num.sqrt(g*h)
            assert num.allclose(pressure_flux, 0.5*g*h**2)
            assert num.allclose(edgeflux, [0., 0., 0.])


    def test_flux1(self):
        # Use data from previous version of abstract_2d_finite_volumes
        normal = num.array([1., 0])
        ql = num.array([-0.2, 2, 3])
        qr = num.array([-0.2, 2, 3])

        zl = zr = ze = -0.5
        hc = hc_n = hl = hr = hle = hre = h = 0.3
        low_froude = 2
        H0 = 0.0
        edgeflux = num.zeros(3, float)
        max_speed, pressure_flux = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux, epsilon, ze, g, H0, hc, hc_n, low_froude)

        assert num.allclose(max_speed, 8.38130948661)
        assert num.allclose(edgeflux, [2., 13.3333333, 20.])


    def test_flux2(self):
        # Use data from previous version of abstract_2d_finite_volumes
        normal = num.array([0., -1.])
        ql = num.array([-0.075, 2, 3])
        qr = num.array([-0.075, 2, 3])
        zl = zr = ze = -0.375

        edgeflux = num.zeros(3, float)
        H0 = 0.0

        hl = hle = ql[0] - zl
        hr = hre = qr[0] - zr
        h = hc = hc_n = (hl + hr)/2
        low_froude = 2
        H0 = 0.0

        max_speed, pressure_flux = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux, epsilon, ze, g, H0, hc, hc_n, low_froude)

        assert num.allclose(edgeflux, [-3., -20.0, -30.0])
        assert num.allclose(max_speed, 11.7146428199)
        assert num.allclose(pressure_flux, 0.441)


    def test_flux3(self):
        # Use data from previous version of abstract_2d_finite_volumes
        normal = num.array([-sqrt(2) / 2, sqrt(2) / 2])
        ql = num.array([-0.075, 2, 3])
        qr = num.array([-0.075, 2, 3])
        zl = zr = ze = -0.375

        edgeflux = num.zeros(3, float)

        hl = hle = ql[0] - zl
        hr = hre = qr[0] - zr
        h = hc = hc_n = (hl + hr)/2
        low_froude = 2
        H0 = 0.0

        max_speed, pressure_flux = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux, epsilon, ze, g, H0, hc, hc_n, low_froude)

        assert num.allclose(edgeflux, [sqrt(2) / 2, 4.71404521, 7.07106781])
        assert num.allclose(max_speed, 4.0716654239)
        assert num.allclose(pressure_flux, 0.441)


    def test_flux4(self):
        # Use data from previous version of abstract_2d_finite_volumes
        normal = num.array([-sqrt(2) / 2, sqrt(2) / 2])
        ql = num.array([-0.34319278, 0.10254161, 0.07273855])
        qr = num.array([-0.30683287, 0.1071986, 0.05930515])
        zl = zr = ze = -0.375

        edgeflux = num.zeros(3, float)
        H0 = 0.0

        hl = hle = ql[0] - zl
        hr = hre = qr[0] - zr
        h = hc = hc_n = (hl + hr)/2
        low_froude = 2
        H0 = 0.0

        max_speed, pressure_flux = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux, epsilon, ze, g, H0, hc, hc_n, low_froude)
        assert num.allclose(max_speed, 1.31414103233)
        assert num.allclose(edgeflux, [-0.04072676, -0.05733579, -0.02967421])
        assert num.allclose(pressure_flux, 0.0192765311)


    def test_flux_computation(self):
        """test flux calculation (actual C implementation)

        This one tests the constant case where only the pressure term
        contributes to each edge and cancels out once the total flux has
        been summed up.
        """

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #              bac,     bce,     ecf,     dbe
        vertices = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)
        domain.check_integrity()

        # The constant case
        domain.set_quantity('elevation', -1)
        domain.set_quantity('stage', 1)

        domain.compute_fluxes()
        # Central triangle
        assert num.allclose(domain.get_quantity('stage').explicit_update[1], 0)

        # The more general case
        def surface(x, y):
            return -x / 2

        domain.set_quantity('elevation', -10)
        domain.set_quantity('stage', surface)
        domain.set_quantity('xmomentum', 1)

        domain.compute_fluxes()

        # FIXME (Ole): TODO the general case
        #assert allclose(domain.get_quantity('stage').explicit_update[1], ...??)


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


    def test_compute_fluxes_structure_0(self):
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
        domain.set_quantity('stage', [[2,4,6], [4,3,2], [2,0,1], [1,2,3]])
        domain.check_integrity()

        assert num.allclose(domain.neighbours,
                            [[-1,1,-2], [2,3,0], [-3,-4,1],[1,-5,-6]])
        assert num.allclose(domain.neighbour_edges,
                            [[-1,2,-1], [2,0,1], [-1,-1,0],[1,-1,-1]])

        ze = zl = zr = 0.     # Assume flat bed

        edgeflux = num.zeros(3, float)
        edgeflux0 = num.zeros(3, float)
        edgeflux1 = num.zeros(3, float)
        edgeflux2 = num.zeros(3, float)
        H0 = 0.0

        # Flux across right edge of volume 1
        normal = domain.get_normal(1, 0)
        ql = domain.get_conserved_quantities(vol_id=1, edge=0)
        qr = domain.get_conserved_quantities(vol_id=2, edge=2)

        hl = hle = ql[0] - zl
        hr = hre = qr[0] - zr
        hc = hc_n = 2
        low_froude = 1

        max_speed, pressure_flux = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux0, epsilon, ze, g, H0, hc, hc_n, low_froude)
        assert(num.any(edgeflux0 != 0))

        # Check that flux seen from other triangles is inverse
        (ql, qr) = (qr, ql)
        normal = domain.get_normal(2, 2)

        max_speed, pressure_flux = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux, epsilon, ze, g, H0, hc, hc_n, low_froude)
        assert(num.any(edgeflux != 0))
        assert num.allclose(edgeflux0 + edgeflux, 0.)

        # Flux across upper edge of volume 1
        normal = domain.get_normal(1, 1)
        ql = domain.get_conserved_quantities(vol_id=1, edge=1)
        qr = domain.get_conserved_quantities(vol_id=3, edge=0)
        max_speed, pressure_flux = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux1, epsilon, ze, g, H0, hc, hc_n, low_froude)
        assert(num.any(edgeflux1 != 0))

        # Check that flux seen from other triangles is inverse
        (ql, qr) = (qr, ql)
        normal = domain.get_normal(3, 0)

        max_speed, pressure_flux = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux, epsilon, ze, g, H0, hc, hc_n, low_froude)
        assert(num.any(edgeflux != 0))
        assert num.allclose(edgeflux1 + edgeflux, 0.)

        # Flux across lower left hypotenuse of volume 1
        normal = domain.get_normal(1, 2)
        ql = domain.get_conserved_quantities(vol_id=1, edge=2)
        qr = domain.get_conserved_quantities(vol_id=0, edge=1)

        max_speed, pressure_flux = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux2, epsilon, ze, g, H0, hc, hc_n, low_froude)
        assert(num.any(edgeflux2 != 0))

        # Check that flux seen from other triangles is inverse
        (ql, qr) = (qr, ql)
        normal = domain.get_normal(0, 1)
        max_speed, pressure_flux = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux, epsilon, ze, g, H0, hc, hc_n, low_froude)
        assert(num.any(edgeflux != 0))
        assert num.allclose(edgeflux2 + edgeflux, 0.)

        # Now check that compute_fluxes is correct as well
        domain.compute_fluxes()
        for name in ['stage', 'xmomentum', 'ymomentum']:
            assert num.allclose(domain.quantities[name].explicit_update[1], 0)


    def test_compute_fluxes_structure_1(self):
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

        domain.set_quantity('height', [[val0, val0, val0], [val1, val1, val1],
                                      [val2, val2, val2], [val3, val3, val3]])

        domain.check_integrity()

        zl = zr = zc = 0.    # Assume flat bed

        edgeflux = num.zeros(3, float)
        edgeflux0 = num.zeros(3, float)
        edgeflux1 = num.zeros(3, float)
        edgeflux2 = num.zeros(3, float)
        H0 = domain.H0

        # Flux across right edge of volume 1
        normal = domain.get_normal(1, 0)    # Get normal 0 of triangle 1
        assert num.allclose(normal, [1, 0])

        ql = domain.get_conserved_quantities(vol_id=1, edge=0)
        assert num.allclose(ql, [val1, 0, 0])

        qr = domain.get_conserved_quantities(vol_id=2, edge=2)
        assert num.allclose(qr, [val2, 0, 0])

        ze = (zl + zr) / 2
        hl = hle = val1 - zl
        hr = hre = val2 - zr
        hc = hc_n = (hl + hr) / 2
        low_froude = 0
        max_speed, pressure_flux0 = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux0, epsilon, ze, g, H0, hc, hc_n, low_froude)

        # Flux across edge in the east direction (as per normal vector)
        assert num.allclose(max_speed, 9.21592824046)
        assert num.allclose(pressure_flux0, 253.71111111)
        assert num.allclose(edgeflux0, [-15.3598804, 0, 0.])

        # Flux across edge in the west direction (opposite sign for xmomentum)
        normal_opposite = domain.get_normal(2, 2)   # Get normal 2 of triangle 2
        assert num.allclose(normal_opposite, [-1, 0])

        max_speed, pressure_flux0 = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux0, epsilon, ze, g, H0, hc, hc_n, low_froude)
        assert num.allclose(max_speed, 9.21592824046)
        assert num.allclose(pressure_flux0, 253.71111111)

        # FIXME (Ole): edgeflux used to have direction - now the two last compenents are zero and pressure_flux has a corresponding scalar value).
        #
        # The test used to be this:
        # assert num.allclose(edgeflux, [-15.3598804, -253.71111111, 0.])

        # Now it is this
        #assert num.allclose(edgeflux0, [-15.3598804, 0, 0.])

        bedslope_work0 = -g*0.5*(hl*hl - hle*hle -(hle+hc)*(zl-zc))+pressure_flux0
        edgeflux0[1] += normal[0]*bedslope_work0
        edgeflux0[2] += normal[1]*bedslope_work0

        assert num.allclose(edgeflux0, [-15.3598804, 253.71111111, 0.])

        # Flux across upper edge of volume 1
        normal = domain.get_normal(1, 1)
        ql = domain.get_conserved_quantities(vol_id=1, edge=1)
        qr = domain.get_conserved_quantities(vol_id=3, edge=0)
        hl = hle = ql[0] - zl
        hr = hre = qr[0] - zr
        hc = hc_n = (hl + hr) / 2
        max_speed, pressure_flux1 = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux1, epsilon, ze, g, H0, hc, hc_n, low_froude)


        # FIXME (Ole): edgeflux used to have direction - now it the two last components are zero and pressure_flux has a corresponding scalar value)                           assert num.allclose(edgeflux1, [2.4098563, 0., 123.04444444])

        # Now it is this instead:
        assert num.allclose(max_speed, 7.22956891292)
        assert num.allclose(pressure_flux1, 123.04444444)
        assert num.allclose(edgeflux1, [2.4098563, 0., 0.])

        bedslope_work1 = -g*0.5*(hl*hl - hle*hle -(hle+hc)*(zl-zc))+pressure_flux1
        edgeflux1[1] += normal[0]*bedslope_work1
        edgeflux1[2] += normal[1]*bedslope_work1

        assert num.allclose(edgeflux1, [2.4098563, 0., 123.04444444])

        # Flux across lower left hypotenuse of volume 1
        normal = domain.get_normal(1, 2)
        ql = domain.get_conserved_quantities(vol_id=1, edge=2)
        qr = domain.get_conserved_quantities(vol_id=0, edge=1)

        hl = hle = ql[0] - zl
        hr = hre = qr[0] - zr
        hc = hc_n = (hl + hr) / 2
        max_speed, pressure_flux2 = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux2, epsilon, ze, g, H0, hc, hc_n, low_froude)

        # FIXME (Ole): The test changed from this
        #assert num.allclose(edgeflux2, [9.63942522, -61.59685738, -61.59685738])
        #assert num.allclose(max_speed, 7.22956891292)

        # To this:
        assert num.allclose(max_speed, 7.22956891292)
        assert num.allclose(pressure_flux2, 87.111111111)
        assert num.allclose(edgeflux2, [9.63942522, 0., 0.])


        bedslope_work2 = -g*0.5*(hl*hl - hle*hle -(hle+hc)*(zl-zc))+pressure_flux2
        edgeflux2[1] += normal[0]*bedslope_work2
        edgeflux2[2] += normal[1]*bedslope_work2

        assert num.allclose(edgeflux2, [9.63942522, -61.59685738, -61.59685738])

        # FIXME (Ole): This type of test was great for internal integrity and should really be reinstated by incorporating pressure_flux

        # Scale, add up and check that compute_fluxes is correct for vol 1
        e0 = domain.edgelengths[1, 0]
        e1 = domain.edgelengths[1, 1]
        e2 = domain.edgelengths[1, 2]

        total_flux = -(e0*edgeflux0 +
                       e1*edgeflux1 +
                       e2*edgeflux2) / domain.areas[1]
        assert num.allclose(total_flux, [-0.68218178, -166.6, -35.93333333])

        # Now check that compute_flux yields the same
        domain.compute_fluxes()

        for i, name in enumerate(['stage', 'xmomentum', 'ymomentum']):

            #                    domain.quantities[name].explicit_update[1])
            assert num.allclose(total_flux[i],
                                domain.quantities[name].explicit_update[1])

        assert num.allclose(domain.quantities['stage'].explicit_update,
                            [0., -0.68218178, -111.77316251, -35.68522449])
        assert num.allclose(domain.quantities['xmomentum'].explicit_update,
                            [-69.68888889, -166.6, 69.68888889, 0])
        assert num.allclose(domain.quantities['ymomentum'].explicit_update,
                            [-69.68888889, -35.93333333, 0., 69.68888889])


    def test_compute_fluxes_structure_2(self):
        # Random values, incl momentum
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

        ze = zl = zr = zc = 0    # Assume flat zero bed

        #ze = (zl + zr) / 2

        edgeflux = num.zeros(3, float)
        edgeflux0 = num.zeros(3, float)
        edgeflux1 = num.zeros(3, float)
        edgeflux2 = num.zeros(3, float)
        ql = num.zeros(3, float)
        qr = num.zeros(3, float)

        H0 = domain.H0
        low_froude = 0

        domain.set_quantity('elevation', zc*num.ones((4, 3), int)) #array default#

        domain.set_quantity('stage', [[val0, val0-1, val0-2],
                                      [val1, val1+1, val1],
                                      [val2, val2-2, val2],
                                      [val3-0.5, val3, val3]])

        domain.set_quantity('height', [[val0, val0-1, val0-2],
                                      [val1, val1+1, val1],
                                      [val2, val2-2, val2],
                                      [val3-0.5, val3, val3]])

        domain.set_quantity('xmomentum',
                            [[1,2,3], [3,4,5], [1,-1,0], [0,-2,2]])

        domain.set_quantity('ymomentum',
                            [[1,-1,0], [0,-3,2], [0,1,0], [-1,2,2]])

        domain.check_integrity()

        Stage  = domain.quantities['stage']
        Height = domain.quantities['height']
        Bed    = domain.quantities['elevation']
        Xmom   = domain.quantities['xmomentum']
        Ymom   = domain.quantities['ymomentum']

        hc = Height.centroid_values[1]
        zc = Bed.centroid_values[1]


        # Flux across right edge of volume 1
        normal = domain.get_normal(1, 0)
        # ql = domain.get_conserved_quantities(vol_id=1, edge=0)
        # qr = domain.get_conserved_quantities(vol_id=2, edge=2)

        # hl = hle = ql[0] - zl
        # hr = hre = qr[0] - zr
        # hc = hc_n = (hl + hr) / 2


        ql[0] = Stage.edge_values[1,0]
        ql[1] = Xmom.edge_values[1,0]
        ql[2] = Ymom.edge_values[1,0]
        zl  = Bed.edge_values[1,0]
        hle = Height.edge_values[1,0]

        hc_n = Height.centroid_values[2]
        zc_n = Bed.centroid_values[2]

        qr[0] = Stage.edge_values[2,2]
        qr[1] = Xmom.edge_values[2,2]
        qr[2] = Ymom.edge_values[2,2]
        zr  = Bed.edge_values[2,2]
        hre = Height.edge_values[2,2]

        z_half = max(zl, zr)

        hl = max(hle+zl-z_half,0.0)
        hr = max(hre+zr-z_half,0.0)

        low_froude = 0

        max_speed0, pressure_flux0 = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux0, epsilon, z_half, g, H0, hc, hc_n, low_froude)

        bedslope_work0 = -g*0.5*(hl*hl - hle*hle -(hle+hc)*(zl-zc))+pressure_flux0
        edgeflux0[1] += normal[0]*bedslope_work0
        edgeflux0[2] += normal[1]*bedslope_work0


        # Flux across upper edge of volume 1
        normal = domain.get_normal(1, 1)
        # ql = domain.get_conserved_quantities(vol_id=1, edge=1)
        # qr = domain.get_conserved_quantities(vol_id=3, edge=0)

        # hl = hle = ql[0] - zl
        # hr = hre = qr[0] - zr
        # hc = hc_n = (hl + hr) / 2


        ql[0] = Stage.edge_values[1,1]
        ql[1] = Xmom.edge_values[1,1]
        ql[2] = Ymom.edge_values[1,1]
        zl  = Bed.edge_values[1,1]
        hl = hle = Height.edge_values[1,1]

        hc_n = Height.centroid_values[3]
        zc_n = Bed.centroid_values[3]

        qr[0] = Stage.edge_values[3,0]
        qr[1] = Xmom.edge_values[3,0]
        qr[2] = Ymom.edge_values[3,0]
        zr  = Bed.edge_values[3,0]
        hr = hre = Height.edge_values[3,0]

        z_half = max(zl, zr)
        hl = max(hle+zl-z_half,0.0)
        hr = max(hre+zr-z_half,0.0)

        max_speed1, pressure_flux1 = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux1, epsilon, z_half, g, H0, hc, hc_n, low_froude)

        bedslope_work1 = -g*0.5*(hl*hl - hle*hle -(hle+hc)*(zl-zc))+pressure_flux1
        edgeflux1[1] += normal[0]*bedslope_work1
        edgeflux1[2] += normal[1]*bedslope_work1

        # Flux across lower left hypotenuse of volume 1
        normal = domain.get_normal(1, 2)
        # ql = domain.get_conserved_quantities(vol_id=1, edge=2)
        # qr = domain.get_conserved_quantities(vol_id=0, edge=1)

        # hl = hle = ql[0] - zl
        # hr = hre = qr[0] - zr
        # hc = hc_n = (hl + hr) / 2

        ql[0] = Stage.edge_values[1,2]
        ql[1] = Xmom.edge_values[1,2]
        ql[2] = Ymom.edge_values[1,2]
        zl  = Bed.edge_values[1,2]
        hl = hle = Height.edge_values[1,2]

        hc_n = Height.centroid_values[0]
        zc_n = Bed.centroid_values[0]

        qr[0] = Stage.edge_values[0,1]
        qr[1] = Xmom.edge_values[0,1]
        qr[2] = Ymom.edge_values[0,1]
        zr  = Bed.edge_values[0,1]
        hr = hre = Height.edge_values[0,1]

        z_half = max(zl, zr)
        hl = max(hle+zl-z_half,0.0)
        hr = max(hre+zr-z_half,0.0)


        max_speed2, pressure_flux2 = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux2, epsilon, z_half, g, H0, hc, hc_n, low_froude)

        bedslope_work2 = -g*0.5*(hl*hl - hle*hle -(hle+hc)*(zl-zc))+pressure_flux2
        edgeflux2[1] += normal[0]*bedslope_work2
        edgeflux2[2] += normal[1]*bedslope_work2

        # Scale, add up and check that compute_fluxes is correct for vol 1
        # FIXME (Ole): This no longer works after the introduction of pressure_flux
        e0 = domain.edgelengths[1, 0]
        e1 = domain.edgelengths[1, 1]
        e2 = domain.edgelengths[1, 2]

        total_flux = -(e0*edgeflux0 +
                       e1*edgeflux1 +
                       e2*edgeflux2) / domain.areas[1]


        domain.compute_fluxes()

        for i, name in enumerate(['stage', 'xmomentum', 'ymomentum']):
            #                    domain.quantities[name].explicit_update[1])
            assert num.allclose(total_flux[i],
                                domain.quantities[name].explicit_update[1])

    # FIXME (Ole): Need test like this for fluxes in very shallow water.


    def test_compute_fluxes_structure_3(self):
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
        ze = (zl + zr) / 2

        domain.set_quantity('elevation', zl*num.ones((4, 3), float)) #array default#

        edgeflux = num.zeros(3, float)
        edgeflux0 = num.zeros(3, float)
        edgeflux1 = num.zeros(3, float)
        edgeflux2 = num.zeros(3, float)
        ql = num.zeros(3, float)
        qr = num.zeros(3, float)

        H0 = 0.0
        low_froude = 0

        stage_vertex = num.array([[val0, val0-1, val0-2],
                                      [val1, val1+1, val1],
                                      [val2, val2-2, val2],
                                      [val3-0.5, val3, val3]])

        domain.set_quantity('stage', stage_vertex )

        domain.set_quantity('height', stage_vertex + zl)

        domain.set_quantity('xmomentum',
                            [[1,2,3], [3,4,5], [1,-1,0], [0,-2,2]])

        domain.set_quantity('ymomentum',
                            [[1,-1,0], [0,-3,2], [0,1,0], [-1,2,2]])

        domain.check_integrity()

        # # Flux across right edge of volume 1
        # normal = domain.get_normal(1, 0)


        # ql = domain.get_conserved_quantities(vol_id=1, edge=0)
        # qr = domain.get_conserved_quantities(vol_id=2, edge=2)

        # hl = hle = ql[0] - zl
        # hr = hre = qr[1] - zr
        # hc = hc_n = (hl + hr) / 2
        # max_speed0, pressure_flux0 = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux0, epsilon, ze, g, H0, hc, hc_n, low_froude)

        # # Flux across upper edge of volume 1
        # normal = domain.get_normal(1, 1)
        # ql = domain.get_conserved_quantities(vol_id=1, edge=1)
        # qr = domain.get_conserved_quantities(vol_id=3, edge=0)

        # hl = hle = ql[0] - zl
        # hr = hre = qr[1] - zr
        # hc = hc_n = (hl + hr) / 2
        # max_speed1, pressure_flux1 = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux1, epsilon, ze, g, H0, hc, hc_n, low_froude)

        # # Flux across lower left hypotenuse of volume 1
        # normal = domain.get_normal(1, 2)
        # ql = domain.get_conserved_quantities(vol_id=1, edge=2)
        # qr = domain.get_conserved_quantities(vol_id=0, edge=1)

        # hl = hle = ql[0] - zl
        # hr = hre = qr[1] - zr
        # hc = hc_n = (hl + hr) / 2
        # max_speed2, pressure_flux2 = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux2, epsilon, ze, g, H0, hc, hc_n, low_froude)

        Stage  = domain.quantities['stage']
        Height = domain.quantities['height']
        Bed    = domain.quantities['elevation']
        Xmom   = domain.quantities['xmomentum']
        Ymom   = domain.quantities['ymomentum']

        hc = Height.centroid_values[1]
        zc = Bed.centroid_values[1]


        # Flux across right edge of volume 1
        normal = domain.get_normal(1, 0)
        # ql = domain.get_conserved_quantities(vol_id=1, edge=0)
        # qr = domain.get_conserved_quantities(vol_id=2, edge=2)

        # hl = hle = ql[0] - zl
        # hr = hre = qr[0] - zr
        # hc = hc_n = (hl + hr) / 2


        ql[0] = Stage.edge_values[1,0]
        ql[1] = Xmom.edge_values[1,0]
        ql[2] = Ymom.edge_values[1,0]
        zl  = Bed.edge_values[1,0]
        hle = Height.edge_values[1,0]

        hc_n = Height.centroid_values[2]
        zc_n = Bed.centroid_values[2]

        qr[0] = Stage.edge_values[2,2]
        qr[1] = Xmom.edge_values[2,2]
        qr[2] = Ymom.edge_values[2,2]
        zr  = Bed.edge_values[2,2]
        hre = Height.edge_values[2,2]

        z_half = max(zl, zr)

        hl = max(hle+zl-z_half,0.0)
        hr = max(hre+zr-z_half,0.0)

        low_froude = 0

        max_speed0, pressure_flux0 = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux0, epsilon, z_half, g, H0, hc, hc_n, low_froude)

        bedslope_work0 = -g*0.5*(hl*hl - hle*hle -(hle+hc)*(zl-zc))+pressure_flux0
        edgeflux0[1] += normal[0]*bedslope_work0
        edgeflux0[2] += normal[1]*bedslope_work0


        # Flux across upper edge of volume 1
        normal = domain.get_normal(1, 1)
        # ql = domain.get_conserved_quantities(vol_id=1, edge=1)
        # qr = domain.get_conserved_quantities(vol_id=3, edge=0)

        # hl = hle = ql[0] - zl
        # hr = hre = qr[0] - zr
        # hc = hc_n = (hl + hr) / 2


        ql[0] = Stage.edge_values[1,1]
        ql[1] = Xmom.edge_values[1,1]
        ql[2] = Ymom.edge_values[1,1]
        zl  = Bed.edge_values[1,1]
        hl = hle = Height.edge_values[1,1]

        hc_n = Height.centroid_values[3]
        zc_n = Bed.centroid_values[3]

        qr[0] = Stage.edge_values[3,0]
        qr[1] = Xmom.edge_values[3,0]
        qr[2] = Ymom.edge_values[3,0]
        zr  = Bed.edge_values[3,0]
        hr = hre = Height.edge_values[3,0]

        z_half = max(zl, zr)
        hl = max(hle+zl-z_half,0.0)
        hr = max(hre+zr-z_half,0.0)

        max_speed1, pressure_flux1 = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux1, epsilon, z_half, g, H0, hc, hc_n, low_froude)

        bedslope_work1 = -g*0.5*(hl*hl - hle*hle -(hle+hc)*(zl-zc))+pressure_flux1
        edgeflux1[1] += normal[0]*bedslope_work1
        edgeflux1[2] += normal[1]*bedslope_work1

        # Flux across lower left hypotenuse of volume 1
        normal = domain.get_normal(1, 2)
        # ql = domain.get_conserved_quantities(vol_id=1, edge=2)
        # qr = domain.get_conserved_quantities(vol_id=0, edge=1)

        # hl = hle = ql[0] - zl
        # hr = hre = qr[0] - zr
        # hc = hc_n = (hl + hr) / 2

        ql[0] = Stage.edge_values[1,2]
        ql[1] = Xmom.edge_values[1,2]
        ql[2] = Ymom.edge_values[1,2]
        zl  = Bed.edge_values[1,2]
        hl = hle = Height.edge_values[1,2]

        hc_n = Height.centroid_values[0]
        zc_n = Bed.centroid_values[0]

        qr[0] = Stage.edge_values[0,1]
        qr[1] = Xmom.edge_values[0,1]
        qr[2] = Ymom.edge_values[0,1]
        zr  = Bed.edge_values[0,1]
        hr = hre = Height.edge_values[0,1]

        z_half = max(zl, zr)
        hl = max(hle+zl-z_half,0.0)
        hr = max(hre+zr-z_half,0.0)


        max_speed2, pressure_flux2 = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux2, epsilon, z_half, g, H0, hc, hc_n, low_froude)

        bedslope_work2 = -g*0.5*(hl*hl - hle*hle -(hle+hc)*(zl-zc))+pressure_flux2
        edgeflux2[1] += normal[0]*bedslope_work2
        edgeflux2[2] += normal[1]*bedslope_work2

        # Scale, add up and check that compute_fluxes is correct for vol 1
        e0 = domain.edgelengths[1, 0]
        e1 = domain.edgelengths[1, 1]
        e2 = domain.edgelengths[1, 2]

        total_flux = -(e0*edgeflux0 +
                       e1*edgeflux1 +
                       e2*edgeflux2) / domain.areas[1]


        # Now check that compute_flux yields the same
        domain.compute_fluxes()

        for i, name in enumerate(['stage', 'xmomentum', 'ymomentum']):
            assert num.allclose(total_flux[i],
                                domain.quantities[name].explicit_update[1])


    def test_compute_fluxes_default_1(self):
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
        domain.set_flow_algorithm('DE0')

        val0 = 2. + 2.0/3
        val1 = 4. + 4.0/3
        val2 = 8. + 2.0/3
        val3 = 2. + 8.0/3

        domain.set_quantity('stage', [[val0, val0, val0], [val1, val1, val1],
                                      [val2, val2, val2], [val3, val3, val3]])
        domain.check_integrity()

        zl = zr = ze = 0.    # Assume flat bed

        edgeflux = num.zeros(3, float)
        edgeflux0 = num.zeros(3, float)
        edgeflux1 = num.zeros(3, float)
        edgeflux2 = num.zeros(3, float)
        H0 = 0.0

        # Flux across right edge of volume 1
        normal = domain.get_normal(1, 0)    # Get normal 0 of triangle 1
        assert num.allclose(normal, [1, 0])

        ql = domain.get_conserved_quantities(vol_id=1, edge=0)
        assert num.allclose(ql, [val1, 0, 0])

        qr = domain.get_conserved_quantities(vol_id=2, edge=2)
        assert num.allclose(qr, [val2, 0, 0])

        hl = hle = val1 - zl
        hr = hre = val2 - zr
        hc = hc_n = (hl + hr) / 2
        low_froude = 1
        domain.set_quantity('height', hc)

        max_speed, pressure_flux = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux0, epsilon, ze, g, H0, hc, hc_n, low_froude)

        assert num.allclose(max_speed, 9.21592824046)
        assert num.allclose(pressure_flux, 253.71111)
        # Flux across edge in the east direction (as per normal vector)
        assert num.allclose(edgeflux0, [-15.3598804, 0, 0.])


        # Flux across edge in the west direction (opposite sign for xmomentum)
        normal_opposite = domain.get_normal(2, 2)   # Get normal 2 of triangle 2
        assert num.allclose(normal_opposite, [-1, 0])

        max_speed, pressure_flux = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux, epsilon, ze, g, H0, hc, hc_n, low_froude)


        assert num.allclose(max_speed, 9.21592824046)
        assert num.allclose(pressure_flux, 253.71111)
        assert num.allclose(edgeflux, [-15.3598804, 0, 0.])

        # Flux across upper edge of volume 1
        normal = domain.get_normal(1, 1)
        ql = domain.get_conserved_quantities(vol_id=1, edge=1)
        qr = domain.get_conserved_quantities(vol_id=3, edge=0)

        hl = hle = ql[0] - zl
        hr = hre = qr[0] - zr
        hc = hc_n = (hl + hr) / 2

        max_speed, pressure_flux = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux1, epsilon, ze, g, H0, hc, hc_n, low_froude)

        assert num.allclose(max_speed, 7.22956891292)
        assert num.allclose(pressure_flux, 123.04444444)
        assert num.allclose(edgeflux1, [2.4098563, 0., 0.])


        # Flux across lower left hypotenuse of volume 1
        normal = domain.get_normal(1, 2)
        ql = domain.get_conserved_quantities(vol_id=1, edge=2)
        qr = domain.get_conserved_quantities(vol_id=0, edge=1)

        hl = hle = ql[0] - zl
        hr = hre = qr[0] - zr
        hc = hc_n = (hl + hr) / 2

        max_speed, pressure_flux = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux2, epsilon, ze, g, H0, hc, hc_n, low_froude)

        assert num.allclose(max_speed, 7.22956891292)
        assert num.allclose(pressure_flux, 87.111111)
        assert num.allclose(edgeflux2, [9.63942522, 0., 0.])

        # Scale, add up and check that compute_fluxes is correct for vol 1
        e0 = domain.edgelengths[1, 0]
        e1 = domain.edgelengths[1, 1]
        e2 = domain.edgelengths[1, 2]

        total_flux = -(e0*edgeflux0 +
                       e1*edgeflux1 +
                       e2*edgeflux2) / domain.areas[1]

        assert num.allclose(total_flux, [-0.68218178, 0., 0.])

        domain.compute_fluxes()

        # FIXME (Ole): Why does flux in 'stage' work, but flux in xmom and ymom have been lumped into pressure_flux?
        # Doesn't that affect direction?
        # Can this be tested?
        #for i, name in enumerate(['stage', 'xmomentum', 'ymomentum']):
        #    assert num.allclose(total_flux[i],
        #                        domain.quantities[name].explicit_update[1])

        msg = 'Got %s' % (str(domain.quantities['stage'].explicit_update))
        assert num.allclose(domain.quantities['stage'].explicit_update,
                            [-6.46904403, -4.5743049, -100.45244512, -43.89591759]), msg
        assert num.allclose(domain.quantities['xmomentum'].explicit_update,
                            [-120.05, 0., 120.05, 0.])
        assert num.allclose(domain.quantities['ymomentum'].explicit_update,
                            [-120.05, 0., 0., 120.05])


    def test_compute_fluxes_DE_1(self):
        # This is a reuse of test_compute_fluxes_old_2 which used the original (now deprecated algorithm).
        # We changed the algorithm to DE and used those values to test it.
        # FIXME (Ole): Need to carefully check how this test can be brought back. It is a very good internal integrity
        # test - but currently not working.
        #

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

        domain.set_compute_fluxes_method('DE')

        val0 = 2. + 2.0/3
        val1 = 4. + 4.0/3
        val2 = 8. + 2.0/3
        val3 = 2. + 8.0/3

        ze = zl = zr = zc = 0    # Assume flat zero bed
        edgeflux = num.zeros(3, float)
        edgeflux0 = num.zeros(3, float)
        edgeflux1 = num.zeros(3, float)
        edgeflux2 = num.zeros(3, float)

        ql = num.zeros(3, float)
        qr = num.zeros(3, float)

        H0 = domain.H0

        domain.set_quantity('elevation', zl*num.ones((4, 3), int)) #array default#

        domain.set_quantity('stage', [[val0, val0-1, val0-2],
                                      [val1, val1+1, val1],
                                      [val2, val2-2, val2],
                                      [val3-0.5, val3, val3]])

        domain.set_quantity('height', [[val0 - zl, val0-1 - zl, val0-2 - zl],
                                      [val1 - zl, val1+1 - zl, val1 - zl],
                                      [val2 - zl, val2-2 - zl, val2 - zl],
                                      [val3-0.5 - zl, val3 - zl, val3 - zl]])

        domain.set_quantity('xmomentum',
                            [[1,2,3], [3,4,5], [1,-1,0], [0,-2,2]])

        domain.set_quantity('ymomentum',
                            [[1,-1,0], [0,-3,2], [0,1,0], [-1,2,2]])

        domain.check_integrity()

        Stage  = domain.quantities['stage']
        Height = domain.quantities['height']
        Bed    = domain.quantities['elevation']
        Xmom   = domain.quantities['xmomentum']
        Ymom   = domain.quantities['ymomentum']

        hc = Height.centroid_values[1]
        zc = Bed.centroid_values[1]

        l0 = domain.edgelengths[1, 0]
        l1 = domain.edgelengths[1, 1]
        l2 = domain.edgelengths[1, 2]

        # Flux across right edge of volume 1 (volume 2)
        normal = domain.get_normal(1, 0)


        ql[0] = Stage.edge_values[1,0]
        ql[1] = Xmom.edge_values[1,0]
        ql[2] = Ymom.edge_values[1,0]
        zl  = Bed.edge_values[1,0]
        hle = Height.edge_values[1,0]

        hc_n = Height.centroid_values[2]
        zc_n = Bed.centroid_values[2]

        qr[0] = Stage.edge_values[2,2]
        qr[1] = Xmom.edge_values[2,2]
        qr[2] = Ymom.edge_values[2,2]
        zr  = Bed.edge_values[2,2]
        hre = Height.edge_values[2,2]

        z_half = max(zl, zr)

        hl = max(hle+zl-z_half,0.0)
        hr = max(hre+zr-z_half,0.0)

        low_froude = 0
        #domain.set_quantity('height', hc)

        max_speed0, pressure_flux0 = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux0, epsilon, z_half, g, H0, hc, hc_n, low_froude)

        bedslope_work0 = -g*0.5*(hl*hl - hle*hle -(hle+hc)*(zl-zc))+pressure_flux0
        edgeflux0[1] += normal[0]*bedslope_work0
        edgeflux0[2] += normal[1]*bedslope_work0
        assert(num.any(edgeflux0 != 0))


        # Flux across upper edge of volume 1 (volume 3)
        normal = domain.get_normal(1, 1)

        ql[0] = Stage.edge_values[1,1]
        ql[1] = Xmom.edge_values[1,1]
        ql[2] = Ymom.edge_values[1,1]
        zl  = Bed.edge_values[1,1]
        hl = hle = Height.edge_values[1,1]

        hc_n = Height.centroid_values[3]
        zc_n = Bed.centroid_values[3]

        qr[0] = Stage.edge_values[3,0]
        qr[1] = Xmom.edge_values[3,0]
        qr[2] = Ymom.edge_values[3,0]
        zr  = Bed.edge_values[3,0]
        hr = hre = Height.edge_values[3,0]

        z_half = max(zl, zr)
        hl = max(hle+zl-z_half,0.0)
        hr = max(hre+zr-z_half,0.0)

        #hc = hc_n = (hl + hr)/2
        #domain.set_quantity('height', hc)



        max_speed, pressure_flux1 = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux1, epsilon, z_half, g, H0, hc, hc_n, low_froude)

        bedslope_work1 = -g*0.5*(hl*hl - hle*hle -(hle+hc)*(zl-zc))+pressure_flux1
        edgeflux1[1] += normal[0]*bedslope_work1
        edgeflux1[2] += normal[1]*bedslope_work1
        assert(num.any(edgeflux1 != 0))

        # Flux across lower left hypotenuse of volume 1 (volume 0)
        normal = domain.get_normal(1, 2)

        ql[0] = Stage.edge_values[1,2]
        ql[1] = Xmom.edge_values[1,2]
        ql[2] = Ymom.edge_values[1,2]
        zl  = Bed.edge_values[1,2]
        hl = hle = Height.edge_values[1,2]

        hc_n = Height.centroid_values[0]
        zc_n = Bed.centroid_values[0]

        qr[0] = Stage.edge_values[0,1]
        qr[1] = Xmom.edge_values[0,1]
        qr[2] = Ymom.edge_values[0,1]
        zr  = Bed.edge_values[0,1]
        hr = hre = Height.edge_values[0,1]

        z_half = max(zl, zr)
        hl = max(hle+zl-z_half,0.0)
        hr = max(hre+zr-z_half,0.0)

        max_speed, pressure_flux2 = flux_function(normal, ql, qr, hl, hr, hle, hre, edgeflux2, epsilon, z_half, g, H0, hc, hc_n, low_froude)

        bedslope_work2 = -g*0.5*(hl*hl - hle*hle -(hle+hc)*(zl-zc))+pressure_flux2
        edgeflux2[1] += normal[0]*bedslope_work2
        edgeflux2[2] += normal[1]*bedslope_work2
        assert(num.any(edgeflux2 != 0))

        # Scale, add up and check that compute_fluxes is correct for vol 1
        # FIXME (Ole): This does not work after the introduction of pressure_flux

        total_flux = -(l0*edgeflux0 +
                       l1*edgeflux1 +
                       l2*edgeflux2) / domain.areas[1]

        domain.compute_fluxes()

        for i, name in enumerate(['stage', 'xmomentum', 'ymomentum']):
            msg = 'Expected %f for %s but got %f' % ((domain.quantities[name].explicit_update[1]),
                                                     name,
                                                     total_flux[i])
            assert num.allclose(total_flux[i],
                                domain.quantities[name].explicit_update[1]), msg



if __name__ == "__main__":
    unittest.main()
