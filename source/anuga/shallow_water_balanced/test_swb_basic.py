#!/usr/bin/env python

import unittest, os
import os.path
from math import pi, sqrt
import tempfile

import anuga

from anuga.config import g, epsilon
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a
from anuga.utilities.numerical_tools import mean
from anuga.geometry.polygon import is_inside_polygon
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.abstract_2d_finite_volumes.quantity import Quantity
from anuga.geospatial_data.geospatial_data import Geospatial_data
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross

from anuga.utilities.system_tools import get_pathname_from_package
from swb_domain import *

import numpy as num

# Get gateway to C implementation of flux function for direct testing
from shallow_water_ext import flux_function_central as flux_function


# For test_fitting_using_shallow_water_domain example
def linear_function(point):
    point = num.array(point)
    return point[:,0]+point[:,1]

class Weir:
    """Set a bathymetry for weir with a hole and a downstream gutter
    x,y are assumed to be in the unit square
    """

    def __init__(self, stage):
        self.inflow_stage = stage

    def __call__(self, x, y):
        N = len(x)
        assert N == len(y)

        z = num.zeros(N, num.float)
        for i in range(N):
            z[i] = -x[i]/2  #General slope

            #Flattish bit to the left
            if x[i] < 0.3:
                z[i] = -x[i]/10

            #Weir
            if x[i] >= 0.3 and x[i] < 0.4:
                z[i] = -x[i]+0.9

            #Dip
            x0 = 0.6
            depth = -1.0
            plateaux = -0.6
            if y[i] < 0.7:
                if x[i] > x0 and x[i] < 0.9:
                    z[i] = depth
                #RHS plateaux
                if x[i] >= 0.9:
                    z[i] = plateaux
            elif y[i] >= 0.7 and y[i] < 1.5:
                #Restrict and deepen
                if x[i] >= x0 and x[i] < 0.8:
                    z[i] = depth - (y[i]/3 - 0.3)
                elif x[i] >= 0.8:
                    #RHS plateaux
                    z[i] = plateaux
            elif y[i] >= 1.5:
                if x[i] >= x0 and x[i] < 0.8 + (y[i]-1.5)/1.2:
                    #Widen up and stay at constant depth
                    z[i] = depth-1.5/5
                elif x[i] >= 0.8 + (y[i]-1.5)/1.2:
                    #RHS plateaux
                    z[i] = plateaux

            #Hole in weir (slightly higher than inflow condition)
            if x[i] >= 0.3 and x[i] < 0.4 and y[i] > 0.2 and y[i] < 0.4:
                z[i] = -x[i]+self.inflow_stage + 0.02

            #Channel behind weir
            x0 = 0.5
            if x[i] >= 0.4 and x[i] < x0 and y[i] > 0.2 and y[i] < 0.4:
                z[i] = -x[i]+self.inflow_stage + 0.02

            if x[i] >= x0 and x[i] < 0.6 and y[i] > 0.2 and y[i] < 0.4:
                #Flatten it out towards the end
                z[i] = -x0+self.inflow_stage + 0.02 + (x0-x[i])/5

            # Hole to the east
            x0 = 1.1
            y0 = 0.35
            if num.sqrt((2*(x[i]-x0))**2 + (2*(y[i]-y0))**2) < 0.2:
                z[i] = num.sqrt(((x[i]-x0))**2 + ((y[i]-y0))**2)-1.0

            #Tiny channel draining hole
            if x[i] >= 1.14 and x[i] < 1.2 and y[i] >= 0.4 and y[i] < 0.6:
                z[i] = -0.9 #North south

            if x[i] >= 0.9 and x[i] < 1.18 and y[i] >= 0.58 and y[i] < 0.65:
                z[i] = -1.0 + (x[i]-0.9)/3 #East west

            # Stuff not in use

            # Upward slope at inlet to the north west
            # if x[i] < 0.0: # and y[i] > 0.5:
            #    #z[i] = -y[i]+0.5  #-x[i]/2
            #    z[i] = x[i]/4 - y[i]**2 + 0.5

            # Hole to the west
            # x0 = -0.4; y0 = 0.35 # center
            # if sqrt((2*(x[i]-x0))**2 + (2*(y[i]-y0))**2) < 0.2:
            #    z[i] = sqrt(((x[i]-x0))**2 + ((y[i]-y0))**2)-0.2

        return z/2



#########################################################

class Test_swb_basic(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


    def test_rotate(self):
        from anuga.shallow_water.shallow_water_ext import rotate
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
        except:
            pass
        else:
            raise Exception('Should have raised an exception')

    # Individual flux tests
    def test_flux_zero_case(self):
        ql = num.zeros(3, num.float)
        qr = num.zeros(3, num.float)
        normal = num.zeros(2, num.float)
        edgeflux = num.zeros(3, num.float)
        zl = zr = 0.
        H0 = 1.0e-3 # As suggested in the manual
        
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux, epsilon, g, H0)

        assert num.allclose(edgeflux, [0,0,0])
        assert max_speed == 0.

    def test_flux_constants(self):
        w = 2.0

        normal = num.array([1.,0])
        ql = num.array([w, 0, 0])
        qr = num.array([w, 0, 0])
        edgeflux = num.zeros(3, num.float)
        zl = zr = 0.
        h = w - (zl+zr)/2
        H0 = 0.0

        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux, epsilon, g, H0)        
        assert num.allclose(edgeflux, [0., 0.5*g*h**2, 0.])
        assert max_speed == num.sqrt(g*h)

    #def test_flux_slope(self):
    #    #FIXME: TODO
    #    w = 2.0
    #
    #    normal = array([1.,0])
    #    ql = array([w, 0, 0])
    #    qr = array([w, 0, 0])
    #    zl = zr = 0.
    #    h = w - (zl+zr)/2
    #
    #    flux, max_speed = flux_function(normal, ql, qr, zl, zr)
    #
    #    assert allclose(flux, [0., 0.5*g*h**2, 0.])
    #    assert max_speed == sqrt(g*h)

    def test_flux1(self):
        # Use data from previous version of abstract_2d_finite_volumes
        normal = num.array([1., 0])
        ql = num.array([-0.2, 2, 3])
        qr = num.array([-0.2, 2, 3])
        zl = zr = -0.5
        edgeflux = num.zeros(3, num.float)

        H0 = 0.0

        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux,
                                  epsilon, g, H0)

        assert num.allclose(edgeflux, [2., 13.77433333, 20.])
        assert num.allclose(max_speed, 8.38130948661)

    def test_flux2(self):
        # Use data from previous version of abstract_2d_finite_volumes
        normal = num.array([0., -1.])
        ql = num.array([-0.075, 2, 3])
        qr = num.array([-0.075, 2, 3])
        zl = zr = -0.375

        edgeflux = num.zeros(3, num.float)
        H0 = 0.0
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux,
                                  epsilon, g, H0)

        assert num.allclose(edgeflux, [-3., -20.0, -30.441])
        assert num.allclose(max_speed, 11.7146428199)

    def test_flux3(self):
        # Use data from previous version of abstract_2d_finite_volumes
        normal = num.array([-sqrt(2)/2, sqrt(2)/2])
        ql = num.array([-0.075, 2, 3])
        qr = num.array([-0.075, 2, 3])
        zl = zr = -0.375

        edgeflux = num.zeros(3, num.float)
        H0 = 0.0
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux,
                                  epsilon, g, H0)

        assert num.allclose(edgeflux, [sqrt(2)/2, 4.40221112, 7.3829019])
        assert num.allclose(max_speed, 4.0716654239)

    def test_flux4(self):
        # Use data from previous version of abstract_2d_finite_volumes
        normal = num.array([-sqrt(2)/2, sqrt(2)/2])
        ql = num.array([-0.34319278, 0.10254161, 0.07273855])
        qr = num.array([-0.30683287, 0.1071986, 0.05930515])
        zl = zr = -0.375

        edgeflux = num.zeros(3, num.float)
        H0 = 0.0
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux,
                                  epsilon, g, H0)

        assert num.allclose(edgeflux, [-0.04072676, -0.07096636, -0.01604364])
        assert num.allclose(max_speed, 1.31414103233)

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
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

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
            return -x/2

        domain.set_quantity('elevation', -10)
        domain.set_quantity('stage', surface)
        domain.set_quantity('xmomentum', 1)

        domain.compute_fluxes()

        #print domain.get_quantity('stage').explicit_update
        # FIXME (Ole): TODO the general case
        #assert allclose(domain.get_quantity('stage').explicit_update[1], ...??)

    def test_sw_domain_simple(self):
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #             bac,     bce,     ecf,     dbe
        vertices = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        #from anuga.abstract_2d_finite_volumes.domain import Domain as Generic_domain
        #msg = 'The class %s is not a subclass of the generic domain class %s'\
        #      %(DomainClass, Domain)
        #assert issubclass(DomainClass, Domain), msg

        domain = Domain(points, vertices)
        domain.check_integrity()

        for name in ['stage', 'xmomentum', 'ymomentum',
                     'elevation', 'friction']:
            assert domain.quantities.has_key(name)

        assert num.alltrue(domain.get_conserved_quantities(0, edge=1) == 0.)

    def xtest_catching_negative_heights(self):
        #OBSOLETE

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

        zl = zr = 4    # Too large
        domain.set_quantity('elevation', zl*num.ones((4, 3), num.int)) #array default#
        domain.set_quantity('stage', [[val0, val0-1, val0-2],
                                      [val1, val1+1, val1],
                                      [val2, val2-2, val2],
                                      [val3-0.5, val3, val3]])

        #Should fail
        try:
            domain.check_integrity()
        except:
            pass

    def test_get_wet_elements(self):
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

        zl = zr = 5
        domain.set_quantity('elevation', zl*num.ones((4, 3), num.int)) #array default#
        domain.set_quantity('stage', [[val0, val0-1, val0-2],
                                      [val1, val1+1, val1],
                                      [val2, val2-2, val2],
                                      [val3-0.5, val3, val3]])

        domain.check_integrity()

        indices = domain.get_wet_elements()
        assert num.allclose(indices, [1, 2])

        indices = domain.get_wet_elements(indices=[0, 1, 3])
        assert num.allclose(indices, [1])

    def test_get_maximum_inundation_1(self):
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

        domain.set_quantity('elevation', lambda x, y: x+2*y)    # 2 4 4 6
        domain.set_quantity('stage', 3)

        domain.check_integrity()

        indices = domain.get_wet_elements()
        assert num.allclose(indices, [0])

        q = domain.get_maximum_inundation_elevation()
        assert num.allclose(q, domain.get_quantity('elevation').\
                                   get_values(location='centroids')[0])

        x, y = domain.get_maximum_inundation_location()
        assert num.allclose([x, y], domain.get_centroid_coordinates()[0])

    def test_get_maximum_inundation_2(self):
        """test_get_maximum_inundation_2(self)

        Test multiple wet cells with same elevation
        """

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

        domain.set_quantity('elevation', lambda x, y: x+2*y)    # 2 4 4 6
        domain.set_quantity('stage', 4.1)

        domain.check_integrity()

        indices = domain.get_wet_elements()
        assert num.allclose(indices, [0, 1, 2])

        q = domain.get_maximum_inundation_elevation()
        assert num.allclose(q, 4)

        x, y = domain.get_maximum_inundation_location()
        assert num.allclose([x, y], domain.get_centroid_coordinates()[1])

    def test_get_maximum_inundation_3(self):
        """test_get_maximum_inundation_3(self)

        Test of real runup example:
        """

        from anuga.abstract_2d_finite_volumes.mesh_factory \
                import rectangular_cross

        initial_runup_height = -0.4
        final_runup_height = -0.3

        #--------------------------------------------------------------
        # Setup computational domain
        #--------------------------------------------------------------
        N = 5
        points, vertices, boundary = rectangular_cross(N, N)
        domain = Domain(points, vertices, boundary)
        domain.set_maximum_allowed_speed(1.0)

        #--------------------------------------------------------------
        # Setup initial conditions
        #--------------------------------------------------------------
        def topography(x, y):
            return -x/2                             # linear bed slope

        # Use function for elevation
        domain.set_quantity('elevation', topography)
        domain.set_quantity('friction', 0.)                # Zero friction
        # Constant negative initial stage
        domain.set_quantity('stage', initial_runup_height)

        #--------------------------------------------------------------
        # Setup boundary conditions
        #--------------------------------------------------------------
        Br = anuga.Reflective_boundary(domain)                 # Reflective wall
        Bd = anuga.Dirichlet_boundary([final_runup_height, 0, 0])# Steady inflow

        # All reflective to begin with (still water)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        #--------------------------------------------------------------
        # Test initial inundation height
        #--------------------------------------------------------------

        indices = domain.get_wet_elements()
        z = domain.get_quantity('elevation').get_values(location='centroids',
                                                        indices=indices)
        assert num.alltrue(z < initial_runup_height)

        q = domain.get_maximum_inundation_elevation()
        # First order accuracy
        assert num.allclose(q, initial_runup_height, rtol=1.0/N)

        x, y = domain.get_maximum_inundation_location()

        qref = domain.get_quantity('elevation').\
                     get_values(interpolation_points=[[x, y]])
        assert num.allclose(q, qref)

        wet_elements = domain.get_wet_elements()
        wet_elevations = domain.get_quantity('elevation').\
                                    get_values(location='centroids',
                                               indices=wet_elements)
        assert num.alltrue(wet_elevations < initial_runup_height)
        assert num.allclose(wet_elevations, z)

        #--------------------------------------------------------------
        # Let triangles adjust
        #--------------------------------------------------------------
        for t in domain.evolve(yieldstep = 0.1, finaltime = 1.0):
            pass

        #--------------------------------------------------------------
        # Test inundation height again
        #--------------------------------------------------------------
        indices = domain.get_wet_elements()
        z = domain.get_quantity('elevation').get_values(location='centroids',
                                                        indices=indices)

        assert num.alltrue(z < initial_runup_height)

        q = domain.get_maximum_inundation_elevation()
        # First order accuracy
        assert num.allclose(q, initial_runup_height, rtol=1.0/N)

        x, y = domain.get_maximum_inundation_location()
        qref = domain.get_quantity('elevation').\
                        get_values(interpolation_points=[[x, y]])
        assert num.allclose(q, qref)

        #--------------------------------------------------------------
        # Update boundary to allow inflow
        #--------------------------------------------------------------
        domain.set_boundary({'right': Bd})

        #--------------------------------------------------------------
        # Evolve system through time
        #--------------------------------------------------------------
        for t in domain.evolve(yieldstep = 0.1, finaltime = 3.0):
            pass

        #--------------------------------------------------------------
        # Test inundation height again
        #--------------------------------------------------------------
        indices = domain.get_wet_elements()
        z = domain.get_quantity('elevation').\
                    get_values(location='centroids', indices=indices)

        assert num.alltrue(z < final_runup_height)

        q = domain.get_maximum_inundation_elevation()
        # First order accuracy
        assert num.allclose(q, final_runup_height, rtol=1.0/N)

        x, y = domain.get_maximum_inundation_location()
        qref = domain.get_quantity('elevation').\
                        get_values(interpolation_points=[[x, y]])
        assert num.allclose(q, qref)

        wet_elements = domain.get_wet_elements()
        wet_elevations = domain.get_quantity('elevation').\
                             get_values(location='centroids',
                                        indices=wet_elements)
        assert num.alltrue(wet_elevations < final_runup_height)
        assert num.allclose(wet_elevations, z)

    def test_get_maximum_inundation_from_sww(self):
        """test_get_maximum_inundation_from_sww(self)

        Test of get_maximum_inundation_elevation()
        and get_maximum_inundation_location() from data_manager.py

        This is based on test_get_maximum_inundation_3(self) but works with the
        stored results instead of with the internal data structure.

        This test uses the underlying get_maximum_inundation_data for tests
        """

        from anuga.abstract_2d_finite_volumes.mesh_factory \
                import rectangular_cross
        from anuga.shallow_water.sww_interrogate import \
            get_maximum_inundation_elevation, get_maximum_inundation_location, \
            get_maximum_inundation_data

        initial_runup_height = -0.4
        final_runup_height = -0.3

        #--------------------------------------------------------------
        # Setup computational domain
        #--------------------------------------------------------------
        N = 10
        points, vertices, boundary = rectangular_cross(N, N)
        domain = Domain(points, vertices, boundary)
        domain.set_name('runup_test')
        #domain.set_maximum_allowed_speed(1.0)

        # FIXME: This works better with old limiters so far
        #domain.tight_slope_limiters = 0

        #--------------------------------------------------------------
        # Setup initial conditions
        #--------------------------------------------------------------
        def topography(x, y):
            return -x/2                             # linear bed slope

        # Use function for elevation
        domain.set_quantity('elevation', topography)
        domain.set_quantity('friction', 0.)                # Zero friction
        # Constant negative initial stage
        domain.set_quantity('stage', initial_runup_height)

        #--------------------------------------------------------------
        # Setup boundary conditions
        #--------------------------------------------------------------
        Br = anuga.Reflective_boundary(domain)                 # Reflective wall
        Bd = anuga.Dirichlet_boundary([final_runup_height, 0, 0])# Steady inflow

        # All reflective to begin with (still water)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        #--------------------------------------------------------------
        # Test initial inundation height
        #--------------------------------------------------------------
        indices = domain.get_wet_elements()
        z = domain.get_quantity('elevation').\
                get_values(location='centroids', indices=indices)
        assert num.alltrue(z < initial_runup_height)

        q_ref = domain.get_maximum_inundation_elevation()
        # First order accuracy
        assert num.allclose(q_ref, initial_runup_height, rtol=1.0/N)

        #--------------------------------------------------------------
        # Let triangles adjust
        #--------------------------------------------------------------
        for t in domain.evolve(yieldstep = 0.1, finaltime = 1.0):
            pass

        #--------------------------------------------------------------
        # Test inundation height again
        #--------------------------------------------------------------
        q_ref = domain.get_maximum_inundation_elevation()
        q = get_maximum_inundation_elevation('runup_test.sww')
        msg = 'We got %f, should have been %f' % (q, q_ref)
        assert num.allclose(q, q_ref, rtol=1.0/N), msg

        q = get_maximum_inundation_elevation('runup_test.sww')
        msg = 'We got %f, should have been %f' % (q, initial_runup_height)
        assert num.allclose(q, initial_runup_height, rtol = 1.0/N), msg

        # Test error condition if time interval is out
        try:
            q = get_maximum_inundation_elevation('runup_test.sww',
                                                 time_interval=[2.0, 3.0])
        except ValueError:
            pass
        else:
            msg = 'should have caught wrong time interval'
            raise Exception(msg)

        # Check correct time interval
        q, loc = get_maximum_inundation_data('runup_test.sww',
                                             time_interval=[0.0, 3.0])
        msg = 'We got %f, should have been %f' % (q, initial_runup_height)
        assert num.allclose(q, initial_runup_height, rtol = 1.0/N), msg
        assert num.allclose(-loc[0]/2, q)    # From topography formula

        #--------------------------------------------------------------
        # Update boundary to allow inflow
        #--------------------------------------------------------------
        domain.set_boundary({'right': Bd})

        #--------------------------------------------------------------
        # Evolve system through time
        #--------------------------------------------------------------
        q_max = None
        for t in domain.evolve(yieldstep = 0.1, finaltime = 3.0,
                               skip_initial_step = True):
            q = domain.get_maximum_inundation_elevation()
            if q > q_max:
                q_max = q

        #--------------------------------------------------------------
        # Test inundation height again
        #--------------------------------------------------------------
        indices = domain.get_wet_elements()
        z = domain.get_quantity('elevation').\
                get_values(location='centroids', indices=indices)

        assert num.alltrue(z < final_runup_height)

        q = domain.get_maximum_inundation_elevation()
        # First order accuracy
        assert num.allclose(q, final_runup_height, rtol=1.0/N)

        q, loc = get_maximum_inundation_data('runup_test.sww',
                                             time_interval=[3.0, 3.0])
        msg = 'We got %f, should have been %f' % (q, final_runup_height)
        assert num.allclose(q, final_runup_height, rtol=1.0/N), msg
        assert num.allclose(-loc[0]/2, q)    # From topography formula

        q = get_maximum_inundation_elevation('runup_test.sww')
        loc = get_maximum_inundation_location('runup_test.sww')
        msg = 'We got %f, should have been %f' % (q, q_max)
        assert num.allclose(q, q_max, rtol=1.0/N), msg
        assert num.allclose(-loc[0]/2, q)    # From topography formula

        q = get_maximum_inundation_elevation('runup_test.sww',
                                             time_interval=[0, 3])
        msg = 'We got %f, should have been %f' % (q, q_max)
        assert num.allclose(q, q_max, rtol=1.0/N), msg

        # Check polygon mode
        # Runup region
        polygon = [[0.3, 0.0], [0.9, 0.0], [0.9, 1.0], [0.3, 1.0]]
        q = get_maximum_inundation_elevation('runup_test.sww',
                                             polygon = polygon,
                                             time_interval=[0, 3])
        msg = 'We got %f, should have been %f' % (q, q_max)
        assert num.allclose(q, q_max, rtol=1.0/N), msg

        # Offshore region
        polygon = [[0.9, 0.0], [1.0, 0.0], [1.0, 1.0], [0.9, 1.0]]
        q, loc = get_maximum_inundation_data('runup_test.sww',
                                             polygon = polygon,
                                             time_interval=[0, 3])
        msg = 'We got %f, should have been %f' % (q, -0.475)
        assert num.allclose(q, -0.475, rtol=1.0/N), msg
        assert is_inside_polygon(loc, polygon)
        assert num.allclose(-loc[0]/2, q)    # From topography formula

        # Dry region
        polygon = [[0.0, 0.0], [0.2, 0.0], [0.2, 1.0], [0.0, 1.0]]
        q, loc = get_maximum_inundation_data('runup_test.sww',
                                             polygon = polygon,
                                             time_interval=[0, 3])
        msg = 'We got %s, should have been None' % (q)
        assert q is None, msg
        msg = 'We got %s, should have been None' % (loc)
        assert loc is None, msg

        # Check what happens if no time point is within interval
        try:
            q = get_maximum_inundation_elevation('runup_test.sww',
                                                 time_interval=[2.75, 2.75])
        except AssertionError:
            pass
        else:
            msg = 'Time interval should have raised an exception'
            raise Exception(msg)

        # Cleanup
        try:
            os.remove(domain.get_name() + '.sww')
        except:
            pass
            #FIXME(Ole): Windows won't allow removal of this





    def test_another_runup_example(self):
        """test_another_runup_example

        Test runup example where actual timeseries at interpolated
        points are tested.
        """

        #-----------------------------------------------------------------
        # Setup computational domain
        #-----------------------------------------------------------------
        points, vertices, boundary = anuga.rectangular_cross(10, 10)# Basic mesh
        domain = anuga.Domain(points, vertices, boundary) # Create domain
        domain.set_default_order(2)
        domain.set_quantities_to_be_stored(None)
        domain.H0 = 1.0e-3

        #-----------------------------------------------------------------
        # Setup initial conditions
        #-----------------------------------------------------------------
        def topography(x, y):
            return -x/2                              # linear bed slope

        domain.set_quantity('elevation', topography)
        domain.set_quantity('friction', 0.0)
        domain.set_quantity('stage', expression='elevation')

        #----------------------------------------------------------------
        # Setup boundary conditions
        #----------------------------------------------------------------
        Br = anuga.Reflective_boundary(domain)           # Solid reflective wall
        Bd = anuga.Dirichlet_boundary([-0.2, 0., 0.])    # Constant
        domain.set_boundary({'left': Br, 'right': Bd, 'top': Br, 'bottom': Br})

        #----------------------------------------------------------------
        # Evolve system through time
        #----------------------------------------------------------------
        interpolation_points = [[0.4,0.5], [0.6,0.5], [0.8,0.5], [0.9,0.5]]
        gauge_values = []
        for _ in interpolation_points:
            gauge_values.append([])

        time = []
        for t in domain.evolve(yieldstep=0.1, finaltime=5.0):
            # Record time series at known points
            time.append(domain.get_time())

            stage = domain.get_quantity('stage')
            w = stage.get_values(interpolation_points=interpolation_points)

            for i, _ in enumerate(interpolation_points):
                gauge_values[i].append(w[i])

        #Reference (nautilus 26/6/2008)
        
        G0 = [-0.20000000000000001, -0.20000000000000001, -0.20000000000000001, -0.1958465301767274, -0.19059602372752493, -0.18448466250700923, -0.16979321333876071, -0.15976372740651074, -0.1575649333345176, -0.15710373731900584, -0.1530922283220747, -0.18836084336565725, -0.19921529311644628, -0.19923451799698919, -0.19923795414410964, -0.19923178806924047, -0.19925157557666154, -0.19930747801697429, -0.1993266428576112, -0.19932004930281799, -0.19929691326931867, -0.19926285267313795, -0.19916645449780995, -0.1988942593296438, -0.19900620256621993, -0.19914327423060865, -0.19918708440899577, -0.19921557252449132, -0.1992404368022069, -0.19925070370697717, -0.19925975477892274, -0.1992671090445659, -0.19927254203777162, -0.19927631910959256, -0.19927843552031504, -0.19927880339239365, -0.19927763204453783, -0.19927545249577633, -0.19927289590622824, -0.19927076261495152, -0.19926974334736983, -0.19927002562760332, -0.19927138236272213, -0.1992734501064522, -0.19927573251318192, -0.19927778936001547, -0.1992793776883893, -0.19928040577720926, -0.19928092586206753, -0.19928110982948721, -0.19928118887248453]

        G1 = [-0.29999999999999993, -0.29999999999999993, -0.29139135018319512, -0.27257394456094503, -0.24141437432643265, -0.22089173942479151, -0.20796171092975532, -0.19874580192293825, -0.19014580508752857, -0.18421165368665365, -0.18020808282748838, -0.17518824759550247, -0.16436633464497749, -0.18714479115225544, -0.2045242886738807, -0.21011244240826329, -0.21151316017424124, -0.21048112933621732, -0.20772920477355789, -0.20489184334204144, -0.20286043930678221, -0.20094305756540246, -0.19948172752345467, -0.19886725178309209, -0.1986680808256765, -0.19860718133373548, -0.19862076543539733, -0.19888949069732539, -0.19932190310819023, -0.19982482967777454, -0.20036045468470615, -0.20064263130562704, -0.2007255389410077, -0.20068815669152493, -0.20057471332984647, -0.20042203940851802, -0.20026779013499779, -0.20015995671464712, -0.2000684005446525, -0.20001486753189174, -0.20000743467898013, -0.20003739771775905, -0.20008784600912621, -0.20013758305093884, -0.20017277554845025, -0.20018629233766885, -0.20018106288462198, -0.20016648079299326, -0.20015155958426531, -0.20014259747382254, -0.20014096648362509]
        
        
        G2 = [-0.40000000000000002, -0.38885199453206343, -0.33425057028323962, -0.30154253721772117, -0.27624597383474103, -0.26037811196890087, -0.24415404585285019, -0.22941383121091052, -0.21613496492144549, -0.20418199946908885, -0.19506212965221825, -0.18851924999737435, -0.18271143344718843, -0.16910750701722474, -0.17963775224176018, -0.19442870510406052, -0.20164216917300118, -0.20467219452758528, -0.20608246104917802, -0.20640259931640445, -0.2054749739152594, -0.20380549124050265, -0.20227296931678532, -0.20095834856297176, -0.20000430919304379, -0.19946673053844086, -0.1990733499211611, -0.19882136174363013, -0.19877442300323914, -0.19905182154377868, -0.19943266521643804, -0.19988524183849191, -0.20018905307631765, -0.20031895675727809, -0.20033991149804931, -0.20031574232920274, -0.20027004750680638, -0.20020472427796293, -0.20013382447039607, -0.2000635018536408, -0.20001515436367642, -0.19998427691514989, -0.19997263083178107, -0.19998545383896535, -0.20000134502238734, -0.2000127243362736, -0.20001564474711939, -0.20001267360809977, -0.20002707579781318, -0.20004087961702843, -0.20004212947389177]
        
        G3 = [-0.45000000000000001, -0.38058172993544187, -0.33756059941741273, -0.31015371357441396, -0.29214769368562965, -0.27545447937118606, -0.25871585649808154, -0.24254276680581988, -0.22758633129006092, -0.21417276895743134, -0.20237184768790789, -0.19369491041576814, -0.18721625333717057, -0.1794243868465818, -0.17052113574042196, -0.18534300640363346, -0.19601184621026671, -0.20185028431829469, -0.20476187496918136, -0.20602933256960082, -0.20598569228739247, -0.20492643756666848, -0.20339695828762758, -0.20196440373022595, -0.20070304052919338, -0.19986227854052355, -0.19933020476408528, -0.19898034831018035, -0.19878317651286193, -0.19886879323961787, -0.19915860801206181, -0.19953675278099042, -0.19992828019602107, -0.20015957043092364, -0.20025268671087426, -0.20028559516444974, -0.20027084877341045, -0.20022991487243985, -0.20016234295579871, -0.20009131445092507, -0.20003149397006523, -0.19998473356473795, -0.19996011913447218, -0.19995647168667186, -0.19996526209120422, -0.19996600297827097, -0.19997268800221216, -0.19998682275066659, -0.20000372259781876, -0.20001628681983963, -0.2000173314086407]
        
        assert num.allclose(gauge_values[0], G0)
        assert num.allclose(gauge_values[1], G1)
        assert num.allclose(gauge_values[2], G2)
        assert num.allclose(gauge_values[3], G3)

    #####################################################


    def test_initial_condition(self):
        """test_initial_condition

        Test that initial condition is output at time == 0 and that
        computed values change as system evolves
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

        domain = anuga.Domain(points, vertices)

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

        domain.set_boundary({'exterior': anuga.Reflective_boundary(domain)})

        domain.optimise_dry_cells = True

        #Evolution
        for t in domain.evolve(yieldstep=0.5, finaltime=2.0):
            stage = domain.quantities['stage'].vertex_values

            if t == 0.0:
                assert num.allclose(stage, initial_stage)
            else:
                assert not num.allclose(stage, initial_stage)

        os.remove(domain.get_name() + '.sww')

    #####################################################

    def test_second_order_flat_bed_onestep(self):

        #Create basic mesh
        points, vertices, boundary = anuga.rectangular(6, 6)

        #Create shallow water domain
        domain = anuga.Domain(points, vertices, boundary)
        domain.set_default_order(2)

        # Boundary conditions
        Br = anuga.Reflective_boundary(domain)
        Bd = anuga.Dirichlet_boundary([0.1, 0., 0.])
        domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})

        domain.check_integrity()

        # Evolution
        for t in domain.evolve(yieldstep=0.05, finaltime=0.05):
            pass


        # Data from earlier version of abstract_2d_finite_volumes
        assert num.allclose(domain.recorded_min_timestep, 0.0396825396825) or \
               num.allclose(domain.recorded_min_timestep, 0.0235282801879)
               
        assert num.allclose(domain.recorded_max_timestep, 0.0396825396825) or \
               num.allclose(domain.recorded_max_timestep, 0.0235282801879)


                
        assert num.allclose(domain.quantities['stage'].centroid_values[:12],
                            [0.00171396, 0.02656103, 0.00241523, 0.02656103,
                             0.00241523, 0.02656103, 0.00241523, 0.02656103,
                             0.00241523, 0.02656103, 0.00241523, 0.0272623], atol=1.0e-3) or \
               num.allclose(domain.quantities['stage'].centroid_values[:12],
                            [ 0.00053119,  0.02900893,  0.00077912,  0.02900893,
                              0.00077912,  0.02900893,  0.00077912,  0.02900893,
                              0.00077912,  0.02900893,  0.00077912,  0.02873746], atol=1.0e-3)

        domain.distribute_to_vertices_and_edges()

 

        assert num.allclose(domain.quantities['stage'].vertex_values[:12,0],
                            [ -1.96794125e-03,   2.65610347e-02,   0.00000000e+00,   2.65610347e-02,
                              -8.67361738e-19,   2.65610347e-02,   4.33680869e-19,   2.65610347e-02,
                              -2.16840434e-18,   2.65610347e-02,  -9.44042339e-05,   2.72623006e-02],
                            atol =1.0e-3) or \
                num.allclose(domain.quantities['stage'].vertex_values[:12,0],
                            [ -5.51381419e-04,   5.74866732e-02,   1.00006808e-15,   5.72387383e-02,
                              9.99851243e-16,   5.72387383e-02,   1.00050176e-15,   5.72387383e-02,
                              9.99417563e-16,   5.72387383e-02,   1.09882029e-05,   5.66957956e-02],
                             atol=1.0e-3)


        assert num.allclose(domain.quantities['stage'].vertex_values[:12,1],
                            [  5.14188587e-03,   2.65610347e-02,   0.00000000e+00,   2.65610347e-02,
                               8.67361738e-19,   2.65610347e-02,  -4.33680869e-19,   2.65610347e-02,
                               1.30104261e-18,   2.65610347e-02,   9.44042339e-05,   2.72623006e-02],
                            atol =1.0e-3) or \
               num.allclose(domain.quantities['stage'].vertex_values[:12,1],
                           [  1.59356551e-03,   5.72387383e-02,   1.00006808e-15,   5.72387383e-02,
                              1.00006808e-15,   5.72387383e-02,   9.99634403e-16,   5.72387383e-02,
                              1.00050176e-15,   5.72387383e-02,  -1.09882029e-05,   1.47582915e-02],
                            atol =1.0e-3)
        
        assert num.allclose(domain.quantities['stage'].vertex_values[:12,2],
                            [ 0.00196794,  0.02656103,  0.00724568,  0.02656103,
                              0.00724568,  0.02656103,  0.00724568,  0.02656103,
                              0.00724568,  0.02656103,  0.00724568,  0.0272623 ], atol =1.0e-3) or \
               num.allclose(domain.quantities['stage'].vertex_values[:12,2],
                            [ 0.00055138, -0.02769862,  0.00233737, -0.02745068,
                              0.00233737, -0.02745068,  0.00233737, -0.02745068,
                              0.00233737, -0.02745068,  0.00233737,  0.01475829], atol =1.0e-3)


        assert num.allclose(domain.quantities['xmomentum'].centroid_values[:12],
                            [0.00113961, 0.01302432, 0.00148672,
                             0.01302432, 0.00148672, 0.01302432,
                             0.00148672, 0.01302432, 0.00148672 ,
                             0.01302432, 0.00148672, 0.01337143], atol=1.0e-3) or \
               num.allclose(domain.quantities['xmomentum'].centroid_values[:12],
                        [ 0.00019529,  0.01425863,  0.00025665,
                          0.01425863,  0.00025665,  0.01425863,
                          0.00025665,  0.01425863,  0.00025665,
                          0.01425863,  0.00025665,  0.014423  ], atol=1.0e-3)
        
        assert num.allclose(domain.quantities['ymomentum'].centroid_values[:12],
                            [-2.91240050e-004 , 1.22721531e-004,
                             -1.22721531e-004,  1.22721531e-004 ,
                             -1.22721531e-004,  1.22721531e-004,
                             -1.22721531e-004 , 1.22721531e-004,
                             -1.22721531e-004,  1.22721531e-004,
                             -1.22721531e-004,  -4.57969873e-005], atol=1.0e-5) or \
               num.allclose(domain.quantities['ymomentum'].centroid_values[:12],
                            [ -6.38239364e-05,   2.16943067e-05,
                              -2.16943067e-05,   2.16943067e-05,
                              -2.16943067e-05,   2.16943067e-05,
                              -2.16943067e-05,   2.16943067e-05,
                              -2.16943067e-05,   2.16943067e-05,
                              -2.16943067e-05,  -4.62796434e-04], atol=1.0e-5)
        
        os.remove(domain.get_name() + '.sww')

    def test_second_order_flat_bed_moresteps(self):
        from mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        # Create shallow water domain
        domain = anuga.Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 2

        # Boundary conditions
        Br = anuga.Reflective_boundary(domain)
        Bd = anuga.Dirichlet_boundary([0.1, 0., 0.])
        domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})

        domain.check_integrity()

        # Evolution
        for t in domain.evolve(yieldstep=0.05, finaltime=0.1):
            pass

        # Data from earlier version of abstract_2d_finite_volumes
        #assert allclose(domain.recorded_min_timestep, 0.0396825396825)
        #assert allclose(domain.recorded_max_timestep, 0.0396825396825)
        #print domain.quantities['stage'].centroid_values

        os.remove(domain.get_name() + '.sww')

    def test_flatbed_first_order(self):
        from mesh_factory import rectangular

        # Create basic mesh
        N = 8
        points, vertices, boundary = rectangular(N, N)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 1
        domain.H0 = 1.0e-3 # As suggested in the manual

        # Boundary conditions
        Br = anuga.Reflective_boundary(domain)
        Bd = anuga.Dirichlet_boundary([0.2, 0., 0.])

        domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})
        domain.check_integrity()

        # Evolution
        for t in domain.evolve(yieldstep=0.02, finaltime=0.5):
            pass

        # FIXME: These numbers were from version before 25/10
        #assert allclose(domain.recorded_min_timestep, 0.0140413643926)
        #assert allclose(domain.recorded_max_timestep, 0.0140947355753)

        for i in range(3):
            #assert allclose(domain.quantities['stage'].edge_values[:4,i],
            #                [0.10730244,0.12337617,0.11200126,0.12605666])
            assert num.allclose(domain.quantities['xmomentum'].\
                                    edge_values[:4,i],
                                [0.07610894,0.06901572,0.07284461,0.06819712])
            #assert allclose(domain.quantities['ymomentum'].edge_values[:4,i],
            #                [-0.0060238, -0.00157404, -0.00309633, -0.0001637])

        os.remove(domain.get_name() + '.sww')

    def test_flatbed_second_order(self):
        from mesh_factory import rectangular

        # Create basic mesh
        N = 8
        points, vertices, boundary = rectangular(N, N)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.set_store_vertices_uniquely(True)
        domain.set_default_order(2)

        # Boundary conditions
        Br = anuga.Reflective_boundary(domain)
        Bd = anuga.Dirichlet_boundary([0.2, 0., 0.])

        domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})
        domain.check_integrity()

        # Evolution
        for t in domain.evolve(yieldstep=0.01, finaltime=0.03):
            pass


        
        msg = 'min step was %f instead of %f' % (domain.recorded_min_timestep,
                                                 0.0155604907816)

        assert num.allclose(domain.recorded_min_timestep, 0.0155604907816), msg
        assert num.allclose(domain.recorded_max_timestep, 0.0155604907816)


        assert num.allclose(domain.quantities['stage'].vertex_values[:4,0],
                            [-0.009, 0.0535,  0.0, 0.0535], atol=1.0e-3) or \
               num.allclose(domain.quantities['stage'].vertex_values[:4,0],
                      [-3.54158995e-03,1.22050959e-01,-2.36227400e-05,1.21501627e-01], atol=1.0e-3)
        
        
        assert num.allclose(domain.quantities['xmomentum'].vertex_values[:4,0],
                            [-0.008, 0.0368, 0.0, 0.0368], atol=1.0e-3) or \
               num.allclose(domain.quantities['xmomentum'].vertex_values[:4,0],
                      [-2.32056226e-03,9.10618822e-02, -1.06135035e-05,9.75175956e-02], atol=1.0e-3)

        assert num.allclose(domain.quantities['ymomentum'].vertex_values[:4,0],
                            [ 0.002 , 6.0e-04, 0.0, 6.0e-04],
                            atol=1.0e-3) or \
               num.allclose(domain.quantities['ymomentum'].vertex_values[:4,0],
                            [ 1.43500775e-03,  6.07102924e-05,   1.59329371e-06,   8.44572599e-03],
                            atol=1.0e-3)

        os.remove(domain.get_name() + '.sww')


    def test_flatbed_second_order_vmax_0(self):
        from mesh_factory import rectangular

        # Create basic mesh
        N = 8
        points, vertices, boundary = rectangular(N, N)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)

        domain.set_store_vertices_uniquely(True)
        domain.set_default_order(2)


        # Boundary conditions
        Br = anuga.Reflective_boundary(domain)
        Bd = anuga.Dirichlet_boundary([0.2, 0., 0.])

        domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})
        domain.check_integrity()

        # Evolution
        for t in domain.evolve(yieldstep=0.01, finaltime=0.03):
            pass


        assert num.allclose(domain.recorded_min_timestep, 0.0210448446782) or \
               num.allclose(domain.recorded_min_timestep, 0.0155604907816)
               
        assert num.allclose(domain.recorded_max_timestep, 0.0210448446782) or \
               num.allclose(domain.recorded_min_timestep, 0.0155604907816)


        assert num.allclose(domain.quantities['xmomentum'].vertex_values[:4,0],
                            [ -2.32056226e-03,   9.10618822e-02,  -1.06135035e-05,   9.75175956e-02],
                            atol=1.0e-3)

        assert num.allclose(domain.quantities['ymomentum'].vertex_values[:4,0],
                            [  1.43500775e-03,   6.07102924e-05,   1.59329371e-06,   8.44572599e-03],
                            atol=1.0e-3)

        os.remove(domain.get_name() + '.sww')

    def test_flatbed_second_order_distribute(self):
        #Use real data from anuga.abstract_2d_finite_volumes 2
        #painfully setup and extracted.

        from mesh_factory import rectangular

        # Create basic mesh
        N = 8
        points, vertices, boundary = rectangular(N, N)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)

        domain.set_store_vertices_uniquely(True)
        domain.set_default_order(2)

        # Boundary conditions
        Br = anuga.Reflective_boundary(domain)
        Bd = anuga.Dirichlet_boundary([0.2, 0., 0.])

        domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})
        domain.check_integrity()

        for V in [False, True]:
            if V:
                # Set centroids as if system had been evolved
                L = num.zeros(2*N*N, num.float)
                L[:32] = [7.21205592e-003, 5.35214298e-002, 1.00910824e-002,
                          5.35439433e-002, 1.00910824e-002, 5.35439433e-002,
                          1.00910824e-002, 5.35439433e-002, 1.00910824e-002,
                          5.35439433e-002, 1.00910824e-002, 5.35439433e-002,
                          1.00910824e-002, 5.35393928e-002, 1.02344264e-002,
                          5.59605058e-002, 0.00000000e+000, 3.31027800e-004,
                          0.00000000e+000, 4.37962142e-005, 0.00000000e+000,
                          4.37962142e-005, 0.00000000e+000, 4.37962142e-005,
                          0.00000000e+000, 4.37962142e-005, 0.00000000e+000,
                          4.37962142e-005, 0.00000000e+000, 4.37962142e-005,
                          0.00000000e+000, 5.57305948e-005]

                X = num.zeros(2*N*N, num.float)
                X[:32] = [6.48351607e-003, 3.68571894e-002, 8.50733285e-003,
                          3.68731327e-002, 8.50733285e-003, 3.68731327e-002,
                          8.50733285e-003, 3.68731327e-002, 8.50733285e-003,
                          3.68731327e-002, 8.50733285e-003, 3.68731327e-002,
                          8.50733285e-003, 3.68693861e-002, 8.65220973e-003,
                          3.85055387e-002, 0.00000000e+000, 2.86060840e-004,
                          0.00000000e+000, 3.58905503e-005, 0.00000000e+000,
                          3.58905503e-005, 0.00000000e+000, 3.58905503e-005,
                          0.00000000e+000, 3.58905503e-005, 0.00000000e+000,
                          3.58905503e-005, 0.00000000e+000, 3.58905503e-005,
                          0.00000000e+000, 4.57662812e-005]

                Y = num.zeros(2*N*N, num.float)
                Y[:32] = [-1.39463104e-003, 6.15600298e-004, -6.03637382e-004,
                          6.18272251e-004, -6.03637382e-004, 6.18272251e-004,
                          -6.03637382e-004, 6.18272251e-004, -6.03637382e-004,
                          6.18272251e-004, -6.03637382e-004, 6.18272251e-004,
                          -6.03637382e-004, 6.18599320e-004, -6.74622797e-004,
                          -1.48934756e-004, 0.00000000e+000, -5.35079969e-005,
                          0.00000000e+000, -2.57264987e-005, 0.00000000e+000,
                          -2.57264987e-005, 0.00000000e+000, -2.57264987e-005,
                          0.00000000e+000, -2.57264987e-005, 0.00000000e+000,
                          -2.57264987e-005, 0.00000000e+000, -2.57264987e-005,
                          0.00000000e+000, -2.57635178e-005]

                domain.set_quantity('stage', L, location='centroids')
                domain.set_quantity('xmomentum', X, location='centroids')
                domain.set_quantity('ymomentum', Y, location='centroids')

                domain.check_integrity()
            else:
                # Evolution
                for t in domain.evolve(yieldstep=0.01, finaltime=0.03):
                    pass

                
                assert num.allclose(domain.recorded_min_timestep, 0.0155604907816)
                assert num.allclose(domain.recorded_max_timestep, 0.0155604907816)

            #print domain.quantities['stage'].centroid_values[:4]
            #print domain.quantities['xmomentum'].centroid_values[:4]
            #print domain.quantities['ymomentum'].centroid_values[:4]                        
                
            #Centroids were correct but not vertices.
            #Hence the check of distribute below.

            if not V:

                assert num.allclose(domain.quantities['stage'].centroid_values[:4],
                               [0.00725574, 0.05350737, 0.01008413, 0.0535293], atol=1.0e-3) or \
                       num.allclose(domain.quantities['stage'].centroid_values[:4],
                               [0.00318259,  0.06261678,  0.00420215,  0.06285189], atol=1.0e-3)
                
                assert num.allclose(domain.quantities['xmomentum'].centroid_values[:4],
                               [0.00654964, 0.03684904, 0.00852561, 0.03686323],atol=1.0e-3) or \
                       num.allclose(domain.quantities['xmomentum'].centroid_values[:4],
                               [0.00218173, 0.04482164, 0.0026334,  0.04491656],atol=1.0e-3)       

                assert num.allclose(domain.quantities['ymomentum'].centroid_values[:4],
                               [-0.00143169, 0.00061027, -0.00062325, 0.00061408],atol=1.0e-3) or \
                       num.allclose(domain.quantities['ymomentum'].centroid_values[:4],
                [-6.46340592e-04,-6.16702557e-05,-2.83424134e-04, 6.48556590e-05],atol=1.0e-3)

                                    
                assert num.allclose(domain.quantities['xmomentum'].centroid_values[17], 0.0,
                                    atol=3.0e-4)                
            else:
                assert num.allclose(domain.quantities['xmomentum'].\
                                        centroid_values[17],
                                    0.00028606084)
                return #FIXME - Bailout for V True

            import copy

            XX = copy.copy(domain.quantities['xmomentum'].centroid_values)
            assert num.allclose(domain.quantities['xmomentum'].centroid_values,
                                XX)

            domain.distribute_to_vertices_and_edges()

            assert num.allclose(domain.quantities['xmomentum'].centroid_values[17], 0.0, atol=3.0e-4)

            assert num.allclose(domain.quantities['ymomentum'].vertex_values[:4,0],
                                [ 1.84104149e-03, 6.05658846e-04, 1.77092716e-07, 6.10687334e-04],
                                atol=1.0e-4) or \
                   num.allclose(domain.quantities['ymomentum'].vertex_values[:4,0],
                                [1.43500775e-03, 6.07102924e-05, 1.59329371e-06, 8.44572599e-03],
                                atol=1.0e-4)             

            assert num.allclose(domain.quantities['xmomentum'].vertex_values[:4,0],
                                [ -8.31184293e-03, 3.68841505e-02, -2.42843889e-06, 3.68900189e-02],
                                atol=1.0e-4) or \
                   num.allclose(domain.quantities['xmomentum'].vertex_values[:4,0],
                                [-2.32056226e-03, 9.10618822e-02, -1.06135035e-05, 9.75175956e-02],
                                rtol=1.0e-2)             


        os.remove(domain.get_name() + '.sww')


    def test_bedslope_problem_second_order_more_steps(self):
        """test_bedslope_problem_second_order_more_step

        Test shallow water balanced finite volumes
        """

        from mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)

        domain.set_store_vertices_uniquely(True)
        domain.set_default_order(2)

        # Bed-slope and friction at vertices (and interpolated elsewhere)
        def x_slope(x, y):
            return -x/3

        domain.set_quantity('elevation', x_slope)

        # Boundary conditions
        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        # Initial condition
        domain.set_quantity('stage', expression='elevation+0.05')
        domain.check_integrity()

        assert num.allclose(domain.quantities['stage'].centroid_values,
                            [ 0.01296296,  0.03148148,  0.01296296,
                              0.03148148,  0.01296296,  0.03148148,
                              0.01296296,  0.03148148,  0.01296296,
                              0.03148148,  0.01296296,  0.03148148,
                             -0.04259259, -0.02407407, -0.04259259,
                             -0.02407407, -0.04259259, -0.02407407,
                             -0.04259259, -0.02407407, -0.04259259,
                             -0.02407407, -0.04259259, -0.02407407,
                             -0.09814815, -0.07962963, -0.09814815,
                             -0.07962963, -0.09814815, -0.07962963,
                             -0.09814815, -0.07962963, -0.09814815,
                             -0.07962963, -0.09814815, -0.07962963,
                             -0.1537037 , -0.13518519, -0.1537037,
                             -0.13518519, -0.1537037,  -0.13518519,
                             -0.1537037 , -0.13518519, -0.1537037,
                             -0.13518519, -0.1537037,  -0.13518519,
                             -0.20925926, -0.19074074, -0.20925926,
                             -0.19074074, -0.20925926, -0.19074074,
                             -0.20925926, -0.19074074, -0.20925926,
                             -0.19074074, -0.20925926, -0.19074074,
                             -0.26481481, -0.2462963,  -0.26481481,
                             -0.2462963,  -0.26481481, -0.2462963,
                             -0.26481481, -0.2462963,  -0.26481481,
                             -0.2462963,  -0.26481481, -0.2462963])

        # Evolution
        for t in domain.evolve(yieldstep = 0.05, finaltime = 0.5):
            pass


        assert num.allclose(domain.quantities['stage'].centroid_values,
     [-0.02901283, -0.01619385, -0.03040423, -0.01564474, -0.02936756, -0.01507953,
      -0.02858108, -0.01491531, -0.02793549, -0.0147037,  -0.02792804, -0.014363,
      -0.07794301, -0.05951952, -0.07675098, -0.06182336, -0.07736607, -0.06079504,
      -0.07696764, -0.06039043, -0.07708793, -0.0601453,  -0.07669911, -0.06020719,
      -0.12223185, -0.10857309, -0.12286676, -0.10837591, -0.12386938, -0.10842744,
      -0.12363769, -0.10790002, -0.12304837, -0.10871278, -0.12543768, -0.10961026,
      -0.15664473, -0.14630207, -0.15838364, -0.14910271, -0.15804002, -0.15029627,
      -0.15829717, -0.1503869,  -0.15852604, -0.14971109, -0.15856346, -0.15205092,
      -0.20900931, -0.19658843, -0.20669607, -0.19558708, -0.20654186, -0.19492423,
      -0.20438765, -0.19492931, -0.20644142, -0.19423147, -0.20237449, -0.19198454,
      -0.13699658, -0.14209126, -0.13600697, -0.14334968, -0.1347657,  -0.14224247,
      -0.13442376, -0.14136926, -0.13501004, -0.14339389, -0.13479263, -0.14304073], atol=1.0e-2) or \
              num.allclose(domain.quantities['stage'].centroid_values,      
     [-0.03393968, -0.0166423,  -0.03253538, -0.01722023, -0.03270405, -0.01728606,
      -0.03277786, -0.0173903,  -0.03333736, -0.01743236, -0.03189526, -0.01463918,
      -0.07951756, -0.06410763, -0.07847973, -0.06350794, -0.07842429, -0.06240852,
      -0.07808697, -0.06255924, -0.07854662, -0.06322442, -0.07867314, -0.06287121,
      -0.11533356, -0.10559238, -0.11971301, -0.10742123, -0.1215759 , -0.10830046,
      -0.12202867, -0.10831703, -0.122214,   -0.10854099, -0.12343779, -0.11035803,
      -0.15725714, -0.14300757, -0.15559898, -0.1447275 , -0.15478568, -0.14483551,
      -0.15461918, -0.14489704, -0.15462074, -0.14516256, -0.15522298, -0.1452902,
      -0.22637615, -0.19192974, -0.20922654, -0.1907441 , -0.20900039, -0.19074809,
      -0.20897969, -0.19073365, -0.209195,   -0.19071396, -0.20922513, -0.19067714,
      -0.11357515, -0.14185801, -0.13224763, -0.14395805, -0.13379438, -0.14497114,
      -0.13437773, -0.14536013, -0.13607796, -0.14799629, -0.13148351, -0.15568502], atol=1.0e-1)



        assert num.allclose(domain.quantities['xmomentum'].centroid_values,
               [ 0.00478273,  0.003297,    0.00471129,  0.00320957,  0.00462171,  0.00320135,
                 0.00458295,  0.00317193,  0.00451704,  0.00314308,  0.00442684,  0.00320466,
                 0.01512907,  0.01150756,  0.01604672,  0.01156605,  0.01583911,  0.01135809,
                 0.01578499,  0.01132479,  0.01543668,  0.01100614,  0.01570445,  0.0120152,
                 0.04019477,  0.02721469,  0.03509982,  0.02735229,  0.03369315,  0.02727871,
                 0.03317931,  0.02706421,  0.03332704,  0.02722779,  0.03170258,  0.02556134,
                 0.07157025,  0.06074271,  0.07249738,  0.05570979,  0.07311261,  0.05428175,
                 0.07316986,  0.05379702,  0.0719581,   0.05230996,  0.07034837,  0.05468702,
                 0.08145001,  0.07753479,  0.08148804,  0.08119069,  0.08247295,  0.08134969,
                 0.0823216,   0.081411,    0.08190964,  0.08151441,  0.08163076,  0.08166174,
                 0.03680205,  0.0768216,   0.03943625,  0.07791183,  0.03930529,  0.07760588,
                 0.03949756,  0.07839929,  0.03992892,  0.08001416,  0.04444335,  0.08628738],
                            atol=1.0e-2) or \
               num.allclose(domain.quantities['xmomentum'].centroid_values,
               [ 0.00178414,  0.00147791,  0.00373636,  0.00169124,  0.00395649,  0.0014468,
                 0.00387617,  0.00135572,  0.00338418,  0.00134554,  0.00404961,  0.00252769,
                 0.01365204,  0.00890416,  0.01381613,  0.00986246,  0.01419385,  0.01145017,
                 0.01465116,  0.01125933,  0.01407359,  0.01055426,  0.01403563,  0.01095544,
                 0.04653827,  0.03018236,  0.03709973,  0.0265533 ,  0.0337694 ,  0.02541724,
                 0.03304266,  0.02535335,  0.03264548,  0.02484769,  0.03047682,  0.02205757,
                 0.07400338,  0.06470583,  0.07756503,  0.06098108,  0.07942593,  0.06086531,
                 0.07977427,  0.06074404,  0.07979513,  0.06019911,  0.07806395,  0.06011152,
                 0.07305045,  0.07883894,  0.08120393,  0.08166623,  0.08180501,  0.08166251,
                 0.0818353 ,  0.08169641,  0.08173762,  0.08174118,  0.08176467,  0.08181817,
                 0.01549926,  0.08259719,  0.01835423,  0.07302656,  0.01672924,  0.07198839,
                 0.01676006,  0.07223233,  0.01775672,  0.07362164,  0.01955846,  0.09361223],
                            atol=1.0e-2)


        assert num.allclose(domain.quantities['ymomentum'].centroid_values,
                            [ -1.09771684e-05,  -2.60328801e-05,  -1.03481959e-05,  -7.75907380e-05,
                              -5.00409090e-05,  -7.83807512e-05,  -3.60509918e-05,  -6.19321031e-05,
                              -1.40041903e-05,  -2.95707259e-05,   3.90296618e-06,   1.87556544e-05,
                              9.27848053e-05,   6.66937557e-07,   1.00653468e-04,   8.24734209e-06,
                              -1.04548672e-05,  -4.40402988e-05,  -2.95549946e-05,  -1.86360736e-05,
                              1.12527016e-04,   1.27240681e-04,   2.02147284e-04,   9.18457482e-05,
                              1.41781748e-03,   7.23407624e-04,   5.09160779e-04,   1.29136939e-04,
                              -4.70131286e-05,  -1.00180290e-04,  -1.76806614e-05,  -4.19421384e-06,
                              -6.17759681e-05,  -3.02124967e-05,   4.32689360e-04,   5.49717934e-04,
                              1.15031101e-03,   1.02737170e-03,   5.77937840e-04,   3.36230967e-04,
                              5.44877516e-04,  -7.28594977e-05,   4.60064858e-04,  -3.94125434e-05,
                              7.48242964e-04,   2.88528341e-04,   6.25148041e-05,  -1.74477175e-04,
                              -5.06603166e-05,   7.07720999e-04,  -2.04937748e-04,   3.38595573e-05,
                              -4.64116229e-05,   1.49325340e-04,  -2.41342281e-05,   1.83817970e-04,
                              -1.44417277e-05,   2.47823834e-04,   7.91185571e-05,   1.71615793e-04,
                              1.56883043e-03,   8.39352974e-04,   3.23353846e-03,   1.70597880e-03,
                              2.27789107e-03,   1.48928169e-03,   2.09854126e-03,   1.50248643e-03,
                              2.83029467e-03,   1.09151499e-03,   6.52455118e-03,  -2.04468968e-03],
                            atol=1.0e-3) or \
               num.allclose(domain.quantities['ymomentum'].centroid_values,
                             [ -1.24810991e-04,  -3.08228767e-04,  -1.56701128e-04,  -1.01904208e-04,
                               -3.36282053e-05,  -1.17956840e-04,  -3.55986664e-05,  -9.38578996e-05,
                               7.13704069e-05,   2.47022380e-05,   1.71121489e-04,   2.65941677e-04,
                               6.90055205e-04,   1.99195585e-04,   1.33804448e-04,  -1.66563316e-04,
                               -2.00962830e-04,  -3.81664130e-05,  -9.50456053e-05,  -3.14620186e-06,
                               1.29388102e-04,   3.16945980e-04,   4.77556581e-04,   2.57217342e-04,
                               1.42300612e-03,   9.60776359e-04,   5.08941026e-04,   1.06939990e-04,
                               6.37673950e-05,  -2.69783047e-04,  -8.55760509e-05,  -2.12987309e-04,
                               -5.86840949e-06,  -9.75751293e-05,   8.25447727e-04,   1.14139065e-03,
                               8.56206468e-04,   3.83113329e-04,   1.75041847e-04,   4.39999200e-04,
                               3.75156469e-04,   2.48774698e-04,   4.09671654e-04,   2.07125615e-04,
                               4.59587647e-04,   2.70581830e-04,  -1.24082302e-06,  -4.29155678e-04,
                               -9.66841218e-03,   4.93278794e-04,  -5.25778806e-06,  -4.90396857e-05,
                               -9.75373988e-06,   7.28023591e-06,  -5.20499868e-06,   3.61013683e-05,
                               -7.54919544e-06,   4.14115771e-05,  -1.35778834e-05,  -2.23991903e-05,
                               3.63635844e-02,   5.29865244e-04,   5.13015379e-03,   1.19233296e-03,
                               4.70681275e-04,   2.62292296e-04,  -1.28084045e-04,   7.04826916e-04,
                               1.50377987e-04,   1.35053814e-03,   1.30710492e-02,   1.93011958e-03],
                            atol=1.0e-1)
                            


        os.remove(domain.get_name() + '.sww')


    def test_temp_play(self):
        from mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(5, 5)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)

        domain.set_store_vertices_uniquely(True)
        domain.set_default_order(2)



        # Bed-slope and friction at vertices (and interpolated elsewhere)
        def x_slope(x, y):
            return -x/3

        domain.set_quantity('elevation', x_slope)

        # Boundary conditions
        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        # Initial condition
        domain.set_quantity('stage', expression='elevation+0.05')
        domain.check_integrity()

        # Evolution
        for t in domain.evolve(yieldstep=0.05, finaltime=0.1):
            pass


        assert num.allclose(domain.quantities['stage'].centroid_values[:4],
                            [ 0.01,   0.015,  0.01,  0.015], atol=1.0e-2)
                            
        assert num.allclose(domain.quantities['xmomentum'].centroid_values[:4],
                            [ 0.015,  0.01,  0.015,  0.01], atol=1.0e-2)
        
        assert num.allclose(domain.quantities['ymomentum'].centroid_values[:4],
                            [  0.0,  0.0,  0.0,   0.0]
                            , atol=1.0e-3)

        os.remove(domain.get_name() + '.sww')

    def test_complex_bed(self):
        # No friction is tested here

        from mesh_factory import rectangular

        N = 12
        points, vertices, boundary = rectangular(N, N/2, len1=1.2, len2=0.6,
                                                 origin=(-0.07, 0))


        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.set_default_order(2)
        domain.set_timestepping_method('rk2')
        domain.set_beta(1.0)

        inflow_stage = 0.1
        Z = Weir(inflow_stage)
        domain.set_quantity('elevation', Z)

        Br = anuga.Reflective_boundary(domain)
        Bd = anuga.Dirichlet_boundary([inflow_stage, 0.0, 0.0])
        domain.set_boundary({'left': Bd, 'right': Br, 'bottom': Br, 'top': Br})

        domain.set_quantity('stage', expression='elevation')

        for t in domain.evolve(yieldstep=0.02, finaltime=0.2):
            pass

        #FIXME: These numbers were from version before 25/10
        #assert allclose(domain.quantities['stage'].centroid_values,
# [3.95822638e-002,  5.61022588e-002,  4.66437868e-002,  5.73081011e-002,
#  4.72394613e-002,  5.74684939e-002,  4.74309483e-002,  5.77458084e-002,
#  4.80628177e-002,  5.85656225e-002,  4.90498542e-002,  6.02609831e-002,
#  1.18470315e-002,  1.75136443e-002,  1.18035266e-002,  2.15565695e-002,
#  1.31620268e-002,  2.14351640e-002,  1.32351076e-002,  2.15450687e-002,
#  1.36414028e-002,  2.24274619e-002,  1.51689511e-002,  2.21789655e-002,
# -7.54337535e-003, -6.86362021e-004, -7.74146760e-003, -1.83756530e-003,
# -8.16773628e-003, -4.49916813e-004, -8.08202599e-003, -3.91118720e-004,
# -8.10292716e-003, -3.88584984e-004, -7.35226124e-003,  2.73985295e-004,
#  1.86166683e-001,  8.74070369e-002,  1.86166712e-001,  8.74035875e-002,
#  6.11666935e-002, -3.76173225e-002, -6.38333276e-002, -3.76147365e-002,
#  6.11666725e-002,  8.73846774e-002,  1.86166697e-001,  8.74171550e-002,
# -4.83333333e-002,  1.18333333e-001, -4.83333333e-002,  1.18333333e-001,
# -4.83333333e-002, -6.66666667e-003, -1.73333333e-001, -1.31666667e-001,
# -1.73333333e-001, -6.66666667e-003, -4.83333333e-002,  1.18333333e-001,
# -2.48333333e-001, -2.31666667e-001, -2.48333333e-001, -2.31666667e-001,
# -2.48333333e-001, -2.31666667e-001, -2.48333333e-001, -2.31666667e-001,
# -2.48333333e-001, -2.31666667e-001, -2.48333333e-001, -2.31666667e-001,
# -4.65000000e-001, -3.65000000e-001, -4.65000000e-001, -3.65000000e-001,
# -4.65000000e-001, -3.65000000e-001, -4.65000000e-001, -3.65000000e-001,
# -4.65000000e-001, -3.65000000e-001, -4.65000000e-001, -3.65000000e-001,
# -5.98333333e-001, -5.81666667e-001, -5.98333333e-001, -5.81666667e-001,
# -5.98333333e-001, -5.81666667e-001, -5.98333333e-001, -5.81666667e-001,
# -5.98333333e-001, -5.81666667e-001, -5.98333333e-001, -5.81666667e-001,
# -6.48333333e-001, -6.31666667e-001, -6.48333333e-001, -6.31666667e-001,
# -6.48333333e-001, -6.31666667e-001, -6.48333333e-001, -6.31666667e-001,
# -6.48333333e-001, -6.31666667e-001, -6.48333333e-001, -6.31666667e-001,
# -5.31666667e-001, -5.98333333e-001, -5.31666667e-001, -5.98333333e-001,
# -5.31666667e-001, -5.98333333e-001, -5.31666667e-001, -5.98333333e-001,
# -5.31666667e-001, -5.98333333e-001, -5.31666667e-001, -5.98333333e-001,
# -4.98333333e-001, -4.81666667e-001, -4.98333333e-001, -4.81666667e-001,
# -4.98333333e-001, -4.81666667e-001, -4.98333333e-001, -4.81666667e-001,
# -4.98333333e-001, -4.81666667e-001, -4.98333333e-001, -4.81666667e-001,
# -5.48333333e-001, -5.31666667e-001, -5.48333333e-001, -5.31666667e-001,
# -5.48333333e-001, -5.31666667e-001, -5.48333333e-001, -5.31666667e-001,
# -5.48333333e-001, -5.31666667e-001, -5.48333333e-001, -5.31666667e-001])

        os.remove(domain.get_name() + '.sww')


    def test_tight_slope_limiters(self):
        """Test that new slope limiters (Feb 2007) don't induce extremely
        small timesteps. This test actually reveals the problem as it
        was in March-April 2007
        """
        import time, os
        from Scientific.IO.NetCDF import NetCDFFile
        from mesh_factory import rectangular_cross

        # Create basic mesh
        points, vertices, boundary = rectangular_cross(2, 2)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.set_default_order(2)
        domain.set_beta(1.0)
        domain.set_timestepping_method('euler')
        #domain.set_CFL(0.5)
        

        # This will pass
        #domain.tight_slope_limiters = 1
        #domain.H0 = 0.01

        # This will fail
        #domain.tight_slope_limiters = 1
        #domain.H0 = 0.001

        # This will pass provided C extension implements limiting of
        # momentum in _compute_speeds
        #domain.tight_slope_limiters = 1
        #domain.H0 = 0.001
        #domain.protect_against_isolated_degenerate_timesteps = True

        # Set some field values
        domain.set_quantity('elevation', lambda x,y: -x)
        domain.set_quantity('friction', 0.03)

        # Boundary conditions
        B = Transmissive_boundary(domain)
        domain.set_boundary({'left': B, 'right': B, 'top': B, 'bottom': B})

        # Initial condition - with jumps
        bed = domain.quantities['elevation'].vertex_values
        stage = num.zeros(bed.shape, num.float)

        h = 0.3
        for i in range(stage.shape[0]):
            if i % 2 == 1:
                stage[i,:] = bed[i,:] + h
            else:
                stage[i,:] = bed[i,:]

        domain.set_quantity('stage', stage)

        domain.distribute_to_vertices_and_edges()

        domain.set_name('tight_limiters')
        domain.smooth = True
        domain.reduction = mean
        domain.set_datadir('.')
        domain.smooth = False
        domain.store = True

        # Evolution
        for t in domain.evolve(yieldstep=0.1, finaltime=0.3):
            #domain.write_time(track_speeds=True)
            stage = domain.quantities['stage'].vertex_values

            # Get NetCDF
            #fid = NetCDFFile(domain.writer.filename, netcdf_mode_r)
            #stage_file = fid.variables['stage']

            #fid.close()

        os.remove(domain.writer.filename)



    def test_pmesh2Domain(self):
         import os
         import tempfile

         fileName = tempfile.mktemp(".tsh")
         file = open(fileName, "w")
         file.write("4 3 # <vertex #> <x> <y> [attributes]\n \
0 0.0 0.0 0.0 0.0 0.01 \n \
1 1.0 0.0 10.0 10.0 0.02  \n \
2 0.0 1.0 0.0 10.0 0.03  \n \
3 0.5 0.25 8.0 12.0 0.04  \n \
# Vert att title  \n \
elevation  \n \
stage  \n \
friction  \n \
2 # <triangle #> [<vertex #>] [<neigbouring triangle #>]  \n\
0 0 3 2 -1  -1  1 dsg\n\
1 0 1 3 -1  0 -1   ole nielsen\n\
4 # <segment #> <vertex #>  <vertex #> [boundary tag] \n\
0 1 0 2 \n\
1 0 2 3 \n\
2 2 3 \n\
3 3 1 1 \n\
3 0 # <x> <y> [attributes] ...Mesh Vertices... \n \
0 216.0 -86.0 \n \
1 160.0 -167.0 \n \
2 114.0 -91.0 \n \
3 # <vertex #>  <vertex #> [boundary tag] ...Mesh Segments... \n \
0 0 1 0 \n \
1 1 2 0 \n \
2 2 0 0 \n \
0 # <x> <y> ...Mesh Holes... \n \
0 # <x> <y> <attribute>...Mesh Regions... \n \
0 # <x> <y> <attribute>...Mesh Regions, area... \n\
#Geo reference \n \
56 \n \
140 \n \
120 \n")
         file.close()

         tags = {}
         b1 =  anuga.Dirichlet_boundary(conserved_quantities = num.array([0.0]))
         b2 =  anuga.Dirichlet_boundary(conserved_quantities = num.array([1.0]))
         b3 =  anuga.Dirichlet_boundary(conserved_quantities = num.array([2.0]))
         tags["1"] = b1
         tags["2"] = b2
         tags["3"] = b3

         domain = Domain(mesh_filename=fileName)
                         # verbose=True, use_cache=True)

         ## check the quantities
         answer = [[0., 8., 0.],
                   [0., 10., 8.]]
         assert num.allclose(domain.quantities['elevation'].vertex_values,
                             answer)

         answer = [[0., 12., 10.],
                   [0., 10., 12.]]
         assert num.allclose(domain.quantities['stage'].vertex_values,
                             answer)

         answer = [[0.01, 0.04, 0.03],
                   [0.01, 0.02, 0.04]]
         assert num.allclose(domain.quantities['friction'].vertex_values,
                             answer)

         tagged_elements = domain.get_tagged_elements()
         assert num.allclose(tagged_elements['dsg'][0], 0)
         assert num.allclose(tagged_elements['ole nielsen'][0], 1)

         msg = "test_tags_to_boundaries failed. Single boundary wasn't added."
         self.failUnless( domain.boundary[(1, 0)]  == '1', msg)
         self.failUnless( domain.boundary[(1, 2)]  == '2', msg)
         self.failUnless( domain.boundary[(0, 1)]  == '3', msg)
         self.failUnless( domain.boundary[(0, 0)]  == 'exterior', msg)
         msg = "test_pmesh2Domain Too many boundaries"
         self.failUnless( len(domain.boundary)  == 4, msg)

         # FIXME change to use get_xllcorner
         msg = 'Bad geo-reference'
         self.failUnless(domain.geo_reference.xllcorner  == 140.0, msg)

         domain = Domain(fileName)

         answer = [[0., 8., 0.],
                   [0., 10., 8.]]
         assert num.allclose(domain.quantities['elevation'].vertex_values,
                             answer)

         answer = [[0., 12., 10.],
                   [0., 10., 12.]]
         assert num.allclose(domain.quantities['stage'].vertex_values,
                             answer)

         answer = [[0.01, 0.04, 0.03],
                   [0.01, 0.02, 0.04]]
         assert num.allclose(domain.quantities['friction'].vertex_values,
                             answer)

         tagged_elements = domain.get_tagged_elements()
         assert num.allclose(tagged_elements['dsg'][0], 0)
         assert num.allclose(tagged_elements['ole nielsen'][0], 1)

         msg = "test_tags_to_boundaries failed. Single boundary wasn't added."
         self.failUnless(domain.boundary[(1, 0)]  == '1', msg)
         self.failUnless(domain.boundary[(1, 2)]  == '2', msg)
         self.failUnless(domain.boundary[(0, 1)]  == '3', msg)
         self.failUnless(domain.boundary[(0, 0)]  == 'exterior', msg)
         msg = "test_pmesh2Domain Too many boundaries"
         self.failUnless(len(domain.boundary)  == 4, msg)

         # FIXME change to use get_xllcorner
         msg = 'Bad geo_reference'
         self.failUnless(domain.geo_reference.xllcorner  == 140.0, msg)

         os.remove(fileName)

    def test_get_lone_vertices(self):
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #             bac,     bce,     ecf,     dbe
        vertices = [[1,0,2], [1,2,4], [4,2,5], [3,1,4] ]
        boundary = {(0, 0): 'Third',
                    (0, 2): 'First',
                    (2, 0): 'Second',
                    (2, 1): 'Second',
                    (3, 1): 'Second',
                    (3, 2): 'Third'}

        domain = Domain(points, vertices, boundary)
        domain.get_lone_vertices()

    def test_fitting_using_shallow_water_domain(self):
        #Mesh in zone 56 (absolute coords)

        x0 = 314036.58727982
        y0 = 6224951.2960092

        a = [x0+0.0, y0+0.0]
        b = [x0+0.0, y0+2.0]
        c = [x0+2.0, y0+0.0]
        d = [x0+0.0, y0+4.0]
        e = [x0+2.0, y0+2.0]
        f = [x0+4.0, y0+0.0]

        points = [a, b, c, d, e, f]

        #             bac,     bce,     ecf,     dbe
        elements = [[1,0,2], [1,2,4], [4,2,5], [3,1,4] ]

        # absolute going in ..
        mesh4 = Domain(points, elements, geo_reference=Geo_reference(56, 0, 0))
        mesh4.check_integrity()
        quantity = Quantity(mesh4)

        # Get (enough) datapoints (relative to georef)
        data_points_rel = [[ 0.66666667, 0.66666667],
                           [ 1.33333333, 1.33333333],
                           [ 2.66666667, 0.66666667],
                           [ 0.66666667, 2.66666667],
                           [ 0.0,        1.0],
                           [ 0.0,        3.0],
                           [ 1.0,        0.0],
                           [ 1.0,        1.0],
                           [ 1.0,        2.0],
                           [ 1.0,        3.0],
                           [ 2.0,        1.0],
                           [ 3.0,        0.0],
                           [ 3.0,        1.0]]

        data_geo_spatial = Geospatial_data(data_points_rel,
                                           geo_reference=Geo_reference(56,
                                                                       x0,
                                                                       y0))
        data_points_absolute = data_geo_spatial.get_data_points(absolute=True)
        attributes = linear_function(data_points_absolute)
        att = 'spam_and_eggs'

        # Create .txt file
        ptsfile = tempfile.mktemp(".txt")
        file = open(ptsfile, "w")
        file.write(" x,y," + att + " \n")
        for data_point, attribute in map(None, data_points_absolute, attributes):
            row = (str(data_point[0]) + ',' +
                   str(data_point[1]) + ',' +
                   str(attribute))
            file.write(row + "\n")
        file.close()

        # Check that values can be set from file
        quantity.set_values(filename=ptsfile, attribute_name=att, alpha=0)
        answer = linear_function(quantity.domain.get_vertex_coordinates())

        assert num.allclose(quantity.vertex_values.flat, answer)

        # Check that values can be set from file using default attribute
        quantity.set_values(filename = ptsfile, alpha = 0)
        assert num.allclose(quantity.vertex_values.flat, answer)

        # Cleanup
        import os
        os.remove(ptsfile)

    def test_fitting_example_that_crashed(self):
        """This unit test has been derived from a real world example
        (the Towradgi '98 rainstorm simulation).

        It shows a condition where fitting as called from set_quantity crashes
        when ANUGA mesh is reused. The test passes in the case where a new mesh
        is created.

        See ticket:314
        """

        verbose = False
        
        # Get path where this test is run
        path = get_pathname_from_package('anuga.shallow_water')        
        
        
        #----------------------------------------------------------------------
        # Create domain
        #--------------------------------------------------------------------
        W = 303400
        N = 6195800
        E = 308640
        S = 6193120
        bounding_polygon = [[W, S], [E, S], [E, N], [W, N]]

        offending_regions = []

        # From culvert 8
        offending_regions.append([[307611.43896231, 6193631.6894806],
                                  [307600.11394969, 6193608.2855474],
                                  [307597.41349586, 6193609.59227963],
                                  [307608.73850848, 6193632.99621282]])
        offending_regions.append([[307633.69143231, 6193620.9216536],
                                  [307622.36641969, 6193597.5177204],
                                  [307625.06687352, 6193596.21098818],
                                  [307636.39188614, 6193619.61492137]])

        # From culvert 9
        offending_regions.append([[306326.69660524, 6194818.62900522],
                                  [306324.67939476, 6194804.37099478],
                                  [306323.75856492, 6194804.50127295],
                                  [306325.7757754, 6194818.7592834]])
        offending_regions.append([[306365.57160524, 6194813.12900522],
                                  [306363.55439476, 6194798.87099478],
                                  [306364.4752246, 6194798.7407166],
                                  [306366.49243508, 6194812.99872705]])

        # From culvert 10
        offending_regions.append([[306955.071019428608, 6194465.704096679576],
                                  [306951.616980571358, 6194457.295903320424],
                                  [306950.044491164153, 6194457.941873183474],
                                  [306953.498530021403, 6194466.350066542625]])
        offending_regions.append([[307002.540019428649, 6194446.204096679576],
                                  [306999.085980571399, 6194437.795903320424],
                                  [307000.658469978604, 6194437.149933457375],
                                  [307004.112508835853, 6194445.558126816526]])

        interior_regions = []
        for polygon in offending_regions:
            interior_regions.append( [polygon, 100] ) 

        meshname = os.path.join(path, 'offending_mesh.msh')
        anuga.create_mesh_from_regions(bounding_polygon,
                                 boundary_tags={'south': [0], 'east': [1],
                                                'north': [2], 'west': [3]},
                                 maximum_triangle_area=1000000,
                                 interior_regions=interior_regions,
                                 filename=meshname,
                                 use_cache=False,
                                 verbose=verbose)

        domain = anuga.Domain(meshname, use_cache=False, verbose=verbose)

        #----------------------------------------------------------------------
        # Fit data point to mesh
        #----------------------------------------------------------------------

        points_file = os.path.join(path, 'offending_point.pts')

        # Offending point
        G = Geospatial_data(data_points=[[306953.344, 6194461.5]],
                            attributes=[1])
        G.export_points_file(points_file)
        
        try:
            domain.set_quantity('elevation', filename=points_file,
                                use_cache=False, verbose=verbose, alpha=0.01)
        except RuntimeError, e:
            msg = 'Test failed: %s' % str(e)
            raise Exception(msg)
            # clean up in case raise fails
            os.remove(meshname)
            os.remove(points_file)
        else:
            os.remove(meshname)
            os.remove(points_file)            
        

    def test_fitting_example_that_crashed_2(self):
        """test_fitting_example_that_crashed_2
        
        This unit test has been derived from a real world example 
        (the JJKelly study, by Petar Milevski).
        
        It shows a condition where set_quantity crashes due to AtA
        not being built properly
        
        See ticket:314
        """

        verbose = False        
        
        # Get path where this test is run
        path = get_pathname_from_package('anuga.shallow_water')        

        meshname = os.path.join(path, 'test_mesh.msh')
        W = 304180
        S = 6185270
        E = 307650
        N = 6189040
        maximum_triangle_area = 1000000

        bounding_polygon = [[W, S], [E, S], [E, N], [W, N]]

        anuga.create_mesh_from_regions(bounding_polygon,
                                 boundary_tags={'south': [0], 
                                                'east': [1], 
                                                'north': [2], 
                                                'west': [3]},
                                 maximum_triangle_area=maximum_triangle_area,
                                 filename=meshname,
                                 use_cache=False,
                                 verbose=verbose)

        domain = anuga.Domain(meshname, use_cache=True, verbose=verbose)
        
        # Large test set revealed one problem
        points_file = os.path.join(path, 'test_points_large.csv')
        try:
            domain.set_quantity('elevation', filename=points_file,
                                use_cache=False, verbose=verbose)
        except AssertionError, e:
            msg = 'Test failed: %s' % str(e)
            raise Exception(msg)
            # Cleanup in case this failed
            os.remove(meshname)

        # Small test set revealed another problem
        points_file = os.path.join(path, 'test_points_small.csv')
        try:    
            domain.set_quantity('elevation', filename=points_file,
                                use_cache=False, verbose=verbose)                            
        except AssertionError, e:
            msg = 'Test failed: %s' % str(e)
            raise Exception(msg)
            # Cleanup in case this failed
            os.remove(meshname)
        else:
            os.remove(meshname)




    def test_variable_elevation(self):            
        """test_variable_elevation

        This will test that elevagtion van be stored in sww files
        as a time dependent quantity.
        
        It will also chck that storage of other quantities 
        can be controlled this way.
        """
        
        #---------------------------------------------------------------------
        # Setup computational domain
        #---------------------------------------------------------------------
        length = 8.
        width = 6.
        dx = dy = 1    # Resolution: Length of subdivisions on both axes
        
        inc = 0.05 # Elevation increment

        points, vertices, boundary = rectangular_cross(int(length/dx), 
                                                       int(width/dy),
                                                       len1=length, 
                                                       len2=width)
        domain = anuga.Domain(points, vertices, boundary)
        domain.set_name('channel_variable_test')  # Output name
        domain.set_quantities_to_be_stored({'elevation': 2,
                                            'stage': 2})

        #---------------------------------------------------------------------
        # Setup initial conditions
        #---------------------------------------------------------------------

        def pole_increment(x,y):
            """This provides a small increment to a pole located mid stream
            For use with variable elevation data
            """
            
            z = 0.0*x

            N = len(x)
            for i in range(N):
                # Pole
                if (x[i] - 4)**2 + (y[i] - 2)**2 < 1.0**2:
                    z[i] += inc
            return z
            
        domain.set_quantity('elevation', 0.0)    # Flat bed initially
        domain.set_quantity('friction', 0.01)    # Constant friction
        domain.set_quantity('stage', 0.0)        # Dry initial condition

        #------------------------------------------------------------------
        # Setup boundary conditions
        #------------------------------------------------------------------
        Bi = anuga.Dirichlet_boundary([0.4, 0, 0])          # Inflow
        Br = anuga.Reflective_boundary(domain)              # Solid reflective wall
        Bo = anuga.Dirichlet_boundary([-5, 0, 0])           # Outflow

        domain.set_boundary({'left': Bi, 'right': Bo, 'top': Br, 'bottom': Br})

        #-------------------------------------------------------------------
        # Evolve system through time
        #-------------------------------------------------------------------

        for t in domain.evolve(yieldstep=1, finaltime=3.0):
            #print domain.timestepping_statistics()

            domain.add_quantity('elevation', pole_increment)
        
            
        # Check that quantities have been stored correctly    
        from Scientific.IO.NetCDF import NetCDFFile
        sww_file = domain.get_name() + '.sww'
        fid = NetCDFFile(sww_file)

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        stage = fid.variables['stage'][:]
        elevation = fid.variables['elevation'][:]
        fid.close()

        os.remove(sww_file)
        
                   
        assert len(stage.shape) == 2
        assert len(elevation.shape) == 2        
        
        M, N = stage.shape
                
        for i in range(M): 
            # For each timestep
            assert num.allclose(max(elevation[i,:]), i * inc) 



#################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_swb_basic, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
