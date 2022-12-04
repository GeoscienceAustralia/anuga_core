#!/usr/bin/env python



from builtins import zip
from builtins import str
from builtins import range
from builtins import object
import unittest, os, time
import os.path
from math import pi, sqrt
import tempfile

from anuga.file.netcdf import NetCDFFile
from anuga.file.sww import extent_sww

from anuga.config import g, epsilon
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a
from anuga.utilities.numerical_tools import mean, ensure_numeric
from anuga.geometry.polygon import is_inside_polygon
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.geospatial_data.geospatial_data import Geospatial_data
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross, \
                                            rectangular
from anuga.abstract_2d_finite_volumes.quantity import Quantity
from anuga.shallow_water.forcing import Inflow, Cross_section
from anuga.geospatial_data.geospatial_data import ensure_geospatial

from anuga.utilities.system_tools import get_pathname_from_package

from anuga.abstract_2d_finite_volumes.generic_boundary_conditions \
        import Dirichlet_boundary
from anuga.shallow_water.forcing import Rainfall, Wind_stress
from anuga.shallow_water.forcing import Inflow, Cross_section
from anuga.shallow_water.sww_interrogate import get_flow_through_cross_section

from anuga.shallow_water.shallow_water_domain import Domain

# boundary functions
from anuga.shallow_water.boundaries import Reflective_boundary, \
            Field_boundary, Transmissive_momentum_set_stage_boundary, \
            Transmissive_stage_zero_momentum_boundary
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions\
     import Transmissive_boundary, Dirichlet_boundary, \
            Time_boundary, File_boundary, AWI_boundary

import numpy as num
from anuga.config import g

from pprint import pprint

# Get gateway to C implementation of flux function for direct testing
from anuga.shallow_water.shallow_water_ext import flux_function_central as flux_function
from anuga.shallow_water.shallow_water_ext import rotate


def set_bottom_friction(tag, elements, domain):
    if tag == "bottom":
        domain.set_quantity('friction', 0.09, indices = elements)

def set_top_friction(tag, elements, domain):
    if tag == "top":
        domain.set_quantity('friction', 1., indices = elements)


def set_all_friction(tag, elements, domain):
    if tag == 'all':
        new_values = domain.get_quantity('friction').get_values(indices = elements) + 10.0

        domain.set_quantity('friction', new_values, indices = elements)


# For test_fitting_using_shallow_water_domain example
def linear_function(point):
    point = num.array(point)
    return point[:,0]+point[:,1]

# for help creating asc and dem files
def axes2points(x, y):
    """Generate all combinations of grid point coordinates from x and y axes

    Args:
        * x: x coordinates (array)
        * y: y coordinates (array)

    Returns:
        * P: Nx2 array consisting of coordinates for all
             grid points defined by x and y axes. The x coordinate
             will vary the fastest to match the way 2D numpy
             arrays are laid out by default ('C' order). That way,
             the x and y coordinates will match a corresponding
             2D array A when flattened (A.flat[:] or A.reshape(-1))

    Note:
        Example

        x = [1, 2, 3]
        y = [10, 20]

        P = [[1, 10],
             [2, 10],
             [3, 10],
             [1, 20],
             [2, 20],
             [3, 20]]
    """
    import numpy

    # Reverse y coordinates to have them start at bottom of array
    y = numpy.flipud(y)

    # Repeat x coordinates for each y (fastest varying)
    X = numpy.kron(numpy.ones(len(y)), x)

    # Repeat y coordinates for each x (slowest varying)
    Y = numpy.kron(y, numpy.ones(len(x)))

    # Check
    N = len(X)
    assert len(Y) == N

    # Create Nx2 array of x and y coordinates
    X = numpy.reshape(X, (N, 1))
    Y = numpy.reshape(Y, (N, 1))
    P = numpy.concatenate((X, Y), axis=1)

    # Return
    return P

class Weir(object):
    """Set a bathymetry for weir with a hole and a downstream gutter
    x,y are assumed to be in the unit square
    """

    def __init__(self, stage):
        self.inflow_stage = stage

    def __call__(self, x, y):
        N = len(x)
        assert N == len(y)

        z = num.zeros(N, float)
        for i in range(N):
            z[i] = -x[i] / 2  # General slope

            # Flattish bit to the left
            if x[i] < 0.3:
                z[i] = -x[i] / 10

            # Weir
            if x[i] >= 0.3 and x[i] < 0.4:
                z[i] = -x[i] + 0.9

            # Dip
            x0 = 0.6
            depth = -1.0
            plateaux = -0.6
            if y[i] < 0.7:
                if x[i] > x0 and x[i] < 0.9:
                    z[i] = depth
                # RHS plateaux
                if x[i] >= 0.9:
                    z[i] = plateaux
            elif y[i] >= 0.7 and y[i] < 1.5:
                # Restrict and deepen
                if x[i] >= x0 and x[i] < 0.8:
                    z[i] = depth - (y[i] / 3 - 0.3)
                elif x[i] >= 0.8:
                    # RHS plateaux
                    z[i] = plateaux
            elif y[i] >= 1.5:
                if x[i] >= x0 and x[i] < 0.8 + (y[i]-1.5)/1.2:
                    # Widen up and stay at constant depth
                    z[i] = depth-1.5/5
                elif x[i] >= 0.8 + (y[i]-1.5)/1.2:
                    # RHS plateaux
                    z[i] = plateaux

            # Hole in weir (slightly higher than inflow condition)
            if x[i] >= 0.3 and x[i] < 0.4 and y[i] > 0.2 and y[i] < 0.4:
                z[i] = -x[i]+self.inflow_stage + 0.02

            # Channel behind weir
            x0 = 0.5
            if x[i] >= 0.4 and x[i] < x0 and y[i] > 0.2 and y[i] < 0.4:
                z[i] = -x[i]+self.inflow_stage + 0.02

            if x[i] >= x0 and x[i] < 0.6 and y[i] > 0.2 and y[i] < 0.4:
                # Flatten it out towards the end
                z[i] = -x0+self.inflow_stage + 0.02 + (x0-x[i]) / 5

            # Hole to the east
            x0 = 1.1
            y0 = 0.35
            if num.sqrt((2*(x[i]-x0))**2 + (2*(y[i]-y0))**2) < 0.2:
                z[i] = num.sqrt(((x[i]-x0))**2 + ((y[i]-y0))**2)-1.0

            # Tiny channel draining hole
            if x[i] >= 1.14 and x[i] < 1.2 and y[i] >= 0.4 and y[i] < 0.6:
                z[i] = -0.9 # North south

            if x[i] >= 0.9 and x[i] < 1.18 and y[i] >= 0.58 and y[i] < 0.65:
                z[i] = -1.0 + (x[i]-0.9) / 3 # East west

            # Stuff not in use

            # Upward slope at inlet to the north west
            # if x[i] < 0.0: # and y[i] > 0.5:
            #    #z[i] = -y[i]+0.5  #-x[i]/2
            #    z[i] = x[i]/4 - y[i]**2 + 0.5

            # Hole to the west
            # x0 = -0.4; y0 = 0.35 # center
            # if sqrt((2*(x[i]-x0))**2 + (2*(y[i]-y0))**2) < 0.2:
            #    z[i] = sqrt(((x[i]-x0))**2 + ((y[i]-y0))**2)-0.2

        return z / 2


class Weir_simple(object):
    """Set a bathymetry for weir with a hole and a downstream gutter

    x,y are assumed to be in the unit square
    """

    def __init__(self, stage):
        self.inflow_stage = stage

    def __call__(self, x, y):
        N = len(x)
        assert N == len(y)

        z = num.zeros(N, float)
        for i in range(N):
            z[i] = -x[i]  # General slope

            # Flat bit to the left
            if x[i] < 0.3:
                z[i] = -x[i] / 10  # General slope

            # Weir
            if x[i] > 0.3 and x[i] < 0.4:
                z[i] = -x[i] + 0.9

            # Dip
            if x[i] > 0.6 and x[i] < 0.9:
                z[i] = -x[i] - 0.5  #-y[i]/5

            # Hole in weir (slightly higher than inflow condition)
            if x[i] > 0.3 and x[i] < 0.4 and y[i] > 0.2 and y[i] < 0.4:
                z[i] = -x[i] + self.inflow_stage + 0.05

        return z / 2



def scalar_func(t, x, y):
    """Function that returns a scalar.

    Used to test error message when numeric array is expected
    """

    return 17.7



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
        except:
            pass
        else:
            raise Exception('Should have raised an exception')

    # Individual flux tests
    def test_flux_zero_case(self):
        ql = num.zeros(3, float)
        qr = num.zeros(3, float)
        normal = num.zeros(2, float)
        edgeflux = num.zeros(3, float)
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
        edgeflux = num.zeros(3, float)
        zl = zr = 0.
        h = w - (zl+zr) / 2
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
        edgeflux = num.zeros(3, float)

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

        edgeflux = num.zeros(3, float)
        H0 = 0.0
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux,
                                  epsilon, g, H0)

        assert num.allclose(edgeflux, [-3., -20.0, -30.441])
        assert num.allclose(max_speed, 11.7146428199)

    def test_flux3(self):
        # Use data from previous version of abstract_2d_finite_volumes
        normal = num.array([-sqrt(2) / 2, sqrt(2) / 2])
        ql = num.array([-0.075, 2, 3])
        qr = num.array([-0.075, 2, 3])
        zl = zr = -0.375

        edgeflux = num.zeros(3, float)
        H0 = 0.0
        max_speed = flux_function(normal, ql, qr, zl, zr, edgeflux,
                                  epsilon, g, H0)

        assert num.allclose(edgeflux, [sqrt(2) / 2, 4.40221112, 7.3829019])
        assert num.allclose(max_speed, 4.0716654239)

    def test_flux4(self):
        # Use data from previous version of abstract_2d_finite_volumes
        normal = num.array([-sqrt(2) / 2, sqrt(2) / 2])
        ql = num.array([-0.34319278, 0.10254161, 0.07273855])
        qr = num.array([-0.30683287, 0.1071986, 0.05930515])
        zl = zr = -0.375

        edgeflux = num.zeros(3, float)
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
            return -x / 2

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
            assert name in domain.quantities

        assert num.alltrue(domain.get_conserved_quantities(0, edge=1) == 0.)



    def test_algorithm_parameters(self):
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
            assert name in domain.quantities



        parameters = domain.get_algorithm_parameters()




        assert parameters['minimum_allowed_height']  == domain.minimum_allowed_height
        assert parameters['maximum_allowed_speed']   == domain.maximum_allowed_speed
        assert parameters['minimum_storable_height'] == domain.minimum_storable_height
        assert parameters['g']                       == domain.g
        assert parameters['alpha_balance']           == domain.alpha_balance
        assert parameters['tight_slope_limiters']    == domain.tight_slope_limiters
        assert parameters['optimise_dry_cells']      == domain.optimise_dry_cells
        assert parameters['use_edge_limiter']        == domain.use_edge_limiter
        assert parameters['use_centroid_velocities'] == domain.use_centroid_velocities
        assert parameters['use_sloped_mannings']     == domain.use_sloped_mannings
        assert parameters['compute_fluxes_method']   == domain.get_compute_fluxes_method()
        assert parameters['flow_algorithm']          == domain.get_flow_algorithm()

        assert parameters['optimised_gradient_limiter']        == domain.optimised_gradient_limiter
        assert parameters['extrapolate_velocity_second_order'] == domain.extrapolate_velocity_second_order

    def test_institution(self):

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

        domain.set_boundary({'exterior': Reflective_boundary(domain)})

        institution = 'Test Institution'
        domain.set_institution(institution)

        #Evolution so as to create sww file
        for t in domain.evolve(yieldstep=0.5, finaltime=0.5):
            pass

        fid = NetCDFFile(domain.get_name() + '.sww')
        sww_institution = fid.institution
        fid.close()

        assert institution == sww_institution

        os.remove(domain.get_name() + '.sww')

    def test_boundary_conditions(self):
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #             bac,     bce,     ecf,     dbe
        vertices = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]
        boundary = {(0, 0): 'Third',
                    (0, 2): 'First',
                    (2, 0): 'Second',
                    (2, 1): 'Second',
                    (3, 1): 'Second',
                    (3, 2): 'Third'}

        domain = Domain(points, vertices, boundary)
        domain.check_integrity()

        domain.set_quantity('stage', [[1,2,3], [5,5,5], [0,0,9], [-6,3,3]])

        domain.set_quantity('xmomentum', [[1,1,1], [2,2,2], [3,3,3], [4,4,4]])

        domain.set_quantity('ymomentum', [[10,10,10], [20,20,20],
                                          [30,30,30], [40,40,40]])

        D = Dirichlet_boundary([5, 2, 1])
        T = Transmissive_boundary(domain)
        R = Reflective_boundary(domain)
        domain.set_boundary({'First': D, 'Second': T, 'Third': R})

        domain.update_boundary()

        # Stage
        assert domain.quantities['stage'].boundary_values[0] == 2.5
        # Reflective (2.5)
        assert (domain.quantities['stage'].boundary_values[0] ==
                domain.get_conserved_quantities(0, edge=0)[0])
        # Dirichlet
        assert domain.quantities['stage'].boundary_values[1] == 5.
        # Transmissive (4.5)
        assert (domain.quantities['stage'].boundary_values[2] ==
                domain.get_conserved_quantities(2, edge=0)[0])
        # Transmissive (4.5)
        assert (domain.quantities['stage'].boundary_values[3] ==
                domain.get_conserved_quantities(2, edge=1)[0])
        # Transmissive (-1.5)
        assert (domain.quantities['stage'].boundary_values[4] ==
                domain.get_conserved_quantities(3, edge=1)[0])
        # Reflective (-1.5)
        assert (domain.quantities['stage'].boundary_values[5] ==
                domain.get_conserved_quantities(3, edge=2)[0])

        # Xmomentum
        #print domain.quantities['xmomentum'].boundary_values
        #print domain.quantities['ymomentum'].boundary_values

        # Reflective
        assert domain.quantities['xmomentum'].boundary_values[0] == 1.0
        # Dirichlet
        assert domain.quantities['xmomentum'].boundary_values[1] == 2.
        # Transmissive
        assert (domain.quantities['xmomentum'].boundary_values[2] ==
                domain.get_conserved_quantities(2, edge=0)[1])
        # Transmissive
        assert (domain.quantities['xmomentum'].boundary_values[3] ==
                domain.get_conserved_quantities(2, edge=1)[1])
        # Transmissive
        assert (domain.quantities['xmomentum'].boundary_values[4] ==
                domain.get_conserved_quantities(3, edge=1)[1])
        # Reflective
        assert domain.quantities['xmomentum'].boundary_values[5] == -4.0

        # Ymomentum
        # Reflective
        assert domain.quantities['ymomentum'].boundary_values[0] == -10.0
        # Dirichlet
        assert domain.quantities['ymomentum'].boundary_values[1] == 1.
        # Transmissive
        assert domain.quantities['ymomentum'].boundary_values[2] == 30.
        # Transmissive
        assert domain.quantities['ymomentum'].boundary_values[3] == 30.
        # Transmissive
        assert domain.quantities['ymomentum'].boundary_values[4] == 40.
        # Reflective
        assert domain.quantities['ymomentum'].boundary_values[5] == 40.

    def test_boundary_conditionsII(self):
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #             bac,     bce,     ecf,     dbe
        vertices = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]
        boundary = {(0, 0): 'Third',
                    (0, 2): 'First',
                    (2, 0): 'Second',
                    (2, 1): 'Second',
                    (3, 1): 'Second',
                    (3, 2): 'Third',
                    (0, 1): 'Internal',
                    (1, 2): 'Internal'}

        domain = Domain(points, vertices, boundary)
        domain.check_integrity()

        domain.set_quantity('stage', [[1,2,3], [5,5,5], [0,0,9], [-6,3,3]])

        domain.set_quantity('xmomentum', [[1,1,1], [2,2,2], [3,3,3], [4,4,4]])

        domain.set_quantity('ymomentum', [[10,10,10], [20,20,20],
                                          [30,30,30], [40,40,40]])

        D = Dirichlet_boundary([5, 2, 1])
        T = Transmissive_boundary(domain)
        R = Reflective_boundary(domain)
        domain.set_boundary({'First': D, 'Second': T,
                             'Third': R, 'Internal': None})

        domain.update_boundary()
        domain.check_integrity()

    def test_boundary_conditions_inconsistent_internal(self):
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #             bac,     bce,     ecf,     dbe
        vertices = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]
        boundary = {(0, 0): 'Third',
                    (0, 2): 'First',
                    (2, 0): 'Second',
                    (2, 1): 'Second',
                    (3, 1): 'Second',
                    (3, 2): 'Third',
                    (0, 1): 'Internal'}

        domain = Domain(points, vertices, boundary)

        try:
            domain.check_integrity()
        except AssertionError:
            pass
        else:
            msg = 'should have thrown assertion exception'
            raise Exception(msg)





    def test_boundary_conditionsIII(self):
        """test_boundary_conditionsIII

        Test Transmissive_stage_zero_momentum_boundary
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
        boundary = {(0, 0): 'Third',
                    (0, 2): 'First',
                    (2, 0): 'Second',
                    (2, 1): 'Second',
                    (3, 1): 'Second',
                    (3, 2): 'Third'}

        domain = Domain(points, vertices, boundary)
        domain.check_integrity()

        domain.set_quantity('stage', [[1,2,3], [5,5,5], [0,0,9], [-6,3,3]])

        domain.set_quantity('xmomentum', [[1,1,1], [2,2,2], [3,3,3], [4,4,4]])

        domain.set_quantity('ymomentum', [[10,10,10], [20,20,20],
                                          [30,30,30], [40,40,40]])

        D = Dirichlet_boundary([5, 2, 1])
        T = Transmissive_stage_zero_momentum_boundary(domain)
        R = Reflective_boundary(domain)
        domain.set_boundary({'First': D, 'Second': T, 'Third': R})

        domain.update_boundary()

        # Stage
        # Reflective (2.5)
        assert domain.quantities['stage'].boundary_values[0] == 2.5
        assert (domain.quantities['stage'].boundary_values[0] ==
                domain.get_conserved_quantities(0, edge=0)[0])
        # Dirichlet
        assert domain.quantities['stage'].boundary_values[1] == 5.
        # Transmissive (4.5)
        assert (domain.quantities['stage'].boundary_values[2] ==
                domain.get_conserved_quantities(2, edge=0)[0])
        # Transmissive (4.5)
        assert (domain.quantities['stage'].boundary_values[3] ==
                domain.get_conserved_quantities(2, edge=1)[0])
        # Transmissive (-1.5)
        assert (domain.quantities['stage'].boundary_values[4] ==
                domain.get_conserved_quantities(3, edge=1)[0])
        # Reflective (-1.5)
        assert (domain.quantities['stage'].boundary_values[5] ==
                domain.get_conserved_quantities(3, edge=2)[0])

        # Xmomentum
        # Reflective
        assert domain.quantities['xmomentum'].boundary_values[0] == 1.0
        # Dirichlet
        assert domain.quantities['xmomentum'].boundary_values[1] == 2.
        assert domain.quantities['xmomentum'].boundary_values[2] == 0.0
        assert domain.quantities['xmomentum'].boundary_values[3] == 0.0
        assert domain.quantities['xmomentum'].boundary_values[4] == 0.0
        # Reflective
        assert domain.quantities['xmomentum'].boundary_values[5] == -4.0

        # Ymomentum
        # Reflective
        assert domain.quantities['ymomentum'].boundary_values[0] == -10.0
        # Dirichlet
        assert domain.quantities['ymomentum'].boundary_values[1] == 1.
        assert domain.quantities['ymomentum'].boundary_values[2] == 0.0
        assert domain.quantities['ymomentum'].boundary_values[3] == 0.0
        assert domain.quantities['ymomentum'].boundary_values[4] == 0.0
        # Reflective
        assert domain.quantities['ymomentum'].boundary_values[5] == 40.

    def test_boundary_condition_time(self):
        """test_boundary_condition_time

        This tests that boundary conditions are evaluated
        using the right time from domain.
        """

        # Setup computational domain
        from anuga.abstract_2d_finite_volumes.mesh_factory \
                import rectangular_cross

        #--------------------------------------------------------------
        # Setup computational domain
        #--------------------------------------------------------------
        N = 5
        points, vertices, boundary = rectangular_cross(N, N)
        domain = Domain(points, vertices, boundary)

        #--------------------------------------------------------------
        # Setup initial conditions
        #--------------------------------------------------------------
        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('friction', 0.0)
        domain.set_quantity('stage', 0.0)

        #--------------------------------------------------------------
        # Setup boundary conditions
        #--------------------------------------------------------------
        # Time dependent boundary
        Bt = Time_boundary(domain=domain, function=lambda t: [t, 0.0, 0.0])

        Br = Reflective_boundary(domain)              # Reflective wall

        domain.set_boundary({'left': Bt, 'right': Br, 'top': Br, 'bottom': Br})

        for t in domain.evolve(yieldstep = 10, finaltime = 20.0):
            q = Bt.evaluate()

            # FIXME (Ole): This test would not have passed in
            # changeset:5846.
            msg = 'Time boundary not evaluated correctly'
            assert num.allclose(t, q[0]), msg

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
        domain.set_quantity('stage', [[2,2,2], [2,2,2], [2,2,2], [2,2,2]])
        domain.check_integrity()

        assert num.allclose(domain.neighbours,
                            [[-1,1,-2], [2,3,0], [-3,-4,1],[1,-5,-6]])
        assert num.allclose(domain.neighbour_edges,
                            [[-1,2,-1], [2,0,1], [-1,-1,0],[1,-1,-1]])

        zl = zr = 0.     # Assume flat bed

        edgeflux = num.zeros(3, float)
        edgeflux0 = num.zeros(3, float)
        edgeflux1 = num.zeros(3, float)
        edgeflux2 = num.zeros(3, float)
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


        from anuga.shallow_water.shallow_water_ext import compute_fluxes_ext_central_structure as compute_fluxes
        # Now check that compute_flux yields zeros as well
        compute_fluxes(domain)

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
        domain.check_integrity()

        zl = zr = 0.    # Assume flat bed

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



        from anuga.shallow_water.shallow_water_ext import compute_fluxes_ext_central_structure as compute_fluxes
        # Now check that compute_flux yields zeros as well
        compute_fluxes(domain)


        for i, name in enumerate(['stage', 'xmomentum', 'ymomentum']):
            assert num.allclose(total_flux[i],
                                domain.quantities[name].explicit_update[1])

        assert num.allclose(domain.quantities['stage'].explicit_update,
                            [0., -0.68218178, -111.77316251, -35.68522449])
        assert num.allclose(domain.quantities['xmomentum'].explicit_update,
                            [-69.68888889, -166.6, 69.68888889, 0])
        assert num.allclose(domain.quantities['ymomentum'].explicit_update,
                            [-69.68888889, -35.93333333, 0., 69.68888889])

    def test_compute_fluxes_structure_2(self):
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
        edgeflux = num.zeros(3, float)
        edgeflux0 = num.zeros(3, float)
        edgeflux1 = num.zeros(3, float)
        edgeflux2 = num.zeros(3, float)
        H0 = 0.0

        domain.set_quantity('elevation', zl*num.ones((4, 3), int)) #array default#

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


        from anuga.shallow_water.shallow_water_ext import compute_fluxes_ext_central_structure as compute_fluxes
        # Now check that compute_flux yields zeros as well
        compute_fluxes(domain)

        for i, name in enumerate(['stage', 'xmomentum', 'ymomentum']):
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
        domain.set_quantity('elevation', zl*num.ones((4, 3), float)) #array default#

        edgeflux = num.zeros(3, float)
        edgeflux0 = num.zeros(3, float)
        edgeflux1 = num.zeros(3, float)
        edgeflux2 = num.zeros(3, float)
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


        from anuga.shallow_water.shallow_water_ext import compute_fluxes_ext_central_structure as compute_fluxes
        # Now check that compute_flux yields zeros as well
        flux_timestep = compute_fluxes(domain)

        #domain.compute_fluxes()

        #print domain.flux_timestep
        assert num.allclose(flux_timestep, 0.0426244319785)



        for i, name in enumerate(['stage', 'xmomentum', 'ymomentum']):
            #print total_flux[i]
            assert num.allclose(total_flux[i],
                                domain.quantities[name].explicit_update[1])


    def test_compute_fluxes_old_0(self):
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
                            [[-1,1,-2], [2,3,0], [-3,-4,1],[1,-5,-6]])
        assert num.allclose(domain.neighbour_edges,
                            [[-1,2,-1], [2,0,1], [-1,-1,0],[1,-1,-1]])

        zl = zr = 0.     # Assume flat bed

        edgeflux = num.zeros(3, float)
        edgeflux0 = num.zeros(3, float)
        edgeflux1 = num.zeros(3, float)
        edgeflux2 = num.zeros(3, float)
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

    def test_compute_fluxes_old_1(self):
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
        domain.set_flow_algorithm('1_5')

        val0 = 2. + 2.0/3
        val1 = 4. + 4.0/3
        val2 = 8. + 2.0/3
        val3 = 2. + 8.0/3

        domain.set_quantity('stage', [[val0, val0, val0], [val1, val1, val1],
                                      [val2, val2, val2], [val3, val3, val3]])
        domain.check_integrity()

        zl = zr = 0.    # Assume flat bed

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

    def test_compute_fluxes_old_2(self):
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

        domain.set_compute_fluxes_method('original')

        val0 = 2. + 2.0/3
        val1 = 4. + 4.0/3
        val2 = 8. + 2.0/3
        val3 = 2. + 8.0/3

        zl = zr = 0    # Assume flat zero bed
        edgeflux = num.zeros(3, float)
        edgeflux0 = num.zeros(3, float)
        edgeflux1 = num.zeros(3, float)
        edgeflux2 = num.zeros(3, float)
        H0 = 0.0

        domain.set_quantity('elevation', zl*num.ones((4, 3), int)) #array default#

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
    def test_compute_fluxes_old_3(self):
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

        domain.set_compute_fluxes_method('original')

        val0 = 2.+2.0/3
        val1 = 4.+4.0/3
        val2 = 8.+2.0/3
        val3 = 2.+8.0/3

        zl = zr = -3.75    # Assume constant bed (must be less than stage)
        domain.set_quantity('elevation', zl*num.ones((4, 3), int)) #array default#

        edgeflux = num.zeros(3, float)
        edgeflux0 = num.zeros(3, float)
        edgeflux1 = num.zeros(3, float)
        edgeflux2 = num.zeros(3, float)
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
        domain.set_quantity('elevation', zl*num.ones((4, 3), int)) #array default#
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
        domain.set_quantity('elevation', zl*num.ones((4, 3), int)) #array default#
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
        final_runup_height = -0.0166

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
            return -x / 2                             # linear bed slope

        # Use function for elevation
        domain.set_quantity('elevation', topography)
        domain.set_quantity('friction', 0.)                # Zero friction
        # Constant negative initial stage
        domain.set_quantity('stage', initial_runup_height)

        #--------------------------------------------------------------
        # Setup boundary conditions
        #--------------------------------------------------------------
        Br = Reflective_boundary(domain)                       # Reflective wall
        Bd = Dirichlet_boundary([final_runup_height, 0, 0])    # Constant inflow

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

    def test_get_flow_through_cross_section_with_geo(self):
        """test_get_flow_through_cross_section(self):

        Test that the total flow through a cross section can be
        correctly obtained at run-time from the ANUGA domain.

        This test creates a flat bed with a known flow through it and tests
        that the function correctly returns the expected flow.

        The specifics are
        e = -1 m
        u = 2 m/s
        h = 2 m
        w = 3 m (width of channel)

        q = u*h*w = 12 m^3/s

        This run tries it with georeferencing and with elevation = -1
        """

        # Create basic mesh (20m x 3m)
        width = 3
        length = 20
        t_end = 1
        points, vertices, boundary = rectangular(length, width, length, width)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary,
                        geo_reference=Geo_reference(56, 308500, 6189000))

        domain.default_order = 2
        domain.set_quantities_to_be_stored(None)

        e = -1.0
        w = 1.0
        h = w-e
        u = 2.0
        uh = u*h

        Br = Reflective_boundary(domain)     # Side walls
        Bd = Dirichlet_boundary([w, uh, 0])  # 2 m/s across the 3 m inlet:


        # Initial conditions
        domain.set_quantity('elevation', e)
        domain.set_quantity('stage', w)
        domain.set_quantity('xmomentum', uh)
        domain.set_boundary({'left': Bd, 'right': Bd, 'top': Br, 'bottom': Br})

        # Interpolation points down the middle
        I = [[0, width/2.],
             [length/2., width/2.],
             [length, width/2.]]
        interpolation_points = domain.geo_reference.get_absolute(I)

        # Shortcuts to quantites
        stage = domain.get_quantity('stage')
        xmomentum = domain.get_quantity('xmomentum')
        ymomentum = domain.get_quantity('ymomentum')

        for t in domain.evolve(yieldstep=0.1, finaltime=t_end):
            # Check that quantities are they should be in the interior
            w_t = stage.get_values(interpolation_points)
            uh_t = xmomentum.get_values(interpolation_points)
            vh_t = ymomentum.get_values(interpolation_points)

            assert num.allclose(w_t, w)
            assert num.allclose(uh_t, uh)
            assert num.allclose(vh_t, 0.0, atol=1.0e-6)

            # Check flows through the middle
            for i in range(5):
                x = length/2. + i*0.23674563    # Arbitrary
                cross_section = [[x, 0], [x, width]]

                cross_section = domain.geo_reference.get_absolute(cross_section)
                Q = domain.get_flow_through_cross_section(cross_section,
                                                          verbose=False)

                assert num.allclose(Q, uh*width)


    def test_get_energy_through_cross_section_with_geo(self):
        """test_get_energy_through_cross_section(self):

        Test that the total and specific energy through a cross section can be
        correctly obtained at run-time from the ANUGA domain.

        This test creates a flat bed with a known flow through it and tests
        that the function correctly returns the expected energies.

        The specifics are
        e = -1 m
        u = 2 m/s
        h = 2 m
        w = 3 m (width of channel)

        q = u*h*w = 12 m^3/s

        This run tries it with georeferencing and with elevation = -1
        """

        import time, os
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh (20m x 3m)
        width = 3
        length = 20
        t_end = 1
        points, vertices, boundary = rectangular(length, width, length, width)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary,
                        geo_reference=Geo_reference(56, 308500, 6189000))

        domain.default_order = 2
        domain.set_quantities_to_be_stored(None)

        e = -1.0
        w = 1.0
        h = w-e
        u = 2.0
        uh = u*h

        Br = Reflective_boundary(domain)       # Side walls
        Bd = Dirichlet_boundary([w, uh, 0])    # 2 m/s across the 3 m inlet:

        # Initial conditions
        domain.set_quantity('elevation', e)
        domain.set_quantity('stage', w)
        domain.set_quantity('xmomentum', uh)
        domain.set_boundary({'left': Bd, 'right': Bd, 'top': Br, 'bottom': Br})

        # Interpolation points down the middle
        I = [[0, width/2.],
             [length/2., width/2.],
             [length, width/2.]]
        interpolation_points = domain.geo_reference.get_absolute(I)

        # Shortcuts to quantites
        stage = domain.get_quantity('stage')
        xmomentum = domain.get_quantity('xmomentum')
        ymomentum = domain.get_quantity('ymomentum')

        for t in domain.evolve(yieldstep=0.1, finaltime=t_end):
            # Check that quantities are they should be in the interior
            w_t = stage.get_values(interpolation_points)
            uh_t = xmomentum.get_values(interpolation_points)
            vh_t = ymomentum.get_values(interpolation_points)

            assert num.allclose(w_t, w)
            assert num.allclose(uh_t, uh)
            assert num.allclose(vh_t, 0.0, atol=1.0e-6)

            # Check energies through the middle
            for i in range(5):
                x = length/2. + i*0.23674563    # Arbitrary
                cross_section = [[x, 0], [x, width]]

                cross_section = domain.geo_reference.get_absolute(cross_section)
                Es = domain.get_energy_through_cross_section(cross_section,
                                                             kind='specific',
                                                             verbose=False)

                assert num.allclose(Es, h + 0.5 * u * u / g)

                Et = domain.get_energy_through_cross_section(cross_section,
                                                             kind='total',
                                                             verbose=False)
                assert num.allclose(Et, w + 0.5 * u * u / g)


    def test_cross_section_class(self):
        """test_cross_section_class(self):

        Test that the total and specific energy through a cross section can be
        correctly obtained at run-time from the ANUGA cross section class.

        This test creates a flat bed with a known flow through it, creates a cross
        section and tests that the correct flow and energies are calculated

        The specifics are
        e = -1 m
        u = 2 m/s
        h = 2 m
        w = 3 m (width of channel)

        q = u*h*w = 12 m^3/s

        This run tries it with georeferencing and with elevation = -1
        """

        import time, os
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh (20m x 3m)
        width = 3
        length = 20
        t_end = 1
        points, vertices, boundary = rectangular(length, width, length, width)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary,
                        geo_reference=Geo_reference(56, 308500, 6189000))

        domain.default_order = 2
        domain.set_quantities_to_be_stored(None)

        e = -1.0
        w = 1.0
        h = w-e
        u = 2.0
        uh = u*h

        Br = Reflective_boundary(domain)       # Side walls
        Bd = Dirichlet_boundary([w, uh, 0])    # 2 m/s across the 3 m inlet:

        # Initial conditions
        domain.set_quantity('elevation', e)
        domain.set_quantity('stage', w)
        domain.set_quantity('xmomentum', uh)
        domain.set_boundary({'left': Bd, 'right': Bd, 'top': Br, 'bottom': Br})

        # Interpolation points down the middle
        I = [[0, width/2.],
             [length/2., width/2.],
             [length, width/2.]]
        interpolation_points = domain.geo_reference.get_absolute(I)

        # Shortcuts to quantites
        stage = domain.get_quantity('stage')
        xmomentum = domain.get_quantity('xmomentum')
        ymomentum = domain.get_quantity('ymomentum')


        # Create some cross sections
        cross_sections = []
        for i in range(5):
            x = length/2. + i*0.23674563    # Arbitrary
            polyline = [[x, 0], [x, width]]

            polyline = domain.geo_reference.get_absolute(polyline)

            cross_sections.append(Cross_section(domain,polyline))



        for t in domain.evolve(yieldstep=0.1, finaltime=t_end):
            # Check that quantities are they should be in the interior
            w_t = stage.get_values(interpolation_points)
            uh_t = xmomentum.get_values(interpolation_points)
            vh_t = ymomentum.get_values(interpolation_points)

            assert num.allclose(w_t, w)
            assert num.allclose(uh_t, uh)
            assert num.allclose(vh_t, 0.0, atol=1.0e-6)


            # Check flows and energies through the middle
            for cross_section in cross_sections:

                Q = cross_section.get_flow_through_cross_section()

                assert num.allclose(Q, uh*width)

                Es = cross_section.get_energy_through_cross_section(kind='specific')

                assert num.allclose(Es, h + 0.5*u*u / g)

                Et = cross_section.get_energy_through_cross_section(kind='total')

                assert num.allclose(Et, w + 0.5*u*u / g)






    def test_another_runup_example(self):
        """test_another_runup_example

        Test runup example where actual timeseries at interpolated
        points are tested.
        """

        from anuga.pmesh.mesh_interface import create_mesh_from_regions
        from anuga.abstract_2d_finite_volumes.mesh_factory \
                import rectangular_cross

        #-----------------------------------------------------------------
        # Setup computational domain
        #-----------------------------------------------------------------
        points, vertices, boundary = rectangular_cross(10, 10) # Basic mesh
        domain = Domain(points, vertices, boundary) # Create domain
        domain.set_low_froude(0)
        domain.set_default_order(2)
        domain.set_quantities_to_be_stored(None)
        domain.H0 = 1.0e-3

        #-----------------------------------------------------------------
        # Setup initial conditions
        #-----------------------------------------------------------------
        def topography(x, y):
            return -x / 2                              # linear bed slope

        domain.set_quantity('elevation', topography)
        domain.set_quantity('friction', 0.0)
        domain.set_quantity('stage', expression='elevation')

        #----------------------------------------------------------------
        # Setup boundary conditions
        #----------------------------------------------------------------
        Br = Reflective_boundary(domain)           # Solid reflective wall
        Bd = Dirichlet_boundary([-0.2, 0., 0.])    # Constant boundary values
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

        # Captured data from code manually inspected for correctness 11/5/2010


#         import pprint
#         pprint.pprint(gauge_values[0])
#         pprint.pprint(gauge_values[1])
#         pprint.pprint(gauge_values[2])
#         pprint.pprint(gauge_values[3])

        # Steve Note: Had to recapture these values when changed to default
        # flow algorithm to DE0

        G0 = [-0.19166666666666665,
             -0.19166666666666665,
             -0.19166666666666665,
             -0.19166666666666665,
             -0.19166666666666665,
             -0.17357498789004924,
             -0.16134698833835073,
             -0.15299819808948953,
             -0.15658886028668945,
             -0.15619519506566443,
             -0.15983820234428089,
             -0.17201558491115593,
             -0.18809602362873767,
             -0.19790107878825061,
             -0.19916520934913592,
             -0.19897656121883669,
             -0.1994526052093108,
             -0.19959391223111991,
             -0.19953478511280054,
             -0.19972300428585374,
             -0.19630771510281139,
             -0.19142529806554692,
             -0.19150376369159233,
             -0.19155465935100272,
             -0.1915480843183047,
             -0.19158487652868542,
             -0.19160491229801035,
             -0.19161005537288817,
             -0.19161941600794635,
             -0.19162712306975391,
             -0.19163250820388589,
             -0.19163717575751205,
             -0.19164123614810449,
             -0.19164461943275748,
             -0.19164746068339567,
             -0.19164986883200527,
             -0.19165191683640856,
             -0.19165366545685461,
             -0.19165516248134548,
             -0.19165644941311713,
             -0.19165755970877371,
             -0.19165852034152259,
             -0.19165935411364624,
             -0.19166008027417775,
             -0.19166071302904672,
             -0.19166126670954906,
             -0.19166176113272784,
             -0.19166219795907555,
             -0.1916625893478662,
             -0.1916629307243472,
             -0.19166323917681413]

        G1 = [-0.29166666666666669,
             -0.29166666666666669,
             -0.29166666666666669,
             -0.2625120524489587,
             -0.23785652592777537,
             -0.22355577785898192,
             -0.21194260352756192,
             -0.20074023971177196,
             -0.1921307806774481,
             -0.18446465937291226,
             -0.17997582762982173,
             -0.17711270247966202,
             -0.17851892444764331,
             -0.18506971624496477,
             -0.19575447168307283,
             -0.20276874461738162,
             -0.20604869378667856,
             -0.20679384168209414,
             -0.20641078072679761,
             -0.20542266250262109,
             -0.2040758427409933,
             -0.2023210482025932,
             -0.20079135595686456,
             -0.200151095127519,
             -0.19933359856480387,
             -0.1987496662304945,
             -0.19838824488166795,
             -0.19840968868426581,
             -0.19875720118469065,
             -0.1992908187455468,
             -0.19978498430904629,
             -0.20013379308921292,
             -0.20033977063045771,
             -0.20042440892886232,
             -0.20040623053779325,
             -0.20034473053740201,
             -0.20025342290445922,
             -0.20014447072563055,
             -0.20003674497651436,
             -0.19995754677956581,
             -0.19990776759716164,
             -0.19988543049560956,
             -0.199888928325714,
             -0.19990718400594476,
             -0.19993251283774235,
             -0.19996431619318056,
             -0.19999093834815995,
             -0.20001191340810734,
             -0.20002510023168069,
             -0.20003044891380231,
             -0.20002952773077484]

        G2 = [-0.42499999999999999,
             -0.40953908664730287,
             -0.33296103674588662,
             -0.30451769824676844,
             -0.28219604345783056,
             -0.26522350354865254,
             -0.25015947587031895,
             -0.23608529438075876,
             -0.22253484154746356,
             -0.20994131123668461,
             -0.19973743965316049,
             -0.19049780733818705,
             -0.18495761801075922,
             -0.18037302557409396,
             -0.18191116803107951,
             -0.18821767150886343,
             -0.19501674197220067,
             -0.19950705780757777,
             -0.20218620150235145,
             -0.20347437823613723,
             -0.20399600876175172,
             -0.20401809537321336,
             -0.20338240306788496,
             -0.20225546137366332,
             -0.20125582827786589,
             -0.20064315945485711,
             -0.19994146252008138,
             -0.19933940195992941,
             -0.19890775656150686,
             -0.19881040880955592,
             -0.19901575680457334,
             -0.19936426014319533,
             -0.1996865767169046,
             -0.19994731519882072,
             -0.20013715580150337,
             -0.20024380974598521,
             -0.20027966500003219,
             -0.20026916410031764,
             -0.20021596069804154,
             -0.20014200328898363,
             -0.20006589807932401,
             -0.20000079735581283,
             -0.19995360469070572,
             -0.19992526935599028,
             -0.19991808248161369,
             -0.19992443120735207,
             -0.19994339357600277,
             -0.19996469380466553,
             -0.19998467955016122,
             -0.20000213621634999,
             -0.20001413684386354]

        G3 = [-0.44166666666666665,
             -0.36356820384297261,
             -0.32677493201134605,
             -0.30732284916772212,
             -0.29038753972867332,
             -0.27270481540957403,
             -0.25782191624274975,
             -0.2432498394927349,
             -0.22936799633074165,
             -0.2163115182211087,
             -0.20471162972237375,
             -0.19474955403226407,
             -0.18722605352323582,
             -0.18186433340685534,
             -0.18066168631825957,
             -0.1849738852801476,
             -0.19187415350512246,
             -0.19744289107562996,
             -0.20094758082578065,
             -0.2028743338186943,
             -0.20375298105641421,
             -0.20404634370333788,
             -0.20385967844158176,
             -0.20281299387611038,
             -0.20172106447516341,
             -0.20097648698314996,
             -0.20029217955013084,
             -0.19962005657018717,
             -0.19908929111834933,
             -0.19882323339870062,
             -0.19886705890507497,
             -0.19919177307550476,
             -0.1995284015021406,
             -0.19982374500648564,
             -0.20005201207002615,
             -0.20019729220453644,
             -0.20026700514017221,
             -0.20027541321283041,
             -0.20024991859335209,
             -0.20018037785313814,
             -0.20010282964013762,
             -0.20003123082670529,
             -0.19997471558256363,
             -0.19993722523364474,
             -0.19992030412618081,
             -0.19992065759609112,
             -0.19993226086341562,
             -0.19995375100532778,
             -0.19997513774473216,
             -0.19999415800173256,
             -0.20000875071328705]





        assert num.allclose(gauge_values[0], G0)
        assert num.allclose(gauge_values[1], G1)
        assert num.allclose(gauge_values[2], G2)
        assert num.allclose(gauge_values[3], G3)

    #####################################################

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

    def test_gravity(self):
        #Assuming no friction

        from anuga.config import g

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
            return slope(x,y) + h

        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage)

        for name in domain.conserved_quantities:
            assert num.allclose(domain.quantities[name].explicit_update, 0)
            assert num.allclose(domain.quantities[name].semi_implicit_update, 0)

        #domain.compute_forcing_terms()
        from anuga.shallow_water.shallow_water_ext import gravity
        gravity(domain)

        #print domain.quantities['xmomentum'].explicit_update
        #print domain.quantities['ymomentum'].explicit_update


        assert num.allclose(domain.quantities['stage'].explicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].explicit_update,
                            -g*h*3)
        assert num.allclose(domain.quantities['ymomentum'].explicit_update, 0)


    def test_gravity_2(self):
        #Assuming no friction

        from anuga.config import g

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

        h = 15
        def stage(x, y):
            return h

        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage)

        for name in domain.conserved_quantities:
            assert num.allclose(domain.quantities[name].explicit_update, 0)
            assert num.allclose(domain.quantities[name].semi_implicit_update, 0)

        #domain.compute_forcing_terms()
        from anuga.shallow_water.shallow_water_ext import gravity
        gravity(domain)

        #print domain.quantities['xmomentum'].explicit_update
        #print domain.quantities['ymomentum'].explicit_update


        assert num.allclose(domain.quantities['stage'].explicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].explicit_update,
                            [-382.2, -323.4, -205.8, -382.2])
        assert num.allclose(domain.quantities['ymomentum'].explicit_update, 0)



    def test_well_balanced_flux_1(self):
        #testing that compute_fluxes_method = wb_1 is well balance

        from anuga.config import g
        from anuga import Reflective_boundary

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

        h = 15
        def stage(x, y):
            return h

        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage)



        for name in domain.conserved_quantities:
            assert num.allclose(domain.quantities[name].explicit_update, 0)
            assert num.allclose(domain.quantities[name].semi_implicit_update, 0)


        Br = Reflective_boundary(domain)      # Solid reflective wall
        domain.set_boundary({'exterior' :Br})
        domain.update_boundary()

        domain.set_compute_fluxes_method('wb_1')
        domain.compute_fluxes()


        #print domain.quantities['xmomentum'].explicit_update
        #print domain.quantities['stage'].vertex_values
        #print domain.quantities['elevation'].vertex_values
        #print domain.quantities['ymomentum'].explicit_update


        assert num.allclose(domain.quantities['stage'].explicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].explicit_update, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].explicit_update, 0.0)


    def test_well_balanced_flux_2(self):
        #esting that compute_fluxes_method = wb_2 is well balance

        from anuga.config import g
        from anuga import Reflective_boundary

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

        h = 15
        def stage(x, y):
            return h

        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage)



        for name in domain.conserved_quantities:
            assert num.allclose(domain.quantities[name].explicit_update, 0)
            assert num.allclose(domain.quantities[name].semi_implicit_update, 0)


        Br = Reflective_boundary(domain)      # Solid reflective wall
        domain.set_boundary({'exterior' :Br})
        domain.update_boundary()

        domain.set_compute_fluxes_method('wb_2')
        domain.compute_fluxes()


        #print domain.quantities['xmomentum'].explicit_update
        #print domain.quantities['stage'].vertex_values
        #print domain.quantities['elevation'].vertex_values
        #print domain.quantities['ymomentum'].explicit_update


        assert num.allclose(domain.quantities['stage'].explicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].explicit_update, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].explicit_update, 0.0)

    def test_gravity_wb(self):
        #Assuming no friction

        from anuga.config import g

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
            return slope(x,y)+h

        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage)

        for name in domain.conserved_quantities:
            assert num.allclose(domain.quantities[name].explicit_update, 0)
            assert num.allclose(domain.quantities[name].semi_implicit_update, 0)

        from anuga.shallow_water.shallow_water_ext import gravity_wb
        gravity_wb(domain)


        #print domain.quantities['xmomentum'].explicit_update
        #print domain.quantities['ymomentum'].explicit_update


        assert num.allclose(domain.quantities['stage'].explicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].explicit_update,
                            -g*h*3)
        assert num.allclose(domain.quantities['ymomentum'].explicit_update, 0)


    def test_gravity_wb_2(self):
        #Assuming no friction

        from anuga.config import g

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

        h = 15.0
        def stage(x, y):
            return h

        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage)

        for name in domain.conserved_quantities:
            assert num.allclose(domain.quantities[name].explicit_update, 0)
            assert num.allclose(domain.quantities[name].semi_implicit_update, 0)

        from anuga.shallow_water.shallow_water_ext import gravity_wb
        gravity_wb(domain)


        #print domain.quantities['xmomentum'].explicit_update
        #print domain.quantities['ymomentum'].explicit_update

        assert num.allclose(domain.quantities['stage'].explicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].explicit_update,
                            [-396.9, -308.7, -220.5, -396.9])
        assert num.allclose(domain.quantities['ymomentum'].explicit_update, 0)



    def test_compute_fluxes_flat_wb_1(self):
        #Assuming no friction

        from anuga.config import g

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
        B = Reflective_boundary(domain)
        domain.set_boundary( {'exterior': B})
        domain.set_compute_fluxes_method('wb_1')

        #Set up for a gradient of (3,0) at mid triangle (bce)
        def slope(x, y):
            return 3*x

        h = 15.0
        def stage(x, y):
            return h

        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage)

        for name in domain.conserved_quantities:
            assert num.allclose(domain.quantities[name].explicit_update, 0)
            assert num.allclose(domain.quantities[name].semi_implicit_update, 0)

        domain.update_boundary()
        domain.compute_fluxes()


        #print domain.quantities['xmomentum'].explicit_update
        #print domain.quantities['ymomentum'].explicit_update

        assert num.allclose(domain.quantities['stage'].explicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].explicit_update, 0)
        assert num.allclose(domain.quantities['ymomentum'].explicit_update, 0)


    def test_compute_fluxes_slope_wb_1(self):
        #Assuming no friction

        from anuga.config import g

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
        B = Reflective_boundary(domain)
        domain.set_boundary( {'exterior': B})
        domain.set_compute_fluxes_method('wb_1')

        #Set up for a gradient of (3,0) at mid triangle (bce)
        def slope(x, y):
            return 3*x

        h = 0.1
        def stage(x, y):
            return slope(x,y)+h

        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage)

        for name in domain.conserved_quantities:
            assert num.allclose(domain.quantities[name].explicit_update, 0)
            assert num.allclose(domain.quantities[name].semi_implicit_update, 0)

        domain.update_boundary()
        domain.compute_fluxes()


        #print domain.quantities['xmomentum'].explicit_update
        #print domain.quantities['ymomentum'].explicit_update

        assert num.allclose(domain.quantities['stage'].explicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].explicit_update,
                                [ -2.94,  -2.94, -10.29, -10.29])
        assert num.allclose(domain.quantities['ymomentum'].explicit_update,
                                [  7.35, 0.0, 0.0, -7.35])



    def test_compute_fluxes_flat_wb_2(self):
        #Assuming no friction

        from anuga.config import g

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
        B = Reflective_boundary(domain)
        domain.set_boundary( {'exterior': B})
        domain.set_compute_fluxes_method('wb_2')

        #Set up for a gradient of (3,0) at mid triangle (bce)
        def slope(x, y):
            return 3*x

        h = 15.0
        def stage(x, y):
            return h

        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage)

        for name in domain.conserved_quantities:
            assert num.allclose(domain.quantities[name].explicit_update, 0)
            assert num.allclose(domain.quantities[name].semi_implicit_update, 0)

        domain.update_boundary()
        domain.compute_fluxes()


        #print domain.quantities['xmomentum'].explicit_update
        #print domain.quantities['ymomentum'].explicit_update

        assert num.allclose(domain.quantities['stage'].explicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].explicit_update, 0)
        assert num.allclose(domain.quantities['ymomentum'].explicit_update, 0)




    def test_compute_fluxes_flat_wb_3(self):
        #Assuming no friction

        from anuga.config import g

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
        B = Reflective_boundary(domain)
        domain.set_boundary( {'exterior': B})
        domain.set_compute_fluxes_method('wb_3')

        #Set up for a gradient of (3,0) at mid triangle (bce)
        def slope(x, y):
            return 3*x

        h = 15.0
        def stage(x, y):
            return h

        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage)

        for name in domain.conserved_quantities:
            assert num.allclose(domain.quantities[name].explicit_update, 0)
            assert num.allclose(domain.quantities[name].semi_implicit_update, 0)

        domain.update_boundary()
        domain.compute_fluxes()


        #print 'should be 0 ', domain.quantities['xmomentum'].explicit_update
        #print 'should be 0 ', domain.quantities['ymomentum'].explicit_update

        assert num.allclose(domain.quantities['stage'].explicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].explicit_update, 0)
        assert num.allclose(domain.quantities['ymomentum'].explicit_update, 0)

    def test_compute_fluxes_slope_wb_3(self):
        #Assuming no friction

        from anuga.config import g

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
        B = Reflective_boundary(domain)
        domain.set_boundary( {'exterior': B})
        domain.set_compute_fluxes_method('wb_3')

        #Set up for a gradient of (3,0) at mid triangle (bce)
        def slope(x, y):
            return 3*x

        h = 0.1
        def stage(x, y):
            return slope(x,y)+h

        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage)

        for name in domain.conserved_quantities:
            assert num.allclose(domain.quantities[name].explicit_update, 0)
            assert num.allclose(domain.quantities[name].semi_implicit_update, 0)

        domain.update_boundary()
        domain.compute_fluxes()


        #print domain.quantities['xmomentum'].explicit_update
        #print domain.quantities['ymomentum'].explicit_update

        assert num.allclose(domain.quantities['stage'].explicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].explicit_update,
                                [ -2.94,  -2.94, -2.94, -2.94])
        assert num.allclose(domain.quantities['ymomentum'].explicit_update,
                                [  0.0, 0.0, 0.0, 0.0])


    def test_manning_friction(self):
        """ Assuming flat manning frinction is default
        """
        from anuga.config import g

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

        domain.set_sloped_mannings_function(False)

        B = Reflective_boundary(domain)
        domain.set_boundary( {'exterior': B})

        #Set up for a gradient of (3,0) at mid triangle (bce)
        def slope(x, y):
            return 3*x

        h = 0.1
        def stage(x, y):
            return slope(x, y) + h

        eta = 0.07
        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage)
        domain.set_quantity('friction', eta)

        for name in domain.conserved_quantities:
            assert num.allclose(domain.quantities[name].explicit_update, 0)
            assert num.allclose(domain.quantities[name].semi_implicit_update, 0)

        domain.compute_forcing_terms()

        assert num.allclose(domain.quantities['stage'].explicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].explicit_update,
                            0)
        assert num.allclose(domain.quantities['ymomentum'].explicit_update, 0)

        assert num.allclose(domain.quantities['stage'].semi_implicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].semi_implicit_update,
                            0)
        assert num.allclose(domain.quantities['ymomentum'].semi_implicit_update,
                            0)

        #Create some momentum for friction to work with
        domain.set_quantity('xmomentum', 1)
        S = -g*eta**2 /  h**(7.0/3)

        domain.compute_forcing_terms()
        assert num.allclose(domain.quantities['stage'].semi_implicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].semi_implicit_update,
                            S)
        assert num.allclose(domain.quantities['ymomentum'].semi_implicit_update,
                            0)

        #A more complex example
        domain.quantities['stage'].semi_implicit_update[:] = 0.0
        domain.quantities['xmomentum'].semi_implicit_update[:] = 0.0
        domain.quantities['ymomentum'].semi_implicit_update[:] = 0.0

        domain.set_quantity('xmomentum', 3)
        domain.set_quantity('ymomentum', 4)

        S = -g*eta**2*5 /  h**(7.0/3)

        domain.compute_forcing_terms()

        assert num.allclose(domain.quantities['stage'].semi_implicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].semi_implicit_update,
                            3*S)
        assert num.allclose(domain.quantities['ymomentum'].semi_implicit_update,
                            4*S)



    def test_flat_manning_friction(self):
        from anuga.config import g

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
        B = Reflective_boundary(domain)
        domain.set_boundary( {'exterior': B})

        # Use the flat function which doesn't takes into account the extra
        # wetted area due to slope of bed
        domain.set_sloped_mannings_function(False)

        #Set up for a gradient of (3,0) at mid triangle (bce)
        def slope(x, y):
            return 3*x

        h = 0.1
        def stage(x, y):
            return slope(x, y) + h

        eta = 0.07
        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage)
        domain.set_quantity('friction', eta)

        for name in domain.conserved_quantities:
            assert num.allclose(domain.quantities[name].explicit_update, 0)
            assert num.allclose(domain.quantities[name].semi_implicit_update, 0)

        domain.compute_forcing_terms()

        assert num.allclose(domain.quantities['stage'].explicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].explicit_update,
                            0)
        assert num.allclose(domain.quantities['ymomentum'].explicit_update, 0)

        assert num.allclose(domain.quantities['stage'].semi_implicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].semi_implicit_update,
                            0)
        assert num.allclose(domain.quantities['ymomentum'].semi_implicit_update,
                            0)

        #Create some momentum for friction to work with
        domain.set_quantity('xmomentum', 1)
        S = -g*eta**2 /  h**(7.0/3)

        domain.compute_forcing_terms()
        assert num.allclose(domain.quantities['stage'].semi_implicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].semi_implicit_update,
                            S)
        assert num.allclose(domain.quantities['ymomentum'].semi_implicit_update,
                            0)

        #A more complex example
        domain.quantities['stage'].semi_implicit_update[:] = 0.0
        domain.quantities['xmomentum'].semi_implicit_update[:] = 0.0
        domain.quantities['ymomentum'].semi_implicit_update[:] = 0.0

        domain.set_quantity('xmomentum', 3)
        domain.set_quantity('ymomentum', 4)

        S = -g*eta**2*5 /  h**(7.0/3)

        domain.compute_forcing_terms()

        assert num.allclose(domain.quantities['stage'].semi_implicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].semi_implicit_update,
                            3*S)
        assert num.allclose(domain.quantities['ymomentum'].semi_implicit_update,
                            4*S)


    def test_sloped_manning_friction(self):
        from anuga.config import g

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
        B = Reflective_boundary(domain)
        domain.set_boundary( {'exterior': B})

        # Use the sloped function which takes into account the extra
        # wetted area due to slope of bed
        domain.set_sloped_mannings_function(True)

        #Set up for a gradient of (3,0) at mid triangle (bce)
        def slope(x, y):
            return 3*x

        h = 0.1
        def stage(x, y):
            return slope(x, y) + h

        eta = 0.07
        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage)
        domain.set_quantity('friction', eta)

        for name in domain.conserved_quantities:
            assert num.allclose(domain.quantities[name].explicit_update, 0)
            assert num.allclose(domain.quantities[name].semi_implicit_update, 0)

        domain.compute_forcing_terms()

        assert num.allclose(domain.quantities['stage'].explicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].explicit_update,
                            0)
        assert num.allclose(domain.quantities['ymomentum'].explicit_update, 0)

        assert num.allclose(domain.quantities['stage'].semi_implicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].semi_implicit_update,
                            0)
        assert num.allclose(domain.quantities['ymomentum'].semi_implicit_update,
                            0)

        #Create some momentum for friction to work with
        domain.set_quantity('xmomentum', 1)
        S = -g*eta**2 /  h**(7.0/3) * sqrt(10)

        domain.compute_forcing_terms()
        assert num.allclose(domain.quantities['stage'].semi_implicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].semi_implicit_update,
                            S)
        assert num.allclose(domain.quantities['ymomentum'].semi_implicit_update,
                            0)

        #A more complex example
        domain.quantities['stage'].semi_implicit_update[:] = 0.0
        domain.quantities['xmomentum'].semi_implicit_update[:] = 0.0
        domain.quantities['ymomentum'].semi_implicit_update[:] = 0.0

        domain.set_quantity('xmomentum', 3)
        domain.set_quantity('ymomentum', 4)

        S = -g*eta**2*5 /  h**(7.0/3) * sqrt(10.0)

        domain.compute_forcing_terms()

        assert num.allclose(domain.quantities['stage'].semi_implicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].semi_implicit_update,
                            3*S)
        assert num.allclose(domain.quantities['ymomentum'].semi_implicit_update,
                            4*S)



    def test_inflow_using_circle(self):
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
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        # Setup only one forcing term, constant inflow of 2 m^3/s
        # on a circle affecting triangles #0 and #1 (bac and bce)
        domain.forcing_terms = []

        I = Inflow(domain, rate=2.0, center=(1,1), radius=1)
        domain.forcing_terms.append(I)
        domain.compute_forcing_terms()


        A = I.exchange_area
        assert num.allclose(A, 4) # Two triangles

        assert num.allclose(domain.quantities['stage'].explicit_update[1], 2.0/A)
        assert num.allclose(domain.quantities['stage'].explicit_update[0], 2.0/A)
        assert num.allclose(domain.quantities['stage'].explicit_update[2:], 0)


    def test_inflow_using_circle_function(self):
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
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        # Setup only one forcing term, time dependent inflow of 2 m^3/s
        # on a circle affecting triangles #0 and #1 (bac and bce)
        domain.forcing_terms = []
        I = Inflow(domain, rate=lambda t: 2., center=(1,1), radius=1)
        domain.forcing_terms.append(I)

        domain.compute_forcing_terms()

        A = I.exchange_area
        assert num.allclose(A, 4) # Two triangles

        assert num.allclose(domain.quantities['stage'].explicit_update[1], 2.0/A)
        assert num.allclose(domain.quantities['stage'].explicit_update[0], 2.0/A)
        assert num.allclose(domain.quantities['stage'].explicit_update[2:], 0)




    def test_inflow_catch_too_few_triangles(self):
        """
        Test that exception is thrown if no triangles are covered
        by the inflow area
        """

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
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        # Setup only one forcing term, constant inflow of 2 m^3/s
        # on a circle affecting triangles #0 and #1 (bac and bce)
        try:
            Inflow(domain, rate=2.0, center=(1,1.1), radius=0.01)
        except:
            pass
        else:
            msg = 'Should have raised exception'
            raise Exception(msg)

    def Xtest_inflow_outflow_conservation(self):
        """
        Test what happens if water is abstracted from one area and
        injected into another - especially if there is not enough
        water to match the abstraction.
        This tests that the total volume is kept constant under a range of
        scenarios.

        This test will fail as the problem was only fixed for culverts.
        """

        from math import pi, cos, sin

        length = 20.
        width = 10.

        dx = dy = 2    # 1 or 2 OK
        points, vertices, boundary = rectangular_cross(int(length / dx),
                                                       int(width / dy),
                                                       len1=length,
                                                       len2=width)
        domain = Domain(points, vertices, boundary)
        domain.set_name('test_inflow_conservation')    # Output name
        domain.set_default_order(2)

        # Flat surface with 1m of water
        stage = 1.0
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', stage)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'bottom': Br, 'top': Br})

        # Setup one forcing term, constant inflow of 2 m^3/s on a circle
        domain.forcing_terms = []
        domain.forcing_terms.append(Inflow(domain, rate=2.0,
                                           center=(5,5), radius=1))

        domain.compute_forcing_terms()

        # Check that update values are correct
        for x in domain.quantities['stage'].explicit_update:
            assert num.allclose(x, 2.0/pi) or num.allclose(x, 0.0)

        # Check volumes without inflow
        domain.forcing_terms = []
        initial_volume = domain.quantities['stage'].get_integral()

        assert num.allclose(initial_volume, width*length*stage)

        for t in domain.evolve(yieldstep = 0.05, finaltime = 5.0):
            volume = domain.quantities['stage'].get_integral()
            assert num.allclose(volume, initial_volume)

        # Now apply the inflow and check volumes for a range of stage values
        for stage in [2.0, 1.0, 0.5, 0.25, 0.1, 0.0]:
            domain.set_time(0.0)
            domain.set_quantity('stage', stage)
            domain.forcing_terms = []
            domain.forcing_terms.append(Inflow(domain, rate=2.0,
                                               center=(5,5), radius=1))
            initial_volume = domain.quantities['stage'].get_integral()
            predicted_volume = initial_volume
            dt = 0.05
            for t in domain.evolve(yieldstep=dt, finaltime=5.0):
                volume = domain.quantities['stage'].get_integral()
                assert num.allclose (volume, predicted_volume)
                predicted_volume = predicted_volume + 2.0 / pi / 100 / dt # Why 100?

        # Apply equivalent outflow only and check volumes
        # for a range of stage values
        for stage in [2.0, 1.0, 0.5, 0.25, 0.1, 0.0]:
            #print stage

            domain.set_time(0.0)
            domain.set_quantity('stage', stage)
            domain.forcing_terms = []
            domain.forcing_terms.append(Inflow(domain, rate=-2.0,
                                               center=(15,5), radius=1))
            initial_volume = domain.quantities['stage'].get_integral()
            predicted_volume = initial_volume
            dt = 0.05
            for t in domain.evolve(yieldstep=dt, finaltime=5.0):
                volume = domain.quantities['stage'].get_integral()
                print(t, volume, predicted_volume)
                assert num.allclose (volume, predicted_volume)
                predicted_volume = predicted_volume - 2.0 / pi / 100 / dt # Why 100?

        # Apply both inflow and outflow and check volumes being constant for a
        # range of stage values
        for stage in [2.0, 1.0, 0.5, 0.25, 0.1, 0.0]:
            #print stage

            domain.set_time(0.0)
            domain.set_quantity('stage', stage)
            domain.forcing_terms = []
            domain.forcing_terms.append(Inflow(domain, rate=2.0,
                                               center=(5,5), radius=1))
            domain.forcing_terms.append(Inflow(domain, rate=-2.0,
                                               center=(15,5), radius=1))
            initial_volume = domain.quantities['stage'].get_integral()

            dt = 0.05
            for t in domain.evolve(yieldstep=dt, finaltime=5.0):
                volume = domain.quantities['stage'].get_integral()

                #print t, volume
                assert num.allclose(volume, initial_volume)

    #####################################################

    def test_first_order_extrapolator_const_z(self):
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

        zl = zr = -3.75    # Assume constant bed (must be less than stage)
        domain.set_quantity('elevation', zl*num.ones((4, 3), int)) #array default#
        domain.set_quantity('stage', [[val0, val0-1, val0-2],
                                      [val1, val1+1, val1],
                                      [val2, val2-2, val2],
                                      [val3-0.5, val3, val3]])

        domain._order_ = 1
        domain.distribute_to_vertices_and_edges()

        #Check that centroid values were distributed to vertices
        C = domain.quantities['stage'].centroid_values
        V = num.sum(domain.quantities['stage'].vertex_values, axis=1)/3.0

        assert num.allclose(V,C)


    def test_first_order_limiter_variable_z(self):
        '''Check that first order limiter follows bed_slope'''

        from anuga.config import epsilon

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]
        #bac, bce, ecf, dbe
        vertices = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices)
        val0 = 2.+2.0/3
        val1 = 4.+4.0/3
        val2 = 8.+2.0/3
        val3 = 2.+8.0/3

        domain.set_quantity('elevation', [[0,0,0], [6,0,0],
                                          [6,6,6], [6,6,6]])
        domain.set_quantity('stage', [[val0, val0, val0],
                                      [val1, val1, val1],
                                      [val2, val2, val2],
                                      [val3, val3, val3]])

        E = domain.quantities['elevation'].vertex_values
        L = domain.quantities['stage'].vertex_values


        #Check that some stages are not above elevation (within eps)
        #- so that the limiter has something to work with
        assert not num.alltrue(num.alltrue(num.greater_equal(L,E-epsilon)))

        domain._order_ = 1
        domain.distribute_to_vertices_and_edges()

        #Check that all stages are above elevation (within eps)
        assert num.alltrue(num.alltrue(num.greater_equal(L,E-epsilon)))

    def test_distribute_basic(self):
        #Using test data generated by abstract_2d_finite_volumes-2
        #Assuming no friction and flat bed (0.0)

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

        val0 = 2.
        val1 = 4.
        val2 = 8.
        val3 = 2.

        domain.set_quantity('stage', [val0, val1, val2, val3],
                            location='centroids')
        L = domain.quantities['stage'].vertex_values

        # First order
        domain._order_ = 1
        domain.distribute_to_vertices_and_edges()

        assert num.allclose(L[1], [ 0.,  6.,  6.])
        assert num.allclose(num.sum(L[1]), 3*val1)

        # Second order
        domain._order_ = 2
        domain.beta_w = 0.9
        domain.beta_w_dry = 0.9
        domain.beta_uh = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh = 0.9
        domain.beta_vh_dry = 0.9
        domain.distribute_to_vertices_and_edges()
        assert num.allclose(L[1], [ 0.,  6.,  6.])
        assert num.allclose(num.sum(L[1]), 3*val1)


    def test_distribute_away_from_bed(self):
        #Using test data generated by abstract_2d_finite_volumes-2
        #Assuming no friction and flat bed (0.0)

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
        #domain.set_flow_algorithm('1_5')
        L = domain.quantities['stage'].vertex_values

        def stage(x, y):
            return x**2

        domain.set_quantity('stage', stage, location='centroids')

        domain.quantities['stage'].compute_gradients()

        a, b = domain.quantities['stage'].get_gradients()

        assert num.allclose(a[1], 3.33333334)
        assert num.allclose(b[1], 0.0)

        domain._order_ = 1
        domain.distribute_to_vertices_and_edges()

        assert num.allclose(L[1], [-0.22222222,  2.77777778,  2.77777778])

        domain._order_ = 2
        domain.beta_w = 0.9
        domain.beta_w_dry = 0.9
        domain.beta_uh = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh = 0.9
        domain.beta_vh_dry = 0.9
        domain.distribute_to_vertices_and_edges()

        assert num.allclose(L[1], [-2.66666667,  4.        ,  4.        ])

    def test_distribute_away_from_bed1(self):
        #Using test data generated by abstract_2d_finite_volumes-2
        #Assuming no friction and flat bed (0.0)

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
        L = domain.quantities['stage'].vertex_values

        def stage(x, y):
            return x**4 + y**2

        domain.set_quantity('stage', stage, location='centroids')

        domain.quantities['stage'].compute_gradients()
        a, b = domain.quantities['stage'].get_gradients()
        assert num.allclose(a[1], 25.18518519)
        assert num.allclose(b[1], 3.33333333)

        domain._order_ = 1
        domain.distribute_to_vertices_and_edges()

        assert num.allclose(L[1], 4.9382716)

        domain._order_ = 2
        domain.beta_w = 0.9
        domain.beta_w_dry = 0.9
        domain.beta_uh = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh = 0.9
        domain.beta_vh_dry = 0.9
        domain.distribute_to_vertices_and_edges()

        assert num.allclose(L[1], [ -7.81670675,   9.95991663,  12.67160494])

    def test_distribute_near_bed(self):
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

        # Set up for a gradient of (10,0) at mid triangle (bce)
        def slope(x, y):
            return 10*x

        h = 0.1
        def stage(x, y):
            return slope(x, y) + h

        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage, location='centroids')

        E = domain.quantities['elevation'].vertex_values
        L = domain.quantities['stage'].vertex_values

        # Get reference values
        volumes = []
        for i in range(len(L)):
            volumes.append(num.sum(L[i]) / 3)
            assert num.allclose(volumes[i],
                                domain.quantities['stage'].centroid_values[i])

        domain._order_ = 1

        domain.tight_slope_limiters = 0
        domain.distribute_to_vertices_and_edges()

        assert num.allclose(L[1], [0.1, 20.1, 20.1])
        for i in range(len(L)):
            assert num.allclose(volumes[i], num.sum(L[i]) / 3)

        # Allow triangle to be flatter (closer to bed)
        domain.tight_slope_limiters = 1

        domain.distribute_to_vertices_and_edges()


        #print L[1]
        #print [0.298, 20.001, 20.001]

        assert num.allclose(L[1], [  0.1,  20.1,  20.1])




        for i in range(len(L)):
            assert num.allclose(volumes[i], num.sum(L[i]) / 3)

        domain._order_ = 2

        domain.tight_slope_limiters = 0
        domain.distribute_to_vertices_and_edges()
        assert num.allclose(L[1], [0.1, 20.1, 20.1])
        for i in range(len(L)):
            assert num.allclose(volumes[i], num.sum(L[i]) / 3)

        # Allow triangle to be flatter (closer to bed)
        domain.tight_slope_limiters = 1

        domain.distribute_to_vertices_and_edges()
        assert num.allclose(L[1], [0.1, 20.1, 20.1])

        for i in range(len(L)):
            assert num.allclose(volumes[i], num.sum(L[i]) / 3)

    def test_distribute_near_bed1(self):
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
        domain.set_flow_algorithm('1_5')

        # Set up for a gradient of (8,2) at mid triangle (bce)
        def slope(x, y):
            return x**4 + y**2

        h = 0.1
        def stage(x, y):
            return slope(x, y) + h

        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage)

        E = domain.quantities['elevation'].vertex_values
        L = domain.quantities['stage'].vertex_values

        # Get reference values
        volumes = []
        for i in range(len(L)):
            volumes.append(num.sum(L[i]) / 3)
            assert num.allclose(volumes[i],
                                domain.quantities['stage'].centroid_values[i])

        domain._order_ = 1

        domain.tight_slope_limiters = 0
        domain.distribute_to_vertices_and_edges()
        assert num.allclose(L[1], [4.1, 16.1, 20.1])
        for i in range(len(L)):
            assert num.allclose(volumes[i], num.sum(L[i]) / 3)


        domain.tight_slope_limiters = 1 # Allow triangle to be flatter (closer to bed)
        domain.distribute_to_vertices_and_edges()
        #print L[1]
        assert num.allclose(L[1], [  4.239986, 16.060004,  20.00001 ], rtol=1.0e-2)
        for i in range(len(L)):
            assert num.allclose(volumes[i], num.sum(L[i]) / 3)


        domain._order_ = 2

        domain.tight_slope_limiters = 0
        domain.distribute_to_vertices_and_edges()
        assert num.allclose(L[1], [4.1, 16.1, 20.1])
        for i in range(len(L)):
            assert num.allclose(volumes[i], num.sum(L[i]) / 3)

        # Allow triangle to be flatter (closer to bed)
        domain.tight_slope_limiters = 1

        domain.distribute_to_vertices_and_edges()
        assert  num.allclose(L[1], [4.19444976, 16.10554024, 20.00001]) or\
                num.allclose(L[1], [4.23370103, 16.06529897, 20.001]) or\
                num.allclose(L[1], [4.18944138, 16.10955862, 20.001]) or\
                num.allclose(L[1], [4.19351461, 16.10548539, 20.001]) or\
                num.allclose(L[1], [4.00001   , 16.15431247, 20.14567753])

        for i in range(len(L)):
            assert num.allclose(volumes[i], num.sum(L[i]) / 3)


    def test_second_order_distribute_real_data(self):
        #Using test data generated by abstract_2d_finite_volumes-2
        #Assuming no friction and flat bed (0.0)

        a = [0.0, 0.0]
        b = [0.0, 1.0/5]
        c = [0.0, 2.0/5]
        d = [1.0/5, 0.0]
        e = [1.0/5, 1.0/5]
        f = [1.0/5, 2.0/5]
        g = [2.0/5, 2.0/5]

        points = [a, b, c, d, e, f, g]
        #             bae,     efb,     cbf,     feg
        vertices = [[1,0,4], [4,5,1], [2,1,5], [5,4,6]]

        domain = Domain(points, vertices)

        def slope(x, y):
            return -x / 3

        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage',
                            [0.01298164, 0.00365611,
                             0.01440365, -0.0381856437096],
                            location='centroids')
        domain.set_quantity('xmomentum',
                            [0.00670439, 0.01263789,
                             0.00647805, 0.0178180740668],
                            location='centroids')
        domain.set_quantity('ymomentum',
                            [-7.23510980e-004, -6.30413883e-005,
                             6.30413883e-005, 0.000200907255866],
                            location='centroids')

        E = domain.quantities['elevation'].vertex_values
        L = domain.quantities['stage'].vertex_values
        X = domain.quantities['xmomentum'].vertex_values
        Y = domain.quantities['ymomentum'].vertex_values

        domain._order_ = 2
        domain.beta_w = 0.9
        domain.beta_w_dry = 0.9
        domain.beta_uh = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh = 0.9
        domain.beta_vh_dry = 0.9

        # FIXME (Ole): Need tests where this is commented out
        domain.tight_slope_limiters = 0       # Backwards compatibility (14/4/7)
        domain.use_centroid_velocities = 0    # Backwards compatibility (7/5/8)

        domain.distribute_to_vertices_and_edges()


        L_EX = [-0.00825736, -0.00801881,  0.0272445 ]
        X_EX = [ 0.0170432 ,  0.01674667,  0.00654035]
        Y_EX = [ -1.56119371e-04,   1.10107449e-04,  -5.74034758e-05]

        #pprint(L[1,:])
        #pprint(X[1,:])
        #pprint(Y[1,:])

        L_EX1 = [-0.01434766, -0.01292565,  0.03824164]
        X_EX1 = [ 0.01631331,  0.01583494,  0.00596076]
        Y_EX1 = [-0.00036261,  0.00068027, -0.00048923]

        assert num.allclose(L[1,:], L_EX) or \
            num.allclose(L[1,:], L_EX1)

        assert num.allclose(X[1,:], X_EX) or \
            num.allclose(X[1,:], X_EX1)

        assert num.allclose(Y[1,:], Y_EX) or \
            num.allclose(Y[1,:], Y_EX1)


    def test_extrapolate_second_order_sw_real_data(self):
        #Using test data generated by abstract_2d_finite_volumes-2
        #Assuming no friction and flat bed (0.0)

        a = [0.0, 0.0]
        b = [0.0, 1.0/5]
        c = [0.0, 2.0/5]
        d = [1.0/5, 0.0]
        e = [1.0/5, 1.0/5]
        f = [1.0/5, 2.0/5]
        g = [2.0/5, 2.0/5]

        points = [a, b, c, d, e, f, g]
        #             bae,     efb,     cbf,     feg
        vertices = [[1,0,4], [4,5,1], [2,1,5], [5,4,6]]

        domain = Domain(points, vertices)

        def slope(x, y):
            return -x / 3

        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage',
                            [0.01298164, 0.00365611,
                             0.01440365, -0.0381856437096],
                            location='centroids')
        domain.set_quantity('xmomentum',
                            [0.00670439, 0.01263789,
                             0.00647805, 0.0178180740668],
                            location='centroids')
        domain.set_quantity('ymomentum',
                            [-7.23510980e-004, -6.30413883e-005,
                             6.30413883e-005, 0.000200907255866],
                            location='centroids')

        E = domain.quantities['elevation'].vertex_values
        L = domain.quantities['stage'].vertex_values
        X = domain.quantities['xmomentum'].vertex_values
        Y = domain.quantities['ymomentum'].vertex_values

        domain._order_ = 2
        domain.beta_w = 0.9
        domain.beta_w_dry = 0.9
        domain.beta_uh = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh = 0.9
        domain.beta_vh_dry = 0.9

        # FIXME (Ole): Need tests where this is commented out
        domain.tight_slope_limiters = 0       # Backwards compatibility (14/4/7)
        domain.use_centroid_velocities = 0    # Backwards compatibility (7/5/8)
        domain.optimised_gradient_limiter = True

        #domain.extrapolate_second_order_sw()


        from  anuga.shallow_water.shallow_water_ext import extrapolate_second_order_sw

        extrapolate_second_order_sw(domain)

        #domain.distribute_to_vertices_and_edges()


        #print L[1,:]
        #print X[1,:]
        #print Y[1,:]

        L_EX = [-0.00137913, -0.00098143,  0.0133289 ]
        X_EX = [ 0.01998251,  0.01948814,  0.00247195]
        Y_EX = [ -2.18223005e-04,   2.25636080e-04,  -5.36417393e-05]


        #L_EX = [-0.00825736, -0.00801881,  0.0272445 ]
        #X_EX = [ 0.0170432 ,  0.01674667,  0.00654035]
        #Y_EX = [ -1.56119371e-04,   1.10107449e-04,  -5.74034758e-05]

        assert num.allclose(L[1,:], L_EX)
        assert num.allclose(X[1,:], X_EX)
        assert num.allclose(Y[1,:], Y_EX)

    def test_balance_deep_and_shallow(self):
        """Test that balanced limiters preserve conserved quantites.
        This test is using old depth based balanced limiters
        """

        import copy

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #             bac,     bce,     ecf,     dbe
        elements = [[1,0,2], [1,2,4], [4,2,5], [3,1,4] ]

        domain = Domain(points, elements)
        domain.check_integrity()

        # Create a deliberate overshoot
        domain.set_quantity('stage', [[3,0,3], [2,2,6], [5,3,8], [8,3,5]])
        domain.set_quantity('elevation', 0)    # Flat bed
        stage = domain.quantities['stage']

        ref_centroid_values = copy.copy(stage.centroid_values[:])    # Copy

        # Limit
        domain.tight_slope_limiters = 0
        domain.distribute_to_vertices_and_edges()

        # Assert that quantities are conserved
        for k in range(len(domain)):
            assert num.allclose(ref_centroid_values[k],
                                num.sum(stage.vertex_values[k,:] / 3))

        # Now try with a non-flat bed - closely hugging initial stage in places
        # This will create alphas in the range [0, 0.478260, 1]
        domain.set_quantity('stage', [[3,0,3], [2,2,6], [5,3,8], [8,3,5]])
        domain.set_quantity('elevation', [[0,0,0],
                                          [1.8,1.9,5.9],
                                          [4.6,0,0],
                                          [0,2,4]])
        stage = domain.quantities['stage']

        ref_centroid_values = copy.copy(stage.centroid_values[:])    # Copy
        ref_vertex_values = copy.copy(stage.vertex_values[:])        # Copy

        # Limit
        domain.tight_slope_limiters = 0
        domain.distribute_to_vertices_and_edges()

        # Assert that all vertex quantities have changed
        for k in range(len(domain)):
            assert not num.allclose(ref_vertex_values[k,:],
                                    stage.vertex_values[k,:])
        # and assert that quantities are still conserved
        for k in range(len(domain)):
            assert num.allclose(ref_centroid_values[k],
                                num.sum(stage.vertex_values[k,:] / 3))


        assert num.allclose(stage.vertex_values,
                            [[ 2.        ,  2.        ,  2.        ],
                             [ 3.33333333,  3.33333333,  3.33333333],
                             [ 5.33333333,  5.33333333,  5.33333333],
                             [ 5.33333333,  5.33333333,  5.33333333]])

    def test_balance_deep_and_shallow_tight_SL(self):
        """Test that balanced limiters preserve conserved quantites.
        This test is using Tight Slope Limiters
        """

        import copy

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #             bac,     bce,     ecf,     dbe
        elements = [[1,0,2], [1,2,4], [4,2,5], [3,1,4] ]

        domain = Domain(points, elements)
        domain.check_integrity()

        # Create a deliberate overshoot
        domain.set_quantity('stage', [[3,0,3], [2,2,6], [5,3,8], [8,3,5]])
        domain.set_quantity('elevation', 0)    # Flat bed
        stage = domain.quantities['stage']

        ref_centroid_values = copy.copy(stage.centroid_values[:])    # Copy

        # Limit
        domain.tight_slope_limiters = 1
        domain.distribute_to_vertices_and_edges()

        # Assert that quantities are conserved
        for k in range(len(domain)):
            assert num.allclose (ref_centroid_values[k],
                                 num.sum(stage.vertex_values[k,:] / 3))

        # Now try with a non-flat bed - closely hugging initial stage in places
        # This will create alphas in the range [0, 0.478260, 1]
        domain.set_quantity('stage', [[3,0,3], [2,2,6], [5,3,8], [8,3,5]])
        domain.set_quantity('elevation', [[0,0,0],
                                          [1.8,1.9,5.9],
                                          [4.6,0,0],
                                          [0,2,4]])
        stage = domain.quantities['stage']

        ref_centroid_values = copy.copy(stage.centroid_values[:])    # Copy
        ref_vertex_values = copy.copy(stage.vertex_values[:])        # Copy

        # Limit
        domain.tight_slope_limiters = 1
        domain.distribute_to_vertices_and_edges()

        # Assert that all vertex quantities have changed
        for k in range(len(domain)):
            assert not num.allclose(ref_vertex_values[k,:],
                                    stage.vertex_values[k,:])
        # and assert that quantities are still conserved
        for k in range(len(domain)):
            assert num.allclose(ref_centroid_values[k],
                                num.sum(stage.vertex_values[k,:] / 3))

    def test_balance_deep_and_shallow_Froude(self):
        """Test that balanced limiters preserve conserved quantites -
        and also that excessive Froude numbers are dealt with.
        This test is using tight slope limiters.
        """

        import copy

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #             bac,     bce,     ecf,     dbe
        elements = [[1,0,2], [1,2,4], [4,2,5], [3,1,4] ]

        domain = Domain(points, elements)
        domain.check_integrity()
        domain.tight_slope_limiters = True
        domain.use_centroid_velocities = True

        # Create non-flat bed - closely hugging initial stage in places
        # This will create alphas in the range [0, 0.478260, 1]
        domain.set_quantity('stage', [[3,0,3], [2,2,6], [5,3,8], [8,3,5]])
        domain.set_quantity('elevation', [[0,0,0],
                                          [1.8,1.999,5.999],
                                          [4.6,0,0],
                                          [0,2,4]])

        # Create small momenta, that nonetheless will generate large speeds
        # due to shallow depth at isolated vertices
        domain.set_quantity('xmomentum', -0.0058)
        domain.set_quantity('ymomentum', 0.0890)

        stage = domain.quantities['stage']
        elevation = domain.quantities['elevation']
        xmomentum = domain.quantities['xmomentum']
        ymomentum = domain.quantities['ymomentum']

        # Setup triangle #1 to mimick real Froude explosion observed
        # in the Onslow example 13 Nov 2007.
        stage.vertex_values[1,:] = [1.6385, 1.6361, 1.2953]
        elevation.vertex_values[1,:] = [1.6375, 1.6336, 0.4647]
        xmomentum.vertex_values[1,:] = [-0.0058, -0.0050, -0.0066]
        ymomentum.vertex_values[1,:] = [0.0890, 0.0890, 0.0890]

        xmomentum.interpolate()
        ymomentum.interpolate()
        stage.interpolate()
        elevation.interpolate()

        # Verify interpolation
        assert num.allclose(stage.centroid_values[1], 1.5233)
        assert num.allclose(elevation.centroid_values[1], 1.2452667)
        assert num.allclose(xmomentum.centroid_values[1], -0.0058)
        assert num.allclose(ymomentum.centroid_values[1], 0.089)

        # Derived quantities
        depth = stage - elevation
        u = xmomentum / depth
        v = ymomentum / depth
        assert num.allclose(u.centroid_values, [-0.0029, -0.02086081, -0.00152632, -0.00174])
        assert num.allclose(v.centroid_values, [0.0445, 0.3201055, 0.02342105, 0.0267])

        denom = (depth*g)**0.5
        Fx = u / denom
        Fy = v / denom

        # Verify against Onslow example (14 Nov 2007)
        assert num.allclose(depth.centroid_values[1], 0.278033)
        assert num.allclose(u.centroid_values[1], -0.0208608)
        assert num.allclose(v.centroid_values[1], 0.3201055)

        assert num.allclose(denom.centroid_values[1],
                            num.sqrt(depth.centroid_values[1]*g))

        assert num.allclose(u.centroid_values[1] / denom.centroid_values[1],
                            -0.012637746977)
        assert num.allclose(Fx.centroid_values[1],
                            u.centroid_values[1] / denom.centroid_values[1])

        # Check that Froude numbers are small at centroids.
        assert num.allclose(Fx.centroid_values[1], -0.012637746977)
        assert num.allclose(Fy.centroid_values[1], 0.193924048435)

        # But Froude numbers are huge at some vertices and edges
        assert num.allclose(Fx.vertex_values[1,:], [-5.85888475e+01,
                                                    -1.27775313e+01,
                                                    -2.78511420e-03])

        assert num.allclose(Fx.edge_values[1,:], [-6.89150773e-03,
                                                  -7.38672488e-03,
                                                  -2.35626238e+01])

        assert num.allclose(Fy.vertex_values[1,:], [8.99035764e+02,
                                                    2.27440057e+02,
                                                    3.75568430e-02])

        assert num.allclose(Fy.edge_values[1,:], [1.05748998e-01,
                                                  1.06035244e-01,
                                                  3.88346947e+02])


        # The task is now to arrange the limiters such that Froude numbers
        # remain under control whil at the same time obeying the conservation
        # laws.
        ref_centroid_values = copy.copy(stage.centroid_values[:])    # Copy
        ref_vertex_values = copy.copy(stage.vertex_values[:])        # Copy

        # Limit (and invoke balance_deep_and_shallow)
        domain.tight_slope_limiters = 1
        domain.distribute_to_vertices_and_edges()

        # Redo derived quantities
        depth = stage - elevation
        u = xmomentum / depth
        v = ymomentum / depth

        # Assert that all vertex velocities stay within one
        # order of magnitude of centroid velocities.
        assert num.alltrue(num.absolute(u.vertex_values[1,:]) <=
                           num.absolute(u.centroid_values[1])*10)
        assert num.alltrue(num.absolute(v.vertex_values[1,:]) <=
                           num.absolute(v.centroid_values[1])*10)

        denom = (depth*g)**0.5
        Fx = u / denom
        Fy = v / denom

        # Assert that Froude numbers are less than max value (TBA)
        # at vertices, edges and centroids.
        from anuga.config import maximum_froude_number

        assert num.alltrue(num.absolute(Fx.vertex_values[1,:]) <
                           maximum_froude_number)
        assert num.alltrue(num.absolute(Fy.vertex_values[1,:]) <
                           maximum_froude_number)

        # Assert that all vertex quantities have changed
        for k in range(len(domain)):
            assert not num.allclose(ref_vertex_values[k,:],
                                    stage.vertex_values[k,:])

        # Assert that quantities are still conserved
        for k in range(len(domain)):
            assert num.allclose(ref_centroid_values[k],
                                num.sum(stage.vertex_values[k,:]) / 3)

        return

        qwidth = 12
        for k in [1]:    # range(len(domain)):
            print('Triangle %d (C, V, E)' % k)

            print(('stage'.ljust(qwidth), stage.centroid_values[k],
                   stage.vertex_values[k,:], stage.edge_values[k,:]))
            print(('elevation'.ljust(qwidth), elevation.centroid_values[k],
                   elevation.vertex_values[k,:], elevation.edge_values[k,:]))
            print(('depth'.ljust(qwidth), depth.centroid_values[k],
                   depth.vertex_values[k,:], depth.edge_values[k,:]))
            print(('xmomentum'.ljust(qwidth), xmomentum.centroid_values[k],
                   xmomentum.vertex_values[k,:], xmomentum.edge_values[k,:]))
            print(('ymomentum'.ljust(qwidth), ymomentum.centroid_values[k],
                   ymomentum.vertex_values[k,:], ymomentum.edge_values[k,:]))
            print(('u'.ljust(qwidth),u.centroid_values[k],
                   u.vertex_values[k,:], u.edge_values[k,:]))
            print(('v'.ljust(qwidth), v.centroid_values[k],
                   v.vertex_values[k,:], v.edge_values[k,:]))
            print(('Fx'.ljust(qwidth), Fx.centroid_values[k],
                   Fx.vertex_values[k,:], Fx.edge_values[k,:]))
            print(('Fy'.ljust(qwidth), Fy.centroid_values[k],
                   Fy.vertex_values[k,:], Fy.edge_values[k,:]))

    def test_evolve_finaltime(self):
        """Test evolve with finaltime set
        """

        from anuga import rectangular_cross_domain

        # Create basic mesh
        domain = rectangular_cross_domain(6, 6)
        domain.set_name('evolve_finaltime')

        # IC
        def x_slope(x, y):
            return x / 3

        domain.set_quantity('elevation', 0)
        domain.set_quantity('friction', 0)
        domain.set_quantity('stage', x_slope)

        # Boundary conditions (reflective everywhere)
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        domain.check_integrity()

        # Evolution
        # Test that t is a float
        tt = 0.0
        for t in domain.evolve(yieldstep=0.05, finaltime=5.0):
            tt += t

        assert num.allclose(tt,252.5)


        os.remove(domain.get_name() + '.sww')

    def test_evolve_duration(self):
        """Test evolve with duration set
        """

        from anuga import rectangular_cross_domain

        # Create basic mesh
        domain = rectangular_cross_domain(6, 6)
        domain.set_name('evolve_duration')

        # IC
        def x_slope(x, y):
            return x / 3

        domain.set_quantity('elevation', 0)
        domain.set_quantity('friction', 0)
        domain.set_quantity('stage', x_slope)

        # Boundary conditions (reflective everywhere)
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        domain.check_integrity()

        # Evolution
        # Test that t is a float
        tt = 0.0
        for t in domain.evolve(yieldstep=0.05, duration=5.0):
            tt += t

        assert num.allclose(tt,252.5)

        os.remove(domain.get_name() + '.sww')

    def test_evolve_and_set_time(self):
        """Test evolve with set_time before evolve
        """

        from anuga import rectangular_cross_domain

        # Create basic mesh
        domain = rectangular_cross_domain(6, 6)
        domain.set_name('evolve_and_set_time')

        # IC
        def x_slope(x, y):
            return x / 3

        domain.set_quantity('elevation', 0)
        domain.set_quantity('friction', 0)
        domain.set_quantity('stage', x_slope)

        # Boundary conditions (reflective everywhere)
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        domain.check_integrity()

        
        domain.set_time(0.5)


        # Evolution
        # Test that t is a float
        tt = 0.0
        for t in domain.evolve(yieldstep=0.05, outputstep=1.0, duration=5.0):
            tt += t

        assert num.allclose(tt,252.5) 

        os.remove(domain.get_name() + '.sww')       


    def test_evolve_outputstep(self):
        """Test evolve with outputstep set
        """

        from anuga import rectangular_cross_domain

        # Create basic mesh
        domain = rectangular_cross_domain(6, 6)
        domain.set_name('evolve_outputstep')

        # IC
        def x_slope(x, y):
            return x / 3

        domain.set_quantity('elevation', 0)
        domain.set_quantity('friction', 0)
        domain.set_quantity('stage', x_slope)

        # Boundary conditions (reflective everywhere)
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        domain.check_integrity()

        # Evolution
        # Test that t is a float
        tt = 0.0
        for t in domain.evolve(yieldstep=0.05, outputstep=1.0, duration=5.0):
            tt += t

        assert num.allclose(tt,252.5)

        # Open sww file to check that only store every second

        # Read results for specific timesteps t=1 and t=2
        fid = NetCDFFile(domain.get_name() + '.sww')
        time = fid.variables['time'][:]
        stage = fid.variables['stage'][:,:]
        fid.close()

        os.remove(domain.get_name() + '.sww')

        timeslices = 6
        msg = f' time.shape[0] = {time.shape[0]}, expected {timeslices}'
        assert time.shape[0] == timeslices, msg

        msg = f' time.shape[0] = {stage.shape[0]}, expected {timeslices}'
        assert stage.shape[0] == timeslices

    def test_evolve_outputstep_integer(self):
        """Test evolve outputstep when it is an integer multipe of yieldstep
        """

        from anuga import rectangular_cross_domain

        # Create basic mesh
        domain = rectangular_cross_domain(6, 6)
        domain.set_name('evolve_outputstep_integer')

        # IC
        def x_slope(x, y):
            return x / 3

        domain.set_quantity('elevation', 0)
        domain.set_quantity('friction', 0)
        domain.set_quantity('stage', x_slope)

        # Boundary conditions (reflective everywhere)
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        domain.check_integrity()

        # Evolution
        # Test that t is a float
        tt = 0.0
        for t in domain.evolve(yieldstep=0.05, outputstep=0.05, duration=5.0):
            tt += t

        assert num.allclose(tt,252.5)

        # Open sww file to check that only store every second

        # Read results for specific timesteps t=1 and t=2
        fid = NetCDFFile(domain.get_name() + '.sww')
        time = fid.variables['time'][:]
        stage = fid.variables['stage'][:,:]
        fid.close()

        os.remove(domain.get_name() + '.sww')

        timeslices = 101
        msg = f' time.shape[0] = {time.shape[0]}, expected {timeslices}'
        assert time.shape[0] == timeslices, msg

        msg = f' time.shape[0] = {stage.shape[0]}, expected {timeslices}'
        assert stage.shape[0] == timeslices

    def test_evolve_outputstep_non_integer(self):
        """Test exception if evolve outputstep is not integer multiple of yieldstep
        """

        from anuga import rectangular_cross_domain

        # Create basic mesh
        domain = rectangular_cross_domain(6, 6)
        domain.set_name('evolve_outputstep_non_integer')

        # IC
        def x_slope(x, y):
            return x / 3

        domain.set_quantity('elevation', 0)
        domain.set_quantity('friction', 0)
        domain.set_quantity('stage', x_slope)

        # Boundary conditions (reflective everywhere)
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        domain.check_integrity()

        # Evolution
        # Test that t is a float
        tt = 0.0

        try:
            for t in domain.evolve(yieldstep=0.05, outputstep=0.12, duration=5.0):
                tt += t
        except AssertionError:
            # Getting here is good as outputstep is not an integer multiple of yieldstep
            return

        # Shouldn't get here
        raise Exception('An AssertionError should have occurred earlier')


    def test_conservation_1(self):
        """Test that stage is conserved globally

        This one uses a flat bed, reflective bdries and a suitable
        initial condition
        """

        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 2

        # IC
        def x_slope(x, y):
            return x / 3

        domain.set_quantity('elevation', 0)
        domain.set_quantity('friction', 0)
        domain.set_quantity('stage', x_slope)

        # Boundary conditions (reflective everywhere)
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        domain.check_integrity()

        initial_volume = domain.quantities['stage'].get_integral()
        initial_xmom = domain.quantities['xmomentum'].get_integral()

        # Evolution
        for t in domain.evolve(yieldstep=0.05, finaltime=5.0):
            volume = domain.quantities['stage'].get_integral()
            assert num.allclose(volume, initial_volume)

            #I don't believe that the total momentum should be the same
            #It starts with zero and ends with zero though
            #xmom = domain.quantities['xmomentum'].get_integral()
            #print xmom
            #assert allclose (xmom, initial_xmom)

        os.remove(domain.get_name() + '.sww')

    def test_conservation_2(self):
        """Test that stage is conserved globally

        This one uses a slopy bed, reflective bdries and a suitable
        initial condition
        """

        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 2

        # IC
        def x_slope(x, y):
            return x / 3

        domain.set_quantity('elevation', x_slope)
        domain.set_quantity('friction', 0)
        domain.set_quantity('stage', 0.4)    # Steady

        # Boundary conditions (reflective everywhere)
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        domain.check_integrity()

        initial_volume = domain.quantities['stage'].get_integral()
        initial_xmom = domain.quantities['xmomentum'].get_integral()

        # Evolution
        for t in domain.evolve(yieldstep=0.05, finaltime=5.0):
            volume = domain.quantities['stage'].get_integral()
            assert num.allclose(volume, initial_volume)

            #FIXME: What would we expect from momentum
            #xmom = domain.quantities['xmomentum'].get_integral()
            #print xmom
            #assert allclose (xmom, initial_xmom)

        os.remove(domain.get_name() + '.sww')

    def test_conservation_3(self):
        """Test that stage is conserved globally

        This one uses a larger grid, convoluted bed, reflective boundaries
        and a suitable initial condition
        """

        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(2, 1)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 2


        # IC
        def x_slope(x, y):
            z = 0*x
            for i in range(len(x)):
                if x[i] < 0.3:
                    z[i] = x[i] / 3
                if 0.3 <= x[i] < 0.5:
                    z[i] = -0.5
                if 0.5 <= x[i] < 0.7:
                    z[i] = 0.39
                if 0.7 <= x[i]:
                    z[i] = x[i] / 3
            return z

        domain.set_quantity('elevation', x_slope)
        domain.set_quantity('friction', 0)
        domain.set_quantity('stage', 0.4) #Steady

        # Boundary conditions (reflective everywhere)
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        domain.check_integrity()

        initial_volume = domain.quantities['stage'].get_integral()
        initial_xmom = domain.quantities['xmomentum'].get_integral()

        import copy

        ref_centroid_values = copy.copy(domain.quantities['stage'].\
                                            centroid_values)

        domain.distribute_to_vertices_and_edges()

        assert num.allclose(domain.quantities['stage'].centroid_values,
                            ref_centroid_values)

        # Check that initial limiter doesn't violate cons quan
        assert num.allclose(domain.quantities['stage'].get_integral(),
                            initial_volume)

        # Evolution
        for t in domain.evolve(yieldstep=0.05, finaltime=10):
            volume =  domain.quantities['stage'].get_integral()
            assert num.allclose (volume, initial_volume)

        os.remove(domain.get_name() + '.sww')

    def test_conservation_4(self):
        """Test that stage is conserved globally

        This one uses a larger grid, convoluted bed, reflective boundaries
        and a suitable initial condition
        """

        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 2

        # IC
        def x_slope(x, y):
            z = 0*x
            for i in range(len(x)):
                if x[i] < 0.3:
                    z[i] = x[i] / 3
                if 0.3 <= x[i] < 0.5:
                    z[i] = -0.5
                if 0.5 <= x[i] < 0.7:
                    #z[i] = 0.3     # OK with beta == 0.2
                    z[i] = 0.34     # OK with beta == 0.0
                    #z[i] = 0.35    # Fails after 80 timesteps with an error
                                    # of the order 1.0e-5
                if 0.7 <= x[i]:
                    z[i] = x[i] / 3
            return z

        domain.set_quantity('elevation', x_slope)
        domain.set_quantity('friction', 0)
        domain.set_quantity('stage', 0.4) #Steady

        # Boundary conditions (reflective everywhere)
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        domain.check_integrity()

        initial_volume = domain.quantities['stage'].get_integral()
        initial_xmom = domain.quantities['xmomentum'].get_integral()

        import copy

        ref_centroid_values = copy.copy(domain.quantities['stage'].\
                                            centroid_values)

        # Test limiter by itself
        domain.distribute_to_vertices_and_edges()

        # Check that initial limiter doesn't violate cons quan
        assert num.allclose(domain.quantities['stage'].get_integral(),
                            initial_volume)
        # NOTE: This would fail if any initial stage was less than the
        # corresponding bed elevation - but that is reasonable.

        #Evolution
        for t in domain.evolve(yieldstep=0.05, finaltime=10.0):
            volume =  domain.quantities['stage'].get_integral()
            assert num.allclose (volume, initial_volume)

        os.remove(domain.get_name() + '.sww')

    def test_conservation_5(self):
        """Test that momentum is conserved globally in steady state scenario

        This one uses a slopy bed, dirichlet and reflective bdries
        """

        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.set_low_froude(0)
        domain.smooth = False
        domain.default_order = 2

        # IC
        def x_slope(x, y):
            return x / 3

        domain.set_quantity('elevation', x_slope)
        domain.set_quantity('friction', 0)
        domain.set_quantity('stage', 0.4) # Steady

        # Boundary conditions (reflective everywhere)
        Br = Reflective_boundary(domain)
        Bleft = Dirichlet_boundary([0.5, 0, 0])
        Bright = Dirichlet_boundary([0.1, 0, 0])
        domain.set_boundary({'left': Bleft, 'right': Bright,
                             'top': Br, 'bottom': Br})

        domain.check_integrity()

        initial_volume = domain.quantities['stage'].get_integral()
        initial_xmom = domain.quantities['xmomentum'].get_integral()

        # Evolution
        for t in domain.evolve(yieldstep=0.05, finaltime=15.0):
            stage =  domain.quantities['stage'].get_integral()
            xmom = domain.quantities['xmomentum'].get_integral()
            ymom = domain.quantities['ymomentum'].get_integral()
            #print stage,xmom,ymom

            if num.allclose(t, 10):    # Steady state reached
                steady_xmom = domain.quantities['xmomentum'].get_integral()
                steady_ymom = domain.quantities['ymomentum'].get_integral()
                steady_stage = domain.quantities['stage'].get_integral()

            if t > 10:
                msg = 'xmom=%.2f, steady_xmom=%.2f' % (xmom, steady_xmom)
                assert num.allclose(xmom, steady_xmom), msg
                assert num.allclose(ymom, steady_ymom)
                assert num.allclose(stage, steady_stage)

        os.remove(domain.get_name() + '.sww')

    def test_conservation_real(self):
        """Test that momentum is conserved globally

        Stephen finally made a test that revealed the problem.
        This test failed with code prior to 25 July 2005
        """

        import sys
        import os.path
        sys.path.append(os.path.join('..', 'abstract_2d_finite_volumes'))
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        yieldstep = 0.01
        finaltime = 0.05
        min_depth = 1.0e-2

        #Create shallow water domain
        points, vertices, boundary = rectangular(10, 10, len1=500, len2=500)
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 1
        domain.minimum_allowed_height = min_depth

        # Set initial condition
        class Set_IC(object):
            """Set an initial condition with a constant value, for x0<x<x1"""

            def __init__(self, x0=0.25, x1=0.5, h=1.0):
                self.x0 = x0
                self.x1 = x1
                self.h  = h

            def __call__(self, x, y):
                return self.h*((x > self.x0) & (x < self.x1))

        domain.set_quantity('stage', Set_IC(200.0, 300.0, 5.0))

        # Boundaries
        R = Reflective_boundary(domain)
        domain.set_boundary({'left': R, 'right': R, 'top':R, 'bottom': R})

        ref = domain.quantities['stage'].get_integral()

        # Evolution
        for t in domain.evolve(yieldstep=yieldstep, finaltime=finaltime):
            pass

        now = domain.quantities['stage'].get_integral()

        msg = 'Stage not conserved: was %f, now %f' % (ref, now)
        assert num.allclose(ref, now), msg

        os.remove(domain.get_name() + '.sww')

    def test_second_order_flat_bed_onestep(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        #Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.set_flow_algorithm('1_5')
        domain.smooth = False
        domain.default_order = 2
        domain.beta_w = 0.9
        domain.beta_w_dry = 0.9
        domain.beta_uh = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh = 0.9
        domain.beta_vh_dry = 0.9
        domain.H0 = 1.0e-3

        # Boundary conditions
        Br = Reflective_boundary(domain)
        Bd = Dirichlet_boundary([0.1, 0., 0.])
        domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})

        domain.check_integrity()

        # Evolution
        for t in domain.evolve(yieldstep=0.05, finaltime=0.05):
            pass

        # Data from earlier version of abstract_2d_finite_volumes
        assert num.allclose(domain.recorded_min_timestep, 0.0396825396825)
        assert num.allclose(domain.recorded_max_timestep, 0.0396825396825)

        msg = ("domain.quantities['stage'].centroid_values[:12]=%s"
               % str(domain.quantities['stage'].centroid_values[:12]))
        assert num.allclose(domain.quantities['stage'].centroid_values[:12],
                            [0.00171396, 0.02656103, 0.00241523, 0.02656103,
                             0.00241523, 0.02656103, 0.00241523, 0.02656103,
                             0.00241523, 0.02656103, 0.00241523, 0.0272623]), msg

        domain.distribute_to_vertices_and_edges()

        assert num.allclose(domain.quantities['stage'].vertex_values[:12,0],
                            [0.001, 0.02656103, 0.001, 0.02656103, 0.001, 0.02656103,
                             0.001, 0.02656103, 0.001, 0.02656103, 0.001, 0.0272623])

        assert num.allclose(domain.quantities['stage'].vertex_values[:12,1],
                            [0.00237867, 0.02656103, 0.001, 0.02656103, 0.001,
                             0.02656103, 0.001, 0.02656103, 0.001, 0.02656103,
                             0.00110647, 0.0272623])

        assert num.allclose(domain.quantities['stage'].vertex_values[:12,2],
                            [0.00176321, 0.02656103, 0.00524568,
                             0.02656103, 0.00524568, 0.02656103,
                             0.00524568, 0.02656103, 0.00524568,
                             0.02656103, 0.00513921, 0.0272623])

        assert num.allclose(domain.quantities['xmomentum'].centroid_values[:12],
                            [0.00113961, 0.01302432, 0.00148672,
                             0.01302432, 0.00148672, 0.01302432,
                             0.00148672, 0.01302432, 0.00148672 ,
                             0.01302432, 0.00148672, 0.01337143])

        assert num.allclose(domain.quantities['ymomentum'].centroid_values[:12],
                            [-2.91240050e-004 , 1.22721531e-004,
                             -1.22721531e-004,  1.22721531e-004 ,
                             -1.22721531e-004,  1.22721531e-004,
                             -1.22721531e-004 , 1.22721531e-004,
                             -1.22721531e-004,  1.22721531e-004,
                             -1.22721531e-004,  -4.57969873e-005])

        os.remove(domain.get_name() + '.sww')

    def test_second_order_flat_bed_moresteps(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 2

        # Boundary conditions
        Br = Reflective_boundary(domain)
        Bd = Dirichlet_boundary([0.1, 0., 0.])
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
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        N = 8
        points, vertices, boundary = rectangular(N, N)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.set_flow_algorithm('1_5')
        domain.smooth = False
        domain.default_order = 1
        domain.H0 = 1.0e-3 # As suggested in the manual

        # Boundary conditions
        Br = Reflective_boundary(domain)
        Bd = Dirichlet_boundary([0.2, 0., 0.])

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
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        N = 8
        points, vertices, boundary = rectangular(N, N)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.set_flow_algorithm('1_5')
        domain.smooth = False
        domain.default_order = 2
        domain.beta_w = 0.9
        domain.beta_w_dry = 0.9
        domain.beta_uh = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh = 0.9
        domain.beta_vh_dry = 0.9

        domain.H0 = 1.0e-3 # As suggested in the manual
        domain.use_centroid_velocities = False # Backwards compatibility (8/5/8)
        domain.set_maximum_allowed_speed(1.0)

        # Boundary conditions
        Br = Reflective_boundary(domain)
        Bd = Dirichlet_boundary([0.2, 0., 0.])

        domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})
        domain.check_integrity()

        # Evolution
        for t in domain.evolve(yieldstep=0.01, finaltime=0.03):
            pass

        msg = 'min step was %f instead of %f' % (domain.recorded_min_timestep,
                                                 0.0210448446782)

        assert num.allclose(domain.recorded_min_timestep, 0.0210448446782), msg
        assert num.allclose(domain.recorded_max_timestep, 0.0210448446782)

        #FIXME: These numbers were from version before 25/10
        #assert allclose(domain.quantities['stage'].vertex_values[:4,0],
        #                [0.00101913,0.05352143,0.00104852,0.05354394])

        # Slight change due to flux limiter optimisation 18/04/2012

        #pprint(domain.quantities['stage'].vertex_values[:4,0])
        #pprint(domain.quantities['xmomentum'].vertex_values[:4,0])
        #pprint(domain.quantities['ymomentum'].vertex_values[:4,0])



        #W_0 =  [ 0.001     ,  0.05350388,  0.001     ,  0.05352525]
        #UH_0 = [ 0.00044246,  0.03684648,  0.0008209 ,  0.03686007]
        #VH_0 = [-0.00142112,  0.00061559, -0.00062362,  0.00061896]

        W_0 =  [ 0.001     ,  0.05350737,  0.00106727,  0.0535293 ]
        UH_0 = [ 0.00090262,  0.03684904,  0.00090267,  0.03686323]
        VH_0 = [ -1.97310289e-04,   6.10268320e-04,  -6.59631326e-05, 6.14082609e-04]



        assert num.allclose(domain.quantities['stage'].vertex_values[:4,0], W_0)

        assert num.allclose(domain.quantities['xmomentum'].vertex_values[:4,0], UH_0)

        assert num.allclose(domain.quantities['ymomentum'].vertex_values[:4,0], VH_0)


        os.remove(domain.get_name() + '.sww')


    def test_flatbed_second_order_vmax_0(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        N = 8
        points, vertices, boundary = rectangular(N, N)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.set_flow_algorithm('1_5')
        domain.smooth = False
        domain.default_order = 2
        domain.beta_w = 0.9
        domain.beta_w_dry = 0.9
        domain.beta_uh = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh = 0.9
        domain.beta_vh_dry = 0.9
        domain.maximum_allowed_speed = 0.0    # Makes it like the 'oldstyle'
        domain.H0 = 1.0e-3 # As suggested in the manual
        domain.use_centroid_velocities = False # Backwards compatibility (8/5/8)

        # Boundary conditions
        Br = Reflective_boundary(domain)
        Bd = Dirichlet_boundary([0.2, 0., 0.])

        domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})
        domain.check_integrity()

        # Evolution
        for t in domain.evolve(yieldstep=0.01, finaltime=0.03):
            pass

        assert num.allclose(domain.recorded_min_timestep, 0.0210448446782)
        assert num.allclose(domain.recorded_max_timestep, 0.0210448446782)

        #pprint(domain.quantities['stage'].vertex_values[:4,0])
        #pprint(domain.quantities['xmomentum'].vertex_values[:4,0])
        #pprint(domain.quantities['ymomentum'].vertex_values[:4,0])

        UH_EX = [ 0.00090262,  0.03684904,  0.00090267,  0.03686323]
        VH_EX = [ -1.97310289e-04,   6.10268320e-04,  -6.59631326e-05,   6.14082609e-04]

        #UH_EX = [ 0.00044246,  0.03684648,  0.0008209 ,  0.03686007]
        #VH_EX = [-0.00142112,  0.00061559, -0.00062362,  0.00061896]

        assert num.allclose(domain.quantities['xmomentum'].vertex_values[:4,0], UH_EX)
        assert num.allclose(domain.quantities['ymomentum'].vertex_values[:4,0], VH_EX)

        os.remove(domain.get_name() + '.sww')

    def test_flatbed_second_order_distribute(self):
        #Use real data from anuga.abstract_2d_finite_volumes 2
        #painfully setup and extracted.

        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        N = 8
        points, vertices, boundary = rectangular(N, N)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.set_flow_algorithm('1_5')
        domain.smooth = False
        domain.default_order=domain._order_ = 2
        domain.beta_w = 0.9
        domain.beta_w_dry = 0.9
        domain.beta_uh = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh = 0.9
        domain.beta_vh_dry = 0.9

        domain.H0 = 1.0e-3 # ANUGA manual 28/5/9

        # Boundary conditions
        Br = Reflective_boundary(domain)
        Bd = Dirichlet_boundary([0.2, 0., 0.])

        domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})
        domain.check_integrity()

        for V in [False, True]:
            if V:
                # Set centroids as if system had been evolved
                L = num.zeros(2*N*N, float)
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

                X = num.zeros(2*N*N, float)
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

                Y = num.zeros(2*N*N, float)
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
                assert num.allclose(domain.recorded_min_timestep, 0.0210448446782)
                assert num.allclose(domain.recorded_max_timestep, 0.0210448446782)

            #print domain.quantities['stage'].centroid_values[:4]
            #print domain.quantities['xmomentum'].centroid_values[:4]
            #print domain.quantities['ymomentum'].centroid_values[:4]

            #Centroids were correct but not vertices.
            #Hence the check of distribute below.

            if not V:



                W_EX = [ 0.00725574,  0.05350737,  0.01008413,  0.0535293 ]
                UH_EX = [ 0.00654919,  0.03684904,  0.00852886,  0.03686323]
                VH_EX = [-0.00143163,  0.00061027, -0.00062325,  0.00061408]


                assert num.allclose(domain.quantities['stage'].centroid_values[:4], W_EX)
                assert num.allclose(domain.quantities['xmomentum'].centroid_values[:4], UH_EX)
                assert num.allclose(domain.quantities['ymomentum'].centroid_values[:4], VH_EX)


                assert num.allclose(domain.quantities['xmomentum'].centroid_values[17], 0.0, atol=1.0e-3)
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

            assert num.allclose(domain.quantities['xmomentum'].centroid_values[17], 0.0, atol=1.0e-3)


            UH_EX = [ 0.00090262,  0.03684904,  0.00090267,  0.03686323]
            VH_EX = [ -1.97310289e-04,   6.10268320e-04,  -6.59631326e-05,   6.14082609e-04]

            #pprint(domain.quantities['stage'].vertex_values[:4,0])
            #pprint(domain.quantities['xmomentum'].vertex_values[:4,0])
            #pprint(domain.quantities['ymomentum'].vertex_values[:4,0])


            assert num.allclose(domain.quantities['xmomentum'].vertex_values[:4,0],
                    UH_EX)

            assert num.allclose(domain.quantities['ymomentum'].vertex_values[:4,0],
                    VH_EX)


        os.remove(domain.get_name() + '.sww')

    def test_bedslope_problem_first_order(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 1

        # Bed-slope and friction
        def x_slope(x, y):
            return -x / 3

        domain.set_quantity('elevation', x_slope)

        # Boundary conditions
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        # Initial condition
        domain.set_quantity('stage', expression='elevation+0.05')
        domain.check_integrity()

        # Evolution
        for t in domain.evolve(yieldstep=0.05, finaltime=0.05):
            pass

        # FIXME (Ole): Need some other assertion here!
        #print domain.min_timestep, domain.max_timestep
        #assert allclose(domain.recorded_min_timestep, 0.050010003001)
        #assert allclose(domain.recorded_max_timestep, 0.050010003001)

        os.remove(domain.get_name() + '.sww')

    def test_bedslope_problem_first_order_moresteps_low_froude_0(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 1

        domain.set_low_froude(0)

        # FIXME (Ole): Need tests where these two are commented out
        domain.H0 = 0                         # Backwards compatibility (6/2/7)
        domain.tight_slope_limiters = 0       # Backwards compatibility (14/4/7)
        domain.use_centroid_velocities = 0    # Backwards compatibility (7/5/8)

        # Bed-slope and friction
        def x_slope(x, y):
            return -x / 3

        domain.set_quantity('elevation', x_slope)

        # Boundary conditions
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        # Initial condition
        domain.set_quantity('stage', expression='elevation+0.05')
        domain.check_integrity()

        # Evolution
        for t in domain.evolve(yieldstep=0.05, finaltime=0.5):
            pass


        W_EX = num.array([-0.02933401, -0.01091747, -0.02782667, -0.01057726, -0.02746728,
       -0.01033124, -0.02724838, -0.01019733, -0.02684617, -0.00965007,
       -0.02495054, -0.00706089, -0.06890704, -0.06303447, -0.07568528,
       -0.06161392, -0.07566147, -0.06125286, -0.07537564, -0.06094254,
       -0.07488436, -0.06058909, -0.07465553, -0.06045194, -0.12205956,
       -0.10701804, -0.12318423, -0.10795678, -0.12302672, -0.10763744,
       -0.12279702, -0.10731376, -0.12214772, -0.1071946 , -0.12283814,
       -0.10858945, -0.16389971, -0.15304369, -0.16428248, -0.15298725,
       -0.16430596, -0.15254043, -0.16407102, -0.15236021, -0.16390702,
       -0.15187689, -0.16460111, -0.15549742, -0.2083898 , -0.19912702,
       -0.20845637, -0.19709535, -0.20648985, -0.19663999, -0.20525695,
       -0.19631453, -0.20446298, -0.19578452, -0.20735693, -0.19565282,
       -0.14022868, -0.14262659, -0.13774131, -0.14132395, -0.13707526,
       -0.14041639, -0.13594765, -0.13910709, -0.13533594, -0.1393996 ,
       -0.1328638 , -0.1363085 ])


        #Data from earlier version of abstract_2d_finite_volumes
        assert num.allclose(domain.quantities['stage'].centroid_values, W_EX)


        os.remove(domain.get_name() + '.sww')


    def test_bedslope_problem_first_order_moresteps_low_froude_1(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 1

        domain.set_low_froude(0)

        # FIXME (Ole): Need tests where these two are commented out
        domain.H0 = 0                         # Backwards compatibility (6/2/7)
        domain.tight_slope_limiters = 0       # Backwards compatibility (14/4/7)
        domain.use_centroid_velocities = 0    # Backwards compatibility (7/5/8)

        # Bed-slope and friction
        def x_slope(x, y):
            return -x / 3

        domain.set_quantity('elevation', x_slope)

        # Boundary conditions
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        # Initial condition
        domain.set_quantity('stage', expression='elevation+0.05')
        domain.check_integrity()

        # Evolution
        for t in domain.evolve(yieldstep=0.05, finaltime=0.5):
            pass


        W_EX = num.array([-0.02933401, -0.01091747, -0.02782667, -0.01057726, -0.02746728,
       -0.01033124, -0.02724838, -0.01019733, -0.02684617, -0.00965007,
       -0.02495054, -0.00706089, -0.06890704, -0.06303447, -0.07568528,
       -0.06161392, -0.07566147, -0.06125286, -0.07537564, -0.06094254,
       -0.07488436, -0.06058909, -0.07465553, -0.06045194, -0.12205956,
       -0.10701804, -0.12318423, -0.10795678, -0.12302672, -0.10763744,
       -0.12279702, -0.10731376, -0.12214772, -0.1071946 , -0.12283814,
       -0.10858945, -0.16389971, -0.15304369, -0.16428248, -0.15298725,
       -0.16430596, -0.15254043, -0.16407102, -0.15236021, -0.16390702,
       -0.15187689, -0.16460111, -0.15549742, -0.2083898 , -0.19912702,
       -0.20845637, -0.19709535, -0.20648985, -0.19663999, -0.20525695,
       -0.19631453, -0.20446298, -0.19578452, -0.20735693, -0.19565282,
       -0.14022868, -0.14262659, -0.13774131, -0.14132395, -0.13707526,
       -0.14041639, -0.13594765, -0.13910709, -0.13533594, -0.1393996 ,
       -0.1328638 , -0.1363085 ])


        #Data from earlier version of abstract_2d_finite_volumes
        assert num.allclose(domain.quantities['stage'].centroid_values, W_EX)


        os.remove(domain.get_name() + '.sww')

    def test_bedslope_problem_second_order_one_step(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 2
        domain.beta_w = 0.9
        domain.beta_w_dry = 0.9
        domain.beta_uh = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh = 0.9
        domain.beta_vh_dry = 0.9


        # FIXME (Ole): Need tests where this is commented out
        domain.tight_slope_limiters = 0       # Backwards compatibility (14/4/7)
        domain.use_centroid_velocities = 0    # Backwards compatibility (7/5/8)

        # Bed-slope and friction at vertices (and interpolated elsewhere)
        def x_slope(x, y):
            return -x / 3

        domain.set_quantity('elevation', x_slope)

        # Boundary conditions
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        #Initial condition
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
                             -0.1537037,  -0.13518519, -0.1537037,
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
        for t in domain.evolve(yieldstep=0.05, finaltime=0.05):
            pass



        W_EX = [ 0.01571282,  0.02678718,  0.01765727,  0.02678718,  0.01765727,
        0.02678718,  0.01765727,  0.02678718,  0.01765727,  0.02678718,
        0.01765727,  0.02873162, -0.04259259, -0.02407407, -0.04259259,
       -0.02407407, -0.04259259, -0.02407407, -0.04259259, -0.02407407,
       -0.04259259, -0.02407407, -0.04259259, -0.02407407, -0.09814815,
       -0.07962963, -0.09814815, -0.07962963, -0.09814815, -0.07962963,
       -0.09814815, -0.07962963, -0.09814815, -0.07962963, -0.09814815,
       -0.07962963, -0.1537037 , -0.13518519, -0.1537037 , -0.13518519,
       -0.1537037 , -0.13518519, -0.1537037 , -0.13518519, -0.1537037 ,
       -0.13518519, -0.1537037 , -0.13518519, -0.20925926, -0.19074074,
       -0.20925926, -0.19074074, -0.20925926, -0.19074074, -0.20925926,
       -0.19074074, -0.20925926, -0.19074074, -0.20925926, -0.19074074,
       -0.26206496, -0.2509906 , -0.26012051, -0.2509906 , -0.26012051,
       -0.2509906 , -0.26012051, -0.2509906 , -0.26012051, -0.2509906 ,
       -0.26012051, -0.24904616]


        assert num.allclose(domain.quantities['stage'].centroid_values, W_EX)


        os.remove(domain.get_name() + '.sww')

    def test_bedslope_problem_second_order_two_steps(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)

        domain.set_low_froude(0) 

        domain.smooth = False
        domain.default_order = 2
        domain.beta_w = 0.9
        domain.beta_w_dry = 0.9
        domain.beta_uh = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh = 0.9
        domain.beta_vh_dry = 0.9

        # FIXME (Ole): Need tests where this is commented out
        domain.tight_slope_limiters = 0    # Backwards compatibility (14/4/7)
        domain.H0 = 0    # Backwards compatibility (6/2/7)
        domain.use_centroid_velocities = 0    # Backwards compatibility (7/5/8)

        # Bed-slope and friction at vertices (and interpolated elsewhere)
        def x_slope(x, y):
            return -x / 3

        domain.set_quantity('elevation', x_slope)

        # Boundary conditions
        Br = Reflective_boundary(domain)
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
        for t in domain.evolve(yieldstep=0.05, finaltime=0.1):
            pass

        # Data from earlier version of abstract_2d_finite_volumes ft=0.1
        assert num.allclose(domain.recorded_min_timestep, 0.0344607459654)
        assert num.allclose(domain.recorded_max_timestep, 0.0391090502542)



        W_EX = [ 0.01308246,  0.02201217,  0.01358687,  0.023575  ,  0.01370201,
        0.0235574 ,  0.01370201,  0.02355753,  0.0136998 ,  0.02361447,
        0.01454146,  0.02507853, -0.04254524, -0.022535  , -0.04260359,
       -0.0225144 , -0.04263851, -0.02252512, -0.04263809, -0.02252512,
       -0.04264303, -0.02249883, -0.04257228, -0.02296247, -0.0981472 ,
       -0.07964098, -0.09814502, -0.07968513, -0.09814337, -0.07968807,
       -0.09814296, -0.07968807, -0.09814296, -0.07968807, -0.09814555,
       -0.07965964, -0.15368138, -0.13518778, -0.15366285, -0.13519037,
       -0.15366285, -0.13519037, -0.15366285, -0.13518996, -0.15366479,
       -0.13518832, -0.15369898, -0.13518613, -0.21012122, -0.19071462,
       -0.21066837, -0.19057453, -0.21064807, -0.19057613, -0.21064807,
       -0.19057519, -0.21067682, -0.19061141, -0.21080033, -0.19066023,
       -0.25755282, -0.24897294, -0.25604086, -0.24820174, -0.25598276,
       -0.24820504, -0.25598272, -0.24820504, -0.2559903 , -0.2480719 ,
       -0.25395141, -0.24799187]


        assert num.allclose(domain.quantities['stage'].centroid_values, W_EX)



        os.remove(domain.get_name() + '.sww')

    def test_bedslope_problem_second_order_two_yieldsteps(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        #Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)

        #domain.set_compute_fluxes_method('original')

        domain.set_low_froude(0)

        domain.smooth = False
        domain.default_order = 2
        domain.beta_w = 0.9
        domain.beta_w_dry = 0.9
        domain.beta_uh = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh = 0.9
        domain.beta_vh_dry = 0.9

        # FIXME (Ole): Need tests where this is commented out
        domain.tight_slope_limiters = 0    # Backwards compatibility (14/4/7)
        domain.H0 = 0    # Backwards compatibility (6/2/7)
        domain.use_centroid_velocities = 0    # Backwards compatibility (7/5/8)


        # Bed-slope and friction at vertices (and interpolated elsewhere)
        def x_slope(x, y):
            return -x / 3

        domain.set_quantity('elevation', x_slope)

        # Boundary conditions
        Br = Reflective_boundary(domain)
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
        for t in domain.evolve(yieldstep=0.05, finaltime=0.1):   #0.05??
            pass




        W_EX = [ 0.01308246,  0.02201217,  0.01358687,  0.023575  ,  0.01370201,
        0.0235574 ,  0.01370201,  0.02355753,  0.0136998 ,  0.02361447,
        0.01454146,  0.02507853, -0.04254524, -0.022535  , -0.04260359,
       -0.0225144 , -0.04263851, -0.02252512, -0.04263809, -0.02252512,
       -0.04264303, -0.02249883, -0.04257228, -0.02296247, -0.0981472 ,
       -0.07964098, -0.09814502, -0.07968513, -0.09814337, -0.07968807,
       -0.09814296, -0.07968807, -0.09814296, -0.07968807, -0.09814555,
       -0.07965964, -0.15368138, -0.13518778, -0.15366285, -0.13519037,
       -0.15366285, -0.13519037, -0.15366285, -0.13518996, -0.15366479,
       -0.13518832, -0.15369898, -0.13518613, -0.21012122, -0.19071462,
       -0.21066837, -0.19057453, -0.21064807, -0.19057613, -0.21064807,
       -0.19057519, -0.21067682, -0.19061141, -0.21080033, -0.19066023,
       -0.25755282, -0.24897294, -0.25604086, -0.24820174, -0.25598276,
       -0.24820504, -0.25598272, -0.24820504, -0.2559903 , -0.2480719 ,
       -0.25395141, -0.24799187]


        assert num.allclose(domain.quantities['stage'].centroid_values, W_EX)


        os.remove(domain.get_name() + '.sww')

    def test_bedslope_problem_second_order_more_steps(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)

        #domain.set_compute_fluxes_method('original')

        domain.smooth = False
        domain.default_order = 2
        domain.beta_w = 0.9
        domain.beta_w_dry = 0.9
        domain.beta_uh = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh = 0.9
        domain.beta_vh_dry = 0.9

        domain.set_low_froude(0)

        # FIXME (Ole): Need tests where these two are commented out
        domain.H0 = 0                      # Backwards compatibility (6/2/7)
        domain.tight_slope_limiters = 0    # Backwards compatibility (14/4/7)
        domain.use_centroid_velocities = 0    # Backwards compatibility (7/5/8)

        # Bed-slope and friction at vertices (and interpolated elsewhere)
        def x_slope(x, y):
            return -x / 3

        domain.set_quantity('elevation', x_slope)

        # Boundary conditions
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        # Initial condition
        domain.set_quantity('stage', expression='elevation+0.05')
        domain.check_integrity()



        W_EX = [ 0.01296296,  0.03148148,  0.01296296,  0.03148148,  0.01296296,
        0.03148148,  0.01296296,  0.03148148,  0.01296296,  0.03148148,
        0.01296296,  0.03148148, -0.04259259, -0.02407407, -0.04259259,
       -0.02407407, -0.04259259, -0.02407407, -0.04259259, -0.02407407,
       -0.04259259, -0.02407407, -0.04259259, -0.02407407, -0.09814815,
       -0.07962963, -0.09814815, -0.07962963, -0.09814815, -0.07962963,
       -0.09814815, -0.07962963, -0.09814815, -0.07962963, -0.09814815,
       -0.07962963, -0.1537037 , -0.13518519, -0.1537037 , -0.13518519,
       -0.1537037 , -0.13518519, -0.1537037 , -0.13518519, -0.1537037 ,
       -0.13518519, -0.1537037 , -0.13518519, -0.20925926, -0.19074074,
       -0.20925926, -0.19074074, -0.20925926, -0.19074074, -0.20925926,
       -0.19074074, -0.20925926, -0.19074074, -0.20925926, -0.19074074,
       -0.26481481, -0.2462963 , -0.26481481, -0.2462963 , -0.26481481,
       -0.2462963 , -0.26481481, -0.2462963 , -0.26481481, -0.2462963 ,
       -0.26481481, -0.2462963 ]

        assert num.allclose(domain.quantities['stage'].centroid_values, W_EX)

        # Evolution
        for t in domain.evolve(yieldstep=0.05, finaltime=0.5):
            # Check that diagnostics works
            msg = domain.timestepping_statistics(track_speeds=True)
            #FIXME(Ole): One might check the contents of msg here.







        W_EX = [-0.0301883 , -0.01127593, -0.02834861, -0.0108968 , -0.02806583,
       -0.01074475, -0.02788852, -0.01065176, -0.02753658, -0.01013457,
       -0.02602057, -0.00697983, -0.07333722, -0.06511481, -0.08002821,
       -0.06234795, -0.07884103, -0.06214974, -0.07867692, -0.06199716,
       -0.07848496, -0.06160417, -0.07782372, -0.06214378, -0.12379366,
       -0.10782217, -0.12402788, -0.10970479, -0.12407233, -0.10911802,
       -0.12375966, -0.10870709, -0.1233677 , -0.10897925, -0.1246598 ,
       -0.10913556, -0.16030637, -0.15059851, -0.1619466 , -0.15150383,
       -0.16222821, -0.15126241, -0.16204491, -0.15117909, -0.16170723,
       -0.15047256, -0.16159564, -0.15306895, -0.20650721, -0.19412466,
       -0.20661286, -0.193611  , -0.20611843, -0.1932463 , -0.20476096,
       -0.19336801, -0.20673157, -0.1928249 , -0.20666551, -0.18982871,
       -0.14177628, -0.13867427, -0.138984  , -0.13970989, -0.13837738,
       -0.13785944, -0.13772734, -0.13619726, -0.13708211, -0.13865006,
       -0.13413676, -0.14008116]




        assert num.allclose(domain.quantities['stage'].centroid_values, W_EX)



        UH_EX = [ 0.00439415,  0.00080998,  0.00359639,  0.00086799,  0.00369947,
        0.00089459,  0.00376748,  0.00091612,  0.00394185,  0.00100331,
        0.00423642,  0.00102681,  0.01866626,  0.00841126,  0.01265179,
        0.00989651,  0.01358962,  0.01007976,  0.01375394,  0.01018847,
        0.01388854,  0.01057534,  0.01516146,  0.01155302,  0.03130996,
        0.02572364,  0.03012732,  0.02348833,  0.02999708,  0.02419562,
        0.03046757,  0.02471555,  0.03088795,  0.02415592,  0.02886438,
        0.02425477,  0.06609542,  0.04946308,  0.06236933,  0.04715456,
        0.06176596,  0.04743353,  0.06204203,  0.04760272,  0.06262964,
        0.04858867,  0.06290163,  0.0447266 ,  0.08245285,  0.07715851,
        0.08533206,  0.07817769,  0.08610624,  0.07901377,  0.0879718 ,
        0.07879576,  0.08550203,  0.07963651,  0.08581205,  0.08376399,
        0.01997027,  0.07948531,  0.02703306,  0.08443647,  0.02694374,
        0.08404407,  0.02627307,  0.08400949,  0.02649757,  0.08397135,
        0.03080026,  0.09688469]




        assert num.allclose(domain.quantities['xmomentum'].centroid_values, UH_EX)





        VH_EX = [ -4.55438497e-04,   5.60449890e-05,  -1.05041759e-04,
         1.09474097e-04,   3.37384149e-05,   1.48217622e-04,
         8.15006870e-05,   1.41778277e-04,  -8.95699377e-05,
        -1.50989506e-05,  -1.13309713e-03,  -6.56027069e-04,
         4.08407333e-04,  -5.37824083e-04,  -6.07324848e-04,
        -1.79169050e-04,  -2.26068995e-04,  -9.40081480e-05,
        -1.92103118e-04,  -7.38055580e-05,  -2.61634534e-04,
        -2.18534710e-04,  -3.97532842e-04,  -4.38793389e-04,
         4.02437237e-04,   1.71996550e-04,  -4.71404571e-04,
        -7.07195798e-04,  -2.92606384e-04,  -2.89749829e-04,
        -2.14551908e-04,  -3.47586675e-04,  -4.51408518e-04,
        -2.85326203e-04,   2.60585517e-04,   1.85996351e-04,
         1.76070321e-03,   4.27131620e-04,  -1.36487970e-04,
        -7.48486546e-04,  -2.25542208e-04,  -4.38919415e-04,
        -1.23010639e-05,  -3.49340425e-04,  -2.36098938e-04,
        -7.63395297e-04,   6.66099385e-05,   6.22983544e-04,
         4.52833391e-04,   2.24808947e-03,   4.85972740e-04,
         7.50193061e-04,   4.53207970e-04,   1.05078134e-03,
         6.32114704e-04,   1.24869321e-03,   4.06577106e-04,
         1.01553335e-03,  -5.26250901e-04,   1.78062345e-04,
        -2.04979836e-03,  -2.70644308e-03,  -3.34067897e-03,
        -1.94716541e-03,  -2.28075505e-03,  -1.61494727e-03,
        -2.15457373e-03,  -1.48932625e-03,  -3.24804437e-03,
        -1.09719715e-03,   3.38706650e-03,  -8.98151209e-04]

        #pprint(domain.quantities['ymomentum'].centroid_values)


        assert num.allclose(domain.quantities['ymomentum'].centroid_values, VH_EX)

        os.remove(domain.get_name() + '.sww')

    def NOtest_bedslope_problem_second_order_more_steps_feb_2007(self):
        """test_bedslope_problem_second_order_more_steps_feb_2007

        Test shallow water finite volumes, using parameters from
        feb 2007 rather than backward compatibility ad infinitum
        """

        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(6, 6)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 2
        domain.beta_w = 0.9
        domain.beta_w_dry = 0.9
        domain.beta_uh = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh = 0.9
        domain.beta_vh_dry = 0.9
        domain.H0 = 0.001
        domain.tight_slope_limiters = 1

        # Bed-slope and friction at vertices (and interpolated elsewhere)
        def x_slope(x, y):
            return -x / 3

        domain.set_quantity('elevation', x_slope)

        # Boundary conditions
        Br = Reflective_boundary(domain)
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
         [-0.03348416, -0.01749303, -0.03299091, -0.01739241, -0.03246447, -0.01732016,
          -0.03205390, -0.01717833, -0.03146383, -0.01699831, -0.03076577, -0.01671795,
          -0.07952656, -0.06684763, -0.07721455, -0.06668388, -0.07632976, -0.06600113,
          -0.07523678, -0.06546373, -0.07447040, -0.06508861, -0.07438723, -0.06359288,
          -0.12526729, -0.11205668, -0.12179433, -0.11068104, -0.12048395, -0.10968948,
          -0.11912023, -0.10862628, -0.11784090, -0.10803744, -0.11790629, -0.10742354,
          -0.16859613, -0.15427413, -0.16664444, -0.15464452, -0.16570816, -0.15327556,
          -0.16409162, -0.15204092, -0.16264608, -0.15102139, -0.16162736, -0.14969205,
          -0.18736511, -0.19874036, -0.18811230, -0.19758289, -0.18590182, -0.19580301,
          -0.18234588, -0.19423215, -0.18100376, -0.19380116, -0.18509710, -0.19501636,
          -0.13982382, -0.14166819, -0.14132775, -0.14528694, -0.14096905, -0.14351126,
          -0.13800356, -0.14027920, -0.13613538, -0.13936795, -0.13621902, -0.14204982])


        assert num.allclose(domain.quantities['xmomentum'].centroid_values,
      [0.00600290,  0.00175780,  0.00591905,  0.00190903,  0.00644462,  0.00203095,
       0.00684561,  0.00225089,  0.00708208,  0.00236235,  0.00649095,  0.00222343,
       0.02068693,  0.01164034,  0.01983343,  0.01159526,  0.02044611,  0.01233252,
       0.02135685,  0.01301289,  0.02161290,  0.01260280,  0.01867612,  0.01133078,
       0.04091313,  0.02668283,  0.03634781,  0.02733469,  0.03767692,  0.02836840,
       0.03906338,  0.02958073,  0.04025669,  0.02953292,  0.03665616,  0.02583565,
       0.06314558,  0.04830935,  0.05663609,  0.04564362,  0.05756200,  0.04739673,
       0.05967379,  0.04919083,  0.06124330,  0.04965808,  0.05879240,  0.04629319,
       0.08220739,  0.06924725,  0.07713556,  0.06782640,  0.07909499,  0.06992544,
       0.08116621,  0.07210181,  0.08281548,  0.07222669,  0.07941059,  0.06755612,
       0.01581588,  0.04533609,  0.02017939,  0.04342565,  0.02073232,  0.04476108,
       0.02117439,  0.04573358,  0.02129473,  0.04694267,  0.02220398,  0.05533458])


        assert num.allclose(domain.quantities['ymomentum'].centroid_values,
     [-7.65882069e-005, -1.46087080e-004, -1.09630102e-004, -7.80950424e-005,
      -1.15922807e-005, -9.09134899e-005, -1.35994542e-004, -1.95673476e-004,
      -4.25779199e-004, -2.95890312e-004, -4.00060341e-004, -9.42021290e-005,
      -3.41372596e-004, -1.54560195e-004, -2.94810038e-004, -1.08844546e-004,
      -6.97240892e-005,  3.50299623e-005, -2.40159184e-004, -2.01805883e-004,
      -7.60732405e-004, -5.10897642e-004, -1.00940001e-003, -1.38037759e-004,
      -1.06169131e-003, -3.12307760e-004, -9.90602307e-004, -4.21634250e-005,
      -6.02424239e-004,  1.52230578e-004, -7.63833035e-004, -1.10273481e-004,
      -1.40187071e-003, -5.57831837e-004, -1.63988285e-003, -2.48018092e-004,
      -1.83309840e-003, -6.19360836e-004, -1.29955242e-003, -3.76237145e-004,
      -1.00613007e-003, -8.63641918e-005, -1.13604124e-003, -3.90589728e-004,
      -1.91457355e-003, -9.43783961e-004, -2.28090840e-003, -5.79107025e-004,
      -1.54091533e-003, -2.39785792e-003, -2.47947427e-003, -2.02694009e-003,
      -2.10441194e-003, -1.82082650e-003, -1.80229336e-003, -2.10418336e-003,
      -1.93104408e-003, -2.23200334e-003, -1.57239706e-003, -1.31486358e-003,
      -1.17564993e-003, -2.85846494e-003, -3.52956754e-003, -5.12658193e-003,
      -6.24238960e-003, -6.01820113e-003, -6.09602201e-003, -5.04787190e-003,
      -4.59373845e-003, -3.01393146e-003,  5.08550095e-004, -4.35896549e-004])

        os.remove(domain.get_name() + '.sww')

    def test_temp_play(self):
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(5, 5)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 2
        domain.beta_w = 0.9
        domain.beta_w_dry = 0.9
        domain.beta_uh = 0.9
        domain.beta_uh_dry = 0.9
        domain.beta_vh = 0.9
        domain.beta_vh_dry = 0.9

        # FIXME (Ole): Need tests where these two are commented out
        domain.H0 = 0                         # Backwards compatibility (6/2/7)
        domain.tight_slope_limiters = False   # Backwards compatibility (14/4/7)
        domain.use_centroid_velocities = False # Backwards compatibility (7/5/8)
        domain.use_edge_limiter = False       # Backwards compatibility (9/5/8)
        domain.low_froude = 0                 # Backwards compatibility (25/6/19)

        # Bed-slope and friction at vertices (and interpolated elsewhere)
        def x_slope(x, y):
            return -x / 3

        domain.set_quantity('elevation', x_slope)

        # Boundary conditions
        Br = Reflective_boundary(domain)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        # Initial condition
        domain.set_quantity('stage', expression='elevation+0.05')
        domain.check_integrity()

        # Evolution
        for t in domain.evolve(yieldstep=0.05, finaltime=0.1):
            pass

        #pprint(domain.quantities['stage'].centroid_values[:4])
        #pprint(domain.quantities['xmomentum'].centroid_values[:4])
        #pprint(domain.quantities['ymomentum'].centroid_values[:4])

        W = domain.quantities['stage'].centroid_values[:4]
        UH = domain.quantities['xmomentum'].centroid_values[:4]
        VH = domain.quantities['ymomentum'].centroid_values[:4]

        W_0  = [ 0.001362,    0.01344294,  0.00308829, 0.01470289]
        UH_0 = [ 0.01300239,  0.00537933,  0.01214676,  0.00515825]
        VH_0 = [ -1.13165691e-03,  -6.55330189e-04, -6.62804076e-05,   5.26313051e-05]

        W_1 = [ 0.00707892,  0.01849914,  0.00783274,  0.01997863]
        UH_1 = [ 0.01512518,  0.00354391,  0.01503765,  0.00326075]
        VH_1 = [  5.36531332e-04,  -6.77297008e-04,  -4.58560426e-05, 2.47714988e-05]


        assert num.allclose(W,W_0) or num.allclose(W,W_1)
        assert num.allclose(UH, UH_0) or num.allclose(UH, UH_1)
        assert num.allclose(VH, VH_0) or num.allclose(VH, VH_1)




        # old values pre revision 8402
#        assert num.allclose(domain.quantities['stage'].centroid_values[:4],
#                            [0.00206836, 0.01296714, 0.00363415, 0.01438924])
#        assert num.allclose(domain.quantities['xmomentum'].centroid_values[:4],
#                            [0.01360154, 0.00671133, 0.01264578, 0.00648503])
#        assert num.allclose(domain.quantities['ymomentum'].centroid_values[:4],
#                            [-1.19201077e-003, -7.23647546e-004,
#                             -6.39083123e-005, 6.29815168e-005])

        os.remove(domain.get_name() + '.sww')

    def test_complex_bed(self):
        # No friction is tested here

        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        N = 12
        points, vertices, boundary = rectangular(N, N // 2, len1=1.2, len2=0.6,
                                                 origin=(-0.07, 0))


        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 2

        inflow_stage = 0.1
        Z = Weir(inflow_stage)
        domain.set_quantity('elevation', Z)

        Br = Reflective_boundary(domain)
        Bd = Dirichlet_boundary([inflow_stage, 0.0, 0.0])
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

    def test_spatio_temporal_boundary_1(self):
        """Test that boundary values can be read from file and interpolated
        in both time and space.

        Verify that the same steady state solution is arrived at and that
        time interpolation works.

        The full solution history is not exactly the same as
        file boundary must read and interpolate from *smoothed* version
        as stored in sww.
        """

        # Create sww file of simple propagation from left to right
        # through rectangular domain
        verbose = False

        import time
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create basic mesh
        points, vertices, boundary = rectangular(3, 3)

        # Create shallow water domain
        domain1 = Domain(points, vertices, boundary)

        domain1.set_low_froude(0)
        domain1.reduction = mean
        domain1.smooth = False    # Exact result

        domain1.default_order = 2
        domain1.store = True
        domain1.set_datadir('.')
        domain1.set_name('spatio_temporal_boundary_source' + str(time.time()))

        # FIXME: This is extremely important!
        # How can we test if they weren't stored?


        # Bed-slope and friction at vertices (and interpolated elsewhere)
        domain1.set_quantity('elevation', 0)
        domain1.set_quantity('friction', 0)

        # Boundary conditions
        Br = Reflective_boundary(domain1)
        Bd = Dirichlet_boundary([0.3,0,0])
        domain1.set_boundary({'left': Bd, 'top': Bd, 'right': Br, 'bottom': Br})

        # Initial condition
        domain1.set_quantity('stage', 0)
        domain1.check_integrity()

        finaltime = 10

        # Evolution  (full domain - large steps)
        for t in domain1.evolve(yieldstep=0.671, finaltime=finaltime):
            pass

        cv1 = domain1.quantities['stage'].centroid_values

        # Create a triangle shaped domain (reusing coordinates from domain 1),
        # formed from the lower and right hand  boundaries and
        # the sw-ne diagonal
        # from domain 1. Call it domain2

        points = [[0,0], [1.0/3,0], [1.0/3,1.0/3],
                  [2.0/3,0], [2.0/3,1.0/3], [2.0/3,2.0/3],
                  [1,0], [1,1.0/3], [1,2.0/3], [1,1]]

        vertices = [[1,2,0], [3,4,1], [2,1,4], [4,5,2],
                    [6,7,3], [4,3,7], [7,8,4], [5,4,8], [8,9,5]]

        boundary = {(0,1): 'bottom',   (1,1): 'bottom',   (4,1): 'bottom',
                    (4,2): 'right',    (6,2): 'right',    (8,2): 'right',
                    (0,0): 'diagonal', (3,0): 'diagonal', (8,0): 'diagonal'}

        domain2 = Domain(points, vertices, boundary)

        domain2.low_froude = domain1.low_froude
        domain2.reduction = domain1.reduction
        domain2.smooth = False
        domain2.default_order = 2

        # Bed-slope and friction at vertices (and interpolated elsewhere)
        domain2.set_quantity('elevation', 0)
        domain2.set_quantity('friction', 0)
        domain2.set_quantity('stage', 0)

        # Boundary conditions
        Br = Reflective_boundary(domain2)
        Bf = Field_boundary(domain1.get_name() + '.sww', domain2)
        domain2.set_boundary({'right':Br, 'bottom':Br, 'diagonal':Bf})
        domain2.check_integrity()

        # Evolution (small steps)
        for t in domain2.evolve(yieldstep=0.0711, finaltime=finaltime):
            pass

        # Use output from domain1 as spatio-temporal boundary for domain2
        # and verify that results at right hand side are close.
        cv2 = domain2.quantities['stage'].centroid_values

        assert num.allclose(num.take(cv1, (0,8,16), axis=0),
                            num.take(cv2, (0,3,8), axis=0), rtol=1.0e-4)      # Diag
        assert num.allclose(num.take(cv1, (0,6,12), axis=0),
                            num.take(cv2, (0,1,4), axis=0), rtol=1.0e-4)      # Bottom
        assert num.allclose(num.take(cv1, (12,14,16), axis=0),
                            num.take(cv2, (4,6,8), axis=0), rtol=1.0e-4)      # RHS

        # Cleanup
        os.remove(domain1.get_name() + '.sww')
        os.remove(domain2.get_name() + '.sww')

    def test_spatio_temporal_boundary_2(self):
        """Test that boundary values can be read from file and interpolated
        in both time and space.
        This is a more basic test, verifying that boundary object
        produces the expected results
        """

        import time
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create sww file of simple propagation from left to right
        # through rectangular domain

        # Create basic mesh
        points, vertices, boundary = rectangular(3, 3)

        #Create shallow water domain
        domain1 = Domain(points, vertices, boundary)

        domain1.reduction = mean
        domain1.smooth = True    # To mimic MOST output

        domain1.default_order = 2
        domain1.store = True
        domain1.set_datadir('.')
        domain1.set_name('spatio_temporal_boundary_source' + str(time.time()))

        # FIXME: This is extremely important!
        # How can we test if they weren't stored?
        domain1.quantities_to_be_stored = {'elevation': 1,
                                           'stage': 2,
                                           'xmomentum': 2,
                                           'ymomentum': 2}

        # Bed-slope and friction at vertices (and interpolated elsewhere)
        domain1.set_quantity('elevation', 0)
        domain1.set_quantity('friction', 0)

        # Boundary conditions
        Br = Reflective_boundary(domain1)
        Bd = Dirichlet_boundary([0.3,0,0])
        domain1.set_boundary({'left': Bd, 'top': Bd, 'right': Br, 'bottom': Br})

        # Initial condition
        domain1.set_quantity('stage', 0)
        domain1.check_integrity()

        finaltime = 5

        # Evolution  (full domain - large steps)
        for t in domain1.evolve(yieldstep=1, finaltime=finaltime):
            pass

        # Create an triangle shaped domain (coinciding with some
        # coordinates from domain 1),
        # formed from the lower and right hand  boundaries and
        # the sw-ne diagonal
        # from domain 1. Call it domain2
        points = [[0,0],
                  [1.0/3,0], [1.0/3,1.0/3],
                  [2.0/3,0], [2.0/3,1.0/3], [2.0/3,2.0/3],
                  [1,0],     [1,1.0/3],     [1,2.0/3],     [1,1]]

        vertices = [[1,2,0], [3,4,1], [2,1,4], [4,5,2],
                    [6,7,3], [4,3,7], [7,8,4], [5,4,8], [8,9,5]]

        boundary = {(0,1): 'bottom',   (1,1): 'bottom',   (4,1): 'bottom',
                    (4,2): 'right',    (6,2): 'right',    (8,2): 'right',
                    (0,0): 'diagonal', (3,0): 'diagonal', (8,0): 'diagonal'}

        domain2 = Domain(points, vertices, boundary)

        domain2.reduction = domain1.reduction
        domain2.smooth = False
        domain2.default_order = 2

        # Bed-slope and friction at vertices (and interpolated elsewhere)
        domain2.set_quantity('elevation', 0)
        domain2.set_quantity('friction', 0)
        domain2.set_quantity('stage', 0)

        # Read results for specific timesteps t=1 and t=2
        fid = NetCDFFile(domain1.get_name() + '.sww')

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        s1 = fid.variables['stage'][1,:]
        s2 = fid.variables['stage'][2,:]
        fid.close()

        shp = (len(x), 1)
        points = num.concatenate((num.reshape(x, shp), num.reshape(y, shp)),
                                 axis=1)

        # The diagonal points of domain 1 are 0, 5, 10, 15
        msg = ('value was\n%s\nshould be\n'
               '[[0,0], [1.0/3, 1.0/3],\n'
               '[2.0/3, 2.0/3], [1,1]]'
               % str(num.take(points, [0,5,10,15], axis=0)))
        assert num.allclose(num.take(points, [0,5,10,15], axis=0),
                            [[0,0], [1.0/3, 1.0/3], [2.0/3, 2.0/3], [1,1]]), msg

        # Boundary conditions
        Br = Reflective_boundary(domain2)
        Bf = Field_boundary(domain1.get_name() + '.sww',
                            domain2, verbose=False)
        domain2.set_boundary({'right':Br, 'bottom':Br, 'diagonal':Bf})
        domain2.check_integrity()

        # Test that interpolation points are the mid points of the all boundary
        # segments
        boundary_midpoints = [[1.0/6, 0], [1.0/2, 0], [5.0/6,0],
                              [1.0, 1.0/6], [1.0, 1.0/2], [1.0, 5.0/6],
                              [1.0/6, 1.0/6], [0.5, 0.5], [5.0/6, 5.0/6]]

        boundary_midpoints.sort()
        R = Bf.F.interpolation_points.tolist()
        R.sort()
        assert num.allclose(boundary_midpoints, R)

        # Check spatially interpolated output at time == 1
        domain2.set_time(1)

        # First diagonal midpoint
        R0 = Bf.evaluate(0, 0)
        assert num.allclose(R0[0], (s1[0] + s1[5]) / 2)

        # Second diagonal midpoint
        R0 = Bf.evaluate(3, 0)
        assert num.allclose(R0[0], (s1[5] + s1[10]) / 2)

        # First diagonal midpoint
        R0 = Bf.evaluate(8, 0)
        assert num.allclose(R0[0], (s1[10] + s1[15]) / 2)

        # Check spatially interpolated output at time == 2
        domain2.set_time(2)

        # First diagonal midpoint
        R0 = Bf.evaluate(0, 0)
        assert num.allclose(R0[0], (s2[0] + s2[5]) / 2)

        # Second diagonal midpoint
        R0 = Bf.evaluate(3, 0)
        assert num.allclose(R0[0], (s2[5] + s2[10]) / 2)

        # First diagonal midpoint
        R0 = Bf.evaluate(8, 0)
        assert num.allclose(R0[0], (s2[10] + s2[15]) / 2)

        # Now check temporal interpolation
        domain2.set_time(1 + 2.0/3)

        # First diagonal midpoint
        R0 = Bf.evaluate(0,0)
        assert num.allclose(R0[0],
                            ((s1[0] + s1[5]) / 2 + 2.0*(s2[0] + s2[5]) / 2) / 3)

        # Second diagonal midpoint
        R0 = Bf.evaluate(3, 0)
        assert num.allclose(R0[0],
                            ((s1[5] + s1[10]) / 2 + 2.0*(s2[5] + s2[10]) / 2) / 3)

        # First diagonal midpoint
        R0 = Bf.evaluate(8, 0)
        assert num.allclose(R0[0],
                            ((s1[10] + s1[15]) / 2 + 2.0*(s2[10] + s2[15]) / 2) / 3)

        # Cleanup
        os.remove(domain1.get_name() + '.sww')

    def test_spatio_temporal_boundary_3(self):
        """Test that boundary values can be read from file and interpolated
        in both time and space.
        This is a more basic test, verifying that boundary object
        produces the expected results

        This tests adjusting using mean_stage
        """

        #Create sww file of simple propagation from left to right
        #through rectangular domain

        import time
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        mean_stage = 5.2    # Adjust stage by this amount in boundary

        # Create basic mesh
        points, vertices, boundary = rectangular(3, 3)

        # Create shallow water domain
        domain1 = Domain(points, vertices, boundary)

        domain1.reduction = mean
        domain1.smooth = True    # To mimic MOST output

        domain1.default_order = 2
        domain1.store = True
        domain1.set_datadir('.')
        domain1.set_name('spatio_temporal_boundary_source' + str(time.time()))

        # FIXME: This is extremely important!
        # How can we test if they weren't stored?

        # Bed-slope and friction at vertices (and interpolated elsewhere)
        domain1.set_quantity('elevation', 0)
        domain1.set_quantity('friction', 0)

        # Boundary conditions
        Br = Reflective_boundary(domain1)
        Bd = Dirichlet_boundary([0.3, 0, 0])
        domain1.set_boundary({'left': Bd, 'top': Bd, 'right': Br, 'bottom': Br})

        # Initial condition
        domain1.set_quantity('stage', 0)
        domain1.check_integrity()

        finaltime = 5

        # Evolution  (full domain - large steps)
        for t in domain1.evolve(yieldstep=1, finaltime=finaltime):
            pass

        # Create an triangle shaped domain (coinciding with some
        # coordinates from domain 1),
        # formed from the lower and right hand  boundaries and
        # the sw-ne diagonal
        # from domain 1. Call it domain2
        points = [[0,0],
                  [1.0/3,0], [1.0/3,1.0/3],
                  [2.0/3,0], [2.0/3,1.0/3], [2.0/3,2.0/3],
                  [1,0],     [1,1.0/3],     [1,2.0/3],     [1,1]]

        vertices = [[1,2,0],
                    [3,4,1], [2,1,4], [4,5,2],
                    [6,7,3], [4,3,7], [7,8,4], [5,4,8], [8,9,5]]

        boundary = {(0,1): 'bottom',   (1,1): 'bottom',   (4,1): 'bottom',
                    (4,2): 'right',    (6,2): 'right',    (8,2): 'right',
                    (0,0): 'diagonal', (3,0): 'diagonal', (8,0): 'diagonal'}

        domain2 = Domain(points, vertices, boundary)

        domain2.reduction = domain1.reduction
        domain2.smooth = False
        domain2.default_order = 2

        # Bed-slope and friction at vertices (and interpolated elsewhere)
        domain2.set_quantity('elevation', 0)
        domain2.set_quantity('friction', 0)
        domain2.set_quantity('stage', 0)

        # Read results for specific timesteps t=1 and t=2
        fid = NetCDFFile(domain1.get_name() + '.sww')

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        s1 = fid.variables['stage'][1,:]
        s2 = fid.variables['stage'][2,:]
        fid.close()

        shp = (len(x), 1)
        points = num.concatenate((num.reshape(x, shp), num.reshape(y, shp)),
                                 axis=1)
        #The diagonal points of domain 1 are 0, 5, 10, 15

        msg = ('values was\n%s\nshould be\n'
               '[[0,0], [1.0/3, 1.0/3],\n'
               '[2.0/3, 2.0/3], [1,1]]'
               % str(num.take(points, [0,5,10,15], axis=0)))
        assert num.allclose(num.take(points, [0,5,10,15], axis=0),
                            [[0,0], [1.0/3, 1.0/3], [2.0/3, 2.0/3], [1,1]]), msg

        # Boundary conditions
        Br = Reflective_boundary(domain2)
        Bf = Field_boundary(domain1.get_name() + '.sww',
                            domain2, mean_stage=mean_stage, verbose=False)

        domain2.set_boundary({'right':Br, 'bottom':Br, 'diagonal':Bf})
        domain2.check_integrity()

        # Test that interpolation points are the mid points of the all boundary
        # segments
        boundary_midpoints = [[1.0/6, 0], [1.0/2, 0], [5.0/6,0],
                              [1.0, 1.0/6], [1.0, 1.0/2], [1.0, 5.0/6],
                              [1.0/6, 1.0/6], [0.5, 0.5], [5.0/6, 5.0/6]]

        boundary_midpoints.sort()
        R = Bf.F.interpolation_points.tolist()
        R.sort()
        assert num.allclose(boundary_midpoints, R)

        # Check spatially interpolated output at time == 1
        domain2.set_time(1)

        # First diagonal midpoint
        R0 = Bf.evaluate(0, 0)
        assert num.allclose(R0[0], (s1[0] + s1[5]) / 2 + mean_stage)

        # Second diagonal midpoint
        R0 = Bf.evaluate(3, 0)
        assert num.allclose(R0[0], (s1[5] + s1[10]) / 2 + mean_stage)

        # First diagonal midpoint
        R0 = Bf.evaluate(8, 0)
        assert num.allclose(R0[0], (s1[10] + s1[15]) / 2 + mean_stage)

        # Check spatially interpolated output at time == 2
        domain2.set_time(2)

        # First diagonal midpoint
        R0 = Bf.evaluate(0, 0)
        assert num.allclose(R0[0], (s2[0] + s2[5]) / 2 + mean_stage)

        # Second diagonal midpoint
        R0 = Bf.evaluate(3, 0)
        assert num.allclose(R0[0], (s2[5] + s2[10]) / 2 + mean_stage)

        # First diagonal midpoint
        R0 = Bf.evaluate(8, 0)
        assert num.allclose(R0[0], (s2[10] + s2[15]) / 2 + mean_stage)

        #Now check temporal interpolation
        domain2.set_time(1 + 2.0/3)

        # First diagonal midpoint
        R0 = Bf.evaluate(0, 0)
        assert num.allclose(R0[0],
                            ((s1[0] + s1[5]) / 2 + 2.0*(s2[0] + s2[5]) / 2) / 3 +
                                mean_stage)

        # Second diagonal midpoint
        R0 = Bf.evaluate(3, 0)
        assert num.allclose(R0[0],
                            ((s1[5] + s1[10]) / 2 + 2.0*(s2[5] + s2[10]) / 2) / 3 +
                                mean_stage)

        # First diagonal midpoint
        R0 = Bf.evaluate(8, 0)
        assert num.allclose(R0[0],
                            ((s1[10] + s1[15]) / 2 + 2.0*(s2[10] + s2[15]) / 2) / 3 +
                                mean_stage)

        # Cleanup
        os.remove(domain1.get_name() + '.sww')

    def test_spatio_temporal_boundary_outside(self):
        """Test that field_boundary catches if a point is outside the sww
        that defines it
        """

        import time
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        # Create sww file of simple propagation from left to right
        # through rectangular domain

        # Create basic mesh
        points, vertices, boundary = rectangular(3, 3)

        # Create shallow water domain
        domain1 = Domain(points, vertices, boundary)

        domain1.reduction = mean
        domain1.smooth = True    # To mimic MOST output

        domain1.default_order = 2
        domain1.store = True
        domain1.set_datadir('.')
        domain1.set_name('spatio_temporal_boundary_source' + str(time.time()))

        # FIXME: This is extremely important!
        # How can we test if they weren't stored?
        domain1.set_quantities_to_be_stored({'elevation': 1,
                                             'stage': 2,
                                             'xmomentum': 2,
                                             'ymomentum': 2})


        # Bed-slope and friction at vertices (and interpolated elsewhere)
        domain1.set_quantity('elevation', 0)
        domain1.set_quantity('friction', 0)

        # Boundary conditions
        Br = Reflective_boundary(domain1)
        Bd = Dirichlet_boundary([0.3, 0, 0])
        domain1.set_boundary({'left': Bd, 'top': Bd, 'right': Br, 'bottom': Br})

        # Initial condition
        domain1.set_quantity('stage', 0)
        domain1.check_integrity()

        finaltime = 5

        # Evolution  (full domain - large steps)
        for t in domain1.evolve(yieldstep=1, finaltime=finaltime):
            pass

        # Create an triangle shaped domain (coinciding with some
        # coordinates from domain 1, but one edge outside!),
        # formed from the lower and right hand  boundaries and
        # the sw-ne diagonal as in the previous test but scaled
        # in the x direction by a factor of 2
        points = [[0,0],
                  [2.0/3,0], [2.0/3,1.0/3],
                  [4.0/3,0], [4.0/3,1.0/3], [4.0/3,2.0/3],
                  [2,0],     [2,1.0/3],     [2,2.0/3],     [2,1]]

        vertices = [[1,2,0],
                    [3,4,1], [2,1,4], [4,5,2],
                    [6,7,3], [4,3,7], [7,8,4], [5,4,8], [8,9,5]]

        boundary = {(0,1): 'bottom',   (1,1): 'bottom',   (4,1): 'bottom',
                    (4,2): 'right',    (6,2): 'right',    (8,2): 'right',
                    (0,0): 'diagonal', (3,0): 'diagonal', (8,0): 'diagonal'}

        domain2 = Domain(points, vertices, boundary)

        domain2.reduction = domain1.reduction
        domain2.smooth = False
        domain2.default_order = 2

        # Bed-slope and friction at vertices (and interpolated elsewhere)
        domain2.set_quantity('elevation', 0)
        domain2.set_quantity('friction', 0)
        domain2.set_quantity('stage', 0)

        # Read results for specific timesteps t=1 and t=2
        fid = NetCDFFile(domain1.get_name() + '.sww')

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        s1 = fid.variables['stage'][1,:]
        s2 = fid.variables['stage'][2,:]
        fid.close()

        shp = (len(x), 1)
        points = num.concatenate((num.reshape(x, shp), num.reshape(y, shp)),
                                 axis=1)
        #The diagonal points of domain 1 are 0, 5, 10, 15

        assert num.allclose(num.take(points, [0,5,10,15], axis=0),
                            [[0,0], [1.0/3,1.0/3], [2.0/3,2.0/3], [1,1]])

        # Boundary conditions
        Br = Reflective_boundary(domain2)
        Bf = Field_boundary(domain1.get_name() + '.sww',
                            domain2, mean_stage=1, verbose=False)

        domain2.set_boundary({'right':Br, 'bottom':Br, 'diagonal':Bf})
        domain2.check_integrity()

        try:
            for t in domain2.evolve(yieldstep=1, finaltime=finaltime):
                pass
        except:
            pass
        else:
            msg = 'This should have caught NAN at boundary'
            raise Exception(msg)

        #Cleanup
        os.remove(domain1.get_name() + '.sww')

    def test_extrema(self):
        """Test that extrema of quantities are computed correctly
        Extrema are updated at every *internal* timestep
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
        domain.set_name('extrema_test')

        #--------------------------------------------------------------
        # Setup initial conditions
        #--------------------------------------------------------------
        def topography(x,y):
            return -x / 2                             # linear bed slope

        domain.set_quantity('elevation', topography)    # function for elevation
        domain.set_quantity('friction', 0.)             # Zero friction
        # Constant negative initial stage
        domain.set_quantity('stage', initial_runup_height)
        domain.set_quantities_to_be_monitored(['stage', 'stage-elevation'],
                                              time_interval=[0.5, 2.7],
                                              polygon=[[0,0], [0,1],
                                                       [1,1], [1,0]])

        assert len(domain.quantities_to_be_monitored) == 2
        assert 'stage' in domain.quantities_to_be_monitored
        assert 'stage-elevation' in domain.quantities_to_be_monitored
        for key in list(domain.quantities_to_be_monitored['stage'].keys()):
            assert domain.quantities_to_be_monitored['stage'][key] is None

        #--------------------------------------------------------------
        # Setup boundary conditions
        #--------------------------------------------------------------
        Br = Reflective_boundary(domain)              # Reflective wall
        # Constant inflow
        Bd = Dirichlet_boundary([final_runup_height, 0, 0])

        # All reflective to begin with (still water)
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

        #--------------------------------------------------------------
        # Let triangles adjust and check extrema
        #--------------------------------------------------------------
        for t in domain.evolve(yieldstep=0.1, finaltime=1.0):
            domain.quantity_statistics() # Run it silently

        #--------------------------------------------------------------
        # Test extrema
        #--------------------------------------------------------------
        stage = domain.quantities_to_be_monitored['stage']
        assert stage['min'] <= stage['max']

        assert num.allclose(stage['min'], initial_runup_height,
                            rtol=1.0/N)    # First order accuracy

        depth = domain.quantities_to_be_monitored['stage-elevation']
        assert depth['min'] <= depth['max']
        assert depth['min'] >= 0.0
        assert depth['max'] >= 0.0

        #--------------------------------------------------------------
        # Update boundary to allow inflow
        #--------------------------------------------------------------
        domain.set_boundary({'right': Bd})

        #--------------------------------------------------------------
        # Evolve system through time
        #--------------------------------------------------------------
        for t in domain.evolve(yieldstep=0.1, finaltime=3.0):
            domain.quantity_statistics()    # Run it silently

        #--------------------------------------------------------------
        # Test extrema again
        #--------------------------------------------------------------
        stage = domain.quantities_to_be_monitored['stage']
        assert stage['min'] <= stage['max']

        assert num.allclose(stage['min'], initial_runup_height,
                            rtol = 1.0/N) # First order accuracy

        depth = domain.quantities_to_be_monitored['stage-elevation']
        assert depth['min'] <= depth['max']
        assert depth['min'] >= 0.0
        assert depth['max'] >= 0.0

        # Cleanup
        os.remove(domain.get_name() + '.sww')

    def test_tight_slope_limiters(self):
        """Test that new slope limiters (Feb 2007) don't induce extremely
        small timesteps. This test actually reveals the problem as it
        was in March-April 2007
        """

        # Create basic mesh
        points, vertices, boundary = rectangular(2, 2)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.default_order = 2

        # This will pass
        #domain.tight_slope_limiters = 1
        #domain.H0 = 0.01

        # This will fail
        #domain.tight_slope_limiters = 1
        #domain.H0 = 0.001

        # This will pass provided C extension implements limiting of
        # momentum in _compute_speeds
        domain.tight_slope_limiters = 1
        domain.H0 = 0.001
        domain.protect_against_isolated_degenerate_timesteps = True

        # Set some field values
        domain.set_quantity('elevation', lambda x,y: -x)
        domain.set_quantity('friction', 0.03)

        # Boundary conditions
        B = Transmissive_boundary(domain)
        domain.set_boundary({'left': B, 'right': B, 'top': B, 'bottom': B})

        # Initial condition - with jumps
        bed = domain.quantities['elevation'].vertex_values
        stage = num.zeros(bed.shape, float)

        h = 0.3
        for i in range(stage.shape[0]):
            if i % 2 == 0:
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
            fid = NetCDFFile(domain.writer.filename, netcdf_mode_r)
            stage_file = fid.variables['stage']

            fid.close()

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
         b1 =  Dirichlet_boundary(dirichlet_values = num.array([0.0]))
         b2 =  Dirichlet_boundary(dirichlet_values = num.array([1.0]))
         b3 =  Dirichlet_boundary(dirichlet_values = num.array([2.0]))
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
         self.assertTrue(domain.boundary[(1, 0)]  == '1', msg)
         self.assertTrue(domain.boundary[(1, 2)]  == '2', msg)
         self.assertTrue(domain.boundary[(0, 1)]  == '3', msg)
         self.assertTrue(domain.boundary[(0, 0)]  == 'exterior', msg)
         msg = "test_pmesh2Domain Too many boundaries"
         self.assertTrue(len(domain.boundary)  == 4, msg)

         # FIXME change to use get_xllcorner
         msg = 'Bad geo-reference'
         self.assertTrue(domain.geo_reference.xllcorner  == 140.0, msg)

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
         self.assertTrue(domain.boundary[(1, 0)]  == '1', msg)
         self.assertTrue(domain.boundary[(1, 2)]  == '2', msg)
         self.assertTrue(domain.boundary[(0, 1)]  == '3', msg)
         self.assertTrue(domain.boundary[(0, 0)]  == 'exterior', msg)
         msg = "test_pmesh2Domain Too many boundaries"
         self.assertTrue(len(domain.boundary)  == 4, msg)

         # FIXME change to use get_xllcorner
         msg = 'Bad geo_reference'
         self.assertTrue(domain.geo_reference.xllcorner  == 140.0, msg)

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
        for data_point, attribute in zip(data_points_absolute, attributes):
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


    def test_fitting_in_hole(self):
        '''
            Make sure we can fit a mesh that has a hole in it.
            This is a regression test for ticket:234

        '''
        verbose = False

        from anuga.shallow_water.shallow_water_domain import Domain
        from anuga.pmesh.mesh_interface import create_mesh_from_regions
        from anuga.geospatial_data.geospatial_data import Geospatial_data


        # Get path where this test is run
        path = get_pathname_from_package('anuga.shallow_water')


        #----------------------------------------------------------------------
        # Create domain
        #--------------------------------------------------------------------
        W = 303400
        N = 6195800
        E = 308640
        S = 6193120
        border = 2000
        bounding_polygon = [[W, S], [E, S], [E, N], [W, N]]
        hole_polygon = [[W+border, S+border], [E-border, S+border], \
                        [E-border, N-border], [W+border, N-border]]

        meshname = 'offending_mesh.msh'
        create_mesh_from_regions(bounding_polygon,
                                 boundary_tags={'south': [0], 'east': [1],
                                                'north': [2], 'west': [3]},
                                 maximum_triangle_area=1000000,
                                 filename=meshname,
                                 interior_holes=[hole_polygon],
                                 use_cache=False,
                                 verbose=verbose)

        domain = Domain(meshname, use_cache=False, verbose=verbose)

        #----------------------------------------------------------------------
        # Fit data point inside hole to mesh
        #----------------------------------------------------------------------

        points_file = 'offending_point.pts'

        # Offending point
        G = Geospatial_data(data_points=[[(E+W) / 2, (N+S) / 2]],
                            attributes=[1])
        G.export_points_file(points_file)

        # fit data using the point within the hole.
        domain.set_quantity('elevation', filename=points_file,
                            use_cache=False, verbose=verbose, alpha=0.01)
        os.remove(meshname)
        os.remove(points_file)


    def test_fitting_example_that_crashed(self):
        """This unit test has been derived from a real world example
        (the Towradgi '98 rainstorm simulation).

        It shows a condition where fitting as called from set_quantity crashes
        when ANUGA mesh is reused. The test passes in the case where a new mesh
        is created.

        See ticket:314
        """

        verbose = False

        from anuga.shallow_water.shallow_water_domain import Domain
        from anuga.pmesh.mesh_interface import create_mesh_from_regions
        from anuga.geospatial_data.geospatial_data import Geospatial_data


        # Get path where thie data file are
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

        meshname = 'offending_mesh_1.msh'
        create_mesh_from_regions(bounding_polygon,
                                 boundary_tags={'south': [0], 'east': [1],
                                                'north': [2], 'west': [3]},
                                 maximum_triangle_area=1000000,
                                 interior_regions=interior_regions,
                                 filename=meshname,
                                 use_cache=False,
                                 verbose=verbose)

        domain = Domain(meshname, use_cache=False, verbose=verbose)

        #----------------------------------------------------------------------
        # Fit data point to mesh
        #----------------------------------------------------------------------

        points_file = 'offending_point_1.pts'

        # Offending point
        G = Geospatial_data(data_points=[[306953.344, 6194461.5]],
                            attributes=[1])
        G.export_points_file(points_file)

        try:
            domain.set_quantity('elevation', filename=points_file,
                                use_cache=False, verbose=verbose, alpha=0.01)
        except RuntimeError as e:
            msg = 'Test failed: %s' % str(e)
            raise Exception(msg)
            # clean up in case raise fails
            os.remove(meshname)
            os.remove(points_file)
        else:
            os.remove(meshname)
            os.remove(points_file)

        #finally:
            # Cleanup regardless
            #FIXME(Ole): Finally does not work like this in python2.4
            #FIXME(Ole): Reinstate this when Python2.4 is out of the way
            #FIXME(Ole): Python 2.6 apparently introduces something called 'with'
            #os.remove(meshname)
            #os.remove(points_file)


    def test_fitting_example_that_crashed_2(self):
        """test_fitting_example_that_crashed_2

        This unit test has been derived from a real world example
        (the JJKelly study, by Petar Milevski).

        It shows a condition where set_quantity crashes due to AtA
        not being built properly

        See ticket:314
        """

        verbose = False

        from anuga.shallow_water.shallow_water_domain import Domain
        from anuga.pmesh.mesh_interface import create_mesh_from_regions
        from anuga.geospatial_data import Geospatial_data

        # Get path where this test is run
        path = get_pathname_from_package('anuga.shallow_water')

        meshname = 'test_mesh_2.msh'

        W = 304180
        S = 6185270
        E = 307650
        N = 6189040
        maximum_triangle_area = 1000000

        bounding_polygon = [[W, S], [E, S], [E, N], [W, N]]

        create_mesh_from_regions(bounding_polygon,
                                 boundary_tags={'south': [0],
                                                'east': [1],
                                                'north': [2],
                                                'west': [3]},
                                 maximum_triangle_area=maximum_triangle_area,
                                 filename=meshname,
                                 use_cache=False,
                                 verbose=verbose)

        domain = Domain(meshname, use_cache=True, verbose=verbose)

        # Large test set revealed one problem
        points_file = os.path.join(path, 'tests', 'data', 'test_points_large.csv')

        domain.set_quantity('elevation', filename=points_file,
                                use_cache=False, verbose=verbose)
        try:
            domain.set_quantity('elevation', filename=points_file,
                                use_cache=False, verbose=verbose)
        except AssertionError as e:
            msg = 'Test failed: %s' % str(e)
            raise Exception(msg)
            # Cleanup in case this failed
            os.remove(meshname)

        # Small test set revealed another problem
        points_file = os.path.join(path, 'tests', 'data', 'test_points_small.csv')
        try:
            domain.set_quantity('elevation', filename=points_file,
                                use_cache=False, verbose=verbose)
        except AssertionError as e:
            msg = 'Test failed: %s' % str(e)
            raise Exception(msg)
            # Cleanup in case this failed
            os.remove(meshname)
        else:
            os.remove(meshname)


    def test_total_volume(self):
        """test_total_volume

        Test that total volume can be computed correctly
        """

        #----------------------------------------------------------------------
        # Import necessary modules
        #----------------------------------------------------------------------
        from anuga.abstract_2d_finite_volumes.mesh_factory \
                import rectangular_cross
        from anuga.shallow_water.shallow_water_domain import Domain

        #----------------------------------------------------------------------
        # Setup computational domain
        #----------------------------------------------------------------------

        length = 100.
        width  = 20.
        dx = dy = 5       # Resolution: of grid on both axes

        points, vertices, boundary = rectangular_cross(int(length / dx),
                                                       int(width / dy),
                                                       len1=length,
                                                       len2=width)
        domain = Domain(points, vertices, boundary)

        #----------------------------------------------------------------------
        # Simple flat bottom bathtub
        #----------------------------------------------------------------------

        d = 1.0
        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('stage', d)

        assert num.allclose(domain.compute_total_volume(), length*width*d)

        #----------------------------------------------------------------------
        # Slope
        #----------------------------------------------------------------------

        slope = 1.0/10          # RHS drops to -10m
        def topography(x, y):
            return -x * slope

        domain.set_quantity('elevation', topography)
        domain.set_quantity('stage', 0.0)       # Domain full

        V = domain.compute_total_volume()
        assert num.allclose(V, length*width*10 / 2)

        domain.set_quantity('stage', -5.0)      # Domain 'half' full

        # IMPORTANT: Adjust stage to match elevation
        domain.distribute_to_vertices_and_edges()

        V = domain.compute_total_volume()
        assert num.allclose(V, width*(length / 2)*5.0/2)


    def test_volumetric_balance_computation(self):
        """test_volumetric_balance_computation

        Test that total in and out flows are computed correctly
        in a steady state situation
        """

        # Set to True if volumetric output is sought
        verbose = False

        #----------------------------------------------------------------------
        # Import necessary modules
        #----------------------------------------------------------------------

        from anuga.abstract_2d_finite_volumes.mesh_factory \
                import rectangular_cross
        from anuga.shallow_water.shallow_water_domain import Domain
        from anuga.shallow_water.forcing import Inflow

        #----------------------------------------------------------------------
        # Setup computational domain
        #----------------------------------------------------------------------

        finaltime = 500.0
        length = 300.
        width  = 20.
        dx = dy = 5       # Resolution: of grid on both axes

        # Input parameters
        uh = 1.0
        vh = 0.0
        d = 1.0

        # 20 m^3/s in the x direction across entire domain
        ref_flow = uh*d*width

        points, vertices, boundary = rectangular_cross(int(length / dx),
                                                       int(width / dy),
                                                       len1=length,
                                                       len2=width)

        domain = Domain(points, vertices, boundary)
        domain.set_name('Inflow_flowline_test')              # Output name
        domain.set_low_froude(0)

        #----------------------------------------------------------------------
        # Setup initial conditions
        #----------------------------------------------------------------------

        domain.set_quantity('elevation', 0.0)  # Flat bed
        domain.set_quantity('friction', 0.0)   # Constant zero friction

        domain.set_quantity('stage', expression='elevation + %d' % d)

        #----------------------------------------------------------------------
        # Setup boundary conditions
        #----------------------------------------------------------------------

        Br = Reflective_boundary(domain)      # Solid reflective wall

        # Constant flow in and out of domain
        # Depth = 1m, uh=1 m/s, i.e. a flow of 20 m^3/s
        Bi = Dirichlet_boundary([d, uh, vh])
        Bo = Dirichlet_boundary([d, uh, vh])

        domain.set_boundary({'left': Bi, 'right': Bo, 'top': Br, 'bottom': Br})

        #----------------------------------------------------------------------
        # Evolve system through time
        #----------------------------------------------------------------------

        for t in domain.evolve(yieldstep=50.0, finaltime=finaltime):
            S = domain.volumetric_balance_statistics()
            if verbose :
                print(domain.timestepping_statistics())
                print(S)

            if t >= 400:
                # Steady state reached

                # Square on flowline at 200m
                q = domain.get_flow_through_cross_section([[200.0,  0.0],  [200.0, 20.0]])

                if verbose:
                    print(q, ref_flow)
                
                assert num.allclose(q, ref_flow)

        os.remove('Inflow_flowline_test.sww')

    def test_volume_conservation_inflow(self):
        """test_volume_conservation

        Test that total volume in domain is as expected, based on questions
        raised by Petar Milevski in May 2009.

        This test adds inflow at a known rate and verifies that the total
        terminal volume is as expected.

        """

        verbose = False


        #---------------------------------------------------------------------
        # Import necessary modules
        #---------------------------------------------------------------------
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
        from anuga.shallow_water.shallow_water_domain import Domain

        #----------------------------------------------------------------------
        # Setup computational domain
        #----------------------------------------------------------------------
        finaltime = 200.0

        length = 300.
        width  = 20.
        dx = dy = 5       # Resolution: of grid on both axes


        points, vertices, boundary = rectangular_cross(int(length / dx),
                                                       int(width / dy),
                                                       len1=length, len2=width)


        domain = Domain(points, vertices, boundary)
        domain.set_name('Inflow_volume_test')              # Output name


        #----------------------------------------------------------------------
        # Setup initial conditions
        #----------------------------------------------------------------------
        slope = 0.0
        def topography(x, y):
            z=-x * slope
            return z

        domain.set_quantity('elevation', topography) # Use function for elevation
        domain.set_quantity('friction', 0.0)         # Constant friction

        domain.set_quantity('stage',
                            expression='elevation')  # Dry initially


        #--------------------------------------------------------------
        # Setup Inflow
        #--------------------------------------------------------------

        # Fixed Flowrate onto Area
        fixed_inflow = Inflow(domain,
                              center=(10.0, 10.0),
                              radius=5.00,
                              rate=10.00)

        domain.forcing_terms.append(fixed_inflow)

        #----------------------------------------------------------------------
        # Setup boundary conditions
        #----------------------------------------------------------------------

        Br = Reflective_boundary(domain) # Solid reflective wall
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})


        #----------------------------------------------------------------------
        # Evolve system through time
        #----------------------------------------------------------------------
        ref_volume = 0.0
        ys = 10.0  # Yieldstep
        for t in domain.evolve(yieldstep=ys, finaltime=finaltime):

            # Check volume
            assert num.allclose(domain.compute_total_volume(), ref_volume)

            if verbose :
                print(domain.timestepping_statistics())
                print(domain.volumetric_balance_statistics())
                print('reference volume', ref_volume)


            # Update reference volume
            ref_volume += ys * fixed_inflow.rate


        os.remove('Inflow_volume_test.sww')



    def test_volume_conservation_rain(self):
        """test_volume_conservation

        Test that total volume in domain is as expected, based on questions
        raised by Petar Milevski in May 2009.

        This test adds rain at a known rate and verifies that the total
        terminal volume is as expected.

        """

        verbose = False


        #----------------------------------------------------------------------
        # Setup computational domain
        #----------------------------------------------------------------------
        finaltime = 200.0

        length = 300.
        width  = 20.
        dx = dy = 5       # Resolution: of grid on both axes


        points, vertices, boundary = rectangular_cross(int(length / dx),
                                                       int(width / dy),
                                                       len1=length, len2=width)


        domain = Domain(points, vertices, boundary)
        domain.set_name('Rain_volume_test')              # Output name


        #----------------------------------------------------------------------
        # Setup initial conditions
        #----------------------------------------------------------------------
        slope = 0.0
        def topography(x, y):
            z=-x * slope
            return z

        domain.set_quantity('elevation', topography) # Use function for elevation
        domain.set_quantity('friction', 0.0)         # Constant friction

        domain.set_quantity('stage',
                            expression='elevation')  # Dry initially


        #--------------------------------------------------------------
        # Setup rain
        #--------------------------------------------------------------

        # Fixed rain onto small circular area
        fixed_rain = Rainfall(domain,
                              center=(10.0, 10.0),
                              radius=5.00,
                              rate=10.00)   # 10 mm/s

        domain.forcing_terms.append(fixed_rain)

        #----------------------------------------------------------------------
        # Setup boundary conditions
        #----------------------------------------------------------------------

        Br = Reflective_boundary(domain) # Solid reflective wall
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})


        #----------------------------------------------------------------------
        # Evolve system through time
        #----------------------------------------------------------------------
        ref_volume = 0.0
        ys = 10.0  # Yieldstep
        for t in domain.evolve(yieldstep=ys, finaltime=finaltime):

            # Check volume
            V = domain.compute_total_volume()
            msg = 'V = %e, Ref = %e' % (V, ref_volume)
            assert num.allclose(V, ref_volume), msg

            if verbose :
                print(domain.timestepping_statistics())
                print(domain.volumetric_balance_statistics())
                print('reference volume', ref_volume)
                print(V)


            # Update reference volume.
            # FIXME: Note that rate has now been redefined
            # as m/s internally. This is a little confusing
            # when it was specfied as mm/s.

            delta_V = fixed_rain.rate*fixed_rain.exchange_area
            ref_volume += ys * delta_V

        os.remove('Rain_volume_test.sww')

    def Xtest_rain_conservation_and_runoff(self):
        """test_rain_conservation_and_runoff

        Test that total volume in domain is as expected, based on questions
        raised by Petar Milevski in May 2009.

        This test adds rain at a known rate and verifies that the total
        volume and outflows are as expected.

        """

        # FIXME (Ole): Does not work yet. Investigate boundary flows

        verbose = True #False


        #---------------------------------------------------------------------
        # Import necessary modules
        #---------------------------------------------------------------------
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
        from anuga.shallow_water.shallow_water_domain import Domain
        from anuga.shallow_water.boundaries import Reflective_boundary
        from anuga import Dirichlet_boundary
        from anuga.shallow_water.forcing import Rainfall
        from anuga.shallow_water.sww_interrogate import get_flow_through_cross_section

        #----------------------------------------------------------------------
        # Setup computational domain
        #----------------------------------------------------------------------
        finaltime = 500.0

        length = 300.
        width  = 20.
        dx = dy = 5       # Resolution: of grid on both axes


        points, vertices, boundary = rectangular_cross(int(length / dx),
                                                       int(width / dy),
                                                       len1=length, len2=width)


        domain = Domain(points, vertices, boundary)
        domain.set_name('Rain_volume_runoff_test')         # Output name


        #----------------------------------------------------------------------
        # Setup initial conditions
        #----------------------------------------------------------------------
        slope = 0.0
        def topography(x, y):
            z=-x * slope
            return z

        domain.set_quantity('elevation', topography) # Use function for elevation
        domain.set_quantity('friction', 0.0)         # Constant friction

        domain.set_quantity('stage',
                            expression='elevation')  # Dry initially


        #--------------------------------------------------------------
        # Setup rain
        #--------------------------------------------------------------

        # Fixed rain onto small circular area
        fixed_rain = Rainfall(domain,
                              center=(10.0, 10.0),
                              radius=5.00,
                              rate=10.00)   # 10 mm/s

        domain.forcing_terms.append(fixed_rain)

        #----------------------------------------------------------------------
        # Setup boundary conditions
        #----------------------------------------------------------------------

        Br = Reflective_boundary(domain) # Solid reflective wall
        Bt = Transmissive_stage_zero_momentum_boundary(domain)
        Bd = Dirichlet_boundary([-10, 0, 0])
        domain.set_boundary({'left': Bt, 'right': Bd, 'top': Bt, 'bottom': Bt})


        #----------------------------------------------------------------------
        # Evolve system through time
        #----------------------------------------------------------------------
        ref_volume = 0.0
        ys = 10.0  # Yieldstep
        for t in domain.evolve(yieldstep=ys, finaltime=finaltime):

            # Check volume
            V = domain.compute_total_volume()
            msg = 'V = %e, Ref = %e' % (V, ref_volume)
            #assert num.allclose(V, ref_volume) or V < ref_volume, msg

            if verbose:
                print(domain.timestepping_statistics())
                print(domain.volumetric_balance_statistics())
                print('reference volume', ref_volume)
                print(V)


            # Update reference volume.
            # FIXME: Note that rate has now been redefined
            # as m/s internally. This is a little confusing
            # when it was specfied as mm/s.

            delta_V = fixed_rain.rate*fixed_rain.exchange_area
            ref_volume += ys * delta_V

            # Compute outflow at right hand downstream boundary
            boundary_flows, inflow , outflow = domain.compute_boundary_flows()
            net_outflow = outflow - inflow

            outflow = boundary_flows['right']
            if verbose:
                print('Outflow', outflow)
                print('Net outflow', net_outflow)

            # Update reference volume
            ref_volume += ys * outflow


    def test_variable_elevation_de0(self):
        """test_variable_elevation

        This will test that elevagtion van be stored in sww files
        as a time dependent quantity.

        It will also chck that storage of other quantities
        can be controlled this way.
        """

        #---------------------------------------------------------------------
        # Import necessary modules
        #---------------------------------------------------------------------
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross

        #---------------------------------------------------------------------
        # Setup computational domain
        #---------------------------------------------------------------------
        length = 8.
        width = 6.
        dx = dy = 1    # Resolution: Length of subdivisions on both axes

        inc = 0.05 # Elevation increment

        points, vertices, boundary = rectangular_cross(int(length / dx),
                                                       int(width / dy),
                                                       len1=length,
                                                       len2=width)
        domain = Domain(points, vertices, boundary)
        domain.set_name('channel_variable_test')  # Output name
        domain.set_quantities_to_be_stored({'elevation': 2,
                                            'stage': 2})

        #pprint(domain.get_algorithm_parameters())
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
        domain.set_quantity('stage', 10.0)        # Dry initial condition

        #------------------------------------------------------------------
        # Setup boundary conditions
        #------------------------------------------------------------------
        Bi = Dirichlet_boundary([10.0, 0, 0])          # Inflow
        Br = Reflective_boundary(domain)              # Solid reflective wall
        Bo = Dirichlet_boundary([-5, 0, 0])           # Outflow

        domain.set_boundary({'left': Bi, 'right': Bo, 'top': Br, 'bottom': Br})

        #-------------------------------------------------------------------
        # Evolve system through time
        #-------------------------------------------------------------------

        for t in domain.evolve(yieldstep=1, finaltime=3.0):
            #print domain.timestepping_statistics()

            domain.add_quantity('elevation', pole_increment, location='centroids')


        # Check that quantities have been stored correctly
        sww_file = domain.get_name() + '.sww'
        fid = NetCDFFile(sww_file)

        stage = fid.variables['stage_c'][:]
        elevation = fid.variables['elevation_c'][:]
        fid.close()

        os.remove(sww_file)


        assert len(stage.shape) == 2
        assert len(elevation.shape) == 2

        M, N = stage.shape

        for i in range(M):
            # For each timestep
            assert num.allclose(max(elevation[i,:]), i * inc)

    def test_variable_elevation_1_5(self):
        """test_variable_elevation

        This will test that elevagtion van be stored in sww files
        as a time dependent quantity.

        It will also chck that storage of other quantities
        can be controlled this way.
        """

        #---------------------------------------------------------------------
        # Import necessary modules
        #---------------------------------------------------------------------
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross

        #---------------------------------------------------------------------
        # Setup computational domain
        #---------------------------------------------------------------------
        length = 8.
        width = 6.
        dx = dy = 1    # Resolution: Length of subdivisions on both axes

        inc = 0.05 # Elevation increment

        points, vertices, boundary = rectangular_cross(int(length / dx),
                                                       int(width / dy),
                                                       len1=length,
                                                       len2=width)
        domain = Domain(points, vertices, boundary)
        domain.set_flow_algorithm('1_5')
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
        domain.set_quantity('stage', 10.0)        # Dry initial condition

        #------------------------------------------------------------------
        # Setup boundary conditions
        #------------------------------------------------------------------
        Bi = Dirichlet_boundary([10.0, 0, 0])          # Inflow
        Br = Reflective_boundary(domain)              # Solid reflective wall
        Bo = Dirichlet_boundary([-5, 0, 0])           # Outflow

        domain.set_boundary({'left': Bi, 'right': Bo, 'top': Br, 'bottom': Br})

        #-------------------------------------------------------------------
        # Evolve system through time
        #-------------------------------------------------------------------

        for t in domain.evolve(yieldstep=1, finaltime=3.0):
            #print domain.timestepping_statistics()

            domain.add_quantity('elevation', pole_increment)


        # Check that quantities have been stored correctly
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


    def test_inflow_using_flowline(self):
        """test_inflow_using_flowline

        Test the ability of a flowline to match inflow above the flowline by
        creating constant inflow onto a circle at the head of a 20m
        wide by 300m long plane dipping at various slopes with a
        perpendicular flowline and gauge downstream of the inflow and
        a 45 degree flowlines at 200m downstream.

        A more substantial version of this test with finer resolution and
        including the depth calculation using Manning's equation is
        available under the validate_all suite in the directory
        anuga_validation/automated_validation_tests/flow_tests.
        """


        verbose = False

        #----------------------------------------------------------------------
        # Setup computational domain
        #----------------------------------------------------------------------
        number_of_inflows = 2 # Number of inflows on top of each other
        finaltime = 1000 #700.0 # If this is too short, steady state will not be achieved

        length = 250.
        width  = 20.
        dx = dy = 5                 # Resolution: of grid on both axes

        points, vertices, boundary = rectangular_cross(int(length / dx),
                                                       int(width / dy),
                                                       len1=length,
                                                       len2=width)
        for mannings_n in [0.1, 0.01]:
            # Loop over a range of roughnesses

            for slope in [1.0/300, 1.0/100]:
                # Loop over a range of bedslopes representing
                # sub to super critical flows


                domain = Domain(points, vertices, boundary)
                domain.set_name('inflow_flowline_test')     # Output name

                #--------------------------------------------------------------
                # Setup initial conditions
                #--------------------------------------------------------------

                def topography(x, y):
                    z = -x * slope
                    return z

                # Use function for elevation
                domain.set_quantity('elevation', topography)
                # Constant friction of conc surface
                domain.set_quantity('friction', mannings_n)
                # Dry initial condition
                domain.set_quantity('stage', expression='elevation')

                #--------------------------------------------------------------
                # Setup Inflow
                #--------------------------------------------------------------

                # Fixed Flowrate onto Area
                fixed_inflow = Inflow(domain,
                                      center=(10.0, 10.0),
                                      radius=5.00,
                                      rate=10.00)

                # Stack this flow
                for i in range(number_of_inflows):
                    domain.forcing_terms.append(fixed_inflow)

                ref_flow = fixed_inflow.rate*number_of_inflows

                # Compute normal depth on plane using Mannings equation
                # v=1/n*(r^2/3)*(s^0.5) or r=(Q*n/(s^0.5*W))^0.6
                normal_depth=(ref_flow*mannings_n / (slope**0.5*width))**0.6
                if verbose:
                    print()
                    print('Slope:', slope, 'Mannings n:', mannings_n)


                #--------------------------------------------------------------
                # Setup boundary conditions
                #--------------------------------------------------------------

                Br = Reflective_boundary(domain)

                # Define downstream boundary based on predicted depth
                def normal_depth_stage_downstream(t):
                    return (-slope*length) + normal_depth

                Bt = Transmissive_momentum_set_stage_boundary(domain=domain,
                                                              function=normal_depth_stage_downstream)




                domain.set_boundary({'left': Br,
                                     'right': Bt,
                                     'top': Br,
                                     'bottom': Br})



                #--------------------------------------------------------------
                # Evolve system through time
                #--------------------------------------------------------------

                for t in domain.evolve(yieldstep=10.0, finaltime=finaltime):
                    pass
                    if verbose :
                        print(domain.timestepping_statistics())
                        print(domain.volumetric_balance_statistics())


                #--------------------------------------------------------------
                # Compute flow thru flowlines ds of inflow
                #--------------------------------------------------------------

                # Square on flowline at 200m
                q = domain.get_flow_through_cross_section([[200.0,  0.0],
                                                           [200.0, 20.0]])
                if verbose:
                    print ('90 degree flowline: ANUGA = %f, Ref = %f'
                           % (q, ref_flow))

                msg = ('Predicted flow was %f, should have been %f'
                       % (q, ref_flow))
                assert num.allclose(q, ref_flow, rtol=5.0e-2), msg


                # 45 degree flowline at 200m
                q = domain.get_flow_through_cross_section([[200.0,  0.0],
                                                           [220.0, 20.0]])
                if verbose:
                    print ('45 degree flowline: ANUGA = %f, Ref = %f'
                           % (q, ref_flow))

                msg = ('Predicted flow was %f, should have been %f'
                       % (q, ref_flow))
                assert num.allclose(q, ref_flow, rtol=5.0e-2), msg

        os.remove('inflow_flowline_test.sww')


    def Xtest_inflow_boundary_using_flowline(self):
        """test_inflow_boundary_using_flowline
        Test the ability of a flowline to match inflow above the flowline by
        creating constant inflow into the boundary at the head of a 20m
        wide by 300m long plane dipping at various slopes with a
        perpendicular flowline and gauge downstream of the inflow and
        a 45 degree flowlines at 200m downstream


        """

        # FIXME (Ole): Work in progress

        verbose = False



        #----------------------------------------------------------------------
        # Setup computational domain
        #----------------------------------------------------------------------

        finaltime = 500 #700.0 # If this is too short, steady state will not be achieved

        length = 250.
        width  = 20.
        dx = dy = 5          # Resolution: of grid on both axes

        points, vertices, boundary = rectangular_cross(int(length / dx), int(width / dy),
                                                       len1=length, len2=width)

        for mannings_n in [0.1, 0.01]:
            # Loop over a range of roughnesses

            for slope in [1.0/300, 1.0/100]:
                # Loop over a range of bedslopes representing sub to super critical flows


                domain = Domain(points, vertices, boundary)
                domain.set_name('inflow_boundary_flowline_test')


                #-------------------------------------------------------------
                # Setup initial conditions
                #-------------------------------------------------------------

                def topography(x, y):
                    z=-x * slope
                    return z

                domain.set_quantity('elevation', topography)
                domain.set_quantity('friction', mannings_n)
                domain.set_quantity('stage',
                                    expression='elevation')



                #--------------------------------------------------------------
                # Setup boundary conditions
                #--------------------------------------------------------------



                ref_flow = 10.00

                # Compute normal depth on plane using Mannings equation
                # v=1/n*(r^2/3)*(s^0.5) or r=(Q*n/(s^0.5*W))^0.6
                normal_depth=(ref_flow*mannings_n / (slope**0.5*width))**0.6
                if verbose:
                    print()
                    print('Slope:', slope, 'Mannings n:', mannings_n)


                from anuga.shallow_water.boundaries import Inflow_boundary

                Bi = Inflow_boundary(domain, rate=ref_flow)

                Br = Reflective_boundary(domain)

                # Define downstream boundary based on predicted depth
                def normal_depth_stage_downstream(t):
                    return (-slope*length) + normal_depth

                Bt = Transmissive_momentum_set_stage_boundary(domain=domain,
                                                              function=normal_depth_stage_downstream)




                domain.set_boundary({'left': Bi,
                                     'right': Bt,
                                     'top': Br,
                                     'bottom': Br})



                #--------------------------------------------------------------
                # Evolve system through time
                #--------------------------------------------------------------


                for t in domain.evolve(yieldstep=100.0, finaltime=finaltime):
                    pass
                    #if verbose :
                    #    print domain.timestepping_statistics()
                    #    print domain.volumetric_balance_statistics()



                #--------------------------------------------------------------
                # Compute flow thru flowlines ds of inflow
                #--------------------------------------------------------------

                # Square on flowline at 200m
                q=domain.get_flow_through_cross_section([[200.0,0.0],[200.0,20.0]])
                msg = 'Predicted flow was %f, should have been %f' % (q, ref_flow)
                if verbose:
                    print('90 degree flowline: ANUGA = %f, Ref = %f' % (q, ref_flow))
                assert num.allclose(q, ref_flow, rtol=1.0e-2), msg


                # 45 degree flowline at 200m
                q=domain.get_flow_through_cross_section([[200.0,0.0],[220.0,20.0]])
                msg = 'Predicted flow was %f, should have been %f' % (q, ref_flow)
                if verbose:
                    print('45 degree flowline: ANUGA = %f, Ref = %f' % (q, ref_flow))

                assert num.allclose(q, ref_flow, rtol=1.0e-2), msg



    def Xtest_friction_dependent_flow_using_flowline(self):
        """test_friction_dependent_flow_using_flowline

        Test the internal flow (using flowline) as a function of
        different values of Mannings n and different slopes.

        Flow is applied in the form of boundary conditions with fixed momentum.
        """

        verbose = True

        #----------------------------------------------------------------------
        # Import necessary modules
        #----------------------------------------------------------------------

        from anuga.abstract_2d_finite_volumes.mesh_factory \
                import rectangular_cross
        from anuga.shallow_water.shallow_water_domain import Domain
        from anuga.shallow_water.boundaries import Reflective_boundary
        from anuga import Dirichlet_boundary
        from anuga.shallow_water.forcing import Inflow

        from anuga.abstract_2d_finite_volumes.util \
                import sww2csv_gauges, csv2timeseries_graphs


        #----------------------------------------------------------------------
        # Setup computational domain
        #----------------------------------------------------------------------

        finaltime = 1000.0

        length = 300.
        width  = 20.
        dx = dy = 5       # Resolution: of grid on both axes

        # Input parameters
        uh = 1.0
        vh = 0.0
        d = 1.0

        ref_flow = uh*d*width # 20 m^3/s in the x direction across entire domain

        points, vertices, boundary = rectangular_cross(int(length / dx),
                                                       int(width / dy),
                                                       len1=length,
                                                       len2=width)

        for mannings_n in [0.035]:          #[0.0, 0.012, 0.035]:
            for slope in [1.0/300]:         #[0.0, 1.0/300, 1.0/150]:
                # Loop over a range of bedslopes representing
                # sub to super critical flows
                if verbose:
                    print()
                    print('Slope:', slope, 'Mannings n:', mannings_n)
                domain = Domain(points, vertices, boundary)
                domain.set_name('Inflow_flowline_test')     # Output name

                #--------------------------------------------------------------
                # Setup initial conditions
                #--------------------------------------------------------------

                def topography(x, y):
                    z = -x * slope
                    return z

                # Use function for elevation
                domain.set_quantity('elevation', topography)
                # Constant friction
                domain.set_quantity('friction', mannings_n)

                #domain.set_quantity('stage', expression='elevation')

                # Set initial flow as depth=1m, uh=1.0 m/s, vh = 0.0
                # making it 20 m^3/s across entire domain
                domain.set_quantity('stage', expression='elevation + %f' % d)
                domain.set_quantity('xmomentum', uh)
                domain.set_quantity('ymomentum', vh)

                #--------------------------------------------------------------
                # Setup boundary conditions
                #--------------------------------------------------------------

                Br = Reflective_boundary(domain)      # Solid reflective wall

                # Constant flow in and out of domain
                # Depth = 1m, uh=1 m/s, i.e. a flow of 20 m^3/s
                # across boundaries
                Bi = Dirichlet_boundary([d, uh, vh])
                Bo = Dirichlet_boundary([-length*slope+d, uh, vh])
                #Bo = Dirichlet_boundary([-100, 0, 0])

                domain.set_boundary({'left': Bi, 'right': Bo,
                                     'top': Br,  'bottom': Br})

                #--------------------------------------------------------------
                # Evolve system through time
                #--------------------------------------------------------------

                for t in domain.evolve(yieldstep=100.0, finaltime=finaltime):
                    if verbose :
                        print(domain.timestepping_statistics())
                        print(domain.volumetric_balance_statistics())

                # 90 degree flowline at 200m
                q = domain.get_flow_through_cross_section([[200.0,  0.0],
                                                           [200.0, 20.0]])
                msg = ('Predicted flow was %f, should have been %f'
                       % (q, ref_flow))
                if verbose:
                    print ('90 degree flowline: ANUGA = %f, Ref = %f'
                           % (q, ref_flow))

                # 45 degree flowline at 200m
                q = domain.get_flow_through_cross_section([[200.0,  0.0],
                                                           [220.0, 20.0]])
                msg = ('Predicted flow was %f, should have been %f'
                       % (q, ref_flow))
                if verbose:
                    print ('45 degree flowline: ANUGA = %f, Ref = %f'
                           % (q, ref_flow))

                # Stage recorder (gauge) in middle of plane at 200m
                x = 200.0
                y = 10.00
                w = domain.get_quantity('stage').\
                        get_values(interpolation_points=[[x, y]])[0]
                z = domain.get_quantity('elevation').\
                        get_values(interpolation_points=[[x, y]])[0]
                domain_depth = w-z

                xmom = domain.get_quantity('xmomentum').\
                        get_values(interpolation_points=[[x, y]])[0]
                ymom = domain.get_quantity('ymomentum').\
                        get_values(interpolation_points=[[x, y]])[0]
                if verbose:
                    print(('At interpolation point (h, uh, vh): ',
                           domain_depth, xmom, ymom))
                    print('uh * d * width = ', xmom*domain_depth*width)

                if slope > 0.0:
                    # Compute normal depth at gauge location using Manning eqn
                    # v=1/n*(r^2/3)*(s^0.5) or r=(Q*n/(s^0.5*W))^0.6
                    normal_depth = (ref_flow*mannings_n / (slope**0.5*width))**0.6
                    if verbose:
                        print ('Depth: ANUGA = %f, Mannings = %f'
                               % (domain_depth, normal_depth))

        os.remove('Inflow_flowline_test.sww')


    def test_track_speeds(self):
        """
        get values based on triangle lists.
        """
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.timestepping_statistics(track_speeds=True)



    def test_tag_region_tags(self):
        """
        get values based on triangle lists.
        """
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.build_tagged_elements_dictionary({'bottom':[0,1],
                                                 'top':[4,5],
                                                 'all':[0,1,2,3,4,5]})


        #Set friction
        manning = 0.07
        domain.set_quantity('friction', manning)

        domain.set_tag_region([set_bottom_friction, set_top_friction])
        assert num.allclose(domain.quantities['friction'].get_values(),\
                            [[ 0.09,  0.09,  0.09],
                             [ 0.09,  0.09,  0.09],
                             [ 0.07,  0.07,  0.07],
                             [ 0.07,  0.07,  0.07],
                             [ 1.0,  1.0,  1.0],
                             [ 1.0,  1.0,  1.0]])

        domain.set_tag_region([set_all_friction])
        assert num.allclose(domain.quantities['friction'].get_values(),
                            [[ 10.09, 10.09, 10.09],
                             [ 10.09, 10.09, 10.09],
                             [ 10.07, 10.07, 10.07],
                             [ 10.07, 10.07, 10.07],
                             [ 11.0,  11.0,  11.0],
                             [ 11.0,  11.0,  11.0]])


    def test_region_tags2(self):
        """
        get values based on triangle lists.
        """
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.build_tagged_elements_dictionary({'bottom':[0,1],
                                                 'top':[4,5],
                                                 'all':[0,1,2,3,4,5]})


        #Set friction
        manning = 0.07
        domain.set_quantity('friction', manning)

        domain.set_tag_region('top', 'friction', 1.0)
        domain.set_tag_region('bottom', 'friction', 0.09)

        msg = ("domain.quantities['friction'].get_values()=\n%s\n"
               'should equal\n'
               '[[ 0.09,  0.09,  0.09],\n'
               ' [ 0.09,  0.09,  0.09],\n'
               ' [ 0.07,  0.07,  0.07],\n'
               ' [ 0.07,  0.07,  0.07],\n'
               ' [ 1.0,  1.0,  1.0],\n'
               ' [ 1.0,  1.0,  1.0]]'
               % str(domain.quantities['friction'].get_values()))
        assert num.allclose(domain.quantities['friction'].get_values(),
                            [[ 0.09,  0.09,  0.09],
                             [ 0.09,  0.09,  0.09],
                             [ 0.07,  0.07,  0.07],
                             [ 0.07,  0.07,  0.07],
                             [ 1.0,  1.0,  1.0],
                             [ 1.0,  1.0,  1.0]]), msg

        domain.set_tag_region([set_bottom_friction, set_top_friction])
        assert num.allclose(domain.quantities['friction'].get_values(),
                            [[ 0.09,  0.09,  0.09],
                             [ 0.09,  0.09,  0.09],
                             [ 0.07,  0.07,  0.07],
                             [ 0.07,  0.07,  0.07],
                             [ 1.0,  1.0,  1.0],
                             [ 1.0,  1.0,  1.0]])

        domain.set_tag_region([set_all_friction])
        assert num.allclose(domain.quantities['friction'].get_values(),
                            [[ 10.09, 10.09, 10.09],
                             [ 10.09, 10.09, 10.09],
                             [ 10.07, 10.07, 10.07],
                             [ 10.07, 10.07, 10.07],
                             [ 11.0,  11.0,  11.0],
                             [ 11.0,  11.0,  11.0]])




    def test_vertex_values_no_smoothing(self):

        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular
        from anuga.utilities.numerical_tools import mean


        #Create basic mesh
        points, vertices, boundary = rectangular(2, 2)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.default_order=2
        domain.reduction = mean


        #Set some field values
        domain.set_quantity('elevation', lambda x,y: x)
        domain.set_quantity('friction', 0.03)


        ######################
        #Initial condition - with jumps

        bed = domain.quantities['elevation'].vertex_values
        stage = num.zeros(bed.shape, float)

        h = 0.03
        for i in range(stage.shape[0]):
            if i % 2 == 0:
                stage[i,:] = bed[i,:] + h
            else:
                stage[i,:] = bed[i,:]

        domain.set_quantity('stage', stage)

        #Get stage
        stage = domain.quantities['stage']
        A, V = stage.get_vertex_values(xy=False, smooth=False)
        Q = stage.vertex_values.flatten()

        for k in range(8):
            assert num.allclose(A[k], Q[k])


        for k in range(8):
            assert V[k, 0] == 3*k
            assert V[k, 1] == 3*k+1
            assert V[k, 2] == 3*k+2



        X, Y, A1, V1 = stage.get_vertex_values(xy=True, smooth=False)


        assert num.allclose(A, A1)
        assert num.allclose(V, V1)

        #Check XY
        assert num.allclose(X[1], 0.5)
        assert num.allclose(Y[1], 0.5)
        assert num.allclose(X[4], 0.0)
        assert num.allclose(Y[4], 0.0)
        assert num.allclose(X[12], 1.0)
        assert num.allclose(Y[12], 0.0)



    #Test smoothing
    def test_smoothing_1_5(self):

        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular
        from anuga.utilities.numerical_tools import mean

        #Create basic mesh
        points, vertices, boundary = rectangular(2, 2)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.set_flow_algorithm('1_5')
        domain.reduction = mean


        #Set some field values
        domain.set_quantity('elevation', lambda x,y: x)
        domain.set_quantity('friction', 0.03)


        ######################
        # Boundary conditions
        B = Transmissive_boundary(domain)
        domain.set_boundary( {'left': B, 'right': B, 'top': B, 'bottom': B})


        ######################
        # Initial condition - with jumps

        bed = domain.quantities['elevation'].vertex_values
        stage = num.zeros(bed.shape, float)

        h = 0.03
        for i in range(stage.shape[0]):
            if i % 2 == 0:
                stage[i,:] = bed[i,:] + h
            else:
                stage[i,:] = bed[i,:]

        domain.set_quantity('stage', stage)

        stage = domain.quantities['stage']

        # Get smoothed stage
        A, V = stage.get_vertex_values(xy=False, smooth=True)
        Q = stage.vertex_values


        assert A.shape[0] == 9
        assert V.shape[0] == 8
        assert V.shape[1] == 3

        # First four points
        assert num.allclose(A[0], (Q[0,2] + Q[1,1]) / 2)
        assert num.allclose(A[1], (Q[1,0] + Q[3,1] + Q[2,2]) / 3)
        assert num.allclose(A[2], Q[3,0])
        assert num.allclose(A[3], (Q[0,0] + Q[5,1] + Q[4,2]) / 3)

        # Center point
        assert num.allclose(A[4], (Q[0,1] + Q[1,2] + Q[2,0] +
                                   Q[5,0] + Q[6,2] + Q[7,1]) / 6)

        # Check V
        assert num.allclose(V[0,:], [3,4,0])
        assert num.allclose(V[1,:], [1,0,4])
        assert num.allclose(V[2,:], [4,5,1])
        assert num.allclose(V[3,:], [2,1,5])
        assert num.allclose(V[4,:], [6,7,3])
        assert num.allclose(V[5,:], [4,3,7])
        assert num.allclose(V[6,:], [7,8,4])
        assert num.allclose(V[7,:], [5,4,8])

        # Get smoothed stage with XY
        X, Y, A1, V1 = stage.get_vertex_values(xy=True, smooth=True)

        assert num.allclose(A, A1)
        assert num.allclose(V, V1)

        # Check XY
        assert num.allclose(X[4], 0.5)
        assert num.allclose(Y[4], 0.5)

        assert num.allclose(X[7], 1.0)
        assert num.allclose(Y[7], 0.5)



    def test_smoothing_de0(self):

        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular
        from anuga.utilities.numerical_tools import mean

        # Create basic mesh
        points, vertices, boundary = rectangular(2, 2)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.set_flow_algorithm('DE0')
        domain.reduction = mean


        # Set some field values
        domain.set_quantity('elevation', lambda x,y: x)
        domain.set_quantity('friction', 0.03)


        ######################
        # Boundary conditions
        B = Transmissive_boundary(domain)
        domain.set_boundary( {'left': B, 'right': B, 'top': B, 'bottom': B})


        ######################
        # Initial condition - with jumps

        bed = domain.quantities['elevation'].vertex_values
        stage = num.zeros(bed.shape, float)

        h = 0.03
        for i in range(stage.shape[0]):
            if i % 2 == 0:
                stage[i,:] = bed[i,:] + h
            else:
                stage[i,:] = bed[i,:]

        domain.set_quantity('stage', stage)

        stage = domain.quantities['stage']

        # Get smoothed stage
        A, V = stage.get_vertex_values(xy=False, smooth=True)
        Q = stage.centroid_values


        assert A.shape[0] == 9
        assert V.shape[0] == 8
        assert V.shape[1] == 3

        # First four points
        assert num.allclose(A[0], (Q[0] + Q[1]) / 2)
        assert num.allclose(A[1], (Q[1] + Q[3] + Q[2]) / 3)
        assert num.allclose(A[2], Q[3])
        assert num.allclose(A[3], (Q[0] + Q[5] + Q[4]) / 3)

        # Center point
        assert num.allclose(A[4], (Q[0] + Q[1] + Q[2] +
                                   Q[5] + Q[6] + Q[7]) / 6)


        # Check V
        assert num.allclose(V[0,:], [3,4,0])
        assert num.allclose(V[1,:], [1,0,4])
        assert num.allclose(V[2,:], [4,5,1])
        assert num.allclose(V[3,:], [2,1,5])
        assert num.allclose(V[4,:], [6,7,3])
        assert num.allclose(V[5,:], [4,3,7])
        assert num.allclose(V[6,:], [7,8,4])
        assert num.allclose(V[7,:], [5,4,8])

        # Get smoothed stage with XY
        X, Y, A1, V1 = stage.get_vertex_values(xy=True, smooth=True)

        assert num.allclose(A, A1)
        assert num.allclose(V, V1)

        # Check XY
        assert num.allclose(X[4], 0.5)
        assert num.allclose(Y[4], 0.5)

        assert num.allclose(X[7], 1.0)
        assert num.allclose(Y[7], 0.5)


    #Test calculating velocities and back to momenta
    # useful for kinematic viscosity calc
    def test_update_centroids_of_velocities_and_height(self):

        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular
        from anuga.utilities.numerical_tools import mean

        #Create basic mesh
        points, vertices, boundary = rectangular(2, 2)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.default_order=2
        domain.reduction = mean


        #Set some field values
        domain.set_quantity('elevation', lambda x,y: y)
        domain.set_quantity('friction', 0.03)


        W  = domain.quantities['stage']
        UH = domain.quantities['xmomentum']
        VH = domain.quantities['ymomentum']
        H  = domain.quantities['height']
        Z  = domain.quantities['elevation']
        U  = domain.quantities['xvelocity']
        V  = domain.quantities['yvelocity']
        X  = domain.quantities['x']
        Y  = domain.quantities['y']

        Wc  = W.centroid_values
        UHc = UH.centroid_values
        VHc = VH.centroid_values
        Hc  = H.centroid_values
        Zc  = Z.centroid_values
        Uc  = U.centroid_values
        Vc  = V.centroid_values
        Xc  = X.centroid_values
        Yc  = Y.centroid_values



        ######################
        # Boundary conditions
        #B = Transmissive_boundary(domain)
        B = Reflective_boundary(domain)
        domain.set_boundary( {'left': B, 'right': B, 'top': B, 'bottom': B})


        domain.set_quantity('stage',expression='elevation - 2*x')
        domain.set_quantity('xmomentum', expression='2*x+3*y')
        domain.set_quantity('ymomentum', expression='5*x+7*y')


        assert num.allclose(Wc, Zc-2*Xc)
        assert num.allclose(UHc, 2*Xc+3*Yc)
        assert num.allclose(VHc, 5*Xc+7*Yc)


#        try:
#            domain.update_centroids_of_velocities_and_height()
#        except AssertionError:
#            pass
#        else:
#            raise Exception('should have caught H<0 error')

        domain.set_quantity('stage',expression='elevation + 2*x')
        assert num.allclose(Wc, Zc+2*Xc)

        domain.update_boundary()
        domain.update_centroids_of_velocities_and_height()


        assert num.allclose(Uc, UHc / Hc)
        assert num.allclose(Vc, VHc / Hc)
        assert num.allclose(Hc, Wc - Zc)

        # Lets change the U and V and change back to UH and VH

        domain.set_quantity('xvelocity', expression='2*x+3*y')
        domain.set_quantity('yvelocity', expression='5*x+7*y')

        domain.set_quantity('height', expression='x+y')

        assert num.allclose(Uc, 2*Xc+3*Yc)
        assert num.allclose(Vc, 5*Xc+7*Yc)
        assert num.allclose(Hc, Xc + Yc)

        domain.update_centroids_of_momentum_from_velocity()

        assert num.allclose(UHc, (2*Xc+3*Yc)*(Xc+Yc))
        assert num.allclose(VHc, (5*Xc+7*Yc)*(Xc+Yc))

    def test_set_quantity_from_file(self):
        '''test the new set_values for the set_quantity procedures. The results of setting quantity values by using set_quantities and set_values for pts, asc and dem files are tested here. They should be equal.'''

        # settup
        x0 = 0.0
        y0 = 0.0

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]

        # bac, bce, ecf, dbe
        elements = [ [1, 0, 2], [1, 2, 4], [4, 2, 5], [3, 1, 4] ]

        # absolute going in ..
        mesh4 = Domain(points, elements, geo_reference=Geo_reference(56, 0, 0))
        mesh4.check_integrity()
        quantity = Quantity(mesh4)

        # Get (enough) datapoints (relative to georef)
        data_points_rel = [[ 0.66666667, 0.66666667],
                           [ 1.33333333, 1.33333333],
                           [ 2.66666667, 0.66666667],
                           [ 0.66666667, 2.66666667],
                           [ 0.0, 1.0],
                           [ 0.0, 3.0],
                           [ 1.0, 0.0],
                           [ 1.0, 1.0],
                           [ 1.0, 2.0],
                           [ 1.0, 3.0],
                           [ 2.0, 1.0],
                           [ 3.0, 0.0],
                           [ 3.0, 1.0]]

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
        for data_point, attribute in zip(data_points_absolute, attributes):
            row = (str(data_point[0]) + ',' +
                   str(data_point[1]) + ',' +
                   str(attribute))
            file.write(row + "\n")
        file.close()

        # exact answer (vertex locations)
        answer_vertex_values = linear_function(quantity.domain.get_vertex_coordinates())

        # Check that values can be set from pts file
        # using set_values directly
        quantity.set_values(filename=ptsfile, alpha=0)
        assert num.allclose(quantity.vertex_values.flat, answer_vertex_values)

        # using set_quantity with quantity name stage
        mesh4.set_quantity(name='stage', filename=ptsfile, alpha=0)
        mesh4_stage = mesh4.get_quantity('stage')
        assert num.allclose(mesh4_stage.vertex_values.flat, answer_vertex_values)

        # check set quantity from asc file

        """ Format of asc file
        ncols         11
        nrows         12
        xllcorner     240000
        yllcorner     7620000
        cellsize      6000
        NODATA_value  -9999
        """
        ncols = 11  # Nx
        nrows = 12  # Ny
        xllcorner = x0
        yllcorner = y0
        cellsize = 1.0
        NODATA_value = -9999

        # Create .asc file
        txt_file = 'test_asc.asc'
        datafile = open(txt_file, "w")
        datafile.write('ncols ' + str(ncols) + "\n")
        datafile.write('nrows ' + str(nrows) + "\n")
        datafile.write('xllcorner ' + str(xllcorner) + "\n")
        datafile.write('yllcorner ' + str(yllcorner) + "\n")
        datafile.write('cellsize ' + str(cellsize) + "\n")
        datafile.write('NODATA_value ' + str(NODATA_value) + "\n")

        x = num.linspace(xllcorner, xllcorner + (ncols - 1) * cellsize, ncols)
        y = num.linspace(yllcorner, yllcorner + (nrows - 1) * cellsize, nrows)
        points = axes2points(x, y)
        datavalues = linear_function(points)

        datavalues = datavalues.reshape(nrows, ncols)

        for row in datavalues:
            datafile.write(" ".join(str(elem) for elem in row) + "\n")
        datafile.close()

        # check set_values from asc file
        quantity.set_values(0.0)
        quantity.set_values(filename=txt_file,
                            location='vertices',
                            indices=None,
                            verbose=False)
        assert num.allclose(quantity.vertex_values.flat, answer_vertex_values)

        quantity.set_values(0.0)
        quantity.set_values(filename=txt_file,
                            location='centroids',
                            indices=None,
                            verbose=False)

        # exact answer for centroid locations
        answer_centroid_values = [ 1.33333333, 2.66666667, 3.33333333, 3.33333333]
        assert num.allclose(quantity.centroid_values, answer_centroid_values)

        # check set_quantity from asc file
        mesh4.set_quantity(name='stage', filename=txt_file,
                            location='vertices', indices=None, verbose=False)
        mesh4_stage = mesh4.get_quantity('stage')
        assert num.allclose(mesh4_stage.vertex_values.flat, answer_vertex_values)
        # reset mesh4 stage values
        mesh4.set_quantity(name='stage', numeric=0.0)
        mesh4.set_quantity(name='stage', filename=txt_file,
                            location='centroids', indices=None, verbose=False)
        assert num.allclose(mesh4_stage.centroid_values, answer_centroid_values)

        # check set quantity values from dem file
        from anuga.file_conversion.asc2dem import asc2dem
        # use the same reference solution used above for testing
        # convert test_asc.asc file to .dem file
        txt_file_prj = 'test_asc.prj'
        fid = open(txt_file_prj, 'w')
        fid.write("""Projection UTM
        Zone 56
        Datum WGS84
        Zunits NO
        Units METERS
        Spheroid WGS84
        Xshift 0.0000000000
        Yshift 10000000.0000000000
        Parameters
        """)
        fid.close()

        txt_file_dem = 'test_asc.dem'
        asc2dem(name_in=txt_file, name_out='test_asc',
                use_cache=False, verbose=False)

        # check set_values from dem file
        quantity.set_values(0.0)
        quantity.set_values(filename=txt_file_dem,
                            location='vertices',
                            indices=None,
                            verbose=False)
        assert num.allclose(quantity.vertex_values.flat, answer_vertex_values)

        quantity.set_values(0.0)
        quantity.set_values(filename=txt_file_dem,
                            location='centroids',
                            indices=None,
                            verbose=False)
        assert num.allclose(quantity.centroid_values, answer_centroid_values)

        # check set_quantity from dem file
        mesh4.set_quantity(name='stage', filename=txt_file_dem,
                            location='vertices', indices=None, verbose=False)
        mesh4_stage = mesh4.get_quantity('stage')
        assert num.allclose(mesh4_stage.vertex_values.flat, answer_vertex_values)
        # reset mesh4 stage values
        mesh4.set_quantity(name='stage', numeric=0.0)
        mesh4.set_quantity(name='stage', filename=txt_file_dem,
                            location='centroids')
        mesh4_stage = mesh4.get_quantity('stage')
        assert num.allclose(mesh4_stage.centroid_values, answer_centroid_values)

        # Cleanup
        import os
        try:
            os.remove(ptsfile)
            os.remove(txt_file)
            os.remove(txt_file_prj)
            os.remove(txt_file_dem)
        except:
            pass

    def test_that_mesh_methods_exist(self):
        """test_that_mesh_methods_exist

        Test that relavent mesh methods are made available in
        domain through composition
        """

        # Create basic mesh
        points, vertices, boundary = rectangular(1, 3)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)


        domain.get_centroid_coordinates()
        domain.get_radii()
        domain.get_areas()
        domain.get_area()
        domain.get_vertex_coordinates()
        domain.get_triangles()
        domain.get_nodes()
        domain.get_number_of_nodes()
        domain.get_normal(0,0)
        domain.get_triangle_containing_point([0.4,0.5])
        domain.get_intersecting_segments([[0.0, 0.0], [0.0, 1.0]])
        domain.get_disconnected_triangles()
        domain.get_boundary_tags()
        domain.get_boundary_polygon()
        #domain.get_number_of_triangles_per_node()
        domain.get_triangles_and_vertices_per_node()
        domain.get_interpolation_object()
        domain.get_tagged_elements()
        domain.get_lone_vertices()
        domain.get_unique_vertices()
        g = domain.get_georeference()
        domain.set_georeference(g)
        domain.build_tagged_elements_dictionary()
        domain.statistics()
        domain.get_extent()



#################################################################################

if __name__ == "__main__":
    #suite = unittest.makeSuite(Test_Shallow_Water, 'test_balance_deep_and_shallow_Froude')
    suite = unittest.makeSuite(Test_Shallow_Water, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
