#!/usr/bin/env python
import pytest
import unittest
import os
import tempfile
from math import pi, sqrt

from anuga.file.netcdf import NetCDFFile

from anuga.config import g, epsilon
from anuga.config import netcdf_mode_r
from anuga.utilities.numerical_tools import mean
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.geospatial_data.geospatial_data import Geospatial_data
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross, \
                                            rectangular
from anuga.abstract_2d_finite_volumes.quantity import Quantity
from anuga.shallow_water.forcing import Cross_section
from anuga.operators.rate_operators import Rate_operator

from anuga.utilities.system_tools import get_pathname_from_package

from anuga.shallow_water.shallow_water_domain import Domain

# boundary functions
from anuga.shallow_water.boundaries import Reflective_boundary, \
            Field_boundary, Transmissive_momentum_set_stage_boundary, \
            Transmissive_stage_zero_momentum_boundary
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions \
     import Transmissive_boundary, Dirichlet_boundary, Time_boundary

import numpy as num

# Get gateway to C implementation of flux function for direct testing
from anuga.shallow_water.sw_domain_openmp_ext import flux_function_central as flux_function
from anuga.shallow_water.sw_domain_openmp_ext import rotate


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

class Weir:
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
                z[i] = num.sqrt(((x[i]-x0))**2 + (y[i]-y0)**2)-1.0

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


class Weir_simple:
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

        assert num.all(domain.get_conserved_quantities(0, edge=1) == 0.)


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
        assert parameters['use_centroid_velocities'] == domain.use_centroid_velocities
        assert parameters['use_sloped_mannings']     == domain.use_sloped_mannings
        assert parameters['compute_fluxes_method']   == domain.get_compute_fluxes_method()
        assert parameters['flow_algorithm']          == domain.get_flow_algorithm()

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

        fid = NetCDFFile(os.path.join(domain.get_datadir(), domain.get_name() + '.sww'))
        sww_institution = fid.institution
        fid.close()

        assert institution == sww_institution

        os.remove(os.path.join(domain.get_datadir(), domain.get_name() + '.sww'))


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
        assert num.all(z < initial_runup_height)

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
        assert num.all(wet_elevations < initial_runup_height)
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

        assert num.all(z < initial_runup_height)

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


        assert num.all(z < final_runup_height)

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
        assert num.all(wet_elevations < final_runup_height)
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

        import time
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

        import time
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

        from anuga.pmesh.mesh_interface import create_pmesh_from_regions
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

        G0_1 = [-0.19166666666666665,
                -0.19166666666666665,
                -0.19166666666666665,
                -0.19166666666666665,
                -0.19166666666666665,
                -0.17381190304128577,
                -0.16132280165665422,
                -0.15298925501851754,
                -0.15659583237258828,
                -0.15621776283704272,
                -0.15986460119926058,
                -0.17203963475377068,
                -0.18810864424882115,
                -0.1979117107239742,
                -0.1991754654915448,
                -0.19897683957117226,
                -0.19945206336142793,
                -0.19959500424874452,
                -0.1995340832395464,
                -0.1997215683621538,
                -0.19621436205094506,
                -0.19142467090766674,
                -0.19150252620861344,
                -0.19155438631994037,
                -0.19154761263501988,
                -0.1915842196209439,
                -0.19160468252140941,
                -0.19160980364199998,
                -0.19161914951094586,
                -0.19162691799721024,
                -0.19163232811516856,
                -0.19163701089018073,
                -0.19164109241305016,
                -0.1916444925786887,
                -0.1916473494518105,
                -0.19164977043095288,
                -0.19165182888531462,
                -0.191653584636309,
                -0.19165508793432617,
                -0.19165637965808152,
                -0.19165749374570998,
                -0.19165845769444367,
                -0.19165929440250323,
                -0.19166002071132715,
                -0.19166065870797933,
                -0.19166122408635325,
                -0.1916617220060678,
                -0.19166216670095917,
                -0.19166254866879748,
                -0.19166290099461689,
                -0.19166321355271052]

        G1_1 = [-0.2916666666666667,
                -0.2916666666666667,
                -0.2916666666666667,
                -0.2621244946477121,
                -0.23779844075687265,
                -0.22348541313746256,
                -0.21191178172453382,
                -0.20069915844821729,
                -0.19210278772405676,
                -0.18445068377648965,
                -0.17997481753201655,
                -0.17713487021319155,
                -0.17854582614019068,
                -0.18509953421213768,
                -0.19577227155174948,
                -0.20277790738092685,
                -0.20605687532003003,
                -0.20679289270229415,
                -0.20640002223540913,
                -0.20541003998870688,
                -0.20406235954847418,
                -0.2023080150140285,
                -0.20078527032438662,
                -0.20014862432859945,
                -0.19933334291426985,
                -0.19875104646941422,
                -0.19839233749317733,
                -0.19841347871193735,
                -0.19875959510049476,
                -0.19929247426493973,
                -0.1997856771715278,
                -0.20013391289621135,
                -0.20033959749257807,
                -0.20042398633857525,
                -0.20040574751161644,
                -0.20034414577241405,
                -0.20025288859338328,
                -0.2001441080714588,
                -0.20003640296307307,
                -0.19995722999299903,
                -0.19990753537419478,
                -0.19988530873507943,
                -0.1998889286213763,
                -0.1999074048319535,
                -0.19993271649470815,
                -0.1999645719271424,
                -0.19999118723027182,
                -0.20001208048044197,
                -0.2000251660101789,
                -0.20003043560247108,
                -0.20002940386795437]

        G2_1 = [-0.425,
                -0.4095390866473029,
                -0.3330261455117935,
                -0.3045361465874291,
                -0.28221732217052314,
                -0.2652340446832371,
                -0.2501595186873826,
                -0.23607452401530452,
                -0.2225138088487703,
                -0.20990925129514534,
                -0.19970519353567634,
                -0.1904796214078515,
                -0.18496708357113806,
                -0.18039267065704584,
                -0.1819332181825171,
                -0.18823319476721406,
                -0.19502973847919994,
                -0.19951923246334066,
                -0.20219706511089913,
                -0.20347744635776302,
                -0.20399086903022218,
                -0.20401097008766697,
                -0.203371464109567,
                -0.20224557572542864,
                -0.20125030323360776,
                -0.20064014855972906,
                -0.19993950144620318,
                -0.19933894788386142,
                -0.19890964495280128,
                -0.19881290106390062,
                -0.1990178518631444,
                -0.19936593683758402,
                -0.1996877428851822,
                -0.19994792126027844,
                -0.20013728564641609,
                -0.20024356225509812,
                -0.20027922030413536,
                -0.2002686205152392,
                -0.20021540846484787,
                -0.20014154341367268,
                -0.20006551220087915,
                -0.20000049130415054,
                -0.19995346292199234,
                -0.1999253105355237,
                -0.1999182807741432,
                -0.1999246810523886,
                -0.199943538416856,
                -0.1999648479151982,
                -0.19998476249379993,
                -0.2000021477692965,
                -0.20001411095765922]

        G3_1 = [-0.44166666666666665,
                -0.3635682038429726,
                -0.32681074017864536,
                -0.3073528849348301,
                -0.29041148666354705,
                -0.27272085688313713,
                -0.2578276174888583,
                -0.24324568472071936,
                -0.22935221574285544,
                -0.21628474243957388,
                -0.20467719246908095,
                -0.19472416162725364,
                -0.18722281112041922,
                -0.18188358783594596,
                -0.1806792087036032,
                -0.18499418192345563,
                -0.19188912245108872,
                -0.19745608471860293,
                -0.20095975518130557,
                -0.20288131481341065,
                -0.20375123993116853,
                -0.20403960012188083,
                -0.2038519542729991,
                -0.20280274669983364,
                -0.20171368524339284,
                -0.20097220959433562,
                -0.20028959247814912,
                -0.19961867521172302,
                -0.1990902275489655,
                -0.19882560524809037,
                -0.19886939995385047,
                -0.19919371901595323,
                -0.19952986166513798,
                -0.19982460177167638,
                -0.200052342822042,
                -0.20019718759186875,
                -0.20026664217667609,
                -0.2002749053895271,
                -0.2002493500781802,
                -0.20017988169955314,
                -0.20010239562313728,
                -0.20003089587642928,
                -0.19997451650641843,
                -0.19993720636696186,
                -0.19992039920081678,
                -0.19992087567798306,
                -0.19993239819139774,
                -0.19995390837570612,
                -0.19997524276142248,
                -0.19999418996113488,
                -0.20000874665437668]


        # import pprint
        # pprint.pprint(gauge_values[0])
        # pprint.pprint(gauge_values[1])
        # pprint.pprint(gauge_values[2])
        # pprint.pprint(gauge_values[3])






        assert num.allclose(gauge_values[0], G0) or num.allclose(gauge_values[0], G0_1)
        assert num.allclose(gauge_values[1], G1) or num.allclose(gauge_values[1], G1_1)
        assert num.allclose(gauge_values[2], G2) or num.allclose(gauge_values[2], G2_1)
        assert num.allclose(gauge_values[3], G3) or num.allclose(gauge_values[3], G3_1)

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

        os.remove(os.path.join(domain.get_datadir(), domain.get_name() + '.sww'))

    #####################################################


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
        I = Rate_operator.inflow(domain, rate=2.0, center=(1,1), radius=1)

        A = I.areas.sum()
        assert num.allclose(A, 4) # Two triangles

        stage0 = domain.quantities['stage'].centroid_values.copy()
        domain.timestep = 1.0
        I()
        dstage = domain.quantities['stage'].centroid_values - stage0

        assert num.allclose(dstage[1], 2.0/A)
        assert num.allclose(dstage[0], 2.0/A)
        assert num.allclose(dstage[2:], 0)


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
        I = Rate_operator.inflow(domain, rate=lambda t: 2., center=(1,1), radius=1)

        A = I.areas.sum()
        assert num.allclose(A, 4) # Two triangles

        stage0 = domain.quantities['stage'].centroid_values.copy()
        domain.timestep = 1.0
        I()
        dstage = domain.quantities['stage'].centroid_values - stage0

        assert num.allclose(dstage[1], 2.0/A)
        assert num.allclose(dstage[0], 2.0/A)
        assert num.allclose(dstage[2:], 0)


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
            Rate_operator.inflow(domain, rate=2.0, center=(1,1.1), radius=0.01)
        except Exception:
            pass
        else:
            msg = 'Should have raised exception'
            raise Exception(msg)


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


        os.remove(os.path.join(domain.get_datadir(), domain.get_name() + '.sww'))


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

        os.remove(os.path.join(domain.get_datadir(), domain.get_name() + '.sww'))


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

        os.remove(os.path.join(domain.get_datadir(), domain.get_name() + '.sww'))


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
        fid = NetCDFFile(os.path.join(domain.get_datadir(), domain.get_name() + '.sww'))
        time = fid.variables['time'][:]
        stage = fid.variables['stage'][:,:]
        fid.close()

        os.remove(os.path.join(domain.get_datadir(), domain.get_name() + '.sww'))

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
        fid = NetCDFFile(os.path.join(domain.get_datadir(), domain.get_name() + '.sww'))
        time = fid.variables['time'][:]
        stage = fid.variables['stage'][:,:]
        fid.close()

        os.remove(os.path.join(domain.get_datadir(), domain.get_name() + '.sww'))

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
            #assert allclose (xmom, initial_xmom)

        os.remove(os.path.join(domain.get_datadir(), domain.get_name() + '.sww'))


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
            #assert allclose (xmom, initial_xmom)

        os.remove(os.path.join(domain.get_datadir(), domain.get_name() + '.sww'))


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

        os.remove(os.path.join(domain.get_datadir(), domain.get_name() + '.sww'))


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

        os.remove(os.path.join(domain.get_datadir(), domain.get_name() + '.sww'))


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

            if num.allclose(t, 10):    # Steady state reached
                steady_xmom = domain.quantities['xmomentum'].get_integral()
                steady_ymom = domain.quantities['ymomentum'].get_integral()
                steady_stage = domain.quantities['stage'].get_integral()

            if t > 10:
                msg = 'xmom=%.2f, steady_xmom=%.2f' % (xmom, steady_xmom)
                assert num.allclose(xmom, steady_xmom), msg
                assert num.allclose(ymom, steady_ymom)
                assert num.allclose(stage, steady_stage)

        os.remove(os.path.join(domain.get_datadir(), domain.get_name() + '.sww'))


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
        class Set_IC:
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

        os.remove(os.path.join(domain.get_datadir(), domain.get_name() + '.sww'))


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

        os.remove(os.path.join(domain.get_datadir(), domain.get_name() + '.sww'))


    def test_pmesh2Domain(self):
         import tempfile

         fd, fileName = tempfile.mkstemp(".tsh")
         os.close(fd)
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
        fd, ptsfile = tempfile.mkstemp(".txt")
        os.close(fd)
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
        os.remove(ptsfile)


    def test_fitting_in_hole(self):
        '''
            Make sure we can fit a mesh that has a hole in it.
            This is a regression test for ticket:234

        '''
        verbose = False

        from anuga.shallow_water.shallow_water_domain import Domain
        from anuga.pmesh.mesh_interface import create_pmesh_from_regions
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
        create_pmesh_from_regions(bounding_polygon,
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
        from anuga.pmesh.mesh_interface import create_pmesh_from_regions
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
        create_pmesh_from_regions(bounding_polygon,
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
        from anuga.pmesh.mesh_interface import create_pmesh_from_regions
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

        create_pmesh_from_regions(bounding_polygon,
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
        inflow_rate = 10.00   # m^3/s
        Rate_operator.inflow(domain,
                             center=(10.0, 10.0),
                             radius=5.00,
                             rate=inflow_rate)

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
            ref_volume += ys * inflow_rate


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
        rain_rate = 10.00   # mm/s
        fixed_rain = Rate_operator.rainfall(domain,
                                            center=(10.0, 10.0),
                                            radius=5.00,
                                            rate=rain_rate*3600)   # mm/hr

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


            # Update reference volume (mm/s -> m/s over the rained area)
            delta_V = rain_rate/1000.0*fixed_rain.areas.sum()
            ref_volume += ys * delta_V

        os.remove('Rain_volume_test.sww')


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

            domain.add_quantity('elevation', pole_increment, location='centroids')


        # Check that quantities have been stored correctly
        sww_file = os.path.join(domain.get_datadir(), domain.get_name() + '.sww')
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


    @pytest.mark.slow


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
                inflow_rate = 10.00   # m^3/s

                # Stack this flow
                for i in range(number_of_inflows):
                    Rate_operator.inflow(domain,
                                         center=(10.0, 10.0),
                                         radius=5.00,
                                         rate=inflow_rate)

                ref_flow = inflow_rate*number_of_inflows

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



    # Test smoothing


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
        fd, ptsfile = tempfile.mkstemp(".txt")
        os.close(fd)
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
        try:
            os.remove(ptsfile)
            os.remove(txt_file)
            os.remove(txt_file_prj)
            os.remove(txt_file_dem)
        except OSError:
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


class Test_OMP_Num_Threads(unittest.TestCase):
    """OpenMP thread count is process-wide, so anuga.set_omp_num_threads(n)
    must be reflected by every domain in the session, including ones already
    constructed (the notebook use-case)."""

    def setUp(self):
        from anuga.shallow_water.shallow_water_domain import get_omp_num_threads
        self._saved = get_omp_num_threads()
        self._saved_env = os.environ.get('OMP_NUM_THREADS')

    def tearDown(self):
        from anuga.shallow_water.shallow_water_domain import set_omp_num_threads
        set_omp_num_threads(self._saved, verbose=False)
        if self._saved_env is None:
            os.environ.pop('OMP_NUM_THREADS', None)
        else:
            os.environ['OMP_NUM_THREADS'] = self._saved_env

    def _make_domain(self):
        points, vertices, boundary = rectangular_cross(4, 4)
        return Domain(points, vertices, boundary)

    def test_module_level_propagates_to_existing_and_new_domains(self):
        from anuga.shallow_water.shallow_water_domain import (
            set_omp_num_threads, get_omp_num_threads)

        set_omp_num_threads(1, verbose=False)
        d1 = self._make_domain()
        self.assertEqual(d1.omp_num_threads, 1)

        # Bump the session count AFTER d1 already exists.
        set_omp_num_threads(3, verbose=False)
        self.assertEqual(get_omp_num_threads(), 3)
        self.assertEqual(d1.omp_num_threads, 3)          # existing domain updates
        self.assertEqual(os.environ['OMP_NUM_THREADS'], '3')

        d2 = self._make_domain()
        self.assertEqual(d2.omp_num_threads, 3)          # new domain inherits it

    def test_backward_compatible_setters_are_process_wide(self):
        from anuga.shallow_water.shallow_water_domain import set_omp_num_threads

        set_omp_num_threads(1, verbose=False)
        d1 = self._make_domain()
        d2 = self._make_domain()

        # Legacy attribute assignment still works and is process-wide.
        d1.omp_num_threads = 2
        self.assertEqual(d1.omp_num_threads, 2)
        self.assertEqual(d2.omp_num_threads, 2)

        # Legacy Domain.set_omp_num_threads() method still works and is process-wide.
        d2.set_omp_num_threads(4, verbose=False)
        self.assertEqual(d1.omp_num_threads, 4)
        self.assertEqual(d2.omp_num_threads, 4)


#################################################################################

if __name__ == "__main__":
    suite = unittest.TestSuite([
    Test_Shallow_Water('test_another_runup_example')
])
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)



if __name__ == "__main__":
    unittest.main()
