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


# Variable windfield implemented using functions
def speed(t, x, y):
    """Large speeds halfway between center and edges

    Low speeds at center and edges
    """

    from math import exp, cos, pi

    x = num.array(x)
    y = num.array(y)

    N = len(x)
    s = 0*x  #New array

    for k in range(N):
        r = num.sqrt(x[k]**2 + y[k]**2)
        factor = exp(-(r-0.15)**2)
        s[k] = 4000 * factor * (cos(t*2*pi/150) + 2)

    return s

def scalar_func(t, x, y):
    """Function that returns a scalar.

    Used to test error message when numeric array is expected
    """

    return 17.7

def scalar_func_list(t, x, y):
    """Function that returns a scalar.

    Used to test error message when numeric array is expected
    """

    return [17.7]


def angle(t, x, y):
    """Rotating field
    """
    from math import atan, pi

    x = num.array(x)
    y = num.array(y)

    N = len(x)
    a = 0 * x    # New array

    for k in range(N):
        r = num.sqrt(x[k]**2 + y[k]**2)

        angle = atan(y[k]/x[k])

        if x[k] < 0:
            angle += pi

        # Take normal direction
        angle -= pi/2

        # Ensure positive radians
        if angle < 0:
            angle += 2*pi

        a[k] = angle/pi*180

    return a


###############################################################################

class Test_swb_forcing_terms(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

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
            return slope(x, y) + h

        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage)

        for name in domain.conserved_quantities:
            assert num.allclose(domain.quantities[name].explicit_update, 0)
            assert num.allclose(domain.quantities[name].semi_implicit_update, 0)

        domain.compute_forcing_terms()

        assert num.allclose(domain.quantities['stage'].explicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].explicit_update,
                            -g*h*3)
        assert num.allclose(domain.quantities['ymomentum'].explicit_update, 0)

    def test_manning_friction(self):
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
                            -g*h*3)
        assert num.allclose(domain.quantities['ymomentum'].explicit_update, 0)

        assert num.allclose(domain.quantities['stage'].semi_implicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].semi_implicit_update,
                            0)
        assert num.allclose(domain.quantities['ymomentum'].semi_implicit_update,
                            0)

        #Create some momentum for friction to work with
        domain.set_quantity('xmomentum', 1)
        dz = sqrt(10.0)
        S = -g*eta**2 *dz / h**(7.0/3) 

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

        S = -g*eta**2*5*dz / h**(7.0/3) 

        domain.compute_forcing_terms()

        assert num.allclose(domain.quantities['stage'].semi_implicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].semi_implicit_update,
                            3*S)
        assert num.allclose(domain.quantities['ymomentum'].semi_implicit_update,
                            4*S)



    def test_manning_friction_old(self):
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

        #Turn old mannings function on
        domain.set_new_mannings_function(False)

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
                            -g*h*3)
        assert num.allclose(domain.quantities['ymomentum'].explicit_update, 0)

        assert num.allclose(domain.quantities['stage'].semi_implicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].semi_implicit_update,
                            0)
        assert num.allclose(domain.quantities['ymomentum'].semi_implicit_update,
                            0)

        #Create some momentum for friction to work with
        domain.set_quantity('xmomentum', 1)
        S = -g*eta**2 / h**(7.0/3)

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

        S = -g*eta**2*5 / h**(7.0/3)

        domain.compute_forcing_terms()

        assert num.allclose(domain.quantities['stage'].semi_implicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].semi_implicit_update,
                            3*S)
        assert num.allclose(domain.quantities['ymomentum'].semi_implicit_update,
                            4*S)


    def test_manning_friction_new(self):
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

        # Use the new function which takes into account the extra
        # wetted area due to slope of bed
        domain.set_new_mannings_function(True)
        
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
                            -g*h*3)
        assert num.allclose(domain.quantities['ymomentum'].explicit_update, 0)

        assert num.allclose(domain.quantities['stage'].semi_implicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].semi_implicit_update,
                            0)
        assert num.allclose(domain.quantities['ymomentum'].semi_implicit_update,
                            0)

        #Create some momentum for friction to work with
        domain.set_quantity('xmomentum', 1)
        S = -g*eta**2 / h**(7.0/3) * sqrt(10)

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

        S = -g*eta**2*5 / h**(7.0/3) * sqrt(10.0)

        domain.compute_forcing_terms()

        assert num.allclose(domain.quantities['stage'].semi_implicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].semi_implicit_update,
                            3*S)
        assert num.allclose(domain.quantities['ymomentum'].semi_implicit_update,
                            4*S)
        
    def test_constant_wind_stress(self):
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

        domain = anuga.Domain(points, vertices)

        #Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        #Setup only one forcing term, constant wind stress
        s = 100
        phi = 135
        domain.forcing_terms = []
        domain.forcing_terms.append(anuga.Wind_stress(s, phi))

        domain.compute_forcing_terms()

        const = eta_w*rho_a / rho_w

        #Convert to radians
        phi = phi*pi / 180

        #Compute velocity vector (u, v)
        u = s*cos(phi)
        v = s*sin(phi)

        #Compute wind stress
        S = const * num.sqrt(u**2 + v**2)

        assert num.allclose(domain.quantities['stage'].explicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].explicit_update, S*u)
        assert num.allclose(domain.quantities['ymomentum'].explicit_update, S*v)

    def test_variable_wind_stress(self):
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

        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        domain.time = 5.54    # Take a random time (not zero)

        #Setup only one forcing term, constant wind stress
        s = 100
        phi = 135
        domain.forcing_terms = []
        domain.forcing_terms.append(anuga.Wind_stress(s=speed, phi=angle))

        domain.compute_forcing_terms()

        #Compute reference solution
        const = eta_w*rho_a / rho_w

        N = len(domain)    # number_of_triangles

        xc = domain.get_centroid_coordinates()
        t = domain.time

        x = xc[:,0]
        y = xc[:,1]
        s_vec = speed(t,x,y)
        phi_vec = angle(t,x,y)

        for k in range(N):
            # Convert to radians
            phi = phi_vec[k]*pi / 180
            s = s_vec[k]

            # Compute velocity vector (u, v)
            u = s*cos(phi)
            v = s*sin(phi)

            # Compute wind stress
            S = const * num.sqrt(u**2 + v**2)

            assert num.allclose(domain.quantities['stage'].explicit_update[k],
                                0)
            assert num.allclose(domain.quantities['xmomentum'].\
                                     explicit_update[k],
                                S*u)
            assert num.allclose(domain.quantities['ymomentum'].\
                                     explicit_update[k],
                                S*v)

    def test_windfield_from_file(self):
        import time
        from anuga.config import rho_a, rho_w, eta_w
        from math import pi, cos, sin
        from anuga.config import time_format
        from anuga.abstract_2d_finite_volumes.util import file_function

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

        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        domain.time = 7    # Take a time that is represented in file (not zero)

        # Write wind stress file (ensure that domain.time is covered)
        # Take x=1 and y=0
        filename = 'test_windstress_from_file'
        start = time.mktime(time.strptime('2000', '%Y'))
        fid = open(filename + '.txt', 'w')
        dt = 1    # One second interval
        t = 0.0
        while t <= 10.0:
            t_string = time.strftime(time_format, time.gmtime(t+start))

            fid.write('%s, %f %f\n' %
                      (t_string, speed(t,[1],[0])[0], angle(t,[1],[0])[0]))
            t += dt

        fid.close()

        anuga.timefile2netcdf(filename+'.txt', filename+'.tms')
        os.remove(filename + '.txt')

        # Setup wind stress
        F = file_function(filename + '.tms',
                          quantities=['Attribute0', 'Attribute1'])
        os.remove(filename + '.tms')

        W = anuga.Wind_stress(F)

        domain.forcing_terms = []
        domain.forcing_terms.append(W)

        domain.compute_forcing_terms()

        # Compute reference solution
        const = eta_w*rho_a / rho_w

        N = len(domain)    # number_of_triangles

        t = domain.time

        s = speed(t, [1], [0])[0]
        phi = angle(t, [1], [0])[0]

        # Convert to radians
        phi = phi*pi / 180

        # Compute velocity vector (u, v)
        u = s*cos(phi)
        v = s*sin(phi)

        # Compute wind stress
        S = const * num.sqrt(u**2 + v**2)

        for k in range(N):
            assert num.allclose(domain.quantities['stage'].explicit_update[k],
                                0)
            assert num.allclose(domain.quantities['xmomentum'].\
                                    explicit_update[k],
                                S*u)
            assert num.allclose(domain.quantities['ymomentum'].\
                                    explicit_update[k],
                                S*v)

    def test_windfield_from_file_seconds(self):
        import time
        from anuga.config import rho_a, rho_w, eta_w
        from math import pi, cos, sin
        from anuga.config import time_format
        from anuga.abstract_2d_finite_volumes.util import file_function

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

        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        domain.time = 7    # Take a time that is represented in file (not zero)

        # Write wind stress file (ensure that domain.time is covered)
        # Take x=1 and y=0
        filename = 'test_windstress_from_file.txt'
        file_out = 'test_windstress_from_file.tms'
        start = time.mktime(time.strptime('2000', '%Y'))
        fid = open(filename, 'w')
        dt = 0.5    # Half second interval
        t = 0.0
        while t <= 10.0:
            fid.write('%s, %f %f\n'
                      % (str(t), speed(t, [1], [0])[0], angle(t, [1], [0])[0]))
            t += dt

        fid.close()

        anuga.timefile2netcdf(filename, file_out, time_as_seconds=True)
        os.remove(filename)

        # Setup wind stress
        F = file_function(file_out,
                          quantities=['Attribute0', 'Attribute1'])
        os.remove(file_out)

        W = anuga.Wind_stress(F)

        domain.forcing_terms = []
        domain.forcing_terms.append(W)

        domain.compute_forcing_terms()

        # Compute reference solution
        const = eta_w*rho_a / rho_w

        N = len(domain)    # number_of_triangles

        t = domain.time

        s = speed(t, [1], [0])[0]
        phi = angle(t, [1], [0])[0]

        # Convert to radians
        phi = phi*pi / 180

        # Compute velocity vector (u, v)
        u = s*cos(phi)
        v = s*sin(phi)

        # Compute wind stress
        S = const * num.sqrt(u**2 + v**2)

        for k in range(N):
            assert num.allclose(domain.quantities['stage'].explicit_update[k],
                                0)
            assert num.allclose(domain.quantities['xmomentum'].\
                                    explicit_update[k],
                                S*u)
            assert num.allclose(domain.quantities['ymomentum'].\
                                    explicit_update[k],
                                S*v)

    def test_wind_stress_error_condition(self):
        """Test that windstress reacts properly when forcing functions
        are wrong - e.g. returns a scalar
        """

        from math import pi, cos, sin
        from anuga.config import rho_a, rho_w, eta_w

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

        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        domain.time = 5.54    # Take a random time (not zero)

        # Setup only one forcing term, bad func
        domain.forcing_terms = []

        try:
            domain.forcing_terms.append(anuga.Wind_stress(s=scalar_func_list,
                                                    phi=angle))
        except AssertionError:
            pass
        else:
            msg = 'Should have raised exception'
            raise Exception, msg

        try:
            domain.forcing_terms.append(Wind_stress(s=speed, phi=scalar_func))
        except Exception:
            pass
        else:
            msg = 'Should have raised exception'
            raise Exception, msg

        try:
            domain.forcing_terms.append(Wind_stress(s=speed, phi='xx'))
        except:
            pass
        else:
            msg = 'Should have raised exception'
            raise Exception, msg

    def test_rainfall(self):
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

        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        # Setup only one forcing term, constant rainfall
        domain.forcing_terms = []
        domain.forcing_terms.append(anuga.Rainfall(domain, rate=2.0))

        domain.compute_forcing_terms()
        assert num.allclose(domain.quantities['stage'].explicit_update,
                            2.0/1000)

    def test_rainfall_restricted_by_polygon(self):
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

        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        # Setup only one forcing term, constant rainfall
        # restricted to a polygon enclosing triangle #1 (bce)
        domain.forcing_terms = []
        R = anuga.Rainfall(domain, rate=2.0, polygon=[[1,1], [2,1], [2,2], [1,2]])

        assert num.allclose(R.exchange_area, 2)

        domain.forcing_terms.append(R)

        domain.compute_forcing_terms()

        assert num.allclose(domain.quantities['stage'].explicit_update[1],
                            2.0/1000)
        assert num.allclose(domain.quantities['stage'].explicit_update[0], 0)
        assert num.allclose(domain.quantities['stage'].explicit_update[2:], 0)

    def test_time_dependent_rainfall_restricted_by_polygon(self):
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

        # Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        # Setup only one forcing term, time dependent rainfall
        # restricted to a polygon enclosing triangle #1 (bce)
        domain.forcing_terms = []
        R = anuga.Rainfall(domain,
                     rate=lambda t: 3*t + 7,
                     polygon = [[1,1], [2,1], [2,2], [1,2]])

        assert num.allclose(R.exchange_area, 2)
        
        domain.forcing_terms.append(R)

        domain.time = 10.

        domain.compute_forcing_terms()

        assert num.allclose(domain.quantities['stage'].explicit_update[1],
                            (3*domain.time + 7)/1000)
        assert num.allclose(domain.quantities['stage'].explicit_update[0], 0)
        assert num.allclose(domain.quantities['stage'].explicit_update[2:], 0)

    def test_time_dependent_rainfall_using_starttime(self):
        rainfall_poly = ensure_numeric([[1,1], [2,1], [2,2], [1,2]], num.float)

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

        # Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        # Setup only one forcing term, time dependent rainfall
        # restricted to a polygon enclosing triangle #1 (bce)
        domain.forcing_terms = []
        R = anuga.Rainfall(domain,
                     rate=lambda t: 3*t + 7,
                     polygon=rainfall_poly)                     

        assert num.allclose(R.exchange_area, 2)
        
        domain.forcing_terms.append(R)

        # This will test that time used in the forcing function takes
        # startime into account.
        domain.starttime = 5.0

        domain.time = 7.

        domain.compute_forcing_terms()

        assert num.allclose(domain.quantities['stage'].explicit_update[1],
                            (3*domain.get_time() + 7)/1000)
        assert num.allclose(domain.quantities['stage'].explicit_update[1],
                            (3*(domain.time + domain.starttime) + 7)/1000)

        # Using internal time her should fail
        assert not num.allclose(domain.quantities['stage'].explicit_update[1],
                                (3*domain.time + 7)/1000)

        assert num.allclose(domain.quantities['stage'].explicit_update[0], 0)
        assert num.allclose(domain.quantities['stage'].explicit_update[2:], 0)

    def test_time_dependent_rainfall_using_georef(self):
        """test_time_dependent_rainfall_using_georef

        This will also test the General forcing term using georef
        """

        # Mesh in zone 56 (absolute coords)
        x0 = 314036.58727982
        y0 = 6224951.2960092

        rainfall_poly = ensure_numeric([[1,1], [2,1], [2,2], [1,2]], num.float)
        rainfall_poly += [x0, y0]

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #             bac,     bce,     ecf,     dbe
        vertices = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices,
                        geo_reference=Geo_reference(56, x0, y0))

        # Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        # Setup only one forcing term, time dependent rainfall
        # restricted to a polygon enclosing triangle #1 (bce)
        domain.forcing_terms = []
        R = anuga.Rainfall(domain,
                     rate=lambda t: 3*t + 7,
                     polygon=rainfall_poly)

        assert num.allclose(R.exchange_area, 2)
        
        domain.forcing_terms.append(R)

        # This will test that time used in the forcing function takes
        # startime into account.
        domain.starttime = 5.0

        domain.time = 7.

        domain.compute_forcing_terms()

        assert num.allclose(domain.quantities['stage'].explicit_update[1],
                            (3*domain.get_time() + 7)/1000)
        assert num.allclose(domain.quantities['stage'].explicit_update[1],
                            (3*(domain.time + domain.starttime) + 7)/1000)

        # Using internal time her should fail
        assert not num.allclose(domain.quantities['stage'].explicit_update[1],
                                (3*domain.time + 7)/1000)

        assert num.allclose(domain.quantities['stage'].explicit_update[0], 0)
        assert num.allclose(domain.quantities['stage'].explicit_update[2:], 0)

    def test_time_dependent_rainfall_restricted_by_polygon_with_default(self):
        """
        Test that default rainfall can be used when given rate runs out of data.
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

        # Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        # Setup only one forcing term, time dependent rainfall
        # that expires at t==20
        from anuga.fit_interpolate.interpolate import Modeltime_too_late

        def main_rate(t):
            if t > 20:
                msg = 'Model time exceeded.'
                raise Modeltime_too_late, msg
            else:
                return 3*t + 7

        domain.forcing_terms = []
        R = anuga.Rainfall(domain,
                     rate=main_rate,
                     polygon = [[1,1], [2,1], [2,2], [1,2]],
                     default_rate=5.0)

        assert num.allclose(R.exchange_area, 2)
        
        domain.forcing_terms.append(R)

        domain.time = 10.

        domain.compute_forcing_terms()

        assert num.allclose(domain.quantities['stage'].explicit_update[1],
                            (3*domain.time+7)/1000)
        assert num.allclose(domain.quantities['stage'].explicit_update[0], 0)
        assert num.allclose(domain.quantities['stage'].explicit_update[2:], 0)

        domain.time = 100.
        domain.quantities['stage'].explicit_update[:] = 0.0     # Reset
        domain.compute_forcing_terms()

        assert num.allclose(domain.quantities['stage'].explicit_update[1],
                            5.0/1000) # Default value
        assert num.allclose(domain.quantities['stage'].explicit_update[0], 0)
        assert num.allclose(domain.quantities['stage'].explicit_update[2:], 0)

    def test_rainfall_forcing_with_evolve(self):
        """test_rainfall_forcing_with_evolve

        Test how forcing terms are called within evolve
        """

        # FIXME(Ole): This test is just to experiment

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

        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        # Setup only one forcing term, time dependent rainfall
        # that expires at t==20
        from anuga.fit_interpolate.interpolate import Modeltime_too_late

        def main_rate(t):
            if t > 20:
                msg = 'Model time exceeded.'
                raise Modeltime_too_late, msg
            else:
                return 3*t + 7

        domain.forcing_terms = []
        R = anuga.Rainfall(domain,
                     rate=main_rate,
                     polygon=[[1,1], [2,1], [2,2], [1,2]],
                     default_rate=5.0)

        assert num.allclose(R.exchange_area, 2)
        
        domain.forcing_terms.append(R)

        for t in domain.evolve(yieldstep=1, finaltime=25):
            pass
            #FIXME(Ole):  A test here is hard because explicit_update also
            # receives updates from the flux calculation.



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

        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        # Setup only one forcing term, constant inflow of 2 m^3/s
        # on a circle affecting triangles #0 and #1 (bac and bce)
        domain.forcing_terms = []
        
        I = anuga.Inflow(domain, rate=2.0, center=(1,1), radius=1)
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

        domain = anuga.Domain(points, vertices)

        # Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        # Setup only one forcing term, time dependent inflow of 2 m^3/s
        # on a circle affecting triangles #0 and #1 (bac and bce)
        domain.forcing_terms = []
        I = anuga.Inflow(domain, rate=lambda t: 2., center=(1,1), radius=1)
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

        Br = anuga.Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        # Setup only one forcing term, constant inflow of 2 m^3/s
        # on a circle affecting triangles #0 and #1 (bac and bce)
        try:
            Inflow(domain, rate=2.0, center=(1,1.1), radius=0.01)
        except:
            pass
        else:
            msg = 'Should have raised exception'
            raise Exception, msg

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
        points, vertices, boundary = rectangular_cross(int(length/dx),
                                                       int(width/dy),
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
            domain.time = 0.0
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
                predicted_volume = predicted_volume + 2.0/pi/100/dt # Why 100?

        # Apply equivalent outflow only and check volumes
        # for a range of stage values
        for stage in [2.0, 1.0, 0.5, 0.25, 0.1, 0.0]:
            domain.time = 0.0
            domain.set_quantity('stage', stage)
            domain.forcing_terms = []
            domain.forcing_terms.append(Inflow(domain, rate=-2.0,
                                               center=(15,5), radius=1))
            initial_volume = domain.quantities['stage'].get_integral()
            predicted_volume = initial_volume
            dt = 0.05
            for t in domain.evolve(yieldstep=dt, finaltime=5.0):
                volume = domain.quantities['stage'].get_integral()
                assert num.allclose (volume, predicted_volume)
                predicted_volume = predicted_volume - 2.0/pi/100/dt # Why 100?

        # Apply both inflow and outflow and check volumes being constant for a
        # range of stage values
        for stage in [2.0, 1.0, 0.5, 0.25, 0.1, 0.0]:
            domain.time = 0.0
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
                assert num.allclose(volume, initial_volume)

    #####################################################

        
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
        finaltime = 500 #700.0 # If this is too short, steady state will not be achieved

        length = 250.
        width  = 20.
        dx = dy = 5                 # Resolution: of grid on both axes

        points, vertices, boundary = rectangular_cross(int(length/dx),
                                                       int(width/dy),
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
                fixed_inflow = anuga.Inflow(domain,
                                      center=(10.0, 10.0),
                                      radius=5.00,
                                      rate=10.00)   

                # Stack this flow
                for i in range(number_of_inflows):
                    domain.forcing_terms.append(fixed_inflow)
                
                ref_flow = fixed_inflow.rate*number_of_inflows

                # Compute normal depth on plane using Mannings equation
                # v=1/n*(r^2/3)*(s^0.5) or r=(Q*n/(s^0.5*W))^0.6
                normal_depth=(ref_flow*mannings_n/(slope**0.5*width))**0.6
                if verbose:
                    print
                    print 'Slope:', slope, 'Mannings n:', mannings_n
                    
                    
                #--------------------------------------------------------------
                # Setup boundary conditions
                #--------------------------------------------------------------

                Br = anuga.Reflective_boundary(domain)
                
                # Define downstream boundary based on predicted depth
                def normal_depth_stage_downstream(t):
                    return (-slope*length) + normal_depth
                
                Bt = anuga.Transmissive_momentum_set_stage_boundary(
                        domain=domain, function=normal_depth_stage_downstream)

                domain.set_boundary({'left': Br,
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
                q = domain.get_flow_through_cross_section([[200.0,  0.0],
                                                           [200.0, 20.0]])
                if verbose:
                    print ('90 degree flowline: ANUGA = %f, Ref = %f'
                           % (q, ref_flow))

                msg = ('Predicted flow was %f, should have been %f'
                       % (q, ref_flow))
                assert num.allclose(q, ref_flow, rtol=1.0e-2), msg         

                           
                # 45 degree flowline at 200m
                q = domain.get_flow_through_cross_section([[200.0,  0.0],
                                                           [220.0, 20.0]])
                if verbose:
                    print ('45 degree flowline: ANUGA = %f, Ref = %f'
                           % (q, ref_flow))
                    
                msg = ('Predicted flow was %f, should have been %f'
                       % (q, ref_flow))
                assert num.allclose(q, ref_flow, rtol=1.0e-2), msg         

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
        # Import necessary modules
        #----------------------------------------------------------------------
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
        from anuga.shallow_water import Domain
        from anuga.shallow_water.shallow_water_domain import Reflective_boundary
        from anuga.shallow_water.shallow_water_domain import Dirichlet_boundary
        from anuga.shallow_water.forcing import Inflow_boundary
        from anuga.shallow_water.data_manager import get_flow_through_cross_section
        from anuga.abstract_2d_finite_volumes.util import sww2csv_gauges, csv2timeseries_graphs


        #----------------------------------------------------------------------
        # Setup computational domain
        #----------------------------------------------------------------------

        finaltime = 500 #700.0 # If this is too short, steady state will not be achieved

        length = 250.
        width  = 20.
        dx = dy = 5          # Resolution: of grid on both axes

        points, vertices, boundary = rectangular_cross(int(length/dx), int(width/dy),
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
                normal_depth=(ref_flow*mannings_n/(slope**0.5*width))**0.6
                if verbose:
                    print
                    print 'Slope:', slope, 'Mannings n:', mannings_n
                    
                
                
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
                    print '90 degree flowline: ANUGA = %f, Ref = %f' % (q, ref_flow)
                assert num.allclose(q, ref_flow, rtol=1.0e-2), msg         

                           
                # 45 degree flowline at 200m
                q=domain.get_flow_through_cross_section([[200.0,0.0],[220.0,20.0]])
                msg = 'Predicted flow was %f, should have been %f' % (q, ref_flow)
                if verbose:
                    print '45 degree flowline: ANUGA = %f, Ref = %f' % (q, ref_flow)
                    
                assert num.allclose(q, ref_flow, rtol=1.0e-2), msg         

        
 ################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_swb_forcing_terms, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
