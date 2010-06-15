"""  Test environmental forcing - rain, wind, etc.
"""

import unittest, os
from anuga.shallow_water.shallow_water_domain import Domain
from boundaries import Reflective_boundary
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.file_conversion.file_conversion import timefile2netcdf
from forcing import *

import numpy as num
import warnings


def scalar_func_list(t, x, y):
    """Function that returns a scalar.

    Used to test error message when numeric array is expected
    """

    return [17.7]


def speed(t, x, y):
    """
    Variable windfield implemented using functions    
    Large speeds halfway between center and edges

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

    

class Test_Forcing(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass
        
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

        domain = Domain(points, vertices)

        #Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        #Setup only one forcing term, constant wind stress
        s = 100
        phi = 135
        domain.forcing_terms = []
        domain.forcing_terms.append(Wind_stress(s, phi))

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

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        domain.time = 5.54    # Take a random time (not zero)

        #Setup only one forcing term, constant wind stress
        s = 100
        phi = 135
        domain.forcing_terms = []
        domain.forcing_terms.append(Wind_stress(s=speed, phi=angle))

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

        Br = Reflective_boundary(domain)
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

        timefile2netcdf(filename + '.txt')
        os.remove(filename + '.txt')

        # Setup wind stress
        F = file_function(filename + '.tms',
                          quantities=['Attribute0', 'Attribute1'])
        os.remove(filename + '.tms')

        W = Wind_stress(F)

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

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        domain.time = 7    # Take a time that is represented in file (not zero)

        # Write wind stress file (ensure that domain.time is covered)
        # Take x=1 and y=0
        filename = 'test_windstress_from_file'
        start = time.mktime(time.strptime('2000', '%Y'))
        fid = open(filename + '.txt', 'w')
        dt = 0.5    # Half second interval
        t = 0.0
        while t <= 10.0:
            fid.write('%s, %f %f\n'
                      % (str(t), speed(t, [1], [0])[0], angle(t, [1], [0])[0]))
            t += dt

        fid.close()

        timefile2netcdf(filename + '.txt', time_as_seconds=True)
        os.remove(filename + '.txt')

        # Setup wind stress
        F = file_function(filename + '.tms',
                          quantities=['Attribute0', 'Attribute1'])
        os.remove(filename + '.tms')

        W = Wind_stress(F)

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

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        domain.time = 5.54    # Take a random time (not zero)

        # Setup only one forcing term, bad func
        domain.forcing_terms = []

        try:
            domain.forcing_terms.append(Wind_stress(s=scalar_func_list,
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

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        # Setup only one forcing term, constant rainfall
        domain.forcing_terms = []
        domain.forcing_terms.append(Rainfall(domain, rate=2.0))

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

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        # Setup only one forcing term, constant rainfall
        # restricted to a polygon enclosing triangle #1 (bce)
        domain.forcing_terms = []
        R = Rainfall(domain, rate=2.0, polygon=[[1,1], [2,1], [2,2], [1,2]])

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

        domain = Domain(points, vertices)

        # Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        # Setup only one forcing term, time dependent rainfall
        # restricted to a polygon enclosing triangle #1 (bce)
        domain.forcing_terms = []
        R = Rainfall(domain,
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

        domain = Domain(points, vertices)

        # Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        # Setup only one forcing term, time dependent rainfall
        # restricted to a polygon enclosing triangle #1 (bce)
        domain.forcing_terms = []
        R = Rainfall(domain,
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

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        # Setup only one forcing term, time dependent rainfall
        # restricted to a polygon enclosing triangle #1 (bce)
        domain.forcing_terms = []
        R = Rainfall(domain,
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

        import warnings
        warnings.simplefilter('ignore', UserWarning)


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
        R = Rainfall(domain,
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
        import warnings
        warnings.simplefilter('ignore', UserWarning)


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
        R = Rainfall(domain,
                     rate=main_rate,
                     polygon=[[1,1], [2,1], [2,2], [1,2]],
                     default_rate=5.0)

        assert num.allclose(R.exchange_area, 2)
        
        domain.forcing_terms.append(R)

        for t in domain.evolve(yieldstep=1, finaltime=25):
            pass
            #FIXME(Ole):  A test here is hard because explicit_update also
            # receives updates from the flux calculation.


    def test_rainfall_forcing_with_evolve_1(self):
        """test_rainfall_forcing_with_evolve

        Test how forcing terms are called within evolve.
        This test checks that proper exception is thrown when no default_rate is set
        """

        import warnings
        warnings.simplefilter('ignore', UserWarning)


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
        R = Rainfall(domain,
                     rate=main_rate,
                     polygon=[[1,1], [2,1], [2,2], [1,2]])


        assert num.allclose(R.exchange_area, 2)
        
        domain.forcing_terms.append(R)
        #for t in domain.evolve(yieldstep=1, finaltime=25):
        #    pass
                
        try:
            for t in domain.evolve(yieldstep=1, finaltime=25):
                pass
        except Modeltime_too_late, e:
            # Test that error message is as expected
            assert 'can specify keyword argument default_rate in the forcing function' in str(e)
        else:
            raise Exception, 'Should have raised exception'

            
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Forcing, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
