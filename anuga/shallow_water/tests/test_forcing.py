"""  Test environmental forcing - rain, wind, etc.
"""

import unittest, os
import anuga
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.shallow_water.boundaries import Reflective_boundary
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.file_conversion.file_conversion import timefile2netcdf
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular
from anuga.abstract_2d_finite_volumes.util import file_function
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a

from anuga.shallow_water.forcing import *

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

def time_varying_speed(t, x, y):
    """
    Variable speed windfield
    """

    from math import exp, cos, pi

    x = num.array(x,num.float)
    y = num.array(y,num.float)

    N = len(x)
    s = 0*x  #New array

    #dx=x[-1]-x[0]; dy = y[-1]-y[0]
    S=100.
    for k in range(N):
        s[k]=S*(1.+t/100.)
    return s


def time_varying_angle(t, x, y):
    """Rotating field
    """
    from math import atan, pi

    x = num.array(x,num.float)
    y = num.array(y,num.float)

    N = len(x)
    a = 0 * x    # New array

    phi=135.
    for k in range(N):
        a[k]=phi*(1.+t/100.)

    return a


def time_varying_pressure(t, x, y):
    """Rotating field
    """
    from math import atan, pi

    x = num.array(x,num.float)
    y = num.array(y,num.float)

    N = len(x)
    p = 0 * x    # New array

    p0=1000.
    for k in range(N):
        p[k]=p0*(1.-t/100.)

    return p

def spatial_linear_varying_speed(t, x, y):
    """
    Variable speed windfield
    """

    from math import exp, cos, pi

    x = num.array(x)
    y = num.array(y)

    N = len(x)
    s = 0*x  #New array

    #dx=x[-1]-x[0]; dy = y[-1]-y[0]
    s0=250.
    ymin=num.min(y)
    xmin=num.min(x)
    a=0.000025; b=0.0000125;
    for k in range(N):
        s[k]=s0*(1+t/100.)+a*x[k]+b*y[k]
    return s


def spatial_linear_varying_angle(t, x, y):
    """Rotating field
    """
    from math import atan, pi

    x = num.array(x)
    y = num.array(y)

    N = len(x)
    a = 0 * x    # New array

    phi=135.
    b1=0.000025; b2=0.00001125;
    for k in range(N):
        a[k]=phi*(1+t/100.)+b1*x[k]+b2*y[k]
    return a

def spatial_linear_varying_pressure(t, x, y):
    p0=1000;
    a=0.000025; b=0.0000125;

    x = num.array(x)
    y = num.array(y)

    N = len(x)
    p = 0 * x    # New array

    for k in range(N):
        p[k]=p0*(1.-t/100.)+a*x[k]+b*y[k]
    return p


def grid_1d(x0,dx,nx):
    x = num.empty(nx,dtype=num.float)
    for i in range(nx):
        x[i]=x0+float(i)*dx
    return x
    

def ndgrid(x,y):
    nx = len(x)
    ny = len(y)
    X = num.empty(nx*ny,dtype=num.float)
    Y = num.empty(nx*ny,dtype=num.float)
    k=0
    for i in range(nx):
        for j in range(ny):
            X[k]=x[i]
            Y[k]=y[j]
            k+=1
    return X,Y

class Test_Forcing(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        for file in ['domain.sww']:
            try:
                os.remove(file)
            except:
                pass
        
    def write_wind_pressure_field_sts(self,
                                      field_sts_filename,
                                      nrows=10,
                                      ncols=10,
                                      cellsize=25,
                                      origin=(0.0,0.0),
                                      refzone=50,
                                      timestep=1,
                                      number_of_timesteps=10,
                                      angle=135.0,
                                      speed=100.0,
                                      pressure=1000.0):

        xllcorner=origin[0]
        yllcorner=origin[1]
        starttime = 0; endtime = number_of_timesteps*timestep;
        no_data = -9999

        time = num.arange(starttime, endtime, timestep, dtype='i')

        x = grid_1d(xllcorner,cellsize,ncols)
        y = grid_1d(yllcorner,cellsize,nrows)
        [X,Y] = ndgrid(x,y)
        number_of_points = nrows*ncols

        wind_speed = num.empty((number_of_timesteps,nrows*ncols),dtype=num.float)
        wind_angle = num.empty((number_of_timesteps,nrows*ncols),dtype=num.float)
        barometric_pressure = num.empty((number_of_timesteps,nrows*ncols),
                                        dtype=num.float)

        if ( callable(speed) and callable(angle) and callable(pressure) ):
            x = num.ones(3, num.float)
            y = num.ones(3, num.float)
            try:
                s = speed(1.0, x=x, y=y)
                a = angle(1.0, x=x, y=y)
                p = pressure(1.0, x=x, y=y)
                use_function=True
            except Exception, e:
                msg = 'Function could not be executed.\n'
                raise Exception, msg
        else:
            try :
                speed=float(speed)
                angle=float(angle)
                pressure=float(pressure)
                use_function=False
            except:
                msg = ('Force fields must be a scalar value coercible to float.')
                raise Exception, msg

        for i,t in enumerate(time):
            if ( use_function ):
                wind_speed[i,:] = speed(t,X,Y)
                wind_angle[i,:] = angle(t,X,Y)
                barometric_pressure[i,:] = pressure(t,X,Y)
            else:
                wind_speed[i,:] = speed
                wind_angle[i,:] = angle
                barometric_pressure[i,:] = pressure

        # "Creating the field STS NetCDF file"

        fid = NetCDFFile(field_sts_filename+'.sts', 'w')
        fid.institution = 'Geoscience Australia'
        fid.description = "description"
        fid.starttime = 0.0
        fid.ncols = ncols
        fid.nrows = nrows
        fid.cellsize = cellsize
        fid.no_data = no_data
        fid.createDimension('number_of_points', number_of_points)
        fid.createDimension('number_of_timesteps', number_of_timesteps)
        fid.createDimension('numbers_in_range', 2)

        fid.createVariable('x', 'd', ('number_of_points',))
        fid.createVariable('y', 'd', ('number_of_points',))
        fid.createVariable('time', 'i', ('number_of_timesteps',))
        fid.createVariable('wind_speed', 'd', ('number_of_timesteps', 
                                               'number_of_points'))
        fid.createVariable('wind_speed_range', 'd', ('numbers_in_range', ))
        fid.createVariable('wind_angle', 'd', ('number_of_timesteps', 
                                               'number_of_points'))
        fid.createVariable('wind_angle_range', 'd', ('numbers_in_range',))
        fid.createVariable('barometric_pressure', 'd', ('number_of_timesteps', 
                                             'number_of_points'))
        fid.createVariable('barometric_pressure_range', 'd', ('numbers_in_range',))


        fid.variables['wind_speed_range'][:] = num.array([1e+036, -1e+036])
        fid.variables['wind_angle_range'][:] = num.array([1e+036, -1e+036])
        fid.variables['barometric_pressure_range'][:] = num.array([1e+036, -1e+036])
        fid.variables['time'][:] = time

        ws = fid.variables['wind_speed']
        wa = fid.variables['wind_angle']
        pr = fid.variables['barometric_pressure']

        for i in xrange(number_of_timesteps):
            ws[i] = wind_speed[i,:]
            wa[i] = wind_angle[i,:]
            pr[i] = barometric_pressure[i,:]

        origin = anuga.coordinate_transforms.geo_reference.Geo_reference(refzone,
                                                                         xllcorner,
                                                                         yllcorner)
        geo_ref = anuga.coordinate_transforms.geo_reference.write_NetCDF_georeference(origin, fid)

        fid.variables['x'][:]=X-geo_ref.get_xllcorner()
        fid.variables['y'][:]=Y-geo_ref.get_yllcorner()


        fid.close()

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

        # This will test that time is set to starttime in set_starttime
        domain.set_starttime(5.0)

        domain.compute_forcing_terms()

        assert num.allclose(domain.quantities['stage'].explicit_update[1],
                            (3*domain.get_time() + 7)/1000)
        assert num.allclose(domain.quantities['stage'].explicit_update[1],
                            (3*domain.get_starttime() + 7)/1000)

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

        # This will test that time is set to starttime in set_starttime
        domain.set_starttime(5.0)

        domain.compute_forcing_terms()

        assert num.allclose(domain.quantities['stage'].explicit_update[1],
                            (3*domain.get_time() + 7)/1000)
        assert num.allclose(domain.quantities['stage'].explicit_update[1],
                            (3*domain.get_starttime() + 7)/1000)

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

    def test_constant_wind_stress_from_file(self):
        from anuga.config import rho_a, rho_w, eta_w
        from math import pi, cos, sin
        from anuga.file_conversion.sts2sww_mesh import sts2sww_mesh

        cellsize = 25
        nrows=5; ncols = 6;
        refzone=50
        xllcorner=366000;yllcorner=6369500;
        number_of_timesteps = 6
        timestep=12*60
        eps=2e-16

        points, vertices, boundary =rectangular(nrows-2,ncols-2,
                                                len1=cellsize*(ncols-1),
                                                len2=cellsize*(nrows-1),
                                                origin=(xllcorner,yllcorner))

        domain = Domain(points, vertices, boundary)
        midpoints = domain.get_centroid_coordinates()

        # Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'top': Br, 'bottom' :Br, 'left': Br, 'right': Br})

        # Setup only one forcing term, constant wind stress
        s = 100
        phi = 135
        pressure=1000
        domain.forcing_terms = []
        field_sts_filename = 'wind_field'
        self.write_wind_pressure_field_sts(field_sts_filename,
                                      nrows=nrows,
                                      ncols=ncols,
                                      cellsize=cellsize,
                                      origin=(xllcorner,yllcorner),
                                      refzone=50,
                                      timestep=timestep,
                                      number_of_timesteps=10,
                                      speed=s,
                                      angle=phi,
                                      pressure=pressure)

        sts2sww_mesh(field_sts_filename,spatial_thinning=1,
                     verbose=False)

        # Setup wind stress
        F = file_function(field_sts_filename+'.sww', domain,
                          quantities=['wind_speed', 'wind_angle'],
                          interpolation_points = midpoints)

        W = Wind_stress(F,use_coordinates=False)
        domain.forcing_terms.append(W)
        domain.compute_forcing_terms()

        const = eta_w*rho_a / rho_w

        # Convert to radians
        phi = phi*pi / 180

        # Compute velocity vector (u, v)
        u = s*cos(phi)
        v = s*sin(phi)

        # Compute wind stress
        S = const * num.sqrt(u**2 + v**2)

        assert num.allclose(domain.quantities['stage'].explicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].explicit_update, S*u)
        assert num.allclose(domain.quantities['ymomentum'].explicit_update, S*v)

    def test_variable_windfield_from_file(self):
        from anuga.config import rho_a, rho_w, eta_w
        from math import pi, cos, sin
        from anuga.config import time_format
        from anuga.file_conversion.sts2sww_mesh import sts2sww_mesh

        cellsize = 25
        #nrows=25; ncols = 25;
        nrows=10; ncols = 10;
        refzone=50
        xllcorner=366000;yllcorner=6369500;
        number_of_timesteps = 10
        timestep=1
        eps=2.e-16
        spatial_thinning=1

        points, vertices, boundary =rectangular(nrows-2,ncols-2,
                                                len1=cellsize*(ncols-1),
                                                len2=cellsize*(nrows-1),
                                                origin=(xllcorner,yllcorner))

        time=num.arange(0,10,1,num.float)
        eval_time=time[7];

        domain = Domain(points, vertices, boundary)
        midpoints = domain.get_centroid_coordinates()
        vertexpoints = domain.get_nodes()

        """
        x=grid_1d(xllcorner,cellsize,ncols)
        y=grid_1d(yllcorner,cellsize,nrows)
        X,Y=num.meshgrid(x,y)
        interpolation_points=num.empty((X.shape[0]*X.shape[1],2),num.float)
        k=0
        for i in range(X.shape[0]):
            for j in range(X.shape[1]):
                interpolation_points[k,0]=X[i,j]
                interpolation_points[k,1]=Y[i,j]
                k+=1

        z=spatial_linear_varying_speed(eval_time,interpolation_points[:,0],
                                       interpolation_points[:,1])

        k=0
        Z=num.empty((X.shape[0],X.shape[1]),num.float)
        for i in range(X.shape[0]):
            for j in range(X.shape[1]):
                Z[i,j]=z[k]
                k+=1

        Q=num.empty((time.shape[0],points.shape[0]),num.float)
        for i, t in enumerate(time):
            Q[i,:]=spatial_linear_varying_speed(t,points[:,0],points[:,1])

        from interpolate import Interpolation_function
        I  = Interpolation_function(time,Q,
                                    vertex_coordinates = points,
                                    triangles = domain.triangles,
                                    #interpolation_points = midpoints,
                                    interpolation_points=interpolation_points,
                                    verbose=False)

        V=num.empty((X.shape[0],X.shape[1]),num.float)
        for k in range(len(interpolation_points)):
            assert num.allclose(I(eval_time,k),z[k])
            V[k/X.shape[1],k%X.shape[1]]=I(eval_time,k)


           import mpl_toolkits.mplot3d.axes3d as p3
           fig=P.figure()
           ax = p3.Axes3D(fig)
           ax.plot_surface(X,Y,V)
           ax.plot_surface(X,Y,Z)
           P.show()


        """

        # Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        domain.time = 7*timestep    # Take a time that is represented in file (not zero)

        # Write wind stress file (ensure that domain.time is covered)

        field_sts_filename = 'wind_field'
        self.write_wind_pressure_field_sts(field_sts_filename,
                                      nrows=nrows,
                                      ncols=ncols,
                                      cellsize=cellsize,
                                      origin=(xllcorner,yllcorner),
                                      refzone=50,
                                      timestep=timestep,
                                      number_of_timesteps=10,
                                      speed=spatial_linear_varying_speed,
                                      angle=spatial_linear_varying_angle,
                                      pressure=spatial_linear_varying_pressure)


        sts2sww_mesh(field_sts_filename,spatial_thinning=spatial_thinning,
                     verbose=False)

        # Setup wind stress
        FW = file_function(field_sts_filename+'.sww', domain,
                          quantities=['wind_speed', 'wind_angle'],
                          interpolation_points = midpoints)

        W = Wind_stress(FW,use_coordinates=False)

        domain.forcing_terms = []
        domain.forcing_terms.append(W)

        domain.compute_forcing_terms()

        # Compute reference solution
        const = eta_w*rho_a / rho_w

        N = len(domain)    # number_of_triangles

        xc = domain.get_centroid_coordinates()
        t = domain.time

        x = xc[:,0]
        y = xc[:,1]
        s_vec = spatial_linear_varying_speed(t,x,y)
        phi_vec = spatial_linear_varying_angle(t,x,y)

        for k in range(N):
            # Convert to radians
            phi = phi_vec[k]*pi / 180
            s = s_vec[k]

            # Compute velocity vector (u, v)
            u = s*cos(phi)
            v = s*sin(phi)

            # Compute wind stress
            S = const * num.sqrt(u**2 + v**2)

            assert num.allclose(domain.quantities['stage'].explicit_update[k],0)

            assert num.allclose(domain.quantities['xmomentum'].\
                                    explicit_update[k],S*u,eps)
            assert num.allclose(domain.quantities['ymomentum'].\
                                     explicit_update[k],S*v,eps)

        os.remove(field_sts_filename+'.sts')
        os.remove(field_sts_filename+'.sww')

    def test_variable_pressurefield_from_file(self):
        from anuga.config import rho_a, rho_w, eta_w
        from math import pi, cos, sin
        from anuga.config import time_format
        from anuga.file_conversion.sts2sww_mesh import sts2sww_mesh

        cellsize = 25
        #nrows=25; ncols = 25;
        nrows=10; ncols = 10;
        refzone=50
        xllcorner=366000;yllcorner=6369500;
        number_of_timesteps = 10
        timestep=1
        eps=2.e-16
        spatial_thinning=1

        points, vertices, boundary =rectangular(nrows-2,ncols-2,
                                                len1=cellsize*(ncols-1),
                                                len2=cellsize*(nrows-1),
                                                origin=(xllcorner,yllcorner))

        time=num.arange(0,10,1,num.float)
        eval_time=time[7];

        domain = Domain(points, vertices, boundary)
        midpoints = domain.get_centroid_coordinates()
        vertexpoints = domain.get_nodes()

        # Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        domain.time = 7*timestep    # Take a time that is represented in file (not zero)

        # Write wind stress file (ensure that domain.time is covered)

        field_sts_filename = 'wind_field'
        self.write_wind_pressure_field_sts(field_sts_filename,
                                      nrows=nrows,
                                      ncols=ncols,
                                      cellsize=cellsize,
                                      origin=(xllcorner,yllcorner),
                                      refzone=50,
                                      timestep=timestep,
                                      number_of_timesteps=10,
                                      speed=spatial_linear_varying_speed,
                                      angle=spatial_linear_varying_angle,
                                      pressure=spatial_linear_varying_pressure)


        sts2sww_mesh(field_sts_filename,spatial_thinning=spatial_thinning,
                     verbose=False)

        # Setup barometric pressure
        FP = file_function(field_sts_filename+'.sww', domain,
                           quantities=['barometric_pressure'],
                           interpolation_points = vertexpoints)

        P = Barometric_pressure(FP,use_coordinates=False)


        domain.forcing_terms = []
        domain.forcing_terms.append(P)

        domain.compute_forcing_terms()

        N = len(domain)    # number_of_triangles

        xc = domain.get_centroid_coordinates()
        t = domain.time

        x = xc[:,0]
        y = xc[:,1]
        p_vec = spatial_linear_varying_pressure(t,x,y)

        h=1 #depth
        px=0.000025  #pressure gradient in x-direction
        py=0.0000125 #pressure gradient in y-direction
        for k in range(N):
            # Convert to radians
            p = p_vec[k]

            assert num.allclose(domain.quantities['stage'].explicit_update[k],0)

            assert num.allclose(domain.quantities['xmomentum'].\
                                    explicit_update[k],h*px/rho_w)

            assert num.allclose(domain.quantities['ymomentum'].\
                                     explicit_update[k],h*py/rho_w)

        os.remove(field_sts_filename+'.sts')
        os.remove(field_sts_filename+'.sww')

    def test_constant_wind_stress_from_file_evolve(self):
        from anuga.config import rho_a, rho_w, eta_w
        from math import pi, cos, sin
        from anuga.config import time_format
        from anuga.file_conversion.sts2sww_mesh import sts2sww_mesh

        cellsize = 25
        nrows=5; ncols = 6;
        refzone=50
        xllcorner=366000;yllcorner=6369500;
        number_of_timesteps = 27
        timestep=1
        eps=2e-16

        points, vertices, boundary =rectangular(nrows-2,ncols-2,
                                                len1=cellsize*(ncols-1),
                                                len2=cellsize*(nrows-1),
                                                origin=(xllcorner,yllcorner))

        domain = Domain(points, vertices, boundary)
        midpoints = domain.get_centroid_coordinates()

        # Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'top': Br, 'bottom' :Br, 'left': Br, 'right': Br})

        # Setup only one forcing term, constant wind stress
        s = 100
        phi = 135
        field_sts_filename = 'wind_field'
        self.write_wind_pressure_field_sts(field_sts_filename,
                                      nrows=nrows,
                                      ncols=ncols,
                                      cellsize=cellsize,
                                      origin=(xllcorner,yllcorner),
                                      refzone=50,
                                      timestep=timestep,
                                      number_of_timesteps=number_of_timesteps,
                                      speed=s,
                                      angle=phi)

        sts2sww_mesh(field_sts_filename,spatial_thinning=1,
                     verbose=False)

        # Setup wind stress
        F = file_function(field_sts_filename+'.sww', domain,
                          quantities=['wind_speed', 'wind_angle'],
                          interpolation_points = midpoints)

        W = Wind_stress(F,use_coordinates=False)
        domain.forcing_terms.append(W)

        valuesUsingFunction=num.empty((3,number_of_timesteps+1,midpoints.shape[0]),
                                      num.float)
        i=0
        for t in domain.evolve(yieldstep=1, finaltime=number_of_timesteps*timestep):
            valuesUsingFunction[0,i]=domain.quantities['stage'].explicit_update
            valuesUsingFunction[1,i]=domain.quantities['xmomentum'].explicit_update
            valuesUsingFunction[2,i]=domain.quantities['ymomentum'].explicit_update
            i+=1


        domain_II = Domain(points, vertices, boundary)

        # Flat surface with 1m of water
        domain_II.set_quantity('elevation', 0)
        domain_II.set_quantity('stage', 1.0)
        domain_II.set_quantity('friction', 0)

        Br = Reflective_boundary(domain_II)
        domain_II.set_boundary({'top': Br, 'bottom' :Br, 'left': Br, 'right': Br})

        s = 100
        phi = 135
        domain_II.forcing_terms = []
        domain_II.forcing_terms.append(Wind_stress(s, phi))

        i=0;
        for t in domain_II.evolve(yieldstep=1, 
                                  finaltime=number_of_timesteps*timestep):
            assert num.allclose(valuesUsingFunction[0,i],domain_II.quantities['stage'].explicit_update), max(valuesUsingFunction[0,i]-domain_II.quantities['stage'].explicit_update)
            assert  num.allclose(valuesUsingFunction[1,i],domain_II.quantities['xmomentum'].explicit_update)
            assert num.allclose(valuesUsingFunction[2,i],domain_II.quantities['ymomentum'].explicit_update)
            i+=1

        os.remove(field_sts_filename+'.sts')
        os.remove(field_sts_filename+'.sww')

    def test_temporally_varying_wind_stress_from_file_evolve(self):
        from anuga.config import rho_a, rho_w, eta_w
        from math import pi, cos, sin
        from anuga.config import time_format
        from anuga.file_conversion.sts2sww_mesh import sts2sww_mesh

        cellsize = 25
        #nrows=20; ncols = 20;
        nrows=10; ncols = 10;
        refzone=50
        xllcorner=366000;yllcorner=6369500;
        number_of_timesteps = 28
        timestep=1.
        eps=2e-16

        #points, vertices, boundary =rectangular(10,10,
        points, vertices, boundary =rectangular(5,5,
                                                len1=cellsize*(ncols-1),
                                                len2=cellsize*(nrows-1),
                                                origin=(xllcorner,yllcorner))

        domain = Domain(points, vertices, boundary)
        midpoints = domain.get_centroid_coordinates()

        # Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'top': Br, 'bottom' :Br, 'left': Br, 'right': Br})

        # Setup only one forcing term, constant wind stress
        field_sts_filename = 'wind_field'
        self.write_wind_pressure_field_sts(field_sts_filename,
                                      nrows=nrows,
                                      ncols=ncols,
                                      cellsize=cellsize,
                                      origin=(xllcorner,yllcorner),
                                      refzone=50,
                                      timestep=timestep,
                                      number_of_timesteps=number_of_timesteps,
                                      speed=time_varying_speed,
                                      angle=time_varying_angle,
                                      pressure=time_varying_pressure)

        sts2sww_mesh(field_sts_filename,spatial_thinning=1,
                     verbose=False)

        # Setup wind stress
        F = file_function(field_sts_filename+'.sww', domain,
                          quantities=['wind_speed', 'wind_angle'],
                          interpolation_points = midpoints)

        #W = Wind_stress(F,use_coordinates=False)
        W = Wind_stress_fast(F,filename=field_sts_filename+'.sww', domain=domain)
        domain.forcing_terms.append(W)

        valuesUsingFunction=num.empty((3,2*number_of_timesteps,midpoints.shape[0]),
                                      num.float)
        i=0
        for t in domain.evolve(yieldstep=timestep/2., finaltime=(number_of_timesteps-1)*timestep):
            valuesUsingFunction[0,i]=domain.quantities['stage'].explicit_update
            valuesUsingFunction[1,i]=domain.quantities['xmomentum'].explicit_update
            valuesUsingFunction[2,i]=domain.quantities['ymomentum'].explicit_update
            i+=1


        domain_II = Domain(points, vertices, boundary)

        # Flat surface with 1m of water
        domain_II.set_quantity('elevation', 0)
        domain_II.set_quantity('stage', 1.0)
        domain_II.set_quantity('friction', 0)

        Br = Reflective_boundary(domain_II)
        domain_II.set_boundary({'top': Br, 'bottom' :Br, 'left': Br, 'right': Br})

        domain_II.forcing_terms.append(Wind_stress(s=time_varying_speed, 
                                                   phi=time_varying_angle))

        i=0;
        for t in domain_II.evolve(yieldstep=timestep/2., 
                                  finaltime=(number_of_timesteps-1)*timestep):
            assert num.allclose(valuesUsingFunction[0,i],
                                domain_II.quantities['stage'].explicit_update,
                                eps)
            #print i,valuesUsingFunction[1,i]
            assert  num.allclose(valuesUsingFunction[1,i],
                                 domain_II.quantities['xmomentum'].explicit_update,
                                 eps),(valuesUsingFunction[1,i]-
                                 domain_II.quantities['xmomentum'].explicit_update)
            assert num.allclose(valuesUsingFunction[2,i],
                                domain_II.quantities['ymomentum'].explicit_update,
                                eps)
            #if i==1: assert-1==1
            i+=1

        os.remove(field_sts_filename+'.sts')
        os.remove(field_sts_filename+'.sww')

    def test_spatially_varying_wind_stress_from_file_evolve(self):
        from anuga.config import rho_a, rho_w, eta_w
        from math import pi, cos, sin
        from anuga.config import time_format
        from anuga.file_conversion.sts2sww_mesh import sts2sww_mesh

        cellsize = 25
        nrows=20; ncols = 20;
        nrows=10; ncols = 10;
        refzone=50
        xllcorner=366000;yllcorner=6369500;
        number_of_timesteps = 28
        timestep=1.
        eps=2e-16

        #points, vertices, boundary =rectangular(10,10,
        points, vertices, boundary =rectangular(5,5,
                                                len1=cellsize*(ncols-1),
                                                len2=cellsize*(nrows-1),
                                                origin=(xllcorner,yllcorner))

        domain = Domain(points, vertices, boundary)
        midpoints = domain.get_centroid_coordinates()

        # Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'top': Br, 'bottom' :Br, 'left': Br, 'right': Br})

        # Setup only one forcing term, constant wind stress
        field_sts_filename = 'wind_field'
        self.write_wind_pressure_field_sts(field_sts_filename,
                                      nrows=nrows,
                                      ncols=ncols,
                                      cellsize=cellsize,
                                      origin=(xllcorner,yllcorner),
                                      refzone=50,
                                      timestep=timestep,
                                      number_of_timesteps=number_of_timesteps,
                                      speed=spatial_linear_varying_speed,
                                      angle=spatial_linear_varying_angle,
                                      pressure=spatial_linear_varying_pressure)

        sts2sww_mesh(field_sts_filename,spatial_thinning=1,
                     verbose=False)

        # Setup wind stress
        F = file_function(field_sts_filename+'.sww', domain,
                          quantities=['wind_speed', 'wind_angle'],
                          interpolation_points = midpoints)

        W = Wind_stress(F,use_coordinates=False)
        domain.forcing_terms.append(W)

        valuesUsingFunction=num.empty((3,number_of_timesteps,midpoints.shape[0]),
                                      num.float)
        i=0
        for t in domain.evolve(yieldstep=timestep, finaltime=(number_of_timesteps-1)*timestep):
            valuesUsingFunction[0,i]=domain.quantities['stage'].explicit_update
            valuesUsingFunction[1,i]=domain.quantities['xmomentum'].explicit_update
            valuesUsingFunction[2,i]=domain.quantities['ymomentum'].explicit_update
            i+=1


        domain_II = Domain(points, vertices, boundary)

        # Flat surface with 1m of water
        domain_II.set_quantity('elevation', 0)
        domain_II.set_quantity('stage', 1.0)
        domain_II.set_quantity('friction', 0)

        Br = Reflective_boundary(domain_II)
        domain_II.set_boundary({'top': Br, 'bottom' :Br, 'left': Br, 'right': Br})

        domain_II.forcing_terms.append(Wind_stress(s=spatial_linear_varying_speed, 
                                                   phi=spatial_linear_varying_angle))

        i=0;
        for t in domain_II.evolve(yieldstep=timestep, 
                                  finaltime=(number_of_timesteps-1)*timestep):
            #print valuesUsingFunction[1,i],domain_II.quantities['xmomentum'].explicit_update
            assert num.allclose(valuesUsingFunction[0,i],
                                domain_II.quantities['stage'].explicit_update,
                                eps)
            assert  num.allclose(valuesUsingFunction[1,i],
                                 domain_II.quantities['xmomentum'].explicit_update,
                                 eps)
            assert num.allclose(valuesUsingFunction[2,i],
                                domain_II.quantities['ymomentum'].explicit_update,
                                eps)
            i+=1

        os.remove(field_sts_filename+'.sts')
        os.remove(field_sts_filename+'.sww')

    def test_temporally_varying_pressure_stress_from_file_evolve(self):
        from anuga.config import rho_a, rho_w, eta_w
        from math import pi, cos, sin
        from anuga.config import time_format
        from anuga.file_conversion.sts2sww_mesh import sts2sww_mesh

        cellsize = 25
        #nrows=20; ncols = 20;
        nrows=10; ncols = 10;
        refzone=50
        xllcorner=366000;yllcorner=6369500;
        number_of_timesteps = 28
        timestep=10.
        eps=2e-16

        #print "Building mesh"
        #points, vertices, boundary =rectangular(10,10,
        points, vertices, boundary =rectangular(5,5,
                                                len1=cellsize*(ncols-1),
                                                len2=cellsize*(nrows-1),
                                                origin=(xllcorner,yllcorner))

        domain = Domain(points, vertices, boundary)
        vertexpoints = domain.get_nodes()

        # Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'top': Br, 'bottom' :Br, 'left': Br, 'right': Br})

        # Setup only one forcing term, constant wind stress
        field_sts_filename = 'wind_field'
        #print 'Writing pressure field sts file'
        self.write_wind_pressure_field_sts(field_sts_filename,
                                      nrows=nrows,
                                      ncols=ncols,
                                      cellsize=cellsize,
                                      origin=(xllcorner,yllcorner),
                                      refzone=50,
                                      timestep=timestep,
                                      number_of_timesteps=number_of_timesteps,
                                      speed=time_varying_speed,
                                      angle=time_varying_angle,
                                      pressure=time_varying_pressure)

        #print "converting sts to sww"
        sts2sww_mesh(field_sts_filename,spatial_thinning=1,
                     verbose=False)

        #print 'initialising file_function'
        # Setup wind stress
        F = file_function(field_sts_filename+'.sww', domain,
                          quantities=['barometric_pressure'],
                          interpolation_points = vertexpoints)

        #P = Barometric_pressure(F,use_coordinates=False)
        #print 'initialising pressure forcing term'
        P = Barometric_pressure_fast(p=F,filename=field_sts_filename+'.sww',domain=domain)
        domain.forcing_terms.append(P)

        valuesUsingFunction=num.empty((3,2*number_of_timesteps,len(domain)),
                                      num.float)
        i=0
        import time as timer
        t0=timer.time()
        for t in domain.evolve(yieldstep=timestep/2., finaltime=(number_of_timesteps-1)*timestep):
            valuesUsingFunction[0,i]=domain.quantities['stage'].explicit_update
            valuesUsingFunction[1,i]=domain.quantities['xmomentum'].explicit_update
            valuesUsingFunction[2,i]=domain.quantities['ymomentum'].explicit_update
            i+=1
            #domain.write_time()
        t1=timer.time()
        #print "That took %fs seconds" %(t1-t0)


        domain_II = Domain(points, vertices, boundary)

        # Flat surface with 1m of water
        domain_II.set_quantity('elevation', 0)
        domain_II.set_quantity('stage', 1.0)
        domain_II.set_quantity('friction', 0)

        Br = Reflective_boundary(domain_II)
        domain_II.set_boundary({'top': Br, 'bottom' :Br, 'left': Br, 'right': Br})

        domain_II.forcing_terms.append(Barometric_pressure(p=time_varying_pressure))

        i=0;
        for t in domain_II.evolve(yieldstep=timestep/2., 
                                  finaltime=(number_of_timesteps-1)*timestep):
            assert num.allclose(valuesUsingFunction[0,i],
                                domain_II.quantities['stage'].explicit_update,
                                eps)
            assert  num.allclose(valuesUsingFunction[1,i],
                                 domain_II.quantities['xmomentum'].explicit_update,
                                 eps)
            assert num.allclose(valuesUsingFunction[2,i],
                                domain_II.quantities['ymomentum'].explicit_update,
                                eps)
            i+=1

        os.remove(field_sts_filename+'.sts')
        os.remove(field_sts_filename+'.sww')

    def test_spatially_varying_pressure_stress_from_file_evolve(self):
        from anuga.config import rho_a, rho_w, eta_w
        from math import pi, cos, sin
        from anuga.config import time_format
        from anuga.file_conversion.sts2sww_mesh import sts2sww_mesh

        cellsize = 25
        #nrows=20; ncols = 20;
        nrows=10; ncols = 10;
        refzone=50
        xllcorner=366000;yllcorner=6369500;
        number_of_timesteps = 28
        timestep=1.
        eps=2e-16

        #points, vertices, boundary =rectangular(10,10,
        points, vertices, boundary =rectangular(5,5,
                                                len1=cellsize*(ncols-1),
                                                len2=cellsize*(nrows-1),
                                                origin=(xllcorner,yllcorner))

        domain = Domain(points, vertices, boundary)
        vertexpoints = domain.get_nodes()

        # Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'top': Br, 'bottom' :Br, 'left': Br, 'right': Br})

        # Setup only one forcing term, constant wind stress
        field_sts_filename = 'wind_field'
        self.write_wind_pressure_field_sts(field_sts_filename,
                                      nrows=nrows,
                                      ncols=ncols,
                                      cellsize=cellsize,
                                      origin=(xllcorner,yllcorner),
                                      refzone=50,
                                      timestep=timestep,
                                      number_of_timesteps=number_of_timesteps,
                                      speed=spatial_linear_varying_speed,
                                      angle=spatial_linear_varying_angle,
                                      pressure=spatial_linear_varying_pressure)

        sts2sww_mesh(field_sts_filename,spatial_thinning=1,
                     verbose=False)

        # Setup wind stress
        F = file_function(field_sts_filename+'.sww', domain,
                          quantities=['barometric_pressure'],
                          interpolation_points = vertexpoints)

        P = Barometric_pressure(F,use_coordinates=False)
        domain.forcing_terms.append(P)

        valuesUsingFunction=num.empty((3,number_of_timesteps,len(domain)),
                                      num.float)
        i=0
        for t in domain.evolve(yieldstep=timestep, finaltime=(number_of_timesteps-1)*timestep):
            valuesUsingFunction[0,i]=domain.quantities['stage'].explicit_update
            valuesUsingFunction[1,i]=domain.quantities['xmomentum'].explicit_update
            valuesUsingFunction[2,i]=domain.quantities['ymomentum'].explicit_update
            i+=1


        domain_II = Domain(points, vertices, boundary)

        # Flat surface with 1m of water
        domain_II.set_quantity('elevation', 0)
        domain_II.set_quantity('stage', 1.0)
        domain_II.set_quantity('friction', 0)

        Br = Reflective_boundary(domain_II)
        domain_II.set_boundary({'top': Br, 'bottom' :Br, 'left': Br, 'right': Br})

        domain_II.forcing_terms.append(Barometric_pressure(p=spatial_linear_varying_pressure))

        i=0;
        for t in domain_II.evolve(yieldstep=timestep, 
                                  finaltime=(number_of_timesteps-1)*timestep):

            assert num.allclose(valuesUsingFunction[0,i],
                                domain_II.quantities['stage'].explicit_update,
                                eps)
            assert  num.allclose(valuesUsingFunction[1,i],
                                 domain_II.quantities['xmomentum'].explicit_update,
                                 eps)
            assert num.allclose(valuesUsingFunction[2,i],
                                domain_II.quantities['ymomentum'].explicit_update,
                                eps)
            i+=1

        os.remove(field_sts_filename+'.sts')
        os.remove(field_sts_filename+'.sww')

    def test_flux_gravity(self):
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

        # fluxes and gravity term are now combined. To ensure zero flux on boundary
        # need to set reflective boundaries
        domain.update_boundary()
        domain.compute_fluxes()

        
        assert num.allclose(domain.quantities['stage'].explicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].explicit_update, -g*h*3)
        assert num.allclose(domain.quantities['ymomentum'].explicit_update, 0)



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

        # Use the old function which doesn't take into account the extra
        # wetted area due to slope of bed
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


        # Only manning friction in the forcing terms (gravity now combined with flux calc)
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
        # sqrt(3^2 +4^2) = 5

        S = -g*eta**2 / h**(7.0/3)  * 5

        domain.compute_forcing_terms()

        assert num.allclose(domain.quantities['stage'].semi_implicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].semi_implicit_update,3*S)
        assert num.allclose(domain.quantities['ymomentum'].semi_implicit_update,4*S)


    def test_manning_friction_new(self):
        from anuga.config import g
        import math

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

        # Use the new function which takes into account the extra
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
        S = -g*eta**2 / h**(7.0/3) * math.sqrt(10)

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

        S = -g*eta**2 *5 / h**(7.0/3) * math.sqrt(10.0)

        domain.compute_forcing_terms()

        #print 'S', S
        #print domain.quantities['xmomentum'].semi_implicit_update
        #print domain.quantities['ymomentum'].semi_implicit_update

        assert num.allclose(domain.quantities['stage'].semi_implicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].semi_implicit_update,3*S)
        assert num.allclose(domain.quantities['ymomentum'].semi_implicit_update,4*S)





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
            raise Exception, msg


            
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Forcing, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
