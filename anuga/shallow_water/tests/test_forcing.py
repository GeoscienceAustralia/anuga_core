"""Tests of the environmental forcing operators: wind stress, barometric
pressure and rainfall, including fields read from time-series (.tms) and
spatial (.sww) files, plus the in-Python gravity / Manning-friction terms.

The legacy forcing-function classes (Wind_stress, Rainfall, Inflow,
Barometric_pressure) were removed in 4.1; see UPGRADING_TO_4.0.md.
"""

import unittest
import os
import anuga
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.shallow_water.boundaries import Reflective_boundary
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.file_conversion.file_conversion import timefile2netcdf
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular
from anuga.abstract_2d_finite_volumes.util import file_function
from anuga.file.netcdf import NetCDFFile
from anuga.utilities.numerical_tools import ensure_numeric

from anuga.operators.rate_operators import Rate_operator
from anuga.operators.wind_stress_operator import Wind_stress_operator
from anuga.operators.barometric_pressure import Barometric_pressure_operator

import numpy as num
import warnings




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

    x = num.array(x,float)
    y = num.array(y,float)

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

    x = num.array(x,float)
    y = num.array(y,float)

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

    x = num.array(x,float)
    y = num.array(y,float)

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
    a=0.000025; b=0.0000125
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
    b1=0.000025; b2=0.00001125
    for k in range(N):
        a[k]=phi*(1+t/100.)+b1*x[k]+b2*y[k]
    return a

def spatial_linear_varying_pressure(t, x, y):
    p0=1000
    a=0.000025; b=0.0000125

    x = num.array(x)
    y = num.array(y)

    N = len(x)
    p = 0 * x    # New array

    for k in range(N):
        p[k]=p0*(1.-t/100.)+a*x[k]+b*y[k]
    return p


def grid_1d(x0,dx,nx):
    x = num.empty(nx,dtype=float)
    for i in range(nx):
        x[i]=x0+float(i)*dx
    return x


def ndgrid(x,y):
    nx = len(x)
    ny = len(y)
    X = num.empty(nx*ny,dtype=float)
    Y = num.empty(nx*ny,dtype=float)
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
            except OSError:
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
        starttime = 0; endtime = number_of_timesteps*timestep
        no_data = -9999

        time = num.arange(starttime, endtime, timestep, dtype='i')

        x = grid_1d(xllcorner,cellsize,ncols)
        y = grid_1d(yllcorner,cellsize,nrows)
        [X,Y] = ndgrid(x,y)
        number_of_points = nrows*ncols

        wind_speed = num.empty((number_of_timesteps,nrows*ncols),dtype=float)
        wind_angle = num.empty((number_of_timesteps,nrows*ncols),dtype=float)
        barometric_pressure = num.empty((number_of_timesteps,nrows*ncols),
                                        dtype=float)

        if ( callable(speed) and callable(angle) and callable(pressure) ):
            x = num.ones(3, float)
            y = num.ones(3, float)
            try:
                s = speed(1.0, x=x, y=y)
                a = angle(1.0, x=x, y=y)
                p = pressure(1.0, x=x, y=y)
                use_function=True
            except Exception as e:
                msg = 'Function could not be executed.\n'
                raise(Exception, msg)
        else:
            try :
                speed=float(speed)
                angle=float(angle)
                pressure=float(pressure)
                use_function=False
            except (ValueError, TypeError):
                msg = ('Force fields must be a scalar value coercible to float.')
                raise(Exception, msg)

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

        for i in range(number_of_timesteps):
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


    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _small_domain(self, geo_reference=None, dt=1.0):
        """6-point, 4-triangle flat domain, 1 m of water, timestep *dt*."""
        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]
        #             bac,     bce,     ecf,     dbe
        vertices = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]

        domain = Domain(points, vertices, geo_reference=geo_reference)

        # Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        domain.timestep = dt
        return domain

    @staticmethod
    def _momentum_delta(domain, op):
        """Apply *op* once; return (d_xmom, d_ymom, d_stage) per triangle."""
        xmom0 = domain.quantities['xmomentum'].centroid_values.copy()
        ymom0 = domain.quantities['ymomentum'].centroid_values.copy()
        stage0 = domain.quantities['stage'].centroid_values.copy()
        op()
        return (domain.quantities['xmomentum'].centroid_values - xmom0,
                domain.quantities['ymomentum'].centroid_values - ymom0,
                domain.quantities['stage'].centroid_values - stage0)

    @staticmethod
    def _wind_stress(s, phi):
        """Reference (S*u, S*v) for wind speed *s* [m/s], direction *phi* [deg]."""
        from anuga.config import rho_a, rho_w, eta_w
        from math import pi
        const = eta_w*rho_a/rho_w
        phi = num.asarray(phi, float)*pi/180.0
        u = s*num.cos(phi)
        v = s*num.sin(phi)
        S = const*num.sqrt(u**2 + v**2)
        return S*u, S*v

    def _rectangular_field_domain(self, nrows, ncols, cellsize,
                                  xllcorner, yllcorner):
        points, vertices, boundary = rectangular(nrows-2, ncols-2,
                                                 len1=cellsize*(ncols-1),
                                                 len2=cellsize*(nrows-1),
                                                 origin=(xllcorner, yllcorner))
        domain = Domain(points, vertices, boundary)

        # Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'top': Br, 'bottom': Br, 'left': Br, 'right': Br})
        return domain

    # ------------------------------------------------------------------
    # Wind stress
    # ------------------------------------------------------------------

    def test_constant_wind_stress(self):
        domain = self._small_domain()

        s = 100
        phi = 135
        W = Wind_stress_operator(domain, s, phi)

        dx, dy, dstage = self._momentum_delta(domain, W)
        Su, Sv = self._wind_stress(s, phi)

        assert num.allclose(dstage, 0)
        assert num.allclose(dx, Su)
        assert num.allclose(dy, Sv)

    def test_variable_wind_stress(self):
        domain = self._small_domain()
        domain.set_time(5.54)   # Take a random time (not zero)

        W = Wind_stress_operator(domain, speed=speed, phi=angle)

        dx, dy, dstage = self._momentum_delta(domain, W)

        # Compute reference solution
        xc = domain.get_centroid_coordinates()
        t = domain.get_time()
        Su, Sv = self._wind_stress(speed(t, xc[:,0], xc[:,1]),
                                   angle(t, xc[:,0], xc[:,1]))

        assert num.allclose(dstage, 0)
        assert num.allclose(dx, Su)
        assert num.allclose(dy, Sv)

    def _check_windfield_from_tms(self, dt, time_as_seconds):
        """Uniform, time-varying wind read from a .tms file."""
        import time
        from anuga.config import time_format

        domain = self._small_domain()
        domain.set_time(7)    # Take a time that is represented in file (not zero)

        # Write wind stress file (ensure that domain time is covered)
        # Take x=1 and y=0
        filename = 'test_windstress_from_file'
        start = time.mktime(time.strptime('2000', '%Y'))
        fid = open(filename + '.txt', 'w')
        t = 0.0
        while t <= 10.0:
            if time_as_seconds:
                t_string = str(t)
            else:
                t_string = time.strftime(time_format, time.gmtime(t+start))
            fid.write('%s, %f %f\n' %
                      (t_string, speed(t,[1],[0])[0], angle(t,[1],[0])[0]))
            t += dt
        fid.close()

        timefile2netcdf(filename + '.txt', time_as_seconds=time_as_seconds)
        os.remove(filename + '.txt')

        # Setup wind stress
        F = file_function(filename + '.tms',
                          quantities=['Attribute0', 'Attribute1'])
        os.remove(filename + '.tms')

        W = Wind_stress_operator(domain, F, use_coordinates=False)

        dx, dy, dstage = self._momentum_delta(domain, W)

        # Compute reference solution
        t = domain.get_time()
        Su, Sv = self._wind_stress(speed(t, [1], [0])[0], angle(t, [1], [0])[0])

        assert num.allclose(dstage, 0)
        assert num.allclose(dx, Su)
        assert num.allclose(dy, Sv)

    def test_windfield_from_file(self):
        self._check_windfield_from_tms(dt=1.0, time_as_seconds=False)

    def test_windfield_from_file_seconds(self):
        self._check_windfield_from_tms(dt=0.5, time_as_seconds=True)

    def test_wind_stress_error_condition(self):
        """Wind file must carry exactly (speed, angle); phi must be numeric."""
        domain = self._small_domain()

        with self.assertRaises(ValueError):
            Wind_stress_operator(domain, speed, use_coordinates=False)

        W = Wind_stress_operator(domain, speed=speed, phi='xx')
        with self.assertRaises((ValueError, TypeError)):
            W()

    # ------------------------------------------------------------------
    # Rainfall (Rate_operator.rainfall takes mm/hr; 1 mm/s == 3600 mm/hr)
    # ------------------------------------------------------------------

    MM_S = 3600.0   # multiply a rate in mm/s by this to get mm/hr

    def test_rainfall(self):
        domain = self._small_domain()

        # Constant rainfall of 2 mm/s over the whole domain
        R = Rate_operator.rainfall(domain, rate=2.0*self.MM_S)

        dx, dy, dstage = self._momentum_delta(domain, R)
        assert num.allclose(dstage, 2.0/1000)

    def test_rainfall_restricted_by_polygon(self):
        domain = self._small_domain()

        # Constant rainfall restricted to a polygon enclosing triangle #1 (bce)
        R = Rate_operator.rainfall(domain, rate=2.0*self.MM_S,
                                   polygon=[[1,1], [2,1], [2,2], [1,2]])

        assert num.allclose(R.areas.sum(), 2)

        dx, dy, dstage = self._momentum_delta(domain, R)

        assert num.allclose(dstage[1], 2.0/1000)
        assert num.allclose(dstage[0], 0)
        assert num.allclose(dstage[2:], 0)

    def test_time_dependent_rainfall_restricted_by_polygon(self):
        domain = self._small_domain()

        # Time dependent rainfall restricted to a polygon enclosing
        # triangle #1 (bce)
        R = Rate_operator.rainfall(domain,
                                   rate=lambda t: (3*t + 7)*self.MM_S,
                                   polygon=[[1,1], [2,1], [2,2], [1,2]])

        assert num.allclose(R.areas.sum(), 2)

        domain.set_time(10.0)

        dx, dy, dstage = self._momentum_delta(domain, R)

        assert num.allclose(dstage[1], (3*domain.get_time() + 7)/1000)
        assert num.allclose(dstage[0], 0)
        assert num.allclose(dstage[2:], 0)

    def test_time_dependent_rainfall_using_starttime(self):
        rainfall_poly = ensure_numeric([[1,1], [2,1], [2,2], [1,2]], float)

        domain = self._small_domain()

        R = Rate_operator.rainfall(domain,
                                   rate=lambda t: (3*t + 7)*self.MM_S,
                                   polygon=rainfall_poly)

        assert num.allclose(R.areas.sum(), 2)

        # This will test that time is set to starttime in set_starttime
        domain.set_starttime(5.0)

        dx, dy, dstage = self._momentum_delta(domain, R)

        assert num.allclose(dstage[1], (3*domain.get_time() + 7)/1000)
        assert num.allclose(dstage[1], (3*domain.get_starttime() + 7)/1000)
        assert num.allclose(dstage[0], 0)
        assert num.allclose(dstage[2:], 0)

    def test_time_dependent_rainfall_using_georef(self):
        """Polygon given in absolute (georeferenced) coordinates."""
        # Mesh in zone 56 (absolute coords)
        x0 = 314036.58727982
        y0 = 6224951.2960092

        rainfall_poly = ensure_numeric([[1,1], [2,1], [2,2], [1,2]], float)
        rainfall_poly += [x0, y0]

        domain = self._small_domain(geo_reference=Geo_reference(56, x0, y0))

        R = Rate_operator.rainfall(domain,
                                   rate=lambda t: (3*t + 7)*self.MM_S,
                                   polygon=rainfall_poly)

        assert num.allclose(R.areas.sum(), 2)

        # This will test that time is set to starttime in set_starttime
        domain.set_starttime(5.0)
        domain.set_time(5.0)

        dx, dy, dstage = self._momentum_delta(domain, R)

        assert num.allclose(dstage[1], (3*domain.get_time() + 7)/1000)
        assert num.allclose(dstage[0], 0)
        assert num.allclose(dstage[2:], 0)

    def test_time_dependent_rainfall_restricted_by_polygon_with_default(self):
        """
        Test that default rainfall can be used when given rate runs out of data.
        """
        from anuga.fit_interpolate.interpolate import Modeltime_too_late

        domain = self._small_domain()

        # Time dependent rainfall that expires at t==20
        def main_rate(t):
            if t > 20:
                msg = 'Model time exceeded.'
                raise Modeltime_too_late(msg)
            else:
                return (3*t + 7)*self.MM_S

        R = Rate_operator.rainfall(domain,
                                   rate=main_rate,
                                   polygon=[[1,1], [2,1], [2,2], [1,2]],
                                   default_rate=5.0*self.MM_S)

        assert num.allclose(R.areas.sum(), 2)

        domain.set_time(10.)
        dx, dy, dstage = self._momentum_delta(domain, R)

        assert num.allclose(dstage[1], (3*domain.get_time()+7)/1000)
        assert num.allclose(dstage[0], 0)
        assert num.allclose(dstage[2:], 0)

        domain.set_time(100.)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            dx, dy, dstage = self._momentum_delta(domain, R)

        assert num.allclose(dstage[1], 5.0/1000) # Default value
        assert num.allclose(dstage[0], 0)
        assert num.allclose(dstage[2:], 0)

    def test_rainfall_forcing_with_evolve(self):
        """Rain that runs out of data falls back to default_rate during evolve."""
        from anuga.fit_interpolate.interpolate import Modeltime_too_late

        domain = self._small_domain()

        def main_rate(t):
            if t > 20:
                msg = 'Model time exceeded.'
                raise Modeltime_too_late(msg)
            else:
                return (3*t + 7)*self.MM_S

        Rate_operator.rainfall(domain,
                               rate=main_rate,
                               polygon=[[1,1], [2,1], [2,2], [1,2]],
                               default_rate=5.0*self.MM_S)

        volume0 = domain.compute_total_volume()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            for t in domain.evolve(yieldstep=1, finaltime=25):
                pass

        # Rain fell only on triangle #1 (area 2): 3t+7 mm/s up to t=20,
        # then 5 mm/s (integrate rate over the 2 m^2 area, convert mm to m).
        # The rate is sampled once per (~1 s) step, so allow the resulting
        # quadrature error; the point is that the fallback rate was used
        # rather than an exception raised.
        expected = 2.0*((1.5*20**2 + 7*20) + 5.0*5)/1000
        assert num.allclose(domain.compute_total_volume() - volume0, expected,
                            rtol=2e-2)

    def test_rainfall_forcing_with_evolve_1(self):
        """Without a default_rate, running out of data raises Modeltime_too_late."""
        from anuga.fit_interpolate.interpolate import Modeltime_too_late

        domain = self._small_domain()

        def main_rate(t):
            if t > 20:
                msg = 'Model time exceeded.'
                raise Modeltime_too_late(msg)
            else:
                return (3*t + 7)*self.MM_S

        Rate_operator.rainfall(domain,
                               rate=main_rate,
                               polygon=[[1,1], [2,1], [2,2], [1,2]],
                               default_rate=None)

        with self.assertRaises(Modeltime_too_late):
            for t in domain.evolve(yieldstep=1, finaltime=25):
                pass

    # ------------------------------------------------------------------
    # Wind and pressure fields read from a file (use_coordinates=False)
    # ------------------------------------------------------------------

    def test_constant_wind_stress_from_file(self):
        from anuga.file_conversion.sts2sww_mesh import sts2sww_mesh

        cellsize = 25
        nrows=5; ncols = 6
        xllcorner=366000;yllcorner=6369500
        timestep=12*60

        domain = self._rectangular_field_domain(nrows, ncols, cellsize,
                                                xllcorner, yllcorner)
        domain.timestep = 1.0
        midpoints = domain.get_centroid_coordinates()

        # Constant wind stress
        s = 100
        phi = 135
        pressure=1000

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

        W = Wind_stress_operator(domain, F, use_coordinates=False)

        dx, dy, dstage = self._momentum_delta(domain, W)
        Su, Sv = self._wind_stress(s, phi)

        assert num.allclose(dstage, 0)
        assert num.allclose(dx, Su)
        assert num.allclose(dy, Sv)

        os.remove(field_sts_filename+'.sts')
        os.remove(field_sts_filename+'.sww')

    def test_variable_windfield_from_file(self):
        from anuga.file_conversion.sts2sww_mesh import sts2sww_mesh

        cellsize = 25
        nrows=10; ncols = 10
        xllcorner=366000;yllcorner=6369500
        timestep=1
        eps=2.e-16
        spatial_thinning=1

        domain = self._rectangular_field_domain(nrows, ncols, cellsize,
                                                xllcorner, yllcorner)
        domain.timestep = 1.0
        midpoints = domain.get_centroid_coordinates()

        domain.set_time(7*timestep)    # Take a time that is represented in file (not zero)

        # Write wind stress file (ensure that domain time is covered)
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

        W = Wind_stress_operator(domain, FW, use_coordinates=False)

        dx, dy, dstage = self._momentum_delta(domain, W)

        # Compute reference solution
        xc = domain.get_centroid_coordinates()
        t = domain.get_time()
        Su, Sv = self._wind_stress(spatial_linear_varying_speed(t, xc[:,0], xc[:,1]),
                                   spatial_linear_varying_angle(t, xc[:,0], xc[:,1]))

        assert num.allclose(dstage, 0)
        assert num.allclose(dx, Su, eps)
        assert num.allclose(dy, Sv, eps)

        os.remove(field_sts_filename+'.sts')
        os.remove(field_sts_filename+'.sww')

    def test_variable_pressurefield_from_file(self):
        from anuga.config import rho_w
        from anuga.file_conversion.sts2sww_mesh import sts2sww_mesh

        cellsize = 25
        nrows=10; ncols = 10
        xllcorner=366000;yllcorner=6369500
        timestep=1
        spatial_thinning=1

        domain = self._rectangular_field_domain(nrows, ncols, cellsize,
                                                xllcorner, yllcorner)
        domain.timestep = 1.0
        vertexpoints = domain.get_nodes()

        domain.set_time(7*timestep)    # Take a time that is represented in file (not zero)

        # Write pressure file (ensure that domain time is covered)
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

        P = Barometric_pressure_operator(domain, FP, use_coordinates=False)

        dx, dy, dstage = self._momentum_delta(domain, P)

        h=1 #depth
        px=0.000025  #pressure gradient in x-direction
        py=0.0000125 #pressure gradient in y-direction

        assert num.allclose(dstage, 0)
        assert num.allclose(dx, h*px/rho_w)
        assert num.allclose(dy, h*py/rho_w)

        os.remove(field_sts_filename+'.sts')
        os.remove(field_sts_filename+'.sww')

    def _evolve_momenta(self, domain, yieldstep, finaltime):
        """Evolve, returning the centroid momenta at each yieldstep."""
        xs = []
        ys = []
        for t in domain.evolve(yieldstep=yieldstep, finaltime=finaltime):
            xs.append(domain.quantities['xmomentum'].centroid_values.copy())
            ys.append(domain.quantities['ymomentum'].centroid_values.copy())
        return num.array(xs), num.array(ys)

    def _check_field_from_file_evolve(self, quantity, nrows, ncols,
                                      speed, angle, pressure,
                                      number_of_timesteps, timestep,
                                      yieldstep):
        """A wind or pressure field read from an sww file must drive the
        same evolution as the function it was written from."""
        from anuga.file_conversion.sts2sww_mesh import sts2sww_mesh

        cellsize = 25
        xllcorner=366000;yllcorner=6369500

        domain = self._rectangular_field_domain(nrows, ncols, cellsize,
                                                xllcorner, yllcorner)

        field_sts_filename = 'wind_field'
        self.write_wind_pressure_field_sts(field_sts_filename,
                                      nrows=nrows,
                                      ncols=ncols,
                                      cellsize=cellsize,
                                      origin=(xllcorner,yllcorner),
                                      refzone=50,
                                      timestep=timestep,
                                      number_of_timesteps=number_of_timesteps,
                                      speed=speed,
                                      angle=angle,
                                      pressure=pressure)

        sts2sww_mesh(field_sts_filename,spatial_thinning=1,
                     verbose=False)

        if quantity == 'wind':
            F = file_function(field_sts_filename+'.sww', domain,
                              quantities=['wind_speed', 'wind_angle'],
                              interpolation_points=domain.get_centroid_coordinates())
            Wind_stress_operator(domain, F, use_coordinates=False)
        else:
            F = file_function(field_sts_filename+'.sww', domain,
                              quantities=['barometric_pressure'],
                              interpolation_points=domain.get_nodes())
            Barometric_pressure_operator(domain, F, use_coordinates=False)

        finaltime = (number_of_timesteps-1)*timestep
        xmom_file, ymom_file = self._evolve_momenta(domain, yieldstep, finaltime)

        # Same evolution driven directly by the functions
        domain_II = self._rectangular_field_domain(nrows, ncols, cellsize,
                                                   xllcorner, yllcorner)
        if quantity == 'wind':
            Wind_stress_operator(domain_II, speed, angle)
        else:
            Barometric_pressure_operator(domain_II, pressure)

        xmom_func, ymom_func = self._evolve_momenta(domain_II, yieldstep, finaltime)

        assert xmom_file.shape == xmom_func.shape
        assert num.allclose(xmom_file, xmom_func), \
            num.abs(xmom_file - xmom_func).max()
        assert num.allclose(ymom_file, ymom_func), \
            num.abs(ymom_file - ymom_func).max()
        # ... and the forcing actually did something
        assert num.abs(xmom_func).max() > 0 or num.abs(ymom_func).max() > 0

        os.remove(field_sts_filename+'.sts')
        os.remove(field_sts_filename+'.sww')

    def test_constant_wind_stress_from_file_evolve(self):
        self._check_field_from_file_evolve('wind', nrows=5, ncols=6,
                                           speed=100.0, angle=135.0,
                                           pressure=1000.0,
                                           number_of_timesteps=27, timestep=1,
                                           yieldstep=1)

    def test_temporally_varying_wind_stress_from_file_evolve(self):
        self._check_field_from_file_evolve('wind', nrows=5, ncols=6,
                                           speed=time_varying_speed,
                                           angle=time_varying_angle,
                                           pressure=time_varying_pressure,
                                           number_of_timesteps=28, timestep=1.,
                                           yieldstep=0.5)

    def test_spatially_varying_wind_stress_from_file_evolve(self):
        self._check_field_from_file_evolve('wind', nrows=10, ncols=10,
                                           speed=spatial_linear_varying_speed,
                                           angle=spatial_linear_varying_angle,
                                           pressure=spatial_linear_varying_pressure,
                                           number_of_timesteps=28, timestep=1.,
                                           yieldstep=1)

    def test_temporally_varying_pressure_stress_from_file_evolve(self):
        self._check_field_from_file_evolve('pressure', nrows=5, ncols=6,
                                           speed=time_varying_speed,
                                           angle=time_varying_angle,
                                           pressure=time_varying_pressure,
                                           number_of_timesteps=28, timestep=1.,
                                           yieldstep=0.5)

    def test_spatially_varying_pressure_stress_from_file_evolve(self):
        self._check_field_from_file_evolve('pressure', nrows=10, ncols=10,
                                           speed=spatial_linear_varying_speed,
                                           angle=spatial_linear_varying_angle,
                                           pressure=spatial_linear_varying_pressure,
                                           number_of_timesteps=28, timestep=1.,
                                           yieldstep=1)

    # ------------------------------------------------------------------
    # Gravity and Manning friction (in-Python forcing machinery)
    # ------------------------------------------------------------------

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
        # Mode-2 ('unified') computes fluxes/forcing on-device and never syncs the
        # host explicit_update/semi_implicit_update arrays this white-box test reads.
        # Pin legacy to exercise the in-Python machinery it is written for.
        domain.set_compute_mode('legacy')
        domain.set_flow_algorithm('DE0')

        B = Reflective_boundary(domain)
        domain.set_boundary({'exterior': B})


        # Set up for a gradient of (3,0) at mid triangle (bce)
        def slope(x, y):
            return 3*x

        h = 0.1
        def stage(x, y):
            return slope(x, y) + h

        domain.set_quantity('elevation', slope)
        domain.set_quantity('stage', stage)
        domain.set_quantity('height', h)

        for name in domain.conserved_quantities:
            assert num.allclose(domain.quantities[name].explicit_update, 0)
            assert num.allclose(domain.quantities[name].semi_implicit_update, 0)

        # Fluxes and gravity term are now combined. To ensure zero flux on boundary
        # need to set reflective boundaries
        domain.update_boundary()
        domain.compute_fluxes()

        assert num.allclose(domain.quantities['stage'].explicit_update, 0)

        msg = 'Got %s expected %f' % (domain.quantities['xmomentum'].explicit_update, -g*h*3)
        assert num.allclose(domain.quantities['xmomentum'].explicit_update, -g*h*3), msg
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
        # Mode-2 ('unified') computes friction on-device and never syncs the host
        # semi_implicit_update array this white-box test reads. Pin legacy to
        # exercise the in-Python forcing machinery it is written for.
        domain.set_compute_mode('legacy')

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
        S = -g*eta**2/ h**(7.0/3)

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

        S = -g*eta**2/ h**(7.0/3)  * 5

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
        # Mode-2 ('unified') computes friction on-device and never syncs the host
        # semi_implicit_update array this white-box test reads. Pin legacy to
        # exercise the in-Python forcing machinery it is written for.
        domain.set_compute_mode('legacy')
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
        S = -g*eta**2/ h**(7.0/3) * math.sqrt(10)

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

        S = -g*eta**2 *5/ h**(7.0/3) * math.sqrt(10.0)

        domain.compute_forcing_terms()

        #print 'S', S
        #print domain.quantities['xmomentum'].semi_implicit_update
        #print domain.quantities['ymomentum'].semi_implicit_update

        assert num.allclose(domain.quantities['stage'].semi_implicit_update, 0)
        assert num.allclose(domain.quantities['xmomentum'].semi_implicit_update,3*S)
        assert num.allclose(domain.quantities['ymomentum'].semi_implicit_update,4*S)






if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Forcing)
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
