"""  Test environmental forcing - rain, wind, etc.
"""

import unittest, os
import anuga
import numpy
from anuga import Domain
from anuga import Reflective_boundary
from anuga import rectangular_cross_domain
from anuga import file_function

from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a
from anuga.file_conversion.file_conversion import timefile2netcdf
from anuga.config import time_format

from anuga.fit_interpolate.interpolate import Modeltime_too_early
from anuga.fit_interpolate.interpolate import Modeltime_too_late

from anuga.operators.rate_operators import *

import numpy as num
import warnings
import time
import os


# Setup to skip test if xarray not available
import sys
try:
    import xarray
except ImportError:
    pass

import pytest

warnings.simplefilter("ignore")

verbose = False
class Test_rate_operators(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        try:
            os.remove('test_file_function.txt')
        except:
            pass

        try:
            os.remove('test_file_function.tms')
        except:
            pass


    def test_rate_operator_simple(self):
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


#        print domain.quantities['stage'].centroid_values
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        # Apply operator to these triangles
        indices = [0,1,3]

        rate = 1.0
        factor = 10.0
        default_rate= 0.0

        operator = Rate_operator(domain, rate=rate, factor=factor, \
                      indices=indices, default_rate = default_rate)

        # Apply Operator
        domain.timestep = 2.0
        operator()

        stage_ex = [ 21.,  21.,   1.,  21.]

#        print domain.quantities['stage'].centroid_values
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)
        assert num.allclose(domain.fractional_step_volume_integral, factor*domain.timestep*(rate*domain.areas[indices]).sum())


        # test timestepping_statistics
        stats = operator.timestepping_statistics()
        import re
        rr = re.findall(r"[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", stats)
        assert num.allclose (float(rr[1]), 1.0)
        assert num.allclose (float(rr[2]), 60.0)



    def test_rate_operator_negative_rate(self):
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

#        print domain.quantities['elevation'].centroid_values
#        print domain.quantities['stage'].centroid_values
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        # Apply operator to these triangles
        indices = [0,1,3]



        #Catchment_Rain_Polygon = read_polygon(join('CatchmentBdy.csv'))
        #rainfall = file_function(join('1y120m.tms'), quantities=['rainfall'])
        rate = -1.0
        factor = 10.0
        default_rate= 0.0

        operator = Rate_operator(domain, rate=rate, factor=factor, \
                      indices=indices, default_rate = default_rate)


        # Apply Operator
        domain.timestep = 2.0
        operator()

        stage_ex = [ 0.,  0.,   1.,  0.]
        step_integral = -6.0

        #print domain.quantities['elevation'].centroid_values
        #print domain.quantities['stage'].centroid_values
        #print domain.quantities['xmomentum'].centroid_values
        #print domain.quantities['ymomentum'].centroid_values
        #print domain.fractional_step_volume_integral
        #print factor*domain.timestep*(rate*domain.areas[indices]).sum()

        #increment = factor*domain.timestep*rate*domain.areas



        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)
        assert num.allclose(domain.fractional_step_volume_integral, step_integral)

        # test timestepping_statistics
        stats = operator.timestepping_statistics()
        import re
        rr = re.findall(r"[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", stats)
        assert num.allclose(float(rr[1]), -1.0)
        assert num.allclose(float(rr[2]), -60.0)

    def test_rate_operator_negative_rate_full(self):
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
        domain.set_quantity('stage', 10.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

#        print domain.quantities['elevation'].centroid_values
#        print domain.quantities['stage'].centroid_values
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        # Apply operator to these triangles
        indices = [0,1,3]



        #Catchment_Rain_Polygon = read_polygon(join('CatchmentBdy.csv'))
        #rainfall = file_function(join('1y120m.tms'), quantities=['rainfall'])
        rate = -1.0
        factor = 10.0
        default_rate= 0.0

        operator = Rate_operator(domain, rate=rate, factor=factor, \
                      indices=None, default_rate = default_rate)


        # Apply Operator
        domain.timestep = 2.0
        operator()

        stage_ex = [ 0.,  0.,   0.,  0.]
        step_integral = -80.0

        #print domain.quantities['elevation'].centroid_values
        #print domain.quantities['stage'].centroid_values
        #print domain.quantities['xmomentum'].centroid_values
        #print domain.quantities['ymomentum'].centroid_values
        #print domain.fractional_step_volume_integral


        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)
        assert num.allclose(domain.fractional_step_volume_integral, step_integral)

        # test timestepping_statistics
        stats = operator.timestepping_statistics()
        import re
        rr = re.findall(r"[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", stats)
        assert num.allclose(float(rr[1]), -1.0)
        assert num.allclose(float(rr[2]), -80.0)

    def test_rate_operator_rate_from_file(self):
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


        #---------------------------------
        #Typical ASCII file
        #---------------------------------
        finaltime = 1200
        filename = 'test_file_function'
        fid = open(filename + '.txt', 'w')
        start = time.mktime(time.strptime('2000', '%Y'))
        dt = 60  #One minute intervals
        t = 0.0
        while t <= finaltime:
            t_string = time.strftime(time_format, time.gmtime(t+start))
            fid.write('%s, %f %f %f\n' %(t_string, 2*t, t**2, sin(t*pi/600)))
            t += dt

        fid.close()

        #Convert ASCII file to NetCDF (Which is what we really like!)
        timefile2netcdf(filename+'.txt')


        #Create file function from time series
        F = file_function(filename + '.tms',
                          quantities = ['Attribute0',
                                        'Attribute1',
                                        'Attribute2'])


        #Now try interpolation
        for i in range(20):
            t = i*10
            q = F(t)

            #Exact linear intpolation
            assert num.allclose(q[0], 2*t)
            if i%6 == 0:
                assert num.allclose(q[1], t**2)
                assert num.allclose(q[2], sin(t*pi/600))

        #Check non-exact

        t = 90 #Halfway between 60 and 120
        q = F(t)
        assert num.allclose( (120**2 + 60**2)/2, q[1] )
        assert num.allclose( (sin(120*pi/600) + sin(60*pi/600))/2, q[2] )


        t = 100 #Two thirds of the way between between 60 and 120
        q = F(t)
        assert num.allclose( 2*120**2/3 + 60**2/3, q[1] )
        assert num.allclose( 2*sin(120*pi/600)/3 + sin(60*pi/600)/3, q[2] )

        #os.remove(filename + '.txt')
        #os.remove(filename + '.tms')


        domain = Domain(points, vertices)

        #Flat surface with 1m of water
        domain.set_quantity('elevation', 0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

#        print domain.quantities['elevation'].centroid_values
#        print domain.quantities['stage'].centroid_values
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        # Apply operator to these triangles
        indices = [0,1,3]


        rate = file_function(filename + '.tms', quantities=['Attribute1'])


        # Make starttime of domain consistent with tms file starttime
        domain.set_starttime(rate.starttime)

        factor = 1000.0
        default_rate= 17.7

        operator = Rate_operator(domain, rate=rate, factor=factor, \
                      indices=indices, default_rate = default_rate)


        # Apply Operator
        domain.set_time(360.0)
        domain.timestep = 1.0

        operator()



        d = domain.get_time()**2 * factor + 1.0
        stage_ex0 = [ d,  d,   1.,  d]

#        print d, domain.get_time(), F(360.0)

#        print domain.quantities['elevation'].centroid_values
#        print domain.quantities['stage'].centroid_values
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex0)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)
        assert num.allclose(domain.fractional_step_volume_integral, ((d-1.)*domain.areas[indices]).sum())

        import warnings
        warnings.simplefilter("ignore")

        domain.set_time(1300.0)
        domain.timestep = 1.0

        operator()

        d = default_rate*factor + d
        stage_ex1 = [ d,  d,   1.,  d]

#         print domain.quantities['elevation'].centroid_values
#         print domain.quantities['stage'].centroid_values
#         print domain.quantities['xmomentum'].centroid_values
#         print domain.quantities['ymomentum'].centroid_values

        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex1)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)
        assert num.allclose(domain.fractional_step_volume_integral, ((d-1.)*domain.areas[indices]).sum())


        tmp = numpy.zeros_like(domain.quantities['stage'].centroid_values)
        tmp[:] = domain.quantities['stage'].centroid_values

        d0 = domain.fractional_step_volume_integral

        domain.set_time(-10.0)
        domain.timestep = 1.0

        operator()

        d = default_rate*factor
        stage_ex2 = numpy.array([ d,  d,   0.,  d]) + numpy.array(stage_ex1)

        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex2)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)
        assert num.allclose(domain.fractional_step_volume_integral, d0+(d*domain.areas[indices]).sum())

        # test timestepping_statistics
        stats = operator.timestepping_statistics()
        import re
        rr = re.findall(r"[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", stats)
        assert num.allclose(float(rr[1]), 17.7)
        assert num.allclose(float(rr[2]), 106200.0)

    def test_rate_operator_functions_rate_default_rate(self):
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

        if verbose:
            print(domain.quantities['elevation'].centroid_values)
            print(domain.quantities['stage'].centroid_values)
            print(domain.quantities['xmomentum'].centroid_values)
            print(domain.quantities['ymomentum'].centroid_values)

        # Apply operator to these triangles
        indices = [0,1,3]
        factor = 10.0


        def main_rate(t):
            if t > 20:
                msg = 'Model time exceeded.'
                raise Modeltime_too_late(msg)
            else:
                return 3.0 * t + 7.0

        default_rate = lambda t: 3*t + 7


        operator = Rate_operator(domain, rate=main_rate, factor=factor, \
                      indices=indices, default_rate = default_rate)


        # Apply Operator
        domain.timestep = 2.0
        operator()

        t = operator.get_time()
        d = operator.get_timestep()*main_rate(t)*factor + 1
        stage_ex = [ d,  d,   1.,  d]

        if verbose:
            print("Time ", operator.get_time())
            print("Rate ", main_rate(t))
            print(domain.quantities['elevation'].centroid_values)
            print(domain.quantities['stage'].centroid_values)
            print(domain.quantities['xmomentum'].centroid_values)
            print(domain.quantities['ymomentum'].centroid_values)

        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)
        assert num.allclose(domain.fractional_step_volume_integral, ((d-1.)*domain.areas[indices]).sum())

        domain.set_time(30.0)
        domain.timestep = 1.0
        import warnings
        warnings.simplefilter("ignore")
        operator()

        t = operator.get_time()
        d = operator.get_timestep()*default_rate(t)*factor + d
        stage_ex = [ d,  d,   1.,  d]

        if verbose:
            print("Time ", operator.get_time())
            print("Rate ", default_rate(t))
            print(domain.quantities['elevation'].centroid_values)
            print(domain.quantities['stage'].centroid_values)
            print(domain.quantities['xmomentum'].centroid_values)
            print(domain.quantities['ymomentum'].centroid_values)

        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)


        # test timestepping_statistics
        stats = operator.timestepping_statistics()


        import re
        rr = re.findall(r"[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", stats)

        if verbose:
            print('Operator Statistics: ',stats)
            print('Extracted values: ',rr)
            print('get_Q: ', operator.get_Q())
            print('Get rate value: ', operator.get_non_spatial_rate())
            print('Areas: ', operator.areas)

        assert num.allclose(float(rr[1]), 97.0)
        assert num.allclose(float(rr[2]), 5820.0)


    def test_rate_operator_functions_spatial(self):
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


        area = numpy.sum(domain.areas)

        #Flat surface with 1m of water
        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0.0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        verbose = False

        if verbose:
            print(domain.quantities['elevation'].centroid_values)
            print(domain.quantities['stage'].centroid_values)
            print(domain.quantities['xmomentum'].centroid_values)
            print(domain.quantities['ymomentum'].centroid_values)

        # Apply operator to these triangles
        factor = 10.0


        def main_spatial_rate(x,y,t):
            # x and y should be an n by 1 array
            return x + y

        default_rate = 0.0


        operator = Rate_operator(domain, rate=main_spatial_rate, factor=factor, \
                      default_rate = default_rate)


        # Apply Operator
        domain.timestep = 2.0
        operator()


        t = operator.get_time()
        Q = operator.get_Q()
        x = operator.coord_c[:,0]
        y = operator.coord_c[:,1]
        rate = main_spatial_rate(x,y,t)*factor
        Q_ex = num.sum(domain.areas*rate)
        d = operator.get_timestep()*rate + 1

        #print "d"
        #print d
        #print area, Q, Q_ex
        stage_ex = num.array([ 1.0,  1.0,   1.0,  1.0])
        stage_ex[:] = d

        if verbose:
            print(domain.quantities['elevation'].centroid_values)
            print(domain.quantities['stage'].centroid_values)
            print(domain.quantities['xmomentum'].centroid_values)
            print(domain.quantities['ymomentum'].centroid_values)

        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)
        assert num.allclose(Q_ex, Q)


        # test timestepping_statistics
        stats = operator.timestepping_statistics()
        import re
        rr = re.findall(r"[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", stats)
        assert num.allclose(float(rr[1]), 1.33333)
        assert num.allclose(float(rr[2]), 3.33333)
        assert num.allclose(float(rr[3]), 213.33333)

#operator_5: Min rate = 1.33333 m/s, Max rate = 3.33333 m/s, Total Q = 213.333 m^3/s


    def test_rate_operator_functions_spatial_with_ghost(self):
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


        area = numpy.sum(domain.areas)

        #Flat surface with 1m of water
        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0.0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        verbose = False

        if verbose:
            print(domain.quantities['elevation'].centroid_values)
            print(domain.quantities['stage'].centroid_values)
            print(domain.quantities['xmomentum'].centroid_values)
            print(domain.quantities['ymomentum'].centroid_values)

        # Apply operator to these triangles
        factor = 10.0


        def main_spatial_rate(x,y,t):
            # x and y should be an n by 1 array
            return x + y

        default_rate = 0.0

        # kludge to make a ghost cell
        domain.tri_full_flag[1] = 0

        operator = Rate_operator(domain, rate=main_spatial_rate, factor=factor, \
                      default_rate = default_rate)


        # Apply Operator
        domain.timestep = 2.0
        operator()


        t = operator.get_time()
        Q_all = operator.get_Q(full_only=False)
        Q_full = operator.get_Q()
        x = operator.coord_c[:,0]
        y = operator.coord_c[:,1]
        rate = main_spatial_rate(x,y,t)*factor
        Q_ex_all = num.sum(domain.areas*rate)
        Q_ex_full = num.sum(num.where(domain.tri_full_flag==1,domain.areas*rate,0.0))
        d = operator.get_timestep()*rate + 1

        #print "d"
        #print d
        #print Q_ex_full, Q_ex_all
        stage_ex = num.array([ 1.0,  1.0,   1.0,  1.0])
        stage_ex[:] = d

        if verbose:
            print(domain.quantities['elevation'].centroid_values)
            print(domain.quantities['stage'].centroid_values)
            print(domain.quantities['xmomentum'].centroid_values)
            print(domain.quantities['ymomentum'].centroid_values)

        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)
        assert num.allclose(Q_ex_all, Q_all)
        assert num.allclose(Q_ex_full, Q_full)
        assert num.allclose(domain.fractional_step_volume_integral, ((d-1.)*domain.areas*domain.tri_full_flag).sum())

        # test timestepping_statistics
        stats = operator.timestepping_statistics()
        import re
        rr = re.findall(r"[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", stats)

        assert num.allclose(float(rr[1]), 1.33333)
        assert num.allclose(float(rr[2]), 3.33333)
        assert num.allclose(float(rr[3]), 160.0)

    def test_rate_operator_functions_spatial_indices(self):
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
        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0.0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        verbose = False

        if verbose:
            print(domain.quantities['elevation'].centroid_values)
            print(domain.quantities['stage'].centroid_values)
            print(domain.quantities['xmomentum'].centroid_values)
            print(domain.quantities['ymomentum'].centroid_values)

        # Apply operator to these triangles
        indices = [0,1,3]
        factor = 10.0


        def main_spatial_rate(x,y,t):
            # x and y should be an n by 1 array
            return x + y

        default_rate = 0.0


        operator = Rate_operator(domain, rate=main_spatial_rate, factor=factor, \
                      indices=indices, default_rate = default_rate)


        # Apply Operator
        domain.timestep = 2.0
        operator()

        t = operator.get_time()
        Q = operator.get_Q()
        x = operator.coord_c[indices,0]
        y = operator.coord_c[indices,1]
        rate = main_spatial_rate(x,y,t)*factor
        Q_ex = num.sum(domain.areas[indices]*rate)
        d = operator.get_timestep()*rate + 1


        #print "d"
        #print d
        stage_ex = num.array([ 1.0,  1.0,   1.0,  1.0])
        stage_ex[indices] = d

        if verbose:
            print(domain.quantities['elevation'].centroid_values)
            print(domain.quantities['stage'].centroid_values)
            print(domain.quantities['xmomentum'].centroid_values)
            print(domain.quantities['ymomentum'].centroid_values)

        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)
        assert num.allclose(Q_ex, Q)
        assert num.allclose(domain.fractional_step_volume_integral, ((d-1.)*domain.areas[indices]).sum())


        # test timestepping_statistics
        stats = operator.timestepping_statistics()
        import re
        rr = re.findall(r"[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", stats)

        assert num.allclose(float(rr[1]), 1.33333)
        assert num.allclose(float(rr[2]), 3.33333)
        assert num.allclose(float(rr[3]), 146.667)


    def test_rate_operator_rate_quantity(self):
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
        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0.0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        verbose = False

        if verbose:
            print(domain.quantities['elevation'].centroid_values)
            print(domain.quantities['stage'].centroid_values)
            print(domain.quantities['xmomentum'].centroid_values)
            print(domain.quantities['ymomentum'].centroid_values)

        # Apply operator to these triangles
        indices = [0,1,3]
        factor = 10.0


        from anuga import Quantity
        rate_Q = Quantity(domain)
        rate_Q.set_values(1.0)

        operator = Rate_operator(domain, rate=rate_Q, factor=factor, \
                                 indices=indices)


        # Apply Operator
        domain.timestep = 2.0
        operator()
        rate = rate_Q.centroid_values[indices]
        t = operator.get_time()
        Q = operator.get_Q()

        rate = rate*factor
        Q_ex = num.sum(domain.areas[indices]*rate)
        d = operator.get_timestep()*rate + 1


        #print "d"
        #print d
        #print Q_ex
        #print Q
        stage_ex = num.array([ 1.0,  1.0,   1.0,  1.0])
        stage_ex[indices] = d

        verbose = False

        if verbose:
            print(domain.quantities['elevation'].centroid_values)
            print(domain.quantities['stage'].centroid_values)
            print(domain.quantities['xmomentum'].centroid_values)
            print(domain.quantities['ymomentum'].centroid_values)

        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)
        assert num.allclose(Q_ex, Q)
        assert num.allclose(domain.fractional_step_volume_integral, ((d-1.)*domain.areas[indices]).sum())

        # test timestepping_statistics
        stats = operator.timestepping_statistics()
        import re
        rr = re.findall(r"[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", stats)

        assert num.allclose(float(rr[1]), 1.0)
        assert num.allclose(float(rr[2]), 1.0)
        assert num.allclose(float(rr[3]), 60.0)

    def test_rate_operator_rate_centroid_array(self):
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
        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0.0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        verbose = False

        if verbose:
            print(domain.quantities['elevation'].centroid_values)
            print(domain.quantities['stage'].centroid_values)
            print(domain.quantities['xmomentum'].centroid_values)
            print(domain.quantities['ymomentum'].centroid_values)

        # Apply operator to these triangles
        indices = [0,1,3]
        factor = 10.0


        rate_array = numpy.ones((domain.number_of_triangles,))

        operator = Rate_operator(domain, rate=rate_array, factor=factor, \
                                 indices=indices)


        # Apply Operator
        domain.timestep = 2.0
        operator()
        rate = rate_array[indices]
        t = operator.get_time()
        Q = operator.get_Q()

        rate = rate*factor
        Q_ex = num.sum(domain.areas[indices]*rate)
        d = operator.get_timestep()*rate + 1


        #print "d"
        #print d
        #print Q_ex
        #print Q
        stage_ex = num.array([ 1.0,  1.0,   1.0,  1.0])
        stage_ex[indices] = d

        verbose = False

        if verbose:
            print(domain.quantities['elevation'].centroid_values)
            print(domain.quantities['stage'].centroid_values)
            print(domain.quantities['xmomentum'].centroid_values)
            print(domain.quantities['ymomentum'].centroid_values)

        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)
        assert num.allclose(Q_ex, Q)
        assert num.allclose(domain.fractional_step_volume_integral, ((d-1.)*domain.areas[indices]).sum())

        # test timestepping_statistics
        stats = operator.timestepping_statistics()
        import re
        rr = re.findall(r"[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", stats)

        assert num.allclose(float(rr[1]), 1.0)
        assert num.allclose(float(rr[2]), 1.0)
        assert num.allclose(float(rr[3]), 60.0)
     

    def test_rate_operator_rate_centroid_array_wrong_shape(self):
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
        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0.0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        verbose = False

        if verbose:
            print(domain.quantities['elevation'].centroid_values)
            print(domain.quantities['stage'].centroid_values)
            print(domain.quantities['xmomentum'].centroid_values)
            print(domain.quantities['ymomentum'].centroid_values)

        # Apply operator to these triangles
        indices = [0,1,3]
        factor = 10.0

        # create array with wrong size, should throw an error
        rate_array = numpy.ones((domain.number_of_triangles,2))

        try:
            Rate_operator(domain, rate=rate_array, factor=factor, \
                                 indices=indices)
        except AssertionError: # this is expected
            pass

    @pytest.mark.skipif('xarray' not in sys.modules,
                    reason="requires the xarray module")
    def test_rate_operator_rate_xarray(self):
        from anuga.config import rho_a, rho_w, eta_w
        from math import pi, cos, sin

        import xarray, pandas

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
        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0.0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        verbose = False

        if verbose:
            print(domain.quantities['elevation'].centroid_values)
            print(domain.quantities['stage'].centroid_values)
            print(domain.quantities['xmomentum'].centroid_values)
            print(domain.quantities['ymomentum'].centroid_values)

        # Apply operator to these triangles
        factor = 10.0


        rate_array = numpy.ones((domain.number_of_triangles,))

        operator = Rate_operator(domain, rate=rate_array, factor=factor)


        # Apply Operator
        domain.timestep = 2.0
        operator()
        rate = rate_array
        t = operator.get_time()
        Q = operator.get_Q()

        rate = rate*factor
        Q_ex = num.sum(domain.areas*rate)
        d = operator.get_timestep()*rate + 1


        #print "d"
        #print d
        #print Q_ex
        #print Q
        stage_ex = num.array([ 1.0,  1.0,   1.0,  1.0])
        stage_ex = d

        verbose = False

        if verbose:
            print(domain.quantities['elevation'].centroid_values)
            print(domain.quantities['stage'].centroid_values)
            print(domain.quantities['xmomentum'].centroid_values)
            print(domain.quantities['ymomentum'].centroid_values)

        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)
        assert num.allclose(Q_ex, Q)
        assert num.allclose(domain.fractional_step_volume_integral, ((d-1.)*domain.areas).sum())

        # test timestepping_statistics
        stats = operator.timestepping_statistics()
        import re
        rr = re.findall(r"[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", stats)

        assert num.allclose(float(rr[1]), 1.0)
        assert num.allclose(float(rr[2]), 1.0)
        assert num.allclose(float(rr[3]), 80.0)


    def test_rate_operator_functions_empty_indices(self):
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
        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0.0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        verbose = False

        if verbose:
            print(domain.quantities['elevation'].centroid_values)
            print(domain.quantities['stage'].centroid_values)
            print(domain.quantities['xmomentum'].centroid_values)
            print(domain.quantities['ymomentum'].centroid_values)

        # Apply operator to these triangles
        indices = []
        factor = 10.0


        def main_spatial_rate(x,y,t):
            # x and y should be an n by 1 array
            return x + y

        default_rate = 0.0

        domain.tri_full_flag[0] = 0
        operator = Rate_operator(domain, rate=main_spatial_rate, factor=factor, \
                      indices=indices, default_rate = default_rate)


        # Apply Operator
        domain.timestep = 2.0
        operator()

        t = operator.get_time()
        Q = operator.get_Q()
        x = operator.coord_c[indices,0]
        y = operator.coord_c[indices,1]
        rate = main_spatial_rate(x,y,t)*factor
        Q_ex = num.sum(domain.areas[indices]*rate)
        d = operator.get_timestep()*rate + 1

        # print Q_ex, Q
        # print indices
        # print "d"
        # print d
        stage_ex = num.array([ 1.0,  1.0,   1.0,  1.0])
        stage_ex[indices] = d

        if verbose:
            print(domain.quantities['elevation'].centroid_values)
            print(domain.quantities['stage'].centroid_values)
            print(domain.quantities['xmomentum'].centroid_values)
            print(domain.quantities['ymomentum'].centroid_values)

        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)
        assert num.allclose(Q_ex, Q)
        assert num.allclose(domain.fractional_step_volume_integral, ((d-1.)*domain.areas[indices]).sum())


        # test timestepping_statistics
        stats = operator.timestepping_statistics()
        import re
        rr = re.findall(r"[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", stats)

        assert num.allclose(float(rr[1]), 0.0)
        assert num.allclose(float(rr[2]), 0.0)
        assert num.allclose(float(rr[3]), 0.0)


    def test_rate_operator_functions_empty_region(self):
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
        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0.0)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})

        verbose = False

        if verbose:
            print(domain.quantities['elevation'].centroid_values)
            print(domain.quantities['stage'].centroid_values)
            print(domain.quantities['xmomentum'].centroid_values)
            print(domain.quantities['ymomentum'].centroid_values)

        # Apply operator to these triangles
        indices = []
        region = anuga.Region(domain,indices=indices)

        factor = 10.0


        def main_spatial_rate(x,y,t):
            # x and y should be an n by 1 array
            return x + y

        default_rate = 0.0

        domain.tri_full_flag[0] = 0
        operator = Rate_operator(domain, rate=main_spatial_rate, factor=factor, \
                      region=region, default_rate = default_rate)


        # Apply Operator
        domain.timestep = 2.0
        operator()

        t = operator.get_time()
        Q = operator.get_Q()
        x = operator.coord_c[indices,0]
        y = operator.coord_c[indices,1]
        rate = main_spatial_rate(x,y,t)*factor
        Q_ex = num.sum(domain.areas[indices]*rate)
        d = operator.get_timestep()*rate + 1

        # print Q_ex, Q
        # print indices
        # print "d"
        # print d
        stage_ex = num.array([ 1.0,  1.0,   1.0,  1.0])
        stage_ex[indices] = d

        if verbose:
            print(domain.quantities['elevation'].centroid_values)
            print(domain.quantities['stage'].centroid_values)
            print(domain.quantities['xmomentum'].centroid_values)
            print(domain.quantities['ymomentum'].centroid_values)

        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)
        assert num.allclose(Q_ex, Q)
        assert num.allclose(domain.fractional_step_volume_integral, ((d-1.)*domain.areas[indices]).sum())


        # test timestepping_statistics
        stats = operator.timestepping_statistics()
        import re
        rr = re.findall(r"[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", stats)

        assert num.allclose(float(rr[1]), 0.0)
        assert num.allclose(float(rr[2]), 0.0)
        assert num.allclose(float(rr[3]), 0.0)



if __name__ == "__main__":
    suite = unittest.makeSuite(Test_rate_operators, 'test_')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
