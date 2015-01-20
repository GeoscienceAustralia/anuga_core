"""  Test erosion operators
"""

import unittest, os
import anuga
from anuga import Domain
from anuga import Reflective_boundary
from anuga import rectangular_cross_domain
from anuga import file_function

from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a
from anuga.file_conversion.file_conversion import timefile2netcdf
from anuga.config import time_format

from anuga.operators.erosion_operators import Erosion_operator

import numpy as num
import warnings
import time



class Test_erosion_operators(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass



    def test_erosion_operator_simple(self):
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
        domain.set_quantity('elevation', 0.5)
        domain.set_quantity('stage', 1.0)
        domain.set_quantity('friction', 0)
        domain.set_quantity('xmomentum',2.0)
        domain.set_quantity('ymomentum',3.0)

        Stage = domain.quantities['stage'].centroid_values
        Elevation = domain.quantities['elevation'].centroid_values

        Height = Stage - Elevation

        sum1 = num.sum(Height)

        Br = Reflective_boundary(domain)
        domain.set_boundary({'exterior': Br})


#        print domain.quantities['stage'].centroid_values
#        print domain.quantities['xmomentum'].centroid_values
#        print domain.quantities['ymomentum'].centroid_values

        # Apply operator to these triangles
        indices = [0,1,3]


        operator = Erosion_operator(domain, indices=indices, logging=True)

        # Apply Operator
        domain.timestep = 2.0
        operator()

        elev_ex  = [ 0.05555556,  0.11111111,  0.27777778,  0.05555556]
        stage_ex = [ 0.55555556,  0.61111111,  0.77777778,  0.55555556]

        Stage = domain.quantities['stage'].centroid_values
        Elevation = domain.quantities['elevation'].centroid_values

        Height = Stage - Elevation

        sum2 = num.sum(Height)
        
        #print domain.quantities['elevation'].centroid_values
        #print domain.quantities['stage'].centroid_values
        #print domain.quantities['xmomentum'].centroid_values
        #print domain.quantities['ymomentum'].centroid_values

        assert sum1 == sum2
        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 2.0)
        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 3.0)



#    def test_rate_operator_rate_from_file(self):
#        from anuga.config import rho_a, rho_w, eta_w
#        from math import pi, cos, sin
#
#        a = [0.0, 0.0]
#        b = [0.0, 2.0]
#        c = [2.0, 0.0]
#        d = [0.0, 4.0]
#        e = [2.0, 2.0]
#        f = [4.0, 0.0]
#
#        points = [a, b, c, d, e, f]
#        #             bac,     bce,     ecf,     dbe
#        vertices = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]
#
#
#        #---------------------------------
#        #Typical ASCII file
#        #---------------------------------
#        finaltime = 1200
#        filename = 'test_file_function'
#        fid = open(filename + '.txt', 'w')
#        start = time.mktime(time.strptime('2000', '%Y'))
#        dt = 60  #One minute intervals
#        t = 0.0
#        while t <= finaltime:
#            t_string = time.strftime(time_format, time.gmtime(t+start))
#            fid.write('%s, %f %f %f\n' %(t_string, 2*t, t**2, sin(t*pi/600)))
#            t += dt
#
#        fid.close()
#
#        #Convert ASCII file to NetCDF (Which is what we really like!)
#        timefile2netcdf(filename+'.txt')
#
#
#        #Create file function from time series
#        F = file_function(filename + '.tms',
#                          quantities = ['Attribute0',
#                                        'Attribute1',
#                                        'Attribute2'])
#
#        #Now try interpolation
#        for i in range(20):
#            t = i*10
#            q = F(t)
#
#            #Exact linear intpolation
#            assert num.allclose(q[0], 2*t)
#            if i%6 == 0:
#                assert num.allclose(q[1], t**2)
#                assert num.allclose(q[2], sin(t*pi/600))
#
#        #Check non-exact
#
#        t = 90 #Halfway between 60 and 120
#        q = F(t)
#        assert num.allclose( (120**2 + 60**2)/2, q[1] )
#        assert num.allclose( (sin(120*pi/600) + sin(60*pi/600))/2, q[2] )
#
#
#        t = 100 #Two thirds of the way between between 60 and 120
#        q = F(t)
#        assert num.allclose( 2*120**2/3 + 60**2/3, q[1] )
#        assert num.allclose( 2*sin(120*pi/600)/3 + sin(60*pi/600)/3, q[2] )
#
#        #os.remove(filename + '.txt')
#        #os.remove(filename + '.tms')
#
#
#        domain = Domain(points, vertices)
#
#        #Flat surface with 1m of water
#        domain.set_quantity('elevation', 0)
#        domain.set_quantity('stage', 1.0)
#        domain.set_quantity('friction', 0)
#
#        Br = Reflective_boundary(domain)
#        domain.set_boundary({'exterior': Br})
#
##        print domain.quantities['elevation'].centroid_values
##        print domain.quantities['stage'].centroid_values
##        print domain.quantities['xmomentum'].centroid_values
##        print domain.quantities['ymomentum'].centroid_values
#
#        # Apply operator to these triangles
#        indices = [0,1,3]
#
#
#        rate = file_function('test_file_function.tms', quantities=['Attribute1'])
#
#        factor = 1000.0
#        default_rate= 17.7
#
#        operator = Rate_operator(domain, rate=rate, factor=factor, \
#                      indices=indices, default_rate = default_rate)
#
#
#        # Apply Operator
#        domain.set_starttime(360.0)
#        domain.timestep = 1.0
#
#        operator()
#
#
#        d = domain.get_time()**2 * factor + 1.0
#        stage_ex0 = [ d,  d,   1.,  d]
#
##        print d, domain.get_time(), F(360.0)
#
##        print domain.quantities['elevation'].centroid_values
##        print domain.quantities['stage'].centroid_values
##        print domain.quantities['xmomentum'].centroid_values
##        print domain.quantities['ymomentum'].centroid_values
#
#        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex0)
#        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
#        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)
#
#
#        domain.set_starttime(-10.0)
#        domain.timestep = 1.0
#
#        try:
#            operator()
#        except:
#            pass
#        else:
#            raise Exception('Should have raised an exception, time too early')
#
#
#        domain.set_starttime(1300.0)
#        domain.timestep = 1.0
#
#        operator()
#
#        d = default_rate*factor + d
#        stage_ex1 = [ d,  d,   1.,  d]
#
##        print domain.quantities['elevation'].centroid_values
##        print domain.quantities['stage'].centroid_values
##        print domain.quantities['xmomentum'].centroid_values
##        print domain.quantities['ymomentum'].centroid_values
#
#        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex1)
#        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
#        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)


#    def test_rate_operator_functions_rate_default_rate(self):
#        from anuga.config import rho_a, rho_w, eta_w
#        from math import pi, cos, sin
#
#        a = [0.0, 0.0]
#        b = [0.0, 2.0]
#        c = [2.0, 0.0]
#        d = [0.0, 4.0]
#        e = [2.0, 2.0]
#        f = [4.0, 0.0]
#
#        points = [a, b, c, d, e, f]
#        #             bac,     bce,     ecf,     dbe
#        vertices = [[1,0,2], [1,2,4], [4,2,5], [3,1,4]]
#
#        domain = Domain(points, vertices)
#
#        #Flat surface with 1m of water
#        domain.set_quantity('elevation', 0)
#        domain.set_quantity('stage', 1.0)
#        domain.set_quantity('friction', 0)
#
#        Br = Reflective_boundary(domain)
#        domain.set_boundary({'exterior': Br})
#
#        verbose = False
#
#        if verbose:
#            print domain.quantities['elevation'].centroid_values
#            print domain.quantities['stage'].centroid_values
#            print domain.quantities['xmomentum'].centroid_values
#            print domain.quantities['ymomentum'].centroid_values
#
#        # Apply operator to these triangles
#        indices = [0,1,3]
#        factor = 10.0
#
#
#        def main_rate(t):
#            if t > 20:
#                msg = 'Model time exceeded.'
#                raise Modeltime_too_late, msg
#            else:
#                return 3.0 * t + 7.0
#
#        default_rate = lambda t: 3*t + 7
#
#
#        operator = Rate_operator(domain, rate=main_rate, factor=factor, \
#                      indices=indices, default_rate = default_rate)
#
#
#        # Apply Operator
#        domain.timestep = 2.0
#        operator()
#
#        t = operator.get_time()
#        d = operator.get_timestep()*main_rate(t)*factor + 1
#        stage_ex = [ d,  d,   1.,  d]
#
#        if verbose:
#            print domain.quantities['elevation'].centroid_values
#            print domain.quantities['stage'].centroid_values
#            print domain.quantities['xmomentum'].centroid_values
#            print domain.quantities['ymomentum'].centroid_values
#
#        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
#        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
#        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)
#
#        domain.set_starttime(30.0)
#        domain.timestep = 1.0
#        operator()
#
#        t = operator.get_time()
#        d = operator.get_timestep()*default_rate(t)*factor + d
#        stage_ex = [ d,  d,   1.,  d]
#
#        if verbose:
#            print domain.quantities['elevation'].centroid_values
#            print domain.quantities['stage'].centroid_values
#            print domain.quantities['xmomentum'].centroid_values
#            print domain.quantities['ymomentum'].centroid_values
#
#        assert num.allclose(domain.quantities['stage'].centroid_values, stage_ex)
#        assert num.allclose(domain.quantities['xmomentum'].centroid_values, 0.0)
#        assert num.allclose(domain.quantities['ymomentum'].centroid_values, 0.0)

            
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_erosion_operators, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
