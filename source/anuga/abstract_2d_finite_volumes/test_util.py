#!/usr/bin/env python


import unittest
from Numeric import zeros, array, allclose, Float
from math import sqrt, pi
import tempfile, os
from os import access, F_OK,sep, removedirs,remove,mkdir,getcwd

from anuga.abstract_2d_finite_volumes.util import *
from anuga.config import epsilon
from anuga.shallow_water.data_manager import timefile2netcdf,del_dir

from anuga.utilities.numerical_tools import NAN

from sys import platform 

from anuga.pmesh.mesh import Mesh
from anuga.shallow_water import Domain, Transmissive_boundary
from anuga.shallow_water.data_manager import get_dataobject
from csv import reader,writer
import time
import string

def test_function(x, y):
    return x+y

class Test_Util(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass




    #Geometric
    #def test_distance(self):
    #    from anuga.abstract_2d_finite_volumes.util import distance#
    #
    #    self.failUnless( distance([4,2],[7,6]) == 5.0,
    #                     'Distance is wrong!')
    #    self.failUnless( allclose(distance([7,6],[9,8]), 2.82842712475),
    #                    'distance is wrong!')
    #    self.failUnless( allclose(distance([9,8],[4,2]), 7.81024967591),
    #                    'distance is wrong!')
    #
    #    self.failUnless( distance([9,8],[4,2]) == distance([4,2],[9,8]),
    #                    'distance is wrong!')


    def test_file_function_time1(self):
        """Test that File function interpolates correctly
        between given times. No x,y dependency here.
        """

        #Write file
        import os, time
        from anuga.config import time_format
        from math import sin, pi

        #Typical ASCII file
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
        timefile2netcdf(filename)


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
            assert allclose(q[0], 2*t)
            if i%6 == 0:
                assert allclose(q[1], t**2)
                assert allclose(q[2], sin(t*pi/600))

        #Check non-exact

        t = 90 #Halfway between 60 and 120
        q = F(t)
        assert allclose( (120**2 + 60**2)/2, q[1] )
        assert allclose( (sin(120*pi/600) + sin(60*pi/600))/2, q[2] )


        t = 100 #Two thirds of the way between between 60 and 120
        q = F(t)
        assert allclose( 2*120**2/3 + 60**2/3, q[1] )
        assert allclose( 2*sin(120*pi/600)/3 + sin(60*pi/600)/3, q[2] )

        os.remove(filename + '.txt')
        os.remove(filename + '.tms')        


        
    def test_spatio_temporal_file_function_basic(self):
        """Test that spatio temporal file function performs the correct
        interpolations in both time and space
        NetCDF version (x,y,t dependency)        
        """
        import time

        #Create sww file of simple propagation from left to right
        #through rectangular domain
        from shallow_water import Domain, Dirichlet_boundary
        from mesh_factory import rectangular
        from Numeric import take, concatenate, reshape

        #Create basic mesh and shallow water domain
        points, vertices, boundary = rectangular(3, 3)
        domain1 = Domain(points, vertices, boundary)

        from anuga.utilities.numerical_tools import mean
        domain1.reduction = mean
        domain1.smooth = True #NOTE: Mimic sww output where each vertex has
                              # only one value.

        domain1.default_order = 2
        domain1.store = True
        domain1.set_datadir('.')
        domain1.set_name('spatio_temporal_boundary_source_%d' %(id(self)))
        domain1.quantities_to_be_stored = ['stage', 'xmomentum', 'ymomentum']

        #Bed-slope, friction and IC at vertices (and interpolated elsewhere)
        domain1.set_quantity('elevation', 0)
        domain1.set_quantity('friction', 0)
        domain1.set_quantity('stage', 0)

        # Boundary conditions
        B0 = Dirichlet_boundary([0,0,0])
        B6 = Dirichlet_boundary([0.6,0,0])
        domain1.set_boundary({'left': B6, 'top': B6, 'right': B0, 'bottom': B0})
        domain1.check_integrity()

        finaltime = 8
        #Evolution
        t0 = -1
        for t in domain1.evolve(yieldstep = 0.1, finaltime = finaltime):
            #print 'Timesteps: %.16f, %.16f' %(t0, t)
            #if t == t0:
            #    msg = 'Duplicate timestep found: %f, %f' %(t0, t)
            #   raise msg
            t0 = t
             
            #domain1.write_time()


        #Now read data from sww and check
        from Scientific.IO.NetCDF import NetCDFFile
        filename = domain1.get_name() + '.' + domain1.format
        fid = NetCDFFile(filename)

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        ymomentum = fid.variables['ymomentum'][:]
        time = fid.variables['time'][:]

        #Take stage vertex values at last timestep on diagonal
        #Diagonal is identified by vertices: 0, 5, 10, 15

        last_time_index = len(time)-1 #Last last_time_index
        d_stage = reshape(take(stage[last_time_index, :], [0,5,10,15]), (4,1))
        d_uh = reshape(take(xmomentum[last_time_index, :], [0,5,10,15]), (4,1))
        d_vh = reshape(take(ymomentum[last_time_index, :], [0,5,10,15]), (4,1))
        D = concatenate( (d_stage, d_uh, d_vh), axis=1)

        #Reference interpolated values at midpoints on diagonal at
        #this timestep are
        r0 = (D[0] + D[1])/2
        r1 = (D[1] + D[2])/2
        r2 = (D[2] + D[3])/2

        #And the midpoints are found now
        Dx = take(reshape(x, (16,1)), [0,5,10,15])
        Dy = take(reshape(y, (16,1)), [0,5,10,15])

        diag = concatenate( (Dx, Dy), axis=1)
        d_midpoints = (diag[1:] + diag[:-1])/2

        #Let us see if the file function can find the correct
        #values at the midpoints at the last timestep:
        f = file_function(filename, domain1,
                          interpolation_points = d_midpoints)

        T = f.get_time()
        msg = 'duplicate timesteps: %.16f and %.16f' %(T[-1], T[-2])
        assert not T[-1] == T[-2], msg
        t = time[last_time_index]
        q = f(t, point_id=0); assert allclose(r0, q)
        q = f(t, point_id=1); assert allclose(r1, q)
        q = f(t, point_id=2); assert allclose(r2, q)


        ##################
        #Now do the same for the first timestep

        timestep = 0 #First timestep
        d_stage = reshape(take(stage[timestep, :], [0,5,10,15]), (4,1))
        d_uh = reshape(take(xmomentum[timestep, :], [0,5,10,15]), (4,1))
        d_vh = reshape(take(ymomentum[timestep, :], [0,5,10,15]), (4,1))
        D = concatenate( (d_stage, d_uh, d_vh), axis=1)

        #Reference interpolated values at midpoints on diagonal at
        #this timestep are
        r0 = (D[0] + D[1])/2
        r1 = (D[1] + D[2])/2
        r2 = (D[2] + D[3])/2

        #Let us see if the file function can find the correct
        #values
        q = f(0, point_id=0); assert allclose(r0, q)
        q = f(0, point_id=1); assert allclose(r1, q)
        q = f(0, point_id=2); assert allclose(r2, q)


        ##################
        #Now do it again for a timestep in the middle

        timestep = 33
        d_stage = reshape(take(stage[timestep, :], [0,5,10,15]), (4,1))
        d_uh = reshape(take(xmomentum[timestep, :], [0,5,10,15]), (4,1))
        d_vh = reshape(take(ymomentum[timestep, :], [0,5,10,15]), (4,1))
        D = concatenate( (d_stage, d_uh, d_vh), axis=1)

        #Reference interpolated values at midpoints on diagonal at
        #this timestep are
        r0 = (D[0] + D[1])/2
        r1 = (D[1] + D[2])/2
        r2 = (D[2] + D[3])/2

        q = f(timestep/10., point_id=0); assert allclose(r0, q)
        q = f(timestep/10., point_id=1); assert allclose(r1, q)
        q = f(timestep/10., point_id=2); assert allclose(r2, q)


        ##################
        #Now check temporal interpolation
        #Halfway between timestep 15 and 16

        timestep = 15
        d_stage = reshape(take(stage[timestep, :], [0,5,10,15]), (4,1))
        d_uh = reshape(take(xmomentum[timestep, :], [0,5,10,15]), (4,1))
        d_vh = reshape(take(ymomentum[timestep, :], [0,5,10,15]), (4,1))
        D = concatenate( (d_stage, d_uh, d_vh), axis=1)

        #Reference interpolated values at midpoints on diagonal at
        #this timestep are
        r0_0 = (D[0] + D[1])/2
        r1_0 = (D[1] + D[2])/2
        r2_0 = (D[2] + D[3])/2

        #
        timestep = 16
        d_stage = reshape(take(stage[timestep, :], [0,5,10,15]), (4,1))
        d_uh = reshape(take(xmomentum[timestep, :], [0,5,10,15]), (4,1))
        d_vh = reshape(take(ymomentum[timestep, :], [0,5,10,15]), (4,1))
        D = concatenate( (d_stage, d_uh, d_vh), axis=1)

        #Reference interpolated values at midpoints on diagonal at
        #this timestep are
        r0_1 = (D[0] + D[1])/2
        r1_1 = (D[1] + D[2])/2
        r2_1 = (D[2] + D[3])/2

        # The reference values are
        r0 = (r0_0 + r0_1)/2
        r1 = (r1_0 + r1_1)/2
        r2 = (r2_0 + r2_1)/2

        q = f((timestep - 0.5)/10., point_id=0); assert allclose(r0, q)
        q = f((timestep - 0.5)/10., point_id=1); assert allclose(r1, q)
        q = f((timestep - 0.5)/10., point_id=2); assert allclose(r2, q)

        ##################
        #Finally check interpolation 2 thirds of the way
        #between timestep 15 and 16

        # The reference values are
        r0 = (r0_0 + 2*r0_1)/3
        r1 = (r1_0 + 2*r1_1)/3
        r2 = (r2_0 + 2*r2_1)/3

        #And the file function gives
        q = f((timestep - 1.0/3)/10., point_id=0); assert allclose(r0, q)
        q = f((timestep - 1.0/3)/10., point_id=1); assert allclose(r1, q)
        q = f((timestep - 1.0/3)/10., point_id=2); assert allclose(r2, q)

        fid.close()
        import os
        os.remove(filename)



    def test_spatio_temporal_file_function_different_origin(self):
        """Test that spatio temporal file function performs the correct
        interpolations in both time and space where space is offset by
        xllcorner and yllcorner
        NetCDF version (x,y,t dependency)        
        """
        import time

        #Create sww file of simple propagation from left to right
        #through rectangular domain
        from shallow_water import Domain, Dirichlet_boundary
        from mesh_factory import rectangular
        from Numeric import take, concatenate, reshape


        from anuga.coordinate_transforms.geo_reference import Geo_reference
        xllcorner = 2048
        yllcorner = 11000
        zone = 2

        #Create basic mesh and shallow water domain
        points, vertices, boundary = rectangular(3, 3)
        domain1 = Domain(points, vertices, boundary,
                         geo_reference = Geo_reference(xllcorner = xllcorner,
                                                       yllcorner = yllcorner))
        

        from anuga.utilities.numerical_tools import mean        
        domain1.reduction = mean
        domain1.smooth = True #NOTE: Mimic sww output where each vertex has
                              # only one value.

        domain1.default_order = 2
        domain1.store = True
        domain1.set_datadir('.')
        domain1.set_name('spatio_temporal_boundary_source_%d' %(id(self)))
        domain1.quantities_to_be_stored = ['stage', 'xmomentum', 'ymomentum']

        #Bed-slope, friction and IC at vertices (and interpolated elsewhere)
        domain1.set_quantity('elevation', 0)
        domain1.set_quantity('friction', 0)
        domain1.set_quantity('stage', 0)

        # Boundary conditions
        B0 = Dirichlet_boundary([0,0,0])
        B6 = Dirichlet_boundary([0.6,0,0])
        domain1.set_boundary({'left': B6, 'top': B6, 'right': B0, 'bottom': B0})
        domain1.check_integrity()

        finaltime = 8
        #Evolution
        for t in domain1.evolve(yieldstep = 0.1, finaltime = finaltime):
            pass
            #domain1.write_time()


        #Now read data from sww and check
        from Scientific.IO.NetCDF import NetCDFFile
        filename = domain1.get_name() + '.' + domain1.format
        fid = NetCDFFile(filename)

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        ymomentum = fid.variables['ymomentum'][:]
        time = fid.variables['time'][:]

        #Take stage vertex values at last timestep on diagonal
        #Diagonal is identified by vertices: 0, 5, 10, 15

        last_time_index = len(time)-1 #Last last_time_index     
        d_stage = reshape(take(stage[last_time_index, :], [0,5,10,15]), (4,1))
        d_uh = reshape(take(xmomentum[last_time_index, :], [0,5,10,15]), (4,1))
        d_vh = reshape(take(ymomentum[last_time_index, :], [0,5,10,15]), (4,1))
        D = concatenate( (d_stage, d_uh, d_vh), axis=1)

        #Reference interpolated values at midpoints on diagonal at
        #this timestep are
        r0 = (D[0] + D[1])/2
        r1 = (D[1] + D[2])/2
        r2 = (D[2] + D[3])/2

        #And the midpoints are found now
        Dx = take(reshape(x, (16,1)), [0,5,10,15])
        Dy = take(reshape(y, (16,1)), [0,5,10,15])

        diag = concatenate( (Dx, Dy), axis=1)
        d_midpoints = (diag[1:] + diag[:-1])/2


        #Adjust for georef - make interpolation points absolute
        d_midpoints[:,0] += xllcorner
        d_midpoints[:,1] += yllcorner                

        #Let us see if the file function can find the correct
        #values at the midpoints at the last timestep:
        f = file_function(filename, domain1,
                          interpolation_points = d_midpoints)

        t = time[last_time_index]                         
        q = f(t, point_id=0); assert allclose(r0, q)
        q = f(t, point_id=1); assert allclose(r1, q)
        q = f(t, point_id=2); assert allclose(r2, q)


        ##################
        #Now do the same for the first timestep

        timestep = 0 #First timestep
        d_stage = reshape(take(stage[timestep, :], [0,5,10,15]), (4,1))
        d_uh = reshape(take(xmomentum[timestep, :], [0,5,10,15]), (4,1))
        d_vh = reshape(take(ymomentum[timestep, :], [0,5,10,15]), (4,1))
        D = concatenate( (d_stage, d_uh, d_vh), axis=1)

        #Reference interpolated values at midpoints on diagonal at
        #this timestep are
        r0 = (D[0] + D[1])/2
        r1 = (D[1] + D[2])/2
        r2 = (D[2] + D[3])/2

        #Let us see if the file function can find the correct
        #values
        q = f(0, point_id=0); assert allclose(r0, q)
        q = f(0, point_id=1); assert allclose(r1, q)
        q = f(0, point_id=2); assert allclose(r2, q)


        ##################
        #Now do it again for a timestep in the middle

        timestep = 33
        d_stage = reshape(take(stage[timestep, :], [0,5,10,15]), (4,1))
        d_uh = reshape(take(xmomentum[timestep, :], [0,5,10,15]), (4,1))
        d_vh = reshape(take(ymomentum[timestep, :], [0,5,10,15]), (4,1))
        D = concatenate( (d_stage, d_uh, d_vh), axis=1)

        #Reference interpolated values at midpoints on diagonal at
        #this timestep are
        r0 = (D[0] + D[1])/2
        r1 = (D[1] + D[2])/2
        r2 = (D[2] + D[3])/2

        q = f(timestep/10., point_id=0); assert allclose(r0, q)
        q = f(timestep/10., point_id=1); assert allclose(r1, q)
        q = f(timestep/10., point_id=2); assert allclose(r2, q)


        ##################
        #Now check temporal interpolation
        #Halfway between timestep 15 and 16

        timestep = 15
        d_stage = reshape(take(stage[timestep, :], [0,5,10,15]), (4,1))
        d_uh = reshape(take(xmomentum[timestep, :], [0,5,10,15]), (4,1))
        d_vh = reshape(take(ymomentum[timestep, :], [0,5,10,15]), (4,1))
        D = concatenate( (d_stage, d_uh, d_vh), axis=1)

        #Reference interpolated values at midpoints on diagonal at
        #this timestep are
        r0_0 = (D[0] + D[1])/2
        r1_0 = (D[1] + D[2])/2
        r2_0 = (D[2] + D[3])/2

        #
        timestep = 16
        d_stage = reshape(take(stage[timestep, :], [0,5,10,15]), (4,1))
        d_uh = reshape(take(xmomentum[timestep, :], [0,5,10,15]), (4,1))
        d_vh = reshape(take(ymomentum[timestep, :], [0,5,10,15]), (4,1))
        D = concatenate( (d_stage, d_uh, d_vh), axis=1)

        #Reference interpolated values at midpoints on diagonal at
        #this timestep are
        r0_1 = (D[0] + D[1])/2
        r1_1 = (D[1] + D[2])/2
        r2_1 = (D[2] + D[3])/2

        # The reference values are
        r0 = (r0_0 + r0_1)/2
        r1 = (r1_0 + r1_1)/2
        r2 = (r2_0 + r2_1)/2

        q = f((timestep - 0.5)/10., point_id=0); assert allclose(r0, q)
        q = f((timestep - 0.5)/10., point_id=1); assert allclose(r1, q)
        q = f((timestep - 0.5)/10., point_id=2); assert allclose(r2, q)

        ##################
        #Finally check interpolation 2 thirds of the way
        #between timestep 15 and 16

        # The reference values are
        r0 = (r0_0 + 2*r0_1)/3
        r1 = (r1_0 + 2*r1_1)/3
        r2 = (r2_0 + 2*r2_1)/3

        #And the file function gives
        q = f((timestep - 1.0/3)/10., point_id=0); assert allclose(r0, q)
        q = f((timestep - 1.0/3)/10., point_id=1); assert allclose(r1, q)
        q = f((timestep - 1.0/3)/10., point_id=2); assert allclose(r2, q)

        fid.close()
        import os
        os.remove(filename)

        


    def qtest_spatio_temporal_file_function_time(self):
        """Test that File function interpolates correctly
        between given times.
        NetCDF version (x,y,t dependency)
        """

        #Create NetCDF (sww) file to be read
        # x: 0, 5, 10, 15
        # y: -20, -10, 0, 10
        # t: 0, 60, 120, ...., 1200
        #
        # test quantities (arbitrary but non-trivial expressions):
        #
        #   stage     = 3*x - y**2 + 2*t
        #   xmomentum = exp( -((x-7)**2 + (y+5)**2)/20 ) * t**2
        #   ymomentum = x**2 + y**2 * sin(t*pi/600)

        #NOTE: Nice test that may render some of the others redundant.

        import os, time
        from anuga.config import time_format
        from Numeric import sin, pi, exp
        from mesh_factory import rectangular
        from shallow_water import Domain
        import anuga.shallow_water.data_manager

        finaltime = 1200
        filename = 'test_file_function'

        #Create a domain to hold test grid
        #(0:15, -20:10)
        points, vertices, boundary =\
                rectangular(4, 4, 15, 30, origin = (0, -20))
        print "points", points

        #print 'Number of elements', len(vertices)
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.default_order = 2
        domain.set_datadir('.')
        domain.set_name(filename)
        domain.store = True
        domain.format = 'sww'   #Native netcdf visualisation format

        #print points
        start = time.mktime(time.strptime('2000', '%Y'))
        domain.starttime = start


        #Store structure
        domain.initialise_storage()

        #Compute artificial time steps and store
        dt = 60  #One minute intervals
        t = 0.0
        while t <= finaltime:
            #Compute quantities
            f1 = lambda x,y: 3*x - y**2 + 2*t + 4
            domain.set_quantity('stage', f1)

            f2 = lambda x,y: x+y+t**2
            domain.set_quantity('xmomentum', f2)

            f3 = lambda x,y: x**2 + y**2 * sin(t*pi/600)
            domain.set_quantity('ymomentum', f3)

            #Store and advance time
            domain.time = t
            domain.store_timestep(domain.conserved_quantities)
            t += dt


        interpolation_points = [[0,-20], [1,0], [0,1], [1.1, 3.14], [10,-12.5]]
      
        #Deliberately set domain.starttime to too early
        domain.starttime = start - 1

        #Create file function
        F = file_function(filename + '.sww', domain,
                          quantities = domain.conserved_quantities,
                          interpolation_points = interpolation_points)

        #Check that FF updates fixes domain starttime
        assert allclose(domain.starttime, start)

        #Check that domain.starttime isn't updated if later
        domain.starttime = start + 1
        F = file_function(filename + '.sww', domain,
                          quantities = domain.conserved_quantities,
                          interpolation_points = interpolation_points)
        assert allclose(domain.starttime, start+1)
        domain.starttime = start


        #Check linear interpolation in time
        F = file_function(filename + '.sww', domain,
                          quantities = domain.conserved_quantities,
                          interpolation_points = interpolation_points)                
        for id in range(len(interpolation_points)):
            x = interpolation_points[id][0]
            y = interpolation_points[id][1]

            for i in range(20):
                t = i*10
                k = i%6

                if k == 0:
                    q0 = F(t, point_id=id)
                    q1 = F(t+60, point_id=id)

                if q0 == NAN:
                    actual = q0
                else:
                    actual = (k*q1 + (6-k)*q0)/6
                q = F(t, point_id=id)
                #print i, k, t, q
                #print ' ', q0
                #print ' ', q1
                print "q",q
                print "actual", actual
                #print
                if q0 == NAN:
                     self.failUnless( q == actual, 'Fail!')
                else:
                    assert allclose(q, actual)


        #Another check of linear interpolation in time
        for id in range(len(interpolation_points)):
            q60 = F(60, point_id=id)
            q120 = F(120, point_id=id)

            t = 90 #Halfway between 60 and 120
            q = F(t, point_id=id)
            assert allclose( (q120+q60)/2, q )

            t = 100 #Two thirds of the way between between 60 and 120
            q = F(t, point_id=id)
            assert allclose(q60/3 + 2*q120/3, q)



        #Check that domain.starttime isn't updated if later than file starttime but earlier
        #than file end time
        delta = 23
        domain.starttime = start + delta
        F = file_function(filename + '.sww', domain,
                          quantities = domain.conserved_quantities,
                          interpolation_points = interpolation_points)
        assert allclose(domain.starttime, start+delta)




        #Now try interpolation with delta offset
        for id in range(len(interpolation_points)):            
            x = interpolation_points[id][0]
            y = interpolation_points[id][1]

            for i in range(20):
                t = i*10
                k = i%6

                if k == 0:
                    q0 = F(t-delta, point_id=id)
                    q1 = F(t+60-delta, point_id=id)

                q = F(t-delta, point_id=id)
                assert allclose(q, (k*q1 + (6-k)*q0)/6)


        os.remove(filename + '.sww')



    def NOtest_spatio_temporal_file_function_time(self):
        # FIXME: This passes but needs some TLC
        # Test that File function interpolates correctly
        # When some points are outside the mesh

        import os, time
        from anuga.config import time_format
        from Numeric import sin, pi, exp
        from mesh_factory import rectangular
        from shallow_water import Domain
        import anuga.shallow_water.data_manager 
        from anuga.pmesh.mesh_interface import create_mesh_from_regions
        finaltime = 1200
        
        filename = tempfile.mktemp()
        #print "filename",filename 
        filename = 'test_file_function'

        meshfilename = tempfile.mktemp(".tsh")

        boundary_tags = {'walls':[0,1],'bom':[2]}
        
        polygon_absolute = [[0,-20],[10,-20],[10,15],[-20,15]]
        
        create_mesh_from_regions(polygon_absolute,
                                 boundary_tags,
                                 10000000,
                                 filename=meshfilename)
        domain = Domain(mesh_filename=meshfilename)
        domain.smooth = False
        domain.default_order = 2
        domain.set_datadir('.')
        domain.set_name(filename)
        domain.store = True
        domain.format = 'sww'   #Native netcdf visualisation format

        #print points
        start = time.mktime(time.strptime('2000', '%Y'))
        domain.starttime = start
        

        #Store structure
        domain.initialise_storage()

        #Compute artificial time steps and store
        dt = 60  #One minute intervals
        t = 0.0
        while t <= finaltime:
            #Compute quantities
            f1 = lambda x,y: 3*x - y**2 + 2*t + 4
            domain.set_quantity('stage', f1)

            f2 = lambda x,y: x+y+t**2
            domain.set_quantity('xmomentum', f2)

            f3 = lambda x,y: x**2 + y**2 * sin(t*pi/600)
            domain.set_quantity('ymomentum', f3)

            #Store and advance time
            domain.time = t
            domain.store_timestep(domain.conserved_quantities)
            t += dt

        interpolation_points = [[1,0]]
        interpolation_points = [[100,1000]]
        
        interpolation_points = [[0,-20], [1,0], [0,1], [1.1, 3.14], [10,-12.5],
                                [78787,78787],[7878,3432]]
            
        #Deliberately set domain.starttime to too early
        domain.starttime = start - 1

        #Create file function
        F = file_function(filename + '.sww', domain,
                          quantities = domain.conserved_quantities,
                          interpolation_points = interpolation_points)

        #Check that FF updates fixes domain starttime
        assert allclose(domain.starttime, start)

        #Check that domain.starttime isn't updated if later
        domain.starttime = start + 1
        F = file_function(filename + '.sww', domain,
                          quantities = domain.conserved_quantities,
                          interpolation_points = interpolation_points)
        assert allclose(domain.starttime, start+1)
        domain.starttime = start


        #Check linear interpolation in time
        # checking points inside and outside the mesh
        F = file_function(filename + '.sww', domain,
                          quantities = domain.conserved_quantities,
                          interpolation_points = interpolation_points)
        
        for id in range(len(interpolation_points)):
            x = interpolation_points[id][0]
            y = interpolation_points[id][1]

            for i in range(20):
                t = i*10
                k = i%6

                if k == 0:
                    q0 = F(t, point_id=id)
                    q1 = F(t+60, point_id=id)

                if q0 == NAN:
                    actual = q0
                else:
                    actual = (k*q1 + (6-k)*q0)/6
                q = F(t, point_id=id)
                #print i, k, t, q
                #print ' ', q0
                #print ' ', q1
                #print "q",q
                #print "actual", actual
                #print
                if q0 == NAN:
                     self.failUnless( q == actual, 'Fail!')
                else:
                    assert allclose(q, actual)

        # now lets check points inside the mesh
        interpolation_points = [[0,-20], [1,0], [0,1], [1.1, 3.14]] #, [10,-12.5]] - this point doesn't work WHY?
        interpolation_points = [[10,-12.5]]
            
        print "len(interpolation_points)",len(interpolation_points) 
        F = file_function(filename + '.sww', domain,
                          quantities = domain.conserved_quantities,
                          interpolation_points = interpolation_points)

        domain.starttime = start


        #Check linear interpolation in time
        F = file_function(filename + '.sww', domain,
                          quantities = domain.conserved_quantities,
                          interpolation_points = interpolation_points)                
        for id in range(len(interpolation_points)):
            x = interpolation_points[id][0]
            y = interpolation_points[id][1]

            for i in range(20):
                t = i*10
                k = i%6

                if k == 0:
                    q0 = F(t, point_id=id)
                    q1 = F(t+60, point_id=id)

                if q0 == NAN:
                    actual = q0
                else:
                    actual = (k*q1 + (6-k)*q0)/6
                q = F(t, point_id=id)
                print "############"
                print "id, x, y ", id, x, y #k, t, q
                print "t", t
                #print ' ', q0
                #print ' ', q1
                print "q",q
                print "actual", actual
                #print
                if q0 == NAN:
                     self.failUnless( q == actual, 'Fail!')
                else:
                    assert allclose(q, actual)


        #Another check of linear interpolation in time
        for id in range(len(interpolation_points)):
            q60 = F(60, point_id=id)
            q120 = F(120, point_id=id)

            t = 90 #Halfway between 60 and 120
            q = F(t, point_id=id)
            assert allclose( (q120+q60)/2, q )

            t = 100 #Two thirds of the way between between 60 and 120
            q = F(t, point_id=id)
            assert allclose(q60/3 + 2*q120/3, q)



        #Check that domain.starttime isn't updated if later than file starttime but earlier
        #than file end time
        delta = 23
        domain.starttime = start + delta
        F = file_function(filename + '.sww', domain,
                          quantities = domain.conserved_quantities,
                          interpolation_points = interpolation_points)
        assert allclose(domain.starttime, start+delta)




        #Now try interpolation with delta offset
        for id in range(len(interpolation_points)):            
            x = interpolation_points[id][0]
            y = interpolation_points[id][1]

            for i in range(20):
                t = i*10
                k = i%6

                if k == 0:
                    q0 = F(t-delta, point_id=id)
                    q1 = F(t+60-delta, point_id=id)

                q = F(t-delta, point_id=id)
                assert allclose(q, (k*q1 + (6-k)*q0)/6)


        os.remove(filename + '.sww')

    def test_file_function_time_with_domain(self):
        """Test that File function interpolates correctly
        between given times. No x,y dependency here.
        Use domain with starttime
        """

        #Write file
        import os, time, calendar
        from anuga.config import time_format
        from math import sin, pi
        from domain import Domain

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
        timefile2netcdf(filename)



        a = [0.0, 0.0]
        b = [4.0, 0.0]
        c = [0.0, 3.0]

        points = [a, b, c]
        vertices = [[0,1,2]]
        domain = Domain(points, vertices)

        #Check that domain.starttime is updated if non-existing
        F = file_function(filename + '.tms', domain)
        
        assert allclose(domain.starttime, start)

        #Check that domain.starttime is updated if too early
        domain.starttime = start - 1
        F = file_function(filename + '.tms', domain)
        assert allclose(domain.starttime, start)

        #Check that domain.starttime isn't updated if later
        domain.starttime = start + 1
        F = file_function(filename + '.tms', domain)        
        assert allclose(domain.starttime, start+1)

        domain.starttime = start
        F = file_function(filename + '.tms', domain,
                          quantities = ['Attribute0', 'Attribute1', 'Attribute2'],
                          use_cache=True)
        

        #print F.precomputed_values
        #print 'F(60)', F(60)
        
        #Now try interpolation
        for i in range(20):
            t = i*10
            q = F(t)

            #Exact linear intpolation
            assert allclose(q[0], 2*t)
            if i%6 == 0:
                assert allclose(q[1], t**2)
                assert allclose(q[2], sin(t*pi/600))

        #Check non-exact

        t = 90 #Halfway between 60 and 120
        q = F(t)
        assert allclose( (120**2 + 60**2)/2, q[1] )
        assert allclose( (sin(120*pi/600) + sin(60*pi/600))/2, q[2] )


        t = 100 #Two thirds of the way between between 60 and 120
        q = F(t)
        assert allclose( 2*120**2/3 + 60**2/3, q[1] )
        assert allclose( 2*sin(120*pi/600)/3 + sin(60*pi/600)/3, q[2] )

        os.remove(filename + '.tms')
        os.remove(filename + '.txt')        

    def test_file_function_time_with_domain_different_start(self):
        """Test that File function interpolates correctly
        between given times. No x,y dependency here.
        Use domain with a starttime later than that of file

        ASCII version
        """

        #Write file
        import os, time, calendar
        from anuga.config import time_format
        from math import sin, pi
        from domain import Domain

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
        timefile2netcdf(filename)        

        a = [0.0, 0.0]
        b = [4.0, 0.0]
        c = [0.0, 3.0]

        points = [a, b, c]
        vertices = [[0,1,2]]
        domain = Domain(points, vertices)

        #Check that domain.starttime isn't updated if later than file starttime but earlier
        #than file end time
        delta = 23
        domain.starttime = start + delta
        F = file_function(filename + '.tms', domain,
                          quantities = ['Attribute0', 'Attribute1', 'Attribute2'])        
        assert allclose(domain.starttime, start+delta)




        #Now try interpolation with delta offset
        for i in range(20):
            t = i*10
            q = F(t-delta)

            #Exact linear intpolation
            assert allclose(q[0], 2*t)
            if i%6 == 0:
                assert allclose(q[1], t**2)
                assert allclose(q[2], sin(t*pi/600))

        #Check non-exact

        t = 90 #Halfway between 60 and 120
        q = F(t-delta)
        assert allclose( (120**2 + 60**2)/2, q[1] )
        assert allclose( (sin(120*pi/600) + sin(60*pi/600))/2, q[2] )


        t = 100 #Two thirds of the way between between 60 and 120
        q = F(t-delta)
        assert allclose( 2*120**2/3 + 60**2/3, q[1] )
        assert allclose( 2*sin(120*pi/600)/3 + sin(60*pi/600)/3, q[2] )


        os.remove(filename + '.tms')
        os.remove(filename + '.txt')                



    def test_apply_expression_to_dictionary(self):

        #FIXME: Division is not expected to work for integers.
        #This must be caught.
        foo = array([[1,2,3],
                     [4,5,6]], Float)

        bar = array([[-1,0,5],
                     [6,1,1]], Float)                  

        D = {'X': foo, 'Y': bar}

        Z = apply_expression_to_dictionary('X+Y', D)        
        assert allclose(Z, foo+bar)

        Z = apply_expression_to_dictionary('X*Y', D)        
        assert allclose(Z, foo*bar)        

        Z = apply_expression_to_dictionary('4*X+Y', D)        
        assert allclose(Z, 4*foo+bar)        

        # test zero division is OK
        Z = apply_expression_to_dictionary('X/Y', D)
        assert allclose(1/Z, 1/(foo/bar)) # can't compare inf to inf

        # make an error for zero on zero
        # this is really an error in Numeric, SciPy core can handle it
        # Z = apply_expression_to_dictionary('0/Y', D)

        #Check exceptions
        try:
            #Wrong name
            Z = apply_expression_to_dictionary('4*X+A', D)        
        except NameError:
            pass
        else:
            msg = 'Should have raised a NameError Exception'
            raise msg


        try:
            #Wrong order
            Z = apply_expression_to_dictionary(D, '4*X+A')        
        except AssertionError:
            pass
        else:
            msg = 'Should have raised a AssertionError Exception'
            raise msg        
        

    def test_multiple_replace(self):
        """Hard test that checks a true word-by-word simultaneous replace
        """
        
        D = {'x': 'xi', 'y': 'eta', 'xi':'lam'}
        exp = '3*x+y + xi'
        
        new = multiple_replace(exp, D)
        
        assert new == '3*xi+eta + lam'
                          


    def test_point_on_line_obsolete(self):
        """Test that obsolete call issues appropriate warning"""

        #Turn warning into an exception
        import warnings
        warnings.filterwarnings('error')

        try:
            assert point_on_line( 0, 0.5, 0,1, 0,0 )
        except DeprecationWarning:
            pass
        else:
            msg = 'point_on_line should have issued a DeprecationWarning'
            raise Exception(msg)    

        warnings.resetwarnings()
    
    def test_get_revision_number(self):
        """test_get_revision_number(self):

        Test that revision number can be retrieved.
        """
        if os.environ.has_key('USER') and os.environ['USER'] == 'dgray':
            # I have a known snv incompatability issue,
            # so I'm skipping this test.
            # FIXME when SVN is upgraded on our clusters
            pass
        else:    
            n = get_revision_number()
            assert n>=0


        
    def test_add_directories(self):
        
        import tempfile
        root_dir = tempfile.mkdtemp('_test_util', 'test_util_')
        directories = ['ja','ne','ke']
        kens_dir = add_directories(root_dir, directories)
        assert kens_dir == root_dir + sep + 'ja' + sep + 'ne' + \
               sep + 'ke'
        assert access(root_dir,F_OK)

        add_directories(root_dir, directories)
        assert access(root_dir,F_OK)
        
        #clean up!
        os.rmdir(kens_dir)
        os.rmdir(root_dir + sep + 'ja' + sep + 'ne')
        os.rmdir(root_dir + sep + 'ja')
        os.rmdir(root_dir)

    def test_add_directories_bad(self):
        
        import tempfile
        root_dir = tempfile.mkdtemp('_test_util', 'test_util_')
        directories = ['/\/!@#@#$%^%&*((*:*:','ne','ke']
        
        try:
            kens_dir = add_directories(root_dir, directories)
        except OSError:
            pass
        else:
            msg = 'bad dir name should give OSError'
            raise Exception(msg)    
            
        #clean up!
        os.rmdir(root_dir)

    def test_check_list(self):

        check_list(['stage','xmomentum'])

        
    def test_add_directories(self):
        
        import tempfile
        root_dir = tempfile.mkdtemp('_test_util', 'test_util_')
        directories = ['ja','ne','ke']
        kens_dir = add_directories(root_dir, directories)
        assert kens_dir == root_dir + sep + 'ja' + sep + 'ne' + \
               sep + 'ke'
        assert access(root_dir,F_OK)

        add_directories(root_dir, directories)
        assert access(root_dir,F_OK)
        
        #clean up!
        os.rmdir(kens_dir)
        os.rmdir(root_dir + sep + 'ja' + sep + 'ne')
        os.rmdir(root_dir + sep + 'ja')
        os.rmdir(root_dir)

    def test_add_directories_bad(self):
        
        import tempfile
        root_dir = tempfile.mkdtemp('_test_util', 'test_util_')
        directories = ['/\/!@#@#$%^%&*((*:*:','ne','ke']
        
        try:
            kens_dir = add_directories(root_dir, directories)
        except OSError:
            pass
        else:
            msg = 'bad dir name should give OSError'
            raise Exception(msg)    
            
        #clean up!
        os.rmdir(root_dir)

    def test_check_list(self):

        check_list(['stage','xmomentum'])
        
    def test_remove_lone_verts_d(self):
        verts = [[0,0],[1,0],[0,1]]
        tris = [[0,1,2]]
        new_verts, new_tris = remove_lone_verts(verts, tris)
        assert new_verts == verts
        assert new_tris == tris
     

    def test_remove_lone_verts_e(self):
        verts = [[0,0],[1,0],[0,1],[99,99]]
        tris = [[0,1,2]]
        new_verts, new_tris = remove_lone_verts(verts, tris)
        assert new_verts == verts[0:3]
        assert new_tris == tris
        
    def test_remove_lone_verts_a(self):
        verts = [[99,99],[0,0],[1,0],[99,99],[0,1],[99,99]]
        tris = [[1,2,4]]
        new_verts, new_tris = remove_lone_verts(verts, tris)
        #print "new_verts", new_verts
        assert new_verts == [[0,0],[1,0],[0,1]]
        assert new_tris == [[0,1,2]]
     
    def test_remove_lone_verts_c(self):
        verts = [[0,0],[1,0],[99,99],[0,1]]
        tris = [[0,1,3]]
        new_verts, new_tris = remove_lone_verts(verts, tris)
        #print "new_verts", new_verts
        assert new_verts == [[0,0],[1,0],[0,1]]
        assert new_tris == [[0,1,2]]
        
    def test_remove_lone_verts_b(self):
        verts = [[0,0],[1,0],[0,1],[99,99],[99,99],[99,99]]
        tris = [[0,1,2]]
        new_verts, new_tris = remove_lone_verts(verts, tris)
        assert new_verts == verts[0:3]
        assert new_tris == tris
     

    def test_remove_lone_verts_e(self):
        verts = [[0,0],[1,0],[0,1],[99,99]]
        tris = [[0,1,2]]
        new_verts, new_tris = remove_lone_verts(verts, tris)
        assert new_verts == verts[0:3]
        assert new_tris == tris
        
    def test_get_min_max_values(self):
        
        list=[8,9,6,1,4]
        min1, max1 = get_min_max_values(list)
        
        assert min1==1 
        assert max1==9
        
    def test_get_min_max_values1(self):
        
        list=[-8,-9,-6,-1,-4]
        min1, max1 = get_min_max_values(list)
        
#        print 'min1,max1',min1,max1
        assert min1==-9 
        assert max1==-1

#    def test_get_min_max_values2(self):
#        '''
#        The min and max supplied are greater than the ones in the 
#        list and therefore are the ones returned
#        '''
#        list=[-8,-9,-6,-1,-4]
#        min1, max1 = get_min_max_values(list,-10,10)
#        
##        print 'min1,max1',min1,max1
#        assert min1==-10 
#        assert max1==10
        
    def bad_test_make_plots_from_csv_files(self):
        
        try: 
            import pylab
        except ImportError:
            #ANUGA don't need pylab to work so the system doesn't 
            #rely on pylab being installed 
            return
        
        if sys.platform == 'win32':  #Windows
        
            current_dir=getcwd()+sep+'abstract_2d_finite_volumes'
            temp_dir = tempfile.mkdtemp('','figures')
    #        print 'temp_dir',temp_dir
            fileName = temp_dir+sep+'time_series_3.csv'
            file = open(fileName,"w")
            file.write("Time,Stage,Speed,Momentum,Elevation\n\
1.0, 0, 0, 0, 10 \n\
2.0, 5, 2, 4, 10 \n\
3.0, 3, 3, 5, 10 \n")
            file.close()
    
            fileName1 = temp_dir+sep+'time_series_4.csv'
            file1 = open(fileName1,"w")
            file1.write("Time,Stage,Speed,Momentum,Elevation\n\
1.0, 0, 0, 0, 5 \n\
2.0, -5, -2, -4, 5 \n\
3.0, -4, -3, -5, 5 \n")
            file1.close()
    
            fileName2 = temp_dir+sep+'time_series_5.csv'
            file2 = open(fileName2,"w")
            file2.write("Time,Stage,Speed,Momentum,Elevation\n\
1.0, 0, 0, 0, 7 \n\
2.0, 4, -0.45, 57, 7 \n\
3.0, 6, -0.5, 56, 7 \n")
            file2.close()
            
            dir, name=os.path.split(fileName)
            make_plots_from_csv_file(directories_dic={dir:['gauge', 0, 0]},
                                output_dir=temp_dir,
                                base_name='time_series_',
                                plot_numbers=['3-5'],
                                quantities=['Speed','Stage','Momentum'],
                                assess_all_csv_files=True,
                                extra_plot_name='test')
            
    #        print 'stage+fileName[:-4]+test.png',dir+sep+'stage_'+name[:-4]+'_test.png'
            assert(access(dir+sep+'stage_'+name[:-4]+'_test.png',F_OK)==True)
            assert(access(dir+sep+'speed_'+name[:-4]+'_test.png',F_OK)==True)
            assert(access(dir+sep+'momentum_'+name[:-4]+'_test.png',F_OK)==True)
    
            dir1, name1=os.path.split(fileName1)
            assert(access(dir+sep+'stage_'+name1[:-4]+'_test.png',F_OK)==True)
            assert(access(dir+sep+'speed_'+name1[:-4]+'_test.png',F_OK)==True)
            assert(access(dir+sep+'momentum_'+name1[:-4]+'_test.png',F_OK)==True)
    
    
            dir2, name2=os.path.split(fileName2)
            assert(access(dir+sep+'stage_'+name2[:-4]+'_test.png',F_OK)==True)
            assert(access(dir+sep+'speed_'+name2[:-4]+'_test.png',F_OK)==True)
            assert(access(dir+sep+'momentum_'+name2[:-4]+'_test.png',F_OK)==True)
    
            del_dir(temp_dir)
        

    def test_sww2csv_gauges(self):

        def elevation_function(x, y):
            return -x
        
        """Most of this test was copied from test_interpolate test_interpole_sww2csv
        
        This is testing the gauge_sww2csv function, by creating a sww file and
        then exporting the gauges and checking the results.
        """
        
        # create mesh
        mesh_file = tempfile.mktemp(".tsh")    
        points = [[0.0,0.0],[6.0,0.0],[6.0,6.0],[0.0,6.0]]
        m = Mesh()
        m.add_vertices(points)
        m.auto_segment()
        m.generate_mesh(verbose=False)
        m.export_mesh_file(mesh_file)
        
        #Create shallow water domain
        domain = Domain(mesh_file)
        os.remove(mesh_file)
        
        domain.default_order=2
        domain.beta_h = 0

        #Set some field values
        domain.set_quantity('elevation', elevation_function)
        domain.set_quantity('friction', 0.03)
        domain.set_quantity('xmomentum', 3.0)
        domain.set_quantity('ymomentum', 4.0)

        ######################
        # Boundary conditions
        B = Transmissive_boundary(domain)
        domain.set_boundary( {'exterior': B})

        # This call mangles the stage values.
        domain.distribute_to_vertices_and_edges()
        domain.set_quantity('stage', 1.0)


        domain.set_name('datatest' + str(time.time()))
        domain.format = 'sww'
        domain.smooth = True
        domain.reduction = mean

        sww = get_dataobject(domain)
        sww.store_connectivity()
        sww.store_timestep(['stage', 'xmomentum', 'ymomentum','elevation'])
        domain.set_quantity('stage', 10.0) # This is automatically limited
        # so it will not be less than the elevation
        domain.time = 2.
        sww.store_timestep(['stage','elevation', 'xmomentum', 'ymomentum'])

        # test the function
        points = [[5.0,1.],[0.5,2.]]

        points_file = tempfile.mktemp(".csv")
#        points_file = 'test_point.csv'
        file_id = open(points_file,"w")
        file_id.write("name, easting, northing, elevation \n\
point1, 5.0, 1.0, 3.0\n\
point2, 0.5, 2.0, 9.0\n")
        file_id.close()

        
        sww2csv_gauges(sww.filename, 
                            points_file,
                            verbose=False,
                            use_cache=False)

#        point1_answers_array = [[0.0,1.0,-5.0,3.0,4.0], [2.0,10.0,-5.0,3.0,4.0]]
        point1_answers_array = [[0.0,1.0,6.0,-5.0,3.0,4.0], [2.0,10.0,15.0,-5.0,3.0,4.0]]
        point1_filename = 'gauge_point1.csv'
        point1_handle = file(point1_filename)
        point1_reader = reader(point1_handle)
        point1_reader.next()

        line=[]
        for i,row in enumerate(point1_reader):
            #print 'i',i,'row',row
            line.append([float(row[0]),float(row[1]),float(row[2]),float(row[3]),float(row[4]),float(row[5])])
            #print 'assert line',line[i],'point1',point1_answers_array[i]
            assert allclose(line[i], point1_answers_array[i])

        point2_answers_array = [[0.0,1.0,1.5,-0.5,3.0,4.0], [2.0,10.0,10.5,-0.5,3.0,4.0]]
        point2_filename = 'gauge_point2.csv' 
        point2_handle = file(point2_filename)
        point2_reader = reader(point2_handle)
        point2_reader.next()
                        
        line=[]
        for i,row in enumerate(point2_reader):
            #print 'i',i,'row',row
            line.append([float(row[0]),float(row[1]),float(row[2]),float(row[3]),float(row[4]),float(row[5])])
            #print 'assert line',line[i],'point1',point1_answers_array[i]
            assert allclose(line[i], point2_answers_array[i])
                         
        # clean up
        point1_handle.close()
        point2_handle.close()
        #print "sww.filename",sww.filename 
        os.remove(sww.filename)
        os.remove(points_file)
        os.remove(point1_filename)
        os.remove(point2_filename)



    def test_sww2csv_gauges1(self):
        from anuga.pmesh.mesh import Mesh
        from anuga.shallow_water import Domain, Transmissive_boundary
        from anuga.shallow_water.data_manager import get_dataobject
        from csv import reader,writer
        import time
        import string

        def elevation_function(x, y):
            return -x
        
        """Most of this test was copied from test_interpolate test_interpole_sww2csv
        
        This is testing the gauge_sww2csv function, by creating a sww file and
        then exporting the gauges and checking the results.
        
        This tests the ablity not to have elevation in the points file and 
        not store xmomentum and ymomentum
        """
        
        # create mesh
        mesh_file = tempfile.mktemp(".tsh")    
        points = [[0.0,0.0],[6.0,0.0],[6.0,6.0],[0.0,6.0]]
        m = Mesh()
        m.add_vertices(points)
        m.auto_segment()
        m.generate_mesh(verbose=False)
        m.export_mesh_file(mesh_file)
        
        #Create shallow water domain
        domain = Domain(mesh_file)
        os.remove(mesh_file)
        
        domain.default_order=2
        domain.beta_h = 0

        #Set some field values
        domain.set_quantity('elevation', elevation_function)
        domain.set_quantity('friction', 0.03)
        domain.set_quantity('xmomentum', 3.0)
        domain.set_quantity('ymomentum', 4.0)

        ######################
        # Boundary conditions
        B = Transmissive_boundary(domain)
        domain.set_boundary( {'exterior': B})

        # This call mangles the stage values.
        domain.distribute_to_vertices_and_edges()
        domain.set_quantity('stage', 1.0)


        domain.set_name('datatest' + str(time.time()))
        domain.format = 'sww'
        domain.smooth = True
        domain.reduction = mean

        sww = get_dataobject(domain)
        sww.store_connectivity()
        sww.store_timestep(['stage', 'xmomentum', 'ymomentum'])
        domain.set_quantity('stage', 10.0) # This is automatically limited
        # so it will not be less than the elevation
        domain.time = 2.
        sww.store_timestep(['stage', 'xmomentum', 'ymomentum'])

        # test the function
        points = [[5.0,1.],[0.5,2.]]

        points_file = tempfile.mktemp(".csv")
#        points_file = 'test_point.csv'
        file_id = open(points_file,"w")
        file_id.write("name,easting,northing \n\
point1, 5.0, 1.0\n\
point2, 0.5, 2.0\n")
        file_id.close()

        sww2csv_gauges(sww.filename, 
                            points_file,
                            quantities=['stage', 'elevation'],
                            use_cache=False,
                            verbose=False)

        point1_answers_array = [[0.0,1.0,-5.0], [2.0,10.0,-5.0]]
        point1_filename = 'gauge_point1.csv'
        point1_handle = file(point1_filename)
        point1_reader = reader(point1_handle)
        point1_reader.next()

        line=[]
        for i,row in enumerate(point1_reader):
#            print 'i',i,'row',row
            line.append([float(row[0]),float(row[1]),float(row[2])])
            #print 'line',line[i],'point1',point1_answers_array[i]
            assert allclose(line[i], point1_answers_array[i])

        point2_answers_array = [[0.0,1.0,-0.5], [2.0,10.0,-0.5]]
        point2_filename = 'gauge_point2.csv' 
        point2_handle = file(point2_filename)
        point2_reader = reader(point2_handle)
        point2_reader.next()
                        
        line=[]
        for i,row in enumerate(point2_reader):
#            print 'i',i,'row',row
            line.append([float(row[0]),float(row[1]),float(row[2])])
#            print 'line',line[i],'point1',point1_answers_array[i]
            assert allclose(line[i], point2_answers_array[i])
                         
        # clean up
        point1_handle.close()
        point2_handle.close()
        #print "sww.filename",sww.filename 
        os.remove(sww.filename)
        os.remove(points_file)
        os.remove(point1_filename)
        os.remove(point2_filename)

    def test_greens_law(self):

        from math import sqrt
        
        d1 = 80.0
        d2 = 20.0
        h1 = 1.0
        h2 = greens_law(d1,d2,h1)

        assert h2==sqrt(2.0)
       

#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Util,'test')
#    suite = unittest.makeSuite(Test_Util,'test_sww2csv')
#    runner = unittest.TextTestRunner(verbosity=2)
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)




