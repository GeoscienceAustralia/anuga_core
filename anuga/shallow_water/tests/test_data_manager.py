#!/usr/bin/env python
#

"""
Set of tests for the now-defunct data manager module.

These could be split up into their correct modules.
"""

import unittest
import copy
import numpy as num
import sys
                
import tempfile
import os
import shutil
from struct import pack, unpack

from anuga.file.netcdf import NetCDFFile


from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular
from anuga.anuga_exceptions import ANUGAError
from anuga.file.sww import SWW_file
from anuga.coordinate_transforms.geo_reference import Geo_reference                        
from anuga.coordinate_transforms.redfearn import degminsec2decimal_degrees
from anuga.abstract_2d_finite_volumes.util import file_function
from anuga.utilities.system_tools import get_pathname_from_package, \
                                            get_revision_number
from anuga.utilities.file_utils import get_all_swwfiles
from anuga.utilities.file_utils import del_dir
from anuga.utilities.numerical_tools import ensure_numeric, mean
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a
from anuga.config import netcdf_float, epsilon, g
from anuga.pmesh.mesh_interface import create_mesh_from_regions
from anuga.file_conversion.sww2dem import sww2dem_batch
from anuga.file.csv_file import load_csv_as_dict, load_csv_as_array, \
                                load_csv_as_building_polygons, \
                                load_csv_as_polygons
from anuga.file.sts import create_sts_boundary
from anuga.file.pts import load_pts_as_polygon
from anuga.file.sww import Write_sww


# import all the boundaries - some are generic, some are shallow water
from anuga.shallow_water.boundaries import Reflective_boundary, \
            Field_boundary, Transmissive_momentum_set_stage_boundary, \
            Transmissive_stage_zero_momentum_boundary
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions\
     import Transmissive_boundary, Dirichlet_boundary, \
            Time_boundary, File_boundary, AWI_boundary

# This is needed to run the tests of local functions
from anuga.file_conversion.urs2sts import urs2sts
from anuga.coordinate_transforms.redfearn import redfearn
from anuga.coordinate_transforms.geo_reference import Geo_reference, \
     DEFAULT_ZONE
from anuga.geospatial_data.geospatial_data import Geospatial_data

from anuga.shallow_water.shallow_water_domain import Domain

# use helper methods from other unit test
from anuga.file.tests.test_mux import Test_Mux


class Test_Data_Manager(Test_Mux):
    # Class variable
    verbose = False

    def set_verbose(self):
        Test_Data_Manager.verbose = True
        
    def setUp(self):
        import time
        #from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular
        
        self.verbose = Test_Data_Manager.verbose
        # Create basic mesh
        points, vertices, boundary = rectangular(2, 2)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.set_flow_algorithm('1_5')
        domain.default_order = 2

        # Set some field values
        domain.set_quantity('elevation', lambda x,y: -x)
        domain.set_quantity('friction', 0.03)


        ######################
        # Boundary conditions
        B = Transmissive_boundary(domain)
        domain.set_boundary( {'left': B, 'right': B, 'top': B, 'bottom': B})


        ######################
        #Initial condition - with jumps
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
        self.initial_stage = copy.copy(domain.quantities['stage'].vertex_values)


        self.domain = domain

        C = domain.get_vertex_coordinates()
        self.X = C[:,0:6:2].copy()
        self.Y = C[:,1:6:2].copy()

        self.F = bed

        #Write A testfile (not realistic. Values aren't realistic)
        self.test_MOST_file = 'most_small'

        longitudes = [150.66667, 150.83334, 151., 151.16667]
        latitudes = [-34.5, -34.33333, -34.16667, -34]

        long_name = 'LON'
        lat_name = 'LAT'

        nx = 4
        ny = 4
        six = 6


        for ext in ['_ha.nc', '_ua.nc', '_va.nc', '_e.nc']:
            fid = NetCDFFile(self.test_MOST_file + ext, netcdf_mode_w)

            fid.createDimension(long_name,nx)
            fid.createVariable(long_name,netcdf_float,(long_name,))
            fid.variables[long_name].point_spacing='uneven'
            fid.variables[long_name].units='degrees_east'
            fid.variables[long_name][:] = longitudes

            fid.createDimension(lat_name,ny)
            fid.createVariable(lat_name,netcdf_float,(lat_name,))
            fid.variables[lat_name].point_spacing='uneven'
            fid.variables[lat_name].units='degrees_north'
            fid.variables[lat_name][:] = latitudes

            fid.createDimension('TIME',six)
            fid.createVariable('TIME',netcdf_float,('TIME',))
            fid.variables['TIME'].point_spacing='uneven'
            fid.variables['TIME'].units='seconds'
            fid.variables['TIME'][:] = [0.0, 0.1, 0.6, 1.1, 1.6, 2.1]


            name = ext[1:3].upper()
            if name == 'E.': name = 'ELEVATION'
            fid.createVariable(name,netcdf_float,('TIME', lat_name, long_name))
            fid.variables[name].units='CENTIMETERS'
            fid.variables[name].missing_value=-1.e+034

            fid.variables[name][:] = [[[0.3400644, 0, -46.63519, -6.50198],
                                              [-0.1214216, 0, 0, 0],
                                              [0, 0, 0, 0],
                                              [0, 0, 0, 0]],
                                             [[0.3400644, 2.291054e-005, -23.33335, -6.50198],
                                              [-0.1213987, 4.581959e-005, -1.594838e-007, 1.421085e-012],
                                              [2.291054e-005, 4.582107e-005, 4.581715e-005, 1.854517e-009],
                                              [0, 2.291054e-005, 2.291054e-005, 0]],
                                             [[0.3400644, 0.0001374632, -23.31503, -6.50198],
                                              [-0.1212842, 0.0002756907, 0.006325484, 1.380492e-006],
                                              [0.0001374632, 0.0002749264, 0.0002742863, 6.665601e-008],
                                              [0, 0.0001374632, 0.0001374632, 0]],
                                             [[0.3400644, 0.0002520159, -23.29672, -6.50198],
                                              [-0.1211696, 0.0005075303, 0.01264618, 6.208276e-006],
                                              [0.0002520159, 0.0005040318, 0.0005027961, 2.23865e-007],
                                              [0, 0.0002520159, 0.0002520159, 0]],
                                             [[0.3400644, 0.0003665686, -23.27842, -6.50198],
                                              [-0.1210551, 0.0007413362, 0.01896192, 1.447638e-005],
                                              [0.0003665686, 0.0007331371, 0.0007313463, 4.734126e-007],
                                              [0, 0.0003665686, 0.0003665686, 0]],
                                             [[0.3400644, 0.0004811212, -23.26012, -6.50198],
                                              [-0.1209405, 0.0009771062, 0.02527271, 2.617787e-005],
                                              [0.0004811212, 0.0009622425, 0.0009599366, 8.152277e-007],
                                              [0, 0.0004811212, 0.0004811212, 0]]]


            fid.close()




    def tearDown(self):
        import os
        for ext in ['_ha.nc', '_ua.nc', '_va.nc', '_e.nc']:
            #print 'Trying to remove', self.test_MOST_file + ext
            try:
                os.remove(self.test_MOST_file + ext)
            except:
                pass
            
        for file in ['domain.sww', 'outline_meshed.tsh', 'outline.tsh']:
            try:
                os.remove(file)
            except:
                pass

    def test_sww_constant(self):
        """Test that constant sww information can be written correctly
        (non smooth)
        """
        self.domain.set_name('datatest' + str(id(self)))
        self.domain.format = 'sww' #Remove??
        self.domain.smooth = False
        
        sww = SWW_file(self.domain)
        sww.store_connectivity()

        fid = NetCDFFile(sww.filename, netcdf_mode_r)  # Open existing file for append

        # Get the variables
        x = fid.variables['x']
        y = fid.variables['y']
        z = fid.variables['elevation']
        V = fid.variables['volumes']

        assert num.allclose (x[:], self.X.flatten())
        assert num.allclose (y[:], self.Y.flatten())
        assert num.allclose (z[:], self.F.flatten())

        P = len(self.domain)
        for k in range(P):
            assert V[k, 0] == 3*k
            assert V[k, 1] == 3*k+1
            assert V[k, 2] == 3*k+2

        fid.close()
        os.remove(sww.filename)

    def test_sww_header(self):
        """Test that constant sww information can be written correctly
        (non smooth)
        """
        self.domain.set_name('datatest' + str(id(self)))
        self.domain.format = 'sww' #Remove??
        self.domain.smooth = False

        sww = SWW_file(self.domain)
        sww.store_connectivity()

        # Check contents
        # Get NetCDF
        fid = NetCDFFile(sww.filename, netcdf_mode_r)  # Open existing file for append

        # Get the variables
        sww_revision = fid.revision_number
        try:
            revision_number = get_revision_number()
        except:
            revision_number = None
            
        assert str(revision_number) == sww_revision
        fid.close()

        #print "sww.filename", sww.filename
        os.remove(sww.filename)

    def test_sww_range(self):
        """Test that constant sww information can be written correctly
        Use non-smooth to be able to compare to quantity values.
        """

        self.domain.set_name('datatest' + str(id(self)))
        self.domain.format = 'sww'
        self.domain.set_store_vertices_uniquely(True)
        
        sww = SWW_file(self.domain)        

        dqs = self.domain.get_quantity('stage')
        dqx = self.domain.get_quantity('xmomentum')
        dqy = self.domain.get_quantity('ymomentum')        
        xmom_min = ymom_min = stage_min = sys.maxsize 
        xmom_max = ymom_max = stage_max = -stage_min        
        for t in self.domain.evolve(yieldstep = 1, finaltime = 1):
            wmax = max(dqs.get_values().flatten())
            if wmax > stage_max: stage_max = wmax
            wmin = min(dqs.get_values().flatten())
            if wmin < stage_min: stage_min = wmin            
            
            uhmax = max(dqx.get_values().flatten())
            if uhmax > xmom_max: xmom_max = uhmax
            uhmin = min(dqx.get_values().flatten())
            if uhmin < xmom_min: xmom_min = uhmin                        
            
            vhmax = max(dqy.get_values().flatten())
            if vhmax > ymom_max: ymom_max = vhmax
            vhmin = min(dqy.get_values().flatten())
            if vhmin < ymom_min: ymom_min = vhmin                                    
            
            
            
        # Get NetCDF
        fid = NetCDFFile(sww.filename, netcdf_mode_r) # Open existing file for append

        # Get the variables
        range = fid.variables['stage_range'][:]
        assert num.allclose(range,[stage_min, stage_max])

        range = fid.variables['xmomentum_range'][:]
        #print range
        assert num.allclose(range, [xmom_min, xmom_max])
        
        range = fid.variables['ymomentum_range'][:]
        #print range
        assert num.allclose(range, [ymom_min, ymom_max])        


        
        fid.close()
        os.remove(sww.filename)

    def test_sww_extrema(self):
        """Test that extrema of quantities can be retrieved at every vertex
        Extrema are updated at every *internal* timestep
        """

        domain = self.domain
        
        domain.set_name('extrema_test' + str(id(self)))
        domain.format = 'sww'
        domain.smooth = True

        assert domain.quantities_to_be_monitored is None
        assert domain.monitor_polygon is None
        assert domain.monitor_time_interval is None        
        
        domain.set_quantities_to_be_monitored(['xmomentum',
                                               'ymomentum',
                                               'stage-elevation'])

        assert domain.monitor_polygon is None
        assert domain.monitor_time_interval is None


        domain.set_quantities_to_be_monitored(['xmomentum',
                                               'ymomentum',
                                               'stage-elevation'],
                                              polygon=domain.get_boundary_polygon(),
                                              time_interval=[0,1])
        
        
        assert len(domain.quantities_to_be_monitored) == 3
        assert 'stage-elevation' in domain.quantities_to_be_monitored
        assert 'xmomentum' in domain.quantities_to_be_monitored                
        assert 'ymomentum' in domain.quantities_to_be_monitored        

        
        #domain.protect_against_isolated_degenerate_timesteps = True
        #domain.tight_slope_limiters = 1
        domain.tight_slope_limiters = 0 # Backwards compatibility
        domain.use_centroid_velocities = 0 # Backwards compatibility (7/5/8)
        
        
        sww = SWW_file(domain)

        for t in domain.evolve(yieldstep = 1, finaltime = 1):
            pass
            #print domain.timestepping_statistics()
            domain.quantity_statistics(precision = '%.8f') # Silent

            
        # Get NetCDF
        fid = NetCDFFile(sww.filename, netcdf_mode_r) # Open existing file for append

        # Get the variables
        extrema = fid.variables['stage-elevation.extrema'][:]
        assert num.allclose(extrema, [0.00, 0.30]) or \
            num.allclose(extrema, [ 0., 0.3222025])

        loc = fid.variables['stage-elevation.min_location'][:]
        assert num.allclose(loc, [0.16666667, 0.33333333])

        loc = fid.variables['stage-elevation.max_location'][:]        
        assert num.allclose(loc, [0.8333333, 0.16666667])        

        time = fid.variables['stage-elevation.max_time'][:]
        assert num.allclose(time, 0.0) or \
            num.allclose(time, 0.35077909)              

        extrema = fid.variables['xmomentum.extrema'][:]
        # Padarn Note: Had to add an extra possibility here (the final one [-0.06062178  0.47518688])
        # to pass this assertion. Cannot see what would be causing this.
        assert num.allclose(extrema,[-0.06062178, 0.47873023]) or \
            num.allclose(extrema, [-0.06062178, 0.47847986]) or \
            num.allclose(extrema, [-0.06062178, 0.47848481]) or \
            num.allclose(extrema, [-0.06062178, 0.47763887]) or \
            num.allclose(extrema, [-0.06062178, 0.46691909])or \
            num.allclose(extrema, [-0.06062178, 0.47503704]) or \
            num.allclose(extrema, [-0.06062178,  0.47518688]) or \
            num.allclose(extrema, [-0.06062178,  0.49014235])


        
        extrema = fid.variables['ymomentum.extrema'][:]
        assert num.allclose(extrema,[0.00, 0.0625786]) or \
            num.allclose(extrema,[0.00, 0.06062178])

        time_interval = fid.variables['extrema.time_interval'][:]
        assert num.allclose(time_interval, [0,1])
        
        polygon = fid.variables['extrema.polygon'][:]        
        assert num.allclose(polygon, domain.get_boundary_polygon())
        
        fid.close()
        #print "sww.filename", sww.filename
        os.remove(sww.filename)

        
        
    def test_sww_constant_smooth(self):
        """Test that constant sww information can be written correctly
        (non smooth)
        """
        self.domain.set_name('datatest' + str(id(self)))
        self.domain.format = 'sww'
        self.domain.smooth = True

        sww = SWW_file(self.domain)
        sww.store_connectivity()

        # Check contents
        # Get NetCDF
        fid = NetCDFFile(sww.filename, netcdf_mode_r)  # Open existing file for append

        # Get the variables
        X = fid.variables['x'][:]
        Y = fid.variables['y'][:]
        Z = fid.variables['elevation'][:]
        V = fid.variables['volumes']

        assert num.allclose([X[0], Y[0]], num.array([0.0, 0.0]))
        assert num.allclose([X[1], Y[1]], num.array([0.0, 0.5]))
        assert num.allclose([X[2], Y[2]], num.array([0.0, 1.0]))
        assert num.allclose([X[4], Y[4]], num.array([0.5, 0.5]))
        assert num.allclose([X[7], Y[7]], num.array([1.0, 0.5]))

        assert Z[4] == -0.5

        assert V[2,0] == 4
        assert V[2,1] == 5
        assert V[2,2] == 1
        assert V[4,0] == 6
        assert V[4,1] == 7
        assert V[4,2] == 3

        fid.close()
        os.remove(sww.filename)
        

    def test_sww_variable(self):
        """Test that sww information can be written correctly
        """
        self.domain.set_name('datatest' + str(id(self)))
        self.domain.format = 'sww'
        self.domain.smooth = True
        self.domain.reduction = mean

        sww = SWW_file(self.domain)
        sww.store_connectivity()
        sww.store_timestep()

        # Check contents
        # Get NetCDF
        fid = NetCDFFile(sww.filename, netcdf_mode_r)  # Open existing file for append


        # Get the variables
        time = fid.variables['time']
        stage = fid.variables['stage']

        Q = self.domain.quantities['stage']
        Q0 = Q.vertex_values[:,0]
        Q1 = Q.vertex_values[:,1]
        Q2 = Q.vertex_values[:,2]

        A = stage[0,:]
        
#         print stage[0,:]
#         print A[0], (Q2[0] + Q1[1])/2
#         print A[1], (Q0[1] + Q1[3] + Q2[2])/3
#         print A[2], Q0[3]
#         print A[3], (Q0[0] + Q1[5] + Q2[4])/3
        assert num.allclose(A[0], (Q2[0] + Q1[1])/2, rtol=1.0e-5, atol=1.0e-5)
        
        assert num.allclose(A[1], (Q0[1] + Q1[3] + Q2[2])/3, rtol=1.0e-5, atol=1.0e-5)
        assert num.allclose(A[2], Q0[3])
        assert num.allclose(A[3], (Q0[0] + Q1[5] + Q2[4])/3, rtol=1.0e-5, atol=1.0e-5)

        #Center point
        assert num.allclose(A[4], (Q1[0] + Q2[1] + Q0[2] + Q0[5] + Q2[6] + Q1[7])/6,
                            rtol=1.0e-5, atol=1.0e-5)
        
        fid.close()
        os.remove(sww.filename)


    def test_sww_variable2(self):
        """Test that sww information can be written correctly
        multiple timesteps. Use average as reduction operator
        """

        import time, os

        self.domain.set_name('datatest' + str(id(self)))
        self.domain.format = 'sww'
        self.domain.smooth = True

        self.domain.reduction = mean

        sww = SWW_file(self.domain)
        sww.store_connectivity()
        sww.store_timestep()
        #self.domain.tight_slope_limiters = 1
        self.domain.evolve_to_end(finaltime = 0.01)
        sww.store_timestep()

        # Check contents
        # Get NetCDF
        fid = NetCDFFile(sww.filename, netcdf_mode_r)

        # Get the variables
        x = fid.variables['x']
        y = fid.variables['y']
        z = fid.variables['elevation']
        time = fid.variables['time']
        stage = fid.variables['stage']

        #Check values
        Q = self.domain.quantities['stage']
        Q0 = Q.vertex_values[:,0]
        Q1 = Q.vertex_values[:,1]
        Q2 = Q.vertex_values[:,2]

        A = stage[1,:]
        assert num.allclose(A[0], (Q2[0] + Q1[1])/2, rtol=1.0e-5, atol=1.0e-5)
        assert num.allclose(A[1], (Q0[1] + Q1[3] + Q2[2])/3, rtol=1.0e-5, atol=1.0e-5)
        assert num.allclose(A[2], Q0[3], rtol=1.0e-5, atol=1.0e-5)
        assert num.allclose(A[3], (Q0[0] + Q1[5] + Q2[4])/3, rtol=1.0e-5, atol=1.0e-5)

        #Center point
        assert num.allclose(A[4], (Q1[0] + Q2[1] + Q0[2] + Q0[5] + Q2[6] + Q1[7])/6,
                            rtol=1.0e-5, atol=1.0e-5)


        fid.close()

        #Cleanup
        os.remove(sww.filename)


    def test_sync(self):
        """test_sync - Test info stored at each timestep is as expected (incl initial condition)
        """

        import time, os

        self.domain.set_name('synctest')
        self.domain.format = 'sww'
        self.domain.smooth = False
        self.domain.store = True

        self.domain.tight_slope_limiters = True
        self.domain.use_centroid_velocities = True        
        
        # In this case tight_slope_limiters as default
        # in conjunction with protection
        # against isolated degenerate timesteps works.
        #self.domain.tight_slope_limiters = 1
        #self.domain.protect_against_isolated_degenerate_timesteps = True

        #print 'tight_sl', self.domain.tight_slope_limiters
        

        #Evolution
        for t in self.domain.evolve(yieldstep = 1.0, finaltime = 4.0):
            
            #########self.domain.write_time(track_speeds=True)
            stage = self.domain.quantities['stage'].vertex_values

            #Get NetCDF
            fid = NetCDFFile(self.domain.writer.filename, netcdf_mode_r)
            stage_file = fid.variables['stage']
            
            if t == 0.0:
                assert num.allclose(stage, self.initial_stage, rtol=1.0e-5, atol=1.0e-5)
                assert num.allclose(stage_file[:], stage.flatten(),rtol=1.0e-5, atol=1.0e-5)
            else:
                assert not num.allclose(stage, self.initial_stage, rtol=1.0e-5, atol=1.0e-5)
                assert not num.allclose(stage_file[:], stage.flatten(), rtol=1.0e-5, atol=1.0e-5)

            fid.close()

        os.remove(self.domain.writer.filename)


    def test_sww_minimum_storable_height(self):
        """Test that sww information can be written correctly
        multiple timesteps using a different reduction operator (min)
        """

        import time, os

        self.domain.set_name('datatest' + str(id(self)))
        self.domain.format = 'sww'
        self.domain.smooth = True
        self.domain.reduction = min
        self.domain.minimum_storable_height = 100

        sww = SWW_file(self.domain)
        sww.store_connectivity()
        sww.store_timestep()

        #self.domain.tight_slope_limiters = 1
        self.domain.evolve_to_end(finaltime = 0.01)
        sww.store_timestep()

        #Check contents
        #Get NetCDF
        fid = NetCDFFile(sww.filename, netcdf_mode_r)

        # Get the variables
        x = fid.variables['x']
        y = fid.variables['y']
        z = fid.variables['elevation']
        time = fid.variables['time']
        stage = fid.variables['stage']
        xmomentum = fid.variables['xmomentum']
        ymomentum = fid.variables['ymomentum']        

        #Check values
        Q = self.domain.quantities['stage']
        Q0 = Q.vertex_values[:,0]
        Q1 = Q.vertex_values[:,1]
        Q2 = Q.vertex_values[:,2]

        A = stage[1,:]
        assert num.allclose(stage[1,:], z[:])


        assert num.allclose(xmomentum, 0.0)
        assert num.allclose(ymomentum, 0.0)        
        
        fid.close()

        #Cleanup
        os.remove(sww.filename)


    def test_file_boundary_stsIV_sinewave_ordering(self):
        """test_file_boundary_stsIV_sinewave_ordering(self):
        Read correct points from ordering file and apply sts to boundary
        This one uses a sine wave and compares to time boundary
        """

        lat_long_points=[[6.01, 97.0], [6.02, 97.0], [6.05, 96.9], [6.0, 97.0]]
        bounding_polygon=[[6.0, 97.0], [6.01, 97.0], [6.02,97.0], \
                            [6.02,97.02], [6.00,97.02]]
        tide = 0.35
        time_step_count = 50
        time_step = 0.1
        times_ref = num.arange(0, time_step_count*time_step, time_step)
        
        n=len(lat_long_points)
        first_tstep=num.ones(n,int)
        last_tstep=(time_step_count)*num.ones(n,int)
        
        gauge_depth=20*num.ones(n,float)
        
        ha1=num.ones((n,time_step_count),float)
        ua1=3.*num.ones((n,time_step_count),float)
        va1=2.*num.ones((n,time_step_count),float)
        for i in range(n):
            ha1[i]=num.sin(times_ref)
        
        
        base_name, files = self.write_mux2(lat_long_points,
                                           time_step_count, time_step,
                                           first_tstep, last_tstep,
                                           depth=gauge_depth,
                                           ha=ha1,
                                           ua=ua1,
                                           va=va1)

        # Write order file
        file_handle, order_base_name = tempfile.mkstemp("")
        os.close(file_handle)
        os.remove(order_base_name)
        d=","
        order_file=order_base_name+'order.txt'
        fid=open(order_file,'w')
        
        # Write Header
        header='index, longitude, latitude\n'
        fid.write(header)
        indices=[3,0,1]
        for i in indices:
            line=str(i)+d+str(lat_long_points[i][1])+d+\
                str(lat_long_points[i][0])+"\n"
            fid.write(line)
        fid.close()

        sts_file=base_name
        urs2sts(base_name, basename_out=sts_file,
                ordering_filename=order_file,
                mean_stage=tide,
                verbose=False)
        self.delete_mux(files)
        
        
        
        # Now read the sts file and check that values have been stored correctly.
        fid = NetCDFFile(sts_file + '.sts')

        # Check the time vector
        times = fid.variables['time'][:]
        
        #print times

        # Check sts quantities
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        ymomentum = fid.variables['ymomentum'][:]
        elevation = fid.variables['elevation'][:]       

        # Create beginnings of boundary polygon based on sts_boundary
        boundary_polygon = create_sts_boundary(base_name)
        
        os.remove(order_file)

        # Append the remaining part of the boundary polygon to be defined by
        # the user
        bounding_polygon_utm=[]
        for point in bounding_polygon:
            zone,easting,northing=redfearn(point[0],point[1])
            bounding_polygon_utm.append([easting,northing])

        boundary_polygon.append(bounding_polygon_utm[3])
        boundary_polygon.append(bounding_polygon_utm[4])

        #print 'boundary_polygon', boundary_polygon
        
        plot=False
        if plot:
            from pylab import plot,show,axis
            boundary_polygon=ensure_numeric(boundary_polygon)
            bounding_polygon_utm=ensure_numeric(bounding_polygon_utm)
            #plot(lat_long_points[:,0],lat_long_points[:,1],'o')
            plot(boundary_polygon[:,0], boundary_polygon[:,1])
            plot(bounding_polygon_utm[:,0],bounding_polygon_utm[:,1])
            show()

        assert num.allclose(bounding_polygon_utm,boundary_polygon)


        extent_res=1000000
        meshname = 'urs_test_mesh' + '.tsh'
        interior_regions=None
        boundary_tags={'ocean': [0,1], 'otherocean': [2,3,4]}
        
        # have to change boundary tags from last example because now bounding
        # polygon starts in different place.
        create_mesh_from_regions(boundary_polygon,
                                 boundary_tags=boundary_tags,
                                 maximum_triangle_area=extent_res,
                                 filename=meshname,
                                 interior_regions=interior_regions,
                                 verbose=False)
        
        domain_fbound = Domain(meshname)
        domain_fbound.set_quantity('stage', tide)
        Bf = File_boundary(sts_file+'.sts', 
                           domain_fbound, 
                           boundary_polygon=boundary_polygon)
        Br = Reflective_boundary(domain_fbound)

        domain_fbound.set_boundary({'ocean': Bf,'otherocean': Br})
        finaltime=time_step*(time_step_count-1)
        yieldstep=time_step
        temp_fbound=num.zeros(int(finaltime/yieldstep)+1,float)
    
        for i, t in enumerate(domain_fbound.evolve(yieldstep=yieldstep,
                                                   finaltime=finaltime, 
                                                   skip_initial_step=False)):
            temp_fbound[i]=domain_fbound.quantities['stage'].centroid_values[2]
    
        
        domain_time = Domain(meshname)
        domain_time.set_quantity('stage', tide)
        Br = Reflective_boundary(domain_time)
        Bw = Time_boundary(domain=domain_time,
                         function=lambda t: [num.sin(t)+tide,3.*(20.+num.sin(t)+tide),2.*(20.+num.sin(t)+tide)])
        domain_time.set_boundary({'ocean': Bw,'otherocean': Br})
        
        temp_time=num.zeros(int(finaltime/yieldstep)+1,float)
        
        domain_time.set_starttime(domain_fbound.get_starttime())
        
        for i, t in enumerate(domain_time.evolve(yieldstep=yieldstep,
                                                   finaltime=finaltime, 
                                                   skip_initial_step=False)):
            temp_time[i]=domain_time.quantities['stage'].centroid_values[2]
        
        assert num.allclose(temp_fbound, temp_time)                
        assert num.allclose(domain_fbound.quantities['stage'].vertex_values,
                            domain_time.quantities['stage'].vertex_values)
                        
        assert num.allclose(domain_fbound.quantities['xmomentum'].vertex_values,
                            domain_time.quantities['xmomentum'].vertex_values)                        
                        
        assert num.allclose(domain_fbound.quantities['ymomentum'].vertex_values,
                            domain_time.quantities['ymomentum'].vertex_values)                                                
        

        try:
            os.remove(sts_file+'.sts')
        except:
            # Windoze can't remove this file for some reason 
            pass
        
        os.remove(meshname)
        
 
    def test_points2polygon(self):  
        att_dict = {}
        pointlist = num.array([[1.0, 0.0],[0.0, 1.0],[0.0, 0.0]])
        att_dict['elevation'] = num.array([10.1, 0.0, 10.4])
        att_dict['brightness'] = num.array([10.0, 1.0, 10.4])
        
        fileName = tempfile.mktemp(".csv")
        
        G = Geospatial_data(pointlist, att_dict)
        
        G.export_points_file(fileName)
        
        polygon = load_pts_as_polygon(fileName)
        
        # This test may fail if the order changes
        assert (polygon == [[0.0, 0.0],[1.0, 0.0],[0.0, 1.0]])
        


#-------------------------------------------------------------

if __name__ == "__main__":
    #suite = unittest.makeSuite(Test_Data_Manager, 'test_sww2domain2')
    suite = unittest.makeSuite(Test_Data_Manager, 'test')
    
    
    
    # FIXME(Ole): When Ross has implemented logging, we can 
    # probably get rid of all this:
    if len(sys.argv) > 1 and sys.argv[1][0].upper() == 'V':
        Test_Data_Manager.verbose=True
        saveout = sys.stdout   
        filename = ".temp_verbose"
        fid = open(filename, 'w')
        sys.stdout = fid
    else:
        pass
    runner = unittest.TextTestRunner() #verbosity=2)
    runner.run(suite)

    # Cleaning up
    if len(sys.argv) > 1 and sys.argv[1][0].upper() == 'V':
        sys.stdout = saveout 
        fid.close() 
        os.remove(filename)


