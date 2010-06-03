#!/usr/bin/env python
#

# This file was reverted from changeset:5484 to changeset:5470 on 10th July 
# by Ole.

import unittest
import copy
import numpy as num
                
import tempfile
import os
import shutil
from struct import pack, unpack

from Scientific.IO.NetCDF import NetCDFFile

from anuga.anuga_exceptions import ANUGAError
from anuga.shallow_water.data_manager import *
from anuga.shallow_water.sww_file import SWW_file
from anuga.coordinate_transforms.geo_reference import Geo_reference                        
from anuga.coordinate_transforms.redfearn import degminsec2decimal_degrees
from anuga.abstract_2d_finite_volumes.util import file_function
from anuga.utilities.system_tools import get_pathname_from_package
from anuga.utilities.file_utils import del_dir
from anuga.utilities.numerical_tools import ensure_numeric, mean
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a
from anuga.config import netcdf_float, epsilon, g

from anuga.file.csv_file import load_csv_as_dict, load_csv_as_array
from anuga.file.sts import create_sts_boundary


# import all the boundaries - some are generic, some are shallow water
from boundaries import Reflective_boundary, \
            Field_boundary, Transmissive_momentum_set_stage_boundary, \
            Transmissive_stage_zero_momentum_boundary
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions\
     import Transmissive_boundary, Dirichlet_boundary, \
            Time_boundary, File_boundary, AWI_boundary

# This is needed to run the tests of local functions
import data_manager 
from anuga.file_conversion.urs2sts import urs2sts
from anuga.coordinate_transforms.redfearn import redfearn
from anuga.coordinate_transforms.geo_reference import Geo_reference, \
     DEFAULT_ZONE
from anuga.geospatial_data.geospatial_data import Geospatial_data

# use helper methods from other unit test
from anuga.file.test_mux import Test_Mux


class Test_Data_Manager(Test_Mux):
    # Class variable
    verbose = False

    def set_verbose(self):
        Test_Data_Manager.verbose = True
        
    def setUp(self):
        import time
        from mesh_factory import rectangular
        
        self.verbose = Test_Data_Manager.verbose
        # Create basic mesh
        points, vertices, boundary = rectangular(2, 2)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
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
        stage = num.zeros(bed.shape, num.float)

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
            fid.variables[long_name].assignValue(longitudes)

            fid.createDimension(lat_name,ny)
            fid.createVariable(lat_name,netcdf_float,(lat_name,))
            fid.variables[lat_name].point_spacing='uneven'
            fid.variables[lat_name].units='degrees_north'
            fid.variables[lat_name].assignValue(latitudes)

            fid.createDimension('TIME',six)
            fid.createVariable('TIME',netcdf_float,('TIME',))
            fid.variables['TIME'].point_spacing='uneven'
            fid.variables['TIME'].units='seconds'
            fid.variables['TIME'].assignValue([0.0, 0.1, 0.6, 1.1, 1.6, 2.1])


            name = ext[1:3].upper()
            if name == 'E.': name = 'ELEVATION'
            fid.createVariable(name,netcdf_float,('TIME', lat_name, long_name))
            fid.variables[name].units='CENTIMETERS'
            fid.variables[name].missing_value=-1.e+034

            fid.variables[name].assignValue([[[0.3400644, 0, -46.63519, -6.50198],
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
                                              [0, 0.0004811212, 0.0004811212, 0]]])


            fid.close()




    def tearDown(self):
        import os
        for ext in ['_ha.nc', '_ua.nc', '_va.nc', '_e.nc']:
            #print 'Trying to remove', self.test_MOST_file + ext
            os.remove(self.test_MOST_file + ext)

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
        xmom_min = ymom_min = stage_min = sys.maxint 
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
        assert domain.quantities_to_be_monitored.has_key('stage-elevation')
        assert domain.quantities_to_be_monitored.has_key('xmomentum')                
        assert domain.quantities_to_be_monitored.has_key('ymomentum')        

        
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
        assert num.allclose(extrema, [0.00, 0.30])

        loc = fid.variables['stage-elevation.min_location'][:]
        assert num.allclose(loc, [0.16666667, 0.33333333])

        loc = fid.variables['stage-elevation.max_location'][:]        
        assert num.allclose(loc, [0.8333333, 0.16666667])        

        time = fid.variables['stage-elevation.max_time'][:]
        assert num.allclose(time, 0.0)                

        extrema = fid.variables['xmomentum.extrema'][:]
        assert num.allclose(extrema,[-0.06062178, 0.47873023]) or \
            num.allclose(extrema, [-0.06062178, 0.47847986]) or \
            num.allclose(extrema, [-0.06062178, 0.47848481]) or \
            num.allclose(extrema, [-0.06062178, 0.47763887]) # 18/09/09
        
        extrema = fid.variables['ymomentum.extrema'][:]
        assert num.allclose(extrema,[0.00, 0.0625786]) or num.allclose(extrema,[0.00, 0.06062178])

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
        #print A[0], (Q2[0,0] + Q1[1,0])/2
        assert num.allclose(A[0], (Q2[0] + Q1[1])/2)
        assert num.allclose(A[1], (Q0[1] + Q1[3] + Q2[2])/3)
        assert num.allclose(A[2], Q0[3])
        assert num.allclose(A[3], (Q0[0] + Q1[5] + Q2[4])/3)

        #Center point
        assert num.allclose(A[4], (Q1[0] + Q2[1] + Q0[2] +\
                                   Q0[5] + Q2[6] + Q1[7])/6)
        
        fid.close()
        os.remove(sww.filename)


    def test_sww_variable2(self):
        """Test that sww information can be written correctly
        multiple timesteps. Use average as reduction operator
        """

        import time, os
        from Scientific.IO.NetCDF import NetCDFFile

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
        fid = NetCDFFile(sww.filename, netcdf_mode_r)  # Open existing file for append

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
        assert num.allclose(A[0], (Q2[0] + Q1[1])/2)
        assert num.allclose(A[1], (Q0[1] + Q1[3] + Q2[2])/3)
        assert num.allclose(A[2], Q0[3])
        assert num.allclose(A[3], (Q0[0] + Q1[5] + Q2[4])/3)

        #Center point
        assert num.allclose(A[4], (Q1[0] + Q2[1] + Q0[2] +\
                                   Q0[5] + Q2[6] + Q1[7])/6)


        fid.close()

        #Cleanup
        os.remove(sww.filename)

    def no_test_sww_variable3(self):
        """Test that sww information can be written correctly
        multiple timesteps using a different reduction operator (min)
        """

        # Different reduction in sww files has been made obsolete.
        
        import time, os
        from Scientific.IO.NetCDF import NetCDFFile

        self.domain.set_name('datatest' + str(id(self)))
        self.domain.format = 'sww'
        self.domain.smooth = True
        self.domain.reduction = min

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

        #Check values
        Q = self.domain.quantities['stage']
        Q0 = Q.vertex_values[:,0]
        Q1 = Q.vertex_values[:,1]
        Q2 = Q.vertex_values[:,2]

        A = stage[1,:]
        assert num.allclose(A[0], min(Q2[0], Q1[1]))
        assert num.allclose(A[1], min(Q0[1], Q1[3], Q2[2]))
        assert num.allclose(A[2], Q0[3])
        assert num.allclose(A[3], min(Q0[0], Q1[5], Q2[4]))

        #Center point
        assert num.allclose(A[4], min(Q1[0], Q2[1], Q0[2],
                                      Q0[5], Q2[6], Q1[7]))


        fid.close()

        #Cleanup
        os.remove(sww.filename)


    def test_sync(self):
        """test_sync - Test info stored at each timestep is as expected (incl initial condition)
        """

        import time, os, config
        from Scientific.IO.NetCDF import NetCDFFile

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
                assert num.allclose(stage, self.initial_stage)
                assert num.allclose(stage_file[:], stage.flatten())
            else:
                assert not num.allclose(stage, self.initial_stage)
                assert not num.allclose(stage_file[:], stage.flatten())

            fid.close()

        os.remove(self.domain.writer.filename)


    def test_sww_minimum_storable_height(self):
        """Test that sww information can be written correctly
        multiple timesteps using a different reduction operator (min)
        """

        import time, os
        from Scientific.IO.NetCDF import NetCDFFile

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


    def Not_a_test_sww_DSG(self):
        """Not a test, rather a look at the sww format
        """

        import time, os
        from Scientific.IO.NetCDF import NetCDFFile

        self.domain.set_name('datatest' + str(id(self)))
        self.domain.format = 'sww'
        self.domain.smooth = True
        self.domain.reduction = mean

        sww = SWW_file(self.domain)
        sww.store_connectivity()
        sww.store_timestep()

        #Check contents
        #Get NetCDF
        fid = NetCDFFile(sww.filename, netcdf_mode_r)

        # Get the variables
        x = fid.variables['x']
        y = fid.variables['y']
        z = fid.variables['elevation']

        volumes = fid.variables['volumes']
        time = fid.variables['time']

        # 2D
        stage = fid.variables['stage']

        X = x[:]
        Y = y[:]
        Z = z[:]
        V = volumes[:]
        T = time[:]
        S = stage[:,:]

#         print "****************************"
#         print "X ",X
#         print "****************************"
#         print "Y ",Y
#         print "****************************"
#         print "Z ",Z
#         print "****************************"
#         print "V ",V
#         print "****************************"
#         print "Time ",T
#         print "****************************"
#         print "Stage ",S
#         print "****************************"


        fid.close()

        #Cleanup
        os.remove(sww.filename)



    def test_export_grid(self):
        """
        test_export_grid(self):
        Test that sww information can be converted correctly to asc/prj
        format readable by e.g. ArcView
        """

        import time, os
        from Scientific.IO.NetCDF import NetCDFFile

        try:
            os.remove('teg*.sww')
        except:
            pass

        #Setup
        self.domain.set_name('teg')

        prjfile = self.domain.get_name() + '_elevation.prj'
        ascfile = self.domain.get_name() + '_elevation.asc'
        swwfile = self.domain.get_name() + '.sww'

        self.domain.set_datadir('.')
        self.domain.smooth = True
        self.domain.set_quantity('elevation', lambda x,y: -x-y)
        self.domain.set_quantity('stage', 1.0)

        self.domain.geo_reference = Geo_reference(56,308500,6189000)

        sww = SWW_file(self.domain)
        sww.store_connectivity()
        sww.store_timestep()
        self.domain.evolve_to_end(finaltime = 0.01)
        sww.store_timestep()

        cellsize = 0.25
        #Check contents
        #Get NetCDF

        fid = NetCDFFile(sww.filename, netcdf_mode_r)

        # Get the variables
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        z = fid.variables['elevation'][:]
        time = fid.variables['time'][:]
        stage = fid.variables['stage'][:]

        fid.close()

        #Export to ascii/prj files
        export_grid(self.domain.get_name(),
                quantities = 'elevation',
                cellsize = cellsize,
                verbose = self.verbose,
                format = 'asc')

        #Check asc file
        ascid = open(ascfile)
        lines = ascid.readlines()
        ascid.close()

        L = lines[2].strip().split()
        assert L[0].strip().lower() == 'xllcorner'
        assert num.allclose(float(L[1].strip().lower()), 308500)

        L = lines[3].strip().split()
        assert L[0].strip().lower() == 'yllcorner'
        assert num.allclose(float(L[1].strip().lower()), 6189000)

        #Check grid values
        for j in range(5):
            L = lines[6+j].strip().split()
            y = (4-j) * cellsize
            for i in range(5):
                assert num.allclose(float(L[i]), -i*cellsize - y)
                
        #Cleanup
        os.remove(prjfile)
        os.remove(ascfile)
        os.remove(swwfile)

    def test_export_gridII(self):
        """
        test_export_gridII(self):
        Test that sww information can be converted correctly to asc/prj
        format readable by e.g. ArcView
        """

        import time, os
        from Scientific.IO.NetCDF import NetCDFFile

        try:
            os.remove('teg*.sww')
        except:
            pass

        #Setup
        self.domain.set_name('tegII')

        swwfile = self.domain.get_name() + '.sww'

        self.domain.set_datadir('.')
        self.domain.smooth = True
        self.domain.set_quantity('elevation', lambda x,y: -x-y)
        self.domain.set_quantity('stage', 1.0)

        self.domain.geo_reference = Geo_reference(56,308500,6189000)

        sww = SWW_file(self.domain)
        sww.store_connectivity()
        sww.store_timestep()
        self.domain.evolve_to_end(finaltime = 0.01)
        sww.store_timestep()

        cellsize = 0.25
        #Check contents
        #Get NetCDF

        fid = NetCDFFile(sww.filename, netcdf_mode_r)

        # Get the variables
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        z = fid.variables['elevation'][:]
        time = fid.variables['time'][:]
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        ymomentum = fid.variables['ymomentum'][:]        

        #print 'stage', stage
        #print 'xmom', xmomentum
        #print 'ymom', ymomentum        

        fid.close()

        #Export to ascii/prj files
        if True:
            export_grid(self.domain.get_name(),
                        quantities = ['elevation', 'depth'],
                        cellsize = cellsize,
                        verbose = self.verbose,
                        format = 'asc')

        else:
            export_grid(self.domain.get_name(),
                quantities = ['depth'],
                cellsize = cellsize,
                verbose = self.verbose,
                format = 'asc')


            export_grid(self.domain.get_name(),
                quantities = ['elevation'],
                cellsize = cellsize,
                verbose = self.verbose,
                format = 'asc')

        prjfile = self.domain.get_name() + '_elevation.prj'
        ascfile = self.domain.get_name() + '_elevation.asc'
        
        #Check asc file
        ascid = open(ascfile)
        lines = ascid.readlines()
        ascid.close()

        L = lines[2].strip().split()
        assert L[0].strip().lower() == 'xllcorner'
        assert num.allclose(float(L[1].strip().lower()), 308500)

        L = lines[3].strip().split()
        assert L[0].strip().lower() == 'yllcorner'
        assert num.allclose(float(L[1].strip().lower()), 6189000)

        #print "ascfile", ascfile
        #Check grid values
        for j in range(5):
            L = lines[6+j].strip().split()
            y = (4-j) * cellsize
            for i in range(5):
                #print " -i*cellsize - y",  -i*cellsize - y
                #print "float(L[i])", float(L[i])
                assert num.allclose(float(L[i]), -i*cellsize - y)

        #Cleanup
        os.remove(prjfile)
        os.remove(ascfile)
        
        #Check asc file
        ascfile = self.domain.get_name() + '_depth.asc'
        prjfile = self.domain.get_name() + '_depth.prj'
        ascid = open(ascfile)
        lines = ascid.readlines()
        ascid.close()

        L = lines[2].strip().split()
        assert L[0].strip().lower() == 'xllcorner'
        assert num.allclose(float(L[1].strip().lower()), 308500)

        L = lines[3].strip().split()
        assert L[0].strip().lower() == 'yllcorner'
        assert num.allclose(float(L[1].strip().lower()), 6189000)

        #Check grid values
        for j in range(5):
            L = lines[6+j].strip().split()
            y = (4-j) * cellsize
            for i in range(5):
                #print " -i*cellsize - y",  -i*cellsize - y
                #print "float(L[i])", float(L[i])                
                assert num.allclose(float(L[i]), 1 - (-i*cellsize - y))

        #Cleanup
        os.remove(prjfile)
        os.remove(ascfile)
        os.remove(swwfile)


    def test_export_gridIII(self):
        """
        test_export_gridIII
        Test that sww information can be converted correctly to asc/prj
        format readable by e.g. ArcView
        """

        import time, os
        from Scientific.IO.NetCDF import NetCDFFile

        try:
            os.remove('teg*.sww')
        except:
            pass

        #Setup
        
        self.domain.set_name('tegIII')

        swwfile = self.domain.get_name() + '.sww'

        self.domain.set_datadir('.')
        self.domain.format = 'sww'
        self.domain.smooth = True
        self.domain.set_quantity('elevation', lambda x,y: -x-y)
        self.domain.set_quantity('stage', 1.0)

        self.domain.geo_reference = Geo_reference(56,308500,6189000)
        
        sww = SWW_file(self.domain)
        sww.store_connectivity()
        sww.store_timestep() #'stage')
        self.domain.evolve_to_end(finaltime = 0.01)
        sww.store_timestep() #'stage')

        cellsize = 0.25
        #Check contents
        #Get NetCDF

        fid = NetCDFFile(sww.filename, netcdf_mode_r)

        # Get the variables
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        z = fid.variables['elevation'][:]
        time = fid.variables['time'][:]
        stage = fid.variables['stage'][:]

        fid.close()

        #Export to ascii/prj files
        extra_name_out = 'yeah'
        if True:
            export_grid(self.domain.get_name(),
                        quantities = ['elevation', 'depth'],
                        extra_name_out = extra_name_out,
                        cellsize = cellsize,
                        verbose = self.verbose,
                        format = 'asc')

        else:
            export_grid(self.domain.get_name(),
                quantities = ['depth'],
                cellsize = cellsize,
                verbose = self.verbose,
                format = 'asc')


            export_grid(self.domain.get_name(),
                quantities = ['elevation'],
                cellsize = cellsize,
                verbose = self.verbose,
                format = 'asc')

        prjfile = self.domain.get_name() + '_elevation_yeah.prj'
        ascfile = self.domain.get_name() + '_elevation_yeah.asc'
        
        #Check asc file
        ascid = open(ascfile)
        lines = ascid.readlines()
        ascid.close()

        L = lines[2].strip().split()
        assert L[0].strip().lower() == 'xllcorner'
        assert num.allclose(float(L[1].strip().lower()), 308500)

        L = lines[3].strip().split()
        assert L[0].strip().lower() == 'yllcorner'
        assert num.allclose(float(L[1].strip().lower()), 6189000)

        #print "ascfile", ascfile
        #Check grid values
        for j in range(5):
            L = lines[6+j].strip().split()
            y = (4-j) * cellsize
            for i in range(5):
                #print " -i*cellsize - y",  -i*cellsize - y
                #print "float(L[i])", float(L[i])
                assert num.allclose(float(L[i]), -i*cellsize - y)
                
        #Cleanup
        os.remove(prjfile)
        os.remove(ascfile)
        
        #Check asc file
        ascfile = self.domain.get_name() + '_depth_yeah.asc'
        prjfile = self.domain.get_name() + '_depth_yeah.prj'
        ascid = open(ascfile)
        lines = ascid.readlines()
        ascid.close()

        L = lines[2].strip().split()
        assert L[0].strip().lower() == 'xllcorner'
        assert num.allclose(float(L[1].strip().lower()), 308500)

        L = lines[3].strip().split()
        assert L[0].strip().lower() == 'yllcorner'
        assert num.allclose(float(L[1].strip().lower()), 6189000)

        #Check grid values
        for j in range(5):
            L = lines[6+j].strip().split()
            y = (4-j) * cellsize
            for i in range(5):
                assert num.allclose(float(L[i]), 1 - (-i*cellsize - y))

        #Cleanup
        os.remove(prjfile)
        os.remove(ascfile)
        os.remove(swwfile)

    def test_export_grid_bad(self):
        """Test that sww information can be converted correctly to asc/prj
        format readable by e.g. ArcView
        """

        try:
            export_grid('a_small_round-egg',
                        quantities = ['elevation', 'depth'],
                        cellsize = 99,
                        verbose = self.verbose,
                        format = 'asc')
        except IOError:
            pass
        else:
            self.failUnless(0 ==1,  'Bad input did not throw exception error!')

    def test_export_grid_parallel(self):
        """Test that sww information can be converted correctly to asc/prj
        format readable by e.g. ArcView
        """

        import time, os
        from Scientific.IO.NetCDF import NetCDFFile

        base_name = 'tegp'
        #Setup
        self.domain.set_name(base_name+'_P0_8')
        swwfile = self.domain.get_name() + '.sww'

        self.domain.set_datadir('.')
        self.domain.format = 'sww'
        self.domain.smooth = True
        self.domain.set_quantity('elevation', lambda x,y: -x-y)
        self.domain.set_quantity('stage', 1.0)

        self.domain.geo_reference = Geo_reference(56,308500,6189000)

        sww = SWW_file(self.domain)
        sww.store_connectivity()
        sww.store_timestep()
        self.domain.evolve_to_end(finaltime = 0.0001)
        #Setup
        self.domain.set_name(base_name+'_P1_8')
        swwfile2 = self.domain.get_name() + '.sww'
        sww = SWW_file(self.domain)
        sww.store_connectivity()
        sww.store_timestep()
        self.domain.evolve_to_end(finaltime = 0.0002)
        sww.store_timestep()

        cellsize = 0.25
        #Check contents
        #Get NetCDF

        fid = NetCDFFile(sww.filename, netcdf_mode_r)

        # Get the variables
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        z = fid.variables['elevation'][:]
        time = fid.variables['time'][:]
        stage = fid.variables['stage'][:]

        fid.close()

        #Export to ascii/prj files
        extra_name_out = 'yeah'
        export_grid(base_name,
                    quantities = ['elevation', 'depth'],
                    extra_name_out = extra_name_out,
                    cellsize = cellsize,
                    verbose = self.verbose,
                    format = 'asc')

        prjfile = base_name + '_P0_8_elevation_yeah.prj'
        ascfile = base_name + '_P0_8_elevation_yeah.asc'       
        #Check asc file
        ascid = open(ascfile)
        lines = ascid.readlines()
        ascid.close()
        #Check grid values
        for j in range(5):
            L = lines[6+j].strip().split()
            y = (4-j) * cellsize
            for i in range(5):
                #print " -i*cellsize - y",  -i*cellsize - y
                #print "float(L[i])", float(L[i])
                assert num.allclose(float(L[i]), -i*cellsize - y)               
        #Cleanup
        os.remove(prjfile)
        os.remove(ascfile)

        prjfile = base_name + '_P1_8_elevation_yeah.prj'
        ascfile = base_name + '_P1_8_elevation_yeah.asc'       
        #Check asc file
        ascid = open(ascfile)
        lines = ascid.readlines()
        ascid.close()
        #Check grid values
        for j in range(5):
            L = lines[6+j].strip().split()
            y = (4-j) * cellsize
            for i in range(5):
                #print " -i*cellsize - y",  -i*cellsize - y
                #print "float(L[i])", float(L[i])
                assert num.allclose(float(L[i]), -i*cellsize - y)               
        #Cleanup
        os.remove(prjfile)
        os.remove(ascfile)
        os.remove(swwfile)

        #Check asc file
        ascfile = base_name + '_P0_8_depth_yeah.asc'
        prjfile = base_name + '_P0_8_depth_yeah.prj'
        ascid = open(ascfile)
        lines = ascid.readlines()
        ascid.close()
        #Check grid values
        for j in range(5):
            L = lines[6+j].strip().split()
            y = (4-j) * cellsize
            for i in range(5):
                assert num.allclose(float(L[i]), 1 - (-i*cellsize - y))
        #Cleanup
        os.remove(prjfile)
        os.remove(ascfile)

        #Check asc file
        ascfile = base_name + '_P1_8_depth_yeah.asc'
        prjfile = base_name + '_P1_8_depth_yeah.prj'
        ascid = open(ascfile)
        lines = ascid.readlines()
        ascid.close()
        #Check grid values
        for j in range(5):
            L = lines[6+j].strip().split()
            y = (4-j) * cellsize
            for i in range(5):
                assert num.allclose(float(L[i]), 1 - (-i*cellsize - y))
        #Cleanup
        os.remove(prjfile)
        os.remove(ascfile)
        os.remove(swwfile2)



    def DISABLEDtest_sww2domain2(self):
        ##################################################################
        #Same as previous test, but this checks how NaNs are handled.
        ##################################################################

        #FIXME: See ticket 223

        from mesh_factory import rectangular

        #Create basic mesh
        points, vertices, boundary = rectangular(2,2)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.smooth = False
        domain.store = True
        domain.set_name('test_file')
        domain.set_datadir('.')
        domain.default_order=2

        domain.set_quantity('elevation', lambda x,y: -x/3)
        domain.set_quantity('friction', 0.1)

        from math import sin, pi
        Br = Reflective_boundary(domain)
        Bt = Transmissive_boundary(domain)
        Bd = Dirichlet_boundary([0.2,0.,0.])
        Bw = Time_boundary(domain=domain,
                           f=lambda t: [(0.1*sin(t*2*pi)), 0.0, 0.0])

        domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})

        h = 0.05
        elevation = domain.quantities['elevation'].vertex_values
        domain.set_quantity('stage', elevation + h)

        domain.check_integrity()

        for t in domain.evolve(yieldstep = 1, finaltime = 2.0):
            pass
            #domain.write_time()



        ##################################
        #Import the file as a new domain
        ##################################
        from data_manager import load_sww_as_domain
        import os

        filename = domain.datadir + os.sep + domain.get_name() + '.sww'

        # Fail because NaNs are present
        #domain2 = sww2domain(filename,
        #                     boundary,
        #                     fail_if_NaN=True,
        #                     verbose=self.verbose)        
        try:
            domain2 = load_sww_as_domain(filename,
                                 boundary,
                                 fail_if_NaN=True,
                                 verbose=self.verbose)
        except DataDomainError:
            # Now import it, filling NaNs to be -9999
            filler = -9999
            domain2 = load_sww_as_domain(filename,
                                 None,
                                 fail_if_NaN=False,
                                 NaN_filler=filler,
                                 verbose=self.verbose)
        else:
            raise Exception, 'should have failed'

            
        # Now import it, filling NaNs to be 0
        filler = -9999
        domain2 = load_sww_as_domain(filename,
                             None,
                             fail_if_NaN=False,
                             NaN_filler=filler,
                             verbose=self.verbose)            
                             
        import sys; sys.exit() 
            
        # Clean up
        os.remove(filename)


        bits = ['geo_reference.get_xllcorner()',
                'geo_reference.get_yllcorner()',
                'vertex_coordinates']

        for quantity in domain.quantities_to_be_stored:
            bits.append('get_quantity("%s").get_integral()' %quantity)
            bits.append('get_quantity("%s").get_values()' %quantity)

        for bit in bits:
        #    print 'testing that domain.'+bit+' has been restored'
            assert num.allclose(eval('domain.'+bit),eval('domain2.'+bit))

        print 
        print
        print domain2.get_quantity('xmomentum').get_values()
        print
        print domain2.get_quantity('stage').get_values()
        print
             
        print 'filler', filler
        print max(domain2.get_quantity('xmomentum').get_values().flat)
        
        assert max(max(domain2.get_quantity('xmomentum').get_values()))==filler
        assert min(min(domain2.get_quantity('xmomentum').get_values()))==filler
        assert max(max(domain2.get_quantity('ymomentum').get_values()))==filler
        assert min(min(domain2.get_quantity('ymomentum').get_values()))==filler



    #FIXME This fails - smooth makes the comparism too hard for allclose
    def ztest_sww2domain3(self):
        ################################################
        #DOMAIN.SMOOTH = TRUE !!!!!!!!!!!!!!!!!!!!!!!!!!
        ################################################
        from mesh_factory import rectangular
        #Create basic mesh

        yiel=0.01
        points, vertices, boundary = rectangular(10,10)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.geo_reference = Geo_reference(56,11,11)
        domain.smooth = True
        domain.store = True
        domain.set_name('bedslope')
        domain.default_order=2
        #Bed-slope and friction
        domain.set_quantity('elevation', lambda x,y: -x/3)
        domain.set_quantity('friction', 0.1)
        # Boundary conditions
        from math import sin, pi
        Br = Reflective_boundary(domain)
        Bt = Transmissive_boundary(domain)
        Bd = Dirichlet_boundary([0.2,0.,0.])
        Bw = Time_boundary(domain=domain,
                           f=lambda t: [(0.1*sin(t*2*pi)), 0.0, 0.0])

        domain.set_boundary({'left': Bd, 'right': Bd, 'top': Bd, 'bottom': Bd})

        domain.quantities_to_be_stored['xmomentum'] = 2
        domain.quantities_to_be_stored['ymomentum'] = 2        
        #Initial condition
        h = 0.05
        elevation = domain.quantities['elevation'].vertex_values
        domain.set_quantity('stage', elevation + h)


        domain.check_integrity()
        #Evolution
        for t in domain.evolve(yieldstep = yiel, finaltime = 0.05):
        #    domain.write_time()
            pass


        filename = domain.datadir + os.sep + domain.get_name() + '.sww'
        domain2 = load_sww_as_domain(filename,None,fail_if_NaN=False,verbose=self.verbose)
        #points, vertices, boundary = rectangular(15,15)
        #domain2.boundary = boundary
        ###################
        ##NOW TEST IT!!!
        ###################

        os.remove(domain.get_name() + '.sww')

        #FIXME smooth domain so that they can be compared


        bits = []
        for quantity in domain.quantities_to_be_stored:
            bits.append('quantities["%s"].get_integral()'%quantity)


        for bit in bits:
            #print 'testing that domain.'+bit+' has been restored'
            #print bit
            #print 'done'
            #print ('domain.'+bit), eval('domain.'+bit)
            #print ('domain2.'+bit), eval('domain2.'+bit)
            assert num.allclose(eval('domain.'+bit),eval('domain2.'+bit),rtol=1.0e-1,atol=1.e-3)
            pass

        ######################################
        #Now evolve them both, just to be sure
        ######################################x
        domain.time = 0.
        from time import sleep

        final = .5
        domain.set_quantity('friction', 0.1)
        domain.store = False
        domain.set_boundary({'left': Bd, 'right': Bd, 'top': Bd, 'bottom': Br})

        for t in domain.evolve(yieldstep = yiel, finaltime = final):
            #domain.write_time()
            pass

        domain2.smooth = True
        domain2.store = False
        domain2.default_order=2
        domain2.set_quantity('friction', 0.1)
        #Bed-slope and friction
        # Boundary conditions
        Bd2=Dirichlet_boundary([0.2,0.,0.])
        Br2 = Reflective_boundary(domain2)
        domain2.boundary = domain.boundary
        #print 'domain2.boundary'
        #print domain2.boundary
        domain2.set_boundary({'left': Bd2, 'right': Bd2, 'top': Bd2, 'bottom': Br2})
        #domain2.boundary = domain.boundary
        #domain2.set_boundary({'exterior': Bd})

        domain2.check_integrity()

        for t in domain2.evolve(yieldstep = yiel, finaltime = final):
            #domain2.write_time()
            pass

        ###################
        ##NOW TEST IT!!!
        ##################

        print '><><><><>>'
        bits = [ 'vertex_coordinates']

        for quantity in ['elevation','xmomentum','ymomentum']:
            #bits.append('quantities["%s"].get_integral()'%quantity)
            bits.append('get_quantity("%s").get_values()' %quantity)

        for bit in bits:
            print bit
            assert num.allclose(eval('domain.'+bit),eval('domain2.'+bit))


    def test_decimate_dem(self):
        """Test decimation of dem file
        """

        import os
        from Scientific.IO.NetCDF import NetCDFFile

        #Write test dem file
        root = 'decdemtest'

        filename = root + '.dem'
        fid = NetCDFFile(filename, netcdf_mode_w)

        fid.institution = 'Geoscience Australia'
        fid.description = 'NetCDF DEM format for compact and portable ' +\
                          'storage of spatial point data'

        nrows = 15
        ncols = 18

        fid.ncols = ncols
        fid.nrows = nrows
        fid.xllcorner = 2000.5
        fid.yllcorner = 3000.5
        fid.cellsize = 25
        fid.NODATA_value = -9999

        fid.zone = 56
        fid.false_easting = 0.0
        fid.false_northing = 0.0
        fid.projection = 'UTM'
        fid.datum = 'WGS84'
        fid.units = 'METERS'

        fid.createDimension('number_of_points', nrows*ncols)

        fid.createVariable('elevation', netcdf_float, ('number_of_points',))

        elevation = fid.variables['elevation']

        elevation[:] = (num.arange(nrows*ncols))

        fid.close()

        #generate the elevation values expected in the decimated file
        ref_elevation = [(  0+  1+  2+ 18+ 19+ 20+ 36+ 37+ 38) / 9.0,
                         (  4+  5+  6+ 22+ 23+ 24+ 40+ 41+ 42) / 9.0,
                         (  8+  9+ 10+ 26+ 27+ 28+ 44+ 45+ 46) / 9.0,
                         ( 12+ 13+ 14+ 30+ 31+ 32+ 48+ 49+ 50) / 9.0,
                         ( 72+ 73+ 74+ 90+ 91+ 92+108+109+110) / 9.0,
                         ( 76+ 77+ 78+ 94+ 95+ 96+112+113+114) / 9.0,
                         ( 80+ 81+ 82+ 98+ 99+100+116+117+118) / 9.0,
                         ( 84+ 85+ 86+102+103+104+120+121+122) / 9.0,
                         (144+145+146+162+163+164+180+181+182) / 9.0,
                         (148+149+150+166+167+168+184+185+186) / 9.0,
                         (152+153+154+170+171+172+188+189+190) / 9.0,
                         (156+157+158+174+175+176+192+193+194) / 9.0,
                         (216+217+218+234+235+236+252+253+254) / 9.0,
                         (220+221+222+238+239+240+256+257+258) / 9.0,
                         (224+225+226+242+243+244+260+261+262) / 9.0,
                         (228+229+230+246+247+248+264+265+266) / 9.0]

        # generate a stencil for computing the decimated values
        stencil = num.ones((3,3), num.float) / 9.0

        decimate_dem(root, stencil=stencil, cellsize_new=100)

        # Open decimated NetCDF file
        fid = NetCDFFile(root + '_100.dem', netcdf_mode_r)

        # Get decimated elevation
        elevation = fid.variables['elevation']

        # Check values
        assert num.allclose(elevation, ref_elevation)

        # Cleanup
        fid.close()

        os.remove(root + '.dem')
        os.remove(root + '_100.dem')

    def test_decimate_dem_NODATA(self):
        """Test decimation of dem file that includes NODATA values
        """

        import os
        from Scientific.IO.NetCDF import NetCDFFile

        # Write test dem file
        root = 'decdemtest'

        filename = root + '.dem'
        fid = NetCDFFile(filename, netcdf_mode_w)

        fid.institution = 'Geoscience Australia'
        fid.description = 'NetCDF DEM format for compact and portable ' +\
                          'storage of spatial point data'

        nrows = 15
        ncols = 18
        NODATA_value = -9999

        fid.ncols = ncols
        fid.nrows = nrows
        fid.xllcorner = 2000.5
        fid.yllcorner = 3000.5
        fid.cellsize = 25
        fid.NODATA_value = NODATA_value

        fid.zone = 56
        fid.false_easting = 0.0
        fid.false_northing = 0.0
        fid.projection = 'UTM'
        fid.datum = 'WGS84'
        fid.units = 'METERS'

        fid.createDimension('number_of_points', nrows*ncols)

        fid.createVariable('elevation', netcdf_float, ('number_of_points',))

        elevation = fid.variables['elevation']

        # Generate initial elevation values
        elevation_tmp = (num.arange(nrows*ncols))

        # Add some NODATA values
        elevation_tmp[0]   = NODATA_value
        elevation_tmp[95]  = NODATA_value
        elevation_tmp[188] = NODATA_value
        elevation_tmp[189] = NODATA_value
        elevation_tmp[190] = NODATA_value
        elevation_tmp[209] = NODATA_value
        elevation_tmp[252] = NODATA_value

        elevation[:] = elevation_tmp

        fid.close()

        # Generate the elevation values expected in the decimated file
        ref_elevation = [NODATA_value,
                         (  4+  5+  6+ 22+ 23+ 24+ 40+ 41+ 42) / 9.0,
                         (  8+  9+ 10+ 26+ 27+ 28+ 44+ 45+ 46) / 9.0,
                         ( 12+ 13+ 14+ 30+ 31+ 32+ 48+ 49+ 50) / 9.0,
                         ( 72+ 73+ 74+ 90+ 91+ 92+108+109+110) / 9.0,
                         NODATA_value,
                         ( 80+ 81+ 82+ 98+ 99+100+116+117+118) / 9.0,
                         ( 84+ 85+ 86+102+103+104+120+121+122) / 9.0,
                         (144+145+146+162+163+164+180+181+182) / 9.0,
                         (148+149+150+166+167+168+184+185+186) / 9.0,
                         NODATA_value,
                         (156+157+158+174+175+176+192+193+194) / 9.0,
                         NODATA_value,
                         (220+221+222+238+239+240+256+257+258) / 9.0,
                         (224+225+226+242+243+244+260+261+262) / 9.0,
                         (228+229+230+246+247+248+264+265+266) / 9.0]

        # Generate a stencil for computing the decimated values
        stencil = num.ones((3,3), num.float) / 9.0

        decimate_dem(root, stencil=stencil, cellsize_new=100)

        # Open decimated NetCDF file
        fid = NetCDFFile(root + '_100.dem', netcdf_mode_r)

        # Get decimated elevation
        elevation = fid.variables['elevation']

        # Check values
        assert num.allclose(elevation, ref_elevation)

        # Cleanup
        fid.close()

        os.remove(root + '.dem')
        os.remove(root + '_100.dem')      

        
    def test_file_boundary_stsIV_sinewave_ordering(self):
        """test_file_boundary_stsIV_sinewave_ordering(self):
        Read correct points from ordering file and apply sts to boundary
        This one uses a sine wave and compares to time boundary
        """
        
        from anuga.shallow_water import Domain
        from anuga.shallow_water import Reflective_boundary
        from anuga.shallow_water import Dirichlet_boundary
        from anuga.shallow_water import File_boundary
        from anuga.pmesh.mesh_interface import create_mesh_from_regions

        lat_long_points=[[6.01,97.0],[6.02,97.0],[6.05,96.9],[6.0,97.0]]
        bounding_polygon=[[6.0,97.0],[6.01,97.0],[6.02,97.0],[6.02,97.02],[6.00,97.02]]
        tide = 0.35
        time_step_count = 50
        time_step = 0.1
        times_ref = num.arange(0, time_step_count*time_step, time_step)
        
        n=len(lat_long_points)
        first_tstep=num.ones(n,num.int)
        last_tstep=(time_step_count)*num.ones(n,num.int)
        
        gauge_depth=20*num.ones(n,num.float)
        
        ha1=num.ones((n,time_step_count),num.float)
        ua1=3.*num.ones((n,time_step_count),num.float)
        va1=2.*num.ones((n,time_step_count),num.float)
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

        #print stage
        #print xmomentum
        #print ymomentum
        #print elevation
        
        

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
        temp_fbound=num.zeros(int(finaltime/yieldstep)+1,num.float)
    
        for i, t in enumerate(domain_fbound.evolve(yieldstep=yieldstep,
                                                   finaltime=finaltime, 
                                                   skip_initial_step=False)):
            temp_fbound[i]=domain_fbound.quantities['stage'].centroid_values[2]
    
        
        domain_time = Domain(meshname)
        domain_time.set_quantity('stage', tide)
        Br = Reflective_boundary(domain_time)
        Bw = Time_boundary(domain=domain_time,
                         f=lambda t: [num.sin(t)+tide,3.*(20.+num.sin(t)+tide),2.*(20.+num.sin(t)+tide)])
        domain_time.set_boundary({'ocean': Bw,'otherocean': Br})
        
        temp_time=num.zeros(int(finaltime/yieldstep)+1,num.float)
        for i, t in enumerate(domain_time.evolve(yieldstep=yieldstep,
                                                   finaltime=finaltime, 
                                                   skip_initial_step=False)):
            temp_time[i]=domain_time.quantities['stage'].centroid_values[2]



        #print temp_fbound
        #print temp_time

        #print domain_fbound.quantities['stage'].vertex_values
        #print domain_time.quantities['stage'].vertex_values
        
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
        
        

        
        
    def test_file_boundary_sts_time_limit(self):
        """test_file_boundary_stsIV_sinewave_ordering(self):
        Read correct points from ordering file and apply sts to boundary
        This one uses a sine wave and compares to time boundary
        
        This one tests that times used can be limited by upper limit
        """
        
        from anuga.shallow_water import Domain
        from anuga.shallow_water import Reflective_boundary
        from anuga.shallow_water import Dirichlet_boundary
        from anuga.shallow_water import File_boundary
        from anuga.pmesh.mesh_interface import create_mesh_from_regions

        lat_long_points=[[6.01,97.0],[6.02,97.0],[6.05,96.9],[6.0,97.0]]
        bounding_polygon=[[6.0,97.0],[6.01,97.0],[6.02,97.0],[6.02,97.02],[6.00,97.02]]
        tide = 0.35
        time_step_count = 50
        time_step = 0.1
        times_ref = num.arange(0, time_step_count*time_step, time_step)
        
        n=len(lat_long_points)
        first_tstep=num.ones(n,num.int)
        last_tstep=(time_step_count)*num.ones(n,num.int)
        
        gauge_depth=20*num.ones(n,num.float)
        
        ha1=num.ones((n,time_step_count),num.float)
        ua1=3.*num.ones((n,time_step_count),num.float)
        va1=2.*num.ones((n,time_step_count),num.float)
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
        starttime = fid.starttime[0]
        #print times
        #print starttime

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
        time_vec = Bf.F.get_time()
        assert num.allclose(Bf.F.starttime, starttime)
        assert num.allclose(time_vec, times_ref)                                   
        
        for time_limit in [0.1, 0.2, 0.5, 1.0, 2.2, 3.0, 4.3, 6.0, 10.0]:
            Bf = File_boundary(sts_file+'.sts', 
                               domain_fbound, 
                               time_limit=time_limit+starttime,
                               boundary_polygon=boundary_polygon)
        
            time_vec = Bf.F.get_time()
            assert num.allclose(Bf.F.starttime, starttime)            
            assert num.alltrue(time_vec < time_limit)
            
            
        try:    
            Bf = File_boundary(sts_file+'.sts', 
                               domain_fbound, 
                               time_limit=-1+starttime,
                               boundary_polygon=boundary_polygon)            
            time_vec = Bf.F.get_time()    
            print time_vec    
        except AssertionError:
            pass
        else:
            raise Exception, 'Should have raised Exception here'

    def test_lon_lat2grid(self):
        lonlatdep = [
            [ 113.06700134  ,  -26.06669998 ,   1.        ] ,
            [ 113.06700134  ,  -26.33329964 ,   3.        ] ,
            [ 113.19999695  ,  -26.06669998 ,   2.        ] ,
            [ 113.19999695  ,  -26.33329964 ,   4.        ] ]
            
        long, lat, quantity = lon_lat2grid(lonlatdep)

        for i, result in enumerate(lat):
            assert lonlatdep [i][1] == result
        assert len(lat) == 2 

        for i, result in enumerate(long):
            assert lonlatdep [i*2][0] == result
        assert len(long) == 2

        for i,q in enumerate(quantity):
            assert q == i+1
            
    def test_lon_lat2grid_bad(self):
        lonlatdep  = [
            [ -26.06669998,  113.06700134,    1.        ],
            [ -26.06669998 , 113.19999695 ,   2.        ],
            [ -26.06669998 , 113.33300018,    3.        ],
            [ -26.06669998 , 113.43299866   , 4.        ],
            [ -26.20000076 , 113.06700134,    5.        ],
            [ -26.20000076 , 113.19999695 ,   6.        ],
            [ -26.20000076 , 113.33300018  ,  7.        ],
            [ -26.20000076 , 113.43299866   , 8.        ],
            [ -26.33329964 , 113.06700134,    9.        ],
            [ -26.33329964 , 113.19999695 ,   10.        ],
            [ -26.33329964 , 113.33300018  ,  11.        ],
            [ -26.33329964 , 113.43299866 ,   12.        ],
            [ -26.43330002 , 113.06700134 ,   13        ],
            [ -26.43330002 , 113.19999695 ,   14.        ],
            [ -26.43330002 , 113.33300018,    15.        ],
            [ -26.43330002 , 113.43299866,    16.        ]]
        try:
            long, lat, quantity = lon_lat2grid(lonlatdep)
        except AssertionError:
            pass
        else:
            msg = 'Should have raised exception'
            raise msg
       
    def test_lon_lat2gridII(self):
        lonlatdep = [
            [ 113.06700134  ,  -26.06669998 ,   1.        ] ,
            [ 113.06700134  ,  -26.33329964 ,   2.        ] ,
            [ 113.19999695  ,  -26.06669998 ,   3.        ] ,
            [ 113.19999695  ,  -26.344329964 ,   4.        ] ]
        try:
            long, lat, quantity = lon_lat2grid(lonlatdep)
        except AssertionError:
            pass
        else:
            msg = 'Should have raised exception'
            raise msg
        
    #### END TESTS FOR URS 2 SWW  ###


    def test_triangulation(self):
        # 
        #  
        
        filename = tempfile.mktemp("_data_manager.sww")
        outfile = NetCDFFile(filename, netcdf_mode_w)
        points_utm = num.array([[0.,0.],[1.,1.], [0.,1.]])
        volumes = (0,1,2)
        elevation = [0,1,2]
        new_origin = None
        new_origin = Geo_reference(56, 0, 0)
        times = [0, 10]
        number_of_volumes = len(volumes)
        number_of_points = len(points_utm)
        sww = Write_sww(['elevation'], ['stage', 'xmomentum', 'ymomentum'])
        sww.store_header(outfile, times, number_of_volumes,
                         number_of_points, description='fully sick testing',
                         verbose=self.verbose,sww_precision=netcdf_float)
        sww.store_triangulation(outfile, points_utm, volumes,
                                elevation,  new_origin=new_origin,
                                verbose=self.verbose)       
        outfile.close()
        fid = NetCDFFile(filename)

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        fid.close()

        assert num.allclose(num.array(map(None, x,y)), points_utm)
        os.remove(filename)

        
    def test_triangulationII(self):
        # 
        #  
        
        filename = tempfile.mktemp("_data_manager.sww")
        outfile = NetCDFFile(filename, netcdf_mode_w)
        points_utm = num.array([[0.,0.],[1.,1.], [0.,1.]])
        volumes = (0,1,2)
        elevation = [0,1,2]
        new_origin = None
        #new_origin = Geo_reference(56, 0, 0)
        times = [0, 10]
        number_of_volumes = len(volumes)
        number_of_points = len(points_utm)
        sww = Write_sww(['elevation'], ['stage', 'xmomentum', 'ymomentum'])        
        sww.store_header(outfile, times, number_of_volumes,
                         number_of_points, description='fully sick testing',
                         verbose=self.verbose,sww_precision=netcdf_float)
        sww.store_triangulation(outfile, points_utm, volumes,
                                new_origin=new_origin,
                                verbose=self.verbose)
        sww.store_static_quantities(outfile, elevation=elevation)                                
                                
        outfile.close()
        fid = NetCDFFile(filename)

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        results_georef = Geo_reference()
        results_georef.read_NetCDF(fid)
        assert results_georef == Geo_reference(DEFAULT_ZONE, 0, 0)
        fid.close()

        assert num.allclose(num.array(map(None, x,y)), points_utm)
        os.remove(filename)

        
    def test_triangulation_new_origin(self):
        # 
        #  
        
        filename = tempfile.mktemp("_data_manager.sww")
        outfile = NetCDFFile(filename, netcdf_mode_w)
        points_utm = num.array([[0.,0.],[1.,1.], [0.,1.]])
        volumes = (0,1,2)
        elevation = [0,1,2]
        new_origin = None
        new_origin = Geo_reference(56, 1, 554354)
        points_utm = new_origin.change_points_geo_ref(points_utm)
        times = [0, 10]
        number_of_volumes = len(volumes)
        number_of_points = len(points_utm)
        sww = Write_sww(['elevation'], ['stage', 'xmomentum', 'ymomentum'])        
        sww.store_header(outfile, times, number_of_volumes,
                         number_of_points, description='fully sick testing',
                         verbose=self.verbose,sww_precision=netcdf_float)
        sww.store_triangulation(outfile, points_utm, volumes,
                                elevation,  new_origin=new_origin,
                                verbose=self.verbose)
        outfile.close()
        fid = NetCDFFile(filename)

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        results_georef = Geo_reference()
        results_georef.read_NetCDF(fid)
        assert results_georef == new_origin
        fid.close()

        absolute = Geo_reference(56, 0,0)
        assert num.allclose(num.array(
            absolute.change_points_geo_ref(map(None, x,y),
                                           new_origin)),points_utm)
        
        os.remove(filename)
        
    def test_triangulation_points_georeference(self):
        # 
        #  
        
        filename = tempfile.mktemp("_data_manager.sww")
        outfile = NetCDFFile(filename, netcdf_mode_w)
        points_utm = num.array([[0.,0.],[1.,1.], [0.,1.]])
        volumes = (0,1,2)
        elevation = [0,1,2]
        new_origin = None
        points_georeference = Geo_reference(56, 1, 554354)
        points_utm = points_georeference.change_points_geo_ref(points_utm)
        times = [0, 10]
        number_of_volumes = len(volumes)
        number_of_points = len(points_utm)
        sww = Write_sww(['elevation'], ['stage', 'xmomentum', 'ymomentum'])        
        sww.store_header(outfile, times, number_of_volumes,
                         number_of_points, description='fully sick testing',
                         verbose=self.verbose,sww_precision=netcdf_float)
        sww.store_triangulation(outfile, points_utm, volumes,
                                elevation,  new_origin=new_origin,
                                points_georeference=points_georeference,
                                verbose=self.verbose)       
        outfile.close()
        fid = NetCDFFile(filename)

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        results_georef = Geo_reference()
        results_georef.read_NetCDF(fid)
        assert results_georef == points_georeference
        fid.close()

        assert num.allclose(num.array(map(None, x,y)), points_utm)
        os.remove(filename)
        
    def test_triangulation_2_geo_refs(self):
        # 
        #  
        
        filename = tempfile.mktemp("_data_manager.sww")
        outfile = NetCDFFile(filename, netcdf_mode_w)
        points_utm = num.array([[0.,0.],[1.,1.], [0.,1.]])
        volumes = (0,1,2)
        elevation = [0,1,2]
        new_origin = Geo_reference(56, 1, 1)
        points_georeference = Geo_reference(56, 0, 0)
        points_utm = points_georeference.change_points_geo_ref(points_utm)
        times = [0, 10]
        number_of_volumes = len(volumes)
        number_of_points = len(points_utm)
        sww = Write_sww(['elevation'], ['stage', 'xmomentum', 'ymomentum'])        
        sww.store_header(outfile, times, number_of_volumes,
                         number_of_points, description='fully sick testing',
                         verbose=self.verbose,sww_precision=netcdf_float)
        sww.store_triangulation(outfile, points_utm, volumes,
                                elevation,  new_origin=new_origin,
                                points_georeference=points_georeference,
                                verbose=self.verbose)       
        outfile.close()
        fid = NetCDFFile(filename)

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        results_georef = Geo_reference()
        results_georef.read_NetCDF(fid)
        assert results_georef == new_origin
        fid.close()


        absolute = Geo_reference(56, 0,0)
        assert num.allclose(num.array(
            absolute.change_points_geo_ref(map(None, x,y),
                                           new_origin)),points_utm)
        os.remove(filename)
        

        
    def test_get_all_swwfiles(self):
        try:
            swwfiles = get_all_swwfiles('','test.txt')  #Invalid
        except IOError:
            pass
        else:
            raise 'Should have raised exception' 
        
    def test_get_all_swwfiles1(self):
        
        temp_dir = tempfile.mkdtemp('','sww_test')
        filename0 = tempfile.mktemp('.sww','test',temp_dir)
        filename1 = tempfile.mktemp('.sww','test',temp_dir)
        filename2 = tempfile.mktemp('.sww','test',temp_dir)
        filename3 = tempfile.mktemp('.sww','test',temp_dir)
       
        #print'filename', filename0,filename1,filename2,filename3
        
        fid0 = open(filename0, 'w')
        fid1 = open(filename1, 'w')
        fid2 = open(filename2, 'w')
        fid3 = open(filename3, 'w')

        fid0.write('hello')
        fid1.write('hello')
        fid2.write('hello')
        fid3.write('hello')
        
        fid0.close()
        fid1.close()
        fid2.close()
        fid3.close()
        
        
        dir, name0 = os.path.split(filename0)
        #print 'dir',dir,name0
        
        iterate=get_all_swwfiles(dir,'test')
        
        del_dir(temp_dir)
#        removeall(temp_dir)

        _, name0 = os.path.split(filename0) 
        #print'name0',name0[:-4],iterate[0]    
        _, name1 = os.path.split(filename1)       
        _, name2 = os.path.split(filename2)       
        _, name3 = os.path.split(filename3)       

        assert name0[:-4] in iterate
        assert name1[:-4] in iterate
        assert name2[:-4] in iterate
        assert name3[:-4] in iterate
        
        assert len(iterate)==4

 
    def test_points2polygon(self):  
        att_dict = {}
        pointlist = num.array([[1.0, 0.0],[0.0, 1.0],[0.0, 0.0]])
        att_dict['elevation'] = num.array([10.1, 0.0, 10.4])
        att_dict['brightness'] = num.array([10.0, 1.0, 10.4])
        
        fileName = tempfile.mktemp(".csv")
        
        G = Geospatial_data(pointlist, att_dict)
        
        G.export_points_file(fileName)
        
        polygon = points2polygon(fileName)
        
        # This test may fail if the order changes
        assert (polygon == [[0.0, 0.0],[1.0, 0.0],[0.0, 1.0]])
        
    
    def test_csv2polygons(self):
        """test_csv2polygons
        """
        
        path = get_pathname_from_package('anuga.shallow_water')                
        testfile = os.path.join(path, 'polygon_values_example.csv')                

        polygons, values = csv2polygons(testfile, 
                                        value_name='floors')

        assert len(polygons) == 7, 'Must have seven polygons'
        assert len(values) == 7, 'Must have seven values'
            
        # Known floor values
        floors = {'1': 2, '2': 0, '3': 1, '4': 2, '5': 0, '8': 1, '9': 1}
        
        # Known polygon values
        known_polys = {'1': [[422681.61,871117.55],
                             [422691.02,871117.60],
                             [422690.87,871084.23],
                             [422649.36,871081.85],
                             [422649.36,871080.39],
                             [422631.86,871079.50],
                             [422631.72,871086.75],
                             [422636.75,871087.20],
                             [422636.75,871091.50],
                             [422649.66,871092.09],
                             [422649.83,871084.91],
                             [422652.94,871084.90],
                             [422652.84,871092.39],
                             [422681.83,871093.73],
                             [422681.61,871117.55]],
                       '2': [[422664.22,870785.46],
                             [422672.48,870780.14],
                             [422668.17,870772.62],
                             [422660.35,870777.17],
                             [422664.22,870785.46]],
                       '3': [[422661.30,871215.06],
                             [422667.50,871215.70],
                             [422668.30,871204.86],
                             [422662.21,871204.33],
                             [422661.30,871215.06]],
                       '4': [[422473.44,871191.22],
                             [422478.33,871192.26],
                             [422479.52,871186.03],
                             [422474.78,871185.14],
                             [422473.44,871191.22]],
                       '5': [[422369.69,871049.29],
                             [422378.63,871053.58],
                             [422383.91,871044.51],
                             [422374.97,871040.32],
                             [422369.69,871049.29]],
                       '8': [[422730.56,871203.13],
                             [422734.10,871204.90],
                             [422735.26,871202.18],
                             [422731.87,871200.58],
                             [422730.56,871203.13]],
                       '9': [[422659.85,871213.80],
                             [422660.91,871210.97],
                             [422655.42,871208.85],
                             [422654.36,871211.68],
                             [422659.85,871213.80]]
                       }        
        

        
                
        for id in ['1', '2', '3', '4', '5' ,'8' ,'9']:
            assert id in polygons.keys()
            assert id in values.keys()            

            assert int(values[id]) == int(floors[id])
            assert len(polygons[id]) == len(known_polys[id])
            assert num.allclose(polygons[id], known_polys[id])


    def test_csv2polygons_with_clipping(self):
        """test_csv2polygons with optional clipping
        """
        #FIXME(Ole): Not Done!!
        
        path = get_pathname_from_package('anuga.shallow_water')                
        testfile = os.path.join(path, 'polygon_values_example.csv')                

        polygons, values = csv2polygons(testfile, 
                                        value_name='floors',
                                        clipping_polygons=None)

        assert len(polygons) == 7, 'Must have seven polygons'
        assert len(values) == 7, 'Must have seven values'
            
        # Known floor values
        floors = {'1': 2, '2': 0, '3': 1, '4': 2, '5': 0, '8': 1, '9': 1}
        
        # Known polygon values
        known_polys = {'1': [[422681.61,871117.55],
                             [422691.02,871117.60],
                             [422690.87,871084.23],
                             [422649.36,871081.85],
                             [422649.36,871080.39],
                             [422631.86,871079.50],
                             [422631.72,871086.75],
                             [422636.75,871087.20],
                             [422636.75,871091.50],
                             [422649.66,871092.09],
                             [422649.83,871084.91],
                             [422652.94,871084.90],
                             [422652.84,871092.39],
                             [422681.83,871093.73],
                             [422681.61,871117.55]],
                       '2': [[422664.22,870785.46],
                             [422672.48,870780.14],
                             [422668.17,870772.62],
                             [422660.35,870777.17],
                             [422664.22,870785.46]],
                       '3': [[422661.30,871215.06],
                             [422667.50,871215.70],
                             [422668.30,871204.86],
                             [422662.21,871204.33],
                             [422661.30,871215.06]],
                       '4': [[422473.44,871191.22],
                             [422478.33,871192.26],
                             [422479.52,871186.03],
                             [422474.78,871185.14],
                             [422473.44,871191.22]],
                       '5': [[422369.69,871049.29],
                             [422378.63,871053.58],
                             [422383.91,871044.51],
                             [422374.97,871040.32],
                             [422369.69,871049.29]],
                       '8': [[422730.56,871203.13],
                             [422734.10,871204.90],
                             [422735.26,871202.18],
                             [422731.87,871200.58],
                             [422730.56,871203.13]],
                       '9': [[422659.85,871213.80],
                             [422660.91,871210.97],
                             [422655.42,871208.85],
                             [422654.36,871211.68],
                             [422659.85,871213.80]]
                       }        
        

        
                
        for id in ['1', '2', '3', '4', '5' ,'8' ,'9']:
            assert id in polygons.keys()
            assert id in values.keys()            

            assert int(values[id]) == int(floors[id])
            assert len(polygons[id]) == len(known_polys[id])
            assert num.allclose(polygons[id], known_polys[id])




    
    def test_csv2building_polygons(self):
        """test_csv2building_polygons
        """
        
        path = get_pathname_from_package('anuga.shallow_water')                
        testfile = os.path.join(path, 'polygon_values_example.csv')                

        polygons, values = csv2building_polygons(testfile, 
                                                 floor_height=3)

        assert len(polygons) == 7, 'Must have seven polygons'
        assert len(values) == 7, 'Must have seven values'
            
        # Known floor values
        floors = {'1': 6, '2': 0, '3': 3, '4': 6, '5': 0, '8': 3, '9': 3}
        
                
        for id in ['1', '2', '3', '4', '5' ,'8' ,'9']:
            assert id in polygons.keys()
            assert id in values.keys()            

            assert float(values[id]) == float(floors[id])


#-------------------------------------------------------------

if __name__ == "__main__":
    #suite = unittest.makeSuite(Test_Data_Manager, 'test_sww2domain2')
    suite = unittest.makeSuite(Test_Data_Manager, 'test_sww')
    
    
    
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


