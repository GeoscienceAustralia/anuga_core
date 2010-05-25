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

from anuga.shallow_water.data_manager import *
from anuga.shallow_water.sww_file import SWW_file
from anuga.file_conversion.file_conversion import tsh2sww, \
                        asc_csiro2sww, pmesh_to_domain_instance
from anuga.coordinate_transforms.geo_reference import Geo_reference                        
from anuga.coordinate_transforms.redfearn import degminsec2decimal_degrees
from anuga.abstract_2d_finite_volumes.util import file_function
from anuga.utilities.system_tools import get_pathname_from_package
from anuga.utilities.file_utils import del_dir, load_csv_as_dict, \
                                        load_csv_as_array
from anuga.anuga_exceptions import ANUGAError
from anuga.utilities.numerical_tools import ensure_numeric, mean
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a
from anuga.config import netcdf_float, epsilon, g

# import all the boundaries - some are generic, some are shallow water
from boundaries import Reflective_boundary, \
            Field_boundary, Transmissive_momentum_set_stage_boundary, \
            Transmissive_stage_zero_momentum_boundary
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions\
     import Transmissive_boundary, Dirichlet_boundary, \
            Time_boundary, File_boundary, AWI_boundary

# This is needed to run the tests of local functions
import data_manager 
from anuga.coordinate_transforms.redfearn import redfearn
from anuga.coordinate_transforms.geo_reference import Geo_reference, \
     DEFAULT_ZONE
from anuga.geospatial_data.geospatial_data import Geospatial_data


class Test_Data_Manager(unittest.TestCase):
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


    def test_hecras_cross_sections2pts(self):
        """Test conversion from HECRAS cross sections in ascii format
        to native NetCDF pts format
        """

        import time, os
        from Scientific.IO.NetCDF import NetCDFFile

        #Write test asc file
        root = 'hecrastest'

        filename = root+'.sdf'
        fid = open(filename, 'w')
        fid.write("""
# RAS export file created on Mon 15Aug2005 11:42
# by HEC-RAS Version 3.1.1

BEGIN HEADER:
  UNITS: METRIC
  DTM TYPE: TIN
  DTM: v:\1\cit\perth_topo\river_tin
  STREAM LAYER: c:\\x_local\hecras\21_02_03\up_canning_cent3d.shp
  CROSS-SECTION LAYER: c:\\x_local\hecras\21_02_03\up_can_xs3d.shp
  MAP PROJECTION: UTM
  PROJECTION ZONE: 50
  DATUM: AGD66
  VERTICAL DATUM:
  NUMBER OF REACHES:  19
  NUMBER OF CROSS-SECTIONS:  2
END HEADER:


BEGIN CROSS-SECTIONS:

  CROSS-SECTION:
    STREAM ID:Southern-Wungong
    REACH ID:Southern-Wungong
    STATION:21410
    CUT LINE:
      407546.08 , 6437277.542
      407329.32 , 6437489.482
      407283.11 , 6437541.232
    SURFACE LINE:
     407546.08,   6437277.54,   52.14
     407538.88,   6437284.58,   51.07
     407531.68,   6437291.62,   50.56
     407524.48,   6437298.66,   49.58
     407517.28,   6437305.70,   49.09
     407510.08,   6437312.74,   48.76
  END:

  CROSS-SECTION:
    STREAM ID:Swan River
    REACH ID:Swan Mouth
    STATION:840.*
    CUT LINE:
      381178.0855 , 6452559.0685
      380485.4755 , 6453169.272
    SURFACE LINE:
     381178.09,   6452559.07,   4.17
     381169.49,   6452566.64,   4.26
     381157.78,   6452576.96,   4.34
     381155.97,   6452578.56,   4.35
     381143.72,   6452589.35,   4.43
     381136.69,   6452595.54,   4.58
     381114.74,   6452614.88,   4.41
     381075.53,   6452649.43,   4.17
     381071.47,   6452653.00,   3.99
     381063.46,   6452660.06,   3.67
     381054.41,   6452668.03,   3.67
  END:
END CROSS-SECTIONS:
""")

        fid.close()


        #Convert to NetCDF pts
        hecras_cross_sections2pts(root)

        #Check contents
        #Get NetCDF
        fid = NetCDFFile(root+'.pts', netcdf_mode_r)

        # Get the variables
        #print fid.variables.keys()
        points = fid.variables['points']
        elevation = fid.variables['elevation']

        #Check values
        ref_points = [[407546.08, 6437277.54],
                      [407538.88, 6437284.58],
                      [407531.68, 6437291.62],
                      [407524.48, 6437298.66],
                      [407517.28, 6437305.70],
                      [407510.08, 6437312.74]]

        ref_points += [[381178.09, 6452559.07],
                       [381169.49, 6452566.64],
                       [381157.78, 6452576.96],
                       [381155.97, 6452578.56],
                       [381143.72, 6452589.35],
                       [381136.69, 6452595.54],
                       [381114.74, 6452614.88],
                       [381075.53, 6452649.43],
                       [381071.47, 6452653.00],
                       [381063.46, 6452660.06],
                       [381054.41, 6452668.03]]


        ref_elevation = [52.14, 51.07, 50.56, 49.58, 49.09, 48.76]
        ref_elevation += [4.17, 4.26, 4.34, 4.35, 4.43, 4.58, 4.41, 4.17, 3.99, 3.67, 3.67]

        #print points[:]
        #print ref_points
        assert num.allclose(points, ref_points)

        #print attributes[:]
        #print ref_elevation
        assert num.allclose(elevation, ref_elevation)

        #Cleanup
        fid.close()


        os.remove(root + '.sdf')
        os.remove(root + '.pts')


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



    def test_sww2pts_centroids(self):
        """Test that sww information can be converted correctly to pts data at specified coordinates
        - in this case, the centroids.
        """

        import time, os
        from Scientific.IO.NetCDF import NetCDFFile
        # Used for points that lie outside mesh
        NODATA_value = 1758323

        # Setup
        self.domain.set_name('datatest')

        ptsfile = self.domain.get_name() + '_elevation.pts'
        swwfile = self.domain.get_name() + '.sww'

        self.domain.set_datadir('.')
        self.domain.format = 'sww'
        self.smooth = True #self.set_store_vertices_uniquely(False)
        self.domain.set_quantity('elevation', lambda x,y: -x-y)

        self.domain.geo_reference = Geo_reference(56,308500,6189000)

        sww = SWW_file(self.domain)
        sww.store_connectivity()
        sww.store_timestep()

        #self.domain.tight_slope_limiters = 1
        self.domain.evolve_to_end(finaltime = 0.01)
        sww.store_timestep()

        # Check contents in NetCDF
        fid = NetCDFFile(sww.filename, netcdf_mode_r)

        # Get the variables
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        elevation = fid.variables['elevation'][:]
        time = fid.variables['time'][:]
        stage = fid.variables['stage'][:]

        volumes = fid.variables['volumes'][:]


        # Invoke interpolation for vertex points       
        points = num.concatenate( (x[:,num.newaxis],y[:,num.newaxis]), axis=1 )
        points = num.ascontiguousarray(points)
        sww2pts(self.domain.get_name(),
                quantity = 'elevation',
                data_points = points,
                NODATA_value = NODATA_value,
                verbose = self.verbose)
        ref_point_values = elevation
        point_values = Geospatial_data(ptsfile).get_attributes()
        #print 'P', point_values
        #print 'Ref', ref_point_values        
        assert num.allclose(point_values, ref_point_values)        



        # Invoke interpolation for centroids
        points = self.domain.get_centroid_coordinates()
        #print points
        sww2pts(self.domain.get_name(),
                quantity = 'elevation',
                data_points = points,
                NODATA_value = NODATA_value,
                verbose = self.verbose)
        ref_point_values = [-0.5, -0.5, -1, -1, -1, -1, -1.5, -1.5]   #At centroids

        
        point_values = Geospatial_data(ptsfile).get_attributes()
        #print 'P', point_values
        #print 'Ref', ref_point_values        
        assert num.allclose(point_values, ref_point_values)        



        fid.close()

        #Cleanup
        os.remove(sww.filename)
        os.remove(ptsfile)



    def test_sww2domain1(self):
        ################################################
        #Create a test domain, and evolve and save it.
        ################################################
        from mesh_factory import rectangular

        #Create basic mesh

        yiel=0.01
        points, vertices, boundary = rectangular(10,10)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.geo_reference = Geo_reference(56,11,11)
        domain.smooth = False
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
        Bw = Time_boundary(domain=domain,f=lambda t: [(0.1*sin(t*2*pi)), 0.0, 0.0])

        #domain.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})
        domain.set_boundary({'left': Bd, 'right': Bd, 'top': Bd, 'bottom': Bd})

        domain.quantities_to_be_stored['xmomentum'] = 2
        domain.quantities_to_be_stored['ymomentum'] = 2
        #Initial condition
        h = 0.05
        elevation = domain.quantities['elevation'].vertex_values
        domain.set_quantity('stage', elevation + h)

        domain.check_integrity()
        #Evolution
        #domain.tight_slope_limiters = 1
        for t in domain.evolve(yieldstep = yiel, finaltime = 0.05):
            #domain.write_time()
            pass


        ##########################################
        #Import the example's file as a new domain
        ##########################################
        from data_manager import sww2domain
        import os

        filename = domain.datadir + os.sep + domain.get_name() + '.sww'
        domain2 = sww2domain(filename,None,fail_if_NaN=False,verbose=self.verbose)
        #points, vertices, boundary = rectangular(15,15)
        #domain2.boundary = boundary
        ###################
        ##NOW TEST IT!!!
        ###################

        os.remove(filename)

        bits = ['vertex_coordinates']
        for quantity in domain.quantities_to_be_stored:
            bits.append('get_quantity("%s").get_integral()' % quantity)
            bits.append('get_quantity("%s").get_values()' % quantity)

        for bit in bits:
            #print 'testing that domain.'+bit+' has been restored'
            #print bit
            #print 'done'
            assert num.allclose(eval('domain.'+bit),eval('domain2.'+bit))

        ######################################
        #Now evolve them both, just to be sure
        ######################################x
        domain.time = 0.
        from time import sleep

        final = .1
        domain.set_quantity('friction', 0.1)
        domain.store = False
        domain.set_boundary({'left': Bd, 'right': Bd, 'top': Bd, 'bottom': Bd})


        for t in domain.evolve(yieldstep = yiel, finaltime = final):
            #domain.write_time()
            pass

        final = final - (domain2.starttime-domain.starttime)
        #BUT since domain1 gets time hacked back to 0:
        final = final + (domain2.starttime-domain.starttime)

        domain2.smooth = False
        domain2.store = False
        domain2.default_order=2
        domain2.set_quantity('friction', 0.1)
        #Bed-slope and friction
        # Boundary conditions
        Bd2=Dirichlet_boundary([0.2,0.,0.])
        domain2.boundary = domain.boundary
        #print 'domain2.boundary'
        #print domain2.boundary
        domain2.set_boundary({'left': Bd, 'right': Bd, 'top': Bd, 'bottom': Bd})
        #domain2.set_boundary({'exterior': Bd})

        domain2.check_integrity()

        for t in domain2.evolve(yieldstep = yiel, finaltime = final):
            #domain2.write_time()
            pass

        ###################
        ##NOW TEST IT!!!
        ##################

        bits = ['vertex_coordinates']

        for quantity in ['elevation','stage', 'ymomentum','xmomentum']:
            bits.append('get_quantity("%s").get_integral()' %quantity)
            bits.append('get_quantity("%s").get_values()' %quantity)

        #print bits
        for bit in bits:
            #print bit
            #print eval('domain.'+bit)
            #print eval('domain2.'+bit)
            
            #print eval('domain.'+bit+'-domain2.'+bit)
            msg = 'Values in the two domains are different for ' + bit
            assert num.allclose(eval('domain.'+bit),eval('domain2.'+bit),
                                rtol=1.e-5, atol=3.e-8), msg


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
        from data_manager import sww2domain
        import os

        filename = domain.datadir + os.sep + domain.get_name() + '.sww'

        # Fail because NaNs are present
        #domain2 = sww2domain(filename,
        #                     boundary,
        #                     fail_if_NaN=True,
        #                     verbose=self.verbose)        
        try:
            domain2 = sww2domain(filename,
                                 boundary,
                                 fail_if_NaN=True,
                                 verbose=self.verbose)
        except DataDomainError:
            # Now import it, filling NaNs to be -9999
            filler = -9999
            domain2 = sww2domain(filename,
                                 None,
                                 fail_if_NaN=False,
                                 NaN_filler=filler,
                                 verbose=self.verbose)
        else:
            raise Exception, 'should have failed'

            
        # Now import it, filling NaNs to be 0
        filler = -9999
        domain2 = sww2domain(filename,
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



    def test_weed(self):
        from data_manager import weed

        coordinates1 = [[0.,0.],[1.,0.],[1.,1.],[1.,0.],[2.,0.],[1.,1.]]
        volumes1 = [[0,1,2],[3,4,5]]
        boundary1= {(0,1): 'external',(1,2): 'not external',(2,0): 'external',(3,4): 'external',(4,5): 'external',(5,3): 'not external'}
        coordinates2,volumes2,boundary2=weed(coordinates1,volumes1,boundary1)

        points2 = {(0.,0.):None,(1.,0.):None,(1.,1.):None,(2.,0.):None}

        assert len(points2)==len(coordinates2)
        for i in range(len(coordinates2)):
            coordinate = tuple(coordinates2[i])
            assert points2.has_key(coordinate)
            points2[coordinate]=i

        for triangle in volumes1:
            for coordinate in triangle:
                assert coordinates2[points2[tuple(coordinates1[coordinate])]][0]==coordinates1[coordinate][0]
                assert coordinates2[points2[tuple(coordinates1[coordinate])]][1]==coordinates1[coordinate][1]


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


        ##########################################
        #Import the example's file as a new domain
        ##########################################
        from data_manager import sww2domain
        import os

        filename = domain.datadir + os.sep + domain.get_name() + '.sww'
        domain2 = sww2domain(filename,None,fail_if_NaN=False,verbose=self.verbose)
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


    def write_mux(self, lat_long_points, time_step_count, time_step,
                  depth=None, ha=None, ua=None, va=None):
        """
        This will write 3 non-gridded mux files, for testing.
        If no quantities are passed in,
        na and va quantities will be the Easting values.
        Depth and ua will be the Northing value.
        
        The mux file format has south as positive so
        this function will swap the sign for va.  
        """

        #print "lat_long_points", lat_long_points
        #print "time_step_count",time_step_count
        #print "time_step",

        
        points_num = len(lat_long_points)
        lonlatdeps = []
        quantities = ['HA','UA','VA']
        
        mux_names = [WAVEHEIGHT_MUX_LABEL,
                     EAST_VELOCITY_LABEL,
                     NORTH_VELOCITY_LABEL]
        quantities_init = [[],[],[]]
        # urs binary is latitude fastest
        for point in lat_long_points:
            lat = point[0]
            lon = point[1]
            _ , e, n = redfearn(lat, lon)
            if depth is None:
                this_depth = n
            else:
                this_depth = depth
            if ha is None:
                this_ha = e
            else:
                this_ha = ha
            if ua is None:
                this_ua = n
            else:
                this_ua = ua
            if va is None:
                this_va = e   
            else:
                this_va = va         
            lonlatdeps.append([lon, lat, this_depth])
            quantities_init[0].append(this_ha) # HA
            quantities_init[1].append(this_ua) # UA
            quantities_init[2].append(this_va) # VA 
                
        file_handle, base_name = tempfile.mkstemp("")
        os.close(file_handle)
        os.remove(base_name)

        files = []        
        for i, q in enumerate(quantities): 
            quantities_init[i] = ensure_numeric(quantities_init[i])
            #print "HA_init", HA_init
            q_time = num.zeros((time_step_count, points_num), num.float)
            for time in range(time_step_count):
                q_time[time,:] = quantities_init[i] #* time * 4
            
            #Write C files
            columns = 3 # long, lat , depth
            file = base_name + mux_names[i]
            #print "base_name file",file 
            f = open(file, 'wb')
            files.append(file)
            f.write(pack('i',points_num))
            f.write(pack('i',time_step_count))
            f.write(pack('f',time_step))

            #write lat/long info
            for lonlatdep in lonlatdeps:
                for float in lonlatdep:
                    f.write(pack('f',float))
                    
            # Write quantity info
            for time in  range(time_step_count):
                for point_i in range(points_num):
                    f.write(pack('f',q_time[time,point_i]))
                    #print " mux_names[i]", mux_names[i] 
                    #print "f.write(pack('f',q_time[time,i]))", q_time[time,point_i]
            f.close()
        return base_name, files
        
    
    def delete_mux(self, files):
        for file in files:
            os.remove(file)
        
    def write_mux2(self, lat_long_points, time_step_count, time_step,
                   first_tstep, last_tstep,
                   depth=None, ha=None, ua=None, va=None):
        """
        This will write 3 non-gridded mux files, for testing.
        If no quantities are passed in,
        na and va quantities will be the Easting values.
        Depth and ua will be the Northing value.
        """
        #print "lat_long_points", lat_long_points
        #print "time_step_count",time_step_count
        #print "time_step",

        #irrelevant header information
        ig=ilon=ilat=0
        mcolat=mcolon=centerlat=centerlon=offset=az=baz=id=0.0

        points_num = len(lat_long_points)
        latlondeps = []
        quantities = ['HA','UA','VA']

        mux_names = [WAVEHEIGHT_MUX2_LABEL,
                     EAST_VELOCITY_MUX2_LABEL,
                     NORTH_VELOCITY_MUX2_LABEL]

        msg='first_tstep and last_step arrays must have same length as number of points'
        assert len(first_tstep)==points_num,msg
        assert len(last_tstep)==points_num,msg

        if depth is not None:
            depth=ensure_numeric(depth)
            assert len(depth)==points_num
        if ha is not None:
            ha=ensure_numeric(ha)
            assert ha.shape==(points_num,time_step_count)
        if ua is not None:
            ua=ensure_numeric(ua)
            assert ua.shape==(points_num,time_step_count)
        if va is not None:
            va=ensure_numeric(va)
            assert va.shape==(points_num,time_step_count)

        quantities_init = [[],[],[]]
        # urs binary is latitude fastest
        for i,point in enumerate(lat_long_points):
            lat = point[0]
            lon = point[1]
            _ , e, n = redfearn(lat, lon)
            if depth is None:
                this_depth = n
            else:
                this_depth = depth[i]
            latlondeps.append([lat, lon, this_depth])

            if ha is None:
                this_ha = e
                quantities_init[0].append(num.ones(time_step_count,num.float)*this_ha) # HA
            else:
                quantities_init[0].append(ha[i])
            if ua is None:
                this_ua = n
                quantities_init[1].append(num.ones(time_step_count,num.float)*this_ua) # UA
            else:
                quantities_init[1].append(ua[i])
            if va is None:
                this_va = e
                quantities_init[2].append(num.ones(time_step_count,num.float)*this_va) #
            else:
                quantities_init[2].append(-va[i]) # South is negative in MUX

        file_handle, base_name = tempfile.mkstemp("write_mux2")
        os.close(file_handle)
        os.remove(base_name)

        files = []        
        for i, q in enumerate(quantities):
            q_time = num.zeros((time_step_count, points_num), num.float)
            quantities_init[i] = ensure_numeric(quantities_init[i])
            for time in range(time_step_count):
                #print i, q, time, quantities_init[i][:,time]
                q_time[time,:] = quantities_init[i][:,time]
                #print i, q, time, q_time[time, :]                

            #Write C files
            columns = 3 # long, lat , depth
            file = base_name + mux_names[i]
            
            #print 'base_name file', file 
            f = open(file, 'wb')
            files.append(file)

            f.write(pack('i',points_num))
            #write mux 2 header
            for latlondep in latlondeps:
                f.write(pack('f',latlondep[0]))
                f.write(pack('f',latlondep[1]))
                f.write(pack('f',mcolat))
                f.write(pack('f',mcolon))
                f.write(pack('i',ig))
                f.write(pack('i',ilon))
                f.write(pack('i',ilat))
                f.write(pack('f',latlondep[2]))
                f.write(pack('f',centerlat))
                f.write(pack('f',centerlon))
                f.write(pack('f',offset))
                f.write(pack('f',az))
                f.write(pack('f',baz))
                f.write(pack('f',time_step))
                f.write(pack('i',time_step_count))
                for j in range(4): # identifier
                    f.write(pack('f',id))    

            #first_tstep=1
            #last_tstep=time_step_count
            for i,latlondep in enumerate(latlondeps):
                f.write(pack('i',first_tstep[i]))
            for i,latlondep in enumerate(latlondeps):
                f.write(pack('i',last_tstep[i]))

            # Find when first station starts recording
            min_tstep = min(first_tstep)
            # Find when all stations have stopped recording
            max_tstep = max(last_tstep)

            #for time in  range(time_step_count):
            for time in range(min_tstep-1,max_tstep):
                f.write(pack('f',time*time_step))                
                for point_i in range(points_num):
                    if time+1>=first_tstep[point_i] and time+1<=last_tstep[point_i]:
                        #print 'writing', time, point_i, q_time[time, point_i]
                        f.write(pack('f', q_time[time, point_i]))
            f.close()

        return base_name, files

    def test_urs2sts_read_mux2_pyI(self):
        """test_urs2sts_read_mux2_pyI(self):
        Constant stage,momentum at each gauge
        """
        tide = 1
        time_step_count = 3
        time_step = 2
        lat_long_points =[(-21.5,114.5),(-21,114.5),(-21.5,115), (-21.,115.)]
        n=len(lat_long_points)
        first_tstep=num.ones(n,num.int)
        last_tstep=time_step_count*num.ones(n,num.int)
        depth=20*num.ones(n,num.float)
        ha=2*num.ones((n,time_step_count),num.float)
        ua=5*num.ones((n,time_step_count),num.float)
        va=-10*num.ones((n,time_step_count),num.float)
        #-ve added to take into account mux file format where south is positive.
        base_name, files = self.write_mux2(lat_long_points,
                                      time_step_count, time_step,
                                      first_tstep, last_tstep,
                                      depth=depth,
                                      ha=ha,
                                      ua=ua,
                                      va=va)

        weights=num.ones(1, num.float)
        #ensure that files are indeed mux2 files
        times, latitudes, longitudes, elevation, stage, starttime = read_mux2_py([files[0]],weights)
        ua_times, ua_latitudes, ua_longitudes, ua_elevation, xvelocity,starttime_ua=read_mux2_py([files[1]],weights)
        msg='ha and ua have different gauge meta data'
        assert num.allclose(times,ua_times) and num.allclose(latitudes,ua_latitudes) and num.allclose(longitudes,ua_longitudes) and num.allclose(elevation,ua_elevation) and num.allclose(starttime,starttime_ua), msg
        va_times, va_latitudes, va_longitudes, va_elevation, yvelocity, starttime_va=read_mux2_py([files[2]],weights)
        msg='ha and va have different gauge meta data'
        assert num.allclose(times,va_times) and num.allclose(latitudes,va_latitudes) and num.allclose(longitudes,va_longitudes) and num.allclose(elevation,va_elevation) and num.allclose(starttime,starttime_va), msg

        self.delete_mux(files)

        msg='time array has incorrect length'
        assert times.shape[0]==time_step_count,msg
        
        msg = 'time array is incorrect'
        #assert allclose(times,time_step*num.arange(1,time_step_count+1)),msg
        assert num.allclose(times,time_step*num.arange(time_step_count)), msg
        
        msg='Incorrect gauge positions returned'
        for i,point in enumerate(lat_long_points):
            assert num.allclose(latitudes[i],point[0]) and num.allclose(longitudes[i],point[1]),msg

        msg='Incorrect gauge depths returned'
        assert num.allclose(elevation,-depth),msg
        msg='incorrect gauge height time series returned'
        assert num.allclose(stage,ha)
        msg='incorrect gauge ua time series returned'
        assert num.allclose(xvelocity,ua)
        msg='incorrect gauge va time series returned'
        assert num.allclose(yvelocity, -va)

    def test_urs2sts_read_mux2_pyII(self):
        """Spatially varing stage
        """
        tide = 1
        time_step_count = 3
        time_step = 2
        lat_long_points =[(-21.5,114.5),(-21,114.5),(-21.5,115), (-21.,115.)]
        n=len(lat_long_points)
        first_tstep=num.ones(n,num.int)
        last_tstep=(time_step_count)*num.ones(n,num.int)
        depth=20*num.ones(n,num.float)
        ha=2*num.ones((n,time_step_count),num.float)
        ha[0]=num.arange(0,time_step_count)+1
        ha[1]=time_step_count-num.arange(1,time_step_count+1)
        ha[1]=num.arange(time_step_count,2*time_step_count)
        ha[2]=num.arange(2*time_step_count,3*time_step_count)
        ha[3]=num.arange(3*time_step_count,4*time_step_count)
        ua=5*num.ones((n,time_step_count),num.float)
        va=-10*num.ones((n,time_step_count),num.float)
        #-ve added to take into account mux file format where south is positive.
        base_name, files = self.write_mux2(lat_long_points,
                                      time_step_count, time_step,
                                      first_tstep, last_tstep,
                                      depth=depth,
                                      ha=ha,
                                      ua=ua,
                                      va=va)

        weights=num.ones(1, num.float)
        #ensure that files are indeed mux2 files
        times, latitudes, longitudes, elevation, stage,starttime=read_mux2_py([files[0]],weights)
        ua_times, ua_latitudes, ua_longitudes, ua_elevation, xvelocity,starttime_ua=read_mux2_py([files[1]],weights)
        msg='ha and ua have different gauge meta data'
        assert num.allclose(times,ua_times) and num.allclose(latitudes,ua_latitudes) and num.allclose(longitudes,ua_longitudes) and num.allclose(elevation,ua_elevation) and num.allclose(starttime,starttime_ua), msg
        va_times, va_latitudes, va_longitudes, va_elevation, yvelocity,starttime_va=read_mux2_py([files[2]],weights)
        msg='ha and va have different gauge meta data'
        assert num.allclose(times,va_times) and num.allclose(latitudes,va_latitudes) and num.allclose(longitudes,va_longitudes) and num.allclose(elevation,va_elevation) and num.allclose(starttime,starttime_va), msg


        self.delete_mux(files)

        msg='time array has incorrect length'
        #assert times.shape[0]==time_step_count,msg
        msg = 'time array is incorrect'
        #assert allclose(times,time_step*num.arange(1,time_step_count+1)),msg
        msg='Incorrect gauge positions returned'
        for i,point in enumerate(lat_long_points):
            assert num.allclose(latitudes[i],point[0]) and num.allclose(longitudes[i],point[1]),msg

        msg='Incorrect gauge depths returned'
        assert num.allclose(elevation, -depth),msg
        msg='incorrect gauge height time series returned'
        assert num.allclose(stage, ha)
        msg='incorrect gauge ua time series returned'
        assert num.allclose(xvelocity, ua)
        msg='incorrect gauge va time series returned'
        assert num.allclose(yvelocity, -va) # South is positive in MUX

    def test_urs2sts_read_mux2_pyIII(self):
        """Varying start and finish times
        """
        tide = 1
        time_step_count = 3
        time_step = 2
        lat_long_points =[(-21.5,114.5),(-21,114.5),(-21.5,115), (-21.,115.)]
        n=len(lat_long_points)
        first_tstep=num.ones(n,num.int)
        first_tstep[0]+=1
        first_tstep[2]+=1
        last_tstep=(time_step_count)*num.ones(n,num.int)
        last_tstep[0]-=1

        depth=20*num.ones(n,num.float)
        ha=2*num.ones((n,time_step_count),num.float)
        ha[0]=num.arange(0,time_step_count)
        ha[1]=num.arange(time_step_count,2*time_step_count)
        ha[2]=num.arange(2*time_step_count,3*time_step_count)
        ha[3]=num.arange(3*time_step_count,4*time_step_count)
        ua=5*num.ones((n,time_step_count),num.float)
        va=-10*num.ones((n,time_step_count),num.float)
        #-ve added to take into account mux file format where south is positive.
        base_name, files = self.write_mux2(lat_long_points,
                                      time_step_count, time_step,
                                      first_tstep, last_tstep,
                                      depth=depth,
                                      ha=ha,
                                      ua=ua,
                                      va=va)

        weights=num.ones(1, num.float)
        #ensure that files are indeed mux2 files
        times, latitudes, longitudes, elevation, stage, starttime=read_mux2_py([files[0]],weights)
        ua_times, ua_latitudes, ua_longitudes, ua_elevation, xvelocity, starttime_ua=read_mux2_py([files[1]],weights)
        msg='ha and ua have different gauge meta data'
        assert num.allclose(times,ua_times) and num.allclose(latitudes,ua_latitudes) and num.allclose(longitudes,ua_longitudes) and num.allclose(elevation,ua_elevation) and num.allclose(starttime,starttime_ua), msg
        va_times, va_latitudes, va_longitudes, va_elevation, yvelocity,starttime_va=read_mux2_py([files[2]],weights)
        msg='ha and va have different gauge meta data'
        assert num.allclose(times,va_times) and num.allclose(latitudes,va_latitudes) and num.allclose(longitudes,va_longitudes) and num.allclose(elevation,va_elevation) and num.allclose(starttime,starttime_va), msg

        self.delete_mux(files)

        msg='time array has incorrect length'
        #assert times.shape[0]==time_step_count,msg
        msg = 'time array is incorrect'
        #assert allclose(times,time_step*num.arange(1,time_step_count+1)),msg
        msg='Incorrect gauge positions returned'
        for i,point in enumerate(lat_long_points):
            assert num.allclose(latitudes[i],point[0]) and num.allclose(longitudes[i],point[1]),msg


        # Set original data used to write mux file to be zero when gauges are 
        #not recdoring
        ha[0][0]=0.0
        ha[0][time_step_count-1]=0.0;
        ha[2][0]=0.0;
        ua[0][0]=0.0
        ua[0][time_step_count-1]=0.0;
        ua[2][0]=0.0;
        va[0][0]=0.0
        va[0][time_step_count-1]=0.0;
        va[2][0]=0.0;
        msg='Incorrect gauge depths returned'
        assert num.allclose(elevation,-depth),msg
        msg='incorrect gauge height time series returned'
        assert num.allclose(stage,ha)
        msg='incorrect gauge ua time series returned'
        assert num.allclose(xvelocity,ua)
        msg='incorrect gauge va time series returned'
        assert num.allclose(yvelocity, -va) # South is positive in mux
        

        
    def test_read_mux_platform_problem1(self):
        """test_read_mux_platform_problem1
        
        This is to test a situation where read_mux returned 
        wrong values Win32

        This test passes on Windows but test_read_mux_platform_problem2
        does not
        """
        
        from urs_ext import read_mux2 
        
        verbose = False
                
        tide = 1.5
        time_step_count = 10
        time_step = 0.2
        times_ref = num.arange(0, time_step_count*time_step, time_step)

        lat_long_points = [(-21.5,114.5), (-21,114.5), (-21.5,115), (-21.,115.), (-22., 117.)]
        n = len(lat_long_points)
        
        # Create different timeseries starting and ending at different times 
        first_tstep=num.ones(n, num.int)
        first_tstep[0]+=2   # Point 0 starts at 2
        first_tstep[1]+=4   # Point 1 starts at 4        
        first_tstep[2]+=3   # Point 2 starts at 3
        
        last_tstep=(time_step_count)*num.ones(n,num.int)
        last_tstep[0]-=1    # Point 0 ends 1 step early
        last_tstep[1]-=2    # Point 1 ends 2 steps early                
        last_tstep[4]-=3    # Point 4 ends 3 steps early        
        
        # Create varying elevation data (positive values for seafloor)
        gauge_depth=20*num.ones(n,num.float)
        for i in range(n):
            gauge_depth[i] += i**2
            
        # Create data to be written to first mux file        
        ha0=2*num.ones((n,time_step_count),num.float)
        ha0[0]=num.arange(0,time_step_count)
        ha0[1]=num.arange(time_step_count,2*time_step_count)
        ha0[2]=num.arange(2*time_step_count,3*time_step_count)
        ha0[3]=num.arange(3*time_step_count,4*time_step_count)
        ua0=5*num.ones((n,time_step_count),num.float)
        va0=-10*num.ones((n,time_step_count),num.float)

        # Ensure data used to write mux file to be zero when gauges are
        # not recording
        for i in range(n):
             # For each point
             for j in range(0, first_tstep[i]-1) + range(last_tstep[i], time_step_count):
                 # For timesteps before and after recording range
                 ha0[i][j] = ua0[i][j] = va0[i][j] = 0.0                                  
        
        # Write first mux file to be combined by urs2sts
        base_nameI, filesI = self.write_mux2(lat_long_points,
                                             time_step_count, time_step,
                                             first_tstep, last_tstep,
                                             depth=gauge_depth,
                                             ha=ha0,
                                             ua=ua0,
                                             va=va0)

        # Create ordering file
        permutation = ensure_numeric([4,0,2])

        _, ordering_filename = tempfile.mkstemp('')
        order_fid = open(ordering_filename, 'w')  
        order_fid.write('index, longitude, latitude\n')
        for index in permutation:
            order_fid.write('%d, %f, %f\n' %(index, 
                                             lat_long_points[index][1], 
                                             lat_long_points[index][0]))
        order_fid.close()
        
        

        # -------------------------------------
        # Now read files back and check values
        weights = ensure_numeric([1.0])

        # For each quantity read the associated list of source mux2 file with 
        # extention associated with that quantity
        file_params=-1*num.ones(3,num.float) #[nsta,dt,nt]
        OFFSET = 5

        for j, file in enumerate(filesI):
            data = read_mux2(1, [file], weights, file_params, permutation, verbose)

            number_of_selected_stations = data.shape[0]

            # Index where data ends and parameters begin
            parameters_index = data.shape[1]-OFFSET          
          
            for i in range(number_of_selected_stations):
                if j == 0: assert num.allclose(data[i][:parameters_index], ha0[permutation[i], :])
                if j == 1: assert num.allclose(data[i][:parameters_index], ua0[permutation[i], :])
                if j == 2: assert num.allclose(data[i][:parameters_index], -va0[permutation[i], :])
        
        self.delete_mux(filesI)


        
        
    def test_read_mux_platform_problem2(self):
        """test_read_mux_platform_problem2
        
        This is to test a situation where read_mux returned 
        wrong values Win32

        This test does not pass on Windows but test_read_mux_platform_problem1
        does
        """
        
        from urs_ext import read_mux2 
        
        from anuga.config import single_precision as epsilon        
        
        verbose = False
                
        tide = 1.5
        time_step_count = 10
        time_step = 0.2
        
        times_ref = num.arange(0, time_step_count*time_step, time_step)
        
        lat_long_points = [(-21.5,114.5), (-21,114.5), (-21.5,115),
                           (-21.,115.), (-22., 117.)]
        n = len(lat_long_points)
        
        # Create different timeseries starting and ending at different times 
        first_tstep=num.ones(n,num.int)
        first_tstep[0]+=2   # Point 0 starts at 2
        first_tstep[1]+=4   # Point 1 starts at 4        
        first_tstep[2]+=3   # Point 2 starts at 3
        
        last_tstep=(time_step_count)*num.ones(n,num.int)
        last_tstep[0]-=1    # Point 0 ends 1 step early
        last_tstep[1]-=2    # Point 1 ends 2 steps early                
        last_tstep[4]-=3    # Point 4 ends 3 steps early        
        
        # Create varying elevation data (positive values for seafloor)
        gauge_depth=20*num.ones(n,num.float)
        for i in range(n):
            gauge_depth[i] += i**2
            
        # Create data to be written to second mux file        
        ha1=num.ones((n,time_step_count),num.float)
        ha1[0]=num.sin(times_ref)
        ha1[1]=2*num.sin(times_ref - 3)
        ha1[2]=5*num.sin(4*times_ref)
        ha1[3]=num.sin(times_ref)
        ha1[4]=num.sin(2*times_ref-0.7)
                
        ua1=num.zeros((n,time_step_count),num.float)
        ua1[0]=3*num.cos(times_ref)        
        ua1[1]=2*num.sin(times_ref-0.7)   
        ua1[2]=num.arange(3*time_step_count,4*time_step_count)
        ua1[4]=2*num.ones(time_step_count)
        
        va1=num.zeros((n,time_step_count),num.float)
        va1[0]=2*num.cos(times_ref-0.87)        
        va1[1]=3*num.ones(time_step_count)
        va1[3]=2*num.sin(times_ref-0.71)        
        
        # Ensure data used to write mux file to be zero when gauges are
        # not recording
        for i in range(n):
             # For each point
             for j in range(0, first_tstep[i]-1) + range(last_tstep[i], time_step_count):
                 # For timesteps before and after recording range
                 ha1[i][j] = ua1[i][j] = va1[i][j] = 0.0 


        #print 'Second station to be written to MUX'
        #print 'ha', ha1[0,:]
        #print 'ua', ua1[0,:]
        #print 'va', va1[0,:]
        
        # Write second mux file to be combined by urs2sts 
        base_nameII, filesII = self.write_mux2(lat_long_points,
                                               time_step_count, time_step,
                                               first_tstep, last_tstep,
                                               depth=gauge_depth,
                                               ha=ha1,
                                               ua=ua1,
                                               va=va1)




        # Read mux file back and verify it's correcness

        ####################################################
        # FIXME (Ole): This is where the test should
        # verify that the MUX files are correct.

        #JJ: It appears as though
        #that certain quantities are not being stored with enough precision
        #inn muxfile or more likely that they are being cast into a
        #lower precision when read in using read_mux2 Time step and q_time
        # are equal but only to approx 1e-7
        ####################################################

        #define information as it should be stored in mus2 files
        points_num=len(lat_long_points)
        depth=gauge_depth
        ha=ha1
        ua=ua1
        va=va1
        
        quantities = ['HA','UA','VA']
        mux_names = [WAVEHEIGHT_MUX2_LABEL,
                     EAST_VELOCITY_MUX2_LABEL,
                     NORTH_VELOCITY_MUX2_LABEL]
        quantities_init = [[],[],[]]
        latlondeps = []
        #irrelevant header information
        ig=ilon=ilat=0
        mcolat=mcolon=centerlat=centerlon=offset=az=baz=id=0.0
        # urs binary is latitude fastest
        for i,point in enumerate(lat_long_points):
            lat = point[0]
            lon = point[1]
            _ , e, n = redfearn(lat, lon)
            if depth is None:
                this_depth = n
            else:
                this_depth = depth[i]
            latlondeps.append([lat, lon, this_depth])

            if ha is None:
                this_ha = e
                quantities_init[0].append(num.ones(time_step_count,num.float)*this_ha) # HA
            else:
                quantities_init[0].append(ha[i])
            if ua is None:
                this_ua = n
                quantities_init[1].append(num.ones(time_step_count,num.float)*this_ua) # UA
            else:
                quantities_init[1].append(ua[i])
            if va is None:
                this_va = e
                quantities_init[2].append(num.ones(time_step_count,num.float)*this_va) #
            else:
                quantities_init[2].append(va[i])

        for i, q in enumerate(quantities):
            #print
            #print i, q
            
            q_time = num.zeros((time_step_count, points_num), num.float)
            quantities_init[i] = ensure_numeric(quantities_init[i])
            for time in range(time_step_count):
                #print i, q, time, quantities_init[i][:,time]
                q_time[time,:] = quantities_init[i][:,time]
                #print i, q, time, q_time[time, :]

            
            filename = base_nameII + mux_names[i]
            f = open(filename, 'rb')
            if self.verbose: print 'Reading' + filename
            assert abs(points_num-unpack('i',f.read(4))[0])<epsilon
            #write mux 2 header
            for latlondep in latlondeps:
                assert abs(latlondep[0]-unpack('f',f.read(4))[0])<epsilon
                assert abs(latlondep[1]-unpack('f',f.read(4))[0])<epsilon
                assert abs(mcolat-unpack('f',f.read(4))[0])<epsilon
                assert abs(mcolon-unpack('f',f.read(4))[0])<epsilon
                assert abs(ig-unpack('i',f.read(4))[0])<epsilon
                assert abs(ilon-unpack('i',f.read(4))[0])<epsilon
                assert abs(ilat-unpack('i',f.read(4))[0])<epsilon
                assert abs(latlondep[2]-unpack('f',f.read(4))[0])<epsilon
                assert abs(centerlat-unpack('f',f.read(4))[0])<epsilon
                assert abs(centerlon-unpack('f',f.read(4))[0])<epsilon
                assert abs(offset-unpack('f',f.read(4))[0])<epsilon
                assert abs(az-unpack('f',f.read(4))[0])<epsilon
                assert abs(baz-unpack('f',f.read(4))[0])<epsilon
                
                x = unpack('f', f.read(4))[0]
                #print time_step
                #print x
                assert abs(time_step-x)<epsilon
                assert abs(time_step_count-unpack('i',f.read(4))[0])<epsilon
                for j in range(4): # identifier
                    assert abs(id-unpack('i',f.read(4))[0])<epsilon 

            #first_tstep=1
            #last_tstep=time_step_count
            for i,latlondep in enumerate(latlondeps):
                assert abs(first_tstep[i]-unpack('i',f.read(4))[0])<epsilon
            for i,latlondep in enumerate(latlondeps):
                assert abs(last_tstep[i]-unpack('i',f.read(4))[0])<epsilon

            # Find when first station starts recording
            min_tstep = min(first_tstep)
            # Find when all stations have stopped recording
            max_tstep = max(last_tstep)

            #for time in  range(time_step_count):
            for time in range(min_tstep-1,max_tstep):
                assert abs(time*time_step-unpack('f',f.read(4))[0])<epsilon
                for point_i in range(points_num):
                    if time+1>=first_tstep[point_i] and time+1<=last_tstep[point_i]:
                        x = unpack('f',f.read(4))[0]
                        #print time, x, q_time[time, point_i]
                        if q == 'VA': x = -x # South is positive in MUX
                        assert abs(q_time[time, point_i]-x)<epsilon

            f.close()
                                               
        # Create ordering file
        permutation = ensure_numeric([4,0,2])

       #  _, ordering_filename = tempfile.mkstemp('')
#         order_fid = open(ordering_filename, 'w')  
#         order_fid.write('index, longitude, latitude\n')
#         for index in permutation:
#             order_fid.write('%d, %f, %f\n' %(index, 
#                                              lat_long_points[index][1], 
#                                              lat_long_points[index][0]))
#         order_fid.close()
        
        # -------------------------------------
        # Now read files back and check values
        weights = ensure_numeric([1.0])

        # For each quantity read the associated list of source mux2 file with 
        # extention associated with that quantity
        file_params=-1*num.ones(3,num.float) # [nsta,dt,nt]
        OFFSET = 5

        for j, file in enumerate(filesII):
            # Read stage, u, v enumerated as j
            #print 'Reading', j, file
            data = read_mux2(1, [file], weights, file_params, permutation, verbose)

            #print 'Data received by Python'
            #print data[1][8]
            number_of_selected_stations = data.shape[0]

            # Index where data ends and parameters begin
            parameters_index = data.shape[1]-OFFSET          
                 
            quantity=num.zeros((number_of_selected_stations, parameters_index), num.float)
            
            
            for i in range(number_of_selected_stations):
        
                #print i, parameters_index
                #print quantity[i][:]
                if j == 0: assert num.allclose(data[i][:parameters_index], ha1[permutation[i], :])
                if j == 1: assert num.allclose(data[i][:parameters_index], ua1[permutation[i], :])
                if j == 2:
                    # FIXME (Ole): This is where the output is wrong on Win32
                    
                    #print
                    #print j, i
                    #print 'Input'
                    #print 'u', ua1[permutation[i], 8]       
                    #print 'v', va1[permutation[i], 8]
                
                    #print 'Output'
                    #print 'v ', data[i][:parameters_index][8]  

                    # South is positive in MUX
                    #print "data[i][:parameters_index]", data[i][:parameters_index]
                    #print "-va1[permutation[i], :]", -va1[permutation[i], :]
                    assert num.allclose(data[i][:parameters_index], -va1[permutation[i], :])
        
        self.delete_mux(filesII)
           
    def test_read_mux_platform_problem3(self):
        
        # This is to test a situation where read_mux returned 
        # wrong values Win32

        
        from urs_ext import read_mux2 
        
        from anuga.config import single_precision as epsilon        
        
        verbose = False
                
        tide = 1.5
        time_step_count = 10
        time_step = 0.02

        '''
        Win results
        time_step = 0.2000001
        This is OK        
        '''
        
        '''
        Win results
        time_step = 0.20000001

        ======================================================================
ERROR: test_read_mux_platform_problem3 (__main__.Test_Data_Manager)
----------------------------------------------------------------------
Traceback (most recent call last):
  File "test_data_manager.py", line 6718, in test_read_mux_platform_problem3
    ha1[0]=num.sin(times_ref)
ValueError: matrices are not aligned for copy

        '''

        '''
        Win results
        time_step = 0.200000001
        FAIL
         assert num.allclose(data[i][:parameters_index],
         -va1[permutation[i], :])
        '''
        times_ref = num.arange(0, time_step_count*time_step, time_step)
        #print "times_ref", times_ref
        
        lat_long_points = [(-21.5,114.5), (-21,114.5), (-21.5,115),
                           (-21.,115.), (-22., 117.)]
        stations = len(lat_long_points)
        
        # Create different timeseries starting and ending at different times 
        first_tstep=num.ones(stations, num.int)
        first_tstep[0]+=2   # Point 0 starts at 2
        first_tstep[1]+=4   # Point 1 starts at 4        
        first_tstep[2]+=3   # Point 2 starts at 3
        
        last_tstep=(time_step_count)*num.ones(stations, num.int)
        last_tstep[0]-=1    # Point 0 ends 1 step early
        last_tstep[1]-=2    # Point 1 ends 2 steps early                
        last_tstep[4]-=3    # Point 4 ends 3 steps early        
        
        # Create varying elevation data (positive values for seafloor)
        gauge_depth=20*num.ones(stations, num.float)
        for i in range(stations):
            gauge_depth[i] += i**2
            
        # Create data to be written to second mux file        
        ha1=num.ones((stations,time_step_count), num.float)
        ha1[0]=num.sin(times_ref)
        ha1[1]=2*num.sin(times_ref - 3)
        ha1[2]=5*num.sin(4*times_ref)
        ha1[3]=num.sin(times_ref)
        ha1[4]=num.sin(2*times_ref-0.7)
                
        ua1=num.zeros((stations,time_step_count),num.float)
        ua1[0]=3*num.cos(times_ref)        
        ua1[1]=2*num.sin(times_ref-0.7)   
        ua1[2]=num.arange(3*time_step_count,4*time_step_count)
        ua1[4]=2*num.ones(time_step_count)
        
        va1=num.zeros((stations,time_step_count),num.float)
        va1[0]=2*num.cos(times_ref-0.87)        
        va1[1]=3*num.ones(time_step_count)
        va1[3]=2*num.sin(times_ref-0.71)        
        #print "va1[0]", va1[0]  # The 8th element is what will go bad.
        # Ensure data used to write mux file to be zero when gauges are
        # not recording
        for i in range(stations):
             # For each point
             for j in range(0, first_tstep[i]-1) + range(last_tstep[i],
                                                         time_step_count):
                 # For timesteps before and after recording range
                 ha1[i][j] = ua1[i][j] = va1[i][j] = 0.0 


        #print 'Second station to be written to MUX'
        #print 'ha', ha1[0,:]
        #print 'ua', ua1[0,:]
        #print 'va', va1[0,:]
        
        # Write second mux file to be combined by urs2sts 
        base_nameII, filesII = self.write_mux2(lat_long_points,
                                               time_step_count, time_step,
                                               first_tstep, last_tstep,
                                               depth=gauge_depth,
                                               ha=ha1,
                                               ua=ua1,
                                               va=va1)
        #print "filesII", filesII




        # Read mux file back and verify it's correcness

        ####################################################
        # FIXME (Ole): This is where the test should
        # verify that the MUX files are correct.

        #JJ: It appears as though
        #that certain quantities are not being stored with enough precision
        #inn muxfile or more likely that they are being cast into a
        #lower precision when read in using read_mux2 Time step and q_time
        # are equal but only to approx 1e-7
        ####################################################

        #define information as it should be stored in mus2 files
        points_num=len(lat_long_points)
        depth=gauge_depth
        ha=ha1
        ua=ua1
        va=va1
        
        quantities = ['HA','UA','VA']
        mux_names = [WAVEHEIGHT_MUX2_LABEL,
                     EAST_VELOCITY_MUX2_LABEL,
                     NORTH_VELOCITY_MUX2_LABEL]
        quantities_init = [[],[],[]]
        latlondeps = []
        #irrelevant header information
        ig=ilon=ilat=0
        mcolat=mcolon=centerlat=centerlon=offset=az=baz=id=0.0
        # urs binary is latitude fastest
        for i,point in enumerate(lat_long_points):
            lat = point[0]
            lon = point[1]
            _ , e, n = redfearn(lat, lon)
            if depth is None:
                this_depth = n
            else:
                this_depth = depth[i]
            latlondeps.append([lat, lon, this_depth])

            if ha is None:
                this_ha = e
                quantities_init[0].append(num.ones(time_step_count,
                                                   num.float)*this_ha) # HA
            else:
                quantities_init[0].append(ha[i])
            if ua is None:
                this_ua = n
                quantities_init[1].append(num.ones(time_step_count,
                                                   num.float)*this_ua) # UA
            else:
                quantities_init[1].append(ua[i])
            if va is None:
                this_va = e
                quantities_init[2].append(num.ones(time_step_count,
                                                   num.float)*this_va) #
            else:
                quantities_init[2].append(va[i])

        for i, q in enumerate(quantities):
            #print
            #print i, q
            
            q_time = num.zeros((time_step_count, points_num), num.float)
            quantities_init[i] = ensure_numeric(quantities_init[i])
            for time in range(time_step_count):
                #print i, q, time, quantities_init[i][:,time]
                q_time[time,:] = quantities_init[i][:,time]
                #print i, q, time, q_time[time, :]

            
            filename = base_nameII + mux_names[i]
            f = open(filename, 'rb')
            if self.verbose: print 'Reading' + filename
            assert abs(points_num-unpack('i',f.read(4))[0])<epsilon
            #write mux 2 header
            for latlondep in latlondeps:
                assert abs(latlondep[0]-unpack('f',f.read(4))[0])<epsilon
                assert abs(latlondep[1]-unpack('f',f.read(4))[0])<epsilon
                assert abs(mcolat-unpack('f',f.read(4))[0])<epsilon
                assert abs(mcolon-unpack('f',f.read(4))[0])<epsilon
                assert abs(ig-unpack('i',f.read(4))[0])<epsilon
                assert abs(ilon-unpack('i',f.read(4))[0])<epsilon
                assert abs(ilat-unpack('i',f.read(4))[0])<epsilon
                assert abs(latlondep[2]-unpack('f',f.read(4))[0])<epsilon
                assert abs(centerlat-unpack('f',f.read(4))[0])<epsilon
                assert abs(centerlon-unpack('f',f.read(4))[0])<epsilon
                assert abs(offset-unpack('f',f.read(4))[0])<epsilon
                assert abs(az-unpack('f',f.read(4))[0])<epsilon
                assert abs(baz-unpack('f',f.read(4))[0])<epsilon
                
                x = unpack('f', f.read(4))[0]
                #print time_step
                #print x
                assert abs(time_step-x)<epsilon
                assert abs(time_step_count-unpack('i',f.read(4))[0])<epsilon
                for j in range(4): # identifier
                    assert abs(id-unpack('i',f.read(4))[0])<epsilon 

            #first_tstep=1
            #last_tstep=time_step_count
            for i,latlondep in enumerate(latlondeps):
                assert abs(first_tstep[i]-unpack('i',f.read(4))[0])<epsilon
            for i,latlondep in enumerate(latlondeps):
                assert abs(last_tstep[i]-unpack('i',f.read(4))[0])<epsilon

            # Find when first station starts recording
            min_tstep = min(first_tstep)
            # Find when all stations have stopped recording
            max_tstep = max(last_tstep)

            #for time in  range(time_step_count):
            for time in range(min_tstep-1,max_tstep):
                assert abs(time*time_step-unpack('f',f.read(4))[0])<epsilon
                for point_i in range(points_num):
                    if time+1>=first_tstep[point_i] and time+1<=last_tstep[point_i]:
                        x = unpack('f',f.read(4))[0]
                        #print time, x, q_time[time, point_i]
                        if q == 'VA': x = -x # South is positive in MUX
                        #print q+" q_time[%d, %d] = %f" %(time, point_i, 
                                                      #q_time[time, point_i])
                        assert abs(q_time[time, point_i]-x)<epsilon

            f.close()
                            
        permutation = ensure_numeric([4,0,2])
                   
        # Create ordering file
#         _, ordering_filename = tempfile.mkstemp('')
#         order_fid = open(ordering_filename, 'w')  
#         order_fid.write('index, longitude, latitude\n')
#         for index in permutation:
#             order_fid.write('%d, %f, %f\n' %(index, 
#                                              lat_long_points[index][1], 
#                                              lat_long_points[index][0]))
#         order_fid.close()
        
        # -------------------------------------
        # Now read files back and check values
        weights = ensure_numeric([1.0])

        # For each quantity read the associated list of source mux2 file with 
        # extention associated with that quantity
        file_params=-1*num.ones(3,num.float) # [nsta,dt,nt]
        OFFSET = 5

        for j, file in enumerate(filesII):
            # Read stage, u, v enumerated as j
            #print 'Reading', j, file
            #print "file", file
            data = read_mux2(1, [file], weights, file_params,
                             permutation, verbose)
            #print str(j) + "data", data

            #print 'Data received by Python'
            #print data[1][8]
            number_of_selected_stations = data.shape[0]
            #print "number_of_selected_stations", number_of_selected_stations
            #print "stations", stations

            # Index where data ends and parameters begin
            parameters_index = data.shape[1]-OFFSET          
                 
            for i in range(number_of_selected_stations):
        
                #print i, parameters_index
                if j == 0:
                    assert num.allclose(data[i][:parameters_index],
                                        ha1[permutation[i], :])
                    
                if j == 1: assert num.allclose(data[i][:parameters_index], ua1[permutation[i], :])
                if j == 2:
                    assert num.allclose(data[i][:parameters_index], -va1[permutation[i], :])
        
        self.delete_mux(filesII)            
        
    def test_urs2sts0(self):
        """
        Test single source
        """
        tide=0
        time_step_count = 3
        time_step = 2
        lat_long_points =[(-21.5,114.5),(-21,114.5),(-21.5,115), (-21.,115.)]
        n=len(lat_long_points)
        first_tstep=num.ones(n,num.int)
        first_tstep[0]+=1
        first_tstep[2]+=1
        last_tstep=(time_step_count)*num.ones(n,num.int)
        last_tstep[0]-=1

        gauge_depth=20*num.ones(n,num.float)
        ha=2*num.ones((n,time_step_count),num.float)
        ha[0]=num.arange(0,time_step_count)
        ha[1]=num.arange(time_step_count,2*time_step_count)
        ha[2]=num.arange(2*time_step_count,3*time_step_count)
        ha[3]=num.arange(3*time_step_count,4*time_step_count)
        ua=5*num.ones((n,time_step_count),num.float)
        va=-10*num.ones((n,time_step_count),num.float)

        base_name, files = self.write_mux2(lat_long_points,
                                      time_step_count, time_step,
                                      first_tstep, last_tstep,
                                      depth=gauge_depth,
                                      ha=ha,
                                      ua=ua,
                                      va=va)

        urs2sts(base_name,
                basename_out=base_name, 
                mean_stage=tide,verbose=False)

        # now I want to check the sts file ...
        sts_file = base_name + '.sts'

        #Let's interigate the sww file
        # Note, the sww info is not gridded.  It is point data.
        fid = NetCDFFile(sts_file)

        # Make x and y absolute
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]

        geo_reference = Geo_reference(NetCDFObject=fid)
        points = geo_reference.get_absolute(map(None, x, y))
        points = ensure_numeric(points)

        x = points[:,0]
        y = points[:,1]

        #Check that first coordinate is correctly represented       
        #Work out the UTM coordinates for first point
        for i in range(4):
            zone, e, n = redfearn(lat_long_points[i][0], lat_long_points[i][1]) 
            assert num.allclose([x[i],y[i]], [e,n])

        #Check the time vector
        times = fid.variables['time'][:]

        times_actual = []
        for i in range(time_step_count):
            times_actual.append(time_step * i)

        assert num.allclose(ensure_numeric(times),
                            ensure_numeric(times_actual))

        #Check first value
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        ymomentum = fid.variables['ymomentum'][:]
        elevation = fid.variables['elevation'][:]

        # Set original data used to write mux file to be zero when gauges are
        #not recdoring
        ha[0][0]=0.0
        ha[0][time_step_count-1]=0.0;
        ha[2][0]=0.0;
        ua[0][0]=0.0
        ua[0][time_step_count-1]=0.0;
        ua[2][0]=0.0;
        va[0][0]=0.0
        va[0][time_step_count-1]=0.0;
        va[2][0]=0.0;

        assert num.allclose(num.transpose(ha),stage)  #Meters

        #Check the momentums - ua
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_ua *(stage+depth)

        depth=num.zeros((len(lat_long_points),time_step_count),num.float)
        for i in range(len(lat_long_points)):
            depth[i]=gauge_depth[i]+tide+ha[i]
        assert num.allclose(num.transpose(ua*depth),xmomentum) 

        #Check the momentums - va
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_va *(stage+depth)

        assert num.allclose(num.transpose(va*depth),ymomentum)

        # check the elevation values.
        # -ve since urs measures depth, sww meshers height,
        assert num.allclose(-elevation, gauge_depth)  #Meters

        fid.close()
        self.delete_mux(files)
        os.remove(sts_file)

    def test_urs2sts_nonstandard_meridian(self):
        """
        Test single source using the meridian from zone 50 as a nonstandard meridian
        """
        tide=0
        time_step_count = 3
        time_step = 2
        lat_long_points =[(-21.,114.5),(-21.,113.5),(-21.,114.), (-21.,115.)]
        n=len(lat_long_points)
        first_tstep=num.ones(n,num.int)
        first_tstep[0]+=1
        first_tstep[2]+=1
        last_tstep=(time_step_count)*num.ones(n,num.int)
        last_tstep[0]-=1

        gauge_depth=20*num.ones(n,num.float)
        ha=2*num.ones((n,time_step_count),num.float)
        ha[0]=num.arange(0,time_step_count)
        ha[1]=num.arange(time_step_count,2*time_step_count)
        ha[2]=num.arange(2*time_step_count,3*time_step_count)
        ha[3]=num.arange(3*time_step_count,4*time_step_count)
        ua=5*num.ones((n,time_step_count),num.float)
        va=-10*num.ones((n,time_step_count),num.float)

        base_name, files = self.write_mux2(lat_long_points,
                                           time_step_count, time_step,
                                           first_tstep, last_tstep,
                                           depth=gauge_depth,
                                           ha=ha,
                                           ua=ua,
                                           va=va)

        urs2sts(base_name,
                basename_out=base_name, 
                central_meridian=123,
                mean_stage=tide,
                verbose=False)

        # now I want to check the sts file ...
        sts_file = base_name + '.sts'

        #Let's interigate the sww file
        # Note, the sww info is not gridded.  It is point data.
        fid = NetCDFFile(sts_file)

        # Make x and y absolute
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]

        geo_reference = Geo_reference(NetCDFObject=fid)
        points = geo_reference.get_absolute(map(None, x, y))
        points = ensure_numeric(points)

        x = points[:,0]
        y = points[:,1]

        # Check that all coordinate are correctly represented       
        # Using the non standard projection (50) 
        for i in range(4):
            zone, e, n = redfearn(lat_long_points[i][0],
                                  lat_long_points[i][1],
                                  central_meridian=123)
            assert num.allclose([x[i],y[i]], [e,n])
            assert zone==-1
        
        self.delete_mux(files)

            
    def test_urs2sts_nonstandard_projection_reverse(self):
        """
        Test that a point not in the specified zone can occur first
        """
        tide=0
        time_step_count = 3
        time_step = 2
        lat_long_points =[(-21.,113.5),(-21.,114.5),(-21.,114.), (-21.,115.)]
        n=len(lat_long_points)
        first_tstep=num.ones(n,num.int)
        first_tstep[0]+=1
        first_tstep[2]+=1
        last_tstep=(time_step_count)*num.ones(n,num.int)
        last_tstep[0]-=1

        gauge_depth=20*num.ones(n,num.float)
        ha=2*num.ones((n,time_step_count),num.float)
        ha[0]=num.arange(0,time_step_count)
        ha[1]=num.arange(time_step_count,2*time_step_count)
        ha[2]=num.arange(2*time_step_count,3*time_step_count)
        ha[3]=num.arange(3*time_step_count,4*time_step_count)
        ua=5*num.ones((n,time_step_count),num.float)
        va=-10*num.ones((n,time_step_count),num.float)

        base_name, files = self.write_mux2(lat_long_points,
                                      time_step_count, time_step,
                                      first_tstep, last_tstep,
                                      depth=gauge_depth,
                                      ha=ha,
                                      ua=ua,
                                      va=va)

        urs2sts(base_name,
                basename_out=base_name, 
                zone=50,
                mean_stage=tide,verbose=False)

        # now I want to check the sts file ...
        sts_file = base_name + '.sts'

        #Let's interigate the sww file
        # Note, the sww info is not gridded.  It is point data.
        fid = NetCDFFile(sts_file)

        # Make x and y absolute
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]

        geo_reference = Geo_reference(NetCDFObject=fid)
        points = geo_reference.get_absolute(map(None, x, y))
        points = ensure_numeric(points)

        x = points[:,0]
        y = points[:,1]

        # Check that all coordinate are correctly represented       
        # Using the non standard projection (50) 
        for i in range(4):
            zone, e, n = redfearn(lat_long_points[i][0], lat_long_points[i][1],
                                  zone=50) 
            assert num.allclose([x[i],y[i]], [e,n])
            assert zone==geo_reference.zone
        
        self.delete_mux(files)

            
    def test_urs2stsII(self):
        """
        Test multiple sources
        """
        tide=0
        time_step_count = 3
        time_step = 2
        lat_long_points =[(-21.5,114.5),(-21,114.5),(-21.5,115), (-21.,115.)]
        n=len(lat_long_points)
        first_tstep=num.ones(n,num.int)
        first_tstep[0]+=1
        first_tstep[2]+=1
        last_tstep=(time_step_count)*num.ones(n,num.int)
        last_tstep[0]-=1

        gauge_depth=20*num.ones(n,num.float)
        ha=2*num.ones((n,time_step_count),num.float)
        ha[0]=num.arange(0,time_step_count)
        ha[1]=num.arange(time_step_count,2*time_step_count)
        ha[2]=num.arange(2*time_step_count,3*time_step_count)
        ha[3]=num.arange(3*time_step_count,4*time_step_count)
        ua=5*num.ones((n,time_step_count),num.float)
        va=-10*num.ones((n,time_step_count),num.float)

        # Create two identical mux files to be combined by urs2sts
        base_nameI, filesI = self.write_mux2(lat_long_points,
                                             time_step_count, time_step,
                                             first_tstep, last_tstep,
                                             depth=gauge_depth,
                                             ha=ha,
                                             ua=ua,
                                             va=va)

        base_nameII, filesII = self.write_mux2(lat_long_points,
                                               time_step_count, time_step,
                                               first_tstep, last_tstep,
                                               depth=gauge_depth,
                                               ha=ha,
                                               ua=ua,
                                               va=va)

        # Call urs2sts with multiple mux files
        urs2sts([base_nameI, base_nameII], 
                basename_out=base_nameI, 
                weights=[1.0, 1.0],
                mean_stage=tide,
                verbose=False)

        # now I want to check the sts file ...
        sts_file = base_nameI + '.sts'

        #Let's interrogate the sts file
        # Note, the sts info is not gridded.  It is point data.
        fid = NetCDFFile(sts_file)

        # Make x and y absolute
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]

        geo_reference = Geo_reference(NetCDFObject=fid)
        points = geo_reference.get_absolute(map(None, x, y))
        points = ensure_numeric(points)

        x = points[:,0]
        y = points[:,1]

        #Check that first coordinate is correctly represented       
        #Work out the UTM coordinates for first point
        zone, e, n = redfearn(lat_long_points[0][0], lat_long_points[0][1]) 
        assert num.allclose([x[0],y[0]], [e,n])

        #Check the time vector
        times = fid.variables['time'][:]

        times_actual = []
        for i in range(time_step_count):
            times_actual.append(time_step * i)

        assert num.allclose(ensure_numeric(times),
                            ensure_numeric(times_actual))

        #Check first value
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        ymomentum = fid.variables['ymomentum'][:]
        elevation = fid.variables['elevation'][:]

        # Set original data used to write mux file to be zero when gauges are
        # not recdoring
        
        ha[0][0]=0.0
        ha[0][time_step_count-1]=0.0
        ha[2][0]=0.0
        ua[0][0]=0.0
        ua[0][time_step_count-1]=0.0
        ua[2][0]=0.0
        va[0][0]=0.0
        va[0][time_step_count-1]=0.0
        va[2][0]=0.0;

        # The stage stored in the .sts file should be the sum of the stage
        # in the two mux2 files because both have weights = 1. In this case
        # the mux2 files are the same so stage == 2.0 * ha
        #print 2.0*num.transpose(ha) - stage 
        assert num.allclose(2.0*num.transpose(ha), stage)  #Meters

        #Check the momentums - ua
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_ua *(stage+depth)

        depth=num.zeros((len(lat_long_points),time_step_count),num.float)
        for i in range(len(lat_long_points)):
            depth[i]=gauge_depth[i]+tide+2.0*ha[i]
            #2.0*ha necessary because using two files with weights=1 are used

        # The xmomentum stored in the .sts file should be the sum of the ua
        # in the two mux2 files multiplied by the depth.
        assert num.allclose(2.0*num.transpose(ua*depth), xmomentum) 

        #Check the momentums - va
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_va *(stage+depth)

        # The ymomentum stored in the .sts file should be the sum of the va
        # in the two mux2 files multiplied by the depth.
        assert num.allclose(2.0*num.transpose(va*depth), ymomentum)

        # check the elevation values.
        # -ve since urs measures depth, sww meshers height,
        assert num.allclose(-elevation, gauge_depth)  #Meters

        fid.close()
        self.delete_mux(filesI)
        self.delete_mux(filesII)
        os.remove(sts_file)

    def test_urs2sts_individual_sources(self):   
        """Test that individual sources compare to actual urs output
           Test that the first recording time is the smallest
           over waveheight, easting and northing velocity
        """
        
        # Get path where this test is run
        path = get_pathname_from_package('anuga.shallow_water')        
        
        testdir = os.path.join(path, 'urs_test_data')
        ordering_filename=os.path.join(testdir, 'thinned_bound_order_test.txt')
        
        sources = ['1-z.grd','2-z.grd','3-z.grd']
        
        # Start times by source and station taken manually from urs header files
        time_start_z = num.array([[10.0,11.5,13,14.5,17.7],
                                  [9.8,11.2,12.7,14.2,17.4],
                                  [9.5,10.9,12.4,13.9,17.1]])

        time_start_e = time_start_n = time_start_z

        # time step in urs output
        delta_t = 0.1
        
        # Make sts file for each source
        for k, source_filename in enumerate(sources):
            source_number = k + 1 # Source numbering starts at 1
            
            urs_filenames = os.path.join(testdir, source_filename)
            weights = [1.]
            sts_name_out = 'test'
            
            urs2sts(urs_filenames,
                    basename_out=sts_name_out,
                    ordering_filename=ordering_filename,
                    weights=weights,
                    mean_stage=0.0,
                    verbose=False)

            # Read in sts file for this source file
            fid = NetCDFFile(sts_name_out+'.sts', netcdf_mode_r) # Open existing file for read
            x = fid.variables['x'][:]+fid.xllcorner    # x-coordinates of vertices
            y = fid.variables['y'][:]+fid.yllcorner    # y-coordinates of vertices
            elevation = fid.variables['elevation'][:]
            time=fid.variables['time'][:]+fid.starttime

            # Get quantity data from sts file
            quantity_names=['stage','xmomentum','ymomentum']
            quantities = {}
            for i, name in enumerate(quantity_names):
                quantities[name] = fid.variables[name][:]
            
            # Make sure start time from sts file is the minimum starttime 
            # across all stations (for this source)
            #print k, time_start_z[k,:]
            starttime = min(time_start_z[k, :])
            sts_starttime = fid.starttime[0]
            msg = 'sts starttime for source %d was %f. Should have been %f'\
                %(source_number, sts_starttime, starttime)
            assert num.allclose(sts_starttime, starttime), msg             

            # For each station, compare urs2sts output to known urs output
            for j in range(len(x)):
                index_start_urs = 0
                index_end_urs = 0
                index_start = 0
                index_end = 0
                count = 0

                # read in urs test data for stage, e and n velocity
                urs_file_name_z = 'z_'+str(source_number)+'_'+str(j)+'.csv'
                dict = load_csv_as_array(os.path.join(testdir, urs_file_name_z))
                urs_stage = dict['urs_stage']
                urs_file_name_e = 'e_'+str(source_number)+'_'+str(j)+'.csv'
                dict = load_csv_as_array(os.path.join(testdir, urs_file_name_e))
                urs_e = dict['urs_e']
                urs_file_name_n = 'n_'+str(source_number)+'_'+str(j)+'.csv'
                dict = load_csv_as_array(os.path.join(testdir, urs_file_name_n))
                urs_n = dict['urs_n']

                # find start and end time for stage             
                for i in range(len(urs_stage)):
                    if urs_stage[i] == 0.0:
                        index_start_urs_z = i+1
                    if int(urs_stage[i]) == 99 and count <> 1:
                        count +=1
                        index_end_urs_z = i

                if count == 0: index_end_urs_z = len(urs_stage)

                start_times_z = time_start_z[source_number-1]

                # start times for easting velocities should be the same as stage
                start_times_e = time_start_e[source_number-1]
                index_start_urs_e = index_start_urs_z
                index_end_urs_e = index_end_urs_z

                # start times for northing velocities should be the same as stage
                start_times_n = time_start_n[source_number-1]
                index_start_urs_n = index_start_urs_z
                index_end_urs_n = index_end_urs_z
                
                # Check that actual start time matches header information for stage
                msg = 'stage start time from urs file is not the same as the '
                msg += 'header file for source %i and station %i' %(source_number,j)
                assert num.allclose(index_start_urs_z,start_times_z[j]/delta_t), msg

                msg = 'e velocity start time from urs file is not the same as the '
                msg += 'header file for source %i and station %i' %(source_number,j)
                assert num.allclose(index_start_urs_e,start_times_e[j]/delta_t), msg

                msg = 'n velocity start time from urs file is not the same as the '
                msg += 'header file for source %i and station %i' %(source_number,j)
                assert num.allclose(index_start_urs_n,start_times_n[j]/delta_t), msg
                
                # get index for start and end time for sts quantities
                index_start_stage = 0
                index_end_stage = 0
                count = 0
                sts_stage = quantities['stage'][:,j]
                for i in range(len(sts_stage)):
                    if sts_stage[i] <> 0.0 and count <> 1:
                        count += 1
                        index_start_stage = i
                    if int(sts_stage[i]) == 99 and count <> 1:
                        count += 1
                        index_end_stage = i

                index_end_stage = index_start_stage + len(urs_stage[index_start_urs_z:index_end_urs_z])

                sts_xmom = quantities['xmomentum'][:,j]
                index_start_x = index_start_stage
                index_end_x = index_start_x + len(urs_e[index_start_urs_e:index_end_urs_e])

                sts_ymom = quantities['ymomentum'][:,j]
                index_start_y = index_start_stage
                index_end_y = index_start_y + len(urs_n[index_start_urs_n:index_end_urs_n])

                # check that urs stage and sts stage are the same
                msg = 'urs stage is not equal to sts stage for for source %i and station %i' %(source_number,j)
                assert num.allclose(urs_stage[index_start_urs_z:index_end_urs_z],
                                    sts_stage[index_start_stage:index_end_stage], 
                                    rtol=1.0e-6, atol=1.0e-5 ), msg                                
                                
                # check that urs e velocity and sts xmomentum are the same
                msg = 'urs e velocity is not equivalent to sts x momentum for for source %i and station %i' %(source_number,j)
                assert num.allclose(urs_e[index_start_urs_e:index_end_urs_e]*(urs_stage[index_start_urs_e:index_end_urs_e]-elevation[j]),
                                sts_xmom[index_start_x:index_end_x], 
                                rtol=1.0e-5, atol=1.0e-4 ), msg
                
                # check that urs n velocity and sts ymomentum are the same
                #print 'urs n velocity', urs_n[index_start_urs_n:index_end_urs_n]*(urs_stage[index_start_urs_n:index_end_urs_n]-elevation[j])
                #print 'sts momentum', sts_ymom[index_start_y:index_end_y]                                                             
                msg = 'urs n velocity is not equivalent to sts y momentum for source %i and station %i' %(source_number,j)
                assert num.allclose(urs_n[index_start_urs_n:index_end_urs_n]*(urs_stage[index_start_urs_n:index_end_urs_n]-elevation[j]),
                                -sts_ymom[index_start_y:index_end_y], 
                                rtol=1.0e-5, atol=1.0e-4 ), msg
                                                
                                
            fid.close()
            
        os.remove(sts_name_out+'.sts')

    def test_urs2sts_combined_sources(self):   
        """Test that combined sources compare to actual urs output
           Test that the first recording time is the smallest
           over waveheight, easting and northing velocity
        """

        # combined
        time_start_z = num.array([9.5,10.9,12.4,13.9,17.1])
        time_start_e = time_start_n = time_start_z
         
        # make sts file for combined sources
        weights = [1., 2., 3.]
        
        path = get_pathname_from_package('anuga.shallow_water')        
                
        testdir = os.path.join(path, 'urs_test_data')        
        ordering_filename=os.path.join(testdir, 'thinned_bound_order_test.txt')

        urs_filenames = [os.path.join(testdir,'1-z.grd'),
                         os.path.join(testdir,'2-z.grd'),
                         os.path.join(testdir,'3-z.grd')]
        sts_name_out = 'test'
        
        urs2sts(urs_filenames,
                basename_out=sts_name_out,
                ordering_filename=ordering_filename,
                weights=weights,
                mean_stage=0.0,
                verbose=False)
        
        # read in sts file for combined source
        fid = NetCDFFile(sts_name_out+'.sts', netcdf_mode_r)    # Open existing file for read
        x = fid.variables['x'][:]+fid.xllcorner   # x-coordinates of vertices
        y = fid.variables['y'][:]+fid.yllcorner   # y-coordinates of vertices
        elevation = fid.variables['elevation'][:]
        time=fid.variables['time'][:]+fid.starttime
        
        
        # Check that stored permutation is as per default
        permutation = range(len(x))
        stored_permutation = fid.variables['permutation'][:]
        msg = 'Permutation was not stored correctly. I got '
        msg += str(stored_permutation)
        assert num.allclose(stored_permutation, permutation), msg        

        # Get quantity data from sts file
        quantity_names=['stage','xmomentum','ymomentum']
        quantities = {}
        for i, name in enumerate(quantity_names):
            quantities[name] = fid.variables[name][:]

        # For each station, compare urs2sts output to known urs output   
        delta_t = 0.1
        
        # Make sure start time from sts file is the minimum starttime 
        # across all stations (for this source)
        starttime = min(time_start_z[:])
        sts_starttime = fid.starttime[0]
        msg = 'sts starttime was %f. Should have been %f'\
            %(sts_starttime, starttime)
        assert num.allclose(sts_starttime, starttime), msg
    
        #stations = [1,2,3]
        #for j in stations: 
        for j in range(len(x)):
            index_start_urs_z = 0
            index_end_urs_z = 0
            index_start_urs_e = 0
            index_end_urs_e = 0
            index_start_urs_n = 0
            index_end_urs_n = 0
            count = 0

            # read in urs test data for stage, e and n velocity
            urs_file_name_z = 'z_combined_'+str(j)+'.csv'
            dict = load_csv_as_array(os.path.join(testdir, urs_file_name_z))
            urs_stage = dict['urs_stage']
            urs_file_name_e = 'e_combined_'+str(j)+'.csv'
            dict = load_csv_as_array(os.path.join(testdir, urs_file_name_e))
            urs_e = dict['urs_e']
            urs_file_name_n = 'n_combined_'+str(j)+'.csv'
            dict = load_csv_as_array(os.path.join(testdir, urs_file_name_n))
            urs_n = dict['urs_n']

            # find start and end time for stage         
            for i in range(len(urs_stage)):
                if urs_stage[i] == 0.0:
                    index_start_urs_z = i+1
                if int(urs_stage[i]) == 99 and count <> 1:
                    count +=1
                    index_end_urs_z = i

            if count == 0: index_end_urs_z = len(urs_stage)

            start_times_z = time_start_z[j]

            start_times_e = time_start_e[j]
            index_start_urs_e = index_start_urs_z

            start_times_n = time_start_n[j]
            index_start_urs_n = index_start_urs_z
               
            # Check that actual start time matches header information for stage
            msg = 'stage start time from urs file is not the same as the '
            msg += 'header file at station %i' %(j)
            assert num.allclose(index_start_urs_z,start_times_z/delta_t), msg

            msg = 'e velocity start time from urs file is not the same as the '
            msg += 'header file at station %i' %(j)
            assert num.allclose(index_start_urs_e,start_times_e/delta_t), msg

            msg = 'n velocity start time from urs file is not the same as the '
            msg += 'header file at station %i' %(j)
            assert num.allclose(index_start_urs_n,start_times_n/delta_t), msg
                
            # get index for start and end time for sts quantities
            index_start_stage = 0
            index_end_stage = 0
            index_start_x = 0
            index_end_x = 0
            index_start_y = 0
            index_end_y = 0
            count = 0
            count1 = 0
            sts_stage = quantities['stage'][:,j]
            for i in range(len(sts_stage)):
                if sts_stage[i] <> 0.0 and count <> 1:
                    count += 1
                    index_start_stage = i
                if int(urs_stage[i]) == 99 and count <> 1:
                    count +=1
                    index_end_stage = i
                
            index_end_stage = index_start_stage + len(urs_stage[index_start_urs_z:index_end_urs_z])

            index_start_x = index_start_stage
            index_end_x = index_start_x + len(urs_stage[index_start_urs_e:index_end_urs_e])
            sts_xmom = quantities['ymomentum'][:,j]

            index_start_y = index_start_stage
            index_end_y = index_start_y + len(urs_stage[index_start_urs_n:index_end_urs_n])
            sts_ymom = quantities['ymomentum'][:,j]

            # check that urs stage and sts stage are the same
            msg = 'urs stage is not equal to sts stage for station %i' %j
            #print 'urs stage', urs_stage[index_start_urs_z:index_end_urs_z]
            #print 'sts stage', sts_stage[index_start_stage:index_end_stage]
            #print 'diff', max(urs_stage[index_start_urs_z:index_end_urs_z]-sts_stage[index_start_stage:index_end_stage])
            #print 'index', index_start_stage, index_end_stage, len(sts_stage)
            assert num.allclose(urs_stage[index_start_urs_z:index_end_urs_z],
                            sts_stage[index_start_stage:index_end_stage], 
                                rtol=1.0e-5, atol=1.0e-4 ), msg                                
                                
            # check that urs e velocity and sts xmomentum are the same          
            msg = 'urs e velocity is not equivalent to sts xmomentum for station %i' %j
            assert num.allclose(urs_e[index_start_urs_e:index_end_urs_e]*(urs_stage[index_start_urs_e:index_end_urs_e]-elevation[j]),
                            sts_xmom[index_start_x:index_end_x], 
                            rtol=1.0e-5, atol=1.0e-4 ), msg
                
            # check that urs n velocity and sts ymomentum are the same                            
            msg = 'urs n velocity is not equivalent to sts ymomentum for station %i' %j
            assert num.allclose(urs_n[index_start_urs_n:index_end_urs_n]*(urs_stage[index_start_urs_n:index_end_urs_n]-elevation[j]),
                            sts_ymom[index_start_y:index_end_y], 
                            rtol=1.0e-5, atol=1.0e-4 ), msg

        fid.close()
        
        os.remove(sts_name_out+'.sts')
        
        
        
    def test_urs2sts_ordering(self):
        """Test multiple sources with ordering file
        """
        
        tide = 0.35
        time_step_count = 6 # I made this a little longer (Ole)
        time_step = 2
        lat_long_points =[(-21.5,114.5),(-21,114.5),(-21.5,115), (-21.,115.)]
        n=len(lat_long_points)
        first_tstep=num.ones(n,num.int)
        first_tstep[0]+=1
        first_tstep[2]+=1
        last_tstep=(time_step_count)*num.ones(n,num.int)
        last_tstep[0]-=1

        gauge_depth=20*num.ones(n,num.float)
        ha=2*num.ones((n,time_step_count),num.float)
        ha[0]=num.arange(0,time_step_count)
        ha[1]=num.arange(time_step_count,2*time_step_count)
        ha[2]=num.arange(2*time_step_count,3*time_step_count)
        ha[3]=num.arange(3*time_step_count,4*time_step_count)
        ua=5*num.ones((n,time_step_count),num.float)
        va=-10*num.ones((n,time_step_count),num.float)

        # Create two identical mux files to be combined by urs2sts
        base_nameI, filesI = self.write_mux2(lat_long_points,
                                             time_step_count, time_step,
                                             first_tstep, last_tstep,
                                             depth=gauge_depth,
                                             ha=ha,
                                             ua=ua,
                                             va=va)

        base_nameII, filesII = self.write_mux2(lat_long_points,
                                               time_step_count, time_step,
                                               first_tstep, last_tstep,
                                               depth=gauge_depth,
                                               ha=ha,
                                               ua=ua,
                                               va=va)

                                               
        # Create ordering file
        permutation = [3,0,2]

        _, ordering_filename = tempfile.mkstemp('')
        order_fid = open(ordering_filename, 'w')  
        order_fid.write('index, longitude, latitude\n')
        for index in permutation:
            order_fid.write('%d, %f, %f\n' %(index, 
                                             lat_long_points[index][1], 
                                             lat_long_points[index][0]))
        order_fid.close()
        
            

                                               
        # Call urs2sts with multiple mux files
        urs2sts([base_nameI, base_nameII], 
                basename_out=base_nameI, 
                ordering_filename=ordering_filename,
                weights=[1.0, 1.0],
                mean_stage=tide,
                verbose=False)

        # now I want to check the sts file ...
        sts_file = base_nameI + '.sts'

        #Let's interrogate the sts file
        # Note, the sts info is not gridded.  It is point data.
        fid = NetCDFFile(sts_file)
        
        # Check that original indices have been stored
        stored_permutation = fid.variables['permutation'][:]
        msg = 'Permutation was not stored correctly. I got '
        msg += str(stored_permutation)
        assert num.allclose(stored_permutation, permutation), msg
        

        # Make x and y absolute
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]

        geo_reference = Geo_reference(NetCDFObject=fid)
        points = geo_reference.get_absolute(map(None, x, y))
        points = ensure_numeric(points)

        x = points[:,0]
        y = points[:,1]

        #print
        #print x
        #print y
        for i, index in enumerate(permutation):
            # Check that STS points are stored in the correct order
            
            # Work out the UTM coordinates sts point i
            zone, e, n = redfearn(lat_long_points[index][0], 
                                  lat_long_points[index][1])             

            #print i, [x[i],y[i]], [e,n]
            assert num.allclose([x[i],y[i]], [e,n])
            
                        
        # Check the time vector
        times = fid.variables['time'][:]

        times_actual = []
        for i in range(time_step_count):
            times_actual.append(time_step * i)

        assert num.allclose(ensure_numeric(times),
                            ensure_numeric(times_actual))
                        

        # Check sts values
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        ymomentum = fid.variables['ymomentum'][:]
        elevation = fid.variables['elevation'][:]

        # Set original data used to write mux file to be zero when gauges are
        # not recdoring
        
        ha[0][0]=0.0
        ha[0][time_step_count-1]=0.0
        ha[2][0]=0.0
        ua[0][0]=0.0
        ua[0][time_step_count-1]=0.0
        ua[2][0]=0.0
        va[0][0]=0.0
        va[0][time_step_count-1]=0.0
        va[2][0]=0.0;

        
        # The stage stored in the .sts file should be the sum of the stage
        # in the two mux2 files because both have weights = 1. In this case
        # the mux2 files are the same so stage == 2.0 * ha
        #print 2.0*num.transpose(ha) - stage 
        
        ha_permutation = num.take(ha, permutation, axis=0) 
        ua_permutation = num.take(ua, permutation, axis=0)
        va_permutation = num.take(va, permutation, axis=0)
        gauge_depth_permutation = num.take(gauge_depth, permutation, axis=0)

        
        assert num.allclose(2.0*num.transpose(ha_permutation)+tide, stage)  # Meters

        #Check the momentums - ua
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_ua *(stage+depth)

        depth=num.zeros((len(lat_long_points),time_step_count),num.float)
        for i in range(len(lat_long_points)):
            depth[i]=gauge_depth[i]+tide+2.0*ha[i]
            #2.0*ha necessary because using two files with weights=1 are used
            
        depth_permutation = num.take(depth, permutation, axis=0)                     
        

        # The xmomentum stored in the .sts file should be the sum of the ua
        # in the two mux2 files multiplied by the depth.
        assert num.allclose(2.0*num.transpose(ua_permutation*depth_permutation), xmomentum) 

        #Check the momentums - va
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_va *(stage+depth)

        # The ymomentum stored in the .sts file should be the sum of the va
        # in the two mux2 files multiplied by the depth.
        assert num.allclose(2.0*num.transpose(va_permutation*depth_permutation), ymomentum)

        # check the elevation values.
        # -ve since urs measures depth, sww meshers height,
        assert num.allclose(-gauge_depth_permutation, elevation)  #Meters

        fid.close()
        self.delete_mux(filesI)
        self.delete_mux(filesII)
        os.remove(sts_file)
        

        
        
        
    def Xtest_urs2sts_ordering_exception(self):
        """Test that inconsistent lats and lons in ordering file are caught.
        """
        
        tide=0
        time_step_count = 3
        time_step = 2
        lat_long_points =[(-21.5,114.5),(-21,114.5),(-21.5,115), (-21.,115.)]
        n=len(lat_long_points)
        first_tstep=num.ones(n,num.int)
        first_tstep[0]+=1
        first_tstep[2]+=1
        last_tstep=(time_step_count)*num.ones(n,num.int)
        last_tstep[0]-=1

        gauge_depth=20*num.ones(n,num.float)
        ha=2*num.ones((n,time_step_count),num.float)
        ha[0]=num.arange(0,time_step_count)
        ha[1]=num.arange(time_step_count,2*time_step_count)
        ha[2]=num.arange(2*time_step_count,3*time_step_count)
        ha[3]=num.arange(3*time_step_count,4*time_step_count)
        ua=5*num.ones((n,time_step_count),num.float)
        va=-10*num.ones((n,time_step_count),num.float)

        # Create two identical mux files to be combined by urs2sts
        base_nameI, filesI = self.write_mux2(lat_long_points,
                                             time_step_count, time_step,
                                             first_tstep, last_tstep,
                                             depth=gauge_depth,
                                             ha=ha,
                                             ua=ua,
                                             va=va)

        base_nameII, filesII = self.write_mux2(lat_long_points,
                                               time_step_count, time_step,
                                               first_tstep, last_tstep,
                                               depth=gauge_depth,
                                               ha=ha,
                                               ua=ua,
                                               va=va)

                                               
        # Create ordering file
        permutation = [3,0,2]

        # Do it wrongly and check that exception is being raised
        _, ordering_filename = tempfile.mkstemp('')
        order_fid = open(ordering_filename, 'w')  
        order_fid.write('index, longitude, latitude\n')
        for index in permutation:
            order_fid.write('%d, %f, %f\n' %(index, 
                                             lat_long_points[index][0], 
                                             lat_long_points[index][1]))
        order_fid.close()
        
        try:
            urs2sts([base_nameI, base_nameII], 
                    basename_out=base_nameI, 
                    ordering_filename=ordering_filename,
                    weights=[1.0, 1.0],
                    mean_stage=tide,
                    verbose=False)  
            os.remove(ordering_filename)            
        except:
            pass
        else:
            msg = 'Should have caught wrong lat longs'
            raise Exception, msg

        
        self.delete_mux(filesI)
        self.delete_mux(filesII)

        

        
    def test_urs2sts_ordering_different_sources(self):
        """Test multiple sources with ordering file, different source files and weights.
           This test also has more variable values than the previous ones
        """
        
        tide = 1.5
        time_step_count = 10
        time_step = 0.2
        
        times_ref = num.arange(0, time_step_count*time_step, time_step)
        #print 'time vector', times_ref
        
        lat_long_points = [(-21.5,114.5), (-21,114.5), (-21.5,115), (-21.,115.), (-22., 117.)]
        n = len(lat_long_points)
        
        # Create non-trivial weights
        #weights = [0.8, 1.5] # OK
        #weights = [0.8, 10.5] # Fail (up to allclose tolerance)
        #weights = [10.5, 10.5] # OK
        #weights = [0.0, 10.5] # OK
        #weights = [0.8, 0.] # OK                
        #weights = [8, 0.1] # OK                        
        #weights = [0.8, 10.0] # OK                                
        #weights = [0.8, 10.6] # OK           
        weights = [3.8, 7.6] # OK                   
        #weights = [0.5, 0.5] # OK                           
        #weights = [2., 2.] # OK                            
        #weights = [0.0, 0.5] # OK                                          
        #weights = [1.0, 1.0] # OK                                                  
        
        
        # Create different timeseries starting and ending at different times 
        first_tstep=num.ones(n,num.int)
        first_tstep[0]+=2   # Point 0 starts at 2
        first_tstep[1]+=4   # Point 1 starts at 4        
        first_tstep[2]+=3   # Point 2 starts at 3
        
        last_tstep=(time_step_count)*num.ones(n,num.int)
        last_tstep[0]-=1    # Point 0 ends 1 step early
        last_tstep[1]-=2    # Point 1 ends 2 steps early                
        last_tstep[4]-=3    # Point 4 ends 3 steps early        
        
        #print
        #print 'time_step_count', time_step_count
        #print 'time_step', time_step
        #print 'first_tstep', first_tstep
        #print 'last_tstep', last_tstep                
        
        
        # Create varying elevation data (positive values for seafloor)
        gauge_depth=20*num.ones(n,num.float)
        for i in range(n):
            gauge_depth[i] += i**2
            
        #print 'gauge_depth', gauge_depth
        
        # Create data to be written to first mux file        
        ha0=2*num.ones((n,time_step_count),num.float)
        ha0[0]=num.arange(0,time_step_count)
        ha0[1]=num.arange(time_step_count,2*time_step_count)
        ha0[2]=num.arange(2*time_step_count,3*time_step_count)
        ha0[3]=num.arange(3*time_step_count,4*time_step_count)
        ua0=5*num.ones((n,time_step_count),num.float)
        va0=-10*num.ones((n,time_step_count),num.float)

        # Ensure data used to write mux file to be zero when gauges are
        # not recording
        for i in range(n):
             # For each point
             
             for j in range(0, first_tstep[i]-1) + range(last_tstep[i], time_step_count):
                 # For timesteps before and after recording range
                 ha0[i][j] = ua0[i][j] = va0[i][j] = 0.0                                  


                 
        #print 
        #print 'using varying start and end time'
        #print 'ha0', ha0[4]
        #print 'ua0', ua0
        #print 'va0', va0        
        
        # Write first mux file to be combined by urs2sts
        base_nameI, filesI = self.write_mux2(lat_long_points,
                                             time_step_count, time_step,
                                             first_tstep, last_tstep,
                                             depth=gauge_depth,
                                             ha=ha0,
                                             ua=ua0,
                                             va=va0)

                                             
                                             
        # Create data to be written to second mux file        
        ha1=num.ones((n,time_step_count),num.float)
        ha1[0]=num.sin(times_ref)
        ha1[1]=2*num.sin(times_ref - 3)
        ha1[2]=5*num.sin(4*times_ref)
        ha1[3]=num.sin(times_ref)
        ha1[4]=num.sin(2*times_ref-0.7)
                
        ua1=num.zeros((n,time_step_count),num.float)
        ua1[0]=3*num.cos(times_ref)        
        ua1[1]=2*num.sin(times_ref-0.7)   
        ua1[2]=num.arange(3*time_step_count,4*time_step_count)
        ua1[4]=2*num.ones(time_step_count)
        
        va1=num.zeros((n,time_step_count),num.float)
        va1[0]=2*num.cos(times_ref-0.87)        
        va1[1]=3*num.ones(time_step_count, num.int)       #array default#
        va1[3]=2*num.sin(times_ref-0.71)        
        
        
        # Ensure data used to write mux file to be zero when gauges are
        # not recording
        for i in range(n):
             # For each point
             
             for j in range(0, first_tstep[i]-1) + range(last_tstep[i], time_step_count):
                 # For timesteps before and after recording range
                 ha1[i][j] = ua1[i][j] = va1[i][j] = 0.0                                  


        #print 
        #print 'using varying start and end time'
        #print 'ha1', ha1[4]
        #print 'ua1', ua1
        #print 'va1', va1        
                                             
                                             
        # Write second mux file to be combined by urs2sts                                             
        base_nameII, filesII = self.write_mux2(lat_long_points,
                                               time_step_count, time_step,
                                               first_tstep, last_tstep,
                                               depth=gauge_depth,
                                               ha=ha1,
                                               ua=ua1,
                                               va=va1)

                                               
        # Create ordering file
        permutation = [4,0,2]

        _, ordering_filename = tempfile.mkstemp('')
        order_fid = open(ordering_filename, 'w')  
        order_fid.write('index, longitude, latitude\n')
        for index in permutation:
            order_fid.write('%d, %f, %f\n' %(index, 
                                             lat_long_points[index][1], 
                                             lat_long_points[index][0]))
        order_fid.close()
        
            

        #------------------------------------------------------------
        # Now read the mux files one by one without weights and test
        
        # Call urs2sts with mux file #0
        urs2sts([base_nameI], 
                basename_out=base_nameI, 
                ordering_filename=ordering_filename,
                mean_stage=tide,
                verbose=False)

        # Now read the sts file and check that values have been stored correctly.
        sts_file = base_nameI + '.sts'
        fid = NetCDFFile(sts_file)
        

        # Check that original indices have been stored
        stored_permutation = fid.variables['permutation'][:]
        msg = 'Permutation was not stored correctly. I got '
        msg += str(stored_permutation)
        assert num.allclose(stored_permutation, permutation), msg
        

        
        
        # Make x and y absolute
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]

        geo_reference = Geo_reference(NetCDFObject=fid)
        points = geo_reference.get_absolute(map(None, x, y))
        points = ensure_numeric(points)

        x = points[:,0]
        y = points[:,1]

        for i, index in enumerate(permutation):
            # Check that STS points are stored in the correct order
            
            # Work out the UTM coordinates sts point i
            zone, e, n = redfearn(lat_long_points[index][0], 
                                  lat_long_points[index][1])             

            #print i, [x[i],y[i]], [e,n]
            assert num.allclose([x[i],y[i]], [e,n])
            
                        
        # Check the time vector
        times = fid.variables['time'][:]
        assert num.allclose(ensure_numeric(times),
                            ensure_numeric(times_ref))
                        

        # Check sts values for mux #0
        stage0 = fid.variables['stage'][:].copy()
        xmomentum0 = fid.variables['xmomentum'][:].copy()
        ymomentum0 = fid.variables['ymomentum'][:].copy()
        elevation0 = fid.variables['elevation'][:].copy()

        
        #print 'stage', stage0
        #print 'xmomentum', xmomentum0
        #print 'ymomentum', ymomentum0        
        #print 'elevation', elevation0
        
        # The quantities stored in the .sts file should be the weighted sum of the 
        # quantities written to the mux2 files subject to the permutation vector.
        
        ha_ref = num.take(ha0, permutation, axis=0)
        ua_ref = num.take(ua0, permutation, axis=0)        
        va_ref = num.take(va0, permutation, axis=0)                

        gauge_depth_ref = num.take(gauge_depth, permutation, axis=0)                      
        
        assert num.allclose(num.transpose(ha_ref)+tide, stage0)  # Meters
        
        
        
        #Check the momentums - ua
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_ua *(stage+depth)

        depth_ref = num.zeros((len(permutation), time_step_count), num.float)
        for i in range(len(permutation)):
            depth_ref[i]=gauge_depth_ref[i]+tide+ha_ref[i]


        # The xmomentum stored in the .sts file should be the sum of the ua
        # in the two mux2 files multiplied by the depth.
        assert num.allclose(num.transpose(ua_ref*depth_ref), xmomentum0) 

        #Check the momentums - va
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_va *(stage+depth)

        # The ymomentum stored in the .sts file should be the sum of the va
        # in the two mux2 files multiplied by the depth.
        
        
        #print transpose(va_ref*depth_ref)
        #print ymomentum
        assert num.allclose(num.transpose(va_ref*depth_ref), ymomentum0)        

        # check the elevation values.
        # -ve since urs measures depth, sww meshers height,
        assert num.allclose(-gauge_depth_ref, elevation0) 

        fid.close()
        os.remove(sts_file)
        
        

        
        # Call urs2sts with mux file #1
        urs2sts([base_nameII], 
                basename_out=base_nameI, 
                ordering_filename=ordering_filename,
                mean_stage=tide,
                verbose=False)

        # Now read the sts file and check that values have been stored correctly.
        sts_file = base_nameI + '.sts'
        fid = NetCDFFile(sts_file)
        
        
        # Check that original indices have been stored
        stored_permutation = fid.variables['permutation'][:]
        msg = 'Permutation was not stored correctly. I got '
        msg += str(stored_permutation)
        assert num.allclose(stored_permutation, permutation), msg
        
        # Make x and y absolute
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]

        geo_reference = Geo_reference(NetCDFObject=fid)
        points = geo_reference.get_absolute(map(None, x, y))
        points = ensure_numeric(points)

        x = points[:,0]
        y = points[:,1]

        for i, index in enumerate(permutation):
            # Check that STS points are stored in the correct order
            
            # Work out the UTM coordinates sts point i
            zone, e, n = redfearn(lat_long_points[index][0], 
                                  lat_long_points[index][1])             

            #print i, [x[i],y[i]], [e,n]
            assert num.allclose([x[i],y[i]], [e,n])
            
                        
        # Check the time vector
        times = fid.variables['time'][:]
        assert num.allclose(ensure_numeric(times),
                            ensure_numeric(times_ref))
                        

        # Check sts values for mux #1 
        stage1 = fid.variables['stage'][:].copy()
        xmomentum1 = fid.variables['xmomentum'][:].copy()
        ymomentum1 = fid.variables['ymomentum'][:].copy()
        elevation1 = fid.variables['elevation'][:].copy()

        
        #print 'stage', stage1
        #print 'xmomentum', xmomentum1
        #print 'ymomentum', ymomentum1       
        #print 'elevation', elevation1
        
        # The quantities stored in the .sts file should be the weighted sum of the 
        # quantities written to the mux2 files subject to the permutation vector.
        
        ha_ref = num.take(ha1, permutation, axis=0)
        ua_ref = num.take(ua1, permutation, axis=0)        
        va_ref = num.take(va1, permutation, axis=0)                

        gauge_depth_ref = num.take(gauge_depth, permutation, axis=0)                         


        #print 
        #print stage1
        #print transpose(ha_ref)+tide - stage1
        

        assert num.allclose(num.transpose(ha_ref)+tide, stage1)  # Meters
        #import sys; sys.exit()

        #Check the momentums - ua
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_ua *(stage+depth)

        depth_ref = num.zeros((len(permutation), time_step_count), num.float)
        for i in range(len(permutation)):
            depth_ref[i]=gauge_depth_ref[i]+tide+ha_ref[i]


        # The xmomentum stored in the .sts file should be the sum of the ua
        # in the two mux2 files multiplied by the depth.
        assert num.allclose(num.transpose(ua_ref*depth_ref), xmomentum1) 

        #Check the momentums - va
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_va *(stage+depth)

        # The ymomentum stored in the .sts file should be the sum of the va
        # in the two mux2 files multiplied by the depth.
        
        
        #print transpose(va_ref*depth_ref)
        #print ymomentum
        assert num.allclose(num.transpose(va_ref*depth_ref), ymomentum1)        

        # check the elevation values.
        # -ve since urs measures depth, sww meshers height,
        assert num.allclose(-gauge_depth_ref, elevation1) 

        fid.close()
        os.remove(sts_file)
        
        #----------------------
        # Then read the mux files together and test
        
                                               
        # Call urs2sts with multiple mux files
        urs2sts([base_nameI, base_nameII], 
                basename_out=base_nameI, 
                ordering_filename=ordering_filename,
                weights=weights,
                mean_stage=tide,
                verbose=False)

        # Now read the sts file and check that values have been stored correctly.
        sts_file = base_nameI + '.sts'
        fid = NetCDFFile(sts_file)

        # Make x and y absolute
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]

        geo_reference = Geo_reference(NetCDFObject=fid)
        points = geo_reference.get_absolute(map(None, x, y))
        points = ensure_numeric(points)

        x = points[:,0]
        y = points[:,1]

        for i, index in enumerate(permutation):
            # Check that STS points are stored in the correct order
            
            # Work out the UTM coordinates sts point i
            zone, e, n = redfearn(lat_long_points[index][0], 
                                  lat_long_points[index][1])             

            #print i, [x[i],y[i]], [e,n]
            assert num.allclose([x[i],y[i]], [e,n])
            
                        
        # Check the time vector
        times = fid.variables['time'][:]
        assert num.allclose(ensure_numeric(times),
                            ensure_numeric(times_ref))
                        

        # Check sts values
        stage = fid.variables['stage'][:].copy()
        xmomentum = fid.variables['xmomentum'][:].copy()
        ymomentum = fid.variables['ymomentum'][:].copy()
        elevation = fid.variables['elevation'][:].copy()

        
        #print 'stage', stage
        #print 'elevation', elevation
        
        # The quantities stored in the .sts file should be the weighted sum of the 
        # quantities written to the mux2 files subject to the permutation vector.
        
        ha_ref = (weights[0]*num.take(ha0, permutation, axis=0)
                  + weights[1]*num.take(ha1, permutation, axis=0))
        ua_ref = (weights[0]*num.take(ua0, permutation, axis=0)
                  + weights[1]*num.take(ua1, permutation, axis=0))
        va_ref = (weights[0]*num.take(va0, permutation, axis=0)
                  + weights[1]*num.take(va1, permutation, axis=0))

        gauge_depth_ref = num.take(gauge_depth, permutation, axis=0)                         


        #print 
        #print stage
        #print transpose(ha_ref)+tide - stage

        assert num.allclose(num.transpose(ha_ref)+tide, stage)  # Meters

        #Check the momentums - ua
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_ua *(stage+depth)

        depth_ref = num.zeros((len(permutation), time_step_count), num.float)
        for i in range(len(permutation)):
            depth_ref[i]=gauge_depth_ref[i]+tide+ha_ref[i]

            
        

        # The xmomentum stored in the .sts file should be the sum of the ua
        # in the two mux2 files multiplied by the depth.
        assert num.allclose(num.transpose(ua_ref*depth_ref), xmomentum) 

        #Check the momentums - va
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_va *(stage+depth)

        # The ymomentum stored in the .sts file should be the sum of the va
        # in the two mux2 files multiplied by the depth.
        
        
        #print transpose(va_ref*depth_ref)
        #print ymomentum

        assert num.allclose(num.transpose(va_ref*depth_ref), ymomentum)

        # check the elevation values.
        # -ve since urs measures depth, sww meshers height,
        assert num.allclose(-gauge_depth_ref, elevation)  #Meters

        fid.close()
        self.delete_mux(filesI)
        self.delete_mux(filesII)
        os.remove(sts_file)
        
        #---------------
        # "Manually" add the timeseries up with weights and test
        # Tide is discounted from individual results and added back in       
        #

        stage_man = weights[0]*(stage0-tide) + weights[1]*(stage1-tide) + tide
        assert num.allclose(stage_man, stage)
                
        
    def test_file_boundary_stsI(self):
        """test_file_boundary_stsI(self):
        """
        
        # FIXME (Ole): These tests should really move to 
        # test_generic_boundaries.py
        
        from anuga.shallow_water import Domain
        from anuga.shallow_water import Reflective_boundary
        from anuga.shallow_water import Dirichlet_boundary
        from anuga.shallow_water import File_boundary
        from anuga.pmesh.mesh_interface import create_mesh_from_regions

        bounding_polygon=[[6.0,97.0],[6.01,97.0],[6.02,97.0],[6.02,97.02],[6.00,97.02]]
        tide = 0.37
        time_step_count = 5
        time_step = 2
        lat_long_points =bounding_polygon[0:3]
        n=len(lat_long_points)
        first_tstep=num.ones(n,num.int)
        last_tstep=(time_step_count)*num.ones(n,num.int)

        h = 20        
        w = 2
        u = 10
        v = -10
        gauge_depth=h*num.ones(n,num.float)
        ha=w*num.ones((n,time_step_count),num.float)
        ua=u*num.ones((n,time_step_count),num.float)
        va=v*num.ones((n,time_step_count),num.float)
        base_name, files = self.write_mux2(lat_long_points,
                                           time_step_count, time_step,
                                           first_tstep, last_tstep,
                                           depth=gauge_depth,
                                           ha=ha,
                                           ua=ua,
                                           va=va)

        sts_file=base_name
        urs2sts(base_name,
                sts_file,
                mean_stage=tide,
                verbose=False)
        self.delete_mux(files)

        #print 'start create mesh from regions'
        for i in range(len(bounding_polygon)):
            zone,bounding_polygon[i][0],bounding_polygon[i][1]=redfearn(bounding_polygon[i][0],bounding_polygon[i][1])
        extent_res=1000000
        meshname = 'urs_test_mesh' + '.tsh'
        interior_regions=None
        boundary_tags={'ocean': [0,1], 'otherocean': [2,3,4]}
        create_mesh_from_regions(bounding_polygon,
                                 boundary_tags=boundary_tags,
                                 maximum_triangle_area=extent_res,
                                 filename=meshname,
                                 interior_regions=interior_regions,
                                 verbose=False)
        
        domain_fbound = Domain(meshname)
        domain_fbound.set_quantity('stage', tide)
        Bf = File_boundary(sts_file+'.sts', 
                           domain_fbound, 
                           boundary_polygon=bounding_polygon)
        Br = Reflective_boundary(domain_fbound)
        Bd = Dirichlet_boundary([w+tide, u*(w+h+tide), v*(w+h+tide)])        

        domain_fbound.set_boundary({'ocean': Bf,'otherocean': Br})
        
        # Check boundary object evaluates as it should
        for i, ((vol_id, edge_id), B) in enumerate(domain_fbound.boundary_objects):
            if B is Bf:
            
                qf = B.evaluate(vol_id, edge_id)  # File boundary
                qd = Bd.evaluate(vol_id, edge_id) # Dirichlet boundary

                assert num.allclose(qf, qd) 
                
        
        # Evolve
        finaltime=time_step*(time_step_count-1)
        yieldstep=time_step
        temp_fbound=num.zeros(int(finaltime/yieldstep)+1,num.float)

        for i, t in enumerate(domain_fbound.evolve(yieldstep=yieldstep,
                                                   finaltime=finaltime, 
                                                   skip_initial_step=False)):
                                                   
            D = domain_fbound
            temp_fbound[i]=D.quantities['stage'].centroid_values[2]

            # Check that file boundary object has populated 
            # boundary array correctly  
            # FIXME (Ole): Do this for the other tests too!
            for j, val in enumerate(D.get_quantity('stage').boundary_values):
            
                (vol_id, edge_id), B = D.boundary_objects[j]
                if isinstance(B, File_boundary):
                    #print j, val
                    assert num.allclose(val, w + tide)


        
        domain_drchlt = Domain(meshname)
        domain_drchlt.set_quantity('stage', tide)
        Br = Reflective_boundary(domain_drchlt)

        domain_drchlt.set_boundary({'ocean': Bd,'otherocean': Br})
        temp_drchlt=num.zeros(int(finaltime/yieldstep)+1,num.float)

        for i, t in enumerate(domain_drchlt.evolve(yieldstep=yieldstep,
                                                   finaltime=finaltime, 
                                                   skip_initial_step=False)):
            temp_drchlt[i]=domain_drchlt.quantities['stage'].centroid_values[2]

        #print domain_fbound.quantities['stage'].vertex_values
        #print domain_drchlt.quantities['stage'].vertex_values
                    
        assert num.allclose(temp_fbound,temp_drchlt)
        
        assert num.allclose(domain_fbound.quantities['stage'].vertex_values,
                            domain_drchlt.quantities['stage'].vertex_values)
                        
        assert num.allclose(domain_fbound.quantities['xmomentum'].vertex_values,
                            domain_drchlt.quantities['xmomentum'].vertex_values) 
                        
        assert num.allclose(domain_fbound.quantities['ymomentum'].vertex_values,
                            domain_drchlt.quantities['ymomentum'].vertex_values)
        
        
        os.remove(sts_file+'.sts')
        os.remove(meshname)
                
        
    def test_file_boundary_stsI_beyond_model_time(self):
        """test_file_boundary_stsI(self):
        
        Test that file_boundary can work when model time
        exceeds available data.
        This is optional and requires the user to specify a default 
        boundary object.
        """
        
        # Don't do warnings in unit test
        import warnings
        warnings.simplefilter('ignore')
        
        from anuga.shallow_water import Domain
        from anuga.shallow_water import Reflective_boundary
        from anuga.shallow_water import Dirichlet_boundary
        from anuga.shallow_water import File_boundary
        from anuga.pmesh.mesh_interface import create_mesh_from_regions

        bounding_polygon=[[6.0,97.0],[6.01,97.0],[6.02,97.0],
                          [6.02,97.02],[6.00,97.02]]
        tide = 0.37
        time_step_count = 5
        time_step = 2
        lat_long_points = bounding_polygon[0:3]
        n=len(lat_long_points)
        first_tstep=num.ones(n,num.int)
        last_tstep=(time_step_count)*num.ones(n,num.int)

        h = 20        
        w = 2
        u = 10
        v = -10
        gauge_depth=h*num.ones(n,num.float)
        ha=w*num.ones((n,time_step_count),num.float)
        ua=u*num.ones((n,time_step_count),num.float)
        va=v*num.ones((n,time_step_count),num.float)
        base_name, files = self.write_mux2(lat_long_points,
                                           time_step_count, time_step,
                                           first_tstep, last_tstep,
                                           depth=gauge_depth,
                                           ha=ha,
                                           ua=ua,
                                           va=va)

        sts_file=base_name
        urs2sts(base_name,
                sts_file,
                mean_stage=tide,
                verbose=False)
        self.delete_mux(files)

        #print 'start create mesh from regions'
        for i in range(len(bounding_polygon)):
            zone,\
            bounding_polygon[i][0],\
            bounding_polygon[i][1]=redfearn(bounding_polygon[i][0],
                                            bounding_polygon[i][1])
                                            
        extent_res=1000000
        meshname = 'urs_test_mesh' + '.tsh'
        interior_regions=None
        boundary_tags={'ocean': [0,1], 'otherocean': [2,3,4]}
        create_mesh_from_regions(bounding_polygon,
                                 boundary_tags=boundary_tags,
                                 maximum_triangle_area=extent_res,
                                 filename=meshname,
                                 interior_regions=interior_regions,
                                 verbose=False)
        
        domain_fbound = Domain(meshname)
        domain_fbound.set_quantity('stage', tide)
        
        Br = Reflective_boundary(domain_fbound)
        Bd = Dirichlet_boundary([w+tide, u*(w+h+tide), v*(w+h+tide)])        
        Bf = File_boundary(sts_file+'.sts', 
                           domain_fbound, 
                           boundary_polygon=bounding_polygon,
                           default_boundary=Bd) # Condition to be used 
                                                # if model time exceeds 
                                                # available data

        domain_fbound.set_boundary({'ocean': Bf,'otherocean': Br})
        
        # Check boundary object evaluates as it should
        for i, ((vol_id, edge_id), B) in enumerate(domain_fbound.boundary_objects):
            if B is Bf:
            
                qf = B.evaluate(vol_id, edge_id)  # File boundary
                qd = Bd.evaluate(vol_id, edge_id) # Dirichlet boundary

                assert num.allclose(qf, qd) 
                
        
        # Evolve
        data_finaltime = time_step*(time_step_count-1)
        finaltime = data_finaltime + 10 # Let model time exceed available data
        yieldstep = time_step
        temp_fbound=num.zeros(int(finaltime/yieldstep)+1, num.float)

        for i, t in enumerate(domain_fbound.evolve(yieldstep=yieldstep,
                                                   finaltime=finaltime, 
                                                   skip_initial_step=False)):
                                                   
            D = domain_fbound
            temp_fbound[i]=D.quantities['stage'].centroid_values[2]

            # Check that file boundary object has populated 
            # boundary array correctly  
            # FIXME (Ole): Do this for the other tests too!
            for j, val in enumerate(D.get_quantity('stage').boundary_values):
            
                (vol_id, edge_id), B = D.boundary_objects[j]
                if isinstance(B, File_boundary):
                    #print j, val
                    assert num.allclose(val, w + tide)


        domain_drchlt = Domain(meshname)
        domain_drchlt.set_quantity('stage', tide)
        Br = Reflective_boundary(domain_drchlt)

        domain_drchlt.set_boundary({'ocean': Bd,'otherocean': Br})
        temp_drchlt=num.zeros(int(finaltime/yieldstep)+1,num.float)

        for i, t in enumerate(domain_drchlt.evolve(yieldstep=yieldstep,
                                                   finaltime=finaltime, 
                                                   skip_initial_step=False)):
            temp_drchlt[i]=domain_drchlt.quantities['stage'].centroid_values[2]

        #print domain_fbound.quantities['stage'].vertex_values
        #print domain_drchlt.quantities['stage'].vertex_values
                    
        assert num.allclose(temp_fbound,temp_drchlt)
        
        assert num.allclose(domain_fbound.quantities['stage'].vertex_values,
                            domain_drchlt.quantities['stage'].vertex_values)
                        
        assert num.allclose(domain_fbound.quantities['xmomentum'].vertex_values,
                            domain_drchlt.quantities['xmomentum'].vertex_values) 
                        
        assert num.allclose(domain_fbound.quantities['ymomentum'].vertex_values,
                            domain_drchlt.quantities['ymomentum'].vertex_values) 
        
        os.remove(sts_file+'.sts')
        os.remove(meshname)
                
        
    def test_field_boundary_stsI_beyond_model_time(self):
        """test_field_boundary(self):
        
        Test that field_boundary can work when model time
        exceeds available data whilst adjusting mean_stage.
        
        """
        
        # Don't do warnings in unit test
        import warnings
        warnings.simplefilter('ignore')
        
        from anuga.shallow_water import Domain
        from anuga.shallow_water import Reflective_boundary
        from anuga.shallow_water import Dirichlet_boundary
        from anuga.shallow_water import File_boundary
        from anuga.pmesh.mesh_interface import create_mesh_from_regions

        bounding_polygon=[[6.0,97.0],[6.01,97.0],[6.02,97.0],
                          [6.02,97.02],[6.00,97.02]]
        tide = 0.37
        time_step_count = 5
        time_step = 2
        lat_long_points = bounding_polygon[0:3]
        n=len(lat_long_points)
        first_tstep=num.ones(n,num.int)
        last_tstep=(time_step_count)*num.ones(n,num.int)

        h = 20        
        w = 2
        u = 10
        v = -10
        gauge_depth=h*num.ones(n,num.float)
        ha=w*num.ones((n,time_step_count),num.float)
        ua=u*num.ones((n,time_step_count),num.float)
        va=v*num.ones((n,time_step_count),num.float)
        base_name, files = self.write_mux2(lat_long_points,
                                           time_step_count, time_step,
                                           first_tstep, last_tstep,
                                           depth=gauge_depth,
                                           ha=ha,
                                           ua=ua,
                                           va=va)

        sts_file=base_name
        urs2sts(base_name,
                sts_file,
                mean_stage=0.0, # Deliberately let Field_boundary do the adjustment
                verbose=False)
        self.delete_mux(files)

        #print 'start create mesh from regions'
        for i in range(len(bounding_polygon)):
            zone,\
            bounding_polygon[i][0],\
            bounding_polygon[i][1]=redfearn(bounding_polygon[i][0],
                                            bounding_polygon[i][1])
                                            
        extent_res=1000000
        meshname = 'urs_test_mesh' + '.tsh'
        interior_regions=None
        boundary_tags={'ocean': [0,1], 'otherocean': [2,3,4]}
        create_mesh_from_regions(bounding_polygon,
                                 boundary_tags=boundary_tags,
                                 maximum_triangle_area=extent_res,
                                 filename=meshname,
                                 interior_regions=interior_regions,
                                 verbose=False)
        
        domain_fbound = Domain(meshname)
        domain_fbound.set_quantity('stage', tide)
        
        Br = Reflective_boundary(domain_fbound)
        Bd = Dirichlet_boundary([w+tide, u*(w+h), v*(w+h)])
        Bdefault = Dirichlet_boundary([w, u*(w+h), v*(w+h)])        
                
        Bf = Field_boundary(sts_file+'.sts', 
                           domain_fbound, 
                           mean_stage=tide, # Field boundary to adjust for tide
                           boundary_polygon=bounding_polygon,
                           default_boundary=Bdefault) # Condition to be used 
                                                      # if model time exceeds 
                                                      # available data

        domain_fbound.set_boundary({'ocean': Bf,'otherocean': Br})
        
        # Check boundary object evaluates as it should
        for i, ((vol_id, edge_id), B) in enumerate(domain_fbound.boundary_objects):
            if B is Bf:
            
                qf = B.evaluate(vol_id, edge_id)  # Field boundary
                qd = Bd.evaluate(vol_id, edge_id) # Dirichlet boundary
                
                msg = 'Got %s, should have been %s' %(qf, qd)
                assert num.allclose(qf, qd), msg 
                
        # Evolve
        data_finaltime = time_step*(time_step_count-1)
        finaltime = data_finaltime + 10 # Let model time exceed available data
        yieldstep = time_step
        temp_fbound=num.zeros(int(finaltime/yieldstep)+1, num.float)

        for i, t in enumerate(domain_fbound.evolve(yieldstep=yieldstep,
                                                   finaltime=finaltime, 
                                                   skip_initial_step=False)):
                                                   
            D = domain_fbound
            temp_fbound[i]=D.quantities['stage'].centroid_values[2]

            # Check that file boundary object has populated 
            # boundary array correctly  
            # FIXME (Ole): Do this for the other tests too!
            for j, val in enumerate(D.get_quantity('stage').boundary_values):
            
                (vol_id, edge_id), B = D.boundary_objects[j]
                if isinstance(B, Field_boundary):
                    msg = 'Got %f should have been %f' %(val, w+tide)
                    assert num.allclose(val, w + tide), msg


    def test_file_boundary_stsII(self):
        """test_file_boundary_stsII(self):
        
         mux2 file has points not included in boundary
         mux2 gauges are not stored with the same order as they are 
         found in bounding_polygon. This does not matter as long as bounding
         polygon passed to file_function contains the mux2 points needed (in
         the correct order).
         """
         
        from anuga.shallow_water import Domain
        from anuga.shallow_water import Reflective_boundary
        from anuga.shallow_water import Dirichlet_boundary
        from anuga.shallow_water import File_boundary
        from anuga.pmesh.mesh_interface import create_mesh_from_regions

        bounding_polygon=[[6.01,97.0],[6.02,97.0],[6.02,97.02],[6.00,97.02],[6.0,97.0]]
        tide = -2.20 
        time_step_count = 20
        time_step = 2
        lat_long_points=bounding_polygon[0:2]
        lat_long_points.insert(0,bounding_polygon[len(bounding_polygon)-1])
        lat_long_points.insert(0,[6.0,97.01])
        n=len(lat_long_points)
        first_tstep=num.ones(n,num.int)
        last_tstep=(time_step_count)*num.ones(n,num.int)
        gauge_depth=20*num.ones(n,num.float)
        ha=2*num.ones((n,time_step_count),num.float)
        ua=10*num.ones((n,time_step_count),num.float)
        va=-10*num.ones((n,time_step_count),num.float)
        base_name, files = self.write_mux2(lat_long_points,
                                           time_step_count,
                                           time_step,
                                           first_tstep,
                                           last_tstep,
                                           depth=gauge_depth,
                                           ha=ha,
                                           ua=ua,
                                           va=va)

        sts_file=base_name
        urs2sts(base_name,sts_file,mean_stage=tide,verbose=False)
        self.delete_mux(files)

        #print 'start create mesh from regions'
        for i in range(len(bounding_polygon)):
            zone,\
            bounding_polygon[i][0],\
            bounding_polygon[i][1]=redfearn(bounding_polygon[i][0],
                                            bounding_polygon[i][1])
            
        extent_res=1000000
        meshname = 'urs_test_mesh' + '.tsh'
        interior_regions=None
        boundary_tags={'ocean': [0,1], 'otherocean': [2,3,4]}
        # have to change boundary tags from last example because now bounding
        # polygon starts in different place.
        create_mesh_from_regions(bounding_polygon,boundary_tags=boundary_tags,
                         maximum_triangle_area=extent_res,filename=meshname,
                         interior_regions=interior_regions,verbose=False)
        
        domain_fbound = Domain(meshname)
        domain_fbound.set_quantity('stage', tide)
        Bf = File_boundary(sts_file+'.sts',
                           domain_fbound,
                           boundary_polygon=bounding_polygon)
        Br = Reflective_boundary(domain_fbound)

        domain_fbound.set_boundary({'ocean': Bf,'otherocean': Br})
        finaltime=time_step*(time_step_count-1)
        yieldstep=time_step
        temp_fbound=num.zeros(int(finaltime/yieldstep)+1,num.float)
        
        for i, t in enumerate(domain_fbound.evolve(yieldstep=yieldstep,
                                                   finaltime=finaltime, 
                                                   skip_initial_step = False)):
            temp_fbound[i]=domain_fbound.quantities['stage'].centroid_values[2]
        
        domain_drchlt = Domain(meshname)
        domain_drchlt.set_quantity('stage', tide)
        Br = Reflective_boundary(domain_drchlt)
        Bd = Dirichlet_boundary([2.0+tide,220+10*tide,-220-10*tide])
        domain_drchlt.set_boundary({'ocean': Bd,'otherocean': Br})
        temp_drchlt=num.zeros(int(finaltime/yieldstep)+1,num.float)

        for i, t in enumerate(domain_drchlt.evolve(yieldstep=yieldstep,
                                                   finaltime=finaltime, 
                                                   skip_initial_step = False)):
            temp_drchlt[i]=domain_drchlt.quantities['stage'].centroid_values[2]


        assert num.allclose(temp_fbound,temp_drchlt)            
            
        #print domain_fbound.quantities['stage'].vertex_values
        #print domain_drchlt.quantities['stage'].vertex_values
                    
            
        assert num.allclose(domain_fbound.quantities['stage'].vertex_values,
                            domain_drchlt.quantities['stage'].vertex_values)
                        
        assert num.allclose(domain_fbound.quantities['xmomentum'].vertex_values,
                            domain_drchlt.quantities['xmomentum'].vertex_values)
                        
        assert num.allclose(domain_fbound.quantities['ymomentum'].vertex_values,
                            domain_drchlt.quantities['ymomentum'].vertex_values)
            
            

        os.remove(sts_file+'.sts')
        os.remove(meshname)

        
        
    def test_file_boundary_stsIII_ordering(self):
        """test_file_boundary_stsIII_ordering(self):
        Read correct points from ordering file and apply sts to boundary
        """
        from anuga.shallow_water import Domain
        from anuga.shallow_water import Reflective_boundary
        from anuga.shallow_water import Dirichlet_boundary
        from anuga.shallow_water import File_boundary
        from anuga.pmesh.mesh_interface import create_mesh_from_regions

        lat_long_points=[[6.01,97.0],[6.02,97.0],[6.05,96.9],[6.0,97.0]]
        bounding_polygon=[[6.0,97.0],[6.01,97.0],[6.02,97.0],
                          [6.02,97.02],[6.00,97.02]]
        tide = 3.0
        time_step_count = 50
        time_step = 2
        n=len(lat_long_points)
        first_tstep=num.ones(n,num.int)
        last_tstep=(time_step_count)*num.ones(n,num.int)
        gauge_depth=20*num.ones(n,num.float)
        ha=2*num.ones((n,time_step_count),num.float)
        ua=10*num.ones((n,time_step_count),num.float)
        va=-10*num.ones((n,time_step_count),num.float)
        base_name, files = self.write_mux2(lat_long_points,
                                           time_step_count,
                                           time_step,
                                           first_tstep,
                                           last_tstep,
                                           depth=gauge_depth,
                                           ha=ha,
                                           ua=ua,
                                           va=va)

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
        urs2sts(base_name,
                basename_out=sts_file,
                ordering_filename=order_file,
                mean_stage=tide,
                verbose=False)
        self.delete_mux(files)

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

        plot=False
        if plot:
            from pylab import plot,show,axis
            boundary_polygon=ensure_numeric(boundary_polygon)
            bounding_polygon_utm=ensure_numeric(bounding_polygon_utm)
            #lat_long_points=ensure_numeric(lat_long_points)
            #plot(lat_long_points[:,0],lat_long_points[:,1],'o')
            plot(boundary_polygon[:,0], boundary_polygon[:,1],'d')
            plot(bounding_polygon_utm[:,0],bounding_polygon_utm[:,1],'o')
            show()

        assert num.allclose(bounding_polygon_utm,boundary_polygon)


        extent_res=1000000
        meshname = 'urs_test_mesh' + '.tsh'
        interior_regions=None
        boundary_tags={'ocean': [0,1], 'otherocean': [2,3,4]}
        
        # have to change boundary tags from last example because now bounding
        # polygon starts in different place.
        create_mesh_from_regions(boundary_polygon,boundary_tags=boundary_tags,
                         maximum_triangle_area=extent_res,filename=meshname,
                         interior_regions=interior_regions,verbose=False)
        
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
                                                   skip_initial_step = False)):
            temp_fbound[i]=domain_fbound.quantities['stage'].centroid_values[2]
    
        
        domain_drchlt = Domain(meshname)
        domain_drchlt.set_quantity('stage', tide)
        Br = Reflective_boundary(domain_drchlt)
        Bd = Dirichlet_boundary([2.0+tide,220+10*tide,-220-10*tide])
        domain_drchlt.set_boundary({'ocean': Bd,'otherocean': Br})
        temp_drchlt=num.zeros(int(finaltime/yieldstep)+1,num.float)
        
        for i, t in enumerate(domain_drchlt.evolve(yieldstep=yieldstep,
                                                   finaltime=finaltime, 
                                                   skip_initial_step = False)):
            temp_drchlt[i]=domain_drchlt.quantities['stage'].centroid_values[2]

        
        #print domain_fbound.quantities['stage'].vertex_values
        #print domain_drchlt.quantities['stage'].vertex_values
                    
        assert num.allclose(temp_fbound,temp_drchlt)

        
        assert num.allclose(domain_fbound.quantities['stage'].vertex_values,
                            domain_drchlt.quantities['stage'].vertex_values)
                        
        assert num.allclose(domain_fbound.quantities['xmomentum'].vertex_values,
                            domain_drchlt.quantities['xmomentum'].vertex_values)                        
                        
        assert num.allclose(domain_fbound.quantities['ymomentum'].vertex_values,
                            domain_drchlt.quantities['ymomentum'].vertex_values)
        
        # Use known Dirichlet condition (if sufficient timesteps have been taken)

        #FIXME: What do these assertions test? Also do they assume tide =0
        #print domain_fbound.quantities['stage'].vertex_values
        #assert allclose(domain_drchlt.quantities['stage'].vertex_values[6], 2)        
        #assert allclose(domain_fbound.quantities['stage'].vertex_values[6], 2)
        
        

        try:
            os.remove(sts_file+'.sts')
        except:
            # Windoze can't remove this file for some reason 
            pass
        
        os.remove(meshname)
        

        
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

    #### TESTS URS UNGRIDDED 2 SWW ###
    def test_URS_points_needed(self):
        
        ll_lat = -21.5
        ll_long = 114.5
        grid_spacing = 1./60.
        lat_amount = 30
        long_amount = 30
        zone = 50

        boundary_polygon = [[250000,7660000],[280000,7660000],
                             [280000,7630000],[250000,7630000]]
        geo=URS_points_needed(boundary_polygon, zone, 
                              ll_lat, ll_long, grid_spacing, 
                              lat_amount, long_amount,
                              verbose=self.verbose)
        # to test this geo, can info from it be transfered onto the boundary
        # poly?
        #Maybe see if I can fit the data to the polygon - have to represent
        # the poly as points though.
        #geo.export_points_file("results.txt", as_lat_long=True)
        results = frozenset(geo.get_data_points(as_lat_long=True))
        #print 'results',results

        # These are a set of points that have to be in results
        points = []
        for i in range(18):
            lat = -21.0 - 8./60 - grid_spacing * i
            points.append((lat,degminsec2decimal_degrees(114,35,0))) 
            points.append((lat,degminsec2decimal_degrees(114,36,0))) 
            points.append((lat,degminsec2decimal_degrees(114,52,0))) 
            points.append((lat,degminsec2decimal_degrees(114,53,0)))
        geo_answer = Geospatial_data(data_points=points,
                                     points_are_lats_longs=True) 
        #geo_answer.export_points_file("answer.txt", as_lat_long=True)  
        answer = frozenset(points)
        
        outs = answer.difference(results)
        #print "outs", outs
        # This doesn't work.  Though visualising the results shows that
        # it is correct.
        #assert answer.issubset(results)
        # this is why;
        #point (-21.133333333333333, 114.58333333333333)
        #result (-21.133333332232368, 114.58333333300342)
        
        for point in points:
            found = False
            for result in results:
                if num.allclose(point, result):
                    found = True
                    break
            if not found:
                assert False
        
    
    def dave_test_URS_points_needed(self):
        ll_lat = -21.51667
        ll_long = 114.51667
        grid_spacing = 2./60.
        lat_amount = 15
        long_amount = 15

       
        boundary_polygon = [[250000,7660000],[280000,7660000],
                             [280000,7630000],[250000,7630000]]
        URS_points_needed_to_file('a_test_example',boundary_polygon,
                                  ll_lat, ll_long, grid_spacing, 
                                  lat_amount, long_amount,
                                  verbose=self.verbose)
        
    def X_test_URS_points_neededII(self):
        ll_lat = -21.5
        ll_long = 114.5
        grid_spacing = 1./60.
        lat_amount = 30
        long_amount = 30

        # change this so lats and longs are inputed, then converted
        
        #boundary_polygon = [[7660000,250000],[7660000,280000],
        #                     [7630000,280000],[7630000,250000]]
        URS_points_needed(boundary_polygon, ll_lat, ll_long, grid_spacing, 
                          lat_amount, long_amount,
                          verbose=self.verbose)
        
    def test_URS_points_northern_hemisphere(self):
               
        LL_LAT = 8.0
        LL_LONG = 97.0
        GRID_SPACING = 2.0/60.0
        LAT_AMOUNT = 2
        LONG_AMOUNT = 2
        ZONE = 47

        # 
        points = []
        for i in range(2):
            for j in range(2):
                points.append((degminsec2decimal_degrees(8,1+i*2,0),
                               degminsec2decimal_degrees(97,1+i*2,0)))
        #print "points", points
        geo_poly = Geospatial_data(data_points=points,
                                     points_are_lats_longs=True)
        poly_lat_long = geo_poly.get_data_points(as_lat_long=False,
                                       isSouthHemisphere=False)
        #print "seg_lat_long",  poly_lat_long
        
      #   geo=URS_points_needed_to_file('test_example_poly3', poly_lat_long,
#                                   ZONE,
#                                   LL_LAT, LL_LONG,
#                                   GRID_SPACING,
#                                   LAT_AMOUNT, LONG_AMOUNT,
#                                   isSouthernHemisphere=False,
#                                   export_csv=True,
#                                   verbose=self.verbose)
        
        geo=URS_points_needed(poly_lat_long,
                                  ZONE,
                                  LL_LAT, LL_LONG,
                                  GRID_SPACING,
                                  LAT_AMOUNT, LONG_AMOUNT,
                                  isSouthHemisphere=False,
                                  verbose=self.verbose)
        
        results = frozenset(geo.get_data_points(as_lat_long=True,
                                                isSouthHemisphere=False))
        #print 'results',results

        # These are a set of points that have to be in results
        points = [] 
        for i in range(2):
            for j in range(2):
                points.append((degminsec2decimal_degrees(8,i*2,0),
                               degminsec2decimal_degrees(97,i*2,0)))
        #print "answer points", points
        answer = frozenset(points)
        
        for point in points:
            found = False
            for result in results:
                if num.allclose(point, result):
                    found = True
                    break
            if not found:
                assert False
        

    def covered_in_other_tests_test_URS_points_needed_poly1(self):
        # Values used for FESA 2007 results
        # domain in southern hemisphere zone 51        
        LL_LAT = -50.0
        LL_LONG = 80.0
        GRID_SPACING = 2.0/60.0
        LAT_AMOUNT = 4800
        LONG_AMOUNT = 3600
        ZONE = 51
        
        poly1 = [[296361.89, 8091928.62],
                 [429495.07,8028278.82],
                 [447230.56,8000674.05],
                 [429661.2,7982177.6],
                 [236945.9,7897453.16],
                 [183493.44,7942782.27],
                 [226583.04,8008058.96]]

        URS_points_needed_to_file('test_example_poly2', poly1,
                                  ZONE,
                                  LL_LAT, LL_LONG,
                                  GRID_SPACING,
                                  LAT_AMOUNT, LONG_AMOUNT,
                                  verbose=self.verbose)
        


    def covered_in_other_tests_test_URS_points_needed_poly2(self):
        # Values used for 2004 validation work
        # domain in northern hemisphere zone 47        
        LL_LAT = 0.0
        LL_LONG = 90.0
        GRID_SPACING = 2.0/60.0
        LAT_AMOUNT = (15-LL_LAT)/GRID_SPACING
        LONG_AMOUNT = (100-LL_LONG)/GRID_SPACING 
        ZONE = 47
        
        poly2 = [[419336.424,810100.845],
                 [342405.0281,711455.8026],
                 [274649.9152,723352.9603],
                 [272089.092,972393.0131],
                 [347633.3754,968551.7784],
                 [427979.2022,885965.2313],
                 [427659.0993,875721.9386],
                 [429259.6138,861317.3083],
                 [436301.8775,840830.723]]
        
        URS_points_needed_to_file('test_example_poly2', poly2,
                                  ZONE,
                                  LL_LAT, LL_LONG,
                                  GRID_SPACING,
                                  LAT_AMOUNT, LONG_AMOUNT,
                                  isSouthernHemisphere=False,
                                  verbose=self.verbose) 
        
    #### END TESTS URS UNGRIDDED 2 SWW ###
    def test_Urs_points(self):
        time_step_count = 3
        time_step = 2
        lat_long_points =[(-21.5,114.5),(-21.5,115),(-21.,115)]
        base_name, files = self.write_mux(lat_long_points,
                                          time_step_count, time_step)
        for file in files:
            urs = Urs_points(file)
            assert time_step_count == urs.time_step_count
            assert time_step == urs.time_step

            for lat_lon, dep in map(None, lat_long_points, urs.lonlatdep):
                    _ , e, n = redfearn(lat_lon[0], lat_lon[1])
                    assert num.allclose(n, dep[2])
                        
            count = 0
            for slice in urs:
                count += 1
                #print slice
                for lat_lon, quantity in map(None, lat_long_points, slice):
                    _ , e, n = redfearn(lat_lon[0], lat_lon[1])
                    #print "quantity", quantity
                    #print "e", e
                    #print "n", n
                    if file[-5:] == WAVEHEIGHT_MUX_LABEL[-5:] or \
                           file[-5:] == NORTH_VELOCITY_LABEL[-5:] :
                        assert num.allclose(e, quantity)
                    if file[-5:] == EAST_VELOCITY_LABEL[-5:]:
                        assert num.allclose(n, quantity)
            assert count == time_step_count
                     
        self.delete_mux(files)

    def test_urs_ungridded2sww (self):
        
        #Zone:   50    
        #Easting:  240992.578  Northing: 7620442.472 
        #Latitude:   -21  30 ' 0.00000 ''  Longitude: 114  30 ' 0.00000 '' 
        lat_long = [[-21.5,114.5],[-21,114.5],[-21,115]]
        time_step_count = 2
        time_step = 400
        tide = 9000000
        base_name, files = self.write_mux(lat_long,
                                          time_step_count, time_step)
        urs_ungridded2sww(base_name, mean_stage=tide,
                          verbose=self.verbose)
        
        # now I want to check the sww file ...
        sww_file = base_name + '.sww'
        
        #Let's interigate the sww file
        # Note, the sww info is not gridded.  It is point data.
        fid = NetCDFFile(sww_file)
        
        # Make x and y absolute
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        geo_reference = Geo_reference(NetCDFObject=fid)
        points = geo_reference.get_absolute(map(None, x, y))
        points = ensure_numeric(points)
        x = points[:,0]
        y = points[:,1]
        
        #Check that first coordinate is correctly represented       
        #Work out the UTM coordinates for first point
        zone, e, n = redfearn(lat_long[0][0], lat_long[0][1]) 
        assert num.allclose([x[0],y[0]], [e,n])

        #Check the time vector
        times = fid.variables['time'][:]
        
        times_actual = []
        for i in range(time_step_count):
            times_actual.append(time_step * i)
        
        assert num.allclose(ensure_numeric(times),
                            ensure_numeric(times_actual))
        
        #Check first value
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        ymomentum = fid.variables['ymomentum'][:]
        elevation = fid.variables['elevation'][:]
        assert num.allclose(stage[0,0], e +tide)  #Meters


        #Check the momentums - ua
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_ua *(stage+depth)
        # = n*(e+tide+n) based on how I'm writing these files
        # 
        answer_x = n*(e+tide+n)
        actual_x = xmomentum[0,0]
        #print "answer_x",answer_x
        #print "actual_x",actual_x 
        assert num.allclose(answer_x, actual_x)  #Meters
        
        #Check the momentums - va
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_va *(stage+depth)
        # = e*(e+tide+n) based on how I'm writing these files
        # 
        answer_y = -1*e*(e+tide+n)
        actual_y = ymomentum[0,0]
        #print "answer_y",answer_y
        #print "actual_y",actual_y 
        assert num.allclose(answer_y, actual_y)  #Meters

        # check the stage values, first time step.
        # These arrays are equal since the Easting values were used as
        # the stage
        assert num.allclose(stage[0], x +tide)  #Meters
        # check the elevation values.
        # -ve since urs measures depth, sww meshers height,
        # these arrays are equal since the northing values were used as
        # the elevation
        assert num.allclose(-elevation, y)  #Meters
        
        fid.close()
        self.delete_mux(files)
        os.remove(sww_file)
  
    def test_urs_ungridded2swwII (self):
        
        #Zone:   50    
        #Easting:  240992.578  Northing: 7620442.472 
        #Latitude:   -21  30 ' 0.00000 ''  Longitude: 114  30 ' 0.00000 '' 
        lat_long = [[-21.5,114.5],[-21,114.5],[-21,115]]
        time_step_count = 2
        time_step = 400
        tide = 9000000
        geo_reference = Geo_reference(50, 3434543,34534543)
        base_name, files = self.write_mux(lat_long,
                                          time_step_count, time_step)
        urs_ungridded2sww(base_name, mean_stage=tide, origin = geo_reference,
                          verbose=self.verbose)
        
        # now I want to check the sww file ...
        sww_file = base_name + '.sww'
        
        #Let's interigate the sww file
        # Note, the sww info is not gridded.  It is point data.
        fid = NetCDFFile(sww_file)
        
        # Make x and y absolute
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        geo_reference = Geo_reference(NetCDFObject=fid)
        points = geo_reference.get_absolute(map(None, x, y))
        points = ensure_numeric(points)
        x = points[:,0]
        y = points[:,1]
        
        #Check that first coordinate is correctly represented       
        #Work out the UTM coordinates for first point
        zone, e, n = redfearn(lat_long[0][0], lat_long[0][1]) 
        assert num.allclose([x[0],y[0]], [e,n])

        #Check the time vector
        times = fid.variables['time'][:]
        
        times_actual = []
        for i in range(time_step_count):
            times_actual.append(time_step * i)
        
        assert num.allclose(ensure_numeric(times),
                            ensure_numeric(times_actual))
        
        #Check first value
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        ymomentum = fid.variables['ymomentum'][:]
        elevation = fid.variables['elevation'][:]
        assert num.allclose(stage[0,0], e +tide)  #Meters

        #Check the momentums - ua
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_ua *(stage+depth)
        # = n*(e+tide+n) based on how I'm writing these files
        # 
        answer_x = n*(e+tide+n)
        actual_x = xmomentum[0,0]
        #print "answer_x",answer_x
        #print "actual_x",actual_x 
        assert num.allclose(answer_x, actual_x)  #Meters
        
        #Check the momentums - va
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_va *(stage+depth)
        # = e*(e+tide+n) based on how I'm writing these files
        # 
        answer_y = -1*e*(e+tide+n)
        actual_y = ymomentum[0,0]
        #print "answer_y",answer_y
        #print "actual_y",actual_y 
        assert num.allclose(answer_y, actual_y)  #Meters

        # check the stage values, first time step.
        # These arrays are equal since the Easting values were used as
        # the stage
        assert num.allclose(stage[0], x +tide)  #Meters
        # check the elevation values.
        # -ve since urs measures depth, sww meshers height,
        # these arrays are equal since the northing values were used as
        # the elevation
        assert num.allclose(-elevation, y)  #Meters
        
        fid.close()
        self.delete_mux(files)
        os.remove(sww_file)
  
    def test_urs_ungridded2swwIII (self):
        
        #Zone:   50    
        #Easting:  240992.578  Northing: 7620442.472 
        #Latitude:   -21  30 ' 0.00000 ''  Longitude: 114  30 ' 0.00000 '' 
        lat_long = [[-21.5,114.5],[-21,114.5],[-21,115]]
        time_step_count = 2
        time_step = 400
        tide = 9000000
        base_name, files = self.write_mux(lat_long,
                                          time_step_count, time_step)
        urs_ungridded2sww(base_name, mean_stage=tide, origin =(50,23432,4343),
                          verbose=self.verbose)
        
        # now I want to check the sww file ...
        sww_file = base_name + '.sww'
        
        #Let's interigate the sww file
        # Note, the sww info is not gridded.  It is point data.
        fid = NetCDFFile(sww_file)
        
        # Make x and y absolute
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        geo_reference = Geo_reference(NetCDFObject=fid)
        points = geo_reference.get_absolute(map(None, x, y))
        points = ensure_numeric(points)
        x = points[:,0]
        y = points[:,1]
        
        #Check that first coordinate is correctly represented       
        #Work out the UTM coordinates for first point
        zone, e, n = redfearn(lat_long[0][0], lat_long[0][1]) 
        assert num.allclose([x[0],y[0]], [e,n])

        #Check the time vector
        times = fid.variables['time'][:]
        
        times_actual = []
        for i in range(time_step_count):
            times_actual.append(time_step * i)
        
        assert num.allclose(ensure_numeric(times),
                            ensure_numeric(times_actual))
        
        #Check first value
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        ymomentum = fid.variables['ymomentum'][:]
        elevation = fid.variables['elevation'][:]
        assert num.allclose(stage[0,0], e +tide)  #Meters

        #Check the momentums - ua
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_ua *(stage+depth)
        # = n*(e+tide+n) based on how I'm writing these files
        # 
        answer_x = n*(e+tide+n)
        actual_x = xmomentum[0,0]
        #print "answer_x",answer_x
        #print "actual_x",actual_x 
        assert num.allclose(answer_x, actual_x)  #Meters
        
        #Check the momentums - va
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_va *(stage+depth)
        # = e*(e+tide+n) based on how I'm writing these files
        # 
        answer_y = -1*e*(e+tide+n)
        actual_y = ymomentum[0,0]
        #print "answer_y",answer_y
        #print "actual_y",actual_y 
        assert num.allclose(answer_y, actual_y)  #Meters

        # check the stage values, first time step.
        # These arrays are equal since the Easting values were used as
        # the stage
        assert num.allclose(stage[0], x +tide)  #Meters
        # check the elevation values.
        # -ve since urs measures depth, sww meshers height,
        # these arrays are equal since the northing values were used as
        # the elevation
        assert num.allclose(-elevation, y)  #Meters
        
        fid.close()
        self.delete_mux(files)
        os.remove(sww_file)

        
    def test_urs_ungridded_hole (self):
        
        #Zone:   50    
        #Easting:  240992.578  Northing: 7620442.472 
        #Latitude:   -21  30 ' 0.00000 ''  Longitude: 114  30 ' 0.00000 '' 
        lat_long = [[-20.5, 114.5],
                    [-20.6, 114.6],
                    [-20.5, 115.],
                    [-20.6, 115.],
                    [-20.5, 115.5],
                    [-20.6, 115.4],
                    
                    [-21., 114.5],
                    [-21., 114.6],
                    [-21., 115.5],
                    [-21., 115.4],
                    
                    [-21.5, 114.5],
                    [-21.4, 114.6],
                    [-21.5, 115.],
                    [-21.4, 115.],
                    [-21.5, 115.5],
                    [-21.4, 115.4]
                    ]
        time_step_count = 6
        time_step = 100
        tide = 9000000
        base_name, files = self.write_mux(lat_long,
                                          time_step_count, time_step)
        #Easting:  292110.784  Northing: 7676551.710 
        #Latitude:   -21  0 ' 0.00000 ''  Longitude: 115  0 ' 0.00000 '' 

        urs_ungridded2sww(base_name, mean_stage=-240992.0,
                          hole_points_UTM=[ 292110.784, 7676551.710 ],
                          verbose=self.verbose)
        
        # now I want to check the sww file ...
        sww_file = base_name + '.sww'
        
        #Let's interigate the sww file
        # Note, the sww info is not gridded.  It is point data.
        fid = NetCDFFile(sww_file)
        
        number_of_volumes = fid.variables['volumes']
        #print "number_of_volumes",len(number_of_volumes) 
        assert num.allclose(16, len(number_of_volumes))
        
        fid.close()
        self.delete_mux(files)
        #print "sww_file", sww_file 
        os.remove(sww_file)
        
    def test_urs_ungridded_holeII(self):

        # Check that if using a hole that returns no triangles,
        # urs_ungridded2sww removes the hole label.
        
        lat_long = [[-20.5, 114.5],
                    [-20.6, 114.6],
                    [-20.5, 115.5],
                    [-20.6, 115.4],
                    
                    
                    [-21.5, 114.5],
                    [-21.4, 114.6],
                    [-21.5, 115.5],
                    [-21.4, 115.4]
                    ]
        time_step_count = 6
        time_step = 100
        tide = 9000000
        base_name, files = self.write_mux(lat_long,
                                          time_step_count, time_step)
        #Easting:  292110.784  Northing: 7676551.710 
        #Latitude:   -21  0 ' 0.00000 ''  Longitude: 115  0 ' 0.00000 '' 

        urs_ungridded2sww(base_name, mean_stage=-240992.0,
                          hole_points_UTM=[ 292110.784, 7676551.710 ],
                          verbose=self.verbose)
        
        # now I want to check the sww file ...
        sww_file = base_name + '.sww'
        fid = NetCDFFile(sww_file)
        
        volumes = fid.variables['volumes']
        #print "number_of_volumes",len(volumes)

        fid.close()
        os.remove(sww_file)
        
        urs_ungridded2sww(base_name, mean_stage=-240992.0)
        
        # now I want to check the sww file ...
        sww_file = base_name + '.sww'
        fid = NetCDFFile(sww_file)
        
        volumes_again = fid.variables['volumes']
        #print "number_of_volumes",len(volumes_again) 
        assert num.allclose(len(volumes_again),
                            len(volumes))
        fid.close()
        os.remove(sww_file)
        self.delete_mux(files) 
        
    def test_urs_ungridded2sww_mint_maxt (self):
        
        #Zone:   50    
        #Easting:  240992.578  Northing: 7620442.472 
        #Latitude:   -21  30 ' 0.00000 ''  Longitude: 114  30 ' 0.00000 '' 
        lat_long = [[-21.5,114.5],[-21,114.5],[-21,115]]
        time_step_count = 6
        time_step = 100
        tide = 9000000
        base_name, files = self.write_mux(lat_long,
                                          time_step_count, time_step)
        urs_ungridded2sww(base_name, mean_stage=tide, origin =(50,23432,4343),
                          mint=101, maxt=500,
                          verbose=self.verbose)
        
        # now I want to check the sww file ...
        sww_file = base_name + '.sww'
        
        #Let's interigate the sww file
        # Note, the sww info is not gridded.  It is point data.
        fid = NetCDFFile(sww_file)
        
        # Make x and y absolute
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        geo_reference = Geo_reference(NetCDFObject=fid)
        points = geo_reference.get_absolute(map(None, x, y))
        points = ensure_numeric(points)
        x = points[:,0]
        y = points[:,1]
        
        #Check that first coordinate is correctly represented       
        #Work out the UTM coordinates for first point
        zone, e, n = redfearn(lat_long[0][0], lat_long[0][1]) 
        assert num.allclose([x[0],y[0]], [e,n])

        #Check the time vector
        times = fid.variables['time'][:]
        
        times_actual = [0,100,200,300]
       
        assert num.allclose(ensure_numeric(times),
                            ensure_numeric(times_actual))
        
               #Check first value
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        ymomentum = fid.variables['ymomentum'][:]
        elevation = fid.variables['elevation'][:]
        assert num.allclose(stage[0,0], e +tide)  #Meters

        #Check the momentums - ua
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_ua *(stage+depth)
        # = n*(e+tide+n) based on how I'm writing these files
        # 
        answer_x = n*(e+tide+n)
        actual_x = xmomentum[0,0]
        #print "answer_x",answer_x
        #print "actual_x",actual_x 
        assert num.allclose(answer_x, actual_x)  #Meters
        
        #Check the momentums - va
        #momentum = velocity*(stage-elevation)
        # elevation = - depth
        #momentum = velocity_va *(stage+depth)
        # = e*(e+tide+n) based on how I'm writing these files
        # 
        answer_y = -1*e*(e+tide+n)
        actual_y = ymomentum[0,0]
        #print "answer_y",answer_y
        #print "actual_y",actual_y 
        assert num.allclose(answer_y, actual_y)  #Meters

        # check the stage values, first time step.
        # These arrays are equal since the Easting values were used as
        # the stage
        assert num.allclose(stage[0], x +tide)  #Meters
        # check the elevation values.
        # -ve since urs measures depth, sww meshers height,
        # these arrays are equal since the northing values were used as
        # the elevation
        assert num.allclose(-elevation, y)  #Meters
        
        fid.close()
        self.delete_mux(files)
        os.remove(sww_file)
        
    def test_urs_ungridded2sww_mint_maxtII (self):
        
        #Zone:   50    
        #Easting:  240992.578  Northing: 7620442.472 
        #Latitude:   -21  30 ' 0.00000 ''  Longitude: 114  30 ' 0.00000 '' 
        lat_long = [[-21.5,114.5],[-21,114.5],[-21,115]]
        time_step_count = 6
        time_step = 100
        tide = 9000000
        base_name, files = self.write_mux(lat_long,
                                          time_step_count, time_step)
        urs_ungridded2sww(base_name, mean_stage=tide, origin =(50,23432,4343),
                          mint=0, maxt=100000,
                          verbose=self.verbose)
        
        # now I want to check the sww file ...
        sww_file = base_name + '.sww'
        
        #Let's interigate the sww file
        # Note, the sww info is not gridded.  It is point data.
        fid = NetCDFFile(sww_file)
        
        # Make x and y absolute
        geo_reference = Geo_reference(NetCDFObject=fid)
        points = geo_reference.get_absolute(map(None, fid.variables['x'][:],
                                                fid.variables['y'][:]))
        points = ensure_numeric(points)
        x = points[:,0]
        
        #Check the time vector
        times = fid.variables['time'][:]
        
        times_actual = [0,100,200,300,400,500]
        assert num.allclose(ensure_numeric(times),
                            ensure_numeric(times_actual))
        
        #Check first value
        stage = fid.variables['stage'][:]
        assert num.allclose(stage[0], x +tide)
        
        fid.close()
        self.delete_mux(files)
        os.remove(sww_file)
        
    def test_urs_ungridded2sww_mint_maxtIII (self):
        
        #Zone:   50    
        #Easting:  240992.578  Northing: 7620442.472 
        #Latitude:   -21  30 ' 0.00000 ''  Longitude: 114  30 ' 0.00000 '' 
        lat_long = [[-21.5,114.5],[-21,114.5],[-21,115]]
        time_step_count = 6
        time_step = 100
        tide = 9000000
        base_name, files = self.write_mux(lat_long,
                                          time_step_count, time_step)
        try:
            urs_ungridded2sww(base_name, mean_stage=tide,
                          origin =(50,23432,4343),
                          mint=301, maxt=399,
                              verbose=self.verbose)
        except: 
            pass
        else:
            self.failUnless(0 ==1, 'Bad input did not throw exception error!')

        self.delete_mux(files)
        
    def test_urs_ungridded2sww_mint_maxt_bad (self):       
        #Zone:   50    
        #Easting:  240992.578  Northing: 7620442.472 
        #Latitude:   -21  30 ' 0.00000 ''  Longitude: 114  30 ' 0.00000 '' 
        lat_long = [[-21.5,114.5],[-21,114.5],[-21,115]]
        time_step_count = 6
        time_step = 100
        tide = 9000000
        base_name, files = self.write_mux(lat_long,
                                          time_step_count, time_step)
        try:
            urs_ungridded2sww(base_name, mean_stage=tide,
                          origin =(50,23432,4343),
                          mint=301, maxt=301,
                              verbose=self.verbose)
        except: 
            pass
        else:
            self.failUnless(0 ==1, 'Bad input did not throw exception error!')

        self.delete_mux(files)

        
    def test_URS_points_needed_and_urs_ungridded2sww(self):
        # This doesn't actually check anything
        #  
        ll_lat = -21.5
        ll_long = 114.5
        grid_spacing = 1./60.
        lat_amount = 30
        long_amount = 30
        time_step_count = 2
        time_step = 400
        tide = -200000
        zone = 50

        boundary_polygon = [[250000,7660000],[280000,7660000],
                             [280000,7630000],[250000,7630000]]
        geo=URS_points_needed(boundary_polygon, zone,
                              ll_lat, ll_long, grid_spacing, 
                              lat_amount, long_amount,
                              verbose=self.verbose)
        lat_long = geo.get_data_points(as_lat_long=True)
        base_name, files = self.write_mux(lat_long,
                                          time_step_count, time_step)
        urs_ungridded2sww(base_name, mean_stage=tide,
                          verbose=self.verbose)
        self.delete_mux(files)
        os.remove( base_name + '.sww')
    
    def cache_test_URS_points_needed_and_urs_ungridded2sww(self):
        
        ll_lat = -21.5
        ll_long = 114.5
        grid_spacing = 1./60.
        lat_amount = 30
        long_amount = 30
        time_step_count = 2
        time_step = 400
        tide = -200000
        zone = 50

        boundary_polygon = [[250000,7660000],[270000,7650000],
                             [280000,7630000],[250000,7630000]]
        geo=URS_points_needed(boundary_polygon, zone,
                              ll_lat, ll_long, grid_spacing, 
                              lat_amount, long_amount, use_cache=True,
                              verbose=True)
        
    def visual_test_URS_points_needed_and_urs_ungridded2sww(self):
        
        ll_lat = -21.5
        ll_long = 114.5
        grid_spacing = 1./60.
        lat_amount = 30
        long_amount = 30
        time_step_count = 2
        time_step = 400
        tide = -200000
        zone = 50

        boundary_polygon = [[250000,7660000],[270000,7650000],
                             [280000,7630000],[250000,7630000]]
        geo=URS_points_needed(boundary_polygon, zone,
                              ll_lat, ll_long, grid_spacing, 
                              lat_amount, long_amount)
        lat_long = geo.get_data_points(as_lat_long=True)
        base_name, files = self.write_mux(lat_long,
                                          time_step_count, time_step)
        urs_ungridded2sww(base_name, mean_stage=tide)
        self.delete_mux(files)
        os.remove( base_name + '.sww')
        # extend this so it interpolates onto the boundary.
        # have it fail if there is NaN

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
        
    def test_get_data_from_file(self):
#    from anuga.abstract_2d_finite_volumes.util import get_data_from_file
        
        import os
       
        fileName = tempfile.mktemp(".txt")
#        print"filename",fileName
        file = open(fileName,"w")
        file.write("elevation, stage\n\
1.0, 3  \n\
0.0, 4 \n\
4.0, 3 \n\
1.0, 6 \n")
        file.close()
        
        header,x = get_data_from_file(fileName)
#        print 'x',x
        os.remove(fileName)
        
        assert num.allclose(x[:,0], [1.0, 0.0,4.0, 1.0])
        
    def test_get_data_from_file1(self):
#    from anuga.abstract_2d_finite_volumes.util import get_data_from_file
        
        import os
       
        fileName = tempfile.mktemp(".txt")
#        print"filename",fileName
        file = open(fileName,"w")
        file.write("elevation stage\n\
1.3 3  \n\
0.0 4 \n\
4.5 3.5 \n\
1.0 6 \n")
        file.close()
        
        header, x = get_data_from_file(fileName,separator_value=' ')
        os.remove(fileName)
#        x = get_data_from_file(fileName)
#        print '1x',x[:,0]
        
        assert num.allclose(x[:,0], [1.3, 0.0,4.5, 1.0])
        
    def test_store_parameters(self):
        """tests store temporary file
        """
        
        from os import sep, getenv
        
        output_dir=''
        file_name='details.csv'
        
        kwargs = {'file_name':'new2.txt',
                  'output_dir':output_dir,
                  'file_name':file_name,
                  'who':'me',
                  'what':'detail',
                  'how':2,
                  'why':241,
#                  'completed':345
                  }
        store_parameters(verbose=False,**kwargs)

        temp='detail_temp.csv'
        fid = open(temp)
        file_header = fid.readline()
        file_line = fid.readline()
        fid.close()
        
        
        keys = kwargs.keys()
        keys.sort()
        line=''
        header=''
        count=0
        #used the sorted keys to create the header and line data
        for k in keys:
#            print "%s = %s" %(k, kwargs[k]) 
            header = header+str(k)
            line = line+str(kwargs[k])
            count+=1
            if count <len(kwargs):
                header = header+','
                line = line+','
        header+='\n'
        line+='\n'
        
        
        #file exists
        assert access(temp,F_OK)
        assert header == file_header
        assert line == file_line
        
        os.remove(temp)
        
    def test_store_parameters1(self):
        """tests store in temporary file and other file 
        """
        
        from os import sep, getenv
        
        output_dir=''
        file_name='details.csv'
        
        kwargs = {'file_name':'new2.txt',
                  'output_dir':output_dir,
                  'file_name':file_name,
                  'who':'me',
                  'what':'detail',
                  'how':2,
                  'why':241,
#                  'completed':345
                  }
        store_parameters(verbose=False,**kwargs)
        
        kwargs['how']=55
        kwargs['completed']=345

        keys = kwargs.keys()
        keys.sort()
        line=''
        header=''
        count=0
        #used the sorted keys to create the header and line data
        for k in keys:
#            print "%s = %s" %(k, kwargs[k]) 
            header = header+str(k)
            line = line+str(kwargs[k])
            count+=1
            if count <len(kwargs):
                header = header+','
                line = line+','
        header+='\n'
        line+='\n'
        
        kwargs['how']=55
        kwargs['completed']=345
        
        store_parameters(verbose=False,**kwargs)
        
#        temp='detail_temp.csv'
        fid = open(file_name)
        file_header = fid.readline()
        file_line1 = fid.readline()
        file_line2 = fid.readline()
        fid.close()
        
        
        #file exists
#        print 'header',header,'line',line
#        print 'file_header',file_header,'file_line1',file_line1,'file_line2',file_line2
        assert access(file_name,F_OK)
        assert header == file_header
        assert line == file_line1
        
        temp='detail_temp.csv'
        os.remove(temp)
        os.remove(file_name)        
        
    def test_store_parameters2(self):
        """tests appending the data to the end of an existing file
        """
        
        from os import sep, getenv
        
        output_dir=''
        file_name='details.csv'
        
        kwargs = {'file_name':'new2.txt',
                  'output_dir':output_dir,
                  'file_name':file_name,
                  'who':'me',
                  'what':'detail',
                  'how':2,
                  'why':241,
                  'completed':345
                  }
        store_parameters(verbose=False,**kwargs)
        
        kwargs['how']=55
        kwargs['completed']=23.54532
        
        store_parameters(verbose=False,**kwargs)
        
        keys = kwargs.keys()
        keys.sort()
        line=''
        header=''
        count=0
        #used the sorted keys to create the header and line data
        for k in keys:
#            print "%s = %s" %(k, kwargs[k]) 
            header = header+str(k)
            line = line+str(kwargs[k])
            count+=1
            if count <len(kwargs):
                header = header+','
                line = line+','
        header+='\n'
        line+='\n'
        
        fid = open(file_name)
        file_header = fid.readline()
        file_line1 = fid.readline()
        file_line2 = fid.readline()
        fid.close()
        
        assert access(file_name,F_OK)
        assert header == file_header
        assert line == file_line2
        
        os.remove(file_name)        
        

    def test_get_maximum_inundation(self):
        """Test that sww information can be converted correctly to maximum
        runup elevation and location (without and with georeferencing)

        This test creates a slope and a runup which is maximal (~11m) at around 10s
        and levels out to the boundary condition (1m) at about 30s.
        """

        import time, os
        from Scientific.IO.NetCDF import NetCDFFile

        #Setup

        from mesh_factory import rectangular

        # Create basic mesh (100m x 100m)
        points, vertices, boundary = rectangular(20, 5, 100, 50)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.default_order = 2
        domain.set_minimum_storable_height(0.01)

        domain.set_name('runuptest')
        swwfile = domain.get_name() + '.sww'

        domain.set_datadir('.')
        domain.format = 'sww'
        domain.smooth = True

        # FIXME (Ole): Backwards compatibility
        # Look at sww file and see what happens when
        # domain.tight_slope_limiters = 1
        domain.tight_slope_limiters = 0
        domain.use_centroid_velocities = 0 # Backwards compatibility (7/5/8)        
        
        Br = Reflective_boundary(domain)
        Bd = Dirichlet_boundary([1.0,0,0])


        #---------- First run without geo referencing
        
        domain.set_quantity('elevation', lambda x,y: -0.2*x + 14) # Slope
        domain.set_quantity('stage', -6)
        domain.set_boundary( {'left': Br, 'right': Bd, 'top': Br, 'bottom': Br})

        for t in domain.evolve(yieldstep=1, finaltime = 50):
            pass


        # Check maximal runup
        runup = get_maximum_inundation_elevation(swwfile)
        location = get_maximum_inundation_location(swwfile)
        #print 'Runup, location', runup, location
        assert num.allclose(runup, 11) or num.allclose(runup, 12) # old limiters
        assert num.allclose(location[0], 15) or num.allclose(location[0], 10)

        # Check final runup
        runup = get_maximum_inundation_elevation(swwfile, time_interval=[45,50])
        location = get_maximum_inundation_location(swwfile, time_interval=[45,50])
        # print 'Runup, location:',runup, location        
        assert num.allclose(runup, 1)
        assert num.allclose(location[0], 65)

        # Check runup restricted to a polygon
        p = [[50,1], [99,1], [99,49], [50,49]]
        runup = get_maximum_inundation_elevation(swwfile, polygon=p)
        location = get_maximum_inundation_location(swwfile, polygon=p)
        #print runup, location        
        assert num.allclose(runup, 4)
        assert num.allclose(location[0], 50)                

        # Check that mimimum_storable_height works
        fid = NetCDFFile(swwfile, netcdf_mode_r) # Open existing file
        
        stage = fid.variables['stage'][:]
        z = fid.variables['elevation'][:]
        xmomentum = fid.variables['xmomentum'][:]
        ymomentum = fid.variables['ymomentum'][:]        

        
        
        for i in range(stage.shape[0]):
            h = stage[i]-z # depth vector at time step i
            
            # Check every node location
            for j in range(stage.shape[1]):
                # Depth being either exactly zero implies
                # momentum being zero.
                # Or else depth must be greater than or equal to
                # the minimal storable height
                if h[j] == 0.0:
                    assert xmomentum[i,j] == 0.0
                    assert ymomentum[i,j] == 0.0                
                else:
                    assert h[j] >= domain.minimum_storable_height
        
        fid.close()

        # Cleanup
        os.remove(swwfile)
        


        #------------- Now the same with georeferencing

        domain.time=0.0
        E = 308500
        N = 6189000
        #E = N = 0
        domain.geo_reference = Geo_reference(56, E, N)

        domain.set_quantity('elevation', lambda x,y: -0.2*x + 14) # Slope
        domain.set_quantity('stage', -6)
        domain.set_boundary( {'left': Br, 'right': Bd, 'top': Br, 'bottom': Br})

        for t in domain.evolve(yieldstep=1, finaltime = 50):
            pass

        # Check maximal runup
        runup = get_maximum_inundation_elevation(swwfile)
        location = get_maximum_inundation_location(swwfile)
        assert num.allclose(runup, 11) or num.allclose(runup, 12) # old limiters
        assert num.allclose(location[0], 15+E) or num.allclose(location[0], 10+E)

        # Check final runup
        runup = get_maximum_inundation_elevation(swwfile, time_interval=[45,50])
        location = get_maximum_inundation_location(swwfile, time_interval=[45,50])
        assert num.allclose(runup, 1)
        assert num.allclose(location[0], 65+E)

        # Check runup restricted to a polygon
        p = num.array([[50,1], [99,1], [99,49], [50,49]], num.int) + num.array([E, N], num.int)      #array default#

        runup = get_maximum_inundation_elevation(swwfile, polygon=p)
        location = get_maximum_inundation_location(swwfile, polygon=p)
        assert num.allclose(runup, 4)
        assert num.allclose(location[0], 50+E)                


        # Cleanup
        os.remove(swwfile)


    def test_get_mesh_and_quantities_from_sww_file(self):
        """test_get_mesh_and_quantities_from_sww_file(self):
        """     
        
        # Generate a test sww file with non trivial georeference
        
        import time, os
        from Scientific.IO.NetCDF import NetCDFFile

        # Setup
        from mesh_factory import rectangular

        # Create basic mesh (100m x 5m)
        width = 5
        length = 50
        t_end = 10
        points, vertices, boundary = rectangular(length, width, 50, 5)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary,
                        geo_reference = Geo_reference(56,308500,6189000))

        domain.set_name('flowtest')
        swwfile = domain.get_name() + '.sww'
        domain.set_datadir('.')

        Br = Reflective_boundary(domain)    # Side walls
        Bd = Dirichlet_boundary([1, 0, 0])  # inflow

        domain.set_boundary( {'left': Bd, 'right': Bd, 'top': Br, 'bottom': Br})

        for t in domain.evolve(yieldstep=1, finaltime = t_end):
            pass

        
        # Read it

        # Get mesh and quantities from sww file
        X = get_mesh_and_quantities_from_file(swwfile,
                                              quantities=['elevation',
                                                          'stage',
                                                          'xmomentum',
                                                          'ymomentum'], 
                                              verbose=False)
        mesh, quantities, time = X
        

        # Check that mesh has been recovered
        assert num.alltrue(mesh.triangles == domain.get_triangles())
        assert num.allclose(mesh.nodes, domain.get_nodes())

        # Check that time has been recovered
        assert num.allclose(time, range(t_end+1))

        # Check that quantities have been recovered
        # (sww files use single precision)
        z=domain.get_quantity('elevation').get_values(location='unique vertices')
        assert num.allclose(quantities['elevation'], z)

        for q in ['stage', 'xmomentum', 'ymomentum']:
            # Get quantity at last timestep
            q_ref=domain.get_quantity(q).get_values(location='unique vertices')

            #print q,quantities[q]
            q_sww=quantities[q][-1,:]

            msg = 'Quantity %s failed to be recovered' %q
            assert num.allclose(q_ref, q_sww, atol=1.0e-6), msg
            
        # Cleanup
        os.remove(swwfile)
        

    def test_get_flow_through_cross_section(self):
        """test_get_flow_through_cross_section(self):

        Test that the total flow through a cross section can be
        correctly obtained from an sww file.
        
        This test creates a flat bed with a known flow through it and tests
        that the function correctly returns the expected flow.

        The specifics are
        u = 2 m/s
        h = 1 m
        w = 3 m (width of channel)

        q = u*h*w = 6 m^3/s

        #---------- First run without geo referencing        
        
        """

        import time, os
        from Scientific.IO.NetCDF import NetCDFFile

        # Setup
        from mesh_factory import rectangular

        # Create basic mesh (20m x 3m)
        width = 3
        length = 20
        t_end = 3
        points, vertices, boundary = rectangular(length, width,
                                                 length, width)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.default_order = 2
        domain.set_minimum_storable_height(0.01)

        domain.set_name('flowtest')
        swwfile = domain.get_name() + '.sww'

        domain.set_datadir('.')
        domain.format = 'sww'
        domain.smooth = True

        h = 1.0
        u = 2.0
        uh = u*h

        Br = Reflective_boundary(domain)     # Side walls
        Bd = Dirichlet_boundary([h, uh, 0])  # 2 m/s across the 3 m inlet: 


        
        domain.set_quantity('elevation', 0.0)
        domain.set_quantity('stage', h)
        domain.set_quantity('xmomentum', uh)
        domain.set_boundary( {'left': Bd, 'right': Bd, 'top': Br, 'bottom': Br})

        for t in domain.evolve(yieldstep=1, finaltime = t_end):
            pass

        # Check that momentum is as it should be in the interior

        I = [[0, width/2.],
             [length/2., width/2.],
             [length, width/2.]]
        
        f = file_function(swwfile,
                          quantities=['stage', 'xmomentum', 'ymomentum'],
                          interpolation_points=I,
                          verbose=False)
        for t in range(t_end+1):
            for i in range(3):
                assert num.allclose(f(t, i), [1, 2, 0], atol=1.0e-6)
            

        # Check flows through the middle
        for i in range(5):
            x = length/2. + i*0.23674563 # Arbitrary
            cross_section = [[x, 0], [x, width]]
            time, Q = get_flow_through_cross_section(swwfile,
                                                     cross_section,
                                                     verbose=False)

            assert num.allclose(Q, uh*width)


       
        # Try the same with partial lines
        x = length/2.
        for i in range(5):
            start_point = [length/2., i*width/5.]
            #print start_point
                            
            cross_section = [start_point, [length/2., width]]
            time, Q = get_flow_through_cross_section(swwfile,
                                                     cross_section,
                                                     verbose=False)

            #print i, Q, (width-start_point[1])
            assert num.allclose(Q, uh*(width-start_point[1]))


        # Verify no flow when line is parallel to flow
        cross_section = [[length/2.-10, width/2.], [length/2.+10, width/2.]]
        time, Q = get_flow_through_cross_section(swwfile,
                                                 cross_section,
                                                 verbose=False)

        #print i, Q
        assert num.allclose(Q, 0, atol=1.0e-5)        


        # Try with lines on an angle (all flow still runs through here)
        cross_section = [[length/2., 0], [length/2.+width, width]]
        time, Q = get_flow_through_cross_section(swwfile,
                                                 cross_section,
                                                 verbose=False)

        assert num.allclose(Q, uh*width)        
        


                                      
    def test_get_flow_through_cross_section_with_geo(self):
        """test_get_flow_through_cross_section(self):

        Test that the total flow through a cross section can be
        correctly obtained from an sww file.
        
        This test creates a flat bed with a known flow through it and tests
        that the function correctly returns the expected flow.

        The specifics are
        u = 2 m/s
        h = 2 m
        w = 3 m (width of channel)

        q = u*h*w = 12 m^3/s


        This run tries it with georeferencing and with elevation = -1
        
        """

        import time, os
        from Scientific.IO.NetCDF import NetCDFFile

        # Setup
        from mesh_factory import rectangular

        # Create basic mesh (20m x 3m)
        width = 3
        length = 20
        t_end = 1
        points, vertices, boundary = rectangular(length, width,
                                                 length, width)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary,
                        geo_reference = Geo_reference(56,308500,6189000))

        domain.default_order = 2
        domain.set_minimum_storable_height(0.01)

        domain.set_name('flowtest')
        swwfile = domain.get_name() + '.sww'

        domain.set_datadir('.')
        domain.format = 'sww'
        domain.smooth = True

        e = -1.0
        w = 1.0
        h = w-e
        u = 2.0
        uh = u*h

        Br = Reflective_boundary(domain)     # Side walls
        Bd = Dirichlet_boundary([w, uh, 0])  # 2 m/s across the 3 m inlet: 



        
        domain.set_quantity('elevation', e)
        domain.set_quantity('stage', w)
        domain.set_quantity('xmomentum', uh)
        domain.set_boundary( {'left': Bd, 'right': Bd, 'top': Br, 'bottom': Br})

        for t in domain.evolve(yieldstep=1, finaltime = t_end):
            pass

        # Check that momentum is as it should be in the interior

        I = [[0, width/2.],
             [length/2., width/2.],
             [length, width/2.]]
        
        I = domain.geo_reference.get_absolute(I)
        f = file_function(swwfile,
                          quantities=['stage', 'xmomentum', 'ymomentum'],
                          interpolation_points=I,
                          verbose=False)

        for t in range(t_end+1):
            for i in range(3):
                #print i, t, f(t, i)            
                assert num.allclose(f(t, i), [w, uh, 0], atol=1.0e-6)
            

        # Check flows through the middle
        for i in range(5):
            x = length/2. + i*0.23674563 # Arbitrary
            cross_section = [[x, 0], [x, width]]

            cross_section = domain.geo_reference.get_absolute(cross_section)            
            time, Q = get_flow_through_cross_section(swwfile,
                                                     cross_section,
                                                     verbose=False)

            assert num.allclose(Q, uh*width)


            
    def test_get_energy_through_cross_section(self):
        """test_get_energy_through_cross_section(self):

        Test that the specific and total energy through a cross section can be
        correctly obtained from an sww file.
        
        This test creates a flat bed with a known flow through it and tests
        that the function correctly returns the expected energies.

        The specifics are
        u = 2 m/s
        h = 1 m
        w = 3 m (width of channel)

        q = u*h*w = 6 m^3/s
        Es = h + 0.5*v*v/g  # Specific energy head [m]
        Et = w + 0.5*v*v/g  # Total energy head [m]        


        This test uses georeferencing
        
        """

        import time, os
        from Scientific.IO.NetCDF import NetCDFFile

        # Setup
        from mesh_factory import rectangular

        # Create basic mesh (20m x 3m)
        width = 3
        length = 20
        t_end = 1
        points, vertices, boundary = rectangular(length, width,
                                                 length, width)

        # Create shallow water domain
        domain = Domain(points, vertices, boundary,
                        geo_reference = Geo_reference(56,308500,6189000))

        domain.default_order = 2
        domain.set_minimum_storable_height(0.01)

        domain.set_name('flowtest')
        swwfile = domain.get_name() + '.sww'

        domain.set_datadir('.')
        domain.format = 'sww'
        domain.smooth = True

        e = -1.0
        w = 1.0
        h = w-e
        u = 2.0
        uh = u*h

        Br = Reflective_boundary(domain)     # Side walls
        Bd = Dirichlet_boundary([w, uh, 0])  # 2 m/s across the 3 m inlet: 

        
        domain.set_quantity('elevation', e)
        domain.set_quantity('stage', w)
        domain.set_quantity('xmomentum', uh)
        domain.set_boundary( {'left': Bd, 'right': Bd, 'top': Br, 'bottom': Br})

        for t in domain.evolve(yieldstep=1, finaltime = t_end):
            pass

        # Check that momentum is as it should be in the interior

        I = [[0, width/2.],
             [length/2., width/2.],
             [length, width/2.]]
        
        I = domain.geo_reference.get_absolute(I)
        f = file_function(swwfile,
                          quantities=['stage', 'xmomentum', 'ymomentum'],
                          interpolation_points=I,
                          verbose=False)

        for t in range(t_end+1):
            for i in range(3):
                #print i, t, f(t, i)
                assert num.allclose(f(t, i), [w, uh, 0], atol=1.0e-6)
            

        # Check energies through the middle
        for i in range(5):
            x = length/2. + i*0.23674563 # Arbitrary
            cross_section = [[x, 0], [x, width]]

            cross_section = domain.geo_reference.get_absolute(cross_section)            
            
            time, Es = get_energy_through_cross_section(swwfile,
                                                       cross_section,
                                                       kind='specific',
                                                       verbose=False)
            assert num.allclose(Es, h + 0.5*u*u/g)
            
            time, Et = get_energy_through_cross_section(swwfile,
                                                        cross_section,
                                                        kind='total',
                                                        verbose=False)
            assert num.allclose(Et, w + 0.5*u*u/g)            

            
        
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

    def test_screen_catcher(self):
    
        stdout_orginal = sys.stdout
        stderr_orginal = sys.stderr
        
        output_dir = tempfile.mkdtemp('','output_test')
        #print output_dir
        start_screen_catcher(output_dir+sep,verbose=False)
        
        print 'hello screen catcher'
        print 'goodbye screen catcher'
        
        sys.stdout = stdout_orginal
        sys.stderr = stderr_orginal
        
        output_file = output_dir+sep+'screen_output.txt'
        assert access(output_file,F_OK)
        
        fid = open(output_file)
        file_line1 = fid.readline()
        file_line2 = fid.readline()
        
        #print 'file contents',file_line1, file_line2
        assert (file_line1 == 'hello screen catcher\n')
        assert (file_line2 == 'goodbye screen catcher\n')
        
 
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

            
    def test_copy_code_files(self):
        '''test that the copy_code_files() function is sane.'''

        def create_file(f):
            fd = open(f, 'w')
            fd.write('%s\n' % f)
            fd.close()

        # create working directories and test files
        work_dir = tempfile.mkdtemp()
        dst_dir = tempfile.mkdtemp(dir=work_dir)
        src_dir = tempfile.mkdtemp(dir=work_dir)

        f1 = 'file1'        
        filename1 = os.path.join(src_dir, f1)
        create_file(filename1)
        f2 = 'file2'        
        filename2 = os.path.join(src_dir, f2)
        create_file(filename2)
        f3 = 'file3'        
        filename3 = os.path.join(src_dir, f3)
        create_file(filename3)
        f4 = 'file4'        
        filename4 = os.path.join(src_dir, f4)
        create_file(filename4)
        f5 = 'file5'        
        filename5 = os.path.join(src_dir, f5)
        create_file(filename5)

        # exercise the copy function
        copy_code_files(dst_dir, filename1)
        copy_code_files(dst_dir, filename1, filename2)
        copy_code_files(dst_dir, (filename4, filename5, filename3))

        # test that files were actually copied
        self.failUnless(access(os.path.join(dst_dir, f1), F_OK))
        self.failUnless(access(os.path.join(dst_dir, f2), F_OK))
        self.failUnless(access(os.path.join(dst_dir, f3), F_OK))
        self.failUnless(access(os.path.join(dst_dir, f4), F_OK))
        self.failUnless(access(os.path.join(dst_dir, f5), F_OK))

        # clean up
        shutil.rmtree(work_dir)
            
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


