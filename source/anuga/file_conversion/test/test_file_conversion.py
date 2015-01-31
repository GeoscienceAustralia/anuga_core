
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.file_conversion.ferret2sww import ferret2sww
from anuga.utilities.numerical_tools import ensure_numeric, mean
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions\
     import Transmissive_boundary
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular
from anuga.file.sww import SWW_file
from anuga.file.sww import extent_sww
from anuga.file_conversion.urs2nc import lon_lat2grid
from anuga.config import netcdf_float, epsilon, g
from anuga.file.netcdf import NetCDFFile

from anuga.file_conversion.file_conversion import tsh2sww
from anuga.file_conversion.file_conversion import timefile2netcdf


from anuga.file.mux import WAVEHEIGHT_MUX_LABEL, EAST_VELOCITY_LABEL, \
                            NORTH_VELOCITY_LABEL

import sys
import unittest
import numpy as num
import copy
import os


class Test_File_Conversion(unittest.TestCase):
    """ A suite of tests to test file conversion functions.
        These tests are quite coarse-grained: converting a file
        and checking that its headers and some of its contents
        are correct.
    """
    verbose = False

    def set_verbose(self):
        Test_File_Conversion.verbose = True
        
    def setUp(self):
        import time
        
        self.verbose = Test_File_Conversion.verbose
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
            os.remove(self.test_MOST_file + ext)

        for file in ['timefile2netcdf_seconds.tms', 'timefile2netcdf.tms']:
            try:
                os.remove(file)
            except:
                pass             

    def test_ferret2sww1(self):
        """Test that georeferencing etc works when converting from
        ferret format (lat/lon) to sww format (UTM)
        """
        import os, sys

        #The test file has
        # LON = 150.66667, 150.83334, 151, 151.16667
        # LAT = -34.5, -34.33333, -34.16667, -34 ;
        # TIME = 0, 0.1, 0.6, 1.1, 1.6, 2.1 ;
        #
        # First value (index=0) in small_ha.nc is 0.3400644 cm,
        # Fourth value (index==3) is -6.50198 cm



        #Read
        from anuga.coordinate_transforms.redfearn import redfearn
        #fid = NetCDFFile(self.test_MOST_file)
        fid = NetCDFFile(self.test_MOST_file + '_ha.nc')
        first_value = fid.variables['HA'][:][0,0,0]
        fourth_value = fid.variables['HA'][:][0,0,3]
        fid.close()


        #Call conversion (with zero origin)
        #ferret2sww('small', verbose=False,
        #           origin = (56, 0, 0))
        ferret2sww(self.test_MOST_file, verbose=self.verbose,
                   origin = (56, 0, 0))

        #Work out the UTM coordinates for first point
        zone, e, n = redfearn(-34.5, 150.66667)
        #print zone, e, n

        #Read output file 'small.sww'
        #fid = NetCDFFile('small.sww')
        fid = NetCDFFile(self.test_MOST_file + '.sww')

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]

        #Check that first coordinate is correctly represented
        assert num.allclose(x[0], e)
        assert num.allclose(y[0], n)

        #Check first value
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        ymomentum = fid.variables['ymomentum'][:]

        #print ymomentum

        assert num.allclose(stage[0,0], first_value/100)  #Meters

        #Check fourth value
        assert num.allclose(stage[0,3], fourth_value/100)  #Meters

        fid.close()

        #Cleanup
        import os
        os.remove(self.test_MOST_file + '.sww')



    def test_ferret2sww_zscale(self):
        """Test that zscale workse
        """
        import os, sys

        #The test file has
        # LON = 150.66667, 150.83334, 151, 151.16667
        # LAT = -34.5, -34.33333, -34.16667, -34 ;
        # TIME = 0, 0.1, 0.6, 1.1, 1.6, 2.1 ;
        #
        # First value (index=0) in small_ha.nc is 0.3400644 cm,
        # Fourth value (index==3) is -6.50198 cm


        #Read
        from anuga.coordinate_transforms.redfearn import redfearn
        fid = NetCDFFile(self.test_MOST_file + '_ha.nc')
        first_value = fid.variables['HA'][:][0,0,0]
        fourth_value = fid.variables['HA'][:][0,0,3]
        fid.close()

        #Call conversion (with no scaling)
        ferret2sww(self.test_MOST_file, verbose=self.verbose,
                   origin = (56, 0, 0))

        #Work out the UTM coordinates for first point
        fid = NetCDFFile(self.test_MOST_file + '.sww')

        #Check values
        stage_1 = fid.variables['stage'][:]
        xmomentum_1 = fid.variables['xmomentum'][:]
        ymomentum_1 = fid.variables['ymomentum'][:]

        assert num.allclose(stage_1[0,0], first_value/100)  #Meters
        assert num.allclose(stage_1[0,3], fourth_value/100)  #Meters

        fid.close()

        #Call conversion (with scaling)
        ferret2sww(self.test_MOST_file,
                   zscale = 5,
                   verbose=self.verbose,
                   origin = (56, 0, 0))

        #Work out the UTM coordinates for first point
        fid = NetCDFFile(self.test_MOST_file + '.sww')

        #Check values
        stage_5 = fid.variables['stage'][:]
        xmomentum_5 = fid.variables['xmomentum'][:]
        ymomentum_5 = fid.variables['ymomentum'][:]
        elevation = fid.variables['elevation'][:]

        assert num.allclose(stage_5[0,0], 5*first_value/100)  #Meters
        assert num.allclose(stage_5[0,3], 5*fourth_value/100)  #Meters

        assert num.allclose(5*stage_1, stage_5)

        # Momentum will also be changed due to new depth

        depth_1 = stage_1-elevation
        depth_5 = stage_5-elevation


        for i in range(stage_1.shape[0]):
            for j in range(stage_1.shape[1]):            
                if depth_1[i,j] > epsilon:

                    scale = depth_5[i,j]/depth_1[i,j]
                    ref_xmomentum = xmomentum_1[i,j] * scale
                    ref_ymomentum = ymomentum_1[i,j] * scale
                    
                    #print i, scale, xmomentum_1[i,j], xmomentum_5[i,j]
                    
                    assert num.allclose(xmomentum_5[i,j], ref_xmomentum)
                    assert num.allclose(ymomentum_5[i,j], ref_ymomentum)
                    
        

        fid.close()


        #Cleanup
        import os
        os.remove(self.test_MOST_file + '.sww')



    def test_ferret2sww_2(self):
        """Test that georeferencing etc works when converting from
        ferret format (lat/lon) to sww format (UTM)
        """

        #The test file has
        # LON = 150.66667, 150.83334, 151, 151.16667
        # LAT = -34.5, -34.33333, -34.16667, -34 ;
        # TIME = 0, 0.1, 0.6, 1.1, 1.6, 2.1 ;
        #
        # First value (index=0) in small_ha.nc is 0.3400644 cm,
        # Fourth value (index==3) is -6.50198 cm


        from anuga.coordinate_transforms.redfearn import redfearn
        #fid = NetCDFFile('small_ha.nc')
        fid = NetCDFFile(self.test_MOST_file + '_ha.nc')

        #Pick a coordinate and a value

        time_index = 1
        lat_index = 0
        lon_index = 2

        test_value = fid.variables['HA'][:][time_index, lat_index, lon_index]
        test_time = fid.variables['TIME'][:][time_index]
        test_lat = fid.variables['LAT'][:][lat_index]
        test_lon = fid.variables['LON'][:][lon_index]

        linear_point_index = lat_index*4 + lon_index
        fid.close()

        #Call conversion (with zero origin)
        ferret2sww(self.test_MOST_file, verbose=self.verbose,
                   origin = (56, 0, 0))


        #Work out the UTM coordinates for test point
        zone, e, n = redfearn(test_lat, test_lon)

        #Read output file 'small.sww'
        fid = NetCDFFile(self.test_MOST_file + '.sww')

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]

        #Check that test coordinate is correctly represented
        assert num.allclose(x[linear_point_index], e)
        assert num.allclose(y[linear_point_index], n)

        #Check test value
        stage = fid.variables['stage'][:]

        assert num.allclose(stage[time_index, linear_point_index], test_value/100)

        fid.close()

        #Cleanup
        import os
        os.remove(self.test_MOST_file + '.sww')


    def test_ferret2sww_lat_long(self):
        # Test that min lat long works

        #The test file has
        # LON = 150.66667, 150.83334, 151, 151.16667
        # LAT = -34.5, -34.33333, -34.16667, -34 ;
        
        #Read
        from anuga.coordinate_transforms.redfearn import redfearn
        fid = NetCDFFile(self.test_MOST_file + '_ha.nc')
        first_value = fid.variables['HA'][:][0,0,0]
        fourth_value = fid.variables['HA'][:][0,0,3]
        fid.close()


        #Call conversion (with zero origin)
        #ferret2sww('small', verbose=self.verbose,
        #           origin = (56, 0, 0))
        ferret2sww(self.test_MOST_file, verbose=self.verbose,
                   origin = (56, 0, 0), minlat=-34.5, maxlat=-34)

        #Work out the UTM coordinates for first point
        zone, e, n = redfearn(-34.5, 150.66667)
        #print zone, e, n

        #Read output file 'small.sww'
        #fid = NetCDFFile('small.sww')
        fid = NetCDFFile(self.test_MOST_file + '.sww')

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        #Check that first coordinate is correctly represented
        assert 16 == len(x)

        fid.close()

        #Cleanup
        import os
        os.remove(self.test_MOST_file + '.sww')


    def test_ferret2sww_lat_longII(self):
        # Test that min lat long works

        #The test file has
        # LON = 150.66667, 150.83334, 151, 151.16667
        # LAT = -34.5, -34.33333, -34.16667, -34 ;
        
        #Read
        from anuga.coordinate_transforms.redfearn import redfearn
        fid = NetCDFFile(self.test_MOST_file + '_ha.nc')
        first_value = fid.variables['HA'][:][0,0,0]
        fourth_value = fid.variables['HA'][:][0,0,3]
        fid.close()


        #Call conversion (with zero origin)
        #ferret2sww('small', verbose=False,
        #           origin = (56, 0, 0))
        ferret2sww(self.test_MOST_file, verbose=False,
                   origin = (56, 0, 0), minlat=-34.4, maxlat=-34.2)

        #Work out the UTM coordinates for first point
        zone, e, n = redfearn(-34.5, 150.66667)
        #print zone, e, n

        #Read output file 'small.sww'
        #fid = NetCDFFile('small.sww')
        fid = NetCDFFile(self.test_MOST_file + '.sww')

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        #Check that first coordinate is correctly represented
        assert 12 == len(x)

        fid.close()

        #Cleanup
        import os
        os.remove(self.test_MOST_file + '.sww')
        
    def test_ferret2sww3(self):
        """Elevation included
        """

        #The test file has
        # LON = 150.66667, 150.83334, 151, 151.16667
        # LAT = -34.5, -34.33333, -34.16667, -34 ;
        # ELEVATION = [-1 -2 -3 -4
        #             -5 -6 -7 -8
        #              ...
        #              ...      -16]
        # where the top left corner is -1m,
        # and the ll corner is -13.0m
        #
        # First value (index=0) in small_ha.nc is 0.3400644 cm,
        # Fourth value (index==3) is -6.50198 cm

        from anuga.coordinate_transforms.redfearn import redfearn
        import os
        fid1 = NetCDFFile('test_ha.nc',netcdf_mode_w)
        fid2 = NetCDFFile('test_ua.nc',netcdf_mode_w)
        fid3 = NetCDFFile('test_va.nc',netcdf_mode_w)
        fid4 = NetCDFFile('test_e.nc',netcdf_mode_w)

        h1_list = [150.66667,150.83334,151.]
        h2_list = [-34.5,-34.33333]

        long_name = 'LON'
        lat_name = 'LAT'
        time_name = 'TIME'

        nx = 3
        ny = 2

        for fid in [fid1,fid2,fid3]:
            fid.createDimension(long_name,nx)
            fid.createVariable(long_name,netcdf_float,(long_name,))
            fid.variables[long_name].point_spacing='uneven'
            fid.variables[long_name].units='degrees_east'
            fid.variables[long_name][:] = h1_list

            fid.createDimension(lat_name,ny)
            fid.createVariable(lat_name,netcdf_float,(lat_name,))
            fid.variables[lat_name].point_spacing='uneven'
            fid.variables[lat_name].units='degrees_north'
            fid.variables[lat_name][:] = h2_list

            fid.createDimension(time_name,2)
            fid.createVariable(time_name,netcdf_float,(time_name,))
            fid.variables[time_name].point_spacing='uneven'
            fid.variables[time_name].units='seconds'
            fid.variables[time_name][:] = [0.,1.]
            #if fid == fid3: break


        for fid in [fid4]:
            fid.createDimension(long_name,nx)
            fid.createVariable(long_name,netcdf_float,(long_name,))
            fid.variables[long_name].point_spacing='uneven'
            fid.variables[long_name].units='degrees_east'
            fid.variables[long_name][:] = h1_list

            fid.createDimension(lat_name,ny)
            fid.createVariable(lat_name,netcdf_float,(lat_name,))
            fid.variables[lat_name].point_spacing='uneven'
            fid.variables[lat_name].units='degrees_north'
            fid.variables[lat_name][:] = h2_list

        name = {}
        name[fid1]='HA'
        name[fid2]='UA'
        name[fid3]='VA'
        name[fid4]='ELEVATION'

        units = {}
        units[fid1]='cm'
        units[fid2]='cm/s'
        units[fid3]='cm/s'
        units[fid4]='m'

        values = {}
        values[fid1]=[[[5., 10.,15.], [13.,18.,23.]],[[50.,100.,150.],[130.,180.,230.]]]
        values[fid2]=[[[1., 2.,3.], [4.,5.,6.]],[[7.,8.,9.],[10.,11.,12.]]]
        values[fid3]=[[[13., 12.,11.], [10.,9.,8.]],[[7.,6.,5.],[4.,3.,2.]]]
        values[fid4]=[[-3000,-3100,-3200],[-4000,-5000,-6000]]

        for fid in [fid1,fid2,fid3]:
          fid.createVariable(name[fid],netcdf_float,(time_name,lat_name,long_name))
          fid.variables[name[fid]].point_spacing='uneven'
          fid.variables[name[fid]].units=units[fid]
          fid.variables[name[fid]][:] = values[fid]
          fid.variables[name[fid]].missing_value = -99999999.
          #if fid == fid3: break

        for fid in [fid4]:
            fid.createVariable(name[fid],netcdf_float,(lat_name,long_name))
            fid.variables[name[fid]].point_spacing='uneven'
            fid.variables[name[fid]].units=units[fid]
            fid.variables[name[fid]][:] = values[fid]
            fid.variables[name[fid]].missing_value = -99999999.


        fid1.sync(); fid1.close()
        fid2.sync(); fid2.close()
        fid3.sync(); fid3.close()
        fid4.sync(); fid4.close()

        fid1 = NetCDFFile('test_ha.nc',netcdf_mode_r)
        fid2 = NetCDFFile('test_e.nc',netcdf_mode_r)
        fid3 = NetCDFFile('test_va.nc',netcdf_mode_r)


        first_amp = fid1.variables['HA'][:][0,0,0]
        third_amp = fid1.variables['HA'][:][0,0,2]
        first_elevation = fid2.variables['ELEVATION'][0,0]
        third_elevation= fid2.variables['ELEVATION'][:][0,2]
        first_speed = fid3.variables['VA'][0,0,0]
        third_speed = fid3.variables['VA'][:][0,0,2]

        fid1.close()
        fid2.close()
        fid3.close()

        #Call conversion (with zero origin)
        ferret2sww('test', verbose=self.verbose,
                   origin = (56, 0, 0), inverted_bathymetry=False)

        os.remove('test_va.nc')
        os.remove('test_ua.nc')
        os.remove('test_ha.nc')
        os.remove('test_e.nc')

        #Read output file 'test.sww'
        fid = NetCDFFile('test.sww')


        #Check first value
        elevation = fid.variables['elevation'][:]
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        ymomentum = fid.variables['ymomentum'][:]

        #print ymomentum
        first_height = first_amp/100 - first_elevation
        third_height = third_amp/100 - third_elevation
        first_momentum=first_speed*first_height/100
        third_momentum=third_speed*third_height/100

        assert num.allclose(ymomentum[0][0],first_momentum)  #Meters
        assert num.allclose(ymomentum[0][2],third_momentum)  #Meters

        fid.close()

        #Cleanup
        os.remove('test.sww')



    def test_ferret2sww4(self):
        """Like previous but with augmented variable names as
        in files produced by ferret as opposed to MOST
        """

        #The test file has
        # LON = 150.66667, 150.83334, 151, 151.16667
        # LAT = -34.5, -34.33333, -34.16667, -34 ;
        # ELEVATION = [-1 -2 -3 -4
        #             -5 -6 -7 -8
        #              ...
        #              ...      -16]
        # where the top left corner is -1m,
        # and the ll corner is -13.0m
        #
        # First value (index=0) in small_ha.nc is 0.3400644 cm,
        # Fourth value (index==3) is -6.50198 cm

        from anuga.coordinate_transforms.redfearn import redfearn
        import os
        fid1 = NetCDFFile('test_ha.nc',netcdf_mode_w)
        fid2 = NetCDFFile('test_ua.nc',netcdf_mode_w)
        fid3 = NetCDFFile('test_va.nc',netcdf_mode_w)
        fid4 = NetCDFFile('test_e.nc',netcdf_mode_w)

        h1_list = [150.66667,150.83334,151.]
        h2_list = [-34.5,-34.33333]

#        long_name = 'LON961_1261'
#        lat_name = 'LAT481_841'
#        time_name = 'TIME1'

        long_name = 'LON'
        lat_name = 'LAT'
        time_name = 'TIME'

        nx = 3
        ny = 2

        for fid in [fid1,fid2,fid3]:
            fid.createDimension(long_name,nx)
            fid.createVariable(long_name,netcdf_float,(long_name,))
            fid.variables[long_name].point_spacing='uneven'
            fid.variables[long_name].units='degrees_east'
            fid.variables[long_name][:] = h1_list

            fid.createDimension(lat_name,ny)
            fid.createVariable(lat_name,netcdf_float,(lat_name,))
            fid.variables[lat_name].point_spacing='uneven'
            fid.variables[lat_name].units='degrees_north'
            fid.variables[lat_name][:] = h2_list

            fid.createDimension(time_name,2)
            fid.createVariable(time_name,netcdf_float,(time_name,))
            fid.variables[time_name].point_spacing='uneven'
            fid.variables[time_name].units='seconds'
            fid.variables[time_name][:] = [0.,1.]
            #if fid == fid3: break


        for fid in [fid4]:
            fid.createDimension(long_name,nx)
            fid.createVariable(long_name,netcdf_float,(long_name,))
            fid.variables[long_name].point_spacing='uneven'
            fid.variables[long_name].units='degrees_east'
            fid.variables[long_name][:] = h1_list

            fid.createDimension(lat_name,ny)
            fid.createVariable(lat_name,netcdf_float,(lat_name,))
            fid.variables[lat_name].point_spacing='uneven'
            fid.variables[lat_name].units='degrees_north'
            fid.variables[lat_name][:] = h2_list

        name = {}
        name[fid1]='HA'
        name[fid2]='UA'
        name[fid3]='VA'
        name[fid4]='ELEVATION'

        units = {}
        units[fid1]='cm'
        units[fid2]='cm/s'
        units[fid3]='cm/s'
        units[fid4]='m'

        values = {}
        values[fid1]=[[[5., 10.,15.], [13.,18.,23.]],[[50.,100.,150.],[130.,180.,230.]]]
        values[fid2]=[[[1., 2.,3.], [4.,5.,6.]],[[7.,8.,9.],[10.,11.,12.]]]
        values[fid3]=[[[13., 12.,11.], [10.,9.,8.]],[[7.,6.,5.],[4.,3.,2.]]]
        values[fid4]=[[-3000,-3100,-3200],[-4000,-5000,-6000]]

        for fid in [fid1,fid2,fid3]:
          fid.createVariable(name[fid],netcdf_float,(time_name,lat_name,long_name))
          fid.variables[name[fid]].point_spacing='uneven'
          fid.variables[name[fid]].units=units[fid]
          fid.variables[name[fid]][:] = values[fid]
          fid.variables[name[fid]].missing_value = -99999999.
          #if fid == fid3: break

        for fid in [fid4]:
            fid.createVariable(name[fid],netcdf_float,(lat_name,long_name))
            fid.variables[name[fid]].point_spacing='uneven'
            fid.variables[name[fid]].units=units[fid]
            fid.variables[name[fid]][:] = values[fid]
            fid.variables[name[fid]].missing_value = -99999999.


        fid1.sync(); fid1.close()
        fid2.sync(); fid2.close()
        fid3.sync(); fid3.close()
        fid4.sync(); fid4.close()

        fid1 = NetCDFFile('test_ha.nc',netcdf_mode_r)
        fid2 = NetCDFFile('test_e.nc',netcdf_mode_r)
        fid3 = NetCDFFile('test_va.nc',netcdf_mode_r)


        first_amp = fid1.variables['HA'][:][0,0,0]
        third_amp = fid1.variables['HA'][:][0,0,2]
        first_elevation = fid2.variables['ELEVATION'][0,0]
        third_elevation= fid2.variables['ELEVATION'][:][0,2]
        first_speed = fid3.variables['VA'][0,0,0]
        third_speed = fid3.variables['VA'][:][0,0,2]

        fid1.close()
        fid2.close()
        fid3.close()

        #Call conversion (with zero origin)
        ferret2sww('test', verbose=self.verbose, origin = (56, 0, 0)
                   , inverted_bathymetry=False)

        os.remove('test_va.nc')
        os.remove('test_ua.nc')
        os.remove('test_ha.nc')
        os.remove('test_e.nc')

        #Read output file 'test.sww'
        fid = NetCDFFile('test.sww')


        #Check first value
        elevation = fid.variables['elevation'][:]
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        ymomentum = fid.variables['ymomentum'][:]

        #print ymomentum
        first_height = first_amp/100 - first_elevation
        third_height = third_amp/100 - third_elevation
        first_momentum=first_speed*first_height/100
        third_momentum=third_speed*third_height/100

        assert num.allclose(ymomentum[0][0],first_momentum)  #Meters
        assert num.allclose(ymomentum[0][2],third_momentum)  #Meters

        fid.close()

        #Cleanup
        os.remove('test.sww')




    def test_ferret2sww_nz_origin(self):
        from anuga.coordinate_transforms.redfearn import redfearn

        #Call conversion (with nonzero origin)
        ferret2sww(self.test_MOST_file, verbose=self.verbose,
                   origin = (56, 100000, 200000))


        #Work out the UTM coordinates for first point
        zone, e, n = redfearn(-34.5, 150.66667)

        #Read output file 'small.sww'
        #fid = NetCDFFile('small.sww', netcdf_mode_r)
        fid = NetCDFFile(self.test_MOST_file + '.sww')

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]

        #Check that first coordinate is correctly represented
        assert num.allclose(x[0], e-100000)
        assert num.allclose(y[0], n-200000)

        fid.close()

        #Cleanup
        os.remove(self.test_MOST_file + '.sww')


    def test_ferret2sww_lat_longII(self):
        # Test that min lat long works

        #The test file has
        # LON = 150.66667, 150.83334, 151, 151.16667
        # LAT = -34.5, -34.33333, -34.16667, -34 ;
        
        #Read
        from anuga.coordinate_transforms.redfearn import redfearn
        fid = NetCDFFile(self.test_MOST_file + '_ha.nc')
        first_value = fid.variables['HA'][:][0,0,0]
        fourth_value = fid.variables['HA'][:][0,0,3]
        fid.close()


        #Call conversion (with zero origin)
        #ferret2sww('small', verbose=self.verbose,
        #           origin = (56, 0, 0))
        try:
            ferret2sww(self.test_MOST_file, verbose=self.verbose,
                   origin = (56, 0, 0), minlat=-34.5, maxlat=-35)
        except AssertionError:
            pass
        else:
            self.assertTrue(0 ==1,  'Bad input did not throw exception error!')

    def test_sww_extent(self):
        """Not a test, rather a look at the sww format
        """

        import time, os


        self.domain.set_name('datatest' + str(id(self)))
        self.domain.format = 'sww'
        self.domain.smooth = True
        self.domain.reduction = mean
        self.domain.set_datadir('.')
        #self.domain.tight_slope_limiters = 1        


        sww = SWW_file(self.domain)
        sww.store_connectivity()
        sww.store_timestep()
        self.domain.time = 2.

        #Modify stage at second timestep
        stage = self.domain.quantities['stage'].vertex_values
        self.domain.set_quantity('stage', stage/2)

        sww.store_timestep()

        file_and_extension_name = self.domain.get_name() + ".sww"
        #print "file_and_extension_name",file_and_extension_name
        [xmin, xmax, ymin, ymax, stagemin, stagemax] = \
               extent_sww(file_and_extension_name )

        assert num.allclose(xmin, 0.0)
        assert num.allclose(xmax, 1.0)
        assert num.allclose(ymin, 0.0)
        assert num.allclose(ymax, 1.0)

        # FIXME (Ole): Revisit these numbers
        #assert num.allclose(stagemin, -0.85), 'stagemin=%.4f' %stagemin
        #assert num.allclose(stagemax, 0.15), 'stagemax=%.4f' %stagemax


        #Cleanup
        os.remove(sww.filename)

    def test_timefile2netcdf_seconds(self):

        pass
        #Write txt time file
        root = 'timefile2netcdf_seconds'

        file_text = root+'.txt'
        fid = open(file_text, 'w')
        fid.write(
"""0.0, 1.328223 0 0
0.1, 1.292912 0
0.2, 1.292912 0 0
""")
        fid.flush()
        fid.close()

        # Expecting error to be raised
        try:
            timefile2netcdf(file_text)
        except:
            pass

        # Should pass
        timefile2netcdf(file_text, time_as_seconds=True)
        
        #os.remove(root+'.tms')
        os.remove(root+'.txt')
        

    def test_timefile2netcdf(self):

        #Write txt time file
        root = 'timefile2netcdf'

        file_text = root+'.txt'
        fid = open(file_text, 'w')
        fid.write(
"""31/08/04 00:00:00, 1.328223 0 0
31/08/04 00:15:00, 1.292912 0 0
31/08/04 00:30:00, 1.292912 0 0
""")
        fid.flush()
        fid.close()

        # Expecting error to be raised
        try:
            timefile2netcdf(file_text,time_as_seconds=True)
        except:
            pass

        # Should pass
        timefile2netcdf(file_text)

        #os.remove(root+'.tms')
        os.remove(root+'.txt')

#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_File_Conversion,'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)

