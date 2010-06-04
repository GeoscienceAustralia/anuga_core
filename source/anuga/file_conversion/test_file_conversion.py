
from anuga.shallow_water.shallow_water_domain import Domain
from ferret2sww import ferret2sww
from anuga.utilities.numerical_tools import ensure_numeric, mean
from anuga.abstract_2d_finite_volumes.generic_boundary_conditions\
     import Transmissive_boundary
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular
from anuga.file.sww import SWW_file
from anuga.file.sww import extent_sww
from anuga.file_conversion.urs2nc import lon_lat2grid
from anuga.config import netcdf_float, epsilon, g
from Scientific.IO.NetCDF import NetCDFFile
from anuga.file_conversion.file_conversion import tsh2sww, \
                        pmesh_to_domain_instance


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

    def test_ferret2sww1(self):
        """Test that georeferencing etc works when converting from
        ferret format (lat/lon) to sww format (UTM)
        """
        from Scientific.IO.NetCDF import NetCDFFile
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
        from Scientific.IO.NetCDF import NetCDFFile
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
        from Scientific.IO.NetCDF import NetCDFFile

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
        from Scientific.IO.NetCDF import NetCDFFile

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
            fid.variables[long_name].assignValue(h1_list)

            fid.createDimension(lat_name,ny)
            fid.createVariable(lat_name,netcdf_float,(lat_name,))
            fid.variables[lat_name].point_spacing='uneven'
            fid.variables[lat_name].units='degrees_north'
            fid.variables[lat_name].assignValue(h2_list)

            fid.createDimension(time_name,2)
            fid.createVariable(time_name,netcdf_float,(time_name,))
            fid.variables[time_name].point_spacing='uneven'
            fid.variables[time_name].units='seconds'
            fid.variables[time_name].assignValue([0.,1.])
            #if fid == fid3: break


        for fid in [fid4]:
            fid.createDimension(long_name,nx)
            fid.createVariable(long_name,netcdf_float,(long_name,))
            fid.variables[long_name].point_spacing='uneven'
            fid.variables[long_name].units='degrees_east'
            fid.variables[long_name].assignValue(h1_list)

            fid.createDimension(lat_name,ny)
            fid.createVariable(lat_name,netcdf_float,(lat_name,))
            fid.variables[lat_name].point_spacing='uneven'
            fid.variables[lat_name].units='degrees_north'
            fid.variables[lat_name].assignValue(h2_list)

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
          fid.variables[name[fid]].assignValue(values[fid])
          fid.variables[name[fid]].missing_value = -99999999.
          #if fid == fid3: break

        for fid in [fid4]:
            fid.createVariable(name[fid],netcdf_float,(lat_name,long_name))
            fid.variables[name[fid]].point_spacing='uneven'
            fid.variables[name[fid]].units=units[fid]
            fid.variables[name[fid]].assignValue(values[fid])
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
        from Scientific.IO.NetCDF import NetCDFFile

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
            fid.variables[long_name].assignValue(h1_list)

            fid.createDimension(lat_name,ny)
            fid.createVariable(lat_name,netcdf_float,(lat_name,))
            fid.variables[lat_name].point_spacing='uneven'
            fid.variables[lat_name].units='degrees_north'
            fid.variables[lat_name].assignValue(h2_list)

            fid.createDimension(time_name,2)
            fid.createVariable(time_name,netcdf_float,(time_name,))
            fid.variables[time_name].point_spacing='uneven'
            fid.variables[time_name].units='seconds'
            fid.variables[time_name].assignValue([0.,1.])
            #if fid == fid3: break


        for fid in [fid4]:
            fid.createDimension(long_name,nx)
            fid.createVariable(long_name,netcdf_float,(long_name,))
            fid.variables[long_name].point_spacing='uneven'
            fid.variables[long_name].units='degrees_east'
            fid.variables[long_name].assignValue(h1_list)

            fid.createDimension(lat_name,ny)
            fid.createVariable(lat_name,netcdf_float,(lat_name,))
            fid.variables[lat_name].point_spacing='uneven'
            fid.variables[lat_name].units='degrees_north'
            fid.variables[lat_name].assignValue(h2_list)

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
          fid.variables[name[fid]].assignValue(values[fid])
          fid.variables[name[fid]].missing_value = -99999999.
          #if fid == fid3: break

        for fid in [fid4]:
            fid.createVariable(name[fid],netcdf_float,(lat_name,long_name))
            fid.variables[name[fid]].point_spacing='uneven'
            fid.variables[name[fid]].units=units[fid]
            fid.variables[name[fid]].assignValue(values[fid])
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
        from Scientific.IO.NetCDF import NetCDFFile
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
            self.failUnless(0 ==1,  'Bad input did not throw exception error!')

    def test_sww_extent(self):
        """Not a test, rather a look at the sww format
        """

        import time, os
        from Scientific.IO.NetCDF import NetCDFFile

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



#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_File_Conversion, 'test_sww')
    
    # FIXME(Ole): When Ross has implemented logging, we can 
    # probably get rid of all this:
    if len(sys.argv) > 1 and sys.argv[1][0].upper() == 'V':
        Test_File_Conversion.verbose=True
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


    def test_asc_csiro2sww(self):
        import tempfile

        bath_dir = tempfile.mkdtemp()
        bath_dir_filename = bath_dir + os.sep +'ba19940524.000'
        #bath_dir = 'bath_data_manager_test'
        #print "os.getcwd( )",os.getcwd( )
        elevation_dir =  tempfile.mkdtemp()
        #elevation_dir = 'elev_expanded'
        elevation_dir_filename1 = elevation_dir + os.sep +'el19940524.000'
        elevation_dir_filename2 = elevation_dir + os.sep +'el19940524.001'

        fid = open(bath_dir_filename, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    9000.000    -1000.000    3000.0
   -1000.000    9000.000  -1000.000
""")
        fid.close()

        fid = open(elevation_dir_filename1, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    9000.000    0.000    3000.0
     0.000     9000.000     0.000
""")
        fid.close()

        fid = open(elevation_dir_filename2, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    9000.000    4000.000    4000.0
    4000.000    9000.000    4000.000
""")
        fid.close()

        ucur_dir =  tempfile.mkdtemp()
        ucur_dir_filename1 = ucur_dir + os.sep +'uc19940524.000'
        ucur_dir_filename2 = ucur_dir + os.sep +'uc19940524.001'

        fid = open(ucur_dir_filename1, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    90.000    60.000    30.0
    10.000    10.000    10.000
""")
        fid.close()
        fid = open(ucur_dir_filename2, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    90.000    60.000    30.0
    10.000    10.000    10.000
""")
        fid.close()

        vcur_dir =  tempfile.mkdtemp()
        vcur_dir_filename1 = vcur_dir + os.sep +'vc19940524.000'
        vcur_dir_filename2 = vcur_dir + os.sep +'vc19940524.001'

        fid = open(vcur_dir_filename1, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    90.000    60.000    30.0
    10.000    10.000    10.000
""")
        fid.close()
        fid = open(vcur_dir_filename2, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    90.000    60.000    30.0
    10.000    10.000    10.000
""")
        fid.close()

        sww_file = 'a_test.sww'
        asc_csiro2sww(bath_dir,elevation_dir, ucur_dir, vcur_dir, sww_file,
                      verbose=self.verbose)

        # check the sww file

        fid = NetCDFFile(sww_file, netcdf_mode_r)    # Open existing file for read
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        z = fid.variables['elevation'][:]
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        geo_ref = Geo_reference(NetCDFObject=fid)
        #print "geo_ref",geo_ref
        x_ref = geo_ref.get_xllcorner()
        y_ref = geo_ref.get_yllcorner()
        self.failUnless(geo_ref.get_zone() == 55,  'Failed')
        assert num.allclose(x_ref, 587798.418) # (-38, 148)
        assert num.allclose(y_ref, 5793123.477)# (-38, 148.5)

        #Zone:   55
        #Easting:  588095.674  Northing: 5821451.722
        #Latitude:   -37  45 ' 0.00000 ''  Longitude: 148 0 ' 0.00000 ''
        assert num.allclose((x[0],y[0]), (588095.674 - x_ref, 5821451.722 - y_ref))

        #Zone:   55
        #Easting:  632145.632  Northing: 5820863.269
        #Latitude:   -37  45 ' 0.00000 ''  Longitude: 148  30 ' 0.00000 ''
        assert num.allclose((x[2],y[2]), (632145.632 - x_ref, 5820863.269 - y_ref))

        #Zone:   55
        #Easting:  609748.788  Northing: 5793447.860
        #Latitude:   -38  0 ' 0.00000 ''  Longitude: 148  15 ' 0.00000 ''
        assert num.allclose((x[4],y[4]), (609748.788  - x_ref, 5793447.86 - y_ref))

        assert num.allclose(z[0],9000.0 )
        assert num.allclose(stage[0][1],0.0 )

        #(4000+1000)*60
        assert num.allclose(xmomentum[1][1],300000.0 )


        fid.close()

        #tidy up
        os.remove(bath_dir_filename)
        os.rmdir(bath_dir)

        os.remove(elevation_dir_filename1)
        os.remove(elevation_dir_filename2)
        os.rmdir(elevation_dir)

        os.remove(ucur_dir_filename1)
        os.remove(ucur_dir_filename2)
        os.rmdir(ucur_dir)

        os.remove(vcur_dir_filename1)
        os.remove(vcur_dir_filename2)
        os.rmdir(vcur_dir)


        # remove sww file
        os.remove(sww_file)

    def test_asc_csiro2sww2(self):
        import tempfile

        bath_dir = tempfile.mkdtemp()
        bath_dir_filename = bath_dir + os.sep +'ba19940524.000'
        #bath_dir = 'bath_data_manager_test'
        #print "os.getcwd( )",os.getcwd( )
        elevation_dir =  tempfile.mkdtemp()
        #elevation_dir = 'elev_expanded'
        elevation_dir_filename1 = elevation_dir + os.sep +'el19940524.000'
        elevation_dir_filename2 = elevation_dir + os.sep +'el19940524.001'

        fid = open(bath_dir_filename, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    9000.000    -1000.000    3000.0
   -1000.000    9000.000  -1000.000
""")
        fid.close()

        fid = open(elevation_dir_filename1, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    9000.000    0.000    3000.0
     0.000     -9999.000     -9999.000
""")
        fid.close()

        fid = open(elevation_dir_filename2, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    9000.000    4000.000    4000.0
    4000.000    9000.000    4000.000
""")
        fid.close()

        ucur_dir =  tempfile.mkdtemp()
        ucur_dir_filename1 = ucur_dir + os.sep +'uc19940524.000'
        ucur_dir_filename2 = ucur_dir + os.sep +'uc19940524.001'

        fid = open(ucur_dir_filename1, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    90.000    60.000    30.0
    10.000    10.000    10.000
""")
        fid.close()
        fid = open(ucur_dir_filename2, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    90.000    60.000    30.0
    10.000    10.000    10.000
""")
        fid.close()

        vcur_dir =  tempfile.mkdtemp()
        vcur_dir_filename1 = vcur_dir + os.sep +'vc19940524.000'
        vcur_dir_filename2 = vcur_dir + os.sep +'vc19940524.001'

        fid = open(vcur_dir_filename1, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    90.000    60.000    30.0
    10.000    10.000    10.000
""")
        fid.close()
        fid = open(vcur_dir_filename2, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    90.000    60.000    30.0
    10.000    10.000    10.000
""")
        fid.close()

        try:
            asc_csiro2sww(bath_dir,elevation_dir, ucur_dir,
                          vcur_dir, sww_file,
                      verbose=self.verbose)
        except:
            #tidy up
            os.remove(bath_dir_filename)
            os.rmdir(bath_dir)

            os.remove(elevation_dir_filename1)
            os.remove(elevation_dir_filename2)
            os.rmdir(elevation_dir)

            os.remove(ucur_dir_filename1)
            os.remove(ucur_dir_filename2)
            os.rmdir(ucur_dir)

            os.remove(vcur_dir_filename1)
            os.remove(vcur_dir_filename2)
            os.rmdir(vcur_dir)
        else:
            #tidy up
            os.remove(bath_dir_filename)
            os.rmdir(bath_dir)

            os.remove(elevation_dir_filename1)
            os.remove(elevation_dir_filename2)
            os.rmdir(elevation_dir)
            raise 'Should raise exception'

            os.remove(ucur_dir_filename1)
            os.remove(ucur_dir_filename2)
            os.rmdir(ucur_dir)

            os.remove(vcur_dir_filename1)
            os.remove(vcur_dir_filename2)
            os.rmdir(vcur_dir)



    def test_asc_csiro2sww3(self):
        import tempfile

        bath_dir = tempfile.mkdtemp()
        bath_dir_filename = bath_dir + os.sep +'ba19940524.000'
        #bath_dir = 'bath_data_manager_test'
        #print "os.getcwd( )",os.getcwd( )
        elevation_dir =  tempfile.mkdtemp()
        #elevation_dir = 'elev_expanded'
        elevation_dir_filename1 = elevation_dir + os.sep +'el19940524.000'
        elevation_dir_filename2 = elevation_dir + os.sep +'el19940524.001'

        fid = open(bath_dir_filename, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    9000.000    -1000.000    3000.0
   -1000.000    9000.000  -1000.000
""")
        fid.close()

        fid = open(elevation_dir_filename1, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    9000.000    0.000    3000.0
     0.000     -9999.000     -9999.000
""")
        fid.close()

        fid = open(elevation_dir_filename2, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    9000.000    4000.000    4000.0
    4000.000    9000.000    4000.000
""")
        fid.close()

        ucur_dir =  tempfile.mkdtemp()
        ucur_dir_filename1 = ucur_dir + os.sep +'uc19940524.000'
        ucur_dir_filename2 = ucur_dir + os.sep +'uc19940524.001'

        fid = open(ucur_dir_filename1, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    90.000    60.000    30.0
    10.000    10.000    10.000
""")
        fid.close()
        fid = open(ucur_dir_filename2, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    90.000    60.000    30.0
    10.000    10.000    10.000
""")
        fid.close()

        vcur_dir =  tempfile.mkdtemp()
        vcur_dir_filename1 = vcur_dir + os.sep +'vc19940524.000'
        vcur_dir_filename2 = vcur_dir + os.sep +'vc19940524.001'

        fid = open(vcur_dir_filename1, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    90.000    60.000    30.0
    10.000    10.000    10.000
""")
        fid.close()
        fid = open(vcur_dir_filename2, 'w')
        fid.write(""" ncols             3
 nrows             2
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
    90.000    60.000    30.0
    10.000    10.000    10.000
""")
        fid.close()

        sww_file = 'a_test.sww'
        asc_csiro2sww(bath_dir,elevation_dir, ucur_dir, vcur_dir,
                      sww_file, fail_on_NaN = False, elevation_NaN_filler = 0,
                      mean_stage = 100,
                      verbose=self.verbose)

        # check the sww file

        fid = NetCDFFile(sww_file, netcdf_mode_r)    # Open existing file for read
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        z = fid.variables['elevation'][:]
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        geo_ref = Geo_reference(NetCDFObject=fid)
        #print "geo_ref",geo_ref
        x_ref = geo_ref.get_xllcorner()
        y_ref = geo_ref.get_yllcorner()
        self.failUnless(geo_ref.get_zone() == 55,  'Failed')
        assert num.allclose(x_ref, 587798.418) # (-38, 148)
        assert num.allclose(y_ref, 5793123.477)# (-38, 148.5)

        #Zone:   55
        #Easting:  588095.674  Northing: 5821451.722
        #Latitude:   -37  45 ' 0.00000 ''  Longitude: 148 0 ' 0.00000 ''
        assert num.allclose((x[0],y[0]), (588095.674 - x_ref, 5821451.722 - y_ref))

        #Zone:   55
        #Easting:  632145.632  Northing: 5820863.269
        #Latitude:   -37  45 ' 0.00000 ''  Longitude: 148  30 ' 0.00000 ''
        assert num.allclose((x[2],y[2]), (632145.632 - x_ref, 5820863.269 - y_ref))

        #Zone:   55
        #Easting:  609748.788  Northing: 5793447.860
        #Latitude:   -38  0 ' 0.00000 ''  Longitude: 148  15 ' 0.00000 ''
        assert num.allclose((x[4],y[4]), (609748.788  - x_ref, 5793447.86 - y_ref))

        assert num.allclose(z[0],9000.0 )
        assert num.allclose(stage[0][4],100.0 )
        assert num.allclose(stage[0][5],100.0 )

        #(100.0 - 9000)*10
        assert num.allclose(xmomentum[0][4], -89000.0 )

        #(100.0 - -1000.000)*10
        assert num.allclose(xmomentum[0][5], 11000.0 )

        fid.close()

        #tidy up
        os.remove(bath_dir_filename)
        os.rmdir(bath_dir)

        os.remove(elevation_dir_filename1)
        os.remove(elevation_dir_filename2)
        os.rmdir(elevation_dir)

        os.remove(ucur_dir_filename1)
        os.remove(ucur_dir_filename2)
        os.rmdir(ucur_dir)

        os.remove(vcur_dir_filename1)
        os.remove(vcur_dir_filename2)
        os.rmdir(vcur_dir)

        # remove sww file
        os.remove(sww_file)


    def test_asc_csiro2sww4(self):
        """
        Test specifying the extent
        """

        import tempfile

        bath_dir = tempfile.mkdtemp()
        bath_dir_filename = bath_dir + os.sep +'ba19940524.000'
        #bath_dir = 'bath_data_manager_test'
        #print "os.getcwd( )",os.getcwd( )
        elevation_dir =  tempfile.mkdtemp()
        #elevation_dir = 'elev_expanded'
        elevation_dir_filename1 = elevation_dir + os.sep +'el19940524.000'
        elevation_dir_filename2 = elevation_dir + os.sep +'el19940524.001'

        fid = open(bath_dir_filename, 'w')
        fid.write(""" ncols             4
 nrows             4
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
   -9000.000    -1000.000   -3000.0 -2000.000
   -1000.000    9000.000  -1000.000 -3000.000
   -4000.000    6000.000   2000.000 -5000.000
   -9000.000    -1000.000   -3000.0 -2000.000
""")
        fid.close()

        fid = open(elevation_dir_filename1, 'w')
        fid.write(""" ncols             4
 nrows             4
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
   -900.000    -100.000   -300.0 -200.000
   -100.000    900.000  -100.000 -300.000
   -400.000    600.000   200.000 -500.000
   -900.000    -100.000   -300.0 -200.000
""")
        fid.close()

        fid = open(elevation_dir_filename2, 'w')
        fid.write(""" ncols             4
 nrows             4
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
   -990.000    -110.000   -330.0 -220.000
   -110.000    990.000  -110.000 -330.000
   -440.000    660.000   220.000 -550.000
   -990.000    -110.000   -330.0 -220.000
""")
        fid.close()

        ucur_dir =  tempfile.mkdtemp()
        ucur_dir_filename1 = ucur_dir + os.sep +'uc19940524.000'
        ucur_dir_filename2 = ucur_dir + os.sep +'uc19940524.001'

        fid = open(ucur_dir_filename1, 'w')
        fid.write(""" ncols             4
 nrows             4
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
   -90.000    -10.000   -30.0 -20.000
   -10.000    90.000  -10.000 -30.000
   -40.000    60.000   20.000 -50.000
   -90.000    -10.000   -30.0 -20.000
""")
        fid.close()
        fid = open(ucur_dir_filename2, 'w')
        fid.write(""" ncols             4
 nrows             4
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
   -90.000    -10.000   -30.0 -20.000
   -10.000    99.000  -11.000 -30.000
   -40.000    66.000   22.000 -50.000
   -90.000    -10.000   -30.0 -20.000
""")
        fid.close()

        vcur_dir =  tempfile.mkdtemp()
        vcur_dir_filename1 = vcur_dir + os.sep +'vc19940524.000'
        vcur_dir_filename2 = vcur_dir + os.sep +'vc19940524.001'

        fid = open(vcur_dir_filename1, 'w')
        fid.write(""" ncols             4
 nrows             4
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
   -90.000    -10.000   -30.0 -20.000
   -10.000    80.000  -20.000 -30.000
   -40.000    50.000   10.000 -50.000
   -90.000    -10.000   -30.0 -20.000
""")
        fid.close()
        fid = open(vcur_dir_filename2, 'w')
        fid.write(""" ncols             4
 nrows             4
 xllcorner    148.00000
 yllcorner    -38.00000
 cellsize       0.25
 nodata_value   -9999.0
   -90.000    -10.000   -30.0 -20.000
   -10.000    88.000  -22.000 -30.000
   -40.000    55.000   11.000 -50.000
   -90.000    -10.000   -30.0 -20.000
""")
        fid.close()

        sww_file = tempfile.mktemp(".sww")
        #sww_file = 'a_test.sww'
        asc_csiro2sww(bath_dir,elevation_dir, ucur_dir, vcur_dir,
                      sww_file, fail_on_NaN = False, elevation_NaN_filler = 0,
                      mean_stage = 100,
                       minlat = -37.6, maxlat = -37.6,
                  minlon = 148.3, maxlon = 148.3,
                      verbose=self.verbose
                      #,verbose = True
                      )

        # check the sww file

        fid = NetCDFFile(sww_file, netcdf_mode_r)    # Open existing file for read
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        z = fid.variables['elevation'][:]
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        ymomentum = fid.variables['ymomentum'][:]
        geo_ref = Geo_reference(NetCDFObject=fid)
        #print "geo_ref",geo_ref
        x_ref = geo_ref.get_xllcorner()
        y_ref = geo_ref.get_yllcorner()
        self.failUnless(geo_ref.get_zone() == 55,  'Failed')

        assert num.allclose(fid.starttime, 0.0) # (-37.45, 148.25)
        assert num.allclose(x_ref, 610120.388) # (-37.45, 148.25)
        assert num.allclose(y_ref,  5820863.269 )# (-37.45, 148.5)

        #Easting:  632145.632  Northing: 5820863.269
        #Latitude:   -37 45 ' 0.00000 ''  Longitude: 148  30 ' 0.00000 ''

        #print "x",x
        #print "y",y
        self.failUnless(len(x) == 4,'failed') # 2*2
        self.failUnless(len(x) == 4,'failed') # 2*2

        #Zone:   55
        #Easting:  632145.632  Northing: 5820863.269
        #Latitude:   -37 45 ' 0.00000 ''  Longitude: 148  30 ' 0.00000 ''
        # magic number - y is close enough for me.
        assert num.allclose(x[3], 632145.63 - x_ref)
        assert num.allclose(y[3], 5820863.269  - y_ref + 5.22155314684e-005)

        assert num.allclose(z[0],9000.0 ) #z is elevation info
        #print "z",z
        # 2 time steps, 4 points
        self.failUnless(xmomentum.shape == (2,4), 'failed')
        self.failUnless(ymomentum.shape == (2,4), 'failed')

        #(100.0 - -1000.000)*10
        #assert num.allclose(xmomentum[0][5], 11000.0 )

        fid.close()

        # is the sww file readable?
        #Lets see if we can convert it to a dem!
        # if you uncomment, remember to delete the file
        #print "sww_file",sww_file
        #dem_file = tempfile.mktemp(".dem")
        domain = load_sww_as_domain(sww_file) ###, dem_file)
        domain.check_integrity()

        #tidy up
        os.remove(bath_dir_filename)
        os.rmdir(bath_dir)

        os.remove(elevation_dir_filename1)
        os.remove(elevation_dir_filename2)
        os.rmdir(elevation_dir)

        os.remove(ucur_dir_filename1)
        os.remove(ucur_dir_filename2)
        os.rmdir(ucur_dir)

        os.remove(vcur_dir_filename1)
        os.remove(vcur_dir_filename2)
        os.rmdir(vcur_dir)

        # remove sww file
        os.remove(sww_file)



    def test_get_min_max_indexes(self):
        latitudes = [3,2,1,0]
        longitudes = [0,10,20,30]

        # k - lat
        # l - lon
        kmin, kmax, lmin, lmax = get_min_max_indices(
            latitudes,longitudes,
            -10,4,-10,31)

        #print "kmin",kmin;print "kmax",kmax
        #print "lmin",lmin;print "lmax",lmax
        latitudes_new = latitudes[kmin:kmax]
        longitudes_news = longitudes[lmin:lmax]
        #print "latitudes_new", latitudes_new
        #print "longitudes_news",longitudes_news
        self.failUnless(latitudes == latitudes_new and \
                        longitudes == longitudes_news,
                         'failed')

        ## 2nd test
        kmin, kmax, lmin, lmax = get_min_max_indices(
            latitudes,longitudes,
            0.5,2.5,5,25)
        #print "kmin",kmin;print "kmax",kmax
        #print "lmin",lmin;print "lmax",lmax
        latitudes_new = latitudes[kmin:kmax]
        longitudes_news = longitudes[lmin:lmax]
        #print "latitudes_new", latitudes_new
        #print "longitudes_news",longitudes_news

        self.failUnless(latitudes == latitudes_new and \
                        longitudes == longitudes_news,
                         'failed')

        ## 3rd test
        kmin, kmax, lmin, lmax = get_min_max_indices(\
            latitudes,
            longitudes,
            1.1,1.9,12,17)
        #print "kmin",kmin;print "kmax",kmax
        #print "lmin",lmin;print "lmax",lmax
        latitudes_new = latitudes[kmin:kmax]
        longitudes_news = longitudes[lmin:lmax]
        #print "latitudes_new", latitudes_new
        #print "longitudes_news",longitudes_news

        self.failUnless(latitudes_new == [2, 1] and \
                        longitudes_news == [10, 20],
                         'failed')


        ## 4th test
        kmin, kmax, lmin, lmax = get_min_max_indices(
            latitudes,longitudes,
                                                      -0.1,1.9,-2,17)
        #print "kmin",kmin;print "kmax",kmax
        #print "lmin",lmin;print "lmax",lmax
        latitudes_new = latitudes[kmin:kmax]
        longitudes_news = longitudes[lmin:lmax]
        #print "latitudes_new", latitudes_new
        #print "longitudes_news",longitudes_news

        self.failUnless(latitudes_new == [2, 1, 0] and \
                        longitudes_news == [0, 10, 20],
                         'failed')
        ## 5th test
        kmin, kmax, lmin, lmax = get_min_max_indices(
            latitudes,longitudes,
            0.1,1.9,2,17)
        #print "kmin",kmin;print "kmax",kmax
        #print "lmin",lmin;print "lmax",lmax
        latitudes_new = latitudes[kmin:kmax]
        longitudes_news = longitudes[lmin:lmax]
        #print "latitudes_new", latitudes_new
        #print "longitudes_news",longitudes_news

        self.failUnless(latitudes_new == [2, 1, 0] and \
                        longitudes_news == [0, 10, 20],
                         'failed')

        ## 6th test

        kmin, kmax, lmin, lmax = get_min_max_indices(
            latitudes,longitudes,
            1.5,4,18,32)
        #print "kmin",kmin;print "kmax",kmax
        #print "lmin",lmin;print "lmax",lmax
        latitudes_new = latitudes[kmin:kmax]
        longitudes_news = longitudes[lmin:lmax]
        #print "latitudes_new", latitudes_new
        #print "longitudes_news",longitudes_news

        self.failUnless(latitudes_new == [3, 2, 1] and \
                        longitudes_news == [10, 20, 30],
                         'failed')


        ## 7th test
        m2d = num.array([[0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15]], num.int)    #array default#
        kmin, kmax, lmin, lmax = get_min_max_indices(
            latitudes,longitudes,
            1.5,1.5,15,15)
        #print "kmin",kmin;print "kmax",kmax
        #print "lmin",lmin;print "lmax",lmax
        latitudes_new = latitudes[kmin:kmax]
        longitudes_news = longitudes[lmin:lmax]
        m2d = m2d[kmin:kmax,lmin:lmax]
        #print "m2d", m2d
        #print "latitudes_new", latitudes_new
        #print "longitudes_news",longitudes_news

        self.failUnless(num.alltrue(latitudes_new == [2, 1]) and 
                        num.alltrue(longitudes_news == [10, 20]),
                        'failed')

        self.failUnless(num.alltrue(m2d == [[5,6],[9,10]]), 'failed')

    def test_get_min_max_indexes_lat_ascending(self):
        latitudes = [0,1,2,3]
        longitudes = [0,10,20,30]

        # k - lat
        # l - lon
        kmin, kmax, lmin, lmax = get_min_max_indices(
            latitudes,longitudes,
            -10,4,-10,31)

        #print "kmin",kmin;print "kmax",kmax
        #print "lmin",lmin;print "lmax",lmax
        latitudes_new = latitudes[kmin:kmax]
        longitudes_news = longitudes[lmin:lmax]
        #print "latitudes_new", latitudes_new
        #print "longitudes_news",longitudes_news
        self.failUnless(latitudes == latitudes_new and \
                        longitudes == longitudes_news,
                         'failed')

        ## 3rd test
        kmin, kmax, lmin, lmax = get_min_max_indices(\
            latitudes,
            longitudes,
            1.1,1.9,12,17)
        #print "kmin",kmin;print "kmax",kmax
        #print "lmin",lmin;print "lmax",lmax
        latitudes_new = latitudes[kmin:kmax]
        longitudes_news = longitudes[lmin:lmax]
        #print "latitudes_new", latitudes_new
        #print "longitudes_news",longitudes_news

        self.failUnless(latitudes_new == [1, 2] and \
                        longitudes_news == [10, 20],
                         'failed')

    def test_get_min_max_indexes2(self):
        latitudes = [-30,-35,-40,-45]
        longitudes = [148,149,150,151]

        m2d = num.array([[0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15]], num.int)    #array default#

        # k - lat
        # l - lon
        kmin, kmax, lmin, lmax = get_min_max_indices(
            latitudes,longitudes,
            -37,-27,147,149.5)

        #print "kmin",kmin;print "kmax",kmax
        #print "lmin",lmin;print "lmax",lmax
        #print "m2d", m2d
        #print "latitudes", latitudes
        #print "longitudes",longitudes
        #print "latitudes[kmax]", latitudes[kmax]
        latitudes_new = latitudes[kmin:kmax]
        longitudes_new = longitudes[lmin:lmax]
        m2d = m2d[kmin:kmax,lmin:lmax]
        #print "m2d", m2d
        #print "latitudes_new", latitudes_new
        #print "longitudes_new",longitudes_new

        self.failUnless(latitudes_new == [-30, -35, -40] and
                        longitudes_new == [148, 149,150],
                        'failed')
        self.failUnless(num.alltrue(m2d == [[0,1,2],[4,5,6],[8,9,10]]),
                        'failed')

    def test_get_min_max_indexes3(self):
        latitudes = [-30,-35,-40,-45,-50,-55,-60]
        longitudes = [148,149,150,151]

        # k - lat
        # l - lon
        kmin, kmax, lmin, lmax = get_min_max_indices(
            latitudes,longitudes,
            -43,-37,148.5,149.5)


        #print "kmin",kmin;print "kmax",kmax
        #print "lmin",lmin;print "lmax",lmax
        #print "latitudes", latitudes
        #print "longitudes",longitudes
        latitudes_new = latitudes[kmin:kmax]
        longitudes_news = longitudes[lmin:lmax]
        #print "latitudes_new", latitudes_new
        #print "longitudes_news",longitudes_news

        self.failUnless(latitudes_new == [-35, -40, -45] and
                        longitudes_news == [148, 149,150],
                         'failed')

    def test_get_min_max_indexes4(self):
        latitudes = [-30,-35,-40,-45,-50,-55,-60]
        longitudes = [148,149,150,151]

        # k - lat
        # l - lon
        kmin, kmax, lmin, lmax = get_min_max_indices(
            latitudes,longitudes)

        #print "kmin",kmin;print "kmax",kmax
        #print "lmin",lmin;print "lmax",lmax
        #print "latitudes", latitudes
        #print "longitudes",longitudes
        latitudes_new = latitudes[kmin:kmax]
        longitudes_news = longitudes[lmin:lmax]
        #print "latitudes_new", latitudes_new
        #print "longitudes_news",longitudes_news

        self.failUnless(latitudes_new == latitudes  and \
                        longitudes_news == longitudes,
                         'failed')

    def test_tsh2sww(self):
        import os
        import tempfile

        tsh_file = tempfile.mktemp(".tsh")
        file = open(tsh_file,"w")
        file.write("4 3 # <vertex #> <x> <y> [attributes]\n \
0 0.0 0.0 0.0 0.0 0.01 \n \
1 1.0 0.0 10.0 10.0 0.02  \n \
2 0.0 1.0 0.0 10.0 0.03  \n \
3 0.5 0.25 8.0 12.0 0.04  \n \
# Vert att title  \n \
elevation  \n \
stage  \n \
friction  \n \
2 # <triangle #> [<vertex #>] [<neigbouring triangle #>]  \n\
0 0 3 2 -1  -1  1 dsg\n\
1 0 1 3 -1  0 -1   ole nielsen\n\
4 # <segment #> <vertex #>  <vertex #> [boundary tag] \n\
0 1 0 2 \n\
1 0 2 3 \n\
2 2 3 \n\
3 3 1 1 \n\
3 0 # <x> <y> [attributes] ...Mesh Vertices... \n \
0 216.0 -86.0 \n \
1 160.0 -167.0 \n \
2 114.0 -91.0 \n \
3 # <vertex #>  <vertex #> [boundary tag] ...Mesh Segments... \n \
0 0 1 0 \n \
1 1 2 0 \n \
2 2 0 0 \n \
0 # <x> <y> ...Mesh Holes... \n \
0 # <x> <y> <attribute>...Mesh Regions... \n \
0 # <x> <y> <attribute>...Mesh Regions, area... \n\
#Geo reference \n \
56 \n \
140 \n \
120 \n")
        file.close()

        #sww_file = tempfile.mktemp(".sww")
        #print "sww_file",sww_file
        #print "sww_file",tsh_file
        tsh2sww(tsh_file,
                verbose=self.verbose)

        os.remove(tsh_file)
        os.remove(tsh_file[:-4] + '.sww')
        
          
    def test_urs2sww_test_fail(self):
        points_num = -100
        time_step_count = 45
        time_step = -7
        file_handle, base_name = tempfile.mkstemp("")        
        os.close(file_handle)
        os.remove(base_name)
        
        files = []
        quantities = ['HA','UA','VA']
        
        mux_names = [WAVEHEIGHT_MUX_LABEL,
                     EAST_VELOCITY_LABEL,
                     NORTH_VELOCITY_LABEL]
        for i,q in enumerate(quantities): 
            #Write C files
            columns = 3 # long, lat , depth
            file = base_name + mux_names[i]
            f = open(file, 'wb')
            files.append(file)
            f.write(pack('i',points_num))
            f.write(pack('i',time_step_count))
            f.write(pack('f',time_step))

            f.close()   
        tide = 1
        try:
            urs2sww(base_name, remove_nc_files=True, mean_stage=tide,
                      verbose=self.verbose)        
        except ANUGAError:
            pass
        else:
            self.delete_mux(files)
            msg = 'Should have raised exception'
            raise msg
        sww_file = base_name + '.sww'
        self.delete_mux(files)
        
    def test_urs2sww_test_fail2(self):
        base_name = 'Harry-high-pants'
        try:
            urs2sww(base_name)        
        except IOError:
            pass
        else:
            self.delete_mux(files)
            msg = 'Should have raised exception'
            raise msg
           
    def test_urs2sww(self):
        tide = 1
        base_name, files = self.create_mux()
        urs2sww(base_name
                #, origin=(0,0,0)
                , mean_stage=tide
                , remove_nc_files=True,
                      verbose=self.verbose
                )
        sww_file = base_name + '.sww'
        
        #Let's interigate the sww file
        # Note, the sww info is not gridded.  It is point data.
        fid = NetCDFFile(sww_file)

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        geo_reference = Geo_reference(NetCDFObject=fid)

        
        #Check that first coordinate is correctly represented       
        #Work out the UTM coordinates for first point
        zone, e, n = redfearn(-34.5, 150.66667)
       
        assert num.allclose(geo_reference.get_absolute([[x[0],y[0]]]), [e,n])

        # Make x and y absolute
        points = geo_reference.get_absolute(map(None, x, y))
        points = ensure_numeric(points)
        x = points[:,0]
        y = points[:,1]
        
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
        # -(-elevation) since elevation is inverted in mux files
        #momentum = velocity_va *(stage+elevation)
        # = e*(e+tide+n) based on how I'm writing these files
        answer_y = e*(e+tide+n) * -1 # talking into account mux file format
        actual_y = ymomentum[0,0]
        #print "answer_y",answer_y
        #print "actual_y",actual_y 
        assert num.allclose(answer_y, actual_y)  #Meters
        
        assert num.allclose(answer_x, actual_x)  #Meters
        
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
        
    def test_urs2sww_momentum(self):
        tide = 1
        time_step_count = 3
        time_step = 2
        #lat_long_points =[(-21.5,114.5),(-21.5,115),(-21.,114.5), (-21.,115.)]
        # This is gridded
        lat_long_points =[(-21.5,114.5),(-21,114.5),(-21.5,115), (-21.,115.)]
        depth=20
        ha=2
        ua=5
        va=-10 #-ve added to take into account mux file format where south
               # is positive.
        base_name, files = self.write_mux(lat_long_points,
                                          time_step_count, time_step,
                                          depth=depth,
                                          ha=ha,
                                          ua=ua,
                                          va=va)
        # write_mux(self,lat_long_points, time_step_count, time_step,
        #          depth=None, ha=None, ua=None, va=None
        urs2sww(base_name
                #, origin=(0,0,0)
                , mean_stage=tide
                , remove_nc_files=True,
                      verbose=self.verbose
                )
        sww_file = base_name + '.sww'
        
        #Let's interigate the sww file
        # Note, the sww info is not gridded.  It is point data.
        fid = NetCDFFile(sww_file)

        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        geo_reference = Geo_reference(NetCDFObject=fid)
        
        #Check first value
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        ymomentum = fid.variables['ymomentum'][:]
        elevation = fid.variables['elevation'][:]
        #assert allclose(stage[0,0], e + tide)  #Meters
        #print "xmomentum", xmomentum
        #print "ymomentum", ymomentum
        #Check the momentums - ua
        #momentum = velocity*water height
        #water height = mux_depth + mux_height +tide
        #water height = mux_depth + mux_height +tide
        #momentum = velocity*(mux_depth + mux_height +tide)
        #
        
        answer = 115
        actual = xmomentum[0,0]
        assert num.allclose(answer, actual)  #Meters^2/ sec
        answer = 230
        actual = ymomentum[0,0]
        #print "answer",answer
        #print "actual",actual 
        assert num.allclose(answer, actual)  #Meters^2/ sec

        # check the stage values, first time step.
        # These arrays are equal since the Easting values were used as
        # the stage

        #assert allclose(stage[0], x +tide)  #Meters

        # check the elevation values.
        # -ve since urs measures depth, sww meshers height,
        # these arrays are equal since the northing values were used as
        # the elevation
        #assert allclose(-elevation, y)  #Meters
        
        fid.close()
        self.delete_mux(files)
        os.remove(sww_file)
        
  
    def test_urs2sww_origin(self):
        tide = 1
        base_name, files = self.create_mux()
        urs2sww(base_name
                , origin=(0,0,0)
                , mean_stage=tide
                , remove_nc_files=True,
                      verbose=self.verbose
                )
        sww_file = base_name + '.sww'
        
        #Let's interigate the sww file
        # Note, the sww info is not gridded.  It is point data.
        fid = NetCDFFile(sww_file)

        #  x and y are absolute
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        geo_reference = Geo_reference(NetCDFObject=fid)

        
        time = fid.variables['time'][:]
        #print "time", time
        assert num.allclose([0.,0.5,1.], time)
        assert fid.starttime == 0.0
        #Check that first coordinate is correctly represented       
        #Work out the UTM coordinates for first point
        zone, e, n = redfearn(-34.5, 150.66667)       
       
        assert num.allclose([x[0],y[0]], [e,n])

        
        #Check first value
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        ymomentum = fid.variables['ymomentum'][:]
        elevation = fid.variables['elevation'][:]
        assert num.allclose(stage[0,0], e +tide)  #Meters

        #Check the momentums - ua
        #momentum = velocity*(stage-elevation)
        #momentum = velocity*(stage+elevation)
        # -(-elevation) since elevation is inverted in mux files
        # = n*(e+tide+n) based on how I'm writing these files
        answer = n*(e+tide+n)
        actual = xmomentum[0,0]
        assert num.allclose(answer, actual)  #Meters

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
 
    def test_urs2sww_minmaxlatlong(self):
        
        #longitudes = [150.66667, 150.83334, 151., 151.16667]
        #latitudes = [-34.5, -34.33333, -34.16667, -34]

        tide = 1
        base_name, files = self.create_mux()
        urs2sww(base_name,
                minlat=-34.5,
                maxlat=-34,
                minlon= 150.66667,
                maxlon= 151.16667,
                mean_stage=tide,
                remove_nc_files=True,
                      verbose=self.verbose
                )
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
        zone, e, n = redfearn(-34.5, 150.66667) 
        assert num.allclose([x[0],y[0]], [e,n])

        
        #Check first value
        stage = fid.variables['stage'][:]
        xmomentum = fid.variables['xmomentum'][:]
        ymomentum = fid.variables['ymomentum'][:]
        elevation = fid.variables['elevation'][:]
        assert num.allclose(stage[0,0], e +tide)  #Meters

        #Check the momentums - ua
        #momentum = velocity*(stage-elevation)
        #momentum = velocity*(stage+elevation)
        # -(-elevation) since elevation is inverted in mux files
        # = n*(e+tide+n) based on how I'm writing these files
        answer = n*(e+tide+n)
        actual = xmomentum[0,0]
        assert num.allclose(answer, actual)  #Meters

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
        
    def test_urs2sww_minmaxmintmaxt(self):
        
        #longitudes = [150.66667, 150.83334, 151., 151.16667]
        #latitudes = [-34.5, -34.33333, -34.16667, -34]

        tide = 1
        base_name, files = self.create_mux()
        
        urs2sww(base_name,
                mint=0.25,
                maxt=0.75,
                mean_stage=tide,
                remove_nc_files=True,
                verbose=self.verbose)
        sww_file = base_name + '.sww'
        
        #Let's interigate the sww file
        # Note, the sww info is not gridded.  It is point data.
        fid = NetCDFFile(sww_file)

        
        time = fid.variables['time'][:]
        assert num.allclose(time, [0.0]) # the time is relative
        assert fid.starttime == 0.5
        
        fid.close()
        self.delete_mux(files)
        #print "sww_file", sww_file
        os.remove(sww_file)
  


    def test_read_asc(self):
        """Test conversion from dem in ascii format to native NetCDF format
        """

        import time, os

        from file_conversion import _read_asc
        #Write test asc file
        filename = tempfile.mktemp(".000")
        fid = open(filename, 'w')
        fid.write("""ncols         7
nrows         4
xllcorner     2000.5
yllcorner     3000.5
cellsize      25
NODATA_value  -9999
    97.921    99.285   125.588   180.830   258.645   342.872   415.836
   473.157   514.391   553.893   607.120   678.125   777.283   883.038
   984.494  1040.349  1008.161   900.738   730.882   581.430   514.980
   502.645   516.230   504.739   450.604   388.500   338.097   514.980
""")
        fid.close()
        bath_metadata, grid = _read_asc(filename, verbose=self.verbose)
        self.failUnless(bath_metadata['xllcorner']  == 2000.5,  'Failed')
        self.failUnless(bath_metadata['yllcorner']  == 3000.5,  'Failed')
        self.failUnless(bath_metadata['cellsize']  == 25,  'Failed')
        self.failUnless(bath_metadata['NODATA_value']  == -9999,  'Failed')
        self.failUnless(grid[0][0]  == 97.921,  'Failed')
        self.failUnless(grid[3][6]  == 514.980,  'Failed')

        os.remove(filename)


    #### TESTS FOR URS 2 SWW  ###     
    
    def create_mux(self, points_num=None):
        # write all the mux stuff.
        time_step_count = 3
        time_step = 0.5
        
        longitudes = [150.66667, 150.83334, 151., 151.16667]
        latitudes = [-34.5, -34.33333, -34.16667, -34]

        if points_num == None:
            points_num = len(longitudes) * len(latitudes)

        lonlatdeps = []
        quantities = ['HA','UA','VA']
        mux_names = [WAVEHEIGHT_MUX_LABEL,
                     EAST_VELOCITY_LABEL,
                     NORTH_VELOCITY_LABEL]
        quantities_init = [[],[],[]]
        # urs binary is latitude fastest
        for i,lon in enumerate(longitudes):
            for j,lat in enumerate(latitudes):
                _ , e, n = redfearn(lat, lon)
                lonlatdeps.append([lon, lat, n])
                quantities_init[0].append(e) # HA
                quantities_init[1].append(n ) # UA
                quantities_init[2].append(e) # VA
        #print "lonlatdeps",lonlatdeps

        file_handle, base_name = tempfile.mkstemp("")
        os.close(file_handle)
        os.remove(base_name)
        
        files = []        
        for i,q in enumerate(quantities): 
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
        

#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_File_Conversion,'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)

