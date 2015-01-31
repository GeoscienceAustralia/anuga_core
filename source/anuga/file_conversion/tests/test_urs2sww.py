
# External modules
from anuga.file.netcdf import NetCDFFile
import sys
import unittest
import numpy as num
import copy
import os

# ANUGA modules
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a, \
                            netcdf_float

from anuga.coordinate_transforms.geo_reference import Geo_reference, \
     write_NetCDF_georeference
from anuga.coordinate_transforms.redfearn import redfearn     

from anuga.file_conversion.urs2sww import urs2sww, urs_ungridded2sww
import anuga.file.urs as urs

from anuga.file.mux import WAVEHEIGHT_MUX_LABEL, EAST_VELOCITY_LABEL, \
                            NORTH_VELOCITY_LABEL

from anuga.geospatial_data.geospatial_data import ensure_absolute
from anuga.utilities.numerical_tools import ensure_numeric  

# use helper methods from other unit test
from anuga.file.tests.test_mux import Test_Mux



class Test_Dem2Pts(Test_Mux):
    """ A suite of tests to test urs2sww file conversion functions.
        These tests are quite coarse-grained: converting a file
        and checking that its headers and some of its contents
        are correct.
    """ 

    def setUp(self):
        self.verbose = False # change this to output more debug info

  
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
        urs_ungridded2sww(base_name, mean_stage=tide, origin =(50,23432,4343))
        
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
                          hole_points_UTM=[ 292110.784, 7676551.710 ])
        
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
                          hole_points_UTM=[ 292110.784, 7676551.710 ])
        
        # now I want to check the sww file ...
        sww_file = base_name + '.sww'
        fid = NetCDFFile(sww_file)
        
        volumes = fid.variables['volumes'][:]
        #print "number_of_volumes",len(volumes)

        fid.close()
        os.remove(sww_file)
        
        urs_ungridded2sww(base_name, mean_stage=-240992.0)
        
        # now I want to check the sww file ...
        sww_file = base_name + '.sww'
        fid = NetCDFFile(sww_file)
        
        volumes_again = fid.variables['volumes'][:]
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
                          mint=101, maxt=500)
        
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
                          mint=0, maxt=100000)
        
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
            urs2sww(base_name, mean_stage=tide,
                          origin =(50,23432,4343),
                          mint=301, maxt=399,
                              verbose=self.verbose)
        except: 
            pass
        else:
            self.assertTrue(0 ==1, 'Bad input did not throw exception error!')

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
            urs2sww(base_name, mean_stage=tide,
                          origin =(50,23432,4343),
                          mint=301, maxt=301,
                              verbose=self.verbose)
        except: 
            pass
        else:
            self.assertTrue(0 ==1, 'Bad input did not throw exception error!')

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
        geo=urs.calculate_boundary_points(boundary_polygon, zone,
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
        urs_ungridded2sww(base_name, mean_stage=tide, origin = geo_reference)
        
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


#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Dem2Pts,'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)
