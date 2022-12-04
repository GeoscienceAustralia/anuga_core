
from builtins import range
import sys
import unittest
import numpy as num
import copy
import os

# ANUGA modules
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a, \
                            netcdf_float

from anuga.file.netcdf import NetCDFFile

from anuga.file_conversion.llasc2pts import llasc2pts
from anuga.coordinate_transforms.redfearn import convert_from_latlon_to_utm

class Test_LLAsc2Pts(unittest.TestCase):
    """ A suite of tests to test file conversion functions.
        These tests are quite coarse-grained: converting a file
        and checking that its headers and some of its contents
        are correct.
    """ 

    def tearDown(self):
        for file in ['demtest2.pts']:
            try:
                os.remove(file)
            except:
                pass 
     
    def test_llasc2pts_bounding_box_v2(self):
        """Test conversion from dem in ascii format to native NetCDF format
        """

        import time, os

        ncols = 5
        nrows = 4
    
        xllcorner = 114.0
        yllcorner = -20.0
        cellsize = 0.25 # degrees
        NODATA = -9999

        #Write test asc file
        root = 'llasctest'

        filename = root+'.asc'
        fid = open(filename, 'w')
        fid.write("""ncols         %g
nrows         %g
xllcorner     %f
yllcorner     %f
cellsize      %f
NODATA_value  %g
"""% (ncols, nrows, xllcorner, yllcorner, cellsize, NODATA))
        #Create linear function
        ref_points = []
        ref_elevation = []
        ref_degrees = []
        ref_zones = []

        z = -1
        for i in range(nrows):
            lat = yllcorner + (nrows-1)*cellsize - i*cellsize 
            for j in range(ncols):
                lon = xllcorner + j*cellsize
                z += 1
                ref_degrees.append ([lat,lon])
                utms, zone = convert_from_latlon_to_utm(points = [lat,lon])
                #print([lat,lon],utms,zone)
                ref_points.append(utms[0])
                ref_elevation.append(z)
                ref_zones.append(zone)
                fid.write('%f ' %z)
            fid.write('\n')

        fid.close()

        utm_corner, zone = convert_from_latlon_to_utm(points = [yllcorner,xllcorner])
        ref_points = num.asarray(ref_points)
        ref_elevation = num.asarray(ref_elevation)
        #print 'sending pts', ref_points
        #print 'sending elev', ref_elevation


        #Convert to NetCDF pts
        llasc2pts(filename, verbose=False)

        #Check contents
        #Get NetCDF
        fid = NetCDFFile(root+'.pts', netcdf_mode_r)

        # Get the variables
        #print fid.variables.keys()
        points = fid.variables['points']
        elevation = fid.variables['elevation']

        #import ipdb; ipdb.set_trace()

        #Check values
        assert num.allclose(fid.xllcorner, 186073.679565389)
        assert num.allclose(fid.yllcorner, 7785705.97374621)
        assert fid.zone == 50
        assert num.allclose(points[:], ref_points-utm_corner)
        assert num.allclose(elevation[:], ref_elevation)

        #Cleanup
        fid.close()

        os.remove(root + '.asc')

        try:
            os.remove(root + '.pts')
            os.remove(root + '.dem')
            os.remove(root + '.asc')
            os.remove(root + '.prj')
        except:
            pass





#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_LLAsc2Pts,'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)
