import sys
import unittest
import numpy as num
import copy
import os

# ANUGA modules
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a, \
                            netcdf_float

from anuga.file.netcdf import NetCDFFile

from anuga.file_conversion.dem2pts import dem2pts
from anuga.file_conversion.asc2dem import asc2dem

class Test_Dem2Pts(unittest.TestCase):
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
     
    def test_dem2pts_bounding_box_v2(self):
        """Test conversion from dem in ascii format to native NetCDF format
        """

        import time, os


        #Write test asc file
        root = 'demtest'

        filename = root+'.asc'
        fid = open(filename, 'w')
        fid.write("""ncols         10
nrows         10
xllcorner     2000
yllcorner     3000
cellsize      1
NODATA_value  -9999
""")
        #Create linear function
        ref_points = []
        ref_elevation = []
        x0 = 2000
        y = 3010
        yvec = range(10)
        xvec = range(10)
        z = -1
        for i in range(10):
            y = y - 1
            for j in range(10):
                x = x0 + xvec[j]
                z += 1
                ref_points.append ([x,y])
                ref_elevation.append(z)
                fid.write('%f ' %z)
            fid.write('\n')

        fid.close()

        #print 'sending pts', ref_points
        #print 'sending elev', ref_elevation

        #Write prj file with metadata
        metafilename = root+'.prj'
        fid = open(metafilename, 'w')


        fid.write("""Projection UTM
Zone 56
Datum WGS84
Zunits NO
Units METERS
Spheroid WGS84
Xshift 0.0000000000
Yshift 10000000.0000000000
Parameters
""")
        fid.close()

        #Convert to NetCDF pts
        dem2pts(filename, easting_min=2002.0, easting_max=2007.0,
                northing_min=3003.0, northing_max=3006.0,
                verbose=False)

        #Check contents
        #Get NetCDF
        fid = NetCDFFile(root+'.pts', netcdf_mode_r)

        # Get the variables
        #print fid.variables.keys()
        points = fid.variables['points']
        elevation = fid.variables['elevation']

        #Check values
        assert fid.xllcorner == 2002.0
        assert fid.yllcorner == 3003.0

        #create new reference points
        newz = []
        newz[0:5] = ref_elevation[32:38]
        newz[6:11] = ref_elevation[42:48]
        newz[12:17] = ref_elevation[52:58]
        newz[18:23] = ref_elevation[62:68]
        ref_elevation = []
        ref_elevation = newz
        ref_points = []
        x0 = 2002
        y = 3007
        yvec = range(4)
        xvec = range(6)
        for i in range(4):
            y = y - 1
            ynew = y - 3003.0
            for j in range(6):
                x = x0 + xvec[j]
                xnew = x - 2002.0
                ref_points.append ([xnew,ynew]) #Relative point values

        assert num.allclose(points, ref_points)

        assert num.allclose(elevation, ref_elevation)

        #Cleanup
        fid.close()


        os.remove(root + '.pts')
        os.remove(root + '.dem')
        os.remove(root + '.asc')
        os.remove(root + '.prj')


    def test_dem2pts_bounding_box_removeNullvalues_v2(self):
        """Test conversion from dem in ascii format to native NetCDF format
        """

        import time, os

        #Write test asc file
        root = 'demtest2'

        filename = root+'.asc'
        fid = open(filename, 'w')
        fid.write("""ncols         10
nrows         10
xllcorner     2000
yllcorner     3000
cellsize      1
NODATA_value  -9999
""")
        #Create linear function
        ref_points = []
        ref_elevation = []
        x0 = 2000
        y = 3010
        yvec = range(10)
        xvec = range(10)
        #z = range(100)
        z = num.zeros(100, num.int)     #array default#
        NODATA_value = -9999
        count = -1
        for i in range(10):
            y = y - 1
            for j in range(10):
                x = x0 + xvec[j]
                ref_points.append ([x,y])
                count += 1
                z[count] = (4*i - 3*j)%13
                if j == 4: z[count] = NODATA_value #column inside clipping region
                if j == 8: z[count] = NODATA_value #column outside clipping region
                if i == 9: z[count] = NODATA_value #row outside clipping region
                if i == 4 and j == 6: z[count] = NODATA_value #arbitrary point inside clipping region
                ref_elevation.append( z[count] )
                fid.write('%f ' %z[count])
            fid.write('\n')

        fid.close()

        #print 'sending elev', ref_elevation

        #Write prj file with metadata
        metafilename = root+'.prj'
        fid = open(metafilename, 'w')


        fid.write("""Projection UTM
Zone 56
Datum WGS84
Zunits NO
Units METERS
Spheroid WGS84
Xshift 0.0000000000
Yshift 10000000.0000000000
Parameters
""")
        fid.close()

        #Convert to NetCDF pts
        dem2pts(filename, easting_min=2002.0, easting_max=2007.0,
                northing_min=3003.0, northing_max=3006.0)

        #Check contents
        #Get NetCDF
        fid = NetCDFFile(root+'.pts', netcdf_mode_r)

        # Get the variables
        #print fid.variables.keys()
        points = fid.variables['points'][:]
        elevation = fid.variables['elevation'][:]

        #Check values
        assert fid.xllcorner == 2002.0
        assert fid.yllcorner == 3003.0

        #create new reference points
        newz = num.zeros(19, num.int)       #array default#
        newz[0:2] = ref_elevation[32:34]
        newz[2:5] = ref_elevation[35:38]
        newz[5:7] = ref_elevation[42:44]
        newz[7] = ref_elevation[45]
        newz[8] = ref_elevation[47]
        newz[9:11] = ref_elevation[52:54]
        newz[11:14] = ref_elevation[55:58]
        newz[14:16] = ref_elevation[62:64]
        newz[16:19] = ref_elevation[65:68]


        ref_elevation = newz
        ref_points = []
        new_ref_points = []
        x0 = 2002
        y = 3007
        yvec = range(4)
        xvec = range(6)
        for i in range(4):
            y = y - 1
            ynew = y - 3003.0
            for j in range(6):
                x = x0 + xvec[j]
                xnew = x - 2002.0
                if j <> 2 and (i<>1 or j<>4):
                    ref_points.append([x,y])
                    new_ref_points.append ([xnew,ynew])



        #print points
        #print new_ref_points

        assert num.allclose(points, new_ref_points)

        #print elevation
        #print ref_elevation
        
        assert num.allclose(elevation, ref_elevation)

        #Cleanup
        fid.close()


        #os.remove(root + '.pts')
        os.remove(root + '.dem')
        os.remove(root + '.asc')
        os.remove(root + '.prj')


    def test_dem2pts_bounding_box_removeNullvalues_v3(self):
        """Test conversion from dem in ascii format to native NetCDF format
        Check missing values on clipping boundary
        """

        import time, os

        #Write test asc file
        root = 'demtest3'

        filename = root+'.asc'
        fid = open(filename, 'w')
        fid.write("""ncols         10
nrows         10
xllcorner     2000
yllcorner     3000
cellsize      1
NODATA_value  -9999
""")
        #Create linear function
        ref_points = []
        ref_elevation = []
        x0 = 2000
        y = 3010
        yvec = range(10)
        xvec = range(10)
        #z = range(100)
        z = num.zeros(100, num.int)     #array default#
        NODATA_value = -9999
        count = -1
        for i in range(10):
            y = y - 1
            for j in range(10):
                x = x0 + xvec[j]
                ref_points.append ([x,y])
                count += 1
                z[count] = (4*i - 3*j)%13
                if j == 4: z[count] = NODATA_value #column inside clipping region
                if j == 8: z[count] = NODATA_value #column outside clipping region
                if i == 6: z[count] = NODATA_value #row on clipping boundary
                if i == 4 and j == 6: z[count] = NODATA_value #arbitrary point inside clipping region
                ref_elevation.append( z[count] )
                fid.write('%f ' %z[count])
            fid.write('\n')

        fid.close()

        #print 'sending elev', ref_elevation

        #Write prj file with metadata
        metafilename = root+'.prj'
        fid = open(metafilename, 'w')


        fid.write("""Projection UTM
Zone 56
Datum WGS84
Zunits NO
Units METERS
Spheroid WGS84
Xshift 0.0000000000
Yshift 10000000.0000000000
Parameters
""")
        fid.close()

        #Convert to NetCDF pts
        dem2pts(filename, easting_min=2002.0, easting_max=2007.0,
                northing_min=3003.0, northing_max=3006.0)

        #Check contents
        #Get NetCDF
        fid = NetCDFFile(root+'.pts', netcdf_mode_r)

        # Get the variables
        #print fid.variables.keys()
        points = fid.variables['points']
        elevation = fid.variables['elevation']

        #Check values
        assert fid.xllcorner == 2002.0
        assert fid.yllcorner == 3003.0

        #create new reference points
        newz = num.zeros(14, num.int)       #array default#
        newz[0:2] = ref_elevation[32:34]
        newz[2:5] = ref_elevation[35:38]
        newz[5:7] = ref_elevation[42:44]
        newz[7] = ref_elevation[45]
        newz[8] = ref_elevation[47]
        newz[9:11] = ref_elevation[52:54]
        newz[11:14] = ref_elevation[55:58]



        ref_elevation = newz
        ref_points = []
        new_ref_points = []
        x0 = 2002
        y = 3007
        yvec = range(4)
        xvec = range(6)
        for i in range(4):
            y = y - 1
            ynew = y - 3003.0
            for j in range(6):
                x = x0 + xvec[j]
                xnew = x - 2002.0
                if j <> 2 and (i<>1 or j<>4) and i<>3:
                    ref_points.append([x,y])
                    new_ref_points.append ([xnew,ynew])


        #print points[:],points[:].shape
        #print new_ref_points, len(new_ref_points)

        assert num.allclose(elevation, ref_elevation)
        assert num.allclose(points, new_ref_points)


        #Cleanup
        fid.close()


        os.remove(root + '.pts')
        os.remove(root + '.dem')
        os.remove(root + '.asc')
        os.remove(root + '.prj')


#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Dem2Pts,'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)
