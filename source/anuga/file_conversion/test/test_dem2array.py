import sys
import unittest
import numpy as num
import copy
import os

# ANUGA modules
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a, \
                            netcdf_float

from anuga.file.netcdf import NetCDFFile

from anuga.file_conversion.dem2array import dem2array
from anuga.file_conversion.asc2dem import asc2dem


class Test_dem2array(unittest.TestCase):
    """ A suite of tests to test file conversion functions.
        These tests are quite coarse-grained: converting a file
        and checking that its headers and some of its contents
        are correct.
    """ 
 
    def test_dem2array(self):
        """Test conversion from dem to array
        """

        import time, os


        #Write test asc file
        root = 'dem2arraytest_1'

        filename = root+'.asc'
        fid = open(filename, 'w')
        fid.write("""ncols         11
nrows         12
xllcorner     2000
yllcorner     3000
cellsize      3
NODATA_value  -9999
""")
        #Create linear function
        ref_points = []
        ref_elevation = []
        x0 = 2000
        y = 3000
        cellsize = 3.0
        xvec = x0 + cellsize*num.arange(11)
        yvec = y + cellsize*(num.arange(12))
        z = -1
        for i in range(12):
            y = y - cellsize
            for j in range(11):
                x = x0 + xvec[j]
                z += 1
                ref_points.append ([x,y])
                ref_elevation.append(z)
                fid.write('%f ' %z)
            fid.write('\n')

        fid.close()

        #print 'sending pts', ref_points
        #print 'sending elev', ref_elevation
        
        Z_ex = num.array(ref_elevation,num.float).reshape(12,11)

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

        #Convert to NetCDF dem
        
        asc2dem(filename)


        x,y, Z = dem2array(root+'.dem')
        
#         print Z
#         print Z_ex
#         print x
#         print xvec
#         print y
#         print yvec        
        
        assert num.allclose(Z,Z_ex)
        assert num.allclose(x,xvec)
        assert num.allclose(y,yvec)

               
        #assert num.allclose(x,)

        try:
            os.remove(root + '.dem')
            os.remove(root + '.asc')
            os.remove(root + '.prj')
        except:
            # Expect error on windows
            pass


#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_dem2array,'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)
