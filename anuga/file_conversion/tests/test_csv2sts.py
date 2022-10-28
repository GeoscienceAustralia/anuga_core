#external modules
from builtins import str
import os
import sys
import unittest
import numpy as num
from anuga.file.netcdf import NetCDFFile
from anuga.utilities.system_tools import get_pathname_from_package



# ANUGA modules
from anuga.config import netcdf_float32, netcdf_float64
from anuga.file.csv_file import load_csv_as_dict

# Local modules
from anuga.file_conversion.csv2sts import csv2sts

# some test files we want to generate
testfile_csv = 'small___.csv'
sts_out = 'sts_out.sts'

lat = 10
lon = 20



class Test_csv2sts(unittest.TestCase):
    """
        Test csv to NetCDFFile conversion functionality.
    """
    def setUp(self):
        """ Setup for all tests. """
        self.verbose = True
        fid = open(testfile_csv, 'w')
        fid.write("""time stage
0 4
1 150.66667
2 150.83334
3 151.
4 151.16667
5 -34.
6 -34.16667
7 -34.33333
8 -34.5
9 -1.
10 -5.
11 -9.
12 -13.
""")
        fid.close()
                  
    def tearDown(self):
        """ Cleanup for all tests. """     
        os.remove(testfile_csv)

    def test_missing_input_file(self):
        """
        Test that a missing csv file raises the correct exception.
        """
        got_except = False
        
        try:
            csv2sts('somename_not_here.csv', sts_out, 10, 20)
        except IOError as e:
            got_except = True
        except:
            assert False, 'Missing file raised wrong exception.'

        assert got_except is True, 'Missing file did not raise an exception.'

    def test_csv2sts_output(self):
        """
        Test that a csv file is correctly rendered to .sts (NetCDF) format.
        """
        csv2sts(testfile_csv, sts_out, latitude = lat, longitude = lon)
        self._check_generated_sts()
        
    def test_run_via_commandline(self):
        """
        Make sure that the python file functions as a command-line tool.
        """

        # Look for script in same dir as this unit test.
        path = os.path.dirname( os.path.realpath( __file__ ) )

        path = get_pathname_from_package( 'anuga.file_conversion' )
        
        cmd = sys.executable + ' ' + os.path.join( path, 'csv2sts.py') + ' --latitude ' 
        cmd += '%s --lon %s %s %s' % (str(lat), str(lon), testfile_csv, sts_out)
        
        os.system(cmd)
        self._check_generated_sts()


    def _check_generated_sts(self):
        """ check that we can read data out of the file """
        sts = NetCDFFile(sts_out,'r')
        
        data, names = load_csv_as_dict(testfile_csv, delimiter=' ', d_type = num.float64)
        
        assert sts.latitude == lat, 'latitude does not match'
        assert sts.longitude == lon, 'longitude does not match'
        
        assert len(sts.variables) == len(data), 'num variables does not match'
        
        # make sure data is returned in exactly the expected format
        for key, values in list(data.items()):
            assert list(sts.variables[key][:]) == values, \
                                        'stored data does not match'

        if not sys.platform == 'win32':
            # Windows cannot delete the file for some reason.
            os.remove(sts_out)           

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_csv2sts,'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)
