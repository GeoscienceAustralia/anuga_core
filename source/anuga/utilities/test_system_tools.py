#!/usr/bin/env python


import unittest
from Numeric import zeros, array, allclose, Float
from Scientific.IO.NetCDF import NetCDFFile

from system_tools import *

class Test_system_tools(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_user_name(self):
        user = get_user_name()

        # print user
        assert isinstance(user, basestring), 'User name should be a string'

    def test_host_name(self):
        host = get_host_name()

        # print host
        assert isinstance(host, basestring), 'User name should be a string'        

    def test_compute_checksum(self):
        """test_compute_checksum(self):

        Check that checksums on files are OK
        """

        import zlib
        from tempfile import NamedTemporaryFile, mktemp

        # Generate a text file
        fid = NamedTemporaryFile(mode='w',
                                 suffix='.tmp',
                                 dir='.')
        string = 'My temp file with textual content. AAAABBBBCCCC1234'
        fid.write(string)
        fid.flush()


        ref_crc = zlib.crc32(string)
        checksum = compute_checksum(fid.name)

        assert checksum == ref_crc

        # Close and remove temporary file
        fid.close()



        # Binary file
        fid = NamedTemporaryFile(mode='w+b',
                                 suffix='.tmp',
                                 dir='.')
        string = 'My temp file with binary content. AAAABBBBCCCC1234'
        fid.write(string)
        fid.flush()


        ref_crc = zlib.crc32(string)
        checksum = compute_checksum(fid.name)

        assert checksum == ref_crc

        # Close and remove temporary file
        fid.close()
        
        
        # Binary NetCDF File X 2 (use mktemp's name)

        test_array = array([[7.0, 3.14], [-31.333, 0.0]])

        # First file
        filename1 = mktemp(suffix='.nc', dir='.')
        fid = NetCDFFile(filename1, 'w')
        fid.createDimension('two', 2)
        fid.createVariable('test_array', Float,
                           ('two', 'two'))
        fid.variables['test_array'][:] = test_array
        fid.close()

        # Second file
        filename2 = mktemp(suffix='.nc', dir='.')
        fid = NetCDFFile(filename2, 'w')
        fid.createDimension('two', 2)
        fid.createVariable('test_array', Float,
                           ('two', 'two'))
        fid.variables['test_array'][:] = test_array
        fid.close()

        
        checksum1 = compute_checksum(filename1)
        checksum2 = compute_checksum(filename2)        
        assert checksum1 == checksum2


        os.remove(filename1)
        os.remove(filename2)

        
#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_system_tools, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)




