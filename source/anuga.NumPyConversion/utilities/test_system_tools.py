#!/usr/bin/env python


import unittest
from numpy import zeros, array, allclose, float, dtype
import zlib
from os.path import join, split, sep


# Please, don't add anuga.utilities to these imports.
# I'm trying to keep this file general, so it works for EQRM and ANUGA
# EQRM also uses this file, but has a different directory structure
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

        from tempfile import mkstemp, mktemp

        # Generate a text file
        tmp_fd , tmp_name = mkstemp(suffix='.tmp', dir='.')
        fid = os.fdopen(tmp_fd, 'w+b')
        string = 'My temp file with textual content. AAAABBBBCCCC1234'
        fid.write(string)
        fid.close()

        # Have to apply the 64 bit fix here since we aren't comparing two
        # files, but rather a string and a file.
        ref_crc = safe_crc(string)

        checksum = compute_checksum(tmp_name)
        assert checksum == ref_crc

        os.remove(tmp_name)
        


        # Binary file
        tmp_fd , tmp_name = mkstemp(suffix='.tmp', dir='.')
        fid = os.fdopen(tmp_fd, 'w+b')

        string = 'My temp file with binary content. AAAABBBBCCCC1234'
        fid.write(string)
        fid.close()

        ref_crc = safe_crc(string)
        checksum = compute_checksum(tmp_name)

        assert checksum == ref_crc

        os.remove(tmp_name)        
        
        # Binary NetCDF File X 2 (use mktemp's name)

        try:
            from Scientific.IO.NetCDF import NetCDFFile
        except ImportError:
            # This code is also used by EQRM which does not require NetCDF
            pass
        else:
            test_array = array([[7.0, 3.14], [-31.333, 0.0]])

            # First file
            filename1 = mktemp(suffix='.nc', dir='.')
            fid = NetCDFFile(filename1, 'w')
            fid.createDimension('two', 2)
            fid.createVariable('test_array', dtype(float).char, ('two', 'two'))
            fid.variables['test_array'][:] = test_array
            fid.close()
            
            # Second file
            filename2 = mktemp(suffix='.nc', dir='.')
            fid = NetCDFFile(filename2, 'w')
            fid.createDimension('two', 2)
            fid.createVariable('test_array', dtype(float).char, ('two', 'two'))
            fid.variables['test_array'][:] = test_array
            fid.close()
            
            
            checksum1 = compute_checksum(filename1)
            checksum2 = compute_checksum(filename2)        
            assert checksum1 == checksum2


            os.remove(filename1)
            os.remove(filename2)


    def test_compute_checksum_real(self):
        """test_compute_checksum(self):

        Check that checksums on a png file is OK
        """

        # Get path where this test is run
        # I'm trying to keep this file general, so it works for EQRM and ANUGA
        path, tail = split(__file__)
        if path == '':
            path = '.' + sep
        
        filename = path + sep +  'crc_test_file.png'

        ref_crc = 1203293305 # Computed on Windows box
        checksum = compute_checksum(filename)

        msg = 'Computed checksum = %s, should have been %s'\
              %(checksum, ref_crc)
        assert checksum == ref_crc, msg
        #print checksum
        
        
#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_system_tools, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)




