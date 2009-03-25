#!/usr/bin/env python


import unittest
import tempfile
import Numeric as num
import zlib
from os.path import join, split, sep
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a


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
            test_array = num.array([[7.0, 3.14], [-31.333, 0.0]])

            # First file
            filename1 = mktemp(suffix='.nc', dir='.')
            fid = NetCDFFile(filename1, netcdf_mode_w)
            fid.createDimension('two', 2)
            fid.createVariable('test_array', num.Float,
                               ('two', 'two'))
            fid.variables['test_array'][:] = test_array
            fid.close()
            
            # Second file
            filename2 = mktemp(suffix='.nc', dir='.')
            fid = NetCDFFile(filename2, netcdf_mode_w)
            fid.createDimension('two', 2)
            fid.createVariable('test_array', num.Float,
                               ('two', 'two'))
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
        

    def test_get_vars_in_expression(self):
        '''Test the 'get vars from expression' code.'''

        def test_it(source, expected):
            result = get_vars_in_expression(source)
            result.sort()
            expected.sort()
            msg = ("Source: '%s'\nResult: %s\nExpected: %s"
                   % (source, str(result), str(expected)))
            self.failUnlessEqual(result, expected, msg)
                
        source = 'fred'
        expected = ['fred']
        test_it(source, expected)

        source = 'tom + dick'
        expected = ['tom', 'dick']
        test_it(source, expected)

        source = 'tom * (dick + harry)'
        expected = ['tom', 'dick', 'harry']
        test_it(source, expected)

        source = 'tom + dick**0.5 / (harry - tom)'
        expected = ['tom', 'dick', 'harry']
        test_it(source, expected)

    def test_tar_untar_files(self):
        '''Test that tarring & untarring files is OK.'''

        # these test files must exist in the current directory
        files = ('test_system_tools.py', 'system_tools.py')

        # name of tar file and test (temp) directory
        tar_filename = 'test.tgz'
        tmp_dir = tempfile.mkdtemp()

        # tar and untar the test files into a temporary directory
        tar_file(files, tar_filename)
        untar_file(tar_filename, tmp_dir)

        # see if original files and untarred ones are the same
        for file in files:
            fd = open(file, 'r')
            orig = fd.readlines()
            fd.close()

            fd = open(os.path.join(tmp_dir, file), 'r')
            copy = fd.readlines()
            fd.close()

            msg = "Original file %s isn't the same as untarred copy?" % file
            self.failUnless(orig == copy, msg)

        os.remove(tar_filename)

#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_system_tools, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)

