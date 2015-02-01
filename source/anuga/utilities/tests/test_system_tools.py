#!/usr/bin/env python


import unittest
import numpy as num
import random
import tempfile
import zlib
import os
from os.path import join, split, sep
from anuga.file.netcdf import NetCDFFile
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a
from anuga.config import netcdf_float, netcdf_char, netcdf_int


# Please, don't add anuga.utilities to these imports.
# I'm trying to keep this file general, so it works for EQRM and ANUGA
# EQRM also uses this file, but has a different directory structure
from anuga.utilities.system_tools import *

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
            from anuga.file.netcdf import NetCDFFile
        except ImportError:
            # This code is also used by EQRM which does not require NetCDF
            pass
        else:
            test_array = num.array([[7.0, 3.14], [-31.333, 0.0]])

            # First file
            filename1 = mktemp(suffix='.nc', dir='.')
            fid = NetCDFFile(filename1, netcdf_mode_w)
            fid.createDimension('two', 2)
            fid.createVariable('test_array', netcdf_float,
                               ('two', 'two'))
            fid.variables['test_array'][:] = test_array
            fid.close()

            # Second file
            filename2 = mktemp(suffix='.nc', dir='.')
            fid = NetCDFFile(filename2, netcdf_mode_w)
            fid.createDimension('two', 2)
            fid.createVariable('test_array', netcdf_float,
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


            
        path = get_pathname_from_package('anuga.utilities')        
                
        filename = os.path.join(path, 'tests', 'data', 'crc_test_file.png') 



        ref_crc = 1203293305 # Computed on Windows box
        checksum = compute_checksum(filename)

        msg = 'Computed checksum = %s, should have been %s'\
              %(checksum, ref_crc)
        assert checksum == ref_crc, msg
        #print checksum

################################################################################
# Test the clean_line() utility function.
################################################################################

    # helper routine to test clean_line()
    def clean_line_helper(self, instr, delim, expected):
        result = clean_line(instr, delim)
        self.assertTrue(result == expected,
                        "clean_line('%s', '%s'), expected %s, got %s"
                        % (str(instr), str(delim), str(expected), str(result)))

    def test_clean_line_01(self):
        self.clean_line_helper('abc, ,,xyz,123', ',', ['abc', '', 'xyz', '123'])

    def test_clean_line_02(self):
        self.clean_line_helper(' abc , ,, xyz  , 123  ', ',',
                             ['abc', '', 'xyz', '123'])

    def test_clean_line_03(self):
        self.clean_line_helper('1||||2', '|', ['1', '2'])

    def test_clean_line_04(self):
        self.clean_line_helper('abc, ,,xyz,123, ', ',',
                             ['abc', '', 'xyz', '123']) 

    def test_clean_line_05(self):
        self.clean_line_helper('abc, ,,xyz,123, ,    ', ',',
                             ['abc', '', 'xyz', '123', ''])

    def test_clean_line_06(self):
        self.clean_line_helper(',,abc, ,,xyz,123, ,    ', ',',
                             ['abc', '', 'xyz', '123', ''])

    def test_clean_line_07(self):
        self.clean_line_helper('|1||||2', '|', ['1', '2'])

    def test_clean_line_08(self):
        self.clean_line_helper(' ,a,, , ,b,c , ,, , ', ',',
                             ['a', '', '', 'b', 'c', '', ''])

    def test_clean_line_09(self):
        self.clean_line_helper('a:b:c', ':', ['a', 'b', 'c'])

    def test_clean_line_10(self):
        self.clean_line_helper('a:b:c:', ':', ['a', 'b', 'c'])

################################################################################
# Test the string_to_char() and char_to_string() utility functions.
################################################################################

    def test_string_to_char(self):
        import random

        MAX_CHARS = 10
        MAX_ENTRIES = 10000
        A_INT = ord('a')
        Z_INT = ord('z')

        # generate some random strings in a list, with guaranteed lengths
        str_list = ['x' * MAX_CHARS]        # make first maximum length
        for entry in xrange(MAX_ENTRIES):
            length = random.randint(1, MAX_CHARS)
            s = ''
            for c in range(length):
                s += chr(random.randint(A_INT, Z_INT))
            str_list.append(s)

        x = string_to_char(str_list)
        new_str_list = char_to_string(x)

        self.failUnlessEqual(new_str_list, str_list)

    # special test - input list is ['']
    def test_string_to_char2(self):
        # generate a special list shown bad in load_mesh testing
        str_list = ['']

        x = string_to_char(str_list)
        new_str_list = char_to_string(x)

        self.failUnlessEqual(new_str_list, str_list)


################################################################################
# Test the raw I/O to NetCDF files of string data encoded/decoded with 
# string_to_char() and char_to_string().
################################################################################

    def helper_write_msh_file(self, filename, l):
        # open the NetCDF file
        fd = NetCDFFile(filename, netcdf_mode_w)
        fd.description = 'Test file - string arrays'

        # convert list of strings to num.array
        al = num.array(string_to_char(l), num.character)

        # write the list
        fd.createDimension('num_of_strings', al.shape[0])
        fd.createDimension('size_of_strings', al.shape[1])

        var = fd.createVariable('strings', netcdf_char,
                                ('num_of_strings', 'size_of_strings'))
        var[:] = al

        fd.close()


    def helper_read_msh_file(self, filename):
        fid = NetCDFFile(filename, netcdf_mode_r)
        mesh = {}

        # Get the 'strings' variable
        strings = fid.variables['strings'][:]

        fid.close()

        return char_to_string(strings)


    # test random strings to a NetCDF file
    def test_string_to_netcdf(self):
        import random

        MAX_CHARS = 10
        MAX_ENTRIES = 10000

        A_INT = ord('a')
        Z_INT = ord('z')

        FILENAME = 'test.msh'

        # generate some random strings in a list, with guaranteed lengths
        str_list = ['x' * MAX_CHARS]        # make first maximum length
        for entry in xrange(MAX_ENTRIES):
            length = random.randint(1, MAX_CHARS)
            s = ''
            for c in range(length):
                s += chr(random.randint(A_INT, Z_INT))
            str_list.append(s)

        self.helper_write_msh_file(FILENAME, str_list)
        new_str_list = self.helper_read_msh_file(FILENAME)

        self.failUnlessEqual(new_str_list, str_list)
        os.remove(FILENAME)

    # special test - list [''] to a NetCDF file
    def test_string_to_netcdf2(self):
        FILENAME = 'test.msh'

        # generate some random strings in a list, with guaranteed lengths
        str_list = ['']

        self.helper_write_msh_file(FILENAME, str_list)
        new_str_list = self.helper_read_msh_file(FILENAME)

        self.failUnlessEqual(new_str_list, str_list)
        os.remove(FILENAME)


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

        num_lines = 100
        line_size = 100

        # these test files must exist in the current directory
        # create them with random data
        files = ('alpha', 'beta', 'gamma')
        for file in files:
            fd = open(file, 'w')
            line = ''
            for i in range(num_lines):
                for j in range(line_size):
                    line += chr(random.randint(ord('A'), ord('Z')))
                line += '\n'
                fd.write(line)
            fd.close()

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
            self.assertTrue(orig == copy, msg)

        # clean up
        for file in files:
            os.remove(file)
        os.remove(tar_filename)


    def test_file_digest(self):
        '''Test that file digest functions give 'correct' answer.
        
        Not a good test as we get 'expected_digest' from a digest file,
        but *does* alert us if the digest algorithm ever changes.
        '''

        # we expect this digest string from the data file
        expected_digest = '831a1dde6edd365ec4163a47871fa21b'

        # prepare test directory and filenames
        tmp_dir = tempfile.mkdtemp()
        data_file = os.path.join(tmp_dir, 'test.data')
        digest_file = os.path.join(tmp_dir, 'test.digest')

        # create the data file
        data_line = 'The quick brown fox jumps over the lazy dog. 0123456789\n'
        fd = open(data_file, 'w')
        for line in range(100):
            fd.write(data_line)
        fd.close()

        # create the digest file
        make_digest_file(data_file, digest_file)

        # get digest string for the data file
        digest = get_file_hexdigest(data_file)

        # check that digest is as expected, string
        msg = ("Digest string wrong, got '%s', expected '%s'"
               % (digest, expected_digest))
        self.assertTrue(expected_digest == digest, msg)

        # check that digest is as expected, file
        msg = ("Digest file wrong, got '%s', expected '%s'"
               % (digest, expected_digest))
        fd = open(digest_file, 'r')
        digest = fd.readline()
        fd.close()
        self.assertTrue(expected_digest == digest, msg)


    def test_file_length_function(self):
        '''Test that file_length() give 'correct' answer.'''

        # prepare test directory and filenames
        tmp_dir = tempfile.mkdtemp()
        test_file1 = os.path.join(tmp_dir, 'test.file1')
        test_file2 = os.path.join(tmp_dir, 'test.file2')
        test_file3 = os.path.join(tmp_dir, 'test.file3')
        test_file4 = os.path.join(tmp_dir, 'test.file4')

        # create files of known length
        fd = open(test_file1, 'w')      # 0 lines
        fd.close
        fd = open(test_file2, 'w')      # 5 lines, all '\n'
        for i in range(5):
            fd.write('\n')
        fd.close()
        fd = open(test_file3, 'w')      # 25 chars, no \n, 1 lines
        fd.write('no newline at end of line')
        fd.close()
        fd = open(test_file4, 'w')      # 1000 lines
        for i in range(1000):
            fd.write('The quick brown fox jumps over the lazy dog.\n')
        fd.close()

        # use file_length() to get and check lengths
        size1 = file_length(test_file1)
        msg = 'Expected file_length() to return 0, but got %d' % size1
        self.assertTrue(size1 == 0, msg)
        size2 = file_length(test_file2)
        msg = 'Expected file_length() to return 5, but got %d' % size2
        self.assertTrue(size2 == 5, msg)
        size3 = file_length(test_file3)
        msg = 'Expected file_length() to return 1, but got %d' % size3
        self.assertTrue(size3 == 1, msg)
        size4 = file_length(test_file4)
        msg = 'Expected file_length() to return 1000, but got %d' % size4
        self.assertTrue(size4 == 1000, msg)
        
        
        
    def test_get_revision_number(self):
        """test_get_revision_number
        
        Test that a revision number is returned.
        This should work both from a sandpit with access to Subversion
        and also in distributions where revision number has been stored
        explicitly in version.py
        """

        x = get_revision_number()

        try:
            y = x.split('.')
            assert int(y[0])
            assert int(y[1])
            assert int(y[2])
        except:
            assert int(x)

################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_system_tools, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)

