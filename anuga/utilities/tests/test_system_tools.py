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


from anuga.utilities.system_tools import *


class Test_system_tools(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_user_name(self):
        user = get_user_name()

        # print user
        assert isinstance(user, str), 'User name should be a string'

    def test_host_name(self):
        host = get_host_name()

        # print host
        assert isinstance(host, str), 'User name should be a string'

    def test_compute_checksum(self):
        """test_compute_checksum(self):

        Check that checksums on files are OK
        """

        from tempfile import mkstemp, mktemp

        # Generate a text file
        tmp_fd, tmp_name = mkstemp(suffix='.tmp', dir='.')
        fid = os.fdopen(tmp_fd, 'w+b')
        string = 'My temp file with textual content. AAAABBBBCCCC1234'
        binary_object = string.encode()  # Make it binary
        fid.write(binary_object)
        fid.close()

        # Have to apply the 64 bit fix here since we aren't comparing two
        # files, but rather a string and a file.
        ref_crc = safe_crc(binary_object)

        checksum = compute_checksum(tmp_name)
        assert checksum == ref_crc

        os.remove(tmp_name)

        # Binary file
        tmp_fd, tmp_name = mkstemp(suffix='.tmp', dir='.')
        fid = os.fdopen(tmp_fd, 'w+b')

        binary_object = b'My temp file with binary content. AAAABBBBCCCC1234'
        fid.write(binary_object)
        fid.close()

        ref_crc = safe_crc(binary_object)
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


    def test_get_pathname_from_package(self):
        """test_get_pathname_from_package(self):

        Check that correct pathname can be derived from package
        """

        path = get_pathname_from_package('anuga')
        assert path.endswith('anuga')


    def test_compute_checksum_real(self):
        """test_compute_checksum(self):

        Check that checksums on a png file is OK
        """

        path = get_pathname_from_package('anuga.utilities')

        filename = os.path.join(path, 'tests', 'data', 'crc_test_file.png')

        ref_crc = 1203293305  # Computed on Windows box
        checksum = compute_checksum(filename)

        msg = 'Computed checksum = %s, should have been %s'\
              % (checksum, ref_crc)
        assert checksum == ref_crc, msg
        # print checksum

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
        self.clean_line_helper('abc, ,,xyz,123', ',', [
                               'abc', '', 'xyz', '123'])

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
        for entry in range(MAX_ENTRIES):
            length = random.randint(1, MAX_CHARS)
            s = ''
            for c in range(length):
                s += chr(random.randint(A_INT, Z_INT))
            str_list.append(s)

        x = string_to_char(str_list)
        new_str_list = char_to_string(x)
        self.assertEqual(new_str_list, str_list)


    # special test - input list is ['']
    def test_string_to_char2(self):
        # generate a special list shown bad in load_mesh testing
        str_list = ['']

        x = string_to_char(str_list)
        new_str_list = char_to_string(x)

        self.assertEqual(new_str_list, str_list)


################################################################################
# Test the raw I/O to NetCDF files of string data encoded/decoded with
# string_to_char() and char_to_string().
################################################################################

# Note that the command num.array(string_to_char(l), num.character) gives
# rise to the following warning:
# DeprecationWarning: Converting `np.character` to a dtype is deprecated.
# The current result is `np.dtype(np.str_)` which is not strictly correct.
# Note that `np.character` is generally deprecated and 'S1' should be used.
# I was not able to find out why why searching, but using 'S1' is working.


    def helper_write_msh_file(self, filename, l):
        # open the NetCDF file
        fd = NetCDFFile(filename, netcdf_mode_w)
        fd.description = 'Test file - string arrays'

        # convert list of strings to num.array
        #al = num.array(string_to_char(l), num.character) # See note above
        al = num.array(string_to_char(l), 'S')

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

    def test_string_to_netcdf1(self):
        import random

        MAX_CHARS = 10
        MAX_ENTRIES = 10000

        A_INT = ord('a')
        Z_INT = ord('z')

        FILENAME = 'test.msh'

        # generate some random strings in a list, with guaranteed lengths
        str_list = ['x' * MAX_CHARS]        # make first maximum length
        for entry in range(MAX_ENTRIES):
            length = random.randint(1, MAX_CHARS)
            s = ''
            for c in range(length):
                s += chr(random.randint(A_INT, Z_INT))
            str_list.append(s)

        self.helper_write_msh_file(FILENAME, str_list)
        new_str_list = self.helper_read_msh_file(FILENAME)
        #print(str_list[:10])
        #print(new_str_list[:10])
        self.assertEqual(new_str_list, str_list)
        os.remove(FILENAME)

    # special test - list [''] to a NetCDF file
    def test_string_to_netcdf2(self):
        FILENAME = 'test.msh'

        # generate some random strings in a list, with guaranteed lengths
        str_list = ['']

        self.helper_write_msh_file(FILENAME, str_list)
        new_str_list = self.helper_read_msh_file(FILENAME)

        self.assertEqual(new_str_list, str_list)
        os.remove(FILENAME)

    def test_get_vars_in_expression(self):
        '''Test the 'get vars from expression' code.'''

        import warnings
        warnings.simplefilter('ignore', DeprecationWarning)

        def test_it(source, expected):
            result = get_vars_in_expression(source)
            result.sort()
            expected.sort()
            msg = ("Source: '%s'\nResult: %s\nExpected: %s"
                   % (source, str(result), str(expected)))
            self.assertEqual(result, expected, msg)

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

        warnings.simplefilter('default', DeprecationWarning)

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
        fd.close()
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
        This should work both from a sandpit with access to Git
        and also in distributions where revision number is returned as 0
        """

        x = get_revision_number()
        assert len(x) >= 0

    def test_get_revision_date(self):
        """test_get_revision_date

        Test that a revision date is returned.
        This should work both from a sandpit with access to Git
        and also in distributions where revision date is returned as 0
        """

        x = get_revision_date()
        assert len(str(x)) >= 0  # FIXME not sure how to test that this is a date (or the default string).

################################################################################

class Test_argparsing(unittest.TestCase):
    """Tests for anuga.utilities.argparsing.create_standard_parser (lines 16-54)."""

    def test_create_standard_parser(self):
        """create_standard_parser returns an ArgumentParser."""
        from anuga.utilities.argparsing import create_standard_parser
        parser = create_standard_parser()
        self.assertIsNotNone(parser)

    def test_parser_defaults(self):
        """Parser should have 'alg' default from parameters (line 30)."""
        from anuga.utilities.argparsing import create_standard_parser
        parser = create_standard_parser()
        args = parser.parse_args([])
        self.assertIsNotNone(args.alg)


class Test_data_audit_wrapper(unittest.TestCase):
    """Tests for anuga.utilities.data_audit_wrapper (lines 14-52)."""

    def test_module_import(self):
        """Importing the module covers module-level assignments (lines 14-41)."""
        import anuga.utilities.data_audit_wrapper as daw
        self.assertIsInstance(daw.extensions_to_ignore, list)

    def test_ip_verified_callable(self):
        """IP_verified function can be called on a directory (lines 45-52)."""
        import tempfile
        from anuga.utilities.data_audit_wrapper import IP_verified
        with tempfile.TemporaryDirectory() as d:
            result = IP_verified(d)
            # Just verify it returns without error; result may be True or False
            self.assertIsNotNone(result)


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_system_tools)
    runner = unittest.TextTestRunner()
    runner.run(suite)


# ---------------------------------------------------------------------------
# pytest functions: memory stats helpers
# ---------------------------------------------------------------------------

import pytest


@pytest.fixture(scope='module')
def simple_domain():
    import anuga
    return anuga.rectangular_cross_domain(4, 4)


def test_memory_stats_returns_mb_string():
    from anuga.utilities.system_tools import memory_stats
    result = memory_stats()
    assert isinstance(result, str)
    assert result.startswith('mem=')


def test_memory_stats_gb_branch():
    from anuga.utilities.system_tools import memory_stats
    from unittest.mock import patch
    with patch('anuga.utilities.system_tools._get_rss_mb', return_value=2048.0):
        result = memory_stats()
    assert 'GB' in result


def test_print_memory_stats(capsys):
    from anuga.utilities.system_tools import print_memory_stats
    print_memory_stats()
    out = capsys.readouterr().out
    assert 'mem=' in out


def test_quantity_memory_stats_returns_table(simple_domain):
    from anuga.utilities.system_tools import quantity_memory_stats
    result = quantity_memory_stats(simple_domain)
    assert isinstance(result, str)
    assert 'GRAND TOTAL' in result
    assert 'Quantity memory' in result


def test_print_quantity_memory_stats(simple_domain, capsys):
    from anuga.utilities.system_tools import print_quantity_memory_stats
    print_quantity_memory_stats(simple_domain)
    out = capsys.readouterr().out
    assert 'GRAND TOTAL' in out


def test_domain_memory_stats_returns_table(simple_domain):
    from anuga.utilities.system_tools import domain_memory_stats
    result = domain_memory_stats(simple_domain)
    assert isinstance(result, str)
    assert 'Domain memory' in result
    assert 'geometry' in result
    assert 'quantities' in result


def test_print_domain_memory_stats(simple_domain, capsys):
    from anuga.utilities.system_tools import print_domain_memory_stats
    print_domain_memory_stats(simple_domain)
    out = capsys.readouterr().out
    assert 'Domain memory' in out


def test_domain_struct_stats_no_gpu(simple_domain):
    from anuga.utilities.system_tools import domain_struct_stats
    result = domain_struct_stats(simple_domain)
    assert isinstance(result, str)
    assert 'C / GPU struct diagnostics' in result
    assert 'not initialised' in result


def test_print_domain_struct_stats(simple_domain, capsys):
    from anuga.utilities.system_tools import print_domain_struct_stats
    print_domain_struct_stats(simple_domain)
    out = capsys.readouterr().out
    assert 'C / GPU struct diagnostics' in out
