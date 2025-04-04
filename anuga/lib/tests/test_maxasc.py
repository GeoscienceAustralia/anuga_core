#!/usr/bin/env python



from anuga.lib.maxasc import *

import anuga.utilities.system_tools as aust
from functools import reduce


import unittest

import sys
import os
import re
import glob

HEADER_SIZE = 6

# pattern string used to split multimax data
SpacesPatternString = ' +'

# generate 're' pattern for 'any number of spaces'
SpacesPattern = re.compile(SpacesPatternString)


def FilesEqual(file1, file2):
    """
    Compare two ASC files for equality.
    """
    
    def alltrue(a, b):
        return a and b

    def do_list(prefix, l):
        print(prefix, '#'*100)
        for (i, x) in enumerate(l):
            print('%05d: %s\t' % (i, str(x)), end=' ')
            if (i+1) % 5 == 0:
                print()
        print()
        print()
    
    # get both files into memory
    fd = open(file1, 'r')
    data1 = fd.readlines()
    fd.close()
    fd = open(file2, 'r')
    data2 = fd.readlines()
    fd.close()

    # check we have same number of lines in each file
    if len(data1) != len(data2):
        print('# lines differs: len(data1)=%d, len(data2)=%d' % (len(data1), len(data2)))
        return False
    
    # read header lines, check identical
    for i in range(HEADER_SIZE):
        if data1[i] != data2[i]:
            print('headers differ:')
            print(data1[i])
            print(data2[i])
            return False

    # read data lines, check same *values*
    for line_num in range(HEADER_SIZE, len(data1)):
        d1 = SpacesPattern.split(data1[line_num].strip())
        d2 = SpacesPattern.split(data2[line_num].strip())

        if len(d1) != len(d2):
            print('# columns differs, len(d1)=%d, len(d2)=%d on line %d' % (len(d1), len(d2), line_num))
            return False

        fd1 = [float(value) for value in d1]
        fd2 = [float(value) for value in d2]

        vec = list(map(lambda a,b: a==b, fd1, fd2))
        
        if not reduce(lambda a,b: a and b, vec):
            print('line number = %d (out of %d)' % (line_num, len(data1)))
            do_list('fd1', fd1)
            do_list('fd2', fd2)
            do_list('vec', vec)
            return False
    
    return True


class Test_MaxAsc(unittest.TestCase):
    def tearDown(self):
        # delete all output files
        files = glob.glob('*.out.asc')
        for file in files:
            os.remove(file)

    def test_unequal_lines(self):
        lib_tests_path = os.path.join(aust.get_pathname_from_package('anuga.lib'), 'tests')
        in_file = os.path.join(lib_tests_path, 'test1.asc')
        expected_file = os.path.join(lib_tests_path, 'test1_bad_num_lines.asc')

        self.assertRaises(RuntimeError, MaxAsc,
                              'test.out.asc',
                              [in_file, expected_file])

    def test_headers_differ(self):
        lib_tests_path = os.path.join(aust.get_pathname_from_package('anuga.lib'), 'tests')
        in_file = os.path.join(lib_tests_path, 'test1.asc')
        expected_file = os.path.join(lib_tests_path, 'test1_bad_num_lines.asc')

        self.assertRaises(RuntimeError, MaxAsc,
                              'test.out.asc',
                              [in_file, expected_file])

    def test_wrong_num_columns(self):
        lib_tests_path = os.path.join(aust.get_pathname_from_package('anuga.lib'), 'tests')
        in_file = os.path.join(lib_tests_path, 'test1.asc')
        expected_file = os.path.join(lib_tests_path, 'test1_wrong_num_columns.asc')

        self.assertRaises(RuntimeError, MaxAsc,
                              'test.out.asc',
                              [in_file, expected_file])

    def test_same_input_equals_output1(self):
        lib_tests_path = os.path.join(aust.get_pathname_from_package('anuga.lib'), 'tests')
        in_file = os.path.join(lib_tests_path, 'test1.asc')

        MaxAsc('test1.out.asc', [in_file])
        self.assertTrue(FilesEqual('test1.out.asc', in_file))


    def test_same_input_equals_outputN(self):
        lib_tests_path = os.path.join(aust.get_pathname_from_package('anuga.lib'), 'tests')
        in_file = os.path.join(lib_tests_path, 'test1.asc')

        MaxAsc('test1.out.asc', [in_file] * 30)
        self.assertTrue(FilesEqual('test1.out.asc', in_file))

    def test_different_input2(self):
        lib_tests_path = os.path.join(aust.get_pathname_from_package('anuga.lib'), 'tests')
        in_file = os.path.join(lib_tests_path, 'test1.asc')
        in_file2 = os.path.join(lib_tests_path, 'test2.asc')
        expected_file = os.path.join(lib_tests_path, 'test2.expected.asc')

        MaxAsc('test2.out.asc', [in_file, in_file2])
        self.assertTrue(FilesEqual('test2.out.asc', expected_file))

    def test_different_input3(self):
        lib_tests_path = os.path.join(aust.get_pathname_from_package('anuga.lib'), 'tests')
        in_file = os.path.join(lib_tests_path, 'test1.asc')
        in_file2 = os.path.join(lib_tests_path, 'test2.asc')
        in_file3 = os.path.join(lib_tests_path, 'test3.asc')
        expected_file = os.path.join(lib_tests_path, 'test3.expected.asc')

        MaxAsc('test3.out.asc', [in_file, in_file2, in_file3])
        self.assertTrue(FilesEqual('test3.out.asc', expected_file))

if __name__ == '__main__':
    suite = unittest.makeSuite(Test_MaxAsc,'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
