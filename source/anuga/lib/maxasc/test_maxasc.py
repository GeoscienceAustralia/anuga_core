#!/usr/bin/env python

from maxasc import *

import exceptions
class TestError(exceptions.Exception): pass
import unittest

import sys
import re

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
        print prefix, '#'*100
        for (i, x) in enumerate(l):
            print '%05d: %s\t' % (i, str(x)),
            if (i+1) % 5 == 0:
                print
        print
        print
    
    # get both files into memory
    fd = open(file1, 'r')
    data1 = fd.readlines()
    fd.close()
    fd = open(file2, 'r')
    data2 = fd.readlines()
    fd.close()

    # check we have same number of lines in each file
    if len(data1) != len(data2):
        print '# lines differs: len(data1)=%d, len(data2)=%d' % (len(data1), len(data2))
        return False
    
    # read header lines, check identical
    for i in range(HEADER_SIZE):
        if data1[i] != data2[i]:
            print 'headers differ:'
            print data1[i]
            print data2[i]
            return False

    # read data lines, check same *values*
    for line_num in range(HEADER_SIZE, len(data1)):
        d1 = SpacesPattern.split(data1[line_num].strip())
        d2 = SpacesPattern.split(data2[line_num].strip())

        if len(d1) != len(d2):
            print '# columns differs, len(d1)=%d, len(d2)=%d on line %d' % (len(d1), len(d2), line_num)
            return False

        fd1 = [float(value) for value in d1]
        fd2 = [float(value) for value in d2]

        vec = map(lambda a,b: a==b, fd1, fd2)
        
        if not reduce(lambda a,b: a and b, vec):
            print 'line number = %d (out of %d)' % (line_num, len(data1))
            do_list('fd1', fd1)
            do_list('fd2', fd2)
            do_list('vec', vec)
            return False
    
    return True


class Test_MaxAsc(unittest.TestCase):
    def test_unequal_lines(self):
        self.failUnlessRaises(RuntimeError, MaxAsc,
                              'test.out.asc',
                              ['test1.asc', 'test1_bad_num_lines.asc'])

    def test_headers_differ(self):
        self.failUnlessRaises(RuntimeError, MaxAsc,
                              'test.out.asc',
                              ['test1.asc', 'test1_bad_num_lines.asc'])

    def test_wrong_num_columns(self):
        self.failUnlessRaises(RuntimeError, MaxAsc,
                              'test.out.asc',
                              ['test1.asc', 'test1_wrong_num_columns.asc'])

    def test_same_input_equals_output1(self):
        MaxAsc('test1.out.asc', ['test1.asc'])
        self.failUnless(FilesEqual('test1.out.asc', 'test1.asc'))

    def test_same_input_equals_bigA(self):
        MaxAsc('perth.out.asc', ['perthAll_stage_250m.asc'])
        self.failUnless(FilesEqual('perth.out.asc', 'perthAll_stage_250m.asc'))

    def test_same_input_equals_bigB(self):
        MaxAsc('perth.out.asc', ['perthAll_stage_250m_all.asc'])
        self.failUnless(FilesEqual('perth.out.asc', 'perthAll_stage_250m_all.asc'))

    def test_same_input_equals_bigC(self):
        MaxAsc('perth.out.asc', ['perthAll_stage_original.asc'])
        self.failUnless(FilesEqual('perth.out.asc', 'perthAll_stage_original.asc'))

    def test_same_input_equals_big3(self):
        MaxAsc('perth.out.asc', ['perthAll_stage_250m.asc',
                                 'perthAll_stage_250m.asc',
                                 'perthAll_stage_250m.asc'])
        self.failUnless(FilesEqual('perth.out.asc', 'perthAll_stage_250m.asc'))

    def test_same_input_equals_outputN(self):
        MaxAsc('test1.out.asc', ['test1.asc'] * 30)
        self.failUnless(FilesEqual('test1.out.asc', 'test1.asc'))

    def test_different_input2(self):
        MaxAsc('test2.out.asc', ['test1.asc', 'test2.asc'])
        self.failUnless(FilesEqual('test2.out.asc', 'test2.expected.asc'))

    def test_different_input3(self):
        MaxAsc('test3.out.asc', ['test1.asc', 'test2.asc', 'test3.asc'])
        self.failUnless(FilesEqual('test3.out.asc', 'test3.expected.asc'))


suite = unittest.makeSuite(Test_MaxAsc,'test')
#suite = unittest.makeSuite(Test_MaxAsc,'test_same_input_equals_output1')
runner = unittest.TextTestRunner(verbosity=1)
runner.run(suite)
