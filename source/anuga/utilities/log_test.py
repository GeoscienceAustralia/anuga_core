#!/usr/bin/env python

import os
import sys
import unittest
import logging
import log

LOGFILE_NAME = 'test.log'
STDOUT_LOG_NAME = 'stdout.log'


class Test_Log(unittest.TestCase):
    def setUp(self):
        # just in case
        if os.path.exists(LOGFILE_NAME):
            os.remove(LOGFILE_NAME)
        if os.path.exists(STDOUT_LOG_NAME):
            os.remove(STDOUT_LOG_NAME)

        # set log module logfile name
        log.log_filename = LOGFILE_NAME

        # set logging levels for this test
        log.console_logging_level = logging.INFO
        log.log_logging_level = logging.DEBUG

    def tearDown(self):
        if os.path.exists(LOGFILE_NAME):
            os.remove(LOGFILE_NAME)
        if os.path.exists(STDOUT_LOG_NAME):
            os.remove(STDOUT_LOG_NAME)
        pass

    ##
    # @brief Test the logging routines.
    # @note HAVE ONE TEST CASE ONLY! Multiple tests would concatenate
    #       multiple test output in one log file.
    def test_simple(self):
        '''Check that logging works in simple case.'''

        # vvvvv  WARNING - DO NOT REFORMAT THESE LINES!  vvvvv
        log_expect = '''2009-03-23 12:16:53,487 CRITICAL                       log:0   |Logfile is 'test.log' with logging level of DEBUG, console logging level is INFO
2009-03-23 12:16:53,488 DEBUG                     test_log:37  |test at level DEBUG
2009-03-23 12:35:26,107 INFO                      test_log:97  |Resource usage: memory=0.0 resident=0.0 stacksize=0.0
2009-03-23 12:16:53,488 INFO                      test_log:38  |test at level INFO'''

        stdout_expect = '''Logfile is 'test.log' with logging level of DEBUG, console logging level is INFO
Resource usage: memory=0.0 resident=0.0 stacksize=0.0
test at level INFO'''
        # ^^^^^  WARNING - DO NOT REFORMAT THESE LINES!  ^^^^^

        # capture stdout to a file
        save_stdout = sys.stdout
        save_stderr = sys.stderr
        sys.stdout = sys.stderr = open(STDOUT_LOG_NAME, 'w')

        # do some logging
        log.debug('test at level DEBUG')
        log.resource_usage(logging.INFO)
        log.info('test at level INFO')

        # put stdout/stderr back to normal
        sys.stderr = save_stderr
        sys.stdout = save_stdout

        # check logfile is as expected
        fd = open(LOGFILE_NAME, 'r')
        lines = fd.readlines()
        fd.close()
        result = []
        for line in lines:
            l = line.strip('\n')
            ndx = l.index('|')
            result.append([l.split()[2], l[ndx:]])
        expected = []
        for line in log_expect.split('\n'):
            ndx = line.index('|')
            expected.append([line.split()[2], line[ndx:]])
        self.failUnlessEqual(result, expected)

        # check that captured stdout is as expected
        fd = open(STDOUT_LOG_NAME, 'r')
        lines = fd.readlines()
        fd.close()
        result = []
        for line in lines:
            l = line.strip('\n')
            result.append(l)
        expected = []
        for line in stdout_expect.split('\n'):
            expected.append(line)
        self.failUnlessEqual(result, expected)

################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Log, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)
