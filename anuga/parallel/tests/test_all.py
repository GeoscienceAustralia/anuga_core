"""Regression testing framework
This module will search for scripts in the same directory named
test_*.py.  Each such script should be a test suite that tests a
module through PyUnit. This script will aggregate all
found test suites into one big test suite and run them all at once.
"""

# Author: Mark Pilgrim
# Modified by Ole Nielsen

import unittest
import os
import sys
import tempfile
import time
import anuga.utilities.system_tools as aust
from anuga.utilities.terminal_width import terminal_width





#List files that should be excluded from the testing process.
#E.g. if they are known to fail and under development
exclude_files = []

# Directories that should not be searched for test files.
exclude_dirs = ['documentation',               # Special requirements
                '.svn',                              # subversion
                'props', 'wcprops', 'prop-base', 'text-base', 'tmp']


##
# @brief List a string sequence on the screen in columns.
# @param names Sequence of strings to list.
# @param func Function to apply to each string in sequence.
# @param col_width Force columns to this width (default calculated).
# @param page_width Set displayable page width to this (default 132).
def list_names(names, func=None, col_width=None, page_width=None):
    # set defaults
    p_width = page_width - 1            # set page width
    if p_width is None:
        p_width = 132                   # default page width

    c_width = col_width                 # set column width
    if c_width is None:
        c_width = 0
        for name in names:
            if func:
                name = func(name)
            c_width = max(c_width, len(name))
    c_width += 2                        # 2 column padding

    # calculate number of columns allowed
    max_columns = int(p_width / c_width)

    # print columns
    column = 0
    for name in names:
        if func:
            name = func(name)
        print '%-*s' % (c_width-1, name),
        column += 1
        if column >= max_columns:
            column = 0
            print

    # if last line not finished, end it here
    if column > 0:
        print


##
# @brief Get 'test_*.py' files and paths to directories.
# @param path Path to directory to start walk in.
# @return A tuple (<files>, <dirs>).
# @note Don't include any files in and below forbidden directories.
def get_unittestfiles(path):
    walk = os.walk(path)

    test_files = []
    path_files = []

    for (dirpath, dirnames, filenames) in walk:
        # exclude forbidden directories
        for e_dir in exclude_dirs:
            try:
                dirnames.remove(e_dir)
            except ValueError:
                pass

        # check for test_*.py files
        for filename in filenames:
            if filename.startswith('test_') and filename.endswith('.py'):
                test_files.append(filename)
                if dirpath not in path_files:
                    path_files.append(dirpath)

    return test_files, path_files


def regressionTest(test_verbose=False):
    # start off with where we are
    path = os.getcwd()
    print
    print 'Testing path: %s' % path

    # get the terminal width
    term_width = terminal_width()

    # explain what we are doing
    print
    print "The following directories will be skipped over:"
    exclude_dirs.sort()
    list_names(exclude_dirs, page_width=term_width)

    # get all test_*.py and enclosing directories
    test_files, path_files = get_unittestfiles(path)
    path_files.sort()

    files = [x for x in test_files if not x == 'test_all.py']
    files.sort()        # Ensure same order on all platforms

    print
    print 'Paths searched:'
    list_names(path_files, os.path.basename, page_width=term_width)

    print    
    print 'Files tested:'
    list_names(files, page_width=term_width)
    print

    # update system path with found paths
    for path in path_files:
        sys.path.append(path)
   
    # exclude files that we can't handle 
    for file in exclude_files:
        print 'WARNING: File '+ file + ' to be excluded from testing'
        try:
            files.remove(file)
        except ValueError, e:
            msg = 'File "%s" was not found in test suite.\n' % file
            msg += 'Original error is "%s"\n' % e
            msg += 'Perhaps it should be removed from exclude list?'
            raise Exception, msg

    # import all test_*.py files
    # NOTE: This implies that test_*.py files MUST HAVE UNIQUE NAMES!
    filenameToModuleName = lambda f: os.path.splitext(f)[0]
    moduleNames = map(filenameToModuleName, files)
    modules = map(__import__, moduleNames)

    # Fix up the system path
    for file in path_files:
        sys.path.remove(file)

    # bundle up all the tests
    load = unittest.defaultTestLoader.loadTestsFromModule
    testCaseClasses = map(load, modules)

    if test_verbose is True:
        # Test the code by setting verbose to True.
        # The test cases have to be set up for this to work.
        # See test data manager for an example.
        for test_suite in testCaseClasses:
            for tests in test_suite._tests:
                # tests is of class TestSuite
                if len(tests._tests) > 1:
                    # these are the test functions
                    try:
                        # Calls class method set_verbose in test case classes
                        tests._tests[0].set_verbose()
                    except:
                        pass                # No all classes have set_verbose

    return unittest.TestSuite(testCaseClasses)


##
# @brief Check that the environment is sane.
# @note Stops here if there is an error.
def check_anuga_import():
    try:
        # importing something that loads quickly
        import anuga.anuga_exceptions
    except ImportError:
        print "Python cannot import ANUGA module."
        print "Check you have followed all steps of its installation."
        import sys
        sys.exit()


if __name__ == '__main__':
    check_anuga_import()

    if len(sys.argv) > 1 and sys.argv[1][0].upper() == 'V':
        test_verbose = True
        saveout = sys.stdout
        filename = ".temp"
        fid = open(filename, 'w')
        sys.stdout = fid
    else:
        test_verbose = False
    suite = regressionTest(test_verbose)
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)

    # timestamp at the end
    timestamp = time.asctime()
    version = aust.get_revision_number()
    print
    print 'Finished at %s, version %s' % (timestamp, version)

    # Cleaning up
    if len(sys.argv) > 1 and sys.argv[1][0].upper() == 'V':
        sys.stdout = saveout
        #fid.close() # This was causing an error in windows
        #os.remove(filename)

    
    if sys.platform == 'win32':
        raw_input('Press the RETURN key')
