"""Regression testing framework
This module will search for scripts in the same directory named
test_*.py.  Each such script should be a test suite that tests a
module through PyUnit. This script will aggregate all
found test suites into one big test suite and run them all at once.

Usage: test_all.py [<options>]

where <options> is zero or more of:
      -q         be quiet, minimal output to screen, write to ./.temp
      --quiet    same as above
      -v         be verbose to screen, more '-v's makes more verbose
      --verbose  same as above
You may do things like "test_all.py -q -v -v".
"""

# Author: Mark Pilgrim
# Modified by Ole Nielsen

import unittest
import os
import sys
import tempfile
import time
import anuga.config as config
from anuga.utilities.terminal_width import terminal_width
import anuga.utilities.system_tools as aust


#List files that should be excluded from the testing process.
#E.g. if they are known to fail and under development
exclude_files = []

# Directories that should not be searched for test files.
exclude_dirs = ['shallow_water_balanced' , # Special requirements
                '.svn',          # subversion
                'props', 'wcprops', 'prop-base', 'text-base', 'tmp']

# name of file to capture stdout in
CaptureFilename = '.temp'


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


def get_unittestfiles(path):
    """ Get 'test_*.py' files and paths to directories.
    param path Path to directory to start walk in.
    @return A tuple (<files>, <dirs>).
    Don't include any files in and below forbidden directories.

    """
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
    
    common_dir = os.path.dirname(path)
    common_dir_length = len(common_dir)
    #print common_dir, common_dir_length
    
    # get the terminal width
    term_width = terminal_width()

    # explain what we are doing
    print
    print "The following directories will be skipped over:"
    exclude_dirs.sort()
    list_names(exclude_dirs, page_width=term_width)

    # get all test_*.py and enclosing directories
    test_files, path_files = get_unittestfiles(path)
    #print path_files
    path_files.sort()

    files = [x for x in test_files if not x == 'test_all.py']
    files.sort()        # Ensure same order on all platforms

    def my_filter(pathname):
        
        return pathname[common_dir_length+1:]
    
    print
    print 'Paths searched:'
    list_names(path_files, my_filter, page_width=term_width)

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
    from optparse import OptionParser

    import os
    for file in os.listdir('.'):
        if file.endswith('.sww') or\
                file.endswith('.msh') or\
                file.endswith('.csv') or\
                file.endswith('.asc') or\
                file.endswith('.prj') or\
                file.endswith('.tsh') or\
                file.endswith('.sts') or\
                file.endswith('.tms') or\
                file.endswith('.pickle'):
            try:
                os.remove(file)
            except:
                pass

    check_anuga_import()

    # check the commandline params
    usage = "usage: %prog [-h|-q|-v]"
    parser = OptionParser(usage)
    parser.add_option("-v", "--verbose",
                      action="count", dest="verbosity", default=1,
                      help='make the ouput even more verbose')
    parser.add_option("-q", "--quiet",
                      action="store_true", dest="quiet",
                      help="capture statistics output in file '%s'"
                           % CaptureFilename)
    (options, args) = parser.parse_args()

    if len(args) > 0:
        parser.error("test_all.py doesn't take any parameters.  "
                     "Do 'test_all.py -h' to see the help.")

    if options.quiet:
        saveout = sys.stdout
        filename = CaptureFilename
        fid = open(filename, 'w')
        sys.stdout = fid

    # run the tests
    suite = regressionTest(options.quiet)
    runner = unittest.TextTestRunner(verbosity=options.verbosity)
    runner.run(suite)

    # timestamp at the end
    timestamp = time.asctime()
    major_revision = config.major_revision
    minor_revision = aust.get_revision_number()
    print '\nFinished at %s, version %s %s' % (timestamp, major_revision, minor_revision)

    # Cleaning up
    if options.verbosity < 1:
        sys.stdout = saveout
        #fid.close() # This was causing an error in windows
        #os.remove(filename)
    


    import os
    for file in os.listdir('.'):
        if file.endswith('.sww') or\
                file.endswith('.msh') or\
                file.endswith('.csv') or\
                file.endswith('.asc') or\
                file.endswith('.prj') or\
                file.endswith('.tsh') or\
                file.endswith('.sts') or\
                file.endswith('.tms') or\
                file.endswith('.pickle'):
            try:
                os.remove(file)
            except Exception as inst:
                pass
            
