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


#List files that should be excluded from the testing process.
#E.g. if they are known to fail and under development

exclude_files = []

#if sys.platform != 'win32':  
#    exclude_files.append('test_advection.py') #Weave doesn't work on Linux

# Exclude test_advection on all platforms for the time being. See ticket:205
exclude_files.append('test_advection.py') #Weave doesn't work on Linux


# Directories that should not be searched for test files.    
exclude_dirs = ['pypar_dist', #Special requirements
                'props', 'wcprops', 'prop-base', 'text-base', '.svn', #Svn
                'tmp']


print "The following directories will be skipped over;"
for dir in exclude_dirs:
    print dir
print ""

def get_test_files(path):


    try:
        files = os.listdir(path)
    except:
        return []

    #Check sub directories
    test_files = []

    #Exclude svn admin dirs
    files = [x for x in files if x not in exclude_dirs]
    path_files = []
    for file in files:

        absolute_filename = path + os.sep + file

        #sys.path.append('pmesh')
        if os.path.isdir(absolute_filename):
            sys.path.append(file) #FIXME: May cause name conflicts between pyvolution\mesh.py and pmesh\mesh.py on some systems
            path_files.append(file)
            print  file + ',', 
            more_test_files, more_path_files =\
                             get_test_files(absolute_filename)
            
            test_files += more_test_files
            path_files += more_path_files
        elif file.startswith('test_') and file.endswith('.py'):
            test_files.append(file)
        else:
            pass
        
    return test_files, path_files



def regressionTest(test_verbose=False):
    path = os.getcwd()
    print 'Recursing into;'
    test_files, path_files = get_test_files(path)

    files = [x for x in test_files if not x == 'test_all.py']

    files.sort() # Ensure same order on all platforms
    
    print
    print
    print 'Testing path %s:' %('...'+path[-50:])
    print
    print 'Files tested;'
    #print_files = []
    for file in files:
        #print_files += file + ' '
        print file + ',',
    print
    print
    if globals().has_key('exclude_files'):
        for file in exclude_files:
            print 'WARNING: File '+ file + ' to be excluded from testing'
            try:    
                files.remove(file)
            except ValueError, e:
                msg = 'File "%s" was not found in test suite.\n' %file
                msg += 'Original error is "%s"\n' %e
                msg += 'Perhaps it should be removed from exclude list?' 
                raise Exception, msg

    filenameToModuleName = lambda f: os.path.splitext(f)[0]
    moduleNames = map(filenameToModuleName, files)
    modules = map(__import__, moduleNames)
    
    # Fix up the system path
    for file in path_files:
        sys.path.remove(file)
        
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
                        # Calls class method set_verbose in the test case classes
                        # print 'Tests', tests._tests[0]
                        # print 'Type', type(tests._tests[0])                        
                        tests._tests[0].set_verbose()
                    except:
                        pass # No all classes have set_verbose
    return unittest.TestSuite(testCaseClasses)

def check_anuga_import():
    try:
        # importing something that loads quickly
        import anuga.utilities.anuga_exceptions
    except ImportError:
        print "Python cannot import ANUGA module."
        print "Check you have followed all steps of its installation."
        import sys; sys.exit() 

    
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
    runner = unittest.TextTestRunner() #verbosity=2
    runner.run(suite)
    
    # Cleaning up
    if len(sys.argv) > 1 and sys.argv[1][0].upper() == 'V':
        sys.stdout = saveout 
        #fid.close() # This was causing an error in windows
        #os.remove(filename)

