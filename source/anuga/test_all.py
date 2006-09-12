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


#List files that should be excluded from the testing process.
#E.g. if they are known to fail and under development
exclude_files = ['test_version.py', #'test_least_squares.py',
                 'test_advection.py', # removing this test for a bit
                 ]
                 #'test_calculate_region.py', 'test_calculate_point.py']
                 #'test_init.py']

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
            print 'Recursing into', file
            more_test_files, more_path_files =get_test_files(absolute_filename)
            test_files += more_test_files
            path_files += more_path_files
        elif file.startswith('test_') and file.endswith('.py'):
            test_files.append(file)
        else:
            pass
    return test_files , path_files



def regressionTest():
    path = os.getcwd()
    test_files, path_files = get_test_files(path)
    files = [x for x in test_files if not x == 'test_all.py']

    print 'Testing path %s:' %('...'+path[-50:])
    for file in files:
        print '  ' + file
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
    return unittest.TestSuite(map(load, modules))

if __name__ == '__main__':
    #unittest.main(defaultTest='regressionTest')

    suite = regressionTest()
    runner = unittest.TextTestRunner() #verbosity=2
    runner.run(suite)
