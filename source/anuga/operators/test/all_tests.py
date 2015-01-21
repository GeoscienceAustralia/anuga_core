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


#List files that should be excluded from the testing process.
#E.g. if they are known to fail and under development
exclude = [] 


def get_test_files(path):

    import sys

    files = os.listdir(path)

    #Check sub directories
    test_files = []
    for file in files:
        if os.path.isdir(file):
            sys.path.append(file)
            #print 'Recursing into', file
            test_files += get_test_files(path + os.sep + file)
        elif file[:5] == 'test_' and file[-2:] == 'py':
            #print 'Appending', file
            test_files.append(file)
        else:
            pass
    return test_files



def regressionTest():
    import sys, os, re, unittest
    path = os.path.split(sys.argv[0])[0] or os.getcwd()


    files = get_test_files(path)



    #test = re.compile('^test_[\w]*.py$', re.IGNORECASE)
    #files = filter(test.search, files)


    try:
        files.remove(__file__)  #Remove self from list (Ver 2.3. or later)
    except:
        files.remove('test_all.py')

    print 'Testing:'
    for file in files:
        print '  ' + file

    if globals().has_key('exclude'):
        for file in exclude:
            files.remove(file)
            print 'WARNING: File '+ file + ' excluded from testing'


    filenameToModuleName = lambda f: os.path.splitext(f)[0]
    #print "files",files 
    moduleNames = map(filenameToModuleName, files)
    #print "moduleNames",moduleNames
    modules = map(__import__, moduleNames)
    load = unittest.defaultTestLoader.loadTestsFromModule
    return unittest.TestSuite(map(load, modules))

################################################################################

if __name__ == '__main__':   
    suite = regressionTest()
    runner = unittest.TextTestRunner() #verbosity=2)
    runner.run(suite)
