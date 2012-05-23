"""
Script to run all the produce_results scripts in the Tests/xxx/xxx/ directories
"""

import os

buildroot = os.getcwd()

Upper_dirs = os.listdir('./Tests')
#print Upper_dirs
try:
    Upper_dirs.remove('.svn')
except ValueError:
    pass
#print Upper_dirs
os.chdir('./Tests')

for dir in Upper_dirs:

    os.chdir(dir)
    #print 'Changing to', os.getcwd()
    Lower_dirs = os.listdir('.')
    try:
        Lower_dirs.remove('.svn')
    except ValueError:
        pass
    #print Lower_dirs
    for l_dir in Lower_dirs:
        os.chdir(l_dir)
        print 'Changing to', os.getcwd()
        try:
            os.system('python produce_results.py')
        except:
            print 'Failed running produce_results in '+os.getcwd()
            pass

        os.chdir('..')
        #print 'Changing to', os.getcwd()

    os.chdir('..')
    #print 'Changing to', os.getcwd()
        


