"""
Script to run all the produce_results scripts in the
validation_tests/xxx/xxx/ directories
"""

import os
import time

import anuga
from anuga import indent
#from anuga.validation_utilities.parameters import alg
#from anuga.validation_utilities.parameters import cfl


args = anuga.get_args()
alg = args.alg
np = args.np
verbose = args.verbose

#---------------------------------
# Get the current svn revision
#---------------------------------
timestamp = time.asctime()
major_revision = anuga.get_version()
try:
    # This fails if using git for version control
    minor_revision = anuga.get_revision_number()
except:
    try:
        # This works when using git on unix
        minor_revision = os.popen("git show-ref --head -s | head -n1").read().strip()
    except:
        # This is a fallback position
        minor_revision = 'unknown'


#----------------------------------
# Now it is ok to create the latex 
# macro file with run parameters
#
# FIXME: THis is a little dangerous as
# this is changed before all the tests
# are run. 
#----------------------------------

f = open('saved_parameters.tex', 'w')
#f.write('\\newcommand{\\cfl}{\\UScore{%s}}\n' % str(cfl))
f.write('\\newcommand{\\alg}{\\UScore{%s}}\n' % str(alg))
f.write('\\newcommand{\\majorR}{\\UScore{%s}}\n' % str(major_revision))
f.write('\\newcommand{\\minorR}{\\UScore{%s}}\n' % str(minor_revision))
f.write('\\newcommand{\\timeR}{{%s}}\n' % str(timestamp))

f.close()

#---------------------------------
# Run the tests
#---------------------------------
os.chdir('..')
buildroot = os.getcwd()

Upper_dirs = os.listdir('.')
dir = '.'
Upper_dirs = [name for name in os.listdir(dir) if os.path.isdir(os.path.join(dir, name))]

try:
    Upper_dirs.remove('.svn')
except ValueError:
    pass

try:
    Upper_dirs.remove('reports')
except ValueError:
    pass

try:
    Upper_dirs.remove('case_studies')
except ValueError:
    pass

#print Upper_dirs
#os.chdir('./Tests')

#print 'Tests'
print(Upper_dirs)

time_total = 0.0
test_number = 1
for dir in Upper_dirs:

    os.chdir(dir)

    print(72 * '=')
    print('Directory: ' + dir)
    print(72 * '=')
    
    #print 'Changing to', os.getcwd()
    dir = '.'
    Lower_dirs =  [name for name in os.listdir(dir) if os.path.isdir(os.path.join(dir, name))]
    try:
        Lower_dirs.remove('.svn')
    except ValueError:
        pass
    #print Lower_dirs




    for l_dir in Lower_dirs:
        os.chdir(l_dir)
        #print os.getcwd()
        print(60 * '=')
        print('Subdirectory %g: '% (test_number)  + l_dir)
        test_number += 1
        print(60 * '=')
        try:
            t0 = time.time()
            if verbose:
                cmd = 'python produce_results.py -alg %s -np %s -v '% (str(alg),str(np))
            else:
                cmd = 'python produce_results.py -alg %s -np %s '% (str(alg),str(np))
            print(2 * indent + 'Running: ' + cmd)
            os.system(cmd)
            t1 = time.time() - t0
            time_total += t1
            print(2 * indent + 'That took ' + str(t1) + ' secs')
        except:
            print(2 * indent + 'Failed running produce_results in ' + os.getcwd())
            pass

        os.chdir('..')
        #print 'Changing to', os.getcwd()

    os.chdir('..')
    #print 'Changing to', os.getcwd()
    
os.chdir(buildroot)

print(72 * '=')
print('That took ' + str(time_total) + ' secs')
print(72 * '=')


# go back to reports directory to typeset report
os.chdir('reports')


os.system('python validations_typeset_report.py')

import subprocess
cmd = 'mv validations_report.pdf validations_report_alg_%s.pdf' % (str(alg))
print(cmd)
subprocess.call([cmd], shell=True)





