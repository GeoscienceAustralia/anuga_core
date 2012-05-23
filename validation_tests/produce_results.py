"""
Script to run all the produce_results scripts in the Tests/xxx/xxx/ directories
"""

import os
import anuga.utilities.system_tools as anugast
import anuga
import time

#--------------------------------
# Setup Default values for basis
# algorithm parameters.
#--------------------------------
import argparse
parser = argparse.ArgumentParser(description='produce results')
parser.add_argument('-cfl', type=float, default=1.0,
                   help='cfl condition')
parser.add_argument('-alg', type=str, default = "1_5",
                   help='flow algorithm')
args = parser.parse_args()

cfl = args.cfl
alg = args.alg

#---------------------------------
# Get the current svn revision
#---------------------------------

timestamp = time.asctime()
major_revision = anuga.config.major_revision
minor_revision = anuga.utilities.system_tools.get_revision_number()


#---------------------------------
# Run the tests
#---------------------------------
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
            cmd = 'python produce_results.py -alg %s -cfl %s '% (alg,cfl)
            print 'Running: '+cmd
            os.system( cmd )
        except:
            print 'Failed running produce_results in '+os.getcwd()
            pass

        os.chdir('..')
        #print 'Changing to', os.getcwd()

    os.chdir('..')
    #print 'Changing to', os.getcwd()
        
#----------------------------------
# Now it is ok to create the latex 
# macro file with run parameters
#----------------------------------

f = open('saved_parameters.tex','w')
f.write('\\newcommand{\\cfl}{\\UScore{%s}}\n' % str(cfl))
f.write('\\newcommand{\\alg}{\\UScore{%s}}\n' % str(alg))
f.write('\\newcommand{\\majorR}{\\UScore{%s}}\n' % str(major_revision))
f.write('\\newcommand{\\minorR}{\\UScore{%s}}\n' % str(minor_revision))
f.write('\\newcommand{\\timeR}{{%s}}\n' % str(timestamp))

f.close()





