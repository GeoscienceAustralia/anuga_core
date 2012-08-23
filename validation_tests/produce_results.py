"""
Script to run all the produce_results scripts in the Tests/xxx/xxx/ directories
"""

import os
import anuga
import time
from anuga import indent

#--------------------------------
# Get Default values for
# algorithm parameters.
#--------------------------------
from parameters import alg
from parameters import cfl

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

try:
    Upper_dirs.remove('.svn')
except ValueError:
    pass
#print Upper_dirs
os.chdir('./Tests')

print 'Tests'
for dir in Upper_dirs:

    os.chdir(dir)

    print indent + dir
    #print 'Changing to', os.getcwd()
    Lower_dirs = os.listdir('.')
    try:
        Lower_dirs.remove('.svn')
    except ValueError:
        pass
    #print Lower_dirs
    for l_dir in Lower_dirs:
        os.chdir(l_dir)
        #print os.getcwd()
        print 2*indent + l_dir
        try:
            cmd = 'python produce_results.py'
            print 3*indent + 'Running: '+cmd
            os.system( cmd )
        except:
            print 3*indent + 'Failed running produce_results in '+os.getcwd()
            pass

        os.chdir('..')
        #print 'Changing to', os.getcwd()

    os.chdir('..')
    #print 'Changing to', os.getcwd()
    
os.chdir('..')
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


cmd = 'pdflatex -shell-escape -interaction=batchmode report.tex'
print cmd
import subprocess
subprocess.call( [cmd], shell=True )
cmd = 'mv report.pdf report_cfl_%s_alg_%s.pdf' % (str(cfl), str(alg))
print cmd
subprocess.call( [cmd] , shell=True )





