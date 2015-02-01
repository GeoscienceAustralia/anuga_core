""" Script to compile all C extensions in ANUGA. """

import os
import subprocess
import sys

BUILDROOT = os.getcwd()

#Attempt to compile all extensions



os.chdir('utilities')
subprocess.call([sys.executable, 'compile.py', 'quad_tree.c'])
subprocess.call([sys.executable, 'compile.py', 'sparse_dok.c'])
subprocess.call([sys.executable, 'compile.py', 'sparse_csr.c'])
execfile('compile.py')

os.chdir('..')
os.chdir('advection')
execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')

os.chdir('..')
os.chdir('operators')
execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')

os.chdir('..')
os.chdir('file_conversion')
execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')

os.chdir('..')
os.chdir('geometry')
execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')

os.chdir('..')
os.chdir('structures')
execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')

os.chdir('..')
os.chdir('abstract_2d_finite_volumes')
execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')

os.chdir('..')
os.chdir('file')
execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')

os.chdir('..')
os.chdir('shallow_water')
execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')


os.chdir('..')
os.chdir('mesh_engine')
execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')

os.chdir('..')
os.chdir('fit_interpolate')
subprocess.call([sys.executable, '..' + os.sep + 'utilities' + os.sep + 'compile.py', 'rand48.c'])
subprocess.call([sys.executable, '..' + os.sep + 'utilities' + os.sep + 'compile.py', 'ptinpoly.c'])
execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')



#====================================================================
#os.chdir('..')
#os.chdir('utilities')
#try:
#    from anuga.utilities.system_tools  import store_version_info
#    store_version_info(verbose=True)
#    print
#    print "---------------------------------"
#    print 'Storing of version info succeeded'
#    print "---------------------------------"
#    print
#except:
#    print
#    print "----------------------------------------------------------------"
#    print 'Storage of version info failed (just means svn is not available)'
#    print "----------------------------------------------------------------"
#    print

os.chdir(BUILDROOT)