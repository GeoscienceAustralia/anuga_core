import os
import time
import sys
import subprocess

buildroot = os.getcwd()

os.chdir('source')
os.chdir('anuga')


print 'Changing to', os.getcwd()        

#entries = listdir('.')

t0 = time.time()

# Attempt to compile all ANUGA extensions

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


os.chdir(buildroot)    



print        
print 'That took %.3fs' %(time.time() - t0)



if sys.platform == 'win32':
    raw_input('Press the RETURN key')
