import os

buildroot = os.getcwd()

#Attempt to compile all extensions

os.chdir('utilities')
execfile('compile.py')

os.chdir('..')
os.chdir('pyvolution')
execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')

os.chdir('..')
os.chdir('mesh_engine')
execfile('compile.py')

os.chdir(buildroot)    
#execfile('test_all.py')
    
