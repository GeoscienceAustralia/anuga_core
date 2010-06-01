import os

buildroot = os.getcwd()

#Attempt to compile all extensions

os.chdir('utilities')
execfile('compile.py')

os.chdir('..')
os.chdir('advection')
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
os.chdir('shallow_water_balanced')
execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')

os.chdir('..')
os.chdir('mesh_engine')
execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')

os.chdir(buildroot)    
#execfile('test_all.py')
    
if sys.platform == 'win32':
    raw_input('Press the RETURN key')
