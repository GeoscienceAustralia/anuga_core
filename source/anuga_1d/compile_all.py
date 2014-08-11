import os

BUILDROOT = os.getcwd()
print 'Changing to', os.getcwd()

#Attempt to compile all extensions

os.chdir('utilities')
execfile('compile.py')

os.chdir('..')
os.chdir('base')
execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')

os.chdir('..')
os.chdir('channel')
execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')


#os.chdir('..')
#os.chdir('sww-sudi')
#execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')
#
#os.chdir('..')
#os.chdir('avalanche-sudi')
#execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')

## os.chdir('..')
## os.chdir('pipe')
## execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')

## os.chdir('..')
## os.chdir('sww')
## execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')

os.chdir('..')
os.chdir('sqpipe')
execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')


os.chdir(BUILDROOT)    
#execfile('test_all.py')
    
#if sys.platform == 'win32':
#    raw_input('Press the RETURN key')
