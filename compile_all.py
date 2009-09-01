import os
import time

buildroot = os.getcwd()

os.chdir('source')
os.chdir('anuga')


print 'Changing to', os.getcwd()        

#entries = listdir('.')

t0 = time.time()

# Attempt to compile all ANUGA extensions

os.chdir('utilities')
execfile('compile.py')

os.chdir('..')
os.chdir('abstract_2d_finite_volumes')
execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')

os.chdir('..')
os.chdir('advection')
execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')


os.chdir('..')
os.chdir('shallow_water')
execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')

os.chdir('..')
os.chdir('mesh_engine')
execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')

os.chdir(buildroot)    


# Attempt to compile Metis for use with anuga_parallel
os.chdir('source')
os.chdir('pymetis')

if sys.platform == 'win32':
    os.system('make for_win32')
else:
    if os.name == 'posix':
        if os.uname()[4] in ['x86_64', 'ia64']:
            os.system('make COPTIONS="-fPIC"')
        else:
            os.system('make')
    else:
        msg = 'Not sure how to complie Metis on platform: %s, %s' % (sys.platform, os.name)
        raise Exception, msg


    
print 'That took %.3fs' %(time.time() - t0)

if sys.platform == 'win32':
    raw_input('Press the RETURN key')
