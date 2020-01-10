from __future__ import print_function
from future.utils import raise_
import os
import time
import sys
import subprocess

buildroot = os.getcwd()
t0 = time.time()

#--------------------------------------------------
# Compiling anuga code 
#--------------------------------------------------
#os.chdir('source')

#print 'Changing to', os.getcwd()        

if sys.platform == 'win32':
    cmd = 'python setup.py build --compiler=mingw32  install --install-lib=source '
else:
    cmd = 'python setup.py install --install-lib=source '
print(cmd)
err = os.system(cmd)
if err != 0:
    msg = 'Could not compile anuga '
    msg += 'on platform %s, %s\n' % (sys.platform, os.name)
    raise_(Exception, msg)
else:
    print(50*"-")
    print()
    msg = 'Compiled anuga successfully.'
    print(msg)

print()        
print('That took %.3fs' %(time.time() - t0))




