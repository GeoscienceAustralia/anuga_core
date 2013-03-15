from anuga.utilities.data_audit_wrapper import IP_verified
from tempfile import mktemp

import os

buildroot = os.getcwd()


os.chdir('source')
os.chdir('anuga_parallel')
print
print '===================== anuga_parallel tests =========================='
print 'Changing to', os.getcwd()
execfile('test_all.py')


# Temporary bail out
import sys; sys.exit() 






    
