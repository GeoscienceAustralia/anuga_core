from __future__ import print_function
from past.builtins import execfile
import os


os.chdir('validation_tests')
print()
print(20*'=' + ' anuga automated validation tests ' + 20*'=')
print('Changing to', os.getcwd()) # This is now different from buildroot   
execfile('run_auto_validation_tests.py')

