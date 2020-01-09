from __future__ import print_function
import os

buildroot = os.getcwd()


os.chdir('anuga')
print()
print('======================= anuga tests =================================')    
print('Changing to', os.getcwd()) # This is now different from buildroot

os.system('python test_all.py')

os.chdir(buildroot)

