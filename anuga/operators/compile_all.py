""" Script to compile all C extensions in ANUGA. """

import os

BUILDROOT = os.getcwd()

#Attempt to compile all extensions

execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')


os.chdir(BUILDROOT)    

    
if sys.platform == 'win32':
    raw_input('Press the RETURN key')
