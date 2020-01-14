""" Script to compile all C extensions in ANUGA. """

from past.builtins import execfile
from builtins import input
import os

BUILDROOT = os.getcwd()

#Attempt to compile all extensions

execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')


os.chdir(BUILDROOT)    

    
if sys.platform == 'win32':
    input('Press the RETURN key')
