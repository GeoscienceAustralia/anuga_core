"""compile.py - compile Python C-extension

   Commandline usage: 
     python compile.py <filename>
 
"""     
from __future__ import print_function
from past.builtins import execfile
import os

buildroot = os.getcwd()
os.chdir('..')


print('Changing to', os.getcwd())        

#entries = listdir('.')

#Attempt to compile mesh_engine extensions


os.chdir('mesh_engine')
execfile('..' + os.sep + 'utilities' + os.sep + 'compile.py')

os.chdir(buildroot)   
    
