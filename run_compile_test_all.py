# Compile and unit test Anuga
import os, sys

buildroot = os.getcwd()
execfile('compile_all.py')

os.chdir(buildroot)
print 'Buildroot is ', buildroot

os.chdir('source')
os.chdir('anuga')
    
print 'Changing to', os.getcwd() #This is now different from buildroot   

print "Test verbose "
os.system("python test_all.py v")

execfile('test_all.py')

    

