# Compile and test Anuga
import os, sys

execfile('compile_all.py')

buildroot = os.getcwd()

os.chdir('source')
os.chdir('anuga')
    
print 'Changing to', os.getcwd() #This is now different from buildroot   

print "Test verbose "
os.system("python test_all.py v")

execfile('test_all.py')

    

