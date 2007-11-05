import os

buildroot = os.getcwd()

os.chdir('source')
os.chdir('anuga')
    
print 'Changing to', os.getcwd() #This is now different from buildroot   

execfile('test_all.py')

    
