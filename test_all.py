import os

buildroot = os.getcwd()

os.chdir('source')
os.chdir('anuga')
    
print 'Changing to', os.getcwd() #This is now different from buildroot   

execfile('test_all.py')

print
print '************************** NOTE *************************************'
print 'If all unit tests passed you should run the suite of validation tests'
print 'Go to the directory anuga_validation/automated_validation_tests'
print 'and run'
print '    python validate_all.py'
print
print 'These tests will take a few hours and will verify that ANUGA'
print 'produces the physical results expected.'
print '*********************************************************************'

    
