import os


os.chdir('validation_tests')
print
print 20*'=' + ' anuga automated validation tests ' + 20*'='
print 'Changing to', os.getcwd() # This is now different from buildroot   
execfile('run_auto_validation_tests.py')







    
