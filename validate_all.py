import os

os.chdir('source')
os.chdir('anuga_validation_tests')
print
print 20*'=' + ' anuga automated validation tests ' + 20*'='
print 'Changing to', os.getcwd() # This is now different from buildroot   
execfile('validate_all.py')







    
