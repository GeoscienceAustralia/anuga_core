import os
import time
import sys
import subprocess

buildroot = os.getcwd()
t0 = time.time()

#--------------------------------------------------
# Compiling anuga code 
#--------------------------------------------------
os.chdir('source')

print 'Changing to', os.getcwd()        

cmd = 'python setup.py build_ext --inplace -f '
print cmd
err = os.system(cmd)
if err != 0:
    msg = 'Could not compile anuga '
    msg += 'on platform %s, %s\n' % (sys.platform, os.name)
    raise Exception, msg
else:
    print 50*"-"
    print
    msg = 'Compiled anuga successfully.'
    print msg

# 
# 
# 
# 
# print 'Changing to', os.getcwd()        
# 
# #entries = listdir('.')
# 
# t0 = time.time()
# 
# # Attempt to compile all ANUGA extensions
# execfile('compile_all.py')
# 
# 
# os.chdir(buildroot)
# 
# 
# #--------------------------------------------------
# # Compiling anuga_1d code
# # excluded for the time being
# #--------------------------------------------------
# #os.chdir('anuga_1d')
# 
# 
# #t0 = time.time()
# 
# # Attempt to compile all ANUGA_1D extensions
# #execfile('compile_all.py')
# 
# 
# os.chdir(buildroot)
# print
# print 'Changing to', os.getcwd()        
# 
# #--------------------------------------------------
# # Compiling partition code
# #--------------------------------------------------
# 
# try:
#     print
#     print '-----------------------------------------------'
#     print 'Attempting to compile Metis'
#     print '-----------------------------------------------'
# 
#     import pypar
# 
#     # Attempt to compile Metis for use with anuga_parallel
#     os.chdir('source')
#     os.chdir('anuga')
#     #os.chdir('parallel')
#     os.chdir('pymetis')
# 
#     print 'Changing to', os.getcwd()        
# 
#     cmd = 'python setup.py build_ext --inplace -f '
#     print cmd
#     err = os.system(cmd)
#     if err != 0:
#         msg = 'Could not compile pymetis '
#         msg += 'on platform %s, %s\n' % (sys.platform, os.name)
#         msg += 'You need to compile metis manually '
#         msg += 'if you want to run ANUGA in parallel.'
#         raise Exception, msg
#     else:
#         msg = 'Compiled pymetis successfully.'
#         print msg
# 
# 
#         
# except:
#     print 
#     print 'pymetis could not compile as pypar not installed'
# 
# 
# os.chdir(buildroot)
# print 'Changing to', os.getcwd() 
# 
# #--------------------------------------------------
# # Compiling pypar_extras
# #--------------------------------------------------
# 
# try:
#     print 
#     print '-----------------------------------------------'
#     print 'Attempting to compile mpiextras'
#     print '-----------------------------------------------'
# 
#     import pypar
#     
#     os.chdir('source')
#     os.chdir('anuga')
#     os.chdir('parallel')
#     #os.chdir('pypar_extras')
# 
#     cmd = 'python setup.py build_ext --inplace -f '
#     print cmd
#     err = os.system(cmd)
#     if err != 0:
#         msg = 'Could not compile mpiextras '
#         msg += 'on platform %s, %s\n' % (sys.platform, os.name)
#         msg += 'You need to compile mpiextras manually '
#         msg += 'if you want to run ANUGA in parallel.'
#         raise Exception, msg
#     else:
#         msg = 'Compiled mpiextras successfully.'
#         print msg
# except:
#     print 'anuga.parallel code not compiled as pypar not installed'
# 
# 

print        
print 'That took %.3fs' %(time.time() - t0)




