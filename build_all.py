from __future__ import print_function
import os
import time
import sys
import subprocess

buildroot = os.getcwd()
t0 = time.time()

#--------------------------------------------------
# Compiling anuga code 
#--------------------------------------------------
#os.chdir('source')

#print 'Changing to', os.getcwd()        

if sys.platform == 'win32':
    cmd = 'python setup.py build --compiler=mingw32'
else:
    cmd = 'python setup.py build'
print(cmd)

log_filename_1 = 'build.log'
log_filename_2 = 'build_err.log'
print("Building, see build.log...")
with open(log_filename_1, 'w') as log1:
    with open(log_filename_2, 'w') as log2:
        p = subprocess.Popen(cmd, stdout=log1, stderr=log2, shell=True)

# Wait for it to finish, and print something to indicate the
# process is alive, but only if the log file has grown (to
# allow continuous integration environments kill a hanging
# process accurately if it produces no output)
last_blip = time.time()
last_log_size = os.stat(log_filename_1).st_size
while p.poll() is None:
    time.sleep(0.5)
    if time.time() - last_blip > 10:
        log_size = os.stat(log_filename_1).st_size
        if log_size > last_log_size:
            print("    ... build in progress")
            last_blip = time.time()
            last_log_size = log_size

ret = p.wait()

if ret == 0:
    print("Build OK")
else:
    with open(log_filename_2, 'r') as f:
        print(f.read())
    print("Build failed!")
    sys.exit(1)




print()        
print('That took %.3fs' %(time.time() - t0))




