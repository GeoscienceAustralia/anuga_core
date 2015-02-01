import os
import string
import sys
import shutil


# Work out the directory the pyd will go.
_s = string.split( sys.version )[0]
#print "_s",_s
version_info = tuple( map(int, string.split(_s, '.') ) )
del _s
#print "version_info", version_info
#print "sys.executable",sys.executable 
pyd_in_dir = 'build/lib.%s-%s.%s/'  %(sys.platform,version_info[0],version_info[1])

# See if the pyd has already been compiled.
try:
    _ = os.listdir(pyd_in_dir)
except WindowsError:
    # The build hasn't got the files we want
    print "Removing the current build directory"
    try:
        shutil.rmtree('build')
    except:
        pass

    

if sys.platform == 'win32':  #Windows
    win32_extra = '-cmingw32'

command = '%s setup.py build %s' %(sys.executable,win32_extra) 
os.system(command)
    

# 

#print "file_in_dir", file_in_dir 
for file in os.listdir(pyd_in_dir):
    if os.path.splitext(file)[1] == ".pyd":
        shutil.copy(os.path.join(pyd_in_dir,file), '.')


