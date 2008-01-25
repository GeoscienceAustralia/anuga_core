"""f2py script -

"""

import os

separation_line = '---------------------------------------'
print separation_line

command = 'python setup.py build_src build_ext --inplace'
print 'Trying to run %s in directory %s' % (command, os.getcwd())
os.system(command)


