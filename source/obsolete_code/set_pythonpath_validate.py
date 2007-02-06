import os

env_var = 'PYTHONPATH'
cwd = os.getcwd()

cwd_source = cwd + os.sep + 'anuga_core'  + os.sep + 'source'

#root = os.sep
#cwd_source = os.path.join(root,'d','cit','1','dgray','validate_inundation','ga','anuga_core','source')

#Not tested
av_dir = os.path.join('..','..','..','anuga_core', 'source')
os.chdir(av_dir)
cwd_source = os.getcwd()
av_dir = os.path.join('..','..','anuga_validation','automated_validation_tests','okushiri_tank_validation')
os.chdir(av_dir)

print 'Old python path  => %s' % os.getenv(env_var)
# This is not a perminent change.  It is only for this execution
os.environ[env_var] = cwd_source
print 'New python path   => %s' % os.getenv(env_var)

import validate_all
#execfile('validate_all.py',{"PYTHONPATH":cwd_source})
