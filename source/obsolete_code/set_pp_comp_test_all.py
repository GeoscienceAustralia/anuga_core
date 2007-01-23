import os

env_var = 'PYTHONPATH'
cwd = os.getcwd()

cwd_source = cwd + os.sep  + os.sep + 'source'

print 'Old python path  => %s' % os.getenv(env_var)
# This is not a perminent change.  It is only for this execution
os.environ[env_var] = cwd_source
print 'New python path   => %s' % os.getenv(env_var)

execfile('compile_all.py')
# THIS DOESN'T WORK
execfile('test_all.py',{"PYTHONPATH":cwd_source})
