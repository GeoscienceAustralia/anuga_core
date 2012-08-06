
import os
anuga_source = os.path.split(os.path.abspath(os.getcwd()))[0]


cmd = 'python setup.py install --install-purelib=%s --install-platlib=%s '%( anuga_source , anuga_source)

print cmd

os.chdir('source')
os.system (cmd)

