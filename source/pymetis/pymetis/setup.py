from distutils.core import setup, Extension;
import os;

metis = Extension('metis',
                  sources = ['metis.c', 'metis_bridge.c'],
                  include_dirs = ['..' + os.sep + os.environ['METIS_DIR'] + os.sep + 'Lib' ,  '/usr/local/lib/python2.5/site-packages/numpy/core/include' ],
                  libraries = ['metis', 'm'],
                  library_dirs = ['..' + os.sep + os.environ['METIS_DIR']]);

setup(name = 'PyMetis',
      version = '1.0',
      description = 'Python interface to metis',
      ext_modules = [metis]);
