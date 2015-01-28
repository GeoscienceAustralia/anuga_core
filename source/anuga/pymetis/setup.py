from distutils.core import setup, Extension;
from numpy.distutils.misc_util import get_numpy_include_dirs;
import os;

metis = Extension('metis',
                  sources = ['metis.c', 'metis_bridge.c'],
                  include_dirs = [os.environ['METIS_DIR'] + os.sep + 'Lib' ] + get_numpy_include_dirs(),
                  libraries = ['metis', 'm'],
                  library_dirs = [os.environ['METIS_DIR']]);

setup(name = 'PyMetis',
      version = '1.0',
      description = 'Python interface to metis',
      ext_modules = [metis]);
