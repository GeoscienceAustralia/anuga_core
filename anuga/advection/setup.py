from __future__ import division, print_function

import os
import sys

from os.path import join
from Cython.Build import cythonize
import Cython.Compiler.Options
Cython.Compiler.Options.annotate = True

def configuration(parent_package='',top_path=None):
    
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    
    config = Configuration('advection', parent_package, top_path)

    config.add_data_dir('tests')

    #util_dir = os.path.abspath(join(os.path.dirname(__file__),'..','utilities'))
    util_dir = join('..','utilities')
            
    config.add_extension('advection_ext',
                         sources=['advection_ext.pyx'],
                         include_dirs=[util_dir])

    config.ext_modules = cythonize(config.ext_modules,annotate=True)

    return config
    
if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(configuration=configuration)
