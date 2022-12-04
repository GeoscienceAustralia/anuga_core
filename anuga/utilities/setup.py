

import os
import sys

from os.path import join
from Cython.Build import cythonize
import Cython.Compiler.Options
Cython.Compiler.Options.annotate = True

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    config = Configuration('utilities', parent_package, top_path)

    config.add_data_dir('tests')
    config.add_data_dir(join('tests','data'))

    config.add_extension('sparse_ext',
                         sources='sparse_ext.pyx')

    config.add_extension('sparse_matrix_ext',
                         sources=['sparse_matrix_ext.pyx'])


    config.add_extension('util_ext',
                         sources='util_ext_c.pyx')

    if sys.platform == 'darwin':
        extra_args = None
    else:
        extra_args = ['-fopenmp']

    config.add_extension('cg_ext',
                         sources='cg_ext.pyx',
                         extra_compile_args=extra_args,
                         extra_link_args=extra_args)

    config.add_extension('quad_tree_ext',
                         sources=['quad_tree_ext.pyx'])
    
    config.ext_modules = cythonize(config.ext_modules,annotate=True)

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(configuration=configuration)
