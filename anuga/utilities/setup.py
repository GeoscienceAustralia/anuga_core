from __future__ import division, print_function

import os
import sys

from os.path import join

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    config = Configuration('utilities', parent_package, top_path)

    config.add_data_dir('tests')
    config.add_data_dir(join('tests','data'))

    config.add_extension('sparse_ext',
                         sources='sparse_ext.c')

    config.add_extension('sparse_matrix_ext',
                         sources=['sparse_matrix_ext.c', 'sparse_dok.c'])


    config.add_extension('util_ext',
                         sources='util_ext.c')

    config.add_extension('cg_ext',
                         sources='cg_ext.c',
                         extra_compile_args=['-fopenmp'],
                         extra_link_args=['-fopenmp'])

    config.add_extension('quad_tree_ext',
                         sources=['quad_tree_ext.c', 'quad_tree.c'])
    

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(configuration=configuration)
