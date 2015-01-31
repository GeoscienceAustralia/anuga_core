from __future__ import division, print_function

import os
import sys

from os.path import join

def configuration(parent_package='',top_path=None):
    
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    
    config = Configuration('fit_interpolate', parent_package, top_path)

    config.add_data_dir('tests')


    util_dir = os.path.abspath(join(os.path.dirname(__file__),'..','utilities'))

    
    util_srcs = [join(util_dir,'quad_tree.c'),
                 join(util_dir,'sparse_dok.c'),
                 join(util_dir,'sparse_csr.c')]
    

    config.add_extension('fitsmooth',
                         sources=['fitsmooth.c']+util_srcs,
                         include_dirs=[util_dir],
                         extra_compile_args=['-fopenmp'],
                         extra_link_args=['-fopenmp'])


    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(configuration=configuration)
