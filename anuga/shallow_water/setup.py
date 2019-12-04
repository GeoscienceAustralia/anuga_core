from __future__ import division, print_function

import os
import sys

from os.path import join

def configuration(parent_package='',top_path=None):
    
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    
    config = Configuration('shallow_water', parent_package, top_path)

    config.add_data_dir('tests')
    config.add_data_dir(join('tests','data'))

    #util_dir = os.path.abspath(join(os.path.dirname(__file__),'..','utilities'))
    util_dir = join('..','utilities')
    
    config.add_extension('shallow_water_ext',
                         sources=['shallow_water_ext.c'],
                         include_dirs=[util_dir])
    
    config.add_extension('swb2_domain_ext',
                         sources=['swb2_domain_ext.c'],
                         include_dirs=[util_dir])

    config.add_extension('swDE1_domain_ext',
                         sources=['swDE1_domain_ext.c'],
                         include_dirs=[util_dir])


    return config
    
if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(configuration=configuration)
