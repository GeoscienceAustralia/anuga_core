from __future__ import division, print_function

import os
import sys

from os.path import join

def configuration(parent_package='',top_path=None):
    
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    
    config = Configuration('tsunami_source', parent_package, top_path)

    config.add_data_dir('tests')
    config.add_data_dir(join('tests','data'))

    if sys.platform != 'win32':
        config.add_extension('okada_tsunami_fortran',
                            sources=['okada_tsunami_fortran.f'])

    
    return config
    
if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(configuration=configuration)
