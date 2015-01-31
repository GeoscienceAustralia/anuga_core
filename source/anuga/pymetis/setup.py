from __future__ import division, print_function

import os
import sys

from os.path import join

def configuration(parent_package='',top_path=None):
    
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    
    config = Configuration('pymetis', parent_package, top_path)

    config.add_data_dir('tests')

    if parent_package is '':
        anuga_dir = '..'
    else:
        anuga_dir = '.'
        
    METIS_DIR = 'metis-4.0'
    
    metis_src = [join(METIS_DIR,'*.c')]
    metis_headers = [join(METIS_DIR,'*.h')]
    include_dir = METIS_DIR 
    
    #print(metis_headers)
    #print(metis_src)
    #print(include_dir)
    
    config.add_include_dirs([include_dir])
    
    config.add_library('metis', sources=metis_src)
    
    
    src_files = ['metis_ext.c', 'metis_bridge.c']
    if sys.platform == 'win32':
        src_files = src_files + ['random.c']
        
    

    config.add_extension('metis_ext',
                         sources=src_files,
                         include_dirs = [include_dir],
                         depends=(metis_src + metis_headers),
                         extra_compile_args=['-I'+include_dir],
                         libraries = ['metis', 'm'])

    return config
    
if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(configuration=configuration)
