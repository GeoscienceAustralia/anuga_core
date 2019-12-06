from __future__ import division, print_function

import os
import sys

from os.path import join
from Cython.Build import cythonize
import Cython.Compiler.Options
Cython.Compiler.Options.annotate = True

def configuration(parent_package='',top_path=None):
    
    from numpy.distutils.misc_util import Configuration
    
    config = Configuration('pymetis', parent_package, top_path)

    config.add_data_dir('tests')

        
    METIS_DIR = 'metis-4.0'
    
    metis_src = [join(METIS_DIR,'*.c')]
    metis_headers = [join(METIS_DIR,'*.h')]
    include_dir = METIS_DIR 
    
    #print(metis_headers)
    #print(metis_src)
    #print(include_dir)
    
    config.add_include_dirs([include_dir])
    
    config.add_library('metis', sources=metis_src)
    
    
    src_files = ['metis_ext.pyx']
    if sys.platform == 'win32':
        src_files = src_files + ['random.c']
        
    

    config.add_extension('metis_ext',
                         sources=src_files,
                         include_dirs = [include_dir],
                         depends=(metis_src + metis_headers),
                         extra_compile_args=['-I'+include_dir],
                         libraries = ['metis', 'm'])

    config.ext_modules = cythonize(config.ext_modules,annotate=True)

    return config
    
if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(configuration=configuration)
