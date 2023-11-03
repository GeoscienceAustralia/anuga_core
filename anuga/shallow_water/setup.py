

import os
import sys

from os.path import join
from Cython.Build import cythonize
import Cython.Compiler.Options
Cython.Compiler.Options.annotate = False

os.environ["CC"] = "nvc -O3 -acc=gpu -Minfo=accel -noswitcherror -lm -I$CUDA_HOME/include/ --device-debug --generate-line-info"
os.environ["CXX"] = "nvc++ -O3 -acc=gpu -Minfo=accel -noswitcherror -lm -I$CUDA_HOME/include/ --device-debug --generate-line-info"
os.environ["FC"] = "nvfortran -O3 -acc=gpu -Minfo=accel -noswitcherror -lm -I$CUDA_HOME/include/ --device-debug --generate-line-info"

def configuration(parent_package='',top_path=None):
    
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    
    config = Configuration('shallow_water', parent_package, top_path)

    config.add_data_dir('tests')
    config.add_data_dir(join('tests', 'data'))

    util_dir = join('..', 'utilities')

    config.add_extension('sw_domain_orig_ext',
                         sources=['sw_domain_orig_ext.pyx'],
                         include_dirs=[util_dir])

    config.add_extension('sw_domain_simd_ext',
                         sources=['sw_domain_simd_ext.pyx'],
                         include_dirs=[util_dir])
      
    config.add_extension('sw_domain_openmp_ext',
                         sources=['sw_domain_openmp_ext.pyx'],
                         include_dirs=[util_dir],
                         extra_compile_args=['-fopenmp'],
                         extra_link_args=['-fopenmp'])

    config.add_extension('sw_domain_openacc_ext',
                         sources=['sw_domain_openacc_ext.pyx'],
                         include_dirs=[util_dir],
                         extra_compile_args=None,
                         extra_link_args=None)

    #config.add_extension('sw_domain_cuda_ext',
    #                     sources=['sw_domain_cuda_ext.pyx'],
    #                     include_dirs=[util_dir],
    #                     extra_compile_args=None,
    #                     extra_link_args=None)



    config.ext_modules = cythonize(config.ext_modules, annotate=False)

    return config
    
if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(configuration=configuration)
