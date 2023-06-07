

import os
import sys

from os.path import join
from Cython.Build import cythonize
import Cython.Compiler.Options
Cython.Compiler.Options.annotate = True


def configuration(parent_package='',top_path=None):
    
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    
    config = Configuration('shallow_water', parent_package, top_path)

    config.add_data_dir('tests')
    config.add_data_dir(join('tests','data'))

    #util_dir = os.path.abspath(join(os.path.dirname(__file__),'..','utilities'))
    util_dir = join('..','utilities')

    config.add_extension('shallow_water_ext',
                         sources=['shallow_water_ext.pyx'],
                         include_dirs=[util_dir])

    config.add_extension('swb2_domain_ext',
                         sources=['swb2_domain_ext.pyx'],
                         include_dirs=[util_dir])

    config.add_extension('swDE_domain_original_ext',
                         sources=['swDE_domain_original_ext.pyx'],
                         include_dirs=[util_dir])

    config.add_extension('swDE_domain_local_timestep_ext',
                         sources=['swDE_domain_local_timestep_ext.pyx'],
                         include_dirs=[util_dir])
    
    # # FIXME SR: come back to getting Mac to run with openmp
    # if sys.platform == 'darwin':
    #     extra_compiler_args = None
    #     extra_link_args = None

    
    config.add_extension('swDE_domain_openmp_ext',
                         sources=['swDE_domain_openmp_ext.pyx'],
                         include_dirs=[util_dir],
                         extra_compile_args=['-fopenmp'],
                         extra_link_args=['-fopenmp'])

    config.add_extension('swDE_domain_openacc_ext',
                         sources=['swDE_domain_openacc_ext.pyx'],
                         include_dirs=[util_dir],
                         extra_compile_args=None,
                         extra_link_args=None)

    config.add_extension('swDE_domain_cuda_ext',
                         sources=['swDE_domain_cuda_ext.pyx'],
                         include_dirs=[util_dir],
                         extra_compile_args=None,
                         extra_link_args=None)



    config.ext_modules = cythonize(config.ext_modules, annotate=True)

    return config
    
if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(configuration=configuration)
