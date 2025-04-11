

import os
import sys

from os.path import join
from Cython.Build import cythonize
import Cython.Compiler.Options
Cython.Compiler.Options.annotate = True

def local_fun():
    pass

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('abstract_2d_finite_volumes', parent_package, top_path)

    config.add_data_dir('tests')

    #util_dir = os.path.abspath(join(os.path.dirname(__file__),'..','utilities'))
    #runtime_dir = os.path.abspath(join(os.path.dirname(__file__),'..','runtime_libs'))

    util_dir = join('..','utilities')

    config.add_extension('neighbour_mesh_ext',
                        sources=['neighbour_mesh_ext.pyx'],
                        include_dirs=[util_dir])

    config.add_extension('mesh_factory_ext',
                         sources=['mesh_factory_ext.pyx'],
                         include_dirs=[util_dir])

    config.add_extension('neighbour_table_ext',
                           sources=['neighbour_table_ext.pyx'],
                           extra_compile_args=["-std=c++11"],
                           language='c++',
                           include_dirs=[util_dir])

    config.add_extension('pmesh2domain_ext',
                           sources=['pmesh2domain_ext.pyx'],
                           include_dirs=[util_dir])

    config.add_extension('quantity_ext',
                           sources=['quantity_ext.pyx'],
                           include_dirs=[util_dir])


    config.add_extension('quantity_openmp_ext',
                         sources=['quantity_openmp_ext.pyx'],
                         include_dirs=[util_dir],
                         extra_compile_args=['-fopenmp'],
                         extra_link_args=['-fopenmp'])

    config.add_extension('quantity_openacc_ext',
                           sources=['quantity_openacc_ext.pyx'],
                           include_dirs=[util_dir],
                           extra_compile_args=None,
                           extra_link_args=None)

    config.add_extension('quantity_cuda_ext',
                           sources=['quantity_cuda_ext.pyx'],
                           include_dirs=[util_dir],
                           extra_compile_args=None,
                           extra_link_args=None)


    config.ext_modules = cythonize(config.ext_modules,annotate=True)

    return config

if __name__ == '__main__':
     from numpy.distutils.core import setup
     setup(configuration=configuration)
