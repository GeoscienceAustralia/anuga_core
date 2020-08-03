from __future__ import division, print_function

from future import standard_library
standard_library.install_aliases()
import os
import sys
import subprocess
import shlex
import string

from os.path import join

#=================================================
# Code taken from pysph
#=================================================    
    
def getoutput_mpicc():
    """Returns the output of the command used to compile using
    mpicc."""
    # LAM/OPENMPI/MPICH2
    output = subprocess.getoutput('mpicc -show') + ' -fPIC'

    if output:
        return output

    # MPICH
    # works with MPICH version 1.2.1 (on Debian)
    output = subprocess.getoutput('mpicc -compile_info -link_info')
    if output:
        return output

def parse_command(output):
    # Now get the include, library dirs and the libs to link.
    flags = shlex.split(output)
    #flags = uniq_arr(flags) # Remove repeated values.
    inc_dirs = []
    lib_dirs = []
    libs = []
    def_macros = []
    undef_macros = []
    for f in flags:
        if f[:2] == '-I':
            inc_dirs.append(f[2:])
        elif f[:2] == '-L':
            lib_dirs.append(f[2:])
        elif f[:2] == '-l' and f[-1] != "'": # Patched by Michael McKerns July 2009
            libs.append(f[2:])
        elif f[:2] == '-U':
            undef_macros.append(f[2:])
        elif f[:2] == '-D':
            tmp = string.split(f[2:], '=')
            if len(tmp) == 1:
                def_macros.append((tmp[0], None))
            else:
                def_macros.append(tuple(tmp))
    return {'inc_dirs': inc_dirs, 'lib_dirs': lib_dirs, 'libs':libs,
            'def_macros': def_macros, 'undef_macros': undef_macros}




def configuration(parent_package='',top_path=None):
    
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    
  
    
    #print(mpi_flags)
    
    config = Configuration('parallel', parent_package, top_path)

    try:
        # Use this import to check if we are in a parallel environment
        from anuga.utilities import parallel_abstraction as pypar

        if pypar.pypar_available:
            #We are parallel!
            mpi_flags = parse_command(getoutput_mpicc())

            config.add_data_dir('tests')
            config.add_data_dir('data')

    except:
        #No parallel support, so just copy over the py files
        pass

    
    return config
    
if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(configuration=configuration)



