#!/usr/bin/env python

# Examples of setup
#
# python setup.py install
# python setup.py install --prefix=/opt/python-2.3
# python setup.py install --home=~
#
# Some times you'll have to add a file called pypar.pth
# containing the word pypar to site-packages
#
# See http://docs.python.org/dist/pure-pkg.html for more on distutils

# FIXME: Now mpiext.c and pypar.py are assumed to be in this directory.
# Maybe, we should put them in the default package directory, pypar.
# The repository structure would then be
# 
# pypar
#     demos
#     documentation
#     source
#          pypar

from distutils.core import setup, Extension

# FIXME (Ole): This works, but I don't know how to use it
# Generate Python EGG if possible.
#try:
#   from setuptools import setup, Extension
#except ImportError:
#   pass

import distutils.sysconfig
import distutils.debug
import os, sys
import popen2
import string
import tempfile
import numpy
from __metadata__ import __version__, __date__, __author__


def setup_compiler():
    distutils.sysconfig.get_config_vars()
    config_vars = distutils.sysconfig._config_vars
    
    if sys.platform == 'sunos5':
        config_vars['LDSHARED'] = "gcc -G"
        config_vars['CCSHARED'] = ""
        

def uniq_arr(arr):
    """Remove repeated values from an array and return new array."""
    ret = []
    for i in arr:
        if i not in ret:
            ret.append(i)
    return ret


def _run_command(cmd):
    out_file, in_file, err_file = popen2.popen3(cmd)
    output = out_file.read() + err_file.read()
    out_file.close()
    in_file.close()
    err_file.close()
    # need this hack to get the exit status
    out_file = os.popen(cmd)
    if out_file.close():
        # close returns exit status of command.
        return ""
    else:
        # no errors, out_file.close() returns None.
        return output


def _get_mpi_cmd():
    """Returns the output of the command used to compile using
    mpicc."""
    # LAM
    output = _run_command("mpicc -showme")
    if output:
        return output

    # MPICH
    # works with MPICH version 1.2.1 (on Debian)
    output = _run_command("mpicc -compile_info -link_info")
    if output:
        return output

    # old version of MPICH needs this hack.
    tmp_base = tempfile.mktemp()
    tmp_c = tmp_base + ".c"
    tmp_o = tmp_base + ".o"
    tmp_file = open(tmp_c, "w")
    tmp_file.write('#include "mpi.h"\nint main(){return 0;}\n')
    tmp_file.close()
    output = _run_command("mpicc -show;"\
                          "mpicc -echo -c %s -o %s"%(tmp_c, tmp_o))
    os.remove(tmp_c)
    if os.path.exists(tmp_o):
        os.remove(tmp_o)
    if output:
        return output
    else:
        return ""


def get_mpi_flags():
    output = _get_mpi_cmd()
    print output
    if not output:
        if sys.platform=='win32': #From Simon Frost
            output = "gcc -L$MPICH_DIR\SDK.gcc\lib -lmpich -I$MPICH_DIR\SDK.gcc\include"
        else:
            output = "cc -L/usr/opt/mpi -lmpi -lelan"


    # now get the include, library dirs and the libs to link with.
    flags = string.split(output)
    flags = uniq_arr(flags) # remove repeated values.
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
        elif f[:2] == '-l':
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


if __name__ == "__main__":
    setup_compiler()
    
    mpi_flags = get_mpi_flags()
    mpi_flags['inc_dirs'].append(numpy.get_include())


    # setting some extra compile flags for AMD64, utilizing
    # distutils.sysconfig to check which compiler to use
    if os.name == 'posix' and os.uname()[4] == 'x86_64':
        #Extra flags for 64 bit architectures
        if 'pgcc' in distutils.sysconfig.get_config_var('CC'):
            extra_compile_args = [' -fPIC -tp amd64'] #Valid for pgcc
        elif 'gcc' in distutils.sysconfig.get_config_var('CC'):
            extra_compile_args = [' -fPIC -m64'] #Valid for gcc
        elif 'icc' in distutils.sysconfig.get_config_var('CC'):
            extra_compile_args = [' -fPIC'] #Valid for icc
        else:
            extra_compile_args = None
    else:
        extra_compile_args = None



    setup(name='Pypar',
          version=__version__,
          description='Pypar - Parallel Python',
          long_description='Pypar - Parallel Python, no-frills MPI interface',
          author=__author__,
          author_email='ole.moller.nielsen@gmail.com',
          url='http://sourceforge.net/projects/pypar',
          package_dir = {'pypar': ''}, # Use files in this dirctory 
          packages  = ['pypar'],
          ext_modules = [Extension('pypar.mpiext',
                                   ['mpiext.c'], 
                                   include_dirs=mpi_flags['inc_dirs'],
                                   library_dirs=mpi_flags['lib_dirs'],
                                   libraries=mpi_flags['libs'],
                                   define_macros=mpi_flags['def_macros'],
                                   undef_macros=mpi_flags['undef_macros'],
                                   extra_compile_args=extra_compile_args)]
         )
