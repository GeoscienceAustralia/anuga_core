#!/bin/bash

# License: 3-clause BSD


set -e

#================================
# Install Ubuntu Packages via apt
#================================

sudo apt install python-is-python3 python3-pip gfortran netcdf-bin libnetcdf-dev libhdf5-serial-dev gdal-bin libgdal-dev libopenmpi-dev openmpi-bin

#=================================
# Install python packages via pip3
#=================================

sudo pip3 install scipy matplotlib nose cython netcdf4 matplotlib dill future gitpython pyproj pymetis triangle Pmw mpi4py ipython

#=================================
# Now install anuga. Should be run 
# from anuga_core directory
#=================================

python setup.py develop --user
python runtests.py
