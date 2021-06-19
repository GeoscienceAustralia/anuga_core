#!/bin/bash

# License: 3-clause BSD


set -e

#================================
# Install Ubuntu Packages via apt
#================================

echo "#==========================="
echo "# Install packages via apt"
echo "#==========================="

sudo apt install -y -q python-is-python3 python3-pip gfortran netcdf-bin libnetcdf-dev libhdf5-serial-dev gdal-bin libgdal-dev libopenmpi-dev openmpi-bin

#=================================
# Install python packages via pip3
#=================================

echo "#==========================="
echo "# Install python packages via pip"
echo "#==========================="

sudo pip3 install scipy matplotlib nose cython netcdf4 matplotlib dill future gitpython pyproj pymetis triangle Pmw mpi4py ipython

#=================================
# Now install anuga. Should be run 
# from anuga_core directory
#=================================

echo "#==========================="
echo "# Install anuga"
echo "#==========================="

python setup.py develop --user

echo "#==========================="
echo "# Run Unit tests"
echo "#==========================="

python runtests.py


