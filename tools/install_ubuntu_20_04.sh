#!/bin/bash

# License: 3-clause BSD


set -e

echo "#==========================="
echo "# Install packages via apt"
echo "#==========================="

sudo apt install -y -q python-is-python3 python3-pip gfortran netcdf-bin libnetcdf-dev libhdf5-serial-dev gdal-bin libgdal-dev libopenmpi-dev openmpi-bin

echo "#==========================="
echo "# Install python packages via pip3"
echo "#==========================="

sudo pip3 install scipy matplotlib nose cython netcdf4 matplotlib dill future gitpython pyproj pymetis triangle Pmw mpi4py ipython

echo "#==========================="
echo "# You should now install anuga from the anuga_core directory"
echo "# via"
echo "cd anuga_core"
echo "python setup.py install --user"
echo "#==========================="

# cd "$(dirname "${BASH_SOURCE[0]}")"/..
# python setup.py develop --user


