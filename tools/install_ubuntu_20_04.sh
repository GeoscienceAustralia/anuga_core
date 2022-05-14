#!/bin/bash

# License: 3-clause BSD


set -e

echo "#==========================="
echo "# Install packages via apt"
echo "#==========================="

sudo apt install -y -q python3-pip gfortran netcdf-bin libnetcdf-dev libhdf5-serial-dev gdal-bin libgdal-dev libopenmpi-dev openmpi-bin

echo "#==========================="
echo "# Create a virtual environment and then"
echo "# install python packages via pip"
echo "#==========================="

python3 -m venv anuga_env
source anuga_env/bin/activate
pip install scipy gdal matplotlib pytest nose cython netcdf4 matplotlib dill future gitpython pyproj pymetis pybind11 meshpy Pmw mpi4py ipython pytz utm

echo "#==========================="
echo "# You should now install anuga from the anuga_core directory"
echo "# via"
echo "cd anuga_core"
echo "pip install -e ."
echo "# and then test via"
echo "# pytest"
echo "#==========================="

# cd "$(dirname "${BASH_SOURCE[0]}")"/..
# python setup.py develop --user


