#!/bin/bash

# License: 3-clause BSD


set -e

echo "#==========================="
echo "# Install packages via apt"
echo "#==========================="

sudo apt install -y -q build-essential python-dev-is-python3 gfortran netcdf-bin libnetcdf-dev libhdf5-serial-dev gdal-bin libgdal-dev libopenmpi-dev openmpi-bin

echo "#==========================="
echo "# Create a virtual environment and then"
echo "# install python packages via pip"
echo "#==========================="

cd "$(dirname "${BASH_SOURCE[0]}")"/..
python3 -m venv anuga_env
source anuga_env/bin/activate
pip install wheel scipy gdal==3.4.1 matplotlib pytest cython netcdf4 matplotlib dill future gitpython pyproj pymetis pybind11 meshpy Pmw ipython utm mpi4py

echo "#==========================="
echo "# Installing anuga from the anuga_core directory"
echo "# and then run unittests"
echo "#==========================="

pip install -e .
pytest --pyargs anuga

echo "#================================================"
echo "# To use anuga you must activate the"
echo "# python environment anuga_env that has been"
echo "# created in your anuga_core directory by"
echo "# changing to anuga_core directory then running:"
echo "# source anuga_env/bin/activate"
echo "#================================================"


