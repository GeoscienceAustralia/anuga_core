#!/bin/bash

# License: 3-clause BSD


set -e

echo "#==========================="
echo "# Install packages via apt"
echo "#==========================="

sudo apt install -y -q build-essential python-dev-is-python3 gfortran netcdf-bin \
    libnetcdf-dev libhdf5-serial-dev gdal-bin libgdal-dev libopenmpi-dev openmpi-bin \
    python3-venv

# sudo apt-get install -y build-essential cmake python3-dev python-dev-is-python3 \
#      netcdf-bin libnetcdf-dev libhdf5-serial-dev libpq-dev libgeos-dev libexpat-dev \
#     libxerces-c-dev libwebp-dev libpng-dev libzstd-dev libssl-dev libopenjp2-7-dev \
#     libspatialite-dev libmuparser-dev autoconf automake sqlite3 bash-completion swig \
#     libopenmpi-dev openmpi-bin python3-venv

#build gdal from source
# cd /tmp
# git clone https://github.com/OSGeo/gdal.git
# cd gdal
# mkdir build
# cd build

# cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON
# cmake --build .
# sudo cmake --install .
# sudo ldconfig
# gdalinfo --version


echo "#======================================="
echo "# Create a virtual environment and then"
echo "# install python packages via pip"
echo "#======================================="

cd "$(dirname "${BASH_SOURCE[0]}")"/..

ANUGA_CORE_PATH=`pwd`
echo "ANUGA_CORE_PATH: $ANUGA_CORE_PATH"

# ensure meson picks the pip installed numpy and not system numpy
export PKG_CONFIG_PATH="${ANUGA_CORE_PATH}anuga_env/lib/python3.12/site-packages/numpy/_core/lib/pkgconfig/:$PKG_CONFIG_PATH"
echo "$PKG_CONFIG_PATH"

python3 -m venv anuga_env
source anuga_env/bin/activate
pip install wheel numpy==2.2 scipy gdal matplotlib pytest cython netcdf4 \
     matplotlib dill future gitpython pyproj pymetis pybind11 meshpy Pmw ipython \
     utm affine mpi4py xarray meson meson-python ninja

echo "#==============================================="
echo "# Installing anuga from the anuga_core directory"
echo "#==============================================="

# get numpy include path
NUMPY_INCLUDE_PATH=$(python3 -c "import numpy; print(numpy.get_include())")
pip install .

echo "#==========================="
echo "# Run unittests"
echo "# "
echo "# At present 22/02/2025 the tests are failing"
echo "# due to a incompatiblity between numpy 2.2 and gdal 3.4.1"
echo "# Suggest installing anuga in a miniforge3 environment"
echo "# using the script install_miniforge.sh"
echo "#==========================="
pytest -q --pyargs anuga

echo "#================================================"
echo "# To use anuga you must activate the"
echo "# python environment anuga_env that has been"
echo "# created in your anuga_core directory by"
echo "# changing to anuga_core directory then running:"
echo "# "
echo "# source anuga_env/bin/activate"
echo "# "
echo "# After activating environment use:"
echo "# "
echo "# pip install mpi4py "
echo "# "
echo "# to enable parallel execution"
echo "#================================================"




