#!/bin/bash
# This script is meant to be called by the "install" step defined in
# .travis.yml. See http://docs.travis-ci.com/ for more details.
# The behavior of the script is controlled by environment variabled defined
# in the .travis.yml in the top level folder of the project.

# License: 3-clause BSD


set -e


PYTHON_VERSION=${PYTHON_VERSION:-"3.10"}
ANUGA_BITS=${ANUGA_BITS:-"64"}


# Deactivate the travis-provided virtual environment and setup a
# conda-based environment instead
deactivate || echo "deactivate failed"

# Use the miniforge installer
wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"


chmod +x Miniforge3.sh && ./Miniforge3.sh -b

export PATH=/home/travis/miniforge3/bin:$PATH 

conda init bash

source ~/.bashrc

conda update --yes conda

# Configure the conda environment
conda create -n anuga_env --yes python=$PYTHON_VERSION compilers numpy scipy cython netcdf4 \
     nose matplotlib gdal dill gitpython mpi4py utm Pmw pymetis meshpy pytest pyproj affine \
     meson-python meson ninja xarray future pkg-config
     
#conda create -n anuga_env --yes python=$PYTHON_VERSION pip numpy scipy meshpy cython netcdf4 pytest matplotlib gdal dill gitpython Pmw pymetis utm pyproj affine xarray mpi4py
conda activate anuga_env

# Useful for debugging any issues with conda
conda info -a

# Build and install anuga
pip install .
