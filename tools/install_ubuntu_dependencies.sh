#!/bin/bash
# This script is meant to be called by the "install" step defined in
# .travis.yml. See http://docs.travis-ci.com/ for more details.
# The behavior of the script is controlled by environment variabled defined
# in the .travis.yml in the top level folder of the project.
#
# But can also be used to install anuga dependencies

# License: 3-clause BSD


set -e

PYTHON_VERSION=${PYTHON_VERSION:-"2.7"}
ANUGA_PARALLEL=${ANUGA_PARALLEL:-"openmpi"}

apt-get update -q

##########################################################
# Use standard ubuntu packages in their default version
    
apt-get install -q -y git gfortran python-dev python-numpy \
                             python-scipy \
                             python-matplotlib netcdf-bin \
                             libnetcdf-dev libhdf5-serial-dev \
                             python-gdal gdal-bin python-pip 


pip install nose netCDF4 pyproj
    
##########################################################
# Setup for various versions of MPI
if [[ "$ANUGA_PARALLEL" == "mpich2" ]]; then
    apt-get install -y mpich2;
fi

if [[ "$ANUGA_PARALLEL" == "openmpi" ]]; then
    apt-get install -y libopenmpi-dev openmpi-bin;
fi

# Install pypar if parallel set
if [[ "$ANUGA_PARALLEL" == "mpich2" || "$ANUGA_PARALLEL" == "openmpi" ]]; then
     git clone https://github.com/daleroberts/pypar.git;
     pushd pypar;
     python setup.py build;
     sudo python setup.py install;
     popd;
fi




