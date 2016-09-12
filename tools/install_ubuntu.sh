#!/bin/bash
# This script is meant to be called by the "install" step defined in
# .travis.yml. See http://docs.travis-ci.com/ for more details.
# The behavior of the script is controlled by environment variabled defined
# in the .travis.yml in the top level folder of the project.

# License: 3-clause BSD


set -e

PYTHON_VERSION=${PYTHON_VERSION:-"2.7"}
ANUGA_PARALLEL=${ANUGA_PARALLEL:-"false"}

###########################################################
# Check if openmpi has been installed
if [ $(dpkg-query -W -f='${Status}' openmpi-bin 2>/dev/null | grep -c "ok installed") -gt 0 ];
then
  ANUGA_PARALLEL="openmpi"
fi

###########################################################
# Check if mpich2 has been installed
if [ $(dpkg-query -W -f='${Status}' mpich2 2>/dev/null | grep -c "ok installed") -gt 0 ];
then
  ANUGA_PARALLEL="mpich2"
fi


sudo apt-get update -q

##########################################################
# Use standard ubuntu packages in their default version

echo "+===============================================+"
echo "|  Using apt-get to install standard packages   |"
echo "+===============================================+"
    
sudo apt-get install -q -y git gfortran python-dev python-numpy \
                             python-scipy \
                             python-matplotlib netcdf-bin \
                             libnetcdf-dev libhdf5-serial-dev \
                             python-gdal gdal-bin python-pip 


echo "+===============================================+"
echo "|  Using pip to install standard packages       |"
echo "+===============================================+"

sudo pip install nose netCDF4 pyproj
    
##########################################################
# Setup for various versions of MPI
if [[ "$ANUGA_PARALLEL" == "mpich2" ]]; then
    echo "+===============================================+"
    echo "|  Using apt-get to install mpich package       |"
    echo "+===============================================+"
    sudo apt-get install -y mpich2;
fi

if [[ "$ANUGA_PARALLEL" == "openmpi" ]]; then
    echo "+===============================================+"
    echo "|  Using apt-get to install openmpi package       |"
    echo "+===============================================+"
    sudo apt-get install -y libopenmpi-dev openmpi-bin;
fi

# Install pypar if parallel set
if [[ "$ANUGA_PARALLEL" == "mpich2" || "$ANUGA_PARALLEL" == "openmpi" ]]; then
    echo "+===============================================+"
    echo "|  Installing pypar from source                 |"
    echo "+===============================================+"
     git clone https://github.com/daleroberts/pypar.git;
     pushd pypar;
     python setup.py build;
     sudo python setup.py install;
     popd;
fi

#########################################################
# Build and install anuga

echo "+===============================================+"
echo "|  Using setup.py to install anuga              |"
echo "+===============================================+"
python setup.py build
sudo python setup.py install 



