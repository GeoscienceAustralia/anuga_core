#!/bin/bash
# This script is meant to be called by the "install" step defined in
# .travis.yml. See http://docs.travis-ci.com/ for more details.
# The behavior of the script is controlled by environment variabled defined
# in the .travis.yml in the top level folder of the project.

# License: 3-clause BSD


set -e

PYTHON_VERSION=${PYTHON_VERSION:-"2.7"}
ANUGA_PARALLEL=${ANUGA_PARALLEL:-"false"}

if [[ "$PYTHON_VERSION" == "3.8" ]]; 
then 
    source ~/virtualenv/python3.8/bin/activate
fi

if [[ "$ANUGA_PARALLEL" == "false" ]];
then 
    PYPAR_AVAILABLE="false"
else
    PYPAR_AVAILABLE=${PYPAR_AVAILABLE:-"pypar"}
fi

###########################################################
# Check if mpich2 has been installed
if [ $(dpkg-query -W -f='${Status}\n' mpich2 2>/dev/null | grep -c "ok installed") -gt 0 ];
then
  ANUGA_PARALLEL="mpich2"
fi

###########################################################
# Check if mpich has been installed
if [ $(dpkg-query -W -f='${Status}' mpich 2>/dev/null | grep -c "ok installed") -gt 0 ];
then
  ANUGA_PARALLEL="mpich"
fi


###########################################################
# Check if openmpi has been installed
if [ $(dpkg-query -W -f='${Status}\n' openmpi-bin 2>/dev/null | grep -c "ok installed") -gt 0 ];
then
  ANUGA_PARALLEL="openmpi"
fi


sudo apt-get update -q

##########################################################
# Use standard ubuntu packages in their default version

echo "+===============================================+"
echo "|  Using apt-get to install standard packages   |"
echo "+===============================================+"

sudo apt-get install -q -y git gfortran netcdf-bin \
                             libnetcdf-dev libhdf5-serial-dev \
                             gdal-bin

echo "+===============================================+"
echo "|  Using pip to install scipy                   |"
echo "+===============================================+"
sudo pip install -q scipy


echo "+===============================================+"
echo "|  Using pip to install matplotlib              |"
echo "+===============================================+"
sudo pip install -q matplotlib


echo "+===============================================+"
echo "|  Using pip to install gdal                    |"
echo "+===============================================+"
sudo pip install -q gdal


echo "+===============================================+"
echo "|  Using pip to install nose                    |"
echo "+===============================================+"
sudo pip install -q nose

echo "+===============================================+"
echo "|  Using pip to install dill                    |"
echo "+===============================================+"
sudo pip install -q dill

echo "+===============================================+"
echo "|  Using pip to install netCDF4                 |"
echo "+===============================================+"
sudo pip install -q netCDF4

echo "+===============================================+"
echo "|  Using pip to install Cython                  |"
echo "+===============================================+"
sudo pip install -q Cython

echo "+===============================================+"
echo "|  Using pip to install future                  |"
echo "+===============================================+"
sudo pip install -q future

echo "+===============================================+"
echo "|  Using pip to install pyproj                  |"
echo "+===============================================+"
sudo pip install -q pyproj


##########################################################
# Setup for various versions of MPI
if [[ "$ANUGA_PARALLEL" == "mpich" ]]; then
    echo "+===============================================+"
    echo "|  Using apt-get to install mpich package       |"
    echo "+===============================================+"
    sudo apt-get install -q -y mpich;
fi

if [[ "$ANUGA_PARALLEL" == "mpich2" ]]; then
    echo "+===============================================+"
    echo "|  Using apt-get to install mpich2 package      |"
    echo "+===============================================+"
    sudo apt-get install -q -y mpich2;
fi

if [[ "$ANUGA_PARALLEL" == "openmpi" ]]; then
    echo "+===============================================+"
    echo "|  Using apt-get to install openmpi package     |"
    echo "+===============================================+"
    sudo apt-get install -q -y libopenmpi-dev openmpi-bin;
fi


# Install pypar if parallel set
if [[ "$PYPAR_AVAILABLE" == "pypar" ]]; then
    echo "+===============================================+"
    echo "|  Installing pypar from source                 |"
    echo "+===============================================+"
    if [ ! -d "pypar" ]; then
	git clone https://github.com/daleroberts/pypar.git;
    fi
    pushd pypar;
    git pull
    python setup.py  build;
    sudo python setup.py  install;
    popd;
fi

if [[ "$PYPAR_AVAILABLE" == "mpi4py" ]]; then 
    echo "+===============================================+"
    echo "|  Using pip to install mpi4py                  |"
    echo "+===============================================+"
    sudo pip install -q mpi4py
fi  

#########################################################
# Build and install anuga

echo "+===============================================+"
echo "|  Build anuga                                  |"
echo "+===============================================+"
python build_all.py

echo "+===============================================+"
echo "|  Install anuga using setup.py                 |"
echo "+===============================================+"
sudo python setup.py -q install
