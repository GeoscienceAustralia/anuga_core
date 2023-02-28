#!/bin/bash
# This script is meant to be called by the "install" step defined in
# .travis.yml. See http://docs.travis-ci.com/ for more details.
# The behavior of the script is controlled by environment variabled defined
# in the .travis.yml in the top level folder of the project.

# License: 3-clause BSD


set -e

PYTHON_VERSION=${PYTHON_VERSION:-"3.8"}
ANUGA_PARALLEL=${ANUGA_PARALLEL:-"false"}

echo $PYTHON_VERSION

if [[ "$PYTHON_VERSION" == "3.8" ]];
then
    echo "+===============================================+"
    echo "|  Activate python 3.8 environment              |"
    echo "+===============================================+"
    source ~/virtualenv/python3.8/bin/activate
    python --version
fi

python --version

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
                             libnetcdf-dev libhdf5-serial-dev


echo "+===============================================+"
echo "|  Using apt-get to install gdal                |"
echo "+===============================================+"
#sudo add-apt-repository -y ppa:ubuntugis/ppa
#sudo apt-get update
sudo apt-get install -y gdal-bin libgdal-dev

#sudo apt-get install -y libgdal1-dev

echo "+===============================================+"
echo "|  GDAL version                                 |"
echo "+===============================================+"
ogrinfo --version
#gdal-config --version
#export CPLUS_INCLUDE_PATH=/usr/include/gdal
#export C_INCLUDE_PATH=/usr/include/gdal


echo "+===============================================+"
echo "|  Using pip to install scipy                   |"
echo "+===============================================+"
python -m pip  install -q scipy

echo "+===============================================+"
echo "|  Using pip to install matplotlib              |"
echo "+===============================================+"
python -m pip  install -q matplotlib

echo "+===============================================+"
echo "|  Using pip to install gdal                    |"
echo "+===============================================+"
python -m pip  install setuptools==58
python -m pip  install -q GDAL==3.0.4
#python -m pip  install -q pygdal==3.0.4.11

echo "+===============================================+"
echo "|  Using pip to install pytest                  |"
echo "+===============================================+"
python -m pip  install -q pytest

echo "+===============================================+"
echo "|  Using pip to install dill                    |"
echo "+===============================================+"
python -m pip  install -q dill

echo "+===============================================+"
echo "|  Using pip to install netCDF4                 |"
echo "+===============================================+"
python -m pip  install -q netCDF4

echo "+===============================================+"
echo "|  Using pip to install Cython                  |"
echo "+===============================================+"
python -m pip  install -q Cython

echo "+===============================================+"
echo "|  Using pip to install future                  |"
echo "+===============================================+"
python -m pip  install -q future

echo "+===============================================+"
echo "|  Using pip to install gitpython               |"
echo "+===============================================+"
python -m pip install -q gitpython

echo "+===============================================+"
echo "|  Using pip to install pyproj                  |"
echo "+===============================================+"
python -m pip  install -q pyproj affine

echo "+===============================================+"
echo "|  Using pip to install pymetis                 |"
echo "+===============================================+"
python -m pip install -q pymetis

echo "+===============================================+"
echo "|  Using pip to install meshpy                  |"
echo "+===============================================+"
python -m pip  install -q meshpy

echo "+===============================================+"
echo "|  Using pip to install Pmw                     |"
echo "+===============================================+"
python -m pip  install -q Pmw

echo "+===============================================+"
echo "|  Using pip to install pytz                    |"
echo "+===============================================+"
python -m pip  install -q pytz

echo "+===============================================+"
echo "|  Using pip to install utm                     |"
echo "+===============================================+"
python -m pip  install -q utm

echo "+===============================================+"
echo "|  Using pip to install nbsphinx                |"
echo "+===============================================+"
python -m pip  install -q nbsphinx
    
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
    python -m pip  install -q mpi4py
fi

#########################################################
# Build and install anuga


echo "+===============================================+"
echo "|  Install anuga using pip                      |"
echo "+===============================================+"
python setup.py develop
