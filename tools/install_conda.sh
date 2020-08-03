#!/bin/bash
# This script is meant to be called by the "install" step defined in
# .travis.yml. See http://docs.travis-ci.com/ for more details.
# The behavior of the script is controlled by environment variabled defined
# in the .travis.yml in the top level folder of the project.

# License: 3-clause BSD


set -e


PYTHON_VERSION=${PYTHON_VERSION:-"2.7"}
ANUGA_BITS=${ANUGA_BITS:-"64"}



sudo apt-get update -q
sudo apt-get install gfortran git
sudo apt-get install -y libopenmpi-dev openmpi-bin;

##########################################################

# Deactivate the travis-provided virtual environment and setup a
# conda-based environment instead
deactivate || echo "deactivate failed"

# Use the miniconda installer for faster download
# install of conda itself
if [[ "$ANUGA_BITS" == "64" ]]; then
    wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh ;
fi
if [[ "$ANUGA_BITS" == "32" ]]; then
    wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86.sh -O miniconda.sh ;
fi
chmod +x miniconda.sh && ./miniconda.sh -b

export PATH=/home/travis/miniconda2/bin:$PATH

ls 

ls ..

echo $PATH

conda update --yes conda

# Configure the conda environment and put it in the path using the
# provided versions
conda create -n anuga_env -c conda-forge --yes python=$PYTHON_VERSION pip numpy scipy cython netcdf4 nose matplotlib gdal dill future gitpython

source activate anuga_env
pip install mpi4py 

# python 2.6 doesn't have argparse by default
if [[ "$PYTHON_VERSION" == "2.6" ]]; then conda install --yes argparse; fi

export GDAL_DATA=`gdal-config --datadir`;

# Install more software to deal with geographical projections
#pip install pyproj

# Install pypar if parallel set
# if [[ "$ANUGA_PARALLEL" == "mpich2" || "$ANUGA_PARALLEL" == "openmpi" ]]; then
#     git clone https://github.com/daleroberts/pypar.git;
#     pushd pypar;
#     python setup.py install;
#     popd;
# fi

# Useful for debugging any issues with conda
conda info -a


#########################################################
# Build and install anuga

python setup.py build
python setup.py install
