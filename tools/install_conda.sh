#!/bin/bash
# This script is meant to be called by the "install" step defined in
# .travis.yml. See http://docs.travis-ci.com/ for more details.
# The behavior of the script is controlled by environment variabled defined
# in the .travis.yml in the top level folder of the project.

# License: 3-clause BSD


set -e


[ -z "$PYTHON_VERSION" ] && PYTHON_VERSION="2.7"


sudo apt-get update -q
sudo apt-get install gfortran

##########################################################
# Setup various versions of MPI
if [[ "$ANUGA_PARALLEL" == "mpich2" ]]; then
    sudo apt-get -y install mpich2;
fi

if [[ "$ANUGA_PARALLEL" == "openmpi" ]]; then
    sudo apt-get install -y libopenmpi-dev openmpi-bin;
fi


##########################################################

# Deactivate the travis-provided virtual environment and setup a
# conda-based environment instead
deactivate || echo "deactivate failed"

# Use the miniconda installer for faster download 
# install of conda itself
wget http://repo.continuum.io/miniconda/Miniconda-3.8.3-Linux-x86_64.sh \
    -O miniconda.sh
chmod +x miniconda.sh && ./miniconda.sh -b

export PATH=$HOME/miniconda/bin:$PATH
conda update --yes conda

# Configure the conda environment and put it in the path using the
# provided versions
conda create -n anuga_env --yes python=$PYTHON_VERSION pip numpy scipy netcdf4 \
    nose matplotlib
source activate anuga_env

# python 2.6 doesn't have argparse by default
if [[ "$PYTHON_VERSION" == "2.6" ]]; then conda install --yes argparse; fi

if [[ "$PYTHON_VERSION" == "2.7" ]]; then
    conda install --yes -c pingucarsti gdal
fi

if [[ "$PYTHON_VERSION" == "2.6" ]]; then
    conda install --yes  gdal geos
fi

export GDAL_DATA=`gdal-config --datadir`;

# Install more software to deal with geographical projections
pip install pyproj

# Install pypar if parallel set
if [[ "$ANUGA_PARALLEL" == "mpich2" || "$ANUGA_PARALLEL" == "openmpi" ]]; then
    git clone https://github.com/daleroberts/pypar.git;
    pushd pypar;
    python setup.py install;
    popd;
fi

# Useful for debugging any issues with conda
conda info -a

########################################################
if [[ "$COVERAGE" == "--coverage" ]]; then
    pip install coverage coveralls
fi

########################################################
# build and install anuga

python setup.py build
python setup.py install
