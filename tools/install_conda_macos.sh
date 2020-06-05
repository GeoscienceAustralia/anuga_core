#!/bin/bash
# This script is meant to be called by the "install" step defined in
# .travis.yml. See http://docs.travis-ci.com/ for more details.
# The behavior of the script is controlled by environment variabled defined
# in the .travis.yml in the top level folder of the project.

# License: 3-clause BSD


set -e


PYTHON_VERSION=${PYTHON_VERSION:-"2.7"}
ANUGA_BITS=${ANUGA_BITS:-"64"}
ANUGA_PARALLEL=${ANUGA_PARALLEL:-"openmpi"}


brew update
brew install openmpi

# if [[ "$ANUGA_PARALLEL" == "mpich2" || "$ANUGA_PARALLEL" == "openmpi" ]]; then
#     brew install openmpi
# fi


# Install miniconda
wget https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh -O miniconda.sh ;

chmod +x miniconda.sh && ./miniconda.sh -b

export PATH=/Users/travis/miniconda2/bin:$PATH

# This leads to an error at the moment (20200605)
#conda update --yes conda

# Configure the conda environment and put it in the path using the
# provided versions

#conda create -n anuga_env -c conda-forge --yes python=$PYTHON_VERSION pip numpy scipy cython netcdf4 \
#    nose matplotlib gdal dill
    
conda create -n anuga_env --yes python=2.7.13 gdal=2.2.2 pip nose numpy scipy netcdf4 matplotlib dill cython future

source activate anuga_env
pip install mpi4py

# Install pypar if parallel set
#if [[ "$ANUGA_PARALLEL" == "mpich2" || "$ANUGA_PARALLEL" == "openmpi" ]]; then
#    git clone https://github.com/daleroberts/pypar.git;
#    pushd pypar;
#    python setup.py install;
#    popd;
#fi


# python 2.6 doesn't have argparse by default
if [[ "$PYTHON_VERSION" == "2.6" ]]; then conda install --yes argparse; fi

export GDAL_DATA=`gdal-config --datadir`;

# Useful for debugging any issues with conda
conda info -a


#########################################################
# Build and install anuga

python setup.py build
python setup.py install
