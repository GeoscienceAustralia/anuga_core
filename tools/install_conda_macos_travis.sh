#!/bin/bash
# This script is meant to be called by the "install" step defined in
# .travis.yml. See http://docs.travis-ci.com/ for more details.
# The behavior of the script is controlled by environment variabled defined
# in the .travis.yml in the top level folder of the project.

# License: 3-clause BSD


set -e

#brew update
#brew install openmpi


# Install miniconda
wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh ;

chmod +x miniconda.sh && ./miniconda.sh -b

export PATH=/Users/travis/miniconda3/bin:$PATH

#conda update --yes conda

# Configure the conda environment and put it in the path using the
# provided versions
    
conda create -n anuga_env --yes python=3.7 gdal pip nose numpy scipy netcdf4 matplotlib dill cython future gitpython pytz

source activate anuga_env

pip install triangle
pip install Pmw
pip install pymetis

# Useful for debugging any issues with conda
conda info -a

#########################################################
# Build and install anuga

python setup.py build
python setup.py install
