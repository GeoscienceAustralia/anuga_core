#!/bin/bash
# This script is meant to be called by the "install" step defined in
# .travis.yml. See http://docs.travis-ci.com/ for more details.
# The behavior of the script is controlled by environment variabled defined
# in the .travis.yml in the top level folder of the project.

# License: 3-clause BSD


set -e


PYTHON_VERSION=${PYTHON_VERSION:-"3.8"}
ANUGA_BITS=${ANUGA_BITS:-"64"}



sudo apt-get update -q
sudo apt-get install gfortran git wget

#sudo apt-get install -y libopenmpi-dev openmpi-bin;

##########################################################

# Deactivate the travis-provided virtual environment and setup a
# conda-based environment instead
deactivate || echo "deactivate failed"


echo "#==========================="
echo "# Install miniforge"
echo "#==========================="
wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3.sh

#export PATH=/home/travis/miniconda2/bin:$PATH

echo $PATH

echo "#==========================="
echo "# Create conda environment anuga_env"
echo "#==========================="

conda update --yes conda

# Configure the conda environment and put it in the path using the
# provided versions
conda create -n anuga_env --yes python=3.8 pip numpy scipy cython netcdf4 \
    nose matplotlib gdal dill future gitpython mpi4py backports.zoneinfo utm Pmw pymetis meshpy

source activate anuga_env

# Useful for debugging any issues with conda
conda info -a

echo "#==========================="
echo "# Installing anuga from the anuga_core directory"
echo "# and then run unittests"
echo "#==========================="

cd "$(dirname "${BASH_SOURCE[0]}")"/..
pip install -e .
pytest --pyargs anuga
