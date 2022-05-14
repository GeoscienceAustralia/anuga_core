#!/bin/bash

# License: 3-clause BSD


set -e

echo "#==========================="
echo "# Install packages via apt"
echo "#==========================="
sudo apt-get update -q
sudo apt-get install git wget


echo "#==========================="
echo "# Install miniconda"
echo "#==========================="
wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3.sh


echo "#==========================="
echo "# Create conda environment anuga_env"
echo "#==========================="
conda update --yes conda
conda create -n anuga_env -c conda-forge --yes python=3.8 pip numpy scipy cython netcdf4 \
    pytest nose matplotlib gdal dill future gitpython mpi4py pytz utm Pmw pymetis meshpy
conda activate anuga_env 


echo "#==========================="
echo "# You should now install anuga"
echo "# from the anuga_core directory"
echo "# via"
echo "# conda activate anuga_env"
echo "# cd anuga_core"
echo "# pip install -e ."
echo "# and then test via"
echo "# pytest"
echo "#==========================="

