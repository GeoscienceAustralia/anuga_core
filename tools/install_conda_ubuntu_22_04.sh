#!/bin/bash

# License: 3-clause BSD


set -e

echo "#==========================="
echo "# Install packages via apt"
echo "#==========================="
sudo apt-get update -q
sudo apt-get install git wget


echo "#==========================="
echo "# Install miniforge"
echo "#==========================="
wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3.sh


echo "#==========================="
echo "# Create conda environment anuga_env"
echo "#==========================="
conda update --yes conda
conda create -n anuga_env python=3.10 gxx pip wheel scipy numpy cython netcdf4 pytest nose matplotlib gdal dill future gitpython utm Pmw pymetis meshpy mpi4py
conda activate anuga_env 


echo "#==========================="
echo "# Installing anuga from the anuga_core directory"
echo "# and then run unittests"
echo "#==========================="

cd "$(dirname "${BASH_SOURCE[0]}")"/..
pip install -e .
pytest --pyargs anuga


