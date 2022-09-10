#!/bin/bash

# License: 3-clause BSD


#!/bin/bash

# License: 3-clause BSD


set -e


echo "#==========================="
echo "# Install miniforge"
echo "#==========================="
if [[ -f "Miniforge3.sh" ]]; then 
  echo "Miniforge3.sh already exists"
else
  echo "Downloading Miniforge3.sh"
  #wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
  #bash Miniforge3.sh
fi


echo "#==========================="
echo "# Create conda environment anuga_env"
echo "#==========================="
conda update --yes conda
conda create -n anuga_env --yes python=3.10 gxx pip wheel scipy numpy cython netcdf4 pytest nose matplotlib gdal dill future gitpython utm Pmw pymetis meshpy mpi4py
conda activate anuga_env 


echo "#==========================="
echo "# Installing anuga from the anuga_core directory"
echo "# and then run unittests"
echo "#==========================="

cd "$(dirname "${BASH_SOURCE[0]}")"/..
pip install -e .
python runtests.py -n

echo "#================================================"
echo "# To use anuga you must activate the conda"
echo "# python environment anuga_env that has been"
echo "# created by running:"
echo "# conda activate anuga_env"
echo "#================================================"
