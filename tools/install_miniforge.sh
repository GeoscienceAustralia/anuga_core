#!/bin/bash


PY=${PY:-"3.10"}

set -e 

SCRIPT=$(realpath "$0")
SCRIPTPATH=$(dirname "$SCRIPT")
ANUGA_CORE_PATH=$(realpath "$SCRIPTPATH/..")


echo "#==========================="
echo "# Install miniforge3"
echo "#==========================="
cd $HOME

if [ -f "$HOME/Miniforge3.sh" ]; then
     echo "Miniforge3.sh already exists."
else
     echo "Miniforge3.sh does not exist. Downloading..."
     wget -O "$HOME/Miniforge3.sh" "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
fi


if [ -d "$HOME/miniforge3" ]; then
     echo "miniforge3 already exists."
else
     echo "Miniforge does not exist. Installing from Miniforge3.sh..."
     bash Miniforge3.sh
fi


echo "#==============================================="
echo "# create conda/mamba environment anuga_env_${PY}"
echo "#==============================================="
echo "..."
mamba create -n anuga_env_${PY} --yes python=${PY} compilers pip numpy scipy cython netcdf4 \
     nose matplotlib gdal dill gitpython mpi4py utm Pmw pymetis meshpy pytest pyproj affine \
     meson-python meson ninja xarray future

echo "#======================================"
echo "# activate environment anuga_env_${PY}"
echo "#======================================"
echo "..."
source ./miniforge3/bin/activate anuga_env_${PY}


echo "#================================================================"
echo "# Installing anuga from the ${ANUGA_CORE_PATH} directory"
echo "#================================================================"
echo "..."

cd ${SCRIPTPATH}
cd ..
pip install .
echo " "

echo "#==========================="
echo "# Run unittests"
echo "#==========================="
echo " "

cd ..
pytest -q --disable-warnings --pyargs anuga

echo " "
echo "#=================================================================="
echo "# Congratulations, Looks like you have successfully installed anuga"
echo "#=================================================================="

echo "#=================================================================="
echo "# To use anuga you must activate the python environment anuga_env_${PY}"
echo "# that has just been created. Run the command"
echo "# "
echo "# source ~/miniforge3/bin/activate anuga_env_${PY}"
echo "# "
echo "#=================================================================="

echo "#=================================================================="
echo "# NOTE: If you run the command"
echo "# "
echo "# mamba init"
echo "# "
echo "# (which will change your .bashrc file) "
echo "# then in new terminals you will be able to use "
echo "# the mamba command"
echo "# "
echo "# mamba activate anuga_env_${PY}"
echo "# "
echo "# to activate the anuga_env_${PY} environment"
echo "#=================================================================="
