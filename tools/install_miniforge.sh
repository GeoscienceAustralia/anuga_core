#!/bin/bash
#
# Run this script as follows:
#
# bash /path/to/anuga_core/tools/install_miniforge.sh
#
# This script will install miniforge3 and create a conda environment called
# anuga_env_${PY} where PY is the python version you want to use.
# The script will then install the anuga package from the anuga_core directory
# and run the unittests.
#
# By default a python version of 3.10 will be installed. If you want to install
# a different version of python, set the PY environment variable before running
# the script. For example, to install python 3.9 run the script as follows:
#
# PY=3.9 bash /path/to/anuga_core/tools/install_miniforge.sh
#
# The script will install python 3.9 and create the anuga_env_3.9 environment.


PY=${PY:-"3.10"}

set -e 

SCRIPT=$(realpath "$0")
SCRIPTPATH=$(dirname "$SCRIPT")
ANUGA_CORE_PATH=$(realpath "$SCRIPTPATH/..")


#test PY>3.8 and <3.13
if [[ "$PY" =~ ^3\.(1[0-3]|[9])$ ]]; then
     echo "Requested python version is $PY"
     echo " "
else
    echo "Python version must be greater than 3.8 and less than 3.14"
    exit 1
fi


echo "#==========================="
echo "# Install miniforge3"
echo "#==========================="
cd $HOME




if [ -d "$HOME/miniforge3" ]; then
     echo "miniforge3 seems to already exist."
else
     echo "miniforge3 does not exist."
     echo "Installing from script Miniforge3.sh..."
     if [ -f "$HOME/Miniforge3.sh" ]; then
          echo "Running Miniforge3.sh first..."
     else
          echo "Miniforge3.sh does not exist. Downloading..."
          wget -O "$HOME/Miniforge3.sh" "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
     fi
     bash Miniforge3.sh -b
fi

echo " "
echo "#==============================================="
echo "# create conda/mamba environment anuga_env_${PY}"
echo "#==============================================="
echo "..."
./miniforge3/bin/mamba env create --file ${SCRIPTPATH}/../environments/environment_${PY}.yml

# ./miniforge3/bin/mamba create -n anuga_env_${PY} --yes python=${PY} compilers numpy scipy cython netcdf4 \
#      nose matplotlib gdal dill gitpython mpi4py utm Pmw pymetis meshpy pytest pyproj affine \
#      meson-python meson ninja xarray future pkg-config

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
