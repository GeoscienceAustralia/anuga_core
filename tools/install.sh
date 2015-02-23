#!/bin/bash
# This script is meant to be called by the "install" step defined in
# .travis.yml. See http://docs.travis-ci.com/ for more details.
# The behavior of the script is controlled by environment variabled defined
# in the .travis.yml in the top level folder of the project.

# License: 3-clause BSD


set -e

sudo apt-get update -qq
sudo apt-get install gfortran

if [[ "$PARALLEL" == "mpich2" ]]; then
    sudo apt-get install mpich2;
    mpirun -np 4 pwd;
fi

if [[ "$DISTRIB" == "conda" ]]; then
    # Deactivate the travis-provided virtual environment and setup a
    # conda-based environment instead
    deactivate

    # Use the miniconda installer for faster download 
    # install of conda itself
    wget http://repo.continuum.io/miniconda/Miniconda-3.8.3-Linux-x86_64.sh \
	-O miniconda.sh
    chmod +x miniconda.sh && ./miniconda.sh -b
    export PATH=/home/travis/miniconda/bin:$PATH

    # Let conda install run with waiting for y/n reply
    conda config --set always_yes yes --set changeps1 no
    conda update --yes conda

    # Useful for debugging any issues with conda
    conda info -a

    # Configure the conda environment and put it in the path using the
    # provided versions
    conda create -n testenv --yes python=$PYTHON_VERSION pip numpy scipy netcdf4 \
	nose matplotlib
    source activate testenv

    if [[ "$PYTHON_VERSION" == "2.7" ]]; then
	conda install -c pingucarsti gdal;
	#conda install -c anaconda gdal;
    fi

    if [[ "$PYTHON_VERSION" == "2.6" ]]; then conda install gdal; fi

    export GDAL_DATA=`gdal-config --datadir`;

    # Install more software to deal with geographical projections
    pip install pyproj

    # python 2.6 doesn't have argparse by default
    if [[ "$PYTHON_VERSION" == "2.6" ]]; then conda install argparse; fi

elif [[ "$DISTRIB" == "ubuntu" ]]; then
    # Use standard ubuntu packages in their default version
    sudo apt-get install -qq python-scipy

    pip install nose
fi


# Install pypar if parallel
if [[ "$PARALLEL" == "mpich2" ]]; then
    svn checkout http://pypar.googlecode.com/svn/ pypar;
    pushd pypar/source;
    python setup.py install;
    popd;
fi


if [[ "$COVERAGE" == "true" ]]; then
    pip install coverage coveralls
fi

python setup.py install
