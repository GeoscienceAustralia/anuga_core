

******************
ANUGA Installation
******************

.. contents::


Introduction
============

ANUGA_ is a python package with some C extensions (and an optional fortran 
extension). This version of ANUGA is run and tested using python 3.8.x


Dependencies
============

AnuGA requires the following python packages:

  numpy scipy matplotlib nose cython netcdf4 dill future gitpython gdal pyproj pymetis triangle Pmw mpi4py ipython



Installing the latest DEVELOPMENT version on Ubuntu 20_04
===================================================

AnuGA is developed on Ubuntu. The preferred way to install the dependencies is 
to use a combination of the standard ubuntu ``apt`` method and python pip install.

From your home directory run the following commands which will download anuga to a directory `anuga_core`, install dependencies, install anuga and run the unit tests::

    git clone https://github.com/anuga-community/anuga_core.git
    bash anuga_core/tools/install_ubuntu_20_04.sh

Note: This will set ``python``  as ``python3`` and part of the bash shell will run as sudo so will ask for a password. If you like you can run the package installs manually, run the commands in the script ``anuga_core/tools/install_ubuntu_20._04.sh``

You should check the installation by running the unit tests via::

  cd anuga_core
  python runtests.py
  

Alternative Ubuntu Install with conda
==========================

An alternative is to install the dependencies using the Anaconda_ or the Miniconda_ Python 
distributions by `Continuum Analytics`_.

Miniconda_ has the advantage of allowing you to create multiple 
python environments and is particularly 
useful if you want to keep multiple versions of AnuGA.

Both Anaconda_ and Miniconda_ do not require administrative rights 
to your computer and do not interfere with the Python installed 
in your system.

Folllow these steps::

    

    sudo apt-get update -q
    sudo apt-get install gfortran git wget
    sudo apt-get install -y libopenmpi-dev openmpi-bin
    
Download and install `Miniconda`::

    wget http://repo.continuum.io/miniconda3/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh
    
Create `anuga_env` conda environment::

    conda update --yes conda
    conda create -n anuga_env -c conda-forge --yes python pip numpy scipy cython netcdf4 nose matplotlib gdal dill future gitpython
    conda activate anuga_env
    pip install mpi4py triangle Pmw pymetis
    
Download, install and test `anuga`::

    git clone https://github.com/anuga-community/anuga_core.git
    cd anuga_core
    python setup.py install
    python runtests.py
    

Windows 10 Install using 'Ubuntu on Windows'
==========================

Starting from Windows 10, it is possible to run an Ubuntu Bash console from Windows. This can greatly simplify the install for Windows users. You'll still need administrator access though. First install an ubuntu 20_04 subsystem. Then just use the whichever ubuntu install described above. 




Native Windows Installation uding Miniconda
===============================

We have installed `anuga` on `windows` using miniconda.  

Run the following powershell instructions to download miniconda and the MPI files (for parallel runs)::

    Start-FileDownload "https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe" C:\Miniconda.exe; echo "Finished downloading miniconda"
    Start-FileDownload "https://download.microsoft.com/download/A/E/0/AE002626-9D9D-448D-8197-1EA510E297CE/msmpisetup.exe" C:\msmpisetup.exe; echo "Finished downloading msmpi"
    Start-FileDownload "https://download.microsoft.com/download/A/E/0/AE002626-9D9D-448D-8197-1EA510E297CE/msmpisdk.msi" C:\msmpisdk.msi; echo "Finished downloading msmpisdk"
    
And then the following cmd instructions::

    msiexec.exe /i "C:\msmpisdk.msi" /qn
    C:\msmpisetup.exe -unattend
    C:\Miniconda.exe /S /D=C:\Py
    C:\Py\Scripts\activate.bat
    set PATH=%PATH%;"C:\Program Files\Microsoft MPI\bin"
    conda config --set always_yes yes
    conda update conda
    conda install python=3.7 gdal nose numpy cython scipy netcdf4 matplotlib dill future gitpython
    pip install Pmw
    conda install -c msys2 libpython m2w64-toolchain
    pip install mpi4py triangle
    python setup.py install

