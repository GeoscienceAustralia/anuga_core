Installation
============

.. contents::



ANUGA_ is a python package with some C extensions (and an optional fortran 
extension). This version of ANUGA is run and tested using python 3.8.x


Dependencies
------------

ANUGA requires the following python packages:

  numpy scipy matplotlib nose cython netcdf4 dill future gitpython gdal pyproj pymetis triangle Pmw mpi4py ipython



Installing the latest version on Ubuntu
---------------------------------------

ANUGA is developed on Ubuntu. The preferred way to install the dependencies is 
to use a combination of the standard ubuntu ``apt`` method and python pip install.

From your home directory run the following commands which will download anuga to a directory `anuga_core`, install dependencies, install anuga and run the unit tests::

    git clone https://github.com/anuga-community/anuga_core.git
    sudo bash anuga_core/tools/install_ubuntu_20_04.sh

Note: This will set ``python``  as ``python3`` and part of the bash shell will run as sudo so will ask for a password. If you like you can run the package installs manually, run the commands in the script ``anuga_core/tools/install_ubuntu_20._04.sh``

You should now install and check the installation of anuga by running the unit tests via::

  cd anuga_core
  python setup.py install --user
  python runtests.py
  
Installing the latest version on Ubuntu as a developer
------------------------------------------------------
  
If you wish to install ANUGA and make changes to the code, the installation procedure is as above, but with the setup step as follows::

  python setup.py develop --user
  

Alternative Ubuntu Install with conda
-------------------------------------

An alternative is to install the dependencies using the Anaconda_ or the Miniconda_ Python 
distributions by `Continuum Analytics`_.

Miniconda_ has the advantage of allowing you to create multiple 
python environments and is particularly 
useful if you want to keep multiple versions of ANUGA.

Both Anaconda_ and Miniconda_ do not require administrative rights 
to your computer and do not interfere with the Python installed 
in your system. But it is necessary to install a few packages via `sudo apt-get` in particular 
a compiler and openmpi for parallel code. 

Follow these steps::

    sudo apt-get update -q
    sudo apt-get install gfortran git wget
    sudo apt-get install libopenmpi-dev openmpi-bin
    
Download and install `Miniconda` if you haven't already::

    wget http://repo.continuum.io/miniconda3/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh
    
Create `anuga_env` conda environment::

    conda update conda
    conda create -n anuga_env -c conda-forge python pip numpy scipy cython netcdf4 nose matplotlib gdal dill future gitpython
    conda activate anuga_env
    pip install mpi4py triangle Pmw pymetis
    
Download, install and test `anuga`::

    git clone https://github.com/anuga-community/anuga_core.git
    cd anuga_core
    python setup.py install
    python runtests.py
    

Windows 10 Install using 'Ubuntu on Windows'
--------------------------------------------

Starting from Windows 10, it is possible to run an Ubuntu Bash console from Windows. This can greatly simplify the install for Windows users. You'll still need administrator access though. First install an ubuntu 20_04 subsystem. Then just use your preferred ubuntu install described above. 




Windows Installation using Miniconda
------------------------------------

We have installed `anuga` on `windows` using miniconda.  

Run the following powershell instruction to download miniconda. You can also just download manually::

    Start-FileDownload "https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe" C:\Miniconda.exe; echo "Finished downloading miniconda"
  
  
From a standard `cmd` prompt then install miniconda via::

    C:\Miniconda.exe /S /D=C:\Py
    C:\Py\Scripts\activate.bat
    
Install conda and pip packages::

    conda install -c conda-forge python=3.7 gdal nose numpy cython scipy netcdf4 matplotlib dill future gitpython mpi4py
    pip install triangle Pmw pymetis
    conda install -c msys2 libpython m2w64-toolchain
    
Download ANUGA and install::

    git clone https://github.com/anuga-community/anuga_core.git
    cd anuga_core
    python setup.py install
    
And finally test the installation:: 

    python runtests.py
