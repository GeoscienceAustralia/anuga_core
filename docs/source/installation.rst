Installation
============

.. contents::


Introduction
------------

ANUGA is a python package with some C extensions (and an optional fortran 
extension). This version of ANUGA is run and tested using python 3.8.x


Dependencies
------------

ANUGA requires python 3.X (X>6) and the following python packages:

.. code-block::

  numpy scipy matplotlib pytest cython netcdf4 dill future gitpython gdal pyproj pymetis triangle Pmw mpi4py pytz ipython  meshpy Pmw pymetis

ANUGA is developed on Ubuntu and so we recommend Ubuntu as your production environment
(though ANUGA can be installed on MacOS and Windows using `Miniconda` or `MiniForge`) 

Ubuntu Install with MiniForge3
------------------------------

A clean way to install the dependencies for ANUGA is to use Anaconda, 
or Miniconda Python distributions by Continuum Analytics, or the `conda-forge`
distribution Miniforge3. 

This installation has the advantage of allowing you to create multiple 
python environments and is particularly 
useful if you want to keep multiple versions of ANUGA

Indeed the most stable install is via the `conda-forge` channel
which is easily available using the Miniforge. In particular the installation of 
the `gdal` module is more stable. 

These conda environments do not require administrative rights 
to your computer and do not interfere with the Python installed in your system. 

But it is necessary to install a few packages via :code:`sudo apt-get`. 

Follow these steps:

.. code-block:: bash

    sudo apt-get update -q
    sudo apt-get install git wget
    
Download and install Miniforge if you haven't already:

.. code-block:: bash

    wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
    bash Miniforge3.sh
    
Create a conda environment `anuga_env` (or what ever name you like):

.. code-block:: bash

    conda update conda
    conda create -n anuga_env python=3.8 pip numpy scipy cython netcdf4 pytest nose matplotlib gdal dill future gitpython pytz mpi4py meshpy Pmw pymetis
    conda activate anuga_env

This will setup an environment using python 3.8. `anuga` has be tested on 3.7, 3.8. 3.9.    

Now Download, install and test `anuga`:

.. code-block:: bash

    git clone https://github.com/anuga-community/anuga_core.git
    cd anuga_core
    pip install -e .
    pytest

Remember, to use ANUGA you will have to activate the `anuga_env` environment 
via the command:

.. code-block:: bash
    
    conda activate anuga_env`

You might even like to set this up in your `.bashrc` file. 

Installing on Ubuntu using apt and pip
--------------------------------------

ANUGA can be installed using `pip`, but a complication arrise when installing 
the `gdal` package. 

First set up a python virtual environment and activate  via:

.. code-block:: bash

    python3 -m venv anuga_env
    course anuga_env/bin/activate

Now we first need to install the `gdal` python package. First install the 
gdal library, via:

.. code-block:: bash
    sudo apt-get install -y gdal-bin libgdal-dev

We need to ascertain the version of  `gdal` installed using the following command: 

.. code-block:: bash

    ogrinfo --version

THe version of `gdal` to install via `pip` should match the version of the library. 
For instance on Ubuntu 20.04 the previous command produces:

.. code-block:: bash

    GDAL 3.0.4, released 2020/01/28

So in this case we install the `gdal` python package as follows

.. code-block:: bash

    pip install gdal==3.0.4

Now we complete the installation of `anuga` simply by:

.. code-block:: bash

    pip install anuga

Installing on Ubuntu using apt and pip
---------------------------------------

You can install the `anuga` dependencies via a  combination of the 
standard ubuntu ``apt`` method and python pip install.

From your home directory run the following commands which will download anuga 
to a directory `anuga_core`, install dependencies, install anuga and run the unit tests:

.. code-block:: bash

    git clone https://github.com/anuga-community/anuga_core.git
    sudo bash anuga_core/tools/install_ubuntu_20_04.sh

Note: This will set ``python``  as ``python3`` and part of the bash shell will run as 
sudo so will ask for a password. If you like you can run the package installs manually, 
run the commands in the script ``anuga_core/tools/install_ubuntu_20._04.sh``

You should now install anuga: 
.. code-block:: bash

  cd anuga_core
  pip install -e .

And finally check the installation by running the unit tests via:
.. code-block:: bash

  pytest
      

Windows 10 Install using 'Ubuntu on Windows'
--------------------------------------------

Starting from Windows 10, it is possible to run an Ubuntu Bash console from Windows. 
This can greatly simplify the install for Windows users. 
You'll still need administrator access though. First install an ubuntu 20_04 subsystem. 
Then just use your preferred ubuntu install described above. 



Windows Installation using MiniForge
------------------------------------

We have installed `anuga` on `windows` using miniforge.  

You can download MiniForge manually 
from the MiniForge site https://github.com/conda-forge/miniforge:

Alternatively you can download and install miniforge via CLI commands:

Run the following powershell instruction to download miniforge. 

.. code-block:: bash

    Start-FileDownload "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Windows-x86_64.exe" C:\Miniforge.exe; echo "Finished downloading miniforge"
  
From a standard `cmd` prompt then install miniconda via:

.. code-block::  bash

    C:\Miniconda.exe /S /D=C:\Py
    C:\Py\Scripts\activate.bat
    
Install conda-forge packages:

.. code-block:: bash

    conda create -n anuga_env python=3.8 gdal pytest nose numpy cython scipy netcdf4 matplotlib dill future gitpython mpi4py meshpy Pmw pymetis
    conda activate anuga_env
    conda install libpython m2w64-toolchain
    
Download ANUGA and install:

.. code-block:: bash

    git clone https://github.com/anuga-community/anuga_core.git
    cd anuga_core
    cp windows_setup.cfg setup.cfg
    pip install -e .
    
And finally test the installation:

.. code-block:: bash

    
