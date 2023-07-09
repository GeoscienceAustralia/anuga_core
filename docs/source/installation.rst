Installation
============

.. contents::


Introduction
------------

ANUGA is a python package with some C extensions (and an optional fortran 
extension). This version of ANUGA is run and tested using python 3.7 - 3.10


Dependencies
------------

ANUGA requires python 3.X (X>6) and the following python packages:

.. code-block::

  numpy scipy matplotlib pytest cython netcdf4 dill future gdal \
  pyproj pymetis Pmw mpi4py ipython meshpy Pmw pymetis utm pyproj affine

and if using `python <= 3.8`

.. code-block::
  
  backport.zoneinfo 


ANUGA is developed on Ubuntu and so we recommend Ubuntu as your production environment using 
`Miniconda` or `Miniforge` python environments. 
(ANUGA can be installed on MacOS and Windows also using `Miniconda` or `MiniForge`) 

.. _Install MiniForge:

Install MiniForge (Ubuntu)
--------------------------

A clean way to install the dependencies necessary for ANUGA is to use the Anaconda, 
or the Miniconda Python distributions by Continuum Analytics. 

Using a `conda` installation has the advantage of allowing you to easily create multiple 
python environments and is particularly 
useful if you want to keep multiple versions of ANUGA

Indeed the most stable install is via the `conda-forge` channel
which is easily available using `Miniforge`. In particular the installation of 
the `gdal` and `mpi4py` modules are more stable using `Miniforge`. We recommend  using `Miniforge`. 

These conda environments do not require administrative rights 
to your computer and do not interfere with the Python installed in your system. 

So first we should install the latest version of `Miniforge` from  https://github.com/conda-forge/miniforge or
use, for instance, `wget` to download the latest version via:

.. code-block:: bash

    wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"


.. note::
    
    If you don't have `wget` you can install it via: 

    .. code-block:: bash

        sudo apt-get update -q
        sudo apt-get install wget


Install `Miniforge`
~~~~~~~~~~~~~~~~~~~ 

Run the installation script:

.. code-block:: bash

    bash Miniforge3.sh

and then activate `miniconda` by running 

.. code-block:: bash

    source miniconda3/bin/activate

.. note::

    During the `Miniforge` installation you will be asked to accept the licence 
    (essentially apache 2.0) and whether to run `conda init` to change your `.bashrc` file to allow activation of the 
    base conda environment when opening a new terminal.
    
    If you choose not to run `conda init` you will need to run the 
    following command every time to activate `miniforge`

    .. code-block:: bash

        source miniconda3/bin/activate 


Once `Miniforge` is installed and activated we can now create an environment to run ANUGA. 


Install ANUGA using MiniForge (Ubuntu)
--------------------------------------

.. note::
    If you want the most recent update of ANUGA or intend to develop ANUGA code you 
    should install ANUGA from source 
    (see section `Install ANUGA from source using MiniForge`_ )

Once you have a working `Miniforge` installation (as described in the previous 
section `Install MiniForge`_ ) you are ready to install a prebuilt 
version of ANUGA from `conda-forge`. 

It is always recommended that you create a separate `conda` environment for 
your ANUGA installation. 

So first create a python 3.9 conda environment called `anuga_env` (or what ever name you like):

.. code-block:: bash

    conda create -n anuga_env python=3.9 anuga mpi4py
    conda activate anuga_env

Note we have also installed `mpi4py` to allow anuga to run in parallel. 
On some systems you may need to manually install `mpi4py` to match the version of `mpi` you are using.


This has setup and activated a `conda` environment `anuga_env` which is using python 3.9. (ANUGA has be tested on 3.7, 3.8. 3.9.)    

We are now ready to use ANUGA. 

You can test your installation via:

.. code-block:: bash

    conda activate anuga_env
    pytest --pyargs anuga


.. _Install ANUGA from source using MiniForge:

Install ANUGA from source (e.g. as a developer) using MiniForge (Ubuntu)
------------------------------------------------------------------------

If you want to use the very latest version of ANUGA (or develop ANUGA code) then you need
to download the `anuga_core` repository from `github` and then `pip` install ANUGA from the source. These steps will require that the following packages are installed: git, pip, conda (via miniforge described at the beginning of this document).


First install the latest version of `Miniforge` as described in section `Install MiniForge`_.

Now we need to download the ANUGA source code from `github`

.. code-block:: bash

    git clone https://github.com/anuga-community/anuga_core.git

This creates a directory `anuga_core`.

.. note::

    If you want to also contribute to the code base, you must have a GitHub account and setup authentication from your developer workstation to GitHub as per these instructions:        https://docs.github.com/en/authentication/managing-commit-signature-verification. The command to clone ANUGA as a developer is then 

    .. code-block:: bash

        git clone git@github.com:anuga-community/anuga_core.git

Now create and activate a `conda` environment with ANUGA's current dependencies as 
defined in the file `environment.yml`

.. code-block:: bash

    cd anuga_core
    conda env create -n anuga_env -f environment.yml
    conda activate anuga_env

and finally install ANUGA. Do a standard `pip` install

.. code-block:: bash

    pip install .

.. note::

    If you intend to develop ANUGA code then you should install ANUGA to be "editable". I.e.:

    .. code-block:: bash

        pip install -e .

    In this case the installation is "inplace" and "editable". You will be able to change and 
    develop code in the `anuga_core` directories. Note that if you change any `cython` or `C` 
    code you will need to run `pip install -e .` again for your changes to take effect.

.. note::

    You may need to install a compiler to complete the `pip install`. 
    You can use the system compilers or use `conda` to install compilers as such (for Linux and OSX):

    .. code-block:: bash

        conda install compilers

    or for win32:

    .. code-block:: bash

        conda install m2w64-gcc libpython 
 

Finally it is sensible to test the installation.

.. code-block:: bash

    pytest --pyargs anuga


Updating
~~~~~~~~

From time to time you might like to update your version of anuga to the latest version on 
github. You can do this by going to the `anuga_core` directory and `pulling` the latest
version and then reinstalling via the following commands:
 
.. code-block:: bash

  conda activate anuga_env
  cd anuga_core
  git pull
  pip install .

And finally check the new installation by running the unit tests via:

.. code-block:: bash

  pytest --pyargs anuga 


Installing on Ubuntu_20.04 and 22.04 using script `install_ubuntu.sh`
---------------------------------------------------------------------

For Ubuntu 20.04 and 22.04 you can install ANUGA and its dependencies into a python virtual environment via 
a simple `bash` script.

First from your home directory run the following command which will download anuga 
to a directory `anuga_core`:

.. code-block:: bash

    git clone https://github.com/anuga-community/anuga_core.git

Then the following will install dependencies, install anuga and run the unit tests:

.. code-block:: bash

    bash anuga_core/tools/install_ubuntu.sh

Note: This script will only work for Ubuntu 20_04 and 22_04.

Note: Part of the bash shell will run as 
sudo so will ask for a password. If you don't like this, you can run the package installs manually, 
see the commands in the scripts ``anuga_core/tools/install_ubuntu_20._04.sh`` 
or ``anuga_core/tools/install_ubuntu_22._04.sh`` as appropriate.  

This script also creates a python3 virtual environment `anuga_env`. You should activate this 
virtual environment when working with ANUGA, via the command:

.. code-block:: bash

    source ~/anuga_core/anuga_env/bin/activate

You might like to add this command to your `.bashrc` file to automatically activate this 
python environment. 

Updating
~~~~~~~~

From time to time you might like to update your version of anuga to the latest version on 
github. You can do this by going to the `anuga_core` directory and `pulling` the latest
version and then reinstalling via the following commands:
 
 Activate the environment if necessary:

.. code-block:: bash

    source ~/anuga_core/anuga_env/bin/activate

Then update ANUGA to latest version:

.. code-block:: bash

  cd anuga_core
  git pull
  pip install -e .

And finally check the new installation by running the unit tests via:

.. code-block:: bash

  pytest --pyargs anuga 
      

Windows 10 Install using 'Ubuntu on Windows'
--------------------------------------------

Starting from Windows 10, it is possible to run an Ubuntu Bash console from Windows. 
This can greatly simplify the install for Windows users. 
You'll still need administrator access though. First install an ubuntu 20_04 subsystem. 
Then just use your preferred ubuntu install described above. 



Windows Installation using MiniForge
------------------------------------

We have installed ANUGA on `windows` using miniforge.  

You can download MiniForge manually 
from the MiniForge site https://github.com/conda-forge/miniforge:

Alternatively you can download and install miniforge via CLI commands:

Run the following powershell instruction to download miniforge. 

.. code-block:: bash

    Start-FileDownload "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Windows-x86_64.exe" C:\Miniforge.exe; 
  
From a standard `cmd` prompt then install miniconda via:

.. code-block::  bash

    C:\Miniconda.exe /S /D=C:\Py
    C:\Py\Scripts\activate.bat
    
Install conda-forge packages:

.. code-block:: bash

    conda create -n anuga_env python=3.8  anuga mpi4py
    conda activate anuga_env
    
You can test your installation via:

.. code-block:: bash

    pytest --pyargs anuga

    
Installing GDAL on Ubuntu using apt and pip
-------------------------------------------

ANUGA can be installed using the python provided by the Ubuntu system and using `pip`. 

First set up a python virtual environment and activate  via:

.. code-block:: bash

    python3 -m venv anuga_env
    source anuga_env/bin/activate

A complication arises when installing  the `gdal` package. 
First install the gdal library, via:

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

Now we complete the installation of ANUGA simply by:

.. code-block:: bash

    pip install anuga

If you obtain errors from `pip` regarding "not installing dependencies", it seems that that can be fixed by just 
running the `pip install anuga` again.
