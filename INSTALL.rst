

Installation
============

.. contents::

Introduction
------------

AnuGA is a python package with some C extensions (and an optional fortran 
extension). At present AnuGA has only been run and tested using python 2.x.
We recommend python 2.7.  

If you have a python 2.7 environment with gcc and gfortran support, 
then you can use easy_install to install the latest release 
version of AnuGA. 


Installing the dependencies
---------------------------

AnuGA requires the following packages:

* `numpy <http://numpy.scipy.org/>`_
* `scipy <http://scipy.org/>`_
* `matplotlib <http://matplotlib.sourceforge.net/>`_
* `gdal <http://gdal.org/>`_
* `nose <http://nose.readthedocs.org/en/latest/>`_
* A C compiler (preferably GCC or MinGW_ on Windows)

The easiest and **preferred** way to get all dependencies in the latest
version is using the Anaconda_ or the Miniconda_ Python 
distributions by `Continuum Analytics`_.
Miniconda_ allows you to create multiple python environments and is particularly 
useful if you want to keep multiple versions of AnuGA.
It does not require administrative rights to your computer and doesn't
interfere with the Python installed in your system.
For Windows users Anaconda_ even comes with MinGW_ so you don't have to worry about
the many, many, many issues of compiling under Windows.


Anaconda
++++++++

Once you have downloaded and installed Anaconda_,
open a terminal (or ``cmd.exe`` on Windows) and run::

    conda install pip nose numpy scipy netcdf4 matplotlib 
    conda install -c pingucarsti gdal 
    
and setup GDAL_DATA environment variable:

    export GDAL_DATA=`gdal-config --datadir` 
    
    
Miniconda
+++++++++

Once you have downloaded and installed Miniconda_, 
open a terminal (or ``cmd.exe`` on Windows), create 
a specific environment for AnuGA, by running::

    conda create -n anuga_env python=2.7
    source activate anuga_env
    
    conda install pip nose numpy scipy netcdf4 matplotlib 
    conda install -c pingucarsti gdal 
    
and setup GDAL_DATA environment variable:

    export GDAL_DATA=`gdal-config --datadir` 
    


Extra dependencies for Windows users
++++++++++++++++++++++++++++++++++++

Unfortunately, the ``gcc`` compiler included in Anaconda MinGW_
doesn't have OpenMP_ support. This is required to compile
some extension modules in AnuGA that have multi-threaded parallel code.

You'll have download and install TDM-GCC_
**after** you've installed Anaconda and **before** you install AnuGA.
Don't forget to mark the ``openmp`` option in the "Choose Components" part of
the installation. See this `excellent documentation for Windows users`_
(they even have screenshots!).

Installing AnuGA
----------------

After you've installed the dependencies you can proceed to install AnuGA
using pip_.
Open a terminal (or ``cmd.exe`` on Windows) and run::

    pip install anuga

and that's it!

If you already have AnuGA installed and want to **upgrade** to a newer
version, use::

    pip install anuga --upgrade

To uninstall simply run::

    pip uninstall anuga


.. note::

    The Windows installer from older versions is no longer supported.

Installing the latest development version
-----------------------------------------

If you want the very latest code and features,
you can install AnuGA directly from Github_.
We try to maintain the *master* branch stable and
`passing all tests <https://travis-ci.org/stoiver/anuga_core/branches>`__,
so it should be safe to use.

First, you'll need to `install git`_.
Then, open a terminal and run::

    git clone --depth=50 --branch=master git://github.com/stoiver/anuga_core.git 

This will fetch the source code from Github_
and place it in a folder called ``anuga_core`` in the directory where you ran the
command.
Then, just ``cd`` into the directory and run ``pip``::

    cd anuga_core
    pip install --upgrade .
    
Testing the install
-------------------


From the source directory run ``python runtests.py``

    python runtests.py
    

If you get an error message or weird result,
please write to the `mailing list`_.
To make it easier for us to debug you problem, please include the following
information:

* Operating system
* Python distribution (Anaconda_, PythonXY_, `ETS/Canopy`_, own install)
* Python version (2.6, 2.7, 3.3, 3.4, etc)
* The script you ran (and gave you an error/weird result)
* The error message (the part that says ``Traceback: ...``) or result (figure,
  numbers, etc)
    
    
.. _install git: http://git-scm.com/
.. _Github: https://github.com/stoiver/anuga_core/
.. _Python: http://www.python.org/
.. _pip: http://www.pip-installer.org
.. _MinGW: http://www.mingw.org/
.. _mailing list: anuga-user@lists.sourceforge.net
.. _Continuum Analytics: http://continuum.io/
.. _Anaconda: http://continuum.io/downloads
.. _Miniconda: http://conda.pydata.org/miniconda.html
.. _PythonXY: http://code.google.com/p/pythonxy/
.. _ETS/Canopy: http://code.enthought.com/projects/index.php
.. _OpenMP: http://openmp.org/
.. _TDM-GCC: http://tdm-gcc.tdragon.net/
.. _excellent documentation for Windows users: http://docs-windows.readthedocs.org/en/latest/devel.html#mingw-with-openmp-support

