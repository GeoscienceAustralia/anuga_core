

******************
ANUGA Installation
******************

.. contents::


Introduction
============

ANUGA_ is a python package with some C extensions (and an optional fortran 
extension). At present ANUGA has only been run and tested using python 2.x.
We recommend python 2.7.
We are currently porting ANUGA to Python 3.x


Dependencies
============

AnuGA requires the following packages:

* `numpy <http://numpy.scipy.org/>`_
* `scipy <http://scipy.org/>`_
* `matplotlib <http://matplotlib.sourceforge.net/>`_
* `gdal <http://gdal.org/>`_
* `netcdf <http://www.unidata.ucar.edu/software/netcdf/>`_
* `nose <http://nose.readthedocs.org/en/latest/>`_
* `dill <https://dill.readthedocs.io/>`_
* A C compiler (preferably GCC or TDM-GCC_ MinGW_ on Windows)  # Possibly deprecated by the use of Cython?
* `Cython <https://cython.org/>`
* `future <https://pypi.org/project/future` # Possible superfluous if we ditch backward compatibility with Python 2.7


Installing the latest DEVELOPMENT version on Ubuntu
===================================================

AnuGA is developed on Ubuntu. The preferred way to install the dependencies is 
to use the standard ubuntu ``apt-get`` method. 

We suggest installing the latest version of ANUGA_ from Github_.
We try to maintain the *master* branch stable and
`passing all tests <https://travis-ci.org/GeoscienceAustralia/anuga_core/branches>`_,
so it should be safe to use.

Follow these instructions to
`Install ANUGA on Ubuntu (using apt-get) <https://github.com/GeoscienceAustralia/anuga_core/wiki/Install-ANUGA-on-Ubuntu-(using-apt-get)>`_




Alternative Ubuntu Install
==========================

An alternative is to install the dependencies using the Anaconda_ or the Miniconda_ Python 
distributions by `Continuum Analytics`_.

Miniconda_ has the advantage of allowing you to create multiple 
python environments and is particularly 
useful if you want to keep multiple versions of AnuGA.

Both Anaconda_ and Miniconda_ do not require administrative rights 
to your computer and do not interfere with the Python installed 
in your system.

Follow these instructions to 
`Install ANUGA on Ubuntu (using Miniconda) <https://github.com/GeoscienceAustralia/anuga_core/wiki/Install-ANUGA-on-Ubuntu-(using-Miniconda)>`_



Windows 10 Install using 'Ubuntu on Windows'
==========================

Starting from Windows 10, it is possible to run an Ubuntu Bash console from Windows. This can greatly simplify the install for Windows users. You'll still need administrator access though. 

Follow the instructions 
`Install ANUGA on Window 10 using Ubuntu for Windows <https://github.com/GeoscienceAustralia/anuga_core/wiki/Install-ANUGA-on-Windows-10-using-Ubuntu-on-Windows'>`_



Native Windows Installation
===============================

We have successfully installed AnuGA 'natively' on windows using Gohlke Binaries and using Miniconda. 
At present we recommend using the Miniconda. 

Follow the instructions  
`install ANUGA on Windows using Miniconda
<https://github.com/GeoscienceAustralia/anuga_core/wiki/Install-ANUGA-on-Windows-(using-Miniconda)>`_

Alternatively follow the instructions 
`install ANUGA on Windows using the Gohlke Binaries
<https://github.com/GeoscienceAustralia/anuga_core/wiki/Install-ANUGA-on-Windows-using-Gohlke-Binaries>`_






GCC dependency for Windows users
--------------------------------

Unfortunately, the ``gcc`` compiler MinGW_ included in Anaconda or 
installable via Miniconda_ doesn't have OpenMP_ support. This is required to compile
some extension modules in ANUGA (those that have multi-threaded parallel code).

We suggest that you download and install the version of MinGW_ provided by TDM-GCC_
**after** you've installed Anaconda and **before** you install ANUGA.
Don't forget to mark the ``openmp`` and ``gfortran`` options in the "Choose Components" part of
the installation. See this `excellent documentation for Windows users`_
(they even have screenshots!). The same applies if you are using Miniconda_.





.. note::

    The Windows installer from older versions is no longer supported.
    

    
.. _AnuGA: http://anuga.anu.edu.au/ 
.. _install git: http://git-scm.com/
.. _Github: https://github.com/GeoscienceAustralia/anuga_core/
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
