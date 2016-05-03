

Installation
============

.. contents::

Introduction
------------

AnuGA_ is a python package with some C extensions (and an optional fortran 
extension). At present AnuGA has only been run and tested using python 2.x.
We recommend python 2.7.  

If you have a python 2.7 environment with gcc and gfortran support, 
then you can use pip to install the latest release 
version of AnuGA::

	pip install anuga

You might want to install the dependencies via ubuntu ``apt-get`` or ``conda install`` 
to download pre-compiled versions of the dependencies. See below. 


Dependencies
------------

AnuGA requires the following packages:

* `numpy <http://numpy.scipy.org/>`_
* `scipy <http://scipy.org/>`_
* `matplotlib <http://matplotlib.sourceforge.net/>`_
* `gdal <http://gdal.org/>`_
* `netcdf <http://www.unidata.ucar.edu/software/netcdf/>`_
* `nose <http://nose.readthedocs.org/en/latest/>`_
* A C compiler (preferably GCC or TDM-GCC_ MinGW_ on Windows)


Ubuntu Install
++++++++++++++

AnuGA is developed on Ubuntu. The preferred way to install the dependencies is 
to use the standard ubuntu ``apt-get`` method. 

We have a script in the ``tools`` directory,
`install_ubuntu.sh <https://github.com/GeoscienceAustralia/anuga_core/blob/master/tools/install_ubuntu.sh>`_
which when run from the ``anuga_core`` directory will install AnuGA and all the dependencies using ``apt-get``.

For extended instructions checkout the wiki page
`install ANUGA on Ubuntu <https://github.com/GeoscienceAustralia/anuga_core/wiki/Install-ANUGA-on-Ubuntu>`_



Alternative Ubuntu Install
++++++++++++++++++++++++++

An alternative is to install the dependencies using the Anaconda_ or the Miniconda_ Python 
distributions by `Continuum Analytics`_.

Miniconda_ has the advantage of allowing you to create multiple 
python environments and is particularly 
useful if you want to keep multiple versions of AnuGA.

Both Anaconda_ and Miniconda_ do not require administrative rights 
to your computer and do not interfere with the Python installed 
in your system.


Anaconda and Miniconda
~~~~~~~~~~~~~~~~~~~~~~

Once you have downloaded and installed Anaconda_ or Miniconda_
open a terminal (or ``cmd.exe`` on Windows).

With  Miniconda_, you can create 
a specific environment for AnuGA, by running::

    conda create -n anuga_env python=2.7
    source activate anuga_env
    
    
With either Anaconda_ or Miniconda_ you can now install the dependencies by running::

    conda install pip nose numpy scipy matplotlib netcdf4
    conda install -c pingucarsti gdal 
    
and setup GDAL_DATA environment variable:

    export GDAL_DATA=`gdal-config --datadir` 
    
(You should add this command to your .bashrc file.)    


Windows Dependency Installation
+++++++++++++++++++++++++++++++

We have successfully install AnuGA on windows using Gohlke Binaries and using Miniconda. 
At present we recommend using the Gohlke Binaries. 

Follow the instructions 
`install ANUGA on Windows using the Gohlke Binaries
<https://github.com/GeoscienceAustralia/anuga_core/wiki/Install-ANUGA-on-Windows-using-Gohlke-Binaries>`_

Alternatively if you want ot use Miniconda, follow the instructions 
`install ANUGA on Windows using Miniconda
<https://github.com/GeoscienceAustralia/anuga_core/wiki/Install-ANUGA-on-Windows-using-Miniconda>`_



GCC dependency for Windows users
++++++++++++++++++++++++++++++++++++

Unfortunately, the ``gcc`` compiler MinGW_ included in Anaconda or 
installable via Miniconda_ doesn't have OpenMP_ support. This is required to compile
some extension modules in AnuGA (those that have multi-threaded parallel code).

We suggest that you download and install the version of MinGW_ provided by TDM-GCC_
**after** you've installed Anaconda and **before** you install AnuGA.
Don't forget to mark the ``openmp`` and ``gfortran`` options in the "Choose Components" part of
the installation. See this `excellent documentation for Windows users`_
(they even have screenshots!). The same applies if you are using Miniconda_.


Installing AnuGA
----------------

After you've installed the dependencies you can proceed to install AnuGA
using pip_.
Open a terminal (or ``cmd.exe`` on Windows) and run::

    pip install anuga


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
`passing all tests <https://travis-ci.org/GeoscienceAustralia/anuga_core/branches>`__,
so it should be safe to use.

First, you'll need to `install git`_.
Then, open a terminal and run::

    git clone git://github.com/GeoscienceAustralia/anuga_core.git 

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
please write to the `mailing list`_ or `create an issue on the github site 
<https://github.com/GeoscienceAustralia/anuga_core/issues>`__.

To make it easier for us to debug you problem, please include the following
information:

* Operating system
* Python distribution (Anaconda_, PythonXY_, `ETS/Canopy`_, own install)
* Python version (2.6, 2.7 etc)
* The script you ran (and gave you an error/weird result)
* The error message (the part that says ``Traceback: ...``) or result (figure,
  numbers, etc)
    
.. _AnuGA http://anuga.anu.edu.au/ 
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

