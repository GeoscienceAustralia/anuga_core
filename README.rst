
.. image:: https://app.travis-ci.com/anuga-community/anuga_core.svg?branch=main
    :target: https://app.travis-ci.com/anuga-community/anuga_core
    :alt: travis ci status
   
.. image:: https://ci.appveyor.com/api/projects/status/x5airjv7eq2u805w/branch/main?svg=true
    :target: https://ci.appveyor.com/project/stoiver/anuga-core-nwgr0
    :alt: appveyor status

.. image:: https://img.shields.io/pypi/v/anuga.svg
    :target: https://pypi.python.org/pypi/anuga/
    :alt: Latest PyPi Version

.. image:: https://img.shields.io/pypi/dm/anuga.svg
    :target: https://pypistats.org/packages/anuga
    :alt: PyPi download statistics

.. image:: https://img.shields.io/conda/vn/conda-forge/anuga.svg
    :target: https://anaconda.org/conda-forge/anuga
    :alt: Latest Conda Version
 
.. image:: https://img.shields.io/conda/dn/conda-forge/anuga.svg
    :target: https://anaconda.org/conda-forge/anuga
    :alt: Conda Forge download statistics

.. image:: https://readthedocs.org/projects/anuga/badge/?version=latest
    :target: https://anuga.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status


=====
ANUGA
=====

.. contents::

What is ANUGA?
--------------

ANUGA (pronounced "AHnooGAH") is open-source software for the simulation of
the shallow water equation, in particular it can be used to model tsunamis
and floods.

ANUGA is a python 3 package with some C and Cython extensions (and an optional
fortran extension). 

ANUGA was created in a collaboration by Geoscience Australia and Mathematical Sciences Institute at the
Australian National University, though now it is further developed and maintained by a community of volunteers.

Copyright Commonwealth of Australia (Geoscience Australia) and Australian National University 2004-Now


Where can I get ANUGA?
----------------------

ANUGA is available from either https://github.com/anuga-community/anuga_core or https://github.com/GeoscienceAustralia/anuga_core.

The Python 3.x version is the default and available in the main branches of both repositories. 

If you rely on the old Python 2.7 version, it is available in branches named anuga_py2.

The python 3 version of ANUGA will continue to be developed and the most up-to-date 
versions will be available from the `anuga-community <https://github.com/anuga-community/anuga_core>`_ repository. 



Installation
------------

If you use `conda` to provide your `python` environment, then you can install `anuga` from the conda-forge channel
as follows. First setup the `conda-forge` channel:

.. code-block::

    conda config --add channels conda-forge
    conda config --set channel_priority strict


Once the conda-forge channel has been enabled, anuga can be installed with conda:

.. code-block::

    conda install anuga


For more installation instructions, see https://anuga.readthedocs.io/en/latest/installation.html


Documentation and Help
----------------------


ANUGA documentation is available via "read the docs" at 

    https://anuga.readthedocs.io 

Also you can obtain help via the old
`user_manual <https://github.com/anuga-community/anuga_core/raw/main/doc/anuga_user_manual.pdf>`_

Also helpful information is available online at

    http://anuga.anu.edu.au

A collection of online jupyter notebooks which can run under google's colab environment can be found at:

    https://github.com/anuga-community/anuga-clinic

Mailing Lists
-------------

You can subscribe to our mailing via:

    https://lists.sourceforge.net/lists/listinfo/anuga-user

and send questions using the address

    anuga-user@lists.sourceforge.net

You can also submit issues to:

    https://github.com/anuga-community/anuga_core/issues


Web sites
---------

* http://anuga.anu.edu.au: Collection of information, talks and papers about ANUGA.
* https://en.wikipedia.org/wiki/ANUGA_Hydro: The Wikipedia site for ANUGA. 
* https://github.com/anuga-community/anuga_core: The active GitHub repository for ANUGA.
* https://github.com/GeoscienceAustralia/anuga_core: Mirror GitHub repository for ANUGA. 
* https://github.com/anuga-community/anuga-viewer: Viewer for animating the ANUGA sww output files.  



Latest source code
------------------

The latest development version of ANUGA's sources are is available at:

    https://github.com/anuga-community/anuga_core

They can be downloaded as a zip file or using the Git client as follows

.. code-block::

    git clone https://github.com/anuga-community/anuga_core #(for read only)
    git clone git@github.com:anuga-community/anuga_core.git #(to contribute)

For the latter option see e.g. https://help.github.com/articles/adding-a-new-ssh-key-to-your-github-account/ for instructions on how to authenticate using ssh keys.

Bug reports
-----------

To search for bugs or report them, please use the ANUGA's Bug Tracker at:

    https://github.com/anuga-community/anuga_core/issues


Developer information
---------------------

If you would like to take part in ANUGA development, take a look
at `CONTRIBUTING.rst <https://github.com/anuga-community/anuga_core/blob/main/CONTRIBUTING.rst>`_.


License information
-------------------

See the file `LICENSE.txt <https://github.com/anuga-community/anuga_core/blob/main/LICENCE.txt>`_
for information on the history of this software, terms & conditions for usage,
and a DISCLAIMER OF ALL WARRANTIES.

Contacts
--------

At the Australian National University:
**Stephen Roberts**
*Lead Developer*
<stephen.roberts@anu.edu.au>

At Geoscience Australia:
**Gareth Davies**
*Developer*
<gareth.davies@ga.gov.au>

ANUGA Community:
**Ole Nielsen**
*Architect and Developer*
<ole.moller.nielsen@gmail.com>
