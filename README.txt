

=====
AnuGA
=====

.. contents::

What is AnuGA?
--------------

AnuGA (pronounced "an uu ga") is open-source software for the simulation of
the shallow water equation, in particular it can be used to model tsunamis
and floods.

Developed at the Risk Assessment Methods Project at Geoscience
Australia and Mathematical Sciences Institute at the Australian
National University.


Copyright 2004, 2005, 2006, 2015 
Ole Nielsen, Stephen Roberts, Duncan Gray, Jane Sexton


Installation
------------

For installation instructions, see ``INSTALL.txt``.


Documentation
-------------

AnuGA documentation is available in the user_manual folder or online

    http://anuga.anu.edu.au

How to generate the HTML documentation, see ``doc/README.txt``.


Web sites
---------

The user's site is:

    http://anuga.anu.edu.au/


Mailing Lists
-------------

Please see anuga user mailing list here:

    anuga-user@lists.sourceforge.net


Latest source code
------------------

The latest development version of AnuGA's sources are always available at:

    https://github.com/stoiver/anuga_core

They can be downloaded as a zip file or using the Git client.


Bug reports
-----------

To search for bugs or report them, please use the AnuGA's Bug Tracker at:

    https://github.com/stoiver/anuga_core/issues


Developer information
---------------------

If you would like to take part in AnuGA development, take a look
at ``HACKING.rst.txt``.


License information
-------------------

See the file ``LICENSE.txt`` for information on the history of this
software, terms & conditions for usage, and a DISCLAIMER OF ALL
WARRANTIES.





This is the root-Subversion directory for the ANUGA Inundation 
Software Development project conducted by






PREREQUISITES

python 2.3 or later
python-numeric
python-dev (Interestingly, numeric seems to be installed as well)
python-scientific (or maybe only python-netcdf)
 (Remember to put netcdf.dll somewhere so the PATH - e.g. in C:\bin)

A C compiler such as gcc (from GNU in case of Linux and MinGW in case of Windows)

RECOMMENDED
psyco
visual python
matplotlib (pylab) for quality 2d plotting
weave?


INSTALLATION

In directory anuga:

  python test_all.py
  

PATH


PATH=C:\Python24;C:\MinGW\bin;C:\Program Files\anuga_viewer

It is assumed that the root directory for all anugu modules is on the pythonpath, e.g.
PYTHONPATH=V:\home\onielsen\anuga_core\source
Alternatively this directory should be copied to the python site_packages directory 
  

  
  
LINKS  
    Python: http://www.python.org
    Visual Python: http://vpython.org/
    Python NetCDF Interface: http://www-md.fsl.noaa.gov/eft/developer/netCDFPythonInterface.html    
    Psyco: http://psyco.sourceforge.net/index.html


LICENSE
    see LICENSE.txt
    

