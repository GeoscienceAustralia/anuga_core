This is the root-Subversion directory for the ANUGA Inundation 
Software Development project conducted by

Stephen Roberts (ANU)
Ole Nielsen (GA)
Duncan Gray (GA)
Jane Sexton (GA)
Nick Bartzis (GA)
Jack Kelly (ANU)
John Jakeman (ANU)


This is the top directory for the inundation project source code
developed at the Risk Assessment Methods Project at Geoscience
Australia and Mathematical Sciences Institute at the Australian
National University.


Copyright 2004, 2005, 2006  
Ole Nielsen, Stephen Roberts, Duncan Gray, Jane Sexton


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

In directory pyvolution:

  python test_all.py
  

PATH


PATH=C:\Python24;C:\MinGW\bin;C:\Program Files\swollen

It is assumed that the root directory for all anuga/inundation modules is on the pythonpath, e.g.
PYTHONPATH=V:\1\cit\risk_assessment_methods_project\inundation\sandpits\onielsen\anuga\inundation
  

  
  
LINKS  
    Python: http://www.python.org
    Visual Python: http://vpython.org/
    Python NetCDF Interface: http://www-md.fsl.noaa.gov/eft/developer/netCDFPythonInterface.html    
    Psyco: http://psyco.sourceforge.net/index.html


LICENSE
    see LICENSE.txt
    

