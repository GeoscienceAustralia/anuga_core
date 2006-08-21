This is the root-Subversion directory for the ANUGA Inundation 
Software Development project conducted by

Stephen Roberts (ANU)
Ole Nielsen (GA)
Duncan Gray (GA)
Christopher Zoppou (GA)


This is the top directory for the inundation project source code
developed at the Risk Assessment Methods Project at Geoscience
Australia and Mathematical Sciences Institute at the Australian
National University.


Copyright 2004, 2005, 2006  
Ole Nielsen, Stephen Roberts, Duncan Gray, Jane Sexton, Christopher Zoppou


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
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License (http://www.gnu.org/copyleft/gpl.html)
    for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307

