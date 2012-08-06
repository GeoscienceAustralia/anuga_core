# =============================================================================
# install.py - Installation of pypar (Parallel Python using MPI)
# Copyright (C) 2001 Ole M. Nielsen
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License (http://www.gnu.org/copyleft/gpl.html)
#    for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
#
#
# Contact address: Ole.Nielsen@anu.edu.au
#
# Version 1.0 October 2001
# =============================================================================

"""Module install.py - Compile C extension for pypar
"""

import os, sys      
    

# Attempt to compile mpi.so 
#
try:
  import compile
  compile.compile('mpi.c', 'mpicc', verbose = 1)
except:
  raise "Could compile C extension mpi.c - please try manually"       


if sys.platform in ['osf1V5', 'sunos5']:  #Compaq AlphaServer or Sun
  sys.exit()
  # Alpha Server or Sun cannot import MPI on sequential processes
  # OK with LAM

#Check if MPI module can be initialised
#
# Attempt to import and initialise mpi.so

error = os.system('python -c "import mpi" > /dev/null') 
if error:
  raise "MPI could not be initialised."
else:  
  import mpi
  print "MPI initialised OK"    

