# =============================================================================
# pypar.py - Parallel Python using MPI
# Copyright (C) 2001, 2002, 2003 Ole M. Nielsen 
#              (Center for Mathematics and its Applications ANU and APAC)
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
# Contact address: Ole.Moller.Nielsen@gmail.com
#
# Version: See pypar.__version__
# =============================================================================

"""The python module pypar.py and the C-extension mpi.c 
implements scalable parallelism on distributed and shared 
memory architectures using an important subset of 
the Message Passing Interface (MPI) standard.
"""

from pypar import *
from __metadata__ import __version__, __date__, __author__

#Add path of package to PYTHONPATH to allow C-extension to be loaded
#import sys
#sys.path += __path__




