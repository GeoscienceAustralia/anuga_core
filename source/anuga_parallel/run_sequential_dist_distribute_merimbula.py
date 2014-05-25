
"""Run parallel shallow water domain.

   run using command like:

   python run_sequential_dist_distribute_merimbula.py

   will produce 4 pickled files corresponding partitioned domains
   
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

import os
import sys
import time
import numpy as num

#------------------------
# ANUGA Modules
#------------------------
	
from anuga import Domain
from anuga import Reflective_boundary
from anuga import Dirichlet_boundary
from anuga import Time_boundary
from anuga import Transmissive_boundary

from anuga import rectangular_cross
from anuga import create_domain_from_file


from anuga_parallel.sequential_distribute import sequential_distribute_dump


#--------------------------------------------------------------------------
# Setup parameters
#--------------------------------------------------------------------------

#mesh_filename = "merimbula_10785_1.tsh" ; x0 = 756000.0 ; x1 = 756500.0
mesh_filename = "merimbula_17156.tsh"   ; x0 = 756000.0 ; x1 = 756500.0
#mesh_filename = "merimbula_43200_1.tsh"   ; x0 = 756000.0 ; x1 = 756500.0
#mesh_filename = "test-100.tsh" ; x0 = 0.25 ; x1 = 0.5
#mesh_filename = "test-20.tsh" ; x0 = 250.0 ; x1 = 350.0
yieldstep = 50
finaltime = 1500
verbose = True

#--------------------------------------------------------------------------
# Setup procedures
#--------------------------------------------------------------------------
class Set_Stage:
    """Set an initial condition with constant water height, for x0<x<x1
    """

    def __init__(self, x0=0.25, x1=0.5, h=1.0):
        self.x0 = x0
        self.x1 = x1
        self.h  = h

    def __call__(self, x, y):
        return self.h*((x>self.x0)&(x<self.x1))+1.0


class Set_Elevation:
    """Set an elevation
    """

    def __init__(self, h=1.0):
        self.x0 = x0
        self.x1 = x1
        self.h  = h

    def __call__(self, x, y):
        return x/self.h
    

#--------------------------------------------------------------------------
# Setup Sequential Domain
#--------------------------------------------------------------------------
domain = create_domain_from_file(mesh_filename)
domain.set_quantity('stage', Set_Stage(x0, x1, 2.0))
#domain.set_datadir('.')
domain.set_name('merimbula_new')
domain.set_store(True)


#--------------------------------------------------------------------------
# Distribute sequential domain on processor 0 to other processors
#--------------------------------------------------------------------------

if verbose: print 'DISTRIBUTING DOMAIN'
sequential_distribute_dump(domain, 4, verbose=True)





