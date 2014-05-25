"""Run parallel shallow water domain.

   run using command like:

   mpirun -np m python run_sequential_dist_all.py

   where m is the number of processors to be used.
   
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


#from anuga_parallel.sequential_distribute import sequential_distribute

from anuga_parallel import distribute, myid, numprocs, finalize, barrier


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
if myid == 0 :
    domain = create_domain_from_file(mesh_filename)
    domain.set_quantity('stage', Set_Stage(x0, x1, 2.0))
    #domain.set_datadir('.')
    domain.set_name('merimbula_new')
    domain.set_store(True)
else:
    domain = None


#--------------------------------------------------------------------------
# Distribute sequential domain on processor 0 to other processors
#--------------------------------------------------------------------------
if verbose: print 'DISTRIBUTING DOMAIN'
domain = distribute(domain)


#--------------------------------------------------------------------------
# On all processors, setup evolve parameters for domains on all processors
# (all called "domain"
#--------------------------------------------------------------------------

domain.set_flow_algorithm('DE0')

#------------------------------------------------------------------------------
# Setup boundary conditions
# This must currently happen *after* domain has been distributed
#------------------------------------------------------------------------------
Br = Reflective_boundary(domain)      # Solid reflective wall

domain.set_boundary({'outflow' :Br, 'inflow' :Br, 'inner' :Br, 'exterior' :Br, 'open' :Br})

#------------------------------------------------------------------------------
# Evolution
#------------------------------------------------------------------------------
if myid == 0 and verbose: print 'EVOLVE'

barrier()
t0 = time.time()
for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
    if myid == 0:
        domain.write_time()


barrier()

for p in range(numprocs):
    if myid == p:
        print 'Processor %g ' %myid
        print 'That took %.2f seconds' %(time.time()-t0)
        print 'Communication time %.2f seconds'%domain.communication_time
        print 'Reduction Communication time %.2f seconds'%domain.communication_reduce_time
        print 'Broadcast time %.2f seconds'%domain.communication_broadcast_time
    else:
        pass

    barrier()


#--------------------------------------------------
# Merge the individual sww files into one file
#--------------------------------------------------
domain.sww_merge(delete_old=False)

finalize()











