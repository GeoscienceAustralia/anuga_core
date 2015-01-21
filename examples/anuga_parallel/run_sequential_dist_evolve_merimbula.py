
"""Run parallel shallow water domain.

   run using command like:

   mpirun -np 4  python run_sequential_dist_evolve_merimbula.py

   Need to have run run_sequential_dist_distribute_merimbula.py

   first to produce the required pickle files encoding the partitioning
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


from anuga_parallel import distribute, myid, numprocs, finalize, barrier

from anuga_parallel.sequential_distribute import sequential_distribute_load

#--------------------------------------------------------------------------
# Setup parameters
#--------------------------------------------------------------------------
yieldstep = 50
finaltime = 1500
verbose = True



#-------------------------------------------------------------------------------
# Read in cPickled distribute mesh data
#-------------------------------------------------------------------------------
domain = sequential_distribute_load(filename='merimbula_new', verbose = True)

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




finalize()




