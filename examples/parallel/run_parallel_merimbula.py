
"""Run parallel shallow water domain.

   run using command like:

   mpiexec -np m python run_parallel_merimbula.py

   where m is the number of processors to be used.
   
   Will produce sww files with names domain_Pn_m.sww where m is number of processors and
   n in [0, m-1] refers to specific processor that owned this part of the partitioned mesh.
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
import anuga
	
from anuga import Domain
from anuga import Reflective_boundary
from anuga import Dirichlet_boundary
from anuga import Time_boundary
from anuga import Transmissive_boundary
from anuga import Transmissive_n_momentum_zero_t_momentum_set_stage_boundary

from anuga import rectangular_cross
from anuga import create_domain_from_file


from anuga import distribute, myid, numprocs, finalize, barrier

from anuga.utilities.system_tools import get_pathname_from_package

#--------------------------------------------------------------------------
# Setup parameters
#--------------------------------------------------------------------------

DATA_DIR = 'data'

mesh_filename = anuga.join(DATA_DIR,"merimbula_10785_1.tsh") ; x0 = 756000.0 ; x1 = 756500.0; yieldstep = 10; finaltime = 100

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
# Setup Domain only on processor 0
#--------------------------------------------------------------------------
if myid == 0:
    domain = create_domain_from_file(mesh_filename)
    domain.set_quantity('stage', Set_Stage(x0, x1, 1.0))
    domain.set_store_vertices_smoothly(False)

else:
    domain = None

#--------------------------------------------------------------------------
# Distribute sequential domain on processor 0 to other processors
#--------------------------------------------------------------------------

if myid == 0 and verbose: print ('DISTRIBUTING DOMAIN')
domain = distribute(domain, verbose=verbose)

#--------------------------------------------------------------------------
# On all processors, setup evolve parameters for domains on all processors
# (all called "domain"
#-------------------------------------------------------------------------

domain.set_quantities_to_be_stored({'elevation':1,
                                    'friction':1,
                                    'stage':2,
                                    'xmomentum':2,
                                    'ymomentum':2})
                                 
#------------------------------------------------------------------------------
# Setup boundary conditions
# This must currently happen *after* domain has been distributed
#------------------------------------------------------------------------------
Br = Reflective_boundary(domain)      # Solid reflective wall

from math import sin
Bts = Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(domain, lambda t: 10*sin(t/60))

#domain.set_boundary({'outflow' :Br, 'inflow' :Br, 'inner' :Br, 'exterior' :Br, 'open' :Bts})
domain.set_boundary({'exterior' :Br, 'open' :Bts})

#------------------------------------------------------------------------------
# Evolution
#------------------------------------------------------------------------------
if myid == 0 and verbose: print ('EVOLVE')

barrier()
t0 = time.time()
for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
    #print 'P:'+ str(myid) + ' '+ domain.timestepping_statistics()
    if myid == 0:
        domain.print_timestepping_statistics()


barrier()

for p in range(numprocs):
    if myid == p:
        print (50*'=')
        print ('Processor %g ' %myid)
        print ('That took %.2f seconds' %(time.time()-t0))
        print ('Communication time %.2f seconds'%domain.communication_time)
        print ('Reduction Communication time %.2f seconds'%domain.communication_reduce_time)
        print ('Broadcast time %.2f seconds'%domain.communication_broadcast_time)
    else:
        pass

    barrier()


#--------------------------------------------------
# Merge the individual sww files into one file
#--------------------------------------------------
domain.sww_merge(delete_old=True)
    
if myid==0:
    print(50*"=")
    print('Number of Tirangles:', domain.number_of_global_triangles)

#domain.dump_triangulation()

finalize()




