
"""Run parallel shallow water domain.

   run using command like:

   mpiexec -np m python run_parallel_sw_rectangular_cross.py

   where m is the number of processors to be used.
   
   Will produce sww files with names domain_Pn_m.sww where n is number of processors and
   m in [0, n-1] refers to specific processor that owned this part of the partitioned mesh.
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

from anuga import rectangular_cross_domain
from anuga import create_domain_from_file


from anuga import distribute, myid, numprocs, finalize, barrier


#--------------------------------------------------------------------------
# Setup parameters
#--------------------------------------------------------------------------

#mesh_filename = "merimbula_10785_1.tsh" ; x0 = 756000.0 ; x1 = 756500.0
#mesh_filename = "merimbula_43200.tsh"   ; x0 = 756000.0 ; x1 = 756500.0
#mesh_filename = "test-100.tsh" ; x0 = 0.25 ; x1 = 0.5


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
        return self.h*((x>self.x0)&(x<self.x1))

#--------------------------------------------------------------------------
# Setup Domain only on processor 0
#--------------------------------------------------------------------------
if myid == 0:
    length = 2.0
    width = 2.0
    dx = dy = 0.01
    domain = rectangular_cross_domain(int(length/dx), int(width/dy),
                                              len1=length, len2=width)

    domain.set_quantity('elevation', -1.0)
    domain.set_quantity('stage', Set_Stage(h=2.0))
else:
    domain = None

#--------------------------------------------------------------------------
# Distribute sequential domain on processor 0 to other processors
#--------------------------------------------------------------------------

if myid == 0 and verbose: print ('DISTRIBUTING DOMAIN')
domain = distribute(domain)

#domain.smooth = False
barrier()
for p in range(numprocs):
    if myid == p:
        print ('Process ID %g' %myid)
        print ('Number of triangles %g ' %domain.get_number_of_triangles())
        sys.stdout.flush()
        
    barrier()


domain.set_flow_algorithm(2.0)

if myid == 0:
    domain.print_algorithm_parameters()
    sys.stdout.flush()
    
barrier()

domain.set_name('rectangular_cross')
domain.set_store(False)

#------------------------------------------------------------------------------
# Setup boundary conditions
# This must currently happen *after* domain has been distributed
#------------------------------------------------------------------------------
Br = Reflective_boundary(domain)      # Solid reflective wall

domain.set_boundary({'top' :Br, 'bottom' :Br, 'left' :Br, 'right' :Br, 'ghost' : None })

"""
barrier()
for p in range(numprocs):
    if myid == p:
        print domain.boundary_statistics()
        sys.stdout.flush()
        
    barrier()
"""


#domain.dump_triangulation()

#------------------------------------------------------------------------------
# Evolution
#------------------------------------------------------------------------------
if myid == 0 and verbose: print ('EVOLVE')

t0 = time.time()
finaltime = 0.25
yieldstep = 0.05

for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
    if myid == 0:
        domain.write_time()



## Profiling
#import cProfile
#prof_file = 'evolve-prof'+ str(numprocs) + '_' + str(myid) + '.dat'
#cProfile.run(s,prof_file)


barrier()

#for id in range(numprocs):
#    if myid == id:
#        import pstats
#        p = pstats.Stats(prof_file)
#        #p.sort_stats('cumulative').print_stats(25)
#        p.sort_stats('time').print_stats(25)
#        sys.stdout.flush
#
#    barrier()
#
#barrier()

for p in range(numprocs):
    if myid == p:
        print ('Process ID %g' %myid)
        print ('Number of processors %g ' %numprocs)
        print ('That took %.2f seconds' %(time.time()-t0))
        print ('Communication time %.2f seconds'%domain.communication_time)
        print ('Reduction Communication time %.2f seconds'%domain.communication_reduce_time)
        print ('Broadcast time %.2f seconds'%domain.communication_broadcast_time)

    barrier()


domain.sww_merge(delete_old=True)

finalize()




