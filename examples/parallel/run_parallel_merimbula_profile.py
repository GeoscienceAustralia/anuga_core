
"""Run parallel shallow water domain.

   run using command like:

   mpiexec -np m python run_parallel_sw_merimbula.py

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

from anuga import rectangular_cross
from anuga import create_domain_from_file


from anuga import distribute, myid, numprocs, finalize, barrier


#--------------------------------------------------------------------------
# Setup parameters
#--------------------------------------------------------------------------

mesh_filename = "data/merimbula_10785_1.tsh" ; x0 = 756000.0 ; x1 = 756500.0
#mesh_filename = "data/merimbula_43200.tsh"   ; x0 = 756000.0 ; x1 = 756500.0
#mesh_filename = "data/test-100.tsh" ; x0 = 0.25 ; x1 = 0.5

finaltime = 500
yieldstep = 50
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
        return self.h*((x>self.x0)&(x<self.x1)) + 1.0

#--------------------------------------------------------------------------
# Setup Domain only on processor 0
#--------------------------------------------------------------------------
if myid == 0:
    domain = create_domain_from_file(mesh_filename)
    domain.set_quantity('stage', Set_Stage(x0, x1, 2.0))
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

    barrier()


domain.set_flow_algorithm(2.0)

if myid == 0:
    domain.print_algorithm_parameters()

barrier()

domain.set_name('meribula')

#------------------------------------------------------------------------------
# Setup boundary conditions
# This must currently happen *after* domain has been distributed
#------------------------------------------------------------------------------
Br = Reflective_boundary(domain)      # Solid reflective wall

domain.set_boundary({'exterior' :Br, 'open' :Br})


barrier()
for p in range(numprocs):
    if myid == p:
        print (domain.boundary_statistics())

    barrier()

#------------------------------------------------------------------------------
# Evolution
#------------------------------------------------------------------------------
if myid == 0 and verbose: print ('EVOLVE')

t0 = time.time()


s = """
for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
    if myid == 0:
        domain.write_time()
"""


# Profiling
import cProfile
prof_file = 'evolve-prof'+ str(numprocs) + '_' + str(myid) + '.dat'
cProfile.run(s,prof_file)


barrier()

if myid == 0:
    import pstats
    p = pstats.Stats(prof_file)
    #p.sort_stats('cumulative').print_stats(25)


    p.sort_stats('time').print_stats(25)

barrier()


if myid == 1:
    import pstats
    p = pstats.Stats(prof_file)
    #p.sort_stats('cumulative').print_stats(25)


    p.sort_stats('time').print_stats(25)

    #p.print_stats()

# Run evolve loop
#result = profiler.runctx(s, globals(), locals())

#result.dump_stats("profile." + str(numprocs) + "." + str(myid) + ".dat")
#profiler.close()

barrier()

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




