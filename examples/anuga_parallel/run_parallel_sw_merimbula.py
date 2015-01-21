
"""Run parallel shallow water domain.

   run using command like:

   mpirun -np m python run_parallel_sw_merimbula.py

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
	
from anuga import Domain
from anuga import Reflective_boundary
from anuga import Dirichlet_boundary
from anuga import Time_boundary
from anuga import Transmissive_boundary
from anuga import Transmissive_n_momentum_zero_t_momentum_set_stage_boundary

from anuga import rectangular_cross
from anuga import create_domain_from_file


from anuga_parallel import distribute, myid, numprocs, finalize, barrier


#--------------------------------------------------------------------------
# Setup parameters
#--------------------------------------------------------------------------

mesh_filename = "data/merimbula_10785_1.tsh" ; x0 = 756000.0 ; x1 = 756500.0; yieldstep = 10; finaltime = 100
#mesh_filename = "data/merimbula_17156.tsh"   ; x0 = 756000.0 ; x1 = 756500.0; yieldstep = 50; finaltime = 500
#mesh_filename = "data/merimbula_43200_1.tsh"   ; x0 = 756000.0 ; x1 = 756500.0; yieldstep = 50; finaltime = 500
#mesh_filename = "data/test-100.tsh" ; x0 = 200.0 ; x1 = 300.0; yieldstep = 1; finaltime = 10
#mesh_filename = "data/test-20.tsh" ; x0 = 250.0 ; x1 = 350.0; yieldstep = 1; finaltime = 50

verbose = False

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
    #domain.set_quantity('stage', 0.0)
    #domain.set_datadir('.')
    domain.set_name()
    domain.set_store(True)
    #domain.set_store_vertices_smoothly(True)
    domain.set_store_vertices_smoothly(False)
    #domain.set_quantity('elevation', Set_Elevation(500.0))

    #print domain.statistics()
    #print domain.get_extent()
    #print domain.get_extent(absolute=True)
    #print domain.geo_reference
else:
    domain = None

#--------------------------------------------------------------------------
# Distribute sequential domain on processor 0 to other processors
#--------------------------------------------------------------------------

if myid == 0 and verbose: print 'DISTRIBUTING DOMAIN'
domain = distribute(domain, verbose=verbose)

#--------------------------------------------------------------------------
# On all processors, setup evolve parameters for domains on all processors
# (all called "domain"
#--------------------------------------------------------------------------

domain.set_flow_algorithm('2_0')
domain.set_store_centroids()

#for p in range(numprocs):
#    if myid == p:
#        print 'P%d'%p
#        print domain.get_extent()
#        print domain.get_extent(absolute=True)
#        print domain.geo_reference
#        print domain.s2p_map
#        print domain.p2s_map
#        print domain.tri_l2g
#        print domain.node_l2g
#    else:
#        pass
#
#    barrier()


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

domain.set_boundary({'outflow' :Br, 'inflow' :Br, 'inner' :Br, 'exterior' :Br, 'open' :Bts})


#bdy_ids = domain.tag_boundary_cells['open']
#vol_ids = domain.boundary_cells[bdy_ids]
#edge_ids = domain.boundary_edges[bdy_ids]

#print domain.mesh.edge_midpoint_coordinates[3*vol_ids+edge_ids]



#------------------------------------------------------------------------------
# Evolution
#------------------------------------------------------------------------------
if myid == 0 and verbose: print 'EVOLVE'

barrier()
t0 = time.time()
for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
    #print 'P:'+ str(myid) + ' '+ domain.timestepping_statistics()
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
domain.sww_merge(delete_old=True)


domain.dump_triangulation()

finalize()




