#!/usr/bin/env python
#########################################################
#
#  Main file for parallel mesh testing. Runs an advection
# flow simulation using a rectangular mesh
#
#  This is a modification of the run_parallel_advection.py
# file
#
#
# *) The (new) files that have been added to manage the
# grid partitioning are
#    +) pmesh_divide_metis.py: subdivide a pmesh
#    +) build_submesh.py: build the submeshes on the host
# processor.
#    +) build_local.py: build the GA mesh datastructure
# on each processor.
#    +) build_commun.py: handle the communication between
# the host and processors
#
#  Authors: Linda Stals, Steve Roberts and Matthew Hardy,
# June 2005
#
#
#
#########################################################

import pypar    # The Python-MPI interface
import time

from sys import path

print path

# Numeric arrays
from Numeric import array, zeros, Float

from print_stats import print_test_stats, build_full_flag

from anuga.abstract_2d_finite_volumes.pmesh2domain\
     import pmesh_to_domain_instance
from anuga.advection.advection import Domain as Advection_Domain
from parallel_advection import Parallel_Domain

from anuga.abstract_2d_finite_volumes.generic_boundary_conditions\
     import Transmissive_boundary

# mesh partition routines

from pmesh_divide  import pmesh_divide_metis
from build_submesh import build_submesh
from build_local   import build_local_mesh
from build_commun  import send_submesh, rec_submesh, extract_submesh



class Set_Stage:
    """Set an initial condition with constant water height, for x<x0
    """

    def __init__(self, x0=0.25, x1=0.5, h=1.0):
        self.x0 = x0
        self.x1 = x1
        self.h  = h

    def __call__(self, x, y):
        return self.h*((x>self.x0)&(x<self.x1))

# read in the processor information

numprocs = pypar.size()
myid = pypar.rank()
processor_name = pypar.Get_processor_name()


#-------
# Domain
rect = zeros( 4, Float) # Buffer for results

if myid == 0:

    # read in the test files

    filename = 'test-100.tsh'
#    filename = 'merimbula_10785.tsh'
    nx = numprocs
    ny = 1
    if nx*ny != numprocs:
        print "WARNING: number of subboxes is not equal to the number of proc"

    domain_full = pmesh_to_domain_instance(filename, Advection_Domain)
    domain_full.set_quantity('stage', Set_Stage(200.0,300.0,1.0))
#    domain_full.set_quantity('stage', Set_Stage(756000.0,756500.0,4.0))

    nodes, triangles, boundary, triangles_per_proc, quantities  =\
            pmesh_divide_metis(domain_full, numprocs)
    
    # subdivide the mesh
    rect = array(domain_full.xy_extent, Float)
#    rect = array(rect, Float)

    submesh = build_submesh(nodes, triangles, boundary, quantities, \
                            triangles_per_proc)

    # send the mesh partition to the appropriate processor

    for p in range(1, numprocs):
      send_submesh(submesh, triangles_per_proc, p)

    # Build the local mesh for processor 0
    
    points, vertices, boundary, quantities, ghost_recv_dict, full_send_dict = \
              extract_submesh(submesh, triangles_per_proc)

# read in the mesh partition that belongs to this
# processor (note that the information is in the
# correct form for the GA data structure

else:
    [points, vertices, boundary, quantities, ghost_recv_dict, full_send_dict] = \
             rec_submesh(0)

pypar.broadcast(rect,0)

domain = Parallel_Domain(points, vertices, boundary,
                                   full_send_dict  = full_send_dict,
                                   ghost_recv_dict = ghost_recv_dict,
                                   velocity = [0.1,0.0])

# Make a notes of which triangles are full and which are ghost

tri_full_flag = build_full_flag(domain, ghost_recv_dict)

#domain.initialise_visualiser(rect=rect)

#Boundaries

T = Transmissive_boundary(domain)
#R = Reflective_boundary(domain)
domain.set_boundary( {'outflow': T, 'inflow': T, 'inner':T, 'exterior': T, 'open':T, 'ghost':None} )



domain.set_quantity('stage', quantities['stage'])

#---------
# Evolution
t0 = time.time()


try:
    domain.initialise_visualiser(rect=rect)
    #domain.visualiser.coloring['stage'] = True
    #domain.visualiser.scale_z['stage'] = 0.2
except:
    print 'No visualiser'


from norms import linf_norm

yieldstep = 1
finaltime = 200

#yieldstep = 1000
#finaltime = 50000

for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
    if myid == 0:
        domain.write_time()
    print_test_stats(domain, tri_full_flag)

if myid == 0:
    print 'That took %.2f seconds' %(time.time()-t0)
