#!/usr/bin/env python
#########################################################
#
#  Main file for parallel mesh testing.
#
#  This is a modification of the run_parallel_advection.py
# file.
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
# *) Things still to do:
#    +) Overlap the communication and computation: The
# communication routines in build_commun.py should be
# interdispersed in the build_submesh.py and build_local.py
# files. This will overlap the communication and
# computation and will be far more efficient. This should
# be done after more testing and there more confidence in
# the subpartioning.
#    +) Much more testing especially with large numbers of
# processors.
#  Authors: Linda Stals, Steve Roberts and Matthew Hardy,
# June 2005
#
#
#
#########################################################
import sys
import pypar    # The Python-MPI interface
import time

# Numeric arrays
import numpy as num
#from numpy import array, zeros, float

# Print debugging information
from print_stats import print_test_stats, build_full_flag

# pmesh
from anuga.shallow_water import Domain
from parallel_shallow_water import Parallel_domain
from anuga.abstract_2d_finite_volumes.pmesh2domain\
     import pmesh_to_domain_instance

# Reuse previous mesh import
from anuga.caching import cache

# Mesh partition routines
from distribute_mesh  import pmesh_divide_metis
from distribute_mesh  import build_submesh
from distribute_mesh  import build_local_mesh
from distribute_mesh  import send_submesh, rec_submesh, extract_submesh


###############################
# Read in processor information
###############################

numprocs = pypar.size()
myid = pypar.rank()
processor_name = pypar.get_processor_name()

############################
# Set the initial conditions
############################

rect = num.zeros( 4, num.float) # Buffer for results

class Set_Stage:
    """Set an initial condition with constant water height, for x<x0
    """

    def __init__(self, x0=0.25, x1=0.5, h=1.0):
        self.x0 = x0
        self.x1 = x1
        self.h  = h

    def __call__(self, x, y):
        return self.h*((x>self.x0)&(x<self.x1))

#######################
# Partition the domain
#######################

if myid == 0:

    # Read in the test files

#    filename = 'test-100.tsh'
#    filename = 'merimbula_10785_1.tsh'
    filename = 'merimbula_43200.tsh'

    # Build the whole domain
    
    domain_full = pmesh_to_domain_instance(filename, Domain)

#    domain_full = cache(pmesh_to_domain_instance,
#               (filename, Domain),
#              dependencies = [filename])

    rect = num.array(domain_full.get_extent(), num.float)
    print rect

    # Initialise the wave

    #domain_full.set_quantity('stage', Set_Stage(200.0,300.0,1.0))
    domain_full.set_quantity('stage', Set_Stage(756000.0,756500.0,2.0))
#    domain_full.set_quantity('stage', Set_Stage(756000.0,756500.0,0.0))

    # Subdivide the domain

    # Note the different arguments compared with pmesh_divide,
    # pmesh_divide_steve etc.
    
    nodes, triangles, boundary, triangles_per_proc, quantities = \
         pmesh_divide_metis(domain_full, numprocs)

    print triangles_per_proc
    
    rect = num.array(domain_full.get_extent(), num.float)

    submesh = build_submesh(nodes, triangles, boundary,\
                            quantities, triangles_per_proc)

    # Send the mesh partition to the appropriate processor

    for p in range(1, numprocs):
      send_submesh(submesh, triangles_per_proc, p)

    # Build the local mesh for processor 0

    points, vertices, boundary, quantities, ghost_recv_dict, full_send_dict = \
             extract_submesh(submesh, triangles_per_proc)

# Read in the mesh partition that belongs to this
# processor (note that the information is in the
# correct form for the GA data structure

else:
    points, vertices, boundary, quantities, ghost_recv_dict, full_send_dict , \
            no_full_nodes, no_full_trigs = rec_submesh(0)


###########################################
# Start the computations on each subpartion
###########################################

#if myid == 0:
#    print 'ghost'
#    print ghost_recv_dict
#processor_name
#if myid == 0:
#    print 'full'
#    print full_send_dict

# The visualiser needs to know the size of the whole domain

pypar.broadcast(rect,0)

domain = Parallel_domain(points, vertices, boundary,
                         full_send_dict  = full_send_dict,
                         ghost_recv_dict = ghost_recv_dict)

# Make a note of which triangles are full and which are ghost

tri_full_flag = build_full_flag(domain, ghost_recv_dict)

try:
    #domain.initialise_visualiser(rect=rect)
    #domain.visualiser.coloring['stage'] = True
    #domain.visualiser.scale_z['stage'] = 0.2
    #domain.visualiser.scale_z['elevation'] = 0.05
    pass
except:
    print 'No visualiser'



domain.default_order = 1

#Boundaries
from anuga.interface import Transmissive_boundary, Reflective_boundary

T = Transmissive_boundary(domain)
R = Reflective_boundary(domain)
domain.set_boundary( {'outflow': R, 'inflow': R, 'inner':R, 'exterior': R, 'open':R, 'ghost':None} )


domain.set_quantity('stage', quantities['stage'])
domain.set_quantity('elevation', quantities['elevation'])

domain.store = False

#---------
# Evolution
t0 = time.time()

print 'Processor %d on %s: No of elements %d'%(domain.processor,processor_name,domain.number_of_elements)
yieldstep = 50.0
finaltime = 500.0

#yieldstep = 1000
#finaltime = 40000

#yieldstep = 1
#finaltime = 1
#processor_name
#for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
#    if myid == 0:
#        domain.write_time()
        #print 'Processor %d, Integral of stage %d'%\
        #       (domain.processor,domain.quantities['stage'].get_integral())
        #    print_test_stats(domain, tri_full_flag)


# Profiling
#import profile
#profiler = profile.Profile()
#result.dump_stats("profile." + str(numprocs) + "." + str(myid) + ".dat")

## #New hotshot profiling
## import hotshot
## profiler = hotshot.Profile("hotshot." + str(numprocs) + "." + str(myid) + ".prof")
## s = '''for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
##   if myid == 0:
##     domain.write_time()
##   print_test_stats(domain, tri_full_flag)

## '''
## result = profiler.runctx(s, globals(), locals())
## profiler.close()

#from vtk_realtime_visualiser import Visualiser
#V = Visualiser(domain,default_scale_z=100.0)
#V.coloring['stage'] = True
#V.coloring['elevation'] = False
#V.setup['elevation']=True
#V.updating['stage']=True
#V.qcolor['stage'] = (0.1,0.4,0.99)


#V.start()
#V.idle.wait()
#V.idle.clear()





for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
    if myid == 0:
        domain.write_time()
        #print_test_stats(domain, tri_full_flag)

#    V.redraw_ready.set()
#    V.idle.wait()
#    V.idle.clear()
#    V.unpaused.wait()


#print 'P%d: That took %.2f seconds' %(myid, time.time()-t0)
#print 'P%d: Communication time %.2f seconds' %(myid, domain.communication_time)
#print 'P%d: Reduction Communication time %.2f seconds' %(myid, domain.communication_reduce_time)
#print 'P%d: Broadcast time %.2f seconds' %(myid, domain.communication_broadcast_time)



if myid == 0:
    print 'That took %.2f seconds' %(time.time()-t0)
    print 'Communication time %.2f seconds'%domain.communication_time
    print 'Reduction Communication time %.2f seconds'%domain.communication_reduce_time
    print 'Broadcast time %.2f seconds'%domain.communication_broadcast_time


pypar.finalize()
