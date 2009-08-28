#########################################################
#
#  Main file for parallel mesh testing.
#
#  This is a modification of the run_parallel_advection.py
# file.
#
#
#  Authors: Linda Stals, Steve Roberts and Matthew Hardy,
# June 2005
#
#
#
#########################################################

import pypar    # The Python-MPI interface
import time

#from Numeric import array
# pmesh
import numpy as num

from print_stats import print_test_stats, build_full_flag

from anuga.shallow_water import Domain
from parallel_shallow_water import Parallel_Domain


# mesh partition routines
from parallel_meshes import parallel_rectangle

###############################
# Read in processor information
###############################
numprocs = pypar.size()
myid = pypar.rank()
processor_name = pypar.get_processor_name()

M = 50
N = M*numprocs


if myid == 0:
    print 'N == %d' %N

points, vertices, boundary, full_send_dict, ghost_recv_dict =\
        parallel_rectangle(N, M, len1_g=1.0*numprocs, len2_g = 1.0)

print "Myid = ", myid, "no points = ", len(points), \
      "no vertices = ", len(vertices), "no boundaries = ", len(boundary)

###########################################
# Start the computations on each subpartion
###########################################

domain = Parallel_Domain(points, vertices, boundary,
                         full_send_dict  = full_send_dict,
                         ghost_recv_dict = ghost_recv_dict)


#Boundaries
from parallel_shallow_water import Transmissive_boundary, Reflective_boundary

T = Transmissive_boundary(domain)
R = Reflective_boundary(domain)


domain.set_boundary( {'left': R, 'right': R, 'bottom': R, 'top': R, 'ghost': None} )




domain.check_integrity()

class Set_Stage:
    """Set an initial condition with constant water height, for x<x0
    """

    def __init__(self, x0=0.25, x1=0.75, y0=0.0, y1=1.0, h=5.0, h0=0.0):
        self.x0 = x0
        self.x1 = x1
        self.y0 = y0
        self.y1 = y1
        self.h  = h
        self.h0 = h0

    def __call__(self, x, y):
        return self.h0 + self.h*((x>self.x0)&(x<self.x1)&(y>self.y0)&(y<self.y1))

domain.set_quantity('stage', Set_Stage(0.2, 0.4, 0.25, 0.75, 1.0, 0.00))


# Set Evolve parameters
domain.set_default_order(2)
domain.set_timestepping_method('rk2')

print domain.get_timestepping_method()

#domain.use_edge_limiter = True
#domain.tight_slope_limiters = True
#domain.use_centroid_velocities = False

domain.CFL = 1.0

domain.set_beta(0.8)



if myid == 0:
    import time
    t0 = time.time()


# Turn on the visualisation
visualise = False
if visualise:
    from anuga.visualiser import RealtimeVisualiser
    vis = RealtimeVisualiser(domain)
    vis.render_quantity_height("elevation", offset=0.001, dynamic=False)
    vis.render_quantity_height("stage", dynamic=True)
    vis.colour_height_quantity('stage', (0.2, 0.2, 0.8))
    vis.start()
    import time
    time.sleep(2.0)



yieldstep = 0.1
finaltime = 1.0

#Check that the boundary value gets propagated to all elements
for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
    if myid == 0:
        domain.write_time()
    #print_test_stats(domain, tri_full_flag)
    if visualise:
        vis.update()						



if visualise: vis.evolveFinished()

if myid == 0:
    print 'That took %.2f seconds' %(time.time()-t0)
    print 'Communication time %.2f seconds'%domain.communication_time
    print 'Reduction Communication time %.2f seconds'%domain.communication_reduce_time
    print 'Broadcast time %.2f seconds'%domain.communication_broadcast_time


if visualise: vis.join()
pypar.finalize()
