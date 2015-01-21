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

import time

import numpy as num

from print_stats import print_test_stats, build_full_flag

#----------------------------
# Sequential interface
#---------------------------
from anuga import Domain
from anuga import Transmissive_boundary, Reflective_boundary

#----------------------------
# Parallel interface
#---------------------------
from anuga_parallel import Parallel_shallow_water_domain
from anuga_parallel import parallel_rectangle
from anuga_parallel import myid, numprocs, finalize, get_processor_name


processor_name = get_processor_name()

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

domain = Parallel_shallow_water_domain(points, vertices, boundary,
                                       full_send_dict  = full_send_dict,
                                       ghost_recv_dict = ghost_recv_dict)


print 'after parallel domain'



domain.set_name('sw_rectangle')

#Boundaries
T = Transmissive_boundary(domain)
R = Reflective_boundary(domain)


domain.set_boundary( {'left': R, 'right': R, 'bottom': R, 'top': R, 'ghost': None} )


print 'after set_boundary'



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

print 'after check_integrity'

domain.set_quantity('stage', Set_Stage(0.2, 0.4, 0.25, 0.75, 1.0, 0.00))


print 'after set quantity'

# Set Evolve parameters
domain.set_default_order(2)
domain.set_timestepping_method('rk2')
domain.set_CFL(0.7)
domain.set_beta(1.5)

#print domain.get_timestepping_method()
#domain.use_edge_limiter = True
#domain.tight_slope_limiters = True
#domain.use_centroid_velocities = False


print 'after evolve parameters'


import pdb; pdb.set_trace()

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



yieldstep = 0.05
finaltime = 2.0

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

finalize()
