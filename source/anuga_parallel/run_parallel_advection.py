#!/usr/bin/env python
#########################################################
#
#  Main file for parallel mesh testing. Runs an advection
# flow simulation using a rectangular mesh
#
#
#  Authors: Steve Roberts June 2005
#  Modified by Linda Stals April 2006
#
#
#########################################################



# Parallel communication routines

import pypar

#  Mesh partition routines

from parallel_meshes import parallel_rectangle

# Parallel Domain
 
from parallel_advection import Parallel_Domain
from parallel_advection import Transmissive_boundary

############################
# Set the initial conditions
############################
class Set_Stage:
    """Set an initial condition with constant water height, for x<x0
    """

    def __init__(self, x0=0.25, x1=0.5, h=1.0):
        self.x0 = x0
        self.x1 = x1
        self.h  = h

    def __call__(self, x, y):
        return self.h*((x>self.x0)&(x<self.x1))

###############################
# Read in processor information
###############################

numprocs = pypar.size()
myid = pypar.rank()
processor_name = pypar.get_processor_name()

N = 5
M = 2

N = 10
M = 5

#######################
# Partition the mesh
#######################

# Build a unit mesh, subdivide it over numproces processors with each
# submesh containing M*N nodes

points, vertices, boundary, full_send_dict, ghost_recv_dict =  \
    parallel_rectangle(N, M, len1_g=1.0)

#print 'ghost ',ghost_recv_dict
#print 'full  ',full_send_dict



print "Myid = ", myid, "no points = ", len(points), \
      "no vertices = ", len(vertices), "no boundaries = ", len(boundary)

###########################################
# Start the computations on each subpartion
###########################################

# Create advection domain with direction (1,-1)
# Initial condition is zero by default

domain = Parallel_Domain(points, vertices, boundary,
                         full_send_dict, ghost_recv_dict, velocity=[1.0, 0.0])

#print 'ghost ',ghost_recv_dict
#print 'full  ',full_send_dict


# Boundaries

T = Transmissive_boundary(domain)




domain.set_boundary( {'left': T, 'right': T, 'bottom': T, 'top': T, \
                      'ghost': None} )


# Set Evolve parameters
domain.set_default_order(2)
domain.set_timestepping_method('rk2')

print domain.get_timestepping_method()

#domain.use_edge_limiter = True
#domain.tight_slope_limiters = True
#domain.use_centroid_velocities = False

domain.CFL = 1.0

domain.set_beta(0.8)


# Ensure that the domain definitions make sense

domain.check_integrity()

# Set the inititial conditions

domain.set_quantity('stage', Set_Stage(0.2,0.4,1.0))

# Let processor 0 output some timing information

visualise = False
if visualise:
    from anuga.visualiser import RealtimeVisualiser
    vis = RealtimeVisualiser(domain)
    vis.render_quantity_height("stage", zScale = 5.0, dynamic=True)
    vis.colour_height_quantity('stage', (0.2, 0.2, 0.8))
    vis.start()
    import time
    time.sleep(2.0)
    


if myid == 0:
    import time
    t0 = time.time()

for t in domain.evolve(yieldstep = 5.0, finaltime = 30.0):
    if myid == 0:
        domain.write_time()
        
    if visualise:
        vis.update()						

        

if visualise: vis.evolveFinished()
    

# Output some computation statistics

if myid == 0:
    print 'That took %.2f seconds' %(time.time()-t0)
    print 'Communication time %.2f seconds'\
          %domain.communication_time
    print 'Reduction Communication time %.2f seconds'\
          %domain.communication_reduce_time
        
  
if visualise: vis.join()
