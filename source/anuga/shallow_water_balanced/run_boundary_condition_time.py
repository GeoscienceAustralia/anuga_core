

import numpy as num
from swb_domain import *


"""test_boundary_condition_time

This tests that boundary conditions are evaluated
using the right time from domain.
"""

# Setup computational domain
from anuga.abstract_2d_finite_volumes.mesh_factory \
     import rectangular_cross

#--------------------------------------------------------------
# Setup computational domain
#--------------------------------------------------------------
N = 20
points, vertices, boundary = rectangular_cross(N, N)
domain = Domain(points, vertices, boundary)

#--------------------------------------------------------------
# Setup initial conditions
#--------------------------------------------------------------
domain.set_quantity('elevation', 0.0)
domain.set_quantity('friction', 0.0)
domain.set_quantity('stage', 0.0)

#--------------------------------------------------------------
# Setup boundary conditions
#--------------------------------------------------------------
# Time dependent boundary
Bt = Time_boundary(domain=domain, f=lambda t: [t, 0.0, 0.0])

Br = Reflective_boundary(domain)              # Reflective wall

domain.set_boundary({'left': Bt, 'right': Br, 'top': Br, 'bottom': Br})



interactive_visualisation = True

#===============================================================================
if interactive_visualisation:
    from anuga.visualiser import RealtimeVisualiser
    vis = RealtimeVisualiser(domain)
    vis.render_quantity_height("elevation", zScale =1.0, dynamic=True)    
    vis.render_quantity_height("stage", zScale =1.0, dynamic=True)
    vis.colour_height_quantity('stage', (1.0, 0.5, 0.5))
    vis.start()
#===============================================================================



for t in domain.evolve(yieldstep = 0.02, finaltime = 20.0):
    #print domain.get_time()
    #q = Bt.evaluate()
    
    # FIXME (Ole): This test would not have passed in
    # changeset:5846.
    #msg = 'Time boundary not evaluated correctly'
    #assert num.allclose(t, q[0]), msg

    domain.write_time()
    if interactive_visualisation:
        vis.update()

if interactive_visualisation:
    vis.evolveFinished()    
    
