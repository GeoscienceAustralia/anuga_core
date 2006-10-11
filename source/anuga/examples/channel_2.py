"""Simple water flow example using ANUGA

Water flowing down a channel
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
from anuga.shallow_water import Domain
from anuga.shallow_water import Reflective_boundary
from anuga.shallow_water import Dirichlet_boundary
from anuga.shallow_water import Time_boundary


#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------
length = 10.
width = 5.
dx = dy = 1           # Resolution: Length of subdivisions on both axes

points, vertices, boundary = rectangular_cross(int(length/dx), int(width/dy), len1=length, len2=width)
domain = Domain(points, vertices, boundary)   
domain.set_name('channel_2')                  # Output name


#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------
def topography(x,y):
    return -x/10                             # linear bed slope

domain.set_quantity('elevation', topography)            # Use function for elevation
domain.set_quantity('friction', 0.01)                   # Constant friction 
domain.set_quantity('stage', expression='elevation')    # Dry


#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
Bi = Dirichlet_boundary([0.4, 0, 0])                            # Inflow
Br = Reflective_boundary(domain)                                # Solid reflective wall
Bo = Dirichlet_boundary([-5, 0, 0])                             # Outflow

domain.set_boundary({'left': Bi, 'right': Br, 'top': Br, 'bottom': Br})


#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep = 0.2, finaltime = 40.0):
    domain.write_time()

    #print 'Stage(10,2.5) = ', domain.get_quantity('stage').get_values(interpolation_points=[[10, 2.5]])
    
    from Numeric import allclose
    #if allclose(t, 20):
    if domain.get_quantity('stage').get_values(interpolation_points=[[10, 2.5]]) > 0:        
        print 'Stage > 0: Changing boundary'
        domain.modify_boundary({'right': Bo})


import sys; sys.exit() 

import time
t0 = time.time()


s = 'for t in domain.evolve(yieldstep = 0.2, finaltime = 40.0): domain.write_time(); domain.get_quantity("stage").get_values(interpolation_points=[[10, 2.5]])'

import profile, pstats
FN = 'profile.dat'

profile.run(s, FN)

print 'That took %.2f seconds' %(time.time()-t0)

S = pstats.Stats(FN)
#S.sort_stats('time').print_stats(20)
s = S.sort_stats('cumulative').print_stats(30)

print s

