"""Simple water flow example using ANUGA

Water driven up a linear slope and time varying boundary,
similar to a beach environment
"""


#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular
from anuga.shallow_water import Domain
from anuga.shallow_water import Reflective_boundary
from anuga.shallow_water import Dirichlet_boundary
from anuga.shallow_water import Time_boundary
from anuga.shallow_water import Transmissive_boundary


#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------

points, vertices, boundary = rectangular(4, 1, len1=4, len2=1)
domain = Domain(points, vertices, boundary)
domain.set_name('limiter_test')
domain.set_store_vertices_uniquely(True)
domain.set_default_order(2)
domain.limit2007 = 1

#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------

def topography(x,y):
    from Numeric import zeros, size, Float
    
    z = zeros(size(x), Float)
    for i in range(len(x)):
        if x[i] < 2:
            z[i] = 1.0
        else:
            z[i] = 0.0
    return z 


domain.set_quantity('elevation', topography) # Use function for elevation
domain.set_quantity('stage', 0.5)           # Constant initial stage


#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------

Br = Reflective_boundary(domain)      # Solid reflective wall
domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------

for t in domain.evolve(yieldstep = 0.1, finaltime = 50.0):
    domain.write_time()

print 'done'    

    

