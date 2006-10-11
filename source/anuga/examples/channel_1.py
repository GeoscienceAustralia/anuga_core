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


#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------
length = 10.
width = 5.
dx = dy = 1           # Resolution: Length of subdivisions on both axes

points, vertices, boundary = rectangular_cross(int(length/dx), int(width/dy), len1=length, len2=width)
domain = Domain(points, vertices, boundary)   # Create domain
domain.set_name('channel_1')                  # Output name


#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------
def topography(x,y):
    return -x/10                             # linear bed slope

domain.set_quantity('elevation', topography)            # Use function for elevation
domain.set_quantity('friction', 0.01)                   # Constant friction 
domain.set_quantity('stage', expression='elevation')    # Dry
#domain.set_quantity('stage', expression='elevation + 0.1')    # Wet


#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
Bi = Dirichlet_boundary([0.4, 0, 0])                            # Inflow
Br = Reflective_boundary(domain)                                # Solid reflective wall

domain.set_boundary({'left': Bi, 'right': Br, 'top': Br, 'bottom': Br})


#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep = 0.2, finaltime = 40.0):
    domain.write_time()


