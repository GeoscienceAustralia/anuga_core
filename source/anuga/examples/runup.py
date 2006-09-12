"""Simple water flow example using ANUGA

Water driven up a linear slope and time varying boundary,
similar to a beach environment
"""


#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

from anuga.pyvolution.mesh_factory import rectangular_cross
from anuga.pyvolution.shallow_water import Domain
from anuga.pyvolution.shallow_water import Reflective_boundary
from anuga.pyvolution.shallow_water import Dirichlet_boundary
from anuga.pyvolution.shallow_water import Time_boundary
from anuga.pyvolution.shallow_water import Transmissive_boundary


#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------

points, vertices, boundary = rectangular_cross(10, 10) # Basic mesh

domain = Domain(points, vertices, boundary) # Create domain
domain.set_name('runup')                    # Output to bedslope.sww


#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------

def topography(x,y):
    return -x/2                             # linear bed slope
    #return x*(-(2.0-x)*.5)                  # curved bed slope

domain.set_quantity('elevation', topography) # Use function for elevation
domain.set_quantity('friction', 0.1)         # Constant friction 
domain.set_quantity('stage', -.4)            # Constant negative initial stage


#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------

from math import sin, pi, exp
Br = Reflective_boundary(domain)      # Solid reflective wall
Bt = Transmissive_boundary(domain)    # Continue all values on boundary 
Bd = Dirichlet_boundary([-0.2,0.,0.]) # Constant boundary values
Bw = Time_boundary(domain=domain,     # Time dependent boundary  
                   f=lambda t: [(.1*sin(t*2*pi)-0.3) * exp(-2*t), 0.0, 0.0])

# Associate boundary tags with boundary objects
domain.set_boundary({'left': Br, 'right': Bd, 'top': Br, 'bottom': Br})


#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------

for t in domain.evolve(yieldstep = 0.1, finaltime = 10.0):
    domain.write_time()
    


