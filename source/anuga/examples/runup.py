"""Simple water flow example using ANUGA

Water driven up a linear slope and time varying boundary,
similar to a beach environment
"""


#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
from anuga.shallow_water import Domain
from anuga.shallow_water import Reflective_boundary
from anuga.shallow_water import Dirichlet_boundary
from anuga.shallow_water import Time_boundary
from anuga.shallow_water import Transmissive_Momentum_Set_Stage_boundary


#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------
N=10
points, vertices, boundary = rectangular_cross(N, N) # Basic mesh

domain = Domain(points, vertices, boundary) # Create domain
domain.set_name('runup')                    # Output to bedslope.sww


#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------
def topography(x,y):
    return -x/2                             # linear bed slope
    

domain.set_quantity('elevation', topography) # Use function for elevation
domain.set_quantity('friction', 0.0)         # Constant friction 
domain.set_quantity('stage', -.4)            # Constant negative initial stage


#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
from math import sin, pi, exp

def waveform(t): 
    return (0.1*sin(t*2*pi)-0.8) * exp(-2*t)

Bw = Transmissive_Momentum_Set_Stage_boundary(domain, waveform)
Br = Reflective_boundary(domain)      # Solid reflective wall
Bd = Dirichlet_boundary([-0.3,0,0])

# Associate boundary tags with boundary objects
domain.set_boundary({'left': Br, 'right': Bd, 'top': Br, 'bottom': Br})


#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep = 0.1, finaltime = 5.0):
    domain.write_time()
    


