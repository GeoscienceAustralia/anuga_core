"""Simple water flow example using ANUGA

Water driven up a linear slope and time varying boundary,
similar to a beach environment
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
import anuga

from math import sin, pi, exp

#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------
points, vertices, boundary = anuga.rectangular_cross(10, 10) # Basic mesh

domain = anuga.Domain(points, vertices, boundary)   # Create domain
domain.set_name('runup')                            # Output to file runup.sww
domain.set_datadir('.')                             # Use current folder

#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------
def topography(x, y):
    return -x/2                              # linear bed slope
    #return x*(-(2.0-x)*.5)                  # curved bed slope

domain.set_quantity('elevation', topography) # Use function for elevation
domain.set_quantity('friction', 0.1)         # Constant friction 
domain.set_quantity('stage', -0.4)           # Constant negative initial stage

#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
Br = anuga.Reflective_boundary(domain)      # Solid reflective wall
Bt = anuga.Transmissive_boundary(domain)    # Continue all values on boundary 
Bd = anuga.Dirichlet_boundary([-0.2,0.,0.]) # Constant boundary values
Bw = anuga.Time_boundary(domain=domain,     # Time dependent boundary  
                   f=lambda t: [(0.1*sin(t*2*pi)-0.3)*exp(-t), 0.0, 0.0])

# Associate boundary tags with boundary objects
domain.set_boundary({'left': Br, 'right': Bw, 'top': Br, 'bottom': Br})

#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep=0.1, finaltime=10.0):
    print domain.timestepping_statistics()
