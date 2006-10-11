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
length = 40.
width = 5.
dx = dy = 1           # Resolution: Length of subdivisions on both axes

points, vertices, boundary = rectangular_cross(int(length/dx), int(width/dy), len1=length, len2=width)
domain = Domain(points, vertices, boundary)   
domain.set_name('channel_3')                  # Output name


#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------
def topography(x,y):
    z = -x/10                             # linear bed slope

    N = len(x)
    for i in range(N):

        #Step
        if 10 < x[i] < 12:
            z[i] += 0.4 - 0.05*y[i]        
        
        #Constriction
        if 27 < x[i] < 29 and y[i] > 3:
            z[i] += 2        
        
        # Pole
        if (x[i] - 34)**2 + (y[i] - 2)**2 < 0.4**2:
            z[i] += 2

    return z


domain.set_quantity('elevation', topography)            # Use function for elevation
domain.set_quantity('friction', 0.01)                   # Constant friction 
domain.set_quantity('stage', expression='elevation')    # Dry


#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
Bi = Dirichlet_boundary([0.4, 0, 0])                            # Inflow
Br = Reflective_boundary(domain)                                # Solid reflective wall
Bo = Dirichlet_boundary([-5, 0, 0])                             # Outflow

domain.set_boundary({'left': Bi, 'right': Bo, 'top': Br, 'bottom': Br})


#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep = 0.2, finaltime = 16.0):
    domain.write_time()


