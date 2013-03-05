"""Simple water flow example using ANUGA

Water flowing down a channel
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

# Import standard shallow water domain and standard boundaries.
import anuga


#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------
points, vertices, boundary = anuga.rectangular_cross(10, 5,
                                               len1=10.0, len2=5.0) # Mesh

domain = anuga.Domain(points, vertices, boundary)  # Create domain
domain.set_name('channel1')                  # Output name

#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------
def topography(x, y):
    return -x/10                             # linear bed slope

domain.set_quantity('elevation', topography) # Use function for elevation
domain.set_quantity('friction', 0.01)        # Constant friction 
domain.set_quantity('stage',                 # Dry bed
                    expression='elevation')  

#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
Bi = anuga.Dirichlet_boundary([0.4, 0, 0])         # Inflow
Br = anuga.Reflective_boundary(domain)             # Solid reflective wall

domain.set_boundary({'left': Bi, 'right': Br, 'top': Br, 'bottom': Br})

#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep=0.2, finaltime=40.0):
    print domain.timestepping_statistics()

