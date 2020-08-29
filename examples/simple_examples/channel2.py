"""Simple water flow example using ANUGA

Water flowing down a channel with changing boundary conditions
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
import anuga

#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------
length = 10.
width = 5.
dx = dy = 1.           # Resolution: Length of subdivisions on both axes

domain = anuga.rectangular_cross_domain(int(length/dx), int(width/dy),
                                        len1=length, len2=width)
                                        

domain.set_name('channel2')                 # Output name

#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------
def topography(x,y):
    return -x/10                             # linear bed slope

domain.set_quantity('elevation', topography) # Use function for elevation
domain.set_quantity('friction', 0.01)        # Constant friction 
domain.set_quantity('stage',
                    expression='elevation')  # Dry initial condition

#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
Bi = anuga.Dirichlet_boundary([0.4, 0, 0])   # Inflow
Br = anuga.Reflective_boundary(domain)       # Solid reflective wall
Bo = anuga.Dirichlet_boundary([-5, 0, 0])    # Outflow

domain.set_boundary({'left': Bi, 'right': Br, 'top': Br, 'bottom': Br})

#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep=0.2, finaltime=40.0):
    domain.print_timestepping_statistics()

    if domain.get_quantity('stage').\
           get_values(interpolation_points=[[10, 2.5]]) > 0:        
        print('Stage > 0: Changing to outflow boundary')
        domain.set_boundary({'right': Bo})
        
