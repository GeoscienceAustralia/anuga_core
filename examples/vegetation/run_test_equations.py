"""Simple water flow example using ANUGA

Water flowing down a channel with a topography that varies with time
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
from anuga import rectangular_cross
from anuga import Domain
from anuga import Reflective_boundary
from anuga import Dirichlet_boundary
from anuga import Transmissive_boundary
from anuga import Time_boundary

from anuga import Quantity

import numpy as num

#===============================================================================
# Setup Functions
#===============================================================================


#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------
def topography(x,y):
    """Complex topography defined by a function of vectors x and y."""
    z = 10.-x/50
    
    return z

def depth(x,y):
    """Complex topography defined by a function of vectors x and y."""
    z = topography(x,y)
            
    return z + 0.1

#===============================================================================
# Setup and Run Model
#===============================================================================


#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------
length = 15.
width = 4.
dx = dy = 0.25 #.1           # Resolution: Length of subdivisions on both axes

points, vertices, boundary = rectangular_cross(int(length/dx), int(width/dy),
                                               len1=length, len2=width)

evolved_quantities = ['stage', 'xmomentum', 'ymomentum', 'elevation']
                                               
domain = Domain(points, vertices, boundary, evolved_quantities=evolved_quantities)
domain.set_flow_algorithm('DE0')
domain.set_name('veg') # Output name
# domain.set_store_vertices_uniquely(True)

# print domain.statistics()

domain.set_quantities_to_be_stored({'elevation': 2,
                                    'stage': 2, 
                                    'xmomentum': 2,
                                    'ymomentum': 2})


domain.set_quantity('elevation', topography)           # elevation is a function
domain.set_quantity('friction', 0.01)                  # Constant friction
domain.set_quantity('stage', topography)   # Dry initial condition

#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
Bi = Dirichlet_boundary([11.5, 0, 0])          # Inflow
Br = Reflective_boundary(domain)              # Solid reflective wall
# Bo = Dirichlet_boundary([0, 0., 0.])           # Outflow
Bo = Transmissive_boundary(domain)

domain.set_boundary({'left': Bi, 'right': Bo, 'top': Br, 'bottom': Br})

#------------------------------------------------------------------------------
# Setup erosion operator in the middle of dam
#------------------------------------------------------------------------------
print 'Set up vegetation to test...'


from anuga.operators.vegetation_operator import Vegetation_operator

# create operator
op1 = Vegetation_operator(domain)
# op1.set_inflow_concentration(0.02)

alpha = 0.4
Quantity(domain, name='veg_diameter', register=True)
domain.set_quantity('veg_diameter', 0.00064)

Quantity(domain, name='veg_spacing', register=True)
domain.set_quantity('veg_spacing', 0.15)

#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep=0.5, finaltime=30.):
    domain.print_timestepping_statistics()
    #domain.print_operator_timestepping_statistics()










