"""

Example of use of sediment transport and vegetation drag operators with simple water flow down a channel

M. Perignon
perignon@colorado.edu
April 2016

"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
import anuga
from anuga import rectangular_cross
from anuga import Domain
from anuga import Dirichlet_boundary
from anuga import Reflective_boundary


#===============================================================================
# Setup Functions
#===============================================================================

#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------
def topography(x,y):
    z = -x/100
    
    return z



#===============================================================================
# Setup and Run Model
#===============================================================================


#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------
length = 5.
width = 5.
dx = dy = 2


"""
Must include the process-specific quantities when creating the domain
"""        

points, vertices, boundary = rectangular_cross(int(length/dx), int(width/dy),
                                               len1=length, len2=width)
   
domain = Domain(points, vertices, boundary)



domain.set_flow_algorithm('DE0')
domain.set_name('run_simple_veg') # Output name
domain.set_store_vertices_uniquely(True)
domain.set_quantity('elevation', topography)           # elevation is a function
domain.set_quantity('stage', expression='elevation')   # Dry initial condition


"""
Store process-specific quantities with same functions
""" 
domain.set_quantities_to_be_stored({'elevation': 2,
                                    'stage': 2,
                                    'xmomentum': 2,
                                    'ymomentum': 2,
                                    'concentration': 2})
                                    
#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
max_elev = domain.quantities['elevation'].vertex_values.max()
min_elev = domain.quantities['elevation'].vertex_values.min()


Bi = Dirichlet_boundary([max_elev + 0.5, 0, 0])
Br = Reflective_boundary(domain)           # Solid reflective wall
Bo = Dirichlet_boundary([min_elev - 1, 0, 0])    # Outflow

domain.set_boundary({'left': Bi, 'right': Bo, 'top': Br, 'bottom': Br})

#------------------------------------------------------------------------------
# Setup operators
#------------------------------------------------------------------------------

from anuga.operators.vegetation_operator import Vegetation_operator

veg_op = Vegetation_operator(domain)

Quantity(domain, name='veg_diameter', register=True)
domain.set_quantity('veg_diameter', 0.00064)

Quantity(domain, name='veg_spacing', register=True)
domain.set_quantity('veg_spacing', 0.15)

#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep = 1, finaltime = 30.0):
    domain.print_timestepping_statistics()










