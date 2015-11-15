"""

Example of use of sediment transport and vegetation drag operators with simple water flow down a channel

M. Perignon
perignon@colorado.edu
July 2014

"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
import anuga
from anuga import rectangular_cross
from anuga import Domain
from anuga import Dirichlet_boundary

"""
Import operators
"""
from anuga.operators.sed_transport.sed_transport_operator import Sed_transport_operator, Vegetation_operator
"""
Import operator-specific boundaries
"""
from anuga.operators.sed_transport.sed_transport_utils import Reflective_boundary_Sed, Dirichlet_boundary_Sed


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

points, vertices, boundary = rectangular_cross(int(length/dx), int(width/dy),
                                               len1=length, len2=width)


"""
Must include the process-specific quantities when creating the domain
"""        
evolved_quantities =  ['stage', 'xmomentum', 'ymomentum', 'concentration']

other_quantities=['elevation', 'friction', 'height', 'xvelocity', \
                  'yvelocity', 'x', 'y', 'vegetation', 'diffusivity']
                  
                  
domain = Domain(points, vertices, boundary, evolved_quantities = evolved_quantities, other_quantities = other_quantities)



domain.set_flow_algorithm('1_75')
domain.set_name('run_simple_sed_transport_veg') # Output name
domain.set_store_vertices_uniquely(True)
domain.set_quantity('elevation', topography)           # elevation is a function
domain.set_quantity('stage', expression='elevation')   # Dry initial condition

print domain.statistics()

"""
Store process-specific quantities with same functions
""" 
domain.set_quantities_to_be_stored({'elevation': 2,
                                    'stage': 2,
                                    'xmomentum': 2,
                                    'ymomentum': 2,
                                    'concentration': 2,
                                    'vegetation': 1})
                                    
#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
max_elev = domain.quantities['elevation'].vertex_values.max()
min_elev = domain.quantities['elevation'].vertex_values.min()
"""
Use operator-specific Reflective and Dirichlet boundaries. Use original
Dirichlet boundary for outlet.

""" 
Bi = Dirichlet_boundary_Sed([max_elev + 0.5, 0, 0, 0.2])  # Inflow, 20% sed
Br = Reflective_boundary_Sed(domain)           # Solid reflective wall
Bo = Dirichlet_boundary([min_elev - 1, 0, 0])    # Outflow

domain.set_boundary({'left': Bi, 'right': Bo, 'top': Br, 'bottom': Br})

#------------------------------------------------------------------------------
# Setup operators
#------------------------------------------------------------------------------
"""
Uniform vegetation cover of type "1" (see vegcodes.txt for corresponding stem spacing and diameter)
""" 
domain.set_quantity('vegetation', 1)


"""
Create operators

Calculate turbulence if turbulence = True
Calculate momentum sinks if momentum_sinks = True
""" 
op1 = Sed_transport_operator(domain,
                             erosion = True,
                             deposition = True,
                             turbulence = True,
                             momentum_sinks = True,
                             verbose = True)
                             
op2 = Vegetation_operator(domain,
                          turbulence = True,
                          momentum_sinks = True,
                          vegfile = 'vegcodes.txt',
                          verbose = True)


#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep = 1, finaltime = 30.0):
    domain.print_timestepping_statistics()










