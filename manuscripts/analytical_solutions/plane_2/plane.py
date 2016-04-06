"""Simple water flow example using ANUGA

Water flowing down a ramp
Sediement transport on
Run to steady state


This run is very wide to minimize issues with interactions with the walls
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

import numpy as num

#===============================================================================
# Setup Functions
#===============================================================================


#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------
def topography(x,y):
    """Complex topography defined by a function of vectors x and y."""
    print ' Create topography....'
    z = 10.-x/50
            
    return z
    
def stage(x,y):
    """Complex topography defined by a function of vectors x and y."""
    print ' Create topography....'
    z = 10.-x/50
            
    return z + 0.4



#===============================================================================
# Setup and Run Model
#===============================================================================


#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------
print ' Set up Domain first...'
length = 15.
width = 10.
dx = dy = 0.5 #.1           # Resolution: Length of subdivisions on both axes

points, vertices, boundary = rectangular_cross(int(length/dx), int(width/dy),
                                               len1=length, len2=width)

evolved_quantities = ['stage', 'xmomentum', 'ymomentum', 'elevation', 'concentration']
                                               
domain = Domain(points, vertices, boundary, evolved_quantities=evolved_quantities)
domain.set_flow_algorithm('DE0')
domain.set_name('plane') # Output name
# domain.set_store_vertices_uniquely(True)

print domain.statistics()

domain.set_quantities_to_be_stored({'elevation': 2,
                                    'stage': 2,# 
                                    'xmomentum': 2,
                                    'concentration': 2})

domain.set_quantity('concentration', 0.01)
domain.set_quantity('elevation', topography)           # elevation is a function
domain.set_quantity('friction', 0.0)                  # Constant friction
domain.set_quantity('stage', stage)   # Dry initial condition
# domain.set_quantity('xmomentum', 0.1)   # Dry initial condition
#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
Bi = Dirichlet_boundary([10.42, 8, 0])          # Inflow
Br = Reflective_boundary(domain)              # Solid reflective wall
# Bo = Dirichlet_boundary([10., 0., 0.])           # Outflow
Bo = Transmissive_boundary(domain)

domain.set_boundary({'left': Bi, 'right': Bo, 'top': Bo, 'bottom': Bo})

#------------------------------------------------------------------------------
# Setup erosion operator in the middle of dam
#------------------------------------------------------------------------------
print 'Set up Erosion Area to test...'

from anuga.operators.sed_transport_operator import Sed_transport_operator

op1 = Sed_transport_operator(domain)

x = domain.quantities['x'].centroid_values
#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep=1, finaltime=60.):

#     if t > 0:
#         print op1.conc.max(), op1.conc.min()
#         print op1.edot.max(), op1.edot.min()
#         print op1.ddot.max(), op1.ddot.min()
#         print '-' * 10

    if t>0:
        
        num.save('dstar/dstar' + str(t) + '.npy', [x, op1.d_star])

    domain.print_timestepping_statistics()
    #domain.print_operator_timestepping_statistics()










