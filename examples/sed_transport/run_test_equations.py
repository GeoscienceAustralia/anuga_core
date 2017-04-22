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
 # 
#     z[y < 0.2] = 13.
#     z[y > width - 0.2] = 13.
    
#     N = len(x)
#     for i in range(N):
#         # Step
# #         if 2 < x[i] < 4:
# #             z[i] += 0.4 - 0.05*y[i]
# 
#         # Permanent pole
#         if (x[i] - 2.5)**2 + (y[i] - 2.5)**2 < 0.3**2:
#             z[i] += 5
#         
#         if (x[i] - 3)**2 + (y[i] - 1.5)**2 < 0.4**2:
#             z[i] += 5
#             
#         if (x[i] - 6)**2 + (y[i] - 2.5)**2 < 0.6**2:
#             z[i] += 5
#             
#         if (x[i] - 10)**2 + (y[i] - 2)**2 < 0.5**2:
#             z[i] += 5
            
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

evolved_quantities = ['stage', 'xmomentum', 'ymomentum', 'elevation', 'concentration']
                                               
domain = Domain(points, vertices, boundary, evolved_quantities=evolved_quantities)
domain.set_flow_algorithm('DE0')
domain.set_name('test_equations') # Output name
# domain.set_store_vertices_uniquely(True)


# print domain.statistics()

domain.set_quantities_to_be_stored({'elevation': 2,
                                    'stage': 2,# 
#                                     'xmomentum': 2,
#                                     'ymomentum': 2,
                                    'concentration': 2})

domain.set_quantity('concentration', 0.01)
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
print 'Set up Erosion Area to test...'

from anuga.operators.sed_transport_operator import Sed_transport_operator

# create operator
op1 = Sed_transport_operator(domain)
# op1.set_inflow_concentration(0.02)


#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep=0.5, finaltime=10.):
    domain.print_timestepping_statistics()
    #domain.print_operator_timestepping_statistics()










