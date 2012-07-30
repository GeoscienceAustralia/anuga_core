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
from anuga import Time_boundary

from anuga.operators.erosion_operators import Polygonal_erosion_operator

#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------
length = 24.
width = 5.
dx = dy = 0.1 #.1           # Resolution: Length of subdivisions on both axes

points, vertices, boundary = rectangular_cross(int(length/dx), int(width/dy),
                                               len1=length, len2=width)
domain = Domain(points, vertices, boundary)
domain.set_name('polygon_erosion') # Output name
print domain.statistics()
domain.set_quantities_to_be_stored({'elevation': 2,
                                    'stage': 2,
                                    'xmomentum': 2,
                                    'ymomentum': 2})

#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------
def topography(x,y):
    """Complex topography defined by a function of vectors x and y."""

    z = -x/100

    N = len(x)
    for i in range(N):
        # Step
        if 2 < x[i] < 4:
            z[i] += 0.4 - 0.05*y[i]

        # Permanent pole
        if (x[i] - 8)**2 + (y[i] - 2)**2 < 0.4**2:
            z[i] += 1

        # Dam
        if 12 < x[i] < 13:
            z[i] += 0.4

            
    return z



domain.set_quantity('elevation', topography)           # elevation is a function
domain.set_quantity('friction', 0.01)                  # Constant friction
domain.set_quantity('stage', expression='elevation')   # Dry initial condition

#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
Bi = Dirichlet_boundary([0.5, 0, 0])          # Inflow
Br = Reflective_boundary(domain)              # Solid reflective wall
Bo = Dirichlet_boundary([-5, 0, 0])           # Outflow

domain.set_boundary({'left': Bi, 'right': Bo, 'top': Br, 'bottom': Br})

#------------------------------------------------------------------------------
# Setup erosion operator in the middle of dam
#------------------------------------------------------------------------------

polygon1 = [ [12., 2.0], [13., 2.0], [13., 3.0], [12., 3.0] ]
op1 = Polygonal_erosion_operator(domain, threshold=0.0, polygon=polygon1)


#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------

for t in domain.evolve(yieldstep=0.2, finaltime=30.0):
    domain.print_timestepping_statistics()
    domain.print_operator_timestepping_statistics()

#    elev_v = domain.get_quantity('elevation').vertex_values
#    stage_v = domain.get_quantity('stage').vertex_values
#    elev_c = domain.get_quantity('elevation').centroid_values
#    stage_c = domain.get_quantity('stage').centroid_values
#    ind = op1.indices
#    from pprint import pprint
#    pprint(ind)
#    print elev_v[ind]
#    print stage_v[ind]
#    print elev_c[ind]
#    print stage_c[ind]








