"""Simple water flow example using ANUGA

Water flowing down a channel with a topography that varies with time
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
import numpy
from anuga import rectangular_cross
from anuga import Domain
from anuga import Reflective_boundary
from anuga import Dirichlet_boundary
from anuga import Time_boundary

#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------
length = 24.
width = 5.
dx = dy = 0.2 #.1           # Resolution: Length of subdivisions on both axes

points, vertices, boundary = rectangular_cross(int(length/dx), int(width/dy),
                                               len1=length, len2=width)
domain = Domain(points, vertices, boundary)
domain.set_name()
print (domain.statistics())
domain.set_quantities_to_be_stored({'elevation': 2,
                                    'stage': 2,
                                    'xmomentum': 2,
                                    'ymomentum': 2})

#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------
def topography_dam(x,y):
    """Complex topography defined by a function of vectors x and y."""

    z = -x/100

    # Step
    id = (2 < x) & (x < 4)
    z[id] += 0.4 - 0.05*y[id]


    # Permanent pole
    #id = (x - 8)**2 + (y - 2)**2 < 0.4**2
    #z[id] += 1


    # Dam
    id = (12 < x) & (x  < 13)
    z[id] += 0.4


            
    return z

def topography_dam_break(x,y):
    """Complex topography defined by a function of vectors x and y."""

    z = -x/100
    
    # Step
    id = (2 < x) & (x < 4)
    z[id] += 0.4 - 0.05*y[id]


    # Permanent pole
    #id = (x - 8)**2 + (y - 2)**2 < 0.4**2
    #z[id] += 1

    # Dam with hole
    id = (12 < x) & (x < 13) and (y > 3.0)
    z[id] += 0.4

    id = (12 < x) & (x < 13) and (y < 2.0)
    z[id] += 0.4


    return z


domain.set_quantity('elevation', topography_dam)       # elevation is a function
domain.set_quantity('friction', 0.01)                  # Constant friction
domain.set_quantity('stage', expression='elevation')   # Dry initial condition

#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
Bi = Dirichlet_boundary([0.4, 0, 0])          # Inflow
Br = Reflective_boundary(domain)              # Solid reflective wall
Bo = Dirichlet_boundary([-5, 0, 0])           # Outflow

domain.set_boundary({'left': Bi, 'right': Bo, 'top': Br, 'bottom': Br})

#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------

from anuga.operators.set_elevation import Set_elevation
op1 = Set_elevation(domain, elevation=lambda x,y : -x/100, radius=1.0, center = (12.5,3.0))


#dam_break = False

for t in domain.evolve(yieldstep=0.1, finaltime=30.0):
    domain.print_timestepping_statistics()
    domain.print_operator_timestepping_statistics()

    if numpy.allclose(t, 10.0):
        print ('changing elevation')
        op1()






