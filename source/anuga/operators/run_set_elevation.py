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

#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------
length = 24.
width = 5.
dx = dy = 0.2 #.1           # Resolution: Length of subdivisions on both axes

points, vertices, boundary = rectangular_cross(int(length/dx), int(width/dy),
                                               len1=length, len2=width)
domain = Domain(points, vertices, boundary)
domain.set_name('set_elevation') # Output name
print domain.statistics()
domain.set_quantities_to_be_stored({'elevation': 2,
                                    'stage': 2})

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
    return z


def pole_increment(x,y,t):
    """This provides a small increment to a pole located mid stream
    For use with variable elevation data
    """


    z = 0.0*x
    

    if t<10.0:
        return z
    

    N = len(x)
    for i in range(N):
        # Pole 1
        if (x[i] - 12)**2 + (y[i] - 3)**2 < 0.4**2:
            z[i] += 0.1

    for i in range(N):
        # Pole 2
        if (x[i] - 14)**2 + (y[i] - 2)**2 < 0.4**2:
            z[i] += 0.05

    return z


def pole(t):

    if t<10:
        return 0.0
    elif t>15:
        return 0.0
    else:
        return 1.0


domain.set_quantity('elevation', topography)           # elevation is a function
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

from anuga.operators.set_elevation_operators import Circular_set_elevation_operator

op1 = Circular_set_elevation_operator(domain, elevation=pole, radius=0.5, center = (12.0,3.0))

growing = False
shrinking = False

done = False
for t in domain.evolve(yieldstep=0.1, finaltime=40.0):
    domain.print_timestepping_statistics()
    domain.print_operator_timestepping_statistics()




