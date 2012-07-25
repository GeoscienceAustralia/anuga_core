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
domain.set_name('set_depth_friction') # Output name
print domain.statistics()


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
            
#        # Pole 2
#        if (x[i] - 14)**2 + (y[i] - 3.5)**2 < 0.4**2:
#            z[i] += 1.0

    return z




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

from anuga.operators.set_friction_operators import Set_depth_friction_operator
from anuga.operators.set_friction_operators import Polygonal_depth_friction_operator

#op1 = Set_depth_friction_operator(domain)

polygon1 = [ [12.0, 2.5], [13.5, 2.5], [13.5, 4.0], [12.0, 4.0] ]
op2 = Polygonal_depth_friction_operator(domain, friction_max = 10000, friction_min = 0.0, polygon=polygon1)



for t in domain.evolve(yieldstep=0.1, finaltime=30.0):
    domain.print_timestepping_statistics()
    domain.print_operator_timestepping_statistics()






