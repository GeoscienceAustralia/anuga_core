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
print (domain.statistics())
domain.set_quantities_to_be_stored({'elevation': 2,
                                    'stage': 2})

#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------
def topography_dam(x,y):
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

def topography_dam_break(x,y):
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
        if 12 < x[i] < 13 and y[i] > 3.0:
            z[i] += 0.4

        if 12 < x[i] < 13 and y[i] < 2.0:
            z[i] += 0.4

    return z



def pole(t):

    if t<10:
        return 0.0
    elif t>15:
        return 0.0
    else:
        return 0.05


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

#from anuga.operators.set_elevation_operators import Circular_set_elevation_operator

#op1 = Circular_set_elevation_operator(domain, elevation=pole, radius=0.5, center = (12.0,3.0))


dam_break = False

for t in domain.evolve(yieldstep=0.1, finaltime=40.0):
    domain.print_timestepping_statistics()
    domain.print_operator_timestepping_statistics()

    if t >= 10 and not dam_break:
        print ('changing elevation')

        stage_c = domain.get_quantity('stage').centroid_values
        elev_c =  domain.get_quantity('elevation').centroid_values
        height_c = stage_c - elev_c
        domain.set_quantity('elevation', topography_dam_break)
        stage_c[:] = elev_c + height_c

        dam_break = True





