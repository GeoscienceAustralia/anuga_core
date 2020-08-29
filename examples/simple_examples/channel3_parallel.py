"""Simple water flow example using ANUGA

Water flowing down a channel with more complex topography
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
import anuga

#------------------------------------------------------------------------------
# Setup some initial info
#------------------------------------------------------------------------------
def topography(x,y):
    """Complex topography defined by a function of vectors x and y."""

    z = -x/10

    N = len(x)
    for i in range(N):
        # Step
        if 10 < x[i] < 12:
            z[i] += 0.4 - 0.05*y[i]

        # Constriction
        if 27 < x[i] < 29 and y[i] > 3:
            z[i] += 2

        # Pole
        if (x[i] - 34)**2 + (y[i] - 2)**2 < 0.4**2:
            z[i] += 2

    return z



#------------------------------------------------------------------------------
# Setup computational domain on one processor
#------------------------------------------------------------------------------
length = 40.
width = 5.
dx = dy = .1           # Resolution: Length of subdivisions on both axes


if anuga.myid == 0:
    points, vertices, boundary = anuga.rectangular_cross(int(length/dx),
                                         int(width/dy), len1=length, len2=width)
    domain = anuga.Domain(points, vertices, boundary)
    domain.set_name('channel3')                  # Output name
    domain.set_flow_algorithm('DE0')
    domain.print_statistics()


    
    domain.set_quantity('elevation', topography)           # elevation is a function
    domain.set_quantity('friction', 0.01)                  # Constant friction
    domain.set_quantity('stage', expression='elevation')   # Dry initial condition
else:
    domain = None

#------------------------------------------------------------------------------
# Distribute domain on processor 0 to to other processors
#------------------------------------------------------------------------------
#parameters = dict(ghost_layer_width=3)
domain = anuga.distribute(domain, verbose= True)


#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
Bi = anuga.Dirichlet_boundary([0.4, 0, 0])          # Inflow
Br = anuga.Reflective_boundary(domain)              # Solid reflective wall
Bo = anuga.Dirichlet_boundary([-5, 0, 0])           # Outflow

domain.set_boundary({'left': Bi, 'right': Bo, 'top': Br, 'bottom': Br})

#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep=0.1, finaltime=16.0):
    if anuga.myid == 0:
        domain.print_timestepping_statistics()


    ## if domain.get_quantity('stage').\
    ##        get_values(interpolation_points=[[10, 2.5]]) > 0:
    ##     print 'Stage > 0: Changing to outflow boundary'
    ##     domain.set_boundary({'right': Bo})


domain.sww_merge(verbose=True)

anuga.finalize()
