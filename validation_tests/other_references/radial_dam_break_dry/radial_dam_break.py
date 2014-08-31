"""Simple water flow example using ANUGA

Water driven up a linear slope and time varying boundary,
similar to a beach environment
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
import sys
import anuga
from math import cos
from numpy import zeros, float

#------------------------------------------------------------------------------
# Setup parameters and utilitiy functions
#------------------------------------------------------------------------------
dx = 1.
dy = dx
L = 200.
W = L

h0 = 10.0
h1 = 0.0
radius = 50.0

def height(x,y):
    z = zeros(len(x), float)
    r2 = radius**2
    for i in range(len(x)):
        rad2 = x[i]**2 + y[i]**2

        if rad2 <= r2:
            z[i] = h0
        else:
            z[i] = h1
    return z

#-------------------------------------------------------------------------
# setup sequential domain
#-------------------------------------------------------------------------
if anuga.myid == 0:
    
    points, vertices, boundary = anuga.rectangular_cross(int(L/dx), int(W/dy), L, W, (-L/2.0, -W/2.0))
    domain = anuga.Domain(points, vertices, boundary)
    domain.set_name('radial_dam_break')                


    #-------------------------------------------------------------------------
    # Setup Algorithm, either using command line arguments
    # or override manually yourself
    #--------------------------------------------------------------------------
    args = anuga.get_args()
    alg = args.alg
    verbose = args.verbose
    domain.set_flow_algorithm(alg)

    #-------------------------------------------------------------------------
    # Setup initial conditions
    #-------------------------------------------------------------------------
    domain.set_quantity('elevation',0.0)
    domain.set_quantity('friction', 0.0)
    domain.set_quantity('stage', height)
else:
    domain = None

domain = anuga.distribute(domain)

#-----------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
Br = anuga.Reflective_boundary(domain)      # Solid reflective wall

# Associate boundary tags with boundary objects
domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

#-------------------------------------------------------------------------
# Produce a documentation of parameters
#-------------------------------------------------------------------------
from anuga.validation_utilities import save_parameters_tex
save_parameters_tex(domain)

#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep = 0.1, finaltime = 2.0):
     if anuga.myid == 0 and verbose:
         print domain.timestepping_statistics()

domain.sww_merge(delete_old=True)

anuga.finalize()
