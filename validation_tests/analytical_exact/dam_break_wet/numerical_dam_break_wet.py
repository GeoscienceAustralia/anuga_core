"""Simple water flow example using ANUGA

Water driven up a linear slope and time varying boundary,
similar to a beach environment
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
import sys
import anuga
from anuga import Domain as Domain
from anuga import myid, finalize, distribute

from math import cos
from numpy import zeros
from time import localtime, strftime, gmtime


#================================================================================
# Setup parameters and globally used functions
#================================================================================
args = anuga.get_args()
alg = args.alg
verbose = args.verbose


time = strftime('%Y%m%d_%H%M%S',localtime())
output_dir = '.'
output_file = 'dam_break'

dx = 1.
dy = dx
L = 1000.
W = 5*dx


h0 = 10.0
h1 = 1.0

def height(x,y):
    z = zeros(len(x), float)
    for i in range(len(x)):
        if x[i]<=0.0:
            z[i] = h0
        else:
            z[i] = h1
    return z

#================================================================================
# create sequential domain
#================================================================================
if myid == 0:

    # structured mesh
    points, vertices, boundary = anuga.rectangular_cross(int(L/dx), int(W/dy), L, W, (-L/2.0, -W/2.0))


    domain = Domain(points, vertices, boundary) 

    domain.set_name(output_file)                
    domain.set_datadir(output_dir) 
    domain.set_flow_algorithm(alg)

    #------------------------------------------------------------------------------
    # Setup initial conditions
    #------------------------------------------------------------------------------
    domain.set_quantity('elevation',0.0)
    domain.set_quantity('friction', 0.0)
    domain.set_quantity('stage', height)
else:
    domain = None


#===================================================================================
# create parallel domain
#===================================================================================
domain = distribute(domain)

#-----------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
from math import sin, pi, exp
Br = anuga.Reflective_boundary(domain)      # Solid reflective wall
Bt = anuga.Transmissive_boundary(domain)    # Continue all values on boundary 
Bd = anuga.Dirichlet_boundary([1,0.,0.]) # Constant boundary values

# Associate boundary tags with boundary objects
domain.set_boundary({'left': Bt, 'right': Bt, 'top': Br, 'bottom': Br})

#-------------------------------------------------------------------------
# Produce a documentation of parameters
#-------------------------------------------------------------------------
from anuga.validation_utilities import save_parameters_tex
save_parameters_tex(domain)

#===================================================================================
# Evolve system through time
#===================================================================================
for t in domain.evolve(yieldstep = 0.5, finaltime = 50.):
    if myid == 0 and verbose:
        print(domain.timestepping_statistics())


domain.sww_merge(delete_old=True)

finalize()
