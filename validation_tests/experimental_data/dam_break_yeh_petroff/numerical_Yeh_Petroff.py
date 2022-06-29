"""Simple water flow example using ANUGA

Water driven up a linear slope and time varying boundary,
similar to a beach environment
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
import sys
import anuga
from anuga import Domain
from math import cos
from numpy import zeros
from time import localtime, strftime, gmtime

from anuga.geometry.polygon import inside_polygon, is_inside_triangle

from anuga import distribute, myid, numprocs, finalize, barrier

#================================================================================
# Setup parameters and globally used functions
#================================================================================

args = anuga.get_args()
alg = args.alg
verbose = args.verbose

output_file = 'dam_break'
output_dir = '.'

time = strftime('%Y%m%d_%H%M%S',localtime())

dx = 0.01
dy = dx
L = 1.6
W = 0.61

#-------------------------------
# Setup initial conditions
#-------------------------------
def stage(x,y):
    h = zeros(len(x), float)
    for i in range(len(x)):
        if x[i] <= 0.4:
            h[i] = 0.3
        elif x[i] <= 0.4 + 0.5:
            h[i] = 0.01
        elif x[i] <= 0.4 + 0.5 + 0.12:
            if 0.25 <= y[i] <= 0.25 + 0.12:
                h[i] = 0.75
            else:
                h[i] = 0.01
        else:
            h[i] = 0.01
    return h

def elevation(x,y):
    z = zeros(len(x), float)
    for i in range(len(x)):
        if 0.4 + 0.5 <= x[i] <= 0.4 + 0.5 + 0.12:
            if 0.25 <= y[i] <= 0.25 + 0.12:
                z[i] = 0.75
    return z



#================================================================================
# create sequential domain
#================================================================================
if myid == 0:
    # structured mesh
    domain = anuga.rectangular_cross_domain(int(L/dx), int(W/dy), L, W, (0.0, 0.0))

    domain.set_name(output_file)                
    domain.set_datadir(output_dir) 

    #------------------------------------------------------------------------------
    # Setup Algorithm, either using command line arguments
    # or override manually yourself
    #------------------------------------------------------------------------------
    domain.set_flow_algorithm(alg)

    domain.set_quantity('stage', stage, location='centroids')
    domain.set_quantity('elevation', elevation, location='centroids')
    domain.set_quantity('friction', 0.03, location='centroids')

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
Bd = anuga.Dirichlet_boundary([1,0.,0.])    # Constant boundary values

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
for t in domain.evolve(yieldstep=0.1, finaltime=3.0):
    if myid == 0 and verbose:
        print(domain.timestepping_statistics())


domain.sww_merge(delete_old=True)

finalize()

