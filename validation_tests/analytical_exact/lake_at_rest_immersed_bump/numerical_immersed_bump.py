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
from math import cos
from numpy import zeros, float
from time import localtime, strftime, gmtime
#from balanced_dev import *
from anuga import myid, finalize, distribute

#-------------------------------------------------------------------------------
# Copy scripts to time stamped output directory and capture screen
# output to file
#-------------------------------------------------------------------------------
time = strftime('%Y%m%d_%H%M%S',localtime())

#output_dir = 'immersed_bump'+time
output_dir = '.'
output_file = 'immersed_bump'

def bed_elevation(x,y):
    z = zeros(len(x), float)
    for i in range(len(x)):
        if 8.0 < x[i] < 12.0:
            z[i] = 0.2 - 0.05*(x[i]-10.)**2
        else:
            z[i] = 0.0
    return z    


args = anuga.get_args()
alg = args.alg
verbose = args.verbose

if myid == 0:
    #------------------------------------------------------------------------------
    # Setup domain
    #------------------------------------------------------------------------------
    dx = 1.
    dy = dx
    L = 25.
    W = 5*dx
    
    # structured mesh
    points, vertices, boundary = anuga.rectangular_cross(int(L/dx), int(W/dy), L, W, (0.0, 0.0))
    
    #domain = anuga.Domain(points, vertices, boundary) 
    domain = Domain(points, vertices, boundary) 
    
    domain.set_name(output_file)                
    domain.set_datadir(output_dir) 

    domain.set_flow_algorithm(alg)

    
    #------------------------------------------------------------------------------
    # Setup initial conditions
    #------------------------------------------------------------------------------
    domain.set_quantity('friction', 0.0)
    domain.set_quantity('stage', 0.5)    
    domain.set_quantity('elevation', bed_elevation)

else:
    
    domain = None
    
#-----------------------------------------------------------------------------
# Parallel Domain
#-----------------------------------------------------------------------------   
domain = distribute(domain)

#-----------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
from math import sin, pi, exp
Br = anuga.Reflective_boundary(domain)      # Solid reflective wall
#Bt = anuga.Transmissive_boundary(domain)    # Continue all values on boundary 
#Bd = anuga.Dirichlet_boundary([1,0.,0.]) # Constant boundary values

# Associate boundary tags with boundary objects
domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})




#------------------------------------------------------------------------------
# Produce a documentation of parameters
#------------------------------------------------------------------------------
if myid == 0:
    parameter_file=open('parameters.tex', 'w')
    parameter_file.write('\\begin{verbatim}\n')
    from pprint import pprint
    pprint(domain.get_algorithm_parameters(),parameter_file,indent=4)
    parameter_file.write('\\end{verbatim}\n')
    parameter_file.close()

#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep = 0.1, finaltime = 5.):
    #print domain.timestepping_statistics(track_speeds=True)
    if myid == 0: print domain.timestepping_statistics()

domain.sww_merge(delete_old=True)

finalize()

