"""
Simple water flow example using ANUGA

Supercritical flow over with 'short transition' decrease 

"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
import sys
import anuga
from anuga import Domain as Domain
from anuga import myid, finalize, distribute
from math import cos
from numpy import zeros, ones, float
from time import localtime, strftime, gmtime
from anuga import Inlet_operator



#-------------------------------------------------------------------------------
# Copy scripts to time stamped output directory and capture screen
# output to file
#-------------------------------------------------------------------------------
time = strftime('%Y%m%d_%H%M%S',localtime())

output_dir = '.'
output_file = 'numerical_depth_expansion'

args = anuga.get_args()
alg = args.alg
verbose = args.verbose


#------------------------------------------------------------------------------
# Setup domain
#------------------------------------------------------------------------------
dx = 0.1
dy = dx
L = 25.
W = 3*dx

if myid == 0:
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
    def elevation(x,y):
        z_b = zeros(len(x))
        e_w=2.0 # Width of expansion region (where depth changes)
        for i in range(len(x)):
            if (x[i] <= 10.0-e_w/2.0):
                z_b[i] = 0.2
            elif (10.0-e_w/2.0 < x[i] < 10.0+e_w/2.0):
                z_b[i] = 0.2 - 0.2*(x[i]-(10.-e_w/2.0))/e_w
            else:
                z_b[i] = 0.0
        return z_b
    domain.set_quantity('elevation',elevation)
    domain.set_quantity('friction', 0.0)
    
    
    def height(x,y):
        return 2.0*ones(len(x))
    domain.set_quantity('stage', height)
else:
    domain = None
    
#---------------------------
# Create Parallel Domain
#---------------------------    
domain = distribute(domain)

#-----------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
Br = anuga.Reflective_boundary(domain)      # Solid reflective wall

outflow_w = 1.0
outflow_q = 1.0
#Bd_out = anuga.Dirichlet_boundary([outflow_w, outflow_q, 0.]) # Constant outflow boundary values
def waveform(t):
    return outflow_w
Bd_out=anuga.Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(domain,waveform)

# Inflow boundary -- use discharge inflow, and reflective wall (so we can control the discharge)
line = [ [ 0., 0.], [0., W] ]
Qin = 1.0*W
Inlet_operator(domain, line, Qin)

# Associate boundary tags with boundary objects
domain.set_boundary({'left': Br, 'right': Bd_out, 'top': Br, 'bottom': Br})




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
for t in domain.evolve(yieldstep = 2.0, finaltime = 200.):
    #print(domain.timestepping_statistics(track_speeds=True))
    if myid == 0 and verbose: print(domain.timestepping_statistics())


domain.sww_merge(delete_old=True)


finalize()
