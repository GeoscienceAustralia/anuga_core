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
from anuga import myid, finalize, distribute


#-------------------------------------------------------------------------------
# Copy scripts to time stamped output directory and capture screen
# output to file
#-------------------------------------------------------------------------------
time = strftime('%Y%m%d_%H%M%S',localtime())

#output_dir = 'steep_island'+time
output_dir = '.'
output_file = 'steep_island'

args = anuga.get_args()
alg = args.alg
verbose = args.verbose

dx = 1.
dy = dx
L = 2000.
W = 5*dx

def stage_flat(x,y):
    w=zeros(len(x))
    for i in range(len(x)):
        w[i]=4.5
    return w

def bed_elevation(x,y):
    z=zeros(len(x))
    for i in range(len(x)):
        if 0 <= x[i] < 200.0:
            z[i] = -0.01*(x[i]-200) + 4.0
        elif 200.0 <= x[i] < 300.0:
            z[i] = -0.02*(x[i]-200) + 4.0
        elif 300.0 <= x[i] < 400.0:
            z[i] = -0.01*(x[i]-300) + 2.0
        elif 400.0 <= x[i] < 550.0:
            z[i] = (-1/75.0)*(x[i]-400.0) + 2.0
        elif 550.0 <= x[i] < 700.0:
            z[i] = (1/11250)*(x[i]-550)*(x[i]-550)
        elif 700.0 <= x[i] < 800.0:
            z[i] = 0.03*(x[i]-700)
        elif 800.0 <= x[i] < 900.0:
            z[i] = -0.03*(x[i]-800) + 3.0
        elif 900.0 <= x[i] < 1000.0:
            z[i] = 6.0
        elif 1000.0 <= x[i] < 1400.0:
            z[i] = (-1.0/20000)*(x[i]-1000)*(x[i]-1400)
        elif 1400.0 <= x[i] < 1500.0:
            z[i] = 0.0
        elif 1500.0 <= x[i] < 1700.0:
            z[i] = 3.0
        elif 1700.0 <= x[i] < 1800.0:
            z[i] = -0.03*(x[i]-1700) + 3.0
        else:
            z[i] = (4.5/40000)*(x[i]-1800)*(x[i]-1800) + 2.0
    return z 

#------------------------------------------------------------------------------
# Setup sequential domain
#------------------------------------------------------------------------------
if myid == 0:
    # structured mesh
    points, vertices, boundary = anuga.rectangular_cross(int(L/dx), int(W/dy), L, W, (0.0, 0.0))
 
    domain = Domain(points, vertices, boundary) 
    
    domain.set_name(output_file)                
    domain.set_datadir(output_dir) 
    domain.set_flow_algorithm(alg)
    
    #------------------------------------------------------------------------------
    # Setup initial conditions
    #------------------------------------------------------------------------------
    domain.set_quantity('friction', 0.0)
    domain.set_quantity('stage', stage_flat)
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
    if myid == 0: print(domain.timestepping_statistics())


domain.sww_merge(delete_old=True)

finalize()


#

