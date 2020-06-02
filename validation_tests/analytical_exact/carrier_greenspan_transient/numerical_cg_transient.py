"""
Transient water flows  using ANUGA,
where water driven up a linear sloping beach and time varying boundary.
Ref: Carrier and Greenspan, Journal of Fluid Mechanics, 1958
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
import sys
import anuga
from anuga import Domain as Domain
from numpy import zeros, array, float, hstack
from time import localtime, strftime, gmtime
from scipy.optimize import fsolve
from math import sin, pi, exp, sqrt, cos
from analytical_cg_transient import analytical_sol

from anuga import myid, finalize, distribute

import warnings
warnings.simplefilter('ignore')

#------------------------------------------------------------------------------
# Setup Algorithm, either using command line arguments
# or override manually yourself
#------------------------------------------------------------------------------
args = anuga.get_args()
alg = args.alg
verbose = args.verbose


#-------------------------------------------------------------------------------
# Copy scripts to time stamped output directory and capture screen
# output to file
#-------------------------------------------------------------------------------
time = strftime('%Y%m%d_%H%M%S',localtime())

#output_dir = 'carrier_greenspan_'+time
output_dir = '.'
output_file = 'carrier_greenspan'

#anuga.copy_code_files(output_dir,__file__)
#start_screen_catcher(output_dir+'_')


#------------------------------------------------------------------------------
# Setup domain
#------------------------------------------------------------------------------
#DIMENSIONAL PARAMETERS
dx  = 500.
dy  = dx
L   = 5e4         # Length of channel (m)
W   = 3*dx        # Width of channel (m)       
h0 = 5e2          # Height at origin when the water is still
g   = 9.81        # Gravity


def elevation(x,y):
    return h0*x/L


def stage(x,y):
    w,z,u = analytical_sol(x/L,0.0)
    return h0*w


#===============================================================================
# Create sequential domain
#===============================================================================
if myid == 0:
    # structured mesh
    points, vertices, boundary = anuga.rectangular_cross(int(0.8*L/dx), int(W/dy), 0.8*L, W, (-0.5*L, -0.5*W))
    domain = Domain(points, vertices, boundary) 
    domain.set_name(output_file)                
    domain.set_datadir(output_dir)

    domain.set_flow_algorithm(alg)

    #------------------------------------------------------------------------------
    # Setup initial conditions
    #------------------------------------------------------------------------------
    domain.set_quantity('friction', 0.0)
    domain.set_quantity('stage', stage)
    domain.set_quantity('elevation', elevation)
    
else:
    domain = None
    
#===============================================================================
# Create Parallel domain
#===============================================================================
domain = distribute(domain)

#-----------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
def f_CG(t): 
    timing = t*sqrt(g*h0)/L      # dimensionless
    w, z, u = analytical_sol(array([-0.5]),timing) # dimensionless
    wB = w*h0                # dimensional
    uB = u*sqrt(g*h0)            # dimensional
    zB = z*h0                     # dimensional
    hB = wB - zB                 # dimensional
    pB = uB * hB                 # dimensional
    #[    'stage', 'Xmomentum', 'Ymomentum']    
    return hstack( (wB,  pB, array([0.0])) )        # dimensional

Br = anuga.Reflective_boundary(domain)      # Solid reflective wall
#Bt = anuga.Transmissive_boundary(domain)    # Continue all values on boundary 
#Bd = anuga.Dirichlet_boundary([1,0.,0.])    # Constant boundary values
BTime = anuga.Time_boundary(domain,f_CG)
# Associate boundary tags with boundary objects
domain.set_boundary({'left': BTime, 'right': Br, 'top': Br, 'bottom': Br})


#===============================================================================
##from anuga.visualiser import RealtimeVisualiser
##vis = RealtimeVisualiser(domain)
##vis.render_quantity_height("stage", zScale =h0*500, dynamic=True)
##vis.colour_height_quantity('stage', (0.0, 0.5, 1.0))
##vis.start()
#===============================================================================


if myid == 0:
    #------------------------------------------------------------------------------
    # Produce a documentation of parameters
    #------------------------------------------------------------------------------
    parameter_file=open('parameters.tex', 'w')
    parameter_file.write('\\begin{verbatim}\n')
    from pprint import pprint
    pprint(domain.get_algorithm_parameters(),parameter_file,indent=4)
    parameter_file.write('\\end{verbatim}\n')
    parameter_file.close()

#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep = 1000., finaltime = 30000.):
    #print domain.timestepping_statistics(track_speeds=True)
    if myid == 0 and verbose:
        print(domain.timestepping_statistics())
    #vis.update()


domain.sww_merge(delete_old=True)

finalize()


