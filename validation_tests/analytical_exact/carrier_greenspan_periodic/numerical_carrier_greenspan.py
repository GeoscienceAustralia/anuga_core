"""
Periodic water flows  using ANUGA,
where water driven up a linear sloping beach and time varying boundary.
Ref1: Carrier and Greenspan, Journal of Fluid Mechanics, 1958
Ref2: Mungkasi and Roberts, Int. J. Numerical Methods in Fluids, 2012
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
import sys
import anuga
from anuga import Domain as Domain
from anuga import myid, finalize, distribute

from math import cos
from numpy import zeros, array
from time import localtime, strftime, gmtime
from scipy.optimize import fsolve
from math import sin, pi, exp, sqrt
from scipy.special import jn


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


args = anuga.get_args()
alg = args.alg
verbose = args.verbose

#------------------------------------------------------------------------------
# Setup domain
#------------------------------------------------------------------------------
#DIMENSIONAL PARAMETERS
dx  = 100.
dy  = dx
L   = 5e4         # Length of channel (m)
W   = 5*dx        # Width of channel (m)       
h0  = 5e2          # Height at origin when the water is still
Tp  = 900.0       # Period of oscillation
a   = 1.0         # Amplitude at origin
g   = 9.81        # Gravity

# Bessel functions
def j0(x):
    return jn(0.0, x)

def j1(x):
    return jn(1.0, x)

def elevation(x,y):
    N = len(x)    
    z = zeros(N, float)
    for i in range(N):
        z[i] = (h0/L)*x[i] - h0
    return z

def height(x,y):
    N = len(x)    
    h = zeros(N, float)
    for i in range(N):
        h[i] = h0 - (h0/L)*x[i]
        if h[i] < 0.0:
            h[i] = 0.0
    return h

def stage(x,y):
    h = height(x,y)
    z = elevation(x,y)
    return h+z

#DIMENSIONLESS PARAMETERS
eps = a/h0
T = Tp*sqrt(g*h0)/L
A = eps/j0(4.0*pi/T)

#===============================================================================
# Create sequential domain
#===============================================================================
if myid == 0:
    # structured mesh

    points, vertices, boundary = \
      anuga.rectangular_cross(int(1.1*L/dx), int(W/dy), 1.1*L, W, (0.0, 0.0))

    domain = Domain(points, vertices, boundary) 
    domain.set_name(output_file)                
    domain.set_datadir(output_dir)

    domain.set_flow_algorithm(alg)


    #------------------------------------------------------------------------------
    # Setup initial conditions
    #------------------------------------------------------------------------------
    domain.set_quantity('friction', 0.0)
    domain.set_quantity('elevation', elevation)
    domain.set_quantity('height', height)
    domain.set_quantity('stage', stage)

else:

    domain = None

#==================================================================================
# Parallel Domain
#==================================================================================
domain = distribute(domain)

#-----------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
##def shore(t):
##    def g(u):
##        return u + 2.0*A*pi/T*sin(2.0*pi/T*(t+u))    
##    u = fsolve(g,0.0)
##    xi = -0.5*u*u + A*cos(2.0*pi/T*(t+u))
##    position = 1.0 + xi
##    return position, u           # dimensionless

def prescribe(x,t): 
    q = zeros(2)
    def fun(q):                  # Here q=(w, u)
        f = zeros(2)
        f[0] = q[0] + 0.5*q[1]**2.0 - A*j0(4.0*pi/T*(1.0+q[0]-x)**0.5)*cos(2.0*pi/T*(t+q[1]))
        f[1] = q[1] + A*j1(4.0*pi/T*(1.0+q[0]-x)**0.5)*sin(2.0*pi/T*(t+q[1]))/(1+q[0]-x)**0.5
        return f
    q = fsolve(fun,q)    
    return q[0], q[1]            # dimensionless

def f_CG(t):
    h0 = 5e2 
    timing = t*sqrt(g*h0)/L      # dimensionless
    w, u = prescribe(0.0,timing) # dimensionless
    w0 = w*h0                    # dimensional
    u0 = u*sqrt(g*h0)            # dimensional
    z0 = -h0                     # dimensional
    h0 = w0 - z0                 # dimensional
    p0 = u0 * h0                 # dimensional
    #[    'stage', 'Xmomentum', 'Ymomentum']    
    return [w0,  p0, 0.0]        # dimensional
    #return w0

Br = anuga.Reflective_boundary(domain)      # Solid reflective wall
Bt = anuga.Transmissive_boundary(domain)    # Continue all values on boundary 
#Bd = anuga.Dirichlet_boundary([1,0.,0.])    # Constant boundary values
#BTime = anuga.Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(domain, f_CG)
BTime = anuga.Time_boundary(domain,f_CG)

# Associate boundary tags with boundary objects
domain.set_boundary({'left': BTime, 'right': Bt, 'top': Br, 'bottom': Br})

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
for t in domain.evolve(yieldstep = Tp/48., finaltime = 30*Tp):
    if myid == 0 and verbose: print(domain.timestepping_statistics())


domain.sww_merge(delete_old=True)

finalize()
