"""
Simple water flow example using ANUGA
Transcritical flow over a bump with a shock
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
import sys
import anuga
from anuga import Domain as Domain
from anuga import myid, finalize, distribute
from anuga import g
from math import cos
from numpy import zeros, ones, array
from time import localtime, strftime, gmtime
from scipy.integrate import quad
#from balanced_dev import *


#-------------------------------------------------------------------------------
# Copy scripts to time stamped output directory and capture screen
# output to file
#-------------------------------------------------------------------------------
time = strftime('%Y%m%d_%H%M%S',localtime())
output_dir = '.'
output_file = 'MacDonald'

args = anuga.get_args()
alg = args.alg
verbose = args.verbose


dx = 0.25
dy = dx
L = 100.
W = 3*dx


L     = 100.  # Channel length
q = q_s = 2.    # Steady discharge
q_up  = q_s   # Upstream discharge
#h_down= h_L   # Downstream water height

#Parameters for water height
a1 = 0.674202
a2 = 21.7112
a3 = 14.492
a4 = 1.4305
x_s=200./3   # Specified position of a shock
n  = 0.0328  # Manning's friction coefficient


def bed(x,y):
    N = len(x)
    Z = zeros((N,2))
    def func1(x): #on the left of the shock
        h = (4./g)**(1./3) * (4./3 - x/100.) - 9.*x/1000*(x/100 - 2./3)
        hx = (4./g)**(1./3)*(-0.01) - (9./1000)*(x/100 - 2./3) - (9.*x/1000)*0.01
        return (q**2/(g*h**3) - 1.)*hx - n**2*q*abs(q)/h**(10./3)
    def func2(x): #on right of the shock
        r = x/100. - 2./3.
        rx = 0.01
        h = (4./g)**(1./3) * (a1*r**4 + a1*r**3 - a2*r**2 + a3*r + a4)
        hx = (4./g)**(1./3)*(a1*4.*r**3*rx + a1*3.*r**2*rx - a2*2.*r*rx + a3*rx)
        return (q**2/(g*h**3) - 1.)*hx - n**2*q*abs(q)/h**(10./3)
    for i in range(N):
        if x[i] <= x_s:
            Z[i] = quad(func1, x[i], x_s)  #on the left of the shock
            Z[i] += quad(func2, x_s, L)    #on right of the shock
        else:
            Z[i] = quad(func2, x[i], L)    #on right of the shock
    #print('Integral = ', -Z[:,0], ' with error = ', Z[:,1])
    elevation = -1.*Z[:,0]
    return elevation


#------------------------------------------------------------------------------
# Setup sequential domain
#------------------------------------------------------------------------------
if myid == 0:
    # structured mesh
    domain = anuga.rectangular_cross_domain(int(L/dx), int(W/dy), L, W, (0.0, 0.0))

    domain.set_name(output_file)
    domain.set_datadir(output_dir)
    domain.set_flow_algorithm(alg)

    #------------------------------------------------------------------------------
    # Setup initial conditions
    #------------------------------------------------------------------------------

    domain.set_quantity('stage', 2.87870797)
    domain.set_quantity('elevation', bed)
    domain.set_quantity('friction', n)
else:

    domain = None

#--------------------------------------------------------------------
# Parallel Domain
#--------------------------------------------------------------------
domain = distribute(domain)

#-----------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
from math import sin, pi, exp
Br = anuga.Reflective_boundary(domain)      # Solid reflective wall
Bt = anuga.Transmissive_boundary(domain)    # Continue all values on boundary
BdL = anuga.Dirichlet_boundary([3.58431872, q_s, 0.]) # Constant boundary values
BdR = anuga.Dirichlet_boundary([2.87870797, q_s, 0.]) # Constant boundary values

# Associate boundary tags with boundary objects
domain.set_boundary({'left': BdL, 'right': BdR, 'top': Br, 'bottom': Br})

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
for t in domain.evolve(yieldstep = 1.0, finaltime = 200.):
    #print(domain.timestepping_statistics(track_speeds=True))
    if myid == 0: print(domain.timestepping_statistics())
    #vis.update()


domain.sww_merge(delete_old=True)

finalize()
