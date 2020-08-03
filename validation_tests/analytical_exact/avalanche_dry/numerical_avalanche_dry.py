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
from numpy import zeros
from time import localtime, strftime, gmtime
from numpy import sin, cos, tan, arctan
from anuga import myid, finalize, distribute

#-------------------------------------------------------------------------------
# Copy scripts to time stamped output directory and capture screen
# output to file
#-------------------------------------------------------------------------------
time = strftime('%Y%m%d_%H%M%S',localtime())

#output_dir = 'avalanche_'+time
output_dir = '.'
output_file = 'avalanche'

#anuga.copy_code_files(output_dir,__file__)
#start_screen_catcher(output_dir+'_')


#------------------------------------------------------------------------------
# Setup domain
#------------------------------------------------------------------------------
g = 9.81           # gravity
h_0 = 20.0         # depth upstream. Note that the depth downstream is 0.0
L = 200.0          # length of stream/domain
dx = 0.5
dy = dx
W = 3*dx

def stage(X,Y):
    N = len(X)
    w = zeros(N)
    for i in range(N):
        if X[i]<=0.0:
            w[i] = bed_slope*X[i]
        else:
            w[i] = bed_slope*X[i] + h_0
    return w

def elevation(X,Y):
    N = len(X)
    y=zeros(N)
    for i in range(N):
        y[i] = bed_slope*X[i]
    return y


class Coulomb_friction:
    
    def __init__(self,
                 friction_slope=0.05,
                 bed_slope=0.1):
        
        self.friction_slope = friction_slope   #tan(delta) # NOTE THAT friction_slope must less than bed_slope
        self.bed_slope = bed_slope             #tan(theta) # bottom slope, positive if it is increasing bottom.

        thet = arctan(bed_slope)
        self.F = g*cos(thet)*cos(thet)*friction_slope
        self.m = -1.0*g*bed_slope + self.F
        
    def __call__(self, domain):        

        w = domain.quantities['stage'].centroid_values
        z = domain.quantities['elevation'].centroid_values
        h = w-z

        xmom_update = domain.quantities['xmomentum'].explicit_update
        ymom_update = domain.quantities['ymomentum'].explicit_update

        xmom_update[:] = xmom_update + self.F*h


bed_slope = 0.1
friction_slope = 0.05
Coulomb_forcing_term = Coulomb_friction(friction_slope, bed_slope)

def f_right(t):
    z_r = bed_slope*(0.5*L)
    h_r = h_0 #+ bed_slope*cell_len
    w_r = z_r + h_r
    u_r = Coulomb_forcing_term.m*t
    #['stage', 'xmomentum', 'ymomentum']
    return [w_r,  u_r*h_r,  0.0]

#------------------------------------------------------------------------------
# Setup Algorithm, either using command line arguments
# or override manually yourself
#------------------------------------------------------------------------------
args = anuga.get_args()
alg = args.alg
verbose = args.verbose

if myid == 0:
    
    # structured mesh
    points, vertices, boundary = anuga.rectangular_cross(int(L/dx), int(W/dy), L, W, (-L/2.0, -W/2.0))

    #domain = anuga.Domain(points, vertices, boundary) 
    domain = Domain(points, vertices, boundary) 

    domain.set_name(output_file)                
    domain.set_datadir(output_dir)
    
    domain.set_flow_algorithm(alg)
    
    #------------------------------------------------------------------------------
    # Setup initial conditions
    #------------------------------------------------------------------------------

    # No Mannings friction, but we introduce Coulomb friction below.
    domain.set_quantity('friction', 0.0)
    domain.set_quantity('stage', stage)
    domain.set_quantity('elevation',elevation)

else:

    domain = None


domain = distribute(domain)


#-----------------------------------------------------------------------------
# Special implementation of Coulomb friction
#-----------------------------------------------------------------------------  
domain.forcing_terms.append(Coulomb_forcing_term)


#-----------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
from math import sin, pi, exp
Br = anuga.Reflective_boundary(domain)      # Solid reflective wall
Bt = anuga.Transmissive_boundary(domain)    # Continue all values on boundary 
#Bd = anuga.Dirichlet_boundary([1,0.,0.]) # Constant boundary values
BTime = anuga.Time_boundary(domain,f_right)

# Associate boundary tags with boundary objects
domain.set_boundary({'left': Bt, 'right': BTime, 'top': Br, 'bottom': Br})



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
for t in domain.evolve(yieldstep = 0.1, finaltime = 3.):
    if myid == 0 and verbose: print(domain.timestepping_statistics())
    

domain.sww_merge(delete_old=True)


finalize()


