"""
Numerical solution to debris avalanche involving a shock wave
Ref: Mungkasi and Roberts, Pure and Applied Geophysics 2012
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
h_0 = 20.0         # depth upstream
h_1 = 10.0         # depth downstream
L = 200.0          # length of stream/domain
dx = 0.5
dy = dx
W = 3*dx

args=anuga.get_args()
alg = args.alg
verbose = args.verbose

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
        #ymom_update = domain.quantities['ymomentum'].explicit_update

        xmom_update[:] = xmom_update + self.F*h


bed_slope = 0.1
friction_slope = 0.05
Coulomb_forcing_term = Coulomb_friction(friction_slope, bed_slope)

def stage(X,Y):
    N = len(X)
    w = zeros(N)
    for i in range(N):
        if X[i]<=0.0:
            w[i] = bed_slope*X[i] + h_1
        else:
            w[i] = bed_slope*X[i] + h_0
    return w

def elevation(X,Y):
    N = len(X)
    y=zeros(N)
    for i in range(N):
        y[i] = bed_slope*X[i]
    return y

def f_right(t):
    z_r = bed_slope*(0.5*L)
    h_r = h_0 #+ bed_slope*cell_len
    w_r = z_r + h_r
    u_r = Coulomb_forcing_term.m*t
    #['stage', 'xmomentum', 'ymomentum']
    return [w_r,  u_r*h_r,  0.0]

def f_left(t):
    z_l = bed_slope*(-0.5*L)
    h_l = h_1 #+ bed_slope*cell_len
    w_l = z_l + h_l
    u_l = Coulomb_forcing_term.m*t
    #['stage', 'xmomentum', 'ymomentum']
    return [w_l,  u_l*h_l,  0.0]



#===============================================================================
# Create sequential domain
#===============================================================================
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
    
    domain.forcing_terms.append(Coulomb_forcing_term)
    
    domain.set_quantity('stage', stage)
    
    domain.set_quantity('elevation',elevation)
    
else:
    
    domain = None

#===============================================================================
# Create parallel domain
#===============================================================================
domain = distribute(domain)

#-----------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
from math import sin, pi, exp
Br = anuga.Reflective_boundary(domain)      # Solid reflective wall
#Bt = anuga.Transmissive_boundary(domain)    # Continue all values on boundary 
#Bd = anuga.Dirichlet_boundary([1,0.,0.]) # Constant boundary values
BTimeR = anuga.Time_boundary(domain,f_right)
BTimeL = anuga.Time_boundary(domain,f_left)

# Associate boundary tags with boundary objects
domain.set_boundary({'left': BTimeL, 'right': BTimeR, 'top': Br, 'bottom': Br})





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
for t in domain.evolve(yieldstep = 0.1, finaltime = 4.):
    #print domain.timestepping_statistics(track_speeds=True)
    if myid == 0 and verbose: print domain.timestepping_statistics()
    #vis.update()


domain.sww_merge(delete_old=True)


finalize()

