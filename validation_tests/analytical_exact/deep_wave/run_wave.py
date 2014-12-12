"""Simple water flow example using ANUGA

Will Powers example of a simple sinusoidal wave. 
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

import sys
import anuga
from anuga import Domain
from anuga import myid, finalize, distribute



from math import cos
import numpy as num
from time import localtime, strftime, gmtime
from os import sep


#--------------------------------
# Get Default values for basic
# algorithm parameters.
#--------------------------------
args = anuga.get_args()
alg = args.alg
verbose = args.verbose

#-------------------------------------------------------------------------------
# Copy scripts to time stamped output directory and capture screen
# output to file
#-------------------------------------------------------------------------------
time = strftime('%Y%m%d_%H%M%S',localtime())

output_dir = '.'
output_file = 'data_wave'

interactive_visualisation = False

if myid == 0:
    #------------------------------------------------------------------------------
    # Setup sequential domain
    #------------------------------------------------------------------------------
    dx = 500.
    dy = dx
    L = 100000.
    W = 10*dx
    
    # structured mesh
    points, vertices, boundary = anuga.rectangular_cross(int(L/dx), int(W/dy), L, W, (0.0, -W/2))
    
    domain = Domain(points, vertices, boundary) 
    
    domain.set_name(output_file)                
    domain.set_datadir(output_dir)  
    
    #------------------------------------------------------------------------------
    # Setup Algorithm, either using command line arguments
    # or override manually yourself
    #------------------------------------------------------------------------------
    domain.set_flow_algorithm(alg)
    
    #------------------------------------------------------------------------------
    # Setup initial conditions
    #------------------------------------------------------------------------------
    domain.set_quantity('elevation',-100.0)
    domain.set_quantity('friction', 0.00)
    domain.set_quantity('stage', 0.0)            

else: 

    domain = None
    
#-----------------------------------------------------------------------------------
# Parallel Domain 
#-----------------------------------------------------------------------------------
domain = distribute(domain)
   
#-----------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
from math import sin, pi, exp
Br = anuga.Reflective_boundary(domain)      # Solid reflective wall
#Bt = anuga.Transmissive_boundary(domain)    # Continue all values on boundary

amplitude = 1
wave_length = 300.0

# Incoming wave

def waveform2(t):
    return amplitude*sin((1./wave_length)*t*2*pi)

Bw2 = anuga.Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(domain, waveform2)


# A radiation type condition on the other side
def zero_fun(t):
    return 0.0                    

Bw5 = anuga.Flather_external_stage_zero_velocity_boundary(domain, zero_fun)

# Associate boundary tags with boundary objects
domain.set_boundary({'left': Bw2, 'right': Bw5, 'top': Br, 'bottom': Br})

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

for t in domain.evolve(yieldstep = 10, finaltime = 2e4):
    domain.write_time()
    if interactive_visualisation:
        vis.update()

domain.sww_merge(delete_old=True)

finalize()

