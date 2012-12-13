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


#-------------------------------------------------------------------------------
# Copy scripts to time stamped output directory and capture screen
# output to file
#-------------------------------------------------------------------------------
time = strftime('%Y%m%d_%H%M%S',localtime())

#output_dir = 'steep_island'+time
output_dir = '.'
output_file = 'steep_island'

#anuga.copy_code_files(output_dir,__file__)
#start_screen_catcher(output_dir+'_')


#------------------------------------------------------------------------------
# Setup domain
#------------------------------------------------------------------------------
dx = 1.
dy = dx
L = 2000.
W = 5*dx

# structured mesh
points, vertices, boundary = anuga.rectangular_cross(int(L/dx), int(W/dy), L, W, (0.0, 0.0))

#domain = anuga.Domain(points, vertices, boundary) 
domain = Domain(points, vertices, boundary) 

domain.set_name(output_file)                
domain.set_datadir(output_dir) 

#------------------------------------------------------------------------------
# Setup Algorithm, either using command line arguments
# or override manually yourself
#------------------------------------------------------------------------------
from anuga.utilities.argparsing import parse_standard_args
alg, cfl = parse_standard_args()
domain.set_flow_algorithm(alg)
domain.set_CFL(cfl)

#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------
domain.set_quantity('friction', 0.0)

def stage_flat(x,y):
    w=zeros(len(x))
    for i in range(len(x)):
        w[i]=4.5
    return w
domain.set_quantity('stage', stage_flat)

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
domain.set_quantity('elevation', bed_elevation)

#-----------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
from math import sin, pi, exp
Br = anuga.Reflective_boundary(domain)      # Solid reflective wall
#Bt = anuga.Transmissive_boundary(domain)    # Continue all values on boundary 
#Bd = anuga.Dirichlet_boundary([1,0.,0.]) # Constant boundary values

# Associate boundary tags with boundary objects
domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})


#===============================================================================
##from anuga.visualiser import RealtimeVisualiser
##vis = RealtimeVisualiser(domain)
##vis.render_quantity_height("stage", zScale =h0*500, dynamic=True)
##vis.colour_height_quantity('stage', (0.0, 0.5, 1.0))
##vis.start()
#===============================================================================


#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep = 0.1, finaltime = 5.):
    #print domain.timestepping_statistics(track_speeds=True)
    print domain.timestepping_statistics()
    #vis.update()


#test against know data
    
#vis.evolveFinished()

