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

#output_dir = 'dam_break_'+time
output_dir = '.'
output_file = 'dam_break'

#anuga.copy_code_files(output_dir,__file__)
#start_screen_catcher(output_dir+'_')


#------------------------------------------------------------------------------
# Setup domain
#------------------------------------------------------------------------------
dx = 1.
dy = dx
L = 1000.
W = 5*dx

# structured mesh
points, vertices, boundary = anuga.rectangular_cross(int(L/dx), int(W/dy), L, W, (-L/2.0, -W/2.0))

#domain = anuga.Domain(points, vertices, boundary) 
domain = Domain(points, vertices, boundary) 

domain.set_name(output_file)                
domain.set_datadir(output_dir)  

#------------------------------------------------------------------------------
# Setup Algorithm
#------------------------------------------------------------------------------
#domain.set_timestepping_method('rk2')
#domain.set_default_order(2)
#
#print domain.get_timestepping_method()
#
#domain.use_edge_limiter = True
##domain.use_edge_limiter = False
#domain.tight_slope_limiters = True
#domain.use_centroid_velocities = False
#
#domain.CFL = 1.0

#domain.beta_w      = 0.6
#domain.beta_uh     = 0.6
#domain.beta_vh     = 0.6


#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------
domain.set_quantity('elevation',0.0)
domain.set_quantity('friction', 0.0)

h0 = 10.0
h1 = 1.0

def height(x,y):
    z = zeros(len(x), float)
    for i in range(len(x)):
        if x[i]<=0.0:
            z[i] = h0
        else:
            z[i] = h1
    return z
domain.set_quantity('stage', height)

#-----------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
from math import sin, pi, exp
Br = anuga.Reflective_boundary(domain)      # Solid reflective wall
Bt = anuga.Transmissive_boundary(domain)    # Continue all values on boundary 
Bd = anuga.Dirichlet_boundary([1,0.,0.]) # Constant boundary values

# Associate boundary tags with boundary objects
domain.set_boundary({'left': Bt, 'right': Bt, 'top': Br, 'bottom': Br})


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
for t in domain.evolve(yieldstep = 0.5, finaltime = 50.):
    #print domain.timestepping_statistics(track_speeds=True)
    print domain.timestepping_statistics()
    #vis.update()


#test against know data
    
#vis.evolveFinished()

