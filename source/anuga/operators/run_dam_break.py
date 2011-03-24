"""Simple water flow example using ANUGA

Water driven up a linear slope and time varying boundary,
similar to a beach environment
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

import sys
import anuga


from math import cos
from numpy import zeros, float, where
from time import localtime, strftime, gmtime



#-------------------------------------------------------------------------------
# Copy scripts to time stamped output directory and capture screen
# output to file
#-------------------------------------------------------------------------------
#time = strftime('%Y%m%d_%H%M%S',localtime())

#output_dir = 'dam_break_'+time
output_file = 'dam_break'

#anuga.copy_code_files(output_dir,__file__)
#start_screen_catcher(output_dir+'_')


#------------------------------------------------------------------------------
# Setup domain
#------------------------------------------------------------------------------
dx = 1000.
dy = dx
L = 100000.
W = 10*dx

# structured mesh
points, vertices, boundary = anuga.rectangular_cross(int(L/dx), int(W/dy), L, W, (0.0, -W/2))

domain = anuga.Domain(points, vertices, boundary) 

domain.set_name(output_file)                
#domain.set_datadir(output_dir)

#------------------------------------------------------------------------------
# Setup Algorithm
#------------------------------------------------------------------------------
domain.set_timestepping_method('rk2')
domain.set_default_order(2)

#------------------------------------------------------------------------------
# Setup Kinematic Viscosity Operator
#------------------------------------------------------------------------------
domain.set_use_kinematic_viscosity(True)

domain.set_beta(1.7)
#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------
h0 = 10.0
h1 = 0.0

def elevation(x,y):

    return where(x<=95000.0, 0.0, h0)

def height(x,y):

    return where(x<=50000.0, h0, h1)

domain.set_quantity('elevation',elevation)
domain.set_quantity('friction', 0.0)
domain.add_quantity('stage', height)
#-----------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
from math import sin, pi, exp
Br = anuga.Reflective_boundary(domain)      # Solid reflective wall
Bt = anuga.Transmissive_boundary(domain)    # Continue all values on boundary 
Bd = anuga.Dirichlet_boundary([1,0.,0.]) # Constant boundary values


# Associate boundary tags with boundary objects
domain.set_boundary({'left': Br, 'right': Bt, 'top': Br, 'bottom': Br})


#===============================================================================
from anuga.visualiser import RealtimeVisualiser
vis = RealtimeVisualiser(domain)
vis.render_quantity_height("stage", zScale = 50000/(h0-h1), dynamic=True)
vis.colour_height_quantity('stage', (0.0, 0.5, 1.0))
vis.start()
#===============================================================================


#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------

for t in domain.evolve(yieldstep = 100.0, finaltime = 3*60*60.):
    #print domain.timestepping_statistics(track_speeds=True)
    domain.print_timestepping_statistics()
    domain.print_operator_timestepping_statistics()

#    if domain.get_time() > 2000.0 and domain.get_time() < 4000.0 :
#        domain.set_use_kinematic_viscosity(False)
#    else:
#        domain.set_use_kinematic_viscosity(True)
        
    vis.update()


#test against know data
    
vis.evolveFinished()

