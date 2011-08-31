"""Simple water flow example using ANUGA

Oscillatory wave from off shore coming into shore
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

import sys
import anuga
from swb_domain import *

from math import cos
import numpy as num
from time import localtime, strftime, gmtime
from os import sep



#-------------------------------------------------------------------------------
# Copy scripts to time stamped output directory and capture screen
# output to file
#-------------------------------------------------------------------------------
time = strftime('%Y%m%d_%H%M%S',localtime())

output_dir = '.'
output_file = 'data_wave_shelf_'+time

#copy_code_files(output_dir,__file__)
#start_screen_catcher(output_dir+sep)

interactive_visualisation = True

#------------------------------------------------------------------------------
# Setup domain
#------------------------------------------------------------------------------
dx = 1000.
dy = dx
L = 150000.
W = 10*dx

# structured mesh
points, vertices, boundary = anuga.rectangular_cross(int(L/dx), int(W/dy), L, W, (0.0, -W/2))

domain = Domain(points, vertices, boundary)

domain.set_name(output_file)                
domain.set_datadir(output_dir)  

#------------------------------------------------------------------------------
# Setup Algorithm
#------------------------------------------------------------------------------
domain.set_timestepping_method('rk2')
domain.set_default_order(2)

print domain.get_timestepping_method()

#domain.use_edge_limiter = True
#domain.tight_slope_limiters = False
#domain.use_centroid_velocities = False


#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------
depth1 = -100.0
depth2 = -50.00
def topography(x,y):
    z = num.zeros_like(x)
    z[:] = depth1

    N = len(x)
    for i in range(N):
        if  40000.0 < x[i] < 60000.0:
            z[i] = depth1 + (depth2-depth1)*((x[i] - 40000.0)/20000.0)

        if x[i] > 60000.0:
            z[i] = depth2

        if x[i] > 100000.0:
            z[i] = depth2*(1.0 - (x[i] - 75000.0)/50000)
    return z

domain.set_quantity('elevation',topography)
domain.set_quantity('friction', 0.00)
domain.set_quantity('stage', 0.0)            

#-----------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
from math import sin, pi, exp
Br = anuga.Reflective_boundary(domain)      # Solid reflective wall
Bt = anuga.Transmissive_boundary(domain)    # Continue all values on boundary
Bd = anuga.Dirichlet_boundary([1,0.,0.]) # Constant boundary values
amplitude = 1
wave_length = 2000.0
Bw = anuga.Time_boundary(domain=domain,     # Time dependent boundary
## Sine wave
                  f=lambda t: [(amplitude*sin((1./wave_length)*t*2*pi)), 0.0, 0.0])
## Sawtooth?
#                   f=lambda t: [(-8.0*(sin((1./180.)*t*2*pi))+(1./2.)*sin((2./180.)*t*2*pi)+(1./3.)*sin((3./180.)*t*2*pi)), 0.0, 0.0])
## Sharp rise, linear fall
#                   f=lambda t: [(5.0*(-((t-0.)/300.)*(t<300.)-cos((t-300.)*2.*pi*(1./240.))*(t>=300. and t<420.)+(1.-(t-420.)/300.)*(t>=420. and t <720.))), 0.0, 0.0])
#                   f=lambda t: [amplitude*(1.-2.*(pi*(1./720.)*(t-720.))**2)/exp((pi*(1./720.)*(t-720.))**2) , 0.0, 0.0])
#                   f=lambda t: [(-8.0*sin((1./720.)*t*2*pi))*((t<720.)-0.5*(t<360.)), 0.0, 0.0])

# Associate boundary tags with boundary objects
domain.set_boundary({'left': Bw, 'right': Br, 'top': Br, 'bottom': Br})


#===============================================================================
if interactive_visualisation:
    from anuga.visualiser import RealtimeVisualiser
    vis = RealtimeVisualiser(domain)
    vis.render_quantity_height("stage", zScale =10000, dynamic=True)
    vis.colour_height_quantity('stage', (1.0, 0.5, 0.5))
    vis.start()
#===============================================================================


#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------

min = 60.0
hr = 60*min

for t in domain.evolve(yieldstep = min, finaltime = 5*hr):
    domain.write_time()
    if interactive_visualisation:
        vis.update()

    if domain.get_time()> 2*hr:
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

if interactive_visualisation:
    vis.evolveFinished()

