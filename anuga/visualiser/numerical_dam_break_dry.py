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
dx = 0.25
dy = dx
L = 20.
W = 5.

# structured mesh
points, vertices, boundary = anuga.rectangular_cross(int(L/dx), int(W/dy), L, W, (-L/2.0, -W/2.0))

#domain = anuga.Domain(points, vertices, boundary) 
domain = Domain(points, vertices, boundary) 

domain.set_name(output_file)                
domain.set_datadir(output_dir) 

#------------------------------------------------------------------------------
# Setup Algorithm, either using command line arguments
# or override manually yourself
#------------------------------------------------------------------------------
from anuga.utilities.argparsing import parse_standard_args
alg = parse_standard_args()
domain.set_flow_algorithm(alg)


#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------
domain.set_quantity('elevation',0.0)
domain.set_quantity('friction', 0.0)

h0 = 5.0
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
vtk_visualiser = True
if vtk_visualiser:
    from anuga.visualiser import RealtimeVisualiser
    vis = RealtimeVisualiser(domain)
    vis.render_quantity_height("height",dynamic=True)
    #vis.render_quantity_height("stage", zScale =1.0, dynamic=True)
    vis.colour_height_quantity('stage', (0.0, 0.0, 1.0))
    vis.start()
#===============================================================================


#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
i = 0
for t in domain.evolve(yieldstep = 0.01, finaltime = 0.25):
    #print domain.timestepping_statistics(track_speeds=True)
    print(domain.timestepping_statistics())
    if vtk_visualiser: 
        vis.update()
        fileName = 'stage_%03d' % i + '.vtk'
        i = i+1
        vis.store_height_quantity('height',fileName)

#test against know data


    
if vtk_visualiser: vis.evolveFinished()

