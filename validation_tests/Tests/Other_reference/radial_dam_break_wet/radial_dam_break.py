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
from numpy import zeros, float
from time import localtime, strftime, gmtime
from anuga.geometry.polygon import inside_polygon, is_inside_triangle


#-------------------------------------------------------------------------------
# Copy scripts to time stamped output directory and capture screen
# output to file
#-------------------------------------------------------------------------------
time = strftime('%Y%m%d_%H%M%S',localtime())

output_dir = '.'#'radial_dam_break_'+time
output_file = 'radial_dam'

#anuga.copy_code_files(output_dir,__file__)
#start_screen_catcher(output_dir+'_')


#------------------------------------------------------------------------------
# Setup domain
#------------------------------------------------------------------------------
dx = 1.
dy = dx
L = 200.
W = L

# structured mesh
points, vertices, boundary = anuga.rectangular_cross(int(L/dx), int(W/dy), L, W, (-L/2.0, -W/2.0))
domain = anuga.Domain(points, vertices, boundary)
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
domain.set_quantity('elevation',0.0)
domain.set_quantity('friction', 0.0)

h0 = 10.0
h1 = 1.0
radius = 50.0

def height(x,y):
    z = zeros(len(x), float)
    r2 = radius**2
    for i in range(len(x)):
        rad2 = x[i]**2 + y[i]**2

        if rad2 <= r2:
            z[i] = h0
        else:
            z[i] = h1
    return z

domain.set_quantity('stage', height)

#-----------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
Br = anuga.Reflective_boundary(domain)      # Solid reflective wall

# Associate boundary tags with boundary objects
domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

#===============================================================================
##from anuga.visualiser import RealtimeVisualiser
##vis = RealtimeVisualiser(domain)
##vis.render_quantity_height("stage", zScale =h0, dynamic=True)
##vis.colour_height_quantity('stage', (0.0, 0.5, 1.0))
##vis.start()
#===============================================================================

#---------------------------------------------
# Find triangle that contains the point Point
# and print to file
#---------------------------------------------
##Point = (0.0, 0.0)
##for n in range(len(domain.triangles)):
##    tri = domain.get_vertex_coordinates(n)
##    if is_inside_triangle(Point,tri):
##        #print 'Point is within triangle with vertices '+'%s'%tri
##        n_point = n
##        break
##print 'The triangle ID containing the point of origin is = ',n_point


#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------

for t in domain.evolve(yieldstep = 0.1, finaltime = 2.0):
    #print domain.timestepping_statistics(track_speeds=True)
    print domain.timestepping_statistics()
    #vis.update()
##test against know data   
##vis.evolveFinished()

