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
from anuga.geometry.polygon import inside_polygon, is_inside_triangle


args = anuga.get_args()
alg = args.alg
verbose = args.verbose

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
dx = 0.01
dy = dx
L = 1.6
W = 0.61

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
domain.set_flow_algorithm(alg)


#------------------------------------------------------------------------------
# Setup initial conditions
#------------------------------------------------------------------------------
def stage(x,y):
    h = zeros(len(x), float)
    for i in range(len(x)):
        if x[i] <= 0.4:
            h[i] = 0.3
        elif x[i] <= 0.4+0.5:
            h[i] = 0.01
        elif x[i] <= 0.4+0.5+0.12:
            if 0.25 <= y[i] <= 0.25+0.12:
                h[i] = 0.75
            else:
                h[i] = 0.01
        else:
            h[i] = 0.01
    return h

def elevation(x,y):
    z = zeros(len(x), float)
    for i in range(len(x)):
        if 0.4+0.5 <= x[i] <= 0.4+0.5+0.12:
            if 0.25 <= y[i] <= 0.25+0.12:
                z[i] = 0.75
    return z
domain.set_quantity('stage', stage)
domain.set_quantity('elevation',elevation)
domain.set_quantity('friction', 0.03)

#-----------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
from math import sin, pi, exp
Br = anuga.Reflective_boundary(domain)      # Solid reflective wall
Bt = anuga.Transmissive_boundary(domain)    # Continue all values on boundary 
Bd = anuga.Dirichlet_boundary([1,0.,0.]) # Constant boundary values

# Associate boundary tags with boundary objects
domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})


#===============================================================================
##from anuga.visualiser import RealtimeVisualiser
##vis = RealtimeVisualiser(domain)
##vis.render_quantity_height("stage", zScale =h0*500, dynamic=True)
##vis.colour_height_quantity('stage', (0.0, 0.5, 1.0))
##vis.start()
#===============================================================================


#---------------------------------------------
# Find triangle that contains the point Point
# and print to file
#---------------------------------------------
##PointL = []
##print "0.12/dx=",int(round(0.12/dx))
##STOP
##PointL=[(0.868, 0.255),
##          (0.868, 0.265),
##          (0.868, 0.275),
##          (0.868, 0.285),
##          (0.868, 0.295),
##          (0.868, 0.305),
##          (0.868, 0.315),
##          (0.868, 0.325),
##          (0.868, 0.335),
##          (0.868, 0.345),
##          (0.868, 0.355),
##          (0.868, 0.365)]
###PointR = []
##PointR=[(1.032, 0.255),
##          (1.032, 0.265),
##          (1.032, 0.275),
##          (1.032, 0.285),
##          (1.032, 0.295),
##          (1.032, 0.305),
##          (1.032, 0.315),
##          (1.032, 0.325),
##          (1.032, 0.335),
##          (1.032, 0.345),
##          (1.032, 0.355),
##          (1.032, 0.365)]
##for i in range(len(PointR)):
##    Point = PointR[i]
##    for n in range(len(domain.triangles)):
##        tri = domain.get_vertex_coordinates(n)
##        if is_inside_triangle(Point,tri):
##            #print 'Point is within triangle with vertices '+'%s'%tri
##            n_point = n
##            break
##    print 'n_point = ',n_point


#------------------------------------------------------------------------------
# Produce a documentation of parameters
#------------------------------------------------------------------------------
parameter_file=open('parameters.tex', 'w')
parameter_file.write('\\begin{verbatim}\n')
from pprint import pprint
pprint(domain.get_algorithm_parameters(),parameter_file,indent=4)
parameter_file.write('\\end{verbatim}\n')
parameter_file.close()

#------------------------------------------------------------------------------
# Evolve system through time
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep = 0.05, finaltime = 3.0):
    #print domain.timestepping_statistics(track_speeds=True)
    print domain.timestepping_statistics()
    #vis.update()


#test against know data
    
#vis.evolveFinished()

