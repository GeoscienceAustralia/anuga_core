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
from numpy import zeros, ones, array, interp, polyval, ones_like, zeros_like
from numpy import where, logical_and
from time import localtime, strftime, gmtime
from scipy.interpolate import interp1d
from anuga.geometry.polygon import inside_polygon, is_inside_triangle
#from balanced_dev import *


#-------------------------------------------------------------------------------
# Copy scripts to time stamped output directory and capture screen
# output to file
#-------------------------------------------------------------------------------
time = strftime('%Y%m%d_%H%M%S',localtime())

#output_dir = 'varying_width'+time
output_dir = '.'
output_file = 'varying_width'

#anuga.copy_code_files(output_dir,__file__)
#start_screen_catcher(output_dir+'_')


#------------------------------------------------------------------------------
# Setup domain
#------------------------------------------------------------------------------
dx = 1.
dy = dx
L = 1500.
W = 60.

# structured mesh
points, vertices, boundary = anuga.rectangular_cross(int(L/dx), int(W/dy), L, W, (0.,-W/2.))

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
domain.set_quantity('stage', 12.0)
XX = array([0.,50.,100.,150.,250.,300.,350.,400.,425.,435.,450.,470.,475.,500.,
            505.,530.,550.,565.,575.,600.,650.,700.,750.,800.,820.,900.,950.,
            1000.,1500.])
ZZ = array([0.,0.,2.5,5.,5.,3.,5.,5.,7.5,8.,9.,9.,9.,9.1,9.,9.,6.,5.5,5.5,5.,
            4.,3.,3.,2.3,2.,1.2,0.4,0.,0.])
WW = array([40.,40.,30.,30.,30.,30.,25.,25.,30.,35.,35.,40.,40.,40.,45.,45.,50.,
            45.,40.,40.,30.,40.,40.,5.,40.,35.,25.,40.,40.])/2.


depth = interp1d(XX, ZZ)
width = interp1d(XX, WW)

def bed_elevation(x,y):
    z = 25.0*ones_like(x)
    wid = width(x)
    dep = depth(x)
    z = where( logical_and(y < wid, y>-wid), dep, z)
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
##print 'n_point = ',n_point
##bla




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
for t in domain.evolve(yieldstep = 0.1, finaltime = 5.0):
    #print domain.timestepping_statistics(track_speeds=True)
    print domain.timestepping_statistics()
    #vis.update()


#test against know data
    
#vis.evolveFinished()

