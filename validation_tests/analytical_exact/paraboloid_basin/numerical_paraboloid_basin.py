"""Example of shallow water wave equation analytical solution
consists of a parabolic profile in a parabolic basin. Analytical
solution to this problem was derived by Thacker
and used by Yoon and Chou.

   Copyright 2005
   Christopher Zoppou, Stephen Roberts, ANU, Geoscience Australia
   Modified by Sudi Mungkasi, ANU, 2011
"""

#-------------------------------------------------------------------------------
# Module imports
#-------------------------------------------------------------------------------
import sys
import anuga
from anuga import g
from math import sqrt, cos, sin, pi
from numpy import asarray
from time import localtime, strftime, gmtime
from anuga.geometry.polygon import inside_polygon, is_inside_triangle


#-------------------------------------------------------------------------------
# Copy scripts to time stamped output directory and capture screen
# output to file
#-------------------------------------------------------------------------------
#timer = strftime('%Y%m%d_%H%M%S',localtime())
output_dir = '.'  #'paraboloid_'+timer
output_file = 'paraboloid'
#anuga.copy_code_files(output_dir,__file__)
#start_screen_catcher(output_dir+'_')


#-------------------------------------------------------------------------------
# Domain
#-------------------------------------------------------------------------------
n = 75#50
m = 75#50
lenx = 8000.0
leny = 8000.0
origin = (-4000.0, -4000.0)

points, elements, boundary = anuga.rectangular_cross(m, n, lenx, leny, origin)
domain = anuga.Domain(points, elements, boundary)


#-------------------------------------------------------------------------------
# Provide file name for storing output
#-------------------------------------------------------------------------------
domain.store = True #False
domain.format = 'sww'
domain.set_name(output_file)                
domain.set_datadir(output_dir)


#-------------------------------------------------------------------------------
# Order of scheme
# Good compromise between
# limiting and CFL
# CUT THIS FOR AUTOMATED TESTS
#-------------------------------------------------------------------------------
#domain.set_default_order(2)
#domain.set_timestepping_method(2)
#domain.set_beta(1.5)
##domain.set_CFL(0.6)
#domain.smooth = True


#-------------------------------------------------------------------------------
# Decide which quantities are to be stored at each timestep
#-------------------------------------------------------------------------------
#domain.quantities_to_be_stored = ['stage', 'xmomentum', 'ymomentum']


#-------------------------------------------------------------------------------
# Reduction operation for get_vertex_values
from anuga.utilities.numerical_tools import mean
domain.reduction = mean #domain.reduction = min  #Looks better near steep slopes


#-------------------------------------------------------------------------------
# Initial condition
#-------------------------------------------------------------------------------
t = 0.0
D0 = 1000.
L = 2500.
R0 = 2000.

A = (L**4 - R0**4)/(L**4 + R0**4)
omega = 2./L*sqrt(2.*g*D0)
T = pi/omega

#------------------
# Set bed elevation
def bed_elevation(x,y):
    n = x.shape[0]
    z = 0*x
    for i in range(n):
        r = sqrt(x[i]*x[i] + y[i]*y[i])
        z[i] = -D0*(1.-r*r/L/L)
    return z
domain.set_quantity('elevation', bed_elevation)

#----------------------------
# Set the initial water level
def stage_init(x,y):
    z = bed_elevation(x,y)
    n = x.shape[0]
    w = 0*x
    for i in range(n):
        r = sqrt(x[i]*x[i] + y[i]*y[i])
        w[i] = D0*((sqrt(1-A*A))/(1.-A*cos(omega*t))
                -1.-r*r/L/L*((1.-A*A)/((1.-A*cos(omega*t))**2)-1.))
    if w[i] < z[i]:
        w[i] = z[i]
    return w
domain.set_quantity('stage', stage_init)

#---------
# Boundary
#R = anuga.Reflective_boundary(domain)
#T = anuga.Transmissive_boundary(domain)
D = anuga.Dirichlet_boundary([0.0, 0.0, 0.0])
domain.set_boundary({'left': D, 'right': D, 'top': D, 'bottom': D})

#---------------------------------------------
# Find triangle that contains the point Point
# and print to file
#---------------------------------------------
Point = (0.0, 0.0)
for n in range(len(domain.triangles)):
    tri = domain.get_vertex_coordinates(n)
    if is_inside_triangle(Point,tri):
        #print 'Point is within triangle with vertices '+'%s'%tri
        n_point = n
        break
print 'n_point = ',n_point
#bla
    
#t = domain.triangles[n_point]
#print 't = ', t
#tri = domain.get_vertex_coordinates(n)

filename=domain.get_name()
file = open(filename,'w')

#----------
# Evolution
import time
t0 = time.time()

time_array = []
stage_array = []
Stage     = domain.quantities['stage']
Xmomentum = domain.quantities['xmomentum']
Ymomentum = domain.quantities['ymomentum']

#interactive_visualisation = True
#===============================================================================
##if interactive_visualisation:
##    from anuga.visualiser import RealtimeVisualiser
##    vis = RealtimeVisualiser(domain)
##    vis.render_quantity_height("stage", zScale =1000, dynamic=True)
##    vis.colour_height_quantity('stage', (1.0, 0.5, 0.5))
##    vis.start()
#===============================================================================


#------------------------------------------------------------------------------
# Produce a documentation of parameters
#------------------------------------------------------------------------------
parameter_file=open('parameters.tex', 'w')
parameter_file.write('\\begin{verbatim}\n')
from pprint import pprint
pprint(domain.get_algorithm_parameters(),parameter_file,indent=4)
parameter_file.write('\\end{verbatim}\n')
parameter_file.close()

for t in domain.evolve(yieldstep = 1.0, finaltime = 200.0):
    domain.write_time()

    #tri_array = asarray(tri)
    #t_array = asarray([[0,1,2]])
    #interp = Interpolation(tri_array,t_array,[points])

    stage     = Stage.get_values(location='centroids',indices=[n_point])
    xmomentum = Xmomentum.get_values(location='centroids',indices=[n_point])
    ymomentum = Ymomentum.get_values(location='centroids',indices=[n_point])
    #print '%10.6f   %10.6f  %10.6f   %10.6f\n'%(t,stage,xmomentum,ymomentum)
    file.write( '%10.6f   %10.6f  %10.6f   %10.6f\n'%(t,stage,xmomentum,ymomentum) )

    #time_array.append(t)
    #stage_array.append(stage)

    #if interactive_visualisation:
    #    vis.update()
"""
file.close()
print 'That took %.2f seconds' %(time.time()-t0)
from pylab import *
#ion()
hold(False)
plot(time_array, stage_array, 'r.-')
#title('Gauge %s' %name)
xlabel('time (s)')
ylabel('stage (m)')
#legend(('Observed', 'Modelled'), shadow=True, loc='upper left')
#savefig(name, dpi = 300)

#raw_input('Next')
show()
"""

#if interactive_visualisation:
#    vis.evolveFinished()
