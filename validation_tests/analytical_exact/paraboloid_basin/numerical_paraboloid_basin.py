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
from anuga import myid, finalize, distribute

#-------------------------------------------------------------------------------
# Copy scripts to time stamped output directory and capture screen
# output to file
#-------------------------------------------------------------------------------
#timer = strftime('%Y%m%d_%H%M%S',localtime())
output_dir = '.'  #'paraboloid_'+timer
output_file = 'paraboloid'
#anuga.copy_code_files(output_dir,__file__)
#start_screen_catcher(output_dir+'_')


args = anuga.get_args()
alg = args.alg
verbose = args.verbose

#-------------------------------------------------------------------------------
# Domain
#-------------------------------------------------------------------------------
n = 75#50
m = 75#50
lenx = 8000.0
leny = 8000.0
origin = (-4000.0, -4000.0)



#===============================================================================
# Create Sequential Domain
#===============================================================================
if myid == 0:
    points, elements, boundary = anuga.rectangular_cross(m, n, lenx, leny, origin)
    domain = anuga.Domain(points, elements, boundary)


    domain.set_name(output_file)                
    domain.set_datadir(output_dir)
    
    domain.set_flow_algorithm(alg)

    
    #-------------------------------------------------------------------------------
    # Initial conditions
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

else:
    
    domain = None
 
#============================================================================
# Create Parallel Domain
#============================================================================
domain = distribute(domain)


#---------
# Boundary
#R = anuga.Reflective_boundary(domain)
#T = anuga.Transmissive_boundary(domain)
D = anuga.Dirichlet_boundary([0.0, 0.0, 0.0])
domain.set_boundary({'left': D, 'right': D, 'top': D, 'bottom': D})


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

import time
t0 = time.time()
for t in domain.evolve(yieldstep = 1.0, finaltime = 200.0):
    if myid == 0 and verbose : domain.write_time()

if myid==0 and verbose: print('That took %s secs' % str(time.time()- t0))

domain.sww_merge(delete_old=True)

finalize()
