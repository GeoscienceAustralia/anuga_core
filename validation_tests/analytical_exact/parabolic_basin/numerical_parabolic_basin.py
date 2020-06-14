"""Parabolic channel oscilation example
"""
#---------
#Import Modules
#--------
import anuga
from anuga import myid, finalize, distribute
import numpy
from math import sqrt, cos, sin, pi
from numpy import asarray
from time import localtime, strftime, gmtime
from anuga.geometry.polygon import inside_polygon, is_inside_triangle

#-------------------------------------------------------------------------------
# Copy scripts to time stamped output directory and capture screen
# output to file
#-------------------------------------------------------------------------------
timer = strftime('%Y%m%d_%H%M%S',localtime())
output_dir = '.'  #'parabola_'+timer
output_file = 'parabola'

#anuga.copy_code_files(output_dir,__file__)
#start_screen_catcher(output_dir+'_')

args = anuga.get_args()
alg = args.alg
verbose = args.verbose

m = 200
n = 10
lenx = 40.0
leny = 2.0
origin = (-20.0, -1.0)


#-------------------------------------------------------------------------------
# Sequential Domain
#-------------------------------------------------------------------------------
if myid == 0:
    points, elements, boundary = anuga.rectangular_cross(m, n, lenx, leny, origin)
    domain = anuga.Domain(points, elements, boundary)
    domain.set_name(output_file)                
    domain.set_datadir(output_dir)
    #domain.set_minimum_allowed_height(0.01)
    
    
    
    domain.set_flow_algorithm(alg)
    
    
    #------------------
    # Define topography
    #------------------
    
    # Parameters for analytical solution
    D0=4.0
    L=10.0
    A = 2.0
    
    def topography(x,y):
            return  (D0/(L**2.))*x**2.
    
    def stage_init(x,y):
        wat= D0 + (2.0*A*D0/L**2.)*(x-A/2.0) # Water elevation inside the parabola
        top=topography(x,y) # Bed elevation
        # Return the maximum of the water elevation and the bed elvation
        return wat*(wat>top) + top*(wat<=top)
    
    
    domain.set_quantity('elevation',topography)     # Use function for elevation
    domain.set_quantity('friction',0.0)            # No friction
    domain.set_quantity('stage', stage_init)        # Constant negative initial stage
    
else:
    
    domain = None
    
#--------------------------------------------------------------------
# Parallel Domain
#--------------------------------------------------------------------
domain = distribute(domain)
    
#--------------------------
# Setup boundary conditions
#--------------------------
Br=anuga.Reflective_boundary(domain)            # Solid reflective wall
#Bt=anuga.Transmissive_boundary(domain)          # Continue all values of boundary -- not used in this example
#Bd=anuga.Dirichlet_boundary([-0.2,0.,0.])       # Constant boundary values -- not used in this example
#Bw=anuga.Time_boundary(domain=domain,
#	f=lambda t: [(0.0*sin(t*2*pi)-0.1)*exp(-t)-0.1,0.0,0.0]) # Time varying boundary -- get rid of the 0.0 to do a runup.

#----------------------------------------------
# Associate boundary tags with boundary objects
#----------------------------------------------
domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom':Br})


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
if myid == 0:
    parameter_file=open('parameters.tex', 'w')
    parameter_file.write('\\begin{verbatim}\n')
    from pprint import pprint
    pprint(domain.get_algorithm_parameters(),parameter_file,indent=4)
    parameter_file.write('\\end{verbatim}\n')
    parameter_file.close()


#------------------------------
#Evolve the system through time
#------------------------------
for t in domain.evolve(yieldstep=0.05,finaltime=10.0):
    if myid == 0: print(domain.timestepping_statistics())

domain.sww_merge(delete_old=True)

finalize()




