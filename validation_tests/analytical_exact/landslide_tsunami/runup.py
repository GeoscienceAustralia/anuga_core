"""Runup example from the manual, slightly modified
"""
#---------
#Import Modules
#--------
import anuga
from anuga import myid, finalize, distribute

import numpy
import scipy.interpolate

from math import sin, pi, exp


args = anuga.get_args()
alg = args.alg
verbose = args.verbose

#timer = strftime('%Y%m%d_%H%M%S',localtime())
#---------
#Setup computational domain
# --- Very inefficient mesh!
#---------

dx = dy = 25.
l1 = 60.*1000.
l2 = dy
nx =int(l1/dx)
ny =int(l2/dy)


# Beach slope of 1/10
def topography(x,y):
	return -(x-200.)/10. 

#--------------------------------------------------------------------------------
# Create Sequential Domain
#--------------------------------------------------------------------------------
if myid == 0:
	print(' Building mesh (alternative non-uniform mesh could be much more efficient)')
	points, vertices, boundary = anuga.rectangular_cross(nx,ny, len1=l1,len2=l2, origin=(-200., 0.))
	
	print('Creating Domain')
	domain=anuga.Domain(points,vertices,boundary)    # Create Domain
	domain.set_name('runup_v2')                         # Output to file runup.sww
	domain.set_datadir('.')                          # Use current folder
	domain.set_quantities_to_be_stored({'stage': 2, 'xmomentum': 2, 'ymomentum': 2, 'elevation': 1})
	

	domain.set_flow_algorithm(alg)

	#------------------
	# Define Initial conditions
	#------------------
	
	
	# Define initial condition using file
	initial_prof=scipy.genfromtxt('./DATA/initial_condition.txt', skip_header=13)
	initial_prof_fun=scipy.interpolate.interp1d(initial_prof[:,0], initial_prof[:,1], fill_value=0., bounds_error=False)
	
	def stagefun(x,y):
	    return initial_prof_fun(x-200.) 
	
	domain.set_quantity('elevation',topography)     # Use function for elevation	
	domain.set_quantity('friction',0.000)             # Constant friction
	domain.set_quantity('stage', stagefun)              

else:
	
	domain = None

#------------------------------------------------------------------------
# Parallel Domain
#------------------------------------------------------------------------	
domain = distribute(domain)
	
#--------------------------
# Setup boundary conditions
#--------------------------
Br=anuga.Reflective_boundary(domain)            # Solid reflective wall

#----------------------------------------------
# Associate boundary tags with boundary objects
#----------------------------------------------
domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom':Br})

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

for t in domain.evolve(yieldstep=5.0,finaltime=350.0):
	if myid == 0: print(domain.timestepping_statistics())
    #uh=domain.quantities['xmomentum'].centroid_values
    #vh=domain.quantities['ymomentum'].centroid_values
    #depth=domain.quantities['stage'].centroid_values - domain.quantities['elevation'].centroid_values
    #depth=depth*(depth>1.0e-06) + 1.0e-06
    #vel=((uh/depth)**2 + (vh/depth)**2)**0.5
    #print 'peak speed is', vel.max()


domain.sww_merge(delete_old=True)

finalize()
