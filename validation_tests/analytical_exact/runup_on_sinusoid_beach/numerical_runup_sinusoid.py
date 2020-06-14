"""
Runup example from the manual, slightly modified
"""
#---------
#Import Modules
#--------
import anuga
import numpy
from math import sin, pi, exp
from anuga import Domain
from anuga import myid, finalize, distribute

args = anuga.get_args()
alg = args.alg
verbose = args.verbose

scale_me=1.0

#-------------------------
# create Sequential domain
#-------------------------
if myid == 0:
	#---------
	#Setup computational domain
	#---------
	domain = anuga.rectangular_cross_domain(40,40, len1=1., len2=1.)
	

	domain.set_flow_algorithm(alg)
	domain.set_name('runup_sinusoid')                         # Output to file runup.sww
	domain.set_datadir('.')                          # Use current folder
	
	#------------------
	# Define topography
	#------------------
	def topography(x,y):
		return (-x/2.0 +0.05*numpy.sin((x+y)*50.0))*scale_me
	
	def stagefun(x,y):
		stge=-0.2*scale_me
		return stge
	
	domain.set_quantity('elevation',topography)     # Use function for elevation
	domain.get_quantity('elevation').smooth_vertex_values() 
	domain.set_quantity('friction',0.0)             # Constant friction
	domain.set_quantity('stage', stagefun)          # Constant negative initial stage
	domain.get_quantity('stage').smooth_vertex_values()

else:
	domain = None
	
#------------------------
# Create parallel domain
#------------------------
domain = distribute(domain)
	

#--------------------------
# Setup boundary conditions
#--------------------------
Br=anuga.Reflective_boundary(domain)            # Solid reflective wall
Bt=anuga.Transmissive_boundary(domain)          # Continue all values of boundary -- not used in this example
Bd=anuga.Dirichlet_boundary([-0.1*scale_me,0.,0.])       # Constant boundary values -- not used in this example

#----------------------------------------------
# Associate boundary tags with boundary objects
#----------------------------------------------
domain.set_boundary({'left': Br, 'right': Bd, 'top': Br, 'bottom':Br})


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
for t in domain.evolve(yieldstep=1.0, finaltime=200.0):
	
	if myid == 0 and verbose: print(domain.timestepping_statistics())
# 	xx = domain.quantities['xmomentum'].centroid_values
# 	yy = domain.quantities['ymomentum'].centroid_values
# 	dd = domain.quantities['stage'].centroid_values - domain.quantities['elevation'].centroid_values 
# 	dd = (dd)*(dd>1.0e-03)+1.0e-06
# 	vv = ( (xx/dd)**2 + (yy/dd)**2 )**0.5
# 	vv = vv*(dd>1.0e-03)
# 	print('Peak velocity is: ', vv.max(), vv.argmax())

domain.sww_merge(delete_old=True)

finalize()
