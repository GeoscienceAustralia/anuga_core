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

#---------
#Setup computational domain
#---------
domain = anuga.rectangular_cross_domain(40,40, len1=1., len2=1.)


#------------------------------------------------------------------------------
# Setup Algorithm, either using command line arguments
# or override manually yourself
#------------------------------------------------------------------------------
from anuga.utilities.argparsing import parse_standard_args
alg, cfl = parse_standard_args()
domain.set_flow_algorithm(alg)
domain.set_CFL(cfl)
domain.set_name('runup_sinusoid')                         # Output to file runup.sww
domain.set_datadir('.')                          # Use current folder
domain.set_quantities_to_be_stored({'stage': 2, 'xmomentum': 2, 'ymomentum': 2, 'elevation': 1})
domain.set_minimum_allowed_height(0.01)


#------------------
# Define topography
#------------------
scale_me=1.0
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

#------------------------------
#Evolve the system through time
#------------------------------
for t in domain.evolve(yieldstep=1.0,finaltime=40.0):
    print domain.timestepping_statistics()
    xx = domain.quantities['xmomentum'].centroid_values
    yy = domain.quantities['ymomentum'].centroid_values
    dd = domain.quantities['stage'].centroid_values - domain.quantities['elevation'].centroid_values 
    dd = (dd)*(dd>1.0e-03)+1.0e-06
    vv = ( (xx/dd)**2 + (yy/dd)**2 )**0.5
    vv = vv*(dd>1.0e-03)
    print 'Peak velocity is: ', vv.max(), vv.argmax()
print 'Finished'
