"""Parabolic channel oscilation example
"""
#---------
#Import Modules
#--------
import anuga
from anuga import *
#from anuga.shallow_water.shallow_water_domain import Domain as Domain
import numpy

#from balanced_dev import *
#from balanced_basic import *
#from anuga.shallow_water_balanced2.swb2_domain import Domain as Domain
#from anuga.shallow_water.shallow_water_domain import Domain as Domain
#---------
#Setup computational domain
#---------
points, vertices, boundary = anuga.rectangular_cross(200,10, len1=40.0,len2=2.0)

domain=Domain(points,vertices,boundary)    # Create Domain
domain.set_name('parabola_v2')                         # Output to file runup.sww
domain.set_datadir('.')                          # Use current folder
domain.set_quantities_to_be_stored({'stage': 2, 'xmomentum': 2, 'ymomentum': 2, 'elevation': 1})
domain.set_minimum_allowed_height(0.01)
# Time stepping
#domain.set_timestepping_method('euler') # Default
#domain.set_timestepping_method('rk2') # 
#domain.beta_w= 1. #0.2
#domain.beta_uh= 1. #0.2
#domain.beta_vh= 1. #0.2
#domain.beta_w_dry= 0.2 #0.
#domain.beta_uh_dry= 0.2 #0.
#domain.beta_vh_dry= 0.2 #0.
#domain.alpha=100.

#------------------
# Define topography
#------------------

# Parameters for analytical solution
D0=4.0
L=10.0
A = 2.0

def topography(x,y):
	return  (D0/(L**2.))*(x-20.0)**2.

def stage_init(x,y):
    wat= ( D0 + (2.0*A*D0/L**2.)*((x-20.0)-A/2.0) ) # Water elevation inside the parabola
    top=topography(x,y) # Bed elevation
    # Return the maximum of the water elevation and the bed elvation
    return wat*(wat>top) + top*(wat<=top)


domain.set_quantity('elevation',topography)     # Use function for elevation

domain.set_quantity('friction',0.00)            # No friction

domain.set_quantity('stage', stage_init)        # Constant negative initial stage



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

#------------------------------
#Evolve the system through time
#------------------------------
for t in domain.evolve(yieldstep=0.1,finaltime=90.0):
    print domain.timestepping_statistics()

print 'Finished'
