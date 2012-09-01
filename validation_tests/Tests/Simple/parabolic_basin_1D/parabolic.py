"""Parabolic channel oscilation example
"""
#---------
#Import Modules
#--------
import anuga
import numpy

#---------
#Setup computational domain
#---------
points, vertices, boundary = anuga.rectangular_cross(200,10, len1=40.0,len2=2.0)

domain=anuga.Domain(points,vertices,boundary)    # Create Domain
domain.set_name('parabola_v2')                         # Output to file runup.sww
domain.set_datadir('.')                          # Use current folder
#domain.set_minimum_allowed_height(0.01)


#------------------------------------------------------------------------------
# Setup Algorithm, either using command line arguments
# or override manually yourself
#------------------------------------------------------------------------------
from anuga.utilities.argparsing import parse_standard_args
alg, cfl = parse_standard_args()
domain.set_flow_algorithm(alg)
domain.set_CFL(cfl)

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
for t in domain.evolve(yieldstep=1.0,finaltime=90.0):
    print domain.timestepping_statistics()

print 'Finished'
