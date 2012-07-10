"""Runup example from the manual, slightly modified
"""
#---------
#Import Modules
#--------
import anuga

import numpy

from math import sin, pi, exp
#from anuga.shallow_water_balanced2.swb2_domain import Domain as Domain
from anuga.shallow_water.shallow_water_domain import Domain as Domain
#from shallow_water_balanced_steve.swb_domain import *
#import shallow_water_balanced_steve.swb_domain 
#from shallow_water_balanced_steve.swb_domain import Domain as Domain
#path.append('/home/gareth/storage/anuga_clean/anuga_jan12/trunk/anuga_work/shallow_water_balanced_steve')
#from swb_domain import *
#path.append('/home/gareth/storage/anuga_clean/anuga_jan12/trunk/anuga_work/development/gareth/balanced_basic')
#from balanced_basic import *
#from balanced_dev import *
#---------
#Setup computational domain
#---------
points, vertices, boundary = anuga.rectangular_cross(40,40)

domain=Domain(points,vertices,boundary)    # Create Domain
domain.set_name('runup_v2')                         # Output to file runup.sww
domain.set_datadir('.')                          # Use current folder
domain.set_quantities_to_be_stored({'stage': 2, 'xmomentum': 2, 'ymomentum': 2, 'elevation': 1})
#domain.set_store_vertices_uniquely(True)
#------------------
# Define topography
#------------------

def topography(x,y):
	return -x/2 #+0.05*numpy.sin((x+y)*50.0) #+0.1*(numpy.random.rand(len(x)) -0.5)       # Linear bed slope + small random perturbation 

def stagefun(x,y):
    #stg=-0.2*(x<0.5) -0.1*(x>=0.5)
    stg=-0.2 # Stage
    #topo=topography(x,y) #Bed
    return stg #*(stg>topo) + topo*(stg<=topo)

domain.set_quantity('elevation',topography)     # Use function for elevation
domain.get_quantity('elevation').smooth_vertex_values() # Steve's fix -- without this, substantial artificial velcities are generated everywhere in the domain. With this fix, there are artificial velocities near the coast, but not elsewhere.

domain.set_quantity('friction',0.00)             # Constant friction

domain.set_quantity('stage', stagefun)              # Constant negative initial stage

#--------------------------
# Setup boundary conditions
#--------------------------
Br=anuga.Reflective_boundary(domain)            # Solid reflective wall
Bt=anuga.Transmissive_boundary(domain)          # Continue all values of boundary -- not used in this example
Bd=anuga.Dirichlet_boundary([-0.2,0.,0.])       # Constant boundary values -- not used in this example
Bw=anuga.Time_boundary(domain=domain,
	f=lambda t: [(0.0*sin(t*2*pi)-0.1)*exp(-t)-0.1,0.0,0.0]) # Time varying boundary -- get rid of the 0.0 to do a runup.

#----------------------------------------------
# Associate boundary tags with boundary objects
#----------------------------------------------
domain.set_boundary({'left': Br, 'right': Bw, 'top': Br, 'bottom':Br})

#------------------------------
#Evolve the system through time
#------------------------------

#xwrite=open("xvel.out","wb")
#ywrite=open("yvel.out","wb")
## Set print options to be compatible with file writing via the 'print' statement 
#numpy.set_printoptions(threshold=numpy.nan, linewidth=numpy.nan)

for t in domain.evolve(yieldstep=0.2,finaltime=30.0):
    print domain.timestepping_statistics()


print 'Finished'
