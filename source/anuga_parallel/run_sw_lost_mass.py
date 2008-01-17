#!/usr/bin/env/python
#########################################################
#
# Shows problem of conservation of mass. Due to material being
# removed if too shallow. Can lead to problems if timestep
# set very small.
#
#
#
#  Authors: Linda Stals, Steve Roberts and Matthew Hardy,
# June 2005
#
#
#
#########################################################


##############################################
# Set a small yieldstep (say 0.005) to see
# lose of conservation
# Set it to a larger value (say 0.1) and have
# conservation "restored"
##############################################
yieldstep = 0.005
#yieldstep = 0.08
#yieldstep = 0.05
finaltime = 100.0

import sys
import time


from Numeric import array, zeros, Float

#from shallow_water import Domain

from anuga.shallow_water import Domain

# mesh partition routines

from anuga.abstract_2d_finite_volumes.pmesh2domain\
     import pmesh_to_domain_instance

# read in the processor information


#-------
# Domain
rect = zeros( 4, Float) # Buffer for results

class Set_Stage:
    """Set an initial condition with constant water height, for x0<x<x1
    """

    def __init__(self, x0=0.25, x1=0.5, h=1.0):
        self.x0 = x0
        self.x1 = x1
        self.h  = h

    def __call__(self, x, y):
        return self.h*((x>self.x0)&(x<self.x1))



    # read in the test files

filename = 'test-100.tsh'

nx = 1
ny = 1

domain = pmesh_to_domain_instance(filename, Domain)

domain.set_quantity('stage', Set_Stage(200.0,300.0,5.0))
rect = array(domain.xy_extent, Float)

try:
    domain.initialise_visualiser(rect=rect)
    #domain.visualiser.coloring['stage'] = True
    domain.visualiser.scale_z['stage'] = 0.2
    domain.visualiser.scale_z['elevation'] = 0.05
except:
    print 'No visualiser'


domain.default_order = 2

#Boundaries
from shallow_water import Transmissive_boundary, Reflective_boundary

T = Transmissive_boundary(domain)
R = Reflective_boundary(domain)
domain.set_boundary( {'outflow': R, 'inflow': R, 'inner':R, 'exterior': R, 'open':R} )



#---------
# Evolution
t0 = time.time()

print 'No of triangles %d'%(domain.number_of_triangles)


for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
    #domain.write_time()
    print '    Integral stage = ', domain.quantities['stage'].get_integral(),' Time = ',domain.time


