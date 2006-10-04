#!/usr/bin/env python
#########################################################
#
#  Example showing improved wet dry limiting.
#
#
#  Authors: Steve Roberts
# October
#
#
#
#########################################################


from Numeric import array
from anuga.shallow_water import Domain


# mesh partition routines
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular


M = 30
points, vertices, boundary = rectangular(M, M, len1 = 1.0, len2 = 1.0)
domain = Domain(points, vertices, boundary)
print 'number of triangles = ', domain.number_of_elements


#---------------------------
#Boundaries
#---------------------------
from anuga.shallow_water import Reflective_boundary

R = Reflective_boundary(domain)


domain.set_boundary( {'left': R, 'right': R, 'bottom': R, 'top': R} )
domain.check_integrity()

class Set_Stage:
    """Set an initial condition with constant water height, for x<x0
    """

    def __init__(self, x0=0.25, x1=0.75, y0=0.0, y1=1.0, h=5.0, h0=0.0):
        self.x0 = x0
        self.x1 = x1
        self.y0 = y0
        self.y1 = y1
        self.h  = h
        self.h0 = h0

    def __call__(self, x, y):
        return self.h0 + self.h*((x>self.x0)&(x<self.x1)&(y>self.y0)&(y<self.y1))

domain.set_quantity('stage', Set_Stage(0.2, 0.4, 0.25, 0.75, 1.0, 0.00))

import time
t0 = time.time()


# Turn on the visualisation

rect = [0.0, 0.0, 1.0, 1.0]
domain.initialise_visualiser()

yieldstep = 0.002
finaltime = 0.05


#===============================================================================
#Old Limiter
#===============================================================================

domain.default_order = 2
domain.beta_w      = 0.9
domain.beta_w_dry  = 0.9
domain.beta_uh     = 0.9
domain.beta_uh_dry = 0.9
domain.beta_vh     = 0.9
domain.beta_vh_dry = 0.9

#Check that the boundary value gets propagated to all elements
for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
    domain.write_time()


print 'That took %.2f seconds' %(time.time()-t0)
print 'Note the small timesteps and the irregular flow'
raw_input('press return to continue')

#===============================================================================
#New Limiter
#===============================================================================

t0 = time.time()
domain.set_quantity('stage', Set_Stage(0.2, 0.4, 0.25, 0.75, 1.0, 0.00))
domain.set_quantity('xmomentum', 0.0)
domain.set_quantity('ymomentum', 0.0)
domain.time = 0.0
domain.default_order = 2
domain.minimum_allowed_height = 0.001
domain.beta_w      = 1.0
domain.beta_w_dry  = 0.2
domain.beta_uh     = 1.0
domain.beta_uh_dry = 0.2
domain.beta_vh     = 1.0
domain.beta_vh_dry = 0.2

#Check that the boundary value gets propagated to all elements
for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
    domain.write_time()


print 'That took %.2f seconds' %(time.time()-t0)
print 'Note more uniform and large timesteps'
raw_input('press return to continue')

#===============================================================================
#First Order
#===============================================================================

t0 = time.time()
domain.set_quantity('stage', Set_Stage(0.2, 0.4, 0.25, 0.75, 1.0, 0.00))
domain.set_quantity('xmomentum', 0.0)
domain.set_quantity('ymomentum', 0.0)
domain.time = 0.0

domain.default_order = 1


#Check that the boundary value gets propagated to all elements
for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
    domain.write_time()


print 'That took %.2f seconds' %(time.time()-t0)

