#!/usr/bin/env python
"""
  Example showing improved wet dry limiting.
  This script runs the same simulation twice and stores the outputs in
  old_limiter.sww and new_limiter.sww respectively.
  
  Authors: Steve Roberts
  October
"""

#-----------------------------------------------------------------
# Common structures
#-----------------------------------------------------------------
import time
from Numeric import array
from anuga.shallow_water import Domain
from anuga.shallow_water import Reflective_boundary
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular

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

M = 30
points, vertices, boundary = rectangular(M, M, len1 = 1.0, len2 = 1.0)

yieldstep = 0.002
finaltime = 0.06
rect = [0.0, 0.0, 1.0, 1.0]


#-----------------------------------------------------------------
# Create domain for "old limiter" scenario
#-----------------------------------------------------------------
domain = Domain(points, vertices, boundary)
domain.set_name('old_limiter_second_order')
print 'Number of triangles =', len(domain)

# Turn on the visualisation
try:
    domain.initialise_visualiser()
except:
    pass


#-----------------------------------------------------------------
# Boundaries and Initial conditions
#-----------------------------------------------------------------
R = Reflective_boundary(domain)
domain.set_boundary( {'left': R, 'right': R, 'bottom': R, 'top': R} )
domain.set_quantity('stage', Set_Stage(0.2, 0.4, 0.25, 0.75, 1.0, 0.00))

#-----------------------------------------------------------------
# Values for old limiter
#-----------------------------------------------------------------
domain.default_order = 2
domain.beta_w      = 0.9
domain.beta_w_dry  = 0.9
domain.beta_uh     = 0.9
domain.beta_uh_dry = 0.9
domain.beta_vh     = 0.9
domain.beta_vh_dry = 0.9

#-----------------------------------------------------------------
# Evolve
#-----------------------------------------------------------------
t0 = time.time()
for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
    domain.write_time()

print 'That took %.2f seconds' %(time.time()-t0)
print 'Note the small timesteps and the irregular flow'
#raw_input('press return to continue')


#-----------------------------------------------------------------
# Create domain for "new limiter" scenario (2 order)
#-----------------------------------------------------------------
domain = Domain(points, vertices, boundary)
domain.set_name('new_limiter_second_order')
print 'Number of triangles =', len(domain)

# Turn on the visualisation
try:
    domain.initialise_visualiser()
except:
    pass

#-----------------------------------------------------------------
# Boundaries and Initial conditions
#-----------------------------------------------------------------
R = Reflective_boundary(domain)
domain.set_boundary( {'left': R, 'right': R, 'bottom': R, 'top': R} )
domain.set_quantity('stage', Set_Stage(0.2, 0.4, 0.25, 0.75, 1.0, 0.00))

#-----------------------------------------------------------------
# Values for new limiter
#-----------------------------------------------------------------
domain.set_default_order(2)
domain.minimum_allowed_height = 0.001
domain.beta_w      = 1.0
domain.beta_w_dry  = 0.2
domain.beta_uh     = 1.0
domain.beta_uh_dry = 0.2
domain.beta_vh     = 1.0
domain.beta_vh_dry = 0.2

#-----------------------------------------------------------------
# Evolve
#-----------------------------------------------------------------
t0 = time.time()
for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
    domain.write_time()

print 'That took %.2f seconds' %(time.time()-t0)
print 'Note more uniform and large timesteps'
#raw_input('press return to continue')


#-----------------------------------------------------------------
# Create domain for "new limiter" scenario (1 order)
#-----------------------------------------------------------------
domain = Domain(points, vertices, boundary)
domain.set_name('new_limiter_first_order')
print 'Number of triangles =', len(domain)

# Turn on the visualisation
try:
    domain.initialise_visualiser()
except:
    pass

#-----------------------------------------------------------------
# Boundaries and Initial conditions
#-----------------------------------------------------------------
R = Reflective_boundary(domain)
domain.set_boundary( {'left': R, 'right': R, 'bottom': R, 'top': R} )
domain.set_quantity('stage', Set_Stage(0.2, 0.4, 0.25, 0.75, 1.0, 0.00))

#-----------------------------------------------------------------
# Values for new limiter first order
#-----------------------------------------------------------------
domain.set_default_order(1)
domain.minimum_allowed_height = 0.001
domain.beta_w      = 1.0
domain.beta_w_dry  = 0.2
domain.beta_uh     = 1.0
domain.beta_uh_dry = 0.2
domain.beta_vh     = 1.0
domain.beta_vh_dry = 0.2

#-----------------------------------------------------------------
# Evolve
#-----------------------------------------------------------------
t0 = time.time()
for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
    domain.write_time()

print 'That took %.2f seconds' %(time.time()-t0)


