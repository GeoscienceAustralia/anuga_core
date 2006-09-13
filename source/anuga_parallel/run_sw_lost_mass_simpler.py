#!/usr/bin/env python
#########################################################
#
# Shows problem of conservation of mass. Due to material being
# removed if too shallow. Can lead to problems if timestep
# set very small.
#
#
#
#  Authors: Linda Stals, Steve Roberts, Matthew Hardy, Ole Nielsen
#  July 2005
#
#
#########################################################


##############################################
# Change min_depth and see how larger values aggravate the problem
yieldstep = 0.1
finaltime = 10.0
#min_depth = 1.0e-4
min_depth = 1.0e-2


from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular
from anuga.shallow_water import Domain, Reflective_boundary



#Create shallow water domain
points, vertices, boundary = rectangular(50, 50, len1=500, len2=500)
domain = Domain(points, vertices, boundary)
domain.smooth = False
domain.visualise = True
domain.default_order = 1
domain.minimum_allowed_height = min_depth
print 'Extent', domain.get_extent()

# Set initial condition
class Set_IC:
    """Set an initial condition with a constant value, for x0<x<x1
    """

    def __init__(self, x0=0.25, x1=0.5, h=1.0):
        self.x0 = x0
        self.x1 = x1
        self.h  = h

    def __call__(self, x, y):
        return self.h*((x>self.x0)&(x<self.x1))


domain.set_quantity('stage', Set_IC(200.0,300.0,5.0))

try:
    domain.initialise_visualiser()
except:
    print 'No visualiser'
else:    
    domain.visualiser.scale_z['stage'] = 0.2
    domain.visualiser.scale_z['elevation'] = 0.05
    

#Boundaries
R = Reflective_boundary(domain)
domain.set_boundary( {'left': R, 'right': R, 'top':R, 'bottom': R} )


# Evolution
print 'Minimal allowed water height = ', domain.minimum_allowed_height 
for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
    print 'Integral stage = ', domain.quantities['stage'].get_integral(),' Time = ',domain.time


