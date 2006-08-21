"""Example of shallow water wave equation.

This version is for profiling of timestepping
"""

######################
# Module imports
from shallow_water import Domain, Reflective_boundary
from mesh_factory import rectangular
from Numeric import array


######################
# Domain

#N = 16   #Should take less than 0.1 s
#N = 64  #Should take less than 3s
N = 128  #Should take less than 20s


print 'Creating domain'
#Create basic mesh
points, vertices, boundary = rectangular(N, N)

#Create shallow water domain
domain = Domain(points, vertices, boundary)
domain.set_default_order(2)        
domain.set_quantities_to_be_stored(None)


print 'Setting initial conditions'

def slope(x, y):
    return -x

domain.set_quantity('elevation', slope)
domain.set_quantity('friction', 0.03)
domain.set_quantity('stage', expression = 'elevation + 0.5')



print 'Setting boundary conditions'
Br = Reflective_boundary(domain)
domain.set_boundary({'left': Br, 'right': Br, 'bottom': Br, 'top': Br})


print 'Checking integrity'
domain.check_integrity()


######################
#Evolution

import time
t0 = time.time()

s = 'for t in domain.evolve(yieldstep = 0.02, finaltime = 0.2): domain.write_time()'


import profile, pstats
FN = 'profile.dat'

profile.run(s, FN)

print 'That took %.2f seconds' %(time.time()-t0)

S = pstats.Stats(FN)
#S.sort_stats('time').print_stats(20)
s = S.sort_stats('cumulative').print_stats(30)

print s
