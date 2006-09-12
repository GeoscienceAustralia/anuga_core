"""Example of shallow water wave equation.

This version is for profiling
(The original)


FIXME: This should really measure something else, such as profiling the set up of domains.
"""

######################
# Module imports 
#
from shallow_water import Domain, Reflective_boundary, Dirichlet_boundary,\
     Transmissive_boundary, Time_boundary,\
     Weir_simple as Weir, Constant_height

from mesh_factory import rectangular
from Numeric import array
    

######################
# Domain
#

N = 128

print 'Creating domain'
#Create basic mesh
points, vertices, boundary = rectangular(N, N)

#Create shallow water domain
domain = Domain(points, vertices, boundary)
domain.default_order = 2
domain.set_quantities_to_be_stored(None)
#domain.visualise = True


print 'Field values'
def slope(x, y):
    return -x
    
domain.set_quantity('elevation', slope)
domain.set_quantity('friction', 0.3)


# Boundary conditions
#
Br = Reflective_boundary(domain)
domain.set_boundary({'left': Br, 'right': Br, 'bottom': Br, 'top': Br})


######################
#Initial condition
domain.set_quantity('stage', Constant_height(slope, 0.3))

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
