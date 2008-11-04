
from anuga.shallow_water import Domain, Reflective_boundary, Dirichlet_boundary 
from eqf_v2 import earthquake_tsunami

res = 10000000.
intres = 1000000.
boundingpoly=[[-100000,-80000],[-100000,80000],[100000,80000],[100000,-80000]]
poly = [[-40000,-20000],[40000,-20000],[40000,30000],[-40000,30000]]
interior_regions = [[poly,intres]]

from anuga.pmesh.mesh_interface import create_mesh_from_regions
from caching import cache
meshname = 'another_test'+'.msh'

_ = cache(create_mesh_from_regions,
          boundingpoly,
           {'boundary_tags': {'e0': [0], 'e1': [1], 'e2': [2], 'e3': [3]},
           'maximum_triangle_area': res,
           'filename': meshname,
           'interior_regions': interior_regions},
           verbose = True)

domain = Domain(meshname, use_cache = True, verbose = True)

basename = 'test'
domain.set_name(basename)
domain.set_quantities_to_be_stored(['stage', 'xmomentum', 'ymomentum'])
domain.set_minimum_storable_height(0.01)

tsunami_source = earthquake_tsunami(length=25.0,
                                    width=5.0,
                                    strike=0.,
                                    depth=3.50,
#                                    depth=3500.0,
                                    dip=15.0,
                                    slip=1.,
                                    x0=0.0,
                                    y0=0.0,
                                    domain=domain,
                                    verbose=True)

stage0 = 1.0

#sets tsunami_source to the domain
domain.set_quantity('stage', tsunami_source)
z2 = domain.get_quantity('stage')
int2 = z2.get_integral()/400000.0
#print 'integral before and after', int1, int2

print 'len(z2)',len(z2)
print 'dir(z2)',dir(z2)

domain.set_quantity('elevation', -10.0, alpha = 0.1)

Br = Reflective_boundary(domain)
Bd = Dirichlet_boundary([0,0,0])
domain.set_boundary( {'e0': Br,  'e1': Br, 'e2': Br, 'e3': Br} )

import time

t0 = time.time()

for t in domain.evolve(yieldstep = 1, finaltime = 1): 
    domain.write_time()
    domain.write_boundary_statistics(tags = 'e2')
    
from anuga.abstract_2d_finite_volumes.util import file_function
from Numeric import array, Float
from pylab import plot, ion, hold,savefig

'''
#points=array([[0,0],[0,10]])
max = 100
skip = 1000

y = 10000
points=array([[0]*2]*max)
print points
half_max_skip=(max*skip)/2
for i in range(max):
    print i*skip-half_max_skip
    points[i]=[i*skip-half_max_skip, y]
#    print i
'''
interval=500
profile_lenght= 20000
number_points = profile_lenght/interval
y = 10000
points=array([[0]*2]*number_points)
print points
half_profile=profile_lenght/2
for i in range(number_points):
    print i*interval-half_profile
    points[i]=[i*interval-half_profile, y]


print "points[0,:]",points
ion()
print 'hello'

F = file_function('test.sww', quantities = 'stage', interpolation_points=points,verbose = True,use_cache=True)

print F.statistics()
t=1
x = 0.0
y = 0.0
#print F(t=1, x=5000, y=5000)
profile=[]

for i in range(number_points):
    profile.append(F(0,point_id=i))
    print i, F(0,point_id=i)
print'profile', profile#,points[:,1]

plot(points[:,0],profile)
savefig("profile",dpi=300)

    
    


