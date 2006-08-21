"""

Test example to demonstrate the "Black Screen of Death"

* Simple mesh set up on square with square internal polygon
* Set elevation using geospatial data
* Duplicate point in internal polygon (or bounding polygon) can
  cause the "Black Screen of Death" only when the elevation is set
  using the geospatial data.

Jane Sexton 2006

"""

######################
# Module imports 
#
from anuga.pmesh.mesh_interface import create_mesh_from_regions
from anuga.pyvolution.pmesh2domain import pmesh_to_domain_instance
from caching import cache
from anuga.pyvolution.shallow_water import Domain, Reflective_boundary,\
     Dirichlet_boundary, Time_boundary, Transmissive_boundary
from anuga.geospatial_data.geospatial_data import *


# bounding polygon
c0 = [331000, 6255000]
c1 = [345000, 6255000]
c2 = [345000, 6267000]
c3 = [331000, 6267000]
bound_poly = [c0, c1, c2, c3]

# set up internal region
p0 = [335000, 6260000]
p1 = [340000, 6260000]
p2 = [340000, 6265000]
p3 = [335000, 6265000]
interior_polygon = [p0, p1, p2, p3]
# duplicate point, this should fail, i.e. conjugate gradient doesn't converge
interior_polygon = [p0, p1, p2, p3, p0] 

interior_region_area = 2500000000 #25000
interior_regions = [[interior_polygon, interior_region_area]]

meshname = 'bedslope_bsod.tsh'

exterior_region_area = 10000000000000000 #100000
# create a mesh
_ = cache(create_mesh_from_regions,
          bound_poly,
          {'boundary_tags': {'bottom': [0], 'right': [1],
                             'top': [2], 'left': [3]},
           'maximum_triangle_area': exterior_region_area,
           'filename': meshname,           
           'interior_regions': interior_regions},
          evaluate = True,
          verbose = True)

#Create shallow water domain
domain = pmesh_to_domain_instance(meshname, Domain,
                                  use_cache = True,
                                  verbose = True)

domain.set_name('bedslope_bsod')
domain.set_datadir('.')                      #Use current directory for output
domain.set_quantities_to_be_stored('stage')  #See shallow_water.py

# set up geospatial data object
test_pts = []
test_elev = []
a = [340000, 6264000]
b = [340000, 6266000]
c = [335000, 6259000]
d = [335000, 6266000]
e = [335000, 6258000]
f = [335000, 6259000]
g = [341000, 6260000]
h = [341000, 6256000]
test_pts = [a, b, c, d, e, f, g, h]
test_elev = [1.0, 4.0, 3.0, 0.1, 5, -100.0, -200, -15]
G = Geospatial_data(test_pts, test_elev)

# set geospatial data object to elevation
domain.set_quantity('elevation', G)
# if set elevation as a constant and duplicate a point in the interior region
# the black screen of death will NOT be seen
#domain.set_quantity('elevation', 1)
domain.set_quantity('friction', 0.0)
domain.set_quantity('stage', 0.0)


######################
# Boundary conditions
from math import pi, sin
Br = Reflective_boundary(domain)
Bt = Transmissive_boundary(domain)
Bd = Dirichlet_boundary([0.2,0.,0.])
Bw = Time_boundary(domain=domain,
                   f=lambda t: [(-1000*sin(t*2*pi)), 0.0, 0.0])

print 'Tags are ', domain.get_boundary_tags()

domain.set_boundary({'left': Bd, 'right': Bw, 'top': Br, 'bottom': Br})

######################
#Evolution
domain.check_integrity()

for t in domain.evolve(yieldstep = 5, finaltime = 50.0):
    domain.write_time()
