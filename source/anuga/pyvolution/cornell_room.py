"""Example of shallow water wave equation.

This is called Netherlands because it shows a dam with a gap in it and
stylised housed behind it and below the water surface.

"""

######################
# Module imports 
#
from shallow_water import Domain, Reflective_boundary, Dirichlet_boundary,\
     Transmissive_boundary, Time_boundary, Constant_height

from mesh_factory import from_polyfile, rectangular
from Numeric import array
from math import sqrt
#from least_squares import Interpolation
from fit_interpolate.interpolate import Interpolate
    

print 'Creating domain'
#data_points, _, data_values = from_polyfile('cornell_room_medres')
#points, triangles, values = from_polyfile('hires2')
data_points, _, data_values = from_polyfile('hires2')


#Regrid onto numerically stable mesh
#
#Compute regular mesh based on resolution and extent of data
data_points = array(data_points)
pmax = max(data_points)
pmin = min(data_points)

M = len(data_points)

N = int(0.8*sqrt(M))

#print N

mesh_points, vertices, boundary = rectangular(N, N, 
					      len1=pmax[0]-pmin[0],
                                              len2=pmax[1]-pmin[1],
					      origin = pmin)
					      

#Compute smooth surface on new mesh based on values from old (regrid)
print 'Interp'
interp = Interpolate(mesh_points, vertices, alpha=10)
mesh_values = interp.fit( data_points, data_values) # this has not been tested
print 'Len mesh values', len(mesh_values)
print 'Len mesh points', len(mesh_points)

    
#Create shallow water domain
print 'Creating domain'
domain = Domain(mesh_points, vertices) #, boundary)

domain.check_integrity()
domain.default_order = 2
domain.smooth = True
domain.reduction = min  #Looks a lot better on top of steep slopes

print "Number of triangles = ", len(domain)

domain.visualise = False
domain.checkpoint = False
domain.store = True    #Store for visualisation purposes
domain.format = 'sww'   #Native netcdf visualisation format
import sys, os
root, ext = os.path.splitext(sys.argv[0])
if domain.smooth is True:
    s = 'smooth'
else:
    s = 'nonsmooth'        
domain.filename = root + '_' + s

#Set bed-slope and friction
manning = 0.0

print 'Field values'
domain.set_quantity('elevation', mesh_values)
domain.set_quantity('friction', manning)


######################
# Boundary conditions
#
print 'Boundaries'
Br = Reflective_boundary(domain)
domain.set_boundary({'exterior': Br})

                  

######################
#Initial condition
#
print 'Initial condition'

#Define water height as a lump in one corner
def height(x, y):
    from Numeric import zeros, Float
    
    N = len(x)
    assert N == len(y)    

    xmin = min(x); xmax = max(x)
    ymin = min(y); ymax = max(y)

    xrange = xmax - xmin
    yrange = ymax - ymin    

    z = zeros(N, Float)
    for i in range(N):
        if x[i] <= xmin + 0.25*xrange and y[i] <= ymin + 0.25*yrange:
            z[i] = 300

    return z

domain.set_quantity('stage', height)

E = domain.quantities['elevation'].vertex_values
L = domain.quantities['stage'].vertex_values
domain.set_quantity('stage', E+L)

#Evolve
for t in domain.evolve(yieldstep = 0.05, finaltime = 5.0):
    domain.write_time()
    

    
