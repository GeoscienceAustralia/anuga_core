"""Example of shallow water wave equation.

Two levels - testing that momentum is not genearted in the upward direction

"""

######################
# Module imports 
#
from os import sep, path
from mesh_factory import rectangular
from shallow_water import Domain, Reflective_boundary, Dirichlet_boundary,\
     Constant_height
from Numeric import array
from anuga.pyvolution.util import Polygon_function, read_polygon


#Create basic mesh
N = 2
points, vertices, boundary = rectangular(3, N, 20, 10)

#Create shallow water domain
domain = Domain(points, vertices, boundary)
domain.store = True
domain.set_name('two_levels')
print "Output being written to " + domain.get_datadir() + sep + \
	      domain.filename + "_size%d." %len(domain) + domain.format
	

domain.default_order=2
domain.smooth = False
domain.visualise = True

#PLAY WITH THIS [0;1]: 
#
# beta_h == 0.0 reveals the problem
# beta_h > 0.2 alleviates it
domain.beta_h = 0.2    

#IC
domain.set_quantity('friction', 0.0)
domain.set_quantity('stage', 3.0)

def twostage(x, y):  
    z = 0*x	 
    for i in range(len(x)):	 	
	if x[i]<10:	
    	    z[i] = 4.

    return z			
#Set elevation  
domain.set_quantity('elevation', twostage)


######################
# Boundary conditions
Br = Reflective_boundary(domain)

domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})
domain.check_integrity()


#print domain.quantities['elevation'].vertex_values    
######################
#Evolution
for t in domain.evolve(yieldstep = 0.01, finaltime = 0.05):
    domain.write_time()
    print domain.quantities['xmomentum'].edge_values    
    #print domain.quantities['stage'].centroid_values    

#print domain.quantities['stage'].vertex_values
#print
#print domain.quantities['xmomentum'].vertex_values
#print
#print domain.quantities['ymomentum'].vertex_values
    
print 'Done'    

