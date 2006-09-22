"""Example of the inundationmodel.

A wave of water is washed up ontp a hypothetical beach.
This one uses the parallel api

To run:

mpirun -c 4 python beach.py
"""

######################
# Module imports 

from math import pi
import time

from anuga.shallow_water import Domain
from anuga.shallow_water import Reflective_boundary
from anuga.shallow_water import Dirichlet_boundary
from anuga.shallow_water import Time_boundary
from anuga.utilities.polygon import Polygon_function

from Numeric import choose, greater, ones, sin, exp

from anuga_parallel.parallel_api import myid, numprocs, distribute



#------------------
# Define geometries
#------------------

def bathymetry(x,y):
    cut = 75
    res = choose( greater(x, cut), ( -(x - 55)/5, -4*ones(x.shape) ))
    res == (100-y)/50 + 1
    return res
    

def topography(x, y):
    z = -4.0*x/25 + 8  + (100-y)/50                         # Beach
    z += 6*exp( -((x-30)/10)**2 ) * exp( -((y-50)/8)**2 )   # Mound
    z += 4*exp( -((x-30)/4)**2 ) * exp( -((y-26)/10)**2 )   # Extra ridge    
    z += 4*exp( -((x-30)/5)**2 ) * exp( -((y-10)/8)**2 )    # Extra ridge
    z -= 4*exp( -((x-15)/6)**2 ) * exp( -((y-20)/12)**2 )   # Depression
    z += 1.2*exp( -((x-88)/7)**2 ) + exp( -((y-20)/25)**2 ) # Seafloor
    return z


def riverbed(x, y):
    return (y-100)/70 - 4.0*x/25 + 8
    

# Polygons
shoreline = [[40,0], [100,0], [100,100], [65,100], 
             [55,90], [55,70], [56,50], [50,25], [40,0]]        

land = [[65,100], [55,90], [55,70], [56,50], [50,25], [40,0],
        [0,0], [0,100]]

water = [[55,0], [100,0], [100,100], [55,100]]
all = [[0,0], [0,100], [100,100], [100,0]]

	     
building1 = [[45,80], [49,78], [52,83], [46,83]]	     
building2 = [[35,75], [40,75], [40,80], [35,80]]	     
building3 = [[42,63], [46,61], [48,65], [44,67]]	     
building4 = [[28,56], [28,51], [33,51], [33,56]]	     
building5 = [[10,70], [10,65], [15,65], [15,70]]	     
building6 = [[10,50], [10,45], [15,45], [15,50]]	     

river = [[20,100], [18,90], [20,80], [20,60], [15,40], [11,20], [2,0], [10,0],
         [14,10], [20,30], [24,45], [27,80], [27,85], [35,90], [39,100]]
	      


#----------------------
# Domain
#----------------------
name = 'beach'
print 'Creating domain from %s.tsh' %name
domain = Domain(mesh_filename = name + '.tsh',
                use_cache=True, verbose=True)
domain.set_name(name)
domain.set_default_order(2)


#----------------------
# Initial conditions
#----------------------
domain.set_quantity('elevation',
                    Polygon_function( [(all, topography),
                                       (building1, 7), (building2, 8), 
                                       (building3, 7), (building4, 13),
                                       (building5, 10), (building6, 11)]))
		      
domain.set_quantity('stage', 
		    Polygon_function( [(water, -1.5), 
                                       (land, -10)] )) 

domain.set_quantity('friction', 0.03)
print domain.statistics()


#----------------------
# Boundary conditions
#----------------------
Br = Reflective_boundary(domain)
Bd = Dirichlet_boundary([-12, 0.0, 0.0])
#Bt = Time_boundary(domain, lambda t: [ 3.0*(1+sin(2*pi*t/100)), 0.0, 0.0])

tags = {}
tags['wall'] = Br 
tags['wall1'] = Br 
tags['outflow'] = Bd 
tags['exterior'] = Br
tags['external'] = Br
tags['land'] = Bd  
tags['westbank'] = None    #Riverbank
tags['eastbank'] = None    #Riverbank
tags['eastbankN'] = None   #Riverbank
tags['ocean'] = None       # Bind this one later  

domain.set_boundary(tags)


#--------------------
# Distribute domain
#--------------------
domain = distribute(domain)

Bt = Time_boundary(domain, lambda t: [ 4.0*(1+sin(2*pi*t/50)), -1.0, 0.0])
domain.modify_boundary({'ocean': Bt})

#----------------------
# Evolve through time
#----------------------
t0 = time.time()
for t in domain.evolve(yieldstep = 0.2, finaltime = 300):
    domain.write_time()
    
print 'Simulation took %.2f seconds' %(time.time()-t0)
    
