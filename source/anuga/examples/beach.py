"""Example of the inundationmodel.

A wave of water is washed up ontp a hypothetical beach.

To run:

python beach.py
"""

######################
# Module imports 

import sys
from os import sep, path

from anuga.shallow_water import Domain, Reflective_boundary,\
     Dirichlet_boundary,\
     Transmissive_boundary, Time_boundary

from anuga.shallow_water.shallow_water_domain import Wind_stress

from anuga.utilities.polygon import read_polygon, Polygon_function
from math import pi
from Numeric import choose, greater, ones, sin, exp
import time

######################
# Domain
name = 'beach'
print 'Creating domain from %s.tsh' %name
domain = Domain(mesh_filename = name + '.tsh',
                use_cache=True, verbose=True)

domain.store = True
#domain.minimum_allowed_height = 0.0
domain.set_name(name + '6')
domain.default_order = 2
print "Output being written to " + domain.get_datadir() + sep + \
      domain.get_name() + domain.format


def bathymetry(x,y):
    cut = 75
    
    res = choose( greater(x, cut), ( -(x - 55)/5, -4*ones(x.shape) ))
    res == (100-y)/50 + 1
    #res += (100-y)/130  #Lift up southern bathymetry
    #res -= y/1500  #Pull down northern bathymetry
    
    return res
    

def topography(x, y):

    import RandomArray
    #z = 4*sin(pi*x/50) + (100-y)/50 + 1 #Ridge
    
    z = -4.0*x/25 + 8  + (100-y)/50 #Beach
    
    #z += sin(pi*x/5)/4 
    #z += RandomArray.normal(0, 1, x.shape)/4 
    
    z += 6*exp( -((x-30)/10)**2 ) * exp( -((y-50)/8)**2 )     #Mound
    z += 4*exp( -((x-30)/4)**2 ) * exp( -((y-26)/10)**2 )     #extra ridge    
    z += 4*exp( -((x-30)/5)**2 ) * exp( -((y-10)/8)**2 )     #extra ridge
    z -= 4*exp( -((x-15)/6)**2 ) * exp( -((y-20)/12)**2 )     #depression

    z += 1.2*exp( -((x-88)/7)**2 ) + exp( -((y-20)/25)**2 )     #Seafloor
    return z

def riverbed(x, y):
    return (y-100)/70 - 4.0*x/25 + 8
    

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
	      
print 'Set elevation'    
t0 = time.time()
domain.set_quantity('elevation',
    Polygon_function( [(all, topography),
                      #(shoreline, bathymetry), (land, topography), 
                      (building1, 7), (building2, 8), 
		      (building3, 7), (building4, 13),
		      (building5, 10), (building6, 11)]))
		      #(river, riverbed)]))
		      
print 'That took %.2f seconds' %(time.time()-t0)
		      
print 'Set stage'    		      
domain.set_quantity('stage', 
		    Polygon_function( [(water, -1.5), 
                                       (land, -10)] )) 
		    
print 'Set friction'    		      		    
domain.set_quantity('friction', 0.03)

print domain.get_extent()

print domain.statistics()


#Add lateral wind gusts bearing 135 degrees
def gust(t,x,y): 
    from math import sin, pi
    from Numeric import zeros, ones, Float

    N = len(x)

    tt = sin(2*pi*t/50)

    if tt > 0.98:
        return 100*tt*ones(N, Float)
    else:
        return zeros(N, Float)
    
#domain.forcing_terms.append(Wind_stress(gust, 135))


#Add lateral wind gusts bearing 90 degrees
def gust2(t,x,y): 
    from math import sin, pi
    from Numeric import zeros, ones, Float

    N = len(x)

    tt = sin(2*pi*t/100)

    if tt > 0.95:
        return 100*tt*ones(N, Float)
    else:
        return zeros(N, Float)
    
domain.forcing_terms.append(Wind_stress(gust2, 90))

#Add lateral wind gusts bearing 255 degrees
def gust3(t,x,y): 
    from math import sin, pi
    from Numeric import zeros, ones, Float

    N = len(x)

    tt = sin(2*pi*(t-30)/55)

    if tt > 0.96:
        return 24000*tt*ones(N, Float)
    else:
        return zeros(N, Float)
    
#domain.forcing_terms.append(Wind_stress(gust3, 255))


######################
# Boundary conditions

print 'Boundaries'
Br = Reflective_boundary(domain)
Bo = Transmissive_boundary(domain)

#Constant outflow
Bd = Dirichlet_boundary([-10, 0.0, 0.0])
Bt = Time_boundary(domain, lambda t: [ 3.0*(1+sin(2*pi*t/100)), 0.0, 0.0])

print 'Available boundary tags are', domain.get_boundary_tags()

#Set boundary conditions
tags = {}
tags['ocean'] = Bt
tags['wall'] = Br 
tags['wall1'] = Br 
tags['outflow'] = Bd 
tags['exterior'] = Br
tags['external'] = Br
#tags['land'] = Bo  
tags['land'] = Bd  
tags['westbank'] = None    #Riverbank
tags['eastbank'] = None    #Riverbank
tags['eastbankN'] = None    #Riverbank
domain.set_boundary(tags)


domain.check_integrity()

######################
#Evolution
t0 = time.time()
for t in domain.evolve(yieldstep = 0.2, finaltime = 300):
    domain.write_time()
print 'Simulation took %.2f seconds' %(time.time()-t0)
    
