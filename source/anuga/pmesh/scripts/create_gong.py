
import os, sys
sys.path.append('..')

#problem, it's loading pyvolution mesh 1st!
from mesh import *
from coordinate_transforms.geo_reference import Geo_reference

#-------------------------------------------------------------
if __name__ == "__main__":

    # from 100m dem
    geo = Geo_reference(xllcorner = 241191.35087962,
                        yllcorner = 6111741.8729311,
                        zone = 56)
    m = Mesh(geo_reference=geo)
    bellambi_northing = 6195000
    north_beach_northing = 6183000
    border = 12000
    west = 301000 # From the map 
    inner_east = 310000 # From the map 
    inner_west_south = 303000 # From the map 
    inner_west_north = 307000  # From the map 
    east = 328421
    center = (east+inner_east)/2
    
    dict = {}
    dict['pointlist'] = [[west, north_beach_northing - border],  #sw
                         [west, bellambi_northing + border],  #nw
                         [center,
                          bellambi_northing + border],  #n_center
                         [east, bellambi_northing + border], #ne
                         [east, north_beach_northing - border],    #se
                         [center,
                          north_beach_northing - border]    #s_center
                         ]
    dict['segmentlist'] = [[0,1],[1,2],[2,3],
                           [3,4],[4,5],[5,0], # the outer boarder
                           [2,5]   # the center line
                           ]
    
    m.addVertsSegs(dict)

    dict = {}
    dict['pointlist'] = [
                         [inner_west_south,north_beach_northing], #sw
                         [inner_west_north,bellambi_northing],   #nw
                         [inner_east,bellambi_northing], #ne
                         [inner_east,north_beach_northing] #se
                         ]
    dict['segmentlist'] = [[0,1],[1,2],[2,3],[3,0]] # the inner boarder
    m.addVertsSegs(dict)

    factor = 10000 #low res 10000, high res 1000
    low = m.addRegionEN(east-0.5, bellambi_northing + border-0.5)
    low.setMaxArea(100*factor)
    
    medium = m.addRegionEN(center-0.5, bellambi_northing + border-0.5)
    medium.setMaxArea(10*factor)
    
    high = m.addRegionEN(inner_east-0.5,bellambi_northing-0.5)
    high.setMaxArea(1*factor)
    m.generateMesh()
    
    m.export_triangulation_file("wollongong_outline.msh")
    
