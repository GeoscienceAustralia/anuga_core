"""Create mesh for University of Queensland dam break flume
"""


from anuga.pmesh.mesh import *
from anuga.coordinate_transforms.geo_reference import Geo_reference


def generate(mesh_filename, maximum_triangle_area=0.01):
    """
    Generate mesh for University of Aberbeen dam break flume.
    The gate position is the distance from the right wall of the wave tank
    to the gate. 
    
    """
    #Basic geometry
    
    xright  = 19.0
    ybottom = 0
    ytop    = 0.45
    xleft = 0.0
    xslope = 4.0 

    #Outline 
    point_sw = [xleft, ybottom]
    point_se = [xright, ybottom]
    point_nw = [xleft, ytop]    
    point_ne = [xright, ytop]


    # slope seperation (middle)
    point_slope_top = [xslope, ytop]
    point_slope_bottom = [xslope, ybottom]    

    m = Mesh()

    #Boundary
    points = [point_sw,   #se
              point_nw,
              point_slope_top,
              point_ne, 
              point_se,
              point_slope_bottom
              ]
    
    segments = [
        [0,1],
        [1,2],
        [2,3],
        [3,4],
        [4,5],
        [5,0],  #The outer border
        [2,5]       # slope separator
        ]    
    
    segment_tags = {'wall':[1,2,3,4,5],
                    'wave':[0]} # '':[6]
        
    m.add_points_and_segments(points, segments, segment_tags)
    
    dam = m.add_region(xslope - 0.0000001,(ytop - ybottom)/2)
    # this is the location of the reservoir region.
    dam.setTag("flat")
    
    slope = m.add_region(xslope + 0.0000001,(ytop - ybottom)/2)
    # this is the location of the slope region.
    slope.setTag("slope")
    
    m.generate_mesh(maximum_triangle_area=maximum_triangle_area)

    m.export_mesh_file(mesh_filename)
    print "mesh created"

#-------------------------------------------------------------
if __name__ == "__main__":
    generate("aa.tsh")
