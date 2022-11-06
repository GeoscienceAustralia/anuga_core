"""Functions for geometries related to culvert flows
"""

# Import necessary modules

from math import sqrt
from anuga.geometry.polygon import inside_polygon, polygon_area

import numpy as num


def create_culvert_polygons(end_point0,
                            end_point1, 
                            width, height=None,
                            enquiry_gap_factor=0.2,
                            number_of_barrels=1):
    """Create polygons at the end of a culvert inlet and outlet.
    At either end two polygons will be created; one for the actual flow to pass through and one a little further away
    for enquiring the total energy at both ends of the culvert and transferring flow.

    Input (mandatory):
        end_point0 - one end of the culvert (x,y)
        end_point1 - other end of the culvert (x,y)        
        width - culvert width

    Input (optional):        
        height - culvert height, defaults to width making a square culvert
        enquiry_gap_factor - sets the distance to the enquiry point as fraction of the height
        number_of_barrels - number of identical pipes.
        
    Output:

        Dictionary of four polygons. The dictionary keys are:
            'exchange_polygon0' - polygon defining the flow area at end_point0
            'exchange_polygon1' - polygon defining the flow area at end_point1
            'enquiry_point0' - point beyond exchange_polygon0
            'enquiry_point1' - point beyond exchange_polygon1            
            'vector'
            'length'
            'normal'
    """    


    # Input check
    if height is None:
        height = width

    # Dictionary for calculated polygons
    culvert_polygons = {}
    

    # Calculate geometry
    x0, y0 = end_point0
    x1, y1 = end_point1

    dx = x1-x0
    dy = y1-y0

    vector = num.array([dx, dy])
    length = sqrt(num.sum(vector**2))

    # Adjust polygon width to number of barrels in this culvert
    width *= number_of_barrels
    
    
    # Unit direction vector and normal 
    vector /= length                 # Unit vector in culvert direction
    normal = num.array([-dy, dx])/length # Normal vector
    
    culvert_polygons['vector'] = vector
    culvert_polygons['length'] = length
    culvert_polygons['normal'] = normal    

    # Short hands
    w = 0.5*width*normal # Perpendicular vector of 1/2 width 
    h = height*vector    # Vector of length=height in the
                         # direction of the culvert
    gap = (1 + enquiry_gap_factor)*h 
                         

    # Build exchange polygon and enquiry point for opening 0
    p0 = end_point0 + w
    p1 = end_point0 - w
    p2 = p1 - h
    p3 = p0 - h
    culvert_polygons['exchange_polygon0'] = num.array([p0,p1,p2,p3])
    culvert_polygons['enquiry_point0'] = end_point0 - gap
    

    # Build exchange polygon and enquiry point for opening 1
    p0 = end_point1 + w
    p1 = end_point1 - w
    p2 = p1 + h
    p3 = p0 + h
    culvert_polygons['exchange_polygon1'] = num.array([p0,p1,p2,p3])
    culvert_polygons['enquiry_point1'] = end_point1 + gap  

    # Check that enquiry polygons are outside exchange polygons
    for key1 in ['exchange_polygon0',
                 'exchange_polygon1']:
        polygon = culvert_polygons[key1]
        area = polygon_area(polygon)
        
        msg = 'Polygon %s ' %(polygon)
        msg += ' has area = %f' % area
        assert area > 0.0, msg

        for key2 in ['enquiry_point0', 'enquiry_point1']:
            point = culvert_polygons[key2]
            msg = 'Enquiry point falls inside exchange polygon.'
            msg += 'Email Ole.Nielsen@ga.gov.au'
            assert len(inside_polygon(point, polygon)) == 0, msg

    # Return results
    return culvert_polygons
