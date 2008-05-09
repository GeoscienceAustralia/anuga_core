"""Functions for geometries related to culvert flows
"""

# Import necessary modules
from math import sqrt
from Numeric import array, sum

def create_culvert_polygons(end_point0,
                            end_point1, 
                            width, height=None,
                            enquiry_gap_factor=1.0,
                            enquiry_shape_factor=2.0):
    """Create polygons at the end of a culvert inlet and outlet.
    At either end two polygons will be created; one for the actual flow to pass through and one a little further away
    for enquiring the total energy at both ends of the culvert and transferring flow.

    Input (mandatory):
        end_point0 - one end of the culvert (x,y)
        end_point1 - other end of the culvert (x,y)        
        width - culvert width

    Input (optional):        
        height - culvert height, defaults to width making a square culvert
        enquiry_gap_factor - sets the distance to the enquiry polygon
        enquiry_shape_factor - sets the shape of the enquiry polygon
                               (large value widens polygon but reduces height
                               to preserve same area as exchange polygon)
        
    Output:

        Dictionary of four polygons. The dictionary keys are:
            'exchange_polygon0' - polygon defining the flow area at end_point0
            'exchange_polygon1' - polygon defining the flow area at end_point1
            'enquiry_polygon0' - polygon defining the enquiry field near end_point0
            'enquiry_polygon1' - polygon defining the enquiry field near end_point1           
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

    vector = array([dx, dy])
    length = sqrt(sum(vector**2))

    # Unit direction vector and normal 
    vector /= length
    normal = array([-dy, dx])/length
    

    # Short hands
    w = width/2*normal # Perpendicular vector of 1/2 width 
    h = height*vector  # Vector of length=height in the
                       # direction of the culvert

    # Build exchange polygon 0
    p0 = end_point0 + w
    p1 = end_point0 - w
    p2 = p1 - h
    p3 = p0 - h
    culvert_polygons['exchange_polygon0'] = array([p0,p1,p2,p3])

    
    # Build exchange polygon 1
    p0 = end_point1 + w
    p1 = end_point1 - w
    p2 = p1 + h
    p3 = p0 + h
    culvert_polygons['exchange_polygon1'] = array([p0,p1,p2,p3])



    # Redefine shorthands for enquiry polygons
    w = w*enquiry_shape_factor
    h = h/enquiry_shape_factor
    gap = (enquiry_gap_factor + h)*vector
    
    # Build enquiry polygon 0
    p0 = end_point0 + w - gap 
    p1 = end_point0 - w - gap
    p2 = p1 - h
    p3 = p0 - h
    culvert_polygons['enquiry_polygon0'] = array([p0,p1,p2,p3])

    # Build enquiry polygon 1
    p0 = end_point1 + w + gap 
    p1 = end_point1 - w + gap
    p2 = p1 + h
    p3 = p0 + h
    culvert_polygons['enquiry_polygon1'] = array([p0,p1,p2,p3])    

    # Return results
    culvert_polygons['vector'] = vector
    culvert_polygons['length'] = length
    culvert_polygons['normal'] = normal    
    return culvert_polygons
