"""Functions for geometries related to culvert flows
"""

# Import necessary modules
from math import atan
from math import sin
from math import cos
from math import pi
from math import degrees


def create_culvert_polygons(center1, center2,
                            width, height=None,
                            enquiry_gap_factor=1.0,
                            enquiry_shape_factor=2.0):
    """Create polygons at the end of a culvert inlet and outlet.
    At either end two polygons will be created; one for the actual flow to pass through and one a little further away
    for enquiring the total energy at both ends of the culvert and transferring flow.

    Input (mandatory):
        center1 - center point at one end of the culvert (x,y)
        center2 - center point at other end of the culvert (x,y)        
        width - culvert width

    Input (optional):        
        height - culvert height, defaults to width making a square culvert
        enquiry_gap_factor - sets the distance to the enquiry polygon
        enquiry_shape_factor - sets the shape of the enquiry polygon
        
    Output:

        Dictionary of four polygons. The dictionary keys are:
            'exchange_polygon1' - polygon defining the flow area at center1
            'exchange_polygon2' - polygon defining the flow area at center2
            'enquiry_polygon1' - polygon defining the enquiry field near center1
            'enquiry_polygon2' - polygon defining the enquiry field near center2            
    
    """    


    # Input check
    if height is None:
        height = width


    # Calculate geometry
    x1, y1 = center1
    x2, y2 = center2

    dx=x2-x1
    dy=y2-y1


    # Calculate slope a in y=a*x+b to project points along the direction of the culvert

    if dx==0.0:
        # vertical
        slope=0.0
    else:
        slope=dy/dx


    #  Intercept    b=Y1-slope*X1


    # Calculate basic trigonometric values
    
    # FIXME (Ole): Suggest using
    # from anuga.utilities.numerical_tools import angle
    # theta = angle(slope)
    # sin_theta=sin(theta)
    # cos_theta=cos(theta)    


    angle_atan=atan(slope)
    angle_sin=sin(angle_atan)
    angle_cos=cos(angle_atan)

    # If you need to establish orientation USE THIS other Wise SKIP !!!
    if dx == 0.0:
        if dy > 0.0:        # vertical up
            angle_rad = 0.0
        else:               # dy < 0 vertical down
    	    angle_rad = pi  # 180 Degrees
    elif dx<0:
        if slope > 0:       # Pt 2 less than Pt 1
            angle_rad=3*pi/2-angle_atan
    else:
        if slope > 0:
            angle_rad=pi/2-angle_atan

    # Calculate derived perpendicular trigonometric values
    angle_per=atan(angle_rad)
    angle_per_sin=sin(angle_rad)
    angle_per_cos=cos(angle_rad)
    
    angle_deg = degrees(angle_rad) # Convert to degrees



    # From the 2 X,Y co-ordinates at each end of the Culvert 
    # establish Polygons
    #   - to abstract flow and/or inject flow, 
    #	- and set up enquiry area for total energy

    culvert_polygons = {}


    # Exchange polygon 1
    pt1x1=x1
    pt1y1=y1
    pt2x1=pt1x1-angle_sin*width/2
    pt2y1=pt1y1+angle_cos*width/2
    pt3x1=pt2x1-angle_per_sin*height
    pt3y1=pt2y1-angle_per_cos*height
    pt4x1=pt3x1+angle_sin*width
    pt4y1=pt3y1-angle_cos*width
    pt5x1=pt4x1+angle_per_sin*height
    pt5y1=pt4y1+angle_per_cos*height
    pt6x1=pt5x1-angle_sin*width/2
    pt6y1=pt5y1+angle_cos*width/2

    culvert_polygons['exchange_polygon1'] = [[pt1x1,pt1y1], [pt2x1,pt2y1],
                                             [pt3x1,pt3y1], [pt4x1,pt4y1],
                                             [pt5x1,pt5y1]]


    # Exchange polygon 2
    pt1x2=x2
    pt1y2=y2
    pt2x2=pt1x2+angle_sin*width/2
    pt2y2=pt1y2-angle_cos*width/2
    pt3x2=pt2x2+angle_per_sin*height
    pt3y2=pt2y2+angle_per_cos*height
    pt4x2=pt3x2-angle_sin*width
    pt4y2=pt3y2+angle_cos*width
    pt5x2=pt4x2-angle_per_sin*height
    pt5y2=pt4y2-angle_per_cos*height
    pt6x2=pt5x2+angle_sin*width/2
    pt6y2=pt5y2-angle_cos*width/2

    culvert_polygons['exchange_polygon2'] = [[pt1x2,pt1y2], [pt2x2,pt2y2],
                                             [pt3x2,pt3y2], [pt4x2,pt4y2],
                                             [pt5x2,pt5y2]]


    # Calculate the the energy enquiry fields
    # using a couple of distances first using factors from above

    enq_dist=height*(1+enquiry_gap_factor)
    enq_width=width*enquiry_shape_factor
    enq_height=height/enquiry_shape_factor



    # Enquiry polygon 1
    pt1qx1=x1-angle_per_sin*enq_dist
    pt1qy1=y1-angle_per_cos*enq_dist
    pt2qx1=pt1qx1-angle_sin*enq_width/2
    pt2qy1=pt1qy1+angle_cos*enq_width/2
    pt3qx1=pt2qx1-angle_per_sin*enq_height
    pt3qy1=pt2qy1-angle_per_cos*enq_height
    pt4qx1=pt3qx1+angle_sin*enq_width
    pt4qy1=pt3qy1-angle_cos*enq_width
    pt5qx1=pt4qx1+angle_per_sin*enq_height
    pt5qy1=pt4qy1+angle_per_cos*enq_height
    pt6qx1=pt5qx1-angle_sin*enq_width/2
    pt6qy1=pt5qy1+angle_cos*enq_width/2


    culvert_polygons['enquiry_polygon1'] = [[pt1qx1,pt1qy1], [pt2qx1,pt2qy1],
                                            [pt3qx1,pt3qy1], [pt4qx1,pt4qy1],
                                            [pt5qx1,pt5qy1]]


    # Enquiry polygon 2
    pt1qx2=x2+angle_per_sin*enq_dist
    pt1qy2=y2+angle_per_cos*enq_dist
    pt2qx2=pt1qx2+angle_sin*enq_width/2
    pt2qy2=pt1qy2-angle_cos*enq_width/2
    pt3qx2=pt2qx2+angle_per_sin*enq_height
    pt3qy2=pt2qy2+angle_per_cos*enq_height
    pt4qx2=pt3qx2-angle_sin*enq_width
    pt4qy2=pt3qy2+angle_cos*enq_width
    pt5qx2=pt4qx2-angle_per_sin*enq_height
    pt5qy2=pt4qy2-angle_per_cos*enq_height
    pt6qx2=pt5qx2+angle_sin*enq_width/2
    pt6qy2=pt5qy2-angle_cos*enq_width/2


    culvert_polygons['enquiry_polygon2'] = [[pt1qx2,pt1qy2], [pt2qx2,pt2qy2],
                                            [pt3qx2,pt3qy2], [pt4qx2,pt4qy2],
                                            [pt5qx2,pt5qy2]]



    # Return dictionary of all four polygons
    return culvert_polygons
