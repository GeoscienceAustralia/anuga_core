from anuga.geometry.polygon import inside_polygon, polygon_area
from anuga.config import g
import anuga.utilities.log as log
import inlet
import numpy as num
import math

class Box_culvert:
    """Culvert flow - transfer water from one rectangular box to another.
    Sets up the geometry of problem
    
    This is the base class for culverts. Inherit from this class (and overwrite
    compute_discharge method for specific subclasses)
    
    Input: Two points, pipe_size (either diameter or width, height), 
    mannings_rougness,
    """	

    def __init__(self,
                 domain,
                 end_points, 
                 width=None,
                 height=None,
                 verbose=False):
        
        # Input check
        
        self.domain = domain

        self.end_points = end_points
 
        self.width = width
        self.height = height
        
        self.verbose=verbose
       
        # Create the fundamental culvert polygons and create inlet objects
        self.__create_culvert_polygons()

        #FIXME (SR) Put this into a foe loop to deal with more inlets
        self.inlets = []
        polygon0 = self.inlet_polygons[0]
        inlet0_vector = self.culvert_vector
        self.inlets.append(inlet.Inlet(self.domain, polygon0))

        polygon1 = self.inlet_polygons[1]
        inlet1_vector = - self.culvert_vector
        self.inlets.append(inlet.Inlet(self.domain, polygon1))
   

    def __create_culvert_polygons(self):

        """Create polygons at the end of a culvert inlet and outlet.
        At either end two polygons will be created; one for the actual flow to pass through and one a little further away
        for enquiring the total energy at both ends of the culvert and transferring flow.
        """

        # Calculate geometry
        x0, y0 = self.end_points[0]
        x1, y1 = self.end_points[1]

        dx = x1 - x0
        dy = y1 - y0

        self.culvert_vector = num.array([dx, dy])
        self.culvert_length = math.sqrt(num.sum(self.culvert_vector**2))
        assert self.culvert_length > 0.0, 'The length of culvert is less than 0'

        # Unit direction vector and normal
        self.culvert_vector /= self.culvert_length                      # Unit vector in culvert direction
        self.culvert_normal = num.array([-dy, dx])/self.culvert_length  # Normal vector

        # Short hands
        w = 0.5*self.width*self.culvert_normal # Perpendicular vector of 1/2 width
        h = self.height*self.culvert_vector    # Vector of length=height in the
                             # direction of the culvert

        self.inlet_polygons = []

        # Build exchange polygon and enquiry points 0 and 1
        for i in [0, 1]:
            i0 = (2*i-1)
            p0 = self.end_points[i] + w
            p1 = self.end_points[i] - w
            p2 = p1 + i0*h
            p3 = p0 + i0*h
            self.inlet_polygons.append(num.array([p0, p1, p2, p3]))

        # Check that enquiry points are outside inlet polygons
        for i in [0,1]:
            polygon = self.inlet_polygons[i]
            # FIXME (SR) Probably should calculate the area of all the triangles
            # associated with this polygon, as there is likely to be some
            # inconsistency between triangles and ploygon
            area = polygon_area(polygon)
            
            msg = 'Polygon %s ' %(polygon)
            msg += ' has area = %f' % area
            assert area > 0.0, msg
            
    
    def get_inlets(self):
        
        return self.inlets
        
        
    def get_culvert_length(self):
        
        return self.culvert_length
        
        
    def get_culvert_width(self):
        
        return self.width
        
        
    def get_culvert_height(self):
    
        return self.height
        
