import anuga
import numpy as num
import math
import inlet

class Structure_operator:
    """Structure Operator - transfer water from one rectangular box to another.
    Sets up the geometry of problem
    
    This is the base class for culverts. Inherit from this class (and overwrite
    compute_discharge method for specific subclasses)
    
    Input: Two points, pipe_size (either diameter or width, height), 
    mannings_rougness,
    """ 

    def __init__(self,
                 domain,
                 end_point0, 
                 end_point1,
                 width,
                 height,
                 apron,
                 manning,
                 enquiry_gap,
                 description,
                 verbose):
        
        self.domain = domain
        self.domain.set_fractional_step_operator(self)
        self.end_points = [end_point0, end_point1]
        
        if height is None:
            height = width

        if apron is None:
            apron = width

        self.width  = width
        self.height = height
        self.apron  = apron
        self.manning = manning
        self.enquiry_gap = enquiry_gap
        self.description = description
        self.verbose = verbose

        self.discharge = 0.0
        self.velocity = 0.0
        self.delta_total_energy = 0.0
        self.driving_energy = 0.0
        
        self.__create_exchange_polygons()

        self.inlets = []
        polygon0 = self.inlet_polygons[0]
        enquiry_point0 = self.inlet_equiry_points[0]
        outward_vector0 = self.culvert_vector
        self.inlets.append(inlet.Inlet(self.domain, polygon0, enquiry_point0, outward_vector0))

        polygon1 = self.inlet_polygons[1]
        exchange_polygon1 = self.inlet_equiry_points[1]
        outward_vector1  = - self.culvert_vector
        self.inlets.append(inlet.Inlet(self.domain, polygon1, exchange_polygon1, outward_vector1))

    def __call__(self):

        pass

    def __create_exchange_polygons(self):

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
        h = self.apron*self.culvert_vector    # Vector of length=height in the
                             # direction of the culvert

        gap = (1 + self.enquiry_gap)*h

        self.inlet_polygons = []
        self.inlet_equiry_points = []

        # Build exchange polygon and enquiry point
        for i in [0, 1]:
            i0 = (2*i-1)
            p0 = self.end_points[i] + w
            p1 = self.end_points[i] - w
            p2 = p1 + i0*h
            p3 = p0 + i0*h
            ep = self.end_points[i] + i0*gap

            self.inlet_polygons.append(num.array([p0, p1, p2, p3]))
            self.inlet_equiry_points.append(ep)

        # Check that enquiry points are outside inlet polygons
        for i in [0,1]:
            polygon = self.inlet_polygons[i]
            ep = self.inlet_equiry_points[i]
           
            area = anuga.polygon_area(polygon)
            
            msg = 'Polygon %s ' %(polygon)
            msg += ' has area = %f' % area
            assert area > 0.0, msg

            msg = 'Enquiry point falls inside an exchange polygon.'
            assert not anuga.inside_polygon(ep, polygon), msg
    
            
        #print '   outflow volume ',outflow.get_total_water_volume()
        

    def print_stats(self):

        print '====================================='
        print 'Generic Culvert Operator'
        print '====================================='

        print 'Culvert'
        print self.culvert

        print 'Culvert Routine'
        print self.routine
        
        for i, inlet in enumerate(self.inlets):
            print '-------------------------------------'
            print 'Inlet %i' % i
            print '-------------------------------------'

            print 'inlet triangle indices and centres'
            print inlet.triangle_indices
            print self.domain.get_centroid_coordinates()[inlet.triangle_indices]
        
            print 'polygon'
            print inlet.polygon

        print '====================================='


    def structure_statistics(self):

        message = '---------------------------\n'
        message += 'Structure report for structure %s:\n' % self.description
        message += '--------------------------\n'
        message += 'Discharge [m^3/s]: %.2f\n' % self.discharge
        message += 'Velocity  [m/s]: %.2f\n' % self.velocity
        message += 'Inlet Driving Energy %.2f\n' % self.driving_energy
        message += 'delta total energy %.2f\n' % self.delta_total_energy
#        message += 'Net boundary flow by tags [m^3/s]\n'
#        for tag in boundary_flows:
#            message += '    %s [m^3/s]: %.2f\n' % (tag, boundary_flows[tag])
#
#        message += 'Total net boundary flow [m^3/s]: %.2f\n' % \
#                    (total_boundary_inflow + total_boundary_outflow)
#        message += 'Total volume in domain [m^3]: %.2f\n' % \
#                    self.compute_total_volume()
#
#        # The go through explicit forcing update and record the rate of change
#        # for stage and
#        # record into forcing_inflow and forcing_outflow. Finally compute
#        # integral of depth to obtain total volume of domain.
#
        # FIXME(Ole): This part is not yet done.

        return message

    def get_inlets(self):
        
        return self.inlets
        
        
    def get_culvert_length(self):
        
        return self.culvert_length
        
        
    def get_culvert_width(self):
        
        return self.width
        
        
    def get_culvert_diameter(self):
    
            return self.width
        
        
    def get_culvert_height(self):
    
        return self.height


    def get_culvert_apron(self):

        return self.apron
