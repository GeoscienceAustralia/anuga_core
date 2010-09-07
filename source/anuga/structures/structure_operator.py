import anuga
import numpy as num
import math
import inlet

class Structure_operator:
    """Structure Operator - transfer water from one rectangular box to another.
    Sets up the geometry of problem
    
    This is the base class for structures (culverts, pipes, bridges etc). Inherit from this class (and overwrite
    discharge_routine method for specific subclasses)
    
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

        timestep = self.domain.get_timestep()
        
        self.__determine_inflow_outflow()
        
        Q, barrel_speed, outlet_depth = self.discharge_routine()

        old_inflow_height = self.inflow.get_average_height()
        old_inflow_xmom = self.inflow.get_average_xmom()
        old_inflow_ymom = self.inflow.get_average_ymom()
            
        if old_inflow_height > 0.0 :
                Q_star = Q/old_inflow_height
        else:
                Q_star = 0.0

        factor = 1.0/(1.0 + Q_star*timestep/self.inflow.get_area())

        new_inflow_height = old_inflow_height*factor
        new_inflow_xmom = old_inflow_xmom*factor
        new_inflow_ymom = old_inflow_ymom*factor
            
        self.inflow.set_heights(new_inflow_height)

        #inflow.set_xmoms(Q/inflow.get_area())
        #inflow.set_ymoms(0.0)

        self.inflow.set_xmoms(new_inflow_xmom)
        self.inflow.set_ymoms(new_inflow_ymom)

        loss = (old_inflow_height - new_inflow_height)*self.inflow.get_area()

        # set outflow
        if old_inflow_height > 0.0 :
                timestep_star = timestep*new_inflow_height/old_inflow_height
        else:
            timestep_star = 0.0

        outflow_extra_height = Q*timestep_star/self.outflow.get_area()
        outflow_direction = - self.outflow.outward_culvert_vector
        outflow_extra_momentum = outflow_extra_height*barrel_speed*outflow_direction
            
        gain = outflow_extra_height*self.outflow.get_area()
            
        #print Q, Q*timestep, barrel_speed, outlet_depth, Qstar, factor, timestep_star
        #print '  ', loss, gain

        # Stats
        self.discharge  = Q#outflow_extra_height*self.outflow.get_area()/timestep
        self.velocity = barrel_speed#self.discharge/outlet_depth/self.width

        new_outflow_height = self.outflow.get_average_height() + outflow_extra_height

        if self.use_momentum_jet :
            # FIXME (SR) Review momentum to account for possible hydraulic jumps at outlet
            #new_outflow_xmom = outflow.get_average_xmom() + outflow_extra_momentum[0]
            #new_outflow_ymom = outflow.get_average_ymom() + outflow_extra_momentum[1]

            new_outflow_xmom = barrel_speed*new_outflow_height*outflow_direction[0]
            new_outflow_ymom = barrel_speed*new_outflow_height*outflow_direction[1]

        else:
            #new_outflow_xmom = outflow.get_average_xmom()
            #new_outflow_ymom = outflow.get_average_ymom()

            new_outflow_xmom = 0.0
            new_outflow_ymom = 0.0

        self.outflow.set_heights(new_outflow_height)
        self.outflow.set_xmoms(new_outflow_xmom)
        self.outflow.set_ymoms(new_outflow_ymom)


    def __determine_inflow_outflow(self):
        # Determine flow direction based on total energy difference

        if self.use_velocity_head:
            self.delta_total_energy = self.inlets[0].get_enquiry_total_energy() - self.inlets[1].get_enquiry_total_energy()
        else:
            self.delta_total_energy = self.inlets[0].get_enquiry_stage() - self.inlets[1].get_enquiry_stage()


        self.inflow  = self.inlets[0]
        self.outflow = self.inlets[1]
        

        if self.delta_total_energy < 0:
            self.inflow  = self.inlets[1]
            self.outflow = self.inlets[0]
            self.delta_total_energy = -self.delta_total_energy


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
            
    
    def discharge_routine(self):
        
        pass
            

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
