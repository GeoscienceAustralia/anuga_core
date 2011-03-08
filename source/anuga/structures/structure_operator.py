import anuga
import numpy as num
import math
import inlet_enquiry

from anuga.utilities.system_tools import log_to_file
from anuga.utilities.numerical_tools import ensure_numeric



class Structure_operator(anuga.Operator):
    """Structure Operator - transfer water from one rectangular box to another.
    Sets up the geometry of problem
    
    This is the base class for structures (culverts, pipes, bridges etc). Inherit from this class (and overwrite
    discharge_routine method for specific subclasses)
    
    Input: Two points, pipe_size (either diameter or width, depth),
    mannings_rougness,
    """ 

    counter = 0

    def __init__(self,
                 domain,
                 end_points,
                 exchange_lines,
                 enquiry_points,
                 width,
                 height,
                 apron,
                 manning,
                 enquiry_gap,
                 description,
                 label,
                 structure_type,
                 logging,
                 verbose):

        anuga.Operator.__init__(self,domain)
        
        self.end_points = ensure_numeric(end_points)
        self.exchange_lines = ensure_numeric(exchange_lines)
        self.enquiry_points = ensure_numeric(enquiry_points)

        
        if height is None:
            height = width

        if apron is None:
            apron = width

        self.width  = width
        self.height = height
        self.apron  = apron
        self.manning = manning
        self.enquiry_gap = enquiry_gap

        if description == None:
            self.description = ' '
        else:
            self.description = description
        
        if label == None:
            self.label = "structure_%g" % Structure_operator.counter
        else:
            self.label = label + '_%g' % Structure_operator.counter

        if structure_type == None:
            self.structure_type = 'generic structure'
        else:
            self.structure_type = structure_type
            
        self.verbose = verbose        
        
        # Keep count of structures
        Structure_operator.counter += 1

        # Slots for recording current statistics
        self.discharge = 0.0
        self.velocity = 0.0
        self.delta_total_energy = 0.0
        self.driving_energy = 0.0
        
        if exchange_lines is not None:
            self.__process_skew_culvert()
        elif end_points is not None:
            self.__process_non_skew_culvert()
        else:
            raise Exception, 'Define either exchange_lines or end_points'
        
        self.inlets = []
        line0 = self.exchange_lines[0] #self.inlet_lines[0]
        enquiry_point0 = self.enquiry_points[0]
        outward_vector0 = self.culvert_vector
        self.inlets.append(inlet_enquiry.Inlet_enquiry(self.domain, line0,
                           enquiry_point0, outward_vector0, self.verbose))

        line1 = self.exchange_lines[1]
        enquiry_point1 = self.enquiry_points[1]
        outward_vector1  = - self.culvert_vector
        self.inlets.append(inlet_enquiry.Inlet_enquiry(self.domain, line1,
                           enquiry_point1, outward_vector1, self.verbose))

        self.set_logging(logging)

    def __call__(self):

        timestep = self.domain.get_timestep()
        
        Q, barrel_speed, outlet_depth = self.discharge_routine()

        old_inflow_depth = self.inflow.get_average_depth()
        old_inflow_stage = self.inflow.get_average_stage()
        old_inflow_xmom = self.inflow.get_average_xmom()
        old_inflow_ymom = self.inflow.get_average_ymom()


        # Implement the update of flow over a timestep by
        # using a semi-implict update. This ensures that
        # the update does not create a negative depth
        if old_inflow_depth > 0.0 :
                Q_star = Q/old_inflow_depth
        else:
                Q_star = 0.0

        factor = 1.0/(1.0 + Q_star*timestep/self.inflow.get_area())

        new_inflow_depth = old_inflow_depth*factor
        new_inflow_xmom = old_inflow_xmom*factor
        new_inflow_ymom = old_inflow_ymom*factor
            
        self.inflow.set_depths(new_inflow_depth)

        #inflow.set_xmoms(Q/inflow.get_area())
        #inflow.set_ymoms(0.0)

        self.inflow.set_xmoms(new_inflow_xmom)
        self.inflow.set_ymoms(new_inflow_ymom)

        loss = (old_inflow_depth - new_inflow_depth)*self.inflow.get_area()

        # set outflow
        if old_inflow_depth > 0.0 :
                timestep_star = timestep*new_inflow_depth/old_inflow_depth
        else:
            timestep_star = 0.0




        outflow_extra_depth = Q*timestep_star/self.outflow.get_area()
        outflow_direction = - self.outflow.outward_culvert_vector
        outflow_extra_momentum = outflow_extra_depth*barrel_speed*outflow_direction
            
        gain = outflow_extra_depth*self.outflow.get_area()
            
        #print Q, Q*timestep, barrel_speed, outlet_depth, Qstar, factor, timestep_star
        #print '  ', loss, gain

        # Stats
        self.discharge  = Q#outflow_extra_depth*self.outflow.get_area()/timestep
        self.velocity = barrel_speed#self.discharge/outlet_depth/self.width

        new_outflow_depth = self.outflow.get_average_depth() + outflow_extra_depth

        if self.use_momentum_jet :
            # FIXME (SR) Review momentum to account for possible hydraulic jumps at outlet
            #new_outflow_xmom = outflow.get_average_xmom() + outflow_extra_momentum[0]
            #new_outflow_ymom = outflow.get_average_ymom() + outflow_extra_momentum[1]

            new_outflow_xmom = barrel_speed*new_outflow_depth*outflow_direction[0]
            new_outflow_ymom = barrel_speed*new_outflow_depth*outflow_direction[1]

        else:
            #new_outflow_xmom = outflow.get_average_xmom()
            #new_outflow_ymom = outflow.get_average_ymom()

            new_outflow_xmom = 0.0
            new_outflow_ymom = 0.0

        self.outflow.set_depths(new_outflow_depth)
        self.outflow.set_xmoms(new_outflow_xmom)
        self.outflow.set_ymoms(new_outflow_ymom)





    def __process_non_skew_culvert(self):

        """Create lines at the end of a culvert inlet and outlet.
        At either end two lines will be created; one for the actual flow to pass through and one a little further away
        for enquiring the total energy at both ends of the culvert and transferring flow.
        """
        
        self.culvert_vector = self.end_points[1] - self.end_points[0]
        self.culvert_length = math.sqrt(num.sum(self.culvert_vector**2))   
        assert self.culvert_length > 0.0, 'The length of culvert is less than 0'
        
        self.culvert_vector /= self.culvert_length
        
        culvert_normal = num.array([-self.culvert_vector[1], self.culvert_vector[0]])  # Normal vector
        w = 0.5*self.width*culvert_normal # Perpendicular vector of 1/2 width

        self.exchange_lines = []

        # Build exchange polyline and enquiry point
        if self.enquiry_points is None:
            
            gap = (self.apron + self.enquiry_gap)*self.culvert_vector
            self.enquiry_points = []
            
            for i in [0, 1]:
                p0 = self.end_points[i] + w
                p1 = self.end_points[i] - w
                self.exchange_lines.append(num.array([p0, p1]))
                ep = self.end_points[i] + (2*i - 1)*gap #(2*i - 1) determines the sign of the points
                self.enquiry_points.append(ep)
            
        else:            
            for i in [0, 1]:
                p0 = self.end_points[i] + w
                p1 = self.end_points[i] - w
                self.exchange_lines.append(num.array([p0, p1]))
            
  
    def __process_skew_culvert(self):    
        
        """Compute skew culvert.
        If exchange lines are given, the enquiry points are determined. This is for enquiring 
        the total energy at both ends of the culvert and transferring flow.
        """
            
        centre_point0 = 0.5*(self.exchange_lines[0][0] + self.exchange_lines[0][1])
        centre_point1 = 0.5*(self.exchange_lines[1][0] + self.exchange_lines[1][1])
        
        if self.end_points is None:
            self.culvert_vector = centre_point1 - centre_point0
        else:
            self.culvert_vector = self.end_points[1] - end_points[0]
        
        self.culvert_length = math.sqrt(num.sum(self.culvert_vector**2))
        assert self.culvert_length > 0.0, 'The length of culvert is less than 0'
        
        if self.enquiry_points is None:
        
            self.culvert_vector /= self.culvert_length
            gap = (self.apron + self.enquiry_gap)*self.culvert_vector
        
            self.enquiry_points = []

            self.enquiry_points.append(centre_point0 - gap)
            self.enquiry_points.append(centre_point1 + gap)
            

    def discharge_routine(self):

        msg = 'Need to impelement '
        raise
            

    def statistics(self):


        message  = '=====================================\n'
        message += 'Structure Operator: %s\n' % self.label
        message += '=====================================\n'

        message += 'Structure Type: %s\n' % self.structure_type

        message += 'Description\n'
        message += '%s' % self.description
        message += '\n'
        
        for i, inlet in enumerate(self.inlets):
            message += '-------------------------------------\n'
            message +=  'Inlet %i\n' % i
            message += '-------------------------------------\n'

            message += 'inlet triangle indices and centres\n'
            message += '%s' % inlet.triangle_indices
            message += '\n'
            
            message += '%s' % self.domain.get_centroid_coordinates()[inlet.triangle_indices]
            message += '\n'

            message += 'line\n'
            message += '%s' % inlet.line
            message += '\n'

        message += '=====================================\n'

        return message


    def print_statistics(self):

        print self.statistics()


    def print_timestepping_statistics(self):

        message = '---------------------------\n'
        message += 'Structure report for %s:\n' % self.label
        message += '--------------------------\n'
        message += 'Type: %s\n' % self.structure_type
        message += 'Discharge [m^3/s]: %.2f\n' % self.discharge
        message += 'Velocity  [m/s]: %.2f\n' % self.velocity
        message += 'Inlet Driving Energy %.2f\n' % self.driving_energy
        message += 'Delta Total Energy %.2f\n' % self.delta_total_energy
        message += 'Control at this instant: %s\n' % self.case

        print message


    def set_logging(self, flag=True):

        self.logging = flag

        # If flag is true open file with mode = "w" to form a clean file for logging
        if self.logging:
            self.log_filename = self.label + '.log'
            log_to_file(self.log_filename, self.statistics(), mode='w')
            log_to_file(self.log_filename, 'time,discharge,velocity,driving_energy,delta_total_energy')

            #log_to_file(self.log_filename, self.culvert_type)


    def timestepping_statistics(self):

        message  = '%.5f, ' % self.domain.get_time()
        message += '%.5f, ' % self.discharge
        message += '%.5f, ' % self.velocity
        message += '%.5f, ' % self.driving_energy
        message += '%.5f' % self.delta_total_energy

        return message

    def log_timestepping_statistics(self):

         if self.logging:
             log_to_file(self.log_filename, self.timestepping_statistics())


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
