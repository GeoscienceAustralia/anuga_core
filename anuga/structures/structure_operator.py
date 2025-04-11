
import anuga
import numpy as num
import math
from . import inlet_enquiry

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
                 end_points=None,
                 exchange_lines=None,
                 enquiry_points=None,
                 invert_elevations=None,
                 width=None,
                 height=None,
                 diameter=None,
                 z1=None,
                 z2=None,
                 blockage=None,
                 barrels=None,
                 #culvert_slope=None,
                 apron=None,
                 manning=None,
                 enquiry_gap=None,
                 use_momentum_jet=False,
                 zero_outflow_momentum=True,
                 use_old_momentum_method=True,
                 always_use_Q_wetdry_adjustment=True,
                 force_constant_inlet_elevations=False,
                 description=None,
                 label=None,
                 structure_type=None,
                 logging=None,
                 verbose=None):
                     
        """
        exchange_lines define the input lines for each inlet.

        If end_points = None, then the culvert_vector is calculated in the
        directions from the centre of echange_line[0] to centre of exchange_line[1}

        If end_points != None, then culvert_vector is unit vector in direction
        end_point[1] - end_point[0]
        """

        anuga.Operator.__init__(self,domain)

        self.master_proc = 0
        self.end_points = ensure_numeric(end_points)
        self.exchange_lines = ensure_numeric(exchange_lines)
        self.enquiry_points = ensure_numeric(enquiry_points)
        self.invert_elevations = ensure_numeric(invert_elevations)

        assert self.end_points is None or self.exchange_lines is None

        
        if height is None:
            height = width

        if width is None:
            width = diameter

        if apron is None:
            apron = width


        assert width is not None


        self.width  = width
        self.height = height
        self.diameter = diameter
        self.z1 = z1 
        self.z2 = z2 
        self.blockage = blockage 
        self.barrels = barrels
        self.apron  = apron
        self.manning = manning
        self.enquiry_gap = enquiry_gap
        self.use_momentum_jet = use_momentum_jet
        self.zero_outflow_momentum = zero_outflow_momentum
        if use_momentum_jet and zero_outflow_momentum:
            msg = "Can't have use_momentum_jet and zero_outflow_momentum both True"
            raise Exception(msg)
        self.use_old_momentum_method = use_old_momentum_method
        self.always_use_Q_wetdry_adjustment = always_use_Q_wetdry_adjustment


        if description is None:
            self.description = ' '
        else:
            self.description = description
        
        if label is None:
            self.label = "structure_%g" % Structure_operator.counter
        else:
            self.label = label + '_%g' % Structure_operator.counter

        if structure_type is None:
            self.structure_type = 'generic structure'
        else:
            self.structure_type = structure_type
            
        self.verbose = verbose        
        
        # Keep count of structures
        Structure_operator.counter += 1

        # Slots for recording current statistics
        self.accumulated_flow = 0.0
        self.discharge = 0.0
        self.discharge_abs_timemean = 0.0
        self.velocity = 0.0
        self.outlet_depth = 0.0
        self.delta_total_energy = 0.0
        self.driving_energy = 0.0
        
        if exchange_lines is not None:
            self.__process_skew_culvert()
        elif end_points is not None:
            self.__process_non_skew_culvert()
        else:
            raise Exception('Define either exchange_lines or end_points')
        

        self.inlets = []
        line0 = self.exchange_lines[0] #self.inlet_lines[0]
        if self.apron is None:
            poly0 = line0
        else:
            offset = -self.apron*self.outward_vector_0 
            #print line0
            #print offset
            poly0 = num.array([ line0[0], line0[1], line0[1]+offset, line0[0]+offset])
            #print poly0
        if self.invert_elevations is None:
            invert_elevation0 = None
        else:
            invert_elevation0 = self.invert_elevations[0]

        enquiry_point0 = self.enquiry_points[0]

        #outward_vector0 = - self.culvert_vector
        self.inlets.append(inlet_enquiry.Inlet_enquiry(
                           self.domain,
                           poly0,
                           enquiry_point0,
                           invert_elevation = invert_elevation0,
                           outward_culvert_vector = self.outward_vector_0,
                           verbose = self.verbose))

        if force_constant_inlet_elevations:
            # Try to enforce a constant inlet elevation 
            inlet_global_elevation = self.inlets[-1].get_average_elevation() 
            self.inlets[-1].set_elevations(inlet_global_elevation)

        tris_0 = self.inlets[0].triangle_indices
        #print tris_0
        #print self.domain.centroid_coordinates[tris_0]

        line1 = self.exchange_lines[1]
        if self.apron is None:
            poly1 = line1
        else:
            offset = -self.apron*self.outward_vector_1
            #print line1
            #print offset
            poly1 = num.array([ line1[0], line1[1], line1[1]+offset, line1[0]+offset])
            #print poly1
            
        if self.invert_elevations is None:
            invert_elevation1 = None
        else:
            invert_elevation1 = self.invert_elevations[1]
        enquiry_point1 = self.enquiry_points[1]

        self.inlets.append(inlet_enquiry.Inlet_enquiry(
                           self.domain,
                           poly1,
                           enquiry_point1,
                           invert_elevation = invert_elevation1,
                           outward_culvert_vector = self.outward_vector_1,
                           verbose = self.verbose))

        if force_constant_inlet_elevations:
            # Try to enforce a constant inlet elevation 
            inlet_global_elevation = self.inlets[-1].get_average_elevation() 
            self.inlets[-1].set_elevations(inlet_global_elevation)

        tris_1 = self.inlets[1].triangle_indices
        
        self.set_logging(logging)

        



    def __call__(self):

        timestep = self.domain.get_timestep()

        # implement different types of structures by redefining 
        # the discharge_routine
        Q, barrel_speed, outlet_depth = self.discharge_routine()

        old_inflow_depth = self.inflow.get_average_depth()
        old_inflow_stage = self.inflow.get_average_stage()
        old_inflow_xmom = self.inflow.get_average_xmom()
        old_inflow_ymom = self.inflow.get_average_ymom()

        # Implement the update of flow over a timestep by
        # using a semi-implict update. This ensures that
        # the update does not create a negative depth
        if old_inflow_depth > 0.0 :
            dt_Q_on_d = timestep*Q/old_inflow_depth
        else:
            dt_Q_on_d = 0.0

        always_use_Q_wetdry_adjustment = self.always_use_Q_wetdry_adjustment
        use_Q_wetdry_adjustment = ((always_use_Q_wetdry_adjustment) |\
            (old_inflow_depth*self.inflow.get_area() <= Q*timestep))
        # If we use the 'Q_adjustment', then the discharge is rescaled so that
        # the depth update is: 
        #    new_inflow_depth*inflow_area = 
        #    old_inflow_depth*inflow_area - 
        #    timestep*Q*(new_inflow_depth/old_inflow_depth)
        # The last term in () is a wet-dry improvement trick (which rescales Q)
        #
        # Before Feb 2015 this rescaling was always done (even if the flow was
        # not near wet-dry). Now if use_old_Q_adjustment=False, the rescaling
        # is only done if required to avoid drying
        #
        #
        factor = 1.0/(1.0 + dt_Q_on_d/self.inflow.get_area())


        if use_Q_wetdry_adjustment:
            # FIXME:
            new_inflow_depth = old_inflow_depth*factor

            if(old_inflow_depth>0.):
                timestep_star = timestep*new_inflow_depth/old_inflow_depth
            else:
                timestep_star = 0.

        else:
            new_inflow_depth = old_inflow_depth - timestep*Q/self.inflow.get_area()
            timestep_star = timestep

        if(self.use_old_momentum_method):
            # This method is here for consistency with the old version of the
            # routine
            new_inflow_xmom = old_inflow_xmom*factor
            new_inflow_ymom = old_inflow_ymom*factor

        else:
            # For the momentum balance, note that Q also advects the momentum,
            # The volumetric momentum flux should be Q*momentum/depth, where
            # momentum has an average value of new_inflow_mom (or
            # old_inflow_mom).  We use old_inflow_depth for depth
            #
            #     new_inflow_xmom*inflow_area = 
            #     old_inflow_xmom*inflow_area - 
            #     [timestep*Q]*new_inflow_xmom/old_inflow_depth
            # and:
            #     new_inflow_ymom*inflow_area = 
            #     old_inflow_ymom*inflow_area - 
            #     [timestep*Q]*new_inflow_ymom/old_inflow_depth
            #
            # The units balance: m^2/s*m^2 = m^2/s*m^2 - s*m^3/s*m^2/s *m^(-1)
            #
            if old_inflow_depth > 0.:
                if use_Q_wetdry_adjustment:
                    # Replace dt*Q with dt*Q*new_inflow_depth/old_inflow_depth = dt_Q_on_d*new_inflow_depth
                    factor2 = 1.0/(1.0 + dt_Q_on_d*new_inflow_depth/(old_inflow_depth*self.inflow.get_area()))
                else:
                    factor2 = 1.0/(1.0 + timestep*Q/(old_inflow_depth*self.inflow.get_area()))
            else:
                factor2 = 0.

            new_inflow_xmom = old_inflow_xmom*factor2
            new_inflow_ymom = old_inflow_ymom*factor2
        
        self.inflow.set_depths(new_inflow_depth)

        #inflow.set_xmoms(Q/inflow.get_area())
        #inflow.set_ymoms(0.0)

        self.inflow.set_xmoms(new_inflow_xmom)
        self.inflow.set_ymoms(new_inflow_ymom)

        loss = (old_inflow_depth - new_inflow_depth)*self.inflow.get_area()
        xmom_loss = (old_inflow_xmom - new_inflow_xmom)*self.inflow.get_area()
        ymom_loss = (old_inflow_ymom - new_inflow_ymom)*self.inflow.get_area()

        # set outflow
        outflow_extra_depth = Q*timestep_star/self.outflow.get_area()
        outflow_direction = - self.outflow.outward_culvert_vector
        #outflow_extra_momentum = outflow_extra_depth*barrel_speed*outflow_direction
            
        gain = outflow_extra_depth*self.outflow.get_area()
        
        #print gain, loss
        assert num.allclose(gain-loss, 0.0)
            
        # Stats

        #print('gain', gain)

        self.accumulated_flow += gain
        self.discharge  = Q*timestep_star/timestep 
        self.discharge_abs_timemean += gain/self.domain.yieldstep
        self.velocity =   barrel_speed
        self.outlet_depth = outlet_depth
        

        new_outflow_depth = self.outflow.get_average_depth() + outflow_extra_depth

        self.outflow.set_depths(new_outflow_depth)

        if self.use_momentum_jet:
            # FIXME (SR) Review momentum to account for possible hydraulic jumps at outlet
            # FIXME (GD) Depending on barrel speed I think this will be either
            # a source or sink of momentum (considering the momentum losses
            # above). Might not always be reasonable.
            #new_outflow_xmom = self.outflow.get_average_xmom() + outflow_extra_momentum[0]
            #new_outflow_ymom = self.outflow.get_average_ymom() + outflow_extra_momentum[1]
            new_outflow_xmom = barrel_speed*new_outflow_depth*outflow_direction[0]
            new_outflow_ymom = barrel_speed*new_outflow_depth*outflow_direction[1]
            
        elif self.zero_outflow_momentum:
            new_outflow_xmom = 0.0
            new_outflow_ymom = 0.0
            #new_outflow_xmom = outflow.get_average_xmom()
            #new_outflow_ymom = outflow.get_average_ymom()

        else:
            # Add the momentum lost from the inflow to the outflow. For
            # structures where barrel_speed is unknown + direction doesn't
            # change from inflow to outflow
            new_outflow_xmom = self.outflow.get_average_xmom() + xmom_loss/self.outflow.get_area()
            new_outflow_ymom = self.outflow.get_average_ymom() + ymom_loss/self.outflow.get_area()

        self.outflow.set_xmoms(new_outflow_xmom)
        self.outflow.set_ymoms(new_outflow_ymom)



    def set_culvert_height(self, height):

        self.culvert_height = height

    def set_culvert_width(self, width):

        self.culvert_width = width
        
    def set_culvert_z1(self, z1): 

        self.culvert_z1 = z1 

    def set_culvert_z2(self, z2):

        self.culvert_z2 = z2
        
    def set_culvert_blockage(self, blockage): 

        self.culvert_blockage = blockage 

    def set_culvert_barrels(self, barrels): 

        self.culvert_barrels = barrels 
        
        
    def __process_non_skew_culvert(self):

        """Create lines at the end of a culvert inlet and outlet.
        At either end two lines will be created; one for the actual flow to pass through and one a little further away
        for enquiring the total energy at both ends of the culvert and transferring flow.
        """
        
        self.culvert_vector = self.end_points[1] - self.end_points[0]
        self.culvert_length = math.sqrt(num.sum(self.culvert_vector**2))   
        assert self.culvert_length > 0.0, 'The length of culvert is less than 0'
        
        self.culvert_vector /= self.culvert_length
        self.outward_vector_0 =   self.culvert_vector
        self.outward_vector_1 = - self.culvert_vector

        
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

        n_exchange_0 = len(self.exchange_lines[0])
        n_exchange_1 = len(self.exchange_lines[1])

        assert n_exchange_0 == n_exchange_1, 'There should be the same number of points in both exchange_lines'

        if n_exchange_0 == 2:
        
            if self.end_points is None:
                self.culvert_vector = centre_point1 - centre_point0
            else:
                self.culvert_vector = self.end_points[1] - self.end_points[0]

            self.outward_vector_0 =   self.culvert_vector
            self.outward_vector_1 = - self.culvert_vector


        elif n_exchange_0 == 4:

            self.outward_vector_0 = self.exchange_lines[0][3] - self.exchange_lines[0][2]
            self.outward_vector_1 = self.exchange_lines[1][3] - self.exchange_lines[1][2]

            self.culvert_vector = centre_point1 - centre_point0

        else:
            raise Exception('n_exchange_0 != 2 or 4')


        self.culvert_length = math.sqrt(num.sum(self.culvert_vector**2))
        assert self.culvert_length > 0.0, 'The length of culvert is less than 0'
        self.culvert_vector /= self.culvert_length

        outward_vector_0_length = math.sqrt(num.sum(self.outward_vector_0**2))
        assert outward_vector_0_length > 0.0, 'The length of outlet_vector_0 is less than 0'
        self.outward_vector_0 /= outward_vector_0_length

        outward_vector_1_length = math.sqrt(num.sum(self.outward_vector_1**2))
        assert outward_vector_1_length > 0.0, 'The length of outlet_vector_1 is less than 0'
        self.outward_vector_1 /= outward_vector_1_length


        if self.enquiry_points is None:
        
            gap = (self.apron + self.enquiry_gap)*self.culvert_vector
        
            self.enquiry_points = []

            self.enquiry_points.append(centre_point0 - gap)
            self.enquiry_points.append(centre_point1 + gap)
            

    def discharge_routine(self):

        msg = 'Need to implement '
        raise
            

    def statistics(self):


        message  = '=====================================\n'
        message += 'Structure Operator: %s\n' % self.label
        message += '=====================================\n'

        message += 'Structure Type: %s\n' % self.structure_type

        message += 'Description\n'
        message += '%s' % self.description
        
        #add the culvert dimensions, blockage factor here
        if self.structure_type == 'boyd_pipe':
            message += 'Culvert Diameter: %s\n'% self.diameter
            message += 'Culvert Blockage: %s\n'% self.blockage
            message += 'No.  of  barrels: %s\n'% self.barrels
        elif self.structure_type == 'boyd_box':
            message += 'Culvert  Height: %s\n'% self.height
            message += 'Culvert    Width: %s\n'% self.width
            message += 'Culvert Blockage: %s\n'% self.blockage
            message += 'No.  of  barrels: %s\n'% self.barrels
        else:
            message += 'Culvert Height: %s\n'% self.height
            message += 'Culvert  Width: %s\n'% self.width
            message += 'Batter Slope 1: %s\n'% self.z1
            message += 'Batter Slope 2: %s\n'% self.z2
            message += 'Culvert Blockage: %s\n'% self.blockage
            message += 'No.  of  barrels: %s\n'% self.barrels
            
        message += '\n'
        
        for i, inlet in enumerate(self.inlets):
            message += '-------------------------------------\n'
            message +=  'Inlet %i\n' % i
            message += '-------------------------------------\n'

            message += 'inlet triangle indices and centres and elevations\n'
            message += '%s' % inlet.triangle_indices
            message += '\n'
            
            message += '%s' % self.domain.get_centroid_coordinates()[inlet.triangle_indices]
            message += '\n'

            elev = self.domain.quantities['elevation'].centroid_values[inlet.triangle_indices]
            message += '%s' % elev
            message += '\n'
           
            elevation_range = elev.max() - elev.min() 
            if not num.allclose(elevation_range, 0.):
                message += 'Warning: non-constant inlet elevation can cause well-balancing problems'

            message += 'region\n'
            message += '%s' % inlet.region
            message += '\n'

        message += '=====================================\n'

        return message


    def print_statistics(self):

        print(self.statistics())


    def print_timestepping_statistics(self):

        message = '---------------------------\n'
        message += 'Structure report for %s:\n' % self.label
        message += '--------------------------\n'
        message += 'Type: %s\n' % self.structure_type
        
        message += 'inlets[0]_enquiry_depth [m]:  %.2f\n' %self.inlets[0].get_enquiry_depth()
        message += 'inlets[0]_enquiry_speed [m/s]:  %.2f\n' %self.inlets[0].get_enquiry_speed()
        message += 'inlets[0]_enquiry_stage [m]:  %.2f\n' %self.inlets[0].get_enquiry_stage()
        message += 'inlets[0]_enquiry_elevation [m]:  %.2f\n' %self.inlets[0].get_enquiry_elevation()
        message += 'inlets[0]_average_depth [m]:  %.2f\n' %self.inlets[0].get_average_depth()
        message += 'inlets[0]_average_speed [m/s]:  %.2f\n' %self.inlets[0].get_average_speed()
        message += 'inlets[0]_average_stage [m]:  %.2f\n' %self.inlets[0].get_average_stage()
        message += 'inlets[0]_average_elevation [m]:  %.2f\n' %self.inlets[0].get_average_elevation()

        message += '\n'
       
        message += 'inlets[1]_enquiry_depth [m]:  %.2f\n' %self.inlets[1].get_enquiry_depth()
        message += 'inlets[1]_enquiry_speed [m/s]:  %.2f\n' %self.inlets[1].get_enquiry_speed()
        message += 'inlets[1]_enquiry_stage [m]:  %.2f\n' %self.inlets[1].get_enquiry_stage()
        message += 'inlets[1]_enquiry_elevation [m]:  %.2f\n' %self.inlets[1].get_enquiry_elevation()

        message += 'inlets[1]_average_depth [m]:  %.2f\n' %self.inlets[1].get_average_depth()
        message += 'inlets[1]_average_speed [m/s]:  %.2f\n' %self.inlets[1].get_average_speed()
        message += 'inlets[1]_average_stage [m]:  %.2f\n' %self.inlets[1].get_average_stage()
        message += 'inlets[1]_average_elevation [m]:  %.2f\n' %self.inlets[1].get_average_elevation()

        
        message += 'Discharge [m^3/s]: %.2f\n' % self.discharge
        message += 'Discharge_function_value [m^3/s]: %.2f\n' % self.discharge_abs_timemean
        message += 'Velocity  [m/s]: %.2f\n' % self.velocity
        message += 'Outlet Depth  [m]: %.2f\n' % self.outlet_depth
        message += 'Accumulated Flow [m^3]: %.2f\n' % self.accumulated_flow
        message += 'Inlet Driving Energy %.2f\n' % self.driving_energy
        message += 'Delta Total Energy %.2f\n' % self.delta_total_energy
        message += 'Control at this instant: %s\n' % self.case
        



        print(message)


    def set_logging(self, flag=True):

        self.logging = flag

        # If flag is true open file with mode = "w" to form a clean file for logging
        if self.logging:
            self.log_filename = self.domain.get_datadir() + '/' + self.label + '.log'
            log_to_file(self.log_filename, self.statistics(), mode='w')
            log_to_file(self.log_filename, 'time, discharge_instantaneous, discharge_abs_timemean, velocity, accumulated_flow, driving_energy_instantaneous, delta_total_energy_instantaneous')

            #log_to_file(self.log_filename, self.culvert_type)


    def timestepping_statistics(self):

        message  = '%.5f, ' % self.domain.get_time()
        message += '%.5f, ' % self.discharge
        message += '%.5f, ' % self.discharge_abs_timemean
        message += '%.5f, ' % self.velocity
        message += '%.5f, ' % self.accumulated_flow
        message += '%.5f, ' % self.driving_energy
        message += '%.5f' % self.delta_total_energy
       
        # Reset discharge_abs_timemean since last time this function was called
        # (FIXME: This assumes that the function is called only just after a
        # yield step)
        self.discharge_abs_timemean = 0.

        return message


    def get_inlets(self):
        
        return self.inlets
        
        
    def get_culvert_length(self):
        
        return self.culvert_length
    
    
    
    def get_culvert_slope(self):
        
        inlet0 = self.inlets[0]
        inlet1 = self.inlets[1]
        
        elev0 = inlet0.get_enquiry_invert_elevation()
        elev1 = inlet1.get_enquiry_invert_elevation()
        
        return (elev1-elev0)/self.get_culvert_length()
                          
                          
        
    def get_culvert_width(self):
        
        return self.width
        
        
    def get_culvert_diameter(self):
    
            return self.diameter
        
        
    def get_culvert_height(self):
    
        return self.height

    def get_culvert_z1(self):
    
        return self.z1 

    def get_culvert_z2(self):
    
        return self.z2

    def get_culvert_blockage(self):
		
        return self.blockage 

    def get_culvert_barrels(self):
		
        return self.barrels
                       
    def get_culvert_apron(self):

        return self.apron


    def get_master_proc(self):

        return 0



    #--------------------------------------------------------
    # Set of enquiry functions so that in the sequential and paralle case
    # we can get equiry info fron the master Proc
    #---------------------------------------------------------

    def get_enquiry_stages(self):

        enq0 = self.inlets[0].get_enquiry_stage()
        enq1 = self.inlets[1].get_enquiry_stage()

        return [enq0, enq1]


    def get_enquiry_depths(self):

        enq0 = self.inlets[0].get_enquiry_depth()
        enq1 = self.inlets[1].get_enquiry_depth()

        return [enq0, enq1]


    def get_enquiry_positions(self):

        enq0 = self.inlets[0].get_enquiry_position()
        enq1 = self.inlets[1].get_enquiry_position()

        return [enq0, enq1]


    def get_enquiry_xmoms(self):

        enq0 = self.inlets[0].get_enquiry_xmom()
        enq1 = self.inlets[1].get_enquiry_xmom()

        return [enq0, enq1]

    def get_enquiry_ymoms(self):

        enq0 = self.inlets[0].get_enquiry_ymom()
        enq1 = self.inlets[1].get_enquiry_ymom()

        return [enq0, enq1]


    def get_enquiry_elevations(self):

        enq0 = self.inlets[0].get_enquiry_elevation()
        enq1 = self.inlets[1].get_enquiry_elevation()

        return [enq0, enq1]



    def get_enquiry_water_depths(self):

        enq0 = self.inlets[0].get_enquiry_water_depth()
        enq1 = self.inlets[1].get_enquiry_water_depth()

        return [enq0, enq1]


    def get_enquiry_invert_elevations(self):

        enq0 = self.inlets[0].get_enquiry_invert_elevation()
        enq1 = self.inlets[1].get_enquiry_invert_elevation()

        return [enq0, enq1]


    def get_enquiry_velocitys(self):

        enq0 = self.inlets[0].get_enquiry_velocity()
        enq1 = self.inlets[1].get_enquiry_velocity()

        return [enq0, enq1]


    def get_enquiry_xvelocitys(self):

        enq0 = self.inlets[0].get_enquiry_xvelocity()
        enq1 = self.inlets[1].get_enquiry_xvelocity()

        return [enq0, enq1]

    def get_enquiry_yvelocitys(self):

        enq0 = self.inlets[0].get_enquiry_yvelocity()
        enq1 = self.inlets[1].get_enquiry_yvelocity()

        return [enq0, enq1]


    def get_enquiry_speeds(self):

        enq0 = self.inlets[0].get_enquiry_speed()
        enq1 = self.inlets[1].get_enquiry_speed()

        return [enq0, enq1]


    def get_enquiry_velocity_heads(self):

        enq0 = self.inlets[0].get_enquiry_velocity_head()
        enq1 = self.inlets[1].get_enquiry_velocity_head()

        return [enq0, enq1]


    def get_enquiry_total_energys(self):

        enq0 = self.inlets[0].get_enquiry_total_energy()
        enq1 = self.inlets[1].get_enquiry_total_energy()

        return [enq0, enq1]


    def get_enquiry_specific_energys(self):

        enq0 = self.inlets[0].get_enquiry_specific_energy()
        enq1 = self.inlets[1].get_enquiry_specific_energy()

        return [enq0, enq1]

