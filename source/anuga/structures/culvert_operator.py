import sys

from anuga.shallow_water.forcing import Inflow, General_forcing
from anuga.utilities.system_tools import log_to_file
from anuga.geometry.polygon import inside_polygon, is_inside_polygon
from anuga.geometry.polygon import plot_polygons, polygon_area

from anuga.utilities.numerical_tools import mean
from anuga.utilities.numerical_tools import ensure_numeric, sign
        
from anuga.config import g, epsilon
from anuga.config import minimum_allowed_height, velocity_protection        
import anuga.utilities.log as log

import Inlet

import numpy as num
import math

class Below_interval(Exception): pass 
class Above_interval(Exception): pass


class Generic_box_culvert:
    """Culvert flow - transfer water from one rectangular box to another.
    Sets up the geometry of problem
    
    This is the base class for culverts. Inherit from this class (and overwrite
    compute_discharge method for specific subclasses)
    
    Input: Two points, pipe_size (either diameter or width, height), 
    mannings_rougness,
    """	

    def __init__(self,
                 domain,
                 end_point0=None, 
                 end_point1=None,
                 enquiry_gap_factor=0.2,
                 width=None,
                 height=None,
                 verbose=False):
        
        # Input check
        
        self.domain = domain

        self.domain.set_fractional_step_operator(self)

        self.end_points = [end_point0, end_point1]
        self.enquiry_gap_factor = enquiry_gap_factor
        
        if height is None:
            height = width

        self.width = width
        self.height = height
        
        self.verbose=verbose
        self.filename = None
       
        # Create the fundamental culvert polygons and create inlet objects
        self.create_culvert_polygons()

        #FIXME (SR) Put this into a foe loop to deal with more inlets
        self.inlets = []
        polygon0 = self.inlet_polygons[0]
        enquiry_pt0 = self.enquiry_points[0]
        inlet0_vector = self.culvert_vector
        self.inlets.append(Inlet.Inlet(self.domain, polygon0, enquiry_pt0, inlet0_vector))

        polygon1 = self.inlet_polygons[1]
        enquiry_pt1 = self.enquiry_points[1]
        inlet1_vector = - self.culvert_vector
        self.inlets.append(Inlet.Inlet(self.domain, polygon1, enquiry_pt1, inlet1_vector))
 
        self.areas = self.domain.areas
        self.stages = self.domain.quantities['stage'].centroid_values
        self.elevations = self.domain.quantities['elevation'].centroid_values
        self.xmoms = self.domain.quantities['xmomentum'].centroid_values
        self.ymoms = self.domain.quantities['ymomentum'].centroid_values
   
        self.print_stats()


    def __call__(self):

        # Time stuff
        time     = self.domain.get_time()
        timestep = self.domain.get_timestep()

        inlet0 = self.inlets[0]
        inlet1 = self.inlets[1]

        # Aliases to cell indices
        inlet0_indices = inlet0.triangle_indices
        inlet1_indices = inlet1.triangle_indices

        # Inlet0 averages
        inlet0_heights = inlet0.get_heights()
        inlet0_areas = inlet0.get_areas()
        inlet0_water_volume = num.sum(inlet0_heights*inlet0_areas)
        average_inlet0_height = inlet0_water_volume/inlet0.area

        # Inlet1 averages
        inlet1_heights = inlet1.get_heights()
        inlet1_areas = inlet1.get_areas()
        inlet1_water_volume = num.sum(inlet1_heights*inlet1_areas)
        average_inlet1_height = inlet1_water_volume/inlet1.area

        # Transfer
        transfer_water = timestep*inlet0_water_volume
        
        self.stages.put(inlet0_indices, inlet0.get_elevations() + average_inlet0_height - transfer_water)
        self.xmoms.put(inlet0_indices, 0.0)
        self.ymoms.put(inlet0_indices, 0.0)

        self.stages.put(inlet1_indices, inlet1.get_elevations() + average_inlet1_height + transfer_water)
        self.xmoms.put(inlet1_indices, 0.0)
        self.ymoms.put(inlet1_indices, 0.0)


    def print_stats(self):

        print '====================================='
        print 'Generic Culvert Operator'
        print '====================================='
        print "enquiry_gap_factor"
        print self.enquiry_gap_factor
        
        for i, inlet in enumerate(self.inlets):
            print '-------------------------------------'
            print 'Inlet %i' % i
            print '-------------------------------------'

            print 'inlet triangle indices and centres'
            print inlet.triangle_indices[i]
            print self.domain.get_centroid_coordinates()[inlet.triangle_indices[i]]
        
            print 'polygon'
            print inlet.polygon

            print 'enquiry_point'
            print inlet.enquiry_point

        print '====================================='


    def set_store_hydrograph_discharge(self, filename=None):

        if filename is None:
            self.filename = 'culvert_discharge_hydrograph'
        else:
            self.filename = filename

        self.discharge_hydrograph = True
        
        self.timeseries_filename = self.filename + '_timeseries.csv'
        fid = open(self.timeseries_filename, 'w')
        fid.write('time, discharge\n')
        fid.close()


    def create_culvert_polygons(self):

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
        gap = (1 + self.enquiry_gap_factor)*h

        self.inlet_polygons = []
        self.enquiry_points = []

        # Build exchange polygon and enquiry points 0 and 1
        for i in [0, 1]:
            i0 = (2*i-1)
            p0 = self.end_points[i] + w
            p1 = self.end_points[i] - w
            p2 = p1 + i0*h
            p3 = p0 + i0*h
            self.inlet_polygons.append(num.array([p0, p1, p2, p3]))
            self.enquiry_points.append(self.end_points[i] + i0*gap)

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

            for j in [0,1]:
                point = self.enquiry_points[j]
                msg = 'Enquiry point falls inside a culvert polygon.'

                assert not inside_polygon(point, polygon), msg

    
    def adjust_flow_for_available_water_at_inlet(self, Q, delta_t):
        """Adjust Q downwards depending on available water at inlet

           This is a critical step in modelling bridges and Culverts
           the predicted flow through a structure based on an abstract
           algorithm can at times request for water that is simply not
           available due to any number of constrictions that limit the
           flow approaching the structure In order to ensure that
           there is adequate flow available certain checks are
           required There needs to be a check using the Static Water
           Volume sitting infront of the structure, In addition if the
           water is moving the available water will be larger than the
           static volume
           
           NOTE To temporarily switch this off for Debugging purposes
           rem out line in function def compute_rates below
        """
    
        if delta_t < epsilon:
            # No need to adjust if time step is very small or zero
            # In this case the possible flow will be very large
            # anyway.
            return Q
        
        # Short hands
        domain = self.domain        
        dq = domain.quantities                
        time = domain.get_time()
        I = self.inlet
        idx = I.exchange_indices    

        # Find triangle with the smallest depth
        stage = dq['stage'].get_values(location='centroids', 
                                               indices=[idx])
        elevation = dq['elevation'].get_values(location='centroids', 
                                               indices=[idx])        
        depth = stage-elevation
        min_depth = min(depth.flat)  # This may lead to errors if edge of area is at a higher level !!!!
        avg_depth = mean(depth.flat) # Yes, but this one violates the conservation unit tests



        # FIXME (Ole): If you want these, use log.critical() and
        # make the statements depend on verbose
        #print I.depth
        #print I.velocity
        #print self.width

        # max_Q Based on Volume Calcs


        depth_term = min_depth*I.exchange_area/delta_t
        if min_depth < 0.2:
            # Only add velocity term in shallow waters (< 20 cm)
            # This is a little ad hoc, but maybe it is reasonable
            velocity_term = self.width*min_depth*I.velocity
        else:
            velocity_term = 0.0

        # This one takes approaching water into account    
        max_Q = max(velocity_term, depth_term)

        # This one preserves Volume
        #max_Q = depth_term


        if self.verbose is True:
            log.critical('Max_Q = %f' % max_Q)            
            msg = 'Width = %.2fm, Depth at inlet = %.2f m, Velocity = %.2f m/s.      ' % (self.width, I.depth, I.velocity)
            msg += 'Max Q = %.2f m^3/s' %(max_Q)
            log.critical(msg)

        if self.log_filename is not None:                
            log_to_file(self.log_filename, msg)
        # New Procedure for assessing the flow available to the Culvert 
        # This routine uses the GET FLOW THROUGH CROSS SECTION
        #   Need to check Several Polyline however
        # Firstly 3 sides of the exchange Poly
        # then only the Line Directly infront of the Polygon
        # Access polygon Points from   self.inlet.polygon
      
        #  The Following computes the flow crossing over 3 sides of the exchange polygon for the structure
        # Clearly the flow in the culvert can not be more than that flowing toward it through the exhange polygon
        
        #q1 = domain.get_flow_through_cross_section(self.culvert_polygons['exchange_polygon0'][1:3]) # First Side Segment
        #q2 = domain.get_flow_through_cross_section(self.culvert_polygons['exchange_polygon0'][2:])   # Second Face Segment
        #q3 =domain.get_flow_through_cross_section(self.culvert_polygons['exchange_polygon0'].take([3,0], axis=0)) # Third Side Segment
        # q4 = domain.get_flow_through_cross_section([self.culvert_polygons['exchange_polygon0'][1:4]][0])
        #q4=max(q1,0.0)+max(q2,0.0)+max(q3,0.0)
        # To use only the Flow crossing the 3 sides of the Exchange Polygon use the following Line Only
        #max_Q=max(q1,q2,q3,q4)
        # Try Simple Smoothing using Average of 2 approaches
        #max_Q=(max(q1,q2,q3,q4)+max_Q)/2.0
        # Calculate the minimum in absolute terms of
        # the requsted flow and the possible flow
        Q_reduced = sign(Q)*min(abs(Q), abs(max_Q))
        if self.verbose is True:
            msg = 'Initial Q Reduced = %.2f m3/s.      ' % (Q_reduced)
            log.critical(msg)

        if self.log_filename is not None:                
            log_to_file(self.log_filename, msg)
        # Now Keep Rolling Average of Computed Discharge to Reduce / Remove Oscillations
        #  can use delta_t if we want to averageover a time frame for example
        # N = 5.0/delta_t  Will provide the average over 5 seconds

        self.i=(self.i+1)%self.N
        self.Q_list[self.i]=Q_reduced
        Q_reduced = sum(self.Q_list)/len(self.Q_list)

        if self.verbose is True:
            msg = 'Final Q Reduced = %.2f m3/s.      ' % (Q_reduced)
            log.critical(msg)

        if self.log_filename is not None:                
            log_to_file(self.log_filename, msg)


        if abs(Q_reduced) < abs(Q): 
            msg = '%.2fs: Requested flow is ' % time
            msg += 'greater than what is supported by the smallest '
            msg += 'depth at inlet exchange area:\n        '
            msg += 'inlet exchange area: %.2f '% (I.exchange_area) 
            msg += 'velocity at inlet :%.2f '% (I.velocity)
            msg += 'Vel* Exch Area = : %.2f '% (I.velocity*avg_depth*self.width)
            msg += 'h_min*inlet_area/delta_t = %.2f*%.2f/%.2f '\
                % (avg_depth, I.exchange_area, delta_t)
            msg += ' = %.2f m^3/s\n        ' % Q_reduced
            msg += 'Q will be reduced from %.2f m^3/s to %.2f m^3/s.' % (Q, Q_reduced)
            msg += 'Note calculate max_Q from V %.2f m^3/s ' % (max_Q)
            if self.verbose is True:
                log.critical(msg)
                
            if self.log_filename is not None:                
                log_to_file(self.log_filename, msg)
        
        return Q_reduced    
                        
            
    def compute_rates(self, delta_t):
        """Compute new rates for inlet and outlet
        """
        # Time stuff
        time = self.domain.get_time()
        self.last_update = time

        if hasattr(self, 'log_filename'):
            log_filename = self.log_filename
            
        # Compute stage, energy and velocity at the 
        # enquiry points at each end of the culvert
        total_energies = []
        specific_energies = []
        velocities = []
        
        for i, inlet in enumerate(self.inlets):
           
            stage = inlet.get_average_stage()
            depth = stage - inlet.get_average_elevation()
                                                                       
            if depth > minimum_allowed_height:
                u, v = inlet.get_average_velocities
            else:
                u = 0.0
                v = 0.0
                
            uv2 =  u**2 + v**2   
                
            if self.use_velocity_head is True:
                velocity_head = 0.5*uv2/g    
            else:
                velocity_head = 0.0
            
            inlet.set_total_energy(velocity_head + stage)
            inlet.set_specific_energy(velocity_head + depth)

        # We now need to deal with each opening individually
        # Determine flow direction based on total energy difference
        delta_total_energy = inlets[0].get_total_energy() - inlets[1].get_total_energy()
        
        if delta_total_energy >= 0:
            inlet = inlets[0]
            outlet = inlets[1]
        else:
            inlet = inlets[1]
            outlet = inlets[2]
            
            #msg = 'Total energy difference is negative'
            #assert delta_total_energy >= 0.0, msg

        # Recompute slope and issue warning if flow is uphill
        # These values do not enter the computation
        delta_z = inlet.get_average_elevation() - outlet.get_average_elevation()
        culvert_slope = (delta_z/self.culvert_length)
        if culvert_slope < 0.0:
            # Adverse gradient - flow is running uphill
            # Flow will be purely controlled by uphill outlet face
            if self.verbose is True:
                log.critical('%.2fs - WARNING: Flow is running uphill.' % time)
            
        if self.log_filename is not None:
            s = 'Time=%.2f, inlet stage = %.2f, outlet stage = %.2f'\
                %(time, self.inlet.stage, self.outlet.stage)
            log_to_file(self.log_filename, s)
            s = 'Delta total energy = %.3f' %(delta_total_energy)
            log_to_file(log_filename, s)
            
        # Determine controlling energy (driving head) for culvert
        if inlet.get_specific_energy() > delta_total_energy:
            # Outlet control
            driving_head = delta_total_energy
        else:
            # Inlet control
            driving_head = inlet.get_specific_energy()
            
        inlet_depth = inlet.get_average_stage() - inlet.get_average_elevation()
        outlet_depth = outlet.get_average_stage() - outlet.get_average_elevation()
        
        if inlet_depth <= self.trigger_depth:
            Q = 0.0
        else:
            # Calculate discharge for one barrel and 
            # set inlet.rate and outlet.rate
            
            if self.culvert_description_filename is not None:
                try:
                    Q = interpolate_linearly(driving_head, 
                                             self.rating_curve[:,0], 
                                             self.rating_curve[:,1]) 
                except Below_interval, e:
                    Q = self.rating_curve[0,1]             
                    msg = '%.2fs: ' % time 
                    msg += 'Delta head smaller than rating curve minimum: '
                    msg += str(e)
                    msg += '\n        '
                    msg += 'I will use minimum discharge %.2f m^3/s ' % Q
                    msg += 'for culvert "%s"' % self.label
                    
                    if hasattr(self, 'log_filename'):                    
                        log_to_file(self.log_filename, msg)
                except Above_interval, e:
                    Q = self.rating_curve[-1,1]          
                    msg = '%.2fs: ' % time                 
                    msg += 'Delta head greater than rating curve maximum: '
                    msg += str(e)
                    msg += '\n        '
                    msg += 'I will use maximum discharge %.2f m^3/s ' % Q
                    msg += 'for culvert "%s"' % self.label 
                    
                    if self.log_filename is not None:                    
                        log_to_file(self.log_filename, msg)
            else:
                # User culvert routine
                Q, barrel_velocity, culvert_outlet_depth =\
                    self.culvert_routine(inlet_depth,
                                         outlet_depth,
                                         inlet.get_average_velocity(),
                                         outlet.get_average_velocity(),
                                         inlet.get_specific_energy(), 
                                         delta_total_energy, 
                                         g,
                                         culvert_length=self.length,
                                         culvert_width=self.width,
                                         culvert_height=self.height,
                                         culvert_type=self.culvert_type,
                                         manning=self.manning,
                                         sum_loss=self.sum_loss,
                                         log_filename=self.log_filename)
            
        # Adjust discharge for multiple barrels
        Q *= self.number_of_barrels

        # Adjust discharge for available water at the inlet
        Q = self.adjust_flow_for_available_water_at_inlet(Q, delta_t)
        
        self.inlet.rate = -Q
        self.outlet.rate = Q

        # Momentum jet stuff
        if self.use_momentum_jet is True:

            # Compute barrel momentum
            barrel_momentum = barrel_velocity*culvert_outlet_depth

            if self.log_filename is not None:                                    
                s = 'Barrel velocity = %f' %barrel_velocity
                log_to_file(self.log_filename, s)

            # Compute momentum vector at outlet
            outlet_mom_x, outlet_mom_y = self.vector * barrel_momentum
                
            if self.log_filename is not None:                
                s = 'Directional momentum = (%f, %f)' %(outlet_mom_x, outlet_mom_y)
                log_to_file(self.log_filename, s)

            # Update momentum        
            if delta_t > 0.0:
                xmomentum_rate = outlet_mom_x - outlet.momentum[0].value
                xmomentum_rate /= delta_t
                    
                ymomentum_rate = outlet_mom_y - outlet.momentum[1].value
                ymomentum_rate /= delta_t
                        
                if self.log_filename is not None:                
                    s = 'X Y MOM_RATE = (%f, %f) ' %(xmomentum_rate, ymomentum_rate)
                    log_to_file(self.log_filename, s)                    
            else:
                xmomentum_rate = ymomentum_rate = 0.0


            # Set momentum rates for outlet jet
            outlet.momentum[0].rate = xmomentum_rate
            outlet.momentum[1].rate = ymomentum_rate

            # Remember this value for next step (IMPORTANT)
            outlet.momentum[0].value = outlet_mom_x
            outlet.momentum[1].value = outlet_mom_y                    

            if int(domain.time*100) % 100 == 0:

                if self.log_filename is not None:                
                    s = 'T=%.5f, Culvert Discharge = %.3f f'\
                        %(time, Q)
                    s += ' Depth= %0.3f  Momentum = (%0.3f, %0.3f)'\
                        %(culvert_outlet_depth, outlet_mom_x,outlet_mom_y)
                    s += ' Momentum rate: (%.4f, %.4f)'\
                        %(xmomentum_rate, ymomentum_rate)                    
                    s+='Outlet Vel= %.3f'\
                        %(barrel_velocity)
                    log_to_file(self.log_filename, s)


            # Execute momentum terms
            # This is where Inflow objects are evaluated and update the domain
                self.outlet.momentum[0](domain)
                self.outlet.momentum[1](domain)        
            

            
        # Log timeseries to file
        try:
            fid = open(self.timeseries_filename, 'a')        
        except:
            pass
        else:    
            fid.write('%.2f, %.2f\n' %(time, Q))
            fid.close()

        # Store value of time
        self.last_time = time


# FIXME(Ole): Write in C and reuse this function by similar code
# in interpolate.py
def interpolate_linearly(x, xvec, yvec):

    msg = 'Input to function interpolate_linearly could not be converted '
    msg += 'to numerical scalar: x = %s' % str(x)
    try:
        x = float(x)
    except:
        raise Exception, msg


    # Check bounds
    if x < xvec[0]:
        msg = 'Value provided = %.2f, interpolation minimum = %.2f.'\
            % (x, xvec[0])
        raise Below_interval, msg

    if x > xvec[-1]:
        msg = 'Value provided = %.2f, interpolation maximum = %.2f.'\
            %(x, xvec[-1])
        raise Above_interval, msg


    # Find appropriate slot within bounds
    i = 0
    while x > xvec[i]: i += 1


    x0 = xvec[i-1]
    x1 = xvec[i]
    alpha = (x - x0)/(x1 - x0)

    y0 = yvec[i-1]
    y1 = yvec[i]
    y = alpha*y1 + (1-alpha)*y0

    return y



def read_culvert_description(culvert_description_filename):

    # Read description file
    fid = open(culvert_description_filename)

    read_rating_curve_data = False
    rating_curve = []
    for i, line in enumerate(fid.readlines()):

        if read_rating_curve_data is True:
            fields = line.split(',')
            head_difference = float(fields[0].strip())
            flow_rate = float(fields[1].strip())
            barrel_velocity = float(fields[2].strip())

            rating_curve.append([head_difference, flow_rate, barrel_velocity])

        if i == 0:
            # Header
            continue
        if i == 1:
            # Metadata
            fields = line.split(',')
            label=fields[0].strip()
            type=fields[1].strip().lower()
            assert type in ['box', 'pipe']

            width=float(fields[2].strip())
            height=float(fields[3].strip())
            length=float(fields[4].strip())
            number_of_barrels=int(fields[5].strip())
            #fields[6] refers to losses
            description=fields[7].strip()

        if line.strip() == '': continue # Skip blanks

        if line.startswith('Rating'):
            read_rating_curve_data = True
            # Flow data follows

    fid.close()

    return label, type, width, height, length, number_of_barrels, description, rating_curve

