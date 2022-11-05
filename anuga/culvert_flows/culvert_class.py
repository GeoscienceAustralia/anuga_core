
import sys

from anuga.shallow_water.forcing import Inflow, General_forcing
from anuga.culvert_flows.culvert_polygons import create_culvert_polygons
from anuga.utilities.system_tools import log_to_file
from anuga.geometry.polygon import inside_polygon, is_inside_polygon, plot_polygons

from anuga.utilities.numerical_tools import mean
from anuga.utilities.numerical_tools import ensure_numeric, sign
        
from anuga.config import g, epsilon
from anuga.config import minimum_allowed_height, velocity_protection        
import anuga.utilities.log as log

import numpy as num
from math import sqrt
from math import sqrt

class Below_interval(Exception): pass 
class Above_interval(Exception): pass

# FIXME(Ole): Take a good hard look at logging here


# FIXME(Ole): Write in C and reuse this function by similar code 
# in interpolate.py
def interpolate_linearly(x, xvec, yvec):

    msg = 'Input to function interpolate_linearly could not be converted '
    msg += 'to numerical scalar: x = %s' % str(x)
    try:
        x = float(x)
    except:
        raise Exception(msg)


    # Check bounds
    if x < xvec[0]: 
        msg = 'Value provided = %.2f, interpolation minimum = %.2f.'\
            % (x, xvec[0])
        raise Below_interval(msg)
        
    if x > xvec[-1]: 
        msg = 'Value provided = %.2f, interpolation maximum = %.2f.'\
            %(x, xvec[-1])
        raise Above_interval(msg)
        
        
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
    

    

class Culvert_flow_general(object):
    """Culvert flow - transfer water from one hole to another
    
    This version will work with either rating curve file or with culvert
    routine.
    
    Input: Two points, pipe_size (either diameter or width, height), 
    mannings_rougness,
    """	

    def __init__(self,
                 domain,
                 culvert_description_filename=None,
                 culvert_routine=None,
                 end_point0=None, 
                 end_point1=None,
                 enquiry_point0=None, 
                 enquiry_point1=None,
                 type='box',
                 width=None,
                 height=None,
                 length=None,
                 number_of_barrels=1,
                 number_of_smoothing_steps=2000,
                 trigger_depth=0.01, # Depth below which no flow happens
                 manning=None,          # Mannings Roughness for Culvert 
                 sum_loss=None,
                 use_velocity_head=False, # FIXME(Ole): Get rid of - always True
                 use_momentum_jet=False, # FIXME(Ole): Not yet implemented
                 label=None,
                 description=None,
                 update_interval=None,
                 log_file=False,
                 discharge_hydrograph=False,
                 verbose=False):
        

        
        # Input check
        
        if height is None: height = width
        
        assert number_of_barrels >= 1
        assert use_velocity_head is True or use_velocity_head is False
        
        #msg = 'Momentum jet not yet moved to general culvert'
        #assert use_momentum_jet is False, msg
        self.use_momentum_jet = use_momentum_jet
 
        self.culvert_routine = culvert_routine        
        self.culvert_description_filename = culvert_description_filename
        if culvert_description_filename is not None:
            label, type, width, height, length, number_of_barrels, description, rating_curve = read_culvert_description(culvert_description_filename)
            self.rating_curve = ensure_numeric(rating_curve)            

        self.height = height
        self.width = width

        
        self.domain = domain
        self.trigger_depth = trigger_depth        
                
        if manning is None: 
            self.manning = 0.012   # Default roughness for pipe
        
        if sum_loss is None:
            self.sum_loss = 0.0
            
        
                        
        # Store culvert information
        self.label = label
        self.description = description
        self.culvert_type = type
        self.number_of_barrels = number_of_barrels
        
        # Store options        
        self.use_velocity_head = use_velocity_head

        if label is None: label = 'culvert_flow'
        label += '_' + str(id(self)) 
        self.label = label
        
        # File for storing discharge_hydrograph
        if discharge_hydrograph is True:
            self.timeseries_filename = label + '_timeseries.csv'
            fid = open(self.timeseries_filename, 'w')
            fid.write('time, discharge\n')
            fid.close()

        # Log file for storing general textual output
        if log_file is True:        
            self.log_filename = label + '.log'         
            log_to_file(self.log_filename, self.label)        
            log_to_file(self.log_filename, description)
            log_to_file(self.log_filename, self.culvert_type)        
        else:
            self.log_filename = None


        # Create the fundamental culvert polygons from polygon
        P = create_culvert_polygons(end_point0,
                                    end_point1,
                                    width=width,   
                                    height=height,
                                    number_of_barrels=number_of_barrels)
        self.culvert_polygons = P
        
        # Select enquiry points
        if enquiry_point0 is None:
            enquiry_point0 = P['enquiry_point0']
            
        if enquiry_point1 is None:
            enquiry_point1 = P['enquiry_point1']            
            
        if verbose is True:
            pass
            #plot_polygons([[end_point0, end_point1],
            #               P['exchange_polygon0'],
            #               P['exchange_polygon1'],
            #               [enquiry_point0, 1.005*enquiry_point0],
            #               [enquiry_point1, 1.005*enquiry_point1]],
            #              figname='culvert_polygon_output')

            
            
        self.enquiry_points = [enquiry_point0, enquiry_point1]
        self.enquiry_indices = self.get_enquiry_indices()                  
        self.check_culvert_inside_domain()
        
                    
        # Create inflow object at each end of the culvert. 
        self.openings = []
        self.openings.append(Inflow(domain,
                                    polygon=P['exchange_polygon0']))
        self.openings.append(Inflow(domain,
                                    polygon=P['exchange_polygon1']))
                                    
        # Assume two openings for now: Referred to as 0 and 1
        assert len(self.openings) == 2

        # Establish initial values at each enquiry point
        dq = domain.quantities 
        for i, opening in enumerate(self.openings):
            idx = self.enquiry_indices[i]
            elevation = dq['elevation'].get_values(location='centroids',
                                                   indices=[idx])[0]
            stage = dq['stage'].get_values(location='centroids',
                                           indices=[idx])[0]
            opening.elevation = elevation
            opening.stage = stage
            opening.depth = stage-elevation

            
                            
        # Determine initial pipe direction.
        # This may change dynamically based on the total energy difference     
        # Consequently, this may be superfluous
        delta_z = self.openings[0].elevation - self.openings[1].elevation
        if delta_z > 0.0:
            self.inlet = self.openings[0]
            self.outlet = self.openings[1]
        else:                
            self.outlet = self.openings[0]
            self.inlet = self.openings[1]        
        
        
        # Store basic geometry 
        self.end_points = [end_point0, end_point1]
        self.vector = P['vector']
        self.length = P['length']; assert self.length > 0.0
        if culvert_description_filename is not None:
            if not num.allclose(self.length, length, rtol=1.0e-2, atol=1.0e-2):
                msg = 'WARNING: barrel length specified in "%s" (%.2f m)'\
                    % (culvert_description_filename, 
                       length)
                msg += ' does not match distance between specified'
                msg += ' end points (%.2f m)' %self.length
                log.critical(msg)
        
        self.verbose = verbose

        # Circular index for flow averaging in culvert
        self.N = N = number_of_smoothing_steps
        self.Q_list = [0]*N
        self.i = i
        
        # For use with update_interval                        
        self.last_update = 0.0
        self.update_interval = update_interval
        
        # Create objects to update momentum (a bit crude at this stage). This is used with the momentum jet.
        xmom0 = General_forcing(domain, 'xmomentum',
                                polygon=P['exchange_polygon0'])

        xmom1 = General_forcing(domain, 'xmomentum',
                                polygon=P['exchange_polygon1'])

        ymom0 = General_forcing(domain, 'ymomentum',
                                polygon=P['exchange_polygon0'])

        ymom1 = General_forcing(domain, 'ymomentum',
                                polygon=P['exchange_polygon1'])

        self.opening_momentum = [[xmom0, ymom0], [xmom1, ymom1]]



        # Print some diagnostics to log if requested
        if self.log_filename is not None:
            s = 'Culvert Effective Length = %.2f m' %(self.length)
            log_to_file(self.log_filename, s)
   
            s = 'Culvert Direction is %s\n' %str(self.vector)
            log_to_file(self.log_filename, s)        
        
        
        
        
        
    def __call__(self, domain):

        # Time stuff
        time = domain.get_time()
        
        
        update = False
        if self.update_interval is None:
            # Use next timestep as has been computed in domain.py        
            delta_t = domain.timestep            
            update = True
        else:    
            # Use update interval 
            delta_t = self.update_interval            
            if time - self.last_update > self.update_interval or time == 0.0:
                update = True
            
        if self.log_filename is not None:        
            s = '\nTime = %.2f, delta_t = %f' %(time, delta_t)
            log_to_file(self.log_filename, s)
                
                                
        if update is True:
            self.compute_rates(delta_t)
            
    
        # Execute flow term for each opening
        # This is where Inflow objects are evaluated using the last rate 
        # that has been calculated
        # 
        # This will take place at every internal timestep and update the domain
        for opening in self.openings:
            opening(domain)
            


    def get_enquiry_indices(self):
        """Get indices for nearest centroids to self.enquiry_points
        """
        
        domain = self.domain
        
        enquiry_indices = []                  
        for point in self.enquiry_points:
            # Find nearest centroid 
            N = len(domain)    
            points = domain.get_centroid_coordinates(absolute=True)

            # Calculate indices in exchange area for this forcing term
            
            triangle_id = min_dist = sys.maxsize
            for k in range(N):
                x, y = points[k,:] # Centroid

                c = point                                
                distance = (x-c[0])**2+(y-c[1])**2
                if distance < min_dist:
                    min_dist = distance
                    triangle_id = k

                    
            if triangle_id < sys.maxsize:
                msg = 'found triangle with centroid (%f, %f)'\
                    %tuple(points[triangle_id, :])
                msg += ' for point (%f, %f)' %tuple(point)
                
                enquiry_indices.append(triangle_id)
            else:
                msg = 'Triangle not found for point (%f, %f)' %point
                raise Exception(msg)
        
        return enquiry_indices

        
    def check_culvert_inside_domain(self):
        """Check that all polygons and enquiry points lie within the mesh.
        """
        bounding_polygon = self.domain.get_boundary_polygon()
        P = self.culvert_polygons
        for key in list(P.keys()):
            if key in ['exchange_polygon0', 
                       'exchange_polygon1']:
                for point in list(P[key]) + self.enquiry_points:
                    msg = 'Point %s in polygon %s for culvert %s did not'\
                        %(str(point), key, self.label)
                    msg += 'fall within the domain boundary.'
                    assert is_inside_polygon(point, bounding_polygon), msg
            

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

        # Short hands
        domain = self.domain        
        dq = domain.quantities                
        
        # Time stuff
        time = domain.get_time()
        self.last_update = time

            
        if hasattr(self, 'log_filename'):
            log_filename = self.log_filename
            
        # Compute stage, energy and velocity at the 
        # enquiry points at each end of the culvert
        openings = self.openings
        for i, opening in enumerate(openings):
            idx = self.enquiry_indices[i]                
            
            stage = dq['stage'].get_values(location='centroids',
                                           indices=[idx])[0]
            depth = h = stage-opening.elevation
                                                           
            
            # Get velocity                                 
            xmomentum = dq['xmomentum'].get_values(location='centroids',
                                                   indices=[idx])[0]
            ymomentum = dq['xmomentum'].get_values(location='centroids',
                                                   indices=[idx])[0]

            if h > minimum_allowed_height:
                u = xmomentum/(h + velocity_protection/h)
                v = ymomentum/(h + velocity_protection/h)
            else:
                u = v = 0.0
                
            v_squared = u*u + v*v
            
            if self.use_velocity_head is True:
                velocity_head = 0.5*v_squared/g    
            else:
                velocity_head = 0.0
            
            opening.total_energy = velocity_head + stage
            opening.specific_energy = velocity_head + depth
            opening.stage = stage
            opening.depth = depth
            opening.velocity = sqrt(v_squared)
            

        # We now need to deal with each opening individually
        # Determine flow direction based on total energy difference
        delta_total_energy = openings[0].total_energy - openings[1].total_energy
        if delta_total_energy > 0:
            inlet = openings[0]
            outlet = openings[1]

            # FIXME: I think this whole momentum jet thing could be a bit more elegant
            inlet.momentum = self.opening_momentum[0]
            outlet.momentum = self.opening_momentum[1]
        else:
            inlet = openings[1]
            outlet = openings[0]
            
            inlet.momentum = self.opening_momentum[1]
            outlet.momentum = self.opening_momentum[0]

            delta_total_energy = -delta_total_energy

        self.inlet = inlet
        self.outlet = outlet
            
        msg = 'Total energy difference is negative'
        assert delta_total_energy >= 0.0, msg

        # Recompute slope and issue warning if flow is uphill
        # These values do not enter the computation
        delta_z = inlet.elevation - outlet.elevation
        culvert_slope = (delta_z/self.length)
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
        if inlet.specific_energy > delta_total_energy:
            # Outlet control
            driving_head = delta_total_energy
        else:
            # Inlet control
            driving_head = inlet.specific_energy
            


        if self.inlet.depth <= self.trigger_depth:
            Q = 0.0
        else:
            # Calculate discharge for one barrel and 
            # set inlet.rate and outlet.rate
            
            if self.culvert_description_filename is not None:
                try:
                    Q = interpolate_linearly(driving_head, 
                                             self.rating_curve[:,0], 
                                             self.rating_curve[:,1]) 
                except Below_interval as e:
                    Q = self.rating_curve[0,1]             
                    msg = '%.2fs: ' % time 
                    msg += 'Delta head smaller than rating curve minimum: '
                    msg += str(e)
                    msg += '\n        '
                    msg += 'I will use minimum discharge %.2f m^3/s ' % Q
                    msg += 'for culvert "%s"' % self.label
                    
                    if hasattr(self, 'log_filename'):                    
                        log_to_file(self.log_filename, msg)
                except Above_interval as e:
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
                    self.culvert_routine(inlet.depth,
                                         outlet.depth,
                                         inlet.velocity,
                                         outlet.velocity,
                                         inlet.specific_energy, 
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

            if int(domain.get_time()*100) % 100 == 0:

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


            



                            
# OBSOLETE (Except for momentum jet in Culvert_flow_energy)    
class Culvert_flow_rating(object):
    """Culvert flow - transfer water from one hole to another
    

    Input: Two points, pipe_size (either diameter or width, height), 
    mannings_rougness,
    inlet/outlet energy_loss_coefficients, internal_bend_coefficent,
    top-down_blockage_factor and bottom_up_blockage_factor
    
    """	

    def __init__(self,
                 domain,
                 culvert_description_filename=None,
                 end_point0=None, 
                 end_point1=None,
                 enquiry_point0=None, 
                 enquiry_point1=None,
                 update_interval=None,
                 log_file=False,
                 discharge_hydrograph=False,
                 verbose=False):
        

        
        label, type, width, height, length, number_of_barrels, description, rating_curve = read_culvert_description(culvert_description_filename)
        
                
        # Store culvert information
        self.label = label
        self.description = description
        self.culvert_type = type
        self.rating_curve = ensure_numeric(rating_curve)
        self.number_of_barrels = number_of_barrels

        if label is None: label = 'culvert_flow'
        label += '_' + str(id(self)) 
        self.label = label
        
        # File for storing discharge_hydrograph
        if discharge_hydrograph is True:
            self.timeseries_filename = label + '_timeseries.csv'
            fid = open(self.timeseries_filename, 'w')
            fid.write('time, discharge\n')
            fid.close()

        # Log file for storing general textual output
        if log_file is True:        
            self.log_filename = label + '.log'         
            log_to_file(self.log_filename, self.label)        
            log_to_file(self.log_filename, description)
            log_to_file(self.log_filename, self.culvert_type)        


        # Create the fundamental culvert polygons from POLYGON
        #if self.culvert_type == 'circle':
        #    # Redefine width and height for use with create_culvert_polygons
        #    width = height = diameter
        
        P = create_culvert_polygons(end_point0,
                                    end_point1,
                                    width=width,   
                                    height=height,
                                    number_of_barrels=number_of_barrels)
        
        # Select enquiry points
        if enquiry_point0 is None:
            enquiry_point0 = P['enquiry_point0']
            
        if enquiry_point1 is None:
            enquiry_point1 = P['enquiry_point1']            
            
        if verbose is True:
            pass
            #plot_polygons([[end_point0, end_point1],
            #               P['exchange_polygon0'],
            #               P['exchange_polygon1'],
            #               [enquiry_point0, 1.005*enquiry_point0],
            #               [enquiry_point1, 1.005*enquiry_point1]],                           
            #              figname='culvert_polygon_output')

            
            
        self.enquiry_points = [enquiry_point0, enquiry_point1]                           

        self.enquiry_indices = []                  
        for point in self.enquiry_points:
            # Find nearest centroid 
            N = len(domain)    
            points = domain.get_centroid_coordinates(absolute=True)

            # Calculate indices in exchange area for this forcing term
            
            triangle_id = min_dist = sys.maxsize
            for k in range(N):
                x, y = points[k,:] # Centroid

                c = point                                
                distance = (x-c[0])**2+(y-c[1])**2
                if distance < min_dist:
                    min_dist = distance
                    triangle_id = k

                    
            if triangle_id < sys.maxsize:
                msg = 'found triangle with centroid (%f, %f)'\
                    %tuple(points[triangle_id, :])
                msg += ' for point (%f, %f)' %tuple(point)
                
                self.enquiry_indices.append(triangle_id)
            else:
                msg = 'Triangle not found for point (%f, %f)' %point
                raise Exception(msg)
        
                          

        # Check that all polygons lie within the mesh.
        bounding_polygon = domain.get_boundary_polygon()
        for key in list(P.keys()):
            if key in ['exchange_polygon0', 
                       'exchange_polygon1']:
                for point in list(P[key]) + self.enquiry_points:
                    msg = 'Point %s in polygon %s for culvert %s did not'\
                        %(str(point), key, self.label)
                    msg += 'fall within the domain boundary.'
                    assert is_inside_polygon(point, bounding_polygon), msg
        
                    
        # Create inflow object at each end of the culvert. 
        self.openings = []
        self.openings.append(Inflow(domain,
                                    polygon=P['exchange_polygon0']))

        self.openings.append(Inflow(domain,
                                    polygon=P['exchange_polygon1']))                                    



        dq = domain.quantities                                            
        for i, opening in enumerate(self.openings):                            
            elevation = dq['elevation'].get_values(location='centroids',
                                                   indices=[self.enquiry_indices[i]])            
            opening.elevation = elevation
            opening.stage = elevation # Simple assumption that culvert is dry initially

        # Assume two openings for now: Referred to as 0 and 1
        assert len(self.openings) == 2
                            
        # Determine pipe direction     
        self.delta_z = delta_z = self.openings[0].elevation - self.openings[1].elevation
        if delta_z > 0.0:
            self.inlet = self.openings[0]
            self.outlet = self.openings[1]
        else:                
            self.outlet = self.openings[0]
            self.inlet = self.openings[1]        
        
        
        # Store basic geometry 
        self.end_points = [end_point0, end_point1]
        self.vector = P['vector']
        self.length = P['length']; assert self.length > 0.0
        if not num.allclose(self.length, length, rtol=1.0e-2, atol=1.0e-2):
            msg = 'WARNING: barrel length specified in "%s" (%.2f m)' %(culvert_description_filename, length)
            msg += ' does not match distance between specified'
            msg += ' end points (%.2f m)' %self.length
            log.critical(msg)
        
        self.verbose = verbose
        self.last_update = 0.0 # For use with update_interval        
        self.last_time = 0.0                
        self.update_interval = update_interval
        

        # Print something
        if hasattr(self, 'log_filename'):
            s = 'Culvert Effective Length = %.2f m' %(self.length)
            log_to_file(self.log_filename, s)
   
            s = 'Culvert Direction is %s\n' %str(self.vector)
            log_to_file(self.log_filename, s)        
        
        
        
        
        
    def __call__(self, domain):

        # Time stuff
        time = domain.get_time()
        
        
        update = False
        if self.update_interval is None:
            update = True
            delta_t = domain.timestep # Next timestep has been computed in domain.py
        else:    
            if time - self.last_update > self.update_interval or time == 0.0:
                update = True
            delta_t = self.update_interval
            
        s = '\nTime = %.2f, delta_t = %f' %(time, delta_t)
        if hasattr(self, 'log_filename'):            
            log_to_file(self.log_filename, s)
                
                                
        if update is True:
            self.last_update = time
       
            dq = domain.quantities
                        
            # Get average water depths at each opening        
            openings = self.openings   # There are two Opening [0] and [1]
            for i, opening in enumerate(openings):
                
                # Compute mean values of selected quantitites in the 
                # enquiry area in front of the culvert
                
                stage = dq['stage'].get_values(location='centroids',
                                               indices=[self.enquiry_indices[i]])
                
                # Store current average stage and depth with each opening object
                opening.depth = stage - opening.elevation
                opening.stage = stage

                

            #################  End of the FOR loop ################################################

            # We now need to deal with each opening individually
                
            # Determine flow direction based on total energy difference

            delta_w = self.inlet.stage - self.outlet.stage
            
            if hasattr(self, 'log_filename'):
                s = 'Time=%.2f, inlet stage = %.2f, outlet stage = %.2f' %(time, 
                                                                           self.inlet.stage,
                                                                           self.outlet.stage)
                log_to_file(self.log_filename, s)


            if self.inlet.depth <= 0.01:
                Q = 0.0
            else:
                # Calculate discharge for one barrel and set inlet.rate and outlet.rate
                
                try:
                    Q = interpolate_linearly(delta_w, self.rating_curve[:,0], self.rating_curve[:,1]) 
                except Below_interval as e:
                    Q = self.rating_curve[0,1]             
                    msg = '%.2fs: Delta head smaller than rating curve minimum: ' %time
                    msg += str(e)
                    msg += '\n        I will use minimum discharge %.2f m^3/s for culvert "%s"'\
                        %(Q, self.label)
                    if hasattr(self, 'log_filename'):                    
                        log_to_file(self.log_filename, msg)
                except Above_interval as e:
                    Q = self.rating_curve[-1,1]          
                    msg = '%.2fs: Delta head greater than rating curve maximum: ' %time
                    msg += str(e)
                    msg += '\n        I will use maximum discharge %.2f m^3/s for culvert "%s"'\
                        %(Q, self.label)
                    if hasattr(self, 'log_filename'):                    
                        log_to_file(self.log_filename, msg)                    

                
                
            
            # Adjust discharge for multiple barrels
            Q *= self.number_of_barrels
            

            # Adjust Q downwards depending on available water at inlet
            stage = self.inlet.get_quantity_values(quantity_name='stage')
            elevation = self.inlet.get_quantity_values(quantity_name='elevation')
            depth = stage-elevation
            
            
            V = 0
            for i, d in enumerate(depth):
                V += d * domain.areas[i]
            
            dt = delta_t            
            if Q*dt > V:
            
                Q_reduced = 0.9*V/dt # Reduce with safety factor
                
                msg = '%.2fs: Computed extraction for this time interval (Q*dt) is ' % time
                msg += 'greater than current volume (V) at inlet.\n'
                msg += ' Q will be reduced from %.2f m^3/s to %.2f m^3/s.' % (Q, Q_reduced)
                
                if self.verbose is True:
                    log.critical(msg)
                if hasattr(self, 'log_filename'):                    
                    log_to_file(self.log_filename, msg)                                        
                
                Q = Q_reduced
        
            self.inlet.rate = -Q
            self.outlet.rate = Q

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

            
    
        # Execute flow term for each opening
        # This is where Inflow objects are evaluated using the last rate that has been calculated
        # 
        # This will take place at every internal timestep and update the domain
        for opening in self.openings:
            opening(domain)



        
        

class Culvert_flow_energy(object):
    """Culvert flow - transfer water from one hole to another
    
    Using Momentum as Calculated by Culvert Flow !!
    Could be Several Methods Investigated to do This !!!

    2008_May_08
    To Ole:
    OK so here we need to get the Polygon Creating code to create 
    polygons for the culvert Based on
    the two input Points (X0,Y0) and (X1,Y1) - need to be passed 
    to create polygon

    The two centers are now passed on to create_polygon.
    

    Input: Two points, pipe_size (either diameter or width, height), 
    mannings_rougness,
    inlet/outlet energy_loss_coefficients, internal_bend_coefficent,
    top-down_blockage_factor and bottom_up_blockage_factor
    
    
    And the Delta H enquiry should be change from Openings in line 412 
    to the enquiry Polygons infront of the culvert
    At the moment this script uses only Depth, later we can change it to 
    Energy...

    Once we have Delta H can calculate a Flow Rate and from Flow Rate 
    an Outlet Velocity
    The Outlet Velocity x Outlet Depth = Momentum to be applied at the Outlet...
	
    Invert levels are optional. If left out they will default to the 
    elevation at the opening.
        
    """	

    def __init__(self,
                 domain,
                 label=None, 
                 description=None,
                 end_point0=None, 
                 end_point1=None,
                 width=None,
                 height=None,
                 diameter=None,
                 manning=None,          # Mannings Roughness for Culvert
                 invert_level0=None,    # Invert level at opening 0
                 invert_level1=None,    # Invert level at opening 1
                 loss_exit=None,
                 loss_entry=None,
                 loss_bend=None,
                 loss_special=None,
                 blockage_topdwn=None,
                 blockage_bottup=None,
                 culvert_routine=None,
                 number_of_barrels=1,
                 enquiry_point0=None, 
                 enquiry_point1=None,
                 update_interval=None,
                 verbose=False):
        
        # Input check
        if diameter is not None:
            self.culvert_type = 'circle'
            self.diameter = diameter
            if height is not None or width is not None:
                msg = 'Either diameter or width&height must be specified, '
                msg += 'but not both.'
                raise Exception(msg)
        else:
            if height is not None:
                if width is None:
                    self.culvert_type = 'square'                                
                    width = height
                else:
                    self.culvert_type = 'rectangle'
            elif width is not None:
                if height is None:
                    self.culvert_type = 'square'                                
                    height = width
            else:
                msg = 'Either diameter or width&height must be specified.'
                raise Exception(msg)
                
            if height == width:
                self.culvert_type = 'square'                                                
                
            self.height = height
            self.width = width

            
        assert self.culvert_type in ['circle', 'square', 'rectangle']
        
        assert number_of_barrels >= 1
        self.number_of_barrels = number_of_barrels
        
        
        # Set defaults
        if manning is None: manning = 0.012   # Default roughness for pipe
        if loss_exit is None: loss_exit = 1.00
        if loss_entry is None: loss_entry = 0.50
        if loss_bend is None: loss_bend=0.00
        if loss_special is None: loss_special=0.00
        if blockage_topdwn is None: blockage_topdwn=0.00
        if blockage_bottup is None: blockage_bottup=0.00
        if culvert_routine is None: 
            culvert_routine=boyd_generalised_culvert_model
            
        if label is None: label = 'culvert_flow'
        label += '_' + str(id(self)) 
        self.label = label
        
        # File for storing culvert quantities
        self.timeseries_filename = label + '_timeseries.csv'
        fid = open(self.timeseries_filename, 'w')
        fid.write('time, E0, E1, Velocity, Discharge\n')
        fid.close()

        # Log file for storing general textual output
        self.log_filename = label + '.log'         
        log_to_file(self.log_filename, self.label)        
        log_to_file(self.log_filename, description)
        log_to_file(self.log_filename, self.culvert_type)        


        # Create the fundamental culvert polygons from POLYGON
        if self.culvert_type == 'circle':
            # Redefine width and height for use with create_culvert_polygons
            width = height = diameter
        
        P = create_culvert_polygons(end_point0,
                                    end_point1,
                                    width=width,   
                                    height=height,
                                    number_of_barrels=number_of_barrels)
        
        # Select enquiry points
        if enquiry_point0 is None:
            enquiry_point0 = P['enquiry_point0']
            
        if enquiry_point1 is None:
            enquiry_point1 = P['enquiry_point1']            
            
        if verbose is True:
            pass
            #plot_polygons([[end_point0, end_point1],
            #               P['exchange_polygon0'],
            #               P['exchange_polygon1'],
            #               [enquiry_point0, 1.005*enquiry_point0],
            #               [enquiry_point1, 1.005*enquiry_point1]],
            #              figname='culvert_polygon_output')


        self.enquiry_points = [enquiry_point0, enquiry_point1]                           
        
        
        self.enquiry_indices = []                  
        for point in self.enquiry_points:
            # Find nearest centroid 
            N = len(domain)    
            points = domain.get_centroid_coordinates(absolute=True)

            # Calculate indices in exchange area for this forcing term
            
            triangle_id = min_dist = sys.maxsize
            for k in range(N):
                x, y = points[k,:] # Centroid

                c = point                                
                distance = (x-c[0])**2+(y-c[1])**2
                if distance < min_dist:
                    min_dist = distance
                    triangle_id = k

                    
            if triangle_id < sys.maxsize:
                msg = 'found triangle with centroid (%f, %f)'\
                    %tuple(points[triangle_id, :])
                msg += ' for point (%f, %f)' %tuple(point)
                
                self.enquiry_indices.append(triangle_id)
            else:
                msg = 'Triangle not found for point (%f, %f)' %point
                raise Exception(msg)
        
                          

            
            

        # Check that all polygons lie within the mesh.
        bounding_polygon = domain.get_boundary_polygon()
        for key in list(P.keys()):
            if key in ['exchange_polygon0', 
                       'exchange_polygon1']:
                for point in P[key]:
                
                    msg = 'Point %s in polygon %s for culvert %s did not'\
                        %(str(point), key, self.label)
                    msg += 'fall within the domain boundary.'
                    assert is_inside_polygon(point, bounding_polygon), msg
        

        # Create inflow object at each end of the culvert. 
        self.openings = []
        self.openings.append(Inflow(domain,
                                    polygon=P['exchange_polygon0']))

        self.openings.append(Inflow(domain,
                                    polygon=P['exchange_polygon1']))                                    


        # Assume two openings for now: Referred to as 0 and 1
        assert len(self.openings) == 2
        
        # Store basic geometry 
        self.end_points = [end_point0, end_point1]
        self.invert_levels = [invert_level0, invert_level1]                
        #self.enquiry_polygons = [P['enquiry_polygon0'], P['enquiry_polygon1']]
        #self.enquiry_polylines = [P['enquiry_polygon0'][:2], 
        #                          P['enquiry_polygon1'][:2]]
        self.vector = P['vector']
        self.length = P['length']; assert self.length > 0.0
        self.verbose = verbose
        self.last_time = 0.0        
        self.last_update = 0.0 # For use with update_interval        
        self.update_interval = update_interval
        

        # Store hydraulic parameters 
        self.manning = manning
        self.loss_exit = loss_exit
        self.loss_entry = loss_entry
        self.loss_bend = loss_bend
        self.loss_special = loss_special
        self.sum_loss = loss_exit + loss_entry + loss_bend + loss_special
        self.blockage_topdwn = blockage_topdwn
        self.blockage_bottup = blockage_bottup

        # Store culvert routine
        self.culvert_routine = culvert_routine

        
        # Create objects to update momentum (a bit crude at this stage)
        xmom0 = General_forcing(domain, 'xmomentum',
                                polygon=P['exchange_polygon0'])

        xmom1 = General_forcing(domain, 'xmomentum',
                                polygon=P['exchange_polygon1'])

        ymom0 = General_forcing(domain, 'ymomentum',
                                polygon=P['exchange_polygon0'])

        ymom1 = General_forcing(domain, 'ymomentum',
                                polygon=P['exchange_polygon1'])

        self.opening_momentum = [ [xmom0, ymom0], [xmom1, ymom1] ]
        

        # Print something
        s = 'Culvert Effective Length = %.2f m' %(self.length)
        log_to_file(self.log_filename, s)
   
        s = 'Culvert Direction is %s\n' %str(self.vector)
        log_to_file(self.log_filename, s)        
        
        
    def __call__(self, domain):

        log_filename = self.log_filename
         
        # Time stuff
        time = domain.get_time()
        
        # Short hand
        dq = domain.quantities
                

        update = False
        if self.update_interval is None:
            update = True
            delta_t = domain.timestep # Next timestep has been computed in domain.py
        else:    
            if time - self.last_update > self.update_interval or time == 0.0:
                update = True
            delta_t = self.update_interval
            
        s = '\nTime = %.2f, delta_t = %f' %(time, delta_t)
        if hasattr(self, 'log_filename'):            
            log_to_file(log_filename, s)
                
                                
        if update is True:
            self.last_update = time
                        
            msg = 'Time did not advance'
            if time > 0.0: assert delta_t > 0.0, msg


            # Get average water depths at each opening        
            openings = self.openings   # There are two Opening [0] and [1]
            for i, opening in enumerate(openings):
                
                # Compute mean values of selected quantitites in the 
                # exchange area in front of the culvert
                     
                stage = opening.get_quantity_values(quantity_name='stage')
                w = mean(stage) # Average stage

                # Use invert level instead of elevation if specified
                invert_level = self.invert_levels[i]
                if invert_level is not None:
                    z = invert_level
                else:
                    elevation = opening.get_quantity_values(quantity_name='elevation')
                    z = mean(elevation) # Average elevation

                # Estimated depth above the culvert inlet
                d = w - z  # Used for calculations involving depth
                if d < 0.0:
                    # This is possible since w and z are taken at different locations
                    #msg = 'D < 0.0: %f' %d
                    #raise Exception(msg)
                    d = 0.0
                

                # Ratio of depth to culvert height.
                # If ratio > 1 then culvert is running full
                if self.culvert_type == 'circle':
                    ratio = d/self.diameter
                else:    
                    ratio = d/self.height  
                opening.ratio = ratio
                    
                    
                # Average measures of energy in front of this opening
                
                id = [self.enquiry_indices[i]]
                stage = dq['stage'].get_values(location='centroids',
                                               indices=id)
                elevation = dq['elevation'].get_values(location='centroids',
                                                       indices=id)                                               
                xmomentum = dq['xmomentum'].get_values(location='centroids',
                                                       indices=id)                                               
                ymomentum = dq['xmomentum'].get_values(location='centroids',
                                                       indices=id)                                                                                              
                depth = stage-elevation
                if depth > 0.0:
                    u = xmomentum/(depth + velocity_protection/depth)
                    v = ymomentum/(depth + velocity_protection/depth)
                else:
                    u = v = 0.0

                    
                opening.total_energy = 0.5*(u*u + v*v)/g + stage

                # Store current average stage and depth with each opening object
                opening.depth = d
                opening.depth_trigger = d # Use this for now
                opening.stage = w
                opening.elevation = z
                

            #################  End of the FOR loop ################################################

            # We now need to deal with each opening individually
                
            # Determine flow direction based on total energy difference
            delta_Et = openings[0].total_energy - openings[1].total_energy

            if delta_Et > 0:
                inlet = openings[0]
                outlet = openings[1]

                inlet.momentum = self.opening_momentum[0]
                outlet.momentum = self.opening_momentum[1]

            else:
                inlet = openings[1]
                outlet = openings[0]

                inlet.momentum = self.opening_momentum[1]
                outlet.momentum = self.opening_momentum[0]
                
                delta_Et = -delta_Et

            self.inlet = inlet
            self.outlet = outlet
                
            msg = 'Total energy difference is negative'
            assert delta_Et > 0.0, msg

            delta_h = inlet.stage - outlet.stage
            delta_z = inlet.elevation - outlet.elevation
            culvert_slope = (delta_z/self.length)

            if culvert_slope < 0.0:
                # Adverse gradient - flow is running uphill
                # Flow will be purely controlled by uphill outlet face
                if self.verbose is True:
                    log.critical('WARNING: Flow is running uphill. Watch Out! '
                                 'inlet.elevation=%s, outlet.elevation%s'
                                 % (str(inlet.elevation), str(outlet.elevation)))


            s = 'Delta total energy = %.3f' %(delta_Et)
            log_to_file(log_filename, s)


            # Calculate discharge for one barrel and set inlet.rate and outlet.rate
            Q, barrel_velocity, culvert_outlet_depth = self.culvert_routine(self, inlet, outlet, delta_Et, g)
        
            # Adjust discharge for multiple barrels
            Q *= self.number_of_barrels

            # Compute barrel momentum
            barrel_momentum = barrel_velocity*culvert_outlet_depth
                    
            s = 'Barrel velocity = %f' %barrel_velocity
            log_to_file(log_filename, s)

            # Compute momentum vector at outlet
            outlet_mom_x, outlet_mom_y = self.vector * barrel_momentum
                
            s = 'Directional momentum = (%f, %f)' %(outlet_mom_x, outlet_mom_y)
            log_to_file(log_filename, s)

            # Log timeseries to file
            fid = open(self.timeseries_filename, 'a')        
            fid.write('%f, %f, %f, %f, %f\n'\
                          %(time, 
                            openings[0].total_energy,
                            openings[1].total_energy,
                            barrel_velocity,
                            Q))
            fid.close()

            # Update momentum        
            if delta_t > 0.0:
                xmomentum_rate = outlet_mom_x - outlet.momentum[0].value
                xmomentum_rate /= delta_t
                    
                ymomentum_rate = outlet_mom_y - outlet.momentum[1].value
                ymomentum_rate /= delta_t
                        
                s = 'X Y MOM_RATE = (%f, %f) ' %(xmomentum_rate, ymomentum_rate)
                log_to_file(log_filename, s)                    
            else:
                xmomentum_rate = ymomentum_rate = 0.0


            # Set momentum rates for outlet jet
            outlet.momentum[0].rate = xmomentum_rate
            outlet.momentum[1].rate = ymomentum_rate

            # Remember this value for next step (IMPORTANT)
            outlet.momentum[0].value = outlet_mom_x
            outlet.momentum[1].value = outlet_mom_y                    

            if int(domain.get_time()*100) % 100 == 0:
                s = 'T=%.5f, Culvert Discharge = %.3f f'\
                    %(time, Q)
                s += ' Depth= %0.3f  Momentum = (%0.3f, %0.3f)'\
                     %(culvert_outlet_depth, outlet_mom_x,outlet_mom_y)
                s += ' Momentum rate: (%.4f, %.4f)'\
                     %(xmomentum_rate, ymomentum_rate)                    
                s+='Outlet Vel= %.3f'\
                    %(barrel_velocity)
                log_to_file(log_filename, s)
            
            # Store value of time
            self.last_time = time
                


        # Execute flow term for each opening
        # This is where Inflow objects are evaluated and update the domain
        for opening in self.openings:
            opening(domain)

        # Execute momentum terms
        # This is where Inflow objects are evaluated and update the domain
        self.outlet.momentum[0](domain)
        self.outlet.momentum[1](domain)        
            


Culvert_flow = Culvert_flow_general        
