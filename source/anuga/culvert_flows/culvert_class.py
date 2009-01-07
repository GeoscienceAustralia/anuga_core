from anuga.shallow_water.shallow_water_domain import Inflow, General_forcing
from anuga.culvert_flows.culvert_polygons import create_culvert_polygons
from anuga.utilities.system_tools import log_to_file
from anuga.utilities.polygon import inside_polygon
from anuga.utilities.polygon import is_inside_polygon
from anuga.utilities.polygon import plot_polygons

from anuga.utilities.numerical_tools import mean
from anuga.utilities.numerical_tools import ensure_numeric
        
from anuga.config import g, epsilon
from Numeric import take, sqrt
from anuga.config import velocity_protection        

        


from Numeric import allclose
from Numeric import sqrt, sum



import sys

class Below_interval(Exception): pass 
class Above_interval(Exception): pass



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
    

class Culvert_flow_rating:
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
            
            triangle_id = min_dist = sys.maxint
            for k in range(N):
                x, y = points[k,:] # Centroid

                c = point                                
                distance = (x-c[0])**2+(y-c[1])**2
                if distance < min_dist:
                    min_dist = distance
                    triangle_id = k

                    
            if triangle_id < sys.maxint:
                msg = 'found triangle with centroid (%f, %f)'\
                    %tuple(points[triangle_id, :])
                msg += ' for point (%f, %f)' %tuple(point)
                
                self.enquiry_indices.append(triangle_id)
            else:
                msg = 'Triangle not found for point (%f, %f)' %point
                raise Exception, msg
        
                          

        # Check that all polygons lie within the mesh.
        bounding_polygon = domain.get_boundary_polygon()
        for key in P.keys():
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
        if not allclose(self.length, length, rtol=1.0e-2, atol=1.0e-2):
            msg = 'WARNING: barrel length specified in "%s" (%.2f m)' %(culvert_description_filename, length)
            msg += ' does not match distance between specified'
            msg += ' end points (%.2f m)' %self.length
            print msg
        
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
                except Below_interval, e:
                    Q = self.rating_curve[0,1]             
                    msg = '%.2fs: Delta head smaller than rating curve minimum: ' %time
                    msg += str(e)
                    msg += '\n        I will use minimum discharge %.2f m^3/s for culvert "%s"'\
                        %(Q, self.label)
                    if hasattr(self, 'log_filename'):                    
                        log_to_file(self.log_filename, msg)
                except Above_interval, e:
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
            
            #Vsimple = mean(depth)*self.inlet.exchange_area # Current volume in exchange area  
            #print 'Q', Q, 'dt', delta_t, 'Q*dt', Q*delta_t, 'V', V, 'Vsimple', Vsimple

            dt = delta_t            
            if Q*dt > V:
            
                Q_reduced = 0.9*V/dt # Reduce with safety factor
                
                msg = '%.2fs: Computed extraction for this time interval (Q*dt) is ' % time
                msg += 'greater than current volume (V) at inlet.\n'
                msg += ' Q will be reduced from %.2f m^3/s to %.2f m^3/s.' % (Q, Q_reduced)
                
                #print msg
                
                if self.verbose is True:
                    print msg
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



        
        

class Culvert_flow_energy:
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
        
        from Numeric import sqrt, sum

        # Input check
        if diameter is not None:
            self.culvert_type = 'circle'
            self.diameter = diameter
            if height is not None or width is not None:
                msg = 'Either diameter or width&height must be specified, '
                msg += 'but not both.'
                raise Exception, msg
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
                raise Exception, msg                
                
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
            
            triangle_id = min_dist = sys.maxint
            for k in range(N):
                x, y = points[k,:] # Centroid

                c = point                                
                distance = (x-c[0])**2+(y-c[1])**2
                if distance < min_dist:
                    min_dist = distance
                    triangle_id = k

                    
            if triangle_id < sys.maxint:
                msg = 'found triangle with centroid (%f, %f)'\
                    %tuple(points[triangle_id, :])
                msg += ' for point (%f, %f)' %tuple(point)
                
                self.enquiry_indices.append(triangle_id)
            else:
                msg = 'Triangle not found for point (%f, %f)' %point
                raise Exception, msg
        
                          

            
            

        # Check that all polygons lie within the mesh.
        bounding_polygon = domain.get_boundary_polygon()
        for key in P.keys():
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
                    #raise msg
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
                #print 'Et = %.3f m' %opening.total_energy

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
                #print 'Flow U/S ---> D/S'
                inlet = openings[0]
                outlet = openings[1]

                inlet.momentum = self.opening_momentum[0]
                outlet.momentum = self.opening_momentum[1]

            else:
                #print 'Flow D/S ---> U/S'
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
                    print 'WARNING: Flow is running uphill. Watch Out!', inlet.elevation, outlet.elevation


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

            if int(domain.time*100) % 100 == 0:
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
            


Culvert_flow = Culvert_flow_rating        
