"""Collection of culvert routines for use with Culvert_flow in culvert_class

Usage:

   

"""

#NOTE:
# Inlet control:  Delta_total_energy > inlet_specific_energy
# Outlet control: Delta_total_energy < inlet_specific_energy
# where total energy is (w + 0.5*v^2/g) and
# specific energy is (h + 0.5*v^2/g)


from math import pi, sqrt, sin, cos


def boyd_generalised_culvert_model(inlet_depth, 
                                   outlet_depth,
                                   inlet_specific_energy, 
                                   delta_total_energy, 
                                   g,
                                   culvert_length=0.0,
                                   culvert_width=0.0,
                                   culvert_height=0.0,
                                   culvert_type='box',
                                   manning=0.0,
                                   sum_loss=0.0,
                                   log_filename=None):

    """Boyd's generalisation of the US department of transportation culvert 
    model
      
    The quantity of flow passing through a culvert is controlled by many factors
    It could be that the culvert is controlled by the inlet only, with it 
    being unsubmerged this is effectively equivalent to the weir Equation
    Else the culvert could be controlled by the inlet, with it being 
    submerged, this is effectively the Orifice Equation
    Else it may be controlled by down stream conditions where depending on 
    the down stream depth, the momentum in the culvert etc. flow is controlled 
    """

    from anuga.utilities.system_tools import log_to_file
    from anuga.config import velocity_protection
    from anuga.utilities.numerical_tools import safe_acos as acos


    if inlet_depth > 0.01:
        # Water has risen above inlet
        
        if log_filename is not None:                            
            s = 'Specific energy  = %f m' % inlet_specific_energy
            log_to_file(log_filename, s)
        
        msg = 'Specific energy at inlet is negative'
        assert inlet_specific_energy >= 0.0, msg
                      

        if culvert_type == 'circle':
            # Round culvert (use width as diameter)
            diameter = culvert_width            

            # Calculate flows for inlet control            
            Q_inlet_unsubmerged = 0.421*g**0.5*diameter**0.87*inlet_specific_energy**1.63 # Inlet Ctrl Inlet Unsubmerged 
            Q_inlet_submerged = 0.530*g**0.5*diameter**1.87*inlet_specific_energy**0.63   # Inlet Ctrl Inlet Submerged

            if log_filename is not None:                                        
                s = 'Q_inlet_unsubmerged = %.6f, Q_inlet_submerged = %.6f' % (Q_inlet_unsubmerged, Q_inlet_submerged)
                log_to_file(log_filename, s)


            # FIXME(Ole): Are these functions really for inlet control?                    
            if Q_inlet_unsubmerged < Q_inlet_submerged:
                Q = Q_inlet_unsubmerged
                alpha = acos(1 - inlet_depth/diameter)
                flow_area = diameter**2 * (alpha - sin(alpha)*cos(alpha))
                outlet_culvert_depth = inlet_depth
                case = 'Inlet unsubmerged'
            else:    
                Q = Q_inlet_submerged
                flow_area = (diameter/2)**2 * pi
                outlet_culvert_depth = diameter
                case = 'Inlet submerged'                    



            if delta_total_energy < inlet_specific_energy:
                # Calculate flows for outlet control
                
                # Determine the depth at the outlet relative to the depth of flow in the Culvert
                if outlet_depth > diameter:        # The Outlet is Submerged
                    outlet_culvert_depth=diameter
                    flow_area = (diameter/2)**2 * pi  # Cross sectional area of flow in the culvert
                    perimeter = diameter * pi
                    case = 'Outlet submerged'
                elif outlet_depth==0.0: 
                    outlet_culvert_depth=inlet_depth # For purpose of calculation assume the outlet depth = the inlet depth
                    alpha = acos(1 - inlet_depth/diameter)
                    flow_area = diameter**2 * (alpha - sin(alpha)*cos(alpha))
                    perimeter = alpha*diameter

                    case = 'Outlet depth is zero'                        
                else:   # Here really should use the Culvert Slope to calculate Actual Culvert Depth & Velocity
                    outlet_culvert_depth=outlet_depth   # For purpose of calculation assume the outlet depth = the inlet depth
                    alpha = acos(1 - outlet_depth/diameter)
                    flow_area = diameter**2 * (alpha - sin(alpha)*cos(alpha))
                    perimeter = alpha*diameter
                    case = 'Outlet is open channel flow'

                hyd_rad = flow_area/perimeter
                
                if log_filename is not None:
                    s = 'hydraulic radius at outlet = %f' %hyd_rad
                    log_to_file(log_filename, s)

                # Outlet control velocity using tail water
                culvert_velocity = sqrt(delta_total_energy/((sum_loss/2*g)+(manning**2*culvert_length)/hyd_rad**1.33333)) 
                Q_outlet_tailwater = flow_area * culvert_velocity
                
                if log_filename is not None:                
                    s = 'Q_outlet_tailwater = %.6f' %Q_outlet_tailwater
                    log_to_file(log_filename, s)
                    
                Q = min(Q, Q_outlet_tailwater)
            else:
                pass
                #FIXME(Ole): What about inlet control?


        else:
            # Box culvert (rectangle or square)

            # Calculate flows for inlet control
            height = culvert_height
            width = culvert_width            
            
            Q_inlet_unsubmerged = 0.540*g**0.5*width*inlet_specific_energy**1.50 # Flow based on Inlet Ctrl Inlet Unsubmerged
            Q_inlet_submerged = 0.702*g**0.5*width*height**0.89*inlet_specific_energy**0.61  # Flow based on Inlet Ctrl Inlet Submerged

            if log_filename is not None:                            
                s = 'Q_inlet_unsubmerged = %.6f, Q_inlet_submerged = %.6f' %(Q_inlet_unsubmerged, Q_inlet_submerged)
                log_to_file(log_filename, s)


            # FIXME(Ole): Are these functions really for inlet control?    
            if Q_inlet_unsubmerged < Q_inlet_submerged:
                Q = Q_inlet_unsubmerged
                flow_area = width*inlet_depth
                outlet_culvert_depth = inlet_depth
                case = 'Inlet unsubmerged'
            else:    
                Q = Q_inlet_submerged
                flow_area = width*height
                outlet_culvert_depth = height
                case = 'Inlet submerged'                    

            if delta_total_energy < inlet_specific_energy:
                # Calculate flows for outlet control
                
                # Determine the depth at the outlet relative to the depth of flow in the Culvert
                if outlet_depth > height:        # The Outlet is Submerged
                    outlet_culvert_depth=height
                    flow_area=width*height       # Cross sectional area of flow in the culvert
                    perimeter=2.0*(width+height)
                    case = 'Outlet submerged'
                elif outlet_depth==0.0: 
                    outlet_culvert_depth=inlet_depth   # For purpose of calculation assume the outlet depth = the inlet depth
                    flow_area=width*inlet_depth
                    perimeter=(width+2.0*inlet_depth)
                    case = 'Outlet depth is zero'                        
                else:   # Here really should use the Culvert Slope to calculate Actual Culvert Depth & Velocity
                    outlet_culvert_depth=outlet_depth
                    flow_area=width*outlet_depth
                    perimeter=(width+2.0*outlet_depth)
                    case = 'Outlet is open channel flow'

                hyd_rad = flow_area/perimeter
                
                if log_filename is not None:                                
                    s = 'hydraulic radius at outlet = %f' % hyd_rad
                    log_to_file(log_filename, s)

                # Outlet control velocity using tail water
                culvert_velocity = sqrt(delta_total_energy/((sum_loss/2*g)+(manning**2*culvert_length)/hyd_rad**1.33333)) 
                Q_outlet_tailwater = flow_area * culvert_velocity

                if log_filename is not None:                            
                    s = 'Q_outlet_tailwater = %.6f' % Q_outlet_tailwater
                    log_to_file(log_filename, s)
                Q = min(Q, Q_outlet_tailwater)
            else:
                pass
                #FIXME(Ole): What about inlet control?

                
        # Common code for circle and square geometries
        if log_filename is not None:
            log_to_file(log_filename, 'Case: "%s"' % case)
            
        if log_filename is not None:                        
            s = 'Flow Rate Control = %f' % Q
            log_to_file(log_filename, s)

            
        culv_froude=sqrt(Q**2*width/(g*flow_area**3))
        
        if log_filename is not None:                            
            s = 'Froude in Culvert = %f' % culv_froude
            log_to_file(log_filename, s)

        # Determine momentum at the outlet 
        barrel_velocity = Q/(flow_area + velocity_protection/flow_area)


    else: # inlet_depth < 0.01:
        Q = barrel_velocity = outlet_culvert_depth = 0.0
        
    return Q, barrel_velocity, outlet_culvert_depth 


