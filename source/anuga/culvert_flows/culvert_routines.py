"""Collection of culvert routines for use with Culvert_flow in culvert_class

Usage:

   

"""

#NOTE:
# Inlet control:  Delta_Et > Es at the inlet
# Outlet control: Delta_Et < Es at the inlet
# where Et is total energy (w + 0.5*v^2/g) and
# Es is the specific energy (h + 0.5*v^2/g)



#   NEW DEFINED CULVERT FLOW---- Flow from INLET 1 ------> INLET 2 (Outlet)
# 
# The First Attempt has a Simple Horizontal Circle as a Hole on the Bed
# Flow Is Removed at a Rate of INFLOW 
# Downstream there is a similar Circular Hole on the Bed where INFLOW effectively Surcharges
#
# This SHould be changed to a Vertical Opening Both BOX and Circular
# There will be several Culvert Routines such as:
# CULVERT_Boyd_Channel
# CULVERT_Orifice_and_Weir
# CULVERT_Simple_FLOOR
# CULVERT_Simple_WALL
# CULVERT_Eqn_Floor
# CULVERT_Eqn_Wall
# CULVERT_Tab_Floor
# CULVERT_Tab_Wall
# BRIDGES.....
# NOTE NEED TO DEVELOP CONCEPT 1D Model for Linked Pipe System !!!!

#  COULD USE EPA SWMM Model !!!!


from math import pi, sqrt, sin, cos


def boyd_generalised_culvert_model(culvert, delta_Et, g):

    """Boyd's generalisation of the US department of transportation culvert model
        # == The quantity of flow passing through a culvert is controlled by many factors
        # == It could be that the culvert is controled by the inlet only, with it being Un submerged this is effectively equivalent to the WEIR Equation
        # == Else the culvert could be controlled by the inlet, with it being Submerged, this is effectively the Orifice Equation
        # == Else it may be controlled by Down stream conditions where depending on the down stream depth, the momentum in the culvert etc. flow is controlled 
    """

    from anuga.utilities.system_tools import log_to_file
    from anuga.config import velocity_protection
    from anuga.utilities.numerical_tools import safe_acos as acos

    inlet = culvert.inlet
    outlet = culvert.outlet        
    Q_outlet_tailwater = 0.0
    inlet.rate = 0.0
    outlet.rate = 0.0
    Q_inlet_unsubmerged = 0.0
    Q_inlet_submerged = 0.0
    Q_outlet_critical_depth = 0.0

    log_filename = culvert.log_filename

    manning = culvert.manning
    sum_loss = culvert.sum_loss
    length = culvert.length

    if inlet.depth_trigger >= 0.01 and inlet.depth >= 0.01:
        # Calculate driving energy
        # FIXME(Ole): Should this be specific energy?
        E = inlet.total_energy

        s = 'Driving energy  = %f m' %E
        log_to_file(log_filename, s)
        
        msg = 'Driving energy is negative'
        assert E >= 0.0, msg
                      
        
        # Water has risen above inlet
        if culvert.culvert_type == 'circle':
            # Round culvert

            # Calculate flows for inlet control            
            diameter = culvert.diameter

            Q_inlet_unsubmerged = 0.421*g**0.5*diameter**0.87*E**1.63    # Inlet Ctrl Inlet Unsubmerged 
            Q_inlet_submerged = 0.530*g**0.5*diameter**1.87*E**0.63      # Inlet Ctrl Inlet Submerged

            s = 'Q_inlet_unsubmerged = %.6f, Q_inlet_submerged = %.6f' %(Q_inlet_unsubmerged, Q_inlet_submerged)
            log_to_file(log_filename, s)

            case = ''
            if Q_inlet_unsubmerged < Q_inlet_submerged:
                Q = Q_inlet_unsubmerged

                alpha = acos(1 - inlet.depth/diameter)
                flow_area = diameter**2 * (alpha - sin(alpha)*cos(alpha))
                outlet_culvert_depth = inlet.depth
                width = diameter*sin(alpha)
                #perimeter = alpha*diameter                
                case = 'Inlet unsubmerged'
            else:    
                Q = Q_inlet_submerged
                flow_area = (diameter/2)**2 * pi
                outlet_culvert_depth = diameter
                width = diameter
                #perimeter = diameter
                case = 'Inlet submerged'                    



            if delta_Et < E:
                # Calculate flows for outlet control
                # Determine the depth at the outlet relative to the depth of flow in the Culvert

                if outlet.depth > diameter:        # The Outlet is Submerged
                    outlet_culvert_depth=diameter
                    flow_area = (diameter/2)**2 * pi  # Cross sectional area of flow in the culvert
                    perimeter = diameter * pi
                    width = diameter
                    case = 'Outlet submerged'
                elif outlet.depth==0.0: 
                    outlet_culvert_depth=inlet.depth   # For purpose of calculation assume the outlet depth = the inlet depth
                    alpha = acos(1 - inlet.depth/diameter)
                    flow_area = diameter**2 * (alpha - sin(alpha)*cos(alpha))
                    perimeter = alpha*diameter
                    width = diameter*sin(alpha)                    

                    case = 'Outlet depth is zero'                        
                else:   # Here really should use the Culvert Slope to calculate Actual Culvert Depth & Velocity
                    outlet_culvert_depth=outlet.depth   # For purpose of calculation assume the outlet depth = the inlet depth
                    alpha = acos(1 - outlet.depth/diameter)
                    flow_area = diameter**2 * (alpha - sin(alpha)*cos(alpha))
                    perimeter = alpha*diameter
                    width = diameter*sin(alpha)                    
                    case = 'Outlet is open channel flow'

                hyd_rad = flow_area/perimeter
                s = 'hydraulic radius at outlet = %f' %hyd_rad
                log_to_file(log_filename, s)

                # Outlet control velocity using tail water
                culvert_velocity = sqrt(delta_Et/((sum_loss/2*g)+(manning**2*length)/hyd_rad**1.33333)) 
                Q_outlet_tailwater = flow_area * culvert_velocity

                s = 'Q_outlet_tailwater = %.6f' %Q_outlet_tailwater
                log_to_file(log_filename, s)
                Q = min(Q, Q_outlet_tailwater)
                    


        else:
            # Box culvert (rectangle or square)

            # Calculate flows for inlet control
            height = culvert.height
            width = culvert.width            
            
            Q_inlet_unsubmerged = 0.540*g**0.5*width*E**1.50 # Flow based on Inlet Ctrl Inlet Unsubmerged
            Q_inlet_submerged = 0.702*g**0.5*width*height**0.89*E**0.61  # Flow based on Inlet Ctrl Inlet Submerged

            s = 'Q_inlet_unsubmerged = %.6f, Q_inlet_submerged = %.6f' %(Q_inlet_unsubmerged, Q_inlet_submerged)
            log_to_file(log_filename, s)

            case = ''
            if Q_inlet_unsubmerged < Q_inlet_submerged:
                Q = Q_inlet_unsubmerged
                flow_area = width*inlet.depth
                outlet_culvert_depth = inlet.depth
                #perimeter=(width+2.0*inlet.depth)                
                case = 'Inlet unsubmerged'
            else:    
                Q = Q_inlet_submerged
                flow_area = width*height
                outlet_culvert_depth = height
                #perimeter=2.0*(width+height)                
                case = 'Inlet submerged'                    

            if delta_Et < E:
                # Calculate flows for outlet control
                # Determine the depth at the outlet relative to the depth of flow in the Culvert

                if outlet.depth > height:        # The Outlet is Submerged
                    outlet_culvert_depth=height
                    flow_area=width*height       # Cross sectional area of flow in the culvert
                    perimeter=2.0*(width+height)
                    case = 'Outlet submerged'
                elif outlet.depth==0.0: 
                    outlet_culvert_depth=inlet.depth   # For purpose of calculation assume the outlet depth = the inlet depth
                    flow_area=width*inlet.depth
                    perimeter=(width+2.0*inlet.depth)
                    case = 'Outlet depth is zero'                        
                else:   # Here really should use the Culvert Slope to calculate Actual Culvert Depth & Velocity
                    outlet_culvert_depth=outlet.depth
                    flow_area=width*outlet.depth
                    perimeter=(width+2.0*outlet.depth)
                    case = 'Outlet is open channel flow'

                hyd_rad = flow_area/perimeter
                s = 'hydraulic radius at outlet = %f' %hyd_rad
                log_to_file(log_filename, s)

                # Outlet control velocity using tail water
                culvert_velocity = sqrt(delta_Et/((sum_loss/2*g)+(manning**2*length)/hyd_rad**1.33333)) 
                Q_outlet_tailwater = flow_area * culvert_velocity

                s = 'Q_outlet_tailwater = %.6f' %Q_outlet_tailwater
                log_to_file(log_filename, s)
                Q = min(Q, Q_outlet_tailwater)


        # Common code for circle and square geometries
        log_to_file(log_filename, 'Case: "%s"' %case)
        flow_rate_control=Q

        s = 'Flow Rate Control = %f' %flow_rate_control
        log_to_file(log_filename, s)

        inlet.rate = -flow_rate_control
        outlet.rate = flow_rate_control                
            
        culv_froude=sqrt(flow_rate_control**2*width/(g*flow_area**3))
        s = 'Froude in Culvert = %f' %culv_froude
        log_to_file(log_filename, s)

        # Determine momentum at the outlet 
        barrel_velocity = Q/(flow_area + velocity_protection/flow_area)


    else: #inlet.depth < 0.01:
        Q = barrel_velocity = outlet_culvert_depth = 0.0
        
    return Q, barrel_velocity, outlet_culvert_depth 


