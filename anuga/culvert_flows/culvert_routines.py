"""Collection of culvert routines for use with Culvert_flow in culvert_class

This module holds various routines to determine FLOW through CULVERTS and SIMPLE BRIDGES




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
                                   inlet_velocity,
                                   outlet_velocity,
                                   inlet_specific_energy, 
                                   delta_total_energy, 
                                   g,
                                   culvert_length=0.0,
                                   culvert_width=0.0,
                                   culvert_height=0.0,
                                   culvert_type='box',
                                   manning=0.0,
                                   sum_loss=0.0,
                                   max_velocity=10.0,
                                   log_filename=None):

    """Boyd's generalisation of the US department of transportation culvert methods
	
	WARNING THIS IS A SIMPLISTIC APPROACH and OUTLET VELOCITIES ARE LIMITED TO EITHER
	FULL PIPE OR CRITICAL DEPTH ONLY
	For Supercritical flow this is UNDERESTIMATING the Outlet Velocity
	The obtain the CORRECT velocity requires an iteration of Depth to Establish 
	the Normal Depth of flow in the pipe.
	
	It is proposed to provide this in a seperate routine called
    boyd_generalised_culvert_model_complex
	
    The Boyd Method is based on methods described by the following:
	1.
	US Dept. Transportation Federal Highway Administration (1965) 
	Hydraulic Chart for Selection of Highway Culverts.
	Hydraulic Engineering Circular No. 5 US Government Printing
	2.
	US Dept. Transportation Federal Highway Administration (1972) 
	Capacity charts for the Hydraulic design of highway culverts.
	Hydraulic Engineering Circular No. 10 US Government Printing
    These documents provide around 60 charts for various configurations of culverts and inlets.
	
	Note these documents have been superceded by:
	2005 Hydraulic Design of Highway Culverts, Hydraulic Design Series No. 5 (HDS-5), 
	Which combines culvert design information previously contained in Hydraulic Engineering Circulars 
	(HEC) No. 5, No. 10, and No. 13 with hydrologic, storage routing, and special culvert design information.
	HEC-5 provides 20 Charts
	HEC-10 Provides an additional 36 Charts
	HEC-13 Discusses the Design of improved more efficient inlets
	HDS-5 Provides 60 sets of Charts

	In 1985 Professor Michael Boyd Published "Head-Discharge Relations for Culverts", and in
	1987 published "Generalised Head Discharge Equations for Culverts".
	These papers reviewed the previous work by the US DOT and provided a simplistic approach for 3 configurations.
	
	It may be possible to extend the same approach for additional charts in the original work, but to date this has not been done.
	The additional charts cover a range of culvert shapes and inlet configurations
	
   
    """

    from anuga.utilities.system_tools import log_to_file
    from anuga.config import velocity_protection
    from anuga.utilities.numerical_tools import safe_acos as acos

    local_debug ='false'
    if inlet_depth > 0.1: #this value was 0.01:
        if local_debug =='true':
            log.critical('Specific E & Deltat Tot E = %s, %s'
                         % (str(inlet_specific_energy),
                            str(delta_total_energy)))
            log.critical('culvert type = %s' % str(culvert_type))
        # Water has risen above inlet
        
        if log_filename is not None:                            
            s = 'Specific energy  = %f m' % inlet_specific_energy
            log_to_file(log_filename, s)
        
        msg = 'Specific energy at inlet is negative'
        assert inlet_specific_energy >= 0.0, msg
                      

        if culvert_type =='circle':
            # Round culvert (use height as diameter)
            diameter = culvert_height   

            """
            For a circular pipe the Boyd method reviews 3 conditions
            1. Whether the Pipe Inlet is Unsubmerged (acting as weir flow into the inlet)
            2. Whether the Pipe Inlet is Fully Submerged (acting as an Orifice)
            3. Whether the energy loss in the pipe results in the Pipe being controlled by Channel Flow of the Pipe 
            
            For these conditions we also would like to assess the pipe flow characteristics as it leaves the pipe
            """
            
            # Calculate flows for inlet control            
            Q_inlet_unsubmerged = 0.421*g**0.5*diameter**0.87*inlet_specific_energy**1.63 # Inlet Ctrl Inlet Unsubmerged 
            Q_inlet_submerged = 0.530*g**0.5*diameter**1.87*inlet_specific_energy**0.63   # Inlet Ctrl Inlet Submerged
            # Note for to SUBMERGED TO OCCUR inlet_specific_energy should be > 1.2 x diameter.... Should Check !!!

            if log_filename is not None:                                        
                s = 'Q_inlet_unsubmerged = %.6f, Q_inlet_submerged = %.6f' % (Q_inlet_unsubmerged, Q_inlet_submerged)
                log_to_file(log_filename, s)
            Q = min(Q_inlet_unsubmerged, Q_inlet_submerged)

            # THE LOWEST Value will Control Calcs From here
            # Calculate Critical Depth Based on the Adopted Flow as an Estimate
            dcrit1 = diameter/1.26*(Q/g**0.5*diameter**2.5)**(1/3.75)
            dcrit2 = diameter/0.95*(Q/g**0.5*diameter**2.5)**(1/1.95)
            # From Boyd Paper ESTIMATE of Dcrit has 2 criteria as
            if dcrit1/diameter  > 0.85:
                outlet_culvert_depth = dcrit2
            else:
                outlet_culvert_depth = dcrit1
            #outlet_culvert_depth = min(outlet_culvert_depth, diameter)
            # Now determine Hydraulic Radius Parameters Area & Wetted Perimeter
            if outlet_culvert_depth >= diameter:
                outlet_culvert_depth = diameter  # Once again the pipe is flowing full not partfull
                flow_area = (diameter/2)**2 * pi  # Cross sectional area of flow in the culvert
                perimeter = diameter * pi
                flow_width= diameter
                case = 'Inlet CTRL Outlet submerged Circular PIPE FULL'
                if local_debug == 'true':
                    log.critical('Inlet CTRL Outlet submerged Circular '
                                 'PIPE FULL')
            else:
                #alpha = acos(1 - outlet_culvert_depth/diameter)    # Where did this Come From ????/
                alpha = acos(1-2*outlet_culvert_depth/diameter)*2
                #flow_area = diameter**2 * (alpha - sin(alpha)*cos(alpha))        # Pipe is Running Partly Full at the INLET   WHRE did this Come From ?????
                flow_area = diameter**2/8*(alpha - sin(alpha))   # Equation from  GIECK 5th Ed. Pg. B3
                flow_width= diameter*sin(alpha/2.0)
                perimeter = alpha*diameter/2.0
                case = 'INLET CTRL Culvert is open channel flow we will for now assume critical depth'
                if local_debug =='true':
                    log.critical('INLET CTRL Culvert is open channel flow '
                                 'we will for now assume critical depth')
                    log.critical('Q Outlet Depth and ALPHA = %s, %s, %s'
                                 % (str(Q), str(outlet_culvert_depth),
                                    str(alpha)))
            if delta_total_energy < inlet_specific_energy:  #  OUTLET CONTROL !!!!
                # Calculate flows for outlet control
                
                # Determine the depth at the outlet relative to the depth of flow in the Culvert
                if outlet_depth > diameter:       # Outlet is submerged Assume the end of the Pipe is flowing FULL
                    outlet_culvert_depth=diameter
                    flow_area = (diameter/2)**2 * pi  # Cross sectional area of flow in the culvert
                    perimeter = diameter * pi
                    flow_width= diameter
                    case = 'Outlet submerged'
                    if local_debug =='true':
                        log.critical('Outlet submerged')
                else:   # Culvert running PART FULL for PART OF ITS LENGTH   Here really should use the Culvert Slope to calculate Actual Culvert Depth & Velocity
                    # IF  outlet_depth < diameter
                    dcrit1 = diameter/1.26*(Q/g**0.5*diameter**2.5)**(1/3.75)
                    dcrit2 = diameter/0.95*(Q/g**0.5*diameter**2.5)**(1/1.95)
                    if dcrit1/diameter >0.85:
                        outlet_culvert_depth= dcrit2
                    else:
                        outlet_culvert_depth = dcrit1
                    if outlet_culvert_depth > diameter:
                        outlet_culvert_depth = diameter  # Once again the pipe is flowing full not partfull
                        flow_area = (diameter/2)**2 * pi  # Cross sectional area of flow in the culvert
                        perimeter = diameter * pi
                        flow_width= diameter
                        case = 'Outlet unsubmerged PIPE FULL'
                        if local_debug =='true':
                            log.critical('Outlet unsubmerged PIPE FULL')
                    else:
                        alpha = acos(1-2*outlet_culvert_depth/diameter)*2
                        flow_area = diameter**2/8*(alpha - sin(alpha))   # Equation from  GIECK 5th Ed. Pg. B3
                        flow_width= diameter*sin(alpha/2.0)
                        perimeter = alpha*diameter/2.0
                        case = 'Outlet is open channel flow we will for now assume critical depth'
                        if local_debug == 'true':
                            log.critical('Q Outlet Depth and ALPHA = %s, %s, %s'
                                         % (str(Q), str(outlet_culvert_depth),
                                            str(alpha)))
                            log.critical('Outlet is open channel flow we '
                                         'will for now assume critical depth')
            if local_debug == 'true':
                log.critical('FLOW AREA = %s' % str(flow_area))
                log.critical('PERIMETER = %s' % str(perimeter))
                log.critical('Q Interim = %s' % str(Q))
            hyd_rad = flow_area/perimeter

            if log_filename is not None:
                s = 'hydraulic radius at outlet = %f' %hyd_rad
                log_to_file(log_filename, s)

            # Outlet control velocity using tail water
            if local_debug =='true':
                log.critical('GOT IT ALL CALCULATING Velocity')
                log.critical('HydRad = %s' % str(hyd_rad))
            culvert_velocity = sqrt(delta_total_energy/((sum_loss/2/g)+(manning**2*culvert_length)/hyd_rad**1.33333)) 
            Q_outlet_tailwater = flow_area * culvert_velocity
            if local_debug =='true':
                log.critical('VELOCITY = %s' % str(culvert_velocity))
                log.critical('Outlet Ctrl Q = %s' % str(Q_outlet_tailwater))
            if log_filename is not None:                
                s = 'Q_outlet_tailwater = %.6f' %Q_outlet_tailwater
                log_to_file(log_filename, s)
            Q = min(Q, Q_outlet_tailwater)
            if local_debug =='true':
                log.critical('%s,%.3f,%.3f'
                             % ('dcrit 1 , dcit2 =',dcrit1,dcrit2))
                log.critical('%s,%.3f,%.3f,%.3f'
                             % ('Q and Velocity and Depth=', Q,
                                culvert_velocity, outlet_culvert_depth))

        else:
            pass
                #FIXME(Ole): What about inlet control?
        # ====  END OF CODE BLOCK FOR "IF" CIRCULAR PIPE 

        #else....
        if culvert_type == 'box':
            if local_debug == 'true':
                log.critical('BOX CULVERT')
            # Box culvert (rectangle or square)   ==========================================================

            # Calculate flows for inlet control
            height = culvert_height
            width = culvert_width
            flow_width=culvert_width            
            
            Q_inlet_unsubmerged = 0.540*g**0.5*width*inlet_specific_energy**1.50 # Flow based on Inlet Ctrl Inlet Unsubmerged
            Q_inlet_submerged = 0.702*g**0.5*width*height**0.89*inlet_specific_energy**0.61  # Flow based on Inlet Ctrl Inlet Submerged

            if log_filename is not None:                            
                s = 'Q_inlet_unsubmerged = %.6f, Q_inlet_submerged = %.6f' %(Q_inlet_unsubmerged, Q_inlet_submerged)
                log_to_file(log_filename, s)

            # FIXME(Ole): Are these functions really for inlet control?    
            if Q_inlet_unsubmerged < Q_inlet_submerged:
                Q = Q_inlet_unsubmerged
                dcrit = (Q**2/g/width**2)**0.333333
                if dcrit > height:
                    dcrit = height
                flow_area = width*dcrit
                outlet_culvert_depth = dcrit
                case = 'Inlet unsubmerged Box Acts as Weir'
            else:    
                Q = Q_inlet_submerged
                flow_area = width*height
                outlet_culvert_depth = height
                case = 'Inlet submerged Box Acts as Orifice'                    

            dcrit = (Q**2/g/width**2)**0.333333

            outlet_culvert_depth = dcrit
            if outlet_culvert_depth > height:
                outlet_culvert_depth = height  # Once again the pipe is flowing full not partfull
                flow_area = width*height  # Cross sectional area of flow in the culvert
                perimeter = 2*(width+height)
                case = 'Inlet CTRL Outlet unsubmerged PIPE PART FULL'
            else:
                flow_area = width * outlet_culvert_depth
                perimeter = width+2*outlet_culvert_depth
                case = 'INLET CTRL Culvert is open channel flow we will for now assume critical depth'
	
            if delta_total_energy < inlet_specific_energy:
                # Calculate flows for outlet control
                
                # Determine the depth at the outlet relative to the depth of flow in the Culvert
                if outlet_depth > height:        # The Outlet is Submerged
                    outlet_culvert_depth=height
                    flow_area=width*height       # Cross sectional area of flow in the culvert
                    perimeter=2.0*(width+height)
                    case = 'Outlet submerged'
                else:   # Here really should use the Culvert Slope to calculate Actual Culvert Depth & Velocity
                    dcrit = (Q**2/g/width**2)**0.333333
                    outlet_culvert_depth=dcrit   # For purpose of calculation assume the outlet depth = Critical Depth
                    if outlet_culvert_depth > height:
                        outlet_culvert_depth=height
                        flow_area=width*height
                        perimeter=2.0*(width+height)
                        case = 'Outlet is Flowing Full'
                    else:
                        flow_area=width*outlet_culvert_depth
                        perimeter=(width+2.0*outlet_culvert_depth)
                        case = 'Outlet is open channel flow'

                hyd_rad = flow_area/perimeter

                if log_filename is not None:                                
                    s = 'hydraulic radius at outlet = %f' % hyd_rad
                    log_to_file(log_filename, s)

                # Outlet control velocity using tail water
                culvert_velocity = sqrt(delta_total_energy/((sum_loss/2/g)+(manning**2*culvert_length)/hyd_rad**1.33333))
                Q_outlet_tailwater = flow_area * culvert_velocity

                if log_filename is not None:                            
                    s = 'Q_outlet_tailwater = %.6f' % Q_outlet_tailwater
                    log_to_file(log_filename, s)
                Q = min(Q, Q_outlet_tailwater)
            else:
                pass
                #FIXME(Ole): What about inlet control?
        # ====  END OF CODE BLOCK FOR "IF" BOX 

                
        # Common code for circle and box geometries ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if log_filename is not None:
            log_to_file(log_filename, 'Case: "%s"' % case)
            
        if log_filename is not None:                        
            s = 'Flow Rate Control = %f' % Q
            log_to_file(log_filename, s)

        culv_froude=sqrt(Q**2*flow_width/(g*flow_area**3))
        if local_debug =='true':
            log.critical('FLOW AREA = %s' % str(flow_area))
            log.critical('PERIMETER = %s' % str(perimeter))
            log.critical('Q final = %s' % str(Q))
            log.critical('FROUDE = %s' % str(culv_froude))
        if log_filename is not None:                            
            s = 'Froude in Culvert = %f' % culv_froude
            log_to_file(log_filename, s)

        # Determine momentum at the outlet 
        barrel_velocity = Q/(flow_area + velocity_protection/flow_area)

    # END CODE BLOCK for DEPTH  > Required depth for CULVERT Flow  

    else: # inlet_depth < 0.01:
        Q = barrel_velocity = outlet_culvert_depth = 0.0

    # Temporary flow limit
    if barrel_velocity > max_velocity:
        if log_filename is not None:                            
            s = 'Barrel velocity was reduced from = %f m/s to %f m/s' % (barrel_velocity, max_velocity)
            log_to_file(log_filename, s)
        
        barrel_velocity = max_velocity
        Q = flow_area * barrel_velocity
        
        

        
    return Q, barrel_velocity, outlet_culvert_depth 


