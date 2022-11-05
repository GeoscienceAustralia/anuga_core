
import anuga
import math

class Boyd_box_operator(anuga.Structure_operator):
    """Culvert flow - transfer water from one rectangular box to another.
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
                 losses,
                 width,
                 height=None,
                 apron=None,
                 manning=0.013,
                 enquiry_gap=0.2,
                 use_momentum_jet=True,
                 use_velocity_head=True,
                 description=None,
                 label=None,
                 structure_type='boyd_box',
                 logging=False,
                 verbose=False):

        raise Exception('(Feb 2015) This code seems broken -- the structure operator call is incorrect -- retire?')

        anuga.Structure_operator.__init__(self,
                                          domain,
                                          end_point0, 
                                          end_point1,
                                          width,
                                          height,
                                          apron,
                                          manning,
                                          enquiry_gap,                                                       
                                          description,
                                          label,
                                          structure_type,
                                          logging,
                                          verbose)     
        
        
        if isinstance(losses, dict):
            self.sum_loss = sum(losses.values())
        elif isinstance(losses, list):
            self.sum_loss = sum(losses)
        else:
            self.sum_loss = losses
        
        self.use_momentum_jet = use_momentum_jet
        self.use_velocity_head = use_velocity_head
        
        self.culvert_length = self.get_culvert_length()
        self.culvert_width = self.get_culvert_width()
        self.culvert_height = self.get_culvert_height()

        self.max_velocity = 10.0

        self.inlets = self.get_inlets()


        # Stats
        
        self.discharge = 0.0
        self.velocity = 0.0
        
        self.case = 'N/A'

    
    def discharge_routine(self):

        local_debug ='false'
        
        if self.inflow.get_enquiry_height() > 0.01: #this value was 0.01:  We should call this Enquiry Depth
            if local_debug =='true':
                anuga.log.critical('Specific E & Deltat Tot E = %s, %s'
                             % (str(self.inflow.get_enquiry_specific_energy()),
                                str(self.delta_total_energy)))
                anuga.log.critical('culvert type = %s' % str(culvert_type))
            # Water has risen above inlet



            msg = 'Specific energy at inlet is negative'
            assert self.inflow.get_enquiry_specific_energy() >= 0.0, msg

            if self.use_velocity_head :
                self.driving_energy = self.inflow.get_enquiry_specific_energy()
            else:
                self.driving_energy = self.inflow.get_enquiry_height()

            height = self.culvert_height
            width = self.culvert_width
            flow_width = self.culvert_width
            
            # intially assume the culvert flow is controlled by the inlet
            # check unsubmerged and submerged condition and use Min Q
            # but ensure the correct flow area and wetted perimeter are used
            if self.driving_energy < height:
                # Inlet Unsubmerged
                self.case = 'Inlet unsubmerged Box Acts as Weir'
                Q_inlet = 0.544*anuga.g**0.5*width*self.driving_energy**1.50 # Flow based on Inlet Ctrl Inlet Unsubmerged
                V_inlet = Q_inlet/(width*self.driving_energy)
            elif (self.driving_energy > height) and (self.driving_energy < 1.93*height):
                # Inlet in Transition Zone by Boyd  New equation
                self.case = 'Inlet in Transition between Weir & Orifice'
                Q_inlet = 0.54286*anuga.g**0.5*width*height**0.5*self.driving_energy
                dcrit = (Q_inlet**2/anuga.g/width**2)**0.333333 # Based on Inlet Control
                if dcrit < height:
                    V_inlet = Q_inlet/(width*dcrit) # Full Box Velocity
                else:
                    V_inlet = Q_inlet/(width*height) # Full Box Velocity
                
            else:
                # Inlet Submerged
                self.case = 'Inlet Submerged Box Acts as Orifice'
                Q_inlet = 0.702*anuga.g**0.5*width*height**0.89*self.driving_energy**0.61  # Flow based on Inlet Ctrl Inlet Submerged
                V_inlet = Q_inlet/(width*height) # Full Box Velocity
            
            Q = Q_inlet
            dcrit = (Q**2/anuga.g/width**2)**0.333333 # Based on Inlet Control
            # 
            # May not need this .... check if same is done above  Might move this block Yet ???
            outlet_culvert_depth = dcrit
            
            if outlet_culvert_depth > height:
                outlet_culvert_depth = height  # Once again the pipe is flowing full not partfull
                flow_area = width*height  # Cross sectional area of flow in the culvert
                perimeter = 2*(width+height)
                self.case = 'Inlet CTRL Outlet unsubmerged PIPE PART FULL'
            else:
                flow_area = width * outlet_culvert_depth
                perimeter = width+2*outlet_culvert_depth
                self.case = 'INLET CTRL Culvert is open channel flow we will for now assume critical depth'
            # Initial Estimate of Flow for Outlet Control using energy slope 
            #( may need to include Culvert Bed Slope Comparison)
            hyd_rad = flow_area/perimeter
            
            
            # NEED TO DETERMINE INTERNAL BARREL VELOCITY  Could you Direct Step Method or Standard Step to compute Profile & Friction loss more accurately
            # For Now Assume Normal Depth in the Pipe if Part full
            #culvert_velocity = math.sqrt(self.delta_total_energy/((self.sum_loss/2/anuga.g)+(self.manning**2*self.culvert_length)/hyd_rad**1.33333))
            # Calculate Culvert velocity Based on the greater of the Pipe Slope, or the Slope of the Energy Line
            if pipe_slope > energy_line_slope:
                slope = pipe_slope
            else:
                slope = energy_line_slope
            
            # Here need to iterate Depth or use Table
            while Qpartfull < Q:
                for i in numpy.arange(0.01, self.get_culvert_height(), 0.01):
                    partfull_Area= partfull_depth*width
                    partfull_Perimeter= width+2*partfull_depth
                    partfull_Hyd_Rad = partfull_Area/Partfull_Perimeter
                    Vpartfull = Partfull_Hyd_Rad**(2/3)*slope**0.5/self.manning
                    Qpartfull = Vpartfull*partfull_Area
                    if partfull_depth < dcrit:
                        flow_type = 'Super-critical in Barrel'
                    elif partfull_depth > dcrit:
                        flow_type = 'Sub-critical in Barrel'
                    else: #partfull = dcrit
                        flow_type = 'Critical Depth in Barrel'
            
            #culvert_velocity = math.sqrt(self.delta_total_energy/((self.manning**2*self.culvert_length)/hyd_rad**1.33333)) # Based only on Friction Loss
            #Q_outlet_tailwater = flow_area * culvert_velocity
            
            
            if self.delta_total_energy < self.driving_energy:  # wE SHOULD DO THIS ANYWAY NOT ONLY FOR THIS CONDITION OTHER WISE MIGHT MISS oUTLET CONTROL SOMETIMES
                # Calculate flows for outlet control

                # Determine the depth at the outlet relative to the depth of flow in the Culvert
                if self.outflow.get_enquiry_height() > height:        # The Outlet is Submerged
                    outlet_culvert_depth=height
                    flow_area=width*height       # Cross sectional area of flow in the culvert
                    perimeter=2.0*(width+height)
                    self.case = 'Outlet submerged Culvert Flowing Full'
                else:   # Here really should use the Culvert Slope to calculate Actual Culvert Depth & Velocity
                    dcrit = (Q**2/anuga.g/width**2)**0.333333
                    outlet_culvert_depth=dcrit   # For purpose of calculation assume the outlet depth = Critical Depth
                    if outlet_culvert_depth > height:
                        outlet_culvert_depth=height
                        flow_area=width*height
                        perimeter=2.0*(width+height)
                        self.case = 'Outlet is Flowing Full'
                    else:
                        flow_area=width*outlet_culvert_depth
                        perimeter=(width+2.0*outlet_culvert_depth)
                        self.case = 'Outlet is open channel flow'

                hyd_rad = flow_area/perimeter



                # Rename this next one  V_outlet
                culvert_velocity = math.sqrt(self.delta_total_energy/((self.sum_loss/2/anuga.g)+(self.manning**2*self.culvert_length)/hyd_rad**1.33333))
                Q_outlet_tailwater = flow_area * culvert_velocity
                
                 # Final Outlet control velocity using tail water
                # Determine HEad Loss based on Inlet, Barrel and Outlet Velocities
                

                #FIXME SR: Is this code used?
                Inlet_Loss = Vinlet**2/(2*g)* Inlet_Loss_Coeff
                Barrel_Loss = self.culvert.length*Vpartfull**2/(partfull_Hyd_Rad**(4/3))
                Bend_Loss = Vpartfull**2/(2*g)* Bend_Loss_Coeff
                #Other_Loss ???
                Exit_Loss = culvert_velocity**2/(2*g)*Exit_Loss_Coeff
                Total_Loss = Inlet_Loss+Barrel_Loss+Bend_Loss+Exit_Loss

 
                Q = min(Q, Q_outlet_tailwater)
            else:
                pass
                #FIXME(Ole): What about inlet control?

            culv_froude=math.sqrt(Q**2*flow_width/(anuga.g*flow_area**3))
            if local_debug =='true':
                anuga.log.critical('FLOW AREA = %s' % str(flow_area))
                anuga.log.critical('PERIMETER = %s' % str(perimeter))
                anuga.log.critical('Q final = %s' % str(Q))
                anuga.log.critical('FROUDE = %s' % str(culv_froude))

            # Determine momentum at the outlet
            barrel_velocity = Q/(flow_area + anuga.velocity_protection/flow_area)

        # END CODE BLOCK for DEPTH  > Required depth for CULVERT Flow

        else: # self.inflow.get_enquiry_height() < 0.01:
            Q = barrel_velocity = outlet_culvert_depth = 0.0

        # Temporary flow limit
        if barrel_velocity > self.max_velocity:
            barrel_velocity = self.max_velocity
            Q = flow_area * barrel_velocity

        return Q, barrel_velocity, outlet_culvert_depth
        
