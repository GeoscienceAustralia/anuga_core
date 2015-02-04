import anuga
import math


#=====================================================================
# The class
#=====================================================================


class Weir_orifice_trapezoid_operator(anuga.Structure_operator):
    """Culvert flow - transfer water from one trapezoidal section to another.
    Sets up the geometry of problem
    
    This is the base class for culverts. Inherit from this class (and overwrite
    compute_discharge method for specific subclasses)
    
    Input: minimum arguments
         domain,
         losses (scalar, list or dictionary of losses),
         width  (= height if height not given)
    """


    def __init__(self,
                 domain,
                 losses,
                 width,
                 z1,
                 z2,
                 height=None,
                 end_points=None,
                 exchange_lines=None,
                 enquiry_points=None,
                 invert_elevations=None,
                 apron=0.1,
                 manning=0.013,
                 enquiry_gap=0.0,
                 use_momentum_jet=True,
                 use_velocity_head=True,
                 description=None,
                 label=None,
                 structure_type='weir_orifice_trapezoid',
                 logging=False,
                 verbose=False):
                     
        anuga.Structure_operator.__init__(self,
                                          domain,
                                          end_points=end_points,
                                          exchange_lines=exchange_lines,
                                          enquiry_points=enquiry_points,
                                          invert_elevations=invert_elevations,
                                          width=width,
                                          height=height,
                                          culvert_slope=culvert_slope,
                                          z1=z1,
                                          z2=z2,
                                          diameter=None,
                                          apron=apron,
                                          manning=manning,
                                          enquiry_gap=enquiry_gap,
                                          description=description,
                                          label=label,
                                          structure_type=structure_type,
                                          logging=logging,
                                          verbose=verbose)

        
        if isinstance(losses, dict):
            self.sum_loss = sum(losses.values())
        elif isinstance(losses, list):
            self.sum_loss = sum(losses)
        else:
            self.sum_loss = losses
        
        self.use_momentum_jet = use_momentum_jet
        # Preserves the old default behaviour of zeroing momentum
        self.zero_outflow_momentum = not self.use_momentum_jet

        self.use_velocity_head = use_velocity_head
        
        self.culvert_length = self.get_culvert_length()
        self.culvert_width = self.get_culvert_width()
        self.culvert_height = self.get_culvert_height()
        self.culvert_slope = self.culvert_slope()
        self.culvert_z1 = self.get_culvert_z1()
        self.culvert_z2 = self.get_culvert_z2()
                
        #FIXME SR: Why is this hard coded!
        self.max_velocity = 10.0

        self.inlets = self.get_inlets()


        # Stats
        
        self.discharge = 0.0
        self.velocity = 0.0
        
        self.case = 'N/A'



    def discharge_routine(self):
        """Procedure to determine the inflow and outflow inlets.
        Then use Weir_Orifice_trapezoid_function to do the actual calculation
        """

        local_debug = False

        # If the cuvert has been closed, then no water gets through
        if self.culvert_height <= 0.0:
            Q = 0.0
            barrel_velocity = 0.0
            outlet_culvert_depth = 0.0
            self.case = "Culvert blocked"
            self.inflow  = self.inlets[0]
            self.outflow = self.inlets[1]
            return Q, barrel_velocity, outlet_culvert_depth



        #  delta_total_energy will determine which inlet is inflow
        if self.use_velocity_head:
            self.delta_total_energy = \
                 self.inlets[0].get_enquiry_total_energy() - self.inlets[1].get_enquiry_total_energy()
        else:
            self.delta_total_energy = \
                 self.inlets[0].get_enquiry_stage() - self.inlets[1].get_enquiry_stage()

        self.inflow  = self.inlets[0]
        self.outflow = self.inlets[1]

        if self.delta_total_energy < 0:
            self.inflow  = self.inlets[1]
            self.outflow = self.inlets[0]
            self.delta_total_energy = -self.delta_total_energy


        # Only calculate flow if there is some water at the inflow inlet.
        if self.inflow.get_enquiry_depth() > 0.01: #this value was 0.01:

            if local_debug:
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
                self.driving_energy = self.inflow.get_enquiry_depth()

            verbose = False
            if verbose:
                print 50*'='
                print 'width ',self.culvert_width
                print 'depth ',self.culvert_height
                print 'left bank slope ',self.culvert_z1
                print 'right bank slope ',self.culvert_z2                
                print 'flow_width ',self.culvert_width
                print 'length ' ,self.culvert_length
                print 'culvert slope', self.culvert_slope
                print 'driving_energy ',self.driving_energy
                print 'delta_total_energy ',self.delta_total_energy
                print 'outlet_enquiry_depth ',self.outflow.get_enquiry_depth()
                print 'inflow_enquiry_depth ',self.inflow.get_enquiry_depth()
                print 'outlet_enquiry_speed ',self.outflow.get_enquiry_speed()
                print 'inflow_enquiry_speed ',self.inflow.get_enquiry_speed()
                print 'sum_loss ',self.sum_loss
                print 'manning ',self.manning

            Q, barrel_velocity, outlet_culvert_depth, flow_area, case = \
                              weir_orifice_trapezoid_function(width =self.culvert_width,
                                                depth               =self.culvert_height,
                                                z1                  =self.culvert_z1,
                                                z2                  =self.culvert_z2, 
                                                culvert_slope       =self.culvert_slope,                                                     
                                                flow_width          =self.culvert_width,
                                                length              =self.culvert_length,
                                                driving_energy      =self.driving_energy,
                                                delta_total_energy  =self.delta_total_energy,
                                                outlet_enquiry_depth=self.outflow.get_enquiry_depth(),
                                                sum_loss            =self.sum_loss,
                                                manning             =self.manning)
        else:
             Q = barrel_velocity = outlet_culvert_depth = 0.0
             case = 'Inlet dry'


        self.case = case


        # Temporary flow limit
        if barrel_velocity > self.max_velocity:
            barrel_velocity = self.max_velocity
            Q = flow_area * barrel_velocity

        return Q, barrel_velocity, outlet_culvert_depth


        

#=============================================================================
# define separately so that can be imported in parallel code.
#=============================================================================
def weir_orifice_trapezoid_function(width, 
                        depth,
                        z1,
                        z2, 
                        flow_width, 
                        length, 
                        culvert_slope,
                        driving_energy, 
                        delta_total_energy, 
                        outlet_enquiry_depth, 
                        sum_loss, 
                        manning):

    # intially assume the culvert flow is controlled by the inlet
    # check unsubmerged and submerged condition and use Min Q
    # but ensure the correct flow area and wetted perimeter are used

    local_debug = False
    
    Q_inlet_unsubmerged = 1.7*((2*width+depth*(z1+z2))/2)*driving_energy**1.50 # Flow based on Inlet Ctrl Inlet Unsubmerged (weir flow equation)
    Q_inlet_submerged = 0.8*anuga.g**0.5*(0.5*depth*(2*width+depth*(z1+z2)))*driving_energy**0.5  # Flow based on Inlet Ctrl Inlet Submerged (orifice equation)

    # FIXME(Ole): Are these functions really for inlet control?
    if Q_inlet_unsubmerged < Q_inlet_submerged:
        Q = Q_inlet_unsubmerged
        # Critical Depths Calculation
        dcrit=0.00001
        ic=0
        dyc=0.001
        while abs(dyc)>0.00001:
            Tc=width+(z1+z2)*dcrit
            Pc=width+((z1**2+1)**0.5+(z2**2+1)**0.5)*dcrit
            Ac=0.5*dcrit*(width+Tc)
            Rc=Ac/Pc
            fc=Ac**1.5*Tc**-0.5-Q/(9.81**0.5)
            ffc=Ac**1.5*-0.5*Tc**-1.5*(z1+z2)+Tc**-0.5*1.5*Ac**0.5*Tc
            dyc=-fc/ffc
            dcrit=dcrit+dyc
            ic=ic+1
        dcrit = dcrit
        if dcrit > depth:
            dcrit = depth
            flow_area = width*dcrit+0.5*(z1+z2)*dcrit**2
            perimeter= 2.0*width+(z1+z2)*dcrit + (dcrit**2+(z1*dcrit)**2)**0.5 + (dcrit**2+(z2*dcrit)**2)**0.5
        else: # dcrit < depth
            flow_area = width*dcrit+0.5*(z1+z2)*dcrit**2
            perimeter= width + (dcrit**2+(z1*dcrit)**2)**0.5 + (dcrit**2+(z2*dcrit)**2)**0.5
        outlet_culvert_depth = dcrit
        case = 'Inlet unsubmerged Box Acts as Weir'
    else: # Inlet Submerged but check internal culvert flow depth
        Q = Q_inlet_submerged
        # Critical Depths Calculation
        dcrit=0.00001
        ic=0
        dyc=0.001
        while abs(dyc)>0.00001:
            Tc=width+(z1+z2)*dcrit
            Pc=width+((z1**2+1)**0.5+(z2**2+1)**0.5)*dcrit
            Ac=0.5*dcrit*(width+Tc)
            Rc=Ac/Pc
            fc=Ac**1.5*Tc**-0.5-Q/(9.81**0.5)
            ffc=Ac**1.5*-0.5*Tc**-1.5*(z1+z2)+Tc**-0.5*1.5*Ac**0.5*Tc
            dyc=-fc/ffc
            dcrit=dcrit+dyc
            ic=ic+1
        dcrit = dcrit
        if dcrit > depth:
            dcrit = depth
            flow_area = width*dcrit+0.5*(z1+z2)*dcrit**2
            perimeter= 2.0*width+(z1+z2)*dcrit + (dcrit**2+(z1*dcrit)**2)**0.5 + (dcrit**2+(z2*dcrit)**2)**0.5
        else: # dcrit < depth
            flow_area = width*dcrit+0.5*(z1+z2)*dcrit**2
            perimeter= width + (dcrit**2+(z1*dcrit)**2)**0.5 + (dcrit**2+(z2*dcrit)**2)**0.5
        outlet_culvert_depth = dcrit
        case = 'Inlet submerged Box Acts as Orifice'
    # Critical Depths Calculation
    dcrit=0.00001
    ic=0
    dyc=0.001
    while abs(dyc)>0.00001:
        Tc=width+(z1+z2)*dcrit
        Pc=width+((z1**2+1)**0.5+(z2**2+1)**0.5)*dcrit
        Ac=0.5*dcrit*(width+Tc)
        Rc=Ac/Pc
        fc=Ac**1.5*Tc**-0.5-Q/(9.81**0.5)
        ffc=Ac**1.5*-0.5*Tc**-1.5*(z1+z2)+Tc**-0.5*1.5*Ac**0.5*Tc
        dyc=-fc/ffc
        dcrit=dcrit+dyc
        ic=ic+1
    dcrit = dcrit
    # May not need this .... check if same is done above
    outlet_culvert_depth = dcrit
    if outlet_culvert_depth > depth:
        outlet_culvert_depth = depth  # Once again the pipe is flowing full not partfull
        flow_area = width*depth+0.5*(z1+z2)*depth**2  # Cross sectional area of flow in the culvert
        perimeter = 2.0*width+(z1+z2)*depth + (depth**2+(z1*depth)**2)**0.5 + (depth**2+(z2*depth)**2)**0.5
        case = 'Inlet CTRL Outlet unsubmerged PIPE PART FULL'
    else:
        flow_area = width*outlet_culvert_depth+0.5*(z1+z2)*outlet_culvert_depth**2
        perimeter = 2.0*width+(z1+z2)*outlet_culvert_depth + (outlet_culvert_depth**2+(z1*outlet_culvert_depth)**2)**0.5 + (outlet_culvert_depth**2+(z2*outlet_culvert_depth)**2)**0.5
        case = 'INLET CTRL Culvert is open channel flow we will for now assume critical depth'
    # Initial Estimate of Flow for Outlet Control using energy slope
    #( may need to include Culvert Bed Slope Comparison)
    hyd_rad = flow_area/perimeter
    culvert_velocity = math.sqrt(delta_total_energy/((sum_loss/2/anuga.g) \
                                                          +(manning**2*length)/hyd_rad**1.33333))
    Q_outlet_tailwater = flow_area * culvert_velocity


    if delta_total_energy < driving_energy:
        # Calculate flows for outlet control

        # Determine the depth at the outlet relative to the depth of flow in the Culvert
        if outlet_enquiry_depth > depth:        # The Outlet is Submerged
            outlet_culvert_depth=depth
            flow_area=width*depth+0.5*(z1+z2)*depth**2  # Cross sectional area of flow in the culvert
            perimeter=width + (depth**2+(z1*depth)**2)**0.5 + (depth**2+(z2*depth)**2)**0.5
            case = 'Outlet submerged'
        else:   # Here we use the Culvert Slope to calculate Actual Culvert Depth & Velocity             
            # Iterate to determine Normal Depth...........
            
            Q = min(Q, Q_outlet_tailwater)
            dnorm=0.00001
            idn=1
            dyn=0.001
            
            while abs(dyn)>0.0001:
                Tn=width+(z1+z2)*dnorm
                An=0.5*dnorm*(width+Tn)
                Pn=width+2*dnorm*((z1**2+z2**2)**0.5)
                Rn=An/Pn
                fn=((culvert_slope**0.5*An*Rn**0.67)/manning) - Q
                ffn=((culvert_slope**0.5)/manning)*(((Rn**0.67)*Tn)+(Tn/Pn)-(2*dnorm*Rn/Pn))
                dnorm=dnorm-fn/ffn        
                dyn=-fn/ffn
                idn=idn+1
            outlet_culvert_depth = dnorm
            if outlet_culvert_depth > depth:
                outlet_culvert_depth=depth
                flow_area=width*depth+0.5*(z1+z2)*depth**2
                perimeter=width + (depth**2+(z1*depth)**2)**0.5 + (depth**2+(z2*depth)**2)**0.5
                case = 'Outlet is Flowing Full'
            else:
                flow_area=width*outlet_culvert_depth+0.5*(z1+z2)*outlet_culvert_depth**2
                perimeter=width + (outlet_culvert_depth**2+(z1*outlet_culvert_depth)**2)**0.5 + (outlet_culvert_depth**2+(z2*outlet_culvert_depth)**2)**0.5
                case = 'Outlet is open channel flow'

        hyd_rad = flow_area/perimeter

        # Final Outlet control velocity using tail water
        culvert_velocity = math.sqrt(delta_total_energy/((sum_loss/2/anuga.g)\
                                                          +(manning**2*length)/hyd_rad**1.33333))
        Q_outlet_tailwater = flow_area * culvert_velocity

        Q = min(Q, Q_outlet_tailwater)
    else:

        pass
        #FIXME(Ole): What about inlet control?

    if  flow_area <= 0.0 :
        culv_froude = 0.0
    else:
        culv_froude=math.sqrt(Q**2*flow_width/(anuga.g*flow_area**3))
        
    if local_debug:
        anuga.log.critical('FLOW AREA = %s' % str(flow_area))
        anuga.log.critical('PERIMETER = %s' % str(perimeter))
        anuga.log.critical('Q final = %s' % str(Q))
        anuga.log.critical('FROUDE = %s' % str(culv_froude))

    # Determine momentum at the outlet
    barrel_velocity = Q/(flow_area + anuga.velocity_protection/flow_area)

    # END CODE BLOCK for DEPTH  > Required depth for CULVERT Flow

    return Q, barrel_velocity, outlet_culvert_depth, flow_area, case
