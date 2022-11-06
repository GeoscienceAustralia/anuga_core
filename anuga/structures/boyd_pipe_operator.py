
import anuga
import math
import numpy
from anuga.structures.boyd_box_operator import total_energy, smooth_discharge


#=====================================================================
# The class
#=====================================================================

class Boyd_pipe_operator(anuga.Structure_operator):
    """Culvert flow - transfer water from one location to another via a circular pipe culvert.
    Sets up the geometry of problem

    This is the base class for culverts. Inherit from this class (and overwrite
    compute_discharge method for specific subclasses)

    Input: Two points, pipe_size (diameter),
    mannings_rougness,
    """


    def __init__(self,
                 domain,
                 losses,
                 diameter=None,
                 barrels=1.0,
                 blockage=0.0,
                 z1=0.0,
                 z2=0.0,
                 end_points=None,
                 exchange_lines=None,
                 enquiry_points=None,
                 invert_elevations=None,
                 apron=0.1,
                 manning=0.013,
                 enquiry_gap=0.2,
                 smoothing_timescale=0.0,
                 use_momentum_jet=True,
                 use_velocity_head=True,
                 description=None,
                 label=None,
                 structure_type='boyd_pipe',
                 logging=False,
                 verbose=False):

        anuga.Structure_operator.__init__(self,
                                          domain,
                                          end_points,
                                          exchange_lines,
                                          enquiry_points,
                                          invert_elevations=invert_elevations,
                                          width=None,
                                          height=None,
                                          diameter=diameter,
                                          blockage=blockage,
                                          barrels=barrels,
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
        self.zero_outflow_momentum = not use_momentum_jet

        self.use_velocity_head = use_velocity_head

        self.culvert_length = self.get_culvert_length()
        self.culvert_diameter = self.get_culvert_diameter()
        self.culvert_blockage = self.get_culvert_blockage()
        self.culvert_barrels = self.get_culvert_barrels()

        #print self.culvert_diameter
        self.max_velocity = 10.0

        self.inlets = self.get_inlets()


        # Stats

        self.discharge = 0.0
        self.velocity = 0.0

        self.case = 'N/A'

        # May/June 2014 -- allow 'smoothing ' of driving_energy, delta total energy, and outflow_enq_depth
        self.smoothing_timescale=0.
        self.smooth_delta_total_energy=0.
        self.smooth_Q=0.
        # Set them based on a call to the discharge routine with smoothing_timescale=0.
        # [values of self.smooth_* are required in discharge_routine, hence dummy values above]
        Qvd=self.discharge_routine()
        self.smooth_delta_total_energy=1.0*self.delta_total_energy
        self.smooth_Q=Qvd[0]
        # Finally, set the smoothing timescale we actually want
        self.smoothing_timescale=smoothing_timescale




    def discharge_routine(self):
        """Procedure to determine the inflow and outflow inlets.
        Then use boyd_pipe_function to do the actual calculation
        """

        local_debug = False

        # If the cuvert has been closed, then no water gets through
        if self.culvert_diameter <= 0.0:
            Q = 0.0
            barrel_velocity = 0.0
            outlet_culvert_depth = 0.0
            self.case = "Culvert blocked"
            self.inflow  = self.inlets[0]
            self.outflow = self.inlets[1]
            return Q, barrel_velocity, outlet_culvert_depth



        #  delta_total_energy will determine which inlet is inflow
        # Update 02/07/2014 -- now using smooth_delta_total_energy
        if self.use_velocity_head:
            self.delta_total_energy = \
                 self.inlets[0].get_enquiry_total_energy() - self.inlets[1].get_enquiry_total_energy()
        else:
            self.delta_total_energy = \
                 self.inlets[0].get_enquiry_stage() - self.inlets[1].get_enquiry_stage()

        forward_Euler_smooth = True
        self.smooth_delta_total_energy, ts = total_energy(self.smooth_delta_total_energy,
                                                        self.delta_total_energy,
                                                        self.domain.timestep,
                                                        self.smoothing_timescale,
                                                        forward_Euler_smooth)

        if self.smooth_delta_total_energy >= 0.:
            self.inflow  = self.inlets[0]
            self.outflow = self.inlets[1]
            self.delta_total_energy=self.smooth_delta_total_energy
        else:
            self.inflow  = self.inlets[1]
            self.outflow = self.inlets[0]
            self.delta_total_energy = -self.smooth_delta_total_energy


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



            Q, barrel_velocity, outlet_culvert_depth, flow_area, case = \
                              boyd_pipe_function(depth               =self.inflow.get_enquiry_depth(),
                                                diameter            =self.culvert_diameter,
                                                blockage            =self.culvert_blockage,
                                                barrels             =self.culvert_barrels,
                                                length              =self.culvert_length,
                                                driving_energy      =self.driving_energy,
                                                delta_total_energy  =self.delta_total_energy,
                                                outlet_enquiry_depth=self.outflow.get_enquiry_depth(),
                                                sum_loss            =self.sum_loss,
                                                manning             =self.manning)

            self.smooth_Q, Q, barrel_velocity = smooth_discharge(self.smooth_delta_total_energy,
                                                                self.smooth_Q,
                                                                Q,
                                                                flow_area,
                                                                ts,
                                                                forward_Euler_smooth)

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
def boyd_pipe_function(depth,
                        diameter,
                        blockage,
                        barrels,
                        length,
                        driving_energy,
                        delta_total_energy,
                        outlet_enquiry_depth,
                        sum_loss,
                        manning):


    local_debug = False


    """
    For a circular pipe the Boyd method reviews 3 conditions
    1. Whether the Pipe Inlet is Unsubmerged (acting as weir flow into the inlet)
    2. Whether the Pipe Inlet is Fully Submerged (acting as an Orifice)
    3. Whether the energy loss in the pipe results in the Pipe being controlled by Channel Flow of the Pipe

    For these conditions we also would like to assess the pipe flow characteristics as it leaves the pipe
    """


    # Note this errors if blockage is set to 1.0 (ie 100% blockaage) and i have no idea how to fix it
    if blockage >= 1.0:
        Q = barrel_velocity = outlet_culvert_depth = 0.0
        flow_area = 0.00001
        case = '100 blocked culvert'
        return Q, barrel_velocity, outlet_culvert_depth, flow_area, case
    if blockage > 0.9:
        bf = 3.333-3.333*blockage
    else:
        bf = 1.0-0.4012316798*blockage-0.3768350138*(blockage**2)

    # Calculate flows for inlet control for circular pipe
    Q_inlet_unsubmerged = barrels * (0.421*anuga.g**0.5*((bf*diameter)**0.87)*driving_energy**1.63) # Inlet Ctrl Inlet Unsubmerged
    Q_inlet_submerged = barrels * (0.530*anuga.g**0.5*((bf*diameter)**1.87)*driving_energy**0.63)   # Inlet Ctrl Inlet Submerged
    # Note for to SUBMERGED TO OCCUR operator.inflow.get_average_specific_energy() should be > 1.2 x diameter.... Should Check !!!


    Q = min(Q_inlet_unsubmerged, Q_inlet_submerged)

    # THE LOWEST Value will Control Calcs From here
    # Calculate Critical Depth Based on the Adopted Flow as an Estimate
    dcrit1 = (bf*diameter)/1.26*(Q/anuga.g**0.5*((bf*diameter)**2.5))**(1/3.75)
    dcrit2 = (bf*diameter)/0.95*(Q/anuga.g**0.5*(bf*diameter)**2.5)**(1/1.95)
    # From Boyd Paper ESTIMATE of Dcrit has 2 criteria as
    if dcrit1/(bf*diameter)  > 0.85:
        outlet_culvert_depth = dcrit2
    else:
        outlet_culvert_depth = dcrit1
    #outlet_culvert_depth = min(outlet_culvert_depth, diameter)
    # Now determine Hydraulic Radius Parameters Area & Wetted Perimeter
    if outlet_culvert_depth >= (bf*diameter):
        outlet_culvert_depth = (bf*diameter)  # Once again the pipe is flowing full not partfull
        flow_area = barrels * (bf*diameter/2)**2 * math.pi  # Cross sectional area of flow in the culvert
        perimeter = barrels * bf * diameter * math.pi
        flow_width= barrels * bf * diameter
        case = 'Inlet CTRL Outlet submerged Circular PIPE FULL'
        if local_debug:
            anuga.log.critical('Inlet CTRL Outlet submerged Circular '
                         'PIPE FULL')
    else:
        #alpha = anuga.acos(1 - outlet_culvert_depth/diameter)    # Where did this Come From ????/
        alpha = anuga.acos(1-2*outlet_culvert_depth/(bf*diameter))*2
        #flow_area = diameter**2 * (alpha - sin(alpha)*cos(alpha))        # Pipe is Running Partly Full at the INLET   WHRE did this Come From ?????
        flow_area = barrels * (bf*diameter)**2/8*(alpha - math.sin(alpha))   # Equation from  GIECK 5th Ed. Pg. B3
        flow_width= barrels * bf*diameter*math.sin(alpha/2.0)
        perimeter = barrels * (alpha*bf*diameter/2.0)
        case = 'INLET CTRL Culvert is open channel flow we will for now assume critical depth'
        if local_debug:
            anuga.log.critical('INLET CTRL Culvert is open channel flow '
                         'we will for now assume critical depth')
            anuga.log.critical('Q Outlet Depth and ALPHA = %s, %s, %s'
                         % (str(Q), str(outlet_culvert_depth),
                            str(alpha)))
    if delta_total_energy < driving_energy:  #  OUTLET CONTROL !!!!
        # Calculate flows for outlet control

        # Determine the depth at the outlet relative to the depth of flow in the Culvert
        if outlet_enquiry_depth > bf*diameter:       # Outlet is submerged Assume the end of the Pipe is flowing FULL
            outlet_culvert_depth=bf*diameter
            flow_area = barrels * (bf*diameter/2)**2 * math.pi  # Cross sectional area of flow in the culvert
            perimeter = barrels * bf*diameter * math.pi
            flow_width= barrels * bf*diameter
            case = 'Outlet submerged'
            if local_debug:
                anuga.log.critical('Outlet submerged')
        else:   # Culvert running PART FULL for PART OF ITS LENGTH   Here really should use the Culvert Slope to calculate Actual Culvert Depth & Velocity
            # IF  operator.outflow.get_average_depth() < diameter
            dcrit1 = (bf*diameter)/1.26*(Q/anuga.g**0.5*((bf*diameter)**2.5))**(1/3.75)
            dcrit2 = (bf*diameter)/0.95*(Q/anuga.g**0.5*((bf*diameter)**2.5))**(1/1.95)
            if dcrit1/(bf*diameter) > 0.85:
                outlet_culvert_depth= dcrit2
            else:
                outlet_culvert_depth = dcrit1
            if outlet_culvert_depth > bf*diameter:
                outlet_culvert_depth = bf*diameter  # Once again the pipe is flowing full not partfull
                flow_area = barrels * (bf*diameter/2)**2 * math.pi  # Cross sectional area of flow in the culvert
                perimeter = barrels * bf*diameter * math.pi
                flow_width= barrels * bf*diameter
                case = 'Outlet unsubmerged PIPE FULL'
                if local_debug:
                    anuga.log.critical('Outlet unsubmerged PIPE FULL')
            else:
                alpha = anuga.acos(1-2*outlet_culvert_depth/(bf*diameter))*2
                flow_area = barrels * (bf*diameter)**2/8*(alpha - math.sin(alpha))   # Equation from  GIECK 5th Ed. Pg. B3
                flow_width= barrels * bf*diameter*math.sin(alpha/2.0)
                perimeter = barrels * alpha*bf*diameter/2.0
                case = 'Outlet is open channel flow we will for now assume critical depth'
                if local_debug:
                    anuga.log.critical('Q Outlet Depth and ALPHA = %s, %s, %s'
                                 % (str(Q), str(outlet_culvert_depth),
                                    str(alpha)))
                    anuga.log.critical('Outlet is open channel flow we '
                                 'will for now assume critical depth')
    if local_debug:
        anuga.log.critical('FLOW AREA = %s' % str(flow_area))
        anuga.log.critical('PERIMETER = %s' % str(perimeter))
        anuga.log.critical('Q Interim = %s' % str(Q))
    hyd_rad = flow_area/perimeter


    # Outlet control velocity using tail water
    if local_debug:
        anuga.log.critical('GOT IT ALL CALCULATING Velocity')
        anuga.log.critical('HydRad = %s' % str(hyd_rad))
    # Calculate Pipe Culvert Outlet Control Velocity.... May need initial Estimate First ??

    culvert_velocity = math.sqrt(delta_total_energy/((sum_loss/2/anuga.g)+\
                                                              (manning**2*length)/hyd_rad**1.33333))
    Q_outlet_tailwater = flow_area * culvert_velocity


    if local_debug:
        anuga.log.critical('VELOCITY = %s' % str(culvert_velocity))
        anuga.log.critical('Outlet Ctrl Q = %s' % str(Q_outlet_tailwater))

    Q = min(Q, Q_outlet_tailwater)
    if local_debug:
        anuga.log.critical('%s,%.3f,%.3f'
                     % ('dcrit 1 , dcit2 =',dcrit1,dcrit2))
        anuga.log.critical('%s,%.3f,%.3f,%.3f'
                     % ('Q and Velocity and Depth=', Q,
                        culvert_velocity, outlet_culvert_depth))

    culv_froude=math.sqrt(Q**2*flow_width*barrels/(anuga.g*flow_area**3))
    if local_debug:
        anuga.log.critical('FLOW AREA = %s' % str(flow_area))
        anuga.log.critical('PERIMETER = %s' % str(perimeter))
        anuga.log.critical('Q final = %s' % str(Q))
        anuga.log.critical('FROUDE = %s' % str(culv_froude))
        anuga.log.critical('Case = %s' % case)

    # Determine momentum at the outlet
    barrel_velocity = Q/(flow_area + anuga.velocity_protection/flow_area)

    # END CODE BLOCK for DEPTH  > Required depth for CULVERT Flow

    return Q, barrel_velocity, outlet_culvert_depth, flow_area, case
