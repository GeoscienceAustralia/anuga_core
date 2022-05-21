from __future__ import print_function
from __future__ import division
from builtins import str
from past.utils import old_div
import anuga
import math
import numpy


#=====================================================================
# The class
#=====================================================================

class Boyd_box_operator(anuga.Structure_operator):
    """Culvert flow - transfer water from one rectangular box to another.
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
                 height=None,
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
                 enquiry_gap=0.0,
                 smoothing_timescale=0.0,
                 use_momentum_jet=True,
                 use_velocity_head=True,
                 description=None,
                 label=None,
                 structure_type='boyd_box',
                 logging=False,
                 verbose=False):

        """Create a box culvert using Boyd flow algorithm


        :param domain: Culvert applied to this domain
        :param losses: Losses 
        :param width: Width of culvert
        :param height: height of culvert
        :param barrels: Number of barrels
        :param blockage: Set between 0.0 - 1.0 Set to 1.0 to close off culvert
        :param z1: Elevation of end of Culvert
        :param z2: Elevation of other end of Culvert
        :param end_points: [[x1,y1], [x2,y2]] of centre of ends of culvert
        :param exchange_lines: [ [[x1,y1], [x2,y2]], [[x1,y1], [x2,y2]] ] list of two lines defining ends of culvert
        :param enquiry_points: [[x1,y1], [x2,y2]] location of enquiry points
        :param invert_elevations: [ e1, e2 ] invert elevations of culvert inlets
        :param apron: 
        :param manning:
        :param enquiry_gap:
        :param smoothing_timescale:
        :param use_momentum_jet:
        :param use_velocity_head:
        :param description:
        :param label:
        :param structure_type:
        :param logging:
        :param verbose:
        
        """

        anuga.Structure_operator.__init__(self,
                                          domain,
                                          end_points=end_points,
                                          exchange_lines=exchange_lines,
                                          enquiry_points=enquiry_points,
                                          invert_elevations=invert_elevations,
                                          width=width,
                                          height=height,
                                          blockage=blockage,
                                          barrels=barrels,
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
        self.culvert_blockage = self.get_culvert_blockage()
        self.culvert_barrels = self.get_culvert_barrels()

        #FIXME SR: Why is this hard coded!
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


        if verbose:
            print(self.get_culvert_slope())



    def discharge_routine(self):
        """Procedure to determine the inflow and outflow inlets.
        Then use boyd_box_function to do the actual calculation
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
                anuga.log.critical('culvert type = %s' % str(self.type))
            # Water has risen above inlet


            msg = 'Specific energy at inlet is negative'
            assert self.inflow.get_enquiry_specific_energy() >= 0.0, msg

            if self.use_velocity_head :
                self.driving_energy = self.inflow.get_enquiry_specific_energy()
            else:
                self.driving_energy = self.inflow.get_enquiry_depth()

            verbose = False
            if verbose:
                print(50*'=')
                print('width ',self.culvert_width)
                print('depth ',self.culvert_height)
                print('culvert blockage ',self.culvert_blockage)
                print('number of culverts ',self.culvert_barrels)
                print('flow_width ',self.culvert_width)
                print('length ' ,self.culvert_length)
                print('driving_energy ',self.driving_energy)
                print('delta_total_energy ',self.delta_total_energy)
                print('outlet_enquiry_depth ',self.outflow.get_enquiry_depth())
                print('inflow_enquiry_depth ',self.inflow.get_enquiry_depth())
                print('outlet_enquiry_speed ',self.outflow.get_enquiry_speed())
                print('inflow_enquiry_speed ',self.inflow.get_enquiry_speed())
                print('sum_loss ',self.sum_loss)
                print('manning ',self.manning)

            Q, barrel_velocity, outlet_culvert_depth, flow_area, case = \
                              boyd_box_function(width               =self.culvert_width,
                                                depth               =self.culvert_height,
                                                blockage            =self.culvert_blockage,
                                                barrels             =self.culvert_barrels,
                                                flow_width          =self.culvert_width,
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
def boyd_box_function(width,
                        depth,
                        blockage,
                        barrels,
                        flow_width,
                        length,
                        driving_energy,
                        delta_total_energy,
                        outlet_enquiry_depth,
                        sum_loss,
                        manning):

    # intially assume the culvert flow is controlled by the inlet
    # check unsubmerged and submerged condition and use Min Q
    # but ensure the correct flow area and wetted perimeter are used

    local_debug = False

    #print(outlet_enquiry_depth)

    bf = 1 - blockage

    if blockage >= 1.0:
        Q = barrel_velocity = outlet_culvert_depth = 0.0
        flow_area = 0.00001
        case = '100 blocked culvert'
        return Q, barrel_velocity, outlet_culvert_depth, flow_area, case
    else:
        Q_inlet_unsubmerged = 0.544*anuga.g**0.5*bf*width*barrels*driving_energy**1.50 # Flow based on Inlet Ctrl Inlet Unsubmerged
        Q_inlet_submerged = 0.702*anuga.g**0.5*bf*width*barrels*depth**0.89*driving_energy**0.61  # Flow based on Inlet Ctrl Inlet Submerged

    #print 'blockage ', blockage
    # FIXME(Ole): Are these functions really for inlet control?
    if Q_inlet_unsubmerged < Q_inlet_submerged:
        Q = Q_inlet_unsubmerged
        dcrit = (old_div(old_div(Q**2,anuga.g),(bf*width*barrels)**2))**0.333333
        if dcrit > depth:
            dcrit = depth
            flow_area = bf*width*dcrit*barrels
            perimeter = 2.0*(bf*width*barrels + dcrit)
        else: # dcrit < depth
            flow_area = bf*width*barrels*dcrit
            perimeter = 2.0*dcrit + bf*width*barrels
        outlet_culvert_depth = dcrit
        case = 'Inlet unsubmerged Box Acts as Weir'
    else: # Inlet Submerged but check internal culvert flow depth
        Q = Q_inlet_submerged
        dcrit = (old_div(old_div(Q**2,anuga.g),(bf*width*barrels)**2))**0.333333
        if dcrit > depth:
            dcrit = depth
            flow_area = bf*width*barrels*dcrit
            perimeter = 2.0*(bf*width*barrels + dcrit)
        else: # dcrit < depth
            flow_area = bf*width*barrels*dcrit
            perimeter = 2.0*dcrit + bf*width*barrels
        outlet_culvert_depth = dcrit
        case = 'Inlet submerged Box Acts as Orifice'

    dcrit = (old_div(old_div(Q**2,anuga.g),(bf*width*barrels)**2))**0.333333

    # May not need this .... check if same is done above
    outlet_culvert_depth = dcrit

    if outlet_culvert_depth > depth:
        outlet_culvert_depth = depth  # Once again the pipe is flowing full not partfull
        flow_area = bf*width*barrels*depth  # Cross sectional area of flow in the culvert
        perimeter = 2*(bf*width*barrels + depth)
        case = 'Inlet CTRL Outlet unsubmerged PIPE PART FULL'
    else:
        flow_area = bf*width*barrels*outlet_culvert_depth
        perimeter = bf*width*barrels + 2*outlet_culvert_depth
        case = 'INLET CTRL Culvert is open channel flow we will for now assume critical depth'
    # Initial Estimate of Flow for Outlet Control using energy slope
    #( may need to include Culvert Bed Slope Comparison)

    hyd_rad = old_div(flow_area,perimeter)

    culvert_velocity = math.sqrt(old_div(delta_total_energy,((old_div(old_div(sum_loss,2),anuga.g)) \
                                                          +old_div((manning**2*length),hyd_rad**1.33333))))
    Q_outlet_tailwater = flow_area * culvert_velocity


    if delta_total_energy < driving_energy:
        # Calculate flows for outlet control

        # Determine the depth at the outlet relative to the depth of flow in the Culvert
        if outlet_enquiry_depth > depth:        # The Outlet is Submerged
            outlet_culvert_depth = depth
            flow_area = bf*width*barrels*depth   # Cross sectional area of flow in the culvert
            perimeter = 2.0*(bf*width*barrels + depth)
            case = 'Outlet submerged'
        else:   # Here really should use the Culvert Slope to calculate Actual Culvert Depth & Velocity
            dcrit = (old_div(old_div(Q**2,anuga.g),(bf*width*barrels)**2))**0.333333
            outlet_culvert_depth = dcrit   # For purpose of calculation assume the outlet depth = Critical Depth
            if outlet_culvert_depth > depth:
                outlet_culvert_depth = depth
                flow_area = bf*width*barrels*depth
                perimeter = 2.0*(bf*width*barrels + depth)
                case = 'Outlet is Flowing Full'
            else:
                flow_area = bf*width*barrels*outlet_culvert_depth
                perimeter = bf*width*barrels + 2.0*outlet_culvert_depth
                case = 'Outlet is open channel flow'

        hyd_rad = old_div(flow_area,perimeter)

        # Final Outlet control velocity using tail water
        culvert_velocity = math.sqrt(old_div(delta_total_energy,((old_div(old_div(sum_loss,2),anuga.g))\
                                                          +old_div((manning**2*length),hyd_rad**1.33333))))
        Q_outlet_tailwater = flow_area * culvert_velocity

        Q = min(Q, Q_outlet_tailwater)
    else:

        pass
        #FIXME(Ole): What about inlet control?

    if  flow_area <= 0.0 :
        culv_froude = 0.0
    else:
        culv_froude=math.sqrt(old_div(Q**2*flow_width*barrels,(anuga.g*flow_area**3)))

    if local_debug:
        anuga.log.critical('FLOW AREA = %s' % str(flow_area))
        anuga.log.critical('PERIMETER = %s' % str(perimeter))
        anuga.log.critical('Q final = %s' % str(Q))
        anuga.log.critical('FROUDE = %s' % str(culv_froude))
        anuga.log.critical('Case = %s' % case)

    # Determine momentum at the outlet
    barrel_velocity = old_div(Q,(flow_area + old_div(anuga.velocity_protection,flow_area)))

    # END CODE BLOCK for DEPTH  > Required depth for CULVERT Flow

    return Q, barrel_velocity, outlet_culvert_depth, flow_area, case

def total_energy(smooth_delta_total_energy,
                delta_total_energy,
                timestep,
                smoothing_timescale,
                forward_Euler_smooth=True):

    if(forward_Euler_smooth):
        # To avoid 'overshoot' we ensure ts<1.
        if(timestep>0.):
            ts=old_div(timestep,max(timestep, smoothing_timescale,1.0e-06))
        else:
            # Without this the unit tests with no smoothing fail [since they have domain.timestep=0.]
            # Note though the discontinuous behaviour as domain.timestep-->0. from above
            ts=1.0
        smooth_delta_total_energy=smooth_delta_total_energy+\
                                ts*(delta_total_energy-smooth_delta_total_energy)
    else:
        # Use backward euler -- the 'sensible' ts limitation is different in this case
        # ts --> Inf is reasonable and corresponds to the 'nosmoothing' case
        ts=old_div(timestep,max(smoothing_timescale, 1.0e-06))
        smooth_delta_total_energy = old_div((smooth_delta_total_energy+ts*(delta_total_energy)),(1.+ts))
    return smooth_delta_total_energy, ts

def smooth_discharge(smooth_delta_total_energy,
                    smooth_Q,
                    Q,
                    flow_area,
                    timestep,
                    forward_Euler_smooth=True):
    #
    # Update 02/07/2014 -- using time-smoothed discharge
    ts = timestep
    Qsign=numpy.sign(smooth_delta_total_energy) #(self.outflow_index-self.inflow_index) # To adjust sign of Q
    if(forward_Euler_smooth):
        smooth_Q = smooth_Q +ts*(Q*Qsign-smooth_Q)
    else:
        # Try implicit euler method
        smooth_Q = old_div((smooth_Q+ts*(Q*Qsign)),(1.+ts))

    if numpy.sign(smooth_Q)!=Qsign:
        # The flow direction of the 'instantaneous Q' based on the
        # 'smoothed delta_total_energy' is not the same as the
        # direction of smooth_Q. To prevent 'jumping around', let's
        # set Q to zero
        Q=0.
    else:
        Q = min(abs(smooth_Q), Q) #abs(self.smooth_Q)
    if flow_area == 0:
        barrel_velocity = 0.0
    else:
        barrel_velocity=old_div(Q,flow_area)
    return smooth_Q, Q, barrel_velocity
