
import anuga
import math
import numpy
from anuga.structures.boyd_box_operator import total_energy, smooth_discharge


# =====================================================================
# Critical depth in a circular conduit
# =====================================================================

def circular_critical_depth(Q, diameter, g):
    """Critical depth of discharge Q in a circular conduit of given diameter.

    Solves the critical-flow condition Q**2/g = A**3/T exactly, where for a
    flow depth y subtending a central angle theta = 2*acos(1 - 2*y/diameter)

        A = diameter**2/8 * (theta - sin(theta))   (flow area)
        T = diameter * sin(theta/2)                (top width)

    A**3/T increases monotonically from 0 (dry) to infinity (full), so a
    fixed-iteration bisection on theta always converges to the unique root.
    The result is strictly less than diameter; a pipe only runs full when an
    outlet/inlet condition forces it to, not through this function.

    Q is the discharge through ONE barrel. The same algorithm is used by the
    C kernel (circular_critical_depth in gpu_culvert_operator.c) so that
    CPU and GPU modes agree.

    References: V. T. Chow, Open-Channel Hydraulics (1959), ch. 4
    (critical flow, Q**2/g = A**3/T); W. O. Straub, "A quick and easy way
    to calculate critical and conjugate depths in circular open channels",
    Civil Engineering (ASCE) 48(12), 1978 (approximation used as a check in
    the tests).

    Example (1 m3/s in a 1 m pipe, g = 9.8):

    >>> round(circular_critical_depth(1.0, 1.0, 9.8), 4)
    0.5732
    >>> circular_critical_depth(0.0, 1.0, 9.8)
    0.0
    """
    if Q <= 0.0 or diameter <= 0.0:
        return 0.0
    target = Q * Q / g
    lo = 0.0
    hi = 2.0 * math.pi
    for _ in range(60):
        theta = 0.5 * (lo + hi)
        area = diameter * diameter / 8.0 * (theta - math.sin(theta))
        top_width = diameter * math.sin(0.5 * theta)
        if area * area * area > target * top_width:
            hi = theta
        else:
            lo = theta
    theta = 0.5 * (lo + hi)
    return 0.5 * diameter * (1.0 - math.cos(0.5 * theta))


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
                 losses=0.0,
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

        """Create a circular pipe culvert using the Boyd (1987) flow algorithm.

        Parameters
        ----------
        domain : anuga.Domain
            Shallow-water domain to which the culvert is attached.
        losses : float or list of float or dict, optional
            Head-loss coefficients. A scalar is used directly; a list is summed;
            a dict maps loss names to coefficients and the values are summed.
            Typical values: 0.5 for a sharp-edged inlet, 1.0 for a re-entrant
            inlet. Default 0.0.
        diameter : float
            Internal diameter (m) of the circular pipe barrel.
        barrels : float, optional
            Number of parallel identical barrels. Total discharge is multiplied
            by this value. Default 1.0.
        blockage : float, optional
            Fractional blockage of the culvert cross-section. Must be in
            [0, 1]; 0.0 is fully open, 1.0 is fully blocked. Default 0.0.
        z1 : float, optional
            Batter slope (rise/run) of the embankment at the first culvert end,
            used when computing the invert elevation from the DEM. Default 0.0.
        z2 : float, optional
            Batter slope (rise/run) of the embankment at the second culvert end.
            Default 0.0.
        end_points : list of [float, float], optional
            ``[[x1, y1], [x2, y2]]`` — centre coordinates (m) of each culvert
            end face. Exactly one of ``end_points`` or ``exchange_lines`` must
            be provided.
        exchange_lines : list of two 2-point lines, optional
            ``[[[x1,y1],[x2,y2]], [[x1,y1],[x2,y2]]]`` — line segments defining
            the inlet and outlet faces. Alternative to ``end_points``.
        enquiry_points : list of [float, float] or None, optional
            ``[[x1, y1], [x2, y2]]`` — explicit locations for the upstream and
            downstream head-enquiry points. Computed automatically from
            ``end_points`` and ``apron`` when None. Default None.
        invert_elevations : list of float or None, optional
            ``[e1, e2]`` — invert (floor) elevations (m AHD) at each culvert
            end. Read from the mesh when None. Default None.
        apron : float, optional
            Width (m) of the approach apron at each culvert end, used to offset
            the enquiry point from the inlet/outlet face. Default 0.1.
        manning : float, optional
            Manning's roughness coefficient for the culvert barrel
            (dimensionless). Typical values: 0.012 for smooth concrete,
            0.024 for corrugated-steel pipe. Default 0.013.
        enquiry_gap : float, optional
            Additional gap (m) beyond the apron edge for the head-enquiry point.
            Default 0.2.
        smoothing_timescale : float, optional
            Timescale (s) for exponential smoothing of computed discharge to
            reduce numerical oscillations. Set to 0.0 to disable. Default 0.0.
        use_momentum_jet : bool, optional
            If True, momentum is injected into the outflow cell along the
            culvert axis (momentum-jet model). If False, outflow momentum is
            zeroed. Default True.
        use_velocity_head : bool, optional
            If True, the velocity head at the inlet is included in the total
            driving energy. Default True.
        description : str or None, optional
            Human-readable description of the culvert, stored in statistics
            output. Default None.
        label : str or None, optional
            Short label used as a filename prefix for the CSV log when
            ``logging=True``. Default None.
        structure_type : str, optional
            Internal type identifier written to output files. Default
            ``'boyd_pipe'``.
        logging : bool, optional
            If True, write per-timestep flow statistics to a CSV file named
            ``<label>_<structure_type>.csv``. Default False.
        verbose : bool, optional
            If True, print diagnostic messages during construction and
            discharge calculations. Default False.
        """

        anuga.Structure_operator.__init__(self,
                                          domain,
                                          end_points=end_points,
                                          exchange_lines=exchange_lines,
                                          enquiry_points=enquiry_points,
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
                anuga.log.info('Specific E & Deltat Tot E = %s, %s'
                             % (str(self.inflow.get_enquiry_specific_energy()),
                                str(self.delta_total_energy)))
                anuga.log.info('culvert type = %s' % self.__class__.__name__)
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
    # Critical depth of the adopted flow, per barrel. (This previously used
    # the empirical pair dcrit1/dcrit2 transcribed as Q/sqrt(g)*D**2.5,
    # which is not dimensionless; the intended group is Q/(sqrt(g)*D**2.5).)
    dcrit = circular_critical_depth(Q/barrels, bf*diameter, anuga.g)
    outlet_culvert_depth = dcrit
    #outlet_culvert_depth = min(outlet_culvert_depth, diameter)
    # Now determine Hydraulic Radius Parameters Area & Wetted Perimeter
    if outlet_culvert_depth >= (bf*diameter):
        outlet_culvert_depth = (bf*diameter)  # Once again the pipe is flowing full not partfull
        flow_area = barrels * (bf*diameter/2)**2 * math.pi  # Cross sectional area of flow in the culvert
        perimeter = barrels * bf * diameter * math.pi
        flow_width= barrels * bf * diameter
        case = 'Inlet CTRL Outlet submerged Circular PIPE FULL'
        if local_debug:
            anuga.log.info('Inlet CTRL Outlet submerged Circular '
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
            anuga.log.info('INLET CTRL Culvert is open channel flow '
                         'we will for now assume critical depth')
            anuga.log.info('Q Outlet Depth and ALPHA = %s, %s, %s'
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
                anuga.log.info('Outlet submerged')
        else:   # Culvert running PART FULL for PART OF ITS LENGTH   Here really should use the Culvert Slope to calculate Actual Culvert Depth & Velocity
            # IF  operator.outflow.get_average_depth() < diameter
            dcrit = circular_critical_depth(Q/barrels, bf*diameter, anuga.g)
            outlet_culvert_depth = dcrit
            if outlet_culvert_depth > bf*diameter:
                outlet_culvert_depth = bf*diameter  # Once again the pipe is flowing full not partfull
                flow_area = barrels * (bf*diameter/2)**2 * math.pi  # Cross sectional area of flow in the culvert
                perimeter = barrels * bf*diameter * math.pi
                flow_width= barrels * bf*diameter
                case = 'Outlet unsubmerged PIPE FULL'
                if local_debug:
                    anuga.log.info('Outlet unsubmerged PIPE FULL')
            else:
                alpha = anuga.acos(1-2*outlet_culvert_depth/(bf*diameter))*2
                flow_area = barrels * (bf*diameter)**2/8*(alpha - math.sin(alpha))   # Equation from  GIECK 5th Ed. Pg. B3
                flow_width= barrels * bf*diameter*math.sin(alpha/2.0)
                perimeter = barrels * alpha*bf*diameter/2.0
                case = 'Outlet is open channel flow we will for now assume critical depth'
                if local_debug:
                    anuga.log.info('Q Outlet Depth and ALPHA = %s, %s, %s'
                                 % (str(Q), str(outlet_culvert_depth),
                                    str(alpha)))
                    anuga.log.info('Outlet is open channel flow we '
                                 'will for now assume critical depth')
    if local_debug:
        anuga.log.info('FLOW AREA = %s' % str(flow_area))
        anuga.log.info('PERIMETER = %s' % str(perimeter))
        anuga.log.info('Q Interim = %s' % str(Q))
    hyd_rad = flow_area/perimeter


    # Outlet control velocity using tail water
    if local_debug:
        anuga.log.info('GOT IT ALL CALCULATING Velocity')
        anuga.log.info('HydRad = %s' % str(hyd_rad))
    # Calculate Pipe Culvert Outlet Control Velocity.... May need initial Estimate First ??

    culvert_velocity = math.sqrt(delta_total_energy/((sum_loss/2/anuga.g)+\
                                                              (manning**2*length)/hyd_rad**1.33333))
    Q_outlet_tailwater = flow_area * culvert_velocity


    if local_debug:
        anuga.log.info('VELOCITY = %s' % str(culvert_velocity))
        anuga.log.info('Outlet Ctrl Q = %s' % str(Q_outlet_tailwater))

    Q = min(Q, Q_outlet_tailwater)
    if local_debug:
        anuga.log.info('%s,%.3f' % ('dcrit =', dcrit))
        anuga.log.info('%s,%.3f,%.3f,%.3f'
                     % ('Q and Velocity and Depth=', Q,
                        culvert_velocity, outlet_culvert_depth))

    # flow_width and flow_area already include all barrels
    culv_froude = math.sqrt(Q**2*flow_width/(anuga.g*flow_area**3))
    if local_debug:
        anuga.log.info('FLOW AREA = %s' % str(flow_area))
        anuga.log.info('PERIMETER = %s' % str(perimeter))
        anuga.log.info('Q final = %s' % str(Q))
        anuga.log.info('FROUDE = %s' % str(culv_froude))
        anuga.log.info('Case = %s' % case)

    # Determine momentum at the outlet
    barrel_velocity = Q/(flow_area + anuga.velocity_protection/flow_area)

    # END CODE BLOCK for DEPTH  > Required depth for CULVERT Flow

    return Q, barrel_velocity, outlet_culvert_depth, flow_area, case
