import anuga
import math
import numpy


#=====================================================================
# The class
#=====================================================================


class Internal_boundary_operator(anuga.Structure_operator):
    """
       The internal_boundary_function must accept 2 input arguments (hw, tw). It 
       returns Q:
       - hw will always be the stage (or energy) at the enquiry_point[0]
       - tw will always be the stage (or energy) at the enquiry_point[1]
       - If flow is from hw to tw, then Q should be positive, otherwise Q
         should be negative

       def internal_boundary_function(hw, tw):
           # Compute Q here from headwater hw and tailwater hw
           return(Q)

       smoothing_timescale>0. can be used to make Q vary more slowly

    """


    def __init__(self,
                 domain,
                 internal_boundary_function,
                 width=1.,
                 height=1.,
                 end_points=None,
                 exchange_lines=None,
                 enquiry_points=None,
                 invert_elevation=None,
                 apron=0.0,
                 enquiry_gap=0.0,
                 use_velocity_head=False,
                 zero_outflow_momentum=False,
                 force_constant_inlet_elevations=True,
                 smoothing_timescale=0.0,
                 description=None,
                 label=None,
                 structure_type='internal_boundary',
                 logging=False,
                 verbose=True):

        if verbose:
            print '########################################'
            print 'INTERNAL BOUNDARY OPERATOR'
            print 'THIS IS EXPERIMENTAL'
            print 'SUBJECT TO CHANGE WITHOUT NOTICE'
            print '########################################'

        # Since no barrel_velocity is computed we cannot use_momentum_jet
        use_momentum_jet = False

        anuga.Structure_operator.__init__(self,
                                          domain,
                                          end_points=end_points,
                                          exchange_lines=exchange_lines,
                                          enquiry_points=enquiry_points,
                                          invert_elevations=[invert_elevation, invert_elevation],
                                          width=width,
                                          height=height,
                                          diameter=None,
                                          apron=apron,
                                          manning=None,
                                          enquiry_gap=enquiry_gap,
                                          use_momentum_jet=use_momentum_jet,
                                          zero_outflow_momentum=zero_outflow_momentum,
                                          use_old_momentum_method=False,
                                          force_constant_inlet_elevations=force_constant_inlet_elevations,
                                          description=description,
                                          label=label,
                                          structure_type=structure_type,
                                          logging=logging,
                                          verbose=verbose)

        self.internal_boundary_function = internal_boundary_function
        self.use_momentum_jet = use_momentum_jet
        self.use_velocity_head = use_velocity_head
        self.zero_outflow_momentum = zero_outflow_momentum

        #FIXME SR: Why is this hard coded!
        self.max_velocity = 99999999999.0

        self.inlets = self.get_inlets()

        # Stats
        self.discharge = 0.0
        self.velocity = 0.0
        self.case = 'N/A'
        self.driving_energy = 0.0
        self.delta_total_energy = 0.0

        # Allow 'smoothing ' of  discharge
        self.smoothing_timescale = 0.
        self.smooth_Q = 0.
        self.smooth_delta_total_energy = 0.
        # Set them based on a call to the discharge routine with smoothing_timescale=0.
        # [values of self.smooth_* are required in discharge_routine, hence dummy values above]
        Qvd = self.discharge_routine()
        self.smooth_Q = Qvd[0]
        # Finally, set the smoothing timescale we actually want
        self.smoothing_timescale = smoothing_timescale
        self.smooth_delta_total_energy = self.delta_total_energy
    
    ###########################################################################
    def discharge_routine(self):
        """Procedure to determine the inflow and outflow inlets.
        Then use self.internal_boundary_function to do the actual calculation
        """

        local_debug = False

        # If the structure has been closed, then no water gets through
        if self.height <= 0.0:
            Q = 0.0
            barrel_velocity = 0.0
            outlet_culvert_depth = 0.0
            self.case = "Structure is blocked"
            self.inflow = self.inlets[0]
            self.outflow = self.inlets[1]
            return Q, barrel_velocity, outlet_culvert_depth

        # Compute energy head or stage at inlets 0 and 1
        if self.use_velocity_head:
            self.inlet0_energy = self.inlets[0].get_enquiry_total_energy()
            self.inlet1_energy = self.inlets[1].get_enquiry_total_energy()
        else:
            self.inlet0_energy = self.inlets[0].get_enquiry_stage()
            self.inlet1_energy = self.inlets[1].get_enquiry_stage()
        
        # Store these variables for anuga's structure output
        self.driving_energy = max(self.inlet0_energy, self.inlet1_energy)
        self.delta_total_energy = self.inlet0_energy - self.inlet1_energy

        # Other variables required by anuga's structure operator are not used
        barrel_velocity = numpy.nan
        outlet_culvert_depth = numpy.nan
        flow_area = numpy.nan
        case = ''

        # ts is used for smoothing discharge and delta_total_energy
        if self.domain.timestep > 0.:
            ts = self.domain.timestep/max(self.domain.timestep, self.smoothing_timescale, 1.0e-30)
        else:
            ts = 1.0

        # Compute a 'smoothed' delta_total_energy
        self.smooth_delta_total_energy = self.smooth_delta_total_energy +\
            ts*(self.delta_total_energy - self.smooth_delta_total_energy)

        if numpy.sign(self.smooth_delta_total_energy) != numpy.sign(self.delta_total_energy):
            self.smooth_delta_total_energy = 0.

        # Compute the 'tailwater' energy from the 'headwater' energy and the smooth_delta_total_energy
        # Note if ts = 1 (no smoothing), then the raw inlet energies are used
        if self.inlet0_energy >= self.inlet1_energy:
            inlet0_energy = 1.0*self.inlet0_energy
            inlet1_energy = inlet0_energy - self.smooth_delta_total_energy

            # Compute discharge
            Q = self.internal_boundary_function(inlet0_energy, inlet1_energy)
        else:
            inlet1_energy = 1.0*self.inlet1_energy
            inlet0_energy = inlet1_energy + self.smooth_delta_total_energy

            # Compute discharge
            Q = self.internal_boundary_function(inlet0_energy, inlet1_energy)

        # Use time-smoothed discharge
        self.smooth_Q = self.smooth_Q + ts*(Q - self.smooth_Q)

        # Define 'inflow' and 'outflow' for anuga's structure operator
        if self.smooth_Q >= 0.:
            self.inflow  = self.inlets[0]
            self.outflow = self.inlets[1]
        else:
            self.inflow  = self.inlets[1]
            self.outflow = self.inlets[0]

        if numpy.sign(self.smooth_Q) != numpy.sign(Q):
            # The flow direction of the 'instantaneous Q' based on the
            # 'smoothed delta_total_energy' is not the same as the
            # direction of smooth_Q. To prevent 'jumping around', let's
            # set Q to zero
            Q = 0.
        else:
            # Make Q positive (for anuga's structure operator)
            Q = min( abs(self.smooth_Q), abs(Q) )

        return Q, barrel_velocity, outlet_culvert_depth
