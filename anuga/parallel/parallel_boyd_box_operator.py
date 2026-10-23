
import anuga

from anuga.structures.boyd_box_operator import boyd_box_function, total_energy, smooth_discharge

from .parallel_structure_operator import Parallel_Structure_operator

class Parallel_Boyd_box_operator(Parallel_Structure_operator):
    """Culvert flow - transfer water from one rectangular box to another.
    Sets up the geometry of problem

    This is the base class for culverts. Inherit from this class (and overwrite
    compute_discharge method for specific subclasses)

    Input: Two points, pipe_size (either diameter or width, height),
    mannings_rougness,
    """

    def __init__(self,
                 domain,
                 losses=0.0,
                 width=None,
                 blockage=0.0,
                 barrels=1.0,
                 z1=0.0,
                 z2=0.0,
                 height=None,
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
                 verbose=False,
                 master_proc = 0,
                 procs = None,
                 inlet_master_proc = None,
                 inlet_procs = None,
                 enquiry_proc = None,
                 max_velocity = 10.0):

        if inlet_master_proc is None:
            inlet_master_proc = [0, 0]
        if enquiry_proc is None:
            enquiry_proc = [0, 0]

        Parallel_Structure_operator.__init__(self,
                                          domain=domain,
                                          end_points=end_points,
                                          exchange_lines=exchange_lines,
                                          enquiry_points=enquiry_points,
                                          invert_elevations=invert_elevations,
                                          width=width,
                                          height=height,
                                          blockage=blockage,
                                          barrels=barrels,
                                          z1=0.0,
                                          z2=0.0,
                                          diameter= None,
                                          apron=apron,
                                          manning=manning,
                                          enquiry_gap=enquiry_gap,
                                          use_momentum_jet=use_momentum_jet,
                                          zero_outflow_momentum=(not use_momentum_jet),
                                          use_old_momentum_method=True,
                                          always_use_Q_wetdry_adjustment=True,
                                          force_constant_inlet_elevations=False,
                                          description=description,
                                          label=label,
                                          structure_type=structure_type,
                                          logging=logging,
                                          verbose=verbose,
                                          master_proc=master_proc,
                                          procs=procs,
                                          inlet_master_proc=inlet_master_proc,
                                          inlet_procs=inlet_procs,
                                          enquiry_proc=enquiry_proc)

        if isinstance(losses, dict):
            self.sum_loss = sum(losses.values())
        elif isinstance(losses, list):
            self.sum_loss = sum(losses)
        else:
            self.sum_loss = losses

        self.use_momentum_jet = use_momentum_jet
        self.zero_outflow_momentum = (not use_momentum_jet)
        self.use_old_momentum_method = True
        self.use_velocity_head = use_velocity_head

        self.culvert_length = self.get_culvert_length()
        self.culvert_width = self.get_culvert_width()
        self.culvert_height = self.get_culvert_height()
        self.culvert_blockage = self.get_culvert_blockage()
        self.culvert_barrels = self.get_culvert_barrels()

        self.max_velocity = max_velocity

        self.inlets = self.get_inlets()


        # Stats

        self.discharge = 0.0
        self.velocity = 0.0

        self.case = 'N/A'

        self.domain=domain

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


        '''
        print "ATTRIBUTES OF PARALLEL BOYD BOX::"
        for attr in dir(self):
            print "obj.%s = %s" % (attr, getattr(self, attr))
        '''


    def parallel_safe(self):
        """
        Set that operator is parallel safe
        """

        return True



    def discharge_routine(self):
        """Get info from inlets and then call sequential function."""

        local_debug = False
        forward_Euler_smooth = True

        # Gather total energy and stage from both enquiry points → master
        enq_data = self._gather_enquiry_stage_and_energy()

        # Master computes delta_total_energy, smooths it, and sets flow direction
        self.inflow_index = 0
        self.outflow_index = 1
        reverse = False
        ts = 1.0

        if self.myid == self.master_proc:
            enq_total_energy0, enq_stage0, enq_total_energy1, enq_stage1 = enq_data
            if self.use_velocity_head:
                self.delta_total_energy = enq_total_energy0 - enq_total_energy1
            else:
                self.delta_total_energy = enq_stage0 - enq_stage1

            self.smooth_delta_total_energy, ts = total_energy(
                self.smooth_delta_total_energy, self.delta_total_energy,
                self.domain.timestep, self.smoothing_timescale, forward_Euler_smooth)

            if self.smooth_delta_total_energy < 0:
                self.inflow_index = 1
                self.outflow_index = 0
                self.delta_total_energy = -self.smooth_delta_total_energy
                reverse = True
            else:
                self.delta_total_energy = self.smooth_delta_total_energy

        # Broadcast direction to non-master procs
        self._broadcast_flow_direction(reverse)

        # Gather inflow depth + specific energy and outflow depth → master
        flow_data = self._gather_inflow_outflow_depths()

        if self.myid != self.master_proc:
            return None, None, None

        inflow_enq_depth, inflow_enq_specific_energy, outflow_enq_depth = flow_data

        if inflow_enq_depth > 0.01:
            if local_debug:
                anuga.log.info('Specific E & Deltat Tot E = %s, %s'
                             % (str(inflow_enq_specific_energy), str(self.delta_total_energy)))
                anuga.log.info('culvert type = %s' % self.__class__.__name__)

            assert inflow_enq_specific_energy >= 0.0, 'Specific energy at inlet is negative'

            self.driving_energy = inflow_enq_specific_energy if self.use_velocity_head else inflow_enq_depth

            Q, barrel_velocity, outlet_culvert_depth, flow_area, case = \
                boyd_box_function(depth=self.culvert_height,
                                  width=self.culvert_width,
                                  flow_width=self.culvert_width,
                                  length=self.culvert_length,
                                  blockage=self.culvert_blockage,
                                  barrels=self.culvert_barrels,
                                  driving_energy=self.driving_energy,
                                  delta_total_energy=self.delta_total_energy,
                                  outlet_enquiry_depth=outflow_enq_depth,
                                  sum_loss=self.sum_loss,
                                  manning=self.manning)

            self.smooth_Q, Q, barrel_velocity = smooth_discharge(
                self.smooth_delta_total_energy, self.smooth_Q, Q,
                flow_area, ts, forward_Euler_smooth)

        else:
            Q = barrel_velocity = outlet_culvert_depth = 0.0
            case = 'Inlet dry'

        self.case = case

        if barrel_velocity > self.max_velocity:
            barrel_velocity = self.max_velocity
            Q = flow_area * barrel_velocity

        return Q, barrel_velocity, outlet_culvert_depth
