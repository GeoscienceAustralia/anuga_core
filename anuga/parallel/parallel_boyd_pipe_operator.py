
import anuga
import math
import numpy

from anuga.structures.boyd_pipe_operator import boyd_pipe_function
from anuga.structures.boyd_box_operator import total_energy, smooth_discharge

from .parallel_inlet_operator import Parallel_Inlet_operator
from .parallel_structure_operator import Parallel_Structure_operator

class Parallel_Boyd_pipe_operator(Parallel_Structure_operator):
    """Culvert flow - transfer water from one rectangular box to another.
    Sets up the geometry of problem

    This is the base class for culverts. Inherit from this class (and overwrite
    compute_discharge method for specific subclasses)

    Input: Two points, pipe_size (either diameter or width, height),
    mannings_rougness,
    """

    def __init__(self,
                 domain,
                 losses,
                 diameter,
                 blockage=0.0,
                 barrels=1.0,
                 z1=0.0,
                 z2=0.0,
                 #height=None,
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
                 structure_type='boyd_pipe',
                 logging=False,
                 verbose=False,
                 master_proc = 0,
                 procs = None,
                 inlet_master_proc = [0,0],
                 inlet_procs = None,
                 enquiry_proc = [0,0]):

        Parallel_Structure_operator.__init__(self,
                                          domain=domain,
                                          end_points=end_points,
                                          exchange_lines=exchange_lines,
                                          enquiry_points=enquiry_points,
                                          invert_elevations=invert_elevations,
                                          width=None,
                                          height=None,
                                          diameter=diameter,
                                          blockage=blockage,
                                          barrels=barrels,
                                          z1=0.0,
                                          z2=0.0,
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
        self.culvert_diameter = self.get_culvert_diameter()
        self.culvert_blockage = self.get_culvert_blockage()
        self.culvert_barrels = self.get_culvert_barrels()

        self.max_velocity = 10.0

        self.inlets = self.get_inlets()


        # Stats

        self.discharge = 0.0
        self.velocity = 0.0

        self.case = 'N/A'


        #print(80*'=')
        #print("DON'T USE BOYD PIPES AS ALGORITHM NOT VERIFIED YET")


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
        print "ATTRIBUTES OF PARALLEL BOYD PIPE::"
        for attr in dir(self):
            print "obj.%s = %s" % (attr, getattr(self, attr))
        '''


    def parallel_safe(self):
        """
        Set that operator is parallel safe
        """

        return True



    def discharge_routine(self):
        """
        Get info from inlets and then call sequential function
        """

        from anuga.utilities import parallel_abstraction as pypar

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

        #Send attributes of both enquiry points to the master proc
        if self.myid == self.master_proc:

            if self.myid == self.enquiry_proc[0]:
                enq_total_energy0 = self.inlets[0].get_enquiry_total_energy()
                enq_stage0 = self.inlets[0].get_enquiry_stage()
            else:
                enq_total_energy0 = pypar.receive(self.enquiry_proc[0])
                enq_stage0 = pypar.receive(self.enquiry_proc[0])


            if self.myid == self.enquiry_proc[1]:
                enq_total_energy1 = self.inlets[1].get_enquiry_total_energy()
                enq_stage1 = self.inlets[1].get_enquiry_stage()
            else:
                enq_total_energy1 = pypar.receive(self.enquiry_proc[1])
                enq_stage1 = pypar.receive(self.enquiry_proc[1])

        else:
            if self.myid == self.enquiry_proc[0]:
                pypar.send(self.inlets[0].get_enquiry_total_energy(), self.master_proc)
                pypar.send(self.inlets[0].get_enquiry_stage(), self.master_proc)

            if self.myid == self.enquiry_proc[1]:
                pypar.send(self.inlets[1].get_enquiry_total_energy(), self.master_proc)
                pypar.send(self.inlets[1].get_enquiry_stage(), self.master_proc)


        # Determine the direction of the flow
        if self.myid == self.master_proc:
            if self.use_velocity_head:
                self.delta_total_energy = enq_total_energy0 - enq_total_energy1
            else:
                self.delta_total_energy = enq_stage0 - enq_stage1

        self.inflow_index = 0
        self.outflow_index = 1
        # master proc orders reversal if applicable
        if self.myid == self.master_proc:
            # May/June 2014 -- change the driving forces gradually, with forward euler timestepping
            #
            forward_Euler_smooth = True
            self.smooth_delta_total_energy, ts = total_energy(self.smooth_delta_total_energy,
                                                            self.delta_total_energy,
                                                            self.domain.timestep,
                                                            self.smoothing_timescale,
                                                            forward_Euler_smooth)

            # Reverse the inflow and outflow direction?
            if self.smooth_delta_total_energy < 0:
                self.inflow_index = 1
                self.outflow_index = 0

                #self.delta_total_energy = -self.delta_total_energy
                self.delta_total_energy = -self.smooth_delta_total_energy

                for i in self.procs:
                    if i == self.master_proc: continue
                    pypar.send(True, i)
            else:
                self.delta_total_energy = self.smooth_delta_total_energy
                for i in self.procs:
                    if i == self.master_proc: continue
                    pypar.send(False, i)

            #print "ZZZZ: Delta total energy = %f" %(self.delta_total_energy)
        else:
            reverse = pypar.receive(self.master_proc)

            if reverse:
                self.inflow_index = 1
                self.outflow_index = 0

        # Get attribute from inflow enquiry point
        if self.myid == self.master_proc:

            if self.myid == self.enquiry_proc[self.inflow_index]:
                    inflow_enq_depth = self.inlets[self.inflow_index].get_enquiry_depth()
                    inflow_enq_specific_energy = self.inlets[self.inflow_index].get_enquiry_specific_energy()
            else:
                    inflow_enq_depth = pypar.receive(self.enquiry_proc[self.inflow_index])
                    inflow_enq_specific_energy = pypar.receive(self.enquiry_proc[self.inflow_index])
        else:
            if self.myid == self.enquiry_proc[self.inflow_index]:
                pypar.send(self.inlets[self.inflow_index].get_enquiry_depth(), self.master_proc)
                pypar.send(self.inlets[self.inflow_index].get_enquiry_specific_energy(), self.master_proc)

        # Get attribute from outflow enquiry point
        if self.myid == self.master_proc:

            if self.myid == self.enquiry_proc[self.outflow_index]:
                outflow_enq_depth = self.inlets[self.outflow_index].get_enquiry_depth()
            else:
                outflow_enq_depth = pypar.receive(self.enquiry_proc[self.outflow_index])

            #print "ZZZZZ: outflow_enq_depth = %f" %(outflow_enq_depth)

        else:
            if self.myid == self.enquiry_proc[self.outflow_index]:
                pypar.send(self.inlets[self.outflow_index].get_enquiry_depth(), self.master_proc)




        # Master proc computes return values
        if self.myid == self.master_proc:

            #inflow_enq_specific_energy


            if inflow_enq_depth > 0.01: #this value was 0.01:
                if local_debug:
                    anuga.log.critical('Specific E & Deltat Tot E = %s, %s'
                                 % (str(inflow_enq_specific_energy),
                                    str(self.delta_total_energy)))

                    anuga.log.critical('culvert type = %s' % str(self.structure_type))

                # Water has risen above inlet


                msg = 'Specific energy at inlet is negative'
                assert inflow_enq_specific_energy >= 0.0, msg

                if self.use_velocity_head :
                    self.driving_energy = inflow_enq_specific_energy
                else:
                    self.driving_energy = inflow_enq_depth


                Q, barrel_velocity, outlet_culvert_depth, flow_area, case = \
                              boyd_pipe_function(depth              =inflow_enq_depth,
                                                diameter            =self.culvert_diameter,
                                                length              =self.culvert_length,
                                                blockage            =self.culvert_blockage,
                                                barrels             =self.culvert_barrels,
                                                driving_energy      =self.driving_energy,
                                                delta_total_energy  =self.delta_total_energy,
                                                outlet_enquiry_depth=outflow_enq_depth,
                                                sum_loss            =self.sum_loss,
                                                manning             =self.manning)

                self.smooth_Q, Q, barrel_velocity = smooth_discharge(self.smooth_delta_total_energy,
                                                                    self.smooth_Q,
                                                                    Q,
                                                                    flow_area,
                                                                    ts,
                                                                    forward_Euler_smooth)

            else: # self.inflow.get_enquiry_depth() < 0.01:
                Q = barrel_velocity = outlet_culvert_depth = 0.0
                case = 'Inlet dry'


            self.case = case

            # Temporary flow limit
            if barrel_velocity > self.max_velocity:
                barrel_velocity = self.max_velocity
                Q = flow_area * barrel_velocity



            return Q, barrel_velocity, outlet_culvert_depth
        else:
            return None, None, None
