
import anuga
import math
import numpy

from anuga.structures.weir_orifice_trapezoid_operator import weir_orifice_trapezoid_function 

from .parallel_inlet_operator import Parallel_Inlet_operator
from .parallel_structure_operator import Parallel_Structure_operator

class Parallel_Weir_orifice_trapezoid_operator(Parallel_Structure_operator):
    """Culvert flow - transfer water from one trapezoid section to another.
    Sets up the geometry of problem
    
    This is the base class for culverts. Inherit from this class (and overwrite
    compute_discharge method for specific subclasses)
    
    Input: Two points, pipe_size (either diameter or width, height),
    mannings_rougness,
    """

    def __init__(self,
                 domain,
                 losses,
                 width,
                 height=None,
                 z1=None,
                 z2=None,
                 blockage=0.0,
                 barrels=1.0,
                 end_points=None,
                 exchange_lines=None,
                 enquiry_points=None,
                 invert_elevations=None,
                 #culvert_slope=None,
                 apron=0.1,
                 manning=0.013,
                 enquiry_gap=0.0,
                 smoothing_timescale=0.0,
                 use_momentum_jet=True,
                 use_velocity_head=True,
                 description=None,
                 label=None,
                 structure_type='weir_orifice_trapezoid',
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
                                          width=width,
                                          height=height,
                                          blockage=blockage,
                                          barrels=barrels,
                                          z1=z1,
                                          z2=z2,
                                          #culvert_slope=culvert_slope,
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
        
        #self.culvert_slope = self.get_culvert_slope()
        
        self.culvert_z1 = self.get_culvert_z1()
        self.culvert_z2 = self.get_culvert_z2()
        
        self.max_velocity = 10.0

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
        print "ATTRIBUTES OF PARALLEL WEIR ORIFICE TRAPEZOID::"
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
            forward_Euler_smooth=True
            if(forward_Euler_smooth):
                # To avoid 'overshoot' we ensure ts<1.
                if(self.domain.timestep>0.):
                    ts=self.domain.timestep/max(self.domain.timestep, self.smoothing_timescale,1.0e-06)
                else:
                    # This case is included in the serial version, which ensures the unit tests pass
                    # even when domain.timestep=0.0. 
                    # Note though the discontinuous behaviour as domain.timestep-->0. from above
                    ts=1.0
                self.smooth_delta_total_energy=self.smooth_delta_total_energy+\
                                        ts*(self.delta_total_energy-self.smooth_delta_total_energy)
            else:
                # Use backward euler -- the 'sensible' ts limitation is different in this case
                # ts --> Inf is reasonable and corresponds to the 'nosmoothing' case
                ts=self.domain.timestep/max(self.smoothing_timescale, 1.0e-06)
                self.smooth_delta_total_energy = (self.smooth_delta_total_energy+ts*(self.delta_total_energy))/(1.+ts)

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

                    anuga.log.critical('culvert type = %s' % str(culvert_type))

                # Water has risen above inlet


                msg = 'Specific energy at inlet is negative'
                assert inflow_enq_specific_energy >= 0.0, msg

                if self.use_velocity_head :
                    self.driving_energy = inflow_enq_specific_energy
                else:
                    self.driving_energy = inflow_enq_depth


                Q, barrel_velocity, outlet_culvert_depth, flow_area, case = \
                              weir_orifice_trapezoid_function(depth =self.culvert_height,
                                                width               =self.culvert_width,
                                                z1                  =self.culvert_z1,
                                                z2                  =self.culvert_z2,                                                
                                                flow_width          =self.culvert_width,
                                                length              =self.culvert_length,
                                                blockage            =self.culvert_blockage,
                                                barrels             =self.culvert_barrels,
                                                #culvert_slope       =self.culvert_slope,
                                                driving_energy      =self.driving_energy,
                                                delta_total_energy  =self.delta_total_energy,
                                                outlet_enquiry_depth=outflow_enq_depth,
                                                sum_loss            =self.sum_loss,
                                                manning             =self.manning)

                ################################################
                # Smooth discharge. This can reduce oscillations
                # 
                # NOTE: The sign of smooth_Q assumes that
                #   self.inflow_index=0 and self.outflow_index=1
                #   , whereas the sign of Q is always positive
                Qsign=(self.outflow_index-self.inflow_index) # To adjust sign of Q
                if(forward_Euler_smooth):
                    self.smooth_Q = self.smooth_Q +ts*(Q*Qsign-self.smooth_Q)
                else: 
                    # Try implicit euler method
                    self.smooth_Q = (self.smooth_Q+ts*(Q*Qsign))/(1.+ts)
                
                if numpy.sign(self.smooth_Q)!=Qsign:
                    # The flow direction of the 'instantaneous Q' based on the
                    # 'smoothed delta_total_energy' is not the same as the
                    # direction of smooth_Q. To prevent 'jumping around', let's
                    # set Q to zero
                    Q=0.
                else:
                    Q = min(abs(self.smooth_Q), Q) #abs(self.smooth_Q)
                barrel_velocity=Q/flow_area
            # END CODE BLOCK for DEPTH  > Required depth for CULVERT Flow

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
        
        
