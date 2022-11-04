
import anuga
import math
import numpy
from numpy.linalg import solve
import scipy
import scipy.optimize as sco

#from anuga.structures.boyd_box_operator import boyd_box_function 

from .parallel_inlet_operator import Parallel_Inlet_operator
from .parallel_structure_operator import Parallel_Structure_operator

class Parallel_Internal_boundary_operator(Parallel_Structure_operator):
    """
        Parallel variant of anuga.structures.Internal_boundary_operator
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
                 compute_discharge_implicitly=True,
                 description=None,
                 label=None,
                 structure_type='internal_boundary',
                 logging=False,
                 verbose=True,
                 master_proc = 0,
                 procs = None,
                 inlet_master_proc = [0,0],
                 inlet_procs = None,
                 enquiry_proc = [0,0]):

        if verbose:
            print('########################################')
            print('PARALLEL INTERNAL BOUNDARY OPERATOR')
            print('THIS IS EXPERIMENTAL')
            print('SUBJECT TO CHANGE WITHOUT NOTICE')
            print('########################################')

        # Since no barrel_velocity is computed we cannot use_momentum_jet
        use_momentum_jet = False
                     
        Parallel_Structure_operator.__init__(self,
                                          domain=domain,
                                          end_points=end_points,
                                          exchange_lines=exchange_lines,
                                          enquiry_points=enquiry_points,
                                          invert_elevations=[invert_elevation,invert_elevation],
                                          width=width,
                                          height=height,
                                          z1=0.0,
                                          z2=0.0,
                                          diameter= None,
                                          apron=apron,
                                          manning=None,
                                          blockage=None,#
                                          barrels=None,#
                                          enquiry_gap=enquiry_gap,
                                          use_momentum_jet=use_momentum_jet,
                                          zero_outflow_momentum=zero_outflow_momentum,
                                          use_old_momentum_method=False,
                                          always_use_Q_wetdry_adjustment=False,
                                          force_constant_inlet_elevations=force_constant_inlet_elevations,
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
       
 
        self.internal_boundary_function = internal_boundary_function
        self.use_momentum_jet = use_momentum_jet
        self.use_velocity_head = use_velocity_head
        self.zero_outflow_momentum = zero_outflow_momentum

        self.compute_discharge_implicitly = compute_discharge_implicitly
        
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

    #def __call__(self):
    #    """
    #        Update for n sub-timesteps 
    #    """

    #    number_of_substeps = 20
    #    original_timestep = self.domain.get_timestep()
    #    self.domain.timestep = original_timestep/(1.0*number_of_substeps)
    #    
    #    for i in range(number_of_substeps): 
    #        anuga.parallel.parallel_structure_operator.Parallel_Structure_operator.__call__(self)

    #    self.domain.timestep = original_timestep

    def parallel_safe(self):

        return True


    def discharge_routine(self):
        """Both implicit and explicit methods available
        The former seems more stable and more accurate (in at least some
        cases). The latter will have less communication in parallel, and
        for some simple internal_boundary_functions there is no benefit to
        the implicit approach
            
        """

        if self.compute_discharge_implicitly:
            Q, barrel_velocity, outlet_culvert_depth = self.discharge_routine_implicit()
        else:
            Q, barrel_velocity, outlet_culvert_depth = self.discharge_routine_explicit()

        return Q, barrel_velocity, outlet_culvert_depth

    def discharge_routine_explicit(self):

        from anuga.utilities import parallel_abstraction as pypar

        local_debug = False
        
        # If the structure has been closed, then no water gets through
        if self.height <= 0.0:
            if self.myid == self.master_proc:
                Q = 0.0
                barrel_velocity = 0.0
                outlet_culvert_depth = 0.0
                self.case = "Structure is blocked"
                self.inflow = self.inlets[0]
                self.outflow = self.inlets[1]
                return Q, barrel_velocity, outlet_culvert_depth
            else:
                return None, None, None

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
            # Variables required by anuga's structure operator which are not
            # used
            barrel_velocity = numpy.nan
            outlet_culvert_depth = numpy.nan
            flow_area = numpy.nan
            case = ''

            # 'Timescale' for smoothed discharge and energy
            ts = self.domain.timestep/max(self.domain.timestep, self.smoothing_timescale, 1.0e-30)

            # Energy or stage as head
            if self.use_velocity_head:
                E0 = enq_total_energy0
                E1 = enq_total_energy1
            else:
                E0 = enq_stage0
                E1 = enq_stage1

            self.delta_total_energy = E0 - E1
            self.driving_energy = max(E0, E1)

            # Compute 'smoothed' versions of key variables
            self.smooth_delta_total_energy += ts*(self.delta_total_energy - self.smooth_delta_total_energy)

            if numpy.sign(self.smooth_delta_total_energy) != numpy.sign(self.delta_total_energy):
                self.smooth_delta_total_energy = 0.

            # Compute the 'tailwater' energy from the 'headwater' energy
            # and the smooth_delta_total_energy. Note if ts = 1 (no
            # smoothing), then the raw inlet energies are used
            if E0 >= E1:
                inlet0_energy = 1.0*E0
                inlet1_energy = inlet0_energy - self.smooth_delta_total_energy

            else:
                inlet1_energy = 1.0*E1
                inlet0_energy = inlet1_energy + self.smooth_delta_total_energy

            # Compute discharge
            Q = self.internal_boundary_function(inlet0_energy, inlet1_energy)
            self.smooth_Q = self.smooth_Q + ts*(Q - self.smooth_Q)

            if numpy.sign(self.smooth_Q) != numpy.sign(Q):
                # The flow direction of the 'instantaneous Q' based on the
                # 'smoothed delta_total_energy' is not the same as the
                # direction of smooth_Q. To prevent 'jumping around', let's
                # set Q to zero
                Q = 0.
            else:
                # Make Q positive (for anuga's structure operator)
                Q = min( abs(self.smooth_Q), abs(Q) )
        else:
            self.delta_total_energy=numpy.nan
            self.driving_energy=numpy.nan


        self.inflow_index = 0
        self.outflow_index = 1
        # master proc orders reversal if applicable
        if self.myid == self.master_proc:

            # Reverse the inflow and outflow direction?
            if self.smooth_Q < 0.:
                self.inflow_index = 1
                self.outflow_index = 0

                for i in self.procs:
                    if i == self.master_proc: continue
                    pypar.send(True, i)
            else:
                for i in self.procs:
                    if i == self.master_proc: continue
                    pypar.send(False, i)

        else:
            reverse = pypar.receive(self.master_proc)

            if reverse:
                self.inflow_index = 1
                self.outflow_index = 0

        # Master proc computes return values
        if self.myid == self.master_proc:
            return Q, barrel_velocity, outlet_culvert_depth
        else:
            return None, None, None
        
        
    def discharge_routine_implicit(self):
        """
            Uses semi-implicit discharge estimation:
              Discharge = (1-theta)*Q(H0, T0) + theta*Q(H0 + delta_H, T0+delta_T))
            where H0 = headwater stage, T0 = tailwater stage, delta_H = change in
            headwater stage over a timestep, delta_T = change in tailwater stage over a
            timestep, and Q is the discharge function, and theta is a constant in
            [0,1] determining the degree of implicitness (currently hardcoded).

            Note this is effectively assuming:
            1) Q is a function of stage, not energy (so we can relate mass change directly to delta_H, delta_T). We
               could generalise it to the energy case ok.
            2) The stage is computed on the exchange line (or the change in
                stage at the enquiry point is effectively the same as that on the exchange
                line)

        """

        from anuga.utilities import parallel_abstraction as pypar

        local_debug = False
        
        # If the structure has been closed, then no water gets through
        if self.height <= 0.0:
            if self.myid == self.master_proc:
                Q = 0.0
                barrel_velocity = 0.0
                outlet_culvert_depth = 0.0
                self.case = "Structure is blocked"
                self.inflow = self.inlets[0]
                self.outflow = self.inlets[1]
                return Q, barrel_velocity, outlet_culvert_depth
            else:
                return None, None, None

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

        # Send inlet areas to the master proc. FIXME: Inlet areas don't change
        # -- perhaps we could just do this once?

        # area0
        if self.myid in self.inlet_procs[0]:
            area0 = self.inlets[0].get_global_area()

        if self.myid == self.master_proc:
            if self.myid != self.inlet_master_proc[0]:
                area0 = pypar.receive(self.inlet_master_proc[0])
        elif self.myid == self.inlet_master_proc[0]:
            pypar.send(area0, self.master_proc)
        
        # area1
        if self.myid in self.inlet_procs[1]:
            area1 = self.inlets[1].get_global_area()

        if self.myid == self.master_proc:
            if self.myid != self.inlet_master_proc[1]:
                area1 = pypar.receive(self.inlet_master_proc[1])
        elif self.myid == self.inlet_master_proc[1]:
            pypar.send(area1, self.master_proc)

        # Compute discharge
        if self.myid == self.master_proc:

            # Energy or stage as head
            if self.use_velocity_head:
                E0 = enq_total_energy0
                E1 = enq_total_energy1
            else:
                E0 = enq_stage0
                E1 = enq_stage1

            # Variables for anuga's structure operator
            self.delta_total_energy = E0 - E1
            self.driving_energy = max(E0, E1)


            Q0 = self.internal_boundary_function(E0, E1)
            dt = self.domain.get_timestep()
            
            if dt > 0.:
                # Key constants for iterative solution
                theta = 1.0
                sol = numpy.array([0., 0.]) # estimate of (delta_H, delta_T)
                areas = numpy.array([area0, area1])

                # Use scipy root finding
                def F_to_solve(sol):
                    Q1 =  self.internal_boundary_function(E0 + sol[0], E1 + sol[1])
                    discharge = (1-theta)*Q0 + theta*Q1
                    output = sol*areas - discharge*dt*numpy.array([-1., 1.])
                    return(output) 

                final_sol = sco.root(F_to_solve, sol, method='lm').x
                Q1 =  self.internal_boundary_function(E0 + final_sol[0], E1 + final_sol[1])
                Q = (1.0-theta)*Q0 + theta*Q1

            else:
                Q = Q0

            # Smooth discharge
            if dt > 0.:
                ts = dt/max(dt, self.smoothing_timescale, 1.0e-30)
            else:
                # No smoothing
                ts = 1.0

            self.smooth_Q = self.smooth_Q + ts*(Q - self.smooth_Q)

        else:
            self.delta_total_energy=numpy.nan
            self.driving_energy=numpy.nan


        self.inflow_index = 0
        self.outflow_index = 1
        # master proc orders reversal if applicable
        if self.myid == self.master_proc:

            # Reverse the inflow and outflow direction?
            if Q < 0.:
                self.inflow_index = 1
                self.outflow_index = 0

                for i in self.procs:
                    if i == self.master_proc: continue
                    pypar.send(True, i)
            else:
                for i in self.procs:
                    if i == self.master_proc: continue
                    pypar.send(False, i)

        else:
            reverse = pypar.receive(self.master_proc)

            if reverse:
                self.inflow_index = 1
                self.outflow_index = 0

        # Master proc computes return values
        if self.myid == self.master_proc:
            # Zero Q if sign's of smooth_Q and Q differ
            if numpy.sign(self.smooth_Q) != numpy.sign(Q):
                Q = 0.
                self.smooth_Q = 0.
            else:
                # Make Q positive (for anuga's structure operator)
                Q = min( abs(self.smooth_Q), abs(Q) )
            # Variables required by anuga's structure operator which are
            # not used
            barrel_velocity = numpy.nan
            outlet_culvert_depth = numpy.nan
            return Q, barrel_velocity, outlet_culvert_depth
        else:
            return None, None, None
