import anuga
import math
import numpy

#from anuga.structures.boyd_box_operator import boyd_box_function 

from parallel_inlet_operator import Parallel_Inlet_operator
from parallel_structure_operator import Parallel_Structure_operator

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
                 smoothing_timescale=0.0,
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
            print '########################################'
            print 'PARALLEL INTERNAL BOUNDARY OPERATOR'
            print 'THIS IS EXPERIMENTAL'
            print 'SUBJECT TO CHANGE WITHOUT NOTICE'
            print '########################################'

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
                                          enquiry_gap=enquiry_gap,
                                          use_momentum_jet=use_momentum_jet,
                                          zero_outflow_momentum=zero_outflow_momentum,
                                          use_old_momentum_method=False,
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
        # Set them based on a call to the discharge routine with smoothing_timescale=0.
        # [values of self.smooth_* are required in discharge_routine, hence dummy values above]
        Qvd = self.discharge_routine()
        self.smooth_Q = Qvd[0]
        # Finally, set the smoothing timescale we actually want
        self.smoothing_timescale = smoothing_timescale

    def parallel_safe(self):

        return True



    def discharge_routine(self):

        import pypar

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
            if self.use_velocity_head:
                self.delta_total_energy = enq_total_energy0 - enq_total_energy1
                self.driving_energy = max(enq_total_energy0, enq_total_energy1)
                # Compute discharge
                Q = self.internal_boundary_function(enq_total_energy0, enq_total_energy1)
            else:
                self.delta_total_energy = enq_stage0 - enq_stage1
                self.driving_energy = max(enq_stage0, enq_stage1)
                # Compute discharge
                Q = self.internal_boundary_function(enq_stage0, enq_stage1)

            # Other variables required by anuga's structure operator are not used
            barrel_velocity = numpy.nan
            outlet_culvert_depth = numpy.nan
            flow_area = numpy.nan
            case = ''

            # Use time-smoothed discharge
            ts = self.domain.timestep/max(self.domain.timestep, self.smoothing_timescale, 1.0e-30)
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
        
        
