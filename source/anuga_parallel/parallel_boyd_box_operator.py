import anuga
import math

from anuga.structures.boyd_box_operator import boyd_box_function 

from parallel_inlet_operator import Parallel_Inlet_operator
from parallel_structure_operator import Parallel_Structure_operator

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
                 losses,
                 width,
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
                 structure_type='boyd_box',
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
                                          diameter= None,
                                          apron=apron,
                                          manning=manning,
                                          enquiry_gap=enquiry_gap,
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
        self.use_velocity_head = use_velocity_head
        
        self.culvert_length = self.get_culvert_length()
        self.culvert_width = self.get_culvert_width()
        self.culvert_height = self.get_culvert_height()

        self.max_velocity = 10.0

        self.inlets = self.get_inlets()


        # Stats
        
        self.discharge = 0.0
        self.velocity = 0.0
        
        self.case = 'N/A'

        '''
        print "ATTRIBUTES OF PARALLEL BOYD BOX::"
        for attr in dir(self):
            print "obj.%s = %s" % (attr, getattr(self, attr))
        '''


    def parallel_safe(self):

        return True


#    def debug_discharge_routine(self):
#        local_debug ='false'
#
#        if self.use_velocity_head:
#            self.delta_total_energy = self.inlets[0].get_enquiry_total_energy() - self.inlets[1].get_enquiry_total_energy()
#        else:
#            self.delta_total_energy = self.inlets[0].get_enquiry_stage() - self.inlets[1].get_enquiry_stage()
#
#        self.inflow  = self.inlets[0]
#        self.outflow = self.inlets[1]
#
#        self.inflow_index = 0
#        self.outflow_index = 1
#
#        if self.delta_total_energy < 0:
#            self.inflow  = self.inlets[1]
#            self.outflow = self.inlets[0]
#            self.delta_total_energy = -self.delta_total_energy
#            self.inflow_index = 1
#            self.outflow_index = 0
#
#
#        if self.inflow.get_enquiry_depth() > 0.01: #this value was 0.01:
#            if local_debug =='true':
#                anuga.log.critical('Specific E & Deltat Tot E = %s, %s'
#                             % (str(self.inflow.get_enquiry_specific_energy()),
#                                str(self.delta_total_energy)))
#                anuga.log.critical('culvert type = %s' % str(culvert_type))
#            # Water has risen above inlet
#
#
#            msg = 'Specific energy at inlet is negative'
#            assert self.inflow.get_enquiry_specific_energy() >= 0.0, msg
#
#            if self.use_velocity_head :
#                self.driving_energy = self.inflow.get_enquiry_specific_energy()
#            else:
#                self.driving_energy = self.inflow.get_enquiry_depth()
#
#            depth = self.culvert_height
#            width = self.culvert_width
#            flow_width = self.culvert_width
#            # intially assume the culvert flow is controlled by the inlet
#            # check unsubmerged and submerged condition and use Min Q
#            # but ensure the correct flow area and wetted perimeter are used
#            Q_inlet_unsubmerged = 0.544*anuga.g**0.5*width*self.driving_energy**1.50 # Flow based on Inlet Ctrl Inlet Unsubmerged
#            Q_inlet_submerged = 0.702*anuga.g**0.5*width*depth**0.89*self.driving_energy**0.61  # Flow based on Inlet Ctrl Inlet Submerged
#
#            # FIXME(Ole): Are these functions really for inlet control?
#            if Q_inlet_unsubmerged < Q_inlet_submerged:
#                Q = Q_inlet_unsubmerged
#                dcrit = (Q**2/anuga.g/width**2)**0.333333
#                if dcrit > depth:
#                    dcrit = depth
#                    flow_area = width*dcrit
#                    perimeter= 2.0*(width+dcrit)
#                else: # dcrit < depth
#                    flow_area = width*dcrit
#                    perimeter= 2.0*dcrit+width
#                outlet_culvert_depth = dcrit
#                self.case = 'Inlet unsubmerged Box Acts as Weir'
#            else: # Inlet Submerged but check internal culvert flow depth
#                Q = Q_inlet_submerged
#                dcrit = (Q**2/anuga.g/width**2)**0.333333
#                if dcrit > depth:
#                    dcrit = depth
#                    flow_area = width*dcrit
#                    perimeter= 2.0*(width+dcrit)
#                else: # dcrit < depth
#                    flow_area = width*dcrit
#                    perimeter= 2.0*dcrit+width
#                outlet_culvert_depth = dcrit
#                self.case = 'Inlet submerged Box Acts as Orifice'
#
#            dcrit = (Q**2/anuga.g/width**2)**0.333333
#            # May not need this .... check if same is done above
#            outlet_culvert_depth = dcrit
#            if outlet_culvert_depth > depth:
#                outlet_culvert_depth = depth  # Once again the pipe is flowing full not partfull
#                flow_area = width*depth  # Cross sectional area of flow in the culvert
#                perimeter = 2*(width+depth)
#                self.case = 'Inlet CTRL Outlet unsubmerged PIPE PART FULL'
#            else:
#                flow_area = width * outlet_culvert_depth
#                perimeter = width+2*outlet_culvert_depth
#                self.case = 'INLET CTRL Culvert is open channel flow we will for now assume critical depth'
#            # Initial Estimate of Flow for Outlet Control using energy slope
#            #( may need to include Culvert Bed Slope Comparison)
#            hyd_rad = flow_area/perimeter
#            culvert_velocity = math.sqrt(self.delta_total_energy/((self.sum_loss/2/anuga.g)+(self.manning**2*self.culvert_length)/hyd_rad**1.33333))
#            Q_outlet_tailwater = flow_area * culvert_velocity
#
#
#            if self.delta_total_energy < self.driving_energy:
#                # Calculate flows for outlet control
#
#                # Determine the depth at the outlet relative to the depth of flow in the Culvert
#                if self.outflow.get_enquiry_depth() > depth:        # The Outlet is Submerged
#                    outlet_culvert_depth=depth
#                    flow_area=width*depth       # Cross sectional area of flow in the culvert
#                    perimeter=2.0*(width+depth)
#                    self.case = 'Outlet submerged'
#                else:   # Here really should use the Culvert Slope to calculate Actual Culvert Depth & Velocity
#                    dcrit = (Q**2/anuga.g/width**2)**0.333333
#                    outlet_culvert_depth=dcrit   # For purpose of calculation assume the outlet depth = Critical Depth
#                    if outlet_culvert_depth > depth:
#                        outlet_culvert_depth=depth
#                        flow_area=width*depth
#                        perimeter=2.0*(width+depth)
#                        self.case = 'Outlet is Flowing Full'
#                    else:
#                        flow_area=width*outlet_culvert_depth
#                        perimeter=(width+2.0*outlet_culvert_depth)
#                        self.case = 'Outlet is open channel flow'
#
#                hyd_rad = flow_area/perimeter
#
#
#
#                # Final Outlet control velocity using tail water
#                culvert_velocity = math.sqrt(self.delta_total_energy/((self.sum_loss/2/anuga.g)+(self.manning**2*self.culvert_length)/hyd_rad**1.33333))
#                Q_outlet_tailwater = flow_area * culvert_velocity
#
#                Q = min(Q, Q_outlet_tailwater)
#            else:
#
#                pass
#                #FIXME(Ole): What about inlet control?
#
#            culv_froude=math.sqrt(Q**2*flow_width/(anuga.g*flow_area**3))
#
#            if local_debug =='true':
#                anuga.log.critical('FLOW AREA = %s' % str(flow_area))
#                anuga.log.critical('PERIMETER = %s' % str(perimeter))
#                anuga.log.critical('Q final = %s' % str(Q))
#                anuga.log.critical('FROUDE = %s' % str(culv_froude))
#
#            # Determine momentum at the outlet
#            barrel_velocity = Q/(flow_area + anuga.velocity_protection/flow_area)
#
#        # END CODE BLOCK for DEPTH  > Required depth for CULVERT Flow
#
#        else: # self.inflow.get_enquiry_depth() < 0.01:
#            Q = barrel_velocity = outlet_culvert_depth = 0.0
#
#        # Temporary flow limit
#        if barrel_velocity > self.max_velocity:
#            barrel_velocity = self.max_velocity
#            Q = flow_area * barrel_velocity
#
#        return Q, barrel_velocity, outlet_culvert_depth

    def discharge_routine(self):

        import pypar

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


            # Reverse the inflow and outflow direction?
            if self.delta_total_energy < 0:
                self.inflow_index = 1
                self.outflow_index = 0

                self.delta_total_energy = -self.delta_total_energy

                for i in self.procs:
                    if i == self.master_proc: continue
                    pypar.send(True, i)
            else:
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
                    
                #print "ZZZZZ: driving energy = %f" %(self.driving_energy)

#                depth = self.culvert_height
#                width = self.culvert_width
#                flow_width = self.culvert_width
#                driving_energy = self.driving_energy
#                sum_loss = self.sum_loss
#                culvert_length= self.culvert_length
#                manning = self.manning
#                delta_total_energy = self.delta_total_energy


            
                Q, barrel_velocity, outlet_culvert_depth, flow_area, case = \
                              boyd_box_function(depth               =self.culvert_height,
                                                width               =self.culvert_width,
                                                flow_width          =self.culvert_width,
                                                length              =self.culvert_length,
                                                driving_energy      =self.driving_energy,
                                                delta_total_energy  =self.delta_total_energy,
                                                outlet_enquiry_depth=outflow_enq_depth,
                                                sum_loss            =self.sum_loss,
                                                manning             =self.manning)

                
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
        
        
