
import anuga
import numpy as num
import math
from . import parallel_inlet_enquiry 
from anuga.utilities import parallel_abstraction as pypar

from anuga.utilities.system_tools import log_to_file
from anuga.utilities.numerical_tools import ensure_numeric
from anuga.structures.inlet_enquiry import Inlet_enquiry


class Parallel_Structure_operator(anuga.Operator):
    """Parallel Structure Operator - transfer water from one rectangular box to another.
    Sets up the geometry of problem
    
    This is the base class for structures (culverts, pipes, bridges etc) that exist across multiple
    parallel shallow water domains. Inherit from this class (and overwrite discharge_routine method 
    for specific subclasses)
    
    Input: Two points, pipe_size (either diameter or width, depth),
    mannings_rougness,
    """ 

    counter = 0

    """
     ===========================================================================================
     PETE: Inputs to this constructor are identical to the serial
     structure operator, except for the following arguments:
     master_proc - the processor that coordinates all processors (with domains) associated with this structure [INT]
     procs - the list of processors associated with thisstructure (List[INT])
     inlet_master_proc - master_proc of the first and second inlet (List[2])
     inlet_procs - list of processors associated with the first and second inlet (LIST[2][INT])
     enquiry_proc - processor associated the first and second enquiry point (List[2])
    """

    def __init__(self,
                 domain,
                 end_points,
                 exchange_lines,
                 enquiry_points,
                 invert_elevations,
                 width,
                 height,
                 diameter,
                 z1,
                 z2,
                 blockage,
                 barrels,
                 apron,
                 manning,
                 enquiry_gap,
                 use_momentum_jet,
                 zero_outflow_momentum,
                 use_old_momentum_method,
                 always_use_Q_wetdry_adjustment,
                 force_constant_inlet_elevations,
                 description,
                 label,
                 structure_type,
                 logging,
                 verbose,
                 master_proc = 0,
                 procs = None,
                 inlet_master_proc = [0,0],
                 inlet_procs = None,
                 enquiry_proc = None):


        self.myid = pypar.rank()
        self.num_procs = pypar.size()
        
        anuga.Operator.__init__(self,domain)

        # Allocate default processor associations if not specified in arguments
        # although we assume that such associations are provided correctly by the 
        # parallel_operator_factory.

        self.master_proc = master_proc
        self.inlet_master_proc = inlet_master_proc
        
        if procs is None:
            self.procs = [master_proc]
        else:
            self.procs = procs

        if inlet_procs is None:
            self.inlet_procs = [[inlet_master_proc[0]],[inlet_master_proc[0]]]
        else:
            self.inlet_procs = inlet_procs

        if enquiry_proc is None:
            self.enquiry_proc = [[inlet_master_proc[0]],[inlet_master_proc[0]]]
        else:
            self.enquiry_proc = enquiry_proc

        self.end_points = ensure_numeric(end_points)
        self.exchange_lines = ensure_numeric(exchange_lines)
        self.enquiry_points = ensure_numeric(enquiry_points)
        self.invert_elevations = ensure_numeric(invert_elevations)

        assert (width is not None and diameter is None) or (width is None and diameter is not None)

        if width is None:
            width = diameter

        if diameter is None:
            diameter = width

        if height is None:
            height = width

        if apron is None:
            apron = width

        self.width  = width
        self.height = height
        self.diameter = diameter
        self.z1 = z1
        self.z2 = z2
        self.blockage = blockage
        self.barrels = barrels
        self.apron  = apron
        self.manning = manning
        self.enquiry_gap = enquiry_gap
        self.use_momentum_jet = use_momentum_jet
        self.zero_outflow_momentum = zero_outflow_momentum
        if use_momentum_jet and zero_outflow_momentum:
            msg = "Can't have use_momentum_jet and zero_outflow_momentum both True"
            raise Exception(msg)
        self.use_old_momentum_method = use_old_momentum_method
        self.always_use_Q_wetdry_adjustment = always_use_Q_wetdry_adjustment

        if description is None:
            self.description = ' '
        else:
            self.description = description
        
        if label is None:
            self.label = "structure_%g" % Parallel_Structure_operator.counter + "_P" + str(self.myid)
        else:
            self.label = label + '_%g' % Parallel_Structure_operator.counter + "_P" + str(self.myid)

        if structure_type is None:
            self.structure_type = 'generic structure'
        else:
            self.structure_type = structure_type
            
        self.verbose = verbose        
        
        # Keep count of structures
        if self.myid == master_proc:
            Parallel_Structure_operator.counter += 1

        # Slots for recording current statistics
        self.accumulated_flow = 0.0
        self.discharge = 0.0
        self.discharge_abs_timemean = 0.0
        self.velocity = 0.0
        self.outlet_depth = 0.0
        self.delta_total_energy = 0.0
        self.driving_energy = 0.0
        
        if exchange_lines is not None:
            self.__process_skew_culvert()
        elif end_points is not None:
            self.__process_non_skew_culvert()
        else:
            raise Exception('Define either exchange_lines or end_points')
        
        self.inlets = []

        # Allocate parallel inlet enquiry, assign None if processor is not associated with particular
        # inlet.

        if self.myid in self.inlet_procs[0]:
            line0 = self.exchange_lines[0]
            if self.apron is None:
                poly0 = line0
            else:
                offset = -self.apron*self.outward_vector_0 
                poly0 = num.array([ line0[0], line0[1], line0[1]+offset, line0[0]+offset])

            if self.invert_elevations is None:
                invert_elevation0 = None
            else:
                invert_elevation0 = self.invert_elevations[0]

            enquiry_point0 = self.enquiry_points[0]
            outward_vector0 = self.culvert_vector

            self.inlets.append(parallel_inlet_enquiry.Parallel_Inlet_enquiry(
                               self.domain,
                               line0,
                               enquiry_point0,
                               invert_elevation = invert_elevation0,
                               outward_culvert_vector = outward_vector0, 
                               master_proc = self.inlet_master_proc[0],
                               procs = self.inlet_procs[0],
                               enquiry_proc = self.enquiry_proc[0],
                               verbose = self.verbose))

            if force_constant_inlet_elevations: 
                # Try to enforce a constant inlet elevation 
                inlet_global_elevation = self.inlets[-1].get_global_average_elevation() 
                self.inlets[-1].set_elevations(inlet_global_elevation)
       
        else:
            self.inlets.append(None)

        if self.myid in self.inlet_procs[1]:
            line1 = self.exchange_lines[1]
            if self.apron is None:
                poly1 = line1
            else:
                offset = -self.apron*self.outward_vector_1
                poly1 = num.array([ line1[0], line1[1], line1[1]+offset, line1[0]+offset])

            if self.invert_elevations is None:
                invert_elevation1 = None
            else:
                invert_elevation1 = self.invert_elevations[1]


            enquiry_point1 = self.enquiry_points[1]
            outward_vector1  = - self.culvert_vector


            self.inlets.append(parallel_inlet_enquiry.Parallel_Inlet_enquiry(
                               self.domain,
                               line1,
                               enquiry_point1,
                               invert_elevation = invert_elevation1,
                               outward_culvert_vector = outward_vector1,
                               master_proc = self.inlet_master_proc[1],
                               procs = self.inlet_procs[1],
                               enquiry_proc = self.enquiry_proc[1],
                               verbose = self.verbose))

            if force_constant_inlet_elevations: 
                # Try to enforce a constant inlet elevation 
                inlet_global_elevation = self.inlets[-1].get_global_average_elevation() 
                self.inlets[-1].set_elevations(inlet_global_elevation)
            

        else:
            self.inlets.append(None)

        self.inflow_index = 0
        self.outflow_index = 1

        self.set_parallel_logging(logging)

    def __call__(self):

        timestep = self.domain.get_timestep()

        Q, barrel_speed, outlet_depth = self.discharge_routine()

        # Get attributes of Inflow inlet, all procs associated with inlet must call
        if self.myid in self.inlet_procs[self.inflow_index]:
            old_inflow_depth = self.inlets[self.inflow_index].get_global_average_depth()
            old_inflow_stage = self.inlets[self.inflow_index].get_global_average_stage()
            old_inflow_xmom = self.inlets[self.inflow_index].get_global_average_xmom()
            old_inflow_ymom = self.inlets[self.inflow_index].get_global_average_ymom()
            inflow_area = self.inlets[self.inflow_index].get_global_area()

        # Master proc of inflow inlet sends attributes to master proc of structure
        if self.myid == self.master_proc:
            if self.myid != self.inlet_master_proc[self.inflow_index]:
                old_inflow_depth = pypar.receive(self.inlet_master_proc[self.inflow_index])
                old_inflow_stage = pypar.receive(self.inlet_master_proc[self.inflow_index])
                old_inflow_xmom = pypar.receive(self.inlet_master_proc[self.inflow_index])
                old_inflow_ymom = pypar.receive(self.inlet_master_proc[self.inflow_index])
                inflow_area = pypar.receive(self.inlet_master_proc[self.inflow_index])
        elif self.myid == self.inlet_master_proc[self.inflow_index]:
            pypar.send(old_inflow_depth, self.master_proc)
            pypar.send(old_inflow_stage, self.master_proc)
            pypar.send(old_inflow_xmom, self.master_proc)
            pypar.send(old_inflow_ymom, self.master_proc)
            pypar.send(inflow_area, self.master_proc)

        # Implement the update of flow over a timestep by
        # using a semi-implict update. This ensures that
        # the update does not create a negative depth
        
        # Master proc of structure only
        if self.myid == self.master_proc:
            if old_inflow_depth > 0.0 :
                dt_Q_on_d = timestep*Q/old_inflow_depth
            else:
                dt_Q_on_d = 0.0

            # Check whether we should use the wet-dry Q adjustment (where Q is
            # multiplied by new_inflow_depth/old_inflow_depth)
            always_use_Q_wetdry_adjustment = self.always_use_Q_wetdry_adjustment
            # Always use it if we are near wet-dry
            use_Q_wetdry_adjustment = ((always_use_Q_wetdry_adjustment) |\
                (old_inflow_depth*inflow_area <= Q*timestep))

            factor = 1.0/(1.0 + dt_Q_on_d/inflow_area)
        
            if use_Q_wetdry_adjustment:
                new_inflow_depth = old_inflow_depth*factor
                if old_inflow_depth > 0.:
                    timestep_star = timestep*new_inflow_depth/old_inflow_depth
                else:
                    timestep_star = 0.
            else:
                new_inflow_depth = old_inflow_depth - timestep*Q/inflow_area
                timestep_star = timestep

            #new_inflow_xmom = old_inflow_xmom*factor
            #new_inflow_ymom = old_inflow_ymom*factor
            if(self.use_old_momentum_method):
                # This method is here for consistency with the old version of the
                # routine
                new_inflow_xmom = old_inflow_xmom*factor
                new_inflow_ymom = old_inflow_ymom*factor

            else:
                # For the momentum balance, note that Q also transports the velocity,
                # which has an average value of new_inflow_mom/depth (or old_inflow_mom/depth). 
                #
                #     new_inflow_xmom*inflow_area = 
                #     old_inflow_xmom*inflow_area - 
                #     timestep*Q*(new_inflow_xmom/old_inflow_depth)
                # and:
                #     new_inflow_ymom*inflow_area = 
                #     old_inflow_ymom*inflow_area - 
                #     timestep*Q*(new_inflow_ymom/old_inflow_depth)
                #
                # The choice of new_inflow_mom in the final term might be
                # replaced with old_inflow_mom.
                #
                # The units balance: (m^2/s)*(m^2) = (m^2/s)*(m^2) - s*(m^3/s)*(m^2/s)*(m^(-1))
                #
                if old_inflow_depth > 0.:
                    if use_Q_wetdry_adjustment:
                        factor2 = 1.0/(1.0 + dt_Q_on_d*new_inflow_depth/(old_inflow_depth*inflow_area))
                    else:
                        factor2 = 1.0/(1.0 + timestep*Q/(old_inflow_depth*inflow_area))
                else:
                    factor2 = 0.

                new_inflow_xmom = old_inflow_xmom*factor2
                new_inflow_ymom = old_inflow_ymom*factor2

        # Master proc of structure sends new inflow attributes to all inflow inlet processors

        if self.myid == self.master_proc:
            for i in self.inlet_procs[self.inflow_index]:
                if i == self.master_proc: continue
                pypar.send(new_inflow_depth, i)
                pypar.send(new_inflow_xmom, i)
                pypar.send(new_inflow_ymom, i)
        elif self.myid in self.inlet_procs[self.inflow_index]:
            new_inflow_depth = pypar.receive(self.master_proc)
            new_inflow_xmom = pypar.receive(self.master_proc)
            new_inflow_ymom = pypar.receive(self.master_proc)

        # Inflow inlet procs sets new attributes
        if self.myid in self.inlet_procs[self.inflow_index]:
            self.inlets[self.inflow_index].set_depths(new_inflow_depth)
            self.inlets[self.inflow_index].set_xmoms(new_inflow_xmom)
            self.inlets[self.inflow_index].set_ymoms(new_inflow_ymom)

        # Get outflow inlet attributes, all processors associated with outflow inlet must call
        if self.myid in self.inlet_procs[self.outflow_index]:
            outflow_area = self.inlets[self.outflow_index].get_global_area()
            outflow_average_depth = self.inlets[self.outflow_index].get_global_average_depth()
            outflow_outward_culvert_vector = self.inlets[self.outflow_index].outward_culvert_vector
            outflow_average_xmom = self.inlets[self.outflow_index].get_global_average_xmom()
            outflow_average_ymom = self.inlets[self.outflow_index].get_global_average_ymom()

        # Master proc of outflow inlet sends attribute to master proc of structure
        if self.myid == self.master_proc:
            if self.myid != self.inlet_master_proc[self.outflow_index]:
                outflow_area = pypar.receive(self.inlet_master_proc[self.outflow_index])
                outflow_average_depth = pypar.receive(self.inlet_master_proc[self.outflow_index])
                outflow_outward_culvert_vector = pypar.receive(self.inlet_master_proc[self.outflow_index])
                outflow_average_xmom = pypar.receive(self.inlet_master_proc[self.outflow_index])
                outflow_average_ymom = pypar.receive(self.inlet_master_proc[self.outflow_index])
        elif self.myid == self.inlet_master_proc[self.outflow_index]:
            pypar.send(outflow_area, self.master_proc)
            pypar.send(outflow_average_depth, self.master_proc)
            pypar.send(outflow_outward_culvert_vector, self.master_proc)
            pypar.send(outflow_average_xmom, self.master_proc)
            pypar.send(outflow_average_ymom, self.master_proc)

        # Master proc of structure computes new outflow attributes
        if self.myid == self.master_proc:
            loss = (old_inflow_depth - new_inflow_depth)*inflow_area
            xmom_loss = (old_inflow_xmom - new_inflow_xmom)*inflow_area
            ymom_loss = (old_inflow_ymom - new_inflow_ymom)*inflow_area

            # set outflow
            outflow_extra_depth = Q*timestep_star/outflow_area
            outflow_direction = - outflow_outward_culvert_vector
            #outflow_extra_momentum = outflow_extra_depth*barrel_speed*outflow_direction
            
            gain = outflow_extra_depth*outflow_area

            #print('gain', gain)

            # Update Stats
            self.accumulated_flow += gain
            self.discharge  = Q*timestep_star/timestep #outflow_extra_depth*self.outflow.get_area()/timestep
            self.discharge_abs_timemean += Q*timestep_star/self.domain.yieldstep
            self.velocity = barrel_speed #self.discharge/outlet_depth/self.width

            new_outflow_depth = outflow_average_depth + outflow_extra_depth

            self.outlet_depth = new_outflow_depth
            #if self.use_momentum_jet :
            #    # FIXME (SR) Review momentum to account for possible hydraulic jumps at outlet
            #    #new_outflow_xmom = outflow.get_average_xmom() + outflow_extra_momentum[0]
            #    #new_outflow_ymom = outflow.get_average_ymom() + outflow_extra_momentum[1]

            #    new_outflow_xmom = barrel_speed*new_outflow_depth*outflow_direction[0]
            #    new_outflow_ymom = barrel_speed*new_outflow_depth*outflow_direction[1]

            #else:
            #    #new_outflow_xmom = outflow.get_average_xmom()
            #    #new_outflow_ymom = outflow.get_average_ymom()

            #    new_outflow_xmom = 0.0
            #    new_outflow_ymom = 0.0
            if self.use_momentum_jet:
                # FIXME (SR) Review momentum to account for possible hydraulic jumps at outlet
                # FIXME (GD) Depending on barrel speed I think this will be either
                # a source or sink of momentum (considering the momentum losses
                # above). Might not always be reasonable.
                #new_outflow_xmom = self.outflow.get_average_xmom() + outflow_extra_momentum[0]
                #new_outflow_ymom = self.outflow.get_average_ymom() + outflow_extra_momentum[1]
                new_outflow_xmom = barrel_speed*new_outflow_depth*outflow_direction[0]
                new_outflow_ymom = barrel_speed*new_outflow_depth*outflow_direction[1]
                
            elif self.zero_outflow_momentum:
                new_outflow_xmom = 0.0
                new_outflow_ymom = 0.0
                #new_outflow_xmom = outflow.get_average_xmom()
                #new_outflow_ymom = outflow.get_average_ymom()

            else:
                # Add the momentum lost from the inflow to the outflow. For
                # structures where barrel_speed is unknown + direction doesn't
                # change from inflow to outflow
                new_outflow_xmom = outflow_average_xmom + xmom_loss/outflow_area
                new_outflow_ymom = outflow_average_ymom + ymom_loss/outflow_area

            # master proc of structure sends outflow attributes to all outflow procs
            for i in self.inlet_procs[self.outflow_index]:
                if i == self.myid: continue
                pypar.send(new_outflow_depth, i)
                pypar.send(new_outflow_xmom, i)
                pypar.send(new_outflow_ymom, i)
        # outflow inlet procs receives new outflow attributes
        elif self.myid in self.inlet_procs[self.outflow_index]:
            new_outflow_depth = pypar.receive(self.master_proc)
            new_outflow_xmom = pypar.receive(self.master_proc)
            new_outflow_ymom = pypar.receive(self.master_proc)

        # outflow inlet procs sets new outflow attributes
        if self.myid in self.inlet_procs[self.outflow_index]:
            self.inlets[self.outflow_index].set_depths(new_outflow_depth)
            self.inlets[self.outflow_index].set_xmoms(new_outflow_xmom)
            self.inlets[self.outflow_index].set_ymoms(new_outflow_ymom)

    def __process_non_skew_culvert(self):
        """Create lines at the end of a culvert inlet and outlet.
        At either end two lines will be created; one for the actual flow to pass through and one a little further away
        for enquiring the total energy at both ends of the culvert and transferring flow.
        """
        
        self.culvert_vector = self.end_points[1] - self.end_points[0]
        self.culvert_length = math.sqrt(num.sum(self.culvert_vector**2))   
        assert self.culvert_length > 0.0, 'The length of culvert is less than 0'
        
        self.culvert_vector /= self.culvert_length
        self.outward_vector_0 =   self.culvert_vector
        self.outward_vector_1 = - self.culvert_vector        
        
        culvert_normal = num.array([-self.culvert_vector[1], self.culvert_vector[0]])  # Normal vector
        w = 0.5*self.width*culvert_normal # Perpendicular vector of 1/2 width

        self.exchange_lines = []

        # Build exchange polyline and enquiry point
        if self.enquiry_points is None:
            
            gap = (self.apron + self.enquiry_gap)*self.culvert_vector
            self.enquiry_points = []
            
            for i in [0, 1]:
                p0 = self.end_points[i] + w
                p1 = self.end_points[i] - w
                self.exchange_lines.append(num.array([p0, p1]))
                ep = self.end_points[i] + (2*i - 1)*gap #(2*i - 1) determines the sign of the points
                self.enquiry_points.append(ep)
            
        else:            
            for i in [0, 1]:
                p0 = self.end_points[i] + w
                p1 = self.end_points[i] - w
                self.exchange_lines.append(num.array([p0, p1]))
            
  
    def __process_skew_culvert(self):    
        
        """Compute skew culvert.
        If exchange lines are given, the enquiry points are determined. This is for enquiring 
        the total energy at both ends of the culvert and transferring flow.
        """
            
        centre_point0 = 0.5*(self.exchange_lines[0][0] + self.exchange_lines[0][1])
        centre_point1 = 0.5*(self.exchange_lines[1][0] + self.exchange_lines[1][1])

        n_exchange_0 = len(self.exchange_lines[0])
        n_exchange_1 = len(self.exchange_lines[1])

        assert n_exchange_0 == n_exchange_1, 'There should be the same number of points in both exchange_lines'

        if n_exchange_0 == 2:
            
            if self.end_points is None:
                self.culvert_vector = centre_point1 - centre_point0
            else:
                self.culvert_vector = self.end_points[1] - self.end_points[0]

            self.outward_vector_0 =   self.culvert_vector
            self.outward_vector_1 = - self.culvert_vector

        elif n_exchange_0 == 4:

            self.outward_vector_0 = self.exchange_lines[0][3] - self.exchange_lines[0][2]
            self.outward_vector_1 = self.exchange_lines[1][3] - self.exchange_lines[1][2]

            self.culvert_vector = centre_point1 - centre_point0

        else:
            raise Exception('n_exchange_0 != 2 or 4')

        self.culvert_length = math.sqrt(num.sum(self.culvert_vector**2))
        assert self.culvert_length > 0.0, 'The length of culvert is less than 0'
        self.culvert_vector /= self.culvert_length

        outward_vector_0_length = math.sqrt(num.sum(self.outward_vector_0**2))
        assert outward_vector_0_length > 0.0, 'The length of outlet_vector_0 is less than 0'
        self.outward_vector_0 /= outward_vector_0_length

        outward_vector_1_length = math.sqrt(num.sum(self.outward_vector_1**2))
        assert outward_vector_1_length > 0.0, 'The length of outlet_vector_1 is less than 0'
        self.outward_vector_1 /= outward_vector_1_length

        
        if self.enquiry_points is None:
        

            gap = (self.apron + self.enquiry_gap)*self.culvert_vector
        
            self.enquiry_points = []

            self.enquiry_points.append(centre_point0 - gap)
            self.enquiry_points.append(centre_point1 + gap)
            

    def discharge_routine(self):

        msg = 'Need to impelement '
        raise
            

    def statistics(self):
        # Warning: requires synchronization, must be called by all procs associated
        # with this structure

        message = ' '

        if self.myid == self.master_proc:

            message  = '===============================================\n'
            message += 'Parallel Structure Operator: %s\n' % self.label
            message += '===============================================\n'

            message += 'Structure Type: %s\n' % self.structure_type

            message += 'Description\n'
            message += '%s' % self.description
            message += '\n'

            #add the culvert dimensions, blockage factor here
            if self.structure_type == 'boyd_pipe':
                message += 'Culvert Diameter: %s\n'% self.diameter
                message += 'Culvert Blockage: %s\n'% self.blockage
                message += 'No.  of  barrels: %s\n'% self.barrels
            elif self.structure_type == 'boyd_box':
                message += 'Culvert   Height: %s\n'% self.height
                message += 'Culvert    Width: %s\n'% self.width
                message += 'Culvert Blockage: %s\n'% self.blockage
                message += 'No.  of  barrels: %s\n'% self.barrels
            else:
                message += 'Culvert Height  : %s\n'% self.height
                message += 'Culvert  Width  : %s\n'% self.width
                message += 'Batter Slope 1  : %s\n'% self.z1
                message += 'Batter Slope 2  : %s\n'% self.z2
                message += 'Culvert Blockage: %s\n'% self.blockage
                message += 'No.  of  barrels: %s\n'% self.barrels
                
        #print "Structure Myids ",self.myid, self.label
        
        for i, inlet in enumerate(self.inlets):
            if self.myid == self.master_proc:
                message += '-------------------------------------\n'
                message +=  'Inlet %i\n' %(i)
                message += '-------------------------------------\n'

            #print "*****",inlet, i,self.myid
            if inlet is not None:
                
                
                stats = inlet.statistics()

            if self.myid == self.master_proc:
                if self.myid != self.inlet_master_proc[i]:
                    stats = pypar.receive(self.inlet_master_proc[i])                    
            elif self.myid == self.inlet_master_proc[i]:
                pypar.send(stats, self.master_proc)

            if self.myid == self.master_proc: message += stats
 

        if self.myid == self.master_proc: message += '=====================================\n'

        return message


    def print_statistics(self):
        # Warning: requires synchronization, must be called by all procs associated
        # with this structure

        print(self.statistics())


    def print_timestepping_statistics(self):
        # Warning: must be called by the master proc of this structure to obtain 
        # meaningful output

        message = ' '

        if self.myid == self.master_proc:
            message = '--------------------------------------------------\n'
            message += 'Parallel Structure report for %s:\n' % self.label
            message += '-------------------------------------------------\n'
            message += 'Type: %s\n' % self.structure_type
            message += 'Discharge [m^3/s]: %.2f\n' % self.discharge
            message += 'Discharge function value [m^3/s]: %.2f\n' % self.discharge_abs_timemean
            message += 'Velocity  [m/s]: %.2f\n' % self.velocity
            message += 'Inlet Driving Energy %.2f\n' % self.driving_energy
            message += 'Delta Total Energy %.2f\n' % self.delta_total_energy
            message += 'Control at this instant: %s\n' % self.case

        print(message)


    def set_parallel_logging(self, flag=True):
        # Warning: requires synchronization, must be called by all procs associated
        # with this structure

        stats = self.statistics()
        self.logging = flag

        # If flag is true open file with mode = "w" to form a clean file for logging
        if self.logging and self.myid == self.master_proc:
            self.log_filename = self.domain.get_datadir() + '/' + self.label + '.log'
            log_to_file(self.log_filename, stats, mode='w')
            log_to_file(self.log_filename, 'time,discharge_instantaneous,discharge_abs_timemean,velocity_instantaneous,driving_energy_instantaneous,delta_total_energy_instantaneous')

            #log_to_file(self.log_filename, self.culvert_type)

    def set_logging(self, flag=True):
        # Overwrite the sequential procedure with a dummy procedure.
        # Need to call set_parallel_logging which needs to be done later
        # after the calculation of master processors

        pass


    def log_timestepping_statistics(self):

        from anuga.utilities.system_tools import log_to_file
        if self.logging and self.myid == self.master_proc:
            log_to_file(self.log_filename, self.timestepping_statistics())



    def timestepping_statistics(self):

        message  = '%.5f, ' % self.domain.get_time()
        message += '%.5f, ' % self.discharge
        message += '%.5f, ' % self.discharge_abs_timemean
        message += '%.5f, ' % self.velocity
        message += '%.5f, ' % self.driving_energy
        message += '%.5f' % self.delta_total_energy

        # Reset discharge_abs_timemean each time there is reporting (FIXME:
        # This assumes this function is only called after each yieldstep)
        self.discharge_abs_timemean = 0.

        return message


    def get_inlets(self):
        return self.inlets
        
        
    def get_culvert_length(self):
        return self.culvert_length


    def get_culvert_width(self):        
        return self.width
        
        
    def get_culvert_diameter(self):
        return self.diameter
        
        
    def get_culvert_height(self):
        return self.height

    def get_culvert_z1(self):
        return self.z1

    def get_culvert_z2(self):   
        return self.z2

    def get_culvert_blockage(self):	
        return self.blockage

    def get_culvert_barrels(self):	
        return self.barrels
                        
    def get_culvert_apron(self):
        return self.apron

    # Get id of master proc of this structure
    def get_master_proc(self):
        return self.master_proc

    # Get id of master proc of first and second inlet
    def get_inlet_master_proc(self):
        return self.inlet_master_proc

    # Get id of processors associated with first and second inlet enquiry points
    def get_enquiry_proc(self, id=None):

        if id is none:
            return self.enquiry_proc
        else:
            return self.enquiry_proc[id]


    def set_culvert_height(self, height):

        self.culvert_height = height

    def set_culvert_width(self, width):

        self.culvert_width = width

    def set_culvert_z1(self, z1):

        self.culvet_z1 = z1        

    def set_culvert_z2(self, z2):

        self.culvert_z2 = z2  

    def set_culvert_blockage(self, blockage): 

        self.culvert_blockage = blockage 

    def set_culvert_blockage(self, barrels): 

        self.culvert_barrels = barrels 
        
                        
    def parallel_safe(self):
        return True


    def get_enquiry_stages(self):
        # Should be called from all processors associated with operator

        enq0 = None
        enq1 = None

        get0 = 'self.inlets[0].get_enquiry_stage()'
        get1 = 'self.inlets[1].get_enquiry_stage()'

        if self.myid == self.master_proc:

            if self.myid == self.enquiry_proc[0]:
                enq0 = eval(get0)
            else:
                enq0 = pypar.receive(self.enquiry_proc[0])


            if self.myid == self.enquiry_proc[1]:
                enq1 = eval(get1)
            else:
                enq1 = pypar.receive(self.enquiry_proc[1])

        else:
            if self.myid == self.enquiry_proc[0]:
                enq0 = eval(get0)
                pypar.send(enq0, self.master_proc)

            if self.myid == self.enquiry_proc[1]:
                enq1 = eval(get1)
                pypar.send(enq1, self.master_proc)


        return [enq0, enq1]

    def get_enquiry_depths(self):
        # Should be called from all processors associated with operator

        enq0 = None
        enq1 = None

        get0 = 'self.inlets[0].get_enquiry_depth()'
        get1 = 'self.inlets[1].get_enquiry_depth()'

        if self.myid == self.master_proc:

            if self.myid == self.enquiry_proc[0]:
                enq0 = eval(get0)
            else:
                enq0 = pypar.receive(self.enquiry_proc[0])


            if self.myid == self.enquiry_proc[1]:
                enq1 = eval(get1)
            else:
                enq1 = pypar.receive(self.enquiry_proc[1])

        else:
            if self.myid == self.enquiry_proc[0]:
                enq0 = eval(get0)
                pypar.send(enq0, self.master_proc)

            if self.myid == self.enquiry_proc[1]:
                enq1 = eval(get1)
                pypar.send(enq1, self.master_proc)


        return [enq0, enq1]



    def get_enquiry_positions(self):

        enq0 = None
        enq1 = None

        get0 = 'self.inlets[0].get_enquiry_position()'
        get1 = 'self.inlets[1].get_enquiry_position()'


        if self.myid == self.master_proc:

            if self.myid == self.enquiry_proc[0]:
                enq0 = eval(get0)
            else:
                enq0 = pypar.receive(self.enquiry_proc[0])


            if self.myid == self.enquiry_proc[1]:
                enq1 = eval(get1)
            else:
                enq1 = pypar.receive(self.enquiry_proc[1])

        else:
            if self.myid == self.enquiry_proc[0]:
                enq0 = eval(get0)
                pypar.send(enq0, self.master_proc)

            if self.myid == self.enquiry_proc[1]:
                enq1 = eval(get1)
                pypar.send(enq1, self.master_proc)


        return [enq0, enq1]


    def get_enquiry_xmoms(self):

        enq0 = None
        enq1 = None

        get0 = 'self.inlets[0].get_enquiry_xmom()'
        get1 = 'self.inlets[1].get_enquiry_xmom()'


        if self.myid == self.master_proc:

            if self.myid == self.enquiry_proc[0]:
                enq0 = eval(get0)
            else:
                enq0 = pypar.receive(self.enquiry_proc[0])


            if self.myid == self.enquiry_proc[1]:
                enq1 = eval(get1)
            else:
                enq1 = pypar.receive(self.enquiry_proc[1])

        else:
            if self.myid == self.enquiry_proc[0]:
                enq0 = eval(get0)
                pypar.send(enq0, self.master_proc)

            if self.myid == self.enquiry_proc[1]:
                enq1 = eval(get1)
                pypar.send(enq1, self.master_proc)


        return [enq0, enq1]

    def get_enquiry_ymoms(self):

        enq0 = None
        enq1 = None

        get0 = 'self.inlets[0].get_enquiry_ymom()'
        get1 = 'self.inlets[1].get_enquiry_ymom()'


        if self.myid == self.master_proc:

            if self.myid == self.enquiry_proc[0]:
                enq0 = eval(get0)
            else:
                enq0 = pypar.receive(self.enquiry_proc[0])


            if self.myid == self.enquiry_proc[1]:
                enq1 = eval(get1)
            else:
                enq1 = pypar.receive(self.enquiry_proc[1])

        else:
            if self.myid == self.enquiry_proc[0]:
                enq0 = eval(get0)
                pypar.send(enq0, self.master_proc)

            if self.myid == self.enquiry_proc[1]:
                enq1 = eval(get1)
                pypar.send(enq1, self.master_proc)


        return [enq0, enq1]


    def get_enquiry_elevations(self):

        enq0 = None
        enq1 = None

        get0 = 'self.inlets[0].get_enquiry_elevation()'
        get1 = 'self.inlets[1].get_enquiry_elevation()'


        if self.myid == self.master_proc:

            if self.myid == self.enquiry_proc[0]:
                enq0 = eval(get0)
            else:
                enq0 = pypar.receive(self.enquiry_proc[0])


            if self.myid == self.enquiry_proc[1]:
                enq1 = eval(get1)
            else:
                enq1 = pypar.receive(self.enquiry_proc[1])

        else:
            if self.myid == self.enquiry_proc[0]:
                enq0 = eval(get0)
                pypar.send(enq0, self.master_proc)

            if self.myid == self.enquiry_proc[1]:
                enq1 = eval(get1)
                pypar.send(enq1, self.master_proc)


        return [enq0, enq1]



    def get_enquiry_water_depths(self):

        enq0 = None
        enq1 = None

        get0 = 'self.inlets[0].get_enquiry_water_depth()'
        get1 = 'self.inlets[1].get_enquiry_water_depth()'


        if self.myid == self.master_proc:

            if self.myid == self.enquiry_proc[0]:
                enq0 = eval(get0)
            else:
                enq0 = pypar.receive(self.enquiry_proc[0])


            if self.myid == self.enquiry_proc[1]:
                enq1 = eval(get1)
            else:
                enq1 = pypar.receive(self.enquiry_proc[1])

        else:
            if self.myid == self.enquiry_proc[0]:
                enq0 = eval(get0)
                pypar.send(enq0, self.master_proc)

            if self.myid == self.enquiry_proc[1]:
                enq1 = eval(get1)
                pypar.send(enq1, self.master_proc)


        return [enq0, enq1]


    def get_enquiry_invert_elevations(self):

        enq0 = None
        enq1 = None

        get0 = 'self.inlets[0].get_enquiry_invert_elevation()'
        get1 = 'self.inlets[1].get_enquiry_invert_elevation()'


        if self.myid == self.master_proc:

            if self.myid == self.enquiry_proc[0]:
                enq0 = eval(get0)
            else:
                enq0 = pypar.receive(self.enquiry_proc[0])


            if self.myid == self.enquiry_proc[1]:
                enq1 = eval(get1)
            else:
                enq1 = pypar.receive(self.enquiry_proc[1])

        else:
            if self.myid == self.enquiry_proc[0]:
                enq0 = eval(get0)
                pypar.send(enq0, self.master_proc)

            if self.myid == self.enquiry_proc[1]:
                enq1 = eval(get1)
                pypar.send(enq1, self.master_proc)


        return [enq0, enq1]


    def get_enquiry_velocitys(self):

        enq0 = None
        enq1 = None

        get0 = 'self.inlets[0].get_enquiry_velocity()'
        get1 = 'self.inlets[1].get_enquiry_velocity()'

        if self.myid == self.master_proc:

            if self.myid == self.enquiry_proc[0]:
                enq0 = eval(get0)
            else:
                enq0 = pypar.receive(self.enquiry_proc[0])


            if self.myid == self.enquiry_proc[1]:
                enq1 = eval(get1)
            else:
                enq1 = pypar.receive(self.enquiry_proc[1])

        else:
            if self.myid == self.enquiry_proc[0]:
                enq0 = eval(get0)
                pypar.send(enq0, self.master_proc)

            if self.myid == self.enquiry_proc[1]:
                enq1 = eval(get1)
                pypar.send(enq1, self.master_proc)


        return [enq0, enq1]


    def get_enquiry_xvelocitys(self):

        enq0 = None
        enq1 = None

        get0 = 'self.inlets[0].get_enquiry_xvelocity()'
        get1 = 'self.inlets[1].get_enquiry_xvelocity()'

        if self.myid == self.master_proc:

            if self.myid == self.enquiry_proc[0]:
                enq0 = eval(get0)
            else:
                enq0 = pypar.receive(self.enquiry_proc[0])


            if self.myid == self.enquiry_proc[1]:
                enq1 = eval(get1)
            else:
                enq1 = pypar.receive(self.enquiry_proc[1])

        else:
            if self.myid == self.enquiry_proc[0]:
                enq0 = eval(get0)
                pypar.send(enq0, self.master_proc)

            if self.myid == self.enquiry_proc[1]:
                enq1 = eval(get1)
                pypar.send(enq1, self.master_proc)


        return [enq0, enq1]

    def get_enquiry_yvelocitys(self):

        enq0 = None
        enq1 = None

        get0 = 'self.inlets[0].get_enquiry_yvelocity()'
        get1 = 'self.inlets[1].get_enquiry_yvelocity()'


        if self.myid == self.master_proc:

            if self.myid == self.enquiry_proc[0]:
                enq0 = eval(get0)
            else:
                enq0 = pypar.receive(self.enquiry_proc[0])


            if self.myid == self.enquiry_proc[1]:
                enq1 = eval(get1)
            else:
                enq1 = pypar.receive(self.enquiry_proc[1])

        else:
            if self.myid == self.enquiry_proc[0]:
                enq0 = eval(get0)
                pypar.send(enq0, self.master_proc)

            if self.myid == self.enquiry_proc[1]:
                enq1 = eval(get1)
                pypar.send(enq1, self.master_proc)


        return [enq0, enq1]


    def get_enquiry_speeds(self):

        enq0 = None
        enq1 = None

        get0 = 'self.inlets[0].get_enquiry_speed()'
        get1 = 'self.inlets[1].get_enquiry_speed()'

        
        if self.myid == self.master_proc:

            if self.myid == self.enquiry_proc[0]:
                enq0 = eval(get0)
            else:
                enq0 = pypar.receive(self.enquiry_proc[0])


            if self.myid == self.enquiry_proc[1]:
                enq1 = eval(get1)
            else:
                enq1 = pypar.receive(self.enquiry_proc[1])

        else:
            if self.myid == self.enquiry_proc[0]:
                enq0 = eval(get0)
                pypar.send(enq0, self.master_proc)

            if self.myid == self.enquiry_proc[1]:
                enq1 = eval(get1)
                pypar.send(enq1, self.master_proc)


        return [enq0, enq1]


    def get_enquiry_velocity_heads(self):

        enq0 = None
        enq1 = None

        get0 = 'self.inlets[0].get_enquiry_velocity_head()'
        get1 = 'self.inlets[1].get_enquiry_velocity_head()'

        if self.myid == self.master_proc:

            if self.myid == self.enquiry_proc[0]:
                enq0 = eval(get0)
            else:
                enq0 = pypar.receive(self.enquiry_proc[0])


            if self.myid == self.enquiry_proc[1]:
                enq1 = eval(get1)
            else:
                enq1 = pypar.receive(self.enquiry_proc[1])

        else:
            if self.myid == self.enquiry_proc[0]:
                enq0 = eval(get0)
                pypar.send(enq0, self.master_proc)

            if self.myid == self.enquiry_proc[1]:
                enq1 = eval(get1)
                pypar.send(enq1, self.master_proc)


        return [enq0, enq1]


    def get_enquiry_total_energys(self):

        enq0 = None
        enq1 = None

        get0 = 'self.inlets[0].get_enquiry_total_energy()'
        get1 = 'self.inlets[1].get_enquiry_total_energy()'

        if self.myid == self.master_proc:

            if self.myid == self.enquiry_proc[0]:
                enq0 = eval(get0)
            else:
                enq0 = pypar.receive(self.enquiry_proc[0])


            if self.myid == self.enquiry_proc[1]:
                enq1 = eval(get1)
            else:
                enq1 = pypar.receive(self.enquiry_proc[1])

        else:
            if self.myid == self.enquiry_proc[0]:
                enq0 = eval(get0)
                pypar.send(enq0, self.master_proc)

            if self.myid == self.enquiry_proc[1]:
                enq1 = eval(get1)
                pypar.send(enq1, self.master_proc)


        return [enq0, enq1]


    def get_enquiry_specific_energys(self):

        enq0 = None
        enq1 = None

        get0 = 'self.inlets[0].get_enquiry_specific_energy()'
        get1 = 'self.inlets[1].get_enquiry_specific_energy()'

        if self.myid == self.master_proc:

            if self.myid == self.enquiry_proc[0]:
                enq0 = eval(get0)
            else:
                enq0 = pypar.receive(self.enquiry_proc[0])


            if self.myid == self.enquiry_proc[1]:
                enq1 = eval(get1)
            else:
                enq1 = pypar.receive(self.enquiry_proc[1])

        else:
            if self.myid == self.enquiry_proc[0]:
                enq0 = eval(get0)
                pypar.send(enq0, self.master_proc)

            if self.myid == self.enquiry_proc[1]:
                enq1 = eval(get1)
                pypar.send(enq1, self.master_proc)


        return [enq0, enq1]


