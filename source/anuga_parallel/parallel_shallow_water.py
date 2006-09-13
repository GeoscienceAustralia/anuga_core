"""Class Parallel_Shallow_Water_Domain -
2D triangular domains for finite-volume computations of
the shallow water equation, with extra structures to allow
communication between other Parallel_Domains and itself

This module contains a specialisation of class Domain
from module shallow_water.py

Ole Nielsen, Stephen Roberts, Duncan Gray, Christopher Zoppou
Geoscience Australia, 2004-2005

"""

import logging, logging.config
logger = logging.getLogger('parallel')
logger.setLevel(logging.WARNING)

try:
    logging.config.fileConfig('log.ini')
except:
    pass


from anuga.shallow_water.shallow_water_domain import *
from Numeric import zeros, Float, Int, ones, allclose, array

import pypar


class Parallel_Domain(Domain):

    def __init__(self, coordinates, vertices, boundary = None,
                 full_send_dict = None, ghost_recv_dict = None):

        Domain.__init__(self,
                        coordinates,
                        vertices,
                        boundary,
                        full_send_dict=full_send_dict,
                        ghost_recv_dict=ghost_recv_dict,
                        processor=pypar.rank(),
                        numproc=pypar.size())

        N = self.number_of_elements

#        self.processor = pypar.rank()
#        self.numproc   = pypar.size()
#
#        # Setup Communication Buffers
#        self.nsys = 3
#        for key in full_send_dict:
#            buffer_shape = full_send_dict[key][0].shape[0]
#            full_send_dict[key].append(zeros( (buffer_shape,self.nsys) ,Float))
#
#
#        for key in ghost_recv_dict:
#            buffer_shape = ghost_recv_dict[key][0].shape[0]
#            ghost_recv_dict[key].append(zeros( (buffer_shape,self.nsys) ,Float))
#
#        self.full_send_dict  = full_send_dict
        self.ghost_recv_dict = ghost_recv_dict

        # Buffers for synchronisation of timesteps
        self.local_timestep = zeros(1, Float)
        self.global_timestep = zeros(1, Float)

        self.local_timesteps = zeros(self.numproc, Float)


        self.communication_time = 0.0
        self.communication_reduce_time = 0.0
        self.communication_broadcast_time = 0.0

        


    def set_name(self, name):
        """Assign name based on processor number 
        """

        # Call parents method with processor number attached.
        Domain.set_name(self, name + '_P%d_%d' %(self.processor, self.numproc))


    def check_integrity(self):
        Domain.check_integrity(self)

        msg = 'Will need to check global and local numbering'
        assert self.conserved_quantities[0] == 'stage', msg
        assert self.conserved_quantities[1] == 'xmomentum', msg
        assert self.conserved_quantities[2] == 'ymomentum', msg


    def update_timestep_1(self, yieldstep, finaltime):
        """Calculate local timestep using broadcasts
        """

        #LINDA:
        # Moved below so timestep is found before doing update
        
        #Domain.update_timestep(self, yieldstep, finaltime)

        import time


        t0 = time.time()

        #Broadcast local timestep from every processor to every other
        for pid in range(self.numproc):
            #print 'P%d calling broadcast from %d' %(self.processor, pid)
            self.local_timestep[0] = self.timestep
            pypar.broadcast(self.local_timestep, pid, bypass=True)
            self.local_timesteps[pid] = self.local_timestep[0]

        self.timestep = min(self.local_timesteps)

        pypar.barrier()
        self.communication_broadcast_time += time.time()-t0

        # LINDA:
        # Moved timestep to here
        
        Domain.update_timestep(self, yieldstep, finaltime)


    def update_timestep(self, yieldstep, finaltime):
        """Calculate local timestep
        """

        # LINDA: Moved below so timestep is updated before
        # calculating statistic
        
        #Compute minimal timestep on local process
        #Domain.update_timestep(self, yieldstep, finaltime)

        pypar.barrier()

        import time
        #Compute minimal timestep across all processes
        self.local_timestep[0] = self.timestep
        use_reduce_broadcast = True
        if use_reduce_broadcast:
            t0 = time.time()
            pypar.reduce(self.local_timestep, pypar.MIN, 0,
                         buffer=self.global_timestep,
                         bypass=True)

        else:
            #Alternative: Try using straight send and receives
            t0 = time.time()
            self.global_timestep[0] = self.timestep

            if self.processor == 0:
                for i in range(1, self.numproc):
                    pypar.receive(i,
                                  buffer=self.local_timestep,
                                  bypass=True)

                    if self.local_timestep[0] < self.global_timestep[0]:
                        self.global_timestep[0] = self.local_timestep[0]
            else:
                pypar.send(self.local_timestep, 0,
                           use_buffer=True, bypass=True)


        self.communication_reduce_time += time.time()-t0


        #Broadcast minimal timestep to all
        t0 = time.time()
        pypar.broadcast(self.global_timestep, 0,
                        bypass=True)

        self.communication_broadcast_time += time.time()-t0


        self.timestep = self.global_timestep[0]
        
        # LINDA:
        # update local stats now
        
        #Compute minimal timestep on local process
        Domain.update_timestep(self, yieldstep, finaltime)

        # FIXME (Ole) We should update the variable min_timestep for use
        # with write_time (or redo write_time) 

    #update_timestep = update_timestep_1

    def update_ghosts(self):

        # We must send the information from the full cells and
        # receive the information for the ghost cells
        # We have a dictionary of lists with ghosts expecting updates from
        # the separate processors


        from Numeric import take,put
        import time
        t0 = time.time()

        # update of non-local ghost cells
        for iproc in range(self.numproc):
            if iproc == self.processor:
                #Send data from iproc processor to other processors
                for send_proc in self.full_send_dict:
                    if send_proc != iproc:

                        Idf  = self.full_send_dict[send_proc][0]
                        Xout = self.full_send_dict[send_proc][2]

                        for i, q in enumerate(self.conserved_quantities):
                            #print 'Send',i,q
                            Q_cv =  self.quantities[q].centroid_values
                            Xout[:,i] = take(Q_cv,     Idf)

                        pypar.send(Xout, send_proc,
                                   use_buffer=True, bypass = True)


            else:
                #Receive data from the iproc processor
                if  self.ghost_recv_dict.has_key(iproc):

                    Idg = self.ghost_recv_dict[iproc][0]
                    X = self.ghost_recv_dict[iproc][2]

                    X = pypar.receive(iproc, buffer=X, bypass = True)

                    for i, q in enumerate(self.conserved_quantities):
                        #print 'Receive',i,q
                        Q_cv =  self.quantities[q].centroid_values
                        put(Q_cv,     Idg, X[:,i])

        #local update of ghost cells
        iproc = self.processor
        if self.full_send_dict.has_key(iproc):

            # LINDA:
            # now store full as local id, global id, value
            Idf  = self.full_send_dict[iproc][0]

            # LINDA:
            # now store ghost as local id, global id, value
            Idg = self.ghost_recv_dict[iproc][0]

            for i, q in enumerate(self.conserved_quantities):
                #print 'LOCAL SEND RECEIVE',i,q
                Q_cv =  self.quantities[q].centroid_values
                put(Q_cv,     Idg, take(Q_cv,     Idf))

        self.communication_time += time.time()-t0


    def write_time(self):
        if self.min_timestep == self.max_timestep:
            print 'Processor %d, Time = %.4f, delta t = %.8f, steps=%d (%d)'\
                  %(self.processor, self.time, self.min_timestep, self.number_of_steps,
                    self.number_of_first_order_steps)
        elif self.min_timestep > self.max_timestep:
            print 'Processor %d, Time = %.4f, steps=%d (%d)'\
                  %(self.processor, self.time, self.number_of_steps,
                    self.number_of_first_order_steps)
        else:
            print 'Processor %d, Time = %.4f, delta t in [%.8f, %.8f], steps=%d (%d)'\
                  %(self.processor, self.time, self.min_timestep,
                    self.max_timestep, self.number_of_steps,
                    self.number_of_first_order_steps)


    def evolve(self, yieldstep = None, finaltime = None):
        """Specialisation of basic evolve method from parent class
        """

        #Initialise real time viz if requested
        if self.time == 0.0:
            pass

        #Call basic machinery from parent class
        for t in Domain.evolve(self, yieldstep, finaltime):

            #Pass control on to outer loop for more specific actions
            yield(t)
