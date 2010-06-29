"""Class Parallel_shallow_water_domain -
2D triangular domains for finite-volume computations of
the shallow water equation, with extra structures to allow
communication between other Parallel_domains and itself

This module contains a specialisation of class Domain
from module shallow_water.py

Ole Nielsen, Stephen Roberts, Duncan Gray, Christopher Zoppou
Geoscience Australia, 2004-2005

"""

from anuga import Domain


import numpy as num

import pypar


class Parallel_domain(Domain):

    def __init__(self, coordinates, vertices,
                 boundary=None,
                 full_send_dict=None,
                 ghost_recv_dict=None,
                 number_of_full_nodes=None,
                 number_of_full_triangles=None):

        Domain.__init__(self,
                        coordinates,
                        vertices,
                        boundary,
                        full_send_dict=full_send_dict,
                        ghost_recv_dict=ghost_recv_dict,
                        processor=pypar.rank(),
                        numproc=pypar.size(),
                        number_of_full_nodes=number_of_full_nodes,
                        number_of_full_triangles=number_of_full_triangles)

        N = len(self) # number_of_triangles


        # Buffers for synchronisation of timesteps
        self.local_timestep = num.zeros(1, num.float)
        self.global_timestep = num.zeros(1, num.float)

        self.local_timesteps = num.zeros(self.numproc, num.float)


        self.communication_time = 0.0
        self.communication_reduce_time = 0.0
        self.communication_broadcast_time = 0.0

        


    def set_name(self, name):
        """Assign name based on processor number 
        """

        if name.endswith('.sww'):
            name = name[:-4]

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
            self.local_timestep[0] = self.flux_timestep
            pypar.broadcast(self.local_timestep, pid, bypass=True)
            self.local_timesteps[pid] = self.local_timestep[0]

        self.flux_timestep = min(self.local_timesteps)

        #print 'Flux Timestep %d P%d_%d' %(self.flux_timestep, self.processor, self.numproc)

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
        self.local_timestep[0] = self.flux_timestep
        use_reduce_broadcast = True
        if use_reduce_broadcast:
            t0 = time.time()
            pypar.reduce(self.local_timestep, pypar.MIN, 0,
                         buffer=self.global_timestep)#,
                         #bypass=True)

        else:
            #Alternative: Try using straight send and receives
            t0 = time.time()
            self.global_timestep[0] = self.flux_timestep

            if self.processor == 0:
                for i in range(1, self.numproc):
                    pypar.receive(i,
                                  buffer=self.local_timestep)

                    if self.local_timestep[0] < self.global_timestep[0]:
                        self.global_timestep[0] = self.local_timestep[0]
            else:
                pypar.send(self.local_timestep, 0,
                           use_buffer=True)


        self.communication_reduce_time += time.time()-t0


        #Broadcast minimal timestep to all
        t0 = time.time()
        pypar.broadcast(self.global_timestep, 0)#,
                        #bypass=True)

        self.communication_broadcast_time += time.time()-t0

        old_timestep = self.flux_timestep
        self.flux_timestep = self.global_timestep[0]
        #print 'Flux Timestep %15.5e %15.5e P%d_%d' %(self.flux_timestep, old_timestep, self.processor, self.numproc)
        
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

        import numpy as num
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
                            Xout[:,i] = num.take(Q_cv, Idf)

                        pypar.send(Xout, int(send_proc), use_buffer=True)


            else:
                #Receive data from the iproc processor
                if  self.ghost_recv_dict.has_key(iproc):

                    Idg = self.ghost_recv_dict[iproc][0]
                    X   = self.ghost_recv_dict[iproc][2]

                    X = pypar.receive(int(iproc), buffer=X)

                    for i, q in enumerate(self.conserved_quantities):
                        #print 'Receive',i,q
                        Q_cv =  self.quantities[q].centroid_values
                        num.put(Q_cv, Idg, X[:,i])

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
                num.put(Q_cv, Idg, num.take(Q_cv, Idf))

        self.communication_time += time.time()-t0

