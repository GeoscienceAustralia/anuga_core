
from builtins import range
import sys


"""Class Parallel_domain -
2D triangular domains for finite-volume computations of
the advection equation, with extra structures to allow
communication between other Parallel_domains and itself

This module contains a specialisation of class Domain from module advection.py

Ole Nielsen, Stephen Roberts, Duncan Gray, Christopher Zoppou
Geoscience Australia, 2004-2005
"""

from anuga.advection import *

from anuga import Domain

import numpy as num

from anuga.utilities import parallel_abstraction as pypar


class Parallel_domain(Domain):

    def __init__(self,
                 coordinates,
                 vertices,
                 boundary = None,
                 full_send_dict = None,
                 ghost_recv_dict = None,
                 velocity = None):

        Domain.__init__(self,
                        coordinates,
                        vertices,
                        boundary,
                        velocity = velocity,
                        full_send_dict=full_send_dict,
                        ghost_recv_dict=ghost_recv_dict,
                        processor=pypar.rank(),
                        numproc=pypar.size()
                        )

        N = self.number_of_elements


        self.communication_time = 0.0
        self.communication_reduce_time = 0.0


        print('processor',self.processor)
        print('numproc',self.numproc)

    def check_integrity(self):
        Domain.check_integrity(self)

        msg = 'Will need to check global and local numbering'
        assert self.conserved_quantities[0] == 'stage', msg

    def update_timestep(self, yieldstep, finaltime):

        #LINDA:
        # moved the calculation so that it is done after timestep
        # has been broadcast
        
#        # Calculate local timestep
#        Domain.update_timestep(self, yieldstep, finaltime)

        import time
        t0 = time.time()

        # For some reason it looks like pypar only reduces numeric arrays
        # hence we need to create some dummy arrays for communication
        ltimestep = num.ones( 1, float )
        ltimestep[0] = self.flux_timestep
        gtimestep = num.zeros( 1, float ) # Buffer for results

        #ltimestep = self.flux_timeste

        #print self.processor, ltimestep, gtimestep
        
        gtimestep = pypar.reduce(ltimestep, pypar.MIN, 0, buffer=gtimestep)

        #print self.processor, ltimestep, gtimestep
        
        pypar.broadcast(gtimestep,0)

        #print self.processor, ltimestep, gtimestep

        self.flux_timestep = gtimestep[0]
        
        self.communication_reduce_time += time.time()-t0

        # LINDA:
        # Now update time stats
        
        # Calculate local timestep
        Domain.update_timestep(self, yieldstep, finaltime)

    def update_ghosts(self):

        # We must send the information from the full cells and
        # receive the information for the ghost cells
        # We have a dictionary of lists with ghosts expecting updates from
        # the separate processors

        #from Numeric import take,put
        import numpy as num
        import time
        t0 = time.time()

        stage_cv = self.quantities['stage'].centroid_values

        # update of non-local ghost cells
        for iproc in range(self.numproc):
            if iproc == self.processor:
                #Send data from iproc processor to other processors
                for send_proc in self.full_send_dict:
                    if send_proc != iproc:

                        Idf  = self.full_send_dict[send_proc][0]
                        Xout = self.full_send_dict[send_proc][2]

                        N = len(Idf)

                        #for i in range(N):
                        #    Xout[i,0] = stage_cv[Idf[i]]
                        Xout[:,0] = num.take(stage_cv, Idf)

                        pypar.send(Xout,send_proc)


            else:
                #Receive data from the iproc processor
                if  iproc in self.ghost_recv_dict:

                    # LINDA:
                    # now store ghost as local id, global id, value
                    Idg = self.ghost_recv_dict[iproc][0]
                    X = self.ghost_recv_dict[iproc][2]

                    X = pypar.receive(iproc,X)
                    N = len(Idg)

                    num.put(stage_cv, Idg, X[:,0])
                    #for i in range(N):
                    #    stage_cv[Idg[i]] = X[i,0]


        #local update of ghost cells
        iproc = self.processor
        if iproc in self.full_send_dict:

            # LINDA:
            # now store full as local id, global id, value
            Idf  = self.full_send_dict[iproc][0]

            # LINDA:
            # now store ghost as local id, global id, value
            Idg = self.ghost_recv_dict[iproc][0]

            N = len(Idg)

            #for i in range(N):
            #    #print i,Idg[i],Idf[i]
            #    stage_cv[Idg[i]] = stage_cv[Idf[i]]

            num.put(stage_cv, Idg, num.take(stage_cv, Idf))


        self.communication_time += time.time()-t0


    ## def write_time(self):
    ##     if self.min_timestep == self.max_timestep:
    ##         print 'Processor %d, Time = %.4f, delta t = %.8f, steps=%d (%d)'\
    ##               %(self.processor, self.time, self.min_timestep, self.number_of_steps,
    ##                 self.number_of_first_order_steps)
    ##     elif self.min_timestep > self.max_timestep:
    ##         print 'Processor %d, Time = %.4f, steps=%d (%d)'\
    ##               %(self.processor, self.time, self.number_of_steps,
    ##                 self.number_of_first_order_steps)
    ##     else:
    ##         print 'Processor %d, Time = %.4f, delta t in [%.8f, %.8f], steps=%d (%d)'\
    ##               %(self.processor, self.time, self.min_timestep,
    ##                 self.max_timestep, self.number_of_steps,
    ##                 self.number_of_first_order_steps)



    ## def evolve(self, yieldstep = None, finaltime = None):
    ##     """Specialisation of basic evolve method from parent class
    ##     """

    ##     #Initialise real time viz if requested
    ##     if self.time == 0.0:
    ##         pass

    ##     #Call basic machinery from parent class
    ##     for t in Domain.evolve(self, yieldstep, finaltime):

    ##         #Pass control on to outer loop for more specific actions
    ##         yield(t)

