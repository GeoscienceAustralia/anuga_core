"""
Generic implementation of update_timestep, update_ghosts and set_name for 
parallel domains (eg shallow_water or advection) 

Ole Nielsen, Stephen Roberts, Duncan Gray, Christopher Zoppou
Geoscience Australia, 2004-2010

"""

import numpy as num

import anuga_parallel.parallel_abstraction as pypar




def setup_buffers(domain):
    """Buffers for synchronisation of timesteps
    """

    domain.local_timestep = num.zeros(1, num.float)
    domain.global_timestep = num.zeros(1, num.float)

    domain.local_timesteps = num.zeros(domain.numproc, num.float)

    domain.communication_time = 0.0
    domain.communication_reduce_time = 0.0
    domain.communication_broadcast_time = 0.0


def communicate_flux_timestep(domain, yieldstep, finaltime):
    """Calculate local timestep
    """

    pypar.barrier()

    import time

    #Compute minimal timestep across all processes
    domain.local_timestep[0] = domain.flux_timestep
    t0 = time.time()

    #Determine minimum timestep across all processors
    pypar.reduce(domain.local_timestep, pypar.MIN, 0,
                      buffer=domain.global_timestep)#,
                     #bypass=True)

    domain.communication_reduce_time += time.time()-t0


    #Broadcast minimal timestep to all processors
    t0 = time.time()
    pypar.broadcast(domain.global_timestep, 0)#,
                        #bypass=True)

    domain.communication_broadcast_time += time.time()-t0

    old_timestep = domain.flux_timestep
    domain.flux_timestep = domain.global_timestep[0]
    
    
    #Compute minimal timestep on local process
    #Domain.update_timestep(domain, yieldstep, finaltime)


def communicate_ghosts(domain):

    # We must send the information from the full cells and
    # receive the information for the ghost cells
    # We have a dictionary of lists with ghosts expecting updates from
    # the separate processors

    import numpy as num
    import time
    t0 = time.time()

    # update of non-local ghost cells
    for iproc in range(domain.numproc):
        if iproc == domain.processor:
            #Send data from iproc processor to other processors
            for send_proc in domain.full_send_dict:
                if send_proc != iproc:

                    Idf  = domain.full_send_dict[send_proc][0]
                    Xout = domain.full_send_dict[send_proc][2]

                    for i, q in enumerate(domain.conserved_quantities):
                        #print 'Send',i,q
                        Q_cv =  domain.quantities[q].centroid_values
                        Xout[:,i] = num.take(Q_cv, Idf)

                    pypar.send(Xout, int(send_proc), use_buffer=True)


        else:
            #Receive data from the iproc processor
            if  domain.ghost_recv_dict.has_key(iproc):

                Idg = domain.ghost_recv_dict[iproc][0]
                X   = domain.ghost_recv_dict[iproc][2]

                X = pypar.receive(int(iproc), buffer=X)

                for i, q in enumerate(domain.conserved_quantities):
                    #print 'Receive',i,q
                    Q_cv =  domain.quantities[q].centroid_values
                    num.put(Q_cv, Idg, X[:,i])

    #local update of ghost cells
    iproc = domain.processor
    if domain.full_send_dict.has_key(iproc):

        # LINDA:
        # now store full as local id, global id, value
        Idf  = domain.full_send_dict[iproc][0]

        # LINDA:
        # now store ghost as local id, global id, value
        Idg = domain.ghost_recv_dict[iproc][0]

        for i, q in enumerate(domain.conserved_quantities):
            #print 'LOCAL SEND RECEIVE',i,q
            Q_cv =  domain.quantities[q].centroid_values
            num.put(Q_cv, Idg, num.take(Q_cv, Idf))

    domain.communication_time += time.time()-t0

