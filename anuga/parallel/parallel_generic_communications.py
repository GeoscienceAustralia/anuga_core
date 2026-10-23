"""
Generic implementation of update_timestep and update_ghosts for
parallel domains (eg shallow_water or advection)

Ole Nielsen, Stephen Roberts, Duncan Gray, Christopher Zoppou
Geoscience Australia, 2004-2010

"""
import numpy as num

import anuga.utilities.parallel_abstraction as pypar




def setup_buffers(domain):
    """Buffers for synchronisation of timesteps
    """

    domain.local_timestep = num.zeros(1, float)
    domain.global_timestep = num.zeros(1, float)

    domain.local_timesteps = num.zeros(domain.numproc, float)

    domain.communication_time = 0.0
    domain.communication_reduce_time = 0.0
    domain.communication_broadcast_time = 0.0

    domain.calls_to_update_ghosts = 0
    domain.calls_to_update_timestep = 0


def communicate_flux_timestep(domain, yieldstep, finaltime):
    """Calculate local timestep
    """

    import time
    import anuga

    if anuga.myid == 0:
        #print('o', end = '')
        domain.calls_to_update_timestep += 1

    # disable allreduce if fixed_flux_timestep is set
    if domain.fixed_flux_timestep is not None:
        domain.flux_timestep = domain.fixed_flux_timestep
        if not domain.test_allreduce:
            return


    #Compute minimal timestep across all processes
    domain.local_timestep[0] = domain.flux_timestep
    t0 = time.time()


    local_timestep = domain.local_timestep
    global_timestep = domain.global_timestep

    #pypar.allreduce(domain.local_timestep, pypar.MIN,
    #                  buffer=domain.global_timestep,
    #                  bypass=True)

    from mpi4py.MPI import MIN
    pypar.comm.Allreduce(local_timestep, global_timestep, op=MIN)

    domain.communication_reduce_time += time.time()-t0

    t0 = time.time()

    domain.communication_broadcast_time += time.time()-t0

    domain.flux_timestep = domain.global_timestep[0]


def communicate_ghosts_blocking(domain, quantities=None):

    # We must send the information from the full cells and
    # receive the information for the ghost cells
    # We have a dictionary of lists with ghosts expecting updates from
    # the separate processors

    import numpy as num
    import time
    t0 = time.time()

    if quantities is None:
        quantities = domain.conserved_quantities

    # update of non-local ghost cells
    for iproc in range(domain.numproc):
        if iproc == domain.processor:
            #Send data from iproc processor to other processors
            for send_proc in domain.full_send_dict:
                if send_proc != iproc:

                    Idf  = domain.full_send_dict[send_proc][0]
                    Xout = domain.full_send_dict[send_proc][2]

                    for i, q in enumerate(quantities):
                        #print 'Send',i,q
                        Q_cv =  domain.quantities[q].centroid_values
                        Xout[:,i] = num.take(Q_cv, Idf)

                    pypar.send(Xout, int(send_proc), use_buffer=True, bypass=True)


        else:
            #Receive data from the iproc processor
            if  iproc in domain.ghost_recv_dict:

                Idg = domain.ghost_recv_dict[iproc][0]
                X   = domain.ghost_recv_dict[iproc][2]

                X = pypar.receive(int(iproc), buffer=X, bypass=True)

                for i, q in enumerate(quantities):
                    #print 'Receive',i,q
                    Q_cv =  domain.quantities[q].centroid_values
                    num.put(Q_cv, Idg, X[:,i])

    #local update of ghost cells
    iproc = domain.processor
    if iproc in domain.full_send_dict:

        # LINDA:
        # now store full as local id, global id, value
        Idf  = domain.full_send_dict[iproc][0]

        # LINDA:
        # now store ghost as local id, global id, value
        Idg = domain.ghost_recv_dict[iproc][0]

        for i, q in enumerate(quantities):
            #print 'LOCAL SEND RECEIVE',i,q
            Q_cv =  domain.quantities[q].centroid_values
            num.put(Q_cv, Idg, num.take(Q_cv, Idf))

    domain.communication_time += time.time()-t0



def communicate_ghosts_non_blocking(domain, quantities=None):

    # We must send the information from the full cells and
    # receive the information for the ghost cells
    # We have a dictionary of lists with ghosts expecting updates from
    # the separate processors
    # Using isend and irecv

    import numpy as num
    import time
    import anuga
    t0 = time.time()

    if anuga.myid == 0:
        #print('.', end = '')
        domain.calls_to_update_ghosts += 1

    sendDict = domain.full_send_dict
    recvDict = domain.ghost_recv_dict

    if quantities is None:
        quantities = domain.conserved_quantities

    # update of non-local ghost cells by copying full cell data into the
    # Xout buffer arrays

    #iproc == domain.processor

    #Setup send buffer arrays for sending full data to other processors
    for send_proc in domain.full_send_dict:
        Idf  = sendDict[send_proc][0]
        Xout = sendDict[send_proc][2]

        for i, q in enumerate(quantities):
            #print 'Store send data',i,q
            Q_cv =  domain.quantities[q].centroid_values
            Xout[:,i] = num.take(Q_cv, Idf)

    #--------------------------------------------
    # Do all the comuunication using isend/irecv
    # via the buffers in the
    # full_send_dict and ghost_recv_dict
    #--------------------------------------------


    #-------------------------
    # Do the Irecvs first
    #-------------------------
    recv_requests = []
    for recv_proc in recvDict:

        Idg = recvDict[recv_proc][0]
        X   = recvDict[recv_proc][2]

        request = pypar.comm.Irecv(X, recv_proc, 123)
        recv_requests.append(request)

    #---------------------
    # Do the Isends second
    #---------------------
    send_requests = []
    for send_proc in sendDict:

        Idg = sendDict[send_proc][0]
        X   = sendDict[send_proc][2]

        request = pypar.comm.Isend(X, send_proc, 123)
        send_requests.append(request)

    #-----------------------------------------
    # Now complete communication.
    # We could put some computation between the
    # communication calls above and this call.
    # Question: Do we need to wait for the sends to complete as well?
    # Answer: Yes, we should wait for the sends to complete as well, otherwise
    # we might be overwriting the send buffers before the data has been sent.
    #-----------------------------------------
    import mpi4py
    mpi4py.MPI.Request.Waitall(recv_requests + send_requests)

    # Now copy data from receive buffers to the domain
    for recv_proc in recvDict:
        Idg  = recvDict[recv_proc][0]
        X    = recvDict[recv_proc][2]

        for i, q in enumerate(quantities):
            #print 'Read receive data',i,q
            Q_cv =  domain.quantities[q].centroid_values
            num.put(Q_cv, Idg, X[:,i])


    domain.communication_time += time.time()-t0

    # Tracers are not Quantity objects, so the loop above cannot see them.
    communicate_tracer_ghosts(domain)



def communicate_tracer_ghosts(domain):
    """Exchange generic tracer values for ghost cells (legacy / mode 1).

    The hydro exchange above iterates ``domain.conserved_quantities`` and reads
    ``domain.quantities[q].centroid_values``. Tracers are raw ``(ns, n)`` arrays
    on the Domain rather than Quantity objects, so they are invisible to it and
    a partitioned run silently transports them with stale ghost values. The
    symptom is subtle: total water volume stays exactly conserved while a
    uniform tracer -- which must satisfy m == h exactly -- drifts away from it.

    This cannot reuse the hydro buffers: ``full_send_dict[proc][2]`` is
    allocated at distribute time with width ``len(conserved_quantities)``. So we
    keep our own, allocated lazily and cached on the domain, and use a distinct
    MPI tag.

    Only the CONSERVED m = h*c is exchanged, matching both the hydro path and
    gpu_halo.c: c is re-derived from m for ghost cells in extrapolation Step 1,
    exactly as height is re-derived from stage.
    """
    import numpy as num
    import time
    import mpi4py

    ns = getattr(domain, 'number_of_tracers', 0)
    if ns == 0:
        return

    t0 = time.time()

    sendDict = domain.full_send_dict
    recvDict = domain.ghost_recv_dict
    m = domain.tracer_conserved_values

    # Cache the buffers, but key the cache on ns: add_tracer() reallocates the
    # tracer arrays at a new width, and a stale buffer would silently exchange
    # the wrong number of slots.
    cache = getattr(domain, '_tracer_halo_buffers', None)
    if cache is None or cache.get('ns') != ns:
        cache = {'ns': ns, 'send': {}, 'recv': {}}
        for proc in sendDict:
            cache['send'][proc] = num.zeros((len(sendDict[proc][0]), ns),
                                            dtype=num.float64)
        for proc in recvDict:
            cache['recv'][proc] = num.zeros((len(recvDict[proc][0]), ns),
                                            dtype=num.float64)
        domain._tracer_halo_buffers = cache

    # Pack owned cells
    for send_proc in sendDict:
        Idf = sendDict[send_proc][0]
        Xout = cache['send'][send_proc]
        for s in range(ns):
            Xout[:, s] = num.take(m[s], Idf)

    # Tag 124: the hydro exchange above uses 123 on the same communicator and
    # the two must not be matched against each other.
    recv_requests = []
    for recv_proc in recvDict:
        recv_requests.append(
            pypar.comm.Irecv(cache['recv'][recv_proc], recv_proc, 124))

    send_requests = []
    for send_proc in sendDict:
        send_requests.append(
            pypar.comm.Isend(cache['send'][send_proc], send_proc, 124))

    mpi4py.MPI.Request.Waitall(recv_requests + send_requests)

    # Unpack into ghost cells
    for recv_proc in recvDict:
        Idg = recvDict[recv_proc][0]
        X = cache['recv'][recv_proc]
        for s in range(ns):
            num.put(m[s], Idg, X[:, s])

    domain.communication_time += time.time() - t0


def communicate_ghosts_asynchronous(domain, quantities=None):

    # We must send the information from the full cells and
    # receive the information for the ghost cells
    # We have a dictionary of lists with ghosts expecting updates from
    # the separate processors
    # Using isend and irecv

    import numpy as num
    import time
    t0 = time.time()

    if quantities is None:
        quantities = domain.conserved_quantities

    # update of non-local ghost cells by copying full cell data into the
    # Xout buffer arrays

    #iproc == domain.processor

    #Setup send buffer arrays for sending full data to other processors
    for send_proc in domain.full_send_dict:
        Idf  = domain.full_send_dict[send_proc][0]
        Xout = domain.full_send_dict[send_proc][2]

        for i, q in enumerate(quantities):
            #print 'Store send data',i,q
            Q_cv =  domain.quantities[q].centroid_values
            Xout[:,i] = num.take(Q_cv, Idf)

    # Do all the comuunication using isend/irecv via the buffers in the
    # full_send_dict and ghost_recv_dict

    pypar.send_recv_via_dicts(domain.full_send_dict,domain.ghost_recv_dict)

    # Now copy data from receive buffers to the domain
    for recv_proc in domain.ghost_recv_dict:
        Idg  = domain.ghost_recv_dict[recv_proc][0]
        X    = domain.ghost_recv_dict[recv_proc][2]

        for i, q in enumerate(quantities):
            #print 'Read receive data',i,q
            Q_cv =  domain.quantities[q].centroid_values
            num.put(Q_cv, Idg, X[:,i])


    domain.communication_time += time.time()-t0

