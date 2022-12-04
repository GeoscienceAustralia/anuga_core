"""Trying to lump parallel stuff into simpler interface


"""



from builtins import range
import numpy as num

# The abstract Python-MPI interface
from anuga.utilities.parallel_abstraction import size, rank, get_processor_name
from anuga.utilities.parallel_abstraction import finalize, send, receive, reduce
from anuga.utilities.parallel_abstraction import pypar_available, barrier

from anuga.parallel.sequential_distribute import sequential_distribute_dump
from anuga.parallel.sequential_distribute import sequential_distribute_load

# ANUGA parallel engine (only load if pypar can)
if pypar_available:
    from anuga.parallel.distribute_mesh  import send_submesh
    from anuga.parallel.distribute_mesh  import rec_submesh
    from anuga.parallel.distribute_mesh  import extract_submesh

    # Mesh partitioning using Metis
    from anuga.parallel.distribute_mesh import build_submesh
    from anuga.parallel.distribute_mesh import pmesh_divide_metis_with_map

    from anuga.parallel.parallel_shallow_water import Parallel_domain



from anuga.abstract_2d_finite_volumes.neighbour_mesh import Mesh

#------------------------------------------------------------------------------
# Read in processor information
#------------------------------------------------------------------------------

numprocs = size()
myid = rank()
processor_name = get_processor_name()
#print 'I am processor %d of %d on node %s' %(myid, numprocs, processor_name)



def collect_value(value):

    value = value

    if myid == 0:
        for i in range(numprocs):
            if i == 0: continue
            val = receive(i)
            value = value + val
    else:
        send(value, 0)


    if myid == 0:
        for i in range(1,numprocs):
            send(value,i)
    else:
        value = receive(0)


    return value




def distribute(domain, verbose=False, debug=False, parameters = None):
    """ Distribute the domain to all processes

    parameters allows user to change size of ghost layer
    """

    if not pypar_available or numprocs == 1 : return domain # Bypass


    if myid == 0:
        from .sequential_distribute import Sequential_distribute
        partition = Sequential_distribute(domain, verbose, debug, parameters)

        partition.distribute(numprocs)

        kwargs, points, vertices, boundary, quantities, boundary_map, \
                domain_name, domain_dir, domain_store, domain_store_centroids, \
                domain_minimum_storable_height, domain_minimum_allowed_height, \
                domain_flow_algorithm, domain_georef, \
                domain_quantities_to_be_stored, domain_smooth, domain_low_froude \
                 = partition.extract_submesh(0)

        for p in range(1, numprocs):

            tostore = partition.extract_submesh(p)

            send(tostore,p)

    else:

        kwargs, points, vertices, boundary, quantities, boundary_map, \
            domain_name, domain_dir, domain_store, domain_store_centroids, \
            domain_minimum_storable_height, domain_minimum_allowed_height, \
            domain_flow_algorithm, domain_georef, \
            domain_quantities_to_be_stored, domain_smooth, domain_low_froude\
             = receive(0)

    #---------------------------------------------------------------------------
    # Now Create parallel domain
    #---------------------------------------------------------------------------
    parallel_domain = Parallel_domain(points, vertices, boundary, **kwargs)


    #------------------------------------------------------------------------
    # Copy in quantity data
    #------------------------------------------------------------------------
    for q in quantities:
        try:
            parallel_domain.set_quantity(q, quantities[q])
        except KeyError:
            #print 'Try to create quantity %s'% q
            from anuga import Quantity
            Q = Quantity(parallel_domain, name=q, register=True)
            parallel_domain.set_quantity(q, quantities[q])

    #------------------------------------------------------------------------
    # Transfer boundary conditions to each subdomain
    #------------------------------------------------------------------------
    boundary_map['ghost'] = None  # Add binding to ghost boundary
    parallel_domain.set_boundary(boundary_map)


    #------------------------------------------------------------------------
    # Transfer other attributes to each subdomain
    #------------------------------------------------------------------------

    parallel_domain.set_flow_algorithm(domain_flow_algorithm)
    parallel_domain.set_name(domain_name)
    parallel_domain.set_datadir(domain_dir)
    parallel_domain.set_store(domain_store)
    parallel_domain.set_low_froude(domain_low_froude)
    parallel_domain.set_store_centroids(domain_store_centroids)
    parallel_domain.set_minimum_storable_height(domain_minimum_storable_height)
    parallel_domain.set_minimum_allowed_height(domain_minimum_allowed_height)
    parallel_domain.geo_reference = domain_georef
    parallel_domain.set_quantities_to_be_stored(domain_quantities_to_be_stored)
    parallel_domain.smooth = domain_smooth

    return parallel_domain







def old_distribute(domain, verbose=False, debug=False, parameters = None):
    """ Distribute the domain to all processes

    parameters values
    """


    if debug:
        verbose = True

    barrier()

    # FIXME: Dummy assignment (until boundaries are refactored to
    # be independent of domains until they are applied)
    if myid == 0:
        bdmap = {}
        for tag in domain.get_boundary_tags():
            bdmap[tag] = None


        domain.set_boundary(bdmap)


    if not pypar_available or numprocs == 1 : return domain # Bypass

    # For some obscure reason this communication must happen prior to
    # the more complex mesh distribution - Oh Well!
    if myid == 0:
        domain_name = domain.get_name()
        domain_dir = domain.get_datadir()
        domain_store = domain.get_store()
        domain_store_centroids = domain.get_store_centroids()
        domain_smooth = domain.smooth
        domain_reduction = domain.reduction
        domain_minimum_storable_height = domain.minimum_storable_height
        domain_flow_algorithm = domain.get_flow_algorithm()
        domain_minimum_allowed_height = domain.get_minimum_allowed_height()
        georef = domain.geo_reference
        number_of_global_triangles = domain.number_of_triangles
        number_of_global_nodes = domain.number_of_nodes

        # FIXME - what other attributes need to be transferred?

        for p in range(1, numprocs):
            # FIXME SR: Creates cPickle dump
            send((domain_name, domain_dir, domain_store, \
                  domain_store_centroids, domain_smooth, domain_reduction, \
                  domain_minimum_storable_height, domain_flow_algorithm, \
                  domain_minimum_allowed_height, georef, \
                  number_of_global_triangles, number_of_global_nodes), p)
    else:
        if verbose: print('P%d: Receiving domain attributes' %(myid))

        domain_name, domain_dir, domain_store, \
                  domain_store_centroids, domain_smooth, domain_reduction, \
                  domain_minimum_storable_height, domain_flow_algorithm, \
                  domain_minimum_allowed_height, georef, \
                  number_of_global_triangles, \
                  number_of_global_nodes = receive(0)



    # Distribute boundary conditions
    # FIXME: This cannot handle e.g. Time_boundaries due to
    # difficulties pickling functions
    if myid == 0:
        boundary_map = domain.boundary_map
        for p in range(1, numprocs):
            # FIXME SR: Creates cPickle dump
            send(boundary_map, p)
    else:
        if verbose: print('P%d: Receiving boundary map' %(myid))

        boundary_map = receive(0)


    send_s2p = False

    if myid == 0:
        # Partition and distribute mesh.
        # Structures returned is in the
        # correct form for the ANUGA data structure


        points, vertices, boundary, quantities,\
                ghost_recv_dict, full_send_dict,\
                number_of_full_nodes, number_of_full_triangles,\
                s2p_map, p2s_map, tri_map, node_map, tri_l2g, node_l2g, \
                ghost_layer_width =\
                distribute_mesh(domain, verbose=verbose, debug=debug, parameters=parameters)

        # Extract l2g maps
        #tri_l2g  = extract_l2g_map(tri_map)
        #node_l2g = extract_l2g_map(node_map)
        #tri_l2g = p2s_map[tri_l2g]

        if debug:
            print('P%d' %myid)
            print('tri_map ',tri_map)
            print('node_map',node_map)
            print('tri_l2g', tri_l2g)
            print('node_l2g', node_l2g)
            print('s2p_map', s2p_map)
            print('p2s_map', p2s_map)


        def protocol(x):
            vanilla=False
            from anuga.utilities import parallel_abstraction as pypar
            control_info, x = pypar.create_control_info(x, vanilla, return_object=True)
            print('protocol', control_info[0])

        # Send serial to parallel (s2p) and parallel to serial (p2s) triangle mapping to proc 1 .. numprocs


        if send_s2p :
            n = len(s2p_map)
            s2p_map_keys_flat = num.reshape(num.array(list(s2p_map.keys()),int), (n,1) )
            s2p_map_values_flat = num.array(list(s2p_map.values()),int)
            s2p_map_flat = num.concatenate( (s2p_map_keys_flat, s2p_map_values_flat), axis=1 )

            n = len(p2s_map)
            p2s_map_keys_flat = num.reshape(num.array(list(p2s_map.keys()),int), (n,2) )
            p2s_map_values_flat = num.reshape(num.array(list(p2s_map.values()),int) , (n,1))
            p2s_map_flat = num.concatenate( (p2s_map_keys_flat, p2s_map_values_flat), axis=1 )

            for p in range(1, numprocs):

                # FIXME SR: Creates cPickle dump
                send(s2p_map_flat, p)
                # FIXME SR: Creates cPickle dump
                #print p2s_map
                send(p2s_map_flat, p)
        else:
            if verbose: print('Not sending s2p_map and p2s_map')
            s2p_map = None
            p2s_map = None

        if verbose: print('Communication done')

    else:
        # Read in the mesh partition that belongs to this
        # processor
        if verbose: print('P%d: Receiving submeshes' %(myid))
        points, vertices, boundary, quantities,\
                ghost_recv_dict, full_send_dict,\
                number_of_full_nodes, number_of_full_triangles, \
                tri_map, node_map, tri_l2g, node_l2g, ghost_layer_width =\
                rec_submesh(0, verbose)



        # Extract l2g maps
        #tri_l2g  = extract_l2g_map(tri_map)
        #node_l2g = extract_l2g_map(node_map)
        #tri_l2g = p2s_map[tri_l2g]

        # Receive serial to parallel (s2p) and parallel to serial (p2s) triangle mapping
        if send_s2p :
            s2p_map_flat = receive(0)
            s2p_map = dict.fromkeys(s2p_map_flat[:,0], s2p_map_flat[:,1:2])

            p2s_map_flat = receive(0)
            p2s_map_keys = [tuple(x) for x in p2s_map_flat[:,0:1]]

            p2s_map = dict.fromkeys(p2s_map_keys, p2s_map_flat[:,2])
        else:
            s2p_map = None
            p2s_map = None

    #------------------------------------------------------------------------
    # Build the domain for this processor using partion structures
    #------------------------------------------------------------------------

    if verbose: print('myid = %g, no_full_nodes = %g, no_full_triangles = %g' % (myid, number_of_full_nodes, number_of_full_triangles))


    domain = Parallel_domain(points, vertices, boundary,
                             full_send_dict=full_send_dict,
                             ghost_recv_dict=ghost_recv_dict,
                             number_of_full_nodes=number_of_full_nodes,
                             number_of_full_triangles=number_of_full_triangles,
                             geo_reference=georef,
                             number_of_global_triangles = number_of_global_triangles,
                             number_of_global_nodes = number_of_global_nodes,
                             s2p_map = s2p_map,
                             p2s_map = p2s_map, ## jj added this
                             tri_l2g = tri_l2g, ## SR added this
                             node_l2g = node_l2g,
                             ghost_layer_width = ghost_layer_width)

    #------------------------------------------------------------------------
    # Transfer initial conditions to each subdomain
    #------------------------------------------------------------------------
    for q in quantities:
        domain.set_quantity(q, quantities[q])


    #------------------------------------------------------------------------
    # Transfer boundary conditions to each subdomain
    #------------------------------------------------------------------------
    boundary_map['ghost'] = None  # Add binding to ghost boundary
    domain.set_boundary(boundary_map)


    #------------------------------------------------------------------------
    # Transfer other attributes to each subdomain
    #------------------------------------------------------------------------
    domain.set_name(domain_name)
    domain.set_datadir(domain_dir)
    domain.set_store(domain_store)
    domain.set_store_centroids(domain_store_centroids)
    domain.set_store_vertices_smoothly(domain_smooth,domain_reduction)
    domain.set_minimum_storable_height(domain_minimum_storable_height)
    domain.set_minimum_allowed_height(domain_minimum_allowed_height)
    domain.set_flow_algorithm(domain_flow_algorithm)
    domain.geo_reference = georef

    #------------------------------------------------------------------------
    # Return parallel domain to all nodes
    #------------------------------------------------------------------------
    return domain






def distribute_mesh(domain, verbose=False, debug=False, parameters=None):
    """ Distribute andsend the mesh info to all the processors.
    Should only be run from processor 0 and will send info to the other
    processors.
    There should be a corresponding  rec_submesh called from all the other
    processors
    """

    if debug:
        verbose = True

    numprocs = size()


    # Subdivide the mesh
    if verbose: print('Subdivide mesh')
    new_nodes, new_triangles, new_boundary, triangles_per_proc, quantities, \
           s2p_map, p2s_map = \
           pmesh_divide_metis_with_map(domain, numprocs)

    #PETE: s2p_map (maps serial domain triangles to parallel domain triangles)
    #      sp2_map (maps parallel domain triangles to domain triangles)



    # Build the mesh that should be assigned to each processor,
    # this includes ghost nodes and the communication pattern
    if verbose: print('Build submeshes')
    submesh = build_submesh(new_nodes, new_triangles, new_boundary, quantities, triangles_per_proc, parameters)

    if verbose:
        for p in range(numprocs):
            N = len(submesh['ghost_nodes'][p])
            M = len(submesh['ghost_triangles'][p])
            print('There are %d ghost nodes and %d ghost triangles on proc %d'\
                  %(N, M, p))

    #if debug:
    #    from pprint import pprint
    #    pprint(submesh)


    # Send the mesh partition to the appropriate processor
    if verbose: print('Distribute submeshes')
    for p in range(1, numprocs):
        send_submesh(submesh, triangles_per_proc, p2s_map, p, verbose)

    # Build the local mesh for processor 0
    points, vertices, boundary, quantities, \
            ghost_recv_dict, full_send_dict, \
            tri_map, node_map, tri_l2g, node_l2g, ghost_layer_width =\
              extract_submesh(submesh, triangles_per_proc, p2s_map, 0)



    # Keep track of the number full nodes and triangles.
    # This is useful later if one needs access to a ghost-free domain
    # Here, we do it for process 0. The others are done in rec_submesh.
    number_of_full_nodes = len(submesh['full_nodes'][0])
    number_of_full_triangles = len(submesh['full_triangles'][0])

    #print
    #for p in range(numprocs):
    #    print 'Process %d:' %(p)
    #
    #    print 'full_triangles:'
    #    print submesh['full_triangles'][p]
    #
    #    print 'full_nodes:'
    #    print submesh['full_nodes'][p]
    #
    #    print 'ghost_triangles:'
    #    print submesh['ghost_triangles'][p]#
    #
    #    print 'ghost_nodes:'
    #   print submesh['ghost_nodes'][p]
    #    print
    #
    #print 'Receive dict'
    #print ghost_recv_dict
    #
    #print 'Send dict'
    #print full_send_dict


    # Return structures necessary for building the parallel domain
    return points, vertices, boundary, quantities,\
           ghost_recv_dict, full_send_dict,\
           number_of_full_nodes, number_of_full_triangles, \
           s2p_map, p2s_map, tri_map, node_map, tri_l2g, node_l2g, \
           ghost_layer_width



## def extract_l2g_map(map):
##     # Extract l2g data  from corresponding map
##     # Maps

##     import numpy as num

##     # FIXME: this is where we loss the original order of
##     # sequential domain
##     b = num.arange(len(map))

##     l_ids = num.extract(map>-1,map)
##     g_ids = num.extract(map>-1,b)

## #    print len(g_ids)
## #    print len(l_ids)
## #    print l_ids

##     l2g = num.zeros_like(g_ids)
##     l2g[l_ids] = g_ids

##     return l2g

def mpicmd(script_name='echo', numprocs=3):

    extra_options = mpi_extra_options()

    return "mpiexec -np %d  %s  python -m mpi4py %s" % (numprocs, extra_options, script_name)  

def mpi_extra_options():

    extra_options = '--oversubscribe'
    cmd = 'mpiexec -np 3 ' + extra_options + ' echo '

    #print(cmd)
    import subprocess
    result = subprocess.run(cmd.split(), capture_output=True)
    if result.returncode != 0:
        extra_options = ' '

    import platform
    if platform.system() == 'Windows':
        extra_options = ' '

    return extra_options