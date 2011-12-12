"""Trying to lump parallel stuff into simpler interface


"""



# The abstract Python-MPI interface
from anuga_parallel.parallel_abstraction import size, rank, get_processor_name
from anuga_parallel.parallel_abstraction import finalize, send, receive
from anuga_parallel.parallel_abstraction import pypar_available, barrier


# ANUGA parallel engine (only load if pypar can)
if pypar_available:
    from anuga_parallel.distribute_mesh  import send_submesh
    from anuga_parallel.distribute_mesh  import rec_submesh
    from anuga_parallel.distribute_mesh  import extract_hostmesh

    # Mesh partitioning using Metis
    from anuga_parallel.distribute_mesh import build_submesh
    from anuga_parallel.distribute_mesh import pmesh_divide_metis_with_map

    from anuga_parallel.parallel_shallow_water import Parallel_domain

#------------------------------------------------------------------------------
# Read in processor information
#------------------------------------------------------------------------------

numprocs = size()
myid = rank()
processor_name = get_processor_name()
#print 'I am processor %d of %d on node %s' %(myid, numprocs, processor_name)




def distribute(domain, verbose=False):
    """ Distribute the domain to all processes
    """

    barrier()

    # FIXME: Dummy assignment (until boundaries are refactored to
    # be independent of domains until they are applied)
    if myid == 0:
        bdmap = {}
        for tag in domain.get_boundary_tags():
            bdmap[tag] = None
    
    
        domain.set_boundary(bdmap)




    if not pypar_available: return domain # Bypass

    # For some obscure reason this communication must happen prior to
    # the more complex mesh distribution - Oh Well!
    if myid == 0:
        domain_name = domain.get_name()
        domain_dir = domain.get_datadir()
        georef = domain.geo_reference
        
        # FIXME - what other attributes need to be transferred?

        for p in range(1, numprocs):
            send((domain_name, domain_dir, georef), p)
    else:
        if verbose: print 'P%d: Receiving domain attributes' %(myid)

        domain_name, domain_dir, georef = receive(0)



    # Distribute boundary conditions
    # FIXME: This cannot handle e.g. Time_boundaries due to
    # difficulties pickling functions
    if myid == 0:
        boundary_map = domain.boundary_map
        for p in range(1, numprocs):
            send(boundary_map, p)
    else:
        if verbose: print 'P%d: Receiving boundary map' %(myid)        

        boundary_map = receive(0)
        



    if myid == 0:
        # Partition and distribute mesh.
        # Structures returned is in the
        # correct form for the ANUGA data structure


        points, vertices, boundary, quantities,\
                ghost_recv_dict, full_send_dict,\
                number_of_full_nodes, number_of_full_triangles,\
                s2p_map, p2s_map, tri_map, node_map =\
                distribute_mesh(domain, verbose=verbose)

        number_of_global_triangles = len(tri_map)
        number_of_global_nodes = len(node_map)

        # Extract l2g maps
        tri_l2g  = extract_l2g_map(tri_map)
        node_l2g = extract_l2g_map(node_map)

        # Send serial to parallel (s2p) and parallel to serial (p2s) triangle mapping to proc 1 .. numprocs
        for p in range(1, numprocs):
            send(s2p_map, p)
            send(p2s_map, p)

        if verbose: print 'Communication done'
        
    else:
        # Read in the mesh partition that belongs to this
        # processor
        if verbose: print 'P%d: Receiving submeshes' %(myid)                
        points, vertices, boundary, quantities,\
                ghost_recv_dict, full_send_dict,\
                number_of_full_nodes, number_of_full_triangles, \
                tri_map, node_map =\
                rec_submesh(0, verbose)

        number_of_global_triangles = len(tri_map)
        number_of_global_nodes = len(node_map)

        # Extract l2g maps
        tri_l2g  = extract_l2g_map(tri_map)
        node_l2g = extract_l2g_map(node_map)
        
        # Recieve serial to parallel (s2p) and parallel to serial (p2s) triangle mapping
        s2p_map = receive(0)
        p2s_map = receive(0)


    #------------------------------------------------------------------------
    # Build the domain for this processor using partion structures
    #------------------------------------------------------------------------

    if verbose: print 'myid = %g, no_full_nodes = %g, no_full_triangles = %g' % (myid, number_of_full_nodes, number_of_full_triangles)

    
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
                             node_l2g = node_l2g)

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
    domain.geo_reference = georef   

    #------------------------------------------------------------------------
    # Return parallel domain to all nodes
    #------------------------------------------------------------------------
    return domain    






def distribute_mesh(domain, verbose=False):

    numprocs = size()

    
    # Subdivide the mesh
    if verbose: print 'Subdivide mesh'
    nodes, triangles, boundary, triangles_per_proc, quantities, \
           s2p_map, p2s_map = \
           pmesh_divide_metis_with_map(domain, numprocs)

    #PETE: s2p_map (maps serial domain triangles to parallel domain triangles)
    #      sp2_map (maps parallel domain triangles to domain triangles)


    # Build the mesh that should be assigned to each processor,
    # this includes ghost nodes and the communication pattern
    if verbose: print 'Build submeshes'    
    submesh = build_submesh(nodes, triangles, boundary,\
                            quantities, triangles_per_proc)

    if verbose:
        for p in range(numprocs):
            N = len(submesh['ghost_nodes'][p])                
            M = len(submesh['ghost_triangles'][p])
            print 'There are %d ghost nodes and %d ghost triangles on proc %d'\
                  %(N, M, p)


    # Send the mesh partition to the appropriate processor
    if verbose: print 'Distribute submeshes'        
    for p in range(1, numprocs):
      send_submesh(submesh, triangles_per_proc, p, verbose)

    # Build the local mesh for processor 0
    points, vertices, boundary, quantities, \
            ghost_recv_dict, full_send_dict, tri_map, node_map =\
              extract_hostmesh(submesh, triangles_per_proc)

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
           s2p_map, p2s_map, tri_map, node_map
    


def extract_l2g_map(map):
    # Extract l2g_map

    import numpy as num
    
    b = num.arange(len(map))

    l_ids = num.extract(map>-1,map)
    g_ids = num.extract(map>-1,b)

#    print len(g_ids)
#    print len(l_ids)
#    print l_ids

    l2g = num.zeros_like(g_ids)
    l2g[l_ids] = g_ids

    return l2g

