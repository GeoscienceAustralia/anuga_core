"""Trying to lump parallel stuff into simpler interface


"""



# The abstract Python-MPI interface
from anuga_parallel.parallel_abstraction import size, rank, get_processor_name
from anuga_parallel.parallel_abstraction import finalize, send, receive
from anuga_parallel.parallel_abstraction import pypar_available, barrier


# ANUGA parallel engine (only load if pypar can)
if pypar_available:
    from anuga_parallel.build_commun  import send_submesh
    from anuga_parallel.build_commun  import rec_submesh
    from anuga_parallel.build_commun  import extract_hostmesh
    from anuga_parallel.parallel_shallow_water import Parallel_Domain

    # Mesh partitioning using Metis
    from anuga_parallel.build_submesh import build_submesh
    from anuga_parallel.pmesh_divide  import pmesh_divide_metis


#------------------------------------------------------------------------------
# Read in processor information
#------------------------------------------------------------------------------

numprocs = size()
myid = rank()
processor_name = get_processor_name()
print 'I am processor %d of %d on node %s' %(myid, numprocs, processor_name)




def distribute(domain, verbose=False):
    """ Distribute the domain to all processes
    """


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
                number_of_full_nodes, number_of_full_triangles =\
                distribute_mesh(domain)


        if verbose: print 'Communication done'
        
    else:
        # Read in the mesh partition that belongs to this
        # processor
        if verbose: print 'P%d: Receiving submeshes' %(myid)                
        points, vertices, boundary, quantities,\
                ghost_recv_dict, full_send_dict,\
                number_of_full_nodes, number_of_full_triangles =\
                rec_submesh(0)


    #------------------------------------------------------------------------
    # Build the domain for this processor using partion structures
    #------------------------------------------------------------------------

    if verbose: print 'myid ',myid, number_of_full_nodes, number_of_full_triangles

    
    domain = Parallel_Domain(points, vertices, boundary,
                             full_send_dict=full_send_dict,
                             ghost_recv_dict=ghost_recv_dict,
                             number_of_full_nodes=number_of_full_nodes,
                             number_of_full_triangles=number_of_full_triangles)

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






def distribute_mesh(domain):

    numprocs = size()

    
    # Subdivide the mesh
    print 'Subdivide mesh'
    nodes, triangles, boundary, triangles_per_proc, quantities = \
           pmesh_divide_metis(domain, numprocs)




    # Build the mesh that should be assigned to each processor,
    # this includes ghost nodes and the communication pattern
    print 'Build submeshes'    
    submesh = build_submesh(nodes, triangles, boundary,\
                            quantities, triangles_per_proc)

    for p in range(numprocs):
        N = len(submesh['ghost_nodes'][p])                
        M = len(submesh['ghost_triangles'][p])
        print 'There are %d ghost nodes and %d ghost triangles on proc %d'\
              %(N, M, p)


    # Send the mesh partition to the appropriate processor
    print 'Distribute submeshes'        
    for p in range(1, numprocs):
      send_submesh(submesh, triangles_per_proc, p)

    # Build the local mesh for processor 0
    points, vertices, boundary, quantities, ghost_recv_dict, full_send_dict =\
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
           number_of_full_nodes, number_of_full_triangles
    



