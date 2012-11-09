"""Trying to lump parallel stuff into simpler interface


"""

import numpy as num

from anuga_parallel.distribute_mesh  import send_submesh
from anuga_parallel.distribute_mesh  import rec_submesh
from anuga_parallel.distribute_mesh  import extract_submesh

# Mesh partitioning using Metis
from anuga_parallel.distribute_mesh import build_submesh
from anuga_parallel.distribute_mesh import pmesh_divide_metis_with_map

from anuga_parallel.parallel_shallow_water import Parallel_domain



def sequential_distribute_dump(domain, numprocs=1, verbose=False, debug=False, parameters = None):
    """ Distribute the domain, create parallel domain and pickle result
    """


    if debug:
        verbose = True



    # FIXME: Dummy assignment (until boundaries are refactored to
    # be independent of domains until they are applied)
    bdmap = {}
    for tag in domain.get_boundary_tags():
        bdmap[tag] = None

    domain.set_boundary(bdmap)


    if numprocs == 1 : return # Bypass


    domain_name = domain.get_name()
    domain_dir = domain.get_datadir()
    domain_store = domain.get_store()
    domain_minimum_storable_height = domain.minimum_storable_height
    domain_flow_algorithm = domain.get_flow_algorithm()
    domain_minimum_allowed_height = domain.get_minimum_allowed_height()
    georef = domain.geo_reference
    number_of_global_triangles = domain.number_of_triangles
    number_of_global_nodes = domain.number_of_nodes
    boundary_map = domain.boundary_map


    #sequential_distribute_mesh(domain, numprocs, verbose=verbose, debug=debug, parameters=parameters)


    # Subdivide the mesh
    if verbose: print 'sequential_distribute: Subdivide mesh'
    new_nodes, new_triangles, new_boundary, triangles_per_proc, quantities, \
           s2p_map, p2s_map = \
           pmesh_divide_metis_with_map(domain, numprocs)

    #PETE: s2p_map (maps serial domain triangles to parallel domain triangles)
    #      sp2_map (maps parallel domain triangles to domain triangles)



    # Build the mesh that should be assigned to each processor,
    # this includes ghost nodes and the communication pattern
    if verbose: print 'sequential_distribute: Build submeshes'
    submesh = build_submesh(new_nodes, new_triangles, new_boundary, quantities, triangles_per_proc, parameters)

    if debug:
        for p in range(numprocs):
            N = len(submesh['ghost_nodes'][p])
            M = len(submesh['ghost_triangles'][p])
            print 'There are %d ghost nodes and %d ghost triangles on proc %d'\
                  %(N, M, p)

    #if debug:
    #    from pprint import pprint
    #    pprint(submesh)


    # extract data to create parallel domain
    if verbose: print 'sequential_distribute: Distribute submeshes'
    for p in range(0, numprocs):

        # Build the local mesh for processor 0
        points, vertices, boundary, quantities, \
            ghost_recv_dict, full_send_dict, tri_map, node_map, ghost_layer_width =\
              extract_submesh(submesh, triangles_per_proc, p)


#        from pprint import pprint
#        print '='*80
#        print p
#        print '='*80
#        pprint(tri_map)
#        print len(tri_map)

        # Keep track of the number full nodes and triangles.
        # This is useful later if one needs access to a ghost-free domain
        # Here, we do it for process 0. The others are done in rec_submesh.
        number_of_full_nodes = len(submesh['full_nodes'][p])
        number_of_full_triangles = len(submesh['full_triangles'][p])

        # Extract l2g maps
        tri_l2g  = extract_l2g_map(tri_map)
        node_l2g = extract_l2g_map(node_map)


        s2p_map = None
        p2s_map = None

        #------------------------------------------------------------------------
        # Build the parallel domain for this processor using partion structures
        #------------------------------------------------------------------------

        if verbose:
            print 'sequential_distribute: P%g, no_full_nodes = %g, no_full_triangles = %g' % (p, number_of_full_nodes, number_of_full_triangles)


        #args = [points, vertices, boundary]

        kwargs = {'full_send_dict': full_send_dict,
                'ghost_recv_dict': ghost_recv_dict,
                'number_of_full_nodes': number_of_full_nodes,
                'number_of_full_triangles': number_of_full_triangles,
                'geo_reference': georef,
                'number_of_global_triangles':  number_of_global_triangles,
                'number_of_global_nodes':  number_of_global_nodes,
                'processor':  p,
                'numproc':  numprocs,
                's2p_map':  s2p_map,
                'p2s_map':  p2s_map, ## jj added this
                'tri_l2g':  tri_l2g, ## SR added this
                'node_l2g':  node_l2g,
                'ghost_layer_width':  ghost_layer_width}

#        parallel_domain = Parallel_domain(points, vertices, boundary, **kwargs)



        #------------------------------------------------------------------------
        # Transfer initial conditions to each subdomain
        #------------------------------------------------------------------------
#        for q in quantities:
#            parallel_domain.set_quantity(q, quantities[q])


        #------------------------------------------------------------------------
        # Transfer boundary conditions to each subdomain
        #------------------------------------------------------------------------
#        boundary_map['ghost'] = None  # Add binding to ghost boundary
#        parallel_domain.set_boundary(boundary_map)


        #------------------------------------------------------------------------
        # Transfer other attributes to each subdomain
        #------------------------------------------------------------------------
#        parallel_domain.set_name(domain_name)
#        parallel_domain.set_datadir(domain_dir)
#        parallel_domain.set_store(domain_store)
#        parallel_domain.set_minimum_storable_height(domain_minimum_storable_height)
#        parallel_domain.set_minimum_allowed_height(domain_minimum_allowed_height)
#        parallel_domain.set_flow_algorithm(domain_flow_algorithm)
#        parallel_domain.geo_reference = georef



        #-----------------------------------------------------------------------
        # Now let's store the parallel_domain via cPickle
        #-----------------------------------------------------------------------
#        import cPickle
#        pickle_name = domain_name + '_P%g_%g.pickle'% (numprocs,p)
#        f = file(pickle_name, 'wb')
#        cPickle.dump(parallel_domain, f, protocol=cPickle.HIGHEST_PROTOCOL)
#        f.close()


        #FIXME SR: Looks like we could reduce storage by a factor of 4 by just
        # storing the data to create the parallel_domain instead of pickling
        # a created domain
        import cPickle
        pickle_name = domain_name + '_P%g_%g.pickle'% (numprocs,p)
        f = file(pickle_name, 'wb')
        tostore = (kwargs, points, vertices, boundary, quantities, boundary_map, domain_name, domain_dir, domain_store, domain_minimum_storable_height, \
                   domain_minimum_allowed_height, domain_flow_algorithm, georef)
        cPickle.dump( tostore, f, protocol=cPickle.HIGHEST_PROTOCOL)

    return


def sequential_distribute_load(filename = 'domain', verbose = False):


    from anuga_parallel import myid, numprocs


    #---------------------------------------------------------------------------
    # Open pickle files
    #---------------------------------------------------------------------------
    import cPickle
    pickle_name = filename+'_P%g_%g.pickle'% (numprocs,myid)
    f = file(pickle_name, 'rb')
    kwargs, points, vertices, boundary, quantities, boundary_map, domain_name, domain_dir, domain_store, domain_minimum_storable_height, \
                   domain_minimum_allowed_height, domain_flow_algorithm, georef = cPickle.load(f)
    f.close()

    #---------------------------------------------------------------------------
    # Create parallel domain
    #---------------------------------------------------------------------------
    parallel_domain = Parallel_domain(points, vertices, boundary, **kwargs)


    #------------------------------------------------------------------------
    # Copy in quantity data
    #------------------------------------------------------------------------
    for q in quantities:
        parallel_domain.set_quantity(q, quantities[q])


    #------------------------------------------------------------------------
    # Transfer boundary conditions to each subdomain
    #------------------------------------------------------------------------
    boundary_map['ghost'] = None  # Add binding to ghost boundary
    parallel_domain.set_boundary(boundary_map)


    #------------------------------------------------------------------------
    # Transfer other attributes to each subdomain
    #------------------------------------------------------------------------
    parallel_domain.set_name(domain_name)
    parallel_domain.set_datadir(domain_dir)
    parallel_domain.set_store(domain_store)
    parallel_domain.set_minimum_storable_height(domain_minimum_storable_height)
    parallel_domain.set_minimum_allowed_height(domain_minimum_allowed_height)
    parallel_domain.set_flow_algorithm(domain_flow_algorithm)
    parallel_domain.geo_reference = georef


    return parallel_domain

def extract_l2g_map(map):
    # Extract l2g data  from corresponding map
    # Maps

    import numpy as num

    b = num.arange(len(map))

    l_ids = num.extract(map>-1,map)
    g_ids = num.extract(map>-1,b)


#    print len(g_ids)
#    print len(l_ids)
#    print l_ids
#    print g_ids

    l2g = num.zeros_like(g_ids)
    l2g[l_ids] = g_ids

    return l2g





