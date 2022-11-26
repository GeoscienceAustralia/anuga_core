"""Trying to lump parallel stuff into simpler interface


"""

import numpy as num

from anuga import Domain

from anuga.parallel.distribute_mesh  import send_submesh
from anuga.parallel.distribute_mesh  import rec_submesh
from anuga.parallel.distribute_mesh  import extract_submesh

# Mesh partitioning using Metis
from anuga.parallel.distribute_mesh import build_submesh
from anuga.parallel.distribute_mesh import pmesh_divide_metis_with_map

from anuga.parallel.parallel_shallow_water import Parallel_domain



class Sequential_distribute(object):

    def __init__(self, domain, verbose=False, debug=False, parameters=None):

        if debug:
            verbose = True

        self.domain = domain
        self.verbose = verbose
        self.debug = debug
        self.parameters = parameters


    def distribute(self, numprocs=1):

        self.numprocs = numprocs

        domain = self.domain
        verbose = self.verbose
        debug = self.debug
        parameters = self.parameters

        # FIXME: Dummy assignment (until boundaries are refactored to
        # be independent of domains until they are applied)
        bdmap = {}
        for tag in domain.get_boundary_tags():
            bdmap[tag] = None

        domain.set_boundary(bdmap)


        self.domain_name = domain.get_name()
        self.domain_dir = domain.get_datadir()
        self.domain_store = domain.get_store()
        self.domain_store_centroids = domain.get_store_centroids()
        self.domain_minimum_storable_height = domain.minimum_storable_height
        self.domain_flow_algorithm = domain.get_flow_algorithm()
        self.domain_minimum_allowed_height = domain.get_minimum_allowed_height()
        self.domain_georef = domain.geo_reference
        self.domain_quantities_to_be_stored = domain.quantities_to_be_stored
        self.domain_smooth = domain.smooth
        self.domain_low_froude = domain.low_froude
        self.number_of_global_triangles = domain.number_of_triangles
        self.number_of_global_nodes = domain.number_of_nodes
        self.boundary_map = domain.boundary_map


        # Subdivide the mesh
        if verbose: print('sequential_distribute: Subdivide mesh')

        new_nodes, new_triangles, new_boundary, triangles_per_proc, quantities, \
               s2p_map, p2s_map = \
               pmesh_divide_metis_with_map(domain, numprocs)


        # Build the mesh that should be assigned to each processor,
        # this includes ghost nodes and the communication pattern
        if verbose: print('sequential_distribute: Build submeshes')
        if verbose: print('sequential_distribute: parameters = ',parameters)

        submesh = build_submesh(new_nodes, new_triangles, new_boundary, \
                                quantities, triangles_per_proc, parameters=parameters)

        if verbose:
            for p in range(numprocs):
                N = len(submesh['ghost_nodes'][p])
                M = len(submesh['ghost_triangles'][p])
                print('There are %d ghost nodes and %d ghost triangles on proc %d'\
                      %(N, M, p))


        self.submesh = submesh
        self.triangles_per_proc = triangles_per_proc
        self.p2s_map =  p2s_map


    def extract_submesh(self, p=0):
        """Build the local mesh for processor p
        """

        submesh = self.submesh
        triangles_per_proc = self.triangles_per_proc
        p2s_map = self.p2s_map
        verbose = self.verbose
        debug = self.debug

        assert p>=0
        assert p<self.numprocs


        points, vertices, boundary, quantities, \
            ghost_recv_dict, full_send_dict, \
            tri_map, node_map, tri_l2g, node_l2g, ghost_layer_width =\
              extract_submesh(submesh, triangles_per_proc, p2s_map, p)


        number_of_full_nodes = len(submesh['full_nodes'][p])
        number_of_full_triangles = len(submesh['full_triangles'][p])


        if debug:
            import pprint
            print(50*"=")
            print('NODE_L2G')
            pprint.pprint(node_l2g)

            pprint.pprint(node_l2g[vertices[:,0]])

            print('VERTICES')
            pprint.pprint(vertices[:,0])
            pprint.pprint(new_triangles[tri_l2g,0])

            assert num.allclose(node_l2g[vertices[:,0]], new_triangles[tri_l2g,0])
            assert num.allclose(node_l2g[vertices[:,1]], new_triangles[tri_l2g,1])
            assert num.allclose(node_l2g[vertices[:,2]], new_triangles[tri_l2g,2])


            print('POINTS')
            pprint.pprint(points)

            assert num.allclose(points[:,0], new_nodes[node_l2g,0])
            assert num.allclose(points[:,1], new_nodes[node_l2g,1])


            print('TRI')
            pprint.pprint(tri_l2g)
            pprint.pprint(p2s_map[tri_l2g])


            assert num.allclose(original_triangles[tri_l2orig,0],node_l2g[vertices[:,0]])
            assert num.allclose(original_triangles[tri_l2orig,1],node_l2g[vertices[:,1]])
            assert num.allclose(original_triangles[tri_l2orig,2],node_l2g[vertices[:,2]])

            print('NODES')
            pprint.pprint(node_map)
            pprint.pprint(node_l2g)

        #tri_l2orig = p2s_map[tri_l2g]

        s2p_map = None
        p2s_map = None

        #------------------------------------------------------------------------
        # Build the parallel domain for this processor using partion structures
        #------------------------------------------------------------------------

        if verbose:
            print('sequential_distribute: P%g, no_full_nodes = %g, no_full_triangles = %g' % (p, number_of_full_nodes, number_of_full_triangles))


        kwargs = {'full_send_dict': full_send_dict,
                'ghost_recv_dict': ghost_recv_dict,
                'number_of_full_nodes': number_of_full_nodes,
                'number_of_full_triangles': number_of_full_triangles,
                'geo_reference': self.domain_georef,
                'number_of_global_triangles':  self.number_of_global_triangles,
                'number_of_global_nodes':  self.number_of_global_nodes,
                'processor':  p,
                'numproc':  self.numprocs,
                's2p_map':  s2p_map,
                'p2s_map':  p2s_map, ## jj added this
                'tri_l2g':  tri_l2g, ## SR added this
                'node_l2g':  node_l2g,
                'ghost_layer_width':  ghost_layer_width}

        boundary_map = self.boundary_map
        domain_name = self.domain_name
        domain_dir = self.domain_dir
        domain_store = self.domain_store
        domain_store_centroids = self.domain_store_centroids
        domain_minimum_storable_height = self.domain_minimum_storable_height
        domain_minimum_allowed_height = self.domain_minimum_allowed_height
        domain_flow_algorithm = self.domain_flow_algorithm
        domain_georef = self.domain_georef
        domain_quantities_to_be_stored = self.domain_quantities_to_be_stored
        domain_smooth = self.domain_smooth
        domain_low_froude = self.domain_low_froude

        tostore = (kwargs, points, vertices, boundary, quantities, \
                   boundary_map, \
                   domain_name, domain_dir, domain_store, domain_store_centroids, \
                   domain_minimum_storable_height, \
                   domain_minimum_allowed_height, domain_flow_algorithm, \
                   domain_georef, domain_quantities_to_be_stored, domain_smooth, \
                   domain_low_froude)

        return tostore






def sequential_distribute_dump(domain, numprocs=1, verbose=False, partition_dir='.', debug=False, parameters = None):
    """ Distribute the domain, create parallel domain and pickle result
    """

    from os.path import join

    partition = Sequential_distribute(domain, verbose, debug, parameters)

    partition.distribute(numprocs)

    # Make sure the partition_dir exists
    if partition_dir == '.' :
        pass
    else:
        import os
        import errno
        try:
            os.makedirs(partition_dir)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

    import pickle
    for p in range(0, numprocs):

        tostore = partition.extract_submesh(p)

        pickle_name = partition.domain_name + '_P%g_%g.pickle'% (numprocs,p)
        pickle_name = join(partition_dir,pickle_name)
        f = open(pickle_name, 'wb')

        lst = list(tostore)

        # Write points and triangles to their own files
        num.save(pickle_name+".np1",tostore[1]) # this append .npy to filename
        lst[1] = pickle_name+".np1.npy"
        num.save(pickle_name+".np2",tostore[2])
        lst[2] = pickle_name+".np2.npy"

        # Write each quantity to it's own file
        for k in tostore[4]:
            num.save(pickle_name+".np4."+k,num.array(tostore[4][k]))
            lst[4][k] = pickle_name+".np4."+k+".npy"

        pickle.dump( tuple(lst), f, protocol=pickle.HIGHEST_PROTOCOL)
    return


def sequential_distribute_load(filename = 'domain', partition_dir = '.', verbose = False):


    from anuga import myid, numprocs

    from os.path import join

    pickle_name = filename+'_P%g_%g.pickle'% (numprocs,myid)
    pickle_name = join(partition_dir,pickle_name)

    return sequential_distribute_load_pickle_file(pickle_name, numprocs, verbose = verbose)


def sequential_distribute_load_pickle_file(pickle_name, np=1, verbose = False):
    """
    Open pickle files
    """

    f = open(pickle_name, 'rb')
    import pickle

    kwargs, points, vertices, boundary, quantities, boundary_map, \
                   domain_name, domain_dir, domain_store, domain_store_centroids, \
                   domain_minimum_storable_height, domain_minimum_allowed_height, \
                   domain_flow_algorithm, domain_georef, \
                   domain_quantities_to_be_stored, domain_smooth, \
                   domain_low_froude = pickle.load(f)
    f.close()

    for k in quantities:
        quantities[k] = num.load(quantities[k])
    points = num.load(points)
    vertices = num.load(vertices)

    #---------------------------------------------------------------------------
    # Create domain (parallel if np>1)
    #---------------------------------------------------------------------------
    if np>1:
        domain = Parallel_domain(points, vertices, boundary, **kwargs)
    else:
        domain = Domain(points, vertices, boundary, **kwargs)

    #------------------------------------------------------------------------
    # Copy in quantity data
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
    domain.set_flow_algorithm(domain_flow_algorithm)
    domain.set_low_froude(domain_low_froude)
    domain.set_store(domain_store)
    domain.set_store_centroids(domain_store_centroids)
    domain.set_minimum_storable_height(domain_minimum_storable_height)
    domain.set_minimum_allowed_height(domain_minimum_allowed_height)
    domain.geo_reference = domain_georef
    domain.set_quantities_to_be_stored(domain_quantities_to_be_stored)
    domain.smooth = domain_smooth

    return domain
