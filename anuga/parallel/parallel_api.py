"""Trying to lump parallel stuff into simpler interface


"""
import numpy as num

# The abstract Python-MPI interface
from anuga.utilities.parallel_abstraction import size, rank, get_processor_name
from anuga.utilities.parallel_abstraction import finalize, send, receive, reduce
from anuga.utilities.parallel_abstraction import isend, waitall
from anuga.utilities.parallel_abstraction import pypar_available, barrier

from anuga.parallel.sequential_distribute import sequential_distribute_dump
from anuga.parallel.sequential_distribute import sequential_distribute_load
from anuga.parallel.sequential_distribute import sequential_mesh_dump
from anuga.parallel.sequential_distribute import sequential_mesh_load
from anuga.parallel.sequential_distribute import uniform_refine_domain
from anuga.parallel.sequential_distribute import sequential_mesh_refine
from anuga.parallel.sequential_distribute import create_parallel_mesh

# ANUGA parallel engine (only load if pypar can)
if pypar_available:
    from anuga.parallel.distribute_mesh  import send_submesh
    from anuga.parallel.distribute_mesh  import rec_submesh
    from anuga.parallel.distribute_mesh  import extract_submesh

    # Mesh partitioning using Metis
    from anuga.parallel.distribute_mesh import build_submesh
    from anuga.parallel.distribute_mesh import partition_mesh
    from anuga.parallel.distribute_mesh import ghost_commun_pattern
    from anuga.parallel.distribute_mesh import build_local_mesh

    from anuga.parallel.parallel_shallow_water import Parallel_domain



from anuga.abstract_2d_finite_volumes.neighbour_mesh import Mesh


def _isend_submesh(comm, tostore, dest):
    """Send an extracted submesh tuple to `dest` using buffer-protocol Isend
    for the large numpy arrays and pickle for the small metadata.

    The large arrays (points, vertices, quantities, tri_l2g, node_l2g) are
    sent via ``comm.Isend()`` (zero-copy on RDMA networks).  Everything else
    is bundled into a single small metadata pickle sent via ``comm.isend()``.

    Returns a list of MPI Request objects.  The caller must call
    ``MPI.Request.Waitall(reqs)`` before allowing the source arrays to be
    garbage-collected.
    """
    from mpi4py import MPI

    (kwargs, points, vertices, boundary, quantities, boundary_map,
     domain_name, domain_dir, domain_store, domain_store_centroids,
     domain_minimum_storable_height, domain_minimum_allowed_height,
     domain_flow_algorithm, domain_georef,
     domain_quantities_to_be_stored, domain_smooth, domain_low_froude) = tostore

    quant_keys = list(quantities.keys())
    nq = len(quant_keys)

    # Pull the large numpy arrays out of kwargs (copy so we don't mutate).
    tri_l2g  = kwargs['tri_l2g']
    node_l2g = kwargs['node_l2g']
    kwargs_small = {k: v for k, v in kwargs.items()
                    if k not in ('tri_l2g', 'node_l2g')}

    # Tag layout:  0=metadata  1=points  2=vertices  3..3+nq-1=quantities
    #              3+nq=tri_l2g  3+nq+1=node_l2g
    meta = (kwargs_small, quant_keys,
            boundary, boundary_map,
            domain_name, domain_dir, domain_store, domain_store_centroids,
            domain_minimum_storable_height, domain_minimum_allowed_height,
            domain_flow_algorithm, domain_georef,
            domain_quantities_to_be_stored, domain_smooth, domain_low_froude,
            points.shape,   points.dtype.str,
            vertices.shape, vertices.dtype.str,
            tri_l2g.shape,  tri_l2g.dtype.str,
            node_l2g.shape, node_l2g.dtype.str)

    reqs = [comm.isend(meta, dest=dest, tag=0)]
    reqs.append(comm.Isend(num.ascontiguousarray(points),   dest=dest, tag=1))
    reqs.append(comm.Isend(num.ascontiguousarray(vertices), dest=dest, tag=2))
    for ki, k in enumerate(quant_keys):
        reqs.append(comm.Isend(num.ascontiguousarray(quantities[k]),
                               dest=dest, tag=3 + ki))
    reqs.append(comm.Isend(num.ascontiguousarray(tri_l2g),  dest=dest, tag=3 + nq))
    reqs.append(comm.Isend(num.ascontiguousarray(node_l2g), dest=dest, tag=3 + nq + 1))
    return reqs


def _recv_submesh(comm, source):
    """Receive an extracted submesh from `source`.

    Counterpart to :func:`_isend_submesh`.  Blocks on the metadata receive
    (needed for array shapes), then issues blocking Recv calls for each
    large array.  Returns the same tuple as ``extract_submesh()``.
    """
    (kwargs_small, quant_keys,
     boundary, boundary_map,
     domain_name, domain_dir, domain_store, domain_store_centroids,
     domain_minimum_storable_height, domain_minimum_allowed_height,
     domain_flow_algorithm, domain_georef,
     domain_quantities_to_be_stored, domain_smooth, domain_low_froude,
     points_shape,   points_dtype,
     vertices_shape, vertices_dtype,
     tri_l2g_shape,  tri_l2g_dtype,
     node_l2g_shape, node_l2g_dtype) = comm.recv(source=source, tag=0)

    nq = len(quant_keys)

    points   = num.empty(points_shape,   dtype=num.dtype(points_dtype))
    vertices = num.empty(vertices_shape, dtype=num.dtype(vertices_dtype))
    comm.Recv(points,   source=source, tag=1)
    comm.Recv(vertices, source=source, tag=2)

    quantities = {}
    for ki, k in enumerate(quant_keys):
        arr = num.empty(vertices_shape[:1] + (1,), dtype=num.float64)
        comm.Recv(arr, source=source, tag=3 + ki)
        quantities[k] = arr

    tri_l2g  = num.empty(tri_l2g_shape,  dtype=num.dtype(tri_l2g_dtype))
    node_l2g = num.empty(node_l2g_shape, dtype=num.dtype(node_l2g_dtype))
    comm.Recv(tri_l2g,  source=source, tag=3 + nq)
    comm.Recv(node_l2g, source=source, tag=3 + nq + 1)

    kwargs = dict(kwargs_small)
    kwargs['tri_l2g']  = tri_l2g
    kwargs['node_l2g'] = node_l2g

    return (kwargs, points, vertices, boundary, quantities, boundary_map,
            domain_name, domain_dir, domain_store, domain_store_centroids,
            domain_minimum_storable_height, domain_minimum_allowed_height,
            domain_flow_algorithm, domain_georef,
            domain_quantities_to_be_stored, domain_smooth, domain_low_froude)


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
    """Distribute the domain to all processes in a parallel computing environment.

    This function partitions a computational domain across multiple processes and
    creates parallel domain instances on each process. The master process (myid=0)
    handles the domain partitioning and distributes submeshes to worker processes,
    while worker processes receive their respective submesh data.

    Parameters
    ----------
    domain : Domain
        The computational domain to be distributed across processes.
    verbose : bool, optional
        If True, print detailed information during distribution process.
        Default is False.
    debug : bool, optional
        If True, enable debug mode for additional diagnostics.
        Default is False.
    parameters : dict, optional
        User-defined parameters to customize distribution behavior,
        particularly the size of ghost layers. Default is None.

    Returns
    -------
    Parallel_domain or Domain
        A Parallel_domain object containing the distributed submesh data and
        quantities specific to the current process. If mpi4py is not available
        or only one process is running, returns the original domain unchanged.

    Notes
    -----
    - Requires mpi4py for parallel communication
    - Only functions with multiple processes (numprocs > 1)
    - The master process (myid=0) partitions the domain and distributes data
    - All other processes receive submesh data from the master process
    - Boundary conditions are transferred as ghost boundaries
    - Domain attributes (flow algorithm, georeferencing, etc.) are preserved
      in the parallel domain instances
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

        from mpi4py import MPI as _MPI
        all_reqs = []
        for p in range(1, numprocs):
            tostore = partition.extract_submesh(p)
            all_reqs.extend(_isend_submesh(_MPI.COMM_WORLD, tostore, p))
        _MPI.Request.Waitall(all_reqs)

    else:

        from mpi4py import MPI as _MPI
        kwargs, points, vertices, boundary, quantities, boundary_map, \
            domain_name, domain_dir, domain_store, domain_store_centroids, \
            domain_minimum_storable_height, domain_minimum_allowed_height, \
            domain_flow_algorithm, domain_georef, \
            domain_quantities_to_be_stored, domain_smooth, domain_low_froude \
             = _recv_submesh(_MPI.COMM_WORLD, 0)

    #---------------------------------------------------------------------------
    # Now Create parallel domain
    #---------------------------------------------------------------------------
    # Tracers are not Quantity objects (#276), so the loop below cannot see
    # them. They travelled as reserved entries in kwargs (names, shared beta)
    # and in quantities (the concentrations, partitioned by the same machinery
    # as everything else); take both out before the constructor and the loop
    # meet them. #278
    from anuga.parallel.distribute_mesh import (
        TRACER_METADATA_KEY, pop_tracer_quantities, restore_tracers)
    tracer_metadata = kwargs.pop(TRACER_METADATA_KEY, None)
    tracer_values = pop_tracer_quantities(quantities)

    parallel_domain = Parallel_domain(points, vertices, boundary, **kwargs)


    #------------------------------------------------------------------------
    # Copy in quantity data
    #------------------------------------------------------------------------
    for q in quantities:
        try:
            parallel_domain.set_quantity(q, quantities[q], location='centroids')
        except KeyError:
            #print 'Try to create quantity %s'% q
            from anuga import Quantity
            Q = Quantity(parallel_domain, name=q, register=True)
            parallel_domain.set_quantity(q, quantities[q], location='centroids')

    #------------------------------------------------------------------------
    # Re-create the tracers on this sub-domain
    #------------------------------------------------------------------------
    # After the quantities: set_tracer derives m = h*c from this sub-domain's
    # own depth, so stage and elevation have to be in place first. The values
    # cover ghost triangles too, so the halo starts correct rather than waiting
    # for the first communicate_tracer_ghosts().
    restore_tracers(parallel_domain, tracer_metadata, tracer_values)

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


def distribute_basic_mesh(basic_mesh, verbose=False, parameters=None):
    """Distribute a Basic_mesh and return a Parallel_domain on each rank.

    This is the mesh-first alternative to distribute().  Only mesh topology
    is distributed -- no quantities.  The caller sets initial conditions,
    boundaries, operators and structures on the returned Parallel_domain.

    Parameters
    ----------
    basic_mesh : Basic_mesh
        The mesh to distribute.  Only rank 0 needs to supply a meaningful
        value; other ranks may pass None.
    verbose : bool, optional
    parameters : dict, optional
        Passed to partition_mesh / build_submesh (e.g. ghost layer width,
        partition_scheme).

    Returns
    -------
    Parallel_domain (numprocs > 1) or a plain Domain (numprocs == 1).

    Example
    -------
    >>> from anuga.parallel import distribute_mesh, finalize, myid
    >>> from anuga.abstract_2d_finite_volumes.basic_mesh import Basic_mesh
    >>> from anuga.abstract_2d_finite_volumes.pmesh2domain import pmesh_to_basic_mesh
    >>>
    >>> if myid == 0:
    ...     pmesh = create_mesh_from_regions(...)
    ...     bm = pmesh_to_basic_mesh(pmesh)
    ... else:
    ...     bm = None
    >>>
    >>> domain = distribute_mesh(bm)
    >>> domain.set_quantity('elevation', lambda x, y: -x / 100)
    >>> domain.set_quantity('stage', 0.0)
    >>> domain.set_boundary({'exterior': Reflective_boundary(domain)})
    >>> for t in domain.evolve(yieldstep=60, finaltime=3600):
    ...     pass
    >>> finalize()
    """
    from .distribute_mesh import partition_mesh, build_submesh, extract_submesh
    from .parallel_shallow_water import Parallel_domain

    if not pypar_available or numprocs == 1:
        # Serial fallback: build a plain Domain directly from the Basic_mesh.
        from anuga.shallow_water.shallow_water_domain import Domain
        return Domain(basic_mesh.nodes, basic_mesh.triangles,
                      boundary=basic_mesh.boundary,
                      geo_reference=basic_mesh.geo_reference)

    if myid == 0:
        domain_georef = basic_mesh.geo_reference
        n_global_tri  = basic_mesh.number_of_triangles
        n_global_nodes = basic_mesh.number_of_nodes

        # Partition the mesh topology (quantities will be empty).
        new_mesh, triangles_per_proc, quantities, s2p_map, p2s_map = \
            partition_mesh(basic_mesh, numprocs,
                           parameters=parameters, verbose=verbose)

        # Build ghost layers and communication patterns (no quantity data).
        submesh = build_submesh(new_mesh, quantities, triangles_per_proc,
                                parameters=parameters, verbose=verbose)

        # Defaults for domain metadata -- the user sets what they need after.
        _defaults = dict(
            boundary_map={},
            domain_name='domain',
            domain_dir='.',
            domain_store=True,
            domain_store_centroids=False,
            domain_minimum_storable_height=0.001,
            domain_minimum_allowed_height=1.0e-05,
            domain_flow_algorithm='1_5_order',
            domain_georef=domain_georef,
            domain_quantities_to_be_stored={
                'stage': 2, 'xmomentum': 2, 'ymomentum': 2},
            domain_smooth=True,
            domain_low_froude=0,
        )

        # Extract rank 0's submesh directly.
        (points, vertices, boundary, quantities0,
         ghost_recv_dict, full_send_dict,
         tri_map, node_map, tri_l2g, node_l2g, ghost_layer_width) = \
            extract_submesh(submesh, triangles_per_proc, p2s_map, 0)

        kwargs0 = {
            'full_send_dict':           full_send_dict,
            'ghost_recv_dict':          ghost_recv_dict,
            'number_of_full_nodes':     len(submesh['full_nodes'][0]),
            'number_of_full_triangles': len(submesh['full_triangles'][0]),
            'geo_reference':            domain_georef,
            'number_of_global_triangles': n_global_tri,
            'number_of_global_nodes':     n_global_nodes,
            'processor':    0,
            'numproc':      numprocs,
            's2p_map':      None,
            'p2s_map':      None,
            'tri_l2g':      tri_l2g,
            'node_l2g':     node_l2g,
            'ghost_layer_width': ghost_layer_width,
        }

        # Send all other ranks using the existing _isend_submesh protocol.
        # Build the full tostore tuple with empty quantities and defaults.
        from mpi4py import MPI as _MPI
        comm = _MPI.COMM_WORLD
        all_reqs = []
        for p in range(1, numprocs):
            (pts_p, verts_p, bnd_p, quant_p,
             ghost_recv_p, full_send_p,
             tri_map_p, node_map_p, tri_l2g_p, node_l2g_p, glw_p) = \
                extract_submesh(submesh, triangles_per_proc, p2s_map, p)

            kwargs_p = {
                'full_send_dict':           full_send_p,
                'ghost_recv_dict':          ghost_recv_p,
                'number_of_full_nodes':     len(submesh['full_nodes'][p]),
                'number_of_full_triangles': len(submesh['full_triangles'][p]),
                'geo_reference':            domain_georef,
                'number_of_global_triangles': n_global_tri,
                'number_of_global_nodes':     n_global_nodes,
                'processor':    p,
                'numproc':      numprocs,
                's2p_map':      None,
                'p2s_map':      None,
                'tri_l2g':      tri_l2g_p,
                'node_l2g':     node_l2g_p,
                'ghost_layer_width': glw_p,
            }
            tostore_p = (
                kwargs_p, pts_p, verts_p, bnd_p, {},
                _defaults['boundary_map'],
                _defaults['domain_name'],
                _defaults['domain_dir'],
                _defaults['domain_store'],
                _defaults['domain_store_centroids'],
                _defaults['domain_minimum_storable_height'],
                _defaults['domain_minimum_allowed_height'],
                _defaults['domain_flow_algorithm'],
                _defaults['domain_georef'],
                _defaults['domain_quantities_to_be_stored'],
                _defaults['domain_smooth'],
                _defaults['domain_low_froude'],
            )
            all_reqs.extend(_isend_submesh(comm, tostore_p, p))
        _MPI.Request.Waitall(all_reqs)

        kwargs   = kwargs0
        boundary_recv = boundary
        domain_name = _defaults['domain_name']

    else:
        from mpi4py import MPI as _MPI
        comm = _MPI.COMM_WORLD
        (kwargs, points, vertices, boundary_recv, quantities_recv,
         boundary_map_recv,
         domain_name, domain_dir, domain_store, domain_store_centroids,
         domain_minimum_storable_height, domain_minimum_allowed_height,
         domain_flow_algorithm, domain_georef,
         domain_quantities_to_be_stored, domain_smooth,
         domain_low_froude) = _recv_submesh(comm, 0)

    parallel_domain = Parallel_domain(points, vertices, boundary_recv,
                                      **kwargs)

    # Apply the rank suffix to the domain name so each rank writes to its own
    # SWW file (e.g. domain_P2_0.sww, domain_P2_1.sww).  The user can call
    # set_name() again afterwards to choose a different base name.
    parallel_domain.set_name(domain_name)

    return parallel_domain


def distribute_basic_mesh_collaborative(basic_mesh, verbose=False, parameters=None):
    """Distribute a Basic_mesh using the collaborative (shared-memory) approach.

    Combines the mesh-first interface of :func:`distribute_basic_mesh` with
    the scalable distribution strategy of :func:`distribute_collaborative`:

    1. Rank 0 partitions the mesh topology only (no quantities).
    2. The full mesh topology (nodes, triangles, neighbours) is broadcast to
       all ranks using ``MPI.Win.Allocate_shared`` — one physical copy per
       node, with a single Bcast between node leaders.
    3. Each rank independently builds its own ghost layer using the BFS
       Cython kernel.
    4. Each rank calls ``build_local_mesh`` locally (no quantity data).
    5. A ``Parallel_domain`` is returned with no quantities set; the caller
       sets initial conditions, boundaries, operators and structures after
       the call, exactly as with :func:`distribute_basic_mesh`.

    This avoids both the rank-0 ``build_submesh`` bottleneck of
    :func:`distribute_basic_mesh` (which builds all P submeshes sequentially
    on rank 0) and the need to distribute quantity arrays (which are set
    per-rank after the call).

    Parameters
    ----------
    basic_mesh : Basic_mesh or None
        The mesh to distribute.  Only rank 0 needs a real value; all other
        ranks should pass ``None``.
    verbose : bool, optional
        Print partitioning progress.  Default ``False``.
    parameters : dict, optional
        Partitioning options, e.g. ``{'partition_scheme': 'metis'}``.
        Also accepts ``'ghost_layer_width'`` (default 2).

    Returns
    -------
    Parallel_domain (numprocs > 1) or plain Domain (numprocs == 1).
    """
    from .distribute_mesh import partition_mesh, ghost_commun_pattern, build_local_mesh
    from .parallel_shallow_water import Parallel_domain

    if not pypar_available or numprocs == 1:
        from anuga.shallow_water.shallow_water_domain import Domain
        return Domain(basic_mesh.nodes, basic_mesh.triangles,
                      boundary=basic_mesh.boundary,
                      geo_reference=basic_mesh.geo_reference)

    from mpi4py import MPI
    comm     = MPI.COMM_WORLD
    myid     = comm.Get_rank()
    nprocs   = comm.Get_size()

    layer_width = (parameters or {}).get('ghost_layer_width', 2)

    node_comm = comm.Split_type(MPI.COMM_TYPE_SHARED)

    # ── Step 1: rank 0 partitions ─────────────────────────────────────────
    if myid == 0:
        new_mesh, tpp, _quantities, _s2p, _p2s = partition_mesh(
            basic_mesh, nprocs, parameters=parameters, verbose=verbose)

        domain_attrs = {
            'geo_reference':            basic_mesh.geo_reference,
            'num_global_tri':           basic_mesh.number_of_triangles,
            'num_global_nodes':         basic_mesh.number_of_nodes,
            'domain_name':              'domain',
        }
    else:
        new_mesh = tpp = domain_attrs = None

    # ── Step 2: broadcast small metadata ─────────────────────────────────
    domain_attrs = comm.bcast(domain_attrs, root=0)
    tpp          = comm.bcast(tpp,          root=0)

    # ── Step 3: shared-memory broadcast of mesh topology ─────────────────
    nodes,      _win_nodes = _shared_bcast_ndarray(
        new_mesh.nodes      if myid == 0 else None, comm, node_comm)
    triangles,  _win_tri   = _shared_bcast_ndarray(
        new_mesh.triangles  if myid == 0 else None, comm, node_comm)
    neighbours, _win_nbrs  = _shared_bcast_ndarray(
        new_mesh.neighbours.astype(num.int64, copy=False) if myid == 0 else None,
        comm, node_comm)
    boundary = comm.bcast(new_mesh.boundary if myid == 0 else None, root=0)

    # ── Step 4: each rank builds its own ghost layer ──────────────────────
    cumsum  = num.concatenate([[0], num.cumsum(tpp)])
    tlower  = int(cumsum[myid])
    tupper  = int(cumsum[myid + 1])

    neighbours_c = num.ascontiguousarray(neighbours, dtype=num.int64)

    from .distribute_mesh_ext import ghost_layer_bfs, ghost_bnd_layer_classify

    ghost_ids         = ghost_layer_bfs(neighbours_c, tlower, tupper, layer_width)
    tri_ids, edge_ids = ghost_bnd_layer_classify(neighbours_c, ghost_ids,
                                                 tlower, tupper)

    # ── Step 5: assemble local submesh ────────────────────────────────────
    full_triangles  = triangles[tlower:tupper].copy()
    full_node_ids   = num.unique(full_triangles.flat)
    full_nodes_col  = num.concatenate(
        [full_node_ids.reshape(-1, 1), nodes[full_node_ids]], axis=1)

    ghost_tri_verts = triangles[ghost_ids]
    ghost_node_ids  = num.setdiff1d(num.unique(ghost_tri_verts.flat), full_node_ids)
    ghost_nodes_col = num.concatenate(
        [ghost_node_ids.reshape(-1, 1), nodes[ghost_node_ids]], axis=1)
    ghost_tris_col  = num.concatenate(
        [ghost_ids.reshape(-1, 1), ghost_tri_verts], axis=1)

    # Release shared-memory windows.
    for _win in (_win_nodes, _win_tri, _win_nbrs):
        if _win is not None:
            _win.Free()
    node_comm.Free()

    # Full boundary: binary-search into the sorted global boundary dict.
    sorted_bnd_keys = sorted(boundary.keys(), key=lambda k: k[0])
    if sorted_bnd_keys:
        sorted_bnd_tids = num.array([k[0] for k in sorted_bnd_keys], dtype=int)
        lo = int(num.searchsorted(sorted_bnd_tids, tlower))
        hi = int(num.searchsorted(sorted_bnd_tids, tupper))
        full_boundary = {sorted_bnd_keys[i]: boundary[sorted_bnd_keys[i]]
                         for i in range(lo, hi)}
    else:
        full_boundary = {}

    ghost_boundary = {(int(t), int(e)): 'ghost'
                      for t, e in zip(tri_ids, edge_ids)}
    ghost_boundary.update((k, boundary[k])
                          for k in ghost_boundary.keys() & boundary.keys())

    tri_per_proc_ranges = cumsum[1:] - 1
    ghost_commun = ghost_commun_pattern(ghost_tris_col, myid, tri_per_proc_ranges)

    # ── Step 6: allgather ghost_commun → derive full_commun ──────────────
    all_ghost_commun = comm.allgather(ghost_commun)

    full_commun = {}
    for q, gc in enumerate(all_ghost_commun):
        for row in gc:
            gid   = int(row[0])
            owner = int(row[1])
            if owner == myid:
                if gid not in full_commun:
                    full_commun[gid] = []
                full_commun[gid].append(q)

    # ── Step 7: build local mesh (no quantities) ──────────────────────────
    submesh_cell = {
        'ghost_layer_width': layer_width,
        'full_nodes':        full_nodes_col,
        'ghost_nodes':       ghost_nodes_col,
        'full_triangles':    full_triangles,
        'ghost_triangles':   ghost_tris_col,
        'full_boundary':     full_boundary,
        'ghost_boundary':    ghost_boundary,
        'ghost_commun':      ghost_commun,
        'full_commun':       full_commun,
        'full_quan':         {},
        'ghost_quan':        {},
    }

    points, vertices, boundary_local, _quantities_local, ghost_recv_dict, \
        full_send_dict, _tri_map, _node_map, tri_l2g, node_l2g, glw = \
        build_local_mesh(submesh_cell, tlower, tupper, nprocs)

    number_of_full_nodes     = len(full_nodes_col)
    number_of_full_triangles = tupper - tlower

    # ── Step 8: create Parallel_domain ───────────────────────────────────
    kwargs = {
        'full_send_dict':             full_send_dict,
        'ghost_recv_dict':            ghost_recv_dict,
        'number_of_full_nodes':       number_of_full_nodes,
        'number_of_full_triangles':   number_of_full_triangles,
        'geo_reference':              domain_attrs['geo_reference'],
        'number_of_global_triangles': domain_attrs['num_global_tri'],
        'number_of_global_nodes':     domain_attrs['num_global_nodes'],
        'processor':                  myid,
        'numproc':                    nprocs,
        's2p_map':                    None,
        'p2s_map':                    None,
        'tri_l2g':                    tri_l2g,
        'node_l2g':                   node_l2g,
        'ghost_layer_width':          glw,
    }

    parallel_domain = Parallel_domain(points, vertices, boundary_local, **kwargs)
    parallel_domain.set_name(domain_attrs['domain_name'])

    return parallel_domain


def _shared_bcast_ndarray(arr_on_root, comm, node_comm):
    """Broadcast a numpy array from global rank 0 using shared memory within nodes.

    Within each node, a single physical copy is held in an
    ``MPI.Win.Allocate_shared`` window; all ranks on the node map a numpy
    view into that window.  Between nodes, a Bcast among node leaders
    (rank 0 within each node) distributes the data.

    Parameters
    ----------
    arr_on_root : ndarray or None
        Array to broadcast.  Only meaningful on global rank 0; pass None
        on all other ranks.
    comm : MPI.Comm
        Global communicator.
    node_comm : MPI.Comm
        Shared-memory communicator obtained via
        ``comm.Split_type(MPI.COMM_TYPE_SHARED)``.

    Returns
    -------
    arr : ndarray
        View into shared memory (valid until ``win.Free()`` is called).
    win : MPI.Win
        Shared memory window.  Caller must call ``win.Free()`` when
        finished with ``arr``.
    """
    from mpi4py import MPI
    myid      = comm.Get_rank()
    node_rank = node_comm.Get_rank()

    # Broadcast shape/dtype from global rank 0.
    meta = (arr_on_root.shape, arr_on_root.dtype.str) if myid == 0 else None
    shape, dtype_str = comm.bcast(meta, root=0)
    dtype  = num.dtype(dtype_str)
    nbytes = int(num.prod(shape)) * dtype.itemsize

    # Attempt shared-memory allocation.  Fall back to a plain private buffer
    # if Win.Allocate_shared is unsupported (e.g. some MPI configurations).
    try:
        # Node leaders allocate the full buffer; non-leaders allocate nothing.
        win = MPI.Win.Allocate_shared(
            nbytes if node_rank == 0 else 0,
            dtype.itemsize, MPI.INFO_NULL, node_comm)

        # All ranks on the node get a numpy view into the leader's buffer.
        buf, _ = win.Shared_query(0)
        arr = num.ndarray(shape, dtype=dtype, buffer=buf)

        # Node leaders form a communicator for inter-node Bcast.
        leader_comm = comm.Split(0 if node_rank == 0 else MPI.UNDEFINED, myid)
        if node_rank == 0:
            if myid == 0:
                arr[...] = arr_on_root
            leader_comm.Bcast(arr, root=0)
            leader_comm.Free()

        # Barrier so all node members see the populated data before returning.
        node_comm.Barrier()

    except (MPI.Exception, NotImplementedError):
        # Shared memory unavailable — fall back to a private copy via Bcast.
        win = None
        arr = num.empty(shape, dtype=dtype)
        if myid == 0:
            arr[...] = arr_on_root
        comm.Bcast(arr, root=0)

    return arr, win


def _scatterv_quantity(arr_on_root, cumsum, comm):
    """Scatter a 1-D float64 quantity array: rank p receives arr[cumsum[p]:cumsum[p+1]].

    Parameters
    ----------
    arr_on_root : 1-D float64 ndarray or None
        Full quantity array on rank 0; None on all other ranks.
    cumsum : 1-D int array, length numprocs+1
        Cumulative triangle counts (``cumsum[p]`` = first triangle of rank p).
    comm : MPI.Comm

    Returns
    -------
    recvbuf : 1-D float64 ndarray, length ``cumsum[myid+1] - cumsum[myid]``
    """
    from mpi4py import MPI
    myid  = comm.Get_rank()
    nproc = comm.Get_size()
    count   = int(cumsum[myid + 1] - cumsum[myid])
    recvbuf = num.empty(count, dtype=num.float64)
    if myid == 0:
        sendcounts = [int(cumsum[p + 1] - cumsum[p]) for p in range(nproc)]
        displs     = [int(cumsum[p])                  for p in range(nproc)]
        comm.Scatterv(
            [arr_on_root.astype(num.float64, copy=False), sendcounts, displs, MPI.DOUBLE],
            recvbuf, root=0)
    else:
        comm.Scatterv(None, recvbuf, root=0)
    return recvbuf


def distribute_collaborative(domain, verbose=False, debug=False, parameters=None):
    """MPI-collaborative domain distribution — eliminates the rank-0 bottleneck.

    The standard :func:`distribute` function builds ALL P submeshes sequentially
    on rank 0 — O(P × N) work.  This function instead:

    1. Rank 0 partitions the mesh (METIS / Morton / Hilbert).
    2. Mesh topology (nodes, triangles, neighbours) is distributed using
       ``MPI.Win.Allocate_shared``: one physical copy per node, shared by
       all ranks on that node.  Between nodes a single Bcast among node
       leaders transfers the data.
    3. Each rank independently builds its own ghost layer using the
       Cython BFS kernel (``ghost_layer_bfs`` / ``ghost_bnd_layer_classify``).
    4. Full-triangle quantities are distributed with ``MPI_Scatterv``
       (each rank receives only its O(N/P) rows).  Ghost quantities are
       sent by rank 0 via targeted ``MPI_Isend`` after an
       ``MPI_Allgather`` of ghost-communication arrays reveals every
       rank's ghost triangle IDs.
    5. Each rank calls ``build_local_mesh`` locally.

    Memory per rank: O(N/nodes × topology) + O(N/P + G) for quantities,
    compared to O(N) per rank in the original ``distribute_collaborative``.
    The only serial bottleneck remaining on rank 0 is the METIS partition
    call, which is O(N) and typically a small fraction of total setup time.

    Parameters
    ----------
    domain : Domain
        The full serial domain (only meaningful on rank 0; other ranks may
        pass any Domain object — it will not be read).
    verbose : bool, optional
    debug : bool, optional
    parameters : dict, optional
        Passed through to ``partition_mesh``.  Recognised extra keys:

        ``ghost_layer_width`` : int, default 2

    Returns
    -------
    Parallel_domain
        The local parallel domain for this MPI rank.
    """

    if not pypar_available:
        return domain

    from mpi4py import MPI
    comm  = MPI.COMM_WORLD
    myid  = comm.Get_rank()
    numprocs = comm.Get_size()

    if numprocs == 1:
        return domain

    if debug:
        verbose = True

    layer_width = (parameters or {}).get('ghost_layer_width', 2)

    # Shared-memory communicator: all ranks on the same node.
    # Used by _shared_bcast_ndarray to avoid P copies of the mesh topology.
    node_comm = comm.Split_type(MPI.COMM_TYPE_SHARED)

    # Tracers ride the quantity machinery under reserved keys (#278); the
    # names and shared beta go in domain_attrs, which is broadcast.
    from anuga.parallel.distribute_mesh import (
        tracer_partition_metadata, pop_tracer_quantities, restore_tracers)

    # ── Step 1: rank 0 partitions ────────────────────────────────────────
    if myid == 0:
        new_mesh, tpp, quantities, _s2p, _p2s = partition_mesh(
            domain, numprocs, parameters=parameters, verbose=verbose)

        # Build dummy boundary_map from boundary tags (same as Sequential_distribute).
        # domain.boundary_map may be None if set_boundary() hasn't been called yet.
        bdmap = domain.boundary_map
        if bdmap is None:
            bdmap = {tag: None for tag in domain.get_boundary_tags()}
            domain.set_boundary(bdmap)

        domain_attrs = {
            'name':                     domain.get_name(),
            'dir':                      domain.get_datadir(),
            'store':                    domain.get_store(),
            'store_centroids':          domain.get_store_centroids(),
            'min_store_height':         domain.minimum_storable_height,
            'min_allowed_height':       domain.get_minimum_allowed_height(),
            'flow_algorithm':           domain.get_flow_algorithm(),
            'geo_reference':            domain.geo_reference,
            'quantities_to_be_stored':  domain.quantities_to_be_stored,
            'smooth':                   domain.smooth,
            'low_froude':               domain.low_froude,
            'boundary_map':             domain.boundary_map,
            'num_global_tri':           domain.number_of_triangles,
            'num_global_nodes':         domain.number_of_nodes,
            # Tracer names and shared beta. The concentrations themselves are
            # already in `quantities` under reserved keys, so they ride the
            # Scatterv below like any other quantity. #278
            'tracer_metadata':          tracer_partition_metadata(domain),
        }
        quant_keys = list(quantities.keys())
    else:
        new_mesh = tpp = quantities = domain_attrs = quant_keys = None

    # ── Step 2: broadcast small metadata ─────────────────────────────────
    domain_attrs = comm.bcast(domain_attrs, root=0)
    tpp          = comm.bcast(tpp,          root=0)
    quant_keys   = comm.bcast(quant_keys,   root=0)

    # ── Step 3: shared-memory broadcast of mesh topology ─────────────────
    # One physical copy per node via MPI.Win.Allocate_shared; inter-node
    # transfer is a Bcast among node leaders only.
    nodes,      _win_nodes = _shared_bcast_ndarray(
        new_mesh.nodes      if myid == 0 else None, comm, node_comm)
    triangles,  _win_tri   = _shared_bcast_ndarray(
        new_mesh.triangles  if myid == 0 else None, comm, node_comm)
    # Cast to int64 before sharing so that ascontiguousarray(neighbours, int64)
    # below is a zero-copy no-op (avoids a full private copy per rank).
    neighbours, _win_nbrs  = _shared_bcast_ndarray(
        new_mesh.neighbours.astype(num.int64, copy=False) if myid == 0 else None,
        comm, node_comm)
    boundary = comm.bcast(new_mesh.boundary if myid == 0 else None, root=0)

    # ── Step 4: Scatterv full-triangle quantities ─────────────────────────
    # Each rank receives only its own N/P rows — O(N/P) per rank vs O(N).
    # Ghost quantities are fetched after the BFS identifies ghost_ids.

    # ── Step 5: each rank builds its own ghost layer ──────────────────────
    cumsum  = num.concatenate([[0], num.cumsum(tpp)])
    tlower  = int(cumsum[myid])
    tupper  = int(cumsum[myid + 1])

    neighbours_c = num.ascontiguousarray(neighbours, dtype=num.int64)

    from .distribute_mesh_ext import ghost_layer_bfs, ghost_bnd_layer_classify

    ghost_ids         = ghost_layer_bfs(neighbours_c, tlower, tupper, layer_width)
    tri_ids, edge_ids = ghost_bnd_layer_classify(neighbours_c, ghost_ids,
                                                 tlower, tupper)

    # ── Step 6: assemble local submesh_cell ──────────────────────────────
    # Extract all needed data as private copies from shared arrays, then
    # free the windows so shared memory is released promptly.

    full_triangles  = triangles[tlower:tupper].copy()   # copy: slice is a view
    full_node_ids   = num.unique(full_triangles.flat)
    full_nodes_col  = num.concatenate(                  # (N_fn, 3): [gid, x, y]
        [full_node_ids.reshape(-1, 1), nodes[full_node_ids]], axis=1)

    ghost_tri_verts = triangles[ghost_ids]              # fancy-index → already a copy
    ghost_node_ids  = num.setdiff1d(num.unique(ghost_tri_verts.flat), full_node_ids)
    ghost_nodes_col = num.concatenate(                  # (G_n, 3): [gid, x, y]
        [ghost_node_ids.reshape(-1, 1), nodes[ghost_node_ids]], axis=1)
    ghost_tris_col  = num.concatenate(                  # (G, 4): [gid, v0, v1, v2]
        [ghost_ids.reshape(-1, 1), ghost_tri_verts], axis=1)

    # Topology arrays no longer needed — release shared-memory windows.
    # (win is None when the shared-memory fallback path was used.)
    for _win in (_win_nodes, _win_tri, _win_nbrs):
        if _win is not None:
            _win.Free()
    node_comm.Free()

    # Full boundary: binary-search into the sorted global boundary dict
    sorted_bnd_keys = sorted(boundary.keys(), key=lambda k: k[0])
    if sorted_bnd_keys:
        sorted_bnd_tids = num.array([k[0] for k in sorted_bnd_keys], dtype=int)
        lo = int(num.searchsorted(sorted_bnd_tids, tlower))
        hi = int(num.searchsorted(sorted_bnd_tids, tupper))
        full_boundary = {sorted_bnd_keys[i]: boundary[sorted_bnd_keys[i]]
                         for i in range(lo, hi)}
    else:
        full_boundary = {}

    # Ghost boundary: start with 'ghost', override with real tags where present
    ghost_boundary = {(int(t), int(e)): 'ghost'
                      for t, e in zip(tri_ids, edge_ids)}
    ghost_boundary.update((k, boundary[k])
                          for k in ghost_boundary.keys() & boundary.keys())

    # Ghost communication pattern: ghost_tris_col[:, 0] → owning processor
    tri_per_proc_ranges = cumsum[1:] - 1
    ghost_commun = ghost_commun_pattern(ghost_tris_col, myid, tri_per_proc_ranges)

    # ── Step 7: allgather ghost_commun → each rank derives full_commun ───
    # ghost_commun has shape (G_p, 2): col0 = global_id, col1 = owner_proc.
    # After allgather, rank myid filters for entries where owner == myid.
    # The allgather result also tells rank 0 every other rank's ghost_ids,
    # which is used below to send targeted ghost quantities.
    all_ghost_commun = comm.allgather(ghost_commun)

    full_commun = {}
    for q, gc in enumerate(all_ghost_commun):
        for row in gc:
            gid   = int(row[0])
            owner = int(row[1])
            if owner == myid:
                if gid not in full_commun:
                    full_commun[gid] = []
                full_commun[gid].append(q)

    # ── Step 4b: full-triangle quantities via Scatterv ────────────────────
    # quantities[k] has shape (N, 1); Scatterv works on flattened rows so we
    # pass the 1-D view and reshape back to (N/P, 1) on receipt.
    quant_full = {k: _scatterv_quantity(
                        quantities[k][:, 0] if myid == 0 else None, cumsum, comm
                     ).reshape(-1, 1)
                  for k in quant_keys}

    # ── Step 4c: ghost quantities via targeted Isend/Recv ─────────────────
    # Rank 0 has full quantities; it knows every rank's ghost_ids from
    # all_ghost_commun (col 0 = global triangle id).
    all_ghost_ids = [gc[:, 0].astype(num.int64) if len(gc) > 0
                     else num.empty(0, dtype=num.int64)
                     for gc in all_ghost_commun]

    quant_ghost = {}
    for ki, k in enumerate(quant_keys):
        tag = ki
        if myid == 0:
            # quantities[k] has shape (N, 1); flatten to 1-D for send.
            q = quantities[k][:, 0].astype(num.float64, copy=False)
            reqs = [comm.Isend(q[all_ghost_ids[p]], dest=p, tag=tag)
                    for p in range(1, numprocs)]
            quant_ghost[k] = q[ghost_ids].reshape(-1, 1).copy()
            MPI.Request.Waitall(reqs)
        else:
            G   = len(ghost_ids)
            buf = num.empty(G, dtype=num.float64)
            if G > 0:
                comm.Recv(buf, source=0, tag=tag)
            quant_ghost[k] = buf.reshape(-1, 1)

    # ── Step 8: build local mesh in GA format ────────────────────────────
    submesh_cell = {
        'ghost_layer_width': layer_width,
        'full_nodes':        full_nodes_col,
        'ghost_nodes':       ghost_nodes_col,
        'full_triangles':    full_triangles,
        'ghost_triangles':   ghost_tris_col,
        'full_boundary':     full_boundary,
        'ghost_boundary':    ghost_boundary,
        'ghost_commun':      ghost_commun,
        'full_commun':       full_commun,
        'full_quan':         quant_full,
        'ghost_quan':        quant_ghost,
    }

    points, vertices, boundary_local, quantities_local, ghost_recv_dict, \
        full_send_dict, tri_map, node_map, tri_l2g, node_l2g, glw = \
        build_local_mesh(submesh_cell, tlower, tupper, numprocs)

    number_of_full_nodes     = len(full_nodes_col)
    number_of_full_triangles = tupper - tlower

    # ── Step 9: create Parallel_domain ───────────────────────────────────
    kwargs = {
        'full_send_dict':             full_send_dict,
        'ghost_recv_dict':            ghost_recv_dict,
        'number_of_full_nodes':       number_of_full_nodes,
        'number_of_full_triangles':   number_of_full_triangles,
        'geo_reference':              domain_attrs['geo_reference'],
        'number_of_global_triangles': domain_attrs['num_global_tri'],
        'number_of_global_nodes':     domain_attrs['num_global_nodes'],
        'processor':                  myid,
        'numproc':                    numprocs,
        's2p_map':                    None,
        'p2s_map':                    None,
        'tri_l2g':                    tri_l2g,
        'node_l2g':                   node_l2g,
        'ghost_layer_width':          glw,
    }

    parallel_domain = Parallel_domain(points, vertices, boundary_local, **kwargs)

    # Out before the loop below, which would otherwise register each one as a
    # Quantity called `__tracer__0`. #278
    tracer_values = pop_tracer_quantities(quantities_local)

    for q in quantities_local:
        try:
            parallel_domain.set_quantity(q, quantities_local[q], location='centroids')
        except KeyError:
            from anuga import Quantity
            Q = Quantity(parallel_domain, name=q, register=True)
            parallel_domain.set_quantity(q, quantities_local[q], location='centroids')

    # After the quantities: set_tracer derives m = h*c from the local depth.
    restore_tracers(parallel_domain,
                    domain_attrs.get('tracer_metadata'), tracer_values)

    boundary_map = domain_attrs['boundary_map'].copy()
    boundary_map['ghost'] = None
    parallel_domain.set_boundary(boundary_map)

    parallel_domain.set_name(domain_attrs['name'])
    parallel_domain.set_datadir(domain_attrs['dir'])
    parallel_domain.set_store(domain_attrs['store'])
    parallel_domain.set_store_centroids(domain_attrs['store_centroids'])
    parallel_domain.set_minimum_storable_height(domain_attrs['min_store_height'])
    parallel_domain.set_minimum_allowed_height(domain_attrs['min_allowed_height'])
    parallel_domain.set_flow_algorithm(domain_attrs['flow_algorithm'])
    parallel_domain.geo_reference = domain_attrs['geo_reference']
    parallel_domain.set_quantities_to_be_stored(domain_attrs['quantities_to_be_stored'])
    parallel_domain.smooth = domain_attrs['smooth']
    parallel_domain.set_low_froude(domain_attrs['low_froude'])

    return parallel_domain


def distribute_mesh(domain, verbose=False, debug=False, parameters=None):
    """ Distribute and send the mesh info to all the processors.
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
    new_mesh, triangles_per_proc, quantities, \
           s2p_map, p2s_map = \
           partition_mesh(domain, numprocs)


    # Build the mesh that should be assigned to each processor,
    # this includes ghost nodes and the communication pattern
    if verbose: print('Build submeshes')
    submesh = build_submesh(new_mesh, quantities, triangles_per_proc, parameters)

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


    # Return structures necessary for building the parallel domain
    return points, vertices, boundary, quantities,\
           ghost_recv_dict, full_send_dict,\
           number_of_full_nodes, number_of_full_triangles, \
           s2p_map, p2s_map, tri_map, node_map, tri_l2g, node_l2g, \
           ghost_layer_width



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
